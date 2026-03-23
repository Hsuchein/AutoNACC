"""
autoprep.normalize - GaussView PDB 规范化

功能:
  1. HETATM → ATOM, 统一残基名/链/编号
  2. 原子重排: N, CA, C, O → 与主链相连的H → CB → 其余原子 (原序)
  3. 原子命名: 保留主链名 (N, CA, C, O, CB), 其余按元素编号
  4. 可选保留 CONECT (序号自动映射)
"""

import math
from collections import defaultdict

BACKBONE_NAMES = {"N", "CA", "C", "O", "CB"}
BACKBONE_ORDER = ["N", "CA", "C", "O"]   # 固定顺序

# 共价半径 (Å), 用于键检测
_COV_RADII = {
    'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05,
    'P': 1.07, 'F': 0.57, 'Cl': 0.99, 'Br': 1.14, 'Se': 1.20,
}


def _dist(c1, c2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))


def _parse_element(line, name):
    """从 PDB 行提取元素符号"""
    elem = line[76:78].strip()
    if not elem:
        elem = name[0]
    return elem


def _parse_coord(line):
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def _is_bonded(elem1, coord1, elem2, coord2, factor=1.3):
    r1 = _COV_RADII.get(elem1, 0.77)
    r2 = _COV_RADII.get(elem2, 0.77)
    return _dist(coord1, coord2) < factor * (r1 + r2)


def _reorder_atoms(atom_lines):
    """
    重排原子顺序:
      1. N, CA, C, O  (主链重原子, 固定顺序)
      2. 与 N/CA/C/O 相连的 H  (按原始索引排)
      3. CB
      4. 其余原子 (按原始索引排)
    """
    # 解析所有原子信息
    atoms = []
    for i, (line, is_orig) in enumerate(atom_lines):
        name = line[12:16].strip()
        elem = _parse_element(line, name)
        coord = _parse_coord(line)
        atoms.append({
            'idx': i, 'line': line, 'is_orig': is_orig,
            'name': name, 'elem': elem, 'coord': coord,
        })

    # 找到主链原子 (仅 orig ATOM 记录中)
    bb_map = {}   # name -> atom dict
    for a in atoms:
        if a['is_orig'] and a['name'] in set(BACKBONE_ORDER):
            bb_map[a['name']] = a

    # 检测与主链 N/CA/C/O 相连的 H
    bb_h_indices = set()
    for bb_name in BACKBONE_ORDER:
        if bb_name not in bb_map:
            continue
        bb = bb_map[bb_name]
        for a in atoms:
            if a['elem'] == 'H' and a['idx'] != bb['idx']:
                if _is_bonded(bb['elem'], bb['coord'], a['elem'], a['coord']):
                    bb_h_indices.add(a['idx'])

    # 找 CB
    cb_atom = None
    for a in atoms:
        if a['is_orig'] and a['name'] == 'CB':
            cb_atom = a
            break

    # 构建有序列表
    ordered = []
    used = set()

    # Group 1: N, CA, C, O
    for name in BACKBONE_ORDER:
        if name in bb_map:
            ordered.append(bb_map[name])
            used.add(bb_map[name]['idx'])

    # Group 2: H connected to backbone (按原始索引排)
    bb_hs = sorted([atoms[i] for i in bb_h_indices], key=lambda a: a['idx'])
    for a in bb_hs:
        if a['idx'] not in used:
            ordered.append(a)
            used.add(a['idx'])

    # Group 3: CB
    if cb_atom and cb_atom['idx'] not in used:
        ordered.append(cb_atom)
        used.add(cb_atom['idx'])

    # Group 4: 其余原子按原始顺序
    for a in atoms:
        if a['idx'] not in used:
            ordered.append(a)

    return [(a['line'], a['is_orig']) for a in ordered]


def normalize_pdb(input_pdb, output_pdb, resname="CYM", chain="A", resnum=293,
                  keep_conect=True, reorder=True):
    with open(input_pdb) as f:
        lines = f.readlines()

    atom_lines = []   # (line, is_original_atom)
    conect_lines = []
    old_serials = []  # 原始序号按出现顺序

    for line in lines:
        rec = line[:6].strip()
        if rec == "ATOM":
            old_serials.append(int(line[6:11]))
            atom_lines.append((line, True))
        elif rec == "HETATM":
            old_serials.append(int(line[6:11]))
            atom_lines.append((line, False))
        elif rec == "CONECT":
            conect_lines.append(line)

    # 重排原子
    if reorder:
        atom_lines = _reorder_atoms(atom_lines)

    # 构建 old_serial → new_serial 映射
    # 先建 line → old_serial 的映射 (用原始行文本匹配)
    line_to_old = {}
    for line in lines:
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM"):
            line_to_old[id(line)] = int(line[6:11])

    # 重排后需要重新建映射: 遍历重排后的 atom_lines, 取其 line 的 old serial
    old_to_new = {}
    elem_count = defaultdict(int)
    out_atoms = []

    for i, (line, is_orig_atom) in enumerate(atom_lines):
        new_serial = i + 1
        old_serial = int(line[6:11])
        old_to_new[old_serial] = new_serial

        old_name = line[12:16].strip()
        element = line[76:78].strip()
        if not element:
            element = old_name[0]

        if is_orig_atom and old_name in BACKBONE_NAMES:
            atom_name = old_name
        else:
            elem_count[element] += 1
            atom_name = f"{element}{elem_count[element]}"

        if len(atom_name) <= 3:
            name_field = f" {atom_name:<3s}"
        else:
            name_field = f"{atom_name:<4s}"

        out = (f"ATOM  {new_serial:5d} {name_field}"
               f" {resname:>3s} {chain:1s}{resnum:4d}    "
               f"{line[30:54]}"
               f"{1.0:6.2f}{0.0:6.2f}"
               f"          {element:>2s}  \n")
        out_atoms.append(out)

    out_conect = []
    if keep_conect and conect_lines:
        for line in conect_lines:
            tokens = line.split()
            new_tokens = ["CONECT"]
            valid = True
            for t in tokens[1:]:
                old_s = int(t)
                if old_s not in old_to_new:
                    valid = False
                    break
                new_tokens.append(f"{old_to_new[old_s]:5d}")
            if valid:
                out_conect.append("".join(new_tokens) + "\n")

    with open(output_pdb, "w") as f:
        for line in out_atoms:
            f.write(line)
        for line in out_conect:
            f.write(line)
        f.write("END\n")

    print(f"Wrote {len(out_atoms)} atoms to {output_pdb}")
    print(f"Residue name: {resname}, Chain: {chain}, ResNum: {resnum}")

    print(f"\n{'Idx':>4s}  {'Name':>4s}  {'Elem':>4s}  {'X':>8s}  {'Y':>8s}  {'Z':>8s}")
    print("-" * 44)
    for line in out_atoms:
        serial = int(line[6:11])
        name = line[12:16].strip()
        elem = line[76:78].strip()
        x = line[30:38].strip()
        y = line[38:46].strip()
        z = line[46:54].strip()
        print(f"{serial:4d}  {name:>4s}  {elem:>4s}  {x:>8s}  {y:>8s}  {z:>8s}")
