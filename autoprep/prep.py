"""
autoprep.prep - 非标准氨基酸自动参数化核心逻辑

工作流:
    1. 读取残基片段(肽键断开,有自由价) + ACE/NME帽基模板
    2. 组装 ACE-RES-NME, 输出 PDB (含 CONECT), 确保肽键平面性
    3. RDKit MMFF94 力场预优化
    4. 生成 Gaussian 输入 (--Link1-- 串联两步):
       Job1: Opt opt_method (如 wB97XD/6-31G*)
       Job2: SP  resp_method Geom=AllCheck (如 HF/6-31G* Pop=MK)
    5. 运行 Gaussian, 每分钟检查是否完成
    6. antechamber RESP 电荷拟合 (单构象)
    7. prepgen 生成 .prepin
    8. parmchk2 生成 .frcmod
"""

import numpy as np
import subprocess as sp
import time
import os
import sys
import json
from pathlib import Path
from collections import defaultdict

# =============================================================================
# 模板路径
# =============================================================================
_TEMPLATE_DIR = Path(__file__).parent / "template"
ACE_TEMPLATE = _TEMPLATE_DIR / "ace.xyz"
NME_TEMPLATE = _TEMPLATE_DIR / "nme.xyz"

# =============================================================================
# 常量
# =============================================================================
COVALENT_RADII = {
    'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05,
    'P': 1.07, 'F': 0.57, 'Cl': 0.99, 'Br': 1.14, 'Se': 1.20,
}
PEPTIDE_CN = 1.335   # 肽键 C-N 键长 (Å)


# =============================================================================
# I/O
# =============================================================================

def read_xyz(path):
    lines = open(path).readlines()
    n = int(lines[0].strip())
    elems, coords = [], []
    for line in lines[2:2 + n]:
        p = line.split()
        elems.append(p[0])
        coords.append([float(x) for x in p[1:4]])
    return elems, np.array(coords)


def read_pdb(path):
    elems, coords, names = [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                names.append(line[12:16].strip())
                coords.append([float(line[30:38]), float(line[38:46]),
                               float(line[46:54])])
                elem = line[76:78].strip()
                if not elem:
                    nm = line[12:16].strip()
                    elem = ''.join(c for c in nm if c.isalpha())[:2]
                    if len(elem) == 2 and elem[0] in 'CH' and elem[1].isupper():
                        elem = elem[0]
                elems.append(elem)
    return elems, np.array(coords), names


def read_structure(path):
    ext = Path(path).suffix.lower()
    if ext == '.pdb':
        return read_pdb(path)
    elif ext == '.xyz':
        elems, coords = read_xyz(path)
        ctr = {}
        names = []
        for e in elems:
            ctr[e] = ctr.get(e, 0) + 1
            names.append(f"{e}{ctr[e]}")
        return elems, coords, names
    else:
        raise ValueError(f"不支持的格式: {ext}")


def write_xyz(path, elems, coords, title=""):
    with open(path, 'w') as f:
        f.write(f"{len(elems)}\n{title}\n")
        for e, (x, y, z) in zip(elems, coords):
            f.write(f" {e:<2s}  {x:14.8f}  {y:14.8f}  {z:14.8f}\n")


def get_all_bonds(coords, elems):
    """基于共价半径检测所有键, 返回 [(i, j), ...] (0-indexed)."""
    bonds = []
    n = len(coords)
    for i in range(n):
        ri = COVALENT_RADII.get(elems[i], 0.77)
        for j in range(i + 1, n):
            rj = COVALENT_RADII.get(elems[j], 0.77)
            if np.linalg.norm(coords[j] - coords[i]) < 1.3 * (ri + rj):
                bonds.append((i, j))
    return bonds


def write_assembled_pdb(path, elems, coords, bonds, rn="RES",
                        seg_labels=None):
    """输出带 CONECT 的 PDB.

    seg_labels: 每个原子属于哪一段 ("ACE", "RES", "NME"), 用于 segment ID.
    """
    # 生成唯一原子名
    ctr = {}
    atom_names = []
    for e in elems:
        ctr[e] = ctr.get(e, 0) + 1
        atom_names.append(f"{e}{ctr[e]}")

    with open(path, 'w') as f:
        f.write(f"REMARK  ACE-{rn}-NME assembled structure\n")
        for i, (e, (x, y, z)) in enumerate(zip(elems, coords)):
            serial = i + 1
            nm = atom_names[i]
            if len(nm) <= 3:
                name_field = f" {nm:<3s}"
            else:
                name_field = f"{nm:<4s}"
            seg = seg_labels[i] if seg_labels else rn
            f.write(f"HETATM{serial:5d} {name_field}"
                    f" {seg:>3s} A   1    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}"
                    f"{1.0:6.2f}{0.0:6.2f}"
                    f"          {e:>2s}  \n")

        # CONECT
        bond_map = defaultdict(list)
        for i, j in bonds:
            bond_map[i].append(j)
            bond_map[j].append(i)
        for atom_idx in sorted(bond_map.keys()):
            serial = atom_idx + 1
            neighbors = sorted([n + 1 for n in bond_map[atom_idx]])
            for chunk_start in range(0, len(neighbors), 4):
                chunk = neighbors[chunk_start:chunk_start + 4]
                line = f"CONECT{serial:5d}"
                for n in chunk:
                    line += f"{n:5d}"
                f.write(line + "\n")
        f.write("END\n")


# =============================================================================
# 几何工具
# =============================================================================

def _norm(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-10 else np.zeros(3)


def find_bonded(coords, elems, idx):
    bonded = []
    ri = COVALENT_RADII.get(elems[idx], 0.77)
    for j in range(len(coords)):
        if j == idx:
            continue
        rj = COVALENT_RADII.get(elems[j], 0.77)
        if np.linalg.norm(coords[j] - coords[idx]) < 1.3 * (ri + rj):
            bonded.append(j)
    return bonded


def find_bonded_elem(coords, elems, idx, target_elem):
    for j in find_bonded(coords, elems, idx):
        if elems[j] == target_elem:
            return j
    return None


def outgoing_direction(coords, elems, idx):
    bonded = find_bonded(coords, elems, idx)
    if not bonded:
        return np.array([1.0, 0.0, 0.0])
    vecs = [_norm(coords[j] - coords[idx]) for j in bonded]
    return _norm(-np.mean(vecs, axis=0))


def rotation_matrix_align(a, b):
    a, b = _norm(a), _norm(b)
    v = np.cross(a, b)
    c = np.dot(a, b)
    if np.linalg.norm(v) < 1e-10:
        if c > 0:
            return np.eye(3)
        perp = np.array([1, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1, 0])
        perp = _norm(perp - np.dot(perp, a) * a)
        return 2 * np.outer(perp, perp) - np.eye(3)
    vx = np.array([[0, -v[2], v[1]],
                    [v[2], 0, -v[0]],
                    [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx / (1 + c)


def dihedral_angle(p1, p2, p3, p4):
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    nn1 = np.linalg.norm(n1)
    nn2 = np.linalg.norm(n2)
    if nn1 < 1e-10 or nn2 < 1e-10:
        return 0.0
    n1, n2 = n1 / nn1, n2 / nn2
    x = np.dot(n1, n2)
    y = np.dot(np.cross(n1, _norm(b2)), n2)
    return np.arctan2(y, x)


def rotate_points_around_axis(points, origin, axis, angle):
    k = _norm(axis)
    ca, sa = np.cos(angle), np.sin(angle)
    centered = points - origin
    rotated = np.zeros_like(centered)
    for i in range(len(centered)):
        v = centered[i]
        rotated[i] = v * ca + np.cross(k, v) * sa + k * np.dot(k, v) * (1 - ca)
    return rotated + origin


def place_cap(cap_coords, cap_elems, connect_idx, target_pos, bond_dir):
    cap_out = outgoing_direction(cap_coords, cap_elems, connect_idx)
    R = rotation_matrix_align(cap_out, bond_dir)
    centered = cap_coords - cap_coords[connect_idx]
    rotated = (R @ centered.T).T
    return rotated + target_pos


# =============================================================================
# 组装 ACE-RES-NME
# =============================================================================

def assemble(res_e, res_c, res_n, ace_e, ace_c, nme_e, nme_c, cfg):
    """组装加帽结构, 确保肽键 O=C-N-H 共面 (trans, ω≈180°).

    ACE 的自由价(羰基C) → 残基主链 N
    残基主链 C(=O)      → NME 的自由价(酰胺N)
    """
    bb_n = cfg["backbone_n_idx"]
    bb_c = cfg["backbone_c_idx"]
    ace_ci = cfg.get("ace_connect_idx", 0)
    ace_oi = cfg.get("ace_o_idx", 1)
    nme_ci = cfg.get("nme_connect_idx", 4)
    nme_hi = cfg.get("nme_h_idx", 5)

    n_pos = res_c[bb_n]
    c_pos = res_c[bb_c]

    res_h_on_n = find_bonded_elem(res_c, res_e, bb_n, 'H')
    res_o_on_c = find_bonded_elem(res_c, res_e, bb_c, 'O')

    v_n = outgoing_direction(res_c, res_e, bb_n)
    v_c = outgoing_direction(res_c, res_e, bb_c)

    # 放置 ACE: ACE_C(羰基碳) → 残基_N
    ace_target = n_pos + PEPTIDE_CN * v_n
    ace_bond = _norm(n_pos - ace_target)
    ace_placed = place_cap(ace_c, ace_e, ace_ci, ace_target, ace_bond)

    if res_h_on_n is not None:
        cur = dihedral_angle(ace_placed[ace_oi], ace_placed[ace_ci],
                             n_pos, res_c[res_h_on_n])
        rot_angle = 0.0 - cur
        axis_dir = n_pos - ace_placed[ace_ci]
        ace_placed = rotate_points_around_axis(
            ace_placed, ace_placed[ace_ci], axis_dir, rot_angle)
        final_dihed = np.degrees(dihedral_angle(
            ace_placed[ace_oi], ace_placed[ace_ci], n_pos, res_c[res_h_on_n]))
        print(f"    ACE肽键平面: dihedral(O,C,N,H) = {final_dihed:.1f}°  (目标 0°)")

    # 放置 NME: 残基_C(羰基碳) → NME_N(酰胺氮)
    nme_target = c_pos + PEPTIDE_CN * v_c
    nme_bond = _norm(c_pos - nme_target)
    nme_placed = place_cap(nme_c, nme_e, nme_ci, nme_target, nme_bond)

    if res_o_on_c is not None:
        cur = dihedral_angle(res_c[res_o_on_c], c_pos,
                             nme_placed[nme_ci], nme_placed[nme_hi])
        rot_angle = 0.0 - cur
        axis_dir = c_pos - nme_placed[nme_ci]
        nme_placed = rotate_points_around_axis(
            nme_placed, nme_placed[nme_ci], axis_dir, rot_angle)
        final_dihed = np.degrees(dihedral_angle(
            res_c[res_o_on_c], c_pos, nme_placed[nme_ci], nme_placed[nme_hi]))
        print(f"    NME肽键平面: dihedral(O,C,N,H) = {final_dihed:.1f}°  (目标 0°)")

    # 合并: ACE + RES + NME
    n_ace = len(ace_e)
    n_res = len(res_e)
    n_nme = len(nme_e)

    all_e = list(ace_e) + list(res_e) + list(nme_e)
    all_c = np.vstack([ace_placed, res_c, nme_placed])

    # 段标签
    seg_labels = (["ACE"] * n_ace + ["RES"] * n_res + ["NME"] * n_nme)

    # 检测所有键 (包括新形成的肽键, 因为距离=1.335Å 在检测范围内)
    bonds = get_all_bonds(all_c, all_e)

    # 验证肽键是否被检测到
    ace_c_global = ace_ci                  # ACE羰基C 的全局索引
    res_n_global = n_ace + bb_n            # 残基N 的全局索引
    res_c_global = n_ace + bb_c            # 残基C 的全局索引
    nme_n_global = n_ace + n_res + nme_ci  # NME酰胺N 的全局索引

    bond_set = set((min(a, b), max(a, b)) for a, b in bonds)
    pep1 = (min(ace_c_global, res_n_global), max(ace_c_global, res_n_global))
    pep2 = (min(res_c_global, nme_n_global), max(res_c_global, nme_n_global))

    if pep1 not in bond_set:
        bonds.append(pep1)
        d = np.linalg.norm(all_c[pep1[0]] - all_c[pep1[1]])
        print(f"    手动添加肽键 ACE_C({ace_c_global})-N({res_n_global}), 距离 {d:.3f}Å")
    if pep2 not in bond_set:
        bonds.append(pep2)
        d = np.linalg.norm(all_c[pep2[0]] - all_c[pep2[1]])
        print(f"    手动添加肽键 C({res_c_global})-NME_N({nme_n_global}), 距离 {d:.3f}Å")

    d1 = np.linalg.norm(all_c[ace_c_global] - all_c[res_n_global])
    d2 = np.linalg.norm(all_c[res_c_global] - all_c[nme_n_global])
    print(f"    肽键: ACE_C-N = {d1:.3f}Å, C-NME_N = {d2:.3f}Å")

    return {
        "elements": all_e,
        "coords": all_c,
        "bonds": bonds,
        "seg_labels": seg_labels,
        "n_ace": n_ace,
        "n_res": n_res,
        "n_nme": n_nme,
        "bb_n_in_res": bb_n,
        "bb_c_in_res": bb_c,
    }


# =============================================================================
# RDKit MMFF 优化
# =============================================================================

def optimize_mmff(xyz_path, out_path, charge=0, max_iters=500):
    """用 RDKit MMFF94 力场优化结构.

    读取 xyz, 优化后写回 xyz.
    返回 True 表示成功.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdDetermineBonds
    except ImportError:
        print("  警告: 未安装 rdkit, 跳过 MMFF 优化")
        return False

    raw_mol = Chem.MolFromXYZFile(xyz_path)
    if raw_mol is None:
        print("  警告: RDKit 无法读取 xyz, 跳过 MMFF 优化")
        return False

    try:
        rdDetermineBonds.DetermineBonds(raw_mol, charge=charge)
    except Exception as e:
        print(f"  警告: RDKit 键检测失败 ({e}), 跳过 MMFF 优化")
        return False

    mol = Chem.RWMol(raw_mol)

    # MMFF 优化
    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
    if mp is None:
        print("  警告: MMFF94 参数化失败, 跳过")
        return False

    ff = AllChem.MMFFGetMoleculeForceField(mol, mp)
    if ff is None:
        print("  警告: MMFF94 力场构建失败, 跳过")
        return False

    e_before = ff.CalcEnergy()
    ret = ff.Minimize(maxIts=max_iters)
    e_after = ff.CalcEnergy()
    print(f"  MMFF94: E = {e_before:.1f} → {e_after:.1f} kcal/mol "
          f"(收敛={'是' if ret == 0 else '否'})")

    # 提取优化后坐标, 写 xyz
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    elems = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(n)]
    coords = np.array([[conf.GetAtomPosition(i).x,
                         conf.GetAtomPosition(i).y,
                         conf.GetAtomPosition(i).z] for i in range(n)])
    write_xyz(out_path, elems, coords, "MMFF94 optimized")
    return True


# =============================================================================
# Gaussian
# =============================================================================

def write_gaussian_input(path, elems, coords, cfg):
    chk = os.path.basename(path).replace('.com', '.chk')
    rn = cfg["residue_name"]
    opt_method = cfg.get("opt_method", "wB97XD/6-31G*")
    resp_method = cfg.get("resp_method", "HF/6-31G*")
    chg = cfg["net_charge"]
    mult = cfg.get("multiplicity", 1)

    with open(path, 'w') as f:
        # Job 1: 结构优化
        f.write(f"%chk={chk}\n")
        f.write(f"%nproc={cfg.get('nproc', 16)}\n")
        f.write(f"%mem={cfg.get('mem', '64GB')}\n")
        f.write(f"# Opt {opt_method} SCF=Tight\n")
        f.write(f"\nACE-{rn}-NME optimization\n\n")
        f.write(f"{chg} {mult}\n")
        for e, (x, y, z) in zip(elems, coords):
            f.write(f" {e:<2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")
        f.write("\n")

        # --Link1-- 单点 ESP
        f.write("--Link1--\n")
        f.write(f"%chk={chk}\n")
        f.write(f"%nproc={cfg.get('nproc', 16)}\n")
        f.write(f"%mem={cfg.get('mem', '64GB')}\n")
        f.write(f"# {resp_method} Geom=AllCheck Guess=Read "
                f"Pop=MK IOp(6/33=2) SCF=Tight\n")
        f.write("\n")


def run_gaussian(com_file, cfg):
    log_file = os.path.splitext(com_file)[0] + ".log"
    cmd = f"{cfg.get('gaussian_cmd', 'g16')} {os.path.basename(com_file)}"
    cwd = os.path.dirname(os.path.abspath(com_file))
    print(f"  [RUN] cd {cwd} && {cmd}")
    sp.Popen(cmd, shell=True, cwd=cwd)
    return log_file


def monitor_gaussian(log_file, interval=60):
    print(f"  监控: {log_file}  (每 {interval}s)")
    out_file = log_file.replace('.log', '.out')
    start = time.time()

    while True:
        time.sleep(interval)
        elapsed = time.time() - start
        mins = int(elapsed // 60)

        target = None
        for f in [log_file, out_file]:
            if os.path.exists(f) and os.path.getsize(f) > 0:
                target = f
                break

        if target is None:
            print(f"    [{mins}min] 等待输出文件 ...")
            continue

        content = open(target).read()
        if "Normal termination" in content:
            count = content.count("Normal termination")
            if count >= 2:
                print(f"    [{mins}min] Gaussian 正常结束! (Opt + ESP 均完成)")
                return True, target
            else:
                print(f"    [{mins}min] 优化完成, ESP 计算中 ...")
                continue
        if "Error termination" in content:
            print(f"    [{mins}min] ERROR: Gaussian 异常终止!")
            return False, target

        nsteps = content.count("Step number")
        print(f"    [{mins}min] 运行中, 优化步数 ~{nsteps}")


# =============================================================================
# AmberTools
# =============================================================================

def _run_cmd(cmd, cwd=None):
    print(f"  [CMD] {cmd}")
    r = sp.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
    if r.returncode != 0:
        print(f"    STDERR: {r.stderr[:800]}")
    return r.returncode, r.stdout, r.stderr


def read_ac_names(ac_path):
    names = []
    with open(ac_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                names.append(line.split()[2])
    return names


def write_residue_pdb(ac_path, out_pdb, n_ace, n_res, rn, cfg):
    """从 AC 文件提取残基原子, 用原始坐标输出 PDB (原子名与 prepin 一致)."""
    # 读取残基输入 PDB 的坐标 (未经 Gaussian 优化的原始坐标)
    from autoprep.normalize import _parse_coord
    res_file = cfg["residue_file"]
    with open(res_file) as f:
        res_lines = [l for l in f if l[:6].strip() in ("ATOM", "HETATM")]
    res_coords = [_parse_coord(l) for l in res_lines]

    # 读取 AC 文件中残基部分的原子名
    ac_atoms = []
    with open(ac_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                name = line.split()[2]
                # 元素从原子名提取: 取前缀字母部分
                elem = ''.join(c for c in name if c.isalpha())
                if len(elem) > 2:
                    elem = elem[:1]
                ac_atoms.append((name, elem.capitalize()))
    res_ac = ac_atoms[n_ace:n_ace + n_res]

    if len(res_ac) != len(res_coords):
        print(f"  [WARNING] AC 残基原子数 {len(res_ac)} != 输入 PDB {len(res_coords)}")
        return

    with open(out_pdb, 'w') as f:
        for i, ((name, elem), (x, y, z)) in enumerate(zip(res_ac, res_coords)):
            serial = i + 1
            if len(name) <= 3:
                name_field = f" {name:<3s}"
            else:
                name_field = f"{name:<4s}"
            f.write(f"ATOM  {serial:5d} {name_field}"
                    f" {rn:>3s} A   1    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}"
                    f"{1.0:6.2f}{0.0:6.2f}"
                    f"          {elem:>2s}  \n")
        f.write("END\n")


def run_antechamber(gout_file, cfg, workdir):
    rn = cfg["residue_name"]
    ac = f"{rn}.ac"
    cmd = (f"{cfg.get('antechamber', 'antechamber')} "
           f"-i {os.path.basename(gout_file)} -fi gout "
           f"-o {ac} -fo ac "
           f"-c resp -s 2 "
           f"-rn {rn} -at gaff2 "
           f"-nc {cfg['net_charge']} -pf y")
    rc, _, err = _run_cmd(cmd, cwd=workdir)
    if rc != 0:
        raise RuntimeError(f"antechamber 失败:\n{err}")
    return os.path.join(workdir, ac)


def write_mainchain_mc(path, ace_names, nme_names,
                       head, tail, mainchain, charge):
    with open(path, 'w') as f:
        f.write(f"HEAD_NAME {head}\n")
        f.write(f"TAIL_NAME {tail}\n")
        for mc in mainchain:
            f.write(f"MAIN_CHAIN {mc}\n")
        for nm in ace_names:
            f.write(f"OMIT_NAME {nm}\n")
        for nm in nme_names:
            f.write(f"OMIT_NAME {nm}\n")
        f.write("PRE_HEAD_TYPE C\n")
        f.write("POST_TAIL_TYPE N\n")
        f.write(f"CHARGE {charge:.1f}\n")


def run_prepgen(ac_file, mc_file, cfg, workdir):
    rn = cfg["residue_name"]
    out = f"{rn}.prepin"
    cmd = (f"{cfg.get('prepgen', 'prepgen')} "
           f"-i {os.path.basename(ac_file)} "
           f"-o {out} "
           f"-m {os.path.basename(mc_file)} "
           f"-rn {rn}")
    rc, _, err = _run_cmd(cmd, cwd=workdir)
    if rc != 0:
        raise RuntimeError(f"prepgen 失败:\n{err}")
    return os.path.join(workdir, out)


def run_parmchk2(prepin_file, cfg, workdir):
    rn = cfg["residue_name"]
    out = f"{rn}.frcmod"
    cmd = (f"{cfg.get('parmchk2', 'parmchk2')} "
           f"-i {os.path.basename(prepin_file)} -f prepi "
           f"-o {out} -s gaff2")
    rc, _, err = _run_cmd(cmd, cwd=workdir)
    if rc != 0:
        raise RuntimeError(f"parmchk2 失败:\n{err}")
    return os.path.join(workdir, out)


# =============================================================================
# prepin 电荷校验
# =============================================================================

def verify_prepin_charge(prepin_path, expected_charge):
    """校验 prepin 文件中残基总电荷是否为整数且等于期望值."""
    total = 0.0
    n_atoms = 0
    has_nan = False
    with open(prepin_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 11:
                continue
            # 跳过 DUMM 行
            if parts[1] == 'DUMM':
                continue
            try:
                int(parts[0])
                charge = float(parts[10])
                total += charge
                n_atoms += 1
            except (ValueError, IndexError):
                continue
            # 检测 nan
            if 'nan' in line.lower():
                has_nan = True

    ok = True
    if has_nan:
        print(f"  [WARNING] prepin 含有 nan, prepgen 内坐标树构建可能有问题!")
        ok = False

    rounded = round(total)
    if abs(total - rounded) > 0.01:
        print(f"  [WARNING] prepin 总电荷 {total:.6f} 不是整数! ({n_atoms} atoms)")
        ok = False
    if rounded != expected_charge:
        print(f"  [WARNING] prepin 电荷 {rounded} != 配置 residue_charge {expected_charge}")
        ok = False
    if ok:
        print(f"  电荷校验通过: {total:.4f} ≈ {expected_charge} ({n_atoms} atoms)")
    return ok


def fix_prepin_charge(prepin_path, expected_charge):
    """修正 prepin 电荷: 将误差均匀分配到主链 (M 标记) 原子上.

    仅在总电荷与期望值差异 > 1e-4 时修正.
    """
    with open(prepin_path) as f:
        lines = f.readlines()

    # 第一遍: 计算总电荷, 找主链原子行号
    total = 0.0
    backbone_indices = []   # line indices of M-flagged atoms
    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) < 11:
            continue
        if parts[1] == 'DUMM':
            continue
        try:
            int(parts[0])
            charge = float(parts[10])
        except (ValueError, IndexError):
            continue
        total += charge
        # parts[3] 是 tree 标志: M=mainchain, S=sidechain, B=branch, E=end
        if parts[3] == 'M':
            backbone_indices.append(i)

    error = expected_charge - total
    if abs(error) < 1e-4:
        return False   # 无需修正

    if not backbone_indices:
        print(f"  [WARNING] 无主链原子 (M 标记), 无法分配电荷修正 {error:.6f}")
        return False

    correction = error / len(backbone_indices)
    print(f"  电荷修正: {total:.6f} → {expected_charge} "
          f"(误差 {error:+.6f}, 分配到 {len(backbone_indices)} 个主链原子, "
          f"每个 {correction:+.6f})")

    # 第二遍: 修正主链原子电荷
    for i in backbone_indices:
        line = lines[i]
        parts = line.split()
        old_charge = float(parts[10])
        new_charge = old_charge + correction

        # prepin 电荷字段在最后一列, 找到其起始位置并替换
        # 格式: 最后一个字段是电荷, 宽度一般为 10 字符 (%10.6f)
        last_space = line.rstrip('\n').rfind(' ')
        lines[i] = line[:last_space + 1] + f"{new_charge:.6f}\n"

    with open(prepin_path, 'w') as f:
        f.writelines(lines)

    return True


# =============================================================================
# ff14SB / gaff2 交叉参数库
# =============================================================================
# 非标准氨基酸 (gaff2) 接入蛋白 (ff14SB) 时, 肽键连接处产生大小写
# 混合的原子类型参数. 此库自动补充缺失的交叉项.

def _bond_key(line):
    if len(line) < 5 or line[2] != '-':
        return None
    t = (line[0:2].strip(), line[3:5].strip())
    if not t[0] or not t[1]:
        return None
    return min(t, t[::-1])

def _angle_key(line):
    if len(line) < 8 or line[2] != '-' or line[5] != '-':
        return None
    t = (line[0:2].strip(), line[3:5].strip(), line[6:8].strip())
    if not all(t):
        return None
    return min(t, t[::-1])

def _dihe_key(line):
    if len(line) < 11 or line[2] != '-' or line[5] != '-' or line[8] != '-':
        return None
    t = (line[0:2].strip(), line[3:5].strip(), line[6:8].strip(), line[9:11].strip())
    if not all(t):
        return None
    return min(t, t[::-1])


# --- 交叉项参数 ---
# HEAD 连接: ff14SB C(=O) → gaff2 ns   (蛋白前一残基 → 非标准残基 N端)
# TAIL 连接: gaff2 c(=O) → ff14SB N(-H) (非标准残基 C端 → 蛋白后一残基)

_CROSSTERM_BOND = [
    "C -ns  490.000   1.335       ff14SB/gaff2 peptide bond\n",
    "c -N   490.000   1.335       ff14SB/gaff2 peptide bond\n",
]

_CROSSTERM_ANGLE = [
    # HEAD: protein C(=O) → residue ns
    "O -C -ns    80.000     122.900   ff14SB(O,C)/gaff2(ns)\n",
    "C -ns-hn    50.000     120.000   ff14SB(C)/gaff2(ns,hn)\n",
    "C -ns-c3    50.000     121.900   ff14SB(C)/gaff2(ns,c3)\n",
    "CX-C -ns    70.000     116.600   ff14SB(CX,C)/gaff2(ns)\n",
    # TAIL: residue c(=O) → protein N(-H)
    "o -c -N     80.000     122.900   gaff2(o,c)/ff14SB(N)\n",
    "c -N -H     80.000     122.900   gaff2(c)/ff14SB(N,H)\n",
    "c -N -CX    50.000     121.900   gaff2(c)/ff14SB(N,CX)\n",
    "c3-c -N     70.000     116.600   gaff2(c3,c)/ff14SB(N)\n",
]

_CROSSTERM_DIHE = [
    # HEAD: central bond C-ns (ff14SB C=O to gaff2 amide N)
    "O -C -ns-hn   1    2.500       180.000          -2.000      ff14SB/gaff2\n",
    "O -C -ns-hn   1    2.000         0.000           1.000      ff14SB/gaff2\n",
    "O -C -ns-c3   4   10.000       180.000           2.000      ff14SB/gaff2\n",
    "CX-C -ns-c3   4   10.000       180.000           2.000      ff14SB/gaff2\n",
    "CX-C -ns-hn   4   10.000       180.000           2.000      ff14SB/gaff2\n",
    # HEAD: central bond ns-c3
    "C -ns-c3-c    6    0.000         0.000           2.000      ff14SB/gaff2\n",
    "C -ns-c3-c3   6    0.000         0.000           2.000      ff14SB/gaff2\n",
    "C -ns-c3-h1   6    0.000         0.000           2.000      ff14SB/gaff2\n",
    # TAIL: central bond c-N (gaff2 C=O to ff14SB amide N)
    "o -c -N -H    1    2.500       180.000          -2.000      ff14SB/gaff2\n",
    "o -c -N -H    1    2.000         0.000           1.000      ff14SB/gaff2\n",
    "o -c -N -CX   4   10.000       180.000           2.000      ff14SB/gaff2\n",
    "c3-c -N -H    4   10.000       180.000           2.000      ff14SB/gaff2\n",
    "c3-c -N -CX   4   10.000       180.000           2.000      ff14SB/gaff2\n",
    # TAIL: central bond N-CX
    "c -N -CX-C    4   10.000       180.000           2.000      ff14SB/gaff2\n",
    "c -N -CX-H1   4   10.000       180.000           2.000      ff14SB/gaff2\n",
]


def patch_frcmod_crossterms(frcmod_path):
    """补充 ff14SB/gaff2 交叉参数到 frcmod."""
    with open(frcmod_path) as f:
        lines = f.readlines()

    SECTIONS = {'MASS', 'BOND', 'ANGLE', 'ANGL', 'DIHE', 'IMPROPER', 'IMPR', 'NONBON'}
    current = None
    existing = {'BOND': set(), 'ANGLE': set(), 'DIHE': set()}
    section_end = {}   # section_name -> line index of terminating blank line

    for i, line in enumerate(lines):
        s = line.strip()
        if s in SECTIONS:
            current = {'ANGL': 'ANGLE', 'IMPR': 'IMPROPER'}.get(s, s)
        elif s == '' and current:
            section_end[current] = i
            current = None
        elif current in existing:
            fn = {'BOND': _bond_key, 'ANGLE': _angle_key, 'DIHE': _dihe_key}[current]
            k = fn(line)
            if k:
                existing[current].add(k)

    # 缺失的 BOND / ANGLE
    missing = {'BOND': [], 'ANGLE': [], 'DIHE': []}

    for line in _CROSSTERM_BOND:
        k = _bond_key(line)
        if k and k not in existing['BOND']:
            missing['BOND'].append(line)
            existing['BOND'].add(k)

    for line in _CROSSTERM_ANGLE:
        k = _angle_key(line)
        if k and k not in existing['ANGLE']:
            missing['ANGLE'].append(line)
            existing['ANGLE'].add(k)

    # 缺失的 DIHE (支持多项式: 同一 key 可能有多行)
    dihe_groups = {}
    for line in _CROSSTERM_DIHE:
        k = _dihe_key(line)
        if k:
            dihe_groups.setdefault(k, []).append(line)
    for k, grp in dihe_groups.items():
        if k not in existing['DIHE']:
            missing['DIHE'].extend(grp)
            existing['DIHE'].add(k)

    total = sum(len(v) for v in missing.values())
    if total == 0:
        return 0

    # 在各 section 的结束空行前插入交叉项
    output = []
    for i, line in enumerate(lines):
        for sec in ('BOND', 'ANGLE', 'DIHE'):
            if i == section_end.get(sec) and missing[sec]:
                for ml in missing[sec]:
                    output.append(ml)
        output.append(line)

    with open(frcmod_path, 'w') as f:
        f.writelines(output)
    return total


# =============================================================================
# 自动检测
# =============================================================================

def auto_find_ca(coords, elems, bb_n, bb_c):
    bn = set(find_bonded(coords, elems, bb_n))
    bc = set(find_bonded(coords, elems, bb_c))
    candidates = [i for i in bn & bc if elems[i] == 'C']
    if candidates:
        return candidates[0]
    for i in bn:
        if elems[i] == 'C' and bb_c in set(find_bonded(coords, elems, i)):
            return i
    return None


def print_atoms(elems, coords, names, label=""):
    if label:
        print(f"\n  --- {label} ---")
    print(f"  {'Idx':>4s}  {'Name':>5s}  {'Elem':>4s}  "
          f"{'X':>10s}  {'Y':>10s}  {'Z':>10s}")
    for i, (e, (x, y, z), nm) in enumerate(zip(elems, coords, names)):
        print(f"  {i:4d}  {nm:>5s}  {e:>4s}  "
              f"{x:10.4f}  {y:10.4f}  {z:10.4f}")


# =============================================================================
# 检查点
# =============================================================================

def _gaussian_done(log_file):
    """检查 Gaussian log 是否有两次 Normal termination."""
    out_file = log_file.replace('.log', '.out')
    for f in [log_file, out_file]:
        if os.path.exists(f) and os.path.getsize(f) > 0:
            content = open(f).read()
            if content.count("Normal termination") >= 2:
                return True, f
    return False, log_file


def _save_meta(workdir, meta):
    with open(os.path.join(workdir, ".autoprep_meta.json"), 'w') as f:
        json.dump(meta, f, indent=2)


def _load_meta(workdir):
    p = os.path.join(workdir, ".autoprep_meta.json")
    if os.path.exists(p):
        with open(p) as f:
            return json.load(f)
    return None


# =============================================================================
# 主流程
# =============================================================================

def run_pipeline(cfg, workdir):
    """执行完整参数化流程, 支持断点续跑."""
    rn = cfg["residue_name"]
    os.makedirs(workdir, exist_ok=True)

    opt_method = cfg.get("opt_method", "wB97XD/6-31G*")
    resp_method = cfg.get("resp_method", "HF/6-31G*")

    print("=" * 60)
    print("  autoprep - 非标准氨基酸自动参数化")
    print(f"  残基: {rn}   电荷: {cfg['net_charge']}")
    print(f"  优化: {opt_method}   RESP: {resp_method}")
    print(f"  工作目录: {workdir}")
    print("=" * 60)

    # 关键文件路径
    pdb_out = os.path.join(workdir, f"{rn}_capped.pdb")
    xyz_raw = os.path.join(workdir, f"{rn}_capped.xyz")
    xyz_mmff = os.path.join(workdir, f"{rn}_mmff.xyz")
    com = os.path.join(workdir, f"{rn}.com")
    log = os.path.join(workdir, f"{rn}.log")
    ac_file = os.path.join(workdir, f"{rn}.ac")
    mc_file = os.path.join(workdir, f"{rn}.mc")
    prepin = os.path.join(workdir, f"{rn}.prepin")
    frcmod = os.path.join(workdir, f"{rn}.frcmod")

    # ==========================================================
    # Step 1-2: 读取输入 + 组装
    # ==========================================================
    meta = _load_meta(workdir)
    if meta and os.path.exists(pdb_out):
        print(f"\n[Step 1-2] 跳过 (已有 {rn}_capped.pdb + 元数据)")
        na = meta["n_ace"]
        nr = meta["n_res"]
        nn = meta["n_nme"]
        res_file = cfg["residue_file"]
        if not os.path.isabs(res_file):
            res_file = os.path.abspath(res_file)
        if cfg.get("precapped", False):
            # In precapped mode, residue_file IS the full ACE-RES-NME PDB
            res_e, res_c, res_n = read_structure(res_file)
            # Slice to RES portion only for downstream code that needs res_e
            res_e = res_e[na:na+nr]
            res_c = res_c[na:na+nr]
            res_n = res_n[na:na+nr] if res_n else None
        else:
            res_e, res_c, res_n = read_structure(res_file)
    elif cfg.get("precapped", False):
        # ★ Precapped mode: input PDB already contains ACE+RES+NME with proper labels
        print("\n[Step 1-2] PRECAPPED 模式 — 跳过 ACE/NME 自动加帽")
        res_file = cfg["residue_file"]
        if not os.path.isabs(res_file):
            res_file = os.path.abspath(res_file)

        all_e, all_c, all_n = read_structure(res_file)

        # Identify segment by residue name (need to re-read with PDB parser)
        seg_labels = []
        atom_names = []
        with open(res_file) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    resname = line[17:20].strip().upper()
                    aname = line[12:16].strip()
                    seg_labels.append(resname)
                    atom_names.append(aname)
        if len(seg_labels) != len(all_e):
            raise RuntimeError(
                f"precapped: parsed {len(seg_labels)} atoms from PDB but "
                f"read_structure gave {len(all_e)} atoms")

        # Group consecutive ACE / RES / NME segments
        # Expected order: ACE first, RES middle, NME last
        ace_idx = [i for i, s in enumerate(seg_labels) if s == "ACE"]
        nme_idx = [i for i, s in enumerate(seg_labels) if s == "NME"]
        res_idx = [i for i, s in enumerate(seg_labels) if s not in ("ACE", "NME")]
        if not ace_idx or not res_idx or not nme_idx:
            raise RuntimeError(
                f"precapped: PDB must contain ACE, RES (any name), NME residues. "
                f"Got: {set(seg_labels)}")
        if not (max(ace_idx) < min(res_idx) and max(res_idx) < min(nme_idx)):
            raise RuntimeError(
                "precapped: atoms must be ordered ACE → RES → NME (group together).")

        na, nr, nn = len(ace_idx), len(res_idx), len(nme_idx)

        # Identify backbone N and C in RES section by atom name
        bb_n_atom = cfg.get("backbone_n_atom_name", "N1")
        bb_c_atom = cfg.get("backbone_c_atom_name", "C4")
        res_atom_names = [atom_names[i] for i in res_idx]
        try:
            bb_n_in_res = res_atom_names.index(bb_n_atom)
            bb_c_in_res = res_atom_names.index(bb_c_atom)
        except ValueError:
            raise RuntimeError(
                f"precapped: backbone_n_atom_name='{bb_n_atom}' or "
                f"backbone_c_atom_name='{bb_c_atom}' not found in RES atoms: "
                f"{res_atom_names}")

        print(f"  解析: ACE={na}  RES={nr}  NME={nn}  总={na+nr+nn} 原子")
        print(f"  backbone N atom: '{bb_n_atom}' (RES idx {bb_n_in_res})")
        print(f"  backbone C atom: '{bb_c_atom}' (RES idx {bb_c_in_res})")

        # Copy input PDB to capped.pdb (autoprep expected name) + write xyz
        import shutil
        shutil.copy2(res_file, pdb_out)
        write_xyz(xyz_raw, all_e, all_c, f"ACE-{rn}-NME (precapped)")
        print(f"  保存: {pdb_out}")
        print(f"  保存: {xyz_raw}")

        _save_meta(workdir, {
            "n_ace": na, "n_res": nr, "n_nme": nn,
            "bb_n_in_res": bb_n_in_res,
            "bb_c_in_res": bb_c_in_res,
        })
        meta = _load_meta(workdir)

        # Slice res_* for downstream
        res_e = [all_e[i] for i in res_idx]
        res_c = [all_c[i] for i in res_idx]
        res_n = [all_n[i] for i in res_idx] if all_n else None
    else:
        print("\n[Step 1] 读取输入 ...")
        res_file = cfg["residue_file"]
        if not os.path.isabs(res_file):
            res_file = os.path.abspath(res_file)
        res_e, res_c, res_n = read_structure(res_file)
        ace_e, ace_c = read_xyz(str(ACE_TEMPLATE))
        nme_e, nme_c = read_xyz(str(NME_TEMPLATE))
        print(f"  残基: {len(res_e)}  ACE: {len(ace_e)}  NME: {len(nme_e)}  原子")
        print_atoms(res_e, res_c, res_n, "残基原子 (请确认 backbone_n/c 索引)")

        print(f"\n[Step 2] 组装 ACE-{rn}-NME ...")
        asm = assemble(res_e, res_c, res_n, ace_e, ace_c, nme_e, nme_c, cfg)
        total = len(asm["elements"])
        na, nr, nn = asm["n_ace"], asm["n_res"], asm["n_nme"]
        print(f"  总原子: {total}  (ACE:{na} + RES:{nr} + NME:{nn})")

        # 输出 PDB (含 CONECT) + XYZ
        write_assembled_pdb(pdb_out, asm["elements"], asm["coords"],
                            asm["bonds"], rn=rn, seg_labels=asm["seg_labels"])
        write_xyz(xyz_raw, asm["elements"], asm["coords"],
                  f"ACE-{rn}-NME capped structure")
        print(f"  保存: {pdb_out}  ← 用分子可视化软件检查 CONECT!")
        print(f"  保存: {xyz_raw}")

        _save_meta(workdir, {
            "n_ace": na, "n_res": nr, "n_nme": nn,
            "bb_n_in_res": asm["bb_n_in_res"],
            "bb_c_in_res": asm["bb_c_in_res"],
        })
        meta = _load_meta(workdir)

    bb_n_in_res = meta["bb_n_in_res"]
    bb_c_in_res = meta["bb_c_in_res"]

    # ==========================================================
    # Step 3: RDKit MMFF94 预优化
    # ==========================================================
    if os.path.exists(xyz_mmff):
        print(f"\n[Step 3] 跳过 (已有 {rn}_mmff.xyz)")
    else:
        print(f"\n[Step 3] RDKit MMFF94 预优化 ...")
        ok = optimize_mmff(xyz_raw, xyz_mmff, charge=cfg["net_charge"])
        if not ok:
            # MMFF 失败, 用原始坐标继续
            import shutil
            shutil.copy2(xyz_raw, xyz_mmff)
            print("  使用原始组装坐标继续")

    # ==========================================================
    # Step 4: Gaussian 输入
    # ==========================================================
    if os.path.exists(com):
        print(f"\n[Step 4] 跳过 (已有 {rn}.com)")
    else:
        print(f"\n[Step 4] Gaussian 输入 (Opt + --Link1-- SP ESP) ...")
        elems, coords = read_xyz(xyz_mmff)
        write_gaussian_input(com, elems, coords, cfg)
        print(f"  文件: {com}")
        print(f"  Job1: Opt {opt_method}")
        print(f"  Job2: SP  {resp_method} Pop=MK IOp(6/33=2)")
        print(f"  nproc={cfg.get('nproc', 16)}  mem={cfg.get('mem', '64GB')}")

    # ==========================================================
    # Step 5: 运行 Gaussian
    # ==========================================================
    done, actual_log = _gaussian_done(log)
    if done:
        print(f"\n[Step 5] 跳过 (Gaussian 已正常结束: {os.path.basename(actual_log)})")
    else:
        print(f"\n[Step 5] 运行 Gaussian ...")
        if os.path.exists(log) and os.path.getsize(log) > 0:
            content = open(log).read()
            if "Error termination" in content:
                print(f"  检测到上次 Error termination, 删除旧 log 重新提交 ...")
                os.remove(log)
                chk = log.replace('.log', '.chk')
                if os.path.exists(chk):
                    os.remove(chk)
        actual_log = run_gaussian(com, cfg)
        ok, actual_log = monitor_gaussian(actual_log, cfg.get("check_interval", 60))
        if not ok:
            print("Gaussian 失败, 请检查后重新运行.")
            sys.exit(1)

    # ==========================================================
    # Step 6: antechamber RESP
    # ==========================================================
    if os.path.exists(ac_file):
        print(f"\n[Step 6] 跳过 (已有 {rn}.ac)")
    else:
        print(f"\n[Step 6] antechamber RESP 拟合 ...")
        if os.path.dirname(os.path.abspath(actual_log)) != workdir:
            import shutil
            shutil.copy2(actual_log, workdir)
            actual_log = os.path.join(workdir, os.path.basename(actual_log))
        run_antechamber(actual_log, cfg, workdir)
        print(f"  .ac: {ac_file}")

    ac_names = read_ac_names(ac_file)
    ace_ac = ac_names[:na]
    res_ac = ac_names[na:na + nr]
    nme_ac = ac_names[na + nr:]

    head = res_ac[bb_n_in_res]
    tail = res_ac[bb_c_in_res]
    print(f"  HEAD={head} (主链N)  TAIL={tail} (主链C)")

    ca_idx = auto_find_ca(res_c, res_e, bb_n_in_res, bb_c_in_res)
    # MAIN_CHAIN 只列 HEAD 和 TAIL 之间的中间原子
    # prepgen 会自动把 HEAD_NAME/TAIL_NAME 加入 mainchain
    mc = []
    if ca_idx is not None:
        mc.append(res_ac[ca_idx])
        print(f"  自动检测 CA: {res_ac[ca_idx]} (残基内索引 {ca_idx})")

    if cfg.get("mainchain_names"):
        mc = cfg["mainchain_names"]
        print(f"  使用用户指定 MAIN_CHAIN: {mc}")

    # ==========================================================
    # Step 7: prepgen
    # ==========================================================
    if os.path.exists(prepin):
        print(f"\n[Step 7] 跳过 (已有 {rn}.prepin)")
    else:
        print(f"\n[Step 7] prepgen → .prepin ...")
        charge = cfg.get("residue_charge", cfg["net_charge"])
        write_mainchain_mc(mc_file, ace_ac, nme_ac, head, tail, mc, charge)
        print(f"  mainchain: {mc_file}")
        run_prepgen(ac_file, mc_file, cfg, workdir)
        print(f"  生成: {prepin}")

    # 校验并修正 prepin 电荷
    expected_q = cfg.get("residue_charge", cfg["net_charge"])
    charge_ok = verify_prepin_charge(prepin, expected_q)
    if not charge_ok:
        fixed = fix_prepin_charge(prepin, expected_q)
        if fixed:
            verify_prepin_charge(prepin, expected_q)

    # ==========================================================
    # Step 8: parmchk2
    # ==========================================================
    if os.path.exists(frcmod):
        print(f"\n[Step 8] 跳过 (已有 {rn}.frcmod)")
    else:
        print(f"\n[Step 8] parmchk2 → .frcmod ...")
        run_parmchk2(prepin, cfg, workdir)
        print(f"  生成: {frcmod}")

    # 补充 ff14SB/gaff2 交叉参数
    n_cross = patch_frcmod_crossterms(frcmod)
    if n_cross:
        print(f"  添加 {n_cross} 条 ff14SB/gaff2 交叉参数")

    # ==========================================================
    # Step 9: 输出残基 PDB (原子名与 prepin 一致)
    # ==========================================================
    res_pdb = os.path.join(workdir, f"{rn}.pdb")
    write_residue_pdb(ac_file, res_pdb, na, nr, rn, cfg)
    print(f"\n[Step 9] 残基 PDB → {res_pdb}")
    print(f"  原子名与 prepin 一致, 可直接替换蛋白 PDB 中对应残基")

    # ==========================================================
    # 完成
    # ==========================================================
    print("\n" + "=" * 60)
    print("  完成! 输出文件:")
    print(f"    {prepin}")
    print(f"    {frcmod}")
    print(f"    {res_pdb}")
    print(f"\n  tleap 加载示例:")
    print(f"    source leaprc.protein.ff14SB")
    print(f"    source leaprc.gaff2")
    print(f"    loadamberprep {rn}.prepin")
    print(f"    loadamberparams {rn}.frcmod")
    print(f"    x = loadpdb your_protein.pdb")
    print(f"    check x")
    print("=" * 60)
