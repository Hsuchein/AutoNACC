"""
autoprep.align - 将蛋白 PDB 中非标准残基的原子名对齐到参考 PDB 命名

输入: 蛋白 PDB + autoprep run 输出的参考残基 .pdb 文件
从参考 PDB 读取残基名, 在蛋白中查找同名残基, 按元素顺序匹配替换原子名.
去质子化多余的 H 自动删除.
"""


def _read_ref_pdb(pdb_path):
    """读取参考残基 PDB, 返回 (resname, [(name, elem), ...])."""
    atoms = []
    resname = None
    with open(pdb_path) as f:
        for line in f:
            if line[:6].strip() not in ("ATOM", "HETATM") or len(line) < 54:
                continue
            name = line[12:16].strip()
            elem = line[76:78].strip() if len(line) >= 78 else name[0]
            if resname is None:
                resname = line[17:20].strip()
            atoms.append((name, elem))
    return resname, atoms


def _build_elem_names(atoms):
    """按元素分组, 保持出现顺序. 返回 {elem: [name1, name2, ...]}."""
    d = {}
    for name, elem in atoms:
        d.setdefault(elem, []).append(name)
    return d


def align_pdb(input_pdb, ref_pdbs, output_pdb):
    """根据参考残基 PDB 重命名蛋白中对应残基的原子名."""
    # 读取参考
    refs = {}  # resname -> {elem: [name, ...]}
    for ref in ref_pdbs:
        resname, ref_atoms = _read_ref_pdb(ref)
        if not resname:
            print(f"  [WARNING] 无法读取残基名: {ref}")
            continue
        refs[resname] = _build_elem_names(ref_atoms)
        print(f"  参考: {resname} ({len(ref_atoms)} 原子) ← {ref}")

    if not refs:
        print("  无可用参考, 退出")
        return

    # 读取蛋白 PDB
    with open(input_pdb) as f:
        lines = f.readlines()

    # 按 (resnum, chain, resname) 分组
    residue_groups = {}
    for i, line in enumerate(lines):
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM") or len(line) < 54:
            continue
        resname = line[17:20].strip()
        if resname not in refs:
            continue
        try:
            resnum = int(line[22:26])
        except ValueError:
            continue
        chain = line[21]
        name = line[12:16].strip()
        elem = line[76:78].strip() if len(line) >= 78 else name[0]
        key = (resnum, chain, resname)
        residue_groups.setdefault(key, []).append((i, name, elem))

    if not residue_groups:
        print("  未找到需要对齐的残基")
        return

    # 对每个残基做元素顺序匹配
    lines_to_rename = {}
    lines_to_remove = set()
    n_aligned = 0

    for key in sorted(residue_groups):
        resnum, chain, resname = key
        ref_elem = refs[resname]
        prot_atoms = residue_groups[key]

        prot_elem = {}
        for idx, name, elem in prot_atoms:
            prot_elem.setdefault(elem, []).append((idx, name))

        removed = []
        for elem, prot_list in prot_elem.items():
            ref_names = ref_elem.get(elem, [])
            for j, (idx, old_name) in enumerate(prot_list):
                if j < len(ref_names):
                    lines_to_rename[idx] = ref_names[j]
                else:
                    lines_to_remove.add(idx)
                    removed.append(old_name)

        n_matched = len([a for a in prot_atoms if a[0] not in lines_to_remove])
        info = f"  {resname} {resnum:4d} {chain}: {n_matched} 原子匹配"
        if removed:
            info += f", 删除 {len(removed)} ({', '.join(removed)})"
        print(info)
        n_aligned += 1

    # 输出
    out = []
    serial = 0
    for i, line in enumerate(lines):
        if i in lines_to_remove:
            continue
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM"):
            serial += 1
            if i in lines_to_rename:
                new_name = lines_to_rename[i]
                if len(new_name) <= 3:
                    name_field = f" {new_name:<3s}"
                else:
                    name_field = f"{new_name:<4s}"
                line = f"{line[:6]}{serial:5d} {name_field}{line[16:]}"
            else:
                line = f"{line[:6]}{serial:5d}{line[11:]}"
        out.append(line)

    with open(output_pdb, "w") as f:
        f.writelines(out)

    print(f"\n  对齐 {n_aligned} 个残基, "
          f"重命名 {len(lines_to_rename)} 原子, "
          f"删除 {len(lines_to_remove)} 原子")
    print(f"  保存: {output_pdb}")
