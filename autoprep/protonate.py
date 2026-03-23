"""
autoprep.protonate - 根据 propka3 pKa 预测修改 PDB 残基名以分配质子化状态

逻辑:
  pKa > pH → 质子化 (保留 H)
  pKa < pH → 去质子化 (失去 H)

残基名映射 (ff14SB):
  ASP: pKa > pH → ASH (COOH, 0)   | pKa < pH → ASP (COO⁻, -1)
  GLU: pKa > pH → GLH (COOH, 0)   | pKa < pH → GLU (COO⁻, -1)
  HIS: pKa > pH → HIP (双质子, +1) | pKa < pH → HIE (ε位质子, 0)
  CYS: pKa > pH → CYS (SH, 0)    | pKa < pH → CYM (S⁻, -1)
  LYS: pKa > pH → LYS (NH₃⁺, +1) | pKa < pH → LYN (NH₂, 0)
  TYR: pKa > pH → TYR (OH, 0)    | pKa < pH → TYM (O⁻, -1)  [需自制参数]
  ARG: pKa > pH → ARG (+1)        | pKa < pH → ARM (去质子, 0) [需自制参数]
"""

# pKa > pH 时的残基名 (质子化), pKa < pH 时的残基名 (去质子化)
_PROT_MAP = {
    "ASP": ("ASH", "ASP"),
    "GLU": ("GLH", "GLU"),
    "HIS": ("HIP", "HIE"),
    "CYS": ("CYS", "CYM"),
    "LYS": ("LYS", "LYN"),
    "TYR": ("TYR", "TYM"),
    "ARG": ("ARG", "ARM"),
}

# 需要自制参数的去质子化残基名
_CUSTOM_PARAM = {"TYM", "ARM"}

# 这些 group 跳过不处理
_SKIP_GROUPS = {"N+", "C-"}


def parse_propka(pka_path):
    """解析 propka3 输出, 返回 [(group, resnum, chain, pka), ...]."""
    results = []
    in_summary = False
    with open(pka_path) as f:
        for line in f:
            if "SUMMARY OF THIS PREDICTION" in line:
                in_summary = True
                continue
            if in_summary and line.startswith("---"):
                break
            if not in_summary:
                continue
            # 跳过表头
            if "Group" in line and "pKa" in line:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            group = parts[0]
            if group in _SKIP_GROUPS:
                continue
            try:
                resnum = int(parts[1])
                chain = parts[2]
                pka = float(parts[3])
            except (ValueError, IndexError):
                continue
            results.append((group, resnum, chain, pka))
    return results


def assign_protonation(pka_results, ph, include_all=False):
    """根据 pKa 和 pH 分配残基名.

    include_all: True 时也替换 TYM/ARM (需自制参数), False 时仅输出 WARNING.
    返回: {(resnum, chain): new_resname}, warnings: [str]
    """
    assignments = {}
    warnings = []
    for group, resnum, chain, pka in pka_results:
        key = (resnum, chain)
        if group not in _PROT_MAP:
            continue
        protonated, deprotonated = _PROT_MAP[group]
        new_name = protonated if pka > ph else deprotonated
        # 需要自制参数的残基
        if new_name in _CUSTOM_PARAM and not include_all:
            warnings.append(
                f"  [WARNING] {group} {resnum} {chain} pKa={pka:.2f} < pH={ph:.1f}, "
                f"需去质子化为 {new_name} 但 ff14SB 无内置参数 "
                f"(使用 --all 强制替换)")
            continue
        # 只记录需要改名的
        if new_name != group:
            assignments[key] = new_name
    return assignments, warnings


def apply_protonation(input_pdb, output_pdb, assignments):
    """修改 PDB 中的残基名, 返回修改计数."""
    count = 0
    with open(input_pdb) as f:
        lines = f.readlines()

    out = []
    for line in lines:
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM") and len(line) >= 26:
            try:
                resnum = int(line[22:26])
                chain = line[21]
            except ValueError:
                out.append(line)
                continue
            key = (resnum, chain)
            if key in assignments:
                new_name = assignments[key]
                line = line[:17] + f"{new_name:>3s}" + line[20:]
                count += 1
        out.append(line)

    with open(output_pdb, "w") as f:
        f.writelines(out)
    return count


def protonate_pdb(input_pdb, pka_path, output_pdb, ph=7.0, include_all=False):
    """主入口: 解析 propka → 分配质子化状态 → 修改 PDB."""
    pka_results = parse_propka(pka_path)
    print(f"  propka 残基数: {len(pka_results)}")
    print(f"  pH = {ph:.1f}")

    assignments, warnings = assign_protonation(pka_results, ph, include_all)

    for w in warnings:
        print(w)

    if assignments:
        print(f"\n  质子化状态修改:")
        for (resnum, chain), new_name in sorted(assignments.items()):
            # 找原始 group name
            for group, rn, ch, pka in pka_results:
                if rn == resnum and ch == chain:
                    print(f"    {group:>3s} {resnum:4d} {chain} "
                          f"(pKa={pka:5.2f}) → {new_name}")
                    break
    else:
        print("  无需修改残基名")

    n = apply_protonation(input_pdb, output_pdb, assignments)
    n_res = len(assignments)
    print(f"\n  修改 {n_res} 个残基 ({n} 行), 保存: {output_pdb}")

    # 输出 summary.txt
    import os
    summary_path = os.path.join(os.path.dirname(os.path.abspath(output_pdb)),
                                "summary.txt")
    with open(summary_path, "w") as sf:
        sf.write(f"autoprep prot summary\n")
        sf.write(f"pH = {ph:.1f}\n")
        sf.write(f"input  = {input_pdb}\n")
        sf.write(f"output = {output_pdb}\n")
        sf.write(f"propka = {pka_path}\n\n")
        sf.write(f"{'Group':>5s} {'ResNum':>6s} {'Chain':>5s} "
                 f"{'pKa':>6s} {'->':>4s} {'NewName':>7s}\n")
        sf.write("-" * 42 + "\n")
        for (resnum, chain), new_name in sorted(assignments.items()):
            for group, rn, ch, pka in pka_results:
                if rn == resnum and ch == chain:
                    sf.write(f"{group:>5s} {resnum:6d} {chain:>5s} "
                             f"{pka:6.2f}   -> {new_name:>7s}\n")
                    break
        if warnings:
            sf.write(f"\nWarnings:\n")
            for w in warnings:
                sf.write(w.strip() + "\n")
        sf.write(f"\nTotal: {n_res} residues modified\n")
    print(f"  摘要: {summary_path}")
