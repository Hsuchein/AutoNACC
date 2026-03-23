#!/usr/bin/env python3
"""
autoprep.py - 非标准氨基酸自动参数化脚本 (Amber)

用法:  python autoprep.py config.json

工作流:
    1. 读取残基片段(肽键断开,有自由价) + ACE/NME帽基模板
    2. 组装 ACE-RES-NME, 确保肽键平面性 (O=C-N-H 共面)
    3. 生成 Gaussian 输入 (--Link1-- 串联两步):
       Job1: Opt opt_method (如 wB97X-D3/6-31G*)
       Job2: SP  resp_method Geom=AllCheck (如 HF/6-31G* Pop=MK)
    4. 运行 Gaussian, 每分钟检查是否完成
    5. antechamber RESP 电荷拟合 (单构象)
    6. prepgen 生成 .prepin
    7. parmchk2 生成 .frcmod
"""

import numpy as np
import subprocess as sp
import time
import os
import sys
import json
from pathlib import Path

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
                    # 处理像 "1H" "CA" 等名字
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


# =============================================================================
# 几何工具
# =============================================================================

def normalize(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-10 else np.zeros(3)


def find_bonded(coords, elems, idx):
    """基于共价半径找与 idx 成键的原子"""
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
    """找与 idx 成键的第一个指定元素原子, 返回索引或 None"""
    for j in find_bonded(coords, elems, idx):
        if elems[j] == target_elem:
            return j
    return None


def outgoing_direction(coords, elems, idx):
    """sp2/sp3 原子的自由价方向: 现有键向量均值的反方向"""
    bonded = find_bonded(coords, elems, idx)
    if not bonded:
        return np.array([1.0, 0.0, 0.0])
    vecs = [normalize(coords[j] - coords[idx]) for j in bonded]
    return normalize(-np.mean(vecs, axis=0))


def rotation_matrix_align(a, b):
    """Rodrigues公式: 将单位向量 a 旋转到 b"""
    a, b = normalize(a), normalize(b)
    v = np.cross(a, b)
    c = np.dot(a, b)
    if np.linalg.norm(v) < 1e-10:
        if c > 0:
            return np.eye(3)
        perp = np.array([1, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1, 0])
        perp = normalize(perp - np.dot(perp, a) * a)
        return 2 * np.outer(perp, perp) - np.eye(3)
    vx = np.array([[0, -v[2], v[1]],
                    [v[2], 0, -v[0]],
                    [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx / (1 + c)


def dihedral_angle(p1, p2, p3, p4):
    """计算二面角 p1-p2-p3-p4, 返回弧度 [-π, π]"""
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
    y = np.dot(np.cross(n1, normalize(b2)), n2)
    return np.arctan2(y, x)


def rotate_points_around_axis(points, origin, axis, angle):
    """将点集绕轴旋转 angle 弧度.
    origin: 轴上一点, axis: 轴方向, angle: 旋转角(弧度)"""
    k = normalize(axis)
    ca, sa = np.cos(angle), np.sin(angle)
    centered = points - origin
    rotated = np.zeros_like(centered)
    for i in range(len(centered)):
        v = centered[i]
        rotated[i] = v * ca + np.cross(k, v) * sa + k * np.dot(k, v) * (1 - ca)
    return rotated + origin


def place_cap(cap_coords, cap_elems, connect_idx, target_pos, bond_dir):
    """放置帽基: connect_idx 原子移到 target_pos, 自由价对准 bond_dir"""
    cap_out = outgoing_direction(cap_coords, cap_elems, connect_idx)
    R = rotation_matrix_align(cap_out, bond_dir)
    centered = cap_coords - cap_coords[connect_idx]
    rotated = (R @ centered.T).T
    return rotated + target_pos


# =============================================================================
# 组装 ACE-RES-NME (带肽键平面性)
# =============================================================================

def assemble(res_e, res_c, res_n, ace_e, ace_c, nme_e, nme_c, cfg):
    """组装加帽结构, 确保肽键 O=C-N-H 共面 (trans, ω≈180°).

    ACE 的自由价(羰基C) → 残基主链 N
    残基主链 C    → NME 的自由价(酰胺N)

    返回 dict 包含组装后的所有信息.
    """
    bb_n = cfg["backbone_n_idx"]   # 残基中主链 N 的索引
    bb_c = cfg["backbone_c_idx"]   # 残基中主链 C(=O) 的索引
    ace_ci = cfg.get("ace_connect_idx", 0)  # ACE 羰基C索引
    ace_oi = cfg.get("ace_o_idx", 1)        # ACE 羰基O索引
    nme_ci = cfg.get("nme_connect_idx", 4)  # NME 酰胺N索引
    nme_hi = cfg.get("nme_h_idx", 5)        # NME N上H索引

    n_pos = res_c[bb_n]
    c_pos = res_c[bb_c]

    # 残基上与肽键平面性相关的原子
    # N 上的 H (酰胺H): 用于 ACE 侧的平面性
    res_h_on_n = find_bonded_elem(res_c, res_e, bb_n, 'H')
    # C 上的 O (羰基O): 用于 NME 侧的平面性
    res_o_on_c = find_bonded_elem(res_c, res_e, bb_c, 'O')

    # --- 自由价方向 ---
    v_n = outgoing_direction(res_c, res_e, bb_n)  # N 的自由价方向 (→ ACE_C)
    v_c = outgoing_direction(res_c, res_e, bb_c)  # C 的自由价方向 (→ NME_N)

    # =====================================================
    # 放置 ACE: ACE_C → 残基_N
    # =====================================================
    ace_target = n_pos + PEPTIDE_CN * v_n       # ACE_C 的目标位置
    ace_bond = normalize(n_pos - ace_target)     # ACE_C→N 方向
    ace_placed = place_cap(ace_c, ace_e, ace_ci, ace_target, ace_bond)

    # 肽键平面性: 旋转 ACE 绕 C-N 轴, 使 dihedral(ACE_O, ACE_C, N, H_on_N) ≈ 0°
    # (trans 肽键中 O 和 H 在同侧, 即 cis 关系, 二面角≈0°)
    if res_h_on_n is not None:
        cur = dihedral_angle(ace_placed[ace_oi], ace_placed[ace_ci],
                             n_pos, res_c[res_h_on_n])
        rot_angle = 0.0 - cur  # 目标: 0°
        axis_dir = n_pos - ace_placed[ace_ci]
        ace_placed = rotate_points_around_axis(
            ace_placed, ace_placed[ace_ci], axis_dir, rot_angle)
        final_dihed = np.degrees(dihedral_angle(
            ace_placed[ace_oi], ace_placed[ace_ci], n_pos, res_c[res_h_on_n]))
        print(f"    ACE肽键平面: dihedral(O,C,N,H) = {final_dihed:.1f}°  (目标 0°)")

    # =====================================================
    # 放置 NME: 残基_C → NME_N
    # =====================================================
    nme_target = c_pos + PEPTIDE_CN * v_c       # NME_N 的目标位置
    nme_bond = normalize(c_pos - nme_target)     # NME_N→C 方向
    nme_placed = place_cap(nme_c, nme_e, nme_ci, nme_target, nme_bond)

    # 肽键平面性: 旋转 NME 绕 C-N 轴, 使 dihedral(O_on_C, C, NME_N, NME_H) ≈ 0°
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

    # =====================================================
    # 合并: ACE + RES + NME (残基原子不做任何删除)
    # =====================================================
    n_ace = len(ace_e)
    n_res = len(res_e)
    n_nme = len(nme_e)

    all_e = list(ace_e) + list(res_e) + list(nme_e)
    all_c = np.vstack([ace_placed, res_c, nme_placed])

    # 统一生成原子名 (Gaussian 只用元素, 但后续 .ac 需要)
    ctr = {}
    all_n = []
    for e in all_e:
        ctr[e] = ctr.get(e, 0) + 1
        all_n.append(f"{e}{ctr[e]}")

    return {
        "elements": all_e,
        "coords": all_c,
        "atom_names": all_n,
        "n_ace": n_ace,
        "n_res": n_res,
        "n_nme": n_nme,
        # 主链 N/C 在残基部分中的索引 (不变, 因为没有删除原子)
        "bb_n_in_res": bb_n,
        "bb_c_in_res": bb_c,
    }


# =============================================================================
# Gaussian
# =============================================================================

def write_gaussian_input(path, elems, coords, cfg):
    """写 Gaussian 输入: Opt + --Link1-- SP ESP, 共享 chk"""
    chk = os.path.basename(path).replace('.com', '.chk')
    rn = cfg["residue_name"]
    opt_method = cfg.get("opt_method", cfg.get("method", "HF/6-31G*"))
    resp_method = cfg.get("resp_method", "HF/6-31G*")
    chg = cfg["net_charge"]
    mult = cfg["multiplicity"]

    with open(path, 'w') as f:
        # ---- Job 1: 结构优化 ----
        f.write(f"%chk={chk}\n")
        f.write(f"%nproc={cfg['nproc']}\n")
        f.write(f"%mem={cfg['mem']}\n")
        f.write(f"# Opt {opt_method} SCF=Tight\n")
        f.write(f"\nACE-{rn}-NME optimization\n\n")
        f.write(f"{chg} {mult}\n")
        for e, (x, y, z) in zip(elems, coords):
            f.write(f" {e:<2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")
        f.write("\n")

        # ---- Link1: 单点 ESP ----
        f.write("--Link1--\n")
        f.write(f"%chk={chk}\n")
        f.write(f"%nproc={cfg['nproc']}\n")
        f.write(f"%mem={cfg['mem']}\n")
        f.write(f"# {resp_method} Geom=AllCheck Guess=Read "
                f"Pop=MK IOp(6/33=2) SCF=Tight\n")
        f.write("\n")


def run_gaussian(com_file, cfg):
    log_file = os.path.splitext(com_file)[0] + ".log"
    cmd = f"{cfg['gaussian_cmd']} {os.path.basename(com_file)}"
    cwd = os.path.dirname(os.path.abspath(com_file))
    print(f"  [RUN] cd {cwd} && {cmd}")
    sp.Popen(cmd, shell=True, cwd=cwd)
    return log_file


def monitor_gaussian(log_file, interval=60):
    """每 interval 秒检查 Gaussian 状态. 返回 (success, actual_file)."""
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
            print(f"    [{mins}min] Gaussian 正常结束!")
            return True, target
        if "Error termination" in content:
            print(f"    [{mins}min] ERROR: Gaussian 异常终止!")
            return False, target

        nsteps = content.count("Step number")
        print(f"    [{mins}min] 运行中, 优化步数 ~{nsteps}")


# =============================================================================
# AmberTools
# =============================================================================

def run_cmd(cmd, cwd=None):
    print(f"  [CMD] {cmd}")
    r = sp.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
    if r.returncode != 0:
        print(f"    STDERR: {r.stderr[:800]}")
    return r.returncode, r.stdout, r.stderr


def read_ac_names(ac_path):
    """从 .ac 文件读取原子名列表 (保持顺序)"""
    names = []
    with open(ac_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                names.append(line.split()[1])
    return names


def run_antechamber(gout_file, cfg, workdir):
    rn = cfg["residue_name"]
    ac = f"{rn}.ac"
    cmd = (f"{cfg['antechamber']} "
           f"-i {os.path.basename(gout_file)} -fi gout "
           f"-o {ac} -fo ac "
           f"-c resp -s 2 "
           f"-rn {rn} -at amber "
           f"-nc {cfg['net_charge']} -pf y")
    rc, _, err = run_cmd(cmd, cwd=workdir)
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
    cmd = (f"{cfg['prepgen']} "
           f"-i {os.path.basename(ac_file)} "
           f"-o {out} "
           f"-m {os.path.basename(mc_file)} "
           f"-rn {rn}")
    rc, _, err = run_cmd(cmd, cwd=workdir)
    if rc != 0:
        raise RuntimeError(f"prepgen 失败:\n{err}")
    return os.path.join(workdir, out)


def run_parmchk2(prepin_file, cfg, workdir):
    rn = cfg["residue_name"]
    out = f"{rn}.frcmod"
    cmd = (f"{cfg['parmchk2']} "
           f"-i {os.path.basename(prepin_file)} -f prepi "
           f"-o {out}")
    rc, _, err = run_cmd(cmd, cwd=workdir)
    if rc != 0:
        raise RuntimeError(f"parmchk2 失败:\n{err}")
    return os.path.join(workdir, out)


# =============================================================================
# 自动检测
# =============================================================================

def auto_find_ca(coords, elems, bb_n, bb_c):
    """自动查找 Cα: 同时与主链 N 和 C 成键的碳"""
    bn = set(find_bonded(coords, elems, bb_n))
    bc = set(find_bonded(coords, elems, bb_c))
    candidates = [i for i in bn & bc if elems[i] == 'C']
    if candidates:
        return candidates[0]
    # 放宽: N-X-C 路径中的 X
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
# 主流程
# =============================================================================

def main():
    if len(sys.argv) < 2:
        print("用法: python autoprep.py config.json")
        sys.exit(1)

    with open(sys.argv[1]) as f:
        cfg = json.load(f)

    rn = cfg["residue_name"]
    workdir = os.path.abspath(cfg.get("workdir", rn + "_prep"))
    os.makedirs(workdir, exist_ok=True)

    print("=" * 60)
    print("  autoprep.py - 非标准氨基酸自动参数化")
    opt_m = cfg.get("opt_method", cfg.get("method", "HF/6-31G*"))
    resp_m = cfg.get("resp_method", "HF/6-31G*")
    print(f"  残基: {rn}   电荷: {cfg['net_charge']}")
    print(f"  优化: {opt_m}   RESP: {resp_m}")
    print(f"  工作目录: {workdir}")
    print("=" * 60)

    # ==================== Step 1: 读取输入 ====================
    print("\n[Step 1] 读取输入 ...")
    res_e, res_c, res_n = read_structure(cfg["residue_file"])
    ace_e, ace_c = read_xyz(cfg["ace_file"])
    nme_e, nme_c = read_xyz(cfg["nme_file"])
    print(f"  残基: {len(res_e)}  ACE: {len(ace_e)}  NME: {len(nme_e)}  原子")
    print_atoms(res_e, res_c, res_n, "残基原子 (请确认 backbone_n/c 索引)")

    # ==================== Step 2: 组装 ====================
    print(f"\n[Step 2] 组装 ACE-{rn}-NME ...")
    asm = assemble(res_e, res_c, res_n, ace_e, ace_c, nme_e, nme_c, cfg)
    total = len(asm["elements"])
    na, nr, nn = asm["n_ace"], asm["n_res"], asm["n_nme"]
    print(f"  总原子: {total}  (ACE:{na} + RES:{nr} + NME:{nn})")

    # 保存 XYZ 供检查
    xyz_out = os.path.join(workdir, f"{rn}_capped.xyz")
    write_xyz(xyz_out, asm["elements"], asm["coords"],
              f"ACE-{rn}-NME capped structure")
    print(f"  保存: {xyz_out}  ← 请用分子可视化软件检查!")

    opt_method = cfg.get("opt_method", cfg.get("method", "HF/6-31G*"))
    resp_method = cfg.get("resp_method", "HF/6-31G*")

    # ==================== Step 3: Gaussian 输入 ====================
    print(f"\n[Step 3] Gaussian 输入 (Opt + --Link1-- SP ESP) ...")
    com = os.path.join(workdir, f"{rn}.com")
    write_gaussian_input(com, asm["elements"], asm["coords"], cfg)
    print(f"  文件: {com}")
    print(f"  Job1: Opt {opt_method}")
    print(f"  Job2: SP  {resp_method} Pop=MK IOp(6/33=2)")
    print(f"  nproc={cfg['nproc']}  mem={cfg['mem']}")

    # ==================== Step 4: 运行 Gaussian ====================
    print(f"\n[Step 4] 运行 Gaussian ...")
    log = run_gaussian(com, cfg)
    ok, actual_log = monitor_gaussian(log, cfg.get("check_interval", 60))
    if not ok:
        print("Gaussian 失败, 请检查后重新运行.")
        sys.exit(1)

    # ==================== Step 5: antechamber RESP ====================
    print(f"\n[Step 5] antechamber RESP 拟合 ...")
    # 确保 log 在 workdir 中
    if os.path.dirname(os.path.abspath(actual_log)) != workdir:
        import shutil
        shutil.copy2(actual_log, workdir)
        actual_log = os.path.join(workdir, os.path.basename(actual_log))

    ac_file = run_antechamber(actual_log, cfg, workdir)
    print(f"  .ac: {ac_file}")

    # 从 .ac 文件读原子名 → 映射 ACE/RES/NME
    ac_names = read_ac_names(ac_file)
    ace_ac = ac_names[:na]
    res_ac = ac_names[na:na + nr]
    nme_ac = ac_names[na + nr:]

    bb_n_in_res = asm["bb_n_in_res"]
    bb_c_in_res = asm["bb_c_in_res"]
    head = res_ac[bb_n_in_res]
    tail = res_ac[bb_c_in_res]
    print(f"  HEAD={head} (主链N)  TAIL={tail} (主链C)")

    # 自动检测 CA
    ca_idx = auto_find_ca(res_c, res_e, bb_n_in_res, bb_c_in_res)
    mc = [head]
    if ca_idx is not None:
        mc.append(res_ac[ca_idx])
        print(f"  自动检测 CA: {res_ac[ca_idx]} (残基内索引 {ca_idx})")
    mc.append(tail)

    if cfg.get("mainchain_names"):
        mc = cfg["mainchain_names"]
        print(f"  使用用户指定 MAIN_CHAIN: {mc}")

    # ==================== Step 6: prepgen ====================
    print(f"\n[Step 6] prepgen → .prepin ...")
    mc_file = os.path.join(workdir, f"{rn}.mc")
    charge = cfg.get("residue_charge", cfg["net_charge"])
    write_mainchain_mc(mc_file, ace_ac, nme_ac, head, tail, mc, charge)
    print(f"  mainchain: {mc_file}")

    prepin = run_prepgen(ac_file, mc_file, cfg, workdir)
    print(f"  生成: {prepin}")

    # ==================== Step 7: parmchk2 ====================
    print(f"\n[Step 7] parmchk2 → .frcmod ...")
    frcmod = run_parmchk2(prepin, cfg, workdir)
    print(f"  生成: {frcmod}")

    # ==================== 完成 ====================
    print("\n" + "=" * 60)
    print("  完成! 输出文件:")
    print(f"    {prepin}")
    print(f"    {frcmod}")
    print(f"\n  tleap 加载示例:")
    print(f"    source leaprc.protein.ff14SB")
    print(f"    loadamberprep {rn}.prepin")
    print(f"    loadamberparams {rn}.frcmod")
    print(f"    x = loadpdb your_protein.pdb")
    print(f"    check x")
    print("=" * 60)


if __name__ == "__main__":
    main()
