"""
autoprep CLI 入口

用法:
    autoprep run config.json -o CYM_prep
    autoprep norm input.pdb -r CYM -o output.pdb
    autoprep prot input.pdb propka.pka -ph 7.0 -o output.pdb
    autoprep align protein.pdb --ref TYM.pdb --ref ARM.pdb -o output.pdb
"""

import argparse
import json
import os
import sys


def _cmd_run(args):
    """运行参数化流水线"""
    with open(args.config) as f:
        cfg = json.load(f)

    rn = cfg.get("residue_name", "RES")
    if args.outdir:
        workdir = os.path.abspath(args.outdir)
    elif cfg.get("workdir"):
        workdir = os.path.abspath(cfg["workdir"])
    else:
        workdir = os.path.abspath(f"{rn}_prep")

    # residue_file 相对路径基于 config 所在目录解析
    res_file = cfg.get("residue_file", "")
    if res_file and not os.path.isabs(res_file):
        cfg_dir = os.path.dirname(os.path.abspath(args.config))
        cfg["residue_file"] = os.path.join(cfg_dir, res_file)

    from autoprep.prep import run_pipeline
    run_pipeline(cfg, workdir)


def _cmd_prot(args):
    """根据 propka3 分配质子化状态"""
    output = args.output
    if output is None:
        base = args.input_pdb.rsplit(".", 1)[0]
        output = f"{base}_prot.pdb"

    from autoprep.protonate import protonate_pdb
    protonate_pdb(args.input_pdb, args.pka_file, output,
                  ph=args.ph, include_all=args.all)


def _cmd_align(args):
    """对齐残基原子名"""
    output = args.output
    if output is None:
        base = args.input_pdb.rsplit(".", 1)[0]
        output = f"{base}_aligned.pdb"

    from autoprep.align import align_pdb
    align_pdb(args.input_pdb, args.ref, output)


def _cmd_norm(args):
    """规范化 PDB 文件"""
    output = args.output
    if output is None:
        base = args.input_pdb.rsplit(".", 1)[0]
        output = f"{base}_norm.pdb"

    from autoprep.normalize import normalize_pdb
    normalize_pdb(args.input_pdb, output,
                  resname=args.resname.upper(),
                  chain=args.chain,
                  resnum=args.resnum,
                  keep_conect=not args.no_conect,
                  reorder=not args.no_sort)


def main():
    p = argparse.ArgumentParser(
        prog="autoprep",
        description="非标准氨基酸自动参数化 (Amber)")
    sub = p.add_subparsers(dest="command")

    # --- autoprep run ---
    p_run = sub.add_parser("run", help="运行参数化流水线")
    p_run.add_argument("config", help="配置文件 (JSON)")
    p_run.add_argument("-o", "--outdir", default=None,
                       help="输出目录 (默认: <residue_name>_prep)")

    # --- autoprep norm ---
    p_norm = sub.add_parser("norm", help="规范化 PDB 文件")
    p_norm.add_argument("input_pdb", help="输入 PDB 文件")
    p_norm.add_argument("-o", "--output", default=None,
                        help="输出文件 (默认: <input>_norm.pdb)")
    p_norm.add_argument("-r", "--resname", default="CYM",
                        help="残基名 (默认: CYM)")
    p_norm.add_argument("-c", "--chain", default="A",
                        help="链 ID (默认: A)")
    p_norm.add_argument("-n", "--resnum", type=int, default=1,
                        help="残基编号 (默认: 1)")
    p_norm.add_argument("--no-conect", action="store_true",
                        help="不保留 CONECT 记录")
    p_norm.add_argument("--no-sort", action="store_true",
                        help="不重排原子顺序 (保持原始顺序)")

    # --- autoprep prot ---
    p_prot = sub.add_parser("prot", help="根据 propka3 分配质子化状态")
    p_prot.add_argument("input_pdb", help="输入 PDB 文件")
    p_prot.add_argument("pka_file", help="propka3 输出文件 (.pka)")
    p_prot.add_argument("-ph", type=float, default=7.0,
                        help="pH 值 (默认: 7.0)")
    p_prot.add_argument("-o", "--output", default=None,
                        help="输出文件 (默认: <input>_prot.pdb)")
    p_prot.add_argument("--all", action="store_true",
                        help="包含 TYM/ARM 等需自制参数的替换")

    # --- autoprep align ---
    p_align = sub.add_parser("align", help="对齐残基原子名到参考 PDB 命名")
    p_align.add_argument("input_pdb", help="输入蛋白 PDB 文件")
    p_align.add_argument("--ref", action="append", required=True,
                         help="参考残基 PDB (autoprep run 输出), 可多次指定")
    p_align.add_argument("-o", "--output", default=None,
                         help="输出文件 (默认: <input>_aligned.pdb)")

    args = p.parse_args()
    if args.command == "run":
        _cmd_run(args)
    elif args.command == "norm":
        _cmd_norm(args)
    elif args.command == "prot":
        _cmd_prot(args)
    elif args.command == "align":
        _cmd_align(args)
    else:
        p.print_help()
        sys.exit(1)
