#!/usr/bin/env python3
"""
normalize_pdb.py
Normalize a GaussView-edited PDB file for use with autoprep.py.
- Convert all HETATM to ATOM
- Assign a single residue name (default CYM)
- Remove TER / REMARK / ANISOU lines
- Assign sequential atom names per element (e.g. C1, C2, H1, H2, ...)
  while preserving backbone atom names (N, CA, C, O, CB)
- Renumber atoms sequentially
- Optionally keep CONECT records with updated numbering
"""

import argparse
import sys
from collections import defaultdict


BACKBONE_NAMES = {"N", "CA", "C", "O", "CB"}


def parse_args():
    p = argparse.ArgumentParser(description="Normalize a GaussView PDB for Amber prep")
    p.add_argument("input_pdb", help="Input PDB file (e.g. cys_pho.pdb)")
    p.add_argument("-o", "--output", default=None,
                   help="Output PDB file (default: <input>_norm.pdb)")
    p.add_argument("-r", "--resname", default="CYM",
                   help="Residue name to assign (default: CYM)")
    p.add_argument("-c", "--chain", default="A",
                   help="Chain ID (default: A)")
    p.add_argument("-n", "--resnum", type=int, default=293,
                   help="Residue number (default: 293)")
    p.add_argument("--keep-conect", action="store_true", default=True,
                   help="Keep CONECT records with updated numbering (default: True)")
    p.add_argument("--no-conect", action="store_true",
                   help="Remove CONECT records")
    return p.parse_args()


def normalize_pdb(input_pdb, output_pdb, resname="CYM", chain="A", resnum=293,
                  keep_conect=True):
    with open(input_pdb) as f:
        lines = f.readlines()

    atom_lines = []   # (line, is_original_atom)
    conect_lines = []
    old_to_new = {}  # old serial -> new serial

    for line in lines:
        rec = line[:6].strip()
        if rec == "ATOM":
            atom_lines.append((line, True))
        elif rec == "HETATM":
            atom_lines.append((line, False))
        elif rec == "CONECT":
            conect_lines.append(line)

    # Track element counts for naming non-backbone atoms
    elem_count = defaultdict(int)

    out_atoms = []
    for i, (line, is_orig_atom) in enumerate(atom_lines):
        new_serial = i + 1
        old_serial = int(line[6:11])
        old_to_new[old_serial] = new_serial

        old_name = line[12:16].strip()
        element = line[76:78].strip()
        if not element:
            # Fallback: guess from atom name
            element = old_name[0]

        # Determine atom name: preserve backbone names only for original ATOM records
        if is_orig_atom and old_name in BACKBONE_NAMES:
            atom_name = old_name
        else:
            elem_count[element] += 1
            atom_name = f"{element}{elem_count[element]}"

        # Format atom name in PDB columns 13-16
        # PDB convention: 1-char elements right-justified in col 13-14,
        # 2+ char names start at col 13
        if len(atom_name) <= 3:
            name_field = f" {atom_name:<3s}"
        else:
            name_field = f"{atom_name:<4s}"

        # Build ATOM line (strict PDB column format)
        # cols: 1-6 record, 7-11 serial, 12 space, 13-16 name,
        #       17 altloc, 18-20 resname, 21 space, 22 chain,
        #       23-26 resseq, 27 icode, 28-30 spaces,
        #       31-38 x, 39-46 y, 47-54 z, 55-60 occ, 61-66 bfac,
        #       77-78 element
        out = (f"ATOM  {new_serial:5d} {name_field}"
               f" {resname:>3s} {chain:1s}{resnum:4d}    "
               f"{line[30:54]}"           # x, y, z
               f"{1.0:6.2f}{0.0:6.2f}"   # occupancy, bfactor
               f"          {element:>2s}  \n")
        out_atoms.append(out)

    # Build CONECT records with updated numbering
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

    # Write output
    with open(output_pdb, "w") as f:
        for line in out_atoms:
            f.write(line)
        for line in out_conect:
            f.write(line)
        f.write("END\n")

    print(f"Wrote {len(out_atoms)} atoms to {output_pdb}")
    print(f"Residue name: {resname}, Chain: {chain}, ResNum: {resnum}")

    # Print atom table for verification
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


if __name__ == "__main__":
    args = parse_args()
    output = args.output
    if output is None:
        base = args.input_pdb.rsplit(".", 1)[0]
        output = f"{base}_norm.pdb"

    keep = args.keep_conect and not args.no_conect
    normalize_pdb(args.input_pdb, output,
                  resname=args.resname.upper(),
                  chain=args.chain,
                  resnum=args.resnum,
                  keep_conect=keep)
