<p align="center">
  <h1 align="center">AutoPrep</h1>
  <p align="center">
    <strong>Automated parameterization of non-standard amino acids for Amber molecular dynamics</strong>
  </p>
  <p align="center">
    <a href="#installation">Installation</a> &bull;
    <a href="#quick-start">Quick Start</a> &bull;
    <a href="#usage">Usage</a> &bull;
    <a href="#configuration">Configuration</a> &bull;
    <a href="docs/tutorial.md">Tutorial</a> &bull;
    <a href="#contributing">Contributing</a>
  </p>
</p>

---

## Overview

**AutoPrep** automates the tedious, multi-step process of generating Amber-compatible force field parameters (`.prepin` + `.frcmod`) for non-standard amino acid residues. It handles the full QM-to-MM pipeline:

```
Input PDB/XYZ  ──►  ACE-RES-NME Assembly  ──►  Gaussian Opt + ESP
                                                       │
                    .prepin + .frcmod  ◄──  RESP Fit  ◄─┘
```

### Key Features

- **End-to-end automation** — From a bare residue structure to production-ready Amber parameters in a single command
- **Smart capping** — Automatic ACE/NME cap placement with enforced trans peptide bond geometry (ω ≈ 180°)
- **MMFF pre-optimization** — RDKit-based pre-optimization before expensive QM calculations
- **Gaussian integration** — Generates linked Opt + ESP jobs (`--Link1--`) for efficiency, with automatic progress monitoring
- **RESP charge fitting** — Antechamber-driven RESP charges with automatic charge verification and correction
- **Cross-term patching** — Intelligent supplementation of ff14SB ↔ GAFF2 cross-term parameters
- **Checkpoint/resume** — Long-running pipelines can be safely interrupted and resumed
- **PDB normalization** — Standardize atom names, residue labels, and connectivity
- **Protonation assignment** — Assign pH-dependent protonation states from PropKa3 predictions
- **Atom name alignment** — Align parameterized residue naming to reference structures

## Prerequisites

| Software | Version | Purpose |
|----------|---------|---------|
| [Gaussian](https://gaussian.com/) | g09 / g16 | QM geometry optimization & ESP |
| [AmberTools](https://ambermd.org/AmberTools.php) | ≥ 20 | antechamber, prepgen, parmchk2 |
| [PropKa3](https://github.com/jensengroup/propka) | ≥ 3.0 | pKa prediction (for `prot` subcommand) |

## Installation

```bash
# Clone the repository
git clone https://github.com/Hsuchein/AutoNACC.git
cd AutoNACC

# Install (editable mode recommended for development)
pip install -e .

# Verify installation
autoprep --help
```

### Dependencies

- Python ≥ 3.8
- NumPy
- RDKit (optional, for MMFF pre-optimization)

## Quick Start

**1. Prepare your residue structure**

Export your non-standard residue as a PDB file (e.g., from GaussView). The structure should contain the bare residue fragment with free valence at backbone N and C atoms.

**2. Normalize the PDB** (if needed)

```bash
autoprep norm input.pdb -o residue.pdb -r PTR -c A -n 100
```

**3. Create a configuration file**

```bash
cp examples/config_template.json my_residue.json
```

Edit `my_residue.json` to specify your residue (see [Configuration](#configuration)).

**4. Run the parameterization pipeline**

```bash
autoprep run my_residue.json
```

AutoPrep will:
1. Assemble the ACE–RES–NME capped structure
2. Generate & run Gaussian Opt + ESP calculation
3. Perform RESP charge fitting via antechamber
4. Generate `.prepin` and `.frcmod` files
5. Output the final residue PDB with Amber-compatible atom names

## Usage

AutoPrep provides four subcommands:

### `autoprep run` — Parameterization Pipeline

Run the full parameterization workflow from a JSON config file.

```bash
autoprep run config.json
```

The pipeline creates a working directory named after the residue (e.g., `PTR/`) containing all intermediate and output files:

```
PTR/
├── assembled.pdb          # ACE-RES-NME structure
├── assembled_mmff.xyz     # MMFF pre-optimized structure
├── residue.com / .log     # Gaussian input / output
├── residue.ac             # Antechamber output
├── PTR.prepin             # ✅ Amber residue template
├── PTR.frcmod             # ✅ Amber force field modifications
├── PTR_residue.pdb        # ✅ Final residue PDB
└── .autoprep_meta.json    # Checkpoint metadata
```

### `autoprep norm` — PDB Normalization

Standardize PDB files: unify residue names, reorder atoms (backbone first), rename side-chain atoms, and update CONECT records.

```bash
autoprep norm input.pdb -o output.pdb -r CYM -c A -n 293
```

| Flag | Description | Default |
|------|-------------|---------|
| `-o` | Output file path | `<input>_norm.pdb` |
| `-r` | Residue name | `CYM` |
| `-c` | Chain ID | `A` |
| `-n` | Residue number | `293` |
| `--no-conect` | Drop CONECT records | keep |
| `--no-reorder` | Skip atom reordering | reorder |

### `autoprep prot` — Protonation Assignment

Assign protonation states based on PropKa3 pKa predictions.

```bash
autoprep prot protein.pdb propka.pka -o protonated.pdb --ph 7.0
```

| Flag | Description | Default |
|------|-------------|---------|
| `-o` | Output file path | `<input>_prot.pdb` |
| `--ph` | Target pH | `7.0` |
| `--include-all` | Include residues requiring custom params (TYM, ARM) | off |

Supported residue mappings:

| Residue | Protonated | Deprotonated |
|---------|------------|--------------|
| ASP | ASH | ASP |
| GLU | GLH | GLU |
| HIS | HIP | HIE |
| CYS | CYS | CYM |
| LYS | LYS | LYN |
| TYR | TYR | TYM* |
| ARG | ARG | ARM* |

*\*Requires custom parameters*

### `autoprep align` — Atom Name Alignment

Align atom names in a protein PDB to match reference residue naming.

```bash
autoprep align protein.pdb -ref ref1.pdb ref2.pdb -o aligned.pdb
```

## Configuration

The parameterization pipeline is configured via a JSON file. See [`examples/config_template.json`](examples/config_template.json) for a full template.

### Required Fields

```jsonc
{
    "residue_file": "residue.pdb",     // Input structure (PDB or XYZ)
    "residue_name": "PTR",             // 3-letter Amber residue name
    "net_charge":   -2,                // Total charge of ACE-RES-NME system
    "backbone_n_idx": 0,               // 0-indexed backbone nitrogen atom
    "backbone_c_idx": 2                // 0-indexed backbone carbonyl carbon
}
```

### Optional Fields

```jsonc
{
    // QM Methods
    "opt_method":   "wB97XD/6-31G*",  // Geometry optimization (default)
    "resp_method":  "HF/6-31G*",      // ESP for RESP fitting (default)
    "multiplicity": 1,                 // Spin multiplicity (default: 1)

    // Computational Resources
    "nproc": 16,                       // Gaussian CPU cores (default: 16)
    "mem":   "64GB",                   // Gaussian memory (default: 64GB)

    // Monitoring
    "check_interval": 60,             // Polling interval in seconds (default: 60)

    // Tool Paths (if not in $PATH)
    "gaussian_cmd": "g16",
    "antechamber":  "antechamber",
    "prepgen":      "prepgen",
    "parmchk2":     "parmchk2",

    // Advanced
    "residue_charge": -2,             // Override residue-only charge
    "mainchain_names": ["N","CA","C"] // Manual mainchain atom specification
}
```

## Project Structure

```
AutoNACC/
├── autoprep/                  # Python package
│   ├── __init__.py
│   ├── cli.py                 # CLI entry point & subcommand dispatch
│   ├── prep.py                # Core parameterization pipeline
│   ├── normalize.py           # PDB normalization utilities
│   ├── align.py               # Atom name alignment
│   ├── protonate.py           # Protonation state assignment
│   └── template/              # ACE/NME capping group templates
│       ├── ace.xyz
│       └── nme.xyz
├── docs/
│   └── tutorial.md            # Comprehensive parameterization tutorial
├── examples/
│   ├── config_template.json   # Example configuration file
│   └── cys_pho_norm.pdb       # Example normalized residue PDB
├── scripts/
│   ├── autoprep.py            # Standalone pipeline script
│   └── normalize_pdb.py       # Standalone normalization script
├── pyproject.toml             # Package configuration
├── LICENSE
├── CONTRIBUTING.md
└── CHANGELOG.md
```

## How It Works

### Pipeline Architecture

```
┌─────────────────────────────────────────────────────────┐
│                    autoprep run                          │
├─────────────────────────────────────────────────────────┤
│                                                         │
│  Step 1   Read residue structure (PDB / XYZ)            │
│     │                                                   │
│  Step 2   Assemble ACE ─ RES ─ NME                      │
│     │     ├─ Place caps at backbone N / C                │
│     │     ├─ Set C─N bond = 1.335 Å                     │
│     │     └─ Enforce ω dihedral ≈ 180° (trans)          │
│     │                                                   │
│  Step 3   MMFF94 pre-optimization (RDKit)               │
│     │                                                   │
│  Step 4   Generate Gaussian input                       │
│     │     ├─ Job 1: Geometry optimization                │
│     │     └─ Job 2: ESP for RESP (via --Link1--)        │
│     │                                                   │
│  Step 5   Run Gaussian & monitor progress               │
│     │                                                   │
│  Step 6   RESP charge fitting (antechamber)             │
│     │                                                   │
│  Step 7   Generate .prepin (prepgen)                    │
│     │     └─ Verify & fix charges                       │
│     │                                                   │
│  Step 8   Generate .frcmod (parmchk2)                   │
│     │     └─ Patch ff14SB ↔ GAFF2 cross-terms           │
│     │                                                   │
│  Step 9   Output final residue PDB                      │
│                                                         │
└─────────────────────────────────────────────────────────┘
```

### Checkpoint & Resume

Each step checks for existing output files. If a step has already completed, it is skipped. This allows safe interruption and resumption of long-running Gaussian calculations:

```bash
# Start the pipeline (Gaussian may take hours)
autoprep run config.json

# Safe to Ctrl+C during Gaussian monitoring

# Resume from where it left off
autoprep run config.json
```

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Citation

If you use AutoPrep in your research, please cite:

```bibtex
@software{autoprep,
  title  = {AutoPrep: Automated Parameterization of Non-Standard Amino Acids for Amber MD},
  author = {Hsuchein},
  url    = {https://github.com/Hsuchein/AutoNACC},
  year   = {2025}
}
```

## Acknowledgments

- [AmberTools](https://ambermd.org/) for antechamber, prepgen, and parmchk2
- [Gaussian](https://gaussian.com/) for QM calculations
- [RDKit](https://www.rdkit.org/) for MMFF pre-optimization
- [PropKa3](https://github.com/jensengroup/propka) for pKa predictions
