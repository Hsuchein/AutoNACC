# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/), and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.1.0] - 2025

### Added

- Core parameterization pipeline (`autoprep run`)
  - ACE/NME capping with trans peptide bond enforcement
  - RDKit MMFF94 pre-optimization
  - Gaussian Opt + ESP job generation (--Link1-- linked jobs)
  - Automatic Gaussian monitoring with checkpoint/resume
  - RESP charge fitting via antechamber
  - Prepin generation and charge verification/correction
  - Frcmod generation with ff14SB/GAFF2 cross-term patching
- PDB normalization utility (`autoprep norm`)
- Protonation state assignment from PropKa3 (`autoprep prot`)
- Atom name alignment to reference structures (`autoprep align`)
- Configuration via JSON files
- Comprehensive tutorial documentation
