# Contributing to AutoPrep

Thank you for your interest in contributing to AutoPrep! This document provides guidelines to make the contribution process smooth.

## Getting Started

1. Fork the repository
2. Clone your fork:
   ```bash
   git clone https://github.com/<your-username>/AutoNACC.git
   cd AutoNACC
   ```
3. Install in development mode:
   ```bash
   pip install -e .
   ```
4. Create a feature branch:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## Development Guidelines

### Code Style

- Follow [PEP 8](https://peps.python.org/pep-0008/) conventions
- Use meaningful variable names, especially for atomic/molecular quantities
- Keep functions focused — one function, one responsibility

### Commit Messages

- Use clear, descriptive commit messages
- Start with a verb in imperative mood (e.g., "Add support for ...", "Fix charge rounding in ...")

### Testing

- Test with real non-standard residues when possible
- Verify that generated `.prepin` and `.frcmod` files load correctly in tLEaP
- Check charge conservation after RESP fitting

## Submitting Changes

1. Ensure your code runs without errors
2. Commit your changes with a clear message
3. Push to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```
4. Open a Pull Request against the `main` branch
5. Describe what your changes do and why

## Reporting Issues

When reporting bugs, please include:

- Python version and OS
- Gaussian / AmberTools version
- The configuration JSON used (sensitive paths redacted)
- Full error output or log file
- Input PDB file if possible

## Feature Requests

Feature requests are welcome! Please open an issue describing:

- The problem you're trying to solve
- Your proposed solution (if any)
- Any relevant references or literature

## Questions?

Open an issue with the **question** label — we're happy to help.
