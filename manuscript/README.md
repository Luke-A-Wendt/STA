# Manuscript

This directory contains the single maintained LaTeX manuscript for the STA project.

## Files

- `main.tex`: primary manuscript source.
- `generated/python_verification.tex`: auto-generated symbolic verification table included by `main.tex`.
- `generated/readme_equation_catalog.tex`, `generated/readme_full_conversion.tex`, and `generated/integrated_sections/*.tex`: build-time preservation artifacts generated from `README.txt` and kept for traceability during manuscript assembly.
- `main.pdf`: compiled PDF output after a successful build.

## Build

From the repository root, run:

```bash
./build_latex.sh
```

The build script regenerates the README conversion files and compiles `manuscript/main.tex` with `lualatex`.
