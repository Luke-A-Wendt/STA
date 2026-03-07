# Manuscript

This directory contains the single maintained LaTeX manuscript for the STA project.

## Files

- `main.tex`: primary manuscript source.
- `generated/readme_equation_catalog.tex`: auto-generated catalog of equation-bearing README lines.
- `generated/readme_full_conversion.tex`: auto-generated full line-by-line README conversion included in the manuscript body.
- `main.pdf`: compiled PDF output after a successful build.

## Build

From the repository root, run:

```bash
./build_latex.sh
```

The build script regenerates the README conversion files and compiles `manuscript/main.tex` with `lualatex`.
