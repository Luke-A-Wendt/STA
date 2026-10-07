"""Shared cyan/magenta PGFPlots styling for the paper's explanatory figures.

Only the standard library and the installed TeX tools are required.
White is the page background; all printed marks use cyan or magenta.
"""

from pathlib import Path
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parent

PREAMBLE = r"""\documentclass[11pt,border=3pt]{standalone}
\usepackage{amsmath,bm}
\usepackage{pgfplots}
\usepgfplotslibrary{groupplots}
\usetikzlibrary{arrows.meta,calc,decorations.markings,patterns}
\pgfplotsset{compat=1.18}
\pgfplotsset{
 every axis/.append style={
  width=7.65cm,height=6.8cm,
  axis line style={magenta},tick style={magenta},
  tick label style={font=\small,text=magenta},
  label style={font=\small,text=magenta},
  title style={font=\small\bfseries,text=magenta},
  grid=major,grid style={cyan!18},
  axis background/.style={fill=white},
  every axis plot/.append style={color=cyan,line width=1.1pt},
  legend style={font=\footnotesize,text=magenta,draw=none,
   fill=white,cells={anchor=west}},
  clip=true,
 }
}
\begin{document}
\color{magenta}
\begin{tikzpicture}[every node/.style={text=magenta}]
"""
END = "\\end{tikzpicture}\n\\end{document}\n"


def samples(lo, hi, n=241):
    return [lo+(hi-lo)*k/(n-1) for k in range(n)]


def xy(p):
    return f"{p[0]:.9f},{p[1]:.9f}"


def curve(points, options="cyan"):
    points = " ".join(f"({xy(p)})" for p in points)
    return f"\\addplot[{options}] coordinates {{{points}}};\n"


def line(start, end, options="cyan"):
    return (rf"\draw[{options}] (axis cs:{xy(start)}) -- "
            rf"(axis cs:{xy(end)});"+"\n")


def note(x, y, text, options=""):
    return (rf"\node[font=\footnotesize,fill=white,inner sep=1.5pt,{options}] "
            rf"at (axis cs:{x:.9f},{y:.9f}) {{{text}}};"+"\n")


def build_figure(name, body):
    """Compile a TikZ body to figures/<name>.pdf and report layout warnings."""
    destination = ROOT / "figures" / f"{name}.pdf"
    destination.parent.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"sta-{name}-") as directory:
        work = Path(directory)
        (work / f"{name}.tex").write_text(PREAMBLE+body+END)
        result = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", f"{name}.tex"],
            cwd=work, capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(result.stdout[-7000:])
        warnings = [line for line in result.stdout.splitlines()
                    if "Warning" in line or "Overfull" in line]
        if warnings:
            print(name+": "+"\n".join(warnings))
        shutil.copyfile(work/f"{name}.pdf", destination)
    print(f"Wrote {destination}")
    return destination
