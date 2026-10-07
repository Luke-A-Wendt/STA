#!/usr/bin/env python3
"""Build the Section 19 vector figures using Python's stdlib and PGFPlots.

Run: python plot_metric_examples.py
Requires pdflatex with standalone, TikZ, and PGFPlots. No Python packages.
The output PDFs go in figures/; temporary TeX build files stay outside the repo.
All black-hole coordinates are divided by Gm. The expansion plots use t0.
"""

from math import cos, hypot, isclose, pi, sin, sqrt
from pathlib import Path
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "figures"
OMEGA = 0.25  # Gm * horizon angular velocity; spin points along +r_3.
ELL_SQUARED = 4 / (1 + 4 * OMEGA**2)
SPIN = ELL_SQUARED * OMEGA  # Specific angular momentum divided by Gm.
HORIZON = ELL_SQUARED / 2  # Ellipsoidal radius epsilon_+ / (Gm).

PREAMBLE = r"""\documentclass[11pt,border=3pt]{standalone}
\usepackage{amsmath}
\usepackage{pgfplots}
\usepgfplotslibrary{groupplots}
\usetikzlibrary{arrows.meta}
\pgfplotsset{compat=1.18}
\colorlet{metriccyan}{cyan}
\colorlet{framemagenta}{magenta}
\pgfplotsset{
 every axis/.append style={
  width=7.65cm,height=6.8cm,
  axis line style={magenta},tick style={magenta},
  tick label style={font=\small,text=magenta},label style={font=\small,text=magenta},
  title style={font=\small\bfseries,text=magenta},
  grid=major,grid style={cyan!18},
  axis background/.style={fill=white},
  every axis plot/.append style={line width=0.9pt},
  legend style={font=\footnotesize,text=magenta,draw=none,fill=white,cells={anchor=west}},
  clip=true,
 }
}
\begin{document}
\color{magenta}
\begin{tikzpicture}[every node/.style={text=magenta}]
\begin{groupplot}[group style={group size=2 by 1,horizontal sep=1.25cm}]
"""
END = "\\end{groupplot}\n\\end{tikzpicture}\n\\end{document}\n"


def samples(lo, hi, n=241):
    return [lo + (hi - lo) * k / (n - 1) for k in range(n)]


def xy(p):
    return f"{p[0]:.8f},{p[1]:.8f}"


def curve(points, options):
    coordinates = " ".join(f"({xy(p)})" for p in points)
    return f"\\addplot[{options}] coordinates {{{coordinates}}};\n"


def circle(radius):
    return [(radius * cos(t), radius * sin(t)) for t in samples(0, 2 * pi)]


def ellipse(epsilon):
    return [(sqrt(epsilon**2 + SPIN**2) * sin(t), epsilon * cos(t))
            for t in samples(0, 2 * pi)]


def note(x, y, text, options=""):
    return (rf"\node[font=\footnotesize,fill=white,inner sep=1.5pt,{options}] "
            rf"at (axis cs:{x:.8f},{y:.8f}) {{{text}}};" + "\n")


def line(start, end, options):
    return (rf"\draw[{options}] (axis cs:{xy(start)}) -- "
            rf"(axis cs:{xy(end)});" + "\n")


def cone(x, t, left_speed, right_speed, dt):
    left = (x + left_speed * dt, t + dt)
    right = (x + right_speed * dt, t + dt)
    tex = (rf"\fill[metriccyan!12] (axis cs:{x},{t}) -- "
           rf"(axis cs:{xy(left)}) -- (axis cs:{xy(right)}) -- cycle;" + "\n")
    for endpoint in (left, right):
        tex += line((x, t), endpoint, "metriccyan,line width=0.9pt,-{Stealth[length=3pt]}")
    return tex


def expansion():
    tex = r"""\nextgroupplot[
 title={Expansion: local light cones},
 xlabel={$r_1/t_0$},ylabel={$t/t_0$},
 xmin=-2,xmax=2,ymin=0,ymax=2.05,
 xtick={-2,-1,0,1,2},ytick={0,0.5,1,1.5,2},
 legend style={at={(0.5,-0.21)},anchor=north},
]
"""
    for x in (-1.5, -0.75, 0.75, 1.5):
        tex += line((x, 0), (x, 2.05), "metriccyan!40,densely dotted")
    tex += curve([(0, 0), (0, 2.05)], "framemagenta,line width=1.2pt")
    tex += "\\addlegendentry{Comoving observer}\n"
    for t in (0.15, 0.75, 1.35):
        for x in (-1.25, 0, 1.25):
            tex += cone(x, t, -1 / (1 + t), 1 / (1 + t), 0.38)
    tex += r"\addlegendimage{metriccyan,-{Stealth[length=3pt]}}" + "\n"
    tex += "\\addlegendentry{Local light directions}\n"
    tex += note(0, 1.92, r"$\lambda(t)=1+t/t_0$")
    tex += r"""\nextgroupplot[
 title={Separation and momentum},
 xlabel={$t/t_0$},ylabel={Ratio to its value at $t=0$},
 xmin=0,xmax=2,ymin=0,ymax=3.2,
 xtick={0,0.5,1,1.5,2},ytick={0,1,2,3},
 legend style={at={(0.5,-0.21)},anchor=north},
]
"""
    tex += curve([(t, 1+t) for t in samples(0, 2)], "metriccyan,line width=1.2pt")
    tex += r"\addlegendentry{Separation $\ell/\ell_0=\lambda$}" + "\n"
    tex += curve([(t, 1/(1+t)) for t in samples(0, 2)], "framemagenta,dashed,line width=1.2pt")
    tex += r"\addlegendentry{Momentum $P/P_0=1/\lambda$}" + "\n"
    return tex


def schwarzschild():
    tex = r"""\nextgroupplot[
 title={Schwarzschild: spatial contours},
 xlabel={$r_1/(Gm)$},ylabel={$r_2/(Gm)$},
 xmin=-5.65,xmax=5.65,ymin=-5.65,ymax=5.65,
 xtick={-4,-2,0,2,4},ytick={-4,-2,0,2,4},axis equal image,
 legend style={at={(0.5,-0.21)},anchor=north},
]
"""
    tex += curve(circle(2), "framemagenta,line width=1.2pt,fill=framemagenta!6")
    tex += r"\addlegendentry{Horizon $r=2Gm$}" + "\n"
    for radius, label, angle in ((3, r"$\lambda=1/3$", pi/3),
                                 (4, r"$\lambda=1/4$", 3*pi/4),
                                 (5, r"$\lambda=1/5$", pi/4)):
        tex += curve(circle(radius), "metriccyan")
        tex += note(radius*cos(angle), radius*sin(angle), label)
    tex += note(0, 0, r"\shortstack{Inside the\\horizon}", "fill=framemagenta!6")
    tex += r"""\nextgroupplot[
 title={Radial light cones},
 xlabel={$r/(Gm)$},ylabel={$t/(Gm)$},
 xmin=0.3,xmax=5.6,ymin=0,ymax=5.25,
 xtick={1,2,3,4,5},ytick={0,1,2,3,4,5},axis equal image,
 legend style={at={(0.5,-0.21)},anchor=north},
]
\fill[framemagenta!6] (axis cs:0.3,0) rectangle (axis cs:2,5.25);
"""
    tex += line((2, 0), (2, 5.25), "framemagenta,line width=1.2pt")
    for t in (0.45, 2.05, 3.65):
        for r in (1.1, 2, 3.45, 4.95):
            tex += cone(r, t, -1, (r-2)/(r+2), 0.62)
    tex += r"\addlegendimage{metriccyan,-{Stealth[length=3pt]}}" + "\n"
    tex += "\\addlegendentry{Local light directions}\n"
    # Exact solution of dr/dt = -1/(r+1), starting at r=3 when t=0.
    observer = lambda t: -1 + sqrt(16 - 2*t)
    tex += curve([(observer(t), t) for t in samples(0, 5.2)], "framemagenta,dashed,line width=1.2pt")
    tex += r"\addlegendentry{Frame observer: $d\boldsymbol r'=0$}" + "\n"
    for t in (1.1, 4.4):
        tex += line((observer(t), t), (observer(t+0.2), t+0.2),
                    "framemagenta,line width=1.2pt,-{Stealth[length=4pt]}")
    tex += note(2.12, 4.94, "Horizon", "anchor=west")
    return tex


def stationary_surface():
    points = []
    for t in samples(0, 2*pi, 481):
        epsilon = 1 + sqrt(1 - SPIN**2*cos(t)**2)
        points.append((sqrt(epsilon**2 + SPIN**2)*sin(t), epsilon*cos(t)))
    return points


def kerr_fields(x, y, z):
    radius2 = x*x + y*y + z*z
    delta = radius2 - SPIN**2
    epsilon2 = (delta + sqrt(delta*delta + 4*SPIN**2*z*z))/2
    epsilon = sqrt(epsilon2)
    lam = epsilon**3 / (epsilon**4 + SPIN**2*z*z)
    denom = epsilon2 + SPIN**2
    u = ((epsilon*x + SPIN*y)/denom,
         (epsilon*y - SPIN*x)/denom,
         (epsilon*z + SPIN**2*z/epsilon)/denom)
    return epsilon, lam, u


def kerr():
    tex = r"""\nextgroupplot[
 title={Kerr: meridional slice},
 xlabel={$r_1/(Gm)$},ylabel={$r_3/(Gm)$},
 xmin=-3.65,xmax=3.65,ymin=-3.65,ymax=3.65,
 xtick={-3,-2,-1,0,1,2,3},ytick={-3,-2,-1,0,1,2,3},axis equal image,
 legend style={at={(0.5,-0.21)},anchor=north},
]
"""
    tex += curve(stationary_surface(), "metriccyan,fill=metriccyan!30,line width=1.8pt")
    tex += r"\addlegendentry{Outer stationary limit: $g_{tt}=0$}" + "\n"
    tex += curve(ellipse(HORIZON), "framemagenta,fill=framemagenta!6,line width=1.2pt")
    tex += r"\addlegendentry{Horizon: $\epsilon/(Gm)=1.6$}" + "\n"
    for epsilon, theta in ((2.4, 0.65), (3.2, -0.65)):
        tex += curve(ellipse(epsilon), "metriccyan,dashed,forget plot")
        tex += note(sqrt(epsilon**2+SPIN**2)*sin(theta), epsilon*cos(theta),
                    rf"$\epsilon={epsilon}\,Gm$")
    tex += line((0, 2.5), (0, 3.35), "framemagenta,-{Stealth[length=4pt]}")
    tex += note(0.15, 2.95, r"$\boldsymbol\omega$", "anchor=west")
    tex += note(0, -0.3, r"\shortstack{Inside the\\horizon}", "fill=framemagenta!6")
    tex += r"""\nextgroupplot[
 title={Kerr: equatorial slice},
 xlabel={$r_1/(Gm)$},ylabel={$r_2/(Gm)$},
 xmin=-3.65,xmax=3.65,ymin=-3.65,ymax=3.65,
 xtick={-3,-2,-1,0,1,2,3},ytick={-3,-2,-1,0,1,2,3},axis equal image,
 legend style={at={(0.5,-0.21)},anchor=north},
]
"""
    tex += curve(circle(sqrt(4 + SPIN**2)), "metriccyan,fill=metriccyan!30,line width=1.8pt")
    tex += r"\addlegendentry{Ergoregion (cyan)}" + "\n"
    tex += curve(circle(sqrt(HORIZON**2 + SPIN**2)), "framemagenta,fill=framemagenta!6,line width=1.2pt,forget plot")
    tex += curve(circle(sqrt(9 + SPIN**2)), "metriccyan,dashed,forget plot")
    tex += note(-2.28, 2.1, r"$\lambda=1/3$")
    # Equal-length arrows show directions only, for the paper's frame observers.
    for radius in (2.0, 2.65, 3.3):
        for theta in samples(0, 2*pi, 13)[:-1]:
            x, y = radius*cos(theta), radius*sin(theta)
            _, lam, u = kerr_fields(x, y, 0)
            wx, wy = (-lam*u[0]/(1+lam), -lam*u[1]/(1+lam))
            norm = hypot(wx, wy)
            tex += line((x, y), (x+0.34*wx/norm, y+0.34*wy/norm),
                        "framemagenta,line width=0.9pt,-{Stealth[length=3.5pt]}")
    tex += r"\addlegendimage{framemagenta,-{Stealth[length=3.5pt]}}" + "\n"
    tex += r"\addlegendentry{Frame observer directions}" + "\n"
    tex += note(0, 0, r"\shortstack{Spin points\\out of the page}", "fill=framemagenta!6")
    return tex


def verify():
    """Check the plotted directions and surfaces against the paper's metric."""
    count = 0

    def near(a, b):
        nonlocal count
        assert isclose(a, b, rel_tol=1e-9, abs_tol=1e-9), (a, b)
        count += 1

    for t in samples(0, 2):
        lam = 1+t
        near(1-lam*lam*(1/lam)**2, 0)
        near(lam*(1/lam), 1)  # Momentum redshift.
    for r in samples(0.31, 5.6):
        for w in (-1, (r-2)/(r+2)):
            near(1-w*w-2/r*(1+w)**2, 0)
        w = -1/(r+1)
        near(1-w*w-2/r*(1+w)**2, (r/(r+1))**2)
    for x, z in stationary_surface():
        epsilon, lam, u = kerr_fields(x, 0, z)
        near(1-2*lam, 0)
        near(sum(v*v for v in u), 1)
    for x, z in ellipse(HORIZON):
        epsilon, lam, u = kerr_fields(x, 0, z)
        near(epsilon, HORIZON)
        # Horizon generator dr/dt = omega cross r.
        w = (0, OMEGA*x, 0)
        near(1-sum(v*v for v in w)-2*lam*(1+sum(a*b for a,b in zip(u,w)))**2, 0)
    for radius in (2.0, 2.65, 3.3):
        for theta in samples(0, 2*pi, 13)[:-1]:
            _, lam, u = kerr_fields(radius*cos(theta), radius*sin(theta), 0)
            w = [-lam*v/(1+lam) for v in u]
            uw = sum(a*b for a,b in zip(u,w))
            for v, direction in zip(w,u):
                near(v+lam*(1+uw)*direction, 0)  # Local dr' vanishes.
    print(f"{count} metric checks passed for the plotted curves and surfaces.")


def main():
    verify()
    OUT.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="sta-metric-plots-") as work:
        work = Path(work)
        for name, body in (("metric_expansion", expansion()),
                           ("metric_schwarzschild", schwarzschild()),
                           ("metric_kerr", kerr())):
            tex = work / f"{name}.tex"
            tex.write_text(PREAMBLE + body + END)
            result = subprocess.run(
                ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex.name],
                cwd=work, capture_output=True, text=True)
            if result.returncode:
                raise RuntimeError(result.stdout[-7000:])
            shutil.copyfile(work / f"{name}.pdf", OUT / f"{name}.pdf")
            print(f"Wrote {OUT / (name + '.pdf')}")


if __name__ == "__main__":
    main()
