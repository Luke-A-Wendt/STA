"""Two explanatory relativity figures, using only black and gray ink.

Run with the standard library and an installed pdflatex.  The figures
illustrate equations already derived in sta_notes.tex; no fitted data
or additional physical assumptions enter the plotted curves.
"""

from math import isclose, log, sqrt

from figure_support import build_figure, curve, line, note, samples


def verify_data():
    """Check invariants and both plotted families of null curves."""
    for beta in samples(0, 0.85):
        gamma = 1 / sqrt(1 - beta * beta)
        electric, magnetic = gamma, -gamma * beta
        assert isclose(electric**2 - magnetic**2, 1, abs_tol=3e-15)
        assert isclose(electric**2 + magnetic**2,
                       (1 + beta * beta) / (1 - beta * beta))
    for time in samples(0, 2):
        scale = 1 + time
        for sign in (-1, 1):
            flat_r = sign * time / scale
            flat_speed = sign / scale**2
            curved_speed = sign / scale
            # Dimensionless t/t0 and r1/t0, with t0 > 0 and c=1.
            assert isclose(1 - (scale * flat_speed + flat_r)**2, 0,
                           abs_tol=1e-14)
            assert isclose(1 - scale**2 * curved_speed**2, 0,
                           abs_tol=1e-14)
            assert isclose((sign - flat_r) / scale, flat_speed)
    print("Verified field invariants and both null-ray families.")


def boosted_fields():
    body = r"""
\begin{groupplot}[group style={group size=2 by 1,horizontal sep=1.2cm},
 height=6.5cm,width=7.5cm]
\nextgroupplot[title={Perpendicular boost},
 xmin=-1,xmax=4,ymin=-1.2,ymax=3.55,hide axis,grid=none,clip=false]
"""
    body += note(0, 3.22, r"Before boost")
    body += note(2.45, 3.22, r"After boost: $\beta=0.6$")
    body += line((0, 0), (0, 2), "plotgray,line width=2pt,-{Stealth[length=3mm]}")
    body += line((2.45, 0), (2.45, 2.5), "plotgray,line width=2pt,-{Stealth[length=3mm]}")
    body += note(0, 2.25, r"$E=E_0$", "text=black")
    body += note(2.45, 2.75, r"$E'=1.25E_0$", "text=black")
    body += note(0, -0.35, r"$\vstyle{B}=0$")
    body += r"""
\draw[black,line width=1.2pt] (axis cs:2.45,0) circle[radius=3.5pt];
\draw[black,line width=1.1pt]
 ([xshift=-2.4pt,yshift=-2.4pt]axis cs:2.45,0) --
 ([xshift=2.4pt,yshift=2.4pt]axis cs:2.45,0);
\draw[black,line width=1.1pt]
 ([xshift=-2.4pt,yshift=2.4pt]axis cs:2.45,0) --
 ([xshift=2.4pt,yshift=-2.4pt]axis cs:2.45,0);
"""
    body += note(2.45, -0.4, r"$B'=0.75E_0$ into page")
    body += line((0.65, 0.8), (1.65, 0.8), "black,line width=1.2pt,-{Stealth}")
    body += note(1.15, 0.43, r"boost along $+r_1$")
    body += note(1.55, -0.95, r"Electric field along $+r_2$; $\vstyle{B}'$ along $-r_3$.")
    body += r"""
\nextgroupplot[title={Energy and invariant},
 xmin=0,xmax=0.85,ymin=0,ymax=7,
 xlabel={Boost speed $\beta$},ylabel={Squared fields divided by $E_0^2$},
 xtick={0,0.2,0.4,0.6,0.8},ytick={0,1,2,3,4,5,6,7},
 legend style={at={(0.04,0.96)},anchor=north west}]
"""
    body += curve([(b, (1+b*b)/(1-b*b)) for b in samples(0, 0.85)],
                  "plotgray,line width=1.5pt")
    body += r"\addlegendentry{$(E'^2+B'^2)/E_0^2$}" + "\n"
    body += curve([(0, 1), (0.85, 1)], "black,dashed,line width=1.3pt")
    body += r"\addlegendentry{$(E'^2-B'^2)/E_0^2=1$}" + "\n"
    body += curve([(0.6, 2.125)], "plotgray,only marks,mark=*,mark size=2pt")
    body += line((0.6, 0), (0.6, 2.125), "black,densely dotted")
    body += note(0.52, 2.7, r"$\beta=0.6:\;2.125$")
    body += "\\end{groupplot}\n"
    return build_figure("boosted_fields", body)


def coordinates_curvature():
    body = r"""
\begin{groupplot}[group style={group size=2 by 1,horizontal sep=1.15cm},
 width=7.5cm,height=7cm,
 xmin=-1.3,xmax=1.3,ymin=0,ymax=2.1,
 xlabel={Chart position $r_1/t_0$},ylabel={Chart time $t/t_0$},
 xtick={-1,-0.5,0,0.5,1},ytick={0,0.5,1,1.5,2}]
\nextgroupplot[title={Flat: exact coordinate change}]
"""
    times = samples(0, 2)
    for sign in (-1, 1):
        body += curve([(sign*time/(1+time), time) for time in times],
                      "plotgray,line width=1.8pt")
    for x in (-0.9, 0.9):
        for time in (0.5, 1.25, 1.85):
            for sign in (-1, 1):
                slope = (sign-x)/(1+time)
                dt = 0.13
                body += line((x-slope*dt/2, time-dt/2),
                             (x+slope*dt/2, time+dt/2),
                             "black,line width=0.8pt,-{Stealth[length=1.2mm]}")
    body += note(0, 1.75, r"$K=0$")
    body += r"""
\nextgroupplot[title={Curved: assigned displacement}]
"""
    for sign in (-1, 1):
        body += curve([(sign*log(1+time), time) for time in times],
                      "plotgray,line width=1.8pt")
    for x in (-0.9, 0.9):
        for time in (0.5, 1.25, 1.85):
            for sign in (-1, 1):
                slope = sign/(1+time)
                dt = 0.13
                body += line((x-slope*dt/2, time-dt/2),
                             (x+slope*dt/2, time+dt/2),
                             "black,line width=0.8pt,-{Stealth[length=1.2mm]}")
    body += note(0, 1.75, r"$K=\dfrac{12}{(t+t_0)^4}>0$")
    body += r"""
\end{groupplot}
\node[font=\small,anchor=north,align=center] at
 ([yshift=-1.15cm]group c1r1.south) {
 $\mathrm ds^2=\mathrm dt^2-\mathrm d\!\left(\lambda\vstyle{r}\right)^2$\\[4pt]
 $\displaystyle\frac{\mathrm dr_1}{\mathrm dt}=\frac{\pm1-r_1/t_0}{\lambda}$};
\node[font=\small,anchor=north,align=center] at
 ([yshift=-1.15cm]group c2r1.south) {
 $\mathrm ds^2=\mathrm dt^2-\lambda^2\mathrm d\vstyle{r}^2$\\[4pt]
 $\displaystyle\frac{\mathrm dr_1}{\mathrm dt}=\frac{\pm1}{\lambda}$};
\node[font=\small,anchor=north,align=center] at
 ([yshift=-2.65cm]$(group c1r1.south)!0.5!(group c2r1.south)$) {
 $\lambda(t)=1+t/t_0,\quad t_0>0$: the same chart and the same scale factor.\\[3pt]
 \textcolor{black}{Thick curves: light rays from the origin.}\quad
 \textcolor{black}{Short arrows: local light directions.}\\[3pt]
 Curved coordinate paths do not prove curvature; the four-dimensional $K$ does.};
"""
    return build_figure("coordinates_curvature", body)


if __name__ == "__main__":
    verify_data()
    boosted_fields()
    coordinates_curvature()
