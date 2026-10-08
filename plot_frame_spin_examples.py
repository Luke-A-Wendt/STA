"""Geometric illustrations for the Frenet and spin sections.

All plotted geometry and probabilities are computed with the standard library.
Run from any directory; figure_support supplies the cyan, black, and gray style.
"""

from math import cos, sin, sqrt, pi, isclose

from figure_support import build_figure, curve, line, note, samples


def frenet_figure():
    # At the vertex of y=x^2/2, k=1 and the osculating circle has
    # centre (0,1), radius 1, and the same tangent and second derivative.
    # The plot is a planar example (tau=0); the frame formula is general.
    assert isclose(1 / (1 + 0**2)**1.5, 1)
    centre, radius = (0, 1), 1
    body = r"""
\begin{axis}[
 at={(0,0)},anchor=south west,width=7.8cm,height=7.2cm,
 axis equal image,axis lines=none,grid=none,
 xmin=-1.7,xmax=1.95,ymin=-.7,ymax=2.65,
 title={The path sets the local axes},clip=false]
"""
    body += curve([(x, x*x/2) for x in samples(-1.6, 1.6)],
                  "black,line width=1.7pt")
    body += curve([(centre[0]+radius*cos(q), centre[1]+radius*sin(q))
                   for q in samples(0, 2*pi)], "black,dashed,line width=.9pt")
    body += line((0, 0), centre, "black,densely dotted,line width=1pt")
    body += r"\fill[black] (axis cs:0,1) circle (1.7pt);" + "\n"
    body += note(-.06, 1.11, "centre", "anchor=east")
    body += note(-.08, .52, r"$1/k$", "anchor=east")
    body += note(-1.23, 1.79, "osculating circle", "anchor=south west")
    body += line((0, 0), (.76, 0), "black,-{Latex[length=2.5mm]},line width=1.5pt")
    body += note(.76, -.07, r"$\bm u$ (tangent)", "anchor=north")
    body += line((0, 0), (0, .76), "black,-{Latex[length=2.5mm]},line width=1.5pt")
    body += note(.075, .72, r"$\bm u'$ (normal)", "anchor=west")
    body += r"""
\draw[black,line width=1.2pt,fill=white] (axis cs:0,0) circle (3.5pt);
\fill[black] (axis cs:0,0) circle (1.2pt);
"""
    body += note(-.1, -.12, r"$\bm u''$ out of the page", "anchor=north east")
    # A second tangent arrow demonstrates that the frame follows the path.
    q = -1.05
    tangent = (1/sqrt(1+q*q), q/sqrt(1+q*q))
    normal = (-q/sqrt(1+q*q), 1/sqrt(1+q*q))
    assert isclose(tangent[0]*normal[0]+tangent[1]*normal[1], 0, abs_tol=1e-15)
    start = (q, q*q/2)
    body += line(start, (start[0]+.52*tangent[0], start[1]+.52*tangent[1]),
                 "black,-{Latex[length=1.9mm]},line width=1pt")
    body += line(start, (start[0]+.52*normal[0], start[1]+.52*normal[1]),
                 "black,-{Latex[length=1.9mm]},line width=1pt")
    body += note(1.15, 1.38, "path", "text=black,anchor=west")
    body += note(.08, -.65, r"Planar example: $\tau=0$", "anchor=south")
    body += r"""
\end{axis}
\begin{axis}[
 at={(8.05cm,0)},anchor=south west,width=7.8cm,height=7.2cm,
 axis equal image,axis lines=none,grid=none,
 xmin=-.65,xmax=2.45,ymin=-.7,ymax=2.65,
 title={Acceleration: speed and turning},clip=false]
"""
    origin, tangential, normal_end, total = (0, 0), (1.6, 0), (0, 1.6), (1.6, 1.6)
    body += line(tangential, total, "black,dashed,line width=.8pt")
    body += line(normal_end, total, "black,dashed,line width=.8pt")
    body += line(origin, tangential, "black,-{Latex[length=2.7mm]},line width=1.5pt")
    body += line(origin, normal_end, "black,-{Latex[length=2.7mm]},line width=1.5pt")
    body += line(origin, total, "black,-{Latex[length=2.8mm]},line width=1.8pt")
    body += note(.85, -.11, r"$\dot v\,\bm u$", "anchor=north,text=black")
    body += note(.85, -.36, "changes speed", "anchor=north")
    body += note(-.09, 1.72, r"$v^2k\,\bm u'$", "anchor=south,text=black")
    body += note(-.09, 2.03, "changes direction", "anchor=south")
    body += note(1.1, 1.36, r"$\dot{\bm v}$", "anchor=south east")
    body += note(.85, -.66, r"$\dot{\bm v}=\dot v\,\bm u+v^2k\,\bm u'$", "anchor=north")
    body += r"\end{axis}" + "\n"
    return build_figure("frenet_frame", body)


def spin_figure():
    # Prepare + spin along u'; measuring along u depends only on their angle.
    for theta in samples(0, pi):
        plus, minus = cos(theta/2)**2, sin(theta/2)**2
        assert isclose(plus+minus, 1, abs_tol=1e-15)
        assert isclose(plus, (1+cos(theta))/2, abs_tol=1e-15)
    for theta, expected in [(0, 1), (pi/2, .5), (pi, 0)]:
        assert isclose(cos(theta/2)**2, expected, abs_tol=1e-15)

    body = r"""
\begin{axis}[
 at={(0,0)},anchor=south west,width=7.3cm,height=6.9cm,
 axis equal image,axis lines=none,grid=none,
 xmin=-.6,xmax=2.6,ymin=-.7,ymax=2.75,
 title={Preparation and measurement},clip=false]
"""
    body += line((0, 0), (0, 2), "black,-{Latex[length=3mm]},line width=1.7pt")
    body += note(-.04, 2.03, r"$\bm u'$", "anchor=south east,text=black")
    body += note(-.08, 1.68, "prepared", "anchor=east")
    body += line((0, 0), (sqrt(3), 1), "black,-{Latex[length=3mm]},line width=1.7pt")
    body += note(sqrt(3)+.05, 1.0, r"$\bm u$", "anchor=west")
    body += note(sqrt(3), .77, "measured", "anchor=north")
    body += curve([(.62*cos(q), .62*sin(q)) for q in samples(pi/6, pi/2, 61)],
                  "black,-{Latex[length=1.8mm]},line width=.9pt")
    body += note(.45, .75, r"$\theta$")
    body += note(.72, -.24, r"Prepared spin: $+\frac12$ along $\bm u'$", "anchor=north")
    body += note(.72, -.58, r"$\bm u\cdot\bm u'=\cos(\theta)$", "anchor=north")
    body += r"""
\end{axis}
\begin{axis}[
 at={(8.1cm,.2cm)},anchor=south west,width=7.9cm,height=6.5cm,
 xmin=0,xmax=3.14159265359,ymin=0,ymax=1,
 title={Probabilities along the measurement axis},
 xlabel={angle $\theta$},ylabel={probability},
 xtick={0,1.57079632679,3.14159265359},
 xticklabels={$0$,$\pi/2$,$\pi$},ytick={0,.5,1},
 minor tick num=0,grid=none,clip=false]
"""
    body += line((pi/2, 0), (pi/2, 1), "black!40,densely dotted,line width=.6pt")
    body += line((0, .5), (pi, .5), "black!40,densely dotted,line width=.6pt")
    body += curve([(q, cos(q/2)**2) for q in samples(0, pi)], "black,line width=1.6pt")
    body += curve([(q, sin(q/2)**2) for q in samples(0, pi)], "black,dashed,line width=1.6pt")
    body += note(.96, .88, r"$\mathrm{prob}(+)=\cos^2(\theta/2)$", "text=black")
    body += note(.96, .13, r"$\mathrm{prob}(-)=\sin^2(\theta/2)$")
    body += r"""
\addplot[only marks,mark=*,mark size=2.5pt,black] coordinates {(0,1) (1.57079632679,.5) (3.14159265359,0)};
\addplot[only marks,mark=o,mark size=3.5pt,black,line width=1pt] coordinates {(0,0) (1.57079632679,.5) (3.14159265359,1)};
"""
    body += note(0, -.31, "aligned", "anchor=north")
    body += note(pi/2, -.31, "perpendicular", "anchor=north")
    body += note(pi, -.31, "opposite", "anchor=north")
    body += note(0, -.41, r"$\mathrm{prob}(+)=1$", "anchor=north")
    body += note(pi/2, -.41, r"$\mathrm{prob}(\pm)=1/2$", "anchor=north")
    body += note(pi, -.41, r"$\mathrm{prob}(-)=1$", "anchor=north")
    body += r"\end{axis}" + "\n"
    return build_figure("spin_probability", body)


if __name__ == "__main__":
    frenet_figure()
    spin_figure()
    print("Verified the Frenet geometry and spin-probability identities.")
