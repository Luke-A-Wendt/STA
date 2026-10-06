#!/usr/bin/env python3
"""Regenerate the boost and midpoint-braking figures with stdlib and PGFPlots.

Run: python plot_relativity_examples.py
Natural units c=1. Only cyan and magenta ink is used on a white background.
"""

from math import asinh, cos, cosh, isclose, pi, sin, sinh, sqrt, tanh

from figure_support import build_figure, curve, line, note, samples


def boost_geometry():
    theta = 0.8
    beta = tanh(theta)
    space, time = sinh(theta), cosh(theta)
    assert isclose(time*time-space*space, 1)
    assert isclose(time*time+space*space, cosh(2*theta))
    assert isclose(time*(time-beta*space), 1)
    assert abs(time*(space-beta*time)) < 1e-14
    tex = r"""\begin{groupplot}[group style={group size=2 by 1,horizontal sep=1.4cm}]
\nextgroupplot[
 title={Boosted axes and invariant hyperbola},
 xlabel={$x=\boldsymbol u\cdot\boldsymbol r$},ylabel={$t$},
 xmin=-0.12,xmax=2.1,ymin=-0.12,ymax=2.25,
 width=7.65cm,height=7.4cm,axis equal image,
 xtick={0,1,2},ytick={0,1,2},grid=none,
 axis lines=middle,clip=false,
 legend style={at={(0.5,-0.17)},anchor=north},
]
"""
    tex += line((0, 0), (2.05, 2.05), "cyan!55,densely dotted,line width=1.2pt")
    tex += note(1.7, 1.91, "Light", "text=cyan,fill=none,rotate=45")
    tex += curve([(sinh(h), cosh(h)) for h in samples(0, 1.43)], "magenta,line width=1.4pt")
    tex += r"\addlegendentry{$t^2-x^2=1$: determinant}"+"\n"
    tex += curve([(sin(h), cos(h)) for h in samples(0, pi/2)], "cyan,dashed,line width=1.2pt")
    tex += r"\addlegendentry{$t^2+x^2=1$: coefficient norm}"+"\n"
    tex += line((0, 0), (1.2, 1.2/beta), "cyan,line width=0.85pt,-{Stealth[length=4pt]}")
    tex += line((0, 0), (1.93, 1.93*beta), "cyan,line width=0.85pt,-{Stealth[length=4pt]}")
    tex += note(1.17, 2.07, r"$t'$ axis", "text=cyan")
    tex += note(1.79, 1.1, r"$x'$ axis", "text=cyan")
    tex += curve([(0, 1), (space, time)], "magenta,only marks,mark=*,mark size=2pt")
    tex += note(0.14, 1.11, r"$\mathcal U(0)$", "anchor=west")
    tex += note(space+0.12, time+0.12, r"$\mathcal U(\theta)$", "anchor=west")
    tex += note(0.17, 2.1, r"$\theta=0.8$", "anchor=west")
    tex += r"""\nextgroupplot[
 title={What a boost preserves},
 xlabel={Rapidity $\theta$},ylabel={Squared quantity},
 xmin=0,xmax=1.22,ymin=0,ymax=6.15,
 width=7.65cm,height=7.4cm,
 xtick={0,0.4,0.8,1.2},ytick={0,1,2,3,4,5,6},
 legend style={at={(0.5,-0.17)},anchor=north},
]
"""
    tex += curve([(h, 1) for h in (0, 1.2)], "magenta,line width=1.4pt")
    tex += r"\addlegendentry{$\det\mathcal U=1$}"+"\n"
    tex += curve([(h, cosh(2*h)) for h in samples(0, 1.2)], "cyan,dashed,line width=1.4pt")
    tex += r"\addlegendentry{$\lVert\mathcal U\rVert^2=\cosh(2\theta)$}"+"\n"
    tex += line((theta, 0), (theta, cosh(2*theta)), "magenta!45,densely dotted")
    tex += curve([(theta, 1)], "magenta,only marks,mark=*,mark size=2pt")
    tex += curve([(theta, cosh(2*theta))], "cyan,only marks,mark=*,mark size=2pt")
    tex += note(0.44, 4.9, r"$\mathcal U=\cosh\theta+\sinh\theta\,\boldsymbol u\cdot\boldsymbol\sigma$")
    tex += note(0.23, 3.65, r"$\beta=\tanh\theta$", "anchor=west")
    tex += "\\end{groupplot}\n"
    return tex


def rocket_trip():
    # alpha means the existing parameter dot v_0 in calculations only.
    # Set alpha*ell=2, so gamma_m=2 and alpha*t_m=sqrt(3).
    tm, sm = sqrt(3), asinh(sqrt(3))

    def state(q):
        """q=t/t_m; return alpha*x, v, ds/dt, alpha*s."""
        half = q if q <= 1 else 2-q
        rapidity = asinh(tm*half)
        x = cosh(rapidity)-1 if q <= 1 else 3-cosh(rapidity)
        proper_time = rapidity if q <= 1 else 2*sm-rapidity
        return x, tanh(rapidity), 1/cosh(rapidity), proper_time

    for q in samples(0, 2):
        x, velocity, rate, proper_time = state(q)
        assert 0 <= x <= 2 and 0 <= velocity < 1
        assert isclose(rate**2+velocity**2, 1)
        assert proper_time <= tm*q+1e-12
    assert isclose(state(1)[0], 1)
    assert isclose(state(1)[1], sqrt(3)/2)
    assert isclose(state(1)[2], .5)
    assert isclose(state(2)[3], 2*sm)
    for q in (.1, .4, .8, 1.2, 1.6, 1.9):
        step = 1e-5
        before, after = state(q-step), state(q+step)
        _, velocity, rate, _ = state(q)
        sign = 1 if q < 1 else -1
        assert isclose((after[0]-before[0])/(2*step*tm), velocity, abs_tol=1e-9)
        assert isclose((after[3]-before[3])/(2*step*tm), rate, abs_tol=1e-9)
        assert isclose((after[1]-before[1])/(2*step*tm), sign*rate**3, abs_tol=1e-9)
    tex = r"""\begin{groupplot}[
 group style={group size=2 by 2,horizontal sep=1.45cm,vertical sep=2.6cm},
 width=7.65cm,height=6.1cm,
]
\nextgroupplot[
 title={One trip; thrust reverses halfway},
 xlabel={$\dot v_0 x$},ylabel={$\dot v_0 t$},
 xmin=-0.12,xmax=2.12,ymin=-0.15,ymax=3.75,
 xtick={0,1,2},ytick={0,1,2,3},
]
"""
    tex += line((0, tm), (2, tm), "magenta!45,densely dotted")
    tex += curve([(state(q)[0], tm*q) for q in samples(0, 1)], "cyan,line width=1.4pt")
    tex += curve([(state(q)[0], tm*q) for q in samples(1, 2)], "magenta,line width=1.4pt")
    tex += curve([(1, tm)], "magenta,only marks,mark=*,mark size=2pt")
    tex += note(.34, .45, "Accelerate", "text=cyan,anchor=west")
    tex += note(.28, 2.95, "Brake", "anchor=west")
    tex += note(.06, tm+.26, "Reversal", "anchor=west")
    tex += note(1.93, 3.58, "Arrive at rest", "anchor=east")
    tex += note(.12, .05, "Start at rest", "anchor=west")
    tex += r"""\nextgroupplot[
 title={Speed stays continuous},
 xlabel={$t/t_m$},ylabel={Speed or clock rate},
 xmin=0,xmax=2,ymin=0,ymax=1.08,
 xtick={0,1,2},ytick={0,0.5,1},
 legend style={at={(0.5,-0.24)},anchor=north,legend columns=2},
]
"""
    tex += line((1, 0), (1, 1.08), "magenta!45,densely dotted")
    tex += curve([(q, state(q)[1]) for q in samples(0, 2)], "cyan,line width=1.4pt")
    tex += r"\addlegendentry{$v$}"+"\n"
    tex += curve([(q, state(q)[2]) for q in samples(0, 2)], "magenta,line width=1.4pt")
    tex += r"\addlegendentry{$\mathrm ds/\mathrm dt$}"+"\n"
    tex += note(1, .98, r"$v_{\max}=\sqrt{3}/2$", "text=cyan")
    tex += note(1, .37, r"$\mathrm ds/\mathrm dt=1/2$")
    tex += r"""\nextgroupplot[
 title={Proper and coordinate acceleration},
 xlabel={$t/t_m$},ylabel={Acceleration divided by $\dot v_0$},
 xmin=0,xmax=2,ymin=-1.28,ymax=1.28,
 xtick={0,1,2},ytick={-1,0,1},
 legend style={at={(0.5,-0.24)},anchor=north},
]
"""
    tex += line((1, -1.25), (1, 1.25), "magenta!45,densely dotted")
    tex += curve([(0, 1), (1, 1)], "cyan,line width=1.4pt")
    tex += r"\addlegendentry{Proper: $\mathrm d\theta/\mathrm ds$}"+"\n"
    tex += curve([(q, state(q)[2]**3) for q in samples(0, 1)], "magenta,line width=1.4pt")
    tex += r"\addlegendentry{Coordinate: $\mathrm dv/\mathrm dt$}"+"\n"
    tex += curve([(1, -1), (2, -1)], "cyan,line width=1.4pt,forget plot")
    tex += curve([(q, -state(q)[2]**3) for q in samples(1, 2)], "magenta,line width=1.4pt,forget plot")
    tex += curve([(1, 1), (1, -1)], "cyan,only marks,mark=o,mark options={fill=white},mark size=2pt,forget plot")
    tex += curve([(1, .125), (1, -.125)], "magenta,only marks,mark=o,mark options={fill=white},mark size=2pt,forget plot")
    tex += r"""\nextgroupplot[
 title={The traveling clock records less time},
 xlabel={$t/t_m$},ylabel={Elapsed time multiplied by $\dot v_0$},
 xmin=0,xmax=2.08,ymin=0,ymax=3.8,
 xtick={0,1,2},ytick={0,1,2,3},
 legend style={at={(0.5,-0.24)},anchor=north},
]
"""
    tex += line((1, 0), (1, 3.8), "magenta!45,densely dotted")
    tex += curve([(q, tm*q) for q in samples(0, 2)], "cyan,dashed,line width=1.4pt")
    tex += r"\addlegendentry{Departure-frame clocks: $t$}"+"\n"
    tex += curve([(q, state(q)[3]) for q in samples(0, 2)], "magenta,line width=1.4pt")
    tex += r"\addlegendentry{Traveling clock: $s$}"+"\n"
    tex += note(1.98, 3.59, r"$3.464$", "anchor=east,text=cyan")
    tex += note(1.98, 2.38, r"$2.634$", "anchor=east")
    tex += "\\end{groupplot}\n"
    return tex


if __name__ == "__main__":
    build_figure("boost_geometry", boost_geometry())
    build_figure("rocket_trip", rocket_trip())
