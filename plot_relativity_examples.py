#!/usr/bin/env python3
"""Regenerate the spacetime and midpoint-braking figures with stdlib and PGFPlots.

Run: python plot_relativity_examples.py
Natural units c=1. All figures use black, gray, and white.
"""

from math import asinh, isclose, sqrt

from figure_support import build_figure, curve, line, note, samples


def lorentz_event(x, t, beta=0.6):
    """Return (x', t') for the document's passive, positive-velocity boost."""
    gamma = 1 / sqrt(1-beta*beta)
    return gamma*(x-beta*t), gamma*(t-beta*x)


def verify_spacetime_examples():
    # Exact textbook values and interval checks, independent of drawing styles.
    from fractions import Fraction as F
    beta, gamma = F(3, 5), F(5, 4)
    def transform(x, t):
        return gamma*(x-beta*t), gamma*(t-beta*x)
    cases = [((3, 5), (0, 4)), ((4, 4), (2, 2)),
             ((0, 4), (-3, 5)), ((4, F(32, 5)), (F(1, 5), 5)),
             ((0, 10), (F(-15, 2), F(25, 2)))]
    for (x, t), expected in cases:
        xp, tp = transform(x, t)
        assert (xp, tp) == expected
        assert tp*tp-xp*xp == t*t-x*x
    # Same-t events P,Q are spacelike and differ by -3 in primed time.
    p, q = transform(0, 4), transform(4, 4)
    assert q[1]-p[1] == -3
    assert (q[1]-p[1])**2-(q[0]-p[0])**2 == -16
    # Each leg of the traveling clock accumulates 4, in either frame.
    o, turn, reunion = (0, 0), (3, 5), (0, 10)
    for first, last in [(o, turn), (turn, reunion)]:
        dx, dt = last[0]-first[0], last[1]-first[1]
        assert dt*dt-dx*dx == 16
        xp0, tp0 = transform(*first)
        xp1, tp1 = transform(*last)
        assert (tp1-tp0)**2-(xp1-xp0)**2 == 16
    for s in range(5):
        assert transform(F(3, 4)*s, F(5, 4)*s) == (0, s)
    print('Verified event coordinates, simultaneity, and proper times in both frames.')


def event_dot(x, t):
    return curve([(x, t)], "black,only marks,mark=*,mark size=2pt,forget plot")


def boost_geometry():
    """One clock and one light signal, drawn in two inertial coordinate charts."""
    verify_spacetime_examples()
    tex = r"""\begin{groupplot}[
 group style={group size=2 by 1,horizontal sep=1.3cm},
 width=7.8cm,height=6.8cm,axis equal image,
 xmin=-4.5,xmax=5.5,ymin=-.5,ymax=6.5,
 xtick={-4,-2,0,2,4},ytick={0,1,2,3,4,5,6},
 grid=none,clip=false]
\nextgroupplot[title={Lab frame},xlabel={$x$},ylabel={$t$}]
"""
    for primed in (False, True):
        if primed:
            tex += r"\nextgroupplot[title={Moving clock's rest frame},xlabel={$x'$},ylabel={$t'$}]"+'\n'
        # A 1+1-dimensional slice of the future cone; equal unit scales.
        tex += r"\fill[black!4] (axis cs:0,0) -- (axis cs:-4.5,4.5) -- (axis cs:-4.5,6.5) -- (axis cs:5.5,6.5) -- (axis cs:5.5,5.5) -- cycle;"+'\n'
        for sign in (-1, 1):
            tex += line((0, 0), (sign*(4.5 if sign<0 else 5.5), 4.5 if sign<0 else 5.5), "black,densely dotted,line width=.8pt")
        # Reference lab origin and the clock worldline share the event O.
        tex += line((0, 0), (-3.6, 6) if primed else (0, 6), "black,line width=1pt")
        points = [(0.75*s, 1.25*s) for s in range(5)]
        if primed:
            points = [lorentz_event(*p) for p in points]
        tex += curve(points, "black,dashed,line width=1.6pt")
        for s, (x, t) in enumerate(points[1:-1], 1):
            tex += curve([(x, t)], "black,only marks,mark=o,mark options={fill=white},mark size=1.7pt,forget plot")
            tex += note(x+.15 if primed else x-.15, t+.06, f'${s}$',
                        'anchor=west' if primed else 'anchor=east')
        a = points[-1]
        light = lorentz_event(4, 4) if primed else (4, 4)
        tex += event_dot(0, 0)+event_dot(*a)+event_dot(*light)
        tex += note(-.18, -.08, '$O$', 'anchor=north east,fill=none')
        tex += note(a[0]+.15, a[1]+.15, r"$A:\ s=4$", 'anchor=south west')
        tex += note(light[0]+.15, light[1], '$L$', 'anchor=west')
        tex += note(2.5 if primed else -2.5, 5.7, 'future light cone', 'fill=none')
    tex += r"""\end{groupplot}
\node[anchor=north,align=center,font=\small] at
 ($(group c1r1.south)!0.5!(group c2r1.south)+(0,-1.05cm)$)
 {Same events: $A:(x,t)=(3,5)\ \longleftrightarrow\ (x',t')=(0,4)$\\[3pt]
 Light signal: $L:(4,4)\ \longleftrightarrow\ (2,2)$;\quad $O:(0,0)$ in both frames.};
"""
    return tex


def spacetime_examples():
    """Relativity of simultaneity and the path dependence of proper time."""
    verify_spacetime_examples()
    tex = r"""\begin{groupplot}[
 group style={group size=2 by 1,horizontal sep=1.8cm},
 width=7.6cm,height=8.5cm,axis equal image,grid=none,clip=false]
\nextgroupplot[title={Which events are simultaneous?},
 xlabel={$x$},ylabel={$t$},xmin=-1.4,xmax=5.4,ymin=-.5,ymax=7.4,
 xtick={0,2,4},ytick={0,2,4,6}]
\fill[black!4] (axis cs:0,0) -- (axis cs:-1.4,1.4) -- (axis cs:-1.4,7.4) -- (axis cs:5.4,7.4) -- (axis cs:5.4,5.4) -- cycle;
"""
    for sign in (-1, 1):
        tex += line((0, 0), (sign*(1.4 if sign<0 else 5.4), 1.4 if sign<0 else 5.4), 'black,densely dotted,line width=.8pt')
    tex += line((0, 0), (4.2, 7), 'black,line width=1pt,-{Stealth[length=3pt]}')
    tex += line((0, 0), (5, 3), 'black,line width=1pt,-{Stealth[length=3pt]}')
    tex += note(4.25, 7.06, "$t'$ axis", 'anchor=south')
    tex += note(4.85, 2.75, "$x'$ axis", 'anchor=north')
    tex += line((-1, 4), (5, 4), 'black,line width=.8pt')
    tex += line((-1, 3.4), (5, 7), 'black,line width=.8pt')
    tex += note(1.7, 4.06, '$t=4$', 'anchor=south')
    tex += note(1.8, 5.18, "$t'=5$", 'anchor=south,rotate=31')
    for (x, t), label, opts in [((0, 0), '$O$', 'anchor=north east'),
                               ((0, 4), '$P$', 'anchor=south east'),
                               ((4, 4), '$Q$', 'anchor=south east'),
                               ((4, 6.4), '$R$', 'anchor=south east')]:
        tex += event_dot(x, t)+note(x-.07, t+.08, label, opts)
    tex += r"""\nextgroupplot[title={Reunited clocks},
 xlabel={$x$},ylabel={$t$},xmin=-4.2,xmax=4.2,ymin=-.5,ymax=10.8,
 xtick={-4,-2,0,2,4},ytick={0,2,4,6,8,10}]
\fill[black!4] (axis cs:0,0) -- (axis cs:4.2,4.2) -- (axis cs:4.2,5.8) -- (axis cs:0,10) -- (axis cs:-4.2,5.8) -- (axis cs:-4.2,4.2) -- cycle;
"""
    for sign in (-1, 1):
        tex += line((0, 0), (sign*4.2, 4.2), 'black,densely dotted,line width=.8pt')
        tex += line((0, 10), (sign*4.2, 5.8), 'black,densely dotted,line width=.8pt')
    tex += line((0, 0), (0, 10), 'black,line width=1.3pt')
    tex += curve([(0, 0), (3, 5), (0, 10)], 'black,dashed,line width=1.5pt')
    for x, t, s in [(1.5, 2.5, 2), (1.5, 7.5, 6)]:
        tex += curve([(x, t)], 'black,only marks,mark=o,mark options={fill=white},mark size=1.7pt,forget plot')
        tex += note(x+.12, t, f'$s={s}$', 'anchor=west')
    tex += event_dot(0, 0)+note(-.15, -.05, '$O$', 'anchor=north east')
    tex += event_dot(3, 5)+note(3, 5.3, r'$T:\ s=4$', 'anchor=south east')
    tex += event_dot(0, 10)+note(0, 10.2, '$D$', 'anchor=south')
    tex += note(-.2, 7, r'\shortstack{Stay:\\$s=t$}', 'anchor=east')
    tex += note(-.2, 9.4, '$s=10$', 'anchor=east')
    tex += note(.5, 9.4, '$s=8$', 'anchor=west')
    tex += note(2.5, 1.5, 'Travel', 'anchor=west')
    tex += r"\end{groupplot}"+'\n'
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
\nextgroupplot[grid style={black!15},
 title={One trip; thrust reverses halfway},
 xlabel={$\dot v_0 x$},ylabel={$\dot v_0 t$},
 xmin=-0.12,xmax=2.12,ymin=-0.15,ymax=3.75,
 xtick={0,1,2},ytick={0,1,2,3},
]
"""
    tex += line((0, tm), (2, tm), "black!45,densely dotted")
    tex += curve([(state(q)[0], tm*q) for q in samples(0, 1)], "black,line width=1.4pt")
    tex += curve([(state(q)[0], tm*q) for q in samples(1, 2)], "black,line width=1.4pt")
    tex += curve([(1, tm)], "black,only marks,mark=*,mark size=2pt")
    tex += note(.34, .45, "Accelerate", "text=black,anchor=west")
    tex += note(.28, 2.95, "Brake", "anchor=west")
    tex += note(.06, tm+.26, "Reversal", "anchor=west")
    tex += note(1.93, 3.58, "Arrive at rest", "anchor=east")
    tex += note(.12, .05, "Start at rest", "anchor=west")
    tex += r"""\nextgroupplot[grid style={black!15},
 title={Speed stays continuous},
 xlabel={$t/t_m$},ylabel={Speed or clock rate},
 xmin=0,xmax=2,ymin=0,ymax=1.08,
 xtick={0,1,2},ytick={0,0.5,1},
 legend style={at={(0.5,-0.24)},anchor=north,legend columns=2},
]
"""
    tex += line((1, 0), (1, 1.08), "black!45,densely dotted")
    tex += curve([(q, state(q)[1]) for q in samples(0, 2)], "black,line width=1.4pt")
    tex += r"\addlegendentry{$v$}"+"\n"
    tex += curve([(q, state(q)[2]) for q in samples(0, 2)], "black,dashed,line width=1.4pt")
    tex += r"\addlegendentry{$\mathrm ds/\mathrm dt$}"+"\n"
    tex += note(1, .98, r"$v_{\max}=\sqrt{3}/2$", "text=black")
    tex += note(1, .37, r"$\mathrm ds/\mathrm dt=1/2$")
    tex += r"""\nextgroupplot[grid style={black!15},
 title={Proper and coordinate acceleration},
 xlabel={$t/t_m$},ylabel={Acceleration divided by $\dot v_0$},
 xmin=0,xmax=2,ymin=-1.28,ymax=1.28,
 xtick={0,1,2},ytick={-1,0,1},
 legend style={at={(0.5,-0.24)},anchor=north},
]
"""
    tex += line((1, -1.25), (1, 1.25), "black!45,densely dotted")
    tex += curve([(0, 1), (1, 1)], "black,dashed,line width=1.4pt")
    tex += r"\addlegendentry{Proper: $\mathrm d\theta/\mathrm ds$}"+"\n"
    tex += curve([(q, state(q)[2]**3) for q in samples(0, 1)], "black,line width=1.4pt")
    tex += r"\addlegendentry{Coordinate: $\mathrm dv/\mathrm dt$}"+"\n"
    tex += curve([(1, -1), (2, -1)], "black,dashed,line width=1.4pt,forget plot")
    tex += curve([(q, -state(q)[2]**3) for q in samples(1, 2)], "black,line width=1.4pt,forget plot")
    tex += curve([(1, 1), (1, -1)], "black,only marks,mark=o,mark options={fill=white},mark size=2pt,forget plot")
    tex += curve([(1, .125), (1, -.125)], "black,only marks,mark=o,mark options={fill=white},mark size=2pt,forget plot")
    tex += r"""\nextgroupplot[grid style={black!15},
 title={The traveling clock records less time},
 xlabel={$t/t_m$},ylabel={Elapsed time multiplied by $\dot v_0$},
 xmin=0,xmax=2.08,ymin=0,ymax=3.8,
 xtick={0,1,2},ytick={0,1,2,3},
 legend style={at={(0.5,-0.24)},anchor=north},
]
"""
    tex += line((1, 0), (1, 3.8), "black!45,densely dotted")
    tex += curve([(q, tm*q) for q in samples(0, 2)], "black,line width=1.4pt")
    tex += r"\addlegendentry{Departure-frame clocks: $t$}"+"\n"
    tex += curve([(q, state(q)[3]) for q in samples(0, 2)], "black,dashed,line width=1.4pt")
    tex += r"\addlegendentry{Traveling clock: $s$}"+"\n"
    tex += note(1.98, 3.59, r"$3.464$", "anchor=east,text=black")
    tex += note(1.98, 2.38, r"$2.634$", "anchor=east")
    tex += "\\end{groupplot}\n"
    return tex


if __name__ == "__main__":
    build_figure("boost_geometry", boost_geometry())
    build_figure("spacetime_examples", spacetime_examples())
    build_figure("rocket_trip", rocket_trip())
