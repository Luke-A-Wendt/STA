#!/usr/bin/env python3
"""Rebuild the observer examples using only black/white vector artwork.

Run: PYTHONDONTWRITEBYTECODE=1 python plot_observer_examples.py
Requires the standard library and the same TeX tools as figure_support.py.
All calculations use natural units. Acceleration and jerk plots use the
dimensionless coordinates printed on their axes, not new paper notation.
"""

from math import cos, cosh, exp, isclose, log, pi, sin, sinh, sqrt

from figure_support import build_figure, curve, line, note, samples


BLACK = "black,line width=1.1pt"
PROPER = BLACK + ",dashed"
LIGHT = "black,densely dotted,line width=1.3pt"
GUIDE = "black,line width=.35pt"
MONO = r"""\pgfplotsset{every axis/.append style={
 grid=none,every axis plot/.append style={color=black},
 tick align=inside,scaled ticks=false}}
"""


def group(options="", rows=1, horizontal_sep=1.55):
    return MONO + rf"""\begin{{groupplot}}[
 group style={{group size=2 by {rows},horizontal sep={horizontal_sep}cm,vertical sep=2.2cm}},
 width=7.5cm,height=6.5cm,{options}]
"""


def dot(x, y, open_mark=False):
    mark = "o,mark options={fill=white}" if open_mark else "*"
    return curve([(x, y)], f"black,only marks,mark={mark},mark size=2pt")


def boost(x, t):
    return 1.25*(x-.6*t), 1.25*(t-.6*x)


def events():
    tex = group("axis equal image,xmin=-1.35,xmax=1.65,ymin=-1,ymax=1.65,"
                "xtick={-1,0,1},ytick={-1,0,1}")
    points = [(0, 1, "A"), (.75, 1.25, "B"), (1, 1, "C"), (1, 0, "D")]
    for primed in (False, True):
        tex += (r"\nextgroupplot[title={Moving frame: $\beta=3/5$},xlabel={$x'$},ylabel={$t'$}]"
                if primed else
                r"\nextgroupplot[title={Unprimed frame},xlabel={$x$},ylabel={$t$}]") + "\n"
        for sign in (-1, 1):
            tex += line((-1.6, sign*-1.6), (1.7, sign*1.7), LIGHT)
        tex += line((-1.35, 0), (1.65, 0), GUIDE)
        tex += line((0, -1), (0, 1.65), GUIDE)
        # Only the timelike segments correspond to clock readings.
        for x, t, label in points:
            if primed:
                x, t = boost(x, t)
            if label in ("A", "B"):
                tex += line((0, 0), (x, t), BLACK if label == "A" else PROPER)
            tex += dot(x, t)
            opts = "anchor=south east" if label == "A" else "anchor=south west"
            tex += note(x+(-.04 if label == "A" else .05), t+.035, f"${label}$", opts)
        tex += dot(0, 0) + note(-.05, -.06, "$O$", "anchor=north east")
        tex += note(-1.15, 1.48, r"$OA:\ \Delta s=1$", "anchor=west")
        tex += note(.65, -.72, r"$OB:\ \Delta s=1$")
    return tex + r"\end{groupplot}" + "\n"


def light_clock():
    tex = group("axis equal image,xmin=-.4,xmax=2.1,ymin=-.4,ymax=1.55,"
                "xtick={0,.75,1.5},ytick={0,1},clip=false")
    tex += r"""\nextgroupplot[title={Rest frame: spatial light path},
 xlabel={$x'/\ell$},ylabel={$r'_2/\ell$}]
"""
    tex += line((0, 0), (0, 1), LIGHT)
    for y in (0, 1):
        tex += line((-.18, y), (.18, y), "black,line width=2pt")
    tex += dot(0, 0) + dot(0, 1)
    tex += note(.15, .06, r"$O,R:\ t'/\ell=0,2$", "anchor=west")
    tex += note(.15, 1.07, r"$M:\ t'/\ell=1$", "anchor=west")
    tex += line((-.18, .2), (-.18, .8), BLACK+",-{Stealth[length=4pt]}")
    tex += line((.18, .8), (.18, .2), BLACK+",-{Stealth[length=4pt]}")
    tex += note(.9, .52, r"$\Delta s=2\ell$")
    tex += r"""\nextgroupplot[title={Unprimed frame: moving clock},
 xlabel={$x/\ell$},ylabel={$r_2/\ell$}]
"""
    # These are spatial projections of the emission, reflection, and return.
    tex += curve([(0, 0), (.75, 1), (1.5, 0)], LIGHT)
    for x, y in [(0, 0), (.75, 1), (1.5, 0)]:
        tex += line((x-.12, y), (x+.12, y), "black,line width=2pt") + dot(x, y)
    tex += note(-.02, -.12, r"$O:\ t=0$", "anchor=north")
    tex += note(.75, 1.1, r"$M:\ t=5\ell/4$", "anchor=south")
    tex += note(1.5, -.12, r"$R:\ t=5\ell/2$", "anchor=north")
    tex += note(1.73, .9, r"$\beta=3/5$")
    tex += note(.75, .22, r"$\Delta t=5\ell/2$")
    return tex + "\\end{groupplot}\n"


def rod():
    tex = group("axis equal image,xmin=-.7,xmax=1.55,ymin=-1,ymax=1.35,"
                "xtick={0,.8,1},ytick={-.75,0,.6,1},clip=true")
    points = [(0, 0, "O"), (1, 0, "A"), (1, .6, "B")]
    for primed in (False, True):
        tex += (r"\nextgroupplot[title={Moving frame: measure at $t'=0$},xlabel={$x'/\ell$},ylabel={$t'/\ell$}]"
                if primed else
                r"\nextgroupplot[title={Rod at rest: measure at $t=0$},xlabel={$x/\ell$},ylabel={$t/\ell$}]") + "\n"
        for x in (0, 1):
            p0, p1 = (x, -1.5), (x, 2)
            if primed:
                p0, p1 = boost(*p0), boost(*p1)
            tex += line(p0, p1, BLACK)
        for x, t in [(1, 0), (1, .6)]:
            end = boost(x, t) if primed else (x, t)
            tex += line((0, 0), end, "black,line width=1.1pt" + (",dash dot" if t else ",dotted"))
        for x, t, label in points:
            if primed:
                x, t = boost(x, t)
            tex += dot(x, t) + note(x+.04, t+.04, f"${label}$", "anchor=south west")
        tex += note(.48 if primed else .48, .89,
                    r"$OB:\ \Delta x'=4\ell/5$" if primed else r"$OA:\ \Delta x=\ell$")
        if not primed:
            tex += note(.44, -.25, r"$\Delta t=0$")
            tex += note(.48, .39, r"$\Delta t'=0$", "rotate=31")
        else:
            tex += note(.42, -.18, r"$\Delta t'=0$")
            tex += note(.79, -.63, r"$\Delta t=0$", "rotate=-31")
    return tex + "\\end{groupplot}\n"


def limits():
    tex = group("xmin=0,xmax=1,xtick={0,.6,1},clip=true")
    tex += r"""\nextgroupplot[title={Clock rate and longitudinal length},
 xlabel={$\beta$},ylabel={Fraction of rest value},ymin=0,ymax=1.18,ytick={0,.5,.8,1}]
"""
    # The clock-rate and length-fraction relations coincide exactly.
    tex += curve([(b, sqrt(1-b*b)) for b in samples(0, 1)], PROPER)
    tex += line((.6, 0), (.6, .8), "black,densely dotted,line width=.5pt")
    tex += line((0, .8), (.6, .8), "black,densely dotted,line width=.5pt")
    tex += dot(.6, .8)
    tex += dot(1, 0, True)
    tex += note(.42, 1.04, r"$\mathrm ds/\mathrm dt=1/\gamma$")
    tex += note(.43, .42, r"Length fraction $=1/\gamma$")
    tex += note(.65, .84, r"$4/5$", "anchor=west")
    tex += r"""\nextgroupplot[title={Energy and invariant mass},
 xlabel={$\beta$},ylabel={$E/m$},ymin=0,ymax=4.2,ytick={0,1,2,3,4}]
"""
    tex += curve([(b, 1/sqrt(1-b*b)) for b in samples(0, .98)], BLACK)
    tex += line((.6, 0), (.6, 1.25), "black,densely dotted,line width=.5pt")
    tex += line((1, 0), (1, 4.2), "black,densely dotted,line width=.6pt")
    tex += dot(.6, 1.25)
    tex += note(.56, 1.65, r"$E/m=5/4$")
    tex += note(.32, 3.25, r"$E^2-\mathbf P^2=m^2$")
    tex += note(.79, 3.65, r"$\gamma\to\infty$")
    return tex + "\\end{groupplot}\n"


def rindler_map(xp, tp):
    """Dimensionless (a*x, a*t), with a=dot v_0 only inside this code."""
    return (1+xp)*cosh(tp)-1, (1+xp)*sinh(tp)


def rindler():
    tex = group(rows=2)
    tex += r"""\nextgroupplot[title={Inertial coordinates},
 axis equal image,xlabel={$\dot v_0x$},ylabel={$\dot v_0t$},
 xmin=-1.2,xmax=2.5,ymin=-.15,ymax=3.55,xtick={-1,0,1,2},ytick={0,1,2,3}]
"""
    for q in (0, .5, 1, 1.5):
        tex += curve([rindler_map(x, q) for x in samples(-1, 2.5)], GUIDE)
    for xp in (-.5, 0, 1):
        tex += curve([rindler_map(xp, q) for q in samples(0, 2.7)],
                     BLACK if xp == 0 else "black,line width=.65pt")
    tex += line((-1, 0), (2.5, 3.5), LIGHT)
    tex += note(.6, 2.03, "horizon", "rotate=45")
    # The same two light rays appear curved in the right panel.
    for sign in (-1, 1):
        tex += line((0, 0), (sign*3.5, 3.5), LIGHT)
    tex += dot(0, 0) + note(.1, .06, "$O$", "anchor=west")
    tex += dot(.25, .75) + note(.33, .72, "$E$", "anchor=west")
    tex += r"""\nextgroupplot[title={Accelerated coordinates},
 xlabel={$\dot v_0x'$},ylabel={$\dot v_0t'$},
 xmin=-1.2,xmax=2.5,ymin=-.15,ymax=2,xtick={-1,0,1,2},ytick={0,.5,1,1.5,2}]
\fill[pattern=north east lines,pattern color=black] (axis cs:-1.2,-.15) rectangle (axis cs:-1,2);
"""
    for xp in (-.5, 0, 1):
        tex += line((xp, 0), (xp, 2), BLACK if xp == 0 else "black,line width=.65pt")
    for q in (0, .5, 1, 1.5):
        tex += line((-1, q), (2.5, q), GUIDE)
    for sign in (-1, 1):
        tex += curve([(exp(sign*q)-1, q) for q in samples(0, 2)], LIGHT)
    for xp in (-.5, .5, 1.5):
        for sign in (-1, 1):
            tex += line((xp, 1.5), (xp+sign*(1+xp)*.13, 1.63), GUIDE)
    tex += dot(0, 0) + note(.05, .06, "$O$", "anchor=west")
    tex += dot(0, log(2)) + note(.08, log(2)+.03, "$E$", "anchor=west")
    tex += note(1.7, 1.78, r"$x'\ \mathrm{constant}$")
    tex += r"""\nextgroupplot[title={Clocks at fixed grid positions},
 xlabel={$\dot v_0t'$},ylabel={$\dot v_0s$},
 xmin=0,xmax=1.55,ymin=0,ymax=3.3,xtick={0,.5,1,1.5},ytick={0,1,2,3}]
"""
    for xp in (-.5, 0, 1):
        tex += curve([(q, (1+xp)*q) for q in samples(0, 1.5)], PROPER)
        tex += note(1.1, .37 if xp == -.5 else (1+xp)*1.1+.18,
                    rf"$\dot v_0x'={xp:g}$")
    tex += r"""\nextgroupplot[title={Proper time between $O$ and $E$},
 axis equal image,xlabel={$\dot v_0x$},ylabel={$\dot v_0t$},
 xmin=-.05,xmax=.85,ymin=-.05,ymax=.85,xtick={0,.25},ytick={0,.75},
 xticklabels={$0$,$1/4$},yticklabels={$0$,$3/4$}]
"""
    tex += line((0, 0), (.25, .75), BLACK)
    tex += curve([rindler_map(0, q) for q in samples(0, log(2))], PROPER)
    tex += dot(0, 0) + dot(.25, .75)
    tex += note(.27, .77, "$E$", "anchor=west")
    tex += note(.025, .035, "$O$", "anchor=west")
    tex += note(.54, .49, r"\shortstack{Accelerated (dashed)\\$\dot v_0\Delta s=\ln(2)$}")
    tex += note(.54, .21, r"\shortstack{Inertial (solid)\\$\dot v_0\Delta s=1/\sqrt{2}$}")
    return tex + "\\end{groupplot}\n"


def integrate(f, upper, n=240):
    """Composite Simpson quadrature; n is even."""
    h = upper/n
    return h/3*(f(0)+f(upper)
                +4*sum(f((2*k+1)*h) for k in range(n//2))
                +2*sum(f(2*k*h) for k in range(1, n//2)))


def jerk_origin(q):
    return integrate(lambda z: sinh(z*z/2), q), integrate(lambda z: cosh(z*z/2), q)


def jerk_map(xp, q):
    x, t = jerk_origin(q)
    return x+xp*cosh(q*q/2), t+xp*sinh(q*q/2)


def jerk_ray(q, sign):
    return sign*exp(sign*q*q/2)*integrate(lambda z: exp(-sign*z*z/2), q)


def jerk():
    tex = group(rows=2, horizontal_sep=2)
    tex += r"""\nextgroupplot[title={Inertial coordinates},
 axis equal image,xlabel={$\sqrt{\ddot v_0}\,x$},ylabel={$\sqrt{\ddot v_0}\,t$},
 xmin=-.8,xmax=2.3,ymin=-.1,ymax=3,xtick={0,1,2},ytick={0,1,2}]
"""
    for q in (0, .5, 1, 1.5):
        tex += curve([jerk_map(xp, q) for xp in samples(-.5, .7, 50)], GUIDE)
    for xp in (-.5, 0, .5):
        tex += curve([jerk_map(xp, q) for q in samples(0, 1.5)], BLACK if xp == 0 else GUIDE)
    for sign in (-1, 1):
        tex += line((0, 0), (sign*3, 3), LIGHT)
    tex += dot(0, 0) + note(.07, .02, "$O$", "anchor=west")
    ex, et = jerk_origin(1)
    tex += dot(ex, et) + note(ex+.07, et, "$E$", "anchor=west")
    tex += note(1.03, 2.71, r"$\sqrt{\ddot v_0}\,t'=1.5$")
    tex += r"""\nextgroupplot[title={Coordinates with constant proper jerk},
 xlabel={$\sqrt{\ddot v_0}\,x'$},ylabel={$\sqrt{\ddot v_0}\,t'$},
 xmin=-1.15,xmax=2.15,ymin=-.1,ymax=1.6,xtick={-1,0,1,2},ytick={0,.5,1,1.5}]
"""
    boundary = [(-1/q, q) for q in samples(1/1.15, 1.6)]
    # Cross-hatching denotes excluded coordinates, never a material object.
    tex += (r"\fill[pattern=north east lines,pattern color=black] (axis cs:-1.15,1.6) -- "
            + " -- ".join(f"(axis cs:{x:.9f},{q:.9f})" for x, q in boundary)
            + " -- cycle;\n")
    tex += curve(boundary, "black,densely dotted,line width=.5pt")
    for xp in (-.5, 0, .5):
        tex += line((xp, 0), (xp, 1.6), BLACK if xp == 0 else GUIDE)
    for sign in (-1, 1):
        tex += curve([(jerk_ray(q, sign), q) for q in samples(0, 1.2)], LIGHT)
    for q in (.1, .65, 1.2):
        for sign in (-1, 1):
            tex += line((1.4, q), (1.4+sign*(1+q*1.4)*.1, q+.1), GUIDE)
    tex += dot(0, 0) + note(.07, .03, "$O$", "anchor=west")
    tex += dot(0, 1) + note(.07, 1.04, "$E$", "anchor=west")
    tex += note(1.35, 1.48, "local cones")
    tex += r"""\nextgroupplot[title={Proper acceleration and jerk},
 xlabel={$\sqrt{\ddot v_0}\,t'$},ylabel={Normalized proper derivatives},
 xmin=0,xmax=1.55,ymin=0,ymax=1.85,xtick={0,.5,1,1.5},ytick={0,.5,1,1.5}]
"""
    tex += curve([(q, q) for q in samples(0, 1.5)], PROPER)
    tex += curve([(0, 1), (1.5, 1)], PROPER)
    tex += note(.78, 1.18, r"$(\mathrm d^2\theta/\mathrm dt'^2)/\ddot v_0=1$")
    tex += note(.8, .34, r"$(\mathrm d\theta/\mathrm dt')/\sqrt{\ddot v_0}$")
    tex += r"""\nextgroupplot[title={Changing clock-rate differences},
 xlabel={$\sqrt{\ddot v_0}\,t'$},ylabel={$\mathrm ds/\mathrm dt'$},
 xmin=0,xmax=1.55,ymin=0,ymax=2,xtick={0,.5,1,1.5},ytick={0,.5,1,1.5,2}]
"""
    for xp in (-.5, 0, .5):
        tex += curve([(q, 1+xp*q) for q in samples(0, 1.5)], PROPER)
        tex += note(.92, 1+xp*.92+.08, rf"$\sqrt{{\ddot v_0}}\,x'={xp:g}$", "anchor=south")
    return tex + "\\end{groupplot}\n"


def rotating():
    tex = group()
    tex += r"""\nextgroupplot[title={Coordinate light velocities},
 axis equal image,xlabel={$\mathrm dr'_1/\mathrm dt'$},ylabel={$\mathrm dr'_2/\mathrm dt'$},
 xmin=-1.25,xmax=1.25,ymin=-1.85,ymax=1.35,xtick={-1,0,1},ytick={-1.6,-1,0,1}]
"""
    for b in (0, .6):
        tex += curve([(cos(p), sin(p)-b) for p in samples(0, 2*pi)],
                     LIGHT if b else "black,dash dot,line width=1pt")
        tex += dot(0, -b, True)
    tex += line((-1.25, 0), (1.25, 0), GUIDE)
    tex += line((0, -1.8), (0, 1.3), GUIDE)
    tex += note(.05, 1.11, r"$\beta=0$", "anchor=south west")
    tex += note(.05, -1.7, r"$\beta=3/5$", "anchor=west")
    tex += note(.03, .51, r"$+2/5$", "anchor=west")
    tex += dot(0, .4) + dot(0, -1.6)
    tex += r"""\nextgroupplot[title={Clock rate on a rotating grid},
 xlabel={Grid speed $\beta$},ylabel={$\mathrm ds/\mathrm dt'$},
 xmin=0,xmax=1,ymin=0,ymax=1.18,xtick={0,.6,1},ytick={0,.5,.8,1}]
"""
    tex += curve([(b, sqrt(1-b*b)) for b in samples(0, 1)], PROPER)
    tex += dot(.6, .8) + dot(1, 0, True)
    tex += line((.6, 0), (.6, .8), "black,densely dotted,line width=.5pt")
    tex += note(.33, .55, r"$\sqrt{1-\beta^2}$")
    tex += note(.68, .87, r"$4/5$")
    tex += note(.67, .16, "null boundary", "anchor=east")
    tex += line((.7, .13), (.97, .025), "black,line width=.5pt,-{Stealth[length=3pt]}")
    return tex + "\\end{groupplot}\n"


def verify():
    """Check event mapping, null rays, proper times, and numerical quadrature."""
    for x, t, expected in [(0, 1, (-.75, 1.25)), (.75, 1.25, (0, 1)),
                           (1, 1, (.5, .5)), (1, 0, (1.25, -.75)),
                           (1, .6, (.8, 0))]:
        xp, tp = boost(x, t)
        assert all(isclose(a, b, abs_tol=1e-14) for a, b in zip((xp, tp), expected))
        assert isclose(t*t-x*x, tp*tp-xp*xp, abs_tol=1e-14)
    # Moving light clock: both spatial legs have length equal to elapsed time.
    assert isclose(sqrt(.75**2+1), 1.25)
    assert all(isclose(a, b) for a, b in zip(rindler_map(0, log(2)), (.25, .75)))
    assert log(2) < 1/sqrt(2)
    for q in samples(0, 1.2, 15):
        for sign in (-1, 1):
            x, t = rindler_map(exp(sign*q)-1, q)
            assert isclose(x, sign*t, abs_tol=1e-13)
            xp = jerk_ray(q, sign)
            assert 1+q*xp > 0  # Every plotted ray point lies inside the chart.
            x, t = jerk_map(xp, q)
            assert isclose(x, sign*t, abs_tol=3e-10)
        # Independent derivative check of the moving-origin map and metric.
        h = 1e-5
        for xp in (-.5, 0, .5):
            x0, t0 = jerk_map(xp, q-h)
            x1, t1 = jerk_map(xp, q+h)
            dx, dt = (x1-x0)/(2*h), (t1-t0)/(2*h)
            assert isclose(dt*dt-dx*dx, (1+q*xp)**2, abs_tol=2e-8)
    fine = integrate(lambda z: cosh(z*z/2), 1.5, n=480)
    assert abs(jerk_origin(1.5)[1]-fine) < 2e-9
    assert isclose(sqrt(1-.6**2), .8)
    for angle in samples(0, 2*pi, 31):
        vx, vy = cos(angle), sin(angle)-.6
        assert isclose(vx*vx+(vy+.6)**2, 1)
    print("Verified inertial events, light clock, rod measurements, accelerated/null maps, and rotating cones.")


if __name__ == "__main__":
    verify()
    for name, drawing in [
        ("observer_events", events), ("observer_light_clock", light_clock),
        ("observer_rod", rod), ("observer_limits", limits),
        ("observer_rindler", rindler), ("observer_jerk", jerk),
        ("observer_rotation", rotating),
    ]:
        build_figure(name, drawing())
