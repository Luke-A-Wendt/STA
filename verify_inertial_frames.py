"""Check changing-frame acceleration independently of the manuscript.

Run with Python and SymPy. Checks use noncommuting Pauli matrices, direct
coordinate differentiation, and the geodesic equation of the derived metric.
No equations are read from the LaTeX source.
"""

import sympy as s


count = 0


def check(name, expression, simplify=s.cancel):
    global count
    entries = list(expression) if isinstance(expression, s.MatrixBase) else [expression]
    residual = [simplify(s.expand(x)) for x in entries]
    assert all(x == 0 for x in residual), (name, residual)
    count += 1
    print("PASS", name, flush=True)


one = s.eye(2)
sigma = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))


def vector(v):
    return sum((v[i]*sigma[i] for i in range(3)), s.zeros(2))


def scalar(z):
    return s.trace(z)/2


def spatial(z):
    return s.Matrix([scalar(b*z) for b in sigma])


# A genuinely noncommuting determinant-one frame and a changing real event.
q = s.symbols("q", real=True)
frame = s.Matrix([[1, s.I*q], [0, 1]])*s.Matrix([[1, 0], [q*q, 1]])
event = (2*q+q*q)*one + vector(s.Matrix([q**3, q-q*q, 1+q**4]))
omega = s.expand(frame.diff(q)*frame.inv())
event_prime = s.expand(frame*event*frame.H)
check("frame has determinant one and scalar-free generator",
      s.Matrix([frame.det()-1, scalar(omega)]))
check("first derivative includes both changing-frame terms",
      event_prime.diff(q)-frame*event.diff(q)*frame.H
      -omega*event_prime-event_prime*omega.H)
check("second derivative of the congruence retains all ordered products",
      event_prime.diff(q, 2)-frame*event.diff(q, 2)*frame.H
      -2*omega*event_prime.diff(q)-2*event_prime.diff(q)*omega.H
      -(omega.diff(q)-omega**2)*event_prime
      -event_prime*(omega.H.diff(q)-omega.H**2)+2*omega*event_prime*omega.H)

# A variable pointwise boost preserves X^2 but not the naive differential
# interval: this is why the observer coordinate map has to be specified.
pointwise = s.diag(s.exp(q/2), s.exp(-q/2))
image = pointwise*(q*one)*pointwise.H
check("pointwise Lorentz congruence preserves the event determinant", image.det()-q*q)
check("variable congruence is not a differential Lorentz isometry",
      image.diff(q).det()-(1-q*q))

# The six physical components of the generator give momentum transport.
a = s.Matrix(s.symbols("a1:4", real=True))
w = s.Matrix(s.symbols("w1:4", real=True))
ad = s.Matrix(s.symbols("ad1:4", real=True))
wd = s.Matrix(s.symbols("wd1:4", real=True))
r = s.Matrix(s.symbols("r1:4", real=True))
v = s.Matrix(s.symbols("v1:4", real=True))
vd = s.Matrix(s.symbols("vd1:4", real=True))
energy = s.symbols("E", real=True)
momentum = s.Matrix(s.symbols("p1:4", real=True))
generator = vector((-a+s.I*w)/2)
P = energy*one+vector(momentum)
check("generator components give the correct energy and momentum terms",
      generator*P+P*generator.H + a.dot(momentum)*one
      + vector(energy*a+w.cross(momentum)))
check("frame-only momentum transport preserves the invariant mass",
      -2*energy*a.dot(momentum)+2*momentum.dot(energy*a+w.cross(momentum)))

# Independently differentiate the moving-origin coordinate map. The origin
# tangent is fixed by the rest-frame condition; all integrations are polynomials.
inverse = frame.inv()
origin_tangent = s.expand(inverse*inverse.H)
origin = origin_tangent.applyfunc(lambda x: s.integrate(x, q))
grid_path = s.Matrix([q+q*q, 1-q+q**3, q*q-q**3])
mapped_path = s.expand(origin+inverse*vector(grid_path)*inverse.H)
local_first = s.expand(frame*mapped_path.diff(q)*frame.H)
local_second = s.expand(frame*mapped_path.diff(q, 2)*frame.H)
actual_a = s.expand(-spatial(omega+omega.H))
actual_w = s.expand(spatial((omega-omega.H)/s.I))
clock = 1+actual_a.dot(grid_path)
space = grid_path.diff(q)+actual_w.cross(grid_path)
check("event-map differentiation gives the tangent and clock factor",
      local_first-clock*one-vector(space))
time_second = (actual_a.diff(q).dot(grid_path)
               +2*actual_a.dot(grid_path.diff(q))
               +actual_a.dot(actual_w.cross(grid_path)))
space_second = (grid_path.diff(q, 2)+2*actual_w.cross(grid_path.diff(q))
                +actual_w.diff(q).cross(grid_path)
                +actual_w.cross(actual_w.cross(grid_path))+clock*actual_a)
check("twice-differentiated event map has the stated scalar component",
      scalar(local_second)-time_second)
check("twice-differentiated event map gives Euler, Coriolis, and centrifugal terms",
      spatial(local_second)-space_second)

# Compute the coordinate geodesic acceleration from the metric, rather than
# from the frame-transport formula. a and w have arbitrary values and first
# derivatives at the chosen time; the spatial axes are unrestricted.
N = 1+a.dot(r)
shift = w.cross(r)
coframe = s.eye(4)
coframe[0, 0] = N
coframe[1:4, 0] = shift
eta = s.diag(1, -1, -1, -1)
metric = coframe.T*eta*coframe
coframe_inverse = coframe.inv()
metric_inverse = coframe_inverse*eta*coframe_inverse.T
check("observer metric inverse is exact", metric_inverse*metric-s.eye(4))


def time_derivative(expression):
    return sum(expression.diff(a[k])*ad[k]+expression.diff(w[k])*wd[k] for k in range(3))


derivatives = [metric.applyfunc(time_derivative)] + [metric.diff(x) for x in r]
coordinate_tangent = s.Matrix([1, *v])
# Contract the Christoffel connection with (1,v) twice. This is equivalent
# to constructing all 64 connection entries, with much smaller expressions.
lower_connection = s.Matrix([
    sum((derivatives[alpha][nu, beta]
         -derivatives[nu][alpha, beta]/2)*coordinate_tangent[alpha]*coordinate_tangent[beta]
        for alpha in range(4) for beta in range(4))
    for nu in range(4)
])
connection = (metric_inverse*lower_connection).applyfunc(s.cancel)
metric_acceleration = s.Matrix([-connection[k+1]+v[k]*connection[0] for k in range(3)])
numerator = ad.dot(r)+2*a.dot(v)+a.dot(shift)
inertial_acceleration = (-N*a-2*w.cross(v)-wd.cross(r)-w.cross(w.cross(r))
                         +(v+shift)*numerator/N)
check("metric geodesics independently reproduce the complete inertial acceleration",
      metric_acceleration-inertial_acceleration)

# An arbitrary physical four-acceleration is first expressed in the local
# orthonormal axes. Convert it to chart components and eliminate the time
# equation in the coordinate equation of motion.
physical = s.Matrix(s.symbols("A0:4", real=True))
proper_rate_squared = N*N-(v+shift).dot(v+shift)
chart_acceleration = coframe_inverse*physical
source_from_coordinates = proper_rate_squared*s.Matrix([
    chart_acceleration[k+1]-v[k]*chart_acceleration[0] for k in range(3)])
source_from_frame = proper_rate_squared*(physical[1:4, 0]-(v+shift)*physical[0]/N)
check("nonzero physical four-force has the stated coordinate correction",
      source_from_coordinates-source_from_frame)

# A directly differentiated passive rotation checks the orientation signs
# and the Coriolis factor of two without using the observer metric.
angle = q*q/2
rotation = s.Matrix([[s.cos(angle), s.sin(angle), 0],
                     [-s.sin(angle), s.cos(angle), 0], [0, 0, 1]])
straight = s.Matrix([1+2*q, 3-q, 2+q/3])
rotated = rotation*straight
angular_velocity = s.Matrix([0, 0, q])
check("rotating coordinates of a straight worldline have the classical three rotation terms",
      rotated.diff(q, 2)+2*angular_velocity.cross(rotated.diff(q))
      +angular_velocity.diff(q).cross(rotated)
      +angular_velocity.cross(angular_velocity.cross(rotated)))

# Rindler coordinates of an arbitrary inertial axial trajectory. This also
# exercises the non-affine time correction, including relativistic speeds.
alpha = s.symbols("alpha", positive=True)
speed = s.symbols("speed", real=True)
rindler_position = (1/(s.cosh(alpha*q)-speed*s.sinh(alpha*q))-1)/alpha
residual = (rindler_position.diff(q, 2)+alpha*(1+alpha*rindler_position)
            -2*alpha*rindler_position.diff(q)**2/(1+alpha*rindler_position))
check("Rindler inertial-worldline acceleration includes the relativistic time term",
      residual, simplify=s.simplify)
for sign in (-1, 1):
    ray = (s.exp(sign*alpha*q)-1)/alpha
    check(f"Rindler null ray {sign:+d} obeys the same coordinate acceleration",
          ray.diff(q, 2)+alpha*(1+alpha*ray)-2*alpha*ray.diff(q)**2/(1+alpha*ray))

# SI dimensional conversion, with natural time and length both measured in
# metres, physical acceleration a/c^2 and angular rate w/c.
c = s.symbols("c", positive=True)
replacement = dict(zip([*a, *ad, *w, *wd, *v],
                       [*(a/c**2), *(ad/c**3), *(w/c), *(wd/c**2), *(v/c)]))
si_from_natural = c*c*inertial_acceleration.subs(replacement, simultaneous=True)
si_acceleration = (-(1+a.dot(r)/c**2)*a-2*w.cross(v)-wd.cross(r)-w.cross(w.cross(r))
                   +(v+shift)*numerator/(c*c+a.dot(r)))
check("SI acceleration restores every factor of c", si_from_natural-si_acceleration)
classical = -a-2*w.cross(v)-wd.cross(r)-w.cross(w.cross(r))
check("Newtonian limit gives translational, Coriolis, Euler, and centrifugal terms",
      si_acceleration.applyfunc(lambda z: s.limit(z, c, s.oo))-classical)
check("pure rotation needs no relativistic correction to coordinate acceleration",
      (inertial_acceleration-classical).subs(dict.fromkeys([*a, *ad], 0)))
si_P = energy/c*one+vector(momentum)
si_generator = vector((-a/c+s.I*w)/2)
check("SI momentum transport weights translational inertia by E/c squared",
      si_generator*si_P+si_P*si_generator.H
      +a.dot(momentum)/c*one+vector(energy/c**2*a+w.cross(momentum)))

print(f"\n{count} checks passed.")
