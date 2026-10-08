"""Check the event-comparison and accelerated-frame examples independently.

Run with Python and SymPy. Coordinate Jacobians and Pauli matrices check
the interval, moving-frame terms, null rays, clock rates, jerk, and SI
conversions without parsing or depending on the manuscript's formulas.
"""

import sympy as s


checks = 0


def check(name, expression):
    global checks
    entries = list(expression) if isinstance(expression, s.MatrixBase) else [expression]
    residuals = [s.simplify(s.expand_trig(value)) for value in entries]
    assert all(value == 0 for value in residuals), (name, residuals)
    checks += 1
    print(f"PASS {name}")


# Exact event pairs, including both clock worldlines and a spacelike pair.
metric = s.diag(1, -1)
boost = s.Matrix([[s.Rational(5, 4), -s.Rational(3, 4)],
                  [-s.Rational(3, 4), s.Rational(5, 4)]])
for old, new, interval in [((1, 0), (s.Rational(5, 4), -s.Rational(3, 4)), 1),
                           ((s.Rational(5, 4), s.Rational(3, 4)), (1, 0), 1),
                           ((1, 1), (s.Rational(1, 2), s.Rational(1, 2)), 0),
                           ((0, 1), (-s.Rational(3, 4), s.Rational(5, 4)), -1)]:
    original, transformed = s.Matrix(old), s.Matrix(new)
    check(f"event pair {old} has the stated coordinates and interval",
          s.Matrix([*(boost * original - transformed),
                    (original.T * metric * original)[0] - interval,
                    (transformed.T * metric * transformed)[0] - interval]))

rapidity, length, mass = s.symbols("theta ell m", real=True)
gamma, beta = s.cosh(rapidity), s.tanh(rapidity)
check("transverse light clock round trip consists of null displacements",
      (gamma * length)**2 - length**2 - (beta * gamma * length)**2)
rod_pair = s.Matrix([beta * length, length])
general_boost = s.Matrix([[s.cosh(rapidity), -s.sinh(rapidity)],
                          [-s.sinh(rapidity), s.cosh(rapidity)]])
check("rod endpoints simultaneous in the measuring frame give contracted length",
      general_boost * rod_pair - s.Matrix([0, length / gamma]))
check("energy increase leaves invariant mass unchanged",
      (mass * gamma)**2 - (mass * s.sinh(rapidity))**2 - mass**2)

# Differentiate a worldline plus its moving rest direction to obtain the
# accelerated coordinate Jacobian. Its coefficients are not assumed.
clock, position = s.symbols("clock position", real=True)
theta = s.Function("theta")(clock)
origin_t, origin_x = s.Function("origin_t")(clock), s.Function("origin_x")(clock)
event_map = s.Matrix([origin_t + position * s.sinh(theta),
                     origin_x + position * s.cosh(theta)])
jacobian = event_map.jacobian([clock, position]).subs({
    origin_t.diff(clock): s.cosh(theta), origin_x.diff(clock): s.sinh(theta)})
clock_factor = 1 + position * theta.diff(clock)
check("accelerated event map yields the proper-time interval",
      jacobian.T * metric * jacobian - s.diag(clock_factor**2, -1))
inverse = s.Matrix([[s.cosh(theta) / clock_factor, -s.sinh(theta) / clock_factor],
                    [-s.sinh(theta), s.cosh(theta)]])
check("the stated increment comparison is the inverse coordinate Jacobian",
      inverse * jacobian - s.eye(2))
check("local boost measures a clock factor rather than the coordinate dt alone",
      general_boost.subs(rapidity, theta) * jacobian - s.diag(clock_factor, 1))
check("fixed grid clock rate changes with origin proper jerk",
      s.diff(clock_factor, clock) - position * theta.diff(clock, 2))
check("fixed grid proper jerk uses that grid clock's own proper time",
      s.diff(theta.diff(clock) / clock_factor, clock) / clock_factor
      - theta.diff(clock, 2) / clock_factor**3)
for sign in (-1, 1):
    null_tangent = jacobian * s.Matrix([1, sign * clock_factor])
    check(f"accelerated-coordinate ray {sign:+d} is inertially null",
          (null_tangent.T * metric * null_tangent)[0])
    check(f"accelerated-coordinate ray {sign:+d} has the inertial light direction",
          null_tangent[1] - sign * null_tangent[0])
    slope = sign * clock_factor
    along_ray = s.diff(slope, clock) + s.diff(slope, position) * slope
    check(f"coordinate light acceleration {sign:+d} includes spatial and jerk terms",
          along_ray - (sign * position * theta.diff(clock, 2)
                       + theta.diff(clock) * clock_factor))

# Constant-acceleration coordinates and their exact light rays.
acceleration = s.symbols("acceleration", positive=True)
rindler_map = s.Matrix([(1 / acceleration + position) * s.sinh(acceleration * clock),
                       (1 / acceleration + position) * s.cosh(acceleration * clock)
                       - 1 / acceleration])
rindler_jacobian = rindler_map.jacobian([clock, position])
check("Rindler map specializes the arbitrary-acceleration Jacobian",
      rindler_jacobian - jacobian.subs({theta: acceleration * clock,
                                       theta.diff(clock): acceleration}))
check("zero-acceleration limit of the full Rindler map is the inertial chart",
      rindler_map.applyfunc(lambda value: s.limit(value, acceleration, 0))
      - s.Matrix([clock, position]))
check("Rindler image satisfies the wedge hyperbola identity",
      (rindler_map[1] + 1 / acceleration)**2 - rindler_map[0]**2
      - (position + 1 / acceleration)**2)
check("accelerated origin event at rapidity log(2) has the stated inertial coordinates",
      rindler_map.subs({position: 0, clock: s.log(2) / acceleration})
      - s.Matrix([s.Rational(3, 4) / acceleration, s.Rational(1, 4) / acceleration]))
initial_position = s.symbols("initial_position", real=True)
for sign in (-1, 1):
    ray = ((1 + acceleration * initial_position)
           * s.exp(sign * acceleration * clock) - 1) / acceleration
    check(f"Rindler null solution {sign:+d} has the stated coordinate slope",
          s.diff(ray, clock) - sign * (1 + acceleration * ray))
    inertial_ray = rindler_map.subs(position, ray)
    check(f"Rindler null solution {sign:+d} maps to a straight inertial light ray",
          s.diff(inertial_ray[1] - sign * inertial_ray[0], clock))

# Initial rest expansion including proper jerk, obtained from a rapidity
# history and integration of its unit tangent.
initial_acceleration, initial_jerk = s.symbols("initial_acceleration initial_jerk", real=True)
theta_series = initial_acceleration * clock + initial_jerk * clock**2 / 2
t_series = s.integrate(s.cosh(theta_series).series(clock, 0, 3).removeO(), clock)
x_series = s.integrate(s.sinh(theta_series).series(clock, 0, 3).removeO(), clock)
check("proper jerk enters the rest-worldline position at cubic order",
      s.Matrix([t_series - clock - initial_acceleration**2 * clock**3 / 6,
                x_series - initial_acceleration * clock**2 / 2 - initial_jerk * clock**3 / 6]))

# A noncommuting complex determinant-one frame exercises acceleration and
# rotation together. The origin tangent is determined by the rest condition.
one = s.eye(2)
sigma = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))


def scalar(z):
    return s.trace(z) / 2


def vector(v):
    return sum((v[k] * sigma[k] for k in range(3)), s.zeros(2))


def spatial(z):
    return s.Matrix([scalar(item * z) for item in sigma])


frame = s.Matrix([[1, s.I * clock], [0, 1]]) * s.Matrix([[1, 0], [clock**2, 1]])
frame_inverse = frame.inv()
omega = (frame.diff(clock) * frame_inverse).expand()
rest_tangent = frame_inverse * frame_inverse.H
spatial_position = s.Matrix(s.symbols("x y z", real=True))
spatial_step = s.Matrix(s.symbols("dx dy dz", real=True))
dt = s.symbols("dt", real=True)
offset = frame_inverse * vector(spatial_position) * frame_inverse.H
inertial_step = ((rest_tangent + offset.diff(clock)) * dt
                 + frame_inverse * vector(spatial_step) * frame_inverse.H)
local_step = (frame * inertial_step * frame.H).expand()
expected_step = ((one - omega * vector(spatial_position)
                  - vector(spatial_position) * omega.H) * dt + vector(spatial_step))
check("general accelerated and rotating event map gives both frame derivative terms",
      local_step - expected_step)
real_omega = (omega + omega.H) / 2
imaginary_omega = (omega - omega.H) / (2 * s.I)
local_time = (1 - 2 * spatial(real_omega).dot(spatial_position)) * dt
local_space = spatial_step + 2 * spatial(imaginary_omega).cross(spatial_position) * dt
check("general observer displacement has the stated acceleration and rotation parts",
      local_step - local_time * one - vector(local_space))
check("general observer interval follows from the inertial-coordinate Jacobian",
      s.expand(inertial_step.det()) - local_time**2 + local_space.dot(local_space))
check("origin proper acceleration is minus twice the Hermitian frame generator",
      frame * rest_tangent.diff(clock) * frame.H + omega + omega.H)

# SI conversions retain s=c*tau and physical acceleration c*dtheta/dtau.
c, si_acceleration, si_jerk = s.symbols("c si_acceleration si_jerk", positive=True)
si_rindler = rindler_map.subs({acceleration: si_acceleration / c**2, clock: c * clock})
expected_si = s.Matrix([(c**2 / si_acceleration + position) * s.sinh(si_acceleration * clock / c),
                        (c**2 / si_acceleration + position) * s.cosh(si_acceleration * clock / c)
                        - c**2 / si_acceleration])
check("Rindler SI map is the proper-length conversion of the natural-unit map",
      si_rindler - expected_si)
si_metric = expected_si.jacobian([clock, position]).T * metric * expected_si.jacobian([clock, position])
check("SI accelerated interval uses acceleration times distance divided by c squared",
      si_metric - s.diag(c**2 * (1 + si_acceleration * position / c**2)**2, -1))
check("SI constant proper jerk rapidity follows from s=c*tau",
      (initial_jerk * clock**2 / 2).subs({initial_jerk: si_jerk / c**3, clock: c * clock})
      - si_jerk * clock**2 / (2 * c))

print(f"\n{checks} checks passed.")
