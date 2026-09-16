"""Exact checks for the paravector treatment of local frames and curvature.

Run with ``python3 verify_curvature.py`` (requires SymPy).  Ordinary Pauli
matrices provide an independent representation for checking the manuscript's
paravector products.  The calculations use symbolic directions and fields,
including noncommuting frame changes; no numerical tolerances are used.
The primary frame law is Z' = T.H * Z * T, with right generator
Omega = T.inverse() * dT.  Curvature therefore has a plus commutator.
"""

from time import perf_counter

import sympy as s


started = perf_counter()
checks = 0
one = s.eye(2)
zero = s.zeros(2)
sigma = (
    s.Matrix([[0, 1], [1, 0]]),
    s.Matrix([[0, -s.I], [s.I, 0]]),
    s.diag(1, -1),
)


def vec(values):
    return sum((value * basis for value, basis in zip(values, sigma)), zero)


def scalar(value):
    return s.trace(value) / 2


def adj(value):
    return s.trace(value) * one - value


def pair(first, second):
    return scalar(adj(first) * second)


def comm(first, second):
    return first * second - second * first


def check(name, expression):
    global checks
    entries = expression if isinstance(expression, s.MatrixBase) else (expression,)
    for entry in entries:
        assert s.cancel(s.expand(entry)) == 0, name
    checks += 1
    print(f"PASS {name}", flush=True)


# The existing one-parameter identities, with noncommuting matrix factors.
parameter = s.symbols("parameter", real=True)
frame = s.Matrix([[1 + s.I * parameter**2, parameter], [s.I * parameter, 1]])
inverse = frame.inv()
connection = inverse * frame.diff(parameter)
field = s.Matrix(
    [[1 + parameter**2, parameter + s.I * parameter**2],
     [parameter - s.I * parameter**2, 2 + 3 * parameter]]
)
similarity = frame.H * field * inverse.H
congruence = frame.H * field * frame
check("one-parameter inverse derivative", inverse.diff(parameter) + connection * inverse)
check(
    "one-parameter similarity derivative",
    similarity.diff(parameter) - frame.H * field.diff(parameter) * inverse.H
    - comm(connection.H, similarity),
)
check(
    "one-parameter second similarity derivative",
    similarity.diff(parameter, 2) - frame.H * field.diff(parameter, 2) * inverse.H
    - 2 * comm(connection.H, similarity.diff(parameter))
    - comm(connection.diff(parameter).H, similarity)
    + comm(connection.H, comm(connection.H, similarity)),
)
check(
    "one-parameter congruence derivative",
    congruence.diff(parameter) - frame.H * field.diff(parameter) * frame
    - connection.H * congruence - congruence * connection,
)
check(
    "one-parameter algebra-basis derivatives",
    s.Matrix.vstack(*(
        (frame.H * basis * inverse.H).diff(parameter)
        - comm(connection.H, frame.H * basis * inverse.H)
        for basis in sigma
    )),
)


# Constant real paravector directions A and B have four arbitrary coefficients.
t, x, y, z = coordinates = s.symbols("t x y z", real=True)
a = s.symbols("a0:4", real=True)
b = s.symbols("b0:4", real=True)


def derivative(field, direction):
    initial = zero if isinstance(field, s.MatrixBase) else s.S.Zero
    return sum(
        (coefficient * s.diff(field, coordinate)
         for coefficient, coordinate in zip(direction, coordinates)),
        initial,
    )


def along(values, direction):
    return sum((coefficient * value for coefficient, value in zip(direction, values)), zero)


def covariant(field, direction, omega):
    value = along(omega, direction)
    return derivative(field, direction) - value.H * field - field * value


def field_covariant(field, direction, omega):
    return derivative(field, direction) - comm(along(omega, direction).H, field)


def curvature(first, second, omega):
    first_value = along(omega, first)
    second_value = along(omega, second)
    return (
        derivative(second_value, first) - derivative(first_value, second)
        + comm(first_value, second_value)
    )


# Products of two independent nilpotent shears give a genuinely noncommuting
# determinant-one frame family.  Its inverse remains polynomial.
frame = s.Matrix([[1 + s.I * t * x, t], [s.I * x, 1]])
inverse = frame.inv()
pure_frame = tuple(inverse * frame.diff(coordinate) for coordinate in coordinates)
check("two-coordinate frame determinant", frame.det() - 1)
check("pure-frame connection scalar part", s.Matrix([scalar(value) for value in pure_frame]))
check("pure-frame curvature in arbitrary directions", curvature(a, b, pure_frame))

field = (1 + t * x) * one + vec((t + y, x + z, t * z))
second_field = (2 + y * z) * one + vec((x * y, t - z, 1 + t * y))
check(
    "pure-frame derivative commutator",
    covariant(covariant(field, b, pure_frame), a, pure_frame)
    - covariant(covariant(field, a, pure_frame), b, pure_frame),
)

# These independent pure-vector connections have nonzero curvature and both
# real boost and imaginary rotation components.
omega = (
    t * sigma[0] + s.I * x * sigma[1],
    x * sigma[1] + s.I * y * sigma[2],
    y * sigma[2] + s.I * z * sigma[0],
    z * sigma[0] + s.I * t * sigma[1],
)
value = curvature(a, b, omega)
check(
    "curved derivative commutator and sign",
    covariant(covariant(field, b, omega), a, omega)
    - covariant(covariant(field, a, omega), b, omega)
    + value.H * field + field * value,
)
complex_field = vec((t + s.I * x, y - s.I * z, t * y + s.I * z**2))
check(
    "similarity-field curvature action",
    field_covariant(field_covariant(complex_field, b, omega), a, omega)
    - field_covariant(field_covariant(complex_field, a, omega), b, omega)
    + comm(value.H, complex_field),
)
check(
    "metric-pair derivative compatibility",
    derivative(pair(field, second_field), a)
    - pair(covariant(field, a, omega), second_field)
    - pair(field, covariant(second_field, a, omega)),
)

# This independent polynomial identity covers arbitrary real paravectors and
# arbitrary complex pure-vector connections at a point.
first_parts = s.symbols("p0:4", real=True)
second_parts = s.symbols("q0:4", real=True)
boost_parts = s.symbols("boost0:3", real=True)
rotation_parts = s.symbols("rotation0:3", real=True)
first = first_parts[0] * one + vec(first_parts[1:])
second = second_parts[0] * one + vec(second_parts[1:])
generator = vec(tuple(real + s.I * imag for real, imag in zip(boost_parts, rotation_parts)))
check(
    "general infinitesimal interval compatibility",
    pair(generator.H * first + first * generator, second)
    + pair(first, generator.H * second + second * generator),
)

transformed_omega = tuple(
    inverse * value * frame + inverse * frame.diff(coordinate)
    for coordinate, value in zip(coordinates, omega)
)
transformed_field = frame.H * field * frame
check(
    "connection transformation preserves pure-vector form",
    s.Matrix([scalar(value) for value in transformed_omega]),
)
check(
    "local-frame covariant derivative transformation",
    covariant(transformed_field, a, transformed_omega)
    - frame.H * covariant(field, a, omega) * frame,
)
check(
    "local-frame similarity-field derivative transformation",
    field_covariant(frame.H * complex_field * inverse.H, a, transformed_omega)
    - frame.H * field_covariant(complex_field, a, omega) * inverse.H,
)
check(
    "local-frame curvature transformation",
    curvature(a, b, transformed_omega) - inverse * curvature(a, b, omega) * frame,
)
check(
    "local-frame interval-pair invariance",
    pair(transformed_field, frame.H * second_field * frame) - pair(field, second_field),
)


# Static lapse example.  Use the value and first two derivatives of an arbitrary
# positive real function f at a point.  This avoids assuming a polynomial f.
# The formulas are checked for arbitrary real u, hence also for every unit u.
f = s.symbols("f", positive=True, real=True)
df, ddf = s.symbols("df ddf", real=True)
u = s.Matrix(s.symbols("u0:3", real=True))
unit_axis = s.Matrix([s.Rational(2, 3), -s.Rational(2, 3), s.Rational(1, 3)])
check("example unit axis with three nonzero components", unit_axis.dot(unit_axis) - 1)


def spatial_dot(direction):
    return u.dot(s.Matrix(direction[1:]))


def frame_map(direction):
    return f * direction[0] * one + vec(direction[1:])


def example_omega(direction):
    return -direction[0] * df * vec(u) / 2


def example_derivative(field, direction):
    return spatial_dot(direction) * (field.diff(f) * df + field.diff(df) * ddf)


def example_covariant(field, direction):
    connection = example_omega(direction)
    return example_derivative(field, direction) - connection.H * field - field * connection


check(
    "static lapse interval",
    frame_map(a).det() - f**2 * a[0]**2 + sum(value**2 for value in a[1:]),
)
check(
    "static lapse connection is pure vector",
    scalar(example_omega(a)),
)
check(
    "static lapse torsion-free condition in arbitrary directions",
    example_covariant(frame_map(b), a) - example_covariant(frame_map(a), b),
)
example_curvature = (
    example_derivative(example_omega(b), a)
    - example_derivative(example_omega(a), b)
    + comm(example_omega(a), example_omega(b))
)
expected_curvature = ddf * (a[0] * spatial_dot(b) - b[0] * spatial_dot(a)) * vec(u) / 2
check("static lapse curvature formula", example_curvature - expected_curvature)
check("linear lapse has zero curvature", example_curvature.subs(ddf, 0))

time_direction = (1, 0, 0, 0)
axis_direction = (0, *unit_axis)
time_axis_curvature = expected_curvature.subs(
    dict(zip(a + b + tuple(u), time_direction + axis_direction + tuple(unit_axis)))
)
check("nonlinear lapse curvature in time-axis plane", time_axis_curvature - ddf * vec(unit_axis) / 2)
check(
    "nonlinear lapse curvature has nontrivial action on real tangents",
    -time_axis_curvature - time_axis_curvature.H + ddf * vec(unit_axis),
)

time_speed, time_acceleration = s.symbols("time_speed time_acceleration", real=True)
velocity = s.Matrix(s.symbols("velocity0:3", real=True))
acceleration = s.Matrix(s.symbols("acceleration0:3", real=True))
tangent = f * time_speed * one + vec(velocity)
tangent_derivative = (
    (df * u.dot(velocity) * time_speed + f * time_acceleration) * one
    + vec(acceleration)
)
connection_on_path = -time_speed * df * vec(u) / 2
transport = tangent_derivative - connection_on_path.H * tangent - tangent * connection_on_path
expected_transport = (
    (f * time_acceleration + 2 * df * time_speed * u.dot(velocity)) * one
    + vec(acceleration + f * df * time_speed**2 * u)
)
check("static lapse geodesic scalar and vector equations", transport - expected_transport)
determinant_derivative = (
    2 * f * df * u.dot(velocity) * time_speed**2
    + 2 * f**2 * time_speed * time_acceleration
    - 2 * velocity.dot(acceleration)
)
geodesic_substitutions = {
    time_acceleration: -2 * df * time_speed * u.dot(velocity) / f,
    **dict(zip(acceleration, -f * df * time_speed**2 * u)),
}
check(
    "geodesic preserves the tangent determinant",
    determinant_derivative.subs(geodesic_substitutions, simultaneous=True),
)

print(f"\n{checks} exact checks passed in {perf_counter() - started:.2f} seconds.")
