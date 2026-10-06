"""Check the index-free paravector transport construction in Section 19.

Run: python3 verify_gr_transport.py

E = M^T maps coordinate tangent columns into physical scalar/vector columns.
dE[i] is its derivative in coordinate direction i.  solve_transport returns
the four parallel-transport generators K_i, so y' = K(q)y.  These have the
opposite sign to the physical covariant-derivative connection.

The construction under test uses only the algebraic torsion equation.  A
separately implemented coordinate metric connection provides its oracle.
This module has no third-party dependencies and can be imported safely.
"""

from math import cos, cosh, exp, sin, sinh
from random import Random


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale(s, a):
    return tuple(s * x for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0])


def transpose(a):
    return tuple(zip(*a))


def mv(a, x):
    return tuple(dot(row, x) for row in a)


def mm(a, b):
    return tuple(tuple(dot(row, col) for col in transpose(b)) for row in a)


def ma(a, b):
    return tuple(add(row, other) for row, other in zip(a, b))


def ms(s, a):
    return tuple(scale(s, row) for row in a)


def eye(n):
    return tuple(tuple(int(i == j) for j in range(n)) for i in range(n))


def zeros(n=4):
    return tuple((0,) * n for _ in range(n))


def solve_linear(a, b):
    """Partial-pivot Gaussian elimination; works with real or complex values."""
    n = len(a)
    rows = [list(row) + [value] for row, value in zip(a, b)]
    for j in range(n):
        pivot = max(range(j, n), key=lambda i: abs(rows[i][j]))
        rows[j], rows[pivot] = rows[pivot], rows[j]
        divisor = rows[j][j]
        if abs(divisor) < 1e-14:
            raise ValueError('Singular or poorly conditioned test matrix')
        rows[j][j:] = [x / divisor for x in rows[j][j:]]
        for i in range(j + 1, n):
            factor = rows[i][j]
            rows[i][j:] = [x - factor*y
                           for x, y in zip(rows[i][j:], rows[j][j:])]
    result = [0] * n
    for i in range(n - 1, -1, -1):
        result[i] = rows[i][-1] - dot(rows[i][i+1:n], result[i+1:n])
    return tuple(result)


def inverse(a):
    return transpose(tuple(solve_linear(a, column) for column in eye(len(a))))


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def flatten(a):
    return tuple(x for row in a for x in row)


def adjoint(a):
    return tuple(tuple(complex(x).conjugate() for x in row)
                 for row in transpose(a))


def paravector(c, z):
    """Independent 2 by 2 Pauli representation."""
    return ((c + z[2], z[0] - 1j*z[1]),
            (z[0] + 1j*z[1], c - z[2]))


ETA = ((1, 0, 0, 0), (0, -1, 0, 0),
       (0, 0, -1, 0), (0, 0, 0, -1))
BASIS = eye(4)


def generator(alpha, beta):
    """Action of Omega = (alpha + i beta).sigma / 2 on a real paravector."""
    return ((0,) + tuple(alpha),
            (alpha[0], 0, beta[2], -beta[1]),
            (alpha[1], -beta[2], 0, beta[0]),
            (alpha[2], beta[1], -beta[0], 0))


def transport_parameters(k):
    return ((k[1][0], k[2][0], k[3][0]),
            (k[2][3], k[3][1], k[1][2]))


def omega_matrix(k):
    alpha, beta = transport_parameters(k)
    return paravector(0, tuple((a + 1j*b)/2 for a, b in zip(alpha, beta)))


def solve_transport(e, de):
    """Solve K_i E_j - K_j E_i = partial_i E_j - partial_j E_i.

    Each K_i is constrained to a six-parameter Lorentz generator, which
    builds metric compatibility into the linear system.  No metric
    derivatives or coordinate connection enter this solver.
    """
    pairs = tuple((i, j) for i in range(4) for j in range(i + 1, 4))
    columns = transpose(e)
    rhs = tuple(value for i, j in pairs for value in
                add(mv(de[i], BASIS[j]), scale(-1, mv(de[j], BASIS[i]))))
    unit_generators = []
    for n in range(6):
        unit = tuple(int(n == k) for k in range(6))
        unit_generators.append(generator(unit[:3], unit[3:]))
    system_columns = []
    for direction in range(4):
        for unit in unit_generators:
            column = []
            for i, j in pairs:
                value = (0, 0, 0, 0)
                if direction == i:
                    value = add(value, mv(unit, columns[j]))
                if direction == j:
                    value = add(value, scale(-1, mv(unit, columns[i])))
                column.extend(value)
            system_columns.append(tuple(column))
    values = solve_linear(transpose(system_columns), rhs)
    return tuple(generator(values[6*i:6*i+3], values[6*i+3:6*i+6])
                 for i in range(4))


def metric_connection(e, de):
    """Independent positive coordinate Levi-Civita connection Gamma_i."""
    g = mm(transpose(e), mm(ETA, e))
    gi = inverse(g)
    dg = tuple(ma(mm(transpose(d), mm(ETA, e)),
                  mm(transpose(e), mm(ETA, d))) for d in de)
    result = []
    for i in range(4):
        columns = []
        for j in range(4):
            covector = tuple((dg[i][k][j] + dg[j][k][i] - dg[k][i][j])/2
                             for k in range(4))
            columns.append(mv(gi, covector))
        result.append(transpose(columns))
    return tuple(result)


def assert_matrix_close(actual, expected, tolerance=2e-11):
    error = max(abs(a - b) for a, b in zip(flatten(actual), flatten(expected)))
    bound = tolerance * (1 + max(abs(b) for b in flatten(expected)))
    assert error <= bound, (error, bound)


def lorentz_frame(event):
    """Noncommuting local boost and rotation: variable frame, flat metric."""
    t, x, _, _ = event
    c, s = cosh(t), sinh(t)
    boost = ((c, s, 0, 0), (s, c, 0, 0),
             (0, 0, 1, 0), (0, 0, 0, 1))
    dboost = ((s, c, 0, 0), (c, s, 0, 0),
              (0, 0, 0, 0), (0, 0, 0, 0))
    c, s = cos(x), sin(x)
    rotation = ((1, 0, 0, 0), (0, c, -s, 0),
                (0, s, c, 0), (0, 0, 0, 1))
    drotation = ((0, 0, 0, 0), (0, -s, -c, 0),
                 (0, c, -s, 0), (0, 0, 0, 0))
    return (mm(boost, rotation),
            (mm(dboost, rotation), mm(boost, drotation), zeros(), zeros()))


def coordinate_frame(event):
    """Jacobian of the nonlinear inertial coordinate map (t, exp(t)x, y, z)."""
    t, x, _, _ = event
    a = exp(t)
    e = ((1, 0, 0, 0), (a*x, a, 0, 0),
         (0, 0, 1, 0), (0, 0, 0, 1))
    dt = ((0, 0, 0, 0), (a*x, a, 0, 0),
          (0, 0, 0, 0), (0, 0, 0, 0))
    dx = ((0, 0, 0, 0), (a, 0, 0, 0),
          (0, 0, 0, 0), (0, 0, 0, 0))
    return e, (dt, dx, zeros(), zeros())


def run_checks():
    rng = Random(381)
    checks = 0
    for trial in range(24):
        e = tuple(tuple(int(i == j) + rng.uniform(-.3, .3)
                        for j in range(4)) for i in range(4))
        if trial % 3 == 1:
            # Invertibility, rather than a chosen orientation, is the
            # algebraic hypothesis of the construction.
            e = (e[1], e[0], e[2], e[3])
        elif trial % 3 == 2:
            e = (scale(-1, e[0]),) + e[1:]
        de = tuple(tuple(tuple(rng.uniform(-1, 1) for j in range(4))
                         for i in range(4)) for _ in range(4))
        ei = inverse(e)
        connection = solve_transport(e, de)
        coordinate = metric_connection(e, de)
        for i, k in enumerate(connection):
            expected = ma(mm(de[i], ei), ms(-1, mm(e, mm(coordinate[i], ei))))
            assert_matrix_close(k, expected)
            assert_matrix_close(ma(mm(ETA, k), mm(transpose(k), ETA)), zeros())
            # A two-by-two product oracle verifies the cross-product sign.
            y = tuple(rng.uniform(-1, 1) for _ in range(4))
            omega = omega_matrix(k)
            py = paravector(y[0], y[1:])
            action = ma(mm(omega, py), mm(py, adjoint(omega)))
            ky = mv(k, y)
            assert_matrix_close(action, paravector(ky[0], ky[1:]))
            assert abs(dot(y, mv(ETA, ky))) < 1e-11
            checks += 4
        for i in range(4):
            for j in range(i + 1, 4):
                lhs = add(mv(connection[i], mv(e, BASIS[j])),
                          scale(-1, mv(connection[j], mv(e, BASIS[i]))))
                rhs = add(mv(de[i], BASIS[j]), scale(-1, mv(de[j], BASIS[i])))
                assert max(abs(a - b) for a, b in zip(lhs, rhs)) < 1e-11
                checks += 1

    for event in ((.3, -.4, .5, .8), (-.6, .8, 0, 0), (.8, .2, 0, 0)):
        e, de = lorentz_frame(event)
        ei = inverse(e)
        connection = solve_transport(e, de)
        assert_matrix_close(mm(transpose(e), mm(ETA, e)), ETA)
        assert max(abs(v) for k in connection for v in flatten(k)) > .5
        for k, derivative, coordinate in zip(connection, de, metric_connection(e, de)):
            assert_matrix_close(k, mm(derivative, ei))
            assert_matrix_close(coordinate, zeros())
            checks += 2
        # For our positive transport convention, covariant curvature is
        # partial_j K_i - partial_i K_j + [K_i,K_j].  This control exercises
        # nonzero derivatives AND a nonzero commutator; the wrong sign fails.
        step = 2e-5
        plus = list(event); plus[0] += step
        minus = list(event); minus[0] -= step
        kp = solve_transport(*lorentz_frame(plus))[1]
        km = solve_transport(*lorentz_frame(minus))[1]
        derivative = ms(1/(2*step), ma(kp, ms(-1, km)))
        commutator = ma(mm(connection[0], connection[1]),
                        ms(-1, mm(connection[1], connection[0])))
        assert max(abs(v) for v in flatten(commutator)) > .5
        assert_matrix_close(ma(ms(-1, derivative), commutator), zeros(), 5e-9)
        checks += 2

        e, de = coordinate_frame(event)
        ei = inverse(e)
        for k, derivative, coordinate in zip(solve_transport(e, de), de,
                                             metric_connection(e, de)):
            assert_matrix_close(k, zeros())
            assert_matrix_close(coordinate, mm(ei, derivative))
            checks += 2

    # Recover the already-derived cosmological transport law for arbitrary
    # tangent and carried vectors, including contraction and stationary cases.
    for a, adot in ((1, 0), (2, 3), (1.5, -.4), (.7, 2)):
        e = ((1, 0, 0, 0), (0, a, 0, 0),
             (0, 0, a, 0), (0, 0, 0, a))
        dt = ((0, 0, 0, 0), (0, adot, 0, 0),
              (0, 0, adot, 0), (0, 0, 0, adot))
        connection = solve_transport(e, (dt, zeros(), zeros(), zeros()))
        for i, k in enumerate(connection):
            expected = generator(scale(-adot, BASIS[i][1:]), (0, 0, 0))
            assert_matrix_close(k, expected)
            checks += 1
    print(f'PASS: {checks} paravector transport checks; 24 arbitrary frames, '
          'nonconstant flat Lorentz and nonlinear-coordinate controls, '
          'Pauli action and curvature signs, and cosmological reduction.')


if __name__ == '__main__':
    run_checks()
