"""Independent curvature and Einstein checks for the paravector construction.

Run: python3 verify_gr_curvature.py (standard library only).
The transport generator comes from the torsion equations in verify_gr_transport.
Differentiating those equations gives its curvature without finite differences.
An independent coordinate-metric calculation uses exact rational first/second
jets. Direct 2x2 Pauli matrices check the proposed curvature and Ricci trace map.
"""

from fractions import Fraction as Q
from itertools import product
from random import Random

from verify_gr_transport import solve_transport, omega_matrix


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def scale(s, a):
    return tuple(s*x for x in a)


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


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
    return tuple(tuple(Q(i == j) for j in range(n)) for i in range(n))


def zero(n=4):
    return tuple((Q(0),)*n for _ in range(n))


def inverse(a):
    n = len(a)
    rows = [list(row)+list(unit) for row, unit in zip(a, eye(n))]
    for j in range(n):
        pivot = max(range(j, n), key=lambda i: abs(rows[i][j]))
        rows[j], rows[pivot] = rows[pivot], rows[j]
        divisor = rows[j][j]
        rows[j] = [x/divisor for x in rows[j]]
        for i in range(n):
            if i != j:
                factor = rows[i][j]
                rows[i] = [x-factor*y for x, y in zip(rows[i], rows[j])]
    return tuple(tuple(row[n:]) for row in rows)


def commutator(a, b):
    return ma(mm(a, b), ms(-1, mm(b, a)))


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def flatten(a):
    return tuple(x for row in a for x in row)


def paravector(column):
    b, x, y, z = column
    return ((b+z, x-1j*y), (x+1j*y, b-z))


def coefficients(a):
    return ((a[0][0]+a[1][1])/2, (a[0][1]+a[1][0])/2,
            (a[1][0]-a[0][1])/(2j), (a[0][0]-a[1][1])/2)


def hermitian(a):
    return tuple(tuple(complex(x).conjugate() for x in row) for row in transpose(a))


def curvature_action(curvature, carried):
    value = paravector(carried)
    return coefficients(ma(mm(curvature, value), mm(value, hermitian(curvature))))


signature = tuple(tuple(Q((1 if i == 0 else -1)*int(i == j))
                        for j in range(4)) for i in range(4))
units = eye(4)
checks = 0


def check(name, actual, expected, tolerance=2e-9):
    global checks
    error = max(abs(x-y) for x, y in zip(actual, expected))
    assert error <= tolerance*(1+max(abs(x) for x in expected)), (name, error)
    checks += 1


def check_matrix(name, actual, expected, tolerance=2e-9):
    check(name, flatten(actual), flatten(expected), tolerance)


def coordinate_curvature(frame, first, second):
    """Exact metric oracle, independent of torsion/Pauli construction.

    frame = E = M^T; first[i] = partial_i E; second[i][j] = partial_i partial_j E.
    Returns R(i,j) as coordinate matrices, the covariant Ricci matrix, and scalar.
    """
    metric = mm(mm(transpose(frame), signature), frame)
    inv = inverse(metric)
    dg = tuple(ma(mm(mm(transpose(d), signature), frame),
                  mm(mm(transpose(frame), signature), d)) for d in first)
    ddg = tuple(tuple(ma(ma(mm(mm(transpose(second[i][j]), signature), frame),
                              mm(mm(transpose(frame), signature), second[i][j])),
                           ma(mm(mm(transpose(first[i]), signature), first[j]),
                              mm(mm(transpose(first[j]), signature), first[i])))
                      for j in range(4)) for i in range(4))
    dinv = tuple(ms(-1, mm(mm(inv, d), inv)) for d in dg)

    def numerator(i, r, s):
        return dg[i][r][s]+dg[s][r][i]-dg[r][i][s]

    gamma = tuple(tuple(tuple(sum(inv[r][l]*numerator(i, l, s) for l in range(4))/2
                              for s in range(4)) for r in range(4)) for i in range(4))
    dgamma = tuple(tuple(tuple(tuple(sum(
        dinv[q][r][l]*numerator(i, l, s)
        + inv[r][l]*(ddg[q][i][l][s]+ddg[q][s][l][i]-ddg[q][l][i][s])
        for l in range(4))/2 for s in range(4)) for r in range(4))
        for i in range(4)) for q in range(4))
    curvature = tuple(tuple(ma(ma(dgamma[i][j], ms(-1, dgamma[j][i])),
                               commutator(gamma[i], gamma[j]))
                            for j in range(4)) for i in range(4))
    ricci = tuple(tuple(sum(curvature[r][j][r][i] for r in range(4))
                        for j in range(4)) for i in range(4))
    scalar = trace(mm(inv, ricci))
    return curvature, ricci, scalar


def paravector_curvature(frame, first, second):
    """Differentiate the torsion solve analytically, then use 2x2 commutators."""
    generators = solve_transport(frame, first)
    derivatives = []
    for q in range(4):
        # Differentiate K_i E_j - K_j E_i = partial_i E_j - partial_j E_i.
        # Move K_i partial_q E_j terms to the right, retaining the same solve.
        effective_first = tuple(ma(second[q][i], ms(-1, mm(generators[i], first[q])))
                                for i in range(4))
        derivatives.append(solve_transport(frame, effective_first))
    omega = tuple(omega_matrix(generator) for generator in generators)
    domega = tuple(tuple(omega_matrix(generator) for generator in derivative)
                   for derivative in derivatives)
    curvature = tuple(tuple(ma(ma(domega[j][i], ms(-1, domega[i][j])),
                               commutator(omega[i], omega[j]))
                            for j in range(4)) for i in range(4))
    return generators, curvature


def combine_curvature(curvature, q, p):
    result = zero(2)
    for i, j in product(range(4), repeat=2):
        result = ma(result, ms(q[i]*p[j], curvature[i][j]))
    return result


def ricci_trace_map(frame, curvature):
    """Literal real trace in the proposed index-free Ricci definition.

    For each physical X,Y, trace Z -> curvature(E^-1 Z,E^-1 X) acting on Y.
    Convert its resulting symmetric bilinear form into its Lorentzian map.
    """
    inv = inverse(frame)
    coordinate_units = tuple(mv(inv, unit) for unit in units)
    covariant = []
    for x in range(4):
        row = []
        for y in range(4):
            value = 0
            for z in range(4):
                generator = combine_curvature(curvature, coordinate_units[z], coordinate_units[x])
                value += curvature_action(generator, units[y])[z]
            row.append(value)
        covariant.append(tuple(row))
    return mm(signature, tuple(covariant))


def check_case(name, frame, first, second, expected_flat=False):
    coordinate, ricci, scalar = coordinate_curvature(frame, first, second)
    generators, curvature = paravector_curvature(frame, first, second)
    inv = inverse(frame)
    for i, j in product(range(4), repeat=2):
        physical = mm(mm(frame, coordinate[i][j]), inv)
        for unit in units:
            check(name+' curvature action', curvature_action(curvature[i][j], unit),
                  mv(physical, unit))
        check(name+' pure curvature generator', (trace(curvature[i][j]),), (0,))
        if expected_flat:
            check_matrix(name+' vanishing curvature', physical, zero())
            check_matrix(name+' vanishing Pauli curvature', curvature[i][j], zero(2))
    result = ricci_trace_map(frame, curvature)
    expected = mm(signature, mm(mm(transpose(inv), ricci), inv))
    check_matrix(name+' Ricci trace map', result, expected)
    check(name+' Ricci scalar', (trace(result),), (scalar,))
    check_matrix(name+' Ricci self-adjointness', mm(signature, result),
                 transpose(mm(signature, result)))
    einstein = ma(result, ms(-trace(result)/2, eye(4)))
    check(name+' Einstein trace', (trace(einstein),), (-scalar,))
    return generators, result, einstein


def run_checks():
    rng = Random(1935)

    def rational():
        return Q(rng.randrange(-3, 4), 10)

    def matrix():
        return tuple(tuple(rational() for _ in range(4)) for _ in range(4))

    # Genuine smooth arbitrary frame jets: only second derivative directions
    # are symmetrized; no symmetry is imposed on the frame matrix itself.
    for case in range(12):
        frame = ma(ms(2, eye(4)), matrix())
        first = tuple(matrix() for _ in range(4))
        seed = tuple(tuple(matrix() for _ in range(4)) for _ in range(4))
        second = tuple(tuple(ms(Q(1, 2), ma(seed[i][j], seed[j][i]))
                             for j in range(4)) for i in range(4))
        check_case('arbitrary frame', frame, first, second)

    # A local Lorentz frame with two noncommuting boosts: all metric jets
    # vanish, but connection and its derivative/commutator terms do not.
    boost_x = tuple(tuple(Q((i, j) in ((0, 1), (1, 0))) for j in range(4))
                    for i in range(4))
    boost_y = tuple(tuple(Q((i, j) in ((0, 2), (2, 0))) for j in range(4))
                    for i in range(4))
    first = (boost_x, boost_y, zero(), zero())
    second = [[zero() for _ in range(4)] for _ in range(4)]
    second[0][0], second[1][1] = mm(boost_x, boost_x), mm(boost_y, boost_y)
    second[0][1] = second[1][0] = mm(boost_x, boost_y)
    generators, _, _ = check_case('varying Lorentz frame', eye(4), first, second, True)
    assert max(abs(x) for x in flatten(commutator(generators[0], generators[1]))) > Q(1, 2)

    # E is the Jacobian of a nonlinear inertial coordinate map:
    # X0=t+x²/3, X1=x+ty/5, X2=y+z²/7, X3=z, at the origin.
    first = [[[Q(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
    first[1][0][1] = Q(2, 3)
    first[0][1][2] = first[2][1][0] = Q(1, 5)
    first[3][2][3] = Q(2, 7)
    second = tuple(tuple(zero() for _ in range(4)) for _ in range(4))
    generators, _, _ = check_case('nonlinear inertial coordinates', eye(4), first, second, True)
    for generator in generators:
        check_matrix('inertial physical transport', generator, zero())

    for factor, rate, accel in ((1, 0, 0), (2, 3, 5), (1, 1, 1), (3, 0, 2), (2, 1, 0)):
        factor, rate, accel = map(Q, (factor, rate, accel))
        frame = tuple(tuple((Q(1) if i == 0 else factor)*int(i == j)
                            for j in range(4)) for i in range(4))
        dt = tuple(tuple(rate*int(i == j and i > 0) for j in range(4)) for i in range(4))
        ddt = tuple(tuple(accel*int(i == j and i > 0) for j in range(4)) for i in range(4))
        first = (dt, zero(), zero(), zero())
        second = tuple(tuple(ddt if i == j == 0 else zero() for j in range(4)) for i in range(4))
        _, ricci, einstein = check_case('Robertson-Walker', frame, first, second)
        hubble, tidal = rate/factor, accel/factor
        expected_ricci = tuple(tuple((-3*tidal if i == 0 else -(tidal+2*hubble*hubble))
                                    *int(i == j) for j in range(4)) for i in range(4))
        check_matrix('Robertson-Walker Ricci map', ricci, expected_ricci)
        energy, pressure = 3*hubble*hubble, -(2*tidal+hubble*hubble)
        # Set 8*pi*G=1 here; restoring it gives the usual Friedmann equations.
        matter = tuple(tuple((energy if i == 0 else -pressure)*int(i == j)
                             for j in range(4)) for i in range(4))
        check_matrix('Friedmann energy and pressure signs', einstein, matter)
        assert 6*hubble*(tidal-hubble*hubble)+3*hubble*(energy+pressure) == 0

    # Weak static Newtonian field, at Phi=0, grad Phi=0. Arbitrary Hessian
    # gives an exact linear curvature test with no neglected products here.
    hessian = ((Q(2, 5), Q(-1, 7), Q(1, 9)),
               (Q(-1, 7), Q(3, 8), Q(2, 11)),
               (Q(1, 9), Q(2, 11), Q(-1, 6)))
    first = (zero(),)*4
    second = tuple(tuple(ms(hessian[i-1][j-1], signature) if i > 0 and j > 0 else zero()
                         for j in range(4)) for i in range(4))
    _, _, einstein = check_case('weak Newtonian field', eye(4), first, second)
    laplacian = trace(hessian)
    expected = tuple(tuple(2*laplacian*int(i == j == 0) for j in range(4)) for i in range(4))
    check_matrix('Poisson equation normalization', einstein, expected)

    # Generic energy density, momentum flux, symmetric spatial stress.
    # Check the proposed physical matter map's Lorentz self-adjointness,
    # trace and energy measured by arbitrary timelike observers.
    for _ in range(12):
        energy, flux = Q(3)+rational(), tuple(rational() for _ in range(3))
        seed = tuple(tuple(rational() for _ in range(3)) for _ in range(3))
        stress = ms(Q(1, 2), ma(seed, transpose(seed)))
        matter = ((energy,)+scale(-1, flux),)+tuple((flux[i],)+scale(-1, stress[i]) for i in range(3))
        check_matrix('matter self-adjointness', mm(signature, matter), transpose(mm(signature, matter)))
        check('matter trace', (trace(matter),), (energy-trace(stress),))
        velocity = tuple(rational() for _ in range(3))
        observer = (Q(1),)+velocity
        expected = energy-2*dot(flux, velocity)+dot(velocity, mv(stress, velocity))
        check('measured matter energy pairing', (dot(observer, mv(signature, mv(matter, observer))),),
              (expected,))

    # Independent Pauli multiplication for electromagnetic matter:
    # F=(E+iB).sigma, T(Y)=F Y F^H/2. Include pure fields and a null field
    # as well as generic fields, checking all four carried-vector components.
    field_cases = [((Q(1), Q(0), Q(0)), (Q(0), Q(0), Q(0))),
                   ((Q(0), Q(0), Q(0)), (Q(0), Q(1), Q(0))),
                   ((Q(1), Q(0), Q(0)), (Q(0), Q(1), Q(0)))]
    field_cases += [(tuple(rational() for _ in range(3)),
                     tuple(rational() for _ in range(3))) for _ in range(8)]
    for electric, magnetic in field_cases:
        field = paravector((0,)+tuple(e+1j*b for e, b in zip(electric, magnetic)))
        energy = (dot(electric, electric)+dot(magnetic, magnetic))/2
        flux = (electric[1]*magnetic[2]-electric[2]*magnetic[1],
                electric[2]*magnetic[0]-electric[0]*magnetic[2],
                electric[0]*magnetic[1]-electric[1]*magnetic[0])
        stress = tuple(tuple(energy*int(i == j)-electric[i]*electric[j]-magnetic[i]*magnetic[j]
                             for j in range(3)) for i in range(3))
        matter = ((energy,)+scale(-1, flux),)+tuple((flux[i],)+scale(-1, stress[i]) for i in range(3))
        for unit in units:
            direct = coefficients(ms(Q(1, 2), mm(mm(field, paravector(unit)), hermitian(field))))
            check('electromagnetic paravector stress map', direct, mv(matter, unit))
        check('electromagnetic matter trace', (trace(matter),), (0,))

    print(f'PASS: {checks} independent curvature, Ricci trace, and Einstein checks; '
          '12 arbitrary frames, two nonconstant flat controls, five exact cosmology cases, '
          'Newtonian Poisson normalization, scalar/vector matter blocks, '
          'and electromagnetic energy, flux, and stress.')


if __name__ == '__main__':
    run_checks()
