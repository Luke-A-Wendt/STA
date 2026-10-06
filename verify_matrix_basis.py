"""Independent checks for the matrix bases in the closing variable-basis section.

Run: python3 verify_matrix_basis.py (standard library only).
Pauli matrices and finite differences check algebra, all three conjugations,
and derivatives for real and complex nonorthogonal frames. Coordinate metric calculations
check the wave and geodesic formulas. Exact rational calculations check the
Robertson--Walker tidal terms and squared Riemann curvature.
Kerr--Schild checks cover the basis, rotating metric, Schwarzschild limit,
and an independent finite-difference vacuum Ricci calculation.
"""

from fractions import Fraction as Q
from itertools import permutations, product
from math import sqrt
from random import Random


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def scale(s, a):
    return tuple(s*x for x in a)


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


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


def inverse(a):
    n = len(a)
    rows = [list(row) + list(unit) for row, unit in zip(a, eye(n))]
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


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def det3(a):
    return dot(a[0], cross(a[1], a[2]))


def paravector(c, z):
    return ((c+z[2], z[0]-1j*z[1]), (z[0]+1j*z[1], c-z[2]))


def flatten(a):
    return tuple(x for row in a for x in row)


checks = 0


def check(name, actual, expected, tolerance=3e-10):
    global checks
    error = max(abs(x-y) for x, y in zip(actual, expected))
    assert error <= tolerance*(1+max(abs(x) for x in expected)), (name, error)
    checks += 1


def check_matrix(name, actual, expected, tolerance=3e-10):
    check(name, flatten(actual), flatten(expected), tolerance)


def connection(metric, derivatives):
    inv = inverse(metric)
    return tuple(tuple(tuple(sum(inv[r][l]*(derivatives[m][l][n]
        + derivatives[n][l][m] - derivatives[l][m][n]) for l in range(4))/2
        for n in range(4)) for m in range(4)) for r in range(4))


rng = Random(1141)


def real():
    return rng.uniform(-0.6, 0.6)


def scalar():
    return complex(real(), real())


def vector():
    return tuple(scalar() for _ in range(3))


for case in range(48):
    frame = tuple(tuple(2*int(i == j)+real() for j in range(3)) for i in range(3))
    if case % 2:
        frame = (scale(-1, frame[0]),) + frame[1:]
    if case % 8 == 0:
        frame = eye(3)
    rate = tuple(tuple(tuple(real() for _ in range(3)) for _ in range(3)) for _ in range(4))
    if case % 8 == 1:
        rate = (((0,)*3,)*3,)*4
    inv, ft = inverse(frame), transpose(frame)
    h = mm(frame, ft)
    hi = inverse(h)
    dh = tuple(ma(mm(d, ft), mm(frame, transpose(d))) for d in rate)
    c, b, x, y, z = scalar(), scalar(), vector(), vector(), vector()
    dc = tuple(scalar() for _ in range(4))
    dx, dy, dz = (tuple(vector() for _ in range(4)) for _ in range(3))
    if case % 8 == 2:
        dc, dx, dy, dz = (0,)*4, ((0,)*3,)*4, ((0,)*3,)*4, ((0,)*3,)*4

    def represent(s, v):
        return paravector(s, mv(ft, v))

    p, q = represent(c, x), represent(b, y)
    cv = cross(x, y)
    product_scalar = c*b + dot(x, mv(h, y))
    product_vector = add(add(scale(c, y), scale(b, x)), scale(1j*det3(frame), mv(hi, cv)))
    check_matrix('full paravector product', mm(p, q), represent(product_scalar, product_vector))
    determinant = c*c-dot(x, mv(h, x))
    check('determinant', (p[0][0]*p[1][1]-p[0][1]*p[1][0],), (determinant,))
    check_matrix('inverse', inverse(p), ms(1/determinant, represent(c, scale(-1, x))))
    check('cross transformation', cross(mv(ft, x), mv(ft, y)), scale(det3(frame), mv(inv, cv)))
    check('reciprocal scalar pairing', (dot(mv(ft, x), mv(inv, y)),), (dot(x, y),))

    field_rates = []
    for n in range(4):
        dt = transpose(rate[n])
        body = mm(rate[n], inv)
        original_rate = add(mv(dt, z), mv(ft, dz[n]))
        actual = paravector(dc[n], original_rate)
        predicted_vector = add(dz[n], mv(transpose(body), z))
        check_matrix('field derivative', actual, represent(dc[n], predicted_vector))
        field_rates.append(actual)

        left, right = represent(0, x), represent(0, y)
        dleft = paravector(0, add(mv(dt, x), mv(ft, dx[n])))
        dright = paravector(0, add(mv(dt, y), mv(ft, dy[n])))
        actual_product = ma(mm(dleft, right), mm(left, dright))
        predicted_scalar = dot(dx[n], mv(h, y)) + dot(x, mv(h, dy[n])) + dot(x, mv(dh[n], y))
        cross_rate = add(add(cross(dx[n], y), cross(x, dy[n])),
                         add(scale(trace(body), cv), scale(-1, mv(body, cv))))
        predicted_product = represent(predicted_scalar, scale(1j*det3(frame), mv(hi, cross_rate)))
        check_matrix('whole product derivative', actual_product, predicted_product)

        def sample(epsilon):
            shifted = transpose(ma(frame, ms(epsilon, rate[n])))
            field = paravector(c+epsilon*dc[n], mv(shifted, add(z, scale(epsilon, dz[n]))))
            xp = paravector(0, mv(shifted, add(x, scale(epsilon, dx[n]))))
            yp = paravector(0, mv(shifted, add(y, scale(epsilon, dy[n]))))
            return field, mm(xp, yp)

        epsilon = 1e-5
        plus, minus = sample(epsilon), sample(-epsilon)
        for j, expected in enumerate((actual, actual_product)):
            central = ms(1/(2*epsilon), ma(plus[j], ms(-1, minus[j])))
            check_matrix('finite difference', central, expected, tolerance=2e-8)

    # Differentiate operands first; then multiply by reciprocal sigmas.
    reciprocal = transpose(inv)
    direct = ((0, 0), (0, 0))
    for n in range(3):
        direct = ma(direct, mm(paravector(0, reciprocal[n]), field_rates[n+1]))
    z_rates = tuple(add(mv(transpose(rate[n+1]), z), mv(ft, dz[n+1])) for n in range(3))
    result_scalar = sum(inv[a][n]*z_rates[n][a] for a in range(3) for n in range(3))
    result_vector = mv(inv, dc[1:])
    for a in range(3):
        for n in range(3):
            result_vector = add(result_vector, scale(1j*inv[a][n], cross(eye(3)[a], z_rates[n])))
    check_matrix('reciprocal contraction', direct, represent(result_scalar, mv(transpose(inv), result_vector)))

    # Independent four-dimensional coordinate-metric calculation.
    metric = ((1, 0, 0, 0),) + tuple((0,)+scale(-1, row) for row in h)
    metric_rates = tuple(((0,)*4,)+tuple((0,)+scale(-1, row) for row in derivative) for derivative in dh)
    metric_inv = inverse(metric)
    gamma = connection(metric, metric_rates)
    tangent = tuple(real() for _ in range(4))
    direct_acceleration = tuple(-sum(gamma[r][m][n]*tangent[m]*tangent[n]
                                     for m in range(4) for n in range(4)) for r in range(4))
    spatial = tangent[1:]
    time_acceleration = -dot(spatial, mv(dh[0], spatial))/2
    total_h = ms(tangent[0], dh[0])
    for n in range(3):
        total_h = ma(total_h, ms(spatial[n], dh[n+1]))
    gradient = tuple(dot(spatial, mv(dh[n+1], spatial))/2 for n in range(3))
    spatial_acceleration = mv(hi, add(gradient, scale(-1, mv(total_h, spatial))))
    check('geodesic Euler-Lagrange equations', direct_acceleration,
          (time_acceleration,)+spatial_acceleration)
    norm_rate = sum(tangent[k]*dot(tangent, mv(metric_rates[k], tangent)) for k in range(4))
    norm_rate += 2*dot(tangent, mv(metric, direct_acceleration))
    check('preserved tangent norm', (norm_rate,), (0,))

    # Polarized law: transport an arbitrary vector, not just the tangent.
    components = tuple(real() for _ in range(4))
    transport_direct = tuple(-sum(gamma[r][m][n]*tangent[m]*components[n]
                                 for m in range(4) for n in range(4)) for r in range(4))
    transport_scalar = -dot(spatial, mv(dh[0], components[1:]))/2
    transport_vector = add(scale(tangent[0], mv(dh[0], components[1:])),
                           scale(components[0], mv(dh[0], spatial)))
    for n in range(3):
        transport_vector = add(transport_vector,
            add(scale(spatial[n], mv(dh[n+1], components[1:])),
                scale(components[n+1], mv(dh[n+1], spatial))))
    gradient_pair = tuple(dot(spatial, mv(dh[n+1], components[1:])) for n in range(3))
    transport_vector = scale(-0.5, mv(hi, add(transport_vector, scale(-1, gradient_pair))))
    check('metric parallel transport', transport_direct, (transport_scalar,)+transport_vector)
    transported_norm_rate = sum(tangent[k]*dot(components, mv(metric_rates[k], components)) for k in range(4))
    transported_norm_rate += 2*dot(components, mv(metric, transport_direct))
    check('transport preserves metric pairing', (transported_norm_rate,), (0,))

    f_rate = tuple(real() for _ in range(4))
    raw_hessian = tuple(tuple(real() for _ in range(4)) for _ in range(4))
    hessian = ma(raw_hessian, transpose(raw_hessian))
    wave_direct = sum(metric_inv[m][n]*(hessian[m][n]
        -sum(gamma[r][m][n]*f_rate[r] for r in range(4))) for m in range(4) for n in range(4))
    density_log_rate = tuple(trace(mm(inv, d)) for d in rate)
    wave_predicted = hessian[0][0]+density_log_rate[0]*f_rate[0]
    for n in range(3):
        inverse_rate = ms(-1, mm(mm(hi, dh[n+1]), hi))
        for m in range(3):
            wave_predicted -= hi[n][m]*hessian[n+1][m+1]
            wave_predicted -= (density_log_rate[n+1]*hi[n][m]+inverse_rate[n][m])*f_rate[m+1]
    check('metric scalar wave operator', (wave_direct,), (wave_predicted,))


# Full four-dimensional basis: temporal and spatial elements can all mix.
signature = tuple(tuple((1 if i == 0 else -1)*int(i == j)
                         for j in range(4)) for i in range(4))


def coefficients(matrix):
    return ((matrix[0][0]+matrix[1][1])/2,
            (matrix[0][1]+matrix[1][0])/2,
            (matrix[1][0]-matrix[0][1])/(2j),
            (matrix[0][0]-matrix[1][1])/2)


for case in range(64):
    entry = real if case < 32 else scalar
    frame = tuple(tuple(2*int(i == j)+entry() for j in range(4)) for i in range(4))
    if case % 2:
        frame = (scale(-1, frame[0]),) + frame[1:]
    ft, inv = transpose(frame), inverse(frame)
    inv_ft = transpose(inv)
    metric = mm(mm(frame, signature), ft)
    q, p = (tuple(scalar() for _ in range(4)) for _ in range(2))

    def full_represent(column):
        original = mv(ft, column)
        return paravector(original[0], original[1:])

    value = full_represent(q)
    other = full_represent(p)
    original_q, original_p = mv(ft, q), mv(ft, p)
    determinant = value[0][0]*value[1][1]-value[0][1]*value[1][0]
    check('full-basis determinant metric', (determinant,), (dot(q, mv(metric, q)),))
    check('coefficient extraction', mv(inv_ft, coefficients(value)), q)
    check_matrix('true algebra identity', full_represent(mv(inv_ft, (1, 0, 0, 0))), eye(2))
    adjugate_q = mv(inv_ft, mv(signature, original_q))
    adjugate_direct = ((value[1][1], -value[0][1]), (-value[1][0], value[0][0]))
    check_matrix('full-basis adjugation', full_represent(adjugate_q), adjugate_direct)
    check_matrix('full-basis inverse', full_represent(scale(1/determinant, adjugate_q)), inverse(value))
    sigma_q = mv(inv_ft, mv(signature, tuple(z.conjugate() for z in original_q)))
    s_matrix = ((0, 1), (-1, 0))
    conjugated_matrix = tuple(tuple(z.conjugate() for z in row) for row in value)
    sigma_direct = mm(mm(ms(-1, s_matrix), conjugated_matrix), s_matrix)
    check_matrix('full-basis sigma conjugation', full_represent(sigma_q), sigma_direct)
    hermitian_direct = transpose(conjugated_matrix)
    hermitian_q = mv(inv_ft, tuple(z.conjugate() for z in original_q))
    check_matrix('full-basis Hermitian conjugation', full_represent(hermitian_q), hermitian_direct)
    check_matrix('sigma then Hermitian is adjugation',
                 transpose(tuple(tuple(z.conjugate() for z in row) for row in sigma_direct)),
                 adjugate_direct)
    # Compare the three basis-column laws to direct 2x2 matrix operations.
    # Sigma and Hermitian conjugation also conjugate the field coefficients;
    # adjugation is complex-linear.
    sigma_from_basis = ((0j, 0j), (0j, 0j))
    hermitian_from_basis = ((0j, 0j), (0j, 0j))
    adjugate_from_basis = ((0j, 0j), (0j, 0j))
    for coefficient, unit, row in zip(q, eye(4), frame):
        basis = full_represent(unit)
        conjugated_basis = tuple(tuple(z.conjugate() for z in row) for row in basis)
        conjugated_row = tuple(z.conjugate() for z in row)
        sigma_basis = paravector(conjugated_row[0], scale(-1, conjugated_row[1:]))
        hermitian_basis = paravector(conjugated_row[0], conjugated_row[1:])
        adjugate_basis = paravector(row[0], scale(-1, row[1:]))
        check_matrix('basis-column sigma conjugation', sigma_basis,
                     mm(mm(ms(-1, s_matrix), conjugated_basis), s_matrix))
        check_matrix('basis-column Hermitian conjugation', hermitian_basis,
                     transpose(conjugated_basis))
        check_matrix('basis-column adjugation', adjugate_basis,
                     ((basis[1][1], -basis[0][1]), (-basis[1][0], basis[0][0])))
        sigma_from_basis = ma(sigma_from_basis, ms(coefficient.conjugate(), sigma_basis))
        hermitian_from_basis = ma(hermitian_from_basis, ms(coefficient.conjugate(), hermitian_basis))
        adjugate_from_basis = ma(adjugate_from_basis, ms(coefficient, adjugate_basis))
    check_matrix('field sigma conjugation in local basis', sigma_from_basis, sigma_direct)
    check_matrix('field Hermitian conjugation in local basis', hermitian_from_basis, hermitian_direct)
    check_matrix('field adjugation in local basis', adjugate_from_basis, adjugate_direct)

    product_scalar = original_q[0]*original_p[0] + dot(original_q[1:], original_p[1:])
    product_vector = add(add(scale(original_q[0], original_p[1:]), scale(original_p[0], original_q[1:])),
                         scale(1j, cross(original_q[1:], original_p[1:])))
    product_column = mv(inv_ft, (product_scalar,)+product_vector)
    check_matrix('full-basis refactored product', full_represent(product_column), mm(value, other))

    # Mixed interval, including a temporal element with a vector part and
    # spatial elements with scalar parts.
    increment = tuple(real() for _ in range(4))
    temporal = full_represent((1, 0, 0, 0))
    spatial = full_represent((0,)+increment[1:])
    temporal_adjugate = ((temporal[1][1], -temporal[0][1]), (-temporal[1][0], temporal[0][0]))
    mixed = coefficients(mm(temporal_adjugate, spatial))[0]
    dt_det = temporal[0][0]*temporal[1][1]-temporal[0][1]*temporal[1][0]
    dr_det = spatial[0][0]*spatial[1][1]-spatial[0][1]*spatial[1][0]
    check('mixed temporal-spatial interval', (dt_det*increment[0]**2+2*mixed*increment[0]+dr_det,),
          (dot(increment, mv(metric, increment)),))

    for direction in range(4):
        rate = tuple(tuple(entry() for _ in range(4)) for _ in range(4))
        # Include fixed coordinate increments as well as varying fields.
        q_rate = (0, 0, 0, 0) if direction % 2 == 0 else tuple(scalar() for _ in range(4))
        actual_column = add(mv(transpose(rate), q), mv(ft, q_rate))
        primed_column = add(q_rate, mv(transpose(mm(rate, inv)), q))
        check('full-basis time/spatial derivative', mv(inv_ft, actual_column), primed_column)
        epsilon = 1e-5

        def sample(e):
            result = mv(transpose(ma(frame, ms(e, rate))), add(q, scale(e, q_rate)))
            return paravector(result[0], result[1:])

        central = ms(1/(2*epsilon), ma(sample(epsilon), ms(-1, sample(-epsilon))))
        check_matrix('full-basis finite difference', central,
                     paravector(actual_column[0], actual_column[1:]), tolerance=2e-8)


# Exact curvature from metric first and second derivatives, independent of
# the paravector calculation. Convention: R^r_smn = d_m Gamma^r_ns - ... .
for a, ap, app in ((1, 0, 0), (2, 3, 5), (1, 1, 1), (3, 0, 2), (2, 1, 0)):
    a, ap, app = map(Q, (a, ap, app))
    g = (Q(1), -a*a, -a*a, -a*a)
    gi = tuple(1/x for x in g)

    def dg(k, i, j):
        return -2*a*ap if k == 0 and i == j and i > 0 else Q(0)

    def ddg(k, l, i, j):
        return -2*(ap*ap+a*app) if k == l == 0 and i == j and i > 0 else Q(0)

    def gamma(r, m, n):
        return gi[r]*(dg(m, r, n)+dg(n, r, m)-dg(r, m, n))/2

    def dgamma(k, r, m, n):
        inverse_term = sum(-gi[r]*gi[l]*dg(k, r, l)*(dg(m, l, n)+dg(n, l, m)-dg(l, m, n)) for l in range(4))
        metric_term = gi[r]*(ddg(k, m, r, n)+ddg(k, n, r, m)-ddg(k, r, m, n))
        return (inverse_term+metric_term)/2

    def riemann(r, s, m, n):
        return dgamma(m, r, n, s)-dgamma(n, r, m, s)+sum(
            gamma(r, m, l)*gamma(l, n, s)-gamma(r, n, l)*gamma(l, m, s) for l in range(4))

    # The concise transport law in the paper, checked against the metric
    # connection on a basis for each of its two vector arguments.
    for tangent, carried in product(eye(4), repeat=2):
        direct = tuple(-sum(gamma(r, m, n)*tangent[m]*carried[n]
                            for m, n in product(range(4), repeat=2))
                       for r in range(4))
        predicted = (-a*ap*dot(tangent[1:], carried[1:]),) + scale(
            -ap/a, add(scale(tangent[0], carried[1:]),
                       scale(carried[0], tangent[1:])))
        assert predicted == direct

        # Convert coordinate components to the measured scalar-vector parts.
        # Differentiating a*Y includes the time-dependent basis scale.
        physical = (carried[0],)+scale(a, carried[1:])
        physical_rate = (direct[0],)+add(scale(a, direct[1:]),
                                       scale(ap*tangent[0], carried[1:]))
        expected_physical = (-ap*dot(tangent[1:], physical[1:]),)+scale(
            -ap*physical[0], tangent[1:])
        assert physical_rate == expected_physical

    # Along a massive geodesic, the measured momentum redshifts as 1/a.
    mass, momentum = 1.2, (2., -1., .5)
    energy = sqrt(mass**2+dot(momentum, momentum))
    time_rate = energy/mass
    spatial_rate = scale(1/(mass*a), momentum)
    time_accel = -a*ap*dot(spatial_rate, spatial_rate)
    spatial_accel = scale(-2*ap*time_rate/a, spatial_rate)
    momentum_rate = scale(mass/time_rate,
                         add(scale(ap*time_rate, spatial_rate), scale(a, spatial_accel)))
    energy_rate = mass*time_accel/time_rate
    check('cosmological momentum redshift', momentum_rate, scale(-ap/a, momentum))
    check('massive measured-energy transport', (energy_rate,),
          (-ap/a*dot(momentum, momentum)/energy,))
    check('mass shell preserved during expansion',
          (2*energy*energy_rate-2*dot(momentum, momentum_rate),), (0,))

    # Arbitrary transported vectors and null photon momenta use t, not
    # proper time, in the new scalar-vector transport equations.
    velocity, carried_scalar, carried_vector = (Q(1, 3), Q(-2, 7), Q(1, 5)), Q(2), (Q(3), Q(-4), Q(5))
    scalar_rate = -ap/a*dot(velocity, carried_vector)
    vector_rate = scale(-ap/a*carried_scalar, velocity)
    assert carried_scalar*scalar_rate-dot(carried_vector, vector_rate) == 0
    photon_direction = (Q(2, 3), Q(-1, 3), Q(2, 3))
    photon_energy = Q(3)
    photon_momentum = scale(photon_energy, photon_direction)
    photon_energy_rate = -ap/a*dot(photon_direction, photon_momentum)
    photon_momentum_rate = scale(-ap/a*photon_energy, photon_direction)
    assert ap*photon_energy+a*photon_energy_rate == 0
    assert photon_energy*photon_energy_rate-dot(photon_momentum, photon_momentum_rate) == 0

    squared = sum((g[r]*riemann(r, s, m, n))**2*gi[r]*gi[s]*gi[m]*gi[n]
                  for r, s, m, n in product(range(4), repeat=4))
    assert squared == 12*((app/a)**2+(ap/a)**4)
    for i, j in product(range(1, 4), repeat=2):
        assert riemann(i, 0, j, 0) == -(app/a)*int(i == j)

print(f'PASS: {checks} algebra, derivative, wave, geodesic, and transport checks '
      'in 48 spatial and 32 real and 32 complex full spacetime matrix frames; '
      'exact transport, tidal, and curvature checks in 5 Robertson-Walker cases.')


def kerr_fields(position, spin, mass_length):
    p2, a2, ap = dot(position, position), dot(spin, spin), dot(spin, position)
    r2 = ((p2-a2)+sqrt((p2-a2)**2+4*ap**2))/2
    r = sqrt(r2)
    direction = scale(1/(r2+a2), add(add(scale(r, position), scale(-1, cross(spin, position))),
                                    scale(ap/r, spin)))
    profile = mass_length*r**3/(r**4+ap**2)
    return r, direction, profile


def kerr_metric(event, spin, mass_length):
    _, direction, profile = kerr_fields(event[1:], spin, mass_length)
    null = (1,)+direction
    return tuple(tuple(signature[i][j]-2*profile*null[i]*null[j]
                       for j in range(4)) for i in range(4))


def omega_kerr_fields(position, omega, mass_length):
    rho = 4*mass_length**2/(1+4*mass_length**2*dot(omega, omega))
    omega2, projection = dot(omega, omega), dot(omega, position)
    radial_term = dot(position, position)-rho**2*omega2
    epsilon = sqrt((radial_term+sqrt(radial_term**2+4*rho**2*projection**2))/2)
    direction = scale(1/(epsilon**2+rho**2*omega2),
        add(add(scale(epsilon, position), scale(-rho, cross(omega, position))),
            scale(rho**2*projection/epsilon, omega)))
    profile = mass_length*epsilon**3/(epsilon**4+rho**2*projection**2)
    return rho, epsilon, direction, profile


def omega_kerr_metric(event, omega, mass_length):
    _, _, direction, profile = omega_kerr_fields(event[1:], omega, mass_length)
    mixed = scale(-2*profile, direction)
    spatial = tuple(tuple(-int(i == j)-2*profile*direction[i]*direction[j]
                          for j in range(3)) for i in range(3))
    return ((1-2*profile,)+mixed,)+tuple((mixed[i],)+spatial[i] for i in range(3))


def determinant(matrix):
    # Independent permutation expansion, used only for 2x2 and 4x4 checks.
    total = 0
    for ordering in permutations(range(len(matrix))):
        term = (-1)**sum(ordering[i] > ordering[j]
                         for i in range(len(matrix)) for j in range(i+1, len(matrix)))
        for i, j in enumerate(ordering):
            term *= matrix[i][j]
        total += term
    return total


ks_start = checks
# Reparameterize Kerr by its outer-horizon angular velocity, including
# the Schwarzschild and extremal endpoints.
for fraction in (0, .01, .2, .8, .99, 1):
    mass_length = .7
    omega = scale(fraction/(2*mass_length), (2/3, -1/3, 2/3))
    rho = 4*mass_length**2/(1+4*mass_length**2*dot(omega, omega))
    spin = scale(rho, omega)
    horizon = rho/(2*mass_length)
    check('Kerr horizon angular velocity', scale(1/(2*mass_length*horizon), spin), omega)
    check('Kerr horizon quadratic',
          (horizon**2-2*mass_length*horizon+dot(spin, spin),), (0,))
    assert horizon >= mass_length-1e-12
    expected_horizon = mass_length if fraction == 1 else (
        mass_length+sqrt(mass_length**2-dot(spin, spin)))
    check('Kerr outer horizon branch', (horizon,), (expected_horizon,))
    position = (2.3, -.8, 1.2)
    radius, direction, profile = kerr_fields(position, spin, mass_length)
    projection = dot(omega, position)
    check('angular-velocity ellipsoidal radius',
          (radius**4-radius**2*(dot(position, position)-rho**2*dot(omega, omega))
           -rho**2*projection**2,), (0,))
    new_direction = scale(1/(radius**2+rho**2*dot(omega, omega)),
        add(add(scale(radius, position), scale(-rho, cross(omega, position))),
            scale(rho**2*projection/radius, omega)))
    new_profile = mass_length*radius**3/(radius**4+rho**2*projection**2)
    check('angular-velocity Kerr direction', new_direction, direction)
    check('angular-velocity Kerr profile', (new_profile,), (profile,))

# Verify the displayed omega-based blocks against the paravector basis,
# without constructing the reference metric from the same block formula.
axis = (2/3, -1/3, 2/3)
transverse = cross(axis, (0, 0, 1))
transverse = scale(1/sqrt(dot(transverse, transverse)), transverse)
for fraction, radial_factor, (sine, cosine) in product(
        (0, .01, .2, .8, .99, 1), (.8, 1, 2), ((0, 1), (1, 0), (.8, .6))):
    mass_length = .7
    omega = scale(fraction/(2*mass_length), axis)
    rho = 4*mass_length**2/(1+fraction**2)
    horizon = rho/(2*mass_length)
    epsilon = radial_factor*horizon
    spin = scale(rho, omega)
    position = add(scale(sine*sqrt(epsilon**2+dot(spin, spin)), transverse),
                   scale(cosine*epsilon, axis))
    actual_rho, actual_epsilon, direction, profile = omega_kerr_fields(position, omega, mass_length)
    check('omega radius on oblate surfaces', (actual_rho, actual_epsilon), (rho, epsilon))
    unit_formula = (dot(position, position)+rho**2*dot(omega, position)**2/epsilon**2)/(
        epsilon**2+rho**2*dot(omega, omega))
    check('omega unit-direction identity', (dot(direction, direction), unit_formula), (1, 1))
    basis = (paravector(1-profile, scale(profile, direction)),)+tuple(
        paravector(-profile*direction[i], add(eye(3)[i], scale(profile*direction[i], direction)))
        for i in range(3))
    # Polarize the 2x2 Pauli determinant to extract all four metric blocks.
    from_basis = tuple(tuple((determinant(ma(left, right))-determinant(left)-determinant(right))/2
                            for right in basis) for left in basis)
    event = (0,)+position
    metric = omega_kerr_metric(event, omega, mass_length)
    check('omega g_tt from Pauli determinant', (metric[0][0],), (from_basis[0][0],))
    check('omega g_tr from Pauli determinant', metric[0][1:], from_basis[0][1:])
    check_matrix('omega g_rr from Pauli determinant', tuple(row[1:] for row in metric[1:]),
                 tuple(row[1:] for row in from_basis[1:]))
    check_matrix('omega blocks equal standard Kerr metric', metric, kerr_metric(event, spin, mass_length))
    check('Kerr metric determinant', (determinant(metric),), (-1,))
    increment = (.4, -.3, .2, .1)
    displacement = ((0j, 0j), (0j, 0j))
    for component, element in zip(increment, basis):
        displacement = ma(displacement, ms(component, element))
    interval = increment[0]**2-dot(increment[1:], increment[1:])-2*profile*(
        increment[0]+dot(direction, increment[1:]))**2
    check('omega ds squared from determinant and blocks',
          (determinant(displacement), dot(increment, mv(metric, increment))), (interval, interval))
    if radial_factor == 1:
        generator = (1,)+cross(omega, position)
        check('horizon angular-velocity generator is null', (dot(generator, mv(metric, generator)),), (0,))
        gradient = scale(1/(epsilon*(2*epsilon**2-dot(position, position)+dot(spin, spin))),
                         add(scale(epsilon**2, position), scale(dot(spin, position), spin)))
        normal = (0,)+gradient
        check('outer horizon normal is null', (dot(normal, mv(inverse(metric), normal)),), (0,))

for case in range(80):
    position = (2+real(), 1+real(), .7+real())
    spin = (0, 0, real()) if case % 2 else tuple(real() for _ in range(3))
    radius, direction, profile = kerr_fields(position, spin, .7)
    check('Kerr unit direction', (dot(direction, direction),), (1,))
    if case % 2:
        x, y, z = position
        a = spin[2]
        # Independently quoted Cartesian form, including the spin orientation.
        check('Kerr Cartesian source formula', direction,
              ((radius*x+a*y)/(radius**2+a*a),
               (radius*y-a*x)/(radius**2+a*a), z/radius))

    radial, spherical_u, spherical_f = kerr_fields(position, (0, 0, 0), .7)
    check('Schwarzschild limit', (radial, spherical_f)+spherical_u,
          (sqrt(dot(position, position)), .7/sqrt(dot(position, position)))
          + scale(1/sqrt(dot(position, position)), position))

    # Include large and negative profiles: invertibility is not a weak-field claim.
    f = profile if case % 2 else 4*real()
    null = (1,)+direction
    raised_null = (1,)+scale(-1, direction)
    frame = tuple(tuple(int(i == j)-f*null[i]*raised_null[j]
                        for j in range(4)) for i in range(4))
    inverse_frame = tuple(tuple(int(i == j)+f*null[i]*raised_null[j]
                                for j in range(4)) for i in range(4))
    check_matrix('Kerr-Schild inverse basis', mm(frame, inverse_frame), eye(4))
    metric = tuple(tuple(signature[i][j]-2*f*null[i]*null[j]
                         for j in range(4)) for i in range(4))
    check_matrix('Kerr-Schild metric from basis', mm(mm(frame, signature), transpose(frame)), metric)
    increment = tuple(real() for _ in range(4))
    original = mv(transpose(frame), increment)
    represented = paravector(original[0], original[1:])
    determinant = represented[0][0]*represented[1][1]-represented[0][1]*represented[1][0]
    expected = increment[0]**2-dot(increment[1:], increment[1:])-2*f*dot(null, increment)**2
    check('Kerr-Schild Pauli determinant', (determinant,), (expected,))


def finite_difference_geometry(metric_at, event, step):
    """Coordinate connection and Ricci tensor, independent of the basis rules."""
    metric = metric_at(event)
    inv = inverse(metric)

    def sample(offsets):
        point = list(event)
        for coordinate, amount in offsets:
            point[coordinate] += step*amount
        return metric_at(tuple(point))

    first = tuple(ms(1/(2*step), ma(sample(((k, 1),)), ms(-1, sample(((k, -1),)))))
                  for k in range(4))
    def mixed_second(k, l):
        if k == l:
            value = ma(ma(sample(((k, 1),)), sample(((k, -1),))), ms(-2, metric))
            return ms(1/step**2, value)
        same_sign = ma(sample(((k, 1), (l, 1))), sample(((k, -1), (l, -1))))
        opposite_sign = ma(sample(((k, 1), (l, -1))), sample(((k, -1), (l, 1))))
        return ms(1/(4*step**2), ma(same_sign, ms(-1, opposite_sign)))

    second = tuple(tuple(mixed_second(k, l) for l in range(4)) for k in range(4))
    gamma = connection(metric, first)
    inverse_rate = tuple(ms(-1, mm(mm(inv, derivative), inv)) for derivative in first)

    def gamma_rate(k, r, m, n):
        return sum(inverse_rate[k][r][l]*(first[m][l][n]+first[n][l][m]-first[l][m][n])
                   + inv[r][l]*(second[k][m][l][n]+second[k][n][l][m]-second[k][l][m][n])
                   for l in range(4))/2

    ricci = tuple(tuple(sum(gamma_rate(r, r, n, m)-gamma_rate(n, r, r, m)
                    + sum(gamma[r][r][l]*gamma[l][n][m]-gamma[r][n][l]*gamma[l][r][m]
                          for l in range(4)) for r in range(4))
                        for n in range(4)) for m in range(4))
    return gamma, ricci


for spin, position in (((0, 0, 0), (3, 1, 2)), ((0, 0, .4), (2, -1, .8)),
                       ((.2, -.1, .4), (1.2, .5, .3)), ((0, 0, .7), (3, 2, -1))):
    event = (0,)+position
    metric_at = lambda event: kerr_metric(event, spin, .7)
    _, coarse = finite_difference_geometry(metric_at, event, .002)
    _, fine = finite_difference_geometry(metric_at, event, .001)
    extrapolated = ms(1/3, ma(ms(4, fine), ms(-1, coarse)))
    check('Kerr/Schwarzschild vacuum Ricci', flatten(extrapolated), (0,)*16, tolerance=3e-7)

for fraction in (0, .6, 1):
    omega = scale(fraction/(2*.7), axis)
    metric_at = lambda event: omega_kerr_metric(event, omega, .7)
    _, coarse = finite_difference_geometry(metric_at, (0, 2, -1, .8), .002)
    _, fine = finite_difference_geometry(metric_at, (0, 2, -1, .8), .001)
    extrapolated = ms(1/3, ma(ms(4, fine), ms(-1, coarse)))
    check('omega metric vacuum Ricci', flatten(extrapolated), (0,)*16, tolerance=3e-7)

# Check the Schwarzschild radial reduction for general profiles, not just
# the solution f=C/r. In spherical coordinates the two Ricci factors are
# f''+2f'/r and 2(f+r*f')/r**2; (r*f)'=0 makes both vanish.
# Compare their Cartesian tensor to independent metric finite differences.
radial_event = (0, 3, 1, 2)
radial_position = radial_event[1:]
radial_radius = sqrt(dot(radial_position, radial_position))
radial_direction = scale(1/radial_radius, radial_position)
for terms in (((.7, -1),), ((.1, 0),), ((.02, 2),),
              ((.1, -2),), ((.7, -1), (.03, 0), (.01, 2))):
    def radial_profile_metric(event):
        radius = sqrt(dot(event[1:], event[1:]))
        profile = sum(coefficient*radius**power for coefficient, power in terms)
        null = (1,)+scale(1/radius, event[1:])
        return tuple(tuple(signature[i][j]-2*profile*null[i]*null[j]
                           for j in range(4)) for i in range(4))

    profile = sum(c*radial_radius**p for c, p in terms)
    radial_factor = sum(c*p*(p+1)*radial_radius**(p-2) for c, p in terms)
    angular_factor = sum(2*c*(p+1)*radial_radius**(p-2) for c, p in terms)
    expected = ((-(1-2*profile)*radial_factor,)
                + scale(2*profile*radial_factor, radial_direction),)
    expected += tuple(
        (2*profile*radial_factor*radial_direction[i],)
        + tuple((1+2*profile)*radial_factor*radial_direction[i]*radial_direction[j]
                + angular_factor*(int(i == j)-radial_direction[i]*radial_direction[j])
                for j in range(3)) for i in range(3))
    _, coarse = finite_difference_geometry(radial_profile_metric, radial_event, .002)
    _, fine = finite_difference_geometry(radial_profile_metric, radial_event, .001)
    extrapolated = ms(1/3, ma(ms(4, fine), ms(-1, coarse)))
    check_matrix('spherical vacuum radial reduction', extrapolated, expected, tolerance=3e-7)

# A non-vacuum radial profile must not pass the same Ricci check.
def constant_profile_metric(event):
    position = event[1:]
    null = (1,)+scale(1/sqrt(dot(position, position)), position)
    return tuple(tuple(signature[i][j]-.2*null[i]*null[j] for j in range(4)) for i in range(4))

_, nonvacuum = finite_difference_geometry(constant_profile_metric, (0, 3, 1, 2), .001)
assert max(abs(value) for value in flatten(nonvacuum)) > 1e-3

# For a particle initially at rest, the exact coordinate acceleration tends
# to the Newtonian inverse-square law as the source strength tends to zero.
event = (0, 3, 1, 2)
position = event[1:]
radius = sqrt(dot(position, position))
errors = []
for strength in (1e-3, 5e-4):
    gamma, _ = finite_difference_geometry(lambda event: kerr_metric(event, (0, 0, 0), strength),
                                         event, .001)
    acceleration = tuple(-gamma[i][0][0]/strength for i in range(1, 4))
    newtonian = scale(-1/radius**3, position)
    errors.append(max(abs(x-y) for x, y in zip(acceleration, newtonian)))
assert errors[1] < .51*errors[0] and errors[1] < 2e-5, errors

print(f'PASS: {checks-ks_start} Kerr-Schild, spin, determinant, and vacuum checks; '
      'spherical radial reduction, non-vacuum control, and Newtonian-limit convergence.')
