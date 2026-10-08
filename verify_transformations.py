"""Independent checks of the direct-action conventions in the paper.

Run ``python3 verify_transformations.py`` (requires SymPy and NumPy).
Pauli matrices provide an independent, faithful implementation of paravectors:

    R = Q = exp(-i theta u.sigma/2), L = exp(-eta u.sigma/2),
    X' = T X T^H, F' = T F T^-1, partial' = (T^-1)^H partial T^-1,
    Psi1' = (T^-1)^H Psi1, Psi2' = T Psi2, psi' = R psi.

Exact polynomial checks cover arbitrary axes; rational noncommuting examples
check compositions, differential covariance, and moving frames. Seeded numeric
checks exercise independent trigonometric and hyperbolic parameter values.
These checks validate the convention, independently of the manuscript source.
"""

import time

import numpy as np
import sympy as s


started = time.perf_counter()
checks = 0
one, zero = s.eye(2), s.zeros(2)
sigma = (s.Matrix(((0, 1), (1, 0))),
         s.Matrix(((0, -s.I), (s.I, 0))), s.diag(1, -1))


def vector(v):
    return sum((v[k] * sigma[k] for k in range(3)), zero)


def sc(z):
    return s.trace(z) / 2


def adj(z):
    return s.trace(z) * one - z


def star(z):
    return adj(z.H)


def norm2(z):
    return sc(z.H * z)


def inner(a, b):
    return sc(a.H * b)


def comm(a, b):
    return a * b - b * a


def coeff(z):
    return s.Matrix((sc(z), *(sc(item * z) for item in sigma)))


def W(z):
    return s.BlockMatrix(((zero, adj(z)), (z, zero))).as_explicit()


def check(name, expression, reduce=s.expand):
    global checks
    entries = list(expression) if isinstance(expression, s.MatrixBase) else [expression]
    residuals = [reduce(item) for item in entries]
    assert all(item == 0 for item in residuals), (name, residuals)
    checks += 1
    print(f"PASS {name}")


def counterexample(name, expression):
    global checks
    entries = list(expression) if isinstance(expression, s.MatrixBase) else [expression]
    assert any(s.cancel(item) != 0 for item in entries), name
    checks += 1
    print(f"PASS {name}")


ux, uy, uz, c, h = s.symbols("ux uy uz c h", real=True)
u = s.Matrix((ux, uy, uz))
t0, x, y, z = s.symbols("t0 x y z", real=True)
r = s.Matrix((x, y, z))
N = vector(u)
X = t0 * one + vector(r)


def constrained(expression, circular=False):
    numerator = s.together(expression).as_numer_denom()[0]
    relation = c**2 + h**2 - 1 if circular else c**2 - h**2 - 1
    result = s.rem(s.expand(numerator), relation, c)
    return s.rem(result, uz**2 + ux**2 + uy**2 - 1, uz).expand()


def circular(expression):
    return constrained(expression, circular=True)


R = c * one - s.I * h * N
L, Li = c * one - h * N, c * one + h * N
gamma, gamma_beta = c**2 + h**2, 2 * c * h
rotated_r = ((c**2 - h**2) * r + 2 * h**2 * u.dot(r) * u
             + 2 * c * h * u.cross(r))
boosted_r = r + (gamma - 1) * u.dot(r) * u - gamma_beta * t0 * u
boosted_X = (gamma * t0 - gamma_beta * u.dot(r)) * one + vector(boosted_r)
check("arbitrary-axis involution", N**2 - one, constrained)
check("rotor is unitary", R.H * R - one, circular)
check("rotor determinant", R.det() - 1, circular)
check("direct rotor Rodrigues has positive axis cross vector", R * X * R.H - t0 * one
      - vector(rotated_r), circular)
check("booster inverse", L * Li - one, constrained)
check("booster determinant", L.det() - 1, constrained)
check("booster is Hermitian", L.H - L)
check("booster positive norm is gamma", norm2(L) - gamma, constrained)
check("direct booster gives the standard negative-velocity frame terms", L * X * L - boosted_X, constrained)
check("boost preserves the spacetime determinant", boosted_X.det() - X.det(), constrained)
check("inverse booster square is the positive-velocity tangent", Li**2 - gamma * one
      - gamma_beta * N, constrained)
check("booster square contains the negative-velocity frame terms", L**2 - gamma * one
      + gamma_beta * N, constrained)
check("direct booster brings that tangent to rest", L * Li**2 * L - one, constrained)
check("transverse event components remain fixed", L * vector(r - u.dot(r) * u) * L
      - vector(r - u.dot(r) * u), constrained)

E = s.Matrix(s.symbols("Ex Ey Ez", real=True))
B = s.Matrix(s.symbols("Bx By Bz", real=True))
F = vector(E + s.I * B)
Ep = gamma * E - (gamma - 1) * u.dot(E) * u + gamma_beta * u.cross(B)
Bp = gamma * B - (gamma - 1) * u.dot(B) * u - gamma_beta * u.cross(E)
check("boosted electromagnetic components have standard signs", L * F * Li
      - vector(Ep + s.I * Bp), constrained)
check("boosted field retains zero scalar part", sc(L * F * Li), constrained)
check("electric minus magnetic square invariant", Ep.dot(Ep) - Bp.dot(Bp)
      - E.dot(E) + B.dot(B), constrained)
check("electric dot magnetic invariant", Ep.dot(Bp) - E.dot(B), constrained)

# Hamilton quaternion units and the positive-angle active quaternion Q = R.
q = tuple(-s.I * item for item in sigma)
check("Hamilton quaternion multiplication", q[0] * q[1] - q[2])
check("Hamilton quaternion units square to minus one", q[0]**2 + one)
a = s.symbols("a", real=True)
A = s.Matrix(s.symbols("Ax Ay Az", real=True))
Q = a * one - s.I * vector(A)
quaternion_vector = ((a**2 - A.dot(A)) * r + 2 * A.dot(r) * A + 2 * a * A.cross(r))
check("quaternion vector sandwich gives positive cross product", Q * (-s.I * vector(r)) * Q.H
      + s.I * vector(quaternion_vector))
check("quaternion Q is exactly the rotor R", Q.subs(dict(zip(A, h * u))).subs(a, c) - R)
check("quaternion Rodrigues equals paravector Rodrigues", quaternion_vector.subs(
      dict(zip(A, h * u))).subs(a, c) - rotated_r, circular)

# Exact noncollinear examples include a composite that is neither Hermitian nor unitary.
axis1 = s.Matrix((s.Rational(2, 3), s.Rational(1, 3), s.Rational(2, 3)))
axis2 = s.Matrix((s.Rational(-2, 7), s.Rational(3, 7), s.Rational(6, 7)))
rotor = s.Rational(4, 5) * one - s.I * s.Rational(3, 5) * vector(axis1)
booster = s.Rational(5, 4) * one - s.Rational(3, 4) * vector(axis2)
composite = (booster * rotor).expand()
Z = (2 + s.I) * one + vector(s.Matrix((1 + 2 * s.I, -3 + s.I, 4 - s.I)))
Y = (1 - 2 * s.I) * one + vector(s.Matrix((2 - s.I, 1 + s.I, -1 + 3 * s.I)))
event = 7 * one + vector(s.Matrix((2, -3, 5)))
field = vector(s.Matrix((2 - s.I, -3 + 4 * s.I, 5 + 2 * s.I)))
psi = s.Matrix(((1 + s.I, 0), (2 - s.I, 0)))
chi = s.Matrix(((2 - 3 * s.I, 0), (-1 + 2 * s.I, 0)))


def event_map(T, Z):
    return (T * Z * T.H).expand()


def field_map(T, Z):
    return (T * Z * adj(T)).expand()


def component_map(T):
    return s.Matrix.hstack(*(coeff(event_map(T, b)) for b in (one, *sigma)))


# Arbitrary complex gradient values and real coordinate increments test the
# differential pairing without assuming a particular scalar field.
gradient_values = s.Matrix(s.symbols("gradient0:4", complex=True))
coordinate_increment = s.Matrix(s.symbols("increment0:4", real=True))
gradient_paravector = gradient_values[0] * one + vector(gradient_values[1:, :])
increment_paravector = coordinate_increment[0] * one + vector(coordinate_increment[1:, :])
real_four_values = s.Matrix(s.symbols("four0:4", real=True))
real_four_paravector = real_four_values[0] * one + vector(real_four_values[1:, :])

for name, T in (("rotation", rotor), ("boost", booster), ("composite", composite)):
    Ti = adj(T)
    Zp, Yp = event_map(T, Z), event_map(T, Y)
    Fp = field_map(T, field)
    M = component_map(T)
    transformed_gradient = Ti.H * gradient_paravector * Ti
    check(f"{name}: derivative coefficients match inverse-transpose coordinate Jacobian",
          coeff(transformed_gradient) - M.inv().T * gradient_values)
    check(f"{name}: complex total differential is unchanged",
          sc(transformed_gradient * event_map(T, increment_paravector))
          - gradient_values.dot(coordinate_increment))
    generic_field_prime = field_map(T, F)
    electric_prime = s.re(coeff(generic_field_prime)[1:, :])
    magnetic_prime = s.im(coeff(generic_field_prime)[1:, :])
    check(f"{name}: energy and Poynting flux from transformed field components",
          generic_field_prime * generic_field_prime.H / 2
          - (electric_prime.dot(electric_prime) + magnetic_prime.dot(magnetic_prime)) * one / 2
          - vector(electric_prime.cross(magnetic_prime)))
    generic_four_prime = event_map(T, real_four_paravector)
    check(f"{name}: Lorentz force covariance for arbitrary real four-vector",
          (generic_field_prime * generic_four_prime + generic_four_prime * generic_field_prime.H) / 2
          - event_map(T, (F * real_four_paravector + real_four_paravector * F.H) / 2))
    check(f"{name}: phase invariance for arbitrary real momentum and event",
          sc(star(generic_four_prime) * event_map(T, X)) - sc(star(real_four_paravector) * X))
    check(f"{name}: event determinant", Zp.det() - Z.det())
    check(f"{name}: complex bilinear Minkowski pairing", sc(adj(Zp) * Yp) - sc(adj(Z) * Y))
    check(f"{name}: complex Hermitian Minkowski pairing", sc(star(Zp) * Yp) - sc(star(Z) * Y))
    check(f"{name}: component Lorentz metric", M.T * s.diag(1, -1, -1, -1) * M
          - s.diag(1, -1, -1, -1))
    check(f"{name}: field determinant", Fp.det() - field.det())
    check(f"{name}: similarity trace and scalar", sc(field_map(T, Z)) - sc(Z))
    check(f"{name}: similarity bilinear trace product", sc(field_map(T, Z) * field_map(T, Y))
          - sc(Z * Y))
    psip, chip = T * psi, T * chi
    weight = T.H * event * T
    check(f"{name}: matrix-element pullback", inner(chip, event * psip) - inner(chi, weight * psi))
    check(f"{name}: normalized expectation pullback", inner(psip, event * psip) / norm2(psip)
          - inner(psi, weight * psi) / inner(psi, T.H * T * psi), s.cancel)
    check(f"{name}: co-transformed weights preserve matrix elements",
          inner(chip, Ti.H * event * Ti * psip) - inner(chi, event * psi))
    if name != "rotation":
        counterexample(f"{name}: event trace need not be invariant", s.trace(Zp) - s.trace(Z))
        counterexample(f"{name}: event Euclidean norm need not be invariant", norm2(Zp) - norm2(Z))
        counterexample(f"{name}: event positive inner product need not be invariant", inner(Zp, Yp) - inner(Z, Y))
        counterexample(f"{name}: field Euclidean norm need not be invariant", norm2(Fp) - norm2(field))
        counterexample(f"{name}: spinor norm need not be invariant", norm2(psip) - norm2(psi))

check("rotation preserves event trace", s.trace(event_map(rotor, Z)) - s.trace(Z))
check("rotation preserves coefficient norm", norm2(event_map(rotor, Z)) - norm2(Z))
check("rotation preserves positive inner product", inner(event_map(rotor, Z), event_map(rotor, Y)) - inner(Z, Y))
check("rotation preserves spinor norm", norm2(rotor * psi) - norm2(psi))
check("active quantum rotor expectation pullback", inner(rotor * psi, event * rotor * psi)
      - inner(psi, rotor.H * event * rotor * psi))
polarization = s.Matrix(tuple(inner(psi, item * psi) / norm2(psi) for item in sigma))
rotated_polarization = (s.Rational(7, 25) * polarization
                       + s.Rational(18, 25) * axis1.dot(polarization) * axis1
                       + s.Rational(24, 25) * axis1.cross(polarization))
check("quantum spin expectation has positive Rodrigues angle", s.Matrix(tuple(
    inner(rotor * psi, item * rotor * psi) / norm2(psi) for item in sigma)) - rotated_polarization)
counterexample("composite event action and expectation pullback have different orders",
               event_map(composite, event) - composite.H * event * composite)
counterexample("a boost similarity cannot transform spacetime", field_map(booster, event).H
               - field_map(booster, event))
counterexample("a field congruence creates a scalar component", sc(event_map(booster, field)))
counterexample("left and right boost multiplication are different", booster * event - event * booster)
counterexample("left boost multiplication does not preserve a real event", (booster * event).H - booster * event)
counterexample("right boost multiplication does not preserve a real event", (event * booster).H - event * booster)
counterexample("left multiplication is not a Lorentz event boost", booster * event - event_map(booster, event))
check("one-sided multiplication nevertheless preserves determinant", (booster * event).det() - event.det())
check("noncommuting transverse vectors prevent combining the boost factors",
      comm(booster, vector(axis1)) - s.Rational(-3, 2) * s.I * vector(axis2.cross(axis1)))
check("a boosted future timelike event has positive scalar part", s.sign(sc(event_map(booster, event))) - 1)
check("its coefficient norm is positive", s.sign(norm2(event_map(booster, event))) - 1)
counterexample("positive coefficient norm differs from the invariant determinant",
               norm2(event_map(booster, event)) - event.det())
check("successive event maps use T_total=T2 T1", event_map(booster, event_map(rotor, event))
      - event_map(composite, event))
check("successive field maps use T_total=T2 T1", field_map(booster, field_map(rotor, field))
      - field_map(composite, field))
check("component maps compose in physical order", component_map(composite)
      - component_map(booster) * component_map(rotor))

boosted_generic_field = field_map(booster, F)
check("explicit boosted energy-flux product retains the inner inverse-boost square",
      boosted_generic_field * boosted_generic_field.H / 2
      - booster * F * adj(booster)**2 * F.H * booster / 2)

# W has adj(Z) above Z; its energy-momentum input is a four-vector.
# Conventional gamma^j = -W(sigma_j), with the physical Dirac blocks unchanged.
T, Ti = composite, adj(composite)
spin = s.diag(Ti.H, T)
gamma_matrices = (W(one), *(-W(item) for item in sigma))
eta = s.diag(1, -1, -1, -1)
check("Dirac gamma anticommutators", s.Matrix.vstack(*(g * k + k * g - 2 * eta[i, j] * s.eye(4)
      for i, g in enumerate(gamma_matrices) for j, k in enumerate(gamma_matrices))))
k = s.symbols("k0:4", real=True)
momentum_symbol = k[0] * one + vector(s.Matrix(k[1:]))
new_symbol = T * momentum_symbol * T.H
mass, charge = s.symbols("m q", real=True)
potential = event
new_potential = event_map(T, potential)
D = W(momentum_symbol) - mass * s.eye(4)
Dp = W(new_symbol) - mass * s.eye(4)
check("Dirac principal symbol covariance in every component", Dp * spin - spin * D)
check("Dirac minimal potential covariance in every component",
      (W(new_symbol - charge * new_potential) - mass * s.eye(4)) * spin
      - spin * (W(momentum_symbol - charge * potential) - mass * s.eye(4)))
check("opposite Dirac mass signs factor Klein-Gordon", (W(s.I * momentum_symbol) - mass * s.eye(4))
      * (W(s.I * momentum_symbol) + mass * s.eye(4))
      + (momentum_symbol.det() + mass**2) * s.eye(4))
check("Dirac spin representation composes with T2 T1", spin
      - s.diag(adj(booster).H, booster) * s.diag(adj(rotor).H, rotor))
check("opposite Dirac entries preserve their mixed inner product",
      inner(Ti.H * psi, T * chi) - inner(psi, chi))
current = adj(psi * psi.H) + chi * chi.H
new_current = adj(Ti.H * psi * psi.H * Ti) + T * chi * chi.H * T.H
check("Dirac probability current is a Lorentz four-vector", new_current - event_map(T, current))

# An explicit chain rule checks derivatives at corresponding spacetime events.
time_coordinate, cx, cy, cz = s.symbols("t cx cy cz", real=True)
coordinates = (time_coordinate, cx, cy, cz)
M = component_map(T)
old_coordinates = M.inv() * s.Matrix(coordinates)
substitution = dict(zip(coordinates, old_coordinates))


def at_old(expression):
    return expression.subs(substitution, simultaneous=True).expand()


def differential(expression, sign=1):
    return expression.diff(time_coordinate) + sign * sum((sigma[j] * expression.diff(coordinates[j + 1])
                                                         for j in range(3)), s.zeros(*expression.shape))


scalar_test = (time_coordinate**4 + time_coordinate * cx * cy + cy * cz**2
               + s.I * (cx**3 * time_coordinate + cy**2 * cz + cz**4)) * one
check("coordinate chain rule partial'=(T^-1)^H partial T^-1", differential(at_old(scalar_test))
      - Ti.H * at_old(differential(scalar_test)) * Ti)


def wave_operator(expression):
    return expression.diff(time_coordinate, 2) - sum(
        (expression.diff(coordinate, 2) for coordinate in coordinates[1:]),
        s.zeros(*expression.shape))


check("wave operator chain rule on a complex field with varying second derivatives",
      wave_operator(at_old(scalar_test)) - at_old(wave_operator(scalar_test)))
polynomial_field = (time_coordinate**2 * field + cx * cy * vector(axis1) + cz**2 * vector(axis2)).expand()
field_prime = field_map(T, at_old(polynomial_field))
check("Maxwell differential covariance", differential(field_prime)
      - Ti.H * at_old(differential(polynomial_field)) * Ti)
psi1 = (time_coordinate**2 + cx * cy) * psi + cz * chi
psi2 = (time_coordinate * cx + cz**2) * chi + cy * psi
psi1p, psi2p = Ti.H * at_old(psi1), T * at_old(psi2)
check("first Dirac equation differential covariance", s.I * differential(psi2p) - mass * psi1p
      - Ti.H * at_old(s.I * differential(psi2) - mass * psi1))
check("second Dirac equation differential covariance", s.I * differential(psi1p, -1) - mass * psi2p
      - T * at_old(s.I * differential(psi1, -1) - mass * psi2))

# Varying real potentials exercise potential derivatives in coupled equations.
variable_potential = ((time_coordinate * cx + cy * cz) * one
                      + vector(s.Matrix((cx * cy, cy * cz + time_coordinate * cy,
                                         cz * time_coordinate + cx**2))))
primed_potential = event_map(T, at_old(variable_potential))
raw_potential_derivative = differential(variable_potential)
raw_primed_derivative = differential(primed_potential)
check("potential scalar divergence is Lorentz invariant",
      sc(raw_primed_derivative) - at_old(sc(raw_potential_derivative)))
derived_field = star(raw_potential_derivative - sc(raw_potential_derivative) * one)
derived_primed_field = star(raw_primed_derivative - sc(raw_primed_derivative) * one)
check("field derived from transformed potential obeys similarity",
      derived_primed_field - field_map(T, at_old(derived_field)))
derived_current = star(differential(derived_field))
check("potential-derived Maxwell source is real", derived_current.H - derived_current)
check("Maxwell equation with transformed potential and source",
      differential(derived_primed_field) - star(event_map(T, at_old(derived_current))))


def mechanical_component(expression, potential, index):
    coefficient = coeff(potential)[index]
    derivative_sign = 1 if index == 0 else -1
    return (derivative_sign * s.I * expression.diff(coordinates[index])
            - charge * coefficient * expression).expand()


def coupled_kg(expression, potential):
    squared = [mechanical_component(mechanical_component(expression, potential, index),
                                    potential, index) for index in range(4)]
    return (-squared[0] + sum(squared[1:], s.zeros(*expression.shape))
            + mass**2 * expression).expand()


coupled_scalar_test = (time_coordinate * cx + cy**2 + s.I * cz * time_coordinate) * one
check("coupled Klein-Gordon covariance including potential derivatives",
      coupled_kg(at_old(coupled_scalar_test), primed_potential)
      - at_old(coupled_kg(coupled_scalar_test, variable_potential)))
check("first coupled Dirac equation with varying potentials",
      s.I * differential(psi2p) - charge * star(primed_potential) * psi2p - mass * psi1p
      - Ti.H * at_old(s.I * differential(psi2)
                     - charge * star(variable_potential) * psi2 - mass * psi1))
check("second coupled Dirac equation with varying potentials",
      s.I * differential(psi1p, -1) - charge * primed_potential * psi1p - mass * psi2p
      - T * at_old(s.I * differential(psi1, -1)
                  - charge * variable_potential * psi1 - mass * psi2))

# Moving frames: the signs follow by differentiating the actual actions.
tau = s.symbols("tau", real=True)
varying = (one + tau * s.Matrix(((0, 1), (0, 0)))) * (one + tau**2 * s.Matrix(((0, 0), (1, 0))))
inverse = adj(varying)
omega = varying.diff(tau) * inverse
moving_field = field + tau * vector(axis1) + tau**2 * s.I * vector(axis2)
moving_event = event + tau * one + tau**2 * vector(axis1)
Fp = field_map(varying, moving_field)
Xp = event_map(varying, moving_event)
check("moving generator has zero scalar", sc(omega))
check("moving inverse derivative", inverse.diff(tau) + inverse * omega)
check("moving field first derivative has plus commutator", Fp.diff(tau)
      - field_map(varying, moving_field.diff(tau)) - comm(omega, Fp))
check("moving field second derivative signs", Fp.diff(tau, 2)
      - field_map(varying, moving_field.diff(tau, 2)) - 2 * comm(omega, Fp.diff(tau))
      - comm(omega.diff(tau), Fp) + comm(omega, comm(omega, Fp)))
check("moving event derivative has two positive frame terms", Xp.diff(tau)
      - event_map(varying, moving_event.diff(tau)) - omega * Xp - Xp * omega.H)

# Numerical checks use arbitrary axes and signed angles/rapidities.
rng = np.random.default_rng(20260923)
pauli = np.array([np.array(item).astype(complex) for item in sigma])
identity = np.eye(2, dtype=complex)


def nv(v):
    return np.einsum("i,ijk->jk", v, pauli)


for index in range(8):
    axis = rng.normal(size=3)
    axis /= np.linalg.norm(axis)
    position, electric, magnetic = rng.normal(size=(3, 3))
    t = rng.normal()
    angle, rapidity = rng.uniform(-3, 3, size=2)
    nr = np.cos(angle / 2) * identity - 1j * np.sin(angle / 2) * nv(axis)
    nl = np.cosh(rapidity / 2) * identity - np.sinh(rapidity / 2) * nv(axis)
    nli = np.linalg.inv(nl)
    cosine, sine = np.cos(angle), np.sin(angle)
    g, gb = np.cosh(rapidity), np.sinh(rapidity)
    want_rotation = (cosine * position + (1 - cosine) * axis.dot(position) * axis
                     + sine * np.cross(axis, position))
    want_event = ((g * t - gb * axis.dot(position)) * identity
                  + nv(position + (g - 1) * axis.dot(position) * axis - gb * t * axis))
    want_e = g * electric - (g - 1) * axis.dot(electric) * axis + gb * np.cross(axis, magnetic)
    want_b = g * magnetic - (g - 1) * axis.dot(magnetic) * axis - gb * np.cross(axis, electric)
    np.testing.assert_allclose(nr @ nv(position) @ nr.conj().T, nv(want_rotation), atol=1e-12)
    np.testing.assert_allclose(nl @ (t * identity + nv(position)) @ nl.conj().T, want_event, atol=1e-12)
    np.testing.assert_allclose(nl @ nv(electric + 1j * magnetic) @ nli, nv(want_e + 1j * want_b), atol=1e-12)
    checks += 3
print("PASS 24 seeded numerical rotation, boost, and electromagnetic component checks")
print(f"\n{checks} checks passed in {time.perf_counter() - started:.2f} s.")
