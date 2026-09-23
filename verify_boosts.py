"""Historical checks of the earlier adjoint-factor parameterization.

Run with ``python3 verify_boosts.py`` (requires SymPy).  Ordinary Pauli
matrices are used only as an independent implementation of the algebra.
The current manuscript uses direct actions, checked independently by
``verify_transformations.py``. This retained suite uses historical factor
names and actions:

    L = exp(-theta*u.sigma/2), R = exp(-i*theta*u.sigma/2),
    X' = T^H X T, F' = T^H F (T^H)^-1,
    partial' = T^-1 partial (T^-1)^H.

Its historical frame factor K is related to the current direct frame factor
T by K = T^H.  Hermitian boost factors agree; the historical rotation R labels
the opposite rotation angle.  The equations below remain valid in that
parameterization; they are not assertions about the current manuscript's
general frame-factor labels or order of multiplication.

The same T^H Z T is the observable pullback for the separate active-state
law psi' = T psi.  Frame congruences therefore compose in the reverse
order to left state actions.  A moving frame uses Omega = T^-1 dT/ds.

Generic-axis identities are polynomial proofs modulo u.u = 1 and
c*c-h*h = 1.  Additional exact rational examples test composite and
varying transformations without choosing a physical preferred axis.
"""

import time

import sympy as s


started = time.perf_counter()
checks = 0
one = s.eye(2)
zero = s.zeros(2)
sigma = (
    s.Matrix(((0, 1), (1, 0))),
    s.Matrix(((0, -s.I), (s.I, 0))),
    s.diag(1, -1),
)


def vector(v):
    return sum((v[k] * sigma[k] for k in range(3)), zero)


def scalar(z):
    return s.trace(z) / 2


def coefficients(z):
    return s.Matrix((scalar(z), *(scalar(sigma[k] * z) for k in range(3))))


def adjugate(z):
    return s.trace(z) * one - z


def sigma_conjugate(z):
    return adjugate(z.H)


def norm_squared(z):
    return scalar(z.H * z)


def commutator(a, b):
    return a * b - b * a


def entries(expression):
    return list(expression) if isinstance(expression, s.MatrixBase) else [expression]


def check(name, expression, reduce=s.expand):
    global checks
    residuals = [reduce(value) for value in entries(expression)]
    assert all(value == 0 for value in residuals), (name, residuals)
    checks += 1
    print(f"PASS {name}")


def check_nonzero(name, expression):
    global checks
    assert any(s.cancel(value) != 0 for value in entries(expression)), name
    checks += 1
    print(f"PASS {name}")


# Independent generators and products.
ux, uy, uz, c, h = s.symbols("ux uy uz c h", real=True)
u = s.Matrix((ux, uy, uz))
a0, ax, ay, az = s.symbols("a0 ax ay az", real=True)
v = s.Matrix((ax, ay, az))
N = vector(u)
Z = a0 * one + vector(v)


def unit_reduce(expression):
    """Exact reduction by independent monic polynomial constraints."""
    numerator = s.together(expression).as_numer_denom()[0]
    numerator = s.rem(s.expand(numerator), c**2 - h**2 - 1, c)
    return s.rem(numerator, uz**2 + ux**2 + uy**2 - 1, uz).expand()


check("arbitrary-axis involution", N**2 - one, unit_reduce)
positive = (one + N) / 2
negative = (one - N) / 2
check("arbitrary-axis spectral projectors", positive**2 - positive, unit_reduce)
check("complementary spectral projectors", positive * negative, unit_reduce)
check("spectral resolution of the axis", positive - negative - N)

L = c * one - h * N
Li = c * one + h * N
gamma = c**2 + h**2
gamma_beta = 2 * c * h
check("boost inverse", L * Li - one, unit_reduce)
check("boost determinant", L.det() - 1, unit_reduce)
check("boost Hermitian conjugation", L.H - L)
check("boost norm is gamma", norm_squared(L) - gamma, unit_reduce)
check("Lorentz gamma identity", gamma**2 - gamma_beta**2 - 1, unit_reduce)
check("half-angle gamma coefficient", c**2 - (gamma + 1) / 2, unit_reduce)
check("half-angle vector coefficient", h * (gamma + 1) - c * gamma_beta, unit_reduce)
check("boost square has negative velocity", L**2 - gamma * one + gamma_beta * N, unit_reduce)
U = gamma * one + gamma_beta * N
check("proper velocity is inverse boost square", Li**2 - U, unit_reduce)
check("proper velocity determinant", U.det() - 1, unit_reduce)
check("boost carries its proper velocity to rest", L * U * L - one, unit_reduce)
check("proper velocity norm is not determinant", norm_squared(U) - (2 * gamma**2 - 1), unit_reduce)
check("boost is not generally unitary", L.H * L - gamma * one + gamma_beta * N, unit_reduce)

# Exponent-to-half-angle construction, including a possible real scalar.
theta, real_scalar = s.symbols("theta real_scalar", real=True)
spectral_exponential = s.exp(-theta / 2) * positive + s.exp(theta / 2) * negative
hyperbolic_exponential = s.cosh(theta / 2) * one - s.sinh(theta / 2) * N
check("real exponential gives hyperbolic half-angle form",
      spectral_exponential - hyperbolic_exponential, lambda q: s.simplify(q.rewrite(s.exp)))
with_scalar = s.exp(real_scalar) * L
check("real scalar exponent changes factor determinant",
      with_scalar.det() - s.exp(2 * real_scalar), unit_reduce)
check("real scalar exponent dilates spacetime",
      with_scalar.H * Z * with_scalar - s.exp(2 * real_scalar) * L.H * Z * L)
check("real scalar exponent scales Minkowski square by exp(4a)",
      (with_scalar.H * Z * with_scalar).det()
      - s.exp(4 * real_scalar) * Z.det(), unit_reduce)
check("scalar factor cancels from a field similarity",
      with_scalar * Z * (s.exp(-real_scalar) * Li) - L * Z * Li)
check("positive-norm normalization changes the boost determinant",
      (L / s.sqrt(gamma)).det() - 1 / gamma, unit_reduce)
check("positive-norm normalization changes the spacetime metric",
      ((L * Z * L) / gamma).det() - Z.det() / gamma**2, unit_reduce)

# Generic-axis congruence and similarity are different representations.
boosted = L.H * Z * L
expected = ((gamma * a0 - gamma_beta * u.dot(v)) * one
            + vector(v + (gamma - 1) * u.dot(v) * u - gamma_beta * a0 * u))
check("generic-axis time and position boost", boosted - expected, unit_reduce)
check("generic-axis real spacetime remains real", boosted.H - boosted)
check("generic-axis Minkowski determinant invariance", boosted.det() - Z.det(), unit_reduce)
check("generic-axis positive norm changes for a time increment",
      norm_squared(L * one * L) - (2 * gamma**2 - 1), unit_reduce)
check("similarity fixes scalar time instead of mixing it", L * one * Li - one, unit_reduce)
transverse = v - u.dot(v) * u
check("boost congruence leaves transverse vectors unchanged",
      L.H * vector(transverse) * L - vector(transverse), unit_reduce)
check("right multiplication by the boost square adds a transverse imaginary part",
      vector(transverse) * L**2 - gamma * vector(transverse)
      - s.I * gamma_beta * vector(u.cross(transverse)), unit_reduce)
parallel = a0 * one + u.dot(v) * N
check("boost factors combine for a commuting longitudinal paravector",
      L.H * parallel * L - parallel * L**2, unit_reduce)

ex, ey, ez, bx, by, bz = s.symbols("ex ey ez bx by bz", real=True)
electric = s.Matrix((ex, ey, ez))
magnetic = s.Matrix((bx, by, bz))
field = vector(electric + s.I * magnetic)
electric_prime = (gamma * electric - (gamma - 1) * u.dot(electric) * u
                  + gamma_beta * u.cross(magnetic))
magnetic_prime = (gamma * magnetic - (gamma - 1) * u.dot(magnetic) * u
                  - gamma_beta * u.cross(electric))
field_prime = L * field * Li
check("generic-axis electric and magnetic field boost",
      field_prime - vector(electric_prime + s.I * magnetic_prime), unit_reduce)
check("field similarity retains zero scalar part", scalar(field_prime), unit_reduce)
check("field E squared minus B squared invariance",
      electric_prime.dot(electric_prime) - magnetic_prime.dot(magnetic_prime)
      - electric.dot(electric) + magnetic.dot(magnetic), unit_reduce)
check("field E dot B invariance",
      electric_prime.dot(magnetic_prime) - electric.dot(magnetic), unit_reduce)
check("field determinant encodes both field invariants",
      field.det() + electric.dot(electric) - magnetic.dot(magnetic)
      + 2 * s.I * electric.dot(magnetic))

# The same arbitrary-axis spacetime calculation applies to momentum.
check("generic-axis energy and momentum mass shell",
      (gamma * a0 - gamma_beta * u.dot(v))**2
      - (v + (gamma - 1) * u.dot(v) * u - gamma_beta * a0 * u).dot(
          v + (gamma - 1) * u.dot(v) * u - gamma_beta * a0 * u)
      - a0**2 + v.dot(v), unit_reduce)
clock_step = one + (gamma_beta / gamma) * N
check("moving clock rate and position",
      L * clock_step * L - one / gamma, unit_reduce)
rod_step = gamma_beta / gamma * u.dot(v) * one + vector(v)
check("rod simultaneous in moving frame",
      scalar(L * rod_step * L), unit_reduce)
check("arbitrary-axis measured rod contraction",
      L * rod_step * L - vector(v + (1 / gamma - 1) * u.dot(v) * u), unit_reduce)
check("relativity of simultaneity", scalar(L * vector(v) * L)
      + gamma_beta * u.dot(v), unit_reduce)
velocity_denominator = gamma - gamma_beta * u.dot(v)
velocity_prime = (v + (gamma - 1) * u.dot(v) * u - gamma_beta * u) / velocity_denominator
check("velocity addition from the transformed tangent",
      L * (one + vector(v)) * L
      - velocity_denominator * (one + vector(velocity_prime)), unit_reduce)
check("velocity transformation preserves timelike and null cones",
      (1 - velocity_prime.dot(velocity_prime)) * velocity_denominator**2
      - 1 + v.dot(v), unit_reduce)

# Exact non-collinear examples exercise fully composite frame factors.
axes = (
    s.Matrix((s.Rational(2, 3), s.Rational(1, 3), s.Rational(2, 3))),
    s.Matrix((s.Rational(-2, 7), s.Rational(3, 7), s.Rational(6, 7))),
    s.Matrix((s.Rational(1, 3), s.Rational(2, 3), s.Rational(-2, 3))),
)


def boost(axis, c_value=s.Rational(5, 4), h_value=s.Rational(3, 4)):
    return c_value * one - h_value * vector(axis)


def rotation(axis):
    return s.Rational(4, 5) * one - s.I * s.Rational(3, 5) * vector(axis)


def four_vector_map(factor):
    basis = (one, *sigma)
    return s.Matrix.hstack(*(coefficients(factor.H * item * factor) for item in basis)).applyfunc(s.cancel)


sample_position = 7 * one + vector(s.Matrix((2, -3, 5)))
sample_field = vector(s.Matrix((2, -3, 5)) + s.I * s.Matrix((-1, 4, 2)))
sample_velocity = s.Rational(5, 4) * one + s.Rational(3, 4) * vector(axes[2])
sample_amplitude = ((one + vector(axes[0]))
                    * (2 * one + s.I * vector(axes[1]))).expand()
other_amplitude = (vector(axes[2]) * sample_amplitude).expand()
metric = s.diag(1, -1, -1, -1)

factors = [boost(axis) for axis in axes]
factors.append(rotation(axes[0]) * boost(axes[1]) * rotation(axes[2]))
factors.append(boost(axes[0]) * boost(axes[1], s.Rational(13, 12), s.Rational(5, 12)))
factors = [factor.applyfunc(s.cancel) for factor in factors]

for number, factor in enumerate(factors, 1):
    inverse = adjugate(factor)
    transformed_position = (factor.H * sample_position * factor).expand()
    transformed_field = (factor.H * sample_field * inverse.H).expand()
    transformed_velocity = (factor.H * sample_velocity * factor).expand()
    force = ((sample_field * sample_velocity + sample_velocity * sample_field.H) / 2).expand()
    transformed_force = (transformed_field * transformed_velocity
                         + transformed_velocity * transformed_field.H).expand() / 2
    transform = four_vector_map(factor)
    prefix = f"exact frame {number}"
    check(f"{prefix}: determinant one", factor.det() - 1)
    check(f"{prefix}: sigma conjugate equals inverse Hermitian", sigma_conjugate(factor) - inverse.H)
    check(f"{prefix}: real congruence", transformed_position.H - transformed_position)
    check(f"{prefix}: interval invariance", transformed_position.det() - sample_position.det())
    check(f"{prefix}: component metric invariance", transform.T * metric * transform - metric)
    check(f"{prefix}: coordinate orientation", transform.det() - 1)
    check(f"{prefix}: field invariants", transformed_field.det() - sample_field.det())
    check(f"{prefix}: Lorentz force covariance", transformed_force - factor.H * force * factor)
    check(f"{prefix}: force tangent to mass shell", scalar(adjugate(transformed_velocity) * transformed_force))
    stress = (sample_field * sample_velocity * sample_field.H).expand() / 2
    check(f"{prefix}: energy momentum flux covariance",
          transformed_field * transformed_velocity * transformed_field.H / 2 - factor.H * stress * factor)
    check(f"{prefix}: invariant energy position pairing",
          scalar(adjugate(transformed_velocity) * transformed_position)
          - scalar(adjugate(sample_velocity) * sample_position))
    check_nonzero(f"{prefix}: Euclidean coefficient norm is not preserved",
                  norm_squared(transformed_position) - norm_squared(sample_position))
    amplitude_prime = (factor * sample_amplitude).expand()
    other_prime = (factor * other_amplitude).expand()
    check(f"{prefix}: active-state matrix-element pullback",
          scalar(other_prime.H * sample_position * amplitude_prime)
          - scalar(other_amplitude.H * transformed_position * sample_amplitude))
    check(f"{prefix}: normalized active-state expectation",
          scalar(amplitude_prime.H * sample_position * amplitude_prime) / norm_squared(amplitude_prime)
          - scalar(sample_amplitude.H * transformed_position * sample_amplitude)
          / scalar(sample_amplitude.H * factor.H * factor * sample_amplitude), s.cancel)
    check(f"{prefix}: state bilinear uses the opposite congruence order",
          amplitude_prime * amplitude_prime.H
          - factor * sample_amplitude * sample_amplitude.H * factor.H)

# Alternative actions fail concretely for the required physical objects.
sample_boost = factors[0]
check_nonzero("similarity fails to preserve real spacetime",
              (sample_boost * sample_position * sample_boost.inv()).H
              - sample_boost * sample_position * sample_boost.inv())
check_nonzero("left multiplication fails to preserve real spacetime",
              (sample_boost * sample_position).H - sample_boost * sample_position)
check_nonzero("congruence wrongly creates a field scalar part",
              scalar(sample_boost * sample_field * sample_boost))
check_nonzero("field positive norm changes under a boost",
              norm_squared(sample_boost * sample_field * sample_boost.inv()) - norm_squared(sample_field))
sample_transverse = vector(axes[0].cross(axes[1]))
check_nonzero("boost congruence cannot generally combine to right multiplication",
              sample_boost.H * sample_transverse * sample_boost
              - sample_transverse * sample_boost**2)


def rotation_reduce(expression):
    numerator = s.together(expression).as_numer_denom()[0]
    numerator = s.rem(s.expand(numerator), c**2 + h**2 - 1, c)
    return s.rem(numerator, uz**2 + ux**2 + uy**2 - 1, uz).expand()


generic_rotation = c * one - s.I * h * N
generic_rodrigues = ((c**2 - h**2) * v + 2 * h**2 * u.dot(v) * u
                    - 2 * c * h * u.cross(v))
check("arbitrary-axis rotation pullback has the negative cross-product sign",
      generic_rotation.H * Z * generic_rotation - a0 * one - vector(generic_rodrigues),
      rotation_reduce)

sample_rotation = rotation(axes[1])
check("rotation unitary", sample_rotation.H * sample_rotation - one)
check("rotation determinant one", sample_rotation.det() - 1)
check("rotation norm invariance", norm_squared(sample_rotation.H * sample_field * sample_rotation)
      - norm_squared(sample_field))
cos_angle, sin_angle = s.Rational(7, 25), s.Rational(24, 25)
sample_spatial = s.Matrix((2, -3, 5))
rodrigues = (cos_angle * sample_spatial + (1 - cos_angle) * axes[1].dot(sample_spatial) * axes[1]
             - sin_angle * axes[1].cross(sample_spatial))
check("preferred rotation congruence gives passive Rodrigues rotation",
      sample_rotation.H * vector(sample_spatial) * sample_rotation - vector(rodrigues))
check("active state rotation co-transforms its observable in the opposite order",
      scalar((sample_rotation * sample_amplitude).H
             * (sample_rotation * sample_position * sample_rotation.H)
             * (sample_rotation * sample_amplitude))
      - scalar(sample_amplitude.H * sample_position * sample_amplitude))

# Passive frame congruences pull back in reverse order to left state actions.
first_frame = rotation(axes[0])
second_frame = boost(axes[1])
frame_product = (first_frame * second_frame).expand()
sequential_position = (second_frame.H * first_frame.H * sample_position
                       * first_frame * second_frame).expand()
check("two frame congruences use T1 T2 for first T1 then T2",
      sequential_position - frame_product.H * sample_position * frame_product)
check("component maps have the same reversed factor composition",
      four_vector_map(frame_product) - four_vector_map(second_frame) * four_vector_map(first_frame))
check("field similarities use the same frame factor order",
      second_frame.H * (first_frame.H * sample_field * adjugate(first_frame).H)
      * adjugate(second_frame).H
      - frame_product.H * sample_field * adjugate(frame_product).H)
check_nonzero("reversing noncommuting factors changes the spacetime map",
              sequential_position
              - (second_frame * first_frame).H * sample_position * (second_frame * first_frame))
check("successive left state actions instead use T2 T1",
      second_frame * (first_frame * sample_amplitude)
      - (second_frame * first_frame) * sample_amplitude)
# Collinear composition, and a non-collinear counterexample to pure boosts.
first = boost(axes[0])
second = boost(axes[0], s.Rational(13, 12), s.Rational(5, 12))
composed = first * second
c_total = s.Rational(5, 4) * s.Rational(13, 12) + s.Rational(3, 4) * s.Rational(5, 12)
h_total = s.Rational(3, 4) * s.Rational(13, 12) + s.Rational(5, 4) * s.Rational(5, 12)
check("collinear rapidity addition", composed - c_total * one + h_total * vector(axes[0]))
beta_first, beta_second = s.Rational(15, 17), s.Rational(65, 97)
beta_total = 2 * c_total * h_total / (c_total**2 + h_total**2)
check("collinear velocity composition",
      beta_total - (beta_first + beta_second) / (1 + beta_first * beta_second))
noncollinear = boost(axes[0]) * boost(axes[1])
check_nonzero("noncollinear boost composition need not be Hermitian", noncollinear.H - noncollinear)

# Explicit coordinate chain rules and Maxwell covariance on polynomial fields.
t, x, y, z = s.symbols("t x y z", real=True)
coordinates = (t, x, y, z)
factor = factors[-1]
inverse = adjugate(factor)
transform = four_vector_map(factor)
old_coordinates = transform.inv() * s.Matrix(coordinates)
substitution = dict(zip(coordinates, old_coordinates))


def evaluate_old(expression):
    return expression.subs(substitution, simultaneous=True).expand()


def differential(expression):
    return expression.diff(t) + sum((sigma[k] * expression.diff(coordinates[k + 1])
                                     for k in range(3)), s.zeros(*expression.shape))


def differential_star(expression):
    return expression.diff(t) - sum((sigma[k] * expression.diff(coordinates[k + 1])
                                     for k in range(3)), s.zeros(*expression.shape))


def dalembertian(expression):
    return expression.diff(t, 2) - sum((expression.diff(coordinates[k + 1], 2)
                                      for k in range(3)), s.zeros(*expression.shape))


scalar_test = (t**2 + t * x + 2 * y * z + x**2) * one
check("coordinate chain rule for the paravector derivative",
      differential(evaluate_old(scalar_test))
      - inverse * evaluate_old(differential(scalar_test)) * inverse.H)
potential = ((x**3 + t * y * z) * one
             + vector(s.Matrix((t * x**2 + y * z, t**2 * y + x * z, t * z**2 + x * y))))
potential_product = differential(potential).expand()
potential_gauge = scalar(potential_product)
polynomial_field = sigma_conjugate(potential_product - potential_gauge * one).expand()
current = sigma_conjugate(differential(polynomial_field)).expand()
check("real potential gives a real Maxwell source", current.H - current)
new_potential = (factor.H * evaluate_old(potential) * factor).expand()
new_potential_product = differential(new_potential).expand()
check("Lorenz gauge scalar is invariant",
      scalar(new_potential_product) - evaluate_old(potential_gauge))
new_field = sigma_conjugate(new_potential_product - scalar(new_potential_product) * one).expand()
check("field similarity follows from potential differentiation",
      new_field - factor.H * evaluate_old(polynomial_field) * inverse.H)
new_current = (factor.H * evaluate_old(current) * factor).expand()
check("Maxwell covariance at transformed spacetime arguments",
      differential(new_field) - sigma_conjugate(new_current))

# Scalar and chiral wave equations transform at corresponding events.
mass = s.symbols("mass", real=True)
scalar_wave = (t**3 * x + t**2 * y**2 + x**2 * y * z + z**4) * one
check("Klein-Gordon covariance at transformed spacetime arguments",
      dalembertian(evaluate_old(scalar_wave)) + mass**2 * evaluate_old(scalar_wave)
      - evaluate_old(dalembertian(scalar_wave) + mass**2 * scalar_wave))
first_amplitude = ((t**2 + x * y) * sample_amplitude + (t * z + y) * other_amplitude).expand()
second_amplitude = ((t * x + z**2) * sample_amplitude + (t**3 - y * z) * other_amplitude).expand()
first_prime = (inverse * evaluate_old(first_amplitude)).expand()
second_prime = (factor.H * evaluate_old(second_amplitude)).expand()
first_residual = s.I * differential(second_amplitude) - mass * first_amplitude
second_residual = s.I * differential_star(first_amplitude) - mass * second_amplitude
check("first free Dirac equation uses Psi1'=T^-1 Psi1",
      s.I * differential(second_prime) - mass * first_prime
      - inverse * evaluate_old(first_residual))
check("second free Dirac equation uses Psi2'=T^H Psi2",
      s.I * differential_star(first_prime) - mass * second_prime
      - factor.H * evaluate_old(second_residual))
check("first potential term has the same Dirac covariance",
      sigma_conjugate(new_potential) * second_prime
      - inverse * evaluate_old(sigma_conjugate(potential) * second_amplitude))
check("second potential term has the same Dirac covariance",
      new_potential * first_prime - factor.H * evaluate_old(potential * first_amplitude))

# Constant transformations commute with world-line differentiation/integration.
parameter = s.symbols("parameter", real=True)
worldline = parameter * sample_velocity + parameter**2 * vector(axes[1]) / 7
worldline_prime = factor.H * worldline * factor
check("constant frame and world-line tangent commute",
      worldline_prime.diff(parameter) - factor.H * worldline.diff(parameter) * factor)
check("constant frame and world-line integral commute",
      worldline_prime - s.integrate(factor.H * worldline.diff(parameter) * factor, (parameter, 0, parameter)))

# Varying determinant-one frames, with no assumption that generators commute.
nilpotent_first = s.Matrix(((0, 1), (0, 0)))
nilpotent_second = s.Matrix(((0, 0), (1, 0)))
varying = (one + parameter * nilpotent_first) * (one + parameter**2 * nilpotent_second)
varying_inverse = adjugate(varying)
omega = varying_inverse * varying.diff(parameter)
varying_field = sample_field + parameter * vector(axes[0]) + parameter**2 * s.I * vector(axes[1])
varying_real = sample_position + parameter * sample_velocity + parameter**2 * vector(axes[2])
similarity = varying.H * varying_field * varying_inverse.H
congruence = varying.H * varying_real * varying
check("varying frame determinant one", varying.det() - 1)
check("varying-frame generator scalar part vanishes", scalar(omega))
check("varying inverse derivative", varying_inverse.diff(parameter) + omega * varying_inverse)
check("varying similarity first derivative",
      similarity.diff(parameter) - varying.H * varying_field.diff(parameter) * varying_inverse.H
      - commutator(omega.H, similarity))
check("varying similarity second derivative",
      similarity.diff(parameter, 2) - varying.H * varying_field.diff(parameter, 2) * varying_inverse.H
      - 2 * commutator(omega.H, similarity.diff(parameter))
      - commutator(omega.H.diff(parameter), similarity)
      + commutator(omega.H, commutator(omega.H, similarity)))
check("varying congruence first derivative",
      congruence.diff(parameter) - varying.H * varying_real.diff(parameter) * varying
      - omega.H * congruence - congruence * omega)
check("varying-frame basis derivative",
      (varying.H * vector(axes[2]) * varying_inverse.H).diff(parameter)
      - commutator(omega.H, varying.H * vector(axes[2]) * varying_inverse.H))
check_nonzero("local physical acceleration excludes derivatives of the frame",
              congruence.diff(parameter) - varying.H * varying_real.diff(parameter) * varying)
varying_rotation = s.cos(parameter / 2) * one - s.I * s.sin(parameter / 2) * vector(axes[0])
rotation_omega = varying_rotation.H * varying_rotation.diff(parameter)
check("varying rotation generator is anti-Hermitian",
      rotation_omega.H + rotation_omega, s.trigsimp)
check("rotation congruence frame terms equal its similarity commutator",
      rotation_omega.H * sample_position + sample_position * rotation_omega
      - commutator(rotation_omega.H, sample_position), s.trigsimp)

print(f"\n{checks} exact checks passed in {time.perf_counter() - started:.2f} s.")
