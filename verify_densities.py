"""Independent verification of the local products in sta_notes.tex.

Uses an independent Pauli-matrix realization, matching the matrix appendix.
Run with ``python3 verify_densities.py`` (requires SymPy and NumPy).
"""
import time

import numpy as np
import sympy as s


started = time.monotonic()
I = s.eye(2)
sigma = (
    s.Matrix(((0, 1), (1, 0))),
    s.Matrix(((0, -s.I), (s.I, 0))),
    s.Matrix(((1, 0), (0, -1))),
)
count = 0


def vector(v):
    return sum((v[j] * sigma[j] for j in range(3)), s.zeros(2))


def sc(z):
    return s.trace(z) / 2


def density(psi, z=I):
    return sc(psi.H * z * psi)


def check(name, expr):
    global count
    entries = list(expr) if isinstance(expr, s.MatrixBase) else (expr,)
    assert all(s.expand(x) == 0 for x in entries), name
    count += 1


def require(name, condition):
    global count
    assert bool(condition), name
    count += 1


# Completely general complex paravector amplitude (eight real components).
a = s.symbols("a0:4", real=True)
b = s.symbols("b0:4", real=True)
psi = (a[0] + s.I * b[0]) * I + vector(
    [a[j + 1] + s.I * b[j + 1] for j in range(3)]
)
rho = density(psi)
check("ordinary norm equals component norm", rho - sum(x*x for x in a+b))

S, vx, vy, vz = s.symbols("S vx vy vz", real=True)
z = S * I + vector((vx, vy, vz))
q = density(psi, z)
check("real multiplicative weight gives real density", q - s.conjugate(q))

# A non-coordinate axis keeps all vector components involved.
u = s.Matrix((s.Rational(2, 3), s.Rational(-1, 3), s.Rational(2, 3)))
N = vector(u)
plus, minus = (I + N) / 2, (I - N) / 2
v = s.symbols("v", nonnegative=True)
spectral_z = S * I + v * N
check("general spectral density", density(psi, spectral_z)
      - (S+v)*density(plus*psi) - (S-v)*density(minus*psi))
check("scalar weight including zero", density(psi, S*I)-S*rho)
check("orthogonal sectors exhaust norm", density(plus*psi)+density(minus*psi)-rho)
check("positive sector projected density", density(psi, plus)-density(plus*psi))
check("negative sector projected density", density(psi, minus)-density(minus*psi))
check("singular nonnegative weight has nonzero null amplitudes", density(minus, I+N))
require("null amplitude really nonzero", density(minus) > 0)
check("zero weight", density(psi, s.zeros(2)))
check("boundary nonnegative density", density(psi, v*(I+N))-2*v*density(plus*psi))
check("negative-eigenvalue witness", density(minus, spectral_z)-(S-v)/2)
check("spin signed expectation positive sector", density(plus, N/2)/density(plus)-s.Rational(1,2))
check("spin signed expectation negative sector", density(minus, N/2)/density(minus)+s.Rational(1,2))

# Verify the chosen nonnegative square root, including both singular cases.
for root_plus, root_minus in ((s.Integer(2),s.Integer(1)),
                              (s.Integer(2),s.Integer(0)),
                              (s.Integer(0),s.Integer(2)),
                              (s.Integer(0),s.Integer(0))):
    root = root_plus*plus + root_minus*minus
    weight = root_plus**2*plus + root_minus**2*minus
    check("spectral square-root square", root*root-weight)
    check("spectral square-root Hermitian", root.H-root)
    check("weighted density equals root norm", density(psi, weight)-density(root*psi))

# General complex transform and general complex observable, without assuming
# determinant one, unitarity, Hermiticity or a particular amplitude sector.
c = s.symbols("c0:4", real=True)
d = s.symbols("d0:4", real=True)
T = s.Matrix(((c[0]+s.I*d[0], c[1]+s.I*d[1]),
              (c[2]+s.I*d[2], c[3]+s.I*d[3])))
zz = z + s.I*(2*I+vector((1,-2,3)))
check("generic complex local pullback", density(T*psi, zz)-density(psi,T.H*zz*T))
check("generic transformed norm weight", density(T*psi)-density(psi,T.H*T))

# Exact invertible examples test normalized ratios and nonnegative transformed
# weights. Rotations and boosts use the signs of the manuscript.
R = s.Rational(4,5)*I - s.I*s.Rational(3,5)*N
L = s.Rational(5,4)*I - s.Rational(3,4)*N
generic_T = s.Matrix(((1+s.I, 2-s.I),(-1+2*s.I,3)))
amplitudes = (s.Matrix(((1+s.I,2-s.I),(-2,3+2*s.I))),
              plus*(I+2*s.I*vector((1,2,-1))))
weights = (I, 2*I+N, plus, minus, s.zeros(2))
observables = weights + (N/2, -I, zz.subs({S:3,vx:1,vy:2,vz:-1}))
check("rotation unitary", R.H*R-I)
check("boost Hermitian", L.H-L)
check("boost determinant one", L.det()-1)
check("rotation preserves arbitrary amplitude density", density(R*psi)-rho)
check("direct boost density uses squared boost weight", density(L*psi)-density(psi,L*L))
require("direct boost can change local density", density(L*plus) != density(plus))

# The generic identity is psi'=A psi with pullback A^H Z A.  For the
# paper's active rotation and second Dirac entry, A is the direct frame factor.
for transform in (R, L, L*R, generic_T):
    require("transform invertible", transform.det() != 0)
    for weight in weights:
        transformed_weight = s.simplify(transform.H*weight*transform)
        check("congruence preserves Hermiticity", transformed_weight.H-transformed_weight)
        transformed_scalar = s.simplify(sc(transformed_weight))
        transformed_vector_square = sum(
            s.simplify(sc(transformed_weight*x))**2 for x in sigma
        )
        require("congruence preserves nonnegative scalar", transformed_scalar >= 0)
        require("congruence preserves S squared >= V squared",
                s.simplify(transformed_scalar**2-transformed_vector_square) >= 0)
    for amplitude in amplitudes:
        original_norm = s.simplify(density(amplitude))
        transformed_amplitude = transform*amplitude
        transformed_norm = s.simplify(density(transformed_amplitude))
        require("nonzero amplitude has positive transformed norm", transformed_norm > 0)
        for observable in observables:
            pullback = transform.H*observable*transform
            original_pullback_mean = density(amplitude,pullback)/original_norm
            transformed_mean = density(transformed_amplitude,observable)/transformed_norm
            check("normalized local expectation pullback", transformed_mean
                  - original_norm/transformed_norm*original_pullback_mean)

print(f"{count} exact algebra checks passed in {time.monotonic()-started:.2f}s")

# Reproducible numerical stress checks vary axes, singularity, amplitudes,
# transform determinants and complex signed observables.
rng = np.random.default_rng(77391)
basis = np.asarray([np.asarray(x,dtype=complex) for x in sigma])
eye = np.eye(2,dtype=complex)
numeric_count = 0


def num_density(amplitude, weight=eye):
    return np.trace(amplitude.conj().T @ weight @ amplitude)/2


def near(left, right):
    global numeric_count
    assert np.allclose(left,right,rtol=2e-11,atol=2e-11), (left,right)
    numeric_count += 1


for case in range(200):
    amplitude = rng.normal(size=(2,2)) + 1j*rng.normal(size=(2,2))
    transform = rng.normal(size=(2,2)) + 1j*rng.normal(size=(2,2))
    assert abs(np.linalg.det(transform)) > 1e-10
    direction = rng.normal(size=3)
    direction /= np.linalg.norm(direction)
    axis = np.einsum("j,jab->ab", direction, basis)
    projection_plus, projection_minus = (eye+axis)/2, (eye-axis)/2
    length = float(rng.uniform(0,3))
    scalar = length if case%2 == 0 else length + float(rng.uniform(.1,3))
    weight = scalar*eye+length*axis
    root = np.sqrt(scalar+length)*projection_plus + np.sqrt(scalar-length)*projection_minus
    rho = num_density(amplitude)
    q = num_density(amplitude,weight)
    near(q,(scalar+length)*num_density(projection_plus@amplitude)
         +(scalar-length)*num_density(projection_minus@amplitude))
    near(q,num_density(root@amplitude))
    assert q.real >= -2e-11 and abs(q.imag)<2e-11
    transformed_amplitude = transform@amplitude
    transformed_weight = transform.conj().T@weight@transform
    near(num_density(transformed_amplitude,weight),num_density(amplitude,transformed_weight))
    assert np.linalg.eigvalsh(transformed_weight).min() >= -2e-11
    observable = rng.normal(size=(2,2)) + 1j*rng.normal(size=(2,2))
    rhop = num_density(transformed_amplitude)
    near(num_density(transformed_amplitude,observable)/rhop,
         (rho/rhop)*(num_density(amplitude,transform.conj().T@observable@transform)/rho))
    probability_plus = num_density(projection_plus@amplitude)/rho
    probability_minus = num_density(projection_minus@amplitude)/rho
    near(probability_plus+probability_minus,1)
    assert -.00000000001 <= probability_plus.real <= 1.00000000001

print(f"{numeric_count} numerical identities passed across 200 random cases")
print("All local-density, positivity, singular-weight and normalization checks passed.")
