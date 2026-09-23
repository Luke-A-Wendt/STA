"""Check the energy-plus-momentum W convention against component Dirac physics.

Run ``python3 verify_dirac.py`` (requires SymPy). Potentials and amplitudes are
arbitrary smooth functions of all four coordinates, so the differential checks
retain their derivatives and do not assume constant fields or commuting energy
and momentum. W(Z) has adj(Z) above Z; gamma^j = -W(sigma_j).
"""

import time

import sympy as s


started = time.perf_counter()
checks = 0
i = s.I
one, zero = s.eye(2), s.zeros(2)
identity = s.eye(4)
sigma = (s.Matrix(((0, 1), (1, 0))),
         s.Matrix(((0, -i), (i, 0))), s.diag(1, -1))
t, x, y, z, q = s.symbols('t x y z q', real=True)
m = s.symbols('m', positive=True)
coordinates = (t, x, y, z)
spatial = coordinates[1:]
V = s.Function('V', real=True)(*coordinates)
A = s.Matrix([s.Function(f'A{k}', real=True)(*coordinates) for k in range(3)])
psi = s.Matrix([s.Function(f'psi{k}')(*coordinates) for k in range(4)])


def check(name, residual):
    global checks
    entries = list(residual) if isinstance(residual, s.MatrixBase) else [residual]
    values = [s.expand(entry) for entry in entries]
    assert all(entry == 0 for entry in values), (name, values)
    checks += 1
    print(f'PASS {name}', flush=True)


def vector(v):
    return sum((v[k] * sigma[k] for k in range(3)), zero)


def adj(Z):
    return s.trace(Z) * one - Z


def W(Z):
    return s.BlockMatrix(((zero, adj(Z)), (Z, zero))).as_explicit()


def energy(field, potential=V):
    return i * field.diff(t) - q * potential * field


def momentum(k, field, potential=A):
    return -i * field.diff(spatial[k]) - q * potential[k] * field


def sigma_momentum(field):
    return sum((sigma[k] * momentum(k, field) for k in range(3)), s.zeros(*field.shape))


def momentum_squared(field):
    return sum((momentum(k, momentum(k, field)) for k in range(3)), s.zeros(*field.shape))


def W_momentum(field, potential=A):
    return sum((W(sigma[k]) * momentum(k, field, potential) for k in range(3)),
               s.zeros(*field.shape))


def dirac(field, mass=m, potential=V, vector_potential=A):
    return W(one) * energy(field, potential) + W_momentum(field, vector_potential) - mass * field


# Algebraic and differential operations are separate: det identities apply to
# commuting paravector coefficients; field commutators enter the operator square.
scalar = s.symbols('S')
v = s.Matrix(s.symbols('v0:3'))
Z = scalar * one + vector(v)
check('new algebraic W squares to the paravector determinant', W(Z)**2 - Z.det() * identity)
check('new W on adj Z equals the previous W on Z',
      W(adj(Z)) - s.BlockMatrix(((zero, Z), (adj(Z), zero))).as_explicit())

# Independent conventional Weyl gamma matrices, rather than deriving gamma from W.
gamma = (s.BlockMatrix(((zero, one), (one, zero))).as_explicit(),
         *(s.BlockMatrix(((zero, item), (-item, zero))).as_explicit() for item in sigma))
check('standard Weyl gamma map has gamma^j = -W(sigma_j)',
      s.Matrix.vstack(W(one) - gamma[0], *(-W(sigma[k]) - gamma[k+1] for k in range(3))))
conventional = gamma[0] * energy(psi) - sum((gamma[k+1] * momentum(k, psi)
                                                          for k in range(3)), s.zeros(4, 1)) - m * psi
check('all four Dirac components match conventional minimal coupling', dirac(psi) - conventional)

upper, lower = psi[:2, :], psi[2:, :]
check('explicit E-minus-p and E-plus-p block signs', dirac(psi) - s.Matrix.vstack(
    -m * upper + energy(lower) - sigma_momentum(lower),
    energy(upper) + sigma_momentum(upper) - m * lower))
Phi = V * one + vector(A)
partial_star = lower.diff(t) - sum((sigma[k] * lower.diff(spatial[k]) for k in range(3)), s.zeros(2, 1))
check('boxed energy-plus-momentum equals i partial-star minus q Phi',
      energy(lower) + sigma_momentum(lower) - i * partial_star + q * Phi * lower)

H_psi = q * V * psi + W(one) * (m * psi - W_momentum(psi))
standard_H_psi = q * V * psi + m * gamma[0] * psi + sum((gamma[0] * gamma[k+1] * momentum(k, psi)
                                                                      for k in range(3)), s.zeros(4, 1))
check('Hamiltonian contains minus W(p.sigma)', H_psi - standard_H_psi)
check('Dirac and Hamiltonian equations are equivalent',
      W(one) * dirac(psi) - i * psi.diff(t) + H_psi)

E = -A.diff(t) - s.Matrix([s.diff(V, coordinate) for coordinate in spatial])
B = s.Matrix([s.diff(A[2], y) - s.diff(A[1], z),
              s.diff(A[0], z) - s.diff(A[2], x),
              s.diff(A[1], x) - s.diff(A[0], y)])
F, Fstar = vector(E + i * B), vector(-E + i * B)
commutator = sigma_momentum(energy(lower)) - energy(sigma_momentum(lower))
check('electric commutator for arbitrary time-dependent potentials',
      commutator + i * q * vector(E) * lower)
check('Pauli momentum square for arbitrary vector potentials',
      sigma_momentum(sigma_momentum(lower)) - momentum_squared(lower) + q * vector(B) * lower)
check('F identity includes the electric commutator',
      momentum_squared(lower) + i * q * F * lower
      - sigma_momentum(sigma_momentum(lower)) + commutator)

squared = energy(energy(psi)) - momentum_squared(psi) - m**2 * psi - i * q * s.diag(Fstar, F) * psi
check('opposite-mass square has F-star above F', dirac(dirac(psi), -m) - squared)
check('opposite-mass operators commute', dirac(dirac(psi), -m) - dirac(dirac(psi, -m)))
box = psi.diff(t, 2) - sum((psi.diff(coordinate, 2) for coordinate in spatial), s.zeros(4, 1))
check('zero charge recovers minus Klein-Gordon', squared.subs(q, 0) + box + m**2 * psi)

# The same operator on the envelope measures the residual mechanical energy.
phase = s.exp(-i * m * t)
check('exact rest-phase energy shift', energy(phase * psi) / phase - m * psi - energy(psi))
check('exact rest-phase energy-square shift',
      energy(energy(phase * psi)) / phase - m**2 * psi - 2 * m * energy(psi) - energy(energy(psi)))
lower_shifted = (energy(energy(phase * lower)) - momentum_squared(phase * lower)
                 - m**2 * phase * lower - i * q * F * phase * lower) / phase
expected_shift = 2 * m * energy(lower) + energy(energy(lower)) - momentum_squared(lower) - i * q * F * lower
check('lower squared block after removing the rest phase', lower_shifted - expected_shift)
recursion = energy(lower) - (sigma_momentum(sigma_momentum(lower))
                            - energy(energy(lower)) - commutator) / (2 * m)
check('exact residual-energy recursion retains its electric sign', 2 * m * recursion - expected_shift)

# C takes (Psi1, Psi2) to (phi1, phi2) with phi2 = (Psi2-Psi1)/sqrt(2).
C = s.BlockMatrix(((one, one), (-one, one))).as_explicit() / s.sqrt(2)
phi = s.Matrix([s.Function(f'phi{k}')(*coordinates) for k in range(4)])
phi1, phi2 = phi[:2, :], phi[2:, :]
envelope_dirac = dirac(phase * (C.T * phi)) / phase
coupled = C * W(one) * envelope_dirac
expected_coupled = s.Matrix.vstack(energy(phi1) - sigma_momentum(phi2),
                                  2 * m * phi2 + energy(phi2) - sigma_momentum(phi1))
check('sum-and-difference amplitudes give the exact two envelope equations', coupled - expected_coupled)
leading_small = sigma_momentum(phi1) / (2 * m)
pauli_residual = i * phi1.diff(t) - (momentum_squared(phi1) / (2 * m) + q * V * phi1
                                    - q * vector(B) * phi1 / (2 * m))
check('first small-component iterate produces the conventional Pauli sign',
      energy(phi1) - sigma_momentum(leading_small) - pauli_residual)

# The manuscript uses Phi' = Phi + partial-star lambda, Psi' = exp(-iq lambda) Psi.
lam = s.Function('lambda', real=True)(*coordinates)
gauge_phase = s.exp(-i * q * lam)
gauge_V = V + lam.diff(t)
gauge_A = A - s.Matrix([lam.diff(coordinate) for coordinate in spatial])
check('Dirac gauge covariance including all potential derivatives',
      dirac(gauge_phase * psi, potential=gauge_V, vector_potential=gauge_A) / gauge_phase - dirac(psi))

print(f'\n{checks} checks passed in {time.perf_counter() - started:.2f} s.')
