"""Independent exact checks for Appendix Conventions.

Run ``python3 verify_conventions.py`` (requires SymPy).  The Pauli matrices
faithfully realize the paravector algebra.  Polynomial identities use arbitrary
axes; rational noncollinear examples check nonunitary/composite actions without
accidental commutation.  No assertions are extracted from manuscript text.
"""

import time

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


def inner(a, b):
    return sc(a.H * b)


def norm2(z):
    return inner(z, z)


def expectation(z, psi):
    return sc(psi.H * z * psi) / norm2(psi)


def minkowski(x, y):
    return sc(adj(x) * y)


def W(z):
    return s.BlockMatrix(((zero, adj(z)), (z, zero))).as_explicit()


def check(name, residual, reduce=s.expand):
    global checks
    entries = list(residual) if isinstance(residual, s.MatrixBase) else [residual]
    values = [reduce(entry) for entry in entries]
    assert all(value == 0 for value in values), (name, values)
    checks += 1
    print(f"PASS {name}")


def counterexample(name, residual):
    global checks
    entries = list(residual) if isinstance(residual, s.MatrixBase) else [residual]
    assert any(s.cancel(entry) != 0 for entry in entries), name
    checks += 1
    print(f"PASS {name}")


ux, uy, uz, c, h = s.symbols("ux uy uz c h", real=True)
u = s.Matrix((ux, uy, uz))
t, x, y, z = s.symbols("t x y z", real=True)
r = s.Matrix((x, y, z))
N = vector(u)
X = t * one + vector(r)


def constrained(value, circular=False):
    numerator = s.together(value).as_numer_denom()[0]
    relation = c**2 + h**2 - 1 if circular else c**2 - h**2 - 1
    remainder = s.rem(s.expand(numerator), relation, c)
    return s.rem(remainder, ux**2 + uy**2 + uz**2 - 1, uz).expand()


def circular(value):
    return constrained(value, circular=True)


# Hamilton's right-handed multiplication fixes the sign of its sigma image.
for sign in (-1, 1):
    q = tuple(sign * s.I * item for item in sigma)
    for j in range(3):
        check(f"quaternion sign {sign:+}: unit {j + 1} squares to -1", q[j]**2 + one)
        k, ell = (j + 1) % 3, (j + 2) % 3
        check(f"quaternion sign {sign:+}: oriented product {j + 1},{k + 1}",
              q[j] * q[k] + sign * q[ell])

Q = c * one - s.I * h * N
for sign in (-1, 1):
    R = c * one + sign * s.I * h * N
    check(f"rotation exponent {sign:+}: unitary", R.H * R - one, circular)
    check(f"rotation exponent {sign:+}: determinant one", R.det() - 1, circular)
    check(f"rotation exponent {sign:+}: adjoint reverses exponent",
          R.H - (c * one - sign * s.I * h * N))
    for adjoint_left in (False, True):
        A = R.H if adjoint_left else R
        angle_sign = sign if adjoint_left else -sign
        rotated = ((c**2 - h**2) * r + 2 * h**2 * u.dot(r) * u
                   + angle_sign * 2 * c * h * u.cross(r))
        side = "left" if adjoint_left else "right"
        check(f"rotation exponent {sign:+}, H on {side}: Rodrigues sign {angle_sign:+}",
              A * X * A.H - t * one - vector(rotated), circular)
        check(f"rotation exponent {sign:+}, H on {side}: quaternion conversion",
              A * (-s.I * vector(r)) * A.H + s.I * vector(rotated), circular)
check("Hamilton Q is R with negative exponent", Q - (c * one - s.I * h * N))
check("positive-exponent rotor is Hamilton Q adjoint", c * one + s.I * h * N - Q.H)

# Hermitian boost conjugation cannot change the velocity sign by changing order.
gamma, gamma_beta = c**2 + h**2, 2 * c * h
for sign in (-1, 1):
    L = c * one + sign * h * N
    Li = c * one - sign * h * N
    tp = gamma * t + sign * gamma_beta * u.dot(r)
    rp = r + (gamma - 1) * u.dot(r) * u + sign * gamma_beta * t * u
    transformed = tp * one + vector(rp)
    check(f"boost exponent {sign:+}: Hermitian", L.H - L)
    check(f"boost exponent {sign:+}: inverse reverses exponent", L * Li - one, constrained)
    check(f"boost exponent {sign:+}: determinant one", L.det() - 1, constrained)
    check(f"boost exponent {sign:+}: either H placement gives same event",
          L.H * X * L - L * X * L.H)
    check(f"boost exponent {sign:+}: coordinate velocity signs", L * X * L.H - transformed,
          constrained)
    check(f"boost exponent {sign:+}: Minkowski determinant", transformed.det() - X.det(),
          constrained)
    check(f"boost exponent {sign:+}: positive norm is gamma", norm2(L) - gamma, constrained)

# Exact rational rotations/boosts with distinct non-coordinate axes.
axis_r = s.Matrix((s.Rational(2, 3), s.Rational(1, 3), s.Rational(2, 3)))
axis_l = s.Matrix((s.Rational(-2, 7), s.Rational(3, 7), s.Rational(6, 7)))
Rminus = s.Rational(4, 5) * one - s.I * s.Rational(3, 5) * vector(axis_r)
Rplus = Rminus.H
Lminus = s.Rational(5, 4) * one - s.Rational(3, 4) * vector(axis_l)
Lplus = adj(Lminus)
Z = (2 + s.I) * one + vector(s.Matrix((1 + 2*s.I, -3 + s.I, 4 - s.I)))
Y = (1 - 2*s.I) * one + vector(s.Matrix((2 - s.I, 1 + s.I, -1 + 3*s.I)))
psi = s.Matrix(((1 + s.I, 0), (2 - s.I, 0)))
chi = s.Matrix(((2 - 3*s.I, 0), (-1 + 2*s.I, 0)))
event = 7 * one + vector(s.Matrix((2, -3, 5)))
event2 = 3 * one + vector(s.Matrix((-1, 4, 2)))
field = vector(s.Matrix((2 - s.I, -3 + 4*s.I, 5 + 2*s.I)))

counterexample("test rotation and boost do not commute", Rminus * Lminus - Lminus * Rminus)
counterexample("rotation adjoint is not itself", Rminus.H - Rminus)
counterexample("boost adjoint is not inverse", Lminus.H - adj(Lminus))

for name, R, A in (("negative/direct", Rminus, Rminus),
                   ("positive/adjoint", Rplus, Rplus.H)):
    check(f"{name}: same active state rotation", A * psi - Rminus * psi)
    check(f"{name}: positive inner product invariant", inner(A*chi, A*psi) - inner(chi, psi))
    check(f"{name}: normalized expectation pullback",
          expectation(Z, A*psi) - expectation(A.H*Z*A, psi), s.cancel)
    check(f"{name}: co-rotated observable expectation invariant",
          expectation(A*Z*A.H, A*psi) - expectation(Z, psi), s.cancel)
    check(f"{name}: conjugation preserves scalar", sc(A*Z*A.H) - sc(Z))
    check(f"{name}: conjugation preserves trace", s.trace(A*Z*A.H) - s.trace(Z))
    check(f"{name}: conjugation preserves determinant", (A*Z*A.H).det() - Z.det(), s.cancel)
    check(f"{name}: conjugation preserves positive norm", norm2(A*Z*A.H) - norm2(Z))

for name, L in (("negative", Lminus), ("positive", Lplus)):
    counterexample(f"{name} boost: state positive norm changes", norm2(L*psi) - norm2(psi))
    counterexample(f"{name} boost: state inner product changes", inner(L*chi, L*psi) - inner(chi, psi))
    counterexample(f"{name} boost: event scalar changes", sc(L*event*L.H) - sc(event))
    counterexample(f"{name} boost: event trace changes", s.trace(L*event*L.H) - s.trace(event))
    counterexample(f"{name} boost: event positive norm changes", norm2(L*event*L.H) - norm2(event))
    check(f"{name} boost: event determinant invariant", (L*event*L.H).det() - event.det(), s.cancel)
    check(f"{name} boost: Minkowski pairing invariant",
          minkowski(L*event*L.H, L*event2*L.H) - minkowski(event, event2))
    check(f"{name} boost: normalized expectation includes density ratio",
          expectation(Z, L*psi) - expectation(L.H*Z*L, psi)/expectation(L.H*L, psi), s.cancel)
    Fp = L * field * adj(L)
    check(f"{name} boost: field similarity preserves scalar", sc(Fp) - sc(field))
    check(f"{name} boost: field similarity preserves trace", s.trace(Fp) - s.trace(field))
    check(f"{name} boost: field similarity preserves determinant", Fp.det() - field.det(), s.cancel)
    counterexample(f"{name} boost: field positive norm changes", norm2(Fp) - norm2(field))

# For X'=T^H X T, the physical frame factor is A=T^H.  This covers a
# non-Hermitian, nonunitary composite, so inverse and adjoint cannot coincide.
for name, T in (("rotation", Rplus), ("boost", Lminus),
                ("composite", (Lminus*Rplus).expand())):
    A, Ti = T.H, adj(T)
    Pp = A * event * A.H
    S = s.diag(Ti, T.H)
    check(f"alternate common {name}: Dirac block covariance", W(Pp)*S - S*W(event))
    check(f"alternate common {name}: mixed chiral pairing invariant",
          inner(Ti*psi, T.H*chi) - inner(psi, chi))
    check(f"alternate common {name}: reverse mixed pairing invariant",
          inner(T.H*chi, Ti*psi) - inner(chi, psi))
    check(f"alternate common {name}: frame determinant invariant", Pp.det()-event.det(), s.cancel)
    check(f"alternate common {name}: direct-factor relabeling", Pp-T.H*event*T)

print(f"\n{checks} convention checks passed in {time.perf_counter()-started:.2f}s.")
