"""Check moving-basis products and derivatives against 2x2 Pauli matrices.

Run: python3 verify_local_basis.py (standard library only).
Random complex fields and nonunitary, noncommuting determinant-one frames
test the body/spatial generator order, time/directional product derivatives,
and both spatial contractions. Central differences independently check the
matrix derivatives along smooth determinant-one frame paths.
"""

import cmath
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


def mul(a, b):
    return (a[0]*b[0] + a[1]*b[2], a[0]*b[1] + a[1]*b[3],
            a[2]*b[0] + a[3]*b[2], a[2]*b[1] + a[3]*b[3])


def inverse(a):
    return scale(1/(a[0]*a[3] - a[1]*a[2]),
                 (a[3], -a[1], -a[2], a[0]))


def commutator(a, b):
    return add(mul(a, b), scale(-1, mul(b, a)))


def paravector(c, z):
    return (c + z[2], z[0] - 1j*z[1], z[0] + 1j*z[1], c - z[2])


def vector(a):
    return ((a[1] + a[2])/2, (a[2] - a[1])/(2j), (a[0] - a[3])/2)


identity = (1, 0, 0, 1)
sigma = tuple(paravector(0, tuple(int(j == k) for j in range(3)))
              for k in range(3))
rng = Random(1104)
checks = 0


def check(name, actual, expected, tolerance=2e-10):
    global checks
    error = max(abs(x-y) for x, y in zip(actual, expected))
    size = 1 + max(abs(x) for x in expected)
    assert error <= tolerance * size, (name, error, actual, expected)
    checks += 1


def random_scalar():
    return complex(rng.uniform(-0.6, 0.6), rng.uniform(-0.6, 0.6))


def random_vector():
    return tuple(random_scalar() for _ in range(3))


def curl(derivatives):
    return (derivatives[1][2]-derivatives[2][1],
            derivatives[2][0]-derivatives[0][2],
            derivatives[0][1]-derivatives[1][0])


for case in range(80):
    # Two shears and a diagonal give a generic det-one complex matrix.
    a, b, d = random_scalar(), random_scalar(), cmath.exp(random_scalar())
    frame = mul(mul((1, a, 0, 1), (1, 0, b, 1)), (d, 0, 0, 1/d))
    inv = inverse(frame)

    def conjugate(matrix):
        return mul(mul(frame, matrix), inv)

    basis = tuple(conjugate(s) for s in sigma)
    c, z, x, y = random_scalar(), random_vector(), random_vector(), random_vector()
    dc = tuple(random_scalar() for _ in range(4))
    dz = tuple(random_vector() for _ in range(4))
    dx = tuple(random_vector() for _ in range(4))
    dy = tuple(random_vector() for _ in range(4))
    # Directions 0..3 are time and the three fixed Cartesian coordinates.
    generators = tuple(paravector(0, random_vector()) for _ in range(4))
    if case % 4 == 0:
        generators = ((0, 0, 0, 0),)*4  # Constant basis.
    if case % 4 == 1:
        dc, dz, dx, dy = (0,)*4, ((0, 0, 0),)*4, ((0, 0, 0),)*4, ((0, 0, 0),)*4

    body = tuple(vector(mul(mul(inv, a), frame)) for a in generators)
    field = conjugate(paravector(c, z))
    xm, ym = conjugate(paravector(0, x)), conjugate(paravector(0, y))
    product_vector = cross(x, y)
    check('pointwise product', mul(xm, ym),
          conjugate(paravector(dot(x, y), scale(1j, product_vector))))

    derivatives, product_derivatives = [], []
    for n in range(4):
        # Differentiate each basis matrix first, then the coefficient fields.
        basis_derivative = tuple(commutator(generators[n], s) for s in basis)
        df = conjugate(paravector(dc[n], dz[n]))
        d_x = conjugate(paravector(0, dx[n]))
        d_y = conjugate(paravector(0, dy[n]))
        for k in range(3):
            df = add(df, scale(z[k], basis_derivative[k]))
            d_x = add(d_x, scale(x[k], basis_derivative[k]))
            d_y = add(d_y, scale(y[k], basis_derivative[k]))
        dp = add(mul(d_x, ym), mul(xm, d_y))
        derivatives.append(df)
        product_derivatives.append(dp)

        check('time/spatial field derivative', df, conjugate(paravector(
            dc[n], add(dz[n], scale(2j, cross(body[n], z))))))
        pv = add(scale(1j, add(cross(dx[n], y), cross(x, dy[n]))),
                 scale(-2, cross(body[n], product_vector)))
        check('time/spatial product derivative', dp, conjugate(paravector(
            dot(dx[n], y) + dot(x, dy[n]), pv)))

        # Smooth path T(h)=exp(h*A)T, with linear coefficient fields.
        av = vector(generators[n])
        q = cmath.sqrt(dot(av, av))

        def sample(h):
            factor = h if abs(q) < 1e-12 else cmath.sinh(h*q)/q
            exp_a = add(scale(cmath.cosh(h*q), identity),
                        scale(factor, generators[n]))
            moving = mul(exp_a, frame)
            moving_inv = inverse(moving)

            def evaluate(scalar, vec, scalar_rate, vec_rate):
                return mul(mul(moving, paravector(scalar+h*scalar_rate,
                           add(vec, scale(h, vec_rate)))), moving_inv)

            return (evaluate(c, z, dc[n], dz[n]),
                    mul(evaluate(0, x, 0, dx[n]), evaluate(0, y, 0, dy[n])))

        h = 1e-5
        plus, minus = sample(h), sample(-h)
        for i, exact in enumerate((df, dp)):
            check('central difference', scale(1/(2*h), add(plus[i], scale(-1, minus[i]))),
                  exact, tolerance=2e-8)

    # A directional derivative is a linear combination of spatial derivatives.
    direction = (0.3, -0.4, 0.5)
    w = tuple(sum(direction[n]*body[n+1][k] for n in range(3)) for k in range(3))
    z_rate = tuple(sum(direction[n]*dz[n+1][k] for n in range(3)) for k in range(3))
    scalar_rate = sum(direction[n]*dc[n+1] for n in range(3))
    direct = tuple(sum(direction[n]*derivatives[n+1][k] for n in range(3)) for k in range(4))
    check('directional field derivative', direct,
          conjugate(paravector(scalar_rate, add(z_rate, scale(2j, cross(w, z))))))
    x_rate = tuple(sum(direction[n]*dx[n+1][k] for n in range(3)) for k in range(3))
    y_rate = tuple(sum(direction[n]*dy[n+1][k] for n in range(3)) for k in range(3))
    direct_product = tuple(sum(direction[n]*product_derivatives[n+1][k]
                               for n in range(3)) for k in range(4))
    product_rate = add(scale(1j, add(cross(x_rate, y), cross(x, y_rate))),
                       scale(-2, cross(w, product_vector)))
    check('directional product derivative', direct_product,
          conjugate(paravector(dot(x_rate, y) + dot(x, y_rate), product_rate)))

    # Contract the actual matrix derivatives, retaining left sigma order.
    contracted = (0,)*4
    contracted_product = (0,)*4
    for n in range(3):
        contracted = add(contracted, mul(basis[n], derivatives[n+1]))
        contracted_product = add(contracted_product, mul(basis[n], product_derivatives[n+1]))

    # Compact form in the closing section: conjugate the fixed-basis
    # derivative plus the commutator with T^{-1} partial T.
    compact, compact_product = (0,)*4, (0,)*4
    fixed_field = paravector(c, z)
    fixed_product = paravector(dot(x, y), scale(1j, product_vector))
    for n in range(3):
        body_matrix = paravector(0, body[n+1])
        field_rate = add(paravector(dc[n+1], dz[n+1]),
                         commutator(body_matrix, fixed_field))
        product_rate = paravector(dot(dx[n+1], y)+dot(x, dy[n+1]),
            scale(1j, add(cross(dx[n+1], y), cross(x, dy[n+1]))))
        product_rate = add(product_rate, commutator(body_matrix, fixed_product))
        compact = add(compact, mul(sigma[n], field_rate))
        compact_product = add(compact_product, mul(sigma[n], product_rate))
    check('compact similarity contraction', contracted, conjugate(compact))
    check('compact product contraction', contracted_product, conjugate(compact_product))

    scalar = sum(dz[n+1][n] + 2j*cross(body[n+1], z)[n] for n in range(3))
    vec = add(dc[1:], scale(1j, curl(dz[1:])))
    for n in range(3):
        vec = add(vec, scale(2, add(scale(body[n+1][n], z), scale(-z[n], body[n+1]))))
    check('contracted field derivative', contracted, conjugate(paravector(scalar, vec)))

    d_cross = tuple(add(cross(dx[n+1], y), cross(x, dy[n+1])) for n in range(3))
    grad_dot = tuple(dot(dx[n+1], y) + dot(x, dy[n+1]) for n in range(3))
    scalar = sum(1j*d_cross[n][n] - 2*cross(body[n+1], product_vector)[n] for n in range(3))
    vec = add(grad_dot, scale(-1, curl(d_cross)))
    for n in range(3):
        vec = add(vec, scale(2j, add(scale(body[n+1][n], product_vector),
                                   scale(-product_vector[n], body[n+1]))))
    check('contracted product derivative', contracted_product, conjugate(paravector(scalar, vec)))
    check('divergence of cross product', (sum(d_cross[n][n] for n in range(3)),),
          (dot(y, curl(dx[1:])) - dot(x, curl(dy[1:])),))
    expanded = add(scale(sum(dy[n+1][n] for n in range(3)), x),
                   scale(-sum(dx[n+1][n] for n in range(3)), y))
    for n in range(3):
        expanded = add(expanded, add(scale(y[n], dx[n+1]), scale(-x[n], dy[n+1])))
    check('curl of cross product', curl(d_cross), expanded)

print(f'PASS: {checks} matrix, finite-difference, and coefficient checks in 80 moving-basis cases.')
