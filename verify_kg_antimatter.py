"""Verify the KG charge-conjugation and particle/antiparticle conventions.

Requires SymPy. Differential identities retain arbitrary real, space- and
time-dependent potentials. Oscillator checks use interior occupation states
so the finite-matrix cutoff cannot affect the identities being checked.
"""

import sympy as s


def verify():
    t, x, y, z, q = s.symbols('t x y z q', real=True)
    m = s.symbols('m', positive=True)
    coordinates = (t, x, y, z)
    spatial = coordinates[1:]
    psi = s.Function('u', real=True)(*coordinates) + s.I*s.Function('v', real=True)(*coordinates)
    V = s.Function('V', real=True)(*coordinates)
    A = [s.Function(f'A{k}', real=True)(*coordinates) for k in range(3)]
    count = 0

    def check(name, residual):
        nonlocal count
        entries = list(residual) if isinstance(residual, s.MatrixBase) else [residual]
        def simplify(value):
            # Conjugation commutes with derivatives in these real coordinates.
            value = s.expand(value).replace(
                lambda term: term.func == s.conjugate and isinstance(term.args[0], s.Derivative),
                lambda term: s.diff(s.conjugate(term.args[0].expr), *term.args[0].variables))
            return s.simplify(s.expand(value.doit()))
        assert all(simplify(value) == 0 for value in entries), name
        count += 1

    def energy(field, charge):
        return s.I*s.diff(field, t)-charge*V*field

    def momentum(field, charge, k):
        return -s.I*s.diff(field, spatial[k])-charge*A[k]*field

    def box(field):
        return s.diff(field, t, 2)-sum(s.diff(field, r, 2) for r in spatial)

    def kg(field, charge):
        return -energy(energy(field, charge), charge)+sum(
            momentum(momentum(field, charge, k), charge, k) for k in range(3)
        )+m*m*field

    def rho(field, charge):
        return s.I*charge*(s.conjugate(field)*s.diff(field, t)-s.conjugate(s.diff(field, t))*field)

    def current(field, charge):
        return s.Matrix([-s.I*charge*(s.conjugate(field)*s.diff(field, r)
                          -s.conjugate(s.diff(field, r))*field) for r in spatial])

    check('conjugate energy', s.conjugate(energy(psi, q))+energy(s.conjugate(psi), -q))
    for k in range(3):
        check('conjugate momentum', s.conjugate(momentum(psi, q, k))+momentum(s.conjugate(psi), -q, k))
    check('coupled KG with arbitrary potentials', s.conjugate(kg(psi, q))-kg(s.conjugate(psi), -q))
    check('continuity identity before imposing KG',
          s.diff(rho(psi, q), t)+sum(s.diff(current(psi, q)[k], r) for k, r in enumerate(spatial))
          -s.I*q*(s.conjugate(psi)*box(psi)-s.conjugate(box(psi))*psi))
    check('real density', rho(psi, q)-s.conjugate(rho(psi, q)))
    check('real current', current(psi, q)-s.conjugate(current(psi, q)))
    check('charge reversed at fixed q', rho(s.conjugate(psi), q)+rho(psi, q))
    check('current reversed at fixed q', current(s.conjugate(psi), q)+current(psi, q))
    check('same source after conjugating amplitude and q', current(s.conjugate(psi), -q)-current(psi, q))

    omega = s.symbols('omega', real=True)
    k = s.symbols('k1:4', real=True)
    c = s.symbols('c', complex=True)
    mode = c*s.exp(s.I*(omega*t+2*s.pi*sum(a*b for a, b in zip(k, spatial))))
    check('mode charge', rho(mode, q)+2*q*omega*s.conjugate(c)*c)
    check('mode current', current(mode, q)-4*s.pi*q*s.Matrix(k)*s.conjugate(c)*c)

    # Derive the mode weights from the classical energy, charge, and
    # momentum integrals, including both signs at the same spatial momentum.
    # The constant spatial modulus makes the cell integral its volume.
    E, volume = s.symbols('E volume', positive=True)
    p = s.symbols('p', real=True)
    a, b = s.symbols('a b', complex=True)
    pair = (a*s.exp(-s.I*E*t)+b*s.exp(s.I*E*t))*s.exp(s.I*p*x)/s.sqrt(2*E*volume)
    dt, dx = s.diff(pair, t), s.diff(pair, x)
    density = s.conjugate(dt)*dt+s.conjugate(dx)*dx+(E**2-p**2)*s.conjugate(pair)*pair
    check('normalized two-branch energy', volume*density-E*(s.conjugate(a)*a+s.conjugate(b)*b))
    check('normalized two-branch charge', volume*rho(pair, q)-q*(s.conjugate(a)*a-s.conjugate(b)*b))
    check('normalized two-branch momentum', -volume*(s.conjugate(dt)*dx+s.conjugate(dx)*dt)
          -p*(s.conjugate(a)*a-s.conjugate(b)*b))

    # Positive mode c is an annihilator; negative mode c is a creator.
    # Check both on several non-vacuum states as well as on the vacuum.
    cutoff = 7
    lower = s.zeros(cutoff)
    for n in range(1, cutoff):
        lower[n-1, n] = s.sqrt(n)
    for frequency_sign in (-1, 1):
        coefficient = lower if frequency_sign < 0 else lower.H
        adjoint = coefficient.H
        occupation = adjoint*coefficient if frequency_sign < 0 else coefficient*adjoint
        create = adjoint if frequency_sign < 0 else coefficient
        annihilate = coefficient if frequency_sign < 0 else adjoint
        check('vacuum annihilation', annihilate*s.eye(cutoff)[:, 0])
        for n in range(4):
            state = s.eye(cutoff)[:, n]
            check('mode commutator sign', (coefficient*adjoint-adjoint*coefficient)*state
                  +frequency_sign*state)
            check('occupation spectrum', occupation*state-n*state)
            check('creation adds one quantum', (occupation*create-create*occupation-create)*state)
            E = s.sqrt(m*m+4*s.pi**2*sum(a*a for a in k))
            check('positive energy increment', (E*occupation*create-create*E*occupation-E*create)*state)
            charge = -frequency_sign*q*occupation
            check('opposite charge increment', (charge*create-create*charge+frequency_sign*q*create)*state)
            for p in k:
                total_momentum = -frequency_sign*2*s.pi*p*occupation
                check('opposite mode momentum', (total_momentum*create-create*total_momentum
                      +frequency_sign*2*s.pi*p*create)*state)
    print(f'{count} KG conjugation, current, and oscillator checks passed.')


if __name__ == '__main__':
    verify()
