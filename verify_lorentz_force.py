"""Exact checks of the Lorentz-force projections and proper acceleration.

Run: python3 verify_lorentz_force.py (standard library only).
Independent field-tensor boosts check the rest-frame electric field. The
derivative of relativistic momentum and its Minkowski norm check the proposed
acceleration, in natural and SI units, including rest and force cancellation.
"""

from fractions import Fraction as Q
from random import Random


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def scale(s, a):
    return tuple(s * x for x in a)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def product(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(4))
                       for j in range(4)) for i in range(4))


def boost(beta, gamma):
    speed2 = dot(beta, beta)
    if not speed2:
        return tuple(tuple(Q(i == j) for j in range(4)) for i in range(4))
    return ((gamma,) + scale(-gamma, beta),) + tuple(
        (-gamma * beta[i],) + tuple(
            Q(i == j) + (gamma-1)*beta[i]*beta[j]/speed2
            for j in range(3)) for i in range(3))


def field_tensor(electric, magnetic):
    # Mixed tensor acting on (gamma, gamma*beta), using (E, cB).
    ex, ey, ez = electric
    bx, by, bz = magnetic
    return ((0, ex, ey, ez), (ex, 0, bz, -by),
            (ey, -bz, 0, bx), (ez, by, -bx, 0))


def project(vector, axis):
    parallel = scale(dot(vector, axis)/dot(axis, axis), axis)
    return parallel, add(vector, scale(-1, parallel))


rng = Random(12022)
# Rational half-rapidity coordinates give exact beta and gamma.
parameters = [(Q(0),)*3, (Q(1, 3), 0, 0), (0, Q(-2, 3), 0),
              (Q(99, 100), 0, 0)]
parameters += [tuple(Q(rng.randint(-8, 8), 20) for _ in range(3))
               for _ in range(12)]
mass = Q(7, 3)
checks = 0
for parameter in parameters:
    r2 = dot(parameter, parameter)
    beta = scale(2/(1+r2), parameter)
    gamma = (1+r2)/(1-r2)
    assert gamma**2 * (1-dot(beta, beta)) == 1
    for light_speed in (Q(1), Q(299792458)):
        velocity = scale(light_speed, beta)
        magnetic = scale(1/light_speed, (Q(2, 3), Q(-1), Q(5, 4)))
        electric = (Q(-4, 5), Q(3, 2), Q(7, 6))
        fields = [(electric, (Q(0),)*3), ((Q(0),)*3, magnetic),
                  (electric, magnetic), (scale(-1, cross(velocity, magnetic)), magnetic),
                  ((Q(0),)*3, (Q(0),)*3)]
        for electric, magnetic in fields:
            rest_tensor = product(product(boost(beta, gamma),
                field_tensor(electric, scale(light_speed, magnetic))),
                boost(scale(-1, beta), gamma))
            rest_electric = tuple(rest_tensor[i][0] for i in range(1, 4))
            for charge in (Q(-5, 2), Q(0), Q(4, 3)):
                lorentz = add(electric, cross(velocity, magnetic))
                if dot(velocity, velocity):
                    e_parallel, e_perpendicular = project(electric, velocity)
                    transverse = add(e_perpendicular, cross(velocity, magnetic))
                    a_parallel = scale(charge/(mass*gamma**3), e_parallel)
                    a_perpendicular = scale(charge/(mass*gamma), transverse)
                    acceleration = add(a_parallel, a_perpendicular)
                    assert project(acceleration, velocity) == (a_parallel, a_perpendicular)
                    assert rest_electric == add(e_parallel, scale(gamma, transverse))
                    alpha2 = charge**2/mass**2 * (
                        dot(e_parallel, e_parallel) + gamma**2*dot(transverse, transverse))
                    assert alpha2 == gamma**6*dot(a_parallel, a_parallel) + gamma**4*dot(a_perpendicular, a_perpendicular)
                else:
                    acceleration = scale(charge/mass, electric)
                    alpha2 = dot(acceleration, acceleration)
                    assert rest_electric == electric

                unsplit = scale(charge/(mass*gamma), add(lorentz,
                    scale(-dot(velocity, electric)/light_speed**2, velocity)))
                assert acceleration == unsplit
                gamma_dot = gamma**3 * dot(velocity, acceleration)/light_speed**2
                momentum_dot = scale(mass, add(scale(gamma, acceleration),
                                               scale(gamma_dot, velocity)))
                assert momentum_dot == scale(charge, lorentz)
                assert mass*light_speed**2*gamma_dot == charge*dot(velocity, electric)
                if dot(velocity, velocity):
                    force_parallel, force_perpendicular = project(momentum_dot, velocity)
                    assert force_parallel == scale(charge, e_parallel)
                    assert force_perpendicular == scale(charge, transverse)

                # U is dimensionless; d/ds = gamma/c d/dt in SI units.
                u_dot = (gamma_dot,) + scale(1/light_speed,
                    add(scale(gamma, acceleration), scale(gamma_dot, velocity)))
                u_s = scale(gamma/light_speed, u_dot)
                determinant = u_s[0]**2 - dot(u_s[1:], u_s[1:])
                assert -light_speed**4*determinant == alpha2
                assert alpha2 == charge**2/mass**2*dot(rest_electric, rest_electric)
                assert alpha2 == charge**2*gamma**2/mass**2 * (
                    dot(lorentz, lorentz) - dot(velocity, electric)**2/light_speed**2)
                checks += 1

print(f'PASS: Lorentz projections, momentum/work, rest-frame field, and proper acceleration ({checks} exact natural/SI cases).')
