"""Independent checks for Space--time kinematics; standard library only.

Run: python3 verify_kinematics.py

Exact rational proper-time series inversion checks the acceleration and jerk
formulas. Pauli matrices check the noncommutative identities independently.
Finite differences check explicit trajectories and coordinate/proper-time
conversions; exact polynomial checks cover constant unprojected four-jerk.
Spatial Frenet frames are checked against polynomial paths and independently
differentiated four-velocities.
Matrix boosts independently check inertial-frame relations; numerical
Lagrangian derivatives check momentum and the energy/kinetic-energy formulas.
Nonlinear coordinate charts and noncommuting matrix-valued frames check the
coordinate derivative, product rule and first/second frame derivatives.
"""
from fractions import Fraction as F
from math import acosh, asin, asinh, atanh, cos, cosh, exp, sin, sinh, sqrt, tanh
from random import Random


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def pair(x, y):
    return x[0] * y[0] - dot(x[1:], y[1:])


def cross(x, y):
    return (x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0])


def close(a, b, tol=1e-8):
    assert abs(a-b) <= tol * max(1, abs(a), abs(b)), (a, b)


def derivative(f, s):
    h = 1e-5
    return (f(s-2*h)-8*f(s-h)+8*f(s+h)-f(s+2*h))/(12*h)


def second_derivative(f, s):
    h = 2e-4
    return (-f(s+2*h)+16*f(s+h)-30*f(s)+16*f(s-h)-f(s-2*h))/(12*h*h)


def matrix(x):
    t, a, b, c = map(complex, x)
    return ((t+c, a-1j*b), (a+1j*b, t-c))


def mm(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)) for i in range(2))


def madj(a):
    return ((a[1][1], -a[0][1]), (-a[1][0], a[0][0]))


def mclose(a, b):
    for i in range(2):
        for j in range(2):
            close(a[i][j], b[i][j])


velocities = [
    ((F(0), F(0), F(0)), F(1)),
    ((F(3,5), F(0), F(0)), F(5,4)),
    ((F(1,3), F(2,3), F(0)), F(3,2)),
    ((F(-2,3), F(0), F(1,3)), F(3,2)),
]
# Timelike, null, and spacelike intervals; adjugates use the interval,
# while real proper time is restricted to the timelike case.
for dx in [(2,1,0,0), (1,1,0,0), (1,2,0,0), (-2,0,1,0)]:
    z = matrix(dx)
    mclose(mm(z,madj(z)),matrix([pair(dx,dx),0,0,0]))
    assert pair(dx,dx) == dx[0]**2-dot(dx[1:],dx[1:])
print('PASS: interval and adjugate identities for all causal types.')
accels = [(F(1,5), F(-1,7), F(1,9)), (F(1,3), F(0), F(0)), (F(0), F(1,4), F(0))]
count = 0
for v, g in velocities:
    for a in accels:
        for j in [(F(0),)*3, (F(-1,3), F(2,5), F(3,8)), tuple(-3*g*g*dot(v,a)*x for x in a)]:
            # Starting from r(t)=v*t+a*t²/2+j*t³/6, expand sqrt(1-v(t)²),
            # integrate to s(t), revert to t(s), and compose X(t(s)).
            w0 = 1-dot(v,v)
            w1 = -2*dot(v,a)
            w2 = -dot(a,a)-dot(v,j)
            b1 = 1/g
            b2 = (w1/w0)/(4*g)
            b3 = (w2/(2*w0)-w1*w1/(8*w0*w0))/(3*g)
            c1 = 1/b1
            c2 = -b2*c1*c1/b1
            c3 = -(2*b2*c1*c2+b3*c1**3)/b1
            U = [c1]+[x*c1 for x in v]
            A = [2*c2]+[2*x*c2+y*c1*c1 for x,y in zip(v,a)]
            J = [6*c3]+[6*x*c3+6*y*c1*c2+z*c1**3 for x,y,z in zip(v,a,j)]
            gd = g**3*dot(v,a)
            gdd = g**3*(dot(a,a)+dot(v,j))+3*gd*gd/g
            Ud = [gd]+[g*y+gd*x for x,y in zip(v,a)]
            Udd = [gdd]+[g*z+2*gd*y+gdd*x for x,y,z in zip(v,a,j)]
            assert A == [g*x for x in Ud]
            assert J == [g*g*x+g*gd*y for x,y in zip(Udd,Ud)]
            assert pair(U,U) == 1 and pair(U,A) == 0
            assert pair(U,J) == -pair(A,A)
            assert pair(U,Udd) == -pair(Ud,Ud)
            assert pair(Ud,Ud) == -g*g*(dot(a,a)+g*g*dot(v,a)**2)
            assert pair(Ud,Ud) == -g**4*(dot(a,a)-dot(cross(v,a),cross(v,a)))
            assert pair(A,A) == g*g*pair(Ud,Ud)
            if dot(v,v):
                parallel = [dot(v,a)*x/dot(v,v) for x in v]
                perpendicular = [x-y for x,y in zip(a,parallel)]
                assert pair(A,A) == -g**4*dot(perpendicular,perpendicular)-g**6*dot(parallel,parallel)
                assert Ud[1:] == [g*x+g**3*y for x,y in zip(perpendicular,parallel)]
            projected = [x+pair(A,A)*y for x,y in zip(J,U)]
            assert pair(U,projected) == 0
            assert pair(projected,projected) == pair(J,J)-pair(A,A)**2 <= 0
            assert all(x == 0 for x in projected) == all(x == -3*gd*y/g for x,y in zip(j,a))
            if not dot(v,v):
                assert projected == [0]+list(j)
            # Independently differentiate the explicit acceleration determinant.
            norm_derivative = g*(-4*g**3*gd*dot(a,a)-2*g**4*dot(a,j)
                -6*g**5*gd*dot(v,a)**2-2*g**6*dot(v,a)*(dot(a,a)+dot(v,j)))
            assert norm_derivative == 2*pair(A,J) == 2*pair(A,projected)
            # Matrix multiplication checks adjugation, scalar-free product and commutator.
            um, udm = matrix(U), matrix(Ud)
            expected_product = [0]+[g*g*(complex(x)-1j*complex(y)) for x,y in zip(a,cross(v,a))]
            mclose(mm(madj(um),udm),matrix(expected_product))
            mclose(madj(udm),tuple(tuple(-z for z in row) for row in mm(mm(madj(um),udm),madj(um))))
            first, second = mm(um,udm), mm(udm,um)
            mclose(tuple(tuple(first[i][k]-second[i][k] for k in range(2)) for i in range(2)),
                   matrix([0]+[2j*g*g*x for x in cross(v,a)]))
            # Mass shell, work rate, Euclidean norm and the rest-only unitarity statement.
            mass = F(7,3)
            P = [mass*x for x in U]
            assert pair(P,P) == mass*mass
            assert mass*gd == dot(v,[mass*x for x in Ud[1:]])
            assert dot(U,U) == 2*g*g-1
            mclose(mm(um,um),matrix([2*g*g-1]+[2*g*g*x for x in v]))
            count += 1
print(f'PASS: acceleration, jerk, determinant/projection, Pauli-product and mass-shell identities ({count} exact trajectory cases).')

# Piecewise inertial paths: clock bound and causal endpoints, including equality.
for v,g in velocities:
    for w,h in velocities:
        elapsed = F(2)+F(3)
        dx = [2*g+3*h]+[2*g*x+3*h*y for x,y in zip(v,w)]
        assert pair(dx,dx) >= elapsed*elapsed
        assert (pair(dx,dx) == elapsed*elapsed) == (v == w)
        assert dx[0] > 0
print('PASS: elapsed-clock bound, equality condition and future causal endpoints (16 piecewise inertial paths).')

# Explicit clock and worldline formulas, checked by differentiation.
for initial_accel in [0.2,0.7,1.3]:
    clock = lambda t:(asin(initial_accel*t)+initial_accel*t*sqrt(1-initial_accel**2*t*t))/(2*initial_accel)
    for fraction in [0,0.2,0.8]:
        t = fraction/initial_accel
        close(derivative(clock,t),sqrt(1-initial_accel**2*t*t))
    for s in [-0.6,0,0.8]:
        t = sinh(initial_accel*s)/initial_accel
        r = (cosh(initial_accel*s)-1)/initial_accel
        close((1+initial_accel*r)**2-(initial_accel*t)**2,1)
        close(asinh(initial_accel*t)/initial_accel,s)
        velocity = lambda t:initial_accel*t/sqrt(1+initial_accel**2*t*t)
        close(derivative(velocity,t),initial_accel/(1+initial_accel**2*t*t)**1.5)
print('PASS: constant coordinate-acceleration clock and hyperbolic worldline (18 samples).')

# Arbitrary fixed-axis rapidity histories, with coordinate derivatives computed
# by numerical differentiation and d/dt=(1/cosh(theta))*d/ds.
for theta0, slope, curvature in [(0,0.4,0),(0.3,-0.4,0.2),(0,0,0.7)]:
    theta = lambda s:theta0+slope*s+curvature*s*s/2
    theta_dot = lambda s:(slope+curvature*s)/cosh(theta(s))
    gamma = lambda s:cosh(theta(s))
    velocity = lambda s:tanh(theta(s))
    acceleration = lambda s:(slope+curvature*s)/cosh(theta(s))**3
    for s in [-0.4,0,0.7]:
        g, td = gamma(s), theta_dot(s)
        gd, tdd = derivative(gamma,s)/g, derivative(theta_dot,s)/g
        close(g*g*tdd+g*gd*td,curvature)
        close(derivative(velocity,s)/g,acceleration(s))
        coordinate_jerk = derivative(acceleration,s)/g
        close(g**4*coordinate_jerk+3*g**6*velocity(s)*acceleration(s)**2,curvature)
        U, N = [cosh(theta(s)),sinh(theta(s))], [sinh(theta(s)),cosh(theta(s))]
        Uss = [(g*td)**2*x+curvature*y for x,y in zip(U,N)]
        Nss = [(g*td)**2*x+curvature*y for x,y in zip(N,U)]
        close(pair(Uss,Uss),-pair(Nss,Nss))
        for sign in [-1,1]:
            f = lambda z:exp(sign*theta(z))
            close(derivative(f,s),sign*g*td*f(s))
            close(second_derivative(f,s),((g*td)**2+sign*curvature)*f(s),1e-6)
        # Recover initial dotted data from the polynomial proper-time solution.
    g0 = gamma(0); td0 = theta_dot(0)
    gd0 = sinh(theta0)*td0; tdd0 = derivative(theta_dot,0)/g0
    close(g0*td0,slope)
    close(g0*g0*tdd0+g0*gd0*td0,curvature)
print('PASS: paired U/N derivatives, rapidity conversions, null combinations and constant proper-jerk initial data (9 samples).')

# Exact polynomial normalization in initially resting and boosted frames.
def boost(x, w, g):
    return [g*(x[0]+w*x[1]),g*(x[1]+w*x[0]),x[2],x[3]]

for initial_accel in [F(1,3),F(1),F(5,2)]:
    for w,g in [(F(0),F(1)),(F(3,5),F(5,4)),(F(-4,5),F(5,3))]:
        U = boost([1,0,0,0],w,g)
        A = boost([0,initial_accel,0,0],w,g)
        J = boost([initial_accel**2,0,initial_accel**2,0],w,g)
        gd = A[0]/g
        Ud = [x/g for x in A]
        Udd = [x/g**2-gd*y/g for x,y in zip(J,Ud)]
        assert pair(U,Ud) == 0 and pair(U,Udd) == -pair(Ud,Ud)
        assert g*pair(Ud,Udd)+gd*pair(Ud,Ud) == 0
        assert J == [g*g*x+g*gd*y for x,y in zip(Udd,Ud)]
        # These are every coefficient in det(U+s A+s² J/2)-1.
        coefficients = [pair(U,U)-1,2*pair(U,A),pair(A,A)+pair(U,J),pair(A,J),pair(J,J)/4]
        assert coefficients == [0]*5
        for s in [F(-1,5),F(0),F(2,3)]:
            current_U = [x+s*y+s*s*z/2 for x,y,z in zip(U,A,J)]
            current_A = [x+s*y for x,y in zip(A,J)]
            projected = [x+pair(current_A,current_A)*y for x,y in zip(J,current_U)]
            assert pair(current_A,current_A) == -initial_accel**2
            assert pair(projected,projected) == -initial_accel**4
print('PASS: cubic constant-four-jerk worldline and all normalization constraints (9 exact boosted cases).')

# Constant proper-acceleration magnitude does not require zero projected jerk.
v, g = (F(3,5),0,0), F(5,4)
a, j = (0,F(3,5),0), (F(-3,5),0,0)  # Uniform circular motion, angular speed 1.
U = [g]+[g*x for x in v]
A = [0]+[g*g*x for x in a]
J = [0]+[g**3*x for x in j]
projected = [x+pair(A,A)*y for x,y in zip(J,U)]
assert pair(A,J) == 0 and pair(projected,projected) < 0
print('PASS: uniform circular motion has constant proper-acceleration magnitude and nonzero proper jerk.')

# Spatial Frenet formulas: analytic derivatives of nonplanar polynomial paths
# check the coordinate expressions. Finite differences of U(t), rather than
# differentiated Frenet coefficients, independently check proper acceleration
# and projected jerk. Speed, curvature and torsion vary in these examples.
path_coefficients = [
    ((0, .15, .035, .004), (0, .09, -.04, .02, .002), (0, .04, .018, -.014, .003)),
    ((0, -.18, .025, .012, -.002), (0, .11, .02, -.008), (0, .06, -.015, .01, .004)),
]


def vector_close(x, y, tol=1e-8):
    for a, b in zip(x, y):
        close(a, b, tol)


def combine(coefficients, basis):
    return [sum(c*b[i] for c, b in zip(coefficients, basis)) for i in range(len(basis[0]))]


frenet_count = 0
for coefficients in path_coefficients:
    for scale in [1, 1.6]:
        def path_derivative(t, order):
            values = []
            for component in coefficients:
                terms = list(component)
                for _ in range(order):
                    terms = [i*c for i, c in enumerate(terms)][1:]
                values.append(scale*sum(c*t**i for i, c in enumerate(terms)))
            return values

        def frame(t):
            velocity, acceleration = path_derivative(t, 1), path_derivative(t, 2)
            speed = sqrt(dot(velocity, velocity))
            tangent = [x/speed for x in velocity]
            turning = cross(velocity, acceleration)
            binormal = [x/sqrt(dot(turning, turning)) for x in turning]
            normal = cross(binormal, tangent)
            gamma = 1/sqrt(1-speed*speed)
            U = [gamma]+[gamma*x for x in velocity]
            N = [gamma*speed]+[gamma*x for x in tangent]
            return speed, gamma, (tangent, normal, binormal), (U, N, [0]+list(normal), [0]+binormal)

        for t in [-.6, .2, .8]:
            velocity, acceleration, jerk = [path_derivative(t, n) for n in [1, 2, 3]]
            speed, gamma, triad, rest_basis = frame(t)
            tangent, normal, binormal = triad
            U, N, normal4, binormal4 = rest_basis
            turning = cross(velocity, acceleration)
            turning_norm = sqrt(dot(turning, turning))
            curvature = turning_norm/speed**3
            torsion = dot(turning, jerk)/turning_norm**2
            speed_dot = dot(velocity, acceleration)/speed
            speed_ddot = (dot(acceleration, acceleration)+dot(velocity, jerk)-speed_dot**2)/speed
            curvature_dot = dot(turning, cross(velocity, jerk))/(turning_norm*speed**3)-3*curvature*speed_dot/speed
            gamma_dot = gamma**3*speed*speed_dot
            assert 0 < speed < 1 and curvature > 0
            assert abs(speed_dot) > 1e-6 and abs(curvature_dot) > 1e-6 and abs(torsion) > 1e-6
            vector_close(acceleration, combine([speed_dot, speed**2*curvature, 0], triad))
            vector_close(jerk, combine([speed_ddot-speed**3*curvature**2,
                3*speed*speed_dot*curvature+speed**2*curvature_dot,
                speed**3*curvature*torsion], triad))
            frame_rates = [
                combine([0, speed*curvature, 0], triad),
                combine([-speed*curvature, 0, speed*torsion], triad),
                combine([0, -speed*torsion, 0], triad),
            ]
            for i in range(3):
                vector_close(frame_rates[i], [derivative(lambda z:frame(z)[2][i][k], t) for k in range(3)])
            for i, x in enumerate(rest_basis):
                for k, y in enumerate(rest_basis):
                    close(pair(x, y), (1 if i == 0 else -1) if i == k else 0)
            U_dot = [derivative(lambda z:frame(z)[3][0][k], t) for k in range(4)]
            N_dot = [derivative(lambda z:frame(z)[3][1][k], t) for k in range(4)]
            vector_close(U_dot, combine([0, gamma**2*speed_dot, gamma*speed**2*curvature, 0], rest_basis))
            vector_close(N_dot, combine([gamma**2*speed_dot, 0, gamma*speed*curvature, 0], rest_basis))
            A = [gamma*x for x in U_dot]
            curvature_squared = gamma**6*speed_dot**2+gamma**4*speed**4*curvature**2
            close(pair(A, A), -curvature_squared)
            U_ddot = [second_derivative(lambda z:frame(z)[3][0][k], t) for k in range(4)]
            J = [gamma**2*x+gamma*gamma_dot*y for x, y in zip(U_ddot, U_dot)]
            projected = [x-curvature_squared*y for x, y in zip(J, U)]
            jerk_coefficients = [
                gamma**4*(speed_ddot-speed**3*curvature**2)+3*gamma**3*gamma_dot*speed_dot,
                3*gamma**5*speed*speed_dot*curvature+gamma**3*speed**2*curvature_dot,
                gamma**3*speed**3*curvature*torsion,
            ]
            vector_close(projected, combine([0]+jerk_coefficients, rest_basis), 2e-7)
            close(pair(projected, projected), -dot(jerk_coefficients, jerk_coefficients), 2e-7)
            frenet_count += 1
print(f'PASS: spatial Frenet frame, varying curvature/torsion, rest basis and proper jerk ({frenet_count} nonplanar samples).')

# Degenerate spatial-curvature limit: no normal is defined by a straight path,
# but the surviving tangent term agrees exactly with fixed-axis proper jerk.
speed, gamma = F(3, 5), F(5, 4)
speed_dot, speed_ddot = F(2, 7), F(-1, 9)
gamma_dot = gamma**3*speed*speed_dot
gamma_ddot = gamma**3*(speed_dot**2+speed*speed_ddot)+3*gamma_dot**2/gamma
U, N = [gamma, gamma*speed], [gamma*speed, gamma]
U_dot = [gamma_dot, gamma_dot*speed+gamma*speed_dot]
U_ddot = [gamma_ddot, gamma_ddot*speed+2*gamma_dot*speed_dot+gamma*speed_ddot]
J = [gamma**2*x+gamma*gamma_dot*y for x, y in zip(U_ddot, U_dot)]
projected = [x-gamma**6*speed_dot**2*y for x, y in zip(J, U)]
assert projected == [(gamma**4*speed_ddot+3*gamma**3*gamma_dot*speed_dot)*x for x in N]

# The circle has constant spatial curvature 1/r and speed, but nonzero proper
# jerk from the rotating normal. Angular speed is one, hence radius = speed.
U, N = [gamma, gamma*speed, 0, 0], [gamma*speed, gamma, 0, 0]
A, J = [0, 0, gamma**2*speed, 0], [0, -gamma**3*speed, 0, 0]
curvature = 1/speed
projected = [x+pair(A, A)*y for x, y in zip(J, U)]
assert projected == [-gamma**4*speed**3*curvature**2*x for x in N]
assert pair(projected, projected) == -gamma**8*speed**6*curvature**4
print('PASS: Frenet straight-path limit and circular proper jerk (exact rational checks).')

# General boosts: independent matrix congruences check the component formulas.
# A fixed seed supplies oblique directions; zero and both signs of frame speed
# are included explicitly, along with resting, moving and null tangents.
def inertial_transform(x, direction, beta):
    gamma = 1/sqrt(1-beta*beta)
    longitudinal = dot(direction, x[1:])
    return [gamma*(x[0]-beta*longitudinal)]+[
        r+((gamma-1)*longitudinal-gamma*beta*x[0])*u
        for r, u in zip(x[1:], direction)]


rng = Random(7316)
frame_count = 40
for index in range(frame_count):
    direction = [rng.uniform(-1, 1) for _ in range(3)]
    direction = [x/sqrt(dot(direction, direction)) for x in direction]
    beta = [-.9, -.4, 0, .3, .85][index % 5]
    gamma = 1/sqrt(1-beta*beta)
    event = [rng.uniform(-3, 3) for _ in range(4)]
    transformed = inertial_transform(event, direction, beta)
    half_cosh = sqrt((gamma+1)/2)
    factor = matrix([half_cosh]+[-gamma*beta*x/(2*half_cosh) for x in direction])
    mclose(mm(mm(factor, matrix(event)), factor), matrix(transformed))
    close(pair(transformed, transformed), pair(event, event))
    vector_close(inertial_transform(transformed, direction, -beta), event)
    longitudinal = dot(direction, event[1:])
    for sign in [-1, 1]:
        close(transformed[0]+sign*dot(direction, transformed[1:]),
              gamma*(1-sign*beta)*(event[0]+sign*longitudinal))

    raw_velocity = [rng.uniform(-1, 1) for _ in range(3)]
    speed = 0 if index % 10 == 0 else rng.uniform(.1, .97)
    velocity = [speed*x/sqrt(dot(raw_velocity, raw_velocity)) for x in raw_velocity]
    particle_gamma = 1/sqrt(1-speed*speed)
    longitudinal_velocity = dot(direction, velocity)
    clock_factor = gamma*(1-beta*longitudinal_velocity)
    boosted_velocity = [
        ((v-longitudinal_velocity*u)/gamma+(longitudinal_velocity-beta)*u)
        /(1-beta*longitudinal_velocity) for v, u in zip(velocity, direction)]
    increment = inertial_transform([1]+velocity, direction, beta)
    close(increment[0], clock_factor)
    vector_close(boosted_velocity, [x/increment[0] for x in increment[1:]])
    close(1-dot(boosted_velocity, boosted_velocity), (1-speed*speed)/clock_factor**2)
    close(1/sqrt(1-dot(boosted_velocity, boosted_velocity)), particle_gamma*clock_factor)
    ray = inertial_transform([1, 0, 1, 0], direction, beta)
    close(dot(ray[1:], ray[1:]), ray[0]**2)
    assert ray[0] > 0

    # A carried clock and a rod measured at equal primed times.
    vector_close(inertial_transform([1]+[beta*u for u in direction], direction, beta),
                 [1/gamma, 0, 0, 0])
    rod = inertial_transform([beta*longitudinal]+event[1:], direction, beta)
    close(rod[0], 0)
    vector_close(rod[1:], [r+(1/gamma-1)*longitudinal*u for r, u in zip(event[1:], direction)])
    close(dot(rod[1:], rod[1:]), dot(event[1:], event[1:])-beta**2*longitudinal**2)

    second_beta = rng.uniform(-.9, .9)
    combined_beta = (beta+second_beta)/(1+beta*second_beta)
    vector_close(inertial_transform(transformed, direction, second_beta),
                 inertial_transform(event, direction, combined_beta))
    close(1/sqrt(1-combined_beta**2), gamma*(1+beta*second_beta)/sqrt(1-second_beta**2))

    # Momentum from finite differences of L, independently of the mass shell.
    mass = rng.uniform(.1, 4)
    energy = mass*particle_gamma
    momentum = [energy*v for v in velocity]
    for component in range(3):
        def lagrangian_with_component(value):
            varied = list(velocity)
            varied[component] = value
            return -mass*sqrt(1-dot(varied, varied))
        close(derivative(lagrangian_with_component, velocity[component]), momentum[component])
    close(dot(velocity, momentum)+mass/particle_gamma, energy)
    close(energy-mass, dot(momentum, momentum)/(energy+mass))
    boosted_momentum = inertial_transform([energy]+momentum, direction, beta)
    close(pair(boosted_momentum, boosted_momentum), mass*mass)
    close(boosted_momentum[0]/energy, clock_factor)
    vector_close([x/boosted_momentum[0] for x in boosted_momentum[1:]], boosted_velocity)
print(f'PASS: inertial boosts, velocity/energy clock factor, light speed, lengths, composition and free-particle action ({frame_count} deterministic cases).')


def msum(*matrices):
    return tuple(tuple(sum(m[i][j] for m in matrices) for j in range(2)) for i in range(2))


def mscale(scalar, a):
    return tuple(tuple(scalar*x for x in row) for row in a)


def mhermitian(a):
    return tuple(tuple(a[j][i].conjugate() for j in range(2)) for i in range(2))


def commutator(a, b):
    return msum(mm(a, b), mscale(-1, mm(b, a)))


def matrix_derivative(function, parameter, order=1):
    differentiate = derivative if order == 1 else second_derivative
    return tuple(tuple(differentiate(lambda s:function(s)[i][j], parameter)
                       for j in range(2)) for i in range(2))


def matrix_close(a, b, tolerance=1e-8):
    for row_a, row_b in zip(a, b):
        vector_close(row_a, row_b, tolerance)


sigmas = [matrix([int(i == j) for i in range(4)]) for j in range(4)]

# Polynomial fields with independent complex coefficients. Their partials are
# evaluated analytically; finite differences of the matrix product provide the
# independent left side of the noncommutative differential product rule.
field_terms = [
    [((1, .2j, -.3, .4), (2, 0, 0, 0)),
     ((.1j, 1, .3, -.2j), (0, 1, 1, 0)),
     ((.2, -.3j, .5j, 1), (0, 0, 0, 2)),
     ((1, 0, .2j, -.3), (0, 1, 0, 0))],
    [((.2, .3j, -.4, 1), (1, 0, 0, 1)),
     ((.5j, 1, -.3j, .2), (0, 0, 2, 0)),
     ((.3, .4, 1, -.5j), (0, 1, 0, 0)),
     ((1, -.2j, .4, .5), (0, 0, 0, 0))],
]


def polynomial_field(coordinates, terms, partial=None):
    result = matrix([0, 0, 0, 0])
    for coefficients, exponents in terms:
        powers = list(exponents)
        weight = 1
        if partial is not None:
            weight = powers[partial]
            if not weight:
                continue
            powers[partial] -= 1
        for coordinate, power in zip(coordinates, powers):
            weight *= coordinate**power
        result = msum(result, mscale(weight, matrix(coefficients)))
    return result


for point in [(-.4, .2, .5, -.1), (.3, -.2, .1, .6), (.7, .4, -.3, .2)]:
    first, second = [polynomial_field(point, terms) for terms in field_terms]
    partials = [[polynomial_field(point, terms, k) for k in range(4)] for terms in field_terms]
    operator_values = [msum(*(mm(sigma, value) for sigma, value in zip(sigmas, values)))
                       for values in partials]
    direct_partials = []
    for component in range(4):
        def product_along_coordinate(value):
            varied = list(point)
            varied[component] = value
            return mm(*(polynomial_field(varied, terms) for terms in field_terms))
        direct_partials.append(matrix_derivative(product_along_coordinate, point[component]))
    direct = msum(*(mm(sigma, value) for sigma, value in zip(sigmas, direct_partials)))
    naive_rule = msum(mm(operator_values[0], second), mm(first, operator_values[1]))
    correction = msum(*(mm(commutator(sigmas[k], first), partials[1][k]) for k in range(1, 4)))
    matrix_close(direct, msum(naive_rule, correction))
    assert max(abs(value) for row in correction for value in row) > .01
print('PASS: noncommutative differential product rule (3 complex polynomial field pairs).')

# A determinant-one, noncommuting boost-x followed by rotation-z product.
# Matrix derivatives are taken directly before the frame identities are used.
def varying_factor(s):
    rapidity = .3*s+.1*s*s
    angle = .2+.4*s+.05*s*s
    boost_factor = matrix([cosh(rapidity/2), sinh(rapidity/2), 0, 0])
    rotation_factor = matrix([cos(angle/2), 0, 0, -1j*sin(angle/2)])
    return mm(boost_factor, rotation_factor)


def frame_omega(s):
    return mm(matrix_derivative(varying_factor, s), madj(varying_factor(s)))


def frame_field(s):
    return polynomial_field((s, .2+s*s, -.3+.4*s, .5-.2*s), field_terms[0])


def rest_tangent(s):
    inverse = madj(varying_factor(s))
    return mm(inverse, mhermitian(inverse))


for s in [-.6, -.2, 0, .3, .8]:
    factor = varying_factor(s)
    inverse, hermitian = madj(factor), mhermitian(factor)
    matrix_close(mm(factor, inverse), sigmas[0])
    omega = frame_omega(s)
    # Differentiate T directly to avoid nested finite-difference cancellation.
    omega_derivative = msum(mm(matrix_derivative(varying_factor, s, 2), inverse),
                           mscale(-1, mm(omega, omega)))
    close(omega[0][0]+omega[1][1], 0)
    value_dot = matrix_derivative(frame_field, s)
    value_ddot = matrix_derivative(frame_field, s, 2)
    for congruence in [False, True]:
        def transformed_field(z):
            transform = varying_factor(z)
            right = mhermitian(transform) if congruence else madj(transform)
            return mm(mm(transform, frame_field(z)), right)
        right = hermitian if congruence else inverse
        transformed = transformed_field(s)
        first = matrix_derivative(transformed_field, s)
        second = matrix_derivative(transformed_field, s, 2)
        physical_first = mm(mm(factor, value_dot), right)
        physical_second = mm(mm(factor, value_ddot), right)
        if congruence:
            omega_h, omega_derivative_h = mhermitian(omega), mhermitian(omega_derivative)
            expected_first = msum(physical_first, mm(omega, transformed), mm(transformed, omega_h))
            expected_second = msum(physical_second,
                mscale(2, mm(omega, first)), mscale(2, mm(first, omega_h)),
                mm(msum(omega_derivative, mscale(-1, mm(omega, omega))), transformed),
                mm(transformed, msum(omega_derivative_h, mscale(-1, mm(omega_h, omega_h)))),
                mscale(-2, mm(mm(omega, transformed), omega_h)))
        else:
            expected_first = msum(physical_first, commutator(omega, transformed))
            expected_second = msum(physical_second, mscale(2, commutator(omega, first)),
                commutator(omega_derivative, transformed),
                mscale(-1, commutator(omega, commutator(omega, transformed))))
        matrix_close(first, expected_first)
        matrix_close(second, expected_second, 2e-6)
    matrix_close(mm(mm(factor, rest_tangent(s)), hermitian), sigmas[0])
    local_acceleration = mm(mm(factor, matrix_derivative(rest_tangent, s)), hermitian)
    matrix_close(local_acceleration, mscale(-1, msum(omega, mhermitian(omega))))
    acceleration_det = (local_acceleration[0][0]*local_acceleration[1][1]
                        -local_acceleration[0][1]*local_acceleration[1][0])
    local_jerk = msum(mm(mm(factor, matrix_derivative(rest_tangent, s, 2)), hermitian),
                      mscale(acceleration_det, sigmas[0]))
    matrix_close(local_jerk, msum(mscale(-1, msum(omega_derivative, mhermitian(omega_derivative))),
                                 commutator(omega, mhermitian(omega))), 2e-6)
print('PASS: varying noncommuting similarity/congruence derivatives and comoving acceleration/jerk (5 frame samples).')


def inverse_real_matrix(a):
    size = len(a)
    augmented = [list(row)+[F(i == j) for j in range(size)] for i, row in enumerate(a)]
    for column in range(size):
        pivot = next(row for row in range(column, size) if augmented[row][column])
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value/divisor for value in augmented[column]]
        for row in range(size):
            if row != column:
                weight = augmented[row][column]
                augmented[row] = [a-weight*b for a, b in zip(augmented[row], augmented[column])]
    return [row[size:] for row in augmented]


# A nonlinear chart mixes inertial time and space. Exact rational Jacobians
# check the reciprocal metric formula, including off-diagonal metric entries.
def nonlinear_chart(q):
    t, x, y, z = q
    return [t+F(2, 25)*x*x, (1+F(1, 5)*t)*x,
            y+F(3, 20)*t*x, z+F(1, 10)*x*y]


for q in [(F(-1, 2), F(1, 3), F(1, 5), F(0)),
          (F(0), F(-2, 5), F(1, 2), F(1, 3)),
          (F(2, 3), F(1, 4), F(-1, 3), F(-1, 5))]:
    t, x, y, z = q
    jacobian = [[F(1), F(4, 25)*x, F(0), F(0)],
                [x/5, 1+t/5, F(0), F(0)],
                [F(3, 20)*x, F(3, 20)*t, F(1), F(0)],
                [F(0), y/10, x/10, F(1)]]
    basis = [list(column) for column in zip(*jacobian)]
    metric = [[pair(a, b) for b in basis] for a in basis]
    inverse_metric = inverse_real_matrix(metric)
    reciprocal = [[sum(inverse_metric[j][k]*basis[j][a]*(1 if a == 0 else -1)
                       for j in range(4)) for a in range(4)] for k in range(4)]
    for a in range(4):
        assert [sum(reciprocal[k][b]*jacobian[a][k] for k in range(4))
                for b in range(4)] == [int(a == b) for b in range(4)]
    partials = [polynomial_field(nonlinear_chart(q), field_terms[0], a) for a in range(4)]
    chart_partials = [msum(*(mscale(jacobian[a][k], partials[a]) for a in range(4))) for k in range(4)]
    matrix_close(msum(*(mm(matrix(reciprocal[k]), chart_partials[k]) for k in range(4))),
                 msum(*(mm(sigmas[a], partials[a]) for a in range(4))))


def chart_worldline(t):
    return nonlinear_chart([t, .2*t+.05*t*t, .1-.1*t, .2*t**3])


def chart_velocity(t):
    rate = [derivative(lambda z:chart_worldline(z)[a], t) for a in range(4)]
    return [x/rate[0] for x in rate[1:]]


for t in [-.5, 0, .6]:
    rate = [derivative(lambda z:chart_worldline(z)[a], t) for a in range(4)]
    second_rate = [second_derivative(lambda z:chart_worldline(z)[a], t) for a in range(4)]
    velocity = chart_velocity(t)
    gamma = 1/sqrt(1-dot(velocity, velocity))
    proper_rate = sqrt(pair(rate, rate))
    assert rate[0] > 0
    close(proper_rate, rate[0]/gamma)
    vector_close([x/proper_rate for x in rate], [gamma]+[gamma*x for x in velocity])
    vector_close([derivative(lambda z:chart_velocity(z)[a], t) for a in range(3)],
                 [(a*rate[0]-v*second_rate[0])/rate[0]**2
                  for a, v in zip(second_rate[1:], rate[1:])], 2e-6)
print('PASS: nonlinear chart reciprocal derivatives, proper-time rates and coordinate acceleration (3 charts and 3 path samples).')

# Fixed Rindler position: the chart-time clock rate is the position, while
# proper acceleration is its reciprocal. Differentiate the worldline directly.
for radius in [.4, 1, 2.5]:
    for t in [-.3, 0, .5]:
        worldline = lambda z:[radius*sinh(z), radius*cosh(z), 0, 0]
        rate = [derivative(lambda z:worldline(z)[a], t) for a in range(4)]
        acceleration = [second_derivative(lambda z:worldline(z)[a], t)/radius**2 for a in range(4)]
        close(sqrt(pair(rate, rate)), radius)
        close(pair(acceleration, acceleration), -1/radius**2, 2e-6)
print('PASS: Rindler clock rates and proper acceleration (9 samples).')

# Two joined hyperbolic burns: differentiate the complete trajectory in proper
# time and independently evaluate coordinate-time derivatives on both halves.
rocket_samples = 0
for initial_acceleration in [.2, .7, 1.3]:
    for distance in [.1, 4, 25]:
        peak_gamma = 1+initial_acceleration*distance/2
        half_proper = acosh(peak_gamma)/initial_acceleration
        half_coordinate = sqrt(distance/initial_acceleration+distance**2/4)

        def rocket_theta(s):
            return initial_acceleration*(s if s <= half_proper else 2*half_proper-s)

        def rocket_worldline(s):
            theta = rocket_theta(s)
            elapsed = sinh(theta)/initial_acceleration
            position = (cosh(theta)-1)/initial_acceleration
            if s > half_proper:
                elapsed, position = 2*half_coordinate-elapsed, distance-position
            return [elapsed, position, 0, 0]

        def rocket_velocity(t):
            remaining = t if t <= half_coordinate else 2*half_coordinate-t
            return tanh(asinh(initial_acceleration*remaining))

        def rocket_u(s):
            theta = rocket_theta(s)
            return [cosh(theta), sinh(theta), 0, 0]

        vector_close(rocket_worldline(0), [0, 0, 0, 0])
        vector_close(rocket_worldline(half_proper), [half_coordinate, distance/2, 0, 0])
        vector_close(rocket_worldline(2*half_proper), [2*half_coordinate, distance, 0, 0])
        vector_close(rocket_u(0), rocket_u(2*half_proper))
        close(rocket_velocity(half_coordinate), sqrt(1-peak_gamma**-2))
        # Evaluate both analytic branches at the join, without differencing
        # across the acceleration jump.
        elapsed = sinh(initial_acceleration*half_proper)/initial_acceleration
        position = (cosh(initial_acceleration*half_proper)-1)/initial_acceleration
        close(elapsed, 2*half_coordinate-elapsed)
        close(position, distance-position)
        assert 2*half_coordinate > distance and half_coordinate > half_proper
        for fraction in [.1, .4, .7, 1.3, 1.6, 1.9]:
            s = fraction*half_proper
            sign = 1 if fraction < 1 else -1
            t, r, _, _ = rocket_worldline(s)
            gamma, spatial_u, _, _ = rocket_u(s)
            speed = spatial_u/gamma
            tangent = [derivative(lambda z:rocket_worldline(z)[a], s) for a in range(4)]
            acceleration = [second_derivative(lambda z:rocket_worldline(z)[a], s) for a in range(4)]
            expected_acceleration = [sign*initial_acceleration*spatial_u,
                                     sign*initial_acceleration*gamma, 0, 0]
            vector_close(tangent, rocket_u(s))
            vector_close(acceleration, expected_acceleration, 3e-6)
            close(pair(rocket_u(s), rocket_u(s)), 1)
            close(pair(expected_acceleration, expected_acceleration), -initial_acceleration**2)
            close(rocket_velocity(t), speed)
            close(derivative(rocket_velocity, t), sign*initial_acceleration/gamma**3)
            close(second_derivative(rocket_velocity, t),
                  -3*initial_acceleration**2*speed/gamma**4, 2e-6)
            for component in range(2):
                coordinate_u = lambda z: (1 if component == 0 else rocket_velocity(z))/sqrt(1-rocket_velocity(z)**2)
                close(derivative(coordinate_u, t),
                      sign*initial_acceleration*(speed if component == 0 else 1), 1e-7)
                close(second_derivative(lambda z:rocket_u(z)[component], s),
                      initial_acceleration**2*rocket_u(s)[component], 3e-6)
            reflected = rocket_worldline(2*half_proper-s)
            close(t+reflected[0], 2*half_coordinate)
            close(r+reflected[1], distance)
            rocket_samples += 1

# One onboard hour with a .999c midpoint peak: solve for acceleration, then
# independently recover both clocks and the distance from the SI trip formulas.
# SI clock conversion uses s=c*tau, and natural-unit t as a length.
c_light = 299792458
standard_gravity = 9.80665
trip_proper_seconds = 3600
peak_fraction = .999
peak_rapidity = atanh(peak_fraction)
si_acceleration = 2*c_light*peak_rapidity/trip_proper_seconds
inverse_length_acceleration = si_acceleration/c_light**2
trip_distance = 2*c_light**2/si_acceleration*(cosh(peak_rapidity)-1)
si_coordinate_time = sqrt(trip_distance**2/c_light**2+4*trip_distance/si_acceleration)
si_proper_time = 2*c_light/si_acceleration*acosh(1+si_acceleration*trip_distance/(2*c_light**2))
close(si_proper_time, trip_proper_seconds)
close(si_coordinate_time, 2*c_light/si_acceleration*sinh(peak_rapidity))
close(sqrt(1-(1+si_acceleration*trip_distance/(2*c_light**2))**-2), peak_fraction)
close(2*sqrt(trip_distance/inverse_length_acceleration+trip_distance**2/4)/c_light, si_coordinate_time)
close(2*acosh(1+inverse_length_acceleration*trip_distance/2)/inverse_length_acceleration/c_light, si_proper_time)
assert round(si_coordinate_time/3600, 2) == 5.88
assert round(trip_distance/(c_light*3600), 2) == 5.62
assert round(si_acceleration/standard_gravity) == 64541
print(f'PASS: accelerate/brake trip, both clocks, endpoint/midpoint matching, derivatives and SI example ({rocket_samples} trajectory samples).')
