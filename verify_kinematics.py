"""Independent checks for Space--time kinematics; standard library only.

Run: python3 verify_kinematics.py

Exact rational proper-time series inversion checks the acceleration and jerk
formulas. Pauli matrices check the noncommutative identities independently.
Finite differences check explicit trajectories and coordinate/proper-time
conversions; exact polynomial checks cover constant unprojected four-jerk.
"""
from fractions import Fraction as F
from math import asin, asinh, cosh, exp, sinh, sqrt, tanh


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
