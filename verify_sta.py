"""Independent exact matrix/operator checks for sta_notes.tex (requires SymPy)."""

import sympy as s

i = s.I
one = s.eye(2)
sigma = (
    s.Matrix([[0, 1], [1, 0]]),
    s.Matrix([[0, -i], [i, 0]]),
    s.diag(1, -1),
)
checks = 0


def vec(v):
    return sum((v[k] * sigma[k] for k in range(3)), s.zeros(2))


def sc(z):
    return s.trace(z) / 2


def adj(z):
    return s.trace(z) * one - z


def check(name, expression):
    global checks
    entries = list(expression) if isinstance(expression, s.MatrixBase) else [expression]
    assert all(s.simplify(v) == 0 for v in entries), name
    checks += 1
    print(f"PASS {name}")


# Matrix arithmetic is independent of the paravector product formulas.
A = s.Matrix([2, -3, 5])
B = s.Matrix([-1, 4, 2])
check("sigma dot/cross product", vec(A) * vec(B) - A.dot(B) * one - i * vec(A.cross(B)))
Z = (2 + 3 * i) * one + vec(A + i * B)
S = i * sigma[1]
check("sigma conjugation", S.inv() * Z.conjugate() * S - ((2 - 3*i)*one - vec(A-i*B)))
W = lambda z: s.BlockMatrix([[s.zeros(2), z], [adj(z), s.zeros(2)]]).as_explicit()
check("Dirac block adjoint", W(Z).H - W(S.inv()*Z.conjugate()*S))
check("algebraic Dirac square", W(Z)**2 - Z.det()*s.eye(4))

# Ideal normalization, Born probabilities and the causal Dirac current.
v = s.Matrix([1 + 2*i, 3-i])
w = s.Matrix([-2+i, 1+3*i])
psi = s.sqrt(2) * v * s.Matrix([[1, 0]])
phi = s.sqrt(2) * w * s.Matrix([[1, 0]])
n = (v.H*v)[0]
state = psi*psi.H/n
u = s.Matrix([s.Rational(2, 3), s.Rational(1, 3), s.Rational(2, 3)])
projector = (one + vec(u))/2
check("minimal ideal norm", sc(psi.H*psi) - n)
check("pure state determinant", state.det())
check("Born pairing", sc(projector*state) - (v.H*projector*v)[0]/n)
J = adj(psi*psi.H) + phi*phi.H
check("Dirac current determinant", J.det() - 4*s.Abs((v.H*w)[0])**2)

# Rational rapidity: cosh(theta/2)=5/4, sinh(theta/2)=3/4.
L = s.Rational(5, 4)*one + s.Rational(3, 4)*sigma[2]
Li = L.inv()
X = 7*one + vec(A)
F = vec(A + i*B)
Xp = Li*X*Li
Fp = Li*F*L
check("real Lorentz congruence", Xp.H-Xp)
check("Lorentz interval", Xp.det()-X.det())
check("boost time component", sc(Xp) - (s.Rational(17,8)*7-s.Rational(15,8)*A[2]))
Jp = adj(L*psi*psi.H*L) + Li*phi*phi.H*Li
check("Dirac current covariance", Jp-Li*J*Li)
force = (F*X + X*F.H)/2
forcep = (Fp*Xp + Xp*Fp.H)/2
check("Lorentz force covariance", forcep-Li*force*Li)
dual = (F*X-X*F.H)/(2*i)
check("electric force components", force - (A.dot(A)*one + vec(7*A+A.cross(B))))
check("magnetic dual force components", dual - (A.dot(B)*one + vec(7*B-A.cross(A))))
check("force mass-shell tangency", sc(adj(X)*force))
check("dual mass-shell tangency", sc(adj(X)*dual))

# Free modes in the standard and chiral bases, with non-collinear spin.
p = s.Matrix([2, -1, 2])
m = s.Integer(4)
energy = s.Integer(5)
K = vec(p)
H = s.BlockMatrix([[-K, m*one], [m*one, K]]).as_explicit()
C = s.BlockMatrix([[one, one], [-one, one]]).as_explicit()/s.sqrt(2)
Hs = s.BlockMatrix([[m*one, K], [K, -m*one]]).as_explicit()
check("standard Dirac basis", C*H*C.H-Hs)
e = s.Matrix([s.sqrt(s.Rational(2,3)), (1+i)/s.sqrt(6)])
positive = s.sqrt((energy+m)/(2*energy))*s.Matrix.vstack(e,K*e/(energy+m))
negative = s.sqrt((energy+m)/(2*energy))*s.Matrix.vstack(-K*e/(energy+m),e)
check("positive free mode", Hs*positive-energy*positive)
check("negative free mode", Hs*negative+energy*negative)
check("free mode normalization", (positive.H*positive)[0]-1)
check("energy orthogonality", (positive.H*negative)[0])
Q = (s.eye(4)+H/energy)/2
check("Hermitian energy projector", Q**2-Q)
check("mass-shell determinant multiplicity", (W(adj(energy*one+K))-m*s.eye(4)-s.Symbol('z')*s.eye(4)).det() - ((m+s.Symbol('z'))**2-m**2)**2)
A = (energy+m)*one+K
T = s.BlockMatrix([[adj(A), A], [A, -adj(A)]]).as_explicit()/(2*s.sqrt(energy*(energy+m)))
check("free FW map unitarity", T.H*T-s.eye(4))
check("free FW diagonalization", T.H*H*T-s.diag(energy,energy,-energy,-energy))
spin = s.diag(projector,projector)
rest_spin = T*spin*T.H
check("rest spin commutes with energy", H*rest_spin-rest_spin*H)
mode = T*s.Matrix.vstack(e,s.zeros(2,1))
bloch = mode[:2,:]*mode[:2,:].H+mode[2:,:]*mode[2:,:].H
rest_bloch = s.Matrix([(e.H*sigma[k]*e)[0] for k in range(3)])
lab_bloch = s.Matrix([s.trace(sigma[k]*bloch) for k in range(3)])
check("rest and rotation spin polarization",lab_bloch-(m/energy)*rest_bloch-p*p.dot(rest_bloch)/(energy*(energy+m)))

# Differential operators act on polynomials, including all product rules.
t, x, y, z, q = s.symbols("t x y z q", real=True)
r = (x, y, z)
V = x*x + y*z + t*y
potential = s.Matrix([t*x+y*z, 2*x*z+t*y, 3*x*y+x*x+t*z])
electric = -potential.diff(t)-s.Matrix([s.diff(V,a) for a in r])
magnetic = s.Matrix([
    s.diff(potential[2],y)-s.diff(potential[1],z),
    s.diff(potential[0],z)-s.diff(potential[2],x),
    s.diff(potential[1],x)-s.diff(potential[0],y),
])
f = s.Matrix([x*y+t*z+i*z*z, x*x+y*z-i*t*y])


def momentum(k, f):
    return -i*f.diff(r[k])-q*potential[k]*f


def energy_op(f):
    return i*f.diff(t)-q*V*f


def spin_momentum(f):
    return sum((sigma[k]*momentum(k,f) for k in range(3)),s.zeros(2,1))


p2 = sum((momentum(k,momentum(k,f)) for k in range(3)),s.zeros(2,1))
check("kinetic momentum squared", spin_momentum(spin_momentum(f))-p2+q*vec(magnetic)*f)
for k in range(3):
    check(f"electric commutator {k+1}", energy_op(momentum(k,f))-momentum(k,energy_op(f))-i*q*electric[k]*f)
plus = lambda f: energy_op(f)+spin_momentum(f)
minus = lambda f: energy_op(f)-spin_momentum(f)
check("ordered square electric sign 1", minus(plus(f))-energy_op(energy_op(f))+p2-q*vec(magnetic)*f-i*q*vec(electric)*f)
check("ordered square electric sign 2", plus(minus(f))-energy_op(energy_op(f))+p2-q*vec(magnetic)*f+i*q*vec(electric)*f)

# Static FW double commutator fixes spin-orbit and Darwin signs.
V = x*x+y*z
potential = potential.subs(t,0)
electric = -s.Matrix([s.diff(V,a) for a in r])
comm = lambda f: spin_momentum(q*V*f)-q*V*spin_momentum(f)
double = spin_momentum(comm(f))-comm(spin_momentum(f))
ep = [sum((electric[a]*momentum(b,f)-electric[b]*momentum(a,f) for a,b in [pair]),s.zeros(2,1)) for pair in [(1,2),(2,0),(0,1)]]
spin_orbit = sum((sigma[k]*ep[k] for k in range(3)),s.zeros(2,1))
darwin = sum(s.diff(electric[k],r[k]) for k in range(3))*f
check("FW spin-orbit and Darwin", double-q*darwin-2*q*spin_orbit)

# Uniform magnetic field: include both signs of charge and the lowest level.
for charge in (-1, 1):
    def landau_k(f):
        return -i*sigma[0]*f.diff(x)-charge*x*sigma[1]*f

    for level in (0,1):
        for spin_sign in (-1,1):
            basis = s.Matrix([1,0]) if spin_sign == 1 else s.Matrix([0,1])
            seed = s.hermite(level,x)*s.exp(-x*x/2)*basis
            squared = 2*level+1-charge*spin_sign
            label = f"q={charge}, n={level}, s={spin_sign}"
            check(f"Landau oscillator {label}",landau_k(landau_k(seed))-squared*seed)
            en = s.sqrt(4+squared)
            first = (en+2)*seed-landau_k(seed)
            second = (en+2)*seed+landau_k(seed)
            check(f"Dirac Landau mode {label}",s.Matrix.vstack(-landau_k(first)+2*second-en*first,2*first+landau_k(second)-en*second))

# SG propagator for a dispersing Gaussian: tests translation and cubic phase.
for force_sign in (-1,1):
    mass = s.Integer(2)
    force = 3*force_sign
    width = 1+i*t/mass
    shifted = z-force*t*t/(2*mass)
    packet = s.exp(i*force*t*z-i*force**2*t**3/(6*mass))*s.exp(-shifted**2/(2*width))/s.sqrt(width)
    check(f"Stern-Gerlach packet {force_sign}",i*s.diff(packet,t)+s.diff(packet,z,2)/(2*mass)+force*z*packet)


# Additional full-vector checks for the reorganized transformation derivations.
I=i; eye=one
C=I*sigma[1]
def star(z):
    return C.inv()*z.conjugate()*C
u=s.Matrix([s.Rational(2,3),s.Rational(-1,3),s.Rational(2,3)])
v=s.Matrix([s.Rational(1,3),s.Rational(2,3),s.Rational(2,3)])
R=s.Rational(3,5)*eye+I*s.Rational(4,5)*vec(u)
L=s.Rational(5,4)*eye+s.Rational(3,4)*vec(u)
gamma=s.Rational(17,8); beta=s.Rational(15,17)
S=s.symbols('S'); V=s.Matrix(s.symbols('V1:4')); Z=S*eye+vec(V)
check('rotor unitary',R.H*R-eye)
check('rotor sigma conjugation',star(R)-R)
check('rotor determinant',R.det()-1)
check('boost Hermitian',L.H-L)
check('boost sigma conjugation inverse',star(L)*L-eye)
check('boost determinant',L.det()-1)
check('boost gamma prefactor',L-s.sqrt((gamma+1)/2)*(eye+gamma*beta/(gamma+1)*vec(u)))
check('boost squared',L**2-gamma*(eye+beta*vec(u)))
check('boost norm',sc(L.H*L)-gamma)
rotated=V*s.Rational(-7,25)+u.cross(V)*s.Rational(24,25)+u*(u.dot(V))*s.Rational(32,25)
check('Rodrigues full vector',R.H*Z*R-S*eye-vec(rotated))
parallel=u*u.dot(V); perp=V-parallel
boosted=gamma*(S-beta*u.dot(V))*eye+vec(perp+gamma*(parallel-beta*S*u))
check('congruence scalar and all spatial components',L.inv()*Z*L.inv()-boosted)
similar=S*eye+vec(parallel+gamma*(perp-I*beta*u.cross(V)))
check('similarity all vector components',L.inv()*Z*L-similar)
check('Lie boost commutator',vec(u)*vec(v)/4-vec(v)*vec(u)/4-I*vec(u.cross(v))/2)
U=L**2
check('proper velocity star normalization',star(U)*U-eye)
check('proper velocity Hermitian',U.H-U)
check('proper velocity positive norm',sc(U.H*U)-(2*gamma**2-1))
assert U.H*U != eye
checks+=1; print('PASS proper velocity is not unitary away from rest')
X=7*eye+vec(s.Matrix([2,-3,5])); P=11*eye+vec(s.Matrix([1,4,-2]))
Xp=L.inv()*X*L.inv(); Pp=L.inv()*P*L.inv()
check('phase scalar pairing',sc(adj(Pp)*Xp)-sc(adj(P)*X))
psi=s.Matrix([1+2*I,3-I]); obs=2*eye+vec(s.Matrix([2,-1,3]))
psip=R.H*psi; obsp=R.H*obs*R
check('simultaneous state observable probability invariance',(psip.H*obsp*psip)[0]/(psip.H*psip)[0]-(psi.H*obs*psi)[0]/(psi.H*psi)[0])

# Spin Casimirs: left multiplication versus the three-vector representation.
spin_one = tuple(i*s.Matrix([
    [0,-axis[2],axis[1]], [axis[2],0,-axis[0]], [-axis[1],axis[0],0]
]) for axis in (s.Matrix([1,0,0]),s.Matrix([0,1,0]),s.Matrix([0,0,1])))
check("spin one-half Casimir",sum((a*a/4 for a in sigma),s.zeros(2))-3*one/4)
check("spin one Casimir",sum((a*a for a in spin_one),s.zeros(3))-2*s.eye(3))
check("spin one commutator",spin_one[0]*spin_one[1]-spin_one[1]*spin_one[0]-i*spin_one[2])
for sign in (-1,1):
    vector = s.Matrix([1,sign*i,0])
    check(f"spin one projection {sign}",spin_one[2]*vector-sign*vector)

# The four free massless modes include the two opposite chiralities.
p=s.Matrix([2,-1,2]); energy=s.Integer(3); momentum_matrix=vec(p)
H=s.diag(-momentum_matrix,momentum_matrix)
for energy_sign in (-1,1):
    for chirality in (-1,1):
        projector=(one+energy_sign*chirality*momentum_matrix/energy)/2
        column=projector*s.Matrix([1,0])
        column/=s.sqrt((column.H*column)[0])
        mode=s.Matrix.vstack(column,s.zeros(2,1)) if chirality==-1 else s.Matrix.vstack(s.zeros(2,1),column)
        check(f"massless mode {energy_sign},{chirality}",H*mode-energy_sign*energy*mode)
        check(f"massless norm {energy_sign},{chirality}",(mode.H*mode)[0]-1)

# The proposed unprojected law must pass the ordinary q'=0 magnetic limit.
tau=s.symbols('tau',real=True)
q,magnetic,mass=s.symbols('q magnetic mass',real=True,nonzero=True)
angle=q*magnetic*tau/mass
complex_momentum=mass*s.cos(angle)*one+i*mass*s.sin(angle)*sigma[2]
complex_position=mass*s.sin(angle)*one/(q*magnetic)+i*mass*(1-s.cos(angle))*sigma[2]/(q*magnetic)
check("complex left evolution",complex_momentum.diff(tau)-q*i*magnetic*sigma[2]*complex_momentum/mass)
check("complex determinant conservation",complex_momentum.det()-mass**2)
check("complex coordinate integration",complex_position.diff(tau)-complex_momentum/mass)
check("imaginary magnetic impulse at zero magnetic charge",(complex_momentum.diff(tau)-complex_momentum.diff(tau).H).subs(tau,0)/(2*i)-q*magnetic*sigma[2])

# Constant-field real evolution with both electric and magnetic charges.
field=(2+3*i)*sigma[2]; coupling=1-2*i; mass=s.Integer(4)
T=s.diag(s.exp((8-i)*tau/8),s.exp(-(8-i)*tau/8))
initial=5*one+3*sigma[0]
real_momentum=T*initial*T.H
check("dyon congruence is real",real_momentum.H-real_momentum)
check("dyon congruence determinant",real_momentum.det()-initial.det())
check("dyon congruence solves force",real_momentum.diff(tau)-(coupling*field*real_momentum+real_momentum*(coupling*field).H)/(2*mass))

# Sequential spin/energy probabilities have a different denominator after conditioning.
sequential_H = s.BlockMatrix([[-sigma[2], one], [one, sigma[2]]]).as_explicit()
energy_projection = (s.eye(4)+sequential_H/s.sqrt(2))/2
spin_projection = s.diag((one+sigma[0])/2, (one+sigma[0])/2)
initial_state = s.Matrix([1, 0, 0, 0])
energy_then_spin = spin_projection*energy_projection*initial_state
spin_then_energy = energy_projection*spin_projection*initial_state
joint_forward = (energy_then_spin.H*energy_then_spin)[0]
joint_reverse = (spin_then_energy.H*spin_then_energy)[0]
energy_probability = (initial_state.H*energy_projection*initial_state)[0]
check("sequential energy then spin joint", joint_forward-(1-1/s.sqrt(2))/4)
check("sequential spin then energy joint", joint_reverse-s.Rational(1,4))
check("conditional spin after energy", joint_forward/energy_probability-s.Rational(1,2))

# A varying frame tests the derivative terms with noncommuting matrices.
frame_parameter = s.symbols('frame_parameter', real=True)
T = s.Matrix([[1, frame_parameter], [frame_parameter, 1+frame_parameter**2]])
Z = s.Matrix([[frame_parameter**2, 1+i*frame_parameter],
              [frame_parameter**3-i, 2-frame_parameter]])
omega = T.inv()*T.diff(frame_parameter)
commutator = lambda a, b: a*b-b*a
transformed = T.inv()*Z*T
check("varying similarity first derivative",
      transformed.diff(frame_parameter)-T.inv()*Z.diff(frame_parameter)*T
      -commutator(transformed, omega))
check("varying similarity second derivative",
      transformed.diff(frame_parameter, 2)-T.inv()*Z.diff(frame_parameter, 2)*T
      -2*commutator(transformed.diff(frame_parameter), omega)
      -commutator(transformed, omega.diff(frame_parameter))
      -commutator(omega, commutator(transformed, omega)))
real_field = Z+Z.H
transformed = T.inv()*real_field*T.inv().H
check("varying congruence derivative",
      transformed.diff(frame_parameter)-T.inv()*real_field.diff(frame_parameter)*T.inv().H
      +omega*transformed+transformed*omega.H)

print(f"{checks} exact symbolic checks passed.")
