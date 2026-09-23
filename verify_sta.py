"""Independent exact matrix/operator checks (requires SymPy).

Run with ``python3 verify_sta.py``. Convention-independent algebra and quantum
checks are retained together with historical inverse-factor Lorentz checks.
Those checks use positive-exponent R and L and X'=K^-1 X (K^-1)^H; their
frame factor K is the inverse of the current manuscript's direct factor T.
The moving-frame generator there is K^-1 dK/ds.  These are equivalent
parameterizations, not the current manuscript's factor labels.  Direct-action
R=Q, L, field, spinor, and moving-frame conventions are tested independently
by ``verify_transformations.py``.
"""

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

# Complex bilinear normalization supports the spectral projectors;
# positive-norm normalization can instead destroy their idempotence.
complex_direction = s.Matrix([s.sqrt(2), i, 0])
complex_involution = vec(complex_direction)
positive_normalized = complex_involution/s.sqrt(3)
bilinear_projector = (one+complex_involution)/2
incorrect_projector = (one+positive_normalized)/2
check("complex bilinear unit vector gives involution", complex_involution**2-one)
check("bilinear unit vector need not have unit norm", sc(complex_involution.H*complex_involution)-3)
check("positive normalization changes bilinear square", positive_normalized**2-one/3)
check("positive normalization cannot replace spectral normalization",
      incorrect_projector**2-incorrect_projector+one/6)
check("bilinear spectral projector is idempotent", bilinear_projector**2-bilinear_projector)
check("bilinear spectral projector need not be Hermitian",
      bilinear_projector.H-bilinear_projector+i*sigma[1])
normalized_null = vec(s.Matrix([1, i, 0]))/s.sqrt(2)
check("nonzero bilinear null vector has positive unit norm", sc(normalized_null.H*normalized_null)-1)
check("positive unit norm permits square zero", normalized_null**2)
singular_normalized = (one+sigma[2])/s.sqrt(2)
check("singular paravector can have positive unit norm", sc(singular_normalized.H*singular_normalized)-1)
check("positive unit norm does not imply invertibility", singular_normalized.det())

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

# Coordinate-chain-rule checks use fields evaluated on the inverse event map.
# A noncollinear rotation and boost test a general constant Lorentz factor,
# rather than assuming the operator transformation that is being verified.
coordinates = s.symbols('coordinate_t coordinate_x coordinate_y coordinate_z', real=True)
primed_coordinates = s.symbols('primed_t primed_x primed_y primed_z', real=True)
coordinate_t, coordinate_x, coordinate_y, coordinate_z = coordinates
frame_rotation = s.Rational(3, 5)*one+i*s.Rational(4, 5)*sigma[0]
frame_boost = s.Rational(5, 4)*one+s.Rational(3, 4)*sigma[2]
frame_factor = frame_rotation*frame_boost
inverse_factor = frame_factor.inv()
primed_event = primed_coordinates[0]*one+vec(s.Matrix(primed_coordinates[1:]))
inverse_event = frame_factor*primed_event*frame_factor.H
inverse_components = [sc(inverse_event)]+[sc(axis*inverse_event) for axis in sigma]
coordinate_substitution = dict(zip(coordinates, inverse_components))


def in_primed_coordinates(field):
    return field.subs(coordinate_substitution, simultaneous=True).applyfunc(s.expand)


def spacetime_derivative(field, chart, spatial_sign=1):
    return field.diff(chart[0])+spatial_sign*sum(
        (sigma[k]*field.diff(chart[k+1]) for k in range(3)), s.zeros(*field.shape))


scalar_field = coordinate_t*coordinate_x+coordinate_y**2-coordinate_t*coordinate_z
scalar_matrix = scalar_field*one
check("Lorentz coordinate chain rule",
      spacetime_derivative(in_primed_coordinates(scalar_matrix), primed_coordinates)
      -frame_factor.H*in_primed_coordinates(spacetime_derivative(scalar_matrix, coordinates))*frame_factor)

coordinate_potential = (coordinate_x**2+coordinate_y*coordinate_z)*one+vec(s.Matrix([
    coordinate_t*coordinate_y+coordinate_x*coordinate_z,
    coordinate_t*coordinate_z+coordinate_y**2,
    coordinate_t*coordinate_x+coordinate_x*coordinate_y]))
primed_potential = inverse_factor*in_primed_coordinates(coordinate_potential)*inverse_factor.H


def field_from_potential(potential_field, chart):
    scalar_potential = sc(potential_field)
    vector_potential = s.Matrix([sc(axis*potential_field) for axis in sigma])
    electric_field = -vector_potential.diff(chart[0])-s.Matrix([
        s.diff(scalar_potential, entry) for entry in chart[1:]])
    magnetic_field = s.Matrix([
        s.diff(vector_potential[2], chart[2])-s.diff(vector_potential[1], chart[3]),
        s.diff(vector_potential[0], chart[3])-s.diff(vector_potential[2], chart[1]),
        s.diff(vector_potential[1], chart[1])-s.diff(vector_potential[0], chart[2])])
    return vec(electric_field+i*magnetic_field)


coordinate_field = field_from_potential(coordinate_potential, coordinates)
primed_field = field_from_potential(primed_potential, primed_coordinates)
check("field similarity from transformed potentials and coordinates",
      primed_field-inverse_factor*in_primed_coordinates(coordinate_field)*frame_factor)
coordinate_source = adj(spacetime_derivative(coordinate_field, coordinates))
primed_source = inverse_factor*in_primed_coordinates(coordinate_source)*inverse_factor.H
check("Maxwell covariance on primed field arguments",
      spacetime_derivative(primed_field, primed_coordinates)-adj(primed_source))

coordinate_spinor = s.Matrix([
    coordinate_t*coordinate_x+i*coordinate_y**2,
    coordinate_x*coordinate_z-i*coordinate_t*coordinate_y])
covariance_charge = s.Integer(2)
primed_second_spinor = inverse_factor*in_primed_coordinates(coordinate_spinor)
primed_first_spinor = frame_factor.H*in_primed_coordinates(coordinate_spinor)
check("first Dirac operator covariance on primed field arguments",
      i*spacetime_derivative(primed_second_spinor, primed_coordinates)
      -covariance_charge*adj(primed_potential)*primed_second_spinor
      -frame_factor.H*in_primed_coordinates(
          i*spacetime_derivative(coordinate_spinor, coordinates)
          -covariance_charge*adj(coordinate_potential)*coordinate_spinor))
check("second Dirac operator covariance on primed field arguments",
      i*spacetime_derivative(primed_first_spinor, primed_coordinates, -1)
      -covariance_charge*primed_potential*primed_first_spinor
      -inverse_factor*in_primed_coordinates(
          i*spacetime_derivative(coordinate_spinor, coordinates, -1)
          -covariance_charge*coordinate_potential*coordinate_spinor))

# Integrate a fixed-axis, constant-proper-acceleration worldline and check
# its tangent, invariant acceleration, coordinate speed, and hyperbola.
proper_acceleration = s.symbols('proper_acceleration', positive=True)
proper_s = s.symbols('proper_s', real=True)
worldline_axis = s.Matrix([s.Rational(2, 3), s.Rational(-1, 3), s.Rational(2, 3)])
worldline = (s.sinh(proper_acceleration*proper_s)*one
             +(s.cosh(proper_acceleration*proper_s)-1)*vec(worldline_axis))/proper_acceleration
worldline_velocity = (s.cosh(proper_acceleration*proper_s)*one
                      +s.sinh(proper_acceleration*proper_s)*vec(worldline_axis))
check("uniform proper acceleration worldline tangent", worldline.diff(proper_s)-worldline_velocity)
check("uniform acceleration tangent normalization", worldline_velocity.det()-1)
check("uniform acceleration invariant magnitude",
      worldline_velocity.diff(proper_s).det()+proper_acceleration**2)
worldline_position = s.Matrix([sc(axis*worldline) for axis in sigma])
check("uniform acceleration hyperbolic trajectory",
      (worldline_axis.dot(worldline_position)+1/proper_acceleration)**2
      -sc(worldline)**2-1/proper_acceleration**2)
check("uniform acceleration coordinate speed",
      worldline_position.diff(proper_s)/sc(worldline).diff(proper_s)
      -s.tanh(proper_acceleration*proper_s)*worldline_axis)

# Coordinate velocity is a ratio of transformed differential components.
coordinate_velocity = s.Matrix(s.symbols('velocity_1:4', real=True))
worldline_boost = s.Rational(5, 4)*one+s.Rational(3, 4)*vec(worldline_axis)
worldline_gamma = s.Rational(17, 8)
worldline_beta = s.Rational(15, 17)
transformed_increment = worldline_boost.inv()*(one+vec(coordinate_velocity))*worldline_boost.inv()
velocity_parallel = worldline_axis*worldline_axis.dot(coordinate_velocity)
expected_velocity = ((coordinate_velocity-velocity_parallel)/worldline_gamma
                     +velocity_parallel-worldline_beta*worldline_axis)/(
                         1-worldline_beta*worldline_axis.dot(coordinate_velocity))
transformed_velocity = s.Matrix([sc(axis*transformed_increment) for axis in sigma])/sc(transformed_increment)
check("velocity transformation from spacetime increments", transformed_velocity-expected_velocity)
check("velocity transformation preserves the light-speed bound",
      1-expected_velocity.dot(expected_velocity)
      -(1-coordinate_velocity.dot(coordinate_velocity))/(
          worldline_gamma**2*(1-worldline_beta*worldline_axis.dot(coordinate_velocity))**2))

momentum_energy = s.symbols('momentum_energy', real=True)
momentum_vector = s.Matrix(s.symbols('momentum_1:4', real=True))
momentum_parallel = worldline_axis*worldline_axis.dot(momentum_vector)
boosted_momentum = worldline_boost.inv()*(momentum_energy*one+vec(momentum_vector))*worldline_boost.inv()
check("energy and momentum component boost",
      boosted_momentum-worldline_gamma*(momentum_energy-worldline_beta*worldline_axis.dot(momentum_vector))*one
      -vec(momentum_vector-momentum_parallel+worldline_gamma*(
          momentum_parallel-worldline_beta*momentum_energy*worldline_axis)))

# Differentiate a trajectory with acceleration both along and across motion.
trajectory_time = s.symbols('trajectory_time', real=True)
trajectory_mass = s.symbols('trajectory_mass', positive=True)
trajectory_acceleration = s.Matrix(s.symbols('acceleration_1:4', real=True))
initial_velocity = s.Matrix([s.Rational(1, 3), s.Rational(2, 3), s.Rational(1, 3)])
trajectory_velocity = initial_velocity+trajectory_time*trajectory_acceleration
trajectory_gamma = 1/s.sqrt(1-trajectory_velocity.dot(trajectory_velocity))
initial_gamma = trajectory_gamma.subs(trajectory_time, 0)
parallel_acceleration = initial_velocity*initial_velocity.dot(trajectory_acceleration)/initial_velocity.dot(initial_velocity)
coordinate_force = (trajectory_mass*trajectory_gamma*trajectory_velocity).diff(trajectory_time).subs(trajectory_time, 0)
check("longitudinal and transverse relativistic inertia",
      coordinate_force-trajectory_mass*initial_gamma*(trajectory_acceleration-parallel_acceleration)
      -trajectory_mass*initial_gamma**3*parallel_acceleration)
check("relativistic work equals mechanical energy derivative",
      (trajectory_mass*trajectory_gamma).diff(trajectory_time).subs(trajectory_time, 0)
      -initial_velocity.dot(coordinate_force))

# Append before verify_sta.py's final print. Existing helpers: s, i, check.
# A non-Cartesian orthonormal frame exercises every spatial component.
sg_axis = s.Matrix([s.Rational(2, 3), s.Rational(1, 3), s.Rational(2, 3)])
sg_transverse = s.Matrix([-s.Rational(2, 3), s.Rational(2, 3), s.Rational(1, 3)])
sg_coordinates = s.Matrix(s.symbols('sg_x sg_y sg_z', real=True))
sg_wavevector = s.Matrix(s.symbols('sg_k1 sg_k2 sg_k3', real=True))
sg_time, sg_bias, sg_gradient, sg_moment = s.symbols(
    'sg_time sg_bias sg_gradient sg_moment', real=True)
sg_mass = s.symbols('sg_mass', positive=True)
sg_field = ((sg_bias+sg_gradient*sg_axis.dot(sg_coordinates))*sg_axis
            -sg_gradient*sg_transverse.dot(sg_coordinates)*sg_transverse)
check('general-axis Stern-Gerlach vacuum divergence',
      sum(s.diff(sg_field[k], sg_coordinates[k]) for k in range(3)))
check('general-axis Stern-Gerlach vacuum curl', s.Matrix([
    s.diff(sg_field[2], sg_coordinates[1])-s.diff(sg_field[1], sg_coordinates[2]),
    s.diff(sg_field[0], sg_coordinates[2])-s.diff(sg_field[2], sg_coordinates[0]),
    s.diff(sg_field[1], sg_coordinates[0])-s.diff(sg_field[0], sg_coordinates[1])]))

# Every free Fourier mode obeys the claimed packet translation and phase.
# Their arbitrary superpositions therefore obey the same evolution law.
for sg_sign in (-1, 1):
    sg_shift = (sg_coordinates
                -sg_sign*sg_moment*sg_gradient*sg_time**2*sg_axis/(2*sg_mass))
    sg_phase = (sg_sign*sg_moment*sg_bias*sg_time
                +sg_sign*sg_moment*sg_gradient*sg_time*sg_axis.dot(sg_coordinates)
                -sg_moment**2*sg_gradient**2*sg_time**3/(6*sg_mass)
                +sg_wavevector.dot(sg_shift)
                -sg_wavevector.dot(sg_wavevector)*sg_time/(2*sg_mass))
    # Divide (i d_t-H) exp(i phase) by exp(i phase), retaining all terms.
    sg_residual = (-s.diff(sg_phase, sg_time)
                   -sum(s.diff(sg_phase, x)**2-i*s.diff(sg_phase, x, 2)
                        for x in sg_coordinates)/(2*sg_mass)
                   +sg_sign*sg_moment*(sg_bias+sg_gradient*sg_axis.dot(sg_coordinates)))
    check(f'general-axis Stern-Gerlach packet branch {sg_sign}', sg_residual)

sg_polarization = (sg_transverse*s.cos(2*sg_moment*sg_bias*sg_time)
                   -sg_axis.cross(sg_transverse)*s.sin(2*sg_moment*sg_bias*sg_time))
check('general-axis spin-precession sign',
      sg_polarization.diff(sg_time)-2*sg_moment*sg_polarization.cross(sg_bias*sg_axis))

# Spin three-halves from three constrained ideal entries, using ordinary blocks.
# The matrix realization identifies a normalized ideal with its normalized column.
def check_paravector_spin_three_halves():
    axes = tuple(s.eye(3)[:, k] for k in range(3))
    axis = s.Matrix([s.Rational(2, 3), s.Rational(1, 3), s.Rational(2, 3)])
    transverse = s.Matrix([-s.Rational(2, 3), s.Rational(2, 3), s.Rational(1, 3)])
    third = axis.cross(transverse)
    plus = s.Matrix([s.sqrt(s.Rational(5, 6)), (2+i)/s.sqrt(30)])
    spinors = {1: plus, -1: vec(transverse)*plus}
    circular = {sign: (transverse+i*sign*third)/s.sqrt(2) for sign in (-1, 1)}
    states = {}
    for sign in (-1, 1):
        states[sign*s.Rational(3, 2)] = s.Matrix.vstack(*[
            sign*circular[sign][k]*spinors[sign] for k in range(3)
        ])
        states[sign*s.Rational(1, 2)] = s.Matrix.vstack(*[
            (sign*circular[sign][k]*spinors[-sign]
             -s.sqrt(2)*axis[k]*spinors[sign])/s.sqrt(3)
            for k in range(3)
        ])
    projections = (s.Rational(3, 2), s.Rational(1, 2),
                   -s.Rational(1, 2), -s.Rational(3, 2))
    basis = s.simplify(s.Matrix.hstack(*(states[value] for value in projections)))
    constraint = s.Matrix.hstack(*sigma)

    def generator(direction):
        return s.BlockMatrix([
            [i*direction.cross(axes[j])[k]*one
             +(vec(direction)/2 if k == j else s.zeros(2))
             for j in range(3)] for k in range(3)
        ]).as_explicit()

    generators = tuple(generator(direction) for direction in axes)
    squared_spin = sum((entry*entry for entry in generators), s.zeros(6))
    projector = s.eye(6)-constraint.H*constraint/3
    check('three-entry spin constraint is two independent conditions',
          constraint*constraint.H-3*one)
    check('spin three-halves arbitrary-axis basis constraint', constraint*basis)
    check('spin three-halves arbitrary-axis basis orthonormality', basis.H*basis-s.eye(4))
    check('spin three-halves full six-component squared spin',
          squared_spin-s.Rational(15, 4)*s.eye(6)+constraint.H*constraint)
    check('spin three-halves orthogonal constraint projector', s.Matrix.vstack(
        projector.H-projector, projector*projector-projector, constraint*projector))
    check('spin three-halves constraint projector rank four', s.trace(projector)-4)
    check('spin three-halves arbitrary-axis basis completeness', basis*basis.H-projector)
    check('excluded spin one-half block',
          (squared_spin-s.Rational(3, 4)*s.eye(6))*(s.eye(6)-projector))

    lowering = s.Matrix([[0, 0, 0, 0], [s.sqrt(3), 0, 0, 0],
                         [0, 2, 0, 0], [0, 0, s.sqrt(3), 0]])
    standard = ((lowering+lowering.H)/2,
                (lowering.H-lowering)/(2*i), s.diag(*projections))
    local_generators = tuple(generator(direction) for direction in (transverse, third, axis))
    check('spin three-halves standard generator matrices on an arbitrary axis', s.Matrix.vstack(*[
        s.simplify(basis.H*entry*basis)-expected
        for entry, expected in zip(local_generators, standard)
    ]))
    check('spin three-halves invariant basis and lowering coefficients', s.Matrix.vstack(*[
        entry*basis-basis*expected for entry, expected in zip(local_generators, standard)
    ]))
    check('spin three-halves rotation commutators', s.Matrix.vstack(*[
        generators[a]*generators[b]-generators[b]*generators[a]-i*generators[c]
        for a, b, c in ((0, 1, 2), (1, 2, 0), (2, 0, 1))
    ]))
    check('spin three-halves squared spin on the constrained block',
          squared_spin*basis-s.Rational(15, 4)*basis)

    # A non-coordinate rotor, independent of the measurement axis.
    rotation_axis = s.Matrix([s.Rational(1, 3), s.Rational(2, 3), s.Rational(2, 3)])
    rotor = s.Rational(3, 5)*one+s.Rational(4, 5)*i*vec(rotation_axis)
    rotation = s.Matrix(3, 3, lambda k, j: s.simplify(sc(sigma[k]*rotor.H*sigma[j]*rotor)))
    transform = s.BlockMatrix([
        [rotation[k, j]*rotor.H for j in range(3)] for k in range(3)
    ]).as_explicit()
    check('spin three-halves finite rotation preserves norm', transform.H*transform-s.eye(6))
    check('spin three-halves finite constraint covariance',
          constraint*transform-rotor.H*constraint)
    check('spin three-halves finite rotation preserves the constrained states',
          transform*projector-projector*transform)
    reduced_rotation = (basis.H*transform*basis).applyfunc(s.expand)
    check('spin three-halves finite four-component rotation is unitary',
          reduced_rotation.H*reduced_rotation-s.eye(4))
    check('spin three-halves full-turn sign',
          (-2*s.pi*i*standard[2]).exp()+s.eye(4))
    check('spin three-halves double-turn restoration',
          (-4*s.pi*i*standard[2]).exp()-s.eye(4))


check_paravector_spin_three_halves()

# Append before verify_sta.py's final print. Uses its s, sigma, vec, check.
# Local symbols prevent interactions with the other exact-check sections.
def verify_low_energy_recursions():
    time, x, y, z, charge = s.symbols(
        "low_time low_x low_y low_z low_charge", real=True
    )
    mass = s.symbols("low_mass", positive=True)
    coordinates = (x, y, z)
    potential = time*x+y*y
    vector_potential = s.Matrix([time*z, x*z, x*y+time*y])
    electric = -vector_potential.diff(time)-s.Matrix(
        [s.diff(potential, coordinate) for coordinate in coordinates]
    )
    magnetic = s.Matrix([
        s.diff(vector_potential[2], y)-s.diff(vector_potential[1], z),
        s.diff(vector_potential[0], z)-s.diff(vector_potential[2], x),
        s.diff(vector_potential[1], x)-s.diff(vector_potential[0], y),
    ])
    amplitude = s.Matrix([x*y+time*z+s.I*z*z, x*x+y*z-s.I*time*y])

    def energy(field):
        return s.I*field.diff(time)-charge*potential*field

    def momentum(axis, field):
        return -s.I*field.diff(coordinates[axis])-charge*vector_potential[axis]*field

    def sigma_momentum(field):
        return sum(
            (sigma[axis]*momentum(axis, field) for axis in range(3)),
            s.zeros(2, 1),
        )

    def momentum_square(field):
        return sum(
            (momentum(axis, momentum(axis, field)) for axis in range(3)),
            s.zeros(2, 1),
        )

    def coupled_plus(field):
        return energy(field)+mass*field+sigma_momentum(field)

    def coupled_minus(field):
        return energy(field)+mass*field-sigma_momentum(field)

    ordered_square = coupled_plus(coupled_minus(amplitude))-mass**2*amplitude
    squared_equation = (
        energy(energy(amplitude))+2*mass*energy(amplitude)
        -momentum_square(amplitude)+charge*vec(magnetic)*amplitude
        -s.I*charge*vec(electric)*amplitude
    )
    momentum_energy_commutator = (
        sigma_momentum(energy(amplitude))-energy(sigma_momentum(amplitude))
    )
    residual_recursion = energy(amplitude)-(
        momentum_square(amplitude)-energy(energy(amplitude))
        -charge*vec(magnetic)*amplitude-momentum_energy_commutator
    )/(2*mass)
    check(
        "Dirac residual-energy recurrence with time-dependent fields",
        s.Matrix.hstack(
            ordered_square-squared_equation,
            2*mass*residual_recursion-squared_equation,
        ),
    )

    # Noncommuting Hermitian symbols represent p.sigma and qV. Conjugating
    # the unnormalized Hamiltonian must give the double commutator.
    sigma_operator, potential_operator = s.symbols(
        "low_sigma_operator low_potential_operator", commutative=False,
        hermitian=True,
    )
    inverse_mass = s.symbols("low_inverse_mass", real=True)

    def commutator(left, right):
        return left*right-right*left

    envelope_hamiltonian = (
        potential_operator+inverse_mass*sigma_operator**2/2
        -inverse_mass**2*sigma_operator*commutator(sigma_operator, potential_operator)/4
        -inverse_mass**3*sigma_operator**4/8
    )
    norm_change = 1+inverse_mass**2*sigma_operator**2/8
    inverse_norm_change = 1-inverse_mass**2*sigma_operator**2/8
    physical_hamiltonian = (
        potential_operator+inverse_mass*sigma_operator**2/2
        -inverse_mass**2*commutator(
            sigma_operator, commutator(sigma_operator, potential_operator)
        )/8
        -inverse_mass**3*sigma_operator**4/8
    )
    difference = s.expand(
        norm_change*envelope_hamiltonian*inverse_norm_change-physical_hamiltonian
    )
    check(
        "Dirac envelope normalization produces Hermitian double commutator",
        s.Matrix([
            s.expand(difference.coeff(inverse_mass, power)) for power in range(4)
        ]+[
            s.expand(s.adjoint(physical_hamiltonian)-physical_hamiltonian)
        ]),
    )

    momentum_squared = s.symbols("low_momentum_squared", nonnegative=True)
    previous_energy = momentum_squared/(2*mass)-momentum_squared**2/(8*mass**3)
    next_energy = previous_energy+momentum_squared**3/(16*mass**5)
    iterated_energy = (momentum_squared-previous_energy**2)/(2*mass)
    exact_energy = s.sqrt(mass**2+momentum_squared)-mass
    check(
        "Free residual-energy recursion yields the positive p6 coefficient",
        s.Matrix([
            next_energy-s.series(iterated_energy, momentum_squared, 0, 4).removeO(),
            next_energy-s.series(exact_energy, momentum_squared, 0, 4).removeO(),
        ]),
    )


verify_low_energy_recursions()

# Polarizations on a non-Cartesian axis test the vector and helicity laws.
def check_general_axis_polarizations():
    axis = s.Matrix([s.Rational(2, 3), s.Rational(1, 3), s.Rational(2, 3)])
    transverse = s.Matrix([-s.Rational(2, 3), s.Rational(2, 3), s.Rational(1, 3)])
    plus = (transverse+i*axis.cross(transverse))/s.sqrt(2)
    minus = (transverse-i*axis.cross(transverse))/s.sqrt(2)
    basis = s.Matrix.hstack(plus, axis, minus)
    generator = i*s.Matrix([[0, -axis[2], axis[1]],
                           [axis[2], 0, -axis[0]],
                           [-axis[1], axis[0], 0]])
    check('general-axis spin-one orthonormal polarizations', basis.H*basis-s.eye(3))
    check('general-axis spin-one projections', generator*basis-basis*s.diag(1, 0, -1))
    check('circular polarizations are bilinearly null', s.Matrix([plus.dot(plus), minus.dot(minus)]))
    coefficients = s.Matrix([2+i, -1+3*i, 4-2*i])
    amplitude = basis*coefficients
    check('general-axis spin-one Born normalization',
          (amplitude.H*amplitude)[0]-(coefficients.H*coefficients)[0])
    # Maxwell's equations for C=E+iB and its conjugate sector have opposite signs.
    check('positive-frequency self-dual Maxwell helicity', i*axis.cross(plus)-plus)
    check('positive-frequency conjugate Maxwell helicity', -i*axis.cross(minus)-minus)


check_general_axis_polarizations()

# Append before verify_sta.py's final print. Uses its existing exact helpers.
def check_arbitrary_reference_projector():
    """No Cartesian projector: independent oblique reference and measurement axes."""
    ref_epsilon = -i * sigma[1]

    def ref_star(z):
        return (ref_epsilon * z.conjugate() * (-ref_epsilon)).applyfunc(s.expand)

    def ref_expand(z):
        return z.applyfunc(s.expand)

    # All three components of the reference direction are nonzero.
    ref_v = s.Matrix([1 + 2*i, 3 - i]) / s.sqrt(15)
    ref_pi = ref_expand(ref_v * ref_v.H)
    ref_c = s.Matrix(s.symbols('reference_c1 reference_c2', complex=True))
    ref_d = s.Matrix(s.symbols('reference_d1 reference_d2', complex=True))
    ref_psi = s.sqrt(2) * ref_c * ref_v.H
    ref_phi = s.sqrt(2) * ref_d * ref_v.H
    ref_rho = (ref_c.H * ref_c)[0]
    check('oblique reference projector defining identities', ref_expand(s.Matrix.vstack(
        ref_pi**2-ref_pi, ref_pi.H-ref_pi, ref_star(ref_pi)+ref_pi-one,
        s.Matrix([[sc(ref_pi)-s.Rational(1, 2), ref_pi.det()]]))))
    check('arbitrary-reference ideal and two-component pairing', ref_expand(s.Matrix.vstack(
        ref_psi*ref_pi-ref_psi,
        s.Matrix([[sc(ref_phi.H*ref_psi)-(ref_d.H*ref_c)[0],
                   sc(ref_psi.H*ref_psi)-ref_rho]]))))
    check('arbitrary-reference right norm is proportional to its projector',
          ref_expand(ref_psi.H*ref_psi-2*ref_rho*ref_pi))
    ref_state = ref_psi * ref_psi.H / ref_rho
    check('arbitrary-reference physical state loses the reference projector',
          ref_expand(ref_state-2*ref_c*ref_c.H/ref_rho))

    # Measurement axis differs from both the right reference and re-encoding axis.
    ref_u = s.Matrix([2, 1, 2])/3
    ref_measurement = (one+vec(ref_u))/2
    ref_probability = (ref_c.H*ref_measurement*ref_c)[0]/ref_rho
    check('arbitrary-reference Born probability in all three forms', s.Matrix([
        s.expand(sc((ref_measurement*ref_psi).H*(ref_measurement*ref_psi))/ref_rho
                 -ref_probability),
        s.expand(sc(ref_measurement*ref_state)-ref_probability)]))

    # Include a central phase: the freedom is full unitarity, not det(T)=1.
    ref_rotor = s.Rational(3, 5)*one+4*i*vec(s.Matrix([2, -2, 1])/3)/5
    ref_encoding = (5+12*i)*ref_rotor/13
    ref_pi_prime = ref_expand(ref_encoding.H*ref_pi*ref_encoding)
    ref_psi_prime = ref_expand(ref_psi*ref_encoding)
    ref_phi_prime = ref_expand(ref_phi*ref_encoding)
    check('unitary right re-encoding preserves the ideal', ref_expand(s.Matrix.vstack(
        ref_encoding.H*ref_encoding-one,
        ref_pi_prime**2-ref_pi_prime,
        ref_psi_prime*ref_pi_prime-ref_psi_prime)))
    check('right re-encoding preserves physical state and pairing', ref_expand(s.Matrix.vstack(
        ref_psi_prime*ref_psi_prime.H-ref_psi*ref_psi.H,
        s.Matrix([[sc(ref_phi_prime.H*ref_psi_prime)-sc(ref_phi.H*ref_psi),
                   sc((ref_measurement*ref_psi_prime).H*(ref_measurement*ref_psi_prime))
                   -sc((ref_measurement*ref_psi).H*(ref_measurement*ref_psi))]]))))

    # Construct a real unitary T without selecting any Cartesian vector.
    # Its relative phase also makes the conventional column charge map exact.
    ref_v_prime = ref_epsilon*ref_v.conjugate()
    ref_charge = ref_expand(ref_v_prime*ref_v.H+ref_v*ref_v_prime.H)
    check('oblique charge factor is real unitary and scalar-free', ref_expand(s.Matrix.vstack(
        ref_charge.H-ref_charge, ref_charge**2-one,
        s.Matrix([[sc(ref_charge), 0]]))))
    check('charge factor connects the complementary right projections', ref_expand(s.Matrix.vstack(
        ref_star(ref_pi)*ref_charge-ref_charge*ref_pi,
        ref_star(ref_charge)+ref_charge)))

    def ref_charge_block(first, second):
        return (ref_expand(ref_star(second)*ref_charge),
                ref_expand(-ref_star(first)*ref_charge))

    ref_first, ref_second = ref_charge_block(ref_psi, ref_phi)
    ref_twice_first, ref_twice_second = ref_charge_block(ref_first, ref_second)
    check('arbitrary-reference Dirac charge conjugation squares to one',
          ref_expand(s.Matrix.vstack(ref_twice_first-ref_psi, ref_twice_second-ref_phi)))
    check('arbitrary-reference charge conjugation preserves norm and ideal',
          ref_expand(s.Matrix.vstack(ref_first*ref_pi-ref_first, ref_second*ref_pi-ref_second,
          s.Matrix([[sc(ref_first.H*ref_first+ref_second.H*ref_second)
                     -sc(ref_psi.H*ref_psi+ref_phi.H*ref_phi), 0]]))))
    check('arbitrary-reference normalized matrix-column extraction',
          ref_expand(ref_psi*ref_v/s.sqrt(2)-ref_c))
    check('arbitrary-reference charge phase gives the standard column map', ref_expand(s.Matrix.vstack(
        ref_charge*ref_v-ref_epsilon*ref_v.conjugate(),
        ref_star(ref_psi)*ref_charge*ref_v/s.sqrt(2)-ref_epsilon*ref_c.conjugate())))

    # Canonical momentum reverses when the spacetime phase is conjugated.
    # This polynomial identity also checks every real potential component.
    ref_q, ref_V, ref_m = s.symbols('reference_q reference_V reference_m', real=True)
    ref_k = s.Matrix(s.symbols('reference_k1:4', real=True))
    ref_A = s.Matrix(s.symbols('reference_A1:4', real=True))
    ref_p = vec(ref_k-ref_q*ref_A)
    ref_H = s.BlockMatrix([[ref_q*ref_V*one-ref_p, ref_m*one],
                           [ref_m*one, ref_q*ref_V*one+ref_p]]).as_explicit()
    ref_H_prime = s.BlockMatrix([[-ref_q*ref_V*one+ref_p, ref_m*one],
                                 [ref_m*one, -ref_q*ref_V*one-ref_p]]).as_explicit()
    ref_C = s.BlockMatrix([[s.zeros(2), ref_epsilon],
                           [-ref_epsilon, s.zeros(2)]]).as_explicit()
    check('charge reversal intertwines opposite charge and conjugated spacetime phase',
          ref_expand(ref_H_prime*ref_C+ref_C*ref_H.conjugate()))


check_arbitrary_reference_projector()

# General-frame component formulas, checked in the paravector matrix image.
def check_arbitrary_axis_components():
    axis = s.Matrix([2, 3, 6])/7
    transverse = s.Matrix([3, -2, 0])/s.sqrt(13)
    across = axis.cross(transverse)
    reference = (one+vec(across))/2
    plus = (one+vec(axis))*reference
    minus = vec(transverse)*plus
    states = (plus, minus)
    check('arbitrary-frame ideal basis normalization and spin', s.Matrix.vstack(*[
        (entry*reference-entry).applyfunc(s.expand)
        for entry in states
    ], *[
        (vec(axis)*entry-sign*entry).applyfunc(s.expand)
        for sign,entry in ((1,plus),(-1,minus))
    ], s.Matrix(2,2,lambda j,k: s.expand(sc(states[j].H*states[k])-(1 if j==k else 0)))))
    check('general-axis ladder paravectors and phase conventions', s.Matrix.vstack(
        (vec(transverse)*minus-plus).applyfunc(s.expand),
        (vec(across)*plus-i*minus).applyfunc(s.expand),
        (vec(across)*minus+i*plus).applyfunc(s.expand)))

    coordinates = s.Matrix(s.symbols('frame_x frame_y frame_z',real=True))
    field = s.symbols('frame_B',real=True)
    potential = field*transverse.dot(coordinates)*across
    curl = s.Matrix([
        s.diff(potential[2],coordinates[1])-s.diff(potential[1],coordinates[2]),
        s.diff(potential[0],coordinates[2])-s.diff(potential[2],coordinates[0]),
        s.diff(potential[1],coordinates[0])-s.diff(potential[0],coordinates[1])])
    check('arbitrary-frame magnetic gauge curl', (curl-field*axis).applyfunc(s.expand))

    h, dh, longitudinal, charge_field, shifted = s.symbols(
        'frame_h frame_dh frame_k frame_qB frame_shift',real=True)
    magnitude = s.symbols('frame_d',positive=True)
    residuals = []
    for sign,entry,other in ((1,plus,minus),(-1,minus,plus)):
        direct = (-i*s.sqrt(magnitude)*dh*vec(transverse)
                  -charge_field*shifted*h*vec(across)/s.sqrt(magnitude)
                  +longitudinal*h*vec(axis))*entry
        coefficients = (sign*longitudinal*h*entry
                        -i*s.sqrt(magnitude)*(dh+sign*charge_field*shifted*h/magnitude)*other)
        residuals.append((direct-coefficients).applyfunc(s.expand))
    check('both arbitrary-axis magnetic Dirac component formulas', s.Matrix.vstack(*residuals))

    scalar = s.symbols('frame_scalar',complex=True)
    vector = s.Matrix(s.symbols('frame_vector1:4',complex=True))
    element = scalar*one+vec(vector)
    coefficient_matrix = s.Matrix(2,2,lambda j,k: sc(states[j].H*element*states[k]))
    expected = s.Matrix([[scalar+axis.dot(vector),(transverse-i*across).dot(vector)],
                         [(transverse+i*across).dot(vector),scalar-axis.dot(vector)]])
    check('arbitrary-axis matrix image for every complex paravector',
          (coefficient_matrix-expected).applyfunc(s.expand))


check_arbitrary_axis_components()


def check_scalar_phase_extraction():
    """Product-rule signs with nonconstant phases, potentials and amplitudes."""
    phase_t, phase_x, phase_y, phase_z = s.symbols(
        'phase_t phase_x phase_y phase_z', real=True)
    phase_q = s.symbols('phase_q', real=True)
    phase_m = s.symbols('phase_m', positive=True)
    phase_coordinates = (phase_x, phase_y, phase_z)
    phase_V = phase_t*phase_x+phase_y*phase_z
    phase_A = s.Matrix([phase_t*phase_y, phase_x*phase_z,
                        phase_t**2+phase_x*phase_y])
    phase_f = phase_t**2*phase_x+phase_t*phase_y*phase_z+phase_x**2*phase_z
    phase_factor = s.exp(i*phase_f)
    phase_amplitude = s.Matrix([
        [phase_t*phase_x+i*phase_y, phase_z**2-i*phase_t*phase_y],
        [phase_x*phase_y+i*phase_t**2, phase_t*phase_z-i*phase_x]])

    def phase_energy(amplitude, potential=phase_V):
        return i*amplitude.diff(phase_t)-phase_q*potential*amplitude

    def phase_momentum(index, amplitude, potential=phase_A):
        return (-i*amplitude.diff(phase_coordinates[index])
                -phase_q*potential[index]*amplitude)

    def phase_expand(amplitude):
        return amplitude.applyfunc(s.expand)

    phase_rate = s.diff(phase_f, phase_t)
    phased = phase_factor*phase_amplitude
    check('scalar phase shifts mechanical energy by minus its time derivative',
          phase_expand(phase_energy(phased)/phase_factor
                       -phase_energy(phase_amplitude)+phase_rate*phase_amplitude))
    check('scalar phase shifts all mechanical momenta by its gradient',
          phase_expand(s.Matrix.vstack(*[
              phase_momentum(k, phased)/phase_factor
              -phase_momentum(k, phase_amplitude)
              -s.diff(phase_f, phase_coordinates[k])*phase_amplitude
              for k in range(3)])))
    check('ordered energy square retains the second phase derivative',
          phase_expand(phase_energy(phase_energy(phased))/phase_factor
                       -phase_energy(phase_energy(phase_amplitude))
                       +2*phase_rate*phase_energy(phase_amplitude)
                       -(phase_rate**2-i*s.diff(phase_rate, phase_t))*phase_amplitude))
    check('scalar phase preserves paravector norms pointwise',
          s.expand(sc(phased.H*phased)-sc(phase_amplitude.H*phase_amplitude)))

    # The gauge function is genuinely spacetime dependent.  The same phase
    # derivative that shifts the amplitude cancels the potential change.
    gauge_factor = s.exp(-i*phase_q*phase_f)
    gauge_V = phase_V+s.diff(phase_f, phase_t)
    gauge_A = phase_A-s.Matrix([s.diff(phase_f, coordinate)
                               for coordinate in phase_coordinates])
    gauge_amplitude = gauge_factor*phase_amplitude
    check('local gauge phase cancels both transformed mechanical operators',
          phase_expand(s.Matrix.vstack(
              phase_energy(gauge_amplitude, gauge_V)/gauge_factor
              -phase_energy(phase_amplitude), *[
                  phase_momentum(k, gauge_amplitude, gauge_A)/gauge_factor
                  -phase_momentum(k, phase_amplitude) for k in range(3)])))

    rest_residuals = []
    for sign in (-1, 1):
        rest_factor = s.exp(-i*sign*phase_m*phase_t)
        rest_amplitude = rest_factor*phase_amplitude
        rest_residuals.extend([
            phase_energy(rest_amplitude)/rest_factor
            -phase_energy(phase_amplitude)-sign*phase_m*phase_amplitude,
            phase_energy(phase_energy(rest_amplitude))/rest_factor
            -phase_energy(phase_energy(phase_amplitude))
            -2*sign*phase_m*phase_energy(phase_amplitude)
            -phase_m**2*phase_amplitude])
        rest_residuals.extend([
            phase_momentum(k, rest_amplitude)/rest_factor
            -phase_momentum(k, phase_amplitude) for k in range(3)])
    check('both rest phase signs shift energy and its square but not momentum',
          phase_expand(s.Matrix.vstack(*rest_residuals)))

    scalar_envelope = phase_amplitude[0, 0]
    scalar_full = s.exp(-i*phase_m*phase_t)*scalar_envelope

    def phase_density(amplitude, potential=phase_V):
        operated = phase_energy(amplitude, potential)
        return (s.conjugate(amplitude)*operated
                +s.conjugate(operated)*amplitude)/(2*phase_m)

    check('KG rest phase density retains the residual energy contribution',
          s.expand(phase_density(scalar_full)
                   -s.conjugate(scalar_envelope)*scalar_envelope
                   -phase_density(scalar_envelope)))
    check('KG density is invariant under the simultaneous gauge transformation',
          s.expand(phase_density(gauge_factor*scalar_envelope, gauge_V)
                   -phase_density(scalar_envelope)))


def check_arbitrary_ideal_coefficient_constraints():
    """The ideal constraints use the projector vector, never a chosen axis."""
    ideal_direction = s.Matrix([2, 3, 6])/7
    ideal_projector_vector = ideal_direction/2
    ideal_projector = one/2+vec(ideal_projector_vector)
    ideal_scalar = s.symbols('ideal_scalar', complex=True)
    ideal_vector = s.Matrix(s.symbols('ideal_vector1:4', complex=True))
    ideal_element = ideal_scalar*one+vec(ideal_vector)
    ideal_amplitude = (ideal_element*ideal_projector).applyfunc(s.expand)
    ideal_a = s.expand(sc(ideal_amplitude))
    ideal_b = s.Matrix([s.expand(sc(sigma[k]*ideal_amplitude))
                       for k in range(3)])
    check('arbitrary ideal amplitudes obey scalar and vector coefficient constraints',
          s.Matrix.vstack(
              s.Matrix([2*ideal_b.dot(ideal_projector_vector)-ideal_a]),
              2*ideal_a*ideal_projector_vector
              +2*i*ideal_b.cross(ideal_projector_vector)-ideal_b).applyfunc(s.expand))
    # The entire matrix residual is exactly half of these two constraints,
    # which verifies sufficiency as well as necessity.
    check('the two coefficient constraints exactly characterize the right ideal',
          (2*(ideal_element*ideal_projector-ideal_element)
           -(2*ideal_vector.dot(ideal_projector_vector)-ideal_scalar)*one
           -vec(2*ideal_scalar*ideal_projector_vector
                +2*i*ideal_vector.cross(ideal_projector_vector)-ideal_vector)
           ).applyfunc(s.expand))


check_scalar_phase_extraction()
check_arbitrary_ideal_coefficient_constraints()

print(f"{checks} exact symbolic checks passed.")
