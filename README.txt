SPACE TIME ALGEBRA (STA)
σₙ² = 1
σₙσₘ = –σₘσₙ  if  n ≠ m
i = σ₁σ₂σ₃                        ←  i² = –1
σ = {σ₁, σ₂, σ₃}
i σ = {σ₂σ₃, σ₃σ₁, σ₁σ₂}
( A • σ ) ( B • σ ) = A • B + i ( A × B ) • σ
Z  = (a+ b i) + (A + B i) • σ ←  paravector
components
    real scalar: a ∈ ℝ
    img scalar: b ∈ ℝ
    real vector: A ∈ ℝ³
    img vector: B ∈ ℝ³
Z = X + i Y
    X = real(Z) = a + A • σ
    Y = img(Z) = b + B • σ
Z = S + V • σ
    S = scalar(Z) = a+ b i 
    V = vector(Z) = A + B i 
( )* ← put a minus sign in front of each σₙ
( )ᴿ ← reverse multiplied order of each σₙ 
S* = Sᴿ = a  – b i
V* = Vᴿ = A – B i
Z*   = S* – V* • σ  = X* – i Y*  
Zᴿ   = S* + V* • σ  = X  – i Y 
Z*ᴿ = S – V • σ     = X* + i Y*
note:  (Zᴿ)ᴿ = Z,  (Z*)* = Z,  (Z*)ᴿ = (Zᴿ)*
Z Z' = ( X X' – Y Y' ) + i ( X Y' + Y X' )
        = ( S S' + V • V' ) + ( S V' + S' V + i ( V × V' ) ) • σ
components of product Z Z'
    real scalar: 
        a a' – b b' + A • A' – B • B'
    img scalar: 
        a b' + b a' + A • B' + B • A'
    real vector: 
        a A' – b B' + a' A – b' B – A × B' – B × A' 
    img vector: 
        a B' + b A' + a' B + b' A + A × A' – B × B'
‖Z‖² = scalar(ZZᴿ)
= |S|² + ‖V‖²
⟪ Z ⟫² = ZZ*ᴿ 
= S² – V • V
______________________________
SUMMARY OF RESULTS BELOW
using natural units, i.e., ℏ = c =1
with metric signature (+, −, −, −)
dX = dt + dr • σ  ← differential space time
ds² = ⟪dX⟫² = dX dX*  ← space-time metric
ɣ = dt/ds = 1 / √( 1 – v² )
dX/ds = ɣ ( 1 + v • σ )  ←proper velocity
m dX/ds = E + P • σ  ← proper momentum
F = ( 𝓔 + i 𝓑 ) • σ  ← electromagnetic field 
∂ = ∂ₜ + σ • ∂ᵣ  ← partial operator
∂ ∂* = ∂ₜ² – ∂ᵣ²  ← d’Alembertian
Φ = V + A • σ  ← potential field 
∂ Φ = S + F*  ← potential gradient
S = ∂ₜ V + ∂ᵣ • A ← Lorenz scalar
∂ F = ρ – J • σ  ← Maxwell's equations
∂*∂ F = ∂* (ρ – J • σ)  ← wave equations
m ddX/ds² = q real( F dX/ds )  ← Lorentz force
(∂∂* + m²) Ψ = 0  ← Klein–Gordon equation
    Ω = ⟪m + σ • ∂ᵣ⟫  ← freq operator
    (∂ₜ ± i Ω) Ψ(±) = 0  ← Ψ = Ψ(+) + Ψ(–)
W(Z) := { 0 , Z ; Z*ᴿ , 0 }  ← 2×2 matrix
(W( i ∂ ) – m 𝐈) Ψ = 0  ← Dirac equation
(–E² + P² + m²) Ψ  = 0  ← KG with field
        E = i ∂ₜ – q V  ← energy operator
        P = –i ∂ᵣ – q A  ← momentum operator
(W( i ∂ – q Φ* ) – m 𝐈) Ψ = 0  ← Dirac with field
2mE Ψ = ( P²  – q 𝓑 • σ) Ψ  ← Pauli equation
    if 𝓑 is a uniform constant
        Ψ = Ψ(+) + Ψ(–)  ← spin up + spin down
        Π(±) = ½ ( 1 ± u • σ )  ← spin projector
        Ψ(±) = Π(±) Ψ
        i ∂ₜ Ψ(±) = ( ½ ( P² ∓ q ‖𝓑‖ ) / m + q V ) Ψ(±)
X' = T⁻¹ X T  ← transforms
    T = exp( ½ θ u • i σ )  ← rotation
    T = exp( ½ θ u • σ )  ← Lorentz boost
______________________________
commutator
[X, Y] := X Y – Y X
[Y, X] = –[X, Y]
[X, Y] = 0 → X Y = Y X → commutative
______________________________
let
A = 3-element vector
B = 3-element vector
dot product
A • B = A₁B₁ + A₂B₂ + A₃B₃
square vector notation
A² = A • A
scalar norm
for  a ∈ ℝ
    |a| = √(a²)
for a ∈ ℝ, b ∈ ℝ, c = a + i b
    |c| = √(cc*)
          = √((a + i b)(a – i b))
          = √(a² + b²)
vector norm
for  A ∈ ℝ³
    ‖A‖ = √(A²)
for A ∈ ℝ³, B ∈ ℝ³, C = A + i B
    ‖C‖ = √(CC*)
          = √((A + i B) • (A – i B))
          = √(A² + B²)
note:  C² ≠ ‖C‖²
cross product
A × B = (A₂B₃ – A₃B₂, 
            –A₁B₃ + A₃B₁, 
              A₁B₂ – A₂B₁)
B × A = –A × B
______________________________
σ algebra
for n and m in {1, 2, 3}
σₙ = algebraic object
σₙ² = +1,   σₙσₘ = –σₘσₙ  if  m ≠ n
all product components can be reduced to
    scalar: 1
    vector: σ₁, σ₂, σ₃
    bi-vector: σ₂σ₃, σ₃σ₁, σ₁σ₂
    pseudo-scalar: σ₁σ₂σ₃
connection to imaginary numbers
i = σ₁σ₂σ₃
i² = –1
(σ₁σ₂σ₃)² = σ₁σ₂(σ₃σ₁)σ₂σ₃ 
= –σ₁(σ₂σ₁)(σ₃σ₂)σ₃
= –σ₁²σ₂²σ₃² 
= –1
bivectors
i σ = vector of bivectors
= {σ₂σ₃, σ₃σ₁, σ₁σ₂}
[i, σₙ] = 0
i σₙ = σₙ i  ←  commutative 
Kronecker delta
δₙₘ = 1 if n = m
       = 0 else
Levi-Civita symbol
εₙₘₖ =  +1 if {1,2,3}, {2,3,1}, {3,1,2} ← even permutation
         = –1 if {3,2,1}, {1,3,2}, {2,1,3}  ← odd permutation 
         = 0  if any two indices are equal
σₙσₘ = δₙₘ + εₙₘₖ i σₖ
= εₙₘₖ i σₖ  if  n ≠ m
σ₁σ₂σ₃ = (δ₁₂ + ε₁₂₃ i σ₃) σ₃ 
= (0 + i σ₃) σ₃ = i
non-commutative algebra if n ≠ m
[σₙ, σₘ] = σₙσₘ – σₘσₙ
= σₙσₘ – (–σₙσₘ)
= 2 σₙσₘ 
= 2 εₙₘₖ i σₖ
[σ₂, σ₃] = 2 σ₁σ₂σ₃σ₁ = 2 i σ₁
[σ₃, σ₁] = 2 σ₁σ₂σ₃σ₂ = 2 i σ₂
[σ₁, σ₂] = 2 σ₁σ₂σ₃σ₃ = 2 i σ₃
______________________________
sigma conjugate operation
put a minus sign in front of each σₙ
σₙ* = –σₙ
i* = –i
(i σₙ )* = i σₙ
(σₙσₘ)* = σₙσₘ
(σ₁σ₂σ₃)* = σ₁* σ₂* σ₃* ← does not flip order
= (–σ₁)(–σ₂)(–σ₃)  
= –σ₁σ₂σ₃
(σₙσₘ)* = σₙ*σₘ*
= (–σₙ)(–σₘ)
= σₙσₘ
( i σₙ )* = i* σₙ* 
= (–i)(–σₙ) 
= i σₙ
______________________________
reverse operation
reverse multiplied order of each σₙ
σₙᴿ = σₙ
iᴿ = –i  ←  (σ₁σ₂σ₃)ᴿ = σ₃σ₂σ₁ = –σ₁σ₂σ₃
(σₙσₘ)ᴿ = σₘσₙ
= –σₙσₘ  if  n ≠ m
(i σₙ)ᴿ = –i σₙ  
    (i σₙ)ᴿ = σₙ ( i )ᴿ 
    = σₙ (–i)
    = –i σₙ
______________________________
sign table
0 = scalar
1 = vector, σₙ
2 = bi-vector, σₙσₘ = εₙₘₖ i σₖ
3 = pseudo-scalar, i = σ₁σ₂σ₃
conjugate operation ( )*
0:  1* = 1
1:  σₙ* = –σₙ
2:  (i σₙ)* = i σₙ
3:  i* = –i
reverse operation ( )ᴿ
0:  1ᴿ = 1
1:  σₙᴿ = σₙ
2:  (i σₙ)ᴿ = –i σₙ
3:  iᴿ = –i
conjugate reverse operation ( )*ᴿ
0:  1*ᴿ = 1
1:  σₙ*ᴿ = –σₙ
2:  (i σₙ)*ᴿ = –i σₙ
3:  i*ᴿ = i
           0  1  2  3
( )*      +  –  +  –
( )ᴿ      +  +  –  –
( )*ᴿ     +  –  –  + 
______________________________
σ vector form
σ = {σ₁, σ₂, σ₃}  ←  symbolic vector 
let
A = 3-element vector
B = 3-element vector
dot product
A • σ = A₁σ₁+ A₂σ₂ + A₃σ₃
addition
A • σ + B • σ = (A + B) • σ
multiplication 
(A • σ) (B • σ) = A • B + i (A × B) • σ    
←    i = σ₁σ₂σ₃
= ΣₙΣₘ AₙBₘ σₙσₘ
= Σₙ AₙBₙσₙ² + Σ_{n<m} (AₙBₘ – AₘBₙ)σₙσₘ
= A • B
    + (A₁B₂ – A₂B₁)σ₁σ₂
    + (A₁B₃ – A₃B₁)σ₁σ₃
    + (A₂B₃ – A₃B₂)σ₂σ₃
= A • B
    + (A₁B₂ – A₂B₁)σ₁σ₂(σ₃²)
    + (A₁B₃ – A₃B₁)σ₁(σ₂²)σ₃
    + (A₂B₃ – A₃B₂)(σ₁²)σ₂σ₃
= A • B
    + (A₁B₂ – A₂B₁)(σ₁σ₂σ₃)σ₃
    – (A₁B₃ – A₃B₁)(σ₁σ₂σ₃)σ₂
    + (A₂B₃ – A₃B₂)(σ₁σ₂σ₃)σ₁
(A • σ)² = A²
(A • σ)²ⁿ = (A²)ⁿ
(A • σ)²ⁿ⁺¹ = (A²)ⁿ (A • σ)
((A • σ)(B • σ))* = (A • σ) (B • σ)
(A • σ) σ = A + i A × σ
σₙ σ = (eₙ • σ) σ 
= eₙ + i eₙ × σ
(A • σ) (B • σ) = A • B + i (A × B) • σ
(B • σ) (A • σ) = A • B – i (A × B) • σ
[A • σ, B • σ] = 2 i (A × B) • σ
(A × B) • σ = [A • σ, B • σ] / (2 i)
A • B = ((A • σ)(B • σ) + (B • σ)(A • σ)) / 2
[A • σ, σ] = –2 i A × σ = 2 i σ × A
[A • (iσ), σ] = 2 A ×σ
[A • σ, σ] • B = 2 i (σ × A) • B
= 2 i (A × B) • σ
= [A • σ, B • σ] 
______________________________
real products
X = a + A • σ ∈ ℝ+ℝ³
Y = b + B • σ ∈ ℝ+ℝ³
XY = (a + A • σ) (b + B • σ)
= a b + a (B • σ) + b (A • σ) + (A • σ) (B • σ)
= ( a b + A • B ) +  (a B + b A + i (A × B)) • σ
= YX  if  A × B = 0
XY ∈ ℝ+ℂ³
XY ∈ ℝ+ℝ³  if  A × B = 0
[ X, Y ] = 2 i (A × B) • σ
= 0 if  A × B = 0  ← conditionally commutative
(XY)ᴿ = ( Y )ᴿ Xᴿ
= ... + iᴿ(A × B) • σ
= ... – i  (A × B) • σ
= ... + i  (B × A) • σ
= Y X
scalar(XY) = a b + A • B
vector(XY) = a B + b A + i (A × B)
scalar(YX) = scalar(XY)
vector(YX) = a B + b A + i (B × A)
= vector(XY)  if  A × B = 0
scalar(X²) = a² + A²
vector(X²) = 2 a A
X² = a² + A² + 2 a A • σ
= b + B • σ
= Y
X = ±√Y
a² = ½ ( b ± √(b² – B²) )
A² = b – a²
minimal polynomial
X² − 2 a X + ( a² − A² ) = 0
⟪XY⟫² = ⟪Y⟫²⟪X⟫²   iff  A × B = 0
X⁻¹ = X*/(XX*)
= (a – A • σ)/(a² – A²)  if  a² ≠ A²
X⁻¹ X = 1
d(X⁻¹) = –X⁻¹ (dX) X⁻¹
X⁻¹ Y = Y X⁻¹  if  A × B = 0
______________________________
real exponent
X = a + A • σ ∈ ℝ+ℝ³
exp(X) = exp(a) exp( A • σ )
exp(X*) = exp(a) exp( –A • σ ) = exp(X)*
exp(–X) = exp(–a) exp( –A • σ ) = exp(X)⁻¹
u = A / ‖A‖  if  A ≠ 0
u = (u₁, u₂, u₃),  u² = 1
exp( A • σ ) = exp( ‖A‖ (u • σ) )
= Σₙ ‖A‖ⁿ (u • σ)ⁿ / n!
= (Σₙ ‖A‖²ⁿ/(2n)!) + (Σₙ ‖A‖²ⁿ⁺¹/(2n+1)!) (u • σ)
= cosh(‖A‖) + sinh(‖A‖) (u • σ)
= ɣ ( 1 + v • σ )
with
ɣ = cosh( ‖A‖ )
v = tanh( ‖A‖ ) u
exp( A • σ )* = ɣ ( 1 – v • σ )
= exp( –A • σ )
= exp( (A • σ)*)
= exp( A • σ )⁻¹
rotation
exp( −i A • σ ) = exp( ‖A‖ (−i u • σ) )
= cos(‖A‖) + sin(‖A‖) (u • (−i σ))
= cos(‖A‖) + sin(‖A‖) (u • q) ← see quaternions
______________________________
real logarithm
X = a + A • σ ∈ ℝ+ℝ³
Y = b + B • σ ∈ ℝ+ℝ³
find (a, A) from (b, B) s.t.
b + B • σ = exp(a + A • σ)
Y = exp(X) 
= exp( a + A • σ )
= exp(a) exp( A • σ ) 
= exp(a) ( cosh(‖A‖) + sinh(‖A‖) (u • σ) )
= b + B • σ
= b + ‖B‖ (u • σ)
let  
‖A‖ = arctanh(‖B‖ / b)  
exp(a) = √( b² − B² )  
u = B / ‖B‖  
log( Y ) = log(exp(X))
= X = a + A • σ
where
a = ½ log( b² − B² )  
A = ‖A‖ u = arctanh( ‖B‖ / b ) ( B / ‖B‖ )
______________________________
paravectors: complex scalar + complex vector
i = σ₁σ₂σ₃
X = a + A • σ ∈ ℝ+ℝ³
Y = b + B • σ ∈ ℝ+ℝ³
Z = c + C • σ ∈ ℂ+ℂ³
  = X + i Y
  = (a+ b i) + (A + i B) • σ  ←  paravector
a = scalar(X) ∈ ℝ
A = vector(X) ∈ ℝ³
b = scalar(Y) ∈ ℝ
B = vector(Y) ∈ ℝ³
c = scalar(Z) ∈ ℂ  ← complex scalar
C = vector(Z) ∈ ℂ³ ← complex vector
Z* = (a – b i) + (–A + i B) • σ
Zᴿ = (a – b i) + (A – i B) • σ
Zᴿ* = (a + b i) – (A + i B) • σ
X = real(Z) = (Z + Zᴿ) / 2
Y = img(Z) = (Z – Zᴿ) / (2i)
c = a + b i
= (Z + Zᴿ*) / 2
C • σ = (A + i B) • σ
= (Z − Zᴿ*) / 2
recovering elements
scalar: 
    a + b i = (Z + Zᴿ*) / 2  
vector: 
    (A + i B) • σ = (Z − Zᴿ*) / 2  
real: 
    a + A • σ = (Z + Zᴿ) / 2
img: 
    b + B • σ = (Z – Zᴿ) / (2i)
real scalar: 
    a = (Z + Z* + Zᴿ + Zᴿ*) / 4
img scalar: 
    b = (Z − Z* − Zᴿ + Zᴿ*) / (4i)
real vector: 
    A • σ = (Z − Z* + Zᴿ − Zᴿ*) / 4
img vector: 
    B • σ = (Z + Z* − Zᴿ − Zᴿ*) / (4i)
______________________________
complex products
X = a + A • σ ∈ ℝ+ℝ³
Y = b + B • σ ∈ ℝ+ℝ³
Z = c + C • σ ∈ ℂ+ℂ³
  = X + i Y
components of Z Z'
    real scalar: a a' – b b' + A • A' – B • B'
    img scalar: a b' + b a' + A • B' + B • A'
    real vector: a A' – b B' + a' A – b' B
       – A × B' – B × A' 
    img vector: a B' + b A' + a' B + b' A
        + A × A' – B × B'
components of Z²
    real scalar: a² – b² + A² – B²
    img scalar: 2 a b + 2 A • B
    real vector: 2 (a A – b B)
    img vector: 2 (a B + b A)
components of Z Z*
    real scalar: a² + b² – A² – B²
    img scalar: 0
    real vector:  –2 A × B
        ←  different than Z*Z
    img vector: 2 ( a B – b A )
components of Z* Z
    real scalar: a² + b² – A² – B²
    img scalar: 0
    real vector: +2 ( A × B )                
        ←  different than ZZ*
    img vector: 2 ( a B – b A )
components of Z Zᴿ
    real scalar: a² + b² + A² + B²
    img scalar: 0
    real vector: 2 (a A + b B + A × B)  
        ← different than ZᴿZ
    img vector: 0
components of Zᴿ Z
    real scalar: a² + b² + A² + B²
    img scalar: 0
    real vector: 2 (a A + b B – A × B)  
        ← different than ZZᴿ
    img vector: 0
components of Z Zᴿ*
    real scalar: a² – b² – A² + B²
    img scalar: 2 a b – 2 A • B
    real vector: 0
    img vector: 0
components of Zᴿ* Z
    real scalar: a² – b² – A² + B²
    img scalar: 2 a b – 2 A • B
    real vector: 0
    img vector: 0
scalar(Z Z') = scalar(Z' Z) ← cyclic property
if Z' = M Z*ᴿ,
    scalar(Z M Z*ᴿ) = scalar(M Z*ᴿ Z)
    = scalar(M) (Z*ᴿ Z)  ←  since Z*ᴿ Z is scalar
scalar(Z Z*)  = c c* – C C* ← real
= (a² + b²) – (A² + B²)
‖Z‖² = |c|² + ‖C‖²
= c c* + C C* 
= (a² + b²) + (A² + B²) ← non-negative real
= scalar(Zᴿ Z) 
= scalar(Z Zᴿ)
= ½ ( Zᴿ Z + (Zᴿ Z)*ᴿ)
= ½ ( Zᴿ Z + Z*ᴿ Z* )
‖ r exp(iθ) Z ‖² = r² ‖Z‖²
Z Z*  = real scalar + complex vector
Z Zᴿ  = real scalar + real vector
Z Zᴿ* = complex scalar
Z* Zᴿ = Zᴿ Z* = (Zᴿ* Z)*
Z* Z = Z Z*  if  A × B = 0 
Zᴿ Z = Z Zᴿ  if  A × B = 0 
( Z Z' )* = Z* ( Z' )*
( Z Z' )ᴿ = ( Z' )ᴿ Zᴿ
( Z Z' )*ᴿ = ( Z' )*ᴿ Z*ᴿ
partial operator on products
∂ = ∂ₜ + σ • ∂ᵣ
∂(ZZ') = (∂Z) Z' + Z (∂Z') + [σ, Z] • ∂ᵣ Z'
______________________________
paravector metric (not a norm)
⟪scalar + vector • σ⟫² := scalar² – vector²
= (scalar + vector • σ)(scalar – vector • σ)
= (scalar – vector • σ)(scalar + vector • σ)
X = a + A • σ ∈ ℝ+ℝ³
Y = b + B • σ ∈ ℝ+ℝ³
Z = c + C • σ ∈ ℂ+ℂ³
= X + i Y 
⟪Z⟫² = c² – C²    ←  complex scalar
= Z Zᴿ* = Zᴿ*Z
= (a + i b)² – (A + i B)²
= (a² – b² – A² + B²) + 2 i (a b – A • B)
⟪ r exp(iθ) Z⟫² = r² exp(2iθ) ⟪Z⟫²
⟪X⟫² = a² – A²  ←  real scalar
= X X* = X* X
______________________________
inverse
Z ∈ ℂ+ℂ³
Z⁻¹ = Zᴿ* / ⟪ Z ⟫² 
= Zᴿ* / (Zᴿ* Z)
note:
if Zᴿ* Z = 0, the inverse does not exist
______________________________
minimal polynomial
Z = c + C • σ ∈ ℂ+ℂ³
Z² − 2 c Z + ( c² − C² ) = 0 
    λ(±) = c ± √( C² )
    ( Z − λ(–) ) ( Z − λ(+) ) = 0
______________________________
complex exponent
X = a + A • σ ∈ ℝ+ℝ³
Y = b + B • σ ∈ ℝ+ℝ³
Z = c + C • σ ∈ ℂ+ℂ³
  = X + i Y 
exp(Z) = exp( c + C • σ )
= exp(c) exp( C • σ )
exp(Z*) = exp(c*) exp( –C* • σ )
exp(Z)⁻¹ = exp(–c) exp( –C • σ )
exp(Z*)⁻¹ = exp(–c*) exp( C* • σ )
(C • σ)² = C • C = C² 
= A² – B² +  2 i A • B
(C • σ)²ⁿ = (C²)ⁿ
z = √C²  ←  complex scalar
exp(C • σ) = Σₙ (C • σ)ⁿ / n!
= Σₙ (C²)ⁿ / (2n)! + Σₙ (C²)ⁿ / (2n+1)! (C • σ)
= cosh(z) + sinh(z) (C • σ) / z
______________________________
complex fields
∂ = ∂ₜ + ∂ᵣ • σ
F = (a + i b) + (A + i B) • σ
a : ℝ×ℝ³ → ℝ
A : ℝ×ℝ³ → ℝ³
b : ℝ×ℝ³ → ℝ
B : ℝ×ℝ³ → ℝ³
F : ℝ×ℝ³ → ℂ+ℂ³
components of ∂ F
  real scalar: 
    ∂ₜ a + ∂ᵣ • A
  img scalar: 
    ∂ₜ b + ∂ᵣ • B
  real vector: 
    ∂ᵣ a + ∂ₜ A − ∂ᵣ × B
  img vector: 
    ∂ᵣ b + ∂ᵣ × A + ∂ₜ B
components of i ∂ F
  real scalar: 
    −∂ₜ b − ∂ᵣ • B
  img scalar: 
    ∂ₜ a  + ∂ᵣ • A
  real vector: 
    −∂ᵣ b – ∂ᵣ × A − ∂ₜ B
  img vector: 
    ∂ᵣ a + ∂ₜ A − ∂ᵣ × B
______________________________
space–time numbers
dX = dt + dr₁σ₁+ dr₂σ₂ + dr₃σ₃
with conjugate
dX* = dt – dr₁σ₁ – dr₂σ₂ – dr₃σ₃
the space-time metric is
ds² = dt² – dr₁² – dr₂² – dr₃²
scalar: dt in ℝ, vector: dr in ℝ³
dX  = dt + dr • σ
dX* = dt – dr • σ
components
(dX + dX*) / 2 = dt
(dX – dX*) / 2 = dr • σ
scalar(dX) = dt
vector(dX) = dr
paravector metric
⟪dX⟫² = dX dX* = dt² – dr²
1.  ⟪dX⟫² > 0 for all dX ≠ 0 if dt² > dr²
2.  ⟪dX + dX'⟫² – ⟪dX⟫² – ⟪dX'⟫²
        = (dt + dt')² – (dr + dr')² 
           – (dt² – dr²) – (dt'² – dr'²)
        = 2 ( dt dt'– dr • dr' )
3.  ⟪a dX⟫² = a² (dt² – dr²) = a² ⟪dX⟫²  
[dX, dX'] = [ dr • σ, dr' • σ ]
= 2i (dr × dr') • σ
ds² = dt² – dr²
= ⟪dX⟫² = dX dX* = dX*dX
ds = √( dt² – dr² ) 
= √( 1 – dr²/dt² ) dt
= √(1 – v²) dt
define
v = dr/dt = (v₁, v₂, v₃)  in ℝ³
ɣ = dt/ds = 1 / √( 1 – v² )
note:
dX⁻¹ = dXᴿ* / ⟪ dX ⟫² 
= dX* / (dX dX*)
= dX* / ds²
______________________________
if ⟪dX⟫² > 0, the interval is "time-like"
    the arc length s equals the proper time τ 
    (measured along a path)
    τ = s = ∫ ds = ∫ dτ
    = ∫ √(1 – v²) dt  ←  time dilation
    =  ∫ 1/ɣ dt
    note: in SI units, ds = c dτ
if ⟪dX⟫² = 0, the interval is "light-like"
    dt² = dr²  →  dt = ± ‖dr‖
    dX⁻¹ does not exist
if ⟪dX⟫² < 0, the interval is "space-like"
______________________________
4 velocity  ← a.k.a  proper velocity
dX/ds = (dt/ds) + (dr/ds) • σ
= (dt + σ • dr) / √( dt² – dr² )
= (dt/dt + σ • dr/dt) / √( (dt/dt)² – (dr/dt)² )
= (1 + v • σ) / √( 1 – v² )
= ɣ ( 1 + v • σ )
components of dX/ds
  scalar: ɣ
  vector: ɣ v
⟪dX/ds⟫² = ɣ² ( 1 + v • σ ) ( 1 – v • σ )
= ɣ² ( 1 – v² ) 
= 1
dX/dt = (dX/ds)(ds/dt) 
= (dX/ds)( 1/ɣ )
define
U := dX/ds
= ɣ + ɣ v • σ
U = Uₛ + Uᵥ • σ
    Uₛ = ɣ
    Uᵥ = ɣ v
ds² = dX dX*,  ds > 0
dX*/ds = dX⁻¹ ds
= (dX/ds)⁻¹    →    U⁻¹ = U*
⟪U⟫² = UU* = 1
d⟪U⟫²/ds
= dU/ds U* + U (dU/ds)*
= ( dU/ds U* ) + ( dU/ds U* )*ᴿ
→ scalar( U* dU/ds ) = 0
⟪dX⟫² = dX dX* 
= ((dX/ds) ds) ((dX/ds)* ds)
= U U* ds²
= ⟪U⟫² ds²
= ds²
___________________________
4 acceleration ← a.k.a  proper acceleration
a = dv/dt = {a₁, a₂, a₃}  in ℝ³
dɣ/dt = ɣ³ (v • a)
dɣ/ds = ɣ⁴ (v • a)
∂ᵥ(–1/ɣ) = ɣv
∂ᵥ( ɣ v ) = ɣ 𝐈 + ɣ³ v vᵀ
d(ɣv)/dt = ∂ᵥ( ɣ v ) dv/dt
= ɣ a + ɣ³ (v • a) v
d(ɣ v)/ds = (d(ɣ v)/dt) (dt/ds)
= ɣ² a + ɣ⁴ (v • a) v
recall
U = Uₛ + Uᵥ • σ
    Uₛ = ɣ
    Uᵥ = ɣ v
dU/ds = (dɣ/ds) + (d(ɣv)/ds) • σ  
= ɣ⁴ (v • a) + (ɣ⁴ (v • a) v + ɣ² a) • σ
u = v / ‖v‖,  u² = 1
par(a) = (a • u) u
prp(a) = a − par(a)
a² = par(a)² + prp(a)²
(v • a)² = v² par(a)²
⟪dU/ds⟫² = –ɣ⁴ ( a² + ɣ² (v • a)² )
= –ɣ⁴ ( prp(a)² + ɣ² par(a)² )
let
A = dU/ds 
= ddX/ds²
A = Aₛ + Aᵥ • σ
    Aₛ = ɣ⁴ (v • a)
    Aᵥ = ɣ⁴ (v • a) v + ɣ² a
U A = ( Uₛ Aₛ + Uᵥ • Aᵥ ) 
         +  (Uₛ Aᵥ + Aₛ Uᵥ + i (Uᵥ × Aᵥ)) • σ
components of U A
  real scalar:  2 ɣ⁵ (v • a)
  real vector: ɣ³ a + 2 ɣ⁵ (v • a) v
  img vector: ɣ³ ( v × a )
U* A = ( Uₛ Aₛ – Uᵥ • Aᵥ ) 
         +  (Uₛ Aᵥ – Aₛ Uᵥ – i (Uᵥ × Aᵥ)) • σ
components of U* A
  real scalar:  0
  img scalar: 0
  real vector: ɣ³ a
  img vector: –ɣ³ ( v × a )
note:  scalar( U* A ) = 0
______________________________
4 momentum ← a.k.a  proper momentum 
m U = E + P • σ
relativistic energy
E := scalar( m U )
= m ɣ
if v = 0 
    E = m  ←  rest energy
        = mc²  ←  in SI units
E – m = kinetic energy
relativistic momentum
P := vector(m U) 
= m ɣ v
⟪m U⟫² = m²
= (E + P • σ)(E – P • σ)
= E² – P²
E² = P² + m²
note: if m = 0,  E² = P²  ←  light
______________________________
Lagrangian
action =  ∫ ℒ(t, r, v) dt
Euler-Lagrange equation
( ∂ᵣ – d/dt ∂ᵥ ) ℒ = 0  
→  d/dt ∂ᵥ ℒ  = ∂ᵣ ℒ
Noether's conserved quantities
    canonical momentum
    ∂ᵥ ℒ = constant if ∂ᵣ ℒ = 0
    Hamiltonian
    ℋ = v • (∂ᵥ ℒ) – ℒ
    = (v • ∂ᵥ – 1) ℒ
    = constant if ∂ₜ ℒ = 0
    d/dt ℋ = d/dt (∂ᵥ ℒ) • v + (∂ᵥ ℒ) • dv/dt – dℒ/dt
    = (∂ᵣ ℒ) • v + (∂ᵥ ℒ) • dv/dt – dℒ/dt  ← EL
    = (∂ᵣ ℒ) • v + (∂ᵥ ℒ) • dv/dt 
        – (∂ₜ ℒ + (∂ᵣ ℒ) • v + ∂ᵥ ℒ • dv/dt)
    = – ∂ₜ ℒ
    d/dt ℋ = 0   if   ∂ₜ ℒ = 0
______________________________
action of a free particle = 
mass-weighted proper time along its world line
∫ ℒ dt = –m τ
= –m ∫ ⟪dX⟫           = –m ∫ ds 
= –m ∫ ⟪dX/dt⟫ dt   = –m ∫ (ds/dt) dt
= –m ∫ ⟪U / ɣ⟫ dt    = –m ∫ (1/ɣ) dt
= –m ∫ √(1 – v²) dt
ℒ = –m / ɣ    ←  Lagrangian
∂ᵥ ℒ = m ɣ v = P  →  ( ∂ᵥ • σ ) ℒ = P •  σ
ℋ = m ɣ v • v + m / ɣ
= m ɣ ( v² + 1 / ɣ² )
= m ɣ = E
∂ₜ ℒ = 0  ← m ɣ is constant
∂ᵣ ℒ = 0  ← m ɣ v is constant
______________________________
non-relativistic approximation
ℒ = –m / ɣ
= –m + ½ m v² + O(v⁴)
momentum
∂ᵥ ℒ = ∂ᵥ (½ m v² + O(v⁴)) 
= m v + O(v³) 
≈ m v
total energy = kinetic + rest mass energy
ℋ = (m v + O(v³)) • v + m – ½ m v² – O(v⁴)
= ½ m v² + m + O(v⁴) 
≈ ½ m v² + m
_____________________________
electricity & magnetism
vector
∂ᵣ = {∂/∂r₁, ∂/∂r₂, ∂/∂r₃}
∂ᵢ = ∂ / ∂rᵢ
∂ᵣ² = ∂₁² + ∂₂² + ∂₃²
∂ = ∂ₜ + σ • ∂ᵣ
= ∂ₜ + σ₁∂₁+ σ₂∂₂ + σ₃∂₃
vectors
𝓔 = {𝓔₁, 𝓔₂, 𝓔₃}    ← electric field 
𝓑 = {𝓑₁, 𝓑₂, 𝓑₃}  ← magnetic field
𝓕 = 𝓔 + i 𝓑         ← electro-magnetic field
F = (𝓔 + i 𝓑) • σ,    i = σ₁σ₂σ₃
= 𝓕 • σ
F* = (−𝓔 + i 𝓑) • σ
= −𝓕* • σ
components of  ½ F²
    real scalar:  ½ (𝓔² – 𝓑²)
    img scalar:  i 𝓔 • 𝓑
    vector: 0
components of  ½ F Fᴿ
    scalar:  ½ (𝓔² + 𝓑²)    ← energy density
    vector:  𝓔 × 𝓑             ← Poynting vector
density paravector
D = ρ + J • σ
D* = ρ − J • σ
ρ = scalar charge density
J = vector current density
∂ F = D*       ← Maxwell's equations
= ρ − J • σ  
∂* F* = D
this expands to
∂F = (∂ₜ + σ • ∂ᵣ) ( (𝓔 + i 𝓑) • σ) 
= ( ∂ᵣ • 𝓔 ) 
+ ( ∂ᵣ • 𝓑 ) i
+ (∂ₜ 𝓔 – ∂ᵣ × 𝓑) • σ 
+ (∂ₜ 𝓑 + ∂ᵣ × 𝓔) • (i σ) 
components of ∂F = ρ − J • σ
    real scalar: ∂ᵣ • 𝓔 = ρ
    img scalar: ∂ᵣ • 𝓑 = 0
    real vector: ∂ₜ 𝓔 – ∂ᵣ × 𝓑 = –J
    img vector: ∂ₜ 𝓑 + ∂ᵣ × 𝓔 = 0
re–arranging the last two
∂ᵣ × 𝓑 = J + ∂ₜ 𝓔
∂ᵣ × 𝓔 = –∂ₜ 𝓑
______________________________
wave operator
∂ = ∂ₜ + σ • ∂ᵣ
∂* = ∂ₜ – σ • ∂ᵣ
d’Alembertian
□:= ∂ₜ² – ∂ᵣ² 
= ∂ ∂* = ∂* ∂ 
= ⟪∂⟫²
note: with σ-algebra
∂² = ∂ₜ² + ∂ᵣ² + 2 ∂ₜ ∂ᵣ • σ
≠ □
note: definitions differ across sources
some use the  (–,+,+,+) signature
some use □²
wave equation
□ F = ∂* ∂ F
= ∂* (ρ − J • σ)
□ F* = ∂ (ρ + J • σ)
( □ 𝓔 + i □ 𝓑 ) • σ
= (∂ₜ ρ + ∂ᵣ • J) + (–∂ₜ J – ∂ᵣ ρ + i (∂ᵣ × J)) • σ
components of □ F = ∂* (ρ − J • σ)
    real scalar: 0 = ∂ₜ ρ + ∂ᵣ • J
    img scalar: 0 = 0
    real vector: □ 𝓔 = –∂ₜ J – ∂ᵣ ρ
    img vector: □ 𝓑 = ∂ᵣ × J
∂ₜ ρ = – ∂ᵣ • J  ←  conservation of charge
if ρ = 0 and J = 0
□ 𝓔 = 0  →  ∂ₜ² 𝓔 = ∂ᵣ² 𝓔
□ 𝓑 = 0  →  ∂ₜ² 𝓑 = ∂ᵣ² 𝓑
______________________________
potential fields
F* = (–𝓔 + i 𝓑) • σ
Φ = V + A • σ
∂Φ = ∂ₜ V + ∂ᵣ • A 
    + ( ∂ₜ A + ∂ᵣ V ) • σ 
    + i ( ∂ᵣ × A ) • σ
S := scalar( ∂Φ ) 
= ∂ₜ V + ∂ᵣ • A 
target relationships
𝓔  = – ∂ₜ A – ∂ᵣ V
𝓑 = ∂ᵣ × A
consider
∂Φ = S + F*
F* = ∂Φ – S
scalar(F*) = scalar(∂Φ – S) = 0
vector(F*) = vector( ∂Φ )
note: always differentiate first
(∂Φ)*ᴿ ≠ Φ*ᴿ ∂*ᴿ
without S
F* = (∂Φ – (∂Φ)*ᴿ) / 2
components of  ∂Φ = S + F*
    real scalar: ∂ₜ V + ∂ᵣ • A  = S
    img scalar: 0 = 0
    real vector:  ∂ₜ A + ∂ᵣ V = –𝓔
    img vector:  ∂ᵣ × A = 𝓑
______________________________
density from potential
recall
F = (𝓔 + i 𝓑) • σ
Φ = V + A • σ
S = ∂ₜ V + ∂ᵣ • A
F* = ∂Φ – S  →  F = ∂* Φ* – S
∂F = ∂∂* Φ* – ∂S
= □ Φ* – ∂S
= ( □ V – ∂ₜ S) + (– □ A – ∂ᵣ S) • σ
= ρ − J • σ
gives
ρ = □ V – ∂ₜ S
J = □ A + ∂ᵣ S
if S = 0,  
    Lorenz gauge condition
     ∂ₜ V = –∂ᵣ • A
    ρ = □ V
    J = □ A
______________________________
gauge transform
scalar( ∂ Φ ) = S
S = ∂ₜ V + ∂ᵣ • A
λ = scalar field
Φ' = Φ + ∂* λ
∂Φ' = ∂(Φ + ∂*λ)
= ∂Φ + ∂(∂*λ)
= ∂Φ + □ λ
note:  □ λ is a scalar field
scalar( ∂ Φ' ) = S + □ λ
vector( ∂ Φ' ) = vector( ∂Φ )
F' = F = (𝓔 + i 𝓑) • σ
______________________________
4 force ← a.k.a  proper force
m U = m dX/ds
= m ɣ + m ɣ v • σ
d/ds = ɣ d/dt
relativistic work 
d/dt (m ɣ) = q v • 𝓔
note: 𝓔 = electric field (not energy)
note: 𝓑 does no work
relativistic Lorentz force law
P := m ɣ v
dP/dt = q 𝓔 + q v × 𝓑
in one formula
d/dt (m U) = q v • 𝓔 + (q 𝓔 + q v × 𝓑) • σ
consider
F U / ɣ = ((𝓔 + i 𝓑) • σ)  ( 1 + v • σ ) 
= 𝓔 • σ  + i 𝓑 • σ + (𝓔 • σ) (v • σ)
    + i (𝓑 • σ)(v • σ)
= 𝓔 • σ  + i 𝓑 • σ + 𝓔 • v 
    + i (𝓔 × v) • σ  + i 𝓑 • v – (𝓑 × v) • σ 
components of F U / ɣ
  real scalar:  𝓔 • v
  img scalar:  𝓑 • v
  real vector:  𝓔 + v × 𝓑
  img vector: 𝓑 + 𝓔 × v
real( F U ) = ½ ( F U + U Fᴿ )
= ɣ (v • 𝓔 + ( 𝓔 + v × 𝓑 ) • σ) 
note: img( F U ) swaps 𝓔 and 𝓑
m dU/ds = q real( F U ) 
=  ½ q ( F U + U Fᴿ )
m ddX/ds² = q real( F dX/ds )  ←  4 force
______________________________
Lagrangian derivation
recall
Φ = V + A • σ
U = dX/ds
ɣ = dt/ds
action = –m ∫ ds – q ∫ scalar( Φ* U ) ds
= ∫ ( –m / ɣ – q scalar( Φ* U / ɣ ))  dt
= ∫ ( –m / ɣ 
    – q scalar( (V – A • σ ) (1 + v • σ )) ) dt
= ∫ ( –m / ɣ – q V + q v • A ) dt
ℒ = –m / ɣ – q ( V – v • A )  ←  Lagrangian
verification
∂ᵥ ℒ = m ɣ v + q A 
= P + q A
d/dt ∂ᵥ ℒ = dP/dt + q ∂ₜ A + q (v • ∂ᵣ) A
∂ᵣ ℒ = − q ∂ᵣ V + q ∂ᵣ (v • A)
vector calc Identity
(v(t) • ∂ᵣ) A(t, r) = ∂ᵣ ( v • A ) – v × ( ∂ᵣ × A )
solving d/dt ∂ᵥ ℒ  = ∂ᵣ ℒ for dP/dt
dP/dt = − q ∂ᵣ V − q ∂ₜ A 
    + q ( ∂ᵣ (v • A) − (v • ∂ᵣ) A )
= q ( –∂ᵣV – ∂ₜA + v × ( ∂ᵣ × A ) )
= q ( 𝓔 + v × 𝓑 )
= m vector( dU/ds )  / ɣ
ℋ = (∂ᵥ ℒ) • v – ℒ
= m ɣ v • v + q v • A + m / ɣ 
    + q V – q v • A
= m ɣ ( v² + 1 / ɣ² ) + q V
= m ɣ ( v² + 1 – v² ) + q V
= m ɣ + q V
recall the Hamiltonian property
dH/dt = – ∂ₜ ℒ
dℋ/dt = d(m ɣ)/dt + q dV/dt
= d(m ɣ)/dt + q ∂ₜV + q v • ∂ᵣ V
∂ₜ H = q ∂ₜ V – q v • ∂ₜ A
d(m ɣ)/dt = q v • ( –∂ᵣ V – ∂ₜ A )
= q v • 𝓔
______________________________
stress–energy of EM?
F = (𝓔 + i 𝓑) • σ
Fᴿ = (𝓔 – i 𝓑) • σ
T = ½ F Fᴿ 
= ½ (𝓔² + 𝓑²) + (𝓔 × 𝓑) • σ
= Tₛ + Tᵥ • σ
∂ = ∂ₜ + σ • ∂ᵣ
components of ∂ T
    real scalar: ∂ₜ Tₛ + ∂ᵣ • Tᵥ
    img scalar: 0
    real vector: ∂ᵣ Tₛ + ∂ₜ Tᵥ
    img vector: ∂ᵣ × Tᵥ
recall for any paravectors Z and Z'
∂(ZZ') = (∂Z) Z' + Z (∂Z') + [σ, Z] • ∂ᵣ Z'
∂ T = ½ ∂(FFᴿ) 
= ½(∂F) Fᴿ + ½F (∂Fᴿ) + ½[σ, F] • ∂ᵣ Fᴿ
D = ρ + J • σ  ←  density paravector
D* = ρ − J • σ
∂ = ∂ₜ + σ • ∂ᵣ
∂ F = D*
F² = (𝓔² − 𝓑²) + 2 i (𝓔 • 𝓑)  ← scalar
components of ∂ F²
    real scalar: ∂ₜ (𝓔² − 𝓑²) 
    img scalar: ∂ₜ 2 (𝓔 • 𝓑)
    real vector: ∂ᵣ (𝓔² − 𝓑²)
    img vector: ∂ᵣ 2 (𝓔 • 𝓑)
components of F D*
    real scalar: −𝓔 • J
    img scalar: −𝓑 • J
    real vector: ρ 𝓔 + 𝓑 × J
    img vector: ρ 𝓑 − 𝓔 × J
components of D* F
    real scalar: −𝓔 • J
    img scalar: −𝓑 • J
    real vector: ρ 𝓔 − 𝓑 × J
    img vector: ρ 𝓑 + 𝓔 × J
components of D* F + F D*
    real scalar: −2 𝓔 • J
    img scalar: −2 𝓑 • J
    real vector: 2 ρ 𝓔
    img vector: 2 ρ 𝓑
recall for any paravectors Z and Z'
∂(ZZ') = (∂Z) Z' + Z (∂Z') + [σ, Z] • ∂ᵣ Z'
∂(F²) = (∂F) F + F (∂F) + [σ, F] • ∂ᵣ F
= D* F + F D* + [σ, F] • ∂ᵣ F
F = Fᵥ • σ
[σ, F] = 2 i Fᵥ × σ
[σ, F] • ∂ᵣ F = −2 i Fᵥ • (∂ᵣ × Fᵥ) 
    −2 ( (∂ᵣ • Fᵥ) Fᵥ − ½ ∂ᵣ Fᵥ² ) • σ
□ F² = ∂* ∂(F²)  ← scalar
______________________________
relativistic quantum mechanics
∂ = ∂ₜ + σ • ∂ᵣ
∂* = ∂ₜ – σ • ∂ᵣ
(i ∂)* = –i ∂*
recall
□:= ∂ₜ² – ∂ᵣ² 
= ∂ ∂* = ∂* ∂ 
= ⟪∂⟫²
note: all differential operators act left to right
bra-ket notation
    normalized
    ⟨ψ | ψ⟩ = expectation of ψᴿ ψ = 1
    expectation
    ⟨ψ | M | ψ⟩ = expectation of ψᴿ M ψ
note: ( )ᴿ is equivalent to the Hermetian operator
______________________________
spin 0 wave mechanics
⟪m U⟫² = m²    →    –E² + P² + m² = 0
for a free particle
E = i ∂ₜ
P = –i ∂ᵣ
note: E + P • σ = i ∂*
(–E² + P² + m²) Ψ
= (–(i ∂ₜ)² – ∂ᵣ² + m²) Ψ
= ( ∂ₜ² – ∂ᵣ² + m² ) Ψ
= (∂ ∂* + m²) Ψ
= (□ + m²) Ψ 
(□ + m²) Ψ = 0  ←  Klein–Gordon equation
note: Ψ is a complex scalar
______________________________
positive and negative frequency
Ω = √(m² – ∂ᵣ²) 
= ⟪m + σ • ∂ᵣ⟫
□ + m² = ∂ₜ² + Ω²
= ∂ₜ² + (m² – ∂ᵣ²)
= (∂ₜ – i Ω) (∂ₜ + i Ω)
= (∂ₜ + i Ω) (∂ₜ – i Ω)
(∂ₜ ± i Ω) Ψ(±) = 0 
Ψ(±) = c(±) exp(∓i Ω t )
Ψ = Ψ(+) + Ψ(–)
______________________________
Interpreting root operators with Fourier
Ψ = Σₙ cₙ fₙ
fₙ = exp(i (ωₙ t + kₙ • r))
□ Ψ = Σₙ ( kₙ² – ωₙ²) cₙ fₙ
±√□ Ψ = Σₙ ±√( kₙ² – ωₙ²) cₙ fₙ
(□ + m²) Ψ = Σₙ (kₙ² + m² – ωₙ² ) cₙ fₙ
kₙ² + m² = ωₙ²
□ + m² = (m ∓ i √□)(m ± i √□)
ωₙ = ±√(kₙ² + m²)
Ω² Ψ =  Σₙ ( m² + kₙ²) cₙ fₙ
Ω Ψ =  Σₙ √( m² + kₙ²) cₙ fₙ
______________________________
spin ½ wave mechanics
consider the 2-element vector
Ψ = { Ψ₁ ; Ψ₂ }
define the 2×2 operator matrix
W(Z) := { 0 , Z ; Z*ᴿ , 0 }  ←  2×2 matrix
W(Z)² = diag( Z Z*ᴿ , Z*ᴿ Z )
if a and b are real scalars,
W((a + b i) Z) = (a + b i) W(Z)
W( i ∂ ) = i { 0 , ∂ ; ∂* , 0 }
= i W( ∂ )
W( i ∂ )² = –diag( ∂∂*, ∂*∂ ) 
= – □ 𝐈
define the Dirac operator matrix
D = W( i ∂ ) – m 𝐈
and its mass conjugate
D' = W( i ∂ ) + m 𝐈
with
D' D  = (W + m 𝐈)(W – m 𝐈)
= (W² + m W – m W – m² 𝐈)
= –( □ + m² ) 𝐈
D Ψ  = 0   ←  Dirac equation
0 = det( D – s 𝐈 ) 
= det( { –(s+m) , i ∂ ; i ∂* , –(s+m) })
= (s + m)² + □  →  s = –m ± i √□
2 equations
i ∂ Ψ₂ – m Ψ₁ = 0
i ∂* Ψ₁ – m Ψ₂ = 0
solve Ψ₁ in first equation
Ψ₁ = i ∂ Ψ₂ / m
plug into second equation
–∂* ∂ Ψ₂ / m – m Ψ₂ = 0    
→         –(□ + m²) Ψ₂ = 0
using the scalar solution Ψ' from KG
with constant complex scalars a and b
Ψ₁ = ( m a + i b ∂ ) Ψ' / m
Ψ₂ = ( i a ∂* + m b ) Ψ' / m
gives
    i ∂ Ψ₂ – m Ψ₁
    = ( –∂ a ∂* + i m ∂ b – m² a – i m b ∂ ) Ψ' / m
    = ( –a / m ) ( □ + m² ) Ψ' = 0
    i ∂* Ψ₁ – m Ψ₂ 
    = ( m i ∂* a – ∂* b ∂ – m i a ∂* – m² b ) Ψ' / m
    = ( –b / m ) ( □ + m² ) Ψ' = 0
with a = 1 and b = 0
    Ψ₁ = Ψ' ← scalar
    Ψ₂ = ( i / m ) ∂* Ψ'
______________________________
Klein–Gordon with potential fields
ℒ = –m / ɣ – q ( V – v • A )
note: E = energy (not electric field)
Hamiltonian
i ∂ₜ = ℋ = E + q V
canonical momentum
–i ∂ᵣ = ∂ᵥ ℒ = P + q A
E = i ∂ₜ – q V ← scalar operator
P = –i ∂ᵣ – q A ←vector operator
note:
∂ₜ V Ψ = (∂ₜ V) Ψ + V ∂ₜ Ψ
∂ₘ Aₙ Ψ = (∂ₘ Aₙ) Ψ + Aₙ ∂ₘ Ψ
E² Ψ = E ( i ∂ₜ Ψ – q V Ψ)
=  –∂ₜ² Ψ – i q ∂ₜ (V Ψ) – i q V ∂ₜ Ψ + q² V² Ψ
= (–∂ₜ² – i q (∂ₜV) – 2 i q V ∂ₜ + q² V² ) Ψ
P² Ψ = P • ( –i ∂ᵣ Ψ – q A Ψ)
=  –∂ᵣ² Ψ + i q ∂ᵣ • (A Ψ) + i q A • ∂ᵣ Ψ + q² A² Ψ
= (–∂ᵣ² + i q (∂ᵣ • A) + 2 i q A • ∂ᵣ + q² A² ) Ψ
recall 
□ = ∂ₜ² – ∂ᵣ² 
Φ = V + A • σ
S = ∂ₜ V + ∂ᵣ • A
⟪∂⟫² = □ ← scalar operator
⟪Φ⟫² = V² − A²  ← scalar
Klein–Gordon equation
0 = (–E² + P² + m²) Ψ 
= ( ∂ₜ² –∂ᵣ² 
    + i q (∂ₜV + ∂ᵣ • A ) 
    + 2 i q (V ∂ₜ + A • ∂ᵣ) 
    – q² ( V² – A² ) + m² ) Ψ
= ( □ + i q S – q² ⟪Φ⟫² + m²
    + 2 i q (V ∂ₜ + A • ∂ᵣ) ) Ψ
recall
𝓔 = –∂ₜ A – ∂ᵣ V
𝓑 = ∂ᵣ × A
F = ( 𝓔 + i 𝓑 ) • σ
commutators
[E, Pₘ] Ψ
= (i ∂ₜ – q V)(–i ∂ₘ Ψ – q Aₘ Ψ) 
    – (–i ∂ₘ – q Aₘ)(i ∂ₜ Ψ – q V Ψ)
= – i q ( (∂ₜ Aₘ) Ψ + Aₘ ∂ₜ Ψ ) 
    – i q ( (∂ₘ V) Ψ + V ∂ₘ Ψ )
    + i q V ∂ₘ Ψ + i q Aₘ ∂ₜ Ψ
= – i q ( ∂ₜ Aₘ + ∂ₘ V ) Ψ
=  i q 𝓔ₘ  →  [E, P] = i q 𝓔
[Pₘ, Pₙ] Ψ
= (–i ∂ₘ – q Aₘ)(–i ∂ₙ Ψ – q Aₙ Ψ) 
  –(–i ∂ₙ – q Aₙ)(–i ∂ₘ Ψ – q Aₘ Ψ) 
= i q ( ∂ₘ Aₙ – ∂ₙ Aₘ ) Ψ
P ⨯ P = ( [P₂, P₃], –[P₁, P₃], [P₁, P₂] )
= i q ( ∂ᵣ × A ) 
= i q 𝓑
(E + σ • P) (E – σ • P) 
    = E² – P² + σ • ( –[E, P] – i P ⨯ P )
    = E² – P² – i q F
(E – σ • P) (E + σ • P) 
    = E² – P² + σ • ( +[E, P] – i P ⨯ P )
    = E² – P² – i q F*
E + σ • P = i ∂* – q Φ
E – σ • P = i ∂ – q Φ*  ←   (E + σ • P)*ᴿ
⟪E + σ • P⟫² 
= E² – P² + σ • ( ±[E, P] – i P ⨯ P )
note: operators are not unique scalars
–(E + σ • P) (E + σ • P)*ᴿ
= (∂* + i q Φ) (∂ + i q Φ*)
–(E + σ • P)*ᴿ (E + σ • P)
= (∂ + i q Φ*) (∂* + i q Φ)
(E + σ • P) (E – σ • P) Ψ
= –(∂* + i q Φ) (∂ + i q Φ*) Ψ
= –□ Ψ – i q ∂* Φ* Ψ – i q Φ ∂ Ψ + q² ⟪Φ⟫²  Ψ
= (E² – P² – i q F) Ψ
= –□ Ψ –i q S Ψ + q² ⟪Φ⟫² Ψ
    – 2 i q (V ∂ₜ + A • ∂ᵣ) Ψ – i q F Ψ
∂* Φ* Ψ + Φ ∂ Ψ 
= S Ψ + 2 (V ∂ₜ + A • ∂ᵣ) Ψ + F Ψ
(E – σ • P) (E + σ • P) Ψ
= – (∂ + i q Φ*) (∂* + i q Φ) Ψ
= –□ Ψ – i q ∂ Φ Ψ – i q Φ* ∂* Ψ + q² ⟪Φ⟫²  Ψ
= (E² – P² – i q F*) Ψ
= –□ Ψ –i q S Ψ + q² ⟪Φ⟫² Ψ
    – 2 i q (V ∂ₜ + A • ∂ᵣ) Ψ – i q F* Ψ
∂ Φ Ψ + Φ* ∂* Ψ 
= S Ψ + 2 (V ∂ₜ + A • ∂ᵣ) Ψ + F* Ψ
(∂* Φ* Ψ + Φ ∂ Ψ) – (∂ Φ Ψ + Φ* ∂* Ψ)
= [Φ, ∂] Ψ + [∂, Φ]* Ψ
= (F – F*) Ψ
 = 2 (𝓔 • σ) Ψ
summary
  E := i ∂ₜ – q V   ← scalar operator
  P := –i ∂ᵣ – q A ← vector operator
  ∂ := ∂ₜ + σ • ∂ᵣ
  ∂* := ∂ₜ – σ • ∂ᵣ
  Φ := V + A • σ
  F := ( 𝓔 + i 𝓑 ) • σ
  S := ∂ₜ V + ∂ᵣ • A
  𝓔 := –∂ₜ A – ∂ᵣ V
  𝓑 := ∂ᵣ × A
( E ± σ • P )( E∓σ • P) 
= E² − (σ • P)² ∓ [E, σ • P]
= E² − P² − i q ( ±𝓔 + i 𝓑) • σ
Klein–Gordon equation
0 = (–E² + P² + m²) Ψ 
= ( □ + i q S – q² ⟪Φ⟫² + m²
    + 2 i q (V ∂ₜ + A • ∂ᵣ) ) Ψ
= (∂* + i q Φ) (∂ + i q Φ*) Ψ – i q F  Ψ + m² Ψ
= (∂ + i q Φ*) (∂* + i q Φ) Ψ – i q F* Ψ + m² Ψ 
= scalar( (∂ + i q Φ*) (∂* + i q Φ) ) Ψ + m² Ψ 
______________________________
Dirac with potential field
Ψ = { Ψ₁ ; Ψ₂ }
define
W(Z) := { 0 , Z ;  Z*ᴿ , 0 }
W(Z)² = { Z Z*ᴿ , 0 ; 0 , Z*ᴿ Z }
W(Z)* = W(Z*)
W(Z)ᴿ = W(Zᴿ)
Dirac operator with potential
D = W( i ∂ – q Φ* ) – m 𝐈
= i W( ∂ ) – q W( Φ* ) – m 𝐈
W( ∂ )² = □ 𝐈 
W( ∂ ) W( Φ* ) 
= { 0 , ∂ ; ∂* , 0 } { 0 , Φ* ; Φ , 0 } 
= diag( ∂ Φ , ∂* Φ* ) 
W( Φ* ) W( ∂ )
= { 0 , Φ* ; Φ , 0 } { 0 , ∂ ; ∂* , 0 } 
= diag( Φ* ∂* , Φ ∂ ) 
W( ∂ ) W( Φ* ) + W( Φ* ) W( ∂ )
= diag( ∂ Φ + Φ* ∂*,  ∂* Φ* + Φ ∂ )
W( Φ* )² 
= { 0 , Φ* ; Φ , 0 } { 0 , Φ* ; Φ , 0 } 
= diag( Φ* Φ , Φ Φ* ) 
= ⟪Φ⟫² 𝐈 
W( i ∂ – q Φ* )² 
= –W( ∂ )² 
    – i q ( W( ∂ ) W( Φ* ) + W( Φ* ) W( ∂ ) ) 
    + q² W( Φ* )² 
= –□ 𝐈 
    – i q diag( ∂ Φ + Φ* ∂*, ∂* Φ* + Φ ∂ ) 
     + q² ⟪Φ⟫² 𝐈 
= –diag(
     (∂ + i q Φ*)(∂* + i q Φ),
     (∂* + i q Φ)(∂ + i q Φ*))
D Ψ = 0  ←  Dirac equation
2 equations 
(–m)Ψ₁ + (i ∂ – q Φ*)Ψ₂ = 0 
(i ∂* – q Φ)Ψ₁ + (–m)Ψ₂ = 0
solve Ψ₁ from first equation
Ψ₁ = (i ∂ – q Φ*) Ψ₂ / m 
plug into the second equation 
–((∂* + i q Φ) (∂ + i q Φ*) + m²) Ψ₂ = 0 
mass conjugate operator
D' = W( i ∂ – q Φ* ) + m 𝐈
0 = D' D Ψ 
= (W( i ∂ – q Φ* )² – m² 𝐈) Ψ 
= –diag( 
     (∂ + i q Φ*)(∂* + i q Φ) + m² ,
     (∂* + i q Φ)(∂ + i q Φ*) + m²
    ) Ψ
0 = ( (∂ + i q Φ*)(∂* + i q Φ) + m² ) Ψ₁ 
= (–E² + P² + m² + i q F*) Ψ₁ 
0 = ( (∂* + i q Φ)(∂ + i q Φ*) + m² ) Ψ₂ 
= (–E² + P² + m² + i q F) Ψ₂
______________________________
gauge transform
λ = real scalar
Φ' = Φ + ∂* λ
∂ λ = (Φ' – Φ)*
Ψ = exp(i q λ) Ψ'
∂ exp( i q λ ) Ψ'
= exp(i q λ) ( i q ( ∂λ ) + ∂ ) Ψ'
0 = (W(i ∂ − q Φ*) − m I) Ψ  
= (W(i ∂ − q Φ*) − m I) exp(i q λ) Ψ'
= exp(i q λ) (W(i ∂ –q ( ∂λ ) − q Φ*) − m I) Ψ'
= exp(i q λ) (W(i ∂ –q (Φ' – Φ)*− q Φ*) − m I) Ψ'
= exp(i q λ) (W(i ∂ –q Φ'*) − m I) Ψ'
(W(i ∂ − q Φ'*) − m I) Ψ' = 0
______________________________
Dirac equation components
∂ = ∂ₜ + σ • ∂ᵣ
Φ = V + A • σ  ∈ ℝ+ℝ³
ψ = ψₛ + ψᵥ • σ  ∈ ℂ+ℂ³
F = Fᵥ • σ
    Fₛ = 0
    Fᵥ = 𝓔 + i 𝓑 ∈ ℂ³
∂* = ∂ₜ – σ • ∂ᵣ
Φ* = V – A • σ
F* = −Fᵥ* • σ
components of ∂ ψ
    scalar: ∂ₜ ψₛ + ∂ᵣ • ψᵥ
    vector: ∂ₜ ψᵥ + ∂ᵣ ψₛ + i (∂ᵣ × ψᵥ)
components of ∂* ψ
    scalar: ∂ₜ ψₛ − ∂ᵣ • ψᵥ
    vector: ∂ₜ ψᵥ − ∂ᵣ ψₛ − i (∂ᵣ × ψᵥ)
components of Φ ψ
    scalar: V ψₛ + A • ψᵥ
    vector: V ψᵥ + ψₛ A + i (A × ψᵥ)
components of Φ* ψ
    scalar: V ψₛ − A • ψᵥ
    vector: V ψᵥ − ψₛ A − i (A × ψᵥ)
components of i q F ψ 
    scalar: i q ( Fᵥ • ψᵥ )
    vector: i q ψₛ Fᵥ − q ( Fᵥ × ψᵥ )
components of i q F* ψ
    scalar: − i q ( Fᵥ* • ψᵥ )
    vector: − i q ψₛ Fᵥ* + q ( Fᵥ* × ψᵥ )
______________________________
low energy spin 0
E = i ∂ₜ – q V
P = –i ∂ᵣ – q A
rest mass oscillation
Ψ' = exp(–i m t) Ψ
E² Ψ' = exp(–i m t) (E + m)² Ψ
P² Ψ' = exp(–i m t) P² Ψ
0 = (–E² + P² + m²) Ψ'
= exp(–i m t) ( –(E + m)² + P² + m²) Ψ
0 = (–E² – 2 m E + P²) Ψ
recursively iterate E
E = ½ ( P² – E² ) / m
E = ½ P² / m – ⅛ P⁴ / m³ + ···
low energy approx
E Ψ ≈ ( ½ P² / m ) Ψ  
Schrödinger equation
i ∂ₜ Ψ = ( ½ (i ∂ᵣ + q A)² / m + q V) Ψ
if A = 0
    i ∂ₜ Ψ = ( –½ ∂ᵣ² / m + q V) Ψ
______________________________
low energy spin ½
Ψ' = {Ψ₁' ; Ψ₂'}
Dirac with potential field
(W( i ∂ – q Φ* ) – m 𝐈 )Ψ' = 0
2 equations
(i ∂ – q Φ*) Ψ₂' – m Ψ₁' = 0
(i ∂* – q Φ) Ψ₁' – m Ψ₂' = 0
∂ = ∂ₜ + σ • ∂ᵣ
Φ = V + A • σ
(i ∂ₜ – q V) Ψ₂' + ( (i ∂ᵣ + q A) • σ ) Ψ₂'
    – m Ψ₁' = 0
(i ∂ₜ – q V) Ψ₁'  + ( –(i ∂ᵣ + q A) • σ ) Ψ₁' 
    – m Ψ₂' = 0
E = i ∂ₜ – q V
P = –i ∂ᵣ – q A
2 equations
(E – P • σ) Ψ₂' – m Ψ₁' = 0
(E + P • σ) Ψ₁' – m Ψ₂' = 0
peel off rest mass
Ψₙ' = exp(–i m t) Ψₙ
E Ψₙ' = exp(–i m t) (E + m) Ψₙ
2 equations without rest mass
(E + m – P • σ) Ψ₂ – m Ψ₁ = 0
(E + m + P • σ) Ψ₁ – m Ψ₂ = 0
solve Ψ₁ from first equation
Ψ₁ = (E + m – P • σ) Ψ₂ / m
plug it into second equation
and multiply by m
0 = ((E + m + P • σ) (E + m – P • σ) – m²) Ψ₂
(E + m + P • σ) (E + m – P • σ) – m²
= (E + m)² – m² – (P • σ)² + [P • σ, E]
–(P • σ)² = –P² – i (P × P) • σ
[Pₘ, Pₙ] = [–i ∂ₘ –q Aₘ, –i ∂ₙ – q Aₙ ]
= i q (∂ₘ Aₙ – ∂ₙ Aₘ)
P × P = { [P₂, P₃], –[P₁, P₃], [P₁, P₂] }
= i q ∂ᵣ × A
= i q 𝓑
–i (P × P) • σ = q 𝓑 • σ
[P • σ, E + m] = [P, E] • σ
= i q ( ∂ᵣ V + ∂ₜ A ) • σ
= – i q 𝓔 • σ
– i q F = q 𝓑 • σ – i q 𝓔 • σ
exact solution
0 = ( (E + m)² – m² – P² + q 𝓑 • σ 
    + [P • σ, E] ) Ψ₂
= ( E² + 2 m E – P² + q 𝓑 • σ – i q 𝓔 • σ) Ψ₂
= ( E² + 2 m E – P² – i q F) Ψ₂
recursively iterate E
E = ( P² – E² – q 𝓑 • σ – [P • σ, E]) / (2m)
= ½ (P² – q 𝓑 • σ) / m + O(1/m² )
low energy approximation
0 ≈ ( 2 m E – P² + q 𝓑 • σ ) Ψ₂
Pauli equation
i ∂ₜ Ψ₂ = ( ( –i ∂ᵣ – q A)² / (2 m) + q V ) Ψ₂
            – q (𝓑 • σ) / (2 m) Ψ₂
______________________________
spin - rotors
i = σ₁σ₂σ₃
i σₖ = σₙσₘ  n≠m≠k
unit vectors
n = (n₁, n₂, n₃),  n² = 1
u = (u₁, u₂, u₃),  u² = 1
R = exp( ½ϕ (n • iσ) )  ←  rotor
= cos( ½ϕ ) + sin( ½ϕ ) (n • iσ)
note:  R(φ + 2π) = −R(φ)  
Rᴿ = cos( ½ϕ ) – sin( ½ϕ ) (n • iσ)  
Rᴿ R = cos( ½ϕ )² + sin( ½ϕ )² = 1
→ R⁻¹ = Rᴿ = R*ᴿ
Rᴿ (u • σ) R = u' • σ
Rodrigues’ formula
u' = u cos(ϕ) + (n × u) sin(ϕ) 
    + (n • u) (1 − cos(ϕ)) n
ψ = ψₛ + ψᵥ • σ
    ψₛ = scalar(ψ) ←  complex
    ψᵥ = vector(ψ) ←  complex
ψᴿ = ψₛᴿ − ψᵥᴿ • σ
   = ψₛ* − ψᵥ* • σ
ψ' = R ψ  ←  spinor
note: ( )ᴿ is the Hermitian operator
⟨ψ' |u • σ | ψ'⟩ = ⟨R ψ |u • σ | R ψ⟩ 
= ψᴿ (Rᴿ (u • σ) R) ψ
= ψᴿ (u' • σ) ψ
= ⟨ψ |u' • σ | ψ⟩ 
ψψᴿ = scalar(ψψᴿ) + vector(ψψᴿ) • σ
scalar(ψᴿ (u' • σ) ψ) 
= scalar((u' • σ) ψ ψᴿ) ← cyclic property
= u'• vector(ψ ψᴿ)
= u • vector(ψ'ψ'ᴿ)
______________________________
spin - projectors
projection operator
Π(±) = ½ ( 1 ± u • σ )
Π(±)² = Π(±)
Π(±) Π(∓) = 0
Π(±) + Π(∓) = 1
(u • σ) Π(±) = ±Π(±)    ←  ±1 eigenvalue
Π(±)ᴿ = Π(±)
Π(±)* = Π(∓)
u • σ = (+1) Π(+) + (−1) Π(–)
Π(±)' = ½ ( 1 ± u' • σ )
u' • σ = Rᴿ (u • σ) R
 = (+1) Π(+)' + (−1) Π(–)'
Ψ(±) = Π(±)Ψ
(u • σ) Ψ(±) = ±Ψ(±)
Ψ = Ψ(+) + Ψ(–)  ← spin up + spin down
Ψ(±)ᴿ = Ψᴿ Π(±)ᴿ
= Ψᴿ Π(±)
Ψ(±)ᴿ Ψ(±) = Ψ*ᴿ Π(±)Π(±) Ψ = Ψ*ᴿ Π(±) Ψ
Ψ(±)ᴿ Ψ(∓) = Ψ*ᴿ Π(±)Π(∓) Ψ = 0
Ψ(∓)ᴿ Ψ(±) = Ψ*ᴿ Π(∓)Π(±) Ψ = 0
⟨ψ |u' • σ | ψ⟩ = ⟨Ψ₊ + Ψ₋  |u' • σ | Ψ₊ + Ψ₋⟩
= Ψ(+)ᴿ (u' • σ) Ψ(+) + Ψ(+)ᴿ (u' • σ) Ψ(–)
    + Ψ(–)ᴿ(u' • σ) Ψ(+) + Ψ(–)ᴿ (u' • σ) Ψ(–)
( u' • σ )( u • σ ) = u' • u + i (u' × u) • σ
( u • σ )( u' • σ ) = u' • u – i (u' × u) • σ
( u • σ )( u' • σ )( u • σ ) = ( 2 ( u' • u ) u – u' ) • σ
par( u' ) := (u • u') u
prp( u' ) := u' − par( u' )
same–side projector
Π(±) ( u' • σ ) Π(±) = ± (u • u') Π(±)
cross–side projector
Π(±) ( u' • σ ) Π(∓)
= ½ ( [ prp( u' ) ± i ( u × u' ) ] • σ ) Π(∓)
u' • u =: cos(θ)
= (u • u) cos(ϕ) + (u • (n × u)) sin(ϕ) 
   + (n • u) (1 − cos(ϕ)) (n • u)
= (n • u)² + (1 − (n • u)²) cos(ϕ)
= cos(ϕ)  if  n • u = 0
½ ( 1 + u' • u ) = ½ (1 + cos(θ)) 
= cos²(θ/2)
= cos²(ϕ/2)  if  n • u = 0
½ ( 1 – u' • u ) = ½ (1 – cos(θ)) 
= sin²(θ/2)
= sin²(ϕ/2)  if  n • u = 0
______________________________
spin uncertainty
for u, u' ∈ ℝ³
u • u = u' • u' = 1
(u • σ) (u' • σ) = u • u' + i ( u × u' ) • σ
[ u • σ , u' • σ ] = 2 i ( u × u' ) • σ
(u • σ) (u' • σ) = u • u' + ½ [ u • σ , u' • σ ] 
using paravector metric ⟪Z⟫² = Z Z*ᴿ
⟪(u • σ) (u' • σ)⟫² = (u • u')² – (i ( u × u' ))²
= (u • u')² + ( u × u' )²
≥ ( u × u' )²
= ‖½ [ u • σ , u' • σ ] ‖²
Robertson inequality
⟪(du • σ) (du' • σ) ⟫ ≥ ½ ‖[ du • σ ,du' • σ ]‖
______________________________
raising and lowering about z
(ladder operators)
u = e₃ = (0, 0, 1)
Π(±) = ½ ( 1 ± σ₃ )
a(±) := ½ ( σ₁ ± i σ₂ ) ← ladder operator
commutator
[ σ₃ , a(±) ] = ± 2 a(±)
projector products
a(±) a(∓) = Π(±)
nilpotence
a(±)² = 0
action with projectors
a(+) Π(+) = 0, 
a(+) Π(–) = a(+)
Π(+) a(+) = a(+)
Π(–) a(+) = 0
a(–) Π(–) = 0
a(–) Π(+) = a(–)
Π(–) a(–) = a(–)
Π(+) a(–) = 0
reconstruction
σ₁ = a(+) + a(–)
σ₂ = ( a(+) − a(–) ) / i
σ₃ = Π(+) − Π(–)
______________________________
spin up and spin down
let
u = 𝓑 / ‖𝓑‖
projection operators
Π(±) = ½ ( 1 ± u • σ )
Π(±)² = Π(±)
Π(±) Π(∓) = 0
Π(±) + Π(∓) = 1
(u • σ) Π(±) = ±Π(±)    ←  ±1 eigenvalue
Π(±)ᴿ = Π(±)
Π(±)* = Π(∓)
u • σ = (+1) Π(+) + (−1) Π(–)
Ψ(±) = Π(±) Ψ
(u • σ) Ψ(±)  = ±Ψ(±)
Ψ = Ψ(+) + Ψ(–)  ← spin up + spin down
Pauli equation
P = –i ∂ᵣ – q A
i ∂ₜ Ψ = ( P² / (2 m) + q V ) Ψ
             – q ‖𝓑‖ (u • σ) / (2 m) Ψ
i ∂ₜ Ψ = H Ψ
i ∂ₜ Ψ(±) = i (∂ₜ Π(±)) Ψ + Π(±) (i ∂ₜ Ψ)
= i (∂ₜ Π(±)) Ψ + Π(±) (H Ψ)
add zero
Π(±) (H Ψ) = Π(±) (H Ψ) + ( H Π(±) Ψ – H Π(±) Ψ)
= H Ψ(±) + [Π(±), H] Ψ
i ∂ₜ Ψ(±) =  H Ψ(±) + ( i (∂ₜ Π(±)) + [Π(±), H]) Ψ
=: H Ψ(±) + S(±)  ←  coupled
i ∂ₜ Ψ(±) = ( ½(1/m) P² + q V ) Ψ(±)
                  –½ (q/m) ‖𝓑‖ (u • σ) Ψ(±) + S(±)
S(±) = ( (i ∂ₜ Π(±)) + [P², Π(±)] / (2m) ) Ψ
if uniform static 𝓑 
S(+) = S(–) = 0  ←  decoupled
(–𝓑 • σ) Ψ(±) = ∓‖𝓑‖ Ψ(±)
Zeeman effect 
i ∂ₜ Ψ(±) = ( ½ ( P² ∓ q ‖𝓑‖ ) / m + q V ) Ψ(±)
______________________________
4-state solution to Dirac equation
unit momentum
U = dX/ds
= ( E + P • σ ) / m
on mass shell
 ⟪U⟫² = (E² − P²) / m² = 1
W(U)² = diag( ⟪U⟫² , ⟪U⟫² ) = 𝐈
energy projector
    positive energy ← matter
    negative energy ← antimatter 
Λ(±) := ½ ( 𝐈 ± W(U) )
properties
    Λ(±)² = Λ(±)
    Λ(+) Λ(−) = 0
    Λ(+) + Λ(−) = 𝐈
spin projector
Π(±) = ½ ( 1 ± u • σ )
unit momentum operator
U → i ∂ / m
( W(U) – 𝐈 ) Ψ = 0  ←  normalized Dirac
solution
Ψ = Ψ(+,+) + Ψ(+,–) + Ψ(–,+) + Ψ(–,–)
where
    Ψ(ε, s) := Λ(ε) Π(s) Ψ
    ε ∈ {+, −}     ← energy sector
    s ∈ {+, −}     ← spin along u 
______________________________
transforms
X' = T⁻¹ X T
= T⁻¹ ( t + r • σ ) T
= t + r • ( T⁻¹σ T )
= t + r • σ'
T⁻¹ = Tᴿ* / (Tᴿ* T)
if boost or rotor,  
    Tᴿ* T = 1  and  T⁻¹ = Tᴿ* 
T⁻¹ T = 1
d(T⁻¹) T = –T⁻¹ dT
dX' = d(T⁻¹) X T + T⁻¹ dX T + T⁻¹ X dT
= T⁻¹ dX T 
    + (–T⁻¹ dT) (T T⁻¹) X T 
    + T⁻¹ X (T T⁻¹) dT
= T⁻¹ dX T + (–T⁻¹ dT) X' + X' (T⁻¹ dT)
= T⁻¹ dX T + [X', T⁻¹ dT]
dΘ := T⁻¹ dT
dX' = T⁻¹ dX T + [ X', dΘ ]
ddX' = d(T⁻¹) dX T + T⁻¹ ddX T + T⁻¹ dX dT
           + [dX', dΘ] + [X', ddΘ]
d(T⁻¹) dX T = (–T⁻¹ dT) (T⁻¹ dX T)
= –dΘ (dX' – [ X', dΘ ])
T⁻¹ dX dT = (T⁻¹ dX T) (T⁻¹ dT)
= (dX' – [ X', dΘ ]) dΘ
ddX' = T⁻¹ ddX T 
            +2 [ dX', dΘ ] 
            + [ X', ddΘ ]
            + [ dΘ, [ X', dΘ ] ]
dΘ/dt = Ω
ddΘ/dt² = dΩ/dt
dX'/dt = T⁻¹ (dX/dt) T + [ X', Ω ] 
ddX'/dt² = T⁻¹ (ddX/dt²) T 
                  +2 [ dX'/dt, Ω ] 
                  + [ X', dΩ/dt ]
                  + [ Ω, [ X', Ω ] ] 
______________________________
transformed basis
T(x) = f(x) + g(x) • σ
T⁻¹ = ( f – g • σ ) / ( f² – g² )
dT = df + dg • σ
σₙ' := T⁻¹ σₙ T
σₙ = T σₙ' T⁻¹
conjugate and reverse
    ( σₙ' )* = (T⁻¹)* σₙ* T*
    = –(T*)⁻¹ σₙ T*
    ( σₙ' )ᴿ = Tᴿ σₙᴿ (T⁻¹)ᴿ
    = Tᴿ σₙ (Tᴿ)⁻¹
    ( σₙ' )*ᴿ = –T*ᴿ σₙ (T⁻¹)*ᴿ
    = –T⁻¹σₙ T  if  T⁻¹ = T*ᴿ
    = –σₙ'          if  T⁻¹ = T*ᴿ
σ' algebra = σ algebra
    (σₙ')² = T⁻¹σₙTT⁻¹σₙT 
     = T⁻¹σₙ² T = 1
    if  n ≠ m
    T⁻¹σₙσₘT = T⁻¹σₙTT⁻¹σₘT = σₙ'σₘ'
         = T⁻¹εₙₘₖ i σₖT = εₙₘₖ i σₖ'
    σ₁'σ₂'σ₃' = T⁻¹σ₁TT⁻¹σ₂TT⁻¹σ₃T 
    = T⁻¹σ₁σ₂σ₃ T = i T⁻¹T = i
dσₙ' = [σₙ', T⁻¹ dT]
dσ' = [σ', T⁻¹ dT]  ←  vector
Z = Zₛ + Zᵥ • σ'  ←  paravector
⟪Z⟫² = ZZ*ᴿ
= (Zₛ + Zᵥ • σ' )(Zₛ – Zᵥ • σ' )  if  T⁻¹ = T*ᴿ
= Zₛ² – Zᵥ²  if  T⁻¹ = T*ᴿ
dZ = dZₛ + dZᵥ • σ' + Zᵥ • dσ'
= dZₛ + dZᵥ • σ' + Zᵥ • [σ', T⁻¹ dT]  
T⁻¹ dT = (f df – g dg) / ( f² – g² )
+ (f dg – df g – i g × dg) • σ / ( f² – g² )
define
Ω := T⁻¹ (dT/dt)
Ω' := T Ω T⁻¹ = (dT/dt) T⁻¹
scalar(σₙ') = scalar(T⁻¹ σₙ T)
= scalar(σₙ TT⁻¹) ←  cyclic property
= scalar(σₙ)
= 0
scalar(Ω σₙ') 
= scalar(T⁻¹ (dT/dt) T⁻¹ σₙ T)
= scalar((dT/dt) T⁻¹ σₙ TT⁻¹) ←  cyclic
= scalar((dT/dt) T⁻¹ σₙ )
= scalar( Ω' σₙ )
kₙ + i ωₙ := scalar( 2 Ω' σₙ )
Ω = ½ (k + i ω) • σ' + scalar
[σ', ω • ( i σ' )] = –2 ω ×σ'
[σ', k • σ'] = 2 i k × σ'
dσ'/dt = [σ', Ω] 
= i ( k + i ω ) × σ'
= ( –ω + i k ) × σ'
Zᵥ • [σ', Ω] = ( –Zᵥ × ω + i Zᵥ × k ) • σ'
dZ/dt  = dZₛ/dt
+ (dZᵥ/dt – Zᵥ × ω + i Zᵥ × k ) • σ'
e.g., non-inertial frame
X = t + r • σ' 
dX/dt  = dt/dt
          + (dr/dt – r × ω + i r × k ) • σ'
energy & momentum
m dX/dt  = E + P' • σ'
E  / m = dt/dt
P' / m = dr/dt + r × (–ω + i k)
if  T⁻¹ = T*ᴿ,
    ⟪dX⟫² = dt² – (dr – r × ω dt + i r × k dt)²
    = (1 – (dr/dr – r × ω + i r × k)²) dt²
    ds' := √(1 – (v – r × ω + i r × k)²) dt
______________________________
Lorentz boost ( T = L )
θ = arctanh(β)
β = tanh(θ)
ɣ = cosh(θ)
ɣ β = sinh(θ)
u = (u₁, u₂, u₃),  u² = 1
L(θ,u) := exp( ½ θ (u • σ) )
= cosh( ½ θ) + sinh( ½ θ ) (u • σ)
L⁻¹ = exp( −½ θ (u • σ) )
note: L⁻¹ = L* = L*ᴿ
commutative over the same axis
[ L(θ,u), L(θ',u) ] = 0
co-linear boost
L(θ', u) = L(θ₁, u) L(θ₂, u)
β' = (β₁ + β₂) / (1 + β₁ β₂)
ɣ' = ɣ₁ ɣ₂ (1 + β₁ β₂)
re-mapping space-time numbers
X  = t + r • σ
X' = L* X L
   = t' + r' • σ
t' = t cosh(θ) − sinh(θ) (u • r)
= ɣ ( t − β (u • r) )
split any vector r into components 
parallel and perpendicular to u
par(r) := (r • u) u  ← parallel to u
prp(r) := r − par(r)   ← perpendicular to u
r = par( r ) + prp( r )
par( r' ) = ɣ ( par( r ) - β t u)
prp( r' ) = prp( r )
r' = par( r' ) + prp( r ' )
= prp( r ) + ɣ ( par( r ) – β t u)
= r +( ɣ ( r • u – β t ) − r • u ) u
if β → 0,
ɣ → 1, t' → t, r' → r
paravector metric invariance
⟪X'⟫² = X' X'* 
= (L* X L) (L* X L)*
= L* X L L* X* L
= L* X X* L
= ⟪X⟫² L* L
= ⟪X⟫²
______________________________
Lorentz transform of fields
F' = L* F L
L = exp( ½ θ u • σ ),   u² = 1
β = tanh(θ),   ɣ = cosh(θ)
F = ( 𝓔 + i 𝓑 ) • σ
par( 𝓔 ' ) = par( 𝓔 )
par( 𝓑 ' ) = par( 𝓑 )
prp( 𝓔 ' ) = ɣ ( prp( 𝓔 ) + β u × 𝓑 )
prp( 𝓑 ') = ɣ ( prp( 𝓑 ) − β u × 𝓔 )
in one line
𝓔 ' = par(𝓔) + ɣ ( prp(𝓔) + β u × 𝓑 )
𝓑 ' = par(𝓑) + ɣ ( prp(𝓑) − β u × 𝓔 )
∂' = L* ∂ L
∂' F' = L* (ρ − J • σ) L
______________________________
rotation ( T = R )
R(θ,u) := exp( ½ θ (u • i σ) )
= cos( ½ θ) + sin( ½ θ ) (u • i σ)  ← rotor
R⁻¹ = exp( –½ θ (u • i σ) )
note: R⁻¹ = Rᴿ = R*ᴿ
X' = Rᴿ X R
   = t' + r' • σ
t' = t
r' = r cos(θ) + (u × r) sin(θ) 
    + (u • r) (1 − cos(θ)) u 
    ← Rodrigues’s formula
differentials
d(Rᴿ) = (dR)ᴿ = dRᴿ
Rᴿ R = 1
dRᴿ R + Rᴿ dR = 0
dR = –R dRᴿ R
reverse symmetry 
Ω := –R (dRᴿ/ds) = (dR/ds) Rᴿ
Ωᴿ = –Ω
Ω' = Rᴿ Ω R = R (dR/ds)
Ω'ᴿ = –Ω'
dR/ds = Ω R
X' = Rᴿ X R  →  R X' = X R
dR X' + R dX' = dX R + X dR
R dX' = dX R + X dR – dR X' 
 = dX R + X dR – dR X'
dX'/ds = Rᴿ (dX/ds) R + [ X', Ω' ]
= Rᴿ ( dX/ds + [ X, Ω ] ) R
______________________________
compositions
T = L R
T⁻¹ = T*ᴿ = R*ᴿ L*ᴿ = Rᴿ L*
T = R L
T⁻¹ = T*ᴿ = L*ᴿ R*ᴿ = L* Rᴿ
T = R₁ R₂
T⁻¹ = T*ᴿ = R₂ᴿ R₁ᴿ
T = L₁ L₂
T⁻¹ = T*ᴿ = L₂*L₁*
______________________________
quaternions
q₁ = –σ₂ σ₃
q₂ = –σ₃ σ₁
q₃ = –σ₁ σ₂
q₁ q₂ q₃ = –1
= q₁² 
= q₂² 
= q₃²
vectors
q = (q₁, q₂, q₃)
b = (b₁, b₂, b₃)
q = –i σ
q* = q
qᴿ = –q
Q = a + b • q 
= a – i b • σ
Q* = Q
Q Q* = Q²
= a² –2 i (b • σ) + (–i)² (b • σ)²
= a² – b² –2 i a b • σ 
Qᴿ = a + i b • σ
Q Qᴿ = a² + b²
rotation
r' • q = Q (r • q) Qᴿ  ← see R(θ,u)
______________________________
functions of u • σ,  u² = 1
let
Π(±) = ½ ( 1 ± u • σ )
Π(±)² = Π(±)
Π(±) Π(∓) = 0
Π(±) + Π(∓) = 1
(u • σ) Π(±) = ±Π(±)
if f( ) has a power series
f( u • σ ) 
= ½ ( f(+1) + f(−1) ) 
+ ½ ( f(+1) − f(−1) )  ( u • σ )
= f(+1) Π(+) + f(−1) Π(–)
f( u • σ ) = f(+1) Π(+) + f(−1) Π(−)
√( u • σ ) = ± Π(+) ± i Π(–)
log( u • σ ) = 2 i π ( n Π(+) + (½ + m) Π(–) )
(z + u • σ )⁻¹ = ( z − u • σ ) / (z² − 1)
( i – u • σ ) / ( i + u • σ ) = i u • σ
exp( θ u • σ ) = cosh(θ) + sinh(θ) (u • σ)
exp( i θ u • σ ) = cos(θ) + i sin(θ) (u • σ)
______________________________
functions of X
X = t + r • σ 
t = scalar( X ) ∈ ℝ
r = vector( X ) ∈ ℝ³
u = r / ‖r‖
minimal polynomial
X² − 2 t X + ( t² − r² ) = 0 
    λ(±) = t ± ‖r‖
    0 = ( X − λ(–) ) ( X − λ(+) )
    = X² – λ(–) X – X λ(+) + λ(–) λ(+)
a  = λ(+)
b  = λ(–)
Π(±) = ½ ( 1 ± u • σ )
X = a Π(+) + b Π(–)
if f( ) has a power series,
X' = f(X) 
= f(a) Π(+) + f(b) Π(–)
= t' + r' • σ
t' = ½ ( f(a) + f(b) )
‖r'‖ = ½ ( f(a) – f(b) )
r' = ‖r'‖ u
summary left as projectors
    eigenvalues
    λ(±) = t ± ‖r‖
    eigen-projectors
    u := r / ‖r‖
    Π(±) = ½ ( 1 ± u • σ )
    eigenvalue equation
    Π(±) X = λ(±) Π(±)
    projector decomposition
    X = λ(+) Π(+) + λ(–) Π(–)
    f(X) = f(λ(+)) Π(+) + f(λ(–)) Π(–)
exp(X) = eᵗ ( cosh(‖r‖) + sinh(‖r‖) ( u • σ ) ) 
√X = ±√a Π(+) ± √b Π(–)
log(X) = log( a ) Π(+) + log( b ) Π(–)
X⁻¹ = Π(+) / a + Π(–) / b
Xⁿ = ½ (aⁿ + bⁿ) + ½ (aⁿ - bⁿ) ( u • σ )
product of any two functions
f(X) g(X') = (f(a) Π(+) + f(b) Π(–))(g(a') Π(+)' + g(b') Π(–)')
= f(a) g(a') Π(+)Π(+)' +f(a) g(b') Π(+)Π(–)' 
  + f(b) g(a') Π(–)Π(+)' +f(b) g(b') Π(–)Π(–)'
Π(+) Π(+)' = ¼ (  1 + u • u' ) 
    + ¼ ( u' + u + i (u × u')) • σ
Π(+) Π(–)' = ¼ ( 1 – u • u' ) 
    + ¼ ( –u' + u – i (u × u')) • σ
Π(–) Π(+)' = ¼ ( 1 – u • u' ) 
    + ¼ ( u' –u – i (u × u')) • σ
Π(–) Π(–)' = ¼ ( 1 + u • u' ) 
    + ¼ ( –u' –u + i (u × u')) • σ
( 1 – r • σ ) / ( 1 + r • σ )
= exp(–2 arctanh( r • σ ))
= exp(–2 arctanh( ‖r‖ )) Π(+) 
+ exp(+2 arctanh( ‖r‖ )) Π(–) 
______________________________
generalized spacetime
X = f(t, r) + g(t, r) • σ
dX = df + dg • σ
= (∂ₜf + ∂ₜg • σ) dt + (∂ᵣf + ∂ᵣ(g • σ)) • dr
=: A dt + B • dr
dX/dt = (∂ₜf + ∂ₜg • σ) + (∂ᵣf + ∂ᵣ(g • σ)) • dr/dt
= A + B • dr/dt
dX* = df – dg • σ
= (∂ₜf – ∂ₜg • σ) dt + (∂ᵣf – ∂ᵣ(g • σ)) • dr
=: A* dt + B* • dr
dX*/dt = (∂ₜf – ∂ₜg • σ) + (∂ᵣf – ∂ᵣ(g • σ)) • dr/dt
= A* + B* • dr/dt
⟪dX⟫² = (dX/dt) (dX*/dt) dt²
= (A + B • dr/dt)(A* + B* • dr/dt) dt²
= AA* dt² + (AB* + BA*) • dr dt 
    + (B • dr)(B* • dr)
= df² – dg²
ds = √( df² – dg² ) 
= √( 1 – (dg/df)² ) df
relativistic velocity
U = dX/ds 
= ( df + dg • σ ) / √( df² – dg² )
= ( 1 + dg/df • σ ) / √( 1 – (dg/df)² )
= ( 1 + v • σ ) / √( 1 – v² )
= ɣ ( 1 + v • σ )
geometric velocity
v = dg/df
= (dg/dt) / (df/dt)
= (dg/ds) / (df/ds)
geometric Lorentz factor
ɣ = df/ds
= df / √( df² – dg² ) 
= (df/dt) / (ds/dt)
= 1 / √( 1 – v² )
note: to remain time-like, v² < 1
note: 4 momentum equation is unchanged
geometric acceleration
a = d/ds v 
= d/ds (dg/df)
= d/df (df/ds) (dg/df)
= ɣ d/df (dg/df)
= ɣ d²g/df²
relativistic acceleration
dU/ds = ɣ⁴ (v • a) + (ɣ⁴ (v • a) v + ɣ² a) • σ
let 
f' = df/dt, g' = dg/dt
f'' = df'/dt, g''= dg'/dt
dv/dt = d/dt ( g' / f' )
= ( g'' f' – g' f'' ) / ( f' )²
    consider a plane curve (x(t), y(t))
    k = |x' y'' – y' x''| / √( x'² + y'² )³
[ U, dU/ds ] = [ U, dU/dt] (dt/ds)
= ɣ [ U, dU/dt ]
= 2 ɣ³ ( v × a ) • ( i σ )
inertial motion: 
no change in energy or momentum
m dU/ds = 0  ←  geodesic equation
i.e., ddX/ds² = 0
scalar component: ɣ³ (v • a) = 0
    → v is perpendicular to a, or a = 0
vector component: ɣ³ (v • a) v + ɣ a = 0
    in matrix form
    M = 𝐈 + ɣ² v vᵀ
    = 𝐈 + v vᵀ / ( 1 – v² ),  v² < 1  
    →  positive definite
    → ɣ M • a = 0  if and only if a = 0
    → d²g/df² = 0  →  g(t, r) = c₁ f(t, r) + c₂
    → linear in (f, g) coordinates
    → v = dg/df → c₁ = v = constant
______________________________
metric
x = {t, r₁, r₂, r₃}
ds² = ⟪dX⟫² 
= dxᵀ • Q(x) • dx  ←  matrix product
eigenvalue decomposition
Q = Qᵀ
det( Q − λₙ 𝐈 ) = 0
( Q − λₙ ​𝐈 ) uₙ ​= 0
D = diag {λ₀​, λ₁​, λ₂​, λ₃​} 
    ← four real eigenvalues
N = {u₀, u₁, u₂, u₃}  
    ← orthonormal matrix
NᵀN = NNᵀ = 𝐈  →  N⁻¹ = Nᵀ
Q⁻¹ = N D⁻¹ Nᵀ
Q = N D Nᵀ
= Σₙ λₙ uₙ uₙᵀ
ds² = dxᵀ • Q • dx 
= Σₙ λₙ (uₙ • dx)²
dxₙ' = √λₙ ( uₙ • dx )
ds² = dx' • dx'
______________________________
matrix form
σ₁ = { 0, 1; 1, 0 }
σ₂ = { 0, –i;  i, 0 }
σ₃ = { 1, 0; 0, –1 }
σₙ² = 𝐈
σₘσₙ = –σₙσₘ  if  n ≠ m
σ₁σ₂σ₃ = i 𝐈
spacetime
  X = t + r₁σ₁+ r₂σ₂ + r₃σ₃
  = { t + r₃, r₁ – i r₂; r₁ + i r₂, t – r₃}
  tr(X) = 2 t
  det(X) = t² – r² 
  λ² − tr(X) λ + det(X) = 0  if  λ = t ± ‖r‖
complex conjugate
  ( i )ᶜ = –i
  σ₁ᶜ = +σ₁
  σ₂ᶜ = –σ₂
  σ₃ᶜ = +σ₃
complex conjugate transpose
  σₙᶜᵀ = σₙ
  for n ≠ m
  (σₙσₘ)ᶜᵀ = σₘᶜᵀσₙᶜᵀ 
  = σₘσₙ 
  = –σₙσₘ
  (i 𝐈)ᶜᵀ = –i 𝐈
  ((σ₁σ₂σ₃)σₙ)ᶜᵀ = σₙᶜᵀ (σ₁σ₂σ₃)ᶜᵀ 
  = σₙ (i 𝐈)ᶜᵀ = –i σₙ
reverse operation
  𝐈ᴿ = 𝐈
  σₙᴿ = σₙ
  (i 𝐈)ᴿ = –i 𝐈
  (i σₙ)ᴿ = –i σₙ
  ( )ᴿ = ( )ᶜᵀ  ←  Hermitian adjoint 
sigma conjugation
  𝐈* = 𝐈
  σₙ* = –σₙ
  (i 𝐈)* =  –i 𝐈
  (i σₙ)* = i σₙ
  M = { 0, 1; –1 , 0 }
  M⁻¹ = –M
  ( )* = M ( )ᶜ M⁻¹
reverse sigma conjugation
  ( )*ᴿ = ( )ᴿ* = M ( )ᵀ M⁻¹