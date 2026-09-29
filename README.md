# Space–Time Algebra

**[Complex scalars](https://en.wikipedia.org/wiki/Complex_number) and [three-vectors](https://en.wikipedia.org/wiki/Euclidean_vector) provide a common [paravector](https://en.wikipedia.org/wiki/Paravector) language for [quaternions](https://en.wikipedia.org/wiki/Quaternion), [rotations](https://en.wikipedia.org/wiki/Rotation_%28mathematics%29), [Lorentz boosts](https://en.wikipedia.org/wiki/Lorentz_transformation), [projections](https://en.wikipedia.org/wiki/Paravector#Null_paravectors_as_projectors), [spacetime geometry](https://en.wikipedia.org/wiki/Minkowski_space), [relativistic particle dynamics](https://en.wikipedia.org/wiki/Relativistic_mechanics), [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations), [electromagnetic waves](https://en.wikipedia.org/wiki/Electromagnetic_radiation) and [forces](https://en.wikipedia.org/wiki/Lorentz_force), [gauge coupling](https://en.wikipedia.org/wiki/Minimal_coupling), the [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) and [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) equations, their [Schrödinger](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation) and [Pauli](https://en.wikipedia.org/wiki/Pauli_equation) low-energy limits, and [spin-one-half](https://en.wikipedia.org/wiki/Spin-1/2) [amplitudes](https://en.wikipedia.org/wiki/Probability_amplitude) and [measurement probabilities](https://en.wikipedia.org/wiki/Born_rule).**

[Read the paper](./sta_notes.pdf) · [LaTeX source](./sta_notes.tex)

[Natural units](https://en.wikipedia.org/wiki/Natural_units) $`\hbar=c=1`$, [rationalized electromagnetic units](https://en.wikipedia.org/wiki/Heaviside%E2%80%93Lorentz_units), and [signature](https://en.wikipedia.org/wiki/Metric_signature) $`(+,-,-,-)`$.

## One complex [paravector](https://en.wikipedia.org/wiki/Paravector)

```math
\boxed{
\begin{aligned}
\sigma_n^2&=1,
&\sigma_n\sigma_m&=-\sigma_m\sigma_n\quad(n\ne m),\\
\mathrm{i}&:=\sigma_1\sigma_2\sigma_3,
&\mathrm{i}^2&=-1,\\
\boldsymbol{\sigma}&:=\{\sigma_1,\sigma_2,\sigma_3\},
&\mathrm{i}\boldsymbol{\sigma}&=\{\sigma_2\sigma_3,\sigma_3\sigma_1,\sigma_1\sigma_2\},\\
Z&:=(a+b\mathrm{i})+(\mathbf A+\mathrm{i}\mathbf B)\cdot\boldsymbol{\sigma},
&a,b&\in\mathbb R,\quad\mathbf A,\mathbf B\in\mathbb R^3.
\end{aligned}}
```

These eight real components form a complex scalar and a complex vector:

```math
Z=S+\mathbf V\cdot\boldsymbol{\sigma}=X+\mathrm{i}Y,
\qquad X:=a+\mathbf A\cdot\boldsymbol{\sigma},
\quad Y:=b+\mathbf B\cdot\boldsymbol{\sigma}.
```

The component extractions are

```math
\begin{aligned}
\mathrm{re}(Z)&:=X,\\
\mathrm{im}(Z)&:=Y,\\
\mathrm{sc}(Z)&:=S=a+\mathrm{i}b,\\
\mathrm{vec}(Z)&:=\mathbf V=\mathbf A+\mathrm{i}\mathbf B.
\end{aligned}
```

## Product

The product carries both the [dot](https://en.wikipedia.org/wiki/Dot_product) and [cross](https://en.wikipedia.org/wiki/Cross_product) products:

```math
(\mathbf A\cdot\boldsymbol{\sigma})(\mathbf B\cdot\boldsymbol{\sigma})
=\mathbf A\cdot\mathbf B+\mathrm{i}(\mathbf A\times\mathbf B)\cdot\boldsymbol{\sigma}.
```

Here we use the convention $`\mathbf A^2:=\mathbf A\cdot\mathbf A`$.

The [commutator](https://en.wikipedia.org/wiki/Commutator) is

```math
[A,B]:=AB-BA.
```

## Core operations

[Determinant](https://en.wikipedia.org/wiki/Determinant), [trace](https://en.wikipedia.org/wiki/Trace_%28linear_algebra%29), [squared norm](https://en.wikipedia.org/wiki/Norm_%28mathematics%29), adjugate, [inverse](https://en.wikipedia.org/wiki/Invertible_matrix), and the two conjugations below.

Write $`Z=S+\mathbf V\cdot\boldsymbol{\sigma}`$, where $`S\in\mathbb C`$ is a complex scalar and $`\mathbf V\in\mathbb C^3`$ is a complex vector.

$`(\,)^{*}:=`$ [complex conjugation](https://en.wikipedia.org/wiki/Complex_conjugate) of complex numbers and a sign flip of $`\boldsymbol{\sigma}`$.

$`(\,)^{\mathsf H}:=`$ [Hermitian conjugation](https://en.wikipedia.org/wiki/Conjugate_transpose), which flips the order of multiplied $`\sigma_k`$.

```math
\begin{gathered}
\begin{aligned}
Z^{*}&=S^{*}-\mathbf V^{*}\cdot\boldsymbol{\sigma},
&Z^{\mathsf H}&=S^{*}+\mathbf V^{*}\cdot\boldsymbol{\sigma},\\
\det Z&:=S^2-\mathbf V^2,
&\mathrm{adj}\,Z&:=S-\mathbf V\cdot\boldsymbol{\sigma} = Z^{*\mathsf H},\\
\mathrm{tr}\,Z&:=2S=Z+\mathrm{adj}\,Z,
&\lVert Z\rVert^2&:=SS^{*}+\mathbf V\cdot\mathbf V^{*} = \mathrm{sc}(Z Z^\mathsf H).
\end{aligned}\\[6pt]
Z^{-1}=\frac{\mathrm{adj}\,Z}{\det Z}\qquad(\det Z\ne0).
\end{gathered}
```

## [Quaternions](https://en.wikipedia.org/wiki/Quaternion)

In this representation, quaternions have a real scalar part and a purely imaginary vector part:

```math
\begin{aligned}
q_k&:=-\mathrm{i}\sigma_k,
\qquad q_1^2=q_2^2=q_3^2=q_1q_2q_3=-1,\\
Q&:=a+\mathbf A\cdot\mathbf q
=a-\mathrm{i}\mathbf A\cdot\boldsymbol{\sigma},
\qquad a\in\mathbb R,\quad\mathbf A\in\mathbb R^3,\\
Q^{\mathsf H}&=a+\mathrm{i}\mathbf A\cdot\boldsymbol{\sigma},\\
QQ^{\mathsf H}&=\det Q=a^2+\mathbf A^2.
\end{aligned}
```

## Rotations and boosts

For real spacetime $`X`$ and real [unit axis](https://en.wikipedia.org/wiki/Unit_vector) $`\mathbf u`$:

```math
X:=t+\mathbf r\cdot\boldsymbol{\sigma},
\qquad X':=TXT^{\mathsf H}=t'+\mathbf r'\cdot\boldsymbol{\sigma}.
```

For $`T=R`$: [Rodrigues rotation](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula) through $`\theta`$, with [unit quaternion](https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation) $`Q=R`$, $`q_k=-\mathrm{i}\sigma_k`$.

```math
\begin{aligned}
R&:=e^{-\mathrm{i}\theta\mathbf u\cdot\boldsymbol{\sigma}/2}
=\cos\frac\theta2-\mathrm{i}\mathbf u\cdot\boldsymbol{\sigma}\sin\frac\theta2,\\
t'&=t,\\
\mathbf r'&=\mathbf r\cos\theta+(\mathbf u\times\mathbf r)\sin\theta
+(\mathbf u\cdot\mathbf r)(1-\cos\theta)\mathbf u.
\end{aligned}
```

For $`T=L`$: [Lorentz boost](https://en.wikipedia.org/wiki/Lorentz_transformation) to a frame moving at $`+\beta\mathbf u`$, with [rapidity](https://en.wikipedia.org/wiki/Rapidity) $`\theta`$, $`\beta:=\tanh\theta`$, and [Lorentz factor](https://en.wikipedia.org/wiki/Lorentz_factor) $`\gamma:=\cosh\theta`$.

Define the parallel and perpendicular components:

```math
\mathbf r_\parallel:=(\mathbf u\cdot\mathbf r)\mathbf u,
\qquad \mathbf r_\perp:=\mathbf r-\mathbf r_\parallel.
```

```math
\begin{aligned}
L&:=e^{-\theta\mathbf u\cdot\boldsymbol{\sigma}/2}
=\cosh\frac\theta2-\mathbf u\cdot\boldsymbol{\sigma}\sinh\frac\theta2
=L^{\mathsf H},\\
t'&=\gamma(t-\beta\mathbf u\cdot\mathbf r),\\
\mathbf r'&=\mathbf r_\perp+\gamma(\mathbf r_\parallel-\beta t\mathbf u).
\end{aligned}
```

## Projections and spectral decomposition

```math
Z=a+b\mathbf u\cdot\boldsymbol{\sigma},
\qquad a,b\in\mathbb C,\quad \mathbf u\in\mathbb C^3,\quad \mathbf u^2=1.
```

The eigenvalues and complementary [projectors](https://en.wikipedia.org/wiki/Paravector#Null_paravectors_as_projectors) are

```math
\begin{gathered}
\lambda_\pm:=a\pm b,
\qquad \Pi_\pm:=\frac12(1\pm\mathbf u\cdot\boldsymbol{\sigma}),\\[4pt]
\Pi_\pm^2=\Pi_\pm,\qquad \Pi_+\Pi_-=0,\qquad \Pi_-+\Pi_+=1.
\end{gathered}
```

For $`f`$ analytic near $`\lambda_\pm`$:

```math
\begin{aligned}
Z&=\lambda_-\Pi_-+\lambda_+\Pi_+,\\[4pt]
f(Z)&=f(\lambda_-)\Pi_-+f(\lambda_+)\Pi_+.
\end{aligned}
```

## Spacetime from the [determinant](https://en.wikipedia.org/wiki/Determinant)

For a future-directed massive particle, [proper time](https://en.wikipedia.org/wiki/Proper_time) and [four-momentum](https://en.wikipedia.org/wiki/Four-momentum) follow from a real paravector, with [four-velocity](https://en.wikipedia.org/wiki/Four-velocity) $`U`$.

```math
\begin{aligned}
\mathrm dX&=\mathrm dt+\mathrm d\mathbf r\cdot\boldsymbol{\sigma},
&\mathrm ds^2&:=\det(\mathrm dX)=\mathrm dt^2-\mathrm d\mathbf r^2,\\
U&:=\frac{\mathrm dX}{\mathrm ds}=\gamma(1+\mathbf v\cdot\boldsymbol{\sigma}),
&\gamma&=(1-\mathbf v^2)^{-1/2},\quad \mathbf v:=\frac{\mathrm d\mathbf r}{\mathrm dt},\\
P&:=mU=\underbrace{E}_{\text{energy}}+
\underbrace{\mathbf p\cdot\boldsymbol{\sigma}}_{\text{momentum}},
&\det P&=E^2-\mathbf p^2=m^2.
\end{aligned}
```

## [Maxwell](https://en.wikipedia.org/wiki/Maxwell%27s_equations) in one equation

Combine the [electric](https://en.wikipedia.org/wiki/Electric_field) and [magnetic](https://en.wikipedia.org/wiki/Magnetic_field) fields and use the [paravector derivative](https://en.wikipedia.org/wiki/Paravector#Paragradient):

```math
F:=(\mathbf E+\mathrm{i}\mathbf B)\cdot\boldsymbol{\sigma},
\qquad \partial:=\partial_t+\partial_{\mathbf r}\cdot\boldsymbol{\sigma},
\qquad \boxed{\partial F=\rho-\mathbf J\cdot\boldsymbol{\sigma}.}
```

Here $`\rho`$ is [charge density](https://en.wikipedia.org/wiki/Charge_density), $`\mathbf J`$ is [current density](https://en.wikipedia.org/wiki/Current_density), and $`\partial_{\mathbf r}`$ is the [gradient](https://en.wikipedia.org/wiki/Gradient). Its [divergence](https://en.wikipedia.org/wiki/Divergence) and [curl](https://en.wikipedia.org/wiki/Curl_%28mathematics%29) give the four components of [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations):

```math
\begin{aligned}
\text{Real scalar:}\quad &\partial_{\mathbf r}\cdot\mathbf E=\rho,\\
\text{Imaginary scalar:}\quad &\partial_{\mathbf r}\cdot\mathbf B=0,\\
\text{Real vector:}\quad &\partial_t\mathbf E-\partial_{\mathbf r}\times\mathbf B=-\mathbf J,\\
\text{Imaginary vector:}\quad &\partial_t\mathbf B+\partial_{\mathbf r}\times\mathbf E=0.
\end{aligned}
```

Applying $`\partial^{*}=\partial_t-\partial_{\mathbf r}\cdot\boldsymbol{\sigma}`$ gives the [d'Alembertian](https://en.wikipedia.org/wiki/D%27Alembert_operator) and the [sourced wave equations](https://en.wikipedia.org/wiki/Inhomogeneous_electromagnetic_wave_equation):

```math
\begin{aligned}
\Box&:=\partial^{*}\partial=\partial_t^2-\partial_{\mathbf r}^2,\\
\partial^{*}(\partial F)&=\Box F=\partial^{*}(\rho-\mathbf J\cdot\boldsymbol{\sigma}).
\end{aligned}
```

Since $`\Box F=(\Box\mathbf E+\mathrm{i}\Box\mathbf B)\cdot\boldsymbol{\sigma}`$, its components give [charge conservation](https://en.wikipedia.org/wiki/Charge_conservation) and the field wave equations

```math
\begin{aligned}
\text{Real scalar:}\quad &0=\partial_t\rho+\partial_{\mathbf r}\cdot\mathbf J
&&\text{(charge conservation)},\\
\text{Imaginary scalar:}\quad &0=0,\\
\text{Real vector:}\quad &\Box\mathbf E=-\partial_t\mathbf J-\partial_{\mathbf r}\rho,\\
\text{Imaginary vector:}\quad &\Box\mathbf B=\partial_{\mathbf r}\times\mathbf J.
\end{aligned}
```

In vacuum, $`\Box F=0`$.

## [Electromagnetic potentials](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) and [gauge invariance](https://en.wikipedia.org/wiki/Electromagnetic_four-potential#Gauge_freedom)

Define the real [four-potential](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) from the [scalar potential](https://en.wikipedia.org/wiki/Electric_potential) $`V`$ and [vector potential](https://en.wikipedia.org/wiki/Magnetic_vector_potential) $`\mathbf A`$, with gauge scalar $`S`$:

```math
\Phi:=V+\mathbf A\cdot\boldsymbol{\sigma},
\qquad S:=\mathrm{sc}(\partial\Phi).
```

**Fields from the potential:**

```math
\boxed{\partial\Phi=S+F^{*}}
```

Its components give

```math
\begin{aligned}
\text{Real scalar:}\quad &S=\partial_tV+\partial_{\mathbf r}\cdot\mathbf A,\\
\text{Imaginary scalar:}\quad &0=0,\\
\text{Real vector:}\quad &\mathbf E=-\partial_t\mathbf A-\partial_{\mathbf r}V,\\
\text{Imaginary vector:}\quad &\mathbf B=\partial_{\mathbf r}\times\mathbf A.
\end{aligned}
```

**Source equations in any [gauge](https://en.wikipedia.org/wiki/Gauge_fixing):**

```math
\Box\Phi=\partial^{*}(\partial\Phi)
=\partial^{*}(S+F^{*})
=\partial^{*}S+(\partial F)^{*}.
```

Conjugating and using Maxwell’s equation gives

```math
\boxed{\Box\Phi^{*}-\partial S=\rho-\mathbf J\cdot\boldsymbol{\sigma}}
```

Its components give

```math
\begin{aligned}
\text{Real scalar:}\quad &\Box V-\partial_tS=\rho,\\
\text{Imaginary scalar:}\quad &0=0,\\
\text{Real vector:}\quad &\Box\mathbf A+\partial_{\mathbf r}S=\mathbf J,\\
\text{Imaginary vector:}\quad &\mathbf 0=\mathbf 0.
\end{aligned}
```

**[Gauge invariance](https://en.wikipedia.org/wiki/Electromagnetic_four-potential#Gauge_freedom):** for real scalar $`\lambda`$,

```math
\begin{aligned}
\Phi'&:=\Phi+\partial^{*}\lambda,\\
S'&=S+\Box\lambda,\\
F'&=F.
\end{aligned}
```

**[Lorenz gauge](https://en.wikipedia.org/wiki/Lorenz_gauge_condition):** $`S=0`$ gives the [potential wave equations](https://en.wikipedia.org/wiki/Inhomogeneous_electromagnetic_wave_equation)

```math
\begin{aligned}
\text{Real scalar:}\quad &\Box V=\rho,\\
\text{Real vector:}\quad &\Box\mathbf A=\mathbf J.
\end{aligned}
```

## [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) as a first-order wave equation

In the [Weyl representation](https://en.wikipedia.org/wiki/Gamma_matrices#Weyl_%28chiral%29_basis), write $`\boldsymbol\Psi:=(\Psi_1,\Psi_2)^{\mathsf T}`$, with two complex components in each entry. Let $`\mathbf I`$ be the block identity.

Hats mark named [energy](https://en.wikipedia.org/wiki/Energy_operator), [momentum](https://en.wikipedia.org/wiki/Momentum_operator), and [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) operators; the derivative symbols $`\partial`$ and $`\Box`$ remain unhatted.

For a free particle of mass $`m`$, define

```math
\begin{aligned}
\hat E&:=\mathrm{i}\partial_t,\\
\hat{\mathbf p}&:=-\mathrm{i}\partial_{\mathbf r},\\
\hat P&:=\hat E+\hat{\mathbf p}\cdot\boldsymbol{\sigma} = \mathrm{i}\partial^{*},\\
\mathbf W(\hat P)&:=
\begin{pmatrix}
0 & \mathrm{adj}\,\hat P \\
\hat P & 0
\end{pmatrix},\\
\hat D(m)&:=\mathbf W(\hat P)-m\mathbf I.
\end{aligned}
```

The Dirac equation is

```math
\boxed{\hat D(m)\boldsymbol\Psi=0.}
```

For commuting free operators,

```math
\hat P\,\mathrm{adj}\,\hat P
=\det\hat P
=\hat E^2-\hat{\mathbf p}^2
=-\Box.
```

Thus the spacetime mass shell $`\det P=m^2`$ becomes

```math
(\det\hat P-m^2)\boldsymbol\Psi=0
\quad\Longleftrightarrow\quad
(\Box+m^2)\boldsymbol\Psi=0.
```

Using $`\mathbf W(\hat P)^2=(\det\hat P)\mathbf I`$ gives the [Klein–Gordon factorization](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation):

```math
\begin{aligned}
\hat D(m)\hat D(-m)&=\hat D(-m)\hat D(m)\\
&=\mathbf W(\hat P)^2-m^2\mathbf I
=-(\Box+m^2)\mathbf I.
\end{aligned}
```

Hence every free Dirac solution satisfies Klein–Gordon.

## [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) with an electromagnetic potential

For constant mass $`m>0`$ and [charge](https://en.wikipedia.org/wiki/Electric_charge) $`q`$, [minimal coupling](https://en.wikipedia.org/wiki/Minimal_coupling) to the real [potential](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) $`\Phi=V+\mathbf A\cdot\boldsymbol{\sigma}`$ replaces the free operators by

```math
\begin{aligned}
\hat E&:=\mathrm{i}\partial_t-qV,\\
\hat{\mathbf p}&:=-\mathrm{i}\partial_{\mathbf r}-q\mathbf A,\\
\hat P&:=\hat E+\hat{\mathbf p}\cdot\boldsymbol{\sigma}.
\end{aligned}
```

```math
\boxed{\hat P=\mathrm{i}\partial^{*}-q\Phi.}
```

```math
\begin{aligned}
\hat D(m)&:=\mathbf W(\hat P)-m\mathbf I\\
&=\begin{pmatrix}
-m & \mathrm{adj}\,\hat P\\
\hat P & -m
\end{pmatrix}.
\end{aligned}
```

The coupled operators obey these field identities:

```math
\begin{aligned}
{}[\hat{\mathbf p}\cdot\boldsymbol{\sigma},\hat E]
&=-\mathrm{i}q\mathbf E\cdot\boldsymbol{\sigma},\\
(\hat{\mathbf p}\cdot\boldsymbol{\sigma})^2
&=\hat{\mathbf p}^2-q\mathbf B\cdot\boldsymbol{\sigma},\\
\hat{\mathbf p}^2+\mathrm{i}qF
&=(\hat{\mathbf p}\cdot\boldsymbol{\sigma})^2
-[\hat{\mathbf p}\cdot\boldsymbol{\sigma},\hat E].
\end{aligned}
```

When derivatives act on the potentials as well as the wavefunction, the opposite-mass product becomes

```math
\hat D(-m)\hat D(m)
=(\hat E^2-\hat{\mathbf p}^2-m^2)\mathbf I
-\mathrm{i}q\begin{pmatrix}F^{*}&0\\0&F\end{pmatrix}.
```

Every Dirac solution satisfies $`\hat D(-m)\hat D(m)\boldsymbol\Psi=0`$. Setting $`q=0`$ recovers the free Klein–Gordon equation.

## Low-energy [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) recovers [Pauli](https://en.wikipedia.org/wiki/Pauli_equation)

Remove the [rest-energy phase](https://en.wikipedia.org/wiki/Pauli_equation#Derivation) to define the envelopes $`\Psi_k'`$. On these envelopes, $`\hat E`$ measures energy relative to the rest mass:

```math
\begin{aligned}
\Psi_k'&:=e^{\mathrm{i}mt}\Psi_k,\\
\hat E\Psi_k&=e^{-\mathrm{i}mt}(m+\hat E)\Psi_k'.
\end{aligned}
```

The two rows of the original Dirac equation give the coupled envelope equations:

```math
\begin{aligned}
(m+\hat E-\hat{\mathbf p}\cdot\boldsymbol{\sigma})\Psi_2'&=m\Psi_1',\\
(m+\hat E+\hat{\mathbf p}\cdot\boldsymbol{\sigma})\Psi_1'&=m\Psi_2'.
\end{aligned}
```

Define the large and small components for the positive-energy branch:

```math
\begin{aligned}
\phi_1&:=\frac{\Psi_1'+\Psi_2'}{\sqrt2},\\
\phi_2&:=\frac{\Psi_2'-\Psi_1'}{\sqrt2}.
\end{aligned}
```

Adding and subtracting the envelope equations gives

```math
\begin{aligned}
\hat E\phi_1&=(\hat{\mathbf p}\cdot\boldsymbol{\sigma})\phi_2,\\
(2m+\hat E)\phi_2&=(\hat{\mathbf p}\cdot\boldsymbol{\sigma})\phi_1.
\end{aligned}
```

In the nonrelativistic limit with weak, slowly varying fields, neglect $`\hat E\phi_2`$ relative to $`2m\phi_2`$. Substitute the resulting small component into the first equation and use the field identity above:

```math
\begin{aligned}
\phi_2&\simeq\frac{\hat{\mathbf p}\cdot\boldsymbol{\sigma}}{2m}\phi_1,\\
\hat E\phi_1&\simeq\frac{(\hat{\mathbf p}\cdot\boldsymbol{\sigma})^2}{2m}\phi_1
=\frac{\hat{\mathbf p}^2-q\mathbf B\cdot\boldsymbol{\sigma}}{2m}\phi_1.
\end{aligned}
```

Restoring $`\hat E=\mathrm{i}\partial_t-qV`$ gives the [Pauli equation](https://en.wikipedia.org/wiki/Pauli_equation), retaining the kinetic and spin terms through order $`1/m`$:

```math
\boxed{\mathrm{i}\partial_t\phi_1=
\left(\frac{\hat{\mathbf p}^2}{2m}
+qV-\frac{q}{2m}\mathbf B\cdot\boldsymbol{\sigma}\right)\phi_1.}
```

## Low-energy [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) recovers [Schrödinger](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation)

For a charged scalar field $`\psi'`$, remove the same rest-energy phase by defining $`\psi:=e^{\mathrm{i}mt}\psi'`$. Substitution into the Klein–Gordon equation gives

```math
\begin{aligned}
0&=(\hat E^2-\hat{\mathbf p}^2-m^2)\psi'\\
&=e^{-\mathrm{i}mt}\bigl((m+\hat E)^2-\hat{\mathbf p}^2-m^2\bigr)\psi\\
&=e^{-\mathrm{i}mt}\bigl(2m\hat E+\hat E^2-\hat{\mathbf p}^2\bigr)\psi.
\end{aligned}
```

Cancel the common phase and solve for the residual-energy term:

```math
\hat E\psi=\frac{\hat{\mathbf p}^2-\hat E^2}{2m}\psi.
```

In the same nonrelativistic limit, neglect $`\hat E^2\psi`$ relative to $`2m\hat E\psi`$. Restoring $`\hat E=\mathrm{i}\partial_t-qV`$ gives the Schrödinger equation, retaining the kinetic term through order $`1/m`$:

```math
\boxed{\mathrm{i}\partial_t\psi=
\left(\frac{(-\mathrm{i}\partial_{\mathbf r}-q\mathbf A)^2}{2m}+qV\right)\psi.}
```

With $`\mathbf A=0`$ (hence $`\mathbf B=0`$), this recovers the conventional Schrödinger equation:

```math
\mathrm{i}\partial_t\psi=
\left(-\frac{\partial_{\mathbf r}^2}{2m}+qV\right)\psi.
```
