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
\mathrm{adj}\,Z&:=S-\mathbf V\cdot\boldsymbol{\sigma} = Z^{*\mathsf H},&\det Z&:=S^2-\mathbf V^2 = Z\,\mathrm{adj}\,Z,\\
\mathrm{tr}\,Z&:=2S=Z+\mathrm{adj}\,Z,
&\lVert Z\rVert^2&:=SS^{*}+\mathbf V\cdot\mathbf V^{*} = \mathrm{sc}(Z Z^\mathsf H).
\end{aligned}\\[6pt]
Z^{-1}=\frac{\mathrm{adj}\,Z}{\det Z}\qquad \mathrm{if} \quad \det Z\ne0.
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

In the [Weyl representation](https://en.wikipedia.org/wiki/Gamma_matrices#Weyl_%28chiral%29_basis), write $`\boldsymbol\Psi:=(\Psi_+,\Psi_-)^{\mathsf T}`$, with two complex components in each entry. Let $`\mathbf I`$ be the block identity.

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

Equivalently,

```math
(\hat E\pm\hat{\mathbf p}\cdot\boldsymbol{\sigma})\Psi_\pm=m\Psi_\mp.
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

For constant mass $`m>0`$ and [charge](https://en.wikipedia.org/wiki/Electric_charge) $`q`$, [minimal coupling](https://en.wikipedia.org/wiki/Minimal_coupling) to the real [potential](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) $`\Phi=V+\mathbf A\cdot\boldsymbol{\sigma}`$ defines the gauge-covariant energy and momentum operators

```math
\begin{aligned}
\hat E_q&:=\hat E-qV=\mathrm{i}\partial_t-qV,\\
\hat{\mathbf p}_q&:=\hat{\mathbf p}-q\mathbf A=-\mathrm{i}\partial_{\mathbf r}-q\mathbf A,\\
\hat P_q&:=\hat E_q+\hat{\mathbf p}_q\cdot\boldsymbol{\sigma}.
\end{aligned}
```

```math
\boxed{\hat P_q=\hat P-q\Phi=\mathrm{i}\partial^{*}-q\Phi.}
```

At $`q=0`$, the coupled operators reduce to the canonical ones: $`\hat E_0=\hat E`$, $`\hat{\mathbf p}_0=\hat{\mathbf p}`$, and $`\hat P_0=\hat P`$. The coupled Dirac operator is

```math
\begin{aligned}
\hat D_q(m)&:=\mathbf W(\hat P_q)-m\mathbf I\\
&=\begin{pmatrix}
-m & \mathrm{adj}\,\hat P_q\\
\hat P_q & -m
\end{pmatrix}.
\end{aligned}
```

The coupled Dirac equation is

```math
\boxed{\hat D_q(m)\boldsymbol\Psi=0.}
```

Equivalently,

```math
(\hat E_q\pm\hat{\mathbf p}_q\cdot\boldsymbol{\sigma})\Psi_\pm=m\Psi_\mp.
```

The coupled operators obey these field identities:

```math
\begin{aligned}
{}[\hat{\mathbf p}_q\cdot\boldsymbol{\sigma},\hat E_q]
&=-\mathrm{i}q\mathbf E\cdot\boldsymbol{\sigma},\\
(\hat{\mathbf p}_q\cdot\boldsymbol{\sigma})^2
&=\hat{\mathbf p}_q^2-q\mathbf B\cdot\boldsymbol{\sigma},\\
\hat{\mathbf p}_q^2+\mathrm{i}qF
&=(\hat{\mathbf p}_q\cdot\boldsymbol{\sigma})^2
-[\hat{\mathbf p}_q\cdot\boldsymbol{\sigma},\hat E_q].
\end{aligned}
```

When derivatives act on the potentials as well as the wavefunction, the opposite-mass product becomes

```math
\hat D_q(-m)\hat D_q(m)
=(\hat E_q^2-\hat{\mathbf p}_q^2-m^2)\mathbf I
-\mathrm{i}q\begin{pmatrix}F^{*}&0\\0&F\end{pmatrix}.
```

Every Dirac solution satisfies $`\hat D_q(-m)\hat D_q(m)\boldsymbol\Psi=0`$. Setting $`q=0`$ recovers the free Klein–Gordon equation.

## Low-energy [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) recovers [Schrödinger](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation)

For constant $`m`$, remove the rest-energy phase, $`\psi'=e^{-\mathrm{i}mt}\psi`$. The product rule gives

```math
\begin{aligned}
\hat E_q\psi'
&=(\mathrm{i}\partial_t-qV)(e^{-\mathrm{i}mt}\psi)\\
&=e^{-\mathrm{i}mt}(m\psi+\mathrm{i}\partial_t\psi-qV\psi)\\
&=e^{-\mathrm{i}mt}(m+\hat E_q)\psi,\\[6pt]
\hat{\mathbf p}_q\psi'&=e^{-\mathrm{i}mt}\hat{\mathbf p}_q\psi.
\end{aligned}
```

Substitute into Klein–Gordon:

```math
\begin{aligned}
0&=(\hat E_q^2-\hat{\mathbf p}_q^2-m^2)\psi'\\
&=e^{-\mathrm{i}mt}\bigl((m+\hat E_q)^2-\hat{\mathbf p}_q^2-m^2\bigr)\psi\\
&=e^{-\mathrm{i}mt}\bigl(2m\hat E_q+\hat E_q^2-\hat{\mathbf p}_q^2\bigr)\psi.
\end{aligned}
```

For any $`q`$ and constant $`m>0`$, canceling the phase gives the exact recursion:

```math
\hat E_q\psi=\frac{\hat{\mathbf p}_q^2-\hat E_q^2}{2m}\psi.
```

If $`\hat E_q`$ and $`\hat{\mathbf p}_q^2`$ commute, iteration on the positive-energy branch gives

```math
\hat E_q\psi=\left(\frac{\hat{\mathbf p}_q^2}{2m}
-\frac{\hat{\mathbf p}_q^4}{8m^3}+\cdots\right)\psi.
```

For weak, slowly varying fields, neglect $`\hat E_q^2\psi`$ relative to $`2m\hat E_q\psi`$ to obtain Schrödinger’s equation:

```math
\boxed{\mathrm{i}\partial_t\psi=
\left(\frac{(-\mathrm{i}\partial_{\mathbf r}-q\mathbf A)^2}{2m}+qV\right)\psi.}
```

With $`\mathbf A=0`$ (hence $`\mathbf B=0`$), this recovers the conventional Schrödinger equation:

```math
\mathrm{i}\partial_t\psi=
\left(-\frac{\partial_{\mathbf r}^2}{2m}+qV\right)\psi.
```

## Low-energy [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) recovers [Pauli](https://en.wikipedia.org/wiki/Pauli_equation)

Using the same phase shift, $`\Psi_\pm':=e^{\mathrm{i}mt}\Psi_\pm`$, gives

```math
(m+\hat E_q\pm\hat{\mathbf p}_q\cdot\boldsymbol{\sigma})\Psi_\pm'=m\Psi_\mp'.
```

Define the large and small components for the positive-energy branch:

```math
\begin{aligned}
\phi_+&:=\frac{\Psi_+'+\Psi_-'}{\sqrt2},\\
\phi_-&:=\frac{\Psi_-'-\Psi_+'}{\sqrt2}.
\end{aligned}
```

Adding and subtracting the envelope equations gives

```math
\begin{aligned}
\hat E_q\phi_+&=(\hat{\mathbf p}_q\cdot\boldsymbol{\sigma})\phi_-,\\
(2m+\hat E_q)\phi_-&=(\hat{\mathbf p}_q\cdot\boldsymbol{\sigma})\phi_+.
\end{aligned}
```

In the same nonrelativistic limit, neglect $`\hat E_q\phi_-`$ relative to $`2m\phi_-`$. Substitute the resulting small component into the first equation and use the field identity above:

```math
\begin{aligned}
\phi_-&\simeq\frac{\hat{\mathbf p}_q\cdot\boldsymbol{\sigma}}{2m}\phi_+,\\
\hat E_q\phi_+&\simeq\frac{(\hat{\mathbf p}_q\cdot\boldsymbol{\sigma})^2}{2m}\phi_+
=\frac{\hat{\mathbf p}_q^2-q\mathbf B\cdot\boldsymbol{\sigma}}{2m}\phi_+.
\end{aligned}
```

Define $`\psi:=\phi_+`$. Restoring $`\hat E_q=\mathrm{i}\partial_t-qV`$ gives the [Pauli equation](https://en.wikipedia.org/wiki/Pauli_equation), retaining the kinetic and spin terms through order $`1/m`$:

```math
\boxed{\mathrm{i}\partial_t\psi=
\left(\frac{\hat{\mathbf p}_q^2}{2m}
+qV-\frac{q}{2m}\mathbf B\cdot\boldsymbol{\sigma}\right)\psi.}
```

## Spin up and spin down in a uniform magnetic field

For a nonzero, constant magnetic field $`\mathbf B`$, define its [unit direction](https://en.wikipedia.org/wiki/Unit_vector) and [spin projectors](https://en.wikipedia.org/wiki/Pauli_matrices#Eigenvectors_and_eigenvalues):

```math
\mathbf u:=\frac{\mathbf B}{\lVert\mathbf B\rVert},
\qquad \Pi_\pm:=\frac12(1\pm\mathbf u\cdot\boldsymbol{\sigma}),
\qquad \psi_\pm:=\Pi_\pm\psi,
\qquad \psi=\psi_++\psi_-.
```

The projector identity gives [spin](https://en.wikipedia.org/wiki/Spin-1/2) $`\pm\tfrac12`$ along $`\mathbf u`$:

```math
(\mathbf u\cdot\boldsymbol{\sigma})\Pi_\pm=\pm\Pi_\pm,
\qquad
\frac{\mathbf u\cdot\boldsymbol{\sigma}}2\psi_\pm=\pm\frac12\psi_\pm.
```

Since the projectors are constant, they commute with $`\hat E_q`$ and $`\hat{\mathbf p}_q^2`$. Projecting the Pauli equation gives

```math
\begin{aligned}
2m\hat E_q\psi_\pm
&=\Pi_\pm(\hat{\mathbf p}_q^2-q\mathbf B\cdot\boldsymbol{\sigma})\psi\\
&=(\hat{\mathbf p}_q^2\mp q\lVert\mathbf B\rVert)\psi_\pm.
\end{aligned}
```

Thus the two spin components evolve independently, with opposite [Zeeman shifts](https://en.wikipedia.org/wiki/Zeeman_effect):

```math
\boxed{\mathrm{i}\partial_t\psi_\pm=
\left(\frac{\hat{\mathbf p}_q^2\mp q\lVert\mathbf B\rVert}{2m}+qV\right)\psi_\pm.}
```

The [kinetic momentum](https://en.wikipedia.org/wiki/Momentum_operator#Electromagnetic_field) $`\hat{\mathbf p}_q=-\mathrm{i}\partial_{\mathbf r}-q\mathbf A`$ retains the orbital coupling.

## Spin measurement probabilities

Prepare spin up along a real unit axis $`\mathbf u`$, so $`\psi=\Pi_+\psi\ne0`$. Measure along another real unit axis $`\mathbf u'`$, with $`\cos\theta:=\mathbf u\cdot\mathbf u'`$. Its projectors satisfy

```math
\Pi_\pm':=\frac12(1\pm\mathbf u'\cdot\boldsymbol{\sigma}),
\qquad \Pi_+\Pi_\pm'\Pi_+=\frac{1\pm\cos\theta}{2}\Pi_+.
```

The [Born rule](https://en.wikipedia.org/wiki/Born_rule) therefore gives the probabilities of measuring [spin](https://en.wikipedia.org/wiki/Spin-1/2) $`\pm\tfrac12`$ along $`\mathbf u'`$:

```math
\begin{gathered}
\operatorname{prob}(\pm)=\frac{\lVert\Pi_\pm'\psi\rVert^2}{\lVert\psi\rVert^2}
=\frac{1\pm\cos\theta}{2},\\[6pt]
\boxed{\operatorname{prob}(+)=\cos^2\frac\theta2,\quad \operatorname{prob}(-)=\sin^2\frac\theta2.}
\end{gathered}
```

Aligned axes give spin up with certainty; perpendicular axes give equal probabilities. Initial spin down swaps the probabilities.
