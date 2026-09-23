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
\mathrm{i}&=\sigma_1\sigma_2\sigma_3,
&\mathrm{i}^2&=-1,\\
\boldsymbol{\sigma}&=\{\sigma_1,\sigma_2,\sigma_3\},
&\mathrm{i}\boldsymbol{\sigma}&=\{\sigma_2\sigma_3,\sigma_3\sigma_1,\sigma_1\sigma_2\},\\
Z&=(a+b\mathrm{i})+(\mathbf A+\mathrm{i}\mathbf B)\cdot\boldsymbol{\sigma},
&a,b&\in\mathbb R,\quad\mathbf A,\mathbf B\in\mathbb R^3.
\end{aligned}}
```

Eight real components, written as a complex scalar and a complex vector:

```math
Z=S+\mathbf V\cdot\boldsymbol{\sigma}=X+\mathrm{i}Y,
\qquad X=a+\mathbf A\cdot\boldsymbol{\sigma},
\quad Y=b+\mathbf B\cdot\boldsymbol{\sigma}.
```

The component extractions are

```math
\begin{aligned}
\mathrm{re}(Z)&=X,\\
\mathrm{im}(Z)&=Y,\\
\mathrm{sc}(Z)&=S=a+\mathrm{i}b,\\
\mathrm{vec}(Z)&=\mathbf V=\mathbf A+\mathrm{i}\mathbf B.
\end{aligned}
```

## Product

The product carries both the [dot](https://en.wikipedia.org/wiki/Dot_product) and [cross](https://en.wikipedia.org/wiki/Cross_product) products:

```math
(\mathbf A\cdot\boldsymbol{\sigma})(\mathbf B\cdot\boldsymbol{\sigma})
=\mathbf A\cdot\mathbf B+\mathrm{i}(\mathbf A\times\mathbf B)\cdot\boldsymbol{\sigma}.
```

Here we use the convention $`\mathbf V^2:=\mathbf V\cdot\mathbf V`$.

The [commutator](https://en.wikipedia.org/wiki/Commutator) is

```math
[A,B]:=AB-BA.
```

## Core operations

[Determinant](https://en.wikipedia.org/wiki/Determinant), [trace](https://en.wikipedia.org/wiki/Trace_%28linear_algebra%29), [squared norm](https://en.wikipedia.org/wiki/Norm_%28mathematics%29), adjugate, [inverse](https://en.wikipedia.org/wiki/Invertible_matrix), and the two conjugations below.

Write $`Z=S+\mathbf V\cdot\boldsymbol{\sigma}`$, where $`S\in\mathbb C`$ is a complex scalar and $`\mathbf V\in\mathbb C^3`$ is a complex vector.

$`(\,)^{*}:=`$ [complex conjugation](https://en.wikipedia.org/wiki/Complex_conjugate) of complex numbers and a sign flip of $`\boldsymbol{\sigma}`$.

$`(\,)^{\mathsf H}:=`$ [Hermitian conjugation](https://en.wikipedia.org/wiki/Conjugate_transpose).

```math
\begin{gathered}
\begin{aligned}
\det Z&=S^2-\mathbf V\cdot\mathbf V,
&\mathrm{adj}\,Z&=S-\mathbf V\cdot\boldsymbol{\sigma},\\
\mathrm{tr}\,Z&=2S,
&\lVert Z\rVert^2&=S^{*}S+\mathbf V^{*}\cdot\mathbf V,\\
Z^{*}&=S^{*}-\mathbf V^{*}\cdot\boldsymbol{\sigma},
&Z^{\mathsf H}&=S^{*}+\mathbf V^{*}\cdot\boldsymbol{\sigma}.
\end{aligned}\\[6pt]
Z^{-1}=\frac{\mathrm{adj}\,Z}{\det Z}\qquad(\det Z\ne0).
\end{gathered}
```

## Rotations and boosts

For real spacetime $`X`$ and real [unit axis](https://en.wikipedia.org/wiki/Unit_vector) $`\mathbf u`$:

```math
X=t+\mathbf r\cdot\boldsymbol{\sigma},
\qquad X'=TXT^{\mathsf H}=t'+\mathbf r'\cdot\boldsymbol{\sigma}.
```

For $`T=R`$: [Rodrigues rotation](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula) through $`\theta`$, with [unit quaternion](https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation) $`Q=R`$, $`q_k=-\mathrm{i}\sigma_k`$.

```math
\begin{aligned}
R&=e^{-\mathrm{i}\theta\mathbf u\cdot\boldsymbol{\sigma}/2}
=\cos\frac\theta2-\mathrm{i}\mathbf u\cdot\boldsymbol{\sigma}\sin\frac\theta2,\\
t'&=t,\\
\mathbf r'&=\mathbf r\cos\theta+(\mathbf u\times\mathbf r)\sin\theta
+(\mathbf u\cdot\mathbf r)(1-\cos\theta)\mathbf u.
\end{aligned}
```

For $`T=L`$: [Lorentz boost](https://en.wikipedia.org/wiki/Lorentz_transformation) to a frame moving at $`+\beta\mathbf u`$, with [rapidity](https://en.wikipedia.org/wiki/Rapidity) $`\theta`$, $`\beta=\tanh\theta`$, and [Lorentz factor](https://en.wikipedia.org/wiki/Lorentz_factor) $`\gamma=\cosh\theta`$.

```math
\begin{aligned}
L&=e^{-\theta\mathbf u\cdot\boldsymbol{\sigma}/2}
=\cosh\frac\theta2-\mathbf u\cdot\boldsymbol{\sigma}\sinh\frac\theta2
=L^{\mathsf H},\\
t'&=\gamma(t-\beta\mathbf u\cdot\mathbf r),\\
\mathbf r'&=\mathbf r+(\gamma-1)(\mathbf u\cdot\mathbf r)\mathbf u
-\gamma\beta t\mathbf u.
\end{aligned}
```

## Projections and spectral decomposition

```math
Z=S+\theta\mathbf u\cdot\boldsymbol{\sigma},
\qquad S,\theta\in\mathbb C,\quad \mathbf u\in\mathbb C^3,\quad \mathbf u^2=1.
```

The eigenvalues and complementary [projectors](https://en.wikipedia.org/wiki/Paravector#Null_paravectors_as_projectors) are

```math
\begin{gathered}
\lambda_\pm=S\pm\theta,
\qquad \Pi_\pm=\frac12(1\pm\mathbf u\cdot\boldsymbol{\sigma}),\\[4pt]
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
&\mathrm ds^2&=\det(\mathrm dX)=\mathrm dt^2-\mathrm d\mathbf r^2,\\
U&=\frac{\mathrm dX}{\mathrm ds}=\gamma(1+\mathbf v\cdot\boldsymbol{\sigma}),
&\gamma&=(1-\mathbf v^2)^{-1/2},\\
P&=mU=\underbrace{E}_{\text{energy}}+
\underbrace{\mathbf p\cdot\boldsymbol{\sigma}}_{\text{momentum}},
&\det P&=E^2-\mathbf p^2=m^2.
\end{aligned}
```

## [Maxwell](https://en.wikipedia.org/wiki/Maxwell%27s_equations) in one equation

Combine the [electric](https://en.wikipedia.org/wiki/Electric_field) and [magnetic](https://en.wikipedia.org/wiki/Magnetic_field) fields and use the [paravector derivative](https://en.wikipedia.org/wiki/Paravector#Paragradient). Hats mark named [energy](https://en.wikipedia.org/wiki/Energy_operator), [momentum](https://en.wikipedia.org/wiki/Momentum_operator), and [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) operators; the derivative symbols $`\partial`$ and $`\Box`$ remain unhatted:

```math
F=(\mathbf E+\mathrm{i}\mathbf B)\cdot\boldsymbol{\sigma},
\qquad \partial=\partial_t+\partial_{\mathbf r}\cdot\boldsymbol{\sigma},
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
\partial^{*}\partial&=\Box=\partial_t^2-\partial_{\mathbf r}^2,\\
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
\Phi=V+\mathbf A\cdot\boldsymbol{\sigma},
\qquad S=\mathrm{sc}(\partial\Phi).
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
\Phi'&=\Phi+\partial^{*}\lambda,\\
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

Package two paravector [amplitudes](https://en.wikipedia.org/wiki/Probability_amplitude) into $`\boldsymbol\Psi=(\Psi_1,\Psi_2)^{\mathsf T}`$, each with two independent complex components. With $`\mathbf I`$ the block identity, define $`\mathbf W`$ in the [Weyl (chiral) representation](https://en.wikipedia.org/wiki/Gamma_matrices#Weyl_%28chiral%29_basis):

```math
\mathbf W(Z)=
\begin{pmatrix}
0 & Z \\
\mathrm{adj}\,Z & 0
\end{pmatrix},
\qquad \mathbf W(Z)^2=\det(Z)\mathbf I.
```

Define the [free Dirac operator](https://en.wikipedia.org/wiki/Dirac_equation) and its equation:

```math
\hat D(m):=\mathbf W(\mathrm{i}\partial)-m\mathbf I,
\qquad \boxed{\hat D(m)\boldsymbol\Psi=0.}
```

Since $`\mathbf W(\mathrm{i}\partial)^2=-\partial^{*}\partial\,\mathbf I=-\Box\mathbf I`$, the opposite mass signs factor the [Klein–Gordon operator](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation):

```math
\begin{aligned}
\hat D(m)\hat D(-m)&=\hat D(-m)\hat D(m)\\
&=\mathbf W(\mathrm{i}\partial)^2-m^2\mathbf I
=-(\Box+m^2)\mathbf I.
\end{aligned}
```

Applying this product to a free Dirac solution gives

```math
\hat D(-m)\hat D(m)\boldsymbol\Psi=-(\Box+m^2)\boldsymbol\Psi=0
\quad\Longrightarrow\quad\boxed{(\Box+m^2)\boldsymbol\Psi=0.}
```

## [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) with an electromagnetic potential

For real [potential](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) $`\Phi=V+\mathbf A\cdot\boldsymbol{\sigma}`$, constant mass $`m>0`$, and [charge](https://en.wikipedia.org/wiki/Electric_charge) $`q`$, [minimal coupling](https://en.wikipedia.org/wiki/Minimal_coupling) gives

```math
\begin{aligned}
\hat D(m)&:=\mathbf W(\mathrm{i}\partial-q\Phi^{*})-m\mathbf I,\\
\hat D(-m)&=\mathbf W(\mathrm{i}\partial-q\Phi^{*})+m\mathbf I.
\end{aligned}
\qquad\boxed{\hat D(m)\boldsymbol\Psi=0.}
```

Define

```math
\hat E=\mathrm{i}\partial_t-qV,\qquad
\hat{\mathbf p}=-\mathrm{i}\partial_{\mathbf r}-q\mathbf A.
```

These obey

```math
\begin{aligned}
[\hat{\mathbf p}\cdot\boldsymbol{\sigma},\hat E]
&=-\mathrm{i}q\mathbf E\cdot\boldsymbol{\sigma},\\
(\hat{\mathbf p}\cdot\boldsymbol{\sigma})^2
&=\hat{\mathbf p}^2-q\mathbf B\cdot\boldsymbol{\sigma},\\
\hat{\mathbf p}^2+\mathrm{i}qF
&=(\hat{\mathbf p}\cdot\boldsymbol{\sigma})^2
-[\hat{\mathbf p}\cdot\boldsymbol{\sigma},\hat E].
\end{aligned}
```

The product, including potential derivatives, is

```math
\hat D(-m)\hat D(m)
=(\hat E^2-\hat{\mathbf p}^2-m^2)\mathbf I
-\mathrm{i}q\begin{pmatrix}F^{*}&0\\0&F\end{pmatrix}.
```

Dirac solutions satisfy $`\hat D(-m)\hat D(m)\boldsymbol\Psi=0`$, recovering Klein–Gordon at $`q=0`$.

## Low-energy [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) recovers [Pauli](https://en.wikipedia.org/wiki/Pauli_equation)

The lower block of $`\hat D(-m)\hat D(m)\boldsymbol\Psi=0`$ is

```math
(\hat E^2-\hat{\mathbf p}^2-m^2-\mathrm{i}qF)\Psi_2=0.
```

Remove the [rest phase](https://en.wikipedia.org/wiki/Pauli_equation#Derivation), $`\Psi_k=e^{-\mathrm{i}mt}\Psi_k'`$. Differentiating it gives

```math
\begin{aligned}
\hat E\Psi_k&=(\mathrm{i}\partial_t-qV)(e^{-\mathrm{i}mt}\Psi_k')
=e^{-\mathrm{i}mt}(m+\hat E)\Psi_k',\\
\hat E^2\Psi_k&=e^{-\mathrm{i}mt}(m+\hat E)^2\Psi_k'.
\end{aligned}
```

The spatial operator and $`F`$ commute with this time-only scalar phase. Cancel it and expand:

```math
\begin{aligned}
0&=\bigl((m+\hat E)^2-\hat{\mathbf p}^2-m^2-\mathrm{i}qF\bigr)\Psi_2'\\
&=\bigl(2m\hat E+\hat E^2-\hat{\mathbf p}^2-\mathrm{i}qF\bigr)\Psi_2'.
\end{aligned}
```

Using the identities above, rearrange for the exact residual-energy recursion:

```math
\begin{aligned}
\hat E\Psi_2'
&=\frac{\hat{\mathbf p}^2+\mathrm{i}qF-\hat E^2}{2m}\Psi_2'\\
&=\frac{(\hat{\mathbf p}\cdot\boldsymbol{\sigma})^2-\hat E^2
-[\hat{\mathbf p}\cdot\boldsymbol{\sigma},\hat E]}{2m}\Psi_2'.
\end{aligned}
```

To retain the first-order coupling, use the envelope equations

```math
\begin{aligned}
(m+\hat E-\hat{\mathbf p}\cdot\boldsymbol{\sigma})\Psi_2'&=m\Psi_1',\\
(m+\hat E+\hat{\mathbf p}\cdot\boldsymbol{\sigma})\Psi_1'&=m\Psi_2'.
\end{aligned}
```

Set $`\phi_1=(\Psi_1'+\Psi_2')/\sqrt2`$ and $`\phi_2=(\Psi_2'-\Psi_1')/\sqrt2`$. Adding and subtracting gives

```math
\hat E\phi_1=(\hat{\mathbf p}\cdot\boldsymbol{\sigma})\phi_2,
\qquad
(2m+\hat E)\phi_2=(\hat{\mathbf p}\cdot\boldsymbol{\sigma})\phi_1.
```

For slow positive-energy states in weak, slowly varying fields, $`\phi_2`$ is small. Neglect its residual energy compared with $`2m`$:

```math
\phi_2=\frac{\hat{\mathbf p}\cdot\boldsymbol{\sigma}}{2m}\phi_1
-\frac{\hat E\phi_2}{2m}
\simeq\frac{\hat{\mathbf p}\cdot\boldsymbol{\sigma}}{2m}\phi_1.
```

Substitute into the $`\phi_1`$ equation and use the sigma product:

```math
\hat E\phi_1\simeq\frac{(\hat{\mathbf p}\cdot\boldsymbol{\sigma})^2}{2m}\phi_1
=\frac{\hat{\mathbf p}^2-q\mathbf B\cdot\boldsymbol{\sigma}}{2m}\phi_1.
```

Restoring $`\hat E=\mathrm{i}\partial_t-qV`$ gives the leading-order [Pauli equation](https://en.wikipedia.org/wiki/Pauli_equation):

```math
\boxed{\mathrm{i}\partial_t\phi_1=
\left[\frac{\hat{\mathbf p}^2}{2m}
+qV-\frac{q}{2m}\mathbf B\cdot\boldsymbol{\sigma}\right]\phi_1.}
```
