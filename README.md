# Space–Time Algebra

**A common [paravector](https://en.wikipedia.org/wiki/Paravector) language is provided by [complex scalars](https://en.wikipedia.org/wiki/Complex_number) and [three-vectors](https://en.wikipedia.org/wiki/Euclidean_vector) for [quaternions](https://en.wikipedia.org/wiki/Quaternion), [rotations](https://en.wikipedia.org/wiki/Rotation_%28mathematics%29), [Lorentz boosts](https://en.wikipedia.org/wiki/Lorentz_transformation), [projections](https://en.wikipedia.org/wiki/Paravector#Null_paravectors_as_projectors), [spacetime geometry](https://en.wikipedia.org/wiki/Minkowski_space), [relativistic particle dynamics](https://en.wikipedia.org/wiki/Relativistic_mechanics), [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations), [electromagnetic waves](https://en.wikipedia.org/wiki/Electromagnetic_radiation) and [forces](https://en.wikipedia.org/wiki/Lorentz_force), [gauge coupling](https://en.wikipedia.org/wiki/Minimal_coupling), the [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) and [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) equations, their [Schrödinger](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation) and [Pauli](https://en.wikipedia.org/wiki/Pauli_equation) low-energy limits, and [spin-one-half](https://en.wikipedia.org/wiki/Spin-1/2) [amplitudes](https://en.wikipedia.org/wiki/Probability_amplitude) and [measurement probabilities](https://en.wikipedia.org/wiki/Born_rule).**

[Paper (PDF)](./sta_notes.pdf) · [LaTeX source](./sta_notes.tex)

[Natural units](https://en.wikipedia.org/wiki/Natural_units) $`\hbar=c=1`$, [rationalized electromagnetic units](https://en.wikipedia.org/wiki/Heaviside%E2%80%93Lorentz_units), and [signature](https://en.wikipedia.org/wiki/Metric_signature) $`(+,-,-,-)`$ are used.

## One complex [paravector](https://en.wikipedia.org/wiki/Paravector)

Space-time algebra is developed using complex paravectors. The sigmas are algebraic objects with order-dependent multiplication rules. The sigmas and the general complex paravector are defined by:

```math
\boxed{
\begin{aligned}
\sigma_n^2&=1,
&\sigma_n\sigma_m&=-\sigma_m\sigma_n\quad(n\ne m),\\
\mathrm{i}&:=\sigma_1\sigma_2\sigma_3,
&\mathrm{i}^2&=-1,\\
𝝈&:=\{\sigma_1,\sigma_2,\sigma_3\},
&\mathrm{i}𝝈&=\{\sigma_2\sigma_3,\sigma_3\sigma_1,\sigma_1\sigma_2\},\\
\mathcal{Z}&:=(a+b\mathrm{i})+(𝐗+\mathrm{i}𝐘)\cdot𝝈,
&a,b&\in\mathbb R,\quad𝐗,𝐘\in\mathbb R^3.
\end{aligned}}
```

These eight real components are grouped into a complex scalar and a complex vector:

```math
\begin{aligned}
\mathcal{X}:&=a+𝐗\cdot𝝈,\\
\mathcal{Y}:&=b+𝐘\cdot𝝈,\\[4pt]
\mathcal{Z}&=\mathcal{X}+\mathrm{i}\mathcal{Y}\\
&=(a+b\mathrm{i})+(𝐗+\mathrm{i}𝐘)\cdot𝝈\\
&=c+𝐙\cdot𝝈.
\end{aligned}
```

The component extractions are:

```math
\begin{aligned}
\mathrm{re}(\mathcal{Z})&:=\mathcal{X},\\
\mathrm{im}(\mathcal{Z})&:=\mathcal{Y},\\
\mathrm{sc}(\mathcal{Z})&:=c=a+\mathrm{i}b,\\
\mathrm{vec}(\mathcal{Z})&:=𝐙=𝐗+\mathrm{i}𝐘.
\end{aligned}
```

## Product

Both the [dot](https://en.wikipedia.org/wiki/Dot_product) and [cross](https://en.wikipedia.org/wiki/Cross_product) products are included in the product:

```math
(𝐗\cdot𝝈)(𝐘\cdot𝝈)
=𝐗\cdot𝐘+\mathrm{i}(𝐗\times𝐘)\cdot𝝈.
```

The convention $`𝐗^2:=𝐗\cdot𝐗`$ is used.

The [commutator](https://en.wikipedia.org/wiki/Commutator) is:

```math
[\mathcal{X},\mathcal{Y}]:=\mathcal{X}\mathcal{Y}-\mathcal{Y}\mathcal{X}.
```

For $`n,m,k\in\{1,2,3\}`$, $`\delta_{nm}`$ is defined to be $`1`$ if $`n=m`$ and $`0`$ otherwise. The value of $`\epsilon_{nmk}`$ is defined to be $`+1`$ for even permutations of $`(1,2,3)`$, $`-1`$ for odd permutations, and $`0`$ if any indices repeat. The ordered multiplication rule and commutator are then expressed as:

```math
\begin{aligned}
\sigma_n\sigma_m&=\delta_{nm}
+\sum_{k=1}^{3}\epsilon_{nmk}\,\mathrm{i}\sigma_k,\\[6pt]
[\sigma_n,\sigma_m]&=2\sum_{k=1}^{3}\epsilon_{nmk}\,\mathrm{i}\sigma_k.
\end{aligned}
```

## Hermitian and sigma conjugation

**Hermitian conjugation** $`(\,\cdot\,)^{\mathsf H}`$ and **sigma conjugation** $`(\,\cdot\,)^{*}`$ are defined by:

```math
\begin{aligned}
(\,\cdot\,)^{\mathsf H}
&:=\text{the multiplication order of the sigmas is reversed},\\[6pt]
(\,\cdot\,)^{*}
&:=\text{every }\sigma_n\text{ is replaced by }-\sigma_n
\text{ with product order preserved}.
\end{aligned}
```

Both operations are applied term by term, with real coefficients left unchanged. For two sigmas:

```math
\begin{aligned}
(\sigma_m\sigma_n)^{\mathsf H}&=\sigma_n\sigma_m,\\[4pt]
(\sigma_m\sigma_n)^{*}&=(-\sigma_m)(-\sigma_n).
\end{aligned}
```

Since $`\mathrm{i}=\sigma_1\sigma_2\sigma_3`$, the sign of $`\mathrm{i}`$ is reversed by both operations:

```math
\begin{aligned}
\mathrm{i}^{\mathsf H}&=\sigma_3\sigma_2\sigma_1=-\mathrm{i},\\[4pt]
\mathrm{i}^{*}&=(-\sigma_1)(-\sigma_2)(-\sigma_3)=-\mathrm{i}.
\end{aligned}
```

For $`\mathcal{Z}=c+𝐙\cdot𝝈`$, with $`c\in\mathbb C`$ and $`𝐙\in\mathbb C^3`$, the complex coefficients are therefore conjugated by both operations, while the vector sign is reversed only by sigma conjugation:

```math
\begin{aligned}
\mathcal{Z}^{\mathsf H}&=c^{*}+𝐙^{*}\cdot𝝈,\\[4pt]
\mathcal{Z}^{*}&=c^{*}-𝐙^{*}\cdot𝝈.
\end{aligned}
```

## Core operations

The paravector [adjugate](https://en.wikipedia.org/wiki/Adjugate_matrix), [determinant](https://en.wikipedia.org/wiki/Determinant), and [inverse](https://en.wikipedia.org/wiki/Invertible_matrix) are given by:

```math
\begin{aligned}
\mathrm{adj}\,\mathcal{Z}&:=c-𝐙\cdot𝝈=\mathcal{Z}^{*\mathsf H},\\[4pt]
\det \mathcal{Z}&:=c^2-𝐙^2=\mathcal{Z}\,\mathrm{adj}\,\mathcal{Z},\\[6pt]
\mathcal{Z}^{-1}&=\frac{\mathrm{adj}\,\mathcal{Z}}{\det \mathcal{Z}}\qquad\text{if }\det \mathcal{Z}\ne0.
\end{aligned}
```

The paravector [trace](https://en.wikipedia.org/wiki/Trace_%28linear_algebra%29) and [squared norm](https://en.wikipedia.org/wiki/Norm_%28mathematics%29) are given by:

```math
\begin{aligned}
\mathrm{tr}\,\mathcal{Z}&:=2c=\mathcal{Z}+\mathrm{adj}\,\mathcal{Z},\\[4pt]
\lVert \mathcal{Z}\rVert^2&:=cc^{*}+𝐙\cdot𝐙^{*}=\mathrm{sc}(\mathcal{Z}\mathcal{Z}^{\mathsf H}).
\end{aligned}
```

## [Quaternions](https://en.wikipedia.org/wiki/Quaternion)

The three imaginary quaternion units are defined by:

```math
\begin{gathered}
\nu_k:=-\mathrm{i}\sigma_k,\qquad 𝝂:=-\mathrm{i}𝝈=\{\nu_1,\nu_2,\nu_3\},\\
\nu_1^2=\nu_2^2=\nu_3^2=\nu_1\nu_2\nu_3=-1.
\end{gathered}
```

For $`a\in\mathbb R`$ and $`𝐐\in\mathbb R^3`$, the quaternion and its conjugate are:

```math
\begin{aligned}
\mathcal{Q}&:=a+𝐐\cdot𝝂=a-\mathrm{i}𝐐\cdot𝝈,\\
\mathcal{Q}^{\mathsf H}&:=a-𝐐\cdot𝝂=a+\mathrm{i}𝐐\cdot𝝈,
\end{aligned}
```

with

```math
\mathcal{Q}\mathcal{Q}^{\mathsf H}=\det \mathcal{Q}=a^2+𝐐^2.
```

Nonzero quaternions are normalized with:

```math
\mathcal{Q}\leftarrow\frac{\mathcal{Q}}{\sqrt{\mathcal{Q}\mathcal{Q}^{\mathsf H}}},
```

Thus, with $`\lVert𝐐\rVert^2=𝐐^2`$ for real $`𝐐`$, the normalization is expressed as:

```math
\mathcal{Q}\mathcal{Q}^{\mathsf H}=a^2+\lVert𝐐\rVert^2=1.
```

For a unit quaternion with $`𝐐\ne0`$, the rotation axis and angle are specified by:

```math
𝐮=\frac{𝐐}{\lVert𝐐\rVert},\qquad
a=\cos\frac{\theta}{2},\qquad
𝐐=𝐮\sin\frac{\theta}{2}.
```

A rotation of $`𝐫\in\mathbb R^3`$ about $`𝐮`$ by angle $`\theta`$ is expressed by:

```math
\mathcal{X}'=𝐫'\cdot𝝂, \qquad
\mathcal{X}=𝐫\cdot𝝂, \qquad
\mathcal{X}'=\mathcal{Q}\mathcal{X}\mathcal{Q}^{\mathsf H}.
```

## Rotations and boosts

For real spacetime $`\mathcal{X}`$ and a real [unit axis](https://en.wikipedia.org/wiki/Unit_vector) $`𝐮`$, the transformation is written as:

```math
\mathcal{X}:=t+𝐫\cdot𝝈,
\qquad \mathcal{X}':=\mathcal{T}\mathcal{X}\mathcal{T}^{\mathsf H}=t'+𝐫'\cdot𝝈.
```

For $`\mathcal{T}=\mathcal{R}`$, a [Rodrigues rotation](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula) through $`\theta`$, with [unit quaternion](https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation) $`\mathcal{Q}=\mathcal{R}`$, is given by:

```math
\begin{aligned}
\mathcal{R}&:=e^{-\mathrm{i}\theta𝐮\cdot𝝈/2}
=\cos\frac\theta2-\mathrm{i}𝐮\cdot𝝈\sin\frac\theta2,\\
\mathcal{X}'&=\mathcal{R}\mathcal{X}\mathcal{R}^{\mathsf H},\\
t'&=t,\\
𝐫'&=𝐫\cos\theta+(𝐮\times𝐫)\sin\theta
+(𝐮\cdot𝐫)(1-\cos\theta)𝐮.
\end{aligned}
```

For $`\mathcal{T}=\mathcal{L}`$, a [Lorentz boost](https://en.wikipedia.org/wiki/Lorentz_transformation) to a frame moving at $`+\beta𝐮`$, with [rapidity](https://en.wikipedia.org/wiki/Rapidity) $`\theta`$, $`\beta:=\tanh\theta`$, and [Lorentz factor](https://en.wikipedia.org/wiki/Lorentz_factor) $`\gamma:=\cosh\theta`$ is considered.

The parallel and perpendicular components are defined by:

```math
𝐫_\parallel:=(𝐮\cdot𝐫)𝐮,
\qquad 𝐫_\perp:=𝐫-𝐫_\parallel.
```

The boost and transformed coordinates are:

```math
\begin{aligned}
\mathcal{L}&:=e^{-\theta𝐮\cdot𝝈/2}
=\cosh\frac\theta2-𝐮\cdot𝝈\sinh\frac\theta2
=\mathcal{L}^{\mathsf H},\\
\mathcal{X}'&=\mathcal{L}\mathcal{X}\mathcal{L}^{\mathsf H},\\
t'&=\gamma(t-\beta𝐮\cdot𝐫),\\
𝐫'&=𝐫_\perp+\gamma(𝐫_\parallel-\beta t𝐮).
\end{aligned}
```

## Projections and spectral decomposition

A paravector of the following form is considered:

```math
\mathcal{Z}=c+z𝐮\cdot𝝈,
\qquad c,z\in\mathbb C,\quad 𝐮\in\mathbb C^3,\quad 𝐮^2=1.
```

The eigenvalues and complementary [projectors](https://en.wikipedia.org/wiki/Paravector#Null_paravectors_as_projectors) are:

```math
\begin{gathered}
\lambda_\pm:=c\pm z,
\qquad \Pi_\pm:=\frac12(1\pm𝐮\cdot𝝈),\\[4pt]
\Pi_\pm^2=\Pi_\pm,\qquad \Pi_+\Pi_-=0,\qquad \Pi_-+\Pi_+=1.
\end{gathered}
```

For real $`𝐮`$, these are also Hermitian projectors. For $`f`$ analytic near $`\lambda_\pm`$, the following spectral decomposition is obtained:

```math
\begin{aligned}
\mathcal{Z}&=\lambda_-\Pi_-+\lambda_+\Pi_+,\\[4pt]
f(\mathcal{Z})&=f(\lambda_-)\Pi_-+f(\lambda_+)\Pi_+.
\end{aligned}
```

## Spacetime from the [determinant](https://en.wikipedia.org/wiki/Determinant)

For a future-directed massive particle, proper time and four-momentum are obtained from a real paravector, with four-velocity $`\mathcal{U}`$.

The spacetime increment and its determinant are:

```math
\begin{aligned}
\mathrm d\mathcal{X}&=\mathrm dt+\mathrm d𝐫\cdot𝝈,\\
\det(\mathrm d\mathcal{X})&=\mathrm dt^2-\mathrm d𝐫^2.
\end{aligned}
```

Positive [proper time](https://en.wikipedia.org/wiki/Proper_time) is chosen along the future-directed path $`\mathrm dt>0`$:

```math
\mathrm ds:=\sqrt{\mathrm dt^2-\mathrm d𝐫^2}=\sqrt{\det(\mathrm d\mathcal{X})}.
```

The proper-time rate is determined by the coordinate velocity:

```math
\begin{aligned}
\frac{\mathrm ds}{\mathrm dt}
&=\sqrt{\frac{\mathrm dt^2-\mathrm d𝐫^2}{\mathrm dt^2}}\\
&=\sqrt{1-\left(\frac{\mathrm d𝐫}{\mathrm dt}\right)^2}
=\sqrt{1-𝐯^2},
\qquad 𝐯:=\frac{\mathrm d𝐫}{\mathrm dt}.
\end{aligned}
```

The [Lorentz factor](https://en.wikipedia.org/wiki/Lorentz_factor) is:

```math
\gamma:=\frac{\mathrm dt}{\mathrm ds}=\frac{1}{\sqrt{1-𝐯^2}}.
```

The [four-velocity](https://en.wikipedia.org/wiki/Four-velocity) is:

```math
\begin{aligned}
\mathcal{U}:=\frac{\mathrm d\mathcal{X}}{\mathrm ds}
&=\frac{\mathrm dt+\mathrm d𝐫\cdot𝝈}{\mathrm ds}\\
&=\left(1+\frac{\mathrm d𝐫}{\mathrm dt}\cdot𝝈\right)
\frac{\mathrm dt}{\mathrm ds}\\
&=\gamma(1+𝐯\cdot𝝈).
\end{aligned}
```

Proper-time normalization is obtained:

```math
\det \mathcal{U}=\det\!\left(\frac{\mathrm d\mathcal{X}}{\mathrm ds}\right)
=\frac{\det(\mathrm d\mathcal{X})}{\mathrm ds^2}=1.
```

For energy $`E`$ and momentum $`𝐏`$, the [four-momentum](https://en.wikipedia.org/wiki/Four-momentum) is:

```math
\mathcal{P}:=m\mathcal{U}=E+𝐏\cdot𝝈.
```

For constant mass $`m>0`$, the determinant is scaled quadratically:

```math
\det \mathcal{P}=\det(m\mathcal{U})=m^2\det \mathcal{U}=m^2.
```

The [mass-shell](https://en.wikipedia.org/wiki/On_shell_and_off_shell#Mass_shell) relation is therefore:

```math
\det \mathcal{P}=E^2-𝐏^2=m^2.
```

## [Maxwell](https://en.wikipedia.org/wiki/Maxwell%27s_equations) in one equation

The [electric field](https://en.wikipedia.org/wiki/Electric_field) $`𝐄`$ and [magnetic field](https://en.wikipedia.org/wiki/Magnetic_field) $`𝐁`$ are combined into the complex electromagnetic field $`𝐅`$ and its paravector $`\mathcal{F}`$:

```math
𝐅:=𝐄+\mathrm{i}𝐁,\qquad \mathcal{F}:=𝐅\cdot𝝈.
```

The [charge density](https://en.wikipedia.org/wiki/Charge_density) $`\rho`$ and [current density](https://en.wikipedia.org/wiki/Current_density) $`𝐉`$ are combined into the [four-current](https://en.wikipedia.org/wiki/Four-current) $`\mathcal{J}`$:

```math
\mathcal{J}:=\rho+𝐉\cdot𝝈,\qquad \mathcal{J}^*=\rho-𝐉\cdot𝝈.
```

The [paravector derivative](https://en.wikipedia.org/wiki/Paravector#Paragradient), with [gradient](https://en.wikipedia.org/wiki/Gradient) $`\partial_{𝐫}`$, is defined by:

```math
𝝏:=\partial_t+\partial_{𝐫}\cdot𝝈,\qquad
𝝏^*=\partial_t-\partial_{𝐫}\cdot𝝈.
```

[Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations) are given by:

```math
\boxed{𝝏 \mathcal{F}=\mathcal{J}^*},
```

with components:

```math
\begin{aligned}
\text{Real scalar:}\quad &\partial_{𝐫}\cdot𝐄=\rho,\\
\text{Imaginary scalar:}\quad &\partial_{𝐫}\cdot𝐁=0,\\
\text{Real vector:}\quad &\partial_t𝐄-\partial_{𝐫}\times𝐁=-𝐉,\\
\text{Imaginary vector:}\quad &\partial_t𝐁+\partial_{𝐫}\times𝐄=0.
\end{aligned}
```

By application of $`𝝏^{*}`$, the [d'Alembertian](https://en.wikipedia.org/wiki/D%27Alembert_operator) and the [sourced wave equations](https://en.wikipedia.org/wiki/Inhomogeneous_electromagnetic_wave_equation) are obtained:

```math
\begin{aligned}
\Box&:=𝝏^{*}𝝏=\partial_t^2-\partial_{𝐫}^2,\\
𝝏^{*}(𝝏 \mathcal{F})&=\Box \mathcal{F}=𝝏^{*}\mathcal{J}^*.
\end{aligned}
```

Equivalently:

```math
\boxed{\Box \mathcal{F}^*=𝝏 \mathcal{J}}.
```

From this boxed equation, [charge conservation](https://en.wikipedia.org/wiki/Charge_conservation) and the field wave equations are obtained component by component:

```math
\begin{aligned}
\text{Real scalar:}\quad &0=\partial_t\rho+\partial_{𝐫}\cdot𝐉
&&\text{(charge conservation)},\\
\text{Imaginary scalar:}\quad &0=0,\\
\text{Real vector:}\quad &\Box𝐄=-\partial_t𝐉-\partial_{𝐫}\rho,\\
\text{Imaginary vector:}\quad &\Box𝐁=\partial_{𝐫}\times𝐉.
\end{aligned}
```

In source-free vacuum, $`\mathcal{J}=0`$ and $`\Box \mathcal{F}=0`$.

## [Electromagnetic potentials](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) and [gauge invariance](https://en.wikipedia.org/wiki/Electromagnetic_four-potential#Gauge_freedom)

The real [four-potential](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) is defined from the [scalar potential](https://en.wikipedia.org/wiki/Electric_potential) $`V`$ and [vector potential](https://en.wikipedia.org/wiki/Magnetic_vector_potential) $`𝐀`$, with four-divergence $`S`$:

```math
\mathcal{A}:=V+𝐀\cdot𝝈,
\qquad S:=\mathrm{sc}(𝝏 \mathcal{A}).
```

The field is obtained from the potential:

```math
\boxed{\mathcal{F}^{*}=𝝏 \mathcal{A}-S}.
```

In components:

```math
\begin{aligned}
\text{Real scalar:}\quad &S=\partial_tV+\partial_{𝐫}\cdot𝐀,\\
\text{Imaginary scalar:}\quad &0=0,\\
\text{Real vector:}\quad &𝐄=-\partial_t𝐀-\partial_{𝐫}V,\\
\text{Imaginary vector:}\quad &𝐁=\partial_{𝐫}\times𝐀.
\end{aligned}
```

The source equations are obtained in any [gauge](https://en.wikipedia.org/wiki/Gauge_fixing) by application of $`𝝏^*`$:

```math
\Box \mathcal{A}=𝝏^{*}(𝝏 \mathcal{A})
=𝝏^{*}(S+\mathcal{F}^{*})
=𝝏^{*}S+(𝝏 \mathcal{F})^{*}.
```

By substitution of [Maxwell’s equation](https://en.wikipedia.org/wiki/Maxwell%27s_equations):

```math
\boxed{\Box \mathcal{A}=𝝏^{*}S+\mathcal{J}}.
```

In components:

```math
\begin{aligned}
\text{Real scalar:}\quad &\Box V-\partial_tS=\rho,\\
\text{Imaginary scalar:}\quad &0=0,\\
\text{Real vector:}\quad &\Box𝐀+\partial_{𝐫}S=𝐉,\\
\text{Imaginary vector:}\quad &𝟎=𝟎.
\end{aligned}
```

Under a [gauge transformation](https://en.wikipedia.org/wiki/Electromagnetic_four-potential#Gauge_freedom), with real scalar function $`\lambda`$, the field is left unchanged:

```math
\begin{aligned}
\mathcal{A}'&:=\mathcal{A}+𝝏^{*}\lambda,\\
S'&=S+\Box\lambda,\\
\mathcal{F}'&=\mathcal{F}.
\end{aligned}
```

In the [Lorenz gauge](https://en.wikipedia.org/wiki/Lorenz_gauge_condition), the [potential wave equations](https://en.wikipedia.org/wiki/Inhomogeneous_electromagnetic_wave_equation) are obtained with $`S=0`$:

```math
\begin{aligned}
\text{Real scalar:}\quad &\Box V=\rho,\\
\text{Real vector:}\quad &\Box𝐀=𝐉.
\end{aligned}
```

## [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) from the [mass shell](https://en.wikipedia.org/wiki/On_shell_and_off_shell#Mass_shell)

Hats are used for named [energy](https://en.wikipedia.org/wiki/Energy_operator), [momentum](https://en.wikipedia.org/wiki/Momentum_operator), and [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) operators; the derivative symbols $`𝝏`$ and $`\Box`$ remain unhatted.

On a complex scalar wavefunction $`\psi`$, the canonical [energy](https://en.wikipedia.org/wiki/Energy_operator) and [momentum](https://en.wikipedia.org/wiki/Momentum_operator) operators are defined by:

```math
\begin{aligned}
\hat E&:=\mathrm{i}\partial_t,\\
\hat{𝐏}&:=-\mathrm{i}\partial_{𝐫},
\end{aligned}
```

The canonical four-momentum operator is formed by:

```math
\hat{\mathcal{P}}:=\mathrm{i}𝝏^*=\hat E+\hat{𝐏}\cdot𝝈.
```

By commutativity of the free operators:

```math
\hat{\mathcal{P}}\,\mathrm{adj}\,\hat{\mathcal{P}}
=\det\hat{\mathcal{P}}
=\hat E^2-\hat{𝐏}^2
=-\Box.
```

With the [mass shell](https://en.wikipedia.org/wiki/On_shell_and_off_shell#Mass_shell) imposed on $`\psi`$:

```math
\begin{aligned}
0&=(\det\hat{\mathcal{P}}-m^2)\psi\\
&=(\hat E^2-\hat{𝐏}^2-m^2)\psi\\
&=-(\Box+m^2)\psi.
\end{aligned}
```

Thus the free [Klein–Gordon equation](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) is:

```math
\boxed{(\Box+m^2)\psi=0.}
```

## [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) as a first-order wave equation

In the [Weyl representation](https://en.wikipedia.org/wiki/Gamma_matrices#Weyl_%28chiral%29_basis), the notation $`𝝍:=(𝝍_+,𝝍_-)^{\mathsf T}`$ is used, with two complex components in each entry. The block identity is denoted by $`𝐈`$.

The following operators are defined using the free momentum operator above:

```math
\begin{aligned}
𝐖(\hat{\mathcal{P}})&:=
\begin{pmatrix}
0 & \mathrm{adj}\,\hat{\mathcal{P}} \\
\hat{\mathcal{P}} & 0
\end{pmatrix},\\
\hat{𝐃}(m)&:=𝐖(\hat{\mathcal{P}})-m𝐈.
\end{aligned}
```

The [Dirac equation](https://en.wikipedia.org/wiki/Dirac_equation) is:

```math
\boxed{\hat{𝐃}(m)𝝍=0.}
```

Equivalently:

```math
(\hat E\pm\hat{𝐏}\cdot𝝈)𝝍_\pm=m𝝍_\mp.
```

By $`𝐖(\hat{\mathcal{P}})^2=(\det\hat{\mathcal{P}})𝐈`$, the [Klein–Gordon factorization](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) is obtained:

```math
\begin{aligned}
\hat{𝐃}(m)\hat{𝐃}(-m)&=\hat{𝐃}(-m)\hat{𝐃}(m)\\
&=𝐖(\hat{\mathcal{P}})^2-m^2𝐈
=-(\Box+m^2)𝐈.
\end{aligned}
```

Hence [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) is satisfied by every free [Dirac solution](https://en.wikipedia.org/wiki/Dirac_equation).

## [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) with an electromagnetic potential

Here, the uncoupled operators defined above are labeled by subscript $`0`$:

```math
\begin{aligned}
\hat E_0&:=\mathrm{i}\partial_t,\\
\hat{𝐏}_0&:=-\mathrm{i}\partial_{𝐫},\\
\hat{\mathcal{P}}_0&:=\mathrm{i}𝝏^*.
\end{aligned}
```

The following are assumed: constant mass $`m>0`$, [charge](https://en.wikipedia.org/wiki/Electric_charge) $`q`$, and a real [potential](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) $`\mathcal{A}=V+𝐀\cdot𝝈`$.

By [minimal coupling](https://en.wikipedia.org/wiki/Minimal_coupling), the [gauge-covariant energy](https://en.wikipedia.org/wiki/Gauge_covariant_derivative) and [kinetic momentum](https://en.wikipedia.org/wiki/Momentum_operator#Definition_(position_space)) operators are obtained:

```math
\begin{aligned}
\hat E_q&:=\hat E_0-qV,\\
\hat{𝐏}_q&:=\hat{𝐏}_0-q𝐀.
\end{aligned}
```

The [kinetic four-momentum](https://en.wikipedia.org/wiki/Minimal_coupling) operator is formed by:

```math
\hat{\mathcal{P}}_q:=\hat{\mathcal{P}}_0-q\mathcal{A}=\hat E_q+\hat{𝐏}_q\cdot𝝈.
```

In $`\hat E_q`$, rest energy is included and the potential energy $`qV`$ is excluded.

For operator arguments, the spatial vector part is reversed by $`\mathrm{adj}`$ while derivative order is preserved.

The coupled [Dirac operator](https://en.wikipedia.org/wiki/Dirac_equation) is:

```math
\begin{aligned}
\hat{𝐃}_q(m)&:=𝐖(\hat{\mathcal{P}}_q)-m𝐈\\
&=\begin{pmatrix}
-m & \mathrm{adj}\,\hat{\mathcal{P}}_q\\
\hat{\mathcal{P}}_q & -m
\end{pmatrix}.
\end{aligned}
```

The coupled [Dirac equation](https://en.wikipedia.org/wiki/Dirac_equation) is:

```math
\boxed{\hat{𝐃}_q(m)𝝍=0.}
```

Equivalently:

```math
(\hat E_q\pm\hat{𝐏}_q\cdot𝝈)𝝍_\pm=m𝝍_\mp.
```

The following field identities are satisfied by the coupled operators:

```math
\begin{aligned}
{}[\hat{𝐏}_q\cdot𝝈,\hat E_q]
&=-\mathrm{i}q𝐄\cdot𝝈,\\
(\hat{𝐏}_q\cdot𝝈)^2
&=\hat{𝐏}_q^2-q𝐁\cdot𝝈,\\
\hat{𝐏}_q^2+\mathrm{i}q\mathcal{F}
&=(\hat{𝐏}_q\cdot𝝈)^2
-[\hat{𝐏}_q\cdot𝝈,\hat E_q].
\end{aligned}
```

With derivatives applied to the potentials as well as the wavefunction, the opposite-mass product is obtained:

```math
\hat{𝐃}_q(-m)\hat{𝐃}_q(m)
=(\hat E_q^2-\hat{𝐏}_q^2-m^2)𝐈
-\mathrm{i}q\begin{pmatrix}\mathcal{F}^{*}&0\\0&\mathcal{F}\end{pmatrix}.
```

For every [Dirac solution](https://en.wikipedia.org/wiki/Dirac_equation), $`\hat{𝐃}_q(-m)\hat{𝐃}_q(m)𝝍=0`$ is satisfied. With $`q=0`$, the free [Klein–Gordon equation](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) is recovered.

## Low-energy [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) to [Schrödinger](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation)

For constant $`m>0`$, the full scalar field is written as $`\psi'=e^{-\mathrm{i}mt}\psi`$, where $`\psi`$ is the envelope with the rest-energy phase removed. By the product rule:

```math
\begin{aligned}
\hat E_q\psi'
&=(\mathrm{i}\partial_t-qV)(e^{-\mathrm{i}mt}\psi)\\
&=e^{-\mathrm{i}mt}(m\psi+\mathrm{i}\partial_t\psi-qV\psi)\\
&=e^{-\mathrm{i}mt}(m+\hat E_q)\psi,\\[6pt]
\hat{𝐏}_q\psi'&=e^{-\mathrm{i}mt}\hat{𝐏}_q\psi.
\end{aligned}
```

By substitution into [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation):

```math
\begin{aligned}
0&=(\hat E_q^2-\hat{𝐏}_q^2-m^2)\psi'\\
&=e^{-\mathrm{i}mt}\bigl((m+\hat E_q)^2-\hat{𝐏}_q^2-m^2\bigr)\psi\\
&=e^{-\mathrm{i}mt}\bigl(2m\hat E_q+\hat E_q^2-\hat{𝐏}_q^2\bigr)\psi.
\end{aligned}
```

On the envelope $`\psi`$, energy above rest, excluding $`qV`$, is measured by $`\hat E_q`$. For any $`q`$ and constant $`m>0`$, the exact recursion is obtained by canceling the phase:

```math
\hat E_q\psi=\frac{\hat{𝐏}_q^2-\hat E_q^2}{2m}\psi.
```

If $`\hat E_q`$ and $`\hat{𝐏}_q^2`$ commute, the following is obtained by iteration on the positive-energy branch:

```math
\hat E_q\psi=\left(\frac{\hat{𝐏}_q^2}{2m}
-\frac{\hat{𝐏}_q^4}{8m^3}+\cdots\right)\psi.
```

For weak, slowly varying fields, $`\hat E_q^2\psi`$ is neglected relative to $`2m\hat E_q\psi`$, and the following is obtained for [Schrödinger’s equation](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation):

```math
\boxed{\mathrm{i}\partial_t\psi=
\left(\frac{(-\mathrm{i}\partial_{𝐫}-q𝐀)^2}{2m}+qV\right)\psi.}
```

With $`𝐀=0`$ (hence $`𝐁=0`$), the following is obtained for the conventional [Schrödinger equation](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation):

```math
\mathrm{i}\partial_t\psi=
\left(-\frac{\partial_{𝐫}^2}{2m}+qV\right)\psi.
```

## Low-energy [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) to [Pauli](https://en.wikipedia.org/wiki/Pauli_equation)

The rest-energy phase is removed from the full Dirac components $`𝝍_\pm`$ by defining the envelopes $`𝝍_\pm':=e^{\mathrm{i}mt}𝝍_\pm`$. The following is obtained:

```math
(m+\hat E_q\pm\hat{𝐏}_q\cdot𝝈)𝝍_\pm'=m𝝍_\mp'.
```

The large and small components for the positive-energy branch are defined by:

```math
\begin{aligned}
𝝓_+&:=\frac{𝝍_+'+𝝍_-'}{\sqrt2},\\
𝝓_-&:=\frac{𝝍_-'-𝝍_+'}{\sqrt2}.
\end{aligned}
```

By addition and subtraction of the envelope equations:

```math
\begin{aligned}
\hat E_q𝝓_+&=(\hat{𝐏}_q\cdot𝝈)𝝓_-,\\
(2m+\hat E_q)𝝓_-&=(\hat{𝐏}_q\cdot𝝈)𝝓_+.
\end{aligned}
```

In the same nonrelativistic limit, $`\hat E_q𝝓_-`$ is neglected relative to $`2m𝝓_-`$. The resulting small component is substituted into the first equation, and the field identity above is applied:

```math
\begin{aligned}
𝝓_-&\simeq\frac{\hat{𝐏}_q\cdot𝝈}{2m}𝝓_+,\\
\hat E_q𝝓_+&\simeq\frac{(\hat{𝐏}_q\cdot𝝈)^2}{2m}𝝓_+
=\frac{\hat{𝐏}_q^2-q𝐁\cdot𝝈}{2m}𝝓_+.
\end{aligned}
```

With $`\hat E_q=\mathrm{i}\partial_t-qV`$ restored, the [Pauli equation](https://en.wikipedia.org/wiki/Pauli_equation) is obtained, with kinetic and spin terms retained through order $`1/m`$:

```math
\boxed{\mathrm{i}\partial_t𝝓_+=
\left(\frac{\hat{𝐏}_q^2}{2m}
+qV-\frac{q}{2m}𝐁\cdot𝝈\right)𝝓_+.}
```

## [Spin up and spin down](https://en.wikipedia.org/wiki/Spin-1/2#Observables) in a uniform magnetic field

For a constant nonzero field $`𝐁`$, its unit direction and [spin projectors](https://en.wikipedia.org/wiki/Pauli_matrices#Eigenvectors_and_eigenvalues) are defined by:

```math
𝐮:=\frac{𝐁}{\lVert𝐁\rVert},
\qquad \Pi_\pm:=\frac12(1\pm𝐮\cdot𝝈).
```

By [spectral decomposition](#projections-and-spectral-decomposition):

```math
𝐁\cdot𝝈
=\lVert𝐁\rVert(\Pi_+-\Pi_-),
\qquad
\Pi_\pm(𝐁\cdot𝝈)
=\pm\lVert𝐁\rVert\Pi_\pm.
```

The two spin amplitudes $`\phi_{++},\phi_{+-}`$ are encoded in a paravector $`𝝓_+`$ aligned with $`𝐮`$:

```math
\begin{aligned}
𝝓_+&=\phi_{++}\Pi_++\phi_{+-}\Pi_-,\\
\Pi_\pm𝝓_+&=𝝓_+\Pi_\pm=\phi_{+\pm}\Pi_\pm.
\end{aligned}
```

Constant $`\Pi_\pm`$ commute with $`\hat E_q`$ and $`\hat{𝐏}_q^2`$, so the following is obtained from the [Pauli equation](https://en.wikipedia.org/wiki/Pauli_equation):

```math
\begin{aligned}
0&=\Pi_\pm\left(2m\hat E_q-\hat{𝐏}_q^2
+q𝐁\cdot𝝈\right)𝝓_+\\
&=\left(\left(2m\hat E_q-\hat{𝐏}_q^2
\pm q\lVert𝐁\rVert\right)\phi_{+\pm}\right)\Pi_\pm.
\end{aligned}
```

Since $`\Pi_\pm\ne0`$, the scalar coefficients vanish. With $`\hat E_q=\mathrm{i}\partial_t-qV`$ restored:

```math
\boxed{\mathrm{i}\partial_t\phi_{+\pm}=
\left(\frac{\hat{𝐏}_q^2\mp q\lVert𝐁\rVert}{2m}+qV\right)\phi_{+\pm}.}
```

For states with [spin projections](https://en.wikipedia.org/wiki/Spin-1/2) $`\pm\tfrac12`$ along $`𝐮`$, opposite [Zeeman shifts](https://en.wikipedia.org/wiki/Zeeman_effect) are obtained. In the [kinetic momentum](https://en.wikipedia.org/wiki/Momentum_operator#Definition_(position_space)) $`\hat{𝐏}_q=-\mathrm{i}\partial_{𝐫}-q𝐀`$, the orbital coupling is retained.

This representation is used for uniform-field [Pauli evolution](https://en.wikipedia.org/wiki/Pauli_equation); the fixed spinor space is used for the full [Dirac equation](https://en.wikipedia.org/wiki/Dirac_equation) and spin measurements.

The four complex amplitudes of the [Dirac spinor](https://en.wikipedia.org/wiki/Dirac_spinor), after removal of the rest-energy phase, are arranged as:

```math
𝝓:=\left(\,\begin{matrix}
\phi_{++}&\phi_{+-}\\
\phi_{-+}&\phi_{--}
\end{matrix}\,\right).
```

[Dirac blocks](https://en.wikipedia.org/wiki/Dirac_spinor) are indexed by rows and spin projections by columns; the upper row is governed by the [Pauli limit](https://en.wikipedia.org/wiki/Pauli_equation).

## [Spin measurement](https://en.wikipedia.org/wiki/Spin-1/2#Rotations_and_Spinors) probabilities

Real unit preparation and measurement axes are denoted by $`𝐮',𝐮`$, with $`\cos\theta:=𝐮\cdot𝐮'`$.

The state $`𝝍=\Pi_+'𝝍\ne0`$ is prepared using projectors with the following overlap:

```math
\begin{gathered}
\Pi_+':=\frac12(1+𝐮'\cdot𝝈),\qquad
\Pi_\pm:=\frac12(1\pm𝐮\cdot𝝈),\\[4pt]
\Pi_+'\Pi_\pm\Pi_+'=\frac{1\pm\cos\theta}{2}\Pi_+'.
\end{gathered}
```

Using the squared norm, $`\lVert \mathcal{Z}\rVert^2=\mathrm{sc}(\mathcal{Z}\mathcal{Z}^{\mathsf H})`$, the denominator is:

```math
\lVert𝝍\rVert^2=\mathrm{sc}(𝝍𝝍^{\mathsf H})>0.
```

For the [prepared state](https://en.wikipedia.org/wiki/Spin-1/2#Bloch_Representation), the projected squared norm is:

```math
\begin{aligned}
\lVert\Pi_\pm𝝍\rVert^2
&=\mathrm{sc}\bigl((\Pi_\pm𝝍)(\Pi_\pm𝝍)^{\mathsf H}\bigr)\\
&=\mathrm{sc}(\Pi_\pm𝝍𝝍^{\mathsf H})\\
&=\mathrm{sc}(\Pi_+'\Pi_\pm\Pi_+'𝝍𝝍^{\mathsf H})\\
&=\frac{1\pm\cos\theta}{2}\mathrm{sc}(𝝍𝝍^{\mathsf H})
=\frac{1\pm\cos\theta}{2}\lVert𝝍\rVert^2.
\end{aligned}
```

The [Born probabilities](https://en.wikipedia.org/wiki/Born_rule) for [spin projections](https://en.wikipedia.org/wiki/Spin-1/2) $`\pm\tfrac12`$ along $`𝐮`$ are:

```math
\begin{gathered}
\mathrm{prob}(\pm)=\frac{\lVert\Pi_\pm𝝍\rVert^2}{\lVert𝝍\rVert^2}
=\frac{1\pm\cos\theta}{2},\\[6pt]
\boxed{\mathrm{prob}(+)=\cos^2\frac\theta2,\quad \mathrm{prob}(-)=\sin^2\frac\theta2.}
\end{gathered}
```

For aligned axes, spin up is obtained with certainty; for perpendicular axes, equal probabilities are obtained. The probabilities are exchanged by reversal of the prepared spin.

## Conventions

Alternative names and comparisons with the literature are collected in [Appendix C: Conventions](./sta_notes.pdf#nameddest=appendix.C) of the paper.
