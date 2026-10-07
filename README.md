# Space–Time Algebra

**A common [paravector](https://en.wikipedia.org/wiki/Paravector) language is provided by [complex scalars](https://en.wikipedia.org/wiki/Complex_number) and [three-vectors](https://en.wikipedia.org/wiki/Euclidean_vector) for [quaternions](https://en.wikipedia.org/wiki/Quaternion), [rotations](https://en.wikipedia.org/wiki/Rotation_%28mathematics%29), [Lorentz boosts](https://en.wikipedia.org/wiki/Lorentz_transformation), [projections](https://en.wikipedia.org/wiki/Paravector#Null_paravectors_as_projectors), [spacetime geometry](https://en.wikipedia.org/wiki/Minkowski_space), [relativistic particle dynamics](https://en.wikipedia.org/wiki/Relativistic_mechanics), [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations), [electromagnetic waves](https://en.wikipedia.org/wiki/Electromagnetic_radiation) and [forces](https://en.wikipedia.org/wiki/Lorentz_force), [gauge coupling](https://en.wikipedia.org/wiki/Minimal_coupling), the [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) and [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) equations, their [Schrödinger](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation) and [Pauli](https://en.wikipedia.org/wiki/Pauli_equation) low-energy limits, and [spin-one-half](https://en.wikipedia.org/wiki/Spin-1/2) [amplitudes](https://en.wikipedia.org/wiki/Probability_amplitude) and [measurement probabilities](https://en.wikipedia.org/wiki/Born_rule).**

[Paper (PDF)](./sta_notes.pdf) · [LaTeX source](./sta_notes.tex)

[Natural units](https://en.wikipedia.org/wiki/Natural_units) $`\hbar=c=1`$, [rationalized electromagnetic units](https://en.wikipedia.org/wiki/Heaviside%E2%80%93Lorentz_units), and [signature](https://en.wikipedia.org/wiki/Metric_signature) $`(+,-,-,-)`$ are used.

Paravectors use ordinary, nonbold math letters ($`X`$, $`Y`$, $`Z`$); vectors and matrices are bold math italic ($`\boldsymbol{X}`$, $`\boldsymbol{r}`$, $`\boldsymbol{M}`$). Greek paravectors are nonbold ($`\psi`$, $`\Pi`$), while vector or matrix symbols are bold ($`\boldsymbol{\sigma}`$, $`\boldsymbol{\psi}`$). Thus $`P=E+\boldsymbol{P}\cdot\boldsymbol{\sigma}`$ distinguishes the paravector momentum from its vector part.

## One complex [paravector](https://en.wikipedia.org/wiki/Paravector)

Complex paravectors are formed from sigmas—algebraic objects with order-dependent multiplication:

```math
\boxed{
\begin{aligned}
\sigma_n^2&=1,
&\sigma_n\sigma_m&=-\sigma_m\sigma_n\quad(n\ne m),\\
\mathrm i&:=\sigma_1\sigma_2\sigma_3,
&\mathrm i^2&=-1,\\
\boldsymbol{\sigma}&:=\{\sigma_1,\sigma_2,\sigma_3\},
&\mathrm i\boldsymbol{\sigma}&=\{\sigma_2\sigma_3,\sigma_3\sigma_1,\sigma_1\sigma_2\},\\
{Z}&:=(a+b\mathrm i)+(\boldsymbol{X}+\mathrm i\boldsymbol{Y})\cdot\boldsymbol{\sigma},
&a,b&\in\mathbb R,\quad\boldsymbol{X},\boldsymbol{Y}\in\mathbb R^3.
\end{aligned}}
```

The eight real components are grouped into a complex scalar and vector:

```math
\begin{aligned}
{X}&=a+\boldsymbol{X}\cdot\boldsymbol{\sigma},\\
{Y}&=b+\boldsymbol{Y}\cdot\boldsymbol{\sigma},\\[4pt]
{Z}&={X}+\mathrm i{Y}\\
&=(a+b\mathrm i)+(\boldsymbol{X}+\mathrm i\boldsymbol{Y})\cdot\boldsymbol{\sigma}\\
&=c+\boldsymbol{Z}\cdot\boldsymbol{\sigma}.
\end{aligned}
```

The component extractions are:

```math
\begin{aligned}
\mathrm{re}({Z})&:={X},\\
\mathrm{im}({Z})&:={Y},\\
\mathrm{sc}({Z})&:=c=a+\mathrm ib,\\
\mathrm{vec}({Z})&:=\boldsymbol{Z}=\boldsymbol{X}+\mathrm i\boldsymbol{Y}.
\end{aligned}
```

## Product

Both the [dot](https://en.wikipedia.org/wiki/Dot_product) and [cross](https://en.wikipedia.org/wiki/Cross_product) products are included in the product:

```math
(\boldsymbol{X}\cdot\boldsymbol{\sigma})(\boldsymbol{Y}\cdot\boldsymbol{\sigma})
=\boldsymbol{X}\cdot\boldsymbol{Y}+\mathrm i(\boldsymbol{X}\times\boldsymbol{Y})\cdot\boldsymbol{\sigma}.
```

The convention $`\boldsymbol{X}^2:=\boldsymbol{X}\cdot\boldsymbol{X}`$ is used.

The [commutator](https://en.wikipedia.org/wiki/Commutator) is:

```math
[{X},{Y}]:={X}{Y}-{Y}{X}.
```

```math
[\boldsymbol{X}\cdot\boldsymbol{\sigma},\boldsymbol{Y}\cdot\boldsymbol{\sigma}]=2\mathrm i(\boldsymbol{X}\times\boldsymbol{Y})\cdot\boldsymbol{\sigma}.
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

Both are applied termwise; real coefficients are unchanged:

```math
\begin{aligned}
(\sigma_m\sigma_n)^{\mathsf H}&=\sigma_n\sigma_m,\\[4pt]
(\sigma_m\sigma_n)^{*}&=(-\sigma_m)(-\sigma_n).
\end{aligned}
```

With $`\mathrm i=\sigma_1\sigma_2\sigma_3`$:

```math
\begin{aligned}
\mathrm i^{\mathsf H}&=\sigma_3\sigma_2\sigma_1=-\mathrm i,\\[4pt]
\mathrm i^{*}&=(-\sigma_1)(-\sigma_2)(-\sigma_3)=-\mathrm i.
\end{aligned}
```

For $`{Z}=c+\boldsymbol{Z}\cdot\boldsymbol{\sigma}`$, $`c\in\mathbb C`$, $`\boldsymbol{Z}\in\mathbb C^3`$, the complex coefficients are conjugated:

```math
\begin{aligned}
{Z}^{\mathsf H}&=c^{*}+\boldsymbol{Z}^{*}\cdot\boldsymbol{\sigma},\\[4pt]
{Z}^{*}&=c^{*}-\boldsymbol{Z}^{*}\cdot\boldsymbol{\sigma}.
\end{aligned}
```

## Core operations

The paravector [adjugate](https://en.wikipedia.org/wiki/Adjugate_matrix), [determinant](https://en.wikipedia.org/wiki/Determinant), and [inverse](https://en.wikipedia.org/wiki/Invertible_matrix) are given by:

```math
\begin{aligned}
\mathrm{adj}\,{Z}&:=c-\boldsymbol{Z}\cdot\boldsymbol{\sigma}={Z}^{*\mathsf H},\\[4pt]
\det {Z}&:=c^2-\boldsymbol{Z}^2={Z}\,\mathrm{adj}\,{Z},\\[6pt]
{Z}^{-1}&=\frac{\mathrm{adj}\,{Z}}{\det {Z}}\qquad\text{if }\det {Z}\ne0.
\end{aligned}
```

The paravector [trace](https://en.wikipedia.org/wiki/Trace_%28linear_algebra%29) and [squared norm](https://en.wikipedia.org/wiki/Norm_%28mathematics%29) are given by:

```math
\begin{aligned}
\mathrm{tr}\,{Z}&:=2c={Z}+\mathrm{adj}\,{Z},\\[4pt]
\lVert {Z}\rVert^2&:=cc^{*}+\boldsymbol{Z}\cdot\boldsymbol{Z}^{*}=\mathrm{sc}({Z}{Z}^{\mathsf H}).
\end{aligned}
```

## [Quaternions](https://en.wikipedia.org/wiki/Quaternion)

The quaternion multiplication rules are satisfied by $`(-\mathrm i\sigma_k)`$:

```math
\begin{gathered}
(-\mathrm i\sigma_1)^2=(-\mathrm i\sigma_2)^2=(-\mathrm i\sigma_3)^2
=(-\mathrm i\sigma_1)(-\mathrm i\sigma_2)(-\mathrm i\sigma_3)=-1,\\[1em]
(-\mathrm i\sigma_1)(-\mathrm i\sigma_2)=(-\mathrm i\sigma_3),\\
(-\mathrm i\sigma_2)(-\mathrm i\sigma_3)=(-\mathrm i\sigma_1),\\
(-\mathrm i\sigma_3)(-\mathrm i\sigma_1)=(-\mathrm i\sigma_2).
\end{gathered}
```

For $`a\in\mathbb R`$ and $`\boldsymbol{Q}\in\mathbb R^3`$, the quaternion and its conjugate are:

```math
\begin{aligned}
{Q}&:=a+\boldsymbol{Q}\cdot(-\mathrm i\boldsymbol{\sigma})=a-\mathrm i\boldsymbol{Q}\cdot\boldsymbol{\sigma},\\
{Q}^{\mathsf H}&:=a-\boldsymbol{Q}\cdot(-\mathrm i\boldsymbol{\sigma})=a+\mathrm i\boldsymbol{Q}\cdot\boldsymbol{\sigma},
\end{aligned}
```

with

```math
{Q}{Q}^{\mathsf H}=\det {Q}=a^2+\boldsymbol{Q}^2.
```

Nonzero quaternions are normalized with:

```math
{Q}\leftarrow\frac{{Q}}{\sqrt{{Q}{Q}^{\mathsf H}}},
```

Thus, with $`\lVert\boldsymbol{Q}\rVert^2=\boldsymbol{Q}^2`$ for real $`\boldsymbol{Q}`$, the normalization is expressed as:

```math
{Q}{Q}^{\mathsf H}=a^2+\lVert\boldsymbol{Q}\rVert^2=1.
```

For $`\boldsymbol{Q}\ne0`$, the rotation axis and angle are specified by:

```math
\boldsymbol{u}=\frac{\boldsymbol{Q}}{\lVert\boldsymbol{Q}\rVert},\qquad
a=\cos\frac{\theta}{2},\qquad
\boldsymbol{Q}=\boldsymbol{u}\sin\frac{\theta}{2}.
```

A rotation of $`\boldsymbol{r}\in\mathbb R^3`$ about $`\boldsymbol{u}`$ by angle $`\theta`$ is expressed by:

```math
{X}'=\boldsymbol{r}'\cdot(-\mathrm i\boldsymbol{\sigma}), \qquad
{X}=\boldsymbol{r}\cdot(-\mathrm i\boldsymbol{\sigma}), \qquad
{X}'={Q}{X}{Q}^{\mathsf H}.
```

## Rotations and boosts

For real spacetime $`{X}`$ and a real [unit axis](https://en.wikipedia.org/wiki/Unit_vector) $`\boldsymbol{u}`$, the transformation is written as:

```math
{X}:=t+\boldsymbol{r}\cdot\boldsymbol{\sigma},
\qquad {X}':={T}{X}{T}^{\mathsf H}=t'+\boldsymbol{r}'\cdot\boldsymbol{\sigma}.
```

For $`{T}={R}`$, a [Rodrigues rotation](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula) through $`\theta`$, with [unit quaternion](https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation) $`{Q}={R}`$, is given by:

```math
\begin{aligned}
{R}&:=e^{-\mathrm i\theta\boldsymbol{u}\cdot\boldsymbol{\sigma}/2}
=\cos\frac\theta2-\mathrm i\boldsymbol{u}\cdot\boldsymbol{\sigma}\sin\frac\theta2,\\
{X}'&={R}{X}{R}^{\mathsf H},\\
t'&=t,\\
\boldsymbol{r}'&=\boldsymbol{r}\cos\theta+(\boldsymbol{u}\times\boldsymbol{r})\sin\theta
+(\boldsymbol{u}\cdot\boldsymbol{r})(1-\cos\theta)\boldsymbol{u}.
\end{aligned}
```

For $`{T}={L}`$, a [Lorentz boost](https://en.wikipedia.org/wiki/Lorentz_transformation) to a frame moving at $`+\beta\boldsymbol{u}`$, with [rapidity](https://en.wikipedia.org/wiki/Rapidity) $`\theta`$, $`\beta:=\tanh\theta`$, and [Lorentz factor](https://en.wikipedia.org/wiki/Lorentz_factor) $`\gamma:=\cosh\theta`$ is considered.

The parallel and perpendicular components are defined by:

```math
\boldsymbol{r}_\parallel:=(\boldsymbol{u}\cdot\boldsymbol{r})\boldsymbol{u},
\qquad \boldsymbol{r}_\perp:=\boldsymbol{r}-\boldsymbol{r}_\parallel.
```

The boost and transformed coordinates are:

```math
\begin{aligned}
{L}&:=e^{-\theta\boldsymbol{u}\cdot\boldsymbol{\sigma}/2}
=\cosh\frac\theta2-\boldsymbol{u}\cdot\boldsymbol{\sigma}\sinh\frac\theta2
={L}^{\mathsf H},\\
{X}'&={L}{X}{L}^{\mathsf H},\\
t'&=\gamma(t-\beta\boldsymbol{u}\cdot\boldsymbol{r}),\\
\boldsymbol{r}'&=\boldsymbol{r}_\perp+\gamma(\boldsymbol{r}_\parallel-\beta t\boldsymbol{u}).
\end{aligned}
```

## Projections and spectral decomposition

A paravector of the following form is considered:

```math
{Z}=c+z\boldsymbol{u}\cdot\boldsymbol{\sigma},
\qquad c,z\in\mathbb C,\quad \boldsymbol{u}\in\mathbb C^3,\quad \boldsymbol{u}^2=1.
```

The eigenvalues and complementary [projectors](https://en.wikipedia.org/wiki/Paravector#Null_paravectors_as_projectors) are:

```math
\begin{gathered}
\lambda_\pm:=c\pm z,
\qquad \Pi _\pm:=\frac12(1\pm\boldsymbol{u}\cdot\boldsymbol{\sigma}),\\[4pt]
\Pi _\pm^2=\Pi _\pm,\qquad \Pi _+\Pi _-=0,\qquad \Pi _-+\Pi _+=1.
\end{gathered}
```

For real $`\boldsymbol{u}`$, these are also Hermitian projectors. For $`f`$ analytic near $`\lambda_\pm`$, the following spectral decomposition is obtained:

```math
\begin{aligned}
{Z}&=\lambda_-\Pi _-+\lambda_+\Pi _+,\\[4pt]
f({Z})&=f(\lambda_-)\Pi _-+f(\lambda_+)\Pi _+.
\end{aligned}
```

## Spacetime from the [determinant](https://en.wikipedia.org/wiki/Determinant)

In natural units, inertial coordinate time is denoted by $`t`$ and onboard proper time (often $`\tau`$) by $`s`$.

The spacetime increment and interval are:

```math
\begin{aligned}
\mathrm d{X}&=\mathrm dt+\mathrm d\boldsymbol{r}\cdot\boldsymbol{\sigma},\\
\det(\mathrm d{X})&=\mathrm dt^2-\mathrm d\boldsymbol{r}^2.
\end{aligned}
```

For a future-directed time-like path ($`\mathrm dt>0`$), positive [proper time](https://en.wikipedia.org/wiki/Proper_time) is defined by:

```math
\mathrm ds:=\sqrt{\mathrm dt^2-\mathrm d\boldsymbol{r}^2}=\sqrt{\det(\mathrm d{X})}.
```

Dividing by coordinate time gives:

```math
\begin{aligned}
\frac{\mathrm ds}{\mathrm dt}
&=\sqrt{\frac{\mathrm dt^2-\mathrm d\boldsymbol{r}^2}{\mathrm dt^2}}\\
&=\sqrt{1-\left(\frac{\mathrm d\boldsymbol{r}}{\mathrm dt}\right)^2}
=\sqrt{1-\boldsymbol{v}^2},
\qquad \boldsymbol{v}:=\frac{\mathrm d\boldsymbol{r}}{\mathrm dt}.
\end{aligned}
```

The [Lorentz factor](https://en.wikipedia.org/wiki/Lorentz_factor) is:

```math
\gamma:=\frac{\mathrm dt}{\mathrm ds}=\frac{1}{\sqrt{1-\boldsymbol{v}^2}}.
```

The [four-velocity](https://en.wikipedia.org/wiki/Four-velocity) is:

```math
\begin{aligned}
{U}:=\frac{\mathrm d{X}}{\mathrm ds}
&=\frac{\mathrm dt+\mathrm d\boldsymbol{r}\cdot\boldsymbol{\sigma}}{\mathrm ds}\\
&=\left(1+\frac{\mathrm d\boldsymbol{r}}{\mathrm dt}\cdot\boldsymbol{\sigma}\right)
\frac{\mathrm dt}{\mathrm ds}\\
&=\gamma(1+\boldsymbol{v}\cdot\boldsymbol{\sigma}).
\end{aligned}
```

Thus:

```math
\det {U}=\det\!\left(\frac{\mathrm d{X}}{\mathrm ds}\right)
=\frac{\det(\mathrm d{X})}{\mathrm ds^2}=1.
```

For energy $`E`$ and momentum $`\boldsymbol{P}`$, the [four-momentum](https://en.wikipedia.org/wiki/Four-momentum) is:

```math
{P}:=m{U}=E+\boldsymbol{P}\cdot\boldsymbol{\sigma}.
```

For constant mass $`m>0`$:

```math
\det {P}=\det(m{U})=m^2\det {U}=m^2.
```

The [mass-shell](https://en.wikipedia.org/wiki/On_shell_and_off_shell#Mass_shell) relation follows:

```math
\det {P}=E^2-\boldsymbol{P}^2=m^2.
```

## [Maxwell](https://en.wikipedia.org/wiki/Maxwell%27s_equations) in one equation

The [electric field](https://en.wikipedia.org/wiki/Electric_field) $`\boldsymbol{E}`$ and [magnetic field](https://en.wikipedia.org/wiki/Magnetic_field) $`\boldsymbol{B}`$ are combined into the complex electromagnetic field $`\boldsymbol{F}`$ and its paravector $`{F}`$:

```math
\boldsymbol{F}:=\boldsymbol{E}+\mathrm i\boldsymbol{B},\qquad {F}:=\boldsymbol{F}\cdot\boldsymbol{\sigma}.
```

The [charge density](https://en.wikipedia.org/wiki/Charge_density) $`\rho`$ and [current density](https://en.wikipedia.org/wiki/Current_density) $`\boldsymbol{J}`$ are combined into the [four-current](https://en.wikipedia.org/wiki/Four-current) $`{J}`$:

```math
{J}:=\rho+\boldsymbol{J}\cdot\boldsymbol{\sigma},\qquad {J}^*=\rho-\boldsymbol{J}\cdot\boldsymbol{\sigma}.
```

The [paravector derivative](https://en.wikipedia.org/wiki/Paravector#Paragradient), with [gradient](https://en.wikipedia.org/wiki/Gradient) $`\partial_{\boldsymbol{r}}`$, is defined by:

```math
\partial :=\partial_t+\partial_{\boldsymbol{r}}\cdot\boldsymbol{\sigma},\qquad
\partial ^*=\partial_t-\partial_{\boldsymbol{r}}\cdot\boldsymbol{\sigma}.
```

[Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations) are given by:

```math
\boxed{\partial  {F}={J}^*},
```

with components:

```math
\begin{aligned}
\text{Real scalar:}\quad &\partial_{\boldsymbol{r}}\cdot\boldsymbol{E}=\rho,\\
\text{Imaginary scalar:}\quad &\partial_{\boldsymbol{r}}\cdot\boldsymbol{B}=0,\\
\text{Real vector:}\quad &\partial_t\boldsymbol{E}-\partial_{\boldsymbol{r}}\times\boldsymbol{B}=-\boldsymbol{J},\\
\text{Imaginary vector:}\quad &\partial_t\boldsymbol{B}+\partial_{\boldsymbol{r}}\times\boldsymbol{E}=0.
\end{aligned}
```

By application of $`\partial ^{*}`$, the [d'Alembertian](https://en.wikipedia.org/wiki/D%27Alembert_operator) and the [sourced wave equations](https://en.wikipedia.org/wiki/Inhomogeneous_electromagnetic_wave_equation) are obtained:

```math
\begin{aligned}
\Box&:=\partial ^{*}\partial =\partial_t^2-\partial_{\boldsymbol{r}}^2,\\
\partial ^{*}(\partial  {F})&=\Box {F}=\partial ^{*}{J}^*.
\end{aligned}
```

Equivalently:

```math
\boxed{\Box {F}^*=\partial  {J}}.
```

From this boxed equation, [charge conservation](https://en.wikipedia.org/wiki/Charge_conservation) and the field wave equations are obtained component by component:

```math
\begin{aligned}
\text{Real scalar:}\quad &0=\partial_t\rho+\partial_{\boldsymbol{r}}\cdot\boldsymbol{J}
&&\text{(charge conservation)},\\
\text{Imaginary scalar:}\quad &0=0,\\
\text{Real vector:}\quad &\Box\boldsymbol{E}=-\partial_t\boldsymbol{J}-\partial_{\boldsymbol{r}}\rho,\\
\text{Imaginary vector:}\quad &\Box\boldsymbol{B}=\partial_{\boldsymbol{r}}\times\boldsymbol{J}.
\end{aligned}
```

In source-free vacuum, $`{J}=0`$ and $`\Box {F}=0`$.

## [Electromagnetic potentials](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) and [gauge invariance](https://en.wikipedia.org/wiki/Electromagnetic_four-potential#Gauge_freedom)

The real [four-potential](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) is defined from the [scalar potential](https://en.wikipedia.org/wiki/Electric_potential) $`V`$ and [vector potential](https://en.wikipedia.org/wiki/Magnetic_vector_potential) $`\boldsymbol{A}`$, with four-divergence $`S`$:

```math
{A}:=V+\boldsymbol{A}\cdot\boldsymbol{\sigma},
\qquad S:=\mathrm{sc}(\partial  {A}).
```

The field is obtained from the potential:

```math
\boxed{{F}^{*}=\partial  {A}-S}.
```

In components:

```math
\begin{aligned}
\text{Real scalar:}\quad &S=\partial_tV+\partial_{\boldsymbol{r}}\cdot\boldsymbol{A},\\
\text{Imaginary scalar:}\quad &0=0,\\
\text{Real vector:}\quad &\boldsymbol{E}=-\partial_t\boldsymbol{A}-\partial_{\boldsymbol{r}}V,\\
\text{Imaginary vector:}\quad &\boldsymbol{B}=\partial_{\boldsymbol{r}}\times\boldsymbol{A}.
\end{aligned}
```

The source equations are obtained in any [gauge](https://en.wikipedia.org/wiki/Gauge_fixing) by application of $`\partial ^*`$:

```math
\Box {A}=\partial ^{*}(\partial  {A})
=\partial ^{*}(S+{F}^{*})
=\partial ^{*}S+(\partial  {F})^{*}.
```

By substitution of [Maxwell’s equation](https://en.wikipedia.org/wiki/Maxwell%27s_equations):

```math
\boxed{\Box {A}=\partial ^{*}S+{J}}.
```

In components:

```math
\begin{aligned}
\text{Real scalar:}\quad &\Box V-\partial_tS=\rho,\\
\text{Imaginary scalar:}\quad &0=0,\\
\text{Real vector:}\quad &\Box\boldsymbol{A}+\partial_{\boldsymbol{r}}S=\boldsymbol{J},\\
\text{Imaginary vector:}\quad &\boldsymbol{0}=\boldsymbol{0}.
\end{aligned}
```

Under a [gauge transformation](https://en.wikipedia.org/wiki/Electromagnetic_four-potential#Gauge_freedom), with real scalar function $`\lambda`$, the field is left unchanged:

```math
\begin{aligned}
{A}'&:={A}+\partial ^{*}\lambda,\\
S'&=S+\Box\lambda,\\
{F}'&={F}.
\end{aligned}
```

In the [Lorenz gauge](https://en.wikipedia.org/wiki/Lorenz_gauge_condition), the [potential wave equations](https://en.wikipedia.org/wiki/Inhomogeneous_electromagnetic_wave_equation) are obtained with $`S=0`$:

```math
\begin{aligned}
\text{Real scalar:}\quad &\Box V=\rho,\\
\text{Real vector:}\quad &\Box\boldsymbol{A}=\boldsymbol{J}.
\end{aligned}
```

## [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) from the [mass shell](https://en.wikipedia.org/wiki/On_shell_and_off_shell#Mass_shell)

Hats are used for named [energy](https://en.wikipedia.org/wiki/Energy_operator), [momentum](https://en.wikipedia.org/wiki/Momentum_operator), and [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) operators; the derivative symbols $`\partial `$ and $`\Box`$ remain unhatted.

On a complex scalar wavefunction $`\psi`$, the canonical [energy](https://en.wikipedia.org/wiki/Energy_operator) and [momentum](https://en.wikipedia.org/wiki/Momentum_operator) operators are defined by:

```math
\begin{aligned}
\hat E&:=\mathrm i\partial_t,\\
\hat{\boldsymbol{P}}&:=-\mathrm i\partial_{\boldsymbol{r}},
\end{aligned}
```

The canonical four-momentum operator is formed by:

```math
\hat{P}:=\mathrm i\partial ^*=\hat E+\hat{\boldsymbol{P}}\cdot\boldsymbol{\sigma}.
```

By commutativity of the free operators:

```math
\hat{P}\,\mathrm{adj}\,\hat{P}
=\det\hat{P}
=\hat E^2-\hat{\boldsymbol{P}}^2
=-\Box.
```

With the [mass shell](https://en.wikipedia.org/wiki/On_shell_and_off_shell#Mass_shell) imposed on $`\psi`$:

```math
\begin{aligned}
0&=(\det\hat{P}-m^2)\psi\\
&=(\hat E^2-\hat{\boldsymbol{P}}^2-m^2)\psi\\
&=-(\Box+m^2)\psi.
\end{aligned}
```

Thus the free [Klein–Gordon equation](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) is:

```math
\boxed{(\Box+m^2)\psi=0.}
```

## [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) as a first-order wave equation

In the [Weyl representation](https://en.wikipedia.org/wiki/Gamma_matrices#Weyl_%28chiral%29_basis), the notation $`\boldsymbol{\psi}:=\{\psi _+,\psi _-\}^{\mathsf T}`$ is used, with two complex components in each entry. The block identity is denoted by $`\boldsymbol{I}`$.

The following operators are defined using the free momentum operator above:

```math
\begin{aligned}
\boldsymbol{W}(\hat{P})&:=
\begin{Bmatrix}
0 & \mathrm{adj}\,\hat{P} \\
\hat{P} & 0
\end{Bmatrix},\\
\hat{\boldsymbol{D}}(m)&:=\boldsymbol{W}(\hat{P})-m\boldsymbol{I}.
\end{aligned}
```

The [Dirac equation](https://en.wikipedia.org/wiki/Dirac_equation) is:

```math
\boxed{\hat{\boldsymbol{D}}(m)\boldsymbol{\psi}=0.}
```

Equivalently:

```math
(\hat E\pm\hat{\boldsymbol{P}}\cdot\boldsymbol{\sigma})\psi _\pm=m\psi _\mp.
```

By $`\boldsymbol{W}(\hat{P})^2=(\det\hat{P})\boldsymbol{I}`$, the [Klein–Gordon factorization](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) is obtained:

```math
\begin{aligned}
\hat{\boldsymbol{D}}(m)\hat{\boldsymbol{D}}(-m)&=\hat{\boldsymbol{D}}(-m)\hat{\boldsymbol{D}}(m)\\
&=\boldsymbol{W}(\hat{P})^2-m^2\boldsymbol{I}
=-(\Box+m^2)\boldsymbol{I}.
\end{aligned}
```

Hence [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) is satisfied by every free [Dirac solution](https://en.wikipedia.org/wiki/Dirac_equation).

## [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) with an electromagnetic potential

Here, the uncoupled operators defined above are labeled by subscript $`0`$:

```math
\begin{aligned}
\hat E_0&:=\mathrm i\partial_t,\\
\hat{\boldsymbol{P}}_0&:=-\mathrm i\partial_{\boldsymbol{r}},\\
\hat{P}_0&:=\mathrm i\partial ^*.
\end{aligned}
```

The following are assumed: constant mass $`m>0`$, [charge](https://en.wikipedia.org/wiki/Electric_charge) $`q`$, and a real [potential](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) $`{A}=V+\boldsymbol{A}\cdot\boldsymbol{\sigma}`$.

By [minimal coupling](https://en.wikipedia.org/wiki/Minimal_coupling), the [gauge-covariant energy](https://en.wikipedia.org/wiki/Gauge_covariant_derivative) and [kinetic momentum](https://en.wikipedia.org/wiki/Momentum_operator#Definition_(position_space)) operators are obtained:

```math
\begin{aligned}
\hat E_q&:=\hat E_0-qV,\\
\hat{\boldsymbol{P}}_q&:=\hat{\boldsymbol{P}}_0-q\boldsymbol{A}.
\end{aligned}
```

The [kinetic four-momentum](https://en.wikipedia.org/wiki/Minimal_coupling) operator is formed by:

```math
\hat{P}_q:=\hat{P}_0-q{A}=\hat E_q+\hat{\boldsymbol{P}}_q\cdot\boldsymbol{\sigma}.
```

In $`\hat E_q`$, rest energy is included and the potential energy $`qV`$ is excluded.

For operator arguments, the spatial vector part is reversed by $`\mathrm{adj}`$ while derivative order is preserved.

The coupled [Dirac operator](https://en.wikipedia.org/wiki/Dirac_equation) is:

```math
\begin{aligned}
\hat{\boldsymbol{D}}_q(m)&:=\boldsymbol{W}(\hat{P}_q)-m\boldsymbol{I}\\
&=\begin{Bmatrix}
-m & \mathrm{adj}\,\hat{P}_q\\
\hat{P}_q & -m
\end{Bmatrix}.
\end{aligned}
```

The coupled [Dirac equation](https://en.wikipedia.org/wiki/Dirac_equation) is:

```math
\boxed{\hat{\boldsymbol{D}}_q(m)\boldsymbol{\psi}=0.}
```

Equivalently:

```math
(\hat E_q\pm\hat{\boldsymbol{P}}_q\cdot\boldsymbol{\sigma})\psi _\pm=m\psi _\mp.
```

The following field identities are satisfied by the coupled operators:

```math
\begin{aligned}
{}[\hat{\boldsymbol{P}}_q\cdot\boldsymbol{\sigma},\hat E_q]
&=-\mathrm iq\boldsymbol{E}\cdot\boldsymbol{\sigma},\\
(\hat{\boldsymbol{P}}_q\cdot\boldsymbol{\sigma})^2
&=\hat{\boldsymbol{P}}_q^2-q\boldsymbol{B}\cdot\boldsymbol{\sigma},\\
\hat{\boldsymbol{P}}_q^2+\mathrm iq{F}
&=(\hat{\boldsymbol{P}}_q\cdot\boldsymbol{\sigma})^2
-[\hat{\boldsymbol{P}}_q\cdot\boldsymbol{\sigma},\hat E_q].
\end{aligned}
```

With derivatives applied to the potentials as well as the wavefunction, the opposite-mass product is obtained:

```math
\hat{\boldsymbol{D}}_q(-m)\hat{\boldsymbol{D}}_q(m)
=(\hat E_q^2-\hat{\boldsymbol{P}}_q^2-m^2)\boldsymbol{I}
-\mathrm iq\begin{Bmatrix}{F}^{*}&0\\0&{F}\end{Bmatrix}.
```

For every [Dirac solution](https://en.wikipedia.org/wiki/Dirac_equation), $`\hat{\boldsymbol{D}}_q(-m)\hat{\boldsymbol{D}}_q(m)\boldsymbol{\psi}=0`$ is satisfied. With $`q=0`$, the free [Klein–Gordon equation](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) is recovered.

## Low-energy [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) to [Schrödinger](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation)

For constant $`m>0`$, the rest-energy phase is factored out as $`\psi'=e^{-\mathrm imt}\psi`$. By the product rule:

```math
\begin{aligned}
\hat E_q\psi'
&=(\mathrm i\partial_t-qV)(e^{-\mathrm imt}\psi)\\
&=e^{-\mathrm imt}(m\psi+\mathrm i\partial_t\psi-qV\psi)\\
&=e^{-\mathrm imt}(m+\hat E_q)\psi,\\[6pt]
\hat{\boldsymbol{P}}_q\psi'&=e^{-\mathrm imt}\hat{\boldsymbol{P}}_q\psi.
\end{aligned}
```

By substitution into [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation):

```math
\begin{aligned}
0&=(\hat E_q^2-\hat{\boldsymbol{P}}_q^2-m^2)\psi'\\
&=e^{-\mathrm imt}\bigl((m+\hat E_q)^2-\hat{\boldsymbol{P}}_q^2-m^2\bigr)\psi\\
&=e^{-\mathrm imt}\bigl(2m\hat E_q+\hat E_q^2-\hat{\boldsymbol{P}}_q^2\bigr)\psi.
\end{aligned}
```

On the envelope $`\psi`$, energy above rest, excluding $`qV`$, is measured by $`\hat E_q`$. After phase cancellation:

```math
\hat E_q\psi=\frac{\hat{\boldsymbol{P}}_q^2-\hat E_q^2}{2m}\psi.
```

For commuting $`\hat E_q`$ and $`\hat{\boldsymbol{P}}_q^2`$, the positive-energy expansion is:

```math
\hat E_q\psi=\left(\frac{\hat{\boldsymbol{P}}_q^2}{2m}
-\frac{\hat{\boldsymbol{P}}_q^4}{8m^3}+\cdots\right)\psi.
```

For weak, slowly varying fields, [Schrödinger’s equation](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation) is obtained by neglecting $`\hat E_q^2\psi`$ relative to $`2m\hat E_q\psi`$:

```math
\boxed{\mathrm i\partial_t\psi=
\left(\frac{(-\mathrm i\partial_{\boldsymbol{r}}-q\boldsymbol{A})^2}{2m}+qV\right)\psi.}
```

With $`\boldsymbol{A}=0`$ (hence $`\boldsymbol{B}=0`$):

```math
\mathrm i\partial_t\psi=
\left(-\frac{\partial_{\boldsymbol{r}}^2}{2m}+qV\right)\psi.
```

## Low-energy [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) to [Pauli](https://en.wikipedia.org/wiki/Pauli_equation)

With the rest-energy phase removed by $`\psi _\pm':=e^{\mathrm imt}\psi _\pm`$:

```math
(m+\hat E_q\pm\hat{\boldsymbol{P}}_q\cdot\boldsymbol{\sigma})\psi _\pm'=m\psi _\mp'.
```

The positive-energy large and small components are defined by:

```math
\begin{aligned}
\phi _+&:=\frac{\psi _+'+\psi _-'}{\sqrt2},\\
\phi _-&:=\frac{\psi _-'-\psi _+'}{\sqrt2}.
\end{aligned}
```

By addition and subtraction:

```math
\begin{aligned}
\hat E_q\phi _+&=(\hat{\boldsymbol{P}}_q\cdot\boldsymbol{\sigma})\phi _-,\\
(2m+\hat E_q)\phi _-&=(\hat{\boldsymbol{P}}_q\cdot\boldsymbol{\sigma})\phi _+.
\end{aligned}
```

In the nonrelativistic limit, $`\hat E_q\phi _-`$ is neglected relative to $`2m\phi _-`$. By substitution and the field identity:

```math
\begin{aligned}
\phi _-&\simeq\frac{\hat{\boldsymbol{P}}_q\cdot\boldsymbol{\sigma}}{2m}\phi _+,\\
\hat E_q\phi _+&\simeq\frac{(\hat{\boldsymbol{P}}_q\cdot\boldsymbol{\sigma})^2}{2m}\phi _+
=\frac{\hat{\boldsymbol{P}}_q^2-q\boldsymbol{B}\cdot\boldsymbol{\sigma}}{2m}\phi _+.
\end{aligned}
```

With $`\hat E_q=\mathrm i\partial_t-qV`$, the [Pauli equation](https://en.wikipedia.org/wiki/Pauli_equation) is obtained through order $`1/m`$:

```math
\boxed{\mathrm i\partial_t\phi _+=
\left(\frac{\hat{\boldsymbol{P}}_q^2}{2m}
+qV-\frac{q}{2m}\boldsymbol{B}\cdot\boldsymbol{\sigma}\right)\phi _+.}
```

## [Spin up and spin down](https://en.wikipedia.org/wiki/Spin-1/2#Observables) in a uniform magnetic field

For constant $`\boldsymbol{B}\ne0`$, the unit direction and [spin projectors](https://en.wikipedia.org/wiki/Pauli_matrices#Eigenvectors_and_eigenvalues) are:

```math
\boldsymbol{u}:=\frac{\boldsymbol{B}}{\lVert\boldsymbol{B}\rVert},
\qquad \Pi _\pm:=\frac12(1\pm\boldsymbol{u}\cdot\boldsymbol{\sigma}).
```

By [spectral decomposition](#projections-and-spectral-decomposition):

```math
\boldsymbol{B}\cdot\boldsymbol{\sigma}
=\lVert\boldsymbol{B}\rVert(\Pi _+-\Pi _-),
\qquad
\Pi _\pm(\boldsymbol{B}\cdot\boldsymbol{\sigma})
=\pm\lVert\boldsymbol{B}\rVert\Pi _\pm.
```

Spin amplitudes along $`\boldsymbol{u}`$ are encoded by:

```math
\begin{aligned}
\phi _+&=\phi_{++}\Pi _++\phi_{+-}\Pi _-,\\
\Pi _\pm\phi _+&=\phi _+\Pi _\pm=\phi_{+\pm}\Pi _\pm.
\end{aligned}
```

The [Pauli equation](https://en.wikipedia.org/wiki/Pauli_equation) is projected using constant $`\Pi _\pm`$, which commute with $`\hat E_q`$ and $`\hat{\boldsymbol{P}}_q^2`$:

```math
\begin{aligned}
0&=\Pi _\pm\left(2m\hat E_q-\hat{\boldsymbol{P}}_q^2
+q\boldsymbol{B}\cdot\boldsymbol{\sigma}\right)\phi _+\\
&=\left(\left(2m\hat E_q-\hat{\boldsymbol{P}}_q^2
\pm q\lVert\boldsymbol{B}\rVert\right)\phi_{+\pm}\right)\Pi _\pm.
\end{aligned}
```

With $`\Pi _\pm\ne0`$ and $`\hat E_q=\mathrm i\partial_t-qV`$:

```math
\boxed{\mathrm i\partial_t\phi_{+\pm}=
\left(\frac{\hat{\boldsymbol{P}}_q^2\mp q\lVert\boldsymbol{B}\rVert}{2m}+qV\right)\phi_{+\pm}.}
```

Opposite [Zeeman shifts](https://en.wikipedia.org/wiki/Zeeman_effect) are obtained for [spin](https://en.wikipedia.org/wiki/Spin-1/2) $`\pm\tfrac12`$ along $`\boldsymbol{u}`$; orbital coupling is retained in [kinetic momentum](https://en.wikipedia.org/wiki/Momentum_operator#Definition_(position_space)) $`\hat{\boldsymbol{P}}_q=-\mathrm i\partial_{\boldsymbol{r}}-q\boldsymbol{A}`$.

For [Dirac dynamics](https://en.wikipedia.org/wiki/Dirac_equation) and measurements, the fixed spinor space is retained.

With the rest-energy phase removed, the complex [Dirac amplitudes](https://en.wikipedia.org/wiki/Dirac_spinor) are arranged as:

```math
\boldsymbol{\phi}:=\begin{Bmatrix}
\phi_{++}&\phi_{+-}\\
\phi_{-+}&\phi_{--}
\end{Bmatrix}.
```

Rows are [Dirac blocks](https://en.wikipedia.org/wiki/Dirac_spinor); columns are spin projections. The upper row is governed by the [Pauli limit](https://en.wikipedia.org/wiki/Pauli_equation).

## [Spin measurement](https://en.wikipedia.org/wiki/Spin-1/2#Rotations_and_Spinors) probabilities

Real unit preparation and measurement axes are denoted by $`\boldsymbol{u}',\boldsymbol{u}`$, with $`\cos\theta:=\boldsymbol{u}\cdot\boldsymbol{u}'`$.

The state $`\psi =\Pi _+'\psi \ne0`$ is prepared using projectors with the following overlap:

```math
\begin{gathered}
\Pi _+':=\frac12(1+\boldsymbol{u}'\cdot\boldsymbol{\sigma}),\qquad
\Pi _\pm:=\frac12(1\pm\boldsymbol{u}\cdot\boldsymbol{\sigma}),\\[4pt]
\Pi _+'\Pi _\pm\Pi _+'=\frac{1\pm\cos\theta}{2}\Pi _+'.
\end{gathered}
```

Using the squared norm, $`\lVert {Z}\rVert^2=\mathrm{sc}({Z}{Z}^{\mathsf H})`$, the denominator is:

```math
\lVert\psi \rVert^2=\mathrm{sc}(\psi \psi ^{\mathsf H})>0.
```

For the [prepared state](https://en.wikipedia.org/wiki/Spin-1/2#Bloch_Representation), the projected squared norm is:

```math
\begin{aligned}
\lVert\Pi _\pm\psi \rVert^2
&=\mathrm{sc}\bigl((\Pi _\pm\psi )(\Pi _\pm\psi )^{\mathsf H}\bigr)\\
&=\mathrm{sc}(\Pi _\pm\psi \psi ^{\mathsf H})\\
&=\mathrm{sc}(\Pi _+'\Pi _\pm\Pi _+'\psi \psi ^{\mathsf H})\\
&=\frac{1\pm\cos\theta}{2}\mathrm{sc}(\psi \psi ^{\mathsf H})
=\frac{1\pm\cos\theta}{2}\lVert\psi \rVert^2.
\end{aligned}
```

The [Born probabilities](https://en.wikipedia.org/wiki/Born_rule) for [spin projections](https://en.wikipedia.org/wiki/Spin-1/2) $`\pm\tfrac12`$ along $`\boldsymbol{u}`$ are:

```math
\begin{gathered}
\mathrm{prob}(\pm)=\frac{\lVert\Pi _\pm\psi \rVert^2}{\lVert\psi \rVert^2}
=\frac{1\pm\cos\theta}{2},\\[6pt]
\boxed{\mathrm{prob}(+)=\cos^2\frac\theta2,\quad \mathrm{prob}(-)=\sin^2\frac\theta2.}
\end{gathered}
```
