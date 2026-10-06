# Space–Time Algebra

**A common [paravector](https://en.wikipedia.org/wiki/Paravector) language is provided by [complex scalars](https://en.wikipedia.org/wiki/Complex_number) and [three-vectors](https://en.wikipedia.org/wiki/Euclidean_vector) for [quaternions](https://en.wikipedia.org/wiki/Quaternion), [rotations](https://en.wikipedia.org/wiki/Rotation_%28mathematics%29), [Lorentz boosts](https://en.wikipedia.org/wiki/Lorentz_transformation), [projections](https://en.wikipedia.org/wiki/Paravector#Null_paravectors_as_projectors), [spacetime geometry](https://en.wikipedia.org/wiki/Minkowski_space), [relativistic particle dynamics](https://en.wikipedia.org/wiki/Relativistic_mechanics), [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations), [electromagnetic waves](https://en.wikipedia.org/wiki/Electromagnetic_radiation) and [forces](https://en.wikipedia.org/wiki/Lorentz_force), [gauge coupling](https://en.wikipedia.org/wiki/Minimal_coupling), the [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) and [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) equations, their [Schrödinger](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation) and [Pauli](https://en.wikipedia.org/wiki/Pauli_equation) low-energy limits, and [spin-one-half](https://en.wikipedia.org/wiki/Spin-1/2) [amplitudes](https://en.wikipedia.org/wiki/Probability_amplitude) and [measurement probabilities](https://en.wikipedia.org/wiki/Born_rule).**

[Paper (PDF)](./sta_notes.pdf) · [LaTeX source](./sta_notes.tex)

[Natural units](https://en.wikipedia.org/wiki/Natural_units) $`\hbar=c=1`$, [rationalized electromagnetic units](https://en.wikipedia.org/wiki/Heaviside%E2%80%93Lorentz_units), and [signature](https://en.wikipedia.org/wiki/Metric_signature) $`(+,-,-,-)`$ are used.

Vectors and matrices are bold cyan; paravectors and paravector operators are magenta, with bold calligraphic Latin letters. Greek glyphs are retained; letter case is unrestricted.

## One complex [paravector](https://en.wikipedia.org/wiki/Paravector)

Complex paravectors are formed from sigmas—algebraic objects with order-dependent multiplication:

```math
\boxed{
\begin{aligned}
\sigma_n^2&=1,
&\sigma_n\sigma_m&=-\sigma_m\sigma_n\quad(n\ne m),\\
\mathrm i&:=\sigma_1\sigma_2\sigma_3,
&\mathrm i^2&=-1,\\
{\color{cyan}𝝈}&:=\{\sigma_1,\sigma_2,\sigma_3\},
&\mathrm i{\color{cyan}𝝈}&=\{\sigma_2\sigma_3,\sigma_3\sigma_1,\sigma_1\sigma_2\},\\
{\color{magenta}\boldsymbol{\mathcal Z}}&:=(a+b\mathrm i)+({\color{cyan}𝐗}+\mathrm i{\color{cyan}𝐘})\cdot{\color{cyan}𝝈},
&a,b&\in\mathbb R,\quad{\color{cyan}𝐗},{\color{cyan}𝐘}\in\mathbb R^3.
\end{aligned}}
```

The eight real components are grouped into a complex scalar and vector:

```math
\begin{aligned}
{\color{magenta}\boldsymbol{\mathcal X}}:&=a+{\color{cyan}𝐗}\cdot{\color{cyan}𝝈},\\
{\color{magenta}\boldsymbol{\mathcal Y}}:&=b+{\color{cyan}𝐘}\cdot{\color{cyan}𝝈},\\[4pt]
{\color{magenta}\boldsymbol{\mathcal Z}}&={\color{magenta}\boldsymbol{\mathcal X}}+\mathrm i{\color{magenta}\boldsymbol{\mathcal Y}}\\
&=(a+b\mathrm i)+({\color{cyan}𝐗}+\mathrm i{\color{cyan}𝐘})\cdot{\color{cyan}𝝈}\\
&=c+{\color{cyan}𝐙}\cdot{\color{cyan}𝝈}.
\end{aligned}
```

The component extractions are:

```math
\begin{aligned}
\mathrm{re}({\color{magenta}\boldsymbol{\mathcal Z}})&:={\color{magenta}\boldsymbol{\mathcal X}},\\
\mathrm{im}({\color{magenta}\boldsymbol{\mathcal Z}})&:={\color{magenta}\boldsymbol{\mathcal Y}},\\
\mathrm{sc}({\color{magenta}\boldsymbol{\mathcal Z}})&:=c=a+\mathrm ib,\\
\mathrm{vec}({\color{magenta}\boldsymbol{\mathcal Z}})&:={\color{cyan}𝐙}={\color{cyan}𝐗}+\mathrm i{\color{cyan}𝐘}.
\end{aligned}
```

## Product

Both the [dot](https://en.wikipedia.org/wiki/Dot_product) and [cross](https://en.wikipedia.org/wiki/Cross_product) products are included in the product:

```math
({\color{cyan}𝐗}\cdot{\color{cyan}𝝈})({\color{cyan}𝐘}\cdot{\color{cyan}𝝈})
={\color{cyan}𝐗}\cdot{\color{cyan}𝐘}+\mathrm i({\color{cyan}𝐗}\times{\color{cyan}𝐘})\cdot{\color{cyan}𝝈}.
```

The convention $`{\color{cyan}𝐗}^2:={\color{cyan}𝐗}\cdot{\color{cyan}𝐗}`$ is used.

The [commutator](https://en.wikipedia.org/wiki/Commutator) is:

```math
[{\color{magenta}\boldsymbol{\mathcal X}},{\color{magenta}\boldsymbol{\mathcal Y}}]:={\color{magenta}\boldsymbol{\mathcal X}}{\color{magenta}\boldsymbol{\mathcal Y}}-{\color{magenta}\boldsymbol{\mathcal Y}}{\color{magenta}\boldsymbol{\mathcal X}}.
```

```math
[{\color{cyan}𝐗}\cdot{\color{cyan}𝝈},{\color{cyan}𝐘}\cdot{\color{cyan}𝝈}]=2\mathrm i({\color{cyan}𝐗}\times{\color{cyan}𝐘})\cdot{\color{cyan}𝝈}.
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

For $`{\color{magenta}\boldsymbol{\mathcal Z}}=c+{\color{cyan}𝐙}\cdot{\color{cyan}𝝈}`$, $`c\in\mathbb C`$, $`{\color{cyan}𝐙}\in\mathbb C^3`$, the complex coefficients are conjugated:

```math
\begin{aligned}
{\color{magenta}\boldsymbol{\mathcal Z}}^{\mathsf H}&=c^{*}+{\color{cyan}𝐙}^{*}\cdot{\color{cyan}𝝈},\\[4pt]
{\color{magenta}\boldsymbol{\mathcal Z}}^{*}&=c^{*}-{\color{cyan}𝐙}^{*}\cdot{\color{cyan}𝝈}.
\end{aligned}
```

## Core operations

The paravector [adjugate](https://en.wikipedia.org/wiki/Adjugate_matrix), [determinant](https://en.wikipedia.org/wiki/Determinant), and [inverse](https://en.wikipedia.org/wiki/Invertible_matrix) are given by:

```math
\begin{aligned}
\mathrm{adj}\,{\color{magenta}\boldsymbol{\mathcal Z}}&:=c-{\color{cyan}𝐙}\cdot{\color{cyan}𝝈}={\color{magenta}\boldsymbol{\mathcal Z}}^{*\mathsf H},\\[4pt]
\det {\color{magenta}\boldsymbol{\mathcal Z}}&:=c^2-{\color{cyan}𝐙}^2={\color{magenta}\boldsymbol{\mathcal Z}}\,\mathrm{adj}\,{\color{magenta}\boldsymbol{\mathcal Z}},\\[6pt]
{\color{magenta}\boldsymbol{\mathcal Z}}^{-1}&=\frac{\mathrm{adj}\,{\color{magenta}\boldsymbol{\mathcal Z}}}{\det {\color{magenta}\boldsymbol{\mathcal Z}}}\qquad\text{if }\det {\color{magenta}\boldsymbol{\mathcal Z}}\ne0.
\end{aligned}
```

The paravector [trace](https://en.wikipedia.org/wiki/Trace_%28linear_algebra%29) and [squared norm](https://en.wikipedia.org/wiki/Norm_%28mathematics%29) are given by:

```math
\begin{aligned}
\mathrm{tr}\,{\color{magenta}\boldsymbol{\mathcal Z}}&:=2c={\color{magenta}\boldsymbol{\mathcal Z}}+\mathrm{adj}\,{\color{magenta}\boldsymbol{\mathcal Z}},\\[4pt]
\lVert {\color{magenta}\boldsymbol{\mathcal Z}}\rVert^2&:=cc^{*}+{\color{cyan}𝐙}\cdot{\color{cyan}𝐙}^{*}=\mathrm{sc}({\color{magenta}\boldsymbol{\mathcal Z}}{\color{magenta}\boldsymbol{\mathcal Z}}^{\mathsf H}).
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

The sign is reversed when two distinct factors are exchanged. For $`a\in\mathbb R`$ and $`{\color{cyan}𝐐}\in\mathbb R^3`$, the quaternion and its conjugate are:

```math
\begin{aligned}
{\color{magenta}\boldsymbol{\mathcal Q}}&:=a+{\color{cyan}𝐐}\cdot(-\mathrm i{\color{cyan}𝝈})=a-\mathrm i{\color{cyan}𝐐}\cdot{\color{cyan}𝝈},\\
{\color{magenta}\boldsymbol{\mathcal Q}}^{\mathsf H}&:=a-{\color{cyan}𝐐}\cdot(-\mathrm i{\color{cyan}𝝈})=a+\mathrm i{\color{cyan}𝐐}\cdot{\color{cyan}𝝈},
\end{aligned}
```

with

```math
{\color{magenta}\boldsymbol{\mathcal Q}}{\color{magenta}\boldsymbol{\mathcal Q}}^{\mathsf H}=\det {\color{magenta}\boldsymbol{\mathcal Q}}=a^2+{\color{cyan}𝐐}^2.
```

Nonzero quaternions are normalized with:

```math
{\color{magenta}\boldsymbol{\mathcal Q}}\leftarrow\frac{{\color{magenta}\boldsymbol{\mathcal Q}}}{\sqrt{{\color{magenta}\boldsymbol{\mathcal Q}}{\color{magenta}\boldsymbol{\mathcal Q}}^{\mathsf H}}},
```

Thus, with $`\lVert{\color{cyan}𝐐}\rVert^2={\color{cyan}𝐐}^2`$ for real $`{\color{cyan}𝐐}`$, the normalization is expressed as:

```math
{\color{magenta}\boldsymbol{\mathcal Q}}{\color{magenta}\boldsymbol{\mathcal Q}}^{\mathsf H}=a^2+\lVert{\color{cyan}𝐐}\rVert^2=1.
```

For a unit quaternion with $`{\color{cyan}𝐐}\ne0`$, the rotation axis and angle are specified by:

```math
{\color{cyan}𝐮}=\frac{{\color{cyan}𝐐}}{\lVert{\color{cyan}𝐐}\rVert},\qquad
a=\cos\frac{\theta}{2},\qquad
{\color{cyan}𝐐}={\color{cyan}𝐮}\sin\frac{\theta}{2}.
```

A rotation of $`{\color{cyan}𝐫}\in\mathbb R^3`$ about $`{\color{cyan}𝐮}`$ by angle $`\theta`$ is expressed by:

```math
{\color{magenta}\boldsymbol{\mathcal X}}'={\color{cyan}𝐫}'\cdot(-\mathrm i{\color{cyan}𝝈}), \qquad
{\color{magenta}\boldsymbol{\mathcal X}}={\color{cyan}𝐫}\cdot(-\mathrm i{\color{cyan}𝝈}), \qquad
{\color{magenta}\boldsymbol{\mathcal X}}'={\color{magenta}\boldsymbol{\mathcal Q}}{\color{magenta}\boldsymbol{\mathcal X}}{\color{magenta}\boldsymbol{\mathcal Q}}^{\mathsf H}.
```

## Rotations and boosts

For real spacetime $`{\color{magenta}\boldsymbol{\mathcal X}}`$ and a real [unit axis](https://en.wikipedia.org/wiki/Unit_vector) $`{\color{cyan}𝐮}`$, the transformation is written as:

```math
{\color{magenta}\boldsymbol{\mathcal X}}:=t+{\color{cyan}𝐫}\cdot{\color{cyan}𝝈},
\qquad {\color{magenta}\boldsymbol{\mathcal X}}':={\color{magenta}\boldsymbol{\mathcal T}}{\color{magenta}\boldsymbol{\mathcal X}}{\color{magenta}\boldsymbol{\mathcal T}}^{\mathsf H}=t'+{\color{cyan}𝐫}'\cdot{\color{cyan}𝝈}.
```

For $`{\color{magenta}\boldsymbol{\mathcal T}}={\color{magenta}\boldsymbol{\mathcal R}}`$, a [Rodrigues rotation](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula) through $`\theta`$, with [unit quaternion](https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation) $`{\color{magenta}\boldsymbol{\mathcal Q}}={\color{magenta}\boldsymbol{\mathcal R}}`$, is given by:

```math
\begin{aligned}
{\color{magenta}\boldsymbol{\mathcal R}}&:=e^{-\mathrm i\theta{\color{cyan}𝐮}\cdot{\color{cyan}𝝈}/2}
=\cos\frac\theta2-\mathrm i{\color{cyan}𝐮}\cdot{\color{cyan}𝝈}\sin\frac\theta2,\\
{\color{magenta}\boldsymbol{\mathcal X}}'&={\color{magenta}\boldsymbol{\mathcal R}}{\color{magenta}\boldsymbol{\mathcal X}}{\color{magenta}\boldsymbol{\mathcal R}}^{\mathsf H},\\
t'&=t,\\
{\color{cyan}𝐫}'&={\color{cyan}𝐫}\cos\theta+({\color{cyan}𝐮}\times{\color{cyan}𝐫})\sin\theta
+({\color{cyan}𝐮}\cdot{\color{cyan}𝐫})(1-\cos\theta){\color{cyan}𝐮}.
\end{aligned}
```

For $`{\color{magenta}\boldsymbol{\mathcal T}}={\color{magenta}\boldsymbol{\mathcal L}}`$, a [Lorentz boost](https://en.wikipedia.org/wiki/Lorentz_transformation) to a frame moving at $`+\beta{\color{cyan}𝐮}`$, with [rapidity](https://en.wikipedia.org/wiki/Rapidity) $`\theta`$, $`\beta:=\tanh\theta`$, and [Lorentz factor](https://en.wikipedia.org/wiki/Lorentz_factor) $`\gamma:=\cosh\theta`$ is considered.

The parallel and perpendicular components are defined by:

```math
{\color{cyan}𝐫}_\parallel:=({\color{cyan}𝐮}\cdot{\color{cyan}𝐫}){\color{cyan}𝐮},
\qquad {\color{cyan}𝐫}_\perp:={\color{cyan}𝐫}-{\color{cyan}𝐫}_\parallel.
```

The boost and transformed coordinates are:

```math
\begin{aligned}
{\color{magenta}\boldsymbol{\mathcal L}}&:=e^{-\theta{\color{cyan}𝐮}\cdot{\color{cyan}𝝈}/2}
=\cosh\frac\theta2-{\color{cyan}𝐮}\cdot{\color{cyan}𝝈}\sinh\frac\theta2
={\color{magenta}\boldsymbol{\mathcal L}}^{\mathsf H},\\
{\color{magenta}\boldsymbol{\mathcal X}}'&={\color{magenta}\boldsymbol{\mathcal L}}{\color{magenta}\boldsymbol{\mathcal X}}{\color{magenta}\boldsymbol{\mathcal L}}^{\mathsf H},\\
t'&=\gamma(t-\beta{\color{cyan}𝐮}\cdot{\color{cyan}𝐫}),\\
{\color{cyan}𝐫}'&={\color{cyan}𝐫}_\perp+\gamma({\color{cyan}𝐫}_\parallel-\beta t{\color{cyan}𝐮}).
\end{aligned}
```

## Projections and spectral decomposition

A paravector of the following form is considered:

```math
{\color{magenta}\boldsymbol{\mathcal Z}}=c+z{\color{cyan}𝐮}\cdot{\color{cyan}𝝈},
\qquad c,z\in\mathbb C,\quad {\color{cyan}𝐮}\in\mathbb C^3,\quad {\color{cyan}𝐮}^2=1.
```

The eigenvalues and complementary [projectors](https://en.wikipedia.org/wiki/Paravector#Null_paravectors_as_projectors) are:

```math
\begin{gathered}
\lambda_\pm:=c\pm z,
\qquad \textcolor{magenta}\Pi _\pm:=\frac12(1\pm{\color{cyan}𝐮}\cdot{\color{cyan}𝝈}),\\[4pt]
\textcolor{magenta}\Pi _\pm^2=\textcolor{magenta}\Pi _\pm,\qquad \textcolor{magenta}\Pi _+\textcolor{magenta}\Pi _-=0,\qquad \textcolor{magenta}\Pi _-+\textcolor{magenta}\Pi _+=1.
\end{gathered}
```

For real $`{\color{cyan}𝐮}`$, these are also Hermitian projectors. For $`f`$ analytic near $`\lambda_\pm`$, the following spectral decomposition is obtained:

```math
\begin{aligned}
{\color{magenta}\boldsymbol{\mathcal Z}}&=\lambda_-\textcolor{magenta}\Pi _-+\lambda_+\textcolor{magenta}\Pi _+,\\[4pt]
f({\color{magenta}\boldsymbol{\mathcal Z}})&=f(\lambda_-)\textcolor{magenta}\Pi _-+f(\lambda_+)\textcolor{magenta}\Pi _+.
\end{aligned}
```

## Spacetime from the [determinant](https://en.wikipedia.org/wiki/Determinant)

In natural units, inertial coordinate time is denoted by $`t`$ and onboard proper time (often $`\tau`$) by $`s`$.

The spacetime increment and interval are:

```math
\begin{aligned}
\mathrm d{\color{magenta}\boldsymbol{\mathcal X}}&=\mathrm dt+\mathrm d{\color{cyan}𝐫}\cdot{\color{cyan}𝝈},\\
\det(\mathrm d{\color{magenta}\boldsymbol{\mathcal X}})&=\mathrm dt^2-\mathrm d{\color{cyan}𝐫}^2.
\end{aligned}
```

For a future-directed time-like path ($`\mathrm dt>0`$), positive [proper time](https://en.wikipedia.org/wiki/Proper_time) is defined by:

```math
\mathrm ds:=\sqrt{\mathrm dt^2-\mathrm d{\color{cyan}𝐫}^2}=\sqrt{\det(\mathrm d{\color{magenta}\boldsymbol{\mathcal X}})}.
```

Dividing by coordinate time gives:

```math
\begin{aligned}
\frac{\mathrm ds}{\mathrm dt}
&=\sqrt{\frac{\mathrm dt^2-\mathrm d{\color{cyan}𝐫}^2}{\mathrm dt^2}}\\
&=\sqrt{1-\left(\frac{\mathrm d{\color{cyan}𝐫}}{\mathrm dt}\right)^2}
=\sqrt{1-{\color{cyan}𝐯}^2},
\qquad {\color{cyan}𝐯}:=\frac{\mathrm d{\color{cyan}𝐫}}{\mathrm dt}.
\end{aligned}
```

The [Lorentz factor](https://en.wikipedia.org/wiki/Lorentz_factor) is:

```math
\gamma:=\frac{\mathrm dt}{\mathrm ds}=\frac{1}{\sqrt{1-{\color{cyan}𝐯}^2}}.
```

The [four-velocity](https://en.wikipedia.org/wiki/Four-velocity) is:

```math
\begin{aligned}
{\color{magenta}\boldsymbol{\mathcal U}}:=\frac{\mathrm d{\color{magenta}\boldsymbol{\mathcal X}}}{\mathrm ds}
&=\frac{\mathrm dt+\mathrm d{\color{cyan}𝐫}\cdot{\color{cyan}𝝈}}{\mathrm ds}\\
&=\left(1+\frac{\mathrm d{\color{cyan}𝐫}}{\mathrm dt}\cdot{\color{cyan}𝝈}\right)
\frac{\mathrm dt}{\mathrm ds}\\
&=\gamma(1+{\color{cyan}𝐯}\cdot{\color{cyan}𝝈}).
\end{aligned}
```

Thus:

```math
\det {\color{magenta}\boldsymbol{\mathcal U}}=\det\!\left(\frac{\mathrm d{\color{magenta}\boldsymbol{\mathcal X}}}{\mathrm ds}\right)
=\frac{\det(\mathrm d{\color{magenta}\boldsymbol{\mathcal X}})}{\mathrm ds^2}=1.
```

For energy $`E`$ and momentum $`{\color{cyan}𝐏}`$, the [four-momentum](https://en.wikipedia.org/wiki/Four-momentum) is:

```math
{\color{magenta}\boldsymbol{\mathcal P}}:=m{\color{magenta}\boldsymbol{\mathcal U}}=E+{\color{cyan}𝐏}\cdot{\color{cyan}𝝈}.
```

For constant mass $`m>0`$:

```math
\det {\color{magenta}\boldsymbol{\mathcal P}}=\det(m{\color{magenta}\boldsymbol{\mathcal U}})=m^2\det {\color{magenta}\boldsymbol{\mathcal U}}=m^2.
```

The [mass-shell](https://en.wikipedia.org/wiki/On_shell_and_off_shell#Mass_shell) relation follows:

```math
\det {\color{magenta}\boldsymbol{\mathcal P}}=E^2-{\color{cyan}𝐏}^2=m^2.
```

## [Maxwell](https://en.wikipedia.org/wiki/Maxwell%27s_equations) in one equation

The [electric field](https://en.wikipedia.org/wiki/Electric_field) $`{\color{cyan}𝐄}`$ and [magnetic field](https://en.wikipedia.org/wiki/Magnetic_field) $`{\color{cyan}𝐁}`$ are combined into the complex electromagnetic field $`{\color{cyan}𝐅}`$ and its paravector $`{\color{magenta}\boldsymbol{\mathcal F}}`$:

```math
{\color{cyan}𝐅}:={\color{cyan}𝐄}+\mathrm i{\color{cyan}𝐁},\qquad {\color{magenta}\boldsymbol{\mathcal F}}:={\color{cyan}𝐅}\cdot{\color{cyan}𝝈}.
```

The [charge density](https://en.wikipedia.org/wiki/Charge_density) $`\rho`$ and [current density](https://en.wikipedia.org/wiki/Current_density) $`{\color{cyan}𝐉}`$ are combined into the [four-current](https://en.wikipedia.org/wiki/Four-current) $`{\color{magenta}\boldsymbol{\mathcal J}}`$:

```math
{\color{magenta}\boldsymbol{\mathcal J}}:=\rho+{\color{cyan}𝐉}\cdot{\color{cyan}𝝈},\qquad {\color{magenta}\boldsymbol{\mathcal J}}^*=\rho-{\color{cyan}𝐉}\cdot{\color{cyan}𝝈}.
```

The [paravector derivative](https://en.wikipedia.org/wiki/Paravector#Paragradient), with [gradient](https://en.wikipedia.org/wiki/Gradient) $`\partial_{{\color{cyan}𝐫}}`$, is defined by:

```math
\textcolor{magenta}\partial :=\partial_t+\partial_{{\color{cyan}𝐫}}\cdot{\color{cyan}𝝈},\qquad
\textcolor{magenta}\partial ^*=\partial_t-\partial_{{\color{cyan}𝐫}}\cdot{\color{cyan}𝝈}.
```

[Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations) are given by:

```math
\boxed{\textcolor{magenta}\partial  {\color{magenta}\boldsymbol{\mathcal F}}={\color{magenta}\boldsymbol{\mathcal J}}^*},
```

with components:

```math
\begin{aligned}
\text{Real scalar:}\quad &\partial_{{\color{cyan}𝐫}}\cdot{\color{cyan}𝐄}=\rho,\\
\text{Imaginary scalar:}\quad &\partial_{{\color{cyan}𝐫}}\cdot{\color{cyan}𝐁}=0,\\
\text{Real vector:}\quad &\partial_t{\color{cyan}𝐄}-\partial_{{\color{cyan}𝐫}}\times{\color{cyan}𝐁}=-{\color{cyan}𝐉},\\
\text{Imaginary vector:}\quad &\partial_t{\color{cyan}𝐁}+\partial_{{\color{cyan}𝐫}}\times{\color{cyan}𝐄}=0.
\end{aligned}
```

By application of $`\textcolor{magenta}\partial ^{*}`$, the [d'Alembertian](https://en.wikipedia.org/wiki/D%27Alembert_operator) and the [sourced wave equations](https://en.wikipedia.org/wiki/Inhomogeneous_electromagnetic_wave_equation) are obtained:

```math
\begin{aligned}
\Box&:=\textcolor{magenta}\partial ^{*}\textcolor{magenta}\partial =\partial_t^2-\partial_{{\color{cyan}𝐫}}^2,\\
\textcolor{magenta}\partial ^{*}(\textcolor{magenta}\partial  {\color{magenta}\boldsymbol{\mathcal F}})&=\Box {\color{magenta}\boldsymbol{\mathcal F}}=\textcolor{magenta}\partial ^{*}{\color{magenta}\boldsymbol{\mathcal J}}^*.
\end{aligned}
```

Equivalently:

```math
\boxed{\Box {\color{magenta}\boldsymbol{\mathcal F}}^*=\textcolor{magenta}\partial  {\color{magenta}\boldsymbol{\mathcal J}}}.
```

From this boxed equation, [charge conservation](https://en.wikipedia.org/wiki/Charge_conservation) and the field wave equations are obtained component by component:

```math
\begin{aligned}
\text{Real scalar:}\quad &0=\partial_t\rho+\partial_{{\color{cyan}𝐫}}\cdot{\color{cyan}𝐉}
&&\text{(charge conservation)},\\
\text{Imaginary scalar:}\quad &0=0,\\
\text{Real vector:}\quad &\Box{\color{cyan}𝐄}=-\partial_t{\color{cyan}𝐉}-\partial_{{\color{cyan}𝐫}}\rho,\\
\text{Imaginary vector:}\quad &\Box{\color{cyan}𝐁}=\partial_{{\color{cyan}𝐫}}\times{\color{cyan}𝐉}.
\end{aligned}
```

In source-free vacuum, $`{\color{magenta}\boldsymbol{\mathcal J}}=0`$ and $`\Box {\color{magenta}\boldsymbol{\mathcal F}}=0`$.

## [Electromagnetic potentials](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) and [gauge invariance](https://en.wikipedia.org/wiki/Electromagnetic_four-potential#Gauge_freedom)

The real [four-potential](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) is defined from the [scalar potential](https://en.wikipedia.org/wiki/Electric_potential) $`V`$ and [vector potential](https://en.wikipedia.org/wiki/Magnetic_vector_potential) $`{\color{cyan}𝐀}`$, with four-divergence $`S`$:

```math
{\color{magenta}\boldsymbol{\mathcal A}}:=V+{\color{cyan}𝐀}\cdot{\color{cyan}𝝈},
\qquad S:=\mathrm{sc}(\textcolor{magenta}\partial  {\color{magenta}\boldsymbol{\mathcal A}}).
```

The field is obtained from the potential:

```math
\boxed{{\color{magenta}\boldsymbol{\mathcal F}}^{*}=\textcolor{magenta}\partial  {\color{magenta}\boldsymbol{\mathcal A}}-S}.
```

In components:

```math
\begin{aligned}
\text{Real scalar:}\quad &S=\partial_tV+\partial_{{\color{cyan}𝐫}}\cdot{\color{cyan}𝐀},\\
\text{Imaginary scalar:}\quad &0=0,\\
\text{Real vector:}\quad &{\color{cyan}𝐄}=-\partial_t{\color{cyan}𝐀}-\partial_{{\color{cyan}𝐫}}V,\\
\text{Imaginary vector:}\quad &{\color{cyan}𝐁}=\partial_{{\color{cyan}𝐫}}\times{\color{cyan}𝐀}.
\end{aligned}
```

The source equations are obtained in any [gauge](https://en.wikipedia.org/wiki/Gauge_fixing) by application of $`\textcolor{magenta}\partial ^*`$:

```math
\Box {\color{magenta}\boldsymbol{\mathcal A}}=\textcolor{magenta}\partial ^{*}(\textcolor{magenta}\partial  {\color{magenta}\boldsymbol{\mathcal A}})
=\textcolor{magenta}\partial ^{*}(S+{\color{magenta}\boldsymbol{\mathcal F}}^{*})
=\textcolor{magenta}\partial ^{*}S+(\textcolor{magenta}\partial  {\color{magenta}\boldsymbol{\mathcal F}})^{*}.
```

By substitution of [Maxwell’s equation](https://en.wikipedia.org/wiki/Maxwell%27s_equations):

```math
\boxed{\Box {\color{magenta}\boldsymbol{\mathcal A}}=\textcolor{magenta}\partial ^{*}S+{\color{magenta}\boldsymbol{\mathcal J}}}.
```

In components:

```math
\begin{aligned}
\text{Real scalar:}\quad &\Box V-\partial_tS=\rho,\\
\text{Imaginary scalar:}\quad &0=0,\\
\text{Real vector:}\quad &\Box{\color{cyan}𝐀}+\partial_{{\color{cyan}𝐫}}S={\color{cyan}𝐉},\\
\text{Imaginary vector:}\quad &{\color{cyan}𝟎}={\color{cyan}𝟎}.
\end{aligned}
```

Under a [gauge transformation](https://en.wikipedia.org/wiki/Electromagnetic_four-potential#Gauge_freedom), with real scalar function $`\lambda`$, the field is left unchanged:

```math
\begin{aligned}
{\color{magenta}\boldsymbol{\mathcal A}}'&:={\color{magenta}\boldsymbol{\mathcal A}}+\textcolor{magenta}\partial ^{*}\lambda,\\
S'&=S+\Box\lambda,\\
{\color{magenta}\boldsymbol{\mathcal F}}'&={\color{magenta}\boldsymbol{\mathcal F}}.
\end{aligned}
```

In the [Lorenz gauge](https://en.wikipedia.org/wiki/Lorenz_gauge_condition), the [potential wave equations](https://en.wikipedia.org/wiki/Inhomogeneous_electromagnetic_wave_equation) are obtained with $`S=0`$:

```math
\begin{aligned}
\text{Real scalar:}\quad &\Box V=\rho,\\
\text{Real vector:}\quad &\Box{\color{cyan}𝐀}={\color{cyan}𝐉}.
\end{aligned}
```

## [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) from the [mass shell](https://en.wikipedia.org/wiki/On_shell_and_off_shell#Mass_shell)

Hats are used for named [energy](https://en.wikipedia.org/wiki/Energy_operator), [momentum](https://en.wikipedia.org/wiki/Momentum_operator), and [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) operators; the derivative symbols $`\textcolor{magenta}\partial `$ and $`\Box`$ remain unhatted.

On a complex scalar wavefunction $`\psi`$, the canonical [energy](https://en.wikipedia.org/wiki/Energy_operator) and [momentum](https://en.wikipedia.org/wiki/Momentum_operator) operators are defined by:

```math
\begin{aligned}
\hat E&:=\mathrm i\partial_t,\\
\hat{{\color{cyan}𝐏}}&:=-\mathrm i\partial_{{\color{cyan}𝐫}},
\end{aligned}
```

The canonical four-momentum operator is formed by:

```math
\hat{{\color{magenta}\boldsymbol{\mathcal P}}}:=\mathrm i\textcolor{magenta}\partial ^*=\hat E+\hat{{\color{cyan}𝐏}}\cdot{\color{cyan}𝝈}.
```

By commutativity of the free operators:

```math
\hat{{\color{magenta}\boldsymbol{\mathcal P}}}\,\mathrm{adj}\,\hat{{\color{magenta}\boldsymbol{\mathcal P}}}
=\det\hat{{\color{magenta}\boldsymbol{\mathcal P}}}
=\hat E^2-\hat{{\color{cyan}𝐏}}^2
=-\Box.
```

With the [mass shell](https://en.wikipedia.org/wiki/On_shell_and_off_shell#Mass_shell) imposed on $`\psi`$:

```math
\begin{aligned}
0&=(\det\hat{{\color{magenta}\boldsymbol{\mathcal P}}}-m^2)\psi\\
&=(\hat E^2-\hat{{\color{cyan}𝐏}}^2-m^2)\psi\\
&=-(\Box+m^2)\psi.
\end{aligned}
```

Thus the free [Klein–Gordon equation](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) is:

```math
\boxed{(\Box+m^2)\psi=0.}
```

## [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) as a first-order wave equation

In the [Weyl representation](https://en.wikipedia.org/wiki/Gamma_matrices#Weyl_%28chiral%29_basis), the notation $`{\color{cyan}𝝍}:=(\textcolor{magenta}\psi _+,\textcolor{magenta}\psi _-)^{\mathsf T}`$ is used, with two complex components in each entry. The block identity is denoted by $`{\color{cyan}𝐈}`$.

The following operators are defined using the free momentum operator above:

```math
\begin{aligned}
{\color{cyan}𝐖}(\hat{{\color{magenta}\boldsymbol{\mathcal P}}})&:=
\begin{pmatrix}
0 & \mathrm{adj}\,\hat{{\color{magenta}\boldsymbol{\mathcal P}}} \\
\hat{{\color{magenta}\boldsymbol{\mathcal P}}} & 0
\end{pmatrix},\\
\hat{{\color{cyan}𝐃}}(m)&:={\color{cyan}𝐖}(\hat{{\color{magenta}\boldsymbol{\mathcal P}}})-m{\color{cyan}𝐈}.
\end{aligned}
```

The [Dirac equation](https://en.wikipedia.org/wiki/Dirac_equation) is:

```math
\boxed{\hat{{\color{cyan}𝐃}}(m){\color{cyan}𝝍}=0.}
```

Equivalently:

```math
(\hat E\pm\hat{{\color{cyan}𝐏}}\cdot{\color{cyan}𝝈})\textcolor{magenta}\psi _\pm=m\textcolor{magenta}\psi _\mp.
```

By $`{\color{cyan}𝐖}(\hat{{\color{magenta}\boldsymbol{\mathcal P}}})^2=(\det\hat{{\color{magenta}\boldsymbol{\mathcal P}}}){\color{cyan}𝐈}`$, the [Klein–Gordon factorization](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) is obtained:

```math
\begin{aligned}
\hat{{\color{cyan}𝐃}}(m)\hat{{\color{cyan}𝐃}}(-m)&=\hat{{\color{cyan}𝐃}}(-m)\hat{{\color{cyan}𝐃}}(m)\\
&={\color{cyan}𝐖}(\hat{{\color{magenta}\boldsymbol{\mathcal P}}})^2-m^2{\color{cyan}𝐈}
=-(\Box+m^2){\color{cyan}𝐈}.
\end{aligned}
```

Hence [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) is satisfied by every free [Dirac solution](https://en.wikipedia.org/wiki/Dirac_equation).

## [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) with an electromagnetic potential

Here, the uncoupled operators defined above are labeled by subscript $`0`$:

```math
\begin{aligned}
\hat E_0&:=\mathrm i\partial_t,\\
\hat{{\color{cyan}𝐏}}_0&:=-\mathrm i\partial_{{\color{cyan}𝐫}},\\
\hat{{\color{magenta}\boldsymbol{\mathcal P}}}_0&:=\mathrm i\textcolor{magenta}\partial ^*.
\end{aligned}
```

The following are assumed: constant mass $`m>0`$, [charge](https://en.wikipedia.org/wiki/Electric_charge) $`q`$, and a real [potential](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) $`{\color{magenta}\boldsymbol{\mathcal A}}=V+{\color{cyan}𝐀}\cdot{\color{cyan}𝝈}`$.

By [minimal coupling](https://en.wikipedia.org/wiki/Minimal_coupling), the [gauge-covariant energy](https://en.wikipedia.org/wiki/Gauge_covariant_derivative) and [kinetic momentum](https://en.wikipedia.org/wiki/Momentum_operator#Definition_(position_space)) operators are obtained:

```math
\begin{aligned}
\hat E_q&:=\hat E_0-qV,\\
\hat{{\color{cyan}𝐏}}_q&:=\hat{{\color{cyan}𝐏}}_0-q{\color{cyan}𝐀}.
\end{aligned}
```

The [kinetic four-momentum](https://en.wikipedia.org/wiki/Minimal_coupling) operator is formed by:

```math
\hat{{\color{magenta}\boldsymbol{\mathcal P}}}_q:=\hat{{\color{magenta}\boldsymbol{\mathcal P}}}_0-q{\color{magenta}\boldsymbol{\mathcal A}}=\hat E_q+\hat{{\color{cyan}𝐏}}_q\cdot{\color{cyan}𝝈}.
```

In $`\hat E_q`$, rest energy is included and the potential energy $`qV`$ is excluded.

For operator arguments, the spatial vector part is reversed by $`\mathrm{adj}`$ while derivative order is preserved.

The coupled [Dirac operator](https://en.wikipedia.org/wiki/Dirac_equation) is:

```math
\begin{aligned}
\hat{{\color{cyan}𝐃}}_q(m)&:={\color{cyan}𝐖}(\hat{{\color{magenta}\boldsymbol{\mathcal P}}}_q)-m{\color{cyan}𝐈}\\
&=\begin{pmatrix}
-m & \mathrm{adj}\,\hat{{\color{magenta}\boldsymbol{\mathcal P}}}_q\\
\hat{{\color{magenta}\boldsymbol{\mathcal P}}}_q & -m
\end{pmatrix}.
\end{aligned}
```

The coupled [Dirac equation](https://en.wikipedia.org/wiki/Dirac_equation) is:

```math
\boxed{\hat{{\color{cyan}𝐃}}_q(m){\color{cyan}𝝍}=0.}
```

Equivalently:

```math
(\hat E_q\pm\hat{{\color{cyan}𝐏}}_q\cdot{\color{cyan}𝝈})\textcolor{magenta}\psi _\pm=m\textcolor{magenta}\psi _\mp.
```

The following field identities are satisfied by the coupled operators:

```math
\begin{aligned}
{}[\hat{{\color{cyan}𝐏}}_q\cdot{\color{cyan}𝝈},\hat E_q]
&=-\mathrm iq{\color{cyan}𝐄}\cdot{\color{cyan}𝝈},\\
(\hat{{\color{cyan}𝐏}}_q\cdot{\color{cyan}𝝈})^2
&=\hat{{\color{cyan}𝐏}}_q^2-q{\color{cyan}𝐁}\cdot{\color{cyan}𝝈},\\
\hat{{\color{cyan}𝐏}}_q^2+\mathrm iq{\color{magenta}\boldsymbol{\mathcal F}}
&=(\hat{{\color{cyan}𝐏}}_q\cdot{\color{cyan}𝝈})^2
-[\hat{{\color{cyan}𝐏}}_q\cdot{\color{cyan}𝝈},\hat E_q].
\end{aligned}
```

With derivatives applied to the potentials as well as the wavefunction, the opposite-mass product is obtained:

```math
\hat{{\color{cyan}𝐃}}_q(-m)\hat{{\color{cyan}𝐃}}_q(m)
=(\hat E_q^2-\hat{{\color{cyan}𝐏}}_q^2-m^2){\color{cyan}𝐈}
-\mathrm iq\begin{pmatrix}{\color{magenta}\boldsymbol{\mathcal F}}^{*}&0\\0&{\color{magenta}\boldsymbol{\mathcal F}}\end{pmatrix}.
```

For every [Dirac solution](https://en.wikipedia.org/wiki/Dirac_equation), $`\hat{{\color{cyan}𝐃}}_q(-m)\hat{{\color{cyan}𝐃}}_q(m){\color{cyan}𝝍}=0`$ is satisfied. With $`q=0`$, the free [Klein–Gordon equation](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) is recovered.

## Low-energy [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) to [Schrödinger](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation)

For constant $`m>0`$, the rest-energy phase is factored out as $`\psi'=e^{-\mathrm imt}\psi`$. By the product rule:

```math
\begin{aligned}
\hat E_q\psi'
&=(\mathrm i\partial_t-qV)(e^{-\mathrm imt}\psi)\\
&=e^{-\mathrm imt}(m\psi+\mathrm i\partial_t\psi-qV\psi)\\
&=e^{-\mathrm imt}(m+\hat E_q)\psi,\\[6pt]
\hat{{\color{cyan}𝐏}}_q\psi'&=e^{-\mathrm imt}\hat{{\color{cyan}𝐏}}_q\psi.
\end{aligned}
```

By substitution into [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation):

```math
\begin{aligned}
0&=(\hat E_q^2-\hat{{\color{cyan}𝐏}}_q^2-m^2)\psi'\\
&=e^{-\mathrm imt}\bigl((m+\hat E_q)^2-\hat{{\color{cyan}𝐏}}_q^2-m^2\bigr)\psi\\
&=e^{-\mathrm imt}\bigl(2m\hat E_q+\hat E_q^2-\hat{{\color{cyan}𝐏}}_q^2\bigr)\psi.
\end{aligned}
```

On the envelope $`\psi`$, energy above rest, excluding $`qV`$, is measured by $`\hat E_q`$. After phase cancellation:

```math
\hat E_q\psi=\frac{\hat{{\color{cyan}𝐏}}_q^2-\hat E_q^2}{2m}\psi.
```

For commuting $`\hat E_q`$ and $`\hat{{\color{cyan}𝐏}}_q^2`$, the positive-energy expansion is:

```math
\hat E_q\psi=\left(\frac{\hat{{\color{cyan}𝐏}}_q^2}{2m}
-\frac{\hat{{\color{cyan}𝐏}}_q^4}{8m^3}+\cdots\right)\psi.
```

For weak, slowly varying fields, [Schrödinger’s equation](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation) is obtained by neglecting $`\hat E_q^2\psi`$ relative to $`2m\hat E_q\psi`$:

```math
\boxed{\mathrm i\partial_t\psi=
\left(\frac{(-\mathrm i\partial_{{\color{cyan}𝐫}}-q{\color{cyan}𝐀})^2}{2m}+qV\right)\psi.}
```

With $`{\color{cyan}𝐀}=0`$ (hence $`{\color{cyan}𝐁}=0`$):

```math
\mathrm i\partial_t\psi=
\left(-\frac{\partial_{{\color{cyan}𝐫}}^2}{2m}+qV\right)\psi.
```

## Low-energy [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) to [Pauli](https://en.wikipedia.org/wiki/Pauli_equation)

With the rest-energy phase removed by $`\textcolor{magenta}\psi _\pm':=e^{\mathrm imt}\textcolor{magenta}\psi _\pm`$:

```math
(m+\hat E_q\pm\hat{{\color{cyan}𝐏}}_q\cdot{\color{cyan}𝝈})\textcolor{magenta}\psi _\pm'=m\textcolor{magenta}\psi _\mp'.
```

The positive-energy large and small components are defined by:

```math
\begin{aligned}
\textcolor{magenta}\phi _+&:=\frac{\textcolor{magenta}\psi _+'+\textcolor{magenta}\psi _-'}{\sqrt2},\\
\textcolor{magenta}\phi _-&:=\frac{\textcolor{magenta}\psi _-'-\textcolor{magenta}\psi _+'}{\sqrt2}.
\end{aligned}
```

By addition and subtraction:

```math
\begin{aligned}
\hat E_q\textcolor{magenta}\phi _+&=(\hat{{\color{cyan}𝐏}}_q\cdot{\color{cyan}𝝈})\textcolor{magenta}\phi _-,\\
(2m+\hat E_q)\textcolor{magenta}\phi _-&=(\hat{{\color{cyan}𝐏}}_q\cdot{\color{cyan}𝝈})\textcolor{magenta}\phi _+.
\end{aligned}
```

In the nonrelativistic limit, $`\hat E_q\textcolor{magenta}\phi _-`$ is neglected relative to $`2m\textcolor{magenta}\phi _-`$. By substitution and the field identity:

```math
\begin{aligned}
\textcolor{magenta}\phi _-&\simeq\frac{\hat{{\color{cyan}𝐏}}_q\cdot{\color{cyan}𝝈}}{2m}\textcolor{magenta}\phi _+,\\
\hat E_q\textcolor{magenta}\phi _+&\simeq\frac{(\hat{{\color{cyan}𝐏}}_q\cdot{\color{cyan}𝝈})^2}{2m}\textcolor{magenta}\phi _+
=\frac{\hat{{\color{cyan}𝐏}}_q^2-q{\color{cyan}𝐁}\cdot{\color{cyan}𝝈}}{2m}\textcolor{magenta}\phi _+.
\end{aligned}
```

With $`\hat E_q=\mathrm i\partial_t-qV`$, the [Pauli equation](https://en.wikipedia.org/wiki/Pauli_equation) is obtained through order $`1/m`$:

```math
\boxed{\mathrm i\partial_t\textcolor{magenta}\phi _+=
\left(\frac{\hat{{\color{cyan}𝐏}}_q^2}{2m}
+qV-\frac{q}{2m}{\color{cyan}𝐁}\cdot{\color{cyan}𝝈}\right)\textcolor{magenta}\phi _+.}
```

## [Spin up and spin down](https://en.wikipedia.org/wiki/Spin-1/2#Observables) in a uniform magnetic field

For a constant nonzero field $`{\color{cyan}𝐁}`$, its unit direction and [spin projectors](https://en.wikipedia.org/wiki/Pauli_matrices#Eigenvectors_and_eigenvalues) are defined by:

```math
{\color{cyan}𝐮}:=\frac{{\color{cyan}𝐁}}{\lVert{\color{cyan}𝐁}\rVert},
\qquad \textcolor{magenta}\Pi _\pm:=\frac12(1\pm{\color{cyan}𝐮}\cdot{\color{cyan}𝝈}).
```

By [spectral decomposition](#projections-and-spectral-decomposition):

```math
{\color{cyan}𝐁}\cdot{\color{cyan}𝝈}
=\lVert{\color{cyan}𝐁}\rVert(\textcolor{magenta}\Pi _+-\textcolor{magenta}\Pi _-),
\qquad
\textcolor{magenta}\Pi _\pm({\color{cyan}𝐁}\cdot{\color{cyan}𝝈})
=\pm\lVert{\color{cyan}𝐁}\rVert\textcolor{magenta}\Pi _\pm.
```

The two spin amplitudes $`\phi_{++},\phi_{+-}`$ are encoded in a paravector $`\textcolor{magenta}\phi _+`$ aligned with $`{\color{cyan}𝐮}`$:

```math
\begin{aligned}
\textcolor{magenta}\phi _+&=\phi_{++}\textcolor{magenta}\Pi _++\phi_{+-}\textcolor{magenta}\Pi _-,\\
\textcolor{magenta}\Pi _\pm\textcolor{magenta}\phi _+&=\textcolor{magenta}\phi _+\textcolor{magenta}\Pi _\pm=\phi_{+\pm}\textcolor{magenta}\Pi _\pm.
\end{aligned}
```

Constant $`\textcolor{magenta}\Pi _\pm`$ commute with $`\hat E_q`$ and $`\hat{{\color{cyan}𝐏}}_q^2`$, so the following is obtained from the [Pauli equation](https://en.wikipedia.org/wiki/Pauli_equation):

```math
\begin{aligned}
0&=\textcolor{magenta}\Pi _\pm\left(2m\hat E_q-\hat{{\color{cyan}𝐏}}_q^2
+q{\color{cyan}𝐁}\cdot{\color{cyan}𝝈}\right)\textcolor{magenta}\phi _+\\
&=\left(\left(2m\hat E_q-\hat{{\color{cyan}𝐏}}_q^2
\pm q\lVert{\color{cyan}𝐁}\rVert\right)\phi_{+\pm}\right)\textcolor{magenta}\Pi _\pm.
\end{aligned}
```

With $`\textcolor{magenta}\Pi _\pm\ne0`$ and $`\hat E_q=\mathrm i\partial_t-qV`$:

```math
\boxed{\mathrm i\partial_t\phi_{+\pm}=
\left(\frac{\hat{{\color{cyan}𝐏}}_q^2\mp q\lVert{\color{cyan}𝐁}\rVert}{2m}+qV\right)\phi_{+\pm}.}
```

Opposite [Zeeman shifts](https://en.wikipedia.org/wiki/Zeeman_effect) are obtained for [spin projections](https://en.wikipedia.org/wiki/Spin-1/2) $`\pm\tfrac12`$ along $`{\color{cyan}𝐮}`$. Orbital coupling is retained in [kinetic momentum](https://en.wikipedia.org/wiki/Momentum_operator#Definition_(position_space)) $`\hat{{\color{cyan}𝐏}}_q=-\mathrm i\partial_{{\color{cyan}𝐫}}-q{\color{cyan}𝐀}`$.

This encoding is used for uniform-field [Pauli evolution](https://en.wikipedia.org/wiki/Pauli_equation); the fixed spinor space is retained for [Dirac dynamics](https://en.wikipedia.org/wiki/Dirac_equation) and measurements.

With the rest-energy phase removed, the four complex [Dirac amplitudes](https://en.wikipedia.org/wiki/Dirac_spinor) are arranged as:

```math
{\color{cyan}𝝓}:=\left(\,\begin{matrix}
\phi_{++}&\phi_{+-}\\
\phi_{-+}&\phi_{--}
\end{matrix}\,\right).
```

Rows are [Dirac blocks](https://en.wikipedia.org/wiki/Dirac_spinor); columns are spin projections. The upper row is governed by the [Pauli limit](https://en.wikipedia.org/wiki/Pauli_equation).

## [Spin measurement](https://en.wikipedia.org/wiki/Spin-1/2#Rotations_and_Spinors) probabilities

Real unit preparation and measurement axes are denoted by $`{\color{cyan}𝐮}',{\color{cyan}𝐮}`$, with $`\cos\theta:={\color{cyan}𝐮}\cdot{\color{cyan}𝐮}'`$.

The state $`\textcolor{magenta}\psi =\textcolor{magenta}\Pi _+'\textcolor{magenta}\psi \ne0`$ is prepared using projectors with the following overlap:

```math
\begin{gathered}
\textcolor{magenta}\Pi _+':=\frac12(1+{\color{cyan}𝐮}'\cdot{\color{cyan}𝝈}),\qquad
\textcolor{magenta}\Pi _\pm:=\frac12(1\pm{\color{cyan}𝐮}\cdot{\color{cyan}𝝈}),\\[4pt]
\textcolor{magenta}\Pi _+'\textcolor{magenta}\Pi _\pm\textcolor{magenta}\Pi _+'=\frac{1\pm\cos\theta}{2}\textcolor{magenta}\Pi _+'.
\end{gathered}
```

Using the squared norm, $`\lVert {\color{magenta}\boldsymbol{\mathcal Z}}\rVert^2=\mathrm{sc}({\color{magenta}\boldsymbol{\mathcal Z}}{\color{magenta}\boldsymbol{\mathcal Z}}^{\mathsf H})`$, the denominator is:

```math
\lVert\textcolor{magenta}\psi \rVert^2=\mathrm{sc}(\textcolor{magenta}\psi \textcolor{magenta}\psi ^{\mathsf H})>0.
```

For the [prepared state](https://en.wikipedia.org/wiki/Spin-1/2#Bloch_Representation), the projected squared norm is:

```math
\begin{aligned}
\lVert\textcolor{magenta}\Pi _\pm\textcolor{magenta}\psi \rVert^2
&=\mathrm{sc}\bigl((\textcolor{magenta}\Pi _\pm\textcolor{magenta}\psi )(\textcolor{magenta}\Pi _\pm\textcolor{magenta}\psi )^{\mathsf H}\bigr)\\
&=\mathrm{sc}(\textcolor{magenta}\Pi _\pm\textcolor{magenta}\psi \textcolor{magenta}\psi ^{\mathsf H})\\
&=\mathrm{sc}(\textcolor{magenta}\Pi _+'\textcolor{magenta}\Pi _\pm\textcolor{magenta}\Pi _+'\textcolor{magenta}\psi \textcolor{magenta}\psi ^{\mathsf H})\\
&=\frac{1\pm\cos\theta}{2}\mathrm{sc}(\textcolor{magenta}\psi \textcolor{magenta}\psi ^{\mathsf H})
=\frac{1\pm\cos\theta}{2}\lVert\textcolor{magenta}\psi \rVert^2.
\end{aligned}
```

The [Born probabilities](https://en.wikipedia.org/wiki/Born_rule) for [spin projections](https://en.wikipedia.org/wiki/Spin-1/2) $`\pm\tfrac12`$ along $`{\color{cyan}𝐮}`$ are:

```math
\begin{gathered}
\mathrm{prob}(\pm)=\frac{\lVert\textcolor{magenta}\Pi _\pm\textcolor{magenta}\psi \rVert^2}{\lVert\textcolor{magenta}\psi \rVert^2}
=\frac{1\pm\cos\theta}{2},\\[6pt]
\boxed{\mathrm{prob}(+)=\cos^2\frac\theta2,\quad \mathrm{prob}(-)=\sin^2\frac\theta2.}
\end{gathered}
```
