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

## Six core operations

[Determinant](https://en.wikipedia.org/wiki/Determinant), [trace](https://en.wikipedia.org/wiki/Trace_%28linear_algebra%29), [squared norm](https://en.wikipedia.org/wiki/Norm_%28mathematics%29), [inverse](https://en.wikipedia.org/wiki/Invertible_matrix), and the two conjugations below.

Write $`Z=S+\mathbf V\cdot\boldsymbol{\sigma}`$, where $`S\in\mathbb C`$ is a complex scalar and $`\mathbf V\in\mathbb C^3`$ is a complex vector.

On coefficients, $`^{*}`$ means ordinary [complex conjugation](https://en.wikipedia.org/wiki/Complex_conjugate). On a paravector, it also reverses the vector sign; $`^{\mathsf H}`$ is [Hermitian conjugation](https://en.wikipedia.org/wiki/Conjugate_transpose).

```math
\begin{aligned}
\det Z&=S^2-\mathbf V\cdot\mathbf V,
&\mathrm{tr}\,Z&=2S,\\
\lVert Z\rVert^2&=S^{*}S+\mathbf V^{*}\cdot\mathbf V,
&Z^{-1}&=\frac{S-\mathbf V\cdot\boldsymbol{\sigma}}{\det Z}\quad(\det Z\ne0),\\
Z^{*}&=S^{*}-\mathbf V^{*}\cdot\boldsymbol{\sigma},
&Z^{\mathsf H}&=S^{*}+\mathbf V^{*}\cdot\boldsymbol{\sigma}.
\end{aligned}
```

## Rotations, boosts, and projections

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

[Projections](https://en.wikipedia.org/wiki/Paravector#Null_paravectors_as_projectors) along $`\mathbf u`$:

```math
\Pi_\pm=\frac12(1\pm\mathbf u\cdot\boldsymbol{\sigma}),
\qquad \Pi_\pm^2=\Pi_\pm,\qquad\Pi_+\Pi_-=0.
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

Equating its four components gives [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations):

```math
\begin{aligned}
\text{Real scalar:}\quad &\partial_{\mathbf r}\cdot\mathbf E=\rho,\\
\text{Imaginary scalar:}\quad &\partial_{\mathbf r}\cdot\mathbf B=0,\\
\text{Real vector:}\quad &\partial_t\mathbf E-\partial_{\mathbf r}\times\mathbf B=-\mathbf J,\\
\text{Imaginary vector:}\quad &\partial_t\mathbf B+\partial_{\mathbf r}\times\mathbf E=0.
\end{aligned}
```

Applying $`\partial^{*}=\partial_t-\partial_{\mathbf r}\cdot\boldsymbol{\sigma}`$ gives the [d'Alembertian](https://en.wikipedia.org/wiki/D%27Alembert_operator) and the [sourced wave equations](https://en.wikipedia.org/wiki/Wave_equation):

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

## [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) as a first-order wave equation

Package two paravector [amplitudes](https://en.wikipedia.org/wiki/Probability_amplitude) into $`\boldsymbol\Psi=(\Psi_1,\Psi_2)^{\mathsf T}`$, each with two independent complex components. With $`\mathbf I`$ the block identity, define $`\mathbf W`$ in the [Weyl (chiral) representation](https://en.wikipedia.org/wiki/Gamma_matrices#Weyl_%28chiral%29_basis):

```math
\mathbf W(Z)=
\begin{pmatrix}
0 & Z \\
Z^{*\mathsf H} & 0
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

**With an [electromagnetic potential](https://en.wikipedia.org/wiki/Electromagnetic_four-potential)**, $`\Phi=V+\mathbf A\cdot\boldsymbol{\sigma}`$ and [charge](https://en.wikipedia.org/wiki/Electric_charge) $`q`$, [minimal coupling](https://en.wikipedia.org/wiki/Minimal_coupling) gives

```math
\boxed{\bigl(\mathbf W(\mathrm{i}\partial-q\Phi^{*})-m\mathbf I\bigr)\boldsymbol\Psi=0.}
```

## Low-energy [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) recovers [Pauli](https://en.wikipedia.org/wiki/Pauli_equation)

Remove the [rest phase](https://en.wikipedia.org/wiki/Pauli_equation#Derivation), $`\Psi_k=e^{-\mathrm{i}mt}\Psi_k'`$, and define

```math
\hat E=\mathrm{i}\partial_t-qV,\qquad
\hat{\mathbf p}=-\mathrm{i}\partial_{\mathbf r}-q\mathbf A.
```

The exact energy recursion includes the [commutator](https://en.wikipedia.org/wiki/Commutator):

```math
\begin{aligned}
\hat E\Psi_2'&=\frac{(\hat{\mathbf p}\cdot\boldsymbol{\sigma})^2-\hat E^2
-[\hat{\mathbf p}\cdot\boldsymbol{\sigma},\hat E]}{2m}\Psi_2',\\
[\hat{\mathbf p}\cdot\boldsymbol{\sigma},\hat E]&=-\mathrm{i}q\mathbf E\cdot\boldsymbol{\sigma},
\qquad [A,B]:=AB-BA.
\end{aligned}
```

With $`\phi_1=(\Psi_1'+\Psi_2')/\sqrt2`$ and $`\phi_2=(\Psi_2'-\Psi_1')/\sqrt2`$,

```math
\hat E\phi_1=(\hat{\mathbf p}\cdot\boldsymbol{\sigma})\phi_2,
\qquad
\phi_2=\frac{\hat{\mathbf p}\cdot\boldsymbol{\sigma}}{2m}\phi_1-\frac{\hat E\phi_2}{2m}.
```

For slow positive-energy envelopes in weak, slowly varying fields, iterate once:

```math
\phi_2\simeq\frac{\hat{\mathbf p}\cdot\boldsymbol{\sigma}}{2m}\phi_1,
\qquad
\hat E\phi_1\simeq\frac{(\hat{\mathbf p}\cdot\boldsymbol{\sigma})^2}{2m}\phi_1
=\frac{\hat{\mathbf p}^2-q\mathbf B\cdot\boldsymbol{\sigma}}{2m}\phi_1.
```

This gives [Pauli's equation](https://en.wikipedia.org/wiki/Pauli_equation) to leading order:

```math
\boxed{\mathrm{i}\partial_t\phi_1=
\left[\frac{\hat{\mathbf p}^2}{2m}
+qV-\frac{q}{2m}\mathbf B\cdot\boldsymbol{\sigma}\right]\phi_1.}
```
