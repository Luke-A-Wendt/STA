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
𝝈&:=\{\sigma_1,\sigma_2,\sigma_3\},
&\mathrm{i}𝝈&=\{\sigma_2\sigma_3,\sigma_3\sigma_1,\sigma_1\sigma_2\},\\
Z&:=(a+b\mathrm{i})+(𝐀+\mathrm{i}𝐁)\cdot𝝈,
&a,b&\in\mathbb R,\quad𝐀,𝐁\in\mathbb R^3.
\end{aligned}}
```

These eight real components form a complex scalar and a complex vector:

```math
Z=S+𝐕\cdot𝝈=X+\mathrm{i}Y,
\qquad X:=a+𝐀\cdot𝝈,
\quad Y:=b+𝐁\cdot𝝈.
```

The component extractions are

```math
\begin{aligned}
\mathrm{re}(Z)&:=X,\\
\mathrm{im}(Z)&:=Y,\\
\mathrm{sc}(Z)&:=S=a+\mathrm{i}b,\\
\mathrm{vec}(Z)&:=𝐕=𝐀+\mathrm{i}𝐁.
\end{aligned}
```

## Product

The product carries both the [dot](https://en.wikipedia.org/wiki/Dot_product) and [cross](https://en.wikipedia.org/wiki/Cross_product) products:

```math
(𝐀\cdot𝝈)(𝐁\cdot𝝈)
=𝐀\cdot𝐁+\mathrm{i}(𝐀\times𝐁)\cdot𝝈.
```

Here we use the convention $`𝐀^2:=𝐀\cdot𝐀`$.

The [commutator](https://en.wikipedia.org/wiki/Commutator) is

```math
[X,Y]:=XY-YX.
```

## Core operations

[Determinant](https://en.wikipedia.org/wiki/Determinant), [trace](https://en.wikipedia.org/wiki/Trace_%28linear_algebra%29), [squared norm](https://en.wikipedia.org/wiki/Norm_%28mathematics%29), adjugate, [inverse](https://en.wikipedia.org/wiki/Invertible_matrix), and the two conjugations below.

Write $`Z=S+𝐕\cdot𝝈`$, where $`S\in\mathbb C`$ is a complex scalar and $`𝐕\in\mathbb C^3`$ is a complex vector.

$`(\,)^{*}:=`$ [complex conjugation](https://en.wikipedia.org/wiki/Complex_conjugate) of complex numbers and a sign flip of $`𝝈`$.

$`(\,)^{\mathsf H}:=`$ [Hermitian conjugation](https://en.wikipedia.org/wiki/Conjugate_transpose), which flips the order of multiplied $`\sigma_k`$.

```math
\begin{gathered}
\begin{aligned}
Z^{*}&=S^{*}-𝐕^{*}\cdot𝝈,
&Z^{\mathsf H}&=S^{*}+𝐕^{*}\cdot𝝈,\\
\mathrm{adj}\,Z&:=S-𝐕\cdot𝝈 = Z^{*\mathsf H},&\det Z&:=S^2-𝐕^2 = Z\,\mathrm{adj}\,Z,\\
\mathrm{tr}\,Z&:=2S=Z+\mathrm{adj}\,Z,
&\lVert Z\rVert^2&:=SS^{*}+𝐕\cdot𝐕^{*} = \mathrm{sc}(Z Z^\mathsf H).
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
Q&:=a+𝐀\cdot𝐪
=a-\mathrm{i}𝐀\cdot𝝈,
\qquad a\in\mathbb R,\quad𝐀\in\mathbb R^3,\\
Q^{\mathsf H}&=a+\mathrm{i}𝐀\cdot𝝈,\\
QQ^{\mathsf H}&=\det Q=a^2+𝐀^2.
\end{aligned}
```

## Rotations and boosts

For real spacetime $`X`$ and real [unit axis](https://en.wikipedia.org/wiki/Unit_vector) $`𝐮`$:

```math
X:=t+𝐫\cdot𝝈,
\qquad X':=TXT^{\mathsf H}=t'+𝐫'\cdot𝝈.
```

For $`T=R`$: [Rodrigues rotation](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula) through $`\theta`$, with [unit quaternion](https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation) $`Q=R`$, $`q_k=-\mathrm{i}\sigma_k`$.

```math
\begin{aligned}
R&:=e^{-\mathrm{i}\theta𝐮\cdot𝝈/2}
=\cos\frac\theta2-\mathrm{i}𝐮\cdot𝝈\sin\frac\theta2,\\
X'&=RXR^{\mathsf H},\\
t'&=t,\\
𝐫'&=𝐫\cos\theta+(𝐮\times𝐫)\sin\theta
+(𝐮\cdot𝐫)(1-\cos\theta)𝐮.
\end{aligned}
```

For $`T=L`$: [Lorentz boost](https://en.wikipedia.org/wiki/Lorentz_transformation) to a frame moving at $`+\beta𝐮`$, with [rapidity](https://en.wikipedia.org/wiki/Rapidity) $`\theta`$, $`\beta:=\tanh\theta`$, and [Lorentz factor](https://en.wikipedia.org/wiki/Lorentz_factor) $`\gamma:=\cosh\theta`$.

Define the parallel and perpendicular components:

```math
𝐫_\parallel:=(𝐮\cdot𝐫)𝐮,
\qquad 𝐫_\perp:=𝐫-𝐫_\parallel.
```

```math
\begin{aligned}
L&:=e^{-\theta𝐮\cdot𝝈/2}
=\cosh\frac\theta2-𝐮\cdot𝝈\sinh\frac\theta2
=L^{\mathsf H},\\
X'&=LXL^{\mathsf H},\\
t'&=\gamma(t-\beta𝐮\cdot𝐫),\\
𝐫'&=𝐫_\perp+\gamma(𝐫_\parallel-\beta t𝐮).
\end{aligned}
```

## Projections and spectral decomposition

```math
Z=a+b𝐮\cdot𝝈,
\qquad a,b\in\mathbb C,\quad 𝐮\in\mathbb C^3,\quad 𝐮^2=1.
```

The eigenvalues and complementary [projectors](https://en.wikipedia.org/wiki/Paravector#Null_paravectors_as_projectors) are

```math
\begin{gathered}
\lambda_\pm:=a\pm b,
\qquad \Pi_\pm:=\frac12(1\pm𝐮\cdot𝝈),\\[4pt]
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

For a future-directed massive particle, proper time and four-momentum follow from a real paravector, with four-velocity $`U`$.

Spacetime increment and proper-time interval:

```math
\begin{aligned}
\mathrm dX&=\mathrm dt+\mathrm d𝐫\cdot𝝈,\\
\det(\mathrm dX)&=\mathrm dt^2-\mathrm d𝐫^2.
\end{aligned}
```

Choose positive [proper time](https://en.wikipedia.org/wiki/Proper_time) along the future-directed path $`\mathrm dt>0`$.

```math
\mathrm ds:=\sqrt{\mathrm dt^2-\mathrm d𝐫^2}=\sqrt{\det(\mathrm dX)}.
```

Coordinate velocity:

```math
\begin{aligned}
\frac{\mathrm ds}{\mathrm dt}
&=\sqrt{\frac{\mathrm dt^2-\mathrm d𝐫^2}{\mathrm dt^2}}\\
&=\sqrt{1-\left(\frac{\mathrm d𝐫}{\mathrm dt}\right)^2}
=\sqrt{1-𝐯^2},
\qquad 𝐯:=\frac{\mathrm d𝐫}{\mathrm dt}.
\end{aligned}
```

```math
\gamma:=\frac{\mathrm dt}{\mathrm ds}=\frac{1}{\sqrt{1-𝐯^2}}.
```

[Four-velocity](https://en.wikipedia.org/wiki/Four-velocity):

```math
\begin{aligned}
U:=\frac{\mathrm dX}{\mathrm ds}
&=\frac{\mathrm dt+\mathrm d𝐫\cdot𝝈}{\mathrm ds}\\
&=\left(1+\frac{\mathrm d𝐫}{\mathrm dt}\cdot𝝈\right)
\frac{\mathrm dt}{\mathrm ds}\\
&=\gamma(1+𝐯\cdot𝝈).
\end{aligned}
```

Proper-time normalization gives

```math
\det U=\det\!\left(\frac{\mathrm dX}{\mathrm ds}\right)
=\frac{\det(\mathrm dX)}{\mathrm ds^2}=1.
```

For energy $`E`$ and momentum $`𝐏`$, the [four-momentum](https://en.wikipedia.org/wiki/Four-momentum) is

```math
P:=mU=E+𝐏\cdot𝝈.
```

For constant mass $`m>0`$, the determinant scales quadratically:

```math
\det P=\det(mU)=m^2\det U=m^2.
```

```math
\det P=E^2-𝐏^2=m^2.
```

## [Maxwell](https://en.wikipedia.org/wiki/Maxwell%27s_equations) in one equation

Combine the [electric field](https://en.wikipedia.org/wiki/Electric_field) $`𝐄`$ and [magnetic field](https://en.wikipedia.org/wiki/Magnetic_field) $`𝐁`$ into the complex electromagnetic field $`𝐅`$ and its paravector $`F`$:

```math
𝐅:=𝐄+\mathrm{i}𝐁,\qquad F:=𝐅\cdot𝝈.
```

Combine the [charge density](https://en.wikipedia.org/wiki/Charge_density) $`\rho`$ and [current density](https://en.wikipedia.org/wiki/Current_density) $`𝐉`$ into the paravector density $`J`$:

```math
J:=\rho+𝐉\cdot𝝈,\qquad J^*=\rho-𝐉\cdot𝝈.
```

Define the [paravector derivative](https://en.wikipedia.org/wiki/Paravector#Paragradient), with [gradient](https://en.wikipedia.org/wiki/Gradient) $`\partial_{𝐫}`$:

```math
\partial:=\partial_t+\partial_{𝐫}\cdot𝝈,\qquad
\partial^*=\partial_t-\partial_{𝐫}\cdot𝝈.
```

[Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations) are given by

```math
\boxed{\partial F=J^*}
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

Applying $`\partial^{*}`$ gives the [d'Alembertian](https://en.wikipedia.org/wiki/D%27Alembert_operator) and the [sourced wave equations](https://en.wikipedia.org/wiki/Inhomogeneous_electromagnetic_wave_equation):

```math
\begin{aligned}
\Box&:=\partial^{*}\partial=\partial_t^2-\partial_{𝐫}^2,\\
\partial^{*}(\partial F)&=\Box F=\partial^{*}J^*.
\end{aligned}
```

```math
\boxed{\Box F^*=\partial J}
```

This boxed equation gives [charge conservation](https://en.wikipedia.org/wiki/Charge_conservation) and the field wave equations component by component:

```math
\begin{aligned}
\text{Real scalar:}\quad &0=\partial_t\rho+\partial_{𝐫}\cdot𝐉
&&\text{(charge conservation)},\\
\text{Imaginary scalar:}\quad &0=0,\\
\text{Real vector:}\quad &\Box𝐄=-\partial_t𝐉-\partial_{𝐫}\rho,\\
\text{Imaginary vector:}\quad &\Box𝐁=\partial_{𝐫}\times𝐉.
\end{aligned}
```

In vacuum, $`\Box F=0`$.

## [Electromagnetic potentials](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) and [gauge invariance](https://en.wikipedia.org/wiki/Electromagnetic_four-potential#Gauge_freedom)

Define the real [four-potential](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) from the [scalar potential](https://en.wikipedia.org/wiki/Electric_potential) $`V`$ and [vector potential](https://en.wikipedia.org/wiki/Magnetic_vector_potential) $`𝐀`$, with gauge scalar $`S`$:

```math
A:=V+𝐀\cdot𝝈,
\qquad S:=\mathrm{sc}(\partial A).
```

Fields from the potential:

```math
\boxed{F^{*}=\partial A-S}
```

Its components give

```math
\begin{aligned}
\text{Real scalar:}\quad &S=\partial_tV+\partial_{𝐫}\cdot𝐀,\\
\text{Imaginary scalar:}\quad &0=0,\\
\text{Real vector:}\quad &𝐄=-\partial_t𝐀-\partial_{𝐫}V,\\
\text{Imaginary vector:}\quad &𝐁=\partial_{𝐫}\times𝐀.
\end{aligned}
```

Source equations in any [gauge](https://en.wikipedia.org/wiki/Gauge_fixing):

```math
\Box A=\partial^{*}(\partial A)
=\partial^{*}(S+F^{*})
=\partial^{*}S+(\partial F)^{*}.
```

Using Maxwell’s equation gives

```math
\boxed{\Box A=\partial^{*}S+J}
```

Its components give

```math
\begin{aligned}
\text{Real scalar:}\quad &\Box V-\partial_tS=\rho,\\
\text{Imaginary scalar:}\quad &0=0,\\
\text{Real vector:}\quad &\Box𝐀+\partial_{𝐫}S=𝐉,\\
\text{Imaginary vector:}\quad &\mathbf 0=\mathbf 0.
\end{aligned}
```

[Gauge invariance](https://en.wikipedia.org/wiki/Electromagnetic_four-potential#Gauge_freedom): for real scalar $`\lambda`$,

```math
\begin{aligned}
A'&:=A+\partial^{*}\lambda,\\
S'&=S+\Box\lambda,\\
F'&=F.
\end{aligned}
```

[Lorenz gauge](https://en.wikipedia.org/wiki/Lorenz_gauge_condition): $`S=0`$ gives the [potential wave equations](https://en.wikipedia.org/wiki/Inhomogeneous_electromagnetic_wave_equation)

```math
\begin{aligned}
\text{Real scalar:}\quad &\Box V=\rho,\\
\text{Real vector:}\quad &\Box𝐀=𝐉.
\end{aligned}
```

## [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) from the mass shell

Hats mark named [energy](https://en.wikipedia.org/wiki/Energy_operator), [momentum](https://en.wikipedia.org/wiki/Momentum_operator), and [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) operators; the derivative symbols $`\partial`$ and $`\Box`$ remain unhatted.

Represent energy and momentum on a complex scalar wavefunction $`\psi`$ by

```math
\begin{aligned}
\hat E&:=\mathrm{i}\partial_t,\\
\hat{𝐏}&:=-\mathrm{i}\partial_{𝐫},\\
\hat{P}&:=\hat E+\hat{𝐏}\cdot𝝈 = \mathrm{i}\partial^{*}.
\end{aligned}
```

The free operators commute, so

```math
\hat{P}\,\mathrm{adj}\,\hat{P}
=\det\hat{P}
=\hat E^2-\hat{𝐏}^2
=-\Box.
```

Imposing the mass shell on $`\psi`$ gives

```math
\begin{aligned}
0&=(\det\hat{P}-m^2)\psi\\
&=(\hat E^2-\hat{𝐏}^2-m^2)\psi\\
&=-(\Box+m^2)\psi.
\end{aligned}
```

Thus the free Klein–Gordon equation is

```math
\boxed{(\Box+m^2)\psi=0.}
```

## [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) as a first-order wave equation

In the [Weyl representation](https://en.wikipedia.org/wiki/Gamma_matrices#Weyl_%28chiral%29_basis), write $`𝝍:=(\psi_+,\psi_-)^{\mathsf T}`$, with two complex components in each entry. Let $`𝐈`$ be the block identity.

Using the free momentum operator above, define

```math
\begin{aligned}
𝐖(\hat{P})&:=
\begin{pmatrix}
0 & \mathrm{adj}\,\hat{P} \\
\hat{P} & 0
\end{pmatrix},\\
\hat{𝐃}(m)&:=𝐖(\hat{P})-m𝐈.
\end{aligned}
```

The Dirac equation is

```math
\boxed{\hat{𝐃}(m)𝝍=0.}
```

Equivalently,

```math
(\hat E\pm\hat{𝐏}\cdot𝝈)\psi_\pm=m\psi_\mp.
```

Using $`𝐖(\hat{P})^2=(\det\hat{P})𝐈`$ gives the [Klein–Gordon factorization](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation):

```math
\begin{aligned}
\hat{𝐃}(m)\hat{𝐃}(-m)&=\hat{𝐃}(-m)\hat{𝐃}(m)\\
&=𝐖(\hat{P})^2-m^2𝐈
=-(\Box+m^2)𝐈.
\end{aligned}
```

Hence every free Dirac solution satisfies Klein–Gordon.

## [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) with an electromagnetic potential

For constant mass $`m>0`$ and [charge](https://en.wikipedia.org/wiki/Electric_charge) $`q`$, [minimal coupling](https://en.wikipedia.org/wiki/Minimal_coupling) to the real [potential](https://en.wikipedia.org/wiki/Electromagnetic_four-potential) $`A=V+𝐀\cdot𝝈`$ defines the gauge-covariant energy and momentum operators

```math
\begin{aligned}
\hat E_q&:=\hat E-qV=\mathrm{i}\partial_t-qV,\\
\hat{𝐏}_q&:=\hat{𝐏}-q𝐀=-\mathrm{i}\partial_{𝐫}-q𝐀,\\
\hat{P}_q&:=\hat E_q+\hat{𝐏}_q\cdot𝝈=\mathrm{i}\partial^{*}-qA.
\end{aligned}
```

```math
\boxed{\hat{P}_q=\hat{P}-qA.}
```

At $`q=0`$, the coupled operators reduce to the canonical ones: $`\hat E_0=\hat E`$, $`\hat{𝐏}_0=\hat{𝐏}`$, and $`\hat{P}_0=\hat{P}`$. The coupled Dirac operator is

```math
\begin{aligned}
\hat{𝐃}_q(m)&:=𝐖(\hat{P}_q)-m𝐈\\
&=\begin{pmatrix}
-m & \mathrm{adj}\,\hat{P}_q\\
\hat{P}_q & -m
\end{pmatrix}.
\end{aligned}
```

The coupled Dirac equation is

```math
\boxed{\hat{𝐃}_q(m)𝝍=0.}
```

Equivalently,

```math
(\hat E_q\pm\hat{𝐏}_q\cdot𝝈)\psi_\pm=m\psi_\mp.
```

The coupled operators obey these field identities:

```math
\begin{aligned}
{}[\hat{𝐏}_q\cdot𝝈,\hat E_q]
&=-\mathrm{i}q𝐄\cdot𝝈,\\
(\hat{𝐏}_q\cdot𝝈)^2
&=\hat{𝐏}_q^2-q𝐁\cdot𝝈,\\
\hat{𝐏}_q^2+\mathrm{i}qF
&=(\hat{𝐏}_q\cdot𝝈)^2
-[\hat{𝐏}_q\cdot𝝈,\hat E_q].
\end{aligned}
```

When derivatives act on the potentials as well as the wavefunction, the opposite-mass product becomes

```math
\hat{𝐃}_q(-m)\hat{𝐃}_q(m)
=(\hat E_q^2-\hat{𝐏}_q^2-m^2)𝐈
-\mathrm{i}q\begin{pmatrix}F^{*}&0\\0&F\end{pmatrix}.
```

Every Dirac solution satisfies $`\hat{𝐃}_q(-m)\hat{𝐃}_q(m)𝝍=0`$. Setting $`q=0`$ recovers the free Klein–Gordon equation.

## Low-energy [Klein–Gordon](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation) recovers [Schrödinger](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation)

For constant $`m`$, remove the rest-energy phase, $`\psi'=e^{-\mathrm{i}mt}\psi`$. The product rule gives

```math
\begin{aligned}
\hat E_q\psi'
&=(\mathrm{i}\partial_t-qV)(e^{-\mathrm{i}mt}\psi)\\
&=e^{-\mathrm{i}mt}(m\psi+\mathrm{i}\partial_t\psi-qV\psi)\\
&=e^{-\mathrm{i}mt}(m+\hat E_q)\psi,\\[6pt]
\hat{𝐏}_q\psi'&=e^{-\mathrm{i}mt}\hat{𝐏}_q\psi.
\end{aligned}
```

Substitute into Klein–Gordon:

```math
\begin{aligned}
0&=(\hat E_q^2-\hat{𝐏}_q^2-m^2)\psi'\\
&=e^{-\mathrm{i}mt}\bigl((m+\hat E_q)^2-\hat{𝐏}_q^2-m^2\bigr)\psi\\
&=e^{-\mathrm{i}mt}\bigl(2m\hat E_q+\hat E_q^2-\hat{𝐏}_q^2\bigr)\psi.
\end{aligned}
```

For any $`q`$ and constant $`m>0`$, canceling the phase gives the exact recursion:

```math
\hat E_q\psi=\frac{\hat{𝐏}_q^2-\hat E_q^2}{2m}\psi.
```

If $`\hat E_q`$ and $`\hat{𝐏}_q^2`$ commute, iteration on the positive-energy branch gives

```math
\hat E_q\psi=\left(\frac{\hat{𝐏}_q^2}{2m}
-\frac{\hat{𝐏}_q^4}{8m^3}+\cdots\right)\psi.
```

For weak, slowly varying fields, neglect $`\hat E_q^2\psi`$ relative to $`2m\hat E_q\psi`$ to obtain Schrödinger’s equation:

```math
\boxed{\mathrm{i}\partial_t\psi=
\left(\frac{(-\mathrm{i}\partial_{𝐫}-q𝐀)^2}{2m}+qV\right)\psi.}
```

With $`𝐀=0`$ (hence $`𝐁=0`$), this recovers the conventional Schrödinger equation:

```math
\mathrm{i}\partial_t\psi=
\left(-\frac{\partial_{𝐫}^2}{2m}+qV\right)\psi.
```

## Low-energy [Dirac](https://en.wikipedia.org/wiki/Dirac_equation) recovers [Pauli](https://en.wikipedia.org/wiki/Pauli_equation)

Using the same phase shift, $`\psi_\pm':=e^{\mathrm{i}mt}\psi_\pm`$, gives

```math
(m+\hat E_q\pm\hat{𝐏}_q\cdot𝝈)\psi_\pm'=m\psi_\mp'.
```

Define the large and small components for the positive-energy branch:

```math
\begin{aligned}
\phi_+&:=\frac{\psi_+'+\psi_-'}{\sqrt2},\\
\phi_-&:=\frac{\psi_-'-\psi_+'}{\sqrt2}.
\end{aligned}
```

Adding and subtracting the envelope equations gives

```math
\begin{aligned}
\hat E_q\phi_+&=(\hat{𝐏}_q\cdot𝝈)\phi_-,\\
(2m+\hat E_q)\phi_-&=(\hat{𝐏}_q\cdot𝝈)\phi_+.
\end{aligned}
```

In the same nonrelativistic limit, neglect $`\hat E_q\phi_-`$ relative to $`2m\phi_-`$. Substitute the resulting small component into the first equation and use the field identity above:

```math
\begin{aligned}
\phi_-&\simeq\frac{\hat{𝐏}_q\cdot𝝈}{2m}\phi_+,\\
\hat E_q\phi_+&\simeq\frac{(\hat{𝐏}_q\cdot𝝈)^2}{2m}\phi_+
=\frac{\hat{𝐏}_q^2-q𝐁\cdot𝝈}{2m}\phi_+.
\end{aligned}
```

Restoring $`\hat E_q=\mathrm{i}\partial_t-qV`$ gives the [Pauli equation](https://en.wikipedia.org/wiki/Pauli_equation), retaining the kinetic and spin terms through order $`1/m`$:

```math
\boxed{\mathrm{i}\partial_t\phi_+=
\left(\frac{\hat{𝐏}_q^2}{2m}
+qV-\frac{q}{2m}𝐁\cdot𝝈\right)\phi_+.}
```

## [Spin up and spin down](https://en.wikipedia.org/wiki/Spin-1/2#Observables) in a uniform magnetic field

For a constant nonzero field $`𝐁`$, define its unit direction and [spin projectors](https://en.wikipedia.org/wiki/Pauli_matrices#Eigenvectors_and_eigenvalues):

```math
𝐮:=\frac{𝐁}{\lVert𝐁\rVert},
\qquad \Pi_\pm:=\frac12(1\pm𝐮\cdot𝝈).
```

[Spectral decomposition](#projections-and-spectral-decomposition) gives

```math
𝐁\cdot𝝈
=\lVert𝐁\rVert(\Pi_+-\Pi_-),
\qquad
\Pi_\pm(𝐁\cdot𝝈)
=\pm\lVert𝐁\rVert\Pi_\pm.
```

For a paravector $`\phi_+`$ with vector part along $`𝐮`$, spectral decomposition gives complex scalar wavefunctions:

```math
\begin{aligned}
\phi_+&=\phi_{++}\Pi_++\phi_{+-}\Pi_-,\\
\Pi_\pm\phi_+&=\phi_+\Pi_\pm=\phi_{+\pm}\Pi_\pm.
\end{aligned}
```

Constant $`\Pi_\pm`$ commute with $`\hat E_q`$ and $`\hat{𝐏}_q^2`$, so the [Pauli equation](https://en.wikipedia.org/wiki/Pauli_equation) gives

```math
\begin{aligned}
0&=\Pi_\pm\left(2m\hat E_q-\hat{𝐏}_q^2
+q𝐁\cdot𝝈\right)\phi_+\\
&=\left(\left(2m\hat E_q-\hat{𝐏}_q^2
\pm q\lVert𝐁\rVert\right)\phi_{+\pm}\right)\Pi_\pm.
\end{aligned}
```

Since $`\Pi_\pm\ne0`$, the scalar coefficients vanish. Restoring $`\hat E_q=\mathrm{i}\partial_t-qV`$ gives

```math
\boxed{\mathrm{i}\partial_t\phi_{+\pm}=
\left(\frac{\hat{𝐏}_q^2\mp q\lVert𝐁\rVert}{2m}+qV\right)\phi_{+\pm}.}
```

[Spin](https://en.wikipedia.org/wiki/Spin-1/2) $`\pm\tfrac12`$ along $`𝐮`$ has opposite [Zeeman shifts](https://en.wikipedia.org/wiki/Zeeman_effect). The [kinetic momentum](https://en.wikipedia.org/wiki/Momentum_operator#Electromagnetic_field) $`\hat{𝐏}_q=-\mathrm{i}\partial_{𝐫}-q𝐀`$ retains the orbital coupling.

This representation applies to uniform-field Pauli evolution; the full Dirac equation and spin measurements use the fixed spinor space.

The full [Dirac spinor](https://en.wikipedia.org/wiki/Dirac_spinor) has four complex scalar amplitudes:

```math
𝝓:=\left(\,\begin{matrix}
\phi_{++}&\phi_{+-}\\
\phi_{-+}&\phi_{--}
\end{matrix}\,\right).
```

Rows label Dirac blocks, columns spin; the Pauli limit governs the upper row.

## [Spin measurement](https://en.wikipedia.org/wiki/Spin-1/2#Rotations_and_Spinors) probabilities

Let $`𝐮',𝐮`$ be real unit preparation and measurement axes, with $`\cos\theta:=𝐮\cdot𝐮'`$.

Prepare $`\psi=\Pi_+'\psi\ne0`$:

```math
\begin{gathered}
\Pi_+':=\frac12(1+𝐮'\cdot𝝈),\qquad
\Pi_\pm:=\frac12(1\pm𝐮\cdot𝝈),\\[4pt]
\Pi_+'\Pi_\pm\Pi_+'=\frac{1\pm\cos\theta}{2}\Pi_+'.
\end{gathered}
```

Using the squared norm, $`\lVert Z\rVert^2=\mathrm{sc}(ZZ^{\mathsf H})`$, the denominator is

```math
\lVert\psi\rVert^2=\mathrm{sc}(\psi\psi^{\mathsf H})>0.
```

For the [prepared state](https://en.wikipedia.org/wiki/Spin-1/2#Bloch_Representation), the projected norm is

```math
\begin{aligned}
\lVert\Pi_\pm\psi\rVert^2
&=\mathrm{sc}\bigl((\Pi_\pm\psi)(\Pi_\pm\psi)^{\mathsf H}\bigr)\\
&=\mathrm{sc}(\Pi_\pm\psi\psi^{\mathsf H})\\
&=\mathrm{sc}(\Pi_+'\Pi_\pm\Pi_+'\psi\psi^{\mathsf H})\\
&=\frac{1\pm\cos\theta}{2}\mathrm{sc}(\psi\psi^{\mathsf H})
=\frac{1\pm\cos\theta}{2}\lVert\psi\rVert^2.
\end{aligned}
```

The [Born probabilities](https://en.wikipedia.org/wiki/Born_rule) for [spin](https://en.wikipedia.org/wiki/Spin-1/2) $`\pm\tfrac12`$ are

```math
\begin{gathered}
\mathrm{prob}(\pm)=\frac{\lVert\Pi_\pm\psi\rVert^2}{\lVert\psi\rVert^2}
=\frac{1\pm\cos\theta}{2},\\[6pt]
\boxed{\mathrm{prob}(+)=\cos^2\frac\theta2,\quad \mathrm{prob}(-)=\sin^2\frac\theta2.}
\end{gathered}
```

Aligned axes give certain spin up; perpendicular axes give equal odds. Reversing the prepared spin swaps the probabilities.
