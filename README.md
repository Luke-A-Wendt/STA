# Space–Time Algebra

**Complex scalars and three-vectors connect spacetime, Maxwell's equations, and Dirac theory.**

[Read the paper](./sta_notes.pdf) · [LaTeX source](./sta_notes.tex)

Natural units $\hbar=c=1$, rationalized electromagnetic units, and signature $(+,-,-,-)$.

## One complex paravector

$$
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
$$

Eight real components, written as a complex scalar and a complex vector:
$Z=S+\mathbf V\cdot\boldsymbol{\sigma}$. The product carries both the dot and cross products:

$$
(\mathbf A\cdot\boldsymbol{\sigma})(\mathbf B\cdot\boldsymbol{\sigma})
=\mathbf A\cdot\mathbf B+\mathrm{i}(\mathbf A\times\mathbf B)\cdot\boldsymbol{\sigma}.
$$

## Six core operations

Write $Z=S+\mathbf V\cdot\boldsymbol{\sigma}$, where $S\in\mathbb C$ is a complex scalar and $\mathbf V\in\mathbb C^3$ is a complex vector.

On coefficients, $^*$ means ordinary complex conjugation. On a paravector, it also reverses the vector sign; $^{\mathsf H}$ is Hermitian conjugation.

$$
\begin{aligned}
\det Z&=S^2-\mathbf V\cdot\mathbf V,
&\operatorname{tr}Z&=2S,\\
\lVert Z\rVert^2&=S^*S+\mathbf V^*\cdot\mathbf V,
&Z^{-1}&=\frac{S-\mathbf V\cdot\boldsymbol{\sigma}}{\det Z}\quad(\det Z\ne0),\\
Z^*&=S^*-\mathbf V^*\cdot\boldsymbol{\sigma},
&Z^{\mathsf H}&=S^*+\mathbf V^*\cdot\boldsymbol{\sigma}.
\end{aligned}
$$

## Rotations, boosts, and projections

For a real unit axis $\mathbf u$, a rotation angle $\theta$, and a boost rapidity $\eta$:

$$
\begin{aligned}
R&=e^{-\mathrm{i}\theta\mathbf u\cdot\boldsymbol{\sigma}/2}
=\cos\frac\theta2-\mathrm{i}\mathbf u\cdot\boldsymbol{\sigma}\sin\frac\theta2,
&Z'&=R^{\mathsf H}ZR,\\
L&=e^{-\eta\mathbf u\cdot\boldsymbol{\sigma}/2}
=\cosh\frac\eta2-\mathbf u\cdot\boldsymbol{\sigma}\sinh\frac\eta2,
&Z'&=L^{\mathsf H}ZL,\\
\Pi_\pm&=\frac12(1\pm\mathbf u\cdot\boldsymbol{\sigma}),
&\Pi_\pm^2&=\Pi_\pm,\qquad\Pi_+\Pi_-=0.
\end{aligned}
$$

These are the paper's frame-component transformations; the boost velocity is $\tanh\eta\,\mathbf u$. Rotations are also unit quaternions with $q_k=-\mathrm{i}\sigma_k$. Their vector sandwich recovers Rodrigues' formula in the same convention:

$$
\mathbf r'=\mathbf r\cos\theta-(\mathbf u\times\mathbf r)\sin\theta
+(\mathbf u\cdot\mathbf r)(1-\cos\theta)\mathbf u.
$$

## Spacetime from the determinant

For a future-directed massive particle, proper time and four-momentum follow from a real paravector. Here $\mathrm d\mathbf r^2:=\mathrm d\mathbf r\cdot\mathrm d\mathbf r$:

$$
\begin{aligned}
\mathrm dX&=\mathrm dt+\mathrm d\mathbf r\cdot\boldsymbol{\sigma},
&\mathrm ds^2&=\det(\mathrm dX)=\mathrm dt^2-\mathrm d\mathbf r^2,\\
U&=\frac{\mathrm dX}{\mathrm ds}=\gamma(1+\mathbf v\cdot\boldsymbol{\sigma}),
&\gamma&=(1-\mathbf v^2)^{-1/2},\\
P&=mU=\underbrace{E}_{\text{energy}}+
\underbrace{\mathbf p\cdot\boldsymbol{\sigma}}_{\text{momentum}},
&\det P&=E^2-\mathbf p^2=m^2.
\end{aligned}
$$

## Maxwell in one equation

Combine the electric and magnetic fields and use the paravector derivative. Hats mark named energy, momentum, and Dirac operators; the derivative symbols $\partial$ and $\Box$ remain unhatted:

$$
F=(\mathbf E+\mathrm{i}\mathbf B)\cdot\boldsymbol{\sigma},
\qquad \partial=\partial_t+\partial_{\mathbf r}\cdot\boldsymbol{\sigma},
\qquad \boxed{\partial F=\rho-\mathbf J\cdot\boldsymbol{\sigma}.}
$$

Equating its four components gives Maxwell's equations:

$$
\begin{aligned}
\text{Real scalar:}\quad &\partial_{\mathbf r}\cdot\mathbf E=\rho,\\
\text{Imaginary scalar:}\quad &\partial_{\mathbf r}\cdot\mathbf B=0,\\
\text{Real vector:}\quad &\partial_t\mathbf E-\partial_{\mathbf r}\times\mathbf B=-\mathbf J,\\
\text{Imaginary vector:}\quad &\partial_t\mathbf B+\partial_{\mathbf r}\times\mathbf E=0.
\end{aligned}
$$

Applying $\partial^*=\partial_t-\partial_{\mathbf r}\cdot\boldsymbol{\sigma}$ gives the d'Alembertian and the sourced wave equations:

$$
\begin{aligned}
\partial^*\partial&=\Box=\partial_t^2-\partial_{\mathbf r}^2,\\
\partial^*(\partial F)&=\Box F=\partial^*(\rho-\mathbf J\cdot\boldsymbol{\sigma}).
\end{aligned}
$$

Since $\Box F=(\Box\mathbf E+\mathrm{i}\Box\mathbf B)\cdot\boldsymbol{\sigma}$, its components give

$$
\begin{aligned}
\text{Real scalar:}\quad &0=\partial_t\rho+\partial_{\mathbf r}\cdot\mathbf J
&&\text{(charge conservation)},\\
\text{Imaginary scalar:}\quad &0=0,\\
\text{Real vector:}\quad &\Box\mathbf E=-\partial_t\mathbf J-\partial_{\mathbf r}\rho,\\
\text{Imaginary vector:}\quad &\Box\mathbf B=\partial_{\mathbf r}\times\mathbf J.
\end{aligned}
$$

In vacuum, $\Box F=0$.

## Dirac as a first-order wave equation

Package two paravector amplitudes into $\boldsymbol\Psi=(\Psi_1,\Psi_2)^{\mathsf T}$, each with two independent complex components. With $\mathbf I$ the block identity, define

$$
\mathbf W(Z)=\begin{pmatrix}0&Z\\Z^{*\mathsf H}&0\end{pmatrix},
\qquad \mathbf W(Z)^2=\det(Z)\mathbf I.
$$

Define the free Dirac operator and its equation:

$$
\hat D(m):=\mathbf W(\mathrm{i}\partial)-m\mathbf I,
\qquad \boxed{\hat D(m)\boldsymbol\Psi=0.}
$$

Since $\mathbf W(\mathrm{i}\partial)^2=-\partial^*\partial\,\mathbf I=-\Box\mathbf I$, the opposite mass signs factor the Klein–Gordon operator:

$$
\begin{aligned}
\hat D(m)\hat D(-m)&=\hat D(-m)\hat D(m)\\
&=\mathbf W(\mathrm{i}\partial)^2-m^2\mathbf I
=-(\Box+m^2)\mathbf I.
\end{aligned}
$$

Applying this product to a free Dirac solution gives

$$
\hat D(-m)\hat D(m)\boldsymbol\Psi=-(\Box+m^2)\boldsymbol\Psi=0
\quad\Longrightarrow\quad\boxed{(\Box+m^2)\boldsymbol\Psi=0.}
$$

**With an electromagnetic potential**, $\Phi=V+\mathbf A\cdot\boldsymbol{\sigma}$ and charge $q$, the equation is simply

$$
\boxed{\bigl(\mathbf W(\mathrm{i}\partial-q\Phi^*)-m\mathbf I\bigr)\boldsymbol\Psi=0.}
$$

## Low-energy Dirac recovers Pauli

Remove the rest phase by writing $\Psi_k=e^{-\mathrm{i}mt}\Psi_k'$. On these envelopes, hats mark the residual-energy and kinetic-momentum operators

$$
\hat E=\mathrm{i}\partial_t-qV,\qquad
\hat{\mathbf p}=-\mathrm{i}\partial_{\mathbf r}-q\mathbf A.
$$

Eliminating one Dirac entry gives the exact residual-energy recursion, including the electric-field commutator:

$$
\begin{aligned}
\hat E\Psi_2'&=\frac{(\hat{\mathbf p}\cdot\boldsymbol{\sigma})^2-\hat E^2
-[\hat{\mathbf p}\cdot\boldsymbol{\sigma},\hat E]}{2m}\Psi_2',\\
[\hat{\mathbf p}\cdot\boldsymbol{\sigma},\hat E]&=-\mathrm{i}q\mathbf E\cdot\boldsymbol{\sigma},
\qquad [A,B]:=AB-BA.
\end{aligned}
$$

To retain both entries and their probability density, set $\phi_1=(\Psi_1'+\Psi_2')/\sqrt2$ and $\phi_2=(\Psi_2'-\Psi_1')/\sqrt2$. Their exact equations give the small-amplitude recursion

$$
\hat E\phi_1=(\hat{\mathbf p}\cdot\boldsymbol{\sigma})\phi_2,
\qquad
\phi_2=\frac{\hat{\mathbf p}\cdot\boldsymbol{\sigma}}{2m}\phi_1-\frac{\hat E\phi_2}{2m}.
$$

For slow positive-energy envelopes in weak, slowly varying fields, the first iterate gives

$$
\phi_2\simeq\frac{\hat{\mathbf p}\cdot\boldsymbol{\sigma}}{2m}\phi_1,
\qquad
\hat E\phi_1\simeq\frac{(\hat{\mathbf p}\cdot\boldsymbol{\sigma})^2}{2m}\phi_1
=\frac{\hat{\mathbf p}^2-q\mathbf B\cdot\boldsymbol{\sigma}}{2m}\phi_1.
$$

Thus the large two-component amplitude obeys Pauli's equation at leading nonrelativistic order:

$$
\boxed{\mathrm{i}\partial_t\phi_1=
\left[\frac{(-\mathrm{i}\partial_{\mathbf r}-q\mathbf A)^2}{2m}
+qV-\frac{q}{2m}\mathbf B\cdot\boldsymbol{\sigma}\right]\phi_1.}
$$

The magnetic spin coupling follows directly from the squared paravector momentum above.
