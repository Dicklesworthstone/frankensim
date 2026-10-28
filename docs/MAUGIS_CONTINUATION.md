# Maugis absolute energy and noncontact continuation

`fs_contact::normal_patch::maugis::continuation` extends the original immutable
`SourceMaugisSphereContact`; it does not fit another normal law or change its
material sources. The parent envelope, source identity and synthetic authority
remain attached. Existing contacting-branch and work APIs are unchanged.

`evaluate_contact_energy(delta)` adds separate elastic storage and surface
potential, with undeformed separated bodies as the common zero. The surface
potential density is `sigma0*(gap-z0)` within the Dugdale range and zero outside;
`sigma0=W/z0`. Its sum with elastic storage obeys `dU/d(delta)=F`.

The calculation uses the equivalent Abel line load: `E*(delta-x*x/R)` for
`0<x<a`, and `-2*sigma0*sqrt(c*c-x*x)` for `a<x<c`. Integrating `q*q/E*`
gives elastic energy. Reciprocal work against a uniform disk gives the surface
part. The resulting closed expressions use cancellation-resistant annulus
widths and angular differences. Evaluation adds constant work to the original
bounded contact inversion, with no allocations or energy quadrature.

`evaluate_noncontact(delta)` retains deformation from a uniform attractive disk:

```
delta = c*c/(2*R) - 4*sigma0*c/(pi*E*) - z0
F = -pi*sigma0*c*c
gap_at_center = -delta - 2*sigma0*c/E*
```

Only the locally stable axisymmetric branch with nonnegative center gap is
admitted. `noncontact_limits()` returns its interval, or `None` when contact
precedes the stable disk branch. A fold indistinguishable from the cutoff in
floating point refuses rather than silently erasing the distinction. The real
solid-contact area remains zero, including in the existing junction-shear
adapter. Radius/depth/layer limits are checked on the actual disk.

`evaluate_detached(delta)` is a separate zero-force, zero-energy branch beyond
the interaction range. Coexisting states are not silently selected. Smooth
ports retain signed action/reaction and withhold tangents at folds and branch
junctions. These evaluations alone neither evolve history nor assign released
energy to caloric heat. They are not nonaxisymmetric stability certificates,
finite-body models, mixed-mode peeling, plasticity, viscoelastic adhesion or
experimental material validation.

Five integration test groups compare energy gradients and the existing branch
work integral, independent line-load energy quadrature, deformed noncontact
geometry and force, zero-area continuity, and source/domain/port behavior.

Theory: Maugis (1992), DOI 10.1016/0021-9797(92)90285-T; Leroy and Willert
(2021), DOI 10.1002/zamm.201900313; the uniform-disk noncontact extension is
also studied by Shi and Polycarpou (2005), DOI 10.1016/j.jcis.2004.08.069.
The energy formulas here are derived from the stated traction potential and
linear-elastic reciprocal work, not fitted from those papers' jump curves.
