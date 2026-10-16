# Radiation equations selected from the physical components

The exterior piano paths now choose their existing BEM formulation from the
actual closed solid components, not `wavenumber * whole_scene_radius < 0.5`.
This applies consistently to raw `response`, loaded `admittance`, ordinary
`render`, and passive-feedback `render-loaded`. It requires no new user flag.
Native section skins, posed rigid parts and supplied closed OBJ surfaces all
use the same rule. Microphone location does not choose the source equation.

## Why the previous lid cases failed

The 28-panel native-board/posed-lid regression crossed the old switch at about
281.735 Hz. Above that switch, the existing Burton--Miller implementation uses
centroid/equivalent-disc hypersingular terms rather than the triangle weak-
kernel integration of the plain equation. In this thin geometry its discrete
radiation resistance became negative, correctly causing the physical power
gate to refuse. Moving the lid changed a global bounding radius even though
empty space between disjoint solids is not an interior resonance cavity.

The change does not clip that negative resistance, add damping, alter the
geometry, narrow the requested band, or relax a failed test. The two original
rigid-assembly playback/CLI regressions and their thresholds remain unchanged.

## Geometric resonance exclusion

`fs_bem::radiation_policy::GeometryPolicy` reuses the exact-coordinate closed-
component admission. Each component has an enclosing axis-aligned box of side
lengths `Lx, Ly, Lz`. Dirichlet domain monotonicity and the box spectrum give

```text
first interior Dirichlet wavenumber >= pi * sqrt(1/Lx^2 + 1/Ly^2 + 1/Lz^2).
```

These are the fictitious frequencies of the conventional exterior equation,
not the piano's physical string, wood, or exterior acoustic resonances. Each
solid has its own interior spectrum; the minimum component bound controls the
joint problem. Box widths are rounded upward and the positive bound arithmetic
downward. The existing triangle plain-CBIE image is selected below 80% of that
minimum bound. The existing Burton--Miller image is retained above it.

The 20% margin is a fixed spectral-selection margin, not an adjustable power
or fit tolerance. A selected solve is performed once; a refusal does not retry
another equation. Both microphone pressure and the full signed radiation-
impedance matrix come from that same selected source equation. All existing
wavelength, power, geometry, fit, energy and output-publication checks remain.
Failed power/conditioning diagnostics now also name the formulation and report
the actual power interval, resolution and conditioning lower bound.

## Scope of the guarantee

The component bound excludes *continuum fictitious resonances* in the selected
band. It does not certify a discrete collocation matrix, solve conditioning,
quadrature, positive radiation, acoustic fit, or full-band instrument fidelity.
Axis-aligned bounds can become less tight after rotation without becoming
invalid. The input contract still requires disjoint, outward, non-self-
intersecting components; no global collision or cavity-accessibility proof is
added. A full Model D mesh may still fail spatial or passive-fit admission.

The public explicit `fs_bem::helmholtz` formulations remain available unchanged
for controlled comparisons and the existing sphere/resonance oracles. The new
policy composes those operators; it is not a replacement acoustic solver.

References: *An improved form of the hypersingular boundary integral equation
for exterior acoustic problems*, Engineering Analysis with Boundary Elements
34 (2010), 189--195, DOI 10.1016/j.enganabound.2009.10.005; V. Ivrii's PDE
textbook, sections 13.2--13.3, for the box spectrum and min--max principle.
