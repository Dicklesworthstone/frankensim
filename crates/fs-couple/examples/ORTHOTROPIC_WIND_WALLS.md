# Directional material walls in the native wind renderer

The shared `DgMaterialWallInput` accepts original isotropic or principal-aligned
orthotropic sources. Both `frankensim-flue-dg-performance-v1` and
`frankensim-material-modal-dg-valve-v1` use this owner, including their physical
`dg_wall regions` assignments. There is no flute-specific or reed-specific
elasticity implementation.

## Source declaration

After the DG bore declaration, a uniform directional wall is written as:

```text
dg_wall orthotropic PACK T OMEGA H AXIAL FMIN FMAX HOOP_AXIS AXIAL_AXIS LOSS
```

`T` is kelvin, `OMEGA` is rad/s, `H` is metres and the applicability band is hertz.
`AXIAL` is `stress-free` (zero incremental axial stress) or `strain-fixed` (zero
incremental axial strain). `HOOP_AXIS` and `AXIAL_AXIS` are distinct one-based
principal material directions in `1..=3`. These directions follow the **local
cylindrical texture**, not one fixed Cartesian grain frame.

`LOSS` may be `elastic`, `proportional-kelvin-voigt TIME_PROPERTY`, or the
independent normal-viscosity declaration below. The proportional form
explicitly declares the normal-plane law `C_viscous = tau * C_elastic`, with a
nonnegative source property `tau` in seconds. A material name never supplies
this property or substitutes an empirical decay coefficient. General anisotropic
Maxwell spectra and full 3-D viscous condensation are not supported by this
declaration; they refuse rather than acquire a scalar approximation. The original
isotropic elastic, Kelvin--Voigt and Maxwell declarations remain unchanged.

For the small synthetic tube in `flue-dg.performance`, an illustrative declaration
is:

```text
dg_wall orthotropic wall.fsmcdpk 293.15 1000 0.00002 strain-fixed 50 2000 1 2 proportional-kelvin-voigt normal_plane_viscosity_time
```

This requires a supplied source pack; it is **not** a bundled wood/material model.
The shared resolver asks that pack for `density`, `young_modulus_1`,
`young_modulus_2`, `poisson_ratio_12`, and `normal_plane_viscosity_time`, at the
same declared `(T, omega)` point. Density is kg/m3, moduli are pascals, the Poisson
ratio is dimensionless, and the viscosity time is seconds. Missing, ambiguous,
wrong-dimension, out-of-domain or unsupported frequency-dependent data refuse.
The unused radial and shear coefficients are not invented or required.

## Physical reduction

The major Poisson ratio `nu_ij` denotes transverse strain in direction `j` under
uniaxial stress in direction `i`, with `i < j`. The reciprocal normal compliance
is `S_ij = S_ji = -nu_ij/E_i`, not a symmetric assignment of Poisson ratios.
For hoop direction `h` and axial direction `a`:

```text
S_hh = 1/E_h
S_aa = 1/E_a
E_radial_response = E_h                         (stress-free)
E_radial_response = 1/(S_hh - S_ha^2/S_aa)       (strain-fixed)
surface_mass = density * H
radial_stiffness_per_area = E_radial_response * H / radius^2
radial_resistance = tau * radial_stiffness_per_area   (proportional loss only)
```

The represented normal compliance block must be positive definite; anisotropic
Poisson ratios are not incorrectly clamped to the isotropic interval. The source
card, original query receipts, direction mapping and loss declaration remain
attached to `LocalMaterialWall`. This is a geometry/constitutive reduction, not
an effective isotropic material card.

The existing DG gas/wall operator owns wall inertia, pressure feedback, stored
energy and loss. Its gas field and viscothermal histories remain separate and
unchanged. A direction or source substitution changes the actual operator and
pressure output, without retuning or changing the PCM scale.

## Independent normal-plane viscosity

For a material with a supplied two-dimensional normal-plane law
`sigma = C * strain + D * strain_rate`, select:

```text
dg_wall orthotropic wall.fsmcdpk 293.15 1000 0.00002 stress-free 50 2000 1 2 normal-kelvin-voigt 1 0.8 wall_eta_hh wall_eta_ha wall_eta_aa
```

The three properties are **stiffness-viscosity matrix entries** `D_hh`, `D_ha`
and `D_aa`, in Pa s, in the ordered local hoop/axial frame. They are not
uniaxial Young viscosities. The cross entry can be negative. The whole symmetric
2-by-2 D matrix must be positive semidefinite; exact zero and rank-deficient
laws are admitted, but negative coefficients are never clipped. Reversing the
source axes requires assigning the viscosity properties to that new frame too.

The two numerical limits before the names retain the shared wall source
contract: `0 < MAX_DT_OVER_TAU <= 2` and `0 < MAX_ANGULAR_STEP <= 1`. The exact DG
exponential does not apply a midpoint step to the generated relaxation. Its
complete-state cap, accuracy admission and nonlinear source clock still apply.

Fixed axial strain uses `C_hh + s*D_hh`. Free axial stress requires the **whole
constitutive Schur complement**, not separate condensation of C and D:

```text
M(s) = C_hh + s*D_hh - (C_ha + s*D_ha)^2 / (C_aa + s*D_aa)
     = E_h + eta*s + delta_E*s*tau / (1 + s*tau)
eta = D_hh - D_ha^2 / D_aa
tau = D_aa / C_aa
delta_E = C_aa * (C_ha/C_aa - D_ha/D_aa)^2
```

The native Maxwell state carries this induced axial memory, while the base wall
retains residual Kelvin resistance. Both contribute to the same DG pressure,
energy and dissipation. A zero `D_aa` requires zero `D_ha` and has no axial memory;
proportional coupling has zero induced strength. No spectrum is fitted. Existing
gas viscothermal states remain independent. The source record retains the
original three names (`normal_viscosity_properties()`), pack and property
receipts. The existing sidecar records the realized Maxwell arm when present;
its `loss: maxwell` label identifies that representation, not a new source pack.

This declaration supplies an already two-dimensional normal constitutive law.
It does **not** take an arbitrary 3-D viscous tensor and independently eliminate
radial stress. Generated memory cannot use the live source-pair rebind API:
there is no direct original modulus/time property pair to identify that history.
It is cold-bound at the authored state, like the other native wind inputs.

## Regional sources and Instrument Lab

Each `dg_wall_region X_BEGIN X_END` row accepts the same `orthotropic ...`
sub-declaration. Regions must cover the main bore contiguously and their physical
interfaces must coincide with DG faces. One original pack may supply several
orientations or states; the loader caches its bytes once and resolves each
assignment independently. Isotropic and orthotropic source regions can coexist.
Refinement must preserve those physical interfaces; no element averages two
materials to admit a misaligned boundary.

Render with the existing command:

```sh
music_render wind source.performance output.wav --decimate
```

Use `--decimate` only when the file's pressure-output clock is above 48 kHz.
Material-modal DG sources keep their original clock rules. In Instrument Lab,
open the performance in the native wind panel and attach the exact original pack
at its declared path. The existing self-contained snapshot route handles these
files; no browser material parameters or new server endpoint are needed.

The native sidecar records `orthotropic-principal-local-ring-membrane`, one-based
hoop/axial directions, the original pack/state identities, and the proportional
viscosity property and time when present. Existing isotropic metadata stays
byte-identical. Mixed/directional regional layouts are identified as
`regional-local-ring-membranes`, with each original region described separately.

## Scope and verification

This remains a frozen-state, axisymmetric **thin cylindrical membrane** with
`H/radius <= 0.05`. It omits axial wave inertia, bending between rings, thick-wall
radial stresses, off-axis shear coupling, end-cap traction, changing geometry,
and exterior wall radiation. Source-frame alignment must describe the specimen:
a straight Cartesian grain direction cannot silently become cylindrical texture.
The normal-plane check does not validate an unspecified full 3-D elastic tensor.

This addition does not by itself qualify a wooden clarinet, oboe or flute as a
measured instrument. Complete specimen data, applicable geometry, spatial/time
refinement and recording comparisons remain separate requirements.

Focused regressions are in `tests/orthotropic_duct_wall.rs`,
`tests/normal_viscous_duct_wall.rs` and `tests/orthotropic_wind_input.rs`: reciprocal compliance, the isotropic limit,
loss, source refusals, gas/solid pressure and energy, both native source families,
regional orientations, exact PCM and attached-source parity. The independent
normal-viscosity cases additionally compare the full two-strain dynamics,
stored energy and dissipation with the native condensed history, then exercise
both wind loaders and regional memory support. Run them with the repository's
pinned build environment:

```sh
cargo test -p fs-couple --test orthotropic_duct_wall --test normal_viscous_duct_wall --test orthotropic_wind_input
```

At initial publication these Rust tests have not been executed in the editing
environment, which lacks a Rust toolchain. Independent numerical matrix-inversion
and full two-strain time-domain checks validate the reduction algebra, not
compilation, native renderer execution or instrument realism.
