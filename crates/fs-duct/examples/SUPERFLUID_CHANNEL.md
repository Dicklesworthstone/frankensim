# Superfluid helium-II: pressure waves, second sound and counterflow

```sh
cargo run -p fs-duct --features superfluid --example superfluid_second_sound
```

This native example starts a one-microkelvin temperature standing wave in a
10 cm closed channel using the shared `fs-material::superfluid` owner. It
prints the resulting pressure, both component velocities, heat flux, wave
energy and exported friction heat. The pressure/temperature coupling follows
the material's finite thermal expansion; no output filter creates the wave.

The material reference is **one point at 1.8 K and saturated vapor pressure**,
from Donnelly and Barenghi, *J. Phys. Chem. Ref. Data* 27, 1217 (1998),
[DOI 10.1063/1.556028](https://doi.org/10.1063/1.556028),
[open NIST paper](https://srd.nist.gov/jpcrdreprint/1.556028.pdf).
The original density, entropy, heat capacity, expansion, viscosity and sound
observations remain available from `HeliumSvpReference::raw()`.
Saturation-path heat capacity and expansion are converted to their isobaric
forms using the saturation-pressure derivative. Measured fast sound is also
distinguished from the thermodynamic adiabatic derivative. This reference
predicts second sound at about 19.7455 m/s, compared with the independently
tabulated 19.90 m/s. That approximately 0.8% difference is retained; there is
no tuning to make the independent comparison exact.

## Public physical interfaces

`TwoFluidChannel` retains one admitted `TwoFluidState`, a uniform length and
area, and explicit perturbation limits. Initialize N pressure/temperature
cell samples and N+1 normal/superfluid face velocities. `step` advances one
declared SI interval with component velocities prescribed at both ends and a
cancellation callback. Positive flow points along +x at both ends.

`cell` returns pressure, temperature, density and specific-entropy
perturbations. `face` returns normal and superfluid velocities, barycentric
velocity, counterflow, mass flux, entropy flux and barycentric-relative heat
flux. The quantities share the same accepted state:

```text
v = f_n v_n + f_s v_s
w = v_n - v_s
j = rho v
entropy flux = rho s v_n
counterflow heat flux = rho_s s T0 w
```

Both sound modes evolve, including their coupling through expansion. The
existing lossless TMM/DG acoustic interfaces also accept the explicitly
selected `first_sound_projection()` of a **zero-expansion** material. That
projection retains its original source but discards the entropy mode;
nonzero expansion refuses. Neither the two-fluid owner nor that projection
implements ordinary Newtonian/Fourier transport or thermoviscous wall loss.

## Friction and thermal handoff

The default example has no assigned mutual friction. To exercise an explicitly
authored linear coefficient (this number is **not** a helium measurement):

```sh
cargo run -p fs-duct --features superfluid --example superfluid_second_sound -- 128 1000
```

The kernel also accepts the declared Gorter–Mellink law. Both use exact local
relative-velocity relaxation, conserve total momentum and export the removed
kinetic energy exactly once. `last_mutual_friction_heat_j()` supplies per-cell
joules for an external thermal owner; interior-face loss is partitioned
equally to its two adjacent cells. A cancelled/refused advance preserves the
previous fields, clock and heat account.

The first-order wave model keeps its reference properties fixed. Its positive
quadratic storage is **wave availability**, not absolute thermal energy. The
boundary wave-work rate is `area * (p' v + T' q/T0)`; it differs from the
physical counterflow heat rate `area * q`. Step reports keep both quantities
separate. Frictional heat is exported rather than added again to the linear
temperature perturbation. Finite-temperature feedback requires a separate
caloric evolution and a new admitted reference state.

## Numerical and physical scope

The staggered spatial incidence and implicit midpoint step preserve quadratic
wave storage up to the returned residual. A 2x2 singular-mode transform and
the existing `fs-la::tridiagonal` factorization give O(N) work and storage.
Successful steps reuse allocated buffers and same-step factors. A Courant
limit of one is an explicit temporal-resolution policy; it is not an error
estimate. Refine space and time for a specific output quantity.

This is a constant-reference, small-signal, uniform-area continuum channel
with free normal bulk motion. It does not resolve normal viscous boundary
layers, longitudinal/second viscosities, nonlinear convection, quantum vortex
formation or reconnection, critical-velocity onset, lambda-transition dynamics,
thin-film third sound, porous fourth sound, Kapitza resistance, boiling or
arbitrary-pressure helium thermodynamics. The retained normal shear viscosity
does not silently activate a classical no-slip wall law.

Focused native checks:

```sh
cargo test -p fs-material --test superfluid
cargo test -p fs-material --lib superfluid::helium
cargo test -p fs-duct --test superfluid_projection
cargo test -p fs-duct --features superfluid --test superfluid_channel
```

The channel tests check independent first/second-sound speeds and spatial/time
refinement, finite-expansion polarization, mass/entropy port balances, passive
friction with one heat account, deterministic retry and atomic cancellation.
The source comparison checks a measured reference point; it does not certify
the complete simulation against an experiment.
