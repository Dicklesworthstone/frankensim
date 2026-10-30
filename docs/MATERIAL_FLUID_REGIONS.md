# Spatial exterior materials without splitting the structural mesh

`MaterialFluidFrame::new_with_regions` consumes `FrameFluidRegion { element,
span, fluid }`. The span is the oriented reference interval `[start,end]` in
an existing frame member, with `0 <= start < end <= 1`. Adjacent regions can
carry different source materials. No structural node, mass, spring or material
section is fabricated at a fluid boundary.

```rust,ignore
use fs_couple::air_path::cylinder::frame::{FrameFluidRegion, MaterialFluidFrame};

let model = MaterialFluidFrame::new_with_regions(
    &frame, retained_modes,
    vec![
        FrameFluidRegion { element: 0, span: [0.0, 0.35], fluid: first_source },
        FrameFluidRegion { element: 0, span: [0.35, 1.0], fluid: second_source },
    ],
    spectrum, frame_budget, fluid_budget, &gate,
)?;
```

The whole-member constructor is the `[0,1]` specialization of the same assembly.
Each region keeps its original fluid object. Radius, mass, elasticity and solid
viscosity still come from the original structural material/geometry owners.
The four Gauss stations move inside each span; their physical weights include
that span's actual length. Cubic displacement traces remain those of the
ORIGINAL beam, including rotational contributions to velocity and nodal moments.

This distinction is physical: half of a guided beam near its root and half near
its tip do not contribute equal added modal inertia or fluid loss. Multiplying
a whole-member damping coefficient by immersed length would lose that effect.
The regional model retains the full added-mass matrix and all causal memory,
not independent retuned frequencies or a length-weighted average fluid.

Region declaration order determines `members()`, every region's memory block,
and `fluid_member_dissipation_j`. Original element IDs can repeat. Inspect
`span()` and `length_m()` on each loaded entry; do not merge its losses merely
because another entry has the same element ID. Existing prepared, liquid-source,
finite-solid-heat and thermoelastic consumers use these same operators and
ordered histories. This does not add a separate integrator or thermal state.

The complete state count is `2*m + 4*regions*spectral_intervals`. All regional
states count before source evaluation or allocation. Overlap checking sorts
indices, not material objects, in O(regions log regions). Shared endpoints are
allowed; overlaps, reversed/nonfinite spans and spans too narrow to represent
four distinct interior stations refuse. Gaps explicitly omit exterior loading;
they do not infer vacuum, air, insulation or another fluid phase.

## Physical boundary

These are frozen, prescribed regions for the existing local 2-D Newtonian
cylinder approximation. Four-point integration is exact for cubic-trace products
in exact arithmetic, not for omitted 3-D flow or arbitrary free-surface effects.
Cross-region/inter-member flow, menisci, interface waves, buoyancy, axial drag,
radially partial cross sections, cavitation, changing immersion/history transfer
and acoustic radiation remain outside this model. A fluid-interface application
must justify that local approximation independently. Moving a boundary requires
a new physical history-transfer model, not reconstruction during a running step.

## Focused qualification

```sh
rch exec -- cargo test -p fs-couple --test frame_fluid_regions
rch exec -- cargo test -p fs-couple --lib air_path::cylinder::frame
```

Six public-API regression groups cover whole-span parity, analytic partial-mass
integrals, heterogeneous cross-mode impedance, partition-invariant native
motion, prepared/dense agreement and pre-source refusal/cancellation. Native
Rust compilation, formatting and test execution are UNVERIFIED in the editing
environment: Rust/Cargo/DSR/RCH are absent and network installation is unavailable.

Independent NumPy/SciPy reference equations were executed. Six one-/four-mode
spatial cases agree with analytic polynomial inertia and impedance elimination
to about 1.1e-14 relative/scaled error. Splitting a uniform region leaves the
100-step structural reference trajectory unchanged within 2.9e-18. An 80-step
heterogeneous case closes its separate material/fluid/work account within
5.8e-19 J. These synthetic equation checks do not compile or execute Rust and
are not measured-material validation, convergence proof or performance evidence.
