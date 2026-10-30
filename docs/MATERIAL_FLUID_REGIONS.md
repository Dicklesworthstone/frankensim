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

## Source geometry and heterogeneous material owners

`FrameFluidRegion::split_at_reference_plane` derives up to two spans from the
original `Frame2d::reference_member_segment`, not a caller-authored duplicate
segment. Supply an explicit unit normal and metre offset in the frame's own
planar coordinates. Both side materials are explicit; there is no air/water
preset, phase selector, or default material. Reversing the member preserves
physical side ownership and reverses the oriented coordinate intervals.

```rust,ignore
use fs_material::fluid_transport::FluidTransport;

let negative: &dyn FluidTransport = &resolved_liquid;
let positive: &dyn FluidTransport = &resolved_gas;
let regions = FrameFluidRegion::split_at_reference_plane(
    &frame, element, [0.0, 1.0], reference_level_m,
    negative, positive, &gate,
)?;
let model = MaterialFluidFrame::new_with_regions(
    &frame, retained_modes, regions,
    spectrum, frame_budget, fluid_budget, &gate,
)?;
```

`FluidTransport for &T` delegates to the original source without making a
transport-only substitute. This supports concrete borrows and heterogeneous
`&dyn FluidTransport` lists in all existing generic consumers. The original
resolved material objects remain alive and immutably borrowed; their admission
errors are not suppressed. `ThermalFluidTransport for &T` likewise delegates
the same temperature and known/unknown expansion rather than assuming 1/T for
liquids. Dynamic dispatch occurs during material preparation, not per spectral
memory coordinate or native mechanical timestep.

The plane helper returns regions in original member traversal order. Endpoint
contact creates no zero-length region; a member wholly on the plane refuses
as ambiguous. Nonfinite/non-unit planes, unsupported reference geometry and
unrepresentable crossings refuse without changing the frame. The structural
reference segment query is read-only and never uses deformed coordinates.

## Physical boundary

These are frozen, prescribed regions for the existing local 2-D Newtonian
cylinder approximation. Four-point integration is exact for cubic-trace products
in exact arithmetic, not for omitted 3-D flow or arbitrary free-surface effects.
Cross-region/inter-member flow, menisci, interface waves, buoyancy, axial drag,
radially partial cross sections, cavitation, changing immersion/history transfer
and acoustic radiation remain outside this model. A fluid-interface application
must justify that local approximation independently. Moving a boundary requires
a new physical history-transfer model, not reconstruction during a running step.
The plane clips a REFERENCE CENTERLINE, not the circular exterior surface; it
is not a finite-radius free-surface intersection or an interface equilibrium.

## Focused qualification

```sh
rch exec -- cargo test -p fs-couple --test frame_fluid_regions
rch exec -- cargo test -p fs-couple --lib air_path::cylinder::frame
rch exec -- cargo test -p fs-solid --lib transverse_trace_tests
rch exec -- cargo test -p fs-material --test borrowed_fluid_transport
```

Six public-API regression groups cover whole-span parity, analytic partial-mass
integrals, heterogeneous cross-mode impedance, partition-invariant native
motion, prepared/dense agreement and pre-source refusal/cancellation. Four
plane regressions additionally cover orientation/translation, endpoint and
coincident cases, cancellation/refusal, and actual prepared motion with original
resolved-liquid/gas source borrows. Three material tests cover shared numerical
and thermal delegation/error preservation; one structural test covers original
endpoint orientation and reference-vs-current geometry. Existing five trace
tests remain unchanged. Native Rust compilation, formatting and test execution
are UNVERIFIED in the editing environment: Rust/Cargo/DSR/RCH are absent and
network installation is unavailable.

Independent NumPy/SciPy reference equations were executed. Six one-/four-mode
spatial cases agree with analytic polynomial inertia and impedance elimination
to about 1.1e-14 relative/scaled error. Splitting a uniform region leaves the
100-step structural reference trajectory unchanged within 2.9e-18. An 80-step
heterogeneous case closes its separate material/fluid/work account within
5.8e-19 J. A separate 90-digit Decimal comparison of 256 manufactured plane
crossings gives maximum absolute coordinate error 1.17e-15, reversal error
2.23e-16 and translated-coordinate discrepancy 4.0e-15. These synthetic equation
checks do not compile or execute Rust and are not measured-material validation,
convergence proof or performance evidence.
