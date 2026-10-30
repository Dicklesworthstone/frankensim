# Material-driven spatial electrothermal simulation

`frankensim electrothermal` loads existing canonical material-card packs and
runs the shared spatial current-balance, Joule-heating and heat-conduction
owners. A changed material pack changes both electrical and thermal response;
there are no copper/aluminum/other-material branches in the solver.

```sh
cargo run -p fs-cli --bin frankensim -- --json electrothermal scene.json > result.json
```

Use the repository's required remote execution policy for Cargo work. This is
an experimental file-driven workflow, like `cooling-cht`, not a ledger-backed
`.fsim` solve. It does not modify the input scene or material packs.

## Material input

Each `pack` is a canonical **FSMCDPK** (`NormalizedMaterialCardPack`) artifact,
the same material-card envelope accepted by `solve --materials`. Paths are
relative to the scene. A plain **FSMATPK** properties pack is not a named
material card and is deliberately refused. The existing `fs-matdb` adapter
`NormalizedMaterialCardPack::new(MaterialStateId, NormalizedPack)` binds an
explicit chemistry/phase/process/revision-0 identity to a compiled claims pack;
`to_bytes()` writes its runtime envelope. No label or filename infers this
identity, and unrelated source conditions must not be silently merged.

A conducting material needs positive `electrical_resistivity` in ohm metres.
Heat transport uses either positive scalar `thermal_conductivity` in W/(m K)
or three explicitly bound principal conductivity properties, described below.
An explicitly insulating material needs only the thermal properties. Only
consumed properties are demanded, not mechanical, optical or acoustic data.
Each conductivity preserves the source's dimension-only descriptor or its
static thermal-conductivity kind. Electrical admission is the existing
`resolve_conductor_state` contract.
Conflicting keys, surviving claim-selection ties and unsupported states refuse.
There is no fallback source, implicit averaging, or inferred temperature slope.

The caller supplies every non-temperature condition in `state`. Values use
the existing quantity parser (`"101325 Pa"`, for example). `kind` is explicit:
`legacy` for an untyped legacy source coordinate, `dimensional` for a typed
unit-only coordinate, or a supported semantic kind such as `pressure`.
These modes are not interchangeable. The separate `temperature_axis` is
updated at each cell's actual temperature; `temperature_kind` must be
`legacy`, `dimensional` or `absolute-temperature`. Temperature is always
absolute kelvin, never a temperature difference or a guessed offset.

## Principal material directions

A material can declare three source property roles and their orientation in
the Cartesian heat grid. This optional object belongs inside a `materials`
entry:

```json
"thermal_conductivity": {
  "principal_properties": [
    "conductivity_along_fibre",
    "conductivity_across_fibre",
    "conductivity_through_stack"
  ],
  "principal_axes": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
}
```

Names are caller-selected keys in the **same material card**. Every key must
carry a positive thermal conductivity response with valid source support.
No canonical scalar conductivity is required for this mode. No missing
direction is averaged or copied from another direction. All three responses
are re-resolved at the same actual cell temperature and complete pressure,
composition and other state coordinates as the electrical response.

Rows of `principal_axes` give each principal direction in grid coordinates.
For example, `[[0,0,-1],[0,1,0],[1,0,0]]` maps the first material direction
to grid z, the second to y, and the third to x. Changing this scene data
changes actual heat flux. The admitted frames are exact signed permutations
of the Cartesian unit axes; off-grid rotations are refused because the
native stencil does not represent their off-diagonal cross-gradient terms.
This does not approximate a rotated tensor by retaining only its diagonal.

The existing spatial heat owner receives `[kx, ky, kz]` per cell. It keeps
its original face-normal half-cell resistances, thermal interfaces,
convection and radiation boundary exchange. The same source path operates
in steady and transient solves, including orthotropic solid regions beside
a prescribed isotropic fluid. The native fluid reference remains isotropic;
declaring this directional role on the fluid itself refuses explicitly.

Omitting the object retains the existing scalar isotropic material contract.

## Example scene structure

This is illustrative geometry and loading, **not calibrated material data**.
Supply `conductor.fsmcdpk` with both required properties and source support over
the computed temperature range. The scene alone does not fabricate that pack.

```json
{
  "schema": "frankensim.electrothermal.v1",
  "cells": [3, 1, 1],
  "cell_size_m": 0.01,
  "materials": [
    {
      "pack": "conductor.fsmcdpk",
      "conducting": true,
      "temperature_axis": "T",
      "temperature_kind": "legacy",
      "state": [
        {"name": "pressure", "value": "101325 Pa", "kind": "legacy"}
      ]
    }
  ],
  "material_indices": [0, 0, 0],
  "initial_temperature_k": 300.0,
  "thermal_faces": [
    {"temperature_k": 300.0}, {"temperature_k": 300.0},
    "adiabatic", "adiabatic", "adiabatic", "adiabatic"
  ],
  "electrodes": [
    {"lo": [0, 0, 0], "hi": [1, 1, 1], "voltage_v": 0.0},
    {"lo": [2, 0, 0], "hi": [3, 1, 1], "current_a": 1.0}
  ],
  "numerics": {
    "max_passes": 100,
    "max_thermal_passes": 100,
    "temperature_tolerance_k": 1e-7,
    "power_tolerance_w": 1e-9,
    "relative_tolerance": 1e-9,
    "conductivity_relative_tolerance": 1e-8,
    "energy_tolerance": 1e-12,
    "energy_max_iterations": 20000
  }
}
```

Assignments are in `(z*ny+y)*nx+x` order; every cell is solid and must name a
material index. Each declared material must be used. Initial temperature may
instead be a complete per-cell array. There are at most 65,536 cells and 64
materials, with a 4 MiB scene/individual-pack cap and a 64 MiB aggregate pack
cap. Field geometry and assignments remain fixed during the solve.

All six thermal faces are explicit, ordered x-, x+, y-, y+, z-, z+.
In addition to the rules shown, a face can declare `{"heat_flux_w_m2": 10}`
(positive into the domain), or `{"h_w_m2_k": 5, "ambient_k": 300}`.
Optional `power_w` is an exactly cell-sized signed source array in watts.
Optional `contacts` contains unique pairs such as
`{"materials": [0, 1], "resistance_m2_k_w": 0.001}`. Contact resistance
remains a separate interface input, not a bulk-material property.

Electrode boxes use inclusive `lo`, exclusive `hi`. Each chooses exactly one
signed `current_a` or `voltage_v`; at least one voltage fixes the potential
reference. The existing owner rejects overlaps, touching electrodes and
unsupported current paths. Electrical electrodes are not automatically
isothermal: the thermal boundary rules remain independent.

## What is coupled

For each proposed temperature field, the common material owners re-query the
same immutable cards and complete state coordinates. Resistivity enters the
existing spatial electrical solve. Its per-cell Joule heat is added to the
authored heat sources. Conductivity is then re-resolved at each cell inside a
bounded nonlinear heat iteration using `EnergyConfig::solve_with_conductivity`.
This field ingress does not relabel voxels or create one artificial material
per cell; all interface and geometric identities stay intact.

The inner solve requires both temperature and relative conductivity agreement.
The outer solve requires temperature and per-cell Joule-source agreement. Both
re-query the proposed final state before publication, so leaving a material
source's support is a refusal, not a constant-endpoint extension. Exhausted
budgets also refuse. There is no stability or convergence theorem for arbitrary
input laws, and no silent relaxation or fallback material is introduced.

JSON output contains temperatures, potentials, Joule heat, solved and final
source-resolved conductivities, electrode currents/voltages, power balance,
coupling defects, work counts, pack/card identities and final per-cell source
state identities. The solved conductivity field is reported separately from
the re-evaluated final field: they agree to the declared tolerance, not by
pretending a finite iteration is exact. The complete fields are
`solved_conductivity_grid_w_m_k` and `resolved_conductivity_grid_w_m_k`,
each an array of per-cell `[kx, ky, kz]` triples. Legacy scalar conductivity
entries are `null` where the three components differ. Material records retain
the principal property names and axes alongside pack/card identities.
Unconnected electrical potentials are `null`, not zero. Failed input or physics admission prints only a diagnostic
on stderr and returns nonzero; no partial result is published on stdout.

## Verification and limits

Four `fs-lbm` field-ingress regressions cover a spatial conductivity profile,
shared contact behavior, fluid eddy composition and refusals. Five command
regressions exercise actual file reads, canonical pack decoding, both material
owners, spatial solvers, source-domain refusal, coordinate semantics and replay.
The synthetic three-cell tests independently predict 303.375 K and 300.84375 K
for data-only constant-material substitutions. With both rho(T) and k(T)
varying, an independent resistor-network reduction predicts a centre
303.0302313928287 K and total Joule heat 4.24286287444763 W.

Those independent Python network calculations were executed. The Rust tests
were authored but **not executed in the implementation environment**, which
lacks Cargo/rustc/rustfmt. Run the focused selectors under the repository's
normal execution policy:

```sh
cargo test -p fs-lbm --lib conjugate::energy::material_field::tests
cargo test -p fs-cli --bin frankensim electrothermal_command::tests
```

Five additional directional regressions exercise source-only principal
properties, material orientation, independently varying temperature laws,
native transverse face fluxes, backward-Euler storage and admission failures.
The independent three-cell reference gives a centre rise of `3.375/kx` K
for a 4.5 W source between 300 K end boundaries. With adiabatic boundaries,
mass 2 kg, Cp 3 J/(kg K) and one 0.1 s step, changing only the orientation
from kx=1 to kx=9 W/(m K) changes the centre rise from 0.072619047619 to
0.059482758621 K while stored energy remains 0.45 J. Those independent
arithmetic checks were executed and the changed Rust parsed with rustfmt;
the new Rust tests could not be compiled or executed in the full-disk
implementation environment.

The command uses scalar DC electrical response and isotropic or grid-aligned
orthotropic thermal response on fixed voxels. Optional transient evolution,
prescribed fluid enthalpy transport, source-backed electrical and thermal
interfaces, and exterior gray radiation use their existing owners. It does
not model melting, changing phase topology, expansion, off-grid thermal
tensors, anisotropic electrical transport, skin effects or inductance.
Nominal source values have no propagated uncertainty; source support is
checked at discrete states, without a certified continuum trajectory or
experimental qualification of a real material, joint or device.
