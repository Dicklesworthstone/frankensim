# File-driven coupled cooling network

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network examples/cooling-network/mixed-slab.json
cargo test -p fs-cli --bin frankensim network_command::tests
cargo test -p fs-cli --test cooling_network
```

> **Transient numbers (2026-09-25).** Commit 05db922bf made the transient
> heat capacity row-sum lumped. The consistent P1 mass let a cold radiative sink
> raise a sampled peak, which violates the comparison principle. On the 12-tet
> teaching meshes, sampled transient peaks moved by up to about 2 K. On a
> refined plate the two rules agree within 0.6% and converge together. Steady
> results are unchanged. Several examples' temperature limits were re-tuned to
> keep their role; see 45af6da3d. TRANSIENT_COOLING, ADAPTIVE_COOLING,
> NONLINEAR_TRANSIENT_COOLING, REPEATED_COOLING, TRANSIENT_DESIGN,
> ENCLOSURE_RADIATION, COMPONENT_POWER_ALLOCATION and MULTI_LIMIT_ALLOCATION now
> state FrankenSim's measured lumped values next to their consistent-mass
> references. Other documents' transient figures (for example
> TIME_CONVERGENCE, ADAPTIVE_MESH_COOLING, RADIATIVE_MESH_STUDIES and
> TRANSIENT_RADIATION) are consistent-mass references until re-measured.

This **experimental binary command** accepts a JSON request containing an actual
linear tetrahedral solid, named exterior cooling faces, a quadratic hydraulic
graph, independent pressure-reservoir temperatures, and explicit numerical
budgets. It runs the existing `LossGraph`, `TransportNetwork`, shared-solid
conjugate solver and optional coupled adjoint. The request is not a `.fsim`
project, and this command does not perform CAD import, create a ledger, or claim
a package/report workflow. Existing `fs_cli::run`/`run_os` library verbs are
unchanged; `cooling-network` is dispatched by the `frankensim` binary.

Every vertex, tetrahedron, surface, coefficient, heat source, hydraulic connection
and source temperature comes from the request. Supply a conforming,
non-overlapping tetrahedral mesh. Admission checks indices, repeated cells, face
incidence and exterior-face ownership; these checks are not a geometric
non-overlap certificate. Materials are declared constant isotropic conductivities,
either uniform or assigned per tetrahedron. Shared vertices imply continuous
solid temperature across material interfaces, not a finite contact resistance.

## Request contract

`schema` is exactly `frankensim.cooling-network.v1`, `units` is exactly `SI`, and
`seed` is a decimal u64 **string**, avoiding JSON floating-point seed rounding.
All quantity keys carry coherent-SI units. Unknown fields and duplicate keys
refuse. Every declared field in the example is required unless stated otherwise.

`budgets` bounds graph sweeps, shared-solid coupling iterations, Krylov iterations
per solve, adjoint interface iterations, and numerical wall time. Wall timing
starts **after** input/mesh admission and applies to the whole numerical invocation;
expiry requests the cancellation gate and prevents publication after draining.
A checkpoint is cooperative, not a hard real-time deadline. The input cap is
16 MiB, with fixed ceilings of 20,000 solid vertices, 100,000 tetrahedra,
4,096 hydraulic nodes and 16,384 branches.

`tolerances` names absolute flow, heat and coupling-temperature gates, a relative
linear/derivative tolerance, and fixed coupling/interface relaxation in `(0,1]`.
The command does not substitute relative allowances for declared absolute heat
or hydraulic tolerances. Excessively tight tolerances may refuse.

Each hydraulic branch obeys `p[from] - p[to] = R Q |Q|`.
`resistance_pa_s2_m6` is in Pa/(m³/s)². Its required `source` identifies the caller's
coefficient declaration, not external validation. `regions` is a stream-wise list
in the declared `from` → `to` direction; reverse flow reverses the thermal order.
An **explicit empty list** declares an adiabatic bypass. No leakage path is inferred.

A boundary's `temperature_k` is optional only when that pressure reservoir does
not supply air. Every actual external supply must declare one, and declarations
on non-supplies refuse after the hydraulic solve. At a reservoir receiving both
incoming branches and external supply, the supplied temperature is the entering
stream's temperature, not an imposed reset of the node mixture.

`solid.surfaces` names triangle vertex triples on the mesh exterior and a uniform
`htc_w_m2_k` per surface. A face cannot belong to two surfaces. Every surface must
belong to exactly one branch, while one branch may contain many surfaces. Wetted
areas are integrated from the mesh. Unlisted exterior faces are insulated only
when `adiabatic_remainder` is true; otherwise every exterior face must be owned.
A zero-flow heat exchanger refuses.

## Heterogeneous materials and component heating

Choose exactly one conductivity representation inside `solid`:

* `conductivity_w_m_k`: one positive scalar, preserving the original input mode.
* `materials` plus `element_materials`: a nonempty table of
  `{"name":"spreader","conductivity_w_m_k":20,"source":"caller declaration"}`
  objects, and exactly one material name for every tetrahedron in input order.

Missing, duplicate and unknown names refuse. Material source strings remain
caller declarations; they do not create matdb receipts. The same per-element
assignment reaches every primal and adjoint solve through `fs-conduction`.

Choose exactly one heat-source representation:

* `source_w_m3`: the original uniform volumetric density, including zero.
* `component_power`: an explicit system total, relative total-power tolerance,
  and nonempty component list. For example:

```json
"component_power": {
  "total_w": 1,
  "relative_tolerance": 1e-12,
  "components": [{"name":"chip","watts":1,"vertices":[4]}]
}
```

A component footprint is a nonempty set of mesh vertices. Component powers must
be nonnegative; repeated vertices within one component and duplicate component
names refuse. Different components may overlap and their sources superpose.
The existing `PowerMap` distributes each component by its footprint's lumped
nodal volume and audits the power injected by the consistent P1 source assembly.
The final solve also checks its independently accumulated source total against
that delivered total.

**This is nodal P1 support, not a sharp cellwise source.** Heating extends over
tetrahedra incident to the selected vertices, including across material
interfaces. The output's `solid_inputs` preserves material declarations and
assignments, component delivered powers and bound volumes. Unstated power
uncertainty stays `null`; numerical total-power agreement is not uncertainty
propagation or experimental evidence.

## Mean and hotspot objectives

The `objective` object requires `gradient` and exactly one selector:

```json
{"mean_wall_region":"first-face","gradient":true}
{"max_wall_region":"first-face","gradient":true}
{"max_solid_temperature":true,"gradient":true}
{"max_vertices":[1,4,7,10],"gradient":true}
```

`mean_wall_region` remains the area-mean objective. `max_wall_region` checks every
vertex of that named exterior surface; `max_solid_temperature` checks every
solid vertex, including non-cooled faces and interior nodes; `max_vertices`
checks the explicit nonempty vertex set. These are exact maxima of the retained
P1 field on the corresponding tetrahedral domain or surface, or of the selected
nodal values, **not bounds on the continuum solution**. A surface mean is never
substituted for a peak.

The `objective` result records kind, value, active vertex and position, exact
maximum tie count, and separation from the next competitor. Optional gradients
use the actual coupled nodal adjoint at the selected hottest vertex. Exact ties
choose the lowest vertex ID and disclose that this is one active-branch
derivative, not a unique differentiable maximum or every directional derivative.
Near-tie stability is not certified. Every design trial reselects the actual
maximum; derivatives guide trials, not passing/failing decisions.

`dobjective_dinlet_k` and each wall's `dobjective_dlog_htc` apply to the selected
objective. Legacy `objective_mean_k` and `dmean_*` fields retain their meaning for
means and are `null` for maxima rather than carrying mislabeled peak results.

## Results and target sizing

Successful output is one JSON document containing actual nodal temperatures,
node mixture temperatures, signed branch flows, branch inlet/outlet temperatures,
surface heat rates and references, the requested objective, optional gradients,
and solid/air energy diagnostics. No result is published on a failed producer,
exhausted iteration budget or expired wall budget. `--json` requests JSON
diagnostics too. Exit classes are 0 success, 2 usage, 3 unreadable input,
4 semantic/numerical refusal, and 6 wall/design-evaluation budget exhaustion.

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network examples/cooling-network/size-mixed-slab.json
cargo run -p fs-cli --bin frankensim -- --json cooling-network examples/cooling-network/size-heterogeneous-hotspot.json
```

An optional `design` selects one surface coefficient, a temperature limit for
the existing objective, coefficient bounds, a full coupled-evaluation budget,
and temperature/log-coefficient tolerances:

```json
"design": {
  "surface": "last-face", "temperature_limit_k": 301.5,
  "min_htc_w_m2_k": 10, "max_htc_w_m2_k": 1000,
  "temperature_tolerance_k": 0.00001, "log_htc_tolerance": 0.0001,
  "max_evaluations": 80
}
```

`objective.gradient=true` is required. The original `mean_temperature_limit_k`
spelling is accepted only with a mean objective; it cannot silently become a
peak limit. Supplying both limit spellings refuses.

Each candidate changes h on both the FEM Robin operator and air exchanger,
solves the coupled field, and computes the total derivative. Hydraulics remain
fixed. Safeguarded Newton proposals outside the central 80% of the bracket use
bisection. A returned design contains the **passing evaluated field**, actual
coefficients, failed lower endpoint, log bracket width, trial history and work
counts. A feasible declared minimum returns `minimum-feasible`; otherwise
`target-bracketed` requires both tolerances. The single wall budget covers the
whole search. No partial field is published on evaluation exhaustion.

There is no global monotonicity theorem for arbitrary networks. If both endpoints
fail, the command reports a **missing passing endpoint**, not proof that no
interior design exists. A returned bracket identifies one local crossing, not a
globally minimum cooling coefficient. This sizes an effective coefficient, not
a fan, fin geometry or validated hardware.

## Reference cases and tests

`mixed-slab.json` is a 50 × 100 × 100 mm, 12-tetrahedron slab: a 3 L/s stream mixes
with a 1 L/s bypass before traversing the opposite surface. Inlets at 330/290 K
have an independent continuous reference of approximately 2.638172672 W through
the solid and an ideal final outlet of 320 K. The 323 K mean-wall sizing inverse
is approximately 194.501441754 W/(m² K).

`size-heterogeneous-hotspot.json` adds two material layers and one localized watt
at vertex 4 with both supplies at 300 K. Independent NumPy P1 FEM calculations
give a baseline first-surface mean of 301.118638 K but a whole-solid peak of
301.690663 K. The mean would pass 301.5 K while the actual discrete peak fails.
A separate bounded numerical search reaches a passing peak near 301.499991 K
at h ≈ 162.459319 W/(m² K). The active vertex changes across the search.
These are independent reference calculations, **not retained Rust execution**.

```bash
cargo test -p fs-cli --bin frankensim network_command::solid_data::tests
cargo test -p fs-cli --bin frankensim network_command::objective::tests
```

No new runtime dependency is introduced. Constant properties, fixed hydraulics,
conforming mesh and declared isotropic materials remain the model boundaries.
Temperature-dependent flow, buoyancy, recirculation, contact, radiation,
mesh-convergence bounds, uncertainty certification and experimental validation
are not inferred. The new Rust tests have not been executed in the authoring
environment, which has no Rust/Cargo, DSR or RCH.
