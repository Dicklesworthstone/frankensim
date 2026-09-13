# File-driven coupled cooling network

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network examples/cooling-network/mixed-slab.json
cargo test -p fs-cli --bin frankensim network_command::tests
cargo test -p fs-cli --test cooling_network
```

This **experimental binary command** accepts a JSON request containing an actual
linear tetrahedral solid, named exterior cooling faces, a quadratic hydraulic
graph, independent pressure-reservoir temperatures, and explicit numerical
budgets. It runs the existing `LossGraph`, `TransportNetwork`, shared-solid
conjugate solver and optional coupled adjoint. The request is not a `.fsim`
project, and this command does not perform CAD import, create a ledger, or claim
a package/report workflow. Existing `fs_cli::run`/`run_os` library verbs are
unchanged; `cooling-network` is dispatched by the `frankensim` binary.

Unlike a fixed example, every vertex, tetrahedron, surface, coefficient, heat
source, hydraulic connection and source temperature comes from the request.
Supply a conforming, non-overlapping tetrahedral mesh. Admission checks indices,
repeated cells, face incidence and exterior-face ownership; those checks are not
a geometric non-overlap certificate. Conductivity is one declared isotropic,
temperature-independent scalar; the source is a uniform volumetric density.

## Request contract

`schema` is exactly `frankensim.cooling-network.v1`, `units` is exactly `SI`, and
`seed` is a decimal u64 **string**, avoiding JSON floating-point seed rounding.
All quantity keys carry coherent-SI units. Unknown fields and duplicate keys
refuse. Every declared field in the example is required unless stated otherwise.

`budgets` explicitly bounds graph sweeps, shared-solid coupling iterations,
Krylov iterations per solve, adjoint interface iterations, and numerical wall
time. Wall timing starts **after** input/mesh admission and applies to the
whole numerical invocation; expiry requests the existing cancellation gate and
prevents publication after draining. A checkpoint is cooperative, not a hard
real-time deadline. The input cap is 16 MiB, with fixed ceilings of 20,000 solid
vertices, 100,000 tetrahedra, 4,096 hydraulic nodes and 16,384 branches.

`tolerances` names absolute flow, heat and coupling-temperature gates, a relative
linear/derivative tolerance, and fixed coupling/interface relaxation in `(0,1]`.
The command does not silently substitute relative allowances for declared
absolute heat or hydraulic tolerances. Excessively tight tolerances may refuse.

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
areas are integrated from the mesh, never copied from a separate input area.
Unlisted exterior faces are insulated only when `adiabatic_remainder` is true;
otherwise every exterior face must be owned. A zero-flow heat exchanger refuses.

`objective.mean_wall_region` selects an **area-mean**, not a nodal maximum.
`objective.gradient` requests total inlet-temperature and log(h) sensitivities,
including the solid matrix/load change and downstream heated-air feedback.
Other geometrical and hydraulic quantities are held fixed.

## Results and example

Successful output is one JSON document, including actual nodal solid temperatures,
node mixture temperatures, signed branch flows, branch inlet/outlet temperatures,
surface heat rates and references, mean-wall objective, optional total gradients,
and independently accumulated solid/air heat-balance diagnostics. `null` means
uncomputed or stagnant/unknown, never zero uncertainty. No result is published
on a failed producer, exhausted iteration budget, or expired wall budget.
`--json` requests JSON diagnostics too. Exit classes reuse the CLI's definitions:
0 success, 2 usage, 3 unreadable input, 4 semantic/numerical refusal, 6 wall budget.

The supplied 12-tetrahedron example is a 50 × 100 × 100 mm slab. A 3 L/s heated
stream mixes with a 1 L/s bypass before traversing the opposite slab surface.
The two inlets are 330 K and 290 K. The continuous constant-property slab oracle
predicts about 2.638172672 W through the solid; the ideal final outlet is 320 K
because the slab has zero net source. These are analytic reference values, not
retained Rust-run measurements. Unit tests exercise the actual input adapter and
FEM APIs; the binary integration test exercises real command dispatch.

All coefficients and mesh data are illustrative declarations. Results remain
nominal estimates: no temperature-dependent flow, buoyancy, recirculation, contact,
radiation, discretization bound, uncertainty certification, material-data authority
or experimental validation is inferred. No new runtime dependency is introduced.
