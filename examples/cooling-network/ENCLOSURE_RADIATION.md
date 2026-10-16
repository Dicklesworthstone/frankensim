# Closed-enclosure radiation in coupled cooling

```bash
frankensim --json cooling-network \
  examples/cooling-network/enclosure-radiation-gap.json

frankensim --json cooling-network \
  examples/cooling-network/enclosure-radiation-pulse.json
```

Two separate solids exchange radiation across a gap while air cools both.
There is no thermal contact or welded node across that gap. Five watts enter
one body. Radiation can redistribute those watts between bodies, but cannot
remove energy from their combined closed enclosure. Only convection exports
heat to the air network. Materials, emissivities and fan data are illustrative,
not measured or experimentally validated.

The examples explicitly declare the ideal infinite-parallel-plate view-factor
limit. Their finite P1 faces represent equal-area samples of that idealization.
Their finite-geometry visibility, edge leakage and occlusion have NOT been
computed. Admitting a reciprocal matrix does not prove it describes the mesh.

For steady/trajectory derivatives, gradient-guided sizing and goal-recovery
marking, see [ENCLOSURE_ADJOINTS.md](ENCLOSURE_ADJOINTS.md) and the runnable
`adjoint-enclosure-pulse.json`. Both models now have total adjoints within the
existing fixed-grid, fixed-cycle derivative domain.

## Declare the radiation model

Choose `radiation.enclosure` instead of the existing `radiation.surfaces` list
of independent fixed-temperature surroundings. Supplying both refuses.

```json
"radiation": {
  "max_iterations": 200,
  "temperature_tolerance_k": 1e-9,
  "relaxation": 0.25,
  "enclosure": {
    "surfaces": [
      {"surface":"emitter","emissivity":0.8,"source":"declared finish"},
      {"surface":"receiver","emissivity":0.6,"source":"declared finish"}
    ],
    "view_factors": [[0,1],[1,0]],
    "row_sum_tolerance": 1e-12,
    "reciprocity_relative_tolerance": 1e-12,
    "evidence": {"kind":"analytic","geometry":"ideal infinite parallel plates"}
  }
}
```

Each named patch must already be a nonempty cooling surface. Areas come from
those mesh traces, not separately entered numbers. This integration requires
convection on the radiating patches; it is not a pure-vacuum solver. Between 2
and 64 patches are admitted. Matrix rows AND columns follow the input surface
order; the complete system is canonicalized by name before use. Coefficients
must lie in [0,1], rows must close to one, and area-weighted reciprocity must
pass the existing fs-conduction admission. Both numerical tolerances must be
explicitly supplied in [0,1e-8]. Supplied factors are not clipped or renormalized
to make an invalid enclosure pass. Emissivity lies in (0,1].

`evidence` may alternatively be an externally generated QMC matrix's retained
coordinates: `{"kind":"external-qmc","seed":"73","samples":"100000",
"generator":"named-generator-and-version"}`. Seed and positive sample count
are decimal u64 strings. These labels preserve the caller's provenance; supplying
this evidence does not run the named external sampler or independently verify
its formula.
No ambient-temperature field exists in a closed enclosure, and missing external
surfaces cannot be replaced by silently discarding unclosed view-factor rows.

## Generate finite-geometry view factors

Run the complete [heated-cavity example](enclosure-radiation-generated.json):

```bash
frankensim --json cooling-network \
  examples/cooling-network/enclosure-radiation-generated.json
```

This example uses a three-metre cube surrounding a one-metre cubic cavity,
represented by 156 tetrahedra. Its 26 m3 of solid receives a uniform
0.2 W/m3 source, for 5.2 W total input. The exterior is adiabatic. Two named
groups of cavity faces exchange radiation and transfer heat to the existing
air network with prescribed coefficients of 1 and 2 W/(m2 K). The geometry,
conductivity, emissivities, heat-transfer coefficients and fan data are
illustrative; they are not measured hardware or a physically validated model.
At steady state, the solid and air balances must account for the 5.2 W input,
while the signed internal radiation transfers sum to zero. No particular
temperature value or sampling-accuracy bound is asserted.

An enclosure can instead sample visibility from its actual tetrahedral mesh.
Keep the named surface declarations and replace all four supplied-matrix fields
(`view_factors`, `row_sum_tolerance`,
`reciprocity_relative_tolerance` and `evidence`) with:

```json
"generate_view_factors": {
  "seed": "41",
  "samples_per_replicate": 512,
  "replicates": 4,
  "max_triangles": 256,
  "max_ray_triangle_tests": "1000000",
  "projection": {
    "max_iterations": 10000,
    "row_sum_tolerance": 1e-12,
    "max_abs_adjustment": 0.25
  }
}
```

The two modes are mutually exclusive. Unknown sampling and projection fields
refuse. The seed and ray/triangle budget are decimal u64 strings;
`samples_per_replicate` is a positive power of two fitting u32;
`replicates` lies in 2..=32. The triangle limit lies in 1..=1,000,000.
These example budgets are illustrative, not a sampling-accuracy guarantee.

Each replicate uses an independent Owen-scrambled Sobol stream. Samples select
positions uniformly by patch area and directions with the diffuse cosine
distribution into the outward hemisphere. The first intersected boundary
triangle determines visibility. Every boundary triangle is opaque, including
triangles outside the named radiating patches. A first hit on an unassigned
triangle or on a back face is recorded separately from a ray that escapes the
mesh; none is reassigned to an unobstructed named patch.

The existing gap examples declare an ideal infinite-plate matrix. Substituting
the generator into those finite, open geometries is expected to refuse when
rays escape. To use generated factors for a closed enclosure, the named patches
must cover the cavity seen by their outward hemispheres. Zero sampled escaped,
unassigned and back-face rays are necessary for this implementation's closed
projection; they do not prove geometric closure or watertightness.

Sampling retains its raw matrix and independent-replicate empirical standard
errors. A separate, explicit projection adjusts area-weighted exchange to
reciprocity and row closure. Its iteration cap is 1..=100,000, its positive
row tolerance is at most 1e-8, and `max_abs_adjustment` is an absolute
view-factor-entry correction budget in [0,1). The projected matrix is admitted
only if the declared correction and closure limits are met. Missing sampled
exchange is not silently repaired. These empirical standard errors are not
calibrated confidence intervals, and the projection is a numerical conservation
correction rather than an error certificate.

Generation takes place during the first physical bind under the command's
existing cancellation context and wall deadline. Before sampling, the product
of patches, replicates, samples per replicate and opaque boundary triangles
must fit the declared ray/triangle budget. Exceeding this bound or the triangle
limit returns the existing radiation-budget exit code 6. A successful record is
retained throughout later steady iterations or transient endpoints for that
parsed request. Failed or cancelled generation is not retained. The retained
record is checked against the complete mesh, selected traces and emissivity
declarations; a newly parsed refined-mesh request generates its own factors.

The result's `radiation.view_factors` is the projected matrix used by the
radiosity solve. Its `radiation.view_factor_evidence` retains:

- the seed, sampling configuration, declared work limits and actual ray work;
- exact mesh and estimate identities, canonical surface names and areas;
- `raw_view_factors`, `empirical_standard_errors` and `hit_counts`;
- per-surface `escaped_counts`, `unassigned_counts` and `backface_counts`;
- declared projection limits, actual iterations and maximum entry adjustment.

Counts and u64 work totals use decimal strings, as declared by
`count_encoding: "decimal-u64-string"`. Matrix axes and count rows follow
the canonical surface-name order. Per emitting surface, hits plus all three
loss counts equal `samples_per_surface`.

This path computes nominal finite-mesh visibility and occlusion. It does not
certify the geometry, calibrate stochastic error, model participating media or
validate physical material data. Existing temperature and inlet derivatives
keep the generated matrix fixed; derivatives of geometry or visibility are not
provided.

## Numerical and energy coupling

The existing `GrayDiffuseEnclosure` owns the radiosity solve

```
J_i - (1-epsilon_i) sum_j F_ij J_j = epsilon_i sigma T_mean_i^4
q_i = J_i - sum_j F_ij J_j.
```

Each q_i is a UNIFORM outward flux over its original patch. Both reflected
and emitted radiation participate. It is not pointwise integration of T(x)^4.
The adapter adds this flux alongside convection by applying the algebraically
equivalent Robin reference `T_air - q_i/h_air`, leaving h_air unchanged.
That shifted number is an assembly device, not a physical air temperature;
it may be negative and is never passed to the air model or reported as an inlet.
Published solid temperatures must remain finite and strictly positive.

Every iteration solves the real conduction problem, including assigned
conductivity laws and matching/nonmatching contact operators. It then recomputes
radiosity at the solved patch means. The original radiation iteration cap and
relaxation apply; the raw mean-temperature residual and EACH patch's applied
versus recomputed watt mismatch must pass. The radiosity equation residual and
closed-enclosure heat sum are separately checked. The independently accumulated
Robin boundary integral must split into convection and the applied patch flux,
and the full source/air/solid energy gates still apply.

`radiation.surfaces` in the RESULT lists emissivity, mean temperature, radiosity,
irradiation, signed outward heat and applied heat. Negative outward heat is
absorption by that solid, not an error or value to clip. `radiative_out_w` is
the near-zero closed-enclosure balance residual, NOT an external energy sink.
In a steady result, `iterations` refers to the last inner radiation solve.
`forward_solid_solves` counts inner solid evaluations across air coupling;
`total_solid_solves` also includes any requested adjoint reconstruction.
All iterations share the original command wall deadline. Exhaustion returns
no partial success.

## Transient heating, cooldown and repeated cycles

The pulse example uses two 30-second cycles, each with ten seconds at five watts
and twenty seconds at zero watts. It has nonlinear conductivity and different
constant heat capacities in the two bodies. The deliberately small declared
capacities (20,000 and 10,000 J/(m3 K)) expose storage effects over a short example;
they are not measured material properties. The second cycle starts from the
first cycle's accepted final field, not from a fresh 300 K state.

No additional radiation option is needed: the ordinary `transient` declaration
selects the existing backward-Euler driver. Every inner radiation and air trial
uses the SAME immutable previous physical temperature field. The old surface
means initialize radiation iteration, but emission is recomputed at the NEW
solved means. Freezing radiation at the old temperature is not this model.
The timestep acceptance equation is

```
stored-energy change = dt * (input power - air heat gain - enclosure residual).
```

Signed patch transfers cancel internally. The existing
`radiative_energy_loss_j` field therefore accumulates only the near-zero closure
residual for this CLOSED model, not the emitter's positive transfer or the sum
of absolute patch transfers. It is not a measure of heat exported by radiation.
The final `radiation` report retains both individual signed patch powers and
`temporal_scope: "final-accepted-endpoint"`. Its `iterations` and
`total_solid_solves` are null: final-field heat recomputation cannot recover
those counts. The trajectory's work counters count every actual solid callback,
including inner radiation iterations, discarded adaptive trials and any requested
adjoint reconstruction. Physical energy belongs only to accepted forward steps.

Adaptive full steps and rejected half-step pairs never enter accepted storage
or energy history. Fixed-count repetition, ordinary periodic/controller runs,
workload/fan sizing, and full-trajectory timestep studies all reach the same
endpoint producer; none substitutes a convection-only trajectory. Their original
timestep, cycle, trial, refinement, solver and wall limits remain in effect.
Complete trajectories retain their usual independent energy checks.

Steady enclosure adjoints and goal-recovery spatial marking are supported, as
are fixed-grid, fixed-count trajectory adjoints and their sizing consumers.
Adaptive-time-grid, automatic timestep-selection, periodic-stop and controller
derivatives still refuse explicitly. There is no frozen-radiosity fallback.
Finite-mesh view-factor generation and opaque first-hit occlusion use the optional
sampling path above. Participating media, changing geometry during a solve,
continuous-time peak certification and native .fsim lowering are not implemented.

## Uncertain surface finish

```bash
frankensim --json cooling-network-uq \
  examples/cooling-network/enclosure-radiation-gap.json \
  examples/cooling-network/uq-enclosure-finish.json \
  --checkpoint enclosure-finish.uqcp

frankensim --json cooling-network-uq \
  examples/cooling-network/enclosure-radiation-pulse.json \
  examples/cooling-network/uq-enclosure-pulse.json \
  --checkpoint enclosure-pulse.uqcp
```

Each checkpoint destination must be new. The existing `radiation-emissivity`
target resolves against the named enclosure patch and changes the actual input
before each cooling solve. Geometry, view factors, power and other emissivities
stay fixed. For transient UQ, one draw applies throughout the entire trajectory
and its repeated cycles. The observable is the all-cycle sampled peak, not the
cooled final temperature. An interrupted trajectory contributes no partial peak
and retries the same ordinal on resume. A surroundings-temperature target
refuses for an enclosure; nonphysical samples are terminal errors, not clipped
or redrawn. Existing checkpoint and confidence policies remain in use. The
example probability distributions are illustrative, not empirical evidence.

## Focused checks and boundaries

```bash
cargo test -p fs-cli --test cooling_enclosure_radiation
cargo test -p fs-cli --test cooling_enclosure_transient
cargo test -p fs-cli --test enclosure_radiosity_adjoint
cargo test -p fs-cli --bin frankensim uq_command::model::radiation
```

The geometry submodule in the enclosure forward test target adds G0 count,
conservation and refusal checks plus G3 surface-order replay. It exercises a
closed cavity formed by 156 tetrahedra through the actual command and checks
the retained raw estimates and correction; the open gap must refuse. The complete
heated example checks positive temperature rise, the 5.2 W solid and air heat
balances, and zero net internal radiative export. Run just these additions with `cargo test -p fs-cli --test cooling_enclosure_radiation geometry::`.

The original sixteen enclosure forward/UQ regressions remain. Former missing-
adjoint checks now exercise the still-unsupported adaptive derivative policy.
Thirteen additional derivative/consumer regressions are described in
ENCLOSURE_ADJOINTS.md. These Rust tests were executed natively on
2026-09-25 and pass.

Independent Python P1 calculations use direct endpoint equations with analytic
two-plate radiation and air elimination, and a separate nested radiosity/solid
iteration. Twelve transient endpoint comparisons had at most 3.65e-11 K field
difference. Two manufactured endpoints had at most 1.71e-13 K field error;
freezing radiation at the old field instead caused 0.09831 and 0.60137 K errors.

For the two-cycle nonlinear pulse example, the independent sampled peak is
315.427903 K at 40 seconds. Its 100 J input splits into 35.006204 J stored and
64.993796 J exported by air. About 15.594296 J is redistributed internally from
emitter to receiver, with zero net radiative export in the analytic reference.
Freezing old-temperature radiation changes the peak by 0.10510 K. Refining from
30 to 60 and 120 total steps gives reference peaks 315.427903, 315.564990 and
315.634980 K, illustrating why a successful algebraic solve is not a timestep
error bound. These are independent numerical references, NOT executions of
FrankenSim, runtime speedups, continuum certificates or physical validation.
FrankenSim itself (measured 2026-09-25, release build, after 05db922bf made the transient capacitance row-sum lumped; the reference above uses the consistent P1 mass, so the two differ on this 12-tet mesh) samples 315.2480936 K at 40 seconds, with 34.9980786 J stored and
65.0019214 J exported by air.
