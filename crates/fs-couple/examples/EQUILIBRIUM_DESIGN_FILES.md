# File-driven physical inverse design

```bash
cargo run -p fs-ascent --features equilibrium-design --bin equilibrium_fit -- \
  crates/fs-couple/examples/equilibrium-design.model \
  crates/fs-couple/examples/equilibrium-design.fit \
  --iterations 128 --evaluations 256 --tolerance 1e-8
```

The command is not tied to the old two-body `contact_fit` rig. The model file
supplies components, modes or physical free masses, port maps, springs, normal
contacts and numerical limits. The design file supplies independent experiments,
physical observations, unknown fields, sharing, scales and bounds. Both files
are loaded once; every candidate is evaluated through `EquilibriumDesign`, and
`EquilibriumStudy` supplies the existing bounded SQP engine. No new optimizer,
mechanical solver, finite-difference production gradient or dependency is added.

The shipped pair describes a supported 40-gram mass and an elastic receiver.
Its three load cases and six displacements are **synthetic**, independently
calculated from quadratic equilibrium with planted support stiffness 600 N/m,
contact coefficient 1.8e8 N/m² and gap 0.2 mm. They are not specimen data or
experimental calibration evidence. The numerical reference recovers these
parameters; the Rust executable and tests remain unexecuted in the authoring
environment because Cargo/rustc are unavailable.

## Model input

Use the existing modal-performance v2, v3 or v4 syntax. Every voice must be
`retain-state` or `free-mass`, with zero supplied position/velocity and zero
initial port loads. End with `events 0`. The design file owns all stationary
loads. Preload voice declarations, existing vibration, nonzero initial loads,
scheduled events and v5 friction models refuse rather than being discarded.
Even a contact-free model uses v2 with explicit `coupling_limits` and
`connections 0`; no numerical limits are invented for an old v1 input.

Each load or observation names a model `port` by zero-based component and port
index. Its original mass-normalized shape is copied exactly. Add a zero-load
port to declare a distinct sensor footprint; an observation is not required to
use an actuated port. A unit translation of a 0.04 kg mass has shape 5, not 1.
Clock, horizon, acoustic transfer and PCM scale records retain the original
model parser's admission, but are not a stationary displacement objective or an
audio-trajectory fit. Use `samples 1` for a static model template.

## Design records

Records are ordered, whitespace-separated and have no ignored fields:

```text
frankensim-equilibrium-design-v1
preload_limits MAX_CONTACTS MAX_SWEEPS MAX_SETUP_TERMS
sensitivity_limits MAX_CONTACTS MAX_SETUP_TERMS MAX_QUERY_TERMS MIN_SWITCH_DISTANCE_M
design_limits MAX_CASES MAX_VARIABLES MAX_BINDINGS MAX_PORTS_PER_CASE
cases COUNT
case NAME LOAD_COUNT TARGET_COUNT
load COMPONENT PORT FORCE_N
... repeat loads ...
target COMPONENT PORT TARGET_M SCALE_M WEIGHT
... repeat targets, then remaining cases ...
variables COUNT
variable NAME REFERENCE SCALE MINIMUM MAXIMUM BINDING_COUNT
bind FIELD INDEX
... repeat bindings, then remaining variables ...
```

`FIELD` is `spring-stiffness`, `spring-rest`, `contact-stiffness`, `contact-gap`,
`contact-weight`, or `actuator-force`. The last form has **two indices**:
`bind actuator-force CASE LOAD_INDEX`; it addresses an authored load record,
not a model port. All other indices refer to model connection/contact order.
Several compatible bindings can share one variable. Duplicate assignments and
mixed physical parameter kinds refuse through the existing design admission.

Physical values are `reference + scale * decision`, starting at decision zero.
Scales are positive and the finite physical bounds contain the reference.
Spring stiffness is N/m, spring rest and contact gaps are metres, and actuator
forces are newtons. Contact stiffness retains its source exponent and weight
convention; no exponent or damping fit is inferred. Each target contributes
`0.5 * weight * ((prediction_m-target_m)/scale_m)^2`. Its positive scale and
nonnegative weight are explicit, not inferred noise or confidence estimates.

The original component, connection, contact and derivative limits remain active.
Design preload limits may tighten, never enlarge, a model's v4 joint-contact
caps. Resource hard ceilings are 4 MiB model input, 1 MiB design input, 64 cases,
128 variables, 1024 total bindings and 1024 load/target ports per case. Expanded
port maps across the **whole** design are independently capped at 65536 shape
coefficients, preventing short reference rows from allocating unbounded arrays.
Each case needs at least one target. Names are unique within their family and
at most 128 bytes. Truncated, extra, nonfinite and unknown records refuse.

## Execution and results

`--iterations` permits 0..512 additional accepted steps. `--evaluations` permits
2..4096 total physics objective attempts, including initialization and a reserved
final complete re-solve. The case-solve ceiling is that total times the number
of independent cases. `--max-kkt-dimension` defaults to 384 and permits 1..384;
the existing small-dense SQP cap counts n decisions plus both n bound faces.
`--tolerance` defaults to 1e-8 and applies to the existing KKT residuals in scaled
decision coordinates. A budget or stall result is not relabelled convergence.

A successful invocation prints one JSON object to stdout: both exact input
hashes, actual stop reason and convergence flag, objective, physical parameter
values, decision-space bound multipliers, all case predictions, original
adjoint residuals and spent work. Input paths do not change result identity.
The accepted physical evaluation is recomputed within the budget and compared
with the optimizer's retained value/derivatives before reporting. Invalid input
or a physical refusal prints an error to stderr and no partial JSON result.
The command does not overwrite source files or create output artifacts itself.

Native API: `fs_couple::render::schedule::force::file::design::EquilibriumDesignFile`.
Its immutable `problem()` can be passed to a caller-owned `EquilibriumStudy` for
bounded segments, cancellation and continuation using the existing work ledger.
The command does not install signal handlers or serialize an interrupted study.

This remains local, fixed-basis, static normal-contact inverse design: no global
optimality, experimental identification, identifiability, mode extraction,
frictional equilibrium, contact-only free-body grounding or transient audio
sensitivity is implied. Contact activity margins can still refuse a trial.
Numerically out-of-bound trials retain the existing SQP adapter's refusal;
this command does not relax domain checks to force a boundary result.

```bash
cargo test -p fs-couple --test equilibrium_design_file
cargo test -p fs-ascent --features equilibrium-design --bin equilibrium_fit --test equilibrium_fit
```
