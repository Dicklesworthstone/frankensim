# Connected modal performances

```bash
cargo run -p fs-couple --bin music_render -- \
  modal crates/fs-couple/examples/coupled-modal.performance \
  /tmp/coupled-modal.wav --block 37
```

The example contains two structural components joined by a spring and dashpot.
Only the first component receives actuator forces; only the second contributes
observer pressure. Its sound therefore requires real mechanical energy transfer.
Replacing `connection 300000 10 0` with `connection 0 0 0` produces silence,
not an independently triggered receiver. Both components start in the settled
**full-network** equilibrium; the load is released before sample 37 and applied
again at sample 1200, then released at sample 1800. The 4801-sample duration
includes a short last callback. No input or output is normalized automatically.

## Input

The existing `frankensim-modal-performance-v1` format remains independent: its
meaning, input hash domain and command provenance are unchanged. Version 2 uses
`frankensim-modal-performance-v2` as its first record and adds the following
records **after all voices and before `events`**:

```text
coupling_limits MAX_CONNECTIONS MAX_SETUP_TERMS NYQUIST_FRACTION MAX_TOTAL_ENERGY MAX_PRESSURE MAX_REACTION_FORCE SOLVE_RTOL ENERGY_ATOL ENERGY_RTOL
connections COUNT
connection STIFFNESS_N_PER_M DAMPING_N_S_PER_M REST_EXTENSION_M
left COMPONENT_INDEX SHAPE_0 ... SHAPE_N_MINUS_1
right COMPONENT_INDEX SHAPE_0 ... SHAPE_N_MINUS_1
```

Repeat the last three records for each connection. Component indices refer to
`voice` records in file order. Each attachment has exactly as many shape values
as its component has modes. Shapes have units `1/sqrt(kg)` and must belong to
that component's existing mass-normalized basis. Signed shapes are important:
the relative displacement is `left^T q_left - right^T q_right - rest`, in metres.
A positive signed reaction acts along the left coordinate and its opposite
along the right coordinate. Two attachments may belong to the same component;
a zero shape vector explicitly declares an immobile attachment in this model.

Stiffness and damping must be finite and nonnegative; either or both may be
zero. Rest extension can be signed. These are fixed **bilateral** links, not a
unilateral obstacle or an infinitely stiff constraint. Force records keep the
v1 physical actuator semantics: newtons on a source component's port, applied
before the named sample. Releasing one actuator does not remove its connections.

`MAX_TOTAL_ENERGY` is joules including the springs, `MAX_PRESSURE` is summed
observer pascals, and `MAX_REACTION_FORCE` is newtons per connection. The setup
work screen is `total_modes * (connections + 1)^2`. File input caps connections
at 64 and setup terms at 16777216. The coupled Nyquist screen uses the algebraic
upper bound `max(omega^2) + sum(k * ||B||^2)`: it is conservative and is not an
interval eigenvalue certificate. Energy admission uses the supplied absolute
joule tolerance plus the supplied relative tolerance times the step's energy
scale. No tolerance is enlarged automatically.

Every component must use `retain-state`, or every component must use
`static-preload` with zero Q/V fields. Mixed initialization is refused. Coupled
preload solves the full spring stiffness system; independently settling each
component would omit the attachment reactions. Dashpots carry no static force.
The pre-window loading work is not reconstructed or claimed.

## Native composition and boundaries

`render::schedule::force::coupled::CoupledModalSystem` accepts existing
`ModalAcousticTimeModel` components and attachment maps. It keeps each
component's damping, pressure transfers, numerical ceilings and initial state.
`ScheduledRenderer::from_coupled_modal_forces` reuses the physical-force compiler,
sample scheduler and streaming output. The direct system exposes accepted
component states and the last complete reaction, storage, work and loss record.
Its `step_under_gate` can cancel before publication; the renderer observes
cancellation at complete host callback boundaries. Failed component/energy trials
publish no coupled state; a failed audio callback poisons its wrapper because an
earlier part of that callback may have advanced.

The component transition is the existing exact held-force modal step. Connection
forces use midpoint spring extension and average relative velocity. Their small
symmetric system is factored once using `fs-la`, and step buffers are reused.
This is a finite-step approximation of the continuously connected network,
**not** its exact matrix exponential. Refine the timestep to evaluate its error;
a small algebraic residual alone does not establish temporal accuracy.

Modal frequency, mass normalization, attachment geometry and acoustic transfers
must already be valid for the caller's reduced model. The reader does not infer
geometry, source material data, contact mechanics, rigid constraints, nonlinear
springs or changing topology. Acoustic transfers retain their original
narrow-band read-only interpretation; no new fluid backreaction or radiation
energy balance is implied. Fixed connection and scratch bounds are not a measured
real-time guarantee. The source and numerical checks for this contribution do
not establish Rust build success or experimental physical validation.

Focused native checks:

```bash
cargo test -p fs-couple --test modal_coupling --test coupled_render --test music_render_coupled
cargo test -p fs-couple --lib render::schedule::force::file
cargo test -p fs-couple --bin music_render
```
