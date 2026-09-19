# Shared-body modal friction

```bash
cargo run -p fs-couple --bin music_render -- \
  modal crates/fs-couple/examples/modal-friction.performance \
  /tmp/modal-friction.wav --block 37
```

The example has three two-mode bodies and two compliant contacts. Only the
center body's initially stationary tangential mode is observed. Normal contact
alone cannot excite that mode; setting both friction coefficients to zero gives
exact receiver silence. The outside bodies carry initial tangential motion and
later authored force events. Negative initial gaps explicitly retain compressed
contact potential, not a silently solved static preload. The 601-sample output
includes a short final callback. These are synthetic reduced-model inputs, not
an experimentally calibrated instrument or interface.

## Why this adjacent capability

The comprehensive plan's power-port composition and the music plan's reuse
ledger already assign the friction law to `fs-tribo` / `StribeckFriction`.
The missing connection was a jointly solved load and tangential motion on the
same shared modal bodies, not another constitutive implementation. This extends
that composition seam. It does not close the broader finite-patch/partial-slip
Beads (`frankensim-b8bxd.7`) or replace the distinct Moreau-Jean/global-impact
lane (`frankensim-ext-contact-nonsmooth-lane-oh0i`). Their original acceptance
requirements and scientific authority boundaries remain unchanged.

## Version 5 input

Use `frankensim-modal-performance-v5`. All version-4 component, bilateral,
normal-contact and work-budget records retain their meanings. After the entire
normal-contact set and before `events`, supply one explicit friction record for
EVERY normal contact, in the same order:

```text
frictions NORMAL_CONTACT_COUNT
friction regularized-coulomb MU RAMP_SPEED_M_S MAX_FORCE_N SOURCE_LABEL
friction_left LEFT_COMPONENT SHAPE_0 ... SHAPE_N_MINUS_1
friction_right RIGHT_COMPONENT SHAPE_0 ... SHAPE_N_MINUS_1
friction none
```

The last line illustrates a second, frictionless contact; do not supply shape
rows for `none`. For a friction law, component indices must match that normal
contact's left and right bodies. Tangential shapes have units `1/sqrt(kg)` and
must use those bodies' existing mass-normalized mode bases. Different normal
and tangential projections may share modes; their Euclidean coefficient dot
product need not be zero. Geometry/frame validity remains the author's duty.
All decoded actuator, connection, normal and tangential shape coefficients share
the existing 65536-coefficient budget. Source labels are bounded to 1024 bytes.

`MU` is finite and nonnegative. `RAMP_SPEED_M_S` and `MAX_FORCE_N` are finite and
strictly positive. The ramp speed is an authored physical-model regularization,
not a solver tolerance. It is used by the existing `StribeckFriction` adapter
with equal static/dynamic coefficients. No static reaction, velocity-weakening
curve, material property, or missing law is inferred. An all-`none` set is legal.
The sidecar records the v5 input hash, total normal contacts, authored non-`none`
friction count (including zero-coefficient laws), and the selected friction
rung. Earlier schemas retain their original compiler, hash and sidecar paths.

## Joint physics and accepted work

Let `B` be the signed normal shapes, `C` the signed tangential shapes, and `D`
the existing bilateral network's condensed held-force displacement response.
For compressive normal reactions `R` and signed tangential reactions `T`:

```text
x1 = x_free - (B^T D B) R - (B^T D C) T
Delta_y = y_free - y0 - (C^T D B) R - (C^T D C) T
T_i = mu_i R_i clamp(Delta_y_i / (dt * ramp_speed_i), -1, 1)
applied modal force = external - B R - C T
friction dissipation = sum_i T_i Delta_y_i
```

All four blocks include shared-body and bilateral responses. Normal loads are
not frozen before friction is evaluated. A deterministic coordinate solve must
pass every normal and tangent residual together. It then stages exactly one
physical sample, rechecks both laws on the actual endpoint motion, and verifies
nonnegative applied tangential work plus whole-system energy closure. Authored
external work excludes solved internal reactions. Mechanical dissipation is
reported; no thermal partition is claimed.

Each tangent inherits its normal contact's root-iteration and force-residual
budgets. `MAX_SWEEPS` bounds their joint solve. The existing `MAX_SETUP_TERMS`
now covers original normal setup plus the friction extension using the
conservative bound `3*p*(n*(k+2)+(k+1)^2)+4*n*p^2`, where `p` is the normal-contact
count, `n` the total mode count and `k` the bilateral-connection count. Force
ceilings and residual tolerances are never enlarged to make a step succeed.

Native entry points are `MultiContactModalSystem::with_friction` and
`ScheduledRenderer::from_frictional_modal_forces`. Accepted states, reports and
sample clocks change together; a refused/cancelled direct step is retryable.
The existing callback host remains poisoned after a failed partial callback,
and the scheduler retains its existing cancellation and atomic-control rules.

## Limits and verification status

This is one-dimensional, continuous regularized Coulomb friction on fixed
compliant contacts. It is **not** set-valued sticking, a 2-D cone, finite-patch
microslip, velocity weakening, contact discovery, free rigid-body integration,
or a rigid-impact law. Scalar monotonicity does not guarantee convergence of
the complete coupled problem; budget exhaustion refuses. Time refinement is
still required. No hard-real-time, alias-free, experimental-validation or
repository-wide quality claim is made.

The following native regression targets are included but were **not executed**
in the authoring environment, which has no Cargo/Rust/RCH installation:

```bash
cargo test -p fs-couple --test modal_contact_friction --test music_render_friction
cargo test -p fs-couple --test modal_multi_contact --test music_render_multi_contact
cargo test -p fs-couple --test modal_contact --test music_render_contact
cargo test -p fs-couple --bin music_render
```

The tests exercise independent analytic/work checks, shared normal-load
feedback, direction reversal, permutation, exact zero-law limits, redundant
maps, budgets, cancellation/retry, file/native equivalence, hostile inputs,
callback partitioning, the real command/encoder, source identity, relocation
and non-overwrite behavior. Separate development-only Python numerical checks
were executed; they do not establish that these Rust targets compile or pass.
