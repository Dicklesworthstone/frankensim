# Constrained elasticity in the native study command

`bracket-projected-stress-2d.fsim` opts the existing `frankensim study` command
into hard numerical material-area equality and a sampled von Mises limit. It
uses the existing level-set/CutFEM optimizer; it is not a density-only study or
a shell wrapper around the separate marquee executable.

```sh
cargo run --release -p fs-cli -- --json study \
  examples/marquee/bracket-projected-stress-2d.fsim ./constrained.db --budget 1
# For this two-update example, a successful first update returns exit 6 and a
# budget-exhausted receipt. Copy its run_id (study- followed by 64 hex digits).

cargo run --release -p fs-cli -- --json study --resume \
  "$RUN_ID" ./constrained.db
# Copy the new run_id before exporting the final or partially completed result.
cargo run --release -p fs-cli -- --json report "$RUN_ID" ./constrained.db
cargo run --release -p fs-cli -- --json package "$RUN_ID" ./constrained.db
```

The example is a **software/numerical exercise**, not an engineering design.
Its 2 Pa Young's modulus, 1 Pa traction, and deliberately loose 1e12 Pa stress
limit are not a calibrated material or a safety recommendation. The model is
still the native 1 m × 1 m, left-clamped, downward right-loaded plane-strain
plate with interior circular holes and the existing unit-thickness convention.
Change the complete explicit source to declare another admitted study. This
mode does not add arbitrary boundary conditions, 3-D physics, external geometry
imports, material/void masks, or stress-feasibility restoration to native `.fsim`.
The separate marquee commands retain their own broader geometry input features.

## Explicit policy

After `:steps` in the optimizer section, `:constraint-mode projected-stress`
requires area tolerance, maximum field shift, maximum area evaluations,
candidates per update, contraction, minimum relative compliance decrease,
CG polling interval, sampled-stress limit and absolute stress tolerance.
The tracked example spells every field in canonical order. Missing, repeated,
unknown, or reordered fields refuse rather than silently selecting defaults.
Without this mode the original native elasticity declarations and behavior stay
unchanged; their augmented-Lagrange volume control is not a hard area guarantee.

The area target is the existing `:volume-fraction` times the unit plate's area.
Projection fixes both boundary nodal traces before independently solving the
baseline. Stress must also be feasible before any accepted baseline is published.
A candidate must decrease compliance relative to the current feasible design
and satisfy both constraints. The reported reduction is against the ORIGINAL
feasible study baseline, never the overfilled input or a freshly reset resume
baseline. All stress claims cover only the deterministic sampled set.

## Retention, continuation and stopping

The existing ledger retains source, exact nodal bits, iteration rows and HTML/
JSON reports. Receipts additionally retain the feasible baseline, each accepted
stress evaluation, candidate counts and final search-refusal reasons. Every
accepted update is committed before more physics. Native resume needs only its
run_id and ledger; it preserves the original policy, iteration target and
lifetime wall charge. It independently re-solves the retained endpoint and
requires exact mechanics/stress agreement under the same executable, without
replaying earlier geometry updates. Earlier receipts are not overwritten.

Exit 0 / `completed` means only the declared update count was reached, NOT KKT
convergence. Exit 6 / `budget-exhausted` retains partial accepted work. A bounded
search stall is `no-feasible-descent` with exit 4 and an exportable retained
result, not an optimum. Cancellation retains the prior admitted state; a stop
before the first feasible baseline publishes no feasible receipt. Numerical or
ledger errors fail explicitly and identify a previously saved state when one
exists. The existing package remains structural evidence, not a physical or
continuous-stress certificate; reports display the actual constraint evidence.

CG and stress-cell checks observe the cooperative wall/cancellation policy.
Assembly, one area evaluation, and ledger I/O are not preemptible. The native
memory declaration is an admission gate, not a measured peak-memory guarantee.

Focused native tests: `cargo test --release -p fs-cli --lib study::elasticity`
and `cargo test --release -p fs-cli --test study_checkpoint_cli`. These new
regressions require native execution; source inspection alone does not establish
that this example accepts a step or that a checkpoint replays successfully.
