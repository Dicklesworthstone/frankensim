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
mode does not add arbitrary boundary conditions, 3-D physics or external geometry
imports. The explicit load-family stress-restoration option is described in
`MULTI_LOAD_STUDY.md`; without that option, initialization remains strictly
stress-feasible as described below.
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

## Protected material and empty regions

`bracket-protected-regions-2d.fsim` adds material that optimization may not
remove and empty regions that it may not fill. Run it with the same native
commands, for example:

```sh
cargo run --release -p fs-cli -- --json study \
  examples/marquee/bracket-protected-regions-2d.fsim ./regions.db --budget 1
```

The optional final optimizer field, after `:stress-tolerance-pa`, is:

```lisp
    :design-regions (
      (region :phase material :lower (0.125 0.125) :upper (0.25 0.25) :phi-margin 0.01)
      (region :phase void :lower (0.3 0.4) :upper (0.32 0.42) :phi-margin 0.01)
    ))
```

The final `))` closes this list and the optimizer. Declare 1..=64 rectangles
with ordered bounds in the unit-square coordinate system. `phi-margin` must
be finite and positive. It is a **field-value margin**, not a certified
distance, wall thickness, or machinability requirement.

The existing region implementation prescribes every corner of every cell
intersecting a rectangle's interior. Thus the entire bilinear cell, not just
its centre, retains the requested sign; coverage may extend by less than one
cell per side. Opposite phases sharing a required node refuse even when their
rectangles do not overlap geometrically. Same-phase overlaps impose the stronger
margin. Neither policy may change the existing support/load boundary traces.

Regions are imposed before area projection and the independent baseline solves.
The resulting baseline must still satisfy the declared area and stress limits;
regions do not relax either gate. Fixed node bits survive every accepted update.
The source and constraint reports retain the declarations, and resume rebuilds
their original prescriptions from the saved source before admitting the saved
geometry. It never derives new prescriptions from the optimized field or repairs
a violating endpoint silently. Resume needs no original input file. A cancelled
region preparation/reconstruction returns no partial field and leaves the prior
receipt available. All existing source files without this field remain valid.

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
