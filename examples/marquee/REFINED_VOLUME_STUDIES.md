# Continue a native volume study on a finer grid

A `projected-volume` study can start from an actual retained endpoint instead
of restarting its original holes. Add `:refine-from "study-..."` immediately
after `:cg-poll-iters` and before any `:mesh-check` or `:design-regions` block.
The identifier must name a sealed native projected-volume receipt in the same
ledger. Set `:mesh-level` to exactly one greater than that source study.

This is a NEW, explicitly funded study, not a budget-resetting `--resume`.
The source is unchanged, including its work charges and terminal status.
Normal `study --resume` on the new receipt continues the fine study exactly.

## Run the supplied example at two levels

First retain an actual coarse update. This command intentionally returns the
budget exit code after one update; retain its JSON output, including `run_id`.

```sh
cargo run -p fs-cli --bin frankensim -- --json study \
  examples/marquee/bracket-projected-volume-2d.fsim studies.db \
  --budget 1 > coarse-result.json
```

The following development-only snippet creates a new source for this specific
supplied example, using the REAL returned receipt. It refuses to overwrite an
existing `fine.fsim`. Python is not a production dependency.

```sh
python3 - <<'PY'
import json
from pathlib import Path
result = json.loads(Path('coarse-result.json').read_text())
assert result['receipt']['iterations_completed'] == 1
pointer = result['run_id']
assert pointer.startswith('study-') and len(pointer) == 70
text = Path('examples/marquee/bracket-projected-volume-2d.fsim').read_text()
assert text.count(':mesh-level 3') == 1
assert text.count('    :cg-poll-iters 1)') == 1
text = text.replace(':mesh-level 3', ':mesh-level 4', 1)
text = text.replace('    :cg-poll-iters 1)',
    '    :cg-poll-iters 1\n    :refine-from "' + pointer + '"\n  )', 1)
with Path('fine.fsim').open('x') as output:
    output.write(text)
PY
cargo run -p fs-cli --bin frankensim -- --json study fine.fsim studies.db
```

Use the new command's actual `run_id` with the usual `study --resume`, `report`,
and `package` commands. A new grid does not guarantee descent or convergence;
the existing no-descent, cancellation, work-budget and optional mesh-check
outcomes still apply.

## What is preserved and what changes

The complete canonical numerical declaration is compared with the source.
Material, physical load, seed, original geometry declarations, protected
regions, material-area target/tolerance and proposal/search policy must match.
Metadata, new wall/memory/update budgets, the one-level refinement, and an
optional mesh-check policy may change. No stress policy is silently dropped:
both ends must explicitly select `projected-volume`.

The coarse field comes from the receipt's exact retained nodal bits, not the
source's original holes. Its current mechanics are independently re-admitted.
The existing bilinear prolongator then transfers it. Prescribed nodes inherit
the nonzero-weight coarse support; protected regions are not re-rasterized in
a way that changes their already-admitted footprint. Refinement chains remain
bounded by the native grid envelope, and ordinary fine-grid resume reconstructs
these prescriptions without repeating any coarse PDE solve.

Fine-grid cut quadrature may measure a different area. The SAME area policy is
therefore restored and the resulting fine geometry is independently solved.
This becomes the new comparison baseline at update zero. Cross-grid compliance
changes and baseline area restoration are never counted as improvement.

`constraints.refinement_origin` records the source receipt, its grid level and
accepted-update count, its measured endpoint, the transferred area, and the
largest nodal projection correction. That correction is a field-value change,
not a certified geometric distance. HTML and JSON reports show the same data.
The canonical fine source is ledger-derived from the coarse receipt through
real input/output edges; old sealed study operations are never rewritten.

Initial transfer/source/fine-solve cancellation returns no child study. The
original coarse state survives. Once a fine baseline exists, standard accepted-
state checkpointing and lifetime wall accounting apply. No hard interruption
latency is claimed inside assembly, prolongation, quadrature or ledger I/O.

## Verification boundary

```sh
cargo test --release -p fs-topols --lib refinement::projected
cargo test --release -p fs-cli --lib study::elasticity::continuation
cargo test --release -p fs-cli --test study_checkpoint_cli
```

The new tests exercise actual coarse/fine mechanics, load scaling, cancellation,
retained-source lineage, multi-level inherited pins, fine-grid continuation,
changed-input refusals and the public study/report commands. Authored tests
are not evidence of a passing build until executed. This remains numerical
2-D CutFEM optimization, not automatic adaptivity, a continuum error/stress
certificate, material calibration, physical validation or an optimum proof.
