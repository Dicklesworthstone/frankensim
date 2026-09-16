# Native elasticity study continuation

The existing native `study` command continues an elasticity study from its last
accepted level-set geometry rather than rerunning every earlier prefix.
It uses `fs_topols::OptimizeCheckpoint` and the existing `.fsim` source,
`study-design`, `study-iterations`, and sealed ledger receipt. No extra
checkpoint-file format or alternate optimizer is introduced.

```bash
frankensim --json study examples/marquee/bracket-2d.fsim study.db --budget 2
```

The returned `run_id` identifies the retained state. Continue it with:

```bash
frankensim --json study --resume study-<receipt-hash> study.db --budget 2
```

`--budget N` allows at most N additional accepted geometry updates for this
invocation. The source's total step count does not reset. The cumulative wall
charge also does not reset: a study whose declared lifetime wall budget is
exhausted cannot obtain more physics merely by resuming it.

## What resumes

The existing design artifact retains every phi node as sixteen hexadecimal
IEEE-754 bits. The final iteration row retains the augmented-Lagrange
multiplier, and its ordinal determines the next global update. Resume checks
the lattice and node count, finite values, contiguous iteration ordinals,
exact row trace hash, and agreement between the final row and geometry.

The multiplier and global ordinal are restored with the geometry. Resetting
just one of them would change the volume-control trajectory or restart the
hole-nucleation schedule. The retained geometry is canonically re-solved before
the next shape sensitivity is used; an old displacement field is not reused on
a new geometry.

New receipts carry `continuation.version=1`, an exact executable fingerprint,
and counts of new updates and legacy prefix updates replayed in this invocation.
Direct continuation requires the same executable and deterministic runtime
profile. This is not a cross-ISA replay claim. A changed executable refuses
without replacing the old artifacts. Older receipts without this binding get
one full-prefix reconstruction and exact trace/geometry verification before
extension; subsequent native continuations do not repeat that migration.

## Work and recovery

The former driver ran prefixes of length 1, 2, ..., N. With an initial solve
and one candidate solve per prefix update, that costs N(N+3)/2 elasticity
solves. The checkpoint implementation instead canonically solves the current
and candidate geometry for each update: 2N solves. Thus an eight-update fresh
run changes from 44 to 16 solves, and a thirty-update run from 495 to 60.
These are code-path operation counts, not measured wall-time speedups. A single
uninterrupted library optimizer remains cheaper at N+1 solves; the native
checkpoint path deliberately pays one canonical reconstruction per update.

The seed and each accepted update are committed through the existing ledger
before another update starts. An actual solver or persistence failure is still
an error, with the last committed run pointer in the diagnostic. A trial that
failed its physics solve is never substituted for the accepted geometry.
Cancellation and time checks occur between updates, not inside an individual
CutFEM solve. No fine-grained interruption mechanism is added here.

A final accepted numerical state is committed before final publication. If
final publication fails, that running checkpoint can be finalized with zero
additional geometry updates, subject to the same admission and budget checks.
A completed receipt is returned unchanged on a redundant resume.

`completed` means all declared geometry updates were performed; it is not a
proof of optimality, volume feasibility, monotone descent, stress compliance,
physical validation, or continuum error control. The original optimizer's
numerical scope is unchanged.

Report and package commands still project retained artifacts without solving:

```bash
frankensim --json report study-<receipt-hash> study.db
frankensim --json package study-<receipt-hash> study.db
```

Exact replay concerns the design and iteration artifact bytes. Receipt hashes
and wall charges can differ between chunked and unchunked invocations because
their execution history and predecessor chain differ.

## Focused regressions

```bash
cargo test -p fs-cli --lib study::elasticity
cargo test -p fs-cli --test study_checkpoint_cli
cargo test -p fs-topols --lib checkpoint
```

The native tests compare a real CutFEM/level-set trajectory against the original
fixed-run optimizer, cross a global nucleation ordinal, restore a cancelled
seed, finalize the last accepted state, and refuse corrupted data, changed
executables and exhausted lifetime budgets. Actual-binary tests compare
on-disk design and iteration bytes across full and 1+2-update runs, then export
the retained results. The authoring environment did not have a Rust toolchain;
these test sources are not a claim of successful execution.
