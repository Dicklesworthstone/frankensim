# Resume a fixed-grid minimum-material search

```sh
# Create a NEW directory; no existing run or checkpoint is overwritten.
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example minimum_material_sdf3 -- 10 250000 0.0012 0.003 \
  --checkpoint-dir run-a

# Use the exact last complete checkpoint path printed by the first command.
# The first positional value is ADDITIONAL accepted updates, not lifetime steps.
CHECKPOINT=run-a/REPLACE_WITH_PRINTED_FILE.fsmcp
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example minimum_material_sdf3 -- 20 250000 0.0012 0.003 \
  --resume "$CHECKPOINT" --checkpoint-dir run-b
```

Use the repository's remote execution lane where available. The original
Krylov budget and reaction bounds must be repeated exactly; restart cannot
increase them. Every rebuild spends the original quadrature box/point limits,
and accepted-endpoint/feasible-incumbent reconstruction spends the remaining
original solver and optimizer-evaluation allowances. The unchanged per-solve
limit still applies. Save/resume is explicit and fixed-grid: do not combine
it with adaptive, binary, estimate or refit modes.

The actual executable bytes are hashed using the existing fs-blake3 owner.
A changed binary, problem bounds, payload checksum or version is refused before
physics. The payload retains exact floating-point bits, original observation
signatures, multiplier/penalty/spectral schedule, consumed optimizer work,
accepted history and a separate least-volume feasible incumbent. The decoder
caps payload bytes at 64 MiB and cumulative decoded entries at four million,
rejecting unrepresentable counts, nonfinite values, truncation and trailing data.
A checksum detects accidental corruption; it is NOT a signature or authority
against an adversary who can rewrite the checkpoint. Use trusted owned files.

Restoration re-solves one accepted endpoint and a second endpoint only when the
feasible incumbent differs. Every observation, signed original limit, physical
density, stiffness and derivative sample is checked. Previous line searches
are not replayed; no stored displacement, adjoint or preconditioner is trusted.
Failed restoration reserves its endpoint evaluations on the mutable checkpoint
and records actual spent solver work in the new directory before returning.
Use that updated file for a retry rather than resetting accounting to an old copy.

Each file uses create-new semantics and is flushed with file sync_all. Baseline,
accepted steps, late failed steps and the final state get separate numbered
files, with paths printed on stderr. The source files remain intact. An abrupt
process death or I/O failure can leave the latest file incomplete; choose an
earlier complete file, whose checksum and endpoints still need admission.
Work since the last complete write cannot be recovered after a crash. No
power-loss atomicity of the parent directory or mid-Krylov recovery is claimed.

Physical output selects the least-volume feasible incumbent when available;
the saved optimizer continues from its latest accepted iterate, not that older
incumbent. Its stop and KKT diagnostics belong to that iterate. An infeasible
run returns nonzero after saving its resumable state. Fixed-grid feasibility
is still not continuum certification, unique recovery or global optimality.

Programmatic users call MinimumMaterialStudy3::checkpoint/restore and
MinimumMaterialCheckpoint3::to_bytes/from_bytes. Bind arbitrary callbacks,
geometry, source identity, executable and original budgets externally. Exact
endpoint checks cannot establish callback equivalence throughout design space.
The command's immutable fixed fixture plus executable digest supplies that
binding for this example; the library byte codec alone does not.

Tests are included in response_limits_sdf3 and the minimum_material_sdf3
example. They exercise same-trajectory reconstruction, original-limit changes,
failed-restore charges, bounded decoding, actual file roundtrips and preserved
source files. Native execution must be checked separately from source checks.
The example adds fs-blake3 only as an existing-workspace dev-dependency; no
production dependency or external version changes. Cargo.lock was not
regenerated in the editing environment because Cargo is unavailable; normal
Cargo resolution must record that dependency edge before --locked validation.
