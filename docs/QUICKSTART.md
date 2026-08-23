# FrankenSim Quickstart

Bootstrap to a first validated project in about 15 minutes on a warm
machine (measured: Apple M4 Pro, pinned nightly toolchain, siblings already
materialized; the first cold build takes longer because it compiles the
workspace). Every command below is executed continuously by
`scripts/ci/examples_freshness_e2e.sh` and
`scripts/ci/solve_stage_producers_e2e.sh`, so if a command here drifts from
what the product actually does, a lane breaks.

## 0. What you are validating

FrankenSim is fail-closed by design. `validate` reports every structural
finding about your project and never guesses. A green validate means the
file is a well-formed project under the frozen schema — not that physics
has been solved. Solve stages beyond the durable producer prefix are typed
gaps that refuse by name; that refusal is the product working, and this
quickstart ends by reading one deliberately.

## 1. Materialize the sibling constellation

The workspace path-depends on seven pinned sibling Franken repositories.
From a fresh checkout:

```bash
rustup toolchain install nightly
rustup component add rustfmt clippy --toolchain nightly
cargo run --manifest-path tools/bootstrap/Cargo.toml
```

The bootstrap is standalone and zero-dependency so Cargo never has to
resolve the missing path dependencies before it builds. See
`docs/BOOTSTRAP.md` for offline and mirror modes.

## 2. Validate your first project

```bash
cargo run -p fs-cli --bin frankensim -- --json validate examples/heated-plate/heated-plate.fsim
```

Expected: exit 0 and JSON with

- `"status":"ok"`,
- `"finding_count":0`,
- `"project_hash"` — the canonical hash of exactly these bytes,
- `authority` and `no_claim` strings stating what this verdict is and is
  not.

Edit any value (try changing `:duty 1.0` to `2.0`) and re-run: the validator
names the code (`project-duty-range`), what happened, and the fix. That
triple — code, what, fix — is the shape of every FrankenSim refusal.

## 3. Break it on purpose

```bash
cargo run -p fs-cli --bin frankensim -- --json validate examples/refusal-loop/broken.fsim
```

Expected: exit 4 (schema/semantic-refusal class) and a finding carrying
code `project-duty-range`. `broken.fsim` is the tracked reference project
with exactly one token changed; restoring `:duty 1.0` makes it
byte-for-byte `data/reference-project/cooling-reference.fsim`. The
freshness lane proves both directions continuously.

## 4. Run the solve producer prefix

```bash
LEDGER="$(mktemp -d)/plate.db"
cargo run -p fs-cli --bin frankensim -- --json import \
  data/reference-project/cooling-reference.fsim data/reference-project/plate.stl \
  "${LEDGER}" --unit m --max-hole-edges 0
cargo run -p fs-cli --bin frankensim -- solve \
  data/reference-project/cooling-reference.fsim "${LEDGER}" \
  --materials data/reference-project/aa6061.fsmcdpk
```

The import admits the quarantined STL into the ledger; the solve executes
the durable producer stages — import-verify, assign, material-resolve, and
an interval-certified flow-network operating point — and then refuses at
the conduction stage naming the bead that owns the gap. Report and package
refuse the same way. Nothing is written that claims more than happened.

## 5. Where to go next

| Path | What it teaches |
|------|-----------------|
| `examples/heated-plate/README.md` | Every mandatory schema section, minimally |
| `examples/refusal-loop/README.md` | The refusal/fix loop as a workflow |
| `examples/heatsink-fan/README.md` | The full cooling contract and the flow-network operating point |
| `examples/cooling-enclosure/README.md` | Reading receipts, budgets, and no-claims |
| `crates/fs-project/src/spec.rs` | The validated schema, section by section |
| Each crate's `CONTRACT.md` | Invariants, determinism class, and no-claim boundaries |
