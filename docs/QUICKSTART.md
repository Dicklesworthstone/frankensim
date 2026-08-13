# Quickstart: first honest solve

This is the shortest path from a fresh checkout to the product's current
end of the road — a four-stage ledgered solve that stops, by name, at the
first stage that does not exist yet. Everything here is executed
continuously by `scripts/ci/solve_stage_producers_e2e.sh`; the walkthrough
version with output-reading guidance is `examples/cooling-enclosure/`.

```bash
# 0. Toolchain + sibling constellation (see docs/BOOTSTRAP.md for modes).
rustup toolchain install nightly
cargo run --manifest-path tools/bootstrap/Cargo.toml

# 1. Validate the tracked reference project (canonical hash, zero findings).
cargo run -p fs-cli --bin frankensim -- --json validate \
  data/reference-project/cooling-reference.fsim

# 2. Import its geometry into a fresh FrankenSQLite ledger.
cargo run -p fs-cli --bin frankensim -- --json import \
  data/reference-project/cooling-reference.fsim \
  data/reference-project/plate.stl \
  /tmp/quickstart.db --unit m --max-hole-edges 0

# 3. Solve: four stages execute with durable receipts, then the pipeline
#    fails closed at conduction, naming the owning bead.
cargo run -p fs-cli --bin frankensim -- --json solve \
  data/reference-project/cooling-reference.fsim \
  /tmp/quickstart.db \
  --materials data/reference-project/aa6061.fsmcdpk
```

What "done" looks like today: stderr shows `import-verify`, `assign`,
`material-resolve`, and `flow-network` completing (the last with an
interval-certified operating point retained in the ledger); stdout reports
`"status":"unavailable"`, `"stage":"conduction"`,
`"dependency":"frankensim-s93ej"`. That is the product refusing to
impersonate a finished simulator — the refusal is the feature. To see the
validation-refusal loop on purpose, run step 1 against
`examples/refusal-loop/broken.fsim`.

Timing honesty: no time-to-first-result number is quoted here because the
dominant cost is the first cold workspace build, which varies by an order
of magnitude across machines and cache states. The three commands
themselves are seconds once built; measure your own cold path if the
metric matters to you.
