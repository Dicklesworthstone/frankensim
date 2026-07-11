# Fleet reference baseline stores (bead c40j)

Governed machine-axis baselines for citable roofline gates (dfh3
design): each JSONL row is a `promote_baseline` record — ≥3 mutually
agreeing quiet probes, named operator, justification, promotion day,
age policy. Gates consume via `FRANKENSIM_BASELINE_STORE=<file>` +
`FRANKENSIM_FIRMWARE_ID=<id>`, or `roofline --baseline <file>`.

Bootstrap/update a machine (re-promotion is the only update path):

    roofline promote --store perf-baselines/<machine>.jsonl \
      --firmware "<os-kernel-id>" --operator "<who>" \
      --justification "<why>" [--probes 3] [--age-days 90]

A loaded host REFUSES promotion (drift bands) — measure quiet.

TRUST BOUNDARY: these stores are operator-trusted and tamper-evident
(content-hashed), NOT independently verified — promotion-authority
signatures are bead fz2.7's layer. Do not present gate results as
third-party-verifiable until that lands.

DEPENDENCY-GRAPH RECEIPTS (bead fz2.6): citable perf builds must also
carry the resolved dependency+feature receipt so tune rows cannot cross
binaries whose dependency codegen differs. Before building the lane:

    export FRANKENSIM_DEPGRAPH_RECEIPT="$(cargo run -p xtask -- \
      depgraph-receipt -- <same cargo selection as the build>)"

and verify after the run with `... depgraph-receipt --verify -- <same
selection>`. Interactive workspace builds instead carry the explicit
`FRANKENSIM_DEPGRAPH_SALT` equivalence class from .cargo/config.toml;
builds with neither fail closed in crates/fs-la/build.rs.
