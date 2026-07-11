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

Promotion reads the existing store through the 1 MiB parser bound, serializes
same-store writers with an OS file lock, durably writes a same-directory staging
generation, and atomically renames it over the prior store. A crash before the
rename leaves the prior authority store intact; a later promotion safely
overwrites the staging generation.

TRUST BOUNDARY: these stores are operator-trusted and tamper-evident
(content-hashed), NOT independently verified — promotion-authority
signatures are bead fz2.7's layer. Do not present gate results as
third-party-verifiable until that lands.

DEPENDENCY-GRAPH RECEIPTS (bead fz2.6): evidence-bearing perf builds must also
carry the resolved dependency+feature receipt so tune rows cannot cross
binaries whose dependency codegen differs. Select exactly one production root
that reaches fs-la through normal/build edges. Before building, for example:

    export FRANKENSIM_DEPGRAPH_RECEIPT="$(cargo run -p xtask -- \
      depgraph-receipt -- --package fs-roofline \
      --target x86_64-unknown-linux-gnu)"

Use the identical root, target, and root-feature flags for verification:

    cargo run -p xtask -- depgraph-receipt --verify -- \
      --package fs-roofline --target x86_64-unknown-linux-gnu

Receipt v1 accepts only `--package`, `--target`, `--features`/`-F`,
`--all-features`, and `--no-default-features`. It refuses workspace,
test/dev/all-target, target-kind, and profile selection rather than silently
approximating them. The build profile remains separately bound by fs-la's full
build fingerprint, so pass `--release`/`--profile` only to the actual Cargo
build, not to the receipt command.

TRUST BOUNDARY: xtask runs locked metadata + `normal,build` tree observations
through one content-addressed/versioned Cargo executable. Tree rows must map to
unique structured metadata package/source IDs. Every local path package in the
fs-la closure is content-addressed over its bounded package-root file tree;
the package-root `.git` and Cargo `target` directories are excluded, while
nested source/data directories with those names remain hashed. Escaping,
non-regular symlinks, unreadable trees, and trees beyond the explicit file,
directory, depth, byte, or manifest bounds fail closed. The canonical receipt
is capped at 1 MiB. `build.rs` strictly validates and binds both its exact bytes
and domain-separated digest, writes the bytes under `OUT_DIR`, and compiles them
with `include_str!` rather than a full-receipt rustc environment variable.
Cargo still cannot prove to a dependency build script that this operator-
supplied receipt is the invoking unit graph. Dynamic build-script inputs from
outside a package root, the environment, network, or generated output require
an explicit `FRANKENSIM_GEMM_CODEGEN_ID` and retained operator protocol. Roots can
inspect/store it through `fs_session::gemm_tune_build_evidence()`. It is
operator-observed evidence, not a signature or independent verification.
Consumers rehash stored bytes with the re-exported
`fs_session::GEMM_DEPGRAPH_RECEIPT_DOMAIN`.
Interactive workspace builds instead carry the explicit
`FRANKENSIM_DEPGRAPH_SALT` development equivalence class from
`.cargo/config.toml`; that salt is never verified graph evidence. Builds with
neither class fail closed in `crates/fs-la/build.rs`.
