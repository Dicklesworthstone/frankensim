# Fleet reference baseline candidates (bead c40j)

These files are operator-promoted machine-axis **candidates**, not authority-
admitted citable baselines. Each plain JSONL row is a `promote_baseline` record
from at least three mutually agreeing quiet probes, with a named operator,
justification, promotion day, and age policy. Plain rows can ground diagnostic
comparison, but they are structurally report-only even when every axis is
inside its band.

Perf lanes read a candidate or attested store via
`FRANKENSIM_BASELINE_STORE=<file>` plus `FRANKENSIM_FIRMWARE_ID=<id>`; the
`roofline` CLI uses `--baseline <file>` plus `--firmware <id>`. A positive
authority-admitted run additionally needs the configured promotion-authority
policy and retained source-receipt inputs. Perf lanes name them
`FRANKENSIM_PROMOTION_AUTHORITY_POLICY=<tsv-file>` and
`FRANKENSIM_RETAINED_SOURCE_RECEIPTS=<sorted-lowerhex-lines-file>`; the CLI uses
`--authority-policy <tsv-file>` and `--retained-receipts
<sorted-lowerhex-lines-file>`. If any authority input is absent, malformed,
denied, revoked, or inconsistent, measurement may continue but the result
remains explicitly report-only.

Bootstrap/update a candidate (re-promotion is the only update path):

    roofline promote --store perf-baselines/<machine>.jsonl \
      --firmware "<os-kernel-id>" --operator "<who>" \
      --justification "<why>" [--probes 3] [--age-days 90]

A loaded host REFUSES promotion (drift bands) — measure quiet. Promotion emits
a plain candidate; it does not mint an attestation or authority decision.

Promotion reads the existing store through the 1 MiB parser bound, serializes
same-store writers with an OS file lock, durably writes a same-directory staging
generation, and atomically renames it over the prior store. A crash before the
rename leaves the prior candidate store intact; a later promotion safely
overwrites the staging generation.

TRUST BOUNDARY: the committed files are operator-trusted and tamper-evident
(content-hashed), but contain neither an attestation nor an atomic authority-
policy decision. Do not describe a run using only one of these files as
`citation_eligible` or citable. An attested envelope is admitted only when the
configured verifier accepts it and every named source hash is declared present
in the protected inventory; the exact decision snapshot must travel with the
run's retained output.

The retained-source file is a protected operator inventory declaration. This
entrypoint checks canonical hash membership only; it does not fetch, rehash, or
independently prove that the named source bytes remain retrievable.

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
