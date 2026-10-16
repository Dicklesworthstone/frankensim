# Retained supplier CAD corpus

This directory is the external-reality input for Bead
`frankensim-extreal-program-f85xj.11.6`. It retains exact files from real
supplier-maintained CAD repositories; it is not a set of FrankenSim-authored
fixtures.

## Authority

`corpus-v1.tsv` binds every local file to:

- an upstream source kind, origin, revision, path, and object identity;
- a BLAKE3 identity of the retained local bytes;
- an SPDX license or public-domain dedication and its upstream record;
- a source-quality tier and expected result governed by the same review state.

Git sources use a 40-hex commit and Git blob SHA-1. Versioned HTTPS snapshots
use a provider revision and SHA-256. In both cases the independently computed
BLAKE3 binds the exact bytes exercised by the scorecard.

The first admitted population intentionally leaves every annotation in
`proposed` state. Agent observations are useful review material, but they are
not human review. The standing scorecard therefore remains fail-closed until an
identified reviewer inspects the per-file receipts, changes each accepted row
to `human-locked`, records a date and positive annotation revision, and explains
any tier or golden update. Changing annotations merely to make a regression
pass is forbidden.

The current files come from Adafruit Industries' official
`Adafruit_CAD_Parts` repository at the exact commit recorded in the manifest.
The pinned upstream `LICENSE` is retained as `LICENSE-Adafruit-MIT.txt`.

The scan-derived stratum contains a Smithsonian 3D Digitization PLY snapshot of
USNM 174698, Cuneiform3 Right. The Smithsonian object record identifies the
media as public domain/CC0; the manifest binds its 3D package revision, resource
path, SHA-256, and local BLAKE3. The population therefore covers clean
parametric exports, tessellation-heavy exports, scan-derived or repaired
meshes, and real files retained as known-broken regression cases without
transforming a FrankenSim-authored fixture into artificial evidence.

## Reproducible acquisition

Run `bash data/cad-import-corpus/fetch-sources.sh` from any directory in the
checkout. The script:

1. refuses to overwrite an existing destination;
2. verifies the remote bytes against the pinned Git blob or HTTPS-snapshot
   SHA-256 before writing;
3. downloads only the exact Git commit or provider resource;
4. verifies the retained file again after writing.

The script never deletes files. A partial or mismatched destination is left in
place and reported so a human can inspect it before authorizing cleanup.

## Growth policy

Every permission-cleared user-reported import failure is a corpus candidate.
Admission requires a new stable case ID, exact byte and upstream identities,
license or written permission, a quality tier, and an initially proposed
annotation. A human then reviews the observed receipt and locks the baseline.
Privacy-sensitive or proprietary inputs must not enter this public directory;
use a separately governed encrypted/private corpus with an evidence-only
summary instead.

`scorecard-summary-v1.json` is the compact tracked projection generated from
the real full sweep. It binds the exact manifest BLAKE3 identity and the full
per-file scorecard identity, then records two deliberately distinct views:

- `population` reports what the importer observed across all retained files.
- `reviewed` counts only files whose tier and expected outcome have independent
  human-locked annotation authority.

Only `reviewed` may supply dashboard rate denominators. Until a human locks at
least one annotation, clean/repaired/refused rates remain `NO-DATA`; proposed
annotations cannot manufacture a favorable metric. The locked-annotation
mismatch count remains a separate regression signal. Do not edit this artifact
by hand: regenerate it from the sweep whenever the manifest, retained bytes,
import implementation, or annotation authority changes.

## Human annotation review runbook

The reviewer must judge the retained bytes, source provenance, importer receipt,
quality tier, and expected outcome independently. The importer's observed
outcome is evidence for that judgment, not authority to approve itself. An AI
agent may prepare receipts or propose manifest text, but it must not identify
itself as the human reviewer or change `review_state` to `human-locked`.

From a stable checkout, replay the deterministic sweep with its canonical
machine-readable output visible:

```bash
RCH_REQUIRE_REMOTE=1 rch exec -- env \
  CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" \
  cargo test -j 2 --locked -p fs-io --test supplier_corpus -- --nocapture
```

The test prints one `supplier_corpus_scorecard=` line containing the full
per-file receipts and one `supplier_corpus_summary=` line containing the exact
compact projection. The output is emitted before the assertions that pin the
current annotation phase, so it remains available when a deliberate authority
transition makes those assertions or the tracked summary stale.

For every row the reviewer accepts:

1. confirm or correct `quality_tier`, `expected_outcome`, and
   `expected_detail` from the evidence;
2. set `review_state` to `human-locked`;
3. record an identified human in `reviewer`, a calendar-valid `YYYY-MM-DD` in
   `reviewed_at`, and a positive `annotation_revision`;
4. replace the proposal note with a justification that explains the reviewed
   judgment and any changed tier or golden outcome.

After review, replay the sweep. Update the phase-specific assertions in
`crates/fs-io/tests/supplier_corpus.rs` only to describe the independently
reviewed state, and replace `scorecard-summary-v1.json` only with the exact JSON
from the emitted `supplier_corpus_summary=` line. Replay until the focused test
passes, then regenerate and check the dashboard consumer:

```bash
cargo run --locked -p xtask -- generate-program-metrics
cargo run --locked -p xtask -- check-program-metrics
```

Review the manifest, focused-test assertion changes, compact summary, and
dashboard diff together. A passing replay cannot substitute for identified
human judgment, and a dashboard update cannot broaden this corpus's authority.

This retained population measures only itself. Its clean, repaired, and refused
rates are not universal supplier-CAD success probabilities, and its sampled
intersection census is diagnostic evidence rather than a certificate of
self-intersection freedom.
