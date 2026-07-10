# Constellation bootstrap (beads huq.17, 1t8i)

The workspace's Franken dependencies resolve through sibling paths
(`../asupersync`, `../franken_numpy`, …). `constellation.lock`
(schema v2) records each repository's semantic identity — version,
git head — plus its remote as TRANSPORT (content identity is the
commit hash, never the URL). A host whose sibling dependencies already exist
can verify or refresh the pinned sources from the lock:

    cargo run -p xtask -- bootstrap-constellation            # workspace parent only
    cargo run -p xtask -- bootstrap-constellation --offline  # verify-only, no network
    cargo run -p xtask -- bootstrap-constellation --from B   # air-gapped mirror base

This is **not yet a clean-clone bootstrap**. Cargo resolves every workspace
manifest and its path dependencies before it can build `xtask`; therefore the
command cannot start when the siblings are absent. The manifests also contain
fixed relative sibling paths, so an arbitrary source cache cannot satisfy the
workspace. Bead `frankensim-1t8i` tracks extraction of a zero-dependency
pre-Cargo entry point and strict lock parsing.

Behavior of an already-built `xtask` binary, per library (the `cargo run`
commands above can exercise only the verification path while every required
sibling resolves):

- **Missing from the workspace parent**: a previously built binary can clone
  the declared remote (transform-free,
  `core.autocrlf=false`), check out the locked revision DETACHED, and
  verify the resulting head equals the lock. An unavailable revision or
  unreachable remote is a structured failure.
- **Present in the workspace parent**: verify head == lock and the tree is clean.
  A wrong head or a dirty tree REFUSES — the bootstrap never silently
  substitutes a nearby working tree.
- **Case-collision artifacts**: paths differing only by case cannot
  coexist on case-insensitive filesystems (macOS/Windows), so such
  checkouts cannot satisfy the clean-tree contract. The currently pinned
  FrankenNumpy revision contains a `seed_M`/`seed_m` collision and therefore
  cannot be verified cleanly on a default macOS filesystem. A case-safe
  upstream commit and deliberate lock update are the current clean-macOS
  blocker; the bootstrap refuses rather than relabeling a changed byte.
- **Offline re-runs** succeed from a verified sibling set with no network.

Provenance: `constellation-bootstrap.json` is written into the workspace parent
with the lock hash and every library's head, remote, and state — the
build-provenance record the acceptance requires.

Re-locking (`lock-constellation`) is a DELIBERATE act: it re-records
live heads and remotes; `check-constellation` gates drift in CI.

Until `frankensim-1t8i` closes, provisioning a genuinely clean host requires a
separately prepared, case-sensitive constellation checkout at the exact sibling
paths. Do not cite the Cargo command above as clean-machine proof.
