# Constellation bootstrap (bead huq.17)

The workspace's Franken dependencies resolve through sibling paths
(`../asupersync`, `../franken_numpy`, …). `constellation.lock`
(schema v2) records each repository's semantic identity — version,
git head — plus its remote as TRANSPORT (content identity is the
commit hash, never the URL). A clean host obtains the pinned sources
FROM the lock:

    cargo run -p xtask -- bootstrap-constellation            # dest = workspace parent
    cargo run -p xtask -- bootstrap-constellation --dest D   # explicit source cache
    cargo run -p xtask -- bootstrap-constellation --offline  # verify-only, no network
    cargo run -p xtask -- bootstrap-constellation --from B   # air-gapped mirror base

Behavior, per library:

- **Missing from the cache**: clone the declared remote (transform-free,
  `core.autocrlf=false`), check out the locked revision DETACHED, and
  verify the resulting head equals the lock. An unavailable revision or
  unreachable remote is a structured failure.
- **Present in the cache**: verify head == lock and the tree is clean.
  A wrong head or a dirty tree REFUSES — the bootstrap never silently
  substitutes a nearby working tree (developer siblings remain an
  explicit fast path: align them deliberately or use a fresh `--dest`).
- **Case-collision artifacts**: paths differing only by case cannot
  coexist on case-insensitive filesystems (macOS/Windows), so such
  checkouts always show phantom dirt. When EVERY dirty path is provably
  one of a case-colliding indexed pair, the library is
  `verified-case-collision` (commit identity intact, noted loudly);
  any other dirt stays a hard failure.
- **Offline re-runs** succeed from a verified cache with no network.

Provenance: `constellation-bootstrap.json` is written into the dest
with the lock hash and every library's head, remote, and state — the
build-provenance record the acceptance requires.

Re-locking (`lock-constellation`) is a DELIBERATE act: it re-records
live heads and remotes; `check-constellation` gates drift in CI.
