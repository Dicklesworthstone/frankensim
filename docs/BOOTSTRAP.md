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

The xtask form requires a workspace that already resolves. A CLEAN CLONE
cannot build xtask (Cargo resolves the fixed relative sibling path
dependencies first), so the clean-machine entry point is the standalone,
zero-dependency tool at `tools/bootstrap` (bead 1t8i) — deliberately NOT a
workspace member, so it builds alone:

    cargo run --manifest-path tools/bootstrap/Cargo.toml            # fetch + verify
    cargo run --manifest-path tools/bootstrap/Cargo.toml -- --offline
    cargo run --manifest-path tools/bootstrap/Cargo.toml -- --from <mirror-base>

It reads `constellation.lock`, clones each missing sibling into the
workspace parent (checkout DETACHED at the pin, then the same pinned-head
and clean-tree verification applied to an existing sibling; no
branches/worktrees created anywhere), verifies existing siblings (pinned
head + clean tree; drift and dirt are refusals, never silent
substitutions — a case-folding checkout collision surfaces as a dirty
tree and refuses), and writes `constellation-bootstrap.json` provenance
beside the siblings. The sibling layout itself is the reproducible Cargo
configuration: no config files are generated or mutated. Idempotent:
re-runs verify. Hermetic offline-cache replay drills (clean-machine
clone from a local bare mirror, idempotent replay, existing- and
newly-cloned-tree dirt, drift, offline-missing, and malformed CLI
refusals) live in `tools/bootstrap/tests/replay.rs`. A bare `--root` or
`--from`, including one immediately followed by another option, is a
structured admission failure. The xtask command remains the in-workspace
verifier once the workspace builds and applies the same post-clone clean-tree
check; its bare, empty, or option-followed `--dest` and `--from` operands also
refuse before any repository operation.

Standalone bootstrap behavior, per library:

- **Missing from the workspace parent**: clone the declared remote
  (transform-free,
  `core.autocrlf=false`), check out the locked revision DETACHED, and
  verify both the resulting head and clean tree. An unavailable revision,
  unreachable remote, or dirty post-checkout result is a structured failure.
- **Present in the workspace parent**: verify head == lock and the tree is clean.
  A wrong head or a dirty tree REFUSES — the bootstrap never silently
  substitutes a nearby working tree.
- **Case-collision artifacts**: paths differing only by case cannot
  coexist on case-insensitive filesystems (macOS/Windows), so such
  checkouts cannot satisfy the clean-tree contract. FrankenNumpy was
  deliberately re-pinned after the colliding corpus paths were renamed, and
  the current pin has clean-checkout evidence on case-insensitive macOS and
  case-sensitive Linux. Any future collision still surfaces as dirt and
  refuses rather than relabeling a changed byte.
- **Offline re-runs** succeed from a verified sibling set with no network.

Provenance: `constellation-bootstrap.json` is written into the workspace parent
with the lock hash and every library's head, remote, and state — the
build-provenance record the acceptance requires.

Re-locking (`lock-constellation`) is a DELIBERATE act: it re-records
live heads and remotes; `check-constellation` gates drift in CI.

`frankensim-1t8i` is closed with the standalone tool, hermetic local-mirror
clone/replay drills, and a real seven-sibling offline verification. The
recorded session did not perform a literal blank-host fetch from every public
remote, so do not cite it as evidence for remote availability or public-network
provisioning; that no-claim does not weaken the lock, checkout, or cleanliness
contracts exercised by the retained replay tests.
