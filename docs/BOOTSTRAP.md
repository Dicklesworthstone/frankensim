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

It reads `constellation.lock` through a 1 MiB bound and accepts only the
canonical v2 grammar, exact seven-library set, unique safe library identities,
lowercase pinned heads, canonical note, and recomputed lock hash before deriving
any destination path. It initializes each missing destination in place with
local marker `frankensim.bootstrapIncomplete=true`, fetches the exact pin,
checks it out DETACHED, applies the same pinned-head and clean-tree verification
as an existing sibling, and clears the marker. A retry resumes only a clean
marked repository or a clean unmarked unborn repository whose `origin` exactly
matches the selected transport. Ordinary non-empty/non-git paths, wrong-origin
unborn repositories, and unmarked wrong-head repositories are refusals (no
branches or worktrees are created anywhere). Existing siblings are verified at
the pinned head with a clean tree; drift and dirt are refusals, never silent
substitutions — a case-folding checkout collision surfaces as a dirty
tree and refuses. Clean verification disables global excludes, independently
enumerates untracked files hidden by repository-local excludes, and rejects
`assume-unchanged`/`skip-worktree` index flags, so local Git configuration
cannot conceal source drift. The tool then writes `constellation-bootstrap.json` provenance
(schema `frankensim-constellation-bootstrap-v2`) beside the siblings. The
sibling layout itself is the reproducible Cargo
configuration: no config files are generated or mutated. Idempotent:
re-runs verify. Hermetic offline-cache replay drills (clean-machine
fetch from a local bare mirror, idempotent replay, interrupted marked and
exact-origin unborn resume, unsafe-destination and lock tamper/path-traversal
refusal, ordinary and ignore/index-hidden dirt, newly-fetched-tree dirt, drift,
offline-missing, and malformed CLI
refusals) live in `tools/bootstrap/tests/replay.rs`. A bare `--root` or
`--from`, including one immediately followed by another option, is a
structured admission failure. The xtask command remains the in-workspace
verifier once the workspace builds and applies the same post-clone clean-tree
check; its bare, empty, or option-followed `--dest` and `--from` operands also
refuse before any repository operation.

Standalone bootstrap behavior, per library:

- **Missing from the workspace parent**: initialize and mark the destination,
  fetch the declared remote (transform-free,
  `core.autocrlf=false`), check out the locked revision DETACHED, and
  verify both the resulting head and clean tree before clearing the marker. An unavailable revision,
  unreachable remote, or dirty post-checkout result is a structured failure.
- **Interrupted destination**: resume only when clean and marked, or when clean,
  unborn, unmarked, and already bound to the exact selected origin. Success
  checks the exact pin and clears the marker; unsafe partial destinations remain
  untouched.
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
with the lock hash and every library's head, canonical lock remote, selected
transport/mirror, whether that transport was actually used, and terminal state.
The two Rust bootstrap implementations emit the same v2 row shape. A verified or
offline-cache row records `transport_used: false`; a fresh clone records the
exact selected transport with `transport_used: true`.

Re-locking (`lock-constellation`) is a DELIBERATE act: it re-records
live heads and remotes; `check-constellation` gates drift in CI.

`frankensim-1t8i` is closed with the standalone tool, hermetic local-mirror
clone/replay drills, and a real seven-sibling offline verification. The
recorded session did not perform a literal blank-host fetch from every public
remote, so do not cite it as evidence for remote availability or public-network
provisioning; that no-claim does not weaken the lock, checkout, or cleanliness
contracts exercised by the retained replay tests.
