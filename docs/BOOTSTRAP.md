# Constellation bootstrap (beads huq.17, 1t8i)

The workspace's Franken dependencies resolve through sibling paths
(`../asupersync`, `../franken_numpy`, …). `constellation.lock` (schema v2) is
owned by identity domain `org.frankensim.xtask.constellation-lock.v1` at
identity version `1`. It records each repository's semantic identity — version
and git head — plus its remote as TRANSPORT (content identity is the commit
hash, never the URL). The row-only `lock_hash` continues to cover
`(lib, version, git_head)`; schema/domain/version, paths, and remotes do not
silently enter that existing identity. A host whose sibling dependencies
already exist can verify or refresh the pinned sources from the lock:

The durable Rust writer is tracked separately under implementation-authority
domain `org.frankensim.xtask.constellation-lock-writer.v2` at authority version
`2`. That coupling epoch binds the live encoder and sink without changing or
serializing a new lock identity: the canonical lock envelope remains at the v1
domain/version above and its row-only hash is unchanged.

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

## Typed admission core (bead sj31i.50.1, first tranche)

`xtask/src/constellation_admission.rs` is the zero-dependency, Rust-2021-compatible
policy core shared by the workspace and standalone bootstrap binaries. Its
versioned `AdmissionContext` makes the command class, offline or pinned-fetch
authority, publication authority, Cx/cancellation/clock identities, deadline,
work, memory, process, file, output, network, retry budgets, path capabilities,
executable capabilities, and trust-anchor state explicit. Omitted values,
booleans, positive budgets, ambient configuration, and `Option` presence cannot
mint authority.

The state machine distinguishes `Diagnostic`, `Unanchored`, `Admitted`,
`Refused`, `Cancelled`, and `Indeterminate`. It transactionally refuses illegal
phase changes and over-budget charges before mutating consumption or the event
ledger. Monotonic clock observations cannot regress, explicit cancellation polls
are bounded by the Cx declaration, and event capacity always reserves a path to
an honest terminal record. Publication requires preflight, an exact stability
recheck, explicit publication authority, evidence that work is quiescent, a
second exact stability recheck with a conditional sink fence, single-use receipt
authorization, then either successful finalization or an indeterminate failure
record. Cancellation requires an identified liveness observation, request,
identified deterministic drain, and identified explicit finalization; failed or
uncertain cleanup is `Indeterminate`, never clean refusal or cancellation. Retry
preserves the request identity and cumulative spend, advances a bounded attempt
ordinal, binds a fresh attempt-scoped child Cx/cancellation source, and cannot
reuse an indeterminate attempt's authority.

Canonical binary records retain only schema/domain/version, opaque fixed-size
identities, policy, budgets, and the complete deterministic transition history.
State, consumption, attempt/Cx disposition, event ordering, terminal rule, and
that rule's ranked closed remedies are replay-derived rather than trusted as
independent fields. Records never retain paths, URLs, credentials, environment
values, or runtime handles. Decoding replays the complete history and returns an
inert `RecordedAdmission`; it cannot reconstruct a live machine or capability.
The canonical audit envelope is bounded separately and does not spend the
operation's external stdout/stderr/provenance output budget.

This first tranche establishes the shared core, canonical receipt, and
G0/G4/G5 transition tests. It does **not** yet claim that every existing
bootstrap, shell, snapshot, verify-only, DSR, or RCH effect is gated by the new
machine. Those adapters remain part of the open bead, and the existing entry
points retain their documented behavior until each effect boundary is migrated
without inventing authority.

It reads `constellation.lock` through a 1 MiB bound and accepts only the
canonical v2 grammar, exact identity domain and integer identity version, exact
seven-library set, unique safe library identities, lowercase pinned heads,
canonical note, and recomputed row-only lock hash before deriving any destination
path. Missing, wrong, duplicate, unknown, or type-invalid identity metadata is a
lock parse failure before repository or destination access. It initializes each
missing destination in place with
local marker `frankensim.bootstrapIncomplete=true`, fetches the exact pin,
checks it out DETACHED, applies the same pinned-head and clean-tree verification
as an existing sibling, and clears the marker. A retry resumes only a clean
marked repository or a clean unmarked unborn repository whose `origin` exactly
matches the selected transport. Ordinary non-empty/non-git paths, directories
that cannot be enumerated as empty, wrong-origin unborn repositories, and
unmarked wrong-head repositories are refusals (no
branches or worktrees are created anywhere). Every Git subprocess sets
`GIT_NO_LAZY_FETCH=1`, so a partial/promisor checkout cannot turn an offline
verification or checkout into an implicit object fetch. Existing siblings are verified at
the pinned head with a clean tree; drift and dirt are refusals, never silent
substitutions. Before hashing, the verifier records every materialized logical
index prefix, requires every intermediate component to be an ordinary directory,
and refuses two distinct logical prefixes that resolve to one filesystem
identity on Unix. This catches case-folding and Unicode-normalization checkout
collapse there; tracked hard-link aliases are conservatively refused by the same
Unix identity rule. Windows instead requires the exact preserved directory-entry
spelling for every materialized prefix. Recursively,
for every initialized nested submodule at every
depth, clean also means its HEAD matches the containing repository's recorded
gitlink and its tracked and untracked worktree is clean. Every verifier
explicitly overrides repository configuration and `.gitmodules` `ignore`
policy at each repository boundary. This does not initialize an absent
submodule or add nested fields to `constellation.lock`; the outer commit already
binds each gitlink. Clean verification disables global excludes, independently
enumerates untracked files hidden by repository-local excludes at every depth,
rejects untracked `.gitignore` authorities plus `assume-unchanged`,
`skip-worktree`, and persisted FSMonitor index authority. Ordinary status/index
probes disable fsmonitor and the untracked cache. A separate raw primary-index
inspection supports index versions 2 and 3 and conservatively refuses any
`FSMN` extension without claiming to decode its per-entry dirty bitmap: Git's
fsmonitor-valid bit is in-memory state that a failed refresh can erase before
`ls-files -f` reports it. The raw reader authenticates the trailing SHA-1 or
SHA-256 checksum and conservatively validates Git's byte-level pathname grammar
before any worktree join. Case-insensitive `.git` components, symlink
`.gitmodules`, and every UTF-8 sequence Git documents as HFS-ignored refuse on
all platforms; Windows additionally refuses native separators, rooted or
drive/stream-qualified paths, and short-name `.git` or symlink `.gitmodules`
aliases. Malformed, truncated, or
checksum-invalid indexes, version 4, split-index `link`, sparse-index `sdir`,
every other unsupported required extension, and layouts whose object-hash
width or checksum boundary is ambiguous also refuse.
An absent primary index is accepted only when the parsed
Git inventory is empty. Clear FSMonitor or split-index state deliberately before
bootstrap verification rather than relying on a monitor refresh. Cleanliness
also compares every stage-zero regular-file and symlink
object ID against raw `hash-object --no-filters` bytes, so clean/smudge, EOL, or
working-tree-encoding transforms cannot redefine source cleanliness. Execute-bit
comparison is enabled only on platforms where filesystem mode is observable;
raw object identity and type checks remain mandatory everywhere. Every Git
subprocess discards inherited repository/worktree/index/object redirection and
standard-handle, attribute-source, and pathspec-mode controls, and disables
inherited trace destinations and replacement objects. It also isolates
global/system Git configuration,
templates, hooks, attributes, and credential helpers. Git transport admission
is default-deny: only `file`, `https`, and `ssh` are enabled, while `http`,
`git`, custom remote helpers, and all other protocols refuse before helper
execution. The shell entrypoint runs each embedded Python program in isolated
mode so `PYTHONPATH` cannot inject code before lock validation;
local or per-worktree include, filter, hook, fsmonitor, SSH/helper, remote-VCS,
bundle/packfile-URI, URL-rewrite, askpass, and other executable or transport
configuration authorities refuse before any filter-aware observation or
mutation. Reference-storage configuration is also refused, and fresh
initialization ignores inherited default object/ref formats plus
`GIT_REFERENCE_BACKEND`. Symlinked repository roots, linked intermediate index
components, local replace refs, and graft authorities
also refuse. Checkout uses `--no-overwrite-ignore`, so a marked resume cannot
replace owner bytes merely because the wrong commit's tracked `.gitignore`
hides them. Fetch disables automatic maintenance and submodule recursion;
checkout also disables submodule recursion. Initialized nested repositories are
recursively verified but never fetched or updated by the outer bootstrap.
Recursive admission reaches every initialized descendant before each parent
status probe forces gitlink visibility with `--ignore-submodules=none`; a
NUL-delimited non-content cached name diff independently detects staged gitlink
rewrites. After materialization, the tool stages and fsyncs a complete
provenance candidate, then performs two whole-constellation verification passes
and an exact lock-byte recheck before the atomic publication rename. It writes
`constellation-bootstrap.json` provenance (schema
`frankensim-constellation-bootstrap-v2`, identity domain
`org.frankensim.xtask.constellation-bootstrap-provenance.v3`, identity version
`3`) beside the siblings. The encoder uses a freshly reserved same-directory
staging file and atomically renames it only after the full document and final
publication barrier are complete. The staging handle remains open across that
long barrier. Immediately before rename, Unix publishers require the visible
ordinary single-link file to match the sealed device/inode and mutation-sensitive
metadata captured after fsync; Windows publishers apply the equivalent
volume-serial/file-index, link-count, attribute, size, and timestamp seal while
opening the visible entry without following a reparse point. A missing path,
symlink/reparse point, hard-link change, metadata mutation, or substituted file
therefore refuses publication.
A failed validation or staging write therefore does not truncate or replace a
verified provenance document; a receipt from an earlier successful invocation
may remain and must not be attributed to the failure. A failed staging or
rename operation deliberately retains its uniquely named same-directory temp
file for diagnosis; that file is never authoritative and cleanup requires an
explicit operator decision. The staged file is fsynced before rename and a
parent-directory fsync is attempted, but filesystems that reject directory
fsync retain an explicit crash-durability no-claim. The identity seal closes
the long validation-window substitution risk on Unix and Windows. Safe `std`
still has no handle-relative rename primitive, so there is a narrow final
identity-check-to-path-rename race and no claim against an untrusted concurrent
writer that can mutate the destination directory. Cooperating concurrent
bootstrap invocations use distinct staging names and remain complete; they are
not ordered, and the last successful rename wins. Targets that are neither
Unix nor Windows fail a support preflight before any sibling repository
mutation because safe `std` cannot prove the staging file's object identity.
Non-UTF-8 provenance destination paths refuse before repository mutation rather
than being serialized lossily. The
sibling layout itself is the reproducible Cargo
configuration: no config files are generated or mutated. Idempotent:
re-runs verify. Hermetic replay drills cover producer-to-consumer parsing of the
tracked xtask lock; missing/wrong/duplicate/unknown/type-invalid identity-field
refusal; clean-machine fetch from a local bare mirror; idempotent offline replay;
interrupted marked and exact-origin unborn resume; unsafe-destination and lock
tamper/path-traversal refusal; ordinary and ignore/index-hidden dirt;
initialized nested-submodule tracked/untracked dirt, recursively nested dirt,
repository-local excludes and hidden index flags inside submodules, gitlink and
HEAD drift, untracked ignore-policy and clean-filter concealment,
`ignore=dirty`/`ignore=all` concealment; newly-fetched-tree dirt; drift;
offline-missing; malformed CLI refusal; checksum-corrupted and forbidden-path
primary indexes plus persisted FSMonitor, v4, and split-index layouts; and the
shell checkout's synthetic checkout, verify, and stable snapshot modes. They live
in `tools/bootstrap/tests/replay.rs`. The standalone `--root` must name an
ordinary directory (not a symlink or reparse point) and is canonicalized before
the lock is read or the sibling destination is derived. A bare `--root` or
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
  case-sensitive Linux. Every logical index prefix is bound to its physical
  filesystem identity on Unix (or exact preserved spelling on Windows); a future
  case/normalization collision refuses before any collapsed or externally
  redirected byte can be hashed. Unix tracked hard-link aliases are also refused;
  no broader Windows hard-link claim is made.
- **Offline re-runs** succeed from a verified sibling set with no network.

Provenance: `constellation-bootstrap.json` is written into the workspace parent
with the exact schema/domain/version header, lock hash, and every library's head,
canonical lock remote, selected transport/mirror, whether that transport was
actually used, and terminal state. The standalone and in-workspace Rust
implementations use one canonical encoder and emit the same top-level and row
shapes. The retained top-level schema string remains v2 while the identity
epoch is v3. A verified or
offline-cache row records `transport_used: false`; a fresh clone records the exact
selected transport with `transport_used: true`.

Re-locking (`lock-constellation`) is a DELIBERATE act: it re-records
live heads and remotes; `check-constellation` gates drift in CI. The writer
fsyncs a uniquely staged same-directory document and atomically replaces the
lock, so a staging or rename failure cannot truncate the previous valid lock.

`frankensim-1t8i` is closed with the standalone tool, hermetic local-mirror
clone/replay drills, and a real seven-sibling offline verification. The
closure also records a literal blank-box run on `ts1`: a directory containing
only a fresh public clone of FrankenSim fetched all seven siblings from their
locked remotes, completed a locked Cargo check, then completed idempotent and
offline replays with retained provenance. That run proves the bootstrap path
and pins as observed in that session; it is not a standing guarantee that every
public remote will remain available later.
