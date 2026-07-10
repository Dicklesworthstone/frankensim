# fs-blake3 CONTRACT

## Purpose and layer

Layer: **UTIL**. The single in-tree owner of the BLAKE3 hash function
(plain and derive-key domain modes, fixed 32-byte output) and the `ContentHash` identity
type. Extracted verbatim from `fs-ledger` (bead 7uq9) so exactly one
BLAKE3 implementation exists in the workspace; `fs-ledger` re-exports
these types unchanged. Zero dependencies by design: solver-free
distribution cones (fs-checker's) may depend on it without gaining any
solver, geometry, or license surface.

## Public types and semantics

- `Blake3` — streaming hasher. `new()`, `update(&[u8])` (any split
  pattern; digest equals hashing the concatenation), `finalize() ->
  ContentHash`. `finalize` borrows immutably: it may be called
  repeatedly and interleaved with further `update`s.
- `hash_bytes(&[u8]) -> ContentHash` — one-shot convenience.
- `hash_domain(domain: &str, payload: &[u8]) -> ContentHash` — the
  canonical domain-separation scheme for every typed 32-byte root in the
  workspace. It uses BLAKE3's standard derive-key construction:
  `DERIVE_KEY_CONTEXT` for the domain and `DERIVE_KEY_MATERIAL` for the
  payload. Mode flags separate typed roots from plain `hash_bytes`
  artifacts as well as from other domains.
- `ContentHash(pub [u8; 32])` — Copy identity value. `as_bytes`,
  `to_hex` (64 lowercase hex chars), `from_hex` (either case, exactly
  64 chars, else `None`), `from_slice` (exactly 32 bytes, else
  `None`). `Display` renders `to_hex`; `Debug` wraps it.

## Invariants

- Output is bit-identical to the official BLAKE3 specification for
  plain-mode hashing (spec/oracle vectors in the unit tests and in
  fs-ledger's `ledger_001` conformance suite).
- Streaming and one-shot hashing of the same byte sequence produce the
  same digest for every split pattern.
- `hash_domain(d1, p1) == hash_domain(d2, p2)` implies (up to a BLAKE3
  collision) `d1 == d2 && p1 == p2`; equality between a tagged hash and a
  plain `hash_bytes` artifact likewise requires a cross-mode BLAKE3 collision
  because the compression flags differ.
- `from_hex(to_hex(h)) == Some(h)` for all `h`.

## Error model

No panics on any input; the fallible constructors (`from_hex`,
`from_slice`) return `Option` and refuse malformed input with `None`.
Inputs longer than 2^54 chunks are outside the supported envelope
(vastly beyond any in-tree artifact size).

## Determinism class

Pure function of the input bytes — bit-stable across runs, thread
counts, platforms, and ISAs. No floating point, no time, no I/O.

## Cancellation behavior

Not applicable: all operations are synchronous, allocation-free (apart
from `to_hex`'s String), and short.

## Unsafe boundary

None. 100% safe Rust; no capsule.

## Feature flags

None.

## Conformance tests

Unit tests in `src/lib.rs`: official empty-input spec vector, oracle
`abc` vector, hex round-trip and rejection, streaming-vs-one-shot
across block/chunk edges, single- and multi-chunk derive-key oracle
vectors, and raw/tagged namespace separation. The historical
multi-chunk / multi-level tree vectors continue to run in fs-ledger's
`ledger_001` conformance suite through the re-exported paths.

## No-claim boundaries

General-purpose keyed hashing, a public KDF API, and extended (XOF) output
are NOT implemented. `hash_domain` uses the standard derive-key modes only
for public, non-secret identity namespaces. No constant-time claim is made
(content addressing, not secret handling). No SIMD/multithreaded throughput
claim is made.
- Domain strings are caller-owned protocol identifiers. Callers must use
  hardcoded, globally unique, versioned contexts; this crate does not register
  or deduplicate names dynamically.
