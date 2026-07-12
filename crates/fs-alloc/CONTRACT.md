# CONTRACT: fs-alloc

## Purpose and layer
Scope arenas with cancellation-time reclaim independent of allocation count,
unconditional 128-byte alignment, hugepage-eligible chunks with the decision
recorded, sharded object pools, and diffable allocation-site accounting
(plan §5.3). Layer: L0. Depends on `fs-obs` (UTIL) only.

## Public types and semantics
- `ALLOC_ALIGN: usize = 128` — the unconditional alignment policy (superset
  of Apple aarch64's 128 B and x86-64's 64 B cache lines). Every arena
  allocation is aligned to at least this.
- `CachePadded<T>` — `#[repr(align(128))]` wrapper sized/aligned so adjacent
  instances never share a cache line.
- `ArenaPool` (`Clone + Send + Sync`) — configuration (`ArenaConfig`),
  chunk free list, budget enforcement, accounting (`PoolStats`), hugepage
  decision (`HugepageDecision`), and site aggregation (`SiteReport`).
  `reservation_bytes_for_slice` computes the fresh-arena first-chunk
  requirement without allocating, using normalized chunk and alignment rules.
- `Arena` (`Send`, deliberately `!Sync`) — one bump arena per unit of scoped
  work. `alloc`, `alloc_slice_fill`, `alloc_slice_with` return references
  borrowing the arena; `scope`/`ArenaPool::scope` run a closure against a
  fresh (child) arena and reclaim it on exit. Types needing `Drop` are
  rejected at compile time (destructors never run).
- `Site`, `SiteStats`, `SiteReport` — static allocation-site tags with
  cumulative padded-byte/count accounting; reports are sorted and
  deterministic (diffable between runs).
- `HugepagePolicy` / `HugepageOutcome` / `HugepageDecision` — intent, probe
  outcome, and recorded reason; `HUGEPAGE_BYTES = 2 MiB`.
- `ShardedPool<T>` / `PoolItem<T>` / `ShardedPoolStats` — per-shard locked
  free lists (shards `CachePadded`), first-touch construction on the
  acquiring thread, RAII return, explicit `detach` for artifacts that
  graduate out of pooled life.
- `PoolStats::to_event_kind`, `SiteReport::to_event_kind`,
  `ShardedPoolStats::to_event_kind`, `HugepageDecision::to_json` — canonical
  JSON payloads riding `fs_obs::EventKind::Custom` toward the ledger.

## Invariants
1. Every allocation is aligned to `max(align_of::<T>(), 128)` and its byte
   range is disjoint from every other allocation's (G0 law, conformance
   alloc-001/alloc-005).
2. Dropping an arena (completion or cancellation) reclaims all its memory
   with cost proportional to its CHUNK count — O(log allocated-bytes) from
   geometric growth — never to its allocation count; chunks recycle through
   the pool free list (alloc-002/alloc-009).
3. Accounting is exact under the stated granularity: per-arena and per-site
   `bytes` equal the sum of 128-padded payload sizes while representable;
   cumulative byte/allocation counters saturate at `u64::MAX` rather than
   wrapping. Dynamic site names are JSON-escaped before entering canonical
   reports/events. `PoolStats::quiescent`
   (no live arenas, all reserved bytes parked in the free list) holds after
   any interleaving of scope completions and cancellations — verified by a
   10^6-cycle G4 storm (alloc-004) and an 8-thread hammer (alloc-006).
4. Cross-scope escape of an allocation, sharing an `Arena` across threads,
   and arena-placing a `Drop` type are COMPILE ERRORS (doctest
   `compile_fail` battery on `ArenaPool::scope`, `Arena::alloc`, `Arena`).
5. The pool budget (`limit_bytes`) bounds OS-reserved bytes (in-use +
   free-listed); on pressure the free list is drained back to the OS before
   refusing. New-chunk bytes are claimed atomically before allocation, so
   concurrent arenas cannot cross the limit through a check-then-increment
   race; counters are released only after the corresponding chunk is
   deallocated, so claimants can observe conservative over-accounting but
   never stale free capacity. Refusal is structured and leaves the pool fully
   usable (alloc-003). After a pressure drain's final failed capacity check,
   acquisition rechecks the free list exactly once so a concurrently published
   suitable chunk is not missed; the retry is bounded and re-enters both gates.
6. `OperationMemoryLease` (bead wf9.16) is the RUN-SCOPED ledger, distinct
   from the pool budget's process-wide one: cloneable, atomic
   reserve/release, canonical v2 receipts (`requested`/`refusals` exact for the
   observed admission trace and saturated at `u64::MAX`; `used` exact;
   `peak_bytes` is the logical high-water; an impossible unmatched release is
   fail-closed and increments `release_invariant_violations` instead of
   disappearing in release builds). Admission precedes recycled-chunk
   removal and fresh allocation, and a failed pool/OS gate rolls the operation
   charge back, so there is no unaccounted hold interval. Cumulative counters
   saturate instead of wrapping, raw
   reserve/release is crate-private, and dynamic refusal labels are JSON
   escaped. A leased arena (`arena_leased`/`scope_leased`) charges every
   chunk — fresh or recycled — for exactly its hold interval and releases on
   drop (unwinds included) before publishing it to the free list; free-list
   inventory belongs to no operation, so recycled chunks are never
   double-charged. Both gates must admit a chunk;
   the refusal names whichever refused (`Exhausted` = pool,
   `LeaseExhausted` = operation). Child scopes inherit the parent's lease.
   Cached acquisition selects the smallest sufficient chunk; when an
   oversized cached chunk cannot fit a bounded operation lease but the normal
   fresh chunk can, the cache remains untouched and the fresh path is tried.
   Identical plans with identical cache state have deterministic cumulative
   demand. Cache history and near-limit concurrent refusal counts,
   first-refusal snapshots, and peaks are schedule/state-sensitive and carry
   no cross-schedule bit-stability claim.

7. `LeasedVec` (bead wf9.16.1) is lease-admitted OWNED storage: its backing
   buffer is admitted through the operation lease BEFORE the allocator is
   asked (admission precedes allocation — no over-limit live interval),
   the RAII charge lives inside the value and releases on drop/unwind,
   and growth admits the full replacement buffer before reallocating (the
   old+new concurrent peak is the honest cost). The detached-empty state
   backs contextless fold identities: charges nothing, refuses growth.
   Payload bytes are tracked exactly; allocator rounding beyond the
   `try_reserve_exact`-sized buffer is a no-claim.

## Error model
All fallible APIs return `Result<_, AllocError>`; `AllocError` is a
structured enum (`Exhausted`, `OutOfMemory`, `LeaseExhausted`,
`LayoutOverflow`,
`ReservationOverflow`) carrying the allocation site, sizes, and budget
context, with teaching `Display` text including ranked fixes (Decalogue P10).
Out-of-memory is a refusal, never an abort (`handle_alloc_error` is not
reachable from this crate). Panics do not
cross the API boundary except caller-supplied closures' own panics
(`alloc_slice_with`), which unwind cleanly: the arena stays usable and the
reserved bytes reclaim with the scope. Lock poisoning from such a panic is
confined; internal locks are held only around plain data structure edits. The
free-list and site-table mutexes are poison-recovering on admission, pressure
drain, recycle, aggregation, and report paths so cleanup during a later unwind
cannot become a second panic; the regression poisons both deliberately and
proves fresh acquisition, recycle, site aggregation, lease release, and
quiescence still work.

## Determinism class
Deterministic: single-threaded workloads with a fixed seed produce
bit-identical `ArenaStats`/`PoolStats`/`SiteReport`/`ShardedPoolStats` JSON
(G5-style check, alloc-007). Reports contain no addresses and no clocks.
Concurrent workloads are deterministic in their FINAL accounting for a fixed
per-thread schedule (per-item counters are exact), but interleaving-derived
quantities (`chunks_created` vs `chunks_recycled` split under contention)
are schedule-dependent — callers wanting bit-stable reports across thread
counts must aggregate per-worker pools keyed by logical identity (fs-exec's
job, plan §5.4).

## Cancellation behavior
No poll points inside this crate (no operation blocks unboundedly; chunk
allocation is one global-allocator call). The crate's cancellation role is
RECLAIM: fs-exec binds one `Arena` per asupersync scope, and cancellation
drops the arena — invariant 2 makes that O(chunks), leak-free (G4 storm),
and independent of how much the cancelled work had allocated. Losing
branches of speculative races reclaim through exactly this path.

## Unsafe boundary
Exactly one registered capsule: `src/raw/mod.rs` (bump-pointer core, chunk
alloc/dealloc; the "arena allocators" zone sanctioned by Decalogue P1),
under 300 lines behind a safe facade, with `src/raw/SAFETY.md` documenting
invariants, the lifetime-erasure argument, and Miri coverage (the capsule is
plain Rust — Miri interprets every path). Registered in
unsafe-capsules.json; enforced by `cargo run -p xtask -- check-unsafe`.

## Feature flags
None. Everything here is `[S]` solid-tier.

## Conformance tests
`tests/conformance.rs`, cases alloc-001..alloc-009 (JSON-line verdicts;
seeded cases carry their seed): unconditional alignment, scope reclaim,
structured budget refusal, the 10^6-cancellation G4 storm (emits a
`storm_assertion` event through fs-obs), shadow-model accounting +
disjointness, concurrent leak-freedom, G5 deterministic reports, recorded
hugepage decisions, and chunk-recycling bounds. Any reimplementation must
pass this suite.
In-module tests additionally verify first-chunk preflight sizing, structured
reservation overflow, and concurrent hard-limit claims.

## No-claim boundaries
- NO claim that hugepages actually back any allocation: without `madvise`
  (P1 forbids FFI; std exposes none) the crate only makes chunks
  THP-*eligible* (2 MiB size + alignment) when Linux THP mode is `always`,
  and RECORDS the decision. Apple platforms: 16 KiB base pages, no THP
  control — recorded as `unsupported_platform`.
- NO NUMA placement claim: first-touch is a hook (pool `make` runs on the
  acquiring worker; arenas are per-worker by design), but node binding and
  CCD-aware placement are fs-exec/fs-substrate territory.
- NO cancel-latency claim (the ≤200 µs tile-boundary target is fs-exec's
  contract); this crate only guarantees reclaim cost and leak-freedom.
- NO claim about accounting granularity finer than 128-byte padding;
  alignments above 128 may consume additional window padding visible only
  as reserved-vs-allocated slack.
- NO claim that an operation receipt is independent of the pool's prior cache
  state: an oversized reusable chunk is charged at its actual size when
  admitted, while a bounded lease may select the smaller fresh path.
- Not yet verified on x86-64/Linux CI (runner pending — same boundary as
  fs-substrate); Linux THP probing is written but exercised only where the
  sysfs file exists.
