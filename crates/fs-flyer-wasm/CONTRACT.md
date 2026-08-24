# fs-flyer-wasm — CONTRACT

Layer: L6. Bead: frankensim-wf-root-guzez.1.3 (E0.3, Wright Flyer program).
Spec: COMPREHENSIVE_PLAN_FOR_REAL_TIME_WRIGHT_FLYER_SIM_WITH_FRANKENSIM.md
(ROUND 6 steady state), §7.1.

## Purpose and layer

Layer L6 (boundary crate). What this crate IS:

The single boundary between the Wright Flyer presentation plane (three.js)
and the simulation plane (FrankenSim physics), compiled to
`wasm32-unknown-unknown`. It follows the `fs-wasm` pattern exactly:

- **Own `[workspace]`** with a nested `Cargo.lock`: browser builds never
  depend on the state of unrelated in-progress native crates. The lock is
  reviewed; CI must fail if a build mutates it (lock-drift gate).
- **`cdylib` + `rlib`**: the same pure functions compile natively (tests,
  goldens) and to wasm (the shipped boundary).
- **wasm-bindgen confined to wasm32** (`[target.'cfg(target_arch =
  "wasm32")'.dependencies]`): native builds stay dependency-clean.
- **asupersync canonical browser profile pinned** (`wasm-browser-prod`):
  the fs-exec cone requires exactly one canonical wasm profile feature;
  this workspace pins it so downstream flyer crates inherit a working
  wasm32 build (plan §11.3 standing caveat).

## Invariants

Boundary contracts, inherited by every future API entry:

1. **Typed-refusal JS envelope.** Every fallible entry returns
   `{"ok": ...}` or `{"refusal": {"code", "message", "ranked_repairs"}}`.
   Codes are stable machine-readable strings; repairs are ranked, most
   likely fix first. No silent clamping; no runtime traps.
2. **Determinism digests.** Trajectory-bearing entries expose an
   `fs-blake3` content digest under a versioned domain
   (`org.frankensim.fs-flyer-wasm.hello-trajectory.v1` for the hello
   kernel). Bit-identical runs produce identical digests; the pinned
   canonical golden seeds the six-lane golden program (E6.2). The
   golden-bump protocol applies to any change.
3. **Admission caps are exact.** Refusal boundaries are tested at cap AND
   cap+1 (workspace law).

## Public types and semantics

### v0 surface (E0.3 scaffold)

| Entry | Kind | Contract |
|---|---|---|
| `hello_spin` / `flyer_hello_spin` | pure / wasm | deterministic free rigid-body CG2 spin (`fs_time::lie::rigid_body_step`); refuses non-finite input, non-positive inertia, non-unit quaternion, out-of-domain dt, steps > 1,000,000 |
| `hello_digest` / `flyer_hello_digest` | pure / wasm | full-trajectory bit-exact content digest (hex) under the versioned domain |
| `hello_spin_json`, `refusal_envelope`, `hello_envelope` | pure | the JS envelope renderers (shared by native tests and the boundary) |

### Archive-loader surface (E0.9c slice, bead guzez.1.9.3)

`archive::` implements the verifying-loader mechanics of the frozen E0.9a
contract (`data/wright-flyer/replay-identity-schema-v1.json`) against a
LOCAL content-addressed store (`data/wright-flyer/archive-fixture/`):

| Entry | Contract |
|---|---|
| `verify_target_bytes` | size check BEFORE hashing, then BLAKE3 content identity vs the targets manifest |
| `verify_dual_publication` | both copies verified independently + byte equality (the read-back rule; never provider-native replication) |
| `parse_hello_envelope` | STRICT fail-closed v1 parser — exact nine keys in canonical order; any deviation refuses (a tolerant line reader is how a gate dies silently) |
| `replay_generation` | backward playback: re-executes the archived generation, compares trajectory digests; divergence is a typed refusal, the "old-exact" contract's teeth |

The fixture archives generation 0 (the canonical hello scenario, dt as the
integer ratio 1/120) with its pinned trajectory digest; the battery drives
corruption, truncation, mirror-divergence, malformed-envelope, and
wrong-kernel-replay twins.

### Engine surface v1 (E5.1, bead guzez.6.2)

`engine::EngineSlot` wraps the REAL `fs_flyer::simloop::SimLoop` (the
E4-certified physics: equilibrated tick-0 closure, coupled prop–airframe
force build-up per tick, canard mechanism, perception/pilot family, OU
gust). Scalars in, envelopes out — no JSON parsing at the boundary, no
serde, no mocks. The wasm exports are 1:1 wrappers over a thread-local
slot; the SAME functions are battle-tested natively (`engine_battery`)
and from node against the actual wasm binary
(`node-harness/engine_harness.mjs`).

| Entry | Contract |
|---|---|
| `flyer_engine_init(seed, rho, headwind, mode, member, rail_m, max_ticks)` | admission (caps at cap AND cap+1) → equilibrate → RunIntentId minted AFTER the tick-0 digest; returns run identity + trim; replaces any prior run (the E5.0 ring epoch bump is the consumer-side guard) |
| `flyer_engine_step(has_input, lever_n, warp_rad)` | one 120 Hz step; mode words: 0=fixed, 1=historical(member), 2=human (input REQUIRED every tick — absent or non-finite input is a typed refusal, never a silent zero-hold); the success envelope carries reduced lateral `p/phi/r/psi` as real engine state |
| `flyer_engine_digest()` | chained per-tick blake3 digest (`org.frankensim.wf.sim-digest.v1`) — bit-identical lifecycles are checkable |

Terminal events arrive IN-BAND as `phase` words (`ended:ground-contact`,
`ended:rail-end-without-lift`, `ended:max-ticks`,
`ended:envelope-exceeded`); a mid-flight aero refusal ENDS the run as
`ended:envelope-exceeded` with `envelope_refusal_code` carried in the
same envelope — the UI gets a receipted flight ending, never a trap.

## No-claim boundaries

- The hello kernel is a torque-free rigid body. It makes **no aerodynamic,
  historical, or Flyer-specific claim**; it exists to prove the toolchain,
  the refusal envelope, and the determinism discipline end to end.
- Native-vs-wasm digest identity is a **goal enforced by E6.2's six-lane
  golden program**, not a claim this crate can certify alone; until those
  lanes run in CI, cross-platform identity is `Estimated`.
- Nothing here allocates authority over the FlyerScenario/ModelId identity
  machinery (E0.9); the envelope carries no identity claims yet.
- The lateral envelope fields inherit `fs-flyer`'s Estimated reduced-model
  ceiling. Transporting them to the scene does not promote them to calibrated
  historical or structural authority.

## Error model

Every fallible entry returns the typed-refusal JS envelope
(`{"refusal": {"code", "message", "ranked_repairs"}}`); codes are
stable machine-readable strings tested at cap AND cap+1. Vocabulary:

`non-finite-input`, `non-positive-inertia`, `non-unit-quaternion`,
`timestep-outside-domain`, `step-budget-exceeded`; archive loader:
`archive-size-mismatch`, `archive-content-digest-mismatch`,
`archive-mirror-divergence`, `archive-envelope-malformed`,
`archive-replay-digest-mismatch`; ring (E5.0): `ring-config-invalid`,
`ring-abi-mismatch`, `ring-publish-invalid`, `ring-lease-torn`,
`ring-epoch-stale`; engine (E5.1): `engine-not-initialized`,
`mode-invalid`, `scenario-invalid`, `control-input-missing`,
`run-ended` (+ native physics refusals pass through verbatim at init;
mid-flight aero refusals become the `ended:envelope-exceeded`
terminal).

## Determinism class

Bit-exact within a lane: identical inputs produce identical trajectories
and identical `fs-blake3` digests under versioned domains (per-lane
SELFTEST goldens: native aarch64 debug+release, wasm-node dev+release).
Cross-lane identity is enforced by the E6.2 six-lane golden program;
the tracked cross-lane divergence class is bead guzez.7.2.1. All
transcendentals route through `fs-math det::`; randomness is philox
counter-addressed (no wall-clock, no thread order).

## Cancellation behavior

The boundary is synchronous and single-shot per call: no entry blocks,
polls, or spawns. Cancellation is the CALLER's concern (the worker drops
the slot; a new `flyer_engine_init` replaces any prior run, with the
E5.0 ring epoch bump as the consumer-side stale-guard). Long-running
loops live sim-side under `fs-exec` checkpoints, not here.

## Unsafe boundary

No `unsafe` in this crate. The wasm ABI surface is wasm-bindgen-generated;
shared-memory transport (SAB seqlock) lives in the JS worker against the
protocol invariants tested by `ring_battery` and the node harness.

## Feature flags

None. The wasm32 target gates wasm-bindgen via
`[target.'cfg(target_arch = "wasm32")'.dependencies]`; native builds are
dependency-clean. No cargo features alter semantics.

## Conformance tests

`tests/engine_battery.rs` (lifecycle vs native goldens),
`tests/ring_battery.rs` (E5.0 seqlock/lease/epoch laws),
`tests/archive_fixture.rs` (verifying-loader twins),
`tests/fieldlease_battery.rs` (E7.1-ii lease/staleness/claims), and
`node-harness/engine_harness.mjs` (the SAME exports exercised against
the actual wasm binary). Repro: `cd crates/fs-flyer-wasm && cargo test`.
