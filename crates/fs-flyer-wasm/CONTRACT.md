# fs-flyer-wasm — CONTRACT

Layer: L6. Bead: frankensim-wf-root-guzez.1.3 (E0.3, Wright Flyer program).
Spec: COMPREHENSIVE_PLAN_FOR_REAL_TIME_WRIGHT_FLYER_SIM_WITH_FRANKENSIM.md
(ROUND 6 steady state), §7.1.

## What this crate IS

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

## Boundary contracts (inherited by every future API entry)

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

## v0 surface (E0.3 scaffold)

| Entry | Kind | Contract |
|---|---|---|
| `hello_spin` / `flyer_hello_spin` | pure / wasm | deterministic free rigid-body CG2 spin (`fs_time::lie::rigid_body_step`); refuses non-finite input, non-positive inertia, non-unit quaternion, out-of-domain dt, steps > 1,000,000 |
| `hello_digest` / `flyer_hello_digest` | pure / wasm | full-trajectory bit-exact content digest (hex) under the versioned domain |
| `hello_spin_json`, `refusal_envelope`, `hello_envelope` | pure | the JS envelope renderers (shared by native tests and the boundary) |

## No-claims

- The hello kernel is a torque-free rigid body. It makes **no aerodynamic,
  historical, or Flyer-specific claim**; it exists to prove the toolchain,
  the refusal envelope, and the determinism discipline end to end.
- Native-vs-wasm digest identity is a **goal enforced by E6.2's six-lane
  golden program**, not a claim this crate can certify alone; until those
  lanes run in CI, cross-platform identity is `Estimated`.
- Nothing here allocates authority over the FlyerScenario/ModelId identity
  machinery (E0.9); the envelope carries no identity claims yet.

## Refusal vocabulary (v0)

`non-finite-input`, `non-positive-inertia`, `non-unit-quaternion`,
`timestep-outside-domain`, `step-budget-exceeded`.
