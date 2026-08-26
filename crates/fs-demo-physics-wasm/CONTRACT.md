# CONTRACT: fs-demo-physics-wasm

Status: **v1 — live surface.** Browser boundary for the CMA-ES explainer
site's parametric wing and bridge demo physics. Own nested workspace;
wasm-bindgen confined to wasm32; native builds stay dependency-clean.

## Purpose

Give the explainer site (cmaes_explainer) real FrankenSim-computed physics
for its two design-optimization demos. The formulas are a 1:1 port of the
site's TypeScript analytic models, so WASM-vs-fallback is provenance, not
behavior: same inputs produce the same displayed numbers either way.

No-claims: teaching/viz surrogates. NOT certified structural or aeroelastic
analysis; no time integration, no FEA.

## Public surface (wasm32 exports)

| Export | Signature | Returns |
|---|---|---|
| `wing_eval` | (aspect_ratio, sweep_deg, thickness_ratio, max_camber, camber_position, taper_ratio, family_id: u32, rib_count, cruise_mach) -> String | JSON envelope |
| `bridge_eval` | (span_m, sag_m, deck_stiffness, topology_id: u32, material_id: u32, suspender_count, tower_aspect, damping, truck_pos_m) -> String | JSON envelope |
| `demo_physics_kernel_version` | () -> String | `"fs-demo-physics-wasm 0.1.0"` |

## Envelope contract (frozen v1)

- wing success: `{"ok":{"kernel","liftCoeffCL","dragCoeffCD","inducedDragCDi",
  "profileDragCD0","waveDragCDw","liftToDragRatio","rootBendingMomentKNm",
  "wingMassKg","criticalMach","costScore"}}`
- bridge success: `{"ok":{"kernel","totalMassTons","maxVonMisesStressMPa",
  "maxDeflectionMm","cableTensionKN","flutterCriticalSpeedKmh",
  "yieldLimitMPa","isCompliant","costScore"}}`
- refusal: `{"refusal":{"code","message","ranked_repairs"}}`

Field names match the site's `WingAnalysisResult` / `BridgeAnalysisResult`,
including the site's display rounding, so the adapter is a straight spread.

## Model summary

Wing: 3D lifting-line lift slope with sweep/compressibility, Oswald induced
drag, form-factor profile drag, Korn drag-divergence Mach + Lock fourth-power
wave drag, half-wing root bending moment at the elliptical lift centroid
(cruise atmosphere ~35 kft), spar/skin/rib mass rollup.

Bridge: parabolic-cable tension H = wL²/8s with T = H·√(1+16(s/L)²), beam
deflection 5wL⁴/384EI with cable relief, point-load deck moment Pa(L−a)/L,
suspender-sized cable area (stress is a design output), simplified Selberg
flutter estimate, mass rollup, quadratic constraint penalties.

## Refusal codes

`input-non-finite` · `aspect-ratio-non-positive` · `span-non-positive` ·
`family-id-out-of-range` (0..=4) · `topology-id-out-of-range` (0..=4) ·
`material-id-out-of-range` (0..=3) · `non-finite-result`

Nothing is silently clamped; no traps cross the boundary.

## Determinism

Pure functions of scalar inputs: no wall-clock, no entropy. Same inputs ⇒
byte-identical envelope, native and wasm. Display fields carry the site's
rounding; `costScore` (the CMA-ES objective) is quantized to 1e-6 on both
sides of the boundary so sub-ULP libm differences between Rust and V8 cannot
flip rankings on near-ties.

## CI gates

`cargo check --locked --target wasm32-unknown-unknown` then
`wasm-pack build --target web --release`. Native tests: `cargo test`
(invariant and behavior checks, including sweep-raises-Mdd and
stiffer-deck-deflects-less).
