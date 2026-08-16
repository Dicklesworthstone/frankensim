# COMPREHENSIVE PLAN: Real-Time Wright Flyer Simulation with FrankenSim

**Working title:** *First Flight — Kitty Hawk, December 17, 1903*
**Document status:** Planning-workflow ROUND 0 (initial comprehensive plan; not yet reviewed).
**Process:** This document follows `/planning-workflow`. It must survive ≥4 review rounds
(GPT Pro Extended Reasoning per the skill's EXACT PROMPT) and reach steady-state before
conversion to beads. A "Review round log" section at the end tracks rounds.
**Repo:** `frankensim` (156-crate Rust workspace, layers L0–L6, Franken-only dependencies,
determinism doctrine, typed refusals, evidence colors Verified/Validated/Estimated).

---

## Table of Contents

1. [Vision](#1-vision)
2. [Product Experience](#2-product-experience)
3. [Historical & Physical Grounding](#3-historical--physical-grounding)
4. [System Architecture](#4-system-architecture)
5. [Physics Core Design](#5-physics-core-design)
6. [Numerics, Determinism, and the Execution Contract](#6-numerics-determinism-and-the-execution-contract)
7. [WASM Engineering & the Real-Time Budget](#7-wasm-engineering--the-real-time-budget)
8. [Rendering & UX (three.js)](#8-rendering--ux-threejs)
9. [Configuration Space, Experiments & KPIs](#9-configuration-space-experiments--kpis)
10. [Validation & Evidence Program](#10-validation--evidence-program)
11. [Crate Reuse Matrix & New Crates](#11-crate-reuse-matrix--new-crates)
12. [Milestones & Dependency-Aware Task Graph](#12-milestones--dependency-aware-task-graph)
13. [Risks & Mitigations](#13-risks--mitigations)
14. [Open Questions for Review Round 1](#14-open-questions-for-review-round-1)
15. [Appendices](#15-appendices)

---

## 1. Vision

### 1.1 What we are building

A browser-native, real-time, physically honest simulation of the 1903 Wright Flyer's
first flights — the airplane, the air, the ground, and the pilot — rendered beautifully
with three.js, computed by FrankenSim physics compiled to WebAssembly, with the invisible
made visible: at any moment the user can flip on the **wind vector field** and watch the
turbulent boundary layer roll over the Kill Devil Hills sand, see the wingtip vortex
sheets peel off and sink toward the ground, and overlay **divergence, gradient, and curl**
of the velocity field in 3-D as living mathematical objects rather than textbook symbols.

Every load-bearing number of the real Flyer is dialed in — span, chord, camber, canard
area, gross weight, engine power, propeller geometry, the December headwind — and every
one of them is a **user-adjustable configuration** with immediate physical consequence.
Stretch the aspect ratio and watch induced drag fall and the wing spar's load rise. Trim
the canard smaller and watch the pitch instability that nearly killed the project become
untamable. Take away the headwind and discover, as the Wrights did, that the launch rail
is suddenly too short. Then read the numbers: lift margin, distance flown, peak speed,
smoothness of flight, control effort — the same quantities of interest the Wrights argued
about in letters, now live on screen.

### 1.2 Why this project, and why FrankenSim

The Euler-disc E2E pipeline proved the pattern: one unified simulation system produced
mechanics, spectral-rendered video, and physically derived audio, all from a single
identity-tracked configuration with honest evidence labels. This project scales the same
pattern up from a spinning disc on a table to an aircraft in a turbulent atmosphere —
and moves the interactive tier into the browser.

This is not a from-scratch flight sim. The plan's central engineering claim, grounded in
the current tree (verified 2026-08-16):

- **`fs-wasm` already exists and already ships a large slice of FrankenSim to wasm32**
  via its own decoupled workspace and `wasm-pack` browser builds (`crates/fs-wasm`,
  CI lane `fs-wasm-build` in `scripts/ci/quality_lanes.sh`). Its dependency list already
  includes `fs-bem`, `fs-vpm`, `fs-lbm`, `fs-exec`, `fs-alloc`, `fs-viz`, `fs-render`,
  `fs-scenario`, `fs-uq`, `fs-qty`, `fs-material` — i.e., the browser toolchain problem
  is *solved infrastructure*, not a research risk.
- **`fs-ornith` (the ornithoid flagship) already flies a parameterized aircraft in the
  browser** at screening fidelity: `fs-bem` Laplace panel aerodynamics (FMM-accelerated
  by `fs-fmm`) plus an `fs-vpm` vortex-particle wake metric, wrapped in a staged
  campaign. The Wright Flyer sim is a sibling flagship with a real-time loop instead of
  a campaign loop.
- **`fs-mbd`** provides deterministic quaternion rigid-body dynamics; **`fs-time`**
  provides structure-preserving (symplectic) integrators; **`fs-contact`** provides
  capability-routed contact for the skid/rail/sand interactions; **`fs-lbm`** provides
  D2Q9/D3Q19 lattice-Boltzmann cores for the offline high-fidelity validation tier;
  **`fs-vpm`** provides the desingularized Biot–Savart vortex core; **`fs-viz`**
  provides scientific-visualization primitives with analytic ground truth.
- The **evidence doctrine** (Verified/Validated/Estimated, typed refusals, no-claims
  blocks, content-identity binding) is exactly the honesty instrument a "real physics,
  real time" product needs so the real-time tier never silently over-claims.

### 1.3 The one-sentence positioning

> A museum-quality interactive: the most historically and physically serious Wright Flyer
> simulation that runs in a browser tab — where the physics is a first-class citizen you
> can *see*, *interrogate*, and *change*, and where every number tells you how sure it is.

### 1.4 Explicit non-goals (v1)

- **Not a general flight simulator.** One aircraft family (Flyer I/II/III lineage), two
  sites (Kill Devil Hills 1903, Huffman Prairie 1904–05), minutes-long flights.
- **Not a CFD product.** The browser tier is a reduced-order model with a wake; the
  full-field CFD (`fs-lbm`) runs offline in the validation program, not per-frame.
- **No multiplayer, no VR in v1** (VR is a natural v2; the architecture must not preclude it).
- **No gamification** (score chasing, unlocks). The "game" is the physics itself.
- **The browser tier never claims better than `Estimated` evidence color.** Validation
  receipts from the native tier can promote *specific pinned configurations* to
  `Validated`; the UI displays the color honestly.

---

## 2. Product Experience

### 2.1 The five-minute journey (first-time user)

1. **Arrival.** The page loads to a dawn scene: the camp at Kill Devil Hills, the wooden
   launch rail pointing into a 24 mph wind, the Flyer on its dolly, sand streaming in low
   ribbons. Ambient wind audio. A single "Fly" button and a date stamp: December 17, 1903,
   10:35 a.m.
2. **First flight, authentic mode.** The user launches. The engine clatters to speed, the
   dolly accelerates down the 60-foot rail, and the Flyer staggers into the air exactly as
   the historical Flyer did: porpoising in pitch (it was statically unstable), 120 feet in
   12 seconds, and a firm arrival in the sand. A results card appears with the historical
   comparison: *"Orville: 120 ft, 12 s. You: 97 ft, 10.4 s."*
3. **Seeing the air.** The user toggles "Show the wind." The world dims slightly; the
   velocity field appears as animated glyphs and streamribbons: the logarithmic boundary
   layer profile, gust structures advecting through, the bound circulation around the
   wings, tip vortices curling off and descending in ground effect. Toggles for
   |curl| (vorticity magnitude), divergence (near-zero — with an explanation of why),
   and pressure-gradient overlays.
4. **Changing history.** The user opens the Design panel. Sliders with the real 1903
   values marked: wingspan, chord, camber, canard area and arm, gross weight, engine
   power, propeller diameter/pitch, headwind speed. They drag aspect ratio from 6.4
   to 9, watch the predicted lift/drag polar redraw live, fly again, and go 240 feet.
   A delta card explains *why* (induced drag ∝ 1/AR).
5. **The hook.** A "Challenges" rail: *Beat Wilbur's 852 feet. Survive a 30 mph gust
   day. Fly Huffman Prairie with no headwind (you'll want the catapult). Tame the pitch
   oscillation by moving the CG.* Each challenge is a scenario config with its own
   leaderboard-free results card and replay.

### 2.2 The advanced journeys

- **The engineer's loop.** Full config editor (typed, unit-checked via `fs-qty`
  semantics), KPI dashboard (lift margin, L/D, static margin, phugoid/short-period
  frequencies from live linearization, control-effort RMS, smoothness index), parameter
  sweeps run in Web Workers with result overlays (e.g., distance-vs-camber curve), and
  export of any run as a replay file with full config identity.
- **The educator's loop.** Curated "lessons": ground effect (fly at 3 m vs 30 m and
  compare induced drag live), why the Flyer's anhedral helped in gusts, why the canard
  configuration is statically unstable but was chosen anyway (stall recovery behavior),
  div/grad/curl as physical objects with the wind field as the worked example.
- **The historian's loop.** Authentic-day presets: the four December 17 flights with
  recorded wind conditions; the 1904 Huffman Prairie problem (light winds → catapult);
  Flyer III 1905 (first practical airplane; 24.5 miles on Oct 5, 1905). Each preset
  cites its sources in-app.
- **The cinematic loop (offline, native).** Any browser replay can be exported and
  re-rendered by the native FrankenSim cinematic pipeline (`fs-render` spectral tracer +
  the existing EXR→ProRes quarantined mux adapter from the Euler-disc program) into a
  film-quality clip — same trajectory identity, same config identity, hero-shot camera
  paths. The browser is the instrument; the native pipeline is the film studio.

### 2.3 Interaction model (flying the Flyer)

The real Flyer's controls were bizarre by modern standards and that is part of the story:

- **Pitch:** front elevator (canard) via lever — sensitive, statically unstable aircraft.
- **Roll+yaw:** wing-warping *coupled to rudder* via the hip cradle — the pilot shifted
  his hips to bank.
- **Throttle:** none in the modern sense; engine ran essentially full-out.

Input mappings (all remappable):

| Mode | Pitch | Roll/yaw (coupled warp+rudder) | Notes |
|---|---|---|---|
| Keyboard | ↑/↓ | ←/→ | default; discrete-feel filtered to lever rate |
| Mouse | vertical drag | horizontal drag | "hip cradle" metaphor |
| Gamepad | left stick Y | left stick X | best feel for instability |

**Stability-assist ladder** (crucial for playability without lying about physics):
- **Authentic** — raw aircraft, raw instability; the historical difficulty.
- **Wright-hands** — a modeled pilot inner loop (the trained reflexes Wilbur/Orville
  earned with ~1,000 glider flights in 1900–1902): a low-gain pitch-rate damper with
  human-plausible latency (~180 ms) and saturation. This is a *modeled historical pilot*,
  not a cheat, and is labeled as such.
- **Modern SAS** — full stability augmentation for accessibility; labeled "not historical."

All three run the same airframe physics; assists only add control inputs through the same
actuator path with the same rate limits, so KPIs remain comparable and honest.

### 2.4 Visualization modes (the differentiator)

1. **Vector glyph field** — GPU-instanced arrows on a user-positioned 3-D probe box or
   plane; length/color = |u|; live animated.
2. **Streamribbons / streaklines** — integrated through the *current* unsteady field;
   emitted from the leading edge, wingtips, or a hand-placed rake.
3. **Vorticity view** — |ω| volume glow + extracted vortex-core polylines for the tip
   vortices and shed wake (the wake data structure *is* vorticity, so this is exact,
   not post-processed).
4. **Divergence view** — near-zero by construction (incompressible model, div-free gust
   synthesis); displayed as a *verification overlay* with a callout: "the model promises
   ∇·u = 0; the residual you see is sampling/discretization error" — a teachable moment
   that doubles as an on-screen self-check.
5. **Gradient/pressure view** — ∇p from the reduced model along surfaces (wing pressure
   bands via section theory) and ∇|u| in the volume probe.
6. **Force overlay** — per-strip lift/drag vectors on the wings, thrust vectors on the
   props, weight/CG marker, net-force and net-moment gnomon at the CG.
7. **Scalar probes** — draggable point probes reading (u, p̂, ω, T?) with strip-chart
   history; a virtual anemometer you can plant on the dune like the Weather Bureau's.

Every mode works while flying, paused, or scrubbing a replay.

---

## 3. Historical & Physical Grounding

> Grounding rule (planning-workflow): numbers below marked **[V]** are load-bearing and
> must be re-verified against primary/secondary sources during E1 (task E1.2) before any
> code pins them. Values marked **[V?]** are believed-approximate and must be sourced or
> demoted to tunable-with-default. Do not let unverified numbers survive into goldens.

### 3.1 The 1903 Wright Flyer — reference configuration

| Property | Value | Status |
|---|---|---|
| Wingspan | 40 ft 4 in = 12.29 m | [V] |
| Chord | 6 ft 6 in = 1.98 m | [V] |
| Wing area (both wings) | ~510 ft² = 47.4 m² | [V] |
| Aspect ratio (per wing) | ~6.2–6.4 | [V] derived |
| Gap between wings | ~6 ft (≈ chord) | [V?] |
| Camber | ~1/20 (flown at Kitty Hawk 1903) | [V] |
| Airfoil | thin, single-surface-ish, ribs + fabric; Wright #12-derived section | [V?] |
| Anhedral (wing droop) | slight, deliberate (gust response at Kitty Hawk) | [V] qualitative |
| Canard ("front rudder") area | ~48 ft², biplane elevator | [V?] |
| Canard arm (CG→canard AC) | ~7 ft ahead | [V?] measure from drawings |
| Twin rudder area (aft) | ~21 ft² | [V?] |
| Empty weight | ~605 lb = 274 kg | [V] |
| Gross weight (with pilot ~145 lb) | ~750 lb = 340 kg | [V] |
| Engine | 4-cyl inline, ~201 in³, ~12 hp peak (~8–9 hp sustained), ~180 lb w/ accessories | [V] |
| Transmission | chain drive, 23:8 reduction; one chain crossed for counter-rotation | [V?] ratio |
| Propellers | 2 × pusher, ~8.5 ft (2.6 m) diameter, ~330–350 rpm, carved spruce | [V] |
| Propeller efficiency | modern reconstructions estimate ~70–80% peak | [V?] cite AIAA/Wright Experience data |
| Cruise/flight airspeed | ~30–34 mph (13–15 m/s) | [V] |
| Stall speed | within a few mph of flight speed (flying near CL,max) | [V?] derive + cite |
| Static pitch stability | **unstable** (canard + aft-ish CG); the famous porpoising | [V] (Culick et al.) |
| Pilot position | prone, head forward, hip cradle | [V] |

### 3.2 The four flights of December 17, 1903 (validation anchors)

| # | Pilot | Distance | Duration | Notes |
|---|---|---|---|---|
| 1 | Orville | 120 ft | 12 s | 10:35 a.m.; undulating pitch; wind ~24–27 mph [V] |
| 2 | Wilbur | ~175 ft | ~12 s | [V?] |
| 3 | Orville | ~200 ft | ~15 s | [V?] |
| 4 | Wilbur | 852 ft | 59 s | pitch oscillation → hard landing, canard broken [V] |

Ground speed on flight 1 was ≈ 6.8 mph against the headwind → airspeed ≈ 31–34 mph.
These four flights, with the recorded wind, are the **primary end-to-end validation
anchors**: the reference config in authentic wind must reproduce distance/duration within
declared tolerance bands (see §10.3).

### 3.3 Sites

**Kill Devil Hills, NC (1903).** Flat sand plain; launch rail ("the Junction Railroad"),
60 ft of 2×4 on edge, iron-capped; dolly/truck riding the rail; the Flyer on skids (no
wheels). Big Kill Devil Hill (~90 ft dune) nearby but the powered flights launched from
the flat. Steady strong Atlantic winds; Dec 17 measured 24–27 mph [V]. Terrain: sand
(aerodynamic roughness z₀ ≈ 10⁻³–10⁻² m), sparse grass, camp buildings (two sheds).

**Huffman Prairie, Dayton OH (1904–05).** ~84-acre cow pasture; light and variable winds
→ from Sept 1904 the **catapult derrick**: ~1,600 lb weights [V?] dropped from ~16 ft,
rope-and-pulley to the dolly — an assisted-takeoff system the sim must include as a
launch option. Rougher air (thermals, tree lines), longer flights, first circles
(Sept 20, 1904, Flyer II) and the 24.5-mile Oct 5, 1905 flight (Flyer III).

Terrain data plan: USGS 3DEP lidar / SRTM for modern topography of both sites, plus
historical shoreline/dune position references for 1903 (the dunes have moved; the
memorial's landscaping postdates the flights). The terrain need not be survey-grade:
a heightfield at ~1–2 m horizontal resolution over ~2×2 km, hand-adjusted to the
historical photographs, with a documented "representative, not surveyed" label.

### 3.4 The aerodynamic character we must capture (fidelity floor)

These are the phenomena that make the Flyer *the Flyer*; the reduced model is designed
around reproducing them (each maps to a validation case in §10):

1. **Pitch instability + pilot-in-the-loop oscillation** (the porpoising). Requires:
   unstable short-period dynamics from the canard configuration, canard lift-curve and
   downwash-free lead, correct-order pitch inertia, elevator effectiveness, and pilot
   latency. This is THE signature behavior.
2. **Ground effect.** Flights at 1–4 m altitude with a 12.3 m span: strong reduction of
   induced drag and downwash. Method: image-vortex system (exact for the reduced model).
3. **Near-stall cruise.** The Flyer flew at high CL close to stall; lift-curve flattening
   and drag rise near CL,max must be modeled (blended stall model per strip), else
   distance/speed KPIs will be fantasy.
4. **Biplane interference.** Two wings a chord apart: mutual induced effects reduce
   effective aspect ratio (Munk factor). Strip model must include gap-dependent
   interference or lift will be ~10–20% optimistic.
5. **Low-airspeed propulsion.** Props at high advance ratio variation during the rail
   run; thrust falls with airspeed; BEMT (blade-element/momentum) with the actual carved
   geometry class, chain-drive RPM ratio, and engine torque curve.
6. **Gusty boundary-layer wind.** Log-law mean profile + anisotropic turbulence with
   realistic length scales; the aircraft's gust response (Küssner-type lag) drives the
   smoothness KPI and the historical difficulty.
7. **Warp-drag / adverse yaw and the coupled warp–rudder control.** The Wrights coupled
   them for a reason; decoupling them (a config toggle!) should reproduce the 1901–02
   adverse-yaw problem.

---

## 4. System Architecture

### 4.1 Three planes

```
┌────────────────────────────────────────────────────────────────────────┐
│ PRESENTATION PLANE (TypeScript + three.js, apps/wright-flyer/)         │
│ scene graph, terrain+sky, aircraft model, field visualization,         │
│ HUD/instruments, config UI, replay scrubber, WebAudio                  │
│         ▲ zero-copy typed-array views + postMessage commands           │
├────────────────────────────────────────────────────────────────────────┤
│ SIMULATION PLANE (Rust→wasm32, crates/fs-flyer-wasm, own workspace     │
│ like fs-wasm)                                                          │
│ fs-flyer (aircraft assembly, 6-DOF, pilot, launch, KPIs)               │
│ fs-wing (lifting-surface aero: strips + unsteady states + VLM tier)    │
│ fs-airscrew (prop BEMT + engine + chain drive)                         │
│ fs-atmo (wind profile + div-free turbulence + gust events)             │
│ on: fs-mbd, fs-time, fs-contact, fs-vpm, fs-bem, fs-exec, fs-math det, │
│     fs-rand philox, fs-blake3 identity, fs-qty units, fs-scenario      │
├────────────────────────────────────────────────────────────────────────┤
│ TRUTH PLANE (native FrankenSim, offline)                               │
│ fs-lbm D3Q19 wind-over-terrain + high-res VPM/BEM aero re-runs;        │
│ cross-fidelity residual receipts (fs-vvreg / vv-scorecard);            │
│ cinematic re-render of replays (fs-render + euler mux adapter);        │
│ golden replay hashes across aarch64/x86/wasm (four-quadrant doctrine)  │
└────────────────────────────────────────────────────────────────────────┘
```

**Why this split.** three.js is the right tool for the *display* problem (mature scene
graph, materials, instancing, ecosystem) and is outside the Rust workspace so Franken-only
dependency doctrine is untouched. The wasm crate is the *only* boundary between planes,
and it follows the proven `fs-wasm` pattern: an isolated `[workspace]` so browser builds
never depend on unrelated in-progress native crates (this decoupling is already documented
in `fs-wasm/Cargo.toml` and battle-tested by the `fs-wasm-build` CI lane). The truth plane
is where "don't dumb down the physics" is *enforced* rather than merely promised: the
same config identities re-run at higher fidelity, and the differences become published
residuals, not vibes.

### 4.2 The multi-fidelity ladder (the honest answer to "real physics, real time")

| Tier | Where | Model | Cost | Role |
|---|---|---|---|---|
| A | browser, always | strip-theory lifting-line w/ unsteady states (Wagner/Küssner), image ground effect, biplane interference, BEMT prop, synthetic div-free turbulence | ~1 ms/step | the flyable sim; KPIs |
| B | browser, capable machines | Tier A + live vortex-lattice bound sheet + shed vortex-particle wake (fs-vpm kernel), wake feedback at reduced rate | ~4–8 ms/step | wake realism, viz truth |
| C | native, offline | high-res unsteady VLM/panel (fs-bem+fs-fmm) + dense VPM wake; long-run KPI recomputation | minutes | validates A/B force models |
| D | native, offline | fs-lbm D3Q19 wind-over-terrain (no aircraft) → gust statistics; optionally immersed aircraft at reduced Re for qualitative field comparison | hours | validates fs-atmo statistics; hero visuals |

Ladder invariants:
- **One geometry, one config identity** flows through all tiers (content-hashed via
  `fs-blake3`, schema'd via `fs-scenario` conventions).
- Tier A/B force coefficients are tabulated/corrected against Tier C runs at pinned
  configs; the correction tables ship with their provenance and residuals.
- The UI evidence badge reflects the ladder: `Estimated` by default; `Validated (vs
  Tier C @ cfg 3f2a…)` when the current config is within the validated envelope
  (defined in §10.4); never `Verified` (no interval-certified claims in v1).

### 4.3 Process/thread topology in the browser

```
Main thread:        three.js render loop @ display Hz; UI; input capture
Worker "sim":       wasm physics @ fixed 120 Hz (Tier A) / +wake 40 Hz (Tier B)
Worker "field":     wasm field-sampling service @ 15–30 Hz (viz grids, probes)
Worker pool (0–N):  parameter sweeps, replay re-runs (background experiments)
SharedArrayBuffer:  state ring buffer (sim→render), field buffers (field→render)
```

- **Fixed-timestep** sim decoupled from render; render interpolates between the last two
  states (standard, but *mandatory* here for determinism and replay identity).
- SharedArrayBuffer requires COOP/COEP headers; the app must also run (degraded to
  `postMessage` transfer, single-threaded wasm) when cross-origin isolation is absent.
  Feature detection with an honest banner ("running in compatibility mode: Tier A").
- wasm SIMD128 used in the Biot–Savart and BEMT inner loops via `fs-simd`'s tier
  discipline (Tier 0 scalar reference always available and always the goldens' referee).

### 4.4 Data flow (one frame)

1. Input events → control lever targets (rate-limited actuator model in `fs-flyer`).
2. Sim worker advances N fixed steps: atmosphere sample → aero strips (+wake) → prop/
   engine → contact (rail/skid/ground) → 6-DOF integrate → KPI accumulators → state
   snapshot to ring buffer (position, quaternion, rates, per-strip forces, wake summary,
   instrument values).
3. Field worker (if a viz mode is on): samples the *same* model (mean wind + turbulence
   + bound/wake induction via shared state) onto the user's probe grid; computes
   div/curl/grad on that grid; writes to field buffer.
4. Main thread: interpolate state; update aircraft pose/control-surface morphs; update
   instanced glyphs/ribbons from field buffer; draw HUD; WebAudio parameter update.

### 4.5 Replay & identity (the Euler-disc lesson, applied)

Every run records: config content-hash, atmosphere seed (philox key), input trace
(timestamped lever positions), and code identity (crate versions + wasm build hash).
A replay file re-runs bit-identically in the browser (same wasm) and *numerically
compatibly* natively (`fs-math det::` routing makes libm differences a non-issue; the
four-quadrant golden program in §10.5 makes "compatibly" a measured claim, not hope).
Replays are the bridge to the truth plane: validation re-runs and cinematic exports
consume replay files, never screen recordings.

---

## 5. Physics Core Design

### 5.1 `fs-flyer` — aircraft assembly & flight dynamics (new crate, L4)

**Owns:** the parametric airframe (geometry+mass), control system model, pilot models,
launch systems, ground interaction, the per-step orchestration, KPI accumulation, and
the scenario schema. Depends on `fs-wing`, `fs-airscrew`, `fs-atmo`, `fs-mbd`,
`fs-time`, `fs-contact`, `fs-qty`, `fs-scenario`, `fs-blake3`, `fs-exec`.

#### 5.1.1 Parametric airframe

The design panel edits a `FlyerDesign` struct; everything downstream derives from it.
High-level levers (user-facing) map to detailed geometry deterministically:

```
FlyerDesign {
  // Lifting system
  span_m,                    // default 12.29
  chord_m,                   // default 1.98  (aspect ratio is derived+displayed)
  gap_over_chord,            // default ~1.0 (biplane gap)
  camber_ratio,              // default 1/20; range 1/25..1/8 (the Wrights' own sweep!)
  stagger_m,                 // default ~0
  anhedral_deg,              // default ~ -1.5 (droop), range -5..+5 (try dihedral!)
  washout_deg,               // twist root→tip
  // Canard & tail
  canard_area_m2, canard_arm_m, canard_span_m,
  rudder_area_m2, rudder_arm_m,
  // Mass & balance
  empty_mass_kg, pilot_mass_kg, ballast_kg, cg_offset_m,   // fore/aft CG shift
  // Propulsion
  engine_power_hp, engine_rpm_max, chain_ratio,
  prop_diameter_m, prop_pitch_m, prop_blade_count, prop_activity_factor,
  // Control system
  warp_rudder_coupled: bool,     // the Wrights' coupling; false = 1901 mode
  elevator_rate_limit, warp_rate_limit,
  // Structure (v1: limits only, no aeroelasticity)
  structural_load_factor_limit,  // display + refusal, default ~2.5 [V?]
}
```

Derived+displayed: wing area, aspect ratio, wing loading, power loading, static margin,
tail volume coefficients, CL required at target speed. **Typed refusals** at admission
for non-physical inputs (negative areas, CG outside airframe, camber beyond section-model
validity) — the same admission-time refusal doctrine the Euler cinematic config uses
(refuse in milliseconds with a ranked-repair message, not after a run).

Mass model: component build-up (wings by area×areal-density, engine, radiators, chains,
pilot at prone station, ballast) → mass, CG, inertia tensor. Inertia derived from the
component distribution, not hand-waved — pitch inertia matters for the instability
signature. Cross-check task: published estimates of Flyer inertias from the AIAA Wright
Flyer Project [V].

#### 5.1.2 Six-DOF dynamics

`fs-mbd` state (position, canonical unit quaternion, linear/angular velocity in its
documented convention) advanced by `fs-time` integrators. Forces/torques assembled each
step from: gravity, per-strip aero (5.2), canard/rudder aero, prop thrust+torque pair
(counter-rotating → net torque ≈ 0, but *config-breakable*: disable one prop and feel
the asymmetry), fuselage/strut/wire parasite drag (flat-plate build-up with the Wrights'
own strut-fairing insight), ground/rail contact.

Integrator: fixed dt = 1/120 s primary; substep to 1/240 s automatically when contact
is active (rail run, landing). Semi-implicit symplectic base from `fs-time`; quaternion
renormalization per `fs-mbd`'s canonical convention. dt-refinement study is a validation
task (V-05), mirroring the Euler-disc "no full-duration dt convergence certificate"
honesty note until it exists.

#### 5.1.3 Control & actuation

Lever dynamics: first-order rate-limited actuators (the pilot's arms), warp→rudder
coupling matrix when enabled, control-surface deflection → geometry deltas consumed by
`fs-wing` (canard incidence; warp = antisymmetric washout twist of outer panels — the
actual mechanism, not an aileron abstraction).

#### 5.1.4 Pilot models (the assist ladder)

- `PilotDirect` — raw user input.
- `PilotWright` — user input + inner-loop pitch-rate damping: `δe += -k_q·q` filtered
  with ~180 ms first-order lag + deflection/rate saturation; gains fit ONCE against the
  historical porpoising amplitude/period (V-02) and then frozen with provenance.
- `PilotSAS` — modern rate+attitude feedback, clearly labeled non-historical.

#### 5.1.5 Launch systems & ground interaction

- **Rail:** dolly constraint (1-DOF along rail) until liftoff detection (net upward
  force > weight AND speed > threshold); rail friction; the wooden-rail bump spectrum
  as a small deterministic perturbation (seeded).
- **Catapult (Huffman):** weight-drop energy → rope tension profile → dolly force;
  configurable drop mass/height (defaults [V?] from 1904 accounts).
- **Skid–sand contact:** `fs-contact` capability-routed pairs; regularized Coulomb
  friction + plastic normal model for sand (tunable "sink" parameter); landing quality
  feeds the smoothness KPI and a damage flag (canard breakage as on flight 4 —
  threshold on nose-down impact energy [V?]).
- **Terrain queries:** heightfield + material map (sand/grass/marsh) sampled by contact
  and by `fs-atmo` (roughness z₀ varies by material).

#### 5.1.6 KPI accumulation (definitions in §9.2)

Accumulated *in the sim plane* (not derived from rendered frames): distance over ground,
airspeed/groundspeed extrema, altitude trace, load factor trace, pitch-rate RMS,
control-effort integrals, energy budget (engine work → prop losses → aero dissipation
→ kinetic/potential), stall-fraction time. Energy closure is the standing sanity
invariant (Euler-disc lesson: energy-balance gates catch nonsense early) — publish the
per-run closure residual in the results card.

### 5.2 `fs-wing` — lifting-surface aerodynamics (new crate, L3)

**The heart.** Three cooperating layers:

1. **Section layer (2-D, per strip):** lift/drag/moment coefficients
   `cl(α,δ), cd(α), cm(α)` for the cambered thin sections, from a blended model:
   thin-airfoil theory baseline + empirical stall blending (sigmoid to flat-plate
   post-stall) + Reynolds correction; **anchored to the Wrights' own 1901 wind-tunnel
   data** (their published tables for section #12 et al. — a delicious historical
   validation set) and to modern re-tests [V]. Unsteady augmentation per strip: 2-state
   Wagner approximation (R.T. Jones coefficients) for circulatory lag + apparent-mass
   terms + Küssner 2-state gust filter. This is the standard state-space unsteady
   aero used in real-time aeroelasticity; it is cheap (4 ODE states/strip) and captures
   the gust-response physics the smoothness KPI depends on.
2. **Planform layer (3-D):** lifting-line with prescribed elliptic-family influence OR
   (Tier B) a vortex-lattice bound sheet. Includes: biplane interference (gap-dependent
   Munk/Prandtl interference factor — validated against the classical biplane theory
   tables [V]), warp twist distribution, canard–wing mutual interference via the same
   induction machinery, and **ground effect by the image system** (every bound/wake
   element mirrored below the terrain plane; exact within the model class; degrades
   gracefully over sloped terrain by mirroring about the local tangent plane).
3. **Wake layer:** Tier A — prescribed helical/flat wake implied by lifting-line (no
   explicit particles) + analytic tip-vortex descriptors *for visualization only*
   (position/strength from Γ distribution + ground images). Tier B — explicit shed
   vorticity into an `fs-vpm`-kernel particle wake (desingularized Biot–Savart, the
   existing 2-D core generalized: see task E4.7 for the 3-D vortex-filament/particle
   extension of fs-vpm, upstreamed as `fs-vpm::three_d` rather than forked), with
   ground images, O(N) via cell-binned near-field + far-field truncation (fs-fmm is
   available if profiling demands it), wake feedback onto strips at 40 Hz, particle
   merge/decay for a capped population (~1,500–2,500).

Outputs per step: per-strip force/moment vectors (for 6-DOF and force-overlay viz),
bound Γ distribution, shed wake state (for viz + field service), and a diagnostics
block (per-strip α_eff, stall flag, local Re).

**Why not "just use fs-bem live"?** The Laplace panel solve (even FMM-accelerated) at
browser rates for an unsteady maneuvering aircraft is Tier C economics, and fs-bem's
own contract labels it inviscid screening. The strip/VLM+states model is the honest
real-time citizen; fs-bem is its offline referee. (fs-bem *is* used live in one place:
precomputing the biplane/canard interference tables per design change, ~100 ms
one-shot on the worker, cached by geometry hash.)

### 5.3 `fs-airscrew` — propulsion (new crate, L3)

- **Blade-element/momentum (BEMT)** with Prandtl tip/root loss, 20–30 radial stations,
  section polars from the same section layer (the Wrights' props were airfoils — their
  key insight); inflow from local airspeed + induced; yields thrust, torque, efficiency
  across the advance-ratio range including the static/rail-run regime (with the standard
  low-J momentum correction).
- **Engine:** torque-vs-rpm curve for the 1903 4-cylinder (peak ~12 hp fading with
  heat soak — a 60-s thermal derate curve reproduces the observed power fade [V?]);
  no throttle (run/off + spark retard as the only controls, per history); fuel mass
  depletion negligible for v1.
- **Drivetrain:** chain ratio, efficiency ~0.95 [V?], counter-rotation bookkeeping,
  optional chain-failure hostile scenario.
- Geometry defaults reconstructed from published measurements of the original props
  [V]; `prop_pitch_m` and `activity_factor` expose the user-facing "what if they'd
  carved differently" lever.
- Validation anchor: modern wind-tunnel data on 1903 prop reproductions (Wright
  Experience / AIAA papers) for η(J) [V]; BEMT curve must sit within declared bands.

### 5.4 `fs-atmo` — wind, turbulence, gusts (new crate, L2/L3)

- **Mean profile:** log-law `U(z) = (u*/κ) ln(z/z₀)` with site material map driving
  z₀ (sand vs grass vs water upwind); direction and reference speed from scenario.
- **Turbulence:** divergence-free synthetic turbulence via **curl of a vector
  potential** whose components are sums of spatially-frozen Fourier modes with von
  Kármán spectral amplitudes and anisotropy per boundary-layer similarity; frozen-field
  Taylor advection by the mean wind + slow phase evolution. Deterministic from a philox
  seed (`fs-rand`). Divergence-free **by construction** (the curl), which is what makes
  the divergence visualization an honest verification overlay. Length scales/intensities
  parameterized by (z, z₀, U_ref) per standard boundary-layer relations; a "Dec 17"
  preset tuned so the 10-m gust statistics match the historical range [V].
- **Gust events:** deterministic scheduled discrete gusts (1-cosine ramps) for
  challenge scenarios + Küssner-consistent sampling by `fs-wing`.
- **Thermals (Huffman, v1.5):** simple convective plume model (position-seeded rising
  columns with entrainment) — off for 1903 scenarios.
- **API:** `sample(x, t) -> (u, ∂u/∂x optional)` — analytic gradients available
  (sum-of-modes differentiates exactly), so the field service can render *exact*
  ∇, ∇·, ∇× of the ambient field, with discretization error only where the wake's
  particle field is involved.
- Validation: Tier D `fs-lbm` D3Q19 wind-over-terrain runs generate reference gust
  statistics (spectra, integral scales, gust factors at anemometer height) that the
  synthetic model must match within bands (V-04). This is exactly the fs-lbm/fs-aeroac
  "honest hybrid" division of labor already established in the workspace: LBM computes
  base-flow truth offline; a cheap model carries it to the consumer.

### 5.5 Field-sampling service (module in `fs-flyer`, exposed via wasm)

One function is the whole visualization backend:

```
sample_field(grid_spec, t) -> { u[], omega[], div[], grad_p_hat[], meta }
u(x) = U_mean(x,t) + u_turb(x,t) + u_bound(x) + u_wake(x) + u_images(x)
```

- Ambient parts: analytic (exact derivatives).
- Bound/wake parts: Biot–Savart over the current filament/particle set (SIMD kernel
  shared with the sim), derivatives by analytic kernel differentiation.
- Emits both the *exact model derivatives* and the *finite-difference-on-grid*
  versions when the verification overlay is on (the difference IS the lesson).
- Budgeted: grid capped (e.g., 32³ or 64×64 plane), refreshed at 15–30 Hz, runs in its
  own worker so it can never steal the sim budget; `fs-exec` cancellation honored
  mid-grid (the Cx checkpoint discipline, same as every FrankenSim kernel).

### 5.6 Sound (stretch, M6): browser WebAudio synthesis driven by physics state —
engine (4-cyl firing at rpm with exhaust character), prop blade-passing hum (2×rpm×2
blades), wind/airframe aeroacoustic noise shaped by airspeed and stall state, rail
rumble, skid-on-sand. Parameters derive from sim state each frame; synthesis stays in
an AudioWorklet. Labeled "sound design informed by physics, not a physical acoustic
claim" (the Euler-disc audio holds a higher bar — physical radiation — that v1 does not
attempt; fs-aeroac/fs-phs offer a v2 path to honest aeroacoustic estimates).

---

## 6. Numerics, Determinism, and the Execution Contract

1. **Determinism doctrine.** All transcendentals through `fs-math det::` (the
   powi/libm hazard class is documented workspace law; `check-powi` and `check-libm`
   lints already gate CI). All randomness through `fs-rand` philox with scenario-seeded
   keys. No wall-clock, no `HashMap` iteration order, no thread-order dependence in
   results (wake binning uses deterministic cell ordering; parallel reductions use
   fixed-shape trees). Target: **bit-identical replays on the same wasm build**, and
   cross-platform (aarch64/x86/wasm) agreement gated by the four-quadrant golden
   program (§10.5) rather than assumed.
2. **Fixed timestep, explicit schedule.** 120 Hz base, 240 Hz contact substep, 40 Hz
   wake update, 15–30 Hz field service — all integer-ratio locked to the sim clock, so
   a replay is a pure function of (config, seed, input trace).
3. **`fs-exec` everywhere.** Every wasm entry that can run long (field grids, sweeps,
   design-change re-derivations) takes a `Cx` with budget + cancellation and
   checkpoints at bounded intervals; the UI cancel button is therefore real, not
   cosmetic. ExecMode::Deterministic in all shipped paths.
4. **Typed refusals across the wasm boundary.** The JS API returns
   `{ok} | {refusal: {code, message, ranked_repairs}}` — mirroring the workspace's
   refusal style (and the Euler config admission lesson: refuse EARLY with repairs).
   No silent clamping: out-of-domain configs refuse with the domain stated.
5. **Units.** `fs-qty` types internally at module boundaries (a lesson paid for in
   the f85xj milliwatt-guard incident: unit bugs at seams are the expensive ones);
   the wasm boundary speaks SI doubles with documented units in the schema.
6. **Energy accounting invariant.** Per-step energy ledger (engine work in; aero
   dissipation, contact dissipation out; ΔKE+ΔPE) with closure residual tracked and
   exposed; regression tests pin the residual envelope per scenario. (Direct
   translation of the Euler-disc energy-balance ingest gate.)
7. **Identity.** Config, terrain, section tables, correction tables, and replay files
   all carry `fs-blake3` domain-hashed identities; new identity constants are
   registered in `identity-authorities.json` at introduction time (the m3h2e
   regression taught: unregistered identity constants across many lanes are how the
   identity gate rots — this project registers as it goes, never in a bulk sweep).

---

## 7. WASM Engineering & the Real-Time Budget

### 7.1 Build & packaging

- New crate `crates/fs-flyer-wasm` with its **own `[workspace]`** exactly like
  `fs-wasm` (documented rationale: browser builds decoupled from unrelated native WIP;
  the pattern already survives CI via wasm-pack with a locked nested Cargo.lock and a
  lock-drift gate — copy the `fs-wasm-build` lane's protections verbatim).
- Q(R1): extend `fs-wasm` vs. new crate. **Plan: new crate.** fs-wasm is a broad demo
  surface tier-organized by campaign; the Flyer needs an app-shaped, size-budgeted
  bundle (< 8 MB wasm gz target) with only the flyer cone. Shared deps come from the
  same paths; no duplication of infra beyond a small build script.
- `wasm-pack --target web`, `wasm-opt -O3`, SIMD128 on with a no-SIMD fallback build
  (two artifacts; loader picks by feature detection).
- Threads build (SharedArrayBuffer + wasm threads) as the *enhanced* artifact; the
  single-threaded artifact is the baseline. Both produced by CI; both self-tested.

### 7.2 The frame budget (Tier A, mid-range laptop, single sim worker)

| Module | Per step (120 Hz) | Per frame @60fps (2 steps) |
|---|---|---|
| fs-atmo sampling (strips+props+probes ~80 pts) | 0.10 ms | 0.20 |
| fs-wing strips (36 strips × states + 3-D solve) | 0.25 ms | 0.50 |
| fs-airscrew (2 props × 24 stations) | 0.06 ms | 0.12 |
| fs-contact + rail | 0.03 ms | 0.06 |
| fs-mbd + fs-time integrate | 0.02 ms | 0.04 |
| KPIs, ring buffer, bookkeeping | 0.05 ms | 0.10 |
| **Sim worker total** | **~0.5 ms** | **~1.0 ms** |

Tier B adds: wake advance+shed @40 Hz with N≈2,000 particles, cell-binned near field:
~2.5 ms per wake step → ~0.8 ms/frame amortized; strip-feedback induction 36×N with
SIMD ≈ 0.6 ms per wake step. Field service (separate worker): 32³ grid × (analytic +
N-particle kernel) ≈ 15–40 ms per refresh at 15 Hz — fits its worker with headroom.
Main-thread three.js budget: 6–10 ms (scene of ~200k tris + instanced glyphs).
**Verdict: comfortably real-time; the risk lives in the wake and field service, both
rate-decoupled and degradable (particle cap, grid cap) without touching flight physics.**

### 7.3 Interop contract

- One `SharedArrayBuffer` ring (triple-buffered) of POD state structs (fixed layout,
  versioned header, ~2 KB/snapshot incl. 36 strip force vectors + instrument block).
- Field buffers: preallocated Float32Array (u:3, ω:3, div:1, |∇u|:1 per cell), written
  by field worker, read by render as raw attribute buffers for instancing (zero-copy
  into GPU upload).
- Command channel: postMessage JSON for config changes (rare) — every config change
  re-derives tables (~50–150 ms, off-thread, with a progress event) and mints a new
  config identity.
- No per-frame JSON anywhere.

### 7.4 Degradation ladder & self-test

Startup self-test (mirrors the Euler E2E runner's `--self-test` spirit, in-browser):
compute a 1-second canonical scenario headless, compare against an embedded golden
state hash for this build; report tier availability (threads/SIMD), chosen tier, and
the golden verdict on the About panel. A mismatch shows a visible "determinism
self-test failed" badge (never silently ship wrong physics).

---

## 8. Rendering & UX (three.js)

### 8.1 App shell

`apps/wright-flyer/` — TypeScript + Vite + three.js (pinned version; renderer wrapped
in a thin interface so a future WebGPU/TSL migration or in-house renderer swap is a
module change, not a rewrite). State management deliberately minimal (the sim owns
truth; UI state is view-only). No backend required to fly (static hosting + correct
COOP/COEP headers); optional tiny share service later.

### 8.2 The aircraft asset

- Source: a high-detail Wright Flyer 3-D model. **Task E2.1 vets candidates**:
  Smithsonian's 3-D digitization of the actual 1903 Flyer (their open program has
  released Flyer scans; license CC0-class [V]), NASA model releases, or a commissioned
  artist model. License must permit web distribution; provenance recorded in-app.
- Pipeline: source → glTF 2.0 → decimated LODs (hero ~150k tris, mid, far) → KTX2
  textures → **rig**: warp deformation (skeleton or morph targets driven by the warp
  state — the fabric visibly twists, which is historically wonderful), canard/rudder
  pivots, prop rotation with motion-blur impostor discs at speed, chain animation
  (texture scroll), pilot figure with hip-cradle pose coupling.
- The *physics* never reads the visual mesh: `fs-wing` geometry is the parametric
  planform. A calibration overlay mode draws the physics panels/strips on top of the
  visual model to make any mismatch visible (and it becomes the debug view).
- When the user drags design sliders, the visual model morphs procedurally (span/chord
  scaling per wing bay, strut/wire re-solve as a tiny geometric constraint pass) —
  imperfect visually at extremes, clearly flagged ("schematic preview beyond ±25%").

### 8.3 Terrain, sky, environment

- Heightfield terrain (2×2 km, ~1–2 m res) with material splat (sand/grass/scrub);
  Kill Devil Hills 1903 layout: rail, hangar+workshop sheds, camera tripod at the
  famous spot (the user can take the John T. Daniels photo themselves — instant-photo
  mode with period grain, a shareable moment), Big Kill Devil Hill in the middle
  distance. Huffman Prairie: pasture, fence lines, derrick, the honey locust tree.
- Sky: physically-plausible sky model with date/time presets (Dec 17 overcast-bright,
  ~34°F — cold clear-grey North Carolina winter light [V]); sun position from date/site.
- Wind made visible even with viz off: sand streamers near the ground, grass/fabric
  flutter amplitude driven by local `fs-atmo` samples, wind sock (they had a flag/
  anemometer; a period-correct cue), streaming camp smoke.
- Ocean/sound in the distance (flat planes with shader waves; not simulated).

### 8.4 Field visualization implementation

- **Glyphs:** InstancedMesh arrows (up to ~30k), attributes streamed from the field
  buffer; length/color transfer functions with legend; probe-box gizmo (drag/scale).
- **Streamribbons:** GPU-side integration in a shader over the sampled grid texture
  (RK2 in texture space) for dense ephemeral lines + a CPU "hero rake" of 20–50
  ribbons integrated accurately in the field worker for teaching views.
- **Vorticity:** wake particles/filaments rendered directly (they carry Γ) as ribbons
  with age-fade; ambient |ω| as volume slices through the grid.
- **Div verification overlay:** grid cells colored by |∇·u| (finite-difference) with
  the analytic-vs-FD toggle described in §5.5.
- **Force overlay:** per-strip arrows from the state snapshot (exact sim values).
- All overlays render into a separate pass composited with depth-aware transparency so
  the aircraft stays legible.

### 8.5 HUD & instruments

Period instruments (the Flyer carried: anemometer, stopwatch, engine-revolution
counter — that triad *is* the historical KPI set [V]) rendered as a wooden panel;
modern overlay (airspeed, altitude, α, load factor, pitch-rate, thrust, L/D, energy
ledger strip) toggleable. Results card after each run; design-diff card comparing two
configs KPI-by-KPI.

### 8.6 Replay UI

Timeline scrubber with event ticks (liftoff, gusts, stall flags, touchdown), camera
presets (chase, wing-tip, Daniels-tripod, onboard prone-pilot view, free orbit),
A/B ghost mode (two replays superimposed with translucent second aircraft — THE tool
for showing what a design change did), and export (replay file download; cinematic
export instructions pointing at the native pipeline).

---

## 9. Configuration Space, Experiments & KPIs

### 9.1 Config schema & governance

`FlyerScenario` = { design: FlyerDesign, site, weather (mean wind, seed, gust events),
launch (rail/catapult params), pilot mode, sim tier }. Serialized as a versioned JSON
schema following the workspace's sidecar-not-IR and schema-freeze doctrines: the
schema version is a registered constant; additive evolution only after freeze; the
schema file + example scenarios are tracked sources (and therefore covered by
`generate-source-manifest` and the identity gate). Preset library: `dec17-flight1..4`,
`huffman-1904-catapult`, `flyer3-1905`, plus challenge presets.

### 9.2 KPI definitions (exact, so beads can pin tests)

| KPI | Definition |
|---|---|
| Distance | ground-track arc length from liftoff to first touchdown (m and period-correct feet) |
| Air distance | ∫ airspeed dt (the Wrights reported both notions' tension on Dec 17) |
| Duration | liftoff→touchdown |
| Max/mean airspeed | over airborne phase |
| Lift margin | min over airborne phase of (CL,max,local − CL,required)/CL,max — how close to stall |
| Smoothness index | 1 / (1 + w₁·RMS(a_z′) + w₂·RMS(q)) over airborne phase; weights fixed & documented |
| Control effort | ∫(|δ̇e| + |δ̇w|)dt |
| Static margin | (x_np − x_cg)/c̄ from live linearization (negative for the 1903 Flyer!) |
| Short-period / phugoid | eigenvalues of the linearized longitudinal model at trim (displayed on the root-locus mini-plot as sliders move — the single best "see the physics" widget for engineers) |
| Energy closure | |ledger residual| / engine work |
| Structural margin | max load factor vs limit; refusal-flag if exceeded (wing "failure" ends flight in v1 with an honest "structural limit exceeded" card, no debris physics) |

### 9.3 Experiment engine

- **Live linearization:** central-difference Jacobian of the reduced model at trim
  (cheap: model is smooth by construction away from stall) → stability derivatives
  (CLα, Cmα, Cmq …) shown in the engineer panel with the historical estimates for
  comparison [V].
- **Sweeps:** worker-pool batch runs over 1-D/2-D config grids (headless sim, no
  render), progress-streamed, results plotted (canvas) and exportable CSV. Determinism
  makes sweep results cacheable by (config-hash, seed).
- **Optimization hook (v1.5):** `fs-bo`/`fs-dfo` compiled in fs-wasm already — "let
  the optimizer find the best camber+AR under a structural constraint" is a one-worker
  task and a spectacular demo ("the browser just rediscovered wing design").

---

## 10. Validation & Evidence Program

The program that keeps "without dumbing down the physics" an audited claim.

### 10.1 Anchor datasets

| ID | Dataset | Validates |
|---|---|---|
| A1 | Wright 1901 wind-tunnel tables (published) | section model cl/cd trends |
| A2 | Modern re-tests of Wright sections & full-scale 1903 replica (AIAA Wright Flyer Project; Langley full-scale tunnel campaign; Wright Experience prop data) [V] | section+planform force model, prop η(J) |
| A3 | Dec 17 flight records (4 flights + wind) | end-to-end KPI bands |
| A4 | Culick et al. stability analyses of the 1903 Flyer | static margin sign/magnitude, short-period character |
| A5 | Classical biplane interference theory tables | gap-interference implementation |
| A6 | Boundary-layer meteorology relations + fs-lbm Tier D runs | fs-atmo statistics |

### 10.2 Verification (math right) — per-crate batteries

- fs-wing: elliptic-wing analytic induced drag (e = 1 recovery), 2-D thin-airfoil
  limit, Wagner/Küssner step responses vs published curves, image-system ground-effect
  vs classical curves, biplane factor vs A5.
- fs-airscrew: momentum-theory limits, Prandtl-loss asymptotics, energy consistency
  (P_shaft = T·V + P_induced + P_profile).
- fs-atmo: spectral content vs prescription (periodogram test), exact ∇·u = 0 at
  machine precision on analytic samples, seed determinism.
- fs-flyer: energy-ledger closure on ballistic + powered fixtures; contact
  restitution/friction fixtures; replay bit-identity.
- All batteries follow workspace law: typed refusals tested at cap AND cap+1
  (CanonicalLimits lesson), no vacuous limit checks (both sides of every comparison
  computed), falsifier-style negative tests for each gate (sum-tests-are-blind and
  metamorphic-blindness lessons from the workspace memory apply — per-strip oracles,
  not just totals).

### 10.3 Validation (model right) — the V-cases

| ID | Case | Pass band (initial; tightened by evidence) |
|---|---|---|
| V-01 | Dec 17 flight 1 replay (authentic pilot Wright-hands, recorded wind band) | distance 120 ft × [0.6, 1.6]; duration 12 s × [0.7, 1.4] |
| V-02 | Porpoising signature | pitch oscillation period 2–4 s, growing without pilot damping; bounded with PilotWright |
| V-03 | Prop η(J) vs A2 | within ±8% over J ∈ [0.4, 0.9] [V bands with data] |
| V-04 | fs-atmo vs Tier D LBM gust stats | spectra within factor 2 across the energy-containing range; gust factor ±20% |
| V-05 | dt/tier convergence | KPIs stable within 2% between 120→240 Hz and Tier A→B on reference scenarios |
| V-06 | Ground effect | induced-drag reduction vs height matches image-theory curve within 5% (verification-grade) and Tier C panel re-run within 15% |
| V-07 | 852-ft flight 4 envelope | reachable with authentic config + recorded wind + Wright-hands pilot |

Receipts land in the `fs-vvreg`/vv-scorecard machinery the workspace already runs
(corpus edits regenerate the scorecard; the check-all gate keeps it honest).

### 10.4 The validated envelope

The set of configs within declared distance (in normalized design space) of Tier-C-
validated pins. Inside: UI badge `Validated (Estimated dynamics, validated forces)`.
Outside: `Estimated — outside validated envelope`, with the nearest pin distance shown.
This is the product-facing face of the evidence doctrine and must be spec'd precisely
in E8 (envelope metric, pin set, refresh protocol).

### 10.5 Determinism goldens (four-quadrant, wasm-extended)

Golden replay state-hashes for canonical scenarios across {aarch64-native, x86-native,
wasm-in-node} × {debug, release} — extending the workspace's four-quadrant golden
doctrine with a wasm column. Gated in CI next to the existing golden-couplings
machinery; the golden-bump protocol (committed tree, both modes, plausible cause,
same-commit goldens) applies verbatim.

### 10.6 E2E runner & hostile twins

`scripts/ci/e2e_wright_flyer.sh` cloned from the (just-hardened) Euler cinematic
runner pattern: `--list/--check/--self-test/--run smoke/--negative CASE/--replay`,
bounded JSONL logging contract (same schema family), and hostile twins: config
tampered vs identity, replay input-trace truncation, seed mismatch, energy-ledger
violation injection, stale correction-table identity, KPI card vs recomputed KPIs
mismatch, wasm/native golden divergence, terrain-hash drift. Runner reuses the
production CLI (`fs-flyer` native binary) — never parallel logic.

---

## 11. Crate Reuse Matrix & New Crates

### 11.1 Existing crates leveraged (verified present)

| Crate | Role here | Notes |
|---|---|---|
| fs-wasm | pattern + infra precedent for browser workspace, CI lane, wasm-pack recipe | fs-flyer-wasm copies its protections |
| fs-bem (+fs-fmm) | Tier C offline force referee; live one-shot interference tables | already wasm32-proven |
| fs-vpm | Biot–Savart core for wake; extend to 3-D filaments/particles upstream | 2-D today (verified header) |
| fs-lbm | Tier D wind-over-terrain truth runs | D3Q19 + boundaries exist |
| fs-mbd | 6-DOF rigid body, canonical quaternions | unconstrained core, by design |
| fs-contact | skid/rail/sand contact routing | |
| fs-time | symplectic/structure-preserving integrators | |
| fs-exec | Cx, budgets, cancellation, deterministic mode | wasm-proven via fs-wasm |
| fs-math | det:: transcendentals (libm doctrine) | determinism backbone |
| fs-rand | philox streams | |
| fs-blake3 | identity domains for config/replay/tables | register identities on introduction |
| fs-qty | unit-checked quantities at seams | |
| fs-scenario | scenario schema conventions | |
| fs-simd | SIMD tiers for Biot–Savart/BEMT kernels | Tier 0 scalar = referee |
| fs-viz | field-viz primitive algorithms w/ analytic ground truth | CPU side of streamlines etc. |
| fs-uq / fs-surrogate / fs-bo / fs-dfo | sweeps, surrogates, optimization hook | already in fs-wasm |
| fs-vvreg (+ vv-scorecard) | validation receipts + reporting | standing infra |
| fs-render + euler mux adapter + cinematic runner | offline cinematic export of replays | just hardened (h7xu5) |
| fs-evidence | evidence colors / no-claims plumbing | |
| fs-ornith | prior art: parameterized-aircraft campaign staging | pattern, some code reuse in screening stage |

### 11.2 New crates (kept to five)

| Crate | Layer | One-line contract |
|---|---|---|
| fs-atmo | L2/L3 | deterministic boundary-layer wind: log-law mean + div-free spectral turbulence + gust events; analytic derivatives; NO acoustic/thermo claims |
| fs-wing | L3 | strip/lifting-line/VLM unsteady lifting-surface aero with ground images and biplane interference; per-strip typed diagnostics; inviscid+empirical honesty labels |
| fs-airscrew | L3 | BEMT propeller + engine torque + drivetrain; energy-consistent; static-thrust caveats declared |
| fs-flyer | L4 | the aircraft: parametric airframe, 6-DOF assembly, pilots, launch, terrain contact, KPIs, scenario schema, field service |
| fs-flyer-wasm | L6 | own-workspace wasm binding: sim loop, field service, sweeps, replay, typed-refusal JS API |

Each ships CONTRACT.md, no-claims block, refusal vocabulary, and registered identity
constants from day one (workspace law; cheaper at birth than at audit).

---

## 12. Milestones & Dependency-Aware Task Graph

Conventions: `E#.#` epics/tasks, `→` blocks. DONE-WHEN clauses are bead-ready.
Estimates assume the established FrankenSim agent workflow (code-first, batteries
executed, E2E runner discipline).

### E0 — Program setup
- **E0.1** Program root bead + this plan converted to beads (after review rounds).
- **E0.2** `apps/wright-flyer` scaffold (Vite+TS+three.js, COOP/COEP dev server,
  CI lint/build lane). DONE-WHEN: blank scene at 60 fps deployed to a static host.
- **E0.3** `fs-flyer-wasm` scaffold on the fs-wasm pattern (own workspace, nested
  lock, wasm-pack CI lane cloned incl. lock-drift gate, hello-kernel exposed and
  called from E0.2's page). → blocks all wasm integration tasks.

### E1 — Historical grounding & data
- **E1.1** Source dossier: assemble A1–A6 datasets, licenses, citations. → E4, E10.
- **E1.2** Verify every [V]/[V?] number in §3; produce `flyer-reference.json`
  (tracked, identity-hashed) as the single source of defaults. → E3.1, E4.
- **E1.3** Terrain data: DEM acquisition both sites, 1903/1904 historical adjustment
  notes, heightfield + material map assets with provenance file. → E2.3, E5.4.

### E2 — Assets & rendering foundation
- **E2.1** Flyer 3-D model: license vetting, acquisition decision. → E2.2.
- **E2.2** Asset pipeline: glTF, LODs, KTX2, rig (warp morph, canard/rudder pivots,
  props, pilot). DONE-WHEN: rigged model animates from a scripted state file in the
  app at 60 fps.
- **E2.3** Terrain+sky+environment scene for Kill Devil Hills (rail, sheds, dune,
  ambient wind cues driven by stub data). DONE-WHEN: the §2.1 "arrival" shot exists.
- **E2.4** Camera system + input mapping + HUD skeleton.

### E3 — Simulation spine (Tier A minimal)
- **E3.1** `fs-flyer` crate: FlyerDesign schema + admission refusals + mass/inertia
  build-up + derived-quantity panel math. Depends E1.2. DONE-WHEN: battery pins
  reference mass/CG/inertia vs dossier values.
- **E3.2** 6-DOF core on fs-mbd/fs-time with gravity + placeholder aero; fixed-dt
  loop; state ring buffer; replay record/playback bit-identity test. → E3.4, E6.1.
- **E3.3** `fs-atmo` v0: log-law mean + frozen div-free turbulence + seeds; battery:
  spectra, div=0, determinism. (Parallel to E3.2.)
- **E3.4** Rail launch + skid contact + terrain heightfield queries (fs-contact).
  DONE-WHEN: dolly run, liftoff hand-off, sliding landing all stable at 240 Hz substep.

### E4 — Aerodynamics & propulsion (the physics heart)
- **E4.1** Section layer: cambered thin-airfoil + stall blend + Re correction,
  anchored to A1/A2; battery incl. falsifier (perturbed-camber must move cl per
  theory). Depends E1.1/E1.2.
- **E4.2** Strip lifting-line w/ biplane interference + canard/rudder surfaces;
  analytic verification battery (elliptic wing, biplane tables). Depends E4.1.
- **E4.3** Unsteady states (Wagner/Küssner) per strip + gust intake from fs-atmo.
  Depends E4.2, E3.3.
- **E4.4** Ground-effect image system (flat + local-tangent-plane). Depends E4.2.
- **E4.5** `fs-airscrew` BEMT + engine curve + drivetrain; battery vs momentum
  limits + A2 η(J). Depends E4.1 (shares section machinery), E1.2.
- **E4.6** Integration: full force build-up into E3.2's 6-DOF; trim solver; live
  linearization + stability derivatives. DONE-WHEN: V-02 porpoising signature
  reproduced (period band), static margin negative per A4.
- **E4.7** fs-vpm 3-D extension (filaments/particles, desingularized kernel, cell
  binning) upstreamed in fs-vpm; Tier B wake shed/feedback in fs-wing @40 Hz with
  particle cap + merge. Depends E4.2; SIMD via fs-simd. DONE-WHEN: Tier A vs Tier B
  KPI deltas within V-05 band on reference scenarios.
- **E4.8** Interference-table one-shot via fs-bem on design change (cached by
  geometry hash). Depends E4.2.

### E5 — Browser integration (playable alpha)
- **E5.1** fs-flyer-wasm API v1: init(scenario), step-loop in sim worker, state ring,
  control input, refusal envelope. Depends E3.*, E4.6, E0.3.
- **E5.2** three.js consumes real state: aircraft pose, control-surface morphs, prop
  spin; instruments live. Depends E5.1, E2.2. **MILESTONE: FIRST FLYABLE BUILD.**
- **E5.3** Pilot-assist ladder + gamepad/mouse tuning passes (feel work). Depends E5.2.
- **E5.4** Both sites + launch options (rail headwind / catapult). Depends E3.4, E1.3.
- **E5.5** Results card + KPIs + historical comparison presets. Depends E5.2.

### E6 — Determinism, replay, E2E harness
- **E6.1** Replay files (config identity + seed + input trace), scrubber UI, ghost
  A/B mode. Depends E5.2, E3.2.
- **E6.2** In-browser startup self-test + golden hash; four-quadrant+wasm golden CI
  lane. Depends E5.1.
- **E6.3** `e2e_wright_flyer.sh` runner + hostile twins + JSONL logging (clone the
  Euler runner skeleton). Depends E6.1, native `fs-flyer` CLI (small task inside).

### E7 — Field visualization (the wow)
- **E7.1** Field service in fs-flyer + wasm API (grids, probes, exact-vs-FD
  derivative duals). Depends E4.3 (ambient) and E4.7 (wake) for full content; ships
  ambient-only first.
- **E7.2** Glyph + ribbon + vorticity + divergence-overlay renderers (instancing,
  transfer functions, legends, probe gizmos). Depends E7.1, E2.3.
- **E7.3** Force overlay + pressure bands + probes with strip-charts. Depends E5.2.
- **E7.4** Lesson mode scaffolding (curated overlay scripts). Depends E7.2/E7.3.

### E8 — Experiments, sweeps, evidence surfacing
- **E8.1** Worker-pool sweep engine + plots + CSV. Depends E5.1.
- **E8.2** Design panel v2: root-locus mini-plot, polar redraw, design-diff cards.
  Depends E4.6, E5.5.
- **E8.3** Validated-envelope spec + UI evidence badges wired to receipts. Depends E10.2.
- **E8.4** (v1.5) fs-bo optimization demo. Depends E8.1.

### E9 — Sound & polish (stretch)
- **E9.1** AudioWorklet engine/prop/wind synthesis from state. Depends E5.2.
- **E9.2** Instant-photo mode, challenges rail, onboarding flow. Depends E5.5, E7.2.

### E10 — Truth plane & validation program
- **E10.1** Tier C referee harness: batch re-run of pinned configs through
  fs-bem/high-res models; correction tables + residual receipts. Depends E4.*, E1.1.
- **E10.2** V-01…V-07 executed and recorded in fs-vvreg/vv-scorecard; pass bands
  ratified. Depends E10.1, E6.3. **MILESTONE: EVIDENCE-BADGED BETA.**
- **E10.3** Tier D fs-lbm wind-over-terrain runs → fs-atmo band validation (V-04).
  Depends E1.3. (Perf note: runs on the Linux perf hosts per workspace practice.)
- **E10.4** Cinematic export path: replay → native trajectory → fs-render scene
  bridge → EXR/ProRes via existing mux adapter; one hero clip produced. Depends E6.1;
  reuses h7xu5 machinery.

### Critical path

E1.2 → E4.1 → E4.2 → E4.6 → E5.1 → E5.2 (first flyable) → E6.x/E7.x in parallel →
E10.2 (beta). Terrain/assets (E1.3, E2.x) parallel the physics spine. The wake tier
(E4.7) and field viz (E7) are protected: they gate the *wow*, not the *flyable*.

---

## 13. Risks & Mitigations

| # | Risk | Likelihood | Mitigation |
|---|---|---|---|
| 1 | Aero model misses the porpoising signature (canard modeling subtle) | med | V-02 gate early (E4.6 DONE-WHEN); A4 stability-derivative cross-check; Tier C referee before feel-tuning |
| 2 | Wake/field costs blow the frame budget on low-end devices | med | tier ladder + hard caps + rate decoupling; ship Tier A as the contract, Tier B as enhancement |
| 3 | SharedArrayBuffer header requirements on hosting | high (known) | dual artifacts + degraded single-thread mode from day 1 (E0.3), not retrofitted |
| 4 | 3-D model licensing (Smithsonian terms, artist rights) | med | E2.1 is a *blocking* vetting task; fallback commissioned low-poly + procedural rig |
| 5 | Historical numbers contested (prop η, engine hp fade, canard area) | med | [V] discipline: dossier task E1.2, tunable-with-provenance defaults, bands not points in validation |
| 6 | wasm/native numerical divergence breaks goldens | med | det:: doctrine + four-quadrant+wasm goldens early (E6.2), before physics churn |
| 7 | Scope creep toward general flight sim | high | §1.4 non-goals; new-aircraft requests become v2 beads, never v1 scope |
| 8 | Instability makes the game feel "broken" to casual users | high | assist ladder default = Wright-hands, authentic opt-in; onboarding explains WHY it porpoises (the instability is the story, told as such) |
| 9 | Determinism vs three.js frame jitter confusion | low | fixed-step sim + interpolation; replay hashes computed in sim plane only |
| 10 | fs-vpm 3-D extension underestimated | med | Tier A ships without it; E4.7 has its own falsifier battery and can slip without gating first-flyable |

---

## 14. Open Questions for Review Round 1

1. Tier B wake: filaments (connectivity, cheaper visuals) vs particles (simpler
   merging) — plan says particles with filament rendering for tips; challenge this.
2. Should the section layer use a small neural/spline fit to A2 data (fs-surrogate)
   instead of analytic+blend? (Provenance and out-of-domain refusals get harder;
   accuracy gets better.)
3. Biplane interference: full VLM in Tier A (costs ~0.3 ms more) vs corrected
   lifting-line — is the simpler model defensible for the canard-wing coupling that
   drives risk #1?
4. Pilot latency model: is 180 ms + saturation enough to reproduce flight-4's growing
   oscillation, or do we need a proper crossover-model human pilot?
5. Terrain: is 2×2 km enough for the 1905 24-mile circling flight (Huffman scenario) —
   or do we tile/stream and accept it in v1?
6. Evidence UX: is the three-color badge legible to lay users, or does it need a
   plain-language layer ("physics checked against wind-tunnel data for this design")?
7. Replay portability across app versions: freeze policy for the replay schema
   (schema-freeze doctrine says decide *now*).
8. Should Tier D LBM ever run *with* the aircraft immersed (moving IB boundary —
   large fs-lbm extension) or stay terrain-only in v1? (Plan says terrain-only.)
9. WebGPU: three.js WebGPURenderer is maturing — pin WebGL2 for v1 with the renderer
   interface, or bet on WebGPU now for the field viz compute passes?
10. Audio scope: is E9.1 sound-design synthesis enough, or does the Euler-disc
    physical-audio lineage demand a v1.5 fs-phs/fs-aeroac-grounded engine note?

---

## 15. Appendices

### A. Model equations (implementation-normative)

**A.1 Strip force build-up.** For strip i with local chord c, span width Δy, unit
span direction; local flow `u_loc = u_atmo(x_i) + u_wake(x_i) + u_images(x_i) −
(v_cg + ω×r_i)`; α_eff from u_loc in the section frame incl. warp twist and canard
incidence; circulatory lift via unsteady states (A.2), profile drag from cd(α_eff,Re),
induced handled by the planform layer (never double-counted — the lifting-line's
induced α is subtracted from the section α input; document the exact bookkeeping in
fs-wing's CONTRACT to keep the classic double-count bug testable).

**A.2 Unsteady states (per strip).** R.T. Jones' Wagner approximation
Φ(s) ≈ 1 − 0.165e^(−0.0455s) − 0.335e^(−0.30s), s = 2Ut/c, realized as 2 LTI states;
Küssner ψ(s) ≈ 1 − 0.5e^(−0.13s) − 0.5e^(−s) likewise; apparent-mass terms
L_am = πρ(c/2)²(ḧ + U α̇ − (c/4)α̈ term per convention chosen and documented). All
constants dimensionless & sourced (Fung / Leishman) — cite in code.

**A.3 Lifting-line with images & biplane.** Bound Γ(y) on N stations per wing;
induced velocity from bound+trailing horseshoe system of BOTH wings + all images
(ground plane at local terrain tangent). Solve the N×N (≤ 80×80) linear system per
step (fs-la dense solve, reuse factorization when geometry unchanged; rank-1-ish
updates on control deflection are an optimization task).

**A.4 BEMT.** Standard per-annulus momentum/blade-element iteration with Prandtl
factor F = (2/π)acos(e^(−f)); low-J regime via momentum correction (Glauert empirical
region); convergence guarded with typed refusal on non-convergence (never NaN).

**A.5 Turbulence synthesis.** u_turb = ∇×A, A_j(x,t) = Σ_k a_{jk} sin(k·x + ω_k t +
φ_{jk}) with amplitudes shaped so the resulting velocity spectrum matches von Kármán
(σ_u, L_u anisotropic per height); modes ~64–128; φ from philox; ω_k = advective +
slow decorrelation. Exact curl/div/grad available analytically (differentiate the
sum) — the basis of the honest viz claims.

**A.6 Ground effect image bookkeeping.** Every vortex element (bound, trailing,
wake particle) has a mirrored partner with sign-flipped normal circulation across
the local ground plane; images regenerate when the reference tangent plane tilts
by > ε (sloped terrain), with the approximation documented.

**A.7 Longitudinal linearization.** Standard (u, w, q, θ) small-perturbation system
about trim; central differences on the full nonlinear RHS; eigenvalues → short
period/phugoid; x_np from ∂Cm/∂CL root-solve.

### B. Performance math (Tier B wake)

N = 2,000 particles; shed 36/step at 40 Hz → cap reached in ~1.4 s of wake age →
merge policy (Γ-weighted pairwise within cells, oldest-first) holds N. Advance:
mutual induction via 32³ cell binning ⇒ ~N·k evals, k≈60 neighbors ⇒ 120k kernel
evals + strip feedback 36×2,000 = 72k ⇒ ~192k evals/wake-step ≈ 5.8 Mevals/s at
40 Hz ⇒ with 4-wide SIMD ≈ 25 flops/eval ⇒ ~36 MFLOP/s — comfortable. Field grid
32³ = 33k points × (2,000 particles via binned far-field truncation + analytic
ambient) at 15 Hz ⇒ budget ~40 ms/refresh in its own worker — fits.

### C. Historical dossier seed (to be completed in E1)

Primary: Wright brothers' diaries/letters (LOC), 1903 photographs (Daniels plate),
McFarland's *The Papers of Wilbur and Orville Wright*. Secondary: AIAA Wright Flyer
Project publications (Culick et al. stability analyses), Langley full-scale tunnel
replica test reports, Wright Experience propeller reconstruction data, NPS Wright
Brothers National Memorial site documentation, Cerwin/Jakab historical studies.
Terrain: USGS 3DEP, NOAA shoreline history, NPS historical base maps.

### D. Glossary

Advance ratio J = V/(nD). Static margin: (x_np−x_cg)/c̄, negative = unstable.
Wagner function: circulatory lift response to step α. Küssner function: response to
sharp-edged gust. BEMT: blade-element/momentum theory. VLM: vortex-lattice method.
Image system: mirrored singularities enforcing ground tangency. Porpoising: coupled
pilot-aircraft pitch oscillation. Warp: the Wrights' roll control via wing twist.

---

## Review round log

| Round | Reviewer | Date | Disposition |
|---|---|---|---|
| 0 | (this draft) NobleLion / Claude | 2026-08-16 | initial comprehensive plan |
| 1 | — | — | pending (GPT Pro Extended Reasoning, EXACT PROMPT) |
| 2–4+ | — | — | pending; convert to beads only at steady-state |
