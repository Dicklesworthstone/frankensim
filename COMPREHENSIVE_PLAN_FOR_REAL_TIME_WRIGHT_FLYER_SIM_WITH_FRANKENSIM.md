# COMPREHENSIVE PLAN: Real-Time Wright Flyer Simulation with FrankenSim

**Working title:** *First Flight — Kitty Hawk, December 17, 1903*
**Document status:** Planning-workflow ROUND 1 (external review integrated; ≥3 further
rounds required before beads conversion).
**Process:** This document follows `/planning-workflow`. The Round 1 review (GPT Pro
Extended Reasoning) accepted the architecture and required a major physics-and-
validation revision, integrated throughout this version. A "Review round log" at the
end tracks rounds and dispositions.
**Repo:** `frankensim` (156-crate Rust workspace, layers L0–L6, Franken-only
dependencies, determinism doctrine, typed refusals, evidence colors
Verified/Validated/Estimated).

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
14. [Open Questions for Review Round 2](#14-open-questions-for-review-round-2)
15. [Appendices](#15-appendices)

---

## 1. Vision

### 1.1 What we are building

A browser-native, real-time, physically honest simulation of the 1903 Wright Flyer's
first flights — the airplane, the air, the ground, the *controls*, and the pilot —
rendered beautifully with three.js, computed by FrankenSim physics compiled to
WebAssembly, with the invisible made visible: at any moment the user can flip on the
**wind vector field** and watch the turbulent boundary layer roll over the Kill Devil
Hills sand, see the wingtip vortex sheets peel off and sink toward the ground, and
overlay **divergence, gradient, and curl** of the velocity field in 3-D as living
mathematical objects rather than textbook symbols.

Every load-bearing number of the real Flyer is dialed in — span, chord, camber, the
canard's geometry *and its hinge mechanics*, gross weight, engine power, propeller
geometry, the December headwind — and every one of them is a **user-adjustable
configuration** with model-determined consequence. Vary canard area, balance, arm,
and CG and watch static stability, free-control stability, trim authority, and
control force trade against one another — *the model decides* which change helps and
which makes the aircraft uncontrollable. Stretch the aspect ratio and inspect the
resulting changes in induced drag, mass, rigging load, stall distribution, and
structural utilization. Take away the headwind and discover, as the Wrights did,
that the launch rail is suddenly too short. Then read the numbers: separation
margin, distance flown, peak speed, ride quality, control activity — the same
quantities of interest the Wrights argued about in letters, now live on screen with
their evidence pedigree attached.

> Product-copy rule (Round 1): marketing and lesson text never predetermines a
> physics outcome. Tradeoffs are advertised; results are computed. No "and you'll go
> 240 feet" in copy — the delta card reports what the model actually produced and
> decomposes why.

### 1.2 Why this project, and why FrankenSim

The Euler-disc E2E pipeline proved the pattern: one unified simulation system
produced mechanics, spectral-rendered video, and physically derived audio, all from a
single identity-tracked configuration with honest evidence labels. This project
scales the same pattern from a spinning disc on a table to an aircraft in a turbulent
atmosphere — and moves the interactive tier into the browser.

This is not a from-scratch flight sim. The central engineering claim, grounded in the
current tree (verified 2026-08-16, §11.3):

- **`fs-wasm` already ships a large slice of FrankenSim to wasm32** via its own
  decoupled workspace and `wasm-pack` CI builds; its dependency list already includes
  `fs-bem`, `fs-vpm`, `fs-lbm`, `fs-exec`, `fs-alloc`, `fs-viz`, `fs-render`,
  `fs-scenario`, `fs-uq`, `fs-qty`, `fs-la`, `fs-material`. The browser toolchain is
  solved infrastructure.
- **`fs-ornith` already flies a parameterized aircraft in the browser** at screening
  fidelity (fs-bem panel aero + fs-vpm wake metric in a staged campaign).
- **`fs-mbd`** (deterministic quaternion rigid body, leaf crate) and **`fs-time`**
  (Lie-group rigid-body + symplectic integrators) — both probe-compiled to wasm32
  clean during the Round 0.5 audit. **`fs-lbm`** provides D2Q9/D3Q19 cores for the
  offline reference tier; **`fs-vpm`** exposes its desingularized Biot–Savart kernel
  functions for reuse; **`fs-viz`** provides scientific-visualization primitives.
- The **evidence doctrine** (Verified/Validated/Estimated, typed refusals, no-claims
  blocks, content-identity binding) is the honesty instrument that keeps "real
  physics, real time" an audited claim, in the UI as rigorously as in the kernels.

### 1.3 The one-sentence positioning

> A museum-quality interactive scientific instrument: the most historically and
> physically serious Wright Flyer simulation that runs in a browser tab — where the
> instability, the pilot, the canard mechanics, the wind, and the model's own
> uncertainty are all visible at the same time.

### 1.4 Explicit non-goals (v1)

- **Not a general flight simulator.** One aircraft family (Flyer I/II/III lineage),
  two sites, minutes-long flights.
- **Not a CFD product.** The browser tier is a reduced-order model with a resolved
  wake; full-field CFD (`fs-lbm`) runs offline in the reference program.
- **No multiplayer, no VR in v1** (architecture must not preclude a v2).
- **No gamification.** The physics is the game.
- **The browser tier never claims better than `Estimated` evidence** except where
  subsystem-specific validation receipts apply (§10.6); the UI displays evidence
  composition honestly, in color AND plain language (§2.5).

---

## 2. Product Experience

### 2.1 The five-minute journey (first-time user, Round-1 revised)

1. **Arrival.** Dawn at Kill Devil Hills: the camp, the 60-ft wooden launch rail
   pointing into a ~25 mph wind, the Flyer on its dolly, sand streaming in low
   ribbons. A date stamp: December 17, 1903, 10:35 a.m.
2. **Ride along.** The first experience is a clearly labeled **historical-pilot
   hypothesis replay**: the modeled pilot flies the authentic configuration in an
   ensemble-drawn December wind. The user watches the porpoising happen — from the
   wing, from the sand, from Daniels' tripod — with a caption explaining that this is
   a modeled hypothesis of the flight, not a recording.
3. **Fly with Training Assist.** The user takes the controls with the accessibility
   assist engaged and succeeds — imperfectly. The results card compares their run to
   the historical record and to the hypothesis replay, and reports the actual
   distance distribution context rather than a single "beat this" number.
4. **Authentic Controls.** Invited, not forced: "Now try it the way Orville did."
   Raw mechanical controls, raw instability. Most users crash within seconds — and
   the card explains *why this is the point*, with a link to…
5. **"Why it porpoises."** The flagship educational view (§2.4.8): open-loop
   divergence, canard hinge behavior, pilot delay and saturation, all animated
   against the flight they just flew. Then the design panel opens: change one
   physical parameter (CG, canard balance, headwind), keep the same atmosphere
   realization, fly the A/B comparison.

### 2.2 The advanced journeys

- **The engineer's loop.** Full config editor (typed, unit-checked, `fs-qty`
  semantics), KPI dashboard with fixed- and free-control static margins, augmented
  eigenmode view, applicability-domain status, parameter sweeps in workers with
  common-random-number ensembles, replay export with full identity envelope.
- **The educator's loop.** Curated lessons: ground effect (fly at 3 m vs 30 m),
  anhedral vs dihedral in gusts, fixed-control vs free-control stability (the canard
  balance slider makes this tangible), adverse yaw with the warp–rudder coupling
  disengaged, div/grad/curl on the wind field.
- **The historian's loop.** The four December 17 flights as ensemble presets with
  documented uncertainty; the 1904 Huffman catapult problem; Flyer III 1905 circuit
  flights. Every preset cites sources and labels reconstruction vs record.
- **The cinematic loop (offline, native).** Browser replays re-rendered by the
  native pipeline (`fs-render` + the existing EXR→ProRes quarantined mux adapter)
  into film-quality clips — same identity envelope, hero cameras.

### 2.3 Interaction model & control modes (Round-1 revised)

The real Flyer's controls are part of the physics, not an input abstraction: the
canard was operated through a lever and cable against its own hinge moment, and
Orville reported the elevator was balanced too near its center and tended to keep
moving once started — overcontrol by mechanism, not only by pilot lag.

Four distinct modes (each a distinct model identity; see §5.1.4):

| Mode | What it is | Label |
|---|---|---|
| **Authentic Controls** | the human user acts through the full historical mechanical control model (lever ratios, hinge moments, friction, stops). No synthetic human delay is added — the user supplies real perceptual/motor latency. | historical mechanics, your reflexes |
| **Historical Pilot** | an autonomous modeled pilot (delay + lead/lag + neuromuscular lag + saturation + deterministic remnant) used for ensemble replays and demonstrations | modeled hypothesis |
| **Training Assist** | user intent + low-authority accessibility controller through the same actuator path | hybrid assist, not historical |
| **Modern SAS** | full stability augmentation | explicitly ahistorical |

Input mappings (keyboard / mouse-drag "hip cradle" / gamepad) are device mappings
onto pilot force/intent, remappable, sampled and deterministically quantized at
simulation ticks (§6.2).

### 2.4 Visualization modes

1. **Vector glyph field** — GPU-instanced arrows on a user-positioned probe volume;
   live animated.
2. **Streamlines / pathlines / streaklines** — with the terms used *precisely*:
   streamlines integrate one frozen field snapshot; pathlines time-integrate markers
   through successive snapshots; streaklines require continuous release and are
   available only when the field-history ring is retained.
3. **Vorticity view** — three honestly distinguished sources: Tier B connected wake
   (direct rendering of model vorticity/circulation), Tier A prescribed wake
   descriptors (labeled schematic/model-derived), ambient turbulence vorticity
   (analytic derivative of the synthetic atmosphere). No shared "exact" label.
4. **Divergence verification overlay** — normalized residual
   ε_div = |∇·u| / (‖∇u‖ + ε), analytic-vs-finite-difference toggle; the model
   promises solenoidal ambient flow and the overlay shows the discretization error —
   a teachable self-check.
5. **Kinematic-gradient view** — ∇u, strain magnitude, Q-criterion/λ₂ (definitions
   in-app), ∇(½|u|²) labeled *kinematic speed gradient*. **No volume pressure** — a
   kinematically synthesized field does not determine p; a Bernoulli proxy appears
   only on explicitly irrotational components, labeled a proxy. Surface pressure
   appears only when a model actually supplies a Cp distribution; otherwise the UI
   shows strip normal loads.
6. **Force overlay** — per-strip lift/drag vectors, prop thrust, weight/CG marker,
   net force/moment gnomon (exact sim values).
7. **Scalar probes** — draggable point probes with strip-chart history; a virtual
   anemometer with declared reference height.
8. **"Why it porpoises"** *(flagship, Round-1)* — synchronized time plots and a
   live loop diagram: open-loop pole and time-to-double; canard command vs actual
   deflection; hinge moment and pilot force; pitch rate, flight-path angle, height;
   delay/phase indicators; saturation and reversal events; the loop component
   currently driving growth highlighted. A/B mode changes one mechanical or pilot
   parameter with the atmosphere realization held fixed.

Every mode works while flying, paused, or scrubbing a replay, and every field output
carries provenance and validity metadata (§5.5).

### 2.5 Evidence UX (Round-1)

Evidence colors stay, and are never the sole channel. Every badge carries: color +
icon; a one-sentence plain-language claim; subsystem breakdown; applicability-domain
status; a receipt link; and the limiting uncertainty. Example: *"Main-wing
attached-flow loads validated against full-scale data at this Re and α; canard
separated-flow loads and the historical gust realization are estimated."* A
provenance inspector opens from any parameter, force, field, or KPI.

---

## 3. Historical & Physical Grounding

> Grounding rule: values marked **[V]** are load-bearing and must be re-verified
> against sources in E1 before entering goldens; **[V?]** are believed-approximate
> and must be sourced or demoted to tunable-with-provenance. Round 1 added source
> anchors: NPS Wright Brothers materials, Library of Congress flight accounts,
> Culick et al. (Caltech) stability analyses, Deters/Broughton/Selig AIAA-2004-0211
> full-scale Flyer simulation, Smithsonian NASM object records.

### 3.1 The 1903 Wright Flyer — reference configuration

| Property | Value | Status |
|---|---|---|
| Wingspan | 40 ft 4 in = 12.29 m | [V] |
| Chord | 6 ft 6 in = 1.98 m | [V] |
| Wing area (both wings) | ~510 ft² = 47.4 m² | [V] |
| Aspect ratio (per wing) | ~6.2 | [V] derived |
| Gap between wings | ~6 ft (≈ chord) | [V?] |
| Camber | ~1/20 (as flown 1903) | [V] |
| Wing section/construction | digitized rib geometry + fabric-covered construction; relationship to the Wright wind-tunnel section identifiers established by E1 (NOT "single-surface", per Smithsonian fabric records) | [V?] |
| Anhedral (droop) | slight, deliberate (gust response) | [V] qualitative |
| Canard | biplane elevator: two planes, own span/chord/gap/stagger, hinge + balance axes, ~48 ft² total | [V?] full geometry from drawings, E1.5 |
| Canard arm | ~7 ft ahead of CG | [V?] |
| Twin rudder (aft) | ~21 ft², movable, warp-coupled | [V?] |
| Empty weight | ~605 lb = 274 kg | [V] |
| Gross weight (pilot ~145 lb) | ~750 lb = 340 kg | [V] |
| Engine | 4-cyl inline, ~201 in³, ~12 hp peak, ~180 lb w/ accessories | [V] |
| Engine controls | topology pending E1 verification — distinguish gas-flow lever, ignition/spark, run/stop, non-pilot-adjustable settings; do NOT assume a modern throttle | [V?] |
| Transmission | chain drive, ratio ~23:8; one chain crossed for counter-rotation | [V?] |
| Propellers | 2 × pusher, ~8.5 ft (2.6 m), ~330–350 rpm, carved spruce; digitized radial geometry from reconstructions | [V] |
| Propeller performance | NPS cites ~66% for the original design; modern reconstruction figures vary by convention and operating point — record convention with every number | [V?] E1.6 |
| Flight airspeed | ~30–34 mph (13–15 m/s) | [V] |
| Static pitch stability | **unstable** (Culick et al.); barely stabilizable by a trained pilot | [V] |
| Canard mechanics | balanced too near center; tended to continue moving once started (Orville) — overcontrol by mechanism | [V] qualitative, E1.5 quantifies |
| Pilot position | prone, hip cradle; ~"almost 1,000" glides of 1902 practice | [V] narrative ranges |

### 3.2 The four flights of December 17, 1903 (validation anchors)

| # | Pilot | Distance | Duration | Notes |
|---|---|---|---|---|
| 1 | Orville | 120 ft | 12 s | 10:35 a.m.; undulating; wind ~24–27 mph [V] |
| 2 | Wilbur | ~175 ft | ~12 s | [V?] |
| 3 | Orville | ~200 ft | ~15 s | [V?] |
| 4 | Wilbur | 852 ft | 59 s | initially undulating, then a long steadier segment, renewed pitching and ground strike; landing damage and the LATER gust destruction of the aircraft must be separately sourced [V?] |

The historical control traces and gust time-series are **unknown**. These flights
are therefore **probabilistic validation anchors**: ensembles over documented
uncertainty, evaluated against pre-registered joint predictive regions (§10.4) —
never a single "correct replay." Wind provenance matters: Orville's diary
distinguishes the Wrights' anemometer reading from the government reading —
`WindReference` records instrument, height, and averaging interval (§5.4).

### 3.3 Sites

**Kill Devil Hills, NC (1903).** Flat sand plain; 60-ft launch rail; dolly; skids
(no wheels). Big Kill Devil Hill (~90 ft) nearby; powered flights from the flat.
Dec 17 wind 24–27 mph [V]. Sand z₀ ≈ 10⁻³–10⁻² m with sparse grass; two camp sheds.

**Huffman Prairie, Dayton OH (1904–05).** ~84-acre pasture; light variable winds →
catapult derrick from Sept 1904 (~1,600 lb drop weights, ~16 ft tower [V?]);
thermals, tree lines; first circle Sept 20 1904; the great Oct 5 1905 Flyer III
flight — Smithsonian: 24.5 mi / 39.2 km; NPS: 24⅕ mi / 38:03 — keep both in
provenance, select one primary record for the preset [V?].

Terrain data: USGS 3DEP/SRTM + historical adjustment; ~1–2 m heightfield over
2×2 km, "representative, not surveyed" label. 2×2 km suffices even for the 1905
flight: 24 miles is accumulated *circuit length* around the small prairie, not
displacement (Round-1 Q5 answer); E1.3 validates the historical circuit stays
inside the tile; a low-detail visual horizon ring covers the distance views.
Sky/lighting: temperature, wind, sun position, cloud state, and artistic lighting
are separate provenance fields; unsourced cloud/lighting choices are labeled visual
reconstruction.

### 3.4 The signature behavior, decomposed (Round-1 P0 revision)

The Flyer's longitudinal behavior is THREE separately modeled, separately validated
mechanisms — the historical porpoising is a *closed-loop trajectory*, not the
definition of the open-loop mode:

1. **(a) Open-loop airframe instability.** The geometry and aerodynamics produce an
   unstable longitudinal pole set. Validation target: the published open-loop
   derivative/pole structure (Culick et al.), with time-to-double reported — not a
   predeclared "2–4 s growing oscillation."
2. **(b) Canard control mechanics.** Hinge moment, aerodynamic balance near center,
   linkage ratio, friction, inertia, stops, backdrivability — the mechanism-level
   overcontrol Orville described.
3. **(c) Closed-loop pilot–aircraft oscillation.** A delayed, saturated feedback
   loop (human or modeled pilot) around (a) via (b). A pure delayed pitch-rate
   damper is NOT assumed sufficient — in the elementary reduction, rate feedback
   changes damping but not unstable static stiffness; the conclusion must come from
   the full augmented pole set.

Other fidelity-floor phenomena (each mapped to validation cases in §10): ground
effect at h/b ≈ 0.1–0.3 (flat-plane image system); near-stall cruise (separation
margin + lagged separation state); biplane and canard–wing interference (coupled
multisurface solve, §5.2); low-J propulsion through the rail run from J = 0
(CT/CQ maps + rotor dynamics, §5.3); gusty boundary-layer wind (wall-compatible
solenoidal synthesis, §5.4); warp adverse yaw and the coupled warp–rudder control
(lateral topology enum, §5.1.1).

---

## 4. System Architecture

### 4.1 Three planes

```
┌────────────────────────────────────────────────────────────────────────┐
│ PRESENTATION PLANE (TypeScript + three.js, apps/wright-flyer/)         │
│ scene graph, terrain+sky, aircraft model, field visualization,         │
│ HUD/instruments, config UI, replay scrubber, WebAudio                  │
├────────────────────────────────────────────────────────────────────────┤
│ SIMULATION PLANE (Rust→wasm32, crates/fs-flyer-wasm, own workspace)    │
│ fs-flyer (aircraft assembly, 6-DOF+added-mass, controls, pilots,       │
│           launch, KPIs, field service)                                 │
│ fs-wing (coupled multisurface unsteady aero + wake)                    │
│ fs-airscrew (BEMT + rotor dynamics + engine + drivetrain)              │
│ fs-airfoil (L2 section models)   fs-atmo (wind + turbulence)           │
│ on: fs-mbd, fs-time, fs-vpm, fs-bem, fs-exec, fs-math det,            │
│     fs-rand philox, fs-blake3 identity, fs-qty units, fs-scenario      │
├────────────────────────────────────────────────────────────────────────┤
│ REFERENCE / REFEREE PLANE (native FrankenSim, offline)                 │
│ fs-lbm D3Q19 wind-over-terrain; high-res panel/VPM aero re-runs;       │
│ cross-fidelity discrepancy receipts (fs-vvreg / vv-scorecard);         │
│ fs-contact certified no-missed-contact replay verification;            │
│ cinematic re-render of replays; four-quadrant+wasm golden hashes       │
└────────────────────────────────────────────────────────────────────────┘
```

**Round-1 rename:** the offline plane is the **reference/referee plane**, not the
"truth plane." fs-bem/VLM/VPM/LBM are higher-fidelity *models*; experimental data
and analytic solutions establish which claims they may referee. Outside those
domains they supply cross-model discrepancy, not truth — this matters most near
stall, where both inviscid prediction and historical experimental coverage are
limited.

### 4.2 The multi-fidelity ladder

| Tier | Where | Model | Role |
|---|---|---|---|
| A | browser, always | **coupled low-order multisurface** lifting-line/Weissinger solve over both main wings, both canard planes, and vertical surfaces; unsteady section states with declared effect ownership; FlatPlaneExact ground images; BEMT propulsion with rotor dynamics; wall-compatible synthetic turbulence | the flyable sim; KPIs |
| B | browser, capable machines | Tier A + resolved hybrid wake (connected near-wake filaments, coarsened mid wake, multipole far wake) with reduced-rate feedback | wake realism, viz truth |
| C | native, offline | high-res unsteady VLM/panel (fs-bem+fs-fmm) + dense VPM; cross-checks A/B attached-flow, interference, and wake predictions inside its validated domain | referee |
| D | native, offline | fs-lbm D3Q19 wind-over-terrain (no aircraft in v1 — Round-1 Q8) | cross-checks atmospheric statistics under matched conditions |

Ladder invariants (Round-1 identity revision, §4.4): one `PhysicalScenarioId` flows
through all tiers; each tier/approximation carries a distinct `ModelId`; residual
**correction tables are disabled by default** — they are explicit identity-bound
modes, fitted on calibration pins only, evaluated on holdout pins, constrained for
symmetry/smoothness/conservation, and applied only inside their own applicability
domain (§10.6).

### 4.3 Process/thread topology & QoS (Round-1 revision)

```
Main thread:        three.js render loop; UI; input capture (tick-quantized)
Worker "sim":       wasm physics @ fixed 120 Hz (+wake 40 Hz Tier B)
Worker "field":     wasm field-sampling service (immutable FieldSourceSnapshots)
Worker pool (0–N):  parameter sweeps, replay re-runs
SharedArrayBuffer:  seqlock state ring, field buffers (when cross-origin isolated)
```

- **Global QoS governor:** sim deadline > audio continuity > render continuity >
  field freshness > sweeps. Field and sweep workers receive explicit
  pause/rate/cancel commands when the sim lag monitor crosses declared thresholds.
  **The physics tier never changes during a run**; if the sim cannot meet its
  fixed-step contract, it pauses with a typed performance refusal and offers a new
  run at a lower tier.
- **Fixed timestep** decoupled from render; render interpolates. All schedules
  (120 Hz sim, 240 Hz contact substep, 40 Hz wake, field refresh) are integer-ratio
  locked to the sim tick.
- **Fallback without cross-origin isolation:** each wasm instance is single-threaded
  but sim/field/sweeps still run in separate Web Workers with pools of transferable
  preallocated ArrayBuffers (no per-frame allocation/cloning). Tier availability is
  selected by *measured capability* (E0.6 benchmarks), not by SharedArrayBuffer
  presence alone; the compatibility banner reports the actual disabled features.
- **Tab suspension:** on visibility loss, execution pauses at a complete simulation
  tick; never unbounded wall-clock catch-up on resume; the pause is presentation
  metadata, not simulation time.

### 4.4 Identity model (Round-1 revision)

Five distinct identities, all `fs-blake3` domain-hashed, registered in
`identity-authorities.json` at introduction:

| Identity | Contents |
|---|---|
| `PhysicalScenarioId` | aircraft design, site, initial conditions, weather *distribution*, pilot hypothesis, launch system |
| `ModelId` | tier, every approximation/fast mode + parameters, correction-table selections, discretizations, timestep, solver modes |
| `ArtifactId` | exact wasm/native executable + embedded data artifacts |
| `InputTraceId` | tick-addressed, deterministically quantized control trace |
| `RunId` | hash of all of the above + realization seed |

"Same physical scenario, different model" is thus expressible without pretending two
runs share a complete identity; cross-fidelity receipts explicitly bind the pair.

### 4.5 Replay envelope

Every replay records: `PhysicalScenarioId`, `ModelId` (with complete mode
parameters), `ArtifactId`, `InputTraceId`, realization seed, terrain/table/
correction identities, schema version, optional checkpoint digests. Bit-identical
re-run on the same `ArtifactId`; numerically compatible cross-platform per the
golden program (§10.7). **Replay schema freeze (Round-1 Q7): freeze the v1 envelope
now**; migrations parse old envelopes, mint new identities, preserve originals, and
state which semantics changed — never silently reinterpret an old replay under new
physics. Old-exact playback uses the archived artifact by `ArtifactId`.

---

## 5. Physics Core Design

### 5.1 `fs-flyer` — aircraft assembly & flight dynamics (new crate, L4)

**Owns:** parametric airframe (geometry+mass+structure modes), mechanical control
system, pilot models, launch systems, ground interaction, per-step orchestration,
KPI accumulation, scenario schema, field service. Depends on `fs-wing`,
`fs-airscrew`, `fs-airfoil`, `fs-atmo`, `fs-mbd`, `fs-time`, `fs-qty`,
`fs-scenario`, `fs-blake3`, `fs-exec` (real-time), with `fs-contact` as an OFFLINE
dev-dependency of the reference-plane replay verifier only.

#### 5.1.1 Parametric airframe (Round-1 expanded)

```
FlyerDesign {
  // Lifting system
  span_m, chord_m, gap_over_chord, camber_ratio, stagger_m,
  anhedral_deg, washout_deg, wing_section_id,

  // Canard: two physical elevator planes, not one equivalent area
  canard_span_m, canard_chord_m, canard_gap_m, canard_stagger_m,
  canard_arm_m, canard_camber_ratio, canard_section_id,
  canard_hinge_axis_fraction, canard_balance_axis_fraction,
  canard_mass_kg, canard_pitch_inertia_kg_m2,
  canard_min_deflection_deg, canard_max_deflection_deg,

  // Rudder and lateral-control topology
  rudder_span_m, rudder_chord_m, rudder_arm_m,
  lateral_control_topology: LateralControlTopology,

  // Mechanical controls
  canard_lever_ratio, canard_cable_compliance, canard_hinge_friction,
  canard_control_force_limit,
  warp_lever_ratio, warp_cable_compliance,
  elevator_rate_limit, warp_rate_limit,      // hard safety caps (typed limits)

  // Mass & balance
  empty_mass_kg, pilot_mass_kg, ballast_kg, cg_offset_m,

  // Propulsion
  engine_power_hp, engine_rpm_max, chain_ratio,
  prop_geometry_id,                          // digitized radial geometry table
  prop_diameter_m, prop_pitch_m, prop_activity_factor,  // counterfactual levers

  // Structure
  warp_structure_mode: WarpStructureMode,
  spar_section_properties, wire_layout_id, wire_pretension,
  cable_compliance, structural_material_set_id,
}

LateralControlTopology {
  NoRearRudder1901,
  FixedTwinTailEarly1902,
  CoupledMovableRudderLate1902And1903,
  IndependentRudderLaterFlyer,
  Custom(LinkageMatrix),
}

WarpStructureMode {
  PrescribedKinematicEstimated,   // v1 default: kinematic warp, schematic deformation
  QuasiStaticBeamAndRigging,      // spanwise beam + wire tension/slack + loaded warp
}
```

Derived+displayed: main-wing and canard areas/aspect ratios; gap and stagger
ratios; wing loading; power loading; canard volume; control effectiveness;
hinge-moment gradient; **fixed-control AND free-control static margins**; trim
authority; CL required at target speed. When `QuasiStaticBeamAndRigging` is active:
spanwise bending moment, spar utilization, wire tension/slack status, loaded warp
distribution, prop/structure clearance. When it is inactive, structural margin and
structural optimization are **unavailable** rather than inferred from a load-factor
scalar (Round-1: the scalar n-limit cannot support structural claims across span/
chord/rigging changes; load-factor *exposure* remains a kinematic KPI, never
relabeled structural margin).

Typed refusals at admission for non-physical inputs with ranked repairs. Mass model:
component build-up → mass, CG, inertia tensor, cross-checked against published Flyer
inertia estimates [V].

#### 5.1.2 Six-DOF dynamics with generalized added mass (Round-1 P0 revision)

The aerodynamic assembly returns (1) acceleration-independent generalized loads,
(2) a symmetric generalized added-mass operator, (3) control-acceleration coupling
terms, (4) diagnostics. Each step solves

```
(M_rigid + M_added(q, δ)) · ν̇ = Q_nonaccel + Q_added_bias
```

with a deterministic `fs-la` factor/solve (6×6; larger when rotor/control
coordinates couple) BEFORE the Lie-group state update. **No finite-difference α̈ or
control acceleration is ever injected as an explicit force** — this kills the
implicit acceleration loop and its timestep-dependent noise.

Integrator: fixed dt = 1/120 s with a second-order predictor/corrector or midpoint
force evaluation wrapped around the `fs-time` Lie-group update
(`lie::rigid_body_step` / `quat_exp_step`); 1/240 s contact substep as an
integer-ratio schedule; the composition's order is validated explicitly (V-05a).
Type adapters between `fs-mbd`'s and `fs-geom`'s vector/quaternion types are owned
here (audited seam).

Parasite drag is a **component ledger** (pilot, engine/radiators, skids, wires,
uprights, struts, chains, misc.), each with area, orientation, Re-dependent
coefficient source, uncertainty, and power loss; the flat-plate aggregate remains a
separately identified fallback mode. At the Flyer's power margin, parasite-drag
error decides whether it flies at all (V-13).

#### 5.1.3 Mechanical control system (Round-1 P0 revision)

The canard control path is a one-DOF mechanical system: pilot force/torque → lever
and cable ratio → canard hinge dynamics, including canard inertia, hinge moment
(from `fs-airfoil` section data at the actual hinge/balance axes), aerodynamic
balance, friction, cable compliance, travel stops, and pilot-arm impedance. Rate
and force limits *emerge* from this model; hard safety caps remain as typed limits.
Wing-warp and rudder actuation use a separate linkage model (cable travel,
compliance, backlash/friction bounds) with the topology selected by
`LateralControlTopology`. This is what makes the canard-overbalance mechanism
(§3.4b) representable at all.

#### 5.1.4 Pilot models (Round-1 P0 revision)

- `PilotDirectHistoricalControls` — raw human input through the mechanical actuator
  model. No synthetic human delay.
- `PilotWrightModel` — autonomous historical-pilot hypothesis: regulates a
  combination of flight-path angle / near-ground height, pitch attitude, and pitch
  rate through transport delay, lead-lag/crossover dynamics, neuromuscular lag,
  force/travel saturation, and deterministic (philox-seeded) remnant. Capable of a
  genuine pilot–aircraft loop, not a fitted q-damper.
- `PilotTrainingAssist` — user intent + low-authority accessibility controller;
  labeled hybrid.
- `PilotSAS` — modern rate/attitude/flight-path feedback; labeled ahistorical.

**Identifiability rule:** no pilot parameter may be fitted to a validation case and
then credited with passing that same case. Calibration and holdout receipts carry
distinct identities (§10.4). Accessibility tuning is isolated from historical
calibration.

#### 5.1.5 Launch systems & ground interaction (Round-1 revision)

- **Rail:** dolly + aircraft as a unilateral constrained system; the constraint
  solver returns the rail-normal reaction; **release occurs when the admissible
  normal reaction reaches zero and the aircraft is separating**, with a small
  documented hysteresis — no empirical speed threshold. Rail geometry-error input is
  an optional *sourced* scenario field, disabled in historical validation unless the
  dossier supports it; it may not be used to tune porpoising.
- **Catapult (Huffman):** weight-drop energy → rope tension profile → dolly force
  (defaults [V?] E1).
- **Skid–sand contact:** fs-flyer-owned real-time model (heightfield penetration
  springs + regularized Coulomb friction + plastic "sink"). Landing produces an
  **impact report**: contact impulses, peak estimated skid/canard loads, penetration
  residual, structural-envelope status. Component breakage is disabled until a
  sourced reduced structural model exists; any cinematic breakage is labeled
  non-physical. `fs-contact` (certified spacetime CCD) runs OFFLINE as the
  replay-verification pass certifying no missed skid/terrain penetration (E10.5).
- **Terrain queries:** heightfield + material map sampled by contact and `fs-atmo`.

#### 5.1.6 Work ledgers (Round-1 revision)

Separated ledgers, not one "energy closure" scalar: (1) rigid-body mechanical
balance (force/moment power vs ΔKE+ΔPE); (2) engine/drivetrain/propeller power
balance (shaft power vs thrust power + induced + profile + chain loss); (3) contact
storage and dissipation; (4) resolved-wake consistency where the wake carries an
impulse/energy diagnostic; (5) numerical residual. "Aerodynamic work transferred
out of the rigid body" is not automatically called dissipation. Residual
normalization uses max(E_reference, E_floor) so engine-off/ballistic fixtures stay
meaningful.

### 5.2 `fs-wing` — coupled multisurface unsteady aerodynamics (new crate, L3)

Consumes `fs-airfoil` (§5.2.1) for all section data. Three layers:

**Planform layer (Round-1 P0 revision).** Tier A performs a **simultaneous coupled
bound-circulation solve** (lifting-line/Weissinger class, N ≤ ~80 unknowns) over
both main wings, both canard planes, and the vertical surfaces — no scalar biplane
factor in the production force path. Classical Munk/Prandtl biplane factors remain
verification fixtures and a separately identified emergency fallback mode. Tier B
raises spatial resolution and adds the resolved wake; it does not introduce
multisurface coupling for the first time. Surfaces enter section regimes
(attached/transitional/separated) independently. Includes warp-twist distribution
and adverse-yaw bookkeeping (V-12).

**Unsteady layer with effect ownership (Round-1 P0 revision).** An explicit
`AeroMemoryMode` selects the owner of each physical effect; incompatible
compositions refuse at scenario admission:

```
AeroMemoryMode:
  IndicialAttached      { motion_kernel: RationalApprox, gust_kernel: RationalApprox }
  ResolvedWake          { wake_model_id, optional_near_field_residual_id }
  IndicialWithSeparatedFlow { motion_kernel, gust_kernel, separation_lag }
```

| Physical effect | Permitted owner |
|---|---|
| Motion-induced circulatory lag | indicial states **or** resolved wake — never both |
| Gust penetration | Küssner/Sears states **or** chord-resolved gust sampling |
| Noncirculatory load | the generalized added-mass operator (§5.1.2) — only |
| 3-D induced velocity | coupled bound solve + explicit/prescribed wake |
| Separation hysteresis | explicit lagged separation state |
| Viscous/profile drag | section model |
| Far-wake closure | bounded multipole/analytic wake mode |

Reduced time integrates ds/dt = 2|U_local|/c (never 2Ut/c). States initialize to
trim steady state unless an impulsive start is requested. Rational kernel order is
selected offline as the smallest meeting a pinned time- and frequency-domain error
battery against exact Wagner/Theodorsen and Küssner/Sears over the Flyer's reduced-
frequency band (provisional bands §10.5) — Jones' 2-state fit is the candidate, not
the axiom. Near stall, a lagged-separation state with hysteresis replaces the
memoryless sigmoid (still `Estimated`, honestly so). A hostile twin enabling both
resolved-wake circulation and full Wagner memory must refuse.

**Ground-effect layer (Round-1 revision).**
- `FlatPlaneExact` (v1 historical default): every bound/wake element reflected
  across one fixed aerodynamic ground plane — exact within the flat-wall inviscid
  model class. Reflection uses the coordinate-free axial-vector rule
  **ω′ = det(R)·R·ω**; images retain stable identities for the life of their
  sources.
- `SmoothedTangentPlane` (additive `Estimated` mode): ONE continuously filtered
  global plane with hysteresis, slope limits, and a paired residual battery; its
  receipts report plane origin/normal/filter constant/max boundary residual; it
  never inherits the word "exact."
- `HeightfieldBoundary` (deferred): reduced-rate source/panel boundary mode.
No per-element or abruptly switching tangent planes, ever.

**Wake layer (Round-1 P1 revision — hybrid, not particle-only).** A 2,000-particle
cap would discard ~98% of a 59-second wake; instead:
1. **Near wake:** connected vortex rings/filaments shed from the spanwise/temporal
   circulation differences (topology-preserving; tip vortices are real filaments —
   also the natural render primitive).
2. **Mid wake:** deterministic age/distance-based conversion to coarsened connected
   filaments or vector particles, preserving circulation, first vorticity moment /
   hydrodynamic impulse, and declared symmetry invariants.
3. **Far wake:** deterministic cell multipoles or analytic trailing-vortex
   descriptors; bounded tree/cell evaluation for aircraft feedback.
4. **Pruning:** only under a computed induced-velocity remainder bound over the
   protected aircraft region — no element is discarded merely because a population
   cap was reached.
Flight-feedback and visualization evaluations may use different approximations, each
separately identified; a visually dense wake never silently feeds a more approximate
force model or vice versa. Kernels build on `fs-vpm` (3-D extension upstreamed,
§11.4).

#### 5.2.1 `fs-airfoil` — section models (new crate, L2; Round-1 P1 revision)

Generic section machinery, consumed by wing, canard, rudder, AND propeller (fixes
the L3→L3 layering violation): section geometry; analytic thin-section and
flat-plate reference modes; provenance-bound coefficient datasets (separate wing /
canard / rudder / propeller tables — the Wright 1901 tunnel tables are a
trend/convention validation source, NOT automatically a full-scale 1903 polar;
full-scale and reconstruction data dominate where available); **shape-constrained
residual splines** over (α, log Re, δ) on the analytic baseline (not a neural fit —
Round-1 Q2: interpretability, derivative control, small-data behavior, OOD refusal);
normal/axial/moment representation through deep stall; explicit
attached/transitional/separated/post-stall regimes; uncertainty intervals on every
query; applicability-domain refusal instead of unconstrained extrapolation; indicial
kernel definitions and their exact-reference batteries.

### 5.3 `fs-airscrew` — propulsion (new crate, L3; Round-1 P1 revision)

- **BEMT** with Prandtl tip/root loss over provenance-bound **digitized radial
  geometry** (chord, twist, section id per station — E1.6); high-level
  diameter/pitch/activity-factor remain counterfactual levers. The design-commit
  path constructs a bounded **CT/CQ map including J = 0**; runtime interpolates the
  map or performs a warm-started bounded station solve. Low-J momentum correction;
  typed refusal on non-convergence.
- **Rotor dynamics:** `I_eq·Ω̇ = Q_engine(Ω, controls) − Q_prop,L − Q_prop,R −
  Q_drivetrain`, with separately diagnosable left/right shaft torques and optional
  failure states — no algebraic RPM jumps during the rail transient.
- **Inflow:** local atmosphere + aircraft motion + coupled main-wing wake, with at
  least disk-average + first azimuthal harmonic distortion (pushers behind the wing;
  counter-rotation cancels only nominal torque, not unequal-inflow effects). Uniform
  inflow remains an identified fast mode.
- **Engine:** torque-vs-rpm curve is load-bearing; control topology per §3.1 pending
  E1; thermal derate is an optional `Estimated` mode DISABLED in historical presets
  until sourced (Round-1 correction of the invented 60-s derate).
- **Prop sections** use `fs-airfoil` interfaces with propeller-specific data,
  Reynolds, and roughness provenance — never silently reusing the wing polar.
- Validation: CT(J), CQ(J), static thrust/torque/RPM from J = 0 through the
  envelope; η(J) only where useful-power efficiency is well-conditioned (V-03).

### 5.4 `fs-atmo` — wind, turbulence, gusts (new crate, L2; Round-1 P0/P1 revision)

- **WindReference** records speed, reference height, averaging interval, direction,
  instrument/source identity, and uncertainty (the Wrights' vs the government
  anemometer is a modeled distinction, not trivia).
- **Mean profile:** log-law over terrain-relative height with a wind-direction-
  dependent **effective roughness from upstream fetch**; material transitions
  blended over a declared length (z₀ never jumps pointwise); neutral stability as
  an explicit declared assumption; other stability classes additive modes.
- **Turbulence:** solenoidal synthesis via a **wall-compatible vector potential**
  (or equivalent mirrored solenoidal modal basis): wall parity enforces zero normal
  velocity at the flat ground; vertical variation lives INSIDE the analytically
  differentiated potential (never a post-hoc envelope that would break ∇·u = 0);
  von Kármán-shaped anisotropic amplitudes; modal amplitudes/phases evolve by
  deterministic stationary OU-like updates projected onto the solenoidal subspace;
  counter-addressed by (scenario identity, mode index, simulation tick); a
  recurrence-time battery covers the ~39-minute 1905 scenario.
- **Gust events:** deterministic 1-cosine ramps for lessons/challenges.
- **Thermals (Huffman, v1.5):** seeded convective plumes; off for 1903.
- **API:** `sample(x, tick)` / `sample_batch(points, tick)` → velocity, analytic
  gradient, component provenance, applicability diagnostics. Optional fast mode
  `AtmosphereSampling::AffineLocalField { anchors, order }` — exact evaluation at a
  few aircraft reference points, local affine reconstruction across the airframe,
  with a paired exact-modal error battery over span-scale fixtures (§11.4).
- **Dec 17 preset:** an **uncertainty-conditioned ensemble** over documented mean-
  wind ranges, surface assumptions, and qualitative gust evidence. It never claims a
  recovered historical gust trace — no historical 10-m turbulence time series
  exists to match (Round-1 correction).
- Validation split: analytic construction (V-04a), statistical targets (V-04b),
  cross-model comparison with Tier D as a *cross-fidelity receipt* (V-04c).

### 5.5 Field-sampling service (module in `fs-flyer`; Round-1 revision)

```
sample_field(grid_spec, tick) -> {
  u[], grad_u[], omega[], div_analytic[], div_finite_difference[],
  strain_magnitude[], q_criterion[], kinetic_speed_gradient[],
  validity_mask[], provenance[],
  meta { source_tick, source_modes, force_coupled_components,
         visualization_only_components, core_radius, export_precision }
}
```

- **No `grad_p_hat`** (Round-1 P0): a kinematically synthesized solenoidal velocity
  field does not determine pressure. Kinematic quantities only; Bernoulli proxy only
  on explicitly irrotational components, labeled; surface Cp only from models that
  actually produce it.
- Consumes an immutable **FieldSourceSnapshot** (tick, physical state, bound
  circulation, wake buffer references, image-plane state, atmosphere seed/tick,
  kernel/core parameters, model identities); results carry the source tick and are
  never presented as current beyond their staleness threshold.
- Ambient parts analytic; wake parts via the shared (identified) induction mode;
  the exact-vs-finite-difference derivative dual powers the divergence overlay.
- Budgeted, worker-isolated, `fs-exec` cancellation-checkpointed, QoS-governed.

### 5.6 Sound (stretch, M6; Round-1 corrected)

AudioWorklet synthesis driven by state: 4-cylinder firing at engine RPM;
**per-propeller blade-passing at f_BPF = blade_count × prop_RPM / 60** — two
spatially separated, partially coherent 2-blade sources with relative phase, loading
harmonics, and unequal-inflow modulation (NOT "2×rpm×2 blades"); airspeed- and
separation-shaped wind/airframe noise; rail and contact events. Labeled "sound
design informed by physics, not an acoustic claim"; an fs-phs/fs-aeroac-grounded
path is v1.5+ (Round-1 Q10).

---

## 6. Numerics, Determinism, and the Execution Contract

1. **Determinism doctrine.** All transcendentals via `fs-math det::`; all randomness
   via `fs-rand` philox keyed by scenario identity + mode index + tick; no
   wall-clock, no map-iteration order, no thread-order dependence; deterministic
   cell ordering in wake binning; fixed-shape reduction trees.
2. **Tick-addressed inputs.** Input traces are keyed by integer simulation tick and
   deterministically quantized after device sampling, before actuator dynamics.
3. **Fixed integer-ratio schedules** (§4.3) make a replay a pure function of its
   identity envelope (§4.5).
4. **Structured determinism checkpoints (Round-1, moved early: E3.5).** Per-tick
   subsystem digests (atmosphere, section loads, circulation, propulsion,
   generalized loads, integrator state) so divergences localize to a subsystem
   instead of a whole-run hash mismatch.
5. **`fs-exec` everywhere**; ExecMode::Deterministic in shipped paths; typed
   refusals across the wasm boundary with ranked repairs; no silent clamping.
6. **Units:** `fs-qty` at module boundaries; SI doubles with documented units at the
   wasm boundary.
7. **Work ledgers** per §5.1.6 with per-scenario pinned residual envelopes.
8. **Identity:** the §4.4 quintuple; identity constants registered at introduction.

---

## 7. WASM Engineering & the Real-Time Budget

### 7.1 Build & packaging

New crate `crates/fs-flyer-wasm` with its own `[workspace]` on the `fs-wasm` pattern
(nested lock + lock-drift CI gate + asupersync canonical wasm profile feature
pinned — audit §11.3). `wasm-pack --target web`, `wasm-opt -O3`; SIMD and no-SIMD
artifacts; threads build as the enhanced artifact; single-threaded baseline always
shipped; loader feature-detects. Bundle target < 8 MB gz. wasm memory is reserved to
its run-time maximum before typed-array views are exported; growth during a run is
refused.

### 7.2 The frame budget (Round-1 corrected arithmetic)

Tier A per 120 Hz step (planning estimates, mid-range laptop, scalar wasm):

| Module | Per step | Per 60 fps frame (2 steps) |
|---|---|---|
| fs-atmo sampling (~80 pts) | 0.10 ms | 0.20 |
| fs-wing coupled solve + strips (~80 unknowns, 36 strips) | 0.30 ms | 0.60 |
| fs-airscrew (2 props, map interpolation + rotor ODE) | 0.06 ms | 0.12 |
| ground contact (fs-flyer) + rail | 0.03 ms | 0.06 |
| added-mass solve + fs-mbd/fs-time integrate | 0.04 ms | 0.08 |
| KPIs, ring buffer, bookkeeping | 0.05 ms | 0.10 |
| **Sim worker Tier A arithmetic total** | **~0.58 ms** | **~1.16 ms** |

Tier B arithmetic (wake at 40 Hz): wake advance 2.5 ms × 40/60 = **1.67 ms/frame**
average; strip feedback 0.6 ms × 40/60 = **0.40 ms/frame**; Tier B sim average
**≈ 3.2 ms/frame**; a frame containing one wake update ≈ **4.2 ms** on the sim
worker. Field service must state its interaction model: dense 32³×2,000 =
65.5 M interactions/refresh vs k-neighbor 32³×k + declared far-field — no field
budget is accepted without an exact interaction count and measured p50/p95/p99.

**These are planning estimates, not acceptance evidence.** FLOP counts are not a
performance proxy for rsqrt/regularization/memory/branching across browsers, and
deterministic transcendentals may make the atmosphere estimate optimistic — hence
E0.6 microbenchmarks precede physics tuning.

### 7.2.1 Performance acceptance contract (Round-1)

Per supported device/browser class: 120 Hz sim hard deadline 8.33 ms/tick; shipping
target sim p99 ≤ 6.0 ms with zero unbounded backlog; render p95 ≤ 16.67 ms at
selected presentation quality; input-to-visible-state latency distribution reported;
field refresh runtime and field-age distribution reported; no sweep may increase sim
deadline misses; the startup benchmark selects presentation/field quality and NEVER
silently changes the physics tier of an in-progress run. Benchmark suite records
scalar/SIMD, isolated/contended, cold/warm, SAB/transferable variants.

### 7.3 Interop contract (Round-1 revision)

- **Versioned seqlock state ring:** header { abi_version, model_id, tick,
  published_slot, sequence }; writer marks sequence odd, writes its owned slot,
  publishes even; readers retry on torn/odd sequence. Checksummed snapshots.
- **FieldSourceSnapshot** double buffer per §5.5.
- Field buffers zero-copy to render when SAB is available; GPU upload remains an
  explicit measured transfer. Non-SAB fallback uses transferable buffer pools.
- Command channel: postMessage JSON for config changes only; every design commit
  re-derives tables off-thread with progress + cancellation and mints new
  identities. No per-frame JSON.

### 7.4 Startup self-test

In-browser: 1-second canonical scenario headless → compare against the embedded
golden for this ArtifactId; report tier availability, chosen tier, and the golden
verdict on the About panel; visible "determinism self-test failed" badge on
mismatch. Never silently ship wrong physics.

---

## 8. Rendering & UX (three.js)

### 8.1 App shell

`apps/wright-flyer/` — TypeScript + Vite + three.js. **Renderer decision
(Round-1 Q9):** use three.js's renderer abstraction with `WebGPURenderer`/TSL where
practical and a tested `forceWebGL` WebGL2 fallback; WebGPU compute may accelerate
field visualization but NO v1 physics or validation feature may require it; field
semantics compared across backends. Static hosting + COOP/COEP; degraded mode works
without.

### 8.2 The aircraft asset

Source vetting (E2.1): Smithsonian 3-D digitization of the 1903 Flyer (license
verified in-task), NASA releases, or commissioned model; provenance in-app.
Pipeline: glTF 2.0, LODs, KTX2, rig (warp morph driven by warp state — under
`QuasiStaticBeamAndRigging` the LOADED warp distribution drives the morph; canard
planes pivot about the modeled hinge axis; props; chain scroll; prone pilot with
hip-cradle pose). Physics never reads the visual mesh; a calibration overlay draws
the physics surfaces/strips over the visual model (debug + honesty view). Procedural
morphs beyond ±25% design change are flagged "schematic preview."

### 8.3 Terrain, sky, environment

Heightfield terrain (2×2 km, 1–2 m res) + material splat; 1903 camp layout: rail,
sheds, Daniels' tripod (instant-photo mode with period grain), Big Kill Devil Hill;
Huffman: pasture, fences, derrick, tree lines, low-detail horizon ring. Sky presets
with separate provenance fields for temperature/wind/sun/cloud/artistic lighting
(§3.3). Ambient wind cues driven by real fs-atmo samples: sand streamers, fabric
flutter, flag, camp smoke.

### 8.4 Field visualization implementation

Instanced glyphs (≤ ~30k) from field buffers; GPU streamline integration through the
latest immutable snapshot; CPU hero-rake pathlines from the field-history ring
(disabled when field age violates declared limits); wake filaments rendered directly
with age-fade (they carry Γ); divergence overlay per §2.4.4; separate composited
overlay pass. Legends and transfer functions always visible.

### 8.5 HUD & instruments

Period triad (anemometer, stopwatch, engine-revolution counter — the Wrights' own
KPI set) + modern overlay (airspeed, altitude, α, load-factor exposure, pitch rate,
thrust, L/D, ledger strip). Engineer mode: augmented eigenmode view distinguishing
airframe / actuator / aero-memory / rotor / pilot-loop modes — not everything
labeled "short-period/phugoid"; the four-state longitudinal projection appears only
as a labeled teaching view. Results card per §9; design-diff card decomposes KPI
changes into induced/profile/parasite/propulsion/trim/structural contributions.

### 8.6 Replay UI

Timeline with event ticks (liftoff by rail-reaction criterion, gusts, separation
flags, reversals, touchdown), camera presets (chase, wingtip, Daniels tripod,
onboard prone view, free), A/B ghost mode with shared atmosphere realization,
export of the full identity envelope.

---

## 9. Configuration Space, Experiments & KPIs

### 9.1 Scenario schema & the three design spaces (Round-1 revision)

```
FlyerScenario = {
  design_family: HistoricalReconstruction | WrightCounterfactual | FreeformTeaching,
  design: FlyerDesign,
  site, weather_distribution, weather_realization,
  launch, pilot_hypothesis, model_selection,
}
```

- **HistoricalReconstruction:** dossier-supported parameter intervals only; geometry
  relationships locked.
- **WrightCounterfactual:** plausible modifications preserving the construction
  system, dependency rules (span change → mass/inertia/strut/wire/clearance
  re-derivation), and structural applicability checks.
- **FreeformTeaching:** broader abstractions, prominently labeled schematic.

Design edits are **transactional**: a running sim's `PhysicalScenarioId` is
immutable; commits create a candidate scenario + derivation job (cancellable,
progress-reporting). Schema versioned and frozen per workspace doctrine; presets
for the four flights, Huffman 1904, Flyer III 1905, challenges.

### 9.2 KPI definitions (Round-1 revision)

| KPI | Definition |
|---|---|
| Historical downrange distance | projection of touchdown displacement onto the launch/wind axis (the 120 ft / 852 ft comparator) |
| Horizontal path length | ∫ horizontal groundspeed dt |
| 3-D trajectory length | ∫ total groundspeed dt |
| Endpoint displacement | liftoff→touchdown Euclidean |
| Air distance | ∫ airspeed dt |
| Duration | liftoff→touchdown |
| Max/mean airspeed | airborne phase |
| Section separation margin | min provenance-bound margin to the local section separation boundary |
| Stall/separation exposure | time + span fraction in transitional/separated regimes |
| Raw ride metrics | RMS vertical specific force, RMS pitch rate, jerk, PIO-band energy, peak-to-peak pitch, reversal rate |
| Smoothness index | optional dimensionless composite from published reference scales; raw components always shown |
| Control activity | travel, rate, reversal count, saturation time per control |
| Estimated pilot work | ∫ pilot force × control velocity (mechanical-control model only) |
| Fixed/free-control static margin | only at a valid declared trim; includes canard mechanical state |
| Mode set | augmented linearization eigenvalues (§9.3) |
| Load-factor exposure | kinematic; never relabeled structural margin |
| Structural margin | only under `QuasiStaticBeamAndRigging`; reports limiting component + uncertainty |
| Ledger residuals | per §5.1.6 |

### 9.3 Experiments

- **Live linearization (Round-1 revision):** two products — the FULL augmented
  linearization over rigid-body + aero-memory + actuator + propulsion (+ pilot when
  displaying closed-loop modes) states, which owns stability claims; and a labeled
  reduced longitudinal projection for teaching. Available only when trim
  continuation succeeds and every subsystem is inside its differentiable
  applicability domain; otherwise a typed refusal naming the limiting subsystem and
  nearest valid operating point.
- **Sweeps:** worker-pool batches over config grids with common-random-number
  ensembles (same realization seeds across design points); progress-streamed;
  cacheable by (ScenarioId, ModelId, seed).
- **Optimization (v1.5, Round-1 gated):** robust multiobjective exploration over
  uncertainty via fs-bo/fs-dfo — requires active structural model,
  applicability-domain enforcement, correction-model holdouts, and CRN ensembles.
  Disabled under `PrescribedKinematicEstimated`.

---

## 10. Validation & Evidence Program

### 10.1 Anchor datasets — partitioned before calibration (Round-1 P0)

| ID | Dataset | Role |
|---|---|---|
| A1 | Wright 1901 wind-tunnel tables | section trends/convention validation |
| A2 | Modern re-tests: full-scale 1903 replica tunnel campaigns, AIAA Wright Flyer Project, Deters/Selig 2004, prop reconstruction data | section+planform forces, prop CT/CQ/static |
| A3 | Dec 17 flight records (4 flights, wind ranges, diaries) | probabilistic end-to-end anchors |
| A4 | Culick et al. stability analyses | open-loop derivative/pole reference |
| A5 | Classical biplane interference theory | verification fixtures |
| A6 | Boundary-layer meteorology relations + Tier D runs | atmosphere statistics |

Each dataset is partitioned BEFORE model calibration into: calibration subset;
held-out validation subset; convention/digitization uncertainty; applicability
domain; quantities for which it is independent evidence. Tier C/D outputs are
cross-model references until independently validated for the exact quantity
compared.

### 10.2 Validation hierarchy (Round-1)

1. Mathematical verification (exact/converged references).
2. Calibration (explicitly identified subsets only).
3. Held-out component validation (withheld α, Re, deflections, geometries, prop
   operating points).
4. Cross-fidelity discrepancy (Tier A/B vs independently validated Tier C/D).
5. Historical posterior-predictive checks (ensembles over documented uncertainty).
6. Performance qualification (device/browser matrix).

Verification batteries follow workspace law: refusals tested at cap AND cap+1; no
vacuous limit checks; falsifier-style negatives per gate; per-strip oracles, not
only totals; the effect-ownership hostile twin (§5.2) must refuse.

### 10.3 The V-cases (Round-1 rebuilt)

| ID | Case | Gate |
|---|---|---|
| V-01 | Section & full-aircraft steady-load holdouts | coefficient/derivative/uncertainty-calibration metrics on data not used for fitting |
| V-02a | Open-loop longitudinal derivatives & poles | full pole/derivative set vs A4 within declared reconstruction uncertainty; time-to-double reported |
| V-02b | Canard control mechanics | hinge-moment sign, self-driving tendency near the historical balance point, control-force/travel response, stop behavior within the sourced mechanical envelope |
| V-02c | Closed-loop pilot–aircraft | held-out PilotWrightModel parameters stabilize only within a finite delay/gain region; oscillation frequency, phase, saturation, and reversal statistics match sources |
| V-03 | Propulsion maps | CT/CQ/static thrust/torque/RPM from J=0 through the envelope within experimental or predeclared discrepancy bands; η only where well-conditioned |
| V-04a/b/c | Atmosphere | analytic construction (div, wall, determinism, stationarity, recurrence) / statistical targets (PSD, TI, Reynolds-stress ratios, integral scales, coherence, gust quantiles) / Tier D cross-fidelity receipt |
| V-05a/b/c | Convergence & discrepancy | 120→240→480 Hz trajectory/KPI/contact convergence / wake 20→40→80 Hz + field-rate sensitivity / Tier A-vs-B differences REPORTED with uncertainty, not required to vanish |
| V-06a/b/c | Ground effect | flat-wall exactness verification / all six load components + derivatives vs h/b, pitch, roll against high-res references / smoothed-tangent residual envelope in its declared slope domain |
| V-07 | Four-flight historical ensemble | all four observed distance-duration pairs inside pre-registered joint predictive regions; sharpness reported, not only coverage |
| V-08 | Unsteady indicial responses | time+frequency exact-reference battery per §10.5 bands |
| V-09 | Coupled biplane/canard loads | gap, stagger, α, control, h/b holdouts vs referee |
| V-10 | Wake invariants & induction | circulation, impulse, topology, induced-field error, pruning bound |
| V-11 | Rail launch | acceleration, release point (unilateral criterion), RPM transient, reaction-force history |
| V-12 | Lateral control | adverse-yaw sign/magnitude, coupled effectiveness, roll mode, turn coordination |
| V-13 | Parasite drag & power balance | component ledger vs full-aircraft drag/required-power envelope |
| V-14 | Browser real-time contract | §7.2.1 p50/p95/p99 incl. contention |

### 10.4 Historical pass logic (Round-1 P0)

Historical cases are probabilistic because control and gust traces are unknown.
A pass requires: pre-registration of every uncertainty distribution (E10.0 freeze);
observed joint outcomes inside declared predictive regions; sharpness reported;
no post-hoc distribution widening under the same validation identity; component
validation staying green. No pilot/atmosphere parameter fitted on a case may be
credited on that case (V-02c holdout rule). Nothing is gated on "reachable."

### 10.5 Provisional numerical bands (ratify from dossier in E10.0)

Indicial approximations: max normalized step-response error ≤ 2%, gain ≤ 2%, phase
≤ 3° over the declared reduced-frequency band. Airborne timestep convergence:
principal KPIs ≤ 1%, trajectory RMS ≤ 0.5%, contact timing/impulse ≤ 2%.
Atmosphere: median PSD ratio within ±3 dB over predeclared energy-containing bins,
no bin beyond ×2; TI and integral scales ±20%; Reynolds-stress ratios ±25%; gust
quantiles ±20%. Wake fast modes: induced-velocity RMS ≤ 3%, max ≤ 10% outside
cores; circulation residual ≤ 0.5%; impulse residual ≤ 1%. Real time: sim p99 ≤
6 ms against the 8.33 ms deadline. Historical endpoints: pre-registered joint
regions from the dossier uncertainty budget — no universal multiplicative band.

### 10.6 Applicability domains replace the "validated envelope" (Round-1 P0)

Validation is claim- and subsystem-specific. Each subsystem publishes an
`ApplicabilityDomain` in physical dimensionless coordinates (section: Re, α, δ, k,
roughness, separation state; biplane/canard: gap/chord, stagger/chord, CL/CL,max,
canard volume; ground effect: h/b, pitch, roll, slope; prop: J, radial Re, inflow
distortion; dynamics: static margin, q̄, control free/fixed; atmosphere: z/z₀, TI,
L/z, stability class; wake: reduced frequency, loading, coarsening error). The
run's claim domain is the intersection; the UI reports which outputs are validated,
which estimated, the LIMITING subsystem and coordinate, and interpolation vs
extrapolation status. No normalized-slider Euclidean distance promotes a run.

### 10.7 Determinism goldens & E2E runner

Four-quadrant golden replay hashes extended with a wasm column ({aarch64, x86,
wasm-in-node} × {debug, release}); golden-bump protocol applies. Subsystem digests
(§6.4) localize divergence. `scripts/ci/e2e_wright_flyer.sh` cloned from the
hardened Euler cinematic runner: `--list/--check/--self-test/--run smoke/
--negative CASE/--replay`, bounded JSONL logging contract, hostile twins (config
vs identity tamper, input-trace truncation, seed mismatch, ledger violation
injection, stale correction-table identity, KPI-vs-recompute mismatch, wasm/native
golden divergence, terrain-hash drift, effect-ownership double-count, post-hoc
distribution widening). Runner reuses production CLIs; never parallel logic.

---

## 11. Crate Reuse Matrix & New Crates

### 11.1 Existing crates leveraged (verified present)

| Crate | Role here | Notes |
|---|---|---|
| fs-wasm | pattern + infra precedent (workspace, CI lane, wasm-pack recipe) | fs-flyer-wasm copies its protections |
| fs-bem (+fs-fmm) | referee-plane force cross-checks; bounded one-shot interference/residual data | wasm32-proven; `panel3d` confirmed; needs screening preset (§11.4) |
| fs-vpm | Biot–Savart kernels; 3-D + hybrid wake extension upstreamed | kernel exposure verified; O(N²) today |
| fs-lbm | Tier D wind-over-terrain | D3Q19 + boundaries exist; terrain-only in v1 |
| fs-mbd | 6-DOF rigid body, canonical quaternions, force/impulse-at-point API | leaf crate; wasm32 PROVEN by audit probe 2026-08-16 |
| fs-contact | OFFLINE certified no-missed-contact replay verification | certified spacetime CCD — never per-step forces |
| fs-time | Lie-group rigid-body + symplectic integrators | wasm32 PROVEN by audit probe 2026-08-16 |
| fs-exec | Cx, budgets, cancellation, deterministic mode | wasm-proven via fs-wasm |
| fs-math | det:: transcendentals | determinism backbone |
| fs-rand | philox streams | |
| fs-blake3 | identity domains (§4.4 quintuple) | register identities at introduction |
| fs-qty | unit-checked quantities at seams | |
| fs-scenario | BC/load-case algebra with dimensions + provenance (verified general) | wasm32-proven |
| fs-simd | SIMD tiers | NO wasm tier today — E0.5 adds SIMD128 capsule; scalar fits budget |
| fs-la | dense factor/GEMM for coupled solve + added-mass solve | reuse-factorization API confirmed-or-added (§11.4) |
| fs-viz | field-viz primitives w/ analytic ground truth | |
| fs-uq / fs-surrogate / fs-bo / fs-dfo | ensembles, sweeps, optimization (gated) | in fs-wasm |
| fs-vvreg (+ vv-scorecard) | validation receipts + reporting | standing infra |
| fs-render + euler mux adapter + cinematic runner | offline cinematic export | just hardened (h7xu5) |
| fs-evidence | evidence colors / no-claims plumbing | |
| fs-ornith | prior art for staged aircraft campaigns | pattern reuse |

### 11.2 New crates (six — Round-1 adds fs-airfoil for layering correctness)

| Crate | Layer | One-line contract |
|---|---|---|
| fs-airfoil | L2 | generic section geometry, analytic baselines, provenance-bound coefficient tables, shape-constrained residual fits, indicial kernels, uncertainty, applicability-domain refusals |
| fs-atmo | L2 | wall-compatible solenoidal boundary-layer wind + gusts; analytic derivatives; ensemble presets; NO acoustic/thermo claims |
| fs-wing | L3 | coupled multisurface unsteady lifting-surface aero with effect-ownership modes, flat-plane images, hybrid wake |
| fs-airscrew | L3 | BEMT + rotor dynamics + engine + drivetrain on fs-airfoil interfaces; CT/CQ maps incl. J=0; static-regime honesty |
| fs-flyer | L4 | aircraft assembly: airframe+structure modes, mechanical controls, pilots, launch, contact, KPIs, ledgers, scenario schema, field service |
| fs-flyer-wasm | L6 | own-workspace wasm binding: sim loop, field service, sweeps, replay, seqlock ABI, typed-refusal JS API |

Each ships CONTRACT.md, no-claims block, refusal vocabulary, and registered
identity constants from day one.

### 11.3 wasm/real-time readiness audit (EXECUTED 2026-08-16, round 0.5)

| Crate | wasm32 evidence | Verdict |
|---|---|---|
| fs-exec, fs-alloc, fs-rand, fs-math, fs-qty, fs-la, fs-viz, fs-uq, fs-lbm, fs-bem, fs-vpm, fs-scenario, fs-render | in `fs-wasm`'s shipping dependency list (CI-built with wasm-pack) | proven by standing CI |
| fs-mbd, fs-time | probe crate compiled `RigidBodyState` + `lie::rigid_body_step` to wasm32-unknown-unknown clean (asupersync `wasm-browser-prod` feature) | proven by audit |
| fs-simd | NEON/AVX2 only | extension required (E0.5); scalar fallback budgeted |
| fs-contact | not probed | offline-only role |
| fs-fmm | not probed | offline referee only |

API shapes confirmed: fs-vpm exposes `induced_velocity`/`advect`; fs-la has dense
factor/GEMM; fs-time has Lie-group steps; fs-mbd exposes force/impulse application;
fs-bem has `panel3d`. Standing caveat: wasm builds of the fs-exec cone require the
asupersync canonical wasm profile feature — pinned by fs-flyer-wasm's workspace
(E0.3 DONE-WHEN).

### 11.4 Optional-fidelity doctrine for shared crates

1. **Additive, never mutative** — approximations are new flagged modes; exact-path
   defaults, semantics, and goldens do not move.
2. **Self-describing lossiness** — each mode names its approximation and appears in
   receipts.
3. **Paired error battery** — each fast mode ships a deviation-bound battery
   against the exact path on pinned fixtures.
4. **Identity discipline** — mode + parameters enter `ModelId`; a Tier-A and
   Tier-C run may share a `PhysicalScenarioId` but never a `ModelId`/`RunId`.
5. **Effect ownership (Round-1)** — every aerodynamic mode declares the physical
   effects it owns; two enabled modes may not own the same effect unless one is an
   explicitly derived residual with a paired no-double-counting battery;
   incompatible combinations are refused at admission.

Planned extensions:

| Crate | New optional mode | Consumer |
|---|---|---|
| fs-vpm | exact connected 3-D filament/ring path; vector-particle and deterministic tree/multipole paths; circulation-and-impulse-preserving coarsening; bounded-remainder pruning; `BinnedTruncated` fast induction (measured mode) | fs-wing hybrid wake |
| fs-bem | coarse-panel screening preset (declared panel budget, fs-exec bounded) + one-shot influence/residual export | interference data (E4.8) |
| fs-simd | wasm32 SIMD128 Tier-1w capsule | inner loops |
| fs-la | reuse-factorization API for repeated same-structure dense solves (if absent) | coupled solve, added-mass solve |
| fs-atmo | `AffineLocalField` sampling | airframe force path |

---

## 12. Milestones & Dependency-Aware Task Graph (Round-1 rebuilt)

### E0 — Program setup & performance ground truth
- **E0.1** Program root bead; plan→beads conversion (after review rounds).
- **E0.2** `apps/wright-flyer` scaffold (Vite+TS+three.js, COOP/COEP dev server, CI).
- **E0.3** `fs-flyer-wasm` scaffold on the fs-wasm pattern; asupersync wasm profile
  feature pinned. → blocks all wasm integration.
- **E0.4** wasm32 CI guard lane over the flyer cone.
- **E0.5** fs-simd SIMD128 Tier-1w capsule (opportunistic; after E4.7 profiling).
- **E0.6** Browser performance microbench suite: det transcendental batches,
  40–100-unknown dense solves, BEMT loops, exact+fast Biot–Savart kernels, bin/tree
  traversal, SAB publication, transferable fallback, Float32 GPU uploads.
  DONE-WHEN: p50/p95/p99 across the device/browser matrix. → informs E4.2/E4.7.
- **E0.7** Worker transport & suspension prototype: seqlock ring,
  FieldSourceSnapshot double buffer, transferable pools, visibility pause,
  no-catch-up, QoS throttling against a synthetic 120 Hz load. → E5.0.

### E1 — Historical grounding & data
- **E1.1** Source dossier A1–A6 with licenses and citations.
- **E1.2** Verify §3 [V]/[V?] values → `flyer-reference.json` (identity-hashed).
- **E1.3** Terrain: DEM both sites, 1903/1904 adjustments, heightfield + materials;
  1905-circuit containment check.
- **E1.4** Coordinate/convention/uncertainty dossier: body axes, moment signs,
  reference areas, control signs, wind reference heights, digitization uncertainty,
  calibration/holdout partitions, historical-input uncertainty distributions.
  → blocks E4 validation claims.
- **E1.5** Canard & control-mechanics dossier: geometry, hinge/balance axes,
  linkage ratios, travel/stops, mass/inertia, sourced qualitative behavior.
  → blocks historical-control claims.
- **E1.6** Propeller radial-geometry + operating-data package (CT/CQ/static
  fixtures with conventions).

### E2 — Assets & rendering foundation
- **E2.1** Flyer 3-D model vetting/acquisition (license task; blocking).
- **E2.2** Asset pipeline + rig. Depends E2.1 + E3.1's stable visual/physics
  geometry contract.
- **E2.3** Kill Devil Hills scene (the §2.1 arrival shot). Depends E1.3.
- **E2.4** Camera system, input mapping (tick-quantized), HUD skeleton.

### E3 — Simulation spine
- **E3.1** fs-flyer design schema + admission + mass/inertia build-up + derived
  panel. Depends E1.2/E1.5 (defaults may stub pending dossier).
- **E3.2** 6-DOF core on fs-mbd/fs-time; ring buffer; replay record/playback
  bit-identity. Includes the type-adapter seam.
- **E3.2a** Generalized added-mass assembly + second-order force/integrator
  coupling. DONE-WHEN: acceleration-dependent fixtures converge without FD
  acceleration noise; the effective-mass solve stays admissible over the reference
  design domain. → blocks E4.6a.
- **E3.3** fs-atmo v0: wall-compatible potential, mean profile + fetch roughness,
  seeds; batteries V-04a. (Parallel.)
- **E3.4** Rail (unilateral release) + fs-flyer contact + terrain queries.
  DONE-WHEN: dolly acceleration and release location converge under timestep
  refinement; no tensile rail reaction; landing impulse/penetration/friction work
  converge in declared bands.
- **E3.5** Structured determinism checkpoints (per-subsystem tick digests) — early,
  before physics churn.

### E4 — Aerodynamics & propulsion
- **E4.0** fs-airfoil crate (analytic modes, tables, constrained residuals,
  uncertainty, regimes, refusals). Depends E1.1/E1.4.
- **E4.1** Wing/canard/rudder/prop section datasets + models on fs-airfoil.
  Depends E4.0, E1.5, E1.6.
- **E4.2** Coupled Tier-A multisurface circulation solve (classical factors =
  fixtures/fallback). Depends E4.1, E0.6.
- **E4.3** Unsteady effect-ownership modes: variable reduced time, exact-reference
  rational kernels, resolved-wake exclusivity, separation lag. Depends E4.2, E3.3.
- **E4.4a** FlatPlaneExact image system (axial-vector rule; V-06a batteries).
- **E4.4b** SmoothedTangentPlane optional mode (only after E4.4a green).
- **E4.5** fs-airscrew: BEMT + CT/CQ map (J=0 up) + rotor dynamics + nonuniform
  inflow + engine + drivetrain + component power ledger. Depends E1.6, E4.0.
- **E4.6a** Open-loop integrated aircraft. Depends E3.2a, E4.2, E4.3, E4.4a, E4.5.
  DONE-WHEN: V-02a open-loop derivative/pole gates pass BEFORE any pilot exists.
- **E4.6b** Canard/warp/rudder mechanical controls. Depends E1.5, E4.6a.
- **E4.6c** PilotWrightModel + TrainingAssist (calibration subset only; assist
  tuning isolated). Depends E4.6b.
- **E4.7** fs-vpm hybrid wake: connected near wake, invariant-preserving
  coarsening, multipole far field, bounded pruning, exact-vs-fast batteries.
  DONE-WHEN: topology/circulation batteries green; coarsening invariants hold;
  induction error bounds met; wake-rate convergence recorded; Tier A/B KPI deltas
  REPORTED (V-05c), not forced to vanish.
- **E4.8** fs-bem screening preset + one-shot interference/residual derivation with
  the full cache key (geometry identity, operating grid, panel preset, ground mode,
  solver version, coefficient convention); slider drag uses schematic preview;
  commit cancels stale jobs. Depends E4.2.
- **E4.9** Early referee fixtures: attached-flow multisurface + prop cases run as
  E4.2/E4.5 land — discrepancy discovered before browser feel work. → feeds E10.1.

### E5 — Browser integration
- **E5.0** Versioned worker ABI + FieldSourceSnapshot. Depends E0.7, E3.5.
- **E5.1** fs-flyer-wasm API v1 (init/step/control/refusals). Depends E5.0, E3.*,
  E4.6a, E0.3.
- **E5.2** three.js consumes real state. Depends E5.1, E2.2.
  **MILESTONE: PHYSICS-SPINE FLYABLE BUILD.**
- **E5.3a** Authentic mechanical controls + device mappings. Depends E4.6b.
- **E5.3b** Historical-pilot demonstration mode. Depends E4.6c.
  **MILESTONE: HISTORICALLY MODELED FLYABLE BUILD.**
- **E5.3c** TrainingAssist + SAS accessibility tuning (cannot alter historical
  parameter identity).
- **E5.4** Both sites + launch options. Depends E3.4, E1.3.
- **E5.5** Results card + KPIs + ensemble-context historical comparisons.
- **E5.6** Runtime QoS governor, field/sweep throttling, field-age UI, immutable
  tier enforcement. Depends E0.6/E0.7/E5.1.

### E6 — Determinism, replay, E2E harness
- **E6.1** Replay envelope (frozen v1 schema), scrubber, ghost A/B.
- **E6.2** Startup self-test + four-quadrant+wasm golden CI (skeleton from E3.5).
- **E6.3** `e2e_wright_flyer.sh` + hostile twins + JSONL logging.

### E7 — Field visualization
- **E7.1** Field service + wasm API (ambient first; wake after E4.7).
- **E7.2** Glyph/streamline/vorticity/divergence renderers + probe gizmos.
- **E7.3** Force overlay + strip loads + probes with strip-charts.
- **E7.4** "Why it porpoises" view (flagship; depends E4.6b/c, E7.3).
- **E7.5** Lesson scaffolding.

### E8 — Experiments & evidence surfacing
- **E8.1** Worker-pool sweep engine with CRN ensembles + plots + CSV.
- **E8.2** Design panel v2: augmented eigenmode view, polar redraw, decomposed
  design-diff cards.
- **E8.3** ApplicabilityDomain UI + plain-language evidence layer + provenance
  inspector. Depends E10.2.
- **E8.4** (v1.5) Robust optimization (gated per §9.3).

### E9 — Sound & polish (stretch)
- **E9.1** AudioWorklet synthesis (corrected BPF math). Depends E5.2.
- **E9.2** Instant-photo, challenges, onboarding journey (§2.1 order).

### E10 — Reference plane & validation program
- **E10.0** Validation registry freeze: dataset partitions, metrics, provisional
  bands (§10.5), uncertainty distributions, historical ensemble protocol — locked
  BEFORE end-to-end results are inspected.
- **E10.1** Referee harness (begins incrementally via E4.9): batch re-runs at
  pinned configs; discrepancy receipts; optional correction tables under §4.2 rules.
- **E10.2** V-01…V-14 executed into fs-vvreg/vv-scorecard.
  **MILESTONE: EVIDENCE-BADGED BETA.**
- **E10.3** Tier D fs-lbm wind-over-terrain runs (Linux perf hosts) → V-04c.
- **E10.4** Cinematic export path (reuses h7xu5 machinery); one hero clip.
- **E10.5** fs-contact certified-contact replay pass receipt.

### Critical path

E1.2/E1.4/E1.5 → E4.0 → E4.1 → (E4.2 + E3.2a) → E4.3/E4.4a/E4.5 → **E4.6a
open-loop validation** → E5.0/E5.1/E5.2 (physics-spine flyable) → E4.6b/E4.6c/E5.3
(historically modeled flyable) → E6/E7/E10 parallel → E10.2 (evidence-badged beta).
Assets (E2) and terrain (E1.3) parallel the spine. Wake (E4.7) and field viz (E7)
gate the wow, not the flyable. E0.6/E0.7 precede physics tuning.

---

## 13. Risks & Mitigations

| # | Risk | L | Mitigation |
|---|---|---|---|
| 1 | Open-loop model misses the A4 pole structure | med | V-02a gates E4.6a before any pilot/feel work; E4.9 early referee fixtures |
| 2 | Canard mechanical model under-sourced (E1.5 thin) | med | qualitative envelope gates (sign/tendency) + explicit uncertainty; V-02b bands widen honestly rather than fake precision |
| 3 | Wake/field costs blow budget on low-end devices | med | tier ladder + QoS governor + induced-error-driven population control; Tier A is the contract |
| 4 | SAB/COOP/COEP hosting constraints | high (known) | dual artifacts + multi-worker transferable fallback from day 1 (E0.7) |
| 5 | 3-D model licensing | med | E2.1 blocking vetting; commissioned fallback |
| 6 | Historical numbers contested | med | [V] discipline, WindReference provenance, ensembles-not-points, E10.0 freeze |
| 7 | wasm/native numerical divergence | med | det:: doctrine + E3.5 subsystem digests + four-quadrant+wasm goldens early |
| 8 | Scope creep toward general flight sim | high | §1.4 non-goals; new-aircraft asks become v2 beads |
| 9 | Round-1 scope additions (mech controls, rotor dynamics, hybrid wake, ensembles) slip v1 | med-high | staged milestones: physics-spine flyable precedes historically-modeled flyable; structural beam model and SmoothedTangentPlane are explicitly droppable to v1.5 without breaking any shipped claim (claims are mode-gated) |
| 10 | Instability frustrates casual users | high | journey order (§2.1): ride-along → assist → authentic; "why it porpoises" turns failure into the lesson |
| 11 | Real-time contact jitter on landings | med | 240 Hz substep + regularized friction + E10.5 certified-contact backstop |
| 12 | Pilot-model identifiability (V-02c data thin) | med | pre-registered calibration/holdout split; if sources cannot support holdouts, V-02c demotes to a declared-hypothesis label rather than a validation claim |
| 13 | Correction-table misuse outside domains | low | §4.2/§10.6 rules; hostile twin for stale/out-of-domain table application |

---

## 14. Open Questions for Review Round 2

1. **Added-mass matrix content:** which terms enter M_added for a fabric biplane
   (wing panels, canard planes, prop disks?) and from what source — analytic flat-
   plate strips, panel-method extraction, or both with a discrepancy battery?
2. **PilotWrightModel structure:** is delay + lead-lag + neuromuscular lag +
   saturation + remnant the right minimal structure, or should we adopt a crossover-
   model formulation with explicit gain adaptation? What can E1's sources actually
   identify?
3. **Canard mechanics data:** if no quantitative hinge-moment/balance data survives,
   what is the honest fallback — parameterized envelope with sensitivity study, or
   fs-bem-derived hinge moments with declared model uncertainty?
4. **Coupled-solve formulation:** Weissinger-L vs nonlinear lifting line with
   decambering for the multisurface solve — which handles the canard's large
   deflections and near-stall behavior better at N ≤ 80?
5. **Hybrid wake conversion policing:** what invariant set is sufficient at the
   near→mid and mid→far conversions (circulation + impulse + ?), and how is the
   aircraft-region remainder bound computed cheaply enough to run per conversion?
6. **Ensemble size & budget:** how many ensemble members can the browser realistically
   run for the historical presets (background workers) vs how many does V-07's
   joint-region test need for stable coverage claims?
7. **fs-airfoil residual basis:** monotone splines vs constrained RBF vs Bernstein
   bases for the shape-preserving residual — which best supports the uncertainty
   intervals and OOD refusals?
8. **Replay schema v1:** exact field list to freeze now (per §4.5) — review the
   envelope for anything missing before the freeze bead lands.
9. **QoS governor policy:** precise thresholds/hysteresis for field-rate and render-
   quality degradation, and how they are surfaced to the user without alarm fatigue.
10. **Terrain aero-coupling:** with FlatPlaneExact as the v1 ground mode, how do we
    present flights near Big Kill Devil Hill (out of domain?) — refuse, warn, or
    SmoothedTangentPlane with prominent labeling?

---

## 15. Appendices

### A. Model equations (implementation-normative)

**A.1 Strip force build-up.** Per strip: local flow = atmosphere + wake induction +
images − body-point velocity; α_eff in the section frame incl. warp twist and
surface deflection; circulatory forces via the owned unsteady mode; profile drag
from fs-airfoil; induced effects via the coupled planform solve — the induced-α
bookkeeping is documented in fs-wing's CONTRACT with a double-count falsifier test.

**A.2 Unsteady kernels.** Rational approximations to Wagner Φ(s) and Küssner ψ(s)
with order selected by battery (Jones' Φ(s) ≈ 1 − 0.165e^(−0.0455s) −
0.335e^(−0.30s) as candidate); reduced time s from ds/dt = 2|U_local|/c;
apparent-mass terms live exclusively in the generalized added-mass operator
(§5.1.2). All constants dimensionless and cited (Fung; Leishman).

**A.3 Coupled planform solve.** Bound Γ on all surfaces (N ≤ ~80); influence from
bound + trailing systems of every surface + all flat-plane images; dense solve via
fs-la with factorization reuse while geometry is unchanged; deflection updates as
low-rank refresh where profiled as worthwhile.

**A.4 BEMT.** Per-annulus momentum/blade-element with Prandtl F = (2/π)acos(e^(−f));
Glauert low-J momentum correction; map construction at design commit (J ∈ [0, J_max]
grid) with runtime interpolation; typed refusal on non-convergence.

**A.5 Wall-compatible turbulence.** u = ∇×A with A built from mirrored/parity modes
such that u_z(z=0) = 0 identically; von Kármán amplitude shaping σ(z), L(z);
OU-evolved mode coefficients projected solenoidal; exact derivatives by term-wise
differentiation. Recurrence battery over the longest scenario duration.

**A.6 Ground images.** FlatPlaneExact: reflection through the fixed plane with
ω′ = det(R)·R·ω; image identities stable. SmoothedTangentPlane: one filtered global
plane; receipts carry origin/normal/τ/slope/max residual; never "exact."

**A.7 Longitudinal analysis.** Full augmented linearization (rigid + aero-memory +
actuator + rotor [+ pilot]) owns stability claims; the (u, w, q, θ) projection is a
labeled teaching view; x_np from ∂Cm/∂CL root-solve at valid trim only.

### B. Performance math (Round-1 corrected)

Wake budgets are specified per representation (near filaments / mid elements / far
multipoles) and controlled by an induced-error budget, not one global cap. Example
arithmetic at 36 strips, 40 Hz wake, near+mid ≈ 2–3k elements with binned
evaluation ≈ 190–260 MFLOP/s scalar — plausible against the 8.33 ms tick with the
measured E0.6 numbers as the real acceptance evidence (planning verdict only).
Tier B average ≈ 3.2 ms/frame; wake-update frames ≈ 4.2 ms (§7.2). Field-service
budgets quoted only with exact interaction counts (§7.2).

### C. Historical dossier seed (E1 completes)

Primary: Wright diaries/letters (LOC), Daniels photograph, McFarland *Papers*.
Secondary: Culick et al. (Caltech) stability analyses; Deters/Broughton/Selig
AIAA-2004-0211; Langley full-scale replica campaigns; Wright Experience propeller
reconstructions; NPS Wright Brothers National Memorial materials (incl. first-
flight accounts and glide-count narratives); Smithsonian NASM object records
(fabric artifacts, Flyer III figures); LOC *Dream of Flight* flight-4 account.
Terrain: USGS 3DEP, NOAA shoreline history, NPS base maps.

### D. Glossary

Advance ratio J = V/(nD). BPF: blade-passing frequency = B·RPM/60. Fixed-control /
free-control static margin: stability with controls held vs released (hinge-moment-
mediated). Backdrivability: a control surface's tendency to move under aerodynamic
hinge moment. Wagner/Küssner functions: circulatory responses to step α / sharp
gust. Effect ownership: the §5.2 rule that one physical effect has one model owner.
Applicability domain: the dimensionless region where a subsystem's claims hold.
Porpoising: the closed-loop pilot–aircraft pitch oscillation. Warp: roll control by
wing twist. Seqlock: sequence-number ring protocol for torn-read-free snapshots.

---

## Review round log

| Round | Reviewer | Date | Disposition |
|---|---|---|---|
| 0 | NobleLion / Claude | 2026-08-16 | initial comprehensive plan |
| 0.5 | self-audit (fresh eyes + executed wasm32 probes) | 2026-08-16 | corrected fs-contact role, fs-simd wasm claim, fs-time integrator choice; added §11.3 audit evidence, §11.4 optional-fidelity doctrine |
| 1 | GPT Pro Extended Reasoning (external) | 2026-08-16 | **architecture accepted; major physics-and-validation revision required and integrated**: longitudinal contract split (open-loop / canard mechanics / closed-loop), effect-ownership graph, coupled multisurface Tier A, flat-plane-exact ground images, generalized added mass, rotor dynamics + CT/CQ maps, wall-compatible atmosphere + ensemble historical presets, hybrid near/mid/far wake, unilateral rail release, structural claims mode-gated, identity quintuple, referee-plane rename, validation rebuilt on identifiability/holdouts/pre-registration, corrected budget arithmetic + perf acceptance contract, worker ABI/QoS protocol, fs-airfoil L2 crate, product-copy neutrality, "why it porpoises" flagship view, BPF audio fix, historical-claims table revisions |
| 2–4+ | — | — | pending; convert to beads only at steady-state |
