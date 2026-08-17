# COMPREHENSIVE PLAN: Real-Time Wright Flyer Simulation with FrankenSim

**Working title:** *First Flight — Kitty Hawk, December 17, 1903*
**Document status:** Planning-workflow ROUND 4 (closure audit integrated;
Round-4 verdict was "ROUND 5 REQUIRED" — six structural closure defects + twelve
polish groups + seven question decisions, all integrated below; Round 5 is a
narrow integration-verification and mechanical name-scan).
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
14. [Round-5 Verification Charter](#14-round-5-verification-charter)
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
resulting changes in induced drag, mass, rigging load, and stall distribution —
structural utilization appears only when the active structure mode supports the
claim. Take away the headwind and test whether the modeled launch rail remains
sufficient under the resulting configuration and atmosphere. Then read the numbers: separation
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
   ensemble-drawn December wind. The user watches the modeled flight unfold — from
   the wing, from the sand, from Daniels' tripod. If undulation grows, the replay
   identifies the open-loop, mechanical, perceptual, and closed-loop contributions;
   if it does not, the result card explains which sampled conditions suppressed it.
   A caption explains this is a modeled hypothesis, not a recording (Round-2 copy
   rule: even the porpoising is a computed outcome, not a promise).
3. **Fly with Training Assist.** The user takes the controls with the
   accessibility assist engaged; the assist is designed to make a completed
   flight substantially more accessible, but the actual outcome remains model-
   and input-determined. The results card compares their run to
   the historical record and to the hypothesis replay, and reports the actual
   distance distribution context rather than a single "beat this" number.
4. **Authentic Controls.** Invited, not forced: "Now try it the way Orville did."
   Raw mechanical controls, raw instability. Authentic control is expected to be
   difficult; the UI reports observed user outcomes only from anonymous aggregate
   evidence once it exists — and the results card explains *why difficulty is the
   point*, with a link to…
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
| **Training Assist** | user intent + an explicitly BOUNDED accessibility controller through the same actuator path; the HUD always shows the assist authority envelope, current contribution, saturation state, and cumulative user-vs-assist control work; results can generate a same USER-INTENT-input no-assist counterfactual replay — the compared source channel is PRE-ASSIST user intent, never the assisted actuator output; common-prefix and post-terminal semantics per §8.6 (labeled a model counterfactual, never "what you would have done") | hybrid assist, not historical |
| **Modern SAS** | full stability augmentation | explicitly ahistorical |

Input mappings (keyboard / mouse-drag "hip cradle" / gamepad) are device mappings
onto pilot force/intent, remappable, sampled and deterministically quantized at
simulation ticks (§6, item 2).

### 2.4 Visualization modes

1. **Vector glyph field** — GPU-instanced arrows on a user-positioned probe volume;
   live animated.
2. **Streamlines / pathlines / streaklines** — with the terms used *precisely*:
   streamlines integrate one frozen field snapshot; pathlines time-integrate markers
   through successive snapshots; streaklines use the persistent
   Lagrangian tracer service with continuous tick-addressed release — no dense
   field-history retention (Round-3 S-14).
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
9. **Pilot perception view** *(Round-2)* — the visual/vestibular cues available to
   the modeled pilot shown beside the simulator's true state: the gap between "what
   the aircraft is doing" and "what the pilot can perceive right now" becomes
   visible during delay, saturation, and low-altitude optic-flow changes.
10. **"Why it rolls and yaws"** *(Round-2)* — the lateral sibling of the flagship:
    warp command, loaded twist, rudder deflection, rolling/yawing moments, sideslip,
    roll rate, spiral mode, adverse-yaw decomposition, pilot lateral cue, and
    control reversals. A/B toggles the historical warp–rudder linkage while
    preserving the atmosphere realization and the common applied-input prefix
    (terminal divergence per ABComparisonReceiptV1, §8.6) — three-axis control
    was the actual invention, and this view teaches it.

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
| Aspect ratio | ALWAYS name the denominator: AR_plane = b²/S_one_plane ≈ 6.4; AR_system = b²/S_both ≈ 3.2 — bare "AR" is banned in UI and dossier (Round-2) | [V] derived |
| Gap between wings | ~6 ft (≈ chord) | [V?] |
| Camber | ~1/20 (as flown 1903) | [V] |
| Wing section/construction | digitized rib geometry + fabric-covered construction; relationship to the Wright wind-tunnel section identifiers established by E1 (NOT "single-surface", per Smithsonian fabric records) | [V?] |
| Anhedral geometry | measured droop from drawings/reconstruction | [V] |
| Anhedral design intent | possibly gust-response motivated; lesson copy states intent only after a direct Wright-era source (Round-2 demotion) | [V?] |
| Canard | biplane elevator: two planes, own span/chord/gap/stagger, hinge + balance axes, ~48 ft² total | [V?] full geometry from drawings, E1.5 |
| Canard arm | ~7 ft ahead of CG | [V?] |
| Twin rudder (aft) | source/convention range 20–21 ft² (Deters/Selig list 20); projected-vs-developed convention resolved by E1.9 | [V?] |
| Empty weight | ~605 lb = 274 kg | [V] |
| Gross weight (pilot ~145 lb) | ~750 lb = 340 kg | [V] |
| Engine | 4-cyl inline, ~201 in³, ~12 hp peak; MASS BY BOUNDARY convention (NPS: 170 lb engine; ~180 lb variously incl. accessories) — dossier defines dry/installed/cooling/ignition/drivetrain boundaries (E1.9); no single ambiguous "engine weight" enters the mass build-up | [V?] |
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
| 1 | Orville | 120 ft | 12 s | 10:35 a.m.; undulating; wind ~24–27 mph reported [V? pending instrument, reference-height, and averaging convention — WindReference/E1.8] |
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
Dec 17 wind 24–27 mph reported [V? pending WindReference disaggregation]. Sand z₀ ≈ 10⁻³–10⁻² m with sparse grass [V? E1.8]; two camp sheds.

**Huffman Prairie, Dayton OH (1904–05).** ~84-acre pasture; light variable winds →
catapult derrick from Sept 1904 (~1,600 lb weight; NPS historic-registration source says a 30-FT derrick — the ~16 ft figure is unsupported; effective drop/travel dimension unresolved [V?] E1);
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
effect at h/b ≈ 0.1–0.3 with the h and b conventions frozen by E1.9 (flat-plane
vortex-image system); near-stall cruise (separation margin + lagged separation
state); biplane and canard–wing interference (coupled multisurface solve, §5.2);
low-J propulsion — noting (Round-2 correction) that the held-on-rail Dec 17 state
begins at J ≈ 0.7–0.8 in the 24–27 mph headwind (n ≈ 5.83 s⁻¹, D ≈ 2.6 m), NOT at
J = 0; J = 0 coverage remains required for calm-air counterfactuals and static
fixtures (§5.3); gusty boundary-layer wind (wall-compatible solenoidal synthesis,
§5.4); and the LATERAL story — warp adverse yaw, the coupled warp–rudder control,
roll response, and the fast spiral instability Culick describes; NPS's flight-3
account records a side gust and lateral overcorrection, so the pilot task and its
validation are two-axis, not pitch-only (§5.1.4, V-12).

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
| A | browser, always | selectable identified modes — **(A0)** coupled nonlinear lifting-surface + 2-D indicial states (cheapest validated screening mode); **(A1)** coupled prescribed-wake attached-flow model with compact wake state / reduced state-space (preferred historical candidate — it is the cheapest model in which cross-surface influence takes TIME). The historical default is selected by V-08b1 + V-08b2 + E0.6 measurements, not by planning preference. Both: FlatPlaneVortexImageExact ground images; propulsion with rotor dynamics + two-way coupling; wall-compatible turbulence | the flyable sim; KPIs |
| B | browser, capable machines | Tier A + resolved hybrid wake (connected near-wake filaments, coarsened mid wake, multipole far wake) with admission-selected feedback schedule | wake realism and model-native resolved-wake visualization |
| C | native, offline | high-res unsteady VLM/panel (fs-bem+fs-fmm) + dense VPM; cross-checks A/B attached-flow, interference, and wake predictions inside its validated domain | referee |
| D | native, offline | fs-lbm D3Q19 wind-over-terrain (no aircraft in v1 — Round-1 Q8) | cross-checks atmospheric statistics under matched conditions |

Ladder invariants (Round-1 identity revision, §4.4): one `PhysicalScenarioId` flows
through all tiers; each tier/approximation carries a distinct `ModelId`; residual
**correction tables are disabled by default** — they are explicit identity-bound
modes, fitted on calibration pins only, evaluated on holdout pins, constrained for
symmetry/smoothness/conservation, and applied only inside their own applicability
domain (§10.6).

**Runtime mode immutability (Round-3 S-16).** After `RunIntentId` is minted, no
`ModelId`-affecting mode may change. Concretely: a `FlatnessCertificate` failure
REVOKES/downgrades the affected evidence claim while the admitted ground model
continues (selecting SmoothedTangentPlane requires a new run); a pruning-bound
audit failure TERMINATES the run with `PruningCertificateFailed` (restart with
pruning disabled = new ModelId); `WeissingerLLinear` is an admission-selected
fallback/debug mode, never an automatic nonlinear-failure path; QoS may never
activate a reduced physics schedule of any kind.

### 4.3 Process/thread topology & QoS (Round-1 revision)

```
Main thread:        three.js render loop; UI; input capture (tick-quantized)
Worker "sim":       wasm physics @ fixed 120 Hz (+admission-selected 40/60/80 Hz
                    full wake-advection schedule in Tier B)
Worker "field":     wasm field-sampling service (immutable FieldSourceSnapshots)
Worker pool (0–N):  parameter sweeps, replay re-runs
SharedArrayBuffer:  seqlock state ring, field buffers (when cross-origin isolated)
```

- **Global QoS governor (Round-2: a hysteretic state machine, not a priority
  list):** `Normal → Constrained → Critical`, fast metric window 2 s, recovery
  window 10 s. Normal→Constrained after 3-of-5 fast windows with sim
  service-time p99 > 5.5 ms or completion-lateness p99 > 1.5 ms → atomically
  activates ONE frozen Constrained presentation profile (sweeps paused, field
  rate halved, field density at the declared reduced profile, render effects
  trimmed) — no undocumented substates inside Constrained (Round-3 S-16).
  Constrained→Critical immediately on backlog > 1 tick, two consecutive deadline
  misses, or repeated transport starvation → freeze field recomputation, retain
  the last overlay, minimum presentation profile; pause physics with a typed
  refusal only if misses persist after all non-physics work is quiesced.
  Recovery one state at a time after 10 s with p99 < 4.5 ms, zero misses, zero
  backlog. UI: a small persistent badge — "visual analysis reduced; physics
  unchanged" — modal only when physics itself pauses. **The physics tier never
  changes during a run.** Thresholds provisional pending E0.6/E0.8; the state
  machine itself is frozen early.
- **Fixed timestep** decoupled from render; render interpolates. All schedules
  are integer-ratio locked to the 120 Hz sim tick: 240 Hz contact penalty
  updates; admission-selected 40/60/80 Hz full-wake advection; 120 Hz
  shedding; admission-selected 120/60/40 Hz propeller coupling; the perception
  schedule; every field/tracer refresh. The complete selected schedule table
  enters ModelId and ReplayEnvelopeV1 (Round-4 P-03).
- **Fallback without cross-origin isolation:** each wasm instance is single-threaded
  but sim/field/sweeps still run in separate Web Workers with pools of transferable
  preallocated ArrayBuffers (no per-frame allocation/cloning). Tier availability is
  selected by *measured capability* (E0.6 benchmarks), not by SharedArrayBuffer
  presence alone; the compatibility banner reports the actual disabled features.
- **Tab suspension:** on visibility loss, execution pauses at a complete simulation
  tick; never unbounded wall-clock catch-up on resume; the pause is presentation
  metadata, not simulation time.

### 4.4 Identity model (Round-1 revision)

Five public top-level identity kinds plus the immutable pre-final
`RunIntentId`, all `fs-blake3` domain-hashed with REGISTERED DOMAIN SEPARATORS
and canonical serialization (`identity-authorities.json` at introduction).
The closed preimage is `RunIdentityBasisV1 { physical_scenario_id, model_id,
artifact_id, physical_uncertainty_realization_id,
model_uncertainty_realization_id, accepted_tick0_state_digest,
input_trace_schema_id }`; `RunIntentId = H("fs-flyer/run-intent/v1",
canonical(RunIdentityBasisV1))` (Round-4 S-01):

| Identity | Contents |
|---|---|
| `PhysicalScenarioId` | aircraft design, site, initial conditions, weather *distribution*, pilot hypothesis, launch system |
| `ModelId` | tier, every approximation/fast mode + parameters, correction-table selections, discretizations, timestep, solver modes |
| `ArtifactId` | hash of the COMPLETE physics execution closure: wasm/native physics binary, physics-affecting JS glue, deterministic host-import implementation, host-ABI + memory-layout contract, schema adapters, deterministic-math tables, complete physics-data manifest, toolchain/build lock (pure rendering assets stay under `PresentationId`) |
| `InputTraceId` | hash of `InputTraceV1` { schema_id, end_tick_exclusive, ordered ResolvedAppliedEvent(channel, applied_tick, ordinal_within_tick, quantized_value) } — the trace EXTENT is part of the domain (two event-free runs stopped at different ticks differ), and acquisition-clock metadata stays in a replay attachment (Round-4 S-01) |
| `RunId` | `H("fs-flyer/run/v1", RunIntentId, InputTraceId)` — a closed two-input definition, never "hash of the quintuple" (which was self-referential) |

"Same physical scenario, different model" is thus expressible without pretending two
runs share a complete identity; cross-fidelity receipts explicitly bind the pair.

**Lifecycle (Round-3 revision).** `PhysicalUncertaintyRealizationId` is a pure
function of (physical distribution, realization seed, dedicated
`PhysicalRealizationAlgorithmId`) — NOT the whole ArtifactId, so an unrelated
binary rebuild cannot change the wind realization and common-random-number
comparisons survive artifact changes. `ModelUncertaintyRealizationId` derives
independently from the model-uncertainty spec + its algorithm id. An interactive
run's final `InputTraceId` cannot exist at launch, so neither can its `RunId`:

`RunIntentId = H("fs-flyer/run-intent/v1", canonical(RunIdentityBasisV1))` —
minted only AFTER prelaunch closes and tick 0 is frozen. When the trace closes,
`InputTraceId` and `RunId = H("fs-flyer/run/v1", RunIntentId, InputTraceId)`
are minted. `ModelUncertaintyRealizationId` derives from (model-uncertainty
spec, realization seed, `ModelRealizationAlgorithmId`) — both algorithm ids and
their canonical manifests are carried in ReplayEnvelopeV1 (Round-4). `RunAnchor = Intent(RunIntentId) | Final(RunId)`:
active subordinate objects bind to the anchor; final persisted objects mint new
final ids — no provisional identifier is ever relabeled. Physical uncertainty (wind, temperature, pilot mass, sourced
geometry variation) and model uncertainty (section residual coefficients, hinge
friction, wake closure) are SEPARATE identity-bound realizations; both enter
`RunId` (via the basis), but only physical uncertainty belongs to the
reconstructed scenario distribution — model deficiency can never masquerade as
weather. **Pilot-uncertainty allocation (Round-4 P-09):**
`PilotPopulationHypothesisId` (distribution artifact + pilot/day/flight
dependency structure) → PhysicalScenarioId; a concrete latent
pilot/day/flight member draw incl. human-variability remnant realization →
PhysicalUncertaintyRealizationId; perception geometry algorithms, controller
form, delay/filter discretization, gain family, saturation law, remnant
ALGORITHM → ModelId; explicit residual model-deficiency draws →
ModelUncertaintyRealizationId. H-07 manifests name each hierarchy level; a
single undifferentiated pilot_effect field is forbidden. Subordinate
identities never alter the quintuple: `FieldSourceSnapshotId = H("fs-flyer/field-source-snapshot/v1", RunAnchor,
source_tick, complete_source_state_digest, snapshot_schema_id)`;
`FieldQueryId = H("fs-flyer/field-query/v1", FieldSourceSnapshotId,
components, grid, induction_approximation, field_backend, precision,
derivative_and_masking_profile)` — bound to the actual sampled STATE, not just
the run (two ticks of one run are distinct — Round-4 S-02);
`PresentationId` (assets/camera/render/audio); `CheckpointId` (RunAnchor +
tick + complete state snapshot).

### 4.5 Replay envelope

`ReplayEnvelopeV1` (frozen by E0.9): the expanded `RunIdentityBasisV1`,
`RunIntentId`, the closed `InputTraceV1` + `InputTraceId`, the final `RunId`,
both realization-algorithm manifests; an artifact manifest (binary
hash, data/terrain/table/prop-geometry/correction hashes, dossier version); a
conventions block (axes, handedness, moment signs, reference areas,
span/semispan definitions, units schema, floating-point and deterministic-math
profiles); all integer-ratio schedules; the complete `CheckpointStateV1` at tick
0; the canonical applied input trace plus an `InputAcquisitionTrace` attachment
(device_sample_time, requested_tick, applied_tick, late_by_ticks,
sequence_number, input_transducer_mode); the
event/refusal log (liftoff, saturations, separation transitions, wake-mode
events, contact, touchdown, QoS attachments); optional full checkpoints with
mandatory digests; evidence receipts; migration history; optional signature.

`CheckpointStateV1` is COMPLETE internal state: rigid body, rail/contact, canard
mechanical, warp/rudder mechanical, aero-memory, separation, bound circulation,
physical wake + core state, rotor/engine/drivetrain, pilot perception+controller,
atmosphere modal OU state, RNG stream/counter manifest, work ledgers, subsystem
digests — PLUS every algorithmic history item capable of changing the next tick
(Round-3 S-05): nonlinear continuation branch + last accepted iterate,
reduced-wake scheduling state, prop-coupling predictor and relaxation history,
wake-feedback predictor, pruning-audit schedule/status, tangent-plane filter and
hysteresis state, integer-ratio schedule phases, admission-certificate states. A
cache may be excluded only after a checkpoint/resume hostile twin proves its
deterministic rebuild reproduces every subsequent subsystem digest. Presentation quality, QoS transitions, device metadata, and field
queries are replay ATTACHMENTS, excluded from `RunId` unless they altered the
applied input trace.

Migrations parse old envelopes, mint new identities, preserve originals, and
state which semantics changed — never silently reinterpret an old replay under
new physics. **Old-exact playback is a contract, not an aspiration (Round-3 S-18/Q4):**
RATIFIED default hosting (Round-4 Q5) — Cloudflare R2 Standard as the
content-addressed serving primary (native prefix-scoped R2 bucket locks after
promotion — R2's S3 API does NOT implement AWS object-lock headers, so the
native lock interface is used) + AWS S3 Standard with versioning + Object Lock
(compliance mode) in a SEPARATELY ADMINISTERED account as the WORM mirror;
git carries only signed trust roots, manifests, indexes, and small schema
fixtures (no LFS for production binaries). Replication is VERIFIER-MEDIATED
dual publication (independent uploads, read-back size+BLAKE3 verification,
retention controls, then TUF targets publication — never provider-native
replication; source deletions never propagate). SSE-S3 on the mirror unless
independently escrowed KMS retention exists (a deleted KMS key makes locked
objects unreadable). ~1 TiB dual-stored ≈ $39/month. E0.9 may replace a vendor
only via a written equivalence receipt (failure/admin domains, immutable
retention, conditional create, multipart integrity, audit export, retrieval,
regions, budget) — S3 API compatibility alone is not equivalence. Signing
uses TUF-style role separation (offline threshold root, delegated targets,
short-lived timestamp, snapshot metadata, rollback/freeze protection) — BLAKE3
identity is content identity, not publisher authenticity. "Exact artifact"
includes the physics wasm/native binary, exact JS glue, the deterministic
host-import ABI implementation, memory-layout/max-memory contract, schema
adapters, det-math tables, all physics data artifacts, convention manifests,
build provenance + toolchain lock, loader compatibility metadata, licenses, and
the startup golden + expected subsystem digests. Old executables run in a
capability-minimized no-network Worker/iframe sandbox (strict import
allowlists, bounded memory/execution, explicit CSP/CORP). Retention classes:
production/replay-referenced/receipt-referenced artifacts are PERMANENT
(deletion requires a signed policy epoch that explicitly retires the affected
old-exact claims — never a lifecycle job); unpromoted release candidates
180 days; CI candidates 90 days; corrupt uploads 14 days past incident closure.
CI replays ≥1 prior schema/artifact generation, plus corruption and
mirror-divergence hostile twins. **Trusted-loader boundary (Round-4 P-10):**
two capability domains — (1) the trusted CURRENT loader may access the
archive, validate TUF timestamp/snapshot metadata for discovery, retrieve
content-addressed packs, verify publisher signatures + BLAKE3 + lengths, and
create immutable read-only byte sources; (2) the ARCHIVED executable runs with
no network and only the frozen host imports needed to consume those verified
sources (a no-network sandbox cannot fetch its own lazy packs). Discovery
freshness is distinct from replay-pinned verification: an expired online
timestamp blocks fresh discovery but never invalidates a replay whose target
metadata, delegation chain, root epoch, and content hash are archived. A
replay without an archived artifact stays readable but cannot claim exact
re-execution.

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
  mechanical_stop_load_limit,
  warp_lever_ratio, warp_cable_compliance,

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
  PrescribedKinematicEstimated,   // kinematic warp, schematic deformation
  ReducedAeroelasticWarp {        // Round-2: the missing middle mode
    compliance_matrix_id, loaded_twist_basis,
    wire_slack_bounds, optional_first_order_lag,
  },
  QuasiStaticBeamAndRigging,      // spanwise beam + wire tension/slack + loaded warp
}

PilotHypothesis { pilot_max_force, pilot_arm_impedance, neuromuscular_parameters }

ModelSafetyLimits { max_canard_rate_guard, max_warp_rate_guard, iteration_caps,
                    condition_thresholds, max_wake_population }
// Numerical guards are part of ModelId and REFUSE when reached; they are never
// user-editable historical design parameters and never silently clamp motion.
```

`ReducedAeroelasticWarp` supports aerodynamic-control and lateral-mode claims
(loaded twist, effective control power, slack-risk diagnostics) but NOT component
structural margins; it is the minimum mode for a `Validated` V-12 lateral claim —
prescribed kinematics can demonstrate adverse yaw educationally but overstates
control effectiveness and shifts the roll/spiral modes.

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

The aerodynamic assembly returns `AeroGeneralizedLoads { q_rigid_nonaccel[6],
q_control_nonaccel[nc], m_added_rr[6,6], m_added_rc[6,nc], m_added_cc[nc,nc],
q_added_bias[6+nc], diagnostics }` where `nc` contains only ACTIVE DYNAMIC
control coordinates (a prescribed kinematic warp coordinate never enters the
mass solve; a dynamic canard hinge or declared dynamic reduced-warp coordinate
does). **Round-3 S-01 correction — the solve is the JOINT generalized system,
not the rigid six-vector with orphaned cross blocks:**

```
[ M_rigid + M_added_rr        M_added_rc            ] [ ν̇      ]
[ M_added_rcᵀ                 M_control + M_added_cc ] [ q̈_ctrl ]
        = [ Q_rigid_nonaccel; Q_control_nonaccel ] + Q_added_bias
```

via one declared deterministic block factorization or fixed-order Schur
complement, covered by goldens, BEFORE the Lie-group update. Cross-coordinate
accelerations are never lagged or finite-differenced. The whole-canard
generalized hinge load in `q_control_nonaccel` EXCLUDES acceleration-dependent
apparent-mass terms — those live only in the added-mass blocks and bias (no
double counting). For `PanelExtracted`, derivatives needed by `q_added_bias`
come analytically from the interpolation basis or a verified derivative
artifact; runtime finite-differencing of the table is forbidden.

**AddedMassMode (Round-2 ladder):** `AnalyticStrip { wing_planes, canard_planes,
vertical_surfaces, control_coordinate_blocks }` is the browser baseline
(cross-surface fluid-inertia omissions DISCLOSED in its no-claims block);
`PanelExtracted { extraction_artifact_id, interpolation_table_id,
ground_mode_domain }` is an additive identity-bound mode and discrepancy referee
over a pinned grid of canard deflection, height/attitude, and warp coordinates —
never a silent patch on the analytic matrix. **Propeller disks are excluded** from
rigid-body added mass: rotor inertia, blade unsteadiness, dynamic inflow, and
gyroscopic terms stay in `fs-airscrew`. Every assembled M_added is checked for
symmetry, positive semidefiniteness, positive definiteness of the total, frame
covariance, and WORK/ENERGY CONSISTENCY of Q_added_bias (a numerically
differentiated table would recreate the noise problem in another form).

**Partitioned integrator (Round-2):** exact discrete transitions for linear
aero-memory states; implicit midpoint/IMEX for cable/hinge/rotor states where
stiff; the energy-consistent effective-mass solve; a Lie-group midpoint
rigid-body update; event-localized ground contact. V-05a/V-05d gates include augmented
pole and closed-loop phase convergence, not only trajectories. Type adapters
between `fs-mbd`'s and `fs-geom`'s types are owned here (audited seam).

Parasite drag is a **component ledger** (pilot, engine/radiators, skids, wires,
uprights, struts, chains, misc.), each with area, orientation, Re-dependent
coefficient source, uncertainty, and power loss; the flat-plate aggregate remains a
separately identified fallback mode. At the Flyer's power margin, parasite-drag
error can decide whether it flies at all (V-13a/V-13b).

#### 5.1.3 Mechanical control system (Round-1 P0 revision)

The canard control path is a one-DOF mechanical system: pilot force/torque → lever
and cable ratio → canard hinge dynamics, including canard inertia, the
**whole-canard generalized hinge load returned by `fs-wing`** (Round-2 ownership
fix: fs-airfoil supplies section data; it cannot own the aircraft-level hinge
load, which depends on canard-plane interference, main-wing upwash, ground
images, gust penetration, unsteady circulation, and the actual hinge axis),
aerodynamic balance, friction, cable compliance, travel stops, and pilot-arm
impedance. Rate and force limits *emerge*; `ModelSafetyLimits` guards refuse.
Wing-warp and rudder use a separate linkage model with the
`LateralControlTopology` topology.

**Free-control stability under stiction (Round-3 refined):**
`ControlConstraintCase ∈ { FixedAtDeflection { deflection },
MechanismReleasedPilotDisconnected, HandsOnImpedance { pilot_arm_model_id } }` —
"free control" means `MechanismReleasedPilotDisconnected`; hands-on pilot-arm
impedance is a different physical boundary condition and is never reported as
the free-control margin. With Coulomb friction, backlash, or cable deadband the
released state is not one smooth scalar: analyses publish
`FrictionlessReleasedEquilibrium`, `SlidingBranch { direction }`, or
`StictionSet { deflection_interval, pilot_force_interval }`; a scalar
free-control static margin is UNAVAILABLE when the branch is not unique — an
interval/set-valued result shows instead (UX in §8.5).

**Evidence honesty:** Orville's account supports the SIGN and self-driving
tendency of the overbalanced elevator, not a quantitative hinge-moment curve;
the literature's ±30° travel figure was inferred from a photograph. Sign/tendency
may become `Validated` via geometry-driven modeling; force levels, friction,
inertia, balance sensitivity, and stops remain `Estimated` unless E1.5 or an
instrumented-replica campaign (candidate dataset A7) produces data — and no
flight endpoint may tune them (§10.4).

#### 5.1.4 Pilot models (Round-1 P0 revision)

- `PilotDirectHistoricalControls` — raw human input through the mechanical
  actuator model via a declared `InputTransducerMode`
  (`VirtualForceFromPosition { impedance_curve }` | `VirtualPositionCommand` |
  `ForceFeedbackDevice`): keyboards and ordinary gamepads cannot claim to measure
  pilot force, so the transducer model is shown in the UI and enters `ModelId`.
  No synthetic human delay is added; the measured software input-latency
  distribution is reported on the results card (§7.2.1).
- `PilotWrightModel` — autonomous historical-pilot hypothesis, **cue-based and
  multiaxis (Round-2)**: a `PilotPerceptionModel` (visual horizon attitude,
  near-ground optic flow, visual vertical motion, vestibular pitch/roll rate,
  specific-force cues, per-cue delays and noise, field-of-view/occlusion) —
  **computed entirely by a deterministic `fs-flyer::perception` service
  (Round-3 S-19 / Round-4 Q4): `PerceptionSchedule::EverySimTick120Hz` is the
  historical baseline (`Divisor2_60Hz` is an optional identified mode admitted
  only after V-16a/V-16b — a 60 Hz zero-order hold adds up to 16.7 ms cue
  delay; no 10 Hz pilot-crossover premise is assumed, the validation sweep
  merely covers through 10 Hz). Frozen causal order per update: sample
  accepted start-of-tick state → construct raw cue geometry → advance cue
  delay/noise/filter states → evaluate the pilot controller → apply pilot
  force through mechanics during n→n+1. Computed from rigid state, declared
  eye geometry, terrain, atmosphere state, and deterministic occlusion; the
  pilot NEVER consumes renderer pixels, GPU depth, presentation timestamps,
  LOD state, or device orientation; the renderer displays the cue-state
  output. `PerceptionModelId`, schedule, causal convention, and filter states
  enter ModelId and CheckpointStateV1** — feeding a
  fixed-gain longitudinal lead/lag controller constrained by a declared crossover
  envelope AND a lateral roll/heading/optic-flow controller acting through the
  hip-cradle warp–rudder linkage (Culick's fast spiral instability and the
  flight-3 lateral overcorrection make the pilot task two-axis), with shared
  neuromuscular lag, shared attention/force limits, and deterministic remnant.
  The historical pilot does NOT receive exact simulator q, γ, h, β, or φ — exact-
  state cues exist only as an explicitly ahistorical diagnostic mode. No online
  gain adaptation in the historical mode; `PilotAdaptiveEstimated` is a separate
  model identity. Gains draw from a pre-registered family at run start.
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
- **Prelaunch initialization (Round-3 S-03 correction):** `PrelaunchPhase ∈
  { HeldOnRailEquilibrated, SourcedTransient, DeliberateImpulsiveStart }`.
  Historical presets use `HeldOnRailEquilibrated` — NOT a point-residual solve
  of the whole state. **Round-4 S-03 revisions:** (a) *Canonical OU path.*
  Atmosphere uses `StationaryOuPathV1`: the exact stationary-covariance draw
  happens ONCE at a fixed anchor set by `PhysicalRealizationAlgorithmId`
  (v1: tick −3840 = −32 s at 120 Hz), then advances by counter-addressed
  exact-discrete innovations through tick 0 — so every permitted pre-roll
  window consumes a SUFFIX of the same physical path and yields the same
  tick-0 wind (changing a numerical pre-roll length can never change the
  physical realization; common-random-number identity survives). (b) *Branch
  awareness.* `PhysicalInitializationSpec { phase, sourced_transient_id?,
  deliberate_state?, control_equilibrium_selection }` is part of
  PhysicalScenarioId, with `ControlEquilibriumSelection ∈
  { UniqueReleasedEquilibriumRequired, DeclaredBranch { branch_id, provenance },
  DrawFromAdmissibleSet { set_id, distribution_id } }` — a user-selected
  stiction branch is a physical initial condition; a draw from the admissible
  set belongs to PhysicalUncertaintyRealizationId; a numerical initial guess
  may never implicitly select the branch. (c) *Partitioned alternate starts.*
  Within the SAME declared branch/member, hostile numerical starts must
  converge to the same accepted tick-0 state; starts deliberately on ANOTHER
  admissible branch must reproduce that branch or emit
  `PrelaunchBranchMismatch` (a set-recognition battery, not same-point
  convergence); all alternates share the identical StationaryOuPathV1 and
  realizations. Coupled pre-roll covers bound circulation, aero-memory, canard
  mechanics, engine, drivetrain, rotors, AND all active pilot/perception/
  controller states; resolved wakes converge only in protected-region
  observables. `PrelaunchInitializationPolicyV1 { canonical_ou_anchor_ticks:
  -3840, initial_pre_roll_ticks: 960, extension_ticks: 480, max_pre_roll_ticks:
  3840, closure_hold_ticks: 120, alternate_start_battery_id,
  closure_envelope_id }` enters ModelId; failure by 32 s emits
  `PrelaunchDidNotClose` (no band widening, no impulsive fallback); the
  accepted tick-0 digest enters RunIntentId/RunId. Provisional closure bands
  (E4.6d may tighten, never post-hoc loosen): normalized force/moment residual
  ≤ 2e-3; protected-region induction RMS ≤ 1e-3 / max ≤ 5e-3; circulation RMS
  ≤ 1e-3; rotor speed ≤ 5e-4 rel; canard ≤ 0.02° / 0.05°/s; rail-constrained
  pose ≤ 2 mm / 0.01°; identical continuation branch, stiction branch,
  admission/saturation/contact/refusal state, wake-topology class. (This is also where the J ≈ 0.7–0.8 held-on-rail advance ratio is
  established.)
- **Swept contact proxies (Round-3 S-17 revision):** the 240 Hz penalty update
  is preceded by deterministic swept point/segment/CAPSULE-vs-heightfield event
  localization for skids, canard frame/leading edge, wingtips, and
  **PHASE-RESOLVED propeller blade capsules** (rotor angle and angular velocity
  are checkpointed states; at ~330–350 rpm blade phase changes materially
  within one 1/120 s step, so a disk proxy false-triggers where neither blade
  is). **`BladeCollisionProxyV1` (Round-4 Q3):** 16 capsules/blade baseline —
  two spanwise rails at 25% and 75% local chord × eight radial intervals;
  deterministic refinement to ≤ 24 when the cover certificate fails; each
  radius = max distance from segment to its assigned provenance-bound
  blade-solid cell + digitization/interpolation/numerical margins (hand-entered
  radii forbidden); no capsule crosses the rotation axis or covers the hub
  void (hub/drivetrain is a separate proxy); acceptance = full cover, excess ≤
  max(5 mm, geometry uncertainty), radial segments ≤ 0.125 R; artifact id +
  cover certificate enter the manifests; > 24 needed → the blade-strike CLAIM
  refuses and only the disk warning remains. The propeller-DISK envelope
  remains a conservative clearance WARNING only. A blade strike, or any impact whose continuation would need an
  unavailable breakage model, emits `TerminalEvent::DamageModelUnavailable`,
  closes the physical run at the localized event time, and produces the impact
  report — cinematic continuation is a presentation attachment that cannot
  alter the physical replay. fs-flyer primitive, NOT fs-contact; the offline
  fs-contact pass stays the certified referee (V-20 covers disk-warning vs
  blade-contact separately, incl. disk-hits-but-blades-miss hostile cases).
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

**Planform layer (Round-2 revision).** Tier A uses a **warm-started NONLINEAR
low-order lifting-surface solve with section closure** over both main wings, both
canard planes, and the vertical surfaces — no scalar biplane factor in the
production force path. Main-wing and canard surfaces carry at least TWO chordwise
control rows in any mode claiming hinge moments or coupled unsteady phase (a
single row cannot give the hinge generalized force a meaningful pressure arm).
N ≈ 80 is a performance target, not an accuracy axiom — discretization is chosen
by the pinned convergence battery. `WeissingerLLinear` remains an exact fixture
and an admission-selected debug/lower-fidelity mode (never entered
automatically after nonlinear failure); decambering is a SEPARATE
`NonlinearDecamberingEstimated` mode with its own battery, never an unmarked
feature of the attached-flow solver. Every solve reports nonlinear residual,
condition estimate, continuation/branch identity, and iteration count;
factorization reuse only while the complete influence operator is unchanged
(canard deflection, warp, attitude/height images, and wake motion invalidate it
unless a verified update formula applies). On branch ambiguity beyond the
declared domain: deterministic continuation from the previous accepted state, or
typed refusal — never a silent jump to the lowest-residual branch. Surfaces enter
section regimes independently; warp-twist and adverse-yaw bookkeeping feed V-12.
Tier B raises resolution and adds the resolved wake; it does not introduce
multisurface coupling for the first time.

**Unsteady layer with a complete effect-owner record (Round-2 P0 revision).**
The Round-1 enum could not express every legal composition nor refuse every
under-owned one; ownership is now a RECORD with one owner per effect, checked at
admission (missing ownership and duplicate ownership are both typed refusals):

```
AeroEffectOwners {
  motion_circulatory: Indicial2d { kernel_id, state_order }
                    | PrescribedWakeReference3d { wake_rows, convection_model,
                                                  full_operator_id }
                    | PrescribedWakeReduced3d { full_operator_id, common_basis_id,
                                                reduced_order, scheduling_table_id }
                    | ResolvedFreeWake { wake_model_id },
  incident_gust:      KussnerSears2d { kernel_id, state_order }
                    | ChordResolvedIncidentField
                    | ResolvedWakeBoundaryCondition,
  noncirculatory:     AddedMassOnly,
  separation:         None | LaggedSeparationState { params },
  three_d_induction:  CoupledBoundSolve + declared wake owner,
  far_wake:           BoundedMultipole | AnalyticTrailing,
}
```

**Reduced time (Round-2 correction).** For a strip under a 2-D indicial owner,
`U_conv` = the POSITIVE CHORDWISE relative-flow component in the section frame,
and `ds/dt = 2·U_conv/c` — never the 3-D speed norm (a vertical gust must not
advance the motion-memory clock). The update uses midpoint U_conv with an EXACT
matrix-exponential transition over Δs; at U_conv = 0 the states freeze; reversed
chordwise flow or out-of-domain crossflow REFUSES the indicial owner rather than
hiding behind an absolute value. States initialize to the trim steady state
unless an impulsive start is explicitly requested (see PrelaunchPhase, §5.1.5).

**Tier-A candidates.** (A0) `Indicial2d` per strip over the coupled solve — the
cheapest screening mode, with the caveat that a small scalar step-response error
does not by itself validate use in a strongly coupled canard–biplane–ground
system; (A1) `PrescribedWakeReduced3d` — reduced from the separately identified
`PrescribedWakeReference3d` full-order operator (Round-3 S-07/Q1): production
reduction is a SHARED-BASIS block/tangential rational-Krylov method (IRKA-like
shift refinement) over the frozen operating grid — independently reduced
per-point models may not be interpolated (their state coordinates are not
comparable); balanced truncation at anchors and Loewner-from-samples remain
independent diagnostics. Candidate orders [4, 8, 12, 16, 24, 32]; smallest
order passing every gate wins; failure at 32 refuses A1 in that domain (bands
never widen). TWO receipts: V-08b1 (physics — full operator vs the independent
Tier-B/C referee) and V-08b2 (reduction — reduced vs full: complex MIMO
transfer incl. hinge-load channel, phase/group delay, poles/causality/
asymptotes, step/pulse/chirp/reversal shadows, grid-interpolation error,
state-continuation smoothness, short nonlinear trajectory shadows, full
identity manifest). Cross-surface convection DELAY is the point: the canard's
changed circulation reaches the main wing when the wake does. V-08b1 + V-08b2 + E0.6
pick the historical default; the plan does not pre-ordain it.
Kernel order is the smallest passing V-08a (lift AND pitching-moment channels,
initial value, asymptote, tail-relative error, causality, stable poles) — Jones'
2-state fit is a candidate, not the axiom. Near stall, a lagged-separation state
with hysteresis replaces the memoryless sigmoid (`Estimated`, honestly so). A
hostile twin enabling both resolved-wake circulation and full indicial motion
memory must refuse.

**Ground-effect layer (Round-2 rename + certificate).**
- `FlatPlaneVortexImageExact` (v1 historical default): every bound/wake element
  reflected across one fixed aerodynamic ground plane — the claim is precisely
  "exact satisfaction of the flat slip-wall condition FOR THE REPRESENTED
  SINGULARITY FIELD." No exactness is inherited by section viscosity,
  separation, body thickness, propeller–ground interaction, added mass (unless
  its ground terms are explicitly represented), or nonplanar terrain. Reflection
  uses the coordinate-free axial-vector rule **ω′ = det(R)·R·ω**; images retain
  stable identities and are EXCLUDED from physical wake circulation/impulse/
  energy ledgers (they are boundary devices, reported separately in
  boundary-work diagnostics). A `FlatnessCertificate` { aerodynamic plane,
  influence footprint, height/slope residuals, min clearance, expected
  ground-effect load error } gates every flat-ground evidence badge; near Big
  Kill Devil Hill the certificate decides: negligible-influence → continue with
  the coupling marked omitted-negligible; planar-enough → images; else
  SmoothedTangentPlane with an unmistakable `Estimated` label, or refusal of the
  ground-coupled CLAIM only (flight itself continues). The aerodynamic plane is
  visible in engineer/debug mode.
- `SmoothedTangentPlane` (additive `Estimated` mode): ONE continuously filtered
  global plane with hysteresis, slope limits, and a paired residual battery; its
  receipts report plane origin/normal/filter constant/max boundary residual PLUS
  (Round-3 S-21) plane translational/angular velocity, filter/hysteresis state,
  and the induced ARTIFICIAL BOUNDARY-POWER residual over the aircraft
  generalized coordinates — a moving image plane can inject/remove modeled
  work, which must be measured and bounded (V-06c), not assumed zero. Filter
  state is in CheckpointStateV1. It never inherits the word "exact."
- `HeightfieldBoundary` (deferred): reduced-rate source/panel boundary mode.
No per-element or abruptly switching tangent planes, ever.

**Wake layer (Round-2 hardened hybrid).** A 2,000-particle cap would discard
~98% of a 59-second wake; instead near/mid/far with an explicit temporal
contract:
- **WakeTemporalSchedule:** Kelvin-consistent shed circulation computed EVERY
  120 Hz sim tick (a fast canard reversal happens inside one 40 Hz wake
  interval); tick-level shed sheets may aggregate into macroelements only under
  a local moment/error battery; full free-wake advection at 40/60/80 Hz selected
  at run admission under a wake CFL/core-displacement certificate; wake-induced
  aircraft velocity between full evaluations via a deterministic first-order
  predictor (zero-order hold is an identified fallback, never the default).
  The schedule is fixed per run and enters ModelId.
- **Near wake:** connected vortex rings/filaments (topology-preserving; tip
  vortices are real filaments and the natural render primitive).
- **Mid wake:** deterministic age/distance conversion preserving Kelvin
  circulation closure per tube, connectivity and legal endpoints, hydrodynamic
  impulse, the far-field multipole moments through the retained order, core
  second moment, and declared reflection/geometry symmetries. Regularized
  kinetic-energy change is DIAGNOSED, not used as the primary conversion
  constraint (Round-2 Q5).
- **Far wake:** deterministic cell multipoles / analytic trailing structures;
  bounded tree/cell evaluation for aircraft feedback.
- **Core evolution:** `WakeCoreEvolutionMode ∈ { FixedCoreInviscid,
  DeterministicCoreSpreading { viscosity_model, growth_law },
  TurbulenceCoupledEstimated }` — a fixed radius over 59 s leaves unphysically
  concentrated old vortices; the mode enters ModelId. Visual age-fade is
  presentation-only and never changes Γ, core radius, or feedback.
- **Pruning:** cell-level bounding volumes + moment norms + source-to-region
  separation give a cheap multipole truncation bound over protected aircraft
  regions; periodic deterministic exact spot-evaluations AUDIT the bound; a
  failed audit TERMINATES the run with `PruningCertificateFailed` (restart
  with pruning disabled = new ModelId; pruning is never disabled in-run). Induction-error
  metrics use the mixed norm ‖u_fast − u_ref‖ / (U_ref + ‖u_ref‖) — pure
  relative error is meaningless where reference induction ≈ 0.
Flight-feedback and visualization evaluations may differ, each separately
identified; Tier B's role is "model-native resolved-wake visualization," never
"viz truth." Kernels build on `fs-vpm` (3-D extension upstreamed, §11.4).

#### 5.2.1 `fs-airfoil` — section models (new crate, L2; Round-1 P1 revision)

Generic section machinery, consumed by wing, canard, rudder, AND propeller (fixes
the L3→L3 layering violation): section geometry; analytic thin-section and
flat-plate reference modes; provenance-bound coefficient datasets (separate wing /
canard / rudder / propeller tables — the Wright 1901 tunnel tables are a
trend/convention validation source, NOT automatically a full-scale 1903 polar;
full-scale and reconstruction data dominate where available); **regime-partitioned tensor-product cubic B-spline residuals** with
coefficient-difference constraints over (α, log Re, δ) on the analytic baseline
(Round-2 Q7 decision: local support, interpretable influence, easy derivatives,
regime-specific shape constraints — global monotonicity is NOT imposed where the
physics is not monotone; Bernstein patches only at small transitions; RBFs
rejected for global support). Uncertainty is a COHERENT draw of spline
coefficients / low-rank function-space realization, never an independent interval
re-sampled per query;
normal/axial/moment representation through deep stall; explicit
attached/transitional/separated/post-stall regimes; uncertainty intervals on every
query; applicability-domain refusal instead of unconstrained extrapolation; indicial
kernel definitions and their exact-reference batteries.

### 5.3 `fs-airscrew` — propulsion (new crate, L3; Round-1 P1 revision)

- **BEMT** with Prandtl tip/root loss over provenance-bound **digitized radial
  geometry** (E1.6); diameter/pitch/activity-factor remain counterfactual levers.
  **Historical default (Round-2): the warm-started bounded station solve** with
  explicit high-loading/static closure and a convergence receipt — because CT/CQ
  are not single-valued in J alone under rotor acceleration, changing air state,
  radial Reynolds variation, and distorted inflow. The optional fast map uses
  `PropMapCoordinates { J, radial_reynolds_descriptor, inflow_harmonics[0..m],
  blade_roughness_id }` with the smallest dimension/harmonic order passing a
  paired error battery; a J-only table is permitted only where that battery
  passes. Maps cover J = 0 (calm-air counterfactuals, static fixtures, some
  Huffman launches) — but the Dec 17 held-on-rail state begins at J ≈ 0.7–0.8
  in the headwind (§3.4). Typed refusal on non-convergence.
- **Rotor dynamics:** `I_eq·Ω̇ = Q_engine(Ω, controls) − Q_prop,L − Q_prop,R −
  Q_drivetrain`, with separately diagnosable left/right shaft torques and optional
  failure states — no algebraic RPM jumps during the rail transient.
- **Two-way propeller–airframe coupling (Round-3 S-02 ownership fix):**
  `fs-wing` and `fs-airscrew` are sibling L3 crates and MUST NOT depend on each
  other — `fs-flyer` (L4) owns `CoupledPropAirframeStep`: (1) call fs-wing with
  the candidate prop-induced field; (2) project the returned flow into
  left/right disk inflow harmonics; (3) call fs-airscrew per propeller; (4)
  build the next actuator-disk/harmonic field; (5) apply the declared
  deterministic relaxation; (6) accept on the joint residual (disk-average +
  first-harmonic inflow, L/R thrust and torque, wing-load redistribution, total
  force/moment). Both L3 crates expose pure request/response interfaces.
  **Iteration scheme (Round-4 Q2 refinement):** previous-two-tick linear
  predictor → one unrelaxed evaluation → deterministic VECTOR Aitken
  relaxation over a dimensionless admission-frozen weighted residual.
  `PropCouplingSolverSpec` (candidate tuple, residual scales, denominator
  floor, growth guard, cap) enters ModelId. The Round-3 tuple {ω₀ 0.5, clamp
  [0.25, 0.80], cap 4} is CANDIDATE A of a pre-registered six-candidate family
  (clamps down to [0.10, 1.00], caps 3–5), selected by: zero nonconvergence on
  E4.5 calibration fixtures → all V-15/V-16a holdouts pass unwidened → min p99
  corrections → min p99 coupling service time → lexicographic tie-break; the
  cap-5 candidate is eligible only if its worst tick meets the
  service/lateness contract. Guards: denominator below floor/nonfinite → reset
  to ω₀; scaled residual grows > 25% → reject the correction, retry once at
  ω_min; second growth or cap exhaustion → typed
  `PropAirframeCouplingDidNotConverge`, NEVER a silent one-way switch. Targets
  (measured, not claimed): median ≤1, p95 ≤2, p99 ≤3 corrections. `PropCouplingSchedule ∈
  { EverySimTick (120 Hz, historical baseline), ReducedRate2 (60 Hz),
  ReducedRate3 (40 Hz), OneWayWingToProp }` — admission-selected, in ModelId,
  immutable in-run; 60 Hz is tried before 40 Hz, and reduced rates require
  V-15/V-16a discrepancy envelopes over rail/climb/reversal/asymmetric/PIO
  fixtures. Inflow ≥ disk-average + first azimuthal harmonic (pushers). The
  rigid body receives drivetrain reaction torque, unequal-prop torque, and
  rotor GYROSCOPIC moment; one gear-constrained rotor coordinate normally,
  separate shafts only in the compliance/failure mode.
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
- **Mean profile (Round-2 solenoidality fix):** a fetch-varying z₀ inserted
  into U(x, z; z₀(x))·e_x makes ∂U/∂x ≠ 0 — smoothing the transition only makes
  the divergence smooth, it does not repair continuity. `MeanWindMode`:
  `FlatSiteLogLaw { scenario_effective_z0, displacement_height,
  reference_height }` is the historical 1903 mode over the certified launch
  region; `FetchAdjustedMassConsistent { roughness_map, adjustment_solver_id,
  boundary_residual_budget }` is the later spatial mode. Pointwise/blended z₀
  variation may never be inserted directly into a horizontal profile while a
  solenoidal claim stands. Neutral stability is an explicit declared assumption;
  other classes are additive modes.
- **Turbulence (Round-2 space–time contract):** solenoidal synthesis via a
  wall-compatible vector potential (wall parity → u_z(0) = 0; vertical variation
  INSIDE the analytically differentiated potential). The COMPLETE post-projection Hermitian
  spectral tensor is fitted DIRECTLY to a pinned neutral-surface-layer target
  of the Mann-tensor class — auto-spectra, complex cross-spectra, two-point
  coherence, phase, Reynolds-stress ratios, surface blocking (Round-3 Q2:
  diagonal-only or selected-coherence fitting may initialize, never accept —
  unconstrained tensor directions get multiaxis gust loading wrong while
  marginal plots pass). v1 fits NEUTRAL conditions only; non-neutral classes
  are future modes with new fit artifacts; each modal truncation / wall-parity
  construction / height set / domain gets its own fit artifact (the finite
  basis changes the realizable tensor). Time structure: deterministic mean-advection phase
  φ_k = k·(x − U_adv·t) PLUS exact-discrete OU amplitude evolution
  a_{k,n+1} = ρ_k·a_{k,n} + σ_k·ξ_{k,n} with philox counter-addressed
  innovations; the OU STATE is sequential, checkpointed, and part of the
  determinism digests — `sample(x, tick)` never hashes a fresh amplitude per
  query. Recurrence battery covers the ~39-minute 1905 scenario.
- **Gust events:** spatially uniform transverse 1-cosine ramps are permitted
  directly; any spatially LOCALIZED gust is generated through the same
  wall-compatible modal basis — multiplying a divergence-free gust by an
  arbitrary spatial envelope is forbidden.
- **Thermals (Huffman, v1.5):** seeded convective plumes; off for 1903.
- **API (Round-2 air-state):** `sample_air_state(x, tick)` /
  `sample_batch(points, tick)` → velocity, analytic gradient, DENSITY, dynamic
  viscosity, temperature, pressure, component provenance, applicability
  diagnostics — Reynolds number and dynamic pressure must derive from the same
  provenance-bound air state as the velocity (E1.8 sources the historical air
  data). Optional fast mode `AtmosphereSampling::AffineLocalField { anchors,
  order }` with a paired exact-modal error battery over span-scale fixtures
  (§11.4).
- **Dec 17 preset:** an **uncertainty-conditioned ensemble** over documented mean-
  wind ranges, surface assumptions, and qualitative gust evidence. It never claims a
  recovered historical gust trace — no historical 10-m turbulence time series
  exists to match (Round-1 correction).
- Validation split: analytic construction (V-04a), statistical targets
  (V-04b1/V-04b2), cross-model comparison with Tier D as a *cross-fidelity
  receipt* (V-04c).

### 5.5 Field-sampling service (module in `fs-flyer`; Round-1 revision)

```
sample_field(snapshot_lease, grid_spec, component_mask) -> {
  u[], grad_u[], omega[], div_analytic[], div_finite_difference[],
  strain_magnitude[], q_criterion[], lambda2[], kinetic_speed_gradient[],
  validity_mask[], singularity_core_mask[], solid_exclusion_mask[],
  component_mask[], provenance[],
  meta { field_source_snapshot_id, source_tick, source_state_digest,
         source_modes, force_coupled_components,
         visualization_only_components, omitted_components,
         core_radius, export_precision }
}
```

Round-4 S-02: callers supply a LEASED immutable `FieldSourceSnapshotV1`, never
a raw tick as a state substitute (tick lookup is an internal replay convenience
that first reconstructs and publishes the identified snapshot). A stale
snapshot yields a valid HISTORICAL result under its own
`FieldSourceSnapshotId`, never mislabeled current. Persistent tracers bind the
ordered PAIR of snapshot ids for temporal interpolation.
```
```

Selectable components: mean atmosphere, turbulent atmosphere, gust event, bound
circulation, physical wake, ground images, propeller induced field,
visualization-only embellishments. The UI never labels a sum "total flow" when a
force-coupled component the active model supports (body displacement flow,
propeller slipstream, viscous wakes) is absent — omissions are named. The
propeller-slipstream visual is the actuator-disk/harmonic induced field; helical
blade-tip vortices are drawn only from a resolved prop-wake model or as
`visualization_only` cinema. The divergence overlay shows BOTH absolute |∇·u|
and normalized ε_div, masking normalized values where ‖∇u‖ is under the floor or
inside singularity-core exclusions. (λ₂ was promised in §2.4 but missing from
the Round-1 API — Round-2 consistency fix.)

- **No `grad_p_hat`** (Round-1 P0): a kinematically synthesized solenoidal velocity
  field does not determine pressure. Kinematic quantities only; Bernoulli proxy only
  on explicitly irrotational components, labeled; surface Cp only from models that
  actually produce it.
- Consumes an immutable **FieldSourceSnapshot** (tick, physical state, bound
  circulation, wake buffer references, image-plane state,
  `AtmosphereStateSnapshot { modal_state_or_immutable_handle, schedule_phase,
  state_digest }`, kernel/core parameters, model identities) — Round-3 S-06:
  seed+tick alone is FORBIDDEN for sequential atmosphere modes; the OU state
  itself (or an immutable handle to it) ships with the snapshot, and in the
  non-SAB path every coefficient needed for exact sampling is in the pooled
  pack (E0.7 measures source-snapshot bytes separately from output bytes).
  Results carry the source tick and are never presented as current beyond
  their staleness threshold.
- Ambient parts analytic; wake parts via the shared (identified) induction mode;
  the exact-vs-finite-difference derivative dual powers the divergence overlay.
- Budgeted, worker-isolated, `fs-exec` cancellation-checkpointed, QoS-governed.

### 5.6 Sound (stretch, E9; Round-1 corrected)

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
   via `fs-rand` philox rooted in the appropriate
   `PhysicalUncertaintyRealizationId` or `ModelUncertaintyRealizationId`, then
   domain-separated by registered stream id, subsystem/mode index, tick, and
   counter — `PhysicalScenarioId` alone is never an RNG root (it names a
   distribution, not a realization; Round-4 S-01); no
   wall-clock, no map-iteration order, no thread-order dependence; deterministic
   cell ordering in wake binning; fixed-shape reduction trees.
2. **Tick-addressed inputs (Round-2 completed contract).** `InputPacket {
   device_sample_time, quantized_value, requested_tick, sequence_number }`;
   worker clock synchronized to the main-thread monotonic clock at startup with
   periodic drift checks; `LateInputPolicy::ApplyNextEligibleTickAndFlag` — NO
   interactive rollback; the trace records requested_tick, applied_tick, and
   late_by_ticks, and replay uses applied_tick (replays reproduce what occurred,
   not what the UI intended). Active-tab pacing follows a monotonic target
   schedule integrating only fixed dt; at most a small declared number of overdue
   ticks may run in one burst — beyond that, pause with a typed performance
   refusal.
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
shipped; loader feature-detects. Bootstrap executable target < 8 MB gz; large terrain/aero/prop/validation/
replay artifacts are LAZY content-addressed packs whose hashes are in
`ArtifactId` — never silently counted inside the bootstrap target (Round-3). wasm memory is reserved to
its run-time maximum before typed-array views are exported; growth during a run is
refused.

### 7.2 The frame budget (Round-1 corrected arithmetic)

Tier A per 120 Hz step (planning estimates, mid-range laptop, scalar wasm):

| Module | Per step | Per 60 fps frame (2 steps) |
|---|---|---|
| fs-atmo sampling (~80 pts) | 0.10 ms | 0.20 |
| fs-wing coupled solve + strips (~80 unknowns, 36 strips) | 0.30 ms | 0.60 |
| fs-airscrew optional fast-map path (2 props, interpolation + rotor ODE) | 0.06 ms | 0.12 |
| historical warm-started station solves (excl. outer coupling corrections) | E0.6/E4.5-measured | — |
| ground contact (fs-flyer) + rail | 0.03 ms | 0.06 |
| added-mass solve + fs-mbd/fs-time integrate | 0.04 ms | 0.08 |
| KPIs, ring buffer, bookkeeping | 0.05 ms | 0.10 |
| **Fast-map-kernel planning SUBTOTAL** | **~0.58 ms** | **~1.16 ms** |
| influence assembly + nonlinear iterations | E0.6-measured | — |
| control/pilot/perception update | E0.6-measured | — |
| prop–airframe coupling iterations | E0.6-measured | — |
| 120 Hz wake shedding + core bookkeeping | E0.6-measured | — |
| snapshot publication + subsystem digests | E0.6-measured | — |
| non-SAB pack/copy (when applicable) | E0.7-measured | — |

The 0.58 ms row is a **fast-map kernel subtotal — neither a closed Tier A
total nor the historical station-solve baseline** (Round-4 P-06: no
historical-mode performance conclusion until station solves + outer coupling
iterations are measured):
it excludes influence assembly, nonlinear iterations, ground-image refresh,
canard generalized-load integration, exact aero-memory transitions, mechanical
controls, pilot/perception states, prop-coupling iterations, shedding
bookkeeping, publication, digests, JS/wasm boundary overhead, and contention.

Tier B arithmetic under the subtotal assumptions (Round-2 corrected rounding):
wake advance 2.5 × 40/60 = 1.6667 ms/frame; strip feedback 0.6 × 40/60 =
0.40 ms/frame; **Tier B average = 3.23 ms per 60 Hz frame**; **a frame containing
one wake update = 4.26 ms** (report 3.23/4.26, or 3.2/4.3 at one decimal — 4.2
was a rounding error). The per-TICK view matters more than frame averages
because the deadline is per tick: ordinary tick 0.58 ms; wake-update tick
0.58 + 2.5 + 0.6 = 3.68 ms; long-run mean tick ≈ 1.61 ms — all subtotals.

Field service must state its interaction model: dense 32³×2,000 = 65.5 M
interactions/refresh vs k-neighbor 32³×60 ≈ 1.97 M + the declared far-field.
Field MEMORY is also budgeted: u + ∇u + ω + div×2 + strain + Q + λ₂ + masks ≈
3 MB/snapshot at 32³ (structure-of-arrays typed buffers; provenance is
DICTIONARY-encoded with compact per-cell indices + bitmasks, never a rich
object per cell; the estimate names the exact enabled channels and index
widths) → ~9–12 MB for a leased ring, plus tens of MB/s of copy and GPU upload
at interactive refresh — bytes and transfers are budget lines, not only kernel
counts.

**These are planning estimates, not acceptance evidence** (FLOPs are a poor
proxy for rsqrt/regularization/memory/branching across browsers; deterministic
transcendentals may make the atmosphere row optimistic). E0.6 precedes physics
tuning.

### 7.2.1 Performance acceptance contract (Round-1)

Round-2: compute service time, scheduling lateness, and backlog are SEPARATE
metrics — a 3 ms tick scheduled 6 ms late still misses. Per device/browser
class: compute-service p99 ≤ 6.0 ms against the 8.33 ms tick; **hard lateness
gates (Round-3 S-15, provisional until E0.8 ratifies — but never report-only):
completion-lateness p99 ≤ 1.5 ms, p99.9 ≤ 4.0 ms, deadline-miss rate ≤ 1e-4
over the qualified active-tab soak, max consecutive misses ≤ 1;** zero backlog
> 1 tick; zero unbounded catch-up; a 10-minute warm/thermal soak, not only a
short benchmark. Render: with T_vsync the
MEASURED refresh interval, render p95 ≤ 0.80·T_vsync and p99 ≤ T_vsync;
skipped-frame rate and longest skip run reported (a fixed 16.67 ms p95 would
tolerate ~1-in-20 missed frames). Input latency decomposed: device sampling →
main→worker delivery → tick quantization → simulation → publication →
presentation. No sweep may increase sim deadline misses; startup benchmarking
selects presentation/field quality and never mutates the physics tier of an
active run. Correctness and capability tests are SEPARATE: (1) deterministic
golden self-test at full speed; (2) cold+warm capability benchmark; (3) a cached
capability profile keyed by ArtifactId, browser version, device class,
isolation mode, and presentation backend. Suite records scalar/SIMD,
isolated/contended, cold/warm, SAB/transferable variants.

### 7.3 Interop contract (Round-1 revision)

- **Versioned seqlock state ring:** header { abi_version,
  payload_layout_hash, payload_bytes, run_epoch, run_anchor_digest_prefix,
  model_id, tick, published_slot, sequence } — `run_epoch` increments on every
  worker reinitialization, and readers REJECT slots whose epoch, anchor
  prefix, layout hash, or size mismatch the admitted run before touching
  payload (a seqlock alone cannot stop a renderer consuming a complete STALE
  slot after restart; V-18 gains restart/slot-reuse ABA twins — Round-4 P-08);
  ALL header/sequence operations use Atomics; the
  writer marks the sequence odd, writes its owned slot, publishes even with the
  payload writes ordered before the publish store; readers acquire through the
  matching atomic load and retry on torn/odd sequence. Checksummed snapshots.
  Wasm relaxed-SIMD is excluded from deterministic tiers.
- **FieldSourceSnapshot LEASED ring, minimum 3 slots (Round-2 — a double buffer
  is unsafe while a reader holds a slot):** FREE → WRITING → PUBLISHED → LEASED
  → FREE via atomic state transitions; payload ownership is immutable while
  LEASED; the sim NEVER blocks on the field worker — with no free slot it skips
  publication and increments `field_snapshot_drop_count` (surfaced by the
  field-age UI).
- Field buffers zero-copy to render when SAB is available; GPU upload remains an
  explicit measured transfer. **Non-SAB fallback (Round-2 corrected):** a buffer
  backed by wasm linear memory cannot be transferred away while the sim uses
  that memory — each publication performs one explicit wasm→standalone pack/copy
  into a pooled transferable buffer, ownership transfers, and the buffer returns
  via an acknowledgement channel; copy bytes/time, pool starvation, and ack
  latency are measured (E0.7) and included in the acceptance budget.
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
practical and a tested `forceWebGL` WebGL2 fallback; WebGPU compute may
accelerate field visualization but NO v1 physics or validation feature may
require it. **Round-3 Q9 decision: the v1 canonical field backend is Rust/wasm;
WebGPU grid computation is strictly a v1.5 experiment (task E7.6)** — promotion
requires ≥2× END-TO-END speedup (incl. encoding, dispatch, sync, upload, any
readback) at the pinned workload, no readback on the ordinary render path,
backend parity through a precision/masking battery, device-loss fallback that
never touches physics, no force/validation/RunId-bearing consumer, and a
separately identified `FieldBackendId` in `FieldQueryId`. Static hosting + COOP/COEP; degraded mode works
without.

### 8.2 The aircraft asset

Source vetting (E2.1): Smithsonian 3-D digitization of the 1903 Flyer (license
verified in-task), NASA releases, or commissioned model; provenance in-app.
Pipeline: glTF 2.0, LODs, KTX2, rig (warp morph driven by warp state — under
`ReducedAeroelasticWarp` or `QuasiStaticBeamAndRigging` the active mode's
LOADED warp distribution drives the morph; canard
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

Instanced glyphs (≤ ~30k) from field buffers; GPU streamline integration
through the latest immutable snapshot; **persistent Lagrangian tracer service
(Round-3 S-14)** for pathlines/streaklines — tick-addressed release, compact
trajectory samples, validity state, and the two source snapshots needed for
temporal interpolation; a dense field-HISTORY ring is explicitly rejected
(10 s at 10 Hz would be ~300 MB) — dense snapshots remain only the short
leased publication ring, and replay scrubbing reconstructs tracer paths
deterministically from checkpoints + field-source state. Wake filaments render
directly with age-fade (they carry Γ); divergence overlay per §2.4.4; separate
composited overlay pass. Legends and transfer functions always visible.

### 8.5 HUD & instruments

Period triad (anemometer, stopwatch, engine-revolution counter — the Wrights' own
KPI set) + modern overlay (airspeed, altitude, α, load-factor exposure, pitch rate,
thrust, L/D, ledger strip). Engineer mode: augmented eigenmode view distinguishing
airframe / actuator / aero-memory / rotor / pilot-loop modes — not everything
labeled "short-period/phugoid"; the four-state longitudinal projection appears
only as a labeled teaching view. **Set-valued free-control UX (Round-3 Q6):**
never collapsed to a midpoint; lay mode shows a textured stability RIBBON over
admissible canard deflection plus one of four plain-language states
(`StableForAllBranches`, `MixedAcrossBranches`, `UnstableForAllBranches`,
`NoReleasedEquilibrium`) with the headline "when the pilot lets go, friction can
hold the canard at more than one angle — a range of behavior, not one margin";
engineer mode exposes the branch diagram, stiction/pilot-force intervals,
branch-conditioned margins/eigenvalues, best/worst case, the active branch
during replay, and whether a design change shrank/enlarged/split the set.
`HandsOnImpedance` is always displayed as distinct from released control. Results card per §9. The design-diff card (Round-2 honesty fix) distinguishes
direct additive ledger deltas, attribution-method-dependent deltas, and a
NONADDITIVE INTERACTION RESIDUAL, and names its method (fixed-order ablation,
symmetric two-order attribution, or Shapley approximation) — in a nonlinear
coupled system no decomposition is uniquely causal, and none is presented as
such.

### 8.6 Replay UI

Timeline with event ticks (liftoff by rail-reaction criterion, gusts, separation
flags, reversals, touchdown), camera presets (chase, wingtip, Daniels tripod,
onboard prone view, free), export of the full identity envelope. **Two A/B modes, never
conflated (Round-4 S-06 terminal-safe semantics):** `SameInputTrace` — the same
SOURCE applied-input schedule and `PhysicalUncertaintyRealizationId` are
replayed into both runs; the `ABComparisonReceiptV1 { comparison_mode,
source_input_trace_id, run_a_id, run_b_id, common_end_instant,
common_prefix_trace_digest, first_divergent_terminal_event,
common_prefix_kpis, whole_run_outcome_refs }` defines `common_end_instant =
min(run_a_end, run_b_end)` and requires bit-identical applied events over
[0, common_end_instant) — if one run hits a `TerminalEvent` first,
common-prefix KPIs stop there, whole-run cards remain shown with the explicit
note that the surviving continuation has no same-applied-trace counterpart,
and each run keeps its own closed `InputTraceId`. `HumanRefly` — identical
scenario seed with a new input trace: complete closed-loop comparison. They
answer different questions and are labeled accordingly.

---

## 9. Configuration Space, Experiments & KPIs

### 9.1 Scenario schema & the three design spaces (Round-1 revision)

```
FlyerScenario = {
  design_family: HistoricalReconstruction | WrightCounterfactual | FreeformTeaching,
  design: FlyerDesign,
  site, weather_distribution,
  realization_spec: SeededRealization { realization_seed },   // Round-2: the
  // physical realization is a pure function of (distribution, seed, algorithm
  // identity); no independently mutable weather_realization object exists
  launch, pilot_hypothesis, model_selection,
}
```

- **HistoricalReconstruction (Round-2):** a dossier-supported JOINT parameter
  model — exact drawing constraints, source-level covariance/dependency graph,
  coherent uncertainty draws. Span, chord, area, mass, inertia, and CG are
  correlated by construction; independent slider endpoints may not combine into
  geometry the dossier does not support.
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
| Fixed/released-control stability | fixed-control scalar margin where valid; released-control branch-conditioned or set-valued result under the declared ControlConstraintCase; NO scalar emitted when the released branch is nonunique |
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
  cacheable BEFORE execution by a distinct `RunSpecId` (the complete
  `RunIdentityBasisV1` preimage + `InputPolicyId`); completed results are
  indexed by final `RunId` — no partial field list is ever described as "all
  RunId ingredients" (Round-4 S-01).
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
| A7a | Instrumented canard/linkage dry-rig + wind-on campaign (lead-partner candidate: The Wright Experience, pending capability/data-rights audit) | quantitative control-mechanics evidence for V-02b |
| A7b | Optional instrumented human-in-loop replica/control fixture | independent generic-pilot evidence for V-02c1 |

A7 protocol, instrumentation schema, independence groups, and calibration/
holdout splits are pre-registered before collection; A7a and A7b are DISTINCT
independence groups, never credited interchangeably (Round-3 Q3/S-12). The v1
rule: no signed feasibility agreement + credible wind-on instrumentation path by
the E1.7 freeze → v1 keeps the existing `Estimated` ceiling, and no historical
endpoint substitutes.

Each dataset is partitioned BEFORE model calibration into: calibration subset;
held-out validation subset; convention/digitization uncertainty; applicability
domain; quantities for which it is independent evidence. **Evidence lineage
(Round-2):** every record carries `EvidenceLineageId`, source artifacts,
derivation steps, an INDEPENDENCE GROUP, permitted claims, forbidden claims,
calibration/holdout role, and convention uncertainty — agreement inside one
independence group is never counted twice (Deters/Selig drew on Ames full-scale
data and CMARC; overlapping-origin "agreements" are one datum). A2 is
partitioned by origin: A2-exp-fullscale, A2-exp-prop, A2-panel-derived,
A2-synthesized-stall — the last may constrain priors and qualitative trends but
may NOT serve as independent Wright-specific deep-stall validation (its own
paper says the stall region was synthesized from generic literature trends).
Tier C/D outputs are cross-model references until independently validated for
the exact quantity compared.

### 10.2 Validation hierarchy (Round-1)

1. Mathematical verification (exact/converged references).
2. Calibration (explicitly identified subsets only).
3. Held-out component validation (withheld α, Re, deflections, geometries, prop
   operating points).
4. Cross-fidelity discrepancy (Tier A/B vs independently validated Tier C/D).
5. Historical joint PRIOR-predictive and MANDATORY leave-one-flight-out
   predictive checks (full-data posterior-predictive views are diagnostic only —
   Round-3 S-11).
6. Performance qualification (device/browser matrix).
7. **Anti-vacuity and discriminative power (Round-2):** each empirical gate has
   ≥1 pre-registered DEFICIENT baseline that must fail or score materially worse
   — uncoupled canard/wing, no unsteady memory, total-speed reduced time, one-way
   prop coupling, non-advecting turbulence, overly broad historical priors,
   fixed kinematic warp, zero hinge feedback. A gate no deficient model can fail
   validates nothing.

Verification batteries follow workspace law: refusals tested at cap AND cap+1; no
vacuous limit checks; falsifier-style negatives per gate; per-strip oracles, not
only totals; the effect-ownership hostile twin (§5.2) must refuse.

### 10.3 The V-cases (Round-1 rebuilt)

Round-2 restructure: **V-cases** are falsifiable verification/validation gates;
**H-cases** are historical-compatibility checks that can never promote a
component to `Validated` (component evidence comes only from its own V-case).

| ID | Case | Gate |
|---|---|---|
| V-01 | Section & full-aircraft steady-load holdouts | coefficient/derivative/uncertainty-calibration metrics on held-out data |
| V-02a | Open-loop longitudinal derivatives & poles | full pole/derivative set vs A4 within declared reconstruction uncertainty; time-to-double reported |
| V-02b | Canard control mechanics | hinge-moment SIGN, self-driving tendency, control-force/travel response, stop behavior within the sourced mechanical envelope; quantitative levels remain Estimated UNLESS an independent A7a holdout receipt supports promotion per §5.1.3 |
| V-02c1 | Generic pilot–vehicle mechanism | crossover behavior, delay/lead-lag, remnant spectrum, saturation, PIO susceptibility vs independent human-control / instrumented-replica evidence |
| H-02c | Wright historical-pilot compatibility | qualitative overcontrol tendency, finite stabilizable region, undulation summaries, endpoint predictions; remains Estimated unless Wright-specific control data appear |
| V-03 | Propulsion maps | CT/CQ/static thrust/torque/RPM from J=0 through the envelope; η only where well-conditioned |
| V-04a | Atmosphere analytic construction | solenoidality, wall parity, exact derivatives, air-state consistency |
| V-04b1 | Atmosphere tensor construction | analytic finite-modal covariance/tensor vs the complete pinned target (scale-normalized Hermitian matrix error + separately gated autospectral log-ratio, complex-coherence, stress-tensor, principal-axis, PSD-violation, wall-residual, two-point coherence/phase, and OU temporal/convection errors) |
| V-04b2 | Atmosphere sampled statistics | the DECLARED finite-record estimator with independent realization seeds as statistical units (overlapping windows of one realization never count as independent), simultaneous 95% max-statistic bootstrap bands, frozen binning/windows/counts/floors |
| V-04c | Atmosphere cross-fidelity | Tier-D discrepancy receipt under matched conditions |
| V-05d | Joint generalized added mass | block symmetry; added-mass PSD + total effective-mass PD; frame covariance; NONZERO rigid↔control cross-block fixtures; analytic/interpolant derivative consistency; Q_added_bias work closure; deterministic block-factorization equivalence; checkpoint/resume equality |
| V-05a/b/c | Convergence & discrepancy | Round-2 unstable-system protocol: local-order fixtures, short-window state/force shadowing, open-loop pole + time-to-double convergence, augmented-loop gain/phase-margin convergence, event timing, KPI convergence, long-horizon ENSEMBLE distribution distance — full-duration trajectory RMS is diagnostic only / wake 20→40→80 Hz and admitted PHYSICS-schedule sensitivity (field rate belongs to V-18 noninterference) / Tier A-vs-B differences reported, not forced to vanish |
| V-06a/b/c | Ground effect | flat-wall implementation exactness / all six load components + derivatives vs h/b, pitch, roll vs high-res references / smoothed-tangent residual envelope in its slope domain |
| H-07 | Four-flight historical compatibility | ONE hierarchical joint model (shared day weather, pilot effects, flight realizations); PRIMARY gates = joint prior-predictive scoring + MANDATORY four-way leave-one-flight-out predictive scoring (a held-out flight never scores a distribution conditioned on itself); full-data posterior-predictive views are diagnostics; report region membership, volume/entropy, energy + declared proper scores vs pre-registered baseline, prior sensitivity, Monte Carlo uncertainty (Round-3 Q7: inference is FULL BAYESIAN OFFLINE with a frozen signed posterior/predictive artifact shipped to the browser — no in-browser MCMC/fitting; profile likelihood is a sensitivity check only) |
| V-08a | 2-D indicial kernels | lift AND moment time/frequency response, initial value, asymptote, tail error, causality, stable poles |
| V-08b1 | Full prescribed-wake physics | `PrescribedWakeReference3d` multisurface transfer matrices, canard↔wing phase, hinge-load channel, flat-wall residual, image symmetry, and h/b scheduling vs the independent E4.9b referee |
| V-08b2 | Prescribed-wake reduction | reduced-vs-full transfer (incl. hinge channel), phase/group delay, poles/stability, IMAGE-AWARE ground-scheduling interpolation + no-double-count proof, state continuation, nonlinear shadows on held-out operating paths |
| V-09 | Coupled biplane/canard loads | gap, stagger, α, control, h/b holdouts vs referee |
| V-10 | Wake invariants & induction | Kelvin closure, connectivity, impulse, retained multipole moments, core second moment, symmetry, MIXED-norm induction error, force/moment sensitivity, pruning-bound calibration, feedback-phase sensitivity |
| V-11a | Rail mechanics verification | unilateral complementarity, no tensile reaction, release criterion, work balance, timestep/event convergence |
| H-11b | Historical launch compatibility | rail position at release, ground-acceleration envelope, RPM and wind-conditioned advance ratio where independently supported (reaction-force HISTORY is a model output, not a historical datum) |
| V-11c | Prelaunch initialization closure | exact stationary OU covariance; canonical-prehistory invariance under window selection; branch/member ownership; within-branch alternate-start convergence; PrelaunchBranchMismatch recognition; protected-region wake convergence; mode-complete tick-0 state; checkpoint + RunIntentId equivalence |
| V-12 | Lateral control | adverse-yaw sign/magnitude, LOADED warp effectiveness, roll and spiral modes, warp–rudder phase, turn coordination; `Validated` requires ≥ ReducedAeroelasticWarp |
| V-13a | Total parasite drag & required power | totals vs independent full-aircraft evidence |
| V-13b | Component ledger closure | arithmetic/work closure, separately measured component fixtures, uncertainty propagation; ledger includes `unresolved_interference_drag`; total agreement never validates every allocation |
| V-14 | Browser real-time contract | §7.2.1 incl. service/lateness/backlog decomposition and contention soak |
| V-15 | Propeller–airframe coupling | thrust/torque, wing-load redistribution, inflow harmonics, iteration residual, one-way-vs-two-way discrepancy across rail/climb/asymmetric fixtures |
| V-16a | Actuator/integrator/latency phase fidelity | rate-limit, saturation, integration-order, and loop-phase reproduction |
| V-16b | Perception determinism & phase | cue geometry, delay/filter phase, checkpoint equivalence, and INVARIANT pilot-control traces across WebGL/WebGPU, render rates, dropped frames, camera LODs, and headless execution |
| V-17 | Atmosphere continuity & space–time structure | mean-field divergence, wall residual, full spectral tensor, two-point coherence/phase, exact OU variance, convection speed, checkpoint-addressed reconstruction equivalence (restore nearest identified checkpoint + deterministic advance == uninterrupted execution; NO stateless seed+tick OU reconstruction is claimed) |
| V-18 | Worker snapshot ownership, FIELD NONINTERFERENCE & input-tick semantics | leased-ring state machine, torn-read impossibility, late-input fidelity, checkpoint/resume identity, and BIT-IDENTICAL physics digests under varied field-query rate/density/backend, cancellation, snapshot drops, and QoS states (field-rate physics-sensitivity moved here from V-05b — the field is presentation and must not touch physics) |
| V-19 | Historical-score anti-vacuity | deficient baselines score materially worse; prior-width baseline comparison |
| V-20 | Online contact proxies vs offline certified contact | localized event instant + contact-feature agreement; blade-capsule cover certificate; disk-warning/blade-contact separation; adversarial disk-hits-but-blades-miss, hub-near-ground, high-Ω, terrain-ridge, and event-bracketing cases |
| V-21 | Replay identity & old-exact execution | RunIntent/final hash-preimage collision battery; trace-extent distinctions; complete artifact-closure hashing; hostile checkpoint/resume twins; signed-target verification; rollback/freeze/corruption/mirror-divergence; no-network sandbox + trusted lazy-pack injection |

### 10.4 Historical pass logic (Round-1 P0)

Historical cases are probabilistic because control and gust traces are unknown.
A pass requires: pre-registration of every uncertainty distribution (E1.7
freeze — Round-2 moved the freeze BEFORE component calibration); observed joint
outcomes inside declared predictive regions; sharpness (region volume/entropy)
and proper-score comparisons against a pre-registered baseline, not coverage
alone; no post-hoc distribution widening under the same validation identity;
component validation staying green; leave-one-flight-out where any flight
informed calibration. H-cases carry the H-prefix precisely because they cannot
promote components. No pilot/atmosphere parameter fitted on a case may be
credited on that case. Nothing anywhere is gated on "reachable."

**Two-stage campaign schema (Round-4 S-04 — the Round-3 single manifest was
circular: a prior-predictive campaign has no inference artifact, and final
RunIds cannot exist pre-execution):** the H-program runs on the Linux
performance hosts under existing deterministic lane protocols, sharded
content-addressably by (campaign_id × fidelity × member_index × flight_index),
resumable, merged in deterministic member-index order.
`HistoricalCampaignIntentManifestV1` (PRE-execution) binds: campaign id +
stage (`PriorPredictive` | `LeaveOneFlightOutPredictive { held_out_flight,
conditioning_artifact_id }` | `FullDataPosteriorPredictiveDiagnostic {
conditioning_artifact_id }`), frozen registry version, campaign model id,
scenario family id, artifact ids, requested member count, deterministic
PER-MEMBER allocations for day-weather, pilot-effect (Orville/Wilbur levels),
and flight-effect realizations (never one campaign-wide realization),
per-member×flight RunIntent SPECS, shard function, merge order, refusal
policy, and the Monte-Carlo stopping rule — `conditioning_artifact_id` MUST be
absent for PriorPredictive. `HistoricalCampaignReceiptV1` (POST-execution)
carries the completed member×flight RunIds, refusal accounting, actual member
count, deterministic merge digest, stopping receipt, and score/region artifact
ids — final RunIds appear only in receipts, never in intent manifests. **Surrogate rule:** no
surrogate-generated endpoint enters a final H-07 region, score, tail, or
coverage number — fs-surrogate may allocate exact runs, locate boundaries,
select fidelity, or serve as a control variate WITH an exact-simulation
correction and auditable error accounting; stratification + common random
numbers + exact reuse come first.

### 10.5 Provisional numerical bands (freeze in E1.7; audit conformance in E10.0)

Indicial kernels: max normalized step-response error ≤ 2%, gain ≤ 2%, phase ≤ 3°
over the declared reduced-frequency band, on BOTH lift and moment channels.
Airborne timestep convergence (unstable-system protocol): principal KPIs ≤ 1%;
open-loop pole real/imaginary parts within ratified bands; SHORT-WINDOW
trajectory norms within ratified local bands; long-horizon energy/Wasserstein
distance between matched ensembles inside a pre-registered Monte Carlo
confidence interval (full-duration RMS is diagnostic only). Contact
timing/impulse ≤ 2%. Atmosphere: spectral-tensor targets with declared
estimator/windows/bins/realizations/CIs (§10.3 V-04b1/V-04b2); TI and integral scales
±20%; Reynolds-stress ratios ±25%; gust quantiles ±20%. Wake fast modes:
mixed-norm induction RMS ≤ 3%, max ≤ 10% outside cores; circulation residual ≤
0.5%; impulse residual ≤ 1%. Real time: per §7.2.1 (service p99 ≤ 6 ms, lateness
and backlog gates). Ensembles (Round-2 Q6): browser presets 32/128/256–512
members with displayed Monte Carlo uncertainty and NO final coverage claims;
the offline H-program starts at 4,096 joint-day simulations, doubling until
pre-registered precision criteria (region-boundary uncertainty, score batch SE,
LOFO stability, tail quantiles) are met — validation counts are never limited by
what the browser can run. Historical endpoints: pre-registered joint regions —
no universal multiplicative band.

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

Native four-quadrant golden replay hashes ({aarch64, x86} × {debug, release})
PLUS two wasm-in-node lanes ({debug, release}) — six lanes total; golden-bump
protocol applies. Subsystem digests
(§6, item 4) localize divergence. `scripts/ci/e2e_wright_flyer.sh` cloned from the
hardened Euler cinematic runner: `--list/--check/--self-test/--run smoke/
--negative CASE/--replay`, bounded JSONL logging contract, hostile twins (config
vs identity tamper, input-trace truncation, seed mismatch, ledger violation
injection, stale correction-table identity, KPI-vs-recompute mismatch, wasm/native
golden divergence, terrain-hash drift, effect-ownership double-count and
UNDER-ownership, post-hoc distribution widening, deficient-baseline
false-passes, leased-ring protocol violations, late-input trace tampering).
Runner reuses production CLIs; never parallel logic.

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
| fs-atmo | L2 | wall-compatible solenoidal boundary-layer wind + gusts; analytic derivatives; ensemble presets; provenance-bound PRESCRIBED air-state fields consumed by aerodynamic and propulsion force paths (Re, q̄), but no prognostic thermodynamics or acoustic claims |
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
- **E0.1** Program root bead + EPIC/LEAF conversion rule (Round-3 Q10 +
  Round-4 Q7): every E-task becomes an EPIC; leaves have exactly one primary
  contract, one bounded surface, one owner, explicit input artifacts + deps,
  refusal tests, and one deterministic DONE-WHEN receipt; integration leaves
  are named explicitly, never hidden in the final child. **Conversion gate
  (structural, not just a list):** an epic may not enter Ready without leaf
  children when ANY of: estimate > 2 focused days; scope crosses a
  crate/layer, native/wasm, worker/main-thread, or model/presentation
  boundary; DONE-WHEN has > 1 independently falsifiable receipt; the task owns
  a device matrix, ensemble campaign, migration, archive, external dataset, or
  cross-fidelity battery; or the task mixes implementation + integration +
  UX/reporting. The named pre-split list is a MINIMUM SEED, not a whitelist:
  E0.6–E0.9, E1.1/E1.3–E1.7/E1.10, E2.2–E2.4, E3.1/E3.2/E3.2a/b/E3.3a/b/
  E3.4/E3.4b/E3.5, E4.0/E4.1/E4.2/E4.2b/c/E4.3/E4.3b1-3/E4.4a/E4.5/E4.6a-d/
  E4.7/E4.7b/E4.8/E4.9a/b, E5.0–E5.3c/E5.5/E5.6, E6.1–E6.3, E7.1/E7.1b/
  E7.2–E7.5, E8.1–E8.3b, E9.2, E10.0/E10.1/E10.2a-c/E10.2/E10.3/E10.5.
- **E0.2** `apps/wright-flyer` scaffold (Vite+TS+three.js, COOP/COEP dev server, CI).
- **E0.3** `fs-flyer-wasm` scaffold on the fs-wasm pattern; asupersync wasm profile
  feature pinned. → blocks all wasm integration.
- **E0.4** wasm32 CI guard lane over the flyer cone.
- **E0.5** fs-simd SIMD128 Tier-1w capsule — **NONBLOCKING, scheduled after
  E4.7 profiling; never on the critical path** (Round-2 tag).
- **E0.6** Browser performance microbench suite: det transcendental batches,
  40–100-unknown dense solves, BEMT loops, exact+fast Biot–Savart kernels, bin/tree
  traversal, SAB publication, transferable fallback, Float32 GPU uploads.
  DONE-WHEN: p50/p95/p99 across the device/browser matrix. → informs E4.2/E4.7.
- **E0.7** Worker transport & suspension prototype: seqlock ring, LEASED
  snapshot ring (≥3 slots), transferable pack/copy pools, visibility pause,
  no-catch-up, QoS throttling against a synthetic 120 Hz load. → E5.0.
- **E0.8** Worker timing/transport semantics: monotonic clock sync,
  requested/applied input ticks, backlog + bounded-catch-up policy, 10-minute
  contention soak. → E5.0, V-18.
- **E0.9** Replay schema + artifact archive/trust freeze: RunIntent/final
  lifecycle, complete CheckpointStateV1 (incl. algorithmic history),
  artifact-object schema + host ABI, primary+WORM-mirror layout + vendor
  ratification, permanent/ephemeral retention classes, TUF-style signing roles
  + key rotation, loader sandbox, migration semantics, backward-playback CI
  fixture + corruption/mirror-divergence twins. → blocks E3.2 replay work and
  E5.0. (PRE-SPLIT per E0.1.)

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
- **E1.7** Validation/evidence registry freeze (the former E10.0, moved BEFORE
  any component calibration; depends the time-boxed E1.10 go/no-go — `A7NoGo`
  is a valid nonblocking completion freezing the Estimated ceiling):
  A1–A6(+A7) lineage graph + independence groups,
  calibration/holdout partitions, claim permissions, uncertainty distributions,
  provisional bands, anti-vacuity baselines, H-case protocol. → blocks E4
  calibration.
- **E1.8** Air-state & weather dossier: temperature, pressure, density,
  viscosity, stability assumption, roughness/displacement-height uncertainty.
- **E1.9** Geometry-convention dossier: system vs per-plane aspect ratio, span
  vs semispan symbols, ground-effect h and b reference conventions,
  rudder/canard area conventions, engine/accessory mass boundaries.
- **E1.10** A7 feasibility + pre-registration (NONBLOCKING): lead-partner
  capability/data-rights audit, instrumentation plan, A7a/A7b separation,
  calibration/holdout freeze, explicit v1 go/no-go — no agreement means the
  Estimated ceiling stands.

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
  bit-identity (replay/checkpoint leaves depend on the E0.9-frozen schemas;
  integrator-only leaves may begin earlier). Includes the type-adapter seam.
- **E3.2a** Generalized added-mass assembly (AnalyticStrip baseline; panel
  cross-terms deferred to E4.2c). DONE-WHEN: acceleration-dependent fixtures
  converge without FD noise; symmetry/PSD/covariance/work batteries green;
  effective-mass solve admissible over the reference domain. → blocks E4.6a.
- **E3.2b** Partitioned integrator + time-scale certificate: exact aero-memory
  transitions, implicit control/rotor updates, augmented-pole phase convergence,
  stiffness admission tests. → blocks E4.5, E4.6a.
- **E3.3a** fs-atmo foundation: FlatSiteLogLaw, wall-compatible potential,
  exact analytic derivatives, sourced air-state API, deterministic stream
  partitioning, V-04a. Depends E1.7/E1.8.
- **E3.3b** Complete Mann-class post-projection tensor fit, exact-discrete
  sequential OU state (StationaryOuPathV1 anchor), checkpoint-addressed
  reconstruction equivalence, recurrence battery, V-04b1/V-04b2, V-17.
  Depends E3.3a. (PRE-SPLIT.)
- **E3.3c** Optional FetchAdjustedMassConsistent mode + boundary-residual
  battery (pointwise/blended z₀ insertion is NOT an implementation task).
- **E3.4** Rail (unilateral release) + fs-flyer contact + terrain queries +
  the `PrelaunchPhase` SCHEMA and deterministic constrained-state scaffold
  (integrated aero/propulsive equilibration is NOT claimed here — that is
  E4.6d). DONE-WHEN: dolly acceleration and release location converge under
  refinement; no tensile reaction; landing impulse/penetration/friction work
  converge. Depends E1.3.
- **E3.4b** Swept critical-feature heightfield event localization (skid,
  canard, wingtip proxies; PHASE-RESOLVED blade capsules + the separate
  prop-disk clearance warning) + the V-20 harness. Depends E1.6, E3.4, E3.2b.
- **E3.5** Structured determinism checkpoints (per-subsystem tick digests) — early,
  before physics churn.

### E4 — Aerodynamics & propulsion
- **E4.0** fs-airfoil crate (analytic modes, tables, constrained residuals,
  uncertainty, regimes, refusals). Depends E1.1/E1.4.
- **E4.1** Wing/canard/rudder/prop section datasets + models on fs-airfoil.
  Depends E4.0, E1.5, E1.6.
- **E4.2** Coupled Tier-A NONLINEAR multisurface lifting-surface solve (≥2
  chordwise rows on wing+canard; classical factors = fixtures/fallback;
  condition/continuation/invalidation contracts). Depends E4.1, E0.6.
- **E4.2b** Whole-canard generalized hinge-load interface (fs-wing →
  mechanical controls). Depends E4.2.
- **E4.2c** Optional panel-extracted added-mass cross terms + discrepancy
  battery. Depends E4.2, E3.2a.
- **E4.3** Complete `AeroEffectOwners` contract: chordwise reduced time with
  exact Δs transitions, exact-reference kernels (V-08a), ownership admission
  refusals, separation lag. Depends E4.2, E3.3.
- **E4.3b1** Full-order `PrescribedWakeReference3d` operator INCLUDING the
  complete FlatPlaneVortexImageExact contribution for every prescribed wake
  row, + a frozen operating/linearization grid explicitly spanning h/b, pitch,
  roll, canard deflection, warp coordinates, and convection parameters.
  Depends E4.3, E4.4a.
- **E4.3b2** Shared-basis rational-Krylov reduction of the IMAGE-AWARE full
  operator (no post-reduction image bolt-on), order ladder, projected-operator
  scheduling incl. ground scheduling, reduction receipt (V-08b2). Depends
  E4.3b1. (PRE-SPLIT.)
- **E4.3b3** V-08b1/V-08b2 campaign + A0/A1 historical-default selection.
  Depends E4.3b2, E4.9b, E0.6.
- **E4.4a** FlatPlaneVortexImageExact COMMON image kernel, identity,
  certificate, and ledger contract — closes the A0 bound-system AND
  analytic-trailing influence paths; every Tier-A influence operator consumes
  the same coordinate-free image transform and stable image identities
  (Round-4 S-05). Depends E1.3, E1.9, E4.2. Prescribed-wake images close in
  E4.3b1; resolved physical-wake images close in E4.7b.
- **E4.4b** (v1.5, NONBLOCKING) SmoothedTangentPlane optional mode, only
  after E4.4a and V-06a/b green — v1 reacts to a failed FlatnessCertificate by
  revoking the ground-coupled claim, never by switching modes.
- **E4.5** fs-airscrew kernels + the fs-flyer `CoupledPropAirframeStep`
  orchestrator (Round-3 S-02): warm-started station solve (historical default),
  multidimensional fast map, rotor/drivetrain dynamics, sibling-crate pure
  APIs, predictor + deterministic Aitken relaxation, admission-fixed coupling
  schedules, reaction/gyroscopic moments, engine, component power ledger,
  air-state inputs. Depends E1.6, E4.0, E3.3a, E3.2b, E4.2. Feeds V-15/V-16a.
  (PRE-SPLIT.)
- **E4.6a** Open-loop integrated aircraft. Depends E3.2a, E3.2b, E4.2, E4.3,
  E4.4a, E4.5. DONE-WHEN: the source-matched FIXED-CONTROL V-02a
  derivative/pole gates pass before any pilot exists (free-control gates
  belong to E4.6b — Round-3 S-09 ordering fix).
- **E4.6b** Canard/warp/rudder mechanical controls incl. free-control
  branch/set semantics. Depends E1.5, E4.6a, E4.2b. DONE-WHEN:
  released-pilot-disconnected, sliding, and stiction-set analyses pass;
  hands-on impedance separately labeled.
- **E4.6b0** ReducedAeroelasticWarp mode — required before any `Validated`
  V-12 lateral claim. Depends E4.6b.
- **E4.6c** PilotPerceptionModel (deterministic fs-flyer::perception service,
  V-16b) + two-axis PilotWrightModel + TrainingAssist (calibration subset only;
  assist tuning isolated; E1.7 decides Validated-generic vs
  Estimated-historical claims). Depends E4.6b; historical lateral completeness
  additionally needs E4.6b0. (PRE-SPLIT.)
- **E4.6d** Integrated `HeldOnRailEquilibrated` closure: StationaryOuPathV1
  draw, aero/control/rotor pre-roll, branch-aware alternate-start battery
  (V-11c), tick-0 digest. Depends E3.3b, E3.4, E4.2, E4.3, E4.5, E4.6b;
  depends E4.3b1/E4.3b2 only for A1; depends E4.7/E4.7b only for Tier B; and
  depends E4.6c for any pilot mode with dynamic perception/controller/remnant
  state (the mode-complete tick-0 state freezes BEFORE RunIntentId mints —
  Round-4 S-03). (PRE-SPLIT.)
- **E4.7** fs-vpm hybrid wake: connected near wake, 120 Hz shedding,
  invariant-preserving coarsening, multipole far field, core-evolution modes,
  bounded pruning + audit, exact-vs-fast batteries. Depends E4.2, E4.3, E4.4a
  + the additive fs-vpm 3-D paths. DONE-WHEN: topology/circulation batteries
  green; coarsening invariants hold; induction bounds met; wake-rate
  convergence recorded; Tier A/B KPI deltas REPORTED (V-05c). (PRE-SPLIT.)
- **E4.7b** Physical-wake IMAGE integration: image identity, ledger exclusion,
  conversion symmetry, V-06a/V-10 cross-battery. Depends E4.7, E4.4a.
- **E4.8** fs-bem screening preset + one-shot interference/residual derivation with
  the full cache key (geometry identity, operating grid, panel preset, ground mode,
  solver version, coefficient convention); slider drag uses schematic preview;
  commit cancels stale jobs. Depends E4.2.
- **E4.9a** Early STEADY referee fixtures: attached-flow multisurface + prop
  cases run as E4.2/E4.5 land — discrepancy discovered before browser feel
  work. → feeds E10.1.
- **E4.9b** Independent UNSTEADY prescribed-wake referee: high-resolution
  Tier-C canard↔wing impulse/step/chirp/reversal + MIMO transfer fixtures,
  incl. hinge-load and flat-ground cases, built on the exact dense 3-D
  filament/panel reference leaf (an early fs-vpm E4.7 leaf) — NEVER the A1 FOM
  implementation or its reduced basis. Emits the V-08b1 receipt (Round-4
  S-05).

### E5 — Browser integration
- **E5.0** Versioned worker ABI + leased FieldSourceSnapshot ring. Depends
  E0.7, E0.8, E0.9, E3.5.
- **E5.1** fs-flyer-wasm API v1 (init/step/control/refusals). Depends E0.3,
  E5.0, E3.1, E3.2, E3.2a, E3.2b, E3.3a, E3.4, E3.5, E4.6a.
- **E5.2** three.js consumes real state. Depends E5.1, E2.2.
  **MILESTONE: PHYSICS-SPINE RUNNING/REPLAYABLE BUILD.**
- **E5.3a** Authentic mechanical controls + device mappings. Depends E5.2,
  E2.4, E4.6b. **MILESTONE: FIRST HUMAN-FLYABLE BUILD.**
- **E5.3b** Historical-pilot demonstration mode. Depends E5.2, E2.4, E4.6b0,
  E4.6c, E4.6d, and the A0/A1 historical-default selection (E4.3b3).
  **MILESTONE: HISTORICALLY MODELED FLYABLE BUILD.**
- **E5.3c** TrainingAssist + SAS accessibility tuning (cannot alter historical
  parameter identity). Depends E5.3a, E4.6c.
- **E5.4** Both sites + launch options. Depends E3.4, E1.3.
- **E5.5** Results card + KPIs + ensemble-context historical comparisons.
- **E5.6** Runtime QoS governor, field/sweep throttling, field-age UI, immutable
  tier enforcement. Depends E0.6/E0.7/E5.1.

### E6 — Determinism, replay, E2E harness
- **E6.1** Implement the E0.9-frozen replay/checkpoint schema: scrubber,
  artifact retrieval, SameInputTrace (ABComparisonReceiptV1) + HumanRefly
  ghost modes. Depends E0.9, E3.2, E5.1, E5.2.
- **E6.2** Startup self-test + six-lane golden CI (skeleton from E3.5).
  Depends E0.9, E3.5, E5.1.
- **E6.3** `e2e_wright_flyer.sh` + hostile twins + JSONL logging. Depends
  E6.1, E6.2.

### E7 — Field visualization
- **E7.1** Field service + wasm API (ambient first). Depends E5.0, E5.1,
  E3.3a, E4.2; the physical-wake branch additionally depends E4.7/E4.7b.
- **E7.1b** Persistent Lagrangian tracer service: pathlines, streakline
  release, snapshot-pair temporal interpolation, compact retention,
  checkpoint/replay reconstruction, cancellation, memory caps. Depends E7.1,
  E6.1.
- **E7.2** Glyph/streamline/vorticity/divergence renderers + probe gizmos.
  Depends E7.1, E5.2.
- **E7.3** Force overlay + strip loads + probes with strip-charts.
- **E7.4** "Why it porpoises" view (flagship; depends E4.6b/c, E7.3).
- **E7.4b** "Why it rolls and yaws" view (depends E4.6b0, lateral pilot model,
  augmented lateral linearization).
- **E7.5** Lesson scaffolding + pilot-perception view.
- **E7.6** (v1.5 candidate) WebGPU field-compute experiment per §8.1 promotion
  gates; CPU/wasm remains canonical.

### E8 — Experiments & evidence surfacing
- **E8.1** Worker-pool sweep engine with CRN ensembles + plots + CSV.
- **E8.2** Design panel v2: augmented eigenmode view, polar redraw, decomposed
  design-diff cards.
- **E8.3a** Evidence/applicability plumbing + empty-receipt UX (early).
- **E8.3b** Populate final validated receipts. Depends E10.2.
- **E8.4** (v1.5) Robust optimization (gated per §9.3).

### E9 — Sound & polish (stretch)
- **E9.1** AudioWorklet synthesis (corrected BPF math). Depends E5.2.
- **E9.2** Instant-photo, challenges, onboarding journey (§2.1 order).

### E10 — Reference plane & validation program
- **E10.0** Registry-CONFORMANCE audit (the freeze itself moved to E1.7):
  confirm no partition, metric, prior, or band changed after protected results
  were observed.
- **E10.1** Referee harness (begins incrementally via E4.9): batch re-runs at
  pinned configs; discrepancy receipts; optional correction tables under §4.2 rules.
- **E10.2a** H-07 Bayesian INFERENCE ENGINE + signed conditioning artifacts:
  the frozen prior/likelihood (or simulation-based inference) contract, four
  LOFO conditioning fits, one full-data diagnostic fit, prior sensitivity,
  sampler diagnostics, deterministic inference manifests — owns NO final
  predictive scores (Round-4 S-04 de-circularization). (PRE-SPLIT.)
- **E10.2b** Generic distributed exact H-campaign EXECUTOR: consume
  `HistoricalCampaignIntentManifestV1`, run deterministic resumable
  member×flight shards, merge in member-index order, emit
  `HistoricalCampaignReceiptV1` — owns no Bayesian fitting and no cross-stage
  orchestration (surrogates allocate, never replace). (PRE-SPLIT.)
- **E10.2c** H-07 stage ORCHESTRATOR + scorer: exact prior-predictive
  campaign → four LOFO conditioning fits → four exact held-out predictive
  campaigns → full-data posterior-predictive diagnostic campaign → proper
  scores, joint regions, prior-sensitivity summaries, final signed browser
  artifacts. Depends E10.2a, E10.2b. (PRE-SPLIT.)
- **E10.2** Aggregate every case in the FROZEN V/H registry + the E10.2a/b/c
  receipts into fs-vvreg/vv-scorecard (component V-cases execute incrementally
  with E3/E4; E10.2 is aggregation, not first execution).
  **MILESTONE: EVIDENCE-BADGED BETA.**
- **E10.3** Tier D fs-lbm wind-over-terrain runs (Linux perf hosts) → V-04c.
- **E10.4** Cinematic export path (reuses h7xu5 machinery); one hero clip.
- **E10.5** fs-contact certified-contact replay pass receipt.

### Critical path

E0.9 + E1.2/E1.4/E1.5/E1.7/E1.8/E1.9 → E4.0/E4.1 → (E4.2/E4.2b + E3.2a/E3.2b)
→ E4.3 + E4.4a + E4.5 → **E4.6a source-matched FIXED-control open-loop gate**
→ E5.0/E5.1/E5.2 (physics spine running with basic deterministic
record/replay; the archival scrubber + old-artifact retrieval land in E6.1) →
**E4.6b free-control/set-valued gates** → E5.3a (first human-flyable) →
E4.6b0/E4.6c/E4.6d → E5.3b (historically modeled flyable) → E6/E7/E10
aggregation → E10.2 (evidence-badged beta). E4.3b1/E4.3b2/E4.9b/E4.3b3 run in
PARALLEL with the fixed-control spine and gate E5.3b, not E4.6a or E5.3a
(Round-4 P-05). Assets (E2) and terrain (E1.3)
parallel the spine. Wake (E4.7/E4.7b) and field viz (E7) gate the wow, not the
flyable. E0.6/E0.7/E0.8 precede physics tuning; E1.7's freeze precedes any
calibration; E0.5 stays off the critical path. Deferred v1.5+ backlog epics
(explicit, not half-inside v1): QuasiStaticBeamAndRigging, SmoothedTangentPlane
(E4.4b), HeightfieldBoundary,
PilotAdaptiveEstimated, non-neutral atmosphere classes, E7.6 WebGPU compute,
fs-phs/fs-aeroac-grounded audio.

---

## 13. Risks & Mitigations

| # | Risk | L | Mitigation |
|---|---|---|---|
| 1 | Open-loop model misses the A4 pole structure | med | V-02a gates E4.6a before any pilot/feel work; E4.9a early referee fixtures |
| 2 | Canard mechanical model under-sourced (E1.5 thin) | med | qualitative envelope gates (sign/tendency) + explicit uncertainty; V-02b bands widen honestly rather than fake precision |
| 3 | Wake/field costs blow budget on low-end devices | med | tier ladder + QoS governor + induced-error-driven population control; Tier A is the contract |
| 4 | SAB/COOP/COEP hosting constraints | high (known) | dual artifacts + multi-worker transferable fallback from day 1 (E0.7) |
| 5 | 3-D model licensing | med | E2.1 blocking vetting; commissioned fallback |
| 6 | Historical numbers contested | med | [V] discipline, WindReference provenance, ensembles-not-points, E1.7 registry freeze |
| 7 | wasm/native numerical divergence | med | det:: doctrine + E3.5 subsystem digests + four-quadrant+wasm goldens early |
| 8 | Scope creep toward general flight sim | high | §1.4 non-goals; new-aircraft asks become v2 beads |
| 9 | Round-1 scope additions (mech controls, rotor dynamics, hybrid wake, ensembles) slip v1 | med-high | staged milestones: physics-spine flyable precedes historically-modeled flyable; structural beam model and SmoothedTangentPlane are explicitly droppable to v1.5 without breaking any shipped claim (claims are mode-gated) |
| 10 | Instability frustrates casual users | high | journey order (§2.1): ride-along → assist → authentic; "why it porpoises" turns failure into the lesson |
| 11 | Real-time contact jitter on landings | med | 240 Hz substep + regularized friction + E10.5 certified-contact backstop |
| 12 | Pilot-model identifiability (historical data thin) | med | Round-2 split: V-02c1 validates the generic mechanism on independent human-control data; H-02c stays a compatibility check that cannot promote; instrumented-replica campaign (A7) is the only path to more |
| 13 | Correction-table misuse outside domains | low | §4.2/§10.6 rules; hostile twin for stale/out-of-domain table application |
| 14 | Artifact-archive scope (signed registry + CI backward-playback) becomes its own program | med | reuse workspace identity/manifest machinery; ratified R2+S3 pair with verifier-mediated dual publication; the archive is REQUIRED for the old-exact-playback claim — weakening that claim is an explicit decision, never silent |
| 15 | Two-way prop coupling iterations blow the tick budget on low-end devices | med | bounded Aitken with candidate-family selection + reported residual; admission-selected reduced schedules gated by V-15/V-16a |

---

## 14. Round-5 Verification Charter

Rounds 1–4 answered 49 posed questions; the architecture and all closure
contracts are integrated. Round 5 is NOT a design round — its charter, per the
Round-4 verdict, is narrow integration verification:

1. **Diff-integration check:** confirm every Round-4 S-01…S-06 and P-01…P-12
   revision was integrated consistently (the review's exact diffs vs this
   text).
2. **Mechanical name scan:** zero undefined case/task/mode names; zero stale
   pre-split identifiers; zero duplicate contract formulations (the same
   contract stated twice in different words).
3. **Preimage audit:** every identity definition (`RunIdentityBasisV1`,
   `RunIntentId`, `InputTraceV1`, `RunId`, `FieldSourceSnapshotId`,
   `FieldQueryId`, `RunSpecId`, campaign manifests) has exactly one closed
   formula with no self-reference and no partial-field alias.
4. **DAG check:** the task graph is acyclic, every dependency names an
   existing task, and the A1/referee lane (E4.3b1/2, E4.9b, E4.3b3) provably
   gates only E5.3b.
5. **Verdict:** "BEADS READY" or "BEADS READY AFTER LISTED POLISH" is expected
   if 1–4 pass; any genuinely structural finding forces Round 6 with the same
   charter.

Upon a beads-ready verdict: execute E0.1 (epic/leaf conversion under the
conversion gate), then begin E0/E1 leaves.

---

## 15. Appendices

### A. Model equations (implementation-normative)

**A.1 Strip force build-up.** Per strip: local flow = atmosphere + wake induction +
images − body-point velocity; α_eff in the section frame incl. warp twist and
surface deflection; circulatory forces via the owned unsteady mode; profile drag
from fs-airfoil; induced effects via the coupled planform solve — the induced-α
bookkeeping is documented in fs-wing's CONTRACT with a double-count falsifier test.

**A.2 Unsteady kernels (Round-2).** Rational approximations to Wagner Φ(s) and
Küssner ψ(s), order selected by V-08a on BOTH lift and moment channels (Jones'
2-state Φ as candidate); reduced time from ds/dt = 2·U_conv/c with U_conv the
positive CHORDWISE relative-flow component, advanced by exact matrix-exponential
transition over midpoint Δs; states freeze at U_conv = 0; reversed/cross flow
refuses the indicial owner. Apparent mass lives exclusively in the generalized
added-mass operator (§5.1.2). Constants dimensionless and cited (Fung; Leishman;
exact-Wagner references).

**A.3 Coupled planform solve (Round-2).** Warm-started nonlinear lifting-surface
solve with section closure; ≥2 chordwise rows on wing and canard where hinge
moments or unsteady phase are claimed; influence from bound + trailing systems of
every surface + all flat-plane images; safeguarded Picard/Newton with
deterministic continuation and branch identity; typed refusal on ambiguity;
condition estimate reported; factorization reuse only while the complete
influence operator is unchanged. `WeissingerLLinear` = fixture/admission-selected debug mode;
decambering = separate Estimated mode.

**A.4 BEMT (Round-3 S-20 — appendix aligned to the accepted contract).**
Per-annulus momentum/blade-element solve with Prandtl tip/root loss and an
explicitly identified high-loading/static closure. The historical RUNTIME path
is the warm-started bounded station solve with convergence receipt. The
optional design-commit map uses `PropMapCoordinates { J,
radial_reynolds_descriptor, inflow_harmonics[0..m], blade_roughness_id }` with
the smallest dimension/harmonic order passing the paired station-solve
battery; a J-only map is legal only where that battery passes; all maps
include J=0 fixtures. Non-convergence refuses.

**A.5 Wall-compatible turbulence (Round-2).** u = ∇×A with mirrored/parity modes
(u_z(0) = 0 identically); spectral tensor fitted to a pinned Mann-class
neutral-surface-layer target (component spectra, cross-spectra, coherence,
phase, stress ratios); phases advect deterministically with U_adv; amplitudes
follow exact-discrete OU with counter-addressed innovations and SEQUENTIAL
checkpointed state; exact derivatives term-wise. Mean flow: FlatSiteLogLaw for
historical 1903 (one effective z₀); fetch adjustment only via a mass-consistent
mode. Recurrence battery over the longest scenario duration.

**A.6 Ground images (Round-2).** FlatPlaneVortexImageExact: reflection through
the fixed plane with ω′ = det(R)·R·ω; image identities stable; images excluded
from physical wake ledgers; claim scope = the represented singularity field's
slip-wall condition only; FlatnessCertificate gates badges. SmoothedTangentPlane:
one filtered global plane; receipts carry origin/normal, translational and
angular velocity, filter + hysteresis state, τ, slope, max wall residual, and
the artificial generalized boundary-power residual; never "exact."

**A.7 Longitudinal analysis.** Full augmented linearization (rigid + aero-memory +
actuator + rotor [+ pilot]) owns stability claims; the (u, w, q, θ) projection is a
labeled teaching view; x_np from ∂Cm/∂CL root-solve at valid trim only.

### B. Performance math (Round-2 corrected)

Wake budgets are per representation (near filaments / mid elements / far
multipoles), controlled by an induced-error budget, never one global cap.
Example arithmetic at 36 strips, 40 Hz wake, near+mid ≈ 2–3k elements with
binned evaluation ≈ 190–260 MFLOP/s scalar — plausible against the 8.33 ms tick,
with E0.6 measurements as the only acceptance evidence. Tier B average =
3.23 ms/60 Hz frame; wake-update frame = 4.26 ms; per-tick view: 0.58 /
3.68 / mean ≈ 1.61 ms (all kernel subtotals, §7.2). Field budgets are quoted
only with exact interaction counts AND bytes (≈3 MB/snapshot at 32³; leased
ring ≈ 9–12 MB; copy/upload traffic measured).

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
Porpoising: the closed-loop pilot–aircraft pitch oscillation. Warp: roll control
by wing twist. Seqlock: sequence-number ring protocol for torn-read-free
snapshots. U_conv: positive chordwise relative-flow component (the indicial
clock). H-case: historical-compatibility check that cannot promote evidence.
Leased ring: FREE→WRITING→PUBLISHED→LEASED slot protocol. RunIntentId:
provisional identity of an interactive run before its input trace closes.
FlatnessCertificate: the residual/slope/clearance record gating flat-ground
claims. Mann tensor: neutral-surface-layer spectral tensor target (spectra +
coherence + phase). AR_plane vs AR_system: b²/S_one_plane vs b²/S_both — never
bare "AR" for a biplane.

---

## Review round log

| Round | Reviewer | Date | Disposition |
|---|---|---|---|
| 0 | NobleLion / Claude | 2026-08-16 | initial comprehensive plan |
| 0.5 | self-audit (fresh eyes + executed wasm32 probes) | 2026-08-16 | corrected fs-contact role, fs-simd wasm claim, fs-time integrator choice; added §11.3 audit evidence, §11.4 optional-fidelity doctrine |
| 1 | GPT Pro Extended Reasoning (external) | 2026-08-16 | **architecture accepted; major physics-and-validation revision required and integrated**: longitudinal contract split (open-loop / canard mechanics / closed-loop), effect-ownership graph, coupled multisurface Tier A, flat-plane-exact ground images, generalized added mass, rotor dynamics + CT/CQ maps, wall-compatible atmosphere + ensemble historical presets, hybrid near/mid/far wake, unilateral rail release, structural claims mode-gated, identity quintuple, referee-plane rename, validation rebuilt on identifiability/holdouts/pre-registration, corrected budget arithmetic + perf acceptance contract, worker ABI/QoS protocol, fs-airfoil L2 crate, product-copy neutrality, "why it porpoises" flagship view, BPF audio fix, historical-claims table revisions |
| 2 | GPT Pro Extended Reasoning (external) | 2026-08-16 | **"NOT BEADS READY" verdict integrated**: AeroEffectOwners record + chordwise reduced time (the |U| clock bug); prescribed-wake Tier-A candidate A1 + V-08b selection; hinge-load ownership moved to fs-wing + branch/set-valued free-control stability + ModelSafetyLimits split; cue-based two-axis pilot + perception model + InputTransducerMode; AddedMassMode ladder + energy-consistent bias + HeldOnRailEquilibrated prelaunch; two-way prop coupling + warm-started station solve default + J≈0.7–0.8 Dec-17 correction + gyroscopic moments; FlatSiteLogLaw solenoidality fix + Mann-class tensor target + exact-discrete OU + air-state API; 120 Hz shedding + core-evolution modes + moment-complete conversions + mixed-norm errors; FlatPlaneVortexImageExact rename + FlatnessCertificate + swept contact proxies; ReducedAeroelasticWarp; corrected 3.23/4.26 budget + kernel-subtotal honesty + tick view + field memory; service/lateness/backlog metrics + vsync-relative render gates; leased ≥3-slot snapshot ring + non-SAB pack/copy; input-tick sync + ApplyNextEligibleTickAndFlag; hysteretic QoS machine; evidence lineage/independence groups; V/H-case split + hierarchical H-07 + unstable-system convergence protocol + anti-vacuity baselines; RunIntentId + derived weather realization + CheckpointStateV1 + artifact archive; joint historical parameter model + attribution honesty; copy fixes; "why it rolls and yaws"; historical flags (30-ft derrick, 20–21 ft² rudder, AR conventions, anhedral intent, engine-mass boundary) |
| 3 | GPT Pro Extended Reasoning (external) | 2026-08-16 | **"ROUND 4 REQUIRED" — 23 structural + 16 polish items, ALL integrated**: joint (6+nc) added-mass solve; fs-flyer-owned prop coupling (L3 siblings decoupled) + Aitken scheme + immutable schedules; stationary-draw prelaunch initialization + E4.6d closure task; PhysicalRealizationAlgorithmId + RunAnchor + applied-trace-only InputTraceId; checkpoint algorithmic history; AtmosphereStateSnapshot in field snapshots; FOM/ROM split with shared-basis rational-Krylov + two receipts; E3.3a/b/c atmosphere workstream + full-tensor fitting; fixed-control-first gate ordering + ControlConstraintCase; E4.4a/E4.7b/E1.10 dependency closure; mandatory LOFO H-07 + offline full-Bayes + campaign manifest + no-surrogate-finals; A7a/A7b split (GO, Wright Experience candidate); persistent tracer service replacing the dense history ring; hard lateness gates + V-18 field-noninterference; runtime mode-immutability triad; phase-resolved blade capsules + TerminalEvent semantics; archive hosting/signing/sandbox/retention contract; deterministic perception service (V-16b); SmoothedTangentPlane boundary-work accounting; A.4 BEMT appendix aligned; WebGPU strictly v1.5 (E7.6); epic/leaf beads rule with pre-split list; 16 polish fixes (copy neutrality, stale refs, milestone renames, status labels, six-lane goldens naming, deferred-mode backlog) |
| 4 | GPT Pro Extended Reasoning (external closure audit) | 2026-08-17 | **"ROUND 5 REQUIRED" — six structural + twelve polish groups + seven decisions, ALL integrated**: closed identity preimages (RunIdentityBasisV1, two-input RunId, trace-extent-bearing InputTraceV1, execution-closure ArtifactId, RNG realization-rooting, RunSpecId cache); state-bound field identities (FieldSourceSnapshotId, snapshot-lease API); branch-aware CRN-preserving prelaunch (StationaryOuPathV1 −32 s anchor, PhysicalInitializationSpec + ControlEquilibriumSelection, partitioned alternate starts, numerical closure bands); de-circularized H-07 (intent-manifest/receipt split, per-member allocations, E10.2a/b/c); A1 image/referee ownership (common image kernel in E4.4a, image-aware FOM/ROM, independent E4.9b referee); terminal-safe A/B (ABComparisonReceiptV1, user-intent counterfactual channel); Aitken candidate family; 16→24 blade capsules with generated radii; 120 Hz perception baseline; ratified R2+ObjectLock archive with verifier-mediated dual publication + trusted-loader boundary; prelaunch policy numbers; conversion gate + expanded pre-split; V-05d/V-11c/V-21 receipts, V-20 expansion, A7a promotion path; seqlock run-epoch; pilot-uncertainty allocation; stale-name sweep; budget-row honesty; schedule inventory; tangent-plane v1.5 deferral; risk reorder |
| 5+ | — | — | pending; narrow verification charter (§14); convert to beads at steady state |
