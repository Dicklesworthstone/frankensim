# First Flight (`apps/wright-flyer`)

A browser front end (TypeScript, three.js) for the Wright Flyer simulation in
[`crates/fs-flyer`](../../crates/fs-flyer). [`crates/fs-flyer-wasm`](../../crates/fs-flyer-wasm) compiles it to
WebAssembly, and a Worker steps it (`src/sim/simWorker.ts`).

## Running it

```bash
cd apps/wright-flyer
npm run wasm       # builds crates/fs-flyer-wasm into src/wasm-pkg (not committed)
npm run dev        # Vite dev server
npm run build      # tsc --noEmit, then the production bundle
npm run test:sim   # the simulation facade tests
```

`src/wasm-pkg` is a local build and is not in git. If it was built before the current `fs-flyer-wasm`, it lacks
exports that the Worker calls, such as `flyer_engine_checkpoint`. Rebuild it before reading the page as the code on
`main`.

## What the page runs

- **Without `?sim=1`**, the page plays a 12 s keyframe loop (`scriptedState` in `src/airframe/applyPose.ts`, and
  `src/flyerScene.ts`). No simulation runs.
- **With `?sim=1`**, the Worker steps fs-flyer's `SimLoop`. The menu, the challenges and the guided journey all link
  here. `src/menu.ts` `scenarioFromQuery` builds the scenario from the URL:
  - seed 1903, or a per-mission seed;
  - an 11 m/s headwind;
  - air density 1.294 kg/m³;
  - an 18.3 m rail.

  No page scenario sets a spatial wind or a material air model.

## Origin of each part

The classes match the other demos:
- **core solve**: a FrankenSim crate advances the part.
- **in-app model**: the page's own TypeScript computes it.
- **scripted**: fixed constants, keyframes, thresholds or placeholder data set it.
- **display only**: it takes no part in the dynamics.

| Part | Origin | Where | FrankenSim crates |
|---|---|---|---|
| Longitudinal flight (distance, height, pitch) | core solve: a rigid body advanced by a midpoint step | `fs-flyer` `simloop/motion.rs` (`midpoint_stage`) | fs-flyer |
| Roll and yaw | core solve, reduced: a roll and yaw response. There is no sideslip and no lateral translation | `fs-flyer` `lateral.rs`, `simloop.rs` | fs-flyer |
| Wing and canard aerodynamics | core solve: a nonlinear vortex lattice, closed per section by a thin-airfoil and flat-plate blend | `fs-flyer` `aircraft.rs` (`InfluenceOperator::build`), `sections.rs` (`section_lift`) | fs-wing |
| Rudder | scripted: a moment of 66 N·m/rad on the warp command, declared as an estimate | `fs-flyer` `simloop.rs` (`WRIGHT_1903_RUDDER_GAIN_NM_PER_RAD`) | — |
| Propellers | core solve: blade-element momentum, coupled to the wing through the slipstream | `fs-flyer` `propcoupling.rs` | fs-airscrew, fs-wing |
| Engine | scripted rating: 8948 W at 1025 rpm. The propeller speed follows from the torque balance | `fs-airscrew` `drive.rs` (`ShaftDrive::wright_reference`) | fs-airscrew |
| Gusts | core: eight Ornstein–Uhlenbeck modes on top of the fixed headwind. Their correlation times are fixed | `fs-flyer` `simloop.rs` | fs-atmo |
| Launch rail | core: the aircraft leaves the rail when the rail's upward force reaches its weight | `fs-flyer` `simloop.rs` | — |
| Landing | threshold: the run ends when height reaches 0 (`GroundContact`). No skid or sand contact model is stepped | `fs-flyer` `simloop.rs` | — |
| Modelled pilot (member 3) | core, with registered fixed gains. Its calibration predates the plant's momentum and gravity corrections (`fs-flyer` CONTRACT.md) | `fs-flyer` `pilot.rs`, `perception.rs` | fs-flyer |
| Wing warping | core, reduced closed form. It drives the roll response, not the vortex lattice | `fs-flyer` `aerowarp.rs`, `simloop.rs` | fs-flyer |
| Canard mechanism | core: one degree of freedom. The aerodynamic hinge moment is set to 0 | `fs-flyer` `canardmech.rs`, `simloop.rs` | fs-flyer |
| Player controls | in-app mapping of keys to lever force (±220 N) and warp (±0.148 rad), with slew limits | `src/input.ts`, `src/sim/humanControls.ts` | — |
| Flight audio | display only: a two-oscillator engine tone driven by propeller speed and airspeed | `src/audio.ts` | — |
| Eigenmodes tile: poles | core solve: fs-flyer's augmented linearization (rigid u, w, q, θ; canard angle and rate; common prop speed; no pilot) at the canonical trim, with structural mode-family attribution, computed once by the wasm export `flyer_modes_json`. The poles belong to the trim, not the live state | `fs-flyer` `augmented.rs` (`canard_control_column`, `build_engine`, `attribute_modes`), `fs-flyer-wasm` `modes_json` | fs-flyer |
| Field instruments: force overlay, tracers, porpoising, the eigenmodes tile's live-state projection, field glyphs, sweeps, lateral view | scripted. Each of these readouts is made up in the page from the live state, and its first line says so ("ILLUSTRATION" or "APPROX"). Example: the glyph field is the aircraft's velocity with a 2 % shear. fs-flyer's field service is not called | `src/instruments.ts` | — |
| Evidence badges, scorecard, applicability | unavailable. No receipt source is connected to the page, so these tiles build no badge, digest or verdict | `src/instruments.ts` (`NO_RECEIPT_SOURCE`) | — |
| Landing score, missions, lessons | scripted thresholds and text. The mission outcomes come from the engine | `src/landingScore.ts`, `src/missions/`, `src/lessons.ts` | — |

## Not modelled

- **Ground effect.** The flight builds its vortex-lattice operator without ground images. Only `fs-flyer`
  `campaign.rs` uses `fs_wing::images::CertifiedGround`.
- **Added mass.** `aircraft.rs` states that the noncirculatory terms are not assembled.
- **Sideslip and banked translation.**
- **Skid contact on landing.**

## The December 17 runs end on the rail

The menu's "1903 PILOT MODEL" and "FIXED CONTROLS" modes, and stage 1 of the guided journey, run seed 1903.
- With fixed controls, or with pilot member 2 or 3, that scenario ends `RailEndWithoutLift` at tick 978.
- The test `original_dec17_reference_reports_the_corrected_rail_shortfall`
  (`crates/fs-flyer/tests/simloop_battery.rs`) pins this.
- The menu and the journey caption say so (`src/menu.ts`, `src/journey.ts`).

## Labels on the page

- **Ground effect.** The ground-effect lesson (`src/lessons.ts`) teaches the certified ground-image model. Its
  first step also says that this page's flight is solved in free air.
- **Sound.** The menu's key card shows `AUDIO_DISCLAIMER_LABEL` (`src/audio.ts`): synthesized sound design, not
  an acoustic simulation.
- **Landing score.** `src/landingScore.ts` grades the sink rate against declared bands, 1.2 and 2.5 m/s, for
  which no source is cited. The sink rate is sin θ·u − cos θ·w from the final frame.
- **Gust seeds.** The flight 1 and flight 4 challenges (`src/challenges.ts`) change only the gust seed. The
  headwind is 11 m/s for every flight.
