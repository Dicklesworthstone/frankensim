# Chronomètre Souverain · Calibre FS-1903 Tourbillon (`apps/pocket-watch`)

This is an interactive 3D and WebAssembly demo of a one-minute tourbillon pocket watch. It draws its values from
sourced material cards in `data/matdb/seed-v1`. You can heat the watch with a torch, pour liquid nitrogen on it,
flood it with oil, bring a magnet near it, shock it and listen to it through a movable microphone.

**Where the behaviour comes from.** The demo's own Rust engine, [`engine-watch`](engine-watch/src), advances every
subsystem below. FrankenSim crates supply:
- leaf values: material properties, friction laws, convection correlations, and duct and plate-layer terms;
- the escapement's contact solve on each step;
- the loudness filter tables;
- the rate e-value.

No subsystem is stepped by a FrankenSim solver or by the composition kernel. Read the demo as real physics laws
coupled inside the app, not as a FrankenSim assembly.

## Quick Start

```bash
cd apps/pocket-watch
npm run wasm      # Rebuilds the Rust -> WASM engine (engine-watch)
npm test          # Runs the Node WASM integration tests
npm run build     # Typechecks (tsc --noEmit) and builds the production Vite bundle
npm run dev       # Launches the interactive 3D + 48 kHz WebAudio laboratory
```

## Origin of each subsystem

The table uses three classes:
- **core solve**: a FrankenSim crate advances or solves the subsystem.
- **in-app model**: `engine-watch`'s own equations advance it, calling the crates listed for leaf values.
- **scripted**: fixed constants, thresholds or waveforms set it.

All files are in [`engine-watch/src`](engine-watch/src).

| Subsystem | Origin | Where | FrankenSim crates it calls |
|---|---|---|---|
| Mainspring torque and let-down | in-app model | `horology.rs`: `torque_at`, `stiffness`, `step` | — |
| Winding ratchet and click | scripted (30 teeth, fixed click energy) | `horology.rs`: `wind_crown` | — |
| Going train, tourbillon cage, moonphase | in-app kinematics: angles follow the escape wheel, and mesh efficiency uses a fixed coefficient | `horology.rs`: `mesh_efficiency`, `train_efficiency`, `step` | — |
| Swiss lever escapement: lock, impulse, drop | in-app model with a FrankenSim contact solve on each step; the restitutions (0.06 to 0.30) and the three pallet-friction levels are fixed | `escapement.rs`: `fine_step`, `detect`; `horology.rs`: `assemble_loads` | `fs_contact` (`solve_normal_impacts`, `solve_isotropic_coulomb_impulse`), `fs_tribo` (`FrictionLaw::Stribeck`) |
| Balance and hairspring | in-app model: an integrated oscillator. The hairspring length is back-solved so that the movement runs at the nominal 2.5 Hz | `horology.rs`: `hairspring_k`, `balance_inertia`, `free_length`, `regulate`, `REGULATION_FACTOR` | — |
| Rate, beat error, COSC verdict | in-app measurement from the lock drops; the e-value is FrankenSim's | `horology.rs`: `fit_rate` | `fs_eproc` (`BettingEProcess`) |
| Heat flow: torch, liquid nitrogen, conduction | in-app model: explicit enthalpy finite volumes on a 64 × 64 in-plane grid (`GRID_N`); the 12 "nodes" (`NODE_COUNT`) are averages for display | `multiphysics.rs` | `fs_convection`, `fs_material` (air state, sublimation) |
| Torch footprint | scripted: a fixed film coefficient and flame temperature over a 2D footprint | `multiphysics.rs` | — |
| Liquid-nitrogen boiling | in-app model: the Rohsenow, Zuber and Berenson pool-boiling curve | `multiphysics.rs`: `ln2_boiling_flux` | — |
| Rim expansion, hairspring modulus E(T) | in-app model with handbook coefficients | `horology.rs`: `balance_inertia`, `hairspring_k`; `materials.rs` | — |
| Oxide films | in-app growth law (parabolic Arrhenius); colour from FrankenSim's thin-film optics | `materials.rs`, `multiphysics.rs` | `fs_render` |
| Incandescence | display only; it takes no part in the heat balance | `materials.rs`, `multiphysics.rs` | `fs_render` |
| Gold case melting and flow | in-app model: enthalpy with the 18k alloy's solidus and liquidus (1158.15 and 1168.15 K), and shallow-layer flow by virtual pipes. The thresholds for releasing the wall and for mobility are fixed | `multiphysics.rs`, `materials.rs` | — |
| Oil viscosity | in-app model: a Vogel–Fulcher–Tammann law fitted to the material rows | `materials.rs` | — |
| Oil drag at pivots and the train | in-app model (Coulomb plus Petroff) | `horology.rs` | `fs_tribo` (`JournalGeometry`) |
| Oil bridges between hairspring coils | scripted stiffening factor | `horology.rs` | — |
| Magnet: field, eddy braking, lever torque, remanence | in-app model: a dipole field and a play-operator remanence. The hairspring attraction is a fixed coefficient (`magnetic_stiffening_at_sat`), and degaussing is a scripted 50 Hz decaying waveform | `multiphysics.rs`, `horology.rs`, `materials.rs` | — |
| Shock | in-app model: a half-sine load integrated in the movement; its direction and the case spin are scripted | `horology.rs`: `apply_mechanical_shock`, `add_shock_loads` | — |
| Beat sound and microphone | timing from the escapement's contact events. The 12-mode resonator bank has fixed frequencies and decays, and the microphone is a clamped 1/r gain | `acoustics.rs` | `fs_material`, `fs_vkplate`, `fs_duct` (small tuning terms), `fs_modal` (a parallel modal bank), `fs_aeroac` (torch noise) |
| Loudness meter | core solve: FrankenSim's filter tables | `acoustics.rs` | `fs_psycho` |

## Not modelled

- The motion works (12:1 hour train) and the power-reserve differential. The moonphase angle comes from time alone.
- A boundary-to-elastohydrodynamic friction transition. Pallet friction is the three-level table above.
- The oil's surface tension in its flow. `surface_tension` exists for the oil, but the flow uses a constant.
- Heat from the movement's friction. The thermal field does not receive it.
