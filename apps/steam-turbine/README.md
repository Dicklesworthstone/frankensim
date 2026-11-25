# Industrial Steam Turbine Simulation (`apps/steam-turbine`)

Interactive WebAssembly + HDR 3D Three.js + 48 kHz WebAudio multiphysics simulation and L4 Design Studio for a **50 MW industrial extraction-condensing steam turbine**, built on **26 FrankenSim crates** across the convergence architecture (`fs-ports` + `fs-system`), L2 MORPH CAD geometry (`fs-rep-frep`, `fs-rep-nurbs`, `fs-rep-mesh`, `fs-topo`), L3 FLUX physics (`fs-thermochem`, `fs-multiphase`, `fs-mbd`, `fs-tribo`, `fs-phs`, `fs-control`, `fs-plate`, `fs-duct`, `fs-modal`, `fs-aeroac`, `fs-psycho`), and L4 ASCENT derivative-free optimization (`fs-dfo`).

The machine is a generic 50 MW extraction-condensing unit: its flow path, masses and ratings are this model's own, not those of any manufacturer's product.

Time runs at physical speed. The metal and oil thermal nodes carry the full `m·cp` of the parts the geometry builds (about 1.1 MJ/K for the HP rotor and 5.3 MJ/K for the IP casing), and the 2-D field steps with the same clock, so a load change moves the HP rotor by tens of kelvin over minutes, not seconds; use the fast-forward buttons to watch it. One simulated second costs 150 to 160 ms of native release CPU (nearly all in `fs_system::System::step`, whose steam components are differenced, having no declared tangent), which is why the page offers at most 3×.

The condenser shell is a saturated control volume: vapour of mass `V·ρ_g(p)` over its condensate, with pressure solved each step from the energy balance `d/dt[V·ρ_g·(u_g − h_f)] = H_in − ṁ_in·h_f − q`. Wet exhaust therefore leaves its liquid as condensate, and at steady state the duty equals `ṁ_in·(h_in − h_f)` exactly. Entropy production of the valves, stage sections and condenser is reported as computed, never clamped at zero.

The two lumped stage sections of the plant graph (HP: Curtis wheel + 5 stages; IP/LP: 11 stages) carry no efficiency of their own. After every 20 ms sub-step the 17-stage row-by-row meanline is solved at the plant's state (29 µs natively), and each section's isentropic efficiency is read off it: HP from the chest to the extraction pressure, IP/LP total-to-static to the condenser, so the last stage's leaving loss counts. Blade heights, chords, shrouds, tip clearance, velocity ratio and wetness therefore reach the shaft power through one loss model, and the per-stage chart adds up to the plant's section power (HP exactly, after labyrinth leakage). The loss model itself (`stage_meanline.rs`) still starts from per-family base efficiencies of 0.835, 0.905 and 0.920; those are its fitted inputs.

The generator is the classical machine on an infinite bus: a constant EMF `E′` behind the transient reactance `x′_d = 0.25 pu` plus a step-up transformer's `x_t = 0.10 pu` (60 MVA, 11 kV), stator transients neglected. Each step solves the stator voltage equations against the bus at the load angle and takes torque, currents and copper loss from `fs_phs::circuit::synchronous::SynchronousMotor::operating_point`, so the excitation slider moves reactive power (about −45 MVAr at 0.8 pu, +45 MVAr at 1.3 pu) while real power stays with the steam. The rotor swings against the bus at `√(K_s/J)` (about 2.3 Hz bare, 2.4–2.6 Hz in the plant, where the governor's droop acting through the steam-chest lag adds stiffness), damped by an amortisseur coefficient of 15 pu, the one lumped machine constant. A pole slip (`|δ| > π`) trips the generator off the bus; open the breaker to clear it, then close it to resynchronize.

## Quick Start

```bash
cd apps/steam-turbine
npm install
npm run dev        # local Vite server with COOP/COEP headers
npm run test       # Node WASM integration & physical emergence tests
npm run build      # TypeScript check + production bundle
npm run wasm       # rebuild release wasm32-unknown-unknown + wasm-bindgen
```

Native Rust verification:

```bash
CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_turbine" \
  cargo ltest --offline --manifest-path apps/steam-turbine/engine-turbine/Cargo.toml -j 4
```

## Origin of Each Subsystem

Per FrankenSim Operating Rule 2, every physical and visual subsystem is classified explicitly by its computational origin:

| Subsystem | Origin | FrankenSim Crates / Implementation |
| --- | --- | --- |
| **Multiphysics Coupling, Newton Loop Solve & Conservation Audit** | `core solve` | `fs-ports::Component` + `PortTangent`, `fs-system::SystemBuilder` / `System::step()` / `System::ledger()` (12 coupled component instances audited every step for First-Law energy, stream mass, carried H₂O, and Second-Law entropy production) |
| **Steam Properties & 17-Stage Row-by-Row Meanline** | `core solve` + `in-app model` | `fs-thermochem::water::if97` (`If97State`, `If97SaturationPair`, `If97Region`), exact safeguarded `(P,h)`, `(P,s)`, and `(u,v)` inversions; 17-stage Euler velocity triangles, Ainley-Mathieson / tip-clearance / Baumann losses, blade-root centrifugal + bending σ_vM & Goodman SF |
| **2D Axisymmetric Thermal & Thermo-Elastic Stress PDE** | `in-app model` | 64×24 cylindrical finite-volume transient heat conduction PDE with `matdb` temperature-dependent `k(T)`, exact Fourier heat flux `q = -k(T) ∇T`, rotating-disc centrifugal + Lamé pressure hoop + thermal-gradient σ_vM, exported as live GPU `DataTexture`s + bilinear surface probe |
| **L4 ASCENT Multi-Variable Machine Design Optimizer** | `core solve` | `fs-dfo::cma::cmaes` (Natural-Gradient CMA-ES) + `fs-dfo::nelder_mead` simplex polish optimizing 6 physical hardware dimensions (tip clearance, HP reaction, LP Stage 16 blade height, bearing clearance & preload, HP casing wall) live in WASM under structural, wetness & 4-DOF stability constraints |
| **Exhaust Hood Structural FEM Modes & Duct Radiation** | `core solve` | `fs-plate::shell::{generate_bell_shell, assemble_shell, modes_shell}` (Sylvester inertia certified S355 conical-to-cylindrical shell eigenmodes) + `fs-duct::levine_schwinger_impedance` (exact Wiener-Hopf open-end radiation impedance) |
| **Condenser Film Condensation & Cooling-Water Heat Transfer** | `core solve` | `fs-multiphase::condensation::CondensationModelV1` (Nusselt + Colburn-Hougen non-condensable air mass-transfer resistance) + `fs-convection::evaluate(CorrelationId::Gnielinski)` |
| **4-DOF Gyroscopic Rotordynamics, Campbell & Stability Poles** | `core solve` | `fs-mbd::rotordynamics::RotorBearingModel` (`at_speed`, `poles`, `campbell_samples`), unconditionally A-stable 4×4 implicit Newmark-β integration with Lund/Childs flexible-shaft impedance reduction |
| **Hydrodynamic Tilting-Pad Bearings & Thrust** | `core solve` | `fs-tribo::journal::JournalGeometry` + Ocvirk short-bearing Sommerfeld `(K_ij, C_ij)` + Shell S2 ISO VG 32 temperature-viscosity curves from `data/matdb/seed-v1` |
| **Synchronous Generator & Core/Winding Losses** | `core solve` + `in-app model` | `fs-phs::circuit::synchronous::SynchronousMotor::operating_point` (torque, d-q currents, flux linkage of the classical `E′`-behind-`x′_d` machine on an infinite bus) + reactive power, swing dynamics and out-of-step trip + M19 Steinmetz core loss + IACS copper `I²R` losses from `matdb` |
| **Governor Speed Control & Thermal Derating Protection** | `core solve` | `fs-control::speed::SpeedController` + `fs-control::thermal::ThermalDeratingPolicy` (4% ISO droop, 110% overspeed ESV latch, bearing/hood thermal derating & trip) |
| **ISO 7919 Anytime-Valid Vibration Monitoring** | `core solve` | `fs-eproc::BettingEProcess` sequential log e-value |
| **L2 MORPH CAD Geometry (Rotor, Labyrinth, 5 Cutaway Housings, 7 Airfoil Families)** | `core solve` | `fs-rep-frep::AxisymmetricChart` (`tessellate_axisymmetric`, Green's-theorem `mass_properties`, `surface_area`), `fs-rep-nurbs::NurbsSurface` (`partials` analytic normals), UV parameterization & exterior/cap material split, certified via `fs-rep-mesh::Soup` + `fs-topo::manifold_certificate` (watertight rotor, Euler χ = 2, 0 defects) |
| **Material Properties & Blackbody Radiometry** | `core solve` | 21 compile-time embedded `data/matdb/seed-v1` packs (`fs-material`) + `fs-render::thermal_surface::GrayDiffuseThermalSurface` + CIE 1931 Planck blackbody LUT |
| **48 kHz Turbine Acoustics & ISO 532-1 Psychoacoustics** | `core solve` + `in-app model` | `fs-modal::time::ModalAcousticTimeModel` (seeded by `fs-plate` FEM shell modes) + `fs-aeroac::curle3d::CompactForceObserver` + `fs-psycho` 28-band third-octave `loudness_stationary` / `sharpness_din` |
| **HDR Studio Rendering, Procedural Machining Maps, 6 GPU Field Shaders & 6 Vector Fields** | `display only` | HalfFloat MSAA `EffectComposer` + `UnrealBloomPass` + `SMAAPass` + `NeutralToneMapping`, procedural CNC lathe/cast/milled normal maps, `onBeforeCompile` GPU per-pixel continuum colormaps + iso-contours, 6 GPU-instanced tapered-arrow vector fields, interactive 3D surface probe, and full desktop mouse/trackpad + mobile touch CAD navigation |
| **Guided Experiment Presets** | `scripted` (inputs only) | Presets inject physical boundary conditions or faults into `fs-system`; all responses are emergent |

## Not Modelled

- Full 3D blade-to-blade compressible Navier-Stokes CFD (RANS/LES) or 3D non-axisymmetric FEM casing contact (the steam flowpath uses 17-stage row-by-row IAPWS-IF97 meanline + Stodola + Baumann models, and metal fields use a 2D axisymmetric 64×24 FVM PDE).
- Two-phase non-equilibrium Wilson-line nucleation kinetics (wetness uses equilibrium IAPWS-IF97 Region 4 quality).
- Sub-transient / negative-sequence 3-phase short-circuit fault currents inside generator stator slots or external boiler drum water-level dynamics.
- Blade Campbell flutter / high-cycle fatigue crack propagation across individual blade root fir-trees.
- Full acoustic room-reverberation BEM (free-field Farassat/Curle dipoles + FEM modal radiator + Levine-Schwinger open-flange load at a movable 3D virtual microphone).
