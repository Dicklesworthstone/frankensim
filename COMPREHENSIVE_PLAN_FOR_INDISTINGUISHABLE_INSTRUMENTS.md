# Comprehensive Plan: Instruments Indistinguishable from Real Recordings

Status: **v5.1, 2026-10-16 (IndigoLeopard: Beads publication sync, building on MagentaBadger's v5 and the Codex reviews).** v1 came out of the reality check. v2 folded in two reviews (physics, execution). v3 re-centred the plan on generic physics modules after inventorying the workspace at HEAD `47a361a87`. v4 folded in four independent reviews of v3:
- numerical methods and architecture;
- musical-acoustics completeness;
- execution, AGENTS.md fit and genericity;
- a Codex (GPT) full review.

A fifth review, by the Gemini CLI, could not run: the client refused to start ("no longer supported for individuals"). v5 adds the owner's explicit real-geometry/material/anatomy requirement and the results of a completed two-family research duel run by IndigoLeopard: Claude Opus 5.5 and GPT-6.1-Sol researchers, 60 candidates, 10 finalists, dual cross-scoring; Gemini was unavailable. Research findings and review status are recorded below as they are integrated; proposed experiments are not passing evidence. The v4 repository inventory was checked at HEAD `9267fd56d`; unverified historical failures remain dated observations.

v5 integrates:
- the owner's directives of 2026-10-15 (evening):
  - real geometry, materials and interacting anatomy, with simplification only under demonstrated acoustic equivalence (§1.4, drafted by IndigoLeopard; geometry and anatomy work in M6.3/M6.4);
  - three instruments beyond the orchestra: duduk, Irish tin whistle, Andean pan pipes (§6.10–6.12);
- a Codex round-2 review of v4. All 13 structural findings are accepted:
  - gates: transitive gate prerequisites; no accelerations in acceptance; binding spectrogram equivalence; Q5 rows and family coverage; claimable-item readiness;
  - physics: compressible thermal flow; multimodal time-domain ducts; nonlinear edge flow; missing coupling edges; restitution in nonsmooth gen-α;
  - method and rules: per-axis error models; construction/setup/performance separation; band-complete bakes;
- the duel's adjudicated deltas and source leads. The synthesis is `DUELING_WIZARDS_REPORT.md` and `ARTIFACT.md` in `/Users/jemanuel/brennerbot_sessions/frankensim-instruments-20261015/`, and the leads are in its `evidence/`;
- IndigoLeopard's bounded review of the integrated v5. The register has 346 unique work items, with references checked independently.
- **harmonization with origin/main through f9c6f1fcb** (49 commits since the inventory HEAD; 36 not yet in the local tree). The "Now" sections, Appendix B anchors and Appendix E rows now credit work that is partially delivered but natively unverified: ModalAperture multimode reeds, shared-support strings, membranes, lumped shared cavities and gas pistons, physical outlet observers, live-piano oversampling. IRP items verify, ladder and consolidate these owners instead of reimplementing them. Commits pushed after f9c6f1fcb (e.g. 39f15c202, 46a55d5) are re-checked when each bead is claimed (§10.6 rule 3).

**Scope.** Every instrument other than piano and drums:
- winds: clarinet, oboe, saxophone, bassoon, flute, recorder, organ flue pipe;
- brass: trumpet, trombone, horn, tuba;
- bowed strings: violin, viola, cello, double bass;
- guitars: classical and steel-string;
- tuned bars: xylophone, marimba, vibraphone;
- voice;
- instruments beyond the Western orchestra (owner, 2026-10-15):
  - **duduk** (Armenian double reed);
  - **Irish tin whistle** (duct flute);
  - **Andean pan pipes** (siku/zampoña and antara);
  - further candidates the owner selects from §6.13 (OD14).

Each is rendered offline at the highest achievable fidelity and graded against real recordings. It then goes on frankensim.org/music, and is finally played live in the browser on the same physics.

**The organizing idea: an instrument is data, and the physics is a small set of generic, abstract, reusable modules.**
- Every capability in this plan is a module that is just as useful outside music: structural acoustics, aero-acoustics, biomechanics, tribology, HVAC, automotive NVH, machine dynamics.
- Every instrument is a description that composes those modules.
- No instrument-shaped code exists anywhere, except fixtures in tests.

---

## 0. How to read this

- §1 sets the owner standard and the operational rules that follow from it. It overrides every earlier music document where they disagree.
- §2 is where we are: measured sound, the inventory's headline findings, and corrections the reviews made to the inventory.
- §3 is the architecture: stack tiers, how modules map to crates, ports, the integrator doctrine, and module contracts.
- §4 is the core: one specification per module. Each covers current state, target physics, numerics and convergence, interfaces, validation, consumers, and work items.
- §5 covers consolidation, defects, breaches of the standard, and registry obligations.
- §6 composes each instrument family from modules. It includes the first audible milestone and named experiments with their falsifiers.
- §7 is acceptance and grading. §8 is the data program. §9 is runtime, site and live.
- §10 is execution: tracks, week one, swarm rules, build discipline, and per-item acceptance requirements.
- §11 is risks, §12 owner decisions, §13 existing programs and overlapping beads.
- Appendix E is the **authoritative work-item register**: ID, size, dependencies and DONE-WHEN for every bead-sized item. Where Appendix E and the §4 text disagree, Appendix E wins.

**Naming.** These names avoid collisions with note names, workspace layers and existing bead titles.

| Kind | Format | Example |
|---|---|---|
| Stack tiers (not workspace layers) | **S0–S7** | — |
| Modules | **Mx.y** | — |
| Work items (bead titles `[IRP] Mx.y-Wn: …`) | **Mx.y-Wn** | — |
| Quality gates | **Q0–Q5**, **QP** (passive chain) | — |
| Experiments | **X-…** | X-Cl-R2 |
| Milestones | **M-…** | M-CL-Q2 |
| Owner decisions | **OD1–OD16** | — |
| Defects | **DEF-nn** | — |
| Breaches of §1.3 | **BR-nn** | — |
| Musical notes: always the word "note" plus the octave | — | note D3, note A4 = 440 Hz |

Sizes: S ≤ 2 agent-days, M ≤ 1 week, L 2–4 weeks. Nothing in Appendix E is larger than L, and most L items are split.

---

## 1. The owner standard

### 1.1 Statement (owner, 2026-10-15)

> "We are NOT interested in 'reduced' engines, EVER, like even a bit. This project is about THE MOST
> ACCURATE POSSIBLE, even if it can't be done in real-time or even close to it. The true test is that a
> spectrogram of our purely synthetic/simulated instrument is basically indistinguishable to a
> well-recorded sample of a real instrument (assuming the instrument is in good condition, properly tuned,
> properly played, etc.)."

On this plan, the owner asked to close the gaps for all instruments and make them sound perfect, "including adding or improving or wiring in any new physics simulation based modules into FrankenSim, trying to make them as generic, abstract and versatile" as possible.

Standing doctrine (2026-08-08):
- bottom-up physics;
- universal parameterized primitives;
- deep integration, with no bolt-ons;
- total parameterization of materials, geometry and air;
- extreme-regime coupling is in scope long-term;
- no instrument modules.

### 1.2 What counts as a model, operationally

**Rejected:**
- fitted loss, dispersion or reflection filters;
- Karplus–Strong and fitted digital waveguides;
- vector-fit "held note" images used as the instrument;
- single-mode stand-ins with hand-typed constants;
- compact one-mode bodies;
- internal pressure presented as sound;
- authored noise;
- any parameter adjusted to make a recording match;
- any model chosen for speed whose error against the fuller model is not measured;
- **any fixed closure (jet-drive law, contraction coefficient, synthetic noise spectrum) used as the shipped render**. Closures may only accelerate, and must be checked against a resolved reference render (OD11).

**Accepted:** a numerical discretization with an explicit refinement path to the **declared governing physical problem**. Numerical convergence and experimental validation are separate requirements, and both are mandatory.
- Complete-basis modal methods and computed operator approximations qualify only when enrichment recovers the declared problem.
- Every physical mechanism the model omits needs either a recoverable refinement axis or a recorded, explicit limitation. A recorded limitation blocks acceptance of any claim it affects.

1. **Discretization ladder** (verification).
   - **Refinable axes.** Every module declares its refinable axes (modes, mesh h, element order p, BEM panels, Δt/substeps, sample rate, line length, FFT length, rational-approximation pole count, lattice spacing, quadrature order, solver tolerance, ensemble size). Each axis has an **applicable convergence/error model**:
     - an algebraic order (h, Δt, p, lattice spacing, quadrature);
     - a spectral or tail bound (modes, FFT/line length, transverse bases);
     - a certified operator error (pole count, BEM/FMM tolerance);
     - a residual-based bound (algebraic solvers);
     - a statistical uncertainty (ensembles).

     It also declares **coupled axis groups** that must be refined together, e.g. FDTD Δx with Δt at fixed CFL, or LBM Δx with Δt under acoustic scaling, which also moves τ.
   - **Rungs.** Convergence is established on **≥ 3 systematically refined rungs** per axis or group, and more when behaviour is scattered or pre-asymptotic.
     - **Algebraic axes** report observed order against formal order. Richardson/GCI estimates (Roache 1998) are used only where their assumptions are demonstrated; Eça & Hoekstra (2014) separate well-behaved from scattered convergence.
     - **Other axes** report their declared bound.
     - **Axis liveness** is proven from the discretization dimensions and tolerances actually used (logged), plus a sensitivity fixture showing that the axis can change a result. An exact answer that legitimately does not change under refinement is permitted. This catches mode caps that silently clamp and FFT overrides that are ignored.
   - **Joint refinement.** A **jointly refined end-to-end run** exposes interactions between axes, such as contact resolution against Δt against acoustic bandwidth.
   - **Error categories** are reported separately: geometry, spatial, temporal, quadrature, operator approximation, algebraic solver, coupling, truncation, statistical. Algebraic tolerances stay below the discretization errors being measured.
   - **Self-sustained oscillators** also report:
     - sounding pitch (cents);
     - oscillation and extinction thresholds (%);
     - regime identity (register, slot, slip fraction).

     Δt can warp pitch, and a rung near a bifurcation can flip regime while the spectral distance stays small.
   - **Two kinds of operating point.** Convergence is qualified at **interior operating points**. Thresholds, hysteresis and intended register transitions are **mapped separately**, because a phrase-capable instrument cannot ignore near-boundary behaviour.
   - **Stochastic and chaotic outputs** converge in distribution (seed ensembles).
2. **Model-form ladder** (validation of reduced dimension).
   - **Comparison.** Every reduced-dimension or lumped component is compared with the next model up across its declared domain, not at two points. The rungs are:

     | Component | Rungs |
     |---|---|
     | duct | 1-D plane wave → continuous multimodal → 3-D linear FE/BEM with viscothermal layers |
     | shell | DKT/Kirchhoff–Love → shear-deformable shell → 3-D p-FEM through thickness |
     | valve | lumped → modally truncated continuum → full continuum with resolved channel flow |
     | flow | quasi-steady → resolved compressible flow |

   - **Refusal band.** The band where model-form error exceeds the gate tolerance is computed, never asserted, and the module refuses outside it.
3. **Operators.** Rational or modal approximations **of computed operators** count as discretizations: a BEM radiation impedance, a body admittance, the exact Zwikker–Kosten viscothermal operator.
   - They are accepted when their order ladder converges and passivity is enforced as a constraint.
   - Any perturbation that passivity enforcement introduces counts in the operator error budget.
   - Fits **to measurements or recordings** are never the instrument.
4. **Source ownership.** Each physical sound source has exactly one owner in the model. A resolved boundary trace must not also receive a jet-noise, dipole or turbulence term that represents the same source (Ffowcs Williams & Hawkings 1969 accounting).

**Validity audit.** It is computed per description, and the module refuses outside it.

- **Transverse modes and evanescent fields.** The first asymmetric mode propagates above 1.841·c/(2πa), and the first axisymmetric mode above 3.832·c/(2πa).
  - Clarinet bore (a = 7.3 mm): 13.8 / 28.7 kHz. The first asymmetric mode therefore propagates **inside** the audio band once holes or bends break symmetry.
  - Clarinet bell rim (a ≈ 30 mm): about 3.4 kHz.
  - Trumpet bell (a ≈ 60 mm): 1.7 / 3.5 kHz.
  - Tuba bell (a ≈ 0.2 m): 0.5 / 1.0 kHz.
  - **Below cutoff, a 1-D model is admissible only with bounded evanescent-field error** near features (mouthpiece, tone holes, bends, closely spaced discontinuities). Local 3-D subdomains (M2.4) supply that bound.
- **Viscothermal loss.** The exact Bessel (Zwikker–Kosten) solution holds at every shear number for **cylindrical linear sections**. Webster–Lokshin fractional terms are only its wide-tube asymptote, with error of order 1/r_v; r_v ≈ 6.5 for a 1 mm oboe staple at 100 Hz. Reed channels, flue windways and pad gaps fail the cylindrical model for different reasons: slit geometry, mean flow and separation. They need resolved models (M2.4 thin-gap linearized Navier–Stokes; M3.2 flow).
- **Edges.** Edges are linear only while acoustic displacement ≪ edge radius. Beyond that, local nonlinear (vortex-shedding) losses at tone holes and bells apply.
- **Propagation.** Propagation is linear only while the shock-formation distance ≫ the bore's cylindrical length. This fails for trombone and trumpet at ff (Hirschberg et al. 1996, trombone).
- **Air drag on strings.** The oscillating-cylinder Stokes solution applies only at Keulegan–Carpenter number ≲ 1. A guitar pluck reaches about 6.

**Acceleration contract (replaces v3's "closure contract").** A closure (jet-drive law, contraction coefficient, flow-noise source statistics) is admissible only as an **acceleration**. Its inputs are minted from resolved simulation or published measurement in the same regime. Its error against a coupled resolved run is measured across the declared domain. It refuses outside that domain. It never qualifies a render for Q3/Q4: the acceptance render uses the resolved physics (OD11).

### 1.3 Data, never fitting; the player and the microphone are declared before grading

**Recordings are for grading only.** Every physical parameter comes from one of:
- geometry or material measurement;
- first principles;
- a resolved higher-fidelity simulation.

**Construction, setup and performance are separate layers.**
- **Construction design** is the instrument as built: geometry, materials and the maker's tuning (bar undercuts, tone-hole placement and undercutting, reed scraping, voicing).
  - It comes from the measured specimen (§1.4).
  - Where no specimen is measured and the plan itself designs a part (e.g. an undercut optimized to 1:4:10), the part is a **declared synthetic design**. Every targeted ratio is disclosed, and the result is excluded from same-specimen predictive-validation claims.
- **Setup tuning** is what a technician or player sets before playing: string tensions, barrel and tuning-slide pull, reed position, pipe stoppers.
  - It follows a frozen protocol to equal temperament at note A4 = 440 Hz.
  - The protocol is defined independently of the graded recordings, and is applied per independently tunable element (each string, each pipe).
  - It is never set to the recording's pitch.
- **Performance controls** are what the player does while playing (next paragraphs).

**Allowed adjustments against a reference, and nothing else:**
- the frozen setup protocol above;
- a performance level inside the published range, chosen by a non-spectral target (pp < mf < ff).

**Player inputs** come from published player measurements or owner-authored performances, never from the graded recordings. They include breath pressure and its ramps, tonguing, embouchure, bow force/speed/position/tilt, pluck force and release, and vibrato.

**Per-note player controls that no measurement fixes are set by a declared non-spectral objective evaluated on the model alone,** and recorded per render. Examples are lip tension and rest aperture, jet length/offset/lip coverage, lip force on a reed, bow–bridge distance, and pluck point. Allowed objectives include:
- equal-temperament pitch with the measured lip-resonance/playing-frequency ratio (Newton, Campbell & Gilbert 2008);
- slot centre from the regime map;
- minimum threshold pressure.

Without this rule, agents would choose these values by listening to the reference.

**Observation nuisances** (microphone distance, azimuth, elevation, player posture, chair/floor) form a **prior declared before grading** and are **marginalized**: the score is the ensemble median, never the best member.

**Hold-out and freeze.**
- Material, geometry and player inputs are frozen independently of the graded recordings.
- Notes, dynamics and phrases are split into development and held-out sets, and acceptance is reported on the held-out set.

**Same-specimen vs population validation.**
- **Same-specimen predictive validation** compares a description of *a characterized physical instrument* with a recording of *that* instrument.
- **Population realism** compares a typical description with recordings of other instruments of the type, using same-type rulers.

The two are reported separately. Where a co-characterized specimen exists (§8, M8-W7), it is the stronger evidence.

### 1.4 Real geometry, materials and the interacting player (owner clarification, 2026-10-15)

**The target is the real physical assembly, including its environment and player.** Use credible online CAD/scans or construct geometry in Blender or equivalent from calibrated photo references, drawings and measurements. Target exact shape and material accuracy; finite acquisition and numerical uncertainty must be stated rather than represented as exact knowledge. An attractive exterior mesh alone does not establish a physically faithful instrument.

- Preserve internal as well as external shape: bores, mouthpiece cavities, undercuts, tone-hole chimneys, wall/thickness fields, soundboard arching, bracing, joints, gaps, seals, supports and contact surfaces. Preserve the source geometry separately from its numerical mesh.
- Assign actual materials, spatial material directions, layers, moisture/temperature state and constitutive evidence to the parts. Colour or a texture is not a material measurement. Geometry cannot identify an unmeasured loss tensor or tissue law.
- Include interacting anatomy and boundary conditions: mouth, tongue, lips, teeth/palate, airway, tissue layers and active prestress, plus hands/fingers, compliant holding and other external bodies where they participate. Anatomical motion and contact are part of the coupled physical experiment; a prescribed breath-pressure curve alone does not replace them.
- Photo reconstruction must retain scale references, camera/view evidence and uncertainty in occluded dimensions. Use sections, maker drawings, CT/MRI or measurement for hidden features when available. Unknown dimensions stay unknown or explicitly bounded; do not choose them by matching the recording. Speech MRI is a source of anatomy, not evidence of a musician's playing posture or tissue properties.

**Simplification requires demonstrated acoustic equivalence.** This applies to geometry decimation/repair, omitted features or bodies, symmetry assumptions, dimensional reduction and substituted material descriptions. The detailed source assembly remains the reference. The existing §1.2 convergence ladders and §7 comparisons provide the test seam; this requirement does not commission a new harness.

1. Freeze the detailed and proposed descriptions, material states, player gestures and physical receivers. First resolve both numerical solutions sufficiently to decide the comparison. Refining the mesh of an incorrectly simplified shape cannot validate that shape.
2. Predeclare the operating domain and equivalence bounds: frequency band, notes/fingerings, dynamics, attacks/releases, articulation, receiver directions/distances and anatomy/contact configurations. Include intended register transitions, onset/extinction thresholds and hysteresis; a matching steady note cannot establish equivalent playing behaviour.
3. Compare raw calibrated pressure and absolute timing/level/pitch, directivity, partial/envelope evolution and the relevant §7 descriptors. Do not retune, normalize gain, pitch-shift, equalize, or choose the best receiver to hide a difference. Chaotic/stochastic cases require distributional comparison with uncertainty, rather than pointwise waveform identity.
4. Report a bound on the difference including numerical and sampling uncertainty. Acceptance requires that bound to lie inside the declared equivalence tolerance throughout the claimed domain; a nonsignificant difference is inconclusive. Test combined simplifications as well as individual changes. Finite experiments establish bounded evidence in that domain, not a universal proof of equality.
5. If evidence is absent or inconclusive, retain the detail. Passing a population ruler or a listening test cannot by itself authorize a geometry/material simplification. Changed geometry, contact, material state or operating domain requires the affected comparison to be repeated.

**Search, sensitivity and uncertainty (IndigoLeopard duel, 2026-10-16).**
- **Search.** A finite or adversarial search can falsify an equivalence claim, but it cannot certify a domain without coverage bounds. A search that finds no counterexample does not authorize removing detail.
- **Untested hypotheses.** Non-negligibility arguments stay untested hypotheses until measured.
  - A small passive load ratio (for example, tract impedance against bore impedance) does not authorize omitting anatomy.
  - No fixed dB heuristic becomes an acceptance rule.
  - Repeatability and JND data may inform a predeclared tolerance; measurement noise never loosens it.
- **Correlated unknowns.** Unknown geometry and material properties keep their **correlations**. Convenience priors (independent, uniform) are not used without independent evidence; genuinely supported distributions are allowed. Before new acquisition, a passive sensitivity screen ranks the three most consequential unknowns. An insensitive response from a coarse grader never justifies omitting one.
- **Reference state.** Every imported geometry declares its **reference and load state**: unloaded, strung, assembled or played. Using a loaded scan as the unloaded reference and then applying its load again double-loads the body. `examples/grand_piano/downbearing.rs` already refuses a loaded crown as its reference. An equilibrated initial-stress formulation may legitimately start from the current, loaded geometry.

The simplification comparison and same-specimen validation are separate: agreement with an inaccurate detailed model does not establish realism. Idealized fixtures and early milestones remain useful for testing modules and locating residuals, but do not qualify the final instrument until this standard is met. Runtime or browser constraints do not relax it.

---

## 2. Where we are

### 2.1 Sound (graded 2026-10-15 against University of Iowa recordings; chain and room per M2.6)

| Instrument | What was graded | Distance (dB) | Real-vs-real ruler (dB) | Nearest real instrument |
|---|---|---|---|---|
| Grand piano note A3 mf | literature-reconstructed composite Model D (authored fields; a population comparison, not a specimen), full chain to a mic | 11.7 (one-note descriptor result, historical; not whole-piano indistinguishability) | 11.1–13.5 (dynamics) | piano |
| Clarinet (shipped fixture) | 2.2 mm × 0.5 m cylinder, quasistatic reed, bore pressure | 23.1–23.5 | 10.9–23.4 | double bass, cello |
| Clarinet (7.3 mm per-note tubes) | diagnostic, bore pressure | 18.4–30.5 | 5–15 | clarinet on 7/8 notes, still farther than a real cello |
| Trumpet | 0.67 m cylinder+cones, 1-DOF lip, mouthpiece pressure | 24.6–32.3 | 3.5–21.5 | alto sax, horn, bassoon |
| Violin / cello | one open 0.65 m string, point bow, one-mode body | 15.8–27.9 | 6.6–19.5 | cello, alto sax, oboe |
| Flute | 0.30 m × 6 mm pipe (12 kHz) | 26.5–33.9 | 10–24 | wrong register (modes 2–4) |
| Guitar | 3 ideal oscillators + 3 resonators | 23.5–27.5 | 8.6–17.9 | cello/violin pizzicato |
| Xylophone | uniform rosewood prism | 28.4 | 7.4–16.8 | untuned 1 : 2.84 : 5.58 |

No orchestral instrument is described with real dimensions anywhere in the repository. The 2026-10-15 grader also has blind spots (§7.2): it skips the attack and the late decay, fits a single decay slope, and misbins stiff partials. These numbers are therefore an upper bound on how close we are, not a precise measurement. Appendix A has the per-render detail. The grader, renders and audit are in the listening pack (Appendix C).

### 2.2 Modules: headline findings (inventory at 47a361a87, corrected by the reviews)

1. **Strong generic foundations exist.**

   | Area | What exists |
   |---|---|
   | Modal methods | certified eigen-solvers (`fs-modal`); exact modal ZOH stepping (`fs-couple/src/modal_acoustic_time.rs`, 87 KB, ≤ 4096 modes) |
   | Port-Hamiltonian core (`fs-phs`, L2) | Gonzalez discrete gradients, DAE descriptors, and a `prepared` path with analytic tangents (`prepared/analytic.rs`), exact condensation (`condensed.rs`) and backtracking |
   | Distributed ports | `fs-ports` (L3): `PortValueShape::Field` and `PowerPairing::FieldDuality` (`lib.rs:202–262, 376`), moving-stream enthalpy accounting (`lib.rs:1208`) |
   | Energy-consistent exchange | held effort × Δx/h port contract (`render/schedule/force/coupled.rs:1–15`) |
   | Plates (`fs-plate`) | orthotropic DKT plates with stiffeners and thickness maps |
   | Shells (`fs-plate::shell::survey`) | flat-facet shells over arbitrary orthotropic midsurfaces |
   | Beams (`fs-plate::stiffened`) | tapered Hermite beams |
   | Solids (`fs-solid::linear3`) | linear volumetric modal mechanics with consistent stiffness and mass (P1) |
   | Exterior acoustics (`fs-bem`) | Burton–Miller exterior Helmholtz BEM with radiation-impedance matrices and SH directivity; **near-field pressure receivers** (`near_field.rs:95`); **first-order pressure/velocity receivers** (omni/cardioid/figure-eight, `near_field_velocity.rs:104`); reciprocal Galerkin radiation (restricted resonance-free band) |
   | Ducts (`fs-duct`) | plane-wave and multimodal transfer matrices with exact Bessel viscothermal losses from `GasState` |
   | Gases (`fs-material`) | mixtures with humidity and CO₂ that reach audio |
   | Contact | distributed unilateral contact (`fs-dcontact`, energy-exact penalty potentials); finite-gap contact BEM (`fs-contact`), with `interface_binding` to `fs-matdb::InterfaceSystemCard` |
   | Periodic regimes | harmonic balance and continuation (`fs-orbit`) |
   | Psychoacoustics | `fs-psycho`, validated on ISO signals and used in five fs-couple tests, but not in grading |
   | Materials | `fs-matdb` (provenance-bearing) |
   | Ladders | `fs-ladder` (`Ladder`, `Refine1d`) |
   | Lattice-Boltzmann | `fs-lbm` with 2-D grid refinement (`refine.rs`) |

2. **Fragmentation.**
   - **Strings live in five places.**
     - The best physics is in `examples/grand_piano`: reciprocal stiff strings, both polarizations, and longitudinal phantom-partial coupling. It is compiled into the browser engine by `#[path]`; engine-piano now `#[path]`-includes 11 grand_piano files (`refmods.rs:21–46`).
     - The fifth place is the percussion impact host, `render/plate/impact/{string,supported,linear/wire}.rs`. Since f9c6f1fcb and 9b9bb95a5 it carries material strings (`bowed_string/runtime/reciprocal/bank/nonlinear{,/supported}.rs`) under fs-phs midpoint stepping.
   - **Valve owners.** The reed family has four:
     - `DynamicAperture`, including single and symmetric DKT plate reductions;
     - `PairedAperture`, with independent blades;
     - `ModalAperture`: 1–32 signed modes from 1–2 material plates (d66962bfa, 4e3d1de55), hosted by `LoadedAperture`/`TwoSidedAperture`;
     - the CLI's `ReedBoreVoice`/`BeatingReed`.

     Brass lip and glottis come on top. The CLI reaches only `ReedBoreVoice` and the `DynamicAperture` performers. There are two copies of the scalar condensation solve (`paired.rs`, `uniform.rs`) and three copies of `with_plate_forces`.
   - **Three hammer-felt adapters** sit around one `WoolFelt` law.
   - **Three baffled-piston quadratures:** `fs-bem/src/helmholtz.rs:1484`, `fs-phs/src/lib.rs:1707`, `fs-couple/src/thin_plate.rs:1748`.
   - **Three modal viscous-loss realizations:**
     - `thin_plate/material_loss.rs`: diagonal isotropic; refuses others at :26;
     - `plate_viscosity.rs`: full tensor projection over the fs-material law owner `plane_stress_viscosity.rs` (174af92bc);
     - `bernoulli_aperture/modal/viscosity.rs` (70ccffa3b).
   - **Two lumped-cavity families on different gas owners:** `cavity_modal` (`GasState`) and `cavity_phs::piston` (`ThermochemicalGasModel`).
   - **Three membrane-energy owners:**
     - `fs-solid::membrane` (d2fcde2b3);
     - `fs-plate::shell::preload::equilibrate_shell_materials` with independent A (98bfffca7);
     - `shell::head::nonlinear::MembraneReduction`.
   - **Three physical-outlet observer adapters** over the shared `BaffledPressure` owner: `OutletPressure` (reed scheduler), the `TwoSidedPressure` environment (6113eca8c), and `AperturePerformance` (cf8b1e2b0).
   - **About 55 per-feature file formats.**
   - **A closed `RenderVoice` enum** (8 variants; `ControlDelta` has 4).
   - **A fixed-slot `AcousticAssembly` realizer** that no CLI reaches.
3. **Finished but unwired:**
   - **Engines:** `fs-ports` coupling schema and audits; IQN-ILS partitioned coupling; `fs-time` integrators; `fs-orbit` (only `tests/hb_gates.rs`); `fs-modalid` (no consumer).
   - **Physics:**
     - `porous.rs`;
     - `ThermochemicalGasModel` for acoustic consumers: it is consumed only by `cavity_phs::piston` and `material_circuit/thermal`;
     - `fs-plate-transient`;
     - BEM fluid loading (examples only);
     - SH directivity in wind/brass audio.
   - **Interfaces and data:**
     - `ModalAperture`, `TwoSidedAperture`/`PairedAperture` and `new_with_exterior_gas` have no CLI path.
     - `OutletPressure` is wired into the reed scheduler and `music_render` wind/ensemble (2ddc6e603), but unverified.
     - Tissue, contact and string material packs.
     - Nine of sixteen typed gesture targets.
4. **Genuinely missing physics, grouped by tier:**

   | Tier | Missing |
   |---|---|
   | Excitation and contact | bow ribbon (width, hair compliance); thermal/plastic third-body friction; tangential contact compliance; finger and plectrum contact; lip collision and lip continuum; player airway coupling |
   | Strings | wound-string mechanics; nonuniform string FE; geometrically exact dynamics; string–fret–string contact |
   | Ducts and acoustics | staple pressure recovery; nonlinear propagation in bores; continuous multimodal propagation (today's is staircased); exact unflanged radiation; 3-D interior acoustics with boundary layers; axisymmetric bell FE; multi-aperture radiation in the oscillation loop; complex/directional fast BEM |
   | Structures | shear-deformable curved shells and junctions; TET10/hex solids with orthotropic viscoelastic tensors; large and nonlinear eigenproblems |
   | Flow | resolved compressible flow at playing Reynolds numbers (D3Q27, transparent refinement); zonal flow–acoustic coupling |
   | Runtime | conservative moving-geometry state transfer; energy-conserving multi-rate coupling; a generic instrument description and renderer; a generic WASM engine |
5. **Defects** (Appendix B gives file:line).
   - **Integrator warping.** For quadratic storage, the Gonzalez discrete gradient is exactly implicit midpoint. **Every** linear subsystem stepped by `fs-phs` is therefore frequency-warped: −58.6 c at 5 kHz at 48 kHz, −3.9 c at 192 kHz. The Kirchhoff–Carrier path is one instance.
   - **Unverified `impedance_line` convention.** It may time-reverse its impulse response. A comment at `driving_point.rs:205–208` asserts the opposite, so this must be adjudicated by causality tests. It also forces Z(0) = 0, which erases the compliance state of closed volumes.
   - **Characteristic-line FFT.** The default 4-round-trip FFT wraps around (measured 92 c on a horn).
   - **Passivity repair.** Global-gain repair is checked only on a sampled grid.
   - **Newton tolerance.** Newton tolerance is relaxed implicitly and scaled inconsistently. The residual is disclosed in `StepRecord::solver_residual`, but the iterate is still accepted.
   - **Gap flow.** `bernoulli_volume_flow` is singular at closure and has no viscous term.
   - **Traction ramp.** A discontinuity in `regularized_traction_1d` at v = vs.
   - **Swallowed refusals.** `.unwrap_or(0.0)` swallows obstacle-friction refusals.
   - **Observer.** The reed observer adds internal bore pressure to radiation plus a dimensionally wrong jet term.
   - **P1 tets.** P1 tetrahedra are used where bending dominates (slow, locking-limited convergence).
   - **Mislabel.** A Helmholtz-resonator end correction is labelled unflanged but is flanged.
   - **Build break.** `examples/material_strings.rs` does not compile.
   - **Stale CONTRACTs.**
6. **Breaches of §1.3:**
   - **Drum heads.** fs-vkplate drum-head tensions and two-point T60 loss fits were tuned to recordings.
   - **Jawari contact pack.** It was calibrated to a fixture's output.
   - **Piano axial projection.** `BRIDGE_AXIAL_PROJECTION = 0.06` is a hand-estimated effective coefficient ("~1.8° downbearing plus a 28 mm bridge-cap rocking lever arm"). It was introduced under a bead (na4ju) whose DONE-WHEN targeted the Iowa recording's phantom-partial levels. It is not demonstrably a fit, but it is authored rather than computed.
   - **Hand-typed constants throughout.** In the wind and bow paths: reed lay K = 1e7·w, face length 0.025 m, ζ = 0.35, a 5 Pa priming push, the brass lip mass/stiffness/damping, the glottis constants, and the bow capture thresholds. In the piano: the felt parameters, the 0.60 drag rate, and the authored bending spectrum.
7. **Data gaps.**
   - **Wood packs** carry ratios at 12 % moisture only.
   - **Loss data** are single longitudinal values, with no frequency- or direction-dependent loss.
   - **Missing materials:** nylon, gut, fluorocarbon, windings, horsehair, rosin, cork, pad leather, ebony, grenadilla, wet cane, fingertip and tongue.
   - **No measured bores, bodies, admittances or impedance curves,** except an Ernoult lab cylinder, Carcagno guitar-body modes and Assaneo tract tables.
   - **No co-characterized specimen** (geometry + materials + passive response + sound of one physical instrument).
8. **Verification debt.**
   - About 220 commits from 10-12 to 10-15 (bowing, double reeds, flue outlets, material strings) state that native Rust tests were never executed, and HEAD stopped compiling on 10-15.
   - **A remote lane without a toolchain is still pushing unverified instrument code**, including commits `6113eca8c`, `2ddc6e603`, `174af92bc`, `3ae3978ac` and `76764a701` ("Local cargo/rustc absent; native tests … unverified"). It touches `bernoulli_aperture`, `plate_viscosity.rs`, `music_render/wind_input.rs` and `OutletPressure` (OD13).
   - **After 9267fd56d** it pushed 36 more commits by 20:38 on 10-15, and more since (e.g. 39f15c202). About 27 of them state no native run:
     - ModalAperture multimode reeds: d66962bfa, 71bf3b183, 4e3d1de55, 70ccffa3b, cf8b1e2b0;
     - shared-support strings: 2d269fd3d, 9b9bb95a5, f9c6f1fcb;
     - membranes: d2fcde2b3, 6d3dee80d, f03a59872, 27e46fd99, ee6fcaba4;
     - shared cavities and gas pistons: f98b823b1–d3e926481, a5fb73382, 708f9882d;
     - engine-piano oversampling and scale preparation, plus `public/live/piano/*.js`: 2fa573a75, 9583f8384, b1efe209e (WASM not rebuilt).
   - **CI is not evidence.** `instrument-materials.yml` has never passed. Bootstrap fails on the franken_numpy `make.bat` line endings, and every cargo step fails `fsci-special = "=0.4.0"` against `constellation.lock`'s frankenscipy 0.2.0.
   - **Body-cavity paths never ran at 47a361a87.** Until 708f9882d, `HelmholtzModalCavity::neck_component` refused at construction.
9. **Geometry and anatomy.**
   - **What exists.** Import and description seams exist.
     - The piano reads physical OBJ, part and section inputs: `examples/grand_piano/mesh_import.rs`, `solid_import.rs` and `crowned_import.rs`; `downbearing.rs` requires an unloaded reference.
     - Percussion has specimen and measured-profile import seams.
     - `fs-io` reads meshes.
   - **What this audit has not established.** For any added family, there is no complete acquired, co-characterized 3-D specimen assembly carrying materials, reference/load state and anatomy together with acceptance evidence.
   - **Current fixtures.** The wind and string fixtures graded in §2.1 use idealized bores (cylinders, cones) and bodies (rectangles, modal charts).
   - **Anatomy.** Interacting anatomy is limited to area-function tract tables (Assaneo) in `tract.rs` and generic contact and source hooks. There is no qualified, acquired playing anatomy: lips, tongue, teeth, hands and fingers in playing postures.
   - **Tools.** Blender 5.3 is installed on the Mac and runs headless.

---

## 3. Architecture: the generic instrument-physics stack

### 3.1 Stack tiers (S0–S7 are plan tiers, not workspace layers)

```
S7  Validation & grading     grader · references · QP gates · module conformance · fs-modalid · fs-psycho
S6  Descriptions & data      Assembly Description · Performance · player data · geometry & material data
S5  Coupling & runtime       port graph · integrators · multi-rate · moving domains · renderer · ladders · regimes
S4  Contact & friction       normal laws · tangential laws · nonsmooth limit rung · driven bodies · contact heat
S3  Flow & excitation        moving-wall flow interfaces (valves) · resolved flow · zonal coupling · accelerations
S2  Acoustic fields          1-D/multimodal ducts · DG time-domain ducts · nonlinear propagation · 3-D interior · exterior · microphones
S1  Structures               strings/rods · shells/plates · 3-D solids · assembly · large eigenproblems
S0  Media & materials        gases · passive anisotropic viscoelastic tensors · tissue · third-body layers/polymers/felt · matdb
```

**Workspace layers (AGENTS.md) govern crate dependencies, not these tiers.** Most physics crates are L3 (FLUX). `fs-phs` and `fs-orbit` are L2; `fs-matdb`, `fs-modal` and `fs-sparse` are L1; `fs-vkplate` is L4. L3 crates may depend on each other only acyclically, and on L0–L2. The consequences are concrete:
- The port-graph **solver lives in `fs-couple` (L3)** and consumes `fs-phs` (L2) systems and `fs-mbd`/`fs-contact` (L3) nonsmooth nodes through traits.
- `fs_orbit::OrbitProblem` is implemented *for* a graph inside `fs-couple`.
- New L3 crates that `fs-couple` consumes (`fs-string`) must not depend on `fs-couple`. The modal ZOH core moves out first (M1.1-W1a).

| Module | Home crate (layer) | Notes |
|---|---|---|
| M0.1 gas | fs-material (L3), fs-thermochem (L3) | |
| M0.2–M0.4 materials | fs-material (L3); packs in fs-matdb (L1) | the matdb.pack schema is frozen; changes go through schema-policy |
| M0.5 resolution | fs-matdb (L1), fs-matdb-store | |
| M1.1 strings | **fs-string (new, L3; OD10)** | below fs-couple; the ZOH core moves first |
| M1.2 shells/plates | fs-plate (L3); `thin_plate*`, `plate_viscosity.rs` live in **fs-couple** today | |
| M1.3 solids | fs-solid (L3) | |
| M1.4 assembly | fs-plate / fs-solid; orchestration in fs-couple | |
| M1.5 eigen/LA | fs-sparse (L1), fs-modal (L1) | |
| M2.1 ducts | fs-duct (L3) | |
| M2.2–M2.3 DG ducts | fs-duct (L3) | new module, not a new crate |
| M2.4 interior | decided by M2.4-W0: **fs-feec (L3)** extension preferred; supersedes bead 0ja4's interior scope | |
| M2.5 exterior | fs-bem (L3), fs-fmm (L2) | |
| M2.6 microphones | fs-bem receivers; observation in fs-couple | |
| M3.1 valves | fs-couple (pattern over fs-phs/fs-contact/fs-duct) | |
| M3.2–M3.3 flow | fs-lbm (L3); interface in fs-couple | |
| M4.1 normal contact | fs-contact (L3), fs-dcontact (L3) | |
| M4.2 tangential | fs-tribo (L3) | |
| M4.3 nonsmooth | fs-mbd (stepper) + fs-contact (cones) per bead oh0i | |
| M4.4 driven bodies | fs-couple | |
| M5.x runtime | fs-couple (L3) | fs-phs (L2) underneath |
| M6.x descriptions | fs-scenario (L3) on the frozen `scenario.ir` schema (migration via schema-policy) | |
| M6.3 geometry acquisition | `tools/geometry/` (Blender headless scripts, photo calibration, import checks); imports through fs-io (L2); small assets under `data/instruments/geometry/<specimen>/`, large scans off-repo with hashes | not a crate and not a separate governance system: it feeds the §8 ledger and the M6.1 descriptions |
| M6.4 player anatomy | parts in the fs-scenario description; physics from M0.3/M1.3/M2.1/M2.4/M4.1 | anatomy is description data, never instrument code |
| M7.1 grader | `apps/instrument-lab/tools/instrgrade/` (new files; the site maintainer is notified) | |

### 3.2 Ports and interfaces

Modules communicate only through the following.

- **Power ports.**
  - **Shapes.** Each port is an effort/flow pair with units, of scalar, vector or **field** shape (`fs-ports` `PortValueShape::Field` with `PowerPairing::FieldDuality`). Typical pairs are force/velocity, moment/angular velocity, pressure/volume flow, and temperature/entropy flow.
  - **Distributed ports** declare their trace space, orientation, frame and integration measure.
  - **Nonmatching meshes** exchange through power-preserving mortar transfer (Cardoso-Ribeiro, Matignon & Lefèvre 2021), so discrete virtual work is preserved.
  - **Flow boundaries** also carry mass, momentum, species and total energy, with entropy production accounted for (`fs-ports` moving-stream enthalpy).
  - **Modulated interconnection.** Coupling vectors may depend on time and state, e.g. a bow position or finger stop moving along a string.
- **Material handles**: resolved `fs-matdb` receipts, never bare numbers.
- **Refinable axes**: every axis has its declared convergence/error model (§1.2) and coupled groups, so the ladder harness (M5.5) can refine any of them without knowing the module.
- **Validity domains**: every module refuses outside its domain with a typed error.
- **Entity paths**: gestures and configuration changes address description parameters by typed path (`EntityId`), not `u8` indices.

### 3.3 Integrator doctrine (applies to every render path)

1. **Linear parts are stepped by exact (ZOH/exponential) or per-mode frequency-matched energy-consistent propagators** (Chabassier & Joly 2010). They are never stepped by implicit midpoint or trapezoid, which warp frequency. Generalized-α is forbidden on render paths: with ρ∞ = 1 it is trapezoidal, and with ρ∞ < 1 it adds algorithmic damping comparable to wood's η.
2. **Nonlinear storage** uses discrete gradients or energy quadratisation (SAV/IEQ; Ducceschi & Bilbao 2022, non-iterative where applicable).
3. **Port exchange** uses held effort × Δx/h (the `coupled.rs` contract), so exchanged work is exactly e·Δx on both sides.
4. **Admission rules.** A step is admitted only if ω·Δt is below the declared bound for every non-exact part, Δt ≤ T_contact/N for every stiff contact, and an anti-aliasing filter is applied before every rate change.
5. **Energy closure** is to the declared algebraic-tolerance budget. Newton tolerance is tied to that budget and **refuses** rather than accepting a stalled iterate.

### 3.4 Contracts every module satisfies

1. **Generic.** No instrument names appear in non-test code; a description supplies every number.
2. **Ported.** The module exposes `fs-ports` schemas (M5.1-W0) and composes through the port graph. Its energy ledger closes to its budget.
3. **Converged.** Refinable axes are declared, and a G1 convergence test uses a manufactured or analytic solution on ≥ 3 rungs. Where the module reduces dimension, it also has a model-form comparison.
4. **Validated.** At least one literature or measured benchmark from Appendix D.
5. **Honest.** It refuses outside its validity domain with a typed error. There are no silent fallbacks: no `unwrap_or(0.0)`, no best-iterate acceptance, no global-gain clipping.
6. **Material-driven.** Parameters come from packs or geometry. Authored constants appear only in test fixtures, and are named as such.
7. **Workspace rules.**
   - **Determinism.** Elementary math goes through `det::`, enforced by check-libm. Stochastic sources use `fs-rand` counter streams keyed by (description, source, block).
   - **Cancellation and checkpoints.** `Cx` is polled per audio block. Renders longer than 10 minutes checkpoint and resume byte-identically.
   - **No new unsafe** (SIMD only via fs-simd).
   - **wasm32.** It is required only for crates reachable from `engine-perform`, which are listed explicitly; each such bead runs `cargo check --target wasm32-unknown-unknown -p <crate>`.
   - **CONTRACT.md** with the ten AGENTS.md sections, before the crate becomes a dependency target.
8. **Natively verified** before push (OD3).

### 3.5 Why this architecture

- **It satisfies "no instrument modules".** Every family needs the same small set of physics (structures, fields, flows, contacts), only composed differently.
- **It reframes the work.** "Make the clarinet perfect" becomes "make the duct, valve, radiation and contact modules exact". That also makes saxophone, oboe and bassoon nearly free.
- **Every module has a non-musical customer**, noted per module, which matches the project's identity.
- **Agents can work in parallel** along tiers without colliding, with instruments as integration milestones.

---

## 4. Module specifications

Each module lists:
- **Purpose**, including non-musical uses;
- **Now**, with file:line where load-bearing;
- **Target**;
- **Numerics & convergence**;
- **Validation**;
- **Consumers**;
- **Work items** (one line each; size, dependencies and DONE-WHEN are in Appendix E).

### S0 — Media and materials

#### M0.1 Gas media

**Purpose.** Thermodynamic and transport properties of any gas mixture at any state, for acoustics and flow (HVAC, combustion acoustics, breath).

**Now.** `GasState`/`GasSpec` with:
- mixtures (Wilke; Wassiljewa–Mason–Saxena);
- moist air (Buck);
- named species with frozen 298 K cp.

It reaches audio (`acoustic_realize.rs:1012`, and `RegionalDuct` per segment). `ThermochemicalGasModel` (NASA-9 cp(T), `fs-material/src/thermochemical_gas.rs:53`) is consumed only by `cavity_phs::piston` and `material_circuit/thermal`. No `GasState` acoustic consumer (ducts, radiation, `cavity_modal`) uses it. `fs_thermal::enthalpy::advection::Advected` (52457eb93, 0403460cc) transports temperature by first-order upwind with prescribed flow. It has no species and no variable density. There is no molecular-relaxation absorption for breath mixtures and no axial acoustic field.

**Target.**
- γ(T) and cp(T) from NASA-9 for every acoustic consumer.
- CO₂/H₂O vibrational-relaxation absorption for arbitrary mixtures (ISO 9613 form generalized by species).
- Continuous axial fields of temperature and composition along a duct, including the **steady played breath field**. Breath enters at about 35 °C, saturated, with a few % CO₂, and cools along the bore. This changes intervals between registers, which the one allowed tuning action cannot absorb (Fletcher & Rossing 1998).

**Validation.**
- NASA-9 against JANAF.
- Relaxation absorption against ISO 9613 and published CO₂ data.
- Humid sound speed against Cramer, pointwise.

**Consumers.** Every duct, cavity, radiation and flow module.

**Work items.**
- **M0.1-W1:** wire `ThermochemicalGasModel` into `GasState` consumers.
- **M0.1-W2:** relaxation absorption for mixtures.
- **M0.1-W3:** axial temperature/composition fields plus a steady played breath-field generator (from a thermal duct model or measured profiles).

#### M0.2 Passive anisotropic viscoelastic solids

**Purpose.** Frequency-, temperature- and moisture-dependent anisotropic stiffness and loss for any solid (composites, wood products, polymers).

**Now.**
- **Elastic laws:** `IsotropicElastic` and `OrthotropicElastic`.
- **Scalar 1-D loss laws:** `FractionalZener`, `GeneralizedMaxwell`, `lower_to_prony` (band-certified), and isotropic `ThermoelasticZener`.
- **Plate material loss** in `thin_plate/material_loss.rs` is diagonal isotropic Kelvin–Voigt and **refuses others** (:26, in fs-couple).
- **Full-tensor Kelvin–Voigt modal projection** with Orthotropic21 exists in the bowed bank: `plate_viscosity.rs` over the law owner `plane_stress_viscosity.rs` (174af92bc, unverified).
- **A third full modal viscous solve** sits in `bernoulli_aperture/modal/viscosity.rs` (70ccffa3b).
- **No frequency dependence.** Loss is frequency-independent, and there is no Prony tensor.
- **Wood packs** store FPL ratios, while the orthotropic resolver needs E1..E3 (`state_point.rs:1995`).
- **No temperature/moisture shifts.** No WLF/Arrhenius or sorption models.

**Target.** A complex-modulus **tensor** law that is **passive by construction**:
- **Prony form.** In Mandel coordinates, equilibrium stiffness C_eq ≻ 0, every branch tensor C_j ⪰ 0, and τ_j > 0. Equivalently, branches relax per Kelvin eigentensor (Mehrabadi & Cowin 1990; Biot 1954; Hanyga 2014). It has an explicit free energy and non-negative dissipation, preserved through rotations, plane-stress elimination and projection.
- **Fractional laws** are lowered to such branches with a certified band error.
- **Shifts and coupling:** temperature (WLF/Arrhenius), moisture (sorption), thermoelastic coupling.
- **Counted once:** intrinsic, thermoelastic, contact, support and radiation losses each appear exactly once in the budget.
- **Certified modal damping.** For lightly damped bodies, per-mode loss comes from the energy quotient, which is exact to O(η²) away from clusters. A certified non-proportionality bound (off-diagonal modal damping over frequency separation) decides when complex modes are required. Such modes fall back to complex-mode or internal-variable integration. This avoids internal variables in 10⁴-mode bodies without being a reduced engine, because its error is bounded and refused when exceeded.

**Validation.**
- Kramers–Kronig consistency.
- Passivity under random multiaxial strain histories.
- Reduction to isotropic Zener.
- Published wood loss against frequency and direction (Ono & Norimoto; Obataya; Brémaud) once M8 data exist.

**Consumers.** Plates, shells, solids, strings, contact laws, valve continua.

**Coordination.** MR03 loss models (ui2if.1.3, in progress) and MR16 wood bundles (ui2if.2.8). Do not fork the loss law.

**Work items.**
- **M0.2-W1:** FPL-ratio → full orthotropic tensor adapter in the state-point resolver.
- **M0.2-W2a:** passive tensor Prony law.
- **M0.2-W2b:** tensor `lower_to_prony` with certified band.
- **M0.2-W3:** orthotropic material loss in plates/shells by extending `plate_viscosity.rs` (no third path), removing the refusal.
- **M0.2-W4:** temperature/moisture shift factors.
- **M0.2-W5:** certified modal-damping reduction with complex-mode fallback.

#### M0.3 Soft tissue

**Purpose.** Lips, fingertips, tongue, vocal folds, hands holding instruments; biomechanics generally.

**Now.** Neo-Hookean and Mooney–Rivlin (unwired for tissue). Tissue packs (perioral lip 45 kPa, vocal-fold moduli and tan δ) are read only by a sanity test. The perioral pack comes from one quasi-static indentation (n = 1, a teeth-backed 9 mm probe, 3.8 N over 3 mm), reduced to a 45 kPa effective stiffness. That is an apparatus response, not a portable intrinsic dynamic law.

**Target.**
- **Finite-strain visco-hyperelastic laws:** Ogden, Holzapfel–Gasser–Ogden fibre-reinforced, and finite-strain Prony branches.
- **Layered tissue** (vocal-fold cover/ligament/muscle; lip mucosa/muscle).
- **Active prestress:** muscle activation as a prescribed active stress, set by player data.
- **Incompressibility** via the mixed formulation (M1.3-W2b).
- **Intrinsic vs apparatus data.** Packs separate intrinsic constitutive data from effective apparatus responses.
  - An apparatus response is used only by forward-modelling the measured apparatus (probe, backing, thickness, prestress), and only when that apparatus is sufficiently described.
  - One scalar cannot identify thickness, backing, active prestress and a viscoelastic spectrum. Admissible law families are kept, and loading conditions are held out for validation.

**Validation.** Treloar data; published indentation curves for lip, fingertip and vocal fold (Zhang 2016 review).

**Work items.**
- **M0.3-W1a:** Ogden and HGO hyperelastic laws.
- **M0.3-W1b:** finite-strain Prony branches.
- **M0.3-W3:** active-stress tissue law.

#### M0.4 Third-body layers, polymers, felt, rubber

**Purpose.** Interfacial layers (rosin, brake transfer films, lubricant films), string polymers, mallet and hammer heads, pads.

**Now.** `WoolFelt` is rate-independent hysteresis and energy-exact, with hand-typed parameters everywhere. `thermovisco.rs` has temperature-dependent Maxwell with no glass-transition plateau. There are no data for rosin, nylon, gut, rubber or yarn.

**Target.**
- **(a) A generic third-body interfacial-layer law:** shear strength and viscosity against temperature across the glass transition (WLF), plus melting and layer thickness. Rosin and brake transfer films are packs.
- **(b) String polymers:** nylon, fluorocarbon and gut as M0.2 tensors with temperature and humidity shifts.
- **(c) Felt** with rate dependence (Prony branches on the hysteretic envelope), and parameters from measured coupons.
- **(d) Rubber and yarn** as visco-hyperelastic (M0.3) for layered mallet heads.

**Work items.**
- **M0.4-W1:** third-body interfacial-layer law.
- **M0.4-W2:** string polymer laws.
- **M0.4-W3:** felt rate dependence and measured-coupon parameter path.

#### M0.5 Material resolution

**Now.** Packs are mostly copied into test constants instead of being resolved. For example, `guitar_body_gate.rs:59`, and `bakeoff_dispersion`'s E = 200 GPa instead of the music-wire pack.

**Target.** Every render path resolves materials through `fs-matdb` queries. Authored test constants become test-only fixture packs.

**Work items.**
- **M0.5-W1:** pack resolution on render paths (after the M8-W2 data).

### S1 — Structures

#### M1.1 Strings and rods (`fs-string`, new crate, OD10)

**Purpose.** Tensioned 1-D continua: musical strings, cables, wire ropes, tethers, power lines, sutures.

**Now.** The physics is spread across four places:
- `string_specimen.rs` (764 LOC);
- the `acoustic_realize` realizer. Its KC path warps frequency, and its "second polarization" is an uncoupled detuned copy (0.85 scale, `acoustic_realize.rs:1422, 2081`);
- the bowed `ReciprocalStringBank`:
  - it has N strings on one shared body with sympathetic coupling, per-string footprints on one DKT eigenbasis (`bank/chart/spatial.rs`, 2d269fd3d), and full Kelvin–Voigt tensor damping (174af92bc);
  - its limits are one plane, small-angle torsion with no torque on any termination, penalty connections, compact pressure, a 128-mode cap, and ZOH stepping;
- the percussion impact host's strings: fixed-end Kirchhoff–Carrier filaments and `SupportedMaterialStrings` (f9c6f1fcb, 9b9bb95a5). The latter carries ≤ 32 strings with both ends on one translating 1-DOF rail, which is not a bridge/nut topology. Both are stepped by Gonzalez midpoint (DEF-04);
- `examples/grand_piano`. It has exact reciprocal string–board inertia, both polarizations with one tension, Kirchhoff–Carrier stretching including the bridge chord term, `fs-nlmodal::LongitudinalString`, bonded sections and finite-width hammers. `linear.rs:30` imports `fs_couple::modal_acoustic_time`.

Four in-progress piano beads edit these files: frankensim-5c0sj, stxqd, t0odx, i5yei.

**Target physics. Governing equations first:** a geometrically exact tensioned rod (Simo & Vu-Quoc 1986) with:
- transverse motion in two planes;
- axial inertia and longitudinal motion;
- torsion with twist–extension coupling;
- bending and shear (Timoshenko);
- eccentric and terminal loads;
- nonuniform and composite sections;
- losses from M0.2 and M0.4;
- air drag from M0.1/M3.2: the oscillating-cylinder boundary-layer solution where KC ≲ 1, and a resolved oscillating-cylinder characterization beyond, with its model-form error measured;
- temperature-dependent tension.

**Wound strings** have a resolved reference: core plus helical wrap as M1.3 solids, with wrapping prestress and core–wrap and wrap–wrap contact (M4.1/M4.2). They also have a homogenized section with an **enrichment path back to that reference**, plus inter-turn friction loss.

**Terminations** are reciprocal multi-DOF ports (3 translations + 3 rotations) at bridge/saddle, nut, tailpiece and peg. All strings on one termination couple, which gives sympathetic ringing of open strings. **Distributed contact** acts along the length: frets, fingerboard, fingers, neighbouring strings (M4.1). **Excitation** comes only through contact (M4.4).

**Numerics & convergence.** The modal and spatial implementations discretize **the same declared rod equations**.
- **(a) Exact modal.** Linear modal part by ZOH; nonlinear stretching/longitudinal coupling by SAV/discrete gradient (§3.3), replacing the warped KC path. Ladders: modes, substeps.
- **(b) Nonuniform FE.** Hermite/Timoshenko elements with tension, then the geometrically exact rod. Ladders: h, p.
- **Cross-checks.** Agreement of (a) and (b) verifies the discretization. Resolved solid specimens test the omitted cross-sectional and local-contact physics.

**Validation.**
- **Stiff strings:** Fletcher B inharmonicity; Chabassier, Chaigne & Joly (2013) for the piano string.
- **Motion and pitch:** Raman/Helmholtz bowed waveforms; Kirchhoff–Carrier pitch glide against analytic Duffing.
- **Partials and torsion:** longitudinal phantom partials (Conklin 1999; Bank & Sujbert 2005); torsional frequencies of violin strings (Woodhouse & Loach; Bavu).
- **Wound strings:** effective properties against measured data (Järveläinen; maker data, M8).
- **Piano regression:** byte-identical WAV.

**Consumers.** Guitar, violin family, double bass, piano (re-plumbed), harp (future).

**Work items:**
- **Lift (M1.1-W1a–d):**
  - **M1.1-W1a:** move the modal ZOH core to an acyclic owner.
  - **M1.1-W1b:** create fs-string; lift specimen and sections with `git mv`.
  - **M1.1-W1c:** lift linear/stretching/longitudinal/polarization; grand_piano consumes it, byte-identical.
  - **M1.1-W1d:** engine-piano off `#[path]`, byte-identical.
- **Generic bank and stepping:**
  - **M1.1-W2:** exciter-agnostic reciprocal string–body bank on fs-string, moved out of `bowed_string/`. Two planes, torsion with bridge torque, N strings per termination.
  - **M1.1-W3:** exact linear + SAV nonlinear stepping (removes DEF-04 on strings).
- **Discretizations and wound strings:**
  - **M1.1-W4a:** nonuniform FE of the rod equations.
  - **M1.1-W4b:** geometrically exact rod dynamics (model-form reference).
  - **M1.1-W5a:** resolved wound-string reference specimen.
  - **M1.1-W5b:** homogenized wound section with enrichment path.
- **Terminations, drag and contact:**
  - **M1.1-W6:** reciprocal 6-DOF termination ports.
  - **M1.1-W7:** BR-01. Axial projection computed from termination geometry (string axis vs support normal, bridge rocking from the bridge model).
  - **M1.1-W8:** air-drag operator (BL solution + resolved characterization).
  - **M1.1-W9:** distributed contact interfaces along the string.

#### M1.2 Plates and shells

**Purpose.** Thin-walled structures of any curvature and material: instrument tops and backs, bells, panels, enclosures, aircraft skins.

**Now.**
- **Plates and shells:** DKT orthotropic plates (validated Navier/Leissa/Olson–Hazell); flat-facet CST+DKT shells (drilling penalty, lumped mass); arbitrary orthotropic midsurfaces (`survey::MeshShell`, refuses non-manifold); eccentric stiffened shells; preload; `ShellReduction`.
- **Not validated:** no curved-shell literature validation, and no T-junctions in admitted paths.
- **Core limits:** Dirichlet-only supports; K and M only in the core.
- **Sections:** independent membrane A and bending D (`shell/section.rs`) through Green–Lagrange preload and DKT loss (`preload.rs:314`; 98bfffca7). There is no B coupling.
- **Membranes:** `fs-solid::membrane` (d2fcde2b3, 6d3dee80d): P1 orthotropic quadratic-strain sheets with Newton–Krylov equilibrium and dense prestressed modes up to 128 free DOFs (`vibration.rs:17, 71`). Zero damping; natively unverified.
- **Wiring:** `fs-plate-transient` is unwired. The bowed path uses a fixed 5×4 rectangle (`reciprocal/plate.rs:176`).

**Target physics.**
- **Elements.**
  - **Shear-deformable elements** are primary: MITC4+/MITC3+ (Lee & Bathe 2014; Ko, Lee & Bathe 2017), or 3-D high-order p-FEM through the thickness (Düster, Bröker & Rank 2001).
  - **Kirchhoff–Love/DKT is only a model-form rung.** Wood's low G_LT, G_RT (E_L/G_LT ≈ 15–20) make transverse shear audible above about 2 kHz: about 4–5 %, or 70–90 c, at 4 kHz for 3 mm spruce.
  - **Selection** is by geometry, h, p, thin-limit and distorted-mesh convergence. Membrane and shear locking are checked separately.
- **Junctions.** Individually oriented component surfaces are joined by **work-conjugate ties, contact or mortar**, with glue-joint compliance as data and explicit ownership of every interior and exterior fluid boundary. This covers T-junctions, closed boxes, ribs and linings. Shell–beam ties (braces, bass bar) and shell–solid ties (bridge feet, soundpost ends, blocks) use the same mechanism.
- **Laminate sections:** varnish, ground and finishes as viscoelastic layers (Schelleng 1968).
- **Supports and holding:** compliant ports to contacts and tissue (chin, shoulder, hand, knees, thigh, chest, forearm, endpin, floor).
- **Loss, nonlinearity and loading:**
  - **Loss** from M0.2 tensors.
  - **Geometric nonlinearity** (existing reduction, validated) where needed.
  - **Prestressed equilibrium** from M1.4.
  - **Fluid loading** from M2.4/M2.5.

**Numerics.**
- h- and p-ladders.
- Modal truncation ladders with residual flexibility.
- Certified eigen windows (M1.5).

**Validation.**
- **Static (MacNeal–Harder):** Scordelis–Lo, pinched cylinder, pinched hemisphere.
- **Vibration and locking:** Leissa (1973) vibration, and the Chapelle–Bathe asymptotic locking tests.
- **Instrument data:** bell partials; violin free-plate modes (Hutchins; Jansson); assembled violin signature modes; guitar top modes (Elejabarrieta et al.).
- **Shared tests:** with bead g42o's NAFEMS shells.

**Work items:**
- **Elements:**
  - **M1.2-W1a:** shear-deformable curved shell element with static benchmarks.
  - **M1.2-W1b:** vibration and locking benchmarks.
  - **M1.2-W1c:** 3-D p-FEM thin-wall cross-check.
- **Junctions:**
  - **M1.2-W2a:** work-conjugate ties and mortar between component surfaces.
  - **M1.2-W2b:** closed boxes and glue-joint compliance.
  - **M1.2-W3:** shell–beam and shell–solid ties.
- **Supports, loss and layers:**
  - **M1.2-W4:** compliant support/holding ports.
  - **M1.2-W5:** M0.2 damping in shell/plate assembly.
  - **M1.2-W7:** laminate sections.
- **M1.2-W6:** outline mesher with holes (f-holes, soundholes) and arching-surface generators.
- **M1.2-W8:** prestressed sheet membranes with one membrane-energy owner (§5.1).

#### M1.3 Three-dimensional solids

**Purpose.** Chunky and soft parts: bars, bridges, mallets, thick reeds, lips, fingertips, mouthpieces, machine parts.

**Now.**
- **Linear:** P1 tets with consistent K and M, and orthotropic or full tensors.
- **Statics:** finite-strain statics (`hyper3`); volumetric locking is not mitigated.
- **Valve card:** `reduce_valve` reduces tets to a modal valve card with an authored loss factor.
- **Missing:** dynamics beyond linear modal; viscoelastic tensor elements; TET10/hex.
- `fs-solid::membrane` is a 2-D elastic sheet; it does not advance M1.3-W2b.

**Target.**
- **Elements:** TET10 and HEX20/27 with consistent mass, with locking tests (Cook's membrane, bending patch). Practical bending accuracy is distinguished from incompressible locking.
- **Viscoelastic elements:**
  - linear orthotropic viscoelastic tensor elements (M0.2);
  - finite-strain visco-hyperelastic elements (M0.3) with **stable mixed (u,p)** or verified stabilization.
- **Dynamics:** exact modal for linear parts (§3.3) and energy-consistent integration for nonlinear parts.
- **Interfaces:** contact surfaces (M4) and acoustic skins (M2.4/M2.5).

**Validation.**
- Beam frequencies.
- Undercut-bar tuning (Bork 1995; Petrolito & Legge).
- Violin-bridge modes (Reinicke; Woodhouse).
- Treloar data, and the incompressibility patch test.

**Work items:**
- **Elements:**
  - **M1.3-W1a:** TET10 elements.
  - **M1.3-W1b:** HEX20/27 elements.
- **Viscoelastic elements:**
  - **M1.3-W2a:** linear orthotropic viscoelastic elements.
  - **M1.3-W2b:** finite-strain visco-hyperelastic mixed (u,p) elements.
- **M1.3-W3:** solid dynamics per §3.3.
- **M1.3-W4:** volumetric meshing from descriptions plus gmsh/STEP import via fs-io.

#### M1.4 Assembly

**Purpose.** Build a body from parts (top, back, ribs, bass bar, soundpost, bridge, neck, frets, supports) or any machine.

**Now.** Hand-composed per example. `vibroacoustic.rs` has a halving truncation check, and the piano has adaptive Ritz. Single-part prestressed tangent modes exist: membrane `solve_vibration`, shell `equilibrate_shell_materials`, `fs-solid::prestress`.

**Target.**
- **Prestressed equilibrium first.** The assembled prestressed equilibrium is solved before tangent modes or playback. String tension, bridge downbearing, soundpost compression, neck relief, supports and existing contacts all take part in the same equilibrium (codex review; Gough 2015).
- **Eigen route.** **Monolithic certified slicing of the assembled body** (M1.5 nested dissection) is primary.
- **CMS for parametric studies.** Component-mode synthesis is for parametric studies: Craig–Bampton with characteristic constraint modes for glue lines (Castanier, Tan & Pierre 2001), or free-interface methods with residual flexibility (MacNeal 1971; Rubin 1975), with certified truncation ladders.
- **Frequency-dependent damping** is preserved through any reduction.

**Validation.**
- Assembled vs monolithic eigenvalues.
- Violin and guitar body modes against published measurements (Bissinger; Strad3D; Carcagno in the corpus).

**Work items:**
- **M1.4-W1:** prestressed assembled equilibrium.
- **M1.4-W2:** monolithic certified eigen of the assembled body.
- **M1.4-W3:** CMS with characteristic constraint modes and ladders.

#### M1.5 Large and nonlinear eigenproblems (fs-sparse, fs-modal)

**Purpose.** Modes of 10⁵–10⁶-DOF bodies, cavities and coupled systems: any structural or acoustic eigenproblem.

**Now.** `fs-modal` needs real symmetric K and SPD M. `fs-sparse` has AMD ordering only.

**Target.**
- Nested-dissection ordering (George 1973).
- Parallel spectral window slicing.
- Nonlinear complex-symmetric eigenproblems (BEM-loaded and boundary-layer-damped modes) by Beyn contour integrals (Beyn 2012) or NLEIGS (Güttel et al. 2014), certified by residuals.

**Work items:**
- **M1.5-W1:** nested dissection plus window slicing.
- **M1.5-W2:** nonlinear eigen (Beyn/NLEIGS) with residual certificates.

### S2 — Acoustic fields

#### M2.1 Ducts in the frequency domain (`fs-duct`)

**Purpose.** Acoustic networks of any geometry: wind instruments, exhausts, ventilation, vocal tracts.

**Now.**
- **Plane-wave TMM** for cylinders and linear cones, with Lossless, WideTube, AllRegime and Bessel loss options. Bessel is exact from `GasState`. AllRegime's narrow-tube branch is a hand interpolation (`lib.rs:525–530`).
- **Lossy cones** in 2–12 slices with no ladder.
- **Tone holes** are compact T-junctions (Dalmont, Nederveen, plus an authored mutual term).
- **Terminations:** Unflanged refuses ka > 1 (`lib.rs:1103`); Flanged uses an 8-ring piston (`lib.rs:1111–1121`).
- **Multimodal:** m = 0, `MAX_MODES = 8`, staircased, holes refused.
- **Missing:** pads, undercut, external hole interaction, nonlinear edge loss, smooth horns and moving geometry.

**Target physics.**
- **Propagation.**
  - **Continuous multimodal admittance-Riccati propagation** on smooth profiles (Pagneux, Amir & Kergomard 1996; Amir et al. 1997), with refinable radial bases and **azimuthal m ≠ 0 bases** for asymmetric features. This replaces staircasing, which converges only at O(Δx) because of mode conversion at the steps.
  - **Plane-wave Riccati/Webster** is the lowest rung; refining its integration step converges to Webster's equation, not to 3-D.
  - **Stable composition.** Scattering/impedance composition is used where transfer matrices become ill-conditioned.
- **Loss.** **Bessel viscothermal loss** applies in admissible cylindrical linear sections, and resolved thermoviscous enrichment (M2.4) elsewhere.
- **Features** as local 3-D subdomains (M2.4) embedded through multiport transfer matrices **including evanescent modes**:
  - tone holes with chimney, undercut, pads/keys above, closed-hole cavities, and the **fingers above open holes**;
  - register vents, mouthpiece chambers, bends.

  Benchmark: Lefebvre & Scavone (2012).
- **Moving pads.** Pads are moving boundaries driven by **key-mechanism dynamics** (spring, pad, cork/felt bumper contact), so key clicks and pad slaps emerge.
- **Terminations.**
  - **Levine–Schwinger** for thin-walled circular pipe ends (organ pipes; a benchmark).
  - **Bells and mouths:** tabulated terminations baked from axisymmetric FE+PML or BEM (M2.4-W1d/M2.5-W1) across the full band.
- **Coupling and fields.**
  - **Open-hole interaction** through the multi-aperture radiation matrix (M2.5).
  - **Generic side-branch networks** (nasal tract, organ-pipe chests).
  - **Continuous gas fields** (M0.1-W3).
  - **Wall vibration** ports to shells (M1.2).

**Validation.**
- Closed-form cylinder and cone.
- Levine–Schwinger tables.
- **Same-specimen measured input impedance** (§1.3): the UNSW instrument description, or an owner-measured specimen (OD4). The Ernoult cylinder is already in the corpus.
- Pad effects (Dalmont; Lefebvre).

**Consumers.** All winds, voice tract, exhausts, HVAC.

**Work items:**
- **Terminations and loss:**
  - **M2.1-W1:** Levine–Schwinger thin-walled termination.
  - **M2.1-W2:** baked full-band terminations in `driving_point` and characteristic lines.
  - **M2.1-W5:** exact Bessel everywhere (replace AllRegime).
- **Propagation:**
  - **M2.1-W3a:** plane-wave Riccati for smooth profiles plus a cone slice ladder.
  - **M2.1-W3b:** continuous multimodal Riccati, m = 0.
  - **M2.1-W3c:** azimuthal m ≠ 0 bases.
  - **M2.1-W9:** stable scattering composition.
- **Features:**
  - **M2.1-W4a:** feature embedding via M2.4 multiports, with evanescent modes and fingers above holes.
  - **M2.1-W4b:** key-mechanism dynamics and moving pads.
  - **M2.1-W6:** multimodal holes, higher-mode losses and mouth radiation matrix.
- **Networks and gas fields:**
  - **M2.1-W7:** continuous gas fields.
  - **M2.1-W8:** generic side-branch networks.

#### M2.2 Time-domain ducts (DG with exact viscothermal operators)

**Purpose.** Run ducts in the time domain inside self-sustained loops, with moving geometry and nonlinearity. Any duct network qualifies (exhaust transients, HVAC, speech).

**Now.** Five realizations exist:
- **(a) Characteristic lines** from TMM reflectance with a default 4-round-trip FFT, which wraps around (DEF-02). They are used by every reed and wind path.
- **(b) `impedance_line`.** Its convention is unadjudicated (DEF-01), and it forces Z(0) = 0 (DEF-14).
- **(c) An 8-cell LC ladder** with 3 Foster branches as the default reed path in `realize_assembly`. This is a reduced engine.
- **(d) A cylinders-only waveguide network** with a fixed 8-term √ω Foster fit.
- **(e) Crossfaded per-fingering FIR banks.**

**Target.**
1. **Characteristic lines** with convergence-driven FFT sizing (a length ladder until the reflectance tail is below tolerance). Static compliance, undamped poles and propagation delays are realized as explicit states, separate from the decaying convolution kernels.
2. **A high-order 1-D DG/FEM scheme** with **passive diffusive (pole) approximations of the exact Zwikker–Kosten operators**, using a pole ladder (Bilbao & Harrison 2016; Thibault & Chabassier 2021; Tournemenne & Chabassier 2019 for convergence). Webster–Lokshin is kept only as a rung. The scheme is energy-consistent, and natively supports:
   - moving geometry (time-varying area and hole opening) by a conservative ALE formulation with the geometric conservation law (Farhat, Geuzaine & Grandmont 2001), including transfer rules for diffusive memory variables when area or hole state changes;
   - continuous gas fields;
   - wall vibration;
   - side branches;
   - nonlinear propagation (M2.3).
3. **Multimodal time-domain extension** for flares, from M2.1-W3b.
4. **The LC-ladder default is removed.**

**Validation.**
- Impulse responses against analytic cylinder and cone solutions.
- Simulated input impedance against M2.1 Bessel within 1 cent at peaks, including narrow tubes.
- Energy balance with moving geometry.
- A slow-opening hole reproduces the quasi-static impedance.

**Work items:**
- **Characteristic lines:**
  - **M2.2-W1:** adjudicate the `impedance_line`/`characteristic_line` conventions with causality tests (energy at t < 0; analytic delayed-wave and closed-volume compliance cases). Change code only if a test fails.
  - **M2.2-W2:** automatic FFT sizing plus length ladder, with compliance/pole/delay states (DEF-02, DEF-14).
- **DG scheme:**
  - **M2.2-W3a:** lossless variable-area DG.
  - **M2.2-W3b:** passive ZK pole approximations.
  - **M2.2-W3c:** side branches plus baked radiation boundary.
  - **M2.2-W3d:** conservative moving geometry and memory transfer.
  - **M2.2-W3e:** valve-port coupling.
  - **M2.2-W3f:** multimodal time-domain propagation.
    - It is a causal, passive multiport realization of the M2.1 multimodal operator, keeping propagating and evanescent modes, losses and port work.
    - It defines the interface to the nonlinear 1-D bore and bounds that interface's error.
    - It is required wherever the validity audit asks for higher modes (bells, flares, the clarinet bell above about 3.4 kHz).
- **M2.2-W4:** remove the LC-ladder default.

#### M2.3 Nonlinear acoustic propagation

**Purpose.** Shock formation in bores and ducts: brass at ff, exhausts, high-intensity ducts.

**Now.** None. DECISION_MUSIC_MEAN_FLOW.md dismissed it on a calculation that fails at ff.

**Target.** Bidirectional nonlinear propagation inside the M2.2 DG, coupled to excitation and radiation:
- entropy-stable fluxes;
- subcell shock capturing (Gassner 2013), with numerical dissipation shown to be ≪ physical dissipation at mf;
- viscothermal loss.

The one-way generalized-Burgers marcher is a **diagnostic only**.

**Validation.**
- Burgers analytic solutions and shock-formation distance.
- Trombone and trumpet ff spectral enrichment (Hirschberg et al. 1996, trombone; Myers et al. 2012 brassiness).

**Work items:**
- **M2.3-W1:** one-way Burgers diagnostic (X-Br-C1).
- **M2.3-W2:** entropy-stable nonlinear flux with shock capturing in the DG.

#### M2.4 Three-dimensional interior acoustics

**Purpose.** Acoustic fields in complex cavities and junctions: mouthpieces, tone-hole junctions, bends, valve ports, instrument body air (violin/guitar with f-holes/soundholes), the hand in a horn bell, resonator tubes, car cabins.

**Now.** No 3-D interior acoustics.
- **Lumped models only.** All are uniform-pressure, single-neck and fixed-geometry, compute no refusal band, and are natively unverified:
  - `cavity_modal.rs`: fs-phs Helmholtz (DEF-09) lowered into `CoupledModalSystem`; this is the string-body path;
  - `cavity_modal/shared.rs:47` (`SharedCavity`): N moving surfaces on one compression coordinate; used by `textile_cavity`;
  - `cavity_phs::piston`: nonlinear closed gas parcels, double-acting pistons, finite-body heat.

  These are model-form rungs, never acceptance.
- **Axisymmetric FE:** an axisymmetric lossless cavity FE exists (`meridian.rs`, kettle only).
- **Bead 0ja4** planned a new fs-acoustics crate that would depend on fs-couple, which would be cyclic.

**Target.**
- **Interior FE.**
  - **P2 (or higher) H1 Helmholtz elements** with viscothermal boundary-layer admittance (Bossart, Joly & Bruneau 2003).
  - **Thin-gap linearized Navier–Stokes** (or sequential LNS; Kampinga 2010) for gaps, slits and pads.
  - **Axisymmetric Fourier-azimuthal FE with PML** for bells and mouths, seeded from `meridian.rs`.
- **Coupling** to structures (M1 skins), exterior radiation (M2.5) and tissue obstacles (a hand in a bell).
- **Extraction** of:
  - multiport transfer matrices with evanescent modes, for embedding into ducts;
  - modal cavity bases (via M1.5) for vibroacoustics.
- **Home:** decided by M2.4-W0, with an fs-feec extension preferred. It is placed below fs-couple and supersedes 0ja4's interior scope (owner confirms, OD10).

**Validation.**
- Rigid-box modes.
- Helmholtz resonator.
- Damped-cylinder impedance against Bessel.
- Narrow-slit damping.
- Tone-hole parameters (Lefebvre FEM).
- Violin A0/A1 modes.

**Work items:**
- **M2.4-W0:** home decision and design.
- **Interior FE:**
  - **M2.4-W1a:** P2 Helmholtz FE.
  - **M2.4-W1b:** boundary-layer admittance.
  - **M2.4-W1c:** thin-gap LNS.
  - **M2.4-W1d:** axisymmetric FE+PML.
- **Extraction and coupling:**
  - **M2.4-W2:** multiport extraction with evanescent modes.
  - **M2.4-W3a:** structure–cavity coupling.
  - **M2.4-W3b:** aperture-to-exterior coupling.

#### M2.5 Exterior radiation (`fs-bem`, `fs-fmm`)

**Purpose.** Sound radiated by vibrating bodies and openings into free space, for any radiator.

**Now.**
- **BEM solver:** Burton–Miller with constant panels, dense LU, `MAX_DENSE_PANELS = 8192` (`helmholtz.rs:82`), and rigid/velocity BCs.
- **Products:** radiation-impedance matrices, SH directivity, near-field and first-order receivers, and reciprocal Galerkin radiation (restricted band).
- **Bell bakes:** lathed solid bells, 150–1446 Hz.
- **Gaps:**
  - fluid loading only in examples, plus planar baffled Rayleigh modal loading with a positive-semidefinite resistance (`rayleigh_modal.rs:175`, 973f27f43; used by the piano example);
  - in-loop compact first-order baffled loads per outlet (`radiation.rs`; ka ≤ 0.5, no mutual reaction; 6113eca8c). These are not a W3b realization;
  - three piston copies;
  - `fs-fmm` is a Chebyshev FMM with a real-valued kernel interface (`Kernel::eval -> f64`, `lib.rs:159`).

**Target.**
1. **Boundary types:** open thin-shell radiators and impedance BCs.
2. **Multi-aperture radiation matrices** (self + mutual) for holes plus bell, used **inside the oscillation loop**. Compact hole arrays use a retarded-multipole ladder (Lefebvre, Scavone & Kergomard 2013).
3. **Fluid–structure loading** for every shell and solid body.
4. **Fast solvers:**
   - complex kernels in fs-fmm, then the **directional Chebyshev FMM** (Messner, Schanz & Darve 2012) for high kD, which reaches about 130 for a violin at 20 kHz and about 370 for a tuba;
   - ACA/H-matrix only where kD is small;
   - GMRES with a preconditioner.

   Scheduling follows measured matrix size and kD.
5. **Operator realization.**
   - **Observation transfers** (body → microphone, outside the loop) are long causal FIRs with an FFT-length ladder; no passivity is needed.
   - **In-loop loading operators** use delay-extracted rational fits with Hamiltonian-matrix passivity tests and residue perturbation (Grivet-Talocia 2004; Gustavsen 2008), after a causality check of the BEM data (Triverio & Grivet-Talocia 2006). Convolution quadrature (Lubich 1988) is the reference.
   - **Closed-loop audit (joint acoustic-realization budget).** A passive local fit percentage does not bound lightly damped closed-loop poles, quiet notches or nonlinear thresholds.
     - Each realization is audited after the loop is closed: complex response, loaded modes and decay, and sound over the admitted band and control domain.
     - D⁻¹ΔD sensitivity is diagnostic only. It does not certify coverage across frequencies or nonlinear regimes.
     - None of this loosens the full-band requirement.
6. **One closed-form piston** (Bessel–Struve).
7. **Unbaffled full-band bell/aperture bakes.**
8. **Scatterers:** the player's body and chair as generic obstacle meshes (Pätynen & Lokki 2010 directivity was measured with players present).

**Validation.**
- Pulsating and oscillating spheres.
- Baffled piston; Levine–Schwinger pipe.
- Multipoles.
- Measured instrument directivity (TU Berlin, licence pending; Meyer).

**Work items:**
- **Bakes, piston and boundary types:**
  - **M2.5-W1:** band-complete unbaffled bakes.
    - The band and geometry are declared before baking.
    - Hitting the dense-panel cap before convergence leaves the item blocked and routes it to M2.5-W5.
    - Partial-band artifacts carry an explicit band refusal.
  - **M2.5-W2:** one closed-form piston.
  - **M2.5-W7:** open thin shells and impedance BCs.
- **Multi-aperture radiation:**
  - **M2.5-W3a:** multi-aperture Z_rad matrix.
  - **M2.5-W3b:** passive delay-extracted rational realization.
  - **M2.5-W3c:** in-loop coupling (X-Cl-R2).
  - **M2.5-W8:** compact hole arrays.
  - **M2.5-W10:** solver-route qualification for open finite-wall tubes and resonators (bar–tube–exterior).
- **Loading and scatterers:**
  - **M2.5-W4a:** solid radiation loading.
  - **M2.5-W4b:** shell radiation loading.
  - **M2.5-W9:** player/chair scatterers.
- **Fast solvers:**
  - **M2.5-W5a:** complex kernels in fs-fmm.
  - **M2.5-W5b:** directional FMM.
  - **M2.5-W5c:** GMRES and preconditioner.
- **M2.5-W6:** causal-FIR observation transfers with a CQ reference.

#### M2.6 Microphones and the recording chain

**Purpose.** Observe fields as a recording would.

**Now.**
- **Receivers.** fs-bem near-field pressure and first-order (cardioid etc.) receivers exist but are not in the music path. The reed scheduler has two observation inputs (`render/schedule/reed/outlet.rs`):
  - a legacy `listener_m` composition (bore pressure plus extras, DEF-12);
  - a physical `baffled_outlet` receiver over `OutletPressure`, landed recently by the remote lane and not yet natively qualified (OD13). It admits only `Termination::FlangedOpen` (`outlet.rs:43`).
  - further physical-outlet adapters over the same `BaffledPressure` owner: the `TwoSidedPressure` environment (6113eca8c), and `AperturePerformance` with an observation-only exterior gas (cf8b1e2b0). None is natively qualified, and the reaction load and the receiver can see different media.
- **No chain emulation.** There is no recording-chain emulation.

**Target.**
- **Receivers** sampling p and ∇p from M2.5 fields, generalizing `FirstOrder`.
  - **Capsule response and pattern** as data. The Iowa MIS pages state:
    - **Until 2011:** one Neumann KM 84 cardioid condenser. The bassoon page, for example, gives 5 feet, a Mackie 1402-VLZ mixer, a Panasonic SV-3800 DAT and 16-bit/44.1 kHz mono.
    - **From 2011:** three Earthworks QTC-40 omnidirectional microphones in a Decca tree 5 feet in front of the player, recorded at 24-bit/96 kHz. The centre microphone was downsampled to 16-bit/44.1 kHz mono for the chromatic-scale files.

    So pre-2011 references need a cardioid (p and ∇p) receiver, and 2011+ references an omni pressure receiver. All instruments except the piano were recorded in the anechoic chamber of the Wendell Johnson Speech and Hearing Center.
  - **Proximity effect** emerges from near-field sampling.
  - **Player and chair** act as scatterers (M2.5-W9).
- **Ensembles.** Distance/angle ensembles come from a declared prior and are marginalized (§1.3).
- **Recording-chain emulation from measured floors:**
  - the session ambient files (e.g. `BbClar.ambient.silence.aiff`), and inter-note gaps minus breath and key events;
  - the graded files' quantization and sample rate (16-bit/44.1 kHz).
- **Flags and caveats.**
  - The chamber low cutoff is undocumented, so descriptors below 100 Hz are flagged until it is sourced.
  - The 2012 `.mono` files are the downsampled centre omni of the Decca tree (MIS post-2012 page). `.stereo` files are a different mix and are not used for grading until their derivation is documented.
- **Room acoustics** only for site playback.

**Work items:**
- **M2.6-W0:** qualify the existing `baffled_outlet` observer natively (receiver ladder, analytic baffled-piston check) as the default for **flanged development fixtures only**.
  - Real unflanged, belled and tone-hole instruments need the physically matching receiver and load route (M2.5, M2.6-W1). Termination or environment is never changed to make an observer admit a render.
  - Retire the legacy `listener_m` composition (DEF-12) from acceptance; keep it only as the X-Obs1 diagnostic. Coordinate with the remote lane.
- **M2.6-W1:** generic p/∇p receivers with capsule data.
- **M2.6-W2:** recording-chain emulation from measured floors.
- **M2.6-W3:** observation prior and marginalization.

### S3 — Flow and excitation

#### M3.1 Moving-wall flow interfaces (valves): a port-graph pattern, not a type

**Purpose.** Any flow-controlling deformable aperture: single and double reeds, brass lips, vocal folds, relief valves, heart valves.

**Now.**
- **Reed family (four owners):**
  - `DynamicAperture`: 1-DOF, lay contact, swept flow, Maxwell memory, inward-striking only; single and symmetric DKT plate reductions.
  - `PairedAperture`: independent blades, spatial and hereditary.
  - `ModalAperture`: 1–32 signed modes, each with its own swept area and slit projection, and one effective uniform slit contact. It is built from 1–2 material DKT plates with Kelvin–Voigt or full viscous C (d66962bfa, 4e3d1de55, 70ccffa3b, 39f15c202), and hosted by `LoadedAperture`/`TwoSidedAperture` with two-sided tract and bore networks (71bf3b183).
  - The CLI's `ReedBoreVoice`/`BeatingReed`.

  Only `ReedBoreVoice` and the `DynamicAperture` performers are reachable from the CLI. All of this is natively unverified.
- **Brass lip:** 1-DOF outward, one-sample lag, no collision.
- **Glottis:** two-mass, quasi-steady.
- **Common limits:**
  - midpoint linear stepping (DEF-04). `ModalAperture` admits ω·Δt ≤ 1, which is up to −131 c (`material.rs:121–122`);
  - quasi-steady √|Δp| flow (`moving.rs:13`);
  - one effective uniform slit contact;
  - no mode, h or Δt ladder (the tests use a 3-mode toy plate);
  - a work/flow audit that cannot fail for `relative_tolerance` ≥ 1/128 (`modal.rs:293–300`, DEF-05 class).

`fs-phs::bernoulli_volume_flow` (`lib.rs:4675`) is pure √|Δp|, singular at closure, with no viscous term (DEF-13).

**Target.** A valve is a **pattern assembled in the port graph**:
- **Wall:** an S1 continuum (plate, shell, or solid; multimode; distributed pressure traction; fluid loading).
- **Contact:** S4 contacts with lays, rims, other blades and tissue.
- **Gap-flow two-port** with a **flow-law hierarchy**.
  - **Rung A** is quasi-steady Bernoulli, plus the laminar channel term Δp = ρU²/(2(wh)²) + 12μLU/(wh³) (Hirschberg et al. 1990; van Zon et al. 1990), unsteady channel inertance, swept volume and contraction.
  - **Rung B** is resolved channel flow: an M3.2 zone, the **reference for acceptance** (OD11). It is a moving-wall coupling:
    - no-slip and traction exchange, swept mass and wall work close conservatively through opening, contact and reopening;
    - coupled oscillation and thresholds converge.

    Flow through a prescribed channel alone does not qualify.
- **Channel/staple ducts** with pressure recovery (M2.1/M2.4).
- **Player airway port upstream:** lungs, trachea and vocal tract as an M2.1/M2.2 duct with yielding walls. Its area functions and playing impedances come from player data (Fritz & Wolfe 2005; Chen, Smith & Wolfe 2008/2009; Scavone, Lefebvre & da Silva 2008; Fréour & Scavone 2013).
- **The bore downstream.**

Striking direction emerges from geometry, with any number of blades. The reed (cane) and the lips (layered active tissue with rim contact) are **two fixtures of one code path**.

**Numerics.** Implicit, energy-consistent coupling (§3.3). Ladders in modes, h and Δt. Regime maps (M5.6) and time-domain threshold sweeps.

**Validation.**
- Reeds: Dalmont et al. 2003 (quasistatic characteristic and thresholds); Almeida et al. 2007 (double reed with staple recovery).
- Lips: Bromage et al. 2010 (open area); Newton et al. 2008 (mechanical response).
- Thresholds: Wilson & Beavers 1974.
- Glottis: replicas (Ruty et al. 2007).

**Work items:**
- **Pattern (M3.1-W1a–d):**
  - **M3.1-W1a:** the pattern, with DynamicAperture as the first implementation (bitwise).
  - **M3.1-W1b:** multimode distributed traction.
  - **M3.1-W1c:** striking direction from geometry.
  - **M3.1-W1d:** N blades plus contact.
- **Other apertures and the flow law:**
  - **M3.1-W2:** brass lip and glottis onto the pattern (collision, swept flow).
  - **M3.1-W3:** rung-A flow law (laminar, inertance, contraction; DEF-13).
  - **M3.1-W4:** staple/channel networks.
- **Fixtures:**
  - **M3.1-W5:** valve continuum from any S1 description (cane fixture).
  - **M3.1-W6:** lips fixture (layered active tissue, rim contact).
  - **M3.1-W7:** player airway port.
- **M3.1-W8:** rung-B resolved channel flow.

#### M3.2 Resolved flow (`fs-lbm`)

**Purpose.** Resolved compressible flow for jets, edges and channels: aero-acoustics, HVAC, biomedical flows.

**Now.**
- **Lattices and collision models:** D2Q9/D3Q19 with BGK, regularized, RR3, TRT, central-moment and reduced cumulant collision models.
- **Boundaries:** sponge; moving and off-lattice walls.
- **Physics:** isothermal weakly compressible acoustics, a passive FV energy equation, and free surface.
- **Refinement:** 2-D only (`refine.rs`, Dupuis–Chopard).
- **Missing:** D3Q27, a resolved energy equation coupled to the flow, and validated jets at playing Re.

**Target.**
- **Qualified discretization.** **Resolved compressible flow with energy transport** is the reference (DNS at Re ≈ 1,000–3,000, which is laminar or transitional). LBM is qualified as a discretization of the compressible Navier–Stokes equations: acoustic dispersion, damping, thermal behaviour, wall geometry, boundary reflections.
- **D3Q27** for rotational invariance (White & Chong 2011).
- **Collision-operator dissipation as a ladder axis** (τ ≈ 0.5015 at Δx = 50 µm).
- **SGS models** only as rungs whose contribution must vanish under refinement.
- **Acoustically transparent 3-D multi-block refinement** (Gendre et al. 2017).
- **Stated cost:** Δt ≈ 84 ns, i.e. about 10⁷ steps per second of sound. This is campaign work on trj/ts1.

**Validation.**
- Round-jet self-similarity.
- Brown edge-tone stages.
- Verge recorder-like flues.
- Grid ladders.
- Comparison with compressible Navier–Stokes recorder references (Yokoyama et al. 2015; Giordano 2014).

**Work items:**
- **Lattice and qualification:**
  - **M3.2-W1:** D3Q27 plus a collision-dissipation ladder.
  - **M3.2-W3a:** coupled compressible Navier–Stokes–Fourier dynamics: equation of state, total energy, pressure work, viscous heating and heat flux. This needs a thermal LBM formulation or a hybrid finite-volume energy solver; D3Q27 alone supplies none of them.
  - **M3.2-W3b:** qualification of W3a: linear acoustics, nonlinear manufactured solutions and thermal-flow cases, with finite-Mach and lattice defects quantified at fixed physical Re, Ma and Pr.
  - **M3.2-W5:** SGS as vanishing rungs (low priority).
- **Refinement:**
  - **M3.2-W2a:** 3-D multi-block refinement.
  - **M3.2-W2b:** transparent-boundary validation.
- **M3.2-W4:** jet/edge-tone validation campaign at playing Re.

#### M3.3 Zonal flow–acoustic coupling

**Purpose.** Couple a resolved-flow zone to acoustic fields: flues, reed and lip channels, glottis, ducts with fans.

**Now.** None.

**Target.**
- **Interface.** A resolved zone exchanges characteristic variables with a DG duct (M2.2) and the exterior (M2.5). The interface's location, transverse basis and projection quadrature are **refined independently**, and near-cutoff evanescent decay sets the distance. Measured reflection must be < −40 dB.
- **Multi-rate exchange** via M5.2.
- **Source ownership** (§1.2): resolved traces own the sound source.
- **Hydrodynamic/acoustic splitting** (Ewert & Schröder 2003) is only a model-form rung.

**Validation.**
- Pipe-mouth impedance.
- Recorder thresholds and regimes (Auvray et al. 2012; Terrien et al. 2014).
- DNS references (Giordano 2014; Yokoyama et al. 2015).

**Work items:**
- **M3.3-W1:** characteristic interface.
- **Recorder/organ-pipe campaign:**
  - **M3.3-W2a:** recorder/organ-pipe resolved zone with a grid ladder.
  - **M3.3-W2b:** threshold and regime sweep.
  - **M3.3-W2c:** DNS comparison.
- **M3.3-W3:** resolved aperture and edge zones for tone holes, register vents, bell rims and pipe/tube rims. They cover localized nonlinear losses and vortex shedding, which can affect register transitions. Zone extent is refined until omitted hydrodynamic structures and nonlinear propagation outside the zone are below the error budget.

#### M3.4 Accelerations from resolved simulation

**Purpose.** Faster validated surrogates for repeated runs: jet drive, contraction, flow noise. They **never qualify an acceptance render** (OD11).

**Now.**
- Jet cards covering Re 144–264.
- `NoiseTable` relative shapes, unconsumed.

**Target.**
- **Accelerations.** Jet-drive and labium accelerations over (Re, h/W, offset), and contraction Cc(Re, geometry), each minted from M3.2/M3.3 campaigns with error measured across the domain.
- **Flow noise** as **physically scaled volume-flow/dipole sources at the separation point inside the coupled step**, so the oscillation modulates them. Noise is never added to the output.
- **Validation of the noise** by cyclostationary coherence at k·f0 against references.

**Work items:**
- **M3.4-W1:** jet-drive acceleration.
- **M3.4-W2:** contraction acceleration.
- **M3.4-W3:** in-loop flow-noise acceleration with cyclostationary validation, for exploration only. Acceptance renders take their noise from the resolved owner (M3.1-W8, M3.3), must show the same cyclostationary statistics on the resolved output, and must disable every synthetic source that represents the same physics.

### S4 — Contact and friction

#### M4.1 Normal contact

**Purpose.** Any unilateral contact between bodies: frets, lays, mallets, fingers, supports, machine parts.

**Now.**
- **Laws and binding:**
  - Hertz, Hunt–Crossley and finite-gap BEM laws;
  - `fs-contact::interface_binding` with `fs-matdb::InterfaceSystemCard` (the Euler-disc program's generic binding);
  - power-law distributed obstacles in `fs-dcontact` (energy-exact penalty potentials, no LCP).
- **Defects and inconsistencies:**
  - obstacle friction uses `.unwrap_or(0.0)` (DEF-07);
  - weight conventions conflict;
  - the jawari pack is fixture-calibrated (BR-03).

**Target.**
- **Extend the existing binding** (no parallel trait).
- **Compliant contact with physically derived stiffness** (Hertz/footprint from geometry and M0 materials) through energy-consistent potentials is **primary**. A rigid nonsmooth formulation is only the limit rung (M4.3).
- **Refusals** propagate.
- **One weight convention.**
- **Laws** for felt (M0.4), tissue footprints (M0.3) and layered heads.

**Validation.**
- Hertz closed forms.
- Hunt–Crossley restitution, quantitatively.
- Footprint against BEM.

**Coordination.** Beads b8bxd.5.2 and .7.x (finite-patch laws into ports).

**Work items:**
- **M4.1-W1:** extended normal-contact binding with compliant potentials, the DEF-07 fix and one convention.

#### M4.2 Tangential interface laws (`fs-tribo`)

**Purpose.** Tribology generally: bowing, squeal, brake noise, stick–slip, fretting.

**Now.**
- **Friction laws:** Coulomb, velocity-dependent, Stribeck (DEF-06 jump), lumped partial slip, isothermal film.
- **Thermal:** flash temperature (no caller).
- **Name clash:** `fs_tribo::FrictionLaw` is an enum (`lib.rs:437`).

**Target.**
- **A stateful `TangentialLaw` trait.** Its state (slip, temperature, layer state) is port-Hamiltonian state with a dissipative port.
- **Render laws:**
  - **Cattaneo–Mindlin tangential compliance** (Mindlin 1949);
  - **thermal third-body friction**: friction set by the M0.4 layer's shear strength at the contact temperature, computed from spatial layer shear with moving-contact heat transport into string and hair. It is implemented together with M4.5.
- **Benchmarks only,** never on render paths:
  - Smith & Woodhouse (2000) lumped thermal model;
  - rate-and-state, LuGre, Stribeck.
- **Falsifier.** Galluzzo, Woodhouse & Mansour (2017) found that no friction law reproduces measured transients. If ribbon + thermal layer + measured player fluctuation also fails Bw2, the friction physics is still incomplete.

**Work items:**
- **M4.2-W1:** `TangentialLaw` trait plus the DEF-06 fix.
- **M4.2-W2:** Cattaneo–Mindlin.
- **Thermal third-body friction (with M4.5):**
  - **M4.2-W3a:** third-body layer in contact.
  - **M4.2-W3b:** moving-heat-source conduction into the bodies.
- **M4.2-W4:** benchmark baselines (rate-and-state, LuGre), test-only.

#### M4.3 Nonsmooth limit rung

**Purpose.** The rigid limit of contact and stick–slip, for any coupled system.

**Now.**
- Step-averaged set-valued projection and block Gauss–Seidel (≤ 128 sweeps).
- Hand-typed capture thresholds (B-04).
- Bead oh0i plans a Moreau–Jean lane (fs-mbd stepper, fs-contact cones).

**Target.**
- **Nonsmooth generalized-α** (Chen, Acary, Virlez & Brüls 2013) as a **diagnostic comparator only**, delivered through oh0i's ownership.
  - Like Moreau–Jean, it requires an explicit impact law (a Newton restitution coefficient) and retains controllable numerical damping, so it never renders acceptance audio.
  - Any rigid-limit impact response is derived from the compliant material and contact reference, with the limiting assumptions disclosed.
- **A disclosed root-selection rule** for non-monotone friction (McIntyre, Schumacher & Woodhouse 1983 hysteresis). The bow-point admittance does not depend on Δt, so the three roots survive refinement.
- **Δt ladders.**
- **Reference implementation:** Desvages & Bilbao (2016).

**Work items:**
- **M4.3-W1:** nonsmooth generalized-α via oh0i.
- **M4.3-W2:** root-selection disclosure, Δt ladders and stick–slip fixture.

#### M4.4 Driven bodies (one generic interface; exciters are descriptions)

**Purpose.** Bodies whose drive follows a performance: bows, fingers, plectra, mallets, hammers, tongues, and robot end-effectors.

**Now.**
- **Bowed strings:** point bow stations with prescribed loads.
- **Plucks:** initial shapes or released forces.
- **Mallets and hammers:** mallets wired to drums only; three piano felt adapters.

**Target.**
- **One generic interface.** Any S1 body (rod, ribbon, plate, solid) has drive DOFs that follow a Performance through a **declared compliant drive**, with player-hand impedance as data, and contacts other bodies through M4.1/M4.2. Incompatible effort and motion are never prescribed at the same boundary.
- **Exciters as fixtures** (descriptions plus validation):
  - **bow:** stick as an M1.1 rod with real camber/taper, plus hair ribbon as a distributed tensioned ribbon in **distributed line contact** with the string. Slipping varies across the ribbon width (Pitteroff & Woodhouse 1998).
  - **finger:** tissue solid plus nail.
  - **plectrum.**
  - **mallets:** layered heads (core plus rubber/yarn/cord wraps), handle flexibility, bounce, dead-stroke.
  - **hammer:** one unified felt hammer.
  - **tongue.**

**Validation.**
- Hair compliance (Pitteroff & Woodhouse).
- Bowed attacks (Guettler 2002; Galluzzo & Woodhouse 2014).
- Guitar plucks (Pavlidou 1997).
- Mallet contact times (Chaigne & Doutaut 1997).

**Work items:**
- **M4.4-W1:** driven-body interface.
- **M4.4-W2:** ribbon–rod distributed line contact.
- **Exciter fixtures (M4.4-W3a–f):**
  - **M4.4-W3a:** bow fixture.
  - **M4.4-W3b:** finger fixture.
  - **M4.4-W3c:** plectrum fixture.
  - **M4.4-W3d:** mallet fixtures.
  - **M4.4-W3e:** unified felt hammer.
  - **M4.4-W3f:** tongue fixture.

#### M4.5 Contact heat

**Now.** Thermal networks exist:
- fs-thermal `DepositionSimulation` (contacts, enthalpy–phase, advection);
- `fs-material::phase::UniformEnthalpyStep`;
- fs-couple `circuit_enthalpy`.

`cavity_phs` gas–wall heat uses authored effective distances (`thermal/layers.rs:46–49`). On acoustic paths, gas-side wall heat belongs to M2.4-W1b only, so it is counted once.

**Target.** Friction power is partitioned into the contacting bodies' thermal networks (the fs-thermal ones; no new store) and the third-body layer, and temperature feeds back into the friction and material laws. Rosin heat is implemented with M4.2-W3; general friction heat is wired here.

**Work items:**
- **M4.5-W1:** general friction-heat partition plus flash-temperature wiring.

### S5 — Coupling and runtime

#### M5.1 Port graph (fs-ports schema + solver in fs-couple)

**Purpose.** Compose any subsystems through power ports with exact energy accounting.

**Now.**
- **fs-phs (L2)** has dense, pairwise, hand-wired descriptors, plus the `prepared` path (analytic tangents, condensation).
- **Newton acceptance is relaxed implicitly** (`lib.rs:878, 1154`; `prepared.rs:376`; DEF-05).
- **Unwired:** `fs-ports` schema/audits and IQN-ILS.
- **Hand-wired N-way couplings** to subsume with parity:
  - `CoupledModalSystem` column connections;
  - `SharedCavity`;
  - `DoubleActingGasPiston` and its thermal/boundary variants (N = 2 hard-coded; nested scalar solves at `thermal/double_acting.rs:113–114`).

**Target.**
- **(a) Early schema adoption.** Every module exposes ports through `fs_ports::PortSchema` early (W0), so modules can adopt the schema before the solver exists.
- **(b) The solver in fs-couple.** The graph solver builds on `fs-phs::prepared`:
  - block residual/operator interfaces, matrix-free solves, local factorizations and Schur complements, with no flattening of continua into one dense pHS system;
  - **modulated interconnection** and mortar (§3.2);
  - nonsmooth nodes from fs-mbd/fs-contact;
  - refusal certificates;
  - monolithic or partitioned execution (IQN-ILS) chosen from measured coupling strength. IQN-ILS convergence is not a passivity or accuracy certificate.
- **(c) Integration, not a prerequisite.** The graph solver is integration; modules do not wait for it.

**Validation.**
- The graph reproduces existing hand-wired results bitwise.
- Partitioned results match monolithic within tolerance.
- The energy ledger closes to its budget.

**Work items:**
- **M5.1-W0:** port-schema adoption contract plus adapters.
- **Graph solver:**
  - **M5.1-W1a:** builder reproducing reed–bore bitwise.
  - **M5.1-W1b:** sparse block assembly with analytic tangents.
  - **M5.1-W1d:** IQN-ILS strategy.
  - **M5.1-W3:** modulated interconnection plus mortar.
- **M5.1-W2:** refusal certificates (DEF-05, all three sites).

#### M5.2 Integrators and multi-rate

**Purpose.** Accurate, structure-preserving time stepping across disparate clocks.

**Now.**
- **Exact modal ZOH:** narrow-band pressure transfer, exact only at modal frequencies.
- **Gonzalez midpoint everywhere else,** with warping (DEF-04).
- **Unwired:** `fs-time`.
- **Multi-rate** only as decimated observation. The live piano runs its whole engine at 1/2/4× through `pcm_wav::decimate::Decimator` (`engine-piano/src/audio.rs`; 2fa573a75); there is no subcycling.
- **Other midpoint and trapezoid paths:**
  - `CoupledModalSystem` connections are trapezoidal (`coupled.rs:5–14`; cavity compression, plate-viscosity dashpots).
  - `ImpactSystem` (`render/plate/impact.rs`) uses Gonzalez midpoint with only an ω·Δt < 0.9π guard (:276).
  - Piston walls are midpoint, and the network-coupled piston is first order.

**Target.**
- **§3.3 implemented across all render paths,** with a **warping audit**: every render path is listed with its integrator, and warping is measured.
- **Broadband modal pressure transfer.**
- **Energy-conserving multi-rate subcycling** (Prakash & Hjelmstad 2004) or a laddered TLM/wave-variable delay. Gravouil–Combescure dissipates at the interface and is not used.
- **Anti-aliased decimation.**
- **Clocks** (mechanics 192–768 kHz, acoustics 48–192 kHz, lattice MHz) chosen from measured phase, dissipation, contact and interface-error budgets. Force reconstruction and exchange order are refined, as well as step size.

**Work items:**
- **M5.2-W1:** integrator doctrine plus warping audit and fixes.
- **M5.2-W2:** broadband modal pressure transfer.
- **M5.2-W3:** energy-conserving multi-rate with decimation and admission rules.

#### M5.3 Configuration changes

**Purpose.** Change topology or geometry while running (fingerings, valves, slides, stops, contacts, terminations) in any switched system.

**Now.**
- **Configuration changes** are crossfaded frozen FIR banks.
- **Unconsumed gesture targets.** Nine typed gesture targets have no consumer: Fret, ValveInsertedLength, ToneHoleSigma, SlideLength, StringTension, HammerStrike, TerminationSwap, PreStress, JetAngle.

**Target.**
- **Continuous inputs.** Physical configuration parameters (hole opening, pad height, slide length, rotor angle) are continuous description inputs, evolved in the conservative moving-domain DG (M2.2-W3d).
- **Discrete topology events** use **conservative projection**, which retains:
  - trapped-volume pressure;
  - contact state;
  - constitutive memory variables;
  - heat;
  - moving-basis derivative terms.

  Projection error is quantified.
- **Gesture addressing.** Gestures address description parameters **by typed path**. The instrument-shaped targets become path aliases, not new consumers.

**Work items:**
- **M5.3-W1** is delivered by M2.2-W3d.
- **M5.3-W2a:** typed-path gesture addressing.
- **M5.3-W2b:** conservative topology events.

#### M5.4 Generic renderer

**Target.**
- `Renderer::from(description, performance)` built from the port graph, with generic controls.
- `PressureRenderer` extended with control input.
- One code path for offline and WASM.
- The CLI `music_render instrument <description> <performance> <out.wav>`.

Reachability through this renderer is the renderer's own DONE-WHEN, **not a gate** for instruments: Q0 accepts any committed CLI.

**Work items:**
- **M5.4-W1:** open voice trait plus generic controls.
- **M5.4-W2:** renderer from description plus CLI.

#### M5.5 Convergence-ladder harness

**Target.**
- Built on `fs-ladder` over **existing** render paths, with no dependency on the generic renderer.
- Every axis has its declared convergence/error model and coupled groups.
- Reports per §1.2:
  - ≥ 3 rungs;
  - observed order where the axis is algebraic, GCI only where its assumptions hold, the declared bound otherwise;
  - axis liveness from logged dimensions plus a sensitivity fixture;
  - a jointly refined run;
  - seed ensembles.
- Results are stored with the render's provenance.

**Work items:**
- **M5.5-W1:** harness.
- **M5.5-W2:** the §1.4 simplification-equivalence comparison on the same seam, with no new harness. It takes:
  - the detailed and the simplified descriptions;
  - a predeclared domain and bounds;
  - a raw calibrated comparison;

  and reports a bound on the difference.

#### M5.6 Regime analysis

**Target.**
- `fs_orbit::OrbitProblem` implemented for graphs inside fs-couple, giving thresholds, register/slot maps and Schelleng-type diagrams for periodic regimes.
- **Q1 is gated by time-domain sweeps on the render path.** Harmonic balance accelerates and cross-checks, and cannot see attacks or aperiodic regimes.

**Work items:**
- **M5.6-W1:** OrbitProblem for graphs.

### S6 — Descriptions, performance, data

#### M6.1 Assembly Description

**Purpose.** One typed, versioned description of any assembled physical system.

**Now.**
- About 55 per-feature formats.
- The `fs-scenario` IR entity model (Assembly → Part → Region/Surface/Interface) is unused by music.
- The `scenario.ir` v2 schema is frozen in `schema-policy.json`.

**Target.** An Assembly Description on the `fs-scenario` IR, introduced through schema-policy migration with a migration test. It contains:
- **parts:** strings, shells, solids, ducts, cavities, valve patterns, driven bodies;
- **geometry:** parametric generators, or imported meshes via fs-io, bound to the real source assembly and reconstruction uncertainty (§1.4);
- **material references:** spatial frames/layers and physical state, with unknown properties explicit;
- **interfaces:** ports, contacts, glue joints;
- **named configuration states → parameter assignments** (fingering, valve, stop, tuning);
- **observers;**
- **a validator** against module validity domains, run before rendering.

Only the formats the reference instruments need are converted, and inside their family description beads. Following AGENTS.md, bad APIs are not preserved through wrappers.

**Work items:**
- **M6.1-W1a:** schema plus migration.
- **M6.1-W1b:** validity validator.
- **M6.1-W1c:** first description, the R13 clarinet.
- **M6.1-W2:** generators (bore, outline with holes, arching).

#### M6.2 Performance

**Target.**
- **Commands vs physical states.** Player commands are distinguished from physical interaction states. Airway, embouchure contact, active tissue prestress and compliant holding are coupled parts (M3.1-W7, M1.2-W4, M0.3-W3), not prescriptions.
- **Score to gestures.** Score/MIDI becomes gesture tracks through dynamics tables (pp/mf/ff → pressure, bow force/speed, pluck force), articulation, and measured micro-variation.
- **The §1.3 per-note control rule** is implemented as a declared objective.

**Work items:**
- **M6.2-W1:** performance schema, MIDI import, dynamics tables and per-note control rule.
- **M6.2-W2:** player models from measured data.

#### M6.3 Geometry acquisition and reconstruction (any physical assembly)

**Purpose.** Turn real objects into solver-usable physical descriptions under §1.4 and the §8 geometry completion rule: instruments, players, rooms, machine parts.

**Now.**
- No acquisition or reconstruction workflow.
- `fs-io` reads meshes; gmsh/STEP paths are unwired to music. `fs-mesh` volumetrizes closed surfaces.
- Blender 5.3 is installed on the Mac (`/Applications/Blender.app/Contents/MacOS/Blender -b -P <script>`).

**Target.**
- **Source hierarchy per part,** recorded in the M8 ledger:
  1. Direct measurement of the specimen: calipers, bore and profile gauges, thickness gauges, casts of windways and channels, CT or structured-light scans through a service (OD15).
  2. Published measured geometry: CT studies, theses, maker drawings.
  3. Credible CAD or scans found online, checked against independent dimensions before use. 3-D-printing models are rarely dimensionally faithful.
  4. Blender reconstruction from calibrated photo references.
- **Blender reconstruction is scripted and reproducible.**
  - Inputs: photographs with scale references and calibrated cameras (photogrammetry such as COLMAP/Meshroom, or single-view calibration), plus a measured dimension table.
  - Scripts under `tools/geometry/blender/` build parametric models with `blender -b -P`, so a model regenerates deterministically from its inputs.
  - Exports are in SI units with hashes.
  - Every surface carries an uncertainty field marking observed against occluded geometry.
  - Hidden geometry comes from sections or CT, or stays a declared bound (§1.4).
- **Outputs.**
  - Source geometry, kept separate from numerical meshes.
  - Part, material, contact and support mapping.
  - Coordinate frames and material directions.
  - Uncertainty.

  These are consumed through fs-io and fs-mesh into M6.1 descriptions.
- **Uncertainty propagation.** Geometry uncertainty enters the simulation ensembles for population claims (§7.2) and bounds same-specimen claims.

**Validation.**
- Reconstruction of an object of independently known geometry (a machined test piece, a measured tube) within a declared tolerance.
- Caliper cross-checks on every reconstructed instrument.
- Geometry and field QA uses `apps/observatory` (VTU/VTR import, `MESH_FIELDS.md`); no new viewer.

**Work items:**
- **M6.3-W1:** geometry ledger rows plus **representation-aware** import checks through fs-io:
  - closed solids are watertight and manifold;
  - midsurfaces declare their boundary edges and thickness field;
  - real openings (tone holes, f-holes, soundholes, bell mouths, pipe ends) are preserved and never "repaired";
  - units, frames and part mapping are checked.
- **M6.3-W2:** Blender headless parametric reconstruction scripts with deterministic, hashed exports.
- **M6.3-W3:** photo calibration and per-surface uncertainty (photogrammetry route).
- **M6.3-W4:** selected-route qualification. The acquisition route a specimen actually uses (measured CAD, CT or scan import through W1, or Blender/photo reconstruction through W2/W3) is validated on a known object. Photo/Blender work is required only when that route is used.

#### M6.4 Player anatomy (interacting human bodies as description data)

**Purpose.** The parts of the human body that take part in playing, described as parts of the coupled physical experiment. The same description serves any human–machine interaction.
- mouth: lips, tongue, teeth and palate, cheeks;
- airway: oral cavity, vocal tract, trachea, lungs, nasal tract;
- hands and fingers;
- holding contacts (chin, shoulder, thigh, chest).

**Now.** Area-function tract tables only (Assaneo, `tract.rs`).

**Target.**
- **Representation.**
  - Lips and tongue are layered visco-hyperelastic solids with active prestress (M0.3, M1.3).
  - Teeth and palate are stiff boundaries.
  - The oral cavity and airway are M2.4 interiors or M2.1/M2.2 ducts with yielding walls.
  - Fingers and hands are tissue solids (M4.4-W3b).
  - Torso, thigh and chin contacts are compliant supports (M1.2-W4).
  - Each articulator has postures (configuration states) and contact surfaces against the instrument.
- **Geometry sources.**
  - Speech MRI: the USC 75-speaker database (Lim et al., *Sci. Data* 8, 187, 2021; CC BY 4.0). These are speech postures, not playing postures.
  - Real-time MRI of brass players (Göttingen group; access to verify).
  - Measured playing-tract impedances to validate the airway descriptions (Fritz & Wolfe 2005; Chen, Smith & Wolfe 2008–2010).
  - The owner's own anatomy if OD16 approves: an intraoral dental scan, a structured-light face and lip scan, MRI if available.
- **Rules.**
  - Playing postures are declared separately from speech postures, and the difference is an uncertainty.
  - Tissue properties come from M8-W2d, never from geometry.

**Validation.**
- Airway impedance computed from the described anatomy against measured playing-tract impedances (a population check).
- Passive geometry-to-acoustics check: the measured acoustics of the Dresden Vocal Tract Dataset's printed tracts.
- Lip open-area waveforms against Bromage et al. (2010).

**Work items:**
- **M6.4-W1:** anatomy parts in the existing Assembly Description schema (articulators, postures, contact surfaces), with no separate registry.
- **M6.4-W2:** airway from MRI segmentation, with the playing posture as declared uncertainty and an airway impedance check.
- **M6.4-W3:** lip and embouchure geometry and contact surfaces for reeds, brass mouthpieces, flute embouchures, pan-pipe rims and the duduk reed.
- **M6.4-W4:** hand and finger geometry for hole covering (including half-holing), stops and holding.

### S7 — Validation and grading

- **M7.1 Grader.**
  - **Location:** `apps/instrument-lab/tools/instrgrade/`. These are new files, and the site maintainer is notified. numpy is a dev-only dependency. `radiate.py` ships as DIAGNOSTIC-ONLY and is never applied to a graded render.
  - **Work items:**
    - M7.1-W1: commit the grader and reproduce the audit.
    - M7.1-W2: blind-spot fixes and descriptors (§7.2).
    - M7.1-W3: discrimination-equivalence listening harness (§7.1 Q4), reusing the listening-pack page.
    - M7.1-W4: per-family phrase fixtures and phrase grading (Q5), with frozen phrase references.
- **M7.2 Module conformance.** Each module's G1 and literature tests (Appendix D) are `#[test]` targets inside its work items.
- **M7.3 Passive-chain gates (QP).** Per family: input impedance per fingering, bridge admittance and signature modes, modal ratios. They are delivered inside the family milestones.
- **M7.4 Modal identification.** `fs-modalid` is wired into QP comparisons (M7.4-W1).
- **M7.5 Psychoacoustics.** `fs-psycho` loudness, sharpness, roughness and tonality become grader descriptors (M7.5-W1).

---

## 5. Consolidation, defects, breaches, registry obligations

### 5.1 Unify (each move or deletion is enumerated under OD10)

| Duplicate | Into | Item | Files moved or removed (OD10 list) |
|---|---|---|---|
| Strings in five places | fs-string | M1.1-W1a–d, M1.1-W2 | `fs-couple/src/modal_acoustic_time.rs` (→ acyclic owner); `string_specimen.rs`; `examples/grand_piano/{linear,dampers,hammer_footprint,string_stretching,string_damping,longitudinal_string,string_polarization,string_sections}.rs` (`git mv`; linear.rs `#[path]`-includes four of them, `:38–48`); for W2: `bowed_string/runtime/reciprocal/bank/**` and the string storage of `render/plate/impact/{string,supported,linear/wire}.rs` |
| Four reed valve owners, two condensation solves, three `with_plate_forces` | valve pattern | M3.1-W1a | `paired.rs`/`uniform.rs` condensation merged; one `with_plate_forces`; old entry points removed only after bitwise parity |
| Three physical-outlet adapters | one qualified `BaffledPressure` owner | M2.6-W0 | adapters stay thin over it |
| Three membrane-energy owners | fs-plate shell with D → 0 | M1.2-W8 | `fs-solid/src/membrane*` and the `shell::head::nonlinear` membrane folded in |
| Two lumped-cavity families on different gas owners | one cavity owner | M2.4-W3a | `cavity_modal` and `cavity_phs::piston` share one gas owner |
| Three valves | valve pattern | M3.1-W1a, W2 | glottis/lip loops become pattern instances; old entry points removed only after bitwise parity |
| Three felt adapters | unified hammer | M4.4-W3e | piano and drum felt adapters |
| Three piston quadratures | one closed form | M2.5-W2 | `fs-bem/helmholtz.rs:1484`, `fs-phs/lib.rs:1707`, `fs-couple/thin_plate.rs:1748` |
| Three modal viscous realizations | `plate_viscosity.rs` | M0.2-W3 | `thin_plate/material_loss.rs` diagonal image and `bernoulli_aperture/modal/viscosity.rs` routed through it; `plane_stress_viscosity.rs` stays the law owner |
| LC-ladder reed default | DG/characteristic line | M2.2-W4 | default path only |
| Closed `RenderVoice`, fixed-slot `AcousticAssembly` | generic renderer | M5.4-W1/W2 | — |
| Three WASM engines | `engine-perform` | L-2 | IRP families only. The piano, kit and drum engines migrate by their owners, with kernel parity at a fixed rung, never by adopting budgets or grids |

### 5.2 Defects (Appendix B gives file:line)

| ID | Defect | Item |
|---|---|---|
| DEF-01 | `impedance_line` convention unadjudicated | M2.2-W1 |
| DEF-02 | 4-round-trip FFT wrap-around | M2.2-W2 |
| DEF-03 | global-gain passivity on a sampled grid | M2.5-W3b; M2.2-W2 |
| DEF-04 | integrator frequency warping across fs-phs linear parts and the KC path | M5.2-W1; M1.1-W3 |
| DEF-05 | implicitly relaxed, inconsistently scaled Newton acceptance at three sites | M5.1-W2 |
| DEF-06 | traction ramp discontinuity | M4.2-W1 |
| DEF-07 | swallowed refusals | M4.1-W1 |
| DEF-08 | P1 tets in bending | M1.3-W1a |
| DEF-09 | Helmholtz neck end correction and loss (now load-bearing: cavity_modal, SharedCavity) | DEF-09 row (M) |
| DEF-10 | `material_strings.rs` compile | F1 |
| DEF-11 | stale CONTRACTs | updated inside each item |
| DEF-12 | reed observer composition | M2.6-W0 |
| DEF-13 | singular, inviscid gap flow | M3.1-W3 |
| DEF-14 | Z(0) = 0 erases compliance | M2.2-W2 |
| DEF-15 | grader blind spots | M7.1-W2 |

### 5.3 Breaches of §1.3

| ID | Breach | Remedy |
|---|---|---|
| BR-01 | Piano axial projection 0.06 (authored, introduced under a recording-targeted DONE-WHEN) | M1.1-W7 |
| BR-02 | fs-vkplate tensions and T60 fits tuned to recordings | recorded on the existing drum-kit bead frankensim-0nxhr (outside this plan's instrument scope); replace with M0.2 losses and measured head tensions |
| BR-03 | Jawari pack calibrated to a fixture | BR-03 bead: measured data or refusal |
| BR-04 | Hand-typed constants on render paths | each replaced by derivation or pack data inside its module item, or kept as a named test fixture; a BR-04 tracking bead lists them |

### 5.4 Registry and golden obligations (same commit as the change)

- **instrument-claims.json:** no new rows. Rows citing moved or removed code get updated evidence refs, or are transitioned to refused; doctrine D21 says rows are never deleted.
- **consolidation-review.json:** refresh the dispositions of newly wired crates (fs-psycho, fs-modalid, fs-thermochem, fs-aeroac, fs-contact).
- **doc-facts:** update for new crates.
- **Manifest:** regenerate the source manifest last.
- **music-baselines.json** WAV hashes follow the golden-bump protocol.
- **Frozen schemas:** any `scenario.ir` or `matdb.pack` change goes through schema-policy migration.
- **Grader evidence:** `data/instruments/grades/<family>/<run-id>.json`.

### 5.5 Verification debt

- **F1** natively builds and runs the instrument targets the first milestones consume.
- **Every work item** includes native execution of its tests, with host, command and exit code in the close reason.
- **OD3** forbids pushing instrument code that was not natively built and tested.
- **OD13** asks the owner to direct the remote no-toolchain lane.

---

## 6. Instrument compositions

Each family lists:
- description data;
- the module chain;
- named experiments X-…, each with its falsifier in brackets (each is an S bead in Appendix E);
- milestones and their dependencies (Appendix E).

### 6.0 First audible milestone M-CL0 (target ≤ 2 weeks)

**Purpose.** Measure, on existing paths, how close a real-dimension clarinet gets. This produces data, not acceptance.

The initial reduced bore/reed experiment does not waive §1.4: missing 3-D details and player anatomy remain explicit model-form limitations until their effects are resolved or bounded.

**Inputs.**
- F1 (native build of the paths);
- M7.1-W1 (grader);
- M2.6-W0 (radiated observation, DEF-12);
- M2.2-W1/W2 (convention adjudication and FFT sizing);
- M8-W3a (R13 + B4 bore, tone-hole lattice, fingering files);
- the existing reed (`DynamicAperture`) and the best existing termination (flanged or BEM-baked), recorded as a model-form limitation. The compact in-loop baffled load is valid only for ka ≤ 0.5. M-CL0 therefore uses the flanged TMM terminal with the qualified `baffled_outlet` observer, and records that.

**DONE-WHEN.**
- R13 fingerings over the Iowa range are rendered natively at a virtual microphone.
- They are graded with null baselines and the grader's blind-spot report.
- A scorecard names the dominant residual.
- Before/after JSON is committed to `data/instruments/grades/clarinet/`.

**Matched-observer diagnostic (X-Obs1).** It attributes the baseline gap between observation and physics.
- **Design.** One accepted reed–bore trajectory is observed by two operators: the legacy internal composition and the physical baffled outlet. Geometry, mechanism and gas are held fixed.
- **One factor at a time.** Factors are varied singly: observer, bore geometry (7.3 mm tube vs R13), tone-hole operator. Geometry and the tone-hole operator are never varied together.
- **No cross-instrument arm.** There is no "internal audio" comparison across instruments (the piano's internal quantities are unlike).

**Expected outcome.** Clarinet recognized on most notes, with the residual pointing at the lattice/radiation (X-Cl-R2) and noise/onset (X-Cl-R3).

### 6.1 Clarinet

**Description.** The Buffet R13 with a Glenn Bowen B4 mouthpiece, i.e. the UIowa reference:
- mouthpiece cavity, baffle and lay profile;
- barrel, joints and bell;
- the full lattice of about 24 holes: chimneys, undercuts, pads, pad heights, key geometry and springs;
- register key and fingering chart;
- cane reed geometry (vamp thickness map), with **wet** cane material;
- lip and tongue tissue;
- player data: blowing pressure, lip force, tonguing, and **the playing airway** (vocal tract).

**Chain.**

```
player airway ↔ M3.1 valve (reed continuum, lay + lip contact, swept flow, rung-A flow law;
                rung-B channel flow for acceptance)
              ↔ M2.2 bore (characteristic lines first, then the DG) + M2.1/M2.4 features
                (tone holes, pads, fingers above holes, key mechanisms)
              ↔ M2.5 multi-aperture radiation in the loop (bell bake)
              → M2.6 receiver
```

Supporting modules:
- M0.1-W3 steady breath field (for Q3);
- M3.4-W3 in-loop flow noise;
- M4.4-W3f tongue;
- M5.3 fingering changes.

**Experiments.**
- **X-Cl-R0 (QP):** input impedance of a description of the **same physical instrument that was measured**, under its measurement conditions (the UNSW instrument, or an owner-measured R13 under OD4). [Falsified if peaks are > 10 c off on a block of fingerings.]
- **X-Cl-R1:** pp→ff growth from flow nonlinearity (Worman) and from beating vs reed dynamics. [Falsified if a massless-reed pressure sweep already reproduces the per-harmonic slopes below 2 kHz.]
- **X-Cl-R2:** lattice cutoff and multi-port radiation explain the loss of odd/even contrast above about 1.5–3 kHz. [Falsified if radiating holes keep 10–20 dB contrast above 2 kHz.]
- **X-Cl-R3:** in-loop flow noise, airway and tonguing dominate residual noise and onset. [Falsified if output-added noise of equal spectrum matches the references' cyclostationary coherence as well as in-loop noise does.]

**Milestones.**
- **M-CL-Q2:** M-CL0, M2.1-W1/W2/W5, M3.1-W1a/W3/W5/W7, M2.5-W1/W3c, M6.2-W1, M7.1-W2.
- **M-CL-Q3:** adds M2.1-W4a/b, M2.2-W3c/W3f, M3.3-W3, M3.1-W8 (the resolved noise owner), M0.1-W3 + M2.1-W7, M6.4-W2/W3.
- **M-CL-Q4:** adds M7.1-W3 and the listening session (OD8).

### 6.2 Saxophone (soprano first, then alto)

**Description.**
- measured conical bore and keyed tone holes with pads;
- mouthpiece chamber, baffle and lay;
- cane reed;
- the player airway, which is essential for altissimo (Chen, Smith & Wolfe 2008; Scavone et al. 2008).

Verify that UIowa has soprano saxophone files; the audit fetched only alto.

**Chain.** As for the clarinet, with conical DG/TMM segments. The mouthpiece volume replaces the missing cone apex.

**Experiments.**
- **X-Dr0 (QP):** same-specimen impedance (UNSW soprano or owner-measured).
- **X-Sx1:** airway resonance moves altissimo register selection [falsified if the airway impedance has no effect above note C6].

**Milestones.**
- **M-SX-Q2:** M-CL-Q2 modules plus M2.1-W3a.
- **M-SX-Q4.**

### 6.3 Oboe and bassoon

**Description.** Measured bores; bocal (bassoon); double-reed blade geometry; staple.

**Chain.** As for the clarinet, with the M3.1 two-blade valve and staple pressure recovery (M3.1-W1d, M3.1-W4, M2.4-W2).

**Experiments.**
- **X-Dr1:** double-reed quasistatic characteristic against Almeida et al. (2007), including staple recovery.
- **X-Dr2:** bore vs reed as the dominant residual.

**Milestones.**
- **M-DR-Q2:** M-CL-Q2 modules plus M3.1-W1d/W4, M2.4-W2.
- **M-DR-Q4.**

### 6.4 Brass (trombone and trumpet first, then horn, then tuba)

**Description.**
- **Trumpet:** Bach 180S37 (the UIowa reference): mouthpiece cup/throat/backbore, leadpipe, tuning slide, valves and slides with toroidal bends and ports, bell flare, bell wall thickness/alloy/temper, and the player's hold.
- **Trombone:** tenor trombone with slide.
- **Horn:** **the hand is in the bell on every note.** This is normal playing, not stopping; it lowers the cutoff and moves the upper resonances (Dell et al. 2015).
- **Tuba.**
- **Lips** as a layered active-tissue continuum, with player airway.

**Chain.**

```
player airway ↔ M3.1 lips (rim contact, collision, active prestress, rung-A/B flow)
              ↔ M2.2 DG bore with M2.3 nonlinear propagation + M2.1-W3b multimodal flare
                (M2.4 bends, valve ports, hand in bell)
              ↔ M2.5 bell radiation (open thin shell, directivity) + M1.2 bell-wall vibration
                (a few dB; Kausel, Zietlow & Moore 2010)
              → receiver
```

M5.3 handles valve and slide changes.

**Tuba.** It is graded from the 2nd harmonic upward until the chamber cutoff is sourced (§2.6).

**Experiments.**
- **X-Br-C1:** a one-way Burgers diagnostic at ff reproduces the measured 4–8 kHz growth (+7.7 dB relative to the fundamental, mf→ff, notes C4–B4) [falsified if the growth is below half].
- **X-Br-C2:** the lip model dominates pp/mf and the upper register [falsified if a lip-parameter sweep inside player-data ranges moves the distance by less than the mf→ff change].
- **X-Br-C3:** linear bell beaming dominates above about 1.5 kHz [falsified if the position-ensemble spread above 1.5 kHz exceeds the mf→ff change].

**Milestones.**
- **M-BR-Q2:** M2.2-W3a–e, M2.3-W2, M2.1-W3b, M3.1-W2/W6/W7, M2.5-W1/W7, M1.2-W1a, M2.4-W1a/W2, M0.1-W3.
- **M-BR-Q4** covers trumpet and trombone only. A brass-wide claim requires the horn and tuba milestones as well.
- **M-HN-Q2:** adds M2.4-W1d and the hand fixture.
- **M-TU-Q2:** adds M2.5-W5a–c.

### 6.5 Tuned bars

**Description.**
- **Bars:** undercut xylophone and marimba bars, and vibraphone bars (aluminium), with **caliper-measured profiles**. Verify the species: the Iowa bars may be Honduras rosewood or synthetic, while the packs are Brazilian/Indian rosewood.
- **Supports:** nodal holes, plus cord and felt supports.
- **Resonators:** tube resonators and the bar–tube gap.
- **Air temperature.** Tube tuning drifts about 3 c/K against the bar, so air temperature is description data.
- **Mallets:** layered heads (yarn, cord, rubber; Iowa variants) with handle flexibility.
- **Vibraphone:** the damper bar and pedal; bowed-bar variants; motor discs only if the reference has the motor on.

**Chain.**

```
M1.3 bar (P1 h-ladder first via linear3, then TET10; cord/felt supports as M4.1 contacts)
  ↔ one coupled acoustic problem: M2.4 gap + tube interior + M2.5 exterior, two-way loading
    (radiation damping and frequency pull)
  ← M4.4-W3d mallet
```

The resonator's 3f lines up with the xylophone's tuned 1:3 partial (Bork 1995). Its detuning sets bloom and beating: the reference `amp_fluct` is 2.0–3.6 dB.

**First slice (Bar-S1).**
- **Model.** An existing `fs-solid::linear3` 3-D modal bar on an h-ladder, with the fs-plate tapered Hermite/Timoshenko beam as a cross-check rung only.
- **Radiation** by the existing fs-bem velocity-BC solve.
- **Excitation** by an existing Hertz/felt mallet.

Use the measured bar's actual profile, supports and mallet geometry/materials when available. A designed undercut or idealized mallet is a declared fixture, not a reconstruction of the recorded specimen; omitted resonator/support/player effects remain subject to §1.4.

**Experiments.**
- **X-Bar1:** a 3-D FE model of a caliper-measured tuned bar predicts that bar's measured partials within 5 c. Without a measured profile, the undercut is optimized to 1:4 or 1:3 at ET. The undercut is then a declared synthetic construction design (§1.3), and Bar1 tests only untargeted quantities: the third partial, torsional and lateral modes, contact time against Chaigne & Doutaut (1997), and decay against Bork (1995).
- **Staged validation (IndigoLeopard duel).** A same-specimen bar reduces description uncertainty but never removes it, and each stage claims only its own scope:
  1. **Bare bar component (X-Bar1).** Measured-force (impact-hammer) excitation of the caliper-measured bar: modes, mode shapes and decay. This stage makes no contact or resonator claim.
  2. **Tube passive response (X-Bar2a).** The resonator tube's input impedance and radiation, measured and computed.
  3. **Assembled resonator (X-Bar2).** The bar over its tube with tuned, detuned and removed-tube controls. The model must predict the quantitative changes in the mechanical pole and decay and in the radiated pressure. `amp_fluct` alone cannot distinguish reciprocal coupling from two independent ringers.
  4. **Complete instrument.** Mallet contact, supports and the player's stroke.
- **Solver route for the coupled stages (M2.5-W10).** The coupled stages need a qualified route for the open, finite-wall tube. Either closed-surface BEM with the finite wall or interior/axisymmetric FE is acceptable, whichever qualifies. fs-bem's plain-CBIE thin-face treatment is evidence for thin opposing faces, not a tube qualification.
- **First specimen.** A metal (aluminium) bar component is the practical first co-characterized specimen: low material uncertainty and cheap measurement. A complete instrument still needs mallet, contact, resonator and supports.

**Milestones.**
- **M-BAR-Q2:** Bar-S1, M1.3-W1a/W2a/W3, M2.5-W4a, M2.4-W1a/W3b, M4.1-W1, M4.4-W1/W3d, M8-W3i.
- **M-BAR-Q4.**
- **M-VB-Q2 (vibraphone):** aluminium bars, the damper bar and pedal as M4.1 contacts, and the motor-driven resonator discs as a time-varying aperture (M5.3-W2b) when the reference has the motor on.

### 6.6 Guitar

**Description.**
- **Instrument:** a classical guitar matching the UIowa reference (Raimundo 118): fan-braced top with thickness map; back with bars; ribs, linings, neck and heel; bridge, saddle and nut; **frets along the full length**; cavity with soundhole.
- **Strings:** nylon trebles, and multifilament-core basses with silver-plated copper windings (gauges, tensions).
- **Player contact:** thigh, chest and forearm against the body as compliant supports.
- **Steel-string X-braced guitar** once a licensed anechoic reference exists.

**Chain.**

```
M1.1 strings (both planes, longitudinal, wound basses; all six strings on one reciprocal bridge
     termination, so open strings ring sympathetically; fret and fingerboard contact along the length)
  ↔ M1.4 body (prestressed equilibrium; M1.2 shells with junctions, braces)
  ↔ M2.4 cavity air with soundhole
  ↔ M2.5 radiation with fluid loading
  → receiver
```

The finger/nail exciter is M4.4-W3b; stops are M5.3.

**First slice (Gt-S1).** `ReciprocalStringBank::from_material_plate_chart` (a flat plate, recorded as a model-form rung) reaching Q1/Q2. This tests the existing coupling; it does not establish that a flat plate is acoustically equivalent to a real guitar body. The final description uses the real body assembly unless §1.4 equivalence is demonstrated.

**Experiments.**
- **X-Gt1:** two-stage partial decay from both polarizations coupled through the bridge [falsified if single-plane matches; requires the §7.2 per-partial envelope].
- **X-Gt2:** intrinsic string loss sets the late decay [analysed to the reference floor, ≥ 5 s].
- **X-Gt3:** fingertip + nail release (Pavlidou 1997) sets the attack spectrum.
- **X-Gt4:** separate the low HNR of the A and E strings into four candidate causes: stiff-partial binning, sympathetic open strings, longitudinal modes, and fret contact.

**Milestones.**
- **M-GT-Q2:** Gt-S1, M1.1-W1c/W2/W5b/W6/W9, M0.4-W2, M1.2-W1a/W2a/W6, M1.4-W1/W2, M2.4-W3a, M2.5-W4b, M4.4-W3b, M4.1-W1.
- **M-GT-Q4.**

### 6.7 Violin family (pizzicato first, then arco)

**Description.** Violin first (the 2012 UIowa reference; the pre-2012 violin is the same-type ruler), then viola, cello and double bass:
- **Body:** arched top and back with thickness maps; ribs, linings and blocks; bass bar and soundpost; f-holes; neck, fingerboard and chin rest; **varnish and ground as laminate layers**.
- **Bridge:** its own modes, with feet in contact.
- **Strings:** tailpiece and afterlengths; G/D/A/E strings with cores and windings.
- **Bow:** stick, hair ribbon, rosin.
- **Holding:** chin, shoulder and hand; for the cello, endpin, knees and floor. Floor and endpin are in the observation ensemble.
- **Finger stops:** fingertip on string on fingerboard. This damps the nut-side segment, so open strings ring longer and brighter.

**The 2012 references are non-vibrato** (f0 s.d. 1.9 c for violin and 1.5 c for cello). Vibrato modelling, including the body-filter effect (Weinreich 1997), is therefore off the Q2–Q4 path. Check the pre-2012 ruler violin for vibrato before using it as a ruler.

**Chain.**

```
M1.1 strings (torsion, both planes, wound; all strings on the bridge termination)
  ↔ M4.4-W3a bow (ribbon line contact; M4.2-W3 thermal third body; M4.3 limit rung)
  ↔ bridge (M1.3 solid, foot contact)
  ↔ M1.4 body (prestressed equilibrium: tension, downbearing, soundpost compression;
     M1.2 arched shells + junctions + laminate)
  ↔ M2.4 cavity air through f-holes
  ↔ M2.5 radiation (+ player scatterer)
  → receiver
```

M4.4-W3b supplies the finger, and M6.2 the bowing gestures.

**Experiments.**
- **X-Pz1:** pizzicato (Iowa has pizzicato scales for every bowed instrument) through the full body with the finger fixture reaches Q2 **before any friction module**.
- **X-Bw0:** bow settings inside the Schelleng range (the old gate point sat at maximum bow force); slip fraction ≈ β; flattening.
- **X-Bw1:** body screen. A clean bridge force through the new body must make the nearest real instrument a violin; if it does not, body data or coupling is wrong first.
- **X-Bw2:** friction law vs contact kinematics vs player micro-variation [falsified if ribbon + thermal third body + measured bow speed/force fluctuation cannot reach reference HNR (violin 16–27 dB, cello 24–34 dB) and amp_fluct].
- **X-Bw3:** onset gestures from Guettler (2002) and Galluzzo & Woodhouse (2014) land in the reference's attack class (Iowa attack times 290–610 ms). Q1 for bowed strings includes pre-Helmholtz transient length.
- **Cello wolf** at the main body resonance.

**Milestones.**
- **M-VNP-Q2 (pizzicato):** M1.1-W1c/W2/W5b/W6/W9, M4.4-W3b, M1.3-W1a, M1.2-W1a/W2a/W3/W4/W6/W7, M1.4-W1/W2, M2.4-W3a, M2.5-W4b, M8-W3c.
- **M-VN-Q2 (arco):** adds M4.4-W2/W3a, M4.2-W3a/b, M4.3-W2, M5.3-W2b, M5.2-W3, M6.2-W1.
- **M-VN-Q4.**
- **Viola, cello and double bass:** each has its own geometry item (M8-W3c2–c4) and milestone rows (M-VA, M-VC, M-DB). Cello and bass additionally need M2.5-W5a–c.

### 6.8 Flue family

**Order.** Recorder and organ pipe come first, because their windway is fixed. They are **QP/Q1 validations of M3.3 against Auvray, Terrien and Verge**, since Iowa has no reference for them. The transverse flute follows.

**Flute description.** A Boehm flute: 19 mm bore, tapered head joint, cork cavity, embouchure hole and chimney, keyed holes. Embouchure control (jet length, offset, lip coverage) comes from player data (Fletcher 1975; de la Cuadra 2005).

**Chain.**

```
M3.3 zonal (M3.2 resolved compressible zone at mouth/labium ↔ M2.2 DG bore)
  → M2.5 radiation
```

M3.4 accelerations serve exploration only.

**Experiments.**
- **X-Fl1:** on the Boehm description, jet length against pressure follows Fletcher (1975) and de la Cuadra (2005). If the fundamental speaks, the register error was the control surface; the 0.30 m × 6 mm pipe is not a flute.
- **X-Rc1:** recorder thresholds and regimes against Auvray et al. (2012) and Terrien et al. (2014).

**Milestones.**
- **M-RC-QP1:** M3.2-W1/W2a/W2b/W3a/W3b/W4, M3.3-W1/W2a–c.
- **M-FL-Q2:** M-RC-QP1, M8-W3g, M2.5-W1.
- **M-FL-Q4.**

### 6.9 Voice (last; a licensed reference is required)

Voice is a composition milestone, not a module. The parts:
- **Folds:** an M3.1 valve pattern with **layered anisotropic active tissue** (M0.3-W3; Zhang 2016): mucosal wave, collision, muscular prestress.
- **Airways:** subglottal tract (trachea/bronchi; Lulich 2010).
- **Tract:** an M2.2 area function with yielding walls (`tract.rs`), moving articulators, piriform sinuses (Dang & Honda 1997) and a nasal branch (M2.1-W8).
- **Radiation:** mouth radiation (M2.5).
- **Noise:** aspiration and fricatives come from the resolved glottal and constriction flow (M3.1-W8, M3.3). M3.4-W3 is exploration only.
- **Variability:** jitter and shimmer from measured physiological tension tremor.
- **Coupling:** source–filter coupling (Titze 2008).

**Milestone.** M-VO-Q2, after licensed reference M8-W6b.

### 6.10 Duduk

**Description.** An Armenian duduk (the measured specimen; key of A first):
- **Body:** apricot wood (*Prunus armeniaca*) with a near-cylindrical bore, eight finger holes, a thumb hole and any tuning holes, all as measured.
- **Reed:** the large flattened cane double reed (*ghamish*/*ehegh*), with its cap and the tuning bridle that sets the opening. Its blade thickness map and internal volume are measured, and it is kept in its soaked material state (M8-W2a wet cane).
- **Player:** the lips cover a long part of the blades, and the airway couples through the large reed.

**What makes it special.** Maugeais & Dalmont (2024) measured a duduk with an impedance sensor and an artificial mouth. They identify the **low reed resonance frequency** as the main specificity of the instrument. Their dimensional description and measurements are usable data. Their calibrated model parameters are not admissible inputs (§1.3).

**Chain.**

```
player airway (M6.4-W2) ↔ M3.1 two-blade valve (reed continuum from the vamp map and wet-cane tensor;
                          distributed lip contact over the blades, M3.1-W6 + M6.4-W3;
                          bridle as a contact constraint; rung-B resolved channel flow)
                        ↔ reed interior as an M2.4 multiport
                        ↔ M2.2 bore + M2.1-W4a tone holes with fingers (M6.4-W4)
                        ↔ M2.5 multi-aperture radiation
                        → receiver
```

Pitch bending by lip force is a performance control (§1.3).

**Experiments.**
- **X-Dk1 (QP):** input impedance of the described specimen against its measured impedance [peaks > 10 c off on a block of fingerings].
- **X-Dk2:** the reed continuum (geometry + wet-cane tensor + lip contact) reproduces the low reed resonance measured on the same reed with an artificial mouth [falsified if it is > 15 % off].
- **X-Dk3:** playing frequency and threshold against blowing pressure and lip force follow the same-reed artificial-mouth measurements [falsified if they fall outside the measurement spread].

**Milestones.** M-DK-QP1, M-DK-Q2, M-DK-Q4, M-DK-Q5 (Appendix E). They depend on:
- the double-reed modules (M3.1-W1d, M3.1-W4);
- the lips fixture (M3.1-W6, M6.4-W3);
- geometry M8-W3l;
- references M8-W6c (Iowa has no duduk).

### 6.11 Irish tin whistle

**Description.** A D whistle.
- **First specimen: a Clare D**, the model Davies, Pinho & English (2006) measured. Their published flow, geometry and acoustic measurements serve as an independent comparator. Their inferred 21.5 mm end correction is an output to compare against, never a physical input.
- **Later specimens** (OD15): a Generation (brass tube, plastic head), a Clarke (conical rolled tinplate, wooden block), a Feadóg.

The description covers:
- **windway:** length, height, width and chamfers;
- **window and labium:** window length, labium edge and angle;
- **tube:** bore and wall;
- **six finger holes,** with the player's fingers, including half-holing;
- **player:** cuts, taps and rolls as finger gestures, the mouth seal on the mouthpiece, breath pressure and tonguing.

**Physics.** A duct flute of the recorder family:
- a resolved jet from the windway at the window (M3.2/M3.3, M3.3-W3 edge zone);
- the bore (M2.2) and tone holes with fingers above them (M2.1-W4a, M6.4-W4);
- the second register by overblowing, with hysteresis;
- the "chiff" onset from jet formation, which comes from the resolved flow, never authored.

Dequand et al. (2003) show that mouth geometry, jet thickness and labium angle change regime stability. Thresholds and regimes are therefore graded, not only steady spectra.

**Why it matters early.** The whistle is cheap, its geometry is fully measurable (calipers, windway casts or CT), and an artificial mouth can blow it. That makes it the leading candidate for the first **wind** co-characterized specimen (M8-W7, OD15). Its physics readiness, however, waits for the resolved jet (M-RC-QP1). So it is the first wind specimen for passive geometry, QP and the resolved-flue regime-map benchmark, not instantly qualified audio. Acquisition starts early, and the computation follows the recorder validation (IndigoLeopard: choose on acquisition cost *and* physics readiness).

**Leads.**
- Davies, Pinho & English (2006): measurements on a Clare D whistle, the same model as our first specimen but not the same object.
- Dequand et al. (2003).

**Experiments.**
- **X-Tw1 (QP):** same-specimen input impedance.
- **X-Tw2:** thresholds and the regime map (first and second register, hysteresis) against artificial-blowing measurements of the same specimen [falsified if thresholds fall outside the measurement spread].
- **X-Tw3:** geometry sensitivity. Changing windway height and labium angle on the described specimen shifts regime boundaries in the measured direction and order before any listening.

**Milestones.** M-TW-QP1, M-TW-Q2, M-TW-Q4, M-TW-Q5. They depend on:
- the recorder validation (M-RC-QP1);
- M2.1-W4a and M6.4-W4;
- geometry M8-W3m;
- references M8-W6c;
- the specimen pair M8-W8c.

### 6.12 Andean pan pipes (siku/zampoña; antara)

**Description.**
- **Siku:** two rows (*ira* and *arka*) of cane or bamboo tubes stopped by natural nodes, traditionally played in hocket by two players.
- **Antara:** a single row. Modern antaras are cane; the pre-Columbian Nasca antaras are ceramic.

The description covers per-tube length, inner diameter, wall thickness, node (stop) shape and rim cut, plus the material (cane or bamboo species identified on the specimen).

Museum records such as Met 1995.541.8 (Bolivian siku) and Horniman M32-1990 (siku-malta) give identity, materials and overall dimensions, but not internal lumen or node geometry. No accessible CT dataset of Andean antaras was found. Interior geometry of ceramic antaras is an acquisition task, and the plan does not promise it.

**Physics.**
- **The jet.** Edge-blown stopped pipes with no windway: the jet is formed by the player's lips across the open end. Lip aperture, lip-to-rim distance and partial lip coverage of the tube end come from M6.4-W3. The jet itself is resolved (M3.2/M3.3, rim edge zone M3.3-W3).
- **Even harmonics.** Fletcher (2005, JASA) shows that the jet offset is set by the jet's own fluid dynamics rather than by lip position. This produces substantial even harmonics in the upper notes of a stopped pipe, which is a falsifiable signature.
- **Noise and bends.** Breath noise comes from resolved turbulence. Partial coverage of the tube end bends the pitch.
- **Hocket.** Two performances alternate notes between the ira and arka rows.

**Experiments.**
- **X-Pp1:** even-harmonic amplitude in the upper tubes against Fletcher (2005) [falsified if the resolved jet gives near odd-only spectra where Fletcher reports near-equal even harmonics].
- **X-Pp2 (QP):** per-tube impedance and sounding pitch against the measured specimen.
- **X-Pp3:** pitch bend against lip coverage, measured on the specimen.

**Milestones.** M-PP-QP1, M-PP-Q2, M-PP-Q4, M-PP-Q5. They depend on:
- M3.3-W1 and M3.2-W4;
- M6.4-W3;
- geometry M8-W3n;
- references M8-W6c.

### 6.13 Candidate further instruments (owner picks, OD14)

Each candidate exercises existing generic modules, and none needs instrument-specific code. Feasibility depends on a reference recording and specimen data.
- **Shakuhachi:** edge-blown open bamboo with nodes (M3.3, M2.1). Shares physics with the pan pipes.
- **Didgeridoo:** a lip valve with strong vocal-tract coupling and an irregular, termite-hollowed bore that needs CT (M3.1-W6, M6.4-W2, M2.2).
- **Uilleann pipes and bagpipes:** double reeds, drones, and the bag as a compliant reservoir (M3.1, M2.1).
- **Sitar and tanpura:** the jawari bridge's distributed contact gives the buzz, plus sympathetic strings (M1.1-W9, M4.1). BR-03 must be resolved first.
- **Erhu:** a bowed string on a snakeskin membrane resonator (M1.1, M4.4-W3a, M1.2-W8, M2.4).
- **Handpan:** a shell with tuned dimples, struck by the hand (M1.2, M4.4-W3b), with geometry from scans.
- **Mbira and kalimba:** tines as M1.3 bars, with buzzers (M4.1).
- **Ocarina:** a vessel flute, i.e. a duct jet on a 3-D cavity (M3.3, M2.4).
- **Bansuri and ney:** open bamboo and cane flutes (M3.3, M2.1).

---

## 7. Acceptance and grading

### 7.1 Gates

| Gate | Meaning | Test |
|---|---|---|
| Q0 Builds and renders | radiated pressure at a virtual receiver via any committed CLI, natively built | native build + render receipt |
| QP Passive chain | linear, player-independent behaviour matches measurements of the **same specimen** where one exists (otherwise population data, flagged) | input impedance per fingering; bridge admittance and signature modes; modal ratios |
| Q1 Speaks | intended regime, right note within ±10 c after the allowed tuning action, threshold/onset within the literature range (first against machine-player data); thresholds/hysteresis mapped | time-domain sweeps on the render path (M5.6 cross-checks) |
| Q2 Alive (registered 2026-10-22, pexfp.357) | **sustained families:** the note's mean spectral flux after the attack, as a level, lies within the spread of its instrument's real recordings within 2 semitones (`instrgrade/identity.py: q2_alive`), on ≥ 80 % of held-out graded notes. In the same run: real recordings from a second library pass at ≥ 80 %, the five null baselines all fail, +20 c moves no margin by ≥ 1 dB. Identity (nearest set) and confusion are reported beside it, never gating; which instrument it is, is Q4's. **Struck and plucked families:** no Q2 claim until a decay-domain measure is registered (pexfp.360). Retired: identity margin ≥ 3 dB (it passed 1 of 46 real recordings from a second library and 28 of 56 frozen nulls) | `instrgrade/gate_run.py` |
| Q3 Inside the ruler | same-type real-vs-real ruler; descriptor changes pp→mf→ff and along scales match; each bound descriptor is inside the real spread; graded as a simulation ensemble marginalized over the observation prior | grader scorecard |
| Q4 Owner blind test | **discrimination-equivalence**, below | listening harness |
| Q5 Phrases | scales, slurs, articulation and dynamics changes pass Q3/Q4 | phrase grading |

**Bound descriptors for Q3:**
- attack-window ERB distance (0–150 ms);
- per-partial onset order;
- per-partial envelopes down to the noise floor, with two-slope fits and beat rate and depth;
- individually tracked partial frequencies with fitted B;
- cyclostationary noise coherence at k·f0;
- intonation fingerprint.

**Q4 discrimination-equivalence:**
- **Design.** A balanced two-alternative forced-choice test with a prespecified **equivalence band** on discrimination accuracy (default 0.40–0.60). Level, onset and pitch are matched, and the recording chain is emulated.
- **Pass rule.** **Both** confidence bounds must lie inside the band (TOST; Lakens 2017). The analysis uses an appropriate binary listener-and-item method (for example, mixed-effects logistic regression with listener and item effects) and a valid stopping rule.
  - Accuracy reliably below 0.40 is systematic reversed-label discrimination, and fails.
  - Non-significance is not equivalence. A confidence interval that crosses either bound is inconclusive, not a pass.
- **Spectrogram test.** A blind spectrogram-pair test is run alongside.

Every gate report includes raw calibrated spectrograms and absolute level, onset and pitch errors alongside the aligned descriptors.

**Physical-description qualification.** Q0–Q2 may report explicitly limited development fixtures. Final Q3–Q5 instrument acceptance additionally requires the real-geometry/material/player description of §1.4 and evidence for every simplification used. Q4's alignment is an experimental listening control, not permission to suppress physical errors in that evidence. Same-specimen validation, population realism and simplification equivalence must be reported separately.

**Binding acceptance (round-2 review).** Final acceptance of a family requires all of the following on the same frozen held-out renders:
- QP, Q1 and Q3;
- the listening-equivalence test;
- a **spectrogram-equivalence** test.

The spectrogram statistic, equivalence margin, confidence procedure, analysis windows and family/note coverage are prespecified before evaluation. Uncorrected pitch, onset time and calibrated level have their own binding tolerances. Alignment is diagnostic and cannot rescue a failed raw tolerance. Where the reference's calibration is unknown (the Iowa files are not level-calibrated), absolute level is reported as unverified, never inferred from normalized audio.

**Gate expansion rule.** These prerequisites block every applicable milestone and are expanded into Appendix E rows. An omitted edge never waives them. Interacting anatomy (M6.4-W2 airway, W3 lips/embouchure, W4 hands/fingers) is required wherever the family exercises it, unless explicit equivalence-based omission evidence exists (M5.5-W2). A speech fixture or a generic player model is not acquired interacting geometry.
- **Every Q2 milestone:** M7.1-W2, M5.5-W1 and the family's QP/Q1 milestone.
- **Every Q3/Q4 milestone:**
  - M2.6-W1/W2/W3, M6.2-W2, M7.1-W3 and M7.5-W1;
  - the §1.4 physical description of the specimen (geometry and the anatomy it uses);
  - for reed, lip, fold and jet families, the resolved flow owner (M3.1-W8 or M3.3-W1/W3) qualified for that aperture family.
- **Every Q5 milestone:** M7.1-W4 phrase fixtures, plus the configuration and contact mechanisms the phrases exercise (M5.3-W2b; M2.2-W3d for winds; the driven-body fixtures for strings and bars).
- **Unavailable references** leave the corresponding acceptance item blocked.

**Claimable-item readiness.** Before a work item can be claimed, it states:
- its required input fields;
- its source or specimen;
- units and uncertainty;
- the fixture or test target and the observable;
- the numeric acceptance band;
- the refusal condition.

A data item is done only when its named consumer consumes it successfully. A missing required field blocks the consumer and is never filled with a default.

### 7.2 Grader essentials and blind spots

**Essentials.**
- **References:** UIowa raw chromatic scales (the per-pitch edits are processed); `novib` until vibrato is modelled. The 2012 strings are non-vibrato.
- **Rulers and baselines:** same-type rulers (pre-2012 vs 2012 Iowa strings; TU Berlin anechoic after licence; RWC flagged non-anechoic); null baselines (sine, 1/n, odd 1/n, 1/n², and the reference with its life removed: `instrgrade/freeze.py: null_signals`) with margins.
- **Uncertainty, by kind:**
  - **Aleatory and population variables** (player micro-variation, observation position, specimen-to-specimen spread) use justified probability distributions. Ensembles are sized by a convergence check on the reported statistic, not by a fixed N. The observation prior is scored by its median.
  - **Epistemic unknowns** (unmeasured geometry or material properties) use admissible, correlated envelopes: bounds or scenario sets. Results are reported as envelopes, never as a median over invented draws. Eight draws are not coverage.
- **Other descriptors:** ERB-rate onset-aligned spectrogram distance split into attack, sustain and release; log-attack time; spectral flux; AM/FM modulation spectrum; jitter and shimmer; band HNR; fs-psycho descriptors.
- **Band split:** at the transition frequency (strings ≈ 1 kHz; winds at the tone-hole cutoff), with statistics above it.
- **Perceptual tolerances** where published (violin: Fritz et al. 2007).
- **Checks:** segmentation; pitch-shift invariance (+20 c ⇒ < 1 dB); baseline ordering; real recordings of the instrument from a second library graded beside every family (`instrgrade/gate_run.py`).

**Blind spots of the 2026-10-15 grader, fixed by M7.1-W2 before any Q2 claim:**
- **Windows** (`steady_window`: 0.30–1.10 s sustained, 0.02–0.52 s decaying) exclude the attack and the late decay.
- **Decay fit.** One straight-line decay fit over 0.05–1.5 s cannot see double decay or beating.
- **HNR bins.** HNR uses ±0.03·f0 bins at k·f0, so stiff partials (B·k³ ≳ 0.06) count as noise.
- **High-pass.** A 30–50 Hz high-pass removes tuba and bass fundamentals.
- **Noise type.** HNR is blind to pitch-synchronous vs stationary noise.
- **Analysis length.** Decaying notes must be analysed to the reference's own noise floor (≥ 5 s for guitar bass strings and low marimba bars).

---

## 8. Data program

**Geometry completion rule (M8-W3a–n including c2–c4, the anatomy items M6.4-W2–W4 and M8-W9).** Deliver a solver-usable physical description, not just a model file: source/specimen identity and rights; SI scale and coordinate frames; external/internal surfaces and thickness; part/material/contact/support mapping; calibrated reconstruction evidence; and uncertainty for unobserved geometry or properties. Preserve the detailed source and identify any simulation simplifications with their §1.4 evidence or an explicit pending limitation. A documented unknown permits an exploratory fixture, not a fidelity claim it affects. Apply this rule to Appendix E's DONE-WHEN entries as well.

**Acquisition routes.** Prefer source CAD, scans and measured sections where available; use Blender/photo reconstruction with independent scale and section constraints otherwise. Start one complete specimen description alongside the first instrument slice, rather than collecting cosmetic models of every family. Bind player anatomy/configuration and material state to the same experiment. This extends the existing ledger and descriptions; do not create a separate geometry-governance system.

| Item | Content | Sources to verify (each gets a licence row) |
|---|---|---|
| M8-W1 Ledger | one row per source: what, licence, uncertainty, digitization, decision | — |
| M8-W2a Woodwind materials | *Arundo donax* cane dry and water-saturated (Obataya & Norimoto 1999); grenadilla/ABS; ebony; cork; pad leather and felt | literature |
| M8-W2b String materials | nylon, gut, fluorocarbon; multifilament cores; silver/copper windings; horsehair; rosin (third-body data) | literature (Smith & Woodhouse 2000; maker data) |
| M8-W2c Wood loss | loss against frequency and direction; moisture and temperature dependence | Ono & Norimoto; Obataya; Brémaud |
| M8-W2d Tissue | lip, fingertip, tongue, vocal fold (layered) | biomechanics literature |
| M8-W2e Heads and felts | wool felt coupons, rubber, yarn, cord | literature, 7sga6.5.7 / .3.12 |
| M8-W2f Brass alloys | against temperature and temper | literature, 7sga6.1.8 |
| M8-W3a–n Geometry | a: R13 + B4 · b: Bach 180S37 · c: violin (c2 viola, c3 cello, c4 double bass) · d: classical guitar · e: soprano/alto sax · f: oboe/bassoon · g: Boehm flute · h: recorder/organ pipe · i: bars (calipers) and resonators · j: trombone/horn/tuba · k: tract MRI sets · l: duduk and reed · m: tin whistle · n: siku/antara | UNSW Music Acoustics, published bores and theses, Strad3D/Bissinger, guitar plans, maker data; owner measurement (OD4) |
| M8-W4 Players | blowing pressure (Fuks & Sundberg 1999; Fletcher & Tarnopolsky 1999 trumpet); lip force and tonguing (Pàmies-Vilà et al.); playing tract impedance (Fritz & Wolfe 2005; Chen et al. 2008/2009; Scavone et al. 2008; Fréour & Scavone 2013); lip mechanical response (Newton et al. 2008); horn hand posture (Dell et al. 2015); bowing (Askenfelt 1989; Schoonderwaldt & Demoucron 2009; Guettler 2002 onsets); guitar finger/nail plucks (Pavlidou 1997); flutist jet (Fletcher 1975; de la Cuadra 2005); mallet stroke velocity | literature |
| M8-W5 Machine players | artificial-mouth thresholds (Dalmont & Frappé; Atig et al.); bowing machine (Galluzzo & Woodhouse); recorder thresholds (Auvray; Terrien); lip open area (Bromage) | literature |
| M8-W6a References | `data/instruments/references.json` (URL, sha256, licence) plus a fetch script; audio uncommitted and off the Mac data volume; Iowa ambient silence files; Iowa equipment page verified (mic, chain, chamber) | UIowa |
| M8-W6b Extra references | TU Berlin anechoic array recordings (rulers and directivity, after licence); licensed recorder, voice and steel-string references | — |
| M8-W6c Unusual-instrument references | duduk, tin whistle and pan-pipe recordings: licensed anechoic recordings if any exist, otherwise owner recordings of the measured specimens in a quiet, characterized room. Room and chain are declared, and a non-anechoic room is flagged in grading. | — |
| M8-W8a–f Same-specimen passive pairs | geometry plus measured complex input impedance (amplitude, linewidth, phase) or admittance of the same physical object: a clarinet, b saxophone, c tin whistle, d duduk, e bar (modes and decay), f guitar/violin bridge admittance | owner measurement or service (OD15); published pairs where the specimen is identified |
| M8-W9 Owner anatomy | intraoral scan, structured-light face and lip scan, MRI if available, with consent (OD16) | — |
| M8-W7 Co-characterized specimens | per family where OD4 approves: one physical instrument with geometry, materials, passive response (impedance or admittance), controlled excitation, microphone arrangement and sound. **Measurement back-action** is part of the measured assembly: pressure taps, sensor cavities and attached transducers are described, or bounded with capped, probed and non-contact comparisons; calibration alone does not bound loading. State drift is logged over repeated takes (humidity, warming, reed water, supports). First specimens: a metal bar component, and a Clare D tin whistle for wind passive geometry and QP | owner budget |

Facts (dimensions, constants) are usable; verbatim tables and figures may not be. The ledger records the decision per source.

Concrete acquisition leads checked on 2026-10-15: [Strad3D CT study](https://strad3d.org/ct.html) for violin geometry; [Brandeis Digital Amati Messiah CT/STL catalogue](https://scholarworks.brandeis.edu/esploro/outputs/dataset/1716-Stradivari-violin-Messiah-Ashmolean-Museum/9924324268601921), explicitly permission-only; [USC volumetric vocal-tract MRI](https://sail.usc.edu/span/morphdb/), speech configurations with an access form; and [Fraunhofer MUSICES](https://www.iis.fraunhofer.de/en/ff/zfp/projects1/musices.html) for instrument CT acquisition guidance. These are leads, not downloaded, licensed, material-characterized or performance-matched assets. Their availability does not establish a matching specimen recording.

Leads for the added instruments and anatomy (IndigoLeopard and MagentaBadger, 2026-10-16). None of these has been downloaded or validated.
- **Duduk:** Maugeais & Dalmont, *What makes the duduk special*, Acta Acustica 8, 60 (2024), DOI 10.1051/aacus/2024037. It gives a dimensional description and artificial-mouth measurements; its calibrated model parameters are not admissible.
- **Tin whistle:**
  - Davies, Pinho & English, *Some acoustic characteristics of the tin whistle*, Proc. IOA 28(1) (2006). It studies a Clare D: the same model as the first specimen, not the same object.
  - Dequand et al., JASA 113 (2003), on mouth geometry, jet thickness and labium angle against regime stability.
- **Pan pipes:**
  - Fletcher, *Stopped-pipe wind instruments: acoustics of the panpipes*, JASA 117 (2005).
  - Met 1995.541.8 and Horniman M32-1990 object records.
  - Swift & Bosworth, Internet Archaeology 56 (2021): scan and photogrammetry methodology on Roman and Egyptian panpipes only.
  - No accessible Andean antara CT dataset was found.
- **Anatomy:**
  - USC 75-speaker speech MRI (Lim et al., Sci. Data 8, 187, 2021; CC BY 4.0).
  - Dresden Vocal Tract Dataset (vocaltractlab.de, DVTD): MRI, printable and FE geometry, and measured and simulated acoustics and aerodynamics for 22 speech sounds of two subjects. It is a passive geometry-to-acoustics validation seam, not playing posture or tissue data.
  - Göttingen real-time MRI of brass players (access to verify).
  - Playing-tract impedance (Fritz & Wolfe 2005; Chen, Smith & Wolfe 2008–2010).

---

## 9. Runtime, site and live

### 9.1 Offline renders

- **Command:** `music_render instrument <description> <performance> <out.wav>` (M5.4-W2). Until then, the existing CLIs are used.
- **Host:** the render host is ts1 (OD7); never ts2 at scale.
- **Ladders and logging:** gate renders carry ladders (M5.5), and every campaign logs command, wall time, host and checkpoint.

### 9.2 Site

- **Section.** One data-driven `<section id="instruments">`, a direct child of `<main>`, between `#piano` (`</section>` at index.html:219) and the PLAYER comment/`#player` (221–222). One nav link ("Orchestra") goes between the `#piano` and `#player` links (34–35), so `initNav` (main.ts:45–63) tracks it. Anchors were checked at f9c6f1fcb.
- **Module and data.**
  - **Mounting.** One lazy `src/instruments.ts`, mounted by one `initInstruments()` that copies `initVerisim` (main.ts:410–416). It renders its own `ClipGroupView` (src/clipgroup.ts) from `public/data/instruments.json`, using the Manifest/Group schema with one fixed gain per group.
  - **Data.** A new `tools/instrgrade/site_data.py` writes `instruments.json`, `audio/instr-*` and `img/instr-*` from scorecards. `build_data.py` is not edited, because it is a top-level script with hard-coded piano and percussion groups. Every URL goes through `util.asset()` (content-hash URLs, 13f89ae0a).
  - **Comparisons.** Real-vs-ours panels read instrgrade JSON. `verisim.ts`/`verisim.py` stay piano-only.
  - **Copy** is computed, never typed.
- **Visibility** per OD1.
- **Ownership.** `apps/instrument-lab/**` belongs to the t712z maintainer (session frankensim-de). IRP adds new files only, and messages the maintainer before touching any of these: `main.ts`, `index.html`, `verisim.*`, `vite.config.ts`, `.gitignore`, `vercel.json`, `cloudflare/_headers`, `cloudflare/wrangler.jsonc`, `build_data.py`, `clipgroup.ts`, `public/live/**`, `engine-*/**`.

### 9.3 Live

- **Engine.** One `engine-perform` crate compiles the generic renderer (M5.4) to wasm32+simd128, with baked linear operators shipped as data. It installs to `public/live/perform/` (the maintainer adds the gitignore line) and is served under `live/<hash>/`.
- **OD2 (changed default): no exception to the owner standard.**
  - **Fixed rung.** A live path runs the offline modules at one declared rung on every M5.5 axis (Δt/substeps, modes, h, taps), and passes the same gates at that rung. The rung is fixed before playback and identical on every device.
  - **No load adaptation.** There is no load governor (cf. `public/live/piano/worker.js:96–106`) and no device-chosen grid (cf. the kit worker).
  - **Substeps.** Substeps (PianoAudio 1/2/4×) are a Δt ladder axis, valid only with their ladder result.
  - **Slower than real time.** If the gated rung cannot run in real time, render-ahead finishes the performance, or a buffer sized from the measured deficit, before playback starts. It reuses fs-live-mix (`lane.js`, `mix-worklet.js`) and requires `skippedFrames = 0`. Otherwise the page uses offline-rendered playback.
  - **Existing engines.** The existing piano and kit live engines declare themselves reduced (`engine-piano/src/lib.rs:1, 62`). They are out of IRP scope and never a parity target.
- **Determinism.** Output is bit-identical to native on the same ISA, reusing fs-live-fma (engine-drum/fma-shim): `relaxed_madd` only behind a start-up fusing probe, correctly rounded software FMA otherwise, bit-identical either way. Cross-ISA runs agree within a declared tolerance until `det::` routing is complete. The music stack's determinism class is one-host today.
- **Budgets.** At most 4 GiB of address space, with a per-device operator byte budget (4096 modes × 4096 taps × 8 B = 134 MB per microphone). Threads need COOP/COEP. Neither deploy sets them today; adding them is a maintainer-coordinated L-2 step.
- **Piano migration.** The piano's engine migrates onto fs-string (M1.1-W1d), then onto engine-perform. It keeps byte-identity on the fp fixtures and on engine-piano's tests (clone, kernel_parts, audio_rate at 1×, scale_preparation) and `examples/validate.rs`, once those first pass natively.

---

## 10. Execution

### 10.1 Tracks (parallel, by tier)

| Track | Modules | First items |
|---|---|---|
| T0 Foundations | F1, M7.1, defects | F1, M7.1-W1, M2.2-W1, M5.1-W2, M2.6-W0, M5.5-W1 |
| TA Acoustics | M2.1–M2.6 | M2.1-W1/W5, M2.2-W2, M2.5-W2, M2.4-W0 |
| TS Structures | M1.1–M1.5 | M1.1-W1a (after OD10), M1.3-W1a, M1.2-W1a, M1.5-W1 |
| TC Contact & friction | M4.1–M4.5 | M4.1-W1, M4.2-W1 |
| TF Flow | M3.1–M3.4 | M3.1-W1a, M3.1-W3, M3.2-W1 |
| TR Runtime | M5.x, M6.x | M5.1-W0, M5.2-W1, M6.1-W1a |
| TD Data | M0.x, M8 | M8-W1, M8-W3a, M8-W6a, M0.2-W1 |
| TI Instruments | milestones, X-* experiments, first slices | M-CL0, Bar-S1, Gt-S1 |

### 10.2 Week one

1. **F1.** Native verification on ts1 with CARGO_BUILD_JOBS ≤ 12 and named targets: fix DEF-10; build and run `music_render`, `flue_contraction`, and the reed/wind/flue/bowed test targets. Every failure is filed as a bead.
2. **M7.1-W1.** Commit the grader and reproduce the audit's clarinet-bore row. Start M7.1-W2.
3. **Quick fixes with regression tests:**
   - M2.2-W1 (adjudicate DEF-01);
   - M2.2-W2 (first slice);
   - M5.1-W2 (DEF-05);
   - M2.6-W0 (DEF-12, coordinated with the remote lane);
   - M2.5-W2.
4. **Data:** M8-W3a R13 bore, lattice and fingering files, with the §8 geometry completion rule and a named route for missing internal geometry/player anatomy; M8-W6a references manifest; ledger rows. Start one complete physical specimen description in parallel with the early renders.
5. **Baselines.** Start M-CL0. Start Bar-S1, the cheapest complete structure + radiation + contact slice.
6. **Integrator audit.** Start M5.2-W1: list every render path's integrator and measure its warping.
7. **Strings.** After OD10 and the piano beads releasing `examples/grand_piano/**`, start M1.1-W1a.
8. **Geometry.**
   - Start M6.3-W1/W2 (ledger and import checks, Blender headless scripts) on the R13 (M8-W3a) and on one bar specimen.
   - Request OD15 specimens and services: a metal (aluminium) bar and a wooden marimba/xylophone bar, a Clare D tin whistle, a duduk with reeds; CT and impedance measurement.

### 10.3 Dependencies

Item-level edges in **Appendix E are authoritative**. The module sketch below is illustrative only.

```
M0.2 ─► M1.1/M1.2/M1.3 ─► M1.4 (+M1.5) ─► bodies ─► M2.4-W3/M2.5-W4
M0.3 ─► M1.3-W2b ─► M3.1-W6 lips / M4.4-W3b finger / voice folds
M0.4 ─► M4.2-W3 ─► M4.4-W3a bow
M2.1 ─► M2.2 DG (W3a–e) ─► M2.3-W2, M2.2-W3d (moving geometry), M3.3-W1
M2.4-W1 ─► M2.4-W2 ─► M2.1-W4a, M3.1-W4 ;  M2.5-W3a/b ─► M2.4-W3b, M2.1-W6
M5.1-W0 ─► every module's ports ;  M5.1-W1a/b ─► M5.4, M5.6, M5.3-W2b
M3.2 ─► M3.3 ─► M3.1-W8, M3.4, flue/recorder
M7.1-W1 ─► M7.1-W2 ─► every Q2+ milestone
```

Known bottlenecks, staffed first:
- **M2.2-W3a–e** (DG duct) gates M2.3-W2, M2.2-W3d, M3.3-W1 and all wind Q3/Q5.
- **M1.1-W1a–c** (fs-string lift) gates guitar and violin.
- **M7.1-W2** (grader) gates every Q2 claim.

### 10.4 Expected order of Q2, per track (no calendar claimed)

- **Winds:** clarinet → saxophone → oboe/bassoon → trombone + trumpet → horn → tuba → flute. Recorder/organ pipe reach QP/Q1 only, until a licensed reference exists.
- **Structures:** xylophone/marimba → guitar → violin-family pizzicato → violin-family arco → vibraphone.
- **Voice** comes last.

Each milestone records measured progress in its scorecard.

### 10.5 Build and verification discipline

- **Builds.** Use DSR for repo gates and `rch exec` for narrow probes. The git-archive export recipe (Appendix C) is for when RCH is contended; record which was used.
- **Jobs.** Use named targets and CARGO_BUILD_JOBS ≤ 12. Never build at scale on ts2 (it hosts production Postgres).
- **Mac.** Run `sbh check --need 20G` before heavy lanes. Campaign outputs, BEM bakes, LBM fields and UIowa audio never live on the Mac data volume (`~/Downloads` is on it).
- **Formatting.** rustfmt only on your own files; no workspace-wide fmt.
- **Commits.** Commit by pathspec.
- **Native verification.** Every item runs its tests natively before closing. Gate renders carry ladders.
- **CI is not evidence.**
  - `instrument-materials.yml` has no green run. Bootstrap fails on the franken_numpy `make.bat` line endings, and cargo fails `fsci-special =0.4.0` against `constellation.lock`'s frankenscipy 0.2.0 (06108aea).
  - No workflow builds `apps/instrument-lab`.
  - Only native DSR/rch receipts count.

### 10.6 Swarm rules (bind every [IRP] bead)

1. **Remote lane first.** Before launch, the owner directs the remote no-toolchain lane to stop touching IRP-reserved paths, or to adopt OD3 (OD13). Until then, IRP beads touching those paths are not claimed.
2. **Epic IDs.** The new epic gets a new ID; never reuse `frankensim-irp-epic-edinr.*`, which is tombstoned. The converter flushes and commits `.beads/issues.jsonl` by pathspec within minutes, and never runs `br sync --merge`. This plan is committed by pathspec before the first bead cites it.
3. **Claiming.** Claim = check `am agents list` task descriptions + take an exclusive Agent Mail reservation on the exact paths + re-verify the item's "Now" anchors at HEAD. Every "Now" anchor is a diagnosis valid at 47a361a87/9267fd56d; the bead re-derives the remedy at HEAD.
4. **fs-couple sub-lanes, one owner each:**
   - valve: `bernoulli_aperture/**`, `reed_bore.rs`, the reed path in `render.rs`;
   - duct-time: `driving_point*`, `wind_line.rs`, `mm_line.rs`, characteristic/impedance lines;
   - string/body: `bowed_string/**`, `string_specimen.rs`, `thin_plate*`, `plate_viscosity.rs`;
   - piano: `examples/grand_piano/**`, shared with 5c0sj/stxqd/t0odx/i5yei;
   - CLI: `src/bin/music_render/**`;
   - cavity/thermal: `cavity_modal*`, `cavity_phs/**`, `fs-thermal/**`.

   `lib.rs` mod lines are append-only.
5. **Site files.** `apps/instrument-lab/**`: new files only; message the t712z maintainer before touching shared files.
6. **Derived artifacts.** One registrar per wave commits the derived artifacts (Cargo.lock hunks, doc-facts, the source manifest last, consolidation-review.json, instrument-claims.json, music-baselines.json, schema-policy.json).
7. **Bead contents.** Each bead carries:
   - **DONE-WHEN:** a named test `path::fn`, its Gauntlet tier, and a numeric tolerance;
   - **native run:** host, command and exit code in the close reason;
   - **audible items:** before/after grader JSON;
   - **e2e:** a fixture added to `scripts/e2e/music_lane.sh` (existing JSONL-per-stage receipts with detailed logging).

   No new harness frameworks (RULE 0.1).

---

## 11. Risks and falsifiers

- **Data.** Data may be unavailable or unlicensed for a family. The response is to measure a physical instrument (OD4, M8-W7), never to invent data.
- **Resolved physics is expensive.** Resolved compressible flow (about 10⁷ steps per second of sound) and thermal third-body friction cost a lot. Offline cost is acceptable (§1.1). Accelerations exist only for exploration (OD11).
- **Full-band fields.** 3-D interior and exterior fields at full band need fast BEM and fine meshes. Ladders bound the error where the band cannot yet be reached, and claims beyond the band are refused.
- **Integrator warping is global** (DEF-04). Any render path not yet migrated to §3.3 carries a measured warping disclosure.
- **The coarse metric can pass sounds the ear rejects.** Q4 is the acceptance. Temporal and noise descriptors, same-type rulers, ensembles and held-out notes reduce the risk.
- **Piano regressions.** Lifting the piano onto fs-string risks regressions, so byte-identical output is required before switching.
- **Unverified pushes.** Code-first commits broke HEAD on the audit day, and a remote lane keeps pushing unverified code (OD3, OD13).
- **Non-unique friction roots.** Non-monotone friction gives a non-unique per-step problem, so the root-selection rule is a disclosed model assumption.
- **Geometry and anatomy acquisition.** Hidden geometry (bracing, internal bores, reed interiors, ceramic antara interiors) needs CT or measured sections. MRI of playing postures needs a facility and consent. Unknowns remain declared bounds and block the fidelity claims they affect (§1.4).
- **Tail hypothesis (cheap test).** Most of the clarinet residual after radiated output is the tone-hole lattice (X-Cl-R2).

---

## 12. Owner decisions (defaults apply if unanswered, except where marked "ask")

| # | Decision | Default |
|---|---|---|
| OD1 | When a family appears on the site: after Q4, or at Q2 labelled "in progress" | after Q4 |
| OD2 | Live engines: any reduced live discretization allowed? | **no**: same gates, or render-ahead/offline playback |
| OD3 | Forbid pushing instrument-crate code that was not natively built and tested | yes |
| OD4 | Approve literature data per family; budget to measure physical instruments (co-characterized specimens) | literature first; ask for specimen budget per family |
| OD5 | Player gestures only from published measurements and owner-authored performances; per-note controls by declared non-spectral objectives | yes |
| OD6 | IRP items may close the instrument clauses of MR beads without their unrelated MR prerequisites | yes |
| OD7 | Render host for multi-hour jobs | ts1 (not ts2) |
| OD8 | About 30-minute blind listening sessions per family; adjudicate the seven old receipts (q61wp.24) in the first | first at clarinet Q2 |
| OD9 | Realism to recordings before extreme-regime material work for instruments (breath gas state stays in) | yes |
| OD10 | Approve new crate fs-string, the M2.4 home (fs-feec extension, superseding 0ja4's interior scope), and **express permission to move/delete the files enumerated in §5.1**, per item | **Decided by the owner 2026-10-21:** (i) fs-string approved; (ii) M2.4 home = fs-feec extension; (iii) moves of the §5.1 files with `git mv` approved, and removal of an old entry point approved only after bitwise parity is proven |
| OD11 | May closures (accelerations) ever be the shipped render, or only accelerations checked against resolved reference renders? | accelerations only |
| OD12 | Supersede music plan v8 (COMPREHENSIVE_PLAN_FOR_OPTIMAL_MUSIC_ORIENTED_BUILDING_BLOCKS.md), DECISION_MUSIC_MEAN_FLOW.md and BAFFLED_INSTRUMENT_ENVIRONMENT.md where they conflict | yes |
| OD13 | **Urgent.** Direct the remote no-toolchain lane to stop touching IRP paths or to adopt OD3. Its commits run from 6113eca8c through b1efe209e and later, and now include `engine-piano/**`, `public/live/piano/**` and `.github/workflows/instrument-materials.yml`. Meanwhile IRP beads verify, ladder and consolidate the lane's existing owners; they never reimplement them | **ask** |
| OD14 | Which further instruments from §6.13 to add (shakuhachi, didgeridoo, uilleann pipes, sitar, erhu, handpan, mbira, ocarina, bansuri/ney) | **ask** |
| OD15 | Budget for specimens and measurement services: metal and wooden bars, Clare D tin whistle, duduk with reeds, siku; CT scanning; impedance and admittance measurement; an artificial mouth; a quiet characterized room for specimen recordings | **ask** |
| OD16 | May the owner's own anatomy be scanned for the player description (intraoral scan, face and lip scan, MRI if available)? | **ask** |
| OD17 | May numbers read off published figures be retained under facts-only redistribution? | **Decided by the owner 2026-10-21: yes**, when the ledger row states the reading error (stated or bounded uncertainty); unstated error stays refused |

---

## 13. Relationship to existing programs and beads

| Bead / program | Relationship |
|---|---|
| frankensim-irp-epic-edinr + .1–.30 | v2 IRP, **tombstoned**; never reuse the IDs; .10 matdb inventory and .7/.8/.9/.18 text mined for v4 |
| Music plan v8 and its registry | superseded where they conflict (OD12); no new registry rows |
| frankensim-material-reality-ui2if (MR program) | MR42 observers, MR68 reed from specimen (ui2if.5.9; audit its orphaned in-progress children against HEAD first), MR76 bowed friction (5.10; audit .5.10.2), MR31/32/33/35/56 (5.1/5.2/5.3/5.5/5.8), MR03 loss models (1.3, in progress), MR16 wood bundles (2.8): delivered by M items, linked as related, never as blockers (OD6) |
| frankensim-ext-contact-nonsmooth-lane-oh0i | M4.3-W1 delivers through it (fs-mbd stepper, fs-contact cones) |
| frankensim-ext-tribo-dry-baseline-tgbj | M4.2 supersedes its friction scope; link M4.5 |
| frankensim-b8bxd.5.2, .7, .7.1–.7.3 | M4.1-W1 builds on these |
| frankensim-b8bxd.9 | related to M1.2/M1.3 dynamics and fs-plate-transient wiring |
| frankensim-ext-acoustics-core-0ja4 | interior/coupling scope superseded by M2.4/M2.5 (OD10); NVH QoI kept |
| frankensim-rc-root-q61wp.37 | parent of §5.1; its no-blanket-migration and no-deletion clauses bind |
| frankensim-rc-root-q61wp.24 / .30 | OD8; .30 superseded by §9.3 engine-perform once that lands |
| frankensim-5c0sj, stxqd, t0odx, i5yei | in-progress piano beads; block M1.1-W1c until they release grand_piano/** |
| frankensim-mkh9r | Cargo cycle; F1 verifies it is resolved |
| frankensim-8w7ya | strict-clippy campaign; coordinate fs-material edits |
| frankensim-t712z | instrument-lab site; §9.2/§9.3; M7.1-W1 location. The maintainer owns asset hashing (13f89ae0a) and `live/<hash>` versioning |
| frankensim-0nxhr | drum kit; BR-02 recorded there |
| frankensim-5c0sj.2, frankensim-h7xu5.7 | one receiver API (M2.6-W1) and one renderer (M5.4) |
| frankensim-g42o | shares M1.2-W1a/b shell benchmarks |
| frankensim-x65w, frankensim-epic-foundations-huq.12 | related to M3.2-W2 |
| frankensim-7sga6.5.7 / .1.8 / .3.12 | feed M8-W2e/f |
| frankensim-b2can | bore extraction; related to M6.1-W2 |
| frankensim-3yr7n | IndigoLeopard's bead for this revision (real geometry, anatomy, cross-model research) |
| research workspace `/Users/jemanuel/brennerbot_sessions/frankensim-instruments-20261015` | IndigoLeopard's completed two-family duel (Opus 5.5, GPT-6.1-Sol): question of record, source leads, geometry patch, DUELING_WIZARDS_REPORT.md, ARTIFACT.md |

---

## Appendix A — Measured baseline (2026-10-15 grader; see §7.2 blind spots)

| Render | f0 (Hz) | Dist (dB) | dynamics ruler | neighbour ruler | nearest real | HNR ours/real | fluct ours/real (dB) |
|---|---|---|---|---|---|---|---|
| piano note A3 mf | 221.2 | 11.7 | 13.5, 11.1 | — | piano 11.7 | 17.1 / 23.6 | 0.78 / 0.78 |
| reed-sustained (committed) | 141.4 | 23.5 | 10.9, 23.4 | 15.1, 11.7 | double bass 15.5 | −5.6 / 30.4 | 3.18 / 0.22 |
| music_render reed (unclipped) | 143.0 | 23.1 | 10.9, 23.4 | 15.1, 11.7 | cello 17.3 | −3.7 / 30.4 | 3.36 / 0.22 |
| brass low | 212.9 | 28.6 | 3.9, 17.2 | 6.2, 4.3 | alto sax 10.4 | 41.1 / 22.1 | 0.0 / 0.09 |
| brass high | 441.8 | 32.3 | 3.5, 21.5 | 6.1, 3.7 | horn 11.1 | 34.0 / 30.6 | 0.0 / 0.13 |
| brass valve | 192.6 | 24.6 | 8.4, 12.4 | 6.8, 8.1 | alto sax 8.1 | 33.7 / 23.5 | 0.0 / 0.11 |
| bowed, old gate point, vs violin | 236.2 | 20.7 | 8.3, 11.0 | 11.4, 13.1 | cello 12.9 | −9.2 / 20.6 | 2.34 / 0.83 |
| bowed light vs violin pp | 231.4 | 15.8 | 15.3, 11.0 | 9.8, 8.0 | cello 11.4 | −3.7 / 21.2 | 1.79 / 0.76 |
| flue 100 Pa vs flute note C6 | 996.8 | 26.5 | 10.1, 10.0 | 16.3, 10.4 | flute 23.9 | 8.5 / 27.6 | 0.05 / 1.45 |
| flue 300 Pa vs flute note G6 | 1533.2 | 33.9 | 20.1, 23.6 | 17.4, 10.6 | flute 22.5 | 12.1 / 21.8 | 0.14 / 0.70 |
| guitar strum (chord) | — | 27.5 | — | 8.6, 9.7 | — | — | 2.09 / 4.77 |
| 3-mode string vs guitar note B3 | 243.0 | 23.5 | — | 17.3, 10.1 | cello pizz 19.3 | 42.1 / 22.0 | 1.13 / 1.00 |
| xylophone strike 1 | 544.8 | 28.4 | 13.0, 16.8 | 7.4, 9.0 | xylophone 21.7 | — | 1.73 / 3.15 |

Reed on 7.3 mm per-note tubes (bore pressure):

| Note | D3 | E3 | F#3 | G3 | A3 | B3 | C#4 | D4 |
|---|---|---|---|---|---|---|---|---|
| Distance (dB) | 26.3 | 30.5 | 28.2 | 22.2 | 21.4 | 18.8 | 20.8 | 18.4 |

All eight notes are 33–41 c flat (lossless length formula), with HNR 33–61 dB.

## Appendix B — Defect and breach register (anchors checked at 47a361a87 / 9267fd56d)

| ID | Where | What |
|---|---|---|
| DEF-01 | `fs-couple/src/driving_point.rs:211` vs `:153`; comment `:205–208` | `impedance_line` places Z without the conjugation `characteristic_line` uses, under the e^{-iωt} convention (`fs-duct/src/lib.rs:11`). The code comment claims the opposite was validated. Adjudicate by causality tests. |
| DEF-02 | `driving_point.rs:134–137`, `:70–74`; `mm_line.rs:294–297` | FFT defaults to 4 round trips (256..4096); 92 c wrap-around on a horn; F2 19 % low on a tract |
| DEF-03 | `fs-vfit/src/discretize.rs:758–774`; `fs-vfit/tests/characteristic_line.rs:75` | global 1/peak scaling, checked on a caller grid (3 frequencies in the test) |
| DEF-04 | fs-phs Gonzalez midpoint on quadratic storage; `acoustic_realize.rs:845–905` (KC, no ω·dt guard); `render/plate/impact.rs:4–10, 276` (impact-host strings, only an ω·dt < 0.9π guard); `bernoulli_aperture/modal.rs:222–225`, `paired.rs`, `uniform.rs` (valves; ω·dt ≤ 1 admitted at `modal/material.rs:121–122`); `CoupledModalSystem` trapezoidal connections (`coupled.rs:5–14`) | ω_num = (2/h)·atan(ωh/2): −58.6 c at 5 kHz/48 kHz; up to −131 c at the valve admission bound |
| DEF-05 | `fs-phs/src/lib.rs:878, 1154`; `fs-phs/src/prepared.rs:376`; `bernoulli_aperture/modal.rs:293–300` | ≤ 1e-6 relative residual accepted with relaxed/inconsistent scaling (disclosed in `StepRecord::solver_residual`); the modal valve's work/flow audit allows 128·`relative_tolerance` with `relative_tolerance` anywhere in (0,1), so it cannot fail at ≥ 1/128 |
| DEF-06 | `fs-tribo/src/lib.rs:620–640` | regularized ramp reaches μs at vs while the kinetic value is lower |
| DEF-07 | `fs-couple/src/unilateral_contact.rs:290–291` | hand-typed 1e-3 m/s ramp; `.unwrap_or(0.0)` |
| DEF-08 | `fs-solid/src/linear3.rs` | P1 tets in bending-dominated parts |
| DEF-09 | `fs-phs/src/lib.rs:1583` (label), `:2161` (formula), `:1978–1991` (radiating variant); consumers `cavity_modal.rs:73–79`, `cavity_modal/shared.rs`, `examples/textile_cavity.rs` | both neck ends use 8/(3π)a (labelled unflanged) whatever `MouthFlange` is; the radiating variant takes the flange resistance but not its end correction; the neck loss is an authored resistance; `RadiationAtResonance` is a resistance frozen at the isolated resonance |
| DEF-10 | `fs-couple/examples/material_strings.rs`; `render.rs:509, 544` | `?` on `RenderError`, which has no `std::error::Error` impl |
| DEF-12 | `fs-couple/src/render.rs:367–387` (`:385`) | reed audio = bore pressure + plate radiation + ρ·dF/dt/(4πrc) |
| DEF-13 | `fs-phs/src/lib.rs:4675` (`bernoulli_volume_flow`) | pure √|Δp|, returns 0 below 1e-12 Pa, no viscous term, singular Jacobian at closure |
| DEF-14 | `driving_point.rs:173` | Z(0) forced to 0: erases closed-volume compliance (ṗ = ρc²U/V) |
| DEF-15 | `instrgrade.py` (`steady_window`, `decay_profile`, HNR bins, high-pass) | grader blind spots, §7.2 |
| BR-01 | `examples/grand_piano/longitudinal_string.rs:18–20`; bead na4ju | 0.06 hand-estimated "~1.8° downbearing + 28 mm lever arm"; DONE-WHEN targeted Iowa phantom levels |
| BR-02 | `crates/fs-vkplate` | head tensions tuned to recordings; s0+s1β² loss from two T60 points |
| BR-03 | `data/matdb/seed-v1/contact-jawari-bone-bridge` | calibrated to a fixture's output |
| BR-04 | render paths | reed lay K = 1e7·w (`render.rs:227`); 5 Pa priming push (`render.rs:238`); face length 0.025 m; ζ 0.35; brass lip parameters (`tests/bakeoff_brass.rs:75`); glottis constants (`glottis.rs:267–290`); felt parameters (`grand_piano/felt.rs:11`); piano drag 0.60 and bending spectrum (`linear.rs:307, 336`); bow capture thresholds (`bowed_string/runtime.rs:184, 361`); `reed-baffled-*.performance` ζ = 0.35 and 5 Pa prime; `two_sided_reed.rs` (E 4e9, ζ 0.02, reflection 0.7, penalty 1e8); FlueVoice environment `displacement_gain`/`seed_amplitude_rel` (acceleration only, OD11); `textile_cavity` neck resistance 2000 Pa·s/m³ (`air.txt`); `cavity_phs` effective conduction distances (`thermal/layers.rs:46–49`) |

## Appendix C — Audit materials and recipes

- **Listening pack and staging page:** `~/Downloads/frankensim-instrument-reality-check-2026-10-15/`.
  - `index.html` (staging page);
  - `audio/`, `wav/`, `figures/`, `results/`;
  - `tools/`: the grader, committed by M7.1-W1 without absolute paths: instrgrade.py, grade_run.py, summarize.py, radiate.py (DIAGNOSTIC-ONLY) and uiowa-links.json. `make_site.py` and `build_pack.py` stay in the pack;
  - `research/`: the audited modes-of-reasoning analysis.
- **Remote build recipe:**
  1. `git archive -o` HEAD and the siblings at local HEAD.
  2. scp the archives to ts1 or trj.
  3. Build with `RCH_SHIM_LOCAL_IDE=1 CARGO_BUILD_JOBS=12` and named targets.

  Two cautions: dcg blocks `git archive | ssh tar` and redirects to `$VAR` paths, so use literal paths; and never run `pkill -f` inside ssh (it kills the session).
- **Clarinet bore:** `music_render wind <perf> out.wav` with `termination flanged`, `cylinder 0.0073 0.57`, pressure `2.8e3`. Flue: `examples/flue_contraction.rs`.

## Appendix D — Validation benchmark index (per module)

| Module | Analytic / manufactured | Literature / measured |
|---|---|---|
| M0.2 visco tensors | Kramers–Kronig; passivity under random multiaxial histories | wood loss (Ono & Norimoto; Obataya; Brémaud) |
| M1.1 strings | stiff-string inharmonicity; KC Duffing glide; torsional frequencies; MMS for the rod equations | Raman waveforms; phantom partials (Conklin 1999; Bank & Sujbert 2005); Chabassier et al. 2013; wound strings (Järveläinen; maker data) |
| M1.2 shells | MacNeal–Harder; Leissa 1973; Chapelle–Bathe locking | violin plates (Hutchins; Jansson); guitar tops (Elejabarrieta); bell partials |
| M1.3 solids | Cook's membrane; patch tests; beam frequencies; Treloar | undercut bars (Bork 1995; Petrolito & Legge); bridges (Reinicke) |
| M1.4 assembly | assembled vs monolithic | violin signature modes (Bissinger; Strad3D); guitar body (Carcagno) |
| M1.5 eigen | known spectra; residual certificates | — |
| M2.1 ducts | cylinder/cone closed forms; Levine–Schwinger; Pagneux multimodal cases | same-specimen impedance (UNSW instrument or owner specimen); Ernoult; pads (Dalmont; Lefebvre); tone holes (Lefebvre & Scavone 2012) |
| M2.2 DG | impulse responses; impedance vs M2.1 < 1 c; moving-geometry energy | — |
| M2.3 nonlinear | Burgers; shock distance; entropy stability | trombone ff (Hirschberg 1996); brassiness (Myers 2012) |
| M2.4 interior | box modes; Helmholtz resonator; damped cylinder; narrow slit | tone-hole parameters (Lefebvre FEM); violin A0/A1 |
| M2.5 radiation | spheres; baffled piston; multipoles; FMM tolerance | directivity (TU Berlin; Meyer; Pätynen & Lokki 2010) |
| M3.1 valves | manufactured flows | Dalmont 2003; Almeida 2007; Bromage 2010; Newton 2008; Wilson & Beavers 1974 |
| M3.2/M3.3 flow | jet self-similarity; edge-tone stages; interface reflection < −40 dB | Auvray 2012; Terrien 2014; Giordano 2014; Yokoyama 2015; Verge |
| M4.2 tangential | Cattaneo–Mindlin; benchmark laws | Smith & Woodhouse 2000; Galluzzo, Woodhouse & Mansour 2017 |
| M4.3 nonsmooth | block on incline; Schelleng limits | Galluzzo & Woodhouse 2014; Desvages & Bilbao 2016 |
| M4.4 driven | Hertz contact times | Pitteroff & Woodhouse; Guettler 2002; Pavlidou 1997; Chaigne & Doutaut 1997 |
| M6.3 geometry | reconstruction of a known object within tolerance | caliper cross-checks on every specimen |
| M6.4 anatomy | — | playing-tract impedance (Fritz & Wolfe 2005; Chen et al.); lip open area (Bromage 2010) |
| Duduk | — | Maugeais & Dalmont 2024 (measurements only) |
| Tin whistle | — | Davies, Pinho & English 2006; Dequand et al. 2003 |
| Pan pipes | — | Fletcher 2005 |

## Appendix E — Work-item register (authoritative)

**Beads refinements (v5.1).** The 346 work items are tracked under
`frankensim-irp-v5-instrument-realism-pexfp`, alongside seven owner-decision
records. Appendix E supplies the core scope and prerequisites. The task records
add stricter acceptance criteria, named tests, numerical bands, falsifier
protocols and additional blocking prerequisites; these additions bind without
weakening the core requirements or the owner standard. The published roadmap
includes these refinements and the union of the blocking prerequisites at build
time, with a dated task snapshot and a
[download of the task refinements](https://frankensim.org/music-plan/bead-refinements.md).
Rebuild and redeploy after changing Beads. Task state is execution bookkeeping,
not permission, a passed quality gate or proof of validated instrument sound.

How to read this register:
- **Size:** S/M/L (§0).
- **Deps:** blocking items. "OD" entries are owner decisions, and "bead:" entries are external beads.
- **DONE-WHEN:** the core acceptance. Every item also satisfies §10.6 rule 7 (named test, tier, tolerance, native run, e2e fixture for audible items) and §3.4.

### E.1 Foundations and validation

| ID | Title | Size | Deps | DONE-WHEN |
|---|---|---|---|---|
| F1 | Native verification of instrument targets | M | bead:mkh9r | DEF-10 fixed; music_render, flue_contraction, reed/wind/flue/bowed targets built and run natively on ts1; every failure filed; also native: fs-couple `--test modal_aperture`, `--test material_supported_strings`, the bank chart/spatial and supported-strings unit tests, `cavity_modal`/`cavity_phs::piston` tests, fs-solid `--test membrane`, examples `textile_sheet`/`textile_cavity`, engine-piano `audio_rate`/`scale_preparation`; the registrar re-pins frankenscipy in constellation.lock to a pushed 0.4.0 head via xtask and ports piano-physics.yml's `make.bat` workaround so instrument-materials.yml is green once |
| M7.1-W1 | Commit grader + reproduce audit row | S | — | `tools/instrgrade/` holds the pack's instrgrade.py, grade_run.py, summarize.py, radiate.py (DIAGNOSTIC-ONLY) and uiowa-links.json (handed to M8-W6a) with no absolute paths; deps via `uv run --with numpy --with matplotlib`; clarinet-bore row reproduced within 0.1 dB; verisim.py and kitgrade.py untouched |
| M7.1-W2 | Grader blind-spot fixes + Q3 descriptors | M | M7.1-W1 | §7.2 list fixed with unit tests on synthetic signals (double decay, beats, stiff partials, 32 Hz fundamental, pitch-synchronous noise) |
| M7.1-W3 | Discrimination-equivalence listening harness | M | M7.1-W1 | balanced 2AFC equivalence procedure with band (0.40, 0.60): both confidence bounds inside, mixed-effects listener/item analysis, valid stopping rule; implemented on the listening-pack page; simulated listeners recover known accuracies including reversed-label (< 0.40) failures |
| M7.1-W4 | Phrase fixtures + phrase grading (Q5) | M | M7.1-W2, M8-W6a | per family: frozen phrase references and performance fixtures (scales, slurs, articulation, dynamics changes); phrase grading scores note transitions, not only per-note descriptors; validated on real-vs-real phrase pairs |
| M7.4-W1 | Wire fs-modalid into QP comparisons | S | — | modal ID of a simulated FRF recovers known modes; used in one QP test |
| M7.5-W1 | fs-psycho descriptors in grader | S | M7.1-W1 | loudness/sharpness/roughness/TNR in grader JSON; ISO signal check |
| M5.5-W1 | Convergence-ladder harness on fs-ladder | M | — | ≥ 3 rungs; per-axis declared error model (observed order where algebraic, GCI only where its assumptions hold, declared bound otherwise); axis liveness from logged dimensions plus a sensitivity fixture; coupled groups; joint run on one existing render path |
| M5.5-W2 | Simplification-equivalence comparison (§1.4) on the ladder seam | M | M5.5-W1 | detailed vs simplified fixture (tone hole with and without undercut) compared over a predeclared domain on raw calibrated output; reports a difference bound; refuses when the bound exceeds tolerance; no new harness |

### E.2 Media and materials

| ID | Title | Size | Deps | DONE-WHEN |
|---|---|---|---|---|
| M0.1-W1 | Wire ThermochemicalGasModel | S | — | γ(T), cp(T) reach duct/radiation/`cavity_modal` consumers; JANAF check |
| M0.1-W2 | Relaxation absorption for mixtures | M | M0.1-W1 | ISO 9613 air; published CO₂ data |
| M0.1-W3 | Axial gas fields + played breath field | M | M0.1-W1 | field input to ducts; breath-field generator on `fs_thermal::enthalpy::advection::Advected` plus a species extension of the same observer (no second transport kernel); cell ladder with observed order; register-interval shift on a test bore |
| M0.2-W1 | FPL ratio → orthotropic tensor adapter | S | — | wood packs resolve to E1..3, G, ν; round-trip test |
| M0.2-W2a | Passive tensor Prony law | M | — | PSD branches in Mandel coordinates; passivity under random multiaxial histories; isotropic reduction |
| M0.2-W2b | Tensor lower_to_prony | M | M0.2-W2a | certified band error vs frequency-domain law |
| M0.2-W3 | Prony tensor in plate/shell modal loss (one owner) | M | M0.2-W2b | `plate_viscosity.rs` consumes M0.2-W2b branches; `thin_plate/material_loss.rs` and reed `modal/viscosity.rs` routed through it (refusal removed; isotropic bitwise parity); native run |
| M0.2-W4 | Temperature/moisture shifts | M | M0.2-W2a, M8-W2c | WLF/Arrhenius + sorption; monotone tests |
| M0.2-W5 | Certified modal damping | M | M0.2-W2a | non-proportionality bound; complex-mode fallback triggers on a clustered test; reference = the `plate_viscosity` full modal matrix; no new projection |
| M0.3-W1a | Ogden + HGO laws | M | — | Treloar; fibre-reinforced uniaxial tests |
| M0.3-W1b | Finite-strain Prony branches | M | M0.3-W1a | relaxation tests; thermodynamic consistency |
| M0.3-W3 | Active-stress tissue law | M | M0.3-W1a | prescribed activation produces published prestress curves |
| M0.4-W1 | Third-body interfacial-layer law | M | M8-W2b | shear strength vs T across Tg; packs for rosin and a brake film |
| M0.4-W2 | String polymer laws | M | M0.2-W2a, M8-W2b, M0.2-W4 | nylon/fluorocarbon/gut packs; T/RH shifts |
| M0.4-W3 | Felt rate dependence + coupon path | M | M8-W2e | Prony-on-envelope; energy exactness kept |
| M0.5-W1 | Pack resolution on render paths | M | M8-W2a, M8-W2b, M0.2-W1, F1 | copied constants replaced; fixture packs for tests |

### E.3 Structures

| ID | Title | Size | Deps | DONE-WHEN |
|---|---|---|---|---|
| M1.1-W1a | Move modal ZOH core to acyclic owner | M | OD10 | callers unchanged; piano WAV byte-identical |
| M1.1-W1b | Create fs-string; lift specimen + sections | M | M1.1-W1a, OD10 | CONTRACT (10 sections), layer, doc-facts; tests moved and green |
| M1.1-W1c | Lift linear/stretching/longitudinal/polarization | M | M1.1-W1b, bead:5c0sj, bead:stxqd, bead:t0odx, bead:i5yei | grand_piano consumes fs-string; byte-identical on six fp fixtures + music_lane string fixture; registry refs updated |
| M1.1-W1d | engine-piano off #[path] | S | M1.1-W1c | all 11 `#[path]` refmods (`refmods.rs:21–46`, incl. stretch_tuning.rs from b1efe209e) become crate deps or are listed with owners (felt → M4.4-W3e, rest → L-2); wasm builds; clone, kernel_parts, audio_rate, scale_preparation and validate pass natively; byte-identical on LivePiano 1×, PianoAudio 2×/4× and a scale-prepared fixture; `tests/piano-oversampling.test.mjs` green |
| M1.1-W2 | Exciter-agnostic reciprocal string–body bank | M | M1.1-W1c, M5.1-W0, F1 | `git mv` `bowed_string/runtime/reciprocal/bank/**` and the string storage of `render/plate/impact/{string,supported,linear/wire}.rs` into fs-string (fs-couple keeps only the ImpactSystem host); two planes; torsion exchanges work with a termination rotation (nonzero torque row, ledger closed); N strings per termination with distinct footprints; parity with the spatial/viscous bank constructors; linear parts ZOH, not midpoint |
| M1.1-W3 | Exact linear + SAV nonlinear string stepping | M | M1.1-W1c, M5.2-W1, M1.1-W2 | no pitch warping beyond ladder tolerance at 5 kHz/48 kHz; KC glide vs Duffing |
| M1.1-W4a | Nonuniform FE of rod equations | M | M1.1-W1b | MMS order; tapered string vs analytic |
| M1.1-W4b | Geometrically exact rod dynamics | M | M1.1-W4a | large-amplitude benchmarks; agreement with modal path in small amplitude |
| M1.1-W5a | Resolved wound-string reference specimen | L | M1.3-W1a, M4.1-W1, M4.2-W2, M8-W2b, M1.3-W2a | core+wrap solid with contact and prestress; effective EI, μ, loss vs maker data |
| M1.1-W5b | Homogenized wound section with enrichment path | M | M1.1-W5a, M1.1-W4a | matches W5a within ladder tolerance; inter-turn loss |
| M1.1-W6 | Reciprocal 6-DOF termination ports | M | M1.1-W2 | work closure; rocking bridge fixture |
| M1.1-W7 | BR-01: geometric axial projection | S | M1.1-W6 | 0.06 replaced by computed projection; phantom partials re-measured |
| M1.1-W8 | String air-drag operator | M | M1.1-W1b, M3.2-W1 | BL solution KC ≲ 1; resolved characterization beyond; model-form error recorded |
| M1.1-W9 | Distributed string contact interfaces | M | M1.1-W4a, M4.1-W1, M1.1-W2 | fret/fingerboard/finger/neighbour contact along length; energy ledger; extends `render/plate/impact/linear/wire.rs` line contact and lifts grand_piano `dampers.rs`/`hammer_footprint.rs`; no second line-contact implementation |
| M1.2-W1a | Shear-deformable curved shell + static benchmarks | M | — | Scordelis–Lo, pinched cylinder, hemisphere within published tolerance |
| M1.2-W1b | Shell vibration + locking benchmarks | M | M1.2-W1a | Leissa 1973; Chapelle–Bathe asymptotic tests |
| M1.2-W1c | 3-D p-FEM thin-wall cross-check | M | M1.3-W1a, M1.2-W1a, M0.2-W1 | agreement with W1a on a spruce plate above 2 kHz shows the shear effect |
| M1.2-W2a | Work-conjugate ties/mortar between surfaces | M | M1.2-W1a | T-junction patch; virtual work preserved |
| M1.2-W2b | Closed boxes + glue compliance | M | M1.2-W2a | box modes vs monolithic; fluid-boundary ownership audit |
| M1.2-W3 | Shell–beam and shell–solid ties | M | M1.2-W2a, M1.3-W1a | brace and soundpost fixtures |
| M1.2-W4 | Compliant support/holding ports | M | M4.1-W1, M4.4-W1 | holding contact changes modes as measured on a test plate |
| M1.2-W5 | M0.2 damping in shell assembly | M | M0.2-W3, M0.2-W4, M1.2-W1a | orthotropic loss per mode vs energy quotient |
| M1.2-W6 | Outline mesher with holes + arching generators | M | — | f-hole/soundhole plates mesh; arching surface from parameters |
| M1.2-W7 | Laminate (layered) sections | M | M0.2-W2a, M1.2-W1a | varnish layer changes loss and frequency per Schelleng 1968 trend |
| M1.2-W8 | Prestressed sheet membranes, one owner (optional consolidation; consumer: §6.13 erhu candidate) | M | F1, M1.2-W1a, OD10, M0.2-W3 | native tests green; one membrane energy (§5.1: fs-solid::membrane and shell::head::nonlinear folded into the shell with D → 0); sparse eigen past the 128-DOF cap; ≥ 3-rung h-ladder vs tensioned circular-membrane Bessel modes; loss from M0.2/M0.3 |
| M1.3-W1a | TET10 elements | M | — | Cook's membrane; bending patch; DEF-08 replaced on bending paths |
| M1.3-W1b | HEX20/27 elements | M | — | beam frequencies within 0.1 % |
| M1.3-W2a | Linear orthotropic viscoelastic elements | M | M1.3-W1a, M0.2-W2b | Prony ladder certified |
| M1.3-W2b | Finite-strain visco-hyperelastic mixed (u,p) | M | M1.3-W1a, M0.3-W1b | Treloar; incompressibility patch |
| M1.3-W3 | Solid dynamics per §3.3 | M | M1.3-W1a, M5.2-W1, M1.1-W1a | exact modal linear; energy-consistent nonlinear; no algorithmic damping |
| M1.3-W4 | Volumetric meshing + gmsh/STEP import | M | — | bar and bridge meshes from descriptions and imports |
| M1.4-W1 | Prestressed assembled equilibrium | M | M1.2-W2a, M1.1-W6, M1.2-W3 | tension, downbearing, soundpost compression; tangent modes shift as published; reuses the existing equilibrium→tangent→fs-modal pattern (membrane, shell preload, fs-solid prestress); no fourth prestress path |
| M1.4-W2 | Monolithic certified eigen of assembled body | M | M1.5-W1, M0.2-W5, M1.4-W1 | eigenproblem on the assembled prestressed-equilibrium tangent and constraints; loaded vs unloaded spectra compared without resetting prestress; 10⁵-DOF window with residual certificates |
| M1.4-W3 | CMS with characteristic constraint modes (optional: exploration, no milestone consumer) | M | M1.4-W2 | assembled vs monolithic within ladder |
| M1.5-W1 | Nested dissection + window slicing | M | — | 10⁵–10⁶-DOF fixtures; timing receipt |
| M1.5-W2 | Nonlinear eigen (Beyn/NLEIGS) | L | M1.5-W1 | BEM-loaded plate modes with residual certificates |

### E.4 Acoustics

| ID | Title | Size | Deps | DONE-WHEN |
|---|---|---|---|---|
| M2.1-W1 | Levine–Schwinger thin-walled termination | M | — | tables within 1e-3 rel, ka 0–3.8 |
| M2.1-W2 | Baked full-band terminations in lines | M | M2.5-W1, M2.2-W2 | bell termination across band in driving_point/characteristic lines |
| M2.1-W3a | Plane-wave Riccati smooth profiles + cone slice ladder | M | — | exponential/Bessel horns vs analytic; slice ladder converges |
| M2.1-W3b | Continuous multimodal Riccati (m = 0) | M | M2.1-W3a | Pagneux cases; no staircase |
| M2.1-W3c | Azimuthal m ≠ 0 bases | M | M2.1-W3b | asymmetric feature test |
| M2.1-W4a | Feature embedding via M2.4 multiports | M | M2.4-W2, M2.1-W3c, M2.1-W9, M2.4-W1c | tone hole with pad and finger vs Lefebvre & Scavone 2012 |
| M2.1-W4b | Key-mechanism dynamics + moving pads | M | M2.1-W4a, M4.1-W1, M2.2-W3d, M5.3-W2b | key click emerges; energy ledger |
| M2.1-W5 | Exact Bessel everywhere | S | — | AllRegime interpolation removed; regression |
| M2.1-W6 | Multimodal holes + mouth radiation matrix | M | M2.1-W3b, M2.5-W3a, M2.1-W3c | holes in multimodal network |
| M2.1-W7 | Continuous gas fields | M | M0.1-W3, M2.1-W3a | axial T field shifts resonances as analytic |
| M2.1-W8 | Generic side-branch networks | S | — | branched tract fixture |
| M2.1-W9 | Stable scattering composition | S | — | ill-conditioned TMM case fixed |
| M2.2-W1 | Adjudicate line conventions (DEF-01) | S | — | causality tests; code changed only if failed |
| M2.2-W2 | Auto FFT sizing + compliance/pole/delay states | M | M2.2-W1 | wrap-around < tolerance; closed-volume compliance retained; DEF-02/14 |
| M2.2-W3a | Lossless variable-area DG | M | — | energy ledger 1e-12; order on cylinder/cone IRs |
| M2.2-W3b | Passive ZK pole approximations | M | M2.2-W3a, M2.1-W5 | pole ladder; impedance peaks within 1 c of Bessel incl. narrow tubes |
| M2.2-W3c | Side branches + baked radiation boundary | M | M2.2-W3b, M2.1-W2 | lattice impedance within 1 c of M2.1 on 3 fingerings |
| M2.2-W3d | Conservative moving geometry + memory transfer | M | M2.2-W3c | GCL; slow-opening hole reproduces quasi-static impedance |
| M2.2-W3e | Valve-port coupling | M | M2.2-W3c, M3.1-W1a | reed loop speaks; matches characteristic-line render within ladder |
| M2.2-W3f | Multimodal time-domain propagation | M | M2.2-W3b, M2.1-W3c, M2.1-W6 | causal multiport response matches M2.1 and M2.4 over the declared band under transverse-basis and time refinement; propagating and evanescent modes, losses and port work retained; nonlinear-to-multimodal interface error bounded |
| M2.2-W4 | Remove LC-ladder default | S | M2.2-W2 | realize_assembly default replaced; regression |
| M2.3-W1 | One-way Burgers diagnostic | M | M8-W3b | Burgers solutions; X-Br-C1 run |
| M2.3-W2 | Entropy-stable nonlinear flux + shock capturing | M | M2.2-W3b | numerical ≪ physical dissipation at mf; shock distance |
| M2.4-W0 | Interior acoustics home + design | S | OD10 | home decided (fs-feec preferred); 0ja4 scope recorded |
| M2.4-W1a | P2 Helmholtz FE | M | M2.4-W0, M1.3-W4 | rigid-box modes within 1e-4 |
| M2.4-W1b | Viscothermal BL admittance | M | M2.4-W1a | damped cylinder vs Bessel |
| M2.4-W1c | Thin-gap LNS elements | M | M2.4-W1a | narrow-slit damping vs analytic |
| M2.4-W1d | Axisymmetric FE + PML | M | M2.4-W1a, M2.1-W1 | unflanged pipe vs Levine–Schwinger; bell fixture |
| M2.4-W2 | Multiport extraction with evanescent modes | M | M2.4-W1b | tone-hole parameters vs Lefebvre FEM |
| M2.4-W3a | Shell/solid–cavity coupling incl. bar–resonator gaps | M | M2.4-W1b, M1.2-W2a, M1.2-W2b, M1.3-W1a, DEF-09 | reciprocal pressure/normal-velocity work verified on shell and solid interfaces; violin A0/A1 trend; bar–tube gap fixture; energy closure; `SharedCavity` reproduced as the uniform-pressure limit; its refusal band computed and enforced in `cavity_modal` |
| M2.4-W3b | Aperture-to-exterior coupling | M | M2.4-W1a, M2.5-W3b | f-hole radiation; through `SharedCavity::step` exterior neck pressure; `RadiationAtResonance` retired from acceptance |
| M2.5-W1 | Band-complete unbaffled bakes | M | — | band and geometry declared before baking; operator tolerances met throughout the band incl. angular/azimuthal content and termination reference planes; dense cap before convergence = blocked and routed to M2.5-W5; partial-band artifacts carry explicit band refusal |
| M2.5-W2 | One closed-form piston | S | — | three copies replaced; Bessel–Struve check; `baffled_pressure.rs` disk quadrature and `radiation.rs` compact load verified and consolidated |
| M2.5-W3a | Multi-aperture Z_rad matrix | M | M2.5-W2 | two-piston mutual impedance vs literature |
| M2.5-W3b | Passive delay-extracted rational realization | M | M2.5-W3a | Hamiltonian passivity; order ladder; DEF-03 removed; closed-loop audit (complex response, loaded modes and decay, sound over the admitted domain); the compact in-loop load (`radiation.rs`) is not a W3b realization |
| M2.5-W3c | In-loop radiation coupling (X-Cl-R2) | M | M2.5-W3b, M2.2-W2, M8-W3a | Cl-R2 falsifier run; extends BaffledEnvironment (mutual reaction absent today) |
| M2.5-W4a | Radiation loading for solids | M | — | sphere/plate loading vs analytic; planar baffled case partially delivered (973f27f43, piano example) |
| M2.5-W4b | Radiation loading for shells | M | M1.2-W1a, M2.5-W4a | loaded shell modes |
| M2.5-W5a | Complex kernels in fs-fmm | M | — | Helmholtz kernel accuracy vs direct |
| M2.5-W5b | Directional Chebyshev FMM | L | M2.5-W5a | high-kD tolerance ladder |
| M2.5-W5c | GMRES + preconditioner for BEM | M | M2.5-W5a | 64k panels solve; time/memory receipt |
| M2.5-W6 | Causal-FIR observation transfers + CQ reference | M | — | FFT-length ladder; CQ agreement |
| M2.5-W7 | Open thin shells + impedance BCs | M | — | thin disc and impedance sphere benchmarks |
| M2.5-W8 | Compact hole arrays | M | M2.5-W3a | retarded multipole ladder vs full BEM |
| M2.5-W9 | Player/chair scatterers | M | M2.5-W5c | directivity change with a body mesh |
| M2.5-W10 | Solver-route qualification: open finite-wall tubes/resonators | M | M2.5-W4a | the selected route (closed-surface BEM with the finite wall, or M2.4 interior/axisymmetric FE) reproduces a finite-wall open tube's radiation impedance and resonances with viscothermal loss within tolerance; the selected route's own prerequisites (e.g. M2.5-W7 impedance walls or M2.4-W1b) and qualification are stated explicitly; cost receipt |
| M2.6-W0 | Radiated reed observation (DEF-12) | S | F1 | BaffledPressure qualified once (analytic baffled piston on-axis plus 2J1(x)/x directivity; ring/angle/band ladder), then each adapter (ReedPerformance baffled_outlet, TwoSidedPressure environment, AperturePerformance); a reaction-load/receiver gas mismatch refuses; default for flanged development fixtures only, real terminations use the matching receiver/load route; legacy listener_m retired from acceptance (kept as X-Obs1 diagnostic); coordinated with the remote lane |
| M2.6-W1 | Generic p/∇p receivers + capsule data | M | M0.1-W2, M2.5-W6 | cardioid proximity effect vs analytic point source |
| M2.6-W2 | Recording-chain emulation from measured floors | S | M7.1-W1, M8-W6a | floor from ambient files and gaps |
| M2.6-W3 | Observation prior + marginalization | S | M7.1-W1, M2.6-W1, M2.5-W9 | aleatory observation prior marginalized (median) with a convergence-sized ensemble; epistemic unknowns reported as correlated envelopes, not medians |

### E.5 Flow and excitation

| ID | Title | Size | Deps | DONE-WHEN |
|---|---|---|---|---|
| M3.1-W1a | Valve pattern with DynamicAperture first | M | M5.1-W0 | existing reed render bitwise; seeded from `LoadedAperture`+`TwoSidedAperture` (two_sided.rs:60, 245); Single, Paired and Modal hosts reproduced bitwise; `paired.rs`/`uniform.rs` condensation merged into one owner; one `with_plate_forces` |
| M3.1-W1b | Multimode distributed traction (verify + ladder the existing ModalAperture) | M | M3.1-W1a, F1, M5.2-W1, M5.1-W2, M5.5-W1 | partially delivered, unverified (d66962bfa, 71bf3b183, 4e3d1de55, 70ccffa3b, 39f15c202): signed multimode basis, uniform-Δp traction, swept volume, full C, two-sided coupling; remaining: native run of `--test modal_aperture`, spatially varying channel traction, distributed lay contact, frequency-matched modal propagator (DEF-04), ≥ 3-rung mode/h/Δt ladders through TwoSidedAperture with pitch, threshold and regime reported |
| M3.1-W1c | Striking direction from geometry | M | M3.1-W1b | inward/outward fixtures show the right threshold shift sign; raw signed A/g can express either striking sense but the material path fixes the frame toward opening; derive the pressure-face sign from geometry |
| M3.1-W1d | N blades + contact | M | M3.1-W1c, M4.1-W1, M8-W2a | two independent blades exist (Paired: spatial contact; Modal: one uniform contact); remaining: distributed blade–blade/lay contact; paired-blade quasistatic curve vs Almeida 2007 |
| M3.1-W2 | Brass lip and glottis onto pattern | M | M3.1-W1a, M4.1-W1, M3.1-W3 | lip collision + swept flow; old loops' parity then removal (OD10) |
| M3.1-W3 | Rung-A gap-flow law (DEF-13) | M | M3.1-W1a | laminar term, inertance, contraction; regular Jacobian at closure; 70ccffa3b (structural C) does not address it; `moving.rs:13` remains inviscid |
| M3.1-W4 | Staple/channel networks | M | M3.1-W1d, M2.4-W2, M2.4-W1c, M3.1-W3 | pressure recovery vs Almeida 2007 |
| M3.1-W5 | Valve continuum from any S1 description (cane) | M | M3.1-W1b, M1.3-W2a, M0.2-W1, M0.2-W4, M8-W2a, M8-W3a | partially delivered, unverified (4e3d1de55; 39f15c202 adds tensor loss) from 1–2 DKT ResolvedPlateCharts; remaining: shell/solid sources, wet-cane tensor, vamp-map fixture, mode/h ladder |
| M3.1-W6 | Lips fixture (layered active tissue, rim) | M | M3.1-W5, M1.3-W2b, M0.3-W3, M4.1-W1, M8-W2d, M3.1-W2, M8-W5 | open-area waveforms vs Bromage 2010 |
| M3.1-W7 | Player airway port | M | M3.1-W1a, M2.1-W8, M8-W1 | tract impedance moves pitch as Chen et al. measured |
| M3.1-W8 | Rung-B resolved channel flow | L | M3.3-W1, M3.1-W1b, M3.1-W3, M4.1-W1 | moving-wall no-slip/traction exchange, swept mass and wall work close conservatively through opening, contact and reopening; coupled oscillation and thresholds converge; model-form error map vs rung A; prescribed-channel flow alone does not close it |
| M3.2-W1 | D3Q27 + collision-dissipation ladder | M | — | rotational invariance test; dissipation ladder |
| M3.2-W2a | 3-D multi-block refinement | M | M3.2-W1 | refined jet matches uniform fine grid |
| M3.2-W2b | Transparent-boundary validation | M | M3.2-W2a | acoustic pulse reflection < −40 dB |
| M3.2-W3a | Coupled compressible Navier–Stokes–Fourier dynamics | L | M3.2-W1 | EOS, total energy, pressure work, viscous heating and heat flux implemented (thermal LBM or hybrid FV energy solver); energy budget closes on a closed-box test |
| M3.2-W3b | Compressible/thermal qualification | M | M3.2-W3a | linear acoustics (dispersion/damping vs linearized NS), nonlinear manufactured solutions and thermal-flow cases pass; finite-Mach and lattice defects quantified at fixed physical Re, Ma and Pr |
| M3.2-W4 | Jet/edge-tone campaign at playing Re | L | M3.2-W2b, M3.2-W3b | self-similarity; Brown stages; grid ladder (trj/ts1) |
| M3.2-W5 | SGS as vanishing rungs (optional: exploration, no milestone consumer) | M | M3.2-W1 | SGS contribution → 0 under refinement |
| M3.3-W1 | Characteristic interface zone ↔ DG duct | L | M3.2-W2b, M2.2-W3b, M5.2-W3, M5.1-W1d, M3.2-W3b, M2.1-W3c | reflection < −40 dB; interface location/basis ladder |
| M3.3-W2a | Recorder/organ-pipe resolved zone | M | M3.3-W1, M6.1-W1a, M8-W3h | grid ladder on jet profile |
| M3.3-W2b | Threshold/regime sweep (X-Rc1) | M | M3.3-W2a | thresholds inside Auvray 2012 bands |
| M3.3-W2c | DNS comparison | M | M3.3-W2a | spectra vs Giordano/Yokoyama |
| M3.3-W3 | Resolved aperture/edge zones (tone holes, vents, bell and pipe rims) | L | M3.3-W1, M2.2-W3c | zone-extent ladder until omitted structures and outside nonlinearity are below budget; open-end nonlinear losses vs Atig, Dalmont & Gilbert 2004; register-transition effect measured on a test bore |
| M3.4-W1 | Jet-drive acceleration (optional: exploration, no milestone consumer) | M | M3.2-W4, M3.3-W2a | error vs resolved across domain |
| M3.4-W2 | Contraction acceleration (optional: exploration, no milestone consumer) | M | M3.1-W8 | error vs resolved across domain |
| M3.4-W3 | In-loop flow-noise acceleration (exploration only) | M | M3.2-W4, M3.1-W8, M7.1-W2 | cyclostationary coherence vs references; never added to output |

### E.6 Contact and friction

| ID | Title | Size | Deps | DONE-WHEN |
|---|---|---|---|---|
| M4.1-W1 | Extended normal-contact binding (DEF-07) | M | — | compliant potentials; refusals propagate; one weight convention; Hertz/HC quantitative |
| M4.2-W1 | TangentialLaw trait (DEF-06) | M | — | pHS state + dissipative port; continuity at v = vs |
| M4.2-W2 | Cattaneo–Mindlin compliance | M | M4.2-W1, M4.1-W1 | Mindlin analytic loops |
| M4.2-W3a | Third-body layer in contact | M | M4.2-W1, M0.4-W1, M4.1-W1 | layer shear at contact T; energy ledger |
| M4.2-W3b | Moving-heat-source conduction into bodies | M | M4.2-W3a | into fs-thermal networks (no new store); flash-temperature analytic; Smith & Woodhouse 2000 loops as benchmark |
| M4.2-W4 | Benchmark baselines (rate-and-state, LuGre) | S | M4.2-W1 | test-only; classic tests |
| M4.3-W1 | Nonsmooth generalized-α comparator (diagnostic only) | M | bead:ext-contact-nonsmooth-lane-oh0i, M4.1-W1 | rigid-limit impact law derived from the compliant reference and disclosed; block-on-incline; never an acceptance integrator |
| M4.3-W2 | Root selection + Δt ladders + stick–slip fixture | S | M4.2-W1 | disclosed rule; slip fraction ≈ β; Schelleng limits |
| M4.4-W1 | Driven-body interface | M | M5.1-W0, M4.1-W1 | prescribed drive reproduces point-force results |
| M4.4-W2 | Ribbon–rod distributed line contact | M | M4.4-W1, M4.2-W1, M1.1-W4a, M5.1-W3, M8-W1 | width ladder; Pitteroff–Woodhouse trend |
| M4.4-W3a | Bow fixture (stick rod + hair ribbon + rosin; feeds X-Bw0) | M | M4.4-W2, M4.2-W3b, M4.3-W2 | Helmholtz motion inside Schelleng range |
| M4.4-W3b | Finger fixture (tissue + nail) | M | M4.4-W1, M1.3-W2b, M0.3-W1b, M8-W2d | pluck release vs Pavlidou 1997 |
| M4.4-W3c | Plectrum fixture | S | M4.4-W1 | release angle/force trend |
| M4.4-W3d | Mallet fixtures (layered heads) | M | M4.4-W1, M0.4-W3, M1.3-W2a, M8-W2e, M1.3-W2b | contact time vs Chaigne & Doutaut 1997 |
| M4.4-W3e | Unified felt hammer | S | M4.4-W1, M0.4-W3 | piano parity; three adapters removed (OD10) |
| M4.4-W3f | Tongue fixture | S | M4.4-W3b, M3.1-W1a | tongued onset fixture |
| M4.5-W1 | Friction-heat partition | M | M4.2-W3b | into fs-thermal networks (no new store); heat to thermal networks; T feedback |

### E.7 Runtime and descriptions

| ID | Title | Size | Deps | DONE-WHEN |
|---|---|---|---|---|
| M5.1-W0 | Port-schema adoption contract + adapters | M | F1 | five existing modules, including `CoupledModalSystem`/`SharedCavity` and `cavity_phs::piston`, expose PortSchema; audit passes; no new N = 2-specific coupling types |
| M5.1-W1a | Graph builder (reed–bore bitwise) | M | M5.1-W0, M5.1-W2 | bitwise reproduction |
| M5.1-W1b | Sparse block assembly + analytic tangents | M | M5.1-W1a | agrees with FD ≤ 1e-8; Schur/matrix-free path |
| M5.1-W1d | IQN-ILS strategy | M | M5.1-W1b | partitioned = monolithic within tolerance |
| M5.1-W2 | Refusal certificates (DEF-05) | S | — | stalled-Newton fixture refuses at all three sites |
| M5.1-W3 | Modulated interconnection + mortar | M | M5.1-W1b | moving bow position conserves energy; nonmatching meshes preserve work |
| M5.2-W1 | Integrator doctrine + warping audit (DEF-04) | M | — | every render path listed with integrator and measured warping; fs-phs linear parts exact; audit includes `CoupledModalSystem` connections, `ImpactSystem` (`render/plate/impact.rs`, only the 0.9π guard at :276), the reed valves and the live piano at each oversampling factor |
| M5.2-W2 | Broadband modal pressure transfer | M | — | off-resonance transfer exact vs analytic |
| M5.2-W3 | Energy-conserving multi-rate | M | M5.2-W1, M5.1-W1b | energy closure across clocks; anti-aliased decimation; admission rules; reuses `pcm_wav::decimate`; policies counted in mechanical steps (piano voice retirement) are converted to physical time before any factor counts as a Δt rung |
| M5.3-W2a | Typed-path gesture addressing | M | M6.1-W1a | gestures bind EntityId; old targets as aliases |
| M5.3-W2b | Conservative topology events | M | M5.3-W2a, M5.1-W1b | trapped volume, contact, memory, heat retained; projection error reported |
| M5.4-W1 | Open voice trait + generic controls | M | M5.1-W1a | existing voices ported; parity |
| M5.4-W2 | Renderer from description + CLI | M | M5.4-W1, M6.1-W1a, M6.1-W1b, M6.1-W1c | `music_render instrument` renders the R13 description |
| M5.6-W1 | OrbitProblem for graphs | M | M5.1-W1a, M5.5-W1 | thresholds match time-domain sweeps on reed fixture |
| M6.1-W1a | Assembly Description schema + migration | M | — | scenario.ir migration test; check-schemas passes |
| M6.1-W1b | Validity-domain validator | M | M6.1-W1a | out-of-domain fixtures refuse with typed errors |
| M6.1-W1c | R13 clarinet description | M | M6.1-W1a, M8-W3a, M6.1-W1b | validates against module domains; converts to the existing wind-CLI inputs and renders |
| M6.1-W2 | Generators (bore, outline+holes, arching) | M | M6.1-W1a | generated parts validate and retain source dimensions/uncertainty under §1.4; generic shapes are explicitly fixtures |
| M6.2-W1 | Performance schema + MIDI + dynamics + per-note rule | M | M6.1-W1a | MIDI file → gesture tracks; declared objectives recorded |
| M6.2-W2 | Player models from measured data | M | M6.2-W1, M8-W4 | dynamics tables from literature |
| M6.3-W1 | Geometry ledger + representation-aware import checks | S | M8-W1 | three fixtures (a closed solid, a midsurface with declared open edges and thickness, a part with a real opening) import through fs-io; closed solids watertight/manifold, midsurface edges and thickness declared, real openings preserved; SI units, frames and part/material mapping checked; ledger rows link source, licence and uncertainty |
| M6.3-W2 | Blender headless parametric reconstruction (route-conditional: required by specimens that use Blender reconstruction) | M | M6.3-W1 | `blender -b -P` scripts regenerate a parametric bore and a plate outline byte-identically from a dimension table under a recorded, fixed Blender version, export settings and seed; SI exports hashed |
| M6.3-W3 | Photo calibration + per-surface uncertainty (route-conditional: required by specimens that use photo reconstruction) | M | M6.3-W1, M6.1-W1a | calibrated multi-view reconstruction of a test object with any qualified photogrammetry tool (Blender, via M6.3-W2, only when the route uses it); per-surface uncertainty map; occluded regions flagged as bounds |
| M6.3-W4 | Selected acquisition-route qualification | S | M6.3-W1 | the route each specimen actually uses (measured CAD/CT/scan import via W1, or Blender/photo via W2/W3) reproduces an independently known object within declared tolerance (bore diameters ≤ 0.1 mm, outlines ≤ 0.5 mm, or tighter where measured); W2/W3 become prerequisites only for specimens that use that route |
| M6.4-W1 | Anatomy parts in the Assembly Description (articulators, postures, contacts) | M | M6.1-W1a, M6.1-W1b | articulators as parts of the existing description schema with postures and contact surfaces (no separate registry); migration test; validator refuses missing tissue properties |
| M6.4-W2 | Airway from MRI segmentation | M | M6.4-W1, M8-W3k, M2.1-W8, M8-W4 | airway built from licensed MRI with declared playing-posture uncertainty; computed tract impedance inside the measured population spread (Fritz & Wolfe 2005; Chen et al.); DVTD printed-tract acoustics reproduced as a passive check |
| M6.4-W3 | Lip/embouchure geometry and contact surfaces | M | M6.4-W1, M8-W2d, M1.3-W2b, M3.1-W6 | contact surfaces for reed, brass rim, flute and pan-pipe rims, duduk reed; lips fixture open-area trend vs Bromage 2010 |
| M6.4-W4 | Hand/finger geometry for holes, stops, holding | M | M6.4-W1, M4.4-W3b, M2.1-W4a, M1.2-W4 | finger covering and half-holing changes a test hole impedance as measured; stop and holding contacts validate |

### E.8 Data

Every M8-W3 geometry item below also satisfies the §8 geometry completion rule. Its size estimate covers a first declared specimen/configuration, not exhaustive anatomy or all instrument variants; split acquisition work if the actual source requires it. Missing hidden geometry or material evidence stays a named blocker to affected fidelity claims.

| ID | Title | Size | Deps | DONE-WHEN |
|---|---|---|---|---|
| M8-W1 | Data ledger | S | — | ledger file + schema (source, licence, uncertainty, digitization, decision); every later data item writes rows; consumed by M6.3-W1 |
| M8-W2a | Woodwind materials | M | M8-W1 | packs: wet and dry cane tensor, grenadilla/ABS, ebony, cork, pad leather/felt with state and uncertainty; consumed by M3.1-W5 and M0.5-W1 without defaults |
| M8-W2b | String materials | M | M8-W1 | packs: nylon, gut, fluorocarbon, multifilament cores, windings, horsehair, rosin third-body data; consumed by M0.4-W2 and M0.4-W1 |
| M8-W2c | Wood loss vs frequency/direction | M | M8-W1 | packs: wood loss vs frequency and direction, moisture and temperature; consumed by M0.2-W3 and M0.2-W4 |
| M8-W2d | Tissue data | S | M8-W1 | sourced layered properties, state and uncertainty for the selected lip/tongue/finger/fold fixture; geometry alone is insufficient; intrinsic constitutive data kept separate from apparatus responses (e.g. the n = 1 perioral indentation), which are used only through a forward model of the measured apparatus |
| M8-W2e | Felts, rubber, yarn, cord | S | M8-W1 | packs: felt coupons, rubber, yarn, cord; consumed by M0.4-W3 and M4.4-W3d |
| M8-W2f | Brass alloys | S | M8-W1 | packs: brass alloys vs temperature and temper; consumed by M2.5-W7 bell-wall fixture |
| M8-W3a | R13 + B4 geometry | M | M8-W1, M6.3-W1 | bore, lattice (chimneys, undercuts, pads, key mechanics), mouthpiece cavity/baffle/lay, reed vamp map, fingering chart; per the §8 completion rule; consumed by M6.1-W1c |
| M8-W3b | Bach 180S37 geometry | M | M8-W1, M6.3-W1 | mouthpiece, leadpipe, valves and slides with bends and ports, bell flare and wall; per the §8 completion rule; consumed by M-BR-QP1 |
| M8-W3c | Violin geometry | M | M8-W1, M6.3-W1 | arching, thickness maps, bracing (bass bar, soundpost), bridge, f-holes, strings; per the §8 completion rule; consumed by M-VNP-QP1; route-dependent: when Blender/photo reconstruction is used, add M6.3-W2/W3 as prerequisites at claim time |
| M8-W3c2 | Viola geometry | M | M8-W1, M6.3-W1 | as M8-W3c for a viola; consumed by M-VA-QP1; route-dependent: when Blender/photo reconstruction is used, add M6.3-W2/W3 as prerequisites at claim time |
| M8-W3c3 | Cello geometry | M | M8-W1, M6.3-W1 | as M8-W3c for a cello (+ endpin, floor contact); consumed by M-VC-QP1; route-dependent: when Blender/photo reconstruction is used, add M6.3-W2/W3 as prerequisites at claim time |
| M8-W3c4 | Double-bass geometry | M | M8-W1, M6.3-W1 | as M8-W3c for a double bass; consumed by M-DB-QP1; route-dependent: when Blender/photo reconstruction is used, add M6.3-W2/W3 as prerequisites at claim time |
| M8-W3d | Classical guitar geometry | M | M8-W1, M6.3-W1 | top with thickness map and fan bracing, back and bars, ribs/linings, bridge, frets, cavity; per the §8 completion rule; consumed by M-GT-QP1; route-dependent: when Blender/photo reconstruction is used, add M6.3-W2/W3 as prerequisites at claim time |
| M8-W3d2 | Steel-string guitar geometry | M | M8-W1, M6.3-W1 | distinct steel-string specimen: X-braced top with thickness map, back and bracing, steel-string set (gauges, tensions, windings), bridge pins and saddle, materials (spruce, rosewood/mahogany, steel/bronze); per the §8 completion rule; consumed by M-SG-QP1 |
| M8-W3e | Saxophone geometry | M | M8-W1, M6.3-W1 | conical bore, keyed holes with pads, mouthpiece chamber and lay; per the §8 completion rule; consumed by M-SX-QP1 |
| M8-W3f | Oboe/bassoon geometry | M | M8-W1, M6.3-W1 | bores, bocal, staple and double-reed blades; per the §8 completion rule; consumed by M-DR-QP1 |
| M8-W3g | Boehm flute geometry | M | M8-W1, M6.3-W1 | bore, head taper, cork cavity, embouchure hole and chimney, keyed holes; per the §8 completion rule; consumed by M-FL-QP1 |
| M8-W3h | Recorder/organ-pipe geometry | S | M8-W1, M6.3-W1 | windway, window, labium, bore; per the §8 completion rule; consumed by M3.3-W2a |
| M8-W3i | Bars (calipers) + resonators | S | M8-W1, M6.3-W1 | xylophone/marimba/vibraphone bar profiles, nodal holes, supports, tubes and gaps; per the §8 completion rule; consumed by X-Bar1 and M-BAR-QP1 |
| M8-W3j | Trombone/horn/tuba geometry | M | M8-W1, M6.3-W1 | mouthpieces, slides/valves, bells; per the §8 completion rule; consumed by M-BR-QP1, M-HN-QP1, M-TU-QP1 |
| M8-W3k | Tract MRI sets | S | M8-W1, M6.3-W1 | licensed geometry with scale, segmentation uncertainty and posture; playing anatomy distinguished from speech, and lips/tongue/airway interfaces usable by the selected modules |
| M8-W3l | Duduk + reed geometry | M | M8-W1, M6.3-W1 | bore, holes, reed blade thickness map, reed interior, cap and bridle; sources incl. Maugeais & Dalmont 2024 dimensions (no fitted parameters); per the §8 completion rule; consumed by M-DK-QP1; route-dependent: when Blender/photo reconstruction is used, add M6.3-W2/W3 as prerequisites at claim time |
| M8-W3m | Tin whistle geometry | S | M8-W1, M6.3-W1 | windway, window, labium, bore, wall, holes of the chosen specimen (calipers, casts or CT); per the §8 completion rule; consumed by M-TW-QP1 |
| M8-W3n | Siku/antara geometry | M | M8-W1, M6.3-W1 | per tube length, inner diameter, wall, node and rim shape, cane species; museum records as identity only; per the §8 completion rule; consumed by M-PP-QP1; route-dependent: when Blender/photo reconstruction is used, add M6.3-W2/W3 as prerequisites at claim time |
| M8-W4 | Player data | M | M8-W1 | tables per §8 (blowing pressure, lip force, tonguing, tract impedance, bowing, plucks, jets, mallet velocity) with sources; consumed by M6.2-W2 |
| M8-W5 | Machine-player data | S | M8-W1 | threshold and regime tables from artificial mouths, bowing machines and recorder studies; consumed by the Q1 milestones |
| M8-W6a | References manifest + fetch + Iowa equipment check | S | M8-W1, M7.1-W1 | references.json (URL, sha256, licence) + fetch script; Iowa chain facts recorded per M2.6; ambient files fetched; consumed by M7.1-W1 and M2.6-W2 |
| M8-W6b | Extra licensed references | S | M8-W1, M8-W6a | licences recorded for TU Berlin, voice, recorder and steel-string references; consumed by M-RC-Q2, M-VO-Q2, M-SG-Q2 |
| M8-W6c | Duduk/whistle/pan-pipe references | M | M8-W1, M8-W6a | licensed anechoic recordings, or owner recordings of the measured specimens in a characterized room with declared chain; consumed by M-DK-Q2, M-TW-Q2, M-PP-Q2 |
| M8-W7 | Co-characterized specimen program | M | OD4 | per approved family: measurement plan and first data for one physical specimen (geometry, materials, passive response, controlled excitation, microphones, sound); consumed by the family QP1 and Q4; measurement back-action (taps, sensor cavities, transducers) described or bounded; state drift logged over repeated takes |
| M8-W8a | Same-specimen passive pair: clarinet | M | M8-W1, M8-W3a, OD15 | geometry plus measured complex input impedance (amplitude, linewidth, phase) or admittance/modes of the same object; consumed by X-Cl-R0; measurement back-action described or bounded and drift logged as in M8-W7 |
| M8-W8b | Same-specimen passive pair: saxophone | M | M8-W1, M8-W3e, OD15 | geometry plus measured complex input impedance (amplitude, linewidth, phase) or admittance/modes of the same object; consumed by X-Dr0; measurement back-action described or bounded and drift logged as in M8-W7 |
| M8-W8c | Same-specimen passive pair: tin whistle | M | M8-W1, M8-W3m, OD15 | geometry plus measured complex input impedance (amplitude, linewidth, phase) or admittance/modes of the same object; consumed by X-Tw1; measurement back-action described or bounded and drift logged as in M8-W7 |
| M8-W8d | Same-specimen passive pair: duduk | M | M8-W1, M8-W3l, OD15 | geometry plus measured complex input impedance (amplitude, linewidth, phase) or admittance/modes of the same object; consumed by X-Dk1; measurement back-action described or bounded and drift logged as in M8-W7 |
| M8-W8e | Same-specimen passive pair: bar | M | M8-W1, M8-W3i, OD15 | geometry plus measured complex input impedance (amplitude, linewidth, phase) or admittance/modes of the same object; consumed by X-Bar1; measurement back-action described or bounded and drift logged as in M8-W7 |
| M8-W8f | Same-specimen passive pairs: classical guitar, steel-string guitar and violin bridges | M | M8-W1, M8-W3d, OD15, M8-W3c, M8-W3d2 | geometry plus measured bridge admittance and body modes of the same objects (classical guitar, steel-string guitar, violin); measurement back-action described or bounded and drift logged as in M8-W7; consumed by M-GT-QP1, M-SG-QP1 and M-VNP-QP1 |
| M8-W9 | Owner anatomy scans (optional, consented; OD16) | M | M8-W1, OD16 | intraoral scan, face/lip scan, MRI if available, with consent records; consumed by M6.4-W2/W3 |

### E.9 Instruments, experiments, breaches, live

Families: CL clarinet · SX saxophone · DR oboe/bassoon · BR trombone/trumpet · HN horn · TU tuba · BAR xylophone/marimba · VB vibraphone · GT classical guitar · SG steel-string guitar · VNP/VN violin (pizzicato/arco) · VA viola · VC cello · DB double bass · RC recorder/organ pipe · FL flute · VO voice · DK duduk · TW tin whistle · PP Andean pan pipes. Each family has QP1 (QP + Q1), Q2, Q4 (Q3 + binding Q4) and Q5 rows, except violin pizzicato (VNP): its QP1 and Q2 rows feed the violin family's shared M-VN-Q4/Q5, which report pizzicato and arco separately. The §7.1 gate expansion rule is already expanded into the dependencies.

| ID | Title | Size | Deps | DONE-WHEN |
|---|---|---|---|---|
| M-CL0 | First audible milestone (clarinet baseline) | M | F1, M7.1-W1, M2.6-W0, M2.2-W1, M2.2-W2, M8-W3a, M5.1-W2 | §6.0; flanged TMM terminal + qualified baffled_outlet; compact in-loop load (ka ≤ 0.5) recorded as a limitation |
| X-Obs1 | Matched-observer diagnostic | S | F1, M7.1-W1, M2.6-W0, M8-W3a | one accepted trajectory, two observers (legacy internal vs physical baffled outlet) with geometry/mechanism/gas fixed; observation share of the baseline gap reported; factors varied singly; seam exists: one ReedBoreVoice under Legacy or Baffled observation (`outlet.rs:49–83`; test `outlet/tests.rs:21`, unexecuted); two renders assert identical reed state |
| X-Cl-R0 | Same-specimen impedance (QP) | S | M8-W3a, M8-W8a, M2.1-W1, M2.1-W5 | peaks within 10 c; amplitudes, linewidths and phase within measurement uncertainty on a block of fingerings |
| X-Cl-R1 | pp→ff nonlinearity vs reed dynamics | S | M-CL0, M8-W4 | falsifier §6.1 |
| X-Cl-R2 | Lattice + radiation tail test | S | M2.5-W3c, M2.6-W0, M7.1-W1 | falsifier §6.1 |
| X-Cl-R3 | In-loop vs output noise (exploration) | S | M3.4-W3, M7.1-W2, M-CL0 | falsifier §6.1 |
| M-CL-QP1 | Clarinet: passive chain + speaks (QP, Q1) | M | X-Cl-R0, M-CL0, M3.1-W1a, M3.1-W3, M5.5-W1, M8-W5, X-Obs1 | QP pass on the same-specimen pair; Q1: intended register, note within ±10 c after the setup protocol, thresholds within machine-player ranges, threshold/hysteresis map |
| M-CL-Q2 | Clarinet recognisable | M | M-CL-QP1, M-CL0, M2.1-W1, M2.1-W2, M2.1-W5, M3.1-W1a, M3.1-W3, M3.1-W5, M3.1-W7, M2.5-W1, M2.5-W3c, M6.2-W1, M7.1-W2, M2.2-W4, M5.5-W1, M8-W5 | Q2 on held-out notes |
| M-CL-Q3 | Clarinet inside the ruler | M | M-CL-Q2, M2.1-W4a, M2.1-W4b, M2.2-W3c, M2.2-W3f, M3.3-W3, M0.1-W3, M2.1-W7, M3.1-W8, M0.5-W1, M2.5-W8, M4.4-W3f, M5.6-W1, M6.1-W1c, M6.4-W2, M6.4-W3, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M6.4-W4 | Q3; §1.4 description complete |
| M-CL-Q4 | Clarinet blind tests (binding Q4) | S | M-CL-Q3, M7.1-W3, M2.5-W9, M8-W7 | binding Q4 (listening + spectrogram equivalence) on frozen held-out renders |
| M-CL-Q5 | Clarinet phrases (Q5) | M | M-CL-Q4, M7.1-W4, M5.3-W2b, M2.2-W3d, M2.1-W4b | Q5 on frozen phrase fixtures |
| X-Dr0 | Saxophone same-specimen impedance | S | M8-W3e, M8-W8b, M2.1-W3a | as X-Cl-R0 |
| X-Sx1 | Airway in altissimo | S | M3.1-W7, M6.4-W2, M8-W3e | falsifier §6.2 |
| M-SX-QP1 | Saxophone: passive chain + speaks (QP, Q1) | M | X-Dr0, M8-W3e, M-CL-QP1, M5.5-W1, M8-W5 | QP vs same-specimen pair where it exists (else flagged population data); Q1: intended regime, note within ±10 c after the setup protocol, thresholds and onset within machine-player ranges, threshold/hysteresis map |
| M-SX-Q2 | Saxophone recognisable | M | M-SX-QP1, M7.1-W2, M5.5-W1, M-CL-Q2, M2.1-W3a, M8-W3e | Q2 on held-out notes |
| M-SX-Q4 | Saxophone inside the ruler + blind tests (Q3, Q4) | M | M-SX-Q2, M2.1-W6, M2.2-W3f, X-Sx1, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M3.1-W8, M8-W3e, M6.4-W2, M6.4-W3, M6.4-W4 | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete |
| M-SX-Q5 | Saxophone phrases (Q5) | M | M-SX-Q4, M7.1-W4, M5.3-W2b, M2.2-W3d, M2.1-W4b | Q5 on frozen phrase fixtures |
| X-Dr1 | Double-reed characteristic | S | M3.1-W1d, M3.1-W4 | vs Almeida 2007 incl. staple recovery |
| M-DR-QP1 | Oboe/bassoon: passive chain + speaks (QP, Q1) | M | X-Dr1, M8-W3f, M-CL-QP1, M5.5-W1, M8-W5 | QP vs same-specimen pair where it exists (else flagged population data); Q1: intended regime, note within ±10 c after the setup protocol, thresholds and onset within machine-player ranges, threshold/hysteresis map |
| M-DR-Q2 | Oboe/bassoon recognisable | M | M-DR-QP1, M7.1-W2, M5.5-W1, M-CL-Q2, M3.1-W1d, M3.1-W4, M2.4-W2, M8-W3f | Q2 on held-out notes |
| M-DR-Q4 | Oboe/bassoon inside the ruler + blind tests (Q3, Q4) | M | M-DR-Q2, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M3.1-W8, M8-W3f, M6.4-W2, M6.4-W3, M6.4-W4 | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete |
| M-DR-Q5 | Oboe/bassoon phrases (Q5) | M | M-DR-Q4, M7.1-W4, M5.3-W2b, M2.2-W3d, M2.1-W4b | Q5 on frozen phrase fixtures |
| X-Dr2 | Bore vs reed residual | S | M-DR-Q2 | recorded with the dominant residual named |
| X-Br-C1 | Burgers ff enrichment diagnostic | S | M2.3-W1, M8-W3b, M8-W4 | falsifier §6.4 |
| X-Br-C2 | Lip model dominance | S | M3.1-W2, M8-W4, M8-W3b, M7.1-W1 | falsifier §6.4 |
| X-Br-C3 | Bell beaming vs position spread | S | M2.6-W3, M2.5-W1, M8-W3b | falsifier §6.4 |
| M-BR-QP1 | Trombone/trumpet: passive chain + speaks (QP, Q1) | M | M8-W3b, M8-W3j, M2.1-W3b, M2.5-W1, M3.1-W2, M5.5-W1, M8-W5 | QP vs same-specimen pair where it exists (else flagged population data); Q1: intended regime, note within ±10 c after the setup protocol, thresholds and onset within machine-player ranges, threshold/hysteresis map |
| M-BR-Q2 | Trombone/trumpet recognisable | M | M-BR-QP1, M7.1-W2, M5.5-W1, M2.2-W3e, M2.2-W3f, M2.3-W2, M2.1-W3b, M2.1-W6, M3.1-W2, M3.1-W6, M3.1-W7, M2.5-W1, M2.5-W7, M1.2-W1a, M2.4-W1a, M2.4-W2, M0.1-W3, M5.6-W1, M8-W2f, M8-W3b, M8-W3j, M6.4-W3 | Q2 on held-out notes |
| M-BR-Q4 | Trombone/trumpet inside the ruler + blind tests (Q3, Q4) | M | M-BR-Q2, M2.5-W9, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M3.1-W8, M8-W3b, M6.4-W2, M6.4-W3, M6.4-W4 | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete |
| M-BR-Q5 | Trombone/trumpet phrases (Q5) | M | M-BR-Q4, M7.1-W4, M5.3-W2b, M2.2-W3d | Q5 on frozen phrase fixtures |
| M-HN-QP1 | Horn: passive chain + speaks (QP, Q1) | M | M-BR-QP1, M2.4-W1d, M6.4-W4, M8-W3j, M5.5-W1, M8-W5 | QP vs same-specimen pair where it exists (else flagged population data); Q1: intended regime, note within ±10 c after the setup protocol, thresholds and onset within machine-player ranges, threshold/hysteresis map |
| M-HN-Q2 | Horn recognisable | M | M-HN-QP1, M7.1-W2, M5.5-W1, M-BR-Q2, M2.4-W1d, M2.1-W3c | Q2 on held-out notes |
| M-HN-Q4 | Horn inside the ruler + blind tests (Q3, Q4) | M | M-HN-Q2, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M3.1-W8, M8-W3j, M6.4-W2, M6.4-W3, M6.4-W4 | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete |
| M-HN-Q5 | Horn phrases (Q5) | M | M-HN-Q4, M7.1-W4, M5.3-W2b, M2.2-W3d | Q5 on frozen phrase fixtures |
| M-TU-QP1 | Tuba: passive chain + speaks (QP, Q1) | M | M-BR-QP1, M8-W3j, M5.5-W1, M8-W5 | QP vs same-specimen pair where it exists (else flagged population data); Q1: intended regime, note within ±10 c after the setup protocol, thresholds and onset within machine-player ranges, threshold/hysteresis map; graded from the 2nd harmonic until the chamber cutoff is sourced |
| M-TU-Q2 | Tuba recognisable | M | M-TU-QP1, M7.1-W2, M5.5-W1, M-BR-Q2, M2.5-W5b, M2.5-W5c | Q2 on held-out notes |
| M-TU-Q4 | Tuba inside the ruler + blind tests (Q3, Q4) | M | M-TU-Q2, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M3.1-W8, M8-W3j, M6.4-W2, M6.4-W3, M6.4-W4 | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete |
| M-TU-Q5 | Tuba phrases (Q5) | M | M-TU-Q4, M7.1-W4, M5.3-W2b, M2.2-W3d | Q5 on frozen phrase fixtures |
| Bar-S1 | Bar first slice (linear3 3-D + BEM + mallet) | M | F1, M7.1-W1 | h-ladder; graded; beam cross-check |
| X-Bar1 | Bare bar component (measured force) | S | Bar-S1, M8-W3i, M8-W8e, M7.4-W1 | caliper-measured bar under impact-hammer force: modes, shapes and decay predicted; no contact or resonator claim |
| X-Bar2a | Resonator tube passive response | S | M8-W3i, M8-W8e, M2.5-W10 | tube input impedance and radiation computed vs measured |
| X-Bar2 | Assembled bar–tube coupling (tuned/detuned/removed) | S | Bar-S1, M2.5-W10, X-Bar2a, X-Bar1, M2.4-W3a | quantitative mechanical pole/decay and radiated-pressure changes predicted across tuned, detuned and removed-tube controls |
| M-BAR-QP1 | Xylophone/marimba: passive chain + speaks (QP, Q1) | M | X-Bar1, M8-W8e, M1.3-W1a, M5.5-W1, M4.1-W1, M8-W5 | QP vs same-specimen pair where it exists (else flagged population data); Q1: pitch within ±10 c after the setup protocol; attack, contact time and decay under machine-strike conditions within measured ranges |
| M-BAR-Q2 | Xylophone/marimba recognisable | M | M-BAR-QP1, M7.1-W2, M5.5-W1, Bar-S1, M1.3-W1a, M1.3-W2a, M1.3-W3, M2.5-W4a, M2.4-W1a, M4.1-W1, M4.4-W1, M4.4-W3d, M8-W3i, M0.5-W1, M5.2-W2, M8-W2c, M8-W2e, M2.5-W10, X-Bar2 | Q2 on held-out notes |
| M-BAR-Q4 | Xylophone/marimba inside the ruler + blind tests (Q3, Q4) | M | M-BAR-Q2, M1.3-W1b, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M8-W3i | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete |
| M-BAR-Q5 | Xylophone/marimba phrases (Q5) | M | M-BAR-Q4, M7.1-W4, M5.3-W2b, M4.4-W3d | Q5 on frozen phrase fixtures |
| M-VB-QP1 | Vibraphone: passive chain + speaks (QP, Q1) | M | M-BAR-QP1, M8-W3i, M5.5-W1 | QP vs same-specimen pair where it exists (else flagged population data); Q1: pitch within ±10 c after the setup protocol; attack, contact time and decay under machine-strike conditions within measured ranges |
| M-VB-Q2 | Vibraphone recognisable | M | M-VB-QP1, M7.1-W2, M5.5-W1, M-BAR-Q2, M4.1-W1, M5.3-W2b | Q2 on held-out notes |
| M-VB-Q4 | Vibraphone inside the ruler + blind tests (Q3, Q4) | M | M-VB-Q2, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M8-W3i | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete |
| M-VB-Q5 | Vibraphone phrases (Q5) | M | M-VB-Q4, M7.1-W4, M5.3-W2b, M4.4-W3d | Q5 on frozen phrase fixtures |
| Gt-S1 | Guitar first slice (plate-chart bank) | M | F1, M7.1-W1 | uses `from_viscous_material_plate_chart_with_force_ports` (`bank/chart/spatial.rs:73`) with per-string footprints; 174af92bc/2d269fd3d natively verified under F1; flat plate and compact pressure recorded as model-form rungs; Q1 + graded |
| X-Gt1 | Two-polarization decay | S | M1.1-W2, M7.1-W2, M1.1-W6, Gt-S1 | falsifier §6.6 |
| X-Gt2 | String loss and late decay | S | M7.1-W2, Gt-S1, M0.4-W2, M1.1-W8 | §6.6 |
| X-Gt3 | Finger + nail release | S | M4.4-W3b, Gt-S1, M7.1-W2, M8-W4 | §6.6 |
| X-Gt4 | A/E-string HNR causes | S | M7.1-W2, Gt-S1, M1.1-W9 | causes separated |
| M-GT-QP1 | Classical guitar: passive chain + speaks (QP, Q1) | M | M8-W3d, M8-W8f, M1.4-W2, M5.5-W1, M7.4-W1, DEF-09, M8-W5 | QP vs same-specimen admittance/modes; Q1: pitch within ±10 c after the setup protocol; attack and decay under machine-pluck (measured release) conditions within measured ranges |
| M-GT-Q2 | Classical guitar recognisable | M | M-GT-QP1, M7.1-W2, M5.5-W1, Gt-S1, M1.1-W1c, M1.1-W2, M1.1-W5b, M1.1-W6, M1.1-W9, M0.4-W2, M1.2-W1a, M1.2-W2a, M1.2-W6, M1.4-W1, M1.4-W2, M2.4-W3a, M2.4-W3b, M2.5-W4b, M4.4-W3b, M4.1-W1, M8-W3d, M0.5-W1, M1.1-W3, M1.1-W7, M1.1-W8, M1.2-W1b, M1.2-W2b, M1.2-W5, M5.2-W2, M6.1-W2, M8-W2c | Q2 on held-out notes |
| M-GT-Q4 | Classical guitar inside the ruler + blind tests (Q3, Q4) | M | M-GT-Q2, M1.1-W4b, M1.2-W1c, M1.5-W2, M2.5-W9, M6.4-W4, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M8-W3d | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete |
| M-GT-Q5 | Classical guitar phrases (Q5) | M | M-GT-Q4, M7.1-W4, M5.3-W2b, M4.4-W3b, M1.1-W9 | Q5 on frozen phrase fixtures |
| M-SG-QP1 | Steel-string guitar (after licensed reference): passive chain + speaks (QP, Q1) | M | M-GT-QP1, M5.5-W1, M8-W3d2, M8-W8f | QP vs same-specimen admittance/modes; Q1: pitch within ±10 c after the setup protocol; attack and decay under machine-pluck (measured release) conditions within measured ranges |
| M-SG-Q2 | Steel-string guitar (after licensed reference) recognisable | M | M-SG-QP1, M7.1-W2, M5.5-W1, M-GT-Q2, M4.4-W3c, M8-W6b, M8-W3d2 | Q2 on held-out notes |
| M-SG-Q4 | Steel-string guitar (after licensed reference) inside the ruler + blind tests (Q3, Q4) | M | M-SG-Q2, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M6.4-W4, M8-W3d2 | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete |
| M-SG-Q5 | Steel-string guitar (after licensed reference) phrases (Q5) | M | M-SG-Q4, M7.1-W4, M5.3-W2b, M4.4-W3c | Q5 on frozen phrase fixtures |
| M-VNP-QP1 | Violin: passive chain + pizzicato speaks (QP, Q1) | M | M8-W3c, M8-W8f, M1.4-W2, M5.5-W1, M7.4-W1 | QP vs same-specimen admittance/modes; pizzicato Q1 (pitch within ±10 c after setup protocol) |
| M-VNP-Q2 | Violin pizzicato recognisable (shared pizzicato path for the family) | L | M-VNP-QP1, M7.1-W2, M5.5-W1, M8-W3c, M1.1-W1c, M1.1-W2, M1.1-W5b, M1.1-W6, M1.1-W9, M4.4-W3b, M1.3-W1a, M1.2-W1a, M1.2-W2a, M1.2-W3, M1.2-W4, M1.2-W6, M1.2-W7, M1.4-W1, M1.4-W2, M2.4-W3a, M2.4-W3b, M2.5-W4b, M0.5-W1, M1.1-W3, M1.1-W7, M1.1-W8, M1.2-W1b, M1.2-W2b, M1.2-W5, M5.2-W2, M6.1-W2, M8-W2c | Q2 on held-out pizzicato notes |
| X-Pz1 | Violin pizzicato | S | M-VNP-Q2 | §6.7 |
| X-Bw0 | Schelleng-range bowing | S | M4.4-W3a, M8-W5, M8-W4 | §6.7 |
| X-Bw1 | Body screen | S | M-VNP-Q2, M7.4-W1 | §6.7 |
| X-Bw2 | Friction vs kinematics vs player | S | M4.4-W3a, M7.1-W2, M8-W4, M4.2-W4, X-Bw1 | falsifier §6.7 |
| X-Bw3 | Onset classes | S | M4.4-W3a, M7.1-W2, M8-W4, M8-W5 | §6.7 |
| M-VN-QP1 | Violin arco speaks (Q1) | M | M-VNP-QP1, X-Bw0, X-Bw3, M8-W5 | Helmholtz motion, slip fraction ≈ β, Schelleng limits, onset classes and pre-Helmholtz transient length within bowing-machine ranges |
| M-VN-Q2 | Violin arco recognisable | L | M-VN-QP1, M-VNP-Q2, M7.1-W2, M5.5-W1, M4.4-W2, M4.4-W3a, M4.2-W3a, M4.2-W3b, M4.3-W2, M5.3-W2b, M5.2-W3, M6.2-W1 | Q2 on held-out arco notes |
| M-VN-Q4 | Violin inside the ruler + blind tests (Q3, Q4) | M | M-VN-Q2, M1.1-W4b, M1.2-W1c, M1.5-W2, M2.5-W9, M4.3-W1, M4.5-W1, M6.4-W4, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M8-W3c, M-VNP-Q2 | Q3 and binding Q4; acceptance audio uses compliant contact (M4.3-W1 is a comparator only); pizzicato and arco each pass on their own held-out set and are reported separately |
| M-VN-Q5 | Violin phrases (Q5) | M | M-VN-Q4, M7.1-W4, M5.3-W2b, M4.4-W3a | Q5 on frozen phrase fixtures; pizzicato and arco phrase results reported separately |
| M-VA-QP1 | Viola: passive chain + speaks (QP, Q1) | M | M8-W3c2, M-VNP-QP1, M-VN-QP1, M5.5-W1 | QP vs same-specimen pair where it exists (else flagged population data); Q1: intended regime, note within ±10 c after the setup protocol, thresholds and onset within machine-player ranges, threshold/hysteresis map |
| M-VA-Q2 | Viola recognisable | M | M-VA-QP1, M7.1-W2, M5.5-W1, M-VN-Q2, M8-W3c2, M-VNP-Q2 | Q2 on held-out notes, pizzicato and arco reported separately |
| M-VA-Q4 | Viola inside the ruler + blind tests (Q3, Q4) | M | M-VA-Q2, M4.3-W1, M4.5-W1, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M8-W3c2, M6.4-W4 | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete; pizzicato and arco each pass on their own held-out set and are reported separately |
| M-VA-Q5 | Viola phrases (Q5) | M | M-VA-Q4, M7.1-W4, M5.3-W2b, M4.4-W3a | Q5 on frozen phrase fixtures; pizzicato and arco phrase results reported separately |
| M-VC-QP1 | Cello: passive chain + speaks (QP, Q1) | M | M8-W3c3, M-VNP-QP1, M-VN-QP1, M5.5-W1 | QP vs same-specimen pair where it exists (else flagged population data); Q1: intended regime, note within ±10 c after the setup protocol, thresholds and onset within machine-player ranges, threshold/hysteresis map |
| M-VC-Q2 | Cello recognisable | M | M-VC-QP1, M7.1-W2, M5.5-W1, M-VN-Q2, M8-W3c3, M2.5-W5a, M2.5-W5c, M-VNP-Q2, M2.5-W5b | Q2 on held-out notes, pizzicato and arco reported separately |
| M-VC-Q4 | Cello inside the ruler + blind tests (Q3, Q4) | M | M-VC-Q2, M4.3-W1, M4.5-W1, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M8-W3c3, M6.4-W4 | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete; pizzicato and arco each pass on their own held-out set and are reported separately |
| M-VC-Q5 | Cello phrases (Q5) | M | M-VC-Q4, M7.1-W4, M5.3-W2b, M4.4-W3a | Q5 on frozen phrase fixtures; pizzicato and arco phrase results reported separately |
| M-DB-QP1 | Double bass: passive chain + speaks (QP, Q1) | M | M8-W3c4, M-VNP-QP1, M-VN-QP1, M5.5-W1 | QP vs same-specimen pair where it exists (else flagged population data); Q1: intended regime, note within ±10 c after the setup protocol, thresholds and onset within machine-player ranges, threshold/hysteresis map |
| M-DB-Q2 | Double bass recognisable | M | M-DB-QP1, M7.1-W2, M5.5-W1, M-VN-Q2, M8-W3c4, M2.5-W5a, M2.5-W5c, M-VNP-Q2, M2.5-W5b | Q2 on held-out notes, pizzicato and arco reported separately |
| M-DB-Q4 | Double bass inside the ruler + blind tests (Q3, Q4) | M | M-DB-Q2, M4.3-W1, M4.5-W1, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M8-W3c4, M6.4-W4 | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete; pizzicato and arco each pass on their own held-out set and are reported separately |
| M-DB-Q5 | Double bass phrases (Q5) | M | M-DB-Q4, M7.1-W4, M5.3-W2b, M4.4-W3a | Q5 on frozen phrase fixtures; pizzicato and arco phrase results reported separately |
| X-Rc1 | Recorder thresholds | S | M3.3-W2b, M8-W5 | §6.8 |
| M-RC-QP1 | Recorder/organ pipe: passive chain + speaks (QP, Q1) | L | M3.2-W4, M3.3-W1, M3.3-W2a, M3.3-W2b, M3.3-W2c, M8-W3h, M5.5-W1, X-Rc1 | Q1 vs Auvray 2012, Terrien 2014, Verge |
| M-RC-Q2 | Recorder/organ pipe recognisable (blocked until a licensed reference exists) | M | M-RC-QP1, M7.1-W2, M5.5-W1, M8-W6b | Q2 on held-out notes |
| M-RC-Q4 | Recorder/organ pipe inside the ruler + blind tests (Q3, Q4) | M | M-RC-Q2, M3.3-W3, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M8-W3h, M6.4-W2, M6.4-W3, M6.4-W4 | Q3 and binding Q4 |
| M-RC-Q5 | Recorder/organ pipe phrases (Q5) | M | M-RC-Q4, M7.1-W4, M5.3-W2b, M2.2-W3d | Q5 |
| X-Fl1 | Boehm flute register | S | M8-W3g, M6.2-W1, M8-W4 | §6.8 |
| M-FL-QP1 | Flute: passive chain + speaks (QP, Q1) | M | M-RC-QP1, M8-W3g, X-Fl1, M6.4-W3, M5.5-W1, M8-W5 | QP vs same-specimen pair where it exists (else flagged population data); Q1: intended regime, note within ±10 c after the setup protocol, thresholds and onset within machine-player ranges, threshold/hysteresis map |
| M-FL-Q2 | Flute recognisable | M | M-FL-QP1, M7.1-W2, M5.5-W1, M8-W3g, M2.5-W1 | Q2 on held-out notes |
| M-FL-Q4 | Flute inside the ruler + blind tests (Q3, Q4) | M | M-FL-Q2, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M3.3-W1, M3.3-W3, M8-W3g, M6.4-W2, M6.4-W3, M6.4-W4 | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete |
| M-FL-Q5 | Flute phrases (Q5) | M | M-FL-Q4, M7.1-W4, M5.3-W2b, M2.2-W3d, M2.1-W4b | Q5 on frozen phrase fixtures |
| M-VO-QP1 | Voice: passive chain + speaks (QP, Q1) | M | M3.1-W2, M6.4-W2, M8-W3k, M5.5-W1, M8-W5 | QP vs same-specimen pair where it exists (else flagged population data); Q1: intended regime, note within ±10 c after the setup protocol, thresholds and onset within machine-player ranges, threshold/hysteresis map |
| M-VO-Q2 | Voice recognisable | M | M-VO-QP1, M7.1-W2, M5.5-W1, M3.1-W2, M3.1-W6, M3.1-W8, M0.3-W3, M2.1-W8, M2.2-W3d, M8-W3k, M8-W6b | Q2 on held-out notes |
| M-VO-Q4 | Voice inside the ruler + blind tests (Q3, Q4) | M | M-VO-Q2, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M3.1-W8, M8-W3k, M6.4-W1, M6.4-W2, M6.4-W3 | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete |
| M-VO-Q5 | Voice phrases (Q5) | M | M-VO-Q4, M7.1-W4, M5.3-W2b, M2.2-W3d | Q5 on frozen phrase fixtures |
| X-Dk1 | Duduk same-specimen impedance (QP) | S | M8-W3l, M8-W8d, M2.1-W5 | as X-Cl-R0 |
| X-Dk2 | Low reed resonance from the reed continuum | S | M3.1-W5, M3.1-W1d, M6.4-W3, M8-W3l | falsifier §6.10 |
| X-Dk3 | Playing frequency/threshold vs pressure and lip force | S | X-Dk2, M8-W8d | falsifier §6.10 |
| M-DK-QP1 | Duduk: passive chain + speaks (QP, Q1) | M | X-Dk1, X-Dk2, M8-W3l, M-DR-QP1, M5.5-W1, M8-W5, X-Dk3 | QP vs same-specimen pair where it exists (else flagged population data); Q1: intended regime, note within ±10 c after the setup protocol, thresholds and onset within machine-player ranges, threshold/hysteresis map |
| M-DK-Q2 | Duduk recognisable | M | M-DK-QP1, M7.1-W2, M5.5-W1, M-DR-Q2, M3.1-W6, M6.4-W3, M8-W3l, M8-W6c, M2.4-W2, X-Dk3 | Q2 on held-out notes |
| M-DK-Q4 | Duduk inside the ruler + blind tests (Q3, Q4) | M | M-DK-Q2, M6.4-W2, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M3.1-W8, M8-W3l, M6.4-W3, M6.4-W4 | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete |
| M-DK-Q5 | Duduk phrases (Q5) | M | M-DK-Q4, M7.1-W4, M5.3-W2b, M2.2-W3d, M2.1-W4b | Q5 on frozen phrase fixtures |
| X-Tw1 | Tin whistle same-specimen impedance (QP) | S | M8-W3m, M8-W8c, M2.1-W5 | as X-Cl-R0 |
| X-Tw2 | Whistle thresholds and regime map | S | M-RC-QP1, M8-W3m, M8-W8c, M8-W7 | falsifier §6.11 |
| X-Tw3 | Windway/labium geometry sensitivity | S | X-Tw2 | §6.11 |
| M-TW-QP1 | Irish tin whistle: passive chain + speaks (QP, Q1) | M | X-Tw1, X-Tw2, M-RC-QP1, M8-W3m, M5.5-W1, M8-W5 | QP vs same-specimen pair where it exists (else flagged population data); Q1: intended regime, note within ±10 c after the setup protocol, thresholds and onset within machine-player ranges, threshold/hysteresis map |
| M-TW-Q2 | Irish tin whistle recognisable | M | M-TW-QP1, M7.1-W2, M5.5-W1, M2.1-W4a, M6.4-W4, M8-W3m, M8-W6c, X-Tw3 | Q2 on held-out notes |
| M-TW-Q4 | Irish tin whistle inside the ruler + blind tests (Q3, Q4) | M | M-TW-Q2, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M3.3-W1, M3.3-W3, M8-W3m, M6.4-W2, M6.4-W3, M6.4-W4 | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete |
| M-TW-Q5 | Irish tin whistle phrases (Q5) | M | M-TW-Q4, M7.1-W4, M5.3-W2b, M2.2-W3d, M6.4-W4 | Q5 on frozen phrase fixtures |
| X-Pp1 | Even harmonics vs Fletcher 2005 | S | M3.3-W1, M3.3-W3, M8-W3n, M6.4-W3, M3.2-W4 | falsifier §6.12 |
| X-Pp2 | Per-tube impedance and pitch (QP) | S | M8-W3n, M2.1-W5 | §6.12 |
| X-Pp3 | Lip-coverage pitch bend | S | X-Pp1 | §6.12 |
| M-PP-QP1 | Andean pan pipes: passive chain + speaks (QP, Q1) | M | X-Pp1, X-Pp2, M8-W3n, M3.2-W4, M5.5-W1, M8-W5 | QP vs same-specimen pair where it exists (else flagged population data); Q1: intended regime, note within ±10 c after the setup protocol, thresholds and onset within machine-player ranges, threshold/hysteresis map |
| M-PP-Q2 | Andean pan pipes recognisable | M | M-PP-QP1, M7.1-W2, M5.5-W1, M6.4-W3, M8-W3n, M8-W6c, X-Pp3 | Q2 on held-out notes |
| M-PP-Q4 | Andean pan pipes inside the ruler + blind tests (Q3, Q4) | M | M-PP-Q2, M2.6-W1, M2.6-W2, M2.6-W3, M6.2-W2, M7.1-W3, M7.5-W1, M6.3-W1, M6.3-W4, M5.5-W2, M3.3-W1, M3.3-W3, M8-W3n, M6.4-W2, M6.4-W3 | Q3 and binding Q4 (listening + spectrogram equivalence) on frozen held-out renders; §1.4 description complete |
| M-PP-Q5 | Andean pan pipes phrases (Q5) | M | M-PP-Q4, M7.1-W4, M5.3-W2b, M6.2-W1, M2.2-W3d | Q5 on frozen phrase fixtures |
| BR-03 | Jawari pack: measured data or refusal | S | M8-W1 | calibration removed |
| BR-04 | Authored-constant register | S | — | every BR-04 constant tracked to its replacing item |
| DEF-09 | Helmholtz neck end correction and loss | M | M2.1-W5 | each end's correction from its declared termination via the same compact-radiation owner as R; label consistent; neck viscothermal loss from fs-duct Bessel on non-fixture paths; resonance vs analytic within 1 %; `cavity_modal` consumers run natively |
| L-1 | Site section #instruments | M | OD1, M-CL-Q4 | §9.2 edits only (index.html section and nav, one main.ts mount, `tools/instrgrade/site_data.py`); `npm run build` passes; Vercel and frankensim.org/music checked in headless Chrome; maintainer ack in Agent Mail |
| L-2 | engine-perform crate | L | M5.4-W2, M1.1-W1d, M4.4-W3e | bit-identical same-ISA; same gates or render-ahead (OD2); absorbs PianoAudio (1/2/4× through `pcm_wav::decimate`) and the `withOversampling`/`options.substeps` API; retires all engine-piano `#[path]`; bit-identical per factor; fixed rung, no governor; render-ahead `skippedFrames = 0`; fs-live-fma reused; `public/live/perform` gitignored; the existing reduced live piano is not an acceptance engine (OD2) |
