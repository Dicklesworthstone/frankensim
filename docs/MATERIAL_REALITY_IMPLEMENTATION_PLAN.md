# Material Reality: Integrated Physical State and Material Substitution

Status: implementation plan; no new physical capability is claimed by this document.
Requested by the owner on 2026-09-04.
Source investigation began at `1be30297741e491357b75dac56874f05476f5d1d` on shared `main`.
The tree contained unrelated staged and unstaged work throughout the investigation.
This is a targeted execution delta to the existing architecture and domain plans.
It does not replace their ownership, scientific contracts, or deferred research gates.

## 1. The outcome

A user describes bodies, their composition and manufactured condition, geometry,
supports, interfaces, environment, and physical excitation.
FrankenSim evolves their physical state and derives motion, sound, and images.
Changing a material or an environmental condition changes that same calculation.
No instrument name, cinematic preset, or desired sound selects hidden physics.

The motivating examples are deliberately demanding:

- A disc may be copper, a specified stainless grade, gold, wood, or ruby.
- Its radius, thickness, edge profile, supports, and base geometry may change.
- Air pressure, temperature, humidity, and composition may change.
- A lead disc may heat, soften, yield, partially melt, flow, and resolidify.
- A string may be fretted by a compliant finger, plucked, bowed, or struck.
- Its supports, winding, diameter, material, initial stress, and temperature matter.
- A duct may be a flute, exhaust, horn, pipe, or machine cavity by its description.
- Its wall elasticity, thermoviscous gas losses, openings, and flow determine response.
- Electrical, magnetic, chemical, optical, and surface effects use the same state.

The success criterion is reusable physical composition, demonstrated in these cases.
It is not a new Euler solver, clarinet module, material-name switch, or audio preset.

### 1.1 What “100% based on physics” means operationally

Every simulated effect must originate in an explicit physical model and its inputs.
Animation follows solved configuration and deformation.
Sound follows physical excitation, structural/fluid response, propagation, and observation.
Appearance follows solved geometry and state-conditioned optical properties.
Temperature follows heat/energy transport, not elapsed-time appearance curves.
State transitions follow thermodynamics and, where needed, kinetics.

Constitutive physics requires measurements as well as theory.
Chemistry names alone do not determine roughness, grain structure, damping, or fatigue.
Continuum descriptions and numerical discretizations are approximations with domains.
“Physics-derived” therefore cannot honestly mean an exact prediction of all reality.
It means no invented causal effects, explicit closures, quantified error where supported,
and visible unknowns where the available data or model cannot support a prediction.

A numerical conservation check is not experimental validation.
A fitted material law is not an exact first-principles derivation.
A plausible sound is not proof of the pressure field.
A visually convincing melt is not proof of conserved enthalpy or mass.
All four can be useful only under their actual evidence boundaries.

### 1.2 Scope and delivery sequence

The first delivery removes actual material and environment shortcuts in existing paths.
The second connects complete material queries to ordinary solid-state substitutions.
The third adds spatial thermal/mechanical evolution and changes in model validity.
The fourth crosses solid/liquid boundaries with conserved state and real free surfaces.
The fifth broadens electromagnetic, chemical, porous, and degradation couplings.

The plan retains the entire ambition without putting every domain on the first release.
Advanced fluid, EM, and research work remains behind the existing domain gates.
No task claims readiness merely because its interface or data schema exists.
No task can close by adding only receipts when its outcome requires a physical solver.

## 2. Current architecture and verified source findings

These are source observations, not newly executed Cargo or experimental results.
README maturity prose is useful orientation but contains older boundaries in places.
Implementers must inspect the named source and live bead before claiming ownership.

### 2.1 Existing material spine

| Owner | Existing responsibility | Reuse and gap |
|---|---|---|
| L1 `fs-matdb` | Immutable claims, observations, cards, ordered interfaces, normalized packs, query receipts, covariance | Extend typed payloads and actual coverage; preserve one evaluator |
| L6 `fs-matdb-store` | FrankenSQLite derived index and canonical pack vault | Already exists; broaden pack families and compound discovery |
| L3 `fs-material::state_point` | Atomic requirement-driven property resolution | Reuse for every consumer; add missing typed bundles |
| L3 `fs-material::graph` | Law registry, declared roles, state schemas, single-pass DAG execution | Reuse local laws; closed physical feedback needs residual assembly above it |
| L3 `fs-material::phase` | Monotone specific-enthalpy solid/liquid equilibrium curve | Existing ingress only; no spatial heat or melt-flow solve |
| L3 `fs-conduction` | Steady conduction, theta-method transient conduction, duty schedules, lumped latent-heat transport | Extend spatial total enthalpy; do not build a duplicate thermal crate |
| L3 `fs-material::visco` | Fractional/Prony laws, fitting, work/storage/dissipation | Wire sourced laws instead of authored damping shortcuts |
| L1 `fs-thermochem` | Species and standard-state thermodynamic closure | Reuse reference conventions; do not duplicate gas/species identities |
| `fs-qty`, `fs-evidence`, `fs-blake3` | Dimensions, shared validity/evidence, identities | Extend current carriers rather than invent parallel ones |

`fs-matdb` has scalar and one-dimensional curve property values today.
Tensor frames, multidimensional responses, and more complete law transport need work.
Covariance is already preserved and queryable; do not recreate a statistics system.
`MaterialStateId` describes an immutable named reference state.
It must not become the mutable integration state of a heated or damaged specimen.

`MaterialStore` currently ingests `NormalizedPack`, not every pack family.
Its SQL index is a discovery aid; evaluations decode verified canonical pack bytes.
Its seal is an integrity check, not authentication.
Its current contract is single-writer and has no `Cx` integration.
The implementation plan expands that existing service rather than adding another database.

### 2.2 Concrete shortcuts and bounded seams

| Source | Finding | Required change |
|---|---|---|
| `fs-couple/src/thin_plate.rs::thermoelastic_for_density` | Chooses steel above a density threshold, aluminum otherwise; fixed reference temperature | Resolve actual expansion, heat capacity, conductivity, modulus, and state |
| `fs-couple/src/acoustic_realize.rs::gas_state` | Checks RH but constructs dry-air gas | Consume the existing moist-air constructor and carry the same gas state throughout |
| `fs-couple/src/acoustic_realize.rs::mode_zeta` | Authored frequency damping and Prony anchoring remain | Resolve physical loss mechanisms with explicit ownership |
| `fs-euler-disc-e2e/src/cinematic_fixture.rs` | Builds a one-temperature material card from caller scalars | Add canonical catalog/card bindings to the real runnable path |
| `cinematic_fixture.rs::build_fixture_production_mechanics` | Orthotropic contact and unequal material query points refuse | Add admitted anisotropic response and local-state interface evaluation |
| `cinematic_fixture.rs::prepare_fixture_disc_acoustics` | Current residual lane refuses admitted in-band structural modes | Wire actual modal plus residual radiation without double counting |
| `fs-euler-disc-e2e/src/structural_acoustics.rs` | Volume structural path admits only bounded cylinder profiles | Generalize through geometry-owned meshing and traces |
| `fs-material/src/phase.rs` and Euler thermal adapters | Enthalpy and volume changes exist; liquid states refuse | Add genuine thermal/deformable/free-surface coupling |

The Euler contract already describes same-card mechanics/contact/optics,
orthotropic structural modes, lumped heating, and uniform free expansion.
Those are valuable foundations and must not be reimplemented as new features.
The cinematic path's restrictions are still real product blockers.

### 2.3 Musical physics already has generic owners

The authoritative musical plan is
`COMPREHENSIVE_PLAN_FOR_OPTIMAL_MUSIC_ORIENTED_BUILDING_BLOCKS.md`.
Its geometry/material ingest law and generic-crate rule remain binding.
`fs-duct` owns duct/horn propagation and thermoviscous relations.
`fs-modal`, `fs-plate`, and `fs-nlmodal` own structural reduction families.
`fs-phs` owns the applicable passive/discrete-gradient machinery.
`fs-vfit` owns rational realizations with their actual approximation limits.
`fs-dcontact` and `fs-tribo` own collision and friction models.
`fs-bem` owns its current exterior radiation formulations.
`fs-aeroac` owns bounded aeroacoustic source models over fluid computations.
`fs-couple` composes these components and owns coupled runtime adapters.

An instrument fixture may name a clarinet, guitar, flute, or bass.
The constitutive or solver API may not dispatch on those names.
A complete object needs its bore/body geometry, material orientations, joints,
supports, exciters, and observation location; a material label is insufficient.

### 2.4 Existing commitments to reuse

| Existing bead | Current scope | Treatment |
|---|---|---|
| `frankensim-oecdy` | Material store, closed | Preserve; add successor tasks for documented gaps |
| `frankensim-0er85` | Common-material breadth, closed | Preserve source tranches; add functional coverage, not duplicate breadth |
| `frankensim-ext-matdb-seed-dataset-1sxe` | Curated source/process corpus | Existing ownership for further data acquisition |
| `frankensim-ext-thermal-domain-je8y` | Thermal storage, transport, latent heat, losses | Reuse charter; recognize existing conduction/lumped work |
| `frankensim-ext-solid-life-ladder-gahl` | Plasticity, damage, creep, fatigue, failure | Reuse for specific constitutive successors |
| `frankensim-ext-tribo-dry-baseline-tgbj` | Dry friction, thermal partition, wear | Reuse interface-law ownership |
| `frankensim-ext-tribo-lubrication-e1ob` | Film, mixed, EHL, cavitation | Reuse; no second thin-film solver |
| `frankensim-ext-active-materials-wccu` | Piezo, magnetostriction, thermoelectric | Retain staged domain ownership |
| `frankensim-euler-disc-emergent-flagship-t6314` | Scientific Euler campaign | Consumer and validation campaign, not generic physics owner |
| `frankensim-h7xu5` | Cinematic Euler production | Runnable consumer and presentation constraints |
| `frankensim-music-v8-root-3ez8g` | Musical building blocks | Reuse generic physics and existing listening gates |
| `frankensim-leapfrog-2026-program-i94v.3.2` | Material-lot passports/substitution | Preserve broader deferred scope; consume its eventual lineage |

New tasks must be residual implementation slices, linked to these owners.
Do not reopen closed beads merely because their scope was narrower than this request.
Do not remove another agent's assignment or reactivate owner-deferred portfolios.
Bead dependencies must target actual prerequisite deliveries, not whole research portfolios.

## 3. Idea-wizard: thirty candidates, then the best fifteen

The evaluation favors causal fidelity, reuse, and demonstrable user value.
Cost is assessed qualitatively; no unmeasured performance or schedule is promised.

| # | Candidate | Decision |
|---|---|---|
| 1 | Same material/environment state drives every consumer | Top five; central missing integration |
| 2 | Remove density/name-based physical shortcuts | Top five; immediate correctness gain |
| 3 | Conservative thermal, mechanical, and phase transition runtime | Top five; enables the lead example |
| 4 | Requirement-driven material and interface discovery | Top five; makes the catalog useful to simulations |
| 5 | Two reusable complete journeys: disc and excited elastic body/duct | Top five; proves composition |
| 6 | Tensor/multiaxis material responses and explicit frame transforms | Selected; wood, crystals, composites |
| 7 | Shared evolving surface state | Selected; joins friction, wetting, heat, and appearance |
| 8 | Physical damping and acoustic radiation ownership | Selected; prevents fabricated sound and duplicate losses |
| 9 | Manufacturing/lot/history and spatial heterogeneity | Selected; nominal chemistry is insufficient |
| 10 | Domain transitions with conservative state transfer | Selected; reduced models must end or escalate honestly |
| 11 | Source-backed complete material bundles | Selected; usable coverage rather than isolated facts |
| 12 | Electromagnetic/thermal/mechanical coupled materials | Selected; motors, sensors, pickups, hum |
| 13 | Chemical/moisture/degradation state | Selected; wood, corrosion, coatings, aging |
| 14 | Calibration and independent holdout measurements | Selected; converts causality into assessed prediction |
| 15 | Multirate execution and state-aware reduced models | Selected; practical cost under explicit error budgets |
| 16 | Material recommendations ranked by application QoI | Later consumer of 4/14/15; avoid a parallel recommender |
| 17 | Natural-language chemistry-name autocomplete | UI convenience later; never selects process/phase silently |
| 18 | A giant material struct with every property field | Reject; sparse data, different domains, and history make it brittle |
| 19 | Generate missing constants with an LLM | Reject; fabricates physical evidence |
| 20 | Atomistic simulation of every engineering material at runtime | Research horizon; cannot be the universal default |
| 21 | One mesh/resolution/time step for every phenomenon | Reject; incompatible scale and computational requirements |
| 22 | Instrument-specific calibrated oscillator presets | Reject as physical truth; fixtures may remain named objects |
| 23 | GPU backend before integrated profiling | Defer to existing accelerator doctrine and measured bottlenecks |
| 24 | Global proof of every coupled nonlinear solver | Preserve research ambitions in existing theorem portfolios |
| 25 | Learned material surrogate | Later, inside measured domain with fallback and uncertainty |
| 26 | Automatic property blending between unrelated materials | Reject; mixture/microstructure laws need explicit assumptions |
| 27 | One friction coefficient per bulk material | Reject; friction belongs to a stateful ordered interface |
| 28 | One “realism score” for the final video | Reject; different observables have different evidence |
| 29 | Comprehensive new provenance framework | Reject; existing receipts/ledger already own this |
| 30 | Universal ingestion crawler for every materials website | Defer; source-specific, license-aware tranches are more useful now |

### 3.1 Ranked top five and rationale

1. **One physical state across consumers.**
   This makes material substitution an actual simulation operation.
   It prevents mechanics, acoustics, and rendering from describing different bodies.
   Reuse current cards, query receipts, field identities, and coupling boundaries.
   The first acceptance test changes one material binding and checks all dependent outputs.

2. **Remove known physical shortcuts.**
   The density classifier and discarded moist-air state are concrete source findings.
   Fixing them gives immediate fidelity gains without a new framework.
   Tests must include a material whose density would select the wrong old model.
   The result is reusable in every plate, duct, and resonator consumer.

3. **Conserved thermal/mechanical/phase evolution.**
   Melting needs heat transport, latent heat, loss of shear support, fluid momentum,
   surface tension, contact/wetting, changing topology, and updated observations.
   The existing enthalpy primitive is the ingress, not the completed capability.
   Start with tractable Stefan and heated-member cases before the full disc.

4. **Material requirements and compound discovery.**
   A catalog of many isolated numbers cannot answer whether a simulation is supportable.
   Query the complete required property/model bundle at the intended conditions.
   Show missing axes and properties rather than substituting approximate chemistry.
   Reuse the existing FrankenSQLite index and immutable pack evaluator.

5. **Two complete, generic consumer journeys.**
   Disc and string/duct scenarios exercise different configurations of common primitives.
   Their joint acceptance prevents a beautiful single-purpose bolt-on.
   The same contact, material, thermal, transfer, and observation operators must participate.
   The scenarios remain parameterized data with deliberate unsupported outcomes.

### 3.2 Next ten and rationale

6. **Tensor and multiaxis responses:** scalar E cannot represent oriented wood or ruby.
7. **Surface state:** finish, contamination, oxide, and wetting couple several domains.
8. **Loss ownership:** audible radiation, mechanical damping, and heat must not debit twice.
9. **Manufactured specimens:** grade, grain, temper, winding, joints, and lot affect behavior.
10. **State transfer:** changes in resolution, topology, or model must conserve physical quantities.
11. **Complete bundles:** a small well-supported material set is more useful than a wide empty catalog.
12. **Active materials:** current, field, stress, and temperature jointly determine machine behavior.
13. **Chemical/moisture history:** wood swelling, corrosion, and charring are evolving physics.
14. **Independent experiments:** calibration must be separated from assessment of predictions.
15. **Multirate reduction:** audio-rate response and slow heating need different resolutions.

These fifteen ideas are implemented by the task catalog below.
The remaining candidates are either subsumed, rejected for physical reasons,
or retained in existing longer-term portfolios without blocking useful delivery.

## 4. Data, state, laws, and observations

### 4.1 Four distinct objects

**Material definition:** immutable source claims, named process/phase, model cards.
**Specimen:** geometry, composition fields, orientations, coatings, joints, residual state.
**Runtime state:** local enthalpy, composition, deformation, internal variables, surfaces.
**Observation:** microphone pressure, sensor response, radiance, displacement, or another QoI.

Do not collapse these into a monolithic “universal material” object.
Immutable card revision and mutable physical time have different meanings.
A runtime state cites the reference cards and model versions that govern its evolution.
A specimen can contain many material states and many interfaces.
An observation cannot silently inject energy or prescribe the dynamics it reports.

### 4.2 Property payloads and query semantics

Keep the existing `PropertyKey`, `PropertyClaim`, and `ValidityDomain` ownership.
Add quantity-kind distinctions where equal dimensions do not establish equal meaning.
Examples include energy versus torque and absolute temperature versus temperature difference.
Hardness must retain scale, indenter/load/dwell, and test context.
Acoustic attenuation must distinguish amplitude from intensity and nepers from decibels.
Optical extinction coefficient must not alias thermal conductivity.
Frequency must distinguish angular and cyclic conventions at adapters.
Strain/stress tensors must declare engineering versus tensor shear and basis ordering.

Multiaxis responses require complete support domains, not only independent bounding boxes.
An unmeasured corner of a temperature/frequency/moisture box is not automatically supported.
Interpolation stays inside an admitted source/model domain.
Extrapolation requires an explicit model, uncertainty, and caller policy; no implicit extension.
Bounded support, missing data, and uncertain correlation remain distinct query outcomes.

Full tensors carry basis/frame identifiers, symmetry class, and transformation receipts.
Rotate second- and fourth-order quantities with their appropriate transformation laws.
Do not rotate a scalar array as if it were a tensor without declaring its convention.
Anisotropic coupled blocks need the correct stability condition for the whole law.
Individual off-diagonal blocks need not be positive definite.

### 4.3 Complete material requirements

A consumer declares exactly the quantities and models needed for its claimed outputs.
Linear vibration can require density and elastic response without requiring yield stress.
Plastic contact adds strength/history requirements.
Thermoelastic damping adds expansion, heat capacity, and heat transport.
Melting adds latent heat/phase data, liquid properties, and interface laws.
Optical observation requires optical properties; it does not invent mechanical constants.

Discovery returns candidates with a supported-domain intersection and a gap explanation.
Evaluation rechecks the chosen exact cards and complete physical query points.
Loading a candidate is not admission of a complete coupled run.
Missing irrelevant data must not demote a claim that does not depend on it.
Missing relevant data must prevent the corresponding physical claim or trigger escalation.

### 4.4 Shared surfaces

Surface state includes material identity, orientation/lay, roughness spectrum and bandwidth,
coating thickness, oxide/contamination, wetness, wear, temperature, and relevant history.
The state can vary spatially and evolves under domain-owned laws.
Friction, contact conductance, electrical contact, wetting, and BSDF models read that state.
They remain distinct constitutive laws; one empirical conversion cannot replace them all.

There is no universal mapping from Ra to optical microfacet roughness or friction.
The roughness measurement scale and slope/spectral statistics matter.
A conversion requires a sourced/calibrated model with an applicability domain.
Contact angle belongs to a solid/liquid/gas system and can have advancing/receding history.
Restitution and modal damping are often system-dependent observations, not bulk constants.

## 5. Physical evolution and conservation

### 5.1 The coupling boundary

Local `ConstitutiveGraph` law evaluation is currently a single-pass DAG.
Closed feedback is not made correct by adding a cycle to that graph.
`fs-couple` assembles component residuals, interface transfers, and iteration schedules.
Strongly coupled blocks may solve monolithically through public residual/Jacobian APIs.
Partitioned blocks need convergence criteria and coupling-error accounting.
No lower-layer solver may read private state from a peer domain.

Each attempted macro-step operates on proposed state.
Admission, nonlinear convergence, transfer, event localization, and conservation checks
precede committing physical state, source debits, and emitted observations together.
Cancellation, budget exhaustion, or failed model admission preserves the last accepted state.
Retry must not debit a source twice or duplicate an emitted sound block.

### 5.2 Energy and entropy

Track kinetic, elastic/internal-memory, thermal, phase, field, surface, and potential energy
only where represented by the active model and with explicit reference conventions.
Internal transfers cancel when the accounting boundary encloses both participants.
External work, heat, radiation, mass/species transport, and control input remain explicit.
Numerical dissipation has its own account; it is not material damping or physical heating.

Frictional work has one debit and named thermal/other physical destinations.
Radiated acoustic energy must be distinguished from an uncoupled diagnostic microphone.
If a one-way radiation approximation is used, its missing back-reaction must be assessed.
Do not add modal radiation damping and separately debit the same outgoing field twice.
Thermoelastic damping must not be added when resolved thermomechanical transport includes it.
Latent heat is already inside enthalpy; do not add it again as an independent heat capacity.

Thermal entropy uses actual boundary temperatures and the chosen thermodynamic chart.
`temperature * heat_flow` is not a power-conjugate port.
Temperature pairs with entropy flow; ordinary heat flow already has units of power.
Open-stream enthalpy already includes pressure-flow work in its declared convention.
Chemical, electrical, and mechanical splits follow the existing extension charter.
A lossless Dirac interconnection alone does not establish discrete coupled passivity.

### 5.3 Heating before melting

Ambient temperature is a reservoir condition, not an instantaneous body temperature.
Conduction, convection, radiation, and contact determine local heating rates.
Density and thermal expansion cannot independently change the mass of a closed specimen.
Integrate expansion over the temperature path using the reference convention.
For nonuniform heating, solve strain compatibility and stress, not global scaling.
Clamped and free members respond differently.
Temperature changes strength, creep, damping, resistivity, magnetization, and optical response.

A lumped thermal model needs a justified internal-gradient approximation.
A low Biot-number check is useful only with its characteristic length and transfer regime.
When that approximation fails, route to a spatial thermal solve or return a named limit.
A small time step cannot fix a physically invalid lumped approximation.

### 5.4 Solid, mixture, and liquid

Use specific enthalpy as primary thermal state through isothermal latent-heat plateaus.
Retain composition/pressure/reference conventions with the equilibrium phase relation.
An alloy may have a solidus/liquidus interval and segregation behavior.
A polymer may soften through glass transition or decompose.
Wood must not be assigned an invented metal-like melting plateau.
Crystals and metals may have solid-state transitions before fusion.

The equilibrium enthalpy primitive does not represent nucleation, supercooling,
segregation, metastability, or reaction kinetics by itself.
Those need distinct admitted model cards and runtime internal variables.
Mass fraction and volume fraction are different when phase densities differ.
Convert with the constituent densities and preserve the corresponding mixture volume.

For melting with flow, solve thermal transport, momentum, phase support, and free surface.
Use the existing geometry/field representations and conservative transfer operators.
An enthalpy/porosity fixed-grid lane is a useful first spatial rung for bounded regimes.
Its mushy drag is a model closure, not a universal derivation from liquid fraction.
Free-surface motion additionally needs surface tension, wetting, and stress balance.
Density jumps require compatible volume flux/divergence; constant-volume incompressibility
must not silently override a phase-change volume difference.

The loss of load-bearing solid can release elastic energy and change contact topology.
Conserve or explicitly account for that release during the model transition.
A broken string's segments retain momentum; they do not vanish when a solver switches.
Droplets and disconnected components preserve mass/species and lineage.
Resolidification needs stress-free reference/locking conventions and shrinkage behavior.

### 5.5 Acoustics and rendering after state changes

Recompute or update structural response when geometry, pre-stress, or material changes.
Reduced bases need a validated update/transfer method, including near mode crossings.
Do not insert an arbitrary crossfade as a physical transition.
Sound generation uses the actual force/velocity or pressure/volume-flow histories.
A fracture or splash needs a resolved or explicitly modeled source mechanism.
It cannot be supplied by a prerecorded effect.

All observations use one accepted physical clock and state sequence.
Rendering may interpolate poses for exposure integration with a documented approximation.
Audio resampling needs explicit bandwidth, antialiasing, and latency treatment.
Microphone distance/direction, acoustic boundaries, and environmental propagation matter.
Camera exposure, spectrum, polarization support, and display conversion are observations.
Display tone mapping must not feed back into thermal emissivity.

## 6. Material coverage strategy

The existing `docs/MATERIAL_PROPERTY_TAXONOMY.md` remains the vocabulary reference.
Its broad categories are retained, with the semantic corrections in this plan.
No duplicate taxonomy or independent materials database will be introduced.

### 6.1 Coverage means complete use cases

Start from material/use/regime bundles rather than a target number of materials.
For each bundle report present measured data, admitted derived values, and exact gaps.
Room-temperature structural data do not imply a continuous high-temperature law.
An optical spectrum from a polished specimen does not describe an oxidized casting.
One steel's hardness or strength cannot fill another grade's missing property.

| Family | Initial named examples | Required emphasis |
|---|---|---|
| Common metals | OFHC copper, 304/316L/430/17-4PH stainless conditions, brass grade, aluminum condition | Temperature curves, process, elastic/plastic behavior, losses, optics |
| Phase demonstration | Lead of named purity/process, then solder/alloy | Solid/liquid enthalpy, density, viscosity, surface tension, oxide/wetting |
| Precious metals | Gold purity/alloy, silver where sourced | Density/elastic/plastic/thermal/optical consistency |
| Woods | Ebony species, balsa species, spruce/maple reference specimens | L/R/T frames, moisture, density basis, creep, damping, charring |
| Crystals/ceramics | Ruby composition/orientation, glass and alumina references | Anisotropic mechanics/optics, brittle failure, absorption |
| Soft/contact media | Finger tissue, cane, wool felt, elastomer | Protocol-dependent nonlinear/viscoelastic/contact response |
| Electrical/magnetic | Copper winding, electrical steel condition, permanent magnet grade | Field/frequency/temperature/history and loss separation |
| Fluids/coatings | Humid gas, water, named lubricant, oxide/varnish | Transport, wetting, chemical compatibility and coupled boundaries |

The examples identify acquisition targets, not assertions of data completeness.
Existing packs must be reused and their condition-specific limitations retained.
No fixed number of tranches is a scientific acceptance criterion.

### 6.2 Source policy

Prefer independently inspectable primary measurements and government technical datasets.
Handbook values can support bounded engineering estimates with their actual provenance.
Preserve method, specimen, process, uncertainty, covariance, censoring, and unit basis.
Store typical, tolerance, confidence, credible, and design-allowable quantities distinctly.
Do not convert an unstated fit error into a statistical confidence interval.
Do not average incompatible specimens, test scales, moisture bases, or processing states.

Dataset discovery does not establish redistribution rights.
Record the actual source terms before embedding data in redistributable packs.
If a source is restricted, retain a citation/gap or an authorized local-only artifact.
Source acquisition tasks may close their research portion with an honest absence finding;
their consuming physical capability remains blocked until required data or an admitted model exists.

## 7. Generic consumer journeys

### 7.1 Parameterized disc on a deformable, supported surface

Input: geometry charts, material cards, spatial frames, interface preparation,
support/joint conditions, gravity, initial configuration/velocities, and environment.
Mass and inertia derive from geometry and density, with an explicit comparison policy.
Comparisons distinguish equal geometry, equal mass, equal initial angular speed,
and equal initial energy; changing one does not secretly enforce another.

Dynamics combines rigid/deformable response, finite-patch contact, slip/rolling,
air/gap effects inside their domains, base compliance, and thermal feedback.
Finite-patch laws consume actual material response rather than a bulk friction constant.
Sound is generated from disc/base/interface/fluid histories and transported to receivers.
Images follow those same shapes, temperatures, interfaces, and optical states.

Acceptance has four successive rungs:

1. Ordinary solid-state parameter/material swaps in supported regimes.
2. Heating with thermal expansion, altered contact/modes, and independent base temperature.
3. Softening/creep/plastic deformation with acoustic/visual feedback.
4. Partial melting, loss of support, free-surface evolution, and continued observations.

No rung inherits experimental validation from a lower numerical rung.
Spin time and rankings need independent measurements with uncertainty and holdouts.
Euler-specific scientific claims remain owned by the existing scientific campaign.

### 7.2 Excited elastic body/string and gas-filled duct

A fretted string is a pre-stressed elastic body with distributed mass and stiffness,
supports/bridge/body coupling, and unilateral contact with frets and compliant tissue.
The finger is a geometry/material/contact description with a force or motion input.
Fret pressure changes contact area, effective constraints, damping, and tension.
Pitch is an outcome or an explicitly posed initialization target with force accounting.
Changing material at fixed end displacement differs from retensioning to fixed pitch.

A wound string requires winding/core geometry or a declared homogenized model.
Pluck, bow, hammer, and airflow are physical excitations of reusable components.
Body radiation and bridge losses remain coupled and separately identifiable.
Thermal expansion, relaxation, damage, or melting updates the same state.

A duct instrument is gas geometry, elastic/rigid walls, openings, and nonlinear excitation.
Reed contact, jet-edge feedback, lip valves, and tone-hole networks retain generic owners.
Humidity changes the admitted gas state and propagation; it is not merely a detune knob.
Wall-material changes may be small in a rigid-wall limit and substantial outside it.
The simulation must predict the appropriate sensitivity rather than amplify every swap.

### 7.3 Lessons from the jazz project

The related repository is `/Users/jemanuel/projects/jazz_chord_progression_editor_html`.
Its score/control scheduling, playable interaction, synthesis comparisons, and listening
workflows are useful consumer requirements and experimental-design inspiration.
Its synthesis approximations are not evidence of bottom-up material fidelity.
Integration must keep score intent separate from excitation and physical response.
No sampled timbre or authored filter becomes a constitutive property through import.
Useful concrete references are `dsp/concert-grand/src/plucked_v2.rs`
(finite-duration pick contact, retained sympathetic state, passive bridge),
`scripts/physical-foundry-plate-modes.ts` (geometry/material-derived plate operators),
and `src/audio/physical-realization.ts` (immutable excitation scheduling).
Its `flute_v2.rs::tuned_sound_speed` and `RESIDUAL_PULL_CENTS` alter sound speed
by note to reach a target pitch; that mechanism is excluded from this physical lane.
Its `libm` and mutable-global runtime patterns are also not implementation templates.
The extraction map in `docs/PHS1_FRANKENSIM_EXTRACTION_MAP.md` describes an older
FrankenSim snapshot; its claims about missing acoustic/viscoelastic owners are stale.

## 8. Implementation task catalog

Task keys below are stable plan-local references; the conversion section records Bead IDs.
Every task includes direct tests of changed behavior and a concrete consumer.
Dependency edges use `task depends on prerequisite` orientation.
Parent-child edges express hierarchy and are not substitutes for blocking edges.
Existing broad owners are related scope anchors, not blanket prerequisites.

Shared rules copied into every implementation bead:

- Work on `main`; preserve unrelated work; reserve only owned paths when mail is healthy.
- Runtime remains safe Rust with the existing Franken-only dependency direction.
- Use existing contracts, quantity/evidence/state carriers, and test seams.
- Source/current-static evidence does not establish executed or experimental proof.
- Poll `Cx` at bounded work boundaries; fail transactionally before publishing partial state.
- No arbitrary material-name inference, substituted data, or output-driven physical tuning.
- Focused G0/G1/G3 tests precede any larger declared G2/G4/G5 campaign.
- Use RCH for narrow heavy Cargo probes and DSR for a coherent repo-level gate.
- No new validator framework, dashboard, or ceremony is an acceptance substitute.

### MR01 — Use moist-air state in acoustic assembly

Group: A — Coherent state and immediate fidelity.
Priority: 1. Ambition: [S].
Depends on: none.
Existing anchor: `frankensim-music-v8-root-3ez8g.3.5` (closed kernel; reuse).
Files: `fs-couple/src/acoustic_realize.rs`, existing assembly tests and contract.
Outcome: humidity affects the same density, sound speed, impedance, and transport
used by excitation, ducts, and outgoing acoustic propagation.
Replace the dry-air constructor at the assembly seam with the existing admitted
`GasState::try_new_moist_air` path; do not reimplement the mixture equations.
Keep each constructor's temperature and vapor-fraction restrictions visible.
Do not add Stokes–Kirchhoff absorption on top of ISO absorption if already included.
Test the dry limit, a supported humid point, invalid RH, and unsupported temperature.
Check derived duct resonance/impedance against independently computed mixture quantities.
A waveform difference from only the final attenuation filter does not satisfy this task.
G0/G3: preserve deterministic repeats and the declared dry-path compatibility boundary.
Acceptance: the real assembly path consumes the moist state with no hidden dry substitution.

### MR02 — Replace density-guessed thermoelastic plate damping

Group: A — Coherent state and immediate fidelity.
Priority: 1. Ambition: [S].
Depends on: none.
Existing anchor: `frankensim-fsim-visco-damping-ybc75` (closed law implementation).
Files: `fs-couple/src/thin_plate.rs`, `fs-material/src/state_point.rs`, plate tests.
Outcome: copper, lead, wood, steel, and other plates never acquire aluminum/steel
thermal properties solely because their density crosses a threshold.
Replace `thermoelastic_for_density` with explicit admitted physical inputs.
Resolve expansion, heat capacity, conductivity, elastic response, and temperature.
Use an existing isotropic thermoelastic law only inside its physical assumptions.
For orthotropic material, require an admitted anisotropic law or report unavailable loss.
Do not substitute the isotropic law just because a scalar modulus is available.
Test two specimens with equal density but different thermal properties.
Test an old-threshold-crossing pair and temperature-dependent response.
G1: compare the admitted beam/plate limit to the independent Zener expression.
Acceptance: both linear and nonlinear consumers have no material-guessing branch.

### MR03 — Resolve physical loss models for vibration

Group: A — Coherent state and immediate fidelity.
Priority: 1. Ambition: [S/F].
Depends on: MR02.
Existing anchor: `frankensim-fsim-visco-damping-ybc75`.
Files: `fs-material/src/visco.rs`, `fs-couple/src/acoustic_realize.rs`, modal adapters.
Outcome: material-dependent decay follows an admitted complex modulus or loss law.
Replace authored `2e-7 * omega` bending loss in the predictive assembly path.
Retain explicit reduced-model inputs only as declared research-model inputs.
Use source-backed temperature/frequency/moisture/amplitude applicability.
Preserve the current distinction between measured in-band fitting error and a proof bound.
Identify material, air, thermoelastic, joint, and radiation losses separately.
Do not infer a full spectrum from one damping ratio without declaring the fitted model.
G1: compare an admitted Prony response to direct convolution and measured ring-down.
G3: change only the loss law and distinguish decay from stiffness-induced frequency shifts.
Test out-of-band query, missing source parameters, and unsupported amplitude.
Acceptance: a real string/plate assembly exposes each active physical loss source.

### MR04 — Compile coherent specimens from geometry and material cards

Group: A — Coherent state and immediate fidelity.
Priority: 1. Ambition: [S/F].
Depends on: none.
Existing anchor: `frankensim-euler-disc-emergent-flagship-t6314.2.11`.
Files: `fs-scenario/src/acoustic.rs`, material resolvers, existing L6 specimen adapters.
Outcome: geometry, density, elastic response, and manufactured condition produce
consistent mass, inertia, EA, EI, plate operators, and optical requirement bindings.
Reference existing cards by exact identity and resolve only consumer-needed properties.
Keep tension/pre-stress and support constraints in the mechanical state.
Distinguish fixed geometry, fixed mass, fixed tension, fixed displacement, and tuned pitch.
Expose explicit research scalar overrides as authored inputs with their own identity.
Never quietly merge those overrides into a source-backed card's claim authority.
Support regional material assignments and declared orientation fields at the seam.
G1: cylinder mass/inertia and uniform string rho*A/EA/EI against independent formulas.
G3: unit rescaling and one-binding material replacement through the same compiler.
Acceptance: one disc and one string consume the same requirement-driven material API.

### MR05 — Commit coupled physical state atomically

Group: A — Coherent state and immediate fidelity.
Priority: 1. Ambition: [F].
Depends on: MR04.
Existing anchor: `frankensim-ext-couple-cosim-lanes-pelj`.
Files: `fs-couple` runtime, existing state/checkpoint adapters, Euler consumer adapters.
Outcome: temperature, material laws, geometry, contact, modes, and observations
cannot commit different versions of the same physical instant.
Reuse existing transactional Euler stepping and common solver-state conventions.
Separate immutable card revisions from mutable per-run quadrature/region history.
Prepare candidate state, resolve laws, solve/transfer, assess, then publish together.
Rejected steps preserve last accepted state, source accounts, and emitted-block position.
Bind local material states separately; a hot body and cool base are valid inputs.
G4: inject cancellation/refusal after each existing transaction boundary.
Check bitwise unchanged accepted state and no duplicate source debit on retry.
G5: resume an accepted prefix and reproduce the suffix under one admitted profile.
Acceptance: two coupled consumers share the same committed state epochs.

### MR06 — Close physical loss and transfer accounting

Group: A — Coherent state and immediate fidelity.
Priority: 1. Ambition: [S/F].
Depends on: MR05.
Existing anchor: `frankensim-ext-couple-cosim-lanes-pelj`.
Files: existing coupling energy accounts and domain source adapters.
Outcome: friction, viscoelasticity, radiation, Joule heat, and phase energy
have explicit ownership and cannot be silently debited or credited twice.
Use the extension charter's energy/entropy conventions and current receipt carriers.
Separate physical dissipation from numerical stabilization and observation-only radiation.
Track stored internal-variable energy when a model includes it.
Do not enforce conservation by an arbitrary post-step energy renormalization.
Use actual boundary work and temperatures for source/reservoir exchange.
G1: two-body friction heat partition, damped oscillator, and driven thermal reservoir.
G3: enclosing both transfer participants cancels their internal exchange.
Mutation tests duplicate a thermal credit and omit a mechanical debit; both must fail.
Acceptance: a coupled mechanical/thermal step reports a residual with a justified scale.

### MR07 — Solve closed material feedback through public residuals

Group: A — Coherent state and immediate fidelity.
Priority: 2. Ambition: [F].
Depends on: MR05, MR06.
Existing anchor: `frankensim-ext-couple-cosim-lanes-pelj`.
Files: `fs-couple`, `fs-material/src/graph.rs` adapters, existing `fs-solver` interfaces.
Outcome: temperature-dependent response and interface feedback converge as a
coupled problem instead of an unacknowledged one-pass lag.
Keep the local ConstitutiveGraph DAG and law registry; add residual/Jacobian adapters above.
Specify unknowns, time levels, state ownership, stopping norms, and nonlinear tolerances.
Use partitioned iteration where justified and monolithic escalation where required.
Reject unsupported feedback cycles rather than silently choosing an evaluation order.
Consistent tangents differentiate the same trial update at the same committed history.
G1: coupled two-variable thermal/resistive or mechanical/thermal reference problem.
G3: compare converged partitioned and monolithic solutions as tolerances shrink.
G4: failed iteration rolls back through MR05 without state-history advancement.
Acceptance: an actual feedback consumer exercises the public residual path.

### MR08 — Invalidate and rebuild state-dependent operators

Group: A — Coherent state and immediate fidelity.
Priority: 1. Ambition: [S/F].
Depends on: MR04, MR05.
Existing anchor: `frankensim-ext-adaptivity-remap-u7ay`.
Files: existing `fs-recompute`, modal/geometry adapters, coupling state identities.
Outcome: material, geometry, support, or environment changes cannot reuse stale operators.
Bind caches to material/model/state/geometry/boundary dependencies actually consumed.
Distinguish immutable topology from changing metric/operator values.
Reuse recompute machinery; do not add a second dependency graph engine.
An update may reuse a basis only under an admitted change/error criterion.
Topology changes require explicit remap/reassembly and invalidation of geometry witnesses.
G3: alter temperature, orientation, support, and optical-only data independently.
Assert required recomputation and justified unaffected-cache reuse.
G4: interrupted rebuild never publishes a half-updated operator set.
Acceptance: a material swap changes all dependent mechanics/audio/render projections.

### MR09 — Make property axes and quantity kinds unambiguous

Group: B — Material data and discovery.
Priority: 1. Ambition: [S].
Depends on: none.
Existing anchor: `frankensim-ext-matdb-core-5hmy` (closed core).
Files: `fs-qty`, `fs-evidence::ValidityDomain`, `fs-matdb`, pack compiler/decoder tests.
Outcome: queries cannot confuse units, absolute/difference quantities, or test conventions.
Extend the existing shared validity carrier; do not introduce a competing domain type.
Require typed axes and quantity kinds where dimensions alone are insufficient.
Cover angular/cyclic frequency, hardness scales, moisture basis, attenuation convention,
tensor shear convention, and optical versus thermal uses of the symbol k.
Version frozen pack changes and provide explicit migration/refusal behavior.
Preserve old canonical bytes under their old decoder; never reinterpret in place.
G0: dimension/kind mismatch, Celsius-offset versus interval-width, Hz/rad/s fixtures.
G3: source-unit normalization produces the same physical query and retained provenance.
Acceptance: representative material and interface queries fail on semantic axis mismatch.

### MR10 — Add tensor and multiaxis material response payloads

Group: B — Material data and discovery.
Priority: 2. Ambition: [F].
Depends on: MR09.
Existing anchor: `frankensim-ext-matdb-core-5hmy` residual payload scope.
Files: `fs-matdb` values/query/packs, compiler, lower-layer frame identifiers.
Outcome: anisotropic and temperature/frequency/moisture-dependent data are representable.
Carry tensor order, symmetry, basis, support domain, uncertainty, and observation identity.
L1 stores frame identifiers; L2/L3 performs actual spatial transformations.
Distinguish sampled support from the enclosing axis-aligned bounding box.
Use an explicit interpolation rule with derivative/uncertainty limitations.
Do not create a full Cartesian grid of invented observations from sparse measurements.
G0: rotations and reciprocal symmetry under independently evaluated tensor formulas.
G3: equivalent unit/basis representations yield consistent consumer responses.
Reject unsupported support holes, nonfinite entries, and incompatible tensor conventions.
Acceptance: an oriented wood/crystal consumer resolves a source-backed tensor response.

### MR11 — Resolve derived property relations consistently

Group: B — Material data and discovery.
Priority: 2. Ambition: [S/F].
Depends on: MR09.
Existing anchor: `frankensim-ext-constitutive-graph-kagp`.
Files: `fs-material/src/state_point.rs`, `fs-matdb` query receipts, existing joint statistics.
Outcome: derived density, moduli, sound speed, heat capacity, and transport quantities
retain their parent data and thermodynamic/mechanical assumptions.
Choose an independent parameter set per law; cross-check redundant measured quantities.
Do not silently enforce E/G/K/nu identities across incompatible anisotropic observations.
Thermal expansion and density must use a common reference and volume convention.
Preserve joint uncertainty when deriving quantities; unknown covariance remains unknown.
Keep algebraic consequence, empirical correlation, and measured property distinct.
G1: isotropic wave-speed/modulus identities and density under free expansion.
G3: selecting a conflicting observation changes the derivation receipt visibly.
Acceptance: derived quantities never masquerade as independent measured constants.

### MR12 — Store all existing material pack families

Group: B — Material data and discovery.
Priority: 1. Ambition: [S].
Depends on: none.
Existing anchor: `frankensim-oecdy` (closed first store).
Files: `fs-matdb-store/src/lib.rs`, contract/tests; existing normalized pack APIs.
Outcome: the FrankenSQLite catalog can retrieve material cards, interfaces, models,
species associations, and ordinary property packs through one typed vault.
Reuse each existing canonical decoder and exact family/version identity.
Maintain links by immutable identities rather than inferred pack-name equivalence.
SQL remains a derived index; actual evaluations use verified canonical bytes.
Ingest a linked set transactionally and retain the current single-writer boundary.
Never add database calls to L1 material kernels or inner physics loops.
G0: family confusion and wrong version/hash refuse.
G4: mid-ingest failure leaves no partially usable linked set.
G5: rebuild parity and direct-pack versus store evaluation/receipt equality.
Acceptance: one body and one ordered interface resolve from the expanded vault.

### MR13 — Discover complete material/model bundles

Group: B — Material data and discovery.
Priority: 1. Ambition: [S/F].
Depends on: MR04, MR09, MR12.
Existing anchor: `frankensim-oecdy` successor and `frankensim-ext-matdb-seed-dataset-1sxe`.
Files: `fs-matdb-store`, scenario admission, CLI material-selection surface.
Outcome: users can ask which named material conditions support a declared simulation.
Intersect the complete requirements over intended state points or state trajectories.
Include model availability, interface counter-material, and curve support in discovery.
Return complete, partial, and unavailable candidates with named missing requirements.
Range discovery over curves must inspect their admitted domain, not just first-knot values.
Candidate discovery never bypasses exact evaluation or source selection at admission.
G0: unknown property versus known-property empty result remain distinct.
G3: adding an irrelevant weak property does not demote an unrelated supported query.
Test a material with room-temperature data but missing high-temperature coverage.
Acceptance: a lead-heating request explains its actual missing physical inputs.

### MR14 — Complete a lead solid/liquid consumer data bundle

Group: B — Material data and discovery.
Priority: 1. Ambition: [S/F].
Depends on: MR09.
Existing anchor: `frankensim-ext-matdb-seed-dataset-1sxe`.
Files: existing lead source tranche, `data/matdb/seed-v1`, compiler and query tests.
Outcome: a named lead specimen can support a bounded heating/melting experiment.
Reuse existing lead density/melting/hardness data; do not pretend that set is complete.
Acquire condition-matched solid/liquid enthalpy, density, k, cp, expansion, strength/creep,
liquid viscosity, surface tension, relevant wetting, and optical data for chosen observers.
Separate directly measured inputs from fitted closures and validation observations.
Use exact source terms and uncertainty conventions; no fabricated missing rows.
G0/G3: compile sources, normalize units, and query all required operating states.
Run dimensional/thermodynamic consistency checks without upgrading measurements to exactness.
Report remaining source gaps explicitly; a data-refusal finding does not unblock melting.
Acceptance: the declared consumer bundle is complete or its dependent journey remains blocked.

### MR15 — Complete ordinary metal substitution bundles

Group: B — Material data and discovery.
Priority: 2. Ambition: [S/F].
Depends on: MR09, MR13.
Existing anchor: `frankensim-ext-matdb-seed-dataset-1sxe` and Euler `.2.7`.
Files: existing copper, brass, stainless, gold, and aluminum source tranches.
Outcome: supported ordinary material swaps carry complete mechanics/loss/thermal/optical inputs.
Name alloy grade, temper, process, specimen condition, and source-supported temperature range.
Keep stainless families separate, including their magnetic and corrosion differences.
Acquire source-backed curves only over the chosen consumer domains.
Do not convert typical strength into a design allowable or fill gold from copper.
Use MR13's gap report to select the smallest useful next acquisition tranche.
G0: source-condition mismatch and unsupported axis queries refuse.
G3: compile/store/direct-query agreement for every delivered bundle.
Acceptance: at least two materially distinct metal configurations run the same solid-state path;
additional named materials retain explicit coverage gaps until their bundles are delivered.

### MR16 — Complete oriented wood, crystal, and soft-contact bundles

Group: B — Material data and discovery.
Priority: 2. Ambition: [F].
Depends on: MR09, MR10, MR13.
Existing anchor: seed dataset owner; closed music cane/tissue packs are reusable inputs.
Files: existing wood/cane/tissue packs plus bounded crystal/ebony acquisition tranches.
Outcome: wood/ruby/tissue substitutions use actual directional and protocol-dependent data.
Name botanical species, moisture basis, L/R/T axes, density convention, and specimen zone.
Ruby needs composition and crystallographic orientation; generic glass is not a substitute.
Finger tissue requires its contact/time/rate protocol, not a single universal hardness.
Preserve documented ebony and other source dead ends until new compatible evidence exists.
Separate small-strain elastic, viscoelastic, indentation, and nonlinear contact measurements.
G0/G3: orientation/units/moisture transformations and unsupported property refusal.
Validate cross-property consistency without averaging incompatible specimen populations.
Acceptance: one oriented solid and one compliant-contact consumer obtain complete admitted inputs.

### MR17 — Acquire stateful interface and surface data

Group: B — Material data and discovery.
Priority: 2. Ambition: [F].
Depends on: MR09, MR12.
Existing anchor: `frankensim-ext-tribo-dry-baseline-tgbj`, seed dataset owner.
Files: interface source packs, `fs-matdb` interface queries, consumer regressions.
Outcome: polish/lubrication/counter-material changes have sourced constitutive consequences.
Retain ordered surfaces, preparation, roughness bandwidth/lay, contamination, medium,
normal load/pressure, velocity, temperature, and named running-in/history condition.
Acquire friction, thermal/electrical contact, wear, adhesion, and wetting where relevant.
Do not infer restitution from a bulk name or hydrophobicity without named fluids/gas.
Keep unsupported surface conditions absent rather than supplying a default coefficient.
G3: reversing a directional interface changes its identity and only admitted symmetric laws commute.
Test two finishes with otherwise equal bulk states and source-backed response differences.
Acceptance: an existing contact consumer evaluates an ordered interface pack end to end.

### MR18 — Exercise corpus compilation through store and physical consumers

Group: B — Material data and discovery.
Priority: 2. Ambition: [S].
Depends on: MR12, MR13, MR14, MR15, MR16, MR17.
Existing anchor: `frankensim-oecdy` documented compiler/store E2E residual.
Files: existing xtask material-pack tests and one existing-style e2e script.
Outcome: acquired source data reaches real physical queries without hand-entered substitutions.
Compile selected source manifests using the actual pack compiler.
Ingest supported family sets into a disposable task-owned database.
Query through MR13 and execute representative material resolvers.
Compare exact receipts and values to direct compiled-pack evaluation.
Use per-stage source IDs, normalized units, selected claim, query point, and refusal in logs.
G4: interruption and a deliberately malformed linked pack leave the accepted corpus intact.
G5: deterministic rebuild from the same selected sources.
Acceptance: this is one corpus-to-consumer test, not a new general test framework.

### MR19 — Add spatial nonlinear total-enthalpy transport

Group: C — Thermal and mechanical evolution.
Priority: 1. Ambition: [F].
Depends on: MR09, MR11.
Existing anchor: `frankensim-ext-thermal-domain-je8y`.
Files: `fs-conduction/src/transient.rs`, `lumped.rs` patterns, `fs-material::phase`, thermal tests.
Outcome: spatial temperature/phase fields evolve with temperature-dependent material response.
Extend existing fs-conduction storage/transport; do not create a parallel thermal solver.
Use primary total enthalpy with explicit mass/volume/reference convention and nonlinear k.
Retain latent plateaus as an enthalpy relation; separately label any apparent-cp regularization.
Include energetic internal-variable changes or explicitly freeze them in the admitted model.
Specify boundary conditions, nonlinear residual/tangent, and accepted-step energy balance.
G1: manufactured single-phase nonlinear diffusion and one-dimensional Stefan front/energy.
Compare against the existing lumped solution in an independently justified uniform limit.
Test latent plateau, zero applied flux, invalid phase data, and timestep/mesh refinement.
Acceptance: a spatial front and heat balance are computed, not prescribed by a fixture.

### MR20 — Couple local thermal contact and boundary reservoirs

Group: C — Thermal and mechanical evolution.
Priority: 1. Ambition: [F].
Depends on: MR05, MR06, MR17, MR19.
Existing anchor: thermal-domain and dry-tribology owners; Euler `.3.13`.
Files: `fs-conduction` boundary/contact adapters, `fs-contact`, `fs-couple`.
Outcome: a hot disc/cool plate or hot string/cool support exchanges heat physically.
Evaluate each bulk at its own state and the interface at explicitly declared interface variables.
Use contact pressure, separation, finish, and medium in the admitted conductance law.
Account for convection/radiation reservoirs and friction heating exactly once.
Do not require the entire QueryPoint for both bodies to be equal.
Do not add a radiation/Robin source atop an existing occupied boundary without composition.
G1: two-body thermal equilibration and a prescribed flux/contact-resistance problem.
G3: surface relabeling preserves energy and temperature-exchange signs.
G4: failed interface query rolls back thermal and mechanical candidates together.
Acceptance: unequal initial body temperatures evolve through the real coupled interface.

### MR21 — Extend nonuniform thermoelastic deformation and pre-stress

Group: C — Thermal and mechanical evolution.
Priority: 1. Ambition: [F].
Depends on: MR04, MR19.
Existing anchor: thermal-domain owner and existing `fs-solid::linear3` thermal strain.
Files: `fs-solid` existing thermal eigenstrain/assembly, material adapters, coupling tests.
Outcome: temperature gradients and support constraints produce actual strain and stress.
Reuse existing 3-D small-strain thermal loading and integrated expansion conventions.
Advance spatial state under changing thermal loads and recalculate compatible deformation.
Mass, reference configuration, and density must be consistent with the deformation map.
Do not replace nonuniform expansion by global geometric scaling.
Represent fixed-end and force-controlled string/support boundaries distinctly.
G1: free expansion, constrained thermal stress, and a thermal-gradient bending fixture.
G3: rotate an orthotropic frame and geometry together; the physical result transforms coherently.
Test inversion/domain exit and return a named finite-deformation escalation requirement.
Acceptance: changed geometry/pre-stress is passed to contact and vibration consumers.

### MR22 — Add temperature/history-dependent creep and plastic softening

Group: C — Thermal and mechanical evolution.
Priority: 2. Ambition: [F].
Depends on: MR05, MR06, MR19, MR21.
Existing anchor: `frankensim-ext-solid-life-ladder-gahl`.
Files: `fs-material` law nodes and `fs-solid` existing nonlinear update paths.
Outcome: a heated loaded member relaxes, yields, and sags before complete melting.
Add one sourced bounded creep/viscoplastic law with consistent algorithmic tangent.
Retain plastic strain, hardening, creep, thermal state, and relevant stored energy.
Represent tension/compression asymmetry only when the selected law/data supports it.
Do not infer creep or yield evolution from melting point alone.
Use a finite-deformation formulation when small-strain kinematics cease to be valid.
G1: constant-stress creep and stress relaxation against analytical law solutions.
G0: objectivity, return-map admissibility, and nonnegative dissipation under declared assumptions.
Acceptance: an actual loaded string/rod fixture changes shape and tension through the law.

### MR23 — Carry heterogeneous manufactured state through physics

Group: C — Thermal and mechanical evolution.
Priority: 2. Ambition: [F].
Depends on: MR04, MR10, MR21.
Existing anchor: Euler `.2.11`, `frankensim-ext-asbuilt-manufacturing-r58b`.
Files: region assignment, `fs-solid`, `fs-modal`, existing specimen/geometry adapters.
Outcome: composites, coatings, wood grain, wound strings, joints, and residual stress matter.
Represent spatial composition/orientation and assembly interfaces using existing regions/traces.
Homogenized members require explicit scale separation and an error/validation domain.
Do not volume-average nonlinear failure or tensor properties without a stated model.
Support one layered or wound specimen and one rotated orthotropic specimen first.
Preserve manufacturing metadata without demanding the complete deferred lot-passport system.
G1: layered conduction/elastic response against analytical series/parallel or laminate limits.
G3: equivalent regional partitioning preserves integrals within discretization error.
Acceptance: shared region data affects structural, thermal, and acoustic response coherently.

### MR24 — Route thermal and mechanical models at validity boundaries

Group: C — Thermal and mechanical evolution.
Priority: 2. Ambition: [F].
Depends on: MR08, MR19, MR21, MR22.
Existing anchor: existing regime/fidelity and coupling owners; Euler `.3.13`.
Files: `fs-regime`/existing selector adapters, coupling admission, contracts/tests.
Outcome: low-cost models escalate or stop before their physics assumptions fail.
Reuse existing selection machinery with explicit model-domain predicates and budgets.
Check lumped versus spatial thermal validity, small versus finite strain, and solid support.
Use predictive event brackets where available; otherwise report bounded numerical localization.
A property-domain boundary may require refinement, a different law, or a named unavailable result.
Do not use threshold-triggered shape morphs or continue a rigid body through liquid fraction.
G3: tighter resolution/steps stabilize transition observations within stated error.
G4: insufficient budget retains the accepted prefix and exact missing next capability.
Acceptance: at least one actual heated-member run triggers a correct model transition/refusal.

### MR25 — Define and implement bounded solid/liquid mixture response

Group: D — Phase transitions and free surfaces.
Priority: 2. Ambition: [F].
Depends on: MR14, MR19, MR22.
Existing anchor: thermal-domain, solid-life, and thermochemical-core owners.
Files: `fs-material::phase`, local constitutive adapters, thermal/fluid coupling tests.
Outcome: partial melting changes thermal and mechanical response using an admitted mixture model.
Distinguish mass and volume fractions and retain phase densities and reference enthalpy.
Use a declared phase-support/permeability relation for a bounded material/regime.
Do not equate liquid fraction with a universal linear stiffness interpolation.
Represent pressure/composition dependence explicitly or restrict the initial case accordingly.
Keep equilibrium and kinetic phase laws separate; no hidden nucleation/supercooling claims.
G0: fraction conversion, endpoints, positivity/admissibility, and latent-heat consistency.
G1: a stationary mushy-region response with independent balance calculation.
Acceptance: mixed phase supplies actual transport/stress closure to a coupled solver.

### MR26 — Implement capillary and wetting boundary forces

Group: D — Phase transitions and free surfaces.
Priority: 2. Ambition: [F].
Depends on: MR17.
Existing anchor: `frankensim-ext-porous-capillary-biz5`.
Files: existing fluid/geometry boundary operators, interface law adapters, coupling tests.
Outcome: droplets and molten material respond to curvature, surface tension, and substrate wetting.
Use orientation-consistent curvature/normal conventions and resolved surface stress balance.
Wetting requires named liquid, solid, gas, temperature, and contact-line assumptions.
Treat contact-angle hysteresis and dynamic contact-line regularization as model closures.
If thermocapillary gradients are supported, include tangential surface-tension stresses.
Do not claim arbitrary contact-line physics from an imposed equilibrium angle.
G1: Laplace pressure, capillary wave/oscillation limit, and sessile-drop equilibrium.
Measure spurious currents and spatial convergence; a stable picture alone does not pass.
Acceptance: capillary forces enter momentum and surface energy accounting coherently.

### MR27 — Evolve a conservative incompressible free surface

Group: D — Phase transitions and free surfaces.
Priority: 2. Ambition: [F].
Depends on: MR06, MR26.
Existing anchor: moving-boundary fluids, existing `fs-flux`/`fs-lbm` domain owners.
Files: existing fluid solver/field/geometry machinery and bounded free-surface tests.
Outcome: liquid moves, slumps, separates, or coalesces as a computed fluid state.
Choose a concrete existing-solver extension with explicit velocity/pressure/interface coupling.
Implement conservation and free-surface traction; do not build another standalone fluid engine.
State the initial density/viscosity and low-Mach assumptions.
Phase volume sources are supplied explicitly by MR25 when the two phases differ in density.
Preserve owner-deferred compressible/turbulence scope; it is unnecessary for the first melt case.
G1/G2: hydrostatics, gravity-driven flow, capillary equilibrium, and a sourced dam-break/drop case.
G3: refine mass and momentum error and identify numerical surface diffusion.
Acceptance: geometry is reconstructed from solved phase/volume fields and conserves mass.

### MR28 — Transfer state conservatively across representations and topology

Group: D — Phase transitions and free surfaces.
Priority: 2. Ambition: [F].
Depends on: MR05, MR08, MR25.
Existing anchor: `frankensim-ext-adaptivity-remap-u7ay`.
Files: existing remap/field/mesh adapters and coupling state transport.
Outcome: remeshing or solid-to-fluid handoff preserves the modeled physical state.
Transfer mass, species, linear/angular momentum, total energy, and internal variables.
Specify reference frames, reference volumes, constraints, and treatment of newly created cells.
Retain unresolved transfer error and applicable positivity/admissibility conditions.
Track split/merge lineage for disconnected components using existing artifact conventions.
Never erase elastic energy or zero velocities when a body stops being rigid.
G1: analytic translating/rotating body remap and enthalpy-bearing split/merge integrals.
G3: representation round trips with independently calculated conserved moments.
G4: failed transfer commits neither geometry nor material history.
Acceptance: an actual solver handoff consumes the transfer, beyond standalone array tests.

### MR29 — Couple heated deformable solids to evolving melt

Group: D — Phase transitions and free surfaces.
Priority: 2. Ambition: [F].
Depends on: MR07, MR20, MR24, MR25, MR27, MR28.
Existing anchor: thermal-domain and coupled-runtime owners.
Files: `fs-couple`, existing thermal/solid/fluid public adapters, phase conformance cases.
Outcome: a heated loaded specimen continues through partial melting and liquid motion.
Solve thermal transport, deformation, phase support, contact, and fluid momentum together.
Update changing interfaces, stress release, mass/inertia, and boundary heat transfer.
Make density-change volume flux compatible with the fluid continuity equation.
Use a bounded pure-material low-Mach first case; broaden alloy kinetics later.
G1: Stefan heat balance coupled to a mechanically simple supported melt case.
G2: a sourced slumping or melt-front experiment with separate calibration inputs.
G4: cancellation/refusal preserves the complete accepted multiphysics prefix.
Acceptance: the result contains computed phase, shape, velocity, stress, and energy histories.

### MR30 — Add resolidification and non-melting transition behavior

Group: D — Phase transitions and free surfaces.
Priority: 3. Ambition: [F].
Depends on: MR25, MR28, MR29.
Existing anchor: thermal-domain, solid-life, and gas-chemistry owners.
Files: phase/history laws and existing solid/fluid/thermochemical adapters.
Outcome: cooling can restore load-bearing solid with a stated reference state;
wood/polymer regimes do not follow an invented metal melting rule.
Specify solid-locking strain, shrinkage, residual stress, and trapped-phase assumptions.
For an initial non-melting material, implement one sourced softening or decomposition law.
Keep decomposition mass/species/energy products explicit and preserve unmodeled products as gaps.
Do not implement combustion merely by changing optical color or deleting mass.
G1: freeze/melt energy cycle and a constrained shrinkage reference.
G0: phase/mass/species bookkeeping and invalid transition/model selection.
Acceptance: one real cooling case and one distinct non-metal transition exercise the laws.

### MR31 — Derive pre-stressed strings and coupled supports from specimens

Group: E — Generic sound-producing systems.
Priority: 1. Ambition: [S/F].
Depends on: MR03, MR04.
Existing anchor: closed music string/body program; generic acoustic owner `0ja4`.
Files: `fs-scenario::acoustic`, `fs-couple` string/bridge adapters, `fs-nlmodal`.
Outcome: string material/diameter/support changes determine dynamics without independent EA/EI guesses.
Derive distributed mass, stretching/bending stiffness, and physical damping from the specimen.
Solve the declared equilibrium/pre-stress under end displacement, force, or explicit tuning target.
Keep actual actuator work when retensioning; fixed pitch is not automatic material-swap behavior.
Support reciprocal bridge/body loading where the selected model claims coupling.
G1: ideal-string and stiff-string limits with independently calculated frequency/tension.
G3: change only diameter/material and check the corresponding causal parameter changes.
Test preserved unplayed-string/body state under repeated excitations.
Acceptance: one generic rod/string fixture supports pluck and another physical excitation.

### MR32 — Simulate compliant finger/fret and pick contact

Group: E — Generic sound-producing systems.
Priority: 2. Ambition: [F].
Depends on: MR16, MR17, MR31.
Existing anchor: generic contact/dcontact/tribology and music gesture owners.
Files: `fs-scenario::gesture`, `fs-dcontact`, `fs-couple` string/contact adapters.
Outcome: fretting/plucking follows contact geometry and material response rather than pitch edits.
Model finger/pick as bounded deformable or admitted reduced contact bodies.
Drive actual force/displacement histories and include actuator work and dissipation.
Resolve fret contact, changing contact area, slip, release, buzzing, and string tension.
Reuse generic unilateral/contact laws so the same interaction can model a cable or seal.
Borrow finite-duration contact and retained-state ideas from Jazz, not its authored timbre.
G1: compliant impact and quasi-static indentation/contact reference.
G3: finger force/stiffness sweeps with physical contact duration and frequency response.
Acceptance: a pluck while fretted creates solved contact forces and observer pressure.

### MR33 — Resolve regional gas and deformable duct-wall response

Group: E — Generic sound-producing systems.
Priority: 2. Ambition: [F].
Depends on: MR01, MR04, MR10, MR16.
Existing anchor: `frankensim-ext-acoustics-core-0ja4`; existing fs-duct owners.
Files: `fs-duct`, `fs-couple` duct/wall/exciter adapters, scenario region bindings.
Outcome: gas and wall material/environment changes reach the same duct calculation.
Admit separate excitation gas, interior gas, and exterior observation environment where supported.
Use real regional state for density/speed/viscosity/thermal losses and impedance.
Couple wall compliance when modeled; otherwise state the rigid-wall limit and its error question.
No note-dependent speed-of-sound tuning or material-name filter selection is allowed.
G1: closed/open duct limits, thermoviscous attenuation, and compliant-wall perturbation.
G3: humidity/temperature/pressure and wall-material sweeps with justified sensitivity.
Preserve the current unbranched bore-extraction limit and active owner of `frankensim-b2can`.
Acceptance: one geometry-defined duct runs without instrument-name dispatch.

### MR34 — Update reduced vibration states under changing physical conditions

Group: E — Generic sound-producing systems.
Priority: 2. Ambition: [F].
Depends on: MR08, MR21, MR22, MR31.
Existing anchor: modal/reduction and conservative-remap owners.
Files: `fs-modal`, `fs-couple` modal-time adapters, existing structural reduction tests.
Outcome: changing stiffness, pre-stress, geometry, or supports updates modes and state coherently.
Track subspaces through near-degenerate modes instead of relying on sorted mode numbers.
Project displacement/momentum with declared mass-inner-product and energy error accounting.
Use residual/tail estimates to decide whether a basis remains adequate for requested observables.
Escalate or rebuild when a topology/material transition invalidates the reduction.
Do not create energy by rescaling amplitudes to preserve a desired sound level.
G1: slowly varying oscillator/string against a direct time-domain reference.
G3: mode crossing, temperature cycle, and refined mode-budget comparison.
Acceptance: a heated string/plate changes sound through transferred physical vibration state.

### MR35 — Compose modal and broadband physical radiation

Group: E — Generic sound-producing systems.
Priority: 1. Ambition: [F].
Depends on: MR03, MR06.
Existing anchor: `frankensim-h7xu5.7.8` and generic acoustic owner `0ja4`.
Files: `fs-couple::modal_acoustic_time`, `broadband_radiation`, Euler structural adapters.
Outcome: a soft/thin disc can radiate in-band structural modes and complementary response.
Remove the cinematic selection restriction by composing actual supported source subspaces.
Use existing solver-neutral radiation tables; avoid a fs-couple→fs-bem dependency cycle.
Define the modal/residual split and check that no forcing or radiation component is counted twice.
Keep failed held-out radiation fits unavailable; do not hide them with a resonator preset.
Where two-way loading is absent, quantify or explicitly bound the approximation's domain.
G1: single-mode radiation and analytical low-frequency source limits.
G3: increasing modal coverage converges with complementary residual accounting.
Acceptance: a physically resonant disc no longer refuses solely for having in-band modes.

### MR36 — Schedule physical gestures and preserve shared resonant state

Group: E — Generic sound-producing systems.
Priority: 2. Ambition: [S/F].
Depends on: MR05, MR31, MR32, MR33.
Existing anchor: existing music gesture/render contracts.
Files: `fs-scenario/src/gesture.rs`, `fs-couple/src/render.rs`, common control adapters.
Outcome: declared fret/bow/pluck/support/valve gestures actually reach physical runtime controls.
Extend the current limited `ControlDelta` boundary with supported physical variables.
Compile score intent to immutable excitation events with exact physical/sample timing.
Maintain string/body/cavity state across notes, releases, rests, and cooperative blocks.
A performer controller may change pressure/embouchure/force, never medium constants for tuning.
Use bounded controls with energy/source accounting and explicit unsupported actions.
G3: block-size independence and identical event timing under cooperative rendering.
G4: cancellation/resume with sustained and sympathetic vibration.
Acceptance: the same generic graph executes a multi-gesture phrase without state resets.

### MR37 — Admit anisotropic and locally non-isothermal contact

Group: F — Interfaces and physical observations.
Priority: 2. Ambition: [F].
Depends on: MR10, MR16, MR17, MR20.
Existing anchor: `frankensim-b8bxd.5.1`, Euler contact/interface owners.
Files: `fs-contact::interface_binding`, normal response laws, Euler production adapter.
Outcome: wood/crystal contacts and hot-body/cool-base interactions can be modeled honestly.
Use the existing bound normal-contact model ingress instead of constructing private parameters.
Add a bounded anisotropic compliance/contact rung with explicit orientation and geometry assumptions.
Keep isotropic Hertz as its own validated limiting case, not an automatic fallback.
Resolve ordered surface and bulk states independently at the current contact region.
Test normal/tangential coupling and applicability outside small-strain half-space assumptions.
G1: anisotropic reference solution and isotropic degeneration.
G3: joint frame rotation and surface exchange under a declared symmetric law.
Acceptance: one anisotropic specimen passes the real production contact path with sourced inputs.

### MR38 — Generalize structural volume and boundary observation geometry

Group: F — Interfaces and physical observations.
Priority: 2. Ambition: [F].
Depends on: MR04, MR23.
Existing anchor: Euler geometry and `fs-mesh` owners.
Files: geometry/mesh producers and `fs-euler-disc-e2e/src/structural_acoustics.rs` adapters.
Outcome: annular/chamfered/tapered and other admitted solids share structural/acoustic geometry.
Extend existing chart-to-volume machinery rather than adding Euler-specific meshes.
Preserve exact region/material assignments, boundary orientation, and support traces.
Geometry integration, structural mass, contact shape, and radiating surface must agree.
Keep numerical meshing error distinct from analytic chart geometry.
G1: volume/mass/inertia and structural convergence for multiple admitted profiles.
G3: rigid transform, scale, and equivalent representation comparisons.
Test holes, thin regions, poor quality, and unsupported shapes with named refusals.
Acceptance: a non-cylinder disc profile reaches structural modes and boundary radiation.

### MR39 — Share evolving surface state across domain laws

Group: F — Interfaces and physical observations.
Priority: 2. Ambition: [F].
Depends on: MR05, MR17.
Existing anchor: dry-tribology, lubrication, and porous-capillary owners.
Files: existing surface/interface descriptors and L3 coupling state adapters.
Outcome: wear, oxide, contamination, wetness, and finish evolve once and affect all relevant laws.
Keep reference interface cards immutable and physical surface history in the run.
Use region-local roughness/lay/coating/temperature state with explicit measurement bandwidth.
Provide thin adapters to friction, heat/electrical contact, wetting, and optical models.
Do not derive every response from Ra or collapse their different calibration domains.
G3: one surface-state change invalidates and updates the consumers that depend on it.
G1: a simple sourced wear or coating-growth law with known integral evolution.
G4: failed update retains the prior state for all consumers.
Acceptance: a finish/wear change affects both contact and optical observation causally.

### MR40 — Resolve physically supported optical and thermal appearance

Group: F — Interfaces and physical observations.
Priority: 2. Ambition: [S/F].
Depends on: MR04, MR09, MR39.
Existing anchor: `frankensim-h7xu5.4` and existing fs-render optical owners.
Files: `fs-material` optical state queries, `fs-render`, scene material adapters.
Outcome: material, temperature, and surface changes determine measured-spectrum optical response.
Reuse spectral complex-index conductor and dielectric paths with correct state/source bindings.
Add only the required coating/absorption/scattering models for admitted specimens.
Thermal emission uses temperature and directional/spectral emissivity in the model's domain.
Kirchhoff reciprocity requires matched direction/wavelength/polarization and equilibrium assumptions.
A visible n+k table is not automatically a broadband thermal-emissivity model.
G1: Fresnel and blackbody/gray-body limiting cases with integrated power checks.
G3: specimen temperature/finish changes leave camera exposure separate from physical radiance.
Acceptance: hot metal appearance and cooling use consistent supported radiative inputs.

### MR41 — Render the actual deformed and phase-evolving surface

Group: F — Interfaces and physical observations.
Priority: 2. Ambition: [F].
Depends on: MR08, MR21, MR38, MR40.
Existing anchor: Euler render-scene bridge and general geometry observation owners.
Files: existing deformed-surface adapters, `fs-render` scene bridge, acoustic boundary extraction.
Outcome: base plate deformation, string sag, and later melt geometry appear in the rendered shape.
Replace the cinematic whole-base scalar translation where a spatial displacement field exists.
Use the same accepted geometry for contact, radiation traces, and rendering.
Recompute normals/visibility consistently with the current configuration.
A phase surface from MR27/MR29 must enter through the generic geometry boundary.
Do not use visual-only displacement maps as proof of simulated deformation.
G3: compare extracted vertices/normals and landmarks against the solved field.
Test time interpolation/exposure without creating physically nonexistent geometry at topology changes.
Acceptance: a spatially bending plate's pixels and acoustic boundary describe the same shape.

### MR42 — Synchronize pressure, radiance, and sensor observations

Group: F — Interfaces and physical observations.
Priority: 2. Ambition: [S/F].
Depends on: MR05, MR35, MR36, MR41.
Existing anchor: existing coupling PCM/WAV and cinematic render owners.
Files: shared accepted-time observation adapters, audio resampling and frame output paths.
Outcome: video and audio observe the same physical trajectory with declared sensor models.
Retain pressure in Pa before an explicit PCM scale; never peak-normalize each material take.
Specify microphone location/direction, propagation delay, sample rate, and antialiasing.
Specify camera exposure interval and radiance/display conversion independently.
Keep artistic soundtrack/taper/chirp paths outside the predictive observation mode.
Unsupported source tails stop or remain marked incomplete rather than being extrapolated musically.
G1: known moving/impulsive source timing and amplitude through observation/resampling.
G3: render block/frame partition changes preserve physical timestamps.
Acceptance: a material/thermal change appears at physically consistent times in both outputs.

### MR43 — Couple electrical transport and heating to material state

Group: G — Broader physical material response.
Priority: 2. Ambition: [F].
Depends on: MR06, MR07, MR09, MR19.
Existing anchor: `frankensim-ext-em-forces-losses-o9im`, electrical/circuit domain owners.
Files: existing or chartered EM/circuit transport owner, material law nodes, coupling adapters.
Outcome: conductivity/resistivity, contact resistance, and temperature drive electrical/thermal response.
First implement a bounded stationary Ohmic conductor/contact case through current field primitives.
Use tensor transport only when admitted; preserve field, frequency, phase, and temperature domains.
Account for Joule loss once; moving conductors later require conductor-frame heat and Lorentz work.
Do not demand the whole deferred EM portfolio to demonstrate stationary electrothermal coupling.
G1: heated resistor with independently calculated R(T), current, and thermal balance.
G0: charge continuity and passive transport dissipation.
G3: boundary voltage/current control and material swap with the same physical geometry.
Acceptance: a conductor heats and changes electrical response through the shared state runtime.

### MR44 — Couple magnetic response, losses, and temperature

Group: G — Broader physical material response.
Priority: 3. Ambition: [F].
Depends on: MR06, MR09, MR10, MR43.
Existing anchor: `frankensim-ext-em-2d-magnetostatics-zcdl`, `...ext-em-forces-losses-o9im`.
Files: chartered magnetic field/loss owners, fs-material laws, thermal coupling adapters.
Outcome: a named magnetic material has field/history/temperature-dependent response and heating.
Implement a bounded current field formulation with its gauge/boundary/force conventions.
Resolve B-H/hysteresis/demagnetization data and preserve irreversible versus reversible changes.
Curie temperature is not a substitute for an entire magnetization-versus-temperature law.
Do not add empirical core-loss maps to already resolved hysteresis/eddy losses.
G1: analytical magnetic limits and a sourced hysteresis loop/thermal case.
G3: material and temperature sweeps with consistent stored/dissipated field energy.
Acceptance: a simple coil/core or pickup-like system changes physical response without name presets.
The existing domain prerequisite must be implemented before this coupling task is claimed complete.

### MR45 — Implement coupled active-material reference components

Group: G — Broader physical material response.
Priority: 3. Ambition: [F].
Depends on: MR06, MR10, MR21, MR43, MR44.
Existing anchor: `frankensim-ext-active-materials-wccu`.
Files: existing constitutive graph coupled blocks, solid/EM/circuit/thermal adapters.
Outcome: piezoelectric, magnetostrictive, and thermoelectric effects connect physical domains.
Implement separable sub-deliveries under the same owner, each with a real consumer.
Declare independent thermodynamic variables, tensor frames, and open/short-circuit conventions.
Check stability/reciprocity for the complete coupled block, including reversible skew terms.
Thermoelectric heat includes the declared Seebeck/Peltier/Thomson convention.
Magnetostriction can excite structural sound; it is not a synthesized hum oscillator.
G0/G1: direct/converse piezo, Kelvin relation, and a magnetostrictive strain reference.
G3: energy exchange and orientation reversal under the correct parity convention.
Acceptance: each supported active component evolves two domains and reports a physical observable.

### MR46 — Evolve moisture, diffusion, and porous response

Group: G — Broader physical material response.
Priority: 2. Ambition: [F].
Depends on: MR10, MR16, MR19, MR23.
Existing anchor: `frankensim-ext-porous-capillary-biz5`.
Files: existing diffusion/porous/solid material adapters and coupling state fields.
Outcome: wood/cane/porous materials absorb moisture, swell, and change response over time.
Distinguish ambient RH, vapor activity, and specimen moisture content on a stated basis.
Use sourced sorption/hysteresis, diffusion, swelling, and state-dependent mechanical laws.
Do not instantaneously assign equilibrium moisture throughout a thick specimen.
Reuse current JCA(L) rigid-frame acoustics within its domain; frame waves are a separate rung.
G1: slab diffusion and directional swelling reference problems.
G3: thickness/time scaling and dry/wet cycles with history where modeled.
Acceptance: a wood/cane geometry changes dimension and modal response under a humidity history.

### MR47 — Couple chemical and surface degradation to physical state

Group: G — Broader physical material response.
Priority: 3. Ambition: [F].
Depends on: MR06, MR09, MR19, MR39, MR43.
Existing anchor: `frankensim-ext-gas-chemistry-ladder-paqh`, solid-life and thermochem owners.
Files: thermochemical closures, surface/history laws, coupled species/energy transport.
Outcome: oxidation, corrosion, or decomposition changes composition/geometry and relevant properties.
Implement one bounded sourced mechanism first, with explicit environment and reaction products.
Retain element/charge conservation, reaction/transport limits, and reference enthalpy conventions.
Galvanic corrosion requires electrolyte/geometry/potential kinetics, not a bulk corrosion score.
Oxide thickness affects heat/contact/optics only through independently admitted laws.
Do not claim service-life prediction from a constitutive formula or one accelerated test.
G1: closed reaction/transport balance and a known growth-law reference.
G3: absent reagent/potential removes the corresponding modeled reaction source.
Acceptance: one evolving chemical surface changes a real physical consumer's response.

### MR48 — Expose declarative material-driven scenarios

Group: H — Integrated journeys and qualification.
Priority: 1. Ambition: [S/F].
Depends on: MR04, MR12, MR13.
Existing anchor: `fs-project`/`fs-scenario` owners and existing Euler/music entry points.
Files: current project/scenario parsing, CLI/adapters, material catalog bindings.
Outcome: users choose geometry, material/process, environment, excitation, and observers in data.
Route the real cinematic and acoustic execution boundaries through exact material/card selection.
Keep low-level research-model scalar inputs explicit, separate from catalog-supported mode.
Version any frozen schema change with the registry's required migration/refusal path.
Emit all missing requirements before expensive work when they are knowable at admission.
Preserve units/seeds/budgets/versions/capabilities and idempotency at mutating boundaries.
G3: changing only a material binding preserves the declared comparison policy.
G4: unsupported combinations refuse before partially writing a run.
Acceptance: a user can run two materially different objects through one parameterized entry path.

### MR49 — Deliver the ordinary material-swap disc journey

Group: H — Integrated journeys and qualification.
Priority: 1. Ambition: [F].
Depends on: MR15, MR17, MR35, MR37, MR38, MR42, MR48.
Existing anchor: Euler scientific/cinematic programs and their current trajectory owner.
Files: existing Euler scenario fixtures, production adapters, e2e script/tests.
Outcome: supported metal/wood/crystal disc/base swaps change solved motion, pressure, and radiance.
Use actual material/interface binding rather than private friction/normal/rolling parameter construction.
Exercise radius, thickness, edge profile, base support, and environment independently.
Declare fixed geometry/mass/energy policies for comparisons; never silently compensate density.
Include at least one audible in-band disc mode and one non-cylinder supported profile.
G1/G3: mass/inertia, limiting mechanics, modal frequency, pressure scale, and state identity checks.
G2: retain independent measured quantities where available; no ranking from fitted spin time alone.
Acceptance: the same physical pipeline handles the admitted swap matrix with explicit gaps elsewhere.
No new scientific Euler claim is promoted solely by a rendered video.

### MR50 — Deliver fretted-string and duct material-swap journeys

Group: H — Integrated journeys and qualification.
Priority: 1. Ambition: [F].
Depends on: MR15, MR16, MR31, MR32, MR33, MR36, MR42, MR48.
Existing anchor: generic acoustic owner and existing music program successors.
Files: geometry/material/gesture fixtures and existing music render e2e lane.
Outcome: users hear physics change under material, dimension, support, force, and gas-state edits.
Use a plucked/fretted string/body case and a geometry-defined reed or jet-driven duct case.
Compare fixed support/pre-stress against explicitly retuned initialization as separate experiments.
Preserve sympathetic body/string state and finite-duration physical excitation.
Test a rigid-wall limit in which wall-material influence is correctly small.
No note-specific sound-speed change, arbitrary spectral filter, or gain normalization is permitted.
G1/G3: pitch/decay/contact/impedance changes trace to independently checked physical parameters.
Acceptance: instrument labels occur only in fixture descriptions and user-facing names.
Listening approval remains distinct from physical/numerical validation.

### MR51 — Deliver heating and melting with recoupled sound and images

Group: H — Integrated journeys and qualification.
Priority: 2. Ambition: [F].
Depends on: MR14, MR29, MR34, MR40, MR41, MR42, MR48, MR49, MR50.
Existing anchor: thermal-domain, remap, Euler, and generic acoustic owners.
Files: existing e2e scenario lanes plus bounded generic heated-member and melt cases.
Outcome: a lead disc and lead string soften, deform, partly melt, and continue physically.
Ambient heating must pass through solved heat transfer and finite latent-heat uptake.
Retain motion, stress, contact, phase, free surface, source pressure, and radiance on one clock.
Stress release, droplet separation, and loss of tension are computed outcomes.
No temperature-triggered animation, prerecorded sound, or disappearing mass is allowed.
G1/G3: Stefan balance, mass/momentum/energy closure, mesh/time refinement, and state transfer.
G2: independent supported experiments bound claims; unavailable liquid acoustics stays explicit.
Acceptance: both journeys use the same thermal/phase/transfer operators and actual emitted sources.
If a required source/model is absent, the complete melting journey remains open.

### MR52 — Calibrate constitutive inputs and assess independent observations

Group: H — Integrated journeys and qualification.
Priority: 2. Ambition: [S/F].
Depends on: MR13, MR48.
Existing anchor: existing material-identifiability, V&V, and Euler experiment owners.
Files: current calibration/UQ/V&V interfaces and bounded experiment fixtures.
Outcome: physical prediction quality is assessed independently of fitted inputs.
Calibrate coupon/interface laws using source/specimen-matched experiments.
Freeze parameters before holdout disc/string/duct comparisons.
Measure multiple observables: temperature/front/shape, force/motion, impedance/pressure, radiance.
Use pressure/mass/radius/ring/support sweeps to distinguish Euler loss mechanisms.
Do not use one stopping time or chirp exponent to identify several confounded dissipation laws.
Preserve uncertainty, covariance, alignment, metrology, and model-form discrepancy.
Acceptance tolerances are preregistered from measurement uncertainty and intended use.
G2: held-out assessment reports intervals and no-data outcomes without automatic validation promotion.
Acceptance: at least one calibration/holdout pair is executed through the real consumer model.

### MR53 — Qualify integrated runtime cost, determinism, and recovery

Group: H — Integrated journeys and qualification.
Priority: 2. Ambition: [S/F].
Depends on: MR49, MR50, MR51, MR52.
Existing anchor: existing DSR/runtime/performance and music budget owners.
Files: existing e2e scripts, runtime budget adapters, scoped contracts.
Outcome: integrated physical runs have measured cost and recover correctly under interruption.
Profile actual coupled journeys before selecting a new optimization or accelerator.
Use multirate stepping only with conservative exchange and observable-error accounting.
Measure same-ISA/thread-count determinism and record cross-ISA differences separately.
Exercise cancellation, budget exhaustion, resume, failed model fit, and failed material query.
Logs need source/state/model IDs, physical time, error norms, iterations, and stop reason.
Do not build a second general harness; extend the existing product lanes.
Run narrow RCH probes while iterating and one coherent DSR quality gate at delivery.
Acceptance: retained runtime receipts describe the tested domain and unresolved limitations.
No real-time performance claim is made until its actual machine/workload band passes.

### MR54 — Extend failure, fatigue, and long-term material response by consumer need

Group: G — Broader physical material response.
Priority: 3. Ambition: [F].
Depends on: MR09, MR10, MR22, MR28, MR39.
Existing anchor: `frankensim-ext-solid-life-ladder-gahl` and seed dataset owner.
Files: existing material/solid damage and history laws; source tranches and component tests.
Outcome: common machine materials can accumulate damage, crack, wear, and lose functionality.
Deliver bounded law/component slices: cyclic plasticity, fatigue, crack growth, brittle fracture,
creep-fatigue, and environmental damage only when the selected model/data supports them.
Retain test method, load spectrum, mean stress, rate, temperature, process, and specimen population.
Strength/hardness/ductility observations parameterize appropriate laws; they are not interchangeable.
Regularize softening/fracture with a stated length/energy scale to avoid mesh-dependent failure.
G1/G2: coupon law reference, fracture-energy balance, and held-out component observations.
G3: mesh/objectivity and transferred-history checks across damage evolution.
Acceptance: each added family changes an actual component outcome; catalog rows alone do not close it.

## 9. Milestones, evidence, and costs

### 9.1 Delivery groups

| Group | Scope | First useful result |
|---|---|---|
| A | MR01–MR08 | Correct environment/material response with coherent runtime state |
| B | MR09–MR18 | Queryable, source-backed bundles usable by physical consumers |
| C | MR19–MR24 | Spatial heating, stress, creep, and validity-aware evolution |
| D | MR25–MR30 | Computed melting/free surface and conservative state transitions |
| E | MR31–MR36 | Generic contact-excited and gas-excited sound production |
| F | MR37–MR42 | General contact and matching physical observations |
| G | MR43–MR47, MR54 | Electrical, magnetic, active, moisture, chemical, and life effects |
| H | MR48–MR53 | Runnable declarative journeys and measured qualification |

The groups organize work; they must not serialize independent tasks unnecessarily.
MR01, MR02, MR04, MR09, and MR12 can start without waiting for melting research.
MR48 can expose supported ordinary scenarios while phase-change work remains underway.
An integrated milestone depends only on the specific source and solver deliveries it uses.
Broad catalog completion, every EM regime, and every deferred theorem are not prerequisites.

### 9.2 Unit and end-to-end evidence

Use existing unit/integration seams for pure query, constitutive, and operator behavior.
Use manufactured/analytical references for implementation accuracy and refinement.
Use independent measured data for physical validation over its actual specimen/domain.
Pair deterministic repeat tests with an independent physical check so a golden cannot freeze a bug.

For each integrated run retain concise, diagnosable records:

- scenario and comparison policy; geometry/material/interface/model identities;
- source claim selection and the actual state points used;
- mesh, time, nonlinear, coupling, reduction, and observation error controls;
- accepted physical time, conserved quantities, loss channels, and event outcomes;
- observer geometry, pressure/radiance units, filtering/exposure conventions;
- compute budget, machine profile, interruption/restart outcome, and artifact locations.

These are fields in existing result/log conventions, not a new evidence platform.
Test-only plots and diagnostic sensors must not mutate the production trajectory.
An approximate pressure waveform cannot be labeled an exact sound-wave solution.
A simulation-only measurement deck cannot replace independent experimental acquisition.

### 9.3 Performance and fidelity

Begin with bounded offline reference journeys and measure them.
Exploit spatial locality, precomputed immutable property resolution, and cache-valid operators.
Resolve temperature/material state at quadrature or region updates, not via SQL per audio sample.
Separate slow transport, mechanical/contact, and acoustic observation time scales explicitly.
Use reduced models when their observables and transition errors have evidence.
Retain a higher-fidelity comparison path and an honest unavailable result.
No speedup, memory footprint, or real-time throughput number is promised by this plan.

### 9.4 Repository verification policy

This planning task changes documentation and Beads, not Rust implementation.
It therefore needs plan/reference/graph checks, not a workspace compile.
Future Rust tasks run relevant formatting and focused tests under repository policy.
Use RCH for narrow heavy probes and DSR for a coherent repository gate.
Before heavy Mac lanes run the required disk-pressure preflight.
Before attributing constellation failures to source, run
`scripts/ci/checkout_constellation.sh --verify-only` and classify pin movement safely.
Never reset a sibling checkout to make a gate pass.
No Simulator action or audio playback is needed for plan verification.

## 10. Scientific reference basis

These sources ground design decisions; they are not licenses to copy entire datasets.
Source acquisition tasks must independently retain exact terms and scientific context.
The initial research consulted these primary sources on 2026-09-04/05 UTC.

| Reference | Design consequence |
|---|---|
| [NIST, uncertainty data models](https://www.nist.gov/programs-projects/data-models-expression-uncertainty-materials-data) | Preserve uncertainty semantics and measurement metadata rather than one confidence scalar |
| [USDA Forest Products Laboratory, Wood Handbook](https://research.fs.usda.gov/fpl/wood-handbook) | Wood structure, moisture, and mechanical properties require named condition and orientation |
| [NIST Chemistry WebBook, lead](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7439921&Mask=6) | Inspect phase-specific thermochemical references; a melting point alone cannot define a melt model |
| [Voller and Prakash, 1987](https://doi.org/10.1016/0017-9310(87)90317-6) | Enthalpy methods can couple latent heat and mushy-region flow through declared closures |
| [Bilbao, Torin, Chatziioannou, collisions](https://arxiv.org/abs/1405.2589) | Energy-based collision treatment is reusable across strings, reeds, and striking interactions |
| [Thorne et al., Euler-disc investigation](https://arxiv.org/html/2603.14520v1) | Different loss mechanisms require controlled physical sweeps; fitted scaling is not a universal law |
| [NASA thermodynamic species coefficients](https://ntrs.nasa.gov/citations/20020085330) | Reuse existing species/reference-state conventions and source-normalization paths |

The Euler paper's rolling/air distinction depends on regime and experimental controls.
Its suggested adhesion mechanism is not adopted as a proven universal interface law.
The enthalpy literature motivates a numerical formulation, not arbitrary alloy calibration.
The collision literature supports numerical energy methods, not a universal tissue/contact exponent.

## 11. Fresh-eyes review and integration record

The review requirement is four plan passes, followed by bounded Bead refinement passes.
Each pass asks the owner's exact question: once again check for blunders, mistakes,
errors, oversights, omissions, problems, misconceptions, and bugs, thoroughly.
Reviews must revise actual defects; repeated statements of confidence do not count.

### Review pass 1 — Architecture and current reality

Completed during source investigation and first draft.
Found existing fs-matdb-store, taxonomy, transient conduction, and lumped enthalpy.
Revised the plan to extend them rather than propose duplicate implementations.
Separated ConstitutiveGraph's local DAG from closed coupled residual solves.
Preserved generic acoustic owners and excluded instrument-name dispatch.

### Review pass 2 — Physics and state-transition audit

Pending review of the complete task graph.
Must specifically inspect latent heat, density/volume, stress release, phase fractions,
energy/entropy ports, radiation loss duplication, damping assumptions, and observation timing.

### Review pass 3 — Independent implementability and integration review

Pending independent reading of the full plan and a standalone obscure task.
Must inspect actual source owners, missing dependencies, complete inputs/outputs,
and whether either flagship could pass through a hidden hand-tuned shortcut.

### Review pass 4 — Dependency and scope review

Pending after review revisions.
Must inspect the DAG, existing owner/deferred overlaps, five architectural rationales,
and whether the final revision introduces only bounded clarifications rather than a redesign.

### Existing-task instruction reconciliation

`frankensim-rc-root-q61wp.37` currently proposes an `fs-instrument` extraction.
The newest owner instruction requires generic physical mechanisms; its destination
must be revised to mechanism-owned decomposition with no functionality loss.
`frankensim-h7xu5.7` permits artistic extrapolation/chirps in its earlier cinematic scope.
Those are excluded from the new predictive lane; current in-progress work is not overwritten.
`frankensim-h7xu5.7.8` retains an old fs-acoustics description despite newer real radiation code.
Add precise scope reconciliation without claiming its current owner completed this new program.

## 12. Beads conversion and handoff

Pending conversion after plan review.
Use only `br` for creation, updates, comments, and dependency mutations.
Use typed parent-child hierarchy and explicit prerequisite `blocks` edges.
Use `bv --robot-*` for diagnosis; never invoke its TUI.
Keep closed antecedents closed and retain assignments on existing active work.
Do not create branches, worktrees, scratch clones, commits, or pushes for this task.

The final handoff will record the root/group/task IDs, reused owner links,
cycle-check outcome, ready implementation entries, and verification limitations.
