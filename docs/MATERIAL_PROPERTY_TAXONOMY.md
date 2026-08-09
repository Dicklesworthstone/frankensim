# Material Property Taxonomy

Normative property vocabulary for the FrankenSim material database
(beads `frankensim-oecdy` store, `frankensim-0er85` breadth; owner
directive 2026-08-08). The goal it serves: **material-swap fidelity** —
keep a simulation identical, change only what something is made of, and
have the dynamics, the sound, the looks, and the operational behavior
of the machine change the way a real experiment would.

Two doctrine rules govern every row:

1. **Properties are functions, not constants.** Nearly everything below
   depends on temperature at minimum, and often on pressure, frequency,
   strain, strain rate, humidity, field, surface state, or history. The
   pack schema's scalar-plus-validity-axes and Curve claims carry this;
   a value stored without its validity window is a defect. The
   extreme-regime premise (melt the strings) means the taxonomy spans
   phase boundaries: solid properties end at solidus, and the liquid
   phase has its own rows.
2. **Populate on demand, from license-compatible sources, honestly.**
   This taxonomy is the vocabulary, not a promise that every material
   carries every row. Absent data is absent (query refusal), never
   interpolated from vibes. Tier and uncertainty metadata are part of
   the value.

## 1. Mechanical — elastic

| Property | Why swap-fidelity needs it | Consumed by |
|---|---|---|
| Elastic constants: E, G, K, ν (isotropic); 9-constant orthotropic; up to 21 anisotropic | Stiffness of everything; wood/composites/single crystals are strongly anisotropic | fs-plate, fs-solid, fs-modal, waveguides |
| Temperature dependence E(T), G(T) | Modulus softening detunes structures long before melting | thermal-structural coupling |
| Pressure dependence of moduli | Deep submergence, geophysics, presses | fs-solid |
| Storage/loss modulus vs frequency E'(f), E''(f) | Viscoelastic stiffness IS frequency-dependent (polymers, wood) | visco models (ybc75), sound |

## 2. Mechanical — strength, plasticity, failure

| Property | Why | Consumed by |
|---|---|---|
| Yield: tension Fty, compression Fcy, shear Fsu, bearing Fbry/Fbru | The tension/compression asymmetry the owner named; joints bear | fs-solid plasticity |
| Ultimate Ftu; elongation at break; reduction of area | Ductility decides bend-vs-snap | failure models |
| Hardening law (Ramberg–Osgood n, isotropic/kinematic split, Bauschinger) | Cyclic machine loads | fs-material plastic |
| Strain-rate sensitivity (Johnson–Cook class), adiabatic-shear tendency | Impacts, crashes, machining | dynamic solid |
| Hardness — SCALE-TYPED: HB, HRC, HV, Shore (polymers), Janka (wood) | Contact, wear, fret/finger indentation; never silently converted between scales | fs-contact |
| Fracture toughness K_IC, J_IC; impact energy (Charpy); ductile–brittle transition temperature | Cold rockets: FERRITIC/carbon steel that is tough at 20 °C shatters at LOX temperature (austenitic 304/316 has no transition — exactly why it is the standard cryo choice; the swap between them is the point) | fracture |
| Fatigue: S–N curves, endurance limit, mean-stress (Goodman) behavior, Paris-law da/dN constants | Every rotating machine's life | life prediction |
| Creep: Norton exponent + activation energy, stress-rupture curves, relaxation | Turbine/exhaust-valve regimes; string tension relaxation | fs-time + solid |
| Damping: loss factor η(f, T, amplitude), specific damping capacity, thermoelastic damping inputs (α, cp, k already elsewhere) | THE decisive sound property; also machine vibration | fs-modal, vibroacoustics |
| Residual-stress state typical of process; springback | As-built behavior differs from handbook coupons | as-built modeling |

## 3. Thermodynamic and thermal

| Property | Why | Consumed by |
|---|---|---|
| Density ρ(T); porosity | Mass, inertia, buoyancy; integrates α(T) | everything |
| Thermal expansion α(T), anisotropic where real | Tolerances, thermal stress, microscope drift, fret buzz on a hot day | fs-conduction coupling |
| Specific heat cp(T) | Transient heating rates | fs-conduction, fs-time |
| Thermal conductivity k(T), anisotropic (wood: L vs R vs T) | Heat paths in machines; also thermoelastic damping | fs-conduction |
| Melting: solidus/liquidus, latent heat of fusion; glass transition Tg; allotropic transformation temperatures + latent heats | The phase boundary the doctrine wants to cross; pure iron's α→γ transition sits at 1185 K (carbon steels begin transforming lower, from the 1000 K eutectoid) | phase-change models |
| Emissivity ε(λ, T, finish) | Radiative cooling AND incandescent appearance (hot steel glows Planck × ε) | fs-conduction radiation, rendering |
| Maximum service temperature; decomposition/charring onset (wood!), pyrolysis kinetics | Wood instruments near a fire; ablatives | degradation |
| Flammability: ignition temperature, heat of combustion, oxygen index | Safety-relevant realism | combustion |
| Vapor pressure, outgassing rate | Vacuum machines — electron microscopes die from outgassing elastomers | vacuum systems |
| Ablation: heat of ablation, char yield | Rocket nozzles/heat shields | rocket thermal |

## 4. Liquid/extreme-state (past the phase boundary)

| Property | Why | Consumed by |
|---|---|---|
| Melt viscosity μ_liq(T), melt surface tension, melt density | The melting lead string SAGS with these; solder joints, casting | fluid-structure |
| Volume change on melting | Cast shrinkage, freezing damage | phase-change |
| Equation of state (Mie–Grüneisen class), Hugoniot Us–up data | Shocks: detonations, hypervelocity impact, rocket transients | shock physics |
| Vaporization: boiling point(p), latent heat; ionization energies | Plasma/exhaust regimes | rocket exhaust |
| For propellants/fuels: heat of combustion, burn-rate law r = a·pⁿ, flame temperature | An ICE or rocket is a materials-and-chemistry machine | combustion |

## 5. Fluids interaction, moisture, diffusion

| Property | Why | Consumed by |
|---|---|---|
| Surface energy; contact angle vs named liquids (hydrophobicity — owner-named); work of adhesion | Wetting, coating, condensation on optics | interface packs |
| Equilibrium moisture content vs RH; hygroscopic expansion coefficients; moisture diffusivity | Wood swells/detunes with weather — first-order for instruments | wood modeling |
| Gas/liquid permeability; diffusion coefficients (incl. hydrogen) | Seals leak; hydrogen embrittles steel | seals, embrittlement |
| Solubility/compatibility vs named fluids (fuels, oils, refrigerants) | Engine seals and hoses | chemical compat |

## 6. Acoustic (mostly derived — stored for validation)

| Property | Why | Consumed by |
|---|---|---|
| Sound speeds: longitudinal, shear, bar/plate; characteristic impedance ρc | Derived from E, G, ρ — stored values serve as cross-check gates (the existing c_L gate pattern) | validation |
| Ultrasonic attenuation α(f) | NDT realism, high-f damping | wave models |
| Radiation ratio c_L/ρ | Soundboard selection metric | instrument descriptions |
| Porous-absorber parameters: flow resistivity, open porosity, tortuosity, viscous + thermal characteristic lengths (Johnson–Champoux–Allard / Biot set) | THE first-order sound properties of felt, foam, textile, wool — piano-hammer felt, case linings, room/enclosure absorbers; a "material" that is mostly air needs its own acoustics | porous-media acoustics |

## 7. Electrical

| Property | Why | Consumed by |
|---|---|---|
| Resistivity/conductivity σ(T); temperature coefficient; RRR (cryo) | Windings heat, sensors drift; copper vs aluminium generator | electromagnetics |
| Dielectric constant εr(f, T); loss tangent tan δ(f); dielectric strength; surface/volume resistivity | Insulation systems, capacitors, RF parts | insulation |
| Piezoelectric d/g coefficients; pyroelectric coefficient | Sensors, quartz, igniters | transducers |
| Seebeck coefficient | Thermocouples inside every machine | instrumentation |
| Work function; secondary-electron yield | Electron sources and detectors — the microscope's heart | e-optics |
| Electrochemical/galvanic-series potential | Dissimilar-metal contact corrodes — material-swap can CREATE corrosion | degradation |
| Superconducting Tc, Hc (where applicable) | Cryo machines | magnets |
| Semiconductor properties: band gap, carrier mobilities, dopant behavior (where applicable) | Sensors and electronics inside machines; distinguishes silicon from everything above | devices |

## 8. Magnetic

| Property | Why | Consumed by |
|---|---|---|
| B–H curve: initial/max permeability, saturation Bs, remanence Br, coercivity Hc | The generator core IS this curve | electromagnetics |
| Core losses: Steinmetz coefficients, eddy/anomalous split (M19 pack precedent) | Efficiency and heating | machines |
| Curie temperature | A hot generator LOSES its magnetism — exactly the owner's swap-fidelity scenario | thermal-magnetic |
| Magnetostriction λs | Transformer hum — a magnetic property that makes SOUND | vibroacoustics |
| Permanent magnets: BHmax, demag curves, temperature coefficients of Br/Hc (N42 pack precedent) | Motor/generator strength vs temperature | machines |

## 9. Optical and appearance ("the looks")

| Property | Why | Consumed by |
|---|---|---|
| Spectral complex refractive index n(λ) + k(λ) | Gold looks gold, copper looks copper because of dispersion — THE first-principles color of conductors | rendering (fs-render/fs-img), optics |
| Spectral transmission/absorption; scattering coefficients (subsurface) | Glass, marble, varnished wood translucency | rendering |
| BRDF parameters: microfacet roughness from surface finish (Ra/Rz), specular/diffuse split, anisotropy (brushed), clearcoat | Same alloy, different polish, different look — finish is a state axis | rendering |
| Stress-optic coefficient (birefringence) | Photoelasticity; polarized-light microscopy | optics |
| Emissivity ε(λ, T) (shared with thermal) | Incandescence and IR imaging | rendering + thermal |
| Oxidation/patina/tarnish kinetics and their optical effect | Copper turns green, silver tarnishes — appearance evolves | aging |
| Fluorescence, thermochromism where real | Special materials | rendering |

## 10. Tribological and surface

| Property | Why | Consumed by |
|---|---|---|
| Friction coefficients static/kinetic vs COUNTER-MATERIAL, finish, lubrication state, T, velocity (Stribeck behavior) | Owner-named; friction is a PAIR + state property, carried by interface packs (existing precedent) | fs-contact, fs-tribo |
| Wear coefficients (Archard), abrasion resistance, galling tendency | Machine life; fret wear | wear models |
| Achievable roughness per finish process (links to BRDF and friction) | The polish axis that couples looks AND friction | shared state axis |
| Contact stiffness/damping of joints | Joints dominate machine vibration | structural dynamics |
| Adhesion/stiction | MEMS, seals, cold-welding in vacuum | contact |
| Coefficient of restitution vs impact velocity (measured, pair + geometry context) | System-dependent, so stored as VALIDATION data for contact models rather than a material constant — but first-order for perceived impact dynamics | contact validation |

## 11. Chemical/environmental degradation

| Property | Why | Consumed by |
|---|---|---|
| Corrosion rates vs environment; pitting resistance (PREN); passivation behavior | Stainless GRADES differ exactly here | life models |
| High-T oxidation kinetics (scale growth) | Exhaust valves, turbine blades | ICE/turbo |
| Stress-corrosion cracking thresholds; hydrogen-embrittlement susceptibility | Sudden failure modes | fracture |
| UV/radiation degradation rates; polymer embrittlement | Space, outdoors | aging |
| Rot/fungal/insect resistance (wood); biocompatibility | Wood service life; medical | wood, bio |

## 12. Process/microstructure metadata (parameterizes everything above)

| Property | Why | Consumed by |
|---|---|---|
| Condition/temper (annealed, H02, T6, cold-drawn %) as an IDENTITY axis | "304 stainless" is not one material; the corpus already distinguishes conditions | pack identity |
| Grain size, texture/anisotropy from processing (rolling, AM build direction) | Directional properties in "isotropic" metals | anisotropy |
| Hardenability (Jominy), weldability, castability/shrinkage, machinability index | Whether the simulated part can exist as designed | manufacturability |
| Hall–Petch constants, precipitation state and aging curves | Strength drifts with thermal history — links to the extreme-regime doctrine | history-dependent props |

## 13. Uncertainty and basis metadata (every row)

| Metadata | Why |
|---|---|
| Source tier (measured-primary / handbook / derived / graph-read), citation, license | Existing doctrine |
| Scatter: COV, or A-basis/B-basis allowables where the source provides them (MIL-HDBK-5 concept) | Design-grade vs typical values are different claims |
| Validity axes per value: T, p, f, strain, strain rate, RH, field, finish, age | Rule 1 above |
| Batch/specimen context (12 % MC wood convention precedent) | Reproducibility |

## Machine-level sanity map

Which rows a material swap exercises in the owner's example machines:

- **ICE engine**: 2 (fatigue, creep, hardness), 3 (α mismatch, k), 4
  (combustion), 5 (oil compatibility), 7 (spark/sensors), 10 (ring/bore
  tribology), 11 (exhaust oxidation).
- **Electric generator**: 8 (B–H, core loss, Curie, magnetostriction →
  hum), 7 (σ(T), insulation), 2 (rotor fatigue), 10 (bearings), 3
  (cooling).
- **Rocket**: 4 (EOS, combustion, ablation), 2 (creep; ductile–brittle
  screening for ferritic alloys at cryo), 3 (outgassing), 11 (LOX
  compatibility).
- **Microscope**: 3 (α → drift, outgassing), 2 (damping → vibration), 7
  (work function), 9 (n+k, birefringence), 10 (stage stiction).
- **Guitar/clarinet**: 1 (orthotropic E), 2 (η(f), hardness of
  fingertip/fret), 3 (hygroscopic + α), 5 (moisture), 6 (c_L, ρc), 9
  (finish look), 10 (bow/string, fret friction).

## Population strategy

Vocabulary-first, data-on-demand: the store (`frankensim-oecdy`) admits
every property name above with typed units and validity axes now; packs
populate rows as license-compatible sources are verified
(`frankensim-0er85` tranches first). A property a simulation needs but
the corpus lacks is a NAMED refusal at query time — which is itself the
signal for the next data tranche. No row is ever fabricated to make a
simulation run.
