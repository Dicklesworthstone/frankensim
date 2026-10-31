import "./native-wind.css";
// Canonical native fixtures, bundled verbatim by Vite; no copied physics defaults.
import reedExampleSource from "../../../crates/fs-couple/examples/reed-two-medium-tone-hole.performance?raw";
import flueExampleSource from "../../../crates/fs-couple/examples/flue-nolle.performance?raw";
import flueDgSource from "../../../crates/fs-couple/examples/flue-dg-baffled.performance?raw";
import resolvedDgSource from "../../../crates/fs-couple/examples/resolved-dg.performance?raw";
import pressureResolvedDgSource from "../../../crates/fs-couple/examples/resolved-dg-pressure.performance?raw";
import warmRampResolvedDgSource from "../../../crates/fs-couple/examples/resolved-dg-warm-ramp.performance?raw";
import splitApertureResolvedDgSource from "../../../crates/fs-couple/examples/resolved-dg-split-aperture.performance?raw";
import sideApertureResolvedDgSource from "../../../crates/fs-couple/examples/resolved-dg-side-aperture.performance?raw";
import radiatingMouthResolvedDgSource from "../../../crates/fs-couple/examples/resolved-dg-radiating-mouth.performance?raw";
import separateOpeningsResolvedDgSource from "../../../crates/fs-couple/examples/resolved-dg-separate-openings.performance?raw";
import multipleEscapesResolvedDgSource from "../../../crates/fs-couple/examples/resolved-dg-multiple-escapes.performance?raw";
import finiteReservoirResolvedDgSource from "../../../crates/fs-couple/examples/resolved-dg-finite-reservoir.performance?raw";
import actuatedReservoirResolvedDgSource from "../../../crates/fs-couple/examples/resolved-dg-actuated-reservoir.performance?raw";
import coupledEscapesResolvedDgSource from "../../../crates/fs-couple/examples/resolved-dg-coupled-escapes.performance?raw";
import slitEdgeResolvedDgSource from "../../../crates/fs-couple/examples/resolved-dg-slit-edge.performance?raw";
import siSlitEdgeResolvedDgSource from "../../../crates/fs-couple/examples/resolved-dg-si-slit-edge.performance?raw";
import upstreamResolvedDgSource from "../../../crates/fs-couple/examples/resolved-dg-upstream.performance?raw";
import molecularFlueSource from "../../../crates/fs-couple/examples/flue-dg-molecular.performance?raw";
import nolleMolecularFlueSource from "../../../crates/fs-couple/examples/flue-dg-nolle-molecular.performance?raw";
import conicalMolecularFlueSource from "../../../crates/fs-couple/examples/flue-dg-conical-molecular.performance?raw";
import fieldMolecularFlueSource from "../../../crates/fs-couple/examples/flue-dg-field-molecular.performance?raw";
import materialReedSource from "../../../crates/fs-couple/examples/material_modal_valve/paired-walls.performance?raw";
import materialReedDgSource from "../../../crates/fs-couple/examples/material_modal_valve/paired-dg-chimney.performance?raw";
import tensorReedDgSource from "../../../crates/fs-couple/examples/material_modal_valve/paired-tensor-dg.performance?raw";
import solidReedDgSource from "../../../crates/fs-couple/examples/material_modal_valve/solid-dg.performance?raw";
import hereditarySolidReedSource from "../../../crates/fs-couple/examples/material_modal_valve/solid-hereditary-dg.performance?raw";
import surfaceSolidReedSource from "../../../crates/fs-couple/examples/material_modal_valve/solid-surface-dg.performance?raw";
import stringExampleSource from "../../../crates/fs-couple/examples/material_strings/pluck-and-bow.material-strings?raw";
import measuredStringSource from "../../../crates/fs-couple/examples/material_strings/measured-cello-g.material-strings?raw";
import nonlinearPreloadStringSource from "../../../crates/fs-couple/examples/material_strings/nonlinear-preload.material-strings?raw";
import stoppedNonlinearStringSource from "../../../crates/fs-couple/examples/material_strings/nonlinear-stopped-preload.material-strings?raw";
import tensorBodyStringSource from "../../../crates/fs-couple/examples/material_strings/tensor-body.material-strings?raw";
import observedStringSource from "../../../crates/fs-couple/examples/material_strings/rayleigh-observed.material-strings?raw";
import regionalCavityStringSource from "../../../crates/fs-couple/examples/material_strings/regional-cavity.material-strings?raw";
import layeredStringSource from "../../../crates/fs-couple/examples/material_strings/relaxing-layered-strings.material-strings?raw";
import shearStringSource from "../../../crates/fs-couple/examples/material_strings/shear-layer-bow.material-strings?raw";
import compliantRibbonSource from "../../../crates/fs-couple/examples/material_strings/compliant-ribbon.material-strings?raw";
import finiteGasStringSource from "../../../crates/fs-couple/examples/finite_gas_string/browser-feedback-flow.string?raw";
import { nativeWindStatus, renderNativePerformance, type NativeWindStatus, type NativeWindAttachment, type NativeWindGasSources, type NativeStringProfiles, type NativeCoatingStrainGuard, type NativeFiniteGasStringSources, type NativeMaterialManifest, type NativeRenderer } from "./native-wind-client.mjs";
import { el } from "./util";
import { createNativeGasControls } from "./native-wind-gases";
import { createNativeMaterialManifests } from "./native-material-manifests";
import { createNativeAtmosphereControls } from "./native-atmosphere";
import { createNativeStringProfiles } from "./native-string-profiles";
import { createNativeFiniteGasStringControls } from "./native-finite-gas-string";
import { createNativeComparisonPanel } from "./native-comparison";
import { MAX_REQUEST_BYTES, editableAttachmentText, editedAttachmentBytes, prepareNativeTextExperiment, serializeNativeExperiment, validateInput } from "./native-experiment.mjs";
import stringStarterSources from "./native-string-starters.json";

// Bundle the original text assets, loading only those used by the chosen
// experiment. The native compiler still owns parsing and material resolution.
const starterSourceFiles = import.meta.glob<string>([
  "../../../crates/fs-couple/examples/material_strings/*.tsv",
  "../../../crates/fs-couple/examples/material_strings/*.body",
  "../../../crates/fs-couple/examples/material_strings/*.stops",
  "../../../crates/fs-couple/examples/material_modal_valve/tensor-blade*.tsv",
  "../../../data/matdb/seed-v1/cello-g-strings-whole-lampis-2024/*.tsv",
  "../../../data/matdb/seed-v1/bow-hair-horse-pitteroff-woodhouse-acustica1998/*.tsv",
], {query: "?raw", import: "default"});
const stringStarters: Readonly<Record<string, {
  sources: Readonly<Record<string, string>>;
  materialManifests: NativeMaterialManifest[];
  stringProfiles?: NativeStringProfiles;
}>> = stringStarterSources;

const EXAMPLES = [
  {
    name: "reed-two-medium-tone-hole.performance",
    label: "Reed — two-medium tone hole",
    renderer: "wind",
    source: reedExampleSource,
    decimate: false,
    intro: "Starter: two-medium tone hole — a short boundary diagnostic with a reed and three tube sections in different gas regions. Change the gas mixtures, temperature, geometry or gestures in the description.",
    help: "ambient sets the inlet and chimney gas; a segment’s gas entry overrides its local medium; exterior sets the surrounding gas. The interior and exterior must share a static pressure because this boundary model omits mean flow. This short numerical fixture is not a measured clarinet or oboe.",
  },
  {
    name: "flue-nolle.performance",
    label: "Flute — Nolle/Price jet",
    renderer: "wind",
    source: flueExampleSource,
    decimate: false,
    intro: "Starter: flute jet — a short blowing-pressure phrase through a cylindrical tube, using the native Nolle/Price jet transfer. Change the gas, tube, mouth geometry or pressure gesture in the description.",
    help: "gas owns the internal medium; exterior sets the surrounding gas. jet_transfer selects the short-channel Nolle/Price source. The native description owns the mechanical and output clocks, and the jet and acoustic mouth are solved together. This synthetic example is not a measured flute.",
  },
  {
    name: "flue-dg-baffled.performance",
    label: "Flute — resolved duct and material atmosphere",
    renderer: "wind",
    source: flueDgSource,
    decimate: true,
    intro: "Starter: a flute jet drives a spatial duct with a radiating outlet and mouth. Change the bore and exterior gas compositions independently; the native pressure observer propagates each source through the declared exterior atmosphere.",
    help: "dg_exterior supplies density, sound speed, viscosity, conductivity and heat capacities to the shared Stokes–Kirchhoff propagation law. observer_gas selects that law and its numerical budgets; lossless is an explicit comparison. Native decimation observes the unchanged 96 kHz mechanics at 48 kHz. This short integration fixture uses an interim jet model and is not a measured flute.",
  },
  {
    name: "resolved-dg.performance",
    label: "Resolved gas — compressible fields and coupled duct",
    renderer: "resolved",
    source: resolvedDgSource,
    decimate: true,
    intro: "Starter: a small resolved gas field coupled to a spatial duct. Edit the original grid, stationary solids, physical initial states and incoming characteristic history, then render through the native flow solver.",
    help: "This canonical source observes pressure at 96 kHz, with native decimation selected for 48 kHz playback. Changing the gas retains grid geometry and may change the physical flow/duct timestep; the observation clock stays fixed. Attached thermochemical sources and duct-wall materials use the same controls as other performers. This exploratory field fixture is not a complete flute or reed instrument.",
  },
  {
    name: "resolved-dg-pressure.performance",
    label: "Resolved flow — physical blowing pressure",
    renderer: "resolved",
    source: pressureResolvedDgSource,
    decimate: true,
    intro: "Starter: a declared pressure history drives the original resolved gas field and spatial bore. The returning wave and actual gas properties determine inlet flow; the native exterior microphone observes the resulting pressure.",
    help: "boundary pressure gives pascals relative to the original source-face pressure. Change the gas, stationary geometry or pressure history in the source, keeping the declared numerical limits. This is an ideal pressure boundary in the existing local-linear characteristic domain, not finite airway storage or a measured flute. Native decimation preserves the physical solver clock.",
  },
  {
    name: "resolved-dg-warm-ramp.performance",
    label: "Resolved flow — warm breath ramps",
    renderer: "resolved",
    source: warmRampResolvedDgSource,
    decimate: true,
    intro: "Starter: pressure and incoming gas temperature rise and fall through the original resolved flow, bore and exterior microphone. Edit the ramp endpoints to explore breath onset and release.",
    help: "Each ramp acts on the physical source at the start of its existing solver interval. Pressure is in pascals relative to the source reference and temperature is in kelvin; output samples are not enveloped or normalized. This short synthetic source is not a measured player or flute.",
  },
  {
    name: "resolved-dg-split-aperture.performance",
    label: "Resolved flow — bore aperture and escaping air",
    renderer: "resolved",
    source: splitApertureResolvedDgSource,
    decimate: true,
    intro: "Starter: a pressure-driven gas field can flow into a declared bore aperture or escape through the surrounding opening. Change the aperture bounds, matching bore area or surrounding pressure to explore their effect on the coupled flow and microphone pressure.",
    help: "coupled_patch selects original boundary cells; its area must match the bore inlet. The remaining cells keep their explicit pressure or incoming-wave condition. The microphone observes bore-outlet radiation; escaped mouth flow has no separate radiation observer yet. This small numerical fixture is not measured flute geometry or a qualified turbulent jet. Native decimation preserves the physical solver clock.",
  },
  {
    name: "resolved-dg-side-aperture.performance",
    label: "Resolved flow — transverse opening to the atmosphere",
    renderer: "resolved",
    source: sideApertureResolvedDgSource,
    decimate: true,
    intro: "Starter: a pressure-driven gas field has a separate opening in its side wall, while the original end face feeds the bore. Change the gas or opening cells to explore ventilation through the same resolved flow; the microphone observes both openings.",
    help: "reservoir_aperture reference binds the quiescent opening to the original gas material, so pressure, temperature and composition changes also update its density. Explicit density, velocity and temperature remain available for an independent reservoir. Cells retain their original wall and characteristic collars; reservoir mass and energy have separate accounts. This coarse prescribed-boundary fixture does not qualify flute geometry or jet turbulence.",
  },
  {
    name: "resolved-dg-radiating-mouth.performance",
    label: "Resolved flow — mouth and bore radiation feedback",
    renderer: "resolved",
    source: radiatingMouthResolvedDgSource,
    decimate: false,
    intro: "Starter: the resolved side opening and bore outlets push back through one shared exterior radiation model. Change the source gas, matching exterior atmosphere, opening geometry or pressure history to change the physical flow and emitted pressure together.",
    help: "reservoir_feedback solves the original mouth pressure and bore flow together, including their mutual radiation. The reference aperture follows the original source gas; the joint exterior must use the same material state. Common pressure edits reach both gas records. This short 8 kHz numerical fixture retains its declared quadrature and work limits and is not a measured flute.",
  },
  {
    name: "resolved-dg-separate-openings.performance",
    label: "Resolved flow — separate breath and acoustic escape",
    renderer: "resolved",
    source: separateOpeningsResolvedDgSource,
    decimate: false,
    intro: "Starter: warm air enters through one transverse opening while a distinct opening exchanges sound with the exterior. The same resolved gas field couples to the bore and a finite duct on the cork side.",
    help: "boundary aperture selects the driven opening by its declared zero-based index; reservoir_observation aperture selects the escape opening that receives reservoir_feedback. Edit the pressure and temperature history, original gas data or opening cells. The microphone combines escape and bore radiation with their original paths. This short numerical integration fixture still requires physical geometry, grid and jet validation.",
  },
  {
    name: "resolved-dg-finite-reservoir.performance",
    label: "Resolved flow — finite breath chamber",
    renderer: "resolved",
    source: finiteReservoirResolvedDgSource,
    decimate: false,
    intro: "Starter: a finite 10 mL nitrogen chamber feeds one resolved opening while a separate opening radiates alongside the bore and cork duct. The chamber's mass and energy evolve with the actual flow, changing its pressure and temperature.",
    help: "Load the supplied nitrogen species, NASA9 and transport packs through Gas source files, with frozen phase and no water species. finite_reservoir declares the chamber volume and exchange limits; its boundary hold carries only an endpoint. Change volume or initial aperture density and temperature to alter the source. This 8 ms numerical fixture has no breathing pump, mean gas transport through the acoustic ducts, or measured flute qualification. See docs/RESOLVED_NSF_DG_COUPLING.md for the source-pack commands.",
  },
  {
    name: "resolved-dg-multiple-escapes.performance",
    label: "Resolved flow — multiple radiating openings",
    renderer: "resolved",
    source: multipleEscapesResolvedDgSource,
    decimate: true,
    intro: "Starter: the same resolved gas field feeds a bore and two atmospheric openings. The microphone combines each original signed tile flow with its declared placement, delay and exterior gas propagation.",
    help: "Repeat reservoir_observation aperture with each original aperture index. Changing the selected gas updates all exterior paths while retaining the original grid. Multiple observers require an explicit aperture index for each reservoir_feedback reaction. This starter observes both openings without adding reaction; the mutually coupled openings starter demonstrates the shared exterior feedback. These numerical fixtures do not establish measured flute geometry.",
  },
  {
    name: "resolved-dg-actuated-reservoir.performance",
    label: "Resolved flow — chamber work and heat",
    renderer: "resolved",
    source: actuatedReservoirResolvedDgSource,
    decimate: false,
    intro: "Starter: compress and expand a finite nitrogen chamber while adding heat. Its original gas law converts the work and heat into pressure and temperature changes that drive the resolved opening, bore, cork duct and separate radiating escape.",
    help: "Attach the same nitrogen species, NASA9 and transport packs as the finite-chamber starter. boundary finite-reservoir-actuation holds or ramps volume rate in m³/s and heat power in watts; positive rate expands and positive power heats. The native solver integrates those controls on its physical clock and checks the gas, work and exchange budgets before accepting each interval. This 8 ms numerical fixture uses a well-mixed source; it does not supply player anatomy or measured flute validation.",
  },
  {
    name: "resolved-dg-coupled-escapes.performance",
    label: "Resolved flow — mutually coupled openings",
    renderer: "resolved",
    source: coupledEscapesResolvedDgSource,
    decimate: false,
    intro: "Starter: a finite actuated nitrogen chamber drives the resolved fluid, bore and cork duct while two separate openings react through one shared exterior radiation field. Each opening's pressure responds to the other openings' flow.",
    help: "Attach the same nitrogen species, NASA9 and transport packs as the chamber starters. Each reservoir_feedback row selects an observed opening; reservoir_feedback_solve supplies the joint solve limits. The microphone retains every original tile's signed flow and material propagation path. All coupled pressure and work checks precede acceptance. This 8 ms numerical fixture needs native execution and physical refinement before instrument claims.",
  },
  {
    name: "resolved-dg-slit-edge.performance",
    label: "Resolved flow — slit jet and solid edge",
    renderer: "resolved",
    source: slitEdgeResolvedDgSource,
    decimate: true,
    intro: "Starter: blowing pressure drives air through a solid slit nozzle, across a free gap and past a solid edge. The same fluid field feeds a bore aperture and a separate atmospheric opening; the microphone observes the bore outlet.",
    help: "This manufactured geometry uses a 0.5 mm grid, a 1 mm slit and an explicit solid wedge. Edit solids, pressure or gas to change the original flow; no jet card or noise generator supplies the sound. The physical solver keeps its own clock while native decimation produces 48 kHz playback. This short coarse example needs grid and interface refinement before edge-tone or instrument claims.",
  },
  {
    name: "resolved-dg-si-slit-edge.performance",
    label: "Resolved flow — slit and edge in physical coordinates",
    renderer: "resolved",
    source: siSlitEdgeResolvedDgSource,
    decimate: true,
    intro: "Starter: edit the nozzle, slit and solid edge in metres while the original gas solver resolves their flow. The atmospheric opening follows the selected gas material and the bore retains its original radiation and microphone.",
    help: "solid_geometry rasterizes the declared boxes and wedge on the selected grid; changing geometry or mesh data changes the original solid mask. reservoir_aperture_m reference binds the opening's physical bounds and original gas state without reauthoring cells or density. This manufactured 0.5 mm example remains subject to grid, interface and jet refinement; native decimation preserves the physical clock.",
  },
  {
    name: "resolved-dg-upstream.performance",
    label: "Resolved flow — upstream acoustic feedback",
    renderer: "resolved",
    source: upstreamResolvedDgSource,
    decimate: true,
    intro: "Starter: a resolved gas source sits between a finite upstream duct and a downstream bore. Change the upstream geometry, gas or pressure-rise history to change the returning waves and the resulting downstream microphone pressure.",
    help: "upstream_dg declares the original upstream acoustic load. boundary upstream-pressure supplies a pressure rise in pascals between that duct and the resolved zone; both face pressures and flows remain coupled unknowns. Material walls and gas relaxation use the existing source records. This acoustic interface does not transport mean mass or enthalpy into the ducts, and this manufactured fixture does not qualify a complete player airway or flute. Native decimation preserves the 96 kHz solver clock.",
  },
  {
    name: "flue-dg-molecular.performance",
    label: "Flue tube — molecular gas relaxation (exploration)",
    renderer: "wind",
    source: molecularFlueSource,
    decimate: true,
    intro: "Starter: a generic flue-driven tube uses one explicitly declared terrestrial-air relaxation law in its spatial bore and exterior microphone paths. Native decimation retains the 96 kHz mechanics. This numerical fixture uses an interim jet model; it is not a measured flute.",
    help: "dg_molecular and observer_gas molecular select the shared native relaxation law. The ISO terrestrial-air declaration has a limited composition domain: enriching oxygen or adding other gases may be refused until a compatible molecular-law declaration is supplied. The atmosphere editor changes the original gas records and preserves the selected law; it never substitutes an air law or invents molecular coefficients.",
  },
  {
    name: "flue-dg-nolle-molecular.performance",
    label: "Flute — Nolle jet, viscous windway and molecular gas",
    renderer: "wind",
    source: nolleMolecularFlueSource,
    decimate: true,
    intro: "Starter: the existing Nolle/Price jet source drives a spatial bore through a viscous windway. The bore and microphone paths retain molecular gas relaxation. Change the original atmosphere, windway, mouth or bore data and render through the shared native owners.",
    help: "jet_transfer selects the small-signal long-channel source; windway slit derives the channel flow from gas density and viscosity. Internal and exterior gases remain separate. The ISO terrestrial-air relaxation law has a bounded composition domain. This numerical fixture is not a resolved turbulent jet or a measured flute; native decimation retains its 96 kHz mechanics.",
  },
  {
    name: "flue-dg-conical-molecular.performance",
    label: "Conical flue — material gas dynamics (exploration)",
    renderer: "wind",
    source: conicalMolecularFlueSource,
    decimate: true,
    intro: "Starter: a flue drives a conical bore with coupled molecular and wall heat exchange. Change the cone geometry or its uniform gas data, then render the resulting pressure through the original exterior atmosphere.",
    help: "The dg_segment record owns the taper; dg_molecular owns the selected gas relaxation law. The atmosphere controls preserve both. The ISO terrestrial-air law may refuse enriched or additional species until a compatible law is supplied. This example uses a local cylindrical wall approximation and an interim jet source; it does not establish measured oboe, saxophone or flute accuracy. Native decimation retains the 96 kHz mechanics.",
  },
  {
    name: "flue-dg-field-molecular.performance",
    label: "Wind bore — temperature and composition gradients",
    renderer: "wind",
    source: fieldMolecularFlueSource,
    decimate: true,
    intro: "Starter: a conical bore with independently editable temperature and composition stations. Attach the N2 and CO2 species, NASA-9 and transport packs using the axial-gas example recipe. Each local gas state enters the native molecular and wall dynamics; the exterior has its own atmosphere.",
    help: "Select a gas_station in Regional atmosphere to edit its temperature and mole fractions. The dg_bore header owns the shared species columns, and gas owns the common static pressure. This exploration uses manufactured relaxation coefficients and an interim jet; it is not a measured breath field or a qualified instrument. Native decimation retains the 96 kHz mechanics.",
  },
  {
    name: "material_modal_valve/paired-walls.performance",
    label: "Material reeds — paired blades and compliant duct walls",
    renderer: "wind",
    source: materialReedSource,
    decimate: false,
    intro: "Starter: two material-derived reed blades coupled to an upstream tract, downstream bore and material walls. Attach blade.matpack and wall.matpack from the canonical material-modal example. Change each blade region or wall source independently.",
    help: "The original meshes, thickness regions, material stiffness and loss determine blade motion, aperture flow and external microphone pressure. This source already observes its 768 kHz mechanics at 48 kHz; leave additional decimation off. Its manufactured materials are not measured cane or a qualified clarinet/oboe specimen.",
  },
  {
    name: "material_modal_valve/paired-dg-chimney.performance",
    label: "Material reeds — resolved duct, chimney and regional gases",
    renderer: "wind",
    source: materialReedDgSource,
    decimate: true,
    intro: "Starter: the same material-derived paired blades drive a spatially resolved duct with a side chimney. Attach blade.matpack. The bore, chimney and exterior have separately declared gas sources and feed the original radiating boundary.",
    help: "This canonical source runs at 96 kHz, with native decimation selected for 48 kHz playback. observer_gas propagates the outlet and chimney sound through dg_exterior using the shared Stokes–Kirchhoff law. Change the material pack or gas data in the source; the same native solver derives the pressure output.",
  },
  {
    name: "material_modal_valve/paired-tensor-dg.performance",
    label: "Material reeds — full directional relaxation",
    renderer: "wind",
    source: tensorReedDgSource,
    decimate: true,
    intro: "Starter: two material reeds carry a complete three-dimensional relaxation tensor through the original plate motion, spatial bore and exterior microphone. Build and attach tensor-blade.matpack using the canonical material-modal source recipe, then change its directional stiffness or relaxation-time data.",
    help: "loss tensor-maxwell names the source strain limit and complete Mandel branch coefficients, including explicit zeros. The native solver retains coupled normal-strain histories and every declared material point. Equilibrium elasticity stays owned by the same card; material, state-count and clock limits may refuse an unsupported change. This coarse manufactured specimen does not establish measured cane or clarinet/oboe accuracy.",
  },
  {
    name: "material_modal_valve/solid-hereditary-dg.performance",
    label: "Material reed — three-dimensional material memory",
    renderer: "wind",
    source: hereditarySolidReedSource,
    decimate: true,
    intro: "Starter: the original tetrahedral solid retains its source-derived relaxation history while driving the same spatial bore and exterior microphone. Build and attach solid-valve.matpack using the material-modal source recipe, then edit its material properties, geometry, supports or loaded faces.",
    help: "solid_memory selects relaxed or unrelaxed initial material history and explicit projection budgets. Original tensor branches, source relaxation times, signed pressure work and every selected mode enter the existing exact structural solver. Source strain, clock and capacity limits remain active. This small manufactured soft solid does not establish measured cane or instrument accuracy.",
  },
  {
    name: "material_modal_valve/solid-surface-dg.performance",
    label: "Material reed — original solid surface contact",
    renderer: "wind",
    source: surfaceSolidReedSource,
    decimate: true,
    intro: "Starter: a three-dimensional material reed touches its lay along selected original solid faces. Build and attach solid-valve.matpack using the canonical source recipe, then change the material, local clearances or physical slit edges to hear their effect through the same bore and exterior microphone.",
    help: "solid_closure spatial retains each original nodal normal displacement, including modes whose average opening is zero. Partial slit closure, contact forces and material memory enter the shared native solve. This path uses fixed reference Tet4 faces with parallel normals; it does not resolve a moving fluid channel or establish measured reed accuracy. Native decimation preserves the 96 kHz mechanics.",
  },
  {
    name: "material_modal_valve/solid-dg.performance",
    label: "Material reed — three-dimensional solid",
    renderer: "wind",
    source: solidReedDgSource,
    decimate: true,
    intro: "Starter: an original tetrahedral solid drives the same spatial bore and exterior microphone as the material-plate reeds. Build and attach solid-valve.matpack using the material-modal source recipe, then edit the solid geometry, supports, loaded faces or source properties.",
    help: "The native solid owner derives mass, stiffness, signed pressure work and the complete coupled viscosity from the declared material tensor. Every selected mode remains, including modes with zero mean opening. This path uses a frozen-frequency material image with explicit dispersion checks; it does not evolve full solid relaxation histories or establish measured reed accuracy. Native decimation preserves the mechanical clock.",
  },
  {
    name: "material_strings/pluck-and-bow.material-strings",
    label: "Strings — pluck and bow with a shared body",
    renderer: "strings",
    source: stringExampleSource,
    decimate: false,
    intro: "Starter: one material-derived string, body and cavity, with a physical preload and bow gesture. The original wire and body sources are included and compile when you render. Edit their attached tables or open your own complete material-string performance.",
    help: "The string performer accepts the same full material cards as the reed and wall solvers, along with legacy claim packs. Replace an attached source at its declared path to change the solved material without browser coefficients. Optional profiles add body geometry, fingers/frets and contact heat. This is a numerical integration fixture, not a measured violin, cello or bass.",
  },
  {
    name: "material_strings/measured-cello-g.material-strings",
    label: "Strings — measured cello G coefficients on a laboratory plate",
    renderer: "strings",
    source: measuredStringSource,
    decimate: false,
    intro: "Starter: published whole-string mass and bending-stiffness estimates drive a manufactured laboratory plate. The original Lampis source table, attribution, plate sources and compilation bindings are included. Render with the native material compiler configured, then edit a source table to compare.",
    help: "The measured setup and property selectors describe one published string on a laboratory receiver. Its temperature_known=0 stays independent of body and gas temperature. No intrinsic loss law or reference-fluid inertial correction is inferred; sampled gas resistance remains. A replacement measured source must supply matching length, radius, fixed tension and explicit selectors. This fixture does not establish complete cello acoustics or measured recording accuracy.",
  },
  {
    name: "material_strings/nonlinear-preload.material-strings",
    label: "Strings — force-loaded nonlinear pizzicato",
    renderer: "strings",
    source: nonlinearPreloadStringSource,
    decimate: false,
    intro: "Starter: a force-loaded material string with bonded layers, two transverse planes, torsion, geometric stretching and a coupled body/cavity. Original wire, sleeve and body sources are included. The force is released at sample zero from the solved coupled equilibrium.",
    help: "Change material sources, geometry or the declared preload force. The native solver derives the initial deformation from the same storage law used during playback, retaining the original relaxed material histories. Longitudinal-wave descriptions use this same preload path. This manufactured specimen does not establish measured violin, cello or bass accuracy.",
  },
  {
    name: "material_strings/nonlinear-stopped-preload.material-strings",
    label: "Strings — force-loaded contact and nonlinear material",
    renderer: "strings",
    source: stoppedNonlinearStringSource,
    decimate: false,
    intro: "Starter: a force-loaded layered string equilibrates against the original compliant stop while retaining geometric stretching, torsion, body and cavity motion. Original wire, sleeve and body sources are included, with stopped-pluck.stops already selected as the contact profile.",
    help: "The native solver balances material, support and normal-contact forces together before releasing the declared load at sample zero. Change material sources or the attached stop profile to change the initial shape and ensuing pressure. The supplied stop is manufactured data; material-stop profiles use source-bound normal contact through the same solver. This is not a measured violin, cello or bass.",
  },
  {
    name: "material_strings/tensor-body.material-strings",
    label: "Strings — full tensor body memory",
    renderer: "strings",
    source: tensorBodyStringSource,
    decimate: false,
    intro: "Starter: a material string drives a body with complete three-dimensional tensor relaxation. Original wire and tensor-body material sources are included, with tensor-body.body selected. Replace the profile with the canonical finer tensor-body-shared.body to explore spatial refinement.",
    help: "The native body retains all six source strain histories, including out-of-plane effects and reciprocal string reactions. The shared profile uses exact equal-time history coordinates so its finer mesh fits the same state cap; every original material-point history remains reconstructible. Change branch_33 or relaxation_time in the source and rebuild its pack to compare pressure. This fixed-state manufactured cantilever is not a measured violin, cello or bass. Its 48 kHz clock needs no decimation.",
  },
  {
    name: "material_strings/rayleigh-observed.material-strings",
    label: "Strings — material plate and exterior microphone",
    renderer: "strings",
    source: observedStringSource,
    decimate: false,
    intro: "Starter: a plucked material string drives a regional orthotropic plate, with reciprocal acoustic loading and a microphone at its declared position. Original wire and directional-body sources are included, with orthotropic.body already selected as the body profile.",
    help: "The native observer integrates signed plate motion over source-to-microphone paths through the original gas. Change material sources or atmosphere data and render again. The body uses an infinite baffle and synthetic specimen data; this flat-plate fixture does not represent a complete measured violin, cello or bass.",
  },
  {
    name: "material_strings/regional-cavity.material-strings",
    label: "Strings — independent cavity and exterior gases",
    renderer: "strings",
    source: regionalCavityStringSource,
    decimate: false,
    intro: "Starter: a material string and closed shell with a helium-filled interior and a separately declared exterior atmosphere. Original wire and shell sources are included, with closed-shell-cavity.body selected so both gas regions reach the native coupled solve.",
    help: "cavity_gas supplies the original interior acoustic field; gas supplies the exterior and string drag. Change either atmosphere independently while keeping their static pressures equal, or apply a common pressure with the atmosphere editor. The body profile is required for the cavity declaration. This manufactured shell fixture is not a measured violin, cello or bass.",
  },
  {
    name: "material_strings/relaxing-layered-strings.material-strings",
    label: "Strings — bonded layers with original relaxation spectra",
    renderer: "strings",
    source: layeredStringSource,
    decimate: false,
    intro: "Starter: layered strings with independent source relaxation spectra. Original relaxing-wire and body sources are included. Change their attached tables or the layer geometry in the performance; the shared constitutive solver retains each layer's history.",
    help: "This perfectly bonded concentric-section model retains original per-layer moduli and relaxation times. It does not add sliding helical windings. Keep the stated end-load or displacement comparison and pressure scale when comparing material substitutions.",
  },
  {
    name: "material_strings/compliant-ribbon.material-strings",
    label: "Strings — material bow-hair compliance",
    renderer: "strings",
    source: compliantRibbonSource,
    decimate: false,
    intro: "Starter: the bow drives a finite ribbon of compliant hair coupled to the original string, body and cavity. Original wire, body and horsehair sources are included. Change the attached hair material table, radius, length, contacting-hair count or contact position in the source to change its mechanical response.",
    help: "The native solver derives hair stiffness and damping from the selected material properties and geometry, retaining spring energy and contact history during gestures and release. The source explicitly selects the published quasistatic spring-dashpot continuation; it omits fiber inertia and uses a prescribed normal load. The fixture retains an authored Coulomb interface, so it does not establish rosin-contact physics or realistic violin/cello sound. Keep the pressure scale fixed when comparing variants.",
  },
  {
    name: "material_strings/shear-layer-bow.material-strings",
    label: "Strings — material shear layer and distributed bow contact",
    renderer: "strings",
    source: shearStringSource,
    decimate: false,
    intro: "Starter: a material shear layer drives distributed bow contact on the original string/body/cavity system. Original wire, body and manufactured shear-layer sources are included; editing the attached material tables changes the native contact law when rendered.",
    help: "Layer shear modulus, yield stress and plastic viscosity control retained elastic/plastic history and the actual pressure output. The native instructions describe geometric stretching, longitudinal waves and optional contact heat. Heat feedback requires the same materials to supply caloric and transport data; this isothermal fixture does not invent them.",
  },
  {
    name: "finite_gas_string/browser-feedback-flow.string",
    label: "Heated string — electrical feedback and finite gas",
    renderer: "finite-gas-string",
    source: finiteGasStringSource,
    decimate: false,
    intro: "Starter: a current heats a vibrating material string while its resistance, material loss and surrounding gas evolve. The same accepted trajectory exchanges heat and transfers gas through a nozzle. Attach the original string pack, nitrogen source packs, initial gas record and nitrogen-reservoir.gas from the native example recipe.",
    help: "electrical_feedback closes deposited Joule heat against endpoint resistance. The 48 kHz source preserves vibration and observes compact pressure while the finite gas temperature, composition and transport state change. Download its synchronized CSV for temperatures, source work, heat transfer and energy residuals. This manufactured short fixture uses first-order heat-then-flow integration and a diagnostic pressure observer; it is not measured wire data or a chamber acoustics model.",
  },
];
const exampleUrl = (name: string): string =>
  `https://github.com/Dicklesworthstone/frankensim/blob/main/crates/fs-couple/examples/${name}`;
const LOCAL_FILE_LIMIT = 1_048_576; // Bound file reading even on a static host without a renderer.
const message = (error: unknown): string => error instanceof Error ? error.message : String(error);

/** A source editor and native render client. Synthesis stays in the Rust renderer. */
export function createNativeWindPanel(): HTMLElement {
  let selectedExample = EXAMPLES[0]!;
  const renderer = el("select", {id: "native-performance-renderer", class: "btn"},
    el("option", {value: "wind"}, "Wind — reeds, jets and material valves"),
    el("option", {value: "resolved"}, "Resolved flow — physical gas fields and spatial duct"),
    el("option", {value: "strings"}, "Strings — materials, bow/contact and shared body"),
    el("option", {value: "finite-gas-string"}, "Heated string — electrical, thermal and finite-gas feedback"));
  const selectedRenderer = (): NativeRenderer =>
    renderer.value === "strings" || renderer.value === "resolved" || renderer.value === "finite-gas-string" ? renderer.value : "wind";
  const rendererAvailable = (next: NativeWindStatus | null): boolean => {
    const selected = selectedRenderer();
    return Boolean(next?.available && (selected === "wind"
      ? next.renderers?.wind !== false : next.renderers?.[selected] === true));
  };
  const starter = el("select", { id: "native-wind-starter", class: "btn" },
    ...EXAMPLES.map((example) => el("option", { value: example.name }, example.label)));
  const introduction = el("p", { class: "nw-intro" }, selectedExample.intro);
  const sourceHelp = el("p", { id: "native-wind-source-help" }, selectedExample.help);
  const canonical = el("a", { href: exampleUrl(selectedExample.name), target: "_blank", rel: "noopener" }, "Canonical example");
  const source = el("textarea", { id: "native-wind-source", class: "nw-source", rows: "15",
    spellcheck: "false", autocapitalize: "off", autocomplete: "off", wrap: "off",
    "aria-describedby": "native-wind-source-help" });
  source.value = selectedExample.source;
  // A textarea normalizes CRLF. Keep the loaded UTF-8 text until the user edits it.
  let loadedSource = selectedExample.source;
  let loadedEditorValue = source.value;
  const currentSource = (): string => source.value === loadedEditorValue ? loadedSource : source.value;
  const sourceInfo = el("span", { class: "small muted" });
  const connection = el("span", { class: "nw-connection", role: "status", "aria-live": "polite" }, "Checking local renderer…");
  const result = el("p", { class: "nw-result small", role: "status", "aria-live": "polite", "aria-atomic": "true" });
  const render = el("button", { type: "button", class: "btn primary" }, "Render WAV");
  const cancel = el("button", { type: "button", class: "btn" }, "Cancel render");
  const check = el("button", { type: "button", class: "btn" }, "Check renderer");
  const open = el("button", { type: "button", class: "btn" }, "Open performance…");
  const openExperiment = el("button", {type:"button", class:"btn"}, "Open experiment…");
  const saveExperiment = el("button", {type:"button", class:"btn"}, "Save experiment");
  const experimentFile = el("input", {type:"file", accept:".json,application/json",
    "aria-label":"Open a complete native experiment with its material sources"});
  experimentFile.hidden = true;
  const reset = el("button", { type: "button", class: "btn" }, "Reset example");
  const attach = el("button", { type: "button", class: "btn" }, "Attach source files…");
  const attachmentFile = el("input", { type: "file", multiple: "",
    "aria-label": "Attach material packs or other referenced native source files" });
  attachmentFile.hidden = true;
  const decimate = el("input", { type: "checkbox", id: "native-wind-decimate" });
  const decimation = el("label", { for: "native-wind-decimate", class: "small nw-decimation" }, decimate,
    " Use native decimation for a performance above 48 kHz (mechanics keep their original clock)");
  const attachmentList = el("div", { class: "nw-attachments" });
  const attachmentInfo = el("p", { class: "small muted" });
  const attachmentPanel = el("section", { "aria-label": "Attached source files" },
    attachmentInfo, attachmentList);
  type Attached = NativeWindAttachment & { sha256: string };
  let attachments: Attached[] = [];
  let replacing: Attached | null = null;
  const replacementFile = el("input", {type:"file", "aria-label":"Replace an attached source without changing its path"});
  replacementFile.hidden = true;
  let editingAttachment: Attached | null = null;
  let loadedAttachmentSource = "", loadedAttachmentEditorValue = "";
  const attachmentEditorTitle = el("label", {for:"native-attachment-source"});
  const attachmentSource = el("textarea", {id:"native-attachment-source", class:"nw-source", rows:"12",
    spellcheck:"false", autocapitalize:"off", autocomplete:"off", wrap:"off"});
  const currentAttachmentSource = (): string => attachmentSource.value === loadedAttachmentEditorValue
    ? loadedAttachmentSource : attachmentSource.value;
  const applyAttachment = el("button", {type:"button", class:"btn primary"}, "Apply source edit");
  const cancelAttachment = el("button", {type:"button", class:"btn"}, "Cancel edit");
  const attachmentEditor = el("section", {"aria-label":"Edit attached physical source"},
    attachmentEditorTitle,
    el("p", {class:"small muted"}, "Edit original manifest, material table, profile, geometry or control records. Apply keeps this file's path and all profile associations; the native renderer checks the physical data. Binary material, interface and model packs use Replace."),
    attachmentSource, el("div", {class:"btn-row"}, applyAttachment, cancelAttachment));
  attachmentEditor.hidden = true;
  const attachmentIdentity = (): string => JSON.stringify(attachments.map(file => [file.path, file.sha256]));
  const materials = createNativeMaterialManifests(() => { refreshSourcePaths(); edited(); });
  const gas = createNativeGasControls(() => edited());
  const profiles = createNativeStringProfiles(() => edited());
  const finiteGasString = createNativeFiniteGasStringControls(() => edited());
  const atmosphere = createNativeAtmosphereControls(currentSource, next => {
    source.value = next; loadedSource = next; loadedEditorValue = source.value;
    edited();
  }, () => controls());
  // Reuse the same record editor for original attached gas files. Applying the
  // gas draft changes this editor; Apply source edit publishes bytes and identity
  // together through replaceAttachedSource, preserving every source association.
  const attachmentAtmosphere = createNativeAtmosphereControls(currentAttachmentSource, next => {
    attachmentSource.value = next;
    loadedAttachmentSource = next; loadedAttachmentEditorValue = attachmentSource.value;
    attachmentAtmosphere.refresh(next);
    controls();
  }, () => controls());
  attachmentEditor.insertBefore(attachmentAtmosphere.element, attachmentSource);
  const file = el("input", { type: "file", accept: ".performance,.material-strings,.string,.txt,text/plain",
    "aria-label": "Open a native physical performance" });
  file.hidden = true;
  const audio = el("audio", { controls: "", preload: "metadata", "aria-label": "Your native physical render" });
  const download = el("a", { class: "btn", download: "wind-render.wav" }, "Download WAV");
  const historyDownload = el("a", {class:"btn", download:"string-history.csv"}, "Download synchronized CSV");
  historyDownload.hidden = true;
  const renderedFrom = el("p", { class: "small muted" });
  const renderedFiles = el("p", { class: "small muted nw-file-identities" });
  const playback = el("div", { class: "nw-playback" }, audio, download, historyDownload, renderedFrom, renderedFiles);
  playback.hidden = true;
  cancel.hidden = true;
  const comparison = createNativeComparisonPanel(() => audio.pause());
  audio.addEventListener("play", () => comparison.stop());

  let status: NativeWindStatus | null = null;
  let checking: AbortController | null = null;
  let rendering: AbortController | null = null;
  let loadingFile = false;
  let previewUrl: string | null = null;
  let historyUrl: string | null = null;
  let renderedSource: string | null = null;
  let renderedAttachmentIdentity: string | null = null;
  let renderedDecimate = false;
  let renderedGasIdentity: string | null = null;
  let renderedProfileIdentity: string | null = null;
  let renderedFiniteGasIdentity: string | null = null;
  let renderedMaterialIdentity: string | null = null;
  let renderedRenderer: string | null = null;
  let sourceName = selectedExample.name;

  const controls = (): void => {
    const active = rendering !== null;
    const otherEdit = active || loadingFile || editingAttachment !== null;
    const locked = otherEdit || atmosphere.hasDraft();
    render.disabled = locked || checking !== null || !rendererAvailable(status) || Boolean(status?.busy) || !source.value.trim();
    cancel.hidden = !active;
    check.disabled = active || checking !== null;
    open.disabled = locked;
    openExperiment.disabled = locked;
    saveExperiment.disabled = locked || !source.value.trim();
    reset.disabled = locked;
    starter.disabled = locked;
    renderer.disabled = locked;
    attach.disabled = locked;
    const heated = selectedRenderer() === "finite-gas-string";
    decimate.disabled = locked || heated;
    decimation.hidden = heated;
    gas.setDisabled(locked);
    materials.setDisabled(locked);
    materials.setAvailable(status?.supportsMaterialManifests === true);
    atmosphere.setDisabled(otherEdit);
    atmosphere.element.hidden = heated;
    profiles.element.hidden = selectedRenderer() !== "strings";
    profiles.setDisabled(locked || selectedRenderer() !== "strings");
    finiteGasString.element.hidden = !heated;
    finiteGasString.setDisabled(locked || !heated);
    attachmentList.querySelectorAll<HTMLInputElement | HTMLButtonElement>("input,button")
      .forEach(control => { control.disabled = locked; });
    source.readOnly = locked;
    source.setAttribute("aria-busy", String(active));
    attachmentAtmosphere.setDisabled(active || loadingFile || editingAttachment === null);
    attachmentSource.readOnly = active || loadingFile || attachmentAtmosphere.hasDraft();
    applyAttachment.disabled = active || loadingFile || attachmentAtmosphere.hasDraft();
    cancelAttachment.disabled = active || loadingFile;
    render.textContent = active ? "Rendering…" : "Render WAV";
  };
  const edited = (): void => {
    atmosphere.refresh(currentSource());
    sourceInfo.textContent = `${sourceName} · ${new TextEncoder().encode(currentSource()).byteLength.toLocaleString()} bytes`;
    if (renderedSource !== null) {
      renderedFrom.textContent = currentSource() === renderedSource && attachmentIdentity() === renderedAttachmentIdentity && decimate.checked === renderedDecimate && gas.identity() === renderedGasIdentity && selectedRenderer() === renderedRenderer && profiles.identity() === renderedProfileIdentity && finiteGasString.identity() === renderedFiniteGasIdentity && materials.identity() === renderedMaterialIdentity
        ? "This WAV was rendered from the performance and source files currently shown."
        : "This WAV is from your previous performance or source files. Render again to hear your edits.";
    }
    controls();
  };
  const clearPreview = (): void => {
    audio.pause();
    audio.removeAttribute("src");
    audio.load();
    download.removeAttribute("href");
    historyDownload.removeAttribute("href");
    historyDownload.hidden = true;
    playback.hidden = true;
    renderedSource = null;
    renderedAttachmentIdentity = null;
    renderedGasIdentity = null;
    renderedProfileIdentity = null;
    renderedFiniteGasIdentity = null;
    renderedMaterialIdentity = null;
    renderedRenderer = null;
    renderedFiles.textContent = "";
    if (previewUrl !== null) URL.revokeObjectURL(previewUrl);
    previewUrl = null;
    if (historyUrl !== null) URL.revokeObjectURL(historyUrl);
    historyUrl = null;
  };
  const showResult = (text: string, refused = false): void => {
    result.textContent = text;
    result.classList.toggle("nw-refused", refused);
  };

  const checkConnection = async (): Promise<void> => {
    checking?.abort();
    const controller = new AbortController();
    checking = controller;
    status = null;
    connection.textContent = "Checking local renderer…";
    controls();
    const timeout = window.setTimeout(() => controller.abort(), 5000);
    try {
      const next = await nativeWindStatus(controller.signal);
      if (checking !== controller) return;
      status = next;
      connection.textContent = !rendererAvailable(next) ? "Selected native renderer unavailable"
        : next.busy ? "Local renderer is busy" : "Local renderer ready";
      if (!rendererAvailable(next)) {
        showResult("Configure the selected native executable in the local Instrument Lab server to render and listen. You can edit the source here.");
      } else if (next.busy) {
        showResult("Another render is running. Check the renderer again when it finishes.");
      } else {
        showResult(`Ready to render. The local render limit is ${Math.ceil(next.timeoutMs / 1000)} seconds.`);
      }
    } catch {
      if (checking !== controller) return;
      connection.textContent = "Native renderer unavailable";
      showResult("This host has no native renderer connection. Open the local Instrument Lab server to render and listen.");
    } finally {
      window.clearTimeout(timeout);
      if (checking === controller) {
        checking = null;
        controls();
      }
    }
  };

  render.addEventListener("click", () => {
    if (!status || !rendererAvailable(status) || status.busy || rendering || checking || loadingFile) return;
    const connected = status;
    const submitted = currentSource();
    const submittedFiles = attachments.map(file => ({ ...file }));
    const submittedAttachmentIdentity = attachmentIdentity();
    const submittedDecimate = decimate.checked;
    const submittedRenderer = selectedRenderer();
    let submittedMaterials: NativeMaterialManifest[] | undefined;
    let submittedGas: NativeWindGasSources | undefined;
    let submittedProfiles: NativeStringProfiles | undefined;
    let submittedCoatingGuard: NativeCoatingStrainGuard | undefined;
    let submittedFiniteGas: NativeFiniteGasStringSources | undefined;
    try {
      submittedMaterials = materials.read();
      submittedGas = gas.read();
      submittedProfiles = submittedRenderer === "strings" ? profiles.read() : undefined;
      submittedCoatingGuard = submittedRenderer === "strings" ? profiles.readCoatingStrainGuard() : undefined;
      submittedFiniteGas = submittedRenderer === "finite-gas-string" ? finiteGasString.read() : undefined;
    }
    catch (error) { showResult(message(error), true); return; }
    const submittedGasIdentity = gas.identity();
    const submittedProfileIdentity = profiles.identity();
    const submittedFiniteGasIdentity = finiteGasString.identity();
    const submittedMaterialIdentity = materials.identity();
    const controller = new AbortController();
    rendering = controller;
    audio.pause();
    comparison.stop();
    showResult("Rendering your performance from its material, gas and geometry data…");
    controls();
    let timedOut = false;
    const timeout = window.setTimeout(() => {
      timedOut = true;
      controller.abort();
    }, connected.timeoutMs + 5000);
    void (async () => {
      try {
        const native = await renderNativePerformance(submitted, submittedFiles, connected, {renderer:submittedRenderer,
          decimate:submittedDecimate, gasSources:submittedGas, stringProfiles:submittedProfiles,
          coatingStrainGuard:submittedCoatingGuard,
          finiteGasString:submittedFiniteGas, materialManifests:submittedMaterials}, controller.signal);
        if (rendering !== controller || controller.signal.aborted) return;
        let details = submittedFiles.length
          ? "Rendered source files: " + submittedFiles.map(file => `${file.path} · ${file.bytes.byteLength.toLocaleString()} bytes · SHA-256 ${file.sha256}`).join("\n")
          : "This render used no attached source files.";
        if (submittedMaterials) details += "\nCompiled material sources: " + submittedMaterials.map(binding =>
          `${binding.packPath} from ${binding.manifestPath}` + (binding.materialState
            ? ` · material ${binding.materialState.chemistry} / ${binding.materialState.phase} / ${binding.materialState.process} / revision ${binding.materialState.revision}`
            : " · original source envelope")).join("; ");
        if (submittedGas) details += "\nGas associations: " + submittedGas.components
          .map(component => `${component.speciesPath} / ${component.nasa9Path} / ${component.transportPath}`).join("; ")
          + ` · phase ${submittedGas.phase} · water ${submittedGas.waterSpecies ?? "none"}`;
        if (submittedProfiles) details += "\nString profiles: " + Object.entries(submittedProfiles)
          .map(([role, path]) => `${role}: ${path}`).join("; ");
        if (submittedCoatingGuard) details += `\nCoating strain checks: maximum bytes ${submittedCoatingGuard.maximumBytes}, maximum work ${submittedCoatingGuard.maximumWork}`;
        if (submittedFiniteGas) details += `\nHeated string: ${submittedFiniteGas.solidPath} · initial gas: ${submittedFiniteGas.gasRecordPath}`;
        await comparison.setCurrent(native, sourceName,
          `Performer: ${submittedRenderer} · native decimation: ${submittedDecimate}\n${details}\n\nOriginal performance:\n${submitted}`, controller.signal);
        if (rendering !== controller || controller.signal.aborted) return;
        const nextUrl = URL.createObjectURL(native.wave);
        clearPreview();
        previewUrl = nextUrl;
        renderedSource = submitted;
        renderedAttachmentIdentity = submittedAttachmentIdentity;
        renderedDecimate = submittedDecimate;
        renderedGasIdentity = submittedGasIdentity;
        renderedProfileIdentity = submittedProfileIdentity;
        renderedFiniteGasIdentity = submittedFiniteGasIdentity;
        renderedMaterialIdentity = submittedMaterialIdentity;
        renderedRenderer = submittedRenderer;
        renderedFiles.textContent = details;
        audio.src = nextUrl;
        download.href = nextUrl;
        download.download = `${sourceName.replace(/\.[^.]*$/, "") || "physical"}-render.wav`;
        if (native.history) {
          historyUrl = URL.createObjectURL(native.history);
          historyDownload.href = historyUrl;
          historyDownload.download = `${sourceName.replace(/\.[^.]*$/, "") || "physical"}-history.csv`;
          historyDownload.hidden = false;
        }
        playback.hidden = false;
        edited();
        showResult("Render complete. Listen above or keep this result as reference A to compare a material or atmosphere change below.");
      } catch (error) {
        if (rendering !== controller) return;
        showResult(controller.signal.aborted
          ? timedOut ? "The render timed out. Check the local renderer before trying again." : "Render cancelled."
          : message(error), !controller.signal.aborted || timedOut);
      } finally {
        window.clearTimeout(timeout);
        if (rendering === controller) {
          rendering = null;
          controls();
        }
      }
    })();
  });
  cancel.addEventListener("click", () => {
    rendering?.abort();
    showResult("Cancelling render…");
  });
  check.addEventListener("click", () => void checkConnection());
  source.addEventListener("input", edited);
  open.addEventListener("click", () => file.click());
  const loadExample = async (example: typeof EXAMPLES[number]): Promise<void> => {
    if (rendering || loadingFile || editingAttachment || atmosphere.hasDraft()) return;
    const previousStarter = canonical.hidden ? "" : selectedExample.name;
    loadingFile = true; controls();
    showResult("Loading the original performance and its physical sources…");
    try {
      const inputs = stringStarters[example.name];
      const experiment = inputs ? await prepareNativeTextExperiment(example.source, inputs.sources, {
        renderer: example.renderer as NativeRenderer, decimate: example.decimate,
        materialManifests: inputs.materialManifests, stringProfiles: inputs.stringProfiles,
      }, async path => {
        const read = starterSourceFiles[`../../../${path}`];
        if (!read) throw new Error(`Original starter source is unavailable: ${path}`);
        return read();
      }) : {attachments: [] as NativeWindAttachment[], materialManifests: undefined, stringProfiles: undefined};
      const nextAttachments = await Promise.all(experiment.attachments.map(async attached => {
        const digest = await crypto.subtle.digest("SHA-256", new Uint8Array(attached.bytes));
        const sha256 = Array.from(new Uint8Array(digest), byte => byte.toString(16).padStart(2, "0")).join("");
        return {...attached, sha256};
      }));
      // Publish the model only after all original sources, paths and bindings
      // have loaded. A failed chunk keeps the user's current experiment intact.
      selectedExample = example;
      source.value = example.source; loadedSource = example.source; loadedEditorValue = source.value;
      sourceName = example.name.split("/").pop()!;
      renderer.value = example.renderer; starter.value = example.name;
      introduction.textContent = example.intro; sourceHelp.textContent = example.help;
      canonical.href = exampleUrl(example.name); canonical.hidden = false;
      attachments = nextAttachments;
      decimate.checked = example.decimate;
      materials.load(experiment.materialManifests);
      gas.reset(); profiles.reset(); finiteGasString.reset();
      drawAttachments();
      profiles.load(experiment.stringProfiles);
      edited();
      await checkConnection();
      showResult(inputs
        ? `${attachments.length} original source files loaded with their material compilation and profile bindings. Edit a table, atmosphere or gesture, then render with the native material compiler configured.`
        : "Example loaded. Edit the source or render it.");
    } catch (error) {
      starter.value = previousStarter;
      showResult(`Example was not loaded: ${message(error)} Your current sources are unchanged.`, true);
    } finally { loadingFile = false; controls(); }
  };
  reset.addEventListener("click", () => { void loadExample(selectedExample); });
  starter.addEventListener("change", () => {
    if (rendering || loadingFile) return;
    const selected = EXAMPLES.find((example) => example.name === starter.value);
    if (!selected) return;
    void loadExample(selected);
  });
  renderer.addEventListener("change", () => {
    if (selectedRenderer() === "finite-gas-string") decimate.checked = false;
    edited(); void checkConnection();
  });
  const refreshSourcePaths = (): void => {
    const originalPaths = attachments.map(file => file.path);
    materials.setAttachments(originalPaths);
    const sourcePaths = [...originalPaths, ...materials.packPaths()];
    gas.setAttachments(sourcePaths);
    profiles.setAttachments(sourcePaths);
    finiteGasString.setAttachments(sourcePaths);
  };
  const drawAttachments = (): void => {
    attachmentInfo.textContent = attachments.length
      ? `${attachments.length} source files · ${attachments.reduce((sum, file) => sum + file.bytes.byteLength, 0).toLocaleString()} bytes. Paths must match the performance exactly.`
      : "Attach the original material packs or artifacts referenced by your performance. Keep their source paths; edit a path below when the description uses a subdirectory.";
    attachmentList.replaceChildren();
    refreshSourcePaths();
    for (const attached of attachments) {
      const path = el("input", { type: "text", value: attached.path, spellcheck: "false",
        "aria-label": `Source path for ${attached.path}` });
      path.addEventListener("input", () => {
        attached.path = path.value; refreshSourcePaths(); edited();
      });
      const editText = el("button", {type:"button", class:"btn"}, "Edit text…");
      editText.addEventListener("click", () => {
        if (rendering || loadingFile || editingAttachment) return;
        try {
          const text = editableAttachmentText(attached.bytes);
          attachmentSource.value = text;
          loadedAttachmentSource = text; loadedAttachmentEditorValue = attachmentSource.value;
          attachmentAtmosphere.refresh(text);
          attachmentAtmosphere.element.hidden = !attachmentAtmosphere.hasRegions();
          editingAttachment = attached;
          attachmentEditorTitle.textContent = `Editing ${attached.path}`;
          attachmentEditor.hidden = false;
          controls(); attachmentSource.focus();
          showResult("Apply or cancel the attached-source edit before rendering or saving.");
        } catch (error) { showResult(message(error), true); }
      });
      const replace = el("button", {type:"button", class:"btn"}, "Replace…");
      replace.addEventListener("click", () => {
        if (rendering || loadingFile) return;
        replacing = attached; replacementFile.click();
      });
      const remove = el("button", { type: "button", class: "btn" }, "Remove");
      remove.addEventListener("click", () => {
        attachments = attachments.filter(file => file !== attached);
        drawAttachments(); edited();
      });
      attachmentList.append(el("div", { class: "nw-attachment" }, path,
        el("span", { class: "small muted", title: `SHA-256 ${attached.sha256}` },
          `${attached.bytes.byteLength.toLocaleString()} bytes · ${attached.sha256.slice(0, 12)}…`), editText, replace, remove));
    }
    controls();
  };
  const replaceAttachedSource = async (target: Attached, bytes: Uint8Array): Promise<void> => {
    const maximumFile = Math.min(8 * 1024 * 1024, status?.maximumAttachmentFileBytes ?? 8 * 1024 * 1024);
    const maximumTotal = Math.min(32 * 1024 * 1024, status?.maximumAttachmentBytes ?? 32 * 1024 * 1024);
    const otherBytes = attachments.reduce((sum, file) => sum + (file === target ? 0 : file.bytes.byteLength), 0);
    if (!bytes.length || bytes.length > maximumFile || otherBytes + bytes.length > maximumTotal)
      throw new Error("Replacement exceeds the source-file limit. The original source remains attached.");
    const digest = await crypto.subtle.digest("SHA-256", new Uint8Array(bytes));
    const sha256 = Array.from(new Uint8Array(digest), byte => byte.toString(16).padStart(2, "0")).join("");
    const index = attachments.indexOf(target);
    if (index === -1) throw new Error("The source is no longer attached.");
    // Publish bytes and identity together; paths, associations and prior audio stay intact.
    attachments[index] = {path:target.path, bytes, sha256};
    drawAttachments(); edited();
  };
  const closeAttachmentEditor = (): void => {
    editingAttachment = null;
    loadedAttachmentSource = ""; loadedAttachmentEditorValue = "";
    attachmentSource.value = "";
    attachmentAtmosphere.refresh("");
    attachmentEditor.hidden = true;
    controls();
  };
  attachmentSource.addEventListener("input", () => {
    attachmentAtmosphere.refresh(currentAttachmentSource());
    attachmentAtmosphere.element.hidden = !attachmentAtmosphere.hasRegions();
    controls();
  });
  cancelAttachment.addEventListener("click", () => {
    if (rendering || loadingFile) return;
    closeAttachmentEditor();
    showResult("Edit cancelled. The original attached source is unchanged.");
  });
  applyAttachment.addEventListener("click", () => {
    const target = editingAttachment;
    if (!target || rendering || loadingFile || attachmentAtmosphere.hasDraft()) return;
    let bytes: Uint8Array;
    try {
      bytes = editedAttachmentBytes(target.bytes, currentAttachmentSource());
    }
    catch (error) { showResult(message(error), true); return; }
    if (bytes === target.bytes) {
      closeAttachmentEditor();
      showResult("No source change. The original bytes and line endings are retained.");
      return;
    }
    loadingFile = true; controls();
    void replaceAttachedSource(target, bytes).then(() => {
      closeAttachmentEditor();
      showResult(`Updated ${target.path}. Render again to hear the new physical configuration; reference A is unchanged.`);
    }).catch(error => {
      showResult(`Edit was not applied: ${message(error)}`, true);
    }).finally(() => { loadingFile = false; controls(); });
  });
  decimate.addEventListener("change", edited);
  attach.addEventListener("click", () => attachmentFile.click());
  openExperiment.addEventListener("click", () => experimentFile.click());
  saveExperiment.addEventListener("click", () => {
    if (rendering || loadingFile) return;
    try {
      const packet = serializeNativeExperiment(currentSource(), attachments, {
        renderer:selectedRenderer(), decimate:decimate.checked, gasSources:gas.read(),
        materialManifests:materials.read(),
        stringProfiles:selectedRenderer() === "strings" ? profiles.read() : undefined,
        coatingStrainGuard:selectedRenderer() === "strings" ? profiles.readCoatingStrainGuard() : undefined,
        finiteGasString:selectedRenderer() === "finite-gas-string" ? finiteGasString.read() : undefined,
      });
      // Use the service's same envelope checks before saving an unusable packet.
      // Native physics admission remains at render time.
      validateInput(new TextEncoder().encode(packet), true);
      const url = URL.createObjectURL(new Blob([packet], {type:"application/json"}));
      const link = el("a", {href:url, download:sourceName.replace(/(?:\.native)?\.json$|\.(?:performance|material-strings|string)$/, "") + ".native.json"});
      link.click();
      window.setTimeout(() => URL.revokeObjectURL(url), 1000);
      showResult("Experiment saved with the original performance, source files, material bindings, profiles and gas associations. Open it to continue this material comparison.");
    } catch (error) { showResult(message(error), true); }
  });
  experimentFile.addEventListener("change", () => {
    const selected = experimentFile.files?.[0];
    experimentFile.value = "";
    if (!selected || rendering || loadingFile) return;
    if (selected.size > MAX_REQUEST_BYTES) {
      showResult(`Choose an experiment no larger than ${MAX_REQUEST_BYTES.toLocaleString()} bytes.`, true);
      return;
    }
    loadingFile = true; controls();
    void selected.arrayBuffer().then(async buffer => {
      const experiment = validateInput(new Uint8Array(buffer), true);
      const text = new TextDecoder("utf-8", {fatal:true, ignoreBOM:true}).decode(experiment.source);
      const nextAttachments = await Promise.all(experiment.attachments.map(async attached => {
        const digest = await crypto.subtle.digest("SHA-256", new Uint8Array(attached.bytes));
        const sha256 = Array.from(new Uint8Array(digest), byte => byte.toString(16).padStart(2, "0")).join("");
        return {...attached, sha256};
      }));
      // Commit the complete imported experiment only after every byte and
      // association was admitted. A bad file leaves the current model intact.
      attachments = nextAttachments;
      materials.load(experiment.materialManifests);
      source.value = text; loadedSource = text; loadedEditorValue = source.value;
      sourceName = selected.name;
      renderer.value = experiment.renderer ?? "wind";
      decimate.checked = experiment.decimate;
      drawAttachments();
      gas.load(experiment.gasSources);
      profiles.load(experiment.stringProfiles, experiment.coatingStrainGuard);
      finiteGasString.load(experiment.finiteGasString);
      starter.value = ""; canonical.hidden = true;
      introduction.textContent = `Loaded ${selected.name}: the complete performance and ${attachments.length} original source files.`;
      sourceHelp.textContent = "This experiment keeps its original performer, mechanical clock, pressure scale, material paths, specimen profiles and gas-source associations. Replace an attached source or edit the performance to run another physical configuration.";
      edited();
      await checkConnection();
      showResult("Experiment loaded. Render to hear it; earlier comparison audio is unchanged.");
    }).catch(error => {
      showResult(`Experiment was not loaded: ${message(error)} Your current sources are unchanged.`, true);
    }).finally(() => { loadingFile = false; controls(); });
  });
  attachmentFile.addEventListener("change", () => {
    const selected = Array.from(attachmentFile.files ?? []);
    attachmentFile.value = "";
    if (!selected.length || rendering || loadingFile) return;
    const maxFiles = Math.min(32, status?.maximumAttachments ?? 32);
    const maxFileBytes = Math.min(8 * 1024 * 1024, status?.maximumAttachmentFileBytes ?? 8 * 1024 * 1024);
    const maxBytes = Math.min(32 * 1024 * 1024, status?.maximumAttachmentBytes ?? 32 * 1024 * 1024);
    const total = attachments.reduce((sum, file) => sum + file.bytes.byteLength, 0)
      + selected.reduce((sum, file) => sum + file.size, 0);
    if (attachments.length + selected.length + materials.packPaths().length > maxFiles || total > maxBytes
        || selected.some(file => !file.size || file.size > maxFileBytes)) {
      showResult(`Choose at most ${maxFiles} nonempty files, ${maxFileBytes.toLocaleString()} bytes per file and ${maxBytes.toLocaleString()} bytes combined.`, true);
      return;
    }
    loadingFile = true; controls();
    void Promise.all(selected.map(async file => {
      const buffer = await file.arrayBuffer();
      const digest = await crypto.subtle.digest("SHA-256", buffer);
      const sha256 = Array.from(new Uint8Array(digest), byte => byte.toString(16).padStart(2, "0")).join("");
      return {path:file.name, bytes:new Uint8Array(buffer), sha256};
    })).then(added => {
      attachments.push(...added);
      drawAttachments(); edited();
      showResult("Source files attached. Check that each path matches the performance, then render.");
    }).catch(() => {
      showResult("Source files could not be read. Your previous attachments are unchanged.", true);
    }).finally(() => { loadingFile = false; controls(); });
  });
  replacementFile.addEventListener("change", () => {
    const selected = replacementFile.files?.[0], target = replacing;
    replacementFile.value = ""; replacing = null;
    if (!selected || !target || !attachments.includes(target) || rendering || loadingFile) return;
    const maximumFile = Math.min(8 * 1024 * 1024, status?.maximumAttachmentFileBytes ?? 8 * 1024 * 1024);
    const maximumTotal = Math.min(32 * 1024 * 1024, status?.maximumAttachmentBytes ?? 32 * 1024 * 1024);
    const otherBytes = attachments.reduce((sum, file) => sum + (file === target ? 0 : file.bytes.byteLength), 0);
    if (!selected.size || selected.size > maximumFile || otherBytes + selected.size > maximumTotal) {
      showResult("Replacement exceeds the source-file limit. The original material remains attached.", true);
      return;
    }
    loadingFile = true; controls();
    void selected.arrayBuffer().then(buffer => replaceAttachedSource(target, new Uint8Array(buffer)))
      .then(() => showResult(`Replaced ${target.path}. Render again to hear the new source; reference A is unchanged.`))
      .catch(error => showResult(`The replacement was not applied: ${message(error)}`, true))
      .finally(() => { loadingFile = false; controls(); });
  });
  file.addEventListener("change", () => {
    const selected = file.files?.[0];
    file.value = "";
    if (!selected || rendering || loadingFile) return;
    const limit = Math.min(LOCAL_FILE_LIMIT, status?.maximumSourceBytes ?? LOCAL_FILE_LIMIT);
    if (selected.size > limit) {
      showResult(`Choose a performance file no larger than ${limit.toLocaleString()} bytes.`, true);
      return;
    }
    loadingFile = true;
    controls();
    void selected.arrayBuffer().then((bytes) => {
      const text = new TextDecoder("utf-8", { fatal: true, ignoreBOM: true }).decode(bytes);
      source.value = text;
      loadedSource = text;
      loadedEditorValue = source.value;
      sourceName = selected.name;
      starter.value = ""; canonical.hidden = true;
      introduction.textContent = `Loaded performance: ${selected.name}. Existing source attachments and profile choices are retained.`;
      showResult("Performance loaded. Edit the source or render it.");
      edited();
    }).catch(() => {
      showResult("This performance could not be read as UTF-8 text. Your previous source is still in the editor.", true);
    }).finally(() => {
      loadingFile = false;
      controls();
    });
  });
  window.addEventListener("pagehide", () => {
    checking?.abort();
    rendering?.abort();
    clearPreview();
    comparison.dispose();
  });
  window.addEventListener("pageshow", (event) => {
    if (event.persisted) void checkConnection();
  });

  const panel = el("section", { class: "panel nw-panel", "aria-labelledby": "native-wind-title" },
    el("div", { class: "panel-head" },
      el("div", {}, el("h3", { id: "native-wind-title" }, "Render a material-driven performance"),
        el("p", { class: "small muted" }, "Edit the source, run the native physics, and listen to your result.")),
      connection),
    el("div", {class: "nw-source-head"}, el("label", {for: "native-performance-renderer"}, "Physical performer"), renderer),
    el("div", { class: "nw-source-head" }, el("label", { for: "native-wind-starter" }, "Starter performance"), starter),
    introduction,
    el("details", { class: "nw-help small" },
      el("summary", {}, "About the source and local renderer"),
      sourceHelp,
      el("p", {}, "Wind and material-string descriptions use their original native loaders. Attach the exact material packs and artifacts named by the source, including subdirectory paths. Changing a pack changes its physical properties; missing data and unsupported states remain native refusals. Keep the same declared pressure scale when comparing renders."),
      el("p", {}, "Rendering uses the local music_render, material_strings or finite_gas_string executable. The heated string performer evolves its material and finite atmosphere together and returns synchronized CSV observations with the WAV. The hosted review page can show and edit the source; audio is available through the local renderer server. Wind and material-string decimation retains the original mechanics; heated strings require a declared 48 kHz clock."),
      el("p", {}, "Open experiment loads a complete native source packet, including its original material files, manifest bindings, gas associations and specimen profiles. Save experiment keeps the current configuration in that same format; neither action renders sound. The material-string launcher can prepare these packets directly from a physical assembly."),
      el("p", {}, canonical,
        " · ", el("a", {href: "https://github.com/Dicklesworthstone/frankensim/tree/main/crates/fs-couple/examples/material_strings", target: "_blank", rel: "noopener"}, "String source-pack instructions"),
        " · ", el("a", {href: "https://github.com/Dicklesworthstone/frankensim/blob/main/crates/fs-couple/examples/material_strings/MEASURED_STRINGS.md", target: "_blank", rel: "noopener"}, "Measured string source recipe"),
        " · ", el("a", {href: "https://github.com/Dicklesworthstone/frankensim/tree/main/crates/fs-couple/examples/finite_gas_string", target: "_blank", rel: "noopener"}, "Heated string source-pack instructions"),
        " · ", el("a", { href: "https://github.com/Dicklesworthstone/frankensim/tree/main/apps/instrument-lab#native-wind-and-string-renderer", target: "_blank", rel: "noopener" }, "Local renderer setup"))),
    atmosphere.element,
    el("div", { class: "nw-source-head" }, el("label", { for: "native-wind-source" }, "Performance source"),
      sourceInfo),
    source, attachmentPanel, attachmentEditor, materials.element, profiles.element, finiteGasString.element, gas.element, decimation,
    el("div", { class: "btn-row nw-actions" }, render, cancel, openExperiment, saveExperiment, open, attach, reset, check,
      file, experimentFile, attachmentFile, replacementFile),
    result, playback, comparison.element);
  drawAttachments();
  edited();
  void checkConnection();
  return panel;
}
