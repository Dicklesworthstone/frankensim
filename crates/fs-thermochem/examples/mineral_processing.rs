//! Evaluate mineral reaction amounts, heat and Gibbs diagnostics from material packs.
//!
//! Usage: `cargo run -p fs-thermochem --example mineral_processing -- \
//! cao.fsmatpk water.fsmatpk caoh2.fsmatpk 56.077`
//! Select another endpoint by prefixing `--reaction magnesite-calcination` or
//! `--reaction magnetite-oxidation`; pack arguments must match that equation's
//! explicit source identities. Copper oxide reduction and calcite bicarbonate
//! neutralization modes take four packs; dolomite bicarbonate neutralization
//! takes five. The carbonate endpoints use compatible USGS aqueous-ion
//! reference data, with bicarbonate retained in solution.
//! Ferric hydrolysis, goethite precipitation with hydrated lime, and hydrated-
//! lime neutralization use compatible USGS goethite/portlandite source phases.
//! Gold cyanidation uses six matched historical NBS/1989-errata references,
//! including the individual aqueous CN-, Au(CN)2- and OH- ions.
//! `--help` lists every supported reaction and required pack order.
//! Prefix a source-pack invocation with `--equilibrium` to also resolve standard
//! entropy and derive Gibbs energy and log equilibrium constants from enthalpy
//! and entropy at the same source state.
//! `--gibbs` additionally requires a formation-Gibbs claim for every participant.
//! Repeated `--activity 'SPECIES=VALUE'` arguments imply `--gibbs` and evaluate
//! the local reaction Gibbs energy only when every active activity is supplied.
//! Both routes can be requested; independently tabulated source values can differ.
//! Copper electrowinning also reports its signed reversible cell potential from
//! the direct Gibbs route; this does not supply an operating cell voltage.
//! Or explicit malachite inventory accounting, with all four inputs in moles:
//! `cargo run -p fs-thermochem --example mineral_processing -- \
//! --malachite-inventory 10 12 100 3`
//!
//! In the source-pack modes, the final argument is kilograms of the first
//! declared chemical species actually reacted. For aqueous ions this is ion
//! mass, not the mass of a solution. Commercial feed purity, concentration and
//! conversion require separate evidence.
//! This is standard-state bookkeeping at 298.15 K and 1 bar, not a plant
//! temperature, reaction-rate, water-dosing or recovery prediction.

use std::{error::Error, path::Path};

use fs_matdb::{MatDbError, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Amount, Dimensionless, Dims, Mass, MolarMass, QuantitySpec};
use fs_thermochem::{
    MassAmountBasis, SpeciesId,
    stoichiometry::{
        FormationEnthalpySource, electrolysis::ElectrolysisReaction, minerals::MineralReaction,
        thermodynamics::FormationGibbsSource,
    },
};

const REFERENCE: &str = "source-elemental-reference-janaf";
const USGS_REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const NBS_REFERENCE: &str = "source-elemental-reference-nbs1989";

struct ReactionCase {
    reaction: MineralReaction,
    reference: &'static str,
    packs: &'static [&'static str],
    species: &'static [&'static str],
    phase_axes: &'static [&'static str],
    equation: &'static str,
}

const REACTIONS: [&str; 27] = [
    "lime-hydration",
    "hydrated-lime-neutralization",
    "magnesite-calcination",
    "magnetite-oxidation",
    "cuprous-oxide-oxidation",
    "cupric-oxide-reduction",
    "cupric-oxide-acid-leach",
    "cuprous-oxide-reduction",
    "copper-electrowinning",
    "gold-cyanidation",
    "pyrite-roasting",
    "pyrite-oxidation-ferrous",
    "ferrous-oxidation-ferric",
    "ferric-pyrite-oxidation",
    "ferric-hydrolysis-goethite",
    "ferric-goethite-precipitation-with-hydrated-lime",
    "sulfur-dioxide-oxidation",
    "sulfur-trioxide-hydration",
    "pyrite-to-acid",
    "calcite-bicarbonate-neutralization",
    "calcite-carbonic-acid-neutralization",
    "malachite-carbonic-acid-leach",
    "dolomite-bicarbonate-neutralization",
    "ferric-chalcopyrite-leach",
    "ferric-chalcocite-leach",
    "ferric-covellite-leach",
    "ferric-bornite-leach",
];
const USAGE: &str = "usage: mineral_processing [--reaction NAME] [--equilibrium] [--gibbs] [--activity 'SPECIES=VALUE']... [--] PACK... REACTED_SPECIES_KG\n       mineral_processing --malachite-inventory MALACHITE_MOL PROTON_MOL WATER_MOL EXTENT_MOL";

fn case(name: &str) -> Result<ReactionCase, Box<dyn Error>> {
    Ok(match name {
        "lime-hydration" => ReactionCase {
            reaction: MineralReaction::LimeHydration,
            reference: REFERENCE,
            packs: &[
                "calcium-oxide-solid-nist-mining",
                "water-liquid-janaf-mining",
                "calcium-hydroxide-solid-nist-mining",
            ],
            species: &["CaO(s)", "H2O(l)", "Ca(OH)2(s)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-liquid",
                "source-phase-solid",
            ],
            equation: "CaO(s) + H2O(l) -> Ca(OH)2(s)",
        },
        "hydrated-lime-neutralization" => ReactionCase {
            reaction: MineralReaction::HydratedLimeNeutralization,
            reference: USGS_REFERENCE,
            packs: &[
                "portlandite-solid-usgs-b2131-thermochemistry",
                "hydrogen-ion-aqueous-usgs-b2131",
                "calcium-ion-aqueous-usgs-b2131",
                "water-liquid-usgs-b2131",
            ],
            species: &["Ca(OH)2(s)", "H+(aq)", "Ca[2+](aq)", "H2O(l)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-liquid",
            ],
            equation: "Ca(OH)2(portlandite) + 2 H+(aq) -> Ca[2+](aq) + 2 H2O(l), hydrated-lime neutralization",
        },
        "magnesite-calcination" => ReactionCase {
            reaction: MineralReaction::MagnesiteCalcination,
            reference: REFERENCE,
            packs: &[
                "magnesite-solid-janaf-mining",
                "magnesium-oxide-solid-janaf-mining",
                "carbon-dioxide-gas-janaf-mining",
            ],
            species: &["MgCO3(s)", "MgO(s)", "CO2(g)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-solid",
                "source-phase-gas",
            ],
            equation: "MgCO3(s) -> MgO(s) + CO2(g)",
        },
        "magnetite-oxidation" => ReactionCase {
            reaction: MineralReaction::MagnetiteOxidationToHematite,
            reference: REFERENCE,
            packs: &[
                "magnetite-solid-janaf-mining",
                "oxygen-gas-janaf-mining",
                "hematite-solid-janaf-mining",
            ],
            species: &["Fe3O4(s)", "O2(g)", "Fe2O3(s)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-gas",
                "source-phase-solid",
            ],
            equation: "4 Fe3O4(s) + O2(g) -> 6 Fe2O3(s), hematite endpoint",
        },
        "cuprous-oxide-oxidation" => ReactionCase {
            reaction: MineralReaction::CuprousOxideOxidation,
            reference: REFERENCE,
            packs: &[
                "copper-i-oxide-solid-janaf-mining",
                "oxygen-gas-janaf-mining",
                "copper-ii-oxide-solid-janaf-mining",
            ],
            species: &["Cu2O(s)", "O2(g)", "CuO(s)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-gas",
                "source-phase-solid",
            ],
            equation: "2 Cu2O(s) + O2(g) -> 4 CuO(s)",
        },
        "cupric-oxide-reduction" => ReactionCase {
            reaction: MineralReaction::CupricOxideCarbonMonoxideReduction,
            reference: REFERENCE,
            packs: &[
                "copper-ii-oxide-solid-janaf-mining",
                "carbon-monoxide-gas-janaf-mining",
                "copper-solid-janaf-mining",
                "carbon-dioxide-gas-janaf-mining",
            ],
            species: &["CuO(s)", "CO(g)", "Cu(s)", "CO2(g)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-gas",
                "source-phase-solid",
                "source-phase-gas",
            ],
            equation: "CuO(s) + CO(g) -> Cu(s) + CO2(g)",
        },
        "cupric-oxide-acid-leach" => ReactionCase {
            reaction: MineralReaction::CupricOxideAcidLeach,
            reference: USGS_REFERENCE,
            packs: &[
                "copper-ii-oxide-solid-usgs-b2131-thermochemistry",
                "hydrogen-ion-aqueous-usgs-b2131",
                "copper-ii-ion-aqueous-usgs-b2131",
                "water-liquid-usgs-b2131",
            ],
            species: &["CuO(s)", "H+(aq)", "Cu[2+](aq)", "H2O(l)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-liquid",
            ],
            equation: "CuO(tenorite) + 2 H+(aq) -> Cu[2+](aq) + H2O(l), dissolved copper(II) endpoint",
        },
        "cuprous-oxide-reduction" => ReactionCase {
            reaction: MineralReaction::CuprousOxideCarbonMonoxideReduction,
            reference: REFERENCE,
            packs: &[
                "copper-i-oxide-solid-janaf-mining",
                "carbon-monoxide-gas-janaf-mining",
                "copper-solid-janaf-mining",
                "carbon-dioxide-gas-janaf-mining",
            ],
            species: &["Cu2O(s)", "CO(g)", "Cu(s)", "CO2(g)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-gas",
                "source-phase-solid",
                "source-phase-gas",
            ],
            equation: "Cu2O(s) + CO(g) -> 2 Cu(s) + CO2(g)",
        },
        "gold-cyanidation" => ReactionCase {
            reaction: MineralReaction::GoldCyanidation,
            reference: NBS_REFERENCE,
            packs: &[
                "gold-solid-nbs1989-thermochemistry",
                "cyanide-ion-aqueous-nbs1989",
                "oxygen-gas-nbs1989-thermochemistry",
                "water-liquid-nbs1989-thermochemistry",
                "dicyanoaurate-ion-aqueous-nbs1989",
                "hydroxide-ion-aqueous-nbs1989",
            ],
            species: &[
                "Au(s)",
                "CN[-](aq)",
                "O2(g)",
                "H2O(l)",
                "Au(CN)2[-](aq)",
                "OH[-](aq)",
            ],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-gas",
                "source-phase-liquid",
                "source-phase-aqueous",
                "source-phase-aqueous",
            ],
            equation: "4 Au(s) + 8 CN[-](aq) + O2(g) + 2 H2O(l) -> 4 Au(CN)2[-](aq) + 4 OH[-](aq), historical NBS/1989-errata Elsner reference",
        },
        "copper-electrowinning" => ReactionCase {
            reaction: MineralReaction::CopperElectrowinning,
            reference: USGS_REFERENCE,
            packs: &[
                "copper-ii-ion-aqueous-usgs-b2131",
                "water-liquid-usgs-b2131",
                "copper-solid-usgs-b2131",
                "oxygen-gas-usgs-b2131",
                "hydrogen-ion-aqueous-usgs-b2131",
            ],
            species: &["Cu[2+](aq)", "H2O(l)", "Cu(s)", "O2(g)", "H+(aq)"],
            phase_axes: &[
                "source-phase-aqueous",
                "source-phase-liquid",
                "source-phase-solid",
                "source-phase-gas",
                "source-phase-aqueous",
            ],
            equation: "2 Cu[2+](aq) + 2 H2O(l) -> 2 Cu(s) + O2(g) + 4 H+(aq), oxygen-evolving electrowinning endpoint",
        },
        "pyrite-roasting" => ReactionCase {
            reaction: MineralReaction::PyriteRoastingToHematite,
            reference: REFERENCE,
            packs: &[
                "pyrite-solid-janaf-thermochemistry",
                "oxygen-gas-janaf-mining",
                "hematite-solid-janaf-mining",
                "sulfur-dioxide-gas-janaf-mining",
            ],
            species: &["FeS2(s)", "O2(g)", "Fe2O3(s)", "SO2(g)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-gas",
                "source-phase-solid",
                "source-phase-gas",
            ],
            equation: "4 FeS2(s) + 11 O2(g) -> 2 Fe2O3(s) + 8 SO2(g), hematite endpoint",
        },
        "pyrite-oxidation-ferrous" => ReactionCase {
            reaction: MineralReaction::PyriteOxidationFerrous,
            reference: USGS_REFERENCE,
            packs: &[
                "pyrite-solid-usgs-b2131-thermochemistry",
                "oxygen-gas-usgs-b2131",
                "water-liquid-usgs-b2131",
                "iron-ii-ion-aqueous-usgs-b2131",
                "sulfate-ion-aqueous-usgs-b2131",
                "hydrogen-ion-aqueous-usgs-b2131",
            ],
            species: &[
                "FeS2(s)",
                "O2(g)",
                "H2O(l)",
                "Fe[2+](aq)",
                "SO4[2-](aq)",
                "H+(aq)",
            ],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-gas",
                "source-phase-liquid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-aqueous",
            ],
            equation: "2 FeS2(s) + 7 O2(g) + 2 H2O(l) -> 2 Fe[2+](aq) + 4 SO4[2-](aq) + 4 H+(aq), ferrous acid-generation endpoint",
        },
        "ferrous-oxidation-ferric" => ReactionCase {
            reaction: MineralReaction::FerrousOxidationFerric,
            reference: USGS_REFERENCE,
            packs: &[
                "iron-ii-ion-aqueous-usgs-b2131",
                "oxygen-gas-usgs-b2131",
                "hydrogen-ion-aqueous-usgs-b2131",
                "iron-iii-ion-aqueous-usgs-b2131",
                "water-liquid-usgs-b2131",
            ],
            species: &["Fe[2+](aq)", "O2(g)", "H+(aq)", "Fe[3+](aq)", "H2O(l)"],
            phase_axes: &[
                "source-phase-aqueous",
                "source-phase-gas",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-liquid",
            ],
            equation: "4 Fe[2+](aq) + O2(g) + 4 H+(aq) -> 4 Fe[3+](aq) + 2 H2O(l), ferrous oxidation endpoint",
        },
        "ferric-pyrite-oxidation" => ReactionCase {
            reaction: MineralReaction::FerricPyriteOxidation,
            reference: USGS_REFERENCE,
            packs: &[
                "pyrite-solid-usgs-b2131-thermochemistry",
                "iron-iii-ion-aqueous-usgs-b2131",
                "water-liquid-usgs-b2131",
                "iron-ii-ion-aqueous-usgs-b2131",
                "sulfate-ion-aqueous-usgs-b2131",
                "hydrogen-ion-aqueous-usgs-b2131",
            ],
            species: &[
                "FeS2(s)",
                "Fe[3+](aq)",
                "H2O(l)",
                "Fe[2+](aq)",
                "SO4[2-](aq)",
                "H+(aq)",
            ],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-liquid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-aqueous",
            ],
            equation: "FeS2(s) + 14 Fe[3+](aq) + 8 H2O(l) -> 15 Fe[2+](aq) + 2 SO4[2-](aq) + 16 H+(aq), ferric attack on pyrite",
        },
        "ferric-hydrolysis-goethite" => ReactionCase {
            reaction: MineralReaction::FerricHydrolysisGoethite,
            reference: USGS_REFERENCE,
            packs: &[
                "iron-iii-ion-aqueous-usgs-b2131",
                "water-liquid-usgs-b2131",
                "goethite-solid-usgs-b2131-thermochemistry",
                "hydrogen-ion-aqueous-usgs-b2131",
            ],
            species: &["Fe[3+](aq)", "H2O(l)", "FeO(OH)(s)", "H+(aq)"],
            phase_axes: &[
                "source-phase-aqueous",
                "source-phase-liquid",
                "source-phase-solid",
                "source-phase-aqueous",
            ],
            equation: "Fe[3+](aq) + 2 H2O(l) -> FeO(OH)(goethite) + 3 H+(aq), ferric hydrolysis endpoint",
        },
        "ferric-goethite-precipitation-with-hydrated-lime" => ReactionCase {
            reaction: MineralReaction::FerricGoethitePrecipitationWithHydratedLime,
            reference: USGS_REFERENCE,
            packs: &[
                "iron-iii-ion-aqueous-usgs-b2131",
                "portlandite-solid-usgs-b2131-thermochemistry",
                "goethite-solid-usgs-b2131-thermochemistry",
                "calcium-ion-aqueous-usgs-b2131",
                "water-liquid-usgs-b2131",
            ],
            species: &[
                "Fe[3+](aq)",
                "Ca(OH)2(s)",
                "FeO(OH)(s)",
                "Ca[2+](aq)",
                "H2O(l)",
            ],
            phase_axes: &[
                "source-phase-aqueous",
                "source-phase-solid",
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-liquid",
            ],
            equation: "2 Fe[3+](aq) + 3 Ca(OH)2(portlandite) -> 2 FeO(OH)(goethite) + 3 Ca[2+](aq) + 2 H2O(l), ferric precipitation with hydrated lime",
        },
        "sulfur-dioxide-oxidation" => ReactionCase {
            reaction: MineralReaction::SulfurDioxideOxidation,
            reference: REFERENCE,
            packs: &[
                "sulfur-dioxide-gas-janaf-mining",
                "oxygen-gas-janaf-mining",
                "sulfur-trioxide-gas-janaf-mining",
            ],
            species: &["SO2(g)", "O2(g)", "SO3(g)"],
            phase_axes: &["source-phase-gas", "source-phase-gas", "source-phase-gas"],
            equation: "2 SO2(g) + O2(g) -> 2 SO3(g)",
        },
        "sulfur-trioxide-hydration" => ReactionCase {
            reaction: MineralReaction::SulfurTrioxideHydration,
            reference: REFERENCE,
            packs: &[
                "sulfur-trioxide-gas-janaf-mining",
                "water-liquid-janaf-mining",
                "sulfuric-acid-liquid-janaf-mining",
            ],
            species: &["SO3(g)", "H2O(l)", "H2SO4(l)"],
            phase_axes: &[
                "source-phase-gas",
                "source-phase-liquid",
                "source-phase-liquid",
            ],
            equation: "SO3(g) + H2O(l) -> H2SO4(l), pure liquid acid endpoint",
        },
        "pyrite-to-acid" => ReactionCase {
            reaction: MineralReaction::PyriteToHematiteAndSulfuricAcid,
            reference: REFERENCE,
            packs: &[
                "pyrite-solid-janaf-thermochemistry",
                "oxygen-gas-janaf-mining",
                "water-liquid-janaf-mining",
                "hematite-solid-janaf-mining",
                "sulfuric-acid-liquid-janaf-mining",
            ],
            species: &["FeS2(s)", "O2(g)", "H2O(l)", "Fe2O3(s)", "H2SO4(l)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-gas",
                "source-phase-liquid",
                "source-phase-solid",
                "source-phase-liquid",
            ],
            equation: "4 FeS2(s) + 15 O2(g) + 8 H2O(l) -> 2 Fe2O3(s) + 8 H2SO4(l), all sulfur assigned to pure liquid acid",
        },
        "calcite-bicarbonate-neutralization" => ReactionCase {
            reaction: MineralReaction::CalciteBicarbonateNeutralization,
            reference: USGS_REFERENCE,
            packs: &[
                "calcite-solid-usgs-b2131-thermochemistry",
                "hydrogen-ion-aqueous-usgs-b2131",
                "calcium-ion-aqueous-usgs-b2131",
                "bicarbonate-ion-aqueous-usgs-b2131",
            ],
            species: &["CaCO3(s)", "H+(aq)", "Ca[2+](aq)", "HCO3[-](aq)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-aqueous",
            ],
            equation: "CaCO3(calcite) + H+(aq) -> Ca[2+](aq) + HCO3[-](aq), bicarbonate endpoint",
        },
        "calcite-carbonic-acid-neutralization" => ReactionCase {
            reaction: MineralReaction::CalciteCarbonicAcidNeutralization,
            reference: USGS_REFERENCE,
            packs: &[
                "calcite-solid-usgs-b2131-thermochemistry",
                "hydrogen-ion-aqueous-usgs-b2131",
                "calcium-ion-aqueous-usgs-b2131",
                "carbonic-acid-aqueous-usgs-b2131",
            ],
            species: &["CaCO3(s)", "H+(aq)", "Ca[2+](aq)", "H2CO3(aq)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-aqueous",
            ],
            equation: "CaCO3(calcite) + 2 H+(aq) -> Ca[2+](aq) + H2CO3(aq), source-defined neutral carbonic endpoint",
        },
        "malachite-carbonic-acid-leach" => ReactionCase {
            reaction: MineralReaction::MalachiteCarbonicAcidLeach,
            reference: USGS_REFERENCE,
            packs: &[
                "malachite-solid-usgs-b2131-thermochemistry",
                "hydrogen-ion-aqueous-usgs-b2131",
                "copper-ii-ion-aqueous-usgs-b2131",
                "carbonic-acid-aqueous-usgs-b2131",
                "water-liquid-usgs-b2131",
            ],
            species: &[
                "Cu2CO3(OH)2(s)",
                "H+(aq)",
                "Cu[2+](aq)",
                "H2CO3(aq)",
                "H2O(l)",
            ],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-liquid",
            ],
            equation: "Cu2CO3(OH)2(malachite) + 4 H+(aq) -> 2 Cu[2+](aq) + H2CO3(aq) + 2 H2O(l), source-defined neutral carbonic endpoint",
        },
        "dolomite-bicarbonate-neutralization" => ReactionCase {
            reaction: MineralReaction::DolomiteBicarbonateNeutralization,
            reference: USGS_REFERENCE,
            packs: &[
                "dolomite-solid-usgs-b2131-thermochemistry",
                "hydrogen-ion-aqueous-usgs-b2131",
                "calcium-ion-aqueous-usgs-b2131",
                "magnesium-ion-aqueous-usgs-b2131",
                "bicarbonate-ion-aqueous-usgs-b2131",
            ],
            species: &[
                "CaMg(CO3)2(s)",
                "H+(aq)",
                "Ca[2+](aq)",
                "Mg[2+](aq)",
                "HCO3[-](aq)",
            ],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-aqueous",
            ],
            equation: "CaMg(CO3)2(dolomite) + 2 H+(aq) -> Ca[2+](aq) + Mg[2+](aq) + 2 HCO3[-](aq), bicarbonate endpoint",
        },
        "ferric-chalcopyrite-leach" => ReactionCase {
            reaction: MineralReaction::FerricChalcopyriteLeach,
            reference: USGS_REFERENCE,
            packs: &[
                "chalcopyrite-solid-usgs-b2131-thermochemistry",
                "iron-iii-ion-aqueous-usgs-b2131",
                "copper-ii-ion-aqueous-usgs-b2131",
                "iron-ii-ion-aqueous-usgs-b2131",
                "sulfur-orthorhombic-solid-usgs-b2131",
            ],
            species: &[
                "CuFeS2(s)",
                "Fe[3+](aq)",
                "Cu[2+](aq)",
                "Fe[2+](aq)",
                "S(s)",
            ],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-solid",
            ],
            equation: "CuFeS2(chalcopyrite) + 4 Fe[3+](aq) -> Cu[2+](aq) + 5 Fe[2+](aq) + 2 S(orthorhombic)",
        },
        "ferric-chalcocite-leach" => ReactionCase {
            reaction: MineralReaction::FerricChalcociteLeach,
            reference: USGS_REFERENCE,
            packs: &[
                "chalcocite-solid-usgs-b2131-thermochemistry",
                "iron-iii-ion-aqueous-usgs-b2131",
                "copper-ii-ion-aqueous-usgs-b2131",
                "iron-ii-ion-aqueous-usgs-b2131",
                "sulfur-orthorhombic-solid-usgs-b2131",
            ],
            species: &["Cu2S(s)", "Fe[3+](aq)", "Cu[2+](aq)", "Fe[2+](aq)", "S(s)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-solid",
            ],
            equation: "Cu2S(chalcocite) + 4 Fe[3+](aq) -> 2 Cu[2+](aq) + 4 Fe[2+](aq) + S(orthorhombic)",
        },
        "ferric-covellite-leach" => ReactionCase {
            reaction: MineralReaction::FerricCovelliteLeach,
            reference: USGS_REFERENCE,
            packs: &[
                "covellite-solid-usgs-b2131-thermochemistry",
                "iron-iii-ion-aqueous-usgs-b2131",
                "copper-ii-ion-aqueous-usgs-b2131",
                "iron-ii-ion-aqueous-usgs-b2131",
                "sulfur-orthorhombic-solid-usgs-b2131",
            ],
            species: &["CuS(s)", "Fe[3+](aq)", "Cu[2+](aq)", "Fe[2+](aq)", "S(s)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-solid",
            ],
            equation: "CuS(covellite) + 2 Fe[3+](aq) -> Cu[2+](aq) + 2 Fe[2+](aq) + S(orthorhombic)",
        },
        "ferric-bornite-leach" => ReactionCase {
            reaction: MineralReaction::FerricBorniteLeach,
            reference: USGS_REFERENCE,
            packs: &[
                "bornite-solid-usgs-b2131-thermochemistry",
                "iron-iii-ion-aqueous-usgs-b2131",
                "copper-ii-ion-aqueous-usgs-b2131",
                "iron-ii-ion-aqueous-usgs-b2131",
                "sulfur-orthorhombic-solid-usgs-b2131",
            ],
            species: &[
                "Cu5FeS4(s)",
                "Fe[3+](aq)",
                "Cu[2+](aq)",
                "Fe[2+](aq)",
                "S(s)",
            ],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-solid",
            ],
            equation: "Cu5FeS4(bornite) + 12 Fe[3+](aq) -> 5 Cu[2+](aq) + 13 Fe[2+](aq) + 4 S(orthorhombic)",
        },
        _ => return Err(format!("unknown reaction {name:?}; use --help for source order").into()),
    })
}

fn load(path: &Path, expected_id: &str) -> Result<NormalizedPack, Box<dyn Error>> {
    let bytes = std::fs::read(path)?;
    let decoded = NormalizedPack::from_bytes(&bytes)?;
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
    if pack.to_bytes() != bytes {
        return Err(format!("{} is not a canonical material pack", path.display()).into());
    }
    // This example explicitly associates known source profiles with phases;
    // arbitrary property packs do not authenticate a chemical formula.
    if pack.pack_id() != expected_id {
        return Err(format!("expected pack {expected_id}, received {}", pack.pack_id()).into());
    }
    Ok(pack)
}

fn reference_point(
    reference: &str,
    phase_axis: &str,
    pack_id: &str,
) -> Result<QueryPoint, MatDbError> {
    let mut point = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            298.15,
        )?
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )?;
    for (axis, value) in [
        ("source-standard-state", 1.0),
        (reference, 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
        (phase_axis, 1.0),
    ] {
        point = point.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    // The neutral USGS carbonic component retains its own formula convention;
    // it is neither an individual ion nor a pure molecular H2CO3 liquid.
    let state_axes: &[&str] = if pack_id == "carbonic-acid-aqueous-usgs-b2131" {
        &[
            "source-aqueous-unit-activity-infinite-dilution",
            "source-carbonic-reference-usgs-b2131",
        ]
    } else if phase_axis == "source-phase-aqueous" {
        &[
            "source-individual-ion",
            "source-aqueous-unit-activity-infinite-dilution",
            "source-ionic-reference-hydrogen-zero",
        ]
    } else {
        &["source-pure-species"]
    };
    for axis in state_axes {
        point = point.with_quantity(*axis, QuantitySpec::dimensional(Dims::NONE), 1.0)?;
    }
    if pack_id == "malachite-solid-usgs-b2131-thermochemistry" {
        point = point.with_quantity(
            "source-phase-malachite",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )?;
    }
    if pack_id == "copper-ii-oxide-solid-usgs-b2131-thermochemistry" {
        point = point.with_quantity(
            "source-phase-tenorite",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )?;
    }
    if pack_id == "calcite-solid-usgs-b2131-thermochemistry" {
        point = point.with_quantity(
            "source-phase-calcite",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )?;
    }
    if pack_id == "dolomite-solid-usgs-b2131-thermochemistry" {
        point = point.with_quantity(
            "source-phase-dolomite",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )?;
    }
    if pack_id == "goethite-solid-usgs-b2131-thermochemistry" {
        point = point.with_quantity(
            "source-phase-goethite",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )?;
    }
    if pack_id == "portlandite-solid-usgs-b2131-thermochemistry" {
        point = point.with_quantity(
            "source-phase-portlandite",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )?;
    }
    if matches!(
        pack_id,
        "pyrite-solid-janaf-thermochemistry" | "pyrite-solid-usgs-b2131-thermochemistry"
    ) {
        point = point.with_quantity(
            "source-phase-pyrite",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )?;
    }
    let mineral_axis = match pack_id {
        "chalcopyrite-solid-usgs-b2131-thermochemistry" => Some("source-phase-chalcopyrite"),
        "chalcocite-solid-usgs-b2131-thermochemistry" => Some("source-phase-chalcocite"),
        "covellite-solid-usgs-b2131-thermochemistry" => Some("source-phase-covellite"),
        "bornite-solid-usgs-b2131-thermochemistry" => Some("source-phase-bornite"),
        "sulfur-orthorhombic-solid-usgs-b2131" => Some("source-phase-sulfur-orthorhombic"),
        _ => None,
    };
    if let Some(axis) = mineral_axis {
        point = point.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), 1.0)?;
    }
    Ok(point)
}

fn malachite_inventory(args: &[String]) -> Result<(), Box<dyn Error>> {
    if args.len() != 5 {
        return Err("usage: mineral_processing --malachite-inventory MALACHITE_MOL PROTON_MOL WATER_MOL EXTENT_MOL".into());
    }
    let reaction = MineralReaction::MalachiteAcidLeach.build()?;
    let initial = [
        (
            SpeciesId::new("Cu2CO3(OH)2(s)")?,
            Amount::new(args[1].parse()?),
        ),
        (SpeciesId::new("H+(aq)")?, Amount::new(args[2].parse()?)),
        (SpeciesId::new("H2O(l)")?, Amount::new(args[3].parse()?)),
        (SpeciesId::new("Cu[2+](aq)")?, Amount::new(0.0)),
        (SpeciesId::new("CO2(aq)")?, Amount::new(0.0)),
    ];
    let extent = Amount::new(args[4].parse()?);
    let bound = reaction.limiting_extent(&initial[..2])?;
    let final_inventory = reaction.advance_inventory(&initial, extent)?;
    println!("Cu2CO3(OH)2(s) + 4 H+(aq) -> 2 Cu[2+](aq) + CO2(aq) + 3 H2O(l)");
    println!(
        "material-availability extent bound [mol]: {}",
        bound.value()
    );
    println!("caller-declared extent [mol]: {}", extent.value());
    for (species, remaining) in final_inventory {
        let before = initial
            .iter()
            .find(|(name, _)| *name == species)
            .ok_or("final species is absent from initial inventory")?
            .1;
        println!(
            "{} [mol]: {} -> {}",
            species.as_str(),
            before.value(),
            remaining.value()
        );
    }
    println!(
        "Extent is an input, not predicted recovery. Proton inventory is not a commercial acid dose."
    );
    Ok(())
}

struct SourceOptions<'a> {
    reaction_name: &'a str,
    arguments: &'a [String],
    equilibrium: bool,
    gibbs: bool,
    activities: Vec<(SpeciesId, Dimensionless)>,
}

fn source_options(args: &[String]) -> Result<SourceOptions<'_>, Box<dyn Error>> {
    let mut reaction_name = None;
    let mut equilibrium = false;
    let mut gibbs = false;
    let mut activities = Vec::new();
    let mut next = 0;
    while let Some(argument) = args.get(next) {
        match argument.as_str() {
            "--reaction" => {
                if reaction_name.is_some() {
                    return Err("--reaction may be supplied only once".into());
                }
                reaction_name = Some(
                    args.get(next + 1)
                        .ok_or("--reaction requires a reaction name")?
                        .as_str(),
                );
                next += 2;
            }
            "--equilibrium" => {
                equilibrium = true;
                next += 1;
            }
            "--gibbs" => {
                gibbs = true;
                next += 1;
            }
            "--activity" => {
                let input = args
                    .get(next + 1)
                    .ok_or("--activity requires 'SPECIES=VALUE'")?;
                let (species, value) = input
                    .split_once('=')
                    .ok_or("--activity requires 'SPECIES=VALUE'")?;
                let value: f64 = value.parse()?;
                if !value.is_finite() || value <= 0.0 {
                    return Err("activity must be finite and strictly positive".into());
                }
                activities.push((SpeciesId::new(species)?, Dimensionless::new(value)));
                gibbs = true;
                next += 2;
            }
            "--" => {
                next += 1;
                break;
            }
            unknown if unknown.starts_with("--") => {
                return Err(format!("unknown option {unknown:?}; use --help").into());
            }
            _ => break,
        }
    }
    Ok(SourceOptions {
        reaction_name: reaction_name.unwrap_or("lime-hydration"),
        arguments: &args[next..],
        equilibrium,
        gibbs,
        activities,
    })
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.as_slice() == ["--help"] {
        println!("{USAGE}");
        println!("Default without --reaction: lime-hydration. Source pack order:");
        for name in REACTIONS {
            let selected = case(name)?;
            println!("{name}: {}", selected.packs.join(" "));
            println!("  activity species: {}", selected.species.join(" "));
        }
        println!("Options precede the pack list; --activity implies --gibbs.");
        println!(
            "--equilibrium derives standard Gibbs energy from source enthalpy and entropy (Delta G = Delta H - T Delta S); it requires both properties for every species."
        );
        println!(
            "--gibbs requires a direct standard-molar-formation-gibbs-energy claim in every supplied pack; unsupported source sets refuse."
        );
        println!(
            "The two routes are reported separately when both are requested; independently tabulated source values can differ. Activities always use the direct formation-Gibbs route."
        );
        println!(
            "Supply one --activity 'SPECIES=VALUE' for every active species, including solids and water. Values are dimensionless relative activities under the retained standard states, not molarity."
        );
        println!(
            "The final mass is of the first declared species actually reacted; an aqueous-ion mass is not solution mass."
        );
        println!(
            "Activity Gibbs energy is dG/dxi at the supplied state, not the total change for a finite reaction extent. No equilibrium composition, reaction rate or conversion is predicted."
        );
        println!(
            "For copper-electrowinning, --gibbs also reports E = -Delta G/(nF), with two electrons per Cu and four per written equation. The signed reversible potential is not an operating voltage or an efficiency-scaled voltage."
        );
        return Ok(());
    }
    if args
        .first()
        .is_some_and(|arg| arg == "--malachite-inventory")
    {
        return malachite_inventory(&args);
    }
    let options = source_options(&args)?;
    let case = case(options.reaction_name)?;
    let args = options.arguments;
    if args.len() != case.packs.len() + 1 {
        return Err(format!(
            "expected {} source packs in this order, then reacted kg: {}",
            case.packs.len(),
            case.packs.join(" ")
        )
        .into());
    }
    let reacted_kg: f64 = args[case.packs.len()].parse()?;
    let packs = args
        .iter()
        .zip(case.packs)
        .map(|(path, expected)| load(Path::new(path), expected))
        .collect::<Result<Vec<_>, _>>()?;
    let species = case
        .species
        .iter()
        .map(|id| SpeciesId::new(*id))
        .collect::<Result<Vec<_>, _>>()?;
    let points = case
        .phase_axes
        .iter()
        .zip(case.packs)
        .map(|(phase, pack)| reference_point(case.reference, phase, pack))
        .collect::<Result<Vec<_>, _>>()?;
    let sources = (0..packs.len())
        .map(|i| FormationEnthalpySource {
            species: &species[i],
            pack: &packs[i],
            at: &points[i],
            reference_axis: case.reference,
        })
        .collect::<Vec<_>>();
    let reaction = case.reaction.build()?;
    let equilibrium = if options.equilibrium {
        Some(reaction.standard_equilibrium(&sources)?)
    } else {
        None
    };
    let heat = match &equilibrium {
        Some(standard) => standard.enthalpy().clone(),
        None => reaction.standard_enthalpy_change(&sources)?,
    };
    let mut molar_masses = Vec::with_capacity(packs.len());
    for (pack, point) in packs.iter().zip(&points) {
        let answer = pack.claims().query_typed(
            &PropertyKey::with_quantity(
                "molar-mass",
                QuantitySpec::dimensional(Dims([0, 1, 0, 0, 0, -1])),
            ),
            point,
            SelectionPolicy::SingleClaimOnly,
        )?;
        pack.claims().verify_receipt(&answer.receipt)?;
        molar_masses.push(MolarMass::new(answer.evidence.value.value));
    }
    let (gibbs, cell_potential) = if options.gibbs {
        let sources = (0..packs.len())
            .map(|i| FormationGibbsSource {
                species: &species[i],
                pack: &packs[i],
                at: &points[i],
                reference_axis: case.reference,
            })
            .collect::<Vec<_>>();
        if case.reaction == MineralReaction::CopperElectrowinning {
            let model = ElectrolysisReaction::new(
                reaction.clone(),
                species[2].clone(),
                molar_masses[2],
                2,
            )?;
            (None, Some(model.standard_potential(&sources)?))
        } else {
            (Some(reaction.standard_gibbs_change(&sources)?), None)
        }
    } else {
        (None, None)
    };
    let standard_gibbs = cell_potential
        .as_ref()
        .map(|standard| standard.gibbs())
        .or(gibbs.as_ref());
    let activity_potential = if options.activities.is_empty() {
        None
    } else {
        cell_potential
            .as_ref()
            .map(|standard| standard.at_activities(&options.activities))
            .transpose()?
    };
    let activity_gibbs = if options.activities.is_empty() || cell_potential.is_some() {
        None
    } else {
        Some(
            standard_gibbs
                .ok_or("activities require a standard Gibbs evaluation")?
                .at_activities(&options.activities)?,
        )
    };
    let activity_evaluation = activity_potential
        .as_ref()
        .map(|evaluation| evaluation.gibbs())
        .or(activity_gibbs.as_ref());
    let extent =
        reaction.extent_from_reacted_mass(&species[0], Mass::new(reacted_kg), molar_masses[0])?;
    let amounts = reaction.amounts(extent)?;
    let enthalpy_change = heat.at_extent(extent)?;
    let mut converted_amounts = Vec::with_capacity(amounts.len());
    let mut consumed_mass = 0.0;
    let mut formed_mass = 0.0;
    for amount in &amounts {
        let index = species
            .iter()
            .position(|id| id == &amount.term.species)
            .ok_or("reaction participant lacks an associated source")?;
        let mass = MassAmountBasis::from_amount(amount.amount, molar_masses[index])?;
        let direction = if amount.term.coefficient < 0 {
            consumed_mass += mass.mass().value();
            "consumed"
        } else {
            formed_mass += mass.mass().value();
            "formed"
        };
        converted_amounts.push((index, amount.amount, mass.mass(), direction));
    }
    if !consumed_mass.is_finite() || !formed_mass.is_finite() {
        return Err("summed source-molar-mass conversion is not representable".into());
    }

    println!("{}, 298.15 K, 100000 Pa", case.equation);
    if case.phase_axes[0] == "source-phase-aqueous" {
        println!(
            "reacted declared {} species mass [kg]: {reacted_kg} (not solution mass)",
            case.species[0]
        );
    } else {
        println!("reacted pure {} mass [kg]: {reacted_kg}", case.species[0]);
    }
    println!(
        "source reactant molar mass [kg/mol]: {}",
        molar_masses[0].value()
    );
    println!(
        "reaction extent [mol of written equation]: {}",
        extent.value()
    );
    for (index, amount, mass, direction) in converted_amounts {
        println!(
            "{} {direction}: {} mol, {} kg",
            case.species[index],
            amount.value(),
            mass.value()
        );
    }
    println!(
        "Source-molar-mass residual, formed minus consumed [kg]: {}",
        formed_mass - consumed_mass
    );
    println!(
        "Atoms and charge close exactly; independently rounded source molar masses are retained without adjustment."
    );
    println!(
        "standard reaction enthalpy [J/mol extent]: {}",
        heat.molar().value()
    );
    println!("standard enthalpy change [J]: {}", enthalpy_change.value());
    if let Some(standard) = equilibrium {
        println!(
            "caloric-derived standard-state route (--equilibrium): Delta G = Delta H - T Delta S"
        );
        println!(
            "standard reaction entropy [J/(mol extent K)]: {}",
            standard.molar_entropy().value()
        );
        println!(
            "standard reaction Gibbs energy [J/mol extent]: {}",
            standard.molar_gibbs().value()
        );
        println!(
            "ln K (dimensionless standard activities): {}",
            standard.log_equilibrium_constant().value()
        );
        println!(
            "log10 K (dimensionless standard activities): {}",
            standard.log_equilibrium_constant().value() / core::f64::consts::LN_10
        );
        println!(
            "These constants use the selected standard states; activities, pH and equilibrium composition require separate inputs and models."
        );
    }
    if let Some(standard) = standard_gibbs {
        println!("tabulated formation-Gibbs route (--gibbs):");
        println!(
            "standard reaction Gibbs energy [J/mol extent]: {}",
            standard.molar().value()
        );
        println!(
            "ln K (dimensionless thermodynamic equilibrium constant): {}",
            standard.log_equilibrium_constant().value()
        );
        println!("direct K representation: {:?}", standard.status());
        if let Some(value) = standard.equilibrium_constant() {
            println!("K [dimensionless]: {}", value.value());
        }
    }
    if let Some(standard) = &cell_potential {
        let model = standard.model();
        println!(
            "declared electrons per {} product: {}; per written equation: {}",
            model.product().as_str(),
            model.electrons_per_product(),
            model.electrons_per_equation()
        );
        println!(
            "source deposited-product molar mass [kg/mol]: {}",
            model.molar_mass().value()
        );
        println!(
            "signed standard reversible cell potential [V]: {}",
            standard.potential().value()
        );
    }
    if options.equilibrium && options.gibbs {
        println!(
            "The two routes retain independently rounded source quantities; their results are reported separately."
        );
    }
    if let Some(evaluation) = activity_evaluation {
        println!("activity diagnostics from tabulated formation-Gibbs sources:");
        for (species, activity) in &evaluation.activities {
            println!(
                "relative activity {} [dimensionless]: {}",
                species.as_str(),
                activity.value()
            );
        }
        println!("ln Q [dimensionless]: {}", evaluation.log_quotient.value());
        println!(
            "ln(Q/K) [dimensionless]: {}",
            evaluation.log_quotient_over_equilibrium.value()
        );
        println!(
            "activity reaction Gibbs energy dG/dxi [J/mol]: {}",
            evaluation.molar.value()
        );
        println!(
            "forward affinity [J/mol]: {}",
            evaluation.affinity().value()
        );
        println!(
            "Negative dG/dxi favors the written forward direction at these activities; this is not the total Gibbs change over the finite reacted extent."
        );
    }
    if let Some(evaluation) = &activity_potential {
        println!(
            "signed reversible cell potential at supplied activities [V]: {}",
            evaluation.potential().value()
        );
    }
    if cell_potential.is_some() {
        println!(
            "A negative signed potential identifies an externally driven written reaction. Reversing its sign gives the reversible opposing voltage, not the actual operating voltage; no efficiency factor is applied."
        );
        println!(
            "Overpotentials, electrical resistance, current efficiency and side reactions require additional models."
        );
    }
    for pack in &packs {
        println!("source pack: {} ({})", pack.pack_id(), pack.content_hash());
    }
    println!("Negative enthalpy is exothermic. Source confidence probabilities are unstated.");
    println!(
        "Reported source uncertainty magnitudes, where available, are not a propagated reaction interval."
    );
    println!(
        "This uses actually reacted declared chemical species; it does not predict conversion, plant temperature, kiln fuel duty or operating reagent doses."
    );
    if options.gibbs || options.equilibrium {
        println!(
            "Gibbs diagnostics retain the source reference states; no solution activity model, equilibrium composition or reaction rate is inferred."
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn arguments(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_owned()).collect()
    }

    #[test]
    fn source_options_preserve_old_defaults_and_named_reactions() {
        let args = arguments(&["cao", "water", "lime", "56.077"]);
        let options = source_options(&args).unwrap();
        assert_eq!(options.reaction_name, "lime-hydration");
        assert_eq!(options.arguments, args);
        assert!(!options.equilibrium);
        assert!(!options.gibbs);
        assert!(options.activities.is_empty());
        let args = arguments(&[
            "--reaction",
            "pyrite-oxidation-ferrous",
            "pyrite",
            "119.979",
        ]);
        let options = source_options(&args).unwrap();
        assert_eq!(options.reaction_name, "pyrite-oxidation-ferrous");
        assert_eq!(options.arguments, &args[2..]);
        assert!(!options.gibbs);
        let args = arguments(&["--equilibrium", "cao", "water", "lime", "56.077"]);
        let options = source_options(&args).unwrap();
        assert_eq!(options.reaction_name, "lime-hydration");
        assert_eq!(options.arguments, &args[1..]);
        assert!(options.equilibrium);
        assert!(!options.gibbs);
        let args = arguments(&[
            "--equilibrium",
            "--reaction",
            "ferrous-oxidation-ferric",
            "--gibbs",
            "iron",
            "55.847",
        ]);
        let options = source_options(&args).unwrap();
        assert_eq!(options.reaction_name, "ferrous-oxidation-ferric");
        assert_eq!(options.arguments, &args[4..]);
        assert!(options.equilibrium);
        assert!(options.gibbs);
    }

    #[test]
    fn explicit_activities_imply_gibbs_without_inventing_missing_values() {
        let args = arguments(&[
            "--activity",
            "Fe[2+](aq)=0.01",
            "--reaction",
            "ferrous-oxidation-ferric",
            "--activity",
            "H+(aq)=1e-3",
            "--gibbs",
            "--",
            "iron",
            "55.847",
        ]);
        let options = source_options(&args).unwrap();
        assert!(options.gibbs);
        assert_eq!(options.reaction_name, "ferrous-oxidation-ferric");
        assert_eq!(options.activities.len(), 2);
        assert_eq!(options.activities[0].0.as_str(), "Fe[2+](aq)");
        assert_eq!(options.activities[0].1.value(), 0.01);
        assert_eq!(options.arguments, &args[8..]);
    }

    #[test]
    fn malformed_options_and_invalid_activities_refuse() {
        for args in [
            vec!["--reaction"],
            vec![
                "--reaction",
                "lime-hydration",
                "--reaction",
                "lime-hydration",
            ],
            vec!["--activity"],
            vec!["--activity", "H+(aq)"],
            vec!["--activity", "H+(aq)=0"],
            vec!["--activity", "H+(aq)=-1"],
            vec!["--activity", "H+(aq)=NaN"],
            vec!["--activity", "H+(aq)=inf"],
            vec!["--unknown"],
        ] {
            assert!(source_options(&arguments(&args)).is_err());
        }
    }
}
