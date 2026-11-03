//! Reference mineral/reagent equations with explicit final species.
//!
//! Pyrite endpoints follow USGS Open-File Report 99-348:
//! <https://pubs.usgs.gov/of/1999/0348/report.pdf>. Its equations 3 and 4 give
//! ferrous oxidation and ferric attack on pyrite, respectively. Ferrous
//! oxidation is multiplied by four here to retain exact integer coefficients.
//! Combining seven such oxidations with two ferric-pyrite attacks recovers
//! the two-pyrite ferrous endpoint; intermediate iron is accounted for rather
//! than treated as a free oxidant. No reaction rate or cycle turnover is inferred.
//! Carbonate neutralizations
//! combine mineral dissolution with CO3(2-) + 2 H+ -> CO2(aq) + H2O from
//! the USGS PHREEQC database:
//! <https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc.v1/html/phqc_66.htm>.
//! Distinct bicarbonate endpoints follow USGS WRIR 95-4097, equations 4-5:
//! <https://pubs.usgs.gov/wri/1995/4097/report.pdf>. A carbonate's acid capacity
//! is therefore endpoint-dependent, not a universal measured ore property.
//! Lime hydration is CaO + H2O -> Ca(OH)2; source-bound formation enthalpies
//! are supplied separately. Hydrated-lime neutralization follows the EPA acid
//! mine drainage design manual:
//! <https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=20007H0I.TXT>.
//! Ferric hydrolysis to goethite follows USGS Open-File Report 2005-1148:
//! <https://pubs.usgs.gov/of/2005/1148/of2005-1148.pdf>. Combining two
//! hydrolysis extents with three hydrated-lime neutralizations gives the
//! declared goethite precipitation endpoint. FeO(OH)(s) here means goethite,
//! not a generic Fe(OH)3 precipitate or a prediction of the stable solid.
//! Malachite's net ionic acid balance is equation 10 in
//! <https://doi.org/10.3390/met10060833>; the pathway here retains dissolved CO2.
//! Separate carbonic-acid endpoints retain the literal un-ionized H2CO3(aq)
//! species in USGS Bulletin 2131, printed pages 15 and 25:
//! <https://pubs.usgs.gov/bul/2131/report.pdf>. They are atom/charge-balanced
//! source-reference reactions, not a prediction of aqueous speciation or a
//! substitution of H2CO3 formation data for a CO2(aq) property.
//! Net ionic Elsner cyanidation follows
//! <https://doi.org/10.1038/s41598-019-43383-4> and the EPA 1994 report
//! <https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=2000EF6J.TXT>.
//! The dissolved-oxygen Elsner endpoint keeps O2(aq) distinct from O2(g).
//! Aqueous HCN dissociation and water dissociation retain the molecular/ionic
//! states in the NBS tables, printed pages 2-37, 2-38 and 2-89, and the
//! standard-state definitions in sections 5.5-5.8:
//! <https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf> and the corrected NIST dataset
//! <https://doi.org/10.18434/M32124>. Molecular HCN uses the ao row; the ai row
//! already represents H+ and CN-. These equations can be composed, but their
//! presence alone does not solve coupled speciation or supply gas transfer.
//! Magnesite calcination and magnetite oxidation retain the crystalline/gas
//! endpoints of JANAF C-079, Mg-008, C-095, Fe-032, Fe-030 and O-029:
//! <https://janaf.nist.gov/>. Their standard reaction heats do not determine
//! a decomposition temperature, kiln duty, oxidation rate or phase stability.
//! Copper oxide reduction and oxidation use explicit CuO(cr), Cu2O(cr),
//! Cu(cr), CO(g), CO2(g) and O2(g) endpoints from the same JANAF reference.
//! A carbon-monoxide reduction equation does not select an industrial reductant
//! mixture, operating atmosphere, temperature, rate or metal recovery.
//! CuO acid leaching instead retains tenorite, H+(aq), Cu(II)(aq) and liquid
//! water from USGS Bulletin 2131, printed pages 6, 7 and 16. It is a declared
//! conserved ionic endpoint; the material references do not supply a rate.
//! Dry pyrite roasting and sulfur-to-acid balances use the explicit JANAF
//! Fe-028, Fe-030, O-029, O-034, O-058, H-063 and H-075 phases. These are
//! atom-balanced aggregate endpoints, not elementary reaction mechanisms.
//! SO3(g) + H2O(l) -> H2SO4(l) selects pure liquid acid; it does not model
//! absorption into aqueous acid, oleum, acid mist or heat of dilution.
//! These equations conserve atoms and charge; they do not select pH/Eh,
//! dissolved speciation, precipitates, gas transfer, rates or conversion.
//! Copper electrowinning combines Cu(II) reduction and water oxidation with
//! an inert oxygen-evolving anode; see Beukes and Badenhorst (2009):
//! <https://www.saimm.co.za/Journal/v109n06p343.pdf>. Its electron transfer is
//! four electrons per written equation; side reactions require separate models.
//! Ferric copper-sulfide leaching retains dissolved Cu(II), Fe(II) and elemental
//! sulfur endpoints. Chalcocite/covellite equations follow Ferron (2003),
//! *Leaching of Secondary Copper Minerals Using Regenerated Ferric Sulphate*,
//! equations 1-2, <https://www.sgs.com/-/media/sgscorp/documents/corporate/brochures/sgs-min-tp-2003-14-leaching-of-copper-using-ferric-sulphate.cdn.en-CA.pdf>.
//! Chalcopyrite/bornite are the net ionic forms of the oxygen-leach balances
//! after subtracting one/three ferrous-regeneration extents in Dreisinger et al.,
//! *The Competitive Position of the Mt Gordon Copper Process in the Copper
//! Industry*, ALTA Copper 2002, section 2.1, p.3 (PDF p.200):
//! <https://d3e2i5nuh73s15.cloudfront.net/wp-content/uploads/2015/04/ALTA-Copper-2002-Proceedings.pdf>.
//! These are complete-dissolution accounting endpoints, not elementary steps:
//! intermediate sulfides, passivation, sulfate formation and incomplete
//! extraction remain separate pathways. `S(s)` selects orthorhombic sulfur,
//! counted per mole of S atoms, not per mole of S8. This reference-state choice
//! does not predict the form of sulfur produced by an actual leach.

use super::{ConservedReaction, ReactionBalanceError};
use crate::{
    ChargeVector, ElementId, ElementalMatrix, ReactionId, SpeciesId, StoichiometricMatrix,
};

/// A declared endpoint for elementary mineral-processing bookkeeping.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MineralReaction {
    /// CaO(s) + H2O(l) -> Ca(OH)2(s).
    LimeHydration,
    /// CaCO3(s) + 2 H+ -> Ca2+ + CO2(aq) + H2O(l).
    /// Acidic dissolved-CO2 endpoint; bicarbonate buffering is a different equation.
    CalciteAcidNeutralization,
    /// CaCO3(s) + 2 H+ -> Ca2+ + H2CO3(aq).
    /// Literal un-ionized carbonic-acid endpoint for the USGS B2131 reference.
    CalciteCarbonicAcidNeutralization,
    /// CaMg(CO3)2(s) + 4 H+ -> Ca2+ + Mg2+ + 2 CO2(aq) + 2 H2O(l).
    DolomiteAcidNeutralization,
    /// CaCO3(s) + H+ -> Ca2+ + HCO3-; bicarbonate remains in solution.
    CalciteBicarbonateNeutralization,
    /// CaMg(CO3)2(s) + 2 H+ -> Ca2+ + Mg2+ + 2 HCO3-.
    DolomiteBicarbonateNeutralization,
    /// Ca(OH)2(s) + 2 H+ -> Ca2+ + 2 H2O(l).
    HydratedLimeNeutralization,
    /// Cu2CO3(OH)2(s) + 4 H+ -> 2 Cu2+ + CO2(aq) + 3 H2O(l).
    /// Stoichiometric dissolved-copper endpoint, not an ore-recovery prediction.
    MalachiteAcidLeach,
    /// Cu2CO3(OH)2(s) + 4 H+ -> 2 Cu2+ + H2CO3(aq) + 2 H2O(l).
    /// The source-defined H2CO3(aq) endpoint remains distinct from CO2(aq).
    MalachiteCarbonicAcidLeach,
    /// 4 Au(s) + 8 CN- + O2(g) + 2 H2O(l) -> 4 Au(CN)2- + 4 OH-.
    /// Net ionic Elsner balance; gas supply does not imply a dissolution rate.
    GoldCyanidation,
    /// 4 Au(s) + 8 CN- + O2(aq) + 2 H2O(l) -> 4 Au(CN)2- + 4 OH-.
    /// Dissolved oxygen is a finite aqueous stock, distinct from a gas reservoir.
    GoldCyanidationDissolvedOxygen,
    /// HCN(aq) -> H+(aq) + CN[-](aq).
    /// HCN is the neutral undissociated solute, not the ionized source component.
    HydrogenCyanideDissociation,
    /// H2O(l) -> H+(aq) + OH[-](aq).
    /// Standard water dissociation; activity laws and equilibrium remain explicit.
    WaterDissociation,
    /// 2 FeS2(s) + 7 O2(g) + 2 H2O(l) -> 2 Fe2+ + 4 SO4(2-) + 4 H+.
    PyriteOxidationFerrous,
    /// 4 Fe2+(aq) + O2(g) + 4 H+(aq) -> 4 Fe3+(aq) + 2 H2O(l).
    /// Oxidation consumes one proton per ferrous ion; rate remains an input.
    FerrousOxidationFerric,
    /// FeS2(s) + 14 Fe3+(aq) + 8 H2O(l) -> 15 Fe2+(aq) + 2 SO4(2-)(aq) + 16 H+(aq).
    /// Ferric iron is an explicit consumed inventory, not an unlimited catalyst.
    FerricPyriteOxidation,
    /// Fe3+(aq) + 2 H2O(l) -> FeO(OH)(s, goethite) + 3 H+(aq).
    /// Phase-specific hydrolysis releases three protons per ferric ion.
    FerricHydrolysisGoethite,
    /// 2 Fe3+(aq) + 3 Ca(OH)2(s) -> 2 FeO(OH)(s, goethite) + 3 Ca2+(aq) + 2 H2O(l).
    /// Combined hydrolysis and portlandite neutralization; extent is explicit.
    FerricGoethitePrecipitationWithHydratedLime,
    /// 4 FeS2(s) + 15 O2(g) + 14 H2O(l) -> 4 Fe(OH)3(s) + 8 SO4(2-) + 16 H+.
    /// Ideal ferric-hydroxide endpoint; no claim that a particular solid precipitates.
    PyriteOxidationFerricHydroxide,
    /// MgCO3(s) -> MgO(s) + CO2(g); crystalline carbonate/oxide endpoints.
    MagnesiteCalcination,
    /// 4 Fe3O4(s) + O2(g) -> 6 Fe2O3(s), with hematite as the product phase.
    /// One mole of reaction extent consumes four moles of magnetite.
    MagnetiteOxidationToHematite,
    /// 2 Cu2O(s) + O2(g) -> 4 CuO(s).
    /// One mole of reaction extent consumes two moles of cuprous oxide.
    CuprousOxideOxidation,
    /// CuO(s) + CO(g) -> Cu(s) + CO2(g).
    CupricOxideCarbonMonoxideReduction,
    /// CuO(s) + 2 H+(aq) -> Cu2+(aq) + H2O(l).
    /// Tenorite acid-leaching endpoint; two proton equivalents per copper.
    CupricOxideAcidLeach,
    /// Cu2O(s) + CO(g) -> 2 Cu(s) + CO2(g).
    CuprousOxideCarbonMonoxideReduction,
    /// 4 FeS2(s) + 11 O2(g) -> 2 Fe2O3(s) + 8 SO2(g).
    /// Dry hematite/SO2 roasting endpoint, distinct from aqueous oxidation.
    PyriteRoastingToHematite,
    /// 2 SO2(g) + O2(g) -> 2 SO3(g).
    /// Aggregate gas-phase balance; no catalytic conversion or rate is inferred.
    SulfurDioxideOxidation,
    /// SO3(g) + H2O(l) -> H2SO4(l), with pure liquid acid as endpoint.
    /// No aqueous dilution, absorber design or operating-dose claim.
    SulfurTrioxideHydration,
    /// 4 FeS2(s) + 15 O2(g) + 8 H2O(l) -> 2 Fe2O3(s) + 8 H2SO4(l).
    /// Aggregate of roasting, SO2 oxidation and hydration with all sulfur
    /// assigned to pure liquid acid. Actual sulfur capture remains a caller input.
    PyriteToHematiteAndSulfuricAcid,
    /// 2 Cu2+(aq) + 2 H2O(l) -> 2 Cu(s) + O2(g) + 4 H+(aq).
    /// Net productive pathway with an oxygen-evolving anode: four electrons
    /// transferred per equation, or two per copper atom deposited.
    CopperElectrowinning,
    /// CuFeS2(s) + 4 Fe3+(aq) -> Cu2+(aq) + 5 Fe2+(aq) + 2 S(s).
    /// Chalcopyrite dissolution to Cu(II) and orthorhombic sulfur; no sulfate.
    FerricChalcopyriteLeach,
    /// Cu2S(s) + 4 Fe3+(aq) -> 2 Cu2+(aq) + 4 Fe2+(aq) + S(s).
    /// Complete chalcocite endpoint, not the distinct first-stage CuS residue.
    FerricChalcociteLeach,
    /// CuS(s) + 2 Fe3+(aq) -> Cu2+(aq) + 2 Fe2+(aq) + S(s).
    /// Covellite dissolution; source properties must name this mineral phase.
    FerricCovelliteLeach,
    /// Cu5FeS4(s) + 12 Fe3+(aq) -> 5 Cu2+(aq) + 13 Fe2+(aq) + 4 S(s).
    /// Complete bornite endpoint; its one Fe atom remains as additional Fe(II).
    FerricBorniteLeach,
}

// One species: stable phase/charge id, element counts, charge, signed coefficient.
type Participant = (&'static str, &'static [i128], i128, i128);

impl MineralReaction {
    /// Build through the existing exact atom/charge checker. An extent is per
    /// integer equation as written, not implicitly per mole of mineral.
    /// Species ids include phase; Fe2O3(s) here denotes the hematite endpoint.
    pub fn build(self) -> Result<ConservedReaction, ReactionBalanceError> {
        let (name, elements, participants): (&str, &[&str], &[Participant]) = match self {
            Self::LimeHydration => (
                "lime-hydration",
                &["Ca", "H", "O"],
                &[
                    ("CaO(s)", &[1, 0, 1], 0, -1),
                    ("H2O(l)", &[0, 2, 1], 0, -1),
                    ("Ca(OH)2(s)", &[1, 2, 2], 0, 1),
                ],
            ),
            Self::CalciteAcidNeutralization => (
                "calcite-acid-co2",
                &["Ca", "C", "H", "O"],
                &[
                    ("CaCO3(s)", &[1, 1, 0, 3], 0, -1),
                    ("H+(aq)", &[0, 0, 1, 0], 1, -2),
                    ("Ca[2+](aq)", &[1, 0, 0, 0], 2, 1),
                    ("CO2(aq)", &[0, 1, 0, 2], 0, 1),
                    ("H2O(l)", &[0, 0, 2, 1], 0, 1),
                ],
            ),
            Self::CalciteCarbonicAcidNeutralization => (
                "calcite-acid-carbonic",
                &["Ca", "C", "H", "O"],
                &[
                    ("CaCO3(s)", &[1, 1, 0, 3], 0, -1),
                    ("H+(aq)", &[0, 0, 1, 0], 1, -2),
                    ("Ca[2+](aq)", &[1, 0, 0, 0], 2, 1),
                    ("H2CO3(aq)", &[0, 1, 2, 3], 0, 1),
                ],
            ),
            Self::DolomiteAcidNeutralization => (
                "dolomite-acid-co2",
                &["Ca", "Mg", "C", "H", "O"],
                &[
                    ("CaMg(CO3)2(s)", &[1, 1, 2, 0, 6], 0, -1),
                    ("H+(aq)", &[0, 0, 0, 1, 0], 1, -4),
                    ("Ca[2+](aq)", &[1, 0, 0, 0, 0], 2, 1),
                    ("Mg[2+](aq)", &[0, 1, 0, 0, 0], 2, 1),
                    ("CO2(aq)", &[0, 0, 1, 0, 2], 0, 2),
                    ("H2O(l)", &[0, 0, 0, 2, 1], 0, 2),
                ],
            ),
            Self::CalciteBicarbonateNeutralization => (
                "calcite-acid-bicarbonate",
                &["Ca", "C", "H", "O"],
                &[
                    ("CaCO3(s)", &[1, 1, 0, 3], 0, -1),
                    ("H+(aq)", &[0, 0, 1, 0], 1, -1),
                    ("Ca[2+](aq)", &[1, 0, 0, 0], 2, 1),
                    ("HCO3[-](aq)", &[0, 1, 1, 3], -1, 1),
                ],
            ),
            Self::DolomiteBicarbonateNeutralization => (
                "dolomite-acid-bicarbonate",
                &["Ca", "Mg", "C", "H", "O"],
                &[
                    ("CaMg(CO3)2(s)", &[1, 1, 2, 0, 6], 0, -1),
                    ("H+(aq)", &[0, 0, 0, 1, 0], 1, -2),
                    ("Ca[2+](aq)", &[1, 0, 0, 0, 0], 2, 1),
                    ("Mg[2+](aq)", &[0, 1, 0, 0, 0], 2, 1),
                    ("HCO3[-](aq)", &[0, 0, 1, 1, 3], -1, 2),
                ],
            ),
            Self::HydratedLimeNeutralization => (
                "hydrated-lime-neutralization",
                &["Ca", "H", "O"],
                &[
                    ("Ca(OH)2(s)", &[1, 2, 2], 0, -1),
                    ("H+(aq)", &[0, 1, 0], 1, -2),
                    ("Ca[2+](aq)", &[1, 0, 0], 2, 1),
                    ("H2O(l)", &[0, 2, 1], 0, 2),
                ],
            ),
            Self::MalachiteAcidLeach => (
                "malachite-acid-leach",
                &["Cu", "C", "H", "O"],
                &[
                    ("Cu2CO3(OH)2(s)", &[2, 1, 2, 5], 0, -1),
                    ("H+(aq)", &[0, 0, 1, 0], 1, -4),
                    ("Cu[2+](aq)", &[1, 0, 0, 0], 2, 2),
                    ("CO2(aq)", &[0, 1, 0, 2], 0, 1),
                    ("H2O(l)", &[0, 0, 2, 1], 0, 3),
                ],
            ),
            Self::MalachiteCarbonicAcidLeach => (
                "malachite-acid-carbonic",
                &["Cu", "C", "H", "O"],
                &[
                    ("Cu2CO3(OH)2(s)", &[2, 1, 2, 5], 0, -1),
                    ("H+(aq)", &[0, 0, 1, 0], 1, -4),
                    ("Cu[2+](aq)", &[1, 0, 0, 0], 2, 2),
                    ("H2CO3(aq)", &[0, 1, 2, 3], 0, 1),
                    ("H2O(l)", &[0, 0, 2, 1], 0, 2),
                ],
            ),
            Self::GoldCyanidation => (
                "gold-cyanidation-elsner",
                &["Au", "C", "N", "O", "H"],
                &[
                    ("Au(s)", &[1, 0, 0, 0, 0], 0, -4),
                    ("CN[-](aq)", &[0, 1, 1, 0, 0], -1, -8),
                    ("O2(g)", &[0, 0, 0, 2, 0], 0, -1),
                    ("H2O(l)", &[0, 0, 0, 1, 2], 0, -2),
                    ("Au(CN)2[-](aq)", &[1, 2, 2, 0, 0], -1, 4),
                    ("OH[-](aq)", &[0, 0, 0, 1, 1], -1, 4),
                ],
            ),
            Self::GoldCyanidationDissolvedOxygen => (
                "gold-cyanidation-elsner-dissolved-oxygen",
                &["Au", "C", "N", "O", "H"],
                &[
                    ("Au(s)", &[1, 0, 0, 0, 0], 0, -4),
                    ("CN[-](aq)", &[0, 1, 1, 0, 0], -1, -8),
                    ("O2(aq)", &[0, 0, 0, 2, 0], 0, -1),
                    ("H2O(l)", &[0, 0, 0, 1, 2], 0, -2),
                    ("Au(CN)2[-](aq)", &[1, 2, 2, 0, 0], -1, 4),
                    ("OH[-](aq)", &[0, 0, 0, 1, 1], -1, 4),
                ],
            ),
            Self::HydrogenCyanideDissociation => (
                "hydrogen-cyanide-dissociation",
                &["H", "C", "N"],
                &[
                    ("HCN(aq)", &[1, 1, 1], 0, -1),
                    ("H+(aq)", &[1, 0, 0], 1, 1),
                    ("CN[-](aq)", &[0, 1, 1], -1, 1),
                ],
            ),
            Self::WaterDissociation => (
                "water-dissociation",
                &["H", "O"],
                &[
                    ("H2O(l)", &[2, 1], 0, -1),
                    ("H+(aq)", &[1, 0], 1, 1),
                    ("OH[-](aq)", &[1, 1], -1, 1),
                ],
            ),
            Self::PyriteOxidationFerrous => (
                "pyrite-oxidation-ferrous",
                &["Fe", "S", "O", "H"],
                &[
                    ("FeS2(s)", &[1, 2, 0, 0], 0, -2),
                    ("O2(g)", &[0, 0, 2, 0], 0, -7),
                    ("H2O(l)", &[0, 0, 1, 2], 0, -2),
                    ("Fe[2+](aq)", &[1, 0, 0, 0], 2, 2),
                    ("SO4[2-](aq)", &[0, 1, 4, 0], -2, 4),
                    ("H+(aq)", &[0, 0, 0, 1], 1, 4),
                ],
            ),
            Self::FerrousOxidationFerric => (
                "ferrous-oxidation-ferric",
                &["Fe", "O", "H"],
                &[
                    ("Fe[2+](aq)", &[1, 0, 0], 2, -4),
                    ("O2(g)", &[0, 2, 0], 0, -1),
                    ("H+(aq)", &[0, 0, 1], 1, -4),
                    ("Fe[3+](aq)", &[1, 0, 0], 3, 4),
                    ("H2O(l)", &[0, 1, 2], 0, 2),
                ],
            ),
            Self::FerricPyriteOxidation => (
                "ferric-pyrite-oxidation",
                &["Fe", "S", "O", "H"],
                &[
                    ("FeS2(s)", &[1, 2, 0, 0], 0, -1),
                    ("Fe[3+](aq)", &[1, 0, 0, 0], 3, -14),
                    ("H2O(l)", &[0, 0, 1, 2], 0, -8),
                    ("Fe[2+](aq)", &[1, 0, 0, 0], 2, 15),
                    ("SO4[2-](aq)", &[0, 1, 4, 0], -2, 2),
                    ("H+(aq)", &[0, 0, 0, 1], 1, 16),
                ],
            ),
            Self::FerricHydrolysisGoethite => (
                "ferric-hydrolysis-goethite",
                &["Fe", "O", "H"],
                &[
                    ("Fe[3+](aq)", &[1, 0, 0], 3, -1),
                    ("H2O(l)", &[0, 1, 2], 0, -2),
                    ("FeO(OH)(s)", &[1, 2, 1], 0, 1),
                    ("H+(aq)", &[0, 0, 1], 1, 3),
                ],
            ),
            Self::FerricGoethitePrecipitationWithHydratedLime => (
                "ferric-goethite-precipitation-with-hydrated-lime",
                &["Fe", "Ca", "O", "H"],
                &[
                    ("Fe[3+](aq)", &[1, 0, 0, 0], 3, -2),
                    ("Ca(OH)2(s)", &[0, 1, 2, 2], 0, -3),
                    ("FeO(OH)(s)", &[1, 0, 2, 1], 0, 2),
                    ("Ca[2+](aq)", &[0, 1, 0, 0], 2, 3),
                    ("H2O(l)", &[0, 0, 1, 2], 0, 2),
                ],
            ),
            Self::PyriteOxidationFerricHydroxide => (
                "pyrite-oxidation-ferric-hydroxide",
                &["Fe", "S", "O", "H"],
                &[
                    ("FeS2(s)", &[1, 2, 0, 0], 0, -4),
                    ("O2(g)", &[0, 0, 2, 0], 0, -15),
                    ("H2O(l)", &[0, 0, 1, 2], 0, -14),
                    ("Fe(OH)3(s)", &[1, 0, 3, 3], 0, 4),
                    ("SO4[2-](aq)", &[0, 1, 4, 0], -2, 8),
                    ("H+(aq)", &[0, 0, 0, 1], 1, 16),
                ],
            ),
            Self::MagnesiteCalcination => (
                "magnesite-calcination",
                &["Mg", "C", "O"],
                &[
                    ("MgCO3(s)", &[1, 1, 3], 0, -1),
                    ("MgO(s)", &[1, 0, 1], 0, 1),
                    ("CO2(g)", &[0, 1, 2], 0, 1),
                ],
            ),
            Self::MagnetiteOxidationToHematite => (
                "magnetite-oxidation-to-hematite",
                &["Fe", "O"],
                &[
                    ("Fe3O4(s)", &[3, 4], 0, -4),
                    ("O2(g)", &[0, 2], 0, -1),
                    ("Fe2O3(s)", &[2, 3], 0, 6),
                ],
            ),
            Self::CuprousOxideOxidation => (
                "cuprous-oxide-oxidation",
                &["Cu", "O"],
                &[
                    ("Cu2O(s)", &[2, 1], 0, -2),
                    ("O2(g)", &[0, 2], 0, -1),
                    ("CuO(s)", &[1, 1], 0, 4),
                ],
            ),
            Self::CupricOxideCarbonMonoxideReduction => (
                "cupric-oxide-carbon-monoxide-reduction",
                &["Cu", "C", "O"],
                &[
                    ("CuO(s)", &[1, 0, 1], 0, -1),
                    ("CO(g)", &[0, 1, 1], 0, -1),
                    ("Cu(s)", &[1, 0, 0], 0, 1),
                    ("CO2(g)", &[0, 1, 2], 0, 1),
                ],
            ),
            Self::CupricOxideAcidLeach => (
                "cupric-oxide-acid-leach",
                &["Cu", "H", "O"],
                &[
                    ("CuO(s)", &[1, 0, 1], 0, -1),
                    ("H+(aq)", &[0, 1, 0], 1, -2),
                    ("Cu[2+](aq)", &[1, 0, 0], 2, 1),
                    ("H2O(l)", &[0, 2, 1], 0, 1),
                ],
            ),
            Self::CuprousOxideCarbonMonoxideReduction => (
                "cuprous-oxide-carbon-monoxide-reduction",
                &["Cu", "C", "O"],
                &[
                    ("Cu2O(s)", &[2, 0, 1], 0, -1),
                    ("CO(g)", &[0, 1, 1], 0, -1),
                    ("Cu(s)", &[1, 0, 0], 0, 2),
                    ("CO2(g)", &[0, 1, 2], 0, 1),
                ],
            ),
            Self::PyriteRoastingToHematite => (
                "pyrite-roasting-to-hematite",
                &["Fe", "S", "O"],
                &[
                    ("FeS2(s)", &[1, 2, 0], 0, -4),
                    ("O2(g)", &[0, 0, 2], 0, -11),
                    ("Fe2O3(s)", &[2, 0, 3], 0, 2),
                    ("SO2(g)", &[0, 1, 2], 0, 8),
                ],
            ),
            Self::SulfurDioxideOxidation => (
                "sulfur-dioxide-oxidation",
                &["S", "O"],
                &[
                    ("SO2(g)", &[1, 2], 0, -2),
                    ("O2(g)", &[0, 2], 0, -1),
                    ("SO3(g)", &[1, 3], 0, 2),
                ],
            ),
            Self::SulfurTrioxideHydration => (
                "sulfur-trioxide-hydration",
                &["S", "O", "H"],
                &[
                    ("SO3(g)", &[1, 3, 0], 0, -1),
                    ("H2O(l)", &[0, 1, 2], 0, -1),
                    ("H2SO4(l)", &[1, 4, 2], 0, 1),
                ],
            ),
            Self::PyriteToHematiteAndSulfuricAcid => (
                "pyrite-to-hematite-and-sulfuric-acid",
                &["Fe", "S", "O", "H"],
                &[
                    ("FeS2(s)", &[1, 2, 0, 0], 0, -4),
                    ("O2(g)", &[0, 0, 2, 0], 0, -15),
                    ("H2O(l)", &[0, 0, 1, 2], 0, -8),
                    ("Fe2O3(s)", &[2, 0, 3, 0], 0, 2),
                    ("H2SO4(l)", &[0, 1, 4, 2], 0, 8),
                ],
            ),
            Self::CopperElectrowinning => (
                "copper-electrowinning",
                &["Cu", "H", "O"],
                &[
                    ("Cu[2+](aq)", &[1, 0, 0], 2, -2),
                    ("H2O(l)", &[0, 2, 1], 0, -2),
                    ("Cu(s)", &[1, 0, 0], 0, 2),
                    ("O2(g)", &[0, 0, 2], 0, 1),
                    ("H+(aq)", &[0, 1, 0], 1, 4),
                ],
            ),
            Self::FerricChalcopyriteLeach => (
                "ferric-chalcopyrite-leach",
                &["Cu", "Fe", "S"],
                &[
                    ("CuFeS2(s)", &[1, 1, 2], 0, -1),
                    ("Fe[3+](aq)", &[0, 1, 0], 3, -4),
                    ("Cu[2+](aq)", &[1, 0, 0], 2, 1),
                    ("Fe[2+](aq)", &[0, 1, 0], 2, 5),
                    ("S(s)", &[0, 0, 1], 0, 2),
                ],
            ),
            Self::FerricChalcociteLeach => (
                "ferric-chalcocite-leach",
                &["Cu", "Fe", "S"],
                &[
                    ("Cu2S(s)", &[2, 0, 1], 0, -1),
                    ("Fe[3+](aq)", &[0, 1, 0], 3, -4),
                    ("Cu[2+](aq)", &[1, 0, 0], 2, 2),
                    ("Fe[2+](aq)", &[0, 1, 0], 2, 4),
                    ("S(s)", &[0, 0, 1], 0, 1),
                ],
            ),
            Self::FerricCovelliteLeach => (
                "ferric-covellite-leach",
                &["Cu", "Fe", "S"],
                &[
                    ("CuS(s)", &[1, 0, 1], 0, -1),
                    ("Fe[3+](aq)", &[0, 1, 0], 3, -2),
                    ("Cu[2+](aq)", &[1, 0, 0], 2, 1),
                    ("Fe[2+](aq)", &[0, 1, 0], 2, 2),
                    ("S(s)", &[0, 0, 1], 0, 1),
                ],
            ),
            Self::FerricBorniteLeach => (
                "ferric-bornite-leach",
                &["Cu", "Fe", "S"],
                &[
                    ("Cu5FeS4(s)", &[5, 1, 4], 0, -1),
                    ("Fe[3+](aq)", &[0, 1, 0], 3, -12),
                    ("Cu[2+](aq)", &[1, 0, 0], 2, 5),
                    ("Fe[2+](aq)", &[0, 1, 0], 2, 13),
                    ("S(s)", &[0, 0, 1], 0, 4),
                ],
            ),
        };
        let species = participants
            .iter()
            .map(|entry| SpeciesId::new(entry.0))
            .collect::<Result<Vec<_>, _>>()?;
        let elements = elements
            .iter()
            .map(|name| ElementId::new(*name))
            .collect::<Result<Vec<_>, _>>()?;
        let rows = (0..elements.len())
            .map(|row| participants.iter().map(|entry| entry.1[row]).collect())
            .collect();
        let elemental = ElementalMatrix::new(elements, species.clone(), rows)?;
        let id = ReactionId::new(name)?;
        let stoichiometry = StoichiometricMatrix::new(
            species.clone(),
            vec![id.clone()],
            participants.iter().map(|entry| vec![entry.3]).collect(),
        )?;
        let charge =
            ChargeVector::new(species, participants.iter().map(|entry| entry.2).collect())?;
        ConservedReaction::new(&elemental, &stoichiometry, &charge, id)
    }
}
