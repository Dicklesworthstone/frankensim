//! Reference mineral/reagent equations with explicit final species.
//!
//! Pyrite endpoints follow USGS Open-File Report 99-348:
//! <https://pubs.usgs.gov/of/1999/0348/report.pdf>. Carbonate neutralizations
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
//! Malachite's net ionic acid balance is equation 10 in
//! <https://doi.org/10.3390/met10060833>; the pathway here retains dissolved CO2.
//! Net ionic Elsner cyanidation follows
//! <https://doi.org/10.1038/s41598-019-43383-4> and the EPA 1994 report
//! <https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=2000EF6J.TXT>.
//! Magnesite calcination and magnetite oxidation retain the crystalline/gas
//! endpoints of JANAF C-079, Mg-008, C-095, Fe-032, Fe-030 and O-029:
//! <https://janaf.nist.gov/>. Their standard reaction heats do not determine
//! a decomposition temperature, kiln duty, oxidation rate or phase stability.
//! These equations conserve atoms and charge; they do not select pH/Eh,
//! dissolved speciation, precipitates, gas transfer, rates or conversion.

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
    /// 4 Au(s) + 8 CN- + O2(g) + 2 H2O(l) -> 4 Au(CN)2- + 4 OH-.
    /// Net ionic Elsner balance; gas supply does not imply a dissolution rate.
    GoldCyanidation,
    /// 2 FeS2(s) + 7 O2(g) + 2 H2O(l) -> 2 Fe2+ + 4 SO4(2-) + 4 H+.
    PyriteOxidationFerrous,
    /// 4 FeS2(s) + 15 O2(g) + 14 H2O(l) -> 4 Fe(OH)3(s) + 8 SO4(2-) + 16 H+.
    /// Ideal ferric-hydroxide endpoint; no claim that a particular solid precipitates.
    PyriteOxidationFerricHydroxide,
    /// MgCO3(s) -> MgO(s) + CO2(g); crystalline carbonate/oxide endpoints.
    MagnesiteCalcination,
    /// 4 Fe3O4(s) + O2(g) -> 6 Fe2O3(s), with hematite as the product phase.
    /// One mole of reaction extent consumes four moles of magnetite.
    MagnetiteOxidationToHematite,
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
