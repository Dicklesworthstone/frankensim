//! Reference mineral/reagent equations with explicit final species.
//!
//! Pyrite endpoints follow USGS Open-File Report 99-348:
//! <https://pubs.usgs.gov/of/1999/0348/report.pdf>. Carbonate neutralizations
//! combine mineral dissolution with CO3(2-) + 2 H+ -> CO2(aq) + H2O from
//! the USGS PHREEQC database:
//! <https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc.v1/html/phqc_66.htm>.
//! Lime hydration is CaO + H2O -> Ca(OH)2; source-bound formation enthalpies
//! are supplied separately.
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
    /// 2 FeS2(s) + 7 O2(g) + 2 H2O(l) -> 2 Fe2+ + 4 SO4(2-) + 4 H+.
    PyriteOxidationFerrous,
    /// 4 FeS2(s) + 15 O2(g) + 14 H2O(l) -> 4 Fe(OH)3(s) + 8 SO4(2-) + 16 H+.
    /// Ideal ferric-hydroxide endpoint; no claim that a particular solid precipitates.
    PyriteOxidationFerricHydroxide,
}

// One species: stable phase/charge id, element counts, charge, signed coefficient.
type Participant = (&'static str, &'static [i128], i128, i128);

impl MineralReaction {
    /// Build through the existing exact atom/charge checker. Pyrite equations
    /// use integer multiples of two/four FeS2; an extent is per equation as
    /// written, not implicitly per mole of pyrite. Species ids include phase.
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
