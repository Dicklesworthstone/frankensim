//! Explicit neutral ion-pair transfer, separate from the adsorption fit.
//!
//! The supported bookkeeping choice moves one Au(CN)2- with one declared
//! monovalent countercation into one neutral adsorbed pair. It conserves the
//! two cyanide ligands and the countercation. It does not assert a microscopic
//! adsorption mechanism or supply standard formation thermochemistry.

use super::{CilError, Result};
use crate::cyanide_leach::{COMPLEX, CyanidationError};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ReactionId, SpeciesId, StoichiometricMatrix,
    stoichiometry::{
        ConservedReaction, ReactionBalanceError, minerals::MineralReaction,
        network::ReactionNetwork,
    },
};

/// Caller-declared, exactly conserved transfer of an intact gold ion pair.
///
/// The original aqueous gold-complex definition is verified against the
/// admitted Elsner equation. Counterion and adsorbed-pair charges come from
/// the reaction artifacts, never from their species labels.
#[derive(Clone, Debug)]
pub struct IonPairTransfer {
    reaction: ConservedReaction,
    counterion: SpeciesId,
    adsorbed_pair: SpeciesId,
}

impl IonPairTransfer {
    /// Au(CN)2[-](aq) + Na[+](aq) -> NaAu(CN)2(ads).
    ///
    /// This is an explicit neutral-pair accounting choice. Sodium stock must
    /// be supplied by the original bath and remains finite during transfer.
    pub fn sodium() -> Result<Self> {
        Self::from_reaction(
            sodium_reaction().map_err(CyanidationError::from)?,
            SpeciesId::new("Na[+](aq)").expect("fixed sodium identity"),
            SpeciesId::new("NaAu(CN)2(ads)").expect("fixed adsorbed-pair identity"),
        )
    }

    /// Admit a caller-supplied 1:1:1 neutral ion-pair equation.
    ///
    /// Exactly one original Au(CN)2- and one monovalent countercation must be
    /// consumed per neutral pair formed. The countercation cannot itself
    /// contain gold. An owning CIL bath additionally checks that this cation
    /// is an explicit, sufficient aqueous spectator and that all shared
    /// chemical definitions agree with its original chemistry.
    pub fn from_reaction(
        reaction: ConservedReaction,
        counterion: SpeciesId,
        adsorbed_pair: SpeciesId,
    ) -> Result<Self> {
        let complex = SpeciesId::new(COMPLEX).expect("fixed gold-complex identity");
        if counterion == complex || adsorbed_pair == complex || counterion == adsorbed_pair {
            return Err(CilError::Invalid("ion-pair species must be distinct"));
        }
        let expected = [(&complex, -1), (&counterion, -1), (&adsorbed_pair, 1)];
        if reaction.terms().len() != expected.len()
            || expected.iter().any(|(id, coefficient)| {
                !reaction
                    .terms()
                    .iter()
                    .any(|term| &term.species == *id && term.coefficient == *coefficient)
            })
        {
            return Err(CilError::Invalid(
                "ion-pair transfer must consume one gold complex and one counterion per pair",
            ));
        }
        if reaction.species_charge(&complex) != Some(-1)
            || reaction.species_charge(&counterion) != Some(1)
            || reaction.species_charge(&adsorbed_pair) != Some(0)
        {
            return Err(CilError::Invalid(
                "ion-pair transfer requires Au(CN)2-, a monovalent cation, and a neutral pair",
            ));
        }
        if reaction
            .species_elements(&counterion)
            .is_none_or(|elements| {
                elements
                    .iter()
                    .any(|(element, count)| element.as_str() == "Au" && *count != 0)
            })
        {
            return Err(CilError::Invalid(
                "ion-pair counterion must have an admitted gold-free definition",
            ));
        }
        ReactionNetwork::new(vec![
            MineralReaction::GoldCyanidationDissolvedOxygen
                .build()
                .map_err(CyanidationError::from)?,
            reaction.clone(),
        ])
        .map_err(CyanidationError::from)?;
        Ok(Self {
            reaction,
            counterion,
            adsorbed_pair,
        })
    }

    /// Original exactly conserved transfer equation.
    #[must_use]
    pub fn reaction(&self) -> &ConservedReaction {
        &self.reaction
    }

    /// Explicit aqueous monovalent countercation.
    #[must_use]
    pub fn counterion(&self) -> &SpeciesId {
        &self.counterion
    }

    /// Explicit neutral adsorbed species.
    #[must_use]
    pub fn adsorbed_pair(&self) -> &SpeciesId {
        &self.adsorbed_pair
    }
}

fn sodium_reaction() -> core::result::Result<ConservedReaction, ReactionBalanceError> {
    let species: Vec<_> = [COMPLEX, "Na[+](aq)", "NaAu(CN)2(ads)"]
        .into_iter()
        .map(|name| SpeciesId::new(name).expect("fixed ion-pair identity"))
        .collect();
    let id = ReactionId::new("gold-cyanide-sodium-ion-pair-transfer")
        .expect("fixed ion-pair reaction identity");
    let elemental = ElementalMatrix::new(
        ["Au", "C", "N", "Na"]
            .into_iter()
            .map(|name| ElementId::new(name).expect("fixed element identity"))
            .collect(),
        species.clone(),
        vec![vec![1, 0, 1], vec![2, 0, 2], vec![2, 0, 2], vec![0, 1, 1]],
    )?;
    let stoichiometric = StoichiometricMatrix::new(
        species.clone(),
        vec![id.clone()],
        vec![vec![-1], vec![-1], vec![1]],
    )?;
    let charge = ChargeVector::new(species, vec![-1, 1, 0])?;
    ConservedReaction::new(&elemental, &stoichiometric, &charge, id)
}

#[cfg(test)]
mod tests;
