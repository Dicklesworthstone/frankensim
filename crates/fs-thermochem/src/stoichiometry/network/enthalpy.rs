//! Standard chemical enthalpy for an admitted finite-inventory reaction stage.
//!
//! One shared source association serves every occurrence of a species. This is
//! reference-state chemical energy, not heat of mixing, a solution activity
//! model, sensible heat, a temperature solve, or a kinetic pathway selection.

use std::collections::BTreeMap;

use fs_qty::{Amount, Energy};

use super::{NetworkAdvance, ReactionNetwork, ReactionNetworkError};
use crate::stoichiometry::{FormationEnthalpySource, ReactionBalanceError, ReactionEnthalpy};
use crate::{ReactionId, SpeciesId};

/// One admitted reaction's source-bound standard enthalpy and actual extent.
#[derive(Debug, Clone)]
pub struct ReactionExtentEnthalpy {
    /// The integer equation whose extent was admitted against shared stocks.
    pub reaction: ReactionId,
    /// Moles of the integer equation as written.
    pub extent: Amount,
    /// Original formation evidence and usage receipts for all active species.
    pub standard: ReactionEnthalpy,
    /// Products-minus-reactants enthalpy change at this extent [J].
    pub enthalpy: Energy,
}

/// Finite material accounts and their compatible standard chemical energy.
#[derive(Debug, Clone)]
pub struct NetworkThermalAdvance {
    /// Complete, canonically ordered material accounts for this stage.
    pub inventory: NetworkAdvance,
    /// Canonically ordered individual contributions with original evidence.
    pub reactions: Vec<ReactionExtentEnthalpy>,
    /// Sum of the individual standard enthalpy changes [J].
    /// Negative is exothermic at constant pressure. Do not add it to a total
    /// energy formulation that already includes chemical formation enthalpy.
    pub enthalpy: Energy,
}

impl ReactionNetwork {
    /// Admit the shared inventory and evaluate its standard enthalpy change.
    ///
    /// Supply exactly one formation source for every network species, even for
    /// zero-extent reactions. Every selected source must declare the same typed
    /// temperature, pressure and elemental reference. This also applies to
    /// disconnected reactions, which otherwise could silently mix conventions.
    /// Phase and ionic-state selectors remain specific to each source. The
    /// caller supplies the species-to-pack association; it is not inferred
    /// from a material name or independently authenticated by this method.
    ///
    /// All stock and source checks succeed before any result is returned. The
    /// input inventories are never mutated and unavailable reagents never
    /// produce a partial heat result. No new uncertainty interval is inferred.
    pub fn advance_with_standard_enthalpy(
        &self,
        initial: &[(SpeciesId, Amount)],
        requested: &[(ReactionId, Amount)],
        sources: &[FormationEnthalpySource<'_>],
    ) -> Result<NetworkThermalAdvance, ReactionNetworkError> {
        let inventory = self.advance_inventory(initial, requested)?;
        if sources.len() != self.species().len() {
            return Err(ReactionBalanceError::Invalid(
                "supply exactly one formation source per network species",
            )
            .into());
        }
        let mut by_species = BTreeMap::new();
        let first = sources.first().ok_or(ReactionBalanceError::Invalid(
            "network has no formation sources",
        ))?;
        for source in sources {
            source.validate_reference()?;
            if self.species().binary_search(source.species).is_err() {
                return Err(ReactionBalanceError::Species {
                    species: source.species.clone(),
                    reason: "formation source is outside the reaction network",
                }
                .into());
            }
            if by_species.insert(source.species, source).is_some() {
                return Err(ReactionBalanceError::Species {
                    species: source.species.clone(),
                    reason: "duplicate network formation source",
                }
                .into());
            }
            if source.reference_axis != first.reference_axis
                || source.at.axes().get("temperature") != first.at.axes().get("temperature")
                || source.at.axes().get("pressure") != first.at.axes().get("pressure")
            {
                return Err(ReactionBalanceError::Invalid(
                    "network formation reference states differ",
                )
                .into());
            }
        }

        let mut reactions = Vec::with_capacity(self.reactions().len());
        let mut sum = 0.0_f64;
        let mut correction = 0.0_f64;
        for reaction in self.reactions() {
            let active = reaction
                .terms()
                .iter()
                .map(|term| {
                    let source = by_species.get(&term.species).ok_or_else(|| {
                        ReactionBalanceError::Species {
                            species: term.species.clone(),
                            reason: "missing network formation source",
                        }
                    })?;
                    Ok(FormationEnthalpySource {
                        species: source.species,
                        pack: source.pack,
                        at: source.at,
                        reference_axis: source.reference_axis,
                    })
                })
                .collect::<Result<Vec<_>, ReactionBalanceError>>()?;
            let standard = reaction.standard_enthalpy_change(&active)?;
            // The inventory owner has already admitted the exact complete
            // reaction axis. Looking up the same caller request cannot select
            // a different extent or invent a default for a missing reaction.
            let extent = requested
                .iter()
                .find(|(id, _)| id == reaction.id())
                .map(|(_, extent)| *extent)
                .ok_or(ReactionBalanceError::Invalid(
                    "missing admitted reaction extent",
                ))?;
            let enthalpy = standard.at_extent(extent)?;
            let value = enthalpy.value();
            let next = sum + value;
            if !next.is_finite() {
                return Err(
                    ReactionBalanceError::Arithmetic("network reaction enthalpy sum").into(),
                );
            }
            correction += if sum.abs() >= value.abs() {
                (sum - next) + value
            } else {
                (value - next) + sum
            };
            sum = next;
            reactions.push(ReactionExtentEnthalpy {
                reaction: reaction.id().clone(),
                extent,
                standard,
                enthalpy,
            });
        }
        let total = sum + correction;
        if !total.is_finite() {
            return Err(ReactionBalanceError::Arithmetic("network reaction enthalpy").into());
        }
        Ok(NetworkThermalAdvance {
            inventory,
            reactions,
            enthalpy: Energy::new(total),
        })
    }
}

#[cfg(test)]
#[path = "enthalpy_tests.rs"]
mod tests;
