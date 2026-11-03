//! Finite simultaneous reaction stages sharing one complete species inventory.
//!
//! Requested extents are decisions supplied by a caller, not kinetics or an
//! equilibrium solution. Every reaction consumes its gross demand from the
//! INITIAL stocks. Products are credited only after every demand passes; they
//! can fund a subsequent explicit call, never another reaction in this stage.
//! Shared species must retain exactly the same elemental counts and charge in
//! every original fs-qty conservation certificate. Phase-specific IDs remain
//! distinct, and absent zero element rows do not create a false identity clash.
//!
//! Equation conservation is exact; reported amounts are binary64 quantities.
//! FMA product residuals and bounded error-free sum expansions protect shared
//! availability and small net changes. No reaction is clipped, reordered into
//! an implicit sequence, or rescaled to repair a failed material balance.

/// Source-bound standard enthalpy accounting for an admitted simultaneous stage.
#[path = "network/enthalpy.rs"]
pub mod enthalpy;

use super::{
    ConservedReaction, ReactionAmount, ReactionBalanceError, SpeciesDefinition, nonnegative,
};
use crate::{ReactionId, SpeciesId};
use fs_math::eft::two_sum;
use fs_qty::Amount;
use std::collections::BTreeMap;

/// Maximum active species in one shared inventory.
pub const MAX_NETWORK_SPECIES: usize = super::MAX_SPECIES;
/// Maximum explicitly requested reactions in one bounded stage.
pub const MAX_NETWORK_REACTIONS: usize = super::MAX_REACTIONS;

/// Refusal returns no partially advanced inventory.
#[derive(Debug, Clone, PartialEq)]
pub enum ReactionNetworkError {
    /// An original single-reaction or amount admission failed.
    Balance(
        /// Original refusal without changed meaning.
        ReactionBalanceError,
    ),
    /// Invalid resource count or incomplete inventory/extent axes.
    Invalid(
        /// Required correction at the call boundary.
        &'static str,
    ),
    /// The same active species ID was assigned different atoms or charge.
    ConflictingSpecies {
        /// Conflicting canonical species identity.
        species: SpeciesId,
    },
    /// A reaction extent cannot be associated with the complete network.
    Reaction {
        /// Missing, duplicated, or unknown reaction identity.
        reaction: ReactionId,
        /// Required correction at the call boundary.
        reason: &'static str,
    },
    /// Gross demand exceeds INITIAL stock; stage products are unavailable.
    InsufficientInventory {
        /// Shared reagent that would be overdrawn.
        species: SpeciesId,
        /// Initial available amount [mol].
        available: Amount,
        /// Rounded gross demand [mol]; may round to available at a tiny overdraw.
        required: Amount,
    },
    /// Overflow or a nonzero inventory change below binary64 resolution.
    Arithmetic(
        /// Failed numerical operation.
        &'static str,
    ),
}

impl core::fmt::Display for ReactionNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reaction network: {self:?}")
    }
}
impl std::error::Error for ReactionNetworkError {}
impl From<ReactionBalanceError> for ReactionNetworkError {
    fn from(error: ReactionBalanceError) -> Self {
        Self::Balance(error)
    }
}
type Result<T> = core::result::Result<T, ReactionNetworkError>;

/// One caller-requested reaction, retaining its written-equation extent basis.
#[derive(Debug, Clone)]
pub struct ReactionAdvance {
    /// Original admitted reaction identity.
    pub reaction: ReactionId,
    /// Requested nonnegative extent [mol of the integer equation as written].
    pub extent: Amount,
    /// Gross nonnegative quantities; term signs distinguish consumption/formation.
    pub amounts: Vec<ReactionAmount>,
}

/// Atomic stage result with complete species accounts in canonical order.
#[derive(Debug, Clone)]
pub struct NetworkAdvance {
    /// Each original reaction, extent and gross participant quantities.
    pub reactions: Vec<ReactionAdvance>,
    /// Gross consumption from the initial stock, including explicit zeroes.
    pub consumed: Vec<(SpeciesId, Amount)>,
    /// Gross formation credited after consumption admission, including zeroes.
    pub formed: Vec<(SpeciesId, Amount)>,
    /// Final complete inventory; directly usable as a later stage's input.
    pub final_inventory: Vec<(SpeciesId, Amount)>,
}

/// Canonical finite network of original, independently conserved reactions.
#[derive(Debug, Clone)]
pub struct ReactionNetwork {
    reactions: Vec<ConservedReaction>,
    species: Vec<SpeciesId>,
}

impl ReactionNetwork {
    /// Join distinct reaction IDs and exactly compatible active species.
    pub fn new(mut reactions: Vec<ConservedReaction>) -> Result<Self> {
        if reactions.is_empty() || reactions.len() > MAX_NETWORK_REACTIONS {
            return Err(ReactionNetworkError::Invalid(
                "network requires 1..=256 reactions",
            ));
        }
        reactions.sort_by(|a, b| a.id().cmp(b.id()));
        for pair in reactions.windows(2) {
            if pair[0].id() == pair[1].id() {
                return Err(ReactionNetworkError::Reaction {
                    reaction: pair[0].id().clone(),
                    reason: "duplicate reaction identity",
                });
            }
        }
        let mut definitions: BTreeMap<SpeciesId, SpeciesDefinition> = BTreeMap::new();
        for reaction in &reactions {
            for (term, definition) in reaction.terms.iter().zip(&reaction.definitions) {
                if let Some(previous) = definitions.get(&term.species) {
                    if previous != definition {
                        return Err(ReactionNetworkError::ConflictingSpecies {
                            species: term.species.clone(),
                        });
                    }
                } else {
                    if definitions.len() == MAX_NETWORK_SPECIES {
                        return Err(ReactionNetworkError::Invalid(
                            "network exceeds 256 active species",
                        ));
                    }
                    definitions.insert(term.species.clone(), definition.clone());
                }
            }
        }
        Ok(Self {
            species: definitions.into_keys().collect(),
            reactions,
        })
    }

    /// Complete canonical species axis, including initially absent products.
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }
    /// Original admitted reactions in canonical reaction-ID order.
    pub fn reactions(&self) -> &[ConservedReaction] {
        &self.reactions
    }

    /// Advance one simultaneous stage using complete, unique species and reaction
    /// axes in any input order. Every amount/extent is finite and nonnegative;
    /// explicit zero extents disable individual pathways for this stage.
    /// Gross reaction demands must fit the initial inventory, even when another
    /// reaction forms the same species. A nonzero NET change that rounds back to
    /// its initial amount refuses. Exact net-zero turnover remains reportable.
    pub fn advance_inventory(
        &self,
        initial: &[(SpeciesId, Amount)],
        requested: &[(ReactionId, Amount)],
    ) -> Result<NetworkAdvance> {
        if initial.len() != self.species.len() || requested.len() != self.reactions.len() {
            return Err(ReactionNetworkError::Invalid(
                "supply the complete species and reaction axes exactly once",
            ));
        }
        let mut inventory = BTreeMap::new();
        for (species, amount) in initial {
            nonnegative(amount.value(), "inventory must be finite and nonnegative")?;
            let reason = if self.species.binary_search(species).is_err() {
                Some("unknown network species")
            } else if inventory.insert(species, amount.value()).is_some() {
                Some("duplicate inventory")
            } else {
                None
            };
            if let Some(reason) = reason {
                return Err(ReactionBalanceError::Species {
                    species: species.clone(),
                    reason,
                }
                .into());
            }
        }
        let mut extents = BTreeMap::new();
        for (reaction, extent) in requested {
            nonnegative(extent.value(), "extent must be finite and nonnegative")?;
            let reason = if self
                .reactions
                .binary_search_by(|r| r.id().cmp(reaction))
                .is_err()
            {
                Some("unknown network reaction")
            } else if extents.insert(reaction, *extent).is_some() {
                Some("duplicate requested extent")
            } else {
                None
            };
            if let Some(reason) = reason {
                return Err(ReactionNetworkError::Reaction {
                    reaction: reaction.clone(),
                    reason,
                });
            }
        }
        let mut consumed = vec![ExactSum::default(); self.species.len()];
        let mut formed = consumed.clone();
        let mut advances = Vec::with_capacity(self.reactions.len());
        for reaction in &self.reactions {
            let extent = extents[reaction.id()];
            let amounts = reaction.amounts(extent)?;
            for amount in &amounts {
                let index = self
                    .species
                    .binary_search(&amount.term.species)
                    .map_err(|_| {
                        ReactionNetworkError::Invalid("reaction species is absent from network")
                    })?;
                let total = if amount.term.coefficient < 0 {
                    &mut consumed[index]
                } else {
                    &mut formed[index]
                };
                let product = amount.amount.value();
                let residual = (amount.term.coefficient as f64)
                    .abs()
                    .mul_add(extent.value(), -product);
                total.add(product)?;
                total.add(residual)?;
            }
            advances.push(ReactionAdvance {
                reaction: reaction.id().clone(),
                extent,
                amounts,
            });
        }
        let mut result = NetworkAdvance {
            reactions: advances,
            consumed: Vec::new(),
            formed: Vec::new(),
            final_inventory: Vec::new(),
        };
        for (index, species) in self.species.iter().enumerate() {
            let before = inventory[species];
            let demand = consumed[index].value()?;
            let formation = formed[index].value()?;
            let mut remainder = ExactSum::scalar(before);
            remainder.subtract(&consumed[index])?;
            if remainder.sign() < 0 {
                return Err(ReactionNetworkError::InsufficientInventory {
                    species: species.clone(),
                    available: Amount::new(before),
                    required: Amount::new(demand),
                });
            }
            for &part in &formed[index].parts {
                remainder.add(part)?;
            }
            let after = remainder.value()?;
            if after < 0.0 {
                return Err(ReactionNetworkError::Arithmetic(
                    "negative final species inventory",
                ));
            }
            if after == before {
                let mut net = formed[index].clone();
                net.subtract(&consumed[index])?;
                if net.sign() != 0 {
                    return Err(ReactionNetworkError::Arithmetic(
                        "nonzero stage change is below inventory resolution",
                    ));
                }
            }
            result.consumed.push((species.clone(), Amount::new(demand)));
            result
                .formed
                .push((species.clone(), Amount::new(formation)));
            result
                .final_inventory
                .push((species.clone(), Amount::new(after)));
        }
        Ok(result)
    }
}

// A bounded expansion of nonoverlapping binary64 components. Grow-expansion
// addition uses fs-math's original error-free transform, so availability tests
// retain the sign of a demand hidden by rounding its reported scalar total.
// Exact integer coefficients ensure FMA product residuals are representable
// multiples of the least subnormal; this does not assert arbitrary-product EFT.
#[derive(Clone, Default)]
struct ExactSum {
    parts: Vec<f64>,
}

impl ExactSum {
    fn scalar(value: f64) -> Self {
        Self {
            parts: if value == 0.0 {
                Vec::new()
            } else {
                vec![value]
            },
        }
    }
    fn add(&mut self, mut value: f64) -> Result<()> {
        if value == 0.0 {
            return Ok(());
        }
        let mut parts = Vec::with_capacity(self.parts.len() + 1);
        for &part in &self.parts {
            let (sum, residual) = two_sum(value, part);
            if !sum.is_finite() {
                return Err(ReactionNetworkError::Arithmetic(
                    "gross or final amount overflow",
                ));
            }
            if residual != 0.0 {
                parts.push(residual);
            }
            value = sum;
        }
        if value != 0.0 {
            parts.push(value);
        }
        if parts.len() > 4 * MAX_NETWORK_REACTIONS + 2 {
            return Err(ReactionNetworkError::Arithmetic(
                "stage expansion exceeds bounded work",
            ));
        }
        self.parts = parts;
        Ok(())
    }
    fn subtract(&mut self, other: &Self) -> Result<()> {
        // Cancel the largest positive demand first; this avoids a transient
        // negative overflow when a large stock and demand nearly coincide.
        for &part in other.parts.iter().rev() {
            self.add(-part)?;
        }
        Ok(())
    }
    fn sign(&self) -> i8 {
        self.parts
            .last()
            .map_or(0, |value| if *value > 0.0 { 1 } else { -1 })
    }
    fn value(&self) -> Result<f64> {
        let value = self.parts.iter().copied().fold(0.0, |sum, part| sum + part);
        if !value.is_finite() || (value == 0.0 && self.sign() != 0) {
            Err(ReactionNetworkError::Arithmetic(
                "stage amount is unrepresentable",
            ))
        } else {
            Ok(value)
        }
    }
}

#[cfg(test)]
#[path = "network/tests.rs"]
mod tests;
