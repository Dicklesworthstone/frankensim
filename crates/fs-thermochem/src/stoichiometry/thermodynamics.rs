//! Source-backed mixed-phase Gibbs energy and explicit activity diagnostics.
//!
//! IUPAC defines ln(K°) = -Delta_r G°/(R T) and
//! mu_i = mu_i° + R T ln(a_i):
//! <https://goldbook.iupac.org/terms/view/S05915> and
//! <https://goldbook.iupac.org/terms/view/A00115>.
//! Direct formation-Gibbs claims retain their elemental and ionic references;
//! a formation enthalpy minus an absolute species entropy is not a substitute.
//! The separate [`caloric`] route combines conserved REACTION enthalpy and
//! entropy sums. It retains both source properties and does not manufacture
//! per-species formation-Gibbs claims or replace directly tabulated values.
//!
//! The direct diagnostic takes caller-supplied dimensionless activities
//! relative to each retained standard state, including solvent and minerals.
//! The separate [`equilibrium`] solver advances finite inventories along ONE
//! reaction under explicit ideal-dilute and pure-condensed activity assumptions.
//! Neither route infers a pH, a calibrated activity model, rates or complete
//! aqueous speciation. A gas activity is a fugacity ratio; identifying it with
//! p/p° needs an additional ideal-gas assumption that these routes do not impose.
//! The separate [`acid_base`] owner solves coupled protonation families and
//! water under a common explicit ideal activity basis and charge balance.

/// Coupled ideal aqueous acid/base families, water and explicit spectator ions.
pub mod acid_base;
/// Standard reaction properties derived from matched enthalpy and entropy sources.
pub mod caloric;
/// Finite-inventory equilibrium of one reaction under explicit ideal activities.
pub mod equilibrium;

use std::collections::BTreeMap;

use fs_matdb::{MaterialAnswer, NormalizedPack, QueryPoint};
use fs_math::det;
use fs_qty::{Dimensionless, Dims, Pressure, QuantitySpec, Temperature};

use super::{
    ConservedReaction, FormationEnthalpySource, ReactionBalanceError, ReactionTerm, Result,
    product, representable,
};
use crate::{
    MolarEnergyQuantityV1, ReactionId, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    equilibrium::EquilibriumConstantStatusV1,
};

/// Directly tabulated standard formation Gibbs energy for a declared species.
/// The property pack does not independently authenticate that chemical identity.
pub struct FormationGibbsSource<'a> {
    /// Caller-declared species, including its intended phase.
    pub species: &'a SpeciesId,
    /// Original normalized property artifact.
    pub pack: &'a NormalizedPack,
    /// Complete typed reference-state and phase coordinates.
    pub at: &'a QueryPoint,
    /// Common `source-elemental-reference-*` dimensionless flag.
    pub reference_axis: &'a str,
}

/// One load-bearing formation-Gibbs claim, with its original uncertainty/receipt.
#[derive(Clone, Debug)]
pub struct FormationGibbsTerm {
    /// Active participant and signed integer stoichiometric coefficient.
    pub term: ReactionTerm,
    /// Normalized source artifact's declared name.
    pub pack_id: String,
    /// Original selected property evidence; no stronger authority is minted.
    pub answer: MaterialAnswer,
}

/// Standard reaction Gibbs energy and dimensionless thermodynamic equilibrium
/// constant at one common source-admitted temperature and reference pressure.
#[derive(Clone, Debug)]
pub struct StandardReactionGibbs {
    reaction: ReactionId,
    temperature: Temperature,
    reference_pressure: Pressure,
    reference_axis: String,
    molar: MolarEnergyQuantityV1,
    log_equilibrium_constant: Dimensionless,
    equilibrium_constant: Option<Dimensionless>,
    status: EquilibriumConstantStatusV1,
    terms: Vec<FormationGibbsTerm>,
}

impl ConservedReaction {
    /// Evaluate sum(nu_i Delta_f G_i°), retaining every original source receipt.
    ///
    /// Exactly one source per active species must declare the same typed T, p
    /// and elemental reference. Aqueous sources must own the unit-activity,
    /// infinite-dilution convention; charged participants must additionally own
    /// the H+(aq)=0 ionic reference. Caller-only flags cannot supply authority.
    /// Standard-state chemical equilibrium is distinct from plant validation.
    pub fn standard_gibbs_change(
        &self,
        sources: &[FormationGibbsSource<'_>],
    ) -> Result<StandardReactionGibbs> {
        if sources.len() != self.terms.len() {
            return Err(ReactionBalanceError::Invalid(
                "supply one formation Gibbs source per active species",
            ));
        }
        let first = sources
            .first()
            .ok_or(ReactionBalanceError::Invalid("no formation Gibbs sources"))?;
        let mut by_species = BTreeMap::new();
        for source in sources {
            if by_species.insert(source.species, source).is_some() {
                return Err(ReactionBalanceError::Species {
                    species: source.species.clone(),
                    reason: "duplicate formation Gibbs source",
                });
            }
            if source.reference_axis != first.reference_axis
                || source.at.axes().get("temperature") != first.at.axes().get("temperature")
                || source.at.axes().get("pressure") != first.at.axes().get("pressure")
            {
                return Err(ReactionBalanceError::Invalid(
                    "formation Gibbs reference states differ",
                ));
            }
        }

        let mut terms = Vec::with_capacity(self.terms.len());
        let mut sum = Sum::default();
        for (term, definition) in self.terms.iter().zip(&self.definitions) {
            let source =
                by_species
                    .get(&term.species)
                    .ok_or_else(|| ReactionBalanceError::Species {
                        species: term.species.clone(),
                        reason: "missing formation Gibbs source",
                    })?;
            // Share reference/phase/receipt admission with the enthalpy path;
            // only the explicitly requested property differs.
            let reference = FormationEnthalpySource {
                species: source.species,
                pack: source.pack,
                at: source.at,
                reference_axis: source.reference_axis,
            };
            let answer =
                reference.select_formation_energy("standard-molar-formation-gibbs-energy")?;
            validate_activity_reference(source, &answer, definition.charge)?;
            sum.add(product(
                term.coefficient as f64,
                answer.evidence.value.value,
                "reaction Gibbs contribution",
            )?)?;
            terms.push(FormationGibbsTerm {
                term: term.clone(),
                pack_id: source.pack.pack_id().to_owned(),
                answer,
            });
        }
        // Every selected source has already admitted positive typed T and p.
        let temperature =
            first
                .at
                .axes()
                .get("temperature")
                .copied()
                .ok_or(ReactionBalanceError::Invalid(
                    "formation Gibbs temperature is absent",
                ))?;
        let reference_pressure =
            first
                .at
                .axes()
                .get("pressure")
                .copied()
                .ok_or(ReactionBalanceError::Invalid(
                    "formation Gibbs pressure is absent",
                ))?;
        let molar = sum.finish()?;
        let rt = product(
            UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
            temperature,
            "Gibbs gas-constant temperature product",
        )?;
        let log_k = -molar / rt;
        representable(log_k, molar != 0.0, "log equilibrium constant")?;
        let direct = det::exp(log_k);
        let (status, equilibrium_constant) = if direct.is_infinite() {
            (EquilibriumConstantStatusV1::LogOnlyOverflow, None)
        } else if direct == 0.0 {
            (EquilibriumConstantStatusV1::LogOnlyUnderflow, None)
        } else if direct.is_finite() && direct > 0.0 {
            (
                EquilibriumConstantStatusV1::FinitePositive,
                Some(Dimensionless::new(direct)),
            )
        } else {
            return Err(ReactionBalanceError::Arithmetic("equilibrium constant"));
        };
        Ok(StandardReactionGibbs {
            reaction: self.id.clone(),
            temperature: Temperature::new(temperature),
            reference_pressure: Pressure::new(reference_pressure),
            reference_axis: first.reference_axis.to_owned(),
            molar: MolarEnergyQuantityV1::new(molar),
            log_equilibrium_constant: Dimensionless::new(log_k),
            equilibrium_constant,
            status,
            terms,
        })
    }
}

fn validate_activity_reference(
    source: &FormationGibbsSource<'_>,
    answer: &MaterialAnswer,
    charge: i128,
) -> Result<()> {
    let claim = source.pack.claims().claim(answer.receipt.selected).ok_or(
        ReactionBalanceError::Invalid("selected formation Gibbs claim is absent"),
    )?;
    let owns_flag = |axis: &str| {
        claim.validity.bound(axis) == Some((1.0, 1.0))
            && claim.validity.axis_quantities().get(axis)
                == Some(&QuantitySpec::dimensional(Dims::NONE))
            && source.at.axes().get(axis) == Some(&1.0)
            && source.at.axis_quantities().get(axis) == Some(&QuantitySpec::dimensional(Dims::NONE))
    };
    if owns_flag("source-phase-aqueous") {
        if !owns_flag("source-aqueous-unit-activity-infinite-dilution") {
            return Err(ReactionBalanceError::Invalid(
                "aqueous Gibbs source must declare its activity standard state",
            ));
        }
        if charge != 0
            && (!owns_flag("source-ionic-reference-hydrogen-zero")
                || !owns_flag("source-individual-ion"))
        {
            return Err(ReactionBalanceError::Invalid(
                "ionic Gibbs source must declare the H+ zero reference",
            ));
        }
    } else if charge != 0 || !owns_flag("source-pure-species") {
        return Err(ReactionBalanceError::Invalid(
            "non-aqueous Gibbs source requires a neutral pure-species reference",
        ));
    }
    Ok(())
}

impl StandardReactionGibbs {
    /// Identity of the admitted integer reaction.
    #[must_use]
    pub fn reaction(&self) -> &ReactionId {
        &self.reaction
    }
    /// Common source-admitted reference temperature [K].
    #[must_use]
    pub fn temperature(&self) -> Temperature {
        self.temperature
    }
    /// Common source reference pressure [Pa]; no pressure correction is inferred.
    #[must_use]
    pub fn reference_pressure(&self) -> Pressure {
        self.reference_pressure
    }
    /// Common original elemental-reference selector.
    #[must_use]
    pub fn reference_axis(&self) -> &str {
        &self.reference_axis
    }
    /// Standard reaction Gibbs energy [J/mol of the equation as written].
    #[must_use]
    pub fn molar(&self) -> MolarEnergyQuantityV1 {
        self.molar
    }
    /// Finite natural logarithm of the dimensionless thermodynamic K°.
    #[must_use]
    pub fn log_equilibrium_constant(&self) -> Dimensionless {
        self.log_equilibrium_constant
    }
    /// Direct K° when positive and finite; the log value remains authoritative.
    #[must_use]
    pub fn equilibrium_constant(&self) -> Option<Dimensionless> {
        self.equilibrium_constant
    }
    /// Whether direct exponentiation was representable.
    #[must_use]
    pub fn status(&self) -> EquilibriumConstantStatusV1 {
        self.status
    }
    /// Canonically ordered original formation claims and their usage receipts.
    #[must_use]
    pub fn terms(&self) -> &[FormationGibbsTerm] {
        &self.terms
    }

    /// Evaluate Delta_r G = R T (ln Q - ln K°) at supplied activities.
    ///
    /// Supply all and only active species, including condensed phases, with
    /// finite strictly positive dimensionless activity. Activities may exceed
    /// one. Zero-activity boundary limits are outside this logarithmic slice.
    /// Activities must use the original source conventions; no activity model
    /// is inferred. The result is dG/dxi at this state, not the free-energy
    /// change over a finite advance whose activities would change.
    pub fn at_activities(
        &self,
        activities: &[(SpeciesId, Dimensionless)],
    ) -> Result<ReactionActivityEvaluation<'_>> {
        if activities.len() != self.terms.len() {
            return Err(ReactionBalanceError::Invalid(
                "supply one dimensionless activity per active species",
            ));
        }
        let mut by_species = BTreeMap::new();
        for (species, activity) in activities {
            if !activity.value().is_finite() || activity.value() <= 0.0 {
                return Err(ReactionBalanceError::Species {
                    species: species.clone(),
                    reason: "activity must be finite and strictly positive",
                });
            }
            if by_species.insert(species, *activity).is_some() {
                return Err(ReactionBalanceError::Species {
                    species: species.clone(),
                    reason: "duplicate activity",
                });
            }
        }
        let mut sum = Sum::default();
        let mut canonical = Vec::with_capacity(self.terms.len());
        for source in &self.terms {
            let term = &source.term;
            let activity =
                by_species
                    .get(&term.species)
                    .ok_or_else(|| ReactionBalanceError::Species {
                        species: term.species.clone(),
                        reason: "missing active species activity",
                    })?;
            sum.add(product(
                term.coefficient as f64,
                det::ln(activity.value()),
                "log reaction quotient contribution",
            )?)?;
            canonical.push((term.species.clone(), *activity));
        }
        let log_q = sum.finish()?;
        let log_ratio = log_q - self.log_equilibrium_constant.value();
        representable(log_ratio, false, "log reaction quotient over equilibrium")?;
        let rt = product(
            UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
            self.temperature.value(),
            "activity gas-constant temperature product",
        )?;
        let molar = product(rt, log_ratio, "activity reaction Gibbs energy")?;
        Ok(ReactionActivityEvaluation {
            standard: self,
            activities: canonical,
            log_quotient: Dimensionless::new(log_q),
            log_quotient_over_equilibrium: Dimensionless::new(log_ratio),
            molar: MolarEnergyQuantityV1::new(molar),
        })
    }
}

/// Local thermodynamic diagnostic retaining the standard result and all inputs.
/// A finite or near-zero residual is not a certified equilibrium or rate claim.
#[derive(Debug)]
pub struct ReactionActivityEvaluation<'a> {
    /// Source-bound standard result at the original reference T and p.
    pub standard: &'a StandardReactionGibbs,
    /// Canonical, explicit activities relative to the selected source states.
    pub activities: Vec<(SpeciesId, Dimensionless)>,
    /// Natural logarithm of Q; Q is never formed as an overflowing product.
    pub log_quotient: Dimensionless,
    /// ln(Q/K°); negative favors the written forward direction at this state.
    pub log_quotient_over_equilibrium: Dimensionless,
    /// dG/dxi [J/mol], products minus reactants at the supplied activities.
    pub molar: MolarEnergyQuantityV1,
}

impl ReactionActivityEvaluation<'_> {
    /// Forward affinity [J/mol] is minus the reaction Gibbs energy.
    #[must_use]
    pub fn affinity(&self) -> MolarEnergyQuantityV1 {
        MolarEnergyQuantityV1::new(-self.molar.value())
    }
}

#[derive(Default)]
struct Sum {
    sum: f64,
    correction: f64,
}

impl Sum {
    fn add(&mut self, value: f64) -> Result<()> {
        let next = self.sum + value;
        representable(next, false, "thermodynamic sum")?;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        representable(self.correction, false, "thermodynamic sum correction")?;
        self.sum = next;
        Ok(())
    }

    fn finish(self) -> Result<f64> {
        let value = self.sum + self.correction;
        representable(value, false, "thermodynamic result")?;
        Ok(value)
    }
}

#[cfg(test)]
mod tests;
