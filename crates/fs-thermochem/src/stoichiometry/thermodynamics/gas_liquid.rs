//! Ideal gas/aqueous partition from the original two formation-Gibbs sources.
//!
//! For the exactly conserved equation `species(g) -> species(aq)`,
//! `a_aq = K * p_partial/p_ref`, with `ln K = -Delta_r G/(R T)`.
//! This owner explicitly adopts ideal gas fugacity and the common aqueous
//! approximation `a_aq=n_aq/n_ref`. The supplied `n_ref` is an amount at unit
//! activity; a unit-molality reference does not become a concentration basis.
//! IUPAC distinguishes molality/pressure and concentration/pressure Henry
//! constants: <https://goldbook.iupac.org/terms/view/09026> and
//! <https://goldbook.iupac.org/terms/view/09024>.
//!
//! A finite ideal headspace also satisfies `p_partial=n_g R T/V_g`.
//! Thus `n_aq/n_g = n_ref K R T/(V_g p_ref)` and a closed two-phase total
//! has an algebraic partition. No transfer coefficient, time evolution,
//! measured salinity correction, solvent-density conversion, dissolved
//! reaction or gas-mixture EOS is inferred. The source temperature is fixed.
//! Standard reference pressure is not an experimentally validated range of
//! actual gas partial pressures; this closure is an explicit ideal model.

use super::{FormationGibbsSource, StandardReactionGibbs, Sum};
use crate::{
    SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::{ConservedReaction, ReactionBalanceError},
};
use fs_math::det;
use fs_qty::{Amount, Dimensionless, Dims, Pressure, QuantitySpec, Volume};

/// Invalid partition declarations, source admission, or numerical range loss.
#[derive(Clone, Debug, PartialEq)]
pub enum GasLiquidError {
    /// Original exact reaction or original selected property source refused.
    Source(ReactionBalanceError),
    /// The equation, phase, basis, stock or volume does not fit this model.
    Invalid(&'static str),
    /// A required positive result or independent residual is unrepresentable.
    Arithmetic(&'static str),
}
impl core::fmt::Display for GasLiquidError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ideal gas/aqueous partition: {self:?}")
    }
}
impl std::error::Error for GasLiquidError {}
impl From<ReactionBalanceError> for GasLiquidError {
    fn from(error: ReactionBalanceError) -> Self {
        Self::Source(error)
    }
}
type Result<T> = core::result::Result<T, GasLiquidError>;

/// Original neutral transfer equation and its matched standard-state evidence.
#[derive(Clone, Debug)]
pub struct IdealGasLiquidPartition {
    reaction: ConservedReaction,
    standard: StandardReactionGibbs,
    gas: SpeciesId,
    aqueous: SpeciesId,
}

/// Algebraic closed-inventory partition at the retained source temperature.
#[derive(Clone, Debug)]
pub struct FiniteGasLiquidEquilibrium {
    model: IdealGasLiquidPartition,
    /// Supplied total of both phase amounts [mol].
    pub total_amount: Amount,
    /// Explicit aqueous amount corresponding to unit activity [mol].
    pub amount_at_unit_activity: Amount,
    /// Fixed ideal-gas headspace volume [m3].
    pub gas_volume: Volume,
    /// Retained gas amount [mol], independently represented even when tiny.
    pub gas_amount: Amount,
    /// Retained aqueous amount [mol], independently represented even when tiny.
    pub aqueous_amount: Amount,
    /// Species partial pressure in the ideal headspace [Pa].
    pub partial_pressure: Pressure,
    /// Independently reduced `gas + aqueous - total` [mol].
    pub component_balance_residual: Amount,
    /// Recomputed `ln(a_aq/(p_partial/p_ref)/K)`; undefined for an empty total.
    /// This is a numerical diagnostic, not a source-uncertainty certificate.
    pub log_equilibrium_residual: Option<Dimensionless>,
}
impl FiniteGasLiquidEquilibrium {
    /// Retained original equation, phase identities and property receipts.
    pub fn model(&self) -> &IdealGasLiquidPartition {
        &self.model
    }
}

impl IdealGasLiquidPartition {
    /// Admit exactly one neutral gas molecule to the same neutral aqueous
    /// formula with coefficients -1/+1. No scaled or dissociating equation is
    /// silently reinterpreted. Original source claims must own their phases,
    /// aqueous infinite-dilution reference, and common T/p/reference convention.
    pub fn from_sources(
        reaction: ConservedReaction,
        gas: SpeciesId,
        aqueous: SpeciesId,
        sources: &[FormationGibbsSource<'_>],
    ) -> Result<Self> {
        if gas == aqueous
            || reaction.terms().len() != 2
            || !reaction
                .terms()
                .iter()
                .any(|t| t.species == gas && t.coefficient == -1)
            || !reaction
                .terms()
                .iter()
                .any(|t| t.species == aqueous && t.coefficient == 1)
            || reaction
                .definitions
                .iter()
                .any(|definition| definition.charge != 0)
        {
            return Err(GasLiquidError::Invalid(
                "require two distinct neutral phase species with gas -1 and aqueous +1",
            ));
        }
        // Exact elemental identity follows from the original admitted conserved
        // two-species -1/+1 equation, not by parsing a species-name suffix.
        let standard = reaction.standard_gibbs_change(sources)?;
        for source in sources {
            let gas_phase = *source.species == gas;
            require_phase(source, &standard, gas_phase)?;
        }
        Ok(Self {
            reaction,
            standard,
            gas,
            aqueous,
        })
    }

    /// Original exactly conserved unit transfer equation.
    pub fn reaction(&self) -> &ConservedReaction {
        &self.reaction
    }
    /// Original standard Gibbs result and each selected source receipt.
    pub fn standard(&self) -> &StandardReactionGibbs {
        &self.standard
    }
    /// Original gas species identity.
    pub fn gas(&self) -> &SpeciesId {
        &self.gas
    }
    /// Original aqueous species identity.
    pub fn aqueous(&self) -> &SpeciesId {
        &self.aqueous
    }

    /// Ideal equilibrium aqueous stock at a declared nonnegative partial
    /// pressure. `n_ref` remains independent of solution volume and density.
    pub fn equilibrium_aqueous_amount(
        &self,
        partial_pressure: Pressure,
        amount_at_unit_activity: Amount,
    ) -> Result<Amount> {
        let pressure = nonnegative(partial_pressure.value(), "invalid species partial pressure")?;
        let basis = positive(
            amount_at_unit_activity.value(),
            "positive aqueous activity amount required",
        )?;
        if pressure == 0.0 {
            return Ok(Amount::new(0.0));
        }
        let log_amount = sum([
            det::ln(basis),
            self.standard.log_equilibrium_constant().value(),
            det::ln(pressure),
            -det::ln(self.standard.reference_pressure().value()),
        ])?;
        Ok(Amount::new(exp_positive(
            log_amount,
            "equilibrium aqueous amount loses range",
        )?))
    }

    /// Ideal species partial pressure `n_g R T/V_g`, at the source temperature.
    /// Other headspace gases do not replace this species' amount or pressure.
    pub fn ideal_partial_pressure(
        &self,
        gas_amount: Amount,
        gas_volume: Volume,
    ) -> Result<Pressure> {
        let amount = nonnegative(gas_amount.value(), "invalid gas amount")?;
        let volume = positive(gas_volume.value(), "positive gas volume required")?;
        if amount == 0.0 {
            return Ok(Pressure::new(0.0));
        }
        let pressure = exp_positive(
            sum([
                det::ln(amount),
                det::ln(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K),
                det::ln(self.standard.temperature().value()),
                -det::ln(volume),
            ])?,
            "ideal partial pressure loses range",
        )?;
        Ok(Pressure::new(pressure))
    }

    /// Close a finite gas-plus-aqueous component total. The gas volume and
    /// aqueous activity-reference amount remain distinct explicit inputs.
    /// Both positive tails are retained independently; unrepresentable positive
    /// stocks or pressure refuse instead of being clipped to zero.
    pub fn finite_equilibrium(
        &self,
        total: Amount,
        amount_at_unit_activity: Amount,
        gas_volume: Volume,
    ) -> Result<FiniteGasLiquidEquilibrium> {
        let amount = nonnegative(total.value(), "invalid combined phase amount")?;
        let basis = positive(
            amount_at_unit_activity.value(),
            "positive aqueous activity amount required",
        )?;
        let volume = positive(gas_volume.value(), "positive gas volume required")?;
        let (gas, aqueous) = if amount == 0.0 {
            (0.0, 0.0)
        } else {
            let log_ratio = sum([
                det::ln(basis),
                self.standard.log_equilibrium_constant().value(),
                det::ln(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K),
                det::ln(self.standard.temperature().value()),
                -det::ln(volume),
                -det::ln(self.standard.reference_pressure().value()),
            ])?;
            // Divide by 1+small, avoiding overflow in the dominant ratio.
            // The minor stock is not obtained by subtracting from the total.
            let small = det::exp(-log_ratio.abs());
            let major = amount / (1.0 + small);
            let product = major * small;
            let minor = if product > 0.0 {
                product
            } else {
                // The fraction itself may underflow while amount*fraction is
                // representable. Keep that tail in logarithmic coordinates.
                exp_positive(
                    sum([det::ln(amount), -log_ratio.abs(), -det::ln(1.0 + small)])?,
                    "positive phase amount loses range",
                )?
            };
            if !major.is_finite() || major <= 0.0 {
                return Err(GasLiquidError::Arithmetic(
                    "dominant phase amount loses range",
                ));
            }
            if log_ratio >= 0.0 {
                (minor, major)
            } else {
                (major, minor)
            }
        };
        let pressure = self.ideal_partial_pressure(Amount::new(gas), gas_volume)?;
        let balance = sum([gas, aqueous, -amount])?;
        let log_residual = if amount == 0.0 {
            None
        } else {
            Some(Dimensionless::new(sum([
                det::ln(aqueous),
                -det::ln(basis),
                -det::ln(pressure.value()),
                det::ln(self.standard.reference_pressure().value()),
                -self.standard.log_equilibrium_constant().value(),
            ])?))
        };
        Ok(FiniteGasLiquidEquilibrium {
            model: self.clone(),
            total_amount: total,
            amount_at_unit_activity,
            gas_volume,
            gas_amount: Amount::new(gas),
            aqueous_amount: Amount::new(aqueous),
            partial_pressure: pressure,
            component_balance_residual: Amount::new(balance),
            log_equilibrium_residual: log_residual,
        })
    }
}

fn require_phase(
    source: &FormationGibbsSource<'_>,
    standard: &StandardReactionGibbs,
    gas: bool,
) -> Result<()> {
    let selected = standard
        .terms()
        .iter()
        .find(|t| t.term.species == *source.species)
        .ok_or(GasLiquidError::Invalid(
            "missing selected phase participant",
        ))?;
    let claim = source
        .pack
        .claims()
        .claim(selected.answer.receipt.selected)
        .ok_or(GasLiquidError::Invalid("missing selected phase source"))?;
    let owns = |axis: &str| {
        claim.validity.bound(axis) == Some((1.0, 1.0))
            && claim.validity.axis_quantities().get(axis)
                == Some(&QuantitySpec::dimensional(Dims::NONE))
            && source.at.axes().get(axis) == Some(&1.0)
            && source.at.axis_quantities().get(axis) == Some(&QuantitySpec::dimensional(Dims::NONE))
    };
    let required = if gas {
        "source-phase-gas"
    } else {
        "source-phase-aqueous"
    };
    let convention = if gas {
        "source-pure-species"
    } else {
        "source-aqueous-unit-activity-infinite-dilution"
    };
    if !owns(required)
        || !owns(convention)
        || [
            "source-phase-solid",
            "source-phase-liquid",
            "source-phase-gas",
            "source-phase-aqueous",
        ]
        .iter()
        .any(|phase| *phase != required && owns(phase))
    {
        return Err(GasLiquidError::Invalid(
            "selected source must own its declared gas or aqueous phase",
        ));
    }
    Ok(())
}
fn nonnegative(value: f64, reason: &'static str) -> Result<f64> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(GasLiquidError::Invalid(reason))
    }
}
fn positive(value: f64, reason: &'static str) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(GasLiquidError::Invalid(reason))
    }
}
fn exp_positive(value: f64, reason: &'static str) -> Result<f64> {
    let result = det::exp(value);
    if result.is_finite() && result > 0.0 {
        Ok(result)
    } else {
        Err(GasLiquidError::Arithmetic(reason))
    }
}
fn sum<const N: usize>(values: [f64; N]) -> Result<f64> {
    let mut sum = Sum::default();
    for value in values {
        sum.add(value)?;
    }
    Ok(sum.finish()?)
}

#[cfg(test)]
#[path = "gas_liquid/tests.rs"]
mod tests;
