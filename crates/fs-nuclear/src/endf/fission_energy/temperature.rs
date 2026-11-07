//! Evaluated prompt fission products with processed temperature feedback.
//!
//! Each knot combines matching product and prompt-energy condensations of the
//! same processed fission cross section and original incident spectrum. The
//! existing product-series owner interpolates partial reaction coefficients:
//! it preserves fission, product-production and prompt-energy rates at the
//! temperature knots and uses their declared linear law between knots.
//!
//! These are outgoing prompt product energies, including incident-neutron
//! kinetic input, as documented in [`super::burnup`]. Delayed release remains
//! the original product-decay owner's work during irradiation and shutdown.
//! This adapter does not broaden cross sections or infer a changing spectrum.

use std::sync::Arc;

use super::burnup::{PromptFissionRouting, same_values};
use super::model::{CollapsedFissionEnergy, PreparedFissionEnergy, ReferenceEnergyTreatment};
use super::{FissionEnergyComponent, FissionEnergyRelease};
use crate::FissionError;
use crate::burnup::BurnupData;
use crate::endf::fission_yields::model::FissionYieldCollapse;
use crate::endf::fission_yields::temperature::{
    FissionYieldTemperatureKnot, FissionYieldTemperatureLimits, append_temperature_series_to_burnup,
};
use crate::endf::reaction_weighted::TabulatedYield;
use crate::endf::{EndfError, budget};
use crate::irradiation::NuclideId;

/// One processed temperature, obtained by condensing both evaluated sources
/// against exactly the same original fission curve and incident spectrum.
#[derive(Clone, Copy, Debug)]
pub struct FissionEnergyTemperatureKnot<'a> {
    pub energy: &'a CollapsedFissionEnergy,
    pub yields: &'a FissionYieldCollapse,
}

/// Append a temperature-dependent fission series on the original burnup owner.
///
/// Temperature knots must be increasing and share the same original MT458
/// evaluation, energy continuation/preparation, product evaluation, incident
/// spectrum and group grid. At every
/// knot, full pointwise source matching, component routing and its explicit
/// J/fission roundoff allowance use the fixed-state bridge's original checks.
/// Existing represented product decays are retained without additional delayed
/// energy at the fission event.
///
/// The delegated product-series assembly interpolates rates, so an interior
/// state's energy and product production are weighted by the knot fission
/// rates. Independently interpolating energy/yield means would generally not
/// preserve those rates when the cross section changes with temperature.
///
/// `maximum_routing_roundoff_j` bounds each knot/group's coefficient routing,
/// not whole-run heat error. Interpolation errors and incident kinetic-energy
/// means remain observable through each supplied `CollapsedFissionEnergy`.
/// The source data are unchanged on success, refusal and cancellation.
#[allow(clippy::too_many_arguments)]
pub fn append_prompt_temperature_series_to_burnup(
    knots: &[FissionEnergyTemperatureKnot<'_>],
    data: &BurnupData,
    routing: PromptFissionRouting,
    exported_products: &[NuclideId],
    limits: FissionYieldTemperatureLimits,
    maximum_routing_roundoff_j: f64,
    mut checkpoint: impl FnMut() -> Result<(), FissionError>,
) -> Result<BurnupData, EndfError> {
    checkpoint()?;
    let first = knots
        .first()
        .ok_or_else(|| invalid("at least one evaluated fission temperature knot required"))?;
    if knots.len() > limits.cross_sections.max_temperatures {
        return Err(budget("evaluated fission temperatures"));
    }
    let groups = first.energy.reaction_weighted().flux.groups().len();
    if groups > limits.cross_sections.max_neutron_groups
        || groups > limits.burnup.max_neutron_groups
        || knots
            .len()
            .checked_mul(groups)
            .is_none_or(|values| values > limits.cross_sections.max_values)
    {
        return Err(budget("evaluated fission temperature values"));
    }
    // Bound all temporary energy rows before allocating any of them. Exact
    // grid/source compatibility is checked by the shared physical adapters.
    for knot in knots {
        checkpoint()?;
        if knot.energy.reaction_weighted().flux.groups().len() != groups {
            return Err(invalid(
                "matching evaluated fission temperature group counts required",
            ));
        }
    }
    let mut channels = Vec::new();
    channels
        .try_reserve_exact(knots.len())
        .map_err(|_| budget("evaluated fission temperature energy rows"))?;
    for knot in knots {
        checkpoint()?;
        if !same_evaluation(
            first.energy.evaluation(),
            knot.energy.evaluation(),
            &mut checkpoint,
        )? {
            return Err(invalid(
                "one original MT458 evaluation required across temperatures",
            ));
        }
        if !same_preparation(
            first.energy.source_preparation(),
            knot.energy.source_preparation(),
            &mut checkpoint,
        )? {
            return Err(invalid(
                "one fission energy continuation and prepared law required across temperatures",
            ));
        }
        channels.push(knot.energy.prompt_channels_for_burnup(
            knot.yields,
            routing,
            limits.burnup,
            maximum_routing_roundoff_j,
            &mut checkpoint,
        )?);
    }
    let mut product_knots = Vec::new();
    product_knots
        .try_reserve_exact(knots.len())
        .map_err(|_| budget("evaluated fission temperature assembly"))?;
    for (knot, channels) in knots.iter().zip(&channels) {
        checkpoint()?;
        product_knots.push(FissionYieldTemperatureKnot {
            yields: knot.yields,
            prompt_energy_per_group: channels,
        });
    }
    let mut assembled = append_temperature_series_to_burnup(
        &product_knots,
        data,
        exported_products,
        limits,
        &mut checkpoint,
    )?;
    for reaction in &mut assembled.reactions[data.reactions.len()..] {
        checkpoint()?;
        reaction.data_id = format!(
            "{} | {} | MF1/MT458 outgoing prompt products",
            reaction.data_id,
            first.energy.evaluation().data_id()
        );
    }
    checkpoint()?;
    Ok(assembled)
}

fn same_evaluation(
    left: &Arc<FissionEnergyRelease>,
    right: &Arc<FissionEnergyRelease>,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<bool, EndfError> {
    checkpoint()?;
    if Arc::ptr_eq(left, right) {
        return Ok(true);
    }
    if left.metadata() != right.metadata()
        || left.polynomial_order() != right.polynomial_order()
        || !same_values(
            left.data_id().as_bytes(),
            right.data_id().as_bytes(),
            checkpoint,
        )?
    {
        return Ok(false);
    }
    for component in FissionEnergyComponent::ALL {
        checkpoint()?;
        if left.reference_value_j(component) != right.reference_value_j(component)
            || !same_values(
                left.coefficients(component),
                right.coefficients(component),
                checkpoint,
            )?
        {
            return Ok(false);
        }
        match (left.tabulation(component), right.tabulation(component)) {
            (None, None) => {}
            (Some(left), Some(right)) => {
                if left.derivation() != right.derivation()
                    || !same_curve(&left.curve, &right.curve, checkpoint)?
                {
                    return Ok(false);
                }
            }
            _ => return Ok(false),
        }
    }
    Ok(true)
}

fn same_preparation(
    left: &PreparedFissionEnergy,
    right: &PreparedFissionEnergy,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<bool, EndfError> {
    checkpoint()?;
    if left.domain_j() != right.domain_j()
        || !same_reference(
            left.reference_treatment(),
            right.reference_treatment(),
            checkpoint,
        )?
    {
        return Ok(false);
    }
    for component in FissionEnergyComponent::ALL {
        checkpoint()?;
        if left.interpolation_error_j(component) != right.interpolation_error_j(component) {
            return Ok(false);
        }
    }
    let left = left.component_curves();
    let right = right.component_curves();
    if left.len() != right.len() {
        return Ok(false);
    }
    if std::ptr::eq(left, right) {
        return Ok(true);
    }
    for (left, right) in left.iter().zip(right) {
        checkpoint()?;
        if !same_curve(left, right, checkpoint)? {
            return Ok(false);
        }
    }
    Ok(true)
}

fn same_reference(
    left: &ReferenceEnergyTreatment,
    right: &ReferenceEnergyTreatment,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<bool, EndfError> {
    checkpoint()?;
    match (left, right) {
        (ReferenceEnergyTreatment::Reject, ReferenceEnergyTreatment::Reject)
        | (
            ReferenceEnergyTreatment::ConstantReference,
            ReferenceEnergyTreatment::ConstantReference,
        ) => Ok(true),
        (
            ReferenceEnergyTreatment::SherBeck2023 {
                neutron_yield: left,
                reference_nu_at_zero: left_zero,
            },
            ReferenceEnergyTreatment::SherBeck2023 {
                neutron_yield: right,
                reference_nu_at_zero: right_zero,
            },
        ) => {
            if left_zero != right_zero {
                return Ok(false);
            }
            if Arc::ptr_eq(left, right) {
                return Ok(true);
            }
            Ok(left.metadata() == right.metadata()
                && left.kind() == right.kind()
                && same_curve(left.curve(), right.curve(), checkpoint)?
                && same_values(
                    left.decay_constants_s_inv(),
                    right.decay_constants_s_inv(),
                    checkpoint,
                )?)
        }
        _ => Ok(false),
    }
}

fn same_curve(
    left: &TabulatedYield,
    right: &TabulatedYield,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<bool, EndfError> {
    Ok(same_values(
        left.data_id().as_bytes(),
        right.data_id().as_bytes(),
        checkpoint,
    )? && same_values(left.energies_j(), right.energies_j(), checkpoint)?
        && same_values(left.values(), right.values(), checkpoint)?
        && same_values(
            left.interpolation_regions(),
            right.interpolation_regions(),
            checkpoint,
        )?)
}

fn invalid(reason: &'static str) -> EndfError {
    FissionError::InvalidInput(reason).into()
}

#[cfg(test)]
#[path = "temperature/tests.rs"]
pub(crate) mod tests;
