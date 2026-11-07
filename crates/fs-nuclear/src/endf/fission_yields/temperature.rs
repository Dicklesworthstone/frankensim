//! Evaluated fission products on the original temperature-dependent burnup path.
//!
//! Each temperature/group channel retains that knot's products and prompt
//! energy. Its microscopic coefficient is a nonnegative temperature-basis view
//! of one shared table. Summing the partial channels therefore interpolates
//! sigma, sigma * product yield and sigma * prompt energy together, rather than
//! independently interpolating conditional yields or energy per event.
//!
//! The original incident spectrum and independent-yield evaluation must be the
//! same at every knot. Between-group flux amplitudes may subsequently change;
//! temperature-dependent within-group spectra are not inferred. Product decay,
//! finite-inventory propagation, transport and material heat retain their
//! existing owners. Source labels and interpolation do not certify nuclear data.

use std::sync::Arc;

use super::model::FissionYieldCollapse;
use super::{FissionYieldEvaluation, invalid, reserved};
use crate::burnup::{BurnupData, BurnupLimits, BurnupRate};
use crate::cross_sections::{
    CrossSectionData, CrossSectionLimits, CrossSectionTable, NeutronReaction,
};
use crate::endf::spectrum::FluxSpectrum;
use crate::endf::{EndfError, EndfMetadata, budget};
use crate::irradiation::{NeutronGroups, NuclideId};
use crate::{EnergyChannels, FissionError, positive};

/// One processed temperature's original condensation and prompt-only event
/// energies, in the original neutron-group order. A nonreacting group still
/// needs an explicit energy entry, but it creates no partial reaction channel.
#[derive(Clone, Copy, Debug)]
pub struct FissionYieldTemperatureKnot<'a> {
    pub yields: &'a FissionYieldCollapse,
    pub prompt_energy_per_group: &'a [EnergyChannels],
}

/// Global output bounds. The cross-section payload counts temperature knots
/// times neutron groups ONCE; partial channels share that immutable storage.
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct FissionYieldTemperatureLimits {
    pub burnup: BurnupLimits,
    pub cross_sections: CrossSectionLimits,
}

/// Append temperature-dependent evaluated fission through the original partial
/// reaction channels. At each knot this recovers its original weighted removal,
/// product and prompt-energy moments. Between knots the existing linear-Kelvin
/// microscopic law interpolates those three moments, without normalizing yields.
///
/// Temperatures must be strictly increasing and cover the model's entire
/// declared domain. A single knot admits only that exact temperature. A zero-
/// reaction knot is retained as a zero table row; neighboring channels cannot
/// interpolate across it. The complete original yield evaluation, processed
/// isotope/library revision, named energy grid and full incident spectrum must
/// agree across knots. There is no temperature or spectral extrapolation.
///
/// The original per-knot assembly owns product/export and duplicate-fission
/// admission. Its input and all accepted inventories remain unchanged on
/// success, refusal and cancellation. Global reaction/product bounds include
/// the original reactions and every new partial channel.
pub fn append_temperature_series_to_burnup(
    knots: &[FissionYieldTemperatureKnot<'_>],
    data: &BurnupData,
    exported_products: &[NuclideId],
    limits: FissionYieldTemperatureLimits,
    mut checkpoint: impl FnMut() -> Result<(), FissionError>,
) -> Result<BurnupData, EndfError> {
    checkpoint()?;
    let first = knots
        .first()
        .ok_or_else(|| invalid("at least one evaluated fission temperature knot required"))?;
    let groups = data.groups.len();
    let product_limit = limits
        .burnup
        .max_products
        .min(limits.burnup.propagation.maximum_edges)
        .min(fs_la::decay::MAX_DECAY_EDGES);
    if knots.len() > limits.cross_sections.max_temperatures
        || groups > limits.cross_sections.max_neutron_groups
        || groups > limits.burnup.max_neutron_groups
        || data.nuclides.len() > limits.burnup.propagation.maximum_nodes
        || data.nuclides.len() > fs_la::decay::MAX_DECAY_NODES
        || data.reactions.len() > limits.burnup.max_reactions
    {
        return Err(budget("fission temperature-series shape"));
    }
    if knots
        .len()
        .checked_mul(groups)
        .is_none_or(|count| count > limits.cross_sections.max_values)
    {
        return Err(budget("fission temperature-series cross-section values"));
    }
    if !positive(data.temperature_min_k)
        || !positive(data.temperature_max_k)
        || data.temperature_min_k > data.temperature_max_k
    {
        return Err(invalid(
            "finite ordered positive burnup temperature domain required",
        ));
    }
    let mut original_products = 0_usize;
    for reaction in &data.reactions {
        checkpoint()?;
        original_products = original_products
            .checked_add(reaction.products.len())
            .filter(|&count| count <= product_limit)
            .ok_or_else(|| budget("fission temperature-series retained products"))?;
    }
    let first_weighted = first.yields.reaction_weighted();
    let mut previous_temperature = 0.0;
    for knot in knots {
        checkpoint()?;
        let weighted = knot.yields.reaction_weighted();
        let temperature = weighted.metadata.temperature_k;
        if !positive(temperature)
            || temperature <= previous_temperature
            || knot.yields.parent() != first.yields.parent()
            || weighted.reaction != NeutronReaction::Fission
            || !same_processed_evaluation(&first_weighted.metadata, &weighted.metadata)
            || !same_evaluation(
                first.yields.source_evaluation(),
                knot.yields.source_evaluation(),
                &mut checkpoint,
            )?
            || !same_groups(&data.groups, weighted.flux.groups(), &mut checkpoint)?
            || !same_values(
                first_weighted.flux.values_m2_s(),
                weighted.flux.values_m2_s(),
                &mut checkpoint,
            )?
            || !same_spectrum(
                first.yields.source_spectrum(),
                knot.yields.source_spectrum(),
                &mut checkpoint,
            )?
        {
            return Err(invalid(
                "fission temperature knots require one original evaluation, processed isotope/revision, grid and incident spectrum",
            ));
        }
        if weighted.cross_sections_m2.len() != groups
            || knot.prompt_energy_per_group.len() != groups
        {
            return Err(invalid(
                "one collapsed cross section and prompt energy per original group required",
            ));
        }
        previous_temperature = temperature;
    }
    if data.temperature_min_k < first_weighted.metadata.temperature_k
        || data.temperature_max_k > previous_temperature
    {
        return Err(invalid(
            "processed temperature knots must cover the entire burnup domain",
        ));
    }

    let mut temperatures = reserved(knots.len(), "fission temperature memory")?;
    let mut rows = reserved(knots.len(), "fission temperature rows")?;
    for knot in knots {
        checkpoint()?;
        let weighted = knot.yields.reaction_weighted();
        temperatures.push(weighted.metadata.temperature_k);
        let mut row = reserved(groups, "fission temperature coefficients")?;
        for &value in &weighted.cross_sections_m2 {
            checkpoint()?;
            row.push(value);
        }
        rows.push(row);
    }
    let table = CrossSectionTable::new_with_checkpoint(
        CrossSectionData {
            data_id: format!(
                "{} | {} | fission temperature basis | {}",
                data.data_id,
                first.yields.evaluation_data_id(),
                first.yields.source_spectrum().data_id()
            ),
            nuclide: first.yields.parent(),
            reaction: NeutronReaction::Fission,
            groups: Arc::clone(&data.groups),
            temperatures_k: temperatures,
            cross_sections_m2: rows,
        },
        limits.cross_sections,
        &mut checkpoint,
    )?;

    let mut basis_data = data.clone();
    let mut additions = Vec::new();
    let mut added_products = 0_usize;
    for (temperature_index, knot) in knots.iter().enumerate() {
        checkpoint()?;
        let weighted = knot.yields.reaction_weighted();
        basis_data.temperature_min_k = weighted.metadata.temperature_k;
        basis_data.temperature_max_k = weighted.metadata.temperature_k;
        let remaining = BurnupLimits {
            max_reactions: limits.burnup.max_reactions - additions.len(),
            max_products: product_limit - added_products,
            ..limits.burnup
        };
        let admitted = knot.yields.append_to_burnup(
            &basis_data,
            knot.prompt_energy_per_group,
            exported_products,
            remaining,
            &mut checkpoint,
        )?;
        let mut channels = admitted.reactions.into_iter().skip(data.reactions.len());
        for (group, &sigma) in weighted.cross_sections_m2.iter().enumerate() {
            checkpoint()?;
            if sigma == 0.0 {
                continue;
            }
            let mut channel = channels
                .next()
                .ok_or_else(|| invalid("original fission assembly omitted a reacting group"))?;
            let partial = table.temperature_component(
                format!(
                    "{} | temperature knot {}",
                    channel.data_id, temperature_index
                ),
                temperature_index,
                group,
            )?;
            channel.rate = BurnupRate::TabulatedFission {
                cross_sections: Arc::new(partial),
            };
            added_products = added_products
                .checked_add(channel.products.len())
                .filter(|&count| count <= product_limit - original_products)
                .ok_or_else(|| budget("fission temperature-series retained products"))?;
            additions
                .try_reserve(1)
                .map_err(|_| budget("fission temperature-series channel memory"))?;
            additions.push(channel);
        }
        if channels.next().is_some() {
            return Err(invalid("original fission assembly returned extra channels"));
        }
    }
    checkpoint()?;
    let mut next = data.clone();
    next.reactions
        .try_reserve_exact(additions.len())
        .map_err(|_| budget("fission temperature-series reaction memory"))?;
    next.reactions.extend(additions);
    checkpoint()?;
    Ok(next)
}

fn same_processed_evaluation(left: &EndfMetadata, right: &EndfMetadata) -> bool {
    left.nuclide == right.nuclide
        && left.atomic_weight_ratio == right.atomic_weight_ratio
        && left.library_id == right.library_id
        && left.library_version == right.library_version
        && left.library_release == right.library_release
        && left.material_revision == right.material_revision
}

fn same_evaluation(
    left: &FissionYieldEvaluation,
    right: &FissionYieldEvaluation,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<bool, EndfError> {
    if std::ptr::eq(left, right) {
        return Ok(true);
    }
    if left.data_id() != right.data_id()
        || left.metadata() != right.metadata()
        || !same_values(left.products(), right.products(), checkpoint)?
    {
        return Ok(false);
    }
    if std::ptr::eq(left.points(), right.points()) {
        return Ok(true);
    }
    if left.points().len() != right.points().len() {
        return Ok(false);
    }
    for (left, right) in left.points().iter().zip(right.points()) {
        checkpoint()?;
        if left.energy_j != right.energy_j
            || left.interpolation_from_previous != right.interpolation_from_previous
            || !same_values(&left.values, &right.values, checkpoint)?
        {
            return Ok(false);
        }
    }
    Ok(true)
}

fn same_groups(
    left: &Arc<NeutronGroups>,
    right: &Arc<NeutronGroups>,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<bool, EndfError> {
    Ok(Arc::ptr_eq(left, right)
        || (left.data_id() == right.data_id()
            && same_values(left.edges_j(), right.edges_j(), checkpoint)?))
}

fn same_spectrum(
    left: &FluxSpectrum,
    right: &FluxSpectrum,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<bool, EndfError> {
    Ok(std::ptr::eq(left, right)
        || (left.data_id() == right.data_id()
            && same_values(left.energies_j(), right.energies_j(), checkpoint)?
            && same_values(
                left.flux_density_m2_s_j(),
                right.flux_density_m2_s_j(),
                checkpoint,
            )?
            && same_values(
                left.interpolation_regions(),
                right.interpolation_regions(),
                checkpoint,
            )?))
}

fn same_values<T: PartialEq>(
    left: &[T],
    right: &[T],
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<bool, EndfError> {
    if std::ptr::eq(left, right) {
        return Ok(true);
    }
    if left.len() != right.len() {
        return Ok(false);
    }
    for (left, right) in left.iter().zip(right) {
        checkpoint()?;
        if left != right {
            return Ok(false);
        }
    }
    Ok(true)
}

#[cfg(test)]
mod tests;
