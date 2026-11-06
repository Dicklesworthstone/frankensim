//! Import-time adapters into the original condensation and burnup owners.
//!
//! Yields use reaction-rate weighting, not flux weighting alone. Each neutron
//! group becomes a disjoint fission channel, so its product distribution still
//! responds to changes in group flux. Cross sections are never duplicated across
//! channels. Subsequent decay and heating remain the original burnup/thermal
//! owners' work; this module creates no additional inventory or propagator.

use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

use super::{FissionYieldEvaluation, invalid, reserved};
use crate::burnup::{BurnupData, BurnupLimits, BurnupProduct, BurnupRate, BurnupReaction};
use crate::cross_sections::NeutronReaction;
use crate::endf::reaction_weighted::{
    ReactionWeightedCollapse, ReactionWeightedLimits, TabulatedYield, TabulatedYieldData,
    TabulatedYieldLimits, collapse_reaction_weighted,
};
use crate::endf::spectrum::FluxSpectrum;
use crate::endf::{
    EndfError, EndfInterpolation, EndfInterpolationRegion, PointwiseCrossSection, budget,
};
use crate::irradiation::{NeutronGroups, NuclideId};
use crate::{EnergyChannels, FissionError};

/// The existing weighted-collapse bounds plus temporary curve-storage bounds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FissionYieldCollapseLimits {
    pub weighting: ReactionWeightedLimits,
    /// Product count times represented incident-energy points (at least two).
    pub max_tabulated_values: usize,
}

impl Default for FissionYieldCollapseLimits {
    fn default() -> Self {
        Self {
            weighting: ReactionWeightedLimits::default(),
            max_tabulated_values: 1_000_000,
        }
    }
}

/// Product identities bound to the original reaction-weighted result. Construction
/// is restricted to independent yields; MT459 cannot enter this physical bridge.
#[derive(Clone, Debug, PartialEq)]
pub struct FissionYieldCollapse {
    evaluation_data_id: Arc<str>,
    parent: NuclideId,
    products: Arc<[NuclideId]>,
    weighted: ReactionWeightedCollapse,
    source_cross_section: PointwiseCrossSection,
    source_spectrum: FluxSpectrum,
}

impl FissionYieldEvaluation {
    /// Condense every independent product using the supplied processed fission
    /// curve and incident spectrum. Both source parents must match, including
    /// isomer. The caller remains responsible for evaluation compatibility.
    ///
    /// Multiple energy rows never extrapolate. A source-declared LE=0 law is
    /// constant on the explicitly supplied group domain. Reported uncertainties
    /// remain available in the evaluation; no covariance is inferred or summed.
    pub fn collapse(
        &self,
        cross_section: &PointwiseCrossSection,
        spectrum: &FluxSpectrum,
        groups: Arc<NeutronGroups>,
        limits: FissionYieldCollapseLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<FissionYieldCollapse, EndfError> {
        checkpoint()?;
        if cross_section.reaction() != NeutronReaction::Fission
            || cross_section.metadata().nuclide != self.metadata.nuclide
        {
            return Err(invalid(
                "independent yields require the same fission parent and isomer",
            ));
        }
        if self.products.len() > limits.weighting.max_curves {
            return Err(budget("independent product curves"));
        }
        let point_count = self.points.len().max(2);
        let values = point_count
            .checked_mul(self.products.len())
            .ok_or_else(|| budget("fission-yield curve values"))?;
        if values > limits.max_tabulated_values {
            return Err(budget("fission-yield curve values"));
        }
        let mut energies = reserved(point_count, "fission-yield energy memory")?;
        let mut regions: Vec<EndfInterpolationRegion> =
            reserved(point_count - 1, "fission-yield regions")?;
        if self.is_energy_independent() {
            let edges = groups.edges_j();
            energies.extend([edges[0], edges[edges.len() - 1]]);
            regions.push(EndfInterpolationRegion {
                end_point: 2,
                interpolation: EndfInterpolation::Linear,
            });
        } else {
            for (index, point) in self.points.iter().enumerate() {
                checkpoint()?;
                energies.push(point.energy_j);
                if let Some(interpolation) = point.interpolation_from_previous {
                    if let Some(last) = regions
                        .last_mut()
                        .filter(|last| last.interpolation == interpolation)
                    {
                        last.end_point = index + 1;
                    } else {
                        regions.push(EndfInterpolationRegion {
                            end_point: index + 1,
                            interpolation,
                        });
                    }
                }
            }
        }
        let mut curves = reserved(self.products.len(), "fission-product curve memory")?;
        for (index, product) in self.products.iter().enumerate() {
            checkpoint()?;
            let mut values = reserved(point_count, "fission-product yield memory")?;
            if self.is_energy_independent() {
                values.resize(2, self.points[0].values[index].yield_per_fission);
            } else {
                for point in self.points.iter() {
                    checkpoint()?;
                    values.push(point.values[index].yield_per_fission);
                }
            }
            curves.push(TabulatedYield::new_with_checkpoint(
                TabulatedYieldData {
                    data_id: format!(
                        "{} | MF8/MT454 Z{} A{} m{}",
                        self.data_id,
                        product.atomic_number(),
                        product.mass_number(),
                        product.isomer()
                    ),
                    energies_j: energies.clone(),
                    values,
                    interpolation_regions: regions.clone(),
                },
                TabulatedYieldLimits {
                    max_points: point_count,
                    max_interpolation_regions: regions.len(),
                },
                &mut checkpoint,
            )?);
        }
        let weighted = collapse_reaction_weighted(
            cross_section,
            spectrum,
            &curves,
            groups,
            limits.weighting,
            &mut checkpoint,
        )?;
        checkpoint()?;
        Ok(FissionYieldCollapse {
            evaluation_data_id: Arc::clone(&self.data_id),
            parent: self.metadata.nuclide,
            products: Arc::clone(&self.products),
            weighted,
            source_cross_section: cross_section.clone(),
            source_spectrum: spectrum.clone(),
        })
    }
}

impl FissionYieldCollapse {
    /// Original pointwise reaction, shared so subsequent energy/data assembly
    /// can verify full source correspondence rather than relying on labels.
    #[must_use]
    pub fn source_cross_section(&self) -> &PointwiseCrossSection {
        &self.source_cross_section
    }

    /// Original incident spectrum, including its within-group shape.
    #[must_use]
    pub fn source_spectrum(&self) -> &FluxSpectrum {
        &self.source_spectrum
    }

    #[must_use]
    pub fn evaluation_data_id(&self) -> &str {
        &self.evaluation_data_id
    }

    #[must_use]
    pub const fn parent(&self) -> NuclideId {
        self.parent
    }

    #[must_use]
    pub fn products(&self) -> &[NuclideId] {
        &self.products
    }

    /// Cross sections, absolute flux, provenance and product-order group yields.
    #[must_use]
    pub fn reaction_weighted(&self) -> &ReactionWeightedCollapse {
        &self.weighted
    }

    /// Append one disjoint channel per reacting group to a candidate data model.
    /// The source model is unchanged on success, refusal and cancellation.
    ///
    /// Prompt energy is explicitly supplied per group, excluding future product
    /// decay. Positive products need a represented nuclide or an explicit export;
    /// neither missing identities nor unknown decay laws silently become stable.
    /// Already represented products cannot also be designated as exports.
    ///
    /// This fixed processed-state bridge requires the exact source temperature
    /// as both model temperature bounds. It does not invent temperature dependence
    /// or broaden the validity range. The within-group spectrum is the one used
    /// for condensation; between-group flux amplitudes may subsequently change.
    pub fn append_to_burnup(
        &self,
        data: &BurnupData,
        prompt_energy_per_group: &[EnergyChannels],
        exported_products: &[NuclideId],
        limits: BurnupLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<BurnupData, EndfError> {
        checkpoint()?;
        let groups = self.weighted.flux.groups();
        if data.groups != *groups
            || prompt_energy_per_group.len() != groups.len()
            || data.temperature_min_k != self.weighted.metadata.temperature_k
            || data.temperature_max_k != self.weighted.metadata.temperature_k
        {
            return Err(invalid(
                "burnup requires the original grid, exact processed temperature and one prompt energy per group",
            ));
        }
        if groups.len() > limits.max_neutron_groups
            || data.nuclides.len() > limits.propagation.maximum_nodes
            || data.reactions.len() > limits.max_reactions
            || exported_products.len() > self.products.len()
        {
            return Err(budget("fission-yield burnup shape"));
        }
        for energy in prompt_energy_per_group {
            checkpoint()?;
            energy.validate()?;
        }
        let mut represented = BTreeMap::new();
        for (index, nuclide) in data.nuclides.iter().enumerate() {
            checkpoint()?;
            if represented.insert(nuclide.id, index).is_some() {
                return Err(invalid("burnup nuclide identities must be unique"));
            }
        }
        let parent = *represented
            .get(&self.parent)
            .ok_or_else(|| invalid("fission parent is absent from the burnup inventory"))?;
        let mut products_used = 0_usize;
        for reaction in &data.reactions {
            checkpoint()?;
            if reaction.parent == parent
                && matches!(
                    reaction.rate,
                    BurnupRate::Fission { .. } | BurnupRate::TabulatedFission { .. }
                )
            {
                return Err(invalid(
                    "existing fission channels would duplicate the imported reaction",
                ));
            }
            products_used = products_used
                .checked_add(reaction.products.len())
                .ok_or_else(|| budget("burnup product count"))?;
        }
        if products_used > limits.max_products {
            return Err(budget("burnup product count"));
        }
        let mut exports = BTreeSet::new();
        for &product in exported_products {
            checkpoint()?;
            if self.products.binary_search(&product).is_err()
                || represented.contains_key(&product)
                || !exports.insert(product)
            {
                return Err(invalid(
                    "exports must be unique evaluated products absent from the tracked inventory",
                ));
            }
        }
        for (product, row) in self.products.iter().zip(&self.weighted.yields) {
            checkpoint()?;
            if !represented.contains_key(product)
                && !exports.contains(product)
                && row.values.iter().flatten().any(|&value| value > 0.0)
            {
                return Err(invalid(
                    "a positive fission product is neither represented nor explicitly exported",
                ));
            }
        }
        let count = self
            .weighted
            .cross_sections_m2
            .iter()
            .filter(|&&value| value > 0.0)
            .count();
        if data
            .reactions
            .len()
            .checked_add(count)
            .is_none_or(|n| n > limits.max_reactions)
        {
            return Err(budget("fission-yield burnup reactions"));
        }
        let mut reactions = reserved(count, "fission-yield reaction memory")?;
        for (group, &sigma) in self.weighted.cross_sections_m2.iter().enumerate() {
            checkpoint()?;
            if sigma == 0.0 {
                continue;
            }
            let mut products = Vec::new();
            for (id, row) in self.products.iter().zip(&self.weighted.yields) {
                checkpoint()?;
                if let Some(&nuclide) = represented.get(id) {
                    let yield_per_event = row.values[group]
                        .ok_or_else(|| invalid("reacting group has no independent yield"))?;
                    if yield_per_event > 0.0 {
                        products_used = products_used
                            .checked_add(1)
                            .ok_or_else(|| budget("burnup product count"))?;
                        if products_used > limits.max_products {
                            return Err(budget("burnup product count"));
                        }
                        products.push(BurnupProduct {
                            nuclide,
                            yield_per_event,
                        });
                    }
                }
            }
            let mut cross_sections_m2 = reserved(groups.len(), "fission-yield channel memory")?;
            cross_sections_m2.resize(groups.len(), 0.0);
            cross_sections_m2[group] = sigma;
            reactions.push(BurnupReaction {
                data_id: format!(
                    "{} | {} | {} | group {}",
                    self.evaluation_data_id,
                    self.weighted.cross_section_data_id,
                    self.weighted.spectrum_data_id,
                    group
                ),
                parent,
                rate: BurnupRate::Fission { cross_sections_m2 },
                products,
                energy_per_event: prompt_energy_per_group[group],
            });
        }
        checkpoint()?;
        let mut next = data.clone();
        next.reactions.extend(reactions);
        checkpoint()?;
        Ok(next)
    }
}
