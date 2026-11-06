//! Reaction-rate-weighted incident-energy yields from supplied pointwise laws.
//!
//! For each group, y_g = integral(sigma*phi*y)/integral(sigma*phi). Together
//! with the original spectrum-weighted sigma_g and integrated flux Phi_g,
//! sigma_g*Phi_g*y_g preserves the pointwise production integral, up to
//! floating-point rounding. Flux weighting alone does not preserve production
//! when cross section and yield vary within a group.
//!
//! Every curve is nonnegative and dimensionless; multiplicities above one are
//! admitted without normalization. Physical meanings, nuclide identities and
//! evaluation compatibility belong to the importer or caller. Each output
//! retains its supplied source identity, and the cross-section metadata is
//! copied unchanged. These declarations are not a data-quality certificate.
//!
//! Histogram/linear segments are integrated on the union of all three knot
//! meshes and group edges. Positive cubic basis terms avoid signed slopes,
//! cancellation and representative-energy sampling. Shared binary scaling
//! keeps narrow-interval products representable. No extrapolation is allowed.
//! A zero reaction denominator has no defined yield and returns None; it never
//! silently acquires an arbitrary multiplicity. A zero-flux group is refused
//! under the existing cross-section collapse contract.

use std::sync::Arc;

use crate::cross_sections::NeutronReaction;
use crate::irradiation::{GroupFlux, NeutronGroups};
use crate::{FissionError, nonnegative};

use super::spectrum::{Cursor, FluxSpectrum, PositiveIntegral, visit, weighted_coefficient};
use super::{
    EndfCollapseLimits, EndfError, EndfInterpolation, EndfInterpolationRegion, EndfMetadata,
    PointwiseCrossSection, budget, convex,
};

/// Authored nonnegative dimensionless values at incident energies in joules.
/// Paired equal energies give left and right limits of a discontinuity.
#[derive(Clone, Debug, PartialEq)]
pub struct TabulatedYieldData {
    pub data_id: String,
    pub energies_j: Vec<f64>,
    pub values: Vec<f64>,
    pub interpolation_regions: Vec<EndfInterpolationRegion>,
}

/// Construction bounds, admitted before retaining the immutable payload.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TabulatedYieldLimits {
    pub max_points: usize,
    pub max_interpolation_regions: usize,
}

impl Default for TabulatedYieldLimits {
    fn default() -> Self {
        Self {
            max_points: 1_000_000,
            max_interpolation_regions: 4096,
        }
    }
}

/// Reusable incident-energy yield. Clone shares the admitted source payload.
#[derive(Clone, Debug, PartialEq)]
pub struct TabulatedYield {
    data_id: Arc<str>,
    energies_j: Arc<[f64]>,
    values: Arc<[f64]>,
    interpolation_regions: Arc<[EndfInterpolationRegion]>,
}

impl TabulatedYield {
    pub fn new(data: TabulatedYieldData, limits: TabulatedYieldLimits) -> Result<Self, EndfError> {
        Self::new_with_checkpoint(data, limits, || Ok(()))
    }

    pub fn new_with_checkpoint(
        data: TabulatedYieldData,
        limits: TabulatedYieldLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, EndfError> {
        checkpoint()?;
        let points = data.energies_j.len();
        let regions = data.interpolation_regions.len();
        if points > limits.max_points {
            return Err(budget("yield curve points"));
        }
        if regions > limits.max_interpolation_regions {
            return Err(budget("yield interpolation regions"));
        }
        if data.data_id.trim().is_empty()
            || points < 2
            || points != data.values.len()
            || regions == 0
            || regions >= points
        {
            return Err(invalid(
                "identified dimensionless yield with matching point values and interpolation regions required",
            ));
        }
        for (index, (&energy, &value)) in data.energies_j.iter().zip(&data.values).enumerate() {
            if index % 64 == 0 {
                checkpoint()?;
            }
            if !nonnegative(energy)
                || !nonnegative(value)
                || (index > 0 && energy < data.energies_j[index - 1])
                || (index > 1
                    && energy == data.energies_j[index - 1]
                    && energy == data.energies_j[index - 2])
            {
                return Err(invalid(
                    "finite nonnegative yields and ascending energies with at most paired discontinuities required",
                ));
            }
        }
        if data.energies_j[0] >= data.energies_j[points - 1] {
            return Err(invalid("nonempty yield energy domain required"));
        }
        let mut previous_end = 1;
        for region in &data.interpolation_regions {
            checkpoint()?;
            if region.end_point <= previous_end || region.end_point > points {
                return Err(invalid(
                    "ascending yield interpolation endpoints bounded by point count required",
                ));
            }
            previous_end = region.end_point;
        }
        if previous_end != points {
            return Err(invalid(
                "last yield interpolation endpoint must equal point count",
            ));
        }
        checkpoint()?;
        Ok(Self {
            data_id: data.data_id.into(),
            energies_j: data.energies_j.into(),
            values: data.values.into(),
            interpolation_regions: data.interpolation_regions.into(),
        })
    }

    #[must_use]
    pub fn data_id(&self) -> &str {
        &self.data_id
    }
    #[must_use]
    pub fn energies_j(&self) -> &[f64] {
        &self.energies_j
    }
    #[must_use]
    pub fn values(&self) -> &[f64] {
        &self.values
    }
    #[must_use]
    pub fn interpolation_regions(&self) -> &[EndfInterpolationRegion] {
        &self.interpolation_regions
    }

    /// Evaluate without extrapolation. Exact jump knots return the right limit.
    pub fn value_at_energy(&self, energy_j: f64) -> Result<f64, EndfError> {
        self.check_domain(energy_j, energy_j)?;
        let cursor = Cursor::new(&self.energies_j, &self.interpolation_regions, energy_j);
        if self.energies_j[cursor.point] == energy_j {
            return Ok(self.values[cursor.point]);
        }
        self.segment_value(cursor, energy_j)
    }

    fn check_domain(&self, lower: f64, upper: f64) -> Result<(), EndfError> {
        if !nonnegative(lower)
            || !nonnegative(upper)
            || lower < self.energies_j[0]
            || upper > self.energies_j[self.energies_j.len() - 1]
        {
            return Err(invalid("energy outside the supplied yield curve"));
        }
        Ok(())
    }

    fn segment_value(&self, cursor: Cursor, energy: f64) -> Result<f64, EndfError> {
        let lower = self.values[cursor.point];
        if self.interpolation_regions[cursor.region].interpolation == EndfInterpolation::Histogram {
            return Ok(lower);
        }
        let low_energy = self.energies_j[cursor.point];
        let high_energy = self.energies_j[cursor.point + 1];
        let width = high_energy - low_energy;
        let lower_weight = (high_energy - energy) / width;
        let upper_weight = (energy - low_energy) / width;
        if (lower_weight == 0.0 && energy != high_energy)
            || (upper_weight == 0.0 && energy != low_energy)
        {
            return Err(FissionError::Arithmetic.into());
        }
        convex(
            lower,
            self.values[cursor.point + 1],
            lower_weight,
            upper_weight,
        )
    }
}

/// One segment-visit budget covers the shared spectrum/reaction passes and
/// every yield numerator. Output bounds include all curve/group combinations.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReactionWeightedLimits {
    pub collapse: EndfCollapseLimits,
    pub max_curves: usize,
    pub max_values: usize,
}

impl Default for ReactionWeightedLimits {
    fn default() -> Self {
        Self {
            collapse: EndfCollapseLimits::default(),
            max_curves: 4096,
            max_values: 1_000_000,
        }
    }
}

/// One curve's yields in the original neutron-group order. None means that
/// group's integral(sigma*phi) is zero. Some(0) is a defined zero production
/// yield in a group whose reaction rate is positive.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionWeightedYield {
    pub data_id: Arc<str>,
    pub values: Vec<Option<f64>>,
}

/// Original spectrum-weighted cross sections and flux, plus reaction-weighted
/// yields in the input curve order. Reusing the coefficients assumes the same
/// within-group spectral shape and processed cross-section temperature.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionWeightedCollapse {
    pub cross_sections_m2: Vec<f64>,
    pub flux: GroupFlux,
    pub yields: Vec<ReactionWeightedYield>,
    pub cross_section_data_id: Arc<str>,
    pub spectrum_data_id: Arc<str>,
    pub metadata: EndfMetadata,
    pub reaction: NeutronReaction,
}

/// Collapse several dimensionless yields against one physical reaction and
/// spectrum. The flux and reaction denominator are integrated only once.
/// Each numerator uses its own union mesh, avoiding a potentially large common
/// knot allocation and repeated interpolation on unrelated yield curves.
pub fn collapse_reaction_weighted(
    cross_section: &PointwiseCrossSection,
    spectrum: &FluxSpectrum,
    curves: &[TabulatedYield],
    groups: Arc<NeutronGroups>,
    limits: ReactionWeightedLimits,
    mut checkpoint: impl FnMut() -> Result<(), FissionError>,
) -> Result<ReactionWeightedCollapse, EndfError> {
    checkpoint()?;
    if curves.is_empty() {
        return Err(invalid("at least one incident-energy yield curve required"));
    }
    if curves.len() > limits.max_curves {
        return Err(budget("reaction-weighted yield curves"));
    }
    if groups.len() > limits.collapse.max_neutron_groups {
        return Err(budget("reaction-weighted neutron groups"));
    }
    let values = curves
        .len()
        .checked_mul(groups.len())
        .ok_or_else(|| budget("reaction-weighted values"))?;
    if values > limits.max_values {
        return Err(budget("reaction-weighted values"));
    }
    let edges = groups.edges_j();
    cross_section.check_domain(edges[0], edges[edges.len() - 1])?;
    for curve in curves {
        checkpoint()?;
        curve.check_domain(edges[0], edges[edges.len() - 1])?;
    }
    let mut visits = 0;
    let integrated =
        spectrum.integrate_groups(groups, limits.collapse, &mut visits, &mut checkpoint)?;
    let reactions = cross_section.weighted_integrals(
        spectrum,
        &integrated,
        limits.collapse,
        &mut visits,
        &mut checkpoint,
    )?;
    let mut cross_sections_m2 =
        reserved(reactions.len(), "reaction-weighted cross-section memory")?;
    for (group, &reaction) in reactions.iter().enumerate() {
        checkpoint()?;
        cross_sections_m2.push(weighted_coefficient(reaction, &integrated, group)?);
    }
    let mut yields = reserved(curves.len(), "reaction-weighted row memory")?;
    for curve in curves {
        checkpoint()?;
        let values = collapse_yield(
            curve,
            cross_section,
            spectrum,
            &integrated.flux,
            &cross_sections_m2,
            &reactions,
            limits.collapse,
            &mut visits,
            &mut checkpoint,
        )?;
        yields.push(ReactionWeightedYield {
            data_id: Arc::clone(&curve.data_id),
            values,
        });
    }
    checkpoint()?;
    Ok(ReactionWeightedCollapse {
        cross_sections_m2,
        flux: integrated.flux,
        yields,
        cross_section_data_id: Arc::clone(&cross_section.data_id),
        spectrum_data_id: Arc::clone(&spectrum.data_id),
        metadata: cross_section.metadata().clone(),
        reaction: cross_section.reaction(),
    })
}

#[allow(clippy::too_many_arguments)]
fn collapse_yield(
    curve: &TabulatedYield,
    cross_section: &PointwiseCrossSection,
    spectrum: &FluxSpectrum,
    flux: &GroupFlux,
    cross_sections: &[f64],
    reactions: &[PositiveIntegral],
    limits: EndfCollapseLimits,
    visits: &mut usize,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<Vec<Option<f64>>, EndfError> {
    let edges = flux.groups().edges_j();
    let mut values = reserved(flux.groups().len(), "reaction-weighted yield memory")?;
    let mut sigma_cursor = Cursor::new(
        cross_section.energies_j(),
        cross_section.interpolation_regions(),
        edges[0],
    );
    let mut flux_cursor = Cursor::new(
        spectrum.energies_j(),
        spectrum.interpolation_regions(),
        edges[0],
    );
    let mut yield_cursor = Cursor::new(&curve.energies_j, &curve.interpolation_regions, edges[0]);
    for (group, bin) in edges.windows(2).enumerate() {
        checkpoint()?;
        if !reactions[group].nonzero {
            values.push(None);
            continue;
        }
        let mut left = bin[0];
        let mut production = PositiveIntegral::new();
        while left < bin[1] {
            visit(visits, limits, checkpoint)?;
            let sigma_advanced = sigma_cursor.advance(
                cross_section.energies_j(),
                cross_section.interpolation_regions(),
                left,
            );
            let flux_advanced = flux_cursor.advance(
                spectrum.energies_j(),
                spectrum.interpolation_regions(),
                left,
            );
            let yield_advanced =
                yield_cursor.advance(&curve.energies_j, &curve.interpolation_regions, left);
            if sigma_advanced || flux_advanced || yield_advanced {
                continue;
            }
            let right = bin[1]
                .min(cross_section.energies_j()[sigma_cursor.point + 1])
                .min(spectrum.energies_j()[flux_cursor.point + 1])
                .min(curve.energies_j[yield_cursor.point + 1]);
            let phi = [
                spectrum.segment_value(flux_cursor, left)?,
                spectrum.segment_value(flux_cursor, right)?,
            ];
            if phi == [0.0, 0.0] {
                left = right;
                continue;
            }
            let sigma = [
                cross_section.segment_value(sigma_cursor.point, sigma_cursor.region, left)?,
                cross_section.segment_value(sigma_cursor.point, sigma_cursor.region, right)?,
            ];
            if sigma == [0.0, 0.0] {
                left = right;
                continue;
            }
            let yield_ = [
                curve.segment_value(yield_cursor, left)?,
                curve.segment_value(yield_cursor, right)?,
            ];
            let width = right - left;
            // Each factor is (1-u)*left + u*right. The two unmixed cubic
            // terms integrate to 1/4; all six mixed terms integrate to 1/12.
            // Histogram intervals have equal endpoint values in this basis.
            for (i, &s) in sigma.iter().enumerate() {
                for (j, &f) in phi.iter().enumerate() {
                    for (k, &y) in yield_.iter().enumerate() {
                        let divisor = if i == j && j == k { 4.0 } else { 12.0 };
                        production.add_triple_product(width, s, f, y, divisor)?;
                    }
                }
            }
            left = right;
        }
        production.value()?;
        let yield_ = production.ratio(reactions[group])?;
        let represented = (cross_sections[group] * flux.values_m2_s()[group]) * yield_;
        if !nonnegative(represented) || (production.nonzero && represented == 0.0) {
            return Err(FissionError::Arithmetic.into());
        }
        values.push(Some(yield_));
    }
    Ok(values)
}

fn reserved<T>(count: usize, name: &'static str) -> Result<Vec<T>, EndfError> {
    let mut values = Vec::new();
    values.try_reserve_exact(count).map_err(|_| budget(name))?;
    Ok(values)
}

fn invalid(reason: &'static str) -> EndfError {
    FissionError::InvalidInput(reason).into()
}

#[cfg(test)]
mod tests;
