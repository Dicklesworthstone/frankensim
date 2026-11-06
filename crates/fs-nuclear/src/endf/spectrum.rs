//! Reaction-preserving collapse under an explicitly supplied neutron spectrum.
//!
//! For each group, Φ_g = ∫φ(E)dE and σ_g = ∫σ(E)φ(E)dE / Φ_g. The
//! immutable input φ is an absolute scalar flux density per unit energy,
//! in neutrons/(m² s J), not a normalized probability or a surface current.
//! Both pointwise laws may be histogram or linear and may contain paired
//! equal-energy points representing left and right limits at a discontinuity.
//!
//! Integration follows the union of the spectrum, cross-section and group
//! knots. Products of linear segments are integrated as positive quadratic
//! basis terms, not sampled at a representative energy. The returned group
//! pair preserves the supplied pointwise reaction integral up to floating-point
//! rounding. Each requested group must have positive integrated flux and lie
//! inside both input domains; no extrapolated or zero-flux coefficients are
//! invented. Reusing a collapsed row assumes the same within-group shape.
//!
//! This is energy condensation of supplied data, not a transport, resonance
//! reconstruction, self-shielding or spectrum-evolution calculation. The
//! scalar-flux weighting definition follows the OpenMC MGXS formulation:
//! <https://docs.openmc.org/en/v0.12.1/examples/mgxs-part-i.html>.

use std::sync::Arc;

use crate::cross_sections::{CrossSectionData, CrossSectionLimits, CrossSectionTable};
use crate::irradiation::{GroupFlux, NeutronGroups};
use crate::{Accumulator, FissionError, nonnegative, positive};

use super::{
    EndfCollapseLimits, EndfError, EndfInterpolation, EndfInterpolationRegion,
    PointwiseCrossSection, budget, convex,
};

/// Supplied absolute spectrum, with energies in joules and density in
/// neutrons/(m² s J). Region endpoints use the same one-based NBT convention
/// as the imported cross sections; the final endpoint must equal point count.
#[derive(Clone, Debug, PartialEq)]
pub struct FluxSpectrumData {
    pub data_id: String,
    pub energies_j: Vec<f64>,
    pub flux_density_m2_s_j: Vec<f64>,
    pub interpolation_regions: Vec<EndfInterpolationRegion>,
}

/// Construction bounds, checked before retaining the shared payload.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FluxSpectrumLimits {
    pub max_points: usize,
    pub max_interpolation_regions: usize,
}

impl Default for FluxSpectrumLimits {
    fn default() -> Self {
        Self {
            max_points: 1_000_000,
            max_interpolation_regions: 4096,
        }
    }
}

/// Identified immutable incident scalar spectrum. Clone shares its payload.
/// The spectrum is supplied independently of the cross-section temperature;
/// no temperature-dependent shape or nuclear-data accuracy is inferred.
#[derive(Clone, Debug, PartialEq)]
pub struct FluxSpectrum {
    pub(super) data_id: Arc<str>,
    energies_j: Arc<[f64]>,
    flux_density_m2_s_j: Arc<[f64]>,
    interpolation_regions: Arc<[EndfInterpolationRegion]>,
}

/// One flux-weighted microscopic row and the absolute group-integrated
/// spectrum that produced it. Both use `flux.groups()` in the original order.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectrumCollapse {
    pub cross_sections_m2: Vec<f64>,
    pub flux: GroupFlux,
}

impl FluxSpectrum {
    pub fn new(data: FluxSpectrumData, limits: FluxSpectrumLimits) -> Result<Self, EndfError> {
        Self::new_with_checkpoint(data, limits, || Ok(()))
    }

    pub fn new_with_checkpoint(
        data: FluxSpectrumData,
        limits: FluxSpectrumLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, EndfError> {
        checkpoint()?;
        let points = data.energies_j.len();
        let regions = data.interpolation_regions.len();
        if points > limits.max_points {
            return Err(budget("spectrum points"));
        }
        if regions > limits.max_interpolation_regions {
            return Err(budget("spectrum interpolation regions"));
        }
        if data.data_id.trim().is_empty()
            || points < 2
            || points != data.flux_density_m2_s_j.len()
            || regions == 0
            || regions >= points
        {
            return Err(invalid(
                "identified spectrum with matching point values and interpolation regions required",
            ));
        }
        for (index, (&energy, &density)) in data
            .energies_j
            .iter()
            .zip(&data.flux_density_m2_s_j)
            .enumerate()
        {
            if index % 64 == 0 {
                checkpoint()?;
            }
            if !nonnegative(energy)
                || !nonnegative(density)
                || (index > 0 && energy < data.energies_j[index - 1])
                || (index > 1
                    && energy == data.energies_j[index - 1]
                    && energy == data.energies_j[index - 2])
            {
                return Err(invalid(
                    "finite nonnegative spectrum values and ascending energies with at most paired discontinuities required",
                ));
            }
        }
        if data.energies_j[0] >= data.energies_j[points - 1] {
            return Err(invalid("nonempty spectrum energy domain required"));
        }
        let mut previous_end = 1;
        for region in &data.interpolation_regions {
            checkpoint()?;
            if region.end_point <= previous_end || region.end_point > points {
                return Err(invalid(
                    "ascending spectrum interpolation endpoints bounded by point count required",
                ));
            }
            previous_end = region.end_point;
        }
        if previous_end != points {
            return Err(invalid(
                "last spectrum interpolation endpoint must equal point count",
            ));
        }
        checkpoint()?;
        Ok(Self {
            data_id: data.data_id.into(),
            energies_j: data.energies_j.into(),
            flux_density_m2_s_j: data.flux_density_m2_s_j.into(),
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
    pub fn flux_density_m2_s_j(&self) -> &[f64] {
        &self.flux_density_m2_s_j
    }
    #[must_use]
    pub fn interpolation_regions(&self) -> &[EndfInterpolationRegion] {
        &self.interpolation_regions
    }

    /// No extrapolation; the right-hand value is returned at a discontinuity.
    pub fn density_at_energy(&self, energy_j: f64) -> Result<f64, EndfError> {
        self.check_domain(energy_j, energy_j)?;
        let cursor = Cursor::new(&self.energies_j, &self.interpolation_regions, energy_j);
        if self.energies_j[cursor.point] == energy_j {
            return Ok(self.flux_density_m2_s_j[cursor.point]);
        }
        self.segment_value(cursor, energy_j)
    }

    pub(super) fn check_domain(&self, lower: f64, upper: f64) -> Result<(), EndfError> {
        if !nonnegative(lower)
            || !nonnegative(upper)
            || lower < self.energies_j[0]
            || upper > self.energies_j[self.energies_j.len() - 1]
        {
            return Err(invalid("energy outside the supplied flux spectrum"));
        }
        Ok(())
    }

    pub(super) fn segment_value(&self, cursor: Cursor, energy: f64) -> Result<f64, EndfError> {
        let lower = self.flux_density_m2_s_j[cursor.point];
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
            self.flux_density_m2_s_j[cursor.point + 1],
            lower_weight,
            upper_weight,
        )
    }

    pub(super) fn integrate_groups(
        &self,
        groups: Arc<NeutronGroups>,
        limits: EndfCollapseLimits,
        visits: &mut usize,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<IntegratedSpectrum, EndfError> {
        checkpoint()?;
        if groups.len() > limits.max_neutron_groups {
            return Err(budget("spectrum collapse groups"));
        }
        let edges = groups.edges_j();
        self.check_domain(edges[0], edges[edges.len() - 1])?;
        let mut integrals = Vec::new();
        let mut values = Vec::new();
        integrals
            .try_reserve_exact(groups.len())
            .map_err(|_| budget("spectrum integral memory"))?;
        values
            .try_reserve_exact(groups.len())
            .map_err(|_| budget("spectrum flux memory"))?;
        let mut cursor = Cursor::new(&self.energies_j, &self.interpolation_regions, edges[0]);
        for bin in edges.windows(2) {
            checkpoint()?;
            let mut left = bin[0];
            let mut integral = PositiveIntegral::new();
            while left < bin[1] {
                visit(visits, limits, checkpoint)?;
                if cursor.advance(&self.energies_j, &self.interpolation_regions, left) {
                    continue;
                }
                let right = bin[1].min(self.energies_j[cursor.point + 1]);
                let width = right - left;
                integral.add_product(width, self.segment_value(cursor, left)?, 1.0, 2.0)?;
                integral.add_product(width, self.segment_value(cursor, right)?, 1.0, 2.0)?;
                left = right;
            }
            if !integral.nonzero {
                return Err(invalid(
                    "each collapsed group requires positive integrated spectral flux",
                ));
            }
            values.push(integral.value()?);
            integrals.push(integral);
        }
        checkpoint()?;
        Ok(IntegratedSpectrum {
            flux: GroupFlux::new(groups, values)?,
            integrals,
        })
    }
}

impl PointwiseCrossSection {
    /// Collapse under an absolute supplied spectrum. Unlike flat-in-energy
    /// averaging, this integrates σφ across both knot meshes. A zero cross
    /// section is valid; a zero-integral flux group has no defined coefficient.
    /// The work budget covers spectrum integration and the merged-grid pass.
    pub fn collapse_with_spectrum(
        &self,
        spectrum: &FluxSpectrum,
        groups: Arc<NeutronGroups>,
        limits: EndfCollapseLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<SpectrumCollapse, EndfError> {
        checkpoint()?;
        let edges = groups.edges_j();
        self.check_domain(edges[0], edges[edges.len() - 1])?;
        let mut visits = 0;
        let integrated = spectrum.integrate_groups(groups, limits, &mut visits, &mut checkpoint)?;
        let cross_sections_m2 =
            self.weighted_row(spectrum, &integrated, limits, &mut visits, &mut checkpoint)?;
        checkpoint()?;
        Ok(SpectrumCollapse {
            cross_sections_m2,
            flux: integrated.flux,
        })
    }

    fn weighted_row(
        &self,
        spectrum: &FluxSpectrum,
        integrated: &IntegratedSpectrum,
        limits: EndfCollapseLimits,
        visits: &mut usize,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<f64>, EndfError> {
        self.weighted_map(
            spectrum,
            integrated,
            limits,
            visits,
            checkpoint,
            |group, reaction| weighted_coefficient(reaction, integrated, group),
        )
    }

    /// Retain the same reaction integrals for multiple reaction-weighted
    /// quantities. Ordinary cross-section collapse retains no extra matrix.
    pub(super) fn weighted_integrals(
        &self,
        spectrum: &FluxSpectrum,
        integrated: &IntegratedSpectrum,
        limits: EndfCollapseLimits,
        visits: &mut usize,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<PositiveIntegral>, EndfError> {
        self.weighted_map(
            spectrum,
            integrated,
            limits,
            visits,
            checkpoint,
            |_, reaction| {
                reaction.value()?;
                Ok(reaction)
            },
        )
    }

    fn weighted_map<T>(
        &self,
        spectrum: &FluxSpectrum,
        integrated: &IntegratedSpectrum,
        limits: EndfCollapseLimits,
        visits: &mut usize,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
        mut value: impl FnMut(usize, PositiveIntegral) -> Result<T, EndfError>,
    ) -> Result<Vec<T>, EndfError> {
        checkpoint()?;
        let groups = integrated.flux.groups();
        let edges = groups.edges_j();
        self.check_domain(edges[0], edges[edges.len() - 1])?;
        let mut row = Vec::new();
        row.try_reserve_exact(groups.len())
            .map_err(|_| budget("spectrum-weighted row memory"))?;
        let mut cross_section_cursor =
            Cursor::new(&self.energies_j, &self.interpolation_regions, edges[0]);
        let mut spectrum_cursor = Cursor::new(
            &spectrum.energies_j,
            &spectrum.interpolation_regions,
            edges[0],
        );
        for (group, bin) in edges.windows(2).enumerate() {
            checkpoint()?;
            let mut left = bin[0];
            let mut reaction = PositiveIntegral::new();
            while left < bin[1] {
                visit(visits, limits, checkpoint)?;
                let advanced_cross_section = cross_section_cursor.advance(
                    &self.energies_j,
                    &self.interpolation_regions,
                    left,
                );
                let advanced_spectrum = spectrum_cursor.advance(
                    &spectrum.energies_j,
                    &spectrum.interpolation_regions,
                    left,
                );
                if advanced_cross_section || advanced_spectrum {
                    continue;
                }
                let right = bin[1]
                    .min(self.energies_j[cross_section_cursor.point + 1])
                    .min(spectrum.energies_j[spectrum_cursor.point + 1]);
                let f0 = spectrum.segment_value(spectrum_cursor, left)?;
                let f1 = spectrum.segment_value(spectrum_cursor, right)?;
                if f0 == 0.0 && f1 == 0.0 {
                    // The admitted spectrum is identically zero on this
                    // interval; no microscopic evaluation can change its rate.
                    left = right;
                    continue;
                }
                let s0 = self.segment_value(
                    cross_section_cursor.point,
                    cross_section_cursor.region,
                    left,
                )?;
                let s1 = self.segment_value(
                    cross_section_cursor.point,
                    cross_section_cursor.region,
                    right,
                )?;
                let width = right - left;
                // Exact integral of the product of two linear endpoint laws:
                // ΔE * (s0*f0/3 + s0*f1/6 + s1*f0/6 + s1*f1/3).
                // Histogram intervals enter with equal endpoint values.
                reaction.add_product(width, s0, f0, 3.0)?;
                reaction.add_product(width, s0, f1, 6.0)?;
                reaction.add_product(width, s1, f0, 6.0)?;
                reaction.add_product(width, s1, f1, 3.0)?;
                left = right;
            }
            row.push(value(group, reaction)?);
        }
        checkpoint()?;
        Ok(row)
    }
}

pub(super) fn weighted_coefficient(
    reaction: PositiveIntegral,
    integrated: &IntegratedSpectrum,
    group: usize,
) -> Result<f64, EndfError> {
    // Admit both the original integral and the returned coefficient/flux pair:
    // an unrepresentable positive hazard cannot become zero.
    reaction.value()?;
    let coefficient = reaction.ratio(integrated.integrals[group])?;
    let rate = coefficient * integrated.flux.values_m2_s()[group];
    if !nonnegative(rate) || (reaction.nonzero && rate == 0.0) {
        return Err(FissionError::Arithmetic.into());
    }
    Ok(coefficient)
}

/// Assemble processed temperature rows using ONE identified spectrum shape
/// for all rows. Flux is integrated once on its own mesh, independent of the
/// cross-section meshes, so the returned GroupFlux is shared exactly by every
/// row. A different shape requires another collapse; interpolation in Kelvin
/// does not infer spectrum evolution. `data_id` identifies the assembled law,
/// whose provenance must include the source curves and `spectrum.data_id()`.
pub fn collapse_temperature_table(
    data_id: String,
    curves: &[PointwiseCrossSection],
    spectrum: &FluxSpectrum,
    groups: Arc<NeutronGroups>,
    collapse_limits: EndfCollapseLimits,
    table_limits: CrossSectionLimits,
    mut checkpoint: impl FnMut() -> Result<(), FissionError>,
) -> Result<(CrossSectionTable, GroupFlux), EndfError> {
    checkpoint()?;
    let first = curves
        .first()
        .ok_or_else(|| invalid("at least one processed ENDF curve required"))?;
    if data_id.trim().is_empty() {
        return Err(invalid(
            "assembled spectrum-weighted source identity required",
        ));
    }
    if curves.len() > table_limits.max_temperatures {
        return Err(budget("cross-section temperatures"));
    }
    if groups.len() > table_limits.max_neutron_groups
        || groups.len() > collapse_limits.max_neutron_groups
    {
        return Err(budget("spectrum collapse groups"));
    }
    let values = curves
        .len()
        .checked_mul(groups.len())
        .ok_or_else(|| budget("cross-section table size"))?;
    if values > table_limits.max_values {
        return Err(budget("cross-section table values"));
    }
    let mut previous_temperature = 0.0;
    let edges = groups.edges_j();
    for curve in curves {
        checkpoint()?;
        let a = &first.metadata;
        let b = &curve.metadata;
        if b.nuclide != a.nuclide
            || curve.reaction != first.reaction
            || (
                b.library_id,
                b.library_version,
                b.library_release,
                b.material_revision,
            ) != (
                a.library_id,
                a.library_version,
                a.library_release,
                a.material_revision,
            )
            || b.temperature_k <= previous_temperature
        {
            return Err(invalid(
                "matching ENDF evaluation/reaction and ascending temperatures required",
            ));
        }
        curve.check_domain(edges[0], edges[edges.len() - 1])?;
        previous_temperature = b.temperature_k;
    }
    let mut visits = 0;
    let integrated = spectrum.integrate_groups(
        Arc::clone(&groups),
        collapse_limits,
        &mut visits,
        &mut checkpoint,
    )?;
    let mut temperatures = Vec::new();
    let mut rows = Vec::new();
    temperatures
        .try_reserve_exact(curves.len())
        .map_err(|_| budget("spectrum temperature memory"))?;
    rows.try_reserve_exact(curves.len())
        .map_err(|_| budget("spectrum row memory"))?;
    for curve in curves {
        checkpoint()?;
        temperatures.push(curve.metadata.temperature_k);
        rows.push(curve.weighted_row(
            spectrum,
            &integrated,
            collapse_limits,
            &mut visits,
            &mut checkpoint,
        )?);
    }
    let table = CrossSectionTable::new_with_checkpoint(
        CrossSectionData {
            data_id,
            nuclide: first.metadata.nuclide,
            reaction: first.reaction,
            groups,
            temperatures_k: temperatures,
            cross_sections_m2: rows,
        },
        table_limits,
        &mut checkpoint,
    )?;
    checkpoint()?;
    Ok((table, integrated.flux))
}

fn invalid(reason: &'static str) -> EndfError {
    FissionError::InvalidInput(reason).into()
}

pub(super) fn visit(
    visits: &mut usize,
    limits: EndfCollapseLimits,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(), EndfError> {
    checkpoint()?;
    if *visits >= limits.max_segment_visits {
        return Err(budget("spectrum collapse segment visits"));
    }
    *visits += 1;
    Ok(())
}

#[derive(Clone, Copy)]
pub(super) struct Cursor {
    pub(super) point: usize,
    pub(super) region: usize,
}

impl Cursor {
    pub(super) fn new(energies: &[f64], regions: &[EndfInterpolationRegion], energy: f64) -> Self {
        let point = energies.partition_point(|&e| e <= energy) - 1;
        let region = regions.partition_point(|r| r.end_point <= point + 1);
        Self { point, region }
    }

    // Called only while left is below the admitted group's upper edge. There
    // must therefore be a subsequent knot in each fully covering input domain.
    pub(super) fn advance(
        &mut self,
        energies: &[f64],
        regions: &[EndfInterpolationRegion],
        left: f64,
    ) -> bool {
        if energies[self.point + 1] > left {
            return false;
        }
        self.point += 1;
        while regions[self.region].end_point <= self.point + 1 {
            self.region += 1;
        }
        true
    }
}

pub(super) struct IntegratedSpectrum {
    pub(super) flux: GroupFlux,
    pub(super) integrals: Vec<PositiveIntegral>,
}

// Keep binary scale separate while integrating positive polynomial terms.
// Raw σ*φ can overflow although multiplying by a narrow ΔE makes the actual
// integral representable; likewise tiny positive basis terms can sum to a
// representable subnormal result. No signed slopes or cancellation are needed.
#[derive(Clone, Copy)]
pub(crate) struct PositiveIntegral {
    sum: Accumulator,
    exponent: i32,
    pub(super) nonzero: bool,
}

impl PositiveIntegral {
    pub(crate) const fn new() -> Self {
        Self {
            sum: Accumulator::new(0.0),
            exponent: 0,
            nonzero: false,
        }
    }

    pub(crate) fn add_product(
        &mut self,
        width: f64,
        a: f64,
        b: f64,
        divisor: f64,
    ) -> Result<(), EndfError> {
        if a == 0.0 || b == 0.0 {
            return Ok(());
        }
        let (w, we) = binary_parts(width);
        let (a, ae) = binary_parts(a);
        let (b, be) = binary_parts(b);
        let (term, te) = binary_parts(w * a * b / divisor);
        let exponent = we + ae + be + te;
        self.add_scaled_term(term, exponent)
    }

    /// Positive cubic basis term, with the same binary scaling as σφ. Keeping
    /// the yield factor separate avoids overflowing σ*φ*y before a small ΔE
    /// makes the physically required production integral representable.
    pub(crate) fn add_triple_product(
        &mut self,
        width: f64,
        a: f64,
        b: f64,
        c: f64,
        divisor: f64,
    ) -> Result<(), EndfError> {
        if a == 0.0 || b == 0.0 || c == 0.0 {
            return Ok(());
        }
        let (w, we) = binary_parts(width);
        let (a, ae) = binary_parts(a);
        let (b, be) = binary_parts(b);
        let (c, ce) = binary_parts(c);
        let (term, te) = binary_parts(w * a * b * c / divisor);
        self.add_scaled_term(term, we + ae + be + ce + te)
    }

    /// One additional independently scaled factor, used by the first-energy
    /// moment of a log-log decay spectrum. E² must not overflow before the
    /// density and normalization make its physical integral representable.
    pub(crate) fn add_four_product(
        &mut self,
        width: f64,
        a: f64,
        b: f64,
        c: f64,
        d: f64,
        divisor: f64,
    ) -> Result<(), EndfError> {
        if a == 0.0 || b == 0.0 || c == 0.0 || d == 0.0 {
            return Ok(());
        }
        let (w, we) = binary_parts(width);
        let (a, ae) = binary_parts(a);
        let (b, be) = binary_parts(b);
        let (c, ce) = binary_parts(c);
        let (d, de) = binary_parts(d);
        let (term, te) = binary_parts(w * a * b * c * d / divisor);
        self.add_scaled_term(term, we + ae + be + ce + de + te)
    }

    fn add_scaled_term(&mut self, term: f64, exponent: i32) -> Result<(), EndfError> {
        if !self.nonzero {
            self.sum = Accumulator::new(term);
            self.exponent = exponent;
            self.nonzero = true;
        } else if exponent > self.exponent {
            let shift = self.exponent - exponent;
            let high = scale_binary(self.sum.hi, shift);
            let low = scale_binary(self.sum.lo, shift);
            self.sum = Accumulator::new(high);
            self.sum.add(low)?;
            self.sum.add(term)?;
            self.exponent = exponent;
        } else {
            self.sum.add(scale_binary(term, exponent - self.exponent))?;
        }
        Ok(())
    }

    pub(crate) fn value(self) -> Result<f64, EndfError> {
        if !self.nonzero {
            return Ok(0.0);
        }
        let value = scale_binary(self.sum.value(), self.exponent);
        if !positive(value) {
            return Err(FissionError::Arithmetic.into());
        }
        Ok(value)
    }

    pub(crate) fn ratio(self, denominator: Self) -> Result<f64, EndfError> {
        if !self.nonzero {
            return Ok(0.0);
        }
        let (numerator, ne) = binary_parts(self.sum.value());
        let (denominator_value, de) = binary_parts(denominator.sum.value());
        let value = scale_binary(
            numerator / denominator_value,
            self.exponent - denominator.exponent + ne - de,
        );
        if !positive(value) {
            return Err(FissionError::Arithmetic.into());
        }
        Ok(value)
    }
}

// Exact finite decomposition; subnormal inputs are normalized by an exact 2^54
// multiplication. This local scaling is only for the positive integral above.
fn binary_parts(value: f64) -> (f64, i32) {
    if value == 0.0 {
        return (value, 0);
    }
    let mut bits = value.to_bits();
    let mut exponent = ((bits >> 52) & 0x7ff) as i32;
    let correction = if exponent == 0 {
        bits = (value * f64::from_bits(1077_u64 << 52)).to_bits();
        exponent = ((bits >> 52) & 0x7ff) as i32;
        -54
    } else {
        0
    };
    let mantissa = f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | (1023_u64 << 52));
    (mantissa, exponent - 1023 + correction)
}

fn scale_binary(value: f64, shift: i32) -> f64 {
    if value == 0.0 {
        return value;
    }
    let (mantissa, own_exponent) = binary_parts(value);
    let exponent = own_exponent + shift;
    if exponent > 1023 {
        f64::INFINITY.copysign(value)
    } else if exponent < -1075 {
        0.0_f64.copysign(value)
    } else if exponent >= -1022 {
        mantissa * f64::from_bits(((exponent + 1023) as u64) << 52)
    } else if exponent == -1075 {
        (mantissa * 0.5) * f64::from_bits(1)
    } else {
        mantissa * f64::from_bits(1_u64 << (exponent + 1074))
    }
}

#[cfg(test)]
mod tests;
