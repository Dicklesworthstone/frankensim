//! Free-gas Doppler heating of complete processed reaction curves.
//!
//! This is an explicit opt-in to the free-gas model: the caller must establish
//! that the input at T0 has the appropriate free-gas temperature treatment.
//! It is not valid to reinterpret arbitrary bound-atom or shielded data this
//! way. Only heating, T >= T0, is defined; no deconvolution is attempted.
//!
//! With theta = k_B (T-T0)/AWR, x = sqrt(E'/theta), y = sqrt(E/theta),
//! the broadened reaction rate g(E) = sqrt(E) sigma(E) is
//! `g_T = integral_0^infinity g(theta*x*x) x/(sqrt(pi)*y)
//!       [exp(-(x-y)^2)-exp(-(x+y)^2)] dx`.
//! This positive kernel integrates to one, so a 1/v cross section is invariant.
//! The negative-speed image term is essential at low energy. The temperature
//! increment, target/neutron mass ratio and relative-speed weighting follow
//! NJOY's BROADR theory (not a Gaussian convolution in incident energy):
//! <https://github.com/njoy/NJOY2016-manual/blob/master/broadx.tex>.
//!
//! Outside the finite source interval, BOTH tail laws must be supplied. They
//! are model assumptions, not evaluated data. Quadrature follows every original
//! histogram/linear segment and jump. The finite kernel cutoff has a separate
//! analytic Maxwell-tail bound; Gauss 8/16 differences estimate integration
//! error but are not interval certificates. A caller-supplied output grid adds
//! its own interpolation error, which these knot estimates do not measure.

use std::sync::Arc;

use fs_math::det;

use super::{EndfError, EndfInterpolation, EndfInterpolationRegion, PointwiseCrossSection, budget};
use crate::{Accumulator, FissionError, nonnegative, positive};

const BOLTZMANN_J_K: f64 = 1.380_649e-23;
const INV_SQRT_PI: f64 = 0.564_189_583_547_756_3;
const MAX_DEPTH: u8 = 32;

/// An authored continuation beyond one endpoint of the supplied energy grid.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DopplerTail {
    /// The reaction is zero throughout this exterior interval.
    Zero,
    /// sigma(E) = sigma(endpoint) sqrt(E_endpoint/E).
    /// A positive endpoint energy is required when the exterior interval exists.
    InverseSpeed,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DopplerTails {
    pub below: DopplerTail,
    pub above: DopplerTail,
}

/// Tolerances concern point evaluations, not interpolation between output knots.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DopplerAccuracy {
    pub relative_tolerance: f64,
    pub absolute_tolerance_m2: f64,
    /// Retain |sqrt(E'/theta)-sqrt(E/theta)| <= this number, in [4, 26].
    pub kernel_half_widths: f64,
}

impl Default for DopplerAccuracy {
    fn default() -> Self {
        Self {
            relative_tolerance: 1.0e-8,
            absolute_tolerance_m2: 0.0,
            kernel_half_widths: 8.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DopplerLimits {
    pub max_source_points: usize,
    pub max_output_points: usize,
    /// Shared by all output knots in one tabulation, including source jumps.
    pub max_segment_visits: usize,
    /// Shared by all quadrature levels and output knots in one tabulation.
    pub max_kernel_evaluations: usize,
}

impl Default for DopplerLimits {
    fn default() -> Self {
        Self {
            max_source_points: 1_000_000,
            max_output_points: 100_000,
            max_segment_visits: 4_000_000,
            max_kernel_evaluations: 16_000_000,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DopplerSample {
    pub cross_section_m2: f64,
    pub estimated_quadrature_error_m2: f64,
    /// Bounds omitted thermal velocities under the declared exterior laws.
    /// It does not bound uncertainty in those physical assumptions.
    pub omitted_kernel_tail_bound_m2: f64,
}

/// The original native curve can go directly to either existing collapse API.
/// Its linear interpolation is the caller's grid approximation. Refining that
/// grid is separate from tightening the quadrature at the retained knots.
#[derive(Clone, Debug, PartialEq)]
pub struct DopplerTabulation {
    pub curve: PointwiseCrossSection,
    pub source_data_id: Arc<str>,
    pub source_temperature_k: f64,
    pub tails: DopplerTails,
    pub maximum_estimated_quadrature_error_m2: f64,
    pub maximum_omitted_kernel_tail_bound_m2: f64,
}

/// Immutable prepared heating operator; cloning shares the original source.
/// The rate bound is prepared once, not rescanned for each requested energy.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedDopplerBroadening {
    source: PointwiseCrossSection,
    target_temperature_k: f64,
    sqrt_theta_j: f64,
    rate_bound_m2_sqrt_j: f64,
    tails: DopplerTails,
}

impl PreparedDopplerBroadening {
    pub fn new(
        source: &PointwiseCrossSection,
        target_temperature_k: f64,
        tails: DopplerTails,
        limits: DopplerLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, EndfError> {
        checkpoint()?;
        if !positive(target_temperature_k) || target_temperature_k < source.metadata.temperature_k {
            return Err(invalid(
                "Doppler heating requires T >= the processed source temperature",
            ));
        }
        if source.energies_j.len() > limits.max_source_points {
            return Err(budget("Doppler source points"));
        }
        let increment = target_temperature_k - source.metadata.temperature_k;
        let sqrt_theta_j = if increment == 0.0 {
            0.0
        } else {
            let theta = (BOLTZMANN_J_K * increment) / source.metadata.atomic_weight_ratio;
            if !positive(theta) {
                return Err(FissionError::Arithmetic.into());
            }
            det::sqrt(theta)
        };
        let mut bound = 0.0_f64;
        if increment > 0.0 {
            // On each original linear/histogram interval sigma is bounded by
            // its endpoints. The sqrt(E) factor increases with E. Inverse-speed
            // exterior laws have constant g, already included in this bound.
            for index in 1..source.energies_j.len() {
                if index % 64 == 0 {
                    checkpoint()?;
                }
                let sigma =
                    source.cross_sections_m2[index - 1].max(source.cross_sections_m2[index]);
                if source.energies_j[index] == 0.0 {
                    continue;
                }
                let rate = det::sqrt(source.energies_j[index]).next_up() * sigma;
                if !nonnegative(rate) || (sigma > 0.0 && rate == 0.0) {
                    return Err(FissionError::Arithmetic.into());
                }
                if rate > 0.0 {
                    bound = bound.max(rate.next_up());
                }
            }
            if !nonnegative(bound) {
                return Err(FissionError::Arithmetic.into());
            }
        }
        checkpoint()?;
        Ok(Self {
            source: source.clone(),
            target_temperature_k,
            sqrt_theta_j,
            rate_bound_m2_sqrt_j: bound,
            tails,
        })
    }

    #[must_use]
    pub const fn source(&self) -> &PointwiseCrossSection {
        &self.source
    }

    #[must_use]
    pub const fn target_temperature_k(&self) -> f64 {
        self.target_temperature_k
    }

    /// Positive incident energy is required for heating: sigma generally has a
    /// physical 1/v singularity at zero. At unchanged temperature the original
    /// domain, discontinuities and bitwise point evaluation are preserved.
    pub fn at_energy(
        &self,
        energy_j: f64,
        accuracy: DopplerAccuracy,
        limits: DopplerLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DopplerSample, EndfError> {
        self.evaluate(energy_j, accuracy, &mut Work::new(limits), &mut checkpoint)
    }

    /// Sample onto an explicitly supplied, strictly increasing positive grid.
    /// No claim is made that this grid resolves the broadened wings. The new
    /// data ID identifies that approximation; the source ID is retained beside
    /// it. Identity heating requires the original grid and preserves its laws,
    /// including duplicate-energy jumps, instead of silently relabeling them.
    pub fn tabulate(
        &self,
        data_id: String,
        energies_j: Vec<f64>,
        accuracy: DopplerAccuracy,
        limits: DopplerLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DopplerTabulation, EndfError> {
        checkpoint()?;
        if data_id.trim().is_empty() || energies_j.len() < 2 {
            return Err(invalid(
                "identified Doppler output with at least two energy knots required",
            ));
        }
        if energies_j.len() > limits.max_output_points {
            return Err(budget("Doppler output points"));
        }
        let mut output = DopplerTabulation {
            curve: self.source.clone(),
            source_data_id: Arc::clone(&self.source.data_id),
            source_temperature_k: self.source.metadata.temperature_k,
            tails: self.tails,
            maximum_estimated_quadrature_error_m2: 0.0,
            maximum_omitted_kernel_tail_bound_m2: 0.0,
        };
        output.curve.data_id = data_id.into();
        if self.sqrt_theta_j == 0.0 {
            if energies_j.as_slice() != self.source.energies_j.as_ref() {
                return Err(invalid(
                    "identity Doppler tabulation retains the original energy grid",
                ));
            }
            checkpoint()?;
            return Ok(output);
        }
        let mut previous = 0.0;
        for (index, &energy) in energies_j.iter().enumerate() {
            if index % 64 == 0 {
                checkpoint()?;
            }
            if !positive(energy) || energy <= previous {
                return Err(invalid(
                    "positive strictly increasing Doppler output energies required",
                ));
            }
            previous = energy;
        }
        let mut values = Vec::new();
        values
            .try_reserve_exact(energies_j.len())
            .map_err(|_| budget("Doppler output memory"))?;
        let mut work = Work::new(limits);
        for &energy in &energies_j {
            let sample = self.evaluate(energy, accuracy, &mut work, &mut checkpoint)?;
            values.push(sample.cross_section_m2);
            output.maximum_estimated_quadrature_error_m2 = output
                .maximum_estimated_quadrature_error_m2
                .max(sample.estimated_quadrature_error_m2);
            output.maximum_omitted_kernel_tail_bound_m2 = output
                .maximum_omitted_kernel_tail_bound_m2
                .max(sample.omitted_kernel_tail_bound_m2);
        }
        output.curve.interpolation_regions = vec![EndfInterpolationRegion {
            end_point: energies_j.len(),
            interpolation: EndfInterpolation::Linear,
        }]
        .into();
        output.curve.energies_j = energies_j.into();
        output.curve.cross_sections_m2 = values.into();
        output.curve.metadata.temperature_k = self.target_temperature_k;
        output.curve.metadata.derivation_flag = 1;
        checkpoint()?;
        Ok(output)
    }

    fn evaluate(
        &self,
        energy: f64,
        accuracy: DopplerAccuracy,
        work: &mut Work,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DopplerSample, EndfError> {
        checkpoint()?;
        if self.sqrt_theta_j == 0.0 {
            return Ok(DopplerSample {
                cross_section_m2: self.source.cross_section_at_energy(energy)?,
                estimated_quadrature_error_m2: 0.0,
                omitted_kernel_tail_bound_m2: 0.0,
            });
        }
        if !positive(energy)
            || !nonnegative(accuracy.relative_tolerance)
            || !nonnegative(accuracy.absolute_tolerance_m2)
            || accuracy.relative_tolerance == 0.0 && accuracy.absolute_tolerance_m2 == 0.0
            || !accuracy.kernel_half_widths.is_finite()
            || !(4.0..=26.0).contains(&accuracy.kernel_half_widths)
        {
            return Err(invalid(
                "positive energy, nonzero finite Doppler tolerances and kernel radius in [4,26] required",
            ));
        }
        if self.rate_bound_m2_sqrt_j == 0.0 {
            return Ok(DopplerSample {
                cross_section_m2: 0.0,
                estimated_quadrature_error_m2: 0.0,
                omitted_kernel_tail_bound_m2: 0.0,
            });
        }
        let sqrt_energy = det::sqrt(energy);
        let y = sqrt_energy / self.sqrt_theta_j;
        let scale = self.rate_bound_m2_sqrt_j / sqrt_energy;
        if !positive(y) || !positive(scale) {
            return Err(FissionError::Arithmetic.into());
        }
        let radius = accuracy.kernel_half_widths;
        let start = (-y).max(-radius);
        let end = radius;
        let parameters = Kernel {
            owner: self,
            energy,
            sqrt_energy,
            y,
            absolute_tolerance: accuracy.absolute_tolerance_m2 / scale,
            relative_tolerance: accuracy.relative_tolerance,
            width: end - start,
        };
        let mut integral = Accumulator::new(0.0);
        let mut error = Accumulator::new(0.0);
        let mut left = start;
        let energies = &self.source.energies_j;
        let first = parameters.coordinate(energies[0]);
        if left < first {
            let right = end.min(first);
            let rate = self.exterior_rate(self.tails.below, 0);
            parameters.integrate_piece(
                left,
                right,
                Law::ConstantRate(rate),
                work,
                checkpoint,
                &mut integral,
                &mut error,
            )?;
            left = right;
        }
        let mut segment = energies
            .partition_point(|&e| parameters.coordinate(e) <= left)
            .saturating_sub(1);
        let mut region = self
            .source
            .interpolation_regions
            .partition_point(|r| r.end_point <= segment + 1);
        while left < end && segment + 1 < energies.len() {
            work.visit(checkpoint)?;
            let right = parameters.coordinate(energies[segment + 1]).min(end);
            if right > left {
                let law = Law::Segment { segment, region };
                parameters.integrate_piece(
                    left,
                    right,
                    law,
                    work,
                    checkpoint,
                    &mut integral,
                    &mut error,
                )?;
                left = right;
            }
            segment += 1;
            while region < self.source.interpolation_regions.len()
                && self.source.interpolation_regions[region].end_point <= segment + 1
            {
                region += 1;
            }
        }
        if left < end {
            let rate = self.exterior_rate(self.tails.above, energies.len() - 1);
            parameters.integrate_piece(
                left,
                end,
                Law::ConstantRate(rate),
                work,
                checkpoint,
                &mut integral,
                &mut error,
            )?;
        }
        // If |relative_speed - incident_speed| > L, the sampled thermal
        // velocity has magnitude > L. Its 3D Maxwell tail therefore bounds
        // all omitted normalized rate, independently of incident speed.
        let tail = (det::erfc(radius) + 2.0 * radius * det::exp(-radius * radius) * INV_SQRT_PI)
            * (1.0 + 64.0 * f64::EPSILON);
        let value = integral.value();
        let estimate = error.value();
        let allowed = parameters.absolute_tolerance + accuracy.relative_tolerance * value;
        if estimate + tail > allowed {
            return Err(budget(
                "Doppler quadrature and omitted kernel tail tolerance",
            ));
        }
        let sample = DopplerSample {
            cross_section_m2: scale * value,
            estimated_quadrature_error_m2: scale * estimate,
            omitted_kernel_tail_bound_m2: (scale * tail).next_up(),
        };
        if !nonnegative(sample.cross_section_m2)
            || !nonnegative(sample.estimated_quadrature_error_m2)
            || !nonnegative(sample.omitted_kernel_tail_bound_m2)
            || value > 0.0 && sample.cross_section_m2 == 0.0
        {
            return Err(FissionError::Arithmetic.into());
        }
        checkpoint()?;
        Ok(sample)
    }

    fn exterior_rate(&self, law: DopplerTail, index: usize) -> f64 {
        match law {
            DopplerTail::Zero => 0.0,
            DopplerTail::InverseSpeed => {
                det::sqrt(self.source.energies_j[index]) * self.source.cross_sections_m2[index]
                    / self.rate_bound_m2_sqrt_j
            }
        }
    }
}

#[derive(Clone, Copy)]
enum Law {
    ConstantRate(f64),
    Segment { segment: usize, region: usize },
}

struct Kernel<'a> {
    owner: &'a PreparedDopplerBroadening,
    energy: f64,
    sqrt_energy: f64,
    y: f64,
    absolute_tolerance: f64,
    relative_tolerance: f64,
    width: f64,
}

impl Kernel<'_> {
    fn coordinate(&self, energy: f64) -> f64 {
        // Rationalized sqrt difference retains tiny thermal windows at high E.
        ((energy - self.energy) / (det::sqrt(energy) + self.sqrt_energy)) / self.owner.sqrt_theta_j
    }

    fn value(&self, t: f64, law: Law) -> Result<f64, EndfError> {
        let x = self.y + t;
        if x <= 0.0 {
            return Ok(0.0);
        }
        let rate = match law {
            Law::ConstantRate(rate) => rate,
            Law::Segment { segment, region } => {
                let source = &self.owner.source;
                let sigma = if source.interpolation_regions[region].interpolation
                    == EndfInterpolation::Histogram
                    || source.cross_sections_m2[segment] == source.cross_sections_m2[segment + 1]
                {
                    source.cross_sections_m2[segment]
                } else {
                    let low = source.energies_j[segment];
                    let high = source.energies_j[segment + 1];
                    let root = self.sqrt_energy + self.owner.sqrt_theta_j * t;
                    let lower_weight = self.owner.sqrt_theta_j
                        * (self.coordinate(high) - t)
                        * (det::sqrt(high) + root)
                        / (high - low);
                    let upper_weight = self.owner.sqrt_theta_j
                        * (t - self.coordinate(low))
                        * (root + det::sqrt(low))
                        / (high - low);
                    let roundoff = 64.0 * f64::EPSILON;
                    if !lower_weight.is_finite()
                        || !upper_weight.is_finite()
                        || lower_weight < -roundoff
                        || upper_weight < -roundoff
                        || lower_weight > 1.0 + roundoff
                        || upper_weight > 1.0 + roundoff
                    {
                        return Err(FissionError::Arithmetic.into());
                    }
                    super::convex(
                        source.cross_sections_m2[segment],
                        source.cross_sections_m2[segment + 1],
                        lower_weight.clamp(0.0, 1.0),
                        upper_weight.clamp(0.0, 1.0),
                    )?
                };
                ((self.sqrt_energy + self.owner.sqrt_theta_j * t) * sigma)
                    / self.owner.rate_bound_m2_sqrt_j
            }
        };
        // x/y * (1-exp(-4xy)) = 4x^2 * (1-exp(-z))/z. The latter has
        // a regular limit at y=0, with no cancellation between two Gaussians.
        let z = 4.0 * self.y * x;
        let factor = if z < 1.0 {
            let phi = if z == 0.0 { 1.0 } else { -det::expm1(-z) / z };
            4.0 * x * x * phi
        } else {
            (x / self.y) * -det::expm1(-z)
        };
        let value = rate * factor * det::exp(-t * t) * INV_SQRT_PI;
        if !nonnegative(value) {
            return Err(FissionError::Arithmetic.into());
        }
        Ok(value)
    }

    #[allow(clippy::too_many_arguments)]
    fn integrate_piece(
        &self,
        mut left: f64,
        right: f64,
        law: Law,
        work: &mut Work,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
        integral: &mut Accumulator,
        error: &mut Accumulator,
    ) -> Result<(), EndfError> {
        if matches!(law, Law::ConstantRate(0.0)) {
            return Ok(());
        }
        if let Law::Segment { segment, region } = law {
            let source = &self.owner.source;
            if source.cross_sections_m2[segment] == 0.0
                && (source.interpolation_regions[region].interpolation
                    == EndfInterpolation::Histogram
                    || source.cross_sections_m2[segment + 1] == 0.0)
            {
                return Ok(());
            }
        }
        // At most one thermal width per initial panel: a narrow Gaussian can
        // never hide between all the nodes on a very wide source segment.
        while left < right {
            work.visit(checkpoint)?;
            let end = (left + 1.0).min(right);
            let (value, estimate) = self.panel(left, end, law, 0, work, checkpoint)?;
            integral.add(value)?;
            error.add(estimate)?;
            left = end;
        }
        Ok(())
    }

    fn panel(
        &self,
        left: f64,
        right: f64,
        law: Law,
        depth: u8,
        work: &mut Work,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(f64, f64), EndfError> {
        checkpoint()?;
        let coarse = self.gauss(left, right, law, &GAUSS8, work, checkpoint)?;
        let fine = self.gauss(left, right, law, &GAUSS16, work, checkpoint)?;
        let error = (8.0 * (fine - coarse).abs()).max(64.0 * f64::EPSILON * fine);
        let allowed = 0.5
            * (self.absolute_tolerance * ((right - left) / self.width)
                + self.relative_tolerance * fine);
        if error <= allowed {
            return Ok((fine, error));
        }
        let middle = left + 0.5 * (right - left);
        if depth == MAX_DEPTH || middle == left || middle == right {
            return Err(budget("Doppler quadrature refinement"));
        }
        let first = self.panel(left, middle, law, depth + 1, work, checkpoint)?;
        let second = self.panel(middle, right, law, depth + 1, work, checkpoint)?;
        Ok((first.0 + second.0, first.1 + second.1))
    }

    fn gauss(
        &self,
        left: f64,
        right: f64,
        law: Law,
        nodes: &[(f64, f64)],
        work: &mut Work,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<f64, EndfError> {
        let half = 0.5 * (right - left);
        let middle = left + half;
        let mut total = Accumulator::new(0.0);
        for &(node, weight) in nodes {
            for t in [middle - half * node, middle + half * node] {
                work.evaluate(checkpoint)?;
                total.add(weight * self.value(t, law)?)?;
            }
        }
        Ok(half * total.value())
    }
}

struct Work {
    limits: DopplerLimits,
    visits: usize,
    evaluations: usize,
}

impl Work {
    fn new(limits: DopplerLimits) -> Self {
        Self {
            limits,
            visits: 0,
            evaluations: 0,
        }
    }

    fn visit(
        &mut self,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(), EndfError> {
        checkpoint()?;
        if self.visits >= self.limits.max_segment_visits {
            return Err(budget("Doppler segment visits"));
        }
        self.visits += 1;
        Ok(())
    }

    fn evaluate(
        &mut self,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(), EndfError> {
        if self.evaluations % 32 == 0 {
            checkpoint()?;
        }
        if self.evaluations >= self.limits.max_kernel_evaluations {
            return Err(budget("Doppler kernel evaluations"));
        }
        self.evaluations += 1;
        Ok(())
    }
}

fn invalid(reason: &'static str) -> EndfError {
    FissionError::InvalidInput(reason).into()
}

// Positive abscissae and weights on [-1, 1], used in symmetric pairs.
const GAUSS8: [(f64, f64); 4] = [
    (0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
];
const GAUSS16: [(f64, f64); 8] = [
    (0.095_012_509_837_637_44, 0.189_450_610_455_068_5),
    (0.281_603_550_779_258_9, 0.182_603_415_044_923_6),
    (0.458_016_777_657_227_4, 0.169_156_519_395_002_54),
    (0.617_876_244_402_643_8, 0.149_595_988_816_576_74),
    (0.755_404_408_355_003, 0.124_628_971_255_533_88),
    (0.865_631_202_387_831_8, 0.095_158_511_682_492_78),
    (0.944_575_023_073_232_6, 0.062_253_523_938_647_89),
    (0.989_400_934_991_649_9, 0.027_152_459_411_754_096),
];

#[cfg(test)]
mod tests;
