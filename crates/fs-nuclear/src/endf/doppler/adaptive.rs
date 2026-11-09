//! Adaptive linear-in-energy tabulation of the original Doppler operator.
//!
//! Accepted intervals pass three interior interpolation probes. A source-knot
//! search prevents narrow broadened features from hiding between those probes:
//! cells whose retained thermal support meets a retained source feature must be at most
//! the admitted width in sqrt(E/theta). Wide energy ratios are split
//! geometrically to resolve the low-energy inverse-speed tail efficiently.
//!
//! This is a deterministic interpolation ESTIMATE, not an interval enclosure
//! or a guarantee for arbitrary spectra. The original quadrature estimates and
//! omitted-kernel bounds remain separate, and their uncertainty is included
//! when comparing each interior probe with its endpoint interpolant.

use super::{
    DopplerAccuracy, DopplerSample, DopplerTabulation, EndfError, EndfInterpolation,
    EndfInterpolationRegion, FissionError, PreparedDopplerBroadening, Work, budget, det, invalid,
    nonnegative, positive,
};

const MAX_GRID_DEPTH: u8 = 64;

/// Controls interpolation between output knots, independently of point
/// integration accuracy. A finite nonzero absolute tolerance is useful where
/// the reaction vanishes or a heated resonance wing becomes negligible.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DopplerGridAccuracy {
    pub relative_tolerance: f64,
    pub absolute_tolerance_m2: f64,
    /// Maximum sqrt(E/theta) span near a retained source feature, in (0, 1].
    /// Endpoints and duplicate-energy jumps are retained. Only interior knots
    /// joining exactly equal constant intervals are omitted from this guard.
    pub maximum_thermal_span: f64,
}

impl Default for DopplerGridAccuracy {
    fn default() -> Self {
        Self {
            relative_tolerance: 1.0e-3,
            absolute_tolerance_m2: 0.0,
            maximum_thermal_span: 0.5,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveDopplerTabulation {
    /// Contains the original pointwise type consumed by the collapse APIs.
    /// Its numerical-error maxima include both retained knots and probes.
    pub tabulation: DopplerTabulation,
    /// Maximum accepted interior discrepancy, plus the integration uncertainty
    /// of that probe and the interpolated endpoint uncertainty. This sampled
    /// estimate is not a bound on unprobed energies or nuclear-data error.
    pub maximum_estimated_interpolation_error_m2: f64,
    /// Counts actual point integrations; shared endpoints and midpoints are
    /// reused during refinement instead of integrated again.
    pub point_evaluations: usize,
}

impl PreparedDopplerBroadening {
    /// Construct a native curve over an explicit closed energy interval.
    ///
    /// Source-knot searches, interval refinements and all point integrations
    /// share one `DopplerLimits` work allowance. Source points are searched,
    /// not copied into the output: a densely sampled input can produce a
    /// smaller heated curve. Output memory is bounded by `max_output_points`.
    ///
    /// Heating requires 0 < lower < upper. At unchanged temperature the
    /// requested interval must be the complete original domain; the original
    /// grid, interpolation laws and jumps are then retained exactly. As with
    /// `tabulate`, a new data ID labels the result without mutating its source.
    pub fn tabulate_adaptive(
        &self,
        data_id: String,
        energy_interval_j: [f64; 2],
        point_accuracy: DopplerAccuracy,
        grid_accuracy: DopplerGridAccuracy,
        limits: super::DopplerLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<AdaptiveDopplerTabulation, EndfError> {
        checkpoint()?;
        if data_id.trim().is_empty() {
            return Err(invalid("identified adaptive Doppler output required"));
        }
        if self.source.energies_j.len() > limits.max_source_points {
            return Err(budget("adaptive Doppler source points"));
        }
        let mut tabulation = DopplerTabulation {
            curve: self.source.clone(),
            source_data_id: self.source.data_id.clone(),
            source_temperature_k: self.source.metadata.temperature_k,
            tails: self.tails,
            maximum_estimated_quadrature_error_m2: 0.0,
            maximum_omitted_kernel_tail_bound_m2: 0.0,
        };
        tabulation.curve.data_id = data_id.into();
        if self.sqrt_theta_j == 0.0 {
            let energies = &self.source.energies_j;
            if energy_interval_j != [energies[0], energies[energies.len() - 1]] {
                return Err(invalid(
                    "identity adaptive Doppler tabulation retains the complete original domain",
                ));
            }
            if energies.len() > limits.max_output_points {
                return Err(budget("adaptive Doppler output points"));
            }
            checkpoint()?;
            return Ok(AdaptiveDopplerTabulation {
                tabulation,
                maximum_estimated_interpolation_error_m2: 0.0,
                point_evaluations: 0,
            });
        }
        let [lower, upper] = energy_interval_j;
        if !positive(lower)
            || !positive(upper)
            || lower >= upper
            || !nonnegative(grid_accuracy.relative_tolerance)
            || grid_accuracy.relative_tolerance > 1.0
            || !nonnegative(grid_accuracy.absolute_tolerance_m2)
            || grid_accuracy.relative_tolerance == 0.0 && grid_accuracy.absolute_tolerance_m2 == 0.0
            || !positive(grid_accuracy.maximum_thermal_span)
            || grid_accuracy.maximum_thermal_span > 1.0
        {
            return Err(invalid(
                "positive increasing adaptive domain, nonzero grid tolerance and thermal span in (0,1] required",
            ));
        }
        if limits.max_output_points < 2 {
            return Err(budget("adaptive Doppler output points"));
        }
        let mut grid = Grid {
            owner: self,
            point_accuracy,
            accuracy: grid_accuracy,
            work: Work::new(limits),
            energies: Vec::new(),
            values: Vec::new(),
            maximum_quadrature_error: 0.0,
            maximum_kernel_tail: 0.0,
            maximum_interpolation_error: 0.0,
            point_evaluations: 0,
        };
        let left = grid.at(lower, &mut checkpoint)?;
        let right = grid.at(upper, &mut checkpoint)?;
        let middle = grid.at(interior(lower, upper)?, &mut checkpoint)?;
        grid.push(left)?;
        grid.refine(left, middle, right, 0, &mut checkpoint)?;

        tabulation.curve.interpolation_regions = vec![EndfInterpolationRegion {
            end_point: grid.energies.len(),
            interpolation: EndfInterpolation::Linear,
        }]
        .into();
        tabulation.curve.energies_j = grid.energies.into();
        tabulation.curve.cross_sections_m2 = grid.values.into();
        tabulation.curve.metadata.temperature_k = self.target_temperature_k;
        tabulation.curve.metadata.derivation_flag = 1;
        tabulation.maximum_estimated_quadrature_error_m2 = grid.maximum_quadrature_error;
        tabulation.maximum_omitted_kernel_tail_bound_m2 = grid.maximum_kernel_tail;
        checkpoint()?;
        Ok(AdaptiveDopplerTabulation {
            tabulation,
            maximum_estimated_interpolation_error_m2: grid.maximum_interpolation_error,
            point_evaluations: grid.point_evaluations,
        })
    }
}

#[derive(Clone, Copy)]
struct Knot {
    energy: f64,
    sample: DopplerSample,
}

struct Grid<'a> {
    owner: &'a PreparedDopplerBroadening,
    point_accuracy: DopplerAccuracy,
    accuracy: DopplerGridAccuracy,
    work: Work,
    energies: Vec<f64>,
    values: Vec<f64>,
    maximum_quadrature_error: f64,
    maximum_kernel_tail: f64,
    maximum_interpolation_error: f64,
    point_evaluations: usize,
}

impl Grid<'_> {
    fn at(
        &mut self,
        energy: f64,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Knot, EndfError> {
        self.work.visit(checkpoint)?;
        let sample =
            self.owner
                .evaluate(energy, self.point_accuracy, &mut self.work, checkpoint)?;
        self.point_evaluations += 1;
        self.maximum_quadrature_error = self
            .maximum_quadrature_error
            .max(sample.estimated_quadrature_error_m2);
        self.maximum_kernel_tail = self
            .maximum_kernel_tail
            .max(sample.omitted_kernel_tail_bound_m2);
        Ok(Knot { energy, sample })
    }

    fn push(&mut self, knot: Knot) -> Result<(), EndfError> {
        let limit = self.work.limits.max_output_points;
        if self.energies.len() >= limit {
            return Err(budget("adaptive Doppler output points"));
        }
        // Exact bounded growth avoids an unbounded doubling reserve. Both
        // vectors are private until the complete output has been admitted.
        for values in [&mut self.energies, &mut self.values] {
            if values.len() == values.capacity() {
                let additional = (limit - values.len()).min(256);
                values
                    .try_reserve_exact(additional)
                    .map_err(|_| budget("adaptive Doppler output memory"))?;
            }
        }
        self.energies.push(knot.energy);
        self.values.push(knot.sample.cross_section_m2);
        Ok(())
    }

    fn refine(
        &mut self,
        left: Knot,
        middle: Knot,
        right: Knot,
        depth: u8,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(), EndfError> {
        self.work.visit(checkpoint)?;
        let quarter = self.at(interior(left.energy, middle.energy)?, checkpoint)?;
        let three_quarters = self.at(interior(middle.energy, right.energy)?, checkpoint)?;
        let mut accepted = true;
        let mut maximum_error = 0.0_f64;
        for probe in [quarter, middle, three_quarters] {
            let upper_weight = (probe.energy - left.energy) / (right.energy - left.energy);
            let lower_weight = (right.energy - probe.energy) / (right.energy - left.energy);
            let prediction = super::super::convex(
                left.sample.cross_section_m2,
                right.sample.cross_section_m2,
                lower_weight,
                upper_weight,
            )?;
            let endpoint_error = super::super::convex(
                point_error(left.sample)?,
                point_error(right.sample)?,
                lower_weight,
                upper_weight,
            )?;
            let error = (probe.sample.cross_section_m2 - prediction).abs()
                + point_error(probe.sample)?
                + endpoint_error;
            let allowed = self.accuracy.absolute_tolerance_m2
                + self.accuracy.relative_tolerance * probe.sample.cross_section_m2;
            if !nonnegative(error) || !nonnegative(allowed) {
                return Err(FissionError::Arithmetic.into());
            }
            if point_error(probe.sample)? > allowed {
                return Err(budget(
                    "Doppler point integration error exceeds adaptive grid tolerance",
                ));
            }
            accepted &= error <= allowed;
            maximum_error = maximum_error.max(error);
        }
        if self.owner.rate_bound_m2_sqrt_j != 0.0 {
            accepted &= right.energy / left.energy <= 4.0;
            if accepted {
                accepted = !self.thermal_guard(left.energy, right.energy, checkpoint)?;
            }
        }
        if accepted {
            self.maximum_interpolation_error = self.maximum_interpolation_error.max(maximum_error);
            return self.push(right);
        }
        if depth == MAX_GRID_DEPTH {
            return Err(budget("adaptive Doppler refinement depth"));
        }
        self.refine(left, quarter, middle, depth + 1, checkpoint)?;
        self.refine(middle, three_quarters, right, depth + 1, checkpoint)
    }

    fn thermal_guard(
        &mut self,
        lower: f64,
        upper: f64,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<bool, EndfError> {
        let low_root = det::sqrt(lower);
        let high_root = det::sqrt(upper);
        let width = ((upper - lower) / (high_root + low_root)) / self.owner.sqrt_theta_j;
        if width <= self.accuracy.maximum_thermal_span {
            return Ok(false);
        }
        // A retained feature affects the integral somewhere in this cell
        // precisely when its speed lies within L thermal widths of the cell.
        // The prepared cache omits only exactly flat redundant interior knots;
        // point quadrature still follows every original source interval.
        // Rationalized square-root differences retain narrow windows at high E.
        let coordinate = |energy: f64, reference: f64, root: f64| {
            ((energy - reference) / (det::sqrt(energy) + root)) / self.owner.sqrt_theta_j
        };
        let energies = &self.owner.thermal_features_j;
        let (mut first, mut end) = (0, energies.len());
        while first < end {
            self.work.visit(checkpoint)?;
            let mid = first + (end - first) / 2;
            if coordinate(energies[mid], lower, low_root) < -self.point_accuracy.kernel_half_widths
            {
                first = mid + 1;
            } else {
                end = mid;
            }
        }
        Ok(first < energies.len()
            && coordinate(energies[first], upper, high_root)
                <= self.point_accuracy.kernel_half_widths)
    }
}

fn point_error(sample: DopplerSample) -> Result<f64, EndfError> {
    let error = sample.estimated_quadrature_error_m2 + sample.omitted_kernel_tail_bound_m2;
    if nonnegative(error) {
        Ok(error)
    } else {
        Err(FissionError::Arithmetic.into())
    }
}

fn interior(lower: f64, upper: f64) -> Result<f64, EndfError> {
    let middle = if upper / lower > 4.0 {
        det::sqrt(lower) * det::sqrt(upper)
    } else {
        lower + 0.5 * (upper - lower)
    };
    if middle > lower && middle < upper {
        Ok(middle)
    } else {
        Err(budget("representable adaptive Doppler refinement"))
    }
}
