//! Number distributions shared by fixed-index and material-backed particles.

use super::{MAX_MIE_POPULATION_BINS, MieError, poll};
use fs_exec::Cx;
use fs_math::det;

/// Discrete radii in metres and normalized particle-number fractions.
///
/// Histogram inputs must be counts per bin. A sampled density needs its bin
/// width or quadrature weight before admission; mass/volume fractions cannot
/// be supplied as number fractions. The declaration does not set concentration.
#[derive(Clone, Debug, PartialEq)]
pub struct ParticleSizeDistribution {
    bins: Vec<(f64, f64)>,
}

impl ParticleSizeDistribution {
    /// Admit 1--256 positive finite radii with finite nonnegative relative
    /// number weights. At least one weight must be positive. Retained bins
    /// preserve input order; duplicate radii are allowed and zeros are removed.
    pub fn try_new(bins: &[(f64, f64)]) -> Result<Self, MieError> {
        if bins.is_empty()
            || bins.len() > MAX_MIE_POPULATION_BINS
            || bins
                .iter()
                .any(|(radius, _)| !radius.is_finite() || *radius <= 0.0)
        {
            return Err(MieError::Parameters);
        }
        Ok(Self {
            bins: normalize_number_bins(bins)?,
        })
    }

    /// Radius in metres paired with normalized number fraction.
    #[must_use]
    pub fn bins(&self) -> &[(f64, f64)] {
        &self.bins
    }

    /// Discretize a lognormal NUMBER distribution on an explicit finite
    /// radius interval. `geometric_std` is `exp(stddev(ln r))`, at least one;
    /// the median is a number median, not a mass or intensity median.
    ///
    /// Cells are uniform in log radius. Their weights are integrated normal
    /// probabilities, including the cell widths, and their radii are geometric
    /// cell midpoints. Weights are conditioned on the admitted finite interval.
    /// A geometric standard deviation of one gives an exact point mass when
    /// the median lies in the interval. Cancellation is polled between cells.
    ///
    /// Increase `bins` to assess quadrature convergence; independently widen
    /// `radius_range_m` to assess omitted tails. Reported number/volume tails
    /// are distribution integrals, NOT bounds on omitted optical scattering.
    /// A small number tail may carry substantial mass or resonant scattering.
    /// Every generated radius still needs admission by the Mie solver at each
    /// wavelength; this constructor never clips unsupported optical sizes.
    ///
    /// Distribution convention: NIST/SEMATECH e-Handbook, §1.3.6.6.9:
    /// <https://www.itl.nist.gov/div898/handbook/eda/section3/eda3669.htm>.
    pub fn lognormal(
        cx: &Cx<'_>,
        median_radius_m: f64,
        geometric_std: f64,
        radius_range_m: [f64; 2],
        bins: usize,
    ) -> Result<(Self, LognormalRadiusReport), MieError> {
        poll(cx)?;
        let [minimum, maximum] = radius_range_m;
        if [median_radius_m, minimum, maximum]
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
            || minimum >= maximum
            || !geometric_std.is_finite()
            || geometric_std < 1.0
            || !(1..=MAX_MIE_POPULATION_BINS).contains(&bins)
        {
            return Err(MieError::Parameters);
        }
        if geometric_std == 1.0 {
            if !(minimum..=maximum).contains(&median_radius_m) {
                return Err(MieError::Parameters);
            }
            return Ok((
                Self::try_new(&[(median_radius_m, 1.0)])?,
                LognormalRadiusReport {
                    radius_range_m,
                    bins: 1,
                    retained_number_fraction: 1.0,
                    omitted_number_fraction: 0.0,
                    omitted_volume_fraction: 0.0,
                },
            ));
        }
        let sigma = det::ln(geometric_std);
        let center = det::ln(median_radius_m);
        let (lower, upper) = (det::ln(minimum), det::ln(maximum));
        let (z_lower, z_upper) = ((lower - center) / sigma, (upper - center) / sigma);
        let retained = normal_interval(z_lower, z_upper);
        if !z_lower.is_finite()
            || !z_upper.is_finite()
            || !retained.is_finite()
            || retained <= 0.0
            || retained > 1.0
        {
            return Err(MieError::Numerical);
        }
        let mut nodes = Vec::new();
        nodes
            .try_reserve_exact(bins)
            .map_err(|_| MieError::Budget)?;
        for index in 0..bins {
            poll(cx)?;
            let lo = (upper - lower).mul_add(index as f64 / bins as f64, lower);
            let hi = if index + 1 == bins {
                upper
            } else {
                (upper - lower).mul_add((index + 1) as f64 / bins as f64, lower)
            };
            let radius = det::exp(0.5 * (lo + hi));
            let weight = normal_interval((lo - center) / sigma, (hi - center) / sigma);
            if !weight.is_finite() || weight <= 0.0 || !(minimum..=maximum).contains(&radius) {
                // Do not silently drop a nonempty positive-probability cell.
                return Err(MieError::Numerical);
            }
            nodes.push((radius, weight));
        }
        // Multiplication by r^3 shifts a lognormal's standardized normal
        // coordinate by 3*sigma. Constants cancel in the retained/full ratio.
        let report = LognormalRadiusReport {
            radius_range_m,
            bins,
            retained_number_fraction: retained,
            omitted_number_fraction: normal_tails(z_lower, z_upper),
            omitted_volume_fraction: normal_tails(z_lower - 3.0 * sigma, z_upper - 3.0 * sigma),
        };
        poll(cx)?;
        Ok((Self::try_new(&nodes)?, report))
    }
}

/// Finite-distribution coverage, separate from angular and wavelength errors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LognormalRadiusReport {
    /// Explicit lower/upper radii, metres.
    pub radius_range_m: [f64; 2],
    /// Number of retained quadrature radii (one for zero width).
    pub bins: usize,
    /// Probability inside the declared radius interval before conditioning.
    pub retained_number_fraction: f64,
    /// Probability excluded by the two radius tails of the full lognormal.
    pub omitted_number_fraction: f64,
    /// Fraction of the full lognormal's particle volume in omitted tails.
    /// Also the mass fraction when all radii share one material density.
    /// This does not bound omitted scattering or absorption.
    pub omitted_volume_fraction: f64,
}

fn normal_interval(lower: f64, upper: f64) -> f64 {
    let scale = core::f64::consts::FRAC_1_SQRT_2;
    if lower >= 0.0 {
        0.5 * (det::erfc(lower * scale) - det::erfc(upper * scale))
    } else if upper <= 0.0 {
        0.5 * (det::erfc(-upper * scale) - det::erfc(-lower * scale))
    } else {
        0.5 * (det::erf(upper * scale) - det::erf(lower * scale))
    }
}

fn normal_tails(lower: f64, upper: f64) -> f64 {
    let scale = core::f64::consts::FRAC_1_SQRT_2;
    0.5 * (det::erfc(-lower * scale) + det::erfc(upper * scale))
}

pub(super) fn normalize_number_bins<T: Copy>(bins: &[(T, f64)]) -> Result<Vec<(T, f64)>, MieError> {
    if bins.is_empty()
        || bins.len() > MAX_MIE_POPULATION_BINS
        || bins
            .iter()
            .any(|(_, value)| !value.is_finite() || *value < 0.0)
    {
        return Err(MieError::Parameters);
    }
    let maximum = bins.iter().map(|(_, value)| *value).fold(0.0_f64, f64::max);
    if maximum == 0.0 {
        return Err(MieError::Parameters);
    }
    let mut scaled = Vec::new();
    scaled
        .try_reserve_exact(bins.len())
        .map_err(|_| MieError::Budget)?;
    let mut total = Compensated::default();
    for &(item, value) in bins {
        if value == 0.0 {
            continue;
        }
        let weight = value / maximum;
        if weight == 0.0 {
            return Err(MieError::Numerical);
        }
        total.add(weight);
        scaled.push((item, weight));
    }
    let total = total.total();
    for (_, weight) in &mut scaled {
        *weight /= total;
        if *weight == 0.0 {
            return Err(MieError::Numerical);
        }
    }
    Ok(scaled)
}

#[derive(Default)]
pub(super) struct Compensated {
    sum: f64,
    correction: f64,
}

impl Compensated {
    pub(super) fn add(&mut self, value: f64) {
        let sum = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - sum) + value
        } else {
            (value - sum) + self.sum
        };
        self.sum = sum;
    }

    pub(super) fn total(self) -> f64 {
        self.sum + self.correction
    }
}
