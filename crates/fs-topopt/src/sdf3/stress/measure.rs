//! Explicit averaging versus a conservative finite-sample strength functional.
use super::{StressError3, StressOptions3, finite, power};
use crate::SimpParams;

/// Stress functional bound to a fixed geometry, load family and sampling rule.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum StressMeasure3 {
    /// Original normalized volume/load-weighted qp aggregate. A low-weight load
    /// or a small-volume hotspot can exceed it. Zero-weight cases are excluded.
    #[default]
    NormalizedAverage,
    /// Unweighted p-norm of BOTH relaxed and physical von Mises stresses at ALL
    /// retained bulk samples of ALL cases, including cases with zero weight:
    /// `A = (sum_lq ((rho^q vm_ref)^p + (k(rho) vm_ref)^p))^(1/p)`.
    ///
    /// In exact arithmetic, `M <= A <= (2*N)^(1/p) M`, where M is the larger
    /// sampled relaxed/physical maximum and N counts retained points over cases.
    /// The implementation also checks A against the actual numerical maxima.
    /// No quadrature-volume or load normalization can dilute a hotspot. The
    /// physical term includes ersatz void stresses rather than erasing them.
    ///
    /// This changes the optimization problem and depends on the sample count;
    /// it is NOT an interval certificate or a bound between quadrature points.
    /// A finite p-norm and a finite tolerance do not certify continuum safety.
    SampledPeakBound,
}
impl StressMeasure3 {
    /// Stable declaration spelling for consumers binding source and recovery.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::NormalizedAverage => "normalized-average",
            Self::SampledPeakBound => "sampled-peak-bound",
        }
    }
    pub(super) fn includes(self, weight: f64) -> bool {
        self == Self::SampledPeakBound || weight > 0.0
    }
    pub(super) fn sample_weight(self, case: f64, point: f64, volume: f64) -> f64 {
        match self {
            Self::NormalizedAverage => case * (point / volume),
            Self::SampledPeakBound => 1.0,
        }
    }
}

// vm is positively homogeneous, so the per-point pair can be evaluated as
// c(rho)*vm_ref, c=(rho^(q*p)+k^p)^(1/p). Differentiate the algebraic norm,
// NOT a frozen maximum scale, and include BOTH direct density contributions.
// The scaling avoids underflow of k^p in a prescribed void at high p.
pub(super) fn peak_scale(rho: f64, relaxed: f64, stiffness: f64,
    params: SimpParams, options: StressOptions3) -> Result<(f64, f64), StressError3> {
    let p = options.aggregation_power;
    let m = relaxed.max(stiffness);
    if !m.is_finite() || m <= 0.0 {
        return Err(StressError3::Invalid("invalid sampled peak material coefficient"));
    }
    let coefficient = finite(m * power(power(relaxed / m, p) + power(stiffness / m, p), 1.0 / p))?;
    let drelaxed = options.relaxation_power * power(rho, options.relaxation_power - 1.0);
    let dstiffness = (1.0 - params.e_min) * params.penal * power(rho, params.penal - 1.0);
    let slope = finite(power(relaxed / coefficient, p - 1.0) * drelaxed
        + power(stiffness / coefficient, p - 1.0) * dstiffness)?;
    Ok((coefficient, slope))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sampled_pair_has_the_full_density_derivative_and_survives_high_p_voids() {
        let params = SimpParams { penal: 3.0, e_min: 1e-6, ..Default::default() };
        for p in [2.0, 8.0, 32.0, 64.0] {
            let options = StressOptions3 { aggregation_power: p, ..Default::default() };
            let sample = |r: f64| peak_scale(r, r,
                params.e_min + (1.0 - params.e_min) * power(r, params.penal), params, options).unwrap();
            let void = sample(0.0);
            assert!((void.0 / params.e_min - 1.0).abs() < 1e-14);
            assert_eq!(void.1, 0.0);
            for r in [0.001, 0.05, 0.3, 0.7, 0.999] {
                let (c, d) = sample(r);
                let h = 1e-6;
                let fd = (sample(r + h).0 - sample(r - h).0) / (2.0 * h);
                assert!((fd - d).abs() < 1e-6 * d.abs().max(1e-8), "p={p}, rho={r}");
                assert!(c >= r);
                assert!(c >= params.e_min + (1.0 - params.e_min) * power(r, params.penal));
            }
        }
    }
    #[test]
    fn sample_weight_cannot_hide_small_volume_or_zero_weight_cases() {
        let peak = StressMeasure3::SampledPeakBound;
        assert!(peak.includes(0.0));
        assert_eq!(peak.sample_weight(0.0, 1e-300, 1e300), 1.0);
        let average = StressMeasure3::NormalizedAverage;
        assert!(!average.includes(0.0));
        assert_eq!(average.sample_weight(0.25, 0.5, 2.0), 0.0625);
    }
}
