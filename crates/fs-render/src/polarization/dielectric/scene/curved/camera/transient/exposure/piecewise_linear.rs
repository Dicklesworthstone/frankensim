//! Exact integration of the declared piecewise-linear source interpolant.
//!
//! Original knot times are seconds relative to the retained DD source epoch.
//! Queries subtract that epoch before any binary64 conversion. Integration
//! clips every crossed original segment; even a narrow interior pulse cannot
//! disappear between temporal quadrature nodes. Work is logarithmic for a
//! point query and linear in crossed segments for an interval, with no query
//! allocation. This is an authored interpolant, not a bandwidth certificate.

use super::{
    Dd, SurfaceEmissionHistory, SurfaceEmissionQuery, SurfaceExposureQuery, SurfaceHistoryError,
};
use fs_exec::Cx;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PiecewiseLinearHistoryError {
    InvalidEpoch,
    InvalidKnots,
    Cancelled,
}

/// Borrowed nonnegative intensity multipliers on a strictly bounded source
/// clock. Times are original local offsets; no resampling or extrapolation is
/// applied. The multiplier is achromatic and spatially uniform on its bound
/// source; the original emitter still owns its full spectral radiance law.
#[derive(Clone, Copy, Debug)]
pub struct PiecewiseLinearHistory<'a> {
    epoch_s: Dd,
    knots: &'a [(f64, f64)],
}

impl<'a> PiecewiseLinearHistory<'a> {
    pub fn try_new(
        cx: &Cx<'_>,
        epoch_s: Dd,
        knots: &'a [(f64, f64)],
    ) -> Result<Self, PiecewiseLinearHistoryError> {
        cx.checkpoint()
            .map_err(|_| PiecewiseLinearHistoryError::Cancelled)?;
        let epoch_s =
            super::admit_time(epoch_s).map_err(|_| PiecewiseLinearHistoryError::InvalidEpoch)?;
        if knots.len() < 2 {
            return Err(PiecewiseLinearHistoryError::InvalidKnots);
        }
        for (index, &(time, scale)) in knots.iter().enumerate() {
            cx.checkpoint()
                .map_err(|_| PiecewiseLinearHistoryError::Cancelled)?;
            if !time.is_finite()
                || !scale.is_finite()
                || scale < 0.0
                || (index > 0 && knots[index - 1].0 >= time)
            {
                return Err(PiecewiseLinearHistoryError::InvalidKnots);
            }
        }
        let extent = Dd::from_f64(knots[knots.len() - 1].0) - Dd::from_f64(knots[0].0);
        if !super::finite_dd(extent) || extent.hi <= 0.0 {
            return Err(PiecewiseLinearHistoryError::InvalidKnots);
        }
        Ok(Self { epoch_s, knots })
    }

    pub const fn epoch_s(&self) -> Dd {
        self.epoch_s
    }
    pub const fn knots(&self) -> &'a [(f64, f64)] {
        self.knots
    }

    pub fn value_at(&self, cx: &Cx<'_>, time_s: Dd) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        let local = self.local_time(time_s)?;
        checked_scale(self.interpolate(self.segment(local), local))
    }

    /// Mean over the entire interval, including every original interior knot.
    /// Normalized segment weights avoid a large energy intermediate; DD sums
    /// retain weak pulse areas next to brighter intervals without subtracting
    /// two nearly equal cumulative integrals.
    pub fn mean_between(
        &self,
        cx: &Cx<'_>,
        open_s: Dd,
        close_s: Dd,
    ) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        let open = self.local_time(open_s)?;
        let close = self.local_time(close_s)?;
        let width = close - open;
        if !open.lt(close) || !super::finite_dd(width) || width.hi <= 0.0 {
            return Err(SurfaceHistoryError::NumericalRange);
        }
        let mut index = self.segment(open);
        let mut left = open;
        let mut sum = Dd::ZERO;
        let mut positive = false;
        loop {
            cx.checkpoint()?;
            let endpoint = Dd::from_f64(self.knots[index + 1].0);
            let right = if endpoint.lt(close) { endpoint } else { close };
            let segment_width = right - left;
            let fraction = segment_width / width;
            let a = self.interpolate(index, left);
            let b = self.interpolate(index, right);
            // A difference preserves equal tiny endpoints and avoids adding
            // two large intensities before division by two.
            let mean = a + (b - a) * Dd::from_f64(0.5);
            positive |= Dd::ZERO.lt(a) || Dd::ZERO.lt(b);
            // A narrow bright pulse can have a subnormal duration fraction
            // while its contribution is ordinary. Prefer its area first. If
            // that product overflows, mean/width is ordinary and the opposite
            // order retains the representable final contribution.
            let contribution = if fraction.hi < f64::MIN_POSITIVE {
                let area = mean * segment_width;
                if super::finite_dd(area) {
                    area / width
                } else {
                    (mean / width) * segment_width
                }
            } else {
                mean * fraction
            };
            sum = sum + contribution;
            if !right.lt(close) {
                break;
            }
            left = right;
            index += 1;
        }
        if positive && sum == Dd::ZERO {
            return Err(SurfaceHistoryError::NumericalRange);
        }
        checked_scale(sum)
    }

    fn local_time(&self, time_s: Dd) -> Result<Dd, SurfaceHistoryError> {
        let time_s = super::admit_time(time_s).map_err(|_| SurfaceHistoryError::NumericalRange)?;
        let local = time_s - self.epoch_s;
        if !super::finite_dd(local) {
            return Err(SurfaceHistoryError::NumericalRange);
        }
        if local.lt(Dd::from_f64(self.knots[0].0))
            || Dd::from_f64(self.knots[self.knots.len() - 1].0).lt(local)
        {
            return Err(SurfaceHistoryError::OutsideTimeDomain);
        }
        Ok(local)
    }

    fn segment(&self, local: Dd) -> usize {
        self.knots
            .partition_point(|&(time, _)| !local.lt(Dd::from_f64(time)))
            .saturating_sub(1)
            .min(self.knots.len() - 2)
    }

    fn interpolate(&self, index: usize, local: Dd) -> Dd {
        let (a, va) = self.knots[index];
        let (b, vb) = self.knots[index + 1];
        if local == Dd::from_f64(a) || va == vb {
            return Dd::from_f64(va);
        }
        if local == Dd::from_f64(b) {
            return Dd::from_f64(vb);
        }
        let fraction = (local - Dd::from_f64(a)) / (Dd::from_f64(b) - Dd::from_f64(a));
        Dd::from_f64(va) * (Dd::ONE - fraction) + Dd::from_f64(vb) * fraction
    }
}

impl SurfaceEmissionHistory for PiecewiseLinearHistory<'_> {
    fn sample(&self, cx: &Cx<'_>, query: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        self.value_at(cx, query.emission_time_s)
    }

    fn average(
        &self,
        cx: &Cx<'_>,
        query: SurfaceExposureQuery,
    ) -> Result<f64, SurfaceHistoryError> {
        self.mean_between(cx, query.opening.emission_time_s, query.emission_close_s)
    }
}

fn checked_scale(value: Dd) -> Result<f64, SurfaceHistoryError> {
    let scale = value.to_f64();
    if !super::finite_dd(value) || !scale.is_finite() || value.lt(Dd::ZERO) {
        Err(SurfaceHistoryError::NumericalRange)
    } else {
        Ok(scale)
    }
}
