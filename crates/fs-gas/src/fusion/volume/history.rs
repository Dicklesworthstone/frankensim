//! Positive, prescribed geometry for homogeneous fusion parcels.

use std::sync::Arc;

use fs_math::det;
use fs_qty::{Time, Volume};

use super::super::{DtFusionError, poll, positive};
use super::log_volume_ratio;

/// One authored physical time [s] and parcel volume [m³].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FusionVolumeKnot {
    pub time: Time,
    pub volume: Volume,
}

/// Immutable, bounded volume history with piecewise linear **log volume**.
///
/// On each interval, `V(t) = V0 * exp(s * ln(V1/V0))`, where `s` is elapsed
/// interval fraction. This prescribes a constant fractional compression rate
/// within that interval; equal consecutive volumes prescribe an exact hold.
/// Knots need not be uniformly spaced. The caller owns this geometry: no
/// piston dynamics, confinement law or linear-volume interpolation is inferred.
///
/// The original knots are retained without resampling. Construction allocates;
/// lookup uses bounded binary search and allocates nothing. Queries outside the
/// authored time range refuse, including for a one-knot, zero-duration history.
#[derive(Clone, Debug, PartialEq)]
pub struct FusionVolumeHistory {
    knots: Arc<[FusionVolumeKnot]>,
    log_ratios: Arc<[f64]>,
}

impl FusionVolumeHistory {
    /// Maximum number of authored knots; includes both interval endpoints.
    pub const MAX_KNOTS: usize = 256;

    /// Admit 1–256 positive finite volumes at strictly increasing finite times.
    /// The first time must equal zero. Cancellation is polled for every knot.
    pub fn new(
        knots: Vec<FusionVolumeKnot>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DtFusionError> {
        poll(&mut cancelled)?;
        if knots.is_empty() || knots.len() > Self::MAX_KNOTS {
            return Err(DtFusionError::Invalid(
                "volume history requires 1..=256 knots",
            ));
        }
        let mut log_ratios = Vec::with_capacity(knots.len() - 1);
        for (index, knot) in knots.iter().enumerate() {
            poll(&mut cancelled)?;
            if !positive(knot.volume.value()) {
                return Err(DtFusionError::Invalid(
                    "volume history requires positive finite volumes",
                ));
            }
            let time = knot.time.value();
            if !time.is_finite()
                || (index == 0 && time != 0.0)
                || (index > 0 && time <= knots[index - 1].time.value())
            {
                return Err(DtFusionError::Invalid(
                    "volume history times must start at zero and increase strictly",
                ));
            }
            if index > 0 {
                log_ratios.push(log_volume_ratio(knots[index - 1].volume, knot.volume)?);
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            knots: knots.into(),
            log_ratios: log_ratios.into(),
        })
    }

    /// Original authored knots, including every hold and turnaround.
    #[must_use]
    pub fn knots(&self) -> &[FusionVolumeKnot] {
        &self.knots
    }

    #[must_use]
    pub fn duration(&self) -> Time {
        self.knots[self.knots.len() - 1].time
    }

    #[must_use]
    pub fn initial_volume(&self) -> Volume {
        self.knots[0].volume
    }

    #[must_use]
    pub fn final_volume(&self) -> Volume {
        self.knots[self.knots.len() - 1].volume
    }

    /// Evaluate the prescribed log-linear volume without extrapolation.
    /// Authored knots and constant segments return their original volume bits.
    pub fn volume_at(&self, time: Time) -> Result<Volume, DtFusionError> {
        let after = self.upper_bound(time)?;
        let index = after - 1;
        let left = self.knots[index];
        if after == self.knots.len() || time == left.time {
            return Ok(left.volume);
        }
        let right = self.knots[after];
        if left.volume == right.volume {
            return Ok(left.volume);
        }
        let elapsed = time.value() - left.time.value();
        let remaining = right.time.value() - time.value();
        let duration = right.time.value() - left.time.value();
        // Anchor at the nearer endpoint. Computing the remaining fraction
        // independently avoids losing tiny changes by subtracting it from one.
        let (anchor, exponent) = if elapsed <= remaining {
            (
                left.volume.value(),
                (elapsed / duration) * self.log_ratios[index],
            )
        } else {
            (
                right.volume.value(),
                -(remaining / duration) * self.log_ratios[index],
            )
        };
        let volume = if exponent.abs() <= 0.5 {
            // expm1 + fused addition retains small relative changes even when
            // exp(exponent) rounds to one before multiplication by the volume.
            det::expm1(exponent).mul_add(anchor, anchor)
        } else {
            let factor = det::exp(exponent);
            let product = anchor * factor;
            if factor.is_normal() && positive(product) {
                product
            } else {
                // A finite volume can have an overflowing or subnormal ratio
                // to its anchor. Combine logarithms before exponentiating.
                det::exp(det::ln(anchor) + exponent)
            }
        };
        if !positive(volume) {
            return Err(DtFusionError::Arithmetic(
                "log-linear volume interpolation is not positive and finite",
            ));
        }
        Ok(Volume::new(volume))
    }

    /// First authored knot strictly after `time`; none at the final endpoint.
    /// Solvers can stop a step here before the fractional compression rate changes.
    pub fn next_knot_after(&self, time: Time) -> Result<Option<Time>, DtFusionError> {
        Ok(self
            .knots
            .get(self.upper_bound(time)?)
            .map(|knot| knot.time))
    }

    fn upper_bound(&self, time: Time) -> Result<usize, DtFusionError> {
        if !time.value().is_finite() {
            return Err(DtFusionError::Invalid(
                "volume history query time must be finite",
            ));
        }
        if time.value() < 0.0 || time > self.duration() {
            return Err(DtFusionError::OutsideDomain);
        }
        Ok(self.knots.partition_point(|knot| knot.time <= time))
    }
}

#[cfg(test)]
mod tests;
