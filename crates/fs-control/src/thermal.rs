//! Authored thermal derating of the sampled shaft-command envelope.
//!
//! Each measured component has its own full-output, zero-output/trip and reset
//! temperatures. Between the first two, an explicitly linear derating law
//! scales the authored torque and motoring/regeneration power limits. The most
//! restrictive component owns the result; temperatures are never averaged.
//! Reaching any trip temperature latches zero torque demand until an explicit
//! reset is requested with EVERY monitored component below its reset limit.
//!
//! This policy consumes actual caller-owned temperatures; it neither predicts
//! cooling nor changes thermal energy. The returned limits feed the existing
//! speed controller before the original machine allocates current. Accept its
//! proposed latch state with the same accepted controller/plant sample.
//!
//! A zero torque request is not an electrical disconnect or a braking command.
//! Stored current, ongoing shaft motion and regeneration remain the native
//! plant's physics. This sampled policy does not certify between-sample peak
//! temperature, component survival, fault isolation or hardware protection.

use crate::ControlError;
use crate::speed::ShaftPowerLimits;

/// Fixed upper bound on work per sample; no sample allocates or sorts.
pub const MAX_THERMAL_DERATING_CHANNELS: usize = 64;

/// One explicitly authored component-temperature law, in kelvin.
///
/// `0 < reset <= full_output < trip`. Equality at `trip` trips; equality at
/// `reset` permits an explicit reset. These are control settings, not material
/// validity limits or automatically inferred damage thresholds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalDeratingBand {
    pub full_output_temperature_k: f64,
    pub trip_temperature_k: f64,
    pub reset_temperature_k: f64,
}

/// Accepted thermal-policy history. Restoring `tripped=true` keeps the latch
/// until a cooled, explicit reset; changing a command never clears it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ThermalDeratingState {
    pub tripped: bool,
}

/// Limits and candidate latch state for the same sampled plant interval.
#[derive(Clone, Copy, Debug, PartialEq)]
#[must_use = "apply limits and accept next_state with the corresponding plant sample"]
pub struct ThermalDeratingProposal {
    pub next_state: ThermalDeratingState,
    /// Feed these limits directly to `SpeedController::demand_controlled`.
    pub limits: ShaftPowerLimits,
    /// Applied common factor in [0,1], including any retained trip latch.
    pub output_fraction: f64,
    /// Most restrictive CURRENT measurement, with lowest index winning ties.
    /// None means all measurements currently permit full output, including
    /// when a previously tripped latch is still withholding the command.
    pub limiting_channel: Option<usize>,
    pub trip_started: bool,
    pub reset_accepted: bool,
    pub reset_blocked: bool,
}

/// Immutable component bands in the caller's declared sensor order.
/// Every sample must supply exactly one temperature for every original band.
#[derive(Clone, Copy, Debug)]
pub struct ThermalDeratingPolicy<'a> {
    bands: &'a [ThermalDeratingBand],
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ControlError> {
    if cancelled() {
        Err(ControlError::Cancelled)
    } else {
        Ok(())
    }
}

impl<'a> ThermalDeratingPolicy<'a> {
    pub fn new(bands: &'a [ThermalDeratingBand]) -> Result<Self, ControlError> {
        if bands.is_empty() || bands.len() > MAX_THERMAL_DERATING_CHANNELS {
            return Err(ControlError::InvalidParameter {
                what: "thermal derating requires 1..=64 component bands",
            });
        }
        if bands.iter().any(|b| {
            !b.full_output_temperature_k.is_finite()
                || !b.trip_temperature_k.is_finite()
                || !b.reset_temperature_k.is_finite()
                || b.reset_temperature_k <= 0.0
                || b.reset_temperature_k > b.full_output_temperature_k
                || b.full_output_temperature_k >= b.trip_temperature_k
        }) {
            return Err(ControlError::InvalidParameter {
                what: "thermal bands require finite 0 < reset <= full_output < trip in kelvin",
            });
        }
        Ok(Self { bands })
    }

    pub fn bands(&self) -> &'a [ThermalDeratingBand] {
        self.bands
    }

    /// Evaluate every component without mutating caller-owned history.
    ///
    /// Each active linear factor is `(trip-T)/(trip-full_output)`. The minimum
    /// factor scales all four ORIGINAL directional torque/power capabilities,
    /// preserving their asymmetry. A trip has priority over a simultaneous
    /// reset request; a retained latch clears only when all channels admit it.
    /// Polling occurs before work, once per channel and before publication.
    pub fn propose_controlled(
        &self,
        state: &ThermalDeratingState,
        temperature_k: &[f64],
        original_limits: ShaftPowerLimits,
        reset_requested: bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ThermalDeratingProposal, ControlError> {
        poll(&mut cancelled)?;
        if temperature_k.len() != self.bands.len() {
            return Err(ControlError::InvalidSample {
                what: "each original thermal channel needs one simultaneous temperature",
            });
        }
        if [
            original_limits.positive_torque_limit_nm,
            original_limits.negative_torque_limit_nm,
            original_limits.motoring_power_w,
            original_limits.regeneration_power_w,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(ControlError::InvalidSample {
                what: "original thermal torque and power limits must be finite nonnegative",
            });
        }
        let mut fraction = 1.0;
        let mut limiting_channel = None;
        let mut trip_now = false;
        let mut cooled = true;
        for (i, (&temperature, band)) in temperature_k.iter().zip(self.bands).enumerate() {
            poll(&mut cancelled)?;
            if !temperature.is_finite() || temperature <= 0.0 {
                return Err(ControlError::InvalidSample {
                    what: "thermal measurements must be finite positive kelvin",
                });
            }
            trip_now |= temperature >= band.trip_temperature_k;
            cooled &= temperature <= band.reset_temperature_k;
            let channel_fraction = if temperature <= band.full_output_temperature_k {
                1.0
            } else if temperature >= band.trip_temperature_k {
                0.0
            } else {
                (band.trip_temperature_k - temperature)
                    / (band.trip_temperature_k - band.full_output_temperature_k)
            };
            if channel_fraction < fraction {
                fraction = channel_fraction;
                limiting_channel = Some(i);
            }
        }
        let reset_accepted = state.tripped && reset_requested && cooled && !trip_now;
        let tripped = trip_now || (state.tripped && !reset_accepted);
        if tripped {
            fraction = 0.0;
        }
        let limits = ShaftPowerLimits {
            positive_torque_limit_nm: original_limits.positive_torque_limit_nm * fraction,
            negative_torque_limit_nm: original_limits.negative_torque_limit_nm * fraction,
            motoring_power_w: original_limits.motoring_power_w * fraction,
            regeneration_power_w: original_limits.regeneration_power_w * fraction,
        };
        let proposal = ThermalDeratingProposal {
            next_state: ThermalDeratingState { tripped },
            limits,
            output_fraction: fraction,
            limiting_channel,
            trip_started: !state.tripped && trip_now,
            reset_accepted,
            reset_blocked: reset_requested && tripped,
        };
        poll(&mut cancelled)?;
        Ok(proposal)
    }
}
