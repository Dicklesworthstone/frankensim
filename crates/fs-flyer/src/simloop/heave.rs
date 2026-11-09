//! Physical COM specific force for the reduced pitch-plane lifecycle.

use super::{Phase, SimLoop, TerminalEvent, det, motion};
use crate::Refusal;
use crate::perception::PerceivedCues;

/// Last acquired heave input, before perceptual delay/filter/gain/remnant.
///
/// Body FRD z is down; units are m/s². The reference is supported rest at the
/// initial trim pitch, `-g cos(theta_trim)`. The acquisition tick names the
/// interval START (`tick / 120` seconds), not its published endpoint.
/// Airborne samples use that accepted start's force and incidence; rail
/// samples use the actual constrained interval impulse, including support,
/// tow, friction and a no-reverse stop. Neither is an impact acceleration.
/// Reduced lateral motion does not bank the translational/sensor frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeaveCueSample {
    /// Physical acquisition interval's start tick.
    pub tick: u64,
    /// Body-down specific force minus the declared reference [m/s²].
    pub relative_specific_force_mps2: f64,
    /// Supported-trim body-down specific force [m/s²].
    pub reference_specific_force_mps2: f64,
}

pub(super) fn rail_specific_force(theta: f64, forward_acceleration: f64) -> f64 {
    det::sin(theta) * forward_acceleration - motion::GRAVITY_MPS2 * det::cos(theta)
}

impl SimLoop {
    pub(super) fn heave_reference(&self) -> f64 {
        -motion::GRAVITY_MPS2 * det::cos(self.tick0.trim.alpha_rad)
    }

    /// Read the most recently acquired physical heave cue without another
    /// force solve or state change. Fresh and envelope-failed ticks have no
    /// accepted sample. The existing perception ring is the sole history.
    #[must_use]
    pub fn heave_cue_sample(&self) -> Option<HeaveCueSample> {
        if self.phase == Phase::Ended(TerminalEvent::EnvelopeExceeded) {
            return None;
        }
        let (tick, raw) = self.perception.latest_heave_input(&self.perception_state)?;
        if tick.checked_add(1) != Some(self.tick) {
            return None;
        }
        Some(HeaveCueSample {
            tick,
            relative_specific_force_mps2: raw,
            reference_specific_force_mps2: self.heave_reference(),
        })
    }

    /// Last cues presented to the pilot, including the exactly six-tick heave
    /// delay and existing filter/remnant. Querying does not consume randomness.
    /// The frozen snapshot payload and registered pilot gains are unchanged.
    ///
    /// # Errors
    /// Perception remnant-stream refusals pass through.
    pub fn perceived_cues(&self) -> Result<Option<PerceivedCues>, Refusal> {
        self.perception.last_output(&self.perception_state)
    }
}

#[cfg(test)]
mod tests;
