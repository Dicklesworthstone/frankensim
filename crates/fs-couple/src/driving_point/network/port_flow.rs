//! Outward volume flow at each further port of a network, per unit wave incident on the active port.
//!
//! The matched source of [`NetworkDrivingPoint::reflection_regional`] (open-circuit pressure 2 behind the real
//! reference `zc`) is an incident wave `p+ = 1` at the active port. The same solve's inflow at another port, negated,
//! is that port's outward volume flow, `T_k = U_k / p+`: finite at the input-impedance poles, as the reflection is.
//! For a one-edge graph this is [`crate::driving_point::outlet`]'s `U_mouth / p+` of the same duct and load. No load,
//! loss or geometry law is added here: the network owner solves every edge with its own gas, and with its wall when
//! the media bind one. The causal FIR comes from the mouth-flow ladder itself, so a branched tract's nostrils and lips
//! each get the line a single duct's mouth gets. The radiated pressure of these flows belongs to the receiver.

use super::*;
use crate::driving_point::outlet::{MouthFlowRealization, transfer_dc_limit, transfer_ladder};
use crate::driving_point::validate_ladder;
use fs_duct::network::{NetworkMedia, network_drive_sources_regional};

#[cfg(test)]
mod tests;

impl NetworkDrivingPoint {
    /// `U_k / p+` at one angular frequency, `e^{-iwt}`: the outward volume flow at `port` per unit pressure wave
    /// incident on the active port, against the real reference `zc`.
    ///
    /// # Errors
    /// Declared coupled ports, the active port itself, a port the graph lacks, or the network owner's refusal.
    pub fn port_flow_per_incident_regional(
        &self,
        media: &NetworkMedia,
        omega: f64,
        loss: LossModel,
        zc: f64,
        port: &str,
    ) -> Result<C64, DuctError> {
        if !self.coupled_labels.is_empty() {
            return Err(DuctError::BadParameter {
                what: "network port flow with coupled ports needs their coupled termination",
            });
        }
        if port == self.active_port {
            return Err(DuctError::BadParameter {
                what: "a port flow is read at a port other than the driven one",
            });
        }
        let sources = self.reflection_source(zc)?;
        let response = network_drive_sources_regional(
            &self.network, media, omega, loss, &sources, &self.terminations, None,
        )?;
        let k = response
            .fields
            .ports
            .iter()
            .position(|label| label == port)
            .ok_or(DuctError::BadParameter {
                what: "a port flow names a port the graph lacks",
            })?;
        let flow = response.fields.port_inflows[k].scale(-1.0);
        if !flow.re.is_finite() || !flow.im.is_finite() {
            return Err(DuctError::Singular);
        }
        Ok(flow)
    }

    /// Causal impulse response of [`Self::port_flow_per_incident_regional`] on the mouth-flow ladder (the same
    /// arithmetic as a single duct's mouth): pushing `p+` through the FIR gives the port's outward flow. The one-way
    /// delay that seeds the ladder is half the regional round trip (to the farthest node).
    ///
    /// # Errors
    /// Non-positive rate or reference, a bad ladder, the port refusals above, or ladder non-convergence.
    #[allow(clippy::too_many_arguments)]
    pub fn realize_port_flow_ir_regional(
        &self,
        media: &NetworkMedia,
        loss: LossModel,
        port: &str,
        sample_rate_hz: u32,
        zc: f64,
        ladder: &LineLadder,
    ) -> Result<(Vec<f64>, MouthFlowRealization), DrivingPointError> {
        if sample_rate_hz == 0 || !(zc.is_finite() && zc > 0.0) {
            return Err(DrivingPointError::Invalid {
                what: "positive sample rate and finite positive port impedance required",
            });
        }
        let max_taps = validate_ladder(ladder)?;
        let dt = 1.0 / f64::from(sample_rate_hz);
        let transit = 0.5 * self.regional_round_trip_seconds(media)? / dt;
        let response = |omega: f64| {
            self.port_flow_per_incident_regional(media, omega, loss, zc, port)
                .map_err(DrivingPointError::Duct)
        };
        let dc = transfer_dc_limit(&response)?;
        transfer_ladder(&response, dc, transit, dt, max_taps, ladder)
    }
}
