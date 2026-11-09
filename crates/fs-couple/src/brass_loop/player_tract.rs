//! The player's tract upstream of the lip (opt-in, [`super::BrassVoice::with_player_tract`]): a measured airway
//! ([`crate::upstream_airway`]) between the card's blowing pressure and the lip's aperture, the lip-side use that module
//! names. The reed side is `crate::render::player_tract`, which owns the pure pieces ([`TractPressure`],
//! [`one_pole_response`], [`admitted_margin`]) this module re-exports; the card rows are its rows, word for word.
//!
//! **Where the card's pressure acts.** [`TractPressure::Lung`]: the card's blowing pressure is the source at the airway's far
//! end and the aperture sees the whole airway, `p = p_oc - Z U`. [`TractPressure::Mouth`]: the card's pressure was measured
//! in the mouth, so it already holds the tract's steady drop (Fletcher and Tarnopolsky's blowing pressures are mouth
//! pressures). The airway's DC resistance `R_dc` ([`AirwayChain::dc_resistance`]) then comes out at DC only: the mouth
//! pressure is the airway's port pressure plus `R_dc s_n`, with `s_n = s_{n-1} + a (U_{n-1} - s_{n-1})` and
//! `a = 1 - exp(-dt / tau)`, a strictly causal one-pole low-pass of the flow into the tract. The aperture sees the load
//! `Z_up(w) - R_dc H(w)` with `H(0) = 1`, so the steady mouth pressure is the card's and the tract keeps its dynamics
//! above about `1 / (2 pi tau)`. `tau` is a numerical regulation constant standing for the player holding a mean mouth
//! pressure, not a measurement; the caller names it. A flat `Z_up - R_dc` is not offered: [`admitted_margin`] reads
//! `Re Z_in(w) - R_dc Re H(w)` on every positive bin of the realized airway's reflectance grid and refuses a negative
//! one by name, so the voice takes no tract that would be an active load.
//!
//! **The junction.** The lip's flow solve is the Bernoulli law against the cup junction's `p = u / a + k` (`a` an
//! admittance); the tract is a Thevenin source `p_mouth = p_oc - Z u` in series with the aperture, so the same quadratic
//! holds with `head = p_oc - k` and `1 / a` replaced by `1 / a + Z` (nothing new is solved). The tract's flow is the
//! aperture flow, the flow the junction takes: this lip has no swept-volume flow of its own. The lip island is explicit, as
//! without a tract; its force takes the mouth pressure `p_oc - Z u_prev`, this sample's open-circuit pressure and
//! impedance against the previous sample's flow, which is the card's pressure exactly when `Z` is zero.
//!
//! Limits, not hidden: the airway's own (cylinders, plane waves, a symmetric tree, a stationary gas, no junction or
//! separation loss where the area steps); one held tract shape over the whole note (a measured posture, not an
//! articulation); the voice's bore gas in every section. Not combined with the jet noise (its bound reads the pressure
//! drop of the ideal-source drive) or the nonlinear bore (a second step); the constriction noise reads only the held pressure
//! and rides on the tract's open-circuit pressure as a series source. Without a tract the voice is
//! the original, bit for bit.

use crate::driving_point::{DrivingPointError, LineLadder};
use crate::mm_line::MmLineError;
use crate::upstream_airway::{AirwayChain, UpstreamAirway};
use fs_material::gas::GasState;

pub use crate::render::player_tract::{TractPressure, admitted_margin, one_pole_response};

/// One sample's tract port, between [`LipTract::port`] and [`LipTract::commit`].
#[derive(Debug, Clone, Copy)]
pub(super) struct LipTractPort {
    /// Mouth pressure with no flow this sample, the one-pole's DC offset included [Pa].
    pub(super) open_circuit_pa: f64,
    /// Series impedance against the flow [Pa s/m^3].
    pub(super) impedance: f64,
}

/// A player's tract realized on the voice's clock, with its one-pole state (module doc).
#[derive(Debug, Clone)]
pub struct LipTract {
    airway: UpstreamAirway,
    pressure: TractPressure,
    dc_resistance: f64,
    pole: f64,
    flow_lp: f64,
    margin: Option<(f64, f64)>,
    mouth_pa: f64,
    flow_m3_s: f64,
}

fn airway_error(_: DrivingPointError) -> MmLineError {
    MmLineError::Invalid {
        what: "player tract airway refused",
    }
}

impl LipTract {
    /// Realize `chain` (aperture side first) filled with `gas` on a `sample_rate_hz` clock whose period is `dt` with
    /// `ladder`, the card's blowing pressure acting where `pressure` says (module doc). A mouth pressure's load is read on
    /// every positive bin of the realized reflectance grid and refused if it is active anywhere ([`admitted_margin`]).
    ///
    /// # Errors
    /// A clock other than `1 / dt`; a `tau_s` that is not positive and finite; the airway's own refusals; a negative
    /// passivity margin.
    pub fn realize(
        chain: &AirwayChain,
        pressure: TractPressure,
        gas: &GasState,
        sample_rate_hz: u32,
        dt: f64,
        ladder: &LineLadder,
    ) -> Result<LipTract, String> {
        if sample_rate_hz == 0 || dt != 1.0 / f64::from(sample_rate_hz) {
            return Err("a player tract is realized on the voice's own clock".to_string());
        }
        let airway = UpstreamAirway::realize(chain, gas, sample_rate_hz, ladder)
            .map_err(|e| format!("player tract airway refused: {e:?}"))?;
        let dc_resistance = airway.receipt().dc_resistance;
        let (pole, margin) = match pressure {
            TractPressure::Lung => (0.0, None),
            TractPressure::Mouth { tau_s } => {
                if !(tau_s > 0.0 && tau_s.is_finite()) {
                    return Err(
                        "a mouth-pressure tract needs a positive finite time constant".to_string(),
                    );
                }
                let pole = 1.0 - fs_math::det::exp(-dt / tau_s);
                if !(pole > 0.0 && pole <= 1.0) {
                    return Err("the tract's one-pole coefficient must lie in (0, 1]".to_string());
                }
                let n = airway.receipt().reflectance.n_fft;
                let margin = admitted_margin(chain, gas, sample_rate_hz, n, |omega| {
                    one_pole_response(pole, omega, dt)
                })
                .map_err(|e| e.to_string())?;
                (pole, Some(margin))
            }
        };
        Ok(LipTract {
            airway,
            pressure,
            dc_resistance,
            pole,
            flow_lp: 0.0,
            margin,
            mouth_pa: 0.0,
            flow_m3_s: 0.0,
        })
    }

    /// The realized airway.
    #[must_use]
    pub const fn airway(&self) -> &UpstreamAirway {
        &self.airway
    }

    /// Where the card's pressure acts.
    #[must_use]
    pub const fn pressure(&self) -> TractPressure {
        self.pressure
    }

    /// The airway's DC resistance `R_dc` [Pa s/m^3].
    #[must_use]
    pub const fn dc_resistance(&self) -> f64 {
        self.dc_resistance
    }

    /// The admitted passivity margin of a mouth-pressure tract, [`admitted_margin`]'s `(value, frequency)`; `None` for a
    /// lung pressure, whose load is the airway's own.
    #[must_use]
    pub const fn passivity_margin(&self) -> Option<(f64, f64)> {
        self.margin
    }

    /// The mouth pressure at the aperture in the last rendered sample [Pa].
    #[must_use]
    pub const fn mouth_pressure_pa(&self) -> f64 {
        self.mouth_pa
    }

    /// The flow through the tract (the aperture flow) in the last rendered sample [m^3/s].
    #[must_use]
    pub const fn flow_m3_s(&self) -> f64 {
        self.flow_m3_s
    }

    fn dc_offset_pa(&self) -> f64 {
        match self.pressure {
            TractPressure::Lung => 0.0,
            TractPressure::Mouth { .. } => self.dc_resistance * self.flow_lp,
        }
    }

    /// This sample's port (module doc): the airway's Thevenin source with `source_pa` at its far end and the one-pole's DC
    /// offset. Nothing moves until [`Self::commit`].
    pub(super) fn port(&mut self, source_pa: f64) -> Result<LipTractPort, MmLineError> {
        let port = self.airway.port(source_pa).map_err(airway_error)?;
        let open_circuit_pa = port.open_circuit_pa + self.dc_offset_pa();
        if !(open_circuit_pa.is_finite() && port.impedance.is_finite()) {
            return Err(MmLineError::Invalid {
                what: "player tract port left the finite set",
            });
        }
        Ok(LipTractPort {
            open_circuit_pa,
            impedance: port.impedance,
        })
    }

    /// The mouth pressure the lip's force reads this sample: the port's open-circuit pressure against the previous sample's
    /// flow (module doc).
    pub(super) fn force_pressure_pa(&self, port: LipTractPort) -> f64 {
        port.open_circuit_pa - port.impedance * self.flow_m3_s
    }

    /// Close the sample at the flow `flow_m3_s` the aperture solved for: the airway takes it, the one-pole advances, and the
    /// mouth pressure is returned [Pa].
    pub(super) fn commit(
        &mut self,
        port: LipTractPort,
        flow_m3_s: f64,
    ) -> Result<f64, MmLineError> {
        let mouth_pa = port.open_circuit_pa - port.impedance * flow_m3_s;
        if !mouth_pa.is_finite() {
            return Err(MmLineError::Invalid {
                what: "player tract mouth pressure left the finite set",
            });
        }
        self.airway.commit(flow_m3_s).map_err(airway_error)?;
        self.flow_lp += self.pole * (flow_m3_s - self.flow_lp);
        self.mouth_pa = mouth_pa;
        self.flow_m3_s = flow_m3_s;
        Ok(mouth_pa)
    }
}

#[cfg(test)]
mod tests;
