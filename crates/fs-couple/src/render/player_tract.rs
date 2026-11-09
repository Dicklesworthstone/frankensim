//! The player's tract upstream of the reed (opt-in, [`ReedBoreVoice::set_player_tract`]): a measured airway
//! ([`crate::upstream_airway`]) between the card's blowing pressure and the reed junction, the reed-side use that module
//! names.
//!
//! **Where the card's pressure acts.** [`TractPressure::Lung`]: the card's pressure is the source at the airway's far
//! end, and the reed sees the whole airway, `p = p_oc - Z_up U`. [`TractPressure::Mouth`]: the card's pressure was
//! measured in the mouth, so it already holds the tract's steady drop. The airway's DC resistance `R_dc`
//! ([`AirwayChain::dc_resistance`]) then comes out at DC only: the mouth pressure is the airway's port pressure plus
//! `R_dc s_n`, with `s_n = s_{n-1} + a (U_{n-1} - s_{n-1})` and `a = 1 - exp(-dt / tau)`, a strictly causal one-pole
//! low-pass of the flow into the tract. The reed sees the load `Z_up(w) - R_dc H(w)` with `H(0) = 1`, so the steady
//! mouth pressure is the card's and the tract keeps its dynamics above about `1 / (2 pi tau)`. `tau` is a numerical
//! regulation constant standing for the player holding a mean mouth pressure, not a measurement; the caller names it.
//!
//! A flat `Z_up - R_dc` is not offered. A measured tract's input resistance falls below its DC resistance between its
//! resonances (participant A of Kaburagi et al. 2021, over about 0.6 to 6 kHz), so the flat subtraction would make the
//! tract an active load there. [`admitted_margin`] reads `Re Z_in(w) - R_dc Re H(w)`, with the chain's analytic
//! `Z_in`, on every positive bin of the realized airway's reflectance grid and refuses a negative one; the voice
//! takes no tract it refuses.
//!
//! **The junction.** The airway's same-sample impedance `Z` is folded into the bore's affine boundary before the reed
//! solve. In its outgoing wave `a`, the bore obeys `p = (1 + d) a + h` and `Z_c U = (1 - d) a - h`. With `eps = Z / Z_c`
//! and `alpha = 1 + eps (1 - d) / 2`, the variable `q = alpha a - eps h / 2` obeys the same two laws with
//! `d' = (d + eps (1 - d) / 2) / alpha` and `h' = h / alpha`, its pressure being the bore's plus `Z U`. The reed solves
//! against the mouth pressure `p_oc + Z u_body`, `u_body` being the plates' flow, which the tract does not carry, so the
//! drop across the reed is `p_oc - Z (U - u_body) - p_bore`. Every reed law (massless, viscous slit, staple recovery,
//! jet contraction, massive) then takes the tract with its own arithmetic unchanged; the massive path's load coordinate
//! sees the impedance `Z_c (1 + d) / (1 - d) + Z` and the same incoming wave. The flow into the tract is the bore's
//! minus the plates': the reed's jet plus the volume its face sweeps.
//!
//! Limits, not hidden: the airway's own (cylinders, plane waves, a symmetric tree, a stationary gas, no junction or
//! separation loss where the area steps); one held tract shape over the whole note (a measured posture, not an
//! articulation); the voice's inlet gas in every section; an opt-in constriction noise stays a series source at the
//! reed. Without a tract the voice is the original, bit for bit.

use super::*;
use crate::driving_point::{DrivingPointError, LineLadder};
use crate::upstream_airway::{AirwayChain, UpstreamAirway};
use fs_math::c64::C64;
use fs_vfit::scattering::ScatteringResponse;

/// Where a card's blowing pressure acts when a player's tract sits upstream of the reed (module doc).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TractPressure {
    /// A lung pressure: the source at the airway's far end, the whole airway at the reed.
    Lung,
    /// A mouth pressure: `R_dc` out at DC only, through the one-pole of time constant `tau_s`.
    Mouth {
        /// The one-pole's time constant [s].
        tau_s: f64,
    },
}

/// A player's tract realized on the voice's clock, with its one-pole state (module doc).
#[derive(Debug, Clone)]
pub struct PlayerTract {
    airway: UpstreamAirway,
    pressure: TractPressure,
    dc_resistance: f64,
    pole: f64,
    flow_lp: f64,
    margin: Option<(f64, f64)>,
    mouth_pa: f64,
    flow_m3_s: f64,
}

/// One sample's tract junction, between [`PlayerTract::fold`] and [`PlayerTract::commit`].
#[derive(Debug, Clone, Copy)]
pub(super) struct TractStep {
    open_circuit_pa: f64,
    impedance: f64,
    half_eps_h: f64,
    alpha: f64,
}

impl TractStep {
    /// The bore's own outgoing wave from the folded solve's outgoing wave `q` (module doc).
    pub(super) fn unfold(self, folded_outgoing: f64) -> f64 {
        (folded_outgoing + self.half_eps_h) / self.alpha
    }
}

fn airway_error(e: DrivingPointError) -> AcousticRealizeError {
    AcousticRealizeError::Nonlinear(format!("player tract airway: {e}"))
}

/// `H(w)` of the one-pole `s_n = s_{n-1} + a (U_{n-1} - s_{n-1})` on a clock of period `dt` (module doc).
#[must_use]
pub fn one_pole_response(pole: f64, omega: f64, dt: f64) -> C64 {
    let (s, c) = fs_math::det::sin_cos(omega * dt);
    let delay = C64::new(c, s);
    delay.scale(pole) * (C64::ONE - delay.scale(1.0 - pole)).recip()
}

/// The smallest `Re Z_in(w) - R_dc Re H(w)` [Pa s/m^3] of `chain` in `gas` over the bins `w_k = 2 pi k rate / n`,
/// `k = 1..=n/2`, with its frequency [Hz]: the load a mouth-pressure tract puts at the reed is passive where this is
/// not negative (module doc). `response(w)` is `H` at `w`. `Z_in` is the chain's ANALYTIC input impedance at each
/// bin ([`AirwayChain::input_impedance`]), not the realized band-limited FIR port, so the refusal is the model's.
///
/// # Errors
/// A zero clock or `n` below 2; the chain's own refusal at a bin; a non-finite margin; a negative margin, refused by
/// name with its value and frequency.
pub fn admitted_margin(
    chain: &AirwayChain,
    gas: &GasState,
    sample_rate_hz: u32,
    n: usize,
    response: impl Fn(f64) -> C64,
) -> Result<(f64, f64), AcousticRealizeError> {
    if sample_rate_hz == 0 || n < 2 {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "a tract passivity read needs a positive clock and at least two bins",
        });
    }
    let r_dc = chain.dc_resistance(gas);
    let mut worst = (f64::INFINITY, 0.0);
    for k in 1..=n / 2 {
        let hz = k as f64 * f64::from(sample_rate_hz) / n as f64;
        let omega = core::f64::consts::TAU * hz;
        let z = chain.input_impedance(gas, omega).map_err(airway_error)?;
        let margin = z.re - r_dc * response(omega).re;
        if !margin.is_finite() {
            return Err(AcousticRealizeError::Reed {
                what: "player tract passivity margin left the finite set",
            });
        }
        if margin < worst.0 {
            worst = (margin, hz);
        }
    }
    if worst.0 < 0.0 {
        return Err(AcousticRealizeError::Nonlinear(format!(
            "player tract refused: Re Z_in - R_dc Re H is {:.6e} Pa s/m^3 at {:.3} Hz, so the mouth-pressure load would be \
             active there",
            worst.0, worst.1
        )));
    }
    Ok(worst)
}

impl PlayerTract {
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

    /// The mouth pressure at the reed in the last rendered sample [Pa].
    #[must_use]
    pub const fn mouth_pressure_pa(&self) -> f64 {
        self.mouth_pa
    }

    /// The flow into the tract in the last rendered sample [m^3/s].
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

    /// This sample's fold (module doc): the airway's port with `source_pa` at its far end, a series source `series_pa`
    /// at the reed, the bore's `boundary` on characteristic impedance `zc`, and the plates' flow. Returns the boundary and
    /// the mouth pressure the reed solves against.
    pub(super) fn fold(
        &mut self,
        source_pa: f64,
        series_pa: f64,
        boundary: ScatteringResponse,
        zc: f64,
        body_flow_m3_s: f64,
    ) -> Result<(ScatteringResponse, f64, TractStep), AcousticRealizeError> {
        let port = self.airway.port(source_pa).map_err(airway_error)?;
        let open_circuit_pa = port.open_circuit_pa + self.dc_offset_pa() + series_pa;
        let impedance = port.impedance;
        let (d, h) = (boundary.direct, boundary.history);
        let eps = impedance / zc;
        let half = 0.5 * eps * (1.0 - d);
        let alpha = 1.0 + half;
        let folded = ScatteringResponse {
            direct: (d + half) / alpha,
            history: h / alpha,
        };
        let mouth_pa = open_circuit_pa + impedance * body_flow_m3_s;
        if !(alpha.is_finite()
            && alpha > 0.0
            && folded.direct.is_finite()
            && folded.history.is_finite()
            && mouth_pa.is_finite())
        {
            return Err(AcousticRealizeError::Reed {
                what: "player tract junction left the finite set",
            });
        }
        Ok((
            folded,
            mouth_pa,
            TractStep {
                open_circuit_pa,
                impedance,
                half_eps_h: 0.5 * eps * h,
                alpha,
            },
        ))
    }

    /// Close the sample at the flow `flow_m3_s` into the tract: the airway takes it, the one-pole advances, and the
    /// mouth pressure at the reed is returned [Pa].
    pub(super) fn commit(
        &mut self,
        step: TractStep,
        flow_m3_s: f64,
    ) -> Result<f64, AcousticRealizeError> {
        let mouth_pa = step.open_circuit_pa - step.impedance * flow_m3_s;
        if !mouth_pa.is_finite() {
            return Err(AcousticRealizeError::Reed {
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

impl ReedBoreVoice {
    /// Put a player's tract upstream of the reed (module doc): `chain`, valve side first, realized in the voice's inlet
    /// gas on its `sample_rate_hz` clock with `ladder`, the card's blowing pressure acting where `pressure` says. A mouth
    /// pressure's load is read on every positive bin of the realized reflectance grid and refused if it is active
    /// anywhere ([`admitted_margin`]). Admission is allowed only before the first rendered sample.
    ///
    /// # Errors
    /// An advanced voice or one that already has a tract; a clock other than the voice's; a `tau_s` that is not positive
    /// and finite; the airway's own refusals; a negative passivity margin.
    pub fn set_player_tract(
        &mut self,
        chain: &AirwayChain,
        pressure: TractPressure,
        sample_rate_hz: u32,
        ladder: &LineLadder,
    ) -> Result<(), AcousticRealizeError> {
        if self.sample_index != 0 || self.player_tract.is_some() {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "a player tract needs a voice that has not advanced and has no tract yet",
            });
        }
        if sample_rate_hz == 0 || self.dt != 1.0 / f64::from(sample_rate_hz) {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "a player tract is realized on the voice's own clock",
            });
        }
        let gas = self.inlet_gas;
        let airway =
            UpstreamAirway::realize(chain, &gas, sample_rate_hz, ladder).map_err(airway_error)?;
        let dc_resistance = airway.receipt().dc_resistance;
        let (pole, margin) = match pressure {
            TractPressure::Lung => (0.0, None),
            TractPressure::Mouth { tau_s } => {
                if !(tau_s > 0.0 && tau_s.is_finite()) {
                    return Err(AcousticRealizeError::InvalidDescription {
                        what: "a mouth-pressure tract needs a positive finite time constant",
                    });
                }
                let pole = 1.0 - fs_math::det::exp(-self.dt / tau_s);
                if !(pole > 0.0 && pole <= 1.0) {
                    return Err(AcousticRealizeError::InvalidDescription {
                        what: "the tract's one-pole coefficient must lie in (0, 1]",
                    });
                }
                let dt = self.dt;
                let n = airway.receipt().reflectance.n_fft;
                let margin = admitted_margin(chain, &gas, sample_rate_hz, n, |omega| {
                    one_pole_response(pole, omega, dt)
                })?;
                (pole, Some(margin))
            }
        };
        self.player_tract = Some(PlayerTract {
            airway,
            pressure,
            dc_resistance,
            pole,
            flow_lp: 0.0,
            margin,
            mouth_pa: 0.0,
            flow_m3_s: 0.0,
        });
        Ok(())
    }

    /// The player's tract upstream of the reed, if one was set.
    #[must_use]
    pub const fn player_tract(&self) -> Option<&PlayerTract> {
        self.player_tract.as_ref()
    }
}

#[cfg(test)]
mod tests;
