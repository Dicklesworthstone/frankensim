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
    gas: GasState,
    ledger: AirwayLedger,
}

/// Relative slack of the airway's passivity certificate ([`PlayerTract::dg_commit`]): the work the source-free
/// airway has taken from the valve may be negative by this fraction of the wave energy through its port. The
/// realized reflectance is refused above `|R| = 1 + 4e-16` on its grid (`UpstreamAirway`), which allows a cumulative
/// deficit of about `4e-16` of that energy; the mouth-pressure load is admitted on the analytic `Re Z_in - R_dc Re H`
/// ([`admitted_margin`]), whose realization error is that of the ladder's tail tolerance. `1e-6` sits above both and
/// far below any active load, which takes an O(1) fraction.
const AIRWAY_PASSIVITY_TOLERANCE: f64 = 1e-6;

/// Cumulative work the source-free airway has taken from the valve, and the wave energy through its port
/// ([`PlayerTract::dg_commit`]). Moved only by an accepted DG sample; the voice route never touches it.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct AirwayLedger {
    absorbed_j: f64,
    wave_energy_j: f64,
}

/// One DG valve sample's airway port, between [`PlayerTract::dg_port`] and [`PlayerTract::dg_commit`].
#[derive(Debug, Clone, Copy)]
pub(crate) struct AirwayPort {
    /// The valve-side open-circuit pressure: the airway's, the regulation's DC offset and the series source [Pa].
    pub(crate) open_circuit_pa: f64,
    /// The airway's same-sample series impedance against the flow into the valve [Pa s/m^3].
    pub(crate) impedance: f64,
    source_pa: f64,
    series_pa: f64,
    /// The open-circuit pressure with the far-end source at rest, regulation included, series source excluded [Pa].
    free_open_circuit_pa: f64,
}

/// What one accepted DG valve sample did to the airway: the separate term the DG frame reports for a tract that is
/// not a DG duct (`DgApertureFrame::airway`). The airway has no stored energy of its own in the DG account: its
/// FIRs are a lumped, band-limited port law, so the work it delivers and takes is the whole of its ledger.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AirwayLedgerFrame {
    /// The far-end source: the card's blowing pressure [Pa].
    pub source_pa: f64,
    /// The series source at the valve (the opt-in constriction noise) [Pa].
    pub series_pa: f64,
    /// The open-circuit pressure at the valve this sample, series source included [Pa].
    pub open_circuit_pa: f64,
    /// The series impedance against the flow [Pa s/m^3].
    pub impedance_pa_s_m3: f64,
    /// The flow from the airway into the valve [m^3/s].
    pub flow_m3_s: f64,
    /// The pressure at the valve face, series source included [Pa].
    pub face_pressure_pa: f64,
    /// Work the airway, its far-end source and its regulation deliver through the valve face, series source
    /// excluded: `(face - series) U dt` [J]. `DgApertureFrame::balance_residual_j` credits it.
    pub delivered_work_j: f64,
    /// Work the source-free airway took from the valve this sample [J]: `-(p_free) U dt`, `p_free` the face pressure
    /// with the far-end source at rest.
    pub absorbed_work_j: f64,
    /// The sum of `absorbed_work_j` over every accepted sample [J]; it is refused below `-1e-6` of `wave_energy_j`.
    pub cumulative_absorbed_work_j: f64,
    /// The sum over samples of `(p_free^2 / Z_s + Z_s U^2) dt / 2`, the wave energy through the port [J].
    pub wave_energy_j: f64,
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
        self.commit_port(step.open_circuit_pa, step.impedance, flow_m3_s)
    }

    /// The airway takes `flow_m3_s`, the one-pole advances, and the face pressure is returned [Pa]; the arithmetic of
    /// [`Self::commit`], which the DG route shares. A refusal leaves the airway as it was.
    fn commit_port(
        &mut self,
        open_circuit_pa: f64,
        impedance: f64,
        flow_m3_s: f64,
    ) -> Result<f64, AcousticRealizeError> {
        let mouth_pa = open_circuit_pa - impedance * flow_m3_s;
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

impl PlayerTract {
    /// Realize `chain` filled with `gas` on a `sample_rate_hz` clock, the card's pressure acting where `pressure`
    /// says, at rest: the tract of [`ReedBoreVoice::set_player_tract`] and of a DG valve's upstream owner. A mouth
    /// pressure's load is read on every positive bin of the realized reflectance grid and refused if it is active
    /// anywhere ([`admitted_margin`]).
    ///
    /// # Errors
    /// A zero clock; a `tau_s` that is not positive and finite; the airway's own refusals; a negative passivity margin.
    pub fn realize(
        chain: &AirwayChain,
        gas: &GasState,
        pressure: TractPressure,
        sample_rate_hz: u32,
        ladder: &LineLadder,
    ) -> Result<Self, AcousticRealizeError> {
        if sample_rate_hz == 0 {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "a player tract needs a positive clock",
            });
        }
        let airway =
            UpstreamAirway::realize(chain, gas, sample_rate_hz, ladder).map_err(airway_error)?;
        let dc_resistance = airway.receipt().dc_resistance;
        let dt = 1.0 / f64::from(sample_rate_hz);
        let (pole, margin) = match pressure {
            TractPressure::Lung => (0.0, None),
            TractPressure::Mouth { tau_s } => {
                if !(tau_s > 0.0 && tau_s.is_finite()) {
                    return Err(AcousticRealizeError::InvalidDescription {
                        what: "a mouth-pressure tract needs a positive finite time constant",
                    });
                }
                let pole = 1.0 - fs_math::det::exp(-dt / tau_s);
                if !(pole > 0.0 && pole <= 1.0) {
                    return Err(AcousticRealizeError::InvalidDescription {
                        what: "the tract's one-pole coefficient must lie in (0, 1]",
                    });
                }
                let n = airway.receipt().reflectance.n_fft;
                let margin = admitted_margin(chain, gas, sample_rate_hz, n, |omega| {
                    one_pole_response(pole, omega, dt)
                })?;
                (pole, Some(margin))
            }
        };
        Ok(Self {
            airway,
            pressure,
            dc_resistance,
            pole,
            flow_lp: 0.0,
            margin,
            mouth_pa: 0.0,
            flow_m3_s: 0.0,
            gas: *gas,
            ledger: AirwayLedger::default(),
        })
    }

    /// The gas the airway was realized in: the inlet gas of the valve it sits upstream of.
    #[must_use]
    pub const fn gas(&self) -> &GasState {
        &self.gas
    }

    /// The airway's same-sample series impedance against the flow into the valve [Pa s/m^3]: the impedance the valve
    /// is declared against, in the DG route a constant of the realization.
    #[must_use]
    pub const fn series_impedance_pa_s_m3(&self) -> f64 {
        self.airway.receipt().series_impedance
    }

    /// The sum of the work the source-free airway took from the valve over every accepted DG sample [J]
    /// ([`AirwayLedgerFrame::cumulative_absorbed_work_j`]); zero in the voice route.
    #[must_use]
    pub const fn cumulative_absorbed_work_j(&self) -> f64 {
        self.ledger.absorbed_j
    }

    /// Whether no sample has been taken: the one-pole, the last flows and the DG ledger are all at rest. The airway's
    /// FIR histories move only with a committed flow, so a zero flow and an empty ledger mean they are empty too.
    #[must_use]
    pub fn is_cold(&self) -> bool {
        self.airway.commit_count() == 0
            && self.flow_lp == 0.0
            && self.flow_m3_s == 0.0
            && self.mouth_pa == 0.0
            && self.ledger == AirwayLedger::default()
    }

    /// This DG valve sample's port with `source_pa` at the airway's far end and a series source `series_pa` at the
    /// valve. Nothing moves until [`Self::dg_commit`]; the open-circuit pressure is [`Self::fold`]'s, term for term.
    pub(crate) fn dg_port(
        &mut self,
        source_pa: f64,
        series_pa: f64,
    ) -> Result<AirwayPort, AcousticRealizeError> {
        if !series_pa.is_finite() {
            return Err(AcousticRealizeError::Reed {
                what: "player tract series source must be finite",
            });
        }
        // The shadow of the flow-only response starts with the first DG sample, on an airway at rest.
        self.airway.track_source_free_response().map_err(airway_error)?;
        let port = self.airway.port(source_pa).map_err(airway_error)?;
        let offset = self.dc_offset_pa();
        let open_circuit_pa = port.open_circuit_pa + offset + series_pa;
        let free_open_circuit_pa = self.airway.source_free_open_circuit_pa().ok_or(AcousticRealizeError::Reed {
            what: "player tract has no source-free response",
        })? + offset;
        if !(open_circuit_pa.is_finite() && free_open_circuit_pa.is_finite() && port.impedance.is_finite()) {
            return Err(AcousticRealizeError::Reed {
                what: "player tract port left the finite set",
            });
        }
        Ok(AirwayPort {
            open_circuit_pa,
            impedance: port.impedance,
            source_pa,
            series_pa,
            free_open_circuit_pa,
        })
    }

    /// Close this DG valve sample at the flow `flow_m3_s` from the airway into the valve: the work the source-free
    /// airway takes is added to its ledger, which is the passivity certificate (a cumulative deficit beyond
    /// `AIRWAY_PASSIVITY_TOLERANCE` of the wave energy refuses by name), then the airway takes the flow. A refusal
    /// leaves the airway, the one-pole and the ledger as they were.
    pub(crate) fn dg_commit(
        &mut self,
        port: AirwayPort,
        flow_m3_s: f64,
    ) -> Result<AirwayLedgerFrame, AcousticRealizeError> {
        let dt = 1.0 / f64::from(self.airway.sample_rate_hz());
        let zs = self.airway.port_impedance();
        let face_pressure_pa = port.open_circuit_pa - port.impedance * flow_m3_s;
        let free_face_pa = port.free_open_circuit_pa - port.impedance * flow_m3_s;
        let absorbed_work_j = -free_face_pa * flow_m3_s * dt;
        let wave_j = 0.5 * (free_face_pa * free_face_pa / zs + zs * flow_m3_s * flow_m3_s) * dt;
        let cumulative_absorbed_work_j = self.ledger.absorbed_j + absorbed_work_j;
        let wave_energy_j = self.ledger.wave_energy_j + wave_j;
        let delivered_work_j = (face_pressure_pa - port.series_pa) * flow_m3_s * dt;
        if ![face_pressure_pa, free_face_pa, absorbed_work_j, wave_j, cumulative_absorbed_work_j, wave_energy_j,
            delivered_work_j].iter().all(|x| x.is_finite())
        {
            return Err(AcousticRealizeError::Reed {
                what: "player tract ledger left the finite set",
            });
        }
        if cumulative_absorbed_work_j < -AIRWAY_PASSIVITY_TOLERANCE * wave_energy_j {
            return Err(AcousticRealizeError::Nonlinear(format!(
                "player tract refused: the source-free airway has taken {cumulative_absorbed_work_j:.6e} J from the \
                 valve against {wave_energy_j:.6e} J of wave energy through its port (limit -{AIRWAY_PASSIVITY_TOLERANCE:e} of it), \
                 so the tract is acting as a source"
            )));
        }
        self.commit_port(port.open_circuit_pa, port.impedance, flow_m3_s)?;
        self.ledger = AirwayLedger { absorbed_j: cumulative_absorbed_work_j, wave_energy_j };
        Ok(AirwayLedgerFrame {
            source_pa: port.source_pa,
            series_pa: port.series_pa,
            open_circuit_pa: port.open_circuit_pa,
            impedance_pa_s_m3: port.impedance,
            flow_m3_s,
            face_pressure_pa,
            delivered_work_j,
            absorbed_work_j,
            cumulative_absorbed_work_j,
            wave_energy_j,
        })
    }
}

#[cfg(test)]
impl PlayerTract {
    /// This tract with a mouth-pressure regulation of time constant `tau_s` and DC resistance `dc_resistance`
    /// [Pa s/m^3] that no admission read: the planted active load of the passivity certificate's test.
    pub(crate) fn with_unadmitted_regulation(mut self, tau_s: f64, dc_resistance: f64) -> Self {
        let dt = 1.0 / f64::from(self.airway.sample_rate_hz());
        self.pressure = TractPressure::Mouth { tau_s };
        self.pole = 1.0 - fs_math::det::exp(-dt / tau_s);
        self.dc_resistance = dc_resistance;
        self.margin = None;
        self
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
        self.player_tract = Some(PlayerTract::realize(
            chain,
            &self.inlet_gas,
            pressure,
            sample_rate_hz,
            ladder,
        )?);
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
