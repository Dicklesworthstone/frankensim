//! Opt-in momentum drive of the jet in the pipe: card row `momentum_drive KAPPA`, after
//! `jet_response`, with KAPPA required and only the two printed values accepted.
//!
//! The forced jet's cut flow `Q_j` enters the pipe as a VOLUME source: the series pressure
//! `p_s = source_coeff dQ_j/dt`, `source_coeff = rho (0.61 a_eq) / S_m` (Fletcher 1976, J. Acoust.
//! Soc. Am. 60, 481: a mouth mass accelerated by the entering flow). The same drive equation has a
//! second, in-phase term, the momentum flux of the jet the labium lets into the pipe. Thwaites and
//! Fletcher (1983, J. Acoust. Soc. Am. 74, 400, Eq. (1) p. 401):
//!
//! ```text
//! Q_p = [ rho (v0 + j w dl) / S_p ] / (Z_p + Z_m) * Q_j
//! ```
//!
//! `v0` is the jet speed at the flue, `dl` the mouth end correction and `S_p` the pipe's cross-section;
//! the two drives are a quarter cycle apart. This module adds the `v0` term to the series source:
//!
//! ```text
//! p_mom = kappa rho v0 (Q_j - Q_rest) / S_p
//! ```
//!
//! * `v0` is the bulk jet speed the loop already uses for `q_max = v0 h W` (the instantaneous
//!   speed, so it follows a pressure gesture); `S_p` is [`FlueGeometry::bore_area_m2`], the bore
//!   the jet enters; `rho` is the voice's gas density.
//! * **The steady value.** `Q_rest = q_max/2 (1 + tanh(-y0/b_cut))` is the cut flow of the QUIESCENT
//!   jet (labium displacement zero), at the instantaneous jet speed. It is exactly the cut flow of
//!   a voice at rest, so the term is exactly zero at rest and through any attack in which the jet
//!   is quiescent: no source enters a closed pipe from the steady jet. It is not a running mean.
//!   Once the loop oscillates, the mean of `Q_j` departs from `Q_rest` (the tanh cut is
//!   asymmetric about the lip), so the term then also delivers a slowly varying pressure offset;
//!   a closed pipe absorbs it as a static pressure with no steady flow, and nothing here removes it.
//! * **The same-sample solve.** The term has no `1/dt`: it adds `kappa rho v0 / S_p` to the slope
//!   of the series source with respect to `Q_j`. The guard `feedback < 1` of `solve_jet_mouth` now
//!   uses `source_coeff/dt + kappa rho v0 / S_p`; it is never loosened, and a voice without the
//!   row evaluates exactly the expressions it did before (nothing is added, not even zero).
//! * **Energy.** The momentum pressure is part of the source pressure the mouth sees, so its work
//!   is in the source work of the ledger, whose closure is unchanged.
//!
//! # The two printed weights
//!
//! * `1`: Eq. (1) p. 401 as printed, derived there and in Fletcher (1976), not fitted.
//! * `2.6`: Thwaites and Fletcher's Sec. II (p. 402, Fig. 3, Eq. (2)-(3)) fit to the pulsating-jet
//!   experiment (3 mm axisymmetric nozzle, 2 m x 10 cm pipe at its 87 Hz resonance, one data set at
//!   `dl` = 8 cm gives C1 = 4.5 and C2 = 3.4 in SI): "a coefficient of 2.6 before the v0 term in
//!   Eq. (1) correctly adjusts the relative magnitudes of the drives for this situation". Their
//!   Sec. IV (p. 406) says the factor is "clearly necessary" for the phase of the jet admittance at
//!   `40 < w/v0 < 100` and might be larger by "perhaps 50%", and their conclusion (p. 408) says
//!   "about 4". The authors could not say whether it reweighs the momentum term or the volume term
//!   through a misread `dl`.
//!
//! # Transfer caveats (also in the provenance sidecar and the site text)
//!
//! * Their jet is turbulent: a 0.5 mm slit, cut-up 1 and 2 cm (cut-up/slit 20 to 40), `Re_h` 430 to
//!   1670. This card's jet is the laminar Bickley kernel of `jet_response`, at the card's own
//!   cut-up/flue-height (a few) and `Re_h`.
//! * Their `dl` (19 to 30 cm) is referred to the pipe cross-section of a 10 cm pipe whose mouth is
//!   about a twentieth of its area; ours, `source_coeff S_p / rho`, is a few millimetres because
//!   the mouth is most of the bore. If 2.6 reweighs through `dl`, the weight equivalent here is far
//!   below 2.6, and below the weight at which the note changes register.
//! * Their measured regime `w/v0` >= 40 per metre; below it the admittance phase went the wrong way
//!   and is unexplained. Eq. (11) (the traveling-wave form) is not what is added here.
//!
//! # Open pipes can latch (observed 2026-11-03; the mechanism is a hypothesis)
//!
//! The steady part of `Q_j - Q_rest` is not removed (no corner frequency for a remover is printed anywhere).
//! A closed pipe cannot respond to it with a steady flow (its steady impedance is infinite), so the
//! pan-pipe renders keep a small static mouth pressure (2 to 8 Pa, 0.4 to 1.8 % of the launched wave's rms) and oscillate normally.
//! An open pipe passes a steady flow, and in the renders the jet then parked at full deflection (`Q_j = q_max`, no oscillation):
//! flute D4 at `kappa` 2.6 (silent, with a `mouth_loss contraction` row) and flute G4 at `kappa` 1 on a card without that row
//! (the first register for about 0.3 s, then silent: mean `Q_j - Q_rest` at 100.0 % of the latch value, +25 Pa at the mouth plane).
//! The same G4 card with only the contraction row added oscillated (65 % of the latch value, +0.01 Pa), one controlled pair.
//! Flute D4 and the measured G4 at `kappa` 1, both with the row, oscillated at 46 to 54 %. The working guess is positive feedback
//! at zero frequency through the open pipe, held back by the mouth's quasi-steady separation loss, so the row raises the weight at
//! which the latch occurs and does not remove it (D4 at 2.6 latched with it); no variant that removes the steady part was run,
//! so the mechanism is not shown.

use super::{FlueError, FlueGeometry, FlueVoice, JetStepCoefficients};

/// Eq. (1) of Thwaites and Fletcher (1983, p. 401) as printed: the momentum term's coefficient.
pub const EQ1_KAPPA: f64 = 1.0;
/// Their Sec. II (p. 402) pulsating-jet fit: "a coefficient of 2.6 before the v0 term".
pub const SEC2_KAPPA: f64 = 2.6;

/// The refusal text of every KAPPA other than the two printed values.
pub const KAPPA_REFUSED: &str = "momentum_drive KAPPA must be 1 (Thwaites and Fletcher 1983 Eq. (1), p. 401) or 2.6 (their Sec. II fit, p. 402); no other value is printed";

/// The drive law, as written in the provenance sidecar.
pub const LAW: &str = "p_mom = kappa rho v0 (Q_j - Q_rest) / S_p in series with the volume drive; v0 = bulk jet speed, S_p = bore area, Q_rest = cut flow of the quiescent jet";

/// Transfer caveats, carried in the sidecar next to the row.
pub const TRANSFER_CAVEATS: &str = "Thwaites and Fletcher's jet is turbulent (0.5 mm slit, cut-up 1 and 2 cm, cut-up/slit 20 to 40, Re_h 430 to 1670); this jet is the laminar Bickley kernel at the card's own cut-up/flue-height and Re_h. Kappa 2.6 is a one-rig fit whose authors could not say whether it reweighs the momentum term or misreads the end correction (theirs 19 to 30 cm referred to the pipe cross-section, ours a few mm): if it reweighs through the end correction the equivalent weight here is far below 2.6. The area is the bore area S_p as printed. The steady part of Q_j - Q_rest is not removed: on an open pipe the jet can park at full deflection and the note goes silent (observed on flute D4 at kappa 2.6 and on a flute G4 card without a mouth_loss row at kappa 1, where adding that row removed it; mechanism unproven), while closed pipes keep a small static mouth pressure and oscillate.";

/// Whether `kappa` is one of the two printed weights (exact values; nothing else is accepted).
#[must_use]
pub fn is_printed_weight(kappa: f64) -> bool {
    kappa.to_bits() == EQ1_KAPPA.to_bits() || kappa.to_bits() == SEC2_KAPPA.to_bits()
}

/// The attached momentum drive: the printed weight and the bore it acts through.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MomentumDrive {
    kappa: f64,
    density_kg_m3: f64,
    bore_area_m2: f64,
}

impl MomentumDrive {
    /// # Errors
    /// A weight other than the two printed ones, or a non-positive or non-finite density or bore.
    pub fn new(kappa: f64, density_kg_m3: f64, bore_area_m2: f64) -> Result<Self, &'static str> {
        if !is_printed_weight(kappa) {
            return Err(KAPPA_REFUSED);
        }
        if !(density_kg_m3.is_finite()
            && density_kg_m3 > 0.0
            && bore_area_m2.is_finite()
            && bore_area_m2 > 0.0)
        {
            return Err("momentum drive requires a positive finite gas density and bore area");
        }
        Ok(Self {
            kappa,
            density_kg_m3,
            bore_area_m2,
        })
    }

    /// The printed weight (1 or 2.6).
    #[must_use]
    pub const fn kappa(&self) -> f64 {
        self.kappa
    }

    /// The bore area `S_p` [m^2] the jet's momentum flux acts over.
    #[must_use]
    pub const fn bore_area_m2(&self) -> f64 {
        self.bore_area_m2
    }

    /// Where the weight is printed.
    #[must_use]
    pub fn source(&self) -> &'static str {
        if self.kappa.to_bits() == EQ1_KAPPA.to_bits() {
            "Thwaites and Fletcher 1983 Eq. (1), p. 401: the momentum term as printed"
        } else {
            "Thwaites and Fletcher 1983 Sec. II, p. 402 (Fig. 3, Eq. (2)-(3)): the pulsating-jet fit, one rig"
        }
    }

    /// Series pressure per unit cut-flow deviation [Pa s / m^3] at jet speed `v0`: `kappa rho v0 / S_p`.
    /// It is the term's direct (no `1/dt`) slope with respect to `Q_j` in the same-sample solve.
    #[must_use]
    pub fn pressure_per_flow(&self, jet_speed_m_s: f64) -> f64 {
        self.kappa * self.density_kg_m3 * jet_speed_m_s / self.bore_area_m2
    }

    /// `p_mom = kappa rho v0 (Q_j - Q_rest) / S_p` [Pa].
    #[must_use]
    pub fn pressure_pa(&self, jet_speed_m_s: f64, cut_flow_m3_s: f64, rest_flow_m3_s: f64) -> f64 {
        self.pressure_per_flow(jet_speed_m_s) * (cut_flow_m3_s - rest_flow_m3_s)
    }
}

/// The cut flow of the quiescent jet (labium displacement zero) at `q_max`: the steady value the
/// momentum drive is measured from. The expression is the one `jet_trial` evaluates at `eta = 0`.
#[must_use]
pub fn rest_cut_flow_m3_s(q_max: f64, offset_m: f64, half_width_m: f64) -> f64 {
    0.5 * q_max * (1.0 + ((-offset_m) / half_width_m).tanh())
}

impl FlueVoice {
    /// Attach the momentum drive after the forced jet: `kappa` is exactly 1 (Thwaites and Fletcher
    /// 1983 Eq. (1), p. 401) or 2.6 (their Sec. II fit, p. 402). See the module doc for the law,
    /// the steady value it acts on and the transfer caveats.
    ///
    /// # Errors
    /// A poisoned or warm voice, a repeated attachment, no forced jet (the drive acts on the
    /// forced jet's cut flow), an attached Nolle/Price transfer, or a weight other than the two
    /// printed ones.
    pub fn with_momentum_drive(mut self, kappa: f64) -> Result<Self, FlueError> {
        if self.poisoned
            || self.samples != 0
            || self.momentum.is_some()
            || self.forced.is_none()
            || self.nolle.is_some()
        {
            return Err(FlueError::Invalid {
                what: "momentum drive requires one cold attachment after the forced jet (jet_response) and no Nolle transfer",
            });
        }
        let drive = MomentumDrive::new(kappa, self.rho, self.geometry.bore_area_m2)
            .map_err(|what| FlueError::Invalid { what })?;
        self.island.provenance.push_str(&format!(
            "; momentum drive: kappa {} ({}), p = kappa rho v0 (Q_j - Q_rest)/S_p with S_p = {:.4e} m^2 added to the volume drive; transfer caveats: {}",
            drive.kappa(),
            drive.source(),
            drive.bore_area_m2(),
            TRANSFER_CAVEATS
        ));
        self.momentum = Some(drive);
        Ok(self)
    }

    /// The attached momentum drive, if any.
    #[must_use]
    pub const fn momentum_drive(&self) -> Option<&MomentumDrive> {
        self.momentum.as_ref()
    }

    /// Bulk jet speed of a trial: the `v0` of `q_max = v0 h W`. Behind a forced jet with a named
    /// spreading law `q_max` is the station's volume flux, `(b(X) / b0)^(1/2)` times the exit's,
    /// and the speed is the station's centre speed, `(b0 / b(X))^(1/2)` times the exit's, so the
    /// momentum flux the cut swings stays the jet's own, which spreading conserves.
    fn momentum_speed_m_s(&self, c: &JetStepCoefficients) -> f64 {
        let bulk = c.q_max / (self.geometry.flue_height_m * self.geometry.flue_width_m);
        match self.forced.as_ref().and_then(super::forced_jet::ForcedJetKernel::station_spread) {
            None => bulk,
            Some((beta, _)) => bulk / beta,
        }
    }

    /// The series pressure the cut jet delivers for cut flow `cut_flow_m3_s`: the volume drive
    /// `source_coeff dQ_j/dt` and, only when the momentum drive is attached, its `p_mom`. Without the row this
    /// is the original expression and nothing is added.
    pub(super) fn jet_source_pa(&self, cut_flow_m3_s: f64, c: JetStepCoefficients) -> f64 {
        let volume = self.source_coeff * (cut_flow_m3_s - self.flow_prev) / self.dt;
        match &self.momentum {
            None => volume,
            Some(drive) => {
                let speed = self.momentum_speed_m_s(&c);
                volume
                    + drive.pressure_pa(
                        speed,
                        cut_flow_m3_s,
                        rest_cut_flow_m3_s(c.q_max, c.offset_m, c.half_width_m),
                    )
            }
        }
    }

    /// The momentum drive's direct slope `kappa rho v0 / S_p` [Pa s / m^3] of the series source with
    /// respect to `Q_j`; exactly zero without the row.
    pub(super) fn momentum_source_slope(&self, c: JetStepCoefficients) -> f64 {
        match &self.momentum {
            None => 0.0,
            Some(drive) => drive.pressure_per_flow(self.momentum_speed_m_s(&c)),
        }
    }
}

#[cfg(test)]
mod tests;
