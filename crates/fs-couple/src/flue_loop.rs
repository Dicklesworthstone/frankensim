//! Flue playing loop (music bead `frankensim-music-v8-root-3ez8g.10.3`):
//! the jet-drive island x the wind line, composed — no `Flute` type, a
//! composition recipe. The island is the aerodynamic source of every
//! flue instrument (recorder, flute, flue organ pipe): a planar jet
//! leaves the flue, crosses the mouth, and is cut by the labium; the
//! pipe's own acoustic field at the mouth deflects the jet at the flue
//! exit, the deflection convects to the labium, and the cut jet's
//! volume flow drives the pipe. Everything the island knows about the
//! JET comes from a jet card minted by the aeroacoustic lab
//! ([`fs_aeroac::jetcard`]); everything it knows about the PIPE comes
//! from the typed fingering chart the wind line already consumes.
//!
//! NO CARD -> REFUSE. A voice cannot be built without a jet card, and
//! only a card whose claim class can drive a tone is admitted (the
//! tonal interim card; a broadband-refusal-boundary card says "no
//! broadband exists below X" and drives nothing). The card's class caps
//! the claims a fixture may mint: a tonal interim card supports
//! edge-tone-class claims only, and [`FlueVoice::claim_check`] refuses
//! any block that ran outside the card's Reynolds validity region.
//!
//! PITCH IS NEVER ASSIGNED. The control surface is blowing pressure
//! and geometry gestures (cut-up, labium offset, jet angle); there is
//! no oscillator, no target note, no frequency input. The note that
//! sounds is the lock the loop finds between the jet's convective delay
//! and the pipe's impedance peaks, and the overblowing ladder emerges
//! from the delay shrinking as the jet speeds up.
//!
//! Model (Fletcher-Rossing jet drive, time domain):
//! - the jet speed is Bernoulli from blowing pressure and the medium,
//!   `U = sqrt(2 p_blow / rho)`; the card is SHAPE/SCALING authority
//!   only (the lab works in lattice units), so every dimensional number
//!   comes from the geometry chart and the medium;
//! - the transverse acoustic velocity at the mouth, `v = q / S_m`,
//!   displaces the jet root by its time integral (receptivity is a
//!   displacement, so its `1/omega` cancels the drive derivative's
//!   `omega` and the pipe's resonance alone selects the mode); the
//!   displacement convects to the labium in
//!   `tau = W / c_p` with `c_p / U` DERIVED from the card's stage-I lock
//!   (`St = f * 2b / U` at `h / delta`, phase condition
//!   `f * h / c_p = n + 1/4`), and grows by an authored spatial gain
//!   inside the sinuous mode's amplification band, which ends at the
//!   tanh shear layer's neutral frequency `omega_c = U / (4 theta)`
//!   (Michalke; `theta` from the card's profile ratio) and opens
//!   linearly from zero at low Strouhal (a high-pass at half the neutral
//!   frequency): the loop cannot lock above or far below the band, and
//!   the band rises with the jet speed;
//! - the jet is cut by the labium: the volume flow entering the pipe is
//!   `Q_j = Q_max * (1 + tanh((eta - y0) / b)) / 2`, the card's
//!   smoothed top-hat profile integrated, `b` its half-width scaled to
//!   the flue height and the card's `theta / b` recorded;
//! - the drive enters the pipe as the series pressure source
//!   `p_j = rho (Delta_m / S_m) dQ_j/dt` in a loop with the mouth
//!   inertance and radiation resistance and the pipe's characteristic
//!   line (the D17 wind line with its fingering lift), solved
//!   implicitly each sample;
//! - the energy ledger is exact for the discretization: source work =
//!   stored mouth inertance + numerical dissipation of the implicit step
//!   + radiation loss + work into the pipe, closure reported per block.
//!
//! Honesty ladder: [`JetIsland::provenance`] states which numbers are
//! card-backed (convection ratio, profile ratio, validity, claim class)
//! and which are authored (spatial gain, mouth radiation resistance,
//! onset seed amplitude); with authored gain the loop is Estimate.
//! Noise injection follows the card: the tonal interim card carries no
//! NoiseTable rows, so no noise SHAPE is injected; the onset seed is
//! the lab's own practice (a symmetric rig never leaves its unstable
//! equilibrium without one) and is keyed by logical voice identity
//! through a Philox stream, so replay is independent of scheduling.

use fs_aeroac::jetcard::{JetCard, JetCardClaim};
use fs_duct::{FingeringTable, Termination};
use fs_material::gas::GasState;
use fs_rand::StreamKey;

use crate::wind_line::{WindLineBank, WindLineError};

/// Philox kernel id of the flue onset seed stream (logical identity;
/// never a scheduling artifact).
pub const FLUE_SEED_KERNEL: u32 = 0x666c_7565; // "flue"

/// Typed refusals of the flue loop.
#[derive(Debug)]
pub enum FlueError {
    /// No jet card was supplied: the island cannot be parameterized.
    NoCard,
    /// The card's claim class cannot drive a tone (a refusal-boundary
    /// card records where broadband does NOT exist).
    CardClassCannotDrive {
        /// The card's claim kind.
        kind: &'static str,
    },
    /// A geometry, medium, or control input is not admissible.
    Invalid {
        /// What was refused.
        what: &'static str,
    },
    /// The wind line refused.
    Line(WindLineError),
    /// A claim was requested for blocks that ran outside the card's
    /// validity region (or with the wrong claim class).
    OutsideCardValidity {
        /// Blocks that ran outside the region.
        blocks: usize,
        /// Lowest and highest jet Reynolds number seen.
        reynolds_seen: (f64, f64),
        /// The card's admitted band.
        reynolds_card: (f64, f64),
    },
}

impl core::fmt::Display for FlueError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            FlueError::NoCard => write!(f, "no jet card: the flue island refuses without a card"),
            FlueError::CardClassCannotDrive { kind } => {
                write!(f, "jet card class `{kind}` cannot drive a tone")
            }
            FlueError::Invalid { what } => write!(f, "invalid flue input: {what}"),
            FlueError::Line(e) => write!(f, "wind line: {e:?}"),
            FlueError::OutsideCardValidity {
                blocks,
                reynolds_seen,
                reynolds_card,
            } => write!(
                f,
                "{blocks} block(s) ran at jet Re in [{:.1}, {:.1}] outside the card's [{:.1}, {:.1}]",
                reynolds_seen.0, reynolds_seen.1, reynolds_card.0, reynolds_card.1
            ),
        }
    }
}

impl From<WindLineError> for FlueError {
    fn from(e: WindLineError) -> Self {
        FlueError::Line(e)
    }
}

/// The jet island's reduced parameters plus their provenance labels.
#[derive(Debug, Clone, PartialEq)]
pub struct JetIsland {
    /// Claim class of the card that minted this island.
    pub claim_kind: &'static str,
    /// Card: locked Strouhal `f * 2b / U` of the recorded stage.
    pub locked_strouhal: f64,
    /// Card: the recorded edge-tone stage.
    pub stage: u8,
    /// Card: labium distance over jet width at the recorded lock.
    pub h_over_delta: f64,
    /// DERIVED from the card: convection speed over centerline speed,
    /// `c_p / U = St * (h / delta) / (stage + 1/4)`.
    pub convection_ratio: f64,
    /// Card: momentum thickness over slot half-width (profile shape).
    pub theta_over_b: f64,
    /// Card: admitted jet Reynolds band (`U * 2b / nu`, dimensionless
    /// so it transfers from the lattice).
    pub reynolds_band: (f64, f64),
    /// Card: amplitude claims require a seeded run.
    pub seeded_amplitude_claims_only: bool,
    /// AUTHORED: spatial amplification of the jet deflection over the
    /// cut-up (the card carries no receptivity gain measurement).
    pub displacement_gain: f64,
    /// AUTHORED: onset seed amplitude relative to the jet speed.
    pub seed_amplitude_rel: f64,
    /// The provenance sentence.
    pub provenance: String,
}

impl JetIsland {
    /// Parameterize the island from a jet card.
    ///
    /// # Errors
    /// [`FlueError::CardClassCannotDrive`] unless the card's claim is
    /// tonal; [`FlueError::Invalid`] on a non-admissible authored gain.
    pub fn from_card(
        card: &JetCard,
        displacement_gain: f64,
        seed_amplitude_rel: f64,
    ) -> Result<JetIsland, FlueError> {
        let JetCardClaim::EdgeToneTonal { feedback } = &card.claim else {
            return Err(FlueError::CardClassCannotDrive {
                kind: card.claim.kind(),
            });
        };
        if !(displacement_gain > 0.0 && displacement_gain.is_finite()) {
            return Err(FlueError::Invalid {
                what: "displacement gain must be positive finite",
            });
        }
        if !(seed_amplitude_rel >= 0.0 && seed_amplitude_rel.is_finite()) {
            return Err(FlueError::Invalid {
                what: "seed amplitude must be finite and non-negative",
            });
        }
        let stage = feedback.stage;
        let phase_cycles = f64::from(stage) + 0.25;
        let convection_ratio = feedback.locked_strouhal * card.validity.h_over_delta / phase_cycles;
        if !(convection_ratio > 0.0 && convection_ratio < 1.0) {
            return Err(FlueError::Invalid {
                what: "card lock implies a convection ratio outside (0, 1)",
            });
        }
        let theta_over_b = card.profile.momentum_thickness / card.profile.slot_half;
        Ok(JetIsland {
            claim_kind: card.claim.kind(),
            locked_strouhal: feedback.locked_strouhal,
            stage,
            h_over_delta: card.validity.h_over_delta,
            convection_ratio,
            theta_over_b,
            reynolds_band: (card.validity.reynolds_lo, card.validity.reynolds_hi),
            seeded_amplitude_claims_only: card.validity.seeded_amplitude_claims_only,
            displacement_gain,
            seed_amplitude_rel,
            provenance: format!(
                "card-backed: claim class {} (stage {stage} lock St {:.5} at h/delta {:.1} -> c_p/U {:.4}), theta/b {:.4} (amplification band omega_c = U/(4 theta), Michalke tanh neutral), Re band [{:.1}, {:.1}], {} receipts; authored: displacement gain {displacement_gain}, onset seed {seed_amplitude_rel} x U, mouth radiation resistance; medium and geometry from the chart",
                card.claim.kind(),
                feedback.locked_strouhal,
                card.validity.h_over_delta,
                convection_ratio,
                theta_over_b,
                card.validity.reynolds_lo,
                card.validity.reynolds_hi,
                card.provenance.receipts.len()
            ),
        })
    }
}

/// Mouth geometry gestures from the chart (all SI).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlueGeometry {
    /// Flue exit width along the labium [m] (jet breadth).
    pub flue_width_m: f64,
    /// Flue exit height [m] (jet thickness `2b`).
    pub flue_height_m: f64,
    /// Cut-up: flue exit to labium distance `W` [m].
    pub cut_up_m: f64,
    /// Labium offset from the jet centerline [m] (positive = the labium
    /// sits on the pipe side of the undeflected jet).
    pub labium_offset_m: f64,
    /// Jet angle towards the pipe side [rad]; tilts the effective
    /// offset by `W tan(angle)`.
    pub jet_angle_rad: f64,
    /// Pipe bore area at the mouth [m^2] (the wind line's `zc` basis).
    pub bore_area_m2: f64,
    /// AUTHORED: mouth radiation resistance as a fraction of the mouth's
    /// characteristic impedance `rho c / S_m`.
    pub mouth_loss_rel: f64,
}

impl FlueGeometry {
    /// Mouth opening area `S_m = width * cut-up`.
    #[must_use]
    pub fn mouth_area_m2(&self) -> f64 {
        self.flue_width_m * self.cut_up_m
    }
    /// Effective labium offset after the jet-angle tilt.
    #[must_use]
    pub fn effective_offset_m(&self) -> f64 {
        self.labium_offset_m + self.cut_up_m * self.jet_angle_rad.tan()
    }
    /// Mouth end correction `0.61 * equivalent radius` (unflanged
    /// opening law; the mouth inertance basis).
    #[must_use]
    pub fn end_correction_m(&self) -> f64 {
        0.61 * (self.mouth_area_m2() / core::f64::consts::PI).sqrt()
    }
    fn validate(&self) -> Result<(), FlueError> {
        let ok = [
            self.flue_width_m,
            self.flue_height_m,
            self.cut_up_m,
            self.bore_area_m2,
        ]
        .iter()
        .all(|v| *v > 0.0 && v.is_finite())
            && self.labium_offset_m.is_finite()
            && self.jet_angle_rad.is_finite()
            && self.jet_angle_rad.abs() < core::f64::consts::FRAC_PI_2
            && self.mouth_loss_rel >= 0.0
            && self.mouth_loss_rel.is_finite();
        if ok {
            Ok(())
        } else {
            Err(FlueError::Invalid {
                what: "flue geometry must be positive finite (angle within +-pi/2)",
            })
        }
    }
}

/// Control-rate inputs (applied between blocks).
#[derive(Debug, Clone, Copy)]
pub enum FlueControl {
    /// Steady blowing pressure [Pa].
    SetBlowingPressure(f64),
    /// Switch fingering on the char image with a crossfade.
    Fingering {
        /// Fingering index in the chart.
        index: usize,
        /// Crossfade length [samples].
        fade_samples: usize,
    },
}

/// Per-block diagnostics (the observer's honesty block).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlueBlockDiag {
    /// Block index.
    pub block: u64,
    /// Blowing pressure [Pa].
    pub blow_pa: f64,
    /// Jet centerline speed [m/s].
    pub jet_speed_m_s: f64,
    /// Jet Reynolds number `U * 2b / nu`.
    pub reynolds: f64,
    /// Whether this block ran inside the card's Reynolds band.
    pub inside_card_validity: bool,
    /// Jet transit delay [s].
    pub transit_s: f64,
    /// Zero-crossing lock estimate of the wave launched into the pipe [Hz].
    pub lock_hz: f64,
    /// RMS of the wave launched into the pipe [Pa].
    pub p_rms_pa: f64,
    /// Fraction of samples the labium cut on the pipe side.
    pub duty_pipe_side: f64,
    /// RMS jet displacement at the labium over the jet half-width
    /// (above about 1 the labium cut saturates: a limit cycle).
    pub eta_rms_over_b: f64,
    /// Work the blowing pressure did on the jet flow [J].
    pub blow_work_j: f64,
    /// Work the jet drive source did on the mouth loop [J].
    pub source_work_j: f64,
    /// Work delivered into the pipe [J].
    pub pipe_work_j: f64,
    /// Radiation loss at the mouth [J].
    pub radiation_loss_j: f64,
    /// Change of stored mouth inertance energy [J].
    pub stored_delta_j: f64,
    /// Numerical dissipation of the implicit mouth step [J].
    pub numerical_dissipation_j: f64,
    /// Ledger closure `source - (stored + numerical + radiation + pipe)`
    /// [J]; roundoff by construction.
    pub ledger_defect_j: f64,
}

/// One flue voice: the island, its mouth loop, and its wind line.
pub struct FlueVoice {
    island: JetIsland,
    geometry: FlueGeometry,
    line: WindLineBank,
    rho: f64,
    nu: f64,
    zc: f64,
    dt: f64,
    blow_pa: f64,
    // Mouth transverse velocity history (ring, oldest overwritten).
    history: Vec<f64>,
    head: usize,
    q_prev: f64,
    flow_prev: f64,
    // Jet amplification band: two cascaded first-order sections at the
    // shear layer's neutral frequency (state of each section).
    band_state: [f64; 2],
    // Low-frequency edge of the band: growth vanishes linearly at small
    // Strouhal, a first-order high-pass (previous input, previous output).
    band_hp: [f64; 2],
    momentum_thickness_m: f64,
    // Jet root displacement: the leaky integral of the band-limited,
    // delayed mouth velocity (receptivity is displacement, so the
    // integral's 1/omega cancels the drive derivative's omega and the
    // pipe's resonance alone selects the mode).
    displacement_state: f64,
    // Mouth loop coefficients.
    inertance: f64,
    resistance: f64,
    source_coeff: f64,
    seed: fs_rand::Stream,
    seed_key: StreamKey,
    block_index: u64,
    diags: Vec<FlueBlockDiag>,
}

impl FlueVoice {
    /// Compose a voice. `card` is REQUIRED: `None` is the typed refusal
    /// the menu law demands.
    ///
    /// # Errors
    /// [`FlueError`] on a missing or non-driving card, an inadmissible
    /// geometry, or a wind line refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        card: Option<&JetCard>,
        displacement_gain: f64,
        seed_amplitude_rel: f64,
        geometry: FlueGeometry,
        table: &FingeringTable,
        gas: &GasState,
        termination: Termination,
        sample_rate_hz: u32,
        voice_seed: u64,
        voice_index: u32,
    ) -> Result<FlueVoice, FlueError> {
        let card = card.ok_or(FlueError::NoCard)?;
        let island = JetIsland::from_card(card, displacement_gain, seed_amplitude_rel)?;
        geometry.validate()?;
        if sample_rate_hz == 0 {
            return Err(FlueError::Invalid {
                what: "sample rate must be positive",
            });
        }
        let rho = gas.density;
        let c = gas.sound_speed;
        let nu = gas.dynamic_viscosity / rho;
        if !(rho > 0.0 && c > 0.0 && nu > 0.0) {
            return Err(FlueError::Invalid {
                what: "medium must have positive density, sound speed, viscosity",
            });
        }
        let zc = rho * c / geometry.bore_area_m2;
        let line = WindLineBank::new(table, gas, termination, sample_rate_hz, zc, false)?;
        let s_m = geometry.mouth_area_m2();
        let inertance = rho * geometry.end_correction_m() / s_m;
        let resistance = geometry.mouth_loss_rel * rho * c / s_m;
        let dt = 1.0 / f64::from(sample_rate_hz);
        // History long enough for the slowest admissible jet (1 m/s).
        let max_delay = geometry.cut_up_m / (island.convection_ratio * 1.0);
        let history_len = ((max_delay / dt).ceil() as usize).max(2) + 2;
        let seed_key = StreamKey {
            seed: voice_seed,
            kernel: FLUE_SEED_KERNEL,
            tile: voice_index,
        };
        let momentum_thickness_m = island.theta_over_b * 0.5 * geometry.flue_height_m;
        Ok(FlueVoice {
            island,
            geometry,
            line,
            rho,
            nu,
            zc,
            dt,
            blow_pa: 0.0,
            history: vec![0.0; history_len],
            head: 0,
            q_prev: 0.0,
            flow_prev: 0.0,
            band_state: [0.0; 2],
            band_hp: [0.0; 2],
            momentum_thickness_m,
            displacement_state: 0.0,
            inertance,
            resistance,
            source_coeff: rho * geometry.end_correction_m() / s_m,
            seed: seed_key.stream(),
            seed_key,
            block_index: 0,
            diags: Vec::new(),
        })
    }

    /// The island's parameters and provenance.
    #[must_use]
    pub fn island(&self) -> &JetIsland {
        &self.island
    }

    /// The seed stream identity (logical voice identity).
    #[must_use]
    pub fn seed_key(&self) -> StreamKey {
        self.seed_key
    }

    /// Per-block diagnostics so far.
    #[must_use]
    pub fn diagnostics(&self) -> &[FlueBlockDiag] {
        &self.diags
    }

    /// Apply a control-rate input.
    ///
    /// # Errors
    /// [`FlueError`] on a non-finite or negative pressure, or an unknown
    /// fingering.
    pub fn apply(&mut self, control: FlueControl) -> Result<(), FlueError> {
        match control {
            FlueControl::SetBlowingPressure(p) => {
                if !(p.is_finite() && p >= 0.0) {
                    return Err(FlueError::Invalid {
                        what: "blowing pressure must be finite and non-negative",
                    });
                }
                self.blow_pa = p;
            }
            FlueControl::Fingering {
                index,
                fade_samples,
            } => self.line.switch_fingering(index, fade_samples)?,
        }
        Ok(())
    }

    /// Jet centerline speed for the current blowing pressure.
    #[must_use]
    pub fn jet_speed_m_s(&self) -> f64 {
        (2.0 * self.blow_pa / self.rho).sqrt()
    }

    /// Jet Reynolds number `U * 2b / nu` for the current blowing pressure.
    #[must_use]
    pub fn reynolds(&self) -> f64 {
        self.jet_speed_m_s() * self.geometry.flue_height_m / self.nu
    }

    /// Jet transit delay `W / c_p` [s] (infinite at zero blowing).
    #[must_use]
    pub fn transit_s(&self) -> f64 {
        self.geometry.cut_up_m / (self.island.convection_ratio * self.jet_speed_m_s())
    }

    fn delayed_velocity(&self, delay_samples: f64) -> f64 {
        let len = self.history.len();
        let clamped = delay_samples.min((len - 2) as f64).max(0.0);
        let whole = clamped.floor();
        let frac = clamped - whole;
        let idx = |back: usize| self.history[(self.head + len - 1 - back) % len];
        let k = whole as usize;
        (1.0 - frac) * idx(k) + frac * idx(k + 1)
    }

    /// Advance one block, writing the wave launched into the pipe
    /// (`p_plus` at the mouth plane) into `out`.
    ///
    /// # Errors
    /// [`FlueError`] on an empty block or a wind line refusal.
    pub fn step_block(&mut self, out: &mut [f64]) -> Result<(), FlueError> {
        if out.is_empty() {
            return Err(FlueError::Invalid {
                what: "empty block",
            });
        }
        let dt = self.dt;
        let u = self.jet_speed_m_s();
        let b = 0.5 * self.geometry.flue_height_m;
        let q_max = u * self.geometry.flue_height_m * self.geometry.flue_width_m;
        let s_m = self.geometry.mouth_area_m2();
        let y0 = self.geometry.effective_offset_m();
        let transit = if u > 0.0 {
            self.transit_s()
        } else {
            f64::INFINITY
        };
        let delay_samples = if transit.is_finite() {
            transit / dt
        } else {
            0.0
        };
        let seed_scale = self.island.seed_amplitude_rel * u;
        // Amplification band of the sinuous jet mode: the tanh shear
        // layer is spatially amplified only below its neutral frequency
        // omega_c = U / (4 theta) (Michalke: neutral wavenumber
        // alpha theta = 1/2 convecting at U/2), so the deflection path is
        // band-limited by two first-order sections at omega_c. The
        // cutoff rises with the jet speed: that, with the shrinking
        // transit, is the overblowing mechanism.
        let omega_c = if u > 0.0 {
            u / (4.0 * self.momentum_thickness_m)
        } else {
            0.0
        };
        let band_alpha = omega_c * dt / (1.0 + omega_c * dt);
        // Low edge of the band: the growth rate vanishes linearly below
        // the most amplified Strouhal (about half the neutral one), so a
        // first-order high-pass at omega_c / 2 carries that slope. Without
        // it the integral receptivity has its largest gain at the lowest
        // frequencies, where the open pipe offers no impedance, and the
        // loop falls into a slow relaxation cycle instead of a note
        // (executed: 6 Hz at every blowing pressure).
        let omega_h = 0.5 * omega_c;
        let hp_alpha = 1.0 / (1.0 + omega_h * dt);
        // Displacement leak two decades below the amplification band.
        let leak_rate = omega_c / 100.0;
        let (l, r, zc) = (self.inertance, self.resistance, self.zc);
        // Implicit mouth loop: p_j = L (q - q_prev)/dt + R q + p_in, with
        // p_in = p_plus + p_minus and q = (p_plus - p_minus)/zc.
        let denom = l / dt + r + zc;
        let mut blow_work = 0.0;
        let mut source_work = 0.0;
        let mut pipe_work = 0.0;
        let mut radiation = 0.0;
        let mut stored_delta = 0.0;
        let mut numerical = 0.0;
        let mut pipe_side = 0usize;
        let mut eta_sq = 0.0;
        for slot in out.iter_mut() {
            // Jet deflection at the labium from the delayed mouth velocity
            // plus the voice-keyed onset seed.
            let v_raw = self.delayed_velocity(delay_samples);
            self.band_state[0] += band_alpha * (v_raw - self.band_state[0]);
            self.band_state[1] += band_alpha * (self.band_state[0] - self.band_state[1]);
            let hp_in = self.band_state[1];
            let hp_out = hp_alpha * (self.band_hp[1] + hp_in - self.band_hp[0]);
            self.band_hp = [hp_in, hp_out];
            let v_delayed = hp_out;
            let seed = if seed_scale > 0.0 {
                seed_scale * (2.0 * self.seed.next_f64() - 1.0)
            } else {
                0.0
            };
            // Receptivity: the jet root is displaced by the integral of the
            // transverse velocity (with a leak far below the band to keep
            // the mean at zero); the seed enters as transverse velocity.
            self.displacement_state +=
                dt * (v_delayed + seed - leak_rate * self.displacement_state);
            let eta = self.island.displacement_gain * self.displacement_state;
            eta_sq += eta * eta;
            let flow = if u > 0.0 {
                0.5 * q_max * (1.0 + ((eta - y0) / b).tanh())
            } else {
                0.0
            };
            if eta > y0 {
                pipe_side += 1;
            }
            blow_work += self.blow_pa * flow * dt;
            let p_source = self.source_coeff * (flow - self.flow_prev) / dt;
            self.flow_prev = flow;
            let p_minus = self.line.incoming();
            let q = (p_source + (l / dt) * self.q_prev - 2.0 * p_minus) / denom;
            let p_plus = p_minus + zc * q;
            let p_in = p_plus + p_minus;
            let _ = self.line.push(p_plus)?;
            // Ledger terms of the implicit step (exact identity).
            source_work += p_source * q * dt;
            pipe_work += p_in * q * dt;
            radiation += r * q * q * dt;
            stored_delta += 0.5 * l * (q * q - self.q_prev * self.q_prev);
            numerical += 0.5 * l * (q - self.q_prev) * (q - self.q_prev);
            self.q_prev = q;
            // Mouth transverse velocity into the history ring.
            self.history[self.head] = q / s_m;
            self.head = (self.head + 1) % self.history.len();
            // The observable is the wave launched into the pipe: the mouth
            // plane is the pipe's pressure node, so the total pressure
            // there says little about the standing wave the pipe carries.
            *slot = p_plus;
        }
        let n = out.len() as f64;
        let mean = out.iter().sum::<f64>() / n;
        let mut crossings = 0usize;
        let mut prev = out[0] - mean;
        let mut sq = 0.0;
        for &p in out.iter() {
            let v = p - mean;
            sq += v * v;
            if prev < 0.0 && v >= 0.0 {
                crossings += 1;
            }
            prev = v;
        }
        let reynolds = self.reynolds();
        let (lo, hi) = self.island.reynolds_band;
        self.diags.push(FlueBlockDiag {
            block: self.block_index,
            blow_pa: self.blow_pa,
            jet_speed_m_s: u,
            reynolds,
            inside_card_validity: reynolds >= lo && reynolds <= hi,
            transit_s: transit,
            lock_hz: crossings as f64 / (n * dt),
            p_rms_pa: (sq / n).sqrt(),
            duty_pipe_side: pipe_side as f64 / n,
            eta_rms_over_b: (eta_sq / n).sqrt() / b,
            blow_work_j: blow_work,
            source_work_j: source_work,
            pipe_work_j: pipe_work,
            radiation_loss_j: radiation,
            stored_delta_j: stored_delta,
            numerical_dissipation_j: numerical,
            ledger_defect_j: source_work - (stored_delta + numerical + radiation + pipe_work),
        });
        self.block_index += 1;
        Ok(())
    }

    /// Whether every block so far ran inside the card's validity region:
    /// the gate a fixture must pass before minting any claim from this
    /// voice's output.
    ///
    /// # Errors
    /// [`FlueError::OutsideCardValidity`] naming the excursion.
    pub fn claim_check(&self) -> Result<(), FlueError> {
        let outside = self
            .diags
            .iter()
            .filter(|d| !d.inside_card_validity)
            .count();
        if outside == 0 {
            return Ok(());
        }
        let seen = self.diags.iter().fold((f64::INFINITY, 0.0f64), |acc, d| {
            (acc.0.min(d.reynolds), acc.1.max(d.reynolds))
        });
        Err(FlueError::OutsideCardValidity {
            blocks: outside,
            reynolds_seen: seen,
            reynolds_card: self.island.reynolds_band,
        })
    }
}

/// An organ rank: N independent voices from one chest chart. A stop is
/// which pipes receive the jet card; kernels (chart, medium, card) are
/// shared, state is per voice, and there is no shared wind stepper.
pub struct FlueRank {
    voices: Vec<FlueVoice>,
}

impl FlueRank {
    /// Build `n` voices with consecutive voice indices under one seed.
    ///
    /// # Errors
    /// [`FlueError`] from any voice.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        card: Option<&JetCard>,
        displacement_gain: f64,
        seed_amplitude_rel: f64,
        geometries: &[FlueGeometry],
        table: &FingeringTable,
        gas: &GasState,
        termination: Termination,
        sample_rate_hz: u32,
        rank_seed: u64,
    ) -> Result<FlueRank, FlueError> {
        let mut voices = Vec::with_capacity(geometries.len());
        for (index, geometry) in geometries.iter().enumerate() {
            voices.push(FlueVoice::new(
                card,
                displacement_gain,
                seed_amplitude_rel,
                *geometry,
                table,
                gas,
                termination,
                sample_rate_hz,
                rank_seed,
                u32::try_from(index).map_err(|_| FlueError::Invalid {
                    what: "rank too large",
                })?,
            )?);
        }
        Ok(FlueRank { voices })
    }

    /// Mutable access to the voices (per-voice control).
    pub fn voices_mut(&mut self) -> &mut [FlueVoice] {
        &mut self.voices
    }

    /// The voices.
    #[must_use]
    pub fn voices(&self) -> &[FlueVoice] {
        &self.voices
    }

    /// Advance every voice one block in the given order, one output
    /// buffer per voice. Order is a scheduling choice and must not
    /// change any voice's output.
    ///
    /// # Errors
    /// [`FlueError`] from any voice.
    pub fn step_block(&mut self, order: &[usize], outs: &mut [Vec<f64>]) -> Result<(), FlueError> {
        if outs.len() != self.voices.len() || order.len() != self.voices.len() {
            return Err(FlueError::Invalid {
                what: "one output buffer and one order slot per voice",
            });
        }
        for &v in order {
            let Some(voice) = self.voices.get_mut(v) else {
                return Err(FlueError::Invalid {
                    what: "order names a voice outside the rank",
                });
            };
            voice.step_block(&mut outs[v])?;
        }
        Ok(())
    }
}
