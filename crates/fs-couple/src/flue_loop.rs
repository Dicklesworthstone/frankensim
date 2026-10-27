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
//! from the typed fingering chart or the explicitly supplied regional tube graph.
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
//! - the default jet speed is Bernoulli from blowing pressure and the medium,
//!   `U = sqrt(2 p_blow / rho)`. An explicitly supplied rectangular windway
//!   adds the shared material-gas viscous slit pressure drop before computing
//!   the bulk jet speed, injected flow and convection. The unsteady windway
//!   option resolves its no-slip momentum profile and convects feedback on an
//!   explicit spatial mesh; source fluid inertia survives pressure release.
//!   The card remains SHAPE/SCALING authority
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
//!   the quasi-steady history is interpolated at the actual retarded time,
//!   with sub-sample transit coupled to the current mouth solve;
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
//!   Optional explicit mouth contraction adds its own same-sample pressure drop
//!   and separate nonnegative separation loss; it is never an output effect.
//!   An explicitly selected rectangular exterior load replaces the legacy
//!   mouth inertance/loss with passive near-field memory and radiation loss.
//!
//! Honesty ladder: [`JetIsland::provenance`] states which numbers are
//! card-backed (convection ratio, profile ratio, validity, claim class)
//! and which are authored (spatial gain, legacy mouth radiation resistance,
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

use crate::wind_line::{WindBoundary, WindLineBank, WindLineError};

/// Physical regional tube networks sharing the same jet-drive loop.
pub mod network;
use network::FlueResonator;

/// Original jet/mouth composition with the shared causal DG port owner.
pub mod dg;

/// Optional geometry/gas-bound separation loss in the same-sample mouth solve.
pub mod mouth_loss;

/// Rectangular exterior-gas reaction sharing the physical mouth receiver.
pub mod mouth_radiation;

/// A wall window's radiation load on the head's outer cylinder.
pub mod window_radiation;

/// Refinable unsteady no-slip source channel and spatial jet feedback transport.
pub mod windway;

/// Opt-in published jet transfer of Nolle/Price (2015) in place of the island's band chain.
pub mod nolle_jet;

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
    /// An existing material, geometry, loss, wave or observation owner refused.
    Physics(crate::acoustic_realize::AcousticRealizeError),
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
            FlueError::Physics(e) => write!(f, "flue physical system: {e:?}"),
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

// A single same-sample boundary for the existing two propagation images.
// The physical graph exposes the already-arrived wave; the TMM bank exposes
// its instantaneous reflection AND the already-determined history remainder.
impl FlueResonator {
    fn mouth_boundary(&self) -> Result<WindBoundary, FlueError> {
        match self {
            Self::Tmm(line) => line.boundary_response().map_err(FlueError::Line),
            Self::Dg(duct) => duct.boundary(),
            Self::Network(_) => Ok(WindBoundary {
                direct_reflection: 0.0, history_pressure_pa: self.incoming(),
            }),
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct MouthStep {
    flow_m3_s: f64,
    outgoing_pa: f64,
    incoming_pa: f64,
    pressure_pa: f64,
}

#[derive(Clone, Copy)]
struct RetardedVelocity {
    past_m_s: f64,
    current_weight: f64,
}

// The newest ring entry is v_(n-1), not v_n. At t_n-tau, a delay below
// one step therefore includes the UNKNOWN current mouth velocity. Return
// that coefficient to the simultaneous port solve instead of adding a step.
fn retarded_velocity(history: &[f64], head: usize, delay_samples: f64)
    -> Result<RetardedVelocity, FlueError>
{
    if history.len() < 2 || head >= history.len() || !delay_samples.is_finite()
        || delay_samples < 0.0 || delay_samples > (history.len()-1) as f64 {
        return Err(FlueError::Invalid { what: "jet convection exceeds its finite retained history" });
    }
    let past = |back: usize| history[(head+history.len()-1-back)%history.len()];
    let value = if delay_samples < 1.0 {
        RetardedVelocity { past_m_s: delay_samples*past(0),
            current_weight: 1.0-delay_samples }
    } else {
        let back = delay_samples-1.0;
        let whole = back.floor() as usize;
        let fraction = back-whole as f64;
        RetardedVelocity { past_m_s: (1.0-fraction)*past(whole)+fraction*past(whole+1),
            current_weight: 0.0 }
    };
    if !value.past_m_s.is_finite() {
        return Err(FlueError::Invalid { what: "retarded mouth velocity is not finite" });
    }
    Ok(value)
}

#[derive(Clone, Copy)]
struct JetStepCoefficients {
    q_max: f64,
    half_width_m: f64,
    offset_m: f64,
    band_alpha: f64,
    hp_alpha: f64,
    leak_rate: f64,
    // The opt-in Nolle/Price transfer supplies the labium displacement directly; None keeps the
    // island's own receptivity chain and every default golden bit-identical.
    eta_override: Option<f64>,
}

#[derive(Clone, Copy)]
struct JetTrial {
    band_state: [f64; 2],
    band_hp: [f64; 2],
    displacement_m: f64,
    eta_m: f64,
    flow_m3_s: f64,
    source_pa: f64,
}

// Series source + inertance + resistance, coupled to an affine wave boundary:
//   source = L*(Q-Qold)/dt + R*Q + Pin
//   Pin = p+ + p-; Zc*Q = p+ - p-; p- = r0*p+ + history.
// All four quantities refer to THIS sample. Never delay r0 or substitute the
// previous reflected sample. Eliminate in wave coordinates so r0=+/-1 remains
// admissible when the complete physical system is nonsingular. This is an
// acoustic port solve, not an instrument-specific oscillator or a second jet law.
fn solve_mouth_port(
    boundary: WindBoundary, source_pa: f64, previous_flow_m3_s: f64,
    inertance: f64, resistance: f64, zc: f64, dt: f64,
) -> Result<MouthStep, FlueError> {
    let refuse = || FlueError::Invalid { what: "mouth port requires finite passive coefficients and a nonsingular wave boundary" };
    if ![boundary.direct_reflection, boundary.history_pressure_pa, source_pa,
        previous_flow_m3_s, inertance, resistance, zc, dt].iter().all(|x| x.is_finite())
        || inertance < 0.0 || resistance < 0.0 || zc <= 0.0 || dt <= 0.0 {
        return Err(refuse());
    }
    let r0 = boundary.direct_reflection;
    let history = boundary.history_pressure_pa;
    let alpha = inertance / dt;
    let impedance = alpha + resistance;
    let denominator = impedance*(1.0-r0) + zc*(1.0+r0);
    if !alpha.is_finite() || !impedance.is_finite()
        || !denominator.is_finite() || denominator <= 0.0 { return Err(refuse()); }
    let drive = source_pa + alpha*previous_flow_m3_s;
    // Preserve the original graph path's arithmetic for its exactly zero
    // direct term; this is the same algebra, not a different physical model.
    let (flow, outgoing) = if r0 == 0.0 {
        let flow = (drive - 2.0*history)/denominator;
        (flow, history + zc*flow)
    } else {
        let outgoing = (zc*drive + (impedance-zc)*history)/denominator;
        let flow = ((1.0-r0)*outgoing-history)/zc;
        (flow, outgoing)
    };
    let incoming = r0*outgoing + history;
    let pressure = outgoing + incoming;
    if ![drive, flow, outgoing, incoming, pressure].iter().all(|x| x.is_finite()) {
        return Err(refuse());
    }
    Ok(MouthStep { flow_m3_s: flow, outgoing_pa: outgoing,
        incoming_pa: incoming, pressure_pa: pressure })
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
                "card-backed: claim class {} (stage {stage} lock St {:.5} at h/delta {:.1} -> c_p/U {:.4}), theta/b {:.4} (amplification band omega_c = U/(4 theta), Michalke tanh neutral), Re band [{:.1}, {:.1}], {} receipts; authored: displacement gain {displacement_gain}, onset seed {seed_amplitude_rel} x U, legacy mouth radiation resistance (unless a geometry/gas load is selected); medium and geometry from the chart",
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
    /// Jet angle towards the pipe side [rad]. The undeflected jet then
    /// reaches the labium `W tan(angle)` further to the pipe side, which
    /// lowers the effective offset by that much.
    pub jet_angle_rad: f64,
    /// Pipe bore area at the mouth [m^2] (the wind line's `zc` basis).
    pub bore_area_m2: f64,
    /// AUTHORED: mouth radiation resistance as a fraction of the mouth's
    /// characteristic impedance `rho c / S_m`. Must be zero when selecting
    /// physical rectangular reaction; the two loads cannot be superimposed.
    pub mouth_loss_rel: f64,
}

impl FlueGeometry {
    /// Mouth opening area `S_m = width * cut-up`.
    #[must_use]
    pub fn mouth_area_m2(&self) -> f64 {
        self.flue_width_m * self.cut_up_m
    }
    /// Effective labium offset after the jet-angle tilt: the labium's
    /// pipe-side distance from the tilted jet's centreline at the labium.
    #[must_use]
    pub fn effective_offset_m(&self) -> f64 {
        self.labium_offset_m - self.cut_up_m * self.jet_angle_rad.tan()
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
    /// Prescribed jet direction [rad] at fixed exit, cut-up and mouth area.
    /// The existing geometric cut law changes the injected flow; all acoustic,
    /// convection and fluid-momentum histories remain live.
    SetJetAngle(f64),
    /// Nonnegative physical closing force at an explicitly bound graph endpoint.
    /// The existing mechanical/contact state evolves on subsequent samples.
    SetBoundaryClosingForce {
        /// Original explicitly driven node address.
        node: usize,
        /// Nonnegative closing force [N].
        force_n: f64,
    },
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
    /// Jet reference speed [m/s], the bulk speed with a viscous windway.
    pub jet_speed_m_s: f64,
    /// Jet Reynolds number `U * 2b / nu`.
    pub reynolds: f64,
    /// Lowest and highest accepted sample Reynolds numbers in this block.
    /// Transient attacks/releases cannot be hidden by the endpoint value.
    pub reynolds_range: (f64, f64),
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
    /// Radiation loss at the mouth [J]. With a rectangular load this is escaped
    /// branch energy, separate from its retained near-field storage.
    pub radiation_loss_j: f64,
    /// Signed exchange into a joint mouth/terminal exterior's mouth channel [J].
    /// Zero on independent self-load paths; this is not dissipated heat.
    pub joint_mouth_work_j: f64,
    /// Irreversible quasi-steady separated-flow loss [J], not radiation.
    /// Zero when the optional mouth contraction is absent.
    pub separation_loss_j: f64,
    /// Change of stored mouth energy [J]: legacy inertance or the selected
    /// rectangular load's retained near-field state.
    pub stored_delta_j: f64,
    /// Numerical dissipation of the implicit mouth step [J].
    pub numerical_dissipation_j: f64,
    /// Ledger closure including the signed joint-mouth exchange [J].
    /// A joint source uses independently checked original DG and MIMO work.
    pub ledger_defect_j: f64,
    /// Total source work into the optional resolved windway [J].
    pub windway_work_j: f64,
    /// Change of the windway fluid's resolved kinetic energy [J].
    pub windway_stored_delta_j: f64,
    /// Windway wall heat [J], separate from the acoustic mouth losses.
    pub windway_viscous_loss_j: f64,
    /// Work through the windway's explicit ideal kinetic exit [J].
    pub windway_exit_work_j: f64,
    /// Windway time-discretization loss [J], not physical wall heat.
    pub windway_numerical_dissipation_j: f64,
    /// Independent windway pressure/storage/loss closure [J].
    pub windway_ledger_defect_j: f64,
}

/// One flue voice: the island, its mouth loop, and its wind line.
pub struct FlueVoice {
    island: JetIsland,
    geometry: FlueGeometry,
    line: FlueResonator,
    gas: GasState,
    mouth_contraction: Option<crate::bernoulli_aperture::contraction::ContractionLoss>,
    mouth_radiation: Option<mouth_radiation::MouthRadiation>,
    window_radiation: Option<window_radiation::WindowRadiation>,
    // The series inertance is the rest of a computed window mass.
    window_mass: bool,
    viscous_windway: Option<crate::bernoulli_aperture::viscous::ViscousSlit>,
    unsteady_windway: Option<windway::UnsteadyWindway>,
    rho: f64,
    nu: f64,
    zc: f64,
    dt: f64,
    blow_pa: f64,
    jet_speed_m_s: f64,
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
    // Opt-in published jet transfer; its source memory commits with the mouth
    // solve, and both feedback flows belong to the same mechanical interval.
    nolle: Option<nolle_jet::NolleJet>,
    // Mouth loop coefficients.
    inertance: f64,
    resistance: f64,
    source_coeff: f64,
    seed: fs_rand::Stream,
    seed_key: StreamKey,
    block_index: u64,
    samples: u64,
    poisoned: bool,
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
        Self::with_resonator(card, displacement_gain, seed_amplitude_rel, geometry,
            gas, sample_rate_hz, voice_seed, voice_index, 1.0, 1_048_576,
            |zc| WindLineBank::new(table, gas, termination, sample_rate_hz, zc, false)
                .map(FlueResonator::Tmm).map_err(FlueError::Line))
    }

    // One jet/mouth initialization for both propagation images. Only the
    // resonator factory changes; no dummy source or replacement jet algorithm.
    #[allow(clippy::too_many_arguments)]
    fn with_resonator(
        card: Option<&JetCard>, displacement_gain: f64, seed_amplitude_rel: f64,
        geometry: FlueGeometry, gas: &GasState, sample_rate_hz: u32,
        voice_seed: u64, voice_index: u32, minimum_jet_speed: f64, max_history: usize,
        make_line: impl FnOnce(f64) -> Result<FlueResonator, FlueError>,
    ) -> Result<Self, FlueError> {
        let card = card.ok_or(FlueError::NoCard)?;
        let island = JetIsland::from_card(card, displacement_gain, seed_amplitude_rel)?;
        geometry.validate()?;
        if sample_rate_hz == 0 {
            return Err(FlueError::Invalid {
                what: "sample rate must be positive",
            });
        }
        gas.validate_coherence().map_err(|_| FlueError::Invalid {
            what: "flue source requires one coherent material gas state",
        })?;
        let rho = gas.density;
        let c = gas.sound_speed;
        let nu = gas.dynamic_viscosity / rho;
        if !(rho > 0.0 && c > 0.0 && nu > 0.0) {
            return Err(FlueError::Invalid {
                what: "medium must have positive density, sound speed, viscosity",
            });
        }
        let zc = rho * c / geometry.bore_area_m2;
        let line = make_line(zc)?;
        let s_m = geometry.mouth_area_m2();
        let inertance = rho * geometry.end_correction_m() / s_m;
        let resistance = geometry.mouth_loss_rel * rho * c / s_m;
        let dt = 1.0 / f64::from(sample_rate_hz);
        // Retain the legacy 1 m/s history for TMM; the graph supplies its
        // explicit minimum and refuses slower positive jets instead of clamping.
        let max_delay = geometry.cut_up_m / (island.convection_ratio * minimum_jet_speed);
        let history_samples = (max_delay / dt).ceil().max(2.0);
        if !history_samples.is_finite() || history_samples > max_history.saturating_sub(2) as f64
            || ![zc, s_m, inertance].iter().all(|x| x.is_finite() && *x > 0.0)
            || !resistance.is_finite() || resistance < 0.0 {
            return Err(FlueError::Invalid { what: "flue geometry exceeds finite coefficient or convection-history limits" });
        }
        let history_len = history_samples as usize + 2;
        let seed_key = StreamKey {
            seed: voice_seed,
            kernel: FLUE_SEED_KERNEL,
            tile: voice_index,
        };
        let momentum_thickness_m = island.theta_over_b * 0.5 * geometry.flue_height_m;
        if !momentum_thickness_m.is_finite() || momentum_thickness_m <= 0.0 {
            return Err(FlueError::Invalid { what: "jet profile must have finite positive momentum thickness" });
        }
        Ok(FlueVoice {
            island,
            geometry,
            line,
            gas: *gas,
            mouth_contraction: None,
            mouth_radiation: None,
            window_radiation: None,
            window_mass: false,
            viscous_windway: None,
            unsteady_windway: None,
            rho,
            nu,
            zc,
            dt,
            blow_pa: 0.0,
            jet_speed_m_s: 0.0,
            history: vec![0.0; history_len],
            head: 0,
            q_prev: 0.0,
            flow_prev: 0.0,
            band_state: [0.0; 2],
            band_hp: [0.0; 2],
            momentum_thickness_m,
            displacement_state: 0.0,
            nolle: None,
            inertance,
            resistance,
            source_coeff: rho * geometry.end_correction_m() / s_m,
            seed: seed_key.stream(),
            seed_key,
            block_index: 0,
            samples: 0,
            poisoned: false,
            diags: Vec::new(),
        })
    }

    /// Accepted source sample count, independent of host block partitioning.
    #[must_use]
    pub const fn samples_rendered(&self) -> u64 { self.samples }

    /// Fixed physical sample period [s]. No retuning or resampling is implied.
    #[must_use]
    pub const fn sample_period_s(&self) -> f64 { self.dt }

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
    /// [`FlueError`] on a non-finite or negative pressure, an inadmissible
    /// jet direction, or an unknown fingering. Refusal leaves held inputs intact.
    pub fn apply(&mut self, control: FlueControl) -> Result<(), FlueError> {
        if self.poisoned { return Err(FlueError::Invalid { what: "flue voice is poisoned by a previous failed sample" }); }
        match control {
            FlueControl::SetBlowingPressure(p) => {
                let speed = self.validated_jet_speed(p)?;
                self.blow_pa = p;
                // A pressure control supplies force; it never overwrites moving gas.
                if self.unsteady_windway.is_none() { self.jet_speed_m_s = speed; }
            }
            FlueControl::SetJetAngle(angle) => {
                self.validate_jet_angle(angle)?;
                self.geometry.jet_angle_rad = angle;
            }
            FlueControl::SetBoundaryClosingForce { node, force_n } => {
                self.set_boundary_closing_force(node, force_n)?;
            }
            FlueControl::Fingering {
                index,
                fade_samples,
            } => self.line.switch_fingering(index, fade_samples)?,
        }
        Ok(())
    }

    /// Validate a direction against the original fixed-window jet geometry.
    /// This is the existing tilted-cut image, not motion of the labium, lip
    /// coverage, a changing mouth area or resolved jet turning. No stored mass,
    /// stiffness or inertance is replaced. Its changed injection enters the
    /// same source-pressure derivative and mouth-work balance when stepping.
    ///
    /// # Errors
    /// Nonfinite direction, an angle outside (-pi/2, pi/2), or a geometric
    /// offset that cannot be represented in the jet's transverse coordinate.
    pub fn validate_jet_angle(&self, angle_rad: f64) -> Result<(), FlueError> {
        let mut geometry = self.geometry;
        geometry.jet_angle_rad = angle_rad;
        geometry.validate()?;
        if !(geometry.effective_offset_m() / (0.5 * geometry.flue_height_m)).is_finite() {
            return Err(FlueError::Invalid { what: "jet angle exceeds the finite transverse geometry" });
        }
        Ok(())
    }

    /// Current prescribed direction [rad], including accepted angle controls.
    #[must_use]
    pub const fn jet_angle_rad(&self) -> f64 { self.geometry.jet_angle_rad }

    /// Accepted jet bulk speed. An unsteady windway retains its velocity field
    /// across pressure controls, including a release to zero blowing pressure.
    #[must_use]
    pub fn jet_speed_m_s(&self) -> f64 {
        self.jet_speed_m_s
    }

    /// Attach a finite, uniform rectangular windway before any accepted sample.
    /// Its actual exit width/height and retained inlet gas determine
    /// `p = 12 mu L v/h^2 + rho v^2/2`; the shared slit owner solves for bulk
    /// speed. That SAME speed sets injected volume, Reynolds number, convective
    /// delay, amplification band and the subsequent acoustic feedback.
    ///
    /// The streamwise length is supplied physical geometry, never inferred
    /// from an instrument name or the mouth cut-up. This is a quasi-steady
    /// parallel-plate resistance plus ideal exit kinetic head. It retains the
    /// existing card profile; it does not resolve entrance/profile development,
    /// windway inertia, turbulent flow or acoustic modulation of the mean jet.
    /// Use only for a declared rectangular channel (e.g. a recorder windway),
    /// not an invented lip-channel length for a transverse flute.
    ///
    /// # Errors
    /// Invalid/unrepresentable geometry, repeated/late attachment, a poisoned
    /// voice, or a held pressure outside the ACTUAL jet-speed/history domain.
    pub fn with_viscous_windway(mut self, length_m: f64) -> Result<Self, FlueError> {
        if self.poisoned || self.samples != 0 || self.viscous_windway.is_some()
            || self.unsteady_windway.is_some() {
            return Err(FlueError::Invalid { what: "viscous windway requires one cold attachment" });
        }
        self.viscous_windway = Some(crate::bernoulli_aperture::viscous::ViscousSlit::from_gas(
            &self.gas, length_m).map_err(FlueError::Physics)?);
        self.jet_speed_m_s = self.validated_jet_speed(self.blow_pa)?;
        Ok(self)
    }

    /// The actual source-gas viscosity and supplied windway length. None retains
    /// the original ideal jet relation. The slit law's quadratic pressure term
    /// represents exit kinetic head here, not a loss in the acoustic mouth loop.
    #[must_use]
    pub const fn viscous_windway(&self) -> Option<&crate::bernoulli_aperture::viscous::ViscousSlit> {
        self.viscous_windway.as_ref()
    }

    /// Replace the legacy `0.61 a_eq` mouth inertance with a computed window
    /// mass before any accepted sample. `branch_mass_per_density` [1/m] is the
    /// whole window branch from the graph's junction to the far field (e.g.
    /// `fs_bem::wall_window`); `riser_mass_per_density` is the part the graph
    /// already carries as its inlet section. The remainder times the inlet gas
    /// density becomes the series mouth inertance, so the static mass between
    /// junction and far field is the computed one. The jet's source
    /// coefficient and the authored linear loss keep their owners. A compact
    /// mouth reaction brings its own outer mass and refuses this; the
    /// window's radiation load ([`FlueVoice::with_window_radiation`]) follows
    /// it and takes its static mass out of this one.
    ///
    /// # Errors
    /// Late attachment, an attached mouth reaction, or a remainder that is
    /// not finite and positive.
    pub fn with_window_mass(mut self, branch_mass_per_density: f64, riser_mass_per_density: f64)
        -> Result<Self, FlueError>
    {
        if self.poisoned || self.samples != 0 || self.mouth_radiation.is_some()
            || self.window_radiation.is_some() {
            return Err(FlueError::Invalid { what: "window mass requires a cold voice without a mouth reaction" });
        }
        let inertance = self.rho * (branch_mass_per_density - riser_mass_per_density);
        if !(inertance.is_finite() && inertance > 0.0) {
            return Err(FlueError::Invalid { what: "window mass must exceed the inlet section's own mass" });
        }
        self.inertance = inertance;
        self.window_mass = true;
        Ok(self)
    }

    /// Series mouth inertance [Pa s^2/m^3]: `rho 0.61 a_eq / S_m`, the rest
    /// of a computed window mass (less a window radiation load's static
    /// mass), or zero under a compact mouth reaction.
    #[must_use]
    pub const fn mouth_inertance_pa_s2_m3(&self) -> f64 { self.inertance }

    /// Attach a material-gas windway with resolved unsteady no-slip momentum.
    ///
    /// The supplied length and transverse mesh determine the evolving velocity
    /// profile. The existing card's feedback is convected across its physical
    /// cut-up on a separate, explicitly refinable spatial mesh. Actual accepted
    /// bulk speed drives source flow, amplification, Reynolds validity and this
    /// transport EACH sample, independently of host callback partitioning.
    ///
    /// The windway starts at rest even when a pressure is already held. A pressure
    /// release leaves the fluid moving; wall friction and the explicit exit remove
    /// its kinetic energy. Zero-speed transport retains the physical field.
    /// This resolves fully developed transverse laminar flow only, with an ideal
    /// uniform kinetic-head terminal. The surrounding jet card, geometry and
    /// observer keep their existing scope; no resolved jet or flute claim follows.
    ///
    /// # Errors
    /// Nonpositive/nonfinite length, either mesh outside 2..=4096, repeated/late
    /// attachment, poisoned voice, or a held pressure above the admitted jet Mach.
    pub fn with_unsteady_windway(mut self, length_m: f64, transverse_cells: usize,
        convection_cells: usize) -> Result<Self, FlueError>
    {
        if self.poisoned || self.samples != 0 || self.viscous_windway.is_some()
            || self.unsteady_windway.is_some() {
            return Err(FlueError::Invalid { what: "unsteady windway requires one cold source-channel attachment" });
        }
        self.unsteady_windway = Some(windway::UnsteadyWindway::new(&self.gas, &self.geometry,
            length_m, transverse_cells, convection_cells, self.dt, self.island.convection_ratio)?);
        self.validated_jet_speed(self.blow_pa)?;
        self.jet_speed_m_s = 0.0;
        Ok(self)
    }

    /// Replace the island's receptivity chain by the published Nolle/Price jet transfer
    /// ([`nolle_jet`]) before any accepted sample. The mouth loop, labium cut, windway and wind line
    /// are unchanged; the drive source becomes `M_d dQ_j/dt` with Price's `M_d = 0.88 rho / H`, and
    /// the labium displacement comes from the transfer's own history. A slit windway implies a
    /// Poiseuille profile (centre speed 1.5 x bulk); no windway keeps the Bernoulli central speed.
    /// The unsteady windway is refused: its evolving profile has no published centre-speed law here.
    /// The authored displacement gain of the card is unused while this transfer is attached.
    ///
    /// # Errors
    /// Late or repeated attachment, a poisoned voice, an unsteady windway, or the transfer's own
    /// geometry/realisation refusal.
    pub fn with_nolle_jet(mut self, spec: nolle_jet::NolleSpec) -> Result<Self, FlueError> {
        if self.poisoned || self.samples != 0 || self.nolle.is_some()
            || self.unsteady_windway.is_some() {
            return Err(FlueError::Invalid {
                what: "nolle jet requires one cold attachment and no unsteady windway",
            });
        }
        let centre_over_bulk = if self.viscous_windway.is_some() {
            nolle_jet::POISEUILLE_CENTRE_OVER_BULK
        } else {
            1.0
        };
        // Slowest bulk jet the retained convection history admits (the voice refuses slower ones).
        let history_s = (self.history.len() - 2) as f64 * self.dt;
        let min_bulk = self.geometry.cut_up_m / (self.island.convection_ratio * history_s);
        let jet = nolle_jet::NolleJet::new(spec, centre_over_bulk, self.geometry.flue_height_m,
            self.geometry.cut_up_m, self.dt, centre_over_bulk * min_bulk)
            .map_err(|what| FlueError::Invalid { what })?;
        self.source_coeff = spec.drive_inertance_coeff * self.rho / self.geometry.flue_width_m;
        self.nolle = Some(jet);
        Ok(self)
    }

    /// The attached published jet transfer, if any.
    #[must_use]
    pub const fn nolle_jet(&self) -> Option<&nolle_jet::NolleJet> {
        self.nolle.as_ref()
    }

    /// Actual source field, material state, energy and both supplied resolutions.
    #[must_use]
    pub const fn unsteady_windway(&self) -> Option<&windway::UnsteadyWindway> {
        self.unsteady_windway.as_ref()
    }

    fn speed_at_pressure(&self, pressure_pa: f64) -> Result<f64, FlueError> {
        if let Some(channel) = &self.unsteady_windway {
            return channel.steady_speed(pressure_pa);
        }
        match self.viscous_windway {
            Some(law) => {
                let flow = law.evaluate(self.geometry.flue_width_m, self.geometry.flue_height_m,
                    self.rho, pressure_pa).map_err(FlueError::Physics)?;
                let speed = flow.flow_m3_s / (self.geometry.flue_width_m * self.geometry.flue_height_m);
                let delay = self.geometry.cut_up_m / (self.island.convection_ratio * speed);
                if !speed.is_finite() || (pressure_pa > 0.0
                    && delay > (self.history.len() - 2) as f64 * self.dt) {
                    return Err(FlueError::Invalid { what: "viscous windway jet exceeds finite speed or convection history" });
                }
                Ok(speed)
            }
            None => Ok((2.0 * pressure_pa / self.rho).sqrt()),
        }
    }

    /// Jet Reynolds number `U * 2b / nu` for the current blowing pressure.
    #[must_use]
    pub fn reynolds(&self) -> f64 {
        self.jet_speed_m_s() * self.geometry.flue_height_m / self.nu
    }

    /// Jet transit delay `W / c_p` [s] (infinite at zero blowing).
    #[must_use]
    pub fn transit_s(&self) -> f64 {
        if let Some(jet) = &self.nolle {
            return jet.transit_s(jet.centre_over_bulk() * self.jet_speed_m_s());
        }
        self.geometry.cut_up_m / (self.island.convection_ratio * self.jet_speed_m_s())
    }

    // Evaluate the original filter, displacement and cut law without advancing
    // any memory. Nonlinear trials all use the same one-per-sample onset seed.
    fn jet_trial(&self, velocity: f64, seed: f64, c: JetStepCoefficients)
        -> Result<JetTrial, FlueError>
    {
        let mut band_state = self.band_state;
        band_state[0] += c.band_alpha*(velocity-band_state[0]);
        band_state[1] += c.band_alpha*(band_state[0]-band_state[1]);
        let hp_in = band_state[1];
        let hp_out = c.hp_alpha*(self.band_hp[1]+hp_in-self.band_hp[0]);
        let displacement_m = self.displacement_state
            + self.dt*(hp_out+seed-c.leak_rate*self.displacement_state);
        let eta_m = match c.eta_override {
            Some(eta) => eta,
            None => self.island.displacement_gain*displacement_m,
        };
        let flow_m3_s = if c.q_max > 0.0 {
            0.5*c.q_max*(1.0+((eta_m-c.offset_m)/c.half_width_m).tanh())
        } else { 0.0 };
        let source_pa = self.source_coeff*(flow_m3_s-self.flow_prev)/self.dt;
        if ![velocity, band_state[0], band_state[1], hp_out, displacement_m,
            eta_m, flow_m3_s, source_pa].iter().all(|v| v.is_finite()) {
            return Err(FlueError::Invalid { what: "jet trial left the finite physical domain" });
        }
        Ok(JetTrial { band_state, band_hp: [hp_in, hp_out], displacement_m,
            eta_m, flow_m3_s, source_pa })
    }

    fn solve_jet_mouth(
        &self,
        boundary: WindBoundary,
        delayed: RetardedVelocity,
        seed: f64,
        c: JetStepCoefficients,
        nolle: Option<nolle_jet::NolleStep>,
    ) -> Result<(MouthStep, JetTrial, f64), FlueError> {
        let (radiation_direct, radiation_history) = self.mouth_radiation.as_ref()
            .map_or((0.0, 0.0), mouth_radiation::MouthRadiation::affine_pressure);
        // At most one of the two is attached; an absent one adds exact zeros.
        let (window_direct, window_history) = self.window_radiation.as_ref()
            .map_or((0.0, 0.0), window_radiation::WindowRadiation::affine_pressure);
        let (radiation_direct, radiation_history) =
            (radiation_direct + window_direct, radiation_history + window_history);
        let (l, r, zc, dt) = (self.inertance, self.resistance + radiation_direct, self.zc, self.dt);
        let area = self.geometry.mouth_area_m2();
        let jet = |q: f64| self.jet_trial(delayed.past_m_s+delayed.current_weight*q/area, seed, c);
        let mouth = |source| match &self.mouth_contraction {
            Some(loss) => mouth_loss::solve(boundary, source-radiation_history, self.q_prev, l, r, zc, dt, loss),
            None => Ok((solve_mouth_port(boundary, source-radiation_history, self.q_prev, l, r, zc, dt)?, 0.0)),
        };
        if let Some(response) = nolle {
            // Condense the exact held-input receptivity and all retained
            // convection history, then solve the cut jet and acoustic mouth
            // together. The feedback uses CURRENT Q and dQj, not a delayed
            // copy that changes the loop phase when the source clock changes.
            let resting_cut = 0.5 * c.q_max * (1.0 + (-c.offset_m / c.half_width_m).tanh());
            let direct = boundary.direct_reflection;
            let denominator = (l / dt + r) * (1.0 - direct) + zc * (1.0 + direct);
            let flow_slope = (1.0 - direct) / denominator * self.source_coeff / dt;
            // A passive contraction can only reduce dQ/dQj. Consequently the
            // largest possible d(eta)/dQj is one of these two endpoints.
            // Negative direct jet feedback is retained rather than discarded
            // in a magnitude bound. tanh' <= 1 then proves a unique cut-flow root.
            let eta_slope = response
                .eta_jet_s_m2
                .max(response.eta_jet_s_m2 + response.eta_mouth_s_m2 * flow_slope);
            let feedback = (0.5 * c.q_max / c.half_width_m) * eta_slope.max(0.0);
            if direct.abs() > 1.0
                || !flow_slope.is_finite()
                || flow_slope < 0.0
                || !feedback.is_finite()
                || feedback >= 1.0
            {
                return Err(FlueError::Invalid {
                    what: "current Nolle jet feedback is underresolved; refine the source time step",
                });
            }
            let evaluate = |cut_flow: f64| -> Result<(f64, MouthStep, JetTrial, f64), FlueError> {
                let source = self.source_coeff * (cut_flow - self.flow_prev) / dt;
                let (step, power) = mouth(source)?;
                let candidate = response
                    .candidate(step.flow_m3_s, cut_flow - resting_cut)
                    .map_err(|what| FlueError::Invalid { what })?;
                let trial = self.jet_trial(
                    0.0,
                    0.0,
                    JetStepCoefficients {
                        eta_override: Some(candidate.eta_m),
                        ..c
                    },
                )?;
                Ok((cut_flow - trial.flow_m3_s, step, trial, power))
            };
            let (mut lo, mut hi) = (0.0, c.q_max);
            if let Some(loss) = &self.mouth_contraction {
                if direct < 1.0 {
                    // Do not ask the constitutive owner to evaluate a trial
                    // beyond its Mach domain just to bracket a valid solution.
                    // Invert the same passive mouth relation at its flow bounds.
                    let mut bound = loss.effective_area_m2()
                        * (loss.gas().sound_speed * loss.spec().maximum_contracted_mach);
                    for _ in 0..32 {
                        bound = bound.next_down();
                    }
                    let cut_at = |q: f64| -> Result<f64, FlueError> {
                        let drop = loss
                            .at_flow(q)
                            .map_err(FlueError::Physics)?
                            .pressure_drop_pa;
                        let pressure = (2.0 * boundary.history_pressure_pa
                            + zc * (1.0 + direct) * q)
                            / (1.0 - direct);
                        let source = l / dt * (q - self.q_prev)
                            + r * q
                            + pressure
                            + drop
                            + radiation_history;
                        let cut = self.flow_prev + dt * source / self.source_coeff;
                        if !cut.is_finite() {
                            return Err(FlueError::Invalid {
                                what: "Nolle mouth flow domain is not representable",
                            });
                        }
                        Ok(cut)
                    };
                    lo = lo.max(cut_at(-bound)?);
                    hi = hi.min(cut_at(bound)?);
                    if lo > hi {
                        return Err(FlueError::Invalid {
                            what: "Nolle jet-mouth root lies outside the admitted flow domain",
                        });
                    }
                }
            }
            let tolerance = 32.0 * f64::EPSILON * c.q_max;
            let ends = [evaluate(lo)?, evaluate(hi)?];
            for (residual, step, trial, power) in ends {
                if residual.abs() <= tolerance {
                    return Ok((step, trial, power));
                }
            }
            if ends[0].0 > 0.0 || ends[1].0 < 0.0 {
                return Err(FlueError::Invalid {
                    what: "Nolle jet-mouth root lies outside the admitted flow domain",
                });
            }
            for _ in 0..96 {
                let flow = 0.5 * lo + 0.5 * hi;
                let (residual, step, trial, power) = evaluate(flow)?;
                if residual.abs() <= tolerance {
                    return Ok((step, trial, power));
                }
                if residual < 0.0 {
                    lo = flow;
                } else {
                    hi = flow;
                }
            }
            return Err(FlueError::Invalid {
                what: "same-sample Nolle jet-mouth flow is unresolved",
            });
        }
        // Traces with no current-sample coefficient are explicit. A closed
        // wave boundary instead fixes q independently of the current jet.
        if delayed.current_weight == 0.0 || c.q_max == 0.0 || boundary.direct_reflection == 1.0 {
            let trial = jet(if boundary.direct_reflection == 1.0 {
                -boundary.history_pressure_pa/zc
            } else { 0.0 })?;
            let (step, power) = mouth(trial.source_pa)?;
            let trial = if delayed.current_weight != 0.0 { jet(step.flow_m3_s)? } else { trial };
            return Ok((step, trial, power));
        }
        // tanh' <= 1 gives a global bound for the current-sample feedback.
        // Together with passive port and separation impedances, feedback < 1
        // makes the scalar residual strictly increasing. Otherwise refine the
        // source clock; choosing an arbitrary root would change its dynamics.
        let direct = boundary.direct_reflection;
        let denominator = (l/dt+r)*(1.0-direct)+zc*(1.0+direct);
        let compliance = (1.0-direct)/denominator;
        let eta_slope = self.island.displacement_gain*dt*c.hp_alpha*c.band_alpha*c.band_alpha
            *delayed.current_weight/area;
        let feedback = compliance*(self.source_coeff/dt)*(0.5*c.q_max)*eta_slope/c.half_width_m;
        if direct.abs() > 1.0 || !compliance.is_finite() || compliance <= 0.0
            || !feedback.is_finite() || !(0.0..1.0).contains(&feedback) {
            return Err(FlueError::Invalid { what: "current jet feedback is underresolved; refine the source time step" });
        }
        // The cut jet lies between zero and Qmax. The lossless responses at
        // those extremes, extended to zero, bracket every passive-loss root.
        let linear = |flow| solve_mouth_port(boundary,
            self.source_coeff*(flow-self.flow_prev)/dt-radiation_history, self.q_prev, l, r, zc, dt);
        let mut lo = linear(0.0)?.flow_m3_s.min(0.0);
        let mut hi = linear(c.q_max)?.flow_m3_s.max(0.0);
        if let Some(loss) = &self.mouth_contraction {
            let area = loss.effective_area_m2();
            let speed = loss.gas().sound_speed;
            let mach = loss.spec().maximum_contracted_mach;
            let mut bound = area*(speed*mach);
            // Multiplication then division can round just outside the strict
            // constitutive domain. Round the search endpoint inward, never
            // weaken the actual Mach admission or clip the accepted solution.
            for _ in 0..4 {
                if (bound/area)/speed <= mach { break; }
                bound = bound.next_down();
            }
            if !bound.is_finite() || bound <= 0.0 {
                return Err(FlueError::Invalid { what: "mouth contraction flow bound is not representable" });
            }
            lo = lo.max(-bound);
            hi = hi.min(bound);
        }
        let evaluate = |q: f64| -> Result<(f64, MouthStep, JetTrial, f64), FlueError> {
            let trial = jet(q)?;
            let (drop, power) = match &self.mouth_contraction {
                Some(loss) => {
                    let frame = loss.at_flow(q).map_err(FlueError::Physics)?;
                    (frame.pressure_drop_pa, frame.dissipated_power_w)
                }
                None => (0.0, 0.0),
            };
            let step = solve_mouth_port(boundary, trial.source_pa-drop-radiation_history,
                self.q_prev, l, r, zc, dt)?;
            let residual = q-step.flow_m3_s;
            if !residual.is_finite() {
                return Err(FlueError::Invalid { what: "jet-mouth residual is not finite" });
            }
            Ok((residual, step, trial, power))
        };
        let low = evaluate(lo)?;
        let high = evaluate(hi)?;
        let accept = |step: MouthStep, trial: JetTrial, power: f64|
            -> Result<(MouthStep, JetTrial, f64), FlueError> {
            // The wave response differs from the trial flow by the stopping
            // residual. Admit and account the ACTUALLY accepted port flow,
            // including the strict constitutive limit at a rounded endpoint.
            let power = if let Some(loss) = &self.mouth_contraction {
                loss.at_flow(step.flow_m3_s).map_err(FlueError::Physics)?.dissipated_power_w
            } else { power };
            Ok((step, trial, power))
        };
        let tolerance = 32.0*f64::EPSILON*lo.abs().max(hi.abs());
        for candidate in [low, high] {
            if candidate.0.abs() <= tolerance { return accept(candidate.1, candidate.2, candidate.3); }
        }
        if low.0 > 0.0 || high.0 < 0.0 {
            return Err(FlueError::Invalid { what: "jet-mouth root lies outside the admitted flow domain" });
        }
        for _ in 0..96 {
            let q = 0.5*lo+0.5*hi;
            let candidate = evaluate(q)?;
            if candidate.0.abs() <= tolerance { return accept(candidate.1, candidate.2, candidate.3); }
            if candidate.0 < 0.0 { lo = q; } else { hi = q; }
        }
        Err(FlueError::Invalid { what: "same-sample jet-mouth flow is unresolved" })
    }

    /// Advance one block, writing the wave launched into the pipe
    /// (`p_plus` at the mouth plane) into `out`.
    ///
    /// # Errors
    /// [`FlueError`] on an empty block or a wind line refusal.
    pub fn step_block(&mut self, out: &mut [f64]) -> Result<(), FlueError> {
        self.step_block_observed(out, |_| Ok(()))
    }

    // Observe accepted physical ports without turning every sample into a
    // separate diagnostic block. The jet's original p+ diagnostics stay intact.
    // A failed sample/observer poisons this voice; partial callbacks cannot resume.
    fn step_block_observed(
        &mut self, out: &mut [f64],
        observe: impl FnMut(&FlueResonator) -> Result<(), FlueError>,
    ) -> Result<(), FlueError> {
        self.validate_blowing_pressure(self.blow_pa)?;
        if out.is_empty() { return Err(FlueError::Invalid { what: "empty block" }); }
        let count = u64::try_from(out.len()).map_err(|_| FlueError::Invalid { what: "flue sample count overflow" })?;
        if self.samples.checked_add(count).is_none() || self.block_index.checked_add(1).is_none() {
            return Err(FlueError::Invalid { what: "flue sample or block clock overflow" });
        }
        self.diags.try_reserve(1).map_err(|_| FlueError::Invalid { what: "cannot reserve flue block diagnostics" })?;
        let result = self.advance_block(out, observe);
        if result.is_err() { self.poisoned = true; }
        result
    }

    fn advance_block(
        &mut self,
        out: &mut [f64],
        mut observe: impl FnMut(&FlueResonator) -> Result<(), FlueError>,
    ) -> Result<(), FlueError> {
        let dt = self.dt;
        let b = 0.5 * self.geometry.flue_height_m;
        let s_m = self.geometry.mouth_area_m2();
        let y0 = self.geometry.effective_offset_m();
        let (l, r) = (self.inertance, self.resistance);
        // The mouth and wave boundary are solved together below. A TMM FIR
        // can reflect instantaneously; a delayed tube graph has r0 exactly zero.
        let mut blow_work = 0.0;
        let mut source_work = 0.0;
        let mut pipe_work = 0.0;
        let mut radiation = 0.0;
        let mut joint_mouth_work = 0.0;
        let mut separation = 0.0;
        let mut stored_delta = 0.0;
        let mut numerical = 0.0;
        let mut pipe_side = 0usize;
        let mut eta_sq = 0.0;
        let mut reynolds_min = f64::INFINITY;
        let mut reynolds_max = 0.0f64;
        let mut windway_work = 0.0;
        let mut windway_stored_delta = 0.0;
        let mut windway_viscous = 0.0;
        let mut windway_exit = 0.0;
        let mut windway_numerical = 0.0;
        let mut windway_defect = 0.0;
        for slot in out.iter_mut() {
            let windway_trial = self.unsteady_windway.as_ref()
                .map(|windway| windway.prepare(self.blow_pa)).transpose()?;
            let u = if let Some(windway) = &windway_trial {
                let frame = windway.frame();
                // Spatial transport admits arbitrarily slow transient motion;
                // the actual source Mach cap remains enforced every sample.
                self.line.validate_jet_speed(frame.bulk_speed_m_s, false)?;
                windway_work += frame.pressure_work_j;
                windway_stored_delta += frame.stored_delta_j;
                windway_viscous += frame.viscous_dissipation_j;
                windway_exit += frame.exit_work_j;
                windway_numerical += frame.numerical_dissipation_j;
                windway_defect += frame.energy_residual_j;
                frame.bulk_speed_m_s
            } else { self.jet_speed_m_s() };
            let q_max = u * self.geometry.flue_height_m * self.geometry.flue_width_m;
            let transit = if u == 0.0 { f64::INFINITY } else if windway_trial.is_some() {
                self.geometry.cut_up_m / (self.island.convection_ratio * u)
            } else { self.transit_s() };
            let seed_scale = self.island.seed_amplitude_rel * u;
            let reynolds = u * self.geometry.flue_height_m / self.nu;
            reynolds_min = reynolds_min.min(reynolds);
            reynolds_max = reynolds_max.max(reynolds);
            // The existing card's growth band follows THIS sample's actual speed.
            // This retains its approximation; the channel does not mint a new card.
            let omega_c = if u > 0.0 { u / (4.0 * self.momentum_thickness_m) } else { 0.0 };
            let band_alpha = omega_c * dt / (1.0 + omega_c * dt);
            let omega_h = 0.5 * omega_c;
            let hp_alpha = 1.0 / (1.0 + omega_h * dt);
            let leak_rate = omega_c / 100.0;
            // The jet's spatial growth exp(alpha_i W) is a MAGNITUDE; its phase is the convective
            // delay W / c_p alone (Fabre 1996; Auvray 2012; the card's own stage law, phase
            // stage + 1/4 cycles, is a pure delay). The causal band filters above are only a
            // magnitude model, but they add phase: two poles at omega_c lag 2/omega_c at low
            // frequency, the high-pass followed by the receptivity integral lags 1/omega_h.
            // That is 16 theta / U of extra loop delay (about 70 percent of the transit for a
            // lip jet), which moved the first-register phase window to Str_W near 0.3 while
            // flutists play it at 0.1 to 0.25 (de la Cuadra 2008, octave jump near 0.1). Tap the
            // mouth-velocity history early by that lag so the low-frequency loop delay is the
            // convective transit. A transit already under one sample, or one that would fall
            // under one sample, keeps its sub-sample current-feedback treatment unchanged.
            let band_lag_s = if u > 0.0 { 2.0 / omega_c + 1.0 / omega_h } else { 0.0 };
            let delay_samples = if transit.is_finite() {
                let transit_samples = transit / dt;
                if transit_samples >= 1.0 {
                    (transit_samples - band_lag_s / dt).max(1.0)
                } else {
                    transit_samples
                }
            } else { 0.0 };
            // The unsteady channel supplies its affine spatial trace. The
            // quasi-steady history is sampled at THIS time minus the physical
            // transit, including a current-mouth term when transit < dt.
            let delayed = if self.nolle.is_some() {
                // The published transfer owns the jet response; the island chain idles at zero.
                RetardedVelocity { past_m_s: 0.0, current_weight: 0.0 }
            } else {
                match &windway_trial {
                    Some(windway) => windway.mouth_response(),
                    None => retarded_velocity(&self.history, self.head, delay_samples)?,
                }
            };
            let seed = if seed_scale > 0.0 {
                seed_scale * (2.0 * self.seed.next_f64() - 1.0)
            } else {
                0.0
            };
            // Prepare the published transfer without advancing its source
            // memory. Both unknown feedback flows enter the current mouth
            // solve; history is published only with the accepted wave step.
            let nolle_step = self.nolle.as_ref().map(|jet|
                jet.prepare_step(jet.centre_over_bulk()*u, self.geometry.flue_width_m, seed)
                    .map_err(|what| FlueError::Invalid { what })).transpose()?;
            let chain_seed = if nolle_step.is_some() { 0.0 } else { seed };
            let coefficients = JetStepCoefficients { q_max, half_width_m: b,
                offset_m: y0, band_alpha, hp_alpha, leak_rate, eta_override: None };
            let boundary = self.line.mouth_boundary()?;
            let (mouth, trial, separation_power) =
                self.solve_jet_mouth(boundary, delayed, chain_seed, coefficients, nolle_step)?;
            let (eta, flow, p_source) = (trial.eta_m, trial.flow_m3_s, trial.source_pa);
            let nolle_candidate = nolle_step.map(|step|
                step.candidate(mouth.flow_m3_s, flow-0.5*q_max*(1.0+(-y0/b).tanh()))
                    .map_err(|what| FlueError::Invalid { what })).transpose()?;
            eta_sq += eta * eta;
            if eta > y0 {
                pipe_side += 1;
            }
            blow_work += self.blow_pa * flow * dt;
            let next_separation = separation + separation_power*dt;
            if !next_separation.is_finite() {
                return Err(FlueError::Invalid { what: "mouth separation energy overflow" });
            }
            let (q, p_plus, p_in) = (mouth.flow_m3_s, mouth.outgoing_pa, mouth.pressure_pa);
            if ![eta, flow, p_source, q, p_plus, p_in].iter().all(|x| x.is_finite()) {
                return Err(FlueError::Invalid { what: "flue trial left the finite physical domain" });
            }
            let radiation_candidate = self.mouth_radiation.as_ref()
                .map(|load| load.candidate(q)).transpose()?;
            let window_candidate = self.window_radiation.as_ref()
                .map(|load| load.candidate(q)).transpose()?;
            // Validate the entire source field with the solved current mouth
            // velocity before either the wave or source history can publish.
            let windway_candidate = windway_trial
                .map(|windway| windway.finish(q / s_m)).transpose()?;
            let committed = self.line.push(mouth)?;
            if matches!(&self.line, FlueResonator::Tmm(_))
                && committed.to_bits() != mouth.incoming_pa.to_bits() {
                return Err(FlueError::Invalid { what: "mouth solve and committed reflection disagree" });
            }
            if let Some(windway) = windway_candidate {
                self.unsteady_windway = Some(windway);
                self.jet_speed_m_s = u;
            }
            if let Some(candidate) = nolle_candidate {
                self.nolle.as_mut().expect("prepared source").commit(candidate);
            }
            self.band_state = trial.band_state;
            self.band_hp = trial.band_hp;
            self.displacement_state = trial.displacement_m;
            self.flow_prev = flow;
            separation = next_separation;
            // Ledger terms of the implicit step (exact identity).
            source_work += p_source * q * dt;
            if let FlueResonator::Dg(duct) = &self.line {
                if let Some(frame) = duct.last_frame().filter(|frame| frame.driven_mouth.is_some()) {
                    pipe_work += frame.inlet_work_j;
                    joint_mouth_work += frame.driven_mouth.as_ref().expect("selected joint mouth").load_work_j;
                } else {
                    pipe_work += p_in * q * dt;
                }
            } else {
                pipe_work += p_in * q * dt;
            }
            radiation += r * q * q * dt;
            stored_delta += 0.5 * l * (q * q - self.q_prev * self.q_prev);
            numerical += 0.5 * l * (q - self.q_prev) * (q - self.q_prev);
            if let Some((next, frame)) = radiation_candidate {
                radiation += frame.port.dissipated_energy_j;
                stored_delta += frame.port.storage_change_j;
                self.mouth_radiation = Some(next);
            }
            if let Some((next, frame)) = window_candidate {
                radiation += frame.dissipated_energy_j;
                stored_delta += frame.storage_change_j;
                self.window_radiation = Some(next);
            }
            self.q_prev = q;
            // Mouth transverse velocity into the history ring.
            self.history[self.head] = q / s_m;
            self.head = (self.head + 1) % self.history.len();
            // The observable is the wave launched into the pipe: the mouth
            // plane is the pipe's pressure node, so the total pressure
            // there says little about the standing wave the pipe carries.
            *slot = match &self.line {
                FlueResonator::Dg(duct) => duct.last_frame()
                    .map_or(p_plus, |frame| frame.outgoing_pressure_pa),
                _ => p_plus,
            };
            self.samples += 1;
            observe(&self.line)?;
        }
        if ![windway_work, windway_stored_delta, windway_viscous, windway_exit,
            windway_numerical, windway_defect].iter().all(|v| v.is_finite()) {
            return Err(FlueError::Invalid { what: "unsteady windway block energy overflow" });
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
            jet_speed_m_s: self.jet_speed_m_s(),
            reynolds,
            reynolds_range: (reynolds_min, reynolds_max),
            inside_card_validity: reynolds_min >= lo && reynolds_max <= hi,
            transit_s: self.transit_s(),
            lock_hz: crossings as f64 / (n * dt),
            p_rms_pa: (sq / n).sqrt(),
            duty_pipe_side: pipe_side as f64 / n,
            eta_rms_over_b: (eta_sq / n).sqrt() / b,
            blow_work_j: blow_work,
            source_work_j: source_work,
            pipe_work_j: pipe_work,
            radiation_loss_j: radiation,
            joint_mouth_work_j: joint_mouth_work,
            separation_loss_j: separation,
            stored_delta_j: stored_delta,
            numerical_dissipation_j: numerical,
            ledger_defect_j: source_work - (stored_delta + numerical + radiation + pipe_work) - separation - joint_mouth_work,
            windway_work_j: windway_work,
            windway_stored_delta_j: windway_stored_delta,
            windway_viscous_loss_j: windway_viscous,
            windway_exit_work_j: windway_exit,
            windway_numerical_dissipation_j: windway_numerical,
            windway_ledger_defect_j: windway_defect,
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
        if self.poisoned { return Err(FlueError::Invalid { what: "failed flue samples cannot support a physical claim" }); }
        let outside = self
            .diags
            .iter()
            .filter(|d| !d.inside_card_validity)
            .count();
        if outside == 0 {
            return Ok(());
        }
        let seen = self.diags.iter().fold((f64::INFINITY, 0.0f64), |acc, d| {
            (acc.0.min(d.reynolds_range.0), acc.1.max(d.reynolds_range.1))
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


#[cfg(test)]
mod convection_tests {
    use super::*;
    use crate::bernoulli_aperture::{contraction::ContractionSpec,
        network::{NetworkNode, TubeNetworkSpec, TubeSection}};
    use fs_exec::CancelGate;
    use fs_material::gas::GasSpec;

    // Short windows deliberately exercise both sides of one sample of transit.
    // These are numerical fixtures, not qualified flute/recorder specimens.
    fn geometry() -> FlueGeometry {
        FlueGeometry { flue_width_m: 0.010, flue_height_m: 0.001, cut_up_m: 6e-5,
            labium_offset_m: 0.0, jet_angle_rad: 0.0,
            bore_area_m2: core::f64::consts::PI*0.006*0.006, mouth_loss_rel: 0.02 }
    }
    fn gas(pressure: f64) -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
    }

    #[test]
    fn a_jet_tilted_towards_the_pipe_lowers_the_effective_offset() {
        // Q_j grows with eta - y0, eta towards the pipe. A tilt towards the
        // pipe moves the undeflected jet W tan(angle) that way at the labium.
        let angle = 0.1_f64;
        let tilted = FlueGeometry { cut_up_m: 0.004, jet_angle_rad: angle, ..geometry() };
        assert_eq!(tilted.effective_offset_m(), -(0.004 * angle.tan()));
        let offset = FlueGeometry { labium_offset_m: 2e-4, ..tilted };
        assert_eq!(offset.effective_offset_m(), 2e-4 - 0.004 * angle.tan());
        assert!(offset.effective_offset_m() < offset.labium_offset_m);
    }
    fn voice(gas: GasState, geometry: FlueGeometry, gain: f64) -> FlueVoice {
        voice_at(gas, geometry, gain, 48_000)
    }
    fn voice_at(gas: GasState, geometry: FlueGeometry, gain: f64, rate: u32) -> FlueVoice {
        let source = network::FlueNetworkSpec {
            graph: TubeNetworkSpec {
                nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection: -0.7 }],
                sections: vec![TubeSection { nodes: [0,1], length_m: 0.3, radius_m: 0.006,
                    max_length_error_m: 0.005 }],
                sound_speed_m_s: gas.sound_speed, max_wave_memory_bytes: 1_000_000,
            },
            gases: vec![gas], viscothermal: vec![], minimum_jet_speed_m_s: 1.0,
            maximum_jet_mach: 0.3, max_jet_history_samples: 4096,
        };
        FlueVoice::from_network(Some(&fs_aeroac::jetcard::mint_tonal_interim_card().unwrap()),
            gain, 1e-6, geometry, source, rate, 3, 0, &CancelGate::new_clock_free()).unwrap()
    }
    fn contraction() -> ContractionSpec {
        ContractionSpec { area_ratio: 0.8, maximum_contracted_mach: 0.3,
            provenance: "authored numerical convection fixture".into() }
    }
    fn assert_work(v: &FlueVoice) {
        for d in v.diagnostics() {
            let scale = d.source_work_j.abs()+d.pipe_work_j.abs()+d.stored_delta_j.abs()
                +d.radiation_loss_j+d.numerical_dissipation_j+d.separation_loss_j;
            assert!(d.ledger_defect_j.abs() <= 2e-12*scale+1e-30, "{d:?}");
            assert!(d.separation_loss_j >= 0.0);
        }
    }

    #[test]
    fn g1_nolle_receptivity_uses_the_current_solved_flows_on_three_clocks() {
        use nolle_jet::{ChannelClass, NolleSpec};
        let geometry = FlueGeometry {
            cut_up_m: 0.004,
            ..geometry()
        };
        for rate in [48_000, 96_000, 192_000] {
            for pressure in [50_000.0, 200_000.0] {
                let spec = NolleSpec::price_2015(ChannelClass::Short, 0.0007).unwrap();
                let mut v = voice_at(gas(pressure), geometry, 1.0, rate)
                    .with_mouth_contraction(contraction())
                    .unwrap()
                    .with_nolle_jet(spec)
                    .unwrap();
                v.island.seed_amplitude_rel = 0.0;
                v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
                let speed = v.jet_speed_m_s();
                let travel = geometry.cut_up_m + spec.chamfer_depth_m;
                let leak = spec.leak_rel * speed / travel;
                let q_max = speed * geometry.flue_height_m * geometry.flue_width_m;
                for _ in 0..17 {
                    let previous = v.nolle_jet().unwrap().receptivity_displacement_m();
                    v.step_block(&mut [0.0]).unwrap();
                    // Exact solution of xi' + leak xi = (alpha_p Q - alpha_j
                    // dQj)/(H W), using the flows actually accepted by the
                    // acoustic solve. A previous-sample coupling fails this.
                    let source = (spec.alpha_p * v.mouth_flow_m3_s()
                        - spec.alpha_j * (v.flow_prev - 0.5 * q_max))
                        / (geometry.flue_width_m * travel);
                    let expected = previous * (-leak * v.dt).exp()
                        + source * (-(-leak * v.dt).exp_m1()) / leak;
                    let actual = v.nolle_jet().unwrap().receptivity_displacement_m();
                    assert!(
                        (actual - expected).abs()
                            < 256.0 * f64::EPSILON * actual.abs().max(expected.abs()).max(1e-30),
                        "{rate} Hz, {pressure} Pa"
                    );
                }
                assert_work(&v);
            }
        }
    }

    #[test]
    fn g3_nolle_source_partition_is_exact_and_a_failed_mouth_does_not_publish_it() {
        use nolle_jet::{ChannelClass, NolleSpec};
        let build = |maximum_contracted_mach| {
            let mut loss = contraction();
            loss.maximum_contracted_mach = maximum_contracted_mach;
            let mut v = voice(
                gas(101_325.0),
                FlueGeometry {
                    cut_up_m: 0.004,
                    ..geometry()
                },
                1.0,
            )
            .with_mouth_contraction(loss)
            .unwrap()
            .with_nolle_jet(NolleSpec::price_2015(ChannelClass::Short, 0.0).unwrap())
            .unwrap();
            v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
            v
        };
        let (mut one, mut many) = (build(0.3), build(0.3));
        let (mut a, mut b) = ([0.0; 257], [0.0; 257]);
        one.step_block(&mut a).unwrap();
        for block in b.chunks_mut(31) {
            many.step_block(block).unwrap();
        }
        assert_eq!(a, b);
        assert_eq!(
            one.nolle_jet().unwrap().receptivity_displacement_m(),
            many.nolle_jet().unwrap().receptivity_displacement_m()
        );
        assert_work(&one);
        assert_work(&many);
        let mut refused = build(1e-10);
        let before = refused.nolle_jet().unwrap().receptivity_displacement_m();
        assert!(refused.step_block(&mut [0.0]).is_err());
        assert_eq!(
            refused.nolle_jet().unwrap().receptivity_displacement_m(),
            before
        );
        assert_eq!(refused.samples, 0);
    }

    #[test]
    fn a_window_radiation_load_takes_its_mass_from_the_window_and_closes_the_ledger() {
        use fs_bem::wall_window::OuterFaceCell;
        use window_radiation::{WindowFlux, WindowRadiationLoad};
        let gas = gas(101_325.0);
        let lossless = FlueGeometry { mouth_loss_rel: 0.0, ..geometry() };
        // A 10 mm by 4 mm piston on an 8.5 mm outer cylinder: a numerical
        // fixture, M0 near 0.39 rho / sqrt(S), about 62 rho per metre.
        let cells = [OuterFaceCell { theta: 0.0, arc: 0.010 / 0.0085, axial: 0.0, length: 0.004, flux: 1.0 }];
        let load = WindowRadiationLoad::new(0.0085, &cells, WindowFlux::CutUpPiston, gas,
            1.0 / 96_000.0, &CancelGate::new_clock_free()).unwrap();
        let window = 300.0;
        assert!(voice_at(gas, lossless, 0.01, 96_000).with_window_radiation(load).is_err(),
            "the radiation load needs the computed window mass first");
        assert!(voice_at(gas, geometry(), 0.01, 96_000).with_window_mass(window, 0.0).unwrap()
            .with_window_radiation(load).is_err(), "it replaces the authored loss");
        assert!(voice(gas, lossless, 0.01).with_window_mass(window, 0.0).unwrap()
            .with_window_radiation(load).is_err(), "it shares the voice's clock");
        let v = voice_at(gas, lossless, 0.01, 96_000).with_window_mass(window, 0.0).unwrap();
        let before = v.mouth_inertance_pa_s2_m3();
        let mut v = v.with_window_radiation(load).unwrap();
        assert_eq!(v.mouth_inertance_pa_s2_m3(), before - load.load().static_inertance_pa_s2_m3());
        v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
        v.step_block(&mut [0.0; 512]).unwrap();
        assert!(v.diagnostics()[0].radiation_loss_j > 0.0);
        assert_work(&v);
        assert!(v.with_window_radiation(load).is_err(), "attachment is cold and single");
    }

    #[test]
    fn g1_retarded_affine_signal_has_no_extra_sample_at_any_ring_position() {
        for head in [0,3,31] {
            let mut history = [0.0; 32];
            for back in 0..32 { history[(head+31-back)%32] = 2.0-0.7*(back+1) as f64; }
            for delay in [0.0,0.25,0.75,1.0,1.25,5.3,31.0] {
                let value = retarded_velocity(&history, head, delay).unwrap();
                let actual = value.past_m_s+value.current_weight*2.0;
                let expected = 2.0-0.7*delay;
                assert!((actual-expected).abs() < 64.0*f64::EPSILON*expected.abs().max(1.0));
            }
            for delay in [-1.0, f64::NAN, f64::INFINITY, 32.0] {
                assert!(retarded_velocity(&history, head, delay).is_err());
            }
        }
    }

    #[test]
    fn g1_retarded_interpolation_has_its_expected_second_order_ladder() {
        // Prescribed v(t)=t^2 has midpoint interpolation error dt^2/4.
        // This grades the delay interpolation, not the full jet integrator.
        let transit = 0.03;
        let mut errors = Vec::new();
        for dt in [0.06,0.02,0.02/3.0] {
            let mut history = [0.0; 32];
            for back in 0..32 { history[31-back] = ((back+1) as f64*dt).powi(2); }
            let value = retarded_velocity(&history, 0, transit/dt).unwrap();
            let error = (value.past_m_s-transit*transit).abs(); // current v(0)=0
            assert!((error-dt*dt/4.0).abs() < 1e-16);
            errors.push(error);
        }
        for pair in errors.windows(2) { assert!((pair[0]/pair[1]-9.0).abs() < 1e-10); }
    }

    #[test]
    fn g3_gas_driven_subsample_jet_preserves_work_and_block_partition() {
        for pressure in [50_000.0,100_000.0,200_000.0] {
            for separated in [false,true] {
                let build = || {
                    let v = voice(gas(pressure), geometry(), 0.01);
                    if separated { v.with_mouth_contraction(contraction()).unwrap() } else { v }
                };
                let (mut one, mut many) = (build(), build());
                for v in [&mut one, &mut many] {
                    v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
                }
                assert_eq!(one.transit_s() < one.dt, pressure < 150_000.0);
                let (mut a, mut b) = ([0.0;257], [0.0;257]);
                one.step_block(&mut a).unwrap();
                for chunk in b.chunks_mut(31) { many.step_block(chunk).unwrap(); }
                assert_eq!(a,b);
                assert!(a.iter().any(|p| *p != 0.0));
                assert_eq!(one.network().unwrap().last_frame(), many.network().unwrap().last_frame());
                assert_eq!(one.mouth_flow_m3_s(), many.mouth_flow_m3_s());
                assert_work(&one); assert_work(&many);
            }
        }
    }

    #[test]
    fn g3_actual_viscous_windway_speed_selects_the_convection_interval() {
        let ordinary = gas(101_325.0);
        let mut spec = GasSpec::dry_air_ussa1976();
        spec.sutherland_beta *= 4.0; // Synthetic viscosity-only substitution.
        let viscous = GasState::try_new(&spec, 293.15, ordinary.pressure).unwrap();
        assert_eq!(ordinary.density, viscous.density);
        assert_eq!(ordinary.sound_speed, viscous.sound_speed);
        let mut a = voice(ordinary, geometry(), 0.01).with_viscous_windway(0.015).unwrap();
        let mut b = voice(viscous, geometry(), 0.01).with_viscous_windway(0.015).unwrap();
        for v in [&mut a, &mut b] { v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap(); }
        assert!(a.transit_s() < a.dt && b.transit_s() > b.dt);
        for v in [&a, &b] {
            assert!((v.transit_s()*v.jet_speed_m_s()*v.island.convection_ratio
                -v.geometry.cut_up_m).abs() < 1e-18);
        }
        let (mut first, mut second) = ([0.0;257], [0.0;257]);
        a.step_block(&mut first).unwrap(); b.step_block(&mut second).unwrap();
        assert_ne!(first,second);
        assert_work(&a); assert_work(&b);
    }

    #[test]
    fn g0_contraction_cap_rounding_cannot_refuse_an_interior_jet_root() {
        let mut geometry = geometry();
        geometry.flue_width_m = 0.001;
        geometry.cut_up_m = 0.00005;
        geometry.mouth_loss_rel = 0.0;
        let mut v = voice(gas(101_325.0), geometry, 0.01)
            .with_mouth_contraction(contraction()).unwrap();
        let loss = v.mouth_contraction().unwrap();
        let raw_bound = loss.effective_area_m2()
            *(loss.gas().sound_speed*loss.spec().maximum_contracted_mach);
        // Arithmetic witness: the nominal endpoint rounds above Mach 0.3.
        assert!(loss.at_flow(raw_bound).is_err());
        v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
        assert!(v.transit_s() < v.dt);
        v.step_block(&mut [0.0;257]).unwrap();
        assert!(v.mouth_flow_m3_s().abs() < raw_bound);
        v.mouth_contraction().unwrap().at_flow(v.mouth_flow_m3_s()).unwrap();
        assert!(v.diagnostics()[0].separation_loss_j > 0.0);
        assert_work(&v);
    }

    #[test]
    fn g0_underresolved_current_jet_feedback_refuses_before_any_wave_commit() {
        let mut v = voice(gas(101_325.0), geometry(), 1e6);
        v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
        assert!(matches!(v.step_block(&mut [0.0]),
            Err(FlueError::Invalid { what }) if what.contains("refine")));
        assert_eq!(v.samples_rendered(), 0);
        assert!(v.network().unwrap().last_frame().is_none());
        assert_eq!(v.mouth_flow_m3_s(), 0.0);
        assert_eq!(v.band_state, [0.0;2]);
        assert_eq!(v.displacement_state, 0.0);
        assert!(v.claim_check().is_err());
    }
}

#[cfg(test)]
mod mouth_boundary_tests {
    use super::*;
    use fs_material::gas::GasSpec;

    #[test]
    fn simultaneous_mouth_constraints_and_work_close_across_gas_states() {
        // Actual material gas states, not instrument-name or timbre presets.
        for temperature in [273.15, 293.15, 333.15] {
            for pressure in [50_000.0, 101_325.0, 180_000.0] {
                let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, pressure).unwrap();
                let geometry = FlueGeometry { flue_width_m: 0.01, flue_height_m: 0.001,
                    cut_up_m: 0.008, labium_offset_m: 0.0, jet_angle_rad: 0.0,
                    bore_area_m2: 3.0e-4, mouth_loss_rel: 0.05 };
                let l = gas.density*geometry.end_correction_m()/geometry.mouth_area_m2();
                let r = geometry.mouth_loss_rel*gas.density*gas.sound_speed/geometry.mouth_area_m2();
                let zc = gas.density*gas.sound_speed/geometry.bore_area_m2;
                let dt = 1.0/48_000.0;
                for direct in [-1.0, -0.7, 0.0, 0.2, 1.0] {
                    for history in [-50.0, 0.0, 70.0] {
                        for source in [-500.0, 0.0, 800.0] {
                            for previous in [-5.0e-5, 0.0, 1.0e-4] {
                                let boundary = WindBoundary { direct_reflection: direct, history_pressure_pa: history };
                                let step = solve_mouth_port(boundary, source, previous, l, r, zc, dt).unwrap();
                                let q = step.flow_m3_s;
                                let pressure_scale = source.abs() + step.pressure_pa.abs()
                                    + (l/dt*(q-previous)).abs() + (r*q).abs() + 1.0;
                                assert!((step.outgoing_pa-step.incoming_pa-zc*q).abs() < 1.0e-12*pressure_scale);
                                assert!((source-l/dt*(q-previous)-r*q-step.pressure_pa).abs() < 1.0e-12*pressure_scale);
                                assert_eq!(step.incoming_pa, direct*step.outgoing_pa+history);
                                let work = source*q*dt;
                                let pipe = step.pressure_pa*q*dt;
                                let radiation = r*q*q*dt;
                                let storage = 0.5*l*(q*q-previous*previous);
                                let numerical = 0.5*l*(q-previous).powi(2);
                                let scale = work.abs()+pipe.abs()+radiation+storage.abs()+numerical+1.0e-30;
                                assert!((work-pipe-radiation-storage-numerical).abs() < 2.0e-12*scale);
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn zero_direct_boundary_preserves_the_existing_delayed_graph_equations() {
        let (l, r, zc, dt) = (50.0, 70_000.0, 1.0e6, 1.0/48_000.0);
        let (source, previous, history) = (250.0, 1.0e-4, -30.0);
        let step = solve_mouth_port(WindBoundary { direct_reflection: 0.0,
            history_pressure_pa: history }, source, previous, l, r, zc, dt).unwrap();
        let flow = (source+(l/dt)*previous-2.0*history)/(l/dt+r+zc);
        assert_eq!(step.flow_m3_s.to_bits(), flow.to_bits());
        assert_eq!(step.outgoing_pa.to_bits(), (history+zc*flow).to_bits());
    }

    #[test]
    fn mouth_and_actual_tmm_commit_use_the_same_sample_through_fingering_lifts() {
        use fs_duct::{Duct, Fingering, HoleState, Segment};
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
        let radius = 0.0075;
        let zc = gas.density*gas.sound_speed/(core::f64::consts::PI*radius*radius);
        let table = FingeringTable::try_new(Duct { segments: vec![
            Segment::Cylinder { radius, length: 0.2 },
            Segment::ToneHole { hole_radius: 0.002, chimney_height: 0.004,
                bore_radius: radius, state: HoleState::Closed },
            Segment::Cylinder { radius, length: 0.2 },
        ] }, vec![
            Fingering { label: "closed".into(), holes: vec![HoleState::Closed] },
            Fingering { label: "open".into(), holes: vec![HoleState::Open] },
        ]).unwrap();
        let mut line = WindLineBank::new(&table, &gas, Termination::FlangedOpen, 48_000, zc, false).unwrap();
        let mut q = 0.0;
        for sample in 0..512 {
            if sample == 200 { line.switch_fingering(1, 32).unwrap(); }
            let boundary = line.boundary_response().unwrap();
            let step = solve_mouth_port(boundary, if sample < 8 { 100.0 } else { 0.0 },
                q, 50.0, 70_000.0, zc, 1.0/48_000.0).unwrap();
            let committed = line.push(step.outgoing_pa).unwrap();
            assert_eq!(committed.to_bits(), step.incoming_pa.to_bits());
            q = step.flow_m3_s;
        }
    }

    #[test]
    fn singular_or_nonfinite_mouth_does_not_invent_a_pressure_release() {
        let boundary = WindBoundary { direct_reflection: -1.0, history_pressure_pa: 0.0 };
        // Zero series impedance and pressure release give an undetermined flow.
        assert!(solve_mouth_port(boundary, 1.0, 0.0, 0.0, 0.0, 1.0, 0.1).is_err());
        assert!(solve_mouth_port(boundary, f64::NAN, 0.0, 1.0, 1.0, 1.0, 0.1).is_err());
        assert!(solve_mouth_port(boundary, 1.0, 0.0, -1.0, 1.0, 1.0, 0.1).is_err());
        assert!(solve_mouth_port(boundary, 1.0, 0.0, 1.0, 1.0, 0.0, 0.1).is_err());
    }
}
