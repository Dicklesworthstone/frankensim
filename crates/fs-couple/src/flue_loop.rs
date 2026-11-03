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
//!   drives a leaky integral, the jet root's displacement. Behind the
//!   band's high-pass at `omega_h = omega_c / 2` that chain is
//!   `1 / (s + omega_h)`: a ZERO-PHASE receptivity of the kind de la Cuadra
//!   (2005, sec. 4.3) measured ("phase difference between Vac and eta0 has
//!   been measured to be zero"), not Fletcher's integral. An integral's
//!   `1/omega` cancels the drive derivative's `omega` only for
//!   `omega >> omega_h`, that is for `x = omega / omega_c` well above 1;
//!   the cards run at x of about 0.4-0.55, so the derivative's quarter cycle
//!   stays in the loop (non-delay loop phase +0.27 to +0.31 cycles for x of
//!   0.1-0.55, bead pexfp.135.6); the
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
//!   with sub-sample transit coupled to the current mouth solve. This chain
//!   is the default. Its gain is authored: the cards' 90 is about 7-18x the
//!   forced jet's response at each note's speaking mode (mail 7491). It
//!   stands until the forced jet passes fg_003, fg_004, fj_001 and
//!   gas_direction_amount (ruling 7518, mail 7859, bead pexfp.135.1);
//! - opt-in by name, the forced Bickley jet's complete causal response
//!   ([`forced_jet`], ruling 7376) replaces that chain:
//!   `eta_n = sum_k g_k v_(n-k)` of the exit transverse velocity, with its
//!   delay, growth and high-Strouhal fall-off from the linearised jet and
//!   `g_0` in the same-sample mouth solve. The published Nolle/Price
//!   transfer ([`nolle_jet`]) and Auvray's jet drive ([`auvray_jet`]: in-phase
//!   receptivity, flat growth `exp(beta W/h)`, delay `W/(gamma U)`, Bickley
//!   `b = 2h/5`, Verge's source length) are the other opt-ins; all three
//!   exclude each other;
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
use crate::turbulence_noise::{BandedGaussian, FLUE_JET_NOISE_KERNEL};

use crate::bernoulli_aperture::developing::{DevelopingSlit, DevelopingSlitState};
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

/// Opt-in forced Bickley jet's causal response in place of the island's band chain.
pub mod forced_jet;

/// Opt-in jet drive of Auvray et al. (2012, 2014) in place of the island's band chain.
pub mod auvray_jet;

/// Philox kernel id of the flue onset seed stream (logical identity;
/// never a scheduling artifact).
pub const FLUE_SEED_KERNEL: u32 = 0x666c_7565; // "flue"

/// The opt-in breath noise of a flue jet (card row `jet_noise SEED GAIN ST_LO ST_HI RE_ON`): a dipole pressure jump
/// in series in the mouth, uncorrelated with the acoustic field, of rms `gain * rho U_j^2 / 2 = gain * P_f` over the band
/// `St_h = f h / U_j` in `[strouhal_lo, strouhal_hi]`, `h` the flue height, and zero while `Re_h = U_j h / nu` is below
/// `reynolds_onset`. `U_j = sqrt(2 P_f / rho)` is VERGE'S jet velocity (the maximum jet velocity of his eq. 2.12) from the
/// driving (foot) pressure `P_f`, the instantaneous blowing pressure upstream of the windway; it is NOT the card speed of
/// the loop (the windway's bulk or centre speed): the measurement defines its gain, its floor and its band on `U_j`, and a
/// windway with a viscous drop leaves `P_f` and so the noise unchanged while it slows the jet the loop sees. This is the
/// FORM of Verge's measured source (PhD thesis, TU Eindhoven 1995, sec. 5.5, Fig. 5.15a: `Delta p / (rho U_j^2 / 2)` per
/// analyser line, the pipe-end values, a recorder-like pipe with W/h 4 and h 1.00 mm, determinable only from 390 Pa); every
/// field of it is the card's. There is NO default for any field: the law has no gain that holds for another geometry
/// (Verge: noise grows with W/h, without a law). The draw is one normal per sample from the driver's own stream
/// ([`FLUE_JET_NOISE_KERNEL`]), shut or open, so the stream index never depends on the gate and the island's onset seed is
/// never touched.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JetNoiseSpec {
    /// Declared seed of the noise stream.
    pub seed: u64,
    /// Dimensionless amplitude `Delta p_rms / (rho U_j^2 / 2)` [-].
    pub gain: f64,
    /// Lower band edge in `St_h = f h / U_j` [-].
    pub strouhal_lo: f64,
    /// Upper band edge in `St_h` [-].
    pub strouhal_hi: f64,
    /// Jet Reynolds number `Re_h = U_j h / nu` [-] below which the source is shut.
    pub reynolds_onset: f64,
}

impl JetNoiseSpec {
    /// # Errors
    /// A named reason: gain finite and positive, `0 < strouhal_lo < strouhal_hi` finite, onset finite and non-negative.
    pub fn validate(&self) -> Result<(), &'static str> {
        if !(self.gain.is_finite() && self.gain > 0.0) {
            return Err("jet_noise gain must be finite and positive");
        }
        if !(self.strouhal_lo.is_finite() && self.strouhal_lo > 0.0
            && self.strouhal_hi.is_finite() && self.strouhal_hi > self.strouhal_lo) {
            return Err("jet_noise needs a Strouhal band 0 < ST_LO < ST_HI");
        }
        if !(self.reynolds_onset.is_finite() && self.reynolds_onset >= 0.0) {
            return Err("jet_noise reynolds onset must be finite and non-negative");
        }
        Ok(())
    }
}

// The attached source: the spec, its own banded stream and the last series pressure it drew (0 while shut).
struct JetNoise {
    spec: JetNoiseSpec,
    noise: BandedGaussian,
    last_pa: f64,
}

impl JetNoise {
    fn new(spec: JetNoiseSpec) -> Result<Self, FlueError> {
        spec.validate().map_err(|what| FlueError::Invalid { what })?;
        Ok(Self { noise: BandedGaussian::new(spec.seed, FLUE_JET_NOISE_KERNEL), spec, last_pa: 0.0 })
    }

    // The series pressure [Pa] of one accepted sample, drawn ONCE before any nonlinear trial; None while the gate is
    // shut (so a shut source leaves the mouth solve and the ledger on their original arithmetic). Everything follows
    // Verge's jet velocity `U_j = sqrt(2 P_f / rho)` of THIS sample's foot pressure `P_f` (never the windway's speed):
    // the amplitude is `gain * rho U_j^2 / 2 = gain * P_f`, the gate `Re_h = U_j h / nu >= reynolds_onset`, the band edges
    // `f = St U_j / h`. The stream advances by exactly one normal whether the gate is open or shut.
    fn pressure_pa(&mut self, foot_pa: f64, flue_height_m: f64, rho: f64, nu: f64, dt: f64) -> Option<f64> {
        let speed_m_s = if foot_pa > 0.0 { (2.0 * foot_pa / rho).sqrt() } else { 0.0 };
        let reynolds = speed_m_s * flue_height_m / nu;
        if !(speed_m_s > 0.0 && reynolds >= self.spec.reynolds_onset) {
            self.noise.skip();
            self.last_pa = 0.0;
            return None;
        }
        let per_hz = speed_m_s / flue_height_m;
        let w = self.noise.step(self.spec.strouhal_lo * per_hz, self.spec.strouhal_hi * per_hz, dt);
        let pressure = self.spec.gain * foot_pa * w;
        self.last_pa = pressure;
        Some(pressure)
    }
}

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

// Whether both points of a retarded lookup `delay_samples` (at least one step) back lie before the
// first nonzero history entry, `moving` entries back. A cold voice's entries before its first
// motion are exactly zero, so such a lookup reads the rest state at any delay (ruling 8132).
fn precedes_motion(delay_samples: f64, moving: u64) -> bool {
    delay_samples >= 1.0 && (delay_samples-1.0).floor() >= moving as f64
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
    // The card's speed is zero: no jet. See `jet_trial`.
    at_rest: bool,
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

/// The labium displacement one accepted sample's cut used: the jet flow was
/// `0.5 q_max (1 + tanh((eta_m - offset_m)/b))` with b half the flue height (the Bickley 2h/5
/// behind the Auvray jet), whichever voice (relay chain, forced jet or published transfer) set eta.
#[derive(Clone, Copy, Debug)]
pub(crate) struct LabiumSample {
    pub(crate) eta_m: f64,
    pub(crate) offset_m: f64,
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
    /// cut-up (the card carries no receptivity gain measurement). The cards'
    /// gain 90 is about 7-18x the forced jet's response at each note's
    /// speaking mode (mail 7491). It stands until the forced jet passes fg_003,
    /// fg_004, fj_001 and gas_direction_amount (ruling 7518, mail 7859). Unused while a jet
    /// response ([`forced_jet`] or [`nolle_jet`]) is attached.
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
    /// Samples of this block in the declared state below the minimum jet speed (ruling 8132):
    /// no deflection, and the windway's steady flow split at zero deflection.
    pub below_speed_edge_samples: usize,
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

/// The jet card's own speed at the held note (opt-in `jet_speed card-centreline`, pexfp.135.5).
/// The card's Strouhal number and band are defined on the centre speed of a near-top-hat jet
/// (`MeanJetProfile::u_centerline`). The windway's bulk speed still carries the jet's flux.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CardCentrelineSpeed {
    /// Held blowing pressure [Pa].
    pub held_pressure_pa: f64,
    /// Windway bulk speed at the held pressure [m/s].
    pub held_bulk_speed_m_s: f64,
    /// Exit centre speed fed to the card at the held pressure [m/s].
    pub held_card_speed_m_s: f64,
    /// Which centre speed the attached windway law supplies.
    pub law: &'static str,
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
    // Shah's developing-flow law and its state at the held note's pressure.
    developing_windway: Option<(DevelopingSlit, DevelopingSlitState)>,
    unsteady_windway: Option<windway::UnsteadyWindway>,
    // Opt-in: the card reads the exit centre speed instead of the bulk speed.
    card_centreline: Option<CardCentrelineSpeed>,
    rho: f64,
    nu: f64,
    zc: f64,
    dt: f64,
    blow_pa: f64,
    jet_speed_m_s: f64,
    // Mouth transverse velocity history (ring, oldest overwritten).
    history: Vec<f64>,
    head: usize,
    // History entries since the first nonzero one; every older entry is exactly zero.
    history_moving: u64,
    // Whether an accepted sample's jet has reached the declared minimum speed; before that the
    // voice is in its onset.
    reached_minimum_speed: bool,
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
    // delayed mouth velocity. Behind the high-pass above it is a zero-phase
    // receptivity, 1/(s + omega_h), not Fletcher's integral: its 1/omega
    // cancels the drive derivative's omega only for omega >> omega_h (see
    // the module doc and bead pexfp.135.6).
    displacement_state: f64,
    // Opt-in published jet transfer; its source memory commits with the mouth
    // solve, and both feedback flows belong to the same mechanical interval.
    nolle: Option<nolle_jet::NolleJet>,
    // Opt-in forced-jet kernel (see `forced_jet`): the labium displacement from the retained exit
    // velocities. None = the island's own chain.
    forced: Option<forced_jet::ForcedJetKernel>,
    // Opt-in Auvray jet drive (see `auvray_jet`): in-phase receptivity, flat growth, pure convection,
    // Bickley cut and Verge source length. None = the island's own chain.
    auvray: Option<auvray_jet::AuvrayJet>,
    // Opt-in breath noise in series in the mouth (see `JetNoiseSpec`); None = the original, byte-identical loop.
    jet_noise: Option<JetNoise>,
    // Declared minimum bulk jet speed, which sizes the retained history, and most retained
    // jet-response input samples.
    minimum_jet_speed_m_s: f64,
    max_jet_history_samples: usize,
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
        // explicit minimum. A slower positive jet never clamps its delay: see
        // `advance_block` for the rest state and the declared state below it.
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
            developing_windway: None,
            unsteady_windway: None,
            card_centreline: None,
            rho,
            nu,
            zc,
            dt,
            blow_pa: 0.0,
            jet_speed_m_s: 0.0,
            history: vec![0.0; history_len],
            head: 0,
            history_moving: 0,
            reached_minimum_speed: false,
            q_prev: 0.0,
            flow_prev: 0.0,
            band_state: [0.0; 2],
            band_hp: [0.0; 2],
            momentum_thickness_m,
            displacement_state: 0.0,
            nolle: None,
            forced: None,
            auvray: None,
            jet_noise: None,
            minimum_jet_speed_m_s: minimum_jet_speed,
            max_jet_history_samples: max_history,
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

    /// The jet card's speed `U` at the accepted state: the bulk speed, or with
    /// [`FlueVoice::with_card_centreline_speed`] the exit's centre speed, or behind a forced Bickley
    /// jet or the Nolle/Price transfer (which own the exit) that jet's centre speed `U0 = (U0/u) u` at
    /// the ratio held at attachment (pexfp.135.5). It sets the convection `c_p = (c_p/U) U`, the
    /// growth band, the onset seed, the Mach admission and the card Reynolds number. The bulk speed
    /// alone sets the jet's flux.
    #[must_use]
    pub fn card_speed_m_s(&self) -> f64 {
        if self.card_centreline.is_some() {
            if let Some(windway) = &self.unsteady_windway {
                return windway.channel().centre_speed_m_s();
            }
        }
        self.card_speed_at(self.blow_pa, self.jet_speed_m_s)
    }

    // A jet that owns the exit (the forced Bickley jet or the Nolle/Price transfer; the opt-in card
    // centreline speed refuses both) runs on its own centre speed `U0 = (U0/u) u` at the ratio held
    // at attachment: the card reads that, not the bulk speed, for the seed, the Mach admission and
    // the Reynolds number. A uniform exit's ratio is exactly 1.
    // Otherwise a uniform exit's bulk speed is its centre speed. A slit's Poiseuille profile has its
    // centre at 1.5 times the bulk speed. A developing slit keeps the inviscid core `sqrt(2P/rho)`
    // until the layers meet, and the Poiseuille centre after that.
    fn card_speed_at(&self, pressure_pa: f64, bulk_m_s: f64) -> f64 {
        if let Some(kernel) = &self.forced {
            return kernel.centre_over_bulk() * bulk_m_s;
        }
        if let Some(jet) = &self.nolle {
            return jet.centre_over_bulk() * bulk_m_s;
        }
        if self.card_centreline.is_none() {
            return bulk_m_s;
        }
        let centre = forced_jet::POISEUILLE_CENTRE_OVER_BULK * bulk_m_s;
        if self.viscous_windway.is_some() {
            centre
        } else if self.developing_windway.is_some() {
            (2.0 * pressure_pa / self.rho).sqrt().min(centre)
        } else {
            bulk_m_s
        }
    }

    /// Attach the opt-in breath noise ([`JetNoiseSpec`]) before any accepted sample. It composes with every jet
    /// model: it adds a series pressure to the mouth solve and changes nothing in the jet's own response. A voice
    /// without it never draws from the noise stream and stays bit-identical.
    ///
    /// # Errors
    /// Repeated or late attachment, a poisoned voice, or a spec [`JetNoiseSpec::validate`] refuses.
    pub fn with_jet_noise(mut self, spec: JetNoiseSpec) -> Result<Self, FlueError> {
        if self.poisoned || self.samples != 0 || self.jet_noise.is_some() {
            return Err(FlueError::Invalid { what: "jet noise requires one cold attachment" });
        }
        self.jet_noise = Some(JetNoise::new(spec)?);
        Ok(self)
    }

    /// The attached breath-noise card, if any.
    #[must_use]
    pub fn jet_noise(&self) -> Option<&JetNoiseSpec> {
        self.jet_noise.as_ref().map(|noise| &noise.spec)
    }

    /// Feed the jet card its own speed (pexfp.135.5) before any accepted sample, after the
    /// windway whose exit sets it. The card's Strouhal number and amplification band are defined
    /// on the centre speed of a near-top-hat jet (`MeanJetProfile::u_centerline`); by default the
    /// voice feeds them the windway's bulk speed. With this, the card reads the exit's centre speed
    /// (see [`FlueVoice::card_speed_m_s`]) and the bulk speed still carries the flux, so the held
    /// flow is unchanged. A slit held as fully developed is refused where its Poiseuille centre
    /// would exceed the inviscid bound `sqrt(2P/rho)`. Both this bound and the declared
    /// low-Mach limit apply to the card's speed at attachment and every pressure control.
    /// An unsteady windway supplies the actual interpolated centre velocity EACH
    /// sample, including after pressure release. Its steady original mesh supplies
    /// only cold pressure admission; the live field is never reset to that state.
    /// The jet transfers, which own their exit profiles, remain incompatible.
    ///
    /// # Errors
    /// Repeated or late attachment, a poisoned voice, a jet transfer, a
    /// non-positive held pressure, or a centre speed above the inviscid or low-Mach bound
    /// at the supplied held pressure or the currently accepted pressure.
    pub fn with_card_centreline_speed(mut self, held_pressure_pa: f64) -> Result<Self, FlueError> {
        if self.poisoned || self.samples != 0 || self.card_centreline.is_some() {
            return Err(FlueError::Invalid { what: "card centreline speed requires one cold attachment" });
        }
        if self.forced.is_some() || self.nolle.is_some() || self.auvray.is_some() {
            return Err(FlueError::Invalid { what: "card centreline speed refuses the jet transfers, which own their exit speed" });
        }
        if !(held_pressure_pa > 0.0 && held_pressure_pa.is_finite()) {
            return Err(FlueError::Invalid { what: "card centreline speed requires a positive held blowing pressure" });
        }
        let bulk = self.speed_at_pressure(held_pressure_pa)?;
        let inviscid = (2.0 * held_pressure_pa / self.rho).sqrt();
        let centre = forced_jet::POISEUILLE_CENTRE_OVER_BULK * bulk;
        let (card, law) = if let Some(windway) = &self.unsteady_windway {
            (windway.channel().steady_centre_speed_m_s(held_pressure_pa)
                .map_err(FlueError::Physics)?,
                "resolved unsteady slit: interpolated original centre velocity")
        } else if self.viscous_windway.is_some() {
            (centre, "fully developed slit: Poiseuille centre 1.5 x bulk")
        } else if self.developing_windway.is_some() {
            if inviscid < centre {
                (inviscid, "developing slit: inviscid core sqrt(2P/rho)")
            } else {
                (centre, "developing slit, developed at the exit: Poiseuille centre 1.5 x bulk")
            }
        } else {
            (bulk, "uniform exit: the bulk speed sqrt(2P/rho)")
        };
        self.card_centreline = Some(CardCentrelineSpeed { held_pressure_pa,
            held_bulk_speed_m_s: bulk, held_card_speed_m_s: card, law });
        // The receipt's held pressure can differ from the already accepted input. Neither
        // may bypass the same source-speed admission used by schedules and live controls.
        self.validated_jet_speed(held_pressure_pa)?;
        self.validated_jet_speed(self.blow_pa)?;
        Ok(self)
    }

    /// The card's centre speed at the held note, when attached.
    #[must_use]
    pub const fn card_centreline_speed(&self) -> Option<&CardCentrelineSpeed> {
        self.card_centreline.as_ref()
    }

    /// Declared minimum bulk jet speed U_min [m/s] (ruling 8132). It sizes the retained history.
    /// Below it, a jet whose lookup would read motion older than that history carries no
    /// deflection, and the windway's steady flow splits at zero deflection. It is a validity
    /// edge, not a stability limit.
    #[must_use]
    pub const fn minimum_jet_speed_m_s(&self) -> f64 { self.minimum_jet_speed_m_s }

    /// Attach a finite, uniform rectangular windway before any accepted sample.
    /// Its actual exit width/height and retained inlet gas determine
    /// `p = 12 mu L v/h^2 + rho v^2/2`; the shared slit owner solves for bulk
    /// speed. That SAME speed sets injected volume, Reynolds number, convective
    /// delay, amplification band and the subsequent acoustic feedback, unless
    /// [`FlueVoice::with_card_centreline_speed`] feeds the card the Poiseuille centre speed.
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
    /// voice, or a held pressure outside the ACTUAL finite low-Mach jet-speed domain.
    pub fn with_viscous_windway(mut self, length_m: f64) -> Result<Self, FlueError> {
        if self.poisoned || self.samples != 0 || self.viscous_windway.is_some()
            || self.developing_windway.is_some() || self.unsteady_windway.is_some() {
            return Err(FlueError::Invalid { what: "viscous windway requires one cold attachment" });
        }
        if self.auvray.is_some() {
            return Err(FlueError::Invalid { what: auvray_jet::WINDWAY_REFUSAL });
        }
        if self.forced.is_some() {
            return Err(FlueError::Invalid { what: "viscous windway must be attached before the forced jet, whose exit it sets" });
        }
        if self.card_centreline.is_some() {
            return Err(FlueError::Invalid { what: "windway must be attached before the card centreline speed, whose exit it sets" });
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

    /// Attach a straight slit windway held by Shah's laminar developing-flow law
    /// ([`crate::bernoulli_aperture::developing`]) before any accepted sample, in
    /// place of the fully developed [`FlueVoice::with_viscous_windway`]. Its exit
    /// width/height and the inlet gas hold the bulk speed at every pressure, and that
    /// speed enters the jet as the slit law's does. The held note's pressure must lie
    /// inside the law's validity (x+, Chen's low-Re term, laminarity); the attack
    /// and release pass through the same quasi-steady law and are not covered by it.
    /// The forced jet or the Nolle/Price transfer, attached after this, takes this law's exit
    /// ([`forced_jet::JetExit::developing`]: its core speed over its bulk speed at the held
    /// pressure); otherwise its core speed reaches only the island's card, through
    /// [`FlueVoice::with_card_centreline_speed`] (pexfp.135.5). A jet transfer attached before
    /// the windway is refused: the windway sets its exit.
    ///
    /// # Errors
    /// Invalid geometry or gas, repeated/late attachment, another windway, an attached
    /// jet transfer, a poisoned voice, or a held pressure outside the law's validity or
    /// the admitted jet-speed domain.
    pub fn with_developing_windway(mut self, length_m: f64, held_pressure_pa: f64)
        -> Result<Self, FlueError>
    {
        if self.poisoned || self.samples != 0 || self.viscous_windway.is_some()
            || self.developing_windway.is_some() || self.unsteady_windway.is_some() {
            return Err(FlueError::Invalid { what: "developing windway requires one cold attachment and no other windway" });
        }
        if self.auvray.is_some() {
            return Err(FlueError::Invalid { what: auvray_jet::WINDWAY_REFUSAL });
        }
        if self.nolle.is_some() {
            return Err(FlueError::Invalid { what: "developing windway must be attached before the Nolle/Price transfer, whose exit it sets" });
        }
        if self.forced.is_some() {
            return Err(FlueError::Invalid { what: "developing windway must be attached before the forced jet, whose exit it sets" });
        }
        if self.card_centreline.is_some() {
            return Err(FlueError::Invalid { what: "windway must be attached before the card centreline speed, whose exit it sets" });
        }
        if !(held_pressure_pa > 0.0 && held_pressure_pa.is_finite()) {
            return Err(FlueError::Invalid { what: "developing windway requires a positive held blowing pressure" });
        }
        let law = DevelopingSlit::from_gas(&self.gas, length_m).map_err(FlueError::Physics)?;
        let held = law.held(self.geometry.flue_width_m, self.geometry.flue_height_m, held_pressure_pa)
            .map_err(FlueError::Physics)?;
        held.validity().map_err(|what| FlueError::Invalid { what })?;
        self.developing_windway = Some((law, held));
        self.jet_speed_m_s = self.validated_jet_speed(self.blow_pa)?;
        Ok(self)
    }

    /// The developing-flow law and its validity state at the held pressure.
    #[must_use]
    pub fn developing_windway(&self) -> Option<(&DevelopingSlit, &DevelopingSlitState)> {
        self.developing_windway.as_ref().map(|(law, held)| (law, held))
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
    /// bulk speed drives source flow EACH sample. Amplification, Reynolds
    /// validity and convection use that speed by default, or the original
    /// profile's centre velocity with `with_card_centreline_speed`. Both images
    /// are independent of host callback partitioning.
    ///
    /// The windway starts at rest even when a pressure is already held. A pressure
    /// release leaves the fluid moving; wall friction and the explicit exit remove
    /// its kinetic energy. Zero-speed transport retains the physical field.
    /// The labium feedback is read at the frozen path's own delay, the transit less
    /// the band filters' low-frequency lag `2/omega_c + 1/omega_h = 16 theta / U`
    /// (with its one-sample floor and sub-sample rule), as a distance along the
    /// field at this sample's speed: `16 theta c_p / U` upstream of the cut-up,
    /// independent of `U`. On both paths the low-frequency loop delay is then the
    /// convective transit `W / c_p`.
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
            || self.developing_windway.is_some() || self.unsteady_windway.is_some() {
            return Err(FlueError::Invalid { what: "unsteady windway requires one cold source-channel attachment" });
        }
        if self.auvray.is_some() {
            return Err(FlueError::Invalid { what: auvray_jet::WINDWAY_REFUSAL });
        }
        if self.forced.is_some() {
            return Err(FlueError::Invalid { what: "the forced jet refuses the unsteady windway" });
        }
        if self.card_centreline.is_some() {
            return Err(FlueError::Invalid { what: "windway must be attached before the card centreline speed, whose exit it sets" });
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
    /// Poiseuille profile (centre speed 1.5 x bulk); no windway keeps the Bernoulli central speed;
    /// a developing windway ([`FlueVoice::with_developing_windway`], attached first) gives its own
    /// exit at the held pressure, the core speed over the bulk speed of Shah's law
    /// ([`forced_jet::JetExit::developing`]), which is exactly the Poiseuille 1.5 once the layers
    /// have met. The developing windway refuses a held pressure outside its validity (x+, Chen's
    /// low-Re term, laminarity, and a developed exit above Re_Dh 1760) when it is attached, so
    /// there is no fallback to the constant 1.5 here, which would exceed the Bernoulli bound
    /// `sqrt(2P/rho)` where the layers have not met (pexfp.135.5).
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
        if self.forced.is_some() {
            return Err(FlueError::Invalid { what: "nolle jet and forced jet exclude each other" });
        }
        if self.auvray.is_some() {
            return Err(FlueError::Invalid { what: "nolle jet and Auvray jet exclude each other" });
        }
        if self.card_centreline.is_some() {
            return Err(FlueError::Invalid { what: "card centreline speed refuses the jet transfers, which own their exit speed" });
        }
        let centre_over_bulk = if let Some((_, held)) = &self.developing_windway {
            forced_jet::JetExit::developing(held, self.geometry.flue_height_m).centre_over_bulk
        } else if self.viscous_windway.is_some() {
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

    /// Replace the island's receptivity chain (band IIR, high-pass, leak and the authored
    /// displacement gain) by the forced Bickley jet's causal response ([`forced_jet`]) before any
    /// accepted sample. The labium displacement becomes `sum_k g_k v_(n-k)` of the exit transverse
    /// velocity (mouth velocity plus onset seed). The tanh cut, source, mouth loop, windway and
    /// wind line are unchanged. A slit windway, attached first, gives the Poiseuille exit; a
    /// developing windway gives its exit at the reference pressure
    /// ([`forced_jet::JetExit::developing`]); no windway gives a uniform lip jet.
    /// `reference_pressure_pa` (the held blowing pressure) sets the Reynolds number of the
    /// integration. The slowest admitted jet and the history budget bound the retained exit
    /// velocities. The unsteady windway is refused: its jet starts from rest and decays toward
    /// zero speed after a release, where this frozen-time response would need unbounded history.
    ///
    /// # Errors
    /// Late or repeated attachment, a poisoned voice, an attached Nolle/Price transfer or
    /// unsteady windway, a non-positive reference pressure, one outside the developing
    /// windway's validity or one whose jet is below the declared minimum speed, a history
    /// longer than the budget at that minimum, or the integrator's own refusal.
    pub fn with_forced_jet(self, spec: forced_jet::ForcedJetSpec, reference_pressure_pa: f64)
        -> Result<Self, FlueError>
    {
        self.with_forced_jet_exit(spec, forced_jet::JetExitLaw::Kutta, reference_pressure_pa)
    }

    /// [`Self::with_forced_jet`] under the named exit law ([`forced_jet::JetExitLaw`]).
    ///
    /// # Errors
    /// As [`Self::with_forced_jet`]. `exit measured` is admitted only where Re_h = U0 h / nu lies
    /// within a factor 2 of 200, the one Re_h at which de la Cuadra measured it (Fig. 2.28);
    /// elsewhere the voice refuses.
    pub fn with_forced_jet_exit(mut self, spec: forced_jet::ForcedJetSpec,
        exit_law: forced_jet::JetExitLaw, reference_pressure_pa: f64) -> Result<Self, FlueError>
    {
        if self.poisoned || self.samples != 0 || self.nolle.is_some() || self.forced.is_some()
            || self.auvray.is_some() || self.unsteady_windway.is_some() {
            return Err(FlueError::Invalid {
                what: "forced jet requires one cold attachment and no other jet transfer or unsteady windway",
            });
        }
        if self.card_centreline.is_some() {
            return Err(FlueError::Invalid { what: "card centreline speed refuses the jet transfers, which own their exit speed" });
        }
        if !(reference_pressure_pa > 0.0 && reference_pressure_pa.is_finite()) {
            return Err(FlueError::Invalid { what: "forced jet reference pressure must be positive" });
        }
        let reference_bulk = self.validated_jet_speed(reference_pressure_pa)?;
        if reference_bulk < self.minimum_jet_speed_m_s {
            return Err(FlueError::Invalid { what: "forced jet reference pressure gives a jet below the declared minimum speed" });
        }
        let exit = if let Some((law, _)) = &self.developing_windway {
            // The response is integrated at the reference pressure, so the exit is held there.
            let held = law.held(self.geometry.flue_width_m, self.geometry.flue_height_m,
                reference_pressure_pa).map_err(FlueError::Physics)?;
            held.validity().map_err(|what| FlueError::Invalid { what })?;
            forced_jet::JetExit::developing(&held, self.geometry.flue_height_m)
        } else if self.viscous_windway.is_some() {
            forced_jet::JetExit::POISEUILLE
        } else {
            forced_jet::JetExit::UNIFORM
        };
        if exit_law == forced_jet::JetExitLaw::Measured {
            let re_h = exit.centre_over_bulk * reference_bulk * self.geometry.flue_height_m / self.nu;
            if !forced_jet::measured_exit_reynolds_admitted(re_h) {
                return Err(FlueError::Invalid { what: "exit measured: U_j h / nu is outside 100 to 400, a numerical neighbourhood (factor 2) chosen in code, not a measured range, around the only Re_h (200) at which the exit receptivity was measured (de la Cuadra 2005 Fig. 2.28)" });
            }
        }
        let kernel = forced_jet::ForcedJetKernel::new_with_exit(spec, exit_law, exit,
            self.geometry.flue_height_m, self.geometry.cut_up_m, self.nu, reference_bulk,
            self.minimum_jet_speed_m_s, self.max_jet_history_samples, self.dt)
            .map_err(|what| FlueError::Invalid { what })?;
        self.island.provenance.push_str(&format!(
            "; jet response: forced Bickley jet (b {:.4e} m, U0/U {}, station {:.3} b, Re {:.1} at {} Pa, spacing {} b, window {} b/U0, held {:.4e} s at that pressure) in place of the band chain and displacement gain",
            kernel.half_width_m(), kernel.centre_over_bulk(), kernel.step_response().station(),
            kernel.reynolds(), reference_pressure_pa, spec.spacing, spec.window,
            kernel.held_gain_s(reference_bulk)));
        if let fs_aeroac::forced_jet::ExitReceptivity::Measured { displacement_per_velocity } =
            kernel.step_response().exit()
        {
            self.island.provenance.push_str(&format!(
                "; exit receptivity: measured, eta(0) = (h/U0) v in phase, R {displacement_per_velocity:.6} b per v/U0, no carried term (Terrien et al. arXiv 1402.3143 Eq. 1 after de la Cuadra 2005 sec. 2.3.5)"));
        }
        self.forced = Some(kernel);
        Ok(self)
    }

    /// The attached forced-jet kernel, if any.
    #[must_use]
    pub const fn forced_jet(&self) -> Option<&forced_jet::ForcedJetKernel> {
        self.forced.as_ref()
    }

    /// Replace the island's receptivity chain (band IIR, high-pass, leak and the authored
    /// displacement gain) by the jet drive of Auvray et al. ([`auvray_jet`]) before any accepted
    /// sample. The labium displacement becomes `h exp(beta W/h) (v + seed) / U` of the mouth
    /// velocity retarded by `W / (gamma U)` alone (no band lag); the cut's half-width becomes the
    /// Bickley `2h/5` with the swing `2 b H U`; the source coefficient becomes `rho delta_d / S_m`
    /// with Verge's `delta_d = (4/pi) sqrt(2 h W)` in place of the mouth end correction. The mouth
    /// inertance and loss, any card `mouth_loss` and the wind line are unchanged. `U = sqrt(2P/rho)`,
    /// so every windway law and the card centreline speed are refused. The onset seed enters as
    /// exit velocity at the current sample, undelayed: it is an input from outside the loop, and a
    /// delay does not change its statistics. While the voice is cold, the retained history grows to
    /// cover `W / (gamma U_min)`.
    ///
    /// # Errors
    /// Late or repeated attachment, a poisoned voice, another jet transfer, any windway law, the
    /// card centreline speed, the kernel's own refusal, or a history beyond the budget at `U_min`.
    pub fn with_auvray_jet(mut self, spec: auvray_jet::AuvraySpec) -> Result<Self, FlueError> {
        if self.poisoned || self.samples != 0 || self.auvray.is_some() {
            return Err(FlueError::Invalid { what: "Auvray jet requires one cold attachment" });
        }
        if self.nolle.is_some() || self.forced.is_some() {
            return Err(FlueError::Invalid { what: "the Auvray jet excludes the other jet transfers" });
        }
        if self.card_centreline.is_some() {
            return Err(FlueError::Invalid { what: "card centreline speed refuses the jet transfers, which own their exit speed" });
        }
        if self.viscous_windway.is_some() || self.developing_windway.is_some()
            || self.unsteady_windway.is_some() {
            return Err(FlueError::Invalid { what: auvray_jet::WINDWAY_REFUSAL });
        }
        let jet = auvray_jet::AuvrayJet::new(spec, self.geometry.flue_height_m, self.geometry.cut_up_m)
            .map_err(|what| FlueError::Invalid { what })?;
        let history_samples = (jet.transit_s(self.minimum_jet_speed_m_s) / self.dt).ceil().max(2.0);
        if !history_samples.is_finite()
            || history_samples > self.max_jet_history_samples.saturating_sub(2) as f64 {
            return Err(FlueError::Invalid {
                what: "Auvray jet convection at the declared minimum speed exceeds the retained history budget",
            });
        }
        let needed = history_samples as usize + 2;
        if needed > self.history.len() {
            // Cold: every entry is the rest state, so a longer ring changes no lookup.
            self.history = vec![0.0; needed];
            self.head = 0;
        }
        let source_coeff = self.rho * jet.source_length_m() / self.geometry.mouth_area_m2();
        if !(source_coeff.is_finite() && source_coeff > 0.0) {
            return Err(FlueError::Invalid { what: "Auvray jet source coefficient is not finite" });
        }
        self.island.provenance.push_str(&format!(
            "; jet transfer: Auvray et al. 2012/2014 jet drive (beta {}, gamma {}, source {}, at W/h {:.3}) in place of the band chain and displacement gain: in-phase receptivity h exp(beta W/h) / U (growth {:.4}), delay W / (gamma U), Bickley b = 2h/5 = {:.4e} m with swing 2bHU, Verge delta_d = (4/pi) sqrt(2hW) = {:.4e} m, U = sqrt(2P/rho)",
            jet.spec().beta, jet.spec().gamma, jet.spec().source, jet.cut_up_over_height(),
            jet.growth(), jet.half_width_m(), jet.source_length_m()));
        self.source_coeff = source_coeff;
        self.auvray = Some(jet);
        Ok(self)
    }

    /// The attached Auvray jet drive, if any.
    #[must_use]
    pub const fn auvray_jet(&self) -> Option<&auvray_jet::AuvrayJet> {
        self.auvray.as_ref()
    }

    /// Half-width `b` of the labium cut `Q_j = (q_max / 2)(1 + tanh((eta - y0) / b))`: half the
    /// flue height, or the Bickley `2h/5` behind the Auvray jet.
    #[must_use]
    pub fn cut_half_width_m(&self) -> f64 {
        match &self.auvray {
            Some(jet) => jet.half_width_m(),
            None => 0.5 * self.geometry.flue_height_m,
        }
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
        if let Some((law, _)) = &self.developing_windway {
            // Validity binds the held note at attachment; other pressures use the same law.
            return Ok(law.held(self.geometry.flue_width_m, self.geometry.flue_height_m, pressure_pa)
                .map_err(FlueError::Physics)?.bulk_speed_m_s);
        }
        match self.viscous_windway {
            Some(law) => {
                let flow = law.evaluate(self.geometry.flue_width_m, self.geometry.flue_height_m,
                    self.rho, pressure_pa).map_err(FlueError::Physics)?;
                let speed = flow.flow_m3_s / (self.geometry.flue_width_m * self.geometry.flue_height_m);
                // A jet whose transit outruns the retained history is admitted: `advance_block`
                // reads the rest state or takes the declared state below the minimum speed.
                if !speed.is_finite() {
                    return Err(FlueError::Invalid { what: "viscous windway jet speed is not finite" });
                }
                Ok(speed)
            }
            None => Ok((2.0 * pressure_pa / self.rho).sqrt()),
        }
    }

    /// Jet Reynolds number `U * 2b / nu` on the card's speed for the current blowing pressure.
    #[must_use]
    pub fn reynolds(&self) -> f64 {
        self.card_speed_m_s() * self.geometry.flue_height_m / self.nu
    }

    /// Jet transit delay `W / c_p` [s] (infinite at zero blowing).
    #[must_use]
    pub fn transit_s(&self) -> f64 {
        if let Some(jet) = &self.nolle {
            return jet.transit_s(jet.centre_over_bulk() * self.jet_speed_m_s());
        }
        if let Some(jet) = &self.auvray {
            return jet.transit_s(self.jet_speed_m_s());
        }
        self.geometry.cut_up_m / (self.island.convection_ratio * self.card_speed_m_s())
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
        // A jet at rest has no band, high-pass or leak corner: all three scale with U / theta and
        // are zero at U = 0, where the high-pass would hold its last output and the integrator,
        // without leak, would ramp on it for as long as the rest lasts (24 m, 0.29 s after the
        // release of the R6 D4 flute render). With no jet there is no deflection and no
        // receptivity memory: the chain takes the rest state the onset starts from, exactly zero,
        // and the next blow starts from it. Before the first blow the chain is already exactly
        // zero there, and a jet at rest carries no flow, so this changes no sample until a jet
        // blows again after a rest.
        let (band_state, hp_in, hp_out, displacement_m) = if c.at_rest {
            ([0.0; 2], 0.0, 0.0, 0.0)
        } else {
            (band_state, hp_in, hp_out, displacement_m)
        };
        let eta_m = match (c.eta_override, &self.forced) {
            (Some(eta), _) => eta,
            // The forced jet's prepared known part plus g_0 times this trial's exit velocity.
            (None, Some(kernel)) => kernel.current().0 + kernel.current().1*velocity,
            // The Auvray drive: in-phase receptivity of this trial's retarded velocity plus the
            // onset seed as exit velocity. Its U is the accepted jet speed: it refuses the
            // unsteady windway, the only owner of a different per-sample speed.
            (None, None) => match &self.auvray {
                Some(jet) => jet.receptivity_s(self.jet_speed_m_s)*(velocity+seed),
                None => self.island.displacement_gain*displacement_m,
            },
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
        series_noise_pa: Option<f64>,
    ) -> Result<(MouthStep, JetTrial, f64), FlueError> {
        let (radiation_direct, radiation_history) = self.mouth_radiation.as_ref()
            .map_or((0.0, 0.0), mouth_radiation::MouthRadiation::affine_pressure);
        // At most one of the two is attached; an absent one adds exact zeros.
        let (window_direct, window_history) = self.window_radiation.as_ref()
            .map_or((0.0, 0.0), window_radiation::WindowRadiation::affine_pressure);
        let (radiation_direct, radiation_history) =
            (radiation_direct + window_direct, radiation_history + window_history);
        // The opt-in breath noise is a pressure in series with the jet's source, the SAME value in every nonlinear trial
        // of this sample. Every use below subtracts `radiation_history` from the jet's source pressure, so a series
        // pressure `p` enters as `radiation_history - p` and the four mouth solves (affine, Nolle, nonlinear, linear
        // bracket) see it identically. Absent noise leaves `radiation_history` untouched, bit for bit.
        let radiation_history = match series_noise_pa {
            Some(p) => radiation_history - p,
            None => radiation_history,
        };
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
                let residual = cut_flow - trial.flow_m3_s;
                // Accept the root's cut flow and the source used by this mouth
                // solve as one state. The constitutive cut is the residual
                // oracle; substituting it after convergence would change the
                // source without re-solving the mouth and break its work identity.
                let accepted = JetTrial {
                    flow_m3_s: cut_flow,
                    source_pa: source,
                    ..trial
                };
                Ok((residual, step, accepted, power))
            };
            let (mut lo, mut hi): (f64, f64) = (0.0, c.q_max);
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
        // A negative current feedback (a signed precursor in the forced kernel,
        // e.g. the measured exit's g_0 < 0) gives the residual a slope >= 1.
        // The bound that matters is feedback < 1.
        let direct = boundary.direct_reflection;
        let denominator = (l/dt+r)*(1.0-direct)+zc*(1.0+direct);
        let compliance = (1.0-direct)/denominator;
        let eta_slope = match (&self.forced, &self.auvray) {
            (Some(kernel), _) => kernel.current().1*delayed.current_weight/area,
            (None, Some(jet)) => jet.receptivity_s(self.jet_speed_m_s)*delayed.current_weight/area,
            (None, None) => self.island.displacement_gain*dt*c.hp_alpha*c.band_alpha*c.band_alpha
                *delayed.current_weight/area,
        };
        let feedback = compliance*(self.source_coeff/dt)*(0.5*c.q_max)*eta_slope/c.half_width_m;
        if direct.abs() > 1.0 || !compliance.is_finite() || compliance <= 0.0
            || !feedback.is_finite() || !(feedback < 1.0) {
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
        mut observe: impl FnMut(&FlueResonator) -> Result<(), FlueError>,
    ) -> Result<(), FlueError> {
        self.step_block_labium_observed(out, |line, _| observe(line))
    }

    // The same, also handing the observer the displacement each accepted cut used.
    fn step_block_labium_observed(
        &mut self, out: &mut [f64],
        observe: impl FnMut(&FlueResonator, LabiumSample) -> Result<(), FlueError>,
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
        mut observe: impl FnMut(&FlueResonator, LabiumSample) -> Result<(), FlueError>,
    ) -> Result<(), FlueError> {
        let dt = self.dt;
        let b = self.cut_half_width_m();
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
        let mut below_edge_samples = 0usize;
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
            let mut windway_trial = self.unsteady_windway.as_ref()
                .map(|windway| windway.prepare(self.blow_pa, self.card_centreline.is_some())).transpose()?;
            let u = if let Some(windway) = &windway_trial {
                let frame = windway.frame();
                // Spatial transport admits arbitrarily slow transient motion;
                // the actual source Mach cap remains enforced every sample.
                self.line.validate_jet_speed(frame.bulk_speed_m_s)?;
                windway_work += frame.pressure_work_j;
                windway_stored_delta += frame.stored_delta_j;
                windway_viscous += frame.viscous_dissipation_j;
                windway_exit += frame.exit_work_j;
                windway_numerical += frame.numerical_dissipation_j;
                windway_defect += frame.energy_residual_j;
                frame.bulk_speed_m_s
            } else { self.jet_speed_m_s() };
            // The bulk speed carries the flux; the card's own speed (the bulk speed unless the
            // card centreline speed is attached) sets its convection, band, seed and Reynolds number.
            let card_u = match &windway_trial {
                Some(windway) if self.card_centreline.is_some() => windway.centre_speed_m_s(),
                _ => self.card_speed_at(self.blow_pa, u),
            };
            // The same candidate centre speed controls transport, instability
            // and Mach admission before any source or bore state is published.
            self.line.validate_jet_speed(card_u)?;
            // The Auvray drive's Bickley jet swings Q_in = H b U (1 + tanh) over 2 b H U.
            let q_max = if self.auvray.is_some() {
                2.0 * b * self.geometry.flue_width_m * u
            } else {
                u * self.geometry.flue_height_m * self.geometry.flue_width_m
            };
            let transit = if card_u == 0.0 { f64::INFINITY } else if windway_trial.is_some() {
                self.geometry.cut_up_m / (self.island.convection_ratio * card_u)
            } else { self.transit_s() };
            let seed_scale = self.island.seed_amplitude_rel * card_u;
            let reynolds = card_u * self.geometry.flue_height_m / self.nu;
            reynolds_min = reynolds_min.min(reynolds);
            reynolds_max = reynolds_max.max(reynolds);
            // The existing card's growth band follows THIS sample's actual speed.
            // This retains its approximation; the channel does not mint a new card.
            let omega_c = if card_u > 0.0 { card_u / (4.0 * self.momentum_thickness_m) } else { 0.0 };
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
            // The Auvray drive has no band filter: its loop delay is the convection W / (gamma U) alone.
            let band_lag_s = if u > 0.0 && self.auvray.is_none() { 2.0 / omega_c + 1.0 / omega_h } else { 0.0 };
            let delay_samples = if transit.is_finite() {
                let transit_samples = transit / dt;
                if transit_samples >= 1.0 {
                    (transit_samples - band_lag_s / dt).max(1.0)
                } else {
                    transit_samples
                }
            } else if u > 0.0 { f64::INFINITY } else { 0.0 };
            // The unsteady channel supplies its affine spatial trace, read at the
            // same delay. The quasi-steady history is sampled at THIS time minus the
            // physical transit, including a current-mouth term when transit < dt.
            // Ruling 8132. The retained history covers every lookup at or above the declared
            // minimum speed U_min. In the onset, before the jet first reaches U_min, a lookup
            // beyond it that precedes the voice's first motion reads the rest state, exactly zero
            // at any delay: unbounded history at zero memory. Any other lookup beyond it below
            // U_min (a release, or a dip mid-note) takes the declared state: no deflection, so
            // the windway's steady flow splits at zero deflection. U_min is a declared validity
            // edge, not a stability limit.
            let below_edge = u > 0.0 && u < self.minimum_jet_speed_m_s;
            let onset = !self.reached_minimum_speed;
            let mut declared = false;
            let rest = RetardedVelocity { past_m_s: 0.0, current_weight: 0.0 };
            let delayed = if self.nolle.is_some() {
                // The published transfer owns the jet response; the island chain idles at zero.
                rest
            } else if self.forced.is_some() {
                // The forced jet's own history is the delay: the trial velocity is this sample's
                // exit velocity, the seed (added below) plus the unknown mouth velocity.
                RetardedVelocity { past_m_s: 0.0, current_weight: 1.0 }
            } else {
                match &mut windway_trial {
                    // At this sample's speed a characteristic covers delay/transit of the
                    // cut-up: W - 16 theta c_p/U in the ordinary case. A sub-sample transit,
                    // and a jet at rest, read at the labium as before.
                    Some(windway) => windway.mouth_response(
                        if transit.is_finite() && transit / dt >= 1.0 {
                            (delay_samples * dt / transit).min(1.0)
                        } else { 1.0 })?,
                    None if delay_samples > (self.history.len() - 1) as f64 && onset
                        && precedes_motion(delay_samples, self.history_moving) => rest,
                    None if delay_samples > (self.history.len() - 1) as f64 && below_edge => {
                        declared = true;
                        rest
                    }
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
            // The forced jet takes the seed as exit velocity, like the mouth's, and prepares this
            // sample's known displacement from its retained intervals: exact when they reach its
            // window or the rest state before its first motion, otherwise the declared state
            // below the minimum speed.
            let (delayed, chain_seed) = match &mut self.forced {
                Some(kernel) => {
                    if kernel.prepare(u).map_err(|what| FlueError::Invalid { what })?.is_some() {
                        (RetardedVelocity { past_m_s: chain_seed, ..delayed }, 0.0)
                    } else if below_edge {
                        declared = true;
                        (rest, 0.0)
                    } else {
                        return Err(FlueError::Invalid { what: "jet is slower than its retained forced-jet history" });
                    }
                }
                None if declared => (delayed, 0.0),
                None => (delayed, chain_seed),
            };
            let coefficients = JetStepCoefficients { q_max, half_width_m: b,
                offset_m: y0, band_alpha, hp_alpha, leak_rate,
                eta_override: declared.then_some(0.0), at_rest: omega_c == 0.0 };
            below_edge_samples += usize::from(declared);
            let boundary = self.line.mouth_boundary()?;
            // The opt-in breath noise: ONE draw per accepted sample, before any nonlinear trial, from its own stream, on
            // Verge's U_j of the instantaneous foot pressure (not on card_u, the loop's own jet speed).
            let (flue_height_m, gas_density, kinematic_viscosity) = (self.geometry.flue_height_m, self.rho, self.nu);
            let foot_pa = self.blow_pa;
            let noise_pa = self.jet_noise.as_mut()
                .and_then(|noise| noise.pressure_pa(foot_pa, flue_height_m, gas_density, kinematic_viscosity, dt));
            let (mouth, trial, separation_power) =
                self.solve_jet_mouth(boundary, delayed, chain_seed, coefficients, nolle_step, noise_pa)?;
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
            // The series breath noise delivers its own work: it is part of the source pressure the mouth sees.
            source_work += match noise_pa {
                Some(p) => (p_source + p) * q * dt,
                None => p_source * q * dt,
            };
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
            if self.history_moving > 0 || q / s_m != 0.0 {
                self.history_moving = self.history_moving.saturating_add(1);
            }
            self.reached_minimum_speed |= u >= self.minimum_jet_speed_m_s;
            // The accepted exit velocity into the forced jet's retained history.
            if let Some(kernel) = &mut self.forced {
                kernel.push(delayed.past_m_s + q / s_m);
            }
            // The observable is the wave launched into the pipe: the mouth
            // plane is the pipe's pressure node, so the total pressure
            // there says little about the standing wave the pipe carries.
            *slot = match &self.line {
                FlueResonator::Dg(duct) => duct.last_frame()
                    .map_or(p_plus, |frame| frame.outgoing_pressure_pa),
                _ => p_plus,
            };
            self.samples += 1;
            observe(&self.line, LabiumSample { eta_m: eta, offset_m: y0 })?;
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
            below_speed_edge_samples: below_edge_samples,
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
                let mut moved = false;
                for _ in 0..17 {
                    let previous = v.nolle_jet().unwrap().receptivity_displacement_m();
                    let response = v
                        .nolle_jet()
                        .unwrap()
                        .prepare_step(speed, geometry.flue_width_m, 0.0)
                        .unwrap();
                    let mut seen = None;
                    v.step_block_labium_observed(&mut [0.0], |_, labium| {
                        seen = Some(labium);
                        Ok(())
                    })
                    .unwrap();
                    let candidate = response
                        .candidate(v.mouth_flow_m3_s(), v.flow_prev - 0.5 * q_max)
                        .unwrap();
                    // The labium tap reads back the displacement the transfer imposed on the cut,
                    // not the idle island chain's.
                    let seen = seen.expect("one observed sample");
                    assert_eq!(seen.offset_m, geometry.effective_offset_m());
                    assert!(
                        (seen.eta_m - candidate.eta_m).abs()
                            <= 256.0 * f64::EPSILON * seen.eta_m.abs().max(candidate.eta_m.abs()),
                        "tapped {} against imposed {} at {rate} Hz, {pressure} Pa",
                        seen.eta_m,
                        candidate.eta_m
                    );
                    moved |= candidate.eta_m != 0.0;
                    let constitutive_cut = 0.5
                        * q_max
                        * (1.0
                            + ((candidate.eta_m - geometry.effective_offset_m())
                                / (0.5 * geometry.flue_height_m))
                                .tanh());
                    assert!(
                        (v.flow_prev - constitutive_cut).abs() <= 32.0 * f64::EPSILON * q_max,
                        "accepted cut-flow residual at {rate} Hz, {pressure} Pa"
                    );
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
                assert!(moved, "the imposed displacement moves at {rate} Hz, {pressure} Pa");
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
    fn g3_developing_windway_holds_shahs_speed_and_refuses_what_it_does_not_cover() {
        use nolle_jet::{ChannelClass, NolleSpec};
        let (air, g) = (gas(101_325.0), geometry());
        let held = DevelopingSlit::from_gas(&air, 0.015).unwrap()
            .held(g.flue_width_m, g.flue_height_m, 100.0).unwrap();
        held.validity().unwrap();
        let mut developing = voice(air, g, 0.01).with_developing_windway(0.015, 100.0).unwrap();
        let mut fully = voice(air, g, 0.01).with_viscous_windway(0.015).unwrap();
        for v in [&mut developing, &mut fully] { v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap(); }
        assert_eq!(developing.jet_speed_m_s(), held.bulk_speed_m_s);
        assert_eq!(developing.developing_windway().unwrap().1, &held);
        assert!(developing.jet_speed_m_s() < fully.jet_speed_m_s(),
            "entrance shear must hold a slower jet than the fully developed law");
        let mut out = [0.0; 257];
        developing.step_block(&mut out).unwrap();
        assert!(out.iter().all(|p| p.is_finite()));
        assert_work(&developing);
        assert!(developing.with_developing_windway(0.015, 100.0).is_err(), "attachment after samples");
        // One windway, attached cold, at a held pressure inside validity.
        assert!(voice(air, g, 0.01).with_viscous_windway(0.015).unwrap()
            .with_developing_windway(0.015, 100.0).is_err());
        assert!(voice(air, g, 0.01).with_developing_windway(0.015, 100.0).unwrap()
            .with_viscous_windway(0.015).is_err());
        assert!(voice(air, g, 0.01).with_developing_windway(0.015, 100.0).unwrap()
            .with_unsteady_windway(0.015, 8, 8).is_err());
        assert!(voice(air, g, 0.01).with_unsteady_windway(0.015, 8, 8).unwrap()
            .with_developing_windway(0.015, 100.0).is_err());
        for pressure in [0.0, -1.0, f64::NAN] {
            assert!(voice(air, g, 0.01).with_developing_windway(0.015, pressure).is_err());
        }
        // A 1 mm channel at 5 kPa holds x+ below 4e-4: refused, not extrapolated.
        assert!(matches!(voice(air, g, 0.01).with_developing_windway(0.001, 5000.0),
            Err(FlueError::Invalid { what }) if what.contains("x+")));
        // The Nolle/Price transfer takes this windway's exit when it is attached after it (see
        // g3_nolle_jet_takes_the_developing_exits_own_centre_over_bulk_not_the_poiseuille_constant);
        // attached before it, the windway cannot set the exit and is refused.
        let nolle = FlueGeometry { cut_up_m: 0.004, ..g };
        let spec = NolleSpec::price_2015(ChannelClass::Short, 0.0007).unwrap();
        assert!(voice(air, nolle, 1.0).with_nolle_jet(spec).unwrap()
            .with_developing_windway(0.015, 100.0).is_err());
    }

    #[test]
    fn g3_nolle_jet_takes_the_developing_exits_own_centre_over_bulk_not_the_poiseuille_constant() {
        use nolle_jet::{ChannelClass, NolleSpec};
        let air = gas(101_325.0);
        let g = FlueGeometry { cut_up_m: 0.004, ..geometry() };
        let spec = NolleSpec::price_2015(ChannelClass::Short, 0.0007).unwrap();
        // Planted exit: a 5 mm slit at 100 Pa holds u_c/U near 1.27 with its layers not yet met, so
        // it is neither the uniform exit (1) nor the Poiseuille exit (1.5) the transfer used to assume.
        let held = DevelopingSlit::from_gas(&air, 0.005).unwrap()
            .held(g.flue_width_m, g.flue_height_m, 100.0).unwrap();
        held.validity().unwrap();
        assert!(!held.developed);
        let ratio = held.core_speed_m_s / held.bulk_speed_m_s;
        assert!(ratio > 1.1 && ratio < 1.4, "planted core over bulk {ratio}");
        let mut v = voice(air, g, 1.0).with_developing_windway(0.005, 100.0).unwrap()
            .with_nolle_jet(spec).unwrap();
        assert_eq!(v.nolle_jet().unwrap().centre_over_bulk(), ratio);
        // The transit follows the centre speed U0 = ratio x bulk. The constant 1.5 would be faster
        // than the exit's own core speed, which is bounded by sqrt(2P/rho).
        v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
        assert_eq!(v.jet_speed_m_s(), held.bulk_speed_m_s);
        let jet = v.nolle_jet().unwrap();
        assert_eq!(v.transit_s(), jet.transit_s(ratio * held.bulk_speed_m_s));
        let constant = jet.transit_s(nolle_jet::POISEUILLE_CENTRE_OVER_BULK * held.bulk_speed_m_s);
        assert!(constant < v.transit_s(), "the constant 1.5 shortens the transit: {constant}");
        assert!(nolle_jet::POISEUILLE_CENTRE_OVER_BULK * held.bulk_speed_m_s > held.core_speed_m_s,
            "the constant centre speed exceeds the exit's Bernoulli core speed");
        let mut out = [0.0; 257];
        v.step_block(&mut out).unwrap();
        assert!(out.iter().all(|p| p.is_finite()));
        assert_work(&v);
        // A slit whose layers met before the exit is exactly the Poiseuille exit.
        let developed = DevelopingSlit::from_gas(&air, 0.03).unwrap()
            .held(g.flue_width_m, g.flue_height_m, 100.0).unwrap();
        assert!(developed.developed);
        let v = voice(air, g, 1.0).with_developing_windway(0.03, 100.0).unwrap()
            .with_nolle_jet(spec).unwrap();
        assert_eq!(v.nolle_jet().unwrap().centre_over_bulk(), nolle_jet::POISEUILLE_CENTRE_OVER_BULK);
        // Where the law refuses the held state the windway refuses, so no transfer is ever built on
        // the constant: a 300 mm slit at 1.5 kPa is developed at Re_Dh near 2500, above 1760.
        assert!(matches!(voice(air, g, 1.0).with_developing_windway(0.3, 1500.0),
            Err(FlueError::Invalid { what }) if what.contains("Re_Dh 1760")));
        // The slit and no-windway exits keep their own constants.
        let slit = voice(air, g, 1.0).with_viscous_windway(0.015).unwrap().with_nolle_jet(spec).unwrap();
        assert_eq!(slit.nolle_jet().unwrap().centre_over_bulk(), nolle_jet::POISEUILLE_CENTRE_OVER_BULK);
        let bare = voice(air, g, 1.0).with_nolle_jet(spec).unwrap();
        assert_eq!(bare.nolle_jet().unwrap().centre_over_bulk(), 1.0);
    }

    #[test]
    fn g3_nolle_mach_admission_reads_the_centre_speed_behind_a_slit() {
        use nolle_jet::{ChannelClass, NolleSpec};
        let air = gas(101_325.0);
        let g = FlueGeometry { cut_up_m: 0.004, ..geometry() };
        let spec = NolleSpec::price_2015(ChannelClass::Short, 0.0007).unwrap();
        let cap = 0.3 * air.sound_speed; // the fixture voice declares maximum_jet_mach 0.3
        // Premise: at 4100 Pa the fully developed slit holds a bulk speed under the cap whose Poiseuille
        // centre (1.5 x bulk, the speed the transfer runs on) is over it.
        let mut bare = voice(air, g, 1.0).with_viscous_windway(0.015).unwrap();
        bare.apply(FlueControl::SetBlowingPressure(4100.0)).unwrap();
        let bulk = bare.jet_speed_m_s();
        assert!(bulk < cap && nolle_jet::POISEUILLE_CENTRE_OVER_BULK * bulk > cap,
            "premise: bulk {bulk} m/s, cap {cap} m/s");
        // With the transfer attached the admission binds the centre speed: a pressure whose centre speed
        // is inside the domain is admitted, 4100 Pa is refused and leaves the accepted input intact.
        let mut v = voice(air, g, 1.0).with_viscous_windway(0.015).unwrap().with_nolle_jet(spec).unwrap();
        v.apply(FlueControl::SetBlowingPressure(2300.0)).unwrap();
        assert!(v.card_speed_m_s() < cap);
        let accepted = (v.jet_speed_m_s(), v.card_speed_m_s());
        assert!(matches!(v.apply(FlueControl::SetBlowingPressure(4100.0)),
            Err(FlueError::Invalid { what }) if what.contains("low-Mach")));
        assert_eq!((v.jet_speed_m_s(), v.card_speed_m_s()), accepted);
    }

    #[test]
    fn g3_card_centreline_speed_feeds_the_card_the_exit_centre_and_keeps_the_flux() {
        use nolle_jet::{ChannelClass, NolleSpec};
        let (air, g) = (gas(101_325.0), geometry());
        // A uniform exit's bulk speed is its centre speed: the opt-in changes no bit.
        let mut ideal = voice(air, g, 0.01);
        let mut centred = voice(air, g, 0.01).with_card_centreline_speed(100.0).unwrap();
        for v in [&mut ideal, &mut centred] { v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap(); }
        assert_eq!(centred.card_speed_m_s(), centred.jet_speed_m_s());
        let (mut a, mut b) = ([0.0; 257], [0.0; 257]);
        ideal.step_block(&mut a).unwrap();
        centred.step_block(&mut b).unwrap();
        assert_eq!(a, b);
        // An undeveloped exit feeds the card its inviscid core; the bulk speed keeps the flux.
        let mut developing = voice(air, g, 0.01).with_developing_windway(0.015, 100.0).unwrap()
            .with_card_centreline_speed(100.0).unwrap();
        developing.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
        let held = *developing.developing_windway().unwrap().1;
        assert!(!held.developed);
        assert_eq!(developing.jet_speed_m_s(), held.bulk_speed_m_s);
        assert_eq!(developing.card_speed_m_s(), held.core_speed_m_s);
        assert_eq!(developing.card_centreline_speed().unwrap().held_card_speed_m_s, held.core_speed_m_s);
        assert!((developing.transit_s()*held.core_speed_m_s*developing.island.convection_ratio
            -g.cut_up_m).abs() < 1e-18);
        assert!((developing.reynolds()-held.core_speed_m_s*g.flue_height_m/developing.nu).abs()
            < 1e-9*developing.reynolds());
        let mut out = [0.0; 257];
        developing.step_block(&mut out).unwrap();
        assert!(out.iter().all(|p| p.is_finite()));
        assert_work(&developing);
        // A fully developed slit feeds the Poiseuille centre while it is under the inviscid
        // bound, and is refused where its own profile would break it.
        let mut slit = voice(air, g, 0.01).with_viscous_windway(0.06).unwrap()
            .with_card_centreline_speed(100.0).unwrap();
        slit.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
        assert_eq!(slit.card_speed_m_s(), 1.5*slit.jet_speed_m_s());
        assert!(slit.card_speed_m_s() < (2.0*100.0/air.density).sqrt());
        assert!(matches!(voice(air, g, 0.01).with_viscous_windway(0.015).unwrap()
            .with_card_centreline_speed(100.0), Err(FlueError::Invalid { what }) if what.contains("inviscid")));
        // One cold attachment, after the windway, at a positive held pressure.
        assert!(developing.with_card_centreline_speed(100.0).is_err(), "attachment after samples");
        assert!(voice(air, g, 0.01).with_card_centreline_speed(100.0).unwrap()
            .with_card_centreline_speed(100.0).is_err());
        for pressure in [0.0, -1.0, f64::NAN] {
            assert!(voice(air, g, 0.01).with_card_centreline_speed(pressure).is_err());
        }
        assert!(voice(air, g, 0.01).with_card_centreline_speed(100.0).unwrap()
            .with_viscous_windway(0.06).is_err());
        assert!(voice(air, g, 0.01).with_card_centreline_speed(100.0).unwrap()
            .with_developing_windway(0.015, 100.0).is_err());
        assert!(voice(air, g, 0.01).with_card_centreline_speed(100.0).unwrap()
            .with_unsteady_windway(0.015, 8, 8).is_err());
        assert!(voice(air, g, 0.01).with_unsteady_windway(0.06, 8, 8).unwrap()
            .with_card_centreline_speed(100.0).is_ok());
        assert!(matches!(voice(air, g, 0.01).with_unsteady_windway(0.015, 8, 8).unwrap()
            .with_card_centreline_speed(100.0), Err(FlueError::Invalid { what }) if what.contains("inviscid")));
        // The jet transfers own their exit speed, in either order.
        let nolle = FlueGeometry { cut_up_m: 0.004, ..g };
        let spec = NolleSpec::price_2015(ChannelClass::Short, 0.0007).unwrap();
        assert!(voice(air, nolle, 1.0).with_card_centreline_speed(100.0).unwrap()
            .with_nolle_jet(spec).is_err());
        assert!(voice(air, nolle, 1.0).with_nolle_jet(spec).unwrap()
            .with_card_centreline_speed(100.0).is_err());
        assert!(voice(air, nolle, 1.0).with_card_centreline_speed(100.0).unwrap()
            .with_forced_jet(forced_jet::COARSE_SPEC, 100.0).is_err());
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
mod forced_jet_tests {
    use super::*;
    use crate::bernoulli_aperture::{contraction::ContractionSpec,
        network::{NetworkNode, TubeNetworkSpec, TubeSection}};
    use fs_exec::CancelGate;
    use fs_material::gas::GasSpec;

    // Plumbing of the opt-in forced jet on the coarse row (forced_jet::COARSE_SPEC, 0.25 b and
    // 40 b/U0: same physics, coarse grid); these are not qualified flute/recorder specimens.
    const PLUMBING: forced_jet::ForcedJetSpec = forced_jet::COARSE_SPEC;
    // A resolved cut-up, and one a 4 kPa jet crosses within one 48 kHz sample. On a 0.1 mm flue
    // that cut-up is 5.3 b, so the jet grows within the sample and its current feedback exceeds
    // one; shorter stations stay below one (measured over flue 0.1-1 mm, cut-up 0.06-0.8 mm,
    // 100-4000 Pa: only 5.3 b and 10.7 b at 4 kPa refuse).
    const CUT_UP_M: f64 = 0.004;
    const SUBSAMPLE_CUT_UP_M: f64 = 4e-4;
    const SUBSAMPLE_FLUE_M: f64 = 1e-4;
    const SUBSAMPLE_PA: f64 = 4000.0;

    fn geometry(cut_up_m: f64) -> FlueGeometry {
        FlueGeometry { flue_width_m: 0.010, flue_height_m: 0.001, cut_up_m,
            labium_offset_m: 0.0, jet_angle_rad: 0.0,
            bore_area_m2: core::f64::consts::PI*0.006*0.006, mouth_loss_rel: 0.02 }
    }
    fn gas(pressure: f64) -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
    }
    // The chain's authored gain is passed and unused once the forced jet attaches.
    fn bare_voice(gas: GasState, geometry: FlueGeometry, rate: u32, history: usize) -> FlueVoice {
        let source = network::FlueNetworkSpec {
            graph: TubeNetworkSpec {
                nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection: -0.7 }],
                sections: vec![TubeSection { nodes: [0,1], length_m: 0.3, radius_m: 0.006,
                    max_length_error_m: 0.005 }],
                sound_speed_m_s: gas.sound_speed, max_wave_memory_bytes: 1_000_000,
            },
            gases: vec![gas], viscothermal: vec![], minimum_jet_speed_m_s: 1.0,
            maximum_jet_mach: 0.3, max_jet_history_samples: history,
        };
        FlueVoice::from_network(Some(&fs_aeroac::jetcard::mint_tonal_interim_card().unwrap()),
            90.0, 1e-6, geometry, source, rate, 3, 0, &CancelGate::new_clock_free()).unwrap()
    }
    fn voice(gas: GasState, geometry: FlueGeometry) -> FlueVoice {
        bare_voice(gas, geometry, 48_000, 4096).with_forced_jet(PLUMBING, 100.0).unwrap()
    }
    fn contraction() -> ContractionSpec {
        ContractionSpec { area_ratio: 0.8, maximum_contracted_mach: 0.3,
            provenance: "authored numerical forced-jet fixture".into() }
    }
    fn assert_work(v: &FlueVoice) {
        for d in v.diagnostics() {
            let scale = d.source_work_j.abs()+d.pipe_work_j.abs()+d.stored_delta_j.abs()
                +d.radiation_loss_j+d.numerical_dissipation_j+d.separation_loss_j;
            assert!(d.ledger_defect_j.abs() <= 2e-12*scale+1e-30, "{d:?}");
            assert!(d.separation_loss_j >= 0.0);
        }
    }

    // Breath noise (JetNoiseSpec): G0 byte-identity, seed, refusals; G3 the rho U^2 law, the gate and the closed ledger.
    fn noise_spec(onset: f64) -> JetNoiseSpec {
        JetNoiseSpec { seed: 11, gain: 1.3e-3, strouhal_lo: 0.045, strouhal_hi: 0.115, reynolds_onset: onset }
    }
    // Eight 257-sample blocks of the relay chain at 100 Pa, as raw bits; the work ledger must close either way.
    fn noise_run(noise: Option<JetNoiseSpec>) -> Vec<u64> {
        let mut v = bare_voice(gas(101_325.0), geometry(CUT_UP_M), 48_000, 4096);
        if let Some(spec) = noise {
            v = v.with_jet_noise(spec).unwrap();
        }
        v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
        let mut bits = Vec::new();
        for _ in 0..8 {
            let mut out = [0.0; 257];
            v.step_block(&mut out).unwrap();
            bits.extend(out.iter().map(|x| x.to_bits()));
        }
        assert_work(&v);
        bits
    }

    #[test]
    fn g0_jet_noise_without_the_row_or_with_a_shut_gate_is_bit_identical_and_an_open_gate_is_seeded() {
        let base = noise_run(None);
        // A source whose gate never opens draws only from its own stream: the loop is the original loop, bit for bit.
        // (A noise drawn from the island's onset-seed stream would move every later seed draw and fail here.)
        assert_eq!(noise_run(Some(noise_spec(1.0e12))), base);
        let open = noise_run(Some(noise_spec(0.0)));
        assert_ne!(open, base, "an open source changes the loop");
        assert_eq!(noise_run(Some(noise_spec(0.0))), open, "a declared seed gives the same bits");
        let other = JetNoiseSpec { seed: 12, ..noise_spec(0.0) };
        assert_ne!(noise_run(Some(other)), open, "another seed gives another noise");
    }

    #[test]
    fn g0_jet_noise_refuses_unusable_fields_and_late_or_repeated_attachment() {
        let fresh = || bare_voice(gas(101_325.0), geometry(CUT_UP_M), 48_000, 4096);
        for bad in [
            JetNoiseSpec { gain: 0.0, ..noise_spec(1700.0) },
            JetNoiseSpec { gain: f64::NAN, ..noise_spec(1700.0) },
            JetNoiseSpec { strouhal_lo: 0.0, ..noise_spec(1700.0) },
            JetNoiseSpec { strouhal_lo: 0.115, strouhal_hi: 0.045, ..noise_spec(1700.0) },
            JetNoiseSpec { strouhal_hi: f64::INFINITY, ..noise_spec(1700.0) },
            JetNoiseSpec { reynolds_onset: -1.0, ..noise_spec(1700.0) },
        ] {
            assert!(fresh().with_jet_noise(bad).is_err(), "{bad:?}");
        }
        let attached = fresh().with_jet_noise(noise_spec(1700.0)).unwrap();
        assert_eq!(attached.jet_noise(), Some(&noise_spec(1700.0)));
        assert!(attached.with_jet_noise(noise_spec(1700.0)).is_err(), "one attachment");
        let mut late = fresh();
        late.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
        late.step_block(&mut [0.0; 16]).unwrap();
        assert!(late.with_jet_noise(noise_spec(1700.0)).is_err(), "cold attachment only");
        assert!(fresh().jet_noise().is_none());
    }

    #[test]
    fn g3_jet_noise_rms_is_gain_times_the_foot_pressure_and_so_scales_with_the_square_of_verges_jet_velocity() {
        // U_j = sqrt(2 P_f / rho): 100 Pa -> 12.9 m/s, 400 Pa -> 25.8 m/s (U_j doubles); rms = gain * rho U_j^2 / 2 = gain * P_f.
        let (h, rho, nu, dt, count) = (1.0e-3, 1.2, 1.5e-5, 1.0 / 48_000.0, 1usize << 16);
        let rms = |foot_pa: f64| {
            let mut n = JetNoise::new(noise_spec(0.0)).unwrap();
            for _ in 0..4096 {
                n.pressure_pa(foot_pa, h, rho, nu, dt);
            }
            let s2: f64 = (0..count).map(|_| n.pressure_pa(foot_pa, h, rho, nu, dt).unwrap().powi(2)).sum();
            (s2 / count as f64).sqrt()
        };
        let (slow, fast) = (rms(100.0), rms(400.0));
        let law = 1.3e-3 * 100.0;
        assert!((slow - law).abs() < 0.05 * law, "rms {slow} against gain P_f = {law}");
        assert!((fast / slow - 4.0).abs() < 0.2, "doubling U_j gave x{}", fast / slow);
    }

    #[test]
    fn g3_jet_noise_gate_sits_at_the_onset_reynolds_number_of_verges_jet_velocity_and_the_ledger_closes_with_it_open() {
        // rho 1.2, nu 1.5e-5, h 1 mm: Re_h = 1700 at U_j = 25.5 m/s, P_f = 390.15 Pa.
        let (h, rho, nu, dt) = (1.0e-3, 1.2, 1.5e-5, 1.0 / 48_000.0);
        let mut n = JetNoise::new(noise_spec(1700.0)).unwrap();
        for _ in 0..64 {
            assert!(n.pressure_pa(390.0, h, rho, nu, dt).is_none(), "Re_h 1699.7 is below his floor");
            assert_eq!(n.last_pa, 0.0);
        }
        assert!(n.pressure_pa(390.5, h, rho, nu, dt).is_some(), "Re_h 1700.8 is above it");
        for dead in [0.0, -5.0, f64::NAN] {
            assert!(n.pressure_pa(dead, h, rho, nu, dt).is_none(), "no jet, no noise ({dead})");
        }
        // `noise_run` asserts the work ledger on its own: with the source open the series pressure's work is booked.
        let _ = noise_run(Some(JetNoiseSpec { gain: 0.05, ..noise_spec(0.0) }));
    }

    #[test]
    fn g3_jet_noise_follows_the_foot_pressure_not_the_windway_loss() {
        // MB 15158: at fixed P_f a card with a different windway loss draws the same noise, bit for bit; the slit slows
        // the jet the LOOP sees (its bulk speed) and nothing else. Block of one sample so the last draw is visible.
        let run = |slit: bool| {
            let mut v = bare_voice(gas(101_325.0), geometry(CUT_UP_M), 48_000, 4096).with_jet_noise(noise_spec(0.0)).unwrap();
            if slit {
                v = v.with_viscous_windway(0.0615).unwrap();
            }
            v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
            let mut draws = Vec::new();
            for _ in 0..512 {
                v.step_block(&mut [0.0; 1]).unwrap();
                draws.push(v.jet_noise.as_ref().unwrap().last_pa.to_bits());
            }
            (v.jet_speed_m_s(), draws)
        };
        let ((free_speed, free_draws), (slit_speed, slit_draws)) = (run(false), run(true));
        assert!(free_speed - slit_speed > 3.0, "the slit must slow the loop's jet: {free_speed} against {slit_speed} m/s");
        assert!(free_draws.iter().any(|b| f64::from_bits(*b) != 0.0), "the source must be open");
        assert_eq!(free_draws, slit_draws, "same P_f, same noise");
    }

    #[test]
    fn g0_forced_jet_attaches_once_after_the_windway_and_refuses_sub_grid_stations() {
        let v = voice(gas(101_325.0), geometry(CUT_UP_M));
        assert!(v.forced_jet().is_some());
        assert!(v.with_viscous_windway(0.015).is_err());
        let v = voice(gas(101_325.0), geometry(CUT_UP_M));
        assert!(v.with_unsteady_windway(0.015, 8, 8).is_err());
        let v = voice(gas(101_325.0), geometry(CUT_UP_M));
        assert!(v.with_forced_jet(PLUMBING, 100.0).is_err());
        let v = bare_voice(gas(101_325.0), geometry(CUT_UP_M), 48_000, 4096)
            .with_unsteady_windway(0.015, 8, 8).unwrap();
        assert!(v.with_forced_jet(PLUMBING, 100.0).is_err());
        // The developing windway sets the forced jet's exit, so it is attached first.
        let v = voice(gas(101_325.0), geometry(CUT_UP_M));
        assert!(matches!(v.with_developing_windway(0.015, 100.0),
            Err(FlueError::Invalid { what }) if what.contains("before the forced jet")));
        // A 60 um cut-up on a 1 mm flue is 0.08 b, inside the first 0.25 b grid cell.
        let v = bare_voice(gas(101_325.0), geometry(6e-5), 48_000, 4096);
        assert!(v.with_forced_jet(PLUMBING, 100.0).is_err());
    }

    #[test]
    fn g3_forced_jet_takes_the_developing_exits_core_speed_and_momentum_flux() {
        // A planted exit whose core speed differs from both the bulk and the Poiseuille centre:
        // a 5 mm slit at 100 Pa holds u_c/U near 1.27. A forced jet that ignored the profile
        // would take the uniform (1, 0.75 h) or Poiseuille (1.5, 0.4 h) exit.
        let (air, g) = (gas(101_325.0), geometry(CUT_UP_M));
        let held = DevelopingSlit::from_gas(&air, 0.005).unwrap()
            .held(g.flue_width_m, g.flue_height_m, 100.0).unwrap();
        held.validity().unwrap();
        assert!(!held.developed);
        let ratio = held.core_speed_m_s / held.bulk_speed_m_s;
        assert!(ratio > 1.1 && ratio < 1.4, "planted core over bulk {ratio}");
        let mut v = bare_voice(air, g, 48_000, 4096).with_developing_windway(0.005, 100.0)
            .unwrap().with_forced_jet(PLUMBING, 100.0).unwrap();
        let kernel = v.forced_jet().unwrap();
        assert_eq!(kernel.centre_over_bulk(), ratio);
        // b carries the momentum flux rho u_c^2 (h - 2 delta* (1 + 1/2.5)) at the core speed.
        let b = 0.75 * (g.flue_height_m - 2.8 * held.displacement_thickness_m);
        assert!((kernel.half_width_m() - b).abs() <= 1e-12 * b);
        assert!(kernel.half_width_m() > 0.45 * g.flue_height_m
            && kernel.half_width_m() < 0.6 * g.flue_height_m);
        let station = g.cut_up_m / b;
        assert!((kernel.step_response().station() - station).abs() <= 1e-12 * station);
        let reynolds = held.core_speed_m_s * b / v.nu;
        assert!((kernel.reynolds() - reynolds).abs() <= 1e-12 * reynolds);
        v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
        let mut out = [0.0; 257];
        v.step_block(&mut out).unwrap();
        assert!(out.iter().all(|p| p.is_finite()) && out.iter().any(|p| *p != 0.0));
        assert_work(&v);
        // A slit whose layers have met before the exit is the Poiseuille exit, exactly.
        let developed = DevelopingSlit::from_gas(&air, 0.03).unwrap()
            .held(g.flue_width_m, g.flue_height_m, 100.0).unwrap();
        assert!(developed.developed);
        let v = bare_voice(air, g, 48_000, 4096).with_developing_windway(0.03, 100.0)
            .unwrap().with_forced_jet(PLUMBING, 100.0).unwrap();
        assert_eq!(v.forced_jet().unwrap().centre_over_bulk(),
            forced_jet::POISEUILLE_CENTRE_OVER_BULK);
        assert_eq!(v.forced_jet().unwrap().half_width_m(),
            forced_jet::POISEUILLE_HALF_WIDTH_OVER_HEIGHT * g.flue_height_m);
        // The exit is held at the reference pressure, inside the law's validity: at 3 kPa the
        // 5 mm slit's x+ falls below 4e-4.
        let v = bare_voice(air, g, 48_000, 4096).with_developing_windway(0.005, 100.0).unwrap();
        assert!(matches!(v.with_forced_jet(PLUMBING, 3000.0),
            Err(FlueError::Invalid { what }) if what.contains("x+")));
    }

    #[test]
    fn g3_forced_jet_card_speed_is_the_exit_centre_speed_for_seed_mach_and_reynolds() {
        // The same planted exit as above: a 5 mm slit at 100 Pa, u_c/U near 1.27. The card's speed
        // (the onset seed, the Mach admission and the Reynolds number all read it in advance_block)
        // is the jet's centre speed, not the bulk speed that carries the flux. The seed has no
        // reading of its own; it is tested through this shared value.
        let (air, g) = (gas(101_325.0), geometry(CUT_UP_M));
        let held = DevelopingSlit::from_gas(&air, 0.005).unwrap()
            .held(g.flue_width_m, g.flue_height_m, 100.0).unwrap();
        let ratio = held.core_speed_m_s / held.bulk_speed_m_s;
        assert!(ratio > 1.1 && ratio < 1.4, "planted core over bulk {ratio}");
        let mut v = bare_voice(air, g, 48_000, 4096).with_developing_windway(0.005, 100.0)
            .unwrap().with_forced_jet(PLUMBING, 100.0).unwrap();
        v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
        assert_eq!(v.jet_speed_m_s(), held.bulk_speed_m_s);
        let centre = ratio * held.bulk_speed_m_s;
        assert_eq!(v.card_speed_m_s(), centre);
        assert!(v.card_speed_m_s() > v.jet_speed_m_s());
        assert_eq!(v.reynolds(), centre * g.flue_height_m / v.nu);
        let mut out = [0.0; 257];
        v.step_block(&mut out).unwrap();
        let expected = centre * g.flue_height_m / v.nu;
        let (lo, hi) = v.diagnostics().last().unwrap().reynolds_range;
        assert!((lo - expected).abs() <= 1e-9 * expected && (hi - expected).abs() <= 1e-9 * expected,
            "per-sample Reynolds range ({lo}, {hi}) against {expected}");
        assert_work(&v);
    }

    #[test]
    fn g3_card_speed_stays_the_bulk_speed_where_no_jet_owns_a_non_uniform_exit() {
        // The relay chain, with and without a slit windway (pexfp.135.5 D, untouched here), and a forced
        // jet on a uniform exit (ratio exactly 1): the card's speed is the bulk speed bit for bit.
        let (air, g) = (gas(101_325.0), geometry(CUT_UP_M));
        for mut v in [
            bare_voice(air, g, 48_000, 4096),
            bare_voice(air, g, 48_000, 4096).with_viscous_windway(0.015).unwrap(),
            bare_voice(air, g, 48_000, 4096).with_forced_jet(PLUMBING, 100.0).unwrap(),
        ] {
            v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
            assert_eq!(v.card_speed_m_s(), v.jet_speed_m_s());
        }
    }

    #[test]
    fn g3_forced_jet_preserves_work_and_block_partition_across_gas_states() {
        for pressure in [50_000.0,100_000.0,200_000.0] {
            for separated in [false,true] {
                let build = || {
                    let v = voice(gas(pressure), geometry(CUT_UP_M));
                    if separated { v.with_mouth_contraction(contraction()).unwrap() } else { v }
                };
                let (mut one, mut many) = (build(), build());
                for v in [&mut one, &mut many] {
                    v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
                }
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
    fn g3_forced_jet_takes_the_actual_viscous_windway_speed() {
        let ordinary = gas(101_325.0);
        let mut spec = GasSpec::dry_air_ussa1976();
        spec.sutherland_beta *= 4.0; // Synthetic viscosity-only substitution.
        let viscous = GasState::try_new(&spec, 293.15, ordinary.pressure).unwrap();
        assert_eq!(ordinary.density, viscous.density);
        assert_eq!(ordinary.sound_speed, viscous.sound_speed);
        let build = |gas| bare_voice(gas, geometry(CUT_UP_M), 48_000, 4096)
            .with_viscous_windway(0.015).unwrap().with_forced_jet(PLUMBING, 100.0).unwrap();
        let (mut a, mut b) = (build(ordinary), build(viscous));
        for v in [&mut a, &mut b] { v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap(); }
        assert!(b.jet_speed_m_s() < a.jet_speed_m_s());
        assert!(b.forced_jet().unwrap().reynolds() < a.forced_jet().unwrap().reynolds());
        let (mut first, mut second) = ([0.0;257], [0.0;257]);
        a.step_block(&mut first).unwrap(); b.step_block(&mut second).unwrap();
        assert_ne!(first,second);
        assert_work(&a); assert_work(&b);
    }

    #[test]
    fn g0_forced_jet_contraction_cap_rounding_cannot_refuse_an_interior_jet_root() {
        // Mouth area 0.25 mm x 0.2 mm, the same bits as 1 mm x 50 um; the station 0.2/0.75 is
        // inside the forced jet's grid. At 400 Pa the full jet flux is 1.6 times the contraction's
        // Mach cap, so the cap bounds the root search.
        let mut geometry = geometry(0.0002);
        geometry.flue_width_m = 0.00025;
        geometry.mouth_loss_rel = 0.0;
        // The jet crosses this cut-up in about a sample at 48 kHz, where the forced jet's
        // current feedback is underresolved; 960 kHz resolves it.
        let mut v = bare_voice(gas(101_325.0), geometry, 960_000, 65_536)
            .with_forced_jet(PLUMBING, 400.0).unwrap()
            .with_mouth_contraction(contraction()).unwrap();
        let loss = v.mouth_contraction().unwrap();
        let raw_bound = loss.effective_area_m2()
            *(loss.gas().sound_speed*loss.spec().maximum_contracted_mach);
        // Arithmetic witness: the nominal endpoint rounds above Mach 0.3.
        assert!(loss.at_flow(raw_bound).is_err());
        v.apply(FlueControl::SetBlowingPressure(400.0)).unwrap();
        assert!(v.jet_speed_m_s()*geometry.flue_height_m*geometry.flue_width_m > raw_bound);
        v.step_block(&mut [0.0;257]).unwrap();
        assert!(v.mouth_flow_m3_s().abs() < raw_bound);
        v.mouth_contraction().unwrap().at_flow(v.mouth_flow_m3_s()).unwrap();
        assert!(v.diagnostics()[0].separation_loss_j > 0.0);
        assert_work(&v);
    }

    #[test]
    fn g0_forced_jet_underresolved_current_feedback_refuses_before_any_wave_commit() {
        // A jet that crosses the cut-up and grows within one sample feeds the current mouth
        // flow back by more than one: the solve has no unique root, so it asks for a finer clock.
        let mut geometry = geometry(SUBSAMPLE_CUT_UP_M);
        geometry.flue_height_m = SUBSAMPLE_FLUE_M;
        let mut v = bare_voice(gas(101_325.0), geometry, 48_000, 4096)
            .with_forced_jet(PLUMBING, SUBSAMPLE_PA).unwrap();
        v.apply(FlueControl::SetBlowingPressure(SUBSAMPLE_PA)).unwrap();
        assert!(v.transit_s() < v.dt);
        assert!(matches!(v.step_block(&mut [0.0]),
            Err(FlueError::Invalid { what }) if what.contains("refine")));
        assert_eq!(v.samples_rendered(), 0);
        assert!(v.network().unwrap().last_frame().is_none());
        assert_eq!(v.mouth_flow_m3_s(), 0.0);
        assert!(v.claim_check().is_err());
    }
}

#[cfg(test)]
mod speed_edge_tests {
    use super::*;
    use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
    use fs_exec::CancelGate;
    use fs_material::gas::GasSpec;

    // Ruling 8132 on a 48 kHz control clock: the onset's rest state and the declared state below
    // the minimum jet speed. Plumbing fixtures, not instrument specimens.
    const RATE: u32 = 48_000;
    const HOLD_PA: f64 = 300.0;

    fn voice(minimum_speed: f64, history: usize, forced: bool) -> FlueVoice {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
        let geometry = FlueGeometry { flue_width_m: 0.010, flue_height_m: 0.001, cut_up_m: 0.004,
            labium_offset_m: 0.0, jet_angle_rad: 0.0,
            bore_area_m2: core::f64::consts::PI*0.006*0.006, mouth_loss_rel: 0.02 };
        let source = network::FlueNetworkSpec {
            graph: TubeNetworkSpec {
                nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection: -0.7 }],
                sections: vec![TubeSection { nodes: [0,1], length_m: 0.3, radius_m: 0.006,
                    max_length_error_m: 0.005 }],
                sound_speed_m_s: gas.sound_speed, max_wave_memory_bytes: 1_000_000,
            },
            gases: vec![gas], viscothermal: vec![], minimum_jet_speed_m_s: minimum_speed,
            maximum_jet_mach: 0.3, max_jet_history_samples: history,
        };
        let v = FlueVoice::from_network(Some(&fs_aeroac::jetcard::mint_tonal_interim_card().unwrap()),
            90.0, 1e-6, geometry, source, RATE, 3, 0, &CancelGate::new_clock_free()).unwrap();
        if forced { v.with_forced_jet(forced_jet::COARSE_SPEC, HOLD_PA).unwrap() } else { v }
    }
    // One control and one sample at a time, so every sample carries its own pressure.
    fn play(v: &mut FlueVoice, pressures: impl Iterator<Item = f64>) -> Vec<f64> {
        pressures.map(|pressure| {
            v.apply(FlueControl::SetBlowingPressure(pressure)).unwrap();
            let mut out = [0.0];
            v.step_block(&mut out).unwrap();
            out[0]
        }).collect()
    }

    #[test]
    fn g1_onset_rest_state_is_unbounded_history_bit_for_bit() {
        // U_min 1 m/s retains a 10 ms transit; U_min 0.01 m/s retains 1 s, more than the 0.1 s
        // render, and its ramp never goes below it. On a smooth 1 s ramp from rest the first 96
        // samples run below 1 m/s, where the short voice's lookups leave its history; they precede
        // the first motion, so both voices give the same bits.
        let ramp = || (0..4800).map(|n| HOLD_PA*f64::from(n)/f64::from(RATE));
        for forced in [false, true] {
            let mut short = voice(1.0, 4096, forced);
            let mut long = voice(0.01, 262_144, forced);
            assert!(short.history.len() < 4800 && long.history.len() > 4800);
            let a = play(&mut short, ramp());
            let b = play(&mut long, ramp());
            assert_eq!(a, b, "forced {forced}");
            assert!(a.iter().any(|p| *p != 0.0));
            for v in [&short, &long] {
                assert!(v.diagnostics().iter().all(|d| d.below_speed_edge_samples == 0));
            }
        }
    }

    // Plays one sample per pressure and counts the declared samples, checking each: no deflection,
    // and the windway's flow split at zero deflection.
    fn declared_samples(v: &mut FlueVoice, pressures: impl Iterator<Item = f64>) -> usize {
        let half_width = 0.5*v.geometry.flue_height_m;
        let offset = v.geometry.effective_offset_m();
        let mut declared = 0;
        for pressure in pressures {
            play(v, core::iter::once(pressure));
            let d = v.diagnostics().last().unwrap();
            if d.below_speed_edge_samples == 1 {
                declared += 1;
                let u = v.jet_speed_m_s();
                assert!(u > 0.0 && u < v.minimum_jet_speed_m_s());
                assert_eq!(d.eta_rms_over_b, 0.0);
                let q_max = u*v.geometry.flue_height_m*v.geometry.flue_width_m;
                assert_eq!(v.flow_prev, 0.5*q_max*(1.0+((0.0-offset)/half_width).tanh()));
            }
        }
        declared
    }

    #[test]
    fn g1_release_below_the_minimum_speed_takes_the_declared_state() {
        // 10 ms up, 20 ms held, a smooth 100 ms release to rest.
        let pressures = || (0..480).map(|n| HOLD_PA*f64::from(n)/480.0)
            .chain(core::iter::repeat_n(HOLD_PA, 960))
            .chain((0..=4800).rev().map(|n| HOLD_PA*f64::from(n)/4800.0));
        // The default chain's last samples run below U_min 1 m/s with lookups beyond its history.
        // The jet is below 1 m/s only for P under 0.6 Pa: the last 9 nonzero samples, each under
        // 4.5 % of the held 22 m/s.
        let mut v = voice(1.0, 4096, false);
        assert_eq!(v.minimum_jet_speed_m_s(), 1.0);
        let declared = declared_samples(&mut v, pressures());
        assert!((1..=9).contains(&declared), "{declared}");
        assert_eq!(v.jet_speed_m_s(), 0.0);
        // The forced jet's retained intervals (a 30 ms window at 1 m/s) still reach its window all
        // the way to rest, so its release declares nothing and stays exact: bit for bit the voice
        // that retains 3 s.
        let mut short = voice(1.0, 4096, true);
        let mut long = voice(0.01, 262_144, true);
        assert_eq!(play(&mut short, pressures()), play(&mut long, pressures()));
        assert!(short.diagnostics().iter().all(|d| d.below_speed_edge_samples == 0));
        assert_eq!(short.jet_speed_m_s(), 0.0);
        // Held at 0.3 Pa (0.71 m/s) for 150 ms, longer than its ring of 1 m/s intervals reaches
        // the window at that speed, it takes the declared state once they fall short.
        let mut slow = voice(1.0, 4096, true);
        let held = (0..480).map(|n| HOLD_PA*f64::from(n)/480.0)
            .chain(core::iter::repeat_n(0.3, 7200));
        assert!(declared_samples(&mut slow, held) > 0);
    }

    #[test]
    fn a_jet_at_rest_holds_no_displacement_and_blows_again_from_rest() {
        // 10 ms up, 20 ms held, a smooth 100 ms release to rest, 200 ms at rest, then the same blow
        // again. At U = 0 the band, high-pass and leak corners are all zero; the chain used to hold
        // its last high-pass output and integrate it without leak, so the displacement grew
        // linearly for as long as the rest lasted. On every sample at rest after the first blow
        // the tapped displacement and the whole chain are exactly zero, so the second blow starts
        // from the rest state the first one did, and it sounds.
        let blow = || (0..480).map(|n| HOLD_PA*f64::from(n)/480.0)
            .chain(core::iter::repeat_n(HOLD_PA, 960));
        let pressures: Vec<f64> = blow()
            .chain((0..=4800).rev().map(|n| HOLD_PA*f64::from(n)/4800.0))
            .chain(core::iter::repeat_n(0.0, 9600))
            .chain(blow()).collect();
        let (hold_end, rest_end) = (1440, 1440 + 4801 + 9600);
        let mut v = voice(1.0, 4096, false);
        let (mut out, mut eta, mut at_rest) = (Vec::new(), Vec::new(), 0);
        for (n, &pressure) in pressures.iter().enumerate() {
            v.apply(FlueControl::SetBlowingPressure(pressure)).unwrap();
            let mut p = [0.0];
            let mut seen = None;
            v.step_block_labium_observed(&mut p, |_, s| { seen = Some(s); Ok(()) }).unwrap();
            let e = seen.expect("one observed sample").eta_m;
            if n >= hold_end && v.jet_speed_m_s() == 0.0 {
                at_rest += 1;
                assert_eq!(e, 0.0, "sample {n}");
                assert_eq!((v.band_state, v.band_hp, v.displacement_state),
                    ([0.0; 2], [0.0; 2], 0.0), "sample {n}");
            }
            out.push(p[0]);
            eta.push(e);
        }
        // The release's last sample, the 9600 at rest and the second blow's first (0 Pa).
        assert_eq!(at_rest, 9602);
        let peak = |r: core::ops::Range<usize>| eta[r].iter().fold(0.0_f64, |m, e| m.max(e.abs()));
        eprintln!("peak |eta|/b: held {:.6e}, release {:.6e}, second blow {:.6e}",
            peak(480..hold_end)/(0.5*v.geometry.flue_height_m),
            peak(hold_end..rest_end)/(0.5*v.geometry.flue_height_m),
            peak(rest_end..eta.len())/(0.5*v.geometry.flue_height_m));
        assert!(eta[rest_end..].iter().any(|e| *e != 0.0));
        assert!(out[rest_end..].iter().any(|p| *p != 0.0));
    }

    #[test]
    fn labium_tap_reads_the_displacement_each_accepted_cut_used() {
        // 10 ms up and 20 ms held. On the relay and the forced jet every tapped (eta, offset) gives
        // back the accepted jet flow bit for bit, and observing changes no sample. A jet never
        // blown carries no displacement: the tap reads exactly 0.
        let blown: Vec<f64> = (0..480).map(|n| HOLD_PA*f64::from(n)/480.0)
            .chain(core::iter::repeat_n(HOLD_PA, 960)).collect();
        for forced in [false, true] {
            for (pressures, rest) in [(blown.clone(), false), (vec![0.0; 480], true)] {
                let (mut plain, mut tapped) = (voice(1.0, 4096, forced), voice(1.0, 4096, forced));
                let (h, w) = (tapped.geometry.flue_height_m, tapped.geometry.flue_width_m);
                let mut moved = false;
                for pressure in pressures {
                    let (mut p, mut q) = ([0.0], [0.0]);
                    plain.apply(FlueControl::SetBlowingPressure(pressure)).unwrap();
                    tapped.apply(FlueControl::SetBlowingPressure(pressure)).unwrap();
                    plain.step_block(&mut p).unwrap();
                    let mut seen = None;
                    tapped.step_block_labium_observed(&mut q, |_, s| { seen = Some(s); Ok(()) }).unwrap();
                    assert_eq!(p[0].to_bits(), q[0].to_bits(), "forced {forced}");
                    let s = seen.expect("one observed sample");
                    assert_eq!(s.offset_m, tapped.geometry.effective_offset_m());
                    let q_max = tapped.jet_speed_m_s()*h*w;
                    let cut = if q_max > 0.0 {
                        0.5*q_max*(1.0+((s.eta_m-s.offset_m)/(0.5*h)).tanh())
                    } else { 0.0 };
                    assert_eq!(tapped.flow_prev.to_bits(), cut.to_bits(), "forced {forced}, {pressure} Pa");
                    if rest { assert_eq!(s.eta_m, 0.0, "forced {forced}"); }
                    moved |= s.eta_m != 0.0;
                }
                assert_eq!(moved, !rest, "forced {forced}");
            }
        }
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
