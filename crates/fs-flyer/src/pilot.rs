//! Two-axis PilotWrightModel (bead wf-root-guzez.5.16.2, E4.6c-ii).
//! Plan §5.1.4: a cue-based crossover-class pilot — FIXED gains drawn
//! from a PRE-REGISTERED family (family id + member index enter the
//! model identity; NO online adaptation), reaction delay in integer
//! perception ticks, first-order neuromuscular lag (exact exponential),
//! output saturation at the declared force/travel limits, and a
//! deterministic tick-addressed remnant on its own philox stream
//! (distinct kernel from perception's).
//!
//! Command chain (longitudinal): perceived pitch cues → desired canard
//! angle → PROPRIOCEPTIVE inner loop (the pilot feels the lever) →
//! lever force. Lateral: perceived roll cues → warp command (the 1903
//! rudder is slaved downstream, control-signs-v1). Both channels pass
//! through delay → lag → saturation, remnant injected pre-saturation.
//!
//! Claim boundary: V-02c1 validates the GENERIC mechanism on generic
//! controlled elements; H-02c against the Wright airframe is a
//! COMPATIBILITY record at Estimated — never a validation claim.

use crate::Refusal;
use crate::perception::{N_CUES, PERCEPTION_HZ, PerceivedCues};
use fs_math::det;
use fs_rand::{Stream, StreamKey};

/// Reaction-delay cap [ticks] (1 s at 120 Hz).
pub const MAX_REACTION_TICKS: usize = 120;

/// Registered family size (members 0..FAMILY_SIZE; index cap tested at
/// cap AND cap+1).
pub const FAMILY_SIZE: u32 = 3;

/// Pilot remnant draws per tick (one per axis).
pub const DRAWS_PER_TICK: u64 = 2;

/// One pre-registered gain member.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PilotGains {
    /// Family identity (enters PilotPopulationHypothesisId).
    pub family_id: &'static str,
    /// Member index within the family.
    pub member: u32,
    /// Desired canard angle per pitch-attitude error [rad/rad].
    pub k_theta: f64,
    /// Pitch-rate lead [rad per rad/s].
    pub k_q: f64,
    /// Proprioceptive lever position gain [N/rad].
    pub k_lever_pos: f64,
    /// Proprioceptive lever rate damping [N per rad/s].
    pub k_lever_rate: f64,
    /// Warp command per roll-attitude error [rad/rad].
    pub k_phi: f64,
    /// Roll-rate lead [rad per rad/s].
    pub k_p: f64,
    /// Reaction delay [perception ticks].
    pub reaction_ticks: usize,
    /// Neuromuscular lag time constant [s].
    pub neuromuscular_tau_s: f64,
    /// Lever force saturation [N].
    pub force_limit_n: f64,
    /// Warp command saturation [rad].
    pub warp_limit_rad: f64,
    /// Remnant one-sigma on the force channel [N].
    pub remnant_sigma_force_n: f64,
    /// Remnant one-sigma on the warp channel [rad].
    pub remnant_sigma_warp: f64,
}

/// The pre-registered v1 family (three members spanning aggressiveness;
/// declared BEFORE any historical calibration — the calibration subset
/// may select among them, never invent new ones).
///
/// # Errors
/// `pilot-member-invalid` (index at the family size — cap and cap+1).
pub fn pilot_family_v1(member: u32) -> Result<PilotGains, Refusal> {
    let base = PilotGains {
        family_id: "wright-pilot-family-v1",
        member,
        k_theta: 1.6,
        k_q: 1.2,
        k_lever_pos: 1200.0,
        k_lever_rate: 110.0,
        k_phi: 0.5,
        k_p: 0.3,
        reaction_ticks: 12, // 100 ms — trained-anticipatory class
        neuromuscular_tau_s: 0.08,
        force_limit_n: 220.0,
        warp_limit_rad: 0.148, // 8.5 deg (dossier warp limit)
        remnant_sigma_force_n: 1.5,
        remnant_sigma_warp: 0.002,
    };
    match member {
        0 => Ok(base),
        1 => Ok(PilotGains {
            k_theta: 2.4,
            k_q: 1.8,
            ..base
        }),
        2 => Ok(PilotGains {
            k_theta: 1.0,
            k_q: 0.8,
            reaction_ticks: 18,
            ..base
        }),
        _ => Err(Refusal {
            code: "pilot-member-invalid",
            message: format!("member {member} outside the registered family (size {FAMILY_SIZE})"),
            ranked_repairs: vec![
                "select a registered member; extending the family is a pre-registration event"
                    .into(),
            ],
        }),
    }
}

/// Pilot dynamic state.
#[derive(Clone, Debug, PartialEq)]
pub struct PilotState {
    // Delay rings on the two raw channel commands.
    ring_long: Vec<f64>,
    ring_lat: Vec<f64>,
    // Neuromuscular filter states.
    nm_long: f64,
    nm_lat: f64,
    tick: u64,
}

/// One tick's commands.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PilotCommand {
    /// Lever force [N] (+ = pull back = nose-up, control-signs-v1).
    pub lever_force_n: f64,
    /// Warp command [rad] (+ = right wing down).
    pub warp_cmd_rad: f64,
    /// Saturation flags (force, warp) — receipts, never silent clamps.
    pub saturated: [bool; 2],
    /// The tick.
    pub tick: u64,
}

/// The pilot model: gains + remnant stream identity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PilotWrightModel {
    /// Gains (registered member).
    pub gains: PilotGains,
    /// Remnant stream key (kernel distinct from perception).
    pub stream_key: StreamKey,
}

impl PilotWrightModel {
    /// Build from a registered member + seed.
    ///
    /// # Errors
    /// `pilot-member-invalid`; `pilot-gains-invalid`.
    pub fn new(member: u32, seed: u64) -> Result<Self, Refusal> {
        let gains = pilot_family_v1(member)?;
        let m = PilotWrightModel {
            gains,
            stream_key: StreamKey {
                seed,
                kernel: 0x4643_504C, // "FCPL" — pilot kernel id
                tile: member,
            },
        };
        m.admit()?;
        Ok(m)
    }

    /// Admit the model.
    ///
    /// # Errors
    /// `pilot-gains-invalid` (non-finite gains, reaction delay above
    /// [`MAX_REACTION_TICKS`] — cap AND cap+1, non-positive lag/limits,
    /// negative remnant).
    pub fn admit(&self) -> Result<(), Refusal> {
        let g = &self.gains;
        let finite = [
            g.k_theta,
            g.k_q,
            g.k_lever_pos,
            g.k_lever_rate,
            g.k_phi,
            g.k_p,
            g.neuromuscular_tau_s,
            g.force_limit_n,
            g.warp_limit_rad,
            g.remnant_sigma_force_n,
            g.remnant_sigma_warp,
        ]
        .iter()
        .all(|v| v.is_finite());
        let ok = finite
            && g.reaction_ticks <= MAX_REACTION_TICKS
            && g.neuromuscular_tau_s > 0.0
            && g.force_limit_n > 0.0
            && g.warp_limit_rad > 0.0
            && g.remnant_sigma_force_n >= 0.0
            && g.remnant_sigma_warp >= 0.0;
        if !ok {
            return Err(Refusal {
                code: "pilot-gains-invalid",
                message: format!("{g:?}"),
                ranked_repairs: vec![format!(
                    "finite gains; reaction <= {MAX_REACTION_TICKS}; positive lag/limits"
                )],
            });
        }
        Ok(())
    }

    /// Fresh state.
    ///
    /// # Errors
    /// Admission refusals.
    pub fn init(&self) -> Result<PilotState, Refusal> {
        self.admit()?;
        Ok(PilotState {
            ring_long: vec![0.0; self.gains.reaction_ticks + 1],
            ring_lat: vec![0.0; self.gains.reaction_ticks + 1],
            nm_long: 0.0,
            nm_lat: 0.0,
            tick: 0,
        })
    }

    /// One 120 Hz step. Cue order is the perception service's frozen
    /// order; `lever` is the proprioceptive lever state (deflection,
    /// rate); references are the commanded attitude targets.
    ///
    /// # Errors
    /// `pilot-input-invalid` (non-finite inputs);
    /// `pilot-stream-invalid`.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &self,
        state: &mut PilotState,
        cues: &PerceivedCues,
        lever_delta_rad: f64,
        lever_rate_rad_s: f64,
        theta_ref_rad: f64,
        phi_ref_rad: f64,
    ) -> Result<PilotCommand, Refusal> {
        if cues.values.iter().any(|v| !v.is_finite())
            || !lever_delta_rad.is_finite()
            || !lever_rate_rad_s.is_finite()
            || !theta_ref_rad.is_finite()
            || !phi_ref_rad.is_finite()
        {
            return Err(Refusal {
                code: "pilot-input-invalid",
                message: format!("cues {:?}", cues.values),
                ranked_repairs: vec!["check the perception/lever plumbing".into()],
            });
        }
        let g = &self.gains;
        let dt = 1.0 / PERCEPTION_HZ;
        // Raw channel commands (before delay/lag).
        let theta = cues.values[0];
        let q = cues.values[1];
        let phi = cues.values[3];
        let p = cues.values[4];
        let dc_desired = g.k_theta * (theta_ref_rad - theta) - g.k_q * q;
        let warp_raw = g.k_phi * (phi_ref_rad - phi) - g.k_p * p;
        // Reaction delay (exactly-D ring, perception.rs semantics).
        let d = g.reaction_ticks;
        let t = state.tick as usize;
        let (dc_delayed, warp_delayed) = if d == 0 {
            (dc_desired, warp_raw)
        } else {
            (
                state.ring_long[(t + 1) % (d + 1)],
                state.ring_lat[(t + 1) % (d + 1)],
            )
        };
        state.ring_long[t % (d + 1)] = dc_desired;
        state.ring_lat[t % (d + 1)] = warp_raw;
        // Neuromuscular lag (exact exponential).
        let a = det::exp(-dt / g.neuromuscular_tau_s);
        state.nm_long = dc_delayed + (state.nm_long - dc_delayed) * a;
        state.nm_lat = warp_delayed + (state.nm_lat - warp_delayed) * a;
        // Remnant (tick-addressed, pre-saturation).
        let mut stream = Stream::resume(fs_rand::StreamCheckpoint {
            checkpoint_version: fs_rand::STREAM_CHECKPOINT_VERSION,
            stream_semantics_version: fs_rand::STREAM_SEMANTICS_VERSION,
            key: self.stream_key,
            index: state.tick * DRAWS_PER_TICK,
        })
        .map_err(|e| Refusal {
            code: "pilot-stream-invalid",
            message: format!("{e:?}"),
            ranked_repairs: vec!["tick out of stream range".into()],
        })?;
        // Longitudinal: proprioceptive inner loop turns the desired
        // canard angle into lever force.
        let force_raw = g.k_lever_pos * (state.nm_long - lever_delta_rad)
            - g.k_lever_rate * lever_rate_rad_s
            + g.remnant_sigma_force_n * stream.next_normal();
        let warp_noisy = state.nm_lat + g.remnant_sigma_warp * stream.next_normal();
        let force = force_raw.clamp(-g.force_limit_n, g.force_limit_n);
        let warp = warp_noisy.clamp(-g.warp_limit_rad, g.warp_limit_rad);
        let out = PilotCommand {
            lever_force_n: force,
            warp_cmd_rad: warp,
            saturated: [
                force_raw.abs() > g.force_limit_n,
                warp_noisy.abs() > g.warp_limit_rad,
            ],
            tick: state.tick,
        };
        state.tick += 1;
        Ok(out)
    }
}

/// Convenience: pack raw physics cues in the frozen perception order.
#[must_use]
pub fn pack_cues(
    theta_rad: f64,
    q_rad_s: f64,
    heave_accel: f64,
    phi_rad: f64,
    p_rad_s: f64,
    r_rad_s: f64,
) -> [f64; N_CUES] {
    [theta_rad, q_rad_s, heave_accel, phi_rad, p_rad_s, r_rad_s]
}
