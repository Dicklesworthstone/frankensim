//! The Wright Flyer lifecycle engine (bead wf-root-guzez.6.2, E5.1 —
//! the native core the wasm surface wraps). init(scenario) →
//! HeldOnRailEquilibrated closure (E4.6d) → RunIntentId minted AFTER
//! the tick-0 digest → 120 Hz step loop:
//!
//!   OnRail: constrained acceleration under REAL thrust/drag from the
//!           coupled force build-up; release when lift carries weight;
//!           running off the rail without lift is a TerminalEvent.
//!   Airborne: perception (E4.6c-i) → control law (fixed hold /
//!           historical pilot / human input) → canard mechanism →
//!           full nonlinear force build-up with the pitch-rate term →
//!           FRD longitudinal kinematics + rotor torque balance + OU
//!           gust coupling + reduced lateral build-up → ground contact.
//!
//! Every tick chains into a running blake3 digest (bit-identical
//! lifecycles are a checkable claim), and the snapshot payload for the
//! E5.0 ring is a frozen layout. Tier notes carried honestly: the
//! mechanism runs the m_aero = 0 stick tier (H-02c convention), the OU
//! gust enters as an angle-of-attack increment (reduced coupling,
//! declared), and roll/yaw use the Estimated ReducedAeroelasticWarp +
//! ReducedLateralBuildUp tier (not a structural-margin claim).

use crate::Refusal;
use crate::aerowarp::ReducedAeroelasticWarp;
use crate::aircraft::{OpenLoopDesign, wright_openloop_v1};
use crate::assist::AssistSpec;
use crate::canardmech::{CANARD_MECH_V1, CanardMechanism, MechState};
use crate::equilibrate::{EquilibrationSpec, PrelaunchBranch, Tick0State, equilibrate};
use crate::lateral::{LateralModel, LateralState, RudderLinkage};
use crate::longitudinal::IYY_KG_M2;
use crate::perception::{PERCEPTION_HZ, PerceptionModelSpec, PerceptionState, perception_v1};
use crate::pilot::{PilotState, PilotWrightModel, pack_cues};
use fs_atmo::ou::{OuMode, STATIONARY_ANCHOR_TICK, StationaryOuPath};
use fs_blake3::hash_domain;
use fs_math::det;

/// Tick cap (absurd-input guard; 10 minutes at 120 Hz).
pub const MAX_TICKS: u64 = 72_000;

/// Headwind admission cap [m/s].
pub const MAX_HEADWIND_MPS: f64 = 20.0;

/// Headwind admission FLOOR [m/s]: the aero core's admitted domain
/// starts at 1 m/s airspeed, and the rail start's airspeed IS the
/// headwind — a scenario below this floor would refuse at tick 1
/// inside the coupled solve instead of at admission. 1.5 keeps a
/// margin above the domain edge; a true dead-calm start needs a
/// low-speed rail tier that does not exist yet (declared).

/// Rail length cap [m] (the 1903 rail was 60 ft ≈ 18.3 m).
pub const MAX_RAIL_M: f64 = 60.0;

/// See [`MAX_HEADWIND_MPS`]'s floor note.
pub const MIN_HEADWIND_MPS: f64 = 1.5;

/// Rail rolling-friction coefficient (declared Estimated).
pub const RAIL_MU: f64 = 0.05;

/// Catapult pull-force cap [N] (absurd-input guard).
pub const MAX_CATAPULT_FORCE_N: f64 = 5_000.0;

/// Catapult pull-length cap [m].
pub const MAX_CATAPULT_PULL_M: f64 = 60.0;

/// The 1904 Huffman catapult (registered variant; flyer-reference
/// 'tunable-with-provenance'). Primary lineage (wright-brothers.org
/// reconstruction): ~1400 lb of weights, 16 ft drop, 3:1 pulley SPEED
/// multiplication — pulleys multiply DISTANCE not force, so the tow is
/// ~1400/3 lb ≈ 466 lb ≈ 2073 N over ~48 ft ≈ 14.6 m of pull. The NPS
/// 1600 lb / 30 ft derrick lineage remains a declared variant axis.
/// Constant-force tow is the declared Estimated tier (rope/pulley
/// dynamics unmodeled).
pub const CATAPULT_1904_V1: CatapultSpec = CatapultSpec {
    pull_force_n: 2_073.0,
    pull_length_m: 14.6,
};

/// Catapult launch assist (rail tow until the drop exhausts).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CatapultSpec {
    /// Constant tow force along the rail [N].
    pub pull_force_n: f64,
    /// Pull length along the rail [m] (drop × pulley ratio).
    pub pull_length_m: f64,
}

/// Snapshot payload length v2. The v1 12-float prefix is unchanged;
/// roll and heading are appended so stored v1 records remain decodable.
pub const SNAPSHOT_LEN: usize = 14;

/// Phase remains at its v1 slot; do not infer it from `SNAPSHOT_LEN`.
pub const SNAPSHOT_PHASE_INDEX: usize = 11;

/// Declared reduced rudder moment gain for the 1903 linked-rudder
/// presentation tier [N·m/rad]. It is Estimated, not calibrated.
const WRIGHT_1903_RUDDER_GAIN_NM_PER_RAD: f64 = 66.0;

/// Control mode.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PilotMode {
    /// Controls held at the trim setting (open-loop — diverges, as the
    /// physics says it must).
    FixedControls,
    /// The E4.6c-ii PilotWrightModel flies (registered member).
    Historical(u32),
    /// A human supplies ControlInput every tick.
    Human,
}

/// Human control input for one tick.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ControlInput {
    /// Lever force [N] (+ pull back = nose up).
    pub lever_force_n: f64,
    /// Warp command [rad].
    pub warp_cmd_rad: f64,
}

/// Lifecycle phase.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Phase {
    /// Constrained on the launch rail.
    OnRail,
    /// Free longitudinal flight.
    Airborne,
    /// Run over.
    Ended(TerminalEvent),
}

/// Why the run ended.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TerminalEvent {
    /// Touched the ground (the historical ending).
    GroundContact,
    /// Ran off the rail without lifting.
    RailEndWithoutLift,
    /// Tick budget exhausted (still flying!).
    MaxTicks,
    /// The flight left the certified aero envelope: a physics solve
    /// refused mid-run (receipt kept via `envelope_refusal`).
    EnvelopeExceeded,
    /// A swept-feature strike (blade/hub) whose continuation would
    /// need an unavailable breakage model: the physical run closes at
    /// the LOCALIZED event time (E3.4b; cinematic continuation is
    /// presentation, never physics).
    DamageModelUnavailable,
}

/// The scenario (enters RunIntentId).
#[derive(Clone, Debug, PartialEq)]
pub struct ScenarioSpec {
    /// Study seed.
    pub seed: u64,
    /// Air density [kg/m³].
    pub rho_kg_m3: f64,
    /// Steady headwind [m/s] (Dec-17 class: 9–12).
    pub headwind_mps: f64,
    /// Control mode.
    pub pilot_mode: PilotMode,
    /// Optional assist (HUD-flagged; calibration-isolated).
    pub assist: Option<AssistSpec>,
    /// Optional catapult tow (Huffman launch option; None = rail
    /// headwind alone, the Dec-17 procedure).
    pub catapult: Option<CatapultSpec>,
    /// Rail length [m].
    pub rail_length_m: f64,
    /// Tick budget.
    pub max_ticks: u64,
}

/// The Dec-17 reference scenario.
#[must_use]
pub fn dec17_scenario(seed: u64, pilot_mode: PilotMode) -> ScenarioSpec {
    ScenarioSpec {
        seed,
        rho_kg_m3: 1.294,
        headwind_mps: 11.0,
        pilot_mode,
        assist: None,
        catapult: None,
        rail_length_m: 18.3,
        max_ticks: 2_400, // 20 s
    }
}

/// The Huffman Prairie 1904-05 reference scenario: near-calm wind
/// (the reason the catapult existed), summer-Dayton density at 815 ft
/// (Estimated — no source publishes the density arithmetic; declared),
/// catapult tow, longer rail.
#[must_use]
pub fn huffman_scenario(seed: u64, pilot_mode: PilotMode) -> ScenarioSpec {
    ScenarioSpec {
        seed,
        rho_kg_m3: 1.17,
        headwind_mps: 2.0, // near-calm, above the admission floor
        pilot_mode,
        assist: None,
        catapult: Some(CATAPULT_1904_V1),
        rail_length_m: 30.0,
        max_ticks: 2_400,
    }
}

/// Published per-tick state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SimStateOut {
    /// Tick index.
    pub tick: u64,
    /// Phase.
    pub phase: Phase,
    /// Ground distance along the flight line [m].
    pub x_m: f64,
    /// Height above ground [m].
    pub h_m: f64,
    /// Body airspeed components (u, w) [m/s].
    pub u_mps: f64,
    /// Body vertical airspeed component [m/s].
    pub w_mps: f64,
    /// Pitch rate [rad/s].
    pub q_rad_s: f64,
    /// Pitch attitude [rad].
    pub theta_rad: f64,
    /// Roll rate [rad/s] from the reduced lateral tier.
    pub p_rad_s: f64,
    /// Roll attitude [rad] (+ = right wing down).
    pub phi_rad: f64,
    /// Yaw rate [rad/s] from the reduced lateral tier.
    pub r_rad_s: f64,
    /// Heading [rad] (+ = nose right).
    pub psi_rad: f64,
    /// Canard deflection [rad].
    pub dc_rad: f64,
    /// Warp command [rad].
    pub warp_rad: f64,
    /// Prop speed [rad/s].
    pub omega_prop_rad_s: f64,
    /// Gust vertical velocity [m/s].
    pub gust_w_mps: f64,
    /// Assist active this tick.
    pub assist_active: bool,
    /// Assist canard contribution this tick [rad] (0 when no assist —
    /// envelope-only visibility; NOT part of the frozen ring payload).
    pub assist_dc_rad: f64,
}

/// The engine.
pub struct SimLoop {
    spec: ScenarioSpec,
    design: OpenLoopDesign,
    mech_spec: CanardMechanism,
    perception: PerceptionModelSpec,
    perception_state: PerceptionState,
    pilot: Option<(PilotWrightModel, PilotState)>,
    warp_model: ReducedAeroelasticWarp,
    lateral_model: LateralModel,
    lateral_state: LateralState,
    ou: StationaryOuPath,
    tick0: Tick0State,
    /// Minted AFTER the tick-0 digest (plan law).
    pub run_intent_id: String,
    // Dynamic state.
    tick: u64,
    phase: Phase,
    x_m: f64,
    h_m: f64,
    u: f64,
    w: f64,
    q: f64,
    theta: f64,
    mech: MechState,
    omega: f64,
    warp: f64,
    warm_slip: Option<[f64; 2]>,
    envelope_refusal: Option<Refusal>,
    digest_acc: [u8; 32],
}

impl SimLoop {
    /// Initialize the lifecycle: admission → equilibrate → RunIntentId.
    ///
    /// # Errors
    /// `scenario-invalid` (caps at cap AND cap+1 on headwind, rail,
    /// ticks, rho); equilibration/trim refusals pass through; assist
    /// admission refusals pass through.
    pub fn init(spec: ScenarioSpec) -> Result<SimLoop, Refusal> {
        let ok = spec.rho_kg_m3.is_finite()
            && spec.rho_kg_m3 > 0.5
            && spec.rho_kg_m3 < 2.0
            && spec.headwind_mps.is_finite()
            && (MIN_HEADWIND_MPS..=MAX_HEADWIND_MPS).contains(&spec.headwind_mps)
            && spec.rail_length_m.is_finite()
            && spec.rail_length_m > 0.0
            && spec.rail_length_m <= MAX_RAIL_M
            && spec.max_ticks > 0
            && spec.max_ticks <= MAX_TICKS
            && spec.catapult.is_none_or(|c| {
                c.pull_force_n.is_finite()
                    && c.pull_force_n > 0.0
                    && c.pull_force_n <= MAX_CATAPULT_FORCE_N
                    && c.pull_length_m.is_finite()
                    && c.pull_length_m > 0.0
                    && c.pull_length_m <= MAX_CATAPULT_PULL_M
            });
        if !ok {
            return Err(Refusal {
                code: "scenario-invalid",
                message: format!("{spec:?}"),
                ranked_repairs: vec![format!(
                    "headwind [{MIN_HEADWIND_MPS}, {MAX_HEADWIND_MPS}]; rail (0, {MAX_RAIL_M}]; ticks [1, {MAX_TICKS}]; rho (0.5, 2.0)"
                )],
            });
        }
        if let Some(a) = &spec.assist {
            a.admit()?;
        }
        // Lifecycle coupling schedule: candidate B (the coupling
        // refusal's sanctioned escape). A full rail→airborne run visits
        // off-trim states where candidate A's cap-4 schedule stalls just
        // above tol; B keeps the SAME tolerance with a deeper schedule,
        // and its spec digest binds into every step's evidence.
        let mut design = wright_openloop_v1();
        design.coupling = crate::propcoupling::CANDIDATE_B;
        let mech_spec = CANARD_MECH_V1;
        let ou_modes: Vec<OuMode> = (0..8)
            .map(|i| {
                OuMode::from_correlation_time(0.35 + 0.05 * f64::from(i), 1.2 + 0.5 * f64::from(i))
            })
            .collect::<Result<_, _>>()
            .map_err(map_atmo)?;
        let eq = EquilibrationSpec {
            branch: PrelaunchBranch::HeldOnRailEquilibrated,
            seed: spec.seed,
            ou_modes,
            anchor_tick: STATIONARY_ANCHOR_TICK,
            preroll_ticks: 960,
            rho_kg_m3: spec.rho_kg_m3,
            trim_start: [13.0, 0.06, 0.1, 45.0],
        };
        let tick0 = equilibrate(&design, &mech_spec, &eq)?;
        // Rebuild the OU path at tick 0 (equilibrate consumed its own).
        let mut ou = StationaryOuPath::stationary_at_anchor(spec.seed, eq.ou_modes.clone())
            .map_err(map_atmo)?;
        ou.advance_to(0).map_err(map_atmo)?;
        // RunIntentId: minted AFTER the tick-0 digest, over digest +
        // intent fields (the one place intent may bind).
        let mut p = Vec::new();
        p.extend_from_slice(tick0.digest.as_bytes());
        p.extend_from_slice(&spec.headwind_mps.to_bits().to_le_bytes());
        p.extend_from_slice(&spec.rail_length_m.to_bits().to_le_bytes());
        p.push(match spec.pilot_mode {
            PilotMode::FixedControls => 0,
            PilotMode::Historical(_) => 1,
            PilotMode::Human => 2,
        });
        if let PilotMode::Historical(m) = spec.pilot_mode {
            p.extend_from_slice(&m.to_le_bytes());
        }
        p.push(u8::from(spec.assist.is_some()));
        // Catapult is intent-bearing scenario state.
        p.push(u8::from(spec.catapult.is_some()));
        if let Some(c) = &spec.catapult {
            p.extend_from_slice(&c.pull_force_n.to_bits().to_le_bytes());
            p.extend_from_slice(&c.pull_length_m.to_bits().to_le_bytes());
        }
        let run_intent_id = hash_domain("org.frankensim.wf.run-intent.v1", &p).to_hex();
        let perception = perception_v1(spec.seed);
        let perception_state = perception.init()?;
        let warp_model = ReducedAeroelasticWarp::wright_v1();
        warp_model.admit()?;
        let lateral_model = LateralModel::wright_v1(RudderLinkage::Linked {
            gain_nm_per_rad: WRIGHT_1903_RUDDER_GAIN_NM_PER_RAD,
        });
        let pilot = match spec.pilot_mode {
            PilotMode::Historical(member) => {
                let m = PilotWrightModel::new(member, spec.seed)?;
                let st = m.init()?;
                Some((m, st))
            }
            _ => None,
        };
        // Rail start: aircraft at rest on the rail (airspeed = headwind).
        let theta0 = tick0.trim.alpha_rad;
        let v_air0 = spec.headwind_mps.max(0.1);
        let digest_acc =
            *hash_domain("org.frankensim.wf.sim-digest.v1", tick0.digest.as_bytes()).as_bytes();
        Ok(SimLoop {
            spec,
            design,
            mech_spec,
            perception,
            perception_state,
            pilot,
            warp_model,
            lateral_model,
            lateral_state: LateralState::default(),
            ou,
            mech: tick0.mech,
            omega: tick0.trim.omega_prop_rad_s,
            tick0,
            run_intent_id,
            tick: 0,
            phase: Phase::OnRail,
            x_m: 0.0,
            h_m: 2.4, // rail height class above the sand (visual datum)
            u: v_air0,
            w: 0.0,
            q: 0.0,
            theta: theta0,
            warp: 0.0,
            warm_slip: None,
            envelope_refusal: None,
            digest_acc,
        })
    }

    /// Sweep-harness twin of `init`: replace the Historical pilot's
    /// gains with an explicit candidate. Calibration driver ONLY — the
    /// registered family is the sole identity-bearing path, so nothing
    /// built here may mint a claim.
    ///
    /// # Errors
    /// As `init`; `pilot-gains-invalid`.
    #[doc(hidden)]
    pub fn init_with_pilot_gains(
        spec: ScenarioSpec,
        gains: crate::pilot::PilotGains,
    ) -> Result<SimLoop, Refusal> {
        let mut sim = SimLoop::init(spec)?;
        if sim.pilot.is_some() {
            let m = PilotWrightModel::from_gains(gains, sim.spec.seed)?;
            let st = m.init()?;
            sim.pilot = Some((m, st));
        }
        Ok(sim)
    }

    /// The tick-0 state (already digest-frozen).
    #[must_use]
    pub fn tick0(&self) -> &Tick0State {
        &self.tick0
    }

    /// The envelope-exit receipt, if the run ended `EnvelopeExceeded`.
    #[must_use]
    pub fn envelope_refusal(&self) -> Option<&Refusal> {
        self.envelope_refusal.as_ref()
    }

    /// Serialize the FULL dynamic state as CheckpointStateV2 bytes
    /// (bead guzez.7.1.1, E6.1-i): bit-exact capture with an embedded
    /// blake3 integrity digest. The OU path is counter-addressed and
    /// reconstructs from the tick alone (philox doctrine).
    ///
    /// # Errors
    /// `checkpoint-after-terminal` (a finished run has its replay
    /// envelope; checkpoints are for live runs).
    pub fn save_checkpoint(&self) -> Result<Vec<u8>, Refusal> {
        if matches!(self.phase, Phase::Ended(_)) {
            return Err(Refusal {
                code: "checkpoint-after-terminal",
                message: format!("terminal at tick {}", self.tick),
                ranked_repairs: vec!["export the replay envelope instead".into()],
            });
        }
        let mut b = Vec::new();
        b.extend_from_slice(&2u32.to_le_bytes()); // version
        b.extend_from_slice(&self.tick.to_le_bytes());
        b.push(match self.phase {
            Phase::OnRail => 0,
            Phase::Airborne => 1,
            Phase::Ended(_) => unreachable!("guarded above"),
        });
        for v in [
            self.x_m,
            self.h_m,
            self.u,
            self.w,
            self.q,
            self.theta,
            self.omega,
            self.warp,
            self.mech.delta_rad,
            self.mech.rate_rad_s,
            self.lateral_state.phi_rad,
            self.lateral_state.p_rad_s,
            self.lateral_state.psi_rad,
            self.lateral_state.r_rad_s,
        ] {
            b.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        match self.warm_slip {
            Some(ws) => {
                b.push(1);
                b.extend_from_slice(&ws[0].to_bits().to_le_bytes());
                b.extend_from_slice(&ws[1].to_bits().to_le_bytes());
            }
            None => {
                b.push(0);
                b.extend_from_slice(&0u64.to_le_bytes());
                b.extend_from_slice(&0u64.to_le_bytes());
            }
        }
        b.extend_from_slice(&self.digest_acc);
        let ps = self.perception_state.to_bytes();
        b.extend_from_slice(&(ps.len() as u32).to_le_bytes());
        b.extend_from_slice(&ps);
        match &self.pilot {
            Some((_, st)) => {
                b.push(1);
                let pb = st.to_bytes();
                b.extend_from_slice(&(pb.len() as u32).to_le_bytes());
                b.extend_from_slice(&pb);
            }
            None => b.push(0),
        }
        let digest = hash_domain("org.frankensim.wf.checkpoint-state.v2", &b);
        b.extend_from_slice(digest.as_bytes());
        Ok(b)
    }

    /// Rebuild a live run from a checkpoint (the spec MUST be the
    /// original scenario — init re-derives tick0/RunIntentId, which is
    /// what makes the bit-identity claim checkable downstream).
    ///
    /// # Errors
    /// `checkpoint-too-large` (cap AND cap+1 at 1 MiB);
    /// `checkpoint-tampered` (embedded digest mismatch);
    /// `checkpoint-version-unsupported` (v1 cannot resume exactly
    /// after lateral state entered the digest/state contract);
    /// `checkpoint-malformed`; init refusals pass through.
    pub fn restore_checkpoint(spec: ScenarioSpec, bytes: &[u8]) -> Result<SimLoop, Refusal> {
        const MAX: usize = 1 << 20;
        if bytes.len() > MAX {
            return Err(Refusal {
                code: "checkpoint-too-large",
                message: format!("{} bytes > {MAX}", bytes.len()),
                ranked_repairs: vec!["a checkpoint is a few KiB".into()],
            });
        }
        let bad = |m: &str| Refusal {
            code: "checkpoint-malformed",
            message: m.into(),
            ranked_repairs: vec!["re-export the checkpoint".into()],
        };
        if bytes.len() < 32 {
            return Err(bad("shorter than the digest"));
        }
        let (body, tail) = bytes.split_at(bytes.len() - 32);
        let digest_v1 = hash_domain("org.frankensim.wf.checkpoint-state.v1", body);
        let digest_v2 = hash_domain("org.frankensim.wf.checkpoint-state.v2", body);
        let matches_v1 = digest_v1.as_bytes() == tail;
        let matches_v2 = digest_v2.as_bytes() == tail;
        if !matches_v1 && !matches_v2 {
            return Err(Refusal {
                code: "checkpoint-tampered",
                message: "embedded digest mismatch".into(),
                ranked_repairs: vec!["the checkpoint bytes were altered; refuse".into()],
            });
        }
        let version = body
            .get(..4)
            .and_then(|v| <[u8; 4]>::try_from(v).ok())
            .map(u32::from_le_bytes)
            .ok_or_else(|| bad("missing version"))?;
        match (version, matches_v1, matches_v2) {
            (2, _, true) => {}
            (1, true, _) => {
                return Err(Refusal {
                    code: "checkpoint-version-unsupported",
                    message: "CheckpointStateV1 lacks lateral state and cannot exact-resume v2"
                        .into(),
                    ranked_repairs: vec![
                        "finish the run with the v1 engine or start a new v2 run".into(),
                    ],
                });
            }
            _ => return Err(bad("version and integrity domain disagree")),
        }
        let mut pos = 0usize;
        let take = |pos: &mut usize, n: usize| -> Result<&[u8], Refusal> {
            let s = body.get(*pos..*pos + n).ok_or_else(|| bad("truncated"))?;
            *pos += n;
            Ok(s)
        };
        let parsed_version = u32::from_le_bytes(take(&mut pos, 4)?.try_into().expect("4"));
        debug_assert_eq!(parsed_version, version);
        let tick = u64::from_le_bytes(take(&mut pos, 8)?.try_into().expect("8"));
        let phase_code = take(&mut pos, 1)?[0];
        let mut f = [0.0f64; 14];
        for v in &mut f {
            *v = f64::from_bits(u64::from_le_bytes(
                take(&mut pos, 8)?.try_into().expect("8"),
            ));
        }
        let warm_flag = take(&mut pos, 1)?[0];
        let w0 = f64::from_bits(u64::from_le_bytes(
            take(&mut pos, 8)?.try_into().expect("8"),
        ));
        let w1 = f64::from_bits(u64::from_le_bytes(
            take(&mut pos, 8)?.try_into().expect("8"),
        ));
        let mut digest_acc = [0u8; 32];
        digest_acc.copy_from_slice(take(&mut pos, 32)?);
        let ps_len = u32::from_le_bytes(take(&mut pos, 4)?.try_into().expect("4")) as usize;
        let ps_bytes = take(&mut pos, ps_len)?.to_vec();
        let (perception_state, used) = PerceptionState::from_bytes(&ps_bytes)?;
        if used != ps_len {
            return Err(bad("perception length mismatch"));
        }
        let pilot_flag = take(&mut pos, 1)?[0];
        let pilot_state = if pilot_flag == 1 {
            let pl = u32::from_le_bytes(take(&mut pos, 4)?.try_into().expect("4")) as usize;
            let pb = take(&mut pos, pl)?.to_vec();
            let (st, used) = crate::pilot::PilotState::from_bytes(&pb)?;
            if used != pl {
                return Err(bad("pilot length mismatch"));
            }
            Some(st)
        } else {
            None
        };
        let mut sim = SimLoop::init(spec)?;
        if pilot_state.is_some() != sim.pilot.is_some() {
            return Err(bad("pilot presence disagrees with the scenario"));
        }
        sim.tick = tick;
        sim.phase = match phase_code {
            0 => Phase::OnRail,
            1 => Phase::Airborne,
            _ => return Err(bad("phase code")),
        };
        sim.x_m = f[0];
        sim.h_m = f[1];
        sim.u = f[2];
        sim.w = f[3];
        sim.q = f[4];
        sim.theta = f[5];
        sim.omega = f[6];
        sim.warp = f[7];
        sim.mech = MechState {
            delta_rad: f[8],
            rate_rad_s: f[9],
        };
        sim.lateral_state = LateralState {
            phi_rad: f[10],
            p_rad_s: f[11],
            psi_rad: f[12],
            r_rad_s: f[13],
        };
        sim.warm_slip = (warm_flag == 1).then_some([w0, w1]);
        sim.digest_acc = digest_acc;
        sim.perception_state = perception_state;
        if let (Some(st), Some((_, slot))) = (pilot_state, sim.pilot.as_mut()) {
            *slot = st;
        }
        // The OU path is counter-addressed: advance to the checkpoint
        // tick (the next step advances to tick+1 exactly as the
        // uninterrupted run would).
        sim.ou.advance_to(tick as i64).map_err(map_atmo)?;
        Ok(sim)
    }

    /// One 120 Hz step.
    ///
    /// # Errors
    /// `run-ended` (stepping past the terminal event);
    /// `control-input-missing` (Human mode without input);
    /// physics refusals pass through.
    pub fn step(&mut self, input: Option<ControlInput>) -> Result<SimStateOut, Refusal> {
        if let Phase::Ended(_) = self.phase {
            return Err(Refusal {
                code: "run-ended",
                message: format!("terminal at tick {}", self.tick),
                ranked_repairs: vec!["init a new run".into()],
            });
        }
        let dt = 1.0 / PERCEPTION_HZ;
        self.tick += 1;
        self.ou.advance_to(self.tick as i64).map_err(map_atmo)?;
        // Reduced gust coupling (declared): vertical gust = weighted OU
        // amplitude sum; enters as an angle-of-attack increment.
        let gust_w: f64 = self
            .ou
            .amplitudes()
            .iter()
            .enumerate()
            .map(|(i, a)| a * 0.12 / (1.0 + i as f64))
            .sum();
        // Statement-split (guzez.7.2.1): no fused mul-add.
        let u2 = self.u * self.u;
        let w2 = self.w * self.w;
        let v_air = (u2 + w2).sqrt().max(0.1);
        // Perception + control.
        let wdot_est = 0.0; // heave cue at the reduced tier (declared)
        let cues = self.perception.step(
            &mut self.perception_state,
            pack_cues(
                self.theta - self.tick0.trim.alpha_rad,
                self.q,
                wdot_est,
                self.lateral_state.phi_rad,
                self.lateral_state.p_rad_s,
                self.lateral_state.r_rad_s,
            ),
        )?;
        let (force, warp_cmd) = match self.spec.pilot_mode {
            PilotMode::FixedControls => {
                // Hold the lever at the trim deflection (the settle law).
                let f = (3000.0 * (self.tick0.trim.delta_canard_rad - self.mech.delta_rad)
                    - 180.0 * self.mech.rate_rad_s)
                    / self.mech_spec.lever_gain_nm_per_n;
                (f.clamp(-220.0, 220.0), 0.0)
            }
            PilotMode::Historical(_) => {
                let (model, st) = self.pilot.as_mut().expect("historical mode has a pilot");
                let cmd = model.step(
                    st,
                    &cues,
                    self.mech.delta_rad - self.tick0.trim.delta_canard_rad,
                    self.mech.rate_rad_s,
                    0.0,
                    0.0,
                )?;
                (cmd.lever_force_n, cmd.warp_cmd_rad)
            }
            PilotMode::Human => {
                let Some(inp) = input else {
                    return Err(Refusal {
                        code: "control-input-missing",
                        message: format!("Human mode requires input at tick {}", self.tick),
                        ranked_repairs: vec![
                            "supply ControlInput every tick (ApplyNextEligibleTickAndFlag)".into(),
                        ],
                    });
                };
                if !inp.lever_force_n.is_finite() || !inp.warp_cmd_rad.is_finite() {
                    return Err(Refusal {
                        code: "control-input-missing",
                        message: "non-finite control input".into(),
                        ranked_repairs: vec!["finite lever force and warp".into()],
                    });
                }
                (inp.lever_force_n.clamp(-220.0, 220.0), inp.warp_cmd_rad)
            }
        };
        self.warp = warp_cmd.clamp(-0.148, 0.148);
        // Mechanism (m_aero = 0 stick tier, H-02c convention, declared).
        self.mech = self.mech_spec.step(self.mech, 0.0, force, dt)?.0;
        // Assist (bounded, flagged).
        let mut dc_eff = self.mech.delta_rad;
        let mut assist_active = false;
        let mut assist_dc = 0.0;
        if let Some(a) = &self.spec.assist {
            let out = a.apply(
                self.q,
                self.theta - self.tick0.trim.alpha_rad,
                self.mech_spec.stop_rad,
            )?;
            dc_eff += out.dc_assist_rad;
            assist_active = out.active;
            assist_dc = out.dc_assist_rad;
        }
        // Aerodynamics at the CURRENT state (+ gust alpha increment).
        let alpha = det::atan2(self.w, self.u) + gust_w / v_air;
        let alpha_clamped = alpha.clamp(-0.55, 0.55);
        let b = match self.design.force_buildup_warm(
            v_air,
            alpha_clamped,
            dc_eff.clamp(-0.55, 0.55),
            self.omega.clamp(10.5, 119.5),
            self.q.clamp(-2.9, 2.9),
            self.spec.rho_kg_m3,
            self.warm_slip,
        ) {
            Ok(b) => b,
            Err(refusal) => {
                // The flight left the certified aero envelope (a wild
                // enough state that the coupled/wing solve refuses). The
                // RUN ends with a typed terminal event and the refusal
                // kept as the receipt — never a mid-flight panic handed
                // to the UI, never a silent fallback model.
                self.phase = Phase::Ended(TerminalEvent::EnvelopeExceeded);
                self.envelope_refusal = Some(refusal);
                let out = self.state_out(gust_w, assist_active, assist_dc);
                let mut bytes = self.digest_acc.to_vec();
                for v in self.snapshot_payload(&out) {
                    bytes.extend_from_slice(&v.to_bits().to_le_bytes());
                }
                self.digest_acc =
                    *hash_domain("org.frankensim.wf.sim-digest.v1", &bytes).as_bytes();
                return Ok(out);
            }
        };
        self.warm_slip = Some(b.coupled.w_slip);
        let m = self.design.gross_mass_kg;
        match self.phase {
            Phase::OnRail => {
                // Constrained: normal force carries weight minus lift;
                // accelerate along the rail under thrust − drag −
                // rolling friction. Gravity is inside force_n already —
                // remove its x-component and use the rail-plane balance.
                let lift = b.lift_n;
                let weight = m * 9.80665;
                let normal = (weight - lift).max(0.0);
                let thrust = b.thrust_n[0] + b.thrust_n[1];
                // Catapult tow: constant force while the pull lasts
                // (x_m measures rail distance; the drop exhausts at
                // pull_length — constant-force tier, declared).
                let tow = match &self.spec.catapult {
                    Some(c) if self.x_m < c.pull_length_m => c.pull_force_n,
                    _ => 0.0,
                };
                let du_ground = (thrust + tow - b.drag_n - RAIL_MU * normal) / m;
                let inc_u = dt * du_ground;
                self.u += inc_u;
                let inc_x = dt * (self.u - self.spec.headwind_mps).max(0.0);
                self.x_m += inc_x;
                let inc_om = dt * b.torque_imbalance_nm / 1.6;
                self.omega += inc_om;
                if lift >= weight {
                    self.phase = Phase::Airborne;
                } else if self.x_m > self.spec.rail_length_m {
                    self.phase = Phase::Ended(TerminalEvent::RailEndWithoutLift);
                }
            }
            Phase::Airborne => {
                // FRD longitudinal kinematics.
                let du = b.force_n[0] / m - self.q * self.w;
                let dw = b.force_n[2] / m + self.q * self.u;
                let dq = b.moment_y_nm / IYY_KG_M2;
                let inc_u = dt * du;
                self.u += inc_u;
                let inc_w = dt * dw;
                self.w += inc_w;
                let inc_q = dt * dq;
                self.q += inc_q;
                let inc_th = dt * self.q;
                self.theta += inc_th;
                let inc_om = dt * b.torque_imbalance_nm / 1.6;
                self.omega += inc_om;
                let climb = self.u * det::sin(self.theta) - self.w * det::cos(self.theta);
                let inc_h = dt * climb;
                self.h_m += inc_h;
                let ground_speed = (self.u * det::cos(self.theta) + self.w * det::sin(self.theta)
                    - self.spec.headwind_mps)
                    .max(0.0);
                let inc_x = dt * ground_speed;
                self.x_m += inc_x;
                // Project the achieved antisymmetric strip twists back
                // onto the registered warp basis. Common trim washout
                // is symmetric and therefore contributes no lateral
                // command; the authority-bearing component drives the
                // reduced roll/yaw state.
                let lateral_result = self
                    .warp_model
                    .evaluate(
                        self.warp,
                        0.5 * self.spec.rho_kg_m3 * v_air * v_air,
                        alpha_clamped,
                    )
                    .and_then(|loaded| {
                        let (projection, norm) = loaded
                            .strips
                            .iter()
                            .zip(&self.warp_model.basis)
                            .fold((0.0, 0.0), |(projection, norm), (strip, basis)| {
                                (projection + strip.loaded_rad * basis, norm + basis * basis)
                            });
                        self.lateral_model.step(
                            &mut self.lateral_state,
                            projection / norm,
                            v_air,
                            self.spec.rho_kg_m3,
                            dt,
                        )
                    });
                if let Err(refusal) = lateral_result {
                    self.phase = Phase::Ended(TerminalEvent::EnvelopeExceeded);
                    self.envelope_refusal = Some(refusal);
                }
                if self.h_m <= 0.0 && !matches!(self.phase, Phase::Ended(_)) {
                    self.h_m = 0.0;
                    self.phase = Phase::Ended(TerminalEvent::GroundContact);
                }
            }
            Phase::Ended(_) => unreachable!("guarded above"),
        }
        if self.tick >= self.spec.max_ticks && !matches!(self.phase, Phase::Ended(_)) {
            self.phase = Phase::Ended(TerminalEvent::MaxTicks);
        }
        // Chain the digest.
        let out = self.state_out(gust_w, assist_active, assist_dc);
        let mut bytes = self.digest_acc.to_vec();
        for v in self.snapshot_payload(&out) {
            bytes.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        self.digest_acc = *hash_domain("org.frankensim.wf.sim-digest.v1", &bytes).as_bytes();
        Ok(out)
    }

    fn state_out(&self, gust_w: f64, assist_active: bool, assist_dc_rad: f64) -> SimStateOut {
        SimStateOut {
            tick: self.tick,
            phase: self.phase,
            x_m: self.x_m,
            h_m: self.h_m,
            u_mps: self.u,
            w_mps: self.w,
            q_rad_s: self.q,
            theta_rad: self.theta,
            p_rad_s: self.lateral_state.p_rad_s,
            phi_rad: self.lateral_state.phi_rad,
            r_rad_s: self.lateral_state.r_rad_s,
            psi_rad: self.lateral_state.psi_rad,
            dc_rad: self.mech.delta_rad,
            warp_rad: self.warp,
            omega_prop_rad_s: self.omega,
            gust_w_mps: gust_w,
            assist_active,
            assist_dc_rad,
        }
    }

    /// The frozen v2 ring payload for a state (`SNAPSHOT_LEN` floats).
    /// Slots 0..12 are the exact v1 prefix; phi/psi are appended.
    #[must_use]
    pub fn snapshot_payload(&self, s: &SimStateOut) -> [f64; SNAPSHOT_LEN] {
        [
            s.x_m,
            s.h_m,
            s.u_mps,
            s.w_mps,
            s.q_rad_s,
            s.theta_rad,
            s.dc_rad,
            s.warp_rad,
            s.omega_prop_rad_s,
            s.gust_w_mps,
            f64::from(u8::from(s.assist_active)),
            match s.phase {
                Phase::OnRail => 0.0,
                Phase::Airborne => 1.0,
                Phase::Ended(TerminalEvent::GroundContact) => 2.0,
                Phase::Ended(TerminalEvent::RailEndWithoutLift) => 3.0,
                Phase::Ended(TerminalEvent::MaxTicks) => 4.0,
                Phase::Ended(TerminalEvent::EnvelopeExceeded) => 5.0,
                Phase::Ended(TerminalEvent::DamageModelUnavailable) => 6.0,
            },
            s.phi_rad,
            s.psi_rad,
        ]
    }

    /// Lifecycle digest so far (hex).
    #[must_use]
    pub fn digest_hex(&self) -> String {
        let mut s = String::with_capacity(64);
        for b in self.digest_acc {
            use core::fmt::Write;
            let _ = write!(s, "{b:02x}");
        }
        s
    }

    /// Current phase.
    #[must_use]
    pub fn phase(&self) -> Phase {
        self.phase
    }
}

fn map_atmo(e: fs_atmo::Refusal) -> Refusal {
    Refusal {
        code: e.code,
        message: e.message,
        ranked_repairs: e.ranked_repairs,
    }
}
