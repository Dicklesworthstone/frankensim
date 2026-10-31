//! ReplayEnvelopeV1 + scrub-to-tick (bead wf-root-guzez.7.1.2,
//! E6.1-ii). The E0.9-frozen replay artifact for a finished run:
//! scenario intent, tick-0 digest, RunIntentId, a checkpoint index
//! every K ticks (E6.1-i CheckpointStateV1 bytes), the event-tick
//! index for the scrubber (liftoff / undulation crossings / terminal),
//! the terminal phase and final chained digest — all under one
//! embedded blake3 integrity digest, parsed fail-closed.
//!
//! Scrub law: `scrub_to_tick(env, t)` restores the NEAREST checkpoint
//! at or below t and re-marches the remainder — the battery proves the
//! result equals uninterrupted execution BITWISE at t.
//!
//! Explicit spatial-wind recorders retain the complete V3/V4 tick-zero
//! checkpoint, so scrubbing before the first periodic checkpoint preserves
//! the original wind model and sampling mode. Legacy recording bytes are
//! unchanged. Initial and restored identities are checked before remarching.
//!
//! v1 records deterministic pilot modes (FixedControls / Historical);
//! a Human run needs its applied input trace bound in — that is the
//! 7.1.3 A/B lane, and recording a Human run here REFUSES rather than
//! silently dropping the inputs.

use crate::Refusal;
use crate::assist::ASSIST_V1;
use crate::atmosphere::MaterialAirProfile;
use crate::simloop::{
    CATAPULT_1904_V1, CatapultSpec, Phase, PilotMode, ScenarioSpec, SimLoop, SimStateOut,
    SpatialWindSpec, TerminalEvent,
};
use fs_blake3::hash_domain;
use fs_material::gas::GasState;

/// Envelope schema id (domain-separated digest).
pub const REPLAY_SCHEMA: &str = "org.frankensim.wf.replay-envelope.v1";

/// Checkpoint-interval cap [ticks].
pub const MAX_CHECKPOINT_INTERVAL: u64 = 7_200;

/// Envelope size cap.
pub const MAX_ENVELOPE_BYTES: usize = 64 << 20;

/// The scrubber's event index.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EventIndex {
    /// First airborne tick (None = never lifted).
    pub liftoff_tick: Option<u64>,
    /// Ticks where the airborne pitch rate flipped sign (undulation
    /// crossings — pairs make porpoise cycles).
    pub undulation_ticks: Vec<u64>,
    /// Terminal tick.
    pub terminal_tick: u64,
    /// Terminal phase code (the snapshot payload convention).
    pub terminal_code: u8,
}

/// The in-memory envelope.
#[derive(Clone, Debug, PartialEq)]
pub struct ReplayEnvelope {
    /// Scenario intent (reconstructable exactly).
    pub spec: ScenarioSpec,
    /// Tick-0 digest (identity anchor).
    pub tick0_digest: String,
    /// RunIntentId.
    pub run_intent_id: String,
    /// Checkpoint interval K.
    pub interval: u64,
    /// (tick, versioned checkpoint bytes) in ascending tick order. Explicit
    /// spatial-wind recordings always begin with their actual tick-zero state.
    pub checkpoints: Vec<(u64, Vec<u8>)>,
    /// Scrubber events.
    pub events: EventIndex,
    /// Final chained digest (hex).
    pub final_digest: String,
}

fn phase_code(phase: Phase) -> u8 {
    match phase {
        Phase::OnRail => 0,
        Phase::Airborne => 1,
        Phase::Ended(TerminalEvent::GroundContact) => 2,
        Phase::Ended(TerminalEvent::RailEndWithoutLift) => 3,
        Phase::Ended(TerminalEvent::MaxTicks) => 4,
        Phase::Ended(TerminalEvent::EnvelopeExceeded) => 5,
        Phase::Ended(TerminalEvent::DamageModelUnavailable) => 6,
    }
}

/// Record a full deterministic run into an envelope.
///
/// # Errors
/// `replay-human-needs-trace` (v1 law); `replay-interval-invalid`
/// (cap AND cap+1); init/step refusals pass through.
pub fn record_replay(spec: ScenarioSpec, interval: u64) -> Result<ReplayEnvelope, Refusal> {
    record_initialized(spec, interval, SimLoop::init, false)
}

/// Record the explicit CG-sampled spatial-wind lifecycle. Its tick-zero V3
/// checkpoint carries the complete declared field and sequential OU state.
///
/// # Errors
/// As [`record_replay`] and [`SimLoop::init_with_spatial_wind`].
pub fn record_replay_with_spatial_wind(
    spec: ScenarioSpec,
    wind: SpatialWindSpec,
    interval: u64,
) -> Result<ReplayEnvelope, Refusal> {
    record_initialized(
        spec,
        interval,
        |scenario| SimLoop::init_with_spatial_wind(scenario, wind),
        true,
    )
}

/// Record distributed wing/canard/rotor wind sampling. The initial V4
/// checkpoint preserves this mode even when the run ends before checkpoint K.
///
/// # Errors
/// As [`record_replay`] and [`SimLoop::init_with_distributed_spatial_wind`].
pub fn record_replay_with_distributed_spatial_wind(
    spec: ScenarioSpec,
    wind: SpatialWindSpec,
    interval: u64,
) -> Result<ReplayEnvelope, Refusal> {
    record_initialized(
        spec,
        interval,
        |scenario| SimLoop::init_with_distributed_spatial_wind(scenario, wind),
        true,
    )
}

/// Record the complete material gas with CG-sampled wind. The retained scenario
/// uses material-derived density and the initial V5 checkpoint preserves every
/// gas coefficient, so scrub never reintroduces the authored scalar density.
pub fn record_replay_with_material_spatial_wind(
    spec: ScenarioSpec,
    wind: SpatialWindSpec,
    gas: GasState,
    interval: u64,
) -> Result<ReplayEnvelope, Refusal> {
    record_initialized(
        spec,
        interval,
        |scenario| SimLoop::init_with_material_spatial_wind(scenario, wind, gas),
        true,
    )
}

/// Record the same complete material gas with distributed strip/disk inflow.
pub fn record_replay_with_material_distributed_spatial_wind(
    spec: ScenarioSpec,
    wind: SpatialWindSpec,
    gas: GasState,
    interval: u64,
) -> Result<ReplayEnvelope, Refusal> {
    record_initialized(
        spec,
        interval,
        |scenario| SimLoop::init_with_material_distributed_spatial_wind(scenario, wind, gas),
        true,
    )
}

/// Record an altitude-dependent material atmosphere on the CG wind lifecycle.
/// V6 checkpoints retain the profile definition and ground datum, so resumed
/// forces continue sampling hydrostatic pressure/density at evolving heights.
pub fn record_replay_with_profile_spatial_wind(
    spec: ScenarioSpec,
    wind: SpatialWindSpec,
    profile: MaterialAirProfile,
    interval: u64,
) -> Result<ReplayEnvelope, Refusal> {
    record_initialized(
        spec,
        interval,
        |scenario| SimLoop::init_with_profile_spatial_wind(scenario, wind, profile),
        true,
    )
}

/// Record the same material profile with distributed strip/disk wind sampling.
pub fn record_replay_with_profile_distributed_spatial_wind(
    spec: ScenarioSpec,
    wind: SpatialWindSpec,
    profile: MaterialAirProfile,
    interval: u64,
) -> Result<ReplayEnvelope, Refusal> {
    record_initialized(
        spec,
        interval,
        |scenario| SimLoop::init_with_profile_distributed_spatial_wind(scenario, wind, profile),
        true,
    )
}

fn record_initialized(
    spec: ScenarioSpec,
    interval: u64,
    initialize: impl FnOnce(ScenarioSpec) -> Result<SimLoop, Refusal>,
    retain_initial: bool,
) -> Result<ReplayEnvelope, Refusal> {
    if matches!(spec.pilot_mode, PilotMode::Human) {
        return Err(Refusal {
            code: "replay-human-needs-trace",
            message: "a Human run replays through its applied input trace (E6.1-iii)".into(),
            ranked_repairs: vec!["record with the A/B lane; never drop inputs silently".into()],
        });
    }
    if interval == 0 || interval > MAX_CHECKPOINT_INTERVAL {
        return Err(Refusal {
            code: "replay-interval-invalid",
            message: format!("interval {interval} outside [1, {MAX_CHECKPOINT_INTERVAL}]"),
            ranked_repairs: vec!["a few hundred ticks is the scrubber class".into()],
        });
    }
    let mut sim = initialize(spec.clone())?;
    let spec = sim.scenario().clone();
    let tick0_digest = sim.tick0().digest.clone();
    let run_intent_id = sim.run_intent_id.clone();
    let mut checkpoints = Vec::new();
    if retain_initial {
        checkpoints.push((0, sim.save_checkpoint()?));
    }
    let mut liftoff = None;
    let mut undulation_ticks = Vec::new();
    let mut last_sign = 0i8;
    let (terminal_tick, terminal_code) = loop {
        let out = sim.step(None)?;
        if matches!(out.phase, Phase::Airborne) {
            if liftoff.is_none() {
                liftoff = Some(out.tick);
            }
            let s = if out.q_rad_s > 1e-3 {
                1i8
            } else if out.q_rad_s < -1e-3 {
                -1i8
            } else {
                0
            };
            if s != 0 && last_sign != 0 && s != last_sign {
                undulation_ticks.push(out.tick);
            }
            if s != 0 {
                last_sign = s;
            }
        }
        if let Phase::Ended(_) = out.phase {
            break (out.tick, phase_code(out.phase));
        }
        if out.tick % interval == 0 {
            checkpoints.push((out.tick, sim.save_checkpoint()?));
        }
    };
    Ok(ReplayEnvelope {
        spec,
        tick0_digest,
        run_intent_id,
        interval,
        checkpoints,
        events: EventIndex {
            liftoff_tick: liftoff,
            undulation_ticks,
            terminal_tick,
            terminal_code,
        },
        final_digest: sim.digest_hex(),
    })
}

impl ReplayEnvelope {
    /// Canonical bytes with the embedded integrity digest.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut b = Vec::new();
        b.extend_from_slice(&1u32.to_le_bytes());
        b.extend_from_slice(&self.spec.seed.to_le_bytes());
        for v in [
            self.spec.rho_kg_m3,
            self.spec.headwind_mps,
            self.spec.rail_length_m,
        ] {
            b.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        b.extend_from_slice(&self.spec.max_ticks.to_le_bytes());
        let (mode, member) = match self.spec.pilot_mode {
            PilotMode::FixedControls => (0u8, 0u32),
            PilotMode::Historical(m) => (1, m),
            PilotMode::Human => (2, 0),
        };
        b.push(mode);
        b.extend_from_slice(&member.to_le_bytes());
        b.push(u8::from(self.spec.assist.is_some()));
        match &self.spec.catapult {
            Some(c) => {
                b.push(1);
                b.extend_from_slice(&c.pull_force_n.to_bits().to_le_bytes());
                b.extend_from_slice(&c.pull_length_m.to_bits().to_le_bytes());
            }
            None => {
                b.push(0);
                b.extend_from_slice(&[0u8; 16]);
            }
        }
        b.extend_from_slice(self.tick0_digest.as_bytes());
        b.extend_from_slice(self.run_intent_id.as_bytes());
        b.extend_from_slice(&self.interval.to_le_bytes());
        b.extend_from_slice(&(self.checkpoints.len() as u32).to_le_bytes());
        for (tick, bytes) in &self.checkpoints {
            b.extend_from_slice(&tick.to_le_bytes());
            b.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
            b.extend_from_slice(bytes);
        }
        match self.events.liftoff_tick {
            Some(t) => {
                b.push(1);
                b.extend_from_slice(&t.to_le_bytes());
            }
            None => {
                b.push(0);
                b.extend_from_slice(&0u64.to_le_bytes());
            }
        }
        b.extend_from_slice(&(self.events.undulation_ticks.len() as u32).to_le_bytes());
        for t in &self.events.undulation_ticks {
            b.extend_from_slice(&t.to_le_bytes());
        }
        b.extend_from_slice(&self.events.terminal_tick.to_le_bytes());
        b.push(self.events.terminal_code);
        b.extend_from_slice(self.final_digest.as_bytes());
        let digest = hash_domain(REPLAY_SCHEMA, &b);
        b.extend_from_slice(digest.as_bytes());
        b
    }

    /// Fail-closed parse (hostile-twin hardened; assist/catapult
    /// reconstruct to the REGISTERED constants — the intent flags are
    /// what the envelope carries).
    ///
    /// # Errors
    /// `replay-too-large`; `replay-tampered`; `replay-malformed`.
    pub fn from_bytes(bytes: &[u8]) -> Result<ReplayEnvelope, Refusal> {
        if bytes.len() > MAX_ENVELOPE_BYTES {
            return Err(Refusal {
                code: "replay-too-large",
                message: format!("{} bytes", bytes.len()),
                ranked_repairs: vec!["64 MiB cap".into()],
            });
        }
        let bad = |m: &str| Refusal {
            code: "replay-malformed",
            message: m.into(),
            ranked_repairs: vec!["re-export the envelope".into()],
        };
        if bytes.len() < 32 {
            return Err(bad("shorter than the digest"));
        }
        let (body, tail) = bytes.split_at(bytes.len() - 32);
        if hash_domain(REPLAY_SCHEMA, body).as_bytes() != tail {
            return Err(Refusal {
                code: "replay-tampered",
                message: "embedded digest mismatch".into(),
                ranked_repairs: vec!["the envelope bytes were altered; refuse".into()],
            });
        }
        let mut pos = 0usize;
        let take = |pos: &mut usize, n: usize| -> Result<&[u8], Refusal> {
            let s = body.get(*pos..*pos + n).ok_or_else(|| bad("truncated"))?;
            *pos += n;
            Ok(s)
        };
        macro_rules! u32le {
            () => {
                u32::from_le_bytes(take(&mut pos, 4)?.try_into().expect("4"))
            };
        }
        macro_rules! u64le {
            () => {
                u64::from_le_bytes(take(&mut pos, 8)?.try_into().expect("8"))
            };
        }
        macro_rules! f64le {
            () => {
                f64::from_bits(u64::from_le_bytes(
                    take(&mut pos, 8)?.try_into().expect("8"),
                ))
            };
        }
        if u32le!() != 1 {
            return Err(bad("unknown version"));
        }
        let seed = u64le!();
        let rho_kg_m3 = f64le!();
        let headwind_mps = f64le!();
        let rail_length_m = f64le!();
        let max_ticks = u64le!();
        let mode = take(&mut pos, 1)?[0];
        let member = u32le!();
        let pilot_mode = match mode {
            0 => PilotMode::FixedControls,
            1 => PilotMode::Historical(member),
            2 => PilotMode::Human,
            _ => return Err(bad("mode code")),
        };
        let assist_flag = take(&mut pos, 1)?[0];
        let cat_flag = take(&mut pos, 1)?[0];
        let cat_force = f64le!();
        let cat_len = f64le!();
        let catapult = (cat_flag == 1).then_some(CatapultSpec {
            pull_force_n: cat_force,
            pull_length_m: cat_len,
        });
        let tick0_digest = String::from_utf8(take(&mut pos, 64)?.to_vec())
            .map_err(|_| bad("tick0 digest not utf8"))?;
        let run_intent_id = String::from_utf8(take(&mut pos, 64)?.to_vec())
            .map_err(|_| bad("intent id not utf8"))?;
        let interval = u64le!();
        let n_ckpts = u32le!() as usize;
        if n_ckpts > 100_000 {
            return Err(bad("checkpoint count"));
        }
        let mut checkpoints = Vec::with_capacity(n_ckpts);
        let mut last_tick = None;
        for _ in 0..n_ckpts {
            let tick = u64le!();
            if last_tick.is_some_and(|previous| tick <= previous) {
                return Err(bad("checkpoint ticks not ascending"));
            }
            last_tick = Some(tick);
            let len = u32le!() as usize;
            if len > (1 << 20) {
                return Err(bad("checkpoint length"));
            }
            checkpoints.push((tick, take(&mut pos, len)?.to_vec()));
        }
        let lo_flag = take(&mut pos, 1)?[0];
        let lo_tick = u64le!();
        let n_und = u32le!() as usize;
        if n_und > 1_000_000 {
            return Err(bad("undulation count"));
        }
        let mut undulation_ticks = Vec::with_capacity(n_und);
        for _ in 0..n_und {
            undulation_ticks.push(u64le!());
        }
        let terminal_tick = u64le!();
        let terminal_code = take(&mut pos, 1)?[0];
        if !(2..=6).contains(&terminal_code) {
            return Err(bad("terminal code"));
        }
        let final_digest = String::from_utf8(take(&mut pos, 64)?.to_vec())
            .map_err(|_| bad("final digest not utf8"))?;
        if pos != body.len() {
            return Err(bad("trailing bytes"));
        }
        Ok(ReplayEnvelope {
            spec: ScenarioSpec {
                seed,
                rho_kg_m3,
                headwind_mps,
                pilot_mode,
                assist: (assist_flag == 1).then_some(ASSIST_V1),
                catapult,
                rail_length_m,
                max_ticks,
            },
            tick0_digest,
            run_intent_id,
            interval,
            checkpoints,
            events: EventIndex {
                liftoff_tick: (lo_flag == 1).then_some(lo_tick),
                undulation_ticks,
                terminal_tick,
                terminal_code,
            },
            final_digest,
        })
    }

    /// Scrub: reconstruct the state AT tick `t` from the nearest
    /// checkpoint strictly before t, re-marching the remainder. Tick zero
    /// returns the actual initial state without advancing physics or wind.
    /// Spatial recordings must retain their original initial checkpoint;
    /// an absent/downgraded initial state never selects legacy physics.
    ///
    /// # Errors
    /// `replay-scrub-out-of-range` (t past the terminal),
    /// `replay-initial-state-missing`, `replay-identity-mismatch`,
    /// `replay-checkpoint-index`; restore/step refusals pass through.
    pub fn scrub_to_tick(&self, t: u64) -> Result<SimStateOut, Refusal> {
        if t > self.events.terminal_tick {
            return Err(Refusal {
                code: "replay-scrub-out-of-range",
                message: format!("tick {t} outside [0, {}]", self.events.terminal_tick),
                ranked_repairs: vec!["scrub inside the recorded run".into()],
            });
        }
        let missing = || Refusal {
            code: "replay-initial-state-missing",
            message: "the recorded lifecycle requires its original tick-zero checkpoint".into(),
            ranked_repairs: vec![
                "retain the initial spatial-wind checkpoint; never replay it as legacy OU".into(),
            ],
        };
        let check_identity = |sim: &SimLoop| -> Result<(), Refusal> {
            if sim.tick0().digest != self.tick0_digest || sim.run_intent_id != self.run_intent_id {
                return Err(Refusal {
                    code: "replay-identity-mismatch",
                    message: "checkpoint/initialization belongs to a different recorded lifecycle"
                        .into(),
                    ranked_repairs: vec![
                        "use the original wind mode, field parameters and scenario".into(),
                    ],
                });
            }
            Ok(())
        };
        let initial = self.checkpoints.first().filter(|(tick, _)| *tick == 0);
        // A step envelope at t>0 comes from stepping t itself. The exact-zero
        // request instead selects checkpoint zero explicitly.
        let base = self
            .checkpoints
            .iter()
            .rev()
            .find(|(ct, _)| *ct < t || (t == 0 && *ct == 0));
        let start = base.map_or(0, |(ct, _)| *ct);
        let mut sim = match base {
            Some((_, bytes)) => SimLoop::restore_checkpoint(self.spec.clone(), bytes)?,
            None => {
                let sim = SimLoop::init(self.spec.clone())?;
                // Old envelopes have no explicit initial checkpoint. Only a
                // matching legacy initialization is a valid fallback.
                check_identity(&sim).map_err(|_| missing())?;
                sim
            }
        };
        check_identity(&sim)?;
        if sim.spatial_ground_velocity_mps().is_some() {
            let (_, bytes) = initial.ok_or_else(missing)?;
            if start != 0 {
                // A later checkpoint must not conceal a missing or downgraded
                // initial field. Reuse the source owner's admission/decoder.
                let initial_sim = SimLoop::restore_checkpoint(self.spec.clone(), bytes)?;
                check_identity(&initial_sim)?;
                initial_sim.replay_initial_state()?;
            }
        }
        if t == 0 {
            return sim.replay_initial_state();
        }
        let mut last: Option<SimStateOut> = None;
        for previous in start..t {
            let out = sim.step(None)?;
            if out.tick != previous + 1 {
                return Err(Refusal {
                    code: "replay-checkpoint-index",
                    message: "checkpoint's actual tick disagrees with the envelope index".into(),
                    ranked_repairs: vec![
                        "retain checkpoint ticks with their original state bytes".into(),
                    ],
                });
            }
            last = Some(out);
        }
        last.ok_or_else(|| Refusal {
            code: "replay-scrub-out-of-range",
            message: "scrub produced no step".into(),
            ranked_repairs: vec!["tick equals a checkpoint tick? scrub one ahead".into()],
        })
    }
}
