//! E5.1 engine surface (bead wf-root-guzez.6.2): the REAL fs-flyer
//! lifecycle engine behind a scalar-argument, JSON-envelope API — the
//! shape the wasm exports wrap 1:1. No JSON parsing at the boundary
//! (scalars in, envelopes out), no serde, no mocks: `EngineSlot` holds
//! an actual `fs_flyer::simloop::SimLoop`.
//!
//! Envelope contract (frozen v1):
//!   success  -> {"ok":{...}}
//!   refusal  -> {"refusal":{"code","message","ranked_repairs"}}
//!
//! Documented refusal codes at this surface:
//!   engine-not-initialized  step/state/digest before a successful init
//!   mode-invalid            mode word outside {0,1,2}
//!   scenario-invalid        scenario caps (from the native engine)
//!   control-input-missing   Human mode without (or with non-finite) input
//!   run-ended               stepping past the terminal event
//!   + native physics refusals pass through verbatim (equilibration/trim
//!     at init; the engine converts mid-flight aero refusals into the
//!     EnvelopeExceeded terminal instead of erroring, by design).

use crate::{Refusal, json_escape, refusal_envelope};
use fs_flyer::simloop::{
    ControlInput, Phase, PilotMode, ScenarioSpec, SimLoop, SimStateOut, TerminalEvent,
};

/// Pilot-mode words at the ABI (frozen v1).
pub const MODE_FIXED: u32 = 0;
/// Historical pilot (member selects the registered family member).
pub const MODE_HISTORICAL: u32 = 1;
/// Human control (input required every tick).
pub const MODE_HUMAN: u32 = 2;

/// One engine slot (a worker owns exactly one lifecycle at a time; a
/// new init replaces the old run — the E5.0 ring's epoch bump is the
/// consumer-side guard).
#[derive(Default)]
pub struct EngineSlot {
    sim: Option<SimLoop>,
}

fn map_refusal(e: fs_flyer::Refusal) -> Refusal {
    Refusal {
        code: e.code,
        message: e.message,
        ranked_repairs: e.ranked_repairs,
    }
}

fn not_initialized() -> String {
    refusal_envelope(&Refusal {
        code: "engine-not-initialized",
        message: "call init before step/state/digest".into(),
        ranked_repairs: vec!["engine_init(...) first".into()],
    })
}

fn phase_word(phase: Phase) -> &'static str {
    match phase {
        Phase::OnRail => "on-rail",
        Phase::Airborne => "airborne",
        Phase::Ended(TerminalEvent::GroundContact) => "ended:ground-contact",
        Phase::Ended(TerminalEvent::RailEndWithoutLift) => "ended:rail-end-without-lift",
        Phase::Ended(TerminalEvent::MaxTicks) => "ended:max-ticks",
        Phase::Ended(TerminalEvent::EnvelopeExceeded) => "ended:envelope-exceeded",
    }
}

fn state_envelope(s: &SimStateOut, envelope_code: Option<&str>) -> String {
    let envelope_field = envelope_code.map_or(String::new(), |c| {
        format!(",\"envelope_refusal_code\":\"{}\"", json_escape(c))
    });
    format!(
        "{{\"ok\":{{\"tick\":{},\"phase\":\"{}\",\"x_m\":{},\"h_m\":{},\"u_mps\":{},\"w_mps\":{},\"q_rad_s\":{},\"theta_rad\":{},\"dc_rad\":{},\"warp_rad\":{},\"omega_prop_rad_s\":{},\"gust_w_mps\":{},\"assist_active\":{}{}}}}}",
        s.tick,
        phase_word(s.phase),
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
        s.assist_active,
        envelope_field,
    )
}

impl EngineSlot {
    /// Initialize a lifecycle (replaces any prior run in this slot).
    /// Returns the init envelope: run intent id, tick-0 digest, trim.
    pub fn init(
        &mut self,
        seed: u64,
        rho_kg_m3: f64,
        headwind_mps: f64,
        mode: u32,
        member: u32,
        rail_length_m: f64,
        max_ticks: u64,
    ) -> String {
        let pilot_mode = match mode {
            MODE_FIXED => PilotMode::FixedControls,
            MODE_HISTORICAL => PilotMode::Historical(member),
            MODE_HUMAN => PilotMode::Human,
            other => {
                return refusal_envelope(&Refusal {
                    code: "mode-invalid",
                    message: format!("mode word {other} not in {{0,1,2}}"),
                    ranked_repairs: vec!["0=fixed, 1=historical(member), 2=human".into()],
                });
            }
        };
        let spec = ScenarioSpec {
            seed,
            rho_kg_m3,
            headwind_mps,
            pilot_mode,
            assist: None,
            rail_length_m,
            max_ticks,
        };
        match SimLoop::init(spec) {
            Ok(sim) => {
                let out = format!(
                    "{{\"ok\":{{\"run_intent_id\":\"{}\",\"tick0_digest\":\"{}\",\"trim_v_mps\":{},\"trim_alpha_rad\":{},\"trim_dc_rad\":{},\"trim_omega_rad_s\":{}}}}}",
                    json_escape(&sim.run_intent_id),
                    json_escape(&sim.tick0().digest),
                    sim.tick0().trim.v_mps,
                    sim.tick0().trim.alpha_rad,
                    sim.tick0().trim.delta_canard_rad,
                    sim.tick0().trim.omega_prop_rad_s,
                );
                self.sim = Some(sim);
                out
            }
            Err(e) => refusal_envelope(&map_refusal(e)),
        }
    }

    /// One 120 Hz step. `has_input` gates whether (lever, warp) is a
    /// `ControlInput` (Human mode requires it every tick).
    pub fn step(&mut self, has_input: bool, lever_force_n: f64, warp_cmd_rad: f64) -> String {
        let Some(sim) = self.sim.as_mut() else {
            return not_initialized();
        };
        let input = has_input.then_some(ControlInput {
            lever_force_n,
            warp_cmd_rad,
        });
        match sim.step(input) {
            Ok(out) => {
                let code = matches!(out.phase, Phase::Ended(TerminalEvent::EnvelopeExceeded))
                    .then(|| sim.envelope_refusal().map_or("unknown", |r| r.code));
                state_envelope(&out, code)
            }
            Err(e) => refusal_envelope(&map_refusal(e)),
        }
    }

    /// The chained lifecycle digest (hex).
    pub fn digest(&self) -> String {
        self.sim.as_ref().map_or_else(not_initialized, |sim| {
            format!("{{\"ok\":{{\"digest\":\"{}\"}}}}", sim.digest_hex())
        })
    }
}
