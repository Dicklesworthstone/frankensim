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
use fs_flyer::assist::ASSIST_V1;
use fs_flyer::simloop::CATAPULT_1904_V1;
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
        Phase::Ended(TerminalEvent::DamageModelUnavailable) => "ended:damage-model-unavailable",
    }
}

fn state_envelope(s: &SimStateOut, envelope_code: Option<&str>) -> String {
    let envelope_field = envelope_code.map_or(String::new(), |c| {
        format!(",\"envelope_refusal_code\":\"{}\"", json_escape(c))
    });
    format!(
        "{{\"ok\":{{\"tick\":{},\"phase\":\"{}\",\"x_m\":{},\"h_m\":{},\"u_mps\":{},\"w_mps\":{},\"q_rad_s\":{},\"theta_rad\":{},\"dc_rad\":{},\"warp_rad\":{},\"omega_prop_rad_s\":{},\"gust_w_mps\":{},\"assist_active\":{},\"assist_dc_rad\":{}{}}}}}",
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
        s.assist_dc_rad,
        envelope_field,
    )
}

/// The 1-second canonical self-test scenario goldens, PER LANE.
/// Cross-lane identity is tracked at bead guzez.7.2.1 (FMA-contraction
/// class, measured 2026-08-21); until it lands each lane pins its own
/// golden and the CI records the cross-lane pair as EXPECTED-DIVERGENT
/// — loudly, never silently.
#[cfg(target_arch = "wasm32")]
pub const SELFTEST_GOLDEN: &str =
    "f088689ae4c60ec33a2034ec7020c85772bfc016968fa9ae5f6d92a308fcbbb6";
/// Native (aarch64/x86 pending the six-lane matrix) canonical golden.
#[cfg(not(target_arch = "wasm32"))]
pub const SELFTEST_GOLDEN: &str =
    "823d9f59dd162c8bc0764e144236d2f00abc48a12142095688a22e59ae95ca9d";

/// Run the canonical 1-second scenario and compare against a golden
/// (the startup self-test core; the shipped entry uses
/// [`SELFTEST_GOLDEN`], the battery's falsifier passes a wrong one).
#[must_use]
pub fn selftest_against(golden: &str) -> String {
    let mut slot = EngineSlot::default();
    let init = slot.init(1903, 1.294, 11.0, MODE_FIXED, 0, 18.3, 120, false, false);
    if !init.starts_with("{\"ok\"") {
        return refusal_envelope(&Refusal {
            code: "determinism-selftest-failed",
            message: format!("canonical init refused: {init}"),
            ranked_repairs: vec!["the build is broken; do not trust results".into()],
        });
    }
    for _ in 0..120 {
        let s = slot.step(false, 0.0, 0.0);
        if !s.starts_with("{\"ok\"") {
            return refusal_envelope(&Refusal {
                code: "determinism-selftest-failed",
                message: format!("canonical step refused: {s}"),
                ranked_repairs: vec!["the build is broken; do not trust results".into()],
            });
        }
    }
    let d = slot.digest();
    let key = "\"digest\":\"";
    let digest = d
        .find(key)
        .map(|i| &d[i + key.len()..i + key.len() + 64])
        .unwrap_or("");
    if digest == golden {
        format!("{{\"ok\":{{\"digest\":\"{digest}\",\"matched\":true}}}}")
    } else {
        refusal_envelope(&Refusal {
            code: "determinism-selftest-failed",
            message: format!("canonical digest {digest} != golden {golden}"),
            ranked_repairs: vec![
                "show the determinism-failure badge; results are untrusted".into(),
                "a perturbed build or platform drift moved the physics".into(),
            ],
        })
    }
}

/// The shipped startup self-test (per-lane golden).
#[must_use]
pub fn selftest() -> String {
    selftest_against(SELFTEST_GOLDEN)
}

impl EngineSlot {
    /// Initialize a lifecycle (replaces any prior run in this slot).
    /// Returns the init envelope: run intent id, tick-0 digest, trim.
    #[allow(clippy::too_many_arguments, clippy::fn_params_excessive_bools)]
    pub fn init(
        &mut self,
        seed: u64,
        rho_kg_m3: f64,
        headwind_mps: f64,
        mode: u32,
        member: u32,
        rail_length_m: f64,
        max_ticks: u64,
        assist: bool,
        catapult: bool,
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
            // E5.3c: the ratified bounded assist (authority 0.3 of the
            // canard stop, HUD-flagged; historical parameter identity
            // is untouched — assist is scenario intent, not model).
            assist: assist.then_some(ASSIST_V1),
            // E5.4: the registered 1904 catapult tow (rail assist;
            // near-calm Huffman remains outside the coupled core's
            // admitted envelope — the run ends with a receipted
            // envelope exit, tracked at guzez.5.7.1).
            catapult: catapult.then_some(CATAPULT_1904_V1),
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
