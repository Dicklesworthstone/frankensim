//! E6.1-iii battery (bead wf-root-guzez.7.1.3): identical twins give
//! the FULL-LENGTH common prefix and no divergence; a perturbed-axis
//! SameInputTrace pair diverges at a receipted tick; a
//! terminal-divergence pair records DIFFERENT TerminalEvent kinds; a
//! HumanRefly pair with differing traces diverges where the traces do;
//! trace/mode mismatches refuse at cap AND cap+1; receipts are
//! deterministic (golden).
//! Repro: cargo test -p fs-flyer --test abcompare_battery --release

use fs_flyer::abcompare::{AB_SCHEMA, MAX_TRACE_LEN, ab_compare};
use fs_flyer::simloop::{ControlInput, PilotMode, dec17_scenario};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-e61iii\",\"case\":\"{case}\",{payload}}}");
}

fn neutral_trace(n: usize) -> Vec<ControlInput> {
    vec![
        ControlInput {
            lever_force_n: 0.0,
            warp_cmd_rad: 0.0,
        };
        n
    ]
}

#[test]
fn identical_twins_share_the_full_prefix() {
    let mut spec = dec17_scenario(21, PilotMode::FixedControls);
    spec.max_ticks = 300; // rail-only, cheap
    let r = ab_compare("same-input-trace", &spec, None, &spec, None).unwrap();
    assert_eq!(r.schema, AB_SCHEMA);
    assert_eq!(r.first_divergence_tick, None, "{r:?}");
    assert_eq!(r.common_prefix_ticks, r.a.terminal_tick);
    assert!(!r.terminal_divergence);
    assert_eq!(r.a.final_digest, r.b.final_digest);
    // Determinism golden (structure-level: digest equality twice).
    let again = ab_compare("same-input-trace", &spec, None, &spec, None).unwrap();
    assert_eq!(r.receipt_digest, again.receipt_digest);
    jlog("twins", &format!("\"prefix\":{}", r.common_prefix_ticks));
}

#[test]
fn perturbed_axis_diverges_at_a_receipted_tick() {
    let mut a = dec17_scenario(21, PilotMode::FixedControls);
    a.max_ticks = 300;
    let mut b = a.clone();
    b.headwind_mps += 0.5; // the modified axis under the same (empty) trace
    let r = ab_compare("same-input-trace", &a, None, &b, None).unwrap();
    let div = r.first_divergence_tick.expect("must diverge");
    assert!(div >= 1 && div <= 300, "receipted tick: {div}");
    assert_eq!(r.common_prefix_ticks, div - 1);
    assert_ne!(r.a.final_digest, r.b.final_digest);
    assert_ne!(r.a.run_intent_id, r.b.run_intent_id, "axis binds intent");
    jlog("perturbed", &format!("\"first_divergence\":{div}"));
}

#[test]
fn terminal_divergence_records_different_kinds() {
    // A: full Dec-17 fixed run ends EnvelopeExceeded (code 5, measured).
    // B: the same run truncated at 1000 ticks ends MaxTicks (code 4).
    let a = dec17_scenario(1903, PilotMode::FixedControls);
    let mut b = a.clone();
    b.max_ticks = 1000;
    let r = ab_compare("same-input-trace", &a, None, &b, None).unwrap();
    assert!(r.terminal_divergence, "different TerminalEvent kinds");
    assert_eq!(r.a.terminal_code, 5, "envelope exit");
    assert_eq!(r.b.terminal_code, 4, "tick budget");
    // Bitwise-identical up to B's terminal-1 (same physics, shorter
    // budget): the divergence tick is B's terminal (phase code slot).
    let div = r.first_divergence_tick.expect("length divergence");
    assert_eq!(div, 1000, "diverges exactly at B's terminal tick");
    jlog(
        "terminal-divergence",
        &format!(
            "\"a_code\":{},\"b_code\":{},\"div\":{div}",
            r.a.terminal_code, r.b.terminal_code
        ),
    );
}

#[test]
fn human_refly_diverges_where_the_traces_do() {
    let mut spec = dec17_scenario(33, PilotMode::Human);
    spec.max_ticks = 200;
    let ghost = neutral_trace(200);
    let mut refly = neutral_trace(200);
    // The re-fly pulls the lever from tick 80 on.
    for c in refly.iter_mut().skip(79) {
        c.lever_force_n = 60.0;
    }
    let r = ab_compare("human-refly", &spec, Some(&ghost), &spec, Some(&refly)).unwrap();
    let div = r.first_divergence_tick.expect("must diverge");
    assert!(
        (80..=95).contains(&div),
        "diverges when the inputs do (mech responds within ticks): {div}"
    );
    jlog("human-refly", &format!("\"first_divergence\":{div}"));
}

#[test]
fn trace_mode_mismatches_and_caps_refuse() {
    let human = {
        let mut s = dec17_scenario(1, PilotMode::Human);
        s.max_ticks = 10;
        s
    };
    let fixed = {
        let mut s = dec17_scenario(1, PilotMode::FixedControls);
        s.max_ticks = 10;
        s
    };
    // Human without trace refuses; deterministic WITH trace refuses.
    assert_eq!(
        ab_compare("x", &human, None, &human, None)
            .unwrap_err()
            .code,
        "ab-trace-mode-mismatch"
    );
    let t = neutral_trace(10);
    assert_eq!(
        ab_compare("x", &fixed, Some(&t), &fixed, Some(&t))
            .unwrap_err()
            .code,
        "ab-trace-mode-mismatch"
    );
    // Trace caps at cap AND cap+1.
    let at_cap = neutral_trace(MAX_TRACE_LEN);
    assert!(ab_compare("x", &human, Some(&at_cap), &human, Some(&at_cap)).is_ok());
    let over = neutral_trace(MAX_TRACE_LEN + 1);
    assert_eq!(
        ab_compare("x", &human, Some(&over), &human, Some(&over))
            .unwrap_err()
            .code,
        "ab-trace-invalid"
    );
    let empty: Vec<ControlInput> = Vec::new();
    assert_eq!(
        ab_compare("x", &human, Some(&empty), &human, Some(&empty))
            .unwrap_err()
            .code,
        "ab-trace-invalid"
    );
    jlog("refusals", "\"mismatch_and_caps\":true");
}
