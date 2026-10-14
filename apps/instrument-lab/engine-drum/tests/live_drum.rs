//! Native proofs for the live drum. These execute the SAME code the browser
//! runs (the wasm surface is a thin shell over `LiveDrumCore`).
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_exec::CancelGate;
use fs_live_drum::engine::{StrikeOutcome, compile_system};
use fs_live_drum::instrument::LAUNCH_GAP_M;
use fs_live_drum::{DrumParams, Instrument, LiveDrumCore};

#[test]
fn head_frequencies_are_the_reference_eigen_solve() {
    // `percussion drum-modal` prints these on stderr for the same cards.
    let inst = Instrument::build(DrumParams::default()).unwrap();
    let batter = inst.heads[0].frequencies_hz();
    let resonant = inst.heads[1].frequencies_hz();
    let expected_batter = [204.61897307668198, 327.0888597986662, 327.0888597986663,
        436.5665022265004, 436.56650222650046, 447.93617286134685];
    let expected_resonant = [264.1200719782797, 422.07135743711854, 422.0713574371189];
    assert_eq!(batter.len(), 6);
    assert_eq!(resonant.len(), 3);
    for (a, b) in batter.iter().zip(expected_batter) { assert!((a - b).abs() < 1e-9, "{a} vs {b}"); }
    for (a, b) in resonant.iter().zip(expected_resonant) { assert!((a - b).abs() < 1e-9, "{a} vs {b}"); }
    assert!((inst.stick.effective_mass_kg - 0.00724314).abs() < 1e-8);
    // The air spring raises the axisymmetric pair and leaves area-free modes alone.
    assert_eq!(inst.coupled_frequencies_hz.len(), 9);
    assert!(inst.coupled_frequencies_hz[0] > batter[0]);
    assert!(inst.coupled_frequencies_hz.iter().any(|f| (f - batter[1]).abs() < 1e-6));
}

#[test]
fn silence_before_the_first_stroke_is_exact_zero() {
    let mut drum = LiveDrumCore::new("{\"substeps\":4}").unwrap();
    let out = drum.render_f32(512);
    assert_eq!(out.len(), 1024);
    assert!(out.iter().all(|v| *v == 0.0));
    let s = drum.snapshot();
    assert_eq!(s.ledger.stored_j, 0.0);
    assert!(!s.stick_present);
}

#[test]
fn the_voice_runs_exactly_the_compiled_reference_image_until_the_catch() {
    let p = DrumParams { substeps: 4, ..DrumParams::default() };
    let inst = Instrument::build(p).unwrap();
    let gate = CancelGate::new_clock_free();
    let (n0, n1) = (inst.heads[0].mode_count(), inst.heads[1].mode_count());
    let row = inst.batter_shapes_at([0.06, 0.01]).unwrap();
    let root = inst.stick.effective_mass_kg.sqrt();
    let tip = ModalAcousticState { displacement_m_sqrt_kg: -LAUNCH_GAP_M * root, velocity_m_sqrt_kg_per_s: 2.0 * root };
    let rest = vec![ModalAcousticState::default(); n0 + n1];
    let mut direct = compile_system(&inst, 192_000, &rest, Some((tip, &row)), &gate).unwrap();
    let force = vec![0.0; direct.mode_count()];

    let mut drum = LiveDrumCore::new(&p.to_json()).unwrap();
    assert_eq!(drum.strike(0.06, 0.01, 2.0).unwrap(), StrikeOutcome::Queued);
    let mut compared = 0;
    for _ in 0..200 {
        drum.render_pa(1);
        for _ in 0..4 { direct.step(&force, &gate).unwrap(); }
        let s = drum.snapshot();
        if !s.stick_present { break; }
        let state = direct.state();
        for i in 0..n0 + n1 {
            assert_eq!(s.q[i].to_bits(), state[2 * (1 + i)].to_bits());
            assert_eq!(s.v[i].to_bits(), state[2 * (1 + i) + 1].to_bits());
        }
        compared += 1;
    }
    assert!(compared > 20, "the stick must stay in the system through contact ({compared} frames)");
}

#[test]
fn ledger_closes_across_launch_contact_catch_and_restrikes() {
    let mut drum = LiveDrumCore::new("{\"substeps\":4}").unwrap();
    let mut launched = 0;
    let mut saw_contact = false;
    let mut saw_catch = false;
    let mut heard = 0.0_f64;
    for frame in 0..24_000 {
        if frame % 3_000 == 0 {
            drum.strike(0.05 - 0.01 * f64::from(launched % 3), 0.02, 1.0 + f64::from(launched)).unwrap();
            launched += 1;
        }
        let (batter, _) = drum.render_pa(1);
        heard = heard.max(batter[0].abs());
        assert!(drum.take_refusal().is_none());
        if frame % 16 == 0 {
            let s = drum.snapshot();
            saw_contact |= s.stick_in_contact;
            saw_catch |= s.strikes_launched > 0 && !s.stick_present;
            let scale = s.ledger.supplied_j.max(1e-6);
            assert!(s.ledger.balance_residual_j.abs() < 1e-7 * scale + 1e-10,
                "frame {frame}: residual {:e}", s.ledger.balance_residual_j);
            assert!(s.ledger.dissipated_j >= -1e-12);
        }
    }
    let s = drum.snapshot();
    assert_eq!(s.strikes_launched, 8);
    assert!(saw_contact && saw_catch && heard > 1.0);
    // Supplied = sum of 1/2 m v^2 over the eight launches.
    let m = drum.instrument().stick.effective_mass_kg;
    let expected: f64 = (0..8).map(|k| 0.5 * m * (1.0 + f64::from(k)).powi(2)).sum();
    assert!((s.ledger.supplied_j - expected).abs() < 1e-12 * expected.max(1.0));
    assert!(s.ledger.recovered_j > 0.0 && s.ledger.recovered_j < s.ledger.supplied_j);
    assert!(s.ledger.dissipated_j > 0.0);
}

#[test]
fn a_stroke_during_contact_waits_for_separation_and_deletes_no_contact_energy() {
    let mut drum = LiveDrumCore::new("{\"substeps\":4}").unwrap();
    drum.strike(0.06, 0.01, 2.0).unwrap();
    let mut deferred = false;
    for _ in 0..2_000 {
        drum.render_pa(1);
        let s = drum.snapshot();
        if s.stick_in_contact {
            assert_eq!(drum.strike(-0.03, 0.04, 3.0).unwrap(), StrikeOutcome::DeferredUntilSeparation);
            deferred = true;
            let pending = drum.snapshot();
            assert!(pending.strike_pending && pending.strikes_launched == 1);
            break;
        }
    }
    assert!(deferred, "contact was never observed");
    let mut relaunched = false;
    for _ in 0..2_000 {
        drum.render_pa(1);
        let s = drum.snapshot();
        assert!(s.ledger.balance_residual_j.abs() < 1e-9, "{:e}", s.ledger.balance_residual_j);
        if s.strikes_launched == 2 {
            assert_eq!(s.stick_xy_m, [-0.03, 0.04]);
            relaunched = true;
            break;
        }
        // Until separation the ORIGINAL stick is still the one in contact.
        assert_eq!(s.stick_xy_m, [0.06, 0.01]);
    }
    assert!(relaunched);
}

#[test]
fn typed_refusals_leave_the_running_instrument_untouched() {
    let mut drum = LiveDrumCore::new("{\"substeps\":2}").unwrap();
    drum.strike(0.06, 0.01, 2.0).unwrap();
    drum.render_pa(480);
    let before = drum.snapshot();
    assert_eq!(drum.strike(0.5, 0.0, 2.0).unwrap_err().code, "strike_outside_head");
    assert_eq!(drum.strike(0.0, 0.0, 0.0).unwrap_err().code, "strike_speed");
    assert_eq!(drum.strike(0.0, 0.0, 25.0).unwrap_err().code, "strike_speed");
    assert_eq!(drum.strike(f64::NAN, 0.0, 1.0).unwrap_err().code, "strike_outside_head");
    assert_eq!(drum.rebuild("{\"batter_tension_n_m\":-1}").unwrap_err().code, "bad_parameter");
    assert_eq!(drum.rebuild("{\"pitch\":440}").unwrap_err().code, "unknown_parameter");
    assert_eq!(drum.rebuild("{\"band_low_hz\":5,\"band_high_hz\":20}").unwrap_err().code, "empty_band");
    // Far over budget: refused by the Weyl estimate before any eigen-solve.
    assert_eq!(drum.rebuild("{\"band_high_hz\":4000,\"radial_intervals\":16,\"azimuths\":64}").unwrap_err().code, "mode_budget");
    // Slightly over budget: the estimate lets it through; the eigen-solve counts and refuses.
    let over = drum.rebuild("{\"band_low_hz\":60,\"band_high_hz\":1250,\"radial_intervals\":14,\"azimuths\":64}").unwrap_err();
    assert_eq!(over.code, "mode_budget", "{over}");
    let after = drum.snapshot();
    assert_eq!(before.q, after.q);
    assert_eq!(before.ledger.stored_j, after.ledger.stored_j);
    assert_eq!(drum.instrument().params.batter.tension_n_m, 3000.0);
}

#[test]
fn rebuild_starts_silent_and_changes_frequencies_only_through_the_eigen_solve() {
    let mut drum = LiveDrumCore::new("{\"substeps\":2}").unwrap();
    drum.strike(0.06, 0.01, 2.0).unwrap();
    drum.render_pa(960);
    assert!(drum.snapshot().ledger.stored_j > 0.0);
    let f0 = drum.instrument().heads[0].frequencies_hz()[0];
    // Membrane limit: frequency scales with sqrt(tension); bending adds a little.
    assert!(drum.rebuild("{\"batter_tension_n_m\":6000}").unwrap());
    let f1 = drum.instrument().heads[0].frequencies_hz()[0];
    assert!((f1 / f0 - 2.0_f64.sqrt()).abs() < 5e-3, "{f0} -> {f1}");
    let s = drum.snapshot();
    assert_eq!(s.ledger.stored_j, 0.0);
    assert_eq!(s.strikes_launched, 0);
    assert!(drum.render_f32(256).iter().all(|v| *v == 0.0));
    // Observation/clock-only changes reuse the eigen-solve.
    assert!(!drum.rebuild("{\"substeps\":4,\"full_scale_pa\":50,\"depth_m\":0.2}").unwrap());
    assert_eq!(drum.voice().substeps(), 4);
    assert_eq!(drum.instrument().heads[0].frequencies_hz()[0], f1);
    // Radius: frequency ~ 1/radius.
    assert!(drum.rebuild("{\"radius_m\":0.12}").unwrap());
    let f2 = drum.instrument().heads[0].frequencies_hz();
    assert!(f2[0] > f1 * 1.3);
}

#[test]
fn tension_variation_splits_the_degenerate_pair_through_the_actual_pencil() {
    let uniform = Instrument::build(DrumParams::default()).unwrap();
    let f = uniform.heads[0].frequencies_hz();
    assert!((f[1] - f[2]).abs() < 1e-6);
    let varied = Instrument::build(DrumParams::default()
        .overlay_json("{\"batter_tension_variation\":[300,-200,80,600,-200,100,400]}").unwrap()).unwrap();
    let g = varied.heads[0].frequencies_hz();
    assert!((g[1] - g[2]).abs() > 1.0, "{g:?}");
}

#[test]
fn info_json_is_parseable_and_declares_scope() {
    let drum = LiveDrumCore::new("{}").unwrap();
    let info = drum.info_json();
    assert!(info.starts_with('{') && info.ends_with('}'));
    for needle in ["\"mechanical_rate_hz\":768000", "\"frequencies_hz\":[204.6", "\"no_claim\":", "\"provenance\":",
        "\"snapshot_layout\":", "\"coupled_frequencies_hz\":[218.0"] {
        assert!(info.contains(needle), "{needle} missing from {info}");
    }
    // Balanced braces/brackets outside strings: a cheap structural check.
    let (mut depth, mut in_string, mut escape) = (0i32, false, false);
    for ch in info.chars() {
        if in_string {
            if escape { escape = false; } else if ch == '\\' { escape = true; } else if ch == '"' { in_string = false; }
        } else {
            match ch { '"' => in_string = true, '{' | '[' => depth += 1, '}' | ']' => depth -= 1, _ => {} }
            assert!(depth >= 0);
        }
    }
    assert_eq!(depth, 0);
}
