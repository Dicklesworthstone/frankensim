//! The live drum's mechanics on the composition kernel (`"mechanics":"kernel"`,
//! cargo feature `kernel`): fs-couple's `ImpactGraph` declares the voice's
//! description on fs-system, against the module (`LinearImpactSystem`) that
//! compiles the same description.
//!
//! Asserted: over a scripted run of three strikes, both microphones' pressure
//! within 1e-6 of the module's peak at every frame, and each modal
//! coordinate within 1e-6 of its full scale; the ledger closing across
//! launches, contact, catches and restrikes on the kernel alone; a stroke
//! during contact waiting for separation; a mechanics change reusing the
//! heads and starting silent, a tension change re-solving them; the kernel's
//! refusals by name, the running drum untouched.
//!
//! Not covered: allocation in the voice's step (this crate forbids unsafe
//! code, so it holds no counting allocator; fs-couple's
//! `tests/impact_graph.rs` gates `ImpactGraph::step`), the browser build
//! (the default and wasm builds have no kernel), timing (`examples/bench.rs`).

use fs_live_drum::engine::StrikeOutcome;
use fs_live_drum::{LiveDrumCore, Mechanics};

const KERNEL: &str = r#"{"substeps":4,"mechanics":"kernel"}"#;

/// Frame, station [m] and speed [m/s] of each scripted strike.
const STRIKES: [(usize, [f64; 2], f64); 3] = [
    (0, [0.06, 0.01], 2.0),
    (3_000, [-0.03, 0.04], 3.0),
    (6_000, [0.02, -0.05], 1.5),
];
const FRAMES: usize = 9_000;

/// Both microphones per frame, and the interleaved modal state per frame.
struct Run {
    pressure: [Vec<f64>; 2],
    state: Vec<Vec<f64>>,
    drum: LiveDrumCore,
}

fn run(params: &str) -> Run {
    let mut drum = LiveDrumCore::new(params).unwrap();
    let mut pressure = [Vec::with_capacity(FRAMES), Vec::with_capacity(FRAMES)];
    let mut state = Vec::with_capacity(FRAMES);
    for frame in 0..FRAMES {
        for (at, xy, speed) in STRIKES {
            if at == frame {
                drum.strike(xy[0], xy[1], speed).unwrap();
            }
        }
        let (batter, resonant) = drum.render_pa(1);
        pressure[0].push(batter[0]);
        pressure[1].push(resonant[0]);
        assert!(drum.take_refusal().is_none(), "frame {frame}");
        let s = drum.snapshot();
        state.push(s.q.iter().zip(&s.v).flat_map(|(&q, &v)| [q, v]).collect());
    }
    Run {
        pressure,
        state,
        drum,
    }
}

#[test]
fn the_kernel_mechanics_render_the_modules_pressure_at_both_microphones() {
    let mut module = run(r#"{"substeps":4}"#);
    let mut kernel = run(KERNEL);
    assert_eq!(module.drum.instrument().params.mechanics, Mechanics::Module);
    assert_eq!(kernel.drum.instrument().params.mechanics, Mechanics::Kernel);

    let full_scale_pa = module.drum.instrument().params.full_scale_pa;
    let mut worst_db = [0.0; 2];
    for (mic, ((m, k), db)) in module
        .pressure
        .iter()
        .zip(&kernel.pressure)
        .zip(&mut worst_db)
        .enumerate()
    {
        let peak = m.iter().fold(0.0_f64, |a, p| a.max(p.abs()));
        let (at, worst) = m
            .iter()
            .zip(k)
            .map(|(a, b)| (a - b).abs())
            .enumerate()
            .fold(
                (0, 0.0_f64),
                |(i, w), (j, d)| if d > w { (j, d) } else { (i, w) },
            );
        assert!(peak > 1.0, "microphone {mic}: peak {peak} Pa");
        assert!(
            worst <= 1e-6 * peak,
            "microphone {mic}: |dp| {worst:e} Pa at frame {at}, peak {peak} Pa"
        );
        *db = 20.0 * (worst / full_scale_pa).log10();
    }

    let width = module.state[0].len();
    assert_eq!(kernel.state[0].len(), width);
    let mut worst_state = 0.0_f64;
    for c in 0..width {
        let scale = module.state.iter().fold(0.0_f64, |a, s| a.max(s[c].abs()));
        let gap = module
            .state
            .iter()
            .zip(&kernel.state)
            .fold(0.0_f64, |a, (m, k)| a.max((m[c] - k[c]).abs()));
        assert!(
            gap <= 1e-6 * scale,
            "coordinate {c}: gap {gap:e}, full scale {scale:e}"
        );
        if scale > 0.0 {
            worst_state = worst_state.max(gap / scale);
        }
    }

    let (m, k) = (module.drum.snapshot(), kernel.drum.snapshot());
    assert_eq!((m.strikes_launched, k.strikes_launched), (3, 3));
    assert_eq!(m.mechanical_steps, k.mechanical_steps);
    assert_eq!(m.ledger.supplied_j.to_bits(), k.ledger.supplied_j.to_bits());
    let scale = m.ledger.supplied_j;
    assert!((m.ledger.dissipated_j - k.ledger.dissipated_j).abs() <= 1e-9 * scale);
    assert!((m.ledger.recovered_j - k.ledger.recovered_j).abs() <= 1e-9 * scale);
    assert!(
        k.ledger.balance_residual_j.abs() < 1e-9 * scale,
        "{:e}",
        k.ledger.balance_residual_j
    );
    println!(
        "{{\"suite\":\"fs-live-drum\",\"case\":\"kernel-vs-module\",\"frames\":{FRAMES},\"strikes\":3,\
         \"worst_db_re_full_scale\":[{:.2},{:.2}],\"worst_state_over_full_scale\":{worst_state:.3e},\
         \"residual_j\":[{:e},{:e}]}}",
        worst_db[0], worst_db[1], m.ledger.balance_residual_j, k.ledger.balance_residual_j
    );
}

#[test]
fn the_kernel_ledger_closes_across_launch_contact_catch_and_restrikes() {
    let mut drum = LiveDrumCore::new(KERNEL).unwrap();
    let mut launched = 0;
    let mut saw_contact = false;
    let mut saw_catch = false;
    let mut heard = 0.0_f64;
    for frame in 0..12_000 {
        if frame % 3_000 == 0 {
            drum.strike(
                0.05 - 0.01 * f64::from(launched % 3),
                0.02,
                1.0 + f64::from(launched),
            )
            .unwrap();
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
            assert!(
                s.ledger.balance_residual_j.abs() < 1e-7 * scale + 1e-10,
                "frame {frame}: residual {:e}",
                s.ledger.balance_residual_j
            );
            assert!(s.ledger.dissipated_j >= -1e-12);
        }
    }
    let s = drum.snapshot();
    assert_eq!(s.strikes_launched, 4);
    assert!(saw_contact && saw_catch && heard > 1.0);
    let m = drum.instrument().stick.effective_mass_kg;
    let expected: f64 = (0..4).map(|k| 0.5 * m * (1.0 + f64::from(k)).powi(2)).sum();
    assert!((s.ledger.supplied_j - expected).abs() < 1e-12 * expected.max(1.0));
    assert!(s.ledger.recovered_j > 0.0 && s.ledger.recovered_j < s.ledger.supplied_j);
    assert!(s.ledger.dissipated_j > 0.0);
}

#[test]
fn a_kernel_stroke_during_contact_waits_for_separation() {
    let mut drum = LiveDrumCore::new(KERNEL).unwrap();
    drum.strike(0.06, 0.01, 2.0).unwrap();
    let mut deferred = false;
    for _ in 0..2_000 {
        drum.render_pa(1);
        if drum.snapshot().stick_in_contact {
            assert_eq!(
                drum.strike(-0.03, 0.04, 3.0).unwrap(),
                StrikeOutcome::DeferredUntilSeparation
            );
            let pending = drum.snapshot();
            assert!(pending.strike_pending && pending.strikes_launched == 1);
            deferred = true;
            break;
        }
    }
    assert!(deferred, "contact was never observed");
    let mut relaunched = false;
    for _ in 0..2_000 {
        drum.render_pa(1);
        let s = drum.snapshot();
        assert!(
            s.ledger.balance_residual_j.abs() < 1e-9,
            "{:e}",
            s.ledger.balance_residual_j
        );
        if s.strikes_launched == 2 {
            assert_eq!(s.stick_xy_m, [-0.03, 0.04]);
            relaunched = true;
            break;
        }
        assert_eq!(s.stick_xy_m, [0.06, 0.01]);
    }
    assert!(relaunched);
}

#[test]
fn a_mechanics_change_reuses_the_heads_and_a_tension_change_resolves_them() {
    let mut drum = LiveDrumCore::new(r#"{"substeps":2}"#).unwrap();
    drum.strike(0.06, 0.01, 2.0).unwrap();
    drum.render_pa(960);
    assert!(drum.snapshot().ledger.stored_j > 0.0);
    let f0 = drum.instrument().heads[0].frequencies_hz();

    assert!(
        !drum.rebuild(r#"{"mechanics":"kernel"}"#).unwrap(),
        "a mechanics change reuses the heads"
    );
    assert_eq!(drum.instrument().params.mechanics, Mechanics::Kernel);
    assert_eq!(drum.instrument().heads[0].frequencies_hz(), f0);
    let s = drum.snapshot();
    assert_eq!((s.ledger.stored_j, s.strikes_launched), (0.0, 0));
    assert!(drum.render_f32(256).iter().all(|v| *v == 0.0));

    assert!(drum.rebuild(r#"{"batter_tension_n_m":6000}"#).unwrap());
    assert_eq!(drum.instrument().params.mechanics, Mechanics::Kernel);
    let f1 = drum.instrument().heads[0].frequencies_hz()[0];
    assert!(
        (f1 / f0[0] - 2.0_f64.sqrt()).abs() < 5e-3,
        "{} -> {f1}",
        f0[0]
    );
    assert_eq!(drum.snapshot().ledger.stored_j, 0.0);
    assert!(drum.render_f32(256).iter().all(|v| *v == 0.0));
    drum.strike(0.06, 0.01, 2.0).unwrap();
    assert!(drum.render_f32(1024).iter().any(|p| p.abs() > 1e-6));
    assert!(drum.take_refusal().is_none());
}

#[test]
fn the_kernel_refuses_contact_loss_by_name_and_leaves_the_running_drum_untouched() {
    let lossy = LiveDrumCore::new(r#"{"mechanics":"kernel","contact_loss_s_m":0.1}"#)
        .err()
        .unwrap();
    assert_eq!(lossy.code, "kernel_refused", "{lossy}");
    let mut drum = LiveDrumCore::new(r#"{"substeps":2,"mechanics":"kernel"}"#).unwrap();
    drum.strike(0.06, 0.01, 2.0).unwrap();
    drum.render_pa(480);
    let before = drum.snapshot();
    assert_eq!(
        drum.rebuild(r#"{"contact_loss_s_m":0.1}"#)
            .unwrap_err()
            .code,
        "kernel_refused"
    );
    assert_eq!(
        drum.rebuild(r#"{"mechanics":"fast"}"#).unwrap_err().code,
        "bad_parameter"
    );
    let after = drum.snapshot();
    assert_eq!(before.q, after.q);
    assert_eq!(before.ledger.stored_j, after.ledger.stored_j);
    assert_eq!(drum.instrument().params.mechanics, Mechanics::Kernel);
    // The module keeps admitting the loss.
    assert!(LiveDrumCore::new(r#"{"contact_loss_s_m":0.1}"#).is_ok());
}
