use super::*;
use crate::modal_acoustic_time::broadband::{BroadbandModalControls, ModalTransferSamples};
use fs_math::{c64::C64, det};
use fs_modal::time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};

const RATE: u32 = 8000;
const NATURAL: f64 = core::f64::consts::TAU * 700.0;
const DRIVE: f64 = core::f64::consts::TAU * 2600.0;
const LIMITS: &str = "tone_limits 8 1048576 10000000 100000\n";

fn observer(damping: f64) -> BroadbandModalObserver {
    let model = ModalAcousticTimeModel::try_new(
        RATE,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: NATURAL,
            damping_ratio: damping,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let samples = ModalTransferSamples {
        source_id: "analytic instantaneous velocity receiver".into(),
        propagation_delay_s: 0.0,
        producer_relative_error_estimate: 0.0,
        omega_rad_s: vec![100.0, 300.0],
        per_mode: vec![vec![C64::from_re(2.0); 2]],
        held_out: vec![(200.0, vec![C64::from_re(2.0)])],
    };
    BroadbandModalObserver::build(
        &model,
        &samples,
        BroadbandModalControls {
            fit_order: 1,
            fit_iterations: 4,
            delay_kernel_half_len: 1,
            max_held_out_rel_error: 1e-9,
        },
        RATE,
    )
    .unwrap()
    .0
}

fn program(text: &str, ports: &[Vec<f64>]) -> Result<ToneProgram, String> {
    ToneProgram::read(text, 0, &[ports], RATE, &mut 10000)
}

fn driver(text: &str, ports: &[Vec<f64>], observer: &BroadbandModalObserver) -> ToneDriver {
    let program = program(text, ports).unwrap();
    program
        .sources
        .into_iter()
        .next()
        .unwrap()
        .unwrap()
        .prepare(observer, program.budget)
        .unwrap()
}

fn velocity(t: f64, omega: f64, cosine: f64, sine: f64, dc: f64) -> f64 {
    let denominator = NATURAL * NATURAL - omega * omega;
    dc * det::sin(NATURAL * t) / NATURAL
        + cosine * (-omega * det::sin(omega * t) + NATURAL * det::sin(NATURAL * t)) / denominator
        + sine * omega * (det::cos(omega * t) - det::cos(NATURAL * t)) / denominator
}

#[test]
fn prescribed_port_tone_matches_exact_motion_and_exposes_sample_hold_bias() {
    let text = format!("{LIMITS}tones 1\ntone 0 0 {DRIVE:.17e} 0.03 -0.01\n");
    let mut exact = observer(0.0);
    let mut sampled = observer(0.0);
    let mut source = driver(&text, &[vec![-2.0]], &exact);
    let mut held_error = 0.0;
    let mut reference_energy = 0.0;
    let mut work = 0.0;
    for i in 0..257 {
        let t = (i + 1) as f64 / f64::from(RATE);
        let expected = 2.0 * velocity(t, DRIVE, -0.06, 0.02, 0.015);
        let observed = source.step(&mut exact, &[0.015]).unwrap();
        assert!(
            (observed - expected).abs() < 2e-15,
            "sample={i}: {observed} != {expected}"
        );
        let before = i as f64 / f64::from(RATE);
        let held = 0.015 - 0.06 * det::cos(DRIVE * before) + 0.02 * det::sin(DRIVE * before);
        let wrong = sampled.step(&[held]).unwrap();
        held_error += (wrong - expected).powi(2);
        reference_energy += expected * expected;
        work += exact.last_modal_frame().input_work_j;
        assert_eq!(exact.last_modal_frame().viscous_dissipation_j, 0.0);
    }
    assert!((held_error / reference_energy).sqrt() > 0.1);
    assert!((work - exact.last_modal_frame().total_modal_energy_j).abs() < 1e-17);
}

#[test]
fn repeated_port_frequencies_cancel_coherently_and_close_tones_keep_phase() {
    let other = DRIVE * 1.0017;
    let text = format!(
        "{LIMITS}tones 3\ntone 0 0 {DRIVE:.17e} 0.03 -0.01\ntone 0 1 {DRIVE:.17e} 0.03 -0.01\ntone 0 0 {other:.17e} 0.002 0.001\n"
    );
    let mut observed = observer(0.0);
    let mut source = driver(&text, &[vec![1.0], vec![-1.0]], &observed);
    for sample in 1..=359 {
        let p = source.step(&mut observed, &[0.0]).unwrap();
        let expected = 2.0 * velocity(sample as f64 / f64::from(RATE), other, 0.002, 0.001, 0.0);
        assert!((p - expected).abs() < 2e-15);
    }
}

#[test]
fn refused_step_preserves_phase_and_accepted_work_heat() {
    let text = format!("{LIMITS}tones 1\ntone 0 0 {DRIVE:.17e} 0.03 -0.01\n");
    let ports = [vec![1.0]];
    let (mut a, mut b) = (observer(0.03), observer(0.03));
    let (mut sa, mut sb) = (driver(&text, &ports, &a), driver(&text, &ports, &b));
    for _ in 0..23 {
        sa.step(&mut a, &[0.01]).unwrap();
        sb.step(&mut b, &[0.01]).unwrap();
    }
    let c = sa.cosine.clone();
    let s = sa.sine.clone();
    let frame = a.last_modal_frame().clone();
    let states = a.model().states().to_vec();
    assert!(sa.step(&mut a, &[f64::NAN]).is_err());
    assert_eq!(&sa.cosine[1..], &c[1..]);
    assert_eq!(&sa.sine[1..], &s[1..]);
    assert_eq!(a.model().states(), states);
    assert_eq!(a.last_modal_frame(), &frame);
    for _ in 0..71 {
        assert_eq!(
            sa.step(&mut a, &[0.01]).unwrap().to_bits(),
            sb.step(&mut b, &[0.01]).unwrap().to_bits()
        );
        assert_eq!(a.last_modal_frame(), b.last_modal_frame());
        assert!(a.last_modal_frame().viscous_dissipation_j >= 0.0);
    }
}

#[test]
fn continuous_record_clock_port_projection_and_preparation_limits_refuse() {
    let good = format!("{LIMITS}tones 1\ntone 0 0 {DRIVE:.17e} 0.03 -0.01\n");
    for bad in [
        good.replace("tones 1", "tones 999999999"),
        good.replace("tone 0 0", "tone 1 0"),
        good.replace("tone 0 0", "tone 0 1"),
        good.replace(" 0.03", " NaN"),
        good.replace(" 0.03", " inf"),
        good.replace("tone_limits 8", "tone_limits 1"),
        good.replace("1048576", "1"),
        format!("{good}ignored\n"),
        format!("{LIMITS}tones 1\ntone 0 0 26000 0.03 -0.01\n"),
    ] {
        assert!(program(&bad, &[vec![1.0]]).is_err(), "accepted {bad}");
    }
    assert!(ToneProgram::read(&good, 0, &[&[vec![1.0]]], RATE, &mut 1).is_err());
    for bad in [
        good.replace("10000000", "1"),
        good.replace("100000\n", "1\n"),
    ] {
        let program = program(&bad, &[vec![1.0]]).unwrap();
        assert!(
            program
                .sources
                .into_iter()
                .next()
                .unwrap()
                .unwrap()
                .prepare(&observer(0.01), program.budget)
                .is_err()
        );
    }
    // Below Nyquist but outside the mechanical owner's original 0.9 guard.
    let guard = format!("{LIMITS}tones 1\ntone 0 0 24000 0.03 -0.01\n");
    let program = program(&guard, &[vec![1.0]]).unwrap();
    assert!(
        program
            .sources
            .into_iter()
            .next()
            .unwrap()
            .unwrap()
            .prepare(&observer(0.01), program.budget)
            .is_err()
    );
}
