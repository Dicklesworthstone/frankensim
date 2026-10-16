//! G0/G1: exact held-force propagation, sampled-delay refinement and the
//! independent pulsating-sphere operator. No recording or fitted oracle.
use fs_couple::broadband_radiation::BroadbandRadiationAuthority;
use fs_couple::broadband_radiation::modal::{
    MAX_MODAL_DELAY_KERNEL_HALF_LEN, MAX_MODAL_TRANSFER_ORDER, admit_modal_samples,
};
use fs_couple::modal_acoustic_time::broadband::{
    BroadbandModalControls, BroadbandModalError, BroadbandModalObserver, ModalTransferSamples,
};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_math::{c64::C64, det};
use std::time::Instant;

const RATE: u32 = 48_000;
const W0: f64 = core::f64::consts::TAU * 500.0;
const ZETA: f64 = 0.01;
const RHO: f64 = 1.2;
const SOUND_SPEED: f64 = 343.0;
const RADIUS: f64 = 0.05;
const DISTANCE: f64 = 1.0;
const DELAY: f64 = (DISTANCE - RADIUS) / SOUND_SPEED;

fn phase(x: f64) -> C64 {
    C64::new(det::cos(x), det::sin(x))
}
fn model(count: usize) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(
        RATE,
        vec![
            ModalAcousticMode {
                angular_frequency_rad_s: W0,
                damping_ratio: ZETA,
                pressure_per_modal_velocity: C64::ZERO,
            };
            count
        ],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap()
}
fn controls(order: usize, half: usize, tolerance: f64) -> BroadbandModalControls {
    BroadbandModalControls {
        fit_order: order,
        fit_iterations: 32,
        delay_kernel_half_len: half,
        max_held_out_rel_error: tolerance,
    }
}
fn sphere(w: f64) -> C64 {
    C64::new(0.0, -w * RHO)
        * phase(w * DELAY)
        * C64::new(1.0, -w * RADIUS / SOUND_SPEED)
            .recip()
            .scale(1.0 / (4.0 * core::f64::consts::PI * DISTANCE))
}
fn samples<F: Fn(f64) -> C64>(n: usize, delay: f64, h: F) -> ModalTransferSamples {
    let low = W0 * 0.125;
    let high = W0 * 9.0;
    let omega: Vec<_> = (0..n)
        .map(|i| low + (high - low) * i as f64 / (n - 1) as f64)
        .collect();
    let held_out = (0..127)
        .map(|i| {
            let w = low + (high - low) * (i as f64 + 0.371) / 127.0;
            (w, vec![h(w)])
        })
        .collect();
    ModalTransferSamples {
        source_id: "analytic:sphere:a=.05:r=1:rho=1.2:c=343".into(),
        per_mode: vec![omega.iter().map(|&w| h(w)).collect()],
        omega_rad_s: omega,
        held_out,
        propagation_delay_s: delay,
        producer_relative_error_estimate: 0.0,
    }
}

// Independent partial fractions of G(s)=C*s^2/[(1+tau*s)*
// (s^2+2*zeta*w0*s+w0^2)]. Integrate EACH exponential exactly over a
// held input. This includes ZOH images; a continuous sinusoid is a different
// forcing problem. The final phase is the declared bandlimited sampled delay,
// not a claim to retime the continuous held-force aliases within a sample.
fn held_sphere_reference(w: f64) -> C64 {
    let tau = RADIUS / SOUND_SPEED;
    let wd = W0 * det::sqrt(1.0 - ZETA * ZETA);
    let poles = [
        C64::from_re(-1.0 / tau),
        C64::new(-ZETA * W0, wd),
        C64::new(-ZETA * W0, -wd),
    ];
    let dt = 1.0 / f64::from(RATE);
    let mut result = C64::ZERO;
    for (i, &p) in poles.iter().enumerate() {
        let denominator = poles
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .fold(C64::from_re(tau), |d, (_, &other)| d * (p - other));
        let residue = p
            * p
            * denominator
                .recip()
                .scale(RHO / (4.0 * core::f64::consts::PI * DISTANCE));
        let a = phase(p.im * dt).scale(det::exp(p.re * dt));
        let b = (a - C64::ONE) * p.recip();
        result = result + residue * b * (phase(-w * dt) - a).recip();
    }
    result * phase(w * DELAY)
}

#[test]
fn off_resonance_transfer_matches_pulsating_sphere() {
    let start = Instant::now();
    let s = samples(256, DELAY, sphere);
    let mut maximum = 0.0_f64;
    for ratio in [0.25, 0.5, 2.0, 4.0, 8.0] {
        let w = ratio * W0;
        let (mut broadband, report) =
            BroadbandModalObserver::build(&model(1), &s, controls(1, 128, 1e-6), RATE).unwrap();
        let mut legacy = model(1);
        let mut measured = C64::ZERO;
        let mut narrow = C64::ZERO;
        // One second removes the slowest modal transient by exp(-31.4).
        // The following second contains an integer number of every drive.
        for n in 0..2 * RATE {
            let force = det::cos(w * f64::from(n) / f64::from(RATE));
            let pressure = broadband.step(&[force]).unwrap();
            legacy.step(&[force]).unwrap();
            let old = legacy
                .observer_pressure_with_transfers_about_static_equilibrium(&[sphere(W0)], &[force])
                .unwrap();
            if n >= RATE {
                let carrier = phase(w * f64::from(n + 1) / f64::from(RATE));
                measured = measured + carrier.scale(pressure * 2.0 / f64::from(RATE));
                narrow = narrow + carrier.scale(old * 2.0 / f64::from(RATE));
            }
        }
        let reference = held_sphere_reference(w);
        let quotient = measured * reference.recip();
        let amplitude = (quotient.abs() - 1.0).abs();
        let phase_error = det::atan2(quotient.im, quotient.re).abs();
        maximum = maximum.max((measured - reference).abs() / reference.abs());
        println!(
            "{{\"suite\":\"broadband-modal\",\"case\":\"off-resonance\",\"ratio\":{ratio},\"producer_id\":\"{}\",\"oracle\":\"analytic exact held-force endpoint plus bandlimited sampled delay\",\"amplitude_rel_error\":{amplitude},\"phase_error_rad\":{phase_error},\"legacy_rel_error\":{},\"held_out_max_rel\":{},\"delay_samples\":{},\"tolerance\":1e-6,\"wall_s\":{}}}",
            report.source_id,
            (narrow - reference).abs() / reference.abs(),
            report.held_out_max_rel_error,
            report.delay_samples,
            start.elapsed().as_secs_f64()
        );
        assert!(amplitude <= 1e-6, "ratio {ratio}: {amplitude}");
        assert!(phase_error <= 1e-6, "ratio {ratio}: {phase_error}");
    }
    assert!(maximum <= 2e-6);
}

#[test]
fn constant_transfer_reduces_to_narrow_band_bitwise() {
    let h = C64::new(2.1, -0.37);
    let s = samples(3, 0.0, |_| h);
    let mut original = model(1);
    original
        .restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: 1e-4,
            velocity_m_sqrt_kg_per_s: -0.03,
        }])
        .unwrap();
    let (mut broadband, report) =
        BroadbandModalObserver::build(&original, &s, controls(1, 1, 1e-6), RATE).unwrap();
    assert!(report.constant_narrow_band_compatibility);
    for n in 0..1024 {
        let f = if n < 317 { 2.3 } else { -0.81 };
        original.step(&[f]).unwrap();
        let expected = original
            .observer_pressure_with_transfers_about_static_equilibrium(&[h], &[f])
            .unwrap();
        assert_eq!(broadband.step(&[f]).unwrap().to_bits(), expected.to_bits());
    }
}

#[test]
fn energy_free_observer_does_not_load_the_mode() {
    let mut original = model(1);
    let s = samples(128, DELAY, sphere);
    let (mut observer, report) =
        BroadbandModalObserver::build(&original, &s, controls(1, 128, 1e-6), RATE).unwrap();
    assert_eq!(report.authority, BroadbandRadiationAuthority::EstimateOnly);
    for n in 0..4096 {
        let f = det::sin(n as f64 * 0.139);
        let frame = original.step(&[f]).unwrap();
        observer.step(&[f]).unwrap();
        assert_eq!(original.states(), observer.model().states());
        let other = observer.last_modal_frame();
        assert_eq!(
            frame.total_modal_energy_j.to_bits(),
            other.total_modal_energy_j.to_bits()
        );
        assert_eq!(frame.input_work_j.to_bits(), other.input_work_j.to_bits());
        assert_eq!(
            frame.viscous_dissipation_j.to_bits(),
            other.viscous_dissipation_j.to_bits()
        );
    }
    let before = observer.model().states().to_vec();
    assert!(observer.step(&[f64::NAN]).is_err());
    assert_eq!(before, observer.model().states());
}

fn maximum_operator_error(observer: &BroadbandModalObserver, s: &ModalTransferSamples) -> f64 {
    s.held_out
        .iter()
        .map(|(w, row)| (observer.transfer_at(0, *w).unwrap() - row[0]).abs() / row[0].abs())
        .fold(0.0_f64, f64::max)
}
fn piston(w: f64) -> C64 {
    let (r, x) = fs_phs::baffled_piston_impedance(RHO, SOUND_SPEED, RADIUS, w).unwrap();
    C64::new(r, x).scale(core::f64::consts::PI * RADIUS * RADIUS)
}
#[test]
fn order_ladder_converges_and_axis_is_live() {
    // The delay-extracted sphere is already first-order rational. It cannot
    // honestly exercise an order ladder. Use the independent nonrational
    // Bessel/Struve piston owner for this axis, as permitted by the bead.
    let s = samples(256, 0.0, piston);
    let mut errors = Vec::new();
    let mut coarse = C64::ZERO;
    let mut fine = C64::ZERO;
    for order in [2, 6, 12] {
        let (observer, report) =
            BroadbandModalObserver::build(&model(1), &s, controls(order, 1, 10.0), RATE).unwrap();
        let error = maximum_operator_error(&observer, &s);
        if order == 2 {
            coarse = observer.transfer_at(0, 4.0 * W0).unwrap();
        }
        fine = observer.transfer_at(0, 4.0 * W0).unwrap();
        println!(
            "{{\"suite\":\"broadband-modal\",\"axis\":\"poles\",\"pole_count\":{order},\"held_out_max_rel\":{error},\"producer_id\":\"analytic:baffled-piston:Bessel-Struve\",\"authority\":\"{:?}\",\"tolerance\":1e-6}}",
            report.authority
        );
        errors.push(error);
    }
    assert!(errors.windows(2).all(|e| e[1] < e[0]), "{errors:?}");
    assert!(errors[2] <= 1e-6, "{errors:?}");
    assert!(
        (coarse - fine).abs() > 1e-8 * fine.abs(),
        "pole axis is not live"
    );
}

#[test]
fn fit_grid_density_converges() {
    // Hold the model order fixed while refining only the sample grid. The
    // independent pole ladder resolves this piston to 7.82e-7 at order six.
    // Order twelve over-resolves it: batch 8 (541935a25, ts1) measured
    // [1.483e-12, 1.825e-13, 3.739e-12] at [16, 64, 256] points, exposing
    // the roundoff floor rather than sample-grid convergence. Start at the
    // fitter's minimum order+2 samples and retain the original finest grid.
    const ORDER: usize = 6;
    let mut errors = Vec::new();
    for density in [8, 16, 256] {
        let mut s = samples(density, 0.0, piston);
        s.source_id = "analytic:baffled-piston:Bessel-Struve:a=.05:rho=1.2:c=343".into();
        let (observer, _) =
            BroadbandModalObserver::build(&model(1), &s, controls(ORDER, 1, 10.0), RATE).unwrap();
        let error = maximum_operator_error(&observer, &s);
        errors.push(error);
        println!(
            "{{\"suite\":\"broadband-modal\",\"axis\":\"fit-grid\",\"points\":{density},\"pole_count\":{ORDER},\"producer_id\":\"{}\",\"held_out_max_rel\":{error},\"tolerance\":1e-6}}",
            s.source_id
        );
    }
    assert!(errors.windows(2).all(|e| e[1] < e[0]), "{errors:?}");
    assert!(errors[2] <= 1e-6, "{errors:?}");
}

#[test]
fn integer_delay_is_exact_fractional_delay_converges() {
    let dt = 1.0 / f64::from(RATE);
    let h = |w: f64| C64::new(1.0, -w * RADIUS / SOUND_SPEED).recip();
    let undelayed = samples(64, 0.0, h);
    let integer = samples(64, 137.0 * dt, |w| h(w) * phase(w * 137.0 * dt));
    let (mut a, _) =
        BroadbandModalObserver::build(&model(1), &undelayed, controls(1, 1, 1e-6), RATE).unwrap();
    let (mut b, _) =
        BroadbandModalObserver::build(&model(1), &integer, controls(1, 1, 1e-6), RATE).unwrap();
    let mut history = Vec::new();
    for n in 0..1024 {
        let f = det::cos(n as f64 * 0.11);
        history.push(a.step(&[f]).unwrap());
        let expected = if n < 137 { 0.0 } else { history[n - 137] };
        // Separate fits see slightly different delay-removal roundoff.
        assert!((b.step(&[f]).unwrap() - expected).abs() < 2e-12);
    }
    let s = samples(128, DELAY, sphere);
    let mut errors = Vec::new();
    for half in [8, 32, 128] {
        let (observer, report) =
            BroadbandModalObserver::build(&model(1), &s, controls(1, half, 10.0), RATE).unwrap();
        let error = maximum_operator_error(&observer, &s);
        errors.push(error);
        println!(
            "{{\"suite\":\"broadband-modal\",\"axis\":\"fractional-delay\",\"half_length\":{half},\"delay_samples\":{},\"held_out_max_rel\":{error},\"tolerance\":1e-6}}",
            report.delay_samples
        );
    }
    assert!(errors.windows(2).all(|e| e[1] < e[0]), "{errors:?}");
    assert!(errors[2] <= 1e-6, "{errors:?}");
}

#[test]
fn mode_cap_and_cap_plus_one() {
    let mut s = samples(2, 0.0, |_| C64::ZERO);
    s.per_mode = vec![vec![C64::ZERO; 2]; 4096];
    for (_, row) in &mut s.held_out {
        *row = vec![C64::ZERO; 4096];
    }
    let (_, report) =
        BroadbandModalObserver::build(&model(4096), &s, controls(128, 4096, 1e-6), RATE).unwrap();
    assert_eq!(report.input_batches, 16);
    s.per_mode.push(vec![C64::ZERO; 2]);
    assert!(matches!(
        BroadbandModalObserver::build(&model(1), &s, controls(1, 1, 1e-6), RATE),
        Err(BroadbandModalError::TooManyModes {
            found: 4097,
            cap: 4096
        })
    ));
    assert!(admit_modal_samples(4096, 512, 512).is_ok());
    assert!(matches!(
        admit_modal_samples(4096, 513, 512),
        Err(BroadbandModalError::Cap { .. })
    ));
    assert!(matches!(
        admit_modal_samples(4096, 512, 513),
        Err(BroadbandModalError::Cap { .. })
    ));
    let s = samples(3, 0.0, |_| C64::ZERO);
    assert!(
        BroadbandModalObserver::build(
            &model(1),
            &s,
            controls(MAX_MODAL_TRANSFER_ORDER + 1, 1, 1e-6),
            RATE
        )
        .is_err()
    );
    assert!(
        BroadbandModalObserver::build(
            &model(1),
            &s,
            controls(1, MAX_MODAL_DELAY_KERNEL_HALF_LEN + 1, 1e-6),
            RATE
        )
        .is_err()
    );
}

#[test]
fn held_out_budget_refuses() {
    let s = samples(128, DELAY, sphere);
    assert!(matches!(
        BroadbandModalObserver::build(&model(1), &s, controls(1, 8, 1e-12), RATE),
        Err(BroadbandModalError::HeldOut { .. })
    ));
    let mut bad = samples(64, 0.0, |w| C64::new(1.0, -w / 1000.0).recip());
    bad.held_out[0].1[0] = bad.held_out[0].1[0].scale(1.1);
    assert!(matches!(
        BroadbandModalObserver::build(&model(1), &bad, controls(1, 1, 1e-6), RATE),
        Err(BroadbandModalError::HeldOut { .. })
    ));
    bad.omega_rad_s[0] = core::f64::consts::PI * f64::from(RATE);
    assert!(BroadbandModalObserver::build(&model(1), &bad, controls(1, 1, 1e-6), RATE).is_err());
}
