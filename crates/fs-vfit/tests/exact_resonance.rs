//! G0/G1 held-input resonance normalization, phase and physical port work.
use fs_math::c64::C64;
use fs_vfit::discretize::DelayedFilter;

const DT: f64 = 1.0 / 48_000.0;
const W: f64 = core::f64::consts::TAU * 1_000.0;

fn line(dt: f64, sigma: f64, residue: C64) -> DelayedFilter {
    DelayedFilter::from_impulse_response(dt, vec![0.0; 4])
        .unwrap()
        .with_exact_poles(0.0, &[(W, sigma, residue)])
        .unwrap()
}
fn exponential(z: C64) -> C64 {
    C64::new(z.im.cos(), z.im.sin()).scale(z.re.exp())
}
fn near(actual: C64, expected: C64, relative: f64) {
    let error = (actual - expected).abs();
    assert!(
        error <= relative * expected.abs().max(1e-20),
        "actual={actual:?}, expected={expected:?}, error={error:e}"
    );
}

#[test]
fn g1_held_force_matches_continuous_state_and_independent_mean_quadrature() {
    for dt in [DT, 1.0 / 8_000.0] {
        for sigma in [0.0, 40.0, 7_000.0] {
            for residue in [C64::from_re(3.0), C64::new(1.4, -0.6)] {
                let mut filter = line(dt, sigma, residue);
                let p = C64::new(-sigma, W);
                let input = 0.7;
                let state =
                    |t: f64| residue * (exponential(p.scale(t)) - C64::ONE) / p.scale(1.0 / input);
                for sample in 0..24 {
                    let t0 = sample as f64 * dt;
                    // Independent quadrature of the continuous trajectory.
                    let subdivisions = 512;
                    let mut mean = 0.0;
                    for i in 0..=subdivisions {
                        let weight = if i == 0 || i == subdivisions {
                            1.0
                        } else if i % 2 == 0 {
                            2.0
                        } else {
                            4.0
                        };
                        mean += weight * 2.0 * state(t0 + dt * i as f64 / subdivisions as f64).re;
                    }
                    mean /= 3.0 * subdivisions as f64;
                    let preview = filter.boundary_response().unwrap().incoming(input).unwrap();
                    let actual = filter.push(input).unwrap();
                    assert!((actual - preview).abs() < 1e-15);
                    assert!(
                        (actual - mean).abs() < 3e-12 * dt,
                        "mean mismatch: dt={dt}, sigma={sigma}, sample={sample}, {actual} != {mean}"
                    );
                    assert!(
                        (filter.exact_poles()[0].state() - state(t0 + dt)).abs()
                            <= 2e-12 * residue.abs() / p.abs()
                    );
                }
            }
        }
    }
}

#[test]
fn g0_lossless_resonance_stores_exact_held_input_work() {
    let mut filter = line(DT, 0.0, C64::from_re(2.0));
    let (mut total_work, mut work_scale) = (0.0, 0.0);
    for i in 0..4_000 {
        let input = if i % 11 < 4 {
            0.0
        } else {
            (i as f64 * 0.17).sin()
        };
        let work = DT * input * filter.push(input).unwrap();
        total_work += work;
        work_scale += work.abs();
        assert!((filter.exact_modal_energy() - total_work).abs() <= 2e-12 * work_scale.max(1e-20));
    }
}

#[test]
fn g1_harmonic_response_matches_driven_playback_for_complex_residue() {
    for frequency_hz in [900.0, 1_000.0, 1_250.0] {
        let omega = core::f64::consts::TAU * frequency_hz;
        let mut filter = line(DT, 40.0, C64::new(2.0, -0.4));
        let expected = filter.resonance_response_e_minus_iwt(omega).unwrap();
        near(
            filter.resonance_response_e_minus_iwt(-omega).unwrap(),
            expected.conj(),
            2e-13,
        );
        let mut observed = C64::ZERO;
        let warmup = 36_000;
        let count = 48_000;
        for i in 0..warmup + count {
            let phase = omega * DT * i as f64;
            let output = filter.push(phase.cos()).unwrap();
            if i >= warmup {
                observed = observed + C64::new(phase.cos(), phase.sin()).scale(output);
            }
        }
        near(observed.scale(2.0 / count as f64), expected, 2e-10);
    }
}

#[test]
fn g1_narrow_resonance_preserves_physical_gain_and_true_poles() {
    let sigma = 1e-12;
    let residue = C64::new(2.0, -0.4);
    let response = line(DT, sigma, residue)
        .resonance_response_e_minus_iwt(W)
        .unwrap();
    let half = 0.5 * W * DT;
    near(
        response.scale(sigma),
        residue.conj().scale((half.sin() / half).powi(2)),
        1e-12,
    );
    assert!(
        line(DT, 0.0, residue)
            .resonance_response_e_minus_iwt(W)
            .is_err()
    );
    assert_eq!(
        line(DT, 0.0, C64::ZERO)
            .resonance_response_e_minus_iwt(W)
            .unwrap(),
        C64::ZERO
    );
}

#[test]
fn g0_state_overflow_refuses_before_any_history_changes() {
    let mut filter = line(1e-6, 0.0, C64::new(0.0, 1e308));
    let before = format!("{filter:?}");
    assert!(filter.push(1e8).is_err());
    assert_eq!(before, format!("{filter:?}"));
    assert!(filter.push(0.1).unwrap().is_finite());
}
