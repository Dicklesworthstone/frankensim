//! G0/G1 held-input resonance normalization, phase and physical port work.
use fs_math::c64::C64;
use fs_vfit::PoleTerm;
use fs_vfit::discretize::DelayedFilter;
use fs_vfit::mimo::delayed::{DelayedMimoLoadLimits, PreparedDelayedMimoLoad};
use fs_vfit::mimo::load::MimoLoadLimits;
use fs_vfit::mimo::{DelayedMimoModel, MimoModel};

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

fn multiport(sigma: f64, omega: f64) -> MimoModel {
    // A positive-definite residue matrix: the mutual response is signed,
    // while the complete impedance is positive real at every frequency.
    MimoModel {
        ports: 2,
        poles: vec![PoleTerm::Pair {
            pole: C64::new(-sigma, omega),
            residue: C64::ZERO,
        }],
        residues: vec![vec![
            C64::from_re(2.0),
            C64::from_re(-0.5),
            C64::from_re(-0.5),
            C64::from_re(1.0),
        ]],
        d: vec![0.0; 4],
        e: vec![0.0; 4],
        reciprocal: true,
    }
}

#[test]
fn g1_multiport_narrow_resonance_retains_gain_phase_and_signed_participation() {
    // exp(-sigma*dt) rounds to 1 for the narrowest resonance. The source
    // damping must still determine its finite frequency response.
    for sigma in [1e-7, 1e-12, 1e-18] {
        let model = multiport(sigma, W);
        let runtime = model.zoh(DT).unwrap();
        let theta = W * DT;
        let sinc = (theta * 0.5).sin() / (theta * 0.5);
        let mean = runtime.eval_mean_discrete_iw(W);
        let start = runtime.eval_discrete_iw(W);
        let negative = runtime.eval_mean_discrete_iw(-W);
        for (index, &residue) in model.residues[0].iter().enumerate() {
            // Independent small-loss limits of the held-input Fourier series.
            near(mean[index].scale(sigma), residue.scale(sinc * sinc), 2e-10);
            near(
                start[index].scale(sigma),
                residue * C64::new((theta * 0.5).cos(), -(theta * 0.5).sin()).scale(sinc),
                2e-10,
            );
            near(negative[index], mean[index].conj(), 3e-14);
        }
        // The scalar acoustic owner observes exactly the same mean map,
        // with its explicitly opposite phasor convention.
        let scalar = line(DT, sigma, C64::from_re(2.0));
        near(
            mean[0],
            scalar.resonance_response_e_minus_iwt(W).unwrap().conj(),
            3e-14,
        );
    }
}

#[test]
fn g1_multiport_slow_relaxation_has_its_exact_finite_dc_compliance() {
    for rate in [1e-8, 1e-12, 1e-18] {
        let model = MimoModel {
            ports: 1,
            poles: vec![PoleTerm::Real {
                pole: -rate,
                residue: 0.0,
            }],
            residues: vec![vec![C64::from_re(rate * 7.0)]],
            d: vec![0.25],
            e: vec![0.0],
            reciprocal: true,
        };
        let runtime = model.zoh(DT).unwrap();
        for response in [
            runtime.eval_discrete_iw(0.0),
            runtime.eval_mean_discrete_iw(0.0),
        ] {
            near(response[0], C64::from_re(7.25), 4e-14);
        }
    }
}

#[test]
fn g0_narrow_passive_multiport_resonance_is_admitted_on_the_sampled_grid() {
    // Put the resonance exactly on a DC-to-Nyquist passivity sample. A
    // nonfinite response here used to refuse this positive-real physical load.
    let omega = core::f64::consts::PI / DT * 0.25;
    let source = DelayedMimoModel {
        model: multiport(1e-12, omega),
        delays: vec![0.0; 4],
    };
    let original = source.model.zoh(DT).unwrap();
    assert!(original.mean_sampled_passivity(17).passive);
    let limits = DelayedMimoLoadLimits {
        source: MimoLoadLimits {
            band_rad_s: [1.0, core::f64::consts::PI / DT],
            sampled_passivity_points: 17,
            max_states: 4,
            max_memory_bytes: 1 << 20,
        },
        max_delayed_states: 8,
        max_history_values: 4,
    };
    let delayed = PreparedDelayedMimoLoad::new(&source, DT, limits).unwrap();
    let actual = delayed.eval_mean_discrete_iw(omega).unwrap();
    let expected = original.eval_mean_discrete_iw(omega);
    for (actual, expected) in actual.into_iter().zip(expected) {
        near(actual, expected, 4e-14);
    }
    assert!(delayed.reports().1.passive);
}
