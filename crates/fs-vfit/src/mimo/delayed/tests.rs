use super::*;
use super::super::{MimoModel, check_mimo_passivity};

const DT: f64 = 0.125;

fn limits() -> DelayedMimoLoadLimits {
    DelayedMimoLoadLimits {
        source: MimoLoadLimits {
            band_rad_s: [0.1, core::f64::consts::PI / DT],
            sampled_passivity_points: 65, max_states: 8, max_memory_bytes: 1 << 20,
        },
        max_delayed_states: 32, max_history_values: 128,
    }
}

fn model(delay: f64) -> DelayedMimoModel {
    DelayedMimoModel {
        delays: vec![0.0, delay, delay, 0.0],
        model: MimoModel {
            ports: 2,
            poles: vec![PoleTerm::Real { pole: -2.0, residue: 0.0 },
                PoleTerm::Pair { pole: C64::new(-3.0, 4.0), residue: C64::ZERO }],
            residues: vec![
                vec![C64::from_re(0.4), C64::from_re(0.05),
                    C64::from_re(0.05), C64::from_re(0.3)],
                vec![C64::new(0.2, 0.03), C64::new(0.04, -0.01),
                    C64::new(0.04, -0.01), C64::new(0.15, 0.02)],
            ],
            d: vec![5.0, 0.2, 0.2, 4.0], e: vec![0.0; 4], reciprocal: true,
        },
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:.17e} != {expected:.17e}");
}

#[test]
fn integer_and_fractional_first_arrivals_retain_their_exact_held_area() {
    // G1: an impulse held for one full interval arrives as two partial areas.
    // The positive diagonal dominates the delayed mutual resistance at all w.
    for (delay, first, weights) in [(2.5 * DT, 2, [0.5, 0.5]),
        (0.25 * DT, 0, [0.75, 0.25]), (3.0 * DT, 3, [1.0, 0.0])] {
        let source = DelayedMimoModel { delays: vec![0.0, delay, delay, 0.0],
            model: MimoModel { ports: 2, poles: vec![], residues: vec![],
                d: vec![1.0, 0.4, 0.4, 1.0], e: vec![0.0; 4], reciprocal: true } };
        let mut runtime = PreparedDelayedMimoLoad::new(&source, DT, limits()).unwrap();
        assert_eq!(runtime.reports().0.class, CertificateClass::GridOnly);
        assert_eq!(runtime.reports().1.class, CertificateClass::GridOnly);
        let mut sum = 0.0;
        for step in 0..8 {
            let input = [if step == 0 { 1.0 } else { 0.0 }, 0.0];
            let mut output = [0.0; 2];
            runtime.step_mean_into(&input, &mut output).unwrap();
            let expected = if step == first { 0.4 * weights[0] }
                else if step == first + 1 { 0.4 * weights[1] } else { 0.0 };
            assert_eq!(output[1], expected);
            sum += output[1];
        }
        assert_eq!(sum, 0.4);
    }
}

#[test]
fn zero_delays_reproduce_the_original_mean_zoh_operator() {
    // G3: only the entrywise state chart differs; the original operator stays.
    let source = model(0.0);
    let mut delayed = PreparedDelayedMimoLoad::new(&source, DT, limits()).unwrap();
    let mut original = source.model.zoh(DT).unwrap();
    for k in 0..80 {
        let input = [(k as f64 * 0.71).sin(), (k as f64 * 0.39).cos()];
        let expected = original.step_mean(&input);
        let mut actual = [0.0; 2];
        delayed.step_mean_into(&input, &mut actual).unwrap();
        for a in 0..2 { close(actual[a], expected[a], 2.0e-14); }
    }
    for fraction in [0.0, 0.13, 0.79, 1.0] {
        let omega = fraction * core::f64::consts::PI / DT;
        let actual = delayed.eval_mean_discrete_iw(omega).unwrap();
        let expected = original.eval_mean_discrete_iw(omega);
        for (a, e) in actual.iter().zip(expected) {
            assert!((*a - e).abs() < 2.0e-14 * (1.0 + e.abs()));
        }
    }
}

// Independently integrate each continuous unit-step response from zero to t.
// Combining these primitives at the actual arrival times integrates any held
// input history without the runtime's two-piece recurrence or phi functions.
fn step_primitive(source: &MimoModel, entry: usize, time: f64) -> f64 {
    if time <= 0.0 { return 0.0; }
    let mut integral = source.d[entry] * time;
    for (pole, residues) in source.poles.iter().zip(&source.residues) {
        let (p, pair) = match pole {
            PoleTerm::Real { pole, .. } => (C64::from_re(*pole), false),
            PoleTerm::Pair { pole, .. } => (*pole, true),
        };
        let argument = p.scale(time);
        let exp = C64::new(argument.im.cos(), argument.im.sin()).scale(argument.re.exp());
        let primitive = ((exp - C64::ONE) / p - C64::from_re(time)) / p;
        let z = residues[entry] * primitive;
        integral += if pair { 2.0 * z.re } else { z.re };
    }
    integral
}

#[test]
fn retarded_real_and_complex_poles_match_independent_continuous_integrals() {
    // G1: abrupt, signed held inputs; exact integer and fractional arrivals.
    for delay in [0.25 * DT, 1.375 * DT, 3.0 * DT] {
        let source = model(delay);
        let mut runtime = PreparedDelayedMimoLoad::new(&source, DT, limits()).unwrap();
        let mut past: Vec<[f64; 2]> = Vec::new();
        for step in 0..40 {
            let input = [((step * 17 % 11) as f64 - 5.0) / 5.0,
                ((step * 7 % 13) as f64 - 6.0) / 7.0];
            past.push(input);
            let mut history = [0.0; 2];
            runtime.mean_history_output_into(&mut history).unwrap();
            let mut actual = [0.0; 2];
            runtime.step_mean_into(&input, &mut actual).unwrap();
            for a in 0..2 {
                let mut expected = 0.0;
                for b in 0..2 {
                    let entry = 2 * a + b;
                    for k in 0..=step {
                        let change = past[k][b] - if k == 0 { 0.0 } else { past[k - 1][b] };
                        let start = (step - k) as f64 * DT - source.delays[entry];
                        expected += change * (step_primitive(&source.model, entry, start + DT)
                            - step_primitive(&source.model, entry, start)) / DT;
                    }
                }
                close(actual[a], expected, 8.0e-13);
                let decomposed = history[a] + runtime.mean_direct_term()[a * 2..a * 2 + 2]
                    .iter().zip(input).map(|(d, u)| d * u).sum::<f64>();
                assert_eq!(actual[a], decomposed);
            }
        }
    }
}

#[test]
fn actual_discrete_response_matches_warm_harmonic_steps() {
    // G1: the query must describe this sampler, including fractional arrivals.
    for fraction in [0.0, 0.13, 0.79, 0.97] {
        let omega = fraction * core::f64::consts::PI / DT;
        for driven in 0..2 {
            let mut runtime = PreparedDelayedMimoLoad::new(&model(1.375 * DT), DT, limits()).unwrap();
            let response = runtime.eval_mean_discrete_iw(omega).unwrap();
            for step in 0..600 {
                let phase = omega * step as f64 * DT;
                let mut input = [0.0; 2];
                input[driven] = phase.cos();
                let mut actual = [0.0; 2];
                runtime.step_mean_into(&input, &mut actual).unwrap();
                if step >= 500 {
                    let wave = C64::new(phase.cos(), phase.sin());
                    for a in 0..2 { close(actual[a], (response[a * 2 + driven] * wave).re, 3.0e-13); }
                }
            }
        }
    }
}

#[test]
fn rejected_and_abandoned_candidates_preserve_poles_and_delay_rings() {
    // G4: a network may preview and abandon a solved boundary without time travel.
    let mut runtime = PreparedDelayedMimoLoad::new(&model(1.375 * DT), DT, limits()).unwrap();
    let mut output = [0.0; 2];
    runtime.step_mean_into(&[0.2, -0.3], &mut output).unwrap();
    let before = runtime.clone();
    assert!(runtime.preview_state(&[0.9, 0.1]));
    assert_eq!(runtime.pole_states(), before.pole_states());
    assert_eq!(runtime.input_history(), before.input_history());
    assert_eq!(runtime.history_cursor(), before.history_cursor());
    let mut unchanged = [42.0; 2];
    for input in [&[f64::NAN, 0.0][..], &[f64::MAX, f64::MAX], &[1.0]] {
        assert!(runtime.step_mean_into(input, &mut unchanged).is_err());
        assert_eq!(unchanged, [42.0; 2]);
    }
    assert!(runtime.step_mean_into(&[0.0, 0.0], &mut unchanged[..1]).is_err());
    assert!(runtime.mean_history_output_into(&mut unchanged[..1]).is_err());
    assert_eq!(runtime.pole_states(), before.pole_states());
    assert_eq!(runtime.input_history(), before.input_history());
    assert_eq!(runtime.history_cursor(), before.history_cursor());
    let mut reference = before;
    for _ in 0..12 {
        let mut expected = [0.0; 2];
        runtime.step_mean_into(&[-0.7, 0.4], &mut output).unwrap();
        reference.step_mean_into(&[-0.7, 0.4], &mut expected).unwrap();
        assert_eq!(output, expected);
        assert_eq!(runtime.pole_states(), reference.pole_states());
        assert_eq!(runtime.input_history(), reference.input_history());
    }
}

#[test]
fn hidden_coordinate_overflow_refuses_before_publishing_a_finite_effort() {
    let source = DelayedMimoModel { delays: vec![0.0], model: MimoModel {
        ports: 1, poles: vec![PoleTerm::Real { pole: -1e-4, residue: 0.0 }],
        residues: vec![vec![C64::from_re(1e-100)]], d: vec![0.0], e: vec![0.0], reciprocal: true,
    } };
    let dt = 1e4;
    let mut bounds = limits();
    bounds.source.band_rad_s = [1e-9, core::f64::consts::PI / dt];
    let mut runtime = PreparedDelayedMimoLoad::new(&source, dt, bounds).unwrap();
    let mut output = [42.0];
    assert!((runtime.mean_direct_term()[0] * 1e307).is_finite());
    assert!(runtime.step_mean_into(&[1e307], &mut output).is_err());
    assert_eq!(output, [42.0]);
    assert_eq!(runtime.pole_states(), [C64::ZERO]);
    assert_eq!(runtime.input_history(), [0.0]);
    runtime.step_mean_into(&[1.0], &mut output).unwrap();
}

#[test]
fn admission_checks_restored_passivity_shapes_and_separate_delay_budgets() {
    // Delay-free [[1,1.2],[-1.2,1]] has Hermitian part 2I. Delaying only
    // the upper mutual path makes the restored load active near phase pi.
    let active = DelayedMimoModel { delays: vec![0.0, DT, 0.0, 0.0],
        model: MimoModel { ports: 2, poles: vec![], residues: vec![],
            d: vec![1.0, 1.2, -1.2, 1.0], e: vec![0.0; 4], reciprocal: false } };
    assert!(check_mimo_passivity(&active.model, (0.1, 30.0)).unwrap().passive);
    assert!(PreparedDelayedMimoLoad::new(&active, DT, limits()).is_err());
    // Conversely, the residual may be active while the restored operator is
    // passive. A quarter-cycle mutual delay suppresses this narrow bandpass's
    // resistive peak; requiring passivity of the residual would reject it.
    let beta = (100.0f64 - 0.25).sqrt();
    let restored = DelayedMimoModel {
        delays: vec![0.0, core::f64::consts::PI / 20.0,
            core::f64::consts::PI / 20.0, 0.0],
        model: MimoModel { ports: 2,
            poles: vec![PoleTerm::Pair { pole: C64::new(-0.5, beta), residue: C64::ZERO }],
            residues: vec![vec![C64::ZERO, C64::new(0.75, 0.375 / beta),
                C64::new(0.75, 0.375 / beta), C64::ZERO]],
            d: vec![1.0, 0.0, 0.0, 1.0], e: vec![0.0; 4], reciprocal: true },
    };
    assert!(!check_mimo_passivity(&restored.model, (0.1, 30.0)).unwrap().passive);
    assert!(PreparedDelayedMimoLoad::new(&restored, DT, limits()).is_ok());
    for delay in [-1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(PreparedDelayedMimoLoad::new(&model(delay), DT, limits()).is_err());
    }
    let source = model(1.375 * DT);
    for field in 0..4 {
        let mut bad = source.clone();
        match field {
            0 => { bad.delays.pop(); },
            1 => { bad.delays[1] = 0.0; },
            2 => { bad.model.e[0] = 1.0; },
            _ => { bad.model.residues[0][0].im = 1.0; },
        }
        assert!(PreparedDelayedMimoLoad::new(&bad, DT, limits()).is_err());
    }
    let mut bounds = limits();
    bounds.max_delayed_states = 11;
    assert!(PreparedDelayedMimoLoad::new(&source, DT, bounds).is_err());
    bounds = limits(); bounds.max_history_values = 3;
    assert!(PreparedDelayedMimoLoad::new(&source, DT, bounds).is_err());
    bounds.max_history_values = 4;
    assert!(PreparedDelayedMimoLoad::new(&source, DT, bounds).is_ok());
    bounds = limits(); bounds.source.max_memory_bytes = 0;
    assert!(PreparedDelayedMimoLoad::new(&source, DT, bounds).is_err());
    let runtime = PreparedDelayedMimoLoad::new(&source, DT, limits()).unwrap();
    for frequency in [-1.0, f64::NAN, core::f64::consts::PI / DT + 1.0] {
        assert!(runtime.eval_mean_discrete_iw(frequency).is_err());
    }
}
