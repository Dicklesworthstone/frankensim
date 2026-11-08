use super::*;
use crate::{AftWorkspace, IslandStateJacobian, OrbitProblem, Packing, hb_residual};
use fs_math::det;

#[test]
fn exact_even_and_odd_grids_preserve_complex_harmonics_and_centered_dc() {
    for (samples, order) in [(256, 32), (512, 64), (1025, 256), (4097, 256)] {
        let pack = Packing {
            d: 2,
            n_harm: order,
            omega_slot: None,
        };
        let mut aft = AftWorkspace::with_samples(&pack, samples, 0).unwrap();
        assert!(aft.transform.is_some());
        let coefficients = |k: usize, n: usize| {
            let phase = (n * (k + 1)) as f64 * 0.371;
            let (s, c) = det::sin_cos(phase);
            C64::new(c, s).scale(0.2 / (n * n) as f64)
        };
        let mut actual = vec![0.0; 2 * samples];
        reconstruct(
            &aft.basis,
            aft.transform.as_mut(),
            samples,
            2,
            order,
            coefficients,
            |j, k, value| actual[j * 2 + k] = value,
        );
        let mut expected = vec![0.0; 2 * samples];
        reconstruct(
            &aft.basis,
            None,
            samples,
            2,
            order,
            coefficients,
            |j, k, value| expected[j * 2 + k] = value,
        );
        for (&a, &b) in actual.iter().zip(&expected) {
            assert!((a - b).abs() < 1e-14, "M={samples}, N={order}: {a} vs {b}");
        }
        aft.gs.copy_from_slice(&actual);
        aft.center_island(2);
        aft.project_island(&pack);
        assert!(aft.transform.as_ref().unwrap().projected);
        for k in 0..2 {
            assert!(aft.island_coefficient(&pack, 0, k).abs() < 2e-15);
            for n in 1..=order {
                assert!((aft.island_coefficient(&pack, n, k) - coefficients(k, n)).abs() < 2e-15);
            }
        }
        // A large constant source must remain exactly constant in every bin;
        // reusing a prepared plan must not leave any previous AC projection.
        aft.gs.fill(1e250);
        aft.center_island(2);
        aft.project_island(&pack);
        for k in 0..2 {
            assert_eq!(aft.island_coefficient(&pack, 0, k), C64::from_re(1e250));
            for n in 1..=order {
                assert_eq!(aft.island_coefficient(&pack, n, k), C64::ZERO);
            }
        }
    }
    assert!(enabled(256, 32));
    assert!(!enabled(513, 64));
    assert!(enabled(4097, 256));
    assert_eq!(payload(513, 64, 2).unwrap(), 0);
}

struct Cubic;

#[test]
fn finite_extreme_source_survives_unrepresentable_fft_intermediate_dc() {
    let pack = Packing {
        d: 1,
        n_harm: 32,
        omega_slot: None,
    };
    let mut aft = AftWorkspace::with_samples(&pack, 256, 0).unwrap();
    aft.gs.fill(-1e308);
    aft.gs[0] = 1e308;
    aft.center_island(1);
    aft.project_island(&pack);
    // The deviations' mean alone is nearly -2e308. Restoring the positive
    // reference in the original compensated owner gives a finite physical DC.
    assert!(!aft.transform.as_ref().unwrap().projected);
    let ac = (1e308 / 256.0) * 2.0;
    assert!((aft.island_coefficient(&pack, 0, 0).re / (ac - 1e308) - 1.0).abs() < 4e-15);
    for n in 1..=32 {
        let value = aft.island_coefficient(&pack, n, 0);
        assert!((value.re / ac - 1.0).abs() < 5e-14);
        assert!((value.im / ac).abs() < 5e-14);
    }
}

impl OrbitProblem for Cubic {
    fn dim(&self) -> usize {
        1
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, _: f64, x: &[f64], out: &mut [f64]) {
        out[0] = -x[0] - 0.3 * x[0].powi(3);
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        Some(|_, _, x, out| {
            out[0] = -1.0 - 0.9 * x[0] * x[0];
            true
        })
    }
}

fn signed_coefficient(u: &[f64], pack: &Packing, n: isize) -> C64 {
    let order = n.unsigned_abs();
    if order > pack.n_harm {
        return C64::ZERO;
    }
    let c = C64::new(
        u[pack.re_index(order, 0)],
        if order == 0 {
            0.0
        } else {
            u[pack.im_index(order, 0)]
        },
    );
    if n < 0 { c.conj() } else { c }
}

#[test]
fn high_order_nonlinear_balance_and_tangent_match_full_convolution() {
    let pack = Packing {
        d: 1,
        n_harm: 64,
        omega_slot: None,
    };
    let mut u = vec![0.0; pack.coeff_len()];
    u[0] = 0.13;
    for (n, re, im) in [(1, 0.2, -0.1), (17, -0.04, 0.02), (32, 0.015, 0.01)] {
        u[pack.re_index(n, 0)] = re;
        u[pack.im_index(n, 0)] = im;
    }
    let mut aft = AftWorkspace::new(&pack).unwrap();
    let mut actual = vec![0.0; pack.coeff_len()];
    hb_residual(&Cubic, &pack, &u, 1.3, &mut aft, &mut actual);
    assert!(aft.transform.is_some());
    let mut tangent = crate::tangent::Workspace::new(&pack).unwrap();
    tangent
        .prepare(
            &Cubic,
            Cubic.island_state_jacobian().unwrap(),
            &pack,
            1.3,
            &aft,
        )
        .unwrap();
    let direction: Vec<_> = (0..pack.coeff_len())
        .map(|i| 0.001 * det::sin(i as f64 * 0.41))
        .collect();
    for n in 0..=64_isize {
        let mut cubic = C64::ZERO;
        let mut derivative = C64::ZERO;
        for a in -32..=32 {
            for b in -32..=32 {
                let product = signed_coefficient(&u, &pack, a) * signed_coefficient(&u, &pack, b);
                cubic = cubic + product * signed_coefficient(&u, &pack, n - a - b);
                derivative =
                    derivative + product * signed_coefficient(&direction, &pack, n - a - b);
            }
        }
        let linear = C64::new(1.0, n as f64 * 1.3);
        let expected = linear * signed_coefficient(&u, &pack, n) + cubic.scale(0.3);
        let expected_derivative =
            linear * signed_coefficient(&direction, &pack, n) + derivative.scale(0.9);
        assert!((signed_coefficient(&actual, &pack, n) - expected).abs() < 2e-14);
        for (row, expected) in [
            (pack.re_index(n as usize, 0), expected_derivative.re),
            (
                if n == 0 {
                    0
                } else {
                    pack.im_index(n as usize, 0)
                },
                if n == 0 {
                    expected_derivative.re
                } else {
                    expected_derivative.im
                },
            ),
        ] {
            let measured: f64 = direction
                .iter()
                .enumerate()
                .map(|(col, &dx)| tangent.coefficient_entry(row, col, &pack) * dx)
                .sum();
            assert!(
                (measured - expected).abs() < 2e-14,
                "order {n}, row {row}: {measured} vs {expected}"
            );
        }
    }
    // The independent odd verification grid keeps its own physical clock and
    // chooses the direct route at this order; both assess the original source.
    let mut finer = aft.finer(&pack).unwrap();
    assert!(finer.transform.is_none());
    let mut checked = vec![0.0; pack.coeff_len()];
    hb_residual(&Cubic, &pack, &u, 1.3, &mut finer, &mut checked);
    assert!(
        actual
            .iter()
            .zip(&checked)
            .all(|(a, b)| (a - b).abs() < 2e-14)
    );
}

#[test]
fn adaptive_refinement_still_resolves_forcing_above_the_initial_fft_order() {
    use crate::{HbAnchor, HbBudget, HbRefinementBudget, solve_hb_adaptive};
    struct UpperTone;
    impl OrbitProblem for UpperTone {
        fn dim(&self) -> usize {
            1
        }
        fn autonomous(&self) -> bool {
            false
        }
        fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
            out[0] = det::cos(65.0 * t + 0.37) - x[0];
        }
        fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
            Some(|_, _, _, out| {
                out[0] = -1.0;
                true
            })
        }
    }
    let result = solve_hb_adaptive(
        &UpperTone,
        HbAnchor::Forced { omega: 1.0 },
        0.0,
        &HbBudget {
            harmonics: 32,
            max_newton: 5,
            tolerance: 1e-10,
            deflate_equilibrium: false,
        },
        &HbRefinementBudget {
            max_harmonics: 129,
            residual_tolerances: vec![1e-8],
            state_tolerances: vec![1e-9],
            frequency_tolerance: 1e-9,
        },
    )
    .unwrap();
    assert_eq!(result.orbit.coeffs[0].len(), 130);
    let (s, c) = det::sin_cos(0.37);
    let expected = C64::new(c, s).scale(0.5) / C64::new(1.0, 65.0);
    assert!((result.orbit.coeffs[0][65] - expected).abs() < 1e-12);
    let check = result.refinement_trace.last().unwrap();
    assert!(check.total_residual_ratio <= 1.0);
    assert!(check.omitted_residual_ratio <= 1.0);
    assert!(check.state_change_ratio.unwrap() <= 1.0);
}

#[test]
fn continuation_admits_base_fft_even_when_finest_odd_grid_stays_direct() {
    use crate::autonomous::{
        AutonomousContinuationBudget, AutonomousContinuationError, AutonomousScales,
        continue_autonomous_branch,
    };
    use crate::{ContinuableProblem, HbOrbit};
    use core::cell::Cell;
    struct Oscillator {
        calls: Cell<usize>,
    }
    impl OrbitProblem for Oscillator {
        fn dim(&self) -> usize {
            2
        }
        fn autonomous(&self) -> bool {
            true
        }
        fn island(&self, _: f64, x: &[f64], out: &mut [f64]) {
            self.calls.set(self.calls.get() + 1);
            out.copy_from_slice(&[x[1], -x[0]]);
        }
    }
    impl ContinuableProblem for Oscillator {
        fn set_parameter(&mut self, _: f64) {}
    }
    let mut problem = Oscillator {
        calls: Cell::new(0),
    };
    let mut coefficients = vec![vec![C64::ZERO; 33]; 2];
    coefficients[0][1] = C64::from_re(0.5);
    coefficients[1][1] = C64::new(0.0, 0.5);
    let seed = HbOrbit {
        omega: 1.0,
        coeffs: coefficients,
        residual: 0.0,
        residual_trace: vec![],
        aft_samples: 0,
        sampling_defect: 0.0,
    };
    let controls = AutonomousContinuationBudget {
        harmonics: 32,
        steps: 0,
        max_newton: 8,
        max_attempts_per_step: 4,
        step_sizes: [0.1, 0.01, 0.2],
        tolerance: 1e-8,
        minimum_oscillation: 1e-5,
        maximum_unknowns: 132,
        maximum_payload_bytes: 900_000,
        maximum_residual_evaluations: 1000,
        scales: AutonomousScales {
            states: vec![1.0; 2],
            equations: vec![1.0; 2],
            angular_frequency: 1.0,
            parameter: 1.0,
        },
    };
    // The previous direct-only allowance was 883,264 bytes. This cap admits
    // that representation but not the coexisting prepared base FFT payload.
    let error = continue_autonomous_branch(
        &mut problem,
        [(0.0, &seed), (0.1, &seed)],
        &controls,
        || false,
    )
    .unwrap_err();
    assert!(
        matches!(error, AutonomousContinuationError::Invalid(what) if what.contains("payload"))
    );
    assert_eq!(problem.calls.get(), 0);
}

#[test]
fn forced_parameter_response_admits_and_reuses_the_original_fft_grids() {
    use crate::sensitivity::{
        ForcedSensitivityControls, ForcedSensitivityDirection, ForcedSensitivityError,
        PreparedForcedSensitivity,
    };
    use crate::{HbOrbit, HbScales};
    use core::cell::Cell;
    struct Forced {
        calls: Cell<usize>,
    }
    impl OrbitProblem for Forced {
        fn dim(&self) -> usize {
            1
        }
        fn autonomous(&self) -> bool {
            false
        }
        fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
            self.calls.set(self.calls.get() + 1);
            out[0] = det::cos(17.0 * t + 0.37) - x[0];
        }
        fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
            Some(|_, _, _, out| {
                out[0] = -1.0;
                true
            })
        }
    }
    let problem = Forced {
        calls: Cell::new(0),
    };
    let (s, c) = det::sin_cos(0.37);
    let expected = C64::new(c, s).scale(0.5) / C64::new(1.0, 17.0);
    let mut coefficients = vec![vec![C64::ZERO; 33]];
    coefficients[0][17] = expected;
    let orbit = HbOrbit {
        omega: 1.0,
        coeffs: coefficients,
        residual: 0.0,
        residual_trace: vec![],
        aft_samples: 0,
        sampling_defect: 0.0,
    };
    let scales = HbScales {
        states: vec![1.0],
        equations: vec![1.0],
        angular_frequency: 1.0,
    };
    let mut controls = ForcedSensitivityControls {
        samples_per_period: 256,
        maximum_unknowns: 65,
        maximum_payload_bytes: 1_700_000,
        maximum_work: 1_000_000_000,
        maximum_base_residual: 1e-10,
        maximum_tangent_sampling_defect: 1e-10,
        maximum_condition_estimate: 1e6,
        maximum_relative_residual: 1e-12,
        maximum_derivative_residual: 1e-10,
    };
    // This admits the previous 1,696,192-byte direct representation, but
    // refuses the additional coarse FFT before either source or port evaluation.
    assert!(matches!(
        PreparedForcedSensitivity::new(&problem, &orbit, controls, &scales, || false),
        Err(ForcedSensitivityError::Budget {
            quantity: "numerical payload bytes",
            ..
        })
    ));
    assert_eq!(problem.calls.get(), 0);
    controls.maximum_payload_bytes = 4_000_000;
    let mut prepared =
        PreparedForcedSensitivity::new(&problem, &orbit, controls, &scales, || false).unwrap();
    let calls = problem.calls.get();
    for factor in [1.0, -0.3, 0.0] {
        let result = prepared
            .solve_island_direction(
                ForcedSensitivityDirection {
                    parameter: 1.0,
                    parameter_scale: 1.0,
                    angular_frequency_derivative: 0.0,
                },
                |t, _, out| {
                    out[0] = factor * det::cos(17.0 * t + 0.37);
                    true
                },
                || false,
            )
            .unwrap();
        assert!((result.coeffs[0][17] - expected.scale(factor)).abs() < 1e-13);
        assert!(
            result.coeffs[0]
                .iter()
                .enumerate()
                .all(|(n, c)| n == 17 || c.abs() < 1e-13)
        );
        assert!(result.derivative_residual < 1e-10);
        assert_eq!(problem.calls.get(), calls);
    }
}
