//! G1/G3 harmonic-balance regressions against analytic forced responses.
//! These verify the retained Fourier balance, not continuum orbit existence.

use fs_math::c64::C64;
use fs_orbit::{HbAnchor, HbBudget, OrbitError, OrbitProblem, solve_hb};

/// x'' + 2 zeta omega_0 x' + omega_0^2 x = dc + F cos(q omega t).
struct ForcedOscillator {
    omega: f64,
    harmonic: usize,
    force: f64,
    dc: f64,
}

impl ForcedOscillator {
    const NATURAL: f64 = 2.3;
    const ZETA: f64 = 0.15;

    fn response(&self) -> C64 {
        let w = self.harmonic as f64 * self.omega;
        C64::new(0.5 * self.force, 0.0)
            / C64::new(
                Self::NATURAL * Self::NATURAL - w * w,
                2.0 * Self::ZETA * Self::NATURAL * w,
            )
    }
}

impl OrbitProblem for ForcedOscillator {
    fn dim(&self) -> usize {
        2
    }

    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        out[0] = x[1];
        out[1] = self.dc + self.force * fs_math::det::cos(self.harmonic as f64 * self.omega * t)
            - 2.0 * Self::ZETA * Self::NATURAL * x[1]
            - Self::NATURAL * Self::NATURAL * x[0];
    }

    fn autonomous(&self) -> bool {
        self.force == 0.0
    }
}

fn assert_complex(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "actual {actual:?}, expected {expected:?}, tolerance {tolerance}"
    );
}

#[test]
fn higher_harmonic_forcing_matches_displacement_and_velocity_phasors() {
    for harmonic in [2, 3] {
        let problem = ForcedOscillator {
            omega: 0.75,
            harmonic,
            force: 0.4,
            dc: 0.0,
        };
        for deflate_equilibrium in [false, true] {
            let orbit = solve_hb(
                &problem,
                HbAnchor::Forced {
                    omega: problem.omega,
                },
                0.1,
                &HbBudget {
                    harmonics: harmonic,
                    tolerance: 1.0e-12,
                    deflate_equilibrium,
                    ..HbBudget::default()
                },
            )
            .expect("a response at a retained higher harmonic is not equilibrium");
            let x = problem.response();
            let v = C64::new(0.0, harmonic as f64 * problem.omega) * x;
            assert_complex(orbit.coeffs[0][harmonic], x, 1.0e-10);
            assert_complex(orbit.coeffs[1][harmonic], v, 1.0e-10);
            for n in 0..=harmonic {
                if n != harmonic {
                    assert!(orbit.coeffs[0][n].abs() < 1.0e-10);
                    assert!(orbit.coeffs[1][n].abs() < 1.0e-10);
                }
            }
        }
    }
}

#[test]
fn small_forced_response_has_no_absolute_amplitude_floor() {
    let problem = ForcedOscillator {
        omega: 0.75,
        harmonic: 2,
        force: 1.0e-11,
        dc: 0.0,
    };
    let orbit = solve_hb(
        &problem,
        HbAnchor::Forced {
            omega: problem.omega,
        },
        0.0,
        &HbBudget {
            harmonics: 2,
            tolerance: 1.0e-22,
            deflate_equilibrium: true,
            ..HbBudget::default()
        },
    )
    .expect("a resolved small forced response must be admitted");
    assert_complex(orbit.coeffs[0][2], problem.response(), 1.0e-20);
    assert!(orbit.coeffs[0][2].abs() > 0.0);
}

#[test]
fn prescribed_period_allows_static_and_zero_forced_responses() {
    for dc in [0.0, 0.7] {
        let problem = ForcedOscillator {
            omega: 1.0,
            harmonic: 1,
            force: 0.0,
            dc,
        };
        let orbit = solve_hb(
            &problem,
            HbAnchor::Forced {
                omega: problem.omega,
            },
            0.0,
            &HbBudget {
                harmonics: 2,
                deflate_equilibrium: true,
                ..HbBudget::default()
            },
        )
        .expect("a prescribed-period balance can be constant");
        let displacement = dc / ForcedOscillator::NATURAL.powi(2);
        assert_complex(orbit.coeffs[0][0], C64::new(displacement, 0.0), 1.0e-10);
        for state in &orbit.coeffs {
            for coefficient in &state[1..] {
                assert!(coefficient.abs() < 1.0e-10);
            }
        }
    }
}

#[test]
fn unknown_period_still_refuses_equilibrium() {
    let problem = ForcedOscillator {
        omega: 1.0,
        harmonic: 1,
        force: 0.0,
        dc: 0.0,
    };
    assert!(matches!(
        solve_hb(
            &problem,
            HbAnchor::Autonomous { omega_guess: 1.0 },
            0.0,
            &HbBudget::default(),
        ),
        Err(OrbitError::TrivialCollapse)
    ));
}

/// Nonlinear algebraic port x + x^7 = F cos(t). For N=1, the exact
/// retained Galerkin equation is z + 35 |z|^6 z = F/2, z = X_1.
/// Its solution at F=99/64 is x=cos(t). Unretained residual harmonics
/// remain: this fixture tests quadrature, not a continuum orbit claim.
struct SeventhOrder {
    force: f64,
}

impl OrbitProblem for SeventhOrder {
    fn dim(&self) -> usize {
        1
    }

    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        out[0] = self.force * fs_math::det::cos(t) - x[0].powi(7);
    }

    fn port(&self, _s: C64) -> Vec<C64> {
        vec![C64::new(1.0, 0.0)]
    }

    fn autonomous(&self) -> bool {
        false
    }
}

impl fs_orbit::ContinuableProblem for SeventhOrder {
    fn set_parameter(&mut self, lambda: f64) {
        self.force = lambda;
    }
}

#[test]
fn nonlinear_aliasing_is_refined_before_a_balance_is_accepted() {
    let problem = SeventhOrder { force: 99.0 / 64.0 };
    let orbit = solve_hb(
        &problem,
        HbAnchor::Forced { omega: 1.0 },
        1.0,
        &HbBudget {
            harmonics: 1,
            tolerance: 1.0e-12,
            ..HbBudget::default()
        },
    )
    .expect("the seventh-degree nonlinearity needs a finer AFT quadrature");

    // On the old eight-point grid, cos(7 theta) aliases into cos(theta)
    // and the amplitude is approximately 0.9968112455659488.
    assert_complex(orbit.coeffs[0][1], C64::new(0.5, 0.0), 1.0e-10);
    assert!(orbit.aft_samples > 8);
    assert!(orbit.sampling_defect < 1.0e-10);

    // Independent odd, shifted grid and platform trig: not the solver's
    // basis table, sample count, or deterministic trigonometric arithmetic.
    let mut g0 = 0.0;
    let mut g1 = C64::new(0.0, 0.0);
    for j in 0..65 {
        let theta = core::f64::consts::TAU * (f64::from(j) + 0.37) / 65.0;
        let x = orbit.sample(0, theta);
        let g = problem.force * theta.cos() - x.powi(7);
        g0 += g / 65.0;
        g1 = g1 + C64::new(theta.cos(), -theta.sin()).scale(g / 65.0);
    }
    assert!((orbit.coeffs[0][0].re - g0).abs() < 1.0e-10);
    assert_complex(orbit.coeffs[0][1], g1, 1.0e-10);
}

#[test]
fn continuation_retains_the_verified_residual_instead_of_zero() {
    let mut problem = SeventhOrder { force: 1.0 };
    let points = fs_orbit::continue_branch(
        &mut problem,
        |_| 1.0,
        1.0,
        1.3,
        0.8,
        &fs_orbit::ContinuationBudget {
            max_steps: 40,
            initial_step: 0.03,
            min_step: 1.0e-6,
            hb: HbBudget {
                harmonics: 1,
                max_newton: 30,
                tolerance: 1.0e-7,
                ..HbBudget::default()
            },
        },
    )
    .expect("smooth scalar force-response branch");
    assert!(points.len() > 2);
    let mut nonzero_residual = false;
    for point in &points {
        let z = point.orbit.coeffs[0][1];
        let analytic =
            (z + z.scale(35.0 * z.abs().powi(6)) - C64::new(0.5 * point.lambda, 0.0)).abs();
        assert!((point.orbit.residual - analytic).abs() < 5.0e-12);
        assert_eq!(
            point.orbit.residual_trace.last().copied(),
            Some(point.orbit.residual),
        );
        assert!(point.orbit.aft_samples > 8);
        assert!(point.orbit.sampling_defect < 1.0e-7);
        nonzero_residual |= point.orbit.residual > 1.0e-12;
    }
    assert!(
        nonzero_residual,
        "this fixture must expose a measurable final defect"
    );
}

/// Smooth unresolved forcing whose retained aliases differ on every grid
/// of the bounded N=1 ladder (8,17,35,71,143). Finite-grid agreement is
/// not a universal proof; persistent disagreement must be a named refusal.
struct UnresolvedForcing;

impl OrbitProblem for UnresolvedForcing {
    fn dim(&self) -> usize {
        1
    }

    fn island(&self, t: f64, _x: &[f64], out: &mut [f64]) {
        out[0] = [7, 16, 34, 70, 142]
            .iter()
            .enumerate()
            .map(|(i, &n)| (i + 1) as f64 * fs_math::det::cos(f64::from(n) * t))
            .sum();
    }

    fn port(&self, _s: C64) -> Vec<C64> {
        vec![C64::new(1.0, 0.0)]
    }

    fn autonomous(&self) -> bool {
        false
    }
}

#[test]
fn persistent_sampling_disagreement_has_a_bounded_named_refusal() {
    let result = solve_hb(
        &UnresolvedForcing,
        HbAnchor::Forced { omega: 1.0 },
        0.0,
        &HbBudget {
            harmonics: 1,
            ..HbBudget::default()
        },
    );
    match result {
        Err(OrbitError::SamplingUnresolved {
            samples,
            residual,
            defect,
        }) => {
            assert_eq!(samples, 143);
            assert!(residual > 0.1);
            assert!(defect > 0.1);
        }
        other => panic!("expected unresolved sampling, got {other:?}"),
    }
}

// The triangular nonlinear ODE has an exact finite Fourier response.
// X1=(1-i)/4, Y1=-3i/32, Y3=(-2+i)/160; the cubic creates the
// third harmonic even when the retained N=1 equations have zero defect.
struct CubicCascade {
    equation_scale: f64,
    dc: f64,
    maximum_port_harmonic: f64,
}

impl OrbitProblem for CubicCascade {
    fn dim(&self) -> usize {
        2
    }
    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        out[0] = fs_math::det::cos(t) - x[0];
        out[1] = self.equation_scale * (self.dc + x[0].powi(3) - x[1]);
    }
    fn port(&self, s: C64) -> Vec<C64> {
        assert!(
            s.im <= self.maximum_port_harmonic,
            "unretained port was queried"
        );
        vec![
            s,
            C64::new(0.0, 0.0),
            C64::new(0.0, 0.0),
            s.scale(self.equation_scale),
        ]
    }
    fn autonomous(&self) -> bool {
        false
    }
}

fn order_budget(max_harmonics: usize) -> fs_orbit::HbRefinementBudget {
    fs_orbit::HbRefinementBudget {
        max_harmonics,
        residual_tolerances: vec![1.0e-9; 2],
        state_tolerances: vec![1.0e-9; 2],
        frequency_tolerance: 1.0e-10,
    }
}

#[test]
fn adaptive_order_resolves_generated_harmonics_and_compares_actual_solutions() {
    let problem = CubicCascade {
        equation_scale: 1.0,
        dc: 0.0,
        maximum_port_harmonic: 7.0,
    };
    let result = fs_orbit::solve_hb_adaptive(
        &problem,
        HbAnchor::Forced { omega: 1.0 },
        0.5,
        &HbBudget {
            harmonics: 1,
            tolerance: 1.0e-13,
            ..HbBudget::default()
        },
        &order_budget(7),
    )
    .expect("cubic cascade resolves after adding its generated third harmonic");
    let orders: Vec<_> = result
        .refinement_trace
        .iter()
        .map(|check| check.harmonics)
        .collect();
    assert_eq!(orders, [1, 3, 7]);
    let first = result.refinement_trace[0];
    assert!(first.retained_residual_ratio < 1.0);
    assert!((first.omitted_residual_ratio * 1.0e-9 - 0.0625).abs() < 1.0e-11);
    assert_eq!(first.state_change_ratio, None);
    let middle = result.refinement_trace[1];
    assert!(middle.omitted_residual_ratio < 1.0);
    assert!(middle.state_change_ratio.expect("actual larger solve") > 1.0e6);
    let last = result.refinement_trace.last().unwrap();
    assert!(last.state_change_ratio.unwrap() < 1.0);
    assert!(last.sampling_defect_ratio < 1.0);
    assert_complex(result.orbit.coeffs[0][1], C64::new(0.25, -0.25), 1.0e-10);
    assert_complex(result.orbit.coeffs[1][1], C64::new(0.0, -0.09375), 1.0e-10);
    assert_complex(
        result.orbit.coeffs[1][3],
        C64::new(-0.0125, 0.00625),
        1.0e-10,
    );
    for n in [0, 2, 4, 5, 6, 7] {
        assert!(result.orbit.coeffs[0][n].abs() < 1.0e-10);
        assert!(result.orbit.coeffs[1][n].abs() < 1.0e-10);
    }
}

#[test]
fn adaptive_order_cap_checks_beyond_the_cap_without_dc_or_unit_masking() {
    for (scale, dc) in [(1.0, 0.0), (1.0e-6, 1000.0)] {
        let problem = CubicCascade {
            equation_scale: scale,
            dc,
            maximum_port_harmonic: 2.0,
        };
        let mut refinement = order_budget(2);
        refinement.residual_tolerances[1] *= scale;
        // Use the exact retained response so this probes harmonic omission,
        // independently of the inner solver's separate global residual scale.
        let result = fs_orbit::solve_hb_adaptive_seeded(
            &problem,
            HbAnchor::Forced { omega: 1.0 },
            vec![0.0, dc, 0.25, 0.0, -0.25, -0.09375],
            &HbBudget {
                harmonics: 1,
                tolerance: 1.0e-14,
                ..HbBudget::default()
            },
            &refinement,
        );
        match result {
            Err(OrbitError::HarmonicsUnresolved { check }) => {
                assert_eq!(check.harmonics, 2);
                assert_eq!(check.checked_harmonics, 5);
                assert!(check.retained_residual_ratio < 1.0);
                assert!((check.omitted_residual_ratio * 1.0e-9 - 0.0625).abs() < 1.0e-8);
                assert!(check.state_change_ratio.is_some());
            }
            other => panic!("unresolved cubic harmonic must survive cap/scaling: {other:?}"),
        }
    }
}

struct LinearBackbone;
impl OrbitProblem for LinearBackbone {
    fn dim(&self) -> usize {
        2
    }
    fn island(&self, _t: f64, x: &[f64], out: &mut [f64]) {
        out[0] = x[1];
        out[1] = -4.0 * x[0];
    }
    fn autonomous(&self) -> bool {
        true
    }
}

#[test]
fn adaptive_order_relocates_the_solved_frequency_in_the_warm_seed() {
    let result = fs_orbit::solve_hb_adaptive(
        &LinearBackbone,
        HbAnchor::Backbone {
            amplitude: 1.0,
            omega_guess: 1.7,
        },
        1.0,
        &HbBudget {
            harmonics: 1,
            tolerance: 1.0e-12,
            ..HbBudget::default()
        },
        &order_budget(3),
    )
    .expect("known-frequency backbone survives a larger canonical packing");
    assert_eq!(result.refinement_trace.len(), 2);
    assert!((result.orbit.omega - 2.0).abs() < 1.0e-10);
    assert!((result.orbit.first_harmonic_amplitude(0) - 1.0).abs() < 1.0e-10);
    assert!(
        result.refinement_trace[1]
            .relative_frequency_change
            .unwrap()
            < 1.0e-10
    );
    assert!(result.refinement_trace[1].state_change_ratio.unwrap() < 1.0);
}

#[test]
fn adaptive_order_distinguishes_loose_newton_from_missing_harmonics() {
    struct UnitEquilibrium;
    impl OrbitProblem for UnitEquilibrium {
        fn dim(&self) -> usize {
            1
        }
        fn island(&self, _: f64, _: &[f64], out: &mut [f64]) {
            out[0] = 1.0;
        }
        fn port(&self, _: C64) -> Vec<C64> {
            vec![C64::new(1.0, 0.0)]
        }
        fn autonomous(&self) -> bool {
            false
        }
    }
    let result = fs_orbit::solve_hb_adaptive(
        &UnitEquilibrium,
        HbAnchor::Forced { omega: 1.0 },
        0.0,
        &HbBudget {
            harmonics: 1,
            tolerance: 2.0,
            ..HbBudget::default()
        },
        &fs_orbit::HbRefinementBudget {
            max_harmonics: 7,
            residual_tolerances: vec![1.0e-9],
            state_tolerances: vec![1.0e-9],
            frequency_tolerance: 1.0e-9,
        },
    );
    match result {
        Err(OrbitError::HarmonicsUnresolved { check }) => {
            assert_eq!(
                check.harmonics, 1,
                "do not spend higher orders on a loose Newton solve"
            );
            assert!(check.retained_residual_ratio > 1.0e8);
            assert!(check.omitted_residual_ratio < 1.0);
            assert_eq!(check.state_change_ratio, None);
        }
        other => panic!("physical retained residual must be checked: {other:?}"),
    }
}

#[test]
fn adaptive_order_admission_precedes_model_evaluation() {
    struct NoEvaluation;
    impl OrbitProblem for NoEvaluation {
        fn dim(&self) -> usize {
            2
        }
        fn island(&self, _: f64, _: &[f64], _: &mut [f64]) {
            panic!("invalid budget evaluated");
        }
        fn autonomous(&self) -> bool {
            true
        }
    }
    let anchor = HbAnchor::Autonomous { omega_guess: 1.0 };
    let hb = HbBudget {
        harmonics: 1,
        ..HbBudget::default()
    };
    let refinement = order_budget(3);
    for omega in [0.0, -1.0, f64::NAN] {
        let mut seed = vec![0.0; 7];
        seed[6] = omega;
        assert!(matches!(
            fs_orbit::solve_hb_adaptive_seeded(&NoEvaluation, anchor, seed, &hb, &refinement),
            Err(OrbitError::BadParameter { .. }),
        ));
    }
    for bad in [
        fs_orbit::HbRefinementBudget {
            max_harmonics: 1,
            ..refinement.clone()
        },
        fs_orbit::HbRefinementBudget {
            residual_tolerances: vec![1.0e-9],
            ..refinement.clone()
        },
        fs_orbit::HbRefinementBudget {
            state_tolerances: vec![0.0; 2],
            ..refinement.clone()
        },
    ] {
        assert!(matches!(
            fs_orbit::solve_hb_adaptive(&NoEvaluation, anchor, 1.0, &hb, &bad),
            Err(OrbitError::BadParameter { .. }),
        ));
    }
    for bad_hb in [
        HbBudget {
            tolerance: f64::NAN,
            ..hb
        },
        HbBudget {
            max_newton: 0,
            ..hb
        },
    ] {
        assert!(matches!(
            fs_orbit::solve_hb_adaptive(&NoEvaluation, anchor, 1.0, &bad_hb, &refinement),
            Err(OrbitError::BadParameter { .. }),
        ));
    }
}

#[test]
fn adaptive_retained_and_omitted_residuals_share_one_physical_budget() {
    struct SplitResidual;
    impl OrbitProblem for SplitResidual {
        fn dim(&self) -> usize {
            1
        }
        fn island(&self, t: f64, _: &[f64], out: &mut [f64]) {
            out[0] = 0.8 + 0.8 * fs_math::det::sqrt(2.0) * fs_math::det::cos(5.0 * t);
        }
        fn port(&self, _: C64) -> Vec<C64> {
            vec![C64::new(1.0, 0.0)]
        }
        fn autonomous(&self) -> bool {
            false
        }
    }
    // A deliberately loose inner solve returns the same zero waveform at N=1
    // and N=3. At N=3 each spectral band has RMS 0.8, but their combined
    // equation defect sqrt(0.8^2 + 0.8^2) exceeds the one physical budget.
    let result = fs_orbit::solve_hb_adaptive(
        &SplitResidual,
        HbAnchor::Forced { omega: 1.0 },
        0.0,
        &HbBudget {
            harmonics: 1,
            tolerance: 2.0,
            ..HbBudget::default()
        },
        &fs_orbit::HbRefinementBudget {
            max_harmonics: 3,
            residual_tolerances: vec![1.0],
            state_tolerances: vec![1.0],
            frequency_tolerance: 1.0e-9,
        },
    );
    match result {
        Err(OrbitError::HarmonicsUnresolved { check }) => {
            assert_eq!(check.harmonics, 3);
            assert!((check.retained_residual_ratio - 0.8).abs() < 1e-12);
            assert!((check.omitted_residual_ratio - 0.8).abs() < 1e-12);
            assert!((check.total_residual_ratio - fs_math::det::sqrt(1.28)).abs() < 1e-12);
            assert!(check.total_residual_ratio > 1.0);
            assert_eq!(check.state_change_ratio, Some(0.0));
        }
        other => panic!("the full measured equation must meet its single budget: {other:?}"),
    }
}
