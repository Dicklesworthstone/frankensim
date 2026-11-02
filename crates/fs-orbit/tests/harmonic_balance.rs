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
