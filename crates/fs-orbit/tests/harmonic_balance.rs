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
