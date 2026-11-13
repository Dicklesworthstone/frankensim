use super::*;
use fs_math::{c64::C64, det};
use fs_phs::{RelaxationBranch, Storage};

fn waveform(omega: f64, coeffs: Vec<Vec<C64>>) -> HbOrbit {
    HbOrbit {
        omega,
        coeffs,
        residual_trace: vec![],
        residual: 0.0,
        aft_samples: 0,
        sampling_defect: 0.0,
    }
}

fn settings(scales: &[f64]) -> (PhsFloquetControls, VariationalShootControls) {
    (
        PhsFloquetControls {
            shooting: ShootBudget {
                steps_per_period: 256,
                max_newton: 20,
                tolerance: 1e-11,
            },
            state_tolerances: scales.iter().map(|s| 2e-6 * s).collect(),
            spectral_radius_tolerance: 1e-6,
            unit_circle_guard: 1e-3,
        },
        VariationalShootControls {
            state_scales: scales.to_vec(),
            period_scale: 1.0,
            minimum_autonomous_excursion: 0.01,
            maximum_dimension: 16,
            maximum_payload_bytes: 2_000_000,
            maximum_work: 10_000_000_000,
            maximum_period_maps: 64,
        },
    )
}

fn amplification(pole: C64, period: f64, steps: usize) -> C64 {
    // Independent RK4 amplification polynomial for a constant physical pole.
    let z = pole.scale(period / steps as f64);
    let z2 = z * z;
    let value =
        C64::ONE + z + z2.scale(0.5) + (z2 * z).scale(1.0 / 6.0) + (z2 * z2).scale(1.0 / 24.0);
    (0..steps).fold(C64::ONE, |product, _| product * value)
}

fn oscillator(scale: f64, damping: f64) -> (PortHamiltonian, PeriodicPortInput, HbOrbit) {
    let (mass, stiffness, omega) = (2.0, 18.0, 2.4);
    let system = fs_phs::mass_spring_damper(mass, stiffness, damping).unwrap();
    let force: Vec<_> = [C64::from_re(0.4), C64::new(0.2, -0.1), C64::new(0.03, 0.02)]
        .iter()
        .map(|z| z.scale(scale))
        .collect();
    let q: Vec<_> = force
        .iter()
        .enumerate()
        .map(|(n, &f)| {
            let w = n as f64 * omega;
            f / C64::new(stiffness - mass * w * w, damping * w)
        })
        .collect();
    let p = q
        .iter()
        .enumerate()
        .map(|(n, &q)| q * C64::new(0.0, mass * n as f64 * omega))
        .collect();
    let input = PeriodicPortInput::new(omega, vec![force]).unwrap();
    (system, input, waveform(omega, vec![q, p]))
}

#[test]
fn linear_resonance_stability_preserves_tiny_physical_scales_and_exact_discrete_maps() {
    let mut reference: Option<Vec<Vec<f64>>> = None;
    for scale in [1e-40, 1.0, 1e40] {
        let (system, input, candidate) = oscillator(scale, 0.8);
        let (controls, variational) = settings(&[0.1 * scale, 0.5 * scale]);
        let result = analyze_forced_orbit_stability_variational(
            &system,
            &input,
            &candidate,
            &controls,
            &variational,
            || false,
        )
        .unwrap();
        assert_eq!(result.classification, PhsNumericalStability::Stable);
        assert_eq!(
            result
                .levels
                .iter()
                .map(|l| l.steps_per_period)
                .collect::<Vec<_>>(),
            [256, 512, 1024]
        );
        assert!((result.dominant_growth_rate_per_s.unwrap() + 0.2).abs() < 1e-8);
        for (index, level) in result.levels.iter().enumerate() {
            let orbit = &level.orbit;
            assert_eq!(
                orbit.orbit.period.to_bits(),
                (core::f64::consts::TAU / 2.4).to_bits()
            );
            assert_eq!(orbit.state_scales, variational.state_scales);
            assert!(orbit.normalized_residual <= controls.shooting.tolerance);
            assert_eq!(orbit.monodromy.len(), 4);
            assert_eq!(
                orbit.rhs_evaluations,
                4 * level.steps_per_period * orbit.period_maps
            );
            // Line-search trials evaluate the original flow without a tangent;
            // accepted Newton maps and the final map retain the exact tangent.
            assert!(orbit.jacobian_evaluations > 0);
            assert!(orbit.jacobian_evaluations <= orbit.rhs_evaluations);
            assert_eq!(orbit.jacobian_evaluations % (4 * level.steps_per_period), 0);
            for sign in [-1.0, 1.0] {
                let expected = amplification(
                    C64::new(-0.2, sign * det::sqrt(8.96)),
                    orbit.orbit.period,
                    level.steps_per_period,
                );
                assert!(
                    orbit
                        .orbit
                        .multipliers
                        .iter()
                        .any(|&(re, im)| (C64::new(re, im) - expected).abs() < 2e-10)
                );
            }
            for (i, &actual) in orbit.orbit.x0.iter().enumerate() {
                assert!((actual - candidate.sample(i, 0.0)).abs() <= controls.state_tolerances[i]);
            }
            if let Some(reference) = &reference {
                for (actual, expected) in orbit.monodromy.iter().zip(&reference[index]) {
                    assert!((actual - expected).abs() < 2e-11);
                }
            }
        }
        reference = Some(
            result
                .levels
                .iter()
                .map(|l| l.orbit.monodromy.clone())
                .collect(),
        );
    }
}

#[test]
fn maxwell_resonance_keeps_internal_memory_and_all_three_decay_multipliers() {
    // s^3 + 6s^2 + 11s + 6 = (s+1)(s+2)(s+3), independently of the solver.
    let (k, kb, tau, w) = (1.5, 1.5_f64, 0.25, core::f64::consts::TAU);
    let system = fs_phs::mass_spring_damper(1.0, k, 2.0)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: kb,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let force = C64::new(1e-18, 3e-19);
    let s = C64::new(0.0, w);
    let arm = C64::ONE + s.scale(tau);
    let q = force / (C64::new(k - w * w, 2.0 * w) + (s.scale(tau) / arm).scale(kb));
    let candidate = waveform(
        w,
        vec![
            vec![C64::ZERO, q],
            vec![C64::ZERO, s * q],
            vec![C64::ZERO, q.scale(det::sqrt(kb)) / arm],
        ],
    );
    let input = PeriodicPortInput::new(w, vec![vec![C64::ZERO, force]]).unwrap();
    let (controls, variational) = settings(&[1e-19, 1e-18, 1e-19]);
    let result = analyze_forced_orbit_stability_variational(
        &system,
        &input,
        &candidate,
        &controls,
        &variational,
        || false,
    )
    .unwrap();
    assert_eq!(result.classification, PhsNumericalStability::Stable);
    for level in &result.levels {
        assert_eq!(level.orbit.monodromy.len(), 9);
        assert_eq!(level.orbit.orbit.x0.len(), 3);
        assert_eq!(level.orbit.orbit.multipliers.len(), 3);
        for pole in [-1.0, -2.0, -3.0] {
            let expected = amplification(C64::from_re(pole), 1.0, level.steps_per_period);
            assert!(
                level
                    .orbit
                    .orbit
                    .multipliers
                    .iter()
                    .any(|&(re, im)| (C64::new(re, im) - expected).abs() < 2e-10)
            );
        }
    }
}

struct CubicSpring;
impl Storage for CubicSpring {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * (x[0] * x[0] + x[1] * x[1]) + 0.25 * x[0].powi(4)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = x[0] + x[0].powi(3);
        out[1] = x[1];
    }
}

#[test]
fn nonlinear_resonance_distinguishes_parametric_growth_from_a_stable_forced_cycle() {
    let damping = 0.05;
    let system = PortHamiltonian::new(
        2,
        1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, damping],
        vec![0.0, 1.0],
        Box::new(CubicSpring),
    )
    .unwrap();
    let hessian = |x: &[f64], v: &[f64], out: &mut [f64]| {
        out[0] = (1.0 + 3.0 * x[0] * x[0]) * v[0];
        out[1] = v[1];
        true
    };
    for (w, expected) in [
        (det::sqrt(2.5), PhsNumericalStability::Unstable),
        (3.0, PhsNumericalStability::Stable),
    ] {
        // q=cos(wt), p=-w*sin(wt) solves the ORIGINAL cubic ODE exactly.
        // Its perturbation obeys u''+.05u'+[2.5+1.5cos(2wt)]u=0.
        // w=sqrt(2.5) lies within its primary parametric-instability tongue;
        // w=3 lies in the stable low Mathieu-a band.
        let candidate = waveform(
            w,
            vec![
                vec![C64::ZERO, C64::from_re(0.5), C64::ZERO, C64::ZERO],
                vec![C64::ZERO, C64::new(0.0, 0.5 * w), C64::ZERO, C64::ZERO],
            ],
        );
        let input = PeriodicPortInput::new(
            w,
            vec![vec![
                C64::ZERO,
                C64::new(0.5 * (1.75 - w * w), 0.5 * damping * w),
                C64::ZERO,
                C64::from_re(0.125),
            ]],
        )
        .unwrap();
        let (controls, variational) = settings(&[1.0, w]);
        assert_eq!(
            analyze_forced_orbit_stability_variational(
                &system,
                &input,
                &candidate,
                &controls,
                &variational,
                || false,
            )
            .unwrap_err(),
            PhsOrbitError::Variational(VariationalShootError::TangentUnavailable),
        );
        let result = analyze_forced_orbit_stability_variational_analytic(
            &system,
            &input,
            &candidate,
            &controls,
            &variational,
            &hessian,
            || false,
        )
        .unwrap();
        assert_eq!(result.classification, expected, "omega={w}: {result:?}");
        assert_eq!(
            result.dominant_growth_rate_per_s.unwrap() > 0.0,
            expected == PhsNumericalStability::Unstable
        );
        for level in &result.levels {
            // Independent midpoint matrix-exponential products at 2048,
            // 4096, 8192 and 16384 panels converge quadratically to these
            // values (Richardson limit); no RK4 or source shooting is used.
            let reference_radius = if expected == PhsNumericalStability::Unstable {
                2.301_406_300_465_6
            } else {
                0.948_987_286_154
            };
            assert!((level.spectral_radius - reference_radius).abs() < 2e-6);
            let map = &level.orbit.monodromy;
            // Liouville's exact determinant is independent of RK4 and the
            // eigenvalue owner: trace of the physical tangent is -damping.
            let determinant = map[0] * map[3] - map[1] * map[2];
            let exact = det::exp(-damping * level.orbit.orbit.period);
            assert!((determinant - exact).abs() < 2e-7);
        }
    }
}

#[test]
fn unresolved_time_grids_and_lossless_neutral_motion_do_not_gain_a_stability_claim() {
    let (system, input, candidate) = oscillator(1.0, 0.8);
    let (mut controls, variational) = settings(&[0.1, 0.5]);
    controls.shooting.steps_per_period = 32;
    controls.state_tolerances = vec![0.1; 2];
    controls.spectral_radius_tolerance = 1e-12;
    let result = analyze_forced_orbit_stability_variational(
        &system,
        &input,
        &candidate,
        &controls,
        &variational,
        || false,
    )
    .unwrap();
    assert_eq!(result.classification, PhsNumericalStability::Unresolved);
    let (system, input, candidate) = oscillator(1.0, 0.0);
    let (controls, variational) = settings(&[0.1, 0.5]);
    let result = analyze_forced_orbit_stability_variational(
        &system,
        &input,
        &candidate,
        &controls,
        &variational,
        || false,
    )
    .unwrap();
    assert_eq!(result.classification, PhsNumericalStability::Unresolved);
    assert!(
        result
            .levels
            .iter()
            .all(|l| l.orbit.orbit.multipliers.len() == 2)
    );
}

#[test]
fn cancellation_budgets_seed_disagreement_and_refusing_tangents_publish_no_assessment() {
    let (system, input, candidate) = oscillator(1.0, 0.8);
    let (controls, mut variational) = settings(&[0.1, 0.5]);
    let original = candidate.coeffs.clone();
    assert_eq!(
        analyze_forced_orbit_stability_variational(
            &system,
            &input,
            &candidate,
            &controls,
            &variational,
            || true,
        )
        .unwrap_err(),
        PhsOrbitError::Variational(VariationalShootError::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        analyze_forced_orbit_stability_variational(
            &system,
            &input,
            &candidate,
            &controls,
            &variational,
            || {
                polls += 1;
                polls == 40
            },
        )
        .unwrap_err(),
        PhsOrbitError::Variational(VariationalShootError::Cancelled)
    );
    variational.maximum_work = 1;
    assert!(matches!(
        analyze_forced_orbit_stability_variational(
            &system,
            &input,
            &candidate,
            &controls,
            &variational,
            || false,
        ),
        Err(PhsOrbitError::Variational(
            VariationalShootError::Budget { .. }
        ))
    ));
    variational.maximum_work = 10_000_000_000;
    let mut wrong = candidate.clone();
    wrong.coeffs[0][0].re += 0.01;
    assert!(matches!(
        analyze_forced_orbit_stability_variational(
            &system,
            &input,
            &wrong,
            &controls,
            &variational,
            || false,
        ),
        Err(PhsOrbitError::SeedDisagreement { .. })
    ));
    let refuse = |_: &[f64], _: &[f64], _: &mut [f64]| false;
    assert!(
        analyze_forced_orbit_stability_variational_analytic(
            &system,
            &input,
            &candidate,
            &controls,
            &variational,
            &refuse,
            || false,
        )
        .is_err()
    );
    assert_eq!(candidate.coeffs, original);
}
