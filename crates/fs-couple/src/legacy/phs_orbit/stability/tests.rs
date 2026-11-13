use super::*;
use crate::phs_orbit::{continue_forced_frequency_response, solve_forced_harmonic_balance};
use fs_orbit::{ContinuationBudget, HbBudget};
use fs_phs::{RelaxationBranch, Storage};

fn hb(harmonics: usize) -> HbBudget {
    HbBudget {
        harmonics,
        max_newton: 40,
        tolerance: 1e-11,
        deflate_equilibrium: false,
    }
}

fn controls(state_tolerances: Vec<f64>) -> PhsFloquetControls {
    PhsFloquetControls {
        shooting: ShootBudget {
            steps_per_period: 256,
            max_newton: 20,
            tolerance: 1e-11,
        },
        state_tolerances,
        spectral_radius_tolerance: 2e-5,
        unit_circle_guard: 1e-3,
    }
}

#[test]
fn forced_linear_resonator_matches_analytic_multipliers_and_authored_phase() {
    let (mass, stiffness, damping, omega) = (2.0, 18.0, 0.8, 2.4);
    let system = fs_phs::mass_spring_damper(mass, stiffness, damping).unwrap();
    let coefficients = vec![C64::from_re(0.4), C64::new(0.2, -0.1), C64::new(0.03, 0.02)];
    let input = PeriodicPortInput::new(omega, vec![coefficients.clone()]).unwrap();
    let candidate = solve_forced_harmonic_balance(&system, &input, 0.1, &hb(3)).unwrap();
    let result =
        analyze_forced_orbit_stability(&system, &input, &candidate, &controls(vec![1e-7, 1e-6]))
            .unwrap();
    let period = core::f64::consts::TAU / omega;
    let decay_rate = damping / (2.0 * mass);
    let damped_omega = det::sqrt(stiffness / mass - decay_rate * decay_rate);
    let magnitude = det::exp(-decay_rate * period);
    let multiplier = C64::new(
        magnitude * det::cos(damped_omega * period),
        magnitude * det::sin(damped_omega * period),
    );
    let mut expected = [coefficients[0].re / stiffness, 0.0];
    for (n, force) in coefficients.iter().enumerate().skip(1) {
        let w = n as f64 * omega;
        let q = *force / C64::new(stiffness - mass * w * w, damping * w);
        expected[0] += 2.0 * q.re;
        expected[1] += 2.0 * (q * C64::new(0.0, mass * w)).re;
    }
    assert_eq!(result.classification, PhsNumericalStability::Stable);
    assert_eq!(
        result
            .levels
            .iter()
            .map(|l| l.steps_per_period)
            .collect::<Vec<_>>(),
        vec![256, 512, 1024]
    );
    for level in &result.levels {
        assert!((level.orbit.x0[0] - expected[0]).abs() < 1e-7);
        assert!((level.orbit.x0[1] - expected[1]).abs() < 1e-6);
        assert_eq!(level.orbit.multipliers.len(), 2);
        for &(re, im) in &level.orbit.multipliers {
            let value = C64::new(re, im);
            assert!(
                (value - multiplier)
                    .abs()
                    .min((value - multiplier.conj()).abs())
                    < 2e-5
            );
        }
        assert!((level.spectral_radius - magnitude).abs() < 2e-5);
    }
    assert!((result.dominant_growth_rate_per_s.unwrap() + decay_rate).abs() < 2e-5);
}

#[test]
fn maxwell_memory_contributes_its_original_third_multiplier() {
    // The original characteristic polynomial is
    // s^3 + (1/tau+c)s^2 + (c/tau+k+kb)s + k/tau
    // = (s+1)(s+2)(s+3). Its exact multipliers over T=1 are exp(-1,-2,-3).
    let omega = core::f64::consts::TAU;
    let system = fs_phs::mass_spring_damper(1.0, 1.5, 2.0)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: 1.5,
            relaxation_time_s: 0.25,
        }])
        .unwrap();
    let input = PeriodicPortInput::new(omega, vec![vec![C64::ZERO, C64::new(0.1, 0.03)]]).unwrap();
    let candidate = solve_forced_harmonic_balance(&system, &input, 0.01, &hb(2)).unwrap();
    let result =
        analyze_forced_orbit_stability(&system, &input, &candidate, &controls(vec![1e-7; 3]))
            .unwrap();
    assert_eq!(result.classification, PhsNumericalStability::Stable);
    for level in &result.levels {
        assert_eq!(level.orbit.x0.len(), 3);
        assert_eq!(level.orbit.multipliers.len(), 3);
        let mut values = level.orbit.multipliers.clone();
        values.sort_by(|a, b| b.0.total_cmp(&a.0));
        for (i, &(re, im)) in values.iter().enumerate() {
            assert!((re - det::exp(-((i + 1) as f64))).abs() < 2e-5);
            assert!(im.abs() < 2e-5);
        }
    }
    assert!((result.dominant_growth_rate_per_s.unwrap() + 1.0).abs() < 2e-5);
}

struct CubicSpring;
impl Storage for CubicSpring {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * (x[0] * x[0] + x[1] * x[1]) + 0.25 * 0.15 * x[0].powi(4)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = x[0] + 0.15 * x[0].powi(3);
        out[1] = x[1];
    }
}

#[test]
fn continued_duffing_branch_distinguishes_stable_outer_and_unstable_middle_orbits() {
    let system = PortHamiltonian::new(
        2,
        1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, 0.1],
        vec![0.0, 1.0],
        Box::new(CubicSpring),
    )
    .unwrap();
    let input = PeriodicPortInput::new(0.6, vec![vec![C64::ZERO, C64::from_re(0.15)]]).unwrap();
    let budget = ContinuationBudget {
        max_steps: 400,
        initial_step: 0.02,
        min_step: 1e-5,
        hb: hb(7),
    };
    let branch =
        continue_forced_frequency_response(&system, &input, 1.8, 0.35, &budget, 128).unwrap();
    let folds: Vec<usize> = branch
        .windows(3)
        .enumerate()
        .filter_map(|(i, p)| {
            let before = p[1].point.lambda - p[0].point.lambda;
            let after = p[2].point.lambda - p[1].point.lambda;
            (before * after < 0.0).then_some(i + 1)
        })
        .collect();
    assert_eq!(folds.len(), 2);
    let segments = [
        &branch[..folds[0]],
        &branch[folds[0] + 1..folds[1]],
        &branch[folds[1] + 1..],
    ];
    for (index, segment) in segments.iter().enumerate() {
        let point = segment
            .iter()
            .min_by(|a, b| {
                (a.point.lambda - 1.22)
                    .abs()
                    .total_cmp(&(b.point.lambda - 1.22).abs())
            })
            .unwrap();
        let forcing = input.at_frequency(point.point.lambda).unwrap();
        let result = analyze_forced_orbit_stability(
            &system,
            &forcing,
            &point.point.orbit,
            &controls(vec![1e-5; 2]),
        )
        .unwrap();
        let expected = if index == 1 {
            PhsNumericalStability::Unstable
        } else {
            PhsNumericalStability::Stable
        };
        assert_eq!(
            result.classification, expected,
            "branch segment {index}: {result:?}"
        );
        assert_eq!(result.dominant_growth_rate_per_s.unwrap() > 0.0, index == 1);
        // Liouville's identity is independent of HB and monodromy differencing:
        // divergence = -0.1, hence det(M) = exp(-0.1*T) on every cubic orbit.
        for level in &result.levels {
            assert_eq!(level.orbit.multipliers.len(), 2);
            let a = level.orbit.multipliers[0];
            let b = level.orbit.multipliers[1];
            let determinant = C64::new(a.0, a.1) * C64::new(b.0, b.1);
            assert!((determinant.re - det::exp(-0.1 * level.orbit.period)).abs() < 2e-5);
            assert!(determinant.im.abs() < 2e-5);
        }
    }
    assert_eq!(input.omega_rad_s(), 0.6);
}

#[test]
fn guarded_unit_circle_and_unsettled_time_grids_remain_unresolved() {
    let system = fs_phs::mass_spring_damper(1.0, 4.0, 0.0002).unwrap();
    let input = PeriodicPortInput::new(2.0, vec![vec![C64::ZERO, C64::from_re(0.5e-5)]]).unwrap();
    let candidate = solve_forced_harmonic_balance(&system, &input, 0.025, &hb(2)).unwrap();
    let result =
        analyze_forced_orbit_stability(&system, &input, &candidate, &controls(vec![1e-4; 2]))
            .unwrap();
    assert_eq!(result.classification, PhsNumericalStability::Unresolved);
    assert!(
        result
            .levels
            .iter()
            .all(|level| (level.spectral_radius - 1.0).abs() < 1e-3)
    );

    let system = fs_phs::mass_spring_damper(1.0, 4.0, 0.4).unwrap();
    let input = PeriodicPortInput::new(1.7, vec![vec![C64::ZERO, C64::from_re(0.1)]]).unwrap();
    let candidate = solve_forced_harmonic_balance(&system, &input, 0.1, &hb(1)).unwrap();
    let mut coarse = controls(vec![1.0; 2]);
    coarse.shooting.steps_per_period = 8;
    coarse.spectral_radius_tolerance = 1e-6;
    let result = analyze_forced_orbit_stability(&system, &input, &candidate, &coarse).unwrap();
    assert_eq!(result.classification, PhsNumericalStability::Unresolved);
    assert!((result.levels[1].spectral_radius - result.levels[0].spectral_radius).abs() > 1e-6);
}

fn scalar_candidate() -> HbOrbit {
    HbOrbit {
        omega: 1.0,
        coeffs: vec![vec![C64::ZERO, C64::from_re(0.1)]],
        residual_trace: vec![0.0],
        residual: 0.0,
        aft_samples: 8,
        sampling_defect: 0.0,
    }
}

#[test]
fn malformed_stability_requests_refuse_before_original_source_evaluation() {
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
    struct Counted(Arc<AtomicUsize>);
    impl Storage for Counted {
        fn hamiltonian(&self, x: &[f64]) -> f64 {
            self.0.fetch_add(1, Ordering::Relaxed);
            0.5 * x[0] * x[0]
        }
        fn gradient(&self, x: &[f64], out: &mut [f64]) {
            self.0.fetch_add(1, Ordering::Relaxed);
            out[0] = x[0];
        }
    }
    let calls = Arc::new(AtomicUsize::new(0));
    let system = PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![1.0],
        vec![1.0],
        Box::new(Counted(calls.clone())),
    )
    .unwrap();
    calls.store(0, Ordering::Relaxed);
    let input = PeriodicPortInput::new(1.0, vec![vec![C64::ZERO, C64::from_re(0.5)]]).unwrap();
    let candidate = scalar_candidate();
    let valid = controls(vec![1e-6]);
    let wrong_clock = input.at_frequency(2.0).unwrap();
    assert!(analyze_forced_orbit_stability(&system, &wrong_clock, &candidate, &valid).is_err());
    let wrong_ports = PeriodicPortInput::new(1.0, vec![]).unwrap();
    assert!(analyze_forced_orbit_stability(&system, &wrong_ports, &candidate, &valid).is_err());
    for mode in 0..5 {
        let mut invalid = valid.clone();
        match mode {
            0 => invalid.shooting.steps_per_period = usize::MAX,
            1 => invalid.shooting.steps_per_period = 2,
            2 => invalid.shooting.max_newton = 0,
            3 => invalid.state_tolerances = vec![f64::NAN],
            _ => invalid.unit_circle_guard = 1.0,
        }
        assert!(analyze_forced_orbit_stability(&system, &input, &candidate, &invalid).is_err());
    }
    let mut missing = candidate.clone();
    missing.coeffs[0] = vec![];
    assert!(analyze_forced_orbit_stability(&system, &input, &missing, &valid).is_err());
    let mut coefficients = vec![C64::ZERO; 10];
    coefficients[9] = C64::from_re(0.5);
    let unresolved = PeriodicPortInput::new(1.0, vec![coefficients]).unwrap();
    assert!(analyze_forced_orbit_stability(&system, &unresolved, &candidate, &valid).is_err());
    assert_eq!(calls.load(Ordering::Relaxed), 0);
    assert!(input.at_frequency(f64::NAN).is_err());
}

#[test]
fn source_refusal_and_large_shooting_corrections_do_not_receive_stability_labels() {
    struct Scalar {
        refuses: bool,
    }
    impl Storage for Scalar {
        fn hamiltonian(&self, x: &[f64]) -> f64 {
            0.5 * x[0] * x[0]
        }
        fn gradient(&self, x: &[f64], out: &mut [f64]) {
            out[0] = if self.refuses { f64::NAN } else { x[0] };
        }
    }
    let input = PeriodicPortInput::new(1.0, vec![vec![C64::ZERO, C64::from_re(0.5)]]).unwrap();
    let candidate = scalar_candidate();
    for refuses in [true, false] {
        let system = PortHamiltonian::new(
            1,
            1,
            vec![0.0],
            vec![1.0],
            vec![1.0],
            Box::new(Scalar { refuses }),
        )
        .unwrap();
        let error =
            analyze_forced_orbit_stability(&system, &input, &candidate, &controls(vec![1e-6]))
                .unwrap_err();
        if refuses {
            assert_eq!(error, PhsOrbitError::SourceRefused);
        } else {
            assert!(
                matches!(error, PhsOrbitError::SeedDisagreement { coordinate: 0, absolute_change, .. } if absolute_change > 0.2)
            );
        }
    }
}
