use super::*;
use fs_math::{c64::C64, det};
use fs_orbit::HbBudget;
use fs_phs::{RelaxationBranch, Storage};

fn budget(harmonics: usize) -> ContinuationBudget {
    ContinuationBudget {
        max_steps: 400,
        initial_step: 0.02,
        min_step: 1e-5,
        hb: HbBudget {
            harmonics,
            max_newton: 30,
            tolerance: 1e-10,
            deflate_equilibrium: false,
        },
    }
}

struct Spring {
    mass: f64,
    stiffness: f64,
    cubic: f64,
}

impl Storage for Spring {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * self.stiffness * x[0] * x[0]
            + 0.25 * self.cubic * x[0].powi(4)
            + 0.5 * x[1] * x[1] / self.mass
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = self.stiffness * x[0] + self.cubic * x[0].powi(3);
        out[1] = x[1] / self.mass;
    }
}

#[test]
fn linear_branch_preserves_complex_multiport_phase_dc_and_physical_work() {
    let (mass, stiffness, damping) = (1.0, 4.0, 0.4);
    let system = PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, damping],
        vec![0.0, 0.0, 1.0, 0.5],
        Box::new(Spring {
            mass,
            stiffness,
            cubic: 0.0,
        }),
    )
    .unwrap();
    let coefficients = vec![
        vec![
            C64::from_re(0.4),
            C64::new(0.2, -0.08),
            C64::new(0.01, 0.015),
        ],
        vec![
            C64::from_re(-0.1),
            C64::new(0.1, 0.04),
            C64::new(-0.008, 0.001),
        ],
    ];
    let input = PeriodicPortInput::new(1.1, coefficients.clone()).unwrap();
    let branch =
        continue_forced_frequency_response(&system, &input, 2.8, 0.2, &budget(2), 256).unwrap();
    assert!(branch.len() > 10);
    assert_eq!(branch[0].point.lambda, input.omega_rad_s());
    assert!(branch.last().unwrap().point.lambda >= 2.8);
    assert!(branch.iter().any(|point| point.point.lambda < 2.0));
    assert!(branch.iter().any(|point| point.point.lambda > 2.0));
    for response in &branch {
        let orbit = &response.point.orbit;
        assert_eq!(response.point.lambda.to_bits(), orbit.omega.to_bits());
        let mut loss = 0.0;
        for n in 0..=2 {
            let omega = n as f64 * orbit.omega;
            let force = coefficients[0][n] + coefficients[1][n].scale(0.5);
            let q = force / C64::new(stiffness - mass * omega * omega, damping * omega);
            let p = q * C64::new(0.0, mass * omega);
            assert!((orbit.coeffs[0][n] - q).abs() < 2e-8 * q.abs().max(1e-4));
            assert!(
                (orbit.coeffs[1][n] - p).abs() < 2e-10 + 2e-8 * p.abs(),
                "omega={omega}, harmonic={n}, measured={:?}, expected={p:?}, residual={}",
                orbit.coeffs[1][n],
                orbit.residual
            );
            if n > 0 {
                loss += 2.0 * damping * omega * omega * q.abs().powi(2);
            }
        }
        assert!((response.audit.mean_resistance_loss_w - loss).abs() < 1e-8 * loss);
        assert!(response.audit.relative_power_residual < 1e-8);
        assert!(orbit.residual.is_finite() && orbit.sampling_defect.is_finite());
    }
    // Continuation must not mutate the caller's frequency or forcing chart.
    assert_eq!(input.omega_rad_s(), 1.1);
    assert_eq!(input.coefficients, coefficients);
}

#[test]
fn maxwell_branch_keeps_the_memory_coordinate_and_frequency_dependent_loss() {
    let (stiffness, branch_stiffness, tau, damping) = (4.0, 3.0, 0.4, 0.02);
    let system = fs_phs::mass_spring_damper(1.0, stiffness, damping)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: branch_stiffness,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let force = C64::new(0.1, -0.04);
    let input = PeriodicPortInput::new(1.2, vec![vec![C64::ZERO, force]]).unwrap();
    let branch =
        continue_forced_frequency_response(&system, &input, 3.2, 0.1, &budget(2), 256).unwrap();
    for response in &branch {
        let orbit = &response.point.orbit;
        let omega = orbit.omega;
        let iw = C64::new(0.0, omega);
        let modulus = (iw.scale(tau) / (C64::ONE + iw.scale(tau))).scale(branch_stiffness);
        let q = force / (C64::new(stiffness - omega * omega, damping * omega) + modulus);
        assert_eq!(orbit.coeffs.len(), 3);
        assert!((orbit.coeffs[0][1] - q).abs() < 2e-8 * q.abs());
        assert!(orbit.first_harmonic_amplitude(2) > 0.0);
        let viscous_only = 2.0 * damping * omega * omega * q.abs().powi(2);
        assert!(response.audit.mean_resistance_loss_w > viscous_only);
        assert!(response.audit.relative_power_residual < 1e-8);
    }
}

// Independent first-harmonic scalar Duffing law. Its fold prediction is an
// approximation to the full cubic-storage branch, not a replacement for it.
fn first_harmonic_fold_frequencies(cubic: f64, damping: f64, force: f64) -> Vec<f64> {
    let mut folds = Vec::new();
    let mut previous = f64::NAN;
    let mut slope = 0.0;
    for i in 200..60_000 {
        let amplitude = i as f64 * 1e-4;
        let a = 1.0 + 0.75 * cubic * amplitude * amplitude;
        let b = a - 0.5 * damping * damping;
        let discriminant = b * b - a * a + (force / amplitude).powi(2);
        if discriminant >= 0.0 {
            let omega = det::sqrt(b + det::sqrt(discriminant));
            if previous.is_finite() {
                let next_slope = omega - previous;
                if next_slope * slope < 0.0 {
                    folds.push(previous);
                }
                slope = next_slope;
            }
            previous = omega;
        }
    }
    folds
}

#[test]
fn nonlinear_branch_traverses_both_folds_with_original_cubic_storage() {
    let (cubic, damping, force) = (0.15, 0.1, 0.3);
    let system = PortHamiltonian::new(
        2,
        1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, damping],
        vec![0.0, 1.0],
        Box::new(Spring {
            mass: 1.0,
            stiffness: 1.0,
            cubic,
        }),
    )
    .unwrap();
    let input =
        PeriodicPortInput::new(0.6, vec![vec![C64::ZERO, C64::from_re(0.5 * force)]]).unwrap();
    let branch =
        continue_forced_frequency_response(&system, &input, 1.8, 0.35, &budget(7), 256).unwrap();
    let mut folds = Vec::new();
    for points in branch.windows(3) {
        let before = points[1].point.lambda - points[0].point.lambda;
        let after = points[2].point.lambda - points[1].point.lambda;
        if before * after < 0.0 {
            folds.push(points[1].point.lambda);
        }
    }
    assert_eq!(folds.len(), 2, "frequency reversals: {folds:?}");
    let reference = first_harmonic_fold_frequencies(cubic, damping, force);
    assert_eq!(reference.len(), 2, "scalar fold scan: {reference:?}");
    // The scalar scan encounters the lower fold first; the continued branch
    // from low frequency encounters the upper fold first.
    assert!((folds[0] / reference[1] - 1.0).abs() < 0.03);
    assert!((folds[1] / reference[0] - 1.0).abs() < 0.03);
    let mut generated_third: f64 = 0.0;
    for response in &branch {
        generated_third = generated_third.max(response.point.orbit.coeffs[0][3].abs());
        assert!(
            response.audit.relative_power_residual < 2e-8,
            "{:?}",
            response.audit
        );
        assert!(response.audit.mean_resistance_loss_w > 0.0);
    }
    assert!(generated_third > 1e-3, "higher harmonics were lost");
    assert!(branch.last().unwrap().point.lambda >= 1.8);
}

#[test]
fn malformed_branch_and_unresolved_forcing_refuse_before_source_evaluation() {
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
    for endpoint in [0.0, -1.0, f64::NAN, f64::INFINITY, 1.0] {
        assert!(
            continue_forced_frequency_response(&system, &input, endpoint, 0.1, &budget(1), 64)
                .is_err()
        );
    }
    for samples in [0, 31, 65_537] {
        assert!(
            continue_forced_frequency_response(&system, &input, 2.0, 0.1, &budget(1), samples)
                .is_err()
        );
    }
    let mut bad_budget = budget(1);
    bad_budget.max_steps = 0;
    assert!(
        continue_forced_frequency_response(&system, &input, 2.0, 0.1, &bad_budget, 64).is_err()
    );
    let mut coefficients = vec![C64::ZERO; 10];
    coefficients[9] = C64::new(0.0, f64::from_bits(1));
    let unresolved = PeriodicPortInput::new(1.0, vec![coefficients]).unwrap();
    assert_eq!(
        continue_forced_frequency_response(&system, &unresolved, 2.0, 0.1, &budget(1), 64)
            .unwrap_err(),
        PhsOrbitError::Invalid {
            what: "nonzero forcing harmonic exceeds the retained HB truncation"
        },
    );
    assert_eq!(calls.load(Ordering::Relaxed), 0);
}

#[test]
fn original_source_refusal_is_not_a_partial_resonance_branch() {
    struct Refusing;
    impl Storage for Refusing {
        fn hamiltonian(&self, _: &[f64]) -> f64 {
            f64::NAN
        }
        fn gradient(&self, _: &[f64], out: &mut [f64]) {
            out.fill(f64::NAN);
        }
    }
    let system =
        PortHamiltonian::new(1, 1, vec![0.0], vec![1.0], vec![1.0], Box::new(Refusing)).unwrap();
    let input = PeriodicPortInput::new(1.0, vec![vec![C64::ZERO, C64::from_re(0.5)]]).unwrap();
    assert_eq!(
        continue_forced_frequency_response(&system, &input, 2.0, 0.1, &budget(1), 64).unwrap_err(),
        PhsOrbitError::SourceRefused,
    );
}
