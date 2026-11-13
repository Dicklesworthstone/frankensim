use super::*;
use fs_math::{c64::C64, det};
use fs_orbit::autonomous::AutonomousScales;
use fs_phs::{RelaxationBranch, Storage};

struct Energy;
impl Storage for Energy {
    fn hamiltonian(&self, x: &[f64]) -> f64 { 0.5 * (x[0]*x[0] + x[1]*x[1]) }
    fn gradient(&self, x: &[f64], out: &mut [f64]) { out.copy_from_slice(x); }
}
fn system() -> PortHamiltonian {
    system_with_storage(Box::new(Energy))
}
fn system_with_storage(storage: Box<dyn Storage>) -> PortHamiltonian {
    PortHamiltonian::new(2, 2, vec![0.0, 1.0, -1.0, 0.0],
        vec![0.2, 0.0, 0.0, 0.2], vec![1.0, 0.0, 0.0, 1.0], storage)
        .unwrap().with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0], stiffness: 0.3, relaxation_time_s: 2.0,
        }]).unwrap()
}

// Explicit actuator law supplies both the native losses and the changing
// angular speed. This is an analytical test excitation, not a material law.
// It never edits original R, H, G or the actual Maxwell state equation.
struct Drive { calls: Cell<usize>, refusal: u8 }
impl AutonomousFeedbackFamily for Drive {
    fn port_count(&self) -> usize { 2 }
    fn input(&self, parameter: f64, x: &[f64], flows: &[f64], out: &mut [f64]) -> bool {
        self.calls.set(self.calls.get() + 1);
        if self.refusal == 1 { return false; }
        if self.refusal == 2 { out[0] = 0.0; return true; }
        let a = x[0]*x[0] + x[1]*x[1];
        let radial = parameter - (a - 1.0).powi(2);
        let omega = 1.0 + 0.2*a;
        out[0] = (omega - 1.0)*x[1] + radial*x[0] + 0.2*flows[0];
        out[1] = flows[0] - omega*x[0] + (radial + 0.2)*x[1];
        true
    }
}
fn seed(radius: f64) -> HbOrbit {
    let omega = 1.0 + 0.2*radius*radius;
    let mut coeffs = vec![vec![C64::ZERO; 4]; 3];
    coeffs[0][1] = C64::from_re(0.5*radius);
    coeffs[1][1] = C64::new(0.0, 0.5*radius);
    coeffs[2][1] = coeffs[0][1].scale(det::sqrt(0.3)) / C64::new(1.0, 2.0*omega);
    HbOrbit { omega, coeffs, residual_trace: vec![0.0], residual: 0.0,
        aft_samples: 0, sampling_defect: 0.0 }
}
fn lambda(radius: f64) -> f64 { (radius*radius - 1.0).powi(2) }
fn controls() -> PhsAutonomousBranchControls {
    PhsAutonomousBranchControls {
        continuation: AutonomousContinuationBudget {
            harmonics: 3, steps: 24, max_newton: 40, max_attempts_per_step: 12,
            step_sizes: [0.045, 1e-5, 0.08], tolerance: 2e-10,
            minimum_oscillation: 1e-5, maximum_unknowns: 64,
            maximum_payload_bytes: 16 << 20, maximum_residual_evaluations: 20_000,
            scales: AutonomousScales { states: vec![1.0;3], equations: vec![1.0;3],
                angular_frequency: 2.0, parameter: 1.0 },
        },
        audit_samples: 64, maximum_audit_evaluations: 10_000,
        maximum_relative_power_residual: 2e-8, maximum_relative_audit_change: 2e-8,
    }
}

#[test]
fn source_driven_branch_keeps_memory_phase_and_original_energy_and_heat_through_the_fold() {
    let system = system();
    let source = Drive { calls: Cell::new(0), refusal: 0 };
    let seeds = [seed(0.65), seed(0.70)];
    let before = system.hamiltonian(&[0.3, -0.2, 0.1]);
    let result = continue_autonomous_feedback_branch(&system, &source,
        [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &controls(), || false).unwrap();
    assert_eq!(result.continuation.points.len(), 26);
    assert_eq!(result.audits.len(), 26);
    assert_eq!(result.audit_evaluations, 26*(64+129));
    assert!(result.continuation.points.windows(2).any(|p| p[1].parameter < p[0].parameter));
    assert!(result.continuation.points.windows(2).any(|p| p[1].parameter > p[0].parameter));
    for (point, audits) in result.continuation.points.iter().zip(&result.audits) {
        let orbit = &point.orbit;
        assert_eq!(orbit.coeffs.len(), 3);
        let radius = orbit.first_harmonic_amplitude(0);
        assert!((point.parameter-lambda(radius)).abs() < 2e-8);
        assert!((orbit.omega-(1.0+0.2*radius*radius)).abs() < 2e-8);
        let q = orbit.coeffs[0][1];
        let z = q.scale(det::sqrt(0.3))/C64::new(1.0, 2.0*orbit.omega);
        assert!((z-orbit.coeffs[2][1]).abs() < 2e-9);
        assert!(z.im.abs() > 1e-3, "original memory lag must survive");
        // Independent complex-coefficient moments for the ORIGINAL H and R.
        let eq = q.scale(1.3) - z.scale(det::sqrt(0.3));
        let ez = z - q.scale(det::sqrt(0.3));
        let p = orbit.coeffs[1][1];
        let expected_heat = 0.4*(eq.norm_sq()+p.norm_sq()) + ez.norm_sq();
        let expected_energy = q.norm_sq()+p.norm_sq()+ez.norm_sq();
        for audit in audits {
            assert!(audit.relative_power_residual < 2e-8);
            assert!((audit.mean_resistance_loss_w-expected_heat).abs() < 2e-8*expected_heat);
            assert!((audit.mean_supplied_power_w-expected_heat).abs() < 2e-8*expected_heat);
            assert!((audit.mean_stored_energy_j-expected_energy).abs() < 2e-8*expected_energy);
        }
    }
    assert_eq!(system.hamiltonian(&[0.3, -0.2, 0.1]), before);
    assert_eq!(system.structure().1, &[0.2,0.0,0.0,0.0,0.2,0.0,0.0,0.0,0.5]);
}

#[test]
fn fixed_source_view_drives_the_existing_independent_stability_and_audit_paths() {
    use crate::phs_orbit::{PhsAutonomousControls, PhsFloquetControls, PhsNumericalStability,
        analyze_autonomous_feedback_stability, audit_autonomous_feedback_orbit};
    let system = system();
    let source = Drive { calls: Cell::new(0), refusal: 0 };
    for radius in [0.8, 1.2] {
        let candidate = seed(radius);
        let law = FixedFeedback::new(&source, lambda(radius)).unwrap();
        let c = PhsAutonomousControls {
            floquet: PhsFloquetControls {
                shooting: fs_orbit::ShootBudget { steps_per_period: 512, max_newton: 40, tolerance: 1e-10 },
                state_tolerances: vec![1e-5;3], spectral_radius_tolerance: 2e-3,
                unit_circle_guard: 1e-2,
            },
            relative_period_tolerance: 1e-5, phase_multiplier_tolerance: 1e-3, audit_samples: 128,
        };
        let answer = analyze_autonomous_feedback_stability(&system, &law, &candidate, &c).unwrap();
        assert_eq!(answer.classification, if radius < 1.0 {
            PhsNumericalStability::Unstable
        } else { PhsNumericalStability::Stable });
        // Radial exponent = -4*a*(a-1); the independent memory direction has -1/2.
        let a = radius*radius;
        let expected_growth = (-4.0*a*(a-1.0)).max(-0.5);
        assert!((answer.dominant_transverse_growth_rate_per_s.unwrap()-expected_growth).abs() < 2e-3);
        let audit = audit_autonomous_feedback_orbit(&system, &law, &candidate, 128).unwrap();
        assert!(audit.relative_power_residual < 1e-12);
        assert_eq!(law.parameter(), lambda(radius));
    }
}

#[test]
fn physical_refusals_and_incomplete_source_efforts_are_not_repaired() {
    let system = system();
    let seeds = [seed(0.65), seed(0.70)];
    let source = Drive { calls: Cell::new(0), refusal: 0 };
    for kind in 0..4 {
        let mut c = controls();
        match kind {
            0 => c.maximum_audit_evaluations = 1,
            1 => c.audit_samples = 0,
            2 => c.maximum_relative_power_residual = f64::NAN,
            _ => c.continuation.maximum_payload_bytes = 1,
        }
        assert!(continue_autonomous_feedback_branch(&system, &source,
            [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &c, || false).is_err());
        assert_eq!(source.calls.get(), 0);
    }
    for refusal in [1, 2] {
        let broken = Drive { calls: Cell::new(0), refusal };
        assert!(matches!(continue_autonomous_feedback_branch(&system, &broken,
            [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &controls(), || false),
            Err(PhsAutonomousBranchError::Source(PhsOrbitError::SourceRefused))));
    }
    assert!(FixedFeedback::new(&source, f64::NAN).is_err());
}

#[test]
fn cancellation_reaches_physical_quadrature_and_preserves_the_original_source() {
    let system = system();
    let source = Drive { calls: Cell::new(0), refusal: 0 };
    let before = system.hamiltonian(&[0.3, -0.2, 0.1]);
    let candidate = seed(0.8);
    let mut polls = 0;
    assert!(matches!(audited(&system, &source, lambda(0.8), &candidate, 128, &mut || {
        polls += 1; polls >= 13
    }), Err(PhsAutonomousBranchError::Cancelled)));
    assert!(source.calls.get() < 128);
    assert_eq!(system.hamiltonian(&[0.3, -0.2, 0.1]), before);
    let audit = audited(&system, &source, lambda(0.8), &candidate, 128, &mut || false).unwrap();
    assert!(audit.relative_power_residual < 1e-12);
    let other = seed(0.9);
    assert!(matches!(continue_autonomous_feedback_branch(&system, &source,
        [(lambda(0.8), &candidate), (lambda(0.9), &other)], &controls(), || true),
        Err(PhsAutonomousBranchError::Cancelled)));
}

fn nonlinear_system(beta: f64) -> PortHamiltonian {
    let storage = fs_nlmodal::SosModalStorage {
        omegas: vec![1.0],
        channels: vec![fs_nlmodal::StressChannel {
            coefficient: beta,
            coupling: vec![1.0],
        }],
    };
    storage.validate().unwrap();
    system_with_storage(Box::new(storage))
}

fn total_feedback_derivative(beta: f64, parameter: f64, x: &[f64], out: &mut [f64]) -> bool {
    let a = x[0] * x[0] + x[1] * x[1];
    let radial = parameter - (a - 1.0).powi(2);
    let omega = 1.0 + 0.2 * a;
    let drdq = -4.0 * (a - 1.0) * x[0];
    let drdp = -4.0 * (a - 1.0) * x[1];
    // The actual flow is q + beta*q^3 + .3*q - sqrt(.3)*z.
    // Its complete Hessian chain belongs in the TOTAL feedback derivative.
    let flow_q = 1.3 + 3.0 * beta * x[0] * x[0];
    let flow_z = -det::sqrt(0.3);
    out.copy_from_slice(&[
        0.4 * x[0] * x[1] + radial + x[0] * drdq + 0.2 * flow_q,
        omega - 1.0 + 0.4 * x[1] * x[1] + x[0] * drdp,
        0.2 * flow_z,
        flow_q - omega - 0.4 * x[0] * x[0] + x[1] * drdq,
        -0.4 * x[0] * x[1] + radial + 0.2 + x[1] * drdp,
        flow_z,
    ]);
    true
}

#[test]
fn analytic_family_retains_native_nonlinear_storage_memory_and_original_cycle_work() {
    let beta = 0.6;
    let system = nonlinear_system(beta);
    assert!(!system.storage().quadratic_part().unwrap().exact);
    assert!(system.storage().continuous_tangent().is_some());
    let source = Drive { calls: Cell::new(0), refusal: 0 };
    let probed_source = Drive { calls: Cell::new(0), refusal: 0 };
    let derivative_calls = Cell::new(0);
    let derivative = |parameter: f64, x: &[f64], _: &[f64], out: &mut [f64]| {
        derivative_calls.set(derivative_calls.get() + 1);
        total_feedback_derivative(beta, parameter, x, out)
    };
    let seeds = [seed(0.65), seed(0.70)];
    let before = system.hamiltonian(&[0.3, -0.2, 0.1]);
    let exact = continue_autonomous_feedback_branch_with_tangent(&system, &source, &derivative,
        [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &controls(), || false).unwrap();
    let probed = continue_autonomous_feedback_branch(&system, &probed_source,
        [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &controls(), || false).unwrap();
    assert_eq!(exact.continuation.points.len(), 26);
    assert!(exact.continuation.tangent_preparations > 0);
    assert_eq!(probed.continuation.tangent_preparations, 0);
    assert!(derivative_calls.get() > 0);
    assert!(3 * exact.continuation.residual_evaluations < probed.continuation.residual_evaluations);
    // These counts include the unchanged physical audits in both routes.
    assert!(2 * source.calls.get() < probed_source.calls.get());
    for ((point, other), audits) in exact.continuation.points.iter()
        .zip(&probed.continuation.points).zip(&exact.audits)
    {
        assert!((point.parameter - other.parameter).abs() < 2e-7);
        assert!((point.orbit.omega - other.orbit.omega).abs() < 2e-7);
        for (row, reference) in point.orbit.coeffs.iter().zip(&other.orbit.coeffs) {
            for (&actual, &expected) in row.iter().zip(reference) {
                assert!((actual - expected).abs() < 2e-7);
            }
        }
        let radius = point.orbit.first_harmonic_amplitude(0);
        assert!((point.parameter - lambda(radius)).abs() < 2e-8);
        assert!((point.orbit.omega - (1.0 + 0.2 * radius * radius)).abs() < 2e-8);
        let q = point.orbit.coeffs[0][1];
        let p = point.orbit.coeffs[1][1];
        let z = q.scale(det::sqrt(0.3)) / C64::new(1.0, 2.0 * point.orbit.omega);
        assert!((z - point.orbit.coeffs[2][1]).abs() < 2e-9);
        let eq1 = q.scale(1.3 + 3.0 * beta * q.norm_sq()) - z.scale(det::sqrt(0.3));
        let eq3 = (q * q * q).scale(beta);
        let ez = z - q.scale(det::sqrt(0.3));
        let heat = 0.4 * (eq1.norm_sq() + eq3.norm_sq() + p.norm_sq()) + ez.norm_sq();
        let energy = q.norm_sq() + p.norm_sq() + ez.norm_sq()
            + 1.5 * beta * q.norm_sq() * q.norm_sq();
        for audit in audits {
            assert!((audit.mean_resistance_loss_w - heat).abs() < 3e-8 * heat);
            assert!((audit.mean_supplied_power_w - heat).abs() < 3e-8 * heat);
            assert!((audit.mean_stored_energy_j - energy).abs() < 3e-8 * energy);
        }
    }
    assert!(exact.continuation.points.windows(2).any(|p| p[1].parameter > p[0].parameter));
    assert!(exact.continuation.points.windows(2).any(|p| p[1].parameter < p[0].parameter));
    assert_eq!(system.hamiltonian(&[0.3, -0.2, 0.1]), before);
}

#[test]
fn analytic_family_refuses_missing_or_incomplete_derivatives_without_substitution() {
    let seeds = [seed(0.65), seed(0.70)];
    let source = Drive { calls: Cell::new(0), refusal: 0 };
    let derivative_calls = Cell::new(0);
    let complete = |parameter: f64, x: &[f64], _: &[f64], out: &mut [f64]| {
        derivative_calls.set(derivative_calls.get() + 1);
        total_feedback_derivative(0.6, parameter, x, out)
    };
    let unavailable = system();
    assert!(matches!(continue_autonomous_feedback_branch_with_tangent(
        &unavailable, &source, &complete,
        [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &controls(), || false),
        Err(PhsAutonomousBranchError::Invalid(_))));
    assert_eq!(source.calls.get(), 0);
    assert_eq!(derivative_calls.get(), 0);
    let system = nonlinear_system(0.6);
    for partial in [false, true] {
        let calls = Cell::new(0);
        let broken = |_: f64, _: &[f64], _: &[f64], out: &mut [f64]| {
            calls.set(calls.get() + 1);
            if partial { out[0] = 0.0; }
            partial
        };
        assert!(matches!(continue_autonomous_feedback_branch_with_tangent(
            &system, &source, &broken,
            [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &controls(), || false),
            Err(PhsAutonomousBranchError::Continuation(
                AutonomousContinuationError::Numerical(fs_orbit::OrbitError::TangentRefused { .. })))));
        assert_eq!(calls.get(), 1);
    }
    let broken_source = Drive { calls: Cell::new(0), refusal: 1 };
    assert!(matches!(continue_autonomous_feedback_branch_with_tangent(
        &system, &broken_source, &complete,
        [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &controls(), || false),
        Err(PhsAutonomousBranchError::Source(PhsOrbitError::SourceRefused))));
}

#[test]
fn analytic_family_honors_payload_limits_and_cancellation_before_publication() {
    let system = nonlinear_system(0.6);
    let seeds = [seed(0.65), seed(0.70)];
    let source = Drive { calls: Cell::new(0), refusal: 0 };
    let derivative_calls = Cell::new(0);
    let derivative = |parameter: f64, x: &[f64], _: &[f64], out: &mut [f64]| {
        derivative_calls.set(derivative_calls.get() + 1);
        total_feedback_derivative(0.6, parameter, x, out)
    };
    for kind in 0..3 {
        let mut c = controls();
        match kind {
            0 => c.continuation.maximum_unknowns = 1,
            1 => c.continuation.maximum_payload_bytes = 1,
            _ => c.maximum_audit_evaluations = 1,
        }
        assert!(continue_autonomous_feedback_branch_with_tangent(&system, &source, &derivative,
            [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &c, || false).is_err());
        assert_eq!(source.calls.get(), 0);
        assert_eq!(derivative_calls.get(), 0);
    }
    let before = system.hamiltonian(&[0.3, -0.2, 0.1]);
    assert!(matches!(continue_autonomous_feedback_branch_with_tangent(&system, &source, &derivative,
        [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &controls(),
        || derivative_calls.get() > 0), Err(PhsAutonomousBranchError::Cancelled)));
    assert!(derivative_calls.get() > 0);
    assert_eq!(system.hamiltonian(&[0.3, -0.2, 0.1]), before);
}
