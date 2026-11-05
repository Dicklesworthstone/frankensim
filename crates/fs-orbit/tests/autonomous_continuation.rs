use fs_math::{c64::C64, det};
use fs_orbit::{ContinuableProblem, HbOrbit, OrbitProblem};
use fs_orbit::autonomous::{AutonomousContinuationBudget, AutonomousScales,
    AutonomousContinuationError, continue_autonomous_branch};
use std::cell::Cell;

#[derive(Clone, Copy)]
enum Tangent { Absent, Exact, Refuse, Incomplete }

// An analytically soluble saddle-node of cycles. The first state is static,
// deliberately invalidating a phase gauge based on its first harmonic.
// a=x1^2+x2^2, lambda=(a-1)^2, omega=1+.2*a.
struct Fold {
    parameter: f64,
    units: [f64; 3],
    equation_units: [f64; 3],
    evaluations: Cell<usize>,
    forced: bool,
    tangent: Tangent,
    tangent_evaluations: Cell<usize>,
    cubic_observation: bool,
}
impl Fold {
    fn new() -> Self {
        Self { parameter: 0.0, units: [1.0; 3], equation_units: [1.0; 3],
            evaluations: Cell::new(0), forced: false, tangent: Tangent::Absent,
            tangent_evaluations: Cell::new(0), cubic_observation: false }
    }
    fn seed(&self, radius: f64, phase: f64, n: usize) -> HbOrbit {
        let mut coeffs = vec![vec![C64::ZERO; n + 1]; 3];
        coeffs[0][0] = C64::from_re(0.3 * self.units[0]);
        let rotation = C64::new(det::cos(phase), det::sin(phase));
        coeffs[1][1] = rotation.scale(0.5 * radius * self.units[1]);
        coeffs[2][1] = (rotation * C64::new(0.0, 1.0)).scale(0.5 * radius * self.units[2]);
        let omega = 1.0 + 0.2 * radius * radius;
        if self.cubic_observation {
            // x0' + x0 = 0.3 + x1^3: a physical cubic generates both
            // fundamental and third harmonic, independently of the HB solver.
            coeffs[0][1] = rotation.scale(3.0 * radius.powi(3) * self.units[0] / 8.0)
                / C64::new(1.0, omega);
            coeffs[0][3] = (rotation * rotation * rotation)
                .scale(radius.powi(3) * self.units[0] / 8.0) / C64::new(1.0, 3.0 * omega);
        }
        HbOrbit { omega, coeffs,
            residual_trace: vec![123.0], residual: 123.0, aft_samples: 0, sampling_defect: 123.0 }
    }
    fn derivative(&self, _: f64, x: &[f64], out: &mut [f64]) -> bool {
        self.tangent_evaluations.set(self.tangent_evaluations.get() + 1);
        if matches!(self.tangent, Tangent::Refuse) { return false; }
        if matches!(self.tangent, Tangent::Incomplete) { out[0] = -1.0; return true; }
        let (x, y) = (x[1] / self.units[1], x[2] / self.units[2]);
        let a = x * x + y * y;
        let radial = self.parameter - (a - 1.0).powi(2);
        let omega = 1.0 + 0.2 * a;
        let radial_x = -4.0 * (a - 1.0) * x;
        let radial_y = -4.0 * (a - 1.0) * y;
        let entries = [
            -1.0, if self.cubic_observation { 3.0 * x * x } else { 0.0 }, 0.0,
            0.0, 0.4 * x * y + radial + x * radial_x,
            omega + 0.4 * y * y + x * radial_y,
            0.0, -omega - 0.4 * x * x + y * radial_x,
            -0.4 * x * y + radial + y * radial_y,
        ];
        for (index, entry) in entries.into_iter().enumerate() {
            out[index] = entry * self.equation_units[index / 3] / self.units[index % 3];
        }
        true
    }
}
impl OrbitProblem for Fold {
    fn dim(&self) -> usize { 3 }
    fn autonomous(&self) -> bool { !self.forced }
    fn island_state_jacobian(&self) -> Option<fs_orbit::IslandStateJacobian<Self>> {
        if matches!(self.tangent, Tangent::Absent) { None } else { Some(Self::derivative) }
    }
    fn island(&self, _: f64, x: &[f64], out: &mut [f64]) {
        self.evaluations.set(self.evaluations.get() + 1);
        let x: [f64; 3] = std::array::from_fn(|i| x[i] / self.units[i]);
        let a = x[1] * x[1] + x[2] * x[2];
        let radial = self.parameter - (a - 1.0).powi(2);
        let w = 1.0 + 0.2 * a;
        out[0] = (-(x[0] - 0.3) + if self.cubic_observation { x[1].powi(3) } else { 0.0 })
            * self.equation_units[0];
        out[1] = (w * x[2] + radial * x[1]) * self.equation_units[1];
        out[2] = (-w * x[1] + radial * x[2]) * self.equation_units[2];
    }
    fn port(&self, s: C64) -> Vec<C64> {
        let mut matrix = vec![C64::ZERO; 9];
        for i in 0..3 { matrix[i * 3 + i] = s.scale(self.equation_units[i] / self.units[i]); }
        matrix
    }
}
impl ContinuableProblem for Fold {
    fn set_parameter(&mut self, value: f64) { self.parameter = value; }
}
fn budget(model: &Fold) -> AutonomousContinuationBudget {
    AutonomousContinuationBudget {
        harmonics: 3, steps: 24, max_newton: 40, max_attempts_per_step: 12,
        step_sizes: [0.045, 1e-5, 0.08], tolerance: 2e-10, minimum_oscillation: 1e-5,
        maximum_unknowns: 64, maximum_payload_bytes: 16 << 20,
        maximum_residual_evaluations: 20_000,
        scales: AutonomousScales { states: model.units.to_vec(),
            equations: model.equation_units.to_vec(), angular_frequency: 2.0, parameter: 1.0 },
    }
}
fn lambda(radius: f64) -> f64 { (radius * radius - 1.0).powi(2) }

#[test]
fn autonomous_cycles_cross_a_parameter_fold_without_freezing_frequency_or_any_state() {
    let mut problem = Fold::new();
    let b = budget(&problem);
    let seeds = [problem.seed(0.65, 0.0, 3), problem.seed(0.70, 0.0, 3)];
    let branch = continue_autonomous_branch(&mut problem,
        [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &b, || false).unwrap();
    assert_eq!(branch.points.len(), 26);
    assert!(branch.residual_evaluations < b.maximum_residual_evaluations);
    let mut falling = false;
    let mut rising = false;
    for pair in branch.points.windows(2) {
        falling |= pair[1].parameter < pair[0].parameter;
        rising |= pair[1].parameter > pair[0].parameter;
    }
    assert!(falling && rising, "continuation must traverse the fold");
    for point in &branch.points {
        let radius = point.orbit.first_harmonic_amplitude(1);
        assert!((point.parameter - lambda(radius)).abs() < 2e-8);
        assert!((point.orbit.omega - (1.0 + 0.2 * radius * radius)).abs() < 2e-8);
        assert!(point.scaled_residual <= b.tolerance);
        assert!(point.scaled_sampling_defect <= b.tolerance);
        assert!(point.phase_residual.abs() <= b.tolerance);
        assert!(point.arclength_residual.abs() <= b.tolerance);
        assert!(point.orbit.aft_samples > 16);
        assert!(point.orbit.residual < 1e-8);
        assert_eq!(point.orbit.residual_trace.last(), Some(&point.orbit.residual));
        assert!(point.orbit.first_harmonic_amplitude(0) < 1e-10);
    }
    assert!(branch.points.last().unwrap().orbit.first_harmonic_amplitude(1) > 1.4);
    assert_eq!(problem.parameter, branch.points.last().unwrap().parameter);
}

#[test]
fn independent_state_and_equation_units_and_common_phase_preserve_the_physical_branch() {
    let mut reference = None;
    for mode in 0..3 {
        let mut problem = Fold::new();
        if mode == 2 {
            problem.units = [1e5, 1e-3, 1e3];
            problem.equation_units = [7.0, 1e4, 1e-4];
        }
        let b = budget(&problem);
        let phase = if mode == 1 { 0.61 } else { 0.0 };
        let seeds = [problem.seed(0.65, phase, 3), problem.seed(0.70, phase, 3)];
        let branch = continue_autonomous_branch(&mut problem,
            [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &b, || false).unwrap();
        let values: Vec<_> = branch.points.iter().map(|p| (p.parameter, p.orbit.omega)).collect();
        if let Some(previous) = &reference {
            let previous: &Vec<(f64, f64)> = previous;
            for ((p, w), (q, v)) in values.iter().zip(previous) {
                assert!((p-q).abs() < 2e-7 && (w-v).abs() < 2e-7);
            }
        } else { reference = Some(values); }
    }
}

#[test]
fn exact_zero_harmonic_padding_does_not_change_the_parseval_arclength_metric() {
    let mut reference: Option<Vec<(f64, f64)>> = None;
    for harmonics in [1, 3] {
        let mut problem = Fold::new();
        let mut b = budget(&problem);
        b.harmonics = harmonics;
        let seeds = [problem.seed(0.65, 0.0, harmonics), problem.seed(0.70, 0.0, harmonics)];
        let branch = continue_autonomous_branch(&mut problem,
            [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &b, || false).unwrap();
        let values: Vec<_> = branch.points.iter().map(|p| (p.parameter, p.orbit.omega)).collect();
        if let Some(previous) = &reference {
            for ((p, w), (q, v)) in values.iter().zip(previous) {
                assert!((p-q).abs() < 2e-7 && (w-v).abs() < 2e-7);
            }
        } else { reference = Some(values); }
    }
}

#[test]
fn resource_and_shape_refusals_precede_physics_and_equilibrium_is_not_a_cycle() {
    let mut problem = Fold::new();
    let base = budget(&problem);
    let seeds = [problem.seed(0.65, 0.0, 3), problem.seed(0.70, 0.0, 3)];
    for case in 0..6 {
        let mut b = base.clone();
        match case {
            0 => b.maximum_payload_bytes = 1,
            1 => b.maximum_unknowns = 1,
            2 => b.step_sizes[1] = 1.0,
            3 => b.scales.equations[2] = f64::NAN,
            4 => b.steps = usize::MAX,
            _ => b.harmonics = 2,
        }
        assert!(matches!(continue_autonomous_branch(&mut problem,
            [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &b, || false),
            Err(AutonomousContinuationError::Invalid(_))));
        assert_eq!(problem.evaluations.get(), 0);
    }
    problem.forced = true;
    assert!(continue_autonomous_branch(&mut problem,
        [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &base, || false).is_err());
    assert_eq!(problem.evaluations.get(), 0);
    problem.forced = false;
    let rest = problem.seed(0.0, 0.0, 3);
    assert!(matches!(continue_autonomous_branch(&mut problem,
        [(0.1, &rest), (0.2, &rest)], &base, || false),
        Err(AutonomousContinuationError::Numerical(fs_orbit::OrbitError::TrivialCollapse))));
}

#[test]
fn cancellation_and_evaluation_exhaustion_restore_the_last_accepted_parameter() {
    let mut problem = Fold::new();
    let mut b = budget(&problem);
    let seeds = [problem.seed(0.65, 0.0, 3), problem.seed(0.70, 0.0, 3)];
    // Each exact seed needs its original residual and one independent grid.
    b.maximum_residual_evaluations = 5;
    assert!(matches!(continue_autonomous_branch(&mut problem,
        [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &b, || false),
        Err(AutonomousContinuationError::EvaluationLimit { evaluations: 5 })));
    assert!((problem.parameter - lambda(0.70)).abs() < 1e-14);
    b.maximum_residual_evaluations = 20_000;
    let mut polls = 0;
    assert!(matches!(continue_autonomous_branch(&mut problem,
        [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])], &b, || {
            polls += 1; polls >= 18
        }), Err(AutonomousContinuationError::Cancelled)));
    assert!((problem.parameter - lambda(0.70)).abs() < 1e-14);
}

#[test]
fn exact_tangents_preserve_generated_harmonics_through_a_fold_with_fewer_source_calls() {
    for mixed_units in [false, true] {
        let mut exact = Fold::new();
        exact.tangent = Tangent::Exact;
        exact.cubic_observation = true;
        if mixed_units {
            exact.units = [1e5, 1e-3, 1e3];
            exact.equation_units = [7.0, 1e4, 1e-4];
        }
        let mut probed = Fold::new();
        probed.units = exact.units;
        probed.equation_units = exact.equation_units;
        probed.cubic_observation = true;
        let b = budget(&exact);
        let seeds = [exact.seed(0.65, 0.37, 3), exact.seed(0.70, 0.37, 3)];
        let initial = [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])];
        let analytic = continue_autonomous_branch(&mut exact, initial, &b, || false).unwrap();
        let finite = continue_autonomous_branch(&mut probed, initial, &b, || false).unwrap();
        assert_eq!(analytic.points.len(), finite.points.len());
        println!(
            "mixed_units={mixed_units}: source calls exact={} finite_difference={}; residual calls exact={} finite_difference={}",
            exact.evaluations.get(), probed.evaluations.get(),
            analytic.residual_evaluations, finite.residual_evaluations
        );
        assert!(exact.tangent_evaluations.get() > 0);
        assert_eq!(probed.tangent_evaluations.get(), 0);
        assert!(
            3 * exact.evaluations.get() < probed.evaluations.get(),
            "source calls: exact {}, probed {}",
            exact.evaluations.get(),
            probed.evaluations.get()
        );
        assert!(3 * analytic.residual_evaluations < finite.residual_evaluations);
        let mut falling = false;
        let mut rising = false;
        for pair in analytic.points.windows(2) {
            falling |= pair[1].parameter < pair[0].parameter;
            rising |= pair[1].parameter > pair[0].parameter;
        }
        assert!(falling && rising);
        for (point, twin) in analytic.points.iter().zip(&finite.points) {
            let radius = point.orbit.first_harmonic_amplitude(1) / exact.units[1];
            assert!((point.parameter - lambda(radius)).abs() < 2e-8);
            assert!((point.orbit.omega - (1.0 + 0.2 * radius * radius)).abs() < 2e-8);
            assert!((point.parameter - twin.parameter).abs() < 2e-7);
            assert!((point.orbit.omega - twin.orbit.omega).abs() < 2e-7);
            for k in 0..3 {
                for n in 0..=3 {
                    assert!(
                        (point.orbit.coeffs[k][n] - twin.orbit.coeffs[k][n])
                            .scale(1.0 / exact.units[k])
                            .abs()
                            < 2e-7
                    );
                }
            }
            let x = point.orbit.coeffs[1][1].scale(1.0 / exact.units[1]);
            let expected_fundamental =
                x.scale(3.0 * x.abs().powi(2)) / C64::new(1.0, point.orbit.omega);
            let expected_third = (x * x * x) / C64::new(1.0, 3.0 * point.orbit.omega);
            assert!(expected_third.abs() > 1e-4);
            assert!(
                (point.orbit.coeffs[0][1].scale(1.0 / exact.units[0]) - expected_fundamental).abs()
                    < 2e-8
            );
            assert!(
                (point.orbit.coeffs[0][3].scale(1.0 / exact.units[0]) - expected_third).abs()
                    < 2e-8
            );
            assert!(point.scaled_residual <= b.tolerance);
            assert!(point.scaled_sampling_defect <= b.tolerance);
        }
    }
}

#[test]
fn declared_autonomous_tangent_refusal_never_becomes_a_residual_probe() {
    for mode in [Tangent::Refuse, Tangent::Incomplete] {
        let mut problem = Fold::new();
        problem.tangent = mode;
        let b = budget(&problem);
        let mut seeds = [problem.seed(0.65, 0.0, 3), problem.seed(0.70, 0.0, 3)];
        seeds[0].omega *= 1.03;
        assert_eq!(
            continue_autonomous_branch(
                &mut problem,
                [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])],
                &b,
                || false
            )
            .unwrap_err(),
            AutonomousContinuationError::Numerical(fs_orbit::OrbitError::TangentRefused {
                sample: 0
            })
        );
        assert_eq!(problem.tangent_evaluations.get(), 1);
        assert_eq!(
            problem.evaluations.get(),
            16,
            "only the original base waveform was evaluated"
        );
        assert_eq!(problem.parameter, lambda(0.65));
    }
}

#[test]
fn exact_autonomous_tangents_respect_payload_and_original_residual_limits() {
    let mut problem = Fold::new();
    problem.tangent = Tangent::Exact;
    let mut b = budget(&problem);
    let seeds = [problem.seed(0.65, 0.0, 3), problem.seed(0.70, 0.0, 3)];
    let initial = [(lambda(0.65), &seeds[0]), (lambda(0.70), &seeds[1])];
    // Enough for the original residual-probe corrector, but not its additional
    // optional tangent tables. Refuse before evaluating or allocating that path.
    b.maximum_payload_bytes = 128 << 10;
    assert!(matches!(
        continue_autonomous_branch(&mut problem, initial, &b, || false),
        Err(AutonomousContinuationError::Invalid(
            "analytic continuation tangent payload exceeds its allowance"
        ))
    ));
    assert_eq!(problem.tangent_evaluations.get(), 0);
    assert_eq!(problem.parameter, lambda(0.70));
    problem.tangent = Tangent::Absent;
    continue_autonomous_branch(&mut problem, initial, &b, || false).unwrap();
    problem.tangent = Tangent::Exact;
    b.maximum_payload_bytes = 16 << 20;
    b.maximum_residual_evaluations = 5;
    assert!(matches!(
        continue_autonomous_branch(&mut problem, initial, &b, || false),
        Err(AutonomousContinuationError::EvaluationLimit { evaluations: 5 })
    ));
    assert!(problem.tangent_evaluations.get() > 0);
    assert_eq!(problem.parameter, lambda(0.70));
}
