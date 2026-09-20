use fs_ascent::equilibrium::{EquilibriumStudy, EquilibriumStudyError};
use fs_ascent::{SqpError, SqpStop};
use fs_couple::modal_acoustic_time::{ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{ModalAttachment, ModalConnection, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::SensitivityBudget;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::DisplacementTarget;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::{
    DesignBudget, DesignControl, DesignError, DesignField, DesignLoad, DesignLoadCase,
    DesignVariable, DesignWork, EquilibriumDesign,
};
use fs_exec::CancelGate;

const FORCES: [f64; 3] = [1.0, 2.0, 4.0];
const OBSERVATION_SCALE: f64 = 0.002;

fn map() -> ModalAttachment { ModalAttachment { component: 0, shapes: vec![1.0] } }

fn problem(target_stiffness: f64, lower: f64, upper: f64) -> EquilibriumDesign {
    let budget = DesignBudget {
        coupling: ModalCouplingConfig {
            max_modes: 4, max_connections: 4, max_setup_terms: 16384,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 1000.0,
            maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 10000.0,
            solve_relative_tolerance: 1e-10, energy_absolute_tolerance_j: 1e-11,
            energy_relative_tolerance: 1e-9,
        },
        contact: MultiContactConfig { max_contacts: 4, max_sweeps: 128, max_setup_terms: 16384 },
        sensitivity: SensitivityBudget { max_contacts: 4, max_setup_terms: 16384,
            max_query_terms: 16384, minimum_contact_margin_m: 1e-9 },
        max_cases: 4, max_variables: 8, max_bindings: 16, max_ports_per_case: 8,
    };
    let body = ModalAcousticTimeModel::try_free_mass(48000, 1.0, 0.0, 0.0,
        ModalAcousticTimeBudget::audible_reference()).unwrap();
    let spring = ModalConnection {
        left: map(), right: ModalAttachment { component: 0, shapes: vec![0.0] },
        stiffness_n_m: 400.0, damping_n_s_m: 0.0, rest_extension_m: 0.0,
    };
    let cases = FORCES.into_iter().enumerate().map(|(i, force)| DesignLoadCase {
        name: format!("force-{i}"),
        loads: vec![DesignLoad { attachment: map(), force_n: force }],
        targets: vec![DisplacementTarget { attachment: map(),
            target_m: force / target_stiffness, scale_m: OBSERVATION_SCALE, weight: 1.0 }],
    }).collect();
    let variables = vec![DesignVariable {
        name: "support-stiffness".into(), reference: 400.0, scale: 400.0,
        minimum: lower, maximum: upper, fields: vec![DesignField::SpringStiffness(0)],
    }];
    EquilibriumDesign::new(vec![body], vec![spring], vec![], cases, variables,
        budget, &CancelGate::new()).unwrap()
}

fn analytic(k: f64, target: f64) -> f64 {
    FORCES.iter().map(|force| {
        0.5 * ((force / k - force / target) / OBSERVATION_SCALE).powi(2)
    }).sum()
}

#[test]
fn fits_shared_physical_stiffness_to_three_independent_experiments() {
    let problem = problem(750.0, 100.0, 1000.0);
    let mut control = DesignControl::new(512, 1536);
    let gate = CancelGate::new();
    let mut study = EquilibriumStudy::new(&problem, &[0.0], &mut control, 3, &gate).unwrap();
    let result = study.run(1e-9, 128, 512, &gate).unwrap();
    assert_eq!(result.stop, SqpStop::Converged);
    assert!(result.solution.kkt.within_tolerance(1e-9));
    assert!((study.accepted().physical_parameters[0] - 750.0).abs() < 1e-5);
    assert_eq!(study.accepted().cases.len(), FORCES.len());
    for (case, force) in study.accepted().cases.iter().zip(FORCES) {
        assert!((case.observations_m[0] - force / 750.0).abs() < 1e-10);
        assert!(case.adjoint_relative_residual <= 1e-10);
    }
    assert_eq!(result.solution.f.to_bits(), study.accepted().value.to_bits());
    assert_eq!(study.work().evaluations, result.solution.evals);
    assert!(study.work().case_solves <= 3 * study.work().evaluations);
    assert!(study.optimizer().history().windows(2).all(|v| v[1] < v[0]));
}

#[test]
fn both_active_bounds_get_kkt_certificates_not_raw_gradient_convergence() {
    for (target, expected, active) in [(1500.0, 1000.0, 1), (50.0, 100.0, 0)] {
        let problem = problem(target, 100.0, 1000.0);
        let mut control = DesignControl::new(512, 1536);
        let gate = CancelGate::new();
        let mut study = EquilibriumStudy::new(&problem, &[0.0], &mut control, 3, &gate).unwrap();
        let result = study.run(1e-8, 128, 512, &gate).unwrap();
        assert_eq!(result.stop, SqpStop::Converged, "target={target}: {result:?}");
        assert!(result.solution.kkt.within_tolerance(1e-8));
        assert!((study.accepted().physical_parameters[0] - expected).abs() < 1e-5);
        assert!(study.accepted().gradient[0].abs() > 1e-3);
        assert!(result.solution.nu[active] > 1e-3);
        assert!(result.solution.nu[1 - active].abs() < 1e-8);
        assert!((result.solution.f - analytic(expected, target)).abs() < 1e-6);
        assert!(study.optimizer().sample().ci.iter().all(|c| *c <= 0.0));
    }
}

#[test]
fn accepted_step_splits_and_pre_cancel_preserve_physical_evidence_and_work() {
    let problem = problem(750.0, 100.0, 1000.0);
    let gate = CancelGate::new();
    let mut full_control = DesignControl::new(512, 1536);
    let mut full = EquilibriumStudy::new(&problem, &[0.0], &mut full_control, 3, &gate).unwrap();
    full.run(1e-9, 128, 512, &gate).unwrap();
    let mut split_control = DesignControl::new(512, 1536);
    let mut split = EquilibriumStudy::new(&problem, &[0.0], &mut split_control, 3, &gate).unwrap();
    split.run(1e-9, 2, 512, &gate).unwrap();
    let before = split.accepted().clone();
    let point = split.optimizer().point().to_vec();
    let work = split.work();
    let cancelled = CancelGate::new(); cancelled.request();
    assert!(matches!(split.run(1e-9, 128, 512, &cancelled), Err(SqpError::Cancelled)));
    assert_eq!(split.accepted(), &before);
    assert_eq!(split.optimizer().point(), point.as_slice());
    assert_eq!(split.work(), work);
    split.run(1e-9, 128, 512, &gate).unwrap();
    assert_eq!(full.optimizer().point(), split.optimizer().point());
    assert_eq!(full.optimizer().history(), split.optimizer().history());
    assert_eq!(full.accepted(), split.accepted());
    assert_eq!(full.work(), split.work());
}

#[test]
fn exhausted_physics_keeps_the_accepted_checkpoint_and_can_extend_without_refunds() {
    let problem = problem(750.0, 100.0, 1000.0);
    let gate = CancelGate::new();
    let mut control = DesignControl::new(512, 3);
    let mut study = EquilibriumStudy::new(&problem, &[0.0], &mut control, 3, &gate).unwrap();
    let accepted = study.accepted().clone();
    assert!(matches!(study.run(1e-9, 128, 512, &gate),
        Err(SqpError::Evaluation(DesignError::Budget { .. }))));
    assert_eq!(study.accepted(), &accepted);
    assert_eq!(study.optimizer().point(), &[0.0]);
    assert_eq!(study.work().case_solves, 3);
    assert!(study.work().evaluations >= 2);
    let spent = study.work();
    assert!(study.extend_physics_budget(1, 1).is_err());
    study.extend_physics_budget(512, 1536).unwrap();
    assert_eq!(study.work(), spent);
    let result = study.run(1e-9, 128, 512, &gate).unwrap();
    assert_eq!(result.stop, SqpStop::Converged);
    assert!(study.work().evaluations > spent.evaluations);
}

#[test]
fn callback_limit_is_hard_and_initialization_failures_leave_work_visible() {
    let problem = problem(750.0, 100.0, 1000.0);
    let gate = CancelGate::new();
    let mut control = DesignControl::new(512, 1536);
    {
        let mut study = EquilibriumStudy::new(&problem, &[0.0], &mut control, 3, &gate).unwrap();
        let result = study.run(1e-9, 128, 1, &gate).unwrap();
        assert_eq!(result.stop, SqpStop::EvaluationLimit);
        assert_eq!(result.solution.evals, 1);
        assert_eq!(study.work(), DesignWork { evaluations: 1, case_solves: 3 });
    }
    let result = EquilibriumStudy::new(&problem, &[-2.0], &mut control, 3, &gate);
    assert!(matches!(result, Err(SqpError::Evaluation(DesignError::OutsideBounds { .. }))));
    assert_eq!(control.work(), DesignWork { evaluations: 2, case_solves: 3 });
}

#[test]
fn dense_admission_and_bad_dimensions_fail_before_any_physics() {
    let problem = problem(750.0, 100.0, 1000.0);
    let gate = CancelGate::new();
    let mut control = DesignControl::new(1, 3);
    assert!(matches!(EquilibriumStudy::new(&problem, &[0.0], &mut control, 2, &gate),
        Err(SqpError::Invalid(_))));
    assert!(matches!(EquilibriumStudy::new(&problem, &[], &mut control, 3, &gate),
        Err(SqpError::Shape { .. })));
    assert_eq!(control.work(), DesignWork::default());
    // The public error keeps its underlying physics type, including sources.
    let _: Option<EquilibriumStudyError> = None;
}
