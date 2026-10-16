//! Forward observations are legal at contact onset without licensing a derivative.
use fs_couple::modal_acoustic_time::{ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{ModalAttachment, ModalConnection, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::SensitivityBudget;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::DisplacementTarget;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::*;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

fn map(shape: f64) -> ModalAttachment { ModalAttachment { component: 0, shapes: vec![shape] } }
fn problem(forces: &[f64], derivative_budget: usize, contacted: bool) -> EquilibriumDesign {
    let coupling = ModalCouplingConfig {
        max_modes: 4, max_connections: 4, max_setup_terms: 16384,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 1000.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 10000.0,
        solve_relative_tolerance: 1e-10, energy_absolute_tolerance_j: 1e-11,
        energy_relative_tolerance: 1e-9,
    };
    let budget = DesignBudget {
        coupling, contact: MultiContactConfig { max_contacts: 4, max_sweeps: 128, max_setup_terms: 16384 },
        sensitivity: SensitivityBudget { max_contacts: 4, max_setup_terms: derivative_budget,
            max_query_terms: derivative_budget, minimum_contact_margin_m: 1e-6 },
        max_cases: 4, max_variables: 4, max_bindings: 8, max_ports_per_case: 8,
    };
    // Nonunit mass/shape: Q = sqrt(m)*q and shape*Q = physical q.
    let body = ModalAcousticTimeModel::try_free_mass(48000, 0.25, 0.0, 0.0,
        ModalAcousticTimeBudget::audible_reference()).unwrap();
    let spring = ModalConnection { left: map(2.0), right: map(0.0),
        stiffness_n_m: 64.0, damping_n_s_m: 0.0, rest_extension_m: 0.0 };
    let contacts = if contacted { vec![(ModalContact { left: map(2.0), right: map(0.0),
        law: Obstacle::new(vec![-1.0], 1, 1, vec![1.0/64.0], vec![1.0], 128.0, 1.0,
            "synthetic-linear-contact".into()).unwrap() },
        ModalContactConfig { max_iterations: 128, maximum_force_n: 10000.0,
            maximum_penetration_m: 0.1, force_absolute_tolerance_n: 1e-12,
            force_relative_tolerance: 1e-13 })] } else { vec![] };
    let cases = forces.iter().enumerate().map(|(i, &force)| DesignLoadCase {
        name: format!("load-{i}"), loads: vec![DesignLoad { attachment: map(2.0), force_n: force }],
        targets: vec![DisplacementTarget { attachment: map(2.0), target_m: 0.0,
            scale_m: 0.01, weight: 1.0 }, DisplacementTarget { attachment: map(-6.0),
            target_m: 0.0, scale_m: 1.0, weight: 0.0 }],
    }).collect();
    let variables = vec![DesignVariable { name: "support".into(), reference: 64.0, scale: 64.0,
        minimum: 32.0, maximum: 128.0, fields: vec![DesignField::SpringStiffness(0)] }];
    EquilibriumDesign::new(vec![body], vec![spring], contacts, cases, variables, budget,
        &CancelGate::new()).unwrap()
}
fn displacement(force: f64, stiffness: f64, contacted: bool) -> f64 {
    if !contacted || force <= stiffness / 64.0 { force / stiffness }
    else { (force + 2.0) / (stiffness + 128.0) }
}

#[test]
fn primal_observations_cross_a_kink_without_weakening_derivative_admission() {
    let forces = [0.5, 1.0, 1.5];
    let p = problem(&forces, 16384, true);
    let gate = CancelGate::new();
    let mut work = DesignControl::new(16, 48);
    let result = p.evaluate_forward(&[0.0], &mut work, &gate).unwrap();
    assert_eq!(result.physical_parameters, vec![64.0]);
    assert_eq!(result.cases.len(), 3);
    for (case, force) in result.cases.iter().zip(forces) {
        let q = displacement(force, 64.0, true);
        assert!((case.observations_m[0] - q).abs() < 1e-11);
        assert!((case.observations_m[1] + 3.0*q).abs() < 1e-11);
        let energy = 0.5*64.0*q*q + 0.5*128.0*(q-1.0/64.0).max(0.0).powi(2);
        assert!((case.stored_energy_j - energy).abs() < 1e-11);
    }
    assert_eq!(work.work(), DesignWork { evaluations: 1, case_solves: 3 });
    assert!(matches!(p.evaluate(&[0.0], &mut work, &gate), Err(DesignError::Case { case: 1, .. })));
    assert_eq!(work.work(), DesignWork { evaluations: 2, case_solves: 5 });
    assert_eq!(p.evaluate_forward(&[0.0], &mut work, &gate).unwrap(), result);
}

#[test]
fn forward_and_adjoint_paths_observe_identical_settled_states_away_from_switches() {
    let gate = CancelGate::new();
    for contacted in [false, true] {
        let p = problem(&[0.5, 1.0, 1.5], 16384, contacted);
        let mut control = DesignControl::new(16, 48);
        for x in [-0.25, 0.25, 0.75] {
            let primal = p.evaluate_forward(&[x], &mut control, &gate).unwrap();
            let adjoint = p.evaluate(&[x], &mut control, &gate).unwrap();
            assert_eq!(primal.physical_parameters, adjoint.physical_parameters);
            for (a, b) in primal.cases.iter().zip(&adjoint.cases) {
                assert_eq!(a.observations_m, b.observations_m);
                assert!(b.adjoint_relative_residual <= 1e-10);
            }
        }
    }
}

#[test]
fn an_observation_does_not_require_a_derivative_setup_or_query_budget() {
    let p = problem(&[0.5], 0, true);
    let mut control = DesignControl::new(4, 4);
    assert!(p.evaluate_forward(&[0.0], &mut control, &CancelGate::new()).is_ok());
    assert!(matches!(p.evaluate(&[0.0], &mut control, &CancelGate::new()), Err(DesignError::Case { .. })));
}

#[test]
fn cancelled_out_of_domain_and_insufficient_case_allowances_preserve_work() {
    let p = problem(&[0.5, 1.0, 1.5], 16384, true);
    let gate = CancelGate::new();
    let mut work = DesignControl::new(8, 2);
    assert!(matches!(p.evaluate_forward(&[0.0], &mut work, &gate), Err(DesignError::Budget { .. })));
    assert_eq!(work.work(), DesignWork { evaluations: 1, case_solves: 0 });
    work.extend(8, 12).unwrap();
    assert!(matches!(p.evaluate_forward(&[-1.0], &mut work, &gate), Err(DesignError::OutsideBounds { .. })));
    let before = work.work();
    let cancelled = CancelGate::new(); cancelled.request();
    assert!(matches!(p.evaluate_forward(&[0.0], &mut work, &cancelled), Err(DesignError::Cancelled)));
    assert_eq!(work.work(), before);
    let result = p.evaluate_forward(&[0.0], &mut work, &gate).unwrap();
    assert_eq!(result, p.evaluate_forward(&[0.0], &mut work, &gate).unwrap());
    assert_eq!(work.work(), DesignWork { evaluations: 4, case_solves: 6 });
}

#[test]
fn late_physical_failure_returns_no_partial_family_and_keeps_attempts_charged() {
    let p = problem(&[0.5, 1.0, 1e9], 16384, true);
    let mut control = DesignControl::new(4, 12);
    for attempt in 1..=2 {
        assert!(matches!(p.evaluate_forward(&[0.0], &mut control, &CancelGate::new()),
            Err(DesignError::Case { case: 2, .. })));
        assert_eq!(control.work(), DesignWork { evaluations: attempt, case_solves: 3*attempt });
    }
}

#[test]
fn forward_observations_do_not_claim_that_declared_response_constraints_pass() {
    use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::constraints::{
        ConstraintSense, ResponseConstraint, ResponseQuantity,
    };
    let gate = CancelGate::new();
    let p = problem(&[0.5], 16384, true).with_constraints(vec![ResponseConstraint {
        name: "nonnegative-left-reaction".into(), case: 0, quantity: ResponseQuantity::SpringForce(0),
        sense: ConstraintSense::AtLeast, bound: 0.0, scale: 1.0,
    }], 1, &gate).unwrap();
    let mut work = DesignControl::new(4, 4);
    let forward = p.evaluate_forward(&[0.0], &mut work, &gate).unwrap();
    assert_eq!(forward.unassessed_response_constraints, 1);
    let full = p.evaluate(&[0.0], &mut work, &gate).unwrap();
    assert_eq!(full.constraints.len(), 1);
    assert!(full.constraints[0].residual > 0.0); // violated is not dropped or declared feasible
    assert_eq!(forward.cases[0].observations_m, full.cases[0].observations_m);
}
