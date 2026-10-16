//! Real preload responses at contact onset, shared load cases and signed maps.
use fs_couple::modal_acoustic_time::{ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{ModalAttachment, ModalConnection, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::SensitivityBudget;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::DisplacementTarget;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::*;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::constraints::*;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

fn map(shape: f64) -> ModalAttachment { ModalAttachment { component: 0, shapes: vec![shape] } }
fn problem(forces: &[f64], derivative_budget: usize) -> EquilibriumDesign {
    let budget = DesignBudget {
        coupling: ModalCouplingConfig { max_modes: 4, max_connections: 4, max_setup_terms: 16384,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 1000.0,
            maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 10000.0,
            solve_relative_tolerance: 1e-10, energy_absolute_tolerance_j: 1e-11, energy_relative_tolerance: 1e-9 },
        contact: MultiContactConfig { max_contacts: 4, max_sweeps: 128, max_setup_terms: 16384 },
        sensitivity: SensitivityBudget { max_contacts: 4, max_setup_terms: derivative_budget,
            max_query_terms: derivative_budget, minimum_contact_margin_m: 1e-6 },
        max_cases: 4, max_variables: 4, max_bindings: 8, max_ports_per_case: 8,
    };
    let body = ModalAcousticTimeModel::try_free_mass(48000, 0.25, 0.0, 0.0,
        ModalAcousticTimeBudget::audible_reference()).unwrap();
    let spring = ModalConnection { left: map(2.0), right: map(0.0),
        stiffness_n_m: 64.0, damping_n_s_m: 0.0, rest_extension_m: 0.0 };
    let contact = (ModalContact { left: map(2.0), right: map(0.0),
        law: Obstacle::new(vec![-1.0], 1, 1, vec![1.0/64.0], vec![1.0], 128.0, 1.0,
            "synthetic-forward-response".into()).unwrap() },
        ModalContactConfig { max_iterations: 128, maximum_force_n: 10000.0,
            maximum_penetration_m: 0.1, force_absolute_tolerance_n: 1e-12, force_relative_tolerance: 1e-13 });
    let cases = forces.iter().enumerate().map(|(i, &force)| DesignLoadCase {
        name: format!("load-{i}"), loads: vec![DesignLoad { attachment: map(2.0), force_n: force }],
        targets: vec![DisplacementTarget { attachment: map(2.0), target_m: 0.0, scale_m: 0.01, weight: 1.0 }],
    }).collect();
    let variables = vec![DesignVariable { name: "support".into(), reference: 64.0, scale: 64.0,
        minimum: 32.0, maximum: 128.0, fields: vec![DesignField::SpringStiffness(0)] }];
    EquilibriumDesign::new(vec![body], vec![spring], vec![contact], cases, variables, budget,
        &CancelGate::new()).unwrap()
}
fn requirements() -> Vec<ResponseConstraint> {
    [
        (2, ResponseQuantity::ContactForce(0), ConstraintSense::AtMost, 0.1, 0.5),
        (0, ResponseQuantity::Displacement(map(-6.0)), ConstraintSense::AtLeast, -0.02, 0.01),
        (1, ResponseQuantity::SpringForce(0), ConstraintSense::AtLeast, -0.9, 2.0),
        (2, ResponseQuantity::ContactPenetration(0), ConstraintSense::AtMost, 0.001, 0.01),
        (0, ResponseQuantity::ContactForce(0), ConstraintSense::AtMost, 0.0, 1.0),
        (1, ResponseQuantity::Displacement(map(2.0)), ConstraintSense::Equal, 1.0/64.0, 1.0),
    ].into_iter().enumerate().map(|(i, (case, quantity, sense, bound, scale))|
        ResponseConstraint { name: format!("response-{i}"), case, quantity, sense, bound, scale }).collect()
}
fn close(a: f64, b: f64) { assert!((a-b).abs() < 2e-10, "{a} != {b}"); }

#[test]
fn all_quantities_senses_and_case_order_work_at_a_contact_kink_without_derivatives() {
    let gate = CancelGate::new();
    let p = problem(&[0.5, 1.0, 1.5], 0).with_constraints(requirements(), 6, &gate).unwrap();
    let mut work = DesignControl::new(2, 6);
    let result = p.evaluate_forward_with_constraints(&[0.0], &mut work, &gate).unwrap();
    assert_eq!(result.constraints.len(), 6);
    assert_eq!(result.observations.unassessed_response_constraints, 0);
    let expected = [1.0/3.0, -3.0/128.0, -1.0, 1.0/384.0, 0.0, 1.0/64.0];
    for ((row, requirement), value) in result.constraints.iter().zip(p.constraints()).zip(expected) {
        close(row.value, value);
        let sign = if requirement.sense == ConstraintSense::AtLeast { -1.0 } else { 1.0 };
        close(row.residual, sign*(value-requirement.bound)/requirement.scale);
    }
    assert!(result.constraints[0].residual > 0.0); // a failed requirement is retained
    assert_eq!(work.work(), DesignWork { evaluations: 1, case_solves: 3 });
    assert!(p.evaluate(&[0.0], &mut work, &gate).is_err()); // no derivative admission weakened
}

#[test]
fn forward_and_adjoint_response_values_agree_away_from_switches() {
    let gate = CancelGate::new();
    let p = problem(&[0.5, 1.0, 1.5], 16384).with_constraints(requirements(), 6, &gate).unwrap();
    let mut work = DesignControl::new(8, 24);
    for x in [-0.25, 0.25, 0.75] {
        let primal = p.evaluate_forward_with_constraints(&[x], &mut work, &gate).unwrap();
        let differentiated = p.evaluate(&[x], &mut work, &gate).unwrap();
        assert_eq!(primal.observations.physical_parameters, differentiated.physical_parameters);
        for (a, b) in primal.constraints.iter().zip(&differentiated.constraints) {
            close(a.value, b.value); close(a.residual, b.residual);
        }
    }
}

#[test]
fn legacy_observations_are_unchanged_and_do_not_claim_constraints_were_assessed() {
    let gate = CancelGate::new();
    let p = problem(&[0.5, 1.0, 1.5], 0).with_constraints(requirements(), 6, &gate).unwrap();
    let mut work = DesignControl::new(8, 24);
    let old = p.evaluate_forward(&[0.0], &mut work, &gate).unwrap();
    let new = p.evaluate_forward_with_constraints(&[0.0], &mut work, &gate).unwrap();
    assert_eq!(old.unassessed_response_constraints, 6);
    assert_eq!(old.cases, new.observations.cases);
    assert_eq!(old.physical_parameters, new.observations.physical_parameters);
    assert_eq!(new, p.evaluate_forward_with_constraints(&[0.0], &mut work, &gate).unwrap());
    let empty = problem(&[0.5], 0).evaluate_forward_with_constraints(&[0.0],
        &mut DesignControl::new(1, 1), &gate).unwrap();
    assert!(empty.constraints.is_empty());
}

#[test]
fn late_physical_failure_cancellation_and_case_budgets_never_publish_partial_responses() {
    let gate = CancelGate::new();
    let p = problem(&[0.5, 1.0, 1e9], 0).with_constraints(requirements(), 6, &gate).unwrap();
    let mut work = DesignControl::new(8, 2);
    assert!(matches!(p.evaluate_forward_with_constraints(&[0.0], &mut work, &gate), Err(DesignError::Budget { .. })));
    assert_eq!(work.work(), DesignWork { evaluations: 1, case_solves: 0 });
    work.extend(8, 24).unwrap();
    assert!(matches!(p.evaluate_forward_with_constraints(&[0.0], &mut work, &gate), Err(DesignError::Case { case: 2, .. })));
    assert_eq!(work.work(), DesignWork { evaluations: 2, case_solves: 3 });
    let cancelled = CancelGate::new(); cancelled.request();
    assert!(matches!(p.evaluate_forward_with_constraints(&[0.0], &mut work, &cancelled), Err(DesignError::Cancelled)));
    assert_eq!(work.work(), DesignWork { evaluations: 2, case_solves: 3 });
}

#[test]
fn unrepresentable_constraint_residual_refuses_without_poisoning_the_physical_template() {
    let gate = CancelGate::new();
    let mut rows = requirements();
    rows[0].bound = f64::MAX;
    rows[0].scale = f64::MIN_POSITIVE;
    let p = problem(&[0.5, 1.0, 1.5], 0).with_constraints(rows, 6, &gate).unwrap();
    let mut work = DesignControl::new(8, 24);
    assert!(matches!(p.evaluate_forward_with_constraints(&[0.0], &mut work, &gate), Err(DesignError::Case { case: 2, .. })));
    let retained = p.evaluate_forward(&[0.0], &mut work, &gate).unwrap();
    assert_eq!(retained.cases.len(), 3);
    assert_eq!(retained.unassessed_response_constraints, 6);
    assert_eq!(work.work(), DesignWork { evaluations: 2, case_solves: 6 });
}
