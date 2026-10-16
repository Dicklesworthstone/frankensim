//! Independent physical-coordinate checks of complete calibration Jacobians.
use fs_couple::modal_acoustic_time::{ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{ModalAttachment, ModalConnection, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::SensitivityBudget;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::DisplacementTarget;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::*;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::observations::*;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

fn map(shape: f64) -> ModalAttachment { ModalAttachment { component: 0, shapes: vec![shape] } }
fn displacement(x: &[f64], force: f64) -> f64 {
    let support = 128.0 * (1.0 + x[0]); // two springs share BOTH stiffness and rest
    let rest = 0.01 * x[1];
    let contact = 128.0 * (1.0 + x[2]) * (1.0 + 0.5 * x[4]);
    let gap = (1.0 + x[3]) / 64.0;
    let free = force / support + rest;
    if free <= gap { free } else { (force + support * rest + contact * gap) / (support + contact) }
}
fn build(forces: [f64; 2], derivative_budget: usize) -> EquilibriumDesign {
    let budget = DesignBudget {
        coupling: ModalCouplingConfig { max_modes: 4, max_connections: 4, max_setup_terms: 16384,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 1000.0, maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 10000.0, solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-11, energy_relative_tolerance: 1e-9 },
        contact: MultiContactConfig { max_contacts: 4, max_sweeps: 128, max_setup_terms: 16384 },
        sensitivity: SensitivityBudget { max_contacts: 4, max_setup_terms: derivative_budget,
            max_query_terms: derivative_budget, minimum_contact_margin_m: 1e-8 },
        max_cases: 4, max_variables: 8, max_bindings: 16, max_ports_per_case: 8,
    };
    let body = ModalAcousticTimeModel::try_free_mass(48000, 0.25, 0.0, 0.0,
        ModalAcousticTimeBudget::audible_reference()).unwrap();
    let spring = ModalConnection { left: map(2.0), right: map(0.0),
        stiffness_n_m: 64.0, damping_n_s_m: 0.0, rest_extension_m: 0.0 };
    let contact = ModalContact { left: map(2.0), right: map(0.0),
        law: Obstacle::new(vec![-1.0], 1, 1, vec![1.0/64.0], vec![1.0], 128.0, 1.0,
            "synthetic-observation-regression".into()).unwrap() };
    let contacts = vec![(contact, ModalContactConfig { max_iterations: 128, maximum_force_n: 10000.0,
        maximum_penetration_m: 0.1, force_absolute_tolerance_n: 1e-12, force_relative_tolerance: 1e-13 })];
    let cases = forces.into_iter().enumerate().map(|(i, force)| {
        let q = displacement(&[0.0; 6], force);
        DesignLoadCase { name: format!("case-{i}"), loads: vec![DesignLoad { attachment: map(2.0), force_n: force }],
            targets: vec![DisplacementTarget { attachment: map(2.0), target_m: q, scale_m: 0.01, weight: 4.0 },
                DisplacementTarget { attachment: map(-6.0), target_m: -3.0*q, scale_m: 0.02, weight: 0.0 }] }
    }).collect();
    let v = |name: &str, reference, scale, minimum, maximum, fields| DesignVariable {
        name: name.into(), reference, scale, minimum, maximum, fields,
    };
    let variables = vec![
        v("support", 64.0, 64.0, 32.0, 128.0, vec![DesignField::SpringStiffness(0), DesignField::SpringStiffness(1)]),
        v("rest", 0.0, 0.01, -0.01, 0.01, vec![DesignField::SpringRest(0), DesignField::SpringRest(1)]),
        v("contact", 128.0, 128.0, 64.0, 256.0, vec![DesignField::ContactStiffness(0)]),
        v("gap", 1.0/64.0, 1.0/64.0, 0.005, 0.04, vec![DesignField::ContactGap(0)]),
        v("weight", 1.0, 0.5, 0.25, 2.0, vec![DesignField::ContactWeight(0)]),
        v("second-load", forces[1], 0.5, 0.01, 1e10, vec![DesignField::ActuatorForce { case: 1, actuator: 0 }]),
    ];
    EquilibriumDesign::new(vec![body], vec![spring.clone(), spring], contacts, cases, variables, budget,
        &CancelGate::new()).unwrap()
}

#[test]
fn all_physical_rows_match_independent_derivatives_with_shared_fields_and_local_loads() {
    let p = build([0.5, 3.0], 16384);
    let gate = CancelGate::new();
    let mut work = DesignControl::new(4, 8);
    let mut queries = ObservationControl::new(4, 16).unwrap();
    let x = [0.1, 0.02, 0.2, 0.05, -0.2, 0.2];
    let result = p.evaluate_observations(&x, &mut work, &mut queries, &gate).unwrap();
    for (i, row) in result.rows().iter().enumerate() {
        assert_eq!((row.case, row.target), (i/2, i%2));
        let factor = if row.target == 0 { 1.0 } else { -3.0 };
        let value = |point: &[f64]| factor * displacement(point, if row.case == 0 { 0.5 } else { 3.0 + 0.5*point[5] });
        assert!((row.value_m - value(&x)).abs() < 1e-11);
        for j in 0..6 {
            let mut plus = x; let mut minus = x; plus[j] += 1e-5; minus[j] -= 1e-5;
            let fd = (value(&plus) - value(&minus)) / 2e-5;
            assert!((row.gradient_m[j] - fd).abs() < 1e-9, "row={i}, column={j}");
        }
        assert!(row.adjoint_relative_residual <= 1e-10);
        if row.target == 1 {
            assert_eq!(row.residual, 0.0);
            assert!(row.residual_gradient.iter().all(|g| *g == 0.0));
            assert!(row.gradient_m.iter().any(|g| g.abs() > 1e-4));
        }
        if row.case == 0 { assert_eq!(row.gradient_m[5], 0.0); }
    }
    let original = p.evaluate(&x, &mut work, &gate).unwrap();
    assert!((result.objective() - original.value).abs() < 1e-10);
    for (a, b) in result.gradient().iter().zip(&original.gradient) { assert!((a-b).abs() < 1e-8); }
    assert_eq!(queries.adjoints_attempted(), 4);
    assert_eq!(work.work(), DesignWork { evaluations: 2, case_solves: 4 });
}

#[test]
fn an_exact_fit_retains_information_and_gauss_newton_action() {
    let p = build([0.5, 3.0], 16384); let gate = CancelGate::new();
    let mut work = DesignControl::new(4, 8); let mut queries = ObservationControl::new(4, 8).unwrap();
    let result = p.evaluate_observations(&[0.0; 6], &mut work, &mut queries, &gate).unwrap();
    assert!(result.objective() < 1e-20);
    assert!(result.gradient().iter().all(|g| g.abs() < 1e-9));
    assert!(result.rows()[0].residual_gradient[0].abs() > 0.1);
    let direction = [0.2, -0.1, 0.3, 0.1, -0.2, 0.4]; let h = 1e-5;
    let plus = direction.map(|d| h*d); let minus = direction.map(|d| -h*d);
    let gp = p.evaluate(&plus, &mut work, &gate).unwrap().gradient;
    let gm = p.evaluate(&minus, &mut work, &gate).unwrap().gradient;
    let product = result.gauss_newton_product(&direction, &gate).unwrap();
    for j in 0..6 { assert!((product[j] - (gp[j]-gm[j])/(2.0*h)).abs() < 1e-6); }
    assert!(direction.iter().zip(&product).map(|(a,b)| a*b).sum::<f64>() >= 0.0);
    assert!(result.gauss_newton_product(&[0.0], &gate).is_err());
    assert!(result.gauss_newton_product(&[f64::NAN; 6], &gate).is_err());
    let cancel = CancelGate::new(); cancel.request();
    assert!(matches!(result.gauss_newton_product(&direction, &cancel), Err(DesignError::Cancelled)));
}

#[test]
fn complete_adjoint_allowance_is_preflighted_and_resume_never_refunds_queries() {
    let p = build([0.5, 3.0], 16384); let gate = CancelGate::new();
    let mut work = DesignControl::new(4, 8); let mut queries = ObservationControl::new(4, 3).unwrap();
    assert!(matches!(p.evaluate_observations(&[0.0; 6], &mut work, &mut queries, &gate), Err(DesignError::Budget { .. })));
    assert_eq!(work.work(), DesignWork::default()); assert_eq!(queries.adjoints_attempted(), 0);
    queries.extend(8).unwrap();
    let first = p.evaluate_observations(&[0.0; 6], &mut work, &mut queries, &gate).unwrap();
    let cancelled = CancelGate::new(); cancelled.request();
    assert!(matches!(p.evaluate_observations(&[0.0; 6], &mut work, &mut queries, &cancelled), Err(DesignError::Cancelled)));
    assert_eq!(queries.adjoints_attempted(), 4);
    assert_eq!(first, p.evaluate_observations(&[0.0; 6], &mut work, &mut queries, &gate).unwrap());
    assert_eq!(queries.adjoints_attempted(), 8); assert!(queries.extend(7).is_err());
    assert!(p.evaluate_observations(&[0.0; 6], &mut work, &mut queries, &gate).is_err());
    assert_eq!(queries.adjoints_attempted(), 8);
    assert!(ObservationControl::new(0, 1).is_err()); assert!(ObservationControl::new(1025, 1).is_err());
}

#[test]
fn contact_switches_and_late_physical_failures_do_not_publish_partial_jacobians() {
    let gate = CancelGate::new();
    for force in [2.0, 1e9] {
        let p = build([0.5, force], 16384);
        let mut work = DesignControl::new(2, 4); let mut queries = ObservationControl::new(4, 8).unwrap();
        assert!(matches!(p.evaluate_observations(&[0.0; 6], &mut work, &mut queries, &gate), Err(DesignError::Case { case: 1, .. })));
        assert_eq!(work.work(), DesignWork { evaluations: 1, case_solves: 2 });
        assert_eq!(queries.adjoints_attempted(), 2);
        if force == 2.0 { assert!(p.evaluate_forward(&[0.0; 6], &mut work, &gate).is_ok()); }
    }
    let p = build([0.5, 3.0], 0); let mut work = DesignControl::new(1, 2);
    let mut queries = ObservationControl::new(4, 4).unwrap();
    assert!(p.evaluate_observations(&[0.0; 6], &mut work, &mut queries, &gate).is_err());
    assert_eq!(queries.adjoints_attempted(), 0);
}

#[test]
fn calibration_rows_do_not_silently_assess_physical_constraints() {
    use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::constraints::{
        ConstraintSense, ResponseConstraint, ResponseQuantity,
    };
    let gate = CancelGate::new();
    let p = build([0.5, 3.0], 16384).with_constraints(vec![ResponseConstraint {
        name: "deliberately-violated".into(), case: 0, quantity: ResponseQuantity::Displacement(map(2.0)),
        sense: ConstraintSense::AtMost, bound: -1.0, scale: 1.0,
    }], 1, &gate).unwrap();
    let mut work = DesignControl::new(2, 4); let mut queries = ObservationControl::new(4, 4).unwrap();
    let rows = p.evaluate_observations(&[0.0; 6], &mut work, &mut queries, &gate).unwrap();
    assert_eq!(rows.unassessed_response_constraints(), 1);
    assert!(p.evaluate(&[0.0; 6], &mut work, &gate).unwrap().constraints[0].residual > 0.0);
}
