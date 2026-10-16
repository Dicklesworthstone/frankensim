//! Physical response Jacobians versus independently solved physical coordinates.
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{ModalAttachment, ModalConnection, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::SensitivityBudget;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::DisplacementTarget;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::*;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::constraints::*;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn map(component: usize) -> ModalAttachment {
    ModalAttachment { component, shapes: vec![if component == 0 { 5.0 } else { 1.0 }] }
}
fn problem(shared: bool) -> EquilibriumDesign {
    let models = vec![
        ModalAcousticTimeModel::try_free_mass(48000, 0.04, 0.0, 0.0, ModalAcousticTimeBudget::audible_reference()).unwrap(),
        ModalAcousticTimeModel::try_new(48000, vec![ModalAcousticMode { angular_frequency_rad_s: 100.0,
            damping_ratio: 0.02, pressure_per_modal_velocity: C64::ZERO }], ModalAcousticTimeBudget::audible_reference()).unwrap(),
    ];
    let spring = ModalConnection { left: map(0), right: ModalAttachment { component: 0, shapes: vec![0.0] },
        stiffness_n_m: 400.0, damping_n_s_m: 0.0, rest_extension_m: 0.0 };
    let springs = if shared { vec![spring.clone(), spring] } else { vec![spring] };
    let contacts = vec![(ModalContact { left: map(0), right: map(1), law: Obstacle::new(vec![-1.0], 1, 1,
        vec![0.0001], vec![1.0], 1e8, 2.0, "authored-constraint-test".into()).unwrap() },
        ModalContactConfig { max_iterations: 128, maximum_force_n: 10000.0, maximum_penetration_m: 0.1,
            force_absolute_tolerance_n: 1e-12, force_relative_tolerance: 1e-13 })];
    let cases = [0.8, 1.6].into_iter().enumerate().map(|(i, force)| DesignLoadCase {
        name: format!("case-{i}"), loads: vec![DesignLoad { attachment: map(0), force_n: force }],
        targets: vec![DisplacementTarget { attachment: map(1), target_m: 0.0002, scale_m: 0.0001, weight: 1.0 }],
    }).collect();
    let fields = [DesignField::SpringStiffness(0), DesignField::SpringRest(0), DesignField::ContactStiffness(0),
        DesignField::ContactGap(0), DesignField::ContactWeight(0), DesignField::ActuatorForce { case: 1, actuator: 0 }];
    let values = [(400.0, 400.0, 100.0, 1200.0), (0.0, 0.0001, -0.001, 0.001),
        (1e8, 1e8, 0.0, 4e8), (0.0001, 0.0002, -0.001, 0.01), (1.0, 0.5, 0.0, 3.0), (1.6, 1.0, 0.0, 4.0)];
    let mut variables: Vec<_> = fields.into_iter().zip(values).enumerate().map(|(i, (field, (reference, scale, minimum, maximum)))|
        DesignVariable { name: format!("p-{i}"), reference, scale, minimum, maximum, fields: vec![field] }).collect();
    if shared { variables[0].fields.push(DesignField::SpringStiffness(1)); variables[1].fields.push(DesignField::SpringRest(1)); }
    EquilibriumDesign::new(models, springs, contacts, cases, variables, DesignBudget {
        coupling: ModalCouplingConfig { max_modes: 8, max_connections: 4, max_setup_terms: 16384,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 1000.0, maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 10000.0, solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-11, energy_relative_tolerance: 1e-9 },
        contact: MultiContactConfig { max_contacts: 4, max_sweeps: 128, max_setup_terms: 16384 },
        sensitivity: SensitivityBudget { max_contacts: 4, max_setup_terms: 16384, max_query_terms: 16384, minimum_contact_margin_m: 1e-9 },
        max_cases: 4, max_variables: 8, max_bindings: 16, max_ports_per_case: 8,
    }, &CancelGate::new()).unwrap()
}
fn requirements() -> Vec<ResponseConstraint> {
    // Deliberately interleave cases; returned rows must retain declaration order.
    vec![
        ResponseConstraint { name: "normal-cap".into(), case: 1, quantity: ResponseQuantity::ContactForce(0), sense: ConstraintSense::AtMost, bound: 0.4, scale: 1.0 },
        ResponseConstraint { name: "motion".into(), case: 0, quantity: ResponseQuantity::Displacement(map(1)), sense: ConstraintSense::AtLeast, bound: 0.0001, scale: 0.0002 },
        ResponseConstraint { name: "support".into(), case: 1, quantity: ResponseQuantity::SpringForce(0), sense: ConstraintSense::Equal, bound: -0.2, scale: 2.0 },
        ResponseConstraint { name: "penetration".into(), case: 1, quantity: ResponseQuantity::ContactPenetration(0), sense: ConstraintSense::AtMost, bound: 0.0001, scale: 0.0002 },
    ]
}
// Closed-form physical-coordinate solution for the quadratic contact, no modal
// matrices, production preload, reaction differentiation or adjoint calls.
fn physical(x: &[f64; 6], case: usize, shared: bool) -> [f64; 4] {
    let support = (400.0 + 400.0*x[0]) * if shared { 2.0 } else { 1.0 };
    let rest = 0.0001*x[1]; let k = 1e8 + 1e8*x[2]; let gap = 0.0001 + 0.0002*x[3];
    let weight = 1.0 + 0.5*x[4]; let load = if case == 0 { 0.8 } else { 1.6 + x[5] };
    let excess = load/support + rest - gap;
    let penetration = if excess <= 0.0 { 0.0 } else {
        2.0*excess / (1.0 + (1.0 + 4.0*weight*k*(1.0/support + 0.0001)*excess).sqrt())
    };
    let reaction = weight*k*penetration*penetration;
    [reaction, reaction/10000.0, (reaction-load)/if shared { 2.0 } else { 1.0 }, penetration]
}
fn oracle(x: &[f64; 6], shared: bool) -> ([f64; 4], [f64; 4]) {
    let a = physical(x, 0, shared); let b = physical(x, 1, shared);
    ([b[0], a[1], b[2], b[3]], [b[0]-0.4, (0.0001-a[1])/0.0002, (b[2]+0.2)/2.0, (b[3]-0.0001)/0.0002])
}
fn evaluate(p: &EquilibriumDesign, x: &[f64;6]) -> DesignEvaluation {
    p.evaluate(x, &mut DesignControl::new(1, 2), &CancelGate::new()).unwrap()
}

#[test]
fn all_response_rows_include_direct_and_implicit_parameter_derivatives() {
    for shared in [false, true] {
        let p = problem(shared).with_constraints(requirements(), 4, &CancelGate::new()).unwrap();
        for x in [[0.0;6], [0.2,0.3,0.4,-0.2,0.2,0.1]] {
            let result = evaluate(&p, &x); let (values, residuals) = oracle(&x, shared);
            for i in 0..4 {
                let row = &result.constraints[i];
                assert!((row.value-values[i]).abs() < 1e-10);
                assert!((row.residual-residuals[i]).abs() < 1e-9);
                assert!(row.adjoint_relative_residual <= 1e-10);
                for j in 0..6 {
                    let h = 1e-5; let mut plus = x; let mut minus = x; plus[j] += h; minus[j] -= h;
                    let fd = (oracle(&plus, shared).1[i]-oracle(&minus, shared).1[i])/(2.0*h);
                    assert!((row.gradient[j]-fd).abs() < 2e-6*fd.abs().max(1e-5), "row {i} field {j}: {} != {fd}", row.gradient[j]);
                }
            }
            assert_eq!(result.constraints[1].gradient[5], 0.0, "case-local force must not affect another experiment");
        }
    }
}

#[test]
fn violated_constraints_are_returned_without_penalizing_the_original_objective() {
    let unconstrained = evaluate(&problem(false), &[0.0;6]);
    let p = problem(false).with_constraints(requirements(), 4, &CancelGate::new()).unwrap();
    let constrained = evaluate(&p, &[0.0;6]);
    assert_eq!(unconstrained.value, constrained.value);
    assert_eq!(unconstrained.gradient, constrained.gradient);
    assert_eq!(unconstrained.cases, constrained.cases);
    assert!(constrained.constraints[0].residual > 0.0, "infeasible but mechanically valid response must survive");
}

#[test]
fn separated_and_zero_coefficient_contacts_keep_the_correct_local_derivatives() {
    let p = problem(false).with_constraints(requirements(), 4, &CancelGate::new()).unwrap();
    let mut separated = [0.0;6]; separated[3] = 20.0;
    let r = evaluate(&p, &separated);
    for i in [0,3] { assert_eq!(r.constraints[i].value, 0.0); assert!(r.constraints[i].gradient.iter().all(|g| *g == 0.0)); }
    let mut disabled = [0.0;6]; disabled[2] = -1.0;
    let r = evaluate(&p, &disabled);
    assert_eq!(r.constraints[0].value, 0.0);
    assert!(r.constraints[0].gradient[2] > 0.0, "K=0 has a nonzero coefficient derivative while penetrating");
    let pen = physical(&disabled, 1, false)[3];
    assert!((r.constraints[0].gradient[2]-1e8*pen*pen).abs() < 1e-8);
}

#[test]
fn constraint_admission_cancellation_and_activity_refusals_do_not_leak_partial_rows() {
    for mutate in 0..5 {
        let mut rows = requirements();
        match mutate {
            0 => rows[1].name = rows[0].name.clone(),
            1 => rows[0].case = 99,
            2 => rows[0].quantity = ResponseQuantity::ContactForce(99),
            3 => rows[0].scale = 0.0,
            _ => rows[1].quantity = ResponseQuantity::Displacement(ModalAttachment { component: 1, shapes: vec![1.0,2.0] }),
        }
        assert!(problem(false).with_constraints(rows, 4, &CancelGate::new()).is_err());
    }
    assert!(problem(false).with_constraints(requirements(), 3, &CancelGate::new()).is_err());
    let p = problem(false).with_constraints(requirements(), 4, &CancelGate::new()).unwrap();
    let mut control = DesignControl::new(4,8); let cancelled = CancelGate::new(); cancelled.request();
    assert!(matches!(p.evaluate(&[0.0;6], &mut control, &cancelled), Err(DesignError::Cancelled)));
    assert_eq!(control.work(), DesignWork::default());
    let mut touching = [0.0;6]; touching[3] = (0.8/400.0-0.0001)/0.0002;
    assert!(p.evaluate(&touching, &mut control, &CancelGate::new()).is_err());
    assert_eq!(control.work().case_solves,1);
    let a = p.evaluate(&[0.0;6], &mut control, &CancelGate::new()).unwrap();
    let b = p.evaluate(&[0.0;6], &mut control, &CancelGate::new()).unwrap();
    assert_eq!(a,b); assert_eq!(control.work(), DesignWork { evaluations:3, case_solves:5 });
}
