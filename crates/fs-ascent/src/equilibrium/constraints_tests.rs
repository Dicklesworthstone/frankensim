use super::*;
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{ModalAttachment, ModalConnection, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::SensitivityBudget;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::DisplacementTarget;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::{
    DesignBudget, DesignField, DesignLoad, DesignLoadCase, DesignVariable,
    constraints::{ResponseQuantity, ResponseConstraint},
};
use fs_dcontact::Obstacle;
use fs_math::c64::C64;

fn budget() -> DesignBudget {
    DesignBudget {
        coupling: ModalCouplingConfig { max_modes: 8, max_connections: 4, max_setup_terms: 16384,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 1000.0, maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 10000.0, solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-11, energy_relative_tolerance: 1e-9 },
        contact: MultiContactConfig { max_contacts: 4, max_sweeps: 128, max_setup_terms: 16384 },
        sensitivity: SensitivityBudget { max_contacts: 4, max_setup_terms: 16384, max_query_terms: 16384, minimum_contact_margin_m: 1e-9 },
        max_cases: 4, max_variables: 4, max_bindings: 8, max_ports_per_case: 8,
    }
}
fn variable() -> DesignVariable {
    DesignVariable { name: "load-N".into(), reference: 1.5, scale: 1.0, minimum: 0.1, maximum: 4.0,
        fields: vec![DesignField::ActuatorForce { case: 0, actuator: 0 }] }
}
fn elastic(omega: f64) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(48000, vec![ModalAcousticMode { angular_frequency_rad_s: omega,
        damping_ratio: 0.02, pressure_per_modal_velocity: C64::ZERO }], ModalAcousticTimeBudget::audible_reference()).unwrap()
}
fn displacement_problem(sense: ConstraintSense) -> EquilibriumDesign {
    let port = ModalAttachment { component: 0, shapes: vec![1.0] };
    let p = EquilibriumDesign::new(vec![elastic(10.0)], vec![], vec![], vec![DesignLoadCase {
        name: "experiment".into(), loads: vec![DesignLoad { attachment: port.clone(), force_n: 1.5 }],
        targets: vec![DisplacementTarget { attachment: port.clone(), target_m: if sense == ConstraintSense::AtLeast { 0.0 } else { 0.02 },
            scale_m: 0.01, weight: 1.0 }],
    }], vec![variable()], budget(), &CancelGate::new()).unwrap();
    p.with_constraints(vec![ResponseConstraint { name: "physical-travel".into(), case: 0, quantity: ResponseQuantity::Displacement(port),
        sense, bound: 0.01, scale: 0.01 }], 1, &CancelGate::new()).unwrap()
}

#[test]
fn physical_upper_lower_and_equality_constraints_have_real_kkt_multipliers() {
    for sense in [ConstraintSense::AtMost, ConstraintSense::AtLeast, ConstraintSense::Equal] {
        let p = displacement_problem(sense); let mut work = DesignControl::new(64,64); let gate = CancelGate::new();
        let mut study = EquilibriumStudy::new(&p, &[0.0], &mut work, 4, &gate).unwrap();
        if sense == ConstraintSense::AtMost { assert!(study.accepted().constraints[0].residual > 0.0); }
        let result = study.run(1e-8,32,64,&gate).unwrap();
        assert_eq!(result.stop, SqpStop::Converged);
        assert!(result.solution.kkt.within_tolerance(1e-8));
        assert!((study.accepted().physical_parameters[0]-1.0).abs() < 1e-7);
        assert!(result.solution.nu[..2].iter().all(|v| v.abs() < 1e-8), "solution is interior to the parameter box");
        let multiplier = if sense == ConstraintSense::Equal { result.solution.lambda[0] } else { result.solution.nu[2] };
        assert!((multiplier-1.0).abs() < 1e-7);
    }
}

#[test]
fn nonlinear_contact_force_limit_is_solved_from_an_infeasible_mechanical_state() {
    let mass = ModalAcousticTimeModel::try_free_mass(48000,0.04,0.0,0.0,ModalAcousticTimeBudget::audible_reference()).unwrap();
    let left = ModalAttachment { component: 0, shapes: vec![5.0] };
    let right = ModalAttachment { component: 1, shapes: vec![1.0] };
    let p = EquilibriumDesign::new(vec![mass,elastic(100.0)],vec![ModalConnection {
        left: left.clone(), right: ModalAttachment { component: 0, shapes: vec![0.0] },
        stiffness_n_m:400.0,damping_n_s_m:0.0,rest_extension_m:0.0,
    }],vec![(ModalContact { left:left.clone(),right:right.clone(),law:Obstacle::new(vec![-1.0],1,1,
        vec![0.0001],vec![1.0],1e8,2.0,"authored-sqp-constraint".into()).unwrap() },
        ModalContactConfig { max_iterations:128,maximum_force_n:10000.0,maximum_penetration_m:0.1,
            force_absolute_tolerance_n:1e-12,force_relative_tolerance:1e-13 })],vec![DesignLoadCase {
        name:"contact-load".into(),loads:vec![DesignLoad {attachment:left,force_n:1.5}],
        targets:vec![DisplacementTarget {attachment:right,target_m:0.0002,scale_m:0.0001,weight:1.0}],
    }],vec![variable()],budget(),&CancelGate::new()).unwrap().with_constraints(vec![ResponseConstraint {
        name:"normal-force".into(),case:0,quantity:ResponseQuantity::ContactForce(0),sense:ConstraintSense::AtMost,bound:0.8,scale:1.0,
    }],1,&CancelGate::new()).unwrap();
    let mut work=DesignControl::new(128,128);let gate=CancelGate::new();
    let mut study=EquilibriumStudy::new(&p,&[0.0],&mut work,4,&gate).unwrap();
    assert!(study.accepted().constraints[0].residual>0.5);
    let report=study.run(1e-8,64,128,&gate).unwrap();
    assert_eq!(report.stop,SqpStop::Converged);
    let expected=0.8+400.0*(0.0001+0.8/10000.0+(0.8_f64/1e8).sqrt());
    assert!((study.accepted().physical_parameters[0]-expected).abs()<1e-7);
    assert!((study.accepted().constraints[0].value-0.8).abs()<1e-8);
    assert!((report.solution.nu[2]-1.2).abs()<1e-7);
}

#[test]
fn physical_rows_count_toward_kkt_admission_and_survive_budgeted_cancellation() {
    let p=displacement_problem(ConstraintSense::AtMost);let gate=CancelGate::new();let mut work=DesignControl::new(32,32);
    assert!(EquilibriumStudy::new(&p,&[0.0],&mut work,3,&gate).is_err());
    assert_eq!(work.work().evaluations,0);
    let mut study=EquilibriumStudy::new(&p,&[0.0],&mut work,4,&gate).unwrap();
    let initial=study.accepted().clone();let before=study.work();
    assert_eq!(study.run(1e-8,32,1,&gate).unwrap().stop,SqpStop::EvaluationLimit);
    assert_eq!(study.accepted(),&initial);
    let cancelled=CancelGate::new();cancelled.request();assert!(study.run(1e-8,1,32,&cancelled).is_err());
    assert_eq!(study.work(),before);assert_eq!(study.accepted(),&initial);
    assert_eq!(study.run(1e-8,32,32,&gate).unwrap().stop,SqpStop::Converged);
    assert!(study.accepted().constraints[0].residual.abs()<1e-8);
}
