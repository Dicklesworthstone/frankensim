use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{ModalAttachment, ModalConnection, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::SensitivityBudget;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::DisplacementTarget;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::*;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

const LOADS: [f64; 3] = [0.8, 1.6, 2.5];
const SCALE: f64 = 0.0002;
fn map(component: usize) -> ModalAttachment {
    ModalAttachment { component, shapes: vec![if component == 0 { 5.0 } else { 1.0 }] }
}
fn budget() -> DesignBudget {
    DesignBudget {
        coupling: ModalCouplingConfig { max_modes:8, max_connections:4, max_setup_terms:16384,
            nyquist_guard_fraction:0.9, maximum_total_energy_j:1000.0, maximum_abs_pressure_pa:1e6,
            maximum_abs_connection_force_n:10000.0, solve_relative_tolerance:1e-10,
            energy_absolute_tolerance_j:1e-11, energy_relative_tolerance:1e-9 },
        contact: MultiContactConfig { max_contacts:4, max_sweeps:128, max_setup_terms:16384 },
        sensitivity: SensitivityBudget { max_contacts:4, max_setup_terms:16384,
            max_query_terms:16384, minimum_contact_margin_m:1e-9 },
        max_cases:4, max_variables:8, max_bindings:16, max_ports_per_case:8,
    }
}
fn variables() -> Vec<DesignVariable> {
    vec![
        DesignVariable { name:"support".into(), reference:400.0, scale:400.0, minimum:200.0, maximum:1000.0,
            fields:vec![DesignField::SpringStiffness(0)] },
        DesignVariable { name:"contact".into(), reference:1e8, scale:1e8, minimum:2e7, maximum:4e8,
            fields:vec![DesignField::ContactStiffness(0)] },
        DesignVariable { name:"gap".into(), reference:0.0001, scale:0.0002, minimum:0.00001, maximum:0.0005,
            fields:vec![DesignField::ContactGap(0)] },
    ]
}
// Independent physical-coordinate quadratic equilibrium, not an iteration,
// modal inverse, discrete gradient or tangent implementation from production.
fn observations(support: f64, contact: f64, gap: f64, force: f64) -> [f64; 2] {
    let excess = force/support-gap;
    if excess <= 0.0 { return [force/support, 0.0]; }
    let compliance = 1.0/support+1.0/10000.0;
    let penetration = 2.0*excess/(1.0+(1.0+4.0*contact*compliance*excess).sqrt());
    let reaction = contact*penetration*penetration;
    [(force-reaction)/support, reaction/10000.0]
}
fn cases() -> Vec<DesignLoadCase> {
    LOADS.into_iter().enumerate().map(|(i, force)| {
        let target = observations(600.0,1.8e8,0.0002,force);
        DesignLoadCase { name:format!("load-{i}"), loads:vec![DesignLoad { attachment:map(0), force_n:force }],
            targets:(0..2).map(|component| DisplacementTarget { attachment:map(component),
                target_m:target[component], scale_m:SCALE, weight:1.0 }).collect() }
    }).collect()
}
fn build(vars: Vec<DesignVariable>, cases: Vec<DesignLoadCase>, shared: bool) -> Result<EquilibriumDesign, DesignError> {
    let models=vec![
        ModalAcousticTimeModel::try_free_mass(48000,0.04,0.0,0.0,ModalAcousticTimeBudget::audible_reference()).unwrap(),
        ModalAcousticTimeModel::try_new(48000,vec![ModalAcousticMode { angular_frequency_rad_s:100.0,
            damping_ratio:0.02,pressure_per_modal_velocity:C64::ZERO }],ModalAcousticTimeBudget::audible_reference()).unwrap(),
    ];
    let spring=ModalConnection { left:map(0),right:ModalAttachment { component:0,shapes:vec![0.0] },
        stiffness_n_m:400.0,damping_n_s_m:0.0,rest_extension_m:0.0 };
    let springs=if shared {vec![spring.clone(),spring]}else{vec![spring]};
    let contacts=vec![(ModalContact { left:map(0),right:map(1),law:Obstacle::new(vec![-1.0],1,1,
        vec![0.0001],vec![1.0],1e8,2.0,"synthetic-design-regression".into()).unwrap() },
        ModalContactConfig { max_iterations:128,maximum_force_n:10000.0,maximum_penetration_m:0.1,
            force_absolute_tolerance_n:1e-12,force_relative_tolerance:1e-13 })];
    EquilibriumDesign::new(models,springs,contacts,cases,vars,budget(),&CancelGate::new())
}
fn analytic(x: &[f64;3], shared: bool, forces: [f64;3]) -> f64 {
    let support=(400.0+400.0*x[0])*if shared {2.0}else{1.0};
    let k=1e8+1e8*x[1];let gap=0.0001+0.0002*x[2];
    LOADS.into_iter().zip(forces).map(|(original,force)| {
        let expected=observations(600.0,1.8e8,0.0002,original);
        let actual=observations(support,k,gap,force);
        (0..2).map(|i|0.5*((actual[i]-expected[i])/SCALE).powi(2)).sum::<f64>()
    }).sum()
}

#[test]
fn complete_multiload_value_and_all_scaled_gradients_match_independent_equilibria() {
    let problem=build(variables(),cases(),false).unwrap();let gate=CancelGate::new();
    let mut control=DesignControl::new(16,48);
    for x in [[0.0,0.0,0.0],[0.25,0.3,0.1],[-0.2,0.8,0.4]] {
        let result=problem.evaluate(&x,&mut control,&gate).unwrap();
        assert!((result.value-analytic(&x,false,LOADS)).abs()<1e-9);
        for j in 0..3 {
            let h=1e-4;let mut plus=x;let mut minus=x;plus[j]+=h;minus[j]-=h;
            let fd=(analytic(&plus,false,LOADS)-analytic(&minus,false,LOADS))/(2.0*h);
            assert!((fd-result.gradient[j]).abs()<1e-6*fd.abs().max(1e-3),"{j}: {fd} vs {:?}",result.gradient);
        }
        for (case,force) in result.cases.iter().zip(LOADS) {
            let expected=observations(400.0+400.0*x[0],1e8+1e8*x[1],0.0001+0.0002*x[2],force);
            for (a,b) in case.observations_m.iter().zip(expected) {assert!((a-b).abs()<1e-12);}
            assert_eq!(case.equilibrium.active_contacts,1);
            assert!(case.adjoint_relative_residual<=1e-10);
        }
    }
    assert_eq!(control.work(),DesignWork {evaluations:3,case_solves:9});
}

#[test]
fn shared_fields_sum_both_pullbacks_instead_of_overwriting_a_gradient() {
    let mut vars=variables();vars[0].fields.push(DesignField::SpringStiffness(1));
    let problem=build(vars,cases(),true).unwrap();let x=[-0.1,0.4,0.2];
    let result=problem.evaluate(&x,&mut DesignControl::new(1,3),&CancelGate::new()).unwrap();
    assert!((result.value-analytic(&x,true,LOADS)).abs()<1e-9);
    let h=1e-4;let mut plus=x;let mut minus=x;plus[0]+=h;minus[0]-=h;
    let expected=(analytic(&plus,true,LOADS)-analytic(&minus,true,LOADS))/(2.0*h);
    assert!((result.gradient[0]-expected).abs()<1e-6*expected.abs());
}

#[test]
fn a_case_specific_load_gradient_is_per_physical_newton_and_local_to_that_case() {
    let mut vars=variables();vars.push(DesignVariable {name:"second-load".into(),reference:1.6,scale:0.3,
        minimum:0.5,maximum:3.0,fields:vec![DesignField::ActuatorForce {case:1,actuator:0}]});
    let problem=build(vars,cases(),false).unwrap();let x=[0.1,0.2,0.2,0.3];
    let result=problem.evaluate(&x,&mut DesignControl::new(1,3),&CancelGate::new()).unwrap();
    let f=1.6+0.3*x[3];let h=1e-5;
    let expected=(analytic(&[x[0],x[1],x[2]],false,[0.8,f+0.3*h,2.5])
        -analytic(&[x[0],x[1],x[2]],false,[0.8,f-0.3*h,2.5]))/(2.0*h);
    assert!((result.gradient[3]-expected).abs()<1e-7*expected.abs());
}

#[test]
fn cancelled_budgeted_and_out_of_domain_candidates_preserve_work_and_replay() {
    let problem=build(variables(),cases(),false).unwrap();let gate=CancelGate::new();
    let mut control=DesignControl::new(2,2);
    assert!(matches!(problem.evaluate(&[0.0;3],&mut control,&gate),Err(DesignError::Budget {..})));
    assert_eq!(control.work(),DesignWork {evaluations:1,case_solves:0});
    control.extend(8,12).unwrap();
    assert!(matches!(problem.evaluate(&[-1.0,0.0,0.0],&mut control,&gate),Err(DesignError::OutsideBounds {..})));
    assert_eq!(control.work(),DesignWork {evaluations:2,case_solves:0});
    let cancelled=CancelGate::new();cancelled.request();
    assert!(matches!(problem.evaluate(&[0.0;3],&mut control,&cancelled),Err(DesignError::Cancelled)));
    assert_eq!(control.work(),DesignWork {evaluations:2,case_solves:0});
    let a=problem.evaluate(&[0.0;3],&mut control,&gate).unwrap();
    let b=problem.evaluate(&[0.0;3],&mut control,&gate).unwrap();
    assert_eq!(a,b);assert_eq!(control.work(),DesignWork {evaluations:4,case_solves:6});
    assert!(control.extend(1,1).is_err());
}

#[test]
fn a_late_case_failure_returns_no_partial_objective_and_consumes_attempted_work() {
    let mut experiments=cases();
    experiments[1].loads[0].force_n=1e9; // finite, but physically inadmissible
    let problem=build(variables(),experiments,false).unwrap();let mut control=DesignControl::new(4,12);
    for evaluation in 1..=2 {
        assert!(matches!(problem.evaluate(&[0.0;3],&mut control,&CancelGate::new()),Err(DesignError::Case {case:1,..})));
        assert_eq!(control.work(),DesignWork {evaluations:evaluation,case_solves:2*evaluation});
    }
}

#[test]
fn duplicate_unknown_and_dimensionally_incompatible_bindings_are_refused() {
    let mut duplicate=variables();duplicate[1].fields.push(DesignField::SpringStiffness(0));
    assert!(build(duplicate,cases(),false).is_err());
    let mut mixed=variables();mixed[0].fields.push(DesignField::ContactGap(0));mixed[2].fields=vec![DesignField::SpringRest(0)];
    assert!(build(mixed,cases(),false).is_err());
    let mut unknown=variables();unknown[1].fields=vec![DesignField::ContactStiffness(8)];
    assert!(build(unknown,cases(),false).is_err());
    let mut bad_scale=variables();bad_scale[0].scale=0.0;assert!(build(bad_scale,cases(),false).is_err());
    let mut bad_name=cases();bad_name[1].name=bad_name[0].name.clone();assert!(build(variables(),bad_name,false).is_err());
    let p=build(variables(),cases(),false).unwrap();
    assert!(p.physical_parameters(&[f64::NAN;3]).is_err());assert!(p.physical_parameters(&[]).is_err());
}
