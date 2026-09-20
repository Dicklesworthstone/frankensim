use super::*;
use fs_couple::render::schedule::force::file::design::EquilibriumDesignFile;

const MODEL: &str = "frankensim-modal-performance-v2\nsample_rate_hz 48000\nsamples 1\nfull_scale_pa 1\nlimits 0.9 10 1000 1000 1000000\ncompile_limits 0 1\nvoices 1\nvoice retain-state 1 1\nmode 2 0.1 0 0 0 0\nport 0 1\ncoupling_limits 0 16384 0.9 1000 1000000 10000 1e-10 1e-11 1e-9\nconnections 0\nevents 0\n";
const DESIGN: &str = "frankensim-equilibrium-design-v1\npreload_limits 0 1 16384\nsensitivity_limits 0 16384 16384 0\ndesign_limits 1 1 1 2\ncases 1\ncase response 1 1\nload 0 0 4\ntarget 0 0 1 1 1\nvariables 1\nvariable load-N 4 4 1 8 1\nbind actuator-force 0 0\n";
fn load(model: &str, design: &str) -> EquilibriumDesignFile {
    EquilibriumDesignFile::from_bytes(model.as_bytes(),design.as_bytes(),&CancelGate::new()).unwrap()
}
fn ensemble<'a>(p: &'a EquilibriumDesign, offsets: &[f64], cap: usize) -> ScenarioProblem<'a> {
    ScenarioProblem::new(p,offsets.iter().enumerate().map(|(i,d)| EquilibriumScenario {
        name:format!("realization-{i}"),physical_offsets:vec![*d],
    }).collect(),32,cap,&CancelGate::new()).unwrap()
}

#[test]
fn empirical_tail_optimizes_fractional_mass_instead_of_rounding_to_the_worst_scenario() {
    // u=F/4; losses are (F-4)^2/32 three times and F^2/32 once.
    // With tail mass m in [2,3], optimum F=4-4/m and risk=(m-1)/(2*m*m).
    let loaded=load(MODEL,DESIGN);let gate=CancelGate::new();
    for alpha in [0.25,0.4] {
        let e=ensemble(loaded.problem(),&[0.0,0.0,0.0,4.0],32).with_cvar(alpha).unwrap();
        let mut work=DesignControl::new(1600,1600);
        let mut study=ScenarioEquilibriumStudy::new(e,&[0.0],&mut work,&gate).unwrap();
        let report=study.run(1e-9,128,400,&gate).unwrap();
        assert!(report.solution.converged,"{:?}: {:?}",report.stop,report.solution.kkt);
        let audit=study.recheck(&gate).unwrap();let m=4.0*(1.0-alpha);
        let expected_force=4.0-4.0/m;let expected_risk=(m-1.0)/(2.0*m*m);
        assert!((audit.nominal_parameters[0]-expected_force).abs()<1e-6);
        let score=study.ensemble().cvar_upper_bound(&audit,report.solution.x[1]).unwrap();
        assert!((score-expected_risk).abs()<1e-8);
        assert!((score-report.solution.f).abs()<1e-8);
        assert!(score<0.125 && audit.worst_objective>0.125);
        assert!((report.solution.nu[2]+report.solution.nu[4]+report.solution.nu[6]+report.solution.nu[8]-1.0).abs()<1e-8);
    }
}

#[test]
fn physical_constraints_remain_binding_even_outside_the_selected_loss_tail() {
    let model=include_str!(concat!(env!("CARGO_MANIFEST_DIR"),"/../fs-couple/examples/equilibrium-design.model"));
    let design=include_str!(concat!(env!("CARGO_MANIFEST_DIR"),"/../fs-couple/examples/equilibrium-response-limits.fit"));
    let loaded=load(model,design);let gate=CancelGate::new();
    let e=ensemble(loaded.problem(),&[-0.1,0.0,0.0,0.1],64).with_cvar(0.4).unwrap();
    let mut work=DesignControl::new(1600,1600);
    let mut study=ScenarioEquilibriumStudy::new(e,&[0.0],&mut work,&gate).unwrap();
    let report=study.run(1e-8,128,400,&gate).unwrap();
    assert!(report.solution.converged,"{:?}: {:?}",report.stop,report.solution.kkt);
    let audit=study.recheck(&gate).unwrap();
    let force=0.8+400.0*(0.0001+(0.8_f64/1e8).sqrt()+0.8/10000.0)-0.1;
    assert!((audit.nominal_parameters[0]-force).abs()<1e-6);
    for scenario in &audit.scenarios { assert!(scenario.constraints.iter().all(|r|r.residual<=1e-7)); }
    // Four scenarios; each block is loss, nonnegative slack, four physical rows.
    let high=2+3*6;
    assert!(report.solution.nu[high].abs()<1e-7);
    assert!(report.solution.nu[high+2]>0.1);
    assert!(audit.scenarios[3].constraints[0].residual.abs()<1e-7);
}

#[test]
fn all_tail_and_physical_jacobian_rows_match_differenced_complete_samples() {
    let design=format!("{DESIGN}constraint_limits 3\nconstraints 3\nconstraint cap 0 displacement 0 0 at-most 1.2 2\nconstraint floor 0 displacement 0 0 at-least 0.1 3\nconstraint exact 0 displacement 0 0 equal 1 4\n");
    let loaded=load(MODEL,&design);let gate=CancelGate::new();
    let e=ensemble(loaded.problem(),&[0.0,0.0,0.0,4.0],64).with_cvar(0.4).unwrap();
    let mut work=DesignControl::new(100,100);
    let point=[-0.1,0.1,0.05,0.02,0.04,0.1];let width=point.len();
    let make=|p:&[f64],w:&mut DesignControl| e.sample(p,&e.evaluate(&p[..1],w,&gate).unwrap());
    let base=make(&point,&mut work);
    assert_eq!(base.ci.len(),18);assert_eq!(base.ce.len(),4);
    for j in 0..width {
        let h=1e-6;let mut plus=point;let mut minus=point;plus[j]+=h;minus[j]-=h;
        let a=make(&plus,&mut work);let b=make(&minus,&mut work);
        assert!(((a.f-b.f)/(2.0*h)-base.gradient[j]).abs()<1e-9);
        for (r,(a,b)) in a.ci.iter().zip(&b.ci).enumerate() {
            assert!(((a-b)/(2.0*h)-base.ji[r*width+j]).abs()<1e-9);
        }
        for (r,(a,b)) in a.ce.iter().zip(&b.ce).enumerate() {
            assert!(((a-b)/(2.0*h)-base.je[r*width+j]).abs()<1e-9);
        }
    }
}

#[test]
fn tail_profile_admission_and_threshold_reporting_do_not_invent_a_quantile() {
    let loaded=load(MODEL,DESIGN);let p=loaded.problem();let offsets=[0.0;4];
    for alpha in [f64::NAN,f64::INFINITY,-1.0,0.0,1.0] {
        assert!(ensemble(p,&offsets,64).with_cvar(alpha).is_err());
    }
    // n=1,s=4,c=0: six decisions plus ten inequality rows need 16 slots.
    assert!(ensemble(p,&offsets,15).with_cvar(0.4).is_err());
    let e=ensemble(p,&offsets,16).with_cvar(0.4).unwrap();assert_eq!(e.cvar_alpha(),Some(0.4));
    let mut result=e.evaluate(&[0.0],&mut DesignControl::new(4,4),&CancelGate::new()).unwrap();
    for (case,value) in result.scenarios.iter_mut().zip([0.0,1.0,2.0,3.0]) { case.value=value; }
    assert!((e.cvar_upper_bound(&result,1.0).unwrap()-2.25).abs()<1e-14);
    assert_eq!(e.cvar_upper_bound(&result,3.0).unwrap(),3.0); // not the minimizing threshold
    assert!(e.cvar_upper_bound(&result,f64::NAN).is_err());
    result.scenarios[0].value=f64::INFINITY;assert!(e.cvar_upper_bound(&result,1.0).is_err());
    result.scenarios.pop();assert!(e.cvar_upper_bound(&result,1.0).is_err());
    assert_eq!(ensemble(p,&offsets,16).cvar_alpha(),None);
}

#[test]
fn cancellation_budget_extension_and_split_runs_keep_slacks_with_the_accepted_physics() {
    let loaded=load(MODEL,DESIGN);let gate=CancelGate::new();
    let e=ensemble(loaded.problem(),&[0.0,0.0,0.0,4.0],32).with_cvar(0.4).unwrap();
    let mut wa=DesignControl::new(1600,1600);let mut wb=DesignControl::new(4,4);
    let mut a=ScenarioEquilibriumStudy::new(e.clone(),&[0.0],&mut wa,&gate).unwrap();
    let mut b=ScenarioEquilibriumStudy::new(e,&[0.0],&mut wb,&gate).unwrap();
    assert_eq!(b.run(1e-9,128,1,&gate).unwrap().stop,SqpStop::EvaluationLimit);
    let point=b.optimizer().point().to_vec();let evidence=b.accepted().clone();let work=b.work();
    let cancelled=CancelGate::new();cancelled.request();
    assert!(matches!(b.run(1e-9,1,400,&cancelled),Err(SqpError::Cancelled)));
    assert_eq!(b.optimizer().point(),point);assert_eq!(b.accepted(),&evidence);assert_eq!(b.work(),work);
    b.extend_physics_budget(1600,1600).unwrap();
    a.run(1e-9,80,400,&gate).unwrap();b.run(1e-9,2,400,&gate).unwrap();b.run(1e-9,78,400,&gate).unwrap();
    assert_eq!(a.optimizer().point(),b.optimizer().point());assert_eq!(a.optimizer().history(),b.optimizer().history());
    assert_eq!(a.accepted(),b.accepted());assert_eq!(a.work(),b.work());
}

#[test]
fn late_physical_refusals_cannot_be_hidden_in_an_unselected_tail() {
    let limited=MODEL.replace("limits 0.9 10 1000 1000 1000000","limits 0.9 10 1000 3 1000000");
    let loaded=load(&limited,DESIGN);let gate=CancelGate::new();
    let e=ensemble(loaded.problem(),&[0.0,0.0,0.0,4.0],32).with_cvar(0.4).unwrap();
    let mut work=DesignControl::new(100,100);
    assert!(matches!(ScenarioEquilibriumStudy::new(e,&[0.0],&mut work,&gate),
        Err(SqpError::Evaluation(ScenarioError::Design {scenario:Some(3),source:DesignError::Case {..}}))));
    assert_eq!(work.work(),DesignWork {evaluations:4,case_solves:4});
}
