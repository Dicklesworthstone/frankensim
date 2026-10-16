use super::*;
use fs_couple::render::schedule::force::file::design::EquilibriumDesignFile;

const MODEL: &str = "frankensim-modal-performance-v2\nsample_rate_hz 48000\nsamples 1\nfull_scale_pa 1\nlimits 0.9 10 1000 1000 1000000\ncompile_limits 0 1\nvoices 1\nvoice retain-state 1 1\nmode 2 0.1 0 0 0 0\nport 0 1\ncoupling_limits 0 16384 0.9 1000 1000000 10000 1e-10 1e-11 1e-9\nconnections 0\nevents 0\n";
const DESIGN: &str = "frankensim-equilibrium-design-v1\npreload_limits 0 1 16384\nsensitivity_limits 0 16384 16384 0\ndesign_limits 1 1 1 2\ncases 1\ncase response 1 1\nload 0 0 4\ntarget 0 0 1 1 1\nvariables 1\nvariable load-N 4 4 1 8 1\nbind actuator-force 0 0\n";
const CONTACT_MODEL: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../fs-couple/examples/equilibrium-design.model"));
const CONTACT_DESIGN: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../fs-couple/examples/equilibrium-response-limits.fit"));
fn load(model: &str, design: &str) -> EquilibriumDesignFile {
    EquilibriumDesignFile::from_bytes(model.as_bytes(), design.as_bytes(), &CancelGate::new()).unwrap()
}
fn scenarios(offsets: &[f64]) -> Vec<EquilibriumScenario> {
    offsets.iter().enumerate().map(|(i, d)| EquilibriumScenario {
        name: format!("scenario-{i}"), physical_offsets: vec![*d],
    }).collect()
}

#[test]
fn minimax_balances_asymmetric_physical_tolerances_without_an_averaged_gradient() {
    // u = F/4 m. The nominal target is F=4 N. Offsets [-0.8,2.4] N
    // give the minimax nominal F=3.2 N, worst loss .08, and equal risk duals.
    let loaded = load(MODEL, DESIGN); let gate = CancelGate::new();
    let ensemble = ScenarioProblem::new(loaded.problem(), scenarios(&[-0.8,2.4]), 2, 16, &gate).unwrap();
    let (lo, hi) = ensemble.decision_bounds();
    assert!((lo[0]+0.55).abs()<1e-14 && (hi[0]-0.4).abs()<1e-14);
    let mut work = DesignControl::new(400,400);
    let mut study = ScenarioEquilibriumStudy::new(ensemble, &[0.0], &mut work, &gate).unwrap();
    let initial = study.accepted().worst_objective;
    let report = study.run(1e-9,128,200,&gate).unwrap();
    assert!(report.solution.converged, "{:?}: {:?}",report.stop,report.solution.kkt);
    let result = study.recheck(&gate).unwrap();
    assert!((result.nominal_parameters[0]-3.2).abs()<1e-7);
    assert!((result.worst_objective-0.08).abs()<1e-9 && result.worst_objective<initial);
    assert!((report.solution.x[1]-result.worst_objective).abs()<1e-9);
    assert!((report.solution.nu[2]-0.5).abs()<1e-7 && (report.solution.nu[3]-0.5).abs()<1e-7);
    for (actual,force) in result.scenarios.iter().zip([2.4,5.6]) {
        assert!((actual.physical_parameters[0]-force).abs()<1e-7);
        assert!((actual.cases[0].observations_m[0]-force/4.0).abs()<1e-8);
    }
}

#[test]
fn the_force_limited_scenario_is_not_required_to_be_the_worst_objective_scenario() {
    let loaded = load(CONTACT_MODEL, CONTACT_DESIGN); let gate=CancelGate::new();
    let ensemble = ScenarioProblem::new(loaded.problem(),scenarios(&[-0.1,0.1]),2,32,&gate).unwrap();
    let mut work=DesignControl::new(800,800);
    let mut study=ScenarioEquilibriumStudy::new(ensemble,&[0.0],&mut work,&gate).unwrap();
    assert!(study.accepted().scenarios.iter().any(|s|s.constraints[0].residual>0.0));
    let report=study.run(1e-8,128,400,&gate).unwrap();
    assert!(report.solution.converged,"{:?}: {:?}",report.stop,report.solution.kkt);
    let result=study.recheck(&gate).unwrap();
    let force=0.8+400.0*(0.0001+(0.8_f64/1e8).sqrt()+0.8/10000.0)-0.1;
    assert!((result.nominal_parameters[0]-force).abs()<1e-6);
    assert_eq!(result.worst_objective,result.scenarios[0].value);
    assert!(result.scenarios[1].value<result.worst_objective);
    for s in &result.scenarios { assert!(s.constraints.iter().all(|r|r.residual<=1e-7)); }
    assert!(result.scenarios[1].constraints[0].residual.abs()<1e-7);
    // Box rows [0,1]; each scenario has an epigraph then four response rows.
    assert!(report.solution.nu[2]>0.9 && report.solution.nu[7].abs()<1e-7);
    assert!(report.solution.nu[8]>0.1, "the high-force physical constraint must bind");
}

#[test]
fn source_constraint_rows_and_scaled_gradients_survive_scenario_and_epigraph_expansion() {
    let design=format!("{DESIGN}constraint_limits 3\nconstraints 3\nconstraint cap 0 displacement 0 0 at-most 1.2 2\nconstraint floor 0 displacement 0 0 at-least 0.1 3\nconstraint exact 0 displacement 0 0 equal 1 4\n");
    let loaded=load(MODEL,&design);let gate=CancelGate::new();
    let ensemble=ScenarioProblem::new(loaded.problem(),scenarios(&[-0.8,2.4]),2,16,&gate).unwrap();
    let mut work=DesignControl::new(2,2);let x=0.07;
    let evaluation=ensemble.evaluate(&[x],&mut work,&gate).unwrap();
    let s=ensemble.sample(&[x,0.5],&evaluation);
    assert_eq!(s.gradient,[0.0,1.0]);assert_eq!(s.f,0.5);
    assert_eq!(s.ce.len(),2);assert_eq!(s.ci.len(),8);
    for (i,delta) in [-0.2,0.6].iter().enumerate() {
        let u=1.0+x+delta;let start=2+i*3;
        assert!((s.ci[start]-(0.5*(u-1.0).powi(2)-0.5)).abs()<1e-14);
        assert!((s.ji[2*start]-(u-1.0)).abs()<1e-14);assert_eq!(s.ji[2*start+1],-1.0);
        assert!((s.ci[start+1]-(u-1.2)/2.0).abs()<1e-14);
        assert!((s.ci[start+2]-(0.1-u)/3.0).abs()<1e-14);
        assert_eq!(&s.ji[2*(start+1)..2*(start+1)+2],&[0.5,0.0]);
        assert_eq!(&s.ji[2*(start+2)..2*(start+2)+2],&[-1.0/3.0,0.0]);
        assert!((s.ce[i]-(u-1.0)/4.0).abs()<1e-14);
        assert_eq!(&s.je[2*i..2*i+2],&[0.25,0.0]);
    }
    // Equality feasibility across these two scenarios is deliberately NOT claimed.
}

#[test]
fn zero_tolerance_is_the_original_physical_evaluation_and_domain_admission_precedes_solves() {
    let loaded=load(MODEL,DESIGN);let gate=CancelGate::new();
    let nominal=loaded.problem().evaluate(&[0.1],&mut DesignControl::new(1,1),&gate).unwrap();
    let ensemble=ScenarioProblem::new(loaded.problem(),scenarios(&[0.0]),1,16,&gate).unwrap();
    let got=ensemble.evaluate(&[0.1],&mut DesignControl::new(1,1),&gate).unwrap();
    assert_eq!(got.scenarios,[nominal]);assert_eq!(got.worst_objective,got.scenarios[0].value);
    let ensemble=ScenarioProblem::new(loaded.problem(),scenarios(&[0.0,2.4]),2,16,&gate).unwrap();
    let mut work=DesignControl::new(100,100);
    assert!(matches!(ensemble.evaluate(&[0.9],&mut work,&gate),Err(ScenarioError::Design {scenario:Some(1),source:DesignError::OutsideBounds {..}})));
    assert_eq!(work.work(),DesignWork::default());
}

#[test]
fn late_failures_publish_no_family_and_retain_attempted_physics_work() {
    let limited=MODEL.replace("limits 0.9 10 1000 1000 1000000","limits 0.9 10 1000 3 1000000");
    let loaded=load(&limited,DESIGN);let gate=CancelGate::new();
    let ensemble=ScenarioProblem::new(loaded.problem(),scenarios(&[-0.8,2.4]),2,16,&gate).unwrap();
    let mut work=DesignControl::new(20,20);
    assert!(matches!(ensemble.evaluate(&[0.0],&mut work,&gate),Err(ScenarioError::Design {scenario:Some(1),source:DesignError::Case {..}})));
    assert_eq!(work.work(),DesignWork {evaluations:2,case_solves:2});
    let resumed=ensemble.evaluate(&[-0.4],&mut work,&gate).unwrap();
    let fresh=ensemble.evaluate(&[-0.4],&mut DesignControl::new(2,2),&gate).unwrap();assert_eq!(resumed,fresh);
    let mut short=DesignControl::new(1,1);
    assert!(matches!(ensemble.evaluate(&[-0.4],&mut short,&gate),Err(ScenarioError::Design {scenario:Some(1),source:DesignError::Budget {..}})));
    assert_eq!(short.work(),DesignWork {evaluations:1,case_solves:1});
}

#[test]
fn accepted_step_splits_cancellation_and_extended_budgets_preserve_the_same_study() {
    let loaded=load(MODEL,DESIGN);let gate=CancelGate::new();
    let ensemble=ScenarioProblem::new(loaded.problem(),scenarios(&[-0.8,2.4]),2,16,&gate).unwrap();
    let mut work_a=DesignControl::new(400,400);let mut work_b=DesignControl::new(2,2);
    let mut a=ScenarioEquilibriumStudy::new(ensemble.clone(),&[0.0],&mut work_a,&gate).unwrap();
    let mut b=ScenarioEquilibriumStudy::new(ensemble,&[0.0],&mut work_b,&gate).unwrap();
    assert_eq!(b.run(1e-9,128,1,&gate).unwrap().stop,SqpStop::EvaluationLimit);
    let before=b.accepted().clone();let work=b.work();let cancelled=CancelGate::new();cancelled.request();
    assert!(matches!(b.run(1e-9,1,200,&cancelled),Err(SqpError::Cancelled)));
    assert_eq!(b.accepted(),&before);assert_eq!(b.work(),work);
    b.extend_physics_budget(400,400).unwrap();
    a.run(1e-9,80,200,&gate).unwrap();b.run(1e-9,2,200,&gate).unwrap();b.run(1e-9,78,200,&gate).unwrap();
    assert_eq!(a.optimizer().point(),b.optimizer().point());assert_eq!(a.optimizer().history(),b.optimizer().history());
    assert_eq!(a.accepted(),b.accepted());assert_eq!(a.work(),b.work());
}

#[test]
fn invalid_scenarios_and_full_expanded_kkt_dimensions_refuse_before_any_evaluation() {
    let loaded=load(MODEL,DESIGN);let p=loaded.problem();let gate=CancelGate::new();
    for entries in [vec![],scenarios(&[f64::NAN]),scenarios(&[-100.0,100.0]),
        vec![EquilibriumScenario {name:"bad-width".into(),physical_offsets:vec![]}],
        vec![EquilibriumScenario {name:"duplicate".into(),physical_offsets:vec![0.0]};2]] {
        assert!(ScenarioProblem::new(p,entries,2,16,&gate).is_err());
    }
    // n=1, s=2, c=0 -> 3+1+2=6, including the auxiliary variable.
    assert!(ScenarioProblem::new(p,scenarios(&[0.0,0.0]),2,5,&gate).is_err());
    assert!(ScenarioProblem::new(p,scenarios(&[0.0,0.0]),2,6,&gate).is_ok());
    assert!(ScenarioProblem::new(p,scenarios(&[0.0]),0,16,&gate).is_err());
    gate.request();assert!(matches!(ScenarioProblem::new(p,scenarios(&[0.0]),1,16,&gate),Err(ScenarioError::Cancelled)));
}
