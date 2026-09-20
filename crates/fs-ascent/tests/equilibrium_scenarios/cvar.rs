use super::*;
const LINEAR_MODEL: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),"/../fs-couple/examples/equilibrium-tail.model"));
const LINEAR_DESIGN: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),"/../fs-couple/examples/equilibrium-tail.fit"));
const TAIL: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),"/../fs-couple/examples/equilibrium-tail.scenarios"));

#[test]
fn authored_fractional_tail_selects_a_different_design_than_default_minimax() {
    let (m,d,s)=inputs(LINEAR_MODEL,LINEAR_DESIGN,TAIL);
    let result=success(run(&m,&d,Some(&s),&["--tolerance","1e-9"]));
    assert!(result.contains("\"schema\":\"frankensim-equilibrium-cvar-fit-v1\""));
    assert!(result.contains("\"converged\":true"));
    assert!((load_parameter(&result)-7.0/3.0).abs()<1e-6);
    assert!((number(&result,"\"cvar_upper_bound\":")-35.0/288.0).abs()<1e-8);
    assert!((number(&result,"\"objective\":")-35.0/288.0).abs()<1e-8);
    assert!((number(&result,"\"tail_mass\":")-2.4).abs()<1e-14);
    assert_eq!(number(&result,"\"scenario_mass\":"),0.25);
    assert_eq!(result.matches("\"tail_multiplier\":").count(),4);
    let high=result.split_once("\"name\":\"high-load\"").unwrap().1;
    assert!((number(high,"\"tail_multiplier\":")-5.0/12.0).abs()<1e-7);
    assert!((number(high,"\"tail_excess\":")-1.0/12.0).abs()<1e-7);
    let maximum=TAIL.replace("risk empirical-cvar 0.4\n","");
    let (mm,dd,ss)=inputs(LINEAR_MODEL,LINEAR_DESIGN,&maximum);
    let worst=success(run(&mm,&dd,Some(&ss),&["--tolerance","1e-9"]));
    assert!(worst.contains("\"scope\":\"local-static-finite-scenario-minimax\""));
    assert!((load_parameter(&worst)-2.0).abs()<1e-7);
    assert!(!worst.contains("cvar_upper_bound"));
}

#[test]
fn a_non_tail_contact_force_limit_keeps_its_own_physical_multiplier() {
    let cases="frankensim-equilibrium-scenarios-v1\nrisk empirical-cvar 0.4\nscenarios 4\nscenario low-load\noffset load-N -0.1\nscenario normal-a\noffset load-N 0\nscenario normal-b\noffset load-N 0\nscenario high-load\noffset load-N 0.1\n";
    let (m,d,s)=inputs(MODEL,DESIGN,cases);
    let result=success(run(&m,&d,Some(&s),&[]));
    assert!(result.contains("\"converged\":true"));
    assert!((load_parameter(&result)-0.8077770876399966).abs()<1e-6);
    assert!((number(&result,"\"cvar_upper_bound\":")-0.8897892300777).abs()<1e-6);
    assert_eq!(result.matches("\"quantity\":\"contact-force\"").count(),4);
    let high=result.split_once("\"name\":\"high-load\"").unwrap().1;
    assert!(number(high,"\"tail_multiplier\":").abs()<1e-7);
    let cap=high.split_once("\"name\":\"normal-cap\"").unwrap().1;
    assert!((number(cap,"\"value\":")-0.8).abs()<1e-7);
    assert!(number(cap,"\"multiplier_normalized\":")>0.5);
    for row in result.split("\"violation\":").skip(1) {
        assert!(row.split(',').next().unwrap().parse::<f64>().unwrap()<=1e-7);
    }
}

#[test]
fn budget_stop_reports_a_recomputed_upper_bound_not_a_fabricated_optimal_quantile() {
    let (m,d,s)=inputs(LINEAR_MODEL,LINEAR_DESIGN,TAIL);
    let result=success(run(&m,&d,Some(&s),&["--evaluations","8"]));
    assert!(result.contains("\"stop\":\"EvaluationLimit\""));assert!(result.contains("\"converged\":false"));
    assert_eq!(number(&result,"\"iterations\":"),0.0);
    assert_eq!(number(&result,"\"physical_evaluations_including_audit\":"),8.0);
    assert_eq!(number(&result,"\"ensemble_evaluations_including_audit\":"),2.0);
    assert_eq!(number(&result,"\"case_solves\":"),8.0);
    // eta=max loss is feasible, but is not a CVaR-minimizing threshold at this design.
    assert_eq!(number(&result,"\"threshold\":"),0.5);
    assert_eq!(number(&result,"\"cvar_upper_bound\":"),0.5);
    assert_eq!(result.matches("\"tail_excess\":").count(),4);
    for flags in [vec!["--evaluations","7"],vec!["--max-kkt-dimension","15"]] {
        let refusal=run(&m,&d,Some(&s),&flags);assert!(!refusal.status.success());assert!(refusal.stdout.is_empty());
    }
}

#[test]
fn risk_inputs_are_strict_relocatable_and_never_modify_the_source_model() {
    let (m,d,s)=inputs(LINEAR_MODEL,LINEAR_DESIGN,TAIL);
    let result=success(run(&m,&d,Some(&s),&["--evaluations","8"]));
    let (mm,dd,ss)=inputs(LINEAR_MODEL,LINEAR_DESIGN,TAIL);
    assert_eq!(result,success(run(&mm,&dd,Some(&ss),&["--evaluations","8"])));
    assert_eq!(std::fs::read_to_string(m).unwrap(),LINEAR_MODEL);
    assert_eq!(std::fs::read_to_string(d).unwrap(),LINEAR_DESIGN);
    assert_eq!(std::fs::read_to_string(s).unwrap(),TAIL);
    for bad in [TAIL.replace("0.4","NaN"),TAIL.replace("0.4","1"),TAIL.replace("empirical-cvar","mean"),
        TAIL.replace("scenarios 4","risk empirical-cvar 0.5\nscenarios 4"),format!("{TAIL}risk empirical-cvar 0.4\n")] {
        let (m,d,s)=inputs(LINEAR_MODEL,LINEAR_DESIGN,&bad);
        let refusal=run(&m,&d,Some(&s),&[]);assert!(!refusal.status.success());assert!(refusal.stdout.is_empty());
    }
}
