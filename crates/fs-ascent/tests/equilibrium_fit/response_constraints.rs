//! Invoke the actual command; reference values are independent physical algebra.
use super::*;
const LIMITED: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../fs-couple/examples/equilibrium-response-limits.fit"));
fn constraint<'a>(text: &'a str, name: &str) -> &'a str {
    text.split_once("\"constraints\":[").unwrap().1
        .split_once(&format!("\"name\":\"{name}\"")).unwrap().1.split('}').next().unwrap()
}

#[test]
fn command_changes_load_to_meet_a_physical_limit_not_a_parameter_bound() {
    let (m,d)=inputs(&directory(),MODEL,LIMITED);let result=text(run(&m,&d,&[]));
    assert!(result.contains("\"stop\":\"Converged\"") && result.contains("\"converged\":true"));
    let expected=0.8+400.0*(0.0001+0.8/10000.0+(0.8_f64/1e8).sqrt());
    assert!((parameter(&result,"load-N")-expected).abs()<1e-7);
    assert!(number_after(&result,"\"upper_multiplier_decision\":").abs()<1e-8);
    assert!(number_after(&result,"\"lower_multiplier_decision\":").abs()<1e-8);
    assert!((number_after(&result,"\"objective\":")-0.72).abs()<1e-7);
    let cap=constraint(&result,"normal-cap");
    assert!((number_after(cap,"\"value\":")-0.8).abs()<1e-8);
    assert!((number_after(cap,"\"multiplier_normalized\":")-1.2).abs()<1e-7);
    for name in ["normal-cap","indentation-cap","support-floor","travel-cap"] {
        assert!(number_after(constraint(&result,name),"\"violation\":")<=1e-8);
    }
    assert!(constraint(&result,"support-floor").contains("\"unit\":\"N\""));
    assert!(constraint(&result,"travel-cap").contains("\"unit\":\"m\""));
    let unconstrained=LIMITED.split_once("constraint_limits").unwrap().0;
    let (m,d)=inputs(&directory(),MODEL,unconstrained);let legacy=text(run(&m,&d,&[]));
    assert!(!legacy.contains("\"constraints\":"));
    assert!(parameter(&legacy,"load-N")>2.0 && number_after(&legacy,"\"objective\":")<1e-8);
}

#[test]
fn equality_and_lower_response_rows_use_their_correct_dual_families() {
    for (sense,target,multiplier) in [("equal","0.0002",1.2),("at-least","0",0.8)] {
        // Retain the other three (inactive) physical limits.
        let design=LIMITED.replace("contact-force 0 at-most",&format!("contact-force 0 {sense}"))
            .replace("target 1 0 0.0002 0.0001",&format!("target 1 0 {target} 0.0001"));
        let (m,d)=inputs(&directory(),MODEL,&design);let result=text(run(&m,&d,&[]));
        assert!(result.contains("\"converged\":true"));
        let row=constraint(&result,"normal-cap");
        assert!(row.contains(&format!("\"sense\":\"{sense}\"")));
        assert!((number_after(row,"\"multiplier_normalized\":")-multiplier).abs()<1e-6);
        assert!(number_after(row,"\"residual\":").abs()<1e-8);
        assert!(number_after(&result,"\"stationarity\":")<=1e-8);
    }
}

#[test]
fn budgeted_infeasibility_and_invalid_rows_are_reported_without_partial_or_relocated_changes() {
    let dir=directory();let (m,d)=inputs(&dir,MODEL,LIMITED);
    let partial=text(run(&m,&d,&["--evaluations","2"]));
    assert!(partial.contains("\"converged\":false") && partial.contains("\"stop\":\"EvaluationLimit\""));
    assert!(number_after(constraint(&partial,"normal-cap"),"\"violation\":")>0.5);
    assert_eq!(number_after(&partial,"\"evaluations_including_audit\":"),2.0);
    let (other_m,other_d)=inputs(&directory(),MODEL,LIMITED);
    assert_eq!(partial,text(run(&other_m,&other_d,&["--evaluations","2"])));
    let bad=run(&m,&d,&["--max-kkt-dimension","6"]);assert!(!bad.status.success() && bad.stdout.is_empty());
    assert_eq!(std::fs::read_to_string(m).unwrap(),MODEL);assert_eq!(std::fs::read_to_string(d).unwrap(),LIMITED);
    for malformed in [LIMITED.replace("contact-force 0","contact-force 99"),
        LIMITED.replace("at-most 0.8 1","at-most 0.8 0"),LIMITED.replace("constraints 4","constraints 5")] {
        let (m,d)=inputs(&directory(),MODEL,&malformed);let bad=run(&m,&d,&[]);
        assert!(!bad.status.success() && bad.stdout.is_empty());
    }
}
