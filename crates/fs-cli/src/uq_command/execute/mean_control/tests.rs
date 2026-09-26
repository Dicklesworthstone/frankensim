use super::*;

const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/fan-correlated-hotspot.json"));
const PLAN: &str = r#"{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":8,"wall_seconds":60,"correlation":{"kind":"independent"},"parameters":[{"target":{"kind":"inlet-temperature","index":0},"distribution":{"kind":"uniform","lo":308,"hi":312}}]}"#;
fn document() -> J {
    J::parse(r#"{"schema":"frankensim.cooling-network.result.v1","objective":{"value_k":315},"adjoint_residual":1e-14,"dobjective_dinlet_k":[1,0.25],"fan_speed_sensitivity":{"status":"available","method":"fan-affinity-coupled-adjoint","dobjective_dlog_speed_ratio_k":-6},"walls":[{"region":"second","htc_w_m2_k":10,"dobjective_dlog_htc":-2},{"region":"first","htc_w_m2_k":50,"dobjective_dlog_htc":-5}],"contact_sensitivities":{"method":"coupled-adjoint-contact-bilinear-form","rows":[{"contact":"bondline","resistance_m2_k_w":0.01,"dobjective_dlog_resistance_k":2}]}}"#).unwrap()
}

#[test]
fn actual_parameter_means_and_gradient_request_are_fixed_before_sampling() {
    let base = J::parse(BASE).unwrap();
    let original = base.clone();
    let config = Config::parse(PLAN, &base).unwrap();
    let run = UqExecution::new(&config.plan()).unwrap();
    let mut calls = 0;
    let control = prepare_with(&base,&config,&run,|request| {
        calls += 1;
        let sent = J::parse(request).unwrap();
        assert_eq!(sent.path(&["hydraulics","fan","temperature_k"]).and_then(J::as_f64),Some(310.0));
        assert_eq!(sent.path(&["objective","gradient"]),Some(&J::Bool(true)));
        assert_eq!(sent.get("solid"),base.get("solid"));
        assert_eq!(sent.get("tolerances"),base.get("tolerances"));
        Ok(document())
    }).unwrap();
    assert_eq!(calls,1);
    assert_eq!(control.frozen.parameter_means(),&[310.0]);
    assert_eq!(control.frozen.gradient(),&[1.0]);
    assert_eq!(run.evaluations_attempted(),0,"nominal call is not a random observation");
    assert_eq!(base,original);
    assert_eq!(J::parse(&config.sample_request(&base,&[311.0]).unwrap()).unwrap()
        .path(&["objective","gradient"]),Some(&J::Bool(false)));
}

#[test]
fn named_total_derivatives_use_linear_coordinates_not_log_coordinates() {
    let d=document();
    assert_eq!(Target::Inlet(1).derivative(&d,300.0).unwrap(),0.25);
    assert_eq!(Target::FanSpeed.derivative(&d,2.0).unwrap(),-3.0);
    assert_eq!(Target::SurfaceHtc("first".into()).derivative(&d,50.0).unwrap(),-0.1);
    assert_eq!(Target::SurfaceHtc("second".into()).derivative(&d,10.0).unwrap(),-0.2);
    assert_eq!(Target::ContactResistance("bondline".into()).derivative(&d,0.01).unwrap(),200.0);
    assert!(Target::SurfaceHtc("first".into()).derivative(&d,40.0).is_err());
    assert!(Target::ContactResistance("bondline".into()).derivative(&d,0.02).is_err());
    assert!(Target::Inlet(2).derivative(&d,300.0).is_err());
    assert!(Target::SurfaceHtc("missing".into()).derivative(&d,50.0).is_err());
    let duplicate=J::parse(r#"[{"region":"same"},{"region":"same"}]"#).unwrap();
    assert!(named_row(&duplicate,"region","same").is_err());
    let missing=J::parse(r#"{"fan_speed_sensitivity":{"status":"unavailable"}}"#).unwrap();
    assert!(Target::FanSpeed.derivative(&missing,1.0).is_err());
    assert!(Target::Inlet(0).derivative(&J::Null,300.0).is_err());
    assert!(linear_derivative(f64::MAX,f64::MIN_POSITIVE).is_err());
}

#[test]
fn unsupported_targets_and_model_failures_never_become_zero_controls() {
    let base=J::parse(BASE).unwrap();
    let text=PLAN.replace(r#"{"kind":"inlet-temperature","index":0}"#,r#"{"kind":"air-density"}"#)
        .replace("308","1.1").replace("312","1.3");
    let config=Config::parse(&text,&base).unwrap();
    let run=UqExecution::new(&config.plan()).unwrap();
    let mut called=false;
    assert!(prepare_with(&base,&config,&run,|_| {called=true;Ok(document())}).is_err());
    assert!(!called);
    let config=Config::parse(PLAN,&base).unwrap();
    let mut run=UqExecution::new(&config.plan()).unwrap();
    let error=prepare_with(&base,&config,&run,|_|Err(model_failure("nominal-refusal-sentinel"))).err().unwrap();
    assert!(error.message.contains("nominal-refusal-sentinel"));
    run.advance(1,||false,|_|Ok::<_,&str>(315.0));
    assert!(prepare_with(&base,&config,&run,|_| {called=true;Ok(document())}).is_err());
    assert!(!called,"a late control must refuse before further physics");
}

#[test]
fn publication_keeps_raw_result_and_requires_complete_work_and_live_budget() {
    let base=J::parse(BASE).unwrap();
    let config=Config::parse(PLAN,&base).unwrap();
    let mut run=UqExecution::new(&config.plan()).unwrap();
    let control=prepare_with(&base,&config,&run,|_|Ok(document())).unwrap();
    let output="{\"unchanged_raw_field\":[1,2,3]}\n".to_string();
    assert!(control.attach(output.clone(),&run,Instant::now()+Duration::from_secs(60)).is_err());
    run.advance(8,||false,|p|Ok::<_,&str>(p[0]+5.0));
    let before=run.observations().to_vec();
    let result=control.attach(output.clone(),&run,Instant::now()+Duration::from_secs(60)).unwrap();
    let parsed=J::parse(&result).unwrap();
    assert_eq!(parsed.get("unchanged_raw_field"),J::parse(&output).unwrap().get("unchanged_raw_field"));
    let adjusted=parsed.get("mean_control_variate").unwrap();
    assert!((adjusted.f64_field("mean_k").unwrap()-315.0).abs()<1e-11);
    assert_eq!(adjusted.f64_field("total_model_evaluations"),Some(9.0));
    assert!(control.attach(output,&run,Instant::now()).is_err());
    assert_eq!(run.observations(),before);
}

#[test]
fn nominal_work_is_admitted_inside_the_total_model_call_cap() {
    let base=J::parse(BASE).unwrap();
    let text=PLAN.replace("\"samples\":8","\"samples\":10000");
    assert_ne!(text,PLAN);
    let config=Config::parse(&text,&base).unwrap();
    let run=UqExecution::new(&config.plan()).unwrap();
    let mut calls=0;
    assert!(prepare_with(&base,&config,&run,|_| {calls+=1;Ok(document())}).is_err());
    assert_eq!(calls,0);
}
