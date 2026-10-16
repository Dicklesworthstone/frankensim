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

#[test]
fn retained_nominal_and_prefix_recover_without_a_second_adjoint_or_refit() {
    use fs_blake3::ContentHash;
    let identity = ContentHash([47; 32]);
    let base = J::parse(BASE).unwrap();
    let config = Config::parse(PLAN, &base).unwrap();
    let mut run = UqExecution::new(&config.plan()).unwrap();
    let mut nominal_calls = 0;
    let control = prepare_with(&base, &config, &run, |_| {
        nominal_calls += 1;
        Ok(document())
    }).unwrap();
    // Persist before the first random observation, not just after some samples.
    let zero = control.checkpoint_bytes(&run, identity).unwrap();
    let (empty, saved) = MeanControl::restore_bytes(&config.plan(), identity, &zero).unwrap();
    assert_eq!(empty.evaluations_attempted(), 0);
    assert_eq!(saved.nominal_temperature.to_bits(), control.nominal_temperature.to_bits());
    assert_eq!(saved.adjoint_residual.to_bits(), control.adjoint_residual.to_bits());
    assert_eq!(saved.frozen.gradient(), control.frozen.gradient());
    assert_eq!(saved.checkpoint_bytes(&empty, identity).unwrap(), zero);
    run.advance(3, || false, |p| Ok::<_, &str>(p[0] + 5.0));
    let bytes = control.checkpoint_bytes(&run, identity).unwrap();
    let raw_before = run.checkpoint(identity).unwrap();
    let (mut restored, frozen) = MeanControl::restore_bytes(&config.plan(), identity, &bytes).unwrap();
    assert_eq!(restored.checkpoint(identity).unwrap(), raw_before);
    let mut sample_calls = 0;
    restored.advance(usize::MAX, || false, |p| { sample_calls += 1; Ok::<_, &str>(p[0] + 5.0) });
    assert_eq!(sample_calls, 5);
    assert_eq!(nominal_calls, 1);
    run.advance(usize::MAX, || false, |p| Ok::<_, &str>(p[0] + 5.0));
    let raw = "{\"raw\":true}\n".to_string();
    let deadline = Instant::now() + Duration::from_secs(60);
    assert_eq!(frozen.attach(raw.clone(), &restored, deadline).unwrap(), control.attach(raw, &run, deadline).unwrap());
    assert_eq!(frozen.checkpoint_bytes(&restored, identity).unwrap(), control.checkpoint_bytes(&run, identity).unwrap());
    // A completed recovered run retains the same work count and raw data.
    let complete = frozen.checkpoint_bytes(&restored, identity).unwrap();
    let (mut terminal, retained) = MeanControl::restore_bytes(&config.plan(), identity, &complete).unwrap();
    terminal.advance(usize::MAX, || false, |_| -> std::result::Result<f64, &str> {
        panic!("terminal controlled execution must not run more physics")
    });
    assert_eq!(retained.checkpoint_bytes(&terminal, identity).unwrap(), complete);
}

#[test]
fn nominal_diagnostics_coefficients_and_control_mode_are_bound_to_the_checkpoint() {
    use fs_blake3::ContentHash;
    let identity = ContentHash([48; 32]);
    let base = J::parse(BASE).unwrap();
    let config = Config::parse(PLAN, &base).unwrap();
    let mut run = UqExecution::new(&config.plan()).unwrap();
    let control = prepare_with(&base, &config, &run, |_| Ok(document())).unwrap();
    run.advance(3, || false, |p| Ok::<_, &str>(p[0] + 5.0));
    let bytes = control.checkpoint_bytes(&run, identity).unwrap();
    for (offset, value) in [(8, 316.0), (16, 2e-14), (16, -1.0), (8, f64::NAN), (40, 2.0)] {
        let mut changed = bytes.clone();
        changed[offset..offset + 8].copy_from_slice(&value.to_bits().to_le_bytes());
        let error = MeanControl::restore_bytes(&config.plan(), identity, &changed).err().unwrap();
        assert_eq!(error.code, "cooling-network-uq-checkpoint");
    }
    for end in [0, 7, 8, 16, 23, 24, bytes.len() - 1] {
        assert!(MeanControl::restore_bytes(&config.plan(), identity, &bytes[..end]).is_err());
    }
    let mut changed_plan = config.plan(); changed_plan.budget_max_samples += 1;
    assert!(MeanControl::restore_bytes(&changed_plan, identity, &bytes).is_err());
    assert!(MeanControl::restore_bytes(&config.plan(), ContentHash([49; 32]), &bytes).is_err());
    assert!(MeanControl::restore_bytes(&config.plan(), identity, &run.checkpoint(identity).unwrap()).is_err());
    assert!(UqExecution::restore(&config.plan(), identity, &bytes).is_err());
    let raw_options = Options::default();
    let controlled_options = Options { adjoint_mean_control: true, ..Options::default() };
    assert_ne!(raw_options.checkpoint_binding(&config).unwrap(), controlled_options.checkpoint_binding(&config).unwrap());
    // Invocation work allowances do not replace the original plan/model binding.
    let longer = Config::parse(&PLAN.replace("\"wall_seconds\":60", "\"wall_seconds\":120"), &base).unwrap();
    assert_eq!(controlled_options.checkpoint_binding(&config).unwrap(), controlled_options.checkpoint_binding(&longer).unwrap());
    assert_eq!(config.plan(), longer.plan());
}

#[test]
fn failed_samples_and_cancelled_assessments_do_not_replace_the_retained_control() {
    use fs_blake3::ContentHash;
    let identity = ContentHash([50; 32]);
    let base = J::parse(BASE).unwrap();
    let config = Config::parse(PLAN, &base).unwrap();
    let mut run = UqExecution::new(&config.plan()).unwrap();
    let control = prepare_with(&base, &config, &run, |_| Ok(document())).unwrap();
    run.advance(8, || false, |p| Ok::<_, &str>(p[0] + 5.0));
    let bytes = control.checkpoint_bytes(&run, identity).unwrap();
    assert!(control.attach("{}\n".into(), &run, Instant::now()).is_err());
    assert_eq!(control.checkpoint_bytes(&run, identity).unwrap(), bytes);
    let (mut failed, frozen) = MeanControl::restore_bytes(&config.plan(), identity,
        &control.checkpoint_bytes(&UqExecution::new(&config.plan()).unwrap(), identity).unwrap()).unwrap();
    failed.advance(1, || false, |_| Err::<f64, _>("sample-refusal"));
    assert!(frozen.checkpoint_bytes(&failed, identity).is_err());
}
