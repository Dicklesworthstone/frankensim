use super::*;
fn args(s:&str)->Vec<String> {s.split_whitespace().map(str::to_owned).collect()}
fn controls()->Controls {Controls {angular_refinement:0,maximum_sweep_work:10_000_000_000}}
#[test]
fn option_preserves_tensioner_and_receiver_controls_and_rejects_atomically() {
    let mut a=args("shape.fsk 100 --edge-force 1 2 --reciprocal-air 1 10000000000 --wav 20 40 500 32 12 0 0 1");
    assert_eq!(option(&mut a).unwrap(),Some(Controls{angular_refinement:1,..controls()}));
    assert_eq!(a,args("shape.fsk 100 --edge-force 1 2 --wav 20 40 500 32 12 0 0 1"));
    for text in ["--reciprocal-air", "--reciprocal-air 3 1000000", "--reciprocal-air 1 0",
        "--reciprocal-air NaN 1000", "--reciprocal-air 1 100 --reciprocal-air 1 100"] {
        let mut a=args(text);let original=a.clone();assert!(option(&mut a).is_err());assert_eq!(a,original);
    }
}
fn sphere()->Boundary {
    let s=SpherePanels::icosphere(1.0,0).unwrap();
    Boundary {triangles:s.triangles().unwrap().to_vec(),
        weights:vec![vec![1.0;s.areas().len()],vec![-0.7;s.areas().len()]],
        radius_m:1.0,minimum_box_width_m:2.0}
}
#[test]
fn actual_reciprocal_samples_keep_impedance_and_acceleration_pressure_units_separate() {
    let medium=Medium{density:1.2,sound_speed:343.0};let gate=CancelGate::new_clock_free();
    let omega=[34.3,102.9,171.5];let receiver=[0.2,0.3,2.0];
    let result=sample(sphere(),medium,&omega,receiver,controls(),&gate).unwrap();
    let b=sphere();let surface=SpherePanels::from_triangles(b.triangles).unwrap();
    let source=Prepared::new(&surface,0.5,Budget{maximum_unknowns:128,maximum_cache_bytes:8*1024*1024,
        maximum_work:1_000_000_000,angular_refinement:0},||false).unwrap();
    let receivers=source.receivers(&[receiver],1e-8).unwrap();
    let rows:Vec<_>=b.weights.iter().map(Vec::as_slice).collect();
    for (i,&w) in omega.iter().enumerate() {
        let batch=source.solve_velocity_ports(w/medium.sound_speed,medium,&rows,||false).unwrap();
        assert_eq!(result.impedances[i],batch.impedance());
        let observed=receivers.prepare(w/medium.sound_speed,medium,Options::default(),false,||false).unwrap();
        let p=observed.evaluate(&batch.solutions()[0]).unwrap().pressure[0];
        let restored=phase(result.pressure[0][i],w*result.delay_s)*C64::new(0.0,-w);
        assert!((restored-p).abs()<1e-10*(1.0+p.abs()));
        assert!((result.pressure[1][i]-result.pressure[0][i].scale(-0.7)).abs()
            <1e-10*(1.0+result.pressure[0][i].abs()));
    }
}
#[test]
fn close_exterior_microphone_is_not_rejected_by_the_enclosing_sphere() {
    let b=sphere();let surface=SpherePanels::from_triangles(b.triangles.clone()).unwrap();
    let p=surface.centroids()[0];let normal=surface.normals()[0];
    let receiver=core::array::from_fn(|j|p[j]+0.001*normal[j]);
    assert!(receiver.iter().map(|v|v*v).sum::<f64>()<1.0);
    let medium=Medium{density:1.2,sound_speed:343.0};let gate=CancelGate::new_clock_free();
    let result=sample(b,medium,&[34.3],receiver,controls(),&gate).unwrap();
    assert_eq!(result.delay_s,0.0);
    assert!(result.pressure.iter().all(|r|r[0].abs()>0.0));
    assert!(sample(sphere(),medium,&[34.3],[0.0;3],controls(),&gate).is_err());
}
#[test]
fn band_budget_and_cancellation_refuse_without_switching_to_collocation() {
    let medium=Medium{density:1.2,sound_speed:343.0};let gate=CancelGate::new_clock_free();
    assert!(sample(sphere(),medium,&[3430.0],[0.0,0.0,3.0],controls(),&gate).is_err());
    assert!(sample(sphere(),medium,&[34.3],[0.0,0.0,3.0],Controls{maximum_sweep_work:1,..controls()},&gate).is_err());
    gate.request();
    assert!(sample(sphere(),medium,&[34.3],[0.0,0.0,3.0],controls(),&gate).is_err());
}
