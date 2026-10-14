use super::*;
use fs_bem::reciprocal_radiation::{Prepared,Budget};
fn sphere()->(SpherePanels,Vec<Vec<f64>>) {
    let s=SpherePanels::icosphere(1.0,0).unwrap();
    let rows=vec![vec![1.0;s.areas().len()],vec![-0.7;s.areas().len()]];
    (s,rows)
}
#[test]
fn velocity_work_and_acceleration_pressure_keep_distinct_physical_units() {
    let medium=Medium{density:1.2,sound_speed:343.0};let gate=CancelGate::new_clock_free();
    let omega=[34.3,102.9,171.5];let receiver=[0.2,0.3,2.0];
    let (surface,rows)=sphere();
    let result=sample(&surface,&rows,1.0,medium,&omega,receiver,10_000_000_000,&gate).unwrap();
    let source=Prepared::new(&surface,0.5,Budget{maximum_unknowns:128,maximum_cache_bytes:8*1024*1024,
        maximum_work:10_000_000_000,angular_refinement:0},||false).unwrap();
    let receivers=source.receivers(&[receiver],1e-8).unwrap();
    let refs:Vec<_>=rows.iter().map(Vec::as_slice).collect();
    for (i,&w) in omega.iter().enumerate() {
        let batch=source.solve_velocity_ports(w/medium.sound_speed,medium,&refs,||false).unwrap();
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
    let (surface,rows)=sphere();
    let p=surface.centroids()[0];let normal=surface.normals()[0];
    let receiver=core::array::from_fn(|j|p[j]+0.001*normal[j]);
    assert!(receiver.iter().map(|v|v*v).sum::<f64>()<1.0);
    let medium=Medium{density:1.2,sound_speed:343.0};let gate=CancelGate::new_clock_free();
    let result=sample(&surface,&rows,1.0,medium,&[34.3],receiver,10_000_000_000,&gate).unwrap();
    assert_eq!(result.delay_s,0.0);
    assert!(result.pressure.iter().all(|r|r[0].abs()>0.0));
    assert!(sample(&surface,&rows,1.0,medium,&[34.3],[0.0;3],10_000_000_000,&gate).is_err());
}
#[test]
fn budgets_band_and_cancelled_inputs_never_switch_equations() {
    let (surface,rows)=sphere();let medium=Medium{density:1.2,sound_speed:343.0};
    let gate=CancelGate::new_clock_free();let receiver=[0.0,0.0,3.0];
    assert!(sample(&surface,&rows,1.0,medium,&[3430.0],receiver,10_000_000_000,&gate).is_err());
    assert!(sample(&surface,&rows,1.0,medium,&[34.3],receiver,1,&gate).is_err());
    assert!(sample(&surface,&rows,1.0,medium,&[34.3,1.0],receiver,10_000_000_000,&gate).is_err());
    gate.request();
    assert!(sample(&surface,&rows,1.0,medium,&[34.3],receiver,10_000_000_000,&gate).is_err());
}
