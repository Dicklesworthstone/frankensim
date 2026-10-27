use super::*;
fn near(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<tol,"{a:e} != {b:e}");}

#[test]
fn one_source_card_drives_mass_decay_and_caloric_rise_against_independent_oscillators() {
    for (rho,young,time,cp) in [(2.,120.,0.01,1.),(4.,120.,0.01,1.),(2.,240.,0.02,2.)] {
        let card=demo_data(rho,young,young,Some(time),cp).unwrap();
        let run=simulate(&card,DEMO_TIME,8,0.03125,&CancelGate::new_clock_free()).unwrap();
        let mass=rho/6.; let m=rho/60.; let k=young/6.;
        let decay=time*k/(2.*m);
        let omega=(k/m-decay*decay).sqrt();
        let energy0=0.5*m*INITIAL_V_M_S.powi(2)+0.5*k*INITIAL_Q_M.powi(2);
        near(run.reference_mass_kg,mass,1e-12);
        assert_eq!(run.card,card.content_hash());
        assert_eq!(run.samples.len(),9);
        for sample in &run.samples {
            let t=sample.time_s;
            let b=(INITIAL_V_M_S+decay*INITIAL_Q_M)/omega;
            let shape=INITIAL_Q_M*(omega*t).cos()+b*(omega*t).sin();
            let q=(-decay*t).exp()*shape;
            let v=(-decay*t).exp()*(-decay*shape-INITIAL_Q_M*omega*(omega*t).sin()+b*omega*(omega*t).cos());
            let energy=0.5*m*v*v+0.5*k*q*q;
            near(sample.displacement_m,q,1e-10);
            near(sample.velocity_m_s,v,1e-9);
            near(sample.mechanical_energy_j,energy,1e-10);
            near(sample.material_heat_j,energy0-energy,1e-10);
            near(sample.temperature_k,300.+(energy0-energy)/(mass*cp),1e-8);
            near(sample.caloric_energy_change_j,energy0-energy,1e-10);
            near(sample.combined_residual_j,0.,1e-9);
        }
    }
}

#[test]
fn heat_capacity_only_swap_changes_temperature_not_the_retained_mechanics() {
    let gate=CancelGate::new_clock_free();
    let a=simulate(&demo(1.,1.).unwrap(),DEMO_TIME,8,0.03125,&gate).unwrap();
    let b=simulate(&demo(1.,2.).unwrap(),DEMO_TIME,8,0.03125,&gate).unwrap();
    let c=simulate(&demo(3.,1.).unwrap(),DEMO_TIME,8,0.03125,&gate).unwrap();
    for (x,y) in a.samples.iter().zip(&b.samples) {
        assert_eq!(x.displacement_m,y.displacement_m);
        assert_eq!(x.velocity_m_s,y.velocity_m_s);
        assert_eq!(x.material_heat_j,y.material_heat_j);
        near(x.temperature_k-300.,2.*(y.temperature_k-300.),1e-9);
    }
    assert!((a.samples[4].displacement_m-c.samples[4].displacement_m).abs()>1e-6);
    assert!((a.samples[4].temperature_k-c.samples[4].temperature_k).abs()>1e-6);
    let repeated=simulate(&demo(1.,1.).unwrap(),DEMO_TIME,8,0.03125,&gate).unwrap();
    assert_eq!(a.samples,repeated.samples);
}

#[test]
fn unprovided_loss_or_changing_mechanical_sources_refuse_instead_of_using_defaults() {
    let gate=CancelGate::new_clock_free();
    for card in [
        demo_data(2.,120.,120.,None,1.).unwrap(),
        demo_data(2.,120.,120.,Some(-0.01),1.).unwrap(),
        demo_data(2.,120.,121.,Some(0.01),1.).unwrap(),
    ] {
        assert!(simulate(&card,DEMO_TIME,8,0.03125,&gate).is_err());
    }
    let card=demo(1.,1.).unwrap();
    assert!(simulate(&card,"density",8,0.03125,&gate).is_err());
    assert!(simulate(&card,DEMO_TIME,0,0.03125,&gate).is_err());
    assert!(simulate(&card,DEMO_TIME,4097,0.03125,&gate).is_err());
    assert!(simulate(&card,DEMO_TIME,1,f64::NAN,&gate).is_err());
    gate.request();
    assert!(simulate(&card,DEMO_TIME,8,0.03125,&gate).is_err());
}
