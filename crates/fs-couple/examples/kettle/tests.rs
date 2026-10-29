use super::{model::{self,Model},spec::Spec};
use fs_exec::CancelGate;
fn spec() -> Spec { Spec::read(include_str!("illustrative-parabolic.fsk")).unwrap() }
fn uniform() -> Spec {
    let mut s=spec();
    s.head.radial_intervals=2; s.head.azimuths=16;
    s.head_band_hz=[0.0,250.0]; s.air_maximum_hz=1.0;
    s.air_mesh=[2,0]; s.creep.clear();
    s
}
#[test]
fn specification_preserves_geometry_and_refuses_incomplete_or_duplicate_data() {
    let s=spec();
    assert_eq!(s.head.radius_m,0.3302); assert_eq!(s.profile.len(),9);
    let text=include_str!("illustrative-parabolic.fsk");
    for changed in [text.replace("head_mesh,4,32",""),format!("{text}\nrate_hz,768000"),
        text.replace("station,0,0","station,0,0.01"), text.replace("air,1.2,343,0","air,1.2,NaN,0"),
        text.replace("rate_hz,768000","rate_hz,0"),format!("{text}\nunknown,2")] {
        assert!(Spec::read(&changed).is_err());
    }
}
#[test]
fn complete_mallet_disk_cannot_cross_polygon_between_sampling_points() {
    let s=uniform(); let gate=CancelGate::new_clock_free();
    let m=Model::new(&s,1,&gate).unwrap();
    assert_eq!(model::footprint(&m.film,[0.0,0.0],0.01).unwrap().len(),4);
    assert!(model::footprint(&m.film,[s.head.radius_m-0.001,0.0],0.01).is_err());
}
#[test]
fn uniform_overlap_is_signed_actual_head_area_and_one_bulk_storage() {
    let s=uniform(); let gate=CancelGate::new_clock_free();
    let m=Model::new(&s,1,&gate).unwrap();
    let (air,overlap)=model::interface(&m.film,&m.modes,&m.air,&gate).unwrap();
    assert_eq!(air.omegas,vec![0.0]);
    assert_eq!(m.cavity.total_modes(),1+m.modes.len()); // no fake zero-mode inertia
    assert_eq!(overlap[0],0.0); // mallet never directly pumps the air
    for (i,mode) in m.modes.iter().enumerate() {
        let area=m.film.modal_area(&mode.phi).unwrap();
        assert!((overlap[i+1]+area).abs()<1e-12);
    }
    let mut state=m.mechanics.state().to_vec(); state[2]=1e-5;
    let expected=-s.air[0]*s.air[1].powi(2)*overlap[1]*state[2]/m.air.volume_m3();
    let mut p=[0.0]; m.cavity.pressures_into(&state,&mut p).unwrap();
    assert!((p[0]-expected).abs()<1e-10*expected.abs().max(1.0));
}
#[test]
fn nonlinear_head_air_and_felt_preserve_history_on_cancel_and_force_refusal() {
    let s=uniform(); let gate=CancelGate::new_clock_free();
    let mut a=Model::new(&s,800,&gate).unwrap(); let mut b=Model::new(&s,800,&gate).unwrap();
    let initial=a.mechanics.stored_energy_j(); let mut loss=0.0; let mut moved=false;
    for tick in 0..800 {
        if tick==200 {
            let before=a.mechanics.state().to_vec(); let history=a.mechanics.felt_history(0).unwrap();
            let cancelled=CancelGate::new_clock_free(); cancelled.request();
            assert!(a.mechanics.step(&a.force,&cancelled).is_err());
            let mut bad=a.force.clone(); bad[0]=f64::NAN;
            assert!(a.mechanics.step(&bad,&gate).is_err());
            assert_eq!(a.mechanics.state(),before);
            assert_eq!(a.mechanics.felt_history(0).unwrap(),history);
        }
        let frame=a.mechanics.step(&a.force,&gate).unwrap();
        b.mechanics.step(&b.force,&gate).unwrap();
        assert_eq!(a.mechanics.state(),b.mechanics.state());
        loss+=frame.dissipated_energy_j;
        assert!((a.mechanics.stored_energy_j()+loss-initial).abs()<1e-6);
        moved|=a.head_displacement().abs()>1e-12;
    }
    assert!(moved); assert!(a.mechanics.felt_history(0).unwrap().eps_max>0.0);
    let before=a.mechanics.state().to_vec();
    assert!(a.mechanics.step(&a.force,&gate).is_err()); assert_eq!(a.mechanics.state(),before);
}
#[test]
fn hysteresis_free_felt_law_strikes_the_head() {
    // F = K x^p unloading on its own curve (crush 0, q = p), as Rhaouti, Chaigne and Joly's mallet.
    let mut s=uniform(); s.felt=[100000.0,0.2,2.2,2.2,0.0,0.7,0.0];
    let gate=CancelGate::new_clock_free();
    let mut m=Model::new(&s,800,&gate).unwrap();
    let initial=m.mechanics.stored_energy_j(); let mut loss=0.0;
    for _ in 0..800 { loss+=m.mechanics.step(&m.force,&gate).unwrap().dissipated_energy_j; }
    assert!(m.mechanics.felt_history(0).unwrap().eps_max>0.0);
    assert!((m.mechanics.stored_energy_j()+loss-initial).abs()<1e-6);
    s.felt[3]=3.0; // crush 0 is the law's limit only when q = p
    assert!(Model::new(&s,800,&gate).is_err());
}
#[test]
fn an_explicit_zero_exterior_load_preserves_the_unloaded_model_exactly() {
    use fs_couple::render::plate::impact::radiation::{Model as Load,Pole};
    let s=uniform();let gate=CancelGate::new_clock_free();
    let mut a=Model::new(&s,40,&gate).unwrap();
    let mut b=Model::with_exterior(&s,40,&gate,|_,modes|Ok(Some(Load {
        ports:modes.len(),poles:vec![Pole {omega:1000.0,zeta:0.2,coupling:vec![0.0;modes.len()]}]
    }))).unwrap();
    for _ in 0..40 {
        a.mechanics.step(&a.force,&gate).unwrap(); b.mechanics.step(&b.force,&gate).unwrap();
        assert_eq!(a.mechanics.state(),b.mechanics.state());
    }
}
#[test]
fn passive_exterior_memory_joins_the_same_impact_and_rollback() {
    use fs_couple::render::plate::impact::radiation::{Model as Load,Pole};
    let s=uniform();let gate=CancelGate::new_clock_free();
    // Declared passive fixture tests composition; it is not a BEM sound model.
    let mut a=Model::with_exterior(&s,800,&gate,|_,modes|Ok(Some(Load {
        ports:modes.len(),poles:vec![Pole {omega:1000.0,zeta:0.2,coupling:vec![100.0;modes.len()]}]
    }))).unwrap();
    let initial=a.mechanics.stored_energy_j();let mut loss=0.0;
    for _ in 0..800 { loss+=a.mechanics.step(&a.force,&gate).unwrap().dissipated_energy_j; }
    assert!((a.mechanics.stored_energy_j()+loss-initial).abs()<1e-6);
    assert!(a.mechanics.radiation_observation().unwrap().stored_energy_j>0.0);
    let state=a.mechanics.state().to_vec();let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(a.mechanics.step(&a.force,&cancelled).is_err());assert_eq!(a.mechanics.state(),state);
}
