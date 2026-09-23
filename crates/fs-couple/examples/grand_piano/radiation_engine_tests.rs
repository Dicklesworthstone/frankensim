use super::*;
fn piano()->Instrument {
    let c=super::super::geometry::demonstration_scale().unwrap()[48];
    Instrument::new(vec![c],&super::super::board::demonstration(),48_000,4,12,true).unwrap()
}
fn load(p:&Instrument)->radiation::Model {
    radiation::Model {ports:p.bank.board_count,poles:vec![radiation::Pole {
        omega:std::f64::consts::TAU*300.,zeta:0.25,
        coupling:(0..p.bank.board_count).map(|j|200./(j+1) as f64).collect(),
    }]}
}
#[test]
fn hammer_contact_reacts_to_air_storage_and_closes_the_combined_energy() {
    let mut loaded=piano();let mut bare=piano();let model=load(&loaded);
    loaded.configure_radiation(&model).unwrap();assert!(loaded.has_radiation());
    loaded.note_on(69,1.).unwrap();bare.note_on(69,1.).unwrap();
    let mut peak_storage=0.0_f64;
    for _ in 0..2400 {
        loaded.step().unwrap();bare.step().unwrap();
        peak_storage=peak_storage.max(loaded.radiation_energy_j());
    }
    assert_ne!(loaded.bank.q,bare.bank.q);
    assert!(peak_storage>1e-15);assert!(loaded.accounting.radiation_loss_j>1e-15);
    assert_eq!(bare.accounting.radiation_loss_j,0.);
    let balance=loaded.accounting.input_work_j-loaded.accounting.dissipated_j()-loaded.energy_j();
    assert!(balance.abs()<1e-7,"{balance:e}");
}
#[test]
fn every_loaded_substep_is_observed_once_and_failed_samples_restore_air_too() {
    let mut a=piano();let mut b=piano();let model=load(&a);
    a.configure_radiation(&model).unwrap();b.configure_radiation(&model).unwrap();
    a.note_on(69,1.).unwrap();b.note_on(69,1.).unwrap();
    let mut trace=vec![0.;a.board_trace_len()];let r=a.bank.board_count;
    for _ in 0..480 {
        a.step_with_board_trace(&mut trace).unwrap();
        for sub in 0..b.substeps {
            b.mechanics_step().unwrap();
            assert_eq!(&trace[sub*r..(sub+1)*r],&b.bank.v[b.bank.modes.len()..]);
        }
        assert_eq!(a.bank.q,b.bank.q);assert_eq!(a.bank.v,b.bank.v);
        assert_eq!(a.radiation_energy_j(),b.radiation_energy_j());
    }
    let energy=a.energy_j();let loss=a.accounting.radiation_loss_j;
    a.damper_drag_ns_m=f64::NAN;assert!(a.step_with_board_trace(&mut trace).is_err());
    assert_eq!(a.energy_j(),energy);assert_eq!(a.bank.q,b.bank.q);assert_eq!(a.bank.v,b.bank.v);
    assert_eq!(a.radiation_energy_j(),b.radiation_energy_j());
    assert_eq!(a.accounting.radiation_loss_j,loss);a.damper_drag_ns_m=0.4;
    a.step().unwrap();b.step().unwrap();
    assert_eq!(a.bank.q,b.bank.q);assert_eq!(a.bank.v,b.bank.v);
    assert_eq!(a.radiation_energy_j(),b.radiation_energy_j());
}
#[test]
fn no_load_and_zero_coupling_keep_the_original_trajectory_and_admission_is_atomic() {
    let mut a=piano();let mut b=piano();let mut model=load(&a);
    model.poles[0].coupling.fill(0.);a.configure_radiation(&model).unwrap();
    assert!(a.configure_radiation(&model).is_err());
    a.note_on(69,0.5).unwrap();b.note_on(69,0.5).unwrap();
    assert!(b.configure_radiation(&model).is_err());assert!(!b.has_radiation());
    for _ in 0..480 {a.step().unwrap();b.step().unwrap();}
    assert_eq!(a.bank.q,b.bank.q);assert_eq!(a.bank.v,b.bank.v);
    assert_eq!(a.accounting.radiation_loss_j,0.);assert_eq!(a.radiation_energy_j(),0.);
}
