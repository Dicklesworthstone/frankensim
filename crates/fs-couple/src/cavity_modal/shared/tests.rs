use super::*;
use fs_material::gas::GasSpec;

fn source(p:f64)->HelmholtzModalCavity {
    HelmholtzModalCavity::new(HelmholtzCavity {volume_m3:0.01,neck_radius_m:0.01,neck_length_m:0.02,
        inner_end:fs_scenario::NeckEnd::Flanged,outer_end:fs_scenario::NeckEnd::Flanged},
        GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,p).unwrap(),CavityNeckLoss::Resistance(2000.0)).unwrap()
}
fn budget()->ModalCouplingConfig {
    ModalCouplingConfig {max_modes:8,max_connections:4,max_setup_terms:1024,nyquist_guard_fraction:0.9,
        maximum_total_energy_j:1.0,maximum_abs_pressure_pa:1e5,maximum_abs_connection_force_n:100.0,
        solve_relative_tolerance:1e-10,energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-8}
}
fn mode(w:f64)->ModalAcousticMode {ModalAcousticMode {angular_frequency_rad_s:w,damping_ratio:0.0,pressure_per_modal_velocity:C64::ZERO}}
fn model(m:Vec<ModalAcousticMode>)->ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(48000,m,ModalAcousticTimeBudget::audible_reference()).unwrap()
}
fn shared(p:f64)->SharedCavity {
    SharedCavity::new(source(p),vec![
        CavitySurface {model:model(vec![mode(500.0)]),outward_volume_shapes:vec![0.01]},
        CavitySurface {model:model(vec![mode(900.0)]),outward_volume_shapes:vec![-0.02]},
    ],48000,ModalAcousticState::default(),ModalAcousticTimeBudget::audible_reference(),budget(),&CancelGate::new_clock_free()).unwrap()
}
fn close(a:f64,b:f64) {assert!((a-b).abs()<1e-13+1e-11*a.abs().max(b.abs()),"{a:e} != {b:e}");}

#[test]
fn multiple_surfaces_match_the_original_single_combined_basis() {
    let mut a=shared(101325.0);let c=source(101325.0);let gate=CancelGate::new_clock_free();
    let mut reference=CoupledModalSystem::new(vec![model(vec![mode(500.0),mode(900.0)]),
        c.neck_component(48000,ModalAcousticTimeBudget::audible_reference()).unwrap()],
        c.connections(0,1,&[0.01,-0.02]).unwrap(),budget(),&gate).unwrap();
    for i in 0..256 {
        let surface=[if i<20 {0.001}else{0.0},if i<30 {-0.002}else{0.0}];
        let exterior=if i<12 {0.5}else{0.0};
        let actual=a.step(&surface,exterior,&gate).unwrap().clone();
        let expected=reference.step(&[surface[0],surface[1],-c.neck_volume_shape_m2_per_sqrt_kg()*exterior]).unwrap();
        close(actual.external_work_j,expected.external_work_j);
        close(actual.connection_energy_j,expected.connection_energy_j);
        for j in 0..2 {
            let x=a.network().components()[j].states()[0];let y=reference.components()[0].states()[j];
            close(x.displacement_m_sqrt_kg,y.displacement_m_sqrt_kg);
            close(x.velocity_m_sqrt_kg_per_s,y.velocity_m_sqrt_kg_per_s);
        }
        let observed=a.observe().unwrap();
        close(observed.compression_energy_j,actual.connection_energy_j);
        close(observed.pressure_increment_pa,
            -c.gas().density*c.gas().sound_speed.powi(2)/c.geometry().volume_m3*observed.volume_change_m3);
    }
    assert!(a.network().components()[1].states()[0].displacement_m_sqrt_kg.abs()>1e-12);
}

#[test]
fn signed_volume_cancellation_does_not_create_per_surface_pressure_springs() {
    let mut m=model(vec![mode(700.0)]);
    m.restore_states(&[ModalAcousticState {displacement_m_sqrt_kg:1e-4,velocity_m_sqrt_kg_per_s:0.0}]).unwrap();
    let gate=CancelGate::new_clock_free();
    let mut a=SharedCavity::new(source(101325.0),vec![
        CavitySurface {model:m.clone(),outward_volume_shapes:vec![0.01]},
        CavitySurface {model:m,outward_volume_shapes:vec![-0.01]},
    ],48000,ModalAcousticState::default(),ModalAcousticTimeBudget::audible_reference(),budget(),&gate).unwrap();
    for _ in 0..100 {
        a.step(&[0.0;2],0.0,&gate).unwrap();
        let observed=a.observe().unwrap();
        assert_eq!(observed.pressure_increment_pa,0.0);
        assert_eq!(observed.compression_energy_j,0.0);
        assert_eq!(a.network().components()[2].states()[0],ModalAcousticState::default());
    }
}

#[test]
fn gas_data_changes_surface_feedback_and_energy_closes_without_extra_heat() {
    let mut a=shared(50000.0);let mut b=shared(150000.0);let gate=CancelGate::new_clock_free();
    assert_eq!(a.network().components()[0].modes(),b.network().components()[0].modes());
    let mut input=[0.0;2];let mut loss=[0.0;2];
    for i in 0..512 {
        let f=[if i<30 {0.004}else{0.0},0.0];
        for (j,s) in [&mut a,&mut b].into_iter().enumerate() {
            let frame=s.step(&f,0.0,&gate).unwrap();
            input[j]+=frame.external_work_j;
            loss[j]+=frame.component_dissipation_j+frame.connection_dissipation_j;
            assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
        }
    }
    for (j,s) in [&a,&b].into_iter().enumerate() {
        close(s.network().total_energy_j().unwrap()+loss[j],input[j]);
        assert!(loss[j]>0.0);
    }
    let qa=a.network().components()[0].states()[0].displacement_m_sqrt_kg;
    let qb=b.network().components()[0].states()[0].displacement_m_sqrt_kg;
    assert!((qa-qb).abs()>1e-12);
}

#[test]
fn bad_projection_pressure_and_cancellation_preserve_all_accepted_components() {
    let gate=CancelGate::new_clock_free();let mut a=shared(101325.0);
    a.step(&[0.001,0.0],0.0,&gate).unwrap();
    let states:Vec<_>=a.network().components().iter().map(|m|m.states().to_vec()).collect();
    let frame=a.network().last_frame().unwrap().clone();let observation=a.observe().unwrap();
    assert!(a.step(&[0.0],0.0,&gate).is_err());
    assert!(a.step(&[0.0;2],f64::NAN,&gate).is_err());
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(matches!(a.step(&[0.0;2],0.0,&cancelled),Err(ModalCouplingError::Cancelled)));
    // Budget failure after a genuine attempted coupled solve also rolls back.
    assert!(a.step(&[1e20,0.0],0.0,&gate).is_err());
    assert_eq!(a.network().last_frame().unwrap(),&frame);
    assert_eq!(a.observe().unwrap(),observation);
    for (m,expected) in a.network().components().iter().zip(states) {assert_eq!(m.states(),expected);}
    assert!(SharedCavity::new(source(101325.0),vec![CavitySurface {model:model(vec![mode(500.0)]),
        outward_volume_shapes:vec![]}],48000,ModalAcousticState::default(),
        ModalAcousticTimeBudget::audible_reference(),budget(),&gate).is_err());
}
