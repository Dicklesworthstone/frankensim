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

fn exterior_cavity(rate: u32, lossless: bool) -> SharedCavity {
    let mut source = HelmholtzModalCavity::new(source(101325.0).geometry,
        GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101325.0).unwrap(),
        CavityNeckLoss::ExteriorMouth).unwrap();
    // The exact conservative limiting fixture isolates radiation coupling;
    // ordinary constructions retain the source-owned Bessel neck resistance.
    if lossless { source.resistance = 0.0; }
    let surface = ModalAcousticTimeModel::try_new(rate, vec![mode(500.0)],
        ModalAcousticTimeBudget::audible_reference()).unwrap();
    let mut config = budget(); config.max_setup_terms = 4096;
    SharedCavity::new(source, vec![CavitySurface { model: surface,
        outward_volume_shapes: vec![if lossless {0.0} else {0.01}] }], rate,
        ModalAcousticState { displacement_m_sqrt_kg:0.0, velocity_m_sqrt_kg_per_s:0.02 },
        ModalAcousticTimeBudget::audible_reference(), config, &CancelGate::new_clock_free()).unwrap()
}

#[test]
fn shared_cavity_exterior_resonances_match_lossless_coupled_roots() {
    use crate::render::plate::impact::radiation::Pole;
    let gate = CancelGate::new_clock_free();
    for rate in [48_000, 96_000] {
        let cavity = exterior_cavity(rate, true);
        let w0 = cavity.source().isolated_omega_rad_s();
        let projection = cavity.source().neck_volume_shape_m2_per_sqrt_kg();
        let (air, coupling, velocity) = (700.0_f64, 300.0_f64, 0.02_f64);
        let load = RadiationModel { ports:1, poles:vec![Pole { omega:air, zeta:0.0,
            coupling:vec![coupling/projection] }] };
        let mut cavity = cavity.with_exterior_radiation(&load, &[1], 1, &gate).unwrap();
        // Independent conservative determinant:
        // (w0^2-x)(wa^2-x)-g^2*x=0. Both roots move; adding only a
        // prescribed neck pressure cannot produce this resonance splitting.
        let trace = w0*w0 + air*air + coupling*coupling;
        let gap = det::sqrt(trace*trace - 4.0*w0*w0*air*air);
        let high2 = 0.5*(trace+gap);
        let low2 = (w0*w0)*(air*air)/high2;
        let (low, high) = (det::sqrt(low2), det::sqrt(high2));
        let a = velocity*(high2-w0*w0-coupling*coupling)/gap;
        let b = velocity-a;
        for step in 1..=384 {
            let frame = cavity.step(&[0.0],0.0,&gate).unwrap().clone();
            let t = f64::from(step)/f64::from(rate);
            let expected_q = a/low*det::sin(low*t) + b/high*det::sin(high*t);
            let expected_v = a*det::cos(low*t) + b*det::cos(high*t);
            let air_v = coupling*velocity/gap * (high*det::sin(high*t)-low*det::sin(low*t));
            let neck = cavity.network().components()[1].states()[0];
            close(neck.displacement_m_sqrt_kg, expected_q);
            close(neck.velocity_m_sqrt_kg_per_s, expected_v);
            close(cavity.neck_exterior_pressure_pa().unwrap(), coupling/projection*air_v);
            close(frame.total_energy_j(),0.5*velocity*velocity);
            assert_eq!(frame.total_dissipation_j(),0.0);
        }
        assert!(low < w0 && high > air);
    }
}

fn mutually_loaded_cavity(rate: u32) -> SharedCavity {
    use crate::render::plate::impact::radiation::Pole;
    let load = RadiationModel { ports:2, poles:vec![Pole { omega:700.0, zeta:0.12,
        coupling:vec![6000.0,-180.0] }] };
    // The neck may be ANY model column, not just the last one. The body
    // shape's opposite sign retains its real mutual reaction.
    exterior_cavity(rate,false).with_exterior_radiation(&load,&[1,0],1,
        &CancelGate::new_clock_free()).unwrap()
}

#[test]
fn shared_cavity_mutual_radiation_preserves_power_and_exact_time_partition() {
    let gate = CancelGate::new_clock_free();
    let mut whole = mutually_loaded_cavity(48_000);
    let mut half = mutually_loaded_cavity(96_000);
    let projection = whole.source().neck_volume_shape_m2_per_sqrt_kg();
    let scaled = whole.network().radiation_model().unwrap();
    close(scaled.poles[0].coupling[0],6000.0*projection);
    assert_eq!(scaled.poles[0].coupling[1],-180.0);
    let initial = whole.network().total_energy_j().unwrap();
    let mut loss = 0.0;
    let mut work = 0.0;
    for _ in 0..64 {
        let full = whole.step(&[0.001],0.2,&gate).unwrap().clone();
        let first = half.step(&[0.001],0.2,&gate).unwrap().clone();
        let second = half.step(&[0.001],0.2,&gate).unwrap().clone();
        for (a,b) in whole.network().components().iter().zip(half.network().components()) {
            close(a.states()[0].displacement_m_sqrt_kg,b.states()[0].displacement_m_sqrt_kg);
            close(a.states()[0].velocity_m_sqrt_kg_per_s,b.states()[0].velocity_m_sqrt_kg_per_s);
        }
        for (&a,&b) in whole.network().radiation_state().iter().zip(half.network().radiation_state()) {
            close(a,b);
        }
        close(full.total_dissipation_j(),first.total_dissipation_j()+second.total_dissipation_j());
        close(full.external_work_j,first.external_work_j+second.external_work_j);
        let pressure = whole.neck_exterior_pressure_pa().unwrap();
        close(pressure,6000.0*whole.network().radiation_state()[1]);
        let flow = whole.observe().unwrap().neck_outflow_m3_s;
        let v = whole.network().components()[1].states()[0].velocity_m_sqrt_kg_per_s;
        close(pressure*flow,(6000.0*projection*whole.network().radiation_state()[1])*v);
        assert!(full.radiation_dissipation_j>0.0);
        assert!(full.energy_residual_j.abs()<=full.energy_tolerance_j);
        loss+=full.total_dissipation_j(); work+=full.external_work_j;
    }
    close(whole.network().total_energy_j().unwrap()+loss,initial+work);
}

#[test]
fn shared_cavity_exterior_load_refuses_double_counting_and_preserves_rejected_state() {
    use crate::render::plate::impact::radiation::Pole;
    let gate = CancelGate::new_clock_free();
    let load = RadiationModel {ports:1,poles:vec![Pole {omega:700.0,zeta:0.1,coupling:vec![6000.0]}]};
    assert!(shared(101325.0).with_exterior_radiation(&load,&[2],1,&gate).is_err());
    assert!(exterior_cavity(48000,false).with_exterior_radiation(&load,&[0],1,&gate).is_err());
    assert!(exterior_cavity(48000,false).with_exterior_radiation(&load,&[1],0,&gate).is_err());
    let duplicate = RadiationModel {ports:2,poles:vec![Pole {omega:700.0,zeta:0.1,coupling:vec![1.0,1.0]}]};
    assert!(exterior_cavity(48000,false).with_exterior_radiation(&duplicate,&[1,1],1,&gate).is_err());
    assert!(mutually_loaded_cavity(48000).with_exterior_radiation(&load,&[1],1,&gate).is_err());
    let mut late = exterior_cavity(48000,false);
    late.step(&[0.0],0.0,&gate).unwrap();
    assert!(late.with_exterior_radiation(&load,&[1],1,&gate).is_err());

    let mut cavity = mutually_loaded_cavity(48000);
    cavity.step(&[0.001],0.2,&gate).unwrap();
    let states: Vec<_> = cavity.network().components().iter().map(|m|m.states().to_vec()).collect();
    let acoustic = cavity.network().radiation_state().to_vec();
    let pressure = cavity.neck_exterior_pressure_pa().unwrap();
    let frame = cavity.network().last_frame().unwrap().clone();
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(matches!(cavity.step(&[0.0],0.0,&cancelled),Err(ModalCouplingError::Cancelled)));
    assert!(cavity.step(&[1e20],0.0,&gate).is_err());
    for (m,expected) in cavity.network().components().iter().zip(states) {assert_eq!(m.states(),expected);}
    assert_eq!(cavity.network().radiation_state(),acoustic);
    assert_eq!(cavity.neck_exterior_pressure_pa().unwrap(),pressure);
    assert_eq!(cavity.network().last_frame().unwrap(),&frame);
}
