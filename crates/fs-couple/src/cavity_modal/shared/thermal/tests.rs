use super::*;
use fs_material::gas::GasSpec;

fn source(loss:CavityNeckLoss)->HelmholtzModalCavity {
    HelmholtzModalCavity::new(HelmholtzCavity {volume_m3:1e-7,neck_radius_m:1e-4,
        neck_length_m:0.1,inner_end:fs_scenario::NeckEnd::Flanged,
        outer_end:fs_scenario::NeckEnd::Flanged},
        GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101325.0).unwrap(),loss).unwrap()
}
fn options()->CylindricalCavityThermalOptions {
    CylindricalCavityThermalOptions {radius_m:0.001,poles:4,band_rad_s:[10.0,1000.0],
        maximum_relative_error:1e-3}
}
fn budget()->ModalCouplingConfig {
    ModalCouplingConfig {max_modes:8,max_connections:6,max_setup_terms:20000,
        nyquist_guard_fraction:0.9,maximum_total_energy_j:1.0,maximum_abs_pressure_pa:1e5,
        maximum_abs_connection_force_n:100.0,solve_relative_tolerance:1e-10,
        energy_absolute_tolerance_j:1e-14,energy_relative_tolerance:1e-8}
}
fn surface(rate:u32)->CavitySurface {
    CavitySurface {model:ModalAcousticTimeModel::try_new(rate,vec![ModalAcousticMode {
        angular_frequency_rad_s:500.0,damping_ratio:0.0,pressure_per_modal_velocity:C64::ZERO}],
        ModalAcousticTimeBudget::audible_reference()).unwrap(),outward_volume_shapes:vec![1e-4]}
}
fn create(rate:u32,viscous:bool,history:Option<Vec<f64>>)->SharedCavity {
    let neck=viscous.then_some(ViscousNeckOptions {poles:4,band_rad_s:[10.0,1000.0],maximum_relative_error:1e-3});
    SharedCavity::new_with_cylindrical_thermal_wall(source(if viscous {CavityNeckLoss::ExteriorMouth}
        else {CavityNeckLoss::Resistance(0.0)}),vec![surface(rate)],rate,
        ModalAcousticState {displacement_m_sqrt_kg:2e-7,velocity_m_sqrt_kg_per_s:0.0},
        ModalAcousticTimeBudget::audible_reference(),budget(),neck,options(),history,
        &CancelGate::new_clock_free()).unwrap()
}
fn close(a:f64,b:f64,absolute:f64) {
    assert!((a-b).abs()<=absolute+2e-10*a.abs().max(b.abs()),"{a:e} != {b:e}");
}

#[test]
fn cylindrical_thermal_compression_converges_to_independent_bessel_compliance() {
    // G2/G3: an independent duct Bessel shunt, three diffusion resolutions,
    // and both physical storage limits. Viscous series impedance cancels in k/Zc.
    let source=source(CavityNeckLoss::Resistance(0.0));let gas=source.gas();
    let gate=CancelGate::new_clock_free();let mut previous=f64::INFINITY;
    for poles in [2,4,8] {
        let law=CylindricalThermalCompression::new(&source,CylindricalCavityThermalOptions {
            poles,band_rad_s:[10.0,5000.0],maximum_relative_error:0.95,..options()},&gate).unwrap();
        close(law.equilibrium_stiffness_n_m(),law.adiabatic_stiffness_n_m()/gas.gamma,1e-18);
        let instantaneous=law.equilibrium_stiffness_n_m()+law.arms().iter().map(|a|a.stiffness_n_m).sum::<f64>();
        close(instantaneous,law.adiabatic_stiffness_n_m(),1e-18);
        assert!(law.arms().iter().all(|a|a.stiffness_n_m>0.0 && a.relaxation_time_s>0.0));
        close(core::f64::consts::PI*law.options().radius_m.powi(2)*law.length_m(),source.geometry().volume_m3,1e-22);
        let area=core::f64::consts::PI*law.options().radius_m.powi(2);
        let mut error=0.0_f64;
        for i in 0..=64 {
            let omega=if i==64 {5000.0}else{10.0*det::pow(500.0,i as f64/64.0)};
            let wave=fs_duct::segment_wave(gas,law.options().radius_m,omega,LossModel::Bessel).unwrap();
            let c_minus=(wave.wavenumber/wave.characteristic_impedance)
                .scale(gas.density*gas.sound_speed.powi(2)/(omega*area));
            let exact=C64::from_re(law.adiabatic_stiffness_n_m())/c_minus.conj();
            let actual=law.stiffness_at(omega).unwrap();
            assert!(actual.re>0.0 && actual.im>0.0);
            error=error.max((actual-exact).abs()/exact.abs());
            let compliance_inverse=C64::new(0.0,omega)
                /law.diffusion().response_at(omega).unwrap().1;
            close((actual.scale(law.adiabatic_stiffness_n_m().recip())-compliance_inverse).abs(),0.0,2e-11);
        }
        assert!(error<previous,"poles={poles}, error={error:e}, previous={previous:e}");
        previous=error;
    }
    assert!(previous<1e-5);
}

#[test]
fn thermal_history_and_runtime_match_physical_pressure_diffusion_and_availability() {
    // G1: pressure plus original thermal diffusion coordinates, integrated by
    // refined RK4, are independent of the production Maxwell transformation.
    let gate=CancelGate::new_clock_free();
    let relaxed=create(48000,false,None);
    let mut actual=create(48000,false,Some(vec![0.0;4]));
    close(actual.observe().unwrap().pressure_increment_pa,
        relaxed.observe().unwrap().pressure_increment_pa*actual.source().gas().gamma,1e-10);
    close(actual.observe().unwrap().compression_energy_j,
        relaxed.observe().unwrap().compression_energy_j*actual.source().gas().gamma,1e-20);
    let source=actual.source().clone();let area=source.neck_area_m2();
    let compliance=area*area/source.compression_stiffness_n_m();
    let inertia=source.neck_mass_kg()/(area*area);
    let terms=actual.cavity_thermal_compression().unwrap().diffusion().shunt_terms().to_vec();
    let mut x=vec![0.0;4+terms.len()]; // volume, neck flow, surface q/v, thermal pressure histories
    x[0]=actual.observe().unwrap().volume_change_m3;
    let pressure=|x:&[f64]|-x[0]/compliance-terms.iter().enumerate()
        .map(|(i,t)|t.gain_per_s/t.rate_per_s*x[i+4]).sum::<f64>();
    let derivative=|x:&[f64],external:f64| {
        let p=pressure(x);let mut d=vec![0.0;x.len()];
        d[0]=x[1]+1e-4*x[3];d[1]=(p-external)/inertia;
        d[2]=x[3];d[3]=-500.0_f64.powi(2)*x[2]+1e-4*p;
        for (i,t) in terms.iter().enumerate() {d[i+4]=t.rate_per_s*(p-x[i+4]);}d
    };
    let initial=actual.network().total_energy_j().unwrap();let (mut work,mut loss)=(0.0,0.0);
    for tick in 0..128 {
        let external=if tick<16 {0.1}else{0.0};let dt=1.0/(48000.0*16.0);
        for _ in 0..16 {
            let stage=|k:&[f64],a:f64|x.iter().zip(k).map(|(x,k)|x+a*dt*k).collect::<Vec<_>>();
            let a=derivative(&x,external);let b=derivative(&stage(&a,0.5),external);
            let c=derivative(&stage(&b,0.5),external);let d=derivative(&stage(&c,1.0),external);
            for i in 0..x.len() {x[i]+=dt*(a[i]+2.0*b[i]+2.0*c[i]+d[i])/6.0;}
        }
        let frame=actual.step(&[0.0],external,&gate).unwrap().clone();
        let observed=actual.observe().unwrap();let p=pressure(&x);
        close(observed.pressure_increment_pa,p,1e-8);
        close(observed.neck_outflow_m3_s,x[1],1e-18);
        close(observed.volume_change_m3,x[0],1e-21);
        close(actual.network().components()[0].states()[0].displacement_m_sqrt_kg,x[2],1e-17);
        close(actual.network().components()[0].states()[0].velocity_m_sqrt_kg_per_s,x[3],1e-14);
        let availability=0.5*compliance*(p*p+terms.iter().enumerate()
            .map(|(i,t)|t.gain_per_s/t.rate_per_s*x[i+4]*x[i+4]).sum::<f64>());
        close(observed.compression_energy_j,availability,1e-19);
        let surface_energy=0.5*(x[3]*x[3]+500.0_f64.powi(2)*x[2]*x[2]);
        close(frame.total_energy_j(),availability+0.5*inertia*x[1]*x[1]+surface_energy,1e-19);
        assert!(frame.material_dissipation_j>=0.0);
        work+=frame.external_work_j;loss+=frame.total_dissipation_j();
    }
    close(actual.network().total_energy_j().unwrap()+loss,initial+work,1e-19);
    assert!(loss>0.001*initial);
}

#[test]
fn thermal_compression_composes_with_viscosity_radiation_and_atomic_rejection() {
    use crate::render::plate::impact::radiation::Pole;
    let gate=CancelGate::new_clock_free();
    let radiation=RadiationModel {ports:1,poles:vec![Pole {omega:700.0,zeta:0.1,coupling:vec![2e6]}]};
    let build=|rate|create(rate,true,Some(vec![0.0;4]))
        .with_exterior_radiation(&radiation,&[1],1,&gate).unwrap();
    let mut whole=build(48000);let mut half=build(96000);
    for _ in 0..12 {
        let full=whole.step(&[0.0],0.03,&gate).unwrap().clone();
        let first=half.step(&[0.0],0.03,&gate).unwrap().clone();
        let second=half.step(&[0.0],0.03,&gate).unwrap().clone();
        close(whole.observe().unwrap().pressure_increment_pa,half.observe().unwrap().pressure_increment_pa,1e-8);
        for (a,b) in whole.network().material_memory().iter().zip(half.network().material_memory()) {
            close(*a,*b,1e-15);
        }
        close(full.total_dissipation_j(),first.total_dissipation_j()+second.total_dissipation_j(),1e-18);
        assert_eq!(whole.ports.schemas.len(),2);
        assert_eq!(whole.ports.transfers().unwrap().len(),2);
    }
    let state:Vec<_>=whole.network().components().iter().map(|m|m.states().to_vec()).collect();
    let thermal=whole.network().material_memory().to_vec();
    let acoustic=whole.network().radiation_state().to_vec();
    let observation=whole.observe().unwrap();let frame=whole.network().last_frame().unwrap().clone();
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(matches!(whole.step(&[0.0],0.0,&cancelled),Err(ModalCouplingError::Cancelled)));
    assert!(whole.step(&[0.0],1e30,&gate).is_err());
    assert_eq!(whole.network().material_memory(),thermal);
    assert_eq!(whole.network().radiation_state(),acoustic);
    assert_eq!(whole.observe().unwrap(),observation);
    assert_eq!(whole.network().last_frame().unwrap(),&frame);
    for (m,expected) in whole.network().components().iter().zip(state) {assert_eq!(m.states(),expected);}
    let attempt=|config,options,history,gate:&CancelGate|SharedCavity::new_with_cylindrical_thermal_wall(
        source(CavityNeckLoss::Resistance(0.0)),vec![surface(48000)],48000,ModalAcousticState::default(),
        ModalAcousticTimeBudget::audible_reference(),config,None,options,history,gate);
    assert!(attempt(budget(),CylindricalCavityThermalOptions {radius_m:0.0,..options()},None,&gate).is_err());
    assert!(attempt(budget(),options(),Some(vec![0.0;2]),&gate).is_err());
    assert!(attempt(ModalCouplingConfig {max_setup_terms:20,..budget()},options(),None,&gate).is_err());
    assert!(attempt(budget(),CylindricalCavityThermalOptions {band_rad_s:[10.0,1e6],..options()},None,&gate).is_err());
    assert!(matches!(attempt(budget(),options(),None,&cancelled),Err(ModalCouplingError::Cancelled)));
}
