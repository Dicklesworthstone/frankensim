use super::*;
use fs_couple::render::plate::impact::ImpactSubstepConfig;
const GAS:&str="frankensim-squeeze-film-v1\nannulus,0.04,0.085,2,4\nviscosity_pa_s,0.000018\nlimits,0.000001,0.04,200000\nboundary,sealed,sealed\nisothermal,100000,293.15,287.05\n";
fn spec()->Spec {Spec::parse(include_str!("estimated-hihat.fshh")).unwrap()}
fn shell()->specimen::Specimen {let mut s=specimen::Specimen::reference();s.azimuths=8;s}
fn stroke()->Stroke {Stroke{speed_m_s:0.0,position_m:None}}

#[test]
fn explicit_gas_selects_storage_and_admits_seals_without_weakening_old_inputs() {
    let c=squeeze::Config::parse(GAS).unwrap();assert!(c.gas().is_some());
    let plain=GAS.replace("isothermal,100000,293.15,287.05\n","");
    assert!(squeeze::Config::parse(&plain).is_err());
    assert!(squeeze::Config::parse(&plain.replace("sealed,sealed","open,open")).unwrap().gas().is_none());
    for bad in [GAS.replace("isothermal,100000","isothermal,0"),
        GAS.replace("293.15","-1"),GAS.replace("287.05","NaN"),
        format!("{GAS}isothermal,100000,293.15,287.05"),
        GAS.replace("isothermal,100000,293.15,287.05","isothermal,100000,293.15"),
        GAS.replace("viscosity_pa_s,0.000018","viscosity_pa_s,0")]
    {assert!(squeeze::Config::parse(&bad).is_err(),"{bad}");}
}

#[test]
fn actual_paired_shells_add_gas_after_histories_without_phantom_source_modes() {
    let c=squeeze::Config::parse(GAS).unwrap();let a=shell();let b=shell();
    let pair=build_with_squeeze(&spec(),&a,&b,stroke(),None,128,2e-6,true,Some(&c)).unwrap();
    let e=pair.experiment;let n=e.force.len();let o=play::gas_observation(&e.system).unwrap().unwrap();
    assert_eq!(o.free_energy_j,0.0);assert_eq!(o.minimum_pressure_pa,100000.0);assert!(o.mass_kg>0.0);
    assert_eq!(e.system.state().len(),2*n+12+8);
    let sources=e.acoustics.as_ref().unwrap().state_modes();
    assert_eq!(sources,pair.upper_modes.clone().chain(pair.lower_modes.clone()).collect::<Vec<_>>());
    assert!(!sources.contains(&pair.pedal.coordinate));assert!(!sources.contains(&0));
    if let Mechanics::Reference(system)=&e.system {
        for i in 0..12 {assert!(system.felt_history(i).is_some());}
    }else{panic!("reference construction must not choose another time image");}
    let plain=build(&spec(),&a,&b,stroke(),None,1,2e-6,false).unwrap();
    assert!(play::gas_observation(&plain.experiment.system).unwrap().is_none());
    assert_eq!(&e.system.state()[..2*n+12],plain.experiment.system.state());
}

#[test]
fn compressed_real_skins_retain_mass_and_publish_one_joint_energy_history() {
    let c=squeeze::Config::parse(GAS).unwrap();let a=shell();let b=shell();
    let pair=build_with_squeeze(&spec(),&a,&b,stroke(),None,128,2e-6,false,Some(&c)).unwrap();
    let mut e=pair.experiment;let n=e.force.len();
    let initial_mass=play::gas_observation(&e.system).unwrap().unwrap().mass_kg;
    e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
        ImpactSubstepConfig{max_depth:8,max_attempts:511}).unwrap();
    // Explicit downward force on the upper shell's retained rigid translation.
    // No direct gas, lower-shell, pedal or stick force is invented.
    let mut force=vec![0.0;n];force[pair.upper_modes.start]=-5.0;
    let gate=CancelGate::new_clock_free();let (mut work,mut loss,mut initial)=(0.0,0.0,None);
    let mut peak=0.0_f64;let mut last=0.0;
    for i in 0..128 {
        if i==32 {let x=e.system.state().to_vec();let cancelled=CancelGate::new_clock_free();cancelled.request();
            assert!(e.system.step(&force,&cancelled).is_err());assert_eq!(e.system.state(),x);}
        let f=e.system.step(&force,&gate).unwrap();let o=play::gas_observation(&e.system).unwrap().unwrap();
        if initial.is_none(){initial=Some(f.stored_energy_j+f.dissipated_energy_j-f.supplied_work_j);}
        work+=f.supplied_work_j;loss+=f.dissipated_energy_j;last=f.stored_energy_j;
        peak=peak.max(o.maximum_pressure_pa);assert!(o.free_energy_j>=0.0);
        assert!((o.mass_kg-initial_mass).abs()<1e-8*initial_mass);
    }
    assert!(peak>100000.001,"actual skin closure must compress retained gas");
    assert!((last+loss-initial.unwrap()-work).abs()<1e-6);
}

#[test]
fn gas_observation_follows_the_existing_player_wrapper_without_a_second_clock() {
    let c=squeeze::Config::parse(GAS).unwrap();let a=shell();let b=shell();
    let pair=build_with_squeeze(&spec(),&a,&b,stroke(),None,4,2e-6,false,Some(&c)).unwrap();
    let e=pair.experiment;let before=play::gas_observation(&e.system).unwrap().unwrap();
    let driven=e.system.with_stick_drive(mechanics::drive::Program::parse("0,0\n0.000008,0").unwrap(),
        2e-6,4,e.stick_weight,e.force.len()).unwrap();
    let after=play::gas_observation(&driven).unwrap().unwrap();
    assert_eq!(before.mass_kg,after.mass_kg);assert_eq!(before.free_energy_j,after.free_energy_j);
}
