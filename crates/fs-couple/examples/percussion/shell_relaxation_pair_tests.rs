use super::*;
use fs_couple::render::plate::impact::ImpactSubstepConfig;

const CARD:&str="frankensim-shell-relaxation-v1\nintrinsic_loss,replace\ninitial,relaxed\nband_hz,0,5000\nshell,lower\nbranch,lower,0.05,0.0002\n";
const GAS:&str="frankensim-squeeze-film-v1\nannulus,0.04,0.085,2,4\nviscosity_pa_s,0.000018\nlimits,0.000001,0.04,200000\nboundary,sealed,sealed\nisothermal,100000,293.15,287.05\n";
fn spec()->Spec {Spec::parse(include_str!("estimated-hihat.fshh")).unwrap()}
fn shell()->specimen::Specimen {let mut s=specimen::Specimen::reference();s.azimuths=8;s}
fn build(material:Option<&shell_relaxation::Spec>,film:bool)->Pair {
    let mut spec=spec();spec.damping[1]=0.;let s=shell();
    let gas=squeeze::Config::parse(GAS).unwrap();
    build_with_material(&spec,&s,&s,Stroke{speed_m_s:0.,position_m:None},None,128,2e-6,true,
        film.then_some(&gas),&shaft_playing::Selection::default(),&mallets::Selection::default(),material).unwrap()
}

#[test]
fn selecting_lower_material_does_not_remove_the_upper_shells_declared_damping() {
    let elastic=shell_relaxation::Spec::read(&CARD.replace("branch,lower,0.05,0.0002\n",""),true).unwrap();
    let material=shell_relaxation::Spec::read(CARD,true).unwrap();let s=shell();
    assert!(!material.selected(0));assert!(material.selected(1));
    assert!(build_with_material(&spec(),&s,&s,Stroke::default(),None,1,2e-6,false,None,
        &shaft_playing::Selection::default(),&mallets::Selection::default(),Some(&material)).is_err());
    let mut original=build(None,false);let mut selected=build(Some(&elastic),false);
    let n=original.experiment.force.len();let mut force=vec![0.;n];force[original.upper_modes.start+1]=0.1;
    original.experiment.system=original.experiment.system.into_analytic_nonlinear().unwrap();
    selected.experiment.system=selected.experiment.system.into_analytic_nonlinear().unwrap();
    let gate=CancelGate::new_clock_free();
    for _ in 0..32 {
        original.experiment.system.step(&force,&gate).unwrap();selected.experiment.system.step(&force,&gate).unwrap();
        assert_eq!(original.experiment.system.state(),selected.experiment.system.state());
    }
    assert_eq!(head_relaxation::observation(&selected.experiment.system).states,0);
}

#[test]
fn paired_shell_memory_follows_gas_and_mounts_without_moving_sources_or_player_coordinates() {
    let material=shell_relaxation::Spec::read(CARD,true).unwrap();
    let pair=build(Some(&material),true);let mut e=pair.experiment;let n=e.force.len();
    let modes=pair.lower_modes.len()-1;
    assert_eq!(head_relaxation::observation(&e.system).states,modes);
    assert_eq!(e.system.state().len(),2*n+12+8+modes);
    assert_eq!(e.acoustics.as_ref().unwrap().state_modes(),
        pair.upper_modes.clone().chain(pair.lower_modes.clone()).collect::<Vec<_>>());
    assert_eq!(pair.pedal.coordinate,pair.lower_modes.end);
    let initial_mass=play::gas_observation(&e.system).unwrap().unwrap().mass_kg;
    let initial=match &e.system {Mechanics::Reference(s)=>s.stored_energy_j(),_=>panic!()};
    e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
        ImpactSubstepConfig{max_depth:8,max_attempts:511}).unwrap();
    // Project a 0.1 N force at the actual off-axis strike station onto upper
    // flexure, in addition to rigid closure. This does not pick an arbitrary
    // member/orientation of a near-degenerate angular mode pair.
    let mut force=vec![0.;n];force[pair.upper_modes.start]=-5.;
    for i in pair.upper_modes.start+1..pair.upper_modes.end {force[i]=0.1*e.observer_a[i];}
    assert!(force[pair.lower_modes.clone()].iter().all(|v|*v==0.));
    let gate=CancelGate::new_clock_free();let (mut net,mut last)=(0.,initial);
    let (mut peak,mut pressure_spread)=(0.0_f64,0.0_f64);
    for tick in 0..128 {
        if tick==32 {
            let state=e.system.state().to_vec();let memory=head_relaxation::observation(&e.system);
            let cancel=CancelGate::new_clock_free();cancel.request();assert!(e.system.step(&force,&cancel).is_err());
            assert_eq!(e.system.state(),state);assert_eq!(head_relaxation::observation(&e.system),memory);
        }
        let frame=e.system.step(&force,&gate).unwrap();net+=frame.supplied_work_j-frame.dissipated_energy_j;
        last=frame.stored_energy_j;peak=peak.max(head_relaxation::observation(&e.system).stored_energy_j);
        let gas=play::gas_observation(&e.system).unwrap().unwrap();assert!((gas.mass_kg-initial_mass).abs()<1e-8*initial_mass);
        pressure_spread=pressure_spread.max(gas.maximum_pressure_pa-gas.minimum_pressure_pa);
    }
    assert!(pressure_spread>1e-10*100000.,"gap pressure variation must exceed roundoff of the supplied atmospheric pressure");
    assert!(peak>0.,"pressure on the unforced lower shell must excite its flexural memory");
    assert!((last-initial-net).abs()<1e-6);
}
