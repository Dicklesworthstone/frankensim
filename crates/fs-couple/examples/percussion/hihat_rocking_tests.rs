use super::*;
use fs_couple::render::plate::impact::ImpactSubstepConfig;
use fs_couple::render::plate::impact::striker::RadiusStation;
use fs_plate::shell::stiffened::beam::RoundBeamSpec;

const GAS:&str="frankensim-squeeze-film-v1\nannulus,0.04,0.085,2,4\nviscosity_pa_s,0.000018\nlimits,0.000001,0.04,200000\nboundary,sealed,sealed\nisothermal,100000,293.15,287.05\n";
fn inner(m:&Mechanics)->&ImpactSystem {
    match m {Mechanics::Reference(s)=>s,Mechanics::Nonlinear(s)=>s,Mechanics::Substepped(s)=>s,
        Mechanics::Driven{inner:s,..}=>inner(s),Mechanics::Prepared(_)=>panic!("paired felt history is required")}
}
fn shaft()->FlexibleStriker {
    FlexibleStriker::new(&[(0.0,0.005),(0.4,0.005)].map(|(position_m,radius_m)|
        RadiusStation{position_m,radius_m}),RoundBeamSpec{young_pa:12e9,density_kg_m3:800.,
        pivot_m:0.1,contact_m:0.39,hand_m:0.16,subdivisions:8,maximum_hz:3000.,maximum_modes:17},0.001).unwrap()
}

#[test]
fn paired_rocking_keeps_both_targets_all_mounts_gap_air_and_the_shared_player_clock() {
    let dt=2e-6;let steps=768;let cap=0.03;
    let make=|rocking| {
        let mut s=specimen::Specimen::reference();s.azimuths=8;
        let spec=Spec::parse(include_str!("estimated-hihat.fshh")).unwrap();
        let film=squeeze::Config::parse(GAS).unwrap();
        let shafts=shaft_playing::Selection{first:Some(shaft()),second:Some(shaft())};
        build_with_rocking(&spec,&s,&s,Stroke{position_m:Some([0.06,0.01]),speed_m_s:0.4},
            Some(Stroke{position_m:Some([-0.06,0.02]),speed_m_s:0.35}),steps,dt,true,Some(&film),
            &shafts,&mallets::Selection::default(),None,None,
            targets::Selection{first:targets::Target::Upper,second:targets::Target::Lower},rocking).unwrap()
    };
    let fixed=make(None);let pair=make(Some(cap));
    assert_eq!(pair.upper_modes.len(),fixed.upper_modes.len()+2);
    assert_eq!(pair.lower_modes.len(),fixed.lower_modes.len()+2);
    assert_eq!(pair.pedal.coordinate,fixed.pedal.coordinate+4);
    assert_eq!(pair.experiment.force.len(),fixed.experiment.force.len()+4);
    assert!(rocking::observation(&fixed.experiment.system,1).is_none());
    let mut e=pair.experiment;let n=e.force.len();
    assert_eq!(e.system.state().len(),2*n+12+8);
    for i in 0..12 {assert!(inner(&e.system).felt_history(i).is_some());}
    assert!(inner(&e.system).felt_history(12).is_none());
    for body in [1,2] {assert_eq!(rocking::observation(&e.system,body).unwrap(),[0.;3]);}
    assert_eq!(e.acoustics.as_ref().unwrap().state_modes(),
        pair.upper_modes.clone().chain(pair.lower_modes.clone()).collect::<Vec<_>>());
    assert_eq!(e.second_stick.unwrap().coordinate,pair.pedal.coordinate+1);
    assert_eq!(e.flexible_sticks[0].as_ref().unwrap().elastic_start(),pair.pedal.coordinate+2);
    assert_eq!(pair.collision.n_points(),fixed.collision.n_points());
    for row in pair.collision.collocation().chunks(n) {
        assert_eq!(row[0],0.);assert_eq!(row[pair.pedal.coordinate],0.);
        assert_eq!(row[e.second_stick.unwrap().coordinate],0.);
        for start in [pair.upper_modes.start,pair.lower_modes.start] {
            assert!(row[start+1..start+3].iter().any(|b|b.abs()>0.),
                "sampled rim closure must include actual rigid tilt");
        }
    }
    let (mut inputs,spatial)=shaft_playing::player_inputs(&e,[
        Some(mechanics::drive::Program::parse("0,0\n0.0003,0.03\n0.0015,0").unwrap()),
        Some(mechanics::drive::Program::parse("0,0\n0.0004,-0.015\n0.0015,0").unwrap())]).unwrap();
    assert_eq!(spatial.len(),2);assert_ne!(spatial[0].weights,spatial[1].weights);
    inputs.push(mechanics::drive::Input{program:mechanics::drive::Program::parse(
        "0,0\n0.0002,0.01\n0.0015,0").unwrap(),coordinate:pair.pedal.coordinate,tip_weight:pair.pedal.weight});
    let initial=inner(&e.system).stored_energy_j();
    e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
        ImpactSubstepConfig{max_depth:8,max_attempts:511}).unwrap()
        .with_player_drives(inputs,spatial,dt,steps,n).unwrap();
    let gate=CancelGate::new_clock_free();let mut peaks=[0.0_f64;2];let mut net=0.;let mut loss=0.;
    for tick in 0..steps {
        if tick==320 {
            let before=e.system.state().to_vec();let histories=(0..12)
                .map(|i|inner(&e.system).felt_history(i)).collect::<Vec<_>>();
            let staged=match &mut e.system {Mechanics::Driven{drive,..}=>drive.forces(&e.force).unwrap().to_vec(),_=>unreachable!()};
            let cancel=CancelGate::new_clock_free();cancel.request();assert!(e.system.step(&e.force,&cancel).is_err());
            let mut bad=e.force.clone();bad[0]=f64::NAN;assert!(e.system.step(&bad,&gate).is_err());
            assert_eq!(e.system.state(),before);
            assert_eq!((0..12).map(|i|inner(&e.system).felt_history(i)).collect::<Vec<_>>(),histories);
            if let Mechanics::Driven{drive,..}=&mut e.system {assert_eq!(drive.forces(&e.force).unwrap(),staged);}
        }
        let f=e.system.step(&e.force,&gate).unwrap();net+=f.supplied_work_j-f.dissipated_energy_j;
        loss+=f.dissipated_energy_j;assert_eq!(f.time_s,(tick+1) as f64*dt);
        assert!((f.stored_energy_j-initial-net).abs()<1e-6);
        for (i,peak) in peaks.iter_mut().enumerate() {
            let r=rocking::observation(&e.system,i+1).unwrap();
            *peak=peak.max(r.iter().map(|a|a*a).sum::<f64>().sqrt());assert!(*peak<=cap);
        }
    }
    assert!(peaks.iter().all(|p|*p>1e-8),"both off-center strikes must excite their shell's angular inertia: {peaks:?}");
    assert!(loss>0.);
    let gas=play::gas_observation(&e.system).unwrap().unwrap();
    assert!(gas.minimum_pressure_pa>0. && gas.maximum_pressure_pa>gas.minimum_pressure_pa);
    // CSV rotations are axial vectors in the WORLD frame on both bodies.
    let mut columns=Vec::new();rocking::row(&e.system,1,false,&mut columns).unwrap();
    rocking::row(&e.system,2,true,&mut columns).unwrap();
    let values=String::from_utf8(columns).unwrap().split(',').skip(1)
        .map(|s|s.parse::<f64>().unwrap()).collect::<Vec<_>>();
    let upper=rocking::observation(&e.system,1).unwrap();let lower=rocking::observation(&e.system,2).unwrap();
    assert_eq!(values,[upper[0],upper[1],upper[2],lower[0],-lower[1],-lower[2]]);
}
