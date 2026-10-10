use super::*;
use fs_couple::render::plate::impact::ImpactSubstepConfig;
use fs_couple::render::plate::impact::striker::RadiusStation;
use fs_plate::shell::stiffened::beam::RoundBeamSpec;

const CARD:&str="frankensim-felt-mallet-v1\ngeometry,0.02,0.003,0.006,0.00002\nfelt,1000000,0.2,2.2,3,0.15,0.7\nconditioning,0\ncreep,3000,6\n";
const GAS:&str="frankensim-squeeze-film-v1\nannulus,0.04,0.085,2,4\nviscosity_pa_s,0.000018\nlimits,0.000001,0.04,200000\nboundary,sealed,sealed\nisothermal,100000,293.15,287.05\n";
fn spec()->Spec {Spec::parse(include_str!("estimated-hihat.fshh")).unwrap()}
fn shell()->specimen::Specimen {let mut s=specimen::Specimen::reference();s.azimuths=8;s}
fn first()->Stroke {Stroke{position_m:Some([0.06,0.01]),speed_m_s:0.5}}
fn second()->Stroke {Stroke{position_m:Some([-0.06,0.01]),speed_m_s:0.3}}
fn tip()->mallets::Spec {mallets::Spec::parse(CARD).unwrap()}
fn loaded_tip()->mallets::Spec {
    mallets::Spec::parse(&CARD.replace("-v1","-v2").replace("conditioning,0",
        "attachment,0.0000012,0\nconditioning,0")).unwrap()
}
fn shaft()->FlexibleStriker {
    FlexibleStriker::new(&[(0.0,0.005),(0.4,0.005)].map(|(position_m,radius_m)|
        RadiusStation{position_m,radius_m}),RoundBeamSpec{young_pa:12e9,density_kg_m3:800.,
        pivot_m:0.1,contact_m:0.39,hand_m:0.16,subdivisions:8,maximum_hz:3000.,maximum_modes:17},0.001).unwrap()
}
fn inner(m:&Mechanics)->&ImpactSystem {
    match m {Mechanics::Reference(s)=>s,Mechanics::Nonlinear(s)=>s,
        Mechanics::Substepped(s)=>s,Mechanics::Driven{inner:s,..}=>inner(s),
        Mechanics::Prepared(_)=>panic!("paired felt contact cannot drop material history")}
}

#[test]
fn second_felt_mallet_retains_both_shells_twelve_mounts_and_compressible_gap_air() {
    let tips=mallets::Selection{first:None,second:Some(tip())};
    let shafts=shaft_playing::Selection{first:Some(shaft()),second:None};
    let film=squeeze::Config::parse(GAS).unwrap();let s=shell();
    let pair=build_with_mallets(&spec(),&s,&s,first(),Some(second()),512,2e-6,true,
        Some(&film),&shafts,&tips).unwrap();
    let e=pair.experiment;let n=e.force.len();let hand=e.second_stick.unwrap();
    let p=e.flexible_sticks[0].as_ref().unwrap();
    assert!(e.flexible_sticks[1].is_none());
    assert_eq!(hand.coordinate,pair.pedal.coordinate+1);
    assert_eq!(p.elastic_start(),hand.coordinate+1);
    assert_eq!(e.system.state().len(),2*n+16+8);
    assert!((hand.weight-1./0.02_f64.sqrt()).abs()<1e-14);
    assert!((e.system.state()[2*hand.coordinate]*hand.weight+0.00002).abs()<1e-17);
    assert!((e.system.state()[2*hand.coordinate+1]*hand.weight-second().speed_m_s).abs()<1e-14);
    for i in 0..16 {assert!(inner(&e.system).felt_history(i).is_some());}
    assert!(inner(&e.system).felt_history(16).is_none());
    for i in 12..16 {assert_eq!(inner(&e.system).felt_observation(i).unwrap().1,0.);}
    assert_eq!(e.acoustics.as_ref().unwrap().state_modes(),
        pair.upper_modes.clone().chain(pair.lower_modes.clone()).collect::<Vec<_>>());
    for row in pair.collision.collocation().chunks(n) {
        assert_eq!(row[0],0.);assert_eq!(row[hand.coordinate],0.);
        assert_eq!(row[pair.pedal.coordinate],0.);
        assert!(row[p.elastic_start()..p.elastic_start()+p.elastic_modes()].iter().all(|v|*v==0.));
    }
    let gas=play::gas_observation(&e.system).unwrap().unwrap();
    assert_eq!(gas.minimum_pressure_pa,100000.);assert_eq!(gas.free_energy_j,0.);
    // Admission never trades away stand history or reinterprets effective mass.
    let two=mallets::Selection{first:Some(tip()),second:Some(tip())};
    let none=shaft_playing::Selection::default();
    let error=build_with_mallets(&spec(),&s,&s,first(),Some(second()),1,2e-6,false,None,&none,&two)
        .err().unwrap().to_string();
    assert!(error.contains("sixteen-pad limit"));
    let physical=mallets::Selection{first:Some(loaded_tip()),second:None};
    assert!(build_with_mallets(&spec(),&s,&s,first(),None,1,2e-6,false,None,&none,&physical).is_err());
    let effective=mallets::Selection{first:Some(tip()),second:None};
    assert!(build_with_mallets(&spec(),&s,&s,first(),None,1,2e-6,false,None,&shafts,&effective).is_err());
}

#[test]
fn loaded_felt_hi_hat_retains_contact_moments_and_one_pedal_and_hand_clock() {
    let make=|| {
        let tips=mallets::Selection{first:Some(loaded_tip()),second:None};
        let shafts=shaft_playing::Selection{first:Some(shaft()),second:Some(shaft())};
        let s=shell();let other=Stroke{speed_m_s:0.,..second()};
        let pair=build_with_mallets(&spec(),&s,&s,first(),Some(other),512,2e-6,false,None,
            &shafts,&tips).unwrap();
        let mut e=pair.experiment;let n=e.force.len();
        let p=e.flexible_sticks[0].as_ref().unwrap();
        let o=p.observe(e.system.state()).unwrap();
        assert!((o.tip_displacement_m+0.00002).abs()<1e-17);
        assert!((o.tip_velocity_m_s-first().speed_m_s).abs()<1e-14);
        assert_eq!(o.flexural_energy_j,0.);
        assert_ne!(p.tip_row(n).unwrap(),p.hand_row(n).unwrap());
        let (mut inputs,spatial)=shaft_playing::player_inputs(&e,[
            Some(mechanics::drive::Program::parse("0,0\n0.0002,0.05\n0.001,0").unwrap()),
            Some(mechanics::drive::Program::parse("0,0\n0.0003,-0.02\n0.001,0").unwrap())]).unwrap();
        inputs.push(mechanics::drive::Input{program:mechanics::drive::Program::parse(
            "0,0\n0.0003,0.2\n0.001,0").unwrap(),coordinate:pair.pedal.coordinate,tip_weight:pair.pedal.weight});
        e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
            ImpactSubstepConfig{max_depth:8,max_attempts:511}).unwrap()
            .with_player_drives(inputs,spatial,2e-6,512,n).unwrap();e
    };
    let mut e=make();let mut clean=make();let initial=inner(&e.system).stored_energy_j();
    let gate=CancelGate::new_clock_free();let (mut net,mut loss)=(0.,0.);
    let (mut contact,mut flexure)=(0.0_f64,0.0_f64);
    for tick in 0..512 {
        if tick==160 {
            let before=e.system.state().to_vec();let histories=(0..16)
                .map(|i|inner(&e.system).felt_history(i)).collect::<Vec<_>>();
            let cancel=CancelGate::new_clock_free();cancel.request();
            assert!(e.system.step(&e.force,&cancel).is_err());
            let mut bad=e.force.clone();bad[0]=f64::NAN;
            assert!(e.system.step(&bad,&gate).is_err());
            assert_eq!(e.system.state(),before);
            assert_eq!((0..16).map(|i|inner(&e.system).felt_history(i)).collect::<Vec<_>>(),histories);
        }
        let f=e.system.step(&e.force,&gate).unwrap();clean.system.step(&clean.force,&gate).unwrap();
        assert_eq!(e.system.state(),clean.system.state());assert_eq!(f.time_s,(tick+1) as f64*2e-6);
        net+=f.supplied_work_j-f.dissipated_energy_j;loss+=f.dissipated_energy_j;
        assert!((f.stored_energy_j-initial-net).abs()<1e-6);
        for i in 12..16 {contact=contact.max(inner(&e.system).felt_observation(i).unwrap().1);}
        flexure=flexure.max(e.flexible_sticks[0].as_ref().unwrap().observe(e.system.state()).unwrap().flexural_energy_j);
    }
    assert!(contact>0. && flexure>1e-12 && loss>0.);
    assert!((12..16).any(|i|inner(&e.system).felt_history(i).unwrap().eps_max>0.));
}

#[test]
fn a_felt_hi_hat_has_no_parallel_hard_tip_before_its_actual_skin_contact() {
    let tips=mallets::Selection{first:Some(mallets::Spec::parse(&CARD.replace("0.003","0.012")).unwrap()),second:None};
    let shafts=shaft_playing::Selection::default();let s=shell();
    let pair=build_with_mallets(&spec(),&s,&s,first(),None,1024,2e-6,false,None,&shafts,&tips).unwrap();
    let mut e=pair.experiment;
    let first_gap=(12..16).map(|i|-0.006*inner(&e.system).felt_observation(i).unwrap().0)
        .fold(f64::INFINITY,f64::min);
    let extra=first_gap-0.00002;
    assert!(extra>16.*2e-6*first().speed_m_s,"fixture needs resolved curved clearance");
    let ticks=((0.00002+0.5*extra)/(2e-6*first().speed_m_s)).floor() as usize;
    assert!(ticks>0 && ticks<1024);
    e.system=e.system.into_analytic_nonlinear().unwrap();
    let gate=CancelGate::new_clock_free();
    for _ in 0..ticks {
        e.system.step(&e.force,&gate).unwrap();
        for i in 12..16 {assert_eq!(inner(&e.system).felt_observation(i).unwrap().1,0.);}
    }
    // A wide, flat face over the curved skin has positive extra site gaps.
    // Its inertial origin has crossed the old point-contact plane, but it is
    // still freely approaching those sites. An accidental Hertz tip brakes it.
    assert!(e.system.state()[0]>0.);
    assert!((e.system.state()[1]*e.stick_weight-first().speed_m_s).abs()<1e-10);
}
