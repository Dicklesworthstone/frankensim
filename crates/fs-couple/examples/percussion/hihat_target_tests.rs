use super::*;
use targets::{Selection,Target};
use fs_couple::render::plate::impact::ImpactSubstepConfig;
use fs_couple::render::plate::impact::striker::RadiusStation;
use fs_plate::shell::stiffened::beam::RoundBeamSpec;

const CARD:&str="frankensim-felt-mallet-v2\ngeometry,0.02,0.003,0.006,0.00002\nattachment,0.0000012,0.61\nfelt,1000000,0.2,2.2,3,0.15,0.7\nconditioning,0\ncreep,3000,6\n";
const GAS:&str="frankensim-squeeze-film-v1\nannulus,0.04,0.085,2,4\nviscosity_pa_s,0.000018\nlimits,0.000001,0.04,200000\nboundary,sealed,sealed\nisothermal,100000,293.15,287.05\n";
fn spec()->Spec {Spec::parse(include_str!("estimated-hihat.fshh")).unwrap()}
fn shell()->specimen::Specimen {let mut s=specimen::Specimen::reference();s.azimuths=8;s}
fn stroke()->Stroke {Stroke{position_m:Some([0.06,0.01]),speed_m_s:0.5}}
fn shaft()->FlexibleStriker {
    FlexibleStriker::new(&[(0.0,0.005),(0.4,0.005)].map(|(position_m,radius_m)|
        RadiusStation{position_m,radius_m}),RoundBeamSpec{young_pa:12e9,density_kg_m3:800.,
        pivot_m:0.1,contact_m:0.39,hand_m:0.16,subdivisions:8,maximum_hz:3000.,maximum_modes:17},0.001).unwrap()
}
fn inner(m:&Mechanics)->&ImpactSystem {
    match m {Mechanics::Reference(s)=>s,Mechanics::Nonlinear(s)=>s,Mechanics::Substepped(s)=>s,
        Mechanics::Driven{inner:s,..}=>inner(s),Mechanics::Prepared(_)=>panic!("paired felt history is required")}
}
fn bounds()->ImpactSubstepConfig {ImpactSubstepConfig{max_depth:8,max_attempts:511}}

#[test]
fn targets_are_independent_exterior_choices_and_keep_legacy_arguments() {
    let args=|s:&str|s.split_whitespace().map(str::to_string).collect::<Vec<_>>();
    let original=args("hihat input.fshh 512 --strike-speed-m-s 0.5");
    let mut plain=original.clone();
    assert_eq!(Selection::options(&mut plain,None).unwrap(),Selection::default());
    assert_eq!(plain,original);
    for (a,b) in [(Target::Upper,Target::Lower),(Target::Lower,Target::Upper),(Target::Lower,Target::Lower)] {
        let mut input=args(&format!("hihat input.fshh 512 --strike-target {} --strike-speed-m-s 0.5 --second-stick-target {}",a.label(),b.label()));
        assert_eq!(Selection::options(&mut input,Some(stroke())).unwrap(),Selection{first:a,second:b});
        assert_eq!(input,original);
    }
    for bad in ["--strike-target","--strike-target inner","--strike-target lower --strike-target upper",
        "--second-stick-target lower --second-stick-target lower"] {
        assert!(Selection::options(&mut args(bad),Some(stroke())).is_err());
    }
    for target in ["upper","lower"] {
        assert!(Selection::options(&mut args(&format!("--second-stick-target {target}")),None).is_err());
    }
}

#[test]
fn lower_loaded_mallet_preserves_world_station_gaps_and_force_moment_power() {
    let s=Shell::new(&shell(),2e-6).unwrap();let tip=mallets::Spec::parse(CARD).unwrap();
    let tips=mallets::Selection{first:Some(tip.clone()),second:None};
    let end=1+s.reduction.mode_count();
    let mut shaft=shaft_playing::Selection{first:Some(shaft()),second:None}
        .build_with_mallets(end,end,stroke(),None,2e-6,&tips).unwrap();
    let raw=Target::Lower.mallet(&tip,&s,stroke(),0,1,shaft.total).unwrap();
    let gaps=raw.pads.iter().map(|p|p.precompression_m+p.weights[0]
        *raw.body.initial[0].displacement_m_sqrt_kg).collect::<Vec<_>>();
    let bound=tip.bind_shaft(raw,0,&mut shaft).unwrap();let p=shaft.ports[0].as_ref().unwrap();
    let n=shaft.total;let mut state=vec![0.;2*n];
    state[0]=bound.body.initial[0].displacement_m_sqrt_kg;
    for i in 0..n {state[2*i+1]=0.003*(i+1) as f64;}
    let o=p.observe(&state).unwrap();let a=0.003/std::f64::consts::SQRT_2;
    let center=stroke().position_m.unwrap();let offsets=[[a,0.],[0.,a],[-a,0.],[0.,-a]];
    let loads=[1.0,0.4,0.2,0.7];let mut generalized=0.;let mut physical=0.;let mut offset_work=0.0_f64;
    for (i,offset) in offsets.into_iter().enumerate() {
        // Independently place the physical world station on the rotated body.
        let world=[center[0]+offset[0],center[1]+offset[1]];
        let surface=s.port([world[0],-world[1]],ShellFace::Positive).unwrap();
        assert!((surface.position_m[0]-world[0]).abs()<1e-12);
        assert!((surface.position_m[1]+world[1]).abs()<1e-12);
        let pad=&bound.pads[i];assert_eq!(&pad.weights[1..end],surface.weights.as_slice());
        assert!((pad.precompression_m+pad.weights[0]*state[0]-gaps[i]).abs()<1e-15);
        let outward=surface.weights.iter().enumerate().map(|(k,b)|b*state[2*(k+1)+1]).sum::<f64>();
        let lever=offset[0]*0.61_f64.cos()+offset[1]*0.61_f64.sin();
        physical+=loads[i]*(outward+o.tip_velocity_m_s+lever*o.tip_angular_velocity_rad_s);
        generalized+=loads[i]*pad.weights.iter().enumerate().map(|(k,b)|b*state[2*k+1]).sum::<f64>();
        let (f,b)=playing::shell_location(&s.outer_nodes,&s.mesh.mesh.tris,[world[0],-world[1]]).unwrap();
        for (skin,middle) in surface.weights.iter().zip(s.reduction.point_port(f,b,[0.,0.,1.]).unwrap()) {
            offset_work=offset_work.max((skin-middle).abs());
        }
    }
    assert!((generalized-physical).abs()<1e-12*physical.abs().max(1.));
    assert!(offset_work>1e-12,"outer skin must retain thickness-offset force moments");
    assert!(shaft.bodies[0].is_none());
}

#[test]
fn second_lower_wood_strike_moves_the_lower_shell_inward_with_independent_hands() {
    let s=shell();let mut spec=spec();spec.separation=0.01;
    let first=Stroke{speed_m_s:0.,..stroke()};
    let second=Stroke{position_m:Some([-0.06,0.01]),speed_m_s:0.8};
    let shafts=shaft_playing::Selection{first:Some(shaft()),second:Some(shaft())};
    let pair=build_with_targets(&spec,&s,&s,first,Some(second),256,2e-6,true,None,&shafts,
        &mallets::Selection::default(),None,None,Selection{first:Target::Upper,second:Target::Lower}).unwrap();
    let mut e=pair.experiment;let n=e.force.len();
    assert_eq!(e.system.state().len(),2*n+12);
    assert_eq!(e.acoustics.as_ref().unwrap().state_modes(),
        pair.upper_modes.clone().chain(pair.lower_modes.clone()).collect::<Vec<_>>());
    let (inputs,spatial)=shaft_playing::player_inputs(&e,[
        Some(mechanics::drive::Program::parse("0,0\n0.0001,-0.05\n0.0005,0").unwrap()),
        Some(mechanics::drive::Program::parse("0,0\n0.00015,0.05\n0.0005,0").unwrap())]).unwrap();
    assert!(inputs.is_empty());assert_eq!(spatial.len(),2);
    assert_ne!(spatial[0].weights,spatial[1].weights);
    let initial=inner(&e.system).stored_energy_j();let initial_first=e.system.state()[0];
    e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(bounds()).unwrap()
        .with_player_drives(inputs,spatial,2e-6,256,n).unwrap();
    let gate=CancelGate::new_clock_free();let mut net=0.;let mut inward=0.0_f64;let mut flexure=0.0_f64;
    for tick in 0..256 {
        let f=e.system.step(&e.force,&gate).unwrap();net+=f.supplied_work_j-f.dissipated_energy_j;
        assert_eq!(f.time_s,(tick+1) as f64*2e-6);
        assert!((f.stored_energy_j-initial-net).abs()<1e-6);
        let x=e.system.state();inward=inward.max(-x[2*pair.lower_modes.start]);
        assert!(pair.upper_modes.clone().all(|i|x[2*i].abs()<1e-10),"no collision or air path may move the unstruck upper shell");
        flexure=flexure.max(e.flexible_sticks[1].as_ref().unwrap().observe(x).unwrap().flexural_energy_j);
    }
    assert!(inward>1e-10 && flexure>1e-12,"lower contact must move its shell upward and bend its actual shaft");
    assert!(e.system.state()[0]<initial_first,"negative first-hand force retracts that hand's rigid coordinate");
}

#[test]
fn first_lower_felt_strike_keeps_gap_air_mount_history_and_tick_rollback() {
    let make=|| {
        let s=shell();let tips=mallets::Selection{first:Some(mallets::Spec::parse(CARD).unwrap()),second:None};
        let shafts=shaft_playing::Selection{first:Some(shaft()),second:Some(shaft())};
        let film=squeeze::Config::parse(GAS).unwrap();
        let pair=build_with_targets(&spec(),&s,&s,stroke(),Some(Stroke{speed_m_s:0.,
            position_m:Some([-0.06,0.01])}),512,2e-6,false,Some(&film),&shafts,&tips,None,None,
            Selection{first:Target::Lower,second:Target::Upper}).unwrap();
        let mut e=pair.experiment;let n=e.force.len();assert_eq!(e.system.state().len(),2*n+16+8);
        for i in 0..16 {assert!(inner(&e.system).felt_history(i).is_some());}
        assert!(inner(&e.system).felt_history(16).is_none());
        e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(bounds()).unwrap();e
    };
    let mut e=make();let mut clean=make();let initial=inner(&e.system).stored_energy_j();
    let gate=CancelGate::new_clock_free();let mut loss=0.;let mut contact=0.0_f64;let mut flexure=0.0_f64;
    for tick in 0..512 {
        if tick==160 {
            let before=e.system.state().to_vec();let histories=(0..16).map(|i|inner(&e.system).felt_history(i)).collect::<Vec<_>>();
            let cancel=CancelGate::new_clock_free();cancel.request();assert!(e.system.step(&e.force,&cancel).is_err());
            let mut bad=e.force.clone();bad[0]=f64::NAN;assert!(e.system.step(&bad,&gate).is_err());
            assert_eq!(e.system.state(),before);
            assert_eq!((0..16).map(|i|inner(&e.system).felt_history(i)).collect::<Vec<_>>(),histories);
        }
        let f=e.system.step(&e.force,&gate).unwrap();clean.system.step(&clean.force,&gate).unwrap();
        assert_eq!(e.system.state(),clean.system.state());loss+=f.dissipated_energy_j;
        assert_eq!(f.supplied_work_j,0.);assert!((f.stored_energy_j+loss-initial).abs()<1e-6);
        contact=contact.max((12..16).map(|i|inner(&e.system).felt_observation(i).unwrap().1).sum());
        flexure=flexure.max(e.flexible_sticks[0].as_ref().unwrap().observe(e.system.state()).unwrap().flexural_energy_j);
    }
    assert!(contact>0. && flexure>1e-12 && loss>0.);
    assert!((12..16).any(|i|inner(&e.system).felt_history(i).unwrap().eps_max>0.));
    assert!(play::gas_observation(&e.system).unwrap().unwrap().minimum_pressure_pa>0.);
}
