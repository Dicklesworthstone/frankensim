use super::*;
use crate::{Stroke, Experiment, config, splash_with_compliant_mute, splash_with_sticks,
    drum_with_compliant_mute, drum_with_cavity_loss};
use fs_couple::render::plate::impact::{ImpactBody, ImpactSystem, ImpactSubstepConfig};
use fs_exec::CancelGate;

// Thin authored test pads resolve contact during a short physical test window.
// These are not hand-tissue measurements or the distributed example's 4mm pad.
fn text(surface:&str, opposed:bool, end:f64) -> String {
    let s=if surface=="resonant" {"below"} else {"above"};
    let mut out=format!("frankensim-compliant-mute-v1\nsurface,{surface}\n\
        site,0.06,0.01,0.0001\nsite,0.065,0.013,0.0001\n\
        jaw,{s},0.025,0.1,0,0.0005,30000,0.2,2.2,3,0.15,0.7,0\n\
        creep,{s},4000,12\nforce,{s},0,0\nforce,{s},{},1\nforce,{s},{end},0\n",end/3.0);
    if opposed {
        out.push_str(&format!("jaw,below,0.035,0.2,0,0.0005,30000,0.2,2.2,3,0.15,0.7,0\n\
            creep,below,5000,15\nforce,below,0,0\nforce,below,{},0.3\nforce,below,{end},0\n",end/2.0));
    }
    out
}
fn bounds() -> ImpactSubstepConfig {ImpactSubstepConfig {max_depth:8,max_attempts:511}}
fn rest(p:[f64;2]) -> Stroke {Stroke {speed_m_s:0.0,position_m:Some(p)}}
fn prepare(mut e:Experiment) -> Experiment {
    e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(bounds()).unwrap();
    e
}

#[test]
fn complete_cards_force_programs_and_exterior_sides_are_admitted_before_geometry() {
    let source=text("shell",true,0.001);let spec=Spec::parse(&source).unwrap();
    assert_eq!(spec.jaw_count(),2);assert_eq!(spec.sites.len(),2);
    for command in ["splash","splash-wav","splash-mic"] {spec.admit_command(command).unwrap();}
    for command in ["drum","snare","drum-modal","unknown"] {assert!(spec.admit_command(command).is_err());}
    Spec::parse(&text("batter",false,0.001)).unwrap().admit_command("drum-stretch-mic").unwrap();
    Spec::parse(&text("resonant",false,0.001)).unwrap().admit_command("drum-wav").unwrap();
    for bad in [source.replace("surface,shell","surface,shell\nsurface,shell"),
        source.replace("0.025,0.1","0,0.1"),source.replace("site,0.06,0.01,0.0001","site,NaN,0.01,0.0001"),
        source.replace("force,above,0,0","force,above,0,1"),
        source.replace("jaw,below,0.035","jaw,above,0.035"),
        text("batter",false,0.001).replace("above","below"),
        text("batter",true,0.001),source+"unknown,1\n",
        text("batter",false,0.001)+"creep,below,100,2\n",
        text("batter",false,0.001)+"force,below,0,0\nforce,below,0.001,0\n"] {
        assert!(Spec::parse(&bad).is_err(),"{bad}");
    }
    assert!(Spec::parse(&"x".repeat(65_537)).is_err());
    for (source,command) in [(include_str!("../estimated-cymbal-mute.fsm"),"splash-mic"),
        (include_str!("../estimated-drum-mute.fsm"),"drum-stretch")] {
        Spec::parse(source).unwrap().admit_command(command).unwrap();
    }
    let file=concat!(env!("CARGO_MANIFEST_DIR"),"/examples/percussion/estimated-cymbal-mute.fsm");
    let mut args=["splash","--compliant-mute",file,"--analytic-newton"].map(String::from).to_vec();
    assert_eq!(option(&mut args).unwrap().unwrap().jaw_count(),2);
    assert_eq!(args,["splash","--analytic-newton"]);

    // Compiled force inputs keep the declared jaw coordinates and physical mass.
    let a=spec.compile(vec![
        fs_couple::render::plate::impact::damping::ViscousDamper {weights:vec![0.0,1.0],damping_n_s_m:0.0};2],3).unwrap();
    let observation=spec.observation(a.ports,0);let inputs=spec.into_inputs(&observation).unwrap();
    assert_eq!(inputs[0].coordinate,3);assert_eq!(inputs[1].coordinate,4);
    assert!(drive::StickDrive::new_inputs(inputs,1e-6,10,5).is_err(),"duration must cover the full force performance");
    let mut args=vec!["splash".into(),"--compliant-mute".into()];assert!(option(&mut args).is_err());
}

#[test]
fn four_physical_force_ports_share_one_accepted_clock_and_exact_retry() {
    let masses:[f64;4]=[0.02,0.03,0.025,0.04];let dt=1e-5;let ticks=40;
    let weights=masses.map(|m|1.0/m.sqrt());
    let inputs=|| (0..4).map(|i|drive::Input {coordinate:i,tip_weight:weights[i],
        program:drive::Program::parse(&format!("0,0\n0.0002,{}\n0.0004,0",[1.0,-0.5,0.7,-0.2][i])).unwrap()}).collect();
    let make=|| Mechanics::Reference(ImpactSystem::new(masses.iter().map(|m|
        ImpactBody::free_mass(*m,0.0,0.0).unwrap().0).collect(),vec![],vec![],vec![],config(ticks,dt)).unwrap())
        .into_analytic_nonlinear().unwrap();
    let mut manual=make();let mut driven=make().with_stick_drives(inputs(),dt,ticks,4).unwrap();
    let mut staging=drive::StickDrive::new_inputs(inputs(),dt,ticks,4).unwrap();
    let gate=CancelGate::new_clock_free();let mut work=0.0;
    for tick in 0..ticks {
        if tick==15 {
            let before=driven.state().to_vec();
            assert!(driven.step(&[0.0,0.0,1e7,0.0],&gate).is_err());assert_eq!(driven.state(),before);
        }
        let f=driven.step(&[0.0;4],&gate).unwrap();
        let expected=manual.step(staging.forces(&[0.0;4]).unwrap(),&gate).unwrap();staging.accept();
        assert_eq!(driven.state(),manual.state());assert_eq!(f.supplied_work_j,expected.supplied_work_j);
        work+=f.supplied_work_j;assert!((f.stored_energy_j-work).abs()<1e-12);
    }
    for i in 0..4 {
        let impulse=0.0002*[1.0,-0.5,0.7,-0.2][i];
        assert!((driven.state()[2*i+1]*weights[i]-impulse/masses[i]).abs()<1e-12);
    }
}

#[test]
fn opposed_jaws_and_two_sticks_react_on_one_shell_without_shifting_stand_or_radiation() {
    let dt=2e-6;let ticks=256;let source=text("shell",true,ticks as f64*dt);
    let spec=Spec::parse(&source).unwrap();
    let mut e=splash_with_compliant_mute(ticks,dt,true,rest([0.06,0.01]),None,&[],
        Some(rest([-0.05,0.02])),Some(&spec)).unwrap();
    let observation=e.mute.as_ref().unwrap();let start=observation.ports[0].coordinate;
    assert_eq!(start,e.second_stick.unwrap().coordinate+1);
    assert_eq!(observation.first_pad,6);assert_eq!(e.force.len(),start+2);
    assert_eq!(e.system.state().len(),2*e.force.len()+6+4); // six stand and four independent pad memories
    assert!(e.observer_a[start..].iter().all(|v|*v==0.0));
    assert!(e.acoustics.as_ref().unwrap().state_modes.iter().all(|i|*i<start-1));
    let Mechanics::Reference(initial)=&e.system else {unreachable!()};
    let energy=initial.stored_energy_j();
    for i in 0..6 {assert!(initial.felt_history(i).is_some());}
    for i in 6..10 {assert_eq!(initial.felt_history(i).unwrap().eps_max,0.0);}
    let mut inputs=spec.into_inputs(observation).unwrap();
    inputs.push(drive::Input {coordinate:0,tip_weight:e.stick_weight,
        program:drive::Program::parse("0,0\n0.0001,0.02\n0.0005,0").unwrap()});
    let second=e.second_stick.unwrap();
    inputs.push(drive::Input {coordinate:second.coordinate,tip_weight:second.weight,
        program:drive::Program::parse("0,0\n0.0002,-0.01\n0.0005,0").unwrap()});
    e=prepare(e);e.system=e.system.with_stick_drives(inputs,dt,ticks,e.force.len()).unwrap();
    let gate=CancelGate::new_clock_free();let mut net=0.0;let mut contact=0.0_f64;
    for tick in 1..=ticks {
        let f=e.system.step(&e.force,&gate).unwrap();net+=f.supplied_work_j-f.dissipated_energy_j;
        assert_eq!(f.time_s,tick as f64*dt);assert!((f.stored_energy_j-energy-net).abs()<1e-6);
        contact=contact.max(felt(&e.system,6).unwrap().1);
    }
    assert!(contact>0.0,"surface contact, not a disconnected force-driven jaw");
    assert!(e.system.state()[2..2*(start-1)].iter().any(|v|v.abs()>1e-14),"jaw reaction must reach the original shell");
    let mut header=Vec::new();let mut row=Vec::new();let o=e.mute.as_ref().unwrap();
    o.header(&mut header).unwrap();o.row(&e.system,&mut row).unwrap();
    assert_eq!(String::from_utf8(header).unwrap().split(',').count(),7);
    assert_eq!(String::from_utf8(row).unwrap().split(',').count(),7);
}

#[test]
fn exterior_head_pad_keeps_cavity_pressure_and_private_material_coordinates_reciprocal() {
    let dt=2e-6;let ticks=128;let source=text("resonant",false,0.00024);
    let spec=Spec::parse(&source).unwrap();
    let mut e=drum_with_compliant_mute(ticks,dt,true,false,None,true,
        Stroke {speed_m_s:4.0,position_m:Some([0.06,0.01])},true,None,None,
        Some(Stroke {speed_m_s:2.5,position_m:Some([-0.05,0.02])}),&[],20.0,Some(&spec)).unwrap();
    let o=e.mute.as_ref().unwrap();let jaw=o.ports[0];let air=e.air.as_ref().unwrap();
    assert_eq!(jaw.coordinate,e.second_stick.unwrap().coordinate+1);
    assert_eq!(air.coupling.structural_modes(),jaw.coordinate+1);
    assert_eq!(e.system.state().len(),2*air.coupling.total_modes()+2);
    assert!(e.acoustics.as_ref().unwrap().state_modes.iter().all(|i|*i<jaw.coordinate));
    assert_eq!(e.pressure.as_ref().unwrap().areas[jaw.coordinate],0.0);
    let mut displaced=e.system.state().to_vec();displaced[2*jaw.coordinate]+=0.001/jaw.inverse_sqrt_mass;
    assert_eq!(air.uniform_pressure(&displaced).unwrap(),0.0,"jaw must not directly compress cavity gas");
    assert_eq!(air.points(&displaced).unwrap(),(0.0,0.0));
    let inputs=spec.into_inputs(o).unwrap();
    let Mechanics::Reference(initial)=&e.system else {unreachable!()};let energy=initial.stored_energy_j();
    e=prepare(e);e.system=e.system.with_stick_drives(inputs,dt,ticks,e.force.len()).unwrap();
    let gate=CancelGate::new_clock_free();let mut net=0.0;let mut pressure=0.0_f64;
    for tick in 0..ticks {
        if tick==32 {
            let before=e.system.state().to_vec();let mut bad=e.force.clone();bad[jaw.coordinate]=1e7;
            assert!(e.system.step(&bad,&gate).is_err());assert_eq!(e.system.state(),before);
        }
        let f=e.system.step(&e.force,&gate).unwrap();net+=f.supplied_work_j-f.dissipated_energy_j;
        assert!((f.stored_energy_j-energy-net).abs()<1e-6);
        pressure=pressure.max(e.air.as_ref().unwrap().uniform_pressure(e.system.state()).unwrap().abs());
    }
    assert!(pressure>0.0);assert!(felt(&e.system,0).is_some());
    assert!(e.system.membrane_observation(1).unwrap().stretching_energy_j>0.0);
}

#[test]
fn no_mute_preserves_legacy_trajectories_and_missing_surface_sites_do_not_snap() {
    let stroke=rest([0.06,0.01]);let gate=CancelGate::new_clock_free();
    let mut old=splash_with_sticks(2,2e-6,false,stroke,None,&[],None).unwrap();
    let mut new=splash_with_compliant_mute(2,2e-6,false,stroke,None,&[],None,None).unwrap();
    assert!(new.mute.is_none());
    for _ in 0..2 {
        let a=old.system.step(&old.force,&gate).unwrap();let b=new.system.step(&new.force,&gate).unwrap();
        assert_eq!(old.system.state(),new.system.state());assert_eq!(a.stored_energy_j,b.stored_energy_j);
    }
    let mut old=drum_with_cavity_loss(2,2e-6,false,false,None,false,stroke,false,None,None,None,&[],0.0).unwrap();
    let mut new=drum_with_compliant_mute(2,2e-6,false,false,None,false,stroke,false,None,None,None,&[],0.0,None).unwrap();
    for _ in 0..2 {old.system.step(&old.force,&gate).unwrap();new.system.step(&new.force,&gate).unwrap();
        assert_eq!(old.system.state(),new.system.state());}
    let source=text("shell",false,0.001).replace("site,0.06,0.01,0.0001","site,0,0,0.0001");
    let bad=Spec::parse(&source).unwrap();
    assert!(splash_with_compliant_mute(1,2e-6,false,stroke,None,&[],None,Some(&bad)).is_err());
}
