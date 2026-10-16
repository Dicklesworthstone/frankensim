use super::*;
use super::super::super::{drum_with_radiation, drum_spec, Experiment, Stroke,
    acoustics, mechanics, nonlinear_snare};
use fs_couple::render::plate::impact::{ImpactSubstepConfig, string::StringStretching};
use fs_exec::CancelGate;

const EXAMPLE:&str=include_str!("../../estimated-snare-carrier.fsc");
const SHORT:&str="frankensim-snare-carrier-v1\nsupport,0.01,12,0.02,0.01,0.2\ninitial,-0.000022,0\nforce,0,0\nforce,0.000256,0.03\nforce,0.001024,0\n";
fn bank(stretching:bool)->SnareSet {
    SnareSet {carrier:Some(Specification::parse(SHORT).unwrap().support),
        stretching:stretching.then_some(StringStretching {axial_rigidity_n:100.0,maximum_slope:0.2}),
        ..SnareSet::reference(false)}
}
fn build(s:SnareSet,heads:bool,air:bool,second:bool,audio:bool)->Experiment {
    drum_with_radiation(512,if audio {acoustics::MECHANICAL_DT}else{2e-6},audio,
        !(heads || s.stretching.is_some() || s.carrier.is_some()),Some(s),heads,
        Stroke {speed_m_s:0.0,position_m:Some([0.06,0.01])},air,None,
        Some(drum_spec::Spec {radial_intervals:2,azimuths:8,..drum_spec::Spec::reference()}),
        second.then_some(Stroke {speed_m_s:0.0,position_m:Some([-0.05,0.02])}),
        &[],if air {25.0}else{0.0},None,false).unwrap()
}
fn energy(s:&Mechanics)->f64 {
    match s {
        Mechanics::Reference(s)=>s.stored_energy_j(),Mechanics::Nonlinear(s)=>s.stored_energy_j(),
        Mechanics::Substepped(s)=>s.stored_energy_j(),Mechanics::Driven {inner,..}=>energy(inner),
        Mechanics::Prepared(_)=>panic!("carrier must use the coupled support dynamics"),
    }
}

#[test]
fn complete_carrier_card_is_required_and_incompatible_commands_refuse_before_loading() {
    let s=Specification::parse(EXAMPLE).unwrap();assert_eq!(s.support.mass_kg,0.02);
    for (a,b) in [("support,0.01","support,0"),("support,0.01","support,NaN"),
        ("12,0.02","-12,0.02"),("12,0.02","12,-0.02"),("0.01,0.2","0,0.2"),
        ("0.01,0.2","0.01,0.4"),("initial,-0.000022,0","initial,0.1,0"),
        ("force,0,0","force,0,1"),("force,0,0","force,-1,0"),
        ("force,0.001024,0","force,0.000256,0"),("force,0.001024,0","force,0.001024,1"),
        ("force,0.000256,0.03","force,NaN,0.03"),("force,0.000256,0.03","force,0.000256,inf"),
        ("support,0.01,12,0.02,0.01,0.2", ""),("initial,-0.000022,0", ""),
        ("frankensim-snare-carrier-v1","frankensim-snare-carrier-v2")] {
        assert!(Specification::parse(&SHORT.replace(a,b)).is_err(),"accepted {b}");
    }
    for row in ["initial,0,0","support,1,0,0,0.1,0.2","unknown,1","force,0.1"] {
        assert!(Specification::parse(&format!("{SHORT}{row}\n")).is_err());
    }
    assert!(Specification::parse(&" ".repeat(MAX_BYTES as usize+1)).is_err());
    assert!(select(Some("must-not-be-opened"),None).is_err());
    let path=format!("{}/examples/percussion/estimated-snare-carrier.fsc",env!("CARGO_MANIFEST_DIR"));
    for command in ["snare","snare-off","snare-wav","snare-off-wav","snare-mic","snare-off-mic"] {
        let s=super::super::spec::select(None,command).unwrap();
        let (s,carrier)=select(Some(&path),s).unwrap();assert!(carrier.is_some());
        assert_eq!(s.unwrap().clearance_m,SnareSet::reference(command.starts_with("snare-off")).clearance_m);
        nonlinear_snare::admit_prepared_command(true,true,command).unwrap();
    }
    let mut args=vec!["snare".into(),"--snare-carrier".into(),path.clone(),"--analytic-newton".into()];
    assert_eq!(option(&mut args).unwrap(),Some(path));assert_eq!(args,["snare","--analytic-newton"]);
    assert!(option(&mut vec!["--snare-carrier".into()]).is_err());
    assert!(option(&mut vec!["--snare-carrier".into(),"--snare-spec".into()]).is_err());
    assert!(option(&mut vec!["--snare-carrier".into(),"a".into(),"--snare-carrier".into(),"b".into()]).is_err());
}

#[test]
fn full_carried_bank_keeps_every_wire_coordinate_and_the_same_radiating_heads() {
    for stretching in [false,true] {for heads in [false,true] {
        let s=bank(stretching);assert_eq!(s.mode_count().unwrap(),161);
        let e=build(s,heads,true,true,true);
        let first=e.second_stick.unwrap().coordinate+1;
        let (port,weight)=force_port(&e.system,4).unwrap();assert_eq!(port,first+160);
        assert_eq!(e.air.as_ref().unwrap().coupling.structural_modes(),first+161);
        assert!(weight.is_finite() && weight>0.0 && e.force.len()<=256);
        assert!(force_port(&e.system,5).is_none());
        let o=observe(&e.system,4).unwrap();assert!((o.position_m+0.000022).abs()<1e-18);
        assert!(o.maximum_wire_slope<1e-14);assert_eq!(o.velocity_m_s,0.0);
        assert_eq!(e.system.membrane_observation(1).is_some(),heads);
        assert!(e.observer_a[first..].iter().chain(&e.observer_b[first..]).all(|v|*v==0.0));
        assert!(e.pressure.as_ref().unwrap().areas[first..].iter().all(|v|*v==0.0));
        assert!(e.acoustics.as_ref().unwrap().state_modes().iter().all(|i|*i<first));
        assert!(drum_with_radiation(1,2e-6,false,true,Some(s),heads,Stroke::default(),false,
            None,None,None,&[],0.0,None,false).is_err(),"linear-only image cannot discard support inertia");
    }}
}

#[test]
fn force_driven_carrier_and_two_sticks_share_cavity_state_work_and_atomic_retry() {
    let s=SnareSet {strands:2,modes_per_strand:2,contact_cells:4,..bank(true)};
    let make=|| {
        let mut e=build(s,true,true,true,false);
        e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
            ImpactSubstepConfig {max_depth:4,max_attempts:31}).unwrap();e
    };
    let inputs=|e:&Experiment| {
        let other=e.second_stick.unwrap();
        vec![Input {coordinate:0,tip_weight:e.stick_weight,program:Program::parse("0,0\n0.000256,0.02\n0.001024,0").unwrap()},
            Input {coordinate:other.coordinate,tip_weight:other.weight,program:Program::parse("0,0\n0.000256,-0.01\n0.001024,0").unwrap()},
            Specification::parse(SHORT).unwrap().into_input(&e.system,4).unwrap()]
    };
    let mut driven=make();let mut manual=make();let initial=energy(&driven.system);
    let mut staging=mechanics::drive::StickDrive::new_inputs(inputs(&manual),2e-6,512,manual.force.len()).unwrap();
    let drive_inputs=inputs(&driven);
    driven.system=driven.system.with_stick_drives(drive_inputs,2e-6,512,driven.force.len()).unwrap();
    let gate=CancelGate::new_clock_free();let cancel=CancelGate::new_clock_free();cancel.request();
    let (mut net,mut work,mut head,mut pressure)=(0.0,0.0,0.0_f64,0.0_f64);
    for tick in 0..128 {
        if tick==64 {let before=driven.system.state().to_vec();assert!(driven.system.step(&driven.force,&cancel).is_err());assert_eq!(driven.system.state(),before);}
        let f=driven.system.step(&driven.force,&gate).unwrap();
        let g=manual.system.step(staging.forces(&manual.force).unwrap(),&gate).unwrap();staging.accept();
        assert_eq!(driven.system.state(),manual.system.state());assert_eq!(f.time_s,g.time_s);
        assert_eq!(f.time_s,(tick+1) as f64*2e-6);assert_eq!(f.supplied_work_j,g.supplied_work_j);
        net+=f.supplied_work_j-f.dissipated_energy_j;work+=f.supplied_work_j.abs();
        assert!((f.stored_energy_j-initial-net).abs()<1e-6);
        head=head.max(driven.system.state()[2].abs());
        pressure=pressure.max(driven.air.as_ref().unwrap().uniform_pressure(driven.system.state()).unwrap().abs());
    }
    let o=observe(&driven.system,4).unwrap();assert!(o.position_m>s.carrier.unwrap().initial_position_m);
    assert!(o.maximum_wire_slope>0.0 && work>0.0 && head>0.0 && pressure>0.0);
    assert_eq!(force_port(&driven.system,4),force_port(&manual.system,4));
}

#[test]
fn unselected_carrier_preserves_original_fixed_bank_and_finite_performance_admission() {
    let stock=SnareSet::reference(false);let (selected,carrier)=select(None,Some(stock)).unwrap();
    assert!(carrier.is_none());assert_eq!(format!("{selected:?}"),format!("{:?}",Some(stock)));
    assert_eq!(selected.unwrap().mode_count().unwrap(),160);
    let mut a=build(selected.unwrap(),false,false,false,false);let mut b=build(stock,false,false,false,false);
    let gate=CancelGate::new_clock_free();
    for _ in 0..16 {a.system.step(&a.force,&gate).unwrap();b.system.step(&b.force,&gate).unwrap();assert_eq!(a.system.state(),b.system.state());}
    let e=build(SnareSet {strands:1,modes_per_strand:2,contact_cells:4,..bank(false)},false,false,false,false);
    let input=Specification::parse(SHORT).unwrap().into_input(&e.system,3).unwrap();
    assert!(e.system.with_stick_drives(vec![input],2e-6,8,e.force.len()).is_err(),"never truncate a carrier performance to render duration");
    let wrong=build(stock,false,false,false,false);
    assert!(Specification::parse(SHORT).unwrap().into_input(&wrong.system,3).is_err());
}
