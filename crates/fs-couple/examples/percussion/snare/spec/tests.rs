use super::*;
use super::super::super::{drum_with_radiation, drum_spec, Mechanics, Stroke, Experiment,
    acoustics, nonlinear_snare, mechanics};
use fs_couple::render::plate::impact::ImpactSubstepConfig;
use fs_exec::CancelGate;

const REFERENCE: &str = include_str!("../../estimated-snare.fsn");
fn supplied(ea: f64) -> SnareSet {
    Specification::parse(&REFERENCE.replace("stretching,off", &format!("stretching,{ea},0.2"))).unwrap().engaged
}
fn build(s: SnareSet, heads: bool, air: bool, second: bool, audio: bool) -> Experiment {
    drum_with_radiation(512, if audio {acoustics::MECHANICAL_DT}else{2e-6}, audio,
        !(heads || s.stretching.is_some()), Some(s), heads,
        Stroke {speed_m_s:0.8,position_m:Some([0.06,0.01])}, air, None,
        Some(drum_spec::Spec {radial_intervals:2,azimuths:8,..drum_spec::Spec::reference()}),
        second.then_some(Stroke {speed_m_s:0.2,position_m:Some([-0.05,0.02])}),
        &[], if air {25.0}else{0.0}, None, false).unwrap()
}
fn energy(s: &Mechanics) -> f64 {
    match s {
        Mechanics::Reference(s)=>s.stored_energy_j(), Mechanics::Nonlinear(s)=>s.stored_energy_j(),
        Mechanics::Substepped(s)=>s.stored_energy_j(), Mechanics::Driven {inner,..}=>energy(inner),
        Mechanics::Prepared(_)=>panic!("this check requires the nonlinear image"),
    }
}

#[test]
fn complete_si_file_reproduces_reference_parameters_and_explicit_disengagement() {
    let parsed=Specification::parse(REFERENCE).unwrap();
    for off in [false,true] {assert_eq!(format!("{:?}",parsed.select(off)),format!("{:?}",SnareSet::reference(off)));}
    let custom=Specification::parse(&REFERENCE.replace("0.00002,0.003", "-0.000001,0.005")).unwrap();
    assert_eq!(custom.select(false).clearance_m,-0.000001);
    assert_eq!(custom.select(true).clearance_m,0.005);
    let p=supplied(200.0);assert_eq!(p.stretching.unwrap(),StringStretching {axial_rigidity_n:200.0,maximum_slope:0.2});
    assert_eq!(p.coil.linear_density_kg_m().unwrap(),parsed.engaged.coil.linear_density_kg_m().unwrap());
    assert_eq!(p.mode_count().unwrap(),160);
    let ordered: Vec<_>=REFERENCE.lines().filter(|r| !r.trim().starts_with('#') && !r.trim().is_empty()).collect();
    let reordered=format!("{}\n{}\n",ordered[0],ordered[1..].iter().rev().copied().collect::<Vec<_>>().join("\n"));
    assert_eq!(format!("{:?}",Specification::parse(&reordered).unwrap()),format!("{parsed:?}"));
}

#[test]
fn malformed_or_incomplete_files_refuse_before_geometry_and_do_not_acquire_defaults() {
    for (from,to) in [("bank,20,8,12","bank,25,8,12"),("bank,20,8,12","bank,20,9999999999999999999999999,12"),
        ("bank,20,8,12","bank,20,8,4"),("coil,0.00015","coil,-0.00015"),
        ("mechanics,0.7","mechanics,NaN"),("mechanics,0.7","mechanics,inf"),
        ("mechanics,0.7,0.000001,4","mechanics,0.7,-1,4"),
        ("0.00002,0.003","0.003,0.00002"),("500000000,1.5,0.05","500000000,0.5,0.05"),
        ("500000000,1.5,0.05","500000000,1.5,-0.05"),
        ("stretching,off","stretching,-1,0.2"),("stretching,off","stretching,10,0.4"),
        ("stretching,off","stretching,10,0"),("stretching,off","stretching,10"),
        ("stretching,off","stretching,off,1"),("stretching,off",""),
        ("frankensim-snare-spec-v1","frankensim-snare-spec-v2")] {
        assert!(Specification::parse(&REFERENCE.replace(from,to)).is_err(),"accepted {to}");
    }
    for row in ["bank,20,8,12,0.30,0.04","unknown,1","coil,0.1","stretching,off"] {
        assert!(Specification::parse(&format!("{REFERENCE}\n{row}\n")).is_err());
    }
    assert!(Specification::parse(&" ".repeat(MAX_BYTES as usize+1)).is_err());
    assert!(select(Some("must-not-be-opened"),"drum").is_err());
    let path=format!("{}/examples/percussion/estimated-snare.fsn",env!("CARGO_MANIFEST_DIR"));
    for c in ["snare","snare-off","snare-wav","snare-off-wav","snare-mic","snare-off-mic"] {
        assert!(select(Some(&path),c).unwrap().is_some());
        nonlinear_snare::admit_prepared_command(true,true,c).unwrap();
    }
    let mut args=vec!["snare".into(),"--snare-spec".into(),path.clone(),"--head-stretching".into()];
    assert_eq!(option(&mut args).unwrap(),Some(path));assert_eq!(args,["snare","--head-stretching"]);
    assert!(option(&mut vec!["--snare-spec".into()]).is_err());
    assert!(option(&mut vec!["--snare-spec".into(),"--analytic-newton".into()]).is_err());
    assert!(option(&mut vec!["--snare-spec".into(),"a".into(),"--snare-spec".into(),"b".into()]).is_err());
}

#[test]
fn file_selected_linear_bank_preserves_the_original_joint_trajectory() {
    let supplied=Specification::parse(REFERENCE).unwrap().engaged;
    let mut file=build(supplied,false,false,false,false);
    let mut stock=build(SnareSet::reference(false),false,false,false,false);
    assert!(matches!(&file.system,Mechanics::Prepared(_)));
    assert_eq!(file.observer_a,stock.observer_a);assert_eq!(file.observer_b,stock.observer_b);
    let gate=CancelGate::new_clock_free();
    for _ in 0..160 {
        let a=file.system.step(&file.force,&gate).unwrap();let b=stock.system.step(&stock.force,&gate).unwrap();
        assert_eq!(file.system.state(),stock.system.state());assert_eq!(a.time_s,b.time_s);
        assert_eq!(a.stored_energy_j,b.stored_energy_j);assert_eq!(a.dissipated_energy_j,b.dissipated_energy_j);
    }
}

#[test]
fn whole_bank_retains_original_addresses_and_radiating_heads_for_either_head_model() {
    for heads in [false,true] {
        let s=supplied(100.0);let e=build(s,heads,true,true,true);
        assert!(matches!(&e.system,Mechanics::Reference(_)));
        assert_eq!(e.system.membrane_observation(1).is_some(),heads);
        assert_eq!(e.system.membrane_observation(2).is_some(),heads);
        let first_wire=e.second_stick.unwrap().coordinate+1;
        assert_eq!(e.air.as_ref().unwrap().coupling.structural_modes(),first_wire+160);
        for body in 4..24 {
            let o=super::super::observe(&e.system,body).unwrap();
            assert_eq!(o.tension_n,s.tension_per_strand_n);assert_eq!(o.stretching_energy_j,0.0);
        }
        assert!(super::super::observe(&e.system,24).is_none()); // appended acoustic inertia is not a wire
        assert!(e.observer_a[first_wire..].iter().chain(&e.observer_b[first_wire..]).all(|v|*v==0.0));
        assert!(e.acoustics.as_ref().unwrap().state_modes().iter().all(|i| *i<first_wire));
        assert!(e.force.len()<=256);
        assert!(drum_with_radiation(2,2e-6,false,true,Some(s),heads,Stroke::default(),false,
            None,None,None,&[],0.0,None,false).is_err(),"linear image cannot erase wire stretching");
    }
}

#[test]
fn supplied_stretching_wires_change_real_head_contact_with_cavity_drive_and_atomic_retry() {
    // Explicit small installation interference puts this runtime fixture into
    // contact immediately; it is not a claim about the reference installation.
    let small=SnareSet {strands:2,modes_per_strand:2,contact_cells:4,clearance_m:-2e-6,..supplied(10000.0)};
    let prepare=|s| {
        let mut e=build(s,true,true,true,false);
        e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
            ImpactSubstepConfig {max_depth:4,max_attempts:31}).unwrap();
        let other=e.second_stick.unwrap();
        let inputs=vec![mechanics::drive::Input {coordinate:0,tip_weight:e.stick_weight,
            program:mechanics::drive::Program::parse("0,0\n0.000256,0.02\n0.001024,0").unwrap()},
            mechanics::drive::Input {coordinate:other.coordinate,tip_weight:other.weight,
            program:mechanics::drive::Program::parse("0,0\n0.000256,-0.01\n0.001024,0").unwrap()}];
        e.system=e.system.with_stick_drives(inputs,2e-6,512,e.force.len()).unwrap();e
    };
    let mut e=prepare(small);let mut retry=prepare(small);
    let mut no_stretch=prepare(SnareSet {stretching:Some(StringStretching {axial_rigidity_n:0.0,..small.stretching.unwrap()}),..small});
    assert_eq!(e.system.state(),no_stretch.system.state());
    let gate=CancelGate::new_clock_free();let cancelled=CancelGate::new_clock_free();cancelled.request();
    let initial=energy(&e.system);let mut net=0.0;let mut peak=small.tension_per_strand_n;let mut change=0.0_f64;
    let mut stretch=0.0_f64;let mut drive_work=0.0;
    for tick in 0..128 {
        if tick==64 {let before=retry.system.state().to_vec();assert!(retry.system.step(&retry.force,&cancelled).is_err());assert_eq!(retry.system.state(),before);}
        let f=e.system.step(&e.force,&gate).unwrap();let retried=retry.system.step(&retry.force,&gate).unwrap();
        no_stretch.system.step(&no_stretch.force,&gate).unwrap();
        assert_eq!(e.system.state(),retry.system.state());assert_eq!(f.time_s,retried.time_s);
        net+=f.supplied_work_j-f.dissipated_energy_j;drive_work+=f.supplied_work_j.abs();
        assert!((f.stored_energy_j-initial-net).abs()<1e-6);
        for body in 4..6 {let o=super::super::observe(&e.system,body).unwrap();peak=peak.max(o.tension_n);stretch=stretch.max(o.stretching_energy_j);}
        change=change.max(e.system.state()[2..].iter().zip(&no_stretch.system.state()[2..])
            .map(|(a,b)|(a-b).abs()).fold(0.0,f64::max));
    }
    assert!(peak>small.tension_per_strand_n && stretch>0.0 && change>1e-14 && drive_work>0.0);
}
