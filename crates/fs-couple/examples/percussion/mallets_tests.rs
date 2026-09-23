use super::*;
use super::super::{drum_spec,drum_with_mallets,drum_with_material,snare::SnareSet,
    Experiment,Mechanics,mechanics::drive,acoustics};
use fs_couple::render::plate::impact::{ImpactSystem,ImpactSubstepConfig};
use fs_exec::CancelGate;

const CARD:&str="frankensim-felt-mallet-v1\ngeometry,0.02,0.012,0.006,0.00002\nfelt,100000,0.2,2.2,3,0.15,0.7\nconditioning,0\n";
fn spec()->Spec {Spec::parse(CARD).unwrap()}
fn first()->Stroke {Stroke {position_m:Some([0.06,0.01]),speed_m_s:0.5}}
fn second()->Stroke {Stroke {position_m:Some([-0.05,0.02]),speed_m_s:0.3}}
fn drum()->drum_spec::Spec {drum_spec::Spec {radial_intervals:2,azimuths:8,..drum_spec::Spec::reference()}}
fn inner(s:&Mechanics)->&ImpactSystem {
    match s {Mechanics::Reference(s)=>s,Mechanics::Nonlinear(s)=>s,
        Mechanics::Substepped(s)=>s,Mechanics::Driven{inner:s,..}=>inner(s),
        Mechanics::Prepared(_)=>panic!("felt cannot be discarded by the linear-only image")}
}
fn build(selected:&Selection,wires:Option<SnareSet>,air:bool,audio:bool)->Experiment {
    drum_with_mallets(1024,if audio {acoustics::MECHANICAL_DT}else{2e-6},audio,
        false,wires,false,first(),air,None,Some(drum()),Some(second()),&[],0.0,None,false,None,selected).unwrap()
}

#[test]
fn mallet_file_has_no_material_or_geometry_defaults_and_admits_only_physical_commands() {
    let s=spec();assert_eq!(s.jaw.mass_kg,0.02);assert!(s.jaw.creep.is_empty());
    for record in ["geometry,0.02,0.012,0.006,0.00002\n","felt,100000,0.2,2.2,3,0.15,0.7\n","conditioning,0\n"] {
        assert!(Spec::parse(&CARD.replace(record,"")).is_err());
        assert!(Spec::parse(&format!("{CARD}{record}")).is_err());
    }
    for text in [CARD.replace("0.02,0.012","0,0.012"),CARD.replace("0.012","-0.012"),
        CARD.replace("0.006","NaN"),format!("{CARD}creep,0,1\n"),format!("{CARD}unknown,1\n")] {
        assert!(Spec::parse(&text).is_err());
    }
    assert!(Spec::parse(&format!("{CARD}creep,3000,6\n")).is_ok());
    let selected=Selection {first:Some(s),second:None};
    for c in ["drum","drum-mic","drum-stretch-wav","snare","snare-off-mic"] {
        selected.admit(c,first(),None).unwrap();
    }
    for c in ["splash","drum-modal","unknown"] {assert!(selected.admit(c,first(),None).is_err());}
    assert!(selected.admit("drum",Stroke::default(),None).is_err());
    let mut args=vec!["snare".into(),"--second-mallet-spec".into(),"right.fsmallet".into(),"32".into()];
    let paths=options(&mut args).unwrap();assert_eq!(paths,[None,Some("right.fsmallet".into())]);
    assert_eq!(args,["snare","32"]);
    let selected=Selection {first:None,second:Some(spec())};
    assert!(selected.admit("snare",first(),None).is_err());
}

#[test]
fn footprint_preserves_disk_area_and_second_moments_and_refuses_unsampled_rim_crossing() {
    let s=spec();let a=std::f64::consts::PI*s.radius_m.powi(2)/4.0;
    let points=s.points([0.0,0.0]);
    assert!((points.iter().map(|p|a*p[0]*p[0]).sum::<f64>()-std::f64::consts::PI*s.radius_m.powi(4)/4.0).abs()<1e-20);
    assert_eq!(points.iter().map(|p|p[0]*p[1]).sum::<f64>(),0.0);
    let (films,modes)=drum().prepare(2e-6,false).unwrap();let n=1+modes.iter().map(Vec::len).sum::<usize>();
    let tip=s.compile(&films[0],&modes[0],first(),0,n).unwrap();
    assert_eq!(tip.pads.len(),4);
    assert!((tip.pads.iter().map(|p|p.area_m2).sum::<f64>()-4.0*a).abs()<1e-18);
    assert!(tip.pads.iter().all(|p|p.weights[1+modes[0].len()..].iter().all(|v|*v==0.0)));
    // The quadrature ring fits inside the circumcircle, but the physical disk
    // crosses the actual polygonal bearing edge and must not be silently clipped.
    let r=films[0].spec.radius_m;let angle=std::f64::consts::PI/8.0;
    let center_radius=r*angle.cos()-0.9*s.radius_m;
    let stroke=Stroke {position_m:Some([center_radius*angle.cos(),center_radius*angle.sin()]),speed_m_s:0.0};
    assert!(s.compile(&films[0],&modes[0],stroke,0,n).is_err());
}

#[test]
fn two_felt_tips_preserve_complete_snare_and_acoustic_addresses_without_duplicate_hertz_contacts() {
    let selected=Selection {first:Some(spec()),second:Some(spec())};
    let e=build(&selected,Some(SnareSet::reference(false)),true,true);
    let s=inner(&e.system);let port=e.second_stick.unwrap();
    for i in 0..8 {assert!(s.felt_history(i).is_some());}
    assert!(s.felt_history(8).is_none());
    assert_eq!(e.air.as_ref().unwrap().coupling.structural_modes(),port.coordinate+1+160);
    assert_eq!(e.acoustics.as_ref().unwrap().state_modes(),&(1..port.coordinate).collect::<Vec<_>>());
    assert!(e.observer_a[port.coordinate..].iter().all(|v|*v==0.0));
    assert!((e.system.state()[0]*e.stick_weight+spec().jaw.initial_gap_m).abs()<1e-18);
    // The dynamic test below exercises actual felt/head/air exchange.
    assert_eq!(e.force.len(),e.air.as_ref().unwrap().coupling.total_modes());
}

#[test]
fn real_mallets_transfer_motion_and_keep_both_force_programs_and_felt_history_on_retry() {
    let make=|| {
        let selected=Selection {first:Some(Spec::parse(&format!("{CARD}creep,3000,6\n")).unwrap()),second:Some(spec())};
        let mut e=build(&selected,None,true,false);let port=e.second_stick.unwrap();
        e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
            ImpactSubstepConfig {max_depth:6,max_attempts:127}).unwrap();
        let programs=vec![drive::Input {program:drive::Program::parse("0,0\n0.0002,0.05\n0.001,0").unwrap(),coordinate:0,tip_weight:e.stick_weight},
            drive::Input {program:drive::Program::parse("0,0\n0.0004,-0.02\n0.001,0").unwrap(),coordinate:port.coordinate,tip_weight:port.weight}];
        e.system=e.system.with_stick_drives(programs,2e-6,1024,e.force.len()).unwrap();e
    };
    let mut e=make();let mut clean=make();let initial=inner(&e.system).stored_energy_j();
    let gate=CancelGate::new_clock_free();let mut net=0.0;let mut crush=0.0;let mut pressure=0.0_f64;
    for tick in 0..512 {
        if tick==128 {
            let before=e.system.state().to_vec();let h=inner(&e.system).felt_history(0).unwrap();
            let cancel=CancelGate::new_clock_free();cancel.request();
            assert!(e.system.step(&e.force,&cancel).is_err());
            assert_eq!(e.system.state(),before);assert_eq!(inner(&e.system).felt_history(0).unwrap(),h);
        }
        let f=e.system.step(&e.force,&gate).unwrap();clean.system.step(&clean.force,&gate).unwrap();
        assert_eq!(e.system.state(),clean.system.state());assert_eq!(f.time_s,(tick+1) as f64*2e-6);
        net+=f.supplied_work_j-f.dissipated_energy_j;crush+=f.felt_crush_loss_j;
        assert!((f.stored_energy_j-initial-net).abs()<1e-6);
        pressure=pressure.max(e.air.as_ref().unwrap().uniform_pressure(e.system.state()).unwrap().abs());
    }
    assert!(crush>0.0 && pressure>0.0 && inner(&e.system).felt_history(0).unwrap().eps_max>0.0);
}

#[test]
fn no_mallet_selection_preserves_the_original_trajectory_and_linear_image_refuses_felt() {
    let empty=Selection::default();let mut a=build(&empty,None,false,false);
    let mut b=drum_with_material(1024,2e-6,false,false,None,false,first(),false,None,
        Some(drum()),Some(second()),&[],0.0,None,false,None).unwrap();
    let gate=CancelGate::new_clock_free();
    for _ in 0..64 {a.system.step(&a.force,&gate).unwrap();b.system.step(&b.force,&gate).unwrap();
        assert_eq!(a.system.state(),b.system.state());}
    let selected=Selection {first:Some(spec()),second:None};
    assert!(drum_with_mallets(1,2e-6,false,true,None,false,first(),false,None,Some(drum()),
        None,&[],0.0,None,false,None,&selected).is_err());
}
