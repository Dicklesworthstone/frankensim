use super::*;
use crate::{Experiment, ImpactSystem, Mechanics, Stroke, mallets,
    shaft_playing, shell_prepare, specimen};
use fs_couple::render::plate::impact::ImpactSubstepConfig;
use fs_exec::CancelGate;

const CARD:&str=include_str!("estimated-splash-mount.fsmount");
fn specimen()->specimen::Specimen {
    let mut s=specimen::Specimen::reference();s.azimuths=8;s
}
fn inner(system:&Mechanics)->&ImpactSystem {
    match system {Mechanics::Reference(s)=>s,Mechanics::Nonlinear(s)=>s,Mechanics::Substepped(s)=>s,
        Mechanics::Driven{inner:s,..}=>inner(s),Mechanics::Prepared(_)=>panic!("mount felt is nonlinear")}
}
fn experiment(mount:Option<&Spec>,second:bool,audio:bool)->Experiment {
    crate::splash_with_mount(1024,2e-6,audio,
        Stroke {speed_m_s:0.4,position_m:Some([0.068,0.01])},Some(specimen()),&[],
        second.then_some(Stroke {speed_m_s:0.2,position_m:Some([-0.055,0.012])}),None,
        &mallets::Selection::default(),&shaft_playing::Selection::default(),None,Some(0.03),mount).unwrap()
}

#[test]
fn complete_mount_input_has_no_silent_physical_defaults() {
    let spec=Spec::parse(CARD).unwrap();assert_eq!(spec.sites.len(),3);assert_eq!(spec.creep.len(),1);
    for record in ["carrier,0.4,3200,8,0,0\n","upper,0.0280,0.006,0.05\n",
        "lower,0.0150,0.006,0.05\n","felt,30000,0.2,2.2,3,0.15,0.7\n"] {
        assert!(Spec::parse(&CARD.replace(record,"")).is_err());
        assert!(Spec::parse(&format!("{CARD}{record}")).is_err());
    }
    for bad in [CARD.replace("carrier,0.4","carrier,0"),CARD.replace("3200,8","-1,8"),
        CARD.replace("0.006,0.05","0,0.05"),CARD.replace("0.006,0.05","0.006,0.8"),
        CARD.replace("upper,0.0280","upper,0.014"),CARD.replace("0.00019137535248117825","-0.001"),
        CARD.replace("site,0.012,0,","site,NaN,0,"),format!("{CARD}creep,0,1\n"),
        format!("{CARD}site,0.012,0,0.001\n"),format!("{CARD}invented,1\n"),
        format!("{CARD}{}","creep,1,1\n".repeat(4)),format!("{CARD}{}","#".repeat(65536))] {
        assert!(Spec::parse(&bad).is_err(),"{bad}");
    }
    assert!(Spec::parse(&CARD.replace("creep,4500,24\n","")).unwrap().creep.is_empty());
    let no_sites=CARD.lines().filter(|line|!line.starts_with("site,")).collect::<Vec<_>>().join("\n");
    assert!(Spec::parse(&no_sites).is_err());
    let too_many=(0..9).fold(no_sites,|s,i|format!("{s}\nsite,{},0,0.001",0.01+i as f64*0.001));
    assert!(Spec::parse(&too_many).is_err());
}

#[test]
fn mount_selection_is_transactional_and_not_silently_ignored_by_other_instruments() {
    let mut args=vec!["splash-mic".into(),"32".into(),"--mount-spec".into(),"stand.fsmount".into()];
    assert_eq!(option(&mut args).unwrap(),Some("stand.fsmount".into()));assert_eq!(args,["splash-mic","32"]);
    for mut bad in [vec!["--mount-spec".into()],vec!["--mount-spec".into(),"--analytic-newton".into()],
        vec!["--mount-spec".into(),"a".into(),"--mount-spec".into(),"b".into()]] {
        let before=bad.clone();assert!(option(&mut bad).is_err());assert_eq!(bad,before);
    }
    for command in ["splash","splash-wav","splash-mic"] {admit_command(true,command).unwrap();}
    for command in ["drum","snare","hihat","kettle","unknown"] {
        assert!(admit_command(true,command).is_err());admit_command(false,command).unwrap();
    }
    // Actual frontend refuses scope before trying to open a nonexistent file.
    let error=crate::run_args(vec!["drum".into(),"--mount-spec".into(),"missing.fsmount".into()]).unwrap_err();
    assert!(error.to_string().contains("single supplied/reference shell"));
}

#[test]
fn finite_skin_geometry_controls_each_gap_force_row_and_local_creep() {
    let spec=Spec::parse(CARD).unwrap();
    let (shell,reduction)=shell_prepare::prepare_rocking(&specimen(),2e-6,0.03).unwrap();
    let n=3+reduction.mode_count();let mount=spec.compile(&shell,&reduction,1,n).unwrap();
    assert_eq!(mount.port.coordinate,n);assert_eq!(mount.pads.len(),6);
    let positive=reduction.surface_positions(&shell.nodal_thickness_m,ShellFace::Positive).unwrap();
    let negative=reduction.surface_positions(&shell.nodal_thickness_m,ShellFace::Negative).unwrap();
    let mut rotary_difference=0.0_f64;
    for (i,site) in spec.sites.iter().enumerate() {
        let mut skins=Vec::new();
        for (nodes,side) in [(&positive,ShellFace::Positive),(&negative,ShellFace::Negative)] {
            let (f,b)=playing::shell_location(nodes,&shell.mesh.tris,site.xy_m).unwrap();
            // Independently interpolate the offset geometry's height.
            let z=(0..3).map(|j|b[j]*nodes[shell.mesh.tris[f][j]][2]).sum::<f64>();
            let surface=reduction.surface_point_port(&shell.nodal_thickness_m,f,b,side,[0.,0.,1.]).unwrap();
            let middle=reduction.point_port(f,b,[0.,0.,1.]).unwrap();
            rotary_difference=rotary_difference.max(surface.weights.iter().zip(&middle)
                .map(|(a,b)|(a-b).abs()).fold(0.0_f64,f64::max));
            skins.push((z,surface.weights));
        }
        skins.sort_by(|a,b|b.0.total_cmp(&a.0));
        for (j,(backing,sign)) in [(spec.upper,1.0),(spec.lower,-1.0)].into_iter().enumerate() {
            let pad=&mount.pads[2*i+j];let (z,weights)=&skins[j];
            assert!((pad.precompression_m-(backing.thickness_m-sign*(backing.z_m-z))).abs()<1e-16);
            for k in 0..reduction.mode_count() {assert_eq!(pad.weights[1+k],sign*weights[k]);}
            assert_eq!(pad.weights[n],-sign*mount.port.inverse_sqrt_mass);
            assert_eq!(pad.weights[0],0.0);assert_eq!(&pad.weights[n-2..n],&[0.0,0.0]);
            assert!((pad.creep[0].stiffness_n_m-1500.0).abs()<1e-10);
            assert!((pad.creep[0].viscosity_n_s_m-8.0).abs()<1e-12);
        }
    }
    assert!(rotary_difference>1e-12,"skin-offset moments must not collapse to midsurface-only contact");
    for parity in [0,1] {
        let pads:Vec<_>=mount.pads.iter().skip(parity).step_by(2).collect();
        assert!((pads.iter().map(|p|p.area_m2).sum::<f64>()-spec.area().unwrap()).abs()<1e-18);
        assert!((pads.iter().map(|p|p.creep[0].stiffness_n_m).sum::<f64>()-4500.0).abs()<1e-10);
    }
}

#[test]
fn surveyed_thickness_changes_compression_and_invalid_mount_locations_do_not_snap() {
    let spec=Spec::parse(CARD).unwrap();let original=specimen();
    let (shell,reduction)=shell_prepare::prepare(&original,2e-6).unwrap();
    let old=spec.compile(&shell,&reduction,1,1+reduction.mode_count()).unwrap();
    let mut changed=original.clone();for station in &mut changed.stations {station.thickness_m*=1.2;}
    let (thick,reduced)=shell_prepare::prepare(&changed,2e-6).unwrap();
    let new=spec.compile(&thick,&reduced,1,1+reduced.mode_count()).unwrap();
    assert!(old.pads.iter().zip(&new.pads).all(|(a,b)|b.precompression_m>a.precompression_m));
    let mut hole=spec.clone();hole.sites[0].xy_m=[0.,0.];
    assert!(hole.compile(&shell,&reduction,1,1+reduction.mode_count()).is_err());
    let mut outside=spec.clone();outside.sites[0].xy_m=[0.2,0.];
    assert!(outside.compile(&shell,&reduction,1,1+reduction.mode_count()).is_err());
    let mut collision=spec.clone();collision.upper.z_m=0.021;
    assert!(collision.compile(&shell,&reduction,1,1+reduction.mode_count()).is_err());
    let mut displaced=spec.clone();displaced.carrier.initial_position_m=0.02;
    assert!(displaced.compile(&shell,&reduction,1,1+reduction.mode_count()).is_err());
    assert!(spec.compile(&shell,&reduction,usize::MAX,4).is_err());
    assert!(spec.compile(&shell,&reduction,1,usize::MAX).is_err());
}

#[test]
fn adding_one_support_preserves_multisite_contact_laws_and_the_old_source_prefix() {
    let contact=Obstacle::new(vec![1.,2.,3.,-4.,5.,6.],2,3,vec![0.1,-0.2],vec![0.3,0.7],
        4e7,1.5,"independent multi-point source".into()).unwrap().with_internal_loss(0.02).unwrap();
    let embedded=extend_contact(&contact,3,4).unwrap();
    assert_eq!(embedded.collocation(),&[1.,2.,3.,0.,-4.,5.,6.,0.]);
    assert_eq!(embedded.gaps(),contact.gaps());assert_eq!(embedded.weights(),contact.weights());
    assert_eq!(embedded.stiffness(),contact.stiffness());assert_eq!(embedded.alpha(),contact.alpha());
    assert_eq!(embedded.provenance(),contact.provenance());assert_eq!(embedded.internal_loss(),0.02);
    assert!(extend_contact(&contact,2,4).is_err());assert!(extend_contact(&contact,3,2).is_err());
    let mount=Spec::parse(CARD).unwrap();let old=experiment(None,true,true);let new=experiment(Some(&mount),true,true);
    assert_eq!(new.force.len(),old.force.len()+1);
    assert_eq!(new.second_stick.unwrap().coordinate,old.second_stick.unwrap().coordinate);
    assert_eq!(new.acoustics.as_ref().unwrap().state_modes(),old.acoustics.as_ref().unwrap().state_modes());
    assert_eq!(&new.observer_a[..old.force.len()],&old.observer_a);
    assert_eq!(&new.observer_b[..old.force.len()],&old.observer_b);
    assert_eq!(new.observer_a[old.force.len()],0.0);assert_eq!(new.observer_b[old.force.len()],0.0);
    assert!(inner(&new.system).felt_history(5).is_some());assert!(inner(&new.system).felt_history(6).is_none());
}

#[test]
fn one_cymbal_and_carrier_ring_together_on_the_same_transactional_clock() {
    let mut spec=Spec::parse(CARD).unwrap();spec.carrier.initial_velocity_m_s=0.002;
    let make=|| {
        let mut e=experiment(Some(&spec),false,false);
        e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
            ImpactSubstepConfig {max_depth:8,max_attempts:511}).unwrap();e
    };
    let mut e=make();let mut clean=make();let initial=inner(&e.system).stored_energy_j();
    let port=spec.port(e.force.len()).unwrap();let gate=CancelGate::new_clock_free();
    let mut loss=0.0;let mut motion=0.0_f64;let mut shell_motion=0.0_f64;
    for tick in 0..1024 {
        if tick==128 {
            let before=e.system.state().to_vec();let history=inner(&e.system).felt_history(0).unwrap();
            let samples=inner(&e.system).samples();let cancel=CancelGate::new_clock_free();cancel.request();
            assert!(e.system.step(&e.force,&cancel).is_err());
            assert_eq!(e.system.state(),before);assert_eq!(inner(&e.system).felt_history(0).unwrap(),history);
            assert_eq!(inner(&e.system).samples(),samples);
        }
        let frame=e.system.step(&e.force,&gate).unwrap();let other=clean.system.step(&clean.force,&gate).unwrap();
        loss+=frame.dissipated_energy_j;
        assert_eq!(e.system.state(),clean.system.state());assert_eq!(frame.stored_energy_j,other.stored_energy_j);
        assert!((frame.stored_energy_j+loss-initial).abs()<1e-6);
        motion=motion.max(port.motion(e.system.state()).unwrap()[0].abs());
        shell_motion=shell_motion.max(e.system.state()[2].abs());
    }
    assert!(motion>1e-9 && shell_motion>1e-10 && loss>0.0);
    let mut out=Vec::new();header(&mut out).unwrap();row(port,e.system.state(),&mut out).unwrap();
    assert!(String::from_utf8(out).unwrap().starts_with(",mount_position_m,mount_velocity_m_s,"));
    let before=e.system.state().to_vec();assert!(e.system.step(&e.force,&gate).is_err());assert_eq!(e.system.state(),before);
}

#[test]
fn mount_inertia_changes_the_actual_shell_trajectory_without_changing_shell_or_stick_materials() {
    let mut light=Spec::parse(CARD).unwrap();light.carrier.mass_kg=0.04;light.carrier.initial_velocity_m_s=0.01;
    let mut heavy=light.clone();heavy.carrier.mass_kg=4.0;
    let make=|spec:&Spec| {
        let mut e=experiment(Some(spec),false,false);
        e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
            ImpactSubstepConfig {max_depth:8,max_attempts:511}).unwrap();e
    };
    let mut a=make(&light);let mut b=make(&heavy);let gate=CancelGate::new_clock_free();
    assert_eq!(a.stick_weight,b.stick_weight);assert_eq!(a.observer_a,b.observer_a);
    for _ in 0..768 {a.system.step(&a.force,&gate).unwrap();b.system.step(&b.force,&gate).unwrap();}
    assert!((a.system.state()[2]-b.system.state()[2]).abs()>1e-11,
        "the shell, not merely support reporting, must respond to changed mounting inertia");
}

#[test]
fn shared_carrier_does_not_reindex_flexible_hands_moving_mutes_or_acoustic_sources() {
    use fs_couple::render::plate::impact::striker::{RadiusStation,flexible::FlexibleStriker};
    use fs_plate::shell::stiffened::beam::RoundBeamSpec;
    let shaft=||FlexibleStriker::new(&[(0.0,0.005),(0.4,0.005)].map(|(position_m,radius_m)|
        RadiusStation {position_m,radius_m}),RoundBeamSpec {young_pa:12e9,density_kg_m3:800.,
        pivot_m:0.1,contact_m:0.39,hand_m:0.16,subdivisions:8,maximum_hz:3000.,maximum_modes:17},0.001).unwrap();
    let shafts=shaft_playing::Selection {first:Some(shaft()),second:Some(shaft())};
    let mute=crate::compliant_mute::Spec::parse(include_str!("estimated-cymbal-mute.fsm")).unwrap();
    let spec=Spec::parse(CARD).unwrap();
    let make=|mount|crate::splash_with_mount(256,2e-6,true,
        Stroke {speed_m_s:0.4,position_m:Some([0.068,0.01])},Some(specimen()),&[],
        Some(Stroke {speed_m_s:0.2,position_m:Some([-0.055,0.012])}),Some(&mute),
        &mallets::Selection::default(),&shafts,None,Some(0.03),mount).unwrap();
    let old=make(None);let new=make(Some(&spec));let at=old.force.len();
    assert_eq!(new.force.len(),at+1);assert_eq!(spec.port(new.force.len()).unwrap().coordinate,at);
    for (a,b) in old.flexible_sticks.iter().zip(&new.flexible_sticks) {
        let (a,b)=(a.as_ref().unwrap(),b.as_ref().unwrap());
        assert_eq!(a.rigid_coordinate(),b.rigid_coordinate());assert_eq!(a.elastic_start(),b.elastic_start());
        assert_eq!(a.elastic_modes(),b.elastic_modes());
        let old_row=a.hand_row(at).unwrap();let new_row=b.hand_row(at+1).unwrap();
        assert_eq!(&new_row[..at],&old_row);assert_eq!(new_row[at],0.0);
        assert_eq!(b.tip_row(at+1).unwrap()[at],0.0);
    }
    assert_eq!(old.acoustics.as_ref().unwrap().state_modes(),new.acoustics.as_ref().unwrap().state_modes());
    let before=mute.clone().into_inputs(old.mute.as_ref().unwrap()).unwrap();
    let after=mute.into_inputs(new.mute.as_ref().unwrap()).unwrap();
    for (a,b) in before.iter().zip(&after) {assert_eq!(a.coordinate,b.coordinate);assert_eq!(a.tip_weight,b.tip_weight);}
    assert!(inner(&new.system).felt_history(13).is_some());assert!(inner(&new.system).felt_history(14).is_none());
    let mut diagnostics=Vec::new();new.mute.as_ref().unwrap().row(&new.system,&mut diagnostics).unwrap();
    assert!(!diagnostics.is_empty());
}
