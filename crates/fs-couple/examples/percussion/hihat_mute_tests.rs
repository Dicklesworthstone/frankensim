use super::*;
use fs_couple::render::plate::impact::ImpactSubstepConfig;
use fs_couple::render::plate::impact::striker::RadiusStation;
use fs_plate::shell::stiffened::beam::RoundBeamSpec;

const GAS:&str="frankensim-squeeze-film-v1\nannulus,0.04,0.085,2,4\nviscosity_pa_s,0.000018\nlimits,0.000001,0.04,200000\nboundary,sealed,sealed\nisothermal,100000,293.15,287.05\n";
fn card(lower:bool)->String {
    let (surface,side)=if lower {("lower","below")}else{("upper","above")};
    format!("frankensim-compliant-mute-v1\nsurface,{surface}\n\
        site,0.06,0.01,0.0001\nsite,0.065,0.013,0.0001\n\
        jaw,{side},0.025,0.1,0,0.0005,30000,0.2,2.2,3,0.15,0.7,0\n\
        creep,{side},4000,12\nforce,{side},0,0\nforce,{side},0.0001,0.5\n\
        force,{side},0.0003,0.5\nforce,{side},0.0004,-2\n\
        force,{side},0.0006,-2\nforce,{side},0.0007,0\n")
}
fn pair_spec()->Spec {Spec::parse(include_str!("estimated-hihat.fshh")).unwrap()}
fn shell()->specimen::Specimen {let mut s=specimen::Specimen::reference();s.azimuths=8;s}
fn rest()->Stroke {Stroke{position_m:Some([-0.06,0.01]),speed_m_s:0.}}
fn shaft()->FlexibleStriker {
    FlexibleStriker::new(&[(0.0,0.005),(0.4,0.005)].map(|(position_m,radius_m)|
        RadiusStation{position_m,radius_m}),RoundBeamSpec{young_pa:12e9,density_kg_m3:800.,
        pivot_m:0.1,contact_m:0.39,hand_m:0.16,subdivisions:8,maximum_hz:3000.,maximum_modes:17},0.001).unwrap()
}
fn inner(m:&Mechanics)->&ImpactSystem {
    match m {Mechanics::Reference(s)=>s,Mechanics::Nonlinear(s)=>s,
        Mechanics::Substepped(s)=>s,Mechanics::Driven{inner:s,..}=>inner(s),
        Mechanics::Prepared(_)=>panic!("moving felt cannot discard material history")}
}

#[test]
fn paired_mute_requires_an_explicit_exterior_surface_and_retains_all_mounts() {
    let s=shell();let tips=mallets::Selection::default();let shafts=shaft_playing::Selection::default();
    for lower in [false,true] {
        let text=card(lower);let mute=compliant_mute::Spec::parse(&text).unwrap();
        for command in ["hihat","hihat-wav","hihat-mic"] {mute.admit_command(command).unwrap();}
        for command in ["splash","snare","drum"] {assert!(mute.admit_command(command).is_err());}
        assert_eq!(mute.lower_shell().unwrap(),lower);
        assert!(compliant_mute::Spec::parse(&text.replace(if lower {"below"}else{"above"},
            if lower {"above"}else{"below"})).is_err());
        let pair=build_with_mute(&pair_spec(),&s,&s,rest(),None,512,2e-6,true,None,
            &shafts,&tips,None,Some(&mute)).unwrap();
        let e=pair.experiment;let n=e.force.len();let jaw=e.mute.as_ref().unwrap().ports[0];
        assert_eq!(jaw.coordinate,pair.pedal.coordinate+1);
        assert_eq!(e.system.state().len(),2*n+14);
        for i in 0..14 {assert!(inner(&e.system).felt_history(i).is_some());}
        assert!(inner(&e.system).felt_history(14).is_none());
        for row in pair.collision.collocation().chunks(n) {assert_eq!(row[jaw.coordinate],0.);}
        assert_eq!(e.acoustics.as_ref().unwrap().state_modes(),
            pair.upper_modes.chain(pair.lower_modes).collect::<Vec<_>>());
    }
    let mute=compliant_mute::Spec::parse(&card(false)).unwrap();
    let tip=mallets::Spec::parse("frankensim-felt-mallet-v1\ngeometry,0.02,0.003,0.006,0.00002\nfelt,1000000,0.2,2.2,3,0.15,0.7\nconditioning,0\n").unwrap();
    let tips=mallets::Selection{first:Some(tip),second:None};
    let error=build_with_mute(&pair_spec(),&s,&s,rest(),None,1,2e-6,false,None,
        &shafts,&tips,None,Some(&mute)).err().unwrap().to_string();
    assert!(error.contains("sixteen-pad limit"));
    let single=compliant_mute::Spec::parse(&card(false).replace("surface,upper","surface,shell")).unwrap();
    assert!(single.admit_command("hihat").is_err());
}

#[test]
fn both_outer_skin_choke_rows_apply_true_force_moments_at_world_stations() {
    let s=Shell::new(&shell(),2e-6).unwrap();let count=s.reduction.mode_count();
    for lower in [false,true] {
        let mute=compliant_mute::Spec::parse(&card(lower)).unwrap();
        let start=if lower {1+count}else{1};let structural=2+2*count;
        let pads=mute.paired_shell(&s.reduction,&s.mesh.nodal_thickness_m,&s.outer_nodes,
            &s.mesh.mesh.tris,start,structural).unwrap();
        assert_eq!(pads.ports[0].coordinate,structural);
        let mut moment_difference=0.0_f64;
        for (pad,p) in pads.pads.iter().zip([[0.06,0.01],[0.065,0.013]]) {
            let local=if lower {[p[0],-p[1]]}else{p};
            let port=s.port(local,ShellFace::Positive).unwrap();
            assert!((port.position_m[0]-local[0]).abs()<1e-12);
            assert!((port.position_m[1]-local[1]).abs()<1e-12);
            // On either exterior, outward shell motion increases compression.
            assert_eq!(&pad.weights[start..start+count],port.weights.as_slice());
            assert!(pad.weights[..start].iter().all(|v|*v==0.));
            assert!(pad.weights[start+count..structural].iter().all(|v|*v==0.));
            assert!((pad.weights[start]-1./s.mesh.mass_kg.sqrt()).abs()<1e-12);
            let (f,b)=playing::shell_location(&s.outer_nodes,&s.mesh.mesh.tris,local).unwrap();
            let middle=s.reduction.point_port(f,b,[0.,0.,1.]).unwrap();
            for (a,b) in port.weights.iter().zip(middle) {moment_difference=moment_difference.max((a-b).abs());}
        }
        assert!(moment_difference>1e-12,"finite-thickness outer contact must retain rotation work");
    }
}

#[test]
fn lower_choke_moves_real_shells_and_retracts_with_one_pedal_hand_and_material_clock() {
    let make=|| {
        let mute=compliant_mute::Spec::parse(&card(true)).unwrap();let s=shell();
        let shafts=shaft_playing::Selection{first:Some(shaft()),second:None};
        let film=squeeze::Config::parse(GAS).unwrap();
        let pair=build_with_mute(&pair_spec(),&s,&s,rest(),None,512,2e-6,false,Some(&film),
            &shafts,&mallets::Selection::default(),None,Some(&mute)).unwrap();
        let mut e=pair.experiment;let n=e.force.len();let a=e.mute.as_ref().unwrap();
        assert_eq!(e.flexible_sticks[0].as_ref().unwrap().elastic_start(),a.ports[0].coordinate+1);
        assert_eq!(e.system.state().len(),2*n+14+8);
        let mut inputs=mute.into_inputs(a).unwrap();
        inputs.push(mechanics::drive::Input{program:mechanics::drive::Program::parse(
            "0,0\n0.0002,0.01\n0.0007,0").unwrap(),coordinate:pair.pedal.coordinate,tip_weight:pair.pedal.weight});
        let hand=mechanics::drive::SpatialInput{program:mechanics::drive::Program::parse(
            "0,0\n0.0002,-0.01\n0.0007,0").unwrap(),weights:e.flexible_sticks[0].as_ref().unwrap().hand_row(n).unwrap()};
        e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
            ImpactSubstepConfig{max_depth:8,max_attempts:511}).unwrap()
            .with_player_drives(inputs,vec![hand],2e-6,512,n).unwrap();e
    };
    let mut e=make();let mut clean=make();let initial=inner(&e.system).stored_energy_j();
    let mut net=0.;let mut peak_contact=0.0_f64;let mut peak_lower=0.0_f64;
    let gate=CancelGate::new_clock_free();
    for tick in 0..512 {
        if tick==128 {
            let before=e.system.state().to_vec();let history=(0..14)
                .map(|i|inner(&e.system).felt_history(i)).collect::<Vec<_>>();
            let cancel=CancelGate::new_clock_free();cancel.request();
            assert!(e.system.step(&e.force,&cancel).is_err());
            let mut bad=e.force.clone();bad[0]=f64::NAN;
            assert!(e.system.step(&bad,&gate).is_err());assert_eq!(e.system.state(),before);
            assert_eq!((0..14).map(|i|inner(&e.system).felt_history(i)).collect::<Vec<_>>(),history);
        }
        let f=e.system.step(&e.force,&gate).unwrap();clean.system.step(&clean.force,&gate).unwrap();
        assert_eq!(e.system.state(),clean.system.state());assert_eq!(f.time_s,(tick+1) as f64*2e-6);
        net+=f.supplied_work_j-f.dissipated_energy_j;
        assert!((f.stored_energy_j-initial-net).abs()<1e-6);
        peak_contact=peak_contact.max((12..14).map(|i|inner(&e.system).felt_observation(i).unwrap().1).sum());
        peak_lower=peak_lower.max(e.observer_b.iter().enumerate().map(|(i,b)|b*e.system.state()[2*i]).sum::<f64>().abs());
    }
    assert!(peak_contact>0. && peak_lower>0.);
    assert!((12..14).any(|i|inner(&e.system).felt_history(i).unwrap().eps_max>0.));
    for i in 12..14 {assert_eq!(inner(&e.system).felt_observation(i).unwrap().1,0.,"retraction must remove actual compression");}
    let jaw=e.mute.as_ref().unwrap().ports[0];assert!(e.system.state()[2*jaw.coordinate]<0.);
}
