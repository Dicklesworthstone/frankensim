use super::*;
use super::super::Selection;
use super::super::super::{shell_prepare,specimen,splash_with_mallets,splash_with_compliant_mute,
    Mechanics,mechanics::drive,muffling,Experiment};
use fs_couple::render::plate::impact::{ImpactSystem,ImpactSubstepConfig};
use fs_exec::CancelGate;

const CARD:&str="frankensim-felt-mallet-v1\ngeometry,0.02,0.003,0.006,0.00002\nfelt,1000000,0.2,2.2,3,0.15,0.7\nconditioning,0\n";
fn spec()->Spec{Spec::parse(CARD).unwrap()}
fn specimen()->specimen::Specimen {let mut s=specimen::Specimen::reference();s.azimuths=8;s}
fn first()->Stroke{Stroke{speed_m_s:0.5,position_m:Some([0.06,0.01])}}
fn second()->Stroke{Stroke{speed_m_s:0.3,position_m:Some([-0.05,0.02])}}
fn inner(m:&Mechanics)->&ImpactSystem {
    match m {Mechanics::Reference(s)|Mechanics::Nonlinear(s)|Mechanics::Substepped(s)=>s,
        Mechanics::Driven{inner:s,..}=>inner(s),Mechanics::Prepared(_)=>panic!("felt-capable image required")}
}
fn build(selection:&Selection,audio:bool)->Experiment {
    splash_with_mallets(1024,2e-6,audio,first(),Some(specimen()),&[
        muffling::Muffler{surface:muffling::Surface::Shell,position_m:[0.075,0.],resistance_n_s_m:0.2}
    ],Some(second()),None,selection).unwrap()
}
#[test]
fn flat_and_sloping_contact_planes_use_the_whole_disk_not_four_samples() {
    let mut nodes=[[-1.,-1.,0.],[1.,-1.,0.],[1.,1.,0.],[-1.,1.,0.]];
    let faces=[[0,1,2],[0,2,3]];let center=[0.1,0.2];let radius=0.3;
    assert_eq!(contact_plane(&nodes,&faces,center,radius).unwrap(),0.);
    for p in &mut nodes {p[2]=1.+0.3*p[0]-0.4*p[1];}
    let expected=1.+0.3*center[0]-0.4*center[1]+0.5*radius;
    let h=contact_plane(&nodes,&faces,center,radius).unwrap();assert!((h-expected).abs()<1e-14);
    let sample_high=(0..4).map(|i| {
        let a=std::f64::consts::FRAC_PI_2*i as f64;
        1.+0.3*(center[0]+radius/2_f64.sqrt()*a.cos())-0.4*(center[1]+radius/2_f64.sqrt()*a.sin())
    }).fold(f64::NEG_INFINITY,f64::max);
    assert!(h>sample_high+0.01);
    // The support point is outside this triangle: an edge/circle intersection wins.
    let t=[[-2.,0.,-2.],[2.,0.,2.],[0.,3.,-6.]];
    assert!((triangle_peak(t,[0.,0.],1.).unwrap().unwrap()-1.).abs()<1e-14);
}
#[test]
fn entire_disk_refuses_hidden_holes_rims_folds_and_multiple_coverage() {
    let mut nodes=Vec::new();for r in [0.2,2.] {for p in [[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.]] {
        nodes.push([r*p[0],r*p[1],0.]);
    }}
    let mut faces=Vec::new();for i in 0..4 {let j=(i+1)%4;faces.push([i,i+4,j+4]);faces.push([i,j+4,j]);}
    let c=[0.7,0.];let r=0.51;let d=r/std::f64::consts::SQRT_2;
    for q in [[c[0]+d,0.],[c[0]-d,0.],[c[0],d],[c[0],-d]] {
        playing::shell_location(&nodes,&faces,q).unwrap(); // all four sites miss the hole
    }
    assert!(contact_plane(&nodes,&faces,c,r).is_err());
    assert!(contact_plane(&nodes,&faces,[1.8,0.],0.3).is_err());
    assert!(contact_plane(&nodes,&faces,[0.,0.],0.1).is_err());
    assert!(contact_plane(&nodes,&faces,[0.8,0.],0.2).is_ok());
    let mut flipped=faces.clone();flipped[0].swap(1,2);
    assert!(contact_plane(&nodes,&flipped,[0.8,0.],0.2).is_err());
    let duplicated=faces.iter().chain(&faces).copied().collect::<Vec<_>>();
    assert!(contact_plane(&nodes,&duplicated,[0.8,0.],0.2).is_err());
}
#[test]
fn curved_mallet_uses_actual_skin_rotational_rows_and_separate_nonnegative_gaps() {
    let (shell,r)=shell_prepare::prepare(&specimen(),2e-6).unwrap();let s=spec();
    let n=1+r.mode_count();let tip=s.compile_shell(&r,&shell,first(),0,1,n).unwrap();
    let nodes=r.surface_positions(&shell.nodal_thickness_m,ShellFace::Positive).unwrap();
    let skin=r.radiation_surface(&shell.nodal_thickness_m,
        fs_plate::shell::reduction::radiation::RadiationSurfaceBudget{max_panels:2048,max_panel_modes:65536}).unwrap();
    for (i,t) in shell.mesh.tris.iter().enumerate(){for j in 0..3 {assert_eq!(nodes[t[j]],skin.triangles()[2*i][j]);}}
    let height=contact_plane(&nodes,&shell.mesh.tris,first().position_m.unwrap(),s.radius_m).unwrap();
    for (pad,q) in tip.pads.iter().zip(s.points(first().position_m.unwrap())) {
        let (f,b)=playing::shell_location(&nodes,&shell.mesh.tris,q).unwrap();
        let p=r.surface_point_port(&shell.nodal_thickness_m,f,b,ShellFace::Positive,[0.,0.,1.]).unwrap();
        assert_eq!(&pad.weights[1..],p.weights.as_slice());
        assert!((pad.precompression_m+p.position_m[2].mul_add(-1.,height)).abs()<1e-14);
        assert!(pad.precompression_m<=0.);
        assert!((pad.weights[0]*tip.body.initial[0].displacement_m_sqrt_kg+pad.precompression_m)<=-s.jaw.initial_gap_m+1e-15);
    }
    assert!(tip.pads.iter().any(|p|p.precompression_m!=tip.pads[0].precompression_m));
    assert!(s.compile_shell(&r,&shell,first(),1,1,n).is_err());
    assert!(r.surface_positions(&vec![0.;shell.mesh.nodes.len()],ShellFace::Positive).is_err());
}
#[test]
fn two_curved_felt_faces_drive_one_shell_and_keep_force_and_material_history_on_retry() {
    let make=|| {
        let selected=Selection{first:Some(spec()),second:Some(spec())};let mut e=build(&selected,true);
        let second=e.second_stick.unwrap();
        assert_eq!(e.acoustics.as_ref().unwrap().state_modes(),&(1..second.coordinate).collect::<Vec<_>>());
        assert!(inner(&e.system).felt_history(13).is_some());assert!(inner(&e.system).felt_history(14).is_none());
        // The six old stand histories remain first; the first and second felt
        // faces follow them. Both new faces begin separated, with zero force.
        for p in 6..14 {assert_eq!(inner(&e.system).felt_observation(p).unwrap().1,0.);}
        assert!((e.stick_weight-1./0.02_f64.sqrt()).abs()<1e-14);
        e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
            ImpactSubstepConfig{max_depth:8,max_attempts:511}).unwrap();
        e.system=e.system.with_stick_drives(vec![
            drive::Input{program:drive::Program::parse("0,0\n0.0002,0.05\n0.001,0").unwrap(),coordinate:0,tip_weight:e.stick_weight},
            drive::Input{program:drive::Program::parse("0,0\n0.0003,-0.02\n0.001,0").unwrap(),coordinate:second.coordinate,tip_weight:second.weight}
        ],2e-6,1024,e.force.len()).unwrap();e
    };
    let mut e=make();let mut clean=make();let initial=inner(&e.system).stored_energy_j();
    let gate=CancelGate::new_clock_free();let mut net=0.;let mut contact=0.0_f64;
    for tick in 0..768 {
        if tick==256 {
            let old=e.system.state().to_vec();let h=inner(&e.system).felt_history(6).unwrap();
            let cancel=CancelGate::new_clock_free();cancel.request();assert!(e.system.step(&e.force,&cancel).is_err());
            assert_eq!(e.system.state(),old);assert_eq!(inner(&e.system).felt_history(6).unwrap(),h);
        }
        let f=e.system.step(&e.force,&gate).unwrap();clean.system.step(&clean.force,&gate).unwrap();
        assert_eq!(e.system.state(),clean.system.state());net+=f.supplied_work_j-f.dissipated_energy_j;
        assert!((f.stored_energy_j-initial-net).abs()<1e-6);
        for p in 6..14 {contact=contact.max(inner(&e.system).felt_observation(p).unwrap().1);}
    }
    assert!(contact>0.);assert!((6..10).any(|p|inner(&e.system).felt_history(p).unwrap().eps_max>0.));
    assert!(e.system.state()[2..2*e.second_stick.unwrap().coordinate].iter().any(|x|*x!=0.));
}
#[test]
fn unselected_shell_preserves_the_original_hard_tip_trajectory() {
    let empty=Selection::default();let mut a=build(&empty,false);
    let mut b=splash_with_compliant_mute(1024,2e-6,false,first(),Some(specimen()),&[
        muffling::Muffler{surface:muffling::Surface::Shell,position_m:[0.075,0.],resistance_n_s_m:0.2}
    ],Some(second()),None).unwrap();
    let gate=CancelGate::new_clock_free();
    for _ in 0..64 {a.system.step(&a.force,&gate).unwrap();b.system.step(&b.force,&gate).unwrap();assert_eq!(a.system.state(),b.system.state());}
}

#[test]
fn a_curved_soft_face_does_not_leave_a_parallel_hertz_tip_at_the_old_plane() {
    let s=Spec::parse(&CARD.replace("0.003","0.012")).unwrap();
    let (shell,r)=shell_prepare::prepare(&specimen(),2e-6).unwrap();
    let tip=s.compile_shell(&r,&shell,first(),0,1,1+r.mode_count()).unwrap();
    let extra=tip.pads.iter().map(|p|-p.precompression_m).fold(f64::INFINITY,f64::min);
    assert!(extra>8.0*2e-6*first().speed_m_s,"fixture needs resolved curved clearance");
    let steps=((s.jaw.initial_gap_m+0.5*extra)/(2e-6*first().speed_m_s)).floor() as u64;
    assert!(steps>0 && steps<1024);
    let selection=Selection{first:Some(s),second:None};
    let mut e=splash_with_mallets(1024,2e-6,false,first(),Some(specimen()),&[],None,None,&selection).unwrap();
    e.system=e.system.into_analytic_nonlinear().unwrap();
    let gate=CancelGate::new_clock_free();
    for _ in 0..steps {
        e.system.step(&e.force,&gate).unwrap();
        for p in 6..10 {assert_eq!(inner(&e.system).felt_observation(p).unwrap().1,0.);}
    }
    // Mallet has crossed q=0, where an accidentally retained old Hertz tip
    // would already strike. It has NOT reached any actual felt contact site.
    assert!(e.system.state()[0]>0.);
    assert!(e.system.state()[2..2*e.force.len()].iter().all(|x|x.abs()<1e-8));
}
