use super::*;
use fs_couple::render::plate::impact::ImpactSubstepConfig;
use fs_phs::Storage;

fn spec()->Spec{Spec::parse(include_str!("estimated-hihat.fshh")).unwrap()}
fn shell(rigid:bool)->specimen::Specimen {
    let mut s=specimen::Specimen::reference();s.azimuths=8;
    if rigid{s.band_hz=[10.,11.];} // Deliberate rigid-only geometry/mass time fixture.
    s
}
fn bounds()->ImpactSubstepConfig{ImpactSubstepConfig{max_depth:8,max_attempts:511}}
fn gap(ob:&Obstacle,x:&[f64],n:usize)->f64{
    ob.collocation().chunks(n).zip(ob.gaps()).map(|(b,g)|g-b.iter().enumerate().map(|(i,w)|w*x[2*i]).sum::<f64>())
        .fold(f64::INFINITY,f64::min)
}
#[test]
fn paired_input_is_complete_and_invalid_physics_has_no_defaults() {
    let text=include_str!("estimated-hihat.fshh");let s=spec();assert_eq!(s.sites.len(),8);
    assert_eq!(s.contact,[1e7,1.5,0.1]);
    for bad in [text.replace("separation_m,0.0008","separation_m,-0.0008"),
        text.replace("carriage,0.05,200,0.5","carriage,0,200,0.5"),
        text.replace("damping_ratio,0.001,0.001","damping_ratio,0.001,-1"),
        text.replace("contact,10000000,1.5,0.1","contact,10000000,1.5,-1"),
        text.replace("site,0.098,0,0.125","site,0.098,0,0"),
        text.replace("pedal,0,0","pedal,0,1"),format!("{text}\nseparation_m,0.1"),
        text.replace("shell,lower,profile,estimated-splash.profile","")] {
        assert!(Spec::parse(&bad).is_err());
    }
    assert!(Spec::parse(&"x".repeat(16385)).is_err());
    assert!(is_command(Some("hihat-mic")));assert!(!is_command(Some("splash")));
}
#[test]
fn both_actual_shells_keep_mass_skin_and_force_moment_reciprocity() {
    let s=spec();let a=Shell::new(&shell(false),2e-6).unwrap();
    let mut material=shell(false);material.density_kg_m3*=1.3;
    let b=Shell::new(&material,2e-6).unwrap();
    let hi=1;let lo=hi+a.reduction.mode_count();let n=lo+b.reduction.mode_count()+1;
    let ob=collision(&s,&a,&b,lo,n).unwrap();assert_eq!(ob.n_points(),s.sites.len());
    assert_ne!(a.mesh.mass_kg,b.mesh.mass_kg);
    assert!(a.reduction.mode_count()>1&&b.reduction.mode_count()>1,"this fixture must retain elastic modes");
    for (i,row) in ob.collocation().chunks(n).enumerate(){
        assert_eq!(row[0],0.);assert_eq!(row[n-1],0.);
        assert!((row[hi]+1./a.mesh.mass_kg.sqrt()).abs()<1e-12);
        assert!((row[lo]+1./b.mesh.mass_kg.sqrt()).abs()<1e-12);
        let (x,y,_)=s.sites[i];let u=a.port([x,y],ShellFace::Negative).unwrap();
        let l=b.port([x,-y],ShellFace::Negative).unwrap();
        assert!((u.position_m[0]-l.position_m[0]).abs()<1e-12);
        assert!((u.position_m[1]+l.position_m[1]).abs()<1e-12);
        assert!((ob.gaps()[i]-s.separation-u.position_m[2]-l.position_m[2]).abs()<1e-12);
        // Arbitrary positive point force and physical opposing translations.
        let force=2.4;let vu=-0.3;let vl=0.2;
        let generalized=-force*(row[hi]*vu*a.mesh.mass_kg.sqrt()+row[lo]*vl*b.mesh.mass_kg.sqrt());
        assert!((generalized-force*(vu+vl)).abs()<1e-12);
    }
    // This collocation owner contains actual distributed potentials, not a mean row.
    struct Zero(usize);
    impl Storage for Zero {fn hamiltonian(&self,_:&[f64])->f64{0.}fn gradient(&self,_:&[f64],out:&mut[f64]){assert_eq!(out.len(),self.0);out.fill(0.);}}
    let contact=fs_dcontact::ContactStorage::new(Box::new(Zero(2*n)),n,vec![ob.clone()]).unwrap();
    let mut x=vec![0.;2*n];x[2*hi]=-(ob.gaps().iter().copied().fold(0.,f64::max)+0.0001)*a.mesh.mass_kg.sqrt();
    let mut g=vec![0.;2*n];contact.gradient(&x,&mut g);
    assert!(contact.hamiltonian(&x)>0.);
    assert!((g[2*hi]*a.mesh.mass_kg.sqrt()-g[2*lo]*b.mesh.mass_kg.sqrt()).abs()<1e-9);
    assert_eq!(g[0],0.);assert_eq!(g[2*(n-1)],0.);
}
#[test]
fn paired_scene_and_two_sticks_keep_source_addresses_and_no_phantom_radiation() {
    let s=spec();let a=shell(false);let b=shell(false);
    let stroke=Stroke{speed_m_s:0.8,position_m:Some([0.06,0.01])};
    let pair=build(&s,&a,&b,stroke,Some(Stroke{speed_m_s:0.,position_m:Some([-0.06,0.01])}),4,2e-6,true).unwrap();
    let e=pair.experiment;let sources=e.acoustics.as_ref().unwrap().state_modes();
    assert_eq!(sources,(pair.upper_modes.clone().chain(pair.lower_modes.clone())).collect::<Vec<_>>());
    assert!(!sources.contains(&0));assert!(!sources.contains(&pair.pedal.coordinate));
    assert!(!sources.contains(&e.second_stick.unwrap().coordinate));
    // Twelve independently retained washer/Kelvin histories, beyond true momenta.
    assert_eq!(e.system.state().len(),2*e.force.len()+12);
    for row in pair.collision.collocation().chunks(e.force.len()){
        assert_eq!(row[0],0.);assert_eq!(row[pair.pedal.coordinate],0.);
        assert_eq!(row[e.second_stick.unwrap().coordinate],0.);
    }
    let mut bad=spec();bad.separation=0.0001;
    assert!(build(&bad,&a,&b,stroke,None,1,2e-6,true).is_err());
    bad=spec();bad.sites[0]=(0.,0.,bad.sites[0].2);
    assert!(build(&bad,&a,&b,stroke,None,1,2e-6,false).is_err());
}
#[test]
fn pedal_closes_contacts_and_retracts_with_one_energy_and_force_clock() {
    let s=spec();let a=shell(true);let b=shell(true);let dt=1e-5;let steps=2400;
    let pair=build(&s,&a,&b,Stroke{speed_m_s:0.,position_m:None},None,steps,dt,false).unwrap();
    assert_eq!(pair.upper_modes.len(),1);assert_eq!(pair.lower_modes.len(),1);
    let mut e=pair.experiment;let n=e.force.len();let initial=e.system.state().to_vec();
    e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(bounds()).unwrap()
        .with_stick_drives(vec![mechanics::drive::Input{program:mechanics::drive::Program::parse(&s.pedal).unwrap(),
            coordinate:pair.pedal.coordinate,tip_weight:pair.pedal.weight}],dt,steps,n).unwrap();
    let gate=CancelGate::new_clock_free();gate.request();assert!(e.system.step(&e.force,&gate).is_err());
    assert_eq!(e.system.state(),initial);
    let gate=CancelGate::new_clock_free();let (mut work,mut loss,mut first,mut last)=(0.,0.,None,0.);
    let (mut penetration,mut lower_motion)=(0.0_f64,0.0_f64);
    for k in 0..steps {
        let f=e.system.step(&e.force,&gate).unwrap();if first.is_none(){first=Some(f.stored_energy_j+f.dissipated_energy_j-f.supplied_work_j);}
        assert!((f.time_s-(k+1) as f64*dt).abs()<1e-14);
        assert!(f.balance_residual_j.abs()<1e-7);work+=f.supplied_work_j;loss+=f.dissipated_energy_j;last=f.stored_energy_j;
        let x=e.system.state();penetration=penetration.max(-gap(&pair.collision,x,n));
        lower_motion=lower_motion.max(x[2*pair.lower_modes.start].abs());
    }
    assert!(penetration>0.,"the physical jaws must actually contact");
    assert!(lower_motion>1e-7,"the unforced lower shell must receive reciprocal collision forces");
    assert!(loss>0.);assert!((last+loss-first.unwrap()-work).abs()<1e-6);
    assert!(gap(&pair.collision,e.system.state(),n)>0.,"signed pedal retraction must reopen the gap");
    assert_ne!(e.system.state(),initial,"release is not a resonator reset");
}
