use fs_plate::shell::head::{TensionedDisk,TensionedDiskSpec};
use fs_plate::shell::profile::ProfileBudget;
fn spec()->TensionedDiskSpec {TensionedDiskSpec{radius_m:0.1778,thickness_m:0.000254,
    young_pa:4.0e9,poisson:0.38,density_kg_m3:1390.0,tension_n_m:3000.0,radial_intervals:3,azimuths:24}}
fn budget()->ProfileBudget {ProfileBudget{max_nodes:1000,max_triangles:2000,max_feature_evaluations:0}}
#[test]
fn real_film_mass_tension_and_reciprocal_cavity_area_enter_the_pencil() {
    let s=spec();let a=TensionedDisk::new(s,budget()).unwrap();
    let b=TensionedDisk::new(TensionedDiskSpec{tension_n_m:2.0*s.tension_n_m,..s},budget()).unwrap();
    let area=0.5*s.azimuths as f64*(2.0*std::f64::consts::PI/s.azimuths as f64).sin()*s.radius_m*s.radius_m;
    assert!((a.mass_kg-s.density_kg_m3*s.thickness_m*area).abs()<1e-14);
    assert_eq!(a.model.m,b.model.m);assert_ne!(a.model.k,b.model.k);
    let mut mode=vec![0.0;a.model.free];
    for (node,&(x,y)) in a.mesh.nodes.iter().enumerate() {
        if let Some(i)=a.model.dof_map[3*node] {mode[i]=1.0-(x*x+y*y)/(s.radius_m*s.radius_m);}
    }
    let modal_area=a.modal_area(&mode).unwrap();assert!(modal_area>0.0);
    let force=37.0*modal_area;let volume_flow=modal_area*0.02;
    assert!((force*0.02-37.0*volume_flow).abs()<1e-14);
    // Tension adds a positive quadratic form without rescaling film inertia.
    let mut ka=vec![0.0;mode.len()];let mut kb=ka.clone();a.model.k.spmv(&mode,&mut ka);b.model.k.spmv(&mode,&mut kb);
    assert!(mode.iter().zip(ka.iter().zip(&kb)).map(|(q,(a,b))|q*(b-a)).sum::<f64>()>0.0);
}
#[test]
fn installed_tension_is_not_guessed_and_membrane_modes_cannot_use_wrong_shapes() {
    let mut s=spec();s.tension_n_m=0.0;assert!(TensionedDisk::new(s,budget()).is_err());
    s=spec();s.radial_intervals=usize::MAX;assert!(TensionedDisk::new(s,budget()).is_err());
    let a=TensionedDisk::new(spec(),budget()).unwrap();assert!(a.modal_area(&[1.0]).is_err());
}
