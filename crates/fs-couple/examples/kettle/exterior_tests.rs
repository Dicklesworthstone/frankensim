use super::*;
use fs_bem::panel3d::SpherePanels;
use fs_plate::shell::{head::TensionedDiskSpec,profile::ProfileBudget};

fn fixture()->(TensionedDisk,Vec<ModePair>,Vec<Station>) {
    let film=TensionedDisk::new(TensionedDiskSpec {radius_m:1.0,thickness_m:0.01,young_pa:1e4,
        poisson:0.3,density_kg_m3:100.0,tension_n_m:10.0,radial_intervals:2,azimuths:8},
        ProfileBudget {max_nodes:100,max_triangles:200,max_feature_evaluations:1000}).unwrap();
    let modes=fs_modal::slice_window(&film.model.k,&film.model.m,(0.01,150.0),
        &fs_plate::SliceOptions::default()).unwrap().modes;
    assert!(!modes.is_empty());
    (film,modes,vec![Station {z_m:0.0,radius_m:0.0},Station {z_m:0.5,radius_m:0.7},Station {z_m:1.0,radius_m:1.0}])
}
fn boundary(f:&TensionedDisk,m:&[ModePair],p:&[Station])->Boundary {
    Boundary::new(f,m,p,1024,32768,&CancelGate::new_clock_free()).unwrap()
}
fn volume(triangles:&[[[f64;3];3]])->f64 {
    triangles.iter().map(|&[a,b,c]| (a[0]*(b[1]*c[2]-b[2]*c[1])
        +a[1]*(b[2]*c[0]-b[0]*c[2])+a[2]*(b[0]*c[1]-b[1]*c[0]))/6.0).sum()
}

#[test]
fn oriented_closed_bowl_preserves_volume_and_actual_film_volume_velocity() {
    let (film,modes,profile)=fixture();let b=boundary(&film,&modes,&profile);
    let panels=SpherePanels::from_triangles(b.triangles.clone()).unwrap();
    // Polygonal frusta, not the circular continuum volume at finite azimuths.
    let factor=4.0*(core::f64::consts::TAU/8.0).sin();
    let expected=profile.windows(2).map(|s|(s[1].z_m-s[0].z_m)*factor/3.0*
        (s[0].radius_m.powi(2)+s[0].radius_m*s[1].radius_m+s[1].radius_m.powi(2))).sum::<f64>();
    assert!((volume(&b.triangles)-expected).abs()<1e-13);
    assert_eq!(b.head_panels,film.mesh.tris.len());
    for (j,m) in modes.iter().enumerate() {
        assert!(b.weights[j][b.head_panels..].iter().all(|v|*v==0.0));
        let flow=b.weights[j].iter().zip(panels.areas()).map(|(v,a)|v*a).sum::<f64>();
        assert!((flow+film.modal_area(&m.phi).unwrap()).abs()<1e-12);
    }
    assert!(panels.normals()[..b.head_panels].iter().all(|n|n[2]>0.999999999999));
    for p in b.triangles.iter().flatten() {
        assert!(p[0].hypot(p[1]).hypot(p[2])<=b.radius_m);
    }
}

#[test]
fn axial_translation_changes_only_the_declared_origin_not_the_acoustic_shape() {
    let (film,modes,profile)=fixture();let a=boundary(&film,&modes,&profile);
    let shifted:Vec<_>=profile.iter().map(|s|Station {z_m:s.z_m+8.0,..*s}).collect();
    let b=boundary(&film,&modes,&shifted);
    assert_eq!(a.weights,b.weights);assert_eq!(a.triangles,b.triangles);
    assert_eq!(b.origin_m[2]-a.origin_m[2],8.0);
}

#[test]
fn invalid_geometry_topology_and_work_refuse_before_bem() {
    let (mut film,modes,profile)=fixture();let gate=CancelGate::new_clock_free();
    assert!(Boundary::new(&film,&modes,&profile,1,32768,&gate).is_err());
    assert!(Boundary::new(&film,&modes,&profile,1024,1,&gate).is_err());
    let mut wrong=profile.clone();wrong[0].radius_m=0.1;
    assert!(Boundary::new(&film,&modes,&wrong,1024,32768,&gate).is_err());
    wrong=profile.clone();wrong[2].radius_m=0.99;
    assert!(Boundary::new(&film,&modes,&wrong,1024,32768,&gate).is_err());
    film.mesh.tris[0].swap(1,2);
    assert!(Boundary::new(&film,&modes,&profile,1024,32768,&gate).is_err());
    film.mesh.tris[0].swap(1,2);gate.request();
    assert!(Boundary::new(&film,&modes,&profile,1024,32768,&gate).is_err());
}

#[test]
fn actual_closed_surface_drives_reciprocal_bem_with_linear_pressure_and_quadratic_power() {
    use fs_bem::{helmholtz::Medium,reciprocal_radiation::{Prepared,Budget}};
    use fs_math::c64::C64;
    let (film,modes,profile)=fixture();let b=boundary(&film,&modes,&profile);
    let panels=SpherePanels::from_triangles(b.triangles).unwrap();
    let prepared=Prepared::new(&panels,0.5,Budget {maximum_unknowns:256,
        maximum_cache_bytes:16*1024*1024,maximum_work:1_000_000_000,angular_refinement:1},||false).unwrap();
    let one:Vec<_>=b.weights[0].iter().map(|&v|C64::from_re(v)).collect();
    let two:Vec<_>=one.iter().map(|v|v.scale(2.0)).collect();let zero=vec![C64::ZERO;one.len()];
    let solutions=prepared.solve_batch(0.5,Medium {density:1.2,sound_speed:343.0},
        &[&one,&two,&zero],||false).unwrap();
    let a=solutions[0].far_field([0.3,0.4,0.8]).unwrap();
    let c=solutions[1].far_field([0.3,0.4,0.8]).unwrap();
    assert!(a.abs()>1e-8);
    assert!((c-a.scale(2.0)).abs()<1e-10*a.abs());
    assert!(solutions[0].outgoing_power_w>0.0);
    assert!((solutions[1].outgoing_power_w-4.0*solutions[0].outgoing_power_w).abs()
        <1e-9*solutions[0].outgoing_power_w);
    assert_eq!(solutions[2].far_field([0.3,0.4,0.8]).unwrap(),C64::ZERO);
}
