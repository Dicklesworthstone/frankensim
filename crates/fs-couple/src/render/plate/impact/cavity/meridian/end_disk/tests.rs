use super::*;
use crate::render::plate::impact::cavity::{CavityCoupling, meridian::{Budget, Station}};
use crate::vibroacoustic::AcousticMedium;
use fs_plate::shell::{head::TensionedDiskSpec, profile::ProfileBudget};

fn film() -> (TensionedDisk, Vec<ModePair>) {
    let disk = TensionedDisk::new(TensionedDiskSpec { radius_m: 1.0, thickness_m: 0.01,
        young_pa: 1e4, poisson: 0.3, density_kg_m3: 100.0, tension_n_m: 10.0,
        radial_intervals: 2, azimuths: 12 }, ProfileBudget {
        max_nodes: 100, max_triangles: 200, max_feature_evaluations: 1000 }).unwrap();
    let modes = fs_modal::slice_window(&disk.model.k, &disk.model.m, (0.01, 150.0),
        &fs_plate::SliceOptions::default()).unwrap().modes;
    assert!(!modes.is_empty());
    (disk,modes)
}
fn air(profile: Vec<Station>) -> MeridianCavity {
    MeridianCavity::new(profile, AcousticMedium {rho0:1.2,c0:1.0}, Budget {
        radial_intervals:4,maximum_azimuthal_order:1,maximum_frequency_hz:0.4,
        maximum_modes:8,maximum_nodes:100,maximum_matrix_entries:30_000,
        relative_tolerance:1e-7 }, &CancelGate::new_clock_free()).unwrap()
}
fn cylinder() -> MeridianCavity {
    air((0..=4).map(|i| Station {z_m:i as f64/4.0,radius_m:1.0}).collect())
}
fn port(n:usize) -> DiskPort { DiskPort {end:EndDisk::Upper,first_mode:2,structural_modes:n+3,positive_z:-1.0} }
fn budget() -> InterfaceBudget { InterfaceBudget {maximum_points:1000,maximum_terms:100_000} }

#[test]
fn constant_pressure_matches_actual_film_volume_work_and_preserves_other_addresses() {
    let (film,modes)=film();let air=cylinder();let p=port(modes.len());
    let interface=air.end_disk_interface(&film,&modes,p,budget(),&CancelGate::new_clock_free()).unwrap();
    let count=interface.cavity.omegas.len();assert!(count>1);
    for (i,mode) in modes.iter().enumerate() {
        let area=film.modal_area(&mode.phi).unwrap();
        assert!((interface.overlaps[(p.first_mode+i)*count]+area).abs()<1e-12);
    }
    for i in [0,1,p.structural_modes-1] {
        assert!(interface.overlaps[i*count..(i+1)*count].iter().all(|v|*v==0.0));
    }
    let coupling=CavityCoupling::new(&interface.cavity,p.structural_modes,&interface.overlaps,&vec![0.0;count]).unwrap();
    let mut state=vec![0.0;2*coupling.total_modes()];state[2*p.first_mode]=0.001;
    let mut pressure=vec![0.0;count];coupling.pressures_into(&state,&mut pressure).unwrap();
    let expected=-1.2*interface.overlaps[p.first_mode*count]*0.001/air.volume_m3();
    assert!((pressure[0]-expected).abs()<1e-14);
}

#[test]
fn reversing_the_physical_axis_reverses_all_pressure_work_not_the_cavity_basis() {
    let (film,modes)=film();let air=cylinder();let gate=CancelGate::new_clock_free();let p=port(modes.len());
    let a=air.end_disk_interface(&film,&modes,p,budget(),&gate).unwrap();
    let b=air.end_disk_interface(&film,&modes,DiskPort {positive_z:1.0,..p},budget(),&gate).unwrap();
    assert_eq!(a.cavity.omegas,b.cavity.omegas);assert_eq!(a.cavity.lambdas,b.cavity.lambdas);
    for (a,b) in a.overlaps.iter().zip(&b.overlaps) {assert_eq!(*a,-b);}
    let lower=air.end_disk_interface(&film,&modes,DiskPort {end:EndDisk::Lower,..p},budget(),&gate).unwrap();
    let count=a.cavity.omegas.len();
    for i in 0..p.structural_modes { assert!((a.overlaps[i*count]+lower.overlaps[i*count]).abs()<1e-12); }
}

#[test]
fn a_real_bowl_cap_uses_its_own_volume_and_does_not_create_a_lower_head() {
    let (film,modes)=film();let p=port(modes.len());let gate=CancelGate::new_clock_free();
    let bowl=air(vec![Station {z_m:0.0,radius_m:0.0},Station {z_m:0.5,radius_m:0.7},Station {z_m:1.0,radius_m:1.0}]);
    let interface=bowl.end_disk_interface(&film,&modes,p,budget(),&gate).unwrap();
    assert_eq!(interface.cavity.lambdas[0],bowl.volume_m3());
    assert!(bowl.volume_m3()<0.6*core::f64::consts::PI);
    assert!(bowl.end_disk_interface(&film,&modes,DiskPort {end:EndDisk::Lower,..p},budget(),&gate).is_err());
}

#[test]
fn mismatched_geometry_addresses_budgets_and_cancelled_assembly_refuse() {
    let (mut film,modes)=film();let air=cylinder();let p=port(modes.len());let gate=CancelGate::new_clock_free();
    for bad in [DiskPort {positive_z:0.0,..p},DiskPort {positive_z:f64::NAN,..p},
        DiskPort {first_mode:usize::MAX,..p},DiskPort {structural_modes:1,..p}] {
        assert!(air.end_disk_interface(&film,&modes,bad,budget(),&gate).is_err());
    }
    for b in [InterfaceBudget {maximum_points:1,..budget()},InterfaceBudget {maximum_terms:1,..budget()}] {
        assert!(air.end_disk_interface(&film,&modes,p,b,&gate).is_err());
    }
    film.spec.radius_m=0.9;
    assert!(air.end_disk_interface(&film,&modes,p,budget(),&gate).is_err());film.spec.radius_m=1.0;
    gate.request();
    assert!(matches!(air.end_disk_interface(&film,&modes,p,budget(),&gate),Err(ImpactError::Cancelled)));
}
