use super::*;
use crate::fixtures::kuhn_cube;

#[test]
fn dirichlet_and_mixed_wall_heat_modes_use_the_original_volume_pencil() {
    let (complex, unit) = kuhn_cube(3);
    let positions: Vec<_> = unit.iter().map(|p| [p[0],0.8*p[1],0.6*p[2]]).collect();
    let fe = HelmholtzFe::new(&complex, &positions, 2,
        Medium { density_kg_m3:1.2, sound_speed_m_s:343.0 }).unwrap();
    let walls = fe.boundary_faces();
    let full = fe.isothermal_diffusion_modes(&walls,100.0,16,
        &SliceOptions::default(),||false).unwrap();
    let exact = core::f64::consts::PI.powi(2)*(1.0+1.0/0.8_f64.powi(2)+1.0/0.6_f64.powi(2));
    assert!((full.eigenvalues_per_m2[0]/exact-1.0).abs() < 0.04);
    let ends: Vec<_> = walls.iter().copied().filter(|face|
        face.iter().all(|&v|positions[v as usize][0]==0.0)
        || face.iter().all(|&v|positions[v as usize][0]==1.0)).collect();
    let mixed = fe.isothermal_diffusion_modes(&ends,15.0,4,
        &SliceOptions::default(),||false).unwrap();
    assert!((mixed.eigenvalues_per_m2[0]/core::f64::consts::PI.powi(2)-1.0).abs() < 0.01);
    assert!(mixed.eigenvalues_per_m2[0] < full.eigenvalues_per_m2[0]);
    for (i,a) in full.shapes.iter().enumerate() {
        for (j,b) in full.shapes.iter().enumerate() {
            let value = dot(a,&spmv_new(fe.mass(),b));
            assert!((value-if i==j {1.0} else {0.0}).abs() < 1e-8);
        }
        for face in &walls { for &vertex in face { assert_eq!(a[vertex as usize],0.0); } }
    }
}

#[test]
fn incomplete_or_duplicate_thermal_domains_refuse() {
    let (complex,positions) = kuhn_cube(2);
    let fe = HelmholtzFe::new(&complex,&positions,2,
        Medium { density_kg_m3:1.2, sound_speed_m_s:343.0 }).unwrap();
    let walls = fe.boundary_faces();
    for (cutoff,count) in [(0.0,8),(1.0,8),(10000.0,1)] {
        assert!(fe.isothermal_diffusion_modes(&walls,cutoff,count,
            &SliceOptions::default(),||false).is_err());
    }
    assert!(fe.isothermal_diffusion_modes(&[walls[0],walls[0]],100.0,8,
        &SliceOptions::default(),||false).is_err());
    assert!(fe.isothermal_diffusion_modes(&walls,100.0,8,
        &SliceOptions::default(),||true).is_err());
}
