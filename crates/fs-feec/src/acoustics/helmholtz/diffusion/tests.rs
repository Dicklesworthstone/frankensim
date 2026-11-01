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

/// Above every mode of the 27-free-dof cube (`kuhn_cube(2)`, order 2): the whole free spectrum is kept.
const WIDE_CEILING: f64 = 1.0e5;

fn cube_modes(ceiling: f64, dense_limit: usize) -> DiffusionModes {
    let (complex, positions) = kuhn_cube(2);
    let fe = HelmholtzFe::new(&complex, &positions, 2,
        Medium { density_kg_m3:1.2, sound_speed_m_s:343.0 }).unwrap();
    fe.isothermal_diffusion_modes_with(&fe.boundary_faces(), ceiling, 128,
        &SliceOptions::default(), dense_limit, ||false).unwrap()
}

#[test]
fn dense_diffusion_modes_do_not_depend_on_the_eigenwindow_bitwise() {
    // The material caller's ceiling is max_decay_rate / diffusivity, so halving the gas pressure doubles it. The
    // unit-diffusivity spectrum is the geometry's: every ceiling must hand back the same bits for the modes it keeps.
    let wide = cube_modes(WIDE_CEILING, DENSE_FULL_PENCIL_DOFS);
    let total = wide.eigenvalues_per_m2.len();
    assert!(wide.dense && total > 12, "{} modes", total);
    // The cube's modes come in degenerate clusters whose members differ by an ulp or two, so a ceiling is only
    // unambiguous inside a real gap (a relative step above 1e-6), not between two members of one cluster.
    let gaps: Vec<usize> = (1..total)
        .filter(|&i| wide.eigenvalues_per_m2[i] > wide.eigenvalues_per_m2[i - 1] * (1.0 + 1e-6))
        .collect();
    assert!(gaps.len() >= 4, "{} gaps in {} modes", gaps.len(), total);
    for count in [gaps[0], gaps[1], gaps[gaps.len() / 2], gaps[gaps.len() - 1]] {
        let ceiling = 0.5 * (wide.eigenvalues_per_m2[count - 1] + wide.eigenvalues_per_m2[count]);
        let part = cube_modes(ceiling, DENSE_FULL_PENCIL_DOFS);
        assert_eq!(part.eigenvalues_per_m2, wide.eigenvalues_per_m2[..count]);
        assert_eq!(part.shapes, wide.shapes[..count]);
        assert_eq!(part.residual, wide.residual[..count]);
        assert_eq!((part.below_low, part.below_high), (0, count));
    }
    let (once, twice) = (cube_modes(40.0, DENSE_FULL_PENCIL_DOFS), cube_modes(80.0, DENSE_FULL_PENCIL_DOFS));
    assert!(twice.eigenvalues_per_m2.len() > once.eigenvalues_per_m2.len());
    assert_eq!(once.eigenvalues_per_m2, twice.eigenvalues_per_m2[..once.eigenvalues_per_m2.len()]);
}

#[test]
fn a_ceiling_keeps_exactly_the_modes_at_or_below_it() {
    // The negative of the test above: the filter is real. The cube's modes come in degenerate clusters, so the
    // expected counts are read off the window-independent spectrum itself: a ceiling one part in 1e12 under mode k
    // keeps the modes strictly below that, a ceiling equal to it keeps every mode up to and including it (the sparse
    // owner's window is (low, high]).
    let wide = cube_modes(WIDE_CEILING, DENSE_FULL_PENCIL_DOFS);
    for k in [2, 5, 9] {
        let lambda = wide.eigenvalues_per_m2[k];
        let under = wide.eigenvalues_per_m2.iter().filter(|&&l| l <= lambda * (1.0 - 1e-12)).count();
        let at = wide.eigenvalues_per_m2.iter().filter(|&&l| l <= lambda).count();
        assert!(under <= k && at > k, "mode {k}: {under} under, {at} at");
        assert_eq!(cube_modes(lambda * (1.0 - 1e-12), DENSE_FULL_PENCIL_DOFS).eigenvalues_per_m2.len(), under);
        assert_eq!(cube_modes(lambda, DENSE_FULL_PENCIL_DOFS).eigenvalues_per_m2.len(), at);
    }
}

#[test]
fn the_dense_and_sparse_owners_agree_on_the_spectrum() {
    // The two routes solve the same pencil; only the dense one promises window-independent bits.
    let dense = cube_modes(150.0, DENSE_FULL_PENCIL_DOFS);
    let sparse = cube_modes(150.0, 0);
    assert!(dense.dense && !sparse.dense);
    assert_eq!(dense.eigenvalues_per_m2.len(), sparse.eigenvalues_per_m2.len());
    for (a, b) in dense.eigenvalues_per_m2.iter().zip(&sparse.eigenvalues_per_m2) {
        assert!((a / b - 1.0).abs() < 1e-9, "{a} against {b}");
    }
}
