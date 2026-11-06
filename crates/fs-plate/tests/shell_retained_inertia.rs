//! Original reference inertia without an elastic section, for retained motion.

use fs_plate::PlateSection;
use fs_plate::shell::inertia::{
    prepare_reference_mass, reference_kinetic_bilinear, reference_kinetic_increment,
};
use fs_plate::shell::{
    ReferenceDensityMoments, ShellMesh, ShellSupport,
    assemble_shell_sections_with_reference_density_moments,
};

fn triangle(oblique: bool) -> ShellMesh {
    let nodes = if oblique {
        vec![[0.2, -0.4, 0.8], [1.4, 0.3, 1.1], [-0.1, 0.8, 1.6]]
    } else {
        vec![[0.0, 0.0, 0.0], [1.3, 0.0, 0.0], [0.2, 0.9, 0.0]]
    };
    ShellMesh::new(nodes, vec![[0, 1, 2]]).unwrap()
}

fn layered(bottom: f64, top: f64) -> ReferenceDensityMoments {
    // Two real density layers: 0 <= z/h <= .4 and .4 <= z/h <= 1.
    ReferenceDensityMoments::new(
        0.4 * bottom + 0.6 * top,
        0.08 * bottom + 0.42 * top,
        (0.064 * bottom + 0.936 * top) / 3.0,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4e-12 * expected.abs().max(1.0),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

fn product(matrix: &[f64; 324], a: &[f64; 18], b: &[f64; 18]) -> f64 {
    (0..18)
        .map(|i| {
            (0..18)
                .map(|j| a[i] * matrix[18 * i + j] * b[j])
                .sum::<f64>()
        })
        .sum()
}

fn field(seed: usize) -> [f64; 18] {
    core::array::from_fn(|i| ((i * 13 + seed * 7) % 23) as f64 * 0.031 - 0.32)
}

#[test]
fn split_mass_matches_original_pencil_in_flat_and_oblique_frames() {
    let uniform = [ReferenceDensityMoments::from_uniform(1200.0).unwrap(); 3];
    let layered = [
        layered(700.0, 2300.0),
        layered(900.0, 1700.0),
        layered(1100.0, 3100.0),
    ];
    for oblique in [false, true] {
        let mesh = triangle(oblique);
        let g = mesh.facet(0).unwrap();
        for moments in [uniform, layered] {
            let h = 0.07;
            let prepared = prepare_reference_mass(&g, [h; 3], moments).unwrap();
            // The existing elastic pencil is an independent public consumer;
            // no elastic section enters the new preparation API.
            let section = PlateSection::isotropic(3e9, 0.25, h, 1200.0).unwrap();
            let original = assemble_shell_sections_with_reference_density_moments(
                &mesh,
                &[section],
                &[moments],
                &[],
                ShellSupport::Free,
            )
            .unwrap();
            for i in 0..18 {
                for j in 0..18 {
                    assert_eq!(prepared.physical[18 * i + j], prepared.physical[18 * j + i]);
                    assert_eq!(prepared.drilling[18 * i + j], prepared.drilling[18 * j + i]);
                    let bi = 6 * (i / 6) + 3 * ((i % 6) / 3);
                    let bj = 6 * (j / 6) + 3 * ((j % 6) / 3);
                    let mut global = 0.0;
                    for p in 0..3 {
                        for q in 0..3 {
                            let index = 18 * (bi + p) + bj + q;
                            global += g.frame[p][i % 3]
                                * (prepared.physical[index] + prepared.drilling[index])
                                * g.frame[q][j % 3];
                        }
                    }
                    close(global, original.m.get(i, j));
                }
            }
            close(
                prepared.reference_mass_kg,
                g.area_m2 * h * moments.iter().map(|q| q.zeroth_kg_m3()).sum::<f64>() / 3.0,
            );
            let a = field(1);
            let b = field(2);
            let (physical, drilling) =
                reference_kinetic_bilinear(&g, [h; 3], moments, &a, &b).unwrap();
            close(physical, product(&prepared.physical, &a, &b));
            close(drilling, product(&prepared.drilling, &a, &b));
            let (physical, drilling) =
                reference_kinetic_increment(&g, [h; 3], moments, &a, &b).unwrap();
            let delta = core::array::from_fn(|i| b[i] - a[i]);
            let mean = core::array::from_fn(|i| 0.5 * a[i] + 0.5 * b[i]);
            for (account, matrix) in [
                (physical, &prepared.physical),
                (drilling, &prepared.drilling),
            ] {
                close(account.kinetic_before_j, 0.5 * product(matrix, &a, &a));
                close(account.kinetic_after_j, 0.5 * product(matrix, &b, &b));
                close(account.kinetic_change_j, product(matrix, &mean, &delta));
                close(
                    account.numerical_loss_j,
                    0.5 * product(matrix, &delta, &delta),
                );
                close(account.inertial_work_j, product(matrix, &b, &delta));
            }
        }
    }
}

#[test]
fn graded_tapered_layers_match_independent_volume_velocity_integral() {
    let mesh = triangle(false);
    let g = mesh.facet(0).unwrap();
    let h = [0.04, 0.08, 0.025];
    let bottom = [700.0, 900.0, 1100.0];
    let top = [2300.0, 1700.0, 3100.0];
    let moments = core::array::from_fn(|i| layered(bottom[i], top[i]));
    let prepared = prepare_reference_mass(&g, h, moments).unwrap();
    let a = field(3);
    let b = field(4);
    let spin = |v: &[f64; 18]| -> f64 {
        (0..3)
            .map(|i| 0.5 * (g.gradient[i][0] * v[6 * i + 1] - g.gradient[i][1] * v[6 * i]))
            .sum()
    };
    // Independent 5x5 surface and three-point physical-depth integration;
    // explicitly evaluate velocities inside both actual density layers.
    let points = [
        0.046910077030668,
        0.230765344947158,
        0.5,
        0.769234655052842,
        0.953089922969332,
    ];
    let weights = [
        0.118463442528095,
        0.239314335249683,
        0.284444444444444,
        0.239314335249683,
        0.118463442528095,
    ];
    let depth = [-(3.0_f64 / 5.0).sqrt(), 0.0, (3.0_f64 / 5.0).sqrt()];
    let depth_weights = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];
    let (mut physical, mut drilling, mut mass) = (0.0, 0.0, 0.0);
    for (i, &u) in points.iter().enumerate() {
        for (j, &v) in points.iter().enumerate() {
            let n = [(1.0 - u) * (1.0 - v), u, (1.0 - u) * v];
            let surface = 2.0 * g.area_m2 * (1.0 - u) * weights[i] * weights[j];
            let at = |values: [f64; 3]| -> f64 { (0..3).map(|k| n[k] * values[k]).sum() };
            let height = at(h);
            let rate = |values: &[f64; 18]| -> [f64; 6] {
                core::array::from_fn(|c| at([values[c], values[6 + c], values[12 + c]]))
            };
            let av = rate(&a);
            let bv = rate(&b);
            for (lo, hi, densities) in [(0.0, 0.4, bottom), (0.4, 1.0, top)] {
                for (&point, &weight) in depth.iter().zip(&depth_weights) {
                    let z = height * (0.5 * (lo + hi) + 0.5 * (hi - lo) * point);
                    let dm = surface * at(densities) * height * 0.5 * (hi - lo) * weight;
                    let va = [av[0] + z * av[4], av[1] - z * av[3], av[2]];
                    let vb = [bv[0] + z * bv[4], bv[1] - z * bv[3], bv[2]];
                    physical += dm * (0..3).map(|c| va[c] * vb[c]).sum::<f64>();
                    drilling += dm * z * z * (av[5] - spin(&a)) * (bv[5] - spin(&b));
                    mass += dm;
                }
            }
        }
    }
    close(product(&prepared.physical, &a, &b), physical);
    close(product(&prepared.drilling, &a, &b), drilling);
    close(prepared.reference_mass_kg, mass);
    let direct = reference_kinetic_bilinear(&g, h, moments, &a, &b).unwrap();
    close(direct.0, physical);
    close(direct.1, drilling);
}

#[test]
fn physical_kinetic_energy_preserves_signed_reference_surface_shift() {
    let mesh = triangle(true);
    let g = mesh.facet(0).unwrap();
    let h = 0.07;
    let shift = 0.09;
    let moments = [layered(800.0, 2400.0); 3];
    let s = shift / h;
    let shifted = moments.map(|q| {
        ReferenceDensityMoments::new(
            q.zeroth_kg_m3(),
            q.first_kg_m3() - s * q.zeroth_kg_m3(),
            q.second_kg_m3() - 2.0 * s * q.first_kg_m3() + s * s * q.zeroth_kg_m3(),
        )
        .unwrap()
    });
    assert!(shifted[0].first_kg_m3() < 0.0);
    let a = field(2);
    let b = field(5);
    let move_surface = |v: [f64; 18]| -> [f64; 18] {
        let mut moved = v;
        for i in 0..3 {
            moved[6 * i] += shift * v[6 * i + 4];
            moved[6 * i + 1] -= shift * v[6 * i + 3];
        }
        moved
    };
    let before = reference_kinetic_bilinear(&g, [h; 3], moments, &a, &b)
        .unwrap()
        .0;
    let after = reference_kinetic_bilinear(&g, [h; 3], shifted, &move_surface(a), &move_surface(b))
        .unwrap()
        .0;
    close(before, after);
}

#[test]
fn rigid_spin_has_zero_gauge_and_pure_drilling_has_zero_physical_energy() {
    for oblique in [false, true] {
        let mesh = triangle(oblique);
        let g = mesh.facet(0).unwrap();
        let moments = [layered(800.0, 2400.0); 3];
        let mut rigid = [0.0; 18];
        let mut drilling = [0.0; 18];
        for i in 0..3 {
            rigid[6 * i] = -0.3 * g.y[i];
            rigid[6 * i + 1] = 0.3 * g.x[i];
            rigid[6 * i + 5] = 0.3;
            drilling[6 * i + 5] = 0.3;
        }
        let (physical, gauge) =
            reference_kinetic_bilinear(&g, [0.07; 3], moments, &rigid, &rigid).unwrap();
        assert!(physical > 0.0);
        assert!(gauge >= 0.0 && gauge < 1e-26 * physical);
        let (physical, gauge) =
            reference_kinetic_bilinear(&g, [0.07; 3], moments, &drilling, &drilling).unwrap();
        assert_eq!(physical, 0.0);
        assert!(gauge > 0.0);
    }
}

#[test]
fn mass_only_ingress_refuses_invalid_geometry_inventory_and_products() {
    let mesh = triangle(false);
    let g = mesh.facet(0).unwrap();
    let moments = [ReferenceDensityMoments::from_uniform(1000.0).unwrap(); 3];
    for thickness in [[0.0; 3], [f64::NAN; 3], [f64::INFINITY; 3], [1e-200; 3]] {
        assert!(prepare_reference_mass(&g, thickness, moments).is_err());
    }
    let mut invalid = g;
    invalid.gradient[0][0] = f64::NAN;
    assert!(prepare_reference_mass(&invalid, [0.07; 3], moments).is_err());
    invalid = g;
    invalid.area_m2 = -1.0;
    assert!(prepare_reference_mass(&invalid, [0.07; 3], moments).is_err());
    assert!(
        reference_kinetic_bilinear(&g, [0.07; 3], moments, &[f64::MAX; 18], &[f64::MAX; 18])
            .is_err()
    );
    assert!(
        reference_kinetic_bilinear(&g, [0.07; 3], moments, &[f64::NAN; 18], &[0.0; 18]).is_err()
    );
}

#[test]
fn low_density_nearly_rigid_increment_retains_the_tiny_gauge_work_identity() {
    let mesh = triangle(true);
    let g = mesh.facet(0).unwrap();
    let moments = [ReferenceDensityMoments::from_uniform(1e-18).unwrap(); 3];
    let mut before = [0.0; 18];
    let mut after = [0.0; 18];
    let old_spin = 0.37;
    let new_spin = old_spin * (1.0 + 1e-12);
    for i in 0..3 {
        before[6 * i] = 0.11 - old_spin * g.y[i];
        before[6 * i + 1] = -0.07 + old_spin * g.x[i];
        before[6 * i + 2] = 0.02;
        before[6 * i + 5] = old_spin;
        after[6 * i] = 0.110000000003 - new_spin * g.y[i];
        after[6 * i + 1] = -0.069999999998 + new_spin * g.x[i];
        after[6 * i + 2] = 0.020000000001;
        after[6 * i + 5] = new_spin;
    }
    after[11] += 2e-16;
    let (physical, gauge) =
        reference_kinetic_increment(&g, [0.07; 3], moments, &before, &after).unwrap();
    assert!(physical.kinetic_after_j > 0.0);
    assert!(gauge.kinetic_after_j < 1e-25 * physical.kinetic_after_j);
    for account in [physical, gauge] {
        assert!(account.kinetic_before_j >= 0.0 && account.kinetic_after_j >= 0.0);
        assert!(account.numerical_loss_j >= 0.0);
        let scale = account
            .kinetic_before_j
            .max(account.kinetic_after_j)
            .max(account.numerical_loss_j)
            .max(account.inertial_work_j.abs());
        assert!(
            (account.inertial_work_j - account.kinetic_change_j - account.numerical_loss_j).abs()
                <= 4096.0 * f64::EPSILON * scale
        );
    }
}
