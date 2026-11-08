//! Original graded inventory through 3-D beam fields, chains and shell ties.
use fs_plate::shell::junction::{
    BeamOnShell, ShellComponent, TieBudget, TieSide, tie_beams_to_shells,
    tie_beams_to_shells_with_reference_density,
    tie_prestressed_beams_to_shells_with_reference_density, tie_shells,
};
use fs_plate::shell::stiffened::timoshenko::{
    BeamKind, ShearAreas, assemble_beam_chain_with_reference_density,
    assemble_prestressed_beam_chain_with_reference_density, segment_matrices_with_mass,
    segment_matrices_with_reference_density,
};
use fs_plate::shell::stiffened::{BeamMass, BeamSection};
use fs_plate::shell::{ShellMesh, ShellModel};
use fs_sparse::{Coo, Csr};

fn section() -> BeamSection {
    BeamSection {
        young_pa: 12e9,
        shear_pa: 0.8e9,
        density_kg_m3: 450.0,
        area_m2: 0.001,
        iy_m4: 1.3e-7,
        iz_m4: 4.5e-8,
        torsion_m4: 8e-8,
    }
}
fn kinds() -> [BeamKind; 2] {
    [
        BeamKind::EulerBernoulli,
        BeamKind::Timoshenko(ShearAreas {
            y_m2: 0.0008,
            z_m2: 0.0006,
        }),
    ]
}
fn near(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() <= 4e-12 * scale,
        "{a:e} != {b:e}; scale {scale:e}"
    );
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn bilinear(m: &[f64], a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .enumerate()
        .map(|(i, x)| {
            x * m[i * b.len()..(i + 1) * b.len()]
                .iter()
                .zip(b)
                .map(|(m, y)| m * y)
                .sum::<f64>()
        })
        .sum()
}
fn sparse_bilinear(m: &Csr, a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .enumerate()
        .map(|(i, x)| {
            let (cols, vals) = m.row(i);
            x * cols.iter().zip(vals).map(|(&j, v)| v * b[j]).sum::<f64>()
        })
        .sum()
}
// Independent integrals of rho(x), x*rho(x), x²*rho(x) on [0,L].
fn moments(l: f64, rho: [f64; 2]) -> [f64; 3] {
    [
        l * (rho[0] + rho[1]) / 2.0,
        l * l * (rho[0] + 2.0 * rho[1]) / 6.0,
        l.powi(3) * (rho[0] + 3.0 * rho[1]) / 12.0,
    ]
}
fn rigid(nodes: &[[f64; 3]], coordinate: usize) -> Vec<f64> {
    let mut axis = [0.0; 3];
    axis[coordinate % 3] = 1.0;
    nodes
        .iter()
        .flat_map(|&p| {
            let (v, w) = if coordinate < 3 {
                (axis, [0.0; 3])
            } else {
                (cross(axis, p), axis)
            };
            [v[0], v[1], v[2], w[0], w[1], w[2]]
        })
        .collect()
}

#[test]
fn graded_beams_retain_complete_rigid_spatial_moments_and_endpoint_reversal() {
    let s = section();
    let l = 1.7;
    let p0 = [0.4, -0.3, 0.2];
    let rho = [200.0, 800.0];
    let [r0, r1, r2] = moments(l, rho);
    let mass = s.area_m2 * r0;
    for axes in [
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [
            [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
            [2.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0],
            [2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0],
        ],
    ] {
        let [x, y, z] = axes;
        let p1 = core::array::from_fn(|i| p0[i] + l * x[i]);
        let first = core::array::from_fn(|i| mass * p0[i] + s.area_m2 * r1 * x[i]);
        let mut expected = [[0.0; 6]; 6];
        for a in 0..3 {
            expected[a][a] = mass;
            let mut axis = [0.0; 3];
            axis[a] = 1.0;
            let mixed = cross(axis, first);
            for b in 0..3 {
                expected[b][a + 3] = mixed[b];
                expected[a + 3][b] = mixed[b];
                let delta = if a == b { 1.0 } else { 0.0 };
                expected[a + 3][b + 3] = s.area_m2
                    * (r0 * (dot(p0, p0) * delta - p0[a] * p0[b])
                        + r1 * (2.0 * dot(p0, x) * delta - p0[a] * x[b] - x[a] * p0[b])
                        + r2 * (delta - x[a] * x[b]))
                    + r0 * ((s.iy_m4 + s.iz_m4) * x[a] * x[b]
                        + s.iy_m4 * y[a] * y[b]
                        + s.iz_m4 * z[a] * z[b]);
            }
        }
        for kind in kinds() {
            for reverse in [false, true] {
                let (nodes, density) = if reverse {
                    ([p1, p0], [rho[1], rho[0]])
                } else {
                    ([p0, p1], rho)
                };
                let (_, m) = segment_matrices_with_reference_density(
                    nodes[0],
                    nodes[1],
                    z,
                    s,
                    kind,
                    density,
                    BeamMass::Consistent,
                )
                .unwrap();
                for a in 0..6 {
                    for b in 0..6 {
                        near(
                            bilinear(&m, &rigid(&nodes, a), &rigid(&nodes, b)),
                            expected[a][b],
                            (expected[a][a] * expected[b][b]).sqrt(),
                        );
                    }
                }
            }
        }
    }
}

// Integrate physical polynomial velocity fields directly; no element shapes
// or quadrature coefficients enter this independent kinetic-energy oracle.
fn weighted_square(c: &[f64], l: f64, rho: [f64; 2]) -> f64 {
    c.iter()
        .enumerate()
        .map(|(i, a)| {
            c.iter()
                .enumerate()
                .map(|(j, b)| {
                    let n = i + j;
                    a * b
                        * l.powi((n + 1) as i32)
                        * (rho[0] / (n + 1) as f64 + (rho[1] - rho[0]) / (n + 2) as f64)
                })
                .sum::<f64>()
        })
        .sum()
}
#[test]
fn density_correlates_with_every_original_velocity_and_independent_section_rotation() {
    let l: f64 = 0.31;
    for kind in kinds() {
        for shear_pa in [8e6, 8e8, 8e10] {
            let s = BeamSection {
                shear_pa,
                ..section()
            };
            let mut energies = Vec::new();
            for rho in [[200.0, 800.0], [800.0, 200.0]] {
                let mut q = [0.0; 12];
                let mut exact = 0.0;
                let tip = |c: &[f64]| c.iter().rev().fold(0.0, |v, x| v * l + x);
                for (slots, sign, inertia, area, theta, v0) in [
                    ([1, 5, 7, 11], 1.0, s.iz_m4, 0.0008, [0.7, -0.8, 0.6], -0.2),
                    ([2, 4, 8, 10], -1.0, s.iy_m4, 0.0006, [-0.4, 0.5, -0.3], 0.1),
                ] {
                    // EI theta''=-V and v'=theta+V/(GA_s).
                    let gamma = match kind {
                        BeamKind::EulerBernoulli => 0.0,
                        BeamKind::Timoshenko(_) => {
                            -2.0 * s.young_pa * inertia * theta[2] / (s.shear_pa * area)
                        }
                    };
                    let v = [v0, theta[0] + gamma, 0.5 * theta[1], theta[2] / 3.0];
                    for (i, value) in
                        slots
                            .into_iter()
                            .zip([v0, sign * theta[0], tip(&v), sign * tip(&theta)])
                    {
                        q[i] = value;
                    }
                    exact += s.area_m2 * weighted_square(&v, l, rho)
                        + inertia * weighted_square(&theta, l, rho);
                }
                for (a, b, inertia, field) in [
                    (0, 6, s.area_m2, [0.4, -0.7]),
                    (3, 9, s.iy_m4 + s.iz_m4, [-0.6, 0.8]),
                ] {
                    q[a] = field[0];
                    q[b] = tip(&field);
                    exact += inertia * weighted_square(&field, l, rho);
                }
                let (_, m) = segment_matrices_with_reference_density(
                    [0.0; 3],
                    [l, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                    s,
                    kind,
                    rho,
                    BeamMass::Consistent,
                )
                .unwrap();
                let energy = bilinear(&m, &q, &q);
                near(energy, exact, exact);
                assert!(energy > 0.0);
                energies.push(energy);
            }
            assert!((energies[0] - energies[1]).abs() > 1e-4 * energies[0]);
        }
    }
}

#[test]
fn uniform_density_preserves_scalar_bits_and_lumping_preserves_graded_inventory() {
    let s = section();
    let l = 0.7;
    let end = [l, 0.0, 0.0];
    let z = [0.0, 0.0, 1.0];
    for kind in kinds() {
        for mass in [BeamMass::Consistent, BeamMass::Lumped] {
            for rho in [s.density_kg_m3, 700.0] {
                let scalar = BeamSection {
                    density_kg_m3: rho,
                    ..s
                };
                let original =
                    segment_matrices_with_mass([0.0; 3], end, z, scalar, kind, mass).unwrap();
                let graded = segment_matrices_with_reference_density(
                    [0.0; 3], end, z, s, kind, [rho; 2], mass,
                )
                .unwrap();
                assert_eq!(original.0.map(f64::to_bits), graded.0.map(f64::to_bits));
                assert_eq!(original.1.map(f64::to_bits), graded.1.map(f64::to_bits));
            }
        }
    }
    let rho = [200.0, 800.0];
    let (_, m) = segment_matrices_with_reference_density(
        [0.0; 3],
        end,
        z,
        s,
        kinds()[1],
        rho,
        BeamMass::Lumped,
    )
    .unwrap();
    for node in 0..2 {
        let partition = l * (2.0 * rho[node] + rho[1 - node]) / 6.0;
        for (dof, shape) in [
            s.area_m2,
            s.area_m2,
            s.area_m2,
            s.iy_m4 + s.iz_m4,
            s.iy_m4,
            s.iz_m4,
        ]
        .into_iter()
        .enumerate()
        {
            let i = 6 * node + dof;
            near(m[12 * i + i], partition * shape, partition * shape);
            for j in 0..12 {
                if i != j {
                    assert_eq!(m[12 * i + j], 0.0);
                }
            }
        }
    }
    let [r0, r1, _] = moments(l, rho);
    near(m[0] + m[78], s.area_m2 * r0, s.area_m2 * r0);
    near(l * m[78], s.area_m2 * r1, s.area_m2 * r1);
}

#[test]
fn native_chain_keeps_density_jumps_and_prestress_does_not_replace_inventory() {
    let nodes = [[0.0; 3], [0.4, 0.0, 0.0], [0.4, 0.7, 0.0]];
    let rho = [[200.0, 800.0], [1200.0, 300.0]];
    let s = section();
    let z = [0.0, 0.0, 1.0];
    let kind = kinds()[1];
    let (k, m) =
        assemble_beam_chain_with_reference_density(&nodes, z, s, kind, &rho, BeamMass::Consistent)
            .unwrap();
    let q: Vec<f64> = (0..18).map(|i| (i as f64 * 0.7).sin()).collect();
    let mut direct = 0.0;
    for e in 0..2 {
        let (_, element) = segment_matrices_with_reference_density(
            nodes[e],
            nodes[e + 1],
            z,
            s,
            kind,
            rho[e],
            BeamMass::Consistent,
        )
        .unwrap();
        direct += bilinear(&element, &q[6 * e..6 * e + 12], &q[6 * e..6 * e + 12]);
    }
    near(sparse_bilinear(&m, &q, &q), direct, direct);
    let translation = rigid(&nodes, 0);
    let total = s.area_m2 * (moments(0.4, rho[0])[0] + moments(0.7, rho[1])[0]);
    near(
        sparse_bilinear(&m, &translation, &translation),
        total,
        total,
    );
    let (loaded, loaded_mass) = assemble_prestressed_beam_chain_with_reference_density(
        &nodes,
        z,
        s,
        kind,
        &rho,
        &[20.0, -5.0],
    )
    .unwrap();
    assert_eq!(m, loaded_mass);
    assert_ne!(k, loaded);
    assert!(
        assemble_beam_chain_with_reference_density(
            &nodes,
            z,
            s,
            kind,
            &rho[..1],
            BeamMass::Consistent
        )
        .is_err()
    );
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            segment_matrices_with_reference_density(
                nodes[0],
                nodes[1],
                z,
                s,
                kind,
                [invalid, 300.0],
                BeamMass::Consistent
            )
            .is_err()
        );
    }
}

#[test]
fn native_off_node_shell_ties_transfer_graded_mass_centroid_and_rotary_work() {
    // A prepared positive shell pencil isolates the actual mortar/offset owner;
    // no surrogate shell constitutive response is part of this inertia check.
    let mesh = ShellMesh {
        nodes: vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
        tris: vec![[0, 1, 3], [0, 3, 2]],
    };
    let mut coo = Coo::new(24, 24);
    for i in 0..24 {
        coo.push(i, i, 0.01);
    }
    let model = ShellModel {
        k: coo.assemble(),
        m: coo.assemble(),
        dof_map: (0..24).map(Some).collect(),
        free: 24,
    };
    let thickness = [0.02; 4];
    let components = [ShellComponent {
        id: 0,
        mesh: &mesh,
        model: &model,
        nodal_thickness_m: &thickness,
    }];
    let budget = TieBudget {
        max_segments: 100,
        max_slave_dofs: 100,
        gap_tolerance_m: 1e-10,
        angle_tolerance_rad: 1e-8,
    };
    let tied = tie_shells(&components, &[], budget, || false).unwrap();
    let nodes = [[0.0, 0.0, 0.05], [0.4, 0.0, 0.05], [1.0, 0.0, 0.05]];
    let beams = [BeamOnShell {
        beam_nodes_m: &nodes,
        section_z: [0.0, 0.0, 1.0],
        section: section(),
        kind: kinds()[1],
        master: TieSide {
            component: 0,
            nodes: vec![0, 1],
        },
        side: 1,
        centroid_height_m: 0.04,
    }];
    let rho = [[200.0, 800.0], [1200.0, 300.0]];
    let graded = tie_beams_to_shells_with_reference_density(
        &tied,
        &components,
        &beams,
        &[&rho],
        budget,
        || false,
    )
    .unwrap();
    let increment = |a: &[f64], b: &[f64]| {
        sparse_bilinear(&graded.model.m, a, b) - sparse_bilinear(&tied.model.m, a, b)
    };
    let mut total = 0.0;
    let mut first = 0.0;
    let mut second = 0.0;
    for (e, &density) in rho.iter().enumerate() {
        let x = nodes[e][0];
        let [r0, r1, r2] = moments(nodes[e + 1][0] - x, density);
        total += r0;
        first += x * r0 + r1;
        second += x * x * r0 + 2.0 * x * r1 + r2;
    }
    let a = section().area_m2;
    let translation = rigid(&mesh.nodes, 2);
    let rotation = rigid(&mesh.nodes, 4);
    near(increment(&translation, &translation), a * total, a * total);
    near(increment(&translation, &rotation), -a * first, a * first);
    let angular = a * (second + 0.05_f64.powi(2) * total) + section().iy_m4 * total;
    near(increment(&rotation, &rotation), angular, angular);
    let lifted = graded.expand_beam(0, &rotation).unwrap();
    let expected = rigid(&nodes, 4);
    for (x, y) in lifted.iter().zip(&expected) {
        near(*x, *y, 1.0);
    }
    let loaded = tie_prestressed_beams_to_shells_with_reference_density(
        &tied,
        &components,
        &beams,
        &[&rho],
        &[&[20.0, -5.0]],
        budget,
        || false,
    )
    .unwrap();
    assert_eq!(loaded.model.m, graded.model.m);
    assert_ne!(loaded.model.k, graded.model.k);
    let uniform = [[section().density_kg_m3; 2]; 2];
    let old = tie_beams_to_shells(&tied, &components, &beams, budget, || false).unwrap();
    let same = tie_beams_to_shells_with_reference_density(
        &tied,
        &components,
        &beams,
        &[&uniform],
        budget,
        || false,
    )
    .unwrap();
    assert_eq!(old.model.m, same.model.m);
    assert_eq!(old.model.k, same.model.k);
    assert!(
        tie_beams_to_shells_with_reference_density(&tied, &components, &beams, &[], budget, || {
            false
        })
        .is_err()
    );
    assert!(
        tie_beams_to_shells_with_reference_density(
            &tied,
            &components,
            &beams,
            &[&rho[..1]],
            budget,
            || false
        )
        .is_err()
    );
    assert!(
        tie_beams_to_shells_with_reference_density(
            &tied,
            &components,
            &beams,
            &[&rho],
            budget,
            || true
        )
        .is_err()
    );
}
