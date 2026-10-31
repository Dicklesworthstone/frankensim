use super::tests::{Part, budget, options, seam};
use super::*;
use crate::shell::junction::TieOffset;
use crate::shell::point::triangle_point_rows;
use crate::shell::preload::tethers::ShellTether;
use crate::shell::preload::{equilibrate_shell_materials, tests::reference_newton};
use fs_sparse::Coo;

fn port(part: &Part, component: usize, node: usize, arm_m: [f64; 3]) -> AssemblyTetherPoint {
    let nodes = *part
        .mesh
        .tris
        .iter()
        .find(|tri| tri.contains(&node))
        .unwrap();
    AssemblyTetherPoint {
        component,
        nodes,
        weights: nodes.map(|n| if n == node { 1. } else { 0. }),
        arm_m,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1. + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}
fn put_force(load: &mut [f64], point: AssemblyTetherPoint, force: [f64; 3]) {
    for c in 0..3 {
        let mut direction = [0.; 3];
        direction[c] = 1.;
        let rows = triangle_point_rows(point.weights, direction, [point.arm_m; 3]).unwrap();
        for a in 0..3 {
            for d in 0..6 {
                load[6 * point.nodes[a] + d] += rows[a][d] * force[c];
            }
        }
    }
}
fn position(part: &Part, point: AssemblyTetherPoint, displacement: &[f64]) -> [f64; 3] {
    std::array::from_fn(|c| {
        let mut direction = [0.; 3];
        direction[c] = 1.;
        let rows = triangle_point_rows(point.weights, direction, [point.arm_m; 3]).unwrap();
        point.arm_m[c]
            + (0..3)
                .map(|a| {
                    point.weights[a] * part.mesh.nodes[point.nodes[a]][c]
                        + (0..6)
                            .map(|d| rows[a][d] * displacement[6 * point.nodes[a] + d])
                            .sum::<f64>()
                })
                .sum::<f64>()
    })
}

#[test]
fn world_anchors_reproduce_the_original_single_shell_tethers_bitwise() {
    let part = Part::rectangle(0, 2, 2, 0., true);
    let bearing = port(&part, 0, 5, [0.005, 0., 0.002]);
    let old: Vec<_> = [(0., 0.), (0.4, 1000.)]
        .into_iter()
        .map(|(x, stiffness)| ShellTether {
            nodes: bearing.nodes,
            weights: bearing.weights,
            arm_m: bearing.arm_m,
            anchor_m: [x, 0.1, -0.04],
            reference_tension_n: 5.,
            axial_stiffness_n_per_m: stiffness,
        })
        .collect();
    let connections: Vec<_> = old
        .iter()
        .map(|t| AssemblyTether {
            bearing,
            anchor: AssemblyTetherAnchor::World(t.anchor_m),
            reference_tension_n: t.reference_tension_n,
            axial_stiffness_n_per_m: t.axial_stiffness_n_per_m,
        })
        .collect();
    let direct = equilibrate_shell_materials(
        &part.mesh,
        &part.sections,
        &part.fixed,
        ShellSupport::Clamped,
        &[],
        &part.loads,
        &old,
        options(),
        &mut reference_newton,
    )
    .unwrap();
    let joint = equilibrate_tethered_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        &connections,
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    assert_eq!(joint.equilibrium.displacements[0], direct.displacement);
    assert_eq!(joint.equilibrium.stored_energy_j, direct.stored_energy_j);
    assert_eq!(joint.equilibrium.residual_force_n, direct.residual_force_n);
    for row in 0..direct.model.free {
        assert_eq!(
            joint.equilibrium.assembly.model.k.row(row),
            direct.model.k.row(row)
        );
        assert_eq!(
            joint.equilibrium.assembly.model.m.row(row),
            direct.model.m.row(row)
        );
    }
    for (a, b) in joint.tether_responses.iter().zip(&direct.tether_responses) {
        assert_eq!(a.force_n, b.force_n);
        assert_eq!(a.tension_n, b.tension_n);
        assert_eq!(a.length_m, b.length_m);
        assert_eq!(a.length_change_m, b.length_change_m);
        assert_eq!(a.potential_change_j, b.potential_change_j);
    }
}

#[test]
fn two_moving_bodies_share_one_tension_energy_and_cross_tangent() {
    let mut left = Part::rectangle(0, 2, 2, 0., true);
    let mut right = Part::rectangle(3, 2, 2, 0.08, true);
    let a = port(&left, 0, 5, [0.01, -0.005, 0.004]);
    let b = port(&right, 1, 5, [-0.01, 0.005, -0.002]);
    let tether = AssemblyTether {
        bearing: a,
        anchor: AssemblyTetherAnchor::Shell(b),
        reference_tension_n: 10.,
        axial_stiffness_n_per_m: 1500.,
    };
    let joint = equilibrate_tethered_shell_assembly(
        &[left.input(), right.input()],
        &[],
        &[],
        budget(),
        &[tether],
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    let response = joint.tether_responses[0];
    let forces = response.endpoint_forces_n();
    assert!(response.tension_n > 0. && response.tension_n < 10.);
    for c in 0..3 {
        assert_eq!(forces[0][c] + forces[1][c], 0.);
    }
    let pa = position(&left, a, &joint.equilibrium.displacements[0]);
    let pb = position(&right, b, &joint.equilibrium.displacements[1]);
    let chord: [f64; 3] = std::array::from_fn(|c| pb[c] - pa[c]);
    let length = chord.iter().fold(0.0_f64, |a, b| a.hypot(*b));
    near(length, response.length_m, 1e-14);
    for c in 0..3 {
        near(forces[0][c], response.tension_n * chord[c] / length, 1e-12);
    }
    near(
        response.potential_change_j,
        response.length_change_m * (10. + 0.5 * 1500. * response.length_change_m),
        1e-13,
    );
    // Freeze the CONVERGED forces through the original reciprocal point owner.
    // They must recover the same structural shape, but cannot reproduce the
    // material/geometric coupling of two moving endpoints in the tangent.
    put_force(&mut left.loads, a, forces[0]);
    put_force(&mut right.loads, b, forces[1]);
    let frozen = equilibrate_shell_assembly(
        &[left.input(), right.input()],
        &[],
        &[],
        budget(),
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    for (a, b) in joint
        .equilibrium
        .displacements
        .iter()
        .flatten()
        .zip(frozen.displacements.iter().flatten())
    {
        near(*a, *b, 2e-10);
    }
    let iz = joint.equilibrium.assembly.map[0][5][2].unwrap();
    let jz = joint.equilibrium.assembly.map[1][5][2].unwrap();
    let geometric = response.tension_n / response.length_m;
    let expected = -(geometric + (1500. - geometric) * (chord[2] / length).powi(2));
    assert_eq!(frozen.assembly.model.k.get(iz, jz), 0.);
    near(
        joint.equilibrium.assembly.model.k.get(iz, jz),
        expected,
        1e-11,
    );
    assert!(joint.equilibrium.assembly.model.k.get(iz, iz) > frozen.assembly.model.k.get(iz, iz));
    for row in 0..frozen.assembly.model.free {
        assert_eq!(
            joint.equilibrium.assembly.model.m.row(row),
            frozen.assembly.model.m.row(row)
        );
    }
}

#[test]
fn shared_component_and_slave_rows_cancel_before_force_and_tangent_assembly() {
    let left = Part::rectangle(0, 2, 2, 0., true);
    let right = Part::rectangle(2, 2, 2, 0., false);
    for tied in [false, true] {
        let components = if tied {
            vec![left.input(), right.input()]
        } else {
            vec![left.input()]
        };
        let ties = if tied {
            vec![seam(2, TieOffset::Midsurface)]
        } else {
            vec![]
        };
        let a = if tied {
            port(&right, 1, 3, [0., 0., 0.01])
        } else {
            port(&left, 0, 5, [0., 0., 0.01])
        };
        let b = port(&left, 0, 5, [0., 0., -0.01]);
        let tether = AssemblyTether {
            bearing: a,
            anchor: AssemblyTetherAnchor::Shell(b),
            reference_tension_n: 8.,
            axial_stiffness_n_per_m: 500.,
        };
        let base = equilibrate_shell_assembly(
            &components,
            &ties,
            &[],
            budget(),
            options(),
            &mut reference_newton,
            || false,
        )
        .unwrap();
        let joint = equilibrate_tethered_shell_assembly(
            &components,
            &ties,
            &[],
            budget(),
            &[tether],
            options(),
            &mut reference_newton,
            || false,
        )
        .unwrap();
        assert!(
            joint
                .equilibrium
                .displacements
                .iter()
                .flatten()
                .all(|v| *v == 0.)
        );
        assert_eq!(joint.tether_responses[0].potential_change_j, 0.);
        let model = &joint.equilibrium.assembly.model;
        for (full, dof) in model.dof_map.iter().enumerate() {
            if full % 6 < 3 {
                if let Some(i) = dof {
                    for j in 0..model.free {
                        near(model.k.get(*i, j), base.assembly.model.k.get(*i, j), 1e-12);
                    }
                }
            }
        }
        // Only the original LINEAR offset-rotation map remains. Its transverse
        // radial curvature is T*delta_arm, without spurious translational K.
        let rx = joint.equilibrium.assembly.map[0][5][3].unwrap();
        near(
            model.k.get(rx, rx) - base.assembly.model.k.get(rx, rx),
            8. * 0.02,
            1e-10,
        );
        let reverse = AssemblyTether {
            bearing: b,
            anchor: AssemblyTetherAnchor::Shell(a),
            ..tether
        };
        let reversed = equilibrate_tethered_shell_assembly(
            &components,
            &ties,
            &[],
            budget(),
            &[reverse],
            options(),
            &mut reference_newton,
            || false,
        )
        .unwrap();
        for i in 0..model.free {
            for j in 0..model.free {
                near(
                    model.k.get(i, j),
                    reversed.equilibrium.assembly.model.k.get(i, j),
                    1e-12,
                );
            }
        }
        for c in 0..3 {
            near(
                joint.tether_responses[0].force_n[c],
                -reversed.tether_responses[0].force_n[c],
                1e-13,
            );
        }
    }
}

#[test]
fn nonmatching_tether_ports_keep_relative_energy_force_and_consistent_tangent() {
    let left = Part::rectangle(0, 2, 2, 0., true);
    let right = Part::rectangle(2, 2, 3, 0., false);
    let a = port(&right, 1, 3, [0.002, -0.001, 0.004]);
    let b = port(&left, 0, 4, [-0.001, 0.002, -0.002]);
    let tether = AssemblyTether {
        bearing: a,
        anchor: AssemblyTetherAnchor::Shell(b),
        reference_tension_n: 0.3,
        axial_stiffness_n_per_m: 20.,
    };
    let mut checked = false;
    equilibrate_tethered_shell_assembly(
        &[left.input(), right.input()],
        &[seam(3, TieOffset::Midsurface)],
        &[],
        budget(),
        &[tether],
        options(),
        &mut |p, initial, limit| {
            if !checked {
                checked = true;
                let connection = &p.tethers[0];
                let x: Vec<_> = (0..p.dimension())
                    .map(|i| ((i * 17 % 19) as f64 - 9.) * 1e-6)
                    .collect();
                let u = p.physical(&x);
                let v = p.physical(
                    &(0..p.dimension())
                        .map(|i| ((i * 11 % 13) as f64 - 6.) / 6.)
                        .collect::<Vec<_>>(),
                );
                let h = 1e-8;
                let plus: Vec<_> = u.iter().zip(&v).map(|(u, v)| u + h * v).collect();
                let minus: Vec<_> = u.iter().zip(&v).map(|(u, v)| u - h * v).collect();
                let mut gradient = vec![0.; p.dimension()];
                connection.add_gradient(&u, &mut gradient, 1.).unwrap();
                let fd = (connection.response(&plus).unwrap().potential_change_j
                    - connection.response(&minus).unwrap().potential_change_j)
                    / (2. * h);
                near(fd, gradient.iter().zip(&v).map(|(a, b)| a * b).sum(), 1e-7);
                let mut matrix = Coo::new(p.dimension(), p.dimension());
                connection.add_tangent(&u, &mut matrix, 1.).unwrap();
                let matrix = matrix.assemble();
                let mut product = vec![0.; p.dimension()];
                matrix.spmv(&v, &mut product);
                let mut action = vec![0.; p.dimension()];
                connection.add_action(&u, &v, &mut action, 1.).unwrap();
                let mut a = vec![0.; p.dimension()];
                let mut b = a.clone();
                connection.add_gradient(&plus, &mut a, 1.).unwrap();
                connection.add_gradient(&minus, &mut b, 1.).unwrap();
                for i in 0..p.dimension() {
                    near(action[i], product[i], 1e-11);
                    near((a[i] - b[i]) / (2. * h), action[i], 1e-7);
                    let (cols, values) = matrix.row(i);
                    for (&j, &value) in cols.iter().zip(values) {
                        assert_eq!(value, matrix.get(j, i));
                    }
                }
            }
            reference_newton(p, initial, limit)
        },
        || false,
    )
    .unwrap();
    assert!(checked);
}

#[test]
fn slack_law_physical_support_requirements_caps_and_final_cancellation_are_retained() {
    let mut part = Part::rectangle(0, 2, 2, 0., true);
    part.loads[6 * 5] = 1.; // Move the bearing toward its anchor: the wire is slack.
    let a = port(&part, 0, 5, [0.; 3]);
    let tether = AssemblyTether {
        bearing: a,
        anchor: AssemblyTetherAnchor::World([0.4, 0.1, 0.]),
        reference_tension_n: 0.,
        axial_stiffness_n_per_m: 1000.,
    };
    let bare = equilibrate_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    let slack = equilibrate_tethered_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        &[tether],
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    assert!(slack.tether_responses[0].length_change_m < 0.);
    assert_eq!(slack.tether_responses[0].force_n, [0.; 3]);
    assert_eq!(slack.tether_responses[0].potential_change_j, 0.);
    assert_eq!(slack.equilibrium.displacements, bare.displacements);
    for i in 0..bare.assembly.model.free {
        for j in 0..bare.assembly.model.free {
            assert_eq!(
                slack.equilibrium.assembly.model.k.get(i, j),
                bare.assembly.model.k.get(i, j)
            );
        }
    }
    let solve = |input: ShellPreloadComponent<'_>, tethers: &[AssemblyTether]| {
        equilibrate_tethered_shell_assembly(
            &[input],
            &[],
            &[],
            budget(),
            tethers,
            options(),
            &mut reference_newton,
            || false,
        )
    };
    let mut free = part.input();
    free.boundary = &[];
    free.support = ShellSupport::Free;
    assert!(solve(free, &[tether]).unwrap_err().contains("rigid modes"));
    let wrong = AssemblyTether {
        anchor: AssemblyTetherAnchor::Shell(AssemblyTetherPoint { component: 1, ..a }),
        ..tether
    };
    assert!(
        solve(part.input(), &[wrong])
            .unwrap_err()
            .contains("unknown tether component")
    );
    let collapsed = AssemblyTether {
        anchor: AssemblyTetherAnchor::Shell(a),
        ..tether
    };
    assert!(
        solve(part.input(), &[collapsed])
            .unwrap_err()
            .contains("coincides")
    );
    assert!(
        solve(part.input(), &vec![tether; 4097])
            .unwrap_err()
            .contains("4096")
    );
    let mut polls = 0;
    equilibrate_tethered_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        &[tether],
        options(),
        &mut reference_newton,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut repeated = 0;
    let error = equilibrate_tethered_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        &[tether],
        options(),
        &mut reference_newton,
        || {
            repeated += 1;
            repeated == polls
        },
    )
    .unwrap_err();
    assert_eq!(repeated, polls);
    assert!(error.contains("cancelled"));
}
