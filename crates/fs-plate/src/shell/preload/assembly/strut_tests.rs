use super::tests::{Part, budget, options, seam};
use super::*;
use crate::shell::junction::TieOffset;
use crate::shell::preload::tests::reference_newton;
use fs_sparse::Coo;

fn port(part: &Part, component: usize, node: usize) -> AssemblyTetherPoint {
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
        arm_m: [0.; 3],
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1. + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}

#[test]
fn compressed_pair_has_analytic_force_and_negative_transverse_geometric_stiffness() {
    let mut lower = Part::rectangle(0, 2, 2, 0., true);
    let mut upper = Part::rectangle(0, 2, 2, 0.1, true);
    let a = port(&lower, 0, 5);
    let b = port(&upper, 1, 5);
    let strut = AssemblyStrut {
        bearing: a,
        anchor: AssemblyTetherAnchor::Shell(b),
        free_length_m: 0.101,
        axial_rigidity_n: 100.,
    };
    let stiffness = strut.axial_rigidity_n / strut.free_length_m;
    let compression = stiffness * (strut.free_length_m - 0.1);
    // Exact opposing dead forces hold the two bearings at their reference
    // positions. This isolates the physical member tangent from shell strain.
    lower.loads[6 * 5 + 2] = compression;
    upper.loads[6 * 5 + 2] = -compression;
    let joint = equilibrate_connected_shell_assembly(
        &[lower.input(), upper.input()],
        &[],
        &[],
        budget(),
        &[],
        &[strut],
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    assert!(joint.tether_responses.is_empty());
    let response = joint.strut_responses[0];
    near(response.compression_n, compression, 1e-12);
    near(response.gap_m, -0.001, 1e-14);
    near(response.stored_energy_j, 0.5 * compression * 0.001, 1e-13);
    near(response.potential_change_j, 0., 1e-13);
    let forces = response.endpoint_forces_n();
    near(forces[0][2], -compression, 1e-12);
    near(forces[1][2], compression, 1e-12);
    for c in 0..3 {
        assert_eq!(forces[0][c] + forces[1][c], 0.);
    }
    assert!(
        joint
            .equilibrium
            .coordinates
            .iter()
            .all(|v| v.abs() < 1e-13)
    );
    assert!(joint.equilibrium.residual_force_n < 1e-9);
    lower.loads.fill(0.);
    upper.loads.fill(0.);
    let bare = equilibrate_shell_assembly(
        &[lower.input(), upper.input()],
        &[],
        &[],
        budget(),
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    let diagonal = [-compression / 0.1, -compression / 0.1, stiffness];
    let model = &joint.equilibrium.assembly.model;
    for (axis, expected) in diagonal.into_iter().enumerate() {
        let i = joint.equilibrium.assembly.map[0][5][axis].unwrap();
        let j = joint.equilibrium.assembly.map[1][5][axis].unwrap();
        near(model.k.get(i, j), -expected, 1e-12);
        near(
            model.k.get(i, i) - bare.assembly.model.k.get(i, i),
            expected,
            1e-8,
        );
        assert_eq!(model.k.get(i, j), model.k.get(j, i));
    }
    for row in 0..model.free {
        assert_eq!(model.m.row(row), bare.assembly.model.m.row(row));
    }
    // A reciprocal tensile connection balances the post in the SAME solve.
    // Their transverse geometric terms cancel; the post's axial stiffness
    // remains. This catches load omission and double-counting at either end.
    let tether = AssemblyTether {
        bearing: a,
        anchor: AssemblyTetherAnchor::Shell(b),
        reference_tension_n: compression,
        axial_stiffness_n_per_m: 0.,
    };
    let mixed = equilibrate_connected_shell_assembly(
        &[lower.input(), upper.input()],
        &[],
        &[],
        budget(),
        &[tether],
        &[strut],
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    assert!(
        mixed
            .equilibrium
            .coordinates
            .iter()
            .all(|v| v.abs() < 1e-13)
    );
    near(
        mixed.tether_responses[0].tension_n,
        mixed.strut_responses[0].compression_n,
        1e-12,
    );
    for axis in 0..3 {
        let i = mixed.equilibrium.assembly.map[0][5][axis].unwrap();
        let j = mixed.equilibrium.assembly.map[1][5][axis].unwrap();
        near(
            mixed.equilibrium.assembly.model.k.get(i, j),
            if axis == 2 { -stiffness } else { 0. },
            1e-11,
        );
    }
}

#[test]
fn separated_short_post_has_no_force_no_tangent_and_can_close() {
    let mut part = Part::rectangle(0, 2, 2, 0., true);
    let bearing = port(&part, 0, 5);
    let short = AssemblyStrut {
        bearing,
        anchor: AssemblyTetherAnchor::World([0.4, 0.1, 0.]),
        free_length_m: 0.19,
        axial_rigidity_n: 2e6,
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
    let open = equilibrate_connected_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        &[],
        &[short],
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    assert_eq!(open.equilibrium.displacements, bare.displacements);
    let response = open.strut_responses[0];
    assert_eq!(response.compression_n, 0.);
    assert_eq!(response.force_n, [0.; 3]);
    assert_eq!(response.potential_change_j, 0.);
    assert_eq!(response.stored_energy_j, 0.);
    near(response.gap_m, 0.01, 1e-14);
    for i in 0..bare.assembly.model.free {
        for j in 0..bare.assembly.model.free {
            assert_eq!(
                open.equilibrium.assembly.model.k.get(i, j),
                bare.assembly.model.k.get(i, j)
            );
        }
    }
    // At first touch the documented open-side tangent is zero. A small load
    // closing that gap must activate the actual EA/L law during Newton.
    let touching = AssemblyStrut {
        free_length_m: 0.2,
        ..short
    };
    part.loads[6 * 5] = 1.;
    let closed = equilibrate_connected_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        &[],
        &[touching],
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    assert!(closed.strut_responses[0].gap_m < 0.);
    assert!(closed.strut_responses[0].compression_n > 0.);
    assert!(closed.strut_responses[0].force_n[0] < 0.);
    assert!(closed.equilibrium.residual_force_n < 1e-8);
    part.loads[6 * 5] = -1.;
    let released = equilibrate_connected_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        &[],
        &[touching],
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    assert!(released.strut_responses[0].gap_m > 0.);
    assert_eq!(released.strut_responses[0].compression_n, 0.);
    assert_eq!(released.strut_responses[0].stored_energy_j, 0.);
    assert_eq!(released.strut_responses[0].force_n, [0.; 3]);
}

#[test]
fn compression_through_nonmatching_ties_preserves_energy_force_and_tangent() {
    let left = Part::rectangle(0, 2, 2, 0., true);
    let right = Part::rectangle(2, 2, 3, 0., false);
    let a = AssemblyTetherPoint {
        arm_m: [0.002, -0.001, 0.004],
        ..port(&right, 1, 3)
    };
    let b = AssemblyTetherPoint {
        arm_m: [-0.001, 0.002, -0.002],
        ..port(&left, 0, 4)
    };
    let chord: [f64; 3] = std::array::from_fn(|c| {
        right.mesh.nodes[3][c] + a.arm_m[c] - left.mesh.nodes[4][c] - b.arm_m[c]
    });
    let length = chord.iter().fold(0.0_f64, |n, c| n.hypot(*c));
    let strut = AssemblyStrut {
        bearing: a,
        anchor: AssemblyTetherAnchor::Shell(b),
        free_length_m: length + 1e-4,
        axial_rigidity_n: 2.,
    };
    let mut checked = false;
    let joint = equilibrate_connected_shell_assembly(
        &[left.input(), right.input()],
        &[seam(3, TieOffset::Midsurface)],
        &[],
        budget(),
        &[],
        &[strut],
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
                assert!(connection.strut_response(&u).unwrap().compression_n > 0.);
                let h = 1e-8;
                let plus: Vec<_> = u.iter().zip(&v).map(|(u, v)| u + h * v).collect();
                let minus: Vec<_> = u.iter().zip(&v).map(|(u, v)| u - h * v).collect();
                let mut gradient = vec![0.; p.dimension()];
                connection.add_gradient(&u, &mut gradient, 1.).unwrap();
                let fd = (connection.strut_response(&plus).unwrap().potential_change_j
                    - connection
                        .strut_response(&minus)
                        .unwrap()
                        .potential_change_j)
                    / (2. * h);
                near(fd, gradient.iter().zip(&v).map(|(a, b)| a * b).sum(), 1e-7);
                let mut matrix = Coo::new(p.dimension(), p.dimension());
                connection.add_tangent(&u, &mut matrix, 1.).unwrap();
                let matrix = matrix.assemble();
                let mut product = vec![0.; p.dimension()];
                matrix.spmv(&v, &mut product);
                let mut action = vec![0.; p.dimension()];
                connection.add_action(&u, &v, &mut action, 1.).unwrap();
                let mut plus_force = vec![0.; p.dimension()];
                let mut minus_force = plus_force.clone();
                connection.add_gradient(&plus, &mut plus_force, 1.).unwrap();
                connection
                    .add_gradient(&minus, &mut minus_force, 1.)
                    .unwrap();
                for i in 0..p.dimension() {
                    near(action[i], product[i], 1e-11);
                    near((plus_force[i] - minus_force[i]) / (2. * h), action[i], 1e-7);
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
    assert!(joint.strut_responses[0].compression_n > 0.);
    assert!(joint.equilibrium.residual_force_n < 1e-8);
}

#[test]
fn invalid_struts_combined_work_caps_and_unstable_compression_refuse() {
    let part = Part::rectangle(0, 2, 2, 0., true);
    let bearing = port(&part, 0, 5);
    let strut = AssemblyStrut {
        bearing,
        anchor: AssemblyTetherAnchor::World([0.4, 0.1, 0.]),
        free_length_m: 0.2,
        axial_rigidity_n: 100.,
    };
    for invalid in [
        AssemblyStrut {
            free_length_m: 0.,
            ..strut
        },
        AssemblyStrut {
            free_length_m: f64::NAN,
            ..strut
        },
        AssemblyStrut {
            axial_rigidity_n: 0.,
            ..strut
        },
        AssemblyStrut {
            axial_rigidity_n: -1.,
            ..strut
        },
        AssemblyStrut {
            axial_rigidity_n: f64::INFINITY,
            ..strut
        },
        AssemblyStrut {
            free_length_m: 1e-9,
            axial_rigidity_n: f64::MAX,
            ..strut
        },
    ] {
        assert!(
            equilibrate_connected_shell_assembly(
                &[part.input()],
                &[],
                &[],
                budget(),
                &[],
                &[invalid],
                options(),
                &mut |_, _, _| panic!("invalid constitutive data reached Newton"),
                || false,
            )
            .is_err()
        );
    }
    // The boundary is shared by BOTH connection kinds, even when the members
    // are open and their supported bearings have zero retained columns.
    let supported = AssemblyStrut {
        bearing: port(&part, 0, 0),
        anchor: AssemblyTetherAnchor::World([0., 0., -0.1]),
        free_length_m: 0.09,
        ..strut
    };
    let at_cap = vec![supported; 4096];
    let admitted = equilibrate_connected_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        &[],
        &at_cap,
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    assert_eq!(admitted.strut_responses.len(), 4096);
    assert!(
        admitted
            .strut_responses
            .iter()
            .all(|s| s.compression_n == 0.)
    );
    let tether = AssemblyTether {
        bearing,
        anchor: strut.anchor,
        reference_tension_n: 0.,
        axial_stiffness_n_per_m: 0.,
    };
    let refused = equilibrate_connected_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        &[tether],
        &at_cap,
        options(),
        &mut |_, _, _| panic!("excess connection count reached Newton"),
        || false,
    )
    .unwrap_err();
    assert!(refused.contains("4096"));
    let free = ShellPreloadComponent {
        boundary: &[],
        support: ShellSupport::Free,
        ..part.input()
    };
    assert!(
        equilibrate_connected_shell_assembly(
            &[free],
            &[],
            &[],
            budget(),
            &[],
            &[strut],
            options(),
            &mut reference_newton,
            || false,
        )
        .unwrap_err()
        .contains("rigid modes")
    );
    // Hold an enormously compressed axial member at exact force balance. The
    // negative transverse curvature MUST survive, so zero residual alone cannot
    // authorize a stable modal pencil. No flexural-member model is claimed.
    let mut lower = Part::rectangle(0, 2, 2, 0., true);
    let mut upper = Part::rectangle(0, 2, 2, 0.1, true);
    let unstable = AssemblyStrut {
        bearing: port(&lower, 0, 5),
        anchor: AssemblyTetherAnchor::Shell(port(&upper, 1, 5)),
        free_length_m: 0.101,
        axial_rigidity_n: 1e12,
    };
    let compression =
        (unstable.axial_rigidity_n / unstable.free_length_m) * (unstable.free_length_m - 0.1);
    lower.loads[6 * 5 + 2] = compression;
    upper.loads[6 * 5 + 2] = -compression;
    let error = equilibrate_connected_shell_assembly(
        &[lower.input(), upper.input()],
        &[],
        &[],
        budget(),
        &[],
        &[unstable],
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap_err();
    assert!(error.contains("unstable"), "{error}");
}

#[test]
fn strut_solver_failure_and_final_cancellation_publish_no_partial_body() {
    let part = Part::rectangle(0, 2, 2, 0., true);
    let original = part.mesh.clone();
    let strut = AssemblyStrut {
        bearing: port(&part, 0, 5),
        anchor: AssemblyTetherAnchor::World([0.4, 0.1, 0.]),
        free_length_m: 0.2001,
        axial_rigidity_n: 100.,
    };
    let mut attempts = 0;
    let error = equilibrate_connected_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        &[],
        &[strut],
        options(),
        &mut |p, x, limit| {
            attempts += 1;
            if attempts == 2 {
                return Err("injected strut Newton failure".into());
            }
            reference_newton(p, x, limit)
        },
        || false,
    )
    .unwrap_err();
    assert!(error.contains("injected strut Newton failure"));
    assert_eq!(part.mesh, original);
    let mut polls = 0;
    equilibrate_connected_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        &[],
        &[strut],
        options(),
        &mut reference_newton,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut repeat = 0;
    let error = equilibrate_connected_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        &[],
        &[strut],
        options(),
        &mut reference_newton,
        || {
            repeat += 1;
            repeat == polls
        },
    )
    .unwrap_err();
    assert_eq!(repeat, polls);
    assert!(error.contains("cancelled"));
    assert_eq!(part.mesh, original);
}
