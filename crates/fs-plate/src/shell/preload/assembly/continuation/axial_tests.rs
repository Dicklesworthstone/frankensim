//! G0/G1/G4: retuning the original axial laws on the same accepted assembly.
use super::*;
use crate::shell::junction::TieOffset;
use crate::shell::preload::assembly::tests::{Part, budget, options, seam};
use crate::shell::preload::assembly::{
    AssemblyTetherAnchor, AssemblyTetherPoint, equilibrate_connected_shell_assembly,
};
use crate::shell::preload::tests::reference_newton;
use fs_sparse::Coo;
use std::cell::Cell;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}; tolerance {tolerance}"
    );
}

fn point(part: &Part, component: usize, node: usize) -> AssemblyTetherPoint {
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

fn panel() -> Part {
    let mut part = Part::rectangle(0, 2, 2, 0., true);
    part.fixed = (0..9).filter(|n| *n != 4).collect();
    part
}

fn connections(part: &Part) -> (AssemblyTether, AssemblyStrut) {
    let bearing = point(part, 0, 4);
    (
        AssemblyTether {
            bearing,
            anchor: AssemblyTetherAnchor::World([0.1, 0.1, 0.15]),
            reference_tension_n: 0.3,
            axial_stiffness_n_per_m: 75.,
        },
        AssemblyStrut {
            bearing,
            anchor: AssemblyTetherAnchor::World([0.1, 0.1, -0.1]),
            free_length_m: 0.10002,
            axial_rigidity_n: 30.,
        },
    )
}

fn make(part: &Part, tether: AssemblyTether, strut: AssemblyStrut) -> ShellAssemblyLoadPath {
    ShellAssemblyLoadPath::new(
        &[part.input()],
        &[],
        &[],
        budget(),
        &[tether],
        &[strut],
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap()
}

fn retained(body: &ShellAssemblyLoadPath) -> (String, Vec<Vec<f64>>, Vec<f64>, Vec<f64>) {
    (
        format!("{:?}", body.accepted()),
        body.component_loads().to_vec(),
        body.tether_reference_tensions_n().to_vec(),
        body.strut_free_lengths_m().to_vec(),
    )
}

fn same_equilibrium(a: &PreloadedConnectedShellAssembly, b: &PreloadedConnectedShellAssembly) {
    for (a, b) in a
        .equilibrium
        .coordinates
        .iter()
        .zip(&b.equilibrium.coordinates)
    {
        near(*a, *b, 2e-10);
    }
    let ak = &a.equilibrium.assembly.model.k;
    let bk = &b.equilibrium.assembly.model.k;
    let n = a.equilibrium.assembly.model.free;
    let mut difference = 0.0_f64;
    let mut scale = 0.0_f64;
    for i in 0..n {
        difference = difference.max((0..n).map(|j| (ak.get(i, j) - bk.get(i, j)).abs()).sum());
        scale = scale.max((0..n).map(|j| bk.get(i, j).abs()).sum());
        assert_eq!(
            a.equilibrium.assembly.model.m.row(i),
            b.equilibrium.assembly.model.m.row(i)
        );
    }
    assert!(difference <= 3e-9 * scale);
    for (a, b) in a.tether_responses.iter().zip(&b.tether_responses) {
        near(a.tension_n, b.tension_n, 2e-8);
        for (a, b) in a.force_n.iter().zip(b.force_n) {
            near(*a, b, 2e-8);
        }
    }
    for (a, b) in a.strut_responses.iter().zip(&b.strut_responses) {
        near(a.compression_n, b.compression_n, 2e-8);
        near(a.stored_energy_j, b.stored_energy_j, 2e-10);
    }
}

#[test]
fn g1_axial_actuation_matches_independent_panel_force_and_tangent_through_release() {
    let mut part = panel();
    let (mut tether, mut strut) = connections(&part);
    let mut body = make(&part, tether, strut);
    let original = format!("{:?}", body.reference_assembly());
    let d = body.reference_assembly().model.dof_map[26].unwrap();
    let linear = body.reference_assembly().model.k.get(d, d);
    // Independently integrate the isotropic Green--Lagrange membrane energy:
    // U(w)=K_b*w²/2 + C*w⁴/4, C=Eh/(2*(1-nu²))*sum(A*|grad N_center|⁴).
    let mut geometric = 0.;
    for tri in &part.mesh.tris {
        let Some(a) = tri.iter().position(|n| *n == 4) else {
            continue;
        };
        let p = part.mesh.nodes[tri[a]];
        let q = part.mesh.nodes[tri[(a + 1) % 3]];
        let r = part.mesh.nodes[tri[(a + 2) % 3]];
        let det = (q[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (q[1] - p[1]);
        let gx = (q[1] - r[1]) / det;
        let gy = (r[0] - q[0]) / det;
        geometric += 0.5 * det.abs() * (gx * gx + gy * gy).powi(2);
    }
    let cubic = 9e9 * 0.004 / (2. * (1. - 0.27_f64.powi(2))) * geometric;
    for (tension, free_length, dead) in
        [(0.8, 0.1004, 0.1), (0., 0.0999, 0.2), (0.6, 0.10005, -0.15)]
    {
        part.loads[26] = dead;
        tether.reference_tension_n = tension;
        strut.free_length_m = free_length;
        let stiffness = strut.axial_rigidity_n / free_length;
        let force = |w: f64| {
            linear * w + cubic * w.powi(3)
                - dead
                - (tension - 75. * w).max(0.)
                - stiffness * (free_length - 0.1 - w).max(0.)
        };
        let (mut lo, mut hi) = (-0.001, 0.001);
        assert!(force(lo) < 0. && force(hi) > 0.);
        for _ in 0..60 {
            let middle = 0.5 * (lo + hi);
            if force(middle) < 0. {
                lo = middle;
            } else {
                hi = middle;
            }
        }
        let w = 0.5 * (lo + hi);
        let expected_tension = (tension - 75. * w).max(0.);
        let shortening = (free_length - 0.1 - w).max(0.);
        let expected_compression = stiffness * shortening;
        let tangent = linear
            + 3. * cubic * w * w
            + if expected_tension > 0. { 75. } else { 0. }
            + if expected_compression > 0. {
                stiffness
            } else {
                0.
            };
        let state = body
            .advance_with_axial_settings(
                &[&part.loads],
                &[tension],
                &[free_length],
                options(),
                &mut reference_newton,
                || false,
            )
            .unwrap();
        near(state.equilibrium.displacements[0][26], w, 2e-12);
        near(state.equilibrium.assembly.model.k.get(d, d), tangent, 2e-5);
        near(
            state.equilibrium.stored_energy_j,
            0.5 * linear * w * w + 0.25 * cubic * w.powi(4),
            2e-10,
        );
        near(state.tether_responses[0].tension_n, expected_tension, 2e-9);
        near(
            state.strut_responses[0].compression_n,
            expected_compression,
            2e-9,
        );
        near(
            state.strut_responses[0].stored_energy_j,
            0.5 * stiffness * shortening * shortening,
            2e-12,
        );
        if tension == 0. {
            assert_eq!(state.tether_responses[0].tension_n, 0.);
            assert_eq!(state.strut_responses[0].compression_n, 0.);
        }
        let cold = equilibrate_connected_shell_assembly(
            &[part.input()],
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
        same_equilibrium(body.accepted(), &cold);
        assert_eq!(body.tether_reference_tensions_n(), &[tension]);
        assert_eq!(body.strut_free_lengths_m(), &[free_length]);
        assert_eq!(format!("{:?}", body.reference_assembly()), original);
    }
}

#[test]
fn g1_tied_parameter_path_uses_original_intermediate_laws_and_full_tangents() {
    let left = Part::rectangle(0, 2, 2, 0., true);
    let mut right = Part::rectangle(2, 2, 4, 0., false);
    let bearing = point(&right, 1, 8);
    let tether = AssemblyTether {
        bearing,
        anchor: AssemblyTetherAnchor::World([0.45, 0.1, -0.04]),
        reference_tension_n: 0.5,
        axial_stiffness_n_per_m: 20.,
    };
    let strut = AssemblyStrut {
        bearing,
        anchor: AssemblyTetherAnchor::World([0.4, 0.1, 0.1]),
        free_length_m: 0.10002,
        axial_rigidity_n: 100.,
    };
    let mut one = ShellAssemblyLoadPath::new(
        &[left.input(), right.input()],
        &[seam(4, TieOffset::Midsurface)],
        &[],
        budget(),
        &[tether],
        &[strut],
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    let mut split = ShellAssemblyLoadPath::new(
        &[left.input(), right.input()],
        &[seam(4, TieOffset::Midsurface)],
        &[],
        budget(),
        &[tether],
        &[strut],
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    let original = format!("{:?}", one.reference_assembly());
    let initial_coordinates = one.accepted().equilibrium.coordinates.clone();
    let initial_loads = right.loads.clone();
    for f in [0.25, 0.5, 0.75, 1.] {
        right.loads[50] = -0.35 * f;
        right.loads[38] = -0.1 * f; // eliminated nonmatching seam coordinate
        split
            .advance_with_axial_settings(
                &[&left.loads, &right.loads],
                &[0.5 + 0.4 * f],
                &[0.10002 + 0.00007 * f],
                options(),
                &mut reference_newton,
                || false,
            )
            .unwrap();
    }
    let mut increments = 0;
    let mut solve_options = options();
    solve_options.load_steps = 4;
    one.advance_with_axial_settings(
        &[&left.loads, &right.loads],
        &[0.9],
        &[0.10009],
        solve_options,
        &mut |p, x, limit| {
            increments += 1;
            let fraction = increments as f64 / 4.;
            near(p.fraction, fraction, 0.);
            let intermediate_tether = AssemblyTether {
                reference_tension_n: (1. - fraction) * 0.5 + fraction * 0.9,
                ..tether
            };
            let intermediate_strut = AssemblyStrut {
                free_length_m: (1. - fraction) * 0.10002 + fraction * 0.10009,
                ..strut
            };
            let expected = prepare_assembly(
                &[left.input(), right.input()],
                &[seam(4, TieOffset::Midsurface)],
                &[],
                budget(),
                &[intermediate_tether],
                &[intermediate_strut],
                &mut || false,
            )?;
            // At the undeformed probe the post is compressed. This distinguishes
            // EA/L_interpolated from a blend of endpoint stiffness/potentials,
            // even if the loaded endpoint subsequently releases the post.
            let u = vec![0.; p.dimension()];
            let v: Vec<_> = (0..p.dimension())
                .map(|i| ((i * 11 % 13) as f64 - 6.) / 6.)
                .collect();
            let mut gradient = vec![0.; p.dimension()];
            let mut original_gradient = gradient.clone();
            let mut action = gradient.clone();
            let mut matrix = Coo::new(p.dimension(), p.dimension());
            let mut original_matrix = Coo::new(p.dimension(), p.dimension());
            for (connection, original) in p.tethers.iter().zip(&expected.tethers) {
                connection.add_gradient_at(&u, &mut gradient, 1., fraction)?;
                original.add_gradient(&u, &mut original_gradient, 1.)?;
                connection.add_action_at(&u, &v, &mut action, 1., fraction)?;
                connection.add_tangent_at(&u, &mut matrix, 1., fraction)?;
                original.add_tangent(&u, &mut original_matrix, 1.)?;
            }
            let matrix = matrix.assemble();
            let original_matrix = original_matrix.assemble();
            let mut product = vec![0.; p.dimension()];
            matrix.spmv(&v, &mut product);
            let h = 1e-8;
            let plus: Vec<_> = v.iter().map(|v| h * v).collect();
            let minus: Vec<_> = v.iter().map(|v| -h * v).collect();
            let mut plus_force = vec![0.; p.dimension()];
            let mut minus_force = plus_force.clone();
            for connection in p.tethers {
                connection.add_gradient_at(&plus, &mut plus_force, 1., fraction)?;
                connection.add_gradient_at(&minus, &mut minus_force, 1., fraction)?;
            }
            for i in 0..p.dimension() {
                near(gradient[i], original_gradient[i], 2e-10);
                near(action[i], product[i], 1e-10 * (1. + action[i].abs()));
                near(
                    (plus_force[i] - minus_force[i]) / (2. * h),
                    action[i],
                    2e-6 * (1. + action[i].abs()),
                );
                for j in 0..p.dimension() {
                    near(
                        matrix.get(i, j),
                        original_matrix.get(i, j),
                        2e-10 * (1. + matrix.get(i, j).abs()),
                    );
                    assert_eq!(matrix.get(i, j), matrix.get(j, i));
                }
            }
            reference_newton(p, x, limit)
        },
        || false,
    )
    .unwrap();
    assert_eq!(increments, 4);
    let cold = equilibrate_connected_shell_assembly(
        &[left.input(), right.input()],
        &[seam(4, TieOffset::Midsurface)],
        &[],
        budget(),
        &[AssemblyTether {
            reference_tension_n: 0.9,
            ..tether
        }],
        &[AssemblyStrut {
            free_length_m: 0.10009,
            ..strut
        }],
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    same_equilibrium(one.accepted(), &cold);
    same_equilibrium(split.accepted(), &cold);
    one.advance_with_axial_settings(
        &[&left.loads, &initial_loads],
        &[0.5],
        &[0.10002],
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    for (a, b) in one
        .accepted()
        .equilibrium
        .coordinates
        .iter()
        .zip(initial_coordinates)
    {
        near(*a, b, 2e-10);
    }
    assert_eq!(format!("{:?}", one.reference_assembly()), original);
}

#[test]
fn g0_setting_noop_and_refusals_preserve_the_complete_accepted_body() {
    let part = panel();
    let (tether, strut) = connections(&part);
    let mut body = make(&part, tether, strut);
    let before = retained(&body);
    body.advance_with_axial_settings(
        &[&part.loads],
        &[0.3],
        &[0.10002],
        options(),
        &mut |_, _, _| panic!("exact unchanged settings and loads must not solve"),
        || false,
    )
    .unwrap();
    assert_eq!(retained(&body), before);
    for (tensions, lengths) in [
        (vec![], vec![0.10002]),
        (vec![0.3, 0.3], vec![0.10002]),
        (vec![0.3], vec![]),
        (vec![0.3], vec![0.10002, 0.10002]),
        (vec![-1.], vec![0.10002]),
        (vec![f64::NAN], vec![0.10002]),
        (vec![0.5], vec![0.]),
        (vec![0.5], vec![f64::INFINITY]),
        (vec![0.5], vec![1e-9_f64.next_down()]),
    ] {
        assert!(
            body.advance_with_axial_settings(
                &[&part.loads],
                &tensions,
                &lengths,
                options(),
                &mut |_, _, _| panic!("invalid settings must refuse before solving"),
                || false
            )
            .is_err()
        );
        assert_eq!(retained(&body), before);
    }
    let mut target = part.loads.clone();
    target[26] = 0.2;
    assert!(
        body.advance_with_axial_settings(
            &[&target],
            &[0.8],
            &[0.1001],
            options(),
            &mut |_, _, _| Err("injected solver failure".into()),
            || false
        )
        .is_err()
    );
    let mut too_small = options();
    too_small.max_evaluations = 1;
    assert!(
        body.advance_with_axial_settings(
            &[&target],
            &[0.8],
            &[0.1001],
            too_small,
            &mut reference_newton,
            || false
        )
        .is_err()
    );
    assert_eq!(retained(&body), before);
    // At the admitted floor the original post simply opens; its EA remains
    // unchanged, and no stiffness is added on the released side.
    body.advance_with_axial_settings(
        &[&part.loads],
        &[0.3],
        &[1e-9],
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    assert_eq!(body.strut_free_lengths_m(), &[1e-9]);
    assert_eq!(body.accepted().strut_responses[0].compression_n, 0.);
    let mut huge = strut;
    huge.free_length_m = 0.09;
    huge.axial_rigidity_n = 1e300;
    let mut huge_body = make(
        &part,
        AssemblyTether {
            reference_tension_n: 0.,
            ..tether
        },
        huge,
    );
    let huge_before = retained(&huge_body);
    assert!(
        huge_body
            .advance_with_axial_settings(
                &[&part.loads],
                &[0.],
                &[1e-9],
                options(),
                &mut |_, _, _| panic!("unrepresentable EA/L must refuse before solving"),
                || false
            )
            .unwrap_err()
            .contains("stiffness is not representable")
    );
    assert_eq!(retained(&huge_body), huge_before);
}

#[test]
fn g4_final_parameter_cancellation_is_atomic_and_retry_keeps_original_sources() {
    let mut part = panel();
    let (tether, strut) = connections(&part);
    let mut body = make(&part, tether, strut);
    let mut control = make(&part, tether, strut);
    let old = retained(&body);
    let mut target = part.loads.clone();
    target[26] = 0.4;
    let polls = Cell::new(0);
    control
        .advance_with_axial_settings(
            &[&target],
            &[0.8],
            &[0.1001],
            options(),
            &mut reference_newton,
            || {
                polls.set(polls.get() + 1);
                false
            },
        )
        .unwrap();
    let last = polls.get();
    polls.set(0);
    assert!(
        body.advance_with_axial_settings(
            &[&target],
            &[0.8],
            &[0.1001],
            options(),
            &mut reference_newton,
            || {
                polls.set(polls.get() + 1);
                polls.get() == last
            }
        )
        .is_err()
    );
    assert_eq!(polls.get(), last);
    assert_eq!(retained(&body), old);
    part.mesh.nodes[4][2] = 1.;
    part.loads.fill(1e9);
    body.advance_with_axial_settings(
        &[&target],
        &[0.8],
        &[0.1001],
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    assert_eq!(retained(&body), retained(&control));
    // Ordinary dead-load continuation must retain the newly accepted settings.
    target[26] = 0.6;
    body.advance(&[&target], options(), &mut reference_newton, || false)
        .unwrap();
    control
        .advance_with_axial_settings(
            &[&target],
            &[0.8],
            &[0.1001],
            options(),
            &mut reference_newton,
            || false,
        )
        .unwrap();
    assert_eq!(retained(&body), retained(&control));
}
