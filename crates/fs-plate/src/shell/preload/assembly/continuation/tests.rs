use super::*;
use crate::shell::junction::TieOffset;
use crate::shell::preload::assembly::tests::{Part, budget, options, seam};
use crate::shell::preload::assembly::{
    AssemblyTetherAnchor, AssemblyTetherPoint, equilibrate_connected_shell_assembly,
};
use crate::shell::preload::tests::reference_newton;
use std::cell::Cell;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} ({tolerance})"
    );
}

fn path(
    parts: &[ShellPreloadComponent<'_>],
    tethers: &[AssemblyTether],
    struts: &[AssemblyStrut],
) -> ShellAssemblyLoadPath {
    ShellAssemblyLoadPath::new(
        parts,
        &[],
        &[],
        budget(),
        tethers,
        struts,
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap()
}

fn center_panel() -> Part {
    let mut part = Part::rectangle(0, 2, 2, 0., true);
    part.fixed = (0..9).filter(|n| *n != 4).collect();
    part
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

#[test]
fn g1_nonlinear_panel_load_reversal_matches_independent_cubic_energy() {
    let mut part = center_panel();
    let mut body = path(&[part.input()], &[], &[]);
    let original = format!("{:?}", body.reference_assembly());
    let d = body.reference_assembly().model.dof_map[26].unwrap();
    let linear = body.reference_assembly().model.k.get(d, d);
    // Symmetry leaves only the centre's transverse translation. Independently
    // integrating isotropic Green--Lagrange membrane storage gives
    // U(w) = K_b*w²/2 + w⁴ * sum(A*|grad N|⁴)*Eh/(8*(1-nu²)).
    // The linear bending coefficient belongs to the unchanged DKT owner.
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
    for w in [0.0002, 0.0006, -0.0004, 0.] {
        part.loads[26] = linear * w + cubic * w * w * w;
        let state = body
            .advance(&[&part.loads], options(), &mut reference_newton, || false)
            .unwrap();
        near(state.equilibrium.displacements[0][26], w, 2e-12);
        near(
            state.equilibrium.stored_energy_j,
            0.5 * linear * w * w + 0.25 * cubic * w.powi(4),
            2e-10,
        );
        near(
            state.equilibrium.assembly.model.k.get(d, d),
            linear + 3. * cubic * w * w,
            2e-5,
        );
        assert!(
            state.equilibrium.residual_force_n
                <= 1e-8 + options().relative_force_tolerance * part.loads[26].abs()
        );
        for (i, value) in state.equilibrium.displacements[0].iter().enumerate() {
            if i != 26 {
                assert!(value.abs() < 2e-12, "unexpected coordinate {i}: {value}");
            }
        }
        assert_eq!(format!("{:?}", body.reference_assembly()), original);
    }
}

#[test]
fn g1_tied_path_partition_matches_cold_equilibrium_with_original_connections() {
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
    let make = || {
        ShellAssemblyLoadPath::new(
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
        .unwrap()
    };
    let mut one = make();
    let mut split = make();
    let initial_tension = one.accepted().tether_responses[0].tension_n;
    for fraction in [0.25, 0.5, 0.75, 1.] {
        right.loads[6 * 8 + 2] = -fraction;
        right.loads[6 * 6 + 2] = -0.25 * fraction; // eliminated seam load
        split
            .advance(
                &[&left.loads, &right.loads],
                options(),
                &mut reference_newton,
                || false,
            )
            .unwrap();
    }
    let accepted_coordinates = one.accepted().equilibrium.coordinates.clone();
    let accepted_tangent = one.accepted().equilibrium.assembly.model.k.clone();
    let mut first = true;
    one.advance(
        &[&left.loads, &right.loads],
        options(),
        &mut |problem, x, limit| {
            if first {
                first = false;
                let u = problem.physical(&x);
                for (a, b) in u.iter().zip(&accepted_coordinates) {
                    near(*a, *b, 1e-16);
                }
                // Changing dead forces has no load stiffness. The first trial
                // must therefore retain the complete previously accepted
                // connection tangent, even at the first fraction of the ramp.
                let tangent = problem.tangent(&u)?;
                for row in 0..problem.dimension() {
                    let (columns, values) = tangent.row(row);
                    for (&column, &value) in columns.iter().zip(values) {
                        near(value, accepted_tangent.get(row, column), 1e-6);
                    }
                }
            }
            reference_newton(problem, x, limit)
        },
        || false,
    )
    .unwrap();
    let cold = equilibrate_connected_shell_assembly(
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
    for result in [one.accepted(), split.accepted()] {
        let maximum_displacement_difference = result
            .equilibrium
            .coordinates
            .iter()
            .zip(&cold.equilibrium.coordinates)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        let mut tangent_difference = 0.0_f64;
        let mut tangent_scale = 0.0_f64;
        for (a, b) in result
            .equilibrium
            .coordinates
            .iter()
            .zip(&cold.equilibrium.coordinates)
        {
            near(*a, *b, 2e-10);
        }
        for i in 0..cold.equilibrium.assembly.model.free {
            let (_, a) = result.equilibrium.assembly.model.k.row(i);
            let (_, b) = cold.equilibrium.assembly.model.k.row(i);
            tangent_difference =
                tangent_difference.max(a.iter().zip(b).map(|(a, b)| (a - b).abs()).sum::<f64>());
            tangent_scale = tangent_scale.max(b.iter().map(|b| b.abs()).sum::<f64>());
            assert_eq!(
                result.equilibrium.assembly.model.m.row(i),
                cold.equilibrium.assembly.model.m.row(i)
            );
        }
        println!(
            "G1 tied continuation: displacement_difference={maximum_displacement_difference:e}, tangent_relative_difference={:e}, warm_force_residual={:e}, cold_force_residual={:e}",
            tangent_difference / tangent_scale,
            result.equilibrium.residual_force_n,
            cold.equilibrium.residual_force_n
        );
        // Independently converged nonlinear endpoints need not have identical
        // last Newton roundoff. Compare the complete tangent in its physical
        // operator scale, while retaining the original unscaled force gate.
        assert!(tangent_difference <= 1e-9 * tangent_scale);
        near(
            result.tether_responses[0].tension_n,
            cold.tether_responses[0].tension_n,
            1e-9,
        );
        near(
            result.strut_responses[0].compression_n,
            cold.strut_responses[0].compression_n,
            1e-9,
        );
    }
    assert!((one.accepted().tether_responses[0].tension_n - initial_tension).abs() > 1e-6);
}

#[test]
fn g0_first_endpoint_matches_one_shot_and_rest_keeps_every_bit() {
    let mut part = center_panel();
    part.loads[26] = 0.2;
    let mut body = path(&[part.input()], &[], &[]);
    let direct = equilibrate_connected_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        &[],
        &[],
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    assert_eq!(format!("{:?}", body.accepted()), format!("{direct:?}"));
    let before = format!("{:?}", body.accepted());
    body.advance(
        &[&part.loads],
        options(),
        &mut |_, _, _| panic!("unchanged load must not solve"),
        || false,
    )
    .unwrap();
    assert_eq!(format!("{:?}", body.accepted()), before);
    assert_eq!(body.component_loads(), &[part.loads]);
}

#[test]
fn g0_failed_load_changes_preserve_all_accepted_fields_and_retry() {
    let mut part = center_panel();
    part.loads[26] = 0.2;
    let mut body = path(&[part.input()], &[], &[]);
    let mut clean = path(&[part.input()], &[], &[]);
    let accepted = format!("{:?}", body.accepted());
    let loads = body.component_loads().to_vec();
    let mut target = part.loads.clone();
    target[26] = 0.4;
    assert!(
        body.advance(&[], options(), &mut reference_newton, || false)
            .is_err()
    );
    let mut invalid = target.clone();
    invalid[26] = f64::NAN;
    assert!(
        body.advance(&[&invalid], options(), &mut reference_newton, || false)
            .is_err()
    );
    assert!(
        body.advance(
            &[&target],
            options(),
            &mut |_, _, _| Err("injected solve failure".into()),
            || false
        )
        .is_err()
    );
    assert!(
        body.advance(
            &[&target],
            PreloadOptions {
                max_evaluations: 1,
                ..options()
            },
            &mut reference_newton,
            || false
        )
        .is_err()
    );
    let mut outside = target.clone();
    outside[26] = 1e12;
    assert!(
        body.advance(&[&outside], options(), &mut reference_newton, || false)
            .is_err()
    );
    assert_eq!(format!("{:?}", body.accepted()), accepted);
    assert_eq!(body.component_loads(), loads);
    body.advance(&[&target], options(), &mut reference_newton, || false)
        .unwrap();
    clean
        .advance(&[&target], options(), &mut reference_newton, || false)
        .unwrap();
    assert_eq!(
        format!("{:?}", body.accepted()),
        format!("{:?}", clean.accepted())
    );
    assert_eq!(body.component_loads(), clean.component_loads());
}

#[test]
fn g4_final_cancellation_is_atomic_and_replay_retains_original_sources() {
    let mut part = center_panel();
    part.loads[26] = 0.2;
    let mut body = path(&[part.input()], &[], &[]);
    let mut control = path(&[part.input()], &[], &[]);
    let old_state = format!("{:?}", body.accepted());
    let old_loads = body.component_loads().to_vec();
    let reference = format!("{:?}", body.reference_assembly());
    let mut target = part.loads.clone();
    target[26] = 0.4;
    let polls = Cell::new(0);
    control
        .advance(&[&target], options(), &mut reference_newton, || {
            polls.set(polls.get() + 1);
            false
        })
        .unwrap();
    let last = polls.get();
    assert!(last > 2);
    polls.set(0);
    assert!(
        body.advance(&[&target], options(), &mut reference_newton, || {
            polls.set(polls.get() + 1);
            polls.get() == last
        })
        .is_err()
    );
    assert_eq!(polls.get(), last);
    assert_eq!(format!("{:?}", body.accepted()), old_state);
    assert_eq!(body.component_loads(), old_loads);
    // Caller geometry and initial load arrays are no longer authoritative:
    // prepared material kernels, tie maps and reference mass remain owned.
    part.mesh.nodes[4][2] = 1.;
    part.loads.fill(1e9);
    body.advance(&[&target], options(), &mut reference_newton, || false)
        .unwrap();
    assert_eq!(format!("{:?}", body.reference_assembly()), reference);
    assert_eq!(
        format!("{:?}", body.accepted()),
        format!("{:?}", control.accepted())
    );
}
