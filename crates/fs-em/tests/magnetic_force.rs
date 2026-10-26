//! G1/G3/G4: field-derived force, material substitution and native mechanics.
//! Synthetic source laws verify the numerical connection, not device accuracy.
mod common;

use common::{config, limits, material, near};
use fs_em::{
    force::MaxwellForceLimits,
    planar::{
        MagneticError, PlanarMagnetostatics,
        eddy::{EddyControls, PlanarEddy},
    },
    transfer::PlanarFrame,
    winding::{PlanarWindingBasis, WindingConfig},
};
use fs_ga::{Mat3, Se3, Vec3, Wrench};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
    free_floating_forward_dynamics,
};
use fs_rep_mesh::{Metric2, TriComplex2};

const DEPTH: f64 = 0.2;
const AMBIENT_MU: f64 = 1e-6;

fn force_limits() -> MaxwellForceLimits {
    MaxwellForceLimits {
        max_cells: 10000,
        max_edges: 20000,
    }
}
fn frame() -> PlanarFrame {
    PlanarFrame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap()
}

// One square material inclusion, with a square integration contour in its
// homogeneous exterior collar. No geometry or solver changes when mu changes.
fn specimen(
    n: usize,
    inner_mu: f64,
    reverse: bool,
    depth: f64,
) -> (PlanarMagnetostatics, Vec<usize>, Vec<usize>) {
    let mut points = Vec::new();
    for y in 0..=n {
        for x in 0..=n {
            points.push([
                -1.0 + 2.0 * x as f64 / n as f64,
                -1.0 + 2.0 * y as f64 / n as f64,
            ]);
        }
    }
    let mut faces = Vec::new();
    for y in 0..n {
        for x in 0..n {
            let a = (y * (n + 1) + x) as u32;
            let b = a + 1;
            let d = a + (n + 1) as u32;
            let c = d + 1;
            faces.extend(if reverse {
                [[a, c, b], [a, d, c]]
            } else {
                [[a, b, c], [a, c, d]]
            });
        }
    }
    let mut enclosed = Vec::new();
    let mut body = Vec::new();
    let ambient = material(AMBIENT_MU);
    let inside = material(inner_mu);
    let mut sources = Vec::new();
    for (cell, face) in faces.iter().enumerate() {
        let center: [f64; 2] =
            [0, 1].map(|axis| face.iter().map(|&v| points[v as usize][axis]).sum::<f64>() / 3.0);
        let radius = center[0].abs().max(center[1].abs());
        if radius < 0.5 {
            enclosed.push(cell);
        }
        if radius < 0.25 {
            body.push(cell);
            sources.push(inside.clone());
        } else {
            sources.push(ambient.clone());
        }
    }
    let mesh = TriComplex2::from_indexed_triangles(
        "magnetic-force-square",
        points,
        faces,
        Metric2::planar(depth).unwrap(),
    )
    .unwrap();
    (
        PlanarMagnetostatics::new(mesh, sources, limits(), &mut || false).unwrap(),
        enclosed,
        body,
    )
}

fn boundary(
    problem: &PlanarMagnetostatics,
    potential: impl Fn([f64; 2]) -> f64,
) -> Vec<(usize, f64)> {
    problem
        .mesh()
        .boundary_trace()
        .unwrap()
        .vertices
        .into_iter()
        .map(|v| (v as usize, potential(problem.mesh().vertices()[v as usize])))
        .collect()
}

fn current_load(n: usize, scale: f64, depth: f64) -> Wrench {
    let (problem, enclosed, body) = specimen(n, AMBIENT_MU, false, depth);
    let mut current = vec![0.0; problem.mesh().faces().len()];
    for cell in body {
        current[cell] = scale * 1000.0;
    }
    let bc = boundary(&problem, |[x, _]| scale * 1e-3 * x);
    let field = problem
        .solve(&current, &bc, config(), &mut || false)
        .unwrap();
    assert_eq!(
        field.impressed_current_density_a_m2(),
        Some(current.as_slice())
    );
    field
        .maxwell_wrench(
            &enclosed,
            &frame(),
            [0.0, 0.0, 0.5 * depth],
            force_limits(),
            &mut || false,
        )
        .unwrap()
}

#[test]
fn uniform_field_has_zero_net_load_for_either_mesh_orientation_and_selection_order() {
    for reverse in [false, true] {
        let (problem, mut enclosed, _) = specimen(8, AMBIENT_MU, reverse, DEPTH);
        let bc = boundary(&problem, |[x, y]| 1e-3 * (2.0 * x + y));
        let field = problem
            .solve(
                &vec![0.0; problem.mesh().faces().len()],
                &bc,
                config(),
                &mut || false,
            )
            .unwrap();
        let first = field
            .maxwell_wrench(
                &enclosed,
                &frame(),
                [0.1, 0.2, 0.3],
                force_limits(),
                &mut || false,
            )
            .unwrap();
        for value in first.to_array() {
            assert!(value.abs() < 1e-12, "{first:?}");
        }
        enclosed.reverse();
        assert_eq!(
            first,
            field
                .maxwell_wrench(
                    &enclosed,
                    &frame(),
                    [0.1, 0.2, 0.3],
                    force_limits(),
                    &mut || false
                )
                .unwrap()
        );
    }
}

#[test]
fn exterior_stress_converges_to_independent_lorentz_force_and_scales_with_depth() {
    // The centered current distribution has zero net self-force by symmetry.
    // Its external B=(0,-1e-3,0) gives J cross B=(1,0,0) N/m3.
    // Conductor area is 0.5^2 m2, so F_x=0.25*depth N, without a field-force call.
    let exact = 0.25 * DEPTH;
    let mut previous = f64::INFINITY;
    for n in [8, 16, 32] {
        let load = current_load(n, 1.0, DEPTH);
        let error = (load.force.x - exact).abs();
        assert!(
            error < 0.6 * previous,
            "n={n}, error={error}, previous={previous}"
        );
        assert!(load.force.y.abs() < 1e-10 && load.torque.z.abs() < 1e-10);
        previous = error;
    }
    assert!(previous < 0.05 * exact);
    let baseline = current_load(16, 1.0, DEPTH);
    let deeper = current_load(16, 1.0, 3.0 * DEPTH);
    let stronger = current_load(16, 2.0, DEPTH);
    near(deeper.force.x, 3.0 * baseline.force.x);
    near(stronger.force.x, 4.0 * baseline.force.x);
}

#[test]
fn material_only_substitution_changes_computed_force_and_native_body_acceleration() {
    let inertia = SpatialInertia::new(
        2.0,
        Vec3::new(0.0, 0.0, 0.0),
        Mat3 {
            m: [0.2, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.2],
        },
    )
    .unwrap();
    let mechanics = ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia,
    )])
    .unwrap();
    let mut acceleration = Vec::new();
    for inner_mu in [AMBIENT_MU, 2.0 * AMBIENT_MU, 4.0 * AMBIENT_MU] {
        let (problem, enclosed, _) = specimen(16, inner_mu, false, DEPTH);
        let bc = boundary(&problem, |[x, y]| 1e-3 * x + 4e-4 * (x * x - y * y));
        let field = problem
            .solve(
                &vec![0.0; problem.mesh().faces().len()],
                &bc,
                config(),
                &mut || false,
            )
            .unwrap();
        let body_frame =
            PlanarFrame::new([0.0, 0.0, -DEPTH / 2.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
        let load = field
            .maxwell_wrench(
                &enclosed,
                &body_frame,
                [0.0; 3],
                force_limits(),
                &mut || false,
            )
            .unwrap();
        // No force coefficient, material switch or alternate integrator intervenes.
        let dynamics = free_floating_forward_dynamics(
            &mechanics,
            FreeFloatingBaseState::stationary(Se3::identity()),
            &[],
            &[],
            &[],
            Vec3::new(0.0, 0.0, 0.0),
            &[load],
        )
        .unwrap();
        near(
            dynamics.base_spatial_acceleration_body.linear.x,
            load.force.x / 2.0,
        );
        near(
            dynamics.base_spatial_acceleration_body.linear.y,
            load.force.y / 2.0,
        );
        near(
            dynamics.base_spatial_acceleration_body.angular.z,
            load.torque.z / 0.2,
        );
        acceleration.push(dynamics.base_spatial_acceleration_body.linear.x);
    }
    assert!(acceleration[0].abs() < 1e-10);
    assert!(acceleration[1] > 1e-3 && acceleration[2] > 1.5 * acceleration[1]);
}

#[test]
fn mechanical_frame_rotation_and_moment_origin_preserve_the_actual_wrench() {
    let (problem, enclosed, body) = specimen(16, AMBIENT_MU, false, DEPTH);
    let mut current = vec![0.0; problem.mesh().faces().len()];
    for cell in body {
        current[cell] = 1000.0;
    }
    let bc = boundary(&problem, |[x, _]| 1e-3 * x);
    let field = problem
        .solve(&current, &bc, config(), &mut || false)
        .unwrap();
    let local = field
        .maxwell_wrench(
            &enclosed,
            &frame(),
            [0.0, 0.0, DEPTH / 2.0],
            force_limits(),
            &mut || false,
        )
        .unwrap();
    let moved = PlanarFrame::new([2.0, 3.0, 4.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]).unwrap();
    let origin = [2.0 + DEPTH / 2.0 + 0.2, 3.3, 4.4];
    let global = field
        .maxwell_wrench(&enclosed, &moved, origin, force_limits(), &mut || false)
        .unwrap();
    let force = Vec3::new(local.force.z, local.force.x, local.force.y);
    let torque = Vec3::new(local.torque.z, local.torque.x, local.torque.y)
        + Vec3::new(-0.2, -0.3, -0.4).cross(force);
    for (got, expected) in global
        .to_array()
        .into_iter()
        .zip(Wrench::new(torque, force).to_array())
    {
        near(got, expected);
    }
}

#[test]
fn unsupported_contours_and_cancellation_cannot_publish_partial_loads() {
    let (problem, enclosed, body) = specimen(8, 2.0 * AMBIENT_MU, false, DEPTH);
    let bc = boundary(&problem, |[x, _]| 1e-3 * x);
    let current = vec![0.0; problem.mesh().faces().len()];
    let field = problem
        .solve(&current, &bc, config(), &mut || false)
        .unwrap();
    let baseline = field.potential_wb_m().to_vec();
    for selected in [
        vec![],
        vec![enclosed[0], enclosed[0]],
        vec![usize::MAX],
        vec![0],
        body,
    ] {
        assert!(
            field
                .maxwell_wrench(&selected, &frame(), [0.0; 3], force_limits(), &mut || false)
                .is_err()
        );
    }
    assert!(
        field
            .maxwell_wrench(
                &enclosed,
                &frame(),
                [f64::NAN, 0.0, 0.0],
                force_limits(),
                &mut || false
            )
            .is_err()
    );
    assert!(
        field
            .maxwell_wrench(
                &enclosed,
                &frame(),
                [0.0; 3],
                MaxwellForceLimits {
                    max_cells: 1,
                    max_edges: 1
                },
                &mut || false
            )
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        field.maxwell_wrench(&enclosed, &frame(), [0.0; 3], force_limits(), &mut || {
            polls += 1;
            polls == 5
        }),
        Err(MagneticError::Cancelled)
    ));
    assert_eq!(field.potential_wb_m(), baseline);

    let (homogeneous, enclosed, _) = specimen(8, AMBIENT_MU, false, DEPTH);
    let bc = boundary(&homogeneous, |[x, _]| 1e-3 * x);
    let field = homogeneous
        .solve(
            &vec![1.0; homogeneous.mesh().faces().len()],
            &bc,
            config(),
            &mut || false,
        )
        .unwrap();
    assert!(
        field
            .maxwell_wrench(&enclosed, &frame(), [0.0; 3], force_limits(), &mut || false)
            .is_err()
    );
}

#[test]
fn a_topological_neighbour_folded_inside_the_contour_is_not_an_exterior_sample() {
    // TriComplex2 admits coherent topological winding and unsigned areas;
    // vertex 3 folds the bottom neighbour over the enclosed triangle.
    let mesh = TriComplex2::from_indexed_triangles(
        "folded-magnetic-contour",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [0.25, 0.25],
            [1.0, 1.0],
            [-1.0, 0.0],
        ],
        vec![[0, 1, 2], [1, 0, 3], [2, 1, 4], [0, 2, 5]],
        Metric2::planar(DEPTH).unwrap(),
    )
    .unwrap();
    let problem =
        PlanarMagnetostatics::new(mesh, vec![material(AMBIENT_MU); 4], limits(), &mut || false)
            .unwrap();
    let bc = boundary(&problem, |[x, _]| 1e-3 * x);
    let field = problem
        .solve(&[0.0; 4], &bc, config(), &mut || false)
        .unwrap();
    assert!(matches!(
        field.maxwell_wrench(&[0], &frame(), [0.0; 3], force_limits(), &mut || false),
        Err(MagneticError::Invalid(
            "Maxwell contour neighbour is folded onto its enclosed side"
        ))
    ));
}

#[test]
fn recovered_induction_history_is_not_silently_treated_as_a_static_field() {
    let (problem, enclosed, _) = specimen(8, AMBIENT_MU, false, DEPTH);
    let bc = boundary(&problem, |[x, _]| 1e-3 * x);
    let field = problem
        .solve(
            &vec![0.0; problem.mesh().faces().len()],
            &bc,
            config(),
            &mut || false,
        )
        .unwrap();
    let induction = PlanarEddy::new(
        &field,
        vec![],
        &(0..problem.mesh().faces().len()).collect::<Vec<_>>(),
        &bc,
        0.0,
        EddyControls {
            linear: config(),
            energy_tolerance_j: 1e-8,
            current_balance_tolerance_a: 1e-8,
        },
        &mut || false,
    )
    .unwrap();
    let recovered = induction.field(&mut || false).unwrap();
    assert!(recovered.impressed_current_density_a_m2().is_none());
    assert!(
        recovered
            .maxwell_wrench(&enclosed, &frame(), [0.0; 3], force_limits(), &mut || false)
            .is_err()
    );
}

#[test]
fn unbiased_unit_current_response_remains_a_valid_static_force_field() {
    let (problem, enclosed, body) = specimen(8, AMBIENT_MU, false, DEPTH);
    let mut profile = vec![0.0; problem.mesh().faces().len()];
    for cell in body {
        let x_sum: f64 = problem.mesh().faces()[cell]
            .iter()
            .map(|&v| problem.mesh().vertices()[v as usize][0])
            .sum();
        profile[cell] = 1000.0 * x_sum.signum();
    }
    let bc = boundary(&problem, |_| 0.0);
    let prescribed: Vec<_> = bc.iter().map(|&(v, _)| v).collect();
    let winding = PlanarWindingBasis::new(
        &problem,
        vec![profile.clone()],
        &prescribed,
        vec![0.0],
        WindingConfig {
            max_windings: 1,
            max_basis_values: 10000,
            pairing_tolerance_h: 1e-8,
            current_balance_tolerance_turns: 1e-8,
            inverse_residual_tolerance: 1e-8,
        },
        config(),
        &mut || false,
    )
    .unwrap();
    let response = &winding.unit_fields()[0];
    assert!(response.is_incremental_response());
    assert!(!problem.has_permanent_magnets());
    assert_eq!(
        response.impressed_current_density_a_m2(),
        Some(profile.as_slice())
    );
    let load = response
        .maxwell_wrench(&enclosed, &frame(), [0.0; 3], force_limits(), &mut || false)
        .unwrap();
    let direct = problem
        .solve(&profile, &bc, config(), &mut || false)
        .unwrap()
        .maxwell_wrench(&enclosed, &frame(), [0.0; 3], force_limits(), &mut || false)
        .unwrap();
    for (got, expected) in load.to_array().into_iter().zip(direct.to_array()) {
        near(got, expected);
    }
    let actual = winding.field_for_currents(&[0.7], &mut || false).unwrap();
    let currents: Vec<_> = profile.iter().map(|j| 0.7 * j).collect();
    assert!(!actual.is_incremental_response());
    assert_eq!(
        actual.impressed_current_density_a_m2(),
        Some(currents.as_slice())
    );
    let reconstructed = actual
        .maxwell_wrench(&enclosed, &frame(), [0.0; 3], force_limits(), &mut || false)
        .unwrap();
    let direct = problem
        .solve(&currents, &bc, config(), &mut || false)
        .unwrap()
        .maxwell_wrench(&enclosed, &frame(), [0.0; 3], force_limits(), &mut || false)
        .unwrap();
    for (got, expected) in reconstructed.to_array().into_iter().zip(direct.to_array()) {
        near(got, expected);
    }
}
