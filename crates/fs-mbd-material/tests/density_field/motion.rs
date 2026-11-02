//! G1/G3/G4: exact material momenta enter the original native motion owner.
use super::*;
use fs_ga::So3Tangent;
use fs_mbd::{RigidBodyIntegrator, Wrench};
use fs_mbd_material::density_field::motion::{RigidMotionProjection, RigidProjectionControls};
use fs_mbd_material::principal::PrincipalFrameError;

fn controls(maximum_unresolved_kinetic_energy_j: f64) -> RigidProjectionControls {
    RigidProjectionControls {
        relative_tolerance: 1e-11,
        maximum_unresolved_kinetic_energy_j,
    }
}
fn placement() -> Se3 {
    Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.3, -0.5, 0.7))).unwrap(),
        Vec3::new(4.0, -2.0, 3.0),
    )
    .unwrap()
}
fn velocity(x: Vec3) -> Vec3 {
    Vec3::new(
        2.0 + 0.7 * x.x - 0.2 * x.y,
        -1.0 + 0.3 * x.x + 0.4 * x.z,
        0.5 - 0.6 * x.y + 0.2 * x.z,
    )
}
// Independent analytic barycentric integrals, without production quadrature.
// Integral rho N_a = V*(sumrho+rho_a)/20; integral rho N_a N_b is
// V*(sumrho+rho_a+rho_b)/120 off diagonal and twice that on diagonal.
fn moments(x: [Vec3; 4], rho: [f64; 4], u: [Vec3; 4], reference_volume: f64) -> (Vec3, Vec3, f64) {
    let sum = rho.into_iter().sum::<f64>();
    let mut p = Vec3::new(0.0, 0.0, 0.0);
    let mut l_origin = Vec3::new(0.0, 0.0, 0.0);
    let mut kinetic = 0.0;
    for a in 0..4 {
        p = p + u[a].scale(reference_volume * (sum + rho[a]) / 20.0);
        for b in 0..4 {
            let mass =
                reference_volume * (sum + rho[a] + rho[b]) / if a == b { 60.0 } else { 120.0 };
            l_origin = l_origin + x[a].cross(u[b]).scale(mass);
            kinetic += 0.5 * mass * u[a].dot(u[b]);
        }
    }
    (p, l_origin, kinetic)
}

fn check_momenta(projection: &RigidMotionProjection<'_>, pose: Se3, p: Vec3, l: Vec3) {
    let expected_p = pose.rotation().rotate(p).unwrap();
    vec_close(ga(projection.state.linear_momentum_world()), expected_p);
    let expected_l = pose.rotation().rotate(l).unwrap() + pose.translation().cross(expected_p);
    let actual = RigidBodyIntegrator::new(Gravity::ZERO)
        .diagnostics(projection.state, projection.frame.properties())
        .unwrap();
    vec_close(ga(actual.angular_momentum_world), expected_l);
    let omega_principal = projection
        .frame
        .properties()
        .angular_velocity_body_checked(projection.state.angular_momentum_body())
        .unwrap();
    let omega = projection
        .frame
        .body_from_principal()
        .rotation()
        .rotate(ga(omega_principal))
        .unwrap();
    vec_close(
        projection.frame.inertia().inertia_com().apply(omega),
        l - projection.frame.inertia().center_of_mass().cross(p),
    );
}

#[test]
fn graded_nodal_motion_preserves_cubic_momenta_and_drives_native_mechanics() {
    let rho = [1000.0, 5000.0, 1600.0, 2700.0];
    let x = vertices();
    let u = x.map(velocity);
    let body = resolve(&x, &rho, &[[0, 1, 2, 3]]);
    let pose = placement();
    let result = body
        .project_rigid_motion(pose, &u, controls(1000.0), || false)
        .unwrap();
    assert!(core::ptr::eq(result.reference(), &body));
    let (p, l, kinetic) = moments(x, rho, u, 1.0 / 6.0);
    check_momenta(&result, pose, p, l);
    close(result.nodal_kinetic_energy_j, kinetic);
    close(
        result.rigid_kinetic_energy_j + result.unresolved_kinetic_energy_j,
        kinetic,
    );
    assert!(result.unresolved_kinetic_energy_j > 0.1);
    assert!(result.relative_energy_residual < 1e-12);

    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let step = integrator
        .step(result.state, result.frame.properties(), Wrench::ZERO, 1e-4)
        .unwrap();
    vec_close(
        ga(step.diagnostics_after.linear_momentum_world),
        ga(step.diagnostics_before.linear_momentum_world),
    );
    vec_close(
        ga(step.diagnostics_after.angular_momentum_world),
        ga(step.diagnostics_before.angular_momentum_world),
    );
    close(
        step.diagnostics_after.mechanical_energy,
        result.rigid_kinetic_energy_j,
    );
    let expected_position = ga(result.state.pose().position_world())
        + pose
            .rotation()
            .rotate(p)
            .unwrap()
            .scale(1e-4 / body.inertia().mass());
    vec_close(
        ga(step.state_after.pose().position_world()),
        expected_position,
    );

    let impulse_world = fs_mbd::Vec3::new(0.7, -0.4, 1.2);
    let impulse = result
        .state
        .apply_impulse_at_body_point(
            result.frame.properties(),
            result.frame.arm_for_body_point(x[1]).unwrap(),
            impulse_world,
        )
        .unwrap();
    let after = integrator
        .diagnostics(impulse.state_after, result.frame.properties())
        .unwrap();
    vec_close(
        ga(after.linear_momentum_world),
        ga(result.state.linear_momentum_world().add(impulse_world)),
    );
    let point_world = pose.rotation().rotate(x[1]).unwrap() + pose.translation();
    vec_close(
        ga(after.angular_momentum_world),
        ga(step.diagnostics_before.angular_momentum_world) + point_world.cross(ga(impulse_world)),
    );

    // Equal total mass with uniform density must not silently replace the field.
    let uniform = resolve(
        &x,
        &[rho.into_iter().sum::<f64>() / 4.0; 4],
        &[[0, 1, 2, 3]],
    );
    let flattened = uniform
        .project_rigid_motion(pose, &u, controls(1000.0), || false)
        .unwrap();
    close(uniform.inertia().mass(), body.inertia().mass());
    assert!(
        (ga(flattened.state.linear_momentum_world()) - ga(result.state.linear_momentum_world()))
            .norm()
            > 1.0
    );
}

#[test]
fn deformed_projection_preserves_reference_inventory_and_partition_equivalence() {
    let rho = [1000.0, 5000.0, 1000.0, 1000.0];
    let body = resolve(&vertices(), &rho, &[[0, 1, 2, 3]]);
    let transform = Mat3 {
        m: [2.0, 0.4, 0.0, 0.0, 3.0, -0.2, 0.0, 0.0, 0.5],
    };
    let current = vertices().map(|v| transform.apply(v) + Vec3::new(10.0, -3.0, 4.0));
    let deformed = body.deform(&current, || false).unwrap();
    let u = current.map(velocity);
    let result = deformed
        .project_rigid_motion(placement(), &u, controls(10000.0), || false)
        .unwrap();
    let (p, l, kinetic) = moments(current, rho, u, 1.0 / 6.0);
    check_momenta(&result, placement(), p, l);
    close(result.nodal_kinetic_energy_j, kinetic);
    assert_eq!(
        result.frame.properties().mass().to_bits(),
        body.inertia().mass().to_bits()
    );
    assert!(core::ptr::eq(result.reference(), &body));

    let mut refined_reference = vertices().to_vec();
    refined_reference.push(Vec3::new(0.25, 0.25, 0.25));
    let refined = resolve(
        &refined_reference,
        &[1000.0, 5000.0, 1000.0, 1000.0, 2000.0],
        &[[4, 1, 2, 3], [0, 4, 2, 3], [0, 1, 4, 3], [0, 1, 2, 4]],
    );
    let positions: Vec<_> = refined_reference
        .iter()
        .map(|&v| transform.apply(v) + Vec3::new(10.0, -3.0, 4.0))
        .collect();
    let velocities: Vec<_> = positions.iter().copied().map(velocity).collect();
    let refined_deformed = refined.deform(&positions, || false).unwrap();
    let refined_result = refined_deformed
        .project_rigid_motion(placement(), &velocities, controls(10000.0), || false)
        .unwrap();
    check_momenta(&refined_result, placement(), p, l);
    inertia_close(result.frame.inertia(), refined_result.frame.inertia());
    close(
        result.nodal_kinetic_energy_j,
        refined_result.nodal_kinetic_energy_j,
    );
    close(
        result.unresolved_kinetic_energy_j,
        refined_result.unresolved_kinetic_energy_j,
    );
}

#[test]
fn unequal_cell_jacobians_keep_each_original_mass_in_the_motion_integral() {
    let mut reference = vertices().to_vec();
    reference.push(Vec3::new(0.0, 0.0, -1.0));
    let rho = [1000.0, 5000.0, 1600.0, 2700.0, 7000.0];
    let body = resolve(&reference, &rho, &[[0, 1, 2, 3], [0, 1, 2, 4]]);
    let mut current = reference.clone();
    current[3].z = 2.0;
    current[4].z = -3.0;
    let deformed = body.deform(&current, || false).unwrap();
    let velocities: Vec<_> = current.iter().copied().map(velocity).collect();
    let result = deformed
        .project_rigid_motion(placement(), &velocities, controls(10000.0), || false)
        .unwrap();
    let mut p = Vec3::new(0.0, 0.0, 0.0);
    let mut l = Vec3::new(0.0, 0.0, 0.0);
    let mut kinetic = 0.0;
    for tet in body.tetrahedra() {
        let (cell_p, cell_l, cell_k) = moments(
            tet.map(|v| current[v]),
            tet.map(|v| rho[v]),
            tet.map(|v| velocities[v]),
            1.0 / 6.0,
        );
        p = p + cell_p;
        l = l + cell_l;
        kinetic += cell_k;
    }
    check_momenta(&result, placement(), p, l);
    close(result.nodal_kinetic_energy_j, kinetic);
    assert_eq!(
        result.frame.properties().mass().to_bits(),
        body.inertia().mass().to_bits()
    );
}

#[test]
fn rigid_motion_is_retained_and_common_boost_cannot_hide_internal_energy() {
    let x = vertices();
    let body = resolve(&x, &[1000.0, 5000.0, 1600.0, 2700.0], &[[0, 1, 2, 3]]);
    let omega = Vec3::new(0.7, -0.4, 1.1);
    let origin_velocity = Vec3::new(2.0, -1.0, 0.5);
    let rigid = x.map(|x| origin_velocity + omega.cross(x));
    let result = body
        .project_rigid_motion(placement(), &rigid, controls(1e-20), || false)
        .unwrap();
    let expected = result
        .frame
        .state_from_kinematics(
            placement(),
            Twist {
                angular: omega,
                linear: origin_velocity,
            },
        )
        .unwrap();
    vec_close(
        ga(result.state.linear_momentum_world()),
        ga(expected.linear_momentum_world()),
    );
    vec_close(
        ga(result.state.angular_momentum_body()),
        ga(expected.angular_momentum_body()),
    );
    close(
        result.nodal_kinetic_energy_j,
        expected.kinetic_energy(result.frame.properties()).unwrap(),
    );
    assert!(result.unresolved_kinetic_energy_j < 1e-20);

    let internal = x.map(|x| Vec3::new(x.x, -2.0 * x.y, 3.0 * x.z));
    let boost = Vec3::new(2.0_f64.powi(40), -2.0_f64.powi(39), 2.0_f64.powi(38));
    let boosted = internal.map(|v| v + boost);
    let unboosted = body
        .project_rigid_motion(placement(), &internal, controls(10000.0), || false)
        .unwrap();
    let boosted_result = body
        .project_rigid_motion(placement(), &boosted, controls(10000.0), || false)
        .unwrap();
    assert_eq!(
        boosted_result.unresolved_kinetic_energy_j.to_bits(),
        unboosted.unresolved_kinetic_energy_j.to_bits()
    );
    assert!(unboosted.unresolved_kinetic_energy_j > 1.0);
    assert!(
        body.project_rigid_motion(
            placement(),
            &boosted,
            controls(0.5 * unboosted.unresolved_kinetic_energy_j),
            || false
        )
        .is_err()
    );
    let translation = body
        .project_rigid_motion(placement(), &[boost; 4], controls(0.0), || false)
        .unwrap();
    assert_eq!(translation.unresolved_kinetic_energy_j, 0.0);
}

#[test]
fn refusal_and_every_cancellation_boundary_publish_no_projection() {
    let body = resolve(
        &vertices(),
        &[1000.0, 5000.0, 1600.0, 2700.0],
        &[[0, 1, 2, 3]],
    );
    let before = body.clone();
    let velocities = vertices().map(velocity);
    assert!(
        body.project_rigid_motion(placement(), &velocities[..3], controls(1000.0), || false)
            .is_err()
    );
    let mut nonfinite = velocities;
    nonfinite[2].y = f64::NAN;
    assert!(
        body.project_rigid_motion(placement(), &nonfinite, controls(1000.0), || false)
            .is_err()
    );
    assert!(
        body.project_rigid_motion(placement(), &velocities, controls(0.0), || false)
            .is_err()
    );
    assert!(
        body.project_rigid_motion(placement(), &velocities, controls(f64::INFINITY), || false)
            .is_err()
    );
    let mut polls = 0;
    body.project_rigid_motion(placement(), &velocities, controls(1000.0), || {
        polls += 1;
        false
    })
    .unwrap();
    for checkpoint in 1..=polls {
        let mut call = 0;
        assert_eq!(
            body.project_rigid_motion(placement(), &velocities, controls(1000.0), || {
                call += 1;
                call == checkpoint
            }),
            Err(PrincipalFrameError::Cancelled)
        );
    }
    assert_eq!(before, body);
}
