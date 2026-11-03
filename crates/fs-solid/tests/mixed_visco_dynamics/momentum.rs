//! Independent volume moments, signed P2 loads and retained Maxwell histories.

use super::*;
use fs_solid::hyper3::mixed::visco::dynamics::SolidMotionMomentum;

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn near_vector(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    for axis in 0..3 {
        near(actual[axis], expected[axis], tolerance);
    }
}

fn closed(momentum: &SolidMotionMomentum) {
    near_vector(momentum.linear_residual_kg_m_s, [0.0; 3], 3e-10);
    near_vector(momentum.angular_residual_kg_m2_s, [0.0; 3], 3e-10);
}

#[test]
fn affine_and_curved_reference_body_loads_have_original_mass_impulse_and_torque() {
    with_cx(|cx| {
        for curvature in [0.0, 0.2] {
            let (mut nodes, elements) = mesh(cx);
            // Original map x = xi + c*xi^2, y = eta, z = zeta has
            // det(Dx) = 1 + 2*c*xi. Exact simplex polynomial integrals give
            // V = (1+c/2)/6 and the physical reference centroid below. Its
            // corner tetrahedron has a DIFFERENT volume whenever c != 0.
            for node in &mut nodes {
                node[0] += curvature * node[0] * node[0];
            }
            let volume = VOLUME * (1.0 + curvature / 2.0);
            let mass = DENSITY * volume;
            let centroid = [
                (0.25 + 0.3 * curvature + 0.1 * curvature * curvature) / (1.0 + curvature / 2.0),
                (0.25 + 0.1 * curvature) / (1.0 + curvature / 2.0),
                (0.25 + 0.1 * curvature) / (1.0 + curvature / 2.0),
            ];
            for treatment in [
                VolumetricTreatment::OriginalMaterial,
                VolumetricTreatment::ExactIncompressible,
            ] {
                let material = law(10.0 * MU, true);
                let before_velocity = [0.2, -0.1, 0.05];
                let acceleration = [0.3, -0.2, 0.1];
                let dt = 0.02;
                let origin = [-0.4, 0.15, -0.2];
                let motion = MixedViscoTetMotion::try_new(
                    solid(&nodes, &elements, &material, treatment, cx),
                    &[DENSITY],
                    &vec![before_velocity; nodes.len()],
                    MASS_BYTES,
                    cx,
                )
                .unwrap();
                assert_eq!(motion.volumetric_treatment(), treatment);
                near(motion.reference_mass_kg(), mass, 3e-13);
                let (step, momentum) = motion
                    .advance_with_momentum(
                        dt,
                        &[],
                        &[],
                        acceleration.map(|a| DENSITY * a),
                        None,
                        settings(),
                        controls(),
                        origin,
                        cx,
                    )
                    .unwrap();
                let after_velocity =
                    std::array::from_fn(|axis| before_velocity[axis] + dt * acceleration[axis]);
                let before_arm = std::array::from_fn(|axis| centroid[axis] - origin[axis]);
                let after_arm =
                    std::array::from_fn(|axis| before_arm[axis] + dt * after_velocity[axis]);
                let impulse = acceleration.map(|a| mass * dt * a);
                let before_angular = cross(before_arm, before_velocity.map(|v| mass * v));
                let after_angular = cross(after_arm, after_velocity.map(|v| mass * v));
                let angular_change =
                    std::array::from_fn(|axis| after_angular[axis] - before_angular[axis]);
                let angular_defect = cross(after_velocity, before_velocity).map(|v| mass * dt * v);
                assert!(angular_defect.iter().any(|v| v.abs() > 1e-6));
                assert_eq!(momentum.change.origin_m, origin);
                assert_eq!(momentum.applied.origin_m, origin);
                assert_eq!(momentum.support.origin_m, origin);
                near_vector(momentum.applied.linear_kg_m_s, impulse, 3e-13);
                near_vector(
                    momentum.applied.angular_kg_m2_s,
                    cross(after_arm, impulse),
                    3e-12,
                );
                near_vector(momentum.change.linear_kg_m_s, impulse, 3e-10);
                near_vector(momentum.change.angular_kg_m2_s, angular_change, 3e-10);
                near_vector(momentum.support.linear_kg_m_s, [0.0; 3], 0.0);
                near_vector(momentum.support.angular_kg_m2_s, [0.0; 3], 0.0);
                near_vector(
                    momentum.angular_integration_defect_kg_m2_s,
                    angular_defect,
                    3e-11,
                );
                for velocity in step.motion.velocity_mps() {
                    near_vector(*velocity, after_velocity, 3e-9);
                }
                closed(&momentum);
                work_identity(&step);
            }
        }
    });
}

#[test]
fn curved_deformation_retains_signed_support_loads_and_original_relaxation_history() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(10.0 * MU, true);
        let shear = 0.08;
        let fixed: Vec<_> = nodes
            .iter()
            .enumerate()
            .flat_map(|(node, x)| {
                [shear * x[1] * x[2], 0.0, 0.0]
                    .into_iter()
                    .enumerate()
                    .map(move |(axis, u)| (3 * node + axis, u))
            })
            .collect();
        let mut initial = solid(
            &nodes,
            &elements,
            &material,
            VolumetricTreatment::OriginalMaterial,
            cx,
        );
        for dt in [0.0, 0.007] {
            initial = initial
                .advance(
                    dt,
                    &fixed,
                    &[],
                    [0.0; 3],
                    None,
                    settings(),
                    controls().material,
                    cx,
                )
                .unwrap()
                .solid;
        }
        assert!(
            initial
                .quadrature_states()
                .iter()
                .any(|s| s.dissipated > 0.0)
        );
        let before_velocity = [0.04, 0.02, -0.03];
        let motion = MixedViscoTetMotion::try_new(
            initial,
            &[DENSITY],
            &vec![before_velocity; nodes.len()],
            MASS_BYTES,
            cx,
        )
        .unwrap();
        let before = (
            motion.solid().displacement_m().to_vec(),
            motion.solid().pressure_pa().to_vec(),
            motion.solid().quadrature_states().to_vec(),
            motion.velocity_mps().to_vec(),
            motion.solid().time_s(),
        );
        let node = elements[0][4 + TET10_EDGES.iter().position(|e| *e == [2, 3]).unwrap()];
        let point_force = [0.3, -0.2, 0.4];
        let nodal = [
            (3 * node, 0.1),
            (3 * node, 0.2),
            (3 * node + 1, point_force[1]),
            (3 * node + 2, point_force[2]),
        ];
        let body_force = [1.2, -0.9, 0.6];
        let origin = [-0.15, 0.1, 0.4];
        let dt = 0.02;
        let baseline = motion
            .advance(dt, &fixed, &[], [0.0; 3], None, settings(), controls(), cx)
            .unwrap();
        let (step, momentum) = motion
            .advance_with_momentum(
                dt,
                &fixed,
                &nodal,
                body_force,
                None,
                settings(),
                controls(),
                origin,
                cx,
            )
            .unwrap();
        // u_x = c*y*z has det(F)=1 and moves the original centroid by c/20.
        // Corner-only geometry misses this shift altogether.
        let centroid = [0.25 + shear / 20.0, 0.25, 0.25];
        let arm = std::array::from_fn(|axis| centroid[axis] - origin[axis]);
        let point = [shear / 4.0, 0.5, 0.5];
        let point_arm = std::array::from_fn(|axis| point[axis] - origin[axis]);
        let body_impulse = body_force.map(|force| dt * VOLUME * force);
        let point_impulse = point_force.map(|force| dt * force);
        let applied = std::array::from_fn(|axis| body_impulse[axis] + point_impulse[axis]);
        let body_torque = cross(arm, body_impulse);
        let point_torque = cross(point_arm, point_impulse);
        let applied_torque = std::array::from_fn(|axis| body_torque[axis] + point_torque[axis]);
        let change = before_velocity.map(|v| -DENSITY * VOLUME * v);
        let angular_change = cross(arm, change);
        near_vector(momentum.applied.linear_kg_m_s, applied, 3e-13);
        near_vector(momentum.applied.angular_kg_m2_s, applied_torque, 3e-13);
        near_vector(momentum.change.linear_kg_m_s, change, 3e-12);
        near_vector(momentum.change.angular_kg_m2_s, angular_change, 3e-12);
        near_vector(
            momentum.support.linear_kg_m_s,
            std::array::from_fn(|axis| change[axis] - applied[axis]),
            3e-11,
        );
        near_vector(
            momentum.support.angular_kg_m2_s,
            std::array::from_fn(|axis| angular_change[axis] - applied_torque[axis]),
            3e-11,
        );
        near_vector(momentum.angular_integration_defect_kg_m2_s, [0.0; 3], 0.0);
        for (&(dof, reaction), &(baseline_dof, baseline_reaction)) in step
            .equilibrium
            .reactions_n
            .iter()
            .zip(&baseline.equilibrium.reactions_n)
        {
            assert_eq!(dof, baseline_dof);
            let signed_volume = if elements[0][..4].contains(&(dof / 3)) {
                -VOLUME / 20.0
            } else {
                VOLUME / 5.0
            };
            let force = signed_volume * body_force[dof % 3]
                + if dof / 3 == node {
                    point_force[dof % 3]
                } else {
                    0.0
                };
            near(reaction - baseline_reaction, -force, 3e-12);
        }
        assert_eq!(step.equilibrium.reactions_n.len(), 3 * nodes.len());
        assert_eq!(
            step.motion.solid().quadrature_states(),
            baseline.motion.solid().quadrature_states()
        );
        assert_eq!(
            step.motion.solid().first_piola_pa(),
            baseline.motion.solid().first_piola_pa()
        );
        assert_eq!(
            step.motion.solid().pressure_pa(),
            baseline.motion.solid().pressure_pa()
        );
        assert_eq!(
            step.motion.solid().displacement_m(),
            baseline.motion.solid().displacement_m()
        );
        assert_eq!(step.motion.velocity_mps(), baseline.motion.velocity_mps());
        assert!(
            step.motion
                .velocity_mps()
                .iter()
                .flatten()
                .all(|v| *v == 0.0)
        );
        assert_eq!(
            step.material_energy.dissipated_j,
            baseline.material_energy.dissipated_j
        );
        assert!(step.material_energy.dissipated_j > 0.0);
        near(step.motion.solid().time_s(), before.4 + dt, 0.0);
        assert_eq!(
            before,
            (
                motion.solid().displacement_m().to_vec(),
                motion.solid().pressure_pa().to_vec(),
                motion.solid().quadrature_states().to_vec(),
                motion.velocity_mps().to_vec(),
                motion.solid().time_s(),
            )
        );
        closed(&momentum);
        work_identity(&step);
        work_identity(&baseline);
    });
}
