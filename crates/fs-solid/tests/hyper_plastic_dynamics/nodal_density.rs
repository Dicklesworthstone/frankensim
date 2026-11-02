//! Graded reference density through the original finite-plastic Newton step.
use super::*;

const NODAL_DENSITY: [f64; 4] = [3.0, 15.0, 3.0, 3.0];

fn nodal_motion<'a>(
    material: &'a FiniteJ2,
    density: &[f64],
    velocity: &[[f64; 3]],
    cx: &Cx<'_>,
) -> Result<FinitePlasticTetMotion<'a>, fs_solid::hyper3::plastic::FinitePlasticTetError> {
    let solid = FinitePlasticTetSolid::try_new(
        &NODES,
        &TETS,
        FinitePlasticLawField::Uniform(material),
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap();
    FinitePlasticTetMotion::try_new_with_nodal_density(solid, density, velocity, cx)
}

#[test]
fn graded_gravity_accelerates_the_native_free_body_without_averaging_its_inventory() {
    with_cx(|cx| {
        let material = law();
        let velocity0 = [0.2, -0.1, 0.05];
        let acceleration = [1.2, -9.81, 0.4];
        let dt = 0.02;
        let mut motion = nodal_motion(&material, &NODAL_DENSITY, &[velocity0; 4], cx).unwrap();
        close(motion.reference_mass_kg(), 1.0, 1e-14);
        assert_eq!(motion.reference_density_kg_m3(), [6.0]);
        for interval in 1..=3 {
            let loads = motion
                .reference_acceleration_load_n(&[acceleration; 4], cx)
                .unwrap();
            let mut actual_load = [0.0; 12];
            for &(dof, force) in &loads {
                actual_load[dof] += force;
            }
            // Independent integral of rho*N_i for rho=3+12*x. The total mass
            // is one, but each point does not carry one quarter of that mass.
            for (node, share) in [0.225, 0.325, 0.225, 0.225].into_iter().enumerate() {
                for i in 0..3 {
                    close(actual_load[3 * node + i], share * acceleration[i], 2e-14);
                }
            }
            let count = f64::from(interval);
            let before: [f64; 3] =
                core::array::from_fn(|i| velocity0[i] + (count - 1.0) * dt * acceleration[i]);
            let after: [f64; 3] =
                core::array::from_fn(|i| velocity0[i] + count * dt * acceleration[i]);
            let displacement: [f64; 3] = core::array::from_fn(|i| {
                count * dt * velocity0[i] + 0.5 * count * (count + 1.0) * dt * dt * acceleration[i]
            });
            let step = motion
                .advance(&[], &loads, [0.0; 3], dt, settings(), controls(), cx)
                .unwrap();
            for (u, v) in step
                .motion
                .solid()
                .displacement_m()
                .iter()
                .zip(step.motion.velocity_mps())
            {
                for i in 0..3 {
                    close(u[i], displacement[i], 1e-10);
                    close(v[i], after[i], 1e-8);
                }
            }
            assert_eq!(
                step.motion.reference_nodal_density_kg_m3(),
                Some(NODAL_DENSITY.as_slice())
            );
            assert_eq!(step.motion.reference_density_kg_m3(), [6.0]);
            assert_eq!(step.motion.reference_mass_kg(), motion.reference_mass_kg());
            assert_eq!(
                step.motion.solid().element_states(),
                &[material.initial_state()]
            );
            close(step.material_energy.stored_after_j, 0.0, 1e-12);
            close(
                step.inertia.kinetic_before_j,
                0.5 * norm_squared(before),
                1e-9,
            );
            close(
                step.inertia.kinetic_after_j,
                0.5 * norm_squared(after),
                1e-9,
            );
            close(
                step.inertia.numerical_loss_j,
                0.5 * dt * dt * norm_squared(acceleration),
                1e-9,
            );
            let work = dt * (0..3).map(|i| acceleration[i] * after[i]).sum::<f64>();
            close(step.inertia.inertial_work_j, work, 1e-9);
            close(step.material_energy.dead_load_endpoint_work_j, work, 1e-9);
            close(step.balance_error_j, 0.0, 1e-9);
            motion = step.motion;
        }
    });
}

#[test]
fn uniform_nodal_inventory_retains_original_native_motion_and_refuses_invalid_fields() {
    with_cx(|cx| {
        let material = law();
        let velocity = [[0.2, -0.1, 0.05]; 4];
        let original = initial_motion(&material, &velocity, cx);
        let nodal = nodal_motion(&material, &[DENSITY[0]; 4], &velocity, cx).unwrap();
        assert_eq!(original.reference_nodal_density_kg_m3(), None);
        assert_eq!(
            original.reference_mass_kg().to_bits(),
            nodal.reference_mass_kg().to_bits()
        );
        let acceleration = [[1.0, -0.5, 0.25]; 4];
        let loads = original
            .reference_acceleration_load_n(&acceleration, cx)
            .unwrap();
        assert_eq!(
            loads,
            nodal
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap()
        );
        let old = original
            .advance(&[], &loads, [0.0; 3], 0.02, settings(), controls(), cx)
            .unwrap();
        let new = nodal
            .advance(&[], &loads, [0.0; 3], 0.02, settings(), controls(), cx)
            .unwrap();
        assert_eq!(
            old.motion.solid().displacement_m(),
            new.motion.solid().displacement_m()
        );
        assert_eq!(old.motion.velocity_mps(), new.motion.velocity_mps());
        assert_eq!(old.inertia, new.inertia);
        for density in [
            &[1.0, 2.0, 3.0][..],
            &[1.0, 0.0, 2.0, 3.0],
            &[1.0, f64::NAN, 2.0, 3.0],
        ] {
            assert!(nodal_motion(&material, density, &velocity, cx).is_err());
        }
        assert!(
            new.motion
                .reference_acceleration_load_n(&[[f64::NAN; 3]; 4], cx)
                .is_err()
        );
        assert!(new.motion.reference_acceleration_load_n(&[], cx).is_err());
        assert_eq!(
            new.motion.reference_nodal_density_kg_m3(),
            Some([6.0; 4].as_slice())
        );
    });
}
