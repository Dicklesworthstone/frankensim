//! G1/G3: source thermal motion retains original mass and accepts complete events.
use super::*;
use fs_couple::finite_solid::thermal::visco::plane::dynamics::ExpansionViscoPlaneMotion;
use fs_solid::hyper2d::visco::dynamics::FiniteViscoMotionControls;

#[test]
fn thermal_motion_keeps_original_inventory_velocity_and_memory_through_atomic_events() {
    with_cx(|_, cx| {
        let source = card(STRETCHES, "T", false);
        let datum = endpoint(&source, Some(shift()));
        let field = ExpansionViscoPlaneField::per_element(vec![datum.clone(), datum], cx).unwrap();
        let target = field
            .at_temperatures(&[&source, &source], &[300.1, 300.2], cx)
            .unwrap();
        let wrong_calibration = ExpansionViscoPlaneField::uniform(endpoint(&source, None));
        let mesh = Mesh2::quads(1.0, 0.7, 2, 1);
        let initial: Vec<_> = mesh.nodes.iter().map(|[x, y]| held(*x, *y)).collect();
        let boundary = walls(&held);
        let velocity = vec![[0.2, -0.1]; mesh.nodes.len()];
        let acceleration = vec![[0.7, -0.3]; mesh.nodes.len()];
        // The existing authored fixture's original density is 6 kg/m³.
        // This explicit inventory is never queried again at a new temperature.
        let density = [6.0, 6.0];
        let controls = FiniteViscoMotionControls {
            material: allowances(),
            maximum_numerical_kinetic_loss_j: 0.1,
        };
        for mode in [
            FiniteViscoPlaneMode::FiniteBulk,
            FiniteViscoPlaneMode::Incompressible,
        ] {
            let placed = ExpansionViscoPlaneSolid::try_new_from(
                &mesh,
                &field,
                mode,
                0.12,
                &initial,
                PlaneStressSettings::default(),
                cx,
            )
            .unwrap();
            let accepted = placed
                .advance(0.08, &boundary, &[], settings(), allowances(), cx)
                .unwrap()
                .solid;
            let motion =
                ExpansionViscoPlaneMotion::try_new(accepted, &density, &velocity, 1 << 20, cx)
                    .unwrap();
            let before_full = motion.native().solid().point_states().map(|s| s.to_vec());
            let before_remainder = motion
                .native()
                .solid()
                .point_remainder_states()
                .map(|s| s.to_vec());
            let before_displacement = motion.native().solid().displacement_m().to_vec();
            near(
                motion.native().reference_mass_kg(),
                1.0 * 0.7 * 0.12 * 6.0,
                1e-13,
            );
            let load = motion
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap();
            for axis in 0..2 {
                let total: f64 = load
                    .iter()
                    .filter(|(dof, _)| dof % 2 == axis)
                    .map(|(_, f)| f)
                    .sum();
                near(
                    total,
                    motion.native().reference_mass_kg() * acceleration[0][axis],
                    1e-13,
                );
            }
            let changed = motion.change_reference(&target, allowances(), cx).unwrap();
            assert_eq!(changed.motion.native().velocity_mps(), velocity);
            assert_eq!(changed.motion.native().reference_density_kg_m3(), density);
            assert_eq!(
                changed.motion.native().reference_mass_kg(),
                motion.native().reference_mass_kg()
            );
            assert_eq!(
                changed
                    .motion
                    .reference_acceleration_load_n(&acceleration, cx)
                    .unwrap(),
                load
            );
            assert_eq!(
                changed.motion.native().solid().displacement_m(),
                before_displacement
            );
            assert_eq!(
                changed.motion.native().solid().time_s(),
                motion.native().solid().time_s()
            );
            if let Some(old) = &before_full {
                for (before, after) in old
                    .iter()
                    .zip(changed.motion.native().solid().point_states().unwrap())
                {
                    assert_eq!(
                        before.inverse_viscous_cauchy_green,
                        after.inverse_viscous_cauchy_green
                    );
                    assert_eq!(before.dissipated, after.dissipated);
                }
            }
            if let Some(old) = &before_remainder {
                for (before, after) in old.iter().zip(
                    changed
                        .motion
                        .native()
                        .solid()
                        .point_remainder_states()
                        .unwrap(),
                ) {
                    assert_eq!(
                        before.inverse_viscous_cauchy_green,
                        after.inverse_viscous_cauchy_green
                    );
                    assert_eq!(before.dissipated, after.dissipated);
                }
            }
            let sequential = changed
                .motion
                .advance(0.02, &boundary, &[], &load, settings(), controls, cx)
                .unwrap();
            let combined = motion
                .change_reference_and_advance(
                    &target,
                    0.02,
                    &boundary,
                    &[],
                    &load,
                    settings(),
                    allowances(),
                    controls,
                    cx,
                )
                .unwrap();
            assert_eq!(combined.reference, changed.energy);
            assert_eq!(combined.material_energy, sequential.material_energy);
            assert_eq!(combined.inertia, sequential.inertia);
            assert_eq!(combined.reactions_n, sequential.reactions_n);
            assert_eq!(
                combined.motion.native().velocity_mps(),
                sequential.motion.native().velocity_mps()
            );
            assert_eq!(
                combined.motion.native().solid().point_states(),
                sequential.motion.native().solid().point_states()
            );
            assert_eq!(
                combined.motion.native().solid().point_remainder_states(),
                sequential.motion.native().solid().point_remainder_states()
            );
            assert!(std::ptr::eq(combined.motion.source_field(), &target));
            assert!(combined.material_energy.dissipated_j > 0.0);
            assert!(combined.inertia.numerical_loss_j > 0.0);
            near(combined.inertia.kinetic_after_j, 0.0, 1e-13);
            near(combined.balance_error_j, 0.0, 1e-8);
            near(
                combined.reference_input_j(),
                combined.reference.reference_work_j
                    + combined.reference.constraint_reference_work_j,
                0.0,
            );
            if mode == FiniteViscoPlaneMode::Incompressible {
                assert!(combined.reference.constraint_reference_work_j.abs() > 1e-9);
            } else {
                near(combined.reference.constraint_reference_work_j, 0.0, 0.0);
            }
            assert!(matches!(
                motion.change_reference(&wrong_calibration, allowances(), cx),
                Err(ExpansionViscoPlaneError::Source { element: 0, .. })
            ));
            // The fully held endpoint removes the admitted initial velocity.
            // Its positive BE kinetic loss is a separate late acceptance gate,
            // so this refusal occurs after a valid source event and native step.
            let strict = FiniteViscoMotionControls {
                maximum_numerical_kinetic_loss_j: 0.0,
                ..controls
            };
            assert!(matches!(
                motion.change_reference_and_advance(
                    &target,
                    0.02,
                    &boundary,
                    &[],
                    &load,
                    settings(),
                    allowances(),
                    strict,
                    cx,
                ),
                Err(ExpansionViscoPlaneError::Native(
                    FiniteViscoPlaneError::EnergyGate {
                        quantity: "backward-Euler numerical kinetic loss",
                        ..
                    }
                ))
            ));
            assert!(std::ptr::eq(motion.source_field(), &field));
            assert_eq!(motion.native().velocity_mps(), velocity);
            assert_eq!(motion.native().reference_density_kg_m3(), density);
            assert_eq!(
                motion.native().solid().displacement_m(),
                before_displacement
            );
            assert_eq!(
                motion.native().solid().point_states(),
                before_full.as_deref()
            );
            assert_eq!(
                motion.native().solid().point_remainder_states(),
                before_remainder.as_deref()
            );
            near(motion.native().solid().time_s(), 0.08, 1e-15);
            near(combined.motion.native().solid().time_s(), 0.1, 1e-15);
        }
    });
}
