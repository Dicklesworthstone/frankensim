//! Retained P2 inertia and pressure admission during physical reference changes.

use super::*;
use fs_material::hyper::reference::ReferenceDeformation;
use fs_solid::hyper3::mixed::visco::MixedViscoReferenceChangeControls;
use fs_solid::hyper3::mixed::visco::dynamics::{
    FiniteViscoContactDriftControls, SolidImpulseControls,
};
use fs_solver::CheckedCgConfig;

#[test]
fn reference_change_retains_graded_quadratic_motion_pressure_and_split_binding() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let bulk = 10.0 * MU;
        let old = law(bulk, true);
        let reference =
            ReferenceDeformation::try_new([1.1, 0.04, 0.0, 0.0, 0.9, 0.03, 0.0, 0.0, 1.05])
                .unwrap();
        let target = old.clone().with_reference_deformation(reference);
        let u: Vec<_> = nodes.iter().map(|x| [0.08 * x[1], 0.0, 0.0]).collect();
        let solid = MixedViscoTetSolid::try_new_from(
            &nodes,
            &elements,
            MixedViscoLawField::Uniform(&old),
            VolumetricTreatment::OriginalMaterial,
            &u,
            &[0.0; 4],
            1e-11,
            TetAssemblyBudget::standard(),
            4,
            64,
            cx,
        )
        .unwrap();
        let density = [[3.0, 15.0, 3.0, 3.0]];
        let velocity: Vec<_> = nodes
            .iter()
            .map(|x| [0.2 + 0.3 * x[0], -0.1 + 0.2 * x[2], 0.05 - 0.2 * x[1]])
            .collect();
        let original = MixedViscoTetMotion::try_new_with_reference_density(
            solid, &density, &velocity, MASS_BYTES, cx,
        )
        .unwrap();
        let impulse = SolidImpulseControls {
            solve: CheckedCgConfig {
                tolerance: 1e-13,
                max_iterations: 64,
                max_corrections: 2,
            },
            maximum_workspace_components: 4096,
            linear_tolerance_kg_m_s: 2e-10,
            angular_tolerance_kg_m2_s: 2e-10,
            work_tolerance_j: 2e-10,
        };
        let (before, kick) = original
            .force_kick_candidate(0.001, &[], &[], [0.0; 3], [0.0; 3], impulse, cx)
            .unwrap();
        let states = before.solid().quadrature_states().to_vec();
        let old_velocity = before.velocity_mps().to_vec();
        let old_density = before.reference_density_kg_m3().to_vec();
        let total_mass = before.reference_mass_kg();
        let time = before.solid().time_s();
        let acceleration: Vec<_> = nodes
            .iter()
            .map(|x| [0.3 + x[0], -0.2 * x[2], 0.1 + 0.4 * x[1]])
            .collect();
        let loads = before
            .reference_acceleration_load_n(&acceleration, cx)
            .unwrap();
        let pressure = [bulk * (1.0 / reference.volume_ratio() - 1.0); 4];
        let admission = MixedViscoReferenceChangeControls {
            maximum_relative_volume_residual: 1e-11,
            maximum_pressure_constraint_work_j: 1e-10,
        };
        let change = before
            .change_reference(
                MixedViscoLawField::Uniform(&target),
                &pressure,
                admission,
                cx,
            )
            .unwrap();
        assert_eq!(before.solid().quadrature_states(), states);
        assert_eq!(before.velocity_mps(), old_velocity);
        near(
            change.stored_after_j - change.stored_before_j,
            change.reference_work_j - change.pressure_constraint_defect_j,
            2e-12,
        );
        near(
            change.reference_work_j,
            change.remainder_reference_work_j + change.bulk_reference_work_j,
            2e-12,
        );
        assert!(change.relative_volume_residual <= admission.maximum_relative_volume_residual);
        let error = change
            .motion
            .contact_drift_candidate(
                &kick,
                settings(),
                FiniteViscoContactDriftControls {
                    material_balance_tolerance_j: 2e-9,
                    maximum_integration_defect_j: 10.0,
                },
                cx,
            )
            .unwrap_err();
        assert!(error.to_string().contains("original kick"));
        drop(before);
        let mut motion = change.motion;
        assert_eq!(motion.velocity_mps(), old_velocity);
        assert_eq!(motion.solid().displacement_m(), u);
        assert_eq!(motion.solid().pressure_pa(), pressure);
        assert_eq!(motion.solid().time_s().to_bits(), time.to_bits());
        assert_eq!(motion.reference_mass_kg().to_bits(), total_mass.to_bits());
        assert_eq!(motion.reference_density_kg_m3(), old_density);
        assert_eq!(
            motion.reference_corner_density_kg_m3(),
            Some(density.as_slice())
        );
        assert_eq!(
            motion
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap(),
            loads
        );
        for (after, before) in motion.solid().quadrature_states().iter().zip(&states) {
            assert_eq!(
                after.inverse_viscous_cauchy_green,
                before.inverse_viscous_cauchy_green
            );
            assert_eq!(after.work.to_bits(), before.work.to_bits());
            assert_eq!(after.dissipated.to_bits(), before.dissipated.to_bits());
        }
        motion = motion
            .change_reference(MixedViscoLawField::Uniform(&old), &[0.0; 4], admission, cx)
            .unwrap()
            .motion;
        assert_eq!(motion.velocity_mps(), old_velocity);
        assert_eq!(
            motion
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap(),
            loads
        );
        let snapshot = motion.solid().quadrature_states().to_vec();
        assert!(
            motion
                .change_reference(
                    MixedViscoLawField::Uniform(&target),
                    &[0.0; 4],
                    admission,
                    cx
                )
                .is_err()
        );
        assert_eq!(motion.solid().quadrature_states(), snapshot);
        assert_eq!(motion.velocity_mps(), old_velocity);
    });
}
