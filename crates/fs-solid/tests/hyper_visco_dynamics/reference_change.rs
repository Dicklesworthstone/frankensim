//! Prescribed eigenstrain changes retain the original graded momentum owner.

use super::*;
use fs_material::hyper::reference::ReferenceDeformation;
use fs_solid::hyper3::visco::dynamics::{FiniteViscoContactDriftControls, SolidImpulseControls};
use fs_solver::CheckedCgConfig;

#[test]
fn reference_change_preserves_graded_mass_and_motion_and_invalidates_old_force_receipts() {
    with_cx(|cx| {
        let old = material(true);
        let target = old.clone().with_reference_deformation(
            ReferenceDeformation::try_new([1.02, 0.03, 0.0, 0.0, 0.99, 0.0, 0.0, 0.0, 1.01])
                .unwrap(),
        );
        let density = [30.0, 150.0, 30.0, 30.0];
        let u = [[0.0; 3], [0.02, 0.0, 0.0], [0.0; 3], [0.0; 3]];
        let prescribed: Vec<_> = u.as_flattened().iter().copied().enumerate().collect();
        let solid = FiniteViscoTetSolid::try_new_from(
            &NODES,
            &TETS,
            ViscoelasticLawField::Uniform(&old),
            &u,
            TetAssemblyBudget::default(),
            cx,
        )
        .unwrap()
        .advance(
            0.007,
            &prescribed,
            &[],
            [0.0; 3],
            settings(),
            controls().material,
            cx,
        )
        .unwrap()
        .solid;
        assert!(solid.element_states()[0].dissipated > 0.0);
        let velocity = [
            [0.2, -0.1, 0.05],
            [-0.1, 0.3, 0.02],
            [0.04, -0.2, 0.1],
            [0.3, 0.02, -0.1],
        ];
        let original =
            FiniteViscoTetMotion::try_new_with_nodal_density(solid, &density, &velocity, cx)
                .unwrap();
        let impulse = SolidImpulseControls {
            solve: CheckedCgConfig {
                tolerance: 1e-13,
                max_iterations: 64,
                max_corrections: 2,
            },
            maximum_workspace_components: 4096,
            linear_tolerance_kg_m_s: 1e-10,
            angular_tolerance_kg_m2_s: 1e-10,
            work_tolerance_j: 1e-10,
        };
        let (before, kick) = original
            .force_kick_candidate(0.001, &[], &[], [0.0; 3], [0.0; 3], impulse, cx)
            .unwrap();
        let old_velocity = before.velocity_mps().to_vec();
        let states = before.solid().element_states().to_vec();
        let time = before.solid().time_s();
        let total_mass = before.reference_mass_kg();
        let acceleration = [
            [0.3, -0.2, 0.1],
            [0.7, 0.05, -0.3],
            [-0.1, 0.4, 0.2],
            [0.2, 0.1, -0.2],
        ];
        let loads = before
            .reference_acceleration_load_n(&acceleration, cx)
            .unwrap();
        let change = before
            .change_reference(ViscoelasticLawField::Uniform(&target), cx)
            .unwrap();
        assert_eq!(before.solid().element_states(), states);
        assert_eq!(before.velocity_mps(), old_velocity);
        close(
            change.stored_after_j - change.stored_before_j,
            change.reference_work_j,
            2e-12,
        );
        assert!(change.reference_work_j.abs() > 1e-8);
        let error = change
            .motion
            .contact_drift_candidate(
                &kick,
                FiniteViscoContactDriftControls {
                    material_balance_tolerance_j: 2e-9,
                    maximum_integration_defect_j: 10.0,
                },
                cx,
            )
            .unwrap_err();
        assert!(error.to_string().contains("original kick"));
        // The successor borrows original geometry/materials, never the wrapper
        // that produced it. This drop is also a compile-time lifetime regression.
        drop(before);
        let mut motion = change.motion;
        assert_eq!(motion.velocity_mps(), old_velocity);
        assert_eq!(motion.solid().time_s().to_bits(), time.to_bits());
        assert_eq!(motion.solid().displacement_m(), u);
        assert_eq!(motion.reference_mass_kg().to_bits(), total_mass.to_bits());
        assert_eq!(motion.reference_density_kg_m3(), &[60.0]);
        assert_eq!(
            motion.reference_nodal_density_kg_m3(),
            Some(density.as_slice())
        );
        assert_eq!(
            motion
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap(),
            loads
        );
        for (after, before) in motion.solid().element_states().iter().zip(&states) {
            assert_eq!(
                after.inverse_viscous_cauchy_green,
                before.inverse_viscous_cauchy_green
            );
            assert_eq!(after.work.to_bits(), before.work.to_bits());
            assert_eq!(after.dissipated.to_bits(), before.dissipated.to_bits());
        }
        // Repeated updates can replace the current owner in a trajectory.
        motion = motion
            .change_reference(ViscoelasticLawField::Uniform(&old), cx)
            .unwrap()
            .motion;
        assert_eq!(motion.velocity_mps(), old_velocity);
        assert_eq!(
            motion
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap(),
            loads
        );
        close(
            motion.solid().element_states()[0].reference_work,
            0.0,
            2e-10,
        );
        let mut incompatible = target.clone();
        incompatible.branches[0].energy = BranchEnergy::NeoHookean {
            mu_pa: 2.0 * BRANCH_MU,
        };
        let snapshot = motion.solid().element_states().to_vec();
        assert!(
            motion
                .change_reference(ViscoelasticLawField::Uniform(&incompatible), cx)
                .is_err()
        );
        assert_eq!(motion.solid().element_states(), snapshot);
        assert_eq!(motion.velocity_mps(), old_velocity);
    });
}
