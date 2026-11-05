//! G1/G3: original hardness selects the spring-loaded plastic response.

use super::*;
use fs_contact::normal_patch::plastic_source::drive::SourcePlasticDriveError;
use fs_contact::normal_patch::{FiniteGapPlasticDrive, FiniteGapPlasticLoadControls};

fn controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1.0e-4,
        maximum_approach_m: 1.0e-4,
        absolute_force_tolerance_n: 1.0e-14,
        relative_force_tolerance: 1.0e-12,
        approach_tolerance_m: 1.0e-18,
        maximum_equilibrium_evaluations: 128,
    }
}

fn compliance(source: &SourcePlasticContact) -> f64 {
    // Integral of the rectangular Boussinesq kernel over one square cell.
    let width = source.model().config().grid.cell_width_m;
    let modulus = source.normal_source().elastic().reduced_modulus_pa();
    4.0 * width * (1.0 + 2.0_f64.sqrt()).ln() / (core::f64::consts::PI * modulus)
}

#[test]
fn source_hardness_and_device_stiffness_determine_retained_indentation() {
    with_cx(|cx, _| {
        let (normal, hardness) = sources([1.0e6, 2.0e6], "spring", 0.0);
        let source = build(normal, hardness, envelope(), cx).unwrap();
        let c = compliance(&source);
        let h = 1.0e6;
        let area = 1.0e-8;
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: area / c,
            position_m: 3.0 * c * h,
        };
        let before = source.initial_state();
        let loaded = source
            .advance_with_drive(&before, drive, controls(), cx)
            .unwrap();
        let native = loaded.candidate.contact_state();
        common::close(loaded.drive.achieved_force_n, area * h);
        common::close(loaded.drive.drive_force_n, area * h);
        common::close(native.approach_m(), 2.0 * c * h);
        common::close(native.plastic_depth_m()[12], c * h);
        common::close(loaded.receipt.plastic_work_j, area * h * c * h);
        common::close(loaded.drive.spring_energy_j, 0.5 * area * c * h * h);
        assert_eq!(loaded.source_identity, source.identity());
        assert_eq!(before, source.initial_state());

        let retained = source
            .advance_with_drive(&loaded.candidate, drive, controls(), cx)
            .unwrap();
        assert_eq!(retained.candidate, loaded.candidate);
        assert_eq!(retained.receipt.plastic_work_j, 0.0);
        let unloaded = source
            .advance_with_drive(
                &loaded.candidate,
                FiniteGapPlasticDrive {
                    position_m: 0.0,
                    ..drive
                },
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(unloaded.drive.achieved_force_n, 0.0);
        assert_eq!(unloaded.drive.spring_energy_j, 0.0);
        assert_eq!(
            unloaded.candidate.contact_state().plastic_depth_m(),
            native.plastic_depth_m()
        );
    });
}

#[test]
fn source_drive_preserves_original_test_owner_and_physical_bounds() {
    with_cx(|cx, gate| {
        let (normal, hardness) = sources([1.0e6, 2.0e6], "spring-owner", 0.0);
        let source = build(normal, hardness, envelope(), cx).unwrap();
        let c = compliance(&source);
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: 1.0e-8 / c,
            position_m: 3.0 * c * 1.0e6,
        };
        let before = source.initial_state();
        let saved = before.clone();
        let (normal, hardness) = sources([1.0e6, 2.0e6], "different-test", 0.0);
        let foreign = build(normal, hardness, envelope(), cx).unwrap();
        assert!(matches!(
            foreign.advance_with_drive(&before, drive, controls(), cx),
            Err(SourcePlasticDriveError::Contact(
                SourcePlasticContactError::Invalid(_)
            ))
        ));
        let excessive = FiniteGapPlasticLoadControls {
            maximum_approach_m: 1.01e-4,
            ..controls()
        };
        assert!(matches!(
            source.advance_with_drive(&before, drive, excessive, cx),
            Err(SourcePlasticDriveError::Contact(
                SourcePlasticContactError::OutsideApplicability { .. }
            ))
        ));
        let (normal, hardness) = sources([1.0e6, 2.0e6], "bounded-depth", 0.0);
        let bounded = build(
            normal,
            hardness,
            PlasticContactEnvelope {
                maximum_plastic_depth_m: 0.0,
                ..envelope()
            },
            cx,
        )
        .unwrap();
        assert!(matches!(
            bounded.advance_with_drive(&bounded.initial_state(), drive, controls(), cx),
            Err(SourcePlasticDriveError::Contact(
                SourcePlasticContactError::OutsideApplicability { .. }
            ))
        ));
        gate.request();
        assert!(
            source
                .advance_with_drive(&before, drive, controls(), cx)
                .is_err()
        );
        assert_eq!(before, saved);
    });
}
