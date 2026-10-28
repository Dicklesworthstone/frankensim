//! Actual building damping direction and gravity-tangent stability.

use fs_frame::{BuildingFrame, BuildingSpec, RcMember};
use fs_solid::frame2d::{FrameMassFormulation, FrameModeStatus};

#[test]
fn horizontal_damping_does_not_use_a_lower_vertical_member_mode() {
    let spec = BuildingSpec {
        stories: 1,
        bays: 1,
        bay_width: 40.0,
        floor_mass: 100.0,
        column: RcMember {
            depth: 0.6,
            width: 0.6,
            steel_area: 1.5e-3,
        },
        beam: RcMember {
            depth: 0.35,
            width: 0.25,
            steel_area: 2e-3,
        },
        ..BuildingSpec::default()
    };
    let building = BuildingFrame::with_member_mass(
        spec,
        864.0,
        210.0,
        FrameMassFormulation::ConsistentReference,
    )
    .unwrap();
    let analysis = building.modal_analysis();
    let floor = 256.0 * f64::EPSILON * analysis.base_influence_mass_kg[0];
    let lowest = &analysis.modes[0];
    assert_eq!(lowest.status(), FrameModeStatus::Oscillatory);
    assert!(
        lowest.effective_mass_kg[0] < floor,
        "lowest mode: {lowest:?}"
    );
    assert!(lowest.effective_mass_kg[1] > 1.0);
    let horizontal = analysis
        .modes
        .iter()
        .find(|m| m.effective_mass_kg[0] > floor)
        .unwrap();
    let omega = horizontal.angular_frequency().unwrap();
    let damping = building.damping();
    let ratio = 0.5 * (damping.a0 / omega + damping.a1 * omega);
    assert!((ratio - spec.zeta).abs() < 1e-12);
    assert!((building.periods()[0] - 2.0 * std::f64::consts::PI / omega).abs() < 1e-12);
    let vertical_omega = lowest.angular_frequency().unwrap();
    let vertical_ratio = 0.5 * (damping.a0 / vertical_omega + damping.a1 * vertical_omega);
    assert!((vertical_ratio - spec.zeta).abs() > 1e-3);
}

#[test]
fn gravity_buckling_cannot_be_discarded_before_damping_fit() {
    let spec = BuildingSpec {
        stories: 1,
        bays: 1,
        story_height: 30.0,
        floor_mass: 10_000.0,
        ..BuildingSpec::default()
    };
    let stable = BuildingFrame::new(spec).unwrap();
    assert!(
        stable
            .modal_analysis()
            .modes
            .iter()
            .all(|m| m.status() == FrameModeStatus::Oscillatory)
    );
    let error = BuildingFrame::new(BuildingSpec {
        floor_mass: 500_000.0,
        ..spec
    })
    .unwrap_err();
    match error {
        fs_solid::SolidError::InvalidInput { what } => {
            assert!(
                what.contains("gravity-loaded tangent") && what.contains("Unstable"),
                "{what}"
            );
        }
        other => panic!("expected signed gravity instability, got {other:?}"),
    }
}
