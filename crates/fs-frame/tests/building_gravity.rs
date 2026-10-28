//! Gravity remains an external load during the full building time history.

use fs_frame::{BuildingFrame, BuildingSpec};

#[test]
fn a_gravity_loaded_building_stays_at_rest_without_base_excitation() {
    let spec = BuildingSpec {
        stories: 2,
        bays: 1,
        ..BuildingSpec::default()
    };
    let mut building = BuildingFrame::new(spec).unwrap();
    let response = building.run(&[0.0; 101], 0.01).unwrap();
    // Dropping the sustained weight releases axial strain energy even when
    // symmetry hides that false vertical vibration from the lateral drift QoI.
    assert!(response.internal_work.abs() < 1e-5, "{:?}", response);
    assert!(response.peak_drift_ratio < 1e-10, "{:?}", response);
}

#[test]
fn member_inertia_enters_loaded_building_modes_and_motion() {
    use fs_solid::frame2d::FrameMassFormulation;
    let spec = BuildingSpec {
        stories: 1,
        bays: 1,
        floor_mass: 20_000.0,
        ..BuildingSpec::default()
    };
    let point_mass_only = BuildingFrame::new(spec).unwrap();
    let loaded = BuildingFrame::with_member_mass(
        spec,
        500.0,
        400.0,
        FrameMassFormulation::ConsistentReference,
    )
    .unwrap();
    // Added self-weight changes the nonlinear RC tangent as well as inertia;
    // a period increase is not required when both stiffness and mass change.
    assert!((loaded.periods()[0] / point_mass_only.periods()[0] - 1.0).abs() > 1e-3);
    println!(
        "floor-only period {}, with member inertia {}",
        point_mass_only.periods()[0],
        loaded.periods()[0]
    );

    let rest = loaded.clone().run(&[0.0; 101], 0.01).unwrap();
    assert!(rest.internal_work.abs() < 1e-4, "{rest:?}");
    let ground: Vec<f64> = (0..=100)
        .map(|i| {
            let t = i as f64 * 0.01;
            0.3 * (4.0 * std::f64::consts::PI * t).sin()
        })
        .collect();
    let response = loaded.clone().run(&ground, 0.01).unwrap();
    assert!(response.peak_drift_ratio > 1e-5 && response.peak_drift_ratio < 0.02);
    assert!(response.energy_balance_error < 1e-6, "{response:?}");
    assert!(response.peak_base_shear.is_finite() && response.peak_base_shear > 0.0);
}

#[test]
fn zero_member_mass_keeps_existing_building_behavior() {
    use fs_solid::frame2d::FrameMassFormulation;
    let spec = BuildingSpec {
        stories: 1,
        bays: 1,
        ..BuildingSpec::default()
    };
    let mut original = BuildingFrame::new(spec).unwrap();
    let mut explicit =
        BuildingFrame::with_member_mass(spec, 0.0, 0.0, FrameMassFormulation::ConsistentReference)
            .unwrap();
    assert_eq!(original.periods(), explicit.periods());
    let ground: Vec<f64> = (0..=30).map(|i| 0.2 * (i as f64 * 0.2).sin()).collect();
    assert_eq!(
        original.run(&ground, 0.01).unwrap(),
        explicit.run(&ground, 0.01).unwrap()
    );
    assert!(
        BuildingFrame::with_member_mass(
            spec,
            f64::NAN,
            0.0,
            FrameMassFormulation::ConsistentReference
        )
        .is_err()
    );
}
