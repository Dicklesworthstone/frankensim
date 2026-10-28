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
