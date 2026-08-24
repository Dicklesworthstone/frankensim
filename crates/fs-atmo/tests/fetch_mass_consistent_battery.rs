//! FetchAdjustedMassConsistent battery (bead `frankensim-wf-root-guzez.4.7`, E3.3c).

use fs_atmo::fetch_mass_consistent::{
    FetchAdjustedMassConsistent, FetchRoughnessProfile, MODEL_ID_FETCH_MASS_CONSISTENT,
};

#[test]
fn fetch_mass_consistent_divergence_cancels_to_machine_precision() {
    let profile = FetchRoughnessProfile::new(0.01, 0.0005, 500.0).expect("valid profile");
    let law = FetchAdjustedMassConsistent::new(profile, 0.0, 10.0, 15.0).expect("valid law");

    assert_eq!(MODEL_ID_FETCH_MASS_CONSISTENT, "FetchAdjustedMassConsistent");

    // Sample across grid of (x, h) points
    for ix in 0..20 {
        let x = ix as f64 * 25.0;
        for ih in 1..20 {
            let h = ih as f64 * 1.5;
            let div = law.divergence(x, h);
            assert!(
                div.abs() < 1e-12,
                "divergence at x={x}, h={h} must cancel to machine precision, got {div}"
            );
        }
    }
}

#[test]
fn fetch_mass_consistent_wall_impermeability() {
    let profile = FetchRoughnessProfile::new(0.02, 0.0002, 300.0).expect("valid profile");
    let law = FetchAdjustedMassConsistent::new(profile, 0.0, 10.0, 12.0).expect("valid law");

    for ix in 0..10 {
        let x = ix as f64 * 30.0;
        let v_ground = law.sample_velocity(x, 0.0);
        assert_eq!(v_ground[2], 0.0, "vertical velocity at ground must be exactly zero");
    }
}

#[test]
fn fetch_mass_consistent_refusals_and_boundaries() {
    // Non-finite roughness refuses
    assert!(FetchRoughnessProfile::new(f64::NAN, 0.001, 100.0).is_err());
    // z0 below minimum domain refuses
    assert!(FetchRoughnessProfile::new(1e-7, 0.001, 100.0).is_err());
    // z0 above maximum domain refuses
    assert!(FetchRoughnessProfile::new(2.5, 0.001, 100.0).is_err());
    // Non-positive max fetch refuses
    assert!(FetchRoughnessProfile::new(0.01, 0.001, 0.0).is_err());

    // Reference height below d + z0 refuses
    let profile = FetchRoughnessProfile::new(0.05, 0.0001, 100.0).expect("valid profile");
    assert!(FetchAdjustedMassConsistent::new(profile, 2.0, 2.02, 10.0).is_err());
}
