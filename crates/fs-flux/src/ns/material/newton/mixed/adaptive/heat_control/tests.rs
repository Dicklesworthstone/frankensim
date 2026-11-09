//! G3/G4: non-cancelling physical errors and accepted absolute-defect accounts.

use super::*;

fn controls(absolute: f64) -> MaterialSdirk2HeatControls {
    MaterialSdirk2HeatControls {
        absolute_temporal_defect_budget_j_per_m: absolute,
        relative_temporal_defect_budget: 0.0,
        ..MaterialSdirk2HeatControls::default()
    }
}
fn trial(defect: f64) -> HeatTrial {
    HeatTrial {
        signed: Dd::from_f64(defect),
        absolute: Dd::from_f64(defect.abs()),
        heat: Dd::ZERO,
        cell: 0.0,
        edge: 0.0,
        bulk_l1: 0.0,
        slip_l1: 0.0,
        local_ratio: 0.0,
        defect_ratio: 0.0,
    }
}

#[test]
fn g3_opposing_original_spatial_heat_discrepancies_cannot_cancel() {
    let tolerance = MaterialSdirk2HeatTolerance {
        absolute_j_per_m: 0.0,
        relative: 0.1,
    };
    let mut error = SpatialHeatError::default();
    // Coarse totals and fine totals are both nine, but each original spatial
    // owner is wrong by three. The second-order estimate is one on each.
    error.add(3.0, 3.0, 3.0, tolerance).unwrap();
    error.add(6.0, 1.0, 2.0, tolerance).unwrap();
    assert_eq!(error.maximum, 1.0);
    assert_eq!(error.l1.to_f64(), 2.0);
    assert!((error.ratio - 1.0 / 0.6).abs() < 1e-14);
}

#[test]
fn g4_time_defect_signs_and_sub_ulp_increments_never_create_budget_credit() {
    let mut cancelling = HeatControl::new(controls(0.015), 1.0).unwrap();
    cancelling.accept(trial(0.01)).unwrap();
    cancelling.accept(trial(-0.01)).unwrap();
    assert_eq!(cancelling.signed.to_f64(), 0.0);
    assert!(cancelling.finish().is_err());

    let mut exhausted = HeatControl::new(controls(0.01), 1.0).unwrap();
    exhausted.accept(trial(0.01)).unwrap();
    exhausted.accept(trial(1e-20)).unwrap();
    assert_eq!(0.01 + 1e-20, 0.01);
    assert!(exhausted.finish().is_err());
    assert!(ratio(1e-20, 0.0).is_infinite());

    let mut quiet = HeatControl::new(controls(0.0), 1.0).unwrap();
    quiet.accept(trial(0.0)).unwrap();
    assert_eq!(
        quiet
            .finish()
            .unwrap()
            .absolute_temporal_kinetic_defect_j_per_m,
        0.0
    );
    assert_eq!(ratio(0.0, 0.0), 0.0);
}
