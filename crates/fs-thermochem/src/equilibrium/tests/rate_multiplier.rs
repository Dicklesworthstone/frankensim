//! Additional absence/range cases on the single native kinetic multiplier.
use super::*;

#[test]
fn zero_multiplier_stops_log_only_rates_without_discarding_activity_or_source_receipts() {
    let model = water_mass_action(7.0);
    let temperature = Temperature::new(900.0);
    let activities = vec![
        activity("H2", 1e200),
        activity("O2", 1e200),
        activity("H2O", 1e200),
    ];
    let base = model.evaluate(temperature, activities.clone()).unwrap();
    assert!(base.forward_activity_product().is_none());
    assert!(base.reverse_activity_product().is_none());
    let scaled = model
        .evaluate_with_kinetic_multiplier(
            temperature,
            activities.clone(),
            Dimensionless::new(1e-300),
        )
        .unwrap();
    // A finite factor cannot recover a base direct rate that is log-only.
    assert_eq!(scaled.forward_status(), base.forward_status());
    assert_eq!(scaled.reverse_status(), base.reverse_status());
    assert!(scaled.forward_progress_rate().is_none());
    assert!(scaled.reverse_progress_rate().is_none());
    assert!(scaled.net_progress_rate().is_none());
    let absent = model
        .evaluate_with_kinetic_multiplier(temperature, activities, Dimensionless::new(0.0))
        .unwrap();
    assert_eq!(absent.reverse_rate(), base.reverse_rate());
    assert_eq!(absent.receipt().terms(), base.receipt().terms());
    assert_eq!(
        absent.log_forward_activity_product(),
        base.log_forward_activity_product()
    );
    assert_eq!(
        absent.log_reverse_activity_product(),
        base.log_reverse_activity_product()
    );
    assert!(absent.forward_activity_product().is_none());
    assert!(absent.reverse_activity_product().is_none());
    assert_eq!(absent.forward_progress_rate().unwrap().value(), 0.0);
    assert_eq!(absent.reverse_progress_rate().unwrap().value(), 0.0);
    assert_eq!(absent.net_progress_rate().unwrap().value(), 0.0);
    assert_eq!(
        absent.forward_status(),
        MassActionDirectionalRateStatusV1::ZeroKineticMultiplier
    );
    assert_eq!(
        absent.reverse_status(),
        MassActionDirectionalRateStatusV1::ZeroKineticMultiplier
    );
    assert_eq!(
        absent.receipt().kinetic_multiplier_bits(),
        0.0_f64.to_bits()
    );
}

#[test]
fn exact_species_absence_keeps_its_own_status_when_the_common_multiplier_is_zero() {
    let model = water_mass_action(7.0);
    let temperature = Temperature::new(900.0);
    let activities = vec![
        activity("H2", 0.0),
        activity("O2", 1.0),
        activity("H2O", 1.0),
    ];
    let base = model.evaluate(temperature, activities.clone()).unwrap();
    let stopped = model
        .evaluate_with_kinetic_multiplier(temperature, activities, Dimensionless::new(0.0))
        .unwrap();
    assert_eq!(
        stopped.forward_status(),
        MassActionDirectionalRateStatusV1::ZeroActivity
    );
    assert_eq!(
        stopped.reverse_status(),
        MassActionDirectionalRateStatusV1::ZeroKineticMultiplier
    );
    assert_eq!(stopped.receipt().terms(), base.receipt().terms());
    assert_eq!(
        stopped.forward_activity_product(),
        base.forward_activity_product()
    );
    assert_eq!(
        stopped.reverse_activity_product(),
        base.reverse_activity_product()
    );
    assert_eq!(stopped.forward_progress_rate().unwrap().value(), 0.0);
    assert_eq!(stopped.reverse_progress_rate().unwrap().value(), 0.0);
    assert_eq!(stopped.net_progress_rate().unwrap().value(), 0.0);
}
