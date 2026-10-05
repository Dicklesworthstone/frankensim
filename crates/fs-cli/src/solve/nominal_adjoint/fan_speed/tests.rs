use super::*;
use fs_convection::ThermalDirection;

#[test]
fn native_fan_reynolds_elasticities_match_the_actual_card_formulas() {
    let step = 1e-5_f64;
    for (card, re, pr) in [
        (CorrelationId::CircularDuctLaminarCwt, 1000.0, 0.72),
        (CorrelationId::CircularDuctLaminarChf, 1000.0, 0.72),
        (CorrelationId::RectangularDuctLaminarCwt, 1000.0, 0.72),
        (CorrelationId::RectangularDuctLaminarChf, 1000.0, 0.72),
        (CorrelationId::CircularDuctHausen, 100.0, 0.72),
        (CorrelationId::CircularDuctHausen, 1000.0, 0.72),
        (CorrelationId::DittusBoelter, 50_000.0, 0.72),
        (CorrelationId::Gnielinski, 20_000.0, 0.72),
        (CorrelationId::Gnielinski, 20_000.0, 7.0),
        (CorrelationId::FlatPlateLaminarAverage, 100_000.0, 0.72),
        (CorrelationId::FlatPlateTurbulentAverage, 1_000_000.0, 0.72),
        (CorrelationId::ChurchillBernsteinCylinder, 10_000.0, 0.72),
    ] {
        let run = |r| evaluate(card, CorrelationInputs::forced(r,pr)
            .with_length_ratio(1000.0).with_aspect_ratio(0.5)
            .with_direction(ThermalDirection::HeatingFluid)).unwrap();
        let nominal = run(re);
        assert!(nominal.evidence().model.in_domain);
        let expected = (run(re*step.exp()).evidence().value.ln()
            - run(re*(-step).exp()).evidence().value.ln())/(2.0*step);
        let actual = reynolds_elasticity(&nominal).unwrap().unwrap();
        assert!((actual-expected).abs() < 1e-7, "{card:?}: {actual} vs {expected}");
    }
}

#[test]
fn native_fan_source_table_does_not_acquire_a_fabricated_smooth_slope() {
    let table = evaluate(CorrelationId::RectangularDuctLaminarCwtDevelopingPr072,
        CorrelationInputs::forced(1000.0,0.72).with_length_ratio(24.0)
            .with_aspect_ratio(0.5)).unwrap();
    assert_eq!(reynolds_elasticity(&table).unwrap(),None);
}
