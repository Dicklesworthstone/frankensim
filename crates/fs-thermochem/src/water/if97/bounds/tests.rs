use super::*;
use crate::water::if97::If97State;

const CONTROLS: If97HeatCapacityBoundSolve = If97HeatCapacityBoundSolve {
    maximum_intervals: 8191,
};

#[test]
fn complete_liquid_and_steam_isobars_have_conservative_useful_caloric_bounds() {
    for (region, p, lo, hi) in [
        (If97Region::Region1, 100_000.0, 273.15, 370.0),
        (If97Region::Region1, 3e6, 300.0, 500.0),
        (If97Region::Region1, 100e6, 273.15, 623.15),
        (If97Region::Region2, 3500.0, 300.0, 1073.15),
        (If97Region::Region2, 100_000.0, 373.2, 1073.15),
        (If97Region::Region2, 30e6, 700.0, 1073.15),
        (If97Region::Region2, 100e6, 863.16, 1073.15),
    ] {
        let bound = isobaric_heat_capacity_bound(
            Pressure::new(p),
            region,
            Temperature::new(lo),
            Temperature::new(hi),
            CONTROLS,
            || false,
        )
        .unwrap();
        assert!(bound.minimum_heat_capacity_j_kg_k() > 0.0);
        assert!(bound.interval_evaluations() <= CONTROLS.maximum_intervals);
        let mut sampled_minimum = f64::INFINITY;
        for i in 0..=128 {
            let t = lo + (hi - lo) * f64::from(i) / 128.0;
            let state =
                If97State::at_pressure_temperature(Pressure::new(p), Temperature::new(t), region)
                    .unwrap();
            let cp = state.cp().value();
            sampled_minimum = sampled_minimum.min(cp);
            assert!(
                cp >= bound.minimum_heat_capacity_j_kg_k(),
                "{region:?} {p} {t}: {cp} {bound:?}"
            );
            assert!(cp <= bound.maximum_heat_capacity_j_kg_k());
            assert!(bound.maximum_temperature_slope() >= 1.0 / cp);
        }
        // Samples validate the implementation and utility; the interval algebra
        // itself proves the bound also at every unsampled temperature.
        assert!(bound.minimum_heat_capacity_j_kg_k() >= sampled_minimum * 0.49);
    }
}

#[test]
fn zero_width_enclosures_cover_the_original_published_reference_states() {
    for (region, p, t) in [
        (If97Region::Region1, 3e6, 300.0),
        (If97Region::Region1, 80e6, 300.0),
        (If97Region::Region1, 3e6, 500.0),
        (If97Region::Region2, 3500.0, 300.0),
        (If97Region::Region2, 3500.0, 700.0),
        (If97Region::Region2, 30e6, 700.0),
    ] {
        let state =
            If97State::at_pressure_temperature(Pressure::new(p), Temperature::new(t), region)
                .unwrap();
        let bound = isobaric_heat_capacity_bound(
            Pressure::new(p),
            region,
            Temperature::new(t),
            Temperature::new(t),
            CONTROLS,
            || false,
        )
        .unwrap();
        assert!(bound.minimum_heat_capacity_j_kg_k() <= state.cp().value());
        assert!(bound.maximum_heat_capacity_j_kg_k() >= state.cp().value());
        assert_eq!(bound.interval_evaluations(), 1);
    }
}

#[test]
fn exhausted_bounds_and_late_cancellation_publish_no_partial_cfl_allowance() {
    let invoke = |controls, cancelled: fn() -> bool| {
        isobaric_heat_capacity_bound(
            Pressure::new(100e6),
            If97Region::Region1,
            Temperature::new(273.15),
            Temperature::new(623.15),
            controls,
            cancelled,
        )
    };
    assert!(matches!(
        invoke(
            If97HeatCapacityBoundSolve {
                maximum_intervals: 1
            },
            || false
        ),
        Err(If97HeatCapacityBoundError::NotProven {
            interval_evaluations: 1,
            ..
        })
    ));
    let receipt = invoke(CONTROLS, || false).unwrap();
    let mut polls = 0;
    let result = isobaric_heat_capacity_bound(
        Pressure::new(100e6),
        If97Region::Region1,
        Temperature::new(273.15),
        Temperature::new(623.15),
        CONTROLS,
        || {
            polls += 1;
            polls == receipt.interval_evaluations() + 2
        },
    );
    assert_eq!(result, Err(If97HeatCapacityBoundError::Cancelled));
    assert_eq!(polls, receipt.interval_evaluations() + 2);
}

#[test]
fn bounds_refuse_invalid_work_and_unsupported_regional_paths() {
    assert!(matches!(
        isobaric_heat_capacity_bound(
            Pressure::new(3e6),
            If97Region::Region1,
            Temperature::new(300.0),
            Temperature::new(500.0),
            If97HeatCapacityBoundSolve {
                maximum_intervals: 0
            },
            || false,
        ),
        Err(If97HeatCapacityBoundError::Invalid(_))
    ));
    assert!(matches!(
        isobaric_heat_capacity_bound(
            Pressure::new(100e6),
            If97Region::Region1,
            Temperature::new(300.0),
            Temperature::new(700.0),
            CONTROLS,
            || false,
        ),
        Err(If97HeatCapacityBoundError::Source(_))
    ));
    assert!(matches!(
        isobaric_heat_capacity_bound(
            Pressure::new(3e6),
            If97Region::Region1,
            Temperature::new(500.0),
            Temperature::new(300.0),
            CONTROLS,
            || false,
        ),
        Err(If97HeatCapacityBoundError::Invalid(_))
    ));
}
