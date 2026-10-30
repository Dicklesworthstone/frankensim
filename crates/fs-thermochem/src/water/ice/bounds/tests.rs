use super::*;

fn bound(p: f64, lo: f64, hi: f64) -> IceIhHeatCapacityBound {
    isobaric_heat_capacity_bound(
        Pressure::new(p),
        Temperature::new(lo),
        Temperature::new(hi),
        IceIhHeatCapacityBoundSolve {
            maximum_intervals: 8191,
        },
        || false,
    )
    .unwrap()
}

#[test]
fn whole_interval_bounds_cover_original_cp_and_caloric_changes() {
    for p in [0.0, 101_325.0, 100e6, 210e6] {
        for (lo, hi) in [(50.0, 273.16), (243.15, 273.16), (273.15, 273.16)] {
            let proof = bound(p, lo, hi);
            let (minimum, maximum) = (
                proof.minimum_heat_capacity_j_kg_k(),
                proof.maximum_heat_capacity_j_kg_k(),
            );
            assert!(minimum > 0.0 && maximum.is_finite());
            for i in 0..=128 {
                let t = lo + (hi - lo) * f64::from(i) / 128.0;
                let state =
                    IceIhState::at_pressure_temperature(Pressure::new(p), Temperature::new(t))
                        .unwrap();
                let cp = state.cp().value();
                assert!(cp >= minimum && cp <= maximum, "p={p}, T={t}, {proof:?}");
                assert!(1.0 / cp <= proof.maximum_temperature_slope());
            }
            let delta = super::super::isobaric_difference(
                Pressure::new(p),
                Temperature::new(lo),
                Temperature::new(hi),
            )
            .unwrap();
            let dh = delta.specific_enthalpy_change().value();
            let ds = delta.specific_entropy_change().value();
            assert!(dh >= minimum * (hi - lo) && dh <= maximum * (hi - lo));
            assert!(ds >= minimum * (hi / lo).ln() && ds <= maximum * (hi / lo).ln());
            assert_eq!(proof, bound(p, lo, hi));
        }
    }
    let point = bound(101_325.0, 250.0, 250.0);
    assert_eq!(point.interval_evaluations(), 1);
}

#[test]
fn unresolved_domain_and_cancelled_bounds_publish_no_partial_result() {
    let run = |lo, hi, work| {
        isobaric_heat_capacity_bound(
            Pressure::new(101_325.0),
            Temperature::new(lo),
            Temperature::new(hi),
            IceIhHeatCapacityBoundSolve {
                maximum_intervals: work,
            },
            || false,
        )
    };
    assert!(matches!(
        run(50.0, 273.16, 1),
        Err(IceIhHeatCapacityBoundError::NotProven { .. })
    ));
    for (lo, hi, work) in [
        (250.0, 273.16, 0),
        (250.0, 273.16, MAX_BOUND_INTERVALS + 1),
        (273.16, 250.0, 8191),
        (f64::NAN, 250.0, 8191),
    ] {
        assert!(matches!(
            run(lo, hi, work),
            Err(IceIhHeatCapacityBoundError::Invalid(_))
        ));
    }
    for (lo, hi) in [(49.0, 273.16), (250.0, 274.0)] {
        assert!(matches!(
            run(lo, hi, 8191),
            Err(IceIhHeatCapacityBoundError::Source(_))
        ));
    }
    for after in [0, 2, 5] {
        let mut polls = 0;
        let result = isobaric_heat_capacity_bound(
            Pressure::new(101_325.0),
            Temperature::new(50.0),
            Temperature::new(273.16),
            IceIhHeatCapacityBoundSolve {
                maximum_intervals: 8191,
            },
            || {
                let cancel = polls >= after;
                polls += 1;
                cancel
            },
        );
        assert_eq!(result, Err(IceIhHeatCapacityBoundError::Cancelled));
    }
}
