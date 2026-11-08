//! Adaptive source closure against a stiff, strictly passive analytical law.
use super::*;
#[path = "../../common/adaptive_conductor.rs"]
mod adaptive_source;

fn adaptation() -> AdaptiveConductionRelaxation {
    AdaptiveConductionRelaxation {
        initial_weight: 0.5,
        minimum_weight: 0.02,
        maximum_weight: 1.0,
    }
}

#[test]
fn adaptive_resistance_closes_stiff_voltage_and_total_current_drives() {
    // In the middle source interval rho(J)=40J-35. At E=5 the unique
    // positive solution is J=1, rho=5. Picard's derivative there is -8;
    // even fixed half relaxation has derivative -3.5 and fails to converge.
    // The conductivity/E dual at J=5 has the same independent root E=1.
    for current_control in [false, true] {
        let law = adaptive_source::law(current_control);
        let original = specimen(law.clone(), law, 0.2);
        let fixed = original.clone().with_relaxation(0.5).unwrap();
        let adaptive = original
            .clone()
            .with_adaptive_relaxation(adaptation())
            .unwrap();
        for sign in [-1.0, 1.0] {
            let voltage = ends(sign * 5.0);
            let solve = |model: &PlanarConductionResponse| {
                if current_control {
                    model.solve_current_candidate(
                        &[2, 5],
                        &[0, 3],
                        sign,
                        limits(),
                        controls(),
                        policy(),
                        &mut || false,
                    )
                } else {
                    model.solve_candidate(
                        &voltage,
                        &[],
                        limits(),
                        controls(),
                        policy(),
                        &mut || false,
                    )
                }
            };
            assert!(matches!(
                solve(&fixed),
                Err(ConductionResponseError::NotConverged { .. })
            ));
            let candidate = solve(&adaptive).unwrap();
            let report = candidate.report();
            assert!(report.field_solves < 32);
            let (e, j) = if current_control {
                (1.0, 5.0)
            } else {
                (5.0, 1.0)
            };
            for cell in 0..4 {
                near(report.electric_field_v_m[cell][0], -sign * e);
                near(report.current_density_a_m2[cell][0], -sign * j);
                near(report.cell_power_w[cell], 0.25);
                near(
                    report.evaluated[cell]
                        .conductor
                        .material()
                        .state_coordinate("operating")
                        .unwrap(),
                    1.0,
                );
                assert_eq!(
                    report.evaluated[cell]
                        .conductor
                        .source_property()
                        .answer()
                        .receipt
                        .selected,
                    original.conductors()[cell]
                        .source_property()
                        .answer()
                        .receipt
                        .selected,
                );
            }
            near(report.voltage_v[2], sign * e);
            near(report.power_w, 1.0);
            near(report.power_balance_residual_w, 0.0);
            near(report.boundary_current_a, 0.0);
            assert_eq!(candidate.problem().adaptive_relaxation(), Some(adaptation()));
            assert_eq!(candidate.problem().relaxation(), adaptation().initial_weight);

            let (native, sources) = candidate.into_parts().0.into_frozen_parts();
            assert_eq!(sources.adaptive_relaxation(), Some(adaptation()));
            let continued =
                PlanarConductionResponse::from_frozen_parts(native, sources, &mut || false)
                    .unwrap();
            assert_eq!(continued.adaptive_relaxation(), Some(adaptation()));
            near(solve(&continued).unwrap().report().power_w, 1.0);
        }
        assert_eq!(adaptive.conductors(), original.conductors());
        assert_eq!(
            adaptive.with_relaxation(0.5).unwrap().adaptive_relaxation(),
            None,
        );
    }
}

#[test]
fn adaptive_source_refusals_and_interruption_leave_no_iteration_history() {
    let law = adaptive_source::law(false);
    let original = specimen(law.clone(), law, 0.2);
    for invalid in [
        AdaptiveConductionRelaxation {
            minimum_weight: 0.0,
            ..adaptation()
        },
        AdaptiveConductionRelaxation {
            initial_weight: 0.01,
            ..adaptation()
        },
        AdaptiveConductionRelaxation {
            maximum_weight: 1.1,
            ..adaptation()
        },
        AdaptiveConductionRelaxation {
            initial_weight: f64::NAN,
            ..adaptation()
        },
    ] {
        assert!(original.clone().with_adaptive_relaxation(invalid).is_err());
    }
    let model = original
        .clone()
        .with_adaptive_relaxation(adaptation())
        .unwrap();
    let short = ConductionResponseControls {
        maximum_iterations: 2,
        ..policy()
    };
    assert!(matches!(
        model.solve_candidate(&ends(5.0), &[], limits(), controls(), short, &mut || false),
        Err(ConductionResponseError::NotConverged { iterations: 2, .. }),
    ));
    let mut calls = 0;
    let accepted = model
        .solve_candidate(&ends(5.0), &[], limits(), controls(), policy(), &mut || {
            calls += 1;
            false
        })
        .unwrap();
    let mut interrupted = 0;
    assert!(matches!(
        model.solve_candidate(&ends(5.0), &[], limits(), controls(), policy(), &mut || {
            interrupted += 1;
            interrupted == calls / 2
        }),
        Err(ConductionResponseError::Cancelled),
    ));
    assert_eq!(model.conductors(), original.conductors());
    let retry = model
        .solve_candidate(&ends(5.0), &[], limits(), controls(), policy(), &mut || false)
        .unwrap();
    assert_eq!(retry.report().voltage_v, accepted.report().voltage_v);
    assert_eq!(retry.report().cell_power_w, accepted.report().cell_power_w);

    // The final root would lie in support, but the first actual trial does not.
    // Adaptive numerical weights cannot clip or skip that physical source query.
    let law = response(
        ElectricalPropertyBasis::Resistivity,
        ConductorMagnitude::CurrentDensity,
        0.01,
        100.0,
        1.0,
    );
    let bounded = specimen(law.clone(), law, 0.2)
        .with_adaptive_relaxation(adaptation())
        .unwrap();
    assert!(matches!(
        bounded.solve_candidate(&ends(1.0), &[], limits(), controls(), policy(), &mut || false),
        Err(ConductionResponseError::Material(_)),
    ));
}
