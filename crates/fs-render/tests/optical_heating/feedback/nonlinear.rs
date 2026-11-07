use super::*;
use fs_render::volumes::heating::feedback::nonlinear::{
    TemperatureTwoPhotonTable, ThermalTwoPhotonCell, TwoPhotonFeedback,
};

fn beta(lo_k: f64, hi_k: f64, lo: f64, hi: f64, ceiling: f64) -> TemperatureTwoPhotonTable {
    TemperatureTwoPhotonTable::try_new(
        ContentHash([83; 32]),
        &[(lo_k, lo), (hi_k, hi)],
        (400.0, 800.0),
        ceiling,
        2,
    )
    .unwrap()
}

#[allow(clippy::too_many_arguments)]
fn nonlinear_heating(
    cx: &Cx<'_>,
    original: &DepositionSimulation<'_>,
    power: f64,
    area: f64,
    linear: &[TemperatureAbsorptionTable],
    nonlinear: &[TemperatureTwoPhotonTable],
    cells: &[ThermalOpticalCell],
    budget: OpticalFeedbackBudget,
) -> Result<TemperatureOpticalHeating<1>, OpticalFeedbackError> {
    let density = vec![1.0; cells.len()];
    let medium = grid(cx, &density, MassAbsorptionSpectrum::constant(1.0).unwrap());
    let bindings: Vec<_> = cells
        .iter()
        .enumerate()
        .map(|(i, cell)| ThermalTwoPhotonCell {
            grid_cell: cell.grid_cell,
            two_photon_table: i,
        })
        .collect();
    TemperatureOpticalHeating::try_new_two_photon(
        cx,
        &medium,
        beam(power, cells.len()),
        linear,
        cells,
        TwoPhotonFeedback {
            beam_area_m2: area,
            tables: nonlinear,
            cells: &bindings,
        },
        original,
        budget,
        1.0e-11,
        || Ok(()),
    )
}

#[test]
fn g1_temperature_two_photon_heating_converges_to_independent_nonlinear_solution() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let mut work = budget();
        work.max_source_evaluations = 4000;
        let optics = nonlinear_heating(
            cx,
            &original,
            4.0,
            1.0,
            &[table(300.0, 400.0, 0.0, 0.0)],
            &[beta(300.0, 400.0, 0.2, 10.2, 100.0)],
            &[binding(0, 11, 0)],
            work,
        )
        .unwrap();
        // alpha=0, beta=.2+.1(T-300), m*cp=4, P=4, L/A=1.
        // Integrating dT/dt=P/4 * (4*beta)/(1+4*beta) gives
        // ln(beta/.2)+4*(beta-.2)=.4*t. Solve that independent implicit
        // closed form, not the production attenuation/time-stepping equation.
        let (mut lo, mut hi) = (0.2_f64, 0.7_f64);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if (mid / 0.2).ln() + 4.0 * (mid - 0.2) < 2.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let exact = 300.0 + (0.5 * (lo + hi) - 0.2) / 0.1;
        let mut errors = Vec::new();
        for dt in [0.5, 0.25, 0.125] {
            let (next, accepted, audit) = original
                .advance_observed(&optics, 5.0, config(dt), || Ok(()))
                .unwrap();
            errors.push((next.temperature_k(CellId(11)).unwrap() - exact).abs());
            close(accepted.energy().thermal_input_j, audit.external_heat_j);
            close(
                accepted.energy().incident_j[0],
                accepted.energy().absorbed_j[0] + accepted.energy().transmitted_j[0],
            );
            assert_eq!(accepted.source_evaluations(), 2 * audit.substeps);
        }
        assert!((1.9..2.1).contains(&(errors[0] / errors[1])));
        assert!((1.9..2.1).contains(&(errors[1] / errors[2])));
        let controls = FeedbackRefinement {
            temperature_absolute_tolerance_k: 0.002,
            enthalpy_absolute_tolerance_j_kg: 0.008,
            optical_energy_absolute_tolerance_j: 0.008,
            ..refinement()
        };
        let (next, accepted, convergence) = optics
            .advance_refined(&original, 5.0, controls, || Ok(()))
            .unwrap();
        assert!((next.temperature_k(CellId(11)).unwrap() - exact).abs() < 0.0011);
        assert!(convergence.normalized_difference <= 1.0);
        assert!((0.9..1.1).contains(&convergence.observed_order.unwrap()));
        close(next.energy_j(), accepted.energy().absorbed_j[0]);
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(optics.energy().absorbed_j, [0.0]);
    });
}

#[test]
fn g1_two_photon_feedback_preserves_linear_weak_loss_area_and_opaque_limits() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let linear = [table(300.0, 400.0, 0.3, 0.3)];
        let cells = [binding(0, 11, 0)];
        let legacy = feedback(cx, &original, 4.0, &linear, &cells);
        let zero = nonlinear_heating(
            cx,
            &original,
            4.0,
            1.0,
            &linear,
            &[beta(300.0, 400.0, 0.0, 0.0, 100.0)],
            &cells,
            budget(),
        )
        .unwrap();
        assert_eq!(legacy.power().absorbed_w, zero.power().absorbed_w);
        assert_eq!(legacy.power().transmitted_w, zero.power().transmitted_w);
        let dark = nonlinear_heating(
            cx,
            &original,
            0.0,
            1.0,
            &[table(300.0, 400.0, 0.0, 0.0)],
            &[beta(300.0, 400.0, 0.5, 0.5, 100.0)],
            &cells,
            budget(),
        )
        .unwrap();
        assert_eq!(dark.power().absorbed_w, [0.0]);
        assert_eq!(dark.power().transmitted_w, [0.0]);
        let (unchanged, dark_after, _) = original
            .advance_observed(&dark, 1.0, config(1.0), || Ok(()))
            .unwrap();
        assert_eq!(unchanged.cells(), original.cells());
        assert_eq!(dark_after.energy().thermal_input_j, 0.0);
        for (beta_value, area) in [(1.0e-20, 1.0), (0.25, 1.0), (0.25, 2.0)] {
            let optics = nonlinear_heating(
                cx,
                &original,
                4.0,
                area,
                &[table(300.0, 400.0, 0.0, 0.0)],
                &[beta(300.0, 400.0, beta_value, beta_value, 100.0)],
                &cells,
                budget(),
            )
            .unwrap();
            let q = beta_value * 4.0 / area;
            let expected = 4.0 * q / (1.0 + q);
            assert!(optics.power().absorbed_w[0] > 0.0);
            assert!((optics.power().absorbed_w[0] / expected - 1.0).abs() < 3.0e-15);
            close(optics.power().transmitted_w[0], 4.0 / (1.0 + q));
        }
        let two = thermal(
            vec![cell(11, 1.0, 1.0, 0.0), cell(22, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let opaque = nonlinear_heating(
            cx,
            &two,
            4.0,
            1.0,
            &[
                table(300.0, 400.0, 1.0e308, 1.0e308),
                table(300.0, 400.0, 0.0, 0.0),
            ],
            &[
                beta(300.0, 400.0, 0.25, 0.25, 100.0),
                beta(300.0, 400.0, 0.25, 0.25, 100.0),
            ],
            &[binding(0, 11, 0), binding(1, 22, 1)],
            budget(),
        )
        .unwrap();
        assert_eq!(opaque.sources_w()[&CellId(11)], 4.0);
        assert_eq!(opaque.sources_w()[&CellId(22)], 0.0);
        assert_eq!(opaque.power().transmitted_w, [0.0]);
    });
}

#[test]
fn g3_local_beta_changes_nonlinear_shielding_on_original_cells() {
    with_cx(|cx| {
        let original = thermal(
            vec![cell(11, 1.0, 1.0, 0.0), cell(22, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let optics = nonlinear_heating(
            cx,
            &original,
            8.0,
            1.0,
            &[table(300.0, 400.0, 0.0, 0.0)],
            &[
                beta(300.0, 400.0, 0.025, 20.025, 100.0),
                beta(300.0, 400.0, 0.4, 0.4, 100.0),
            ],
            &[binding(0, 11, 0), binding(1, 22, 0)],
            budget(),
        )
        .unwrap();
        let (next, accepted, audit) = original
            .advance_observed(&optics, 3.0, config(0.01), || Ok(()))
            .unwrap();
        assert!(accepted.sources_w()[&CellId(11)] > optics.sources_w()[&CellId(11)] * 1.5);
        assert!(accepted.sources_w()[&CellId(22)] < optics.sources_w()[&CellId(22)]);
        close(next.energy_j(), accepted.energy().thermal_input_j);
        close(audit.external_heat_j, accepted.energy().absorbed_j[0]);
        close(
            accepted.energy().incident_j[0],
            accepted.energy().absorbed_j[0] + accepted.energy().transmitted_j[0],
        );
    });
}

#[test]
fn g4_final_endpoint_intensity_temperature_and_work_refuse_transactionally() {
    with_cx(|cx| {
        let original = thermal(
            vec![cell(11, 1.0, 1.0, 0.0), cell(22, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let optics = nonlinear_heating(
            cx,
            &original,
            8.0,
            1.0,
            &[table(300.0, 400.0, 0.0, 0.0)],
            &[
                beta(300.0, 302.0, 1.0, 0.0, 100.0),
                beta(300.0, 400.0, 0.0, 0.0, 1.0),
            ],
            &[binding(0, 11, 0), binding(1, 22, 0)],
            budget(),
        )
        .unwrap();
        // Initial shielding admits only 8/9 W/m² to the rear material.
        // The front heats and becomes clearer. Its accepted endpoint would
        // violate the rear ceiling even though no subsequent step is requested.
        assert!(matches!(
            original.advance_observed(&optics, 1.0, config(1.0), || Ok(())),
            Err(OpticalFeedbackError::IrradianceDomain { grid_cell: 1, .. })
        ));
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(optics.energy().thermal_input_j, 0.0);
        let one = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let narrow = nonlinear_heating(
            cx,
            &one,
            4.0,
            1.0,
            &[table(300.0, 400.0, 0.0, 0.0)],
            &[beta(300.0, 300.01, 0.5, 0.5, 100.0)],
            &[binding(0, 11, 0)],
            budget(),
        )
        .unwrap();
        assert!(matches!(
            one.advance_observed(&narrow, 1.0, config(1.0), || Ok(())),
            Err(OpticalFeedbackError::TemperatureDomain { .. })
        ));
        let mut limited = budget();
        limited.max_source_evaluations = 1;
        let bounded = nonlinear_heating(
            cx,
            &one,
            4.0,
            1.0,
            &[table(300.0, 400.0, 0.0, 0.0)],
            &[beta(300.0, 400.0, 0.5, 0.5, 100.0)],
            &[binding(0, 11, 0)],
            limited,
        )
        .unwrap();
        assert!(matches!(
            one.advance_observed(&bounded, 1.0, config(1.0), || Ok(())),
            Err(OpticalFeedbackError::Budget)
        ));
        assert!(matches!(
            one.advance_observed(&bounded, 1.0, config(1.0), || Err(
                DepositionError::Cancelled
            )),
            Err(OpticalFeedbackError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(one.time_s(), 0.0);
        assert_eq!(bounded.energy().absorbed_j, [0.0]);
    });
}

#[test]
fn g3_nonlinear_feedback_retains_native_latent_events_and_restart() {
    with_cx(|cx| {
        let original = thermal(
            vec![cell(11, 1.0, 1.0, 0.0)],
            &LATENT,
            vec![TimedEvent {
                time_s: 0.5,
                event: Event::Boundary(CellId(11), boundary(4.0)),
            }],
        );
        let optics = nonlinear_heating(
            cx,
            &original,
            40.0,
            1.0,
            &[table(300.0, 400.0, 0.0, 0.0)],
            &[beta(300.0, 400.0, 0.5, 1.5, 100.0)],
            &[binding(0, 11, 0)],
            budget(),
        )
        .unwrap();
        let (whole, accepted, audit) = original
            .advance_observed(&optics, 1.0, config(0.125), || Ok(()))
            .unwrap();
        let (middle, middle_optics, _) = original
            .advance_observed(&optics, 0.5, config(0.125), || Ok(()))
            .unwrap();
        let (split, split_optics, _) = middle
            .advance_observed(&middle_optics, 1.0, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(whole.temperature_k(CellId(11)).unwrap(), 310.0);
        assert_eq!(whole.cells(), split.cells());
        assert_eq!(accepted.energy(), split_optics.energy());
        close(
            audit.external_heat_j,
            accepted.energy().thermal_input_j + 2.0,
        );
        assert_eq!(audit.events, 1);
    });
}

#[test]
fn g0_nonlinear_sources_are_explicit_independent_and_jointly_bounded() {
    for rows in [
        vec![(300.0, -1.0), (400.0, 1.0)],
        vec![(300.0, 1.0), (300.0, 1.0)],
    ] {
        assert!(
            TemperatureTwoPhotonTable::try_new(
                ContentHash([83; 32]),
                &rows,
                (400.0, 800.0),
                100.0,
                2
            )
            .is_err()
        );
    }
    assert!(
        TemperatureTwoPhotonTable::try_new(
            ContentHash([0; 32]),
            &[(300.0, 1.0), (400.0, 1.0)],
            (400.0, 800.0),
            100.0,
            2
        )
        .is_err()
    );
    assert!(
        TemperatureTwoPhotonTable::try_new(
            ContentHash([83; 32]),
            &[(300.0, 1.0), (400.0, 1.0)],
            (400.0, 800.0),
            0.0,
            2
        )
        .is_err()
    );
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let linear = [table(300.0, 400.0, 0.0, 0.0)];
        let sources = [beta(300.0, 400.0, 0.5, 0.5, 100.0)];
        assert_eq!(sources[0].source_identity(), ContentHash([83; 32]));
        assert_eq!(sources[0].temperature_bounds_k(), (300.0, 400.0));
        for constrained in [
            OpticalFeedbackBudget {
                max_tables: 1,
                ..budget()
            },
            OpticalFeedbackBudget {
                max_total_temperature_knots: 3,
                ..budget()
            },
        ] {
            assert!(matches!(
                nonlinear_heating(
                    cx,
                    &original,
                    4.0,
                    1.0,
                    &linear,
                    &sources,
                    &[binding(0, 11, 0)],
                    constrained
                ),
                Err(OpticalFeedbackError::Budget)
            ));
        }
        assert!(
            nonlinear_heating(
                cx,
                &original,
                4.0,
                0.0,
                &linear,
                &sources,
                &[binding(0, 11, 0)],
                budget()
            )
            .is_err()
        );
        let off_band = TemperatureTwoPhotonTable::try_new(
            ContentHash([83; 32]),
            &[(300.0, 0.5), (400.0, 0.5)],
            (600.0, 800.0),
            100.0,
            2,
        )
        .unwrap();
        assert!(matches!(
            nonlinear_heating(
                cx,
                &original,
                4.0,
                1.0,
                &linear,
                &[off_band],
                &[binding(0, 11, 0)],
                budget()
            ),
            Err(OpticalFeedbackError::Transport(_))
        ));
    });
}
