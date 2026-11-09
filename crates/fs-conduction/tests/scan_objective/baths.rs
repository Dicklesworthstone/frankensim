//! Prescribed bath histories through the same scan, response and replay adjoints.
use super::*;
use fs_conduction::bc::history::ThermalRegionHistory;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::{
    checkpoint::{ScanCheckpointBudget, evaluate_scan_objective_checkpointed_with_peaks},
    responses::{
        ScanResponse, ScanResponseBudget,
        peaks::{ScanPeakResponse, evaluate_scan_objective_with_peaks},
    },
    value::{
        ScanValueBudget, evaluate_scan_values,
        responses::{ScanValueResponseBudget, evaluate_scan_values_with_responses},
    },
};
use fs_scenario::bc::{BcKind, BcValue, BoundaryCondition, Physics};
use fs_scenario::signal::LinearChirp;

fn row(region: &str, value: BcValue) -> ThermalRegionHistory {
    ThermalRegionHistory {
        row: BoundaryCondition {
            region: region.into(),
            physics: Physics::Thermal,
            kind: BcKind::Dirichlet,
            value: Some(value),
            compatibility: None,
            frame: 0,
        },
        reference_temperature: None,
    }
}

fn scheduled<'m>(cx: &Cx<'_>, m: &'m ConductionMesh, value: BcValue) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(
        cx,
        ThermalBoundaryBuilder::new(m)
            .region(
                "bath",
                |face| face.centroid[2] == 0.,
                ThermalBc::dirichlet(299.).unwrap(),
            )
            .unwrap()
            .remainder("surface", ThermalBc::adiabatic())
            .unwrap(),
        &Environment::earth_lab(),
        vec![row("bath", value)],
        0.,
    )
    .unwrap()
}

fn values(kind: usize, values: &[f64]) -> BcValue {
    if kind == 0 {
        return BcValue::Uniform(QtyAny::new(values[0], TEMPERATURE_DIMS));
    }
    BcValue::Signal(match kind {
        1 => constant(values[0], TEMPERATURE_DIMS),
        2 => TimeSignal::Ramp {
            t_start: 0.,
            t_end: 0.006,
            from: QtyAny::new(values[0], TEMPERATURE_DIMS),
            to: QtyAny::new(values[1], TEMPERATURE_DIMS),
        },
        3 | 4 => TimeSignal::Table {
            times: vec![0., 0.002, 0.006],
            values: values.to_vec(),
            dims: TEMPERATURE_DIMS,
            interp: if kind == 3 {
                Interp::Linear
            } else {
                Interp::Hold
            },
        },
        _ => unreachable!("fixture family"),
    })
}

fn response_budget() -> ScanResponseBudget {
    ScanResponseBudget {
        max_responses: 2,
        max_gradient_entries: 128,
        max_transpose_solves: 10,
    }
}
fn checkpoints() -> ScanCheckpointBudget {
    ScanCheckpointBudget {
        segment_steps: 2,
        max_checkpoint_vertices: 12,
        max_forward_solves: 8,
    }
}
fn requests() -> ([ScanResponse; 1], [ScanPeakResponse; 1]) {
    (
        [ScanResponse {
            endpoint: 1,
            vertex: 3,
            observable: ScanObservable::Temperature,
        }],
        [ScanPeakResponse {
            first_endpoint: 0,
            last_endpoint: 3,
            vertices: vec![0, 3],
            observable: ScanObservable::Temperature,
            smoothing_scale: 0.5,
        }],
    )
}
fn measurements() -> [ScanObservation; 2] {
    [
        observation(0, 1, ScanObservable::Temperature, 300., 2.),
        observation(3, 3, ScanObservable::SpecificEnthalpy, 42., 4.),
    ]
}

#[test]
fn g3_authored_baths_match_physical_differences_and_checkpointed_responses() {
    let m = mesh();
    let curve = chart();
    let material = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(280., 1.), (380., 3.)]).unwrap(),
    );
    with_cx(|cx| {
        let old = [10., 20., 30., 40.];
        let p = [0.2, 0.25, 0.3, 1.4, 1.3, 1.2, 3.];
        let times = [0.001, 0.002, 0.004, 0.006];
        let (points, peaks) = requests();
        for (kind, sampling) in [
            (0, ScanObjectiveSampling::RightEndpoint),
            (1, ScanObjectiveSampling::EventLeft),
            (2, ScanObjectiveSampling::EventLeft),
            (3, ScanObjectiveSampling::EventLeft),
            (4, ScanObjectiveSampling::EventLeft),
            (4, ScanObjectiveSampling::RightEndpoint),
        ] {
            let controls = match kind {
                0 | 1 => vec![299.],
                2 => vec![297., 301.],
                _ => vec![297., 299., 301.],
            };
            let hist = scheduled(cx, &m, values(kind, &controls));
            let make = || run(cx, &hist, &material, &curve, &old, vec![heater(p)], true);
            let (objective, rows, envelopes) = evaluate_scan_objective_with_peaks(
                cx,
                make(),
                &times,
                &measurements(),
                sampling,
                config(),
                linear(),
                objective_budget(),
                &points,
                &peaks,
                response_budget(),
                8,
            )
            .unwrap();
            let replay = evaluate_scan_objective_checkpointed_with_peaks(
                cx,
                make(),
                &times,
                &measurements(),
                sampling,
                config(),
                linear(),
                objective_budget(),
                &points,
                &peaks,
                response_budget(),
                8,
                checkpoints(),
            )
            .unwrap();
            assert_eq!(
                objective.prescribed_temperatures,
                replay.objective.prescribed_temperatures
            );
            assert_eq!(
                rows[0].prescribed_temperatures,
                replay.responses[0].prescribed_temperatures
            );
            assert_eq!(
                envelopes[0].prescribed_temperatures,
                replay.peaks[0].prescribed_temperatures
            );
            assert_eq!(objective.sources, replay.objective.sources);
            assert_eq!(objective.radiation, replay.objective.radiation);
            assert_eq!(objective.prescribed_temperatures[0].region, "bath");
            let analytic = [
                &objective.prescribed_temperatures[0].values_k,
                &rows[0].prescribed_temperatures[0].values_k,
                &envelopes[0].prescribed_temperatures[0].values_k,
            ];
            let evaluate = |controls: &[f64]| {
                let hist = scheduled(cx, &m, values(kind, controls));
                let (objective, rows, peaks) = evaluate_scan_values_with_responses(
                    cx,
                    run(cx, &hist, &material, &curve, &old, vec![heater(p)], true),
                    &times,
                    &measurements(),
                    sampling,
                    config(),
                    ScanValueBudget {
                        max_steps: 4,
                        max_observations: 2,
                        max_schedule_controls: 7,
                    },
                    &points,
                    &peaks,
                    ScanValueResponseBudget {
                        max_responses: 2,
                        max_peak_samples: 8,
                    },
                )
                .unwrap();
                [objective.value, rows[0].value, peaks[0].value]
            };
            assert_eq!(
                evaluate(&controls),
                [objective.value, rows[0].value, envelopes[0].value]
            );
            for coordinate in 0..controls.len() {
                let (mut plus, mut minus) = (controls.clone(), controls.clone());
                plus[coordinate] += 1e-3;
                minus[coordinate] -= 1e-3;
                let (a, b) = (evaluate(&plus), evaluate(&minus));
                for row in 0..3 {
                    close(analytic[row][coordinate], (a[row] - b[row]) / 2e-3);
                }
            }
            if kind == 4 && sampling == ScanObjectiveSampling::EventLeft {
                assert_eq!(objective.prescribed_temperatures[0].values_k[2], 0.);
                assert_eq!(rows[0].prescribed_temperatures[0].values_k[1..], [0., 0.]);
                assert_eq!(envelopes[0].prescribed_temperatures[0].values_k[2], 0.);
            }
        }
    });
}

#[test]
fn g1_shared_vertices_use_original_primal_region_precedence_once() {
    let m = mesh();
    let curve = chart();
    let material = ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx| {
        for first_scheduled in [false, true] {
            // The two faces share vertices 0 and 1. The first partition region
            // owns them even though its HISTORY is declared after the second.
            let mut histories = vec![row(
                "second",
                BcValue::Uniform(QtyAny::new(301., TEMPERATURE_DIMS)),
            )];
            if first_scheduled {
                histories.push(row(
                    "first",
                    BcValue::Uniform(QtyAny::new(299., TEMPERATURE_DIMS)),
                ));
            }
            let hist = ThermalBoundaryHistory::new(
                cx,
                ThermalBoundaryBuilder::new(&m)
                    .region(
                        "first",
                        |f| f.centroid[2] == 0.,
                        ThermalBc::dirichlet(299.).unwrap(),
                    )
                    .unwrap()
                    .region(
                        "second",
                        |f| f.centroid[1] == 0.,
                        ThermalBc::dirichlet(301.).unwrap(),
                    )
                    .unwrap()
                    .adiabatic_remainder(),
                &Environment::earth_lab(),
                histories,
                0.,
            )
            .unwrap();
            let observed: Vec<_> = (0..4)
                .map(|v| observation(0, v, ScanObservable::Temperature, 298., 1.))
                .collect();
            let result = evaluate_scan_objective(
                cx,
                run(
                    cx,
                    &hist,
                    &material,
                    &curve,
                    &[0.; 4],
                    vec![heater([0.2, 0.25, 0.3, 1.4, 1.3, 1.2, 0.])],
                    false,
                ),
                &[0.001],
                &observed,
                ScanObjectiveSampling::RightEndpoint,
                config(),
                linear(),
                objective_budget(),
            )
            .unwrap();
            assert_eq!(result.predictions, [299., 299., 299., 301.]);
            assert_eq!(
                result.prescribed_temperatures.len(),
                if first_scheduled { 2 } else { 1 }
            );
            assert_eq!(result.prescribed_temperatures[0].region, "second");
            close(result.prescribed_temperatures[0].values_k[0], 301. - 298.);
            if first_scheduled {
                assert_eq!(result.prescribed_temperatures[1].region, "first");
                close(
                    result.prescribed_temperatures[1].values_k[0],
                    3. * (299. - 298.),
                );
            }
            assert_eq!(result.initial_specific_enthalpy, [0.; 4]);
        }
    });
}

#[test]
fn g0_bath_controls_charge_returned_gradient_shapes_but_do_not_restrict_value_only_signals() {
    let m = mesh();
    let curve = chart();
    let material = ConductivityModel::isotropic_declared(2.).unwrap();
    with_cx(|cx| {
        let hist = scheduled(cx, &m, values(4, &[297., 299., 301.]));
        let p = [0.2, 0.25, 0.3, 1.4, 1.3, 1.2, 3.];
        let times = [0.001, 0.002, 0.004, 0.006];
        let (points, peaks) = requests();
        let make = || {
            run(
                cx,
                &hist,
                &material,
                &curve,
                &[10., 20., 30., 40.],
                vec![heater(p)],
                true,
            )
        };
        let call = |budget, rows| {
            evaluate_scan_objective_with_peaks(
                cx,
                make(),
                &times,
                &measurements(),
                ScanObjectiveSampling::EventLeft,
                config(),
                linear(),
                budget,
                &points,
                &peaks,
                rows,
                8,
            )
        };
        // Seven heater values + two radiation values + three bath values.
        let budget = ScanObjectiveBudget {
            max_schedule_controls: 12,
            ..objective_budget()
        };
        let rows = ScanResponseBudget {
            max_gradient_entries: 32,
            ..response_budget()
        };
        assert!(
            call(
                ScanObjectiveBudget {
                    max_schedule_controls: 11,
                    ..budget
                },
                rows
            )
            .is_err()
        );
        assert!(
            call(
                budget,
                ScanResponseBudget {
                    max_gradient_entries: 31,
                    ..rows
                }
            )
            .is_err()
        );
        assert!(call(budget, rows).is_ok());
        assert!(
            evaluate_scan_objective_checkpointed_with_peaks(
                cx,
                make(),
                &times,
                &measurements(),
                ScanObjectiveSampling::EventLeft,
                config(),
                linear(),
                budget,
                &points,
                &peaks,
                ScanResponseBudget {
                    max_gradient_entries: 31,
                    ..rows
                },
                8,
                checkpoints()
            )
            .is_err()
        );
        let chirp = TimeSignal::Chirp(
            LinearChirp::new(
                0.,
                0.006,
                1.,
                2.,
                0.,
                QtyAny::new(299., TEMPERATURE_DIMS),
                QtyAny::new(1., TEMPERATURE_DIMS),
            )
            .unwrap(),
        );
        let chirped = scheduled(cx, &m, BcValue::Signal(chirp));
        let make = || {
            run(
                cx,
                &chirped,
                &material,
                &curve,
                &[10., 20., 30., 40.],
                vec![heater(p)],
                true,
            )
        };
        let observed = [observation(0, 3, ScanObservable::Temperature, 300., 1.)];
        assert!(
            evaluate_scan_values(
                cx,
                make(),
                &[0.001],
                &observed,
                ScanObjectiveSampling::RightEndpoint,
                config(),
                ScanValueBudget {
                    max_steps: 1,
                    max_observations: 1,
                    max_schedule_controls: 7
                }
            )
            .is_ok()
        );
        assert!(
            evaluate_scan_objective(
                cx,
                make(),
                &[0.001],
                &observed,
                ScanObjectiveSampling::RightEndpoint,
                config(),
                linear(),
                objective_budget()
            )
            .is_err()
        );
    });
}
