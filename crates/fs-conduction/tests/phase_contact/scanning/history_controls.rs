//! Shared contact-law controls through objectives, limits and checkpoint replay.
use super::*;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::{
    ScanObjectiveBudget, ScanObjectiveSampling, ScanObservable, ScanObservation,
    checkpoint::{ScanCheckpointBudget, evaluate_scan_objective_checkpointed_with_peaks},
    responses::{
        ScanResponse, ScanResponseBudget,
        peaks::{ScanPeakResponse, evaluate_scan_objective_with_peaks},
    },
    value::{
        ScanValueBudget,
        responses::{ScanValueResponseBudget, evaluate_scan_values_with_responses},
    },
};

const TIMES: [f64; 4] = [0.025, 0.05, 0.075, 0.1];

fn controls() -> ScanObjectiveBudget {
    ScanObjectiveBudget {
        max_steps: 4,
        max_observations: 3,
        max_retained_vertices: 32,
        max_retained_source_vertices: 32,
        // x/y constants, three z table values, three widths, one power,
        // and ONE original interface multiplier, even for delegated traces.
        max_schedule_controls: 10,
        max_retained_radiation_entries: 0,
    }
}

fn limits() -> ScanResponseBudget {
    ScanResponseBudget {
        max_responses: 2,
        max_gradient_entries: 2 * (8 + 10),
        max_transpose_solves: 4 + 2 + 4,
    }
}

fn observations() -> [ScanObservation; 3] {
    [
        ScanObservation {
            endpoint: 0,
            vertex: 0,
            observable: ScanObservable::SpecificEnthalpy,
            target: 30.0,
            scale: 4.0,
        },
        ScanObservation {
            endpoint: 2,
            vertex: 5,
            observable: ScanObservable::LiquidMassFraction,
            target: 0.5,
            scale: 0.1,
        },
        ScanObservation {
            endpoint: 3,
            vertex: 1,
            observable: ScanObservable::SpecificEnthalpy,
            target: 40.0,
            scale: 5.0,
        },
    ]
}

fn responses() -> ([ScanResponse; 1], [ScanPeakResponse; 1]) {
    (
        [ScanResponse {
            endpoint: 1,
            vertex: 0,
            observable: ScanObservable::SpecificEnthalpy,
        }],
        [ScanPeakResponse {
            first_endpoint: 0,
            last_endpoint: 3,
            vertices: vec![0, 1, 4, 5],
            observable: ScanObservable::LiquidMassFraction,
            smoothing_scale: 0.02,
        }],
    )
}

fn path() -> TimeSignal {
    TimeSignal::Table {
        times: vec![0.0, 0.05, 0.1],
        values: vec![-0.2, 0.1, 0.2],
        dims: LENGTH,
        interp: Interp::Hold,
    }
}

#[test]
fn g3_contact_law_controls_survive_whole_history_limits_and_replay() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = scan_history(cx, &mesh);
        let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [common::bulk("cold"), common::bulk("hot")];
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source_card = source(&cards, 1.0, false, false, 1.0);
        let query_point = point();
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let (points, peaks) = responses();
        let transpose = LinearConfig {
            tolerance: 1e-11,
            max_iterations: 128,
            restart: 8,
        };
        for mode in 0..3 {
            for sampling in [
                ScanObjectiveSampling::EventLeft,
                ScanObjectiveSampling::RightEndpoint,
            ] {
                let make = || {
                    let laws = contacts(
                        cx,
                        &mesh,
                        &initial,
                        &cards,
                        &curves,
                        &source_card,
                        &query_point,
                        sides(&mesh, 4, 2, 0.0),
                        mode,
                        false,
                    );
                    assert_eq!(laws.surface_count(), 1);
                    ScanningEnthalpyRun::new(
                        cx,
                        scan_base(cx, &boundary, &conductivity, &curves, 4),
                        vec![heater(path())],
                        1,
                        None,
                    )
                    .unwrap()
                    .with_contacts_newton(cx, laws, config())
                    .unwrap()
                };
                let (objective, rows, envelopes) = evaluate_scan_objective_with_peaks(
                    cx,
                    make(),
                    &TIMES,
                    &observations(),
                    sampling,
                    config().step,
                    transpose,
                    controls(),
                    &points,
                    &peaks,
                    limits(),
                    16,
                )
                .unwrap();
                let mut short_tape = controls();
                short_tape.max_retained_vertices = 16;
                short_tape.max_retained_source_vertices = 16;
                let replay = evaluate_scan_objective_checkpointed_with_peaks(
                    cx,
                    make(),
                    &TIMES,
                    &observations(),
                    sampling,
                    config().step,
                    transpose,
                    short_tape,
                    &points,
                    &peaks,
                    limits(),
                    16,
                    ScanCheckpointBudget {
                        segment_steps: 2,
                        max_checkpoint_vertices: 24,
                        max_forward_solves: 8,
                    },
                )
                .unwrap();
                assert_eq!(
                    objective.contact_resistance_multipliers,
                    replay.objective.contact_resistance_multipliers
                );
                assert_eq!(
                    rows[0].contact_resistance_multipliers,
                    replay.responses[0].contact_resistance_multipliers
                );
                assert_eq!(
                    envelopes[0].contact_resistance_multipliers,
                    replay.peaks[0].contact_resistance_multipliers
                );
                assert_eq!(objective.sources, replay.objective.sources);
                assert_eq!(
                    objective.initial_specific_enthalpy,
                    replay.objective.initial_specific_enthalpy
                );
                assert_eq!(replay.work.original_forward_solves, 4);
                assert_eq!(replay.work.recomputed_forward_solves, 4);
                let controls = [
                    &objective.contact_resistance_multipliers,
                    &rows[0].contact_resistance_multipliers,
                    &envelopes[0].contact_resistance_multipliers,
                ];
                for gradient in controls {
                    assert_eq!(gradient.len(), 1);
                    assert_eq!(gradient[0].interface, "joint");
                    assert_eq!(gradient[0].card_identity, source_card.content_hash());
                    assert!(gradient[0].face_pairs > 0);
                    assert!(gradient[0].derivative.is_finite());
                }
                // Perturb the ENTIRE original law in an independent value-only
                // trajectory, never a frozen endpoint resistance or its tangent.
                let values = |scale| {
                    let card = source(&cards, scale, false, false, 1.0);
                    let laws = contacts(
                        cx,
                        &mesh,
                        &initial,
                        &cards,
                        &curves,
                        &card,
                        &query_point,
                        sides(&mesh, 4, 2, 0.0),
                        mode,
                        false,
                    );
                    let run = ScanningEnthalpyRun::new(
                        cx,
                        scan_base(cx, &boundary, &conductivity, &curves, 4),
                        vec![heater(path())],
                        1,
                        None,
                    )
                    .unwrap()
                    .with_contacts_newton(cx, laws, config())
                    .unwrap();
                    let (objective, rows, peaks) = evaluate_scan_values_with_responses(
                        cx,
                        run,
                        &TIMES,
                        &observations(),
                        sampling,
                        config().step,
                        ScanValueBudget {
                            max_steps: 4,
                            max_observations: 3,
                            max_schedule_controls: 9,
                        },
                        &points,
                        &peaks,
                        ScanValueResponseBudget {
                            max_responses: 2,
                            max_peak_samples: 16,
                        },
                    )
                    .unwrap();
                    [objective.value, rows[0].value, peaks[0].value]
                };
                assert_eq!(
                    values(1.0),
                    [objective.value, rows[0].value, envelopes[0].value]
                );
                let epsilon = 1e-4;
                let (plus, minus) = (values(1.0 + epsilon), values(1.0 - epsilon));
                for (i, gradient) in controls.into_iter().enumerate() {
                    let difference = (plus[i] - minus[i]) / (2.0 * epsilon);
                    near(
                        gradient[0].derivative,
                        difference,
                        2e-5 * (1.0 + difference.abs()),
                    );
                }
                assert!(rows[0].contact_resistance_multipliers[0].derivative.abs() > 1e-3);
                if sampling == ScanObjectiveSampling::EventLeft {
                    assert_eq!(objective.sources[0].center_values_m[2][2], 0.0);
                    assert_eq!(rows[0].sources[0].center_values_m[2][1..], [0.0, 0.0]);
                }
            }
        }
    });
}

#[test]
fn g4_contact_controls_count_once_in_objective_and_response_admission() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = scan_history(cx, &mesh);
        let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [common::bulk("cold"), common::bulk("hot")];
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source_card = source(&cards, 1.0, false, false, 1.0);
        let query_point = point();
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let make = || {
            let laws = contacts(
                cx,
                &mesh,
                &initial,
                &cards,
                &curves,
                &source_card,
                &query_point,
                sides(&mesh, 4, 2, 0.0),
                2,
                false,
            );
            ScanningEnthalpyRun::new(
                cx,
                scan_base(cx, &boundary, &conductivity, &curves, 4),
                vec![heater(path())],
                1,
                None,
            )
            .unwrap()
            .with_contacts_newton(cx, laws, config())
            .unwrap()
        };
        let (points, peaks) = responses();
        let transpose = LinearConfig {
            tolerance: 1e-11,
            max_iterations: 128,
            restart: 8,
        };
        let mut missing_control = controls();
        missing_control.max_schedule_controls -= 1;
        let mut missing_row_entry = limits();
        missing_row_entry.max_gradient_entries -= 1;
        for (budget, rows) in [(missing_control, limits()), (controls(), missing_row_entry)] {
            let error = evaluate_scan_objective_with_peaks(
                cx,
                make(),
                &TIMES,
                &observations(),
                ScanObjectiveSampling::EventLeft,
                config().step,
                transpose,
                budget,
                &points,
                &peaks,
                rows,
                16,
            )
            .unwrap_err();
            assert!(error.to_string().contains("allowance"), "{error}");
        }
    });
}
