//! Complete scan objectives retain the original contact-law control derivatives.
use super::*;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::{
    ScanObjectiveBudget, ScanObjectiveSampling, ScanObservable, ScanObservation,
    checkpoint::{ScanCheckpointBudget, evaluate_scan_objective_checkpointed_with_peaks},
    evaluate_scan_objective,
    responses::{
        ScanResponse, ScanResponseBudget, evaluate_scan_objective_with_responses,
        peaks::{ScanPeakResponse, evaluate_scan_objective_with_peaks},
    },
    value::{
        ScanValueBudget,
        responses::{ScanValueResponseBudget, evaluate_scan_values_with_responses},
    },
};
use fs_conduction::interface::InterfaceSurface;

const TIMES: [f64; 2] = [0.03, 0.08];
const CONTROLS: [f64; 7] = [0.15, 0.2, 0.3, 1.6, 1.5, 1.2, 3.];

fn observations() -> [ScanObservation; 2] {
    [
        ScanObservation {
            endpoint: 0,
            vertex: 3,
            observable: ScanObservable::Temperature,
            target: 309.,
            scale: 2.,
        },
        ScanObservation {
            endpoint: 1,
            vertex: 1,
            observable: ScanObservable::SpecificEnthalpy,
            target: 110.,
            scale: 8.,
        },
    ]
}

fn requests() -> ([ScanResponse; 1], [ScanPeakResponse; 1]) {
    (
        [ScanResponse {
            endpoint: 0,
            vertex: 0,
            observable: ScanObservable::Temperature,
        }],
        [ScanPeakResponse {
            first_endpoint: 0,
            last_endpoint: 1,
            vertices: vec![0, 1, 2, 3],
            observable: ScanObservable::SpecificEnthalpy,
            smoothing_scale: 2.,
        }],
    )
}

fn objective_budget(controls: usize) -> ScanObjectiveBudget {
    ScanObjectiveBudget {
        max_steps: 2,
        max_observations: 2,
        max_retained_vertices: 16,
        max_retained_source_vertices: 16,
        max_schedule_controls: controls,
        max_retained_radiation_entries: 2048,
    }
}

fn response_budget(controls: usize) -> ScanResponseBudget {
    ScanResponseBudget {
        max_responses: 2,
        max_gradient_entries: 2 * (8 + controls),
        max_transpose_solves: 5,
    }
}

fn checkpoints() -> ScanCheckpointBudget {
    ScanCheckpointBudget {
        segment_steps: 1,
        max_checkpoint_vertices: 24,
        max_forward_solves: 4,
    }
}

#[test]
fn g3_contact_objective_point_and_peak_gradients_match_fresh_whole_law_trajectories() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let history = control_history(cx, &mesh);
        let card = card();
        let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic(
            ConductivityTable::declared_curve(vec![(280., 1.6), (380., 3.6)]).unwrap(),
        );
        let (points, peaks) = requests();
        for radiative in [false, true] {
            let step = control_step(radiative);
            // Seven heater coordinates, one complete contact multiplier, and
            // optionally the existing reservoir temperature and emissivity.
            let controls = if radiative { 10 } else { 8 };
            let make =
                || control_scan(cx, &history, &material, &curves, &card, CONTROLS, radiative);
            let (objective, rows, envelopes) = evaluate_scan_objective_with_peaks(
                cx,
                make(),
                &TIMES,
                &observations(),
                ScanObjectiveSampling::RightEndpoint,
                step,
                control_linear(),
                objective_budget(controls),
                &points,
                &peaks,
                response_budget(controls),
                8,
            )
            .unwrap();
            let replay = evaluate_scan_objective_checkpointed_with_peaks(
                cx,
                make(),
                &TIMES,
                &observations(),
                ScanObjectiveSampling::RightEndpoint,
                step,
                control_linear(),
                ScanObjectiveBudget {
                    max_retained_vertices: 8,
                    max_retained_source_vertices: 8,
                    max_retained_radiation_entries: 1024,
                    ..objective_budget(controls)
                },
                &points,
                &peaks,
                response_budget(controls),
                8,
                checkpoints(),
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
            assert_eq!(objective.radiation, replay.objective.radiation);
            assert_eq!(
                objective.initial_specific_enthalpy,
                replay.objective.initial_specific_enthalpy
            );
            assert_eq!(
                objective.transpose_iterations,
                replay.objective.transpose_iterations
            );
            assert_eq!(
                rows[0].transpose_iterations,
                replay.responses[0].transpose_iterations
            );
            assert_eq!(
                envelopes[0].transpose_iterations,
                replay.peaks[0].transpose_iterations
            );
            assert_eq!(replay.work.original_forward_solves, 2);
            assert_eq!(replay.work.recomputed_forward_solves, 2);
            assert_eq!(replay.work.peak_retained_linearizations, 1);
            let native = coupled(&objective.final_step)
                .physical_contacts
                .frozen_interfaces();
            let gradients = [
                &objective.contact_resistance_multipliers,
                &rows[0].contact_resistance_multipliers,
                &envelopes[0].contact_resistance_multipliers,
            ];
            for gradient in gradients {
                assert_eq!(gradient.len(), 1);
                assert_eq!(gradient[0].interface, "joint");
                assert_eq!(gradient[0].card_identity, card.content_hash());
                assert_eq!(Some(gradient[0].mapped), native.surface_is_mapped("joint"));
                assert_eq!(gradient[0].face_pairs, 1);
            }
            assert!(rows[0].contact_resistance_multipliers[0].derivative.abs() > 1e-5);
            assert!(
                envelopes[0].contact_resistance_multipliers[0]
                    .derivative
                    .abs()
                    > 1e-5
            );
            let sample = |scale: f64| {
                // Scale R(T) and its slope together, then evolve both bodies
                // again through every original nonlinear contact/radiation solve.
                let trial_card = card_with_curve(vec![(280., 0.1 * scale), (380., 0.7 * scale)]);
                let run = control_scan(
                    cx,
                    &history,
                    &material,
                    &curves,
                    &trial_card,
                    CONTROLS,
                    radiative,
                );
                let (objective, rows, peaks) = evaluate_scan_values_with_responses(
                    cx,
                    run,
                    &TIMES,
                    &observations(),
                    ScanObjectiveSampling::RightEndpoint,
                    step,
                    ScanValueBudget {
                        max_steps: 2,
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
                sample(1.),
                [objective.value, rows[0].value, envelopes[0].value]
            );
            let epsilon = 1e-4;
            let (plus, minus) = (sample(1. + epsilon), sample(1. - epsilon));
            for (i, gradient) in gradients.into_iter().enumerate() {
                close_control(
                    gradient[0].derivative,
                    (plus[i] - minus[i]) / (2. * epsilon),
                );
            }
        }
    });
}

#[test]
fn g0_contact_multipliers_charge_exact_objective_and_response_scalar_allowances() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let history = control_history(cx, &mesh);
        let card = card();
        let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let (points, peaks) = requests();
        for radiative in [false, true] {
            let controls = if radiative { 10 } else { 8 };
            let step = control_step(radiative);
            let make =
                || control_scan(cx, &history, &material, &curves, &card, CONTROLS, radiative);
            let objective = objective_budget(controls);
            let missing_contact = objective_budget(controls - 1);
            assert!(
                evaluate_scan_objective(
                    cx,
                    make(),
                    &TIMES,
                    &observations(),
                    ScanObjectiveSampling::RightEndpoint,
                    step,
                    control_linear(),
                    missing_contact
                )
                .is_err()
            );
            assert!(
                evaluate_scan_objective_checkpointed_with_peaks(
                    cx,
                    make(),
                    &TIMES,
                    &observations(),
                    ScanObjectiveSampling::RightEndpoint,
                    step,
                    control_linear(),
                    missing_contact,
                    &points,
                    &peaks,
                    response_budget(controls),
                    8,
                    checkpoints()
                )
                .is_err()
            );
            let point_budget = ScanResponseBudget {
                max_responses: 1,
                max_gradient_entries: 8 + controls,
                max_transpose_solves: 3,
            };
            assert!(
                evaluate_scan_objective_with_responses(
                    cx,
                    make(),
                    &TIMES,
                    &observations(),
                    ScanObjectiveSampling::RightEndpoint,
                    step,
                    control_linear(),
                    objective,
                    &points,
                    ScanResponseBudget {
                        max_gradient_entries: 7 + controls,
                        ..point_budget
                    }
                )
                .is_err()
            );
            assert!(
                evaluate_scan_objective_with_responses(
                    cx,
                    make(),
                    &TIMES,
                    &observations(),
                    ScanObjectiveSampling::RightEndpoint,
                    step,
                    control_linear(),
                    objective,
                    &points,
                    point_budget
                )
                .is_ok()
            );
            let short = ScanResponseBudget {
                max_gradient_entries: 2 * (8 + controls) - 1,
                ..response_budget(controls)
            };
            assert!(
                evaluate_scan_objective_with_peaks(
                    cx,
                    make(),
                    &TIMES,
                    &observations(),
                    ScanObjectiveSampling::RightEndpoint,
                    step,
                    control_linear(),
                    objective,
                    &points,
                    &peaks,
                    short,
                    8
                )
                .is_err()
            );
            assert!(
                evaluate_scan_objective_checkpointed_with_peaks(
                    cx,
                    make(),
                    &TIMES,
                    &observations(),
                    ScanObjectiveSampling::RightEndpoint,
                    step,
                    control_linear(),
                    objective,
                    &points,
                    &peaks,
                    short,
                    8,
                    checkpoints()
                )
                .is_err()
            );
        }
    });
}

#[test]
fn g0_unattached_contact_models_return_no_invented_resistance_controls() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let history = control_history(cx, &mesh);
        let card = card();
        let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let initial = history.at(cx, 0.).unwrap();
        // The original fixed-interface path remains an explicit fixed design
        // input. It must not silently acquire a live contact control model.
        let fixed = InterfaceSurface::new(
            "joint",
            ThermalInterfaces::coincident_face_pairs(&mesh).unwrap(),
            law(&mesh, &card).reference_resistance().clone(),
        )
        .unwrap();
        let interfaces = ThermalInterfaces::new(&mesh, &initial, vec![fixed]).unwrap();
        let base = EnthalpyHistoryRun::new(
            cx,
            ThermalHistoryProblem {
                boundary: &history,
                material: &material,
                element_materials: None,
                interfaces: Some(&interfaces),
                source: HeatSourceHistory::Fixed(ScalarField::Uniform(2.)),
            },
            EnthalpyHistoryMaterials {
                table: &[
                    ReferenceEnthalpyMaterial {
                        curve: &curves[0],
                        reference_density_kg_m3: 1.,
                    },
                    ReferenceEnthalpyMaterial {
                        curve: &curves[1],
                        reference_density_kg_m3: 2.,
                    },
                ],
                of_element: &[0, 1],
                phase_conductivity: None,
            },
            &OLD,
            EnthalpyHistoryBudget {
                spatial: EnthalpyBudget {
                    max_vertices: 8,
                    max_elements: 2,
                },
                max_steps: 2,
            },
        )
        .unwrap();
        let run = ScanningEnthalpyRun::new(
            cx,
            base,
            vec![heater(
                constant(0.3, LENGTH),
                constant(3., COMPONENT_POWER_DIMS),
            )],
            1,
            None,
        )
        .unwrap();
        assert!(run.contacts().is_none());
        let (points, peaks) = requests();
        let (objective, rows, envelopes) = evaluate_scan_objective_with_peaks(
            cx,
            run,
            &[0.001, 0.002],
            &observations(),
            ScanObjectiveSampling::RightEndpoint,
            config(),
            control_linear(),
            objective_budget(7),
            &points,
            &peaks,
            response_budget(7),
            8,
        )
        .unwrap();
        assert!(objective.contact_resistance_multipliers.is_empty());
        assert!(rows[0].contact_resistance_multipliers.is_empty());
        assert!(envelopes[0].contact_resistance_multipliers.is_empty());
    });
}
