//! Stage-resolved material emission, independent baths, and retained accepted history.
use super::*;
use fs_conduction::transient::enthalpy::sdirk::{EnthalpySdirk2Stage, EnthalpySdirk2Step, SDIRK2_GAMMA as G};

fn heat(h: f64, bath: f64) -> f64 {
    0.5 * STEFAN_BOLTZMANN_W_M2_K4
        * (400.0_f64.powi(4) - bath.powi(4)) * (0.2 + 0.6 * h / LATENT)
}

// Exact independent scalar stage equations on three one-kilogram face nodes.
fn stages(old: f64, dt: f64, baths: [f64; 2]) -> [f64; 2] {
    let first = endpoint(old, G * dt, baths[0]);
    let carry = -(1.0 - G) * dt * heat(first, baths[0]) / 3.0;
    [first, endpoint(old + carry, G * dt, baths[1])]
}

fn admit_stage(stage: &EnthalpySdirk2Stage, body: usize, expected: f64, datum: f64, bath: f64) {
    assert!(stage.temperature.iter().all(|&temperature| temperature == 400.0));
    for local in 0..3 {
        let vertex = 4 * body + local;
        close(stage.specific_enthalpy_j_kg[vertex] - datum, expected);
        close(stage.liquid_mass_fraction[vertex], expected / LATENT);
    }
    close(stage.specific_enthalpy_j_kg[4 * body + 3] - datum, 500.0);
    let radiation = stage.pointwise_radiation.as_ref().unwrap();
    let patch = &radiation.patches[body];
    close(patch.patch.ambient_temperature_k(), bath);
    close(patch.outward_heat_w, heat(expected, bath));
    assert_eq!(patch.material_states.len(), 3);
    for state in &patch.material_states {
        close(state.emissivity.value(), 0.2 + 0.6 * expected / LATENT);
    }
    assert!(stage.stage_energy_residual_j.abs() <= config().energy_tolerance_j);
}

fn admit_single(step: &EnthalpySdirk2Step, old: f64, dt: f64, datum: f64, baths: [f64; 2]) -> f64 {
    let expected = stages(old, dt, baths);
    admit_stage(&step.first_stage, 0, expected[0], datum, baths[0]);
    admit_stage(&step.second_stage, 0, expected[1], datum, baths[1]);
    let outward = dt * ((1.0 - G) * heat(expected[0], baths[0]) + G * heat(expected[1], baths[1]));
    close(step.radiation_out_j, outward);
    close(step.stored_energy_change_j, 3.0 * (expected[1] - old));
    close(step.stored_energy_change_j + step.radiation_out_j, 0.0);
    assert_eq!(step.source_j, 0.0);
    assert_eq!(step.robin_out_j, 0.0);
    assert!(step.energy_residual_j.abs() <= config().energy_tolerance_j);
    expected[1]
}

#[test]
fn advanced_material_history_retains_phase_mass_and_claims_when_switching_to_sdirk() {
    let mesh = mesh(false);
    let card = source_card(1.0);
    let context = context();
    let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    for datum in [0.0, 4096.0] {
        let curve = phase_chart(&card, datum);
        with_cx(|cx| {
            let boundary = face_history(cx, &mesh);
            let signal = TimeSignal::Ramp { t_start: 0.0, t_end: 1.0,
                from: QtyAny::new(400.0, TEMPERATURE_DIMS), to: QtyAny::new(500.0, TEMPERATURE_DIMS) };
            let declaration = declaration(&card, &curve, &context, signal);
            let mut original = MaterialRadiativeEnthalpyHistoryRun::new(cx,
                base(cx, &boundary, &conductivity, &curve, datum + 500.0, 3),
                vec![declaration.clone()]).unwrap();
            let selected = original.reference_emissivities()[0].receipt().selected;
            let first = original.advance_to(cx, 0.05, config()).unwrap();
            let saved = original.conduction().specific_enthalpy_j_kg().to_vec();
            let old = first.specific_enthalpy_j_kg[0] - datum;
            let mut run = original.into_sdirk2();
            assert_eq!(run.conduction().specific_enthalpy_j_kg(), saved);
            assert_eq!(run.conduction().time_s(), 0.05);
            assert_eq!(run.conduction().accepted_steps(), 1);
            assert_eq!(run.conduction().storage().reference_nodal_masses_kg(), &[1.0; 4]);
            let retained = &run.material_ambient_patches().unwrap()[0];
            assert!(std::ptr::eq(retained.curve, &curve));
            assert!(std::ptr::eq(retained.material.context, &context));
            assert_eq!(retained.ambient_temperature, declaration.ambient_temperature);
            assert!(run.material_enclosure_patches().is_none());

            let end = 0.25;
            let dt = end - 0.05;
            let baths = [400.0 + 100.0 * G.mul_add(dt, 0.05), 425.0];
            let step = run.advance_to(cx, end, config()).unwrap();
            admit_single(&step, old, dt, datum, baths);
            for stage in [&step.first_stage, &step.second_stage] {
                for state in &stage.pointwise_radiation.as_ref().unwrap().patches[0].material_states {
                    assert_eq!(state.emissivity.receipt().selected, selected);
                    assert_eq!(state.emissivity.card_identity(), card.content_hash());
                }
            }
            assert_ne!(step.first_stage.specific_enthalpy_j_kg[0],
                step.second_stage.specific_enthalpy_j_kg[0]);
            assert_eq!(run.conduction().time_s(), end);
            assert_eq!(run.conduction().accepted_steps(), 2);
        });
    }
}

#[test]
fn independent_bath_events_drive_original_fine_stages_and_one_entropy_window() {
    let mesh = mesh(true);
    let card = source_card(1.0);
    let context = context();
    let curves = [phase_chart(&card, 0.0), phase_chart(&card, 4096.0)];
    let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    with_cx(|cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
            .region("left", |face| face.centroid[2] == 0.0 && face.centroid[0] < 2.0,
                ThermalBc::adiabatic()).unwrap()
            .region("right", |face| face.centroid[2] == 0.0 && face.centroid[0] > 2.0,
                ThermalBc::adiabatic()).unwrap().adiabatic_remainder(),
            &Environment::earth_lab(), Vec::new(), 0.0).unwrap();
        let original = EnthalpyHistoryRun::new(cx, problem(&boundary, &conductivity,
            HeatSourceHistory::Fixed(ScalarField::Uniform(0.0))), EnthalpyHistoryMaterials {
                table: &[
                    ReferenceEnthalpyMaterial { curve: &curves[0], reference_density_kg_m3: 24.0 },
                    ReferenceEnthalpyMaterial { curve: &curves[1], reference_density_kg_m3: 24.0 },
                ], of_element: &[0, 1], phase_conductivity: None,
            }, &[500.0, 500.0, 500.0, 500.0, 4596.0, 4596.0, 4596.0, 4596.0],
            budget(&mesh, 8)).unwrap();
        let schedules = [
            TimeSignal::Table { times: vec![0.0, 0.125, 0.375], values: vec![500.0, 400.0, 300.0],
                dims: TEMPERATURE_DIMS, interp: Interp::Hold },
            TimeSignal::Table { times: vec![0.0, 0.25, 0.375], values: vec![300.0, 500.0, 400.0],
                dims: TEMPERATURE_DIMS, interp: Interp::Hold },
        ];
        let patches = [0, 1].map(|body| {
            let mut patch = declaration(&card, &curves[body], &context, schedules[body].clone());
            patch.region = if body == 0 { "left" } else { "right" }.into();
            patch
        }).to_vec();
        let mut run = original.into_sdirk2().with_material_ambient(cx, patches).unwrap();
        let control = AdaptiveHistoryConfig { absolute_temperature_k: 1e-10,
            relative_temperature: 0.0, min_step_s: 1e-5, max_step_s: 1.0, max_attempts: 12 };
        let phase = PhaseErrorTolerance { absolute_enthalpy_j_kg: 1000.0,
            relative_enthalpy: 0.0, absolute_liquid_fraction: 1.0 };
        let mut old = [500.0; 2];
        let mut previous = 0.0;
        let mut whole_heat = 0.0;
        for (end, baths) in [
            (0.125, [500.0, 300.0]), (0.25, [400.0, 300.0]),
            (0.375, [400.0, 500.0]), (0.5, [300.0, 400.0]),
        ] {
            let report = run.advance_adaptive_step_with_entropy(cx, 0.5, 1.0,
                config(), control, phase, 1e-8).unwrap();
            let physical = &report.physical;
            assert_eq!(physical.start_time_s, previous);
            assert_eq!(physical.end_time_s, end);
            assert_eq!(physical.rejected_attempts, 0);
            assert_eq!(physical.endpoint_solves, 3);
            let dt = (end - previous) / 2.0;
            let mut heat_by_bath = [0.0; 2];
            let initial = old;
            for half in [&physical.first_half, &physical.second_half] {
                for body in 0..2 {
                    let expected = stages(old[body], dt, [baths[body]; 2]);
                    let datum = if body == 0 { 0.0 } else { 4096.0 };
                    admit_stage(&half.first_stage, body, expected[0], datum, baths[body]);
                    admit_stage(&half.second_stage, body, expected[1], datum, baths[body]);
                    heat_by_bath[body] += dt * ((1.0 - G) * heat(expected[0], baths[body])
                        + G * heat(expected[1], baths[body]));
                    old[body] = expected[1];
                }
            }
            let supplied = -heat_by_bath.iter().sum::<f64>();
            let entropy_input: f64 = (0..2).map(|body| -heat_by_bath[body] / baths[body]).sum();
            let generation: f64 = (0..2).map(|body|
                heat_by_bath[body] * (1.0 / baths[body] - 1.0 / 400.0)).sum();
            close(physical.first_half.radiation_out_j + physical.second_half.radiation_out_j, -supplied);
            close(supplied, 3.0 * (old[0] - initial[0] + old[1] - initial[1]));
            close(report.entropy.thermal_storage_change_j_k, supplied / 400.0);
            close(report.entropy.radiation_entropy_input_j_k, entropy_input);
            close(report.entropy.radiation_generation_j_k, generation);
            close(report.entropy.temporal_defect_j_k, 0.0);
            assert!(!report.entropy.has_unclassified_inputs);
            report.entropy.admit(1e-8).unwrap();
            whole_heat += supplied;
            previous = end;
        }
        close(whole_heat, 3.0 * (old[0] - 500.0 + old[1] - 500.0));
        assert_eq!(run.conduction().time_s(), 0.5);
        assert_eq!(run.conduction().accepted_steps(), 8);
    });
}

#[test]
fn adaptive_bath_ramp_resamples_rejected_trials_and_accepts_only_fine_phase_heat() {
    let mesh = mesh(false);
    let card = source_card(1.0);
    let curve = phase_chart(&card, 0.0);
    let context = context();
    let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    with_cx(|cx| {
        let boundary = face_history(cx, &mesh);
        let signal = TimeSignal::Ramp { t_start: 0.0, t_end: 0.25,
            from: QtyAny::new(400.0, TEMPERATURE_DIMS), to: QtyAny::new(500.0, TEMPERATURE_DIMS) };
        let mut run = base(cx, &boundary, &conductivity, &curve, 500.0, 2).into_sdirk2()
            .with_material_ambient(cx, vec![declaration(&card, &curve, &context, signal)]).unwrap();
        let report = run.advance_adaptive_step_with_entropy(cx, 0.25, 0.25, config(),
            AdaptiveHistoryConfig { absolute_temperature_k: 1e-10, relative_temperature: 0.0,
                min_step_s: 0.001, max_step_s: 0.25, max_attempts: 8 },
            PhaseErrorTolerance { absolute_enthalpy_j_kg: 1.0,
                relative_enthalpy: 0.0, absolute_liquid_fraction: 1e-5 }, 1e-8).unwrap();
        assert_eq!(report.physical.rejected_attempts, 2);
        assert_eq!(report.physical.endpoint_solves, 9);
        assert_eq!(report.physical.end_time_s, 0.0625);
        assert_eq!(run.conduction().accepted_steps(), 2);
        assert_eq!(run.conduction().temperature(), &[400.0; 4]);
        let mut old = 500.0;
        let mut from = 0.0;
        let mut supplied = 0.0;
        let mut entropy_input = 0.0;
        for (physical, to) in [
            (&report.physical.first_half, report.physical.midpoint_time_s),
            (&report.physical.second_half, report.physical.end_time_s),
        ] {
            let dt = to - from;
            let baths = [400.0 + 400.0 * G.mul_add(dt, from), 400.0 + 400.0 * to];
            let expected = stages(old, dt, baths);
            old = admit_single(physical, old, dt, 0.0, baths);
            supplied -= physical.radiation_out_j;
            entropy_input -= dt * ((1.0 - G) * heat(expected[0], baths[0]) / baths[0]
                + G * heat(expected[1], baths[1]) / baths[1]);
            from = to;
        }
        close(3.0 * (old - 500.0), supplied);
        close(report.entropy.thermal_storage_change_j_k, supplied / 400.0);
        close(report.entropy.radiation_entropy_input_j_k, entropy_input);
        close(report.entropy.radiation_generation_j_k, supplied / 400.0 - entropy_input);
        close(report.entropy.temporal_defect_j_k, 0.0);
        report.entropy.admit(1e-8).unwrap();
    });
}

#[test]
fn refused_second_bath_stage_and_cancellation_preserve_the_accepted_prefix() {
    let mesh = mesh(false);
    let card = source_card(1.0);
    let curve = phase_chart(&card, 0.0);
    let cloned = curve.clone();
    let context = context();
    let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    let boundary = with_cx(|cx| face_history(cx, &mesh));
    with_cx(|cx| {
        assert!(base(cx, &boundary, &conductivity, &curve, 500.0, 4).into_sdirk2()
            .with_material_ambient(cx, vec![declaration(&card, &cloned, &context, constant(500.0))]).is_err());
        let attached = base(cx, &boundary, &conductivity, &curve, 500.0, 4).into_sdirk2()
            .with_material_ambient(cx, vec![declaration(&card, &curve, &context, constant(500.0))]).unwrap();
        assert!(attached.with_material_ambient(cx, Vec::new()).is_err());
    });
    let mut run = with_cx(|cx| {
        let patch = declaration(&card, &curve, &context, TimeSignal::Table {
            times: vec![0.0, 0.25], values: vec![500.0, -100.0],
            dims: TEMPERATURE_DIMS, interp: Interp::Hold });
        base(cx, &boundary, &conductivity, &curve, 500.0, 4).into_sdirk2()
            .with_material_ambient(cx, vec![patch]).unwrap()
    });
    with_cx(|cx| { run.advance_to(cx, 0.1, config()).unwrap(); });
    let saved = run.conduction().specific_enthalpy_j_kg().to_vec();
    with_cancelled_cx(|cx| assert!(run.advance_to(cx, 0.2, config()).is_err()));
    with_cx(|cx| {
        // The first stage's 500 K bath is valid, but the second sample is not.
        assert!(run.advance_to(cx, 0.3, config()).is_err());
        assert_eq!(run.conduction().specific_enthalpy_j_kg(), saved);
        assert_eq!(run.conduction().time_s(), 0.1);
        assert_eq!(run.conduction().accepted_steps(), 1);
        let old = saved[0];
        let step = run.advance_event_step(cx, 0.5, 1.0, config()).unwrap();
        assert_eq!(step.end_time_s, 0.25);
        admit_single(&step.report, old, 0.15, 0.0, [500.0; 2]);
        assert_eq!(run.conduction().accepted_steps(), 2);
        assert!(run.advance_event_step(cx, 0.5, 1.0, config()).is_err());
        assert_eq!(run.conduction().time_s(), 0.25);
        assert_eq!(run.conduction().accepted_steps(), 2);
    });
}
