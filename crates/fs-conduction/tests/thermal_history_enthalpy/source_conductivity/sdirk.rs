//! Second-order scheduled heat with exact original nonlinear conductivity moments.
use super::*;
use fs_conduction::transient::enthalpy::sdirk::SDIRK2_GAMMA as G;

// Algebraic stage carry, not an admitted or queried material state.
fn stages(old: f64, dt: f64, mass: f64, q1: f64, q2: f64) -> [f64; 2] {
    let first = endpoint(old, G * dt, mass, q1);
    let carry = dt * (1.0 - G) * ((10.0 / 6.0) * mean_k(first) + q1 / 24.0) / mass;
    [first, endpoint(old + carry, G * dt, mass, q2)]
}

#[test]
fn source_owned_sdirk_history_samples_both_stage_times_without_resetting_mass_or_phase() {
    let mesh = mesh(false); let card = card();
    for datum in [0.0, 4096.0] {
        let curve = chart(&card, datum);
        let source = source(&card, &curve, [datum - 100.0, datum + 500.0]);
        with_cx(|cx| {
            let materials = [MaterialPhaseTransport::try_new_isotropic(cx, &source, 24.0).unwrap()];
            let boundary = heated(cx, &mesh);
            let initial = [datum + 110.0, datum + 110.0, datum + 110.0, datum + 20.0];
            // Positive P1 rho0 gives exact free-node mass1.4kg, unlike either
            // the source's equilibrium density or uniform reference24kg/m3.
            let mut original = EnthalpyHistoryRun::from_phase_materials_with_nodal_reference_density(cx,
                histories(&boundary, HeatSourceHistory::Uniform(TimeSignal::Ramp {
                    t_start: 0.0, t_end: 2.0, from: QtyAny::new(240.0, VOLUMETRIC_SOURCE_DIMS),
                    to: QtyAny::new(480.0, VOLUMETRIC_SOURCE_DIMS) })),
                &materials, &[0], &[12.0, 24.0, 36.0, 48.0], &initial, budget(&mesh, 4)).unwrap();
            let first = original.advance_to(cx, 0.1, config()).unwrap();
            let mut old = first.specific_enthalpy_j_kg[3] - datum;
            let saved = original.specific_enthalpy_j_kg().to_vec();
            let mut run = original.into_sdirk2();
            assert_eq!(run.conduction().time_s(), 0.1);
            assert_eq!(run.conduction().accepted_steps(), 1);
            assert_eq!(run.conduction().specific_enthalpy_j_kg(), saved);
            assert!(std::ptr::eq(run.conduction().storage().phase_conductivity_source(0).unwrap(), &source));
            assert!(std::ptr::eq(run.conduction().storage().phase_curve_at(3).unwrap(), &curve));
            close(run.conduction().storage().reference_nodal_masses_kg()[3], 1.4);
            let mut from = 0.1;
            for to in [0.3, 0.5, 0.7] {
                let dt = to - from;
                let q = |t: f64| 240.0 + 120.0 * t;
                let expected = stages(old, dt, 1.4, q(from + G * dt), q(to));
                let step = run.advance_to(cx, to, config()).unwrap();
                for (stage, h) in [(&step.first_stage, expected[0]), (&step.second_stage, expected[1])] {
                    close(stage.specific_enthalpy_j_kg[3] - datum, h);
                    close(stage.liquid_mass_fraction[3], h / 100.0);
                    assert_eq!(stage.temperature, [410.0, 410.0, 410.0, 400.0]);
                }
                let integrated_source = dt * (q(from) + q(to)) / 12.0; // Exact linear-in-time integral, V=1/6.
                close(step.source_j, integrated_source);
                close(step.stored_energy_change_j, 1.4 * (expected[1] - old));
                close(step.dirichlet_in_j + step.source_j, step.stored_energy_change_j);
                assert!(step.energy_residual_j.abs() <= config().energy_tolerance_j);
                let heat = run.conduction().storage().phase_transport_heat_flux(cx,
                    run.conduction().specific_enthalpy_j_kg()).unwrap();
                close(heat[0][2], 10.0 * mean_k(expected[1]));
                old = expected[1]; from = to;
            }
            assert_eq!(run.conduction().accepted_steps(), 4);
        });
    }
}

#[test]
fn source_sdirk_events_and_phase_adaptation_keep_original_conductivity_and_heat_accounts() {
    let mesh = mesh(false); let card = card(); let curve = chart(&card, 0.0);
    let source = source(&card, &curve, [-100.0, 500.0]);
    with_cx(|cx| {
        let materials = [MaterialPhaseTransport::try_new_isotropic(cx, &source, 24.0).unwrap()];
        let boundary = heated(cx, &mesh);
        let create = || EnthalpyHistoryRun::from_phase_materials(cx, histories(&boundary,
            HeatSourceHistory::Uniform(TimeSignal::Table { times: vec![0.0, 0.5],
                values: vec![480.0, -480.0], dims: VOLUMETRIC_SOURCE_DIMS, interp: Interp::Hold })),
            &materials, &[0], &[110.0, 110.0, 110.0, 20.0], budget(&mesh, 4)).unwrap().into_sdirk2();
        let mut events = create();
        let first = events.advance_event_step(cx, 1.0, 1.0, config()).unwrap();
        assert_eq!(first.end_time_s, 0.5);
        let h1 = stages(20.0, 0.5, 1.0, 480.0, 480.0);
        close(first.report.second_stage.specific_enthalpy_j_kg[3], h1[1]);
        close(first.report.source_j, 40.0);
        let second = events.advance_event_step(cx, 1.0, 1.0, config()).unwrap();
        let h2 = stages(h1[1], 0.5, 1.0, -480.0, -480.0);
        assert_eq!(second.end_time_s, 1.0);
        close(second.report.second_stage.specific_enthalpy_j_kg[3], h2[1]);
        close(second.report.source_j, -40.0);
        close(first.report.stored_energy_change_j + second.report.stored_energy_change_j,
            first.report.dirichlet_in_j + second.report.dirichlet_in_j);

        let mut adaptive = create();
        let accepted = adaptive.advance_adaptive_step(cx, 0.5, 0.5, config(),
            AdaptiveHistoryConfig { absolute_temperature_k: 1e-10, relative_temperature: 0.0,
                min_step_s: 1e-6, max_step_s: 0.5, max_attempts: 20 },
            PhaseErrorTolerance { absolute_enthalpy_j_kg: 100.0, relative_enthalpy: 0.0,
                absolute_liquid_fraction: 1e-8 }).unwrap();
        assert!(accepted.rejected_attempts > 0 && accepted.end_time_s < 0.5);
        assert!(accepted.normalized_error <= 1.0);
        let hm = stages(20.0, accepted.midpoint_time_s, 1.0, 480.0, 480.0)[1];
        let hf = stages(hm, accepted.end_time_s - accepted.midpoint_time_s, 1.0, 480.0, 480.0)[1];
        close(accepted.first_half.second_stage.specific_enthalpy_j_kg[3], hm);
        close(accepted.second_half.second_stage.specific_enthalpy_j_kg[3], hf);
        assert_eq!(accepted.second_half.second_stage.temperature[3], 400.0);
        assert_eq!(adaptive.conduction().accepted_steps(), 2);
        close(adaptive.conduction().specific_enthalpy_j_kg()[3], hf);
    });
}

#[test]
fn failed_second_source_stage_and_unrepresentable_clock_leave_sdirk_history_resumable() {
    let mesh = mesh(false); let card = card(); let curve = chart(&card, 0.0);
    let source = source(&card, &curve, [40.0, 60.0]);
    with_cx(|cx| {
        let materials = [MaterialPhaseTransport::try_new_isotropic(cx, &source, 24.0).unwrap()];
        let boundary = natural(cx, &mesh);
        let create = || EnthalpyHistoryRun::from_phase_materials(cx, histories(&boundary,
            HeatSourceHistory::Fixed(ScalarField::Uniform(480.0))),
            &materials, &[0], &[50.0; 4], budget(&mesh, 2)).unwrap().into_sdirk2();
        let mut run = create();
        // First stage h=55.857... is supported; second stage h=70 is not.
        assert!(run.advance_to(cx, 1.0, config()).is_err());
        with_cancelled_cx(|cancelled| assert!(run.advance_to(cancelled, 0.25, config()).is_err()));
        assert_eq!(run.conduction().specific_enthalpy_j_kg(), &[50.0; 4]);
        assert_eq!(run.conduction().accepted_steps(), 0);
        assert_eq!(run.conduction().time_s(), 0.0);
        let step = run.advance_to(cx, 0.25, config()).unwrap();
        for &h in &step.second_stage.specific_enthalpy_j_kg { close(h, 55.0); }
        close(step.source_j, 20.0); close(step.stored_energy_change_j, 20.0);
        assert!(run.advance_to(cx, 0.25_f64.next_up(), config()).is_err());
        assert_eq!(run.conduction().time_s(), 0.25);
        assert_eq!(run.conduction().accepted_steps(), 1);
        assert_eq!(run.conduction().specific_enthalpy_j_kg(), step.second_stage.specific_enthalpy_j_kg);
        let replay = create().advance_to(cx, 0.25, config()).unwrap();
        assert_eq!(replay.second_stage.specific_enthalpy_j_kg, step.second_stage.specific_enthalpy_j_kg);
    });
}
