//! G1/G3/G4 dynamic wetting: independent geometry, work, rate and source checks.
#[path = "contact_line_bridge/common.rs"]
mod common;
use common::*;
use fs_contact::capillary_bridge::contact_line::*;
use fs_exec::CancelGate;

#[test]
fn g1_equilibrium_recovers_original_bridge_and_actual_young_force() {
    let model = source(1.0, 10.0, 0.9, 1.1);
    let state = model.initialize(GAP, 1.0).unwrap();
    let p = model.evaluate(state).unwrap();
    let original = model
        .original()
        .evaluate(model.original().initialize(GAP, true).unwrap())
        .unwrap();
    close(
        p.surface_free_energy_j,
        original.surface_free_energy_j,
        2e-22,
    );
    close(p.normal_force_n, original.normal_force_n, 2e-18);
    close(p.wetted_span_m, original.wetted_span_m, 2e-18);
    close(p.footprint_energy_gradient_n, 0.0, 2e-20);
    let angle = 1.15;
    let p = model
        .evaluate(model.initialize(GAP, angle).unwrap())
        .unwrap();
    let radius = GAP / (2.0 * angle.cos());
    let n = 8192;
    let dy = GAP / f64::from(n);
    let mut area = 0.0;
    for i in 0..n {
        let y = -0.5 * GAP + (f64::from(i) + 0.5) * dy;
        let half_width =
            0.5 * p.wetted_span_m + radius * angle.sin() - (radius * radius - y * y).sqrt();
        area += 2.0 * half_width * dy;
    }
    close(WIDTH * area, VOLUME, 3e-19);
    close(angle_at(GAP, p.wetted_span_m), angle, 2e-14);
    let force = p.liquid_minus_gas_pressure_pa * p.wetted_span_m * WIDTH
        - 2.0 * SIGMA * WIDTH * angle.sin();
    close(force, p.normal_force_n, 2e-18);
    close(
        p.footprint_energy_gradient_n,
        2.0 * SIGMA * WIDTH * (angle.cos() - 1.0_f64.cos()),
        2e-20,
    );
}

#[test]
fn g1_fixed_gap_relaxation_converts_original_surface_energy_into_line_loss() {
    let model = source(1.0, 10.0, 1.0, 1.0);
    let mut state = model.initialize(GAP, 1.16).unwrap();
    let initial = model.evaluate(state).unwrap();
    let mut loss = 0.0;
    for _ in 0..100 {
        let step = model
            .advance(state, GAP, 0.001, &CancelGate::new())
            .unwrap();
        assert_eq!(step.regime, ContactLineRegime::Advancing);
        assert!(step.profile.contact_angle_radians < state.contact_angle_radians());
        assert!(step.profile.contact_angle_radians >= 1.0);
        assert!(step.viscous_line_dissipation_j > 0.0);
        close(step.pinning_dissipation_j, 0.0, 0.0);
        close(step.work_on_liquid_j, 0.0, 0.0);
        close(
            step.surface_energy_change_j + step.contact_line_dissipation_j,
            0.0,
            3e-19,
        );
        assert!(
            step.estimated_force_error_n
                <= config().force_absolute_tolerance_n
                    + config().force_relative_tolerance * step.mean_normal_force_n.abs()
        );
        loss += step.contact_line_dissipation_j;
        state = step.state();
    }
    let final_profile = model.evaluate(state).unwrap();
    close(
        final_profile.surface_free_energy_j + loss,
        initial.surface_free_energy_j,
        2e-18,
    );
    assert!(state.contact_angle_radians() < 1.03);
    assert_eq!(final_profile.liquid_volume_m3, VOLUME);
}

#[test]
fn g1_pinned_footprint_depins_with_separate_advancing_and_receding_losses() {
    let model = source(1.0, 1.0, 0.95, 1.05);
    let state = model.initialize(GAP, 1.0).unwrap();
    let pinned = model
        .advance(state, GAP + 2e-6, 0.01, &CancelGate::new())
        .unwrap();
    assert_eq!(pinned.regime, ContactLineRegime::Pinned);
    assert_eq!(pinned.state().wetted_span_m(), state.wetted_span_m());
    assert_eq!(pinned.contact_line_dissipation_j, 0.0);
    close(
        pinned.surface_energy_change_j,
        pinned.work_on_liquid_j,
        3e-19,
    );
    for (gap, regime) in [
        (GAP - 2e-5, ContactLineRegime::Advancing),
        (GAP + 2e-5, ContactLineRegime::Receding),
    ] {
        let step = model.advance(state, gap, 0.01, &CancelGate::new()).unwrap();
        assert_eq!(step.regime, regime);
        assert!(step.pinning_dissipation_j > 0.0);
        assert!(step.viscous_line_dissipation_j > 0.0);
        close(
            step.surface_energy_change_j + step.contact_line_dissipation_j,
            step.work_on_liquid_j,
            4e-19,
        );
        assert!(if regime == ContactLineRegime::Advancing {
            step.state().wetted_span_m() > state.wetted_span_m()
        } else {
            step.state().wetted_span_m() < state.wetted_span_m()
        });
    }
}

#[test]
fn g3_original_line_friction_changes_relaxation_and_implicit_force_tangent() {
    for friction in [1.0, 10.0] {
        let model = source(1.0, friction, 0.98, 1.02);
        for angle in [0.86, 1.0, 1.14] {
            let state = model.initialize(GAP, angle).unwrap();
            let h = GAP + 1e-6;
            let step = model.advance(state, h, 0.003, &CancelGate::new()).unwrap();
            let epsilon = 1e-9;
            let plus = model
                .advance(state, h + epsilon, 0.003, &CancelGate::new())
                .unwrap();
            let minus = model
                .advance(state, h - epsilon, 0.003, &CancelGate::new())
                .unwrap();
            let derivative =
                (plus.mean_normal_force_n - minus.mean_normal_force_n) / (2.0 * epsilon);
            close(
                step.mean_force_gap_derivative_n_per_m,
                derivative,
                3e-6 * derivative.abs().max(1.0),
            );
        }
    }
    let fast = source(1.0, 1.0, 1.0, 1.0);
    let slow = source(1.0, 10.0, 1.0, 1.0);
    let a = fast
        .advance(
            fast.initialize(GAP, 1.16).unwrap(),
            GAP,
            0.001,
            &CancelGate::new(),
        )
        .unwrap();
    let b = slow
        .advance(
            slow.initialize(GAP, 1.16).unwrap(),
            GAP,
            0.001,
            &CancelGate::new(),
        )
        .unwrap();
    assert!(a.profile.contact_angle_radians < b.profile.contact_angle_radians);
    assert!(a.contact_line_dissipation_j > b.contact_line_dissipation_j);
}

#[test]
fn g3_wetting_trajectory_converges_to_independent_ode_under_timestep_refinement() {
    let model = source(1.0, 10.0, 1.0, 1.0);
    let initial = model.initialize(GAP, 1.14).unwrap();
    let end = 0.02;
    let rhs = |span: f64| 2.0 * SIGMA / 10.0 * (1.0_f64.cos() - angle_at(GAP, span).cos());
    let mut reference = initial.wetted_span_m();
    let oracle_dt = end / 4096.0;
    for _ in 0..4096 {
        let k1 = rhs(reference);
        let k2 = rhs(reference + 0.5 * oracle_dt * k1);
        let k3 = rhs(reference + 0.5 * oracle_dt * k2);
        let k4 = rhs(reference + oracle_dt * k3);
        reference += oracle_dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
    let errors = [10, 20, 40].map(|steps| {
        let mut state = initial;
        for _ in 0..steps {
            state = model
                .advance(state, GAP, end / f64::from(steps), &CancelGate::new())
                .unwrap()
                .state();
        }
        (state.wetted_span_m() - reference).abs()
    });
    assert!(
        errors[0] / errors[1] > 3.5 && errors[1] / errors[2] > 3.5,
        "{errors:?}"
    );
}

#[test]
fn g0_g4_source_geometry_accuracy_and_cancel_refusals_preserve_exact_state() {
    let model = source(1.0, 10.0, 0.95, 1.05);
    let state = model.initialize(GAP, 1.12).unwrap();
    assert!(matches!(
        source(1.0, 11.0, 0.95, 1.05).evaluate(state),
        Err(ContactLineError::StateMismatch)
    ));
    assert!(
        SourceContactLineBridge::new(original(1.0, [10.0, 11.0], 0.95, 1.05), config()).is_err()
    );
    assert!(model.initialize(GAP, 1.5).is_err());
    assert!(
        model
            .advance(state, GAP + 1e-3, 0.01, &CancelGate::new())
            .is_err()
    );
    assert!(model.advance(state, GAP, 0.0, &CancelGate::new()).is_err());
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        model.advance(state, GAP, 0.001, &gate),
        Err(ContactLineError::Cancelled)
    );
    let mut limited = config();
    limited.maximum_evaluations = 2;
    let bounded =
        SourceContactLineBridge::new(original(1.0, [10.0; 2], 0.95, 1.05), limited).unwrap();
    let saved = bounded.initialize(GAP, 1.12).unwrap();
    assert_eq!(
        bounded.advance(saved, GAP, 0.001, &CancelGate::new()),
        Err(ContactLineError::Budget)
    );
    assert_eq!(saved.contact_angle_radians(), 1.12);
    let mut unresolved = config();
    unresolved.force_absolute_tolerance_n = 1e-30;
    unresolved.force_relative_tolerance = 1e-16;
    let strict =
        SourceContactLineBridge::new(original(1.0, [10.0; 2], 0.95, 1.05), unresolved).unwrap();
    let strict_state = strict.initialize(GAP, 1.12).unwrap();
    assert!(
        strict
            .advance(strict_state, GAP, 0.001, &CancelGate::new())
            .is_err()
    );
    let mut speed = config();
    speed.maximum_contact_line_speed_m_s = 1e-12;
    let slow = SourceContactLineBridge::new(original(1.0, [10.0; 2], 0.95, 1.05), speed).unwrap();
    assert!(matches!(
        slow.advance(
            slow.initialize(GAP, 1.12).unwrap(),
            GAP,
            0.001,
            &CancelGate::new()
        ),
        Err(ContactLineError::Domain("contact-line speed"))
    ));
    assert_eq!(
        model.advance(state, GAP, 0.001, &CancelGate::new()),
        model.advance(state, GAP, 0.001, &CancelGate::new())
    );
}
