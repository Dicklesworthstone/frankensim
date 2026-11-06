use super::*;
fn liquid(nu: f64, sigma: f64) -> SlenderJetParameters {
    SlenderJetParameters {
        density_kg_m3: 1.0,
        surface_tension_n_m: sigma,
        kinematic_viscosity_m2_s: nu,
    }
}
fn limits() -> SlenderJetLimits {
    SlenderJetLimits {
        minimum_radius_m: 0.1,
        maximum_abs_slope: 0.75,
        minimum_radius_cells: 2.0,
    }
}
fn budget() -> SlenderJetStepBudget {
    SlenderJetStepBudget {
        maximum_step_s: 0.01,
        cfl: 0.2,
        maximum_substeps: 100_000,
    }
}
fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-30),
        "{actual:e} != {expected:e}"
    );
}
fn eigen_profile(n: usize, amplitude: f64, nu: f64) -> (f64, Vec<f64>, Vec<f64>, f64) {
    let k = 0.5;
    let length = core::f64::consts::TAU / k;
    let dx = length / n as f64;
    let half_damping = 1.5 * nu * k * k;
    let drive = 0.5 * k * k * (1.0 - k * k);
    let growth = drive / ((half_damping * half_damping + drive).sqrt() + half_damping);
    let radius = (0..n)
        .map(|i| 1.0 + amplitude * fs_math::det::cos(k * i as f64 * dx))
        .collect();
    let velocity = (0..n)
        .map(|i| -2.0 * growth * amplitude / k * fs_math::det::sin(k * (i as f64 + 0.5) * dx))
        .collect();
    (length, radius, velocity, growth)
}

fn choose_implicit(jet: &mut PeriodicSlenderJet, gate: &CancelGate) {
    jet.set_integration_method(SlenderJetIntegrationMethod::ImplicitViscosity, gate)
        .unwrap();
}

#[test]
fn g1_cyclic_viscosity_solve_matches_a_nonuniform_manufactured_state() {
    let gate = CancelGate::new_clock_free();
    let n = 17;
    let mut state = Fields::zeros(n, &gate).unwrap();
    let u: Vec<f64> = (0..n)
        .map(|i| 0.7 + 0.2 * fs_math::det::sin(core::f64::consts::TAU * i as f64 / n as f64))
        .collect();
    for i in 0..n {
        state.a[i] = 1.0 + 0.25 * fs_math::det::cos(core::f64::consts::TAU * i as f64 / n as f64);
    }
    let (nu, dt, dx) = (4.0, 0.3, 0.2);
    let coefficient = 3.0 * nu * dt / dx / dx;
    for i in 0..n {
        let l = previous(i, n);
        let r = next(i, n);
        let mass = average(state.a[i], state.a[r]);
        state.q[i] =
            mass * u[i] - coefficient * (state.a[r] * (u[r] - u[i]) - state.a[i] * (u[i] - u[l]));
    }
    let initial_momentum: f64 = state.q.iter().sum();
    let mut workspace = Workspace::new(n, &gate).unwrap();
    workspace
        .solve
        .apply(&mut state, nu, dx, dt, &gate)
        .unwrap();
    for i in 0..n {
        let actual = state.q[i] / average(state.a[i], state.a[next(i, n)]);
        near(actual, u[i], 5e-12);
    }
    near(state.q.iter().sum(), initial_momentum, 3e-13);
}

#[test]
fn g1_high_viscosity_finishes_with_less_work_and_agrees_with_explicit_reference() {
    let gate = CancelGate::new_clock_free();
    let (length, radius, mut velocity, _) = eigen_profile(64, 0.04, 50.0);
    for value in &mut velocity {
        *value += 0.05;
    }
    let mut explicit = PeriodicSlenderJet::new(
        liquid(50.0, 1.0),
        length,
        &radius,
        &velocity,
        limits(),
        &gate,
    )
    .unwrap();
    let mut implicit = explicit.clone();
    choose_implicit(&mut implicit, &gate);
    let short_budget = SlenderJetStepBudget {
        maximum_substeps: 1000,
        ..budget()
    };
    assert_eq!(
        explicit.advance(0.04, short_budget, &gate),
        Err(NonlinearJetError::SubstepBudget)
    );
    let result = implicit.advance(0.04, short_budget, &gate).unwrap();
    let reference = explicit.advance(0.04, budget(), &gate).unwrap();
    assert_eq!(result.outcome, SlenderJetOutcome::Completed);
    assert_eq!(
        result.integration_method,
        SlenderJetIntegrationMethod::ImplicitViscosity
    );
    assert!(result.accepted_substeps * 20 < reference.accepted_substeps);
    assert!(reference.temporal_error.is_none());
    let error = result.temporal_error.unwrap();
    assert!(error.maximum_radius_error_m.is_finite());
    assert!(error.maximum_velocity_error_m_s.is_finite());
    for (actual, expected) in implicit.radii_m().zip(explicit.radii_m()) {
        near(actual, expected, 2e-7);
    }
    for (actual, expected) in implicit.velocities_m_s().zip(explicit.velocities_m_s()) {
        assert!((actual - expected).abs() < 2e-7);
    }
    near(result.after.volume_m3, result.before.volume_m3, 3e-12);
    near(
        result.after.axial_momentum_kg_m_s,
        result.before.axial_momentum_kg_m_s,
        3e-12,
    );
    assert!(result.viscous_dissipation_j >= 0.0);
    assert!(result.upwind_dissipation_j >= 0.0);
    assert!(result.energy_balance_residual_j.is_finite());
}

#[test]
fn g1_zero_viscosity_reduces_to_explicit_ssprk3() {
    let gate = CancelGate::new_clock_free();
    let (length, radius, velocity, _) = eigen_profile(64, 0.04, 0.0);
    let mut explicit = PeriodicSlenderJet::new(
        liquid(0.0, 1.0),
        length,
        &radius,
        &velocity,
        limits(),
        &gate,
    )
    .unwrap();
    let mut implicit = explicit.clone();
    choose_implicit(&mut implicit, &gate);
    let a = explicit.advance(0.05, budget(), &gate).unwrap();
    let b = implicit.advance(0.05, budget(), &gate).unwrap();
    assert_eq!(a.outcome, b.outcome);
    for (x, y) in implicit.radii_m().zip(explicit.radii_m()) {
        near(x, y, 5e-13);
    }
    for (x, y) in implicit.velocities_m_s().zip(explicit.velocities_m_s()) {
        assert!((x - y).abs() < 1e-12);
    }
    assert_eq!(b.viscous_dissipation_j, 0.0);
}

#[test]
fn g1_time_refinement_is_third_order_for_resolved_viscous_decay() {
    let gate = CancelGate::new_clock_free();
    let n = 64;
    let length = 4.0 * core::f64::consts::PI;
    let dx = length / n as f64;
    let nu = 2.0;
    let amplitude = 1e-8;
    let duration = 0.2;
    let wave = |i: usize| fs_math::det::sin(core::f64::consts::TAU * (i as f64 + 0.5) / n as f64);
    let velocity: Vec<f64> = (0..n).map(|i| amplitude * wave(i)).collect();
    let discrete_k = 2.0 * fs_math::det::sin(0.25 * dx) / dx;
    let expected = amplitude * fs_math::det::exp(-3.0 * nu * discrete_k * discrete_k * duration);
    let mut errors = Vec::new();
    let mut indicators = Vec::new();
    for maximum_step_s in [0.04, 0.02, 0.01] {
        let mut jet = PeriodicSlenderJet::new(
            liquid(nu, 0.0),
            length,
            &vec![1.0; n],
            &velocity,
            limits(),
            &gate,
        )
        .unwrap();
        choose_implicit(&mut jet, &gate);
        let result = jet
            .advance(
                duration,
                SlenderJetStepBudget {
                    maximum_step_s,
                    ..budget()
                },
                &gate,
            )
            .unwrap();
        assert_eq!(result.outcome, SlenderJetOutcome::Completed);
        let observed = 2.0 / n as f64
            * jet
                .velocities_m_s()
                .enumerate()
                .map(|(i, u)| u * wave(i))
                .sum::<f64>();
        errors.push((observed - expected).abs());
        indicators.push(result.temporal_error.unwrap().maximum_velocity_error_m_s);
        assert!(result.viscous_dissipation_j > 0.0);
    }
    assert!(
        errors[1] < 0.2 * errors[0] && errors[2] < 0.25 * errors[1],
        "{errors:?}"
    );
    assert!(
        indicators[1] < 0.2 * indicators[0] && indicators[2] < 0.2 * indicators[1],
        "{indicators:?}"
    );
}

#[test]
fn g4_policy_selection_cancellation_and_work_refusal_are_transactional() {
    let gate = CancelGate::new_clock_free();
    let (length, radius, velocity, _) = eigen_profile(64, 0.04, 1.0);
    let mut jet = PeriodicSlenderJet::new(
        liquid(1.0, 1.0),
        length,
        &radius,
        &velocity,
        limits(),
        &gate,
    )
    .unwrap();
    let a = jet.radius_squared_m2().to_vec();
    let q = jet.radius_squared_flux_m3_s().to_vec();
    let before = jet.diagnostics();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        jet.set_integration_method(SlenderJetIntegrationMethod::ImplicitViscosity, &cancelled),
        Err(NonlinearJetError::Cancelled)
    );
    assert_eq!(
        jet.integration_method(),
        SlenderJetIntegrationMethod::ExplicitSspRk3
    );
    choose_implicit(&mut jet, &gate);
    assert_eq!(jet.radius_squared_m2(), a);
    assert_eq!(jet.radius_squared_flux_m3_s(), q);
    assert_eq!(jet.diagnostics(), before);
    let mut fresh = jet.clone();
    assert_eq!(
        jet.advance(
            0.1,
            SlenderJetStepBudget {
                maximum_substeps: 1,
                ..budget()
            },
            &gate
        ),
        Err(NonlinearJetError::SubstepBudget)
    );
    assert_eq!(
        jet.advance(0.1, budget(), &cancelled),
        Err(NonlinearJetError::Cancelled)
    );
    assert_eq!(jet.radius_squared_m2(), a);
    assert_eq!(jet.radius_squared_flux_m3_s(), q);
    assert_eq!(jet.diagnostics(), before);
    assert_eq!(jet.time_s(), 0.0);
    assert_eq!(
        jet.advance(0.02, budget(), &gate),
        fresh.advance(0.02, budget(), &gate)
    );
    assert_eq!(jet.radius_squared_m2(), fresh.radius_squared_m2());
    assert_eq!(
        jet.radius_squared_flux_m3_s(),
        fresh.radius_squared_flux_m3_s()
    );
    jet.set_integration_method(SlenderJetIntegrationMethod::ExplicitSspRk3, &gate)
        .unwrap();
    assert_eq!(
        jet.integration_method(),
        SlenderJetIntegrationMethod::ExplicitSspRk3
    );
    assert_eq!(jet.diagnostics(), fresh.diagnostics());
}

#[test]
fn g0_implicit_domain_stop_keeps_only_admitted_substeps() {
    let gate = CancelGate::new_clock_free();
    let (length, radius, velocity, _) = eigen_profile(64, 0.05, 0.05);
    let mut jet = PeriodicSlenderJet::new(
        liquid(0.05, 1.0),
        length,
        &radius,
        &velocity,
        SlenderJetLimits {
            minimum_radius_m: 0.94,
            ..limits()
        },
        &gate,
    )
    .unwrap();
    choose_implicit(&mut jet, &gate);
    let result = jet.advance(1.0, budget(), &gate).unwrap();
    let SlenderJetOutcome::DomainLimit(violation) = result.outcome else {
        panic!("expected a neck stop: {:?}", result.outcome);
    };
    assert_eq!(violation.kind, SlenderJetLimitKind::MinimumRadius);
    assert!(result.accepted_substeps > 0);
    assert!(result.end_time_s > 0.0 && result.end_time_s < 1.0);
    assert!(result.after.minimum_radius_m >= 0.94);
    near(result.after.volume_m3, result.before.volume_m3, 3e-12);
    assert_eq!(result.after, jet.diagnostics());
    assert_eq!(result.end_time_s, jet.time_s());
}
