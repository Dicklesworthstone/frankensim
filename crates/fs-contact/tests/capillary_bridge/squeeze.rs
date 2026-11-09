//! G0/G1/G3/G4: independent Reynolds pressure integration, path work and domains.
use super::*;

const SMALL_VOLUME: f64 = 1e-12;
const GAP: f64 = 2e-6;

fn flow_limits() -> PlanarSqueezeFlowLimits {
    PlanarSqueezeFlowLimits {
        maximum_gap_to_half_span_ratio: 0.1,
        maximum_reynolds_number: 0.1,
        maximum_unsteady_stokes_number: 0.1,
        ambient_pressure_pa: 100_000.0,
        minimum_liquid_pressure_pa: 2_000.0,
        maximum_liquid_pressure_pa: 200_000.0,
    }
}
fn small_limits() -> BridgeLimits {
    BridgeLimits {
        minimum_gap_m: 1e-6,
        maximum_capillary_number: 0.01,
        maximum_weber_number: 0.01,
        ..limits()
    }
}
fn small_model() -> SourcePlanarCapillaryBridge {
    let (liquid, wetting) = sources(PI / 3.0, 0.001, 300.0, true);
    SourcePlanarCapillaryBridge::new(liquid, wetting, SMALL_VOLUME, WIDTH, small_limits()).unwrap()
}
fn viscous_model() -> SourcePlanarCapillaryBridge {
    small_model()
        .with_no_slip_squeeze_flow(flow_limits())
        .unwrap()
}

// Independent spatial integration of the parabolic Reynolds pressure, followed
// by time quadrature; it never uses the production resistance polynomial.
fn pressure_force(m: &SourcePlanarCapillaryBridge, h: f64, velocity: f64) -> f64 {
    let p = m.evaluate(m.initialize(h, true).unwrap()).unwrap();
    let span = p.wetted_span_m;
    let n = 16;
    let mut integral = 0.0;
    for i in 0..=n {
        let x = span * (i as f64 / n as f64 - 0.5);
        let pressure =
            6.0 * m.liquid().dynamic_viscosity_pa_s() * velocity * (x * x - 0.25 * span * span)
                / h.powi(3);
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight * pressure;
    }
    WIDTH * span * integral / (3.0 * n as f64)
}

#[test]
fn g1_bulk_force_matches_pressure_and_time_integrals_with_exact_path_work() {
    let base = small_model();
    let m = viscous_model();
    let gate = CancelGate::new_clock_free();
    for (h0, h1) in [(GAP, 1.1 * GAP), (1.1 * GAP, GAP)] {
        let dt = 0.05;
        let old = m.initialize(h0, true).unwrap();
        let result = m.advance(old, h1, dt, &gate).unwrap();
        let cap = base
            .advance(base.initialize(h0, true).unwrap(), h1, dt, &gate)
            .unwrap();
        let n = 256;
        let mut force_integral = 0.0;
        for i in 0..=n {
            let h = h0 + (h1 - h0) * i as f64 / n as f64;
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            force_integral += weight * pressure_force(&base, h, (h1 - h0) / dt);
        }
        let bulk_force = force_integral / (3.0 * n as f64);
        near(
            result.mean_normal_force_n - cap.mean_normal_force_n,
            bulk_force,
            2e-11,
            1e-9,
        );
        near(
            result.bulk_viscous_dissipation_j,
            -bulk_force * (h1 - h0),
            2e-11,
            1e-15,
        );
        assert!(result.bulk_viscous_dissipation_j > 0.0);
        assert_eq!(result.surface_energy_change_j, cap.surface_energy_change_j);
        assert_eq!(result.irreversible_release_j, 0.0);
        near(
            result.work_on_liquid_j,
            result.surface_energy_change_j + result.bulk_viscous_dissipation_j,
            2e-12,
            1e-12,
        );
        near(result.energy_balance_residual_j, 0.0, 2e-12, 1e-12);
        assert_eq!(cap.bulk_viscous_dissipation_j, 0.0);
        assert_eq!(m.advance(old, h1, dt, &gate).unwrap(), result);
    }
}

#[test]
fn g1_mean_tangent_includes_bulk_resistance_at_rest_and_during_both_directions() {
    let m = viscous_model();
    let gate = CancelGate::new_clock_free();
    let old = m.initialize(GAP, true).unwrap();
    let dt = 0.05;
    for endpoint in [0.9 * GAP, GAP, 1.1 * GAP] {
        let h = 1e-5 * GAP;
        let step = m.advance(old, endpoint, dt, &gate).unwrap();
        let minus = m.advance(old, endpoint - h, dt, &gate).unwrap();
        let plus = m.advance(old, endpoint + h, dt, &gate).unwrap();
        near(
            (plus.mean_normal_force_n - minus.mean_normal_force_n) / (2.0 * h),
            step.mean_force_gap_derivative_n_per_m,
            2e-8,
            1.0,
        );
    }
    let rest = m.advance(old, GAP, dt, &gate).unwrap();
    assert_eq!(rest.bulk_viscous_dissipation_j, 0.0);
    assert_eq!(rest.mean_normal_force_n, rest.profile.normal_force_n);
    let resistance = -pressure_force(&m, GAP, 1.0);
    near(
        rest.mean_force_gap_derivative_n_per_m,
        0.5 * rest.profile.force_gap_derivative_n_per_m - resistance / dt,
        2e-14,
        1.0,
    );
}

#[test]
fn g3_constant_speed_subdivision_and_reversal_preserve_dissipation() {
    let m = viscous_model();
    let gate = CancelGate::new_clock_free();
    let initial = m.initialize(GAP, true).unwrap();
    let end = 1.15 * GAP;
    let dt = 0.2;
    let full = m.advance(initial, end, dt, &gate).unwrap();
    let reverse = m.advance(full.state(), GAP, dt, &gate).unwrap();
    near(
        full.bulk_viscous_dissipation_j,
        reverse.bulk_viscous_dissipation_j,
        2e-14,
        1e-15,
    );
    near(
        full.work_on_liquid_j + reverse.work_on_liquid_j,
        full.bulk_viscous_dissipation_j + reverse.bulk_viscous_dissipation_j,
        1e-12,
        1e-14,
    );
    let mut state = initial;
    let mut loss = 0.0;
    let mut work = 0.0;
    for i in 1..=17 {
        let next = GAP + (end - GAP) * i as f64 / 17.0;
        let step = m.advance(state, next, dt / 17.0, &gate).unwrap();
        loss += step.bulk_viscous_dissipation_j;
        work += step.work_on_liquid_j;
        state = step.state();
    }
    near(loss, full.bulk_viscous_dissipation_j, 2e-13, 1e-15);
    near(work, full.work_on_liquid_j, 2e-13, 1e-15);
}

#[test]
fn g3_original_source_refresh_preserves_squeeze_selection_and_uses_new_viscosity() {
    let gate = CancelGate::new_clock_free();
    let (card, interfaces) = temperature_cards(true, 1.0, 1.0);
    let (liquid, wetting) = resolved_sources(&card, &interfaces, 300.0, true);
    let source =
        SourcePlanarCapillaryBridge::new(liquid, wetting, SMALL_VOLUME, WIDTH, small_limits())
            .unwrap()
            .with_no_slip_squeeze_flow(flow_limits())
            .unwrap();
    let state = source.initialize(GAP, true).unwrap();
    let cold = source.advance(state, 1.1 * GAP, 0.05, &gate).unwrap();
    let refresh = source
        .refresh_at_query_point(
            state,
            &card,
            interfaces.each_ref(),
            &QueryPoint::new().with("T", 320.0).unwrap(),
            &gate,
        )
        .unwrap();
    assert_eq!(refresh.source().squeeze_flow_limits(), Some(flow_limits()));
    assert_eq!(refresh.state().gap_m(), state.gap_m());
    let hot = refresh
        .source()
        .advance(refresh.state(), 1.1 * GAP, 0.05, &gate)
        .unwrap();
    near(
        hot.bulk_viscous_dissipation_j,
        2.0 * cold.bulk_viscous_dissipation_j,
        2e-13,
        1e-15,
    );
    assert_eq!(source.liquid().dynamic_viscosity_pa_s(), 1e-3);
    assert_eq!(refresh.source().liquid().dynamic_viscosity_pa_s(), 2e-3);
}

#[test]
fn g0_g4_pressure_geometry_rate_and_identity_refusals_leave_history_unchanged() {
    let m = viscous_model();
    let old = m.initialize(GAP, true).unwrap();
    let gate = CancelGate::new_clock_free();
    let valid = m.advance(old, 1.1 * GAP, 0.05, &gate).unwrap();
    assert_ne!(m.identity(), small_model().identity());
    assert_eq!(
        m.evaluate(small_model().initialize(GAP, true).unwrap()),
        Err(BridgeError::StateMismatch)
    );
    assert!(m.initialize(5e-5, true).is_err());
    for limits in [
        PlanarSqueezeFlowLimits {
            minimum_liquid_pressure_pa: 99_450.0,
            ..flow_limits()
        },
        PlanarSqueezeFlowLimits {
            maximum_reynolds_number: 1e-4,
            ..flow_limits()
        },
        PlanarSqueezeFlowLimits {
            maximum_unsteady_stokes_number: 1e-5,
            ..flow_limits()
        },
    ] {
        let narrow = small_model().with_no_slip_squeeze_flow(limits).unwrap();
        let state = narrow.initialize(GAP, true).unwrap();
        assert!(narrow.advance(state, 1.1 * GAP, 0.05, &gate).is_err());
        assert_eq!(state.gap_m(), GAP);
    }
    let compressed = small_model()
        .with_no_slip_squeeze_flow(PlanarSqueezeFlowLimits {
            maximum_liquid_pressure_pa: 100_010.0,
            ..flow_limits()
        })
        .unwrap();
    let state = compressed.initialize(1.1 * GAP, true).unwrap();
    assert!(compressed.advance(state, GAP, 0.01, &gate).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        m.advance(old, 1.1 * GAP, 0.05, &cancelled),
        Err(BridgeError::Cancelled)
    );
    assert_eq!(m.advance(old, 1.1 * GAP, 0.05, &gate).unwrap(), valid);
    assert!(
        small_model()
            .with_no_slip_squeeze_flow(PlanarSqueezeFlowLimits {
                maximum_gap_to_half_span_ratio: 0.2,
                ..flow_limits()
            })
            .is_err()
    );
}
