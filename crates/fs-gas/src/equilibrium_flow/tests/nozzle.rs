//! G1 pressure-driven reacting flow against analytic and mass-action oracles.
use super::super::nozzle::{
    EquilibriumNozzleControl, EquilibriumNozzleError, equilibrium_nozzle_flow,
};
use super::*;
use fs_qty::Area;

fn nozzle_controls() -> EquilibriumNozzleControl {
    let mut c = EquilibriumNozzleControl::new(controls().isentropic, Pressure::new(100.0));
    c.relative_sonic_tolerance = 1e-10;
    c.relative_pressure_tolerance = 1e-10;
    c
}

#[test]
fn g1_frozen_limit_resolves_zero_subsonic_and_choked_mass_flow() {
    let chemistry = single();
    let feed = [(id("A"), Amount::new(2.0))];
    let gamma = 4.0_f64 / 3.0;
    let gas_constant = R / 0.01;
    let critical_ratio = (2.0 / (gamma + 1.0)).powf(gamma / (gamma - 1.0));
    for ratio in [
        1.0,
        0.98,
        0.75,
        critical_ratio * (1.0 + 1e-5),
        critical_ratio * (1.0 - 1e-5),
        0.1,
        0.0,
    ] {
        let result = equilibrium_nozzle_flow(
            &chemistry,
            Temperature::new(T0),
            Pressure::new(P0),
            &feed,
            Pressure::new(P0 * ratio),
            Area::new(0.002),
            nozzle_controls(),
            || false,
        )
        .unwrap();
        let pressure = P0 * ratio.max(critical_ratio);
        let t = T0 * (pressure / P0).powf(0.25);
        let velocity = (2.0 * 4.0 * gas_constant * (T0 - t)).sqrt();
        let density = pressure / (gas_constant * t);
        if ratio == 1.0 {
            assert_eq!(result.outlet, result.inlet);
            assert_eq!(result.mass_flow.value(), 0.0);
            assert_eq!(result.velocity_m_s, 0.0);
            assert_eq!(result.pressure_evaluations, 1);
        } else {
            close(
                result.outlet.equilibrium().pressure().value(),
                pressure,
                3e-5,
            );
            close(result.outlet.equilibrium().temperature().value(), t, 2e-7);
            close(result.velocity_m_s, velocity, 2e-6);
            close(result.mass_flow.value(), 0.002 * density * velocity, 2e-9);
            close(
                result.sample.response.sound_speed_m_s.powi(2),
                gamma * gas_constant * t,
                3e-4,
            );
        }
        assert_eq!(result.choked, ratio <= critical_ratio);
        assert_eq!(result.sonic_pressure_bracket_pa.is_some(), result.choked);
        assert!(result.stagnation_enthalpy_residual_j_kg.abs() < 1e-7);
        assert!(
            result.sample.entropy_change_j_k.abs()
                <= nozzle_controls().isentropic.entropy_tolerance_j_k
        );
        assert!(result.outlet.equilibrium().is_bound_to(&chemistry));
        if result.choked {
            close(result.mach_number, 1.0, 2e-9);
        }
    }
}

// Independent A2 <=> 2A mass action gives all local derivatives explicitly.
// These formulas do not call the native equilibrium response or nozzle root.
fn reacting_station(pressure: f64) -> (f64, f64, f64, f64) {
    let t = bisect(dissociation_s(T0, P0), |t| dissociation_s(t, pressure));
    let a = alpha(t, pressure);
    let amount = 1.0 + a;
    let density = 0.002 * pressure / (amount * R * t);
    let cp = R * (3.5 + 1.5 * a + 1.125 * a * (1.0 - a * a));
    let expansion = 1.0 + 0.75 * a * (1.0 - a);
    let compressibility = 1.0 + 0.5 * a * (1.0 - a);
    let sound_squared =
        (pressure / density) / (compressibility - amount * R * expansion * expansion / cp);
    let speed_squared = 2.0 * (dissociation_h(T0, P0) - dissociation_h(t, pressure)) / 0.002;
    (t, density, speed_squared, sound_squared)
}

fn reacting_critical_pressure() -> f64 {
    let (mut lo, mut hi) = (0.1 * P0, P0);
    for _ in 0..70 {
        let middle = 0.5 * (lo + hi);
        let (_, _, v2, c2) = reacting_station(middle);
        if c2 - v2 < 0.0 {
            lo = middle;
        } else {
            hi = middle;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn g1_reacting_nozzle_chokes_at_the_mass_action_maximum_and_preserves_sources() {
    let chemistry = dissociation();
    let feed = [(id("A2"), Amount::new(1.0))];
    let critical = reacting_critical_pressure();
    let (_, rho_critical, v2_critical, _) = reacting_station(critical);
    let maximum_flux = rho_critical * v2_critical.sqrt();
    for factor in [0.99, 1.01] {
        let (_, rho, v2, _) = reacting_station(factor * critical);
        assert!(rho * v2.sqrt() < maximum_flux);
    }
    for back in [0.9 * P0, 1.01 * critical, 0.1 * P0, 0.0] {
        let result = equilibrium_nozzle_flow(
            &chemistry,
            Temperature::new(T0),
            Pressure::new(P0),
            &feed,
            Pressure::new(back),
            Area::new(0.003),
            nozzle_controls(),
            || false,
        )
        .unwrap();
        let expected_pressure = back.max(critical);
        let (t, density, v2, c2) = reacting_station(expected_pressure);
        close(
            result.outlet.equilibrium().pressure().value(),
            expected_pressure,
            4e-5,
        );
        close(result.outlet.equilibrium().temperature().value(), t, 2e-7);
        close(result.sample.response.density_kg_m3, density, 2e-10);
        close(result.sample.response.sound_speed_m_s.powi(2), c2, 5e-3);
        close(result.mass_flow.value(), 0.003 * density * v2.sqrt(), 2e-9);
        assert_eq!(result.choked, back <= critical);
        assert!(
            result.sample.entropy_change_j_k.abs()
                <= nozzle_controls().isentropic.entropy_tolerance_j_k
        );
        assert!(
            result.sample.relative_mass_residual.abs()
                <= controls().isentropic.equilibrium.relative_balance_tolerance
        );
        assert!(result.stagnation_enthalpy_residual_j_kg.abs() <= 2e-8);
        check_products(&result.inlet, &chemistry);
        check_products(&result.outlet, &chemistry);
        if result.choked {
            close(result.mach_number, 1.0, 2e-9);
        }
    }
    // A frozen inlet mixture cannot supply this reacting critical pressure.
    let a0 = alpha(T0, P0);
    let frozen_gamma = (3.5 + 1.5 * a0) / (3.5 + 1.5 * a0 - 1.0 - a0);
    let frozen_critical =
        P0 * (2.0 / (frozen_gamma + 1.0)).powf(frozen_gamma / (frozen_gamma - 1.0));
    assert!((critical - frozen_critical).abs() > 500.0);
}

#[test]
fn g3_feed_normalization_area_scaling_and_lower_back_pressure_preserve_the_same_law() {
    let chemistry = dissociation();
    let mut reference = None;
    for scale in [0.125, 1.0, 8.0] {
        let feed = [(id("A2"), Amount::new(scale))];
        let mut c = nozzle_controls();
        c.isentropic.entropy_tolerance_j_k *= scale;
        for area in [0.002, 0.004] {
            let flow = equilibrium_nozzle_flow(
                &chemistry,
                Temperature::new(T0),
                Pressure::new(P0),
                &feed,
                Pressure::new(0.0),
                Area::new(area),
                c,
                || false,
            )
            .unwrap();
            let normalized = flow.mass_flow.value() / area;
            if let Some(expected) = reference {
                close(normalized, expected, 5e-7);
            } else {
                reference = Some(normalized);
            }
            close(flow.feed_mass_kg, 0.002 * scale, 1e-16);
            let low_back = equilibrium_nozzle_flow(
                &chemistry,
                Temperature::new(T0),
                Pressure::new(P0),
                &feed,
                Pressure::new(0.1 * P0),
                Area::new(area),
                c,
                || false,
            )
            .unwrap();
            assert_eq!(flow.outlet, low_back.outlet);
            assert_eq!(flow.mass_flow, low_back.mass_flow);
        }
    }
}

#[test]
fn g3_weak_drive_retains_source_relative_energy_under_large_reference_datums() {
    let back = P0 * (1.0 - 1e-7);
    let ratio = back / P0;
    let gas_constant = R / 0.01;
    let expected_v2 = -2.0 * 4.0 * gas_constant * T0 * (0.25 * ratio.ln()).exp_m1();
    let expected_t = T0 * ratio.powf(0.25);
    let expected_flux = back / (gas_constant * expected_t) * expected_v2.sqrt();
    for datum in [-1e11, 0.0, 1e11] {
        let chemistry = owner(vec![model("A", 0.01, 4.0, datum, datum / T0)], vec![1]);
        let flow = equilibrium_nozzle_flow(
            &chemistry,
            Temperature::new(T0),
            Pressure::new(P0),
            &[(id("A"), Amount::new(1.0))],
            Pressure::new(back),
            Area::new(0.01),
            nozzle_controls(),
            || false,
        )
        .unwrap();
        assert!(!flow.choked);
        assert!(flow.mass_flow.value() > 0.0);
        assert_eq!(flow.pressure_evaluations, 2);
        close(flow.mass_flux_kg_m2_s, expected_flux, expected_flux * 2e-7);
        assert!(flow.stagnation_enthalpy_residual_j_kg.abs() < 1e-12);
    }
}

#[test]
fn g0_closed_restrictions_and_shallow_source_paths_do_not_query_an_unused_throat() {
    let chemistry = single();
    let feed = [(id("A"), Amount::new(1.0))];
    let mut c = nozzle_controls();
    c.maximum_pressure_evaluations = 1;
    for (back, area) in [(0.0, 0.0), (P0, 0.1)] {
        let flow = equilibrium_nozzle_flow(
            &chemistry,
            Temperature::new(T0),
            Pressure::new(P0),
            &feed,
            Pressure::new(back),
            Area::new(area),
            c,
            || false,
        )
        .unwrap();
        assert_eq!(flow.inlet, flow.outlet);
        assert_eq!(flow.mass_flow.value(), 0.0);
        assert_eq!(flow.sample.specific_enthalpy_change_j_kg, 0.0);
        assert_eq!(flow.sample.entropy_change_j_k, 0.0);
        assert_eq!(flow.pressure_evaluations, 1);
        assert_eq!(flow.sonic_residual_m2_s2, None);
    }
    c = nozzle_controls();
    c.isentropic.minimum_temperature_k = 950.0;
    c.isentropic.maximum_temperature_k = 1100.0;
    let shallow = equilibrium_nozzle_flow(
        &chemistry,
        Temperature::new(T0),
        Pressure::new(P0),
        &feed,
        Pressure::new(0.99 * P0),
        Area::new(0.01),
        c,
        || false,
    )
    .unwrap();
    assert!(!shallow.choked);
    assert_eq!(shallow.pressure_evaluations, 2);
    assert!(matches!(
        equilibrium_nozzle_flow(
            &chemistry,
            Temperature::new(T0),
            Pressure::new(P0),
            &feed,
            Pressure::new(0.0),
            Area::new(0.01),
            c,
            || false,
        ),
        Err(EquilibriumNozzleError::Chemistry { .. })
    ));
}

#[test]
fn g4_original_pressure_work_and_unsupported_search_floor_refuse_without_invented_flow() {
    let chemistry = dissociation();
    let feed = [(id("A2"), Amount::new(1.0))];
    let call = |back, area, c| {
        equilibrium_nozzle_flow(
            &chemistry,
            Temperature::new(T0),
            Pressure::new(P0),
            &feed,
            Pressure::new(back),
            Area::new(area),
            c,
            || false,
        )
    };
    let mut c = nozzle_controls();
    c.maximum_pressure_evaluations = 1;
    assert!(matches!(
        call(0.0, 0.01, c),
        Err(EquilibriumNozzleError::WorkLimit {
            pressure_evaluations: 1
        })
    ));
    c = nozzle_controls();
    c.minimum_pressure_pa = 0.9 * P0;
    assert!(matches!(
        call(0.0, 0.01, c),
        Err(EquilibriumNozzleError::OutsidePressureBracket { .. })
    ));
    assert!(matches!(
        call(1.01 * P0, 0.01, nozzle_controls()),
        Err(EquilibriumNozzleError::Invalid(_))
    ));
    assert!(matches!(
        call(0.0, -1.0, nozzle_controls()),
        Err(EquilibriumNozzleError::Invalid(_))
    ));
    c = nozzle_controls();
    c.isentropic.maximum_temperature_iterations = 1;
    assert!(matches!(
        call(0.7 * P0, 0.01, c),
        Err(EquilibriumNozzleError::Chemistry { .. })
    ));
}

#[test]
fn g4_final_nozzle_cancellation_is_atomic_and_retry_is_deterministic() {
    let chemistry = dissociation();
    let feed = [(id("A2"), Amount::new(1.0))];
    let before = feed.clone();
    let call = |cancelled: &mut dyn FnMut() -> bool| {
        equilibrium_nozzle_flow(
            &chemistry,
            Temperature::new(T0),
            Pressure::new(P0),
            &feed,
            Pressure::new(0.0),
            Area::new(0.001),
            nozzle_controls(),
            cancelled,
        )
    };
    let mut total = 0;
    let expected = call(&mut || {
        total += 1;
        false
    })
    .unwrap();
    assert!(total > 30);
    for stop in [1, total / 2, total] {
        let mut polls = 0;
        assert!(matches!(
            call(&mut || {
                polls += 1;
                polls == stop
            }),
            Err(EquilibriumNozzleError::Cancelled)
        ));
        assert_eq!(polls, stop);
    }
    assert_eq!(call(&mut || false).unwrap(), expected);
    assert_eq!(feed, before);
}
