use super::*;

fn config() -> WallHeatConfig {
    WallHeatConfig { turbulent_prandtl: 0.9, relative_tolerance: 1e-9, max_refinements: 14 }
}

fn fluid(pr: f64) -> FluidProperties {
    FluidProperties { density_kg_m3: 1.0, specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0, kinematic_viscosity_m2_s: pr }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a / b - 1.0).abs() < tolerance, "{a:e} != {b:e}");
}

#[test]
fn rest_and_viscous_sublayer_recover_molecular_conduction() {
    let gate = CancelGate::new_clock_free();
    let air = FluidProperties::dry_air_300k();
    let rest = equilibrium_wall_heat_transfer(0.0, 0.003, &air, config(), &gate).unwrap();
    assert_eq!(rest.resistance_m2_k_w, 0.003 / air.conductivity_w_m_k);
    assert_eq!(rest.friction_velocity_m_s, 0.0);
    assert_eq!(rest.temperature_plus, 0.0);
    assert_eq!(rest.quadrature_evaluations, 0);
    let weak = equilibrium_wall_heat_transfer(7e-11, 1.0, &fluid(0.7), config(), &gate).unwrap();
    close(weak.resistance_m2_k_w, 1.0, 1e-12);
    close(weak.u_plus, weak.y_plus, 1e-12);
}

#[test]
fn equal_prandtl_numbers_obey_exact_reynolds_analogy() {
    let gate = CancelGate::new_clock_free();
    for re in [1e-8, 1.0, 100.0, 1e4, 1e8] {
        let wall = equilibrium_wall_heat_transfer(re * 0.9, 1.0, &fluid(0.9), config(), &gate).unwrap();
        close(wall.u_plus * wall.y_plus, re, 1e-11);
        assert_eq!(wall.temperature_plus, 0.9 * wall.u_plus);
        close(wall.resistance_m2_k_w, wall.u_plus / wall.y_plus, 1e-14);
        close(wall.friction_velocity_m_s, 0.9 * wall.y_plus, 1e-11);
        assert_eq!(wall.quadrature_evaluations, 0);
    }
}

#[test]
fn thermal_resistance_matches_independent_wall_normal_quadrature() {
    // Reference values use adaptive integration in log(1+y+) and a Brent
    // inversion of y+(u+) at every point, not the production u+ quadrature.
    let gate = CancelGate::new_clock_free();
    for (re, pr, expected) in [
        (100.0, 0.7, 0.8282457154560492),
        (100.0, 7.0, 0.5120197051613329),
        (1e4, 0.7, 0.0478636758450352),
        (1e4, 7.0, 0.0149917993638034),
    ] {
        let wall = equilibrium_wall_heat_transfer(re * pr, 1.0, &fluid(pr), config(), &gate).unwrap();
        close(wall.resistance_m2_k_w, expected, 2e-9);
        assert!(wall.resistance_m2_k_w > 0.0 && wall.resistance_m2_k_w < 1.0);
        assert!(wall.quadrature_relative_change <= config().relative_tolerance);
        assert!(wall.quadrature_evaluations <= (1 << config().max_refinements) + 1);
    }
}

#[test]
fn wall_law_respects_physical_unit_scaling() {
    let gate = CancelGate::new_clock_free();
    let a = equilibrium_wall_heat_transfer(7000.0, 1.0, &fluid(0.7), config(), &gate).unwrap();
    let scaled = FluidProperties { density_kg_m3: 2.0, specific_heat_j_kg_k: 4.0,
        kinematic_viscosity_m2_s: 15.0 * 0.7, conductivity_w_m_k: 120.0 };
    let b = equilibrium_wall_heat_transfer(21000.0, 5.0, &scaled, config(), &gate).unwrap();
    close(b.u_plus, a.u_plus, 1e-12);
    close(b.y_plus, a.y_plus, 1e-12);
    close(b.temperature_plus, a.temperature_plus, 1e-12);
    close(b.friction_velocity_m_s, 3.0 * a.friction_velocity_m_s, 1e-12);
    close(b.resistance_m2_k_w, a.resistance_m2_k_w / 24.0, 1e-12);
}

#[test]
fn invalid_physics_budget_and_cancellation_refuse() {
    let gate = CancelGate::new_clock_free();
    for speed in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(equilibrium_wall_heat_transfer(speed, 1.0, &fluid(0.7), config(), &gate).is_err());
    }
    for distance in [0.0, -1.0, f64::NAN] {
        assert!(equilibrium_wall_heat_transfer(1.0, distance, &fluid(0.7), config(), &gate).is_err());
    }
    let tight = WallHeatConfig { relative_tolerance: 1e-12, max_refinements: 3, ..config() };
    assert!(equilibrium_wall_heat_transfer(7e4, 1.0, &fluid(7.0), tight, &gate).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(equilibrium_wall_heat_transfer(1.0, 1.0, &fluid(0.7), config(), &cancelled),
        Err(ChtError::Cancelled)));
}
