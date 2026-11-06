//! G1/G3: original source gas, native wall motion, energy and entropy together.

use super::*;
use super::super::tests::{controls, drive, gas};
use fs_ports::PortedModule;

fn entropy_change(initial: &ClosedGasParcel, final_gas: &ClosedGasParcel) -> f64 {
    initial.mass_kg() * initial.model()
        .specific_entropy_change_j_kg_k(initial.point(), final_gas.point())
        .unwrap()
}

#[test]
fn g1_adiabatic_motion_preserves_source_entropy_and_combined_energy() {
    for velocity in [-0.2, 0.2] {
        // Temperature-dependent heat capacity and a nontrivial mixture.
        let mut piston = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, velocity).unwrap();
        let initial = piston.clone();
        let mut external_work = 0.0;
        let steps = 48;
        for index in 0..steps {
            let input = PistonInput {
                external_force_n: if index < 20 { -0.5 } else { 0.0 },
                ..drive(0.0)
            };
            let report = piston.step_discrete_gradient(input, controls(), || false).unwrap();
            external_work += report.external_force_work_j - report.exterior_pressure_work_j;
            assert!(report.momentum_residual_n_s.abs() <= controls().momentum_tolerance_n_s);
            assert!(report.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
            assert!(report.gas.entropy_residual_j_k.abs()
                <= 8.0 * initial.gas().mass_kg() * controls().gas.energy_tolerance_j_kg
                    / controls().gas.minimum_temperature_k);
        }
        let stored = initial.gas().mass_kg()
            * initial.gas().model().thermodynamics()
                .specific_internal_energy_change(
                    fs_qty::Temperature::new(initial.gas().state().state().temperature),
                    fs_qty::Temperature::new(piston.gas().state().state().temperature),
                ).unwrap().value();
        let kinetic = 0.5 * initial.mass_kg()
            * (piston.velocity_m_s().powi(2) - initial.velocity_m_s().powi(2));
        assert!((stored + kinetic - external_work).abs()
            <= f64::from(steps) * controls().energy_tolerance_j);
        let entropy_bound = 8.0 * f64::from(steps) * initial.gas().mass_kg()
            * controls().gas.energy_tolerance_j_kg / controls().gas.minimum_temperature_k;
        assert!(entropy_change(initial.gas(), piston.gas()).abs() <= entropy_bound);
        assert_eq!(piston.gas().mass_kg(), initial.gas().mass_kg());
        assert_eq!(piston.gas().point().axes().get("specimen-lot"), Some(&42.0));
    }
}

#[test]
fn g1_constant_capacity_adiabatic_endpoint_matches_the_analytic_isentrope() {
    for velocity in [-8.0, 8.0] {
        let mut piston = GasPiston::new(gas(0.75, 0.0), 1.0, 0.01, velocity).unwrap();
        let initial = piston.gas().clone();
        let report = piston.step_discrete_gradient(drive(0.0), controls(), || false).unwrap();
        let properties = initial.state().thermodynamics().properties();
        let exponent = properties.mass_specific_gas_constant().value()
            / properties.mass_specific_cv().value();
        let exact_temperature = initial.state().state().temperature
            * (initial.volume_m3() / piston.gas().volume_m3()).powf(exponent);
        assert!((piston.gas().state().state().temperature / exact_temperature - 1.0).abs() < 2e-12);
        assert_ne!(report.displacement_m, 0.0);
        let pressure = report.thermodynamic_mean.state().pressure;
        assert!((report.gas.work_by_gas_j
            - pressure * (piston.gas().volume_m3() - initial.volume_m3())).abs() < 1e-10);
        let impulse = drive(0.0).duration_s
            * piston.area_m2() * (pressure - drive(0.0).exterior_pressure_pa);
        assert!((piston.mass_kg() * (piston.velocity_m_s() - velocity) - impulse).abs()
            <= controls().momentum_tolerance_n_s);
    }
}

#[test]
fn g1_heat_port_uses_the_material_entropy_conjugate_temperature() {
    let mut piston = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap();
    let initial = piston.gas().clone();
    let input = drive(20.0);
    let report = piston.step_discrete_gradient(input, controls(), || false).unwrap();
    assert!(piston.velocity_m_s() > 0.0);
    let heat_entropy = input.heat_j / report.thermodynamic_mean.state().temperature;
    let actual = entropy_change(&initial, piston.gas());
    assert!((actual - heat_entropy).abs() <= 1e-11);
    assert_eq!(report.gas.entropy_transfer_j_k, heat_entropy);
    assert!(report.gas.gibbs_residual_j.abs() <= controls().energy_tolerance_j);
    assert!(piston.audit_last_window().is_ok());
    assert!(piston.last_window_transfers().is_some());
}

#[test]
fn g1_fixed_wall_matches_original_isochoric_heat_with_zero_work() {
    let initial = gas(0.75, 0.001);
    let mut reference = initial.clone();
    reference.add_heat(50.0, controls().gas, || false).unwrap();
    let mut piston = GasPiston::new(initial.clone(), 1.0, 0.01, 0.0).unwrap()
        .with_fixed_support().unwrap();
    let config = PistonSolve {
        volume_bracket_m3: [initial.volume_m3(); 2],
        ..controls()
    };
    let report = piston.step_discrete_gradient(drive(50.0), config, || false).unwrap();
    assert_eq!(piston.gas().volume_m3(), initial.volume_m3());
    assert_eq!(piston.velocity_m_s(), 0.0);
    assert_eq!(report.gas.work_by_gas_j, 0.0);
    assert_eq!(report.support_work_j, 0.0);
    assert_eq!(report.material_trials, 1);
    assert_eq!(report.iterations, 0);
    assert!((piston.gas().state().state().temperature
        - reference.state().state().temperature).abs() < 1e-8);
    assert!((report.gas.stored_energy_change_j - 50.0).abs()
        <= controls().energy_tolerance_j);
    assert!(report.support_impulse_n_s < 0.0);
}

#[test]
fn g3_stationary_state_and_failure_preserve_accepted_physics_and_ports() {
    let mut piston = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap();
    let untouched = piston.clone();
    let stationary = piston.step_discrete_gradient(drive(0.0), controls(), || false).unwrap();
    assert_eq!(piston, untouched);
    assert_eq!(stationary.total_energy_residual_j, 0.0);
    assert_eq!(stationary.gas.entropy_change_j_k, 0.0);
    let before = piston.clone();
    let transfers = piston.last_window_transfers().map(<[_]>::to_vec);
    let duration = piston.last_window_duration_s();
    for input in [
        PistonInput { duration_s: 0.0, ..drive(1.0) },
        PistonInput { exterior_pressure_pa: -1.0, ..drive(1.0) },
        drive(1e9),
    ] {
        assert!(piston.step_discrete_gradient(input, controls(), || false).is_err());
        assert_eq!(piston, before);
        assert_eq!(piston.last_window_duration_s(), duration);
        assert_eq!(piston.last_window_transfers().map(<[_]>::to_vec), transfers);
    }
    let mut calls = 0;
    let mut probe = piston.clone();
    probe.step_discrete_gradient(drive(20.0), controls(), || {
        calls += 1;
        false
    }).unwrap();
    let mut current = 0;
    assert!(matches!(
        piston.step_discrete_gradient(drive(20.0), controls(), || {
            current += 1;
            current == calls
        }),
        Err(PistonError::Cancelled)
    ));
    assert_eq!(piston, before);
    assert_eq!(piston.last_window_duration_s(), duration);
    assert_eq!(piston.last_window_transfers().map(<[_]>::to_vec), transfers);
}

mod real_gas;
