//! Actual gas, material and port owners under an ideal fixed constraint.
use super::*;
use crate::cavity_phs::piston::tests::{controls, drive, gas};
use crate::cavity_phs::piston::double_acting::DoubleActingGasPiston;
use fs_ports::{ModuleWindowAudit, PortedModule};

#[test]
fn heating_changes_pressure_transport_and_sound_at_exactly_fixed_volume() {
    let initial = gas(0.75, 0.001);
    let mut expected = initial.clone();
    expected.add_heat(20.0, controls().gas, || false).unwrap();
    let input = PistonInput { duration_s: 1.0, external_force_n: 3.0, ..drive(20.0) };
    let fixed = PistonSolve { volume_bracket_m3: [initial.volume_m3(); 2], ..controls() };
    let mut endpoints = Vec::new();
    for mass in [0.125, 800.0] {
        let mut wall = GasPiston::new(initial.clone(), mass, 0.01, 0.0).unwrap()
            .with_fixed_support().unwrap();
        let report = wall.step(input, fixed, || false).unwrap();
        assert_eq!(wall.support(), PistonSupport::Fixed);
        assert_eq!(wall.mass_kg(), mass);
        assert_eq!(wall.gas().volume_m3().to_bits(), initial.volume_m3().to_bits());
        assert_eq!(wall.velocity_m_s(), 0.0);
        assert_eq!(report.velocity_m_s, [0.0; 2]);
        assert_eq!(report.displacement_m, 0.0);
        assert_eq!(report.kinetic_energy_change_j, 0.0);
        assert_eq!(report.gas.work_by_gas_j, 0.0);
        assert_eq!(report.support_work_j, 0.0);
        assert_eq!(report.external_force_work_j, 0.0);
        assert_eq!(report.exterior_pressure_work_j, 0.0);
        let applied = input.duration_s * (wall.area_m2()
            * (report.midpoint.state().pressure-input.exterior_pressure_pa)+input.external_force_n);
        assert_eq!(report.support_impulse_n_s, -applied);
        assert_eq!(report.momentum_residual_n_s, 0.0);
        assert!(report.total_energy_residual_j.abs() <= fixed.energy_tolerance_j);
        assert!(report.mechanical_energy_residual_j.abs() <= fixed.energy_tolerance_j);
        assert_eq!(report.iterations, 0);
        assert_eq!(report.material_trials, 1);
        assert!((wall.gas().state().state().temperature-expected.state().state().temperature).abs() < 1e-7);
        assert!(wall.gas().state().state().pressure > initial.state().state().pressure);
        assert!(wall.gas().state().state().sound_speed > initial.state().state().sound_speed);
        assert_ne!(wall.gas().state().state().dynamic_viscosity, initial.state().state().dynamic_viscosity);
        assert_eq!(wall.gas().point().axes().get("specimen-lot"), Some(&42.0));
        assert!(!matches!(wall.audit_last_window().unwrap(), ModuleWindowAudit::Unavailable { .. }));
        endpoints.push(wall.gas().clone());
    }
    // A stationary support is independent of the specimen's finite inertia.
    assert_eq!(endpoints[0], endpoints[1]);
    let mut free = GasPiston::new(initial, 0.125, 0.01, 0.0).unwrap();
    assert_eq!(free.support(), PistonSupport::Free);
    assert!(matches!(free.step(input, controls(), || false), Err(PistonError::TimeResolution(_))));
}

#[test]
fn attaching_a_support_cannot_discard_motion_or_hide_an_incompatible_linkage() {
    let moving = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 2.0).unwrap();
    assert!(matches!(moving.clone().with_fixed_support(), Err(PistonError::Invalid(_))));
    assert_eq!(moving.velocity_m_s(), 2.0);
    assert_eq!(moving.support(), PistonSupport::Free);
    let fixed = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap()
        .with_fixed_support().unwrap();
    assert!(matches!(DoubleActingGasPiston::new(fixed.clone(), gas(0.25, 0.001), 0.01),
        Err(PistonError::Invalid(_))));
    assert_eq!(fixed.support(), PistonSupport::Fixed);
}

#[test]
fn failed_sources_and_final_cancellation_leave_the_supported_owner_unchanged() {
    let original = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap()
        .with_fixed_support().unwrap();
    let mut count = 0;
    original.clone().step(drive(20.0), controls(), || { count += 1; false }).unwrap();
    for stop in [1, count] {
        let mut wall = original.clone();
        let mut polls = 0;
        assert_eq!(wall.step(drive(20.0), controls(), || { polls += 1; polls == stop }),
            Err(PistonError::Cancelled));
        assert_eq!(wall, original);
        assert_eq!(wall.ports, original.ports);
    }
    let mut wall = original.clone();
    assert!(wall.step(drive(1e20), controls(), || false).is_err());
    assert_eq!(wall, original);
    assert_eq!(wall.ports, original.ports);
}
