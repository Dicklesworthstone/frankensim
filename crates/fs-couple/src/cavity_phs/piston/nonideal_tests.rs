//! Dense manufactured PR mixture through the original native piston mechanics.
//! Source cards are synthetic; no empirical real-fluid validation is asserted.

use super::*;
use fs_matdb::QueryPoint;
use fs_material::thermochemical_gas::nonideal::energy::parcel::IsochoricRealGasParcel;
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
};
use fs_ports::{ModuleWindowAudit, PortedModule};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Pressure, QuantitySpec, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseV1;

pub(super) fn parcel() -> IsochoricRealGasParcel {
    // Reuse the original piston fixture's NASA and source-backed transport.
    let ideal = super::tests::gas(0.75, 0.001);
    let model = ThermochemicalPengRobinsonModel::new(
        ideal.model().clone(),
        vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A").unwrap(),
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B").unwrap(),
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6),
                acentric_factor: -0.00219,
            },
        ],
        vec![0.0, 0.03, 0.03, 0.0],
    )
    .unwrap();
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            300.0,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 5.0e6)
        .unwrap()
        .with("specimen-lot", 42.0)
        .unwrap();
    IsochoricRealGasParcel::new(model, point, 0.001, PengRobinsonPhaseV1::Vapor).unwrap()
}

pub(super) fn controls() -> PistonSolve {
    PistonSolve {
        volume_bracket_m3: [0.00095, 0.00105],
        gas: GasEnergySolve {
            minimum_temperature_k: 280.0,
            maximum_temperature_k: 450.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
        momentum_tolerance_n_s: 1e-9,
        energy_tolerance_j: 5e-7,
        maximum_iterations: 128,
    }
}

pub(super) fn drive(heat_j: f64) -> PistonInput {
    PistonInput {
        duration_s: 2e-5,
        heat_j,
        exterior_pressure_pa: 5e6,
        external_force_n: 50.0,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn dense_pr_gas_drives_native_motion_and_the_same_complete_energy_account() {
    let gas = parcel();
    let ideal = ClosedGasParcel::new(
        gas.model().source().clone(),
        gas.point().clone(),
        gas.volume_m3(),
    )
    .unwrap();
    assert!((gas.state().density_kg_m3() / ideal.state().state().density - 1.0).abs() > 0.05);
    let mut piston = GasPiston::new(gas.clone(), 2.0, 0.01, 0.03).unwrap();
    let mut ideal_piston = GasPiston::new(ideal, 2.0, 0.01, 0.03).unwrap();
    let initial_energy = gas.state().specific_internal_energy_j_kg();
    let initial_kinetic = 0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2);
    let (mut heat, mut exterior_work, mut force_work) = (0.0, 0.0, 0.0);
    for _ in 0..8 {
        let old_volume = piston.gas().volume_m3();
        let input = drive(25.0);
        let step = piston.step(input, controls(), || false).unwrap();
        ideal_piston.step(input, controls(), || false).unwrap();
        let dv = piston.gas().volume_m3() - old_volume;
        let force = piston.area_m2() * (step.midpoint.pressure_pa() - input.exterior_pressure_pa)
            + input.external_force_n;
        near(
            piston.mass_kg() * (step.velocity_m_s[1] - step.velocity_m_s[0]),
            input.duration_s * force,
            controls().momentum_tolerance_n_s,
        );
        near(
            dv,
            piston.area_m2()
                * input.duration_s
                * f64::midpoint(step.velocity_m_s[0], step.velocity_m_s[1]),
            2e-19,
        );
        near(
            step.gas.work_by_gas_j,
            step.midpoint.pressure_pa() * dv,
            1e-12,
        );
        assert!(step.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
        assert_eq!(piston.gas().mass_kg(), gas.mass_kg());
        assert_eq!(piston.gas().point().axes().get("specimen-lot"), Some(&42.0));
        assert_eq!(step.midpoint.source().components().len(), 2);
        assert_eq!(piston.last_window_duration_s(), Some(input.duration_s));
        assert!(piston.last_window_transfers().is_some());
        assert!(matches!(
            piston.audit_last_window().unwrap(),
            ModuleWindowAudit::Available(_)
        ));
        heat += input.heat_j;
        exterior_work += step.exterior_pressure_work_j;
        force_work += step.external_force_work_j;
    }
    // Independent accepted-state energy difference; this manufactured source
    // has no large common formation offset, so direct subtraction is resolved.
    let stored =
        gas.mass_kg() * (piston.gas().state().specific_internal_energy_j_kg() - initial_energy);
    let kinetic = 0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2) - initial_kinetic;
    near(
        stored + kinetic + exterior_work - force_work,
        heat,
        8.0 * controls().energy_tolerance_j,
    );
    assert!(piston.gas().volume_m3() > gas.volume_m3());
    assert!((piston.velocity_m_s() - ideal_piston.velocity_m_s()).abs() > 1e-8);
    assert_ne!(
        piston.gas().state().pressure_pa(),
        ideal_piston.gas().state().state().pressure
    );
}

#[test]
fn fixed_real_gas_wall_heats_with_native_support_impulse_and_zero_mechanical_work() {
    let gas = parcel();
    let mut direct = gas.clone();
    let mut piston = GasPiston::new(gas.clone(), 2.0, 0.01, 0.0)
        .unwrap()
        .with_fixed_support()
        .unwrap();
    let input = PistonInput {
        duration_s: 0.5,
        ..drive(200.0)
    };
    let fixed = PistonSolve {
        volume_bracket_m3: [gas.volume_m3(); 2],
        ..controls()
    };
    let step = piston.step(input, fixed, || false).unwrap();
    direct.add_heat(input.heat_j, fixed.gas, || false).unwrap();
    assert_eq!(piston.support(), PistonSupport::Fixed);
    assert_eq!(piston.velocity_m_s(), 0.0);
    assert_eq!(piston.gas().volume_m3(), gas.volume_m3());
    assert_eq!(step.displacement_m, 0.0);
    assert_eq!(step.gas.work_by_gas_j, 0.0);
    assert_eq!(step.support_work_j, 0.0);
    near(
        step.gas.stored_energy_change_j,
        input.heat_j,
        fixed.energy_tolerance_j,
    );
    near(
        piston.gas().state().temperature_k(),
        direct.state().temperature_k(),
        1e-8,
    );
    near(
        piston.gas().state().pressure_pa(),
        direct.state().pressure_pa(),
        1e-3,
    );
    let force = piston.area_m2() * (step.midpoint.pressure_pa() - input.exterior_pressure_pa)
        + input.external_force_n;
    assert_eq!(step.support_impulse_n_s, -input.duration_s * force);
    assert!(step.support_impulse_n_s < 0.0);
    assert!(piston.gas().state().pressure_pa() > step.midpoint.pressure_pa());
    assert_eq!(step.material_trials, 1);
    assert!(matches!(
        piston.audit_last_window().unwrap(),
        ModuleWindowAudit::Available(_)
    ));
}

#[test]
fn real_gas_acoustic_clock_uses_original_eos_compressibility_not_the_ideal_reference() {
    let gas = parcel();
    let mut piston = GasPiston::new(gas.clone(), 2.0, 0.01, 0.0).unwrap();
    let state = gas.state();
    let rho = state.density_kg_m3();
    let c = state.sound_speed_m_s();
    let omega = (rho * c * c * piston.area_m2() * piston.area_m2()
        / (piston.mass_kg() * gas.volume_m3()))
    .sqrt();
    let ideal = state.source().state();
    let wrong = (ideal.density
        * ideal.sound_speed
        * ideal.sound_speed
        * piston.area_m2()
        * piston.area_m2()
        / (piston.mass_kg() * gas.volume_m3()))
    .sqrt();
    assert!((omega / wrong - 1.0).abs() > 1e-3);
    let initial = piston.clone();
    let input = PistonInput {
        duration_s: 0.041 / omega,
        heat_j: 0.0,
        exterior_pressure_pa: state.pressure_pa(),
        external_force_n: 0.0,
    };
    let error = piston.step(input, controls(), || false).unwrap_err();
    match error {
        PistonError::TimeResolution(fs_phs::linear_exact::ResolvedBandError::Outside {
            rate_rad_s,
            ..
        }) => {
            assert_eq!(rate_rad_s, omega);
        }
        other => panic!("expected original acoustic clock refusal, got {other:?}"),
    }
    assert_eq!(piston, initial);
    assert_eq!(
        piston.audit_last_window().unwrap(),
        ModuleWindowAudit::NotStepped
    );
    let step = piston
        .step(
            PistonInput {
                duration_s: 0.039 / omega,
                ..input
            },
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(piston, initial);
    assert_eq!(step.total_energy_residual_j, 0.0);
}

#[test]
fn real_gas_source_failure_and_final_cancellation_preserve_gas_wall_and_port_window() {
    let mut piston = GasPiston::new(parcel(), 2.0, 0.01, 0.0).unwrap();
    piston.step(drive(10.0), controls(), || false).unwrap();
    assert!(piston.last_window_transfers().is_some());
    let before = piston.clone();
    let window = |owner: &GasPiston<IsochoricRealGasParcel>| {
        (
            owner.last_window_duration_s(),
            owner.last_window_transfers().map(<[_]>::to_vec),
            owner.last_window_transfer_error().cloned(),
            owner.audit_last_window().unwrap(),
        )
    };
    let before_window = window(&before);
    assert!(matches!(
        piston.step(drive(1e12), controls(), || false),
        Err(PistonError::Gas(_))
    ));
    assert_eq!(piston, before);
    assert_eq!(window(&piston), before_window);
    let mut probe = before.clone();
    let mut calls = 0;
    let expected = probe
        .step(drive(25.0), controls(), || {
            calls += 1;
            false
        })
        .unwrap();
    assert!(calls > 2);
    let mut observed = 0;
    let refused = piston.step(drive(25.0), controls(), || {
        observed += 1;
        observed == calls
    });
    assert_eq!(refused, Err(PistonError::Cancelled));
    assert_eq!(piston, before);
    assert_eq!(window(&piston), before_window);
    let retried = piston.step(drive(25.0), controls(), || false).unwrap();
    assert_eq!(retried, expected);
    assert_eq!(piston, probe);
    assert_eq!(window(&piston), window(&probe));
}
