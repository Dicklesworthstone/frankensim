//! Source-preserving continuation through the original dense-gas power solve.

use super::*;
use crate::cavity_phs::piston::material::tests::body as specimen;
use crate::cavity_phs::piston::nonideal_tests::{controls, drive, parcel};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::thermochemical_gas::nonideal::ResolvedThermochemicalRealGasState;
use fs_ports::PortedModule;
use std::convert::Infallible;

fn input() -> PistonInput {
    PistonInput {
        duration_s: 2.0e-6,
        exterior_pressure_pa: 5.2e6,
        external_force_n: 0.0,
        ..drive(0.0)
    }
}
fn exchange() -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-0.1, 0.1],
        source_tolerance_j: 1.0e-9,
        energy_tolerance_j: 1.0e-6,
        maximum_iterations: 128,
    }
}
fn curve() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        fs_blake3::hash_bytes(b"synthetic finite thermal store for sourced piston power"),
        SolidLiquidPhase::Solid,
        [200.0, 800.0]
            .into_iter()
            .map(|t| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 500.0 * (t - 200.0),
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 4000.0,
            })
            .collect(),
    )
    .unwrap()
}
fn thermal_body(curve: &EquilibriumEnthalpyPhaseCurve) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve,
        initial: curve.state_at_specific_enthalpy(150000.0).unwrap(),
        mass_kg: 0.001,
        volume_m3: 0.001 / 4000.0,
        surface_area_m2: 0.02,
        internal_heat_j: 0.0,
        duration_s: input().duration_s,
    }
}
fn conductance(
    gas: &ResolvedThermochemicalRealGasState,
    wall: EquilibriumPhaseState,
) -> Result<ConductanceSample<(f64, f64)>, Infallible> {
    Ok(ConductanceSample {
        conductance_w_k: 20.0,
        source: (gas.temperature_k(), wall.temperature_k()),
    })
}
fn power(
    gas: &ResolvedThermochemicalRealGasState,
    wall: EquilibriumPhaseState,
) -> Result<PowerSample<(f64, f64)>, Infallible> {
    Ok(PowerSample {
        power_w: 10.0 * (gas.pressure_pa() / 1.0e6) * (wall.temperature_k() - 280.0),
        source: (gas.pressure_pa(), wall.temperature_k()),
    })
}

#[test]
fn g3_powered_compression_keeps_original_inertia_sources_and_native_continuation() {
    let original_body = specimen(8000.0);
    let mut sourced = MaterialGasPiston::new(parcel(), original_body.clone(), 0.01, -0.03).unwrap();
    // Existing accepted motion and port history must survive the adapter.
    sourced.step(input(), controls(), || false).unwrap();
    let initial = sourced.clone();
    let curve = curve();
    let mut thermal = thermal_body(&curve);
    let mut supplied = 0.0;
    for _ in 0..4 {
        let before = sourced.clone();
        let native = sourced
            .piston()
            .powered_thermal_candidate(
                input(),
                thermal,
                controls(),
                exchange(),
                exchange(),
                conductance,
                power,
                || false,
            )
            .unwrap();
        let candidate = sourced
            .powered_thermal_candidate(
                input(),
                thermal,
                controls(),
                exchange(),
                exchange(),
                conductance,
                power,
                || false,
            )
            .unwrap();
        assert_eq!(sourced, before);
        assert_eq!(candidate.coupled().piston().body(), &original_body);
        assert_eq!(
            candidate.coupled().piston().piston(),
            native.coupled().piston()
        );
        assert_eq!(candidate.coupled().thermal(), native.coupled().thermal());
        assert_eq!(candidate.coupled().step(), native.coupled().step());
        assert_eq!(candidate.power(), native.power());
        assert_eq!(candidate.step(), native.step());
        supplied += candidate.step().evaluated_energy_j;
        let (joint, retained_power, feedback) = candidate.into_parts();
        assert!(retained_power.power_w > 0.0);
        assert!(feedback.total_energy_residual_j.abs() <= exchange().energy_tolerance_j);
        let (next, body, report) = joint.into_parts();
        assert_eq!(next.piston().mass_kg(), original_body.inertia().mass());
        assert_eq!(
            next.piston().last_window_duration_s(),
            Some(input().duration_s)
        );
        assert_eq!(
            next.piston().last_window_transfers(),
            native.coupled().piston().last_window_transfers()
        );
        assert!(report.piston.gas.work_by_gas_j < 0.0);
        thermal.initial = body.state;
        sourced = next;
    }
    let gas_energy = initial.gas().mass_kg()
        * initial
            .gas()
            .model()
            .specific_internal_energy_change_between_densities(
                initial.gas().density_kg_m3(),
                initial.gas().state().temperature_k(),
                sourced.gas().density_kg_m3(),
                sourced.gas().state().temperature_k(),
            )
            .unwrap();
    let kinetic = 0.5
        * original_body.inertia().mass()
        * (sourced.piston().velocity_m_s().powi(2) - initial.piston().velocity_m_s().powi(2));
    let wall_energy = thermal.mass_kg * (thermal.initial.specific_enthalpy_j_kg() - 150000.0);
    let exterior =
        input().exterior_pressure_pa * (sourced.gas().volume_m3() - initial.gas().volume_m3());
    assert!((gas_energy + kinetic + wall_energy + exterior - supplied).abs() < 4.0e-6);
    assert!(sourced.gas().volume_m3() < initial.gas().volume_m3());
    assert!(sourced.gas().state().temperature_k() > initial.gas().state().temperature_k());
    assert_eq!(sourced.body(), initial.body());
}

#[test]
fn g4_power_refusal_and_last_poll_cancellation_preserve_material_and_port_history() {
    let mut sourced = MaterialGasPiston::new(parcel(), specimen(8000.0), 0.01, -0.03).unwrap();
    sourced.step(input(), controls(), || false).unwrap();
    let saved = sourced.clone();
    let curve = curve();
    let thermal = thermal_body(&curve);
    let refused = sourced.powered_thermal_candidate(
        input(),
        thermal,
        controls(),
        exchange(),
        exchange(),
        |_, wall| {
            Ok::<_, &'static str>(ConductanceSample {
                conductance_w_k: 20.0,
                source: wall,
            })
        },
        |_, _| Err::<PowerSample<()>, _>("power source refused"),
        || false,
    );
    assert!(matches!(
        refused,
        Err(ThermalPistonError::Interface("power source refused"))
    ));
    let mut polls = 0;
    sourced
        .powered_thermal_candidate(
            input(),
            thermal,
            controls(),
            exchange(),
            exchange(),
            conductance,
            power,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut current = 0;
    assert!(matches!(
        sourced.powered_thermal_candidate(
            input(),
            thermal,
            controls(),
            exchange(),
            exchange(),
            conductance,
            power,
            || {
                current += 1;
                current == polls
            }
        ),
        Err(ThermalPistonError::Cancelled)
    ));
    assert_eq!(sourced, saved);
    assert_eq!(thermal.initial, thermal_body(&curve).initial);
    assert_eq!(
        sourced.piston().last_window_duration_s(),
        saved.piston().last_window_duration_s()
    );
    assert_eq!(
        sourced.piston().last_window_transfers(),
        saved.piston().last_window_transfers()
    );
}
