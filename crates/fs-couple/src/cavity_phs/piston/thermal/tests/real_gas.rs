//! The same reciprocal body, gas-energy and momentum solve with dense PR gas.
//! Manufactured material laws; these are conservation and source-view tests.

use super::{body, curve, exchange};
use crate::cavity_phs::piston::material::{MaterialGasPiston, tests::body as inertial_body};
use crate::cavity_phs::piston::tests::real_gas::{real_controls, real_parcel};
use crate::cavity_phs::piston::thermal::{ConductanceSample, ThermalPistonError};
use crate::cavity_phs::piston::{GasPiston, PistonInput};
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::fluid_transport::FluidTransport;
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStepInput};
use fs_material::thermochemical_gas::nonideal::{
    ResolvedThermochemicalRealGasState, energy::parcel::IsochoricRealGasParcel,
};
use std::convert::Infallible;

fn piston() -> GasPiston<IsochoricRealGasParcel> {
    GasPiston::new(real_parcel(300.0, 5e6), 2.0, 0.01, 0.4).unwrap()
}

fn input() -> PistonInput {
    PistonInput {
        duration_s: 0.00005,
        heat_j: 2.0,
        exterior_pressure_pa: 4.5e6,
        external_force_n: -100.0,
    }
}

fn link(
    conductance: f64,
) -> impl FnMut(
    &ResolvedThermochemicalRealGasState,
    EquilibriumPhaseState,
) -> Result<ConductanceSample<ResolvedThermochemicalRealGasState>, Infallible> {
    move |gas, _| {
        Ok(ConductanceSample {
            conductance_w_k: conductance,
            source: gas.clone(),
        })
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_dense_moving_piston_and_finite_body_close_signed_reciprocal_energy() {
    let curve = curve(500.0);
    for temperature in [500.0, 250.0] {
        let old = piston();
        let saved = old.clone();
        let thermal = body(&curve, 500.0 * (temperature - 200.0));
        let candidate = old
            .thermal_candidate(
                input(),
                thermal,
                real_controls(),
                exchange(),
                link(50.0),
                || false,
            )
            .unwrap();
        let next = candidate.piston();
        let report = candidate.step();
        let midpoint = &report.piston.midpoint;
        let heat = report.exchanged_heat_j;
        assert_eq!(heat > 0.0, temperature > 300.0);
        assert!(next.gas().volume_m3() > old.gas().volume_m3());
        assert!(next.velocity_m_s() > old.velocity_m_s());
        assert_eq!(next.gas().mass_kg(), old.gas().mass_kg());
        assert_eq!(
            next.gas().density_kg_m3(),
            next.gas().mass_kg() / next.gas().volume_m3()
        );
        assert_eq!(candidate.thermal().report, *midpoint);
        assert_eq!(candidate.thermal().external_heat_j, -heat);
        assert_ne!(midpoint.density_kg_m3(), midpoint.source().state().density);
        assert_eq!(
            midpoint.transport_state().unwrap().density_kg_m3(),
            midpoint.density_kg_m3()
        );
        assert_eq!(
            midpoint.acoustic_state().unwrap().density_kg_m3(),
            midpoint.density_kg_m3()
        );
        assert_eq!(next.gas().point().axes().get("specimen-lot"), Some(&42.0));
        assert_eq!(next.gas().model(), old.gas().model());
        near(
            heat,
            input().duration_s
                * 50.0
                * (report.thermal_midpoint.temperature_k() - midpoint.temperature_k()),
            exchange().source_tolerance_j,
        );

        // Independent endpoint accounting uses the stored source energies,
        // wall velocities and volumes, not the coupling's reported residual.
        let gas_energy = old.gas().mass_kg()
            * (next.gas().state().specific_internal_energy_j_kg()
                - old.gas().state().specific_internal_energy_j_kg());
        let body_energy = thermal.mass_kg
            * (candidate.thermal().state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        let kinetic = 0.5
            * old.mass_kg()
            * (next.velocity_m_s() - old.velocity_m_s())
            * (next.velocity_m_s() + old.velocity_m_s());
        let dv = next.gas().volume_m3() - old.gas().volume_m3();
        let exterior = input().exterior_pressure_pa * dv;
        let applied = input().external_force_n * dv / old.area_m2();
        near(body_energy + heat, 0.0, exchange().energy_tolerance_j);
        near(
            gas_energy + body_energy + kinetic + exterior - applied,
            input().heat_j,
            exchange().energy_tolerance_j,
        );
        assert!(report.total_energy_residual_j.abs() <= exchange().energy_tolerance_j);
        assert!(
            report.piston.momentum_residual_n_s.abs() <= real_controls().momentum_tolerance_n_s
        );
        assert_eq!(old, saved);
        assert_eq!(
            thermal.initial,
            curve
                .state_at_specific_enthalpy(500.0 * (temperature - 200.0))
                .unwrap()
        );
    }
}

#[test]
fn g4_real_thermal_source_refusal_and_final_cancellation_return_no_joint_state() {
    let old = piston();
    let saved = old.clone();
    let curve = curve(500.0);
    let thermal = body(&curve, 150000.0);
    let initial_body = thermal.initial;
    let failure =
        |_: &ResolvedThermochemicalRealGasState, _| -> Result<ConductanceSample<()>, &str> {
            Err("missing real-gas interface source")
        };
    assert!(matches!(
        old.thermal_candidate(
            input(),
            thermal,
            real_controls(),
            exchange(),
            failure,
            || false
        ),
        Err(ThermalPistonError::Interface(
            "missing real-gas interface source"
        ))
    ));
    let mut total = 0;
    old.thermal_candidate(
        input(),
        thermal,
        real_controls(),
        exchange(),
        link(50.0),
        || {
            total += 1;
            false
        },
    )
    .unwrap();
    let mut polls = 0;
    assert!(matches!(
        old.thermal_candidate(
            input(),
            thermal,
            real_controls(),
            exchange(),
            link(50.0),
            || {
                polls += 1;
                polls == total
            }
        ),
        Err(ThermalPistonError::Cancelled)
    ));
    assert_eq!(old, saved);
    assert_eq!(thermal.initial, initial_body);
}

#[test]
fn g3_source_density_and_inertia_survive_joint_thermal_acceptance_and_continuation() {
    let curve = curve(500.0);
    let drive = PistonInput {
        duration_s: 5e-6,
        ..input()
    };
    let thermal = UniformEnthalpyStepInput {
        duration_s: drive.duration_s,
        ..body(&curve, 150000.0)
    };
    let initial_gas = real_parcel(300.0, 5e6);
    let mut velocities = [0.0; 2];
    for (index, density) in [1000.0, 8000.0].into_iter().enumerate() {
        let specimen = inertial_body(density);
        let old = MaterialGasPiston::new(initial_gas.clone(), specimen.clone(), 0.01, 0.0).unwrap();
        let saved = old.clone();
        let mut total_polls = 0;
        let candidate = old
            .thermal_candidate(
                drive,
                thermal,
                real_controls(),
                exchange(),
                link(50.0),
                || {
                    total_polls += 1;
                    false
                },
            )
            .unwrap();
        assert_eq!(candidate.piston().body(), &specimen);
        assert_eq!(
            candidate.piston().piston().mass_kg(),
            specimen.inertia().mass()
        );
        assert_eq!(candidate.piston().gas().mass_kg(), initial_gas.mass_kg());
        assert_eq!(candidate.thermal().report, candidate.step().piston.midpoint);
        assert_eq!(
            fs_ports::PortedModule::last_window_duration_s(candidate.piston().piston()),
            Some(drive.duration_s)
        );
        assert!(candidate.step().exchanged_heat_j > 0.0);
        assert!(candidate.step().total_energy_residual_j.abs() <= exchange().energy_tolerance_j);
        let mut polls = 0;
        assert!(matches!(
            old.thermal_candidate(
                drive,
                thermal,
                real_controls(),
                exchange(),
                link(50.0),
                || {
                    polls += 1;
                    polls == total_polls
                }
            ),
            Err(ThermalPistonError::Cancelled)
        ));
        assert_eq!(old, saved);

        // Continue from the complete returned wrapper and the jointly accepted
        // thermal state. Reconstructing a numeric piston would lose its body.
        let (next, accepted_body, _) = candidate.into_parts();
        let continued_body = UniformEnthalpyStepInput {
            initial: accepted_body.state,
            ..thermal
        };
        let continued = next
            .thermal_candidate(
                drive,
                continued_body,
                real_controls(),
                exchange(),
                link(50.0),
                || false,
            )
            .unwrap();
        let (next, accepted_body, report) = continued.into_parts();
        assert_eq!(next.body(), &specimen);
        assert_eq!(next.piston().mass_kg(), specimen.inertia().mass());
        assert_eq!(next.gas().mass_kg(), initial_gas.mass_kg());
        assert_eq!(next.gas().model(), initial_gas.model());
        assert_eq!(
            fs_ports::PortedModule::last_window_duration_s(next.piston()),
            Some(drive.duration_s)
        );
        assert_eq!(
            next.gas().density_kg_m3(),
            next.gas().mass_kg() / next.gas().volume_m3()
        );
        assert!(report.total_energy_residual_j.abs() <= exchange().energy_tolerance_j);
        let gas_energy = initial_gas.mass_kg()
            * (next.gas().state().specific_internal_energy_j_kg()
                - initial_gas.state().specific_internal_energy_j_kg());
        let body_energy = thermal.mass_kg
            * (accepted_body.state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        let kinetic = 0.5 * specimen.inertia().mass() * next.piston().velocity_m_s().powi(2);
        let dv = next.gas().volume_m3() - initial_gas.volume_m3();
        let outside =
            drive.exterior_pressure_pa * dv - drive.external_force_n * dv / next.piston().area_m2();
        near(
            gas_energy + body_energy + kinetic + outside,
            2.0 * drive.heat_j,
            2.0 * exchange().energy_tolerance_j,
        );
        velocities[index] = next.piston().velocity_m_s();
    }
    assert!(velocities[0] > velocities[1]);
    assert!(velocities[1] > 0.0);
}
