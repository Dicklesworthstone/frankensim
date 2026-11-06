//! Real-gas work and finite-body heat share the same energy/entropy closure.

use super::*;
use crate::cavity_phs::piston::tests::real_gas::{real_controls, real_parcel};
use fs_material::fluid_transport::FluidTransport;
use fs_material::thermochemical_gas::nonideal::ResolvedThermochemicalRealGasState;

fn drive() -> PistonInput {
    PistonInput {
        duration_s: 0.00002,
        heat_j: 0.0,
        exterior_pressure_pa: 4.5e6,
        external_force_n: -100.0,
    }
}
fn contact() -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-25.0, 25.0],
        source_tolerance_j: 1e-10,
        energy_tolerance_j: 1e-7,
        maximum_iterations: 128,
    }
}
fn resolve(
    gas: &ResolvedThermochemicalRealGasState,
    wall: EquilibriumPhaseState,
) -> Result<ConductanceSample<(f64, f64, f64)>, Infallible> {
    let transport = gas.transport_state().unwrap();
    assert_eq!(transport.density_kg_m3(), gas.density_kg_m3());
    assert_eq!(
        transport.specific_heat_j_kg_k(),
        gas.specific_heat_cp_j_kg_k()
    );
    // Synthetic finite path length, using the original temperature-dependent
    // gas conductivity. The sampled pressure remains the actual PR pressure.
    Ok(ConductanceSample {
        conductance_w_k: 0.02 * transport.thermal_conductivity_w_m_k() / 4e-7,
        source: (gas.temperature_k(), gas.pressure_pa(), wall.temperature_k()),
    })
}

#[test]
fn moving_dense_gas_and_body_close_independent_first_and_second_laws() {
    let curve = sensible();
    for velocity in [-0.4, 0.4] {
        let mut piston = GasPiston::new(real_parcel(300.0, 5e6), 2.0, 0.01, velocity).unwrap();
        let initial = piston.clone();
        let mut thermal = UniformEnthalpyStepInput {
            duration_s: drive().duration_s,
            ..body(&curve, 225000.0)
        };
        let h0 = thermal.initial.specific_enthalpy_j_kg();
        let mut produced = 0.0;
        let mut work = 0.0;
        for _ in 0..12 {
            let candidate = piston
                .thermal_discrete_gradient_candidate(
                    drive(),
                    thermal,
                    real_controls(),
                    contact(),
                    resolve,
                    || false,
                )
                .unwrap();
            let (next, accepted, report) = candidate.into_parts();
            assert!(report.source_residual_j.abs() <= contact().source_tolerance_j);
            assert!(report.total_energy_residual_j.abs() <= contact().energy_tolerance_j);
            assert!(report.contact_entropy_production_j_k >= 0.0);
            assert!(report.entropy_residual_j_k.abs() <= report.entropy_tolerance_j_k);
            assert!(report.total_entropy_change_j_k >= -report.entropy_tolerance_j_k);
            assert_eq!(report.prescribed_source_entropy_j_k, 0.0);
            assert_eq!(
                accepted.report,
                (
                    report.piston.thermodynamic_mean.temperature_k(),
                    report.piston.thermodynamic_mean.pressure_pa(),
                    report.thermal_midpoint.temperature_k(),
                )
            );
            assert_eq!(accepted.external_heat_j, -report.exchanged_heat_j);
            produced += report.contact_entropy_production_j_k;
            work += report.piston.external_force_work_j - report.piston.exterior_pressure_work_j
                + report.piston.support_work_j;
            thermal.initial = accepted.state;
            piston = next;
        }
        let old = initial.gas();
        let new = piston.gas();
        let gas_energy = old.mass_kg()
            * old
                .model()
                .specific_internal_energy_change_between_densities(
                    old.density_kg_m3(),
                    old.state().temperature_k(),
                    new.density_kg_m3(),
                    new.state().temperature_k(),
                )
                .unwrap();
        let gas_entropy = old.mass_kg()
            * old
                .model()
                .specific_entropy_change_between_densities(
                    old.density_kg_m3(),
                    old.state().temperature_k(),
                    new.density_kg_m3(),
                    new.state().temperature_k(),
                )
                .unwrap();
        let body_energy = thermal.mass_kg * (thermal.initial.specific_enthalpy_j_kg() - h0);
        let body_entropy = thermal.mass_kg
            * curve
                .specific_entropy_change_j_kg_k(h0, thermal.initial.specific_enthalpy_j_kg())
                .unwrap();
        let kinetic = 0.5
            * piston.mass_kg()
            * (piston.velocity_m_s().powi(2) - initial.velocity_m_s().powi(2));
        assert!(
            (gas_energy + body_energy + kinetic - work).abs() < 12.0 * contact().energy_tolerance_j
        );
        assert!((gas_entropy + body_entropy - produced).abs() < 1e-8);
        assert!(produced > 0.0);
        assert!(thermal.initial.temperature_k() < 650.0);
        assert_eq!(new.mass_kg(), old.mass_kg());
        assert!(piston.audit_last_window().is_ok());
    }
}

#[test]
fn dense_fixed_gas_receives_latent_heat_without_inventing_a_mean_phase_state() {
    let curve = latent(350.0);
    let piston = GasPiston::new(real_parcel(300.0, 5e6), 2.0, 0.01, 0.0)
        .unwrap()
        .with_fixed_support()
        .unwrap();
    let thermal = UniformEnthalpyStepInput {
        duration_s: drive().duration_s,
        ..body(&curve, 130000.0)
    };
    let controls = PistonSolve {
        volume_bracket_m3: [piston.gas().volume_m3(); 2],
        ..real_controls()
    };
    let candidate = piston
        .thermal_discrete_gradient_candidate(drive(), thermal, controls, contact(), resolve, || {
            false
        })
        .unwrap();
    let report = candidate.step();
    assert_eq!(report.thermal_mean_temperature_k, 350.0);
    assert_eq!(report.thermal_midpoint.temperature_k(), 350.0);
    assert_eq!(candidate.thermal().state.temperature_k(), 350.0);
    assert!(
        candidate.thermal().state.liquid_mass_fraction() < thermal.initial.liquid_mass_fraction()
    );
    assert_eq!(report.piston.gas.work_by_gas_j, 0.0);
    assert!(report.exchanged_heat_j > 0.0);
    assert!(report.contact_entropy_production_j_k > 0.0);
    assert!(report.entropy_residual_j_k.abs() <= report.entropy_tolerance_j_k);
    assert!(report.total_energy_residual_j.abs() <= contact().energy_tolerance_j);
    let mut polls = 0;
    piston
        .thermal_discrete_gradient_candidate(drive(), thermal, controls, contact(), resolve, || {
            polls += 1;
            false
        })
        .unwrap();
    let before = piston.clone();
    let mut current = 0;
    assert!(matches!(
        piston.thermal_discrete_gradient_candidate(
            drive(),
            thermal,
            controls,
            contact(),
            resolve,
            || {
                current += 1;
                current == polls
            },
        ),
        Err(ThermalPistonError::Cancelled)
    ));
    assert_eq!(piston, before);
    assert_eq!(
        piston.last_window_transfers(),
        before.last_window_transfers()
    );
}
