//! Three-owner conservative coupling, using the original synthetic source laws.

use super::{close, controls, drive, exchange, hotter, nozzle, parcel};
use crate::cavity_phs::piston::exchange::NozzleDirection;
use crate::cavity_phs::piston::material::{MaterialGasPiston, tests::body};
use crate::cavity_phs::piston::thermal::{
    ConductanceSample, HeatExchangeSolve, ThermalPistonError,
};
use crate::cavity_phs::piston::{GasPiston, PistonInput};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, SolidLiquidPhase,
    UniformEnthalpyStepInput,
};
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzle;
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ResolvedThermochemicalRealGasState,
    ThermochemicalPengRobinsonModel, energy::parcel::IsochoricRealGasParcel,
};
use fs_ports::{ModuleWindowAudit, PortedModule};
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseV1;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Nasa9RegionV1, Nasa9StandardStateModelV1};
use std::convert::Infallible;

fn curve() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        fs_blake3::hash_bytes(b"synthetic real-gas nozzle lining"),
        SolidLiquidPhase::Solid,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 200.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 4000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 300_000.0,
                temperature_k: 800.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 4000.0,
            },
        ],
    )
    .unwrap()
}

fn thermal_body(
    curve: &EquilibriumEnthalpyPhaseCurve,
    temperature_k: f64,
    input: PistonInput,
) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve,
        initial: curve
            .state_at_specific_enthalpy(500.0 * (temperature_k - 200.0))
            .unwrap(),
        mass_kg: 0.001,
        volume_m3: 0.001 / 4000.0,
        surface_area_m2: 0.02,
        internal_heat_j: 0.02,
        duration_s: input.duration_s,
    }
}

fn heat_solve() -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-1.0, 1.0],
        source_tolerance_j: 1e-9,
        energy_tolerance_j: 2e-6,
        maximum_iterations: 128,
    }
}

fn interface(
    gas: &ResolvedThermochemicalRealGasState,
    _: EquilibriumPhaseState,
) -> Result<ConductanceSample<ResolvedThermochemicalRealGasState>, Infallible> {
    Ok(ConductanceSample {
        conductance_w_k: 50.0,
        source: gas.clone(),
    })
}

#[test]
fn real_flow_in_both_directions_and_signed_wall_heat_close_all_three_stores() {
    for into_piston in [true, false] {
        let (chamber, reservoir, wall_temperature) = if into_piston {
            (parcel(), hotter(parcel()), 500.0)
        } else {
            (hotter(parcel()), parcel(), 280.0)
        };
        let piston = GasPiston::new(chamber.clone(), 2.0, 0.01, 0.03).unwrap();
        let before = piston.clone();
        let old_reservoir = reservoir.clone();
        let input = drive(0.03);
        let curve = curve();
        let thermal = thermal_body(&curve, wall_temperature, input);

        // Compare against ONE original transfer, independently of heat trials.
        let (native, native_flow) = if into_piston {
            reservoir.through_nozzle_to(&chamber, nozzle(), input.duration_s, exchange(), || false)
        } else {
            chamber.through_nozzle_to(&reservoir, nozzle(), input.duration_s, exchange(), || false)
        }
        .unwrap();
        let (donor, receiver, native_transfer) = native.into_parts();
        let expected_reservoir = if into_piston { donor } else { receiver };
        let accepted = piston
            .real_gas_nozzle_thermal_candidate(
                &reservoir,
                nozzle(),
                input,
                thermal,
                controls(),
                exchange(),
                heat_solve(),
                interface,
                || false,
            )
            .unwrap();
        let step = accepted.step();
        let next = accepted.piston();
        assert_eq!(
            step.direction,
            Some(if into_piston {
                NozzleDirection::IntoPiston
            } else {
                NozzleDirection::OutOfPiston
            })
        );
        assert_eq!(step.flow, native_flow);
        assert_eq!(step.transfer, native_transfer);
        assert_eq!(accepted.reservoir(), &expected_reservoir);
        assert!(step.transfer.transferred_mass_kg > 0.0);
        assert!(step.thermal.exchange_trials > 1);
        assert_eq!(step.thermal.exchanged_heat_j > 0.0, into_piston);
        assert_eq!(
            accepted.thermal().external_heat_j,
            -step.thermal.exchanged_heat_j
        );
        assert_eq!(accepted.thermal().report, step.thermal.piston.midpoint);
        close(
            step.thermal.exchanged_heat_j,
            input.duration_s
                * 50.0
                * (step.thermal.thermal_midpoint.temperature_k()
                    - step.thermal.piston.midpoint.temperature_k()),
            heat_solve().source_tolerance_j,
        );
        close(
            next.gas().state().density_kg_m3(),
            next.gas().mass_kg() / next.gas().volume_m3(),
            1e-12,
        );
        assert_eq!(next.gas().model(), chamber.model());
        assert_eq!(next.gas().point().axes().get("specimen-lot"), Some(&42.0));
        assert!(next.gas().state().density_kg_m3() != next.gas().state().source().state().density);

        // Independent endpoint balance is resolved for this zero-offset fixture.
        // Formation-offset robustness is checked separately through increments.
        let gas_energy = |gas: &IsochoricRealGasParcel| {
            gas.mass_kg() * gas.state().specific_internal_energy_j_kg()
        };
        let stored = gas_energy(next.gas()) - gas_energy(&chamber)
            + gas_energy(accepted.reservoir())
            - gas_energy(&reservoir);
        let body_stored = thermal.mass_kg
            * (accepted.thermal().state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        let kinetic =
            0.5 * piston.mass_kg() * (next.velocity_m_s().powi(2) - piston.velocity_m_s().powi(2));
        let dv = next.gas().volume_m3() - chamber.volume_m3();
        close(
            stored + body_stored + kinetic + input.exterior_pressure_pa * dv
                - input.external_force_n * dv / piston.area_m2()
                - input.heat_j
                - thermal.internal_heat_j,
            0.0,
            heat_solve().energy_tolerance_j,
        );
        let mass = (next.gas().mass_kg() - chamber.mass_kg())
            + (accepted.reservoir().mass_kg() - reservoir.mass_kg());
        assert_eq!(mass, step.total_mass_residual_kg);
        assert!(mass.abs() <= exchange().transfer.mass_tolerance_kg);
        assert!(step.total_energy_residual_j.abs() <= heat_solve().energy_tolerance_j);
        assert_eq!(piston, before);
        assert_eq!(reservoir, old_reservoir);
    }
}

#[test]
fn shut_nozzle_and_insulated_wall_reduce_to_the_existing_native_owners() {
    let piston = GasPiston::new(parcel(), 2.0, 0.01, 0.0).unwrap();
    let reservoir = hotter(parcel());
    let curve = curve();
    let input = drive(0.03);
    let thermal = thermal_body(&curve, 500.0, input);
    let expected = piston
        .thermal_candidate(input, thermal, controls(), heat_solve(), interface, || {
            false
        })
        .unwrap();
    let shut = piston
        .real_gas_nozzle_thermal_candidate(
            &reservoir,
            GasNozzle {
                area_m2: 0.0,
                ..nozzle()
            },
            input,
            thermal,
            controls(),
            exchange(),
            heat_solve(),
            interface,
            || false,
        )
        .unwrap();
    assert_eq!(shut.piston(), expected.piston());
    assert_eq!(shut.reservoir(), &reservoir);
    assert_eq!(shut.thermal(), expected.thermal());
    assert_eq!(&shut.step().thermal, expected.step());
    assert_eq!(shut.step().direction, None);
    assert_eq!(shut.step().transfer.transferred_mass_kg, 0.0);

    let thermal = UniformEnthalpyStepInput {
        internal_heat_j: 0.0,
        ..thermal
    };
    let mut expected_piston = piston.clone();
    let mut expected_reservoir = reservoir.clone();
    let flow_only = expected_piston
        .step_with_real_gas_nozzle(
            &mut expected_reservoir,
            nozzle(),
            input,
            exchange(),
            controls(),
            heat_solve().energy_tolerance_j,
            || false,
        )
        .unwrap();
    let insulated = piston
        .real_gas_nozzle_thermal_candidate(
            &reservoir,
            nozzle(),
            input,
            thermal,
            controls(),
            exchange(),
            heat_solve(),
            |gas, _| {
                Ok::<_, Infallible>(ConductanceSample {
                    conductance_w_k: 0.0,
                    source: gas.clone(),
                })
            },
            || false,
        )
        .unwrap();
    assert_eq!(insulated.piston(), &expected_piston);
    assert_eq!(insulated.reservoir(), &expected_reservoir);
    assert_eq!(insulated.step().thermal.piston, flow_only.piston);
    assert_eq!(insulated.step().transfer, flow_only.transfer);
    assert_eq!(insulated.thermal().state, thermal.initial);
    assert_eq!(insulated.step().thermal.exchanged_heat_j, 0.0);
}

#[test]
fn source_bound_wall_density_changes_the_continuing_flow_heat_motion_trajectory() {
    let curve = curve();
    let input = PistonInput {
        duration_s: 5e-6,
        ..drive(0.3)
    };
    let run = |density| {
        let initial = parcel();
        let specimen = body(density);
        let mut piston =
            MaterialGasPiston::new(initial.clone(), specimen.clone(), 0.01, 0.0).unwrap();
        let mut reservoir = hotter(initial.clone());
        let total_mass = piston.gas().mass_kg() + reservoir.mass_kg();
        let mut thermal = thermal_body(&curve, 500.0, input);
        for _ in 0..3 {
            let candidate = piston
                .real_gas_nozzle_thermal_candidate(
                    &reservoir,
                    nozzle(),
                    input,
                    thermal,
                    controls(),
                    exchange(),
                    heat_solve(),
                    interface,
                    || false,
                )
                .unwrap();
            let (next, next_reservoir, next_thermal, report) = candidate.into_parts();
            assert_eq!(next.body(), &specimen);
            assert_eq!(next.piston().mass_kg(), specimen.inertia().mass());
            assert_eq!(next.gas().model(), initial.model());
            assert!(next.gas().mass_kg() > piston.gas().mass_kg());
            assert!(report.thermal.exchanged_heat_j > 0.0);
            assert!(report.total_energy_residual_j.abs() <= heat_solve().energy_tolerance_j);
            close(
                next.gas().mass_kg() + next_reservoir.mass_kg(),
                total_mass,
                exchange().transfer.mass_tolerance_kg,
            );
            assert_eq!(
                next.piston().last_window_duration_s(),
                Some(input.duration_s)
            );
            assert!(matches!(
                next.piston().audit_last_window().unwrap(),
                ModuleWindowAudit::Unavailable { .. }
            ));
            piston = next;
            reservoir = next_reservoir;
            thermal.initial = next_thermal.state;
        }
        piston
    };
    let light = run(1000.0);
    let dense = run(8000.0);
    close(light.piston().mass_kg(), 0.1, 1e-14);
    close(dense.piston().mass_kg(), 0.8, 1e-14);
    assert!(light.piston().velocity_m_s() > dense.piston().velocity_m_s());
    assert!(dense.piston().velocity_m_s() > 0.0);
    assert!(light.gas().volume_m3() > dense.gas().volume_m3());
}

#[test]
fn source_refusal_and_last_wrapper_poll_publish_no_partial_endpoint_or_port_window() {
    let mut piston = MaterialGasPiston::new(parcel(), body(8000.0), 0.01, 0.0).unwrap();
    piston.step(drive(0.1), controls(), || false).unwrap();
    let before = piston.clone();
    let reservoir = hotter(parcel());
    let old_reservoir = reservoir.clone();
    let curve = curve();
    let input = drive(0.03);
    let thermal = thermal_body(&curve, 500.0, input);
    let mut observed_midpoint = false;
    let refused = piston.real_gas_nozzle_thermal_candidate(
        &reservoir,
        nozzle(),
        input,
        thermal,
        controls(),
        exchange(),
        heat_solve(),
        |gas, _| {
            assert_eq!(gas.source().components().len(), 2);
            observed_midpoint = true;
            Err::<ConductanceSample<ResolvedThermochemicalRealGasState>, _>(
                "wall source unavailable",
            )
        },
        || false,
    );
    assert!(observed_midpoint);
    assert_eq!(
        refused,
        Err(ThermalPistonError::Interface("wall source unavailable"))
    );
    assert_eq!(piston, before);
    assert_eq!(reservoir, old_reservoir);

    let mut polls = 0;
    piston
        .real_gas_nozzle_thermal_candidate(
            &reservoir,
            nozzle(),
            input,
            thermal,
            controls(),
            exchange(),
            heat_solve(),
            interface,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 5);
    let mut observed = 0;
    let cancelled = piston.real_gas_nozzle_thermal_candidate(
        &reservoir,
        nozzle(),
        input,
        thermal,
        controls(),
        exchange(),
        heat_solve(),
        interface,
        || {
            observed += 1;
            observed == polls
        },
    );
    assert_eq!(cancelled, Err(ThermalPistonError::Cancelled));
    assert_eq!(piston, before);
    assert_eq!(reservoir, old_reservoir);
    assert_eq!(
        piston.piston().last_window_transfers(),
        before.piston().last_window_transfers()
    );
}

// Change only the NASA enthalpy datum in the existing manufactured cards.
// Density, caloric derivatives, entropy, transport and PR bindings stay original.
fn with_reference(reference_k: f64) -> IsochoricRealGasParcel {
    let original = parcel();
    let source = original.model().source();
    let components = source
        .thermodynamics()
        .components()
        .iter()
        .map(|component| {
            let regions = component
                .regions()
                .iter()
                .map(|region| {
                    let mut card = region.card().clone();
                    card.parameters.get_mut("a7").unwrap().value = reference_k;
                    Nasa9RegionV1::from_card(card).unwrap()
                })
                .collect();
            Nasa9StandardStateModelV1::new(
                component.species().clone(),
                component.molar_mass(),
                component.convention().clone(),
                regions,
            )
            .unwrap()
        })
        .collect();
    let thermodynamics = FrozenIdealGasMixtureModelV1::new(
        components,
        source.thermodynamics().declared_composition().clone(),
    )
    .unwrap();
    let bindings = source
        .thermodynamics()
        .components()
        .iter()
        .zip(original.model().mechanics().components())
        .map(|(component, pr)| PengRobinsonSpeciesBinding {
            species: component.species().clone(),
            critical_temperature: pr.critical_temperature,
            critical_pressure: pr.critical_pressure,
            acentric_factor: pr.acentric_factor,
        })
        .collect();
    let model = ThermochemicalPengRobinsonModel::new(
        ThermochemicalGasModel::new(
            thermodynamics,
            source.transport_materials().to_vec(),
            source.water_species(),
        )
        .unwrap(),
        bindings,
        original.model().mechanics().binary_interactions().to_vec(),
    )
    .unwrap();
    IsochoricRealGasParcel::new(
        model,
        original.point().clone(),
        original.volume_m3(),
        PengRobinsonPhaseV1::Vapor,
    )
    .unwrap()
}

#[test]
fn enormous_formation_reference_does_not_change_flow_wall_heat_or_motion() {
    let curve = curve();
    let input = drive(0.03);
    let thermal = thermal_body(&curve, 500.0, input);
    let run = |reference| {
        let chamber = with_reference(reference);
        let reservoir = hotter(chamber.clone());
        GasPiston::new(chamber, 2.0, 0.01, 0.03)
            .unwrap()
            .real_gas_nozzle_thermal_candidate(
                &reservoir,
                nozzle(),
                input,
                thermal,
                controls(),
                exchange(),
                heat_solve(),
                interface,
                || false,
            )
            .unwrap()
    };
    let reference = run(0.0);
    for datum in [-1e24, 1e24] {
        let shifted = run(datum);
        close(
            shifted.step().flow.mass_flow_kg_s,
            reference.step().flow.mass_flow_kg_s,
            1e-12,
        );
        close(
            shifted.step().thermal.exchanged_heat_j,
            reference.step().thermal.exchanged_heat_j,
            1e-9,
        );
        close(
            shifted.piston().velocity_m_s(),
            reference.piston().velocity_m_s(),
            1e-12,
        );
        close(
            shifted.thermal().state.temperature_k(),
            reference.thermal().state.temperature_k(),
            1e-9,
        );
        for (actual, expected) in [
            (shifted.piston().gas(), reference.piston().gas()),
            (shifted.reservoir(), reference.reservoir()),
        ] {
            close(actual.mass_kg(), expected.mass_kg(), 1e-15);
            close(
                actual.state().temperature_k(),
                expected.state().temperature_k(),
                1e-9,
            );
            close(
                actual.state().pressure_pa(),
                expected.state().pressure_pa(),
                1e-4,
            );
            close(
                actual.state().density_kg_m3(),
                expected.state().density_kg_m3(),
                1e-10,
            );
        }
        assert!(shifted.step().transfer.transferred_enthalpy_j.abs() > 1e12);
        assert!(shifted.step().total_energy_residual_j.abs() <= heat_solve().energy_tolerance_j);
        assert!(
            shifted.step().total_mass_residual_kg.abs() <= exchange().transfer.mass_tolerance_kg
        );
    }
}
