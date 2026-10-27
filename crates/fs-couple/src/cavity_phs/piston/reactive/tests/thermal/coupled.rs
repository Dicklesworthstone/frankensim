//! Material heat changes the SAME chemical endpoint, without staged heat release.

use super::super::coupled::{coupled_controls, temperature_reaction};
use super::*;
use crate::cavity_phs::piston::reactive::thermal::coupled::CoupledReactiveThermalError;
use crate::cavity_phs::piston::thermal::boundary::ports::ThermalPort;
use fs_ports::PortedModule;

#[test]
fn wall_material_changes_endpoint_chemistry_and_continues_the_same_history() {
    let run = |conductivity: f64| {
        // Equal-mass, equal-cp isomers isolate thermal/rate feedback. Formation
        // energies are zero; the reference is an implementation fixture.
        let law = temperature_reaction(0.0, 400.0);
        let (card, key, caloric) = wall(conductivity);
        let context = point(850.0);
        let input = PistonInput {
            duration_s: 0.004,
            heat_j: 0.125,
            ..drive()
        };
        let mut body = UniformEnthalpyStepInput {
            initial: caloric
                .curve()
                .state_at_specific_enthalpy(325000.0)
                .unwrap(),
            duration_s: input.duration_s,
            internal_heat_j: 0.375,
            ..thermal_body(&caloric)
        };
        let mut piston = GasPiston::new(gas(&law, [1.0, 0.0, 0.0]), 1.0, 0.01, 0.0)
            .unwrap()
            .with_fixed_support()
            .unwrap();
        let initial = piston.clone();
        let initial_body = body.initial;
        let bounds = HeatExchangeSolve {
            heat_bracket_j: [-10.0, 10.0],
            ..exchange()
        };
        let staged = law
            .parcel_candidate(piston.gas(), input.duration_s, chemical_controls(), || {
                false
            })
            .unwrap()
            .step()
            .extent_mol;
        let mut first_extent = 0.0;
        let mut sum_extent = 0.0;
        let mut total_heat = 0.0;
        for interval in 0..2 {
            let before = piston.clone();
            let before_body = body.initial;
            let candidate = piston
                .coupled_reactive_thermal_layers_candidate(
                    &law,
                    input,
                    body,
                    chemical_controls(),
                    coupled_controls(),
                    bounds,
                    layers(&card, &key, &context),
                    || false,
                )
                .unwrap();
            let step = candidate.step();
            let boundary = &step.boundary;
            let reaction = &boundary.piston.reaction.reaction;
            let gas_end = candidate.piston().gas();
            let gas_temperature = gas_end.state().state().temperature;
            let wall_end = candidate.thermal().state;
            let wall_temperature = wall_end.temperature_k();
            let q = boundary.exchanged_heat_j;
            assert!(q > 0.0);
            assert!(gas_temperature > before.gas().state().state().temperature);
            assert!(wall_temperature > gas_temperature);

            // Independently reduce the original NASA-9 fixture and finite
            // caloric body, without using the coupling's energy report.
            let total_moles = before.gas().mass_kg() / 0.028;
            let gas_capacity = total_moles * 2.5 * R;
            near(
                gas_capacity * (gas_temperature - before.gas().state().state().temperature),
                input.heat_j + q,
                bounds.energy_tolerance_j,
            );
            near(
                body.mass_kg * 500.0 * (wall_temperature - before_body.temperature_k()),
                body.internal_heat_j - q,
                bounds.energy_tolerance_j,
            );

            // The wall's linear conductivity curve and the gas's original
            // mixture transport must describe the accepted endpoints.
            let wall_k = conductivity * (1.0 + (wall_temperature - 200.0) / 800.0);
            near(
                candidate
                    .thermal()
                    .report
                    .property_by_key(&key)
                    .unwrap()
                    .value_si(),
                wall_k,
                1e-12,
            );
            assert_eq!(
                candidate.thermal().report.card_identity(),
                card.content_hash()
            );
            let conductance =
                0.02 / (0.0001 / gas_end.state().state().thermal_conductivity + 0.002 / wall_k);
            near(
                q,
                input.duration_s * conductance * (wall_temperature - gas_temperature),
                bounds.source_tolerance_j,
            );
            near(candidate.thermal().external_heat_j, -q, 1e-14);

            // A <-> B, K=1: backward Euler reduces exactly to this scalar
            // expression. Its source scale is k(T1)=400*T1/600, not k(T0).
            let z = input.duration_s * (400.0 * gas_temperature / 600.0) * R * gas_temperature
                / 100000.0;
            let difference = reaction.species[0].initial_mol - reaction.species[1].initial_mol;
            near(reaction.extent_mol, z * difference / (1.0 + 2.0 * z), 2e-11);
            near(
                reaction.extent_mol,
                input.duration_s
                    * gas_end.volume_m3()
                    * reaction.rate.net_progress_rate().unwrap().value(),
                2e-11,
            );
            near(
                reaction
                    .rate_material
                    .as_ref()
                    .unwrap()
                    .property(ACTIVITY_RATE_SCALE_PROPERTY)
                    .unwrap()
                    .value_si(),
                400.0 * gas_temperature / 600.0,
                1e-9,
            );
            near(
                gas_end.state().state().sound_speed,
                (1.4 * R * gas_temperature / 0.028).sqrt(),
                1e-8,
            );
            assert_eq!(gas_end.mass_kg(), initial.gas().mass_kg());
            assert_eq!(gas_end.volume_m3(), initial.gas().volume_m3());
            assert_eq!(
                candidate.piston().ports.next_tick(),
                before.ports.next_tick() + 1
            );
            assert_eq!(
                candidate.piston().last_window_duration_s(),
                Some(input.duration_s)
            );
            let transfer = candidate.piston().last_window_transfers().unwrap();
            near(transfer[1].energy_into().value(), input.heat_j + q, 1e-12);
            assert!(step.total_energy_residual_j.abs() <= bounds.energy_tolerance_j);
            assert!(step.thermal_energy_residual_j.abs() <= bounds.energy_tolerance_j);
            assert!(boundary.source_residual_j.abs() <= bounds.source_tolerance_j);
            if interval == 0 {
                first_extent = reaction.extent_mol;
                assert!(first_extent > staged);
            }
            sum_extent += reaction.extent_mol;
            total_heat += q;
            let (next, wall, _) = candidate.into_parts();
            piston = next;
            body.initial = wall.state;
        }
        let gas_delta = piston.gas().mass_kg()
            * (piston.gas().state().specific_internal_energy_j_kg()
                - initial.gas().state().specific_internal_energy_j_kg());
        let wall_delta = body.mass_kg
            * (body.initial.specific_enthalpy_j_kg() - initial_body.specific_enthalpy_j_kg());
        near(
            gas_delta + wall_delta,
            2.0 * (input.heat_j + body.internal_heat_j),
            2.0 * bounds.energy_tolerance_j,
        );
        (staged, first_extent, sum_extent, total_heat)
    };
    let resistive = run(0.02);
    let conductive = run(2.0);
    assert_eq!(resistive.0, conductive.0);
    assert!(conductive.1 > resistive.1);
    assert!(conductive.2 > resistive.2);
    assert!(conductive.3 > resistive.3);
}

#[test]
fn late_cancellation_and_wrong_wall_source_leave_both_physical_states_intact() {
    let law = temperature_reaction(0.0, 400.0);
    let (card, key, caloric) = wall(0.2);
    let context = point(500.0);
    let body = thermal_body(&caloric);
    let piston = GasPiston::new(gas(&law, [1.0, 0.0, 0.0]), 1.0, 0.01, 0.0)
        .unwrap()
        .with_fixed_support()
        .unwrap();
    let before = piston.clone();
    let mut polls = 0;
    piston
        .coupled_reactive_thermal_layers_candidate(
            &law,
            drive(),
            body,
            chemical_controls(),
            coupled_controls(),
            exchange(),
            layers(&card, &key, &context),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 4);
    let mut observed = 0;
    let cancelled = piston.coupled_reactive_thermal_layers_candidate(
        &law,
        drive(),
        body,
        chemical_controls(),
        coupled_controls(),
        exchange(),
        layers(&card, &key, &context),
        || {
            observed += 1;
            observed == polls
        },
    );
    assert!(matches!(
        cancelled,
        Err(CoupledReactiveThermalError::Cancelled)
    ));
    assert_eq!(observed, polls);
    let other = wall(2.0).0;
    let wrong_card = piston.coupled_reactive_thermal_layers_candidate(
        &law,
        drive(),
        body,
        chemical_controls(),
        coupled_controls(),
        exchange(),
        layers(&other, &key, &context),
        || false,
    );
    assert!(wrong_card.is_err());
    assert_eq!(piston, before);
    assert_eq!(piston.ports, before.ports);
    assert_eq!(body.initial, thermal_body(&caloric).initial);
}

#[test]
fn opposing_ports_retain_each_transfer_at_one_reacted_endpoint() {
    let law = temperature_reaction(0.0, 400.0);
    let piston = GasPiston::new(gas(&law, [1.0, 0.0, 0.0]), 1.0, 0.01, 0.0)
        .unwrap()
        .with_fixed_support()
        .unwrap();
    let input = PistonInput {
        duration_s: 0.004,
        ..drive()
    };
    let ports = [
        ThermalPort {
            id: 101,
            temperature_k: 700.0,
        },
        ThermalPort {
            id: 7,
            temperature_k: 500.0,
        },
    ];
    let candidate = piston
        .coupled_reactive_thermal_ports_candidate(
            &law,
            input,
            chemical_controls(),
            coupled_controls(),
            exchange(),
            &ports,
            2,
            |port, gas| {
                Ok::<_, Infallible>(ConductanceSample {
                    conductance_w_k: 0.25,
                    source: (port.id, gas.clone()),
                })
            },
            || false,
        )
        .unwrap();
    let samples = candidate.source();
    assert_eq!(
        samples
            .iter()
            .map(|sample| sample.port.id)
            .collect::<Vec<_>>(),
        vec![7, 101]
    );
    near(samples[0].heat_j, -0.1, exchange().source_tolerance_j);
    near(samples[1].heat_j, 0.1, exchange().source_tolerance_j);
    near(
        samples.iter().map(|sample| sample.heat_j).sum(),
        0.0,
        exchange().source_tolerance_j,
    );
    assert_eq!(candidate.step().exchanged_heat_j, 0.0);
    near(
        candidate.step().evaluated_heat_j,
        0.0,
        exchange().source_tolerance_j,
    );
    assert!(candidate.step().piston.reaction.reaction.extent_mol > 0.0);
    for sample in samples {
        assert_eq!(sample.source.0, sample.port.id);
        assert_eq!(&sample.source.1, candidate.piston().gas().state());
        assert_ne!(&sample.source.1, piston.gas().state());
        near(sample.source.1.state().temperature, 600.0, 1e-7);
        near(
            sample.heat_j,
            input.duration_s
                * sample.conductance_w_k
                * (sample.port.temperature_k - sample.source.1.state().temperature),
            1e-14,
        );
    }
    assert_eq!(
        candidate.piston().ports.next_tick(),
        piston.ports.next_tick() + 1
    );
    assert_eq!(
        candidate.piston().last_window_duration_s(),
        Some(input.duration_s)
    );
    let transfer = candidate.piston().last_window_transfers().unwrap();
    assert_eq!(transfer[1].energy_into().value(), 0.0);
    assert!(candidate.step().total_energy_residual_j.abs() <= exchange().energy_tolerance_j);
}
