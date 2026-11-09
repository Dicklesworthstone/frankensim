//! Original nozzle transport, regional phase/heat, and native pressure act once.
use super::*;
use fs_couple::liquid_gas::evaporation::{
    multiple::regions::pressure::{
        RegionalGasPressure, advance_liquid_regions_with_gas_and_pressure,
    },
    transport::regions::{
        advance_liquid_regions_with_gas_transport_and_heat_and_pressure,
        advance_liquid_regions_with_gas_transport_and_pressure,
    },
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;

const REFERENCE: f64 = 100000.0;
const MAXIMUM: f64 = 0.05;

fn membership() -> BTreeMap<[usize; 3], u64> {
    let mut result = BTreeMap::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if x > 3 || y == 0 || y == 7 || z == 0 || z == 7 {
                    result.insert([x, y, z], if z >= 4 { 0 } else { 1 });
                }
            }
        }
    }
    result
}

fn pressure_parcels(source: &ThermochemicalGasModel) -> BTreeMap<u64, ClosedGasParcel> {
    BTreeMap::from([
        (0, node(source, 300.0, 100000.0, 0.01, 3e-7)),
        (1, node(source, 300.0, 100001.0, 0.06, 3e-7)),
        (9, node(source, 300.0, 90000.0, 0.05, 1e-8)),
    ])
}

fn pressure_from_inventory(gas: &ClosedGasParcel) -> f64 {
    let amount = gas.species_mass_kg(&carrier()).unwrap() / 0.028
        + gas.species_mass_kg(&water()).unwrap() / M;
    amount * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * gas.state().state().temperature
        / gas.volume_m3()
}

#[test]
fn transport_then_regional_phase_and_heat_load_actual_final_pressures_in_one_native_step() {
    for with_heat in [false, true] {
        with_liquid_interfaces_pressure(0.0, &CELLS, Some(REFERENCE), |liquid, gate| {
            let source = source(-13.4e6);
            let parcels = pressure_parcels(&source);
            let preserved = parcels.clone();
            let old_bits = liquid.flow().lattice().grid().state_bits();
            let old_h = liquid.total_enthalpy_j().unwrap();
            let binding = binding(&source, 0.0);
            let rows = regional_faces(&binding);
            let warm_rows = heated(&rows, 50.0);
            let bounds = endpoints(&parcels);
            let links = [nozzle(1e-10)];
            let paths = GasTransportPaths { nozzles: &links, ..Default::default() };
            let moved = plumbing::step_candidate_with_paths(
                &parcels, paths, DT, &bounds, transport_controls(), || false,
            )
            .unwrap();
            assert!(moved.step().nozzles[0].transferred_mass_kg > 0.0);
            let gas_regions = membership();
            let pressure = RegionalGasPressure {
                gas_regions: &gas_regions,
                maximum_density_excursion: MAXIMUM,
            };
            let thermal = GasLiquidChamberHeatSolve { maximum_iterations: 128 };
            let next = if with_heat {
                advance_liquid_regions_with_gas_transport_and_heat_and_pressure(
                    liquid, &parcels, &source, paths, &bounds, &warm_rows,
                    transport_controls(), joint_controls(), thermal, pressure, gate,
                )
            } else {
                advance_liquid_regions_with_gas_transport_and_pressure(
                    liquid, &parcels, &source, paths, &bounds, &rows,
                    transport_controls(), joint_controls(), pressure, gate,
                )
            }
            .unwrap();
            let unloaded = if with_heat {
                advance_liquid_regions_with_gas_transport_and_heat(
                    liquid, &parcels, &source, paths, &bounds, &warm_rows,
                    transport_controls(), joint_controls(), thermal, gate,
                )
            } else {
                advance_liquid_regions_with_gas_transport(
                    liquid, &parcels, &source, paths, &bounds, &rows,
                    transport_controls(), joint_controls(), gate,
                )
            }
            .unwrap();
            let report = next.report();
            let receipt = report.coupled.pressure.as_ref().unwrap();
            assert_eq!(&report.transport, moved.step());
            assert_eq!(next.parcels(), unloaded.parcels());
            assert_eq!(report.coupled.mass_exchange, unloaded.report().coupled.mass_exchange);
            assert!(unloaded.report().coupled.pressure.is_none());
            for (index, (cell, _, _)) in CELLS.iter().enumerate() {
                let id = index as u64;
                let load = receipt.regions[&id];
                // The pressure stage starts AFTER actual nozzle transport,
                // then applies the EOS state after regional heat and phase.
                assert_eq!(load.initial_pressure_pa, moved.parcels()[&id].state().state().pressure);
                assert_ne!(load.initial_pressure_pa, parcels[&id].state().state().pressure);
                assert_eq!(load.applied_pressure_pa, next.parcels()[&id].state().state().pressure);
                assert_eq!(load.volume_m3, parcels[&id].volume_m3());
                assert_eq!(load.maximum_density_excursion, MAXIMUM);
                assert!(receipt.sampled_links[&id] > 0);
                assert_eq!(
                    report.coupled.regions[&id].exchange,
                    unloaded.report().coupled.regions[&id].exchange
                );
                assert_eq!(
                    report.coupled.regions[&id].heat,
                    unloaded.report().coupled.regions[&id].heat
                );
                assert_eq!(report.coupled.regions[&id].heat.is_some(), with_heat);
                // Independent inventory EOS and original five-link flat-front
                // momentum oracle; all those links lie in the same region.
                let expected_pressure = pressure_from_inventory(&next.parcels()[&id]);
                near(load.applied_pressure_pa, expected_pressure, 2e-10);
                let d = (expected_pressure - REFERENCE) / 1000.0;
                let velocity = next.liquid().flow()
                    .sample_at_si(cell[0], cell[1], cell[2])
                    .unwrap().unwrap().velocity_m_s[0];
                near(velocity, -d / (1.0 + d), 2e-11);
            }
            assert_ne!(
                next.liquid().flow().lattice().grid().state_bits(),
                unloaded.liquid().flow().lattice().grid().state_bits()
            );
            assert_eq!(next.parcels()[&9], parcels[&9]);
            assert_eq!(next.liquid().accepted_steps(), 1);
            assert_eq!(next.liquid().physical_time_s(), DT);
            assert_eq!(report.duration_s, DT);
            assert_eq!(report.coupled.duration_s, DT);
            assert_eq!(report.coupled.flow.external_heat_j, 0.0);
            let dm = next.liquid().flow().ledger_mass_kg().unwrap()
                - liquid.flow().ledger_mass_kg().unwrap();
            let mut gas_mass = 0.0;
            let mut gas_energy = 0.0;
            let mut volatile_mass = 0.0;
            let mut carrier_mass = 0.0;
            for (id, old) in &parcels {
                let new = &next.parcels()[id];
                gas_mass += new.mass_kg() - old.mass_kg();
                gas_energy += new.mass_kg() * new.state().specific_internal_energy_j_kg()
                    - old.mass_kg() * old.state().specific_internal_energy_j_kg();
                volatile_mass += new.species_mass_kg(&water()).unwrap()
                    - old.species_mass_kg(&water()).unwrap();
                carrier_mass += new.species_mass_kg(&carrier()).unwrap()
                    - old.species_mass_kg(&carrier()).unwrap();
                assert_eq!(new.volume_m3(), old.volume_m3());
            }
            near(dm + gas_mass, 0.0, joint_controls().mass_tolerance_kg);
            near(dm + volatile_mass, 0.0, joint_controls().mass_tolerance_kg);
            near(carrier_mass, 0.0, joint_controls().mass_tolerance_kg);
            near(
                next.liquid().total_enthalpy_j().unwrap() - old_h
                    + binding.liquid_enthalpy_offset_j_kg * dm + gas_energy,
                report.coupled.flow.viscous_heat_j + report.coupled.flow.body_heat_j,
                joint_controls().energy_tolerance_j,
            );
            assert!(report.absolute_mass_defect_kg <= joint_controls().mass_tolerance_kg);
            assert!(report.absolute_energy_defect_j <= joint_controls().energy_tolerance_j);
            assert_eq!(parcels, preserved);
            assert_eq!(liquid.flow().lattice().grid().state_bits(), old_bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), old_h);
            assert_eq!(liquid.accepted_steps(), 0);
        });
    }
}

#[test]
fn missing_final_pressure_link_discards_transport_and_empty_paths_retain_native_pressure() {
    with_liquid_interfaces_pressure(0.0, &CELLS, Some(REFERENCE), |liquid, gate| {
        let source = source(-13.4e6);
        let parcels = pressure_parcels(&source);
        let preserved = parcels.clone();
        let old_bits = liquid.flow().lattice().grid().state_bits();
        let old_h = liquid.total_enthalpy_j().unwrap();
        let binding = binding(&source, 0.0);
        let rows = regional_faces(&binding);
        let warm_rows = heated(&rows, 50.0);
        let bounds = endpoints(&parcels);
        let links = [nozzle(1e-10)];
        let paths = GasTransportPaths { nozzles: &links, ..Default::default() };
        let complete = membership();
        let mut missing = complete.clone();
        missing.remove(&[4, 3, 5]);
        assert!(matches!(
            advance_liquid_regions_with_gas_transport_and_pressure(
                liquid, &parcels, &source, paths, &bounds, &rows,
                transport_controls(), joint_controls(),
                RegionalGasPressure {
                    gas_regions: &missing,
                    maximum_density_excursion: MAXIMUM,
                },
                gate,
            ),
            Err(TransportingLiquidGasError::Coupled(_))
        ));
        assert_eq!(parcels, preserved);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), old_bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), old_h);
        assert_eq!(liquid.accepted_steps(), 0);
        assert_eq!(liquid.physical_time_s(), 0.0);
        let pressure = RegionalGasPressure {
            gas_regions: &complete,
            maximum_density_excursion: MAXIMUM,
        };
        let retry = advance_liquid_regions_with_gas_transport_and_pressure(
            liquid, &parcels, &source, paths, &bounds, &rows,
            transport_controls(), joint_controls(), pressure, gate,
        )
        .unwrap();
        assert_eq!(retry.liquid().accepted_steps(), 1);
        assert_eq!(retry.liquid().physical_time_s(), DT);
        let direct = advance_liquid_regions_with_gas_and_pressure(
            liquid, &parcels, &source, &rows, joint_controls(), pressure, gate,
        )
        .unwrap();
        let empty = advance_liquid_regions_with_gas_transport_and_pressure(
            liquid, &parcels, &source, GasTransportPaths::default(), &bounds, &rows,
            transport_controls(), joint_controls(), pressure, gate,
        )
        .unwrap();
        assert_eq!(empty.parcels(), direct.parcels());
        assert_eq!(empty.report().coupled.pressure, direct.report().pressure);
        assert_eq!(
            empty.liquid().flow().lattice().grid().state_bits(),
            direct.liquid().flow().lattice().grid().state_bits()
        );
        assert_eq!(
            empty.liquid().total_enthalpy_j().unwrap(),
            direct.liquid().total_enthalpy_j().unwrap()
        );
        gate.request();
        assert!(matches!(
            advance_liquid_regions_with_gas_transport_and_pressure(
                liquid, &parcels, &source, paths, &bounds, &rows,
                transport_controls(), joint_controls(), pressure, gate,
            ),
            Err(TransportingLiquidGasError::Cancelled)
        ));
        assert!(matches!(
            advance_liquid_regions_with_gas_transport_and_heat_and_pressure(
                liquid, &parcels, &source, paths, &bounds, &warm_rows,
                transport_controls(), joint_controls(),
                GasLiquidChamberHeatSolve { maximum_iterations: 128 },
                pressure, gate,
            ),
            Err(TransportingLiquidGasError::Cancelled)
        ));
        assert_eq!(parcels, preserved);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), old_bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), old_h);
        assert_eq!(liquid.accepted_steps(), 0);
    });
}
