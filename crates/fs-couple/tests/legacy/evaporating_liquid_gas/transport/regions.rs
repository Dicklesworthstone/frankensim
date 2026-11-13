//! Original gas transport feeds local native phase and sensible-heat owners.
use super::*;

#[path = "regions/pressure.rs"]
mod pressure;
use fs_couple::liquid_gas::evaporation::{
    multiple::{
        LiquidGasHeatSurface,
        regions::{
            RegionalLiquidGasHeatSurface, RegionalLiquidGasSurface,
            advance_liquid_regions_with_gas, advance_liquid_regions_with_gas_and_heat,
        },
    },
    transport::regions::{
        advance_liquid_regions_with_gas_transport,
        advance_liquid_regions_with_gas_transport_and_heat,
    },
};
use fs_material::evaporation::gas::chamber::heat::GasLiquidChamberHeatSolve;

const CELLS: [([usize; 3], f64, f64); 2] = [([3, 2, 5], 0.4, 295.0), ([3, 4, 2], 0.4, 305.0)];
fn regional_faces(binding: &GasLiquidBinding) -> [RegionalLiquidGasSurface<'_>; 2] {
    [0, 1].map(|index| {
        let coordinates = CELLS[index].0;
        RegionalLiquidGasSurface {
            region_id: index as u64,
            surface: LiquidGasEvaporationSurface {
                liquid: coordinates,
                gas: [coordinates[0] + 1, coordinates[1], coordinates[2]],
                ..surface(&Phase, binding)
            },
        }
    })
}
fn heated<'a>(
    rows: &[RegionalLiquidGasSurface<'a>],
    coefficient: f64,
) -> Vec<RegionalLiquidGasHeatSurface<'a>> {
    rows.iter()
        .map(|row| RegionalLiquidGasHeatSurface {
            region_id: row.region_id,
            surface: LiquidGasHeatSurface {
                surface: row.surface,
                heat_transfer_coefficient_w_m2_k: coefficient,
                heat_controls: heat_controls(),
            },
        })
        .collect()
}

#[test]
fn nozzle_transport_changes_local_phase_heat_and_closes_all_original_inventories_once() {
    with_liquid_interfaces(0.0, &CELLS, |liquid, gate| {
        let source = source(-13.4e6);
        let parcels = BTreeMap::from([
            (0, node(&source, 300.0, 100000.0, 0.01, 1e-8)),
            (1, node(&source, 300.0, 120000.0, 0.5, 1e-6)),
            (9, node(&source, 300.0, 90000.0, 0.05, 1e-8)),
        ]);
        let original = parcels.clone();
        let bits = liquid.flow().lattice().grid().state_bits();
        let old_h = liquid.total_enthalpy_j().unwrap();
        let binding = binding(&source, 0.0);
        let rows = heated(&regional_faces(&binding), 20.0);
        let bounds = endpoints(&parcels);
        let links = [nozzle(1e-10)];
        let paths = GasTransportPaths {
            nozzles: &links,
            ..Default::default()
        };
        let moved = plumbing::step_candidate_with_paths(
            &parcels,
            paths,
            DT,
            &bounds,
            transport_controls(),
            || false,
        )
        .unwrap();
        assert!(moved.step().nozzles[0].transferred_mass_kg > 0.0);
        assert!(
            moved.parcels()[&0].species_mass_kg(&water()).unwrap()
                > parcels[&0].species_mass_kg(&water()).unwrap()
        );
        let isolated = advance_liquid_regions_with_gas_and_heat(
            liquid,
            &parcels,
            &source,
            &rows,
            joint_controls(),
            GasLiquidChamberHeatSolve {
                maximum_iterations: 128,
            },
            gate,
        )
        .unwrap();
        let next = advance_liquid_regions_with_gas_transport_and_heat(
            liquid,
            &parcels,
            &source,
            paths,
            &bounds,
            &rows,
            transport_controls(),
            joint_controls(),
            GasLiquidChamberHeatSolve {
                maximum_iterations: 128,
            },
            gate,
        )
        .unwrap();
        let report = next.report();
        assert_eq!(report.duration_s, DT);
        assert_eq!(report.transport.duration_s, DT);
        assert_eq!(report.coupled.duration_s, DT);
        assert_eq!(next.liquid().accepted_steps(), 1);
        assert_eq!(next.liquid().physical_time_s(), DT);
        assert_eq!(report.coupled.mass_exchange.samples.len(), 2);
        assert_eq!(next.parcels()[&9], parcels[&9]);
        assert!(std::ptr::eq(next.catalog(), &source));
        for id in [0, 1] {
            let region = &report.coupled.regions[&id];
            let heat = region.heat.as_ref().unwrap();
            assert_eq!(
                heat.gas_temperature_k[0],
                moved.parcels()[&id].state().state().temperature
            );
            assert_eq!(heat.gas_temperature_k[1], region.phase.gas_temperature_k[0]);
            assert_eq!(region.exchange.transferred_mass_kg > 0.0, id == 0);
            assert_eq!(heat.heat_into_liquid_j > 0.0, id == 0);
            near(
                next.parcels()[&id].species_mass_kg(&water()).unwrap()
                    - moved.parcels()[&id].species_mass_kg(&water()).unwrap(),
                region.exchange.transferred_mass_kg,
                controls().phase.gas.mass_tolerance_kg,
            );
        }
        assert_ne!(
            report.coupled.regions[&0].exchange.transferred_mass_kg,
            isolated.report().regions[&0].exchange.transferred_mass_kg
        );
        assert!(
            report.coupled.regions[&0].surfaces[0]
                .phase
                .vapor_partial_pressure_pa[0]
                > isolated.report().regions[&0].surfaces[0]
                    .phase
                    .vapor_partial_pressure_pa[0]
        );
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
            volatile_mass +=
                new.species_mass_kg(&water()).unwrap() - old.species_mass_kg(&water()).unwrap();
            carrier_mass +=
                new.species_mass_kg(&carrier()).unwrap() - old.species_mass_kg(&carrier()).unwrap();
            assert_eq!(new.volume_m3(), old.volume_m3());
        }
        near(dm + gas_mass, 0.0, joint_controls().mass_tolerance_kg);
        near(dm + volatile_mass, 0.0, joint_controls().mass_tolerance_kg);
        near(carrier_mass, 0.0, joint_controls().mass_tolerance_kg);
        near(
            next.liquid().total_enthalpy_j().unwrap() - old_h
                + binding.liquid_enthalpy_offset_j_kg * dm
                + gas_energy,
            report.coupled.flow.external_heat_j
                + report.coupled.flow.viscous_heat_j
                + report.coupled.flow.body_heat_j,
            joint_controls().energy_tolerance_j,
        );
        assert!(report.absolute_mass_defect_kg <= joint_controls().mass_tolerance_kg);
        assert!(report.absolute_energy_defect_j <= joint_controls().energy_tolerance_j);
        assert_eq!(parcels, original);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), old_h);
        assert_eq!(liquid.accepted_steps(), 0);
    });
}

#[test]
fn final_regional_heat_pressure_reversal_refuses_all_nodes_and_a_closed_nozzle_retries() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let source = source(-13.4e6);
        let parcels = BTreeMap::from([
            (0, node(&source, 295.0, 100000.0, 0.01, 1e-8)),
            (1, node(&source, 295.0, 100000.01, 0.01, 1e-6)),
        ]);
        let original = parcels.clone();
        let bits = liquid.flow().lattice().grid().state_bits();
        let old_h = liquid.total_enthalpy_j().unwrap();
        let binding = binding(&source, 0.0);
        let rows = heated(
            &[RegionalLiquidGasSurface {
                region_id: 0,
                surface: surface(&Phase, &binding),
            }],
            100.0,
        );
        let bounds = endpoints(&parcels);
        let links = [nozzle(1e-13)];
        let paths = GasTransportPaths {
            nozzles: &links,
            ..Default::default()
        };
        let moved = plumbing::step_candidate_with_paths(
            &parcels,
            paths,
            DT,
            &bounds,
            transport_controls(),
            || false,
        )
        .unwrap();
        assert!(moved.step().nozzles[0].transferred_mass_kg > 0.0);
        let direct = advance_liquid_regions_with_gas_and_heat(
            liquid,
            moved.parcels(),
            &source,
            &rows,
            joint_controls(),
            GasLiquidChamberHeatSolve {
                maximum_iterations: 128,
            },
            gate,
        )
        .unwrap();
        assert!(
            direct.parcels()[&0].state().state().pressure
                > moved.parcels()[&1].state().state().pressure
        );
        assert!((direct.parcels()[&0].state().state().temperature - 295.0).abs() < 1.0);
        let run = |paths| {
            advance_liquid_regions_with_gas_transport_and_heat(
                liquid,
                &parcels,
                &source,
                paths,
                &bounds,
                &rows,
                transport_controls(),
                joint_controls(),
                GasLiquidChamberHeatSolve {
                    maximum_iterations: 128,
                },
                gate,
            )
        };
        assert!(matches!(
            run(paths),
            Err(TransportingLiquidGasError::Transport(
                GasTransportError::Nozzle(SnapshotNozzleError::EquilibriumCrossed { id: 31 })
            ))
        ));
        assert_eq!(parcels, original);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), old_h);
        assert_eq!(liquid.accepted_steps(), 0);
        assert_eq!(liquid.physical_time_s(), 0.0);
        let closed = [nozzle(0.0)];
        let retry = run(GasTransportPaths {
            nozzles: &closed,
            ..Default::default()
        })
        .unwrap();
        assert_eq!(retry.liquid().accepted_steps(), 1);
        assert_eq!(retry.liquid().physical_time_s(), DT);
        assert_eq!(retry.parcels()[&1], parcels[&1]);
        assert!(
            retry.parcels()[&0].state().state().pressure > parcels[&1].state().state().pressure
        );
    });
}

#[test]
fn empty_paths_preserve_direct_regional_phase_and_heat_endpoints_and_cancellation_is_atomic() {
    with_liquid_interfaces(0.0, &CELLS, |liquid, gate| {
        let source = source(-13.4e6);
        let parcels = BTreeMap::from([
            (0, gas(&source, 300.0, 1000.0)),
            (1, gas(&source, 310.0, 6000.0)),
            (9, gas(&source, 305.0, 3000.0)),
        ]);
        let original = parcels.clone();
        let bits = liquid.flow().lattice().grid().state_bits();
        let old_h = liquid.total_enthalpy_j().unwrap();
        let binding = binding(&source, 0.0);
        let rows = regional_faces(&binding);
        let bounds = endpoints(&parcels);
        let direct = advance_liquid_regions_with_gas(
            liquid,
            &parcels,
            &source,
            &rows,
            joint_controls(),
            gate,
        )
        .unwrap();
        let through = advance_liquid_regions_with_gas_transport(
            liquid,
            &parcels,
            &source,
            GasTransportPaths::default(),
            &bounds,
            &rows,
            transport_controls(),
            joint_controls(),
            gate,
        )
        .unwrap();
        assert_eq!(through.parcels(), direct.parcels());
        assert_eq!(
            through.liquid().flow().lattice().grid().state_bits(),
            direct.liquid().flow().lattice().grid().state_bits()
        );
        assert_eq!(
            through.liquid().total_enthalpy_j().unwrap(),
            direct.liquid().total_enthalpy_j().unwrap()
        );
        assert_eq!(
            through.report().coupled.mass_exchange,
            direct.report().mass_exchange
        );
        assert_eq!(
            through.report().transport.exchange.total_energy_residual_j,
            0.0
        );
        let warm_rows = heated(&rows, 20.0);
        let direct_heat = advance_liquid_regions_with_gas_and_heat(
            liquid,
            &parcels,
            &source,
            &warm_rows,
            joint_controls(),
            GasLiquidChamberHeatSolve {
                maximum_iterations: 128,
            },
            gate,
        )
        .unwrap();
        let through_heat = advance_liquid_regions_with_gas_transport_and_heat(
            liquid,
            &parcels,
            &source,
            GasTransportPaths::default(),
            &bounds,
            &warm_rows,
            transport_controls(),
            joint_controls(),
            GasLiquidChamberHeatSolve {
                maximum_iterations: 128,
            },
            gate,
        )
        .unwrap();
        assert_eq!(through_heat.parcels(), direct_heat.parcels());
        assert_eq!(
            through_heat.liquid().flow().lattice().grid().state_bits(),
            direct_heat.liquid().flow().lattice().grid().state_bits()
        );
        assert_eq!(
            through_heat.liquid().total_enthalpy_j().unwrap(),
            direct_heat.liquid().total_enthalpy_j().unwrap()
        );
        for id in [0, 1] {
            assert_eq!(
                through.report().coupled.regions[&id].exchange,
                direct.report().regions[&id].exchange
            );
            assert_eq!(
                through_heat.report().coupled.regions[&id].exchange,
                direct_heat.report().regions[&id].exchange
            );
            assert_eq!(
                through_heat.report().coupled.regions[&id].heat,
                direct_heat.report().regions[&id].heat
            );
        }
        assert_eq!(through_heat.parcels()[&9], parcels[&9]);
        assert_eq!(through_heat.liquid().accepted_steps(), 1);
        assert_eq!(through_heat.liquid().physical_time_s(), DT);
        gate.request();
        assert!(matches!(
            advance_liquid_regions_with_gas_transport(
                liquid,
                &parcels,
                &source,
                GasTransportPaths::default(),
                &bounds,
                &rows,
                transport_controls(),
                joint_controls(),
                gate,
            ),
            Err(TransportingLiquidGasError::Cancelled)
        ));
        assert!(matches!(
            advance_liquid_regions_with_gas_transport_and_heat(
                liquid,
                &parcels,
                &source,
                GasTransportPaths::default(),
                &bounds,
                &warm_rows,
                transport_controls(),
                joint_controls(),
                GasLiquidChamberHeatSolve {
                    maximum_iterations: 128
                },
                gate,
            ),
            Err(TransportingLiquidGasError::Cancelled)
        ));
        assert_eq!(parcels, original);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), old_h);
        assert_eq!(liquid.accepted_steps(), 0);
    });
}
