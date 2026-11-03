//! Actual finite-gas plumbing changes native phase transfer and thermal state.
use super::*;
use crate::conduction::{film, joint_controls, thermal_controls, with_wall};
use fs_couple::liquid_gas::evaporation::{
    conduction::advance_liquid_surfaces_with_gas_and_radiative_conduction,
    transport::{TransportingLiquidGasError, advance_liquid_surfaces_with_gas_transport},
};
use fs_lbm::material::free_surface::enthalpy::conduction::FreeSurfaceRadiationFace;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::{
        GasNozzle,
        network::snapshot::{SnapshotNozzleError, SnapshotNozzleLink},
    },
    species::network::{
        SpeciesNetworkSolve,
        transport::{self as plumbing, GasTransportError, GasTransportPaths, GasTransportSolve},
    },
};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::collections::BTreeMap;

fn node(
    source: &ThermochemicalGasModel,
    t: f64,
    pressure: f64,
    vapor_fraction: f64,
    volume: f64,
) -> ClosedGasParcel {
    let original = gas(source, t, 100000.0 * vapor_fraction);
    let point = original
        .point()
        .clone()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            pressure,
        )
        .unwrap();
    ClosedGasParcel::new(original.model().clone(), point, volume).unwrap()
}
fn transport_controls() -> GasTransportSolve {
    GasTransportSolve {
        exchange: SpeciesNetworkSolve {
            maximum_species_fraction: 0.2,
            mass_tolerance_kg: 1e-18,
            energy_tolerance_j: 1e-11,
        },
        entropy_tolerance_j_k: 1e-12,
    }
}
fn endpoints(parcels: &BTreeMap<u64, ClosedGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    parcels
        .keys()
        .map(|&id| (id, controls().phase.gas.gas))
        .collect()
}
fn nozzle(area_m2: f64) -> SnapshotNozzleLink {
    SnapshotNozzleLink {
        id: 31,
        parcels: [1, 0],
        nozzle: GasNozzle {
            area_m2,
            discharge_coefficient: 0.8,
        },
        nozzle_solve: controls().phase.gas.gas,
    }
}
fn radiation() -> FreeSurfaceRadiationFace {
    FreeSurfaceRadiationFace {
        liquid: INTERFACE,
        gas: [4, 3, 3],
        area_m2: 1e-6,
        emissivity: 0.6.into(),
        surroundings_k: 450.0,
    }
}

#[test]
fn g1_original_nozzle_inventory_reaches_phase_radiation_and_finite_walls_once() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        with_wall(|body, contact| {
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
            let surfaces = [surface(&Phase, &binding)];
            let bounds = endpoints(&parcels);
            let links = [nozzle(1e-10)];
            let paths = GasTransportPaths {
                nozzles: &links,
                ..Default::default()
            };
            let bodies = [body];
            let old_wall_h = bodies[0].enthalpy_j();
            let exposed = [film(INTERFACE)];
            let radiant = [radiation()];
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
            let original_coupling = |gas: &ClosedGasParcel| {
                advance_liquid_surfaces_with_gas_and_radiative_conduction(
                    liquid,
                    gas,
                    &source,
                    &surfaces,
                    &bodies,
                    &[],
                    &[],
                    &[contact],
                    &[],
                    &exposed,
                    &radiant,
                    joint_controls(),
                    thermal_controls(),
                    controls().phase.gas.gas,
                    None,
                    false,
                    gate,
                )
                .unwrap()
            };
            let direct = original_coupling(&moved.parcels()[&0]);
            let isolated = original_coupling(&parcels[&0]);
            let next = advance_liquid_surfaces_with_gas_transport(
                liquid,
                &parcels,
                0,
                &source,
                paths,
                &bounds,
                &surfaces,
                &bodies,
                &[],
                &[],
                &[contact],
                &[],
                &exposed,
                &radiant,
                transport_controls(),
                joint_controls(),
                thermal_controls(),
                None,
                false,
                gate,
            )
            .unwrap();
            assert_eq!(&next.report().transport, moved.step());
            assert_eq!(next.gas(), direct.gas());
            assert_eq!(next.parcels()[&1], moved.parcels()[&1]);
            assert_eq!(next.parcels()[&9], parcels[&9]);
            assert_eq!(
                next.liquid().flow().lattice().grid().state_bits(),
                direct.liquid().flow().lattice().grid().state_bits()
            );
            assert_eq!(
                next.liquid().total_enthalpy_j().unwrap(),
                direct.liquid().total_enthalpy_j().unwrap()
            );
            assert_eq!(
                next.bodies()[0].enthalpy_j(),
                direct.bodies()[0].enthalpy_j()
            );
            assert_ne!(
                next.report().coupled.phase.phase.transferred_mass_kg,
                isolated.report().phase.phase.transferred_mass_kg
            );
            assert!(core::ptr::eq(next.catalog(), &source));
            assert_eq!(next.report().duration_s, DT);
            assert_eq!(next.liquid().physical_time_s(), DT);
            assert_eq!(next.liquid().accepted_steps(), 1);

            let dm = next.liquid().flow().ledger_mass_kg().unwrap()
                - liquid.flow().ledger_mass_kg().unwrap();
            let dh = next.liquid().total_enthalpy_j().unwrap() - old_h;
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
                carrier_mass += new.species_mass_kg(&carrier()).unwrap()
                    - old.species_mass_kg(&carrier()).unwrap();
                assert_eq!(new.volume_m3(), old.volume_m3());
            }
            let thermal = &next.report().coupled.thermal.conduction;
            let liquid_radiation = thermal.radiation.as_ref().unwrap().heat_to_liquid_j;
            let wall_radiation = thermal.body_radiation.as_ref().unwrap().heat_to_solids_j;
            assert!(liquid_radiation > 0.0 && wall_radiation > 0.0);
            near(dm + gas_mass, 0.0, joint_controls().mass_tolerance_kg);
            near(dm + volatile_mass, 0.0, joint_controls().mass_tolerance_kg);
            near(carrier_mass, 0.0, joint_controls().mass_tolerance_kg);
            near(
                dh + binding.liquid_enthalpy_offset_j_kg * dm
                    + gas_energy
                    + next.bodies()[0].enthalpy_j()
                    - old_wall_h,
                liquid_radiation + wall_radiation,
                joint_controls().energy_tolerance_j,
            );
            assert!(next.report().absolute_mass_defect_kg <= joint_controls().mass_tolerance_kg);
            assert!(next.report().absolute_energy_defect_j <= joint_controls().energy_tolerance_j);
            assert_eq!(parcels, original);
            assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
            assert_eq!(liquid.accepted_steps(), 0);
            assert_eq!(bodies[0].enthalpy_j(), old_wall_h);
        })
    });
}

#[test]
fn g3_post_phase_thermal_pressure_reversal_refuses_the_complete_network_atomically() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        with_wall(|body, contact| {
            let source = source(-13.4e6);
            let parcels = BTreeMap::from([
                (0, node(&source, 295.0, 100000.0, 0.01, 1e-8)),
                (1, node(&source, 295.0, 100000.01, 0.01, 1e-6)),
            ]);
            let original = parcels.clone();
            let bits = liquid.flow().lattice().grid().state_bits();
            let h = liquid.total_enthalpy_j().unwrap();
            let binding = binding(&source, 0.0);
            let surfaces = [surface(&Phase, &binding)];
            let bodies = [body];
            let old_wall_h = bodies[0].enthalpy_j();
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
            let direct = advance_liquid_surfaces_with_gas_and_radiative_conduction(
                liquid,
                &moved.parcels()[&0],
                &source,
                &surfaces,
                &bodies,
                &[],
                &[],
                &[contact],
                &[],
                &[film(INTERFACE)],
                &[],
                joint_controls(),
                thermal_controls(),
                controls().phase.gas.gas,
                None,
                false,
                gate,
            )
            .unwrap();
            assert!(
                direct.gas().state().state().pressure
                    > moved.parcels()[&1].state().state().pressure
            );
            let run = |paths| {
                advance_liquid_surfaces_with_gas_transport(
                    liquid,
                    &parcels,
                    0,
                    &source,
                    paths,
                    &bounds,
                    &surfaces,
                    &bodies,
                    &[],
                    &[],
                    &[contact],
                    &[],
                    &[film(INTERFACE)],
                    &[],
                    transport_controls(),
                    joint_controls(),
                    thermal_controls(),
                    None,
                    false,
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
            assert_eq!(liquid.total_enthalpy_j().unwrap(), h);
            assert_eq!(liquid.accepted_steps(), 0);
            assert_eq!(bodies[0].enthalpy_j(), old_wall_h);
            // Closing the physical restriction permits the same thermal/phase
            // candidate; no failed flow, phase debit or wall heat survived.
            let closed = [nozzle(0.0)];
            let valid = run(GasTransportPaths {
                nozzles: &closed,
                ..Default::default()
            })
            .unwrap();
            assert_eq!(valid.liquid().accepted_steps(), 1);
            assert_eq!(valid.parcels()[&1], parcels[&1]);
            assert!(valid.gas().state().state().pressure > parcels[&1].state().state().pressure);
        })
    });
}

#[test]
fn g3_empty_plumbing_recovers_original_radiative_phase_conduction_and_retains_idle_nodes() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let source = source(-13.4e6);
        let parcels = BTreeMap::from([
            (0, gas(&source, 300.0, 1000.0)),
            (9, gas(&source, 305.0, 3000.0)),
        ]);
        let bounds = endpoints(&parcels);
        let binding = binding(&source, 0.0);
        let surfaces = [surface(&Phase, &binding)];
        let exposed = [film(INTERFACE)];
        let radiant = [radiation()];
        let direct = advance_liquid_surfaces_with_gas_and_radiative_conduction(
            liquid,
            &parcels[&0],
            &source,
            &surfaces,
            &[],
            &[],
            &[],
            &[],
            &[],
            &exposed,
            &radiant,
            joint_controls(),
            thermal_controls(),
            controls().phase.gas.gas,
            None,
            false,
            gate,
        )
        .unwrap();
        let next = advance_liquid_surfaces_with_gas_transport(
            liquid,
            &parcels,
            0,
            &source,
            GasTransportPaths::default(),
            &bounds,
            &surfaces,
            &[],
            &[],
            &[],
            &[],
            &[],
            &exposed,
            &radiant,
            transport_controls(),
            joint_controls(),
            thermal_controls(),
            None,
            false,
            gate,
        )
        .unwrap();
        assert_eq!(next.gas(), direct.gas());
        assert_eq!(
            next.liquid().flow().lattice().grid().state_bits(),
            direct.liquid().flow().lattice().grid().state_bits()
        );
        assert_eq!(
            next.liquid().total_enthalpy_j().unwrap(),
            direct.liquid().total_enthalpy_j().unwrap()
        );
        assert_eq!(
            next.report().transport.exchange.total_energy_residual_j,
            0.0
        );
        assert_eq!(
            next.report().energy_residual_j,
            direct.report().energy_residual_j
        );
        assert_eq!(next.parcels()[&9], parcels[&9]);
        assert_eq!(next.liquid().accepted_steps(), 1);
        let (_, complete, catalog, walls, report) = next.into_parts();
        assert_eq!(complete.len(), 2);
        assert!(core::ptr::eq(catalog, &source));
        assert!(walls.is_empty());
        assert_eq!(report.atmosphere_id, 0);
        assert_eq!(report.duration_s, DT);
    });
}
