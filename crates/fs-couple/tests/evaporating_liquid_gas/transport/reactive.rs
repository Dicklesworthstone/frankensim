//! Original gas plumbing, chemistry, phase and heat share one native clock.
use super::*;
use crate::reactive::{reaction, reactive_controls};
use fs_couple::liquid_gas::evaporation::reactive::ReactiveLiquidPhaseError;
use fs_lbm::material::free_surface::enthalpy::FreeSurfaceCellHeat;

#[test]
fn g1_transported_species_react_and_exchange_phase_heat_with_one_complete_inventory() {
    for maximum_density_excursion in [None, Some(0.3)] {
        with_liquid_pressure(
            0.0,
            0.4,
            maximum_density_excursion.map(|_| 100000.0),
            |liquid, gate| {
                with_wall(|body, contact| {
                    let law = reaction(500.0, 75.0);
                    let source = law.catalog();
                    let parcels = BTreeMap::from([
                        (0, node(source, 300.0, 100000.0, 0.01, 1e-8)),
                        (1, node(source, 300.0, 120000.0, 0.5, 1e-6)),
                        (9, node(source, 300.0, 90000.0, 0.05, 1e-8)),
                    ]);
                    let original = parcels.clone();
                    let old_bits = liquid.flow().lattice().grid().state_bits();
                    let old_h = liquid.total_enthalpy_j().unwrap();
                    let binding = binding(source, 0.0);
                    let surfaces = [surface(&Phase, &binding)];
                    let bodies = [body];
                    let old_wall_h = bodies[0].enthalpy_j();
                    let exposed = [film(INTERFACE)];
                    let radiant = [radiation()];
                    let heat = [FreeSurfaceCellHeat {
                        coordinates: INTERFACE,
                        energy_j: 1e-6,
                    }];
                    let bounds = endpoints(&parcels);
                    let links = [nozzle(1e-12)];
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
                    let direct = law
                        .radiative_conducting_evaporating_liquid_candidate(
                            liquid,
                            &moved.parcels()[&0],
                            &surfaces,
                            &bodies,
                            &[],
                            &heat,
                            &[contact],
                            &[],
                            &exposed,
                            &radiant,
                            joint_controls(),
                            thermal_controls(),
                            maximum_density_excursion,
                            false,
                            reactive_controls(),
                            gate,
                        )
                        .unwrap();
                    let isolated_reaction = law
                        .parcel_candidate(&parcels[&0], DT, reactive_controls().reaction, || false)
                        .unwrap();
                    let next = law
                        .transporting_evaporating_liquid_candidate(
                            liquid,
                            &parcels,
                            0,
                            paths,
                            &bounds,
                            &surfaces,
                            &bodies,
                            &[],
                            &heat,
                            &[contact],
                            &[],
                            &exposed,
                            &radiant,
                            transport_controls(),
                            joint_controls(),
                            thermal_controls(),
                            maximum_density_excursion,
                            false,
                            reactive_controls(),
                            gate,
                        )
                        .unwrap();
                    let report = next.report();
                    assert_eq!(&report.transport, moved.step());
                    assert_eq!(report.coupled.reaction, direct.report().reaction);
                    assert!(report.coupled.reaction.extent_mol > 0.0);
                    assert_ne!(
                        report.coupled.reaction.extent_mol,
                        isolated_reaction.step().extent_mol,
                        "the original kinetic law must consume the transported composition"
                    );
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
                    assert!(core::ptr::eq(next.catalog(), source));
                    assert_eq!(report.duration_s, DT);
                    assert_eq!(report.coupled.duration_s, DT);
                    assert_eq!(next.liquid().physical_time_s(), DT);
                    assert_eq!(next.liquid().accepted_steps(), 1);

                    // Independent inventory differences include EVERY gas node.
                    // Chemistry changes constituents; only the original reaction
                    // stoichiometry can supply those changes. Film/contact heat
                    // cancels inside the enclosing inventory, and reaction heat is
                    // already in the original chemical formation energies.
                    let dm = next.liquid().flow().ledger_mass_kg().unwrap()
                        - liquid.flow().ledger_mass_kg().unwrap();
                    let dh = next.liquid().total_enthalpy_j().unwrap() - old_h;
                    let mut gas_mass = 0.0;
                    let mut gas_energy = 0.0;
                    let mut vapor_mass = 0.0;
                    let mut carrier_mass = 0.0;
                    for (id, old) in &parcels {
                        let new = &next.parcels()[id];
                        gas_mass += new.mass_kg() - old.mass_kg();
                        gas_energy += new.mass_kg() * new.state().specific_internal_energy_j_kg()
                            - old.mass_kg() * old.state().specific_internal_energy_j_kg();
                        vapor_mass += new.species_mass_kg(&water()).unwrap()
                            - old.species_mass_kg(&water()).unwrap();
                        carrier_mass += new.species_mass_kg(&carrier()).unwrap()
                            - old.species_mass_kg(&carrier()).unwrap();
                        assert_eq!(new.volume_m3(), old.volume_m3());
                    }
                    let chemistry = &report.coupled.reaction;
                    let vapor = chemistry
                        .species
                        .iter()
                        .find(|row| row.species == water())
                        .unwrap();
                    let carrier = chemistry
                        .species
                        .iter()
                        .find(|row| row.species == carrier())
                        .unwrap();
                    assert!(vapor.stoichiometric_change_mol > 0.0);
                    assert!(carrier.stoichiometric_change_mol < 0.0);
                    near(dm + gas_mass, 0.0, reactive_controls().mass_tolerance_kg);
                    near(
                        dm + vapor_mass,
                        M * vapor.stoichiometric_change_mol,
                        reactive_controls().mass_tolerance_kg,
                    );
                    near(
                        carrier_mass,
                        M * carrier.stoichiometric_change_mol,
                        reactive_controls().mass_tolerance_kg,
                    );
                    let thermal = &report.coupled.phase.thermal.conduction;
                    let radiation_j = thermal.radiation.as_ref().unwrap().heat_to_liquid_j
                        + thermal.body_radiation.as_ref().unwrap().heat_to_solids_j;
                    assert!(radiation_j > 0.0);
                    near(
                        dh + binding.liquid_enthalpy_offset_j_kg * dm
                            + gas_energy
                            + next.bodies()[0].enthalpy_j()
                            - old_wall_h,
                        heat[0].energy_j + radiation_j,
                        reactive_controls().energy_tolerance_j,
                    );
                    assert!(
                        report.absolute_mass_defect_kg <= reactive_controls().mass_tolerance_kg
                    );
                    assert!(
                        report.absolute_energy_defect_j <= reactive_controls().energy_tolerance_j
                    );
                    if let Some(maximum) = maximum_density_excursion {
                        let pressure = thermal.finite_gas.as_ref().unwrap().pressure.unwrap();
                        assert_eq!(
                            pressure.applied_pressure_pa,
                            next.gas().state().state().pressure
                        );
                        assert_eq!(pressure.maximum_density_excursion, maximum);
                    }
                    assert_eq!(parcels, original);
                    assert_eq!(liquid.flow().lattice().grid().state_bits(), old_bits);
                    assert_eq!(liquid.total_enthalpy_j().unwrap(), old_h);
                    assert_eq!(liquid.accepted_steps(), 0);
                    assert_eq!(bodies[0].enthalpy_j(), old_wall_h);
                })
            },
        );
    }
}

#[test]
fn g3_chemical_pressure_reversal_rechecks_the_final_network_and_rolls_back() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let law = reaction(500.0, 75.0);
        let source = law.catalog();
        let parcels = BTreeMap::from([
            (0, node(source, 300.0, 100000.0, 0.01, 1e-8)),
            (1, node(source, 300.0, 100000.01, 0.01, 1e-6)),
        ]);
        let original = parcels.clone();
        let old_bits = liquid.flow().lattice().grid().state_bits();
        let old_h = liquid.total_enthalpy_j().unwrap();
        let binding = binding(source, 0.0);
        let mut patch = surface(&Phase, &binding);
        // Shut phase transfer and direct gas heat to isolate chemical pressure.
        patch.mass_transfer_coefficient_m_s = 0.0;
        let surfaces = [patch];
        let mut reflection = radiation();
        reflection.emissivity = 0.0.into();
        let radiant = [reflection];
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
        let direct = law
            .radiative_conducting_evaporating_liquid_candidate(
                liquid,
                &moved.parcels()[&0],
                &surfaces,
                &[],
                &[],
                &[],
                &[],
                &[],
                &[],
                &radiant,
                joint_controls(),
                thermal_controls(),
                None,
                false,
                reactive_controls(),
                gate,
            )
            .unwrap();
        assert!(direct.report().reaction.extent_mol > 0.0);
        assert!(
            direct.gas().state().state().pressure > moved.parcels()[&1].state().state().pressure
        );
        let run = |paths| {
            law.transporting_evaporating_liquid_candidate(
                liquid,
                &parcels,
                0,
                paths,
                &bounds,
                &surfaces,
                &[],
                &[],
                &[],
                &[],
                &[],
                &[],
                &radiant,
                transport_controls(),
                joint_controls(),
                thermal_controls(),
                None,
                false,
                reactive_controls(),
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
        assert_eq!(liquid.flow().lattice().grid().state_bits(), old_bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), old_h);
        assert_eq!(liquid.accepted_steps(), 0);
        let closed = [nozzle(0.0)];
        let valid = run(GasTransportPaths {
            nozzles: &closed,
            ..Default::default()
        })
        .unwrap();
        assert_eq!(valid.liquid().accepted_steps(), 1);
        assert_eq!(valid.parcels()[&1], parcels[&1]);
        assert!(valid.gas().state().state().pressure > parcels[&1].state().state().pressure);
    });
}

#[test]
fn g3_empty_paths_preserve_reactive_owner_and_its_explicit_caloric_controls() {
    with_liquid(0.0, 0.4, |liquid, gate| {
        let law = reaction(500.0, 75.0);
        let source = law.catalog();
        let parcels = BTreeMap::from([
            (0, gas(source, 300.0, 1000.0)),
            (9, gas(source, 305.0, 3000.0)),
        ]);
        let original = parcels.clone();
        let old_bits = liquid.flow().lattice().grid().state_bits();
        let binding = binding(source, 0.0);
        let surfaces = [surface(&Phase, &binding)];
        let exposed = [film(INTERFACE)];
        let radiant = [radiation()];
        let bounds = endpoints(&parcels);
        let direct = law
            .radiative_conducting_evaporating_liquid_candidate(
                liquid,
                &parcels[&0],
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
                None,
                false,
                reactive_controls(),
                gate,
            )
            .unwrap();
        let run = |controls| {
            law.transporting_evaporating_liquid_candidate(
                liquid,
                &parcels,
                0,
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
                controls,
                gate,
            )
        };
        let next = run(reactive_controls()).unwrap();
        assert_eq!(next.gas(), direct.gas());
        assert_eq!(next.report().coupled.reaction, direct.report().reaction);
        assert_eq!(
            next.liquid().flow().lattice().grid().state_bits(),
            direct.liquid().flow().lattice().grid().state_bits()
        );
        assert_eq!(
            next.liquid().total_enthalpy_j().unwrap(),
            direct.liquid().total_enthalpy_j().unwrap()
        );
        assert_eq!(
            next.report().energy_residual_j,
            direct.report().energy_residual_j
        );
        assert_eq!(
            next.report().absolute_energy_defect_j,
            direct.report().absolute_energy_defect_j
        );
        assert_eq!(next.parcels()[&9], parcels[&9]);
        let mut restricted = reactive_controls();
        restricted.reaction.gas.maximum_temperature_k = 300.0;
        assert!(
            matches!(
                run(restricted),
                Err(TransportingLiquidGasError::Reactive(
                    ReactiveLiquidPhaseError::Reaction(_)
                ))
            ),
            "network endpoint tolerances must not overwrite the declared chemical caloric bracket"
        );
        assert_eq!(parcels, original);
        assert_eq!(liquid.flow().lattice().grid().state_bits(), old_bits);
        assert_eq!(liquid.accepted_steps(), 0);
        let (_, complete, catalog, bodies, report) = next.into_parts();
        assert_eq!(complete.len(), 2);
        assert!(core::ptr::eq(catalog, source));
        assert!(bodies.is_empty());
        assert_eq!(report.duration_s, DT);
        assert_eq!(report.atmosphere_id, 0);
    });
}
