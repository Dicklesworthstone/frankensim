//! Manufactured compartment equations on the original thermal fixtures.
use super::*;
use crate::burnup::production::{BurnupProductionChange, BurnupProductionSchedule};
use crate::burnup::transport::{
    BurnupTransportBudget, BurnupTransportLimits, BurnupTransportState, NuclideTransfer,
    NuclideTransport,
};
use crate::burnup::{
    BurnupData, BurnupInventory, BurnupLimits, BurnupModel, BurnupNuclide, BurnupProduct,
    BurnupRate, BurnupReaction,
};
use crate::cross_sections::{
    CrossSectionData, CrossSectionLimits, CrossSectionTable, NeutronReaction,
};
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

fn transport_limits() -> BurnupTransportLimits {
    BurnupTransportLimits {
        burnup: BurnupLimits::default(),
        max_compartments: 4,
        max_transfers: 4,
        max_outflows: 2,
        max_rate_coefficients: 32,
    }
}

fn transport_budget() -> BurnupTransportBudget {
    BurnupTransportBudget {
        limits: transport_limits(),
        nuclei_tolerance: 1e-7,
        energy_tolerance_j: 1e-7,
    }
}

fn flow(rate: f64) -> NuclideTransport {
    NuclideTransport {
        transfers: vec![NuclideTransfer {
            from_compartment: 0,
            to_compartment: 1,
            rates_per_s: vec![rate, rate],
        }],
        outflows: vec![],
    }
}

fn inventory_model(capture: Option<(f64, f64)>, maximum_k: f64) -> Arc<BurnupModel> {
    let groups = Arc::new(
        NeutronGroups::new("manufactured-mobile-isotope-grid".into(), vec![0.0, 1.0]).unwrap(),
    );
    let parent = NuclideId::new(10, 25, 0).unwrap();
    let energy = EnergyChannels {
        deposited_j: 8.0,
        radiation_j: 2.0,
        escaped_j: 1.0,
    };
    let products = vec![BurnupProduct {
        nuclide: 1,
        yield_per_event: 1.0,
    }];
    let mut reactions = vec![BurnupReaction {
        data_id: "manufactured-moving-parent-decay".into(),
        parent: 0,
        rate: BurnupRate::Decay {
            branch_fraction: 1.0,
        },
        products: products.clone(),
        energy_per_event: energy,
    }];
    if let Some((lower, upper)) = capture {
        reactions.push(BurnupReaction {
            data_id: "manufactured-moving-parent-capture".into(),
            parent: 0,
            rate: BurnupRate::TabulatedCapture {
                cross_sections: Arc::new(
                    CrossSectionTable::new(
                        CrossSectionData {
                            data_id: "manufactured-compartment-temperature-fit".into(),
                            nuclide: parent,
                            reaction: NeutronReaction::Capture,
                            groups: Arc::clone(&groups),
                            temperatures_k: vec![300.0, maximum_k],
                            cross_sections_m2: vec![vec![lower], vec![upper]],
                        },
                        CrossSectionLimits::default(),
                    )
                    .unwrap(),
                ),
            },
            products,
            energy_per_event: energy,
        });
    }
    Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-compartment-decay-and-capture".into(),
                groups,
                temperature_min_k: 300.0,
                temperature_max_k: maximum_k,
                nuclides: vec![
                    BurnupNuclide {
                        id: parent,
                        law: DecayLawV1::from_decay_constant(Frequency::new(0.5)).unwrap(),
                    },
                    BurnupNuclide {
                        id: NuclideId::new(11, 25, 0).unwrap(),
                        law: DecayLawV1::stable(),
                    },
                ],
                reactions,
            },
            transport_limits().burnup,
        )
        .unwrap(),
    )
}

fn cohort(
    model: Arc<BurnupModel>,
    populations: [[f64; 2]; 2],
    transport: NuclideTransport,
) -> Arc<BurnupTransportState> {
    Arc::new(
        BurnupTransportState::new(
            populations
                .into_iter()
                .map(|population| {
                    BurnupInventory::new(Arc::clone(&model), population.to_vec()).unwrap()
                })
                .collect(),
            vec![0.01, 0.01],
            transport,
            transport_limits(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn two_cells(boundary_source: f64) -> DepositionSimulation<'static> {
    thermal(
        vec![
            cell(20, 1, 2.0, 1200.0, boundary_source),
            cell(10, 2, 3.0, 600.0, 0.0),
        ],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    )
}

fn cells() -> Vec<(CellId, MaterialId)> {
    // Deliberately opposite the thermal BTreeMap's order.
    vec![(CellId(20), MaterialId(1)), (CellId(10), MaterialId(2))]
}

fn neutron_flux(model: &BurnupModel, value: f64) -> GroupFlux {
    GroupFlux::new(Arc::clone(model.groups()), vec![value]).unwrap()
}

fn fluxes(model: &BurnupModel, value: f64) -> Vec<NeutronFluxSchedule> {
    (0..2)
        .map(|_| NeutronFluxSchedule::new(neutron_flux(model, value), vec![]).unwrap())
        .collect()
}

fn bind(
    thermal: &DepositionSimulation<'_>,
    state: Arc<BurnupTransportState>,
    value: f64,
) -> BurnupTransportHeating {
    let sources =
        BurnupTransportCellSource::bind(Arc::clone(&state), cells(), fluxes(state.model(), value))
            .unwrap();
    BurnupTransportHeating::new(thermal, sources, transport_budget()).unwrap()
}

fn phi(rate: f64, time: f64) -> f64 {
    -(-rate * time).exp_m1() / rate
}

#[test]
fn g1_migration_and_daughter_ingrowth_release_heat_at_the_actual_residence_location() {
    let initial = two_cells(3.0);
    let state = cohort(
        inventory_model(None, 900.0),
        [[100.0, 0.0], [0.0, 0.0]],
        flow(0.8),
    );
    let heating = bind(&initial, state, 0.0);
    let (heated, accepted, audit) = heating
        .advance_to(&initial, 1.0, controls(2.0), || Ok(()))
        .unwrap();
    let left = accepted.source(CellId(20)).unwrap();
    let right = accepted.source(CellId(10)).unwrap();
    assert!(Arc::ptr_eq(left.state(), right.state()));
    assert_eq!(left.compartment_index(), 0);
    assert_eq!(right.compartment_index(), 1);
    assert!(std::ptr::eq(
        accepted.last_step(CellId(20)).unwrap().transport(),
        accepted.last_step(CellId(10)).unwrap().transport(),
    ));
    let inventories = left.state().inventories();
    let survives = (-0.5_f64).exp();
    let stays = (-0.8_f64).exp();
    for (actual, expected) in [
        (inventories[0].populations()[0], 100.0 * survives * stays),
        (
            inventories[1].populations()[0],
            100.0 * survives * (1.0 - stays),
        ),
        (
            inventories[0].populations()[1],
            100.0 * (1.0 - survives) * stays,
        ),
        (
            inventories[1].populations()[1],
            100.0 * (1.0 - survives) * (1.0 - stays),
        ),
    ] {
        close(actual, expected, 1e-9);
    }
    let decays_left = 50.0 * phi(1.3, 1.0);
    let total_decays = 100.0 * (1.0 - survives);
    let decays_right = total_decays - decays_left;
    close(
        inventories[0].cumulative_reaction_counts()[0],
        decays_left,
        1e-9,
    );
    close(
        inventories[1].cumulative_reaction_counts()[0],
        decays_right,
        1e-9,
    );
    close(
        heated.temperature_k(CellId(20)).unwrap(),
        300.0 + decays_left + 3.0 / 8.0,
        1e-9,
    );
    close(
        heated.temperature_k(CellId(10)).unwrap(),
        300.0 + 8.0 * decays_right / 6.0,
        1e-9,
    );
    close(accepted.audit().thermal_credit_j, 8.0 * total_decays, 1e-8);
    close(audit.external_heat_j, 3.0 + 8.0 * total_decays, 1e-8);
    close(
        accepted.audit().released.radiation_j,
        2.0 * total_decays,
        1e-8,
    );
    close(accepted.audit().decays, total_decays, 1e-9);
    assert_eq!(accepted.audit().produced_nuclei, 0.0);
    assert_eq!(accepted.audit().substeps, 1);
    assert_eq!(accepted.time_s(), heated.time_s());
    assert_eq!(heated.cells()[&CellId(20)].mass_kg, 2.0);
    assert_eq!(heated.cells()[&CellId(10)].mass_kg, 3.0);
}

#[test]
fn g1_continuous_feed_moves_and_decays_with_its_actual_age_before_heating_each_cell() {
    let initial = two_cells(0.0);
    let state = cohort(inventory_model(None, 900.0), [[0.0, 0.0]; 2], flow(0.8));
    let movement = Arc::new(
        NuclideTransportSchedule::new(&state, flow(0.8), vec![], transport_limits(), 0).unwrap(),
    );
    let production = [40.0, 0.0]
        .into_iter()
        .map(|rate| {
            BurnupProductionSchedule::new(Arc::clone(state.model()), vec![rate, 0.0], vec![], 0)
                .unwrap()
        })
        .collect();
    let sources = BurnupTransportCellSource::bind_with_schedules(
        Arc::clone(&state),
        cells(),
        fluxes(state.model(), 0.0),
        movement,
        production,
    )
    .unwrap();
    let heating = BurnupTransportHeating::new(&initial, sources, transport_budget()).unwrap();
    let (heated, accepted, _) = heating
        .advance_to(&initial, 1.0, controls(2.0), || Ok(()))
        .unwrap();
    let state = accepted.source(CellId(20)).unwrap().state();
    let left_parent = 40.0 * phi(1.3, 1.0);
    let total_parent = 40.0 * phi(0.5, 1.0);
    let left_all = 40.0 * phi(0.8, 1.0);
    for (actual, expected) in [
        (state.inventories()[0].populations()[0], left_parent),
        (
            state.inventories()[1].populations()[0],
            total_parent - left_parent,
        ),
        (
            state.inventories()[0].populations()[1],
            left_all - left_parent,
        ),
        (
            state.inventories()[1].populations()[1],
            40.0 - left_all - total_parent + left_parent,
        ),
    ] {
        close(actual, expected, 1e-9);
    }
    let left_decays = 20.0 / 1.3 * (1.0 - phi(1.3, 1.0));
    let all_decays = 40.0 - total_parent;
    close(
        heated.temperature_k(CellId(20)).unwrap(),
        300.0 + left_decays,
        1e-9,
    );
    close(
        heated.temperature_k(CellId(10)).unwrap(),
        300.0 + 8.0 * (all_decays - left_decays) / 6.0,
        1e-9,
    );
    close(accepted.audit().decays, all_decays, 1e-9);
    close(accepted.audit().produced_nuclei, 40.0, 1e-9);
    assert_eq!(
        state.inventories()[0].cumulative_produced_nuclei(),
        [40.0, 0.0]
    );
    assert_eq!(
        state.inventories()[1].cumulative_produced_nuclei(),
        [0.0, 0.0]
    );
    close(
        state.cumulative_transfer_counts()[0].iter().sum(),
        40.0 - left_all,
        1e-9,
    );
    assert_eq!(initial.temperature_k(CellId(20)).unwrap(), 300.0);
}

#[test]
fn g3_independent_feed_flux_and_flow_shutdowns_preserve_absolute_time_on_rebinding() {
    let initial = two_cells(0.0);
    let state = cohort(
        inventory_model(Some((0.25, 0.25)), 900.0),
        [[50.0, 0.0], [0.0, 0.0]],
        flow(0.7),
    );
    let movement = Arc::new(
        NuclideTransportSchedule::new(
            &state,
            flow(0.7),
            vec![NuclideTransportChange {
                time_s: 0.75,
                transport: flow(0.0),
            }],
            transport_limits(),
            1,
        )
        .unwrap(),
    );
    let neutron_histories: Vec<_> = (0..2)
        .map(|_| {
            NeutronFluxSchedule::new(
                neutron_flux(state.model(), 1.0),
                vec![NeutronFluxChange {
                    time_s: 0.5,
                    flux: neutron_flux(state.model(), 0.0),
                }],
            )
            .unwrap()
        })
        .collect();
    let production: Vec<_> = [10.0, 0.0]
        .into_iter()
        .map(|rate| {
            BurnupProductionSchedule::new(
                Arc::clone(state.model()),
                vec![rate, 0.0],
                vec![BurnupProductionChange {
                    time_s: 0.25,
                    rates_nuclei_s: vec![0.0, 0.0],
                }],
                1,
            )
            .unwrap()
        })
        .collect();
    let make_sources = |state: Arc<BurnupTransportState>| {
        BurnupTransportCellSource::bind_with_schedules(
            state,
            cells(),
            neutron_histories.clone(),
            Arc::clone(&movement),
            production.clone(),
        )
        .unwrap()
    };
    let heating =
        BurnupTransportHeating::new(&initial, make_sources(state), transport_budget()).unwrap();
    let (heated, accepted, _) = heating
        .advance_to(&initial, 1.0, controls(2.0), || Ok(()))
        .unwrap();
    assert_eq!(accepted.audit().substeps, 4);
    let parent_at_feed_stop = 50.0 * (-0.75_f64 * 0.25).exp() + 10.0 * phi(0.75, 0.25);
    let parent_at_flux_stop = parent_at_feed_stop * (-0.75_f64 * 0.25).exp();
    let final_parent = parent_at_flux_stop * (-0.5_f64 * 0.5).exp();
    let captures = 0.25
        * (50.0 * phi(0.75, 0.25)
            + 10.0 / 0.75 * (0.25 - phi(0.75, 0.25))
            + parent_at_feed_stop * phi(0.75, 0.25));
    let final_state = accepted.source(CellId(20)).unwrap().state();
    close(
        final_state
            .inventories()
            .iter()
            .map(|inventory| inventory.populations()[0])
            .sum(),
        final_parent,
        1e-9,
    );
    close(accepted.audit().captures, captures, 1e-9);
    close(
        accepted.audit().decays,
        52.5 - final_parent - captures,
        1e-9,
    );
    close(accepted.audit().produced_nuclei, 2.5, 1e-9);
    let (middle_heat, middle, _) = heating
        .advance_to(&initial, 0.5, controls(2.0), || Ok(()))
        .unwrap();
    let rebound = BurnupTransportHeating::new(
        &middle_heat,
        make_sources(Arc::clone(middle.source(CellId(20)).unwrap().state())),
        transport_budget(),
    )
    .unwrap();
    let (restarted_heat, restarted, _) = rebound
        .advance_to(&middle_heat, 1.0, controls(2.0), || Ok(()))
        .unwrap();
    assert_eq!(restarted.audit().produced_nuclei, 0.0);
    assert_eq!(restarted.audit().captures, 0.0);
    assert!(restarted.audit().decays > 0.0);
    for cell in [CellId(20), CellId(10)] {
        close(
            restarted_heat.temperature_k(cell).unwrap(),
            heated.temperature_k(cell).unwrap(),
            1e-9,
        );
    }
    let restarted_state = restarted.source(CellId(20)).unwrap().state();
    for (a, b) in restarted_state
        .inventories()
        .iter()
        .zip(final_state.inventories())
    {
        for (&actual, &expected) in a.populations().iter().zip(b.populations()) {
            close(actual, expected, 1e-9);
        }
        assert_eq!(
            a.cumulative_produced_nuclei(),
            b.cumulative_produced_nuclei()
        );
    }
}

#[test]
fn g0_transport_binding_requires_complete_shared_history_and_matching_physical_cells() {
    let initial = two_cells(0.0);
    let state = cohort(
        inventory_model(None, 900.0),
        [[10.0, 0.0], [0.0, 0.0]],
        flow(0.8),
    );
    let make_sources = |state: Arc<BurnupTransportState>| {
        BurnupTransportCellSource::bind(Arc::clone(&state), cells(), fluxes(state.model(), 0.0))
            .unwrap()
    };
    let sources = make_sources(Arc::clone(&state));
    let complete =
        BurnupTransportHeating::new(&initial, sources.clone(), transport_budget()).unwrap();
    let (_, idle, _) = complete
        .advance_to(&initial, 0.0, controls(1.0), || Ok(()))
        .unwrap();
    assert_eq!(idle.audit().substeps, 0);
    assert!(Arc::ptr_eq(
        idle.source(CellId(20)).unwrap().state(),
        idle.source(CellId(10)).unwrap().state(),
    ));
    assert!(
        BurnupTransportHeating::new(&initial, vec![sources[0].clone()], transport_budget())
            .is_err()
    );
    let rebound_same_state = make_sources(Arc::clone(&state));
    assert!(
        BurnupTransportHeating::new(
            &initial,
            vec![sources[0].clone(), rebound_same_state[1].clone()],
            transport_budget(),
        )
        .is_err()
    );
    let independent = make_sources(Arc::new(state.as_ref().clone()));
    assert!(
        BurnupTransportHeating::new(
            &initial,
            vec![sources[0].clone(), independent[1].clone()],
            transport_budget()
        )
        .is_err()
    );
    let mut wrong_cell = cell(20, 1, 2.0, 1200.0, 0.0);
    wrong_cell.geometry.volume_m3 = 0.02;
    let wrong_volume = thermal(
        vec![wrong_cell, cell(10, 2, 3.0, 600.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    );
    assert!(matches!(
        BurnupTransportHeating::new(&wrong_volume, sources.clone(), transport_budget()),
        Err(FissionThermalError::BindingChanged(CellId(20)))
    ));
    let wrong_material = BurnupTransportCellSource::bind(
        Arc::clone(&state),
        vec![(CellId(20), MaterialId(2)), (CellId(10), MaterialId(2))],
        fluxes(state.model(), 0.0),
    )
    .unwrap();
    assert!(matches!(
        BurnupTransportHeating::new(&initial, wrong_material, transport_budget()),
        Err(FissionThermalError::MaterialMismatch { .. })
    ));
    assert!(
        BurnupTransportCellSource::bind(
            Arc::clone(&state),
            vec![(CellId(20), MaterialId(1)); 2],
            fluxes(state.model(), 0.0),
        )
        .is_err()
    );
    let other = inventory_model(None, 900.0);
    let production = (0..2)
        .map(|_| {
            BurnupProductionSchedule::new(Arc::clone(&other), vec![0.0; 2], vec![], 0).unwrap()
        })
        .collect();
    let movement = Arc::new(
        NuclideTransportSchedule::new(&state, flow(0.8), vec![], transport_limits(), 0).unwrap(),
    );
    assert!(
        BurnupTransportCellSource::bind_with_schedules(
            Arc::clone(&state),
            cells(),
            fluxes(state.model(), 0.0),
            movement,
            production,
        )
        .is_err()
    );
}

#[test]
fn g4_work_cancellation_and_material_refusals_leave_all_compartments_and_heat_unchanged() {
    let initial = two_cells(0.0);
    let state = cohort(
        inventory_model(None, 900.0),
        [[100.0, 0.0], [0.0, 0.0]],
        flow(0.8),
    );
    let sources =
        BurnupTransportCellSource::bind(Arc::clone(&state), cells(), fluxes(state.model(), 0.0))
            .unwrap();
    let heating =
        BurnupTransportHeating::new(&initial, sources.clone(), transport_budget()).unwrap();
    let mut strict = transport_budget();
    strict.limits.burnup.propagation.maximum_work = 0;
    let refused = BurnupTransportHeating::new(&initial, sources, strict).and_then(|source| {
        source
            .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
            .map(|_| ())
    });
    assert!(refused.is_err());
    let mut calls = 0;
    assert!(matches!(
        heating.advance_to(&initial, 1.0, controls(1.0), || {
            calls += 1;
            if calls == 20 {
                Err(DepositionError::MaterialDomain)
            } else {
                Ok(())
            }
        }),
        Err(FissionThermalError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(calls, 20);
    let limited = bind(
        &initial,
        cohort(
            inventory_model(None, 300.00001),
            [[100.0, 0.0], [0.0, 0.0]],
            flow(0.8),
        ),
        0.0,
    );
    assert!(
        limited
            .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
            .is_err()
    );
    for source in [&heating, &limited] {
        assert_eq!(source.time_s(), 0.0);
        assert_eq!(source.audit().thermal_credit_j, 0.0);
        let checkpoint = source.source(CellId(20)).unwrap().state();
        assert_eq!(checkpoint.inventories()[0].populations(), &[100.0, 0.0]);
        assert_eq!(checkpoint.inventories()[1].populations(), &[0.0, 0.0]);
        assert_eq!(
            checkpoint.cumulative_transfer_counts(),
            vec![vec![0.0, 0.0]]
        );
        assert_eq!(
            checkpoint.inventories()[0].cumulative_reaction_counts(),
            [0.0]
        );
    }
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.temperature_k(CellId(20)).unwrap(), 300.0);
    assert_eq!(initial.temperature_k(CellId(10)).unwrap(), 300.0);
}

#[test]
fn g1_adaptive_mobile_capture_heat_converges_to_the_independent_nonlinear_solution() {
    let initial = thermal(
        vec![cell(20, 1, 2.0, 1200.0, 0.0), cell(10, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let mut transport = flow(0.7);
    transport.transfers.push(NuclideTransfer {
        from_compartment: 1,
        to_compartment: 0,
        rates_per_s: vec![0.7, 0.7],
    });
    let state = cohort(
        inventory_model(Some((0.5, 2.5)), 900.0),
        [[80.0, 0.0]; 2],
        transport,
    );
    let sources = BurnupTransportCellSource::bind(
        Arc::clone(&state),
        vec![(CellId(20), MaterialId(1)), (CellId(10), MaterialId(1))],
        fluxes(state.model(), 1.0),
    )
    .unwrap();
    let heating = BurnupTransportHeating::new(&initial, sources, transport_budget()).unwrap();
    // Symmetry preserves equal inventories despite nonzero transfers. C=8 J/K,
    // event heat=8 J, so T=300+80-N. The total removal rate is 1+(80-N)/300:
    // N'=-a*N+b*N^2, a=1+80/300 and b=1/300, solved independently below.
    let a: f64 = 1.0 + 80.0 / 300.0;
    let survival = (-a).exp();
    let exact_parent = 80.0 * a * survival / (1.0 + 80.0 / 300.0 * survival);
    let integrated_parent = 300.0 * (a - exact_parent / 300.0).ln();
    let mut errors = Vec::new();
    let mut comparison_counts = Vec::new();
    for tolerance in [1e-2, 1e-5] {
        let settings = NuclearAdaptiveControls {
            step: StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-6,
            },
            initial_dt_s: 1.0,
            minimum_dt_s: 1e-9,
            maximum_trials: 4096,
            maximum_substeps: 4096,
            tolerances: NuclearTimeTolerances {
                nuclei: tolerance,
                nuclear_energy_j: 8.0 * tolerance,
                specific_energy_j_kg: 4.0 * tolerance,
                temperature_k: tolerance,
                relative: 0.0,
                ..NuclearTimeTolerances::default()
            },
        };
        let (heated, accepted, report) = heating
            .advance_adaptive_to(&initial, 1.0, settings, || Ok(()))
            .unwrap();
        let left = accepted.source(CellId(20)).unwrap();
        let right = accepted.source(CellId(10)).unwrap();
        assert!(Arc::ptr_eq(left.state(), right.state()));
        let actual_parent = left.state().inventories()[0].populations()[0];
        errors.push((actual_parent - exact_parent).abs());
        comparison_counts.push(report.accepted_comparisons);
        close(
            right.state().inventories()[1].populations()[0],
            actual_parent,
            1e-8,
        );
        close(
            heated.temperature_k(CellId(20)).unwrap(),
            380.0 - actual_parent,
            1e-8,
        );
        close(
            heated.temperature_k(CellId(10)).unwrap(),
            380.0 - actual_parent,
            1e-8,
        );
        close(
            accepted.audit().released.deposited_j,
            16.0 * (80.0 - actual_parent),
            1e-7,
        );
        assert!(report.accuracy_rejections > 0);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        assert_eq!(accepted.audit().substeps, report.thermal.substeps);
        assert_eq!(report.thermal.substeps, 2 * report.accepted_comparisons);
        if tolerance < 1e-3 {
            close(
                left.state().cumulative_transfer_counts()[0][0],
                0.7 * integrated_parent,
                0.02,
            );
            close(accepted.audit().decays, integrated_parent, 0.02);
        }
    }
    assert!(
        errors[1] < errors[0] / 5.0,
        "nonlinear population errors {errors:?}"
    );
    assert!(errors[1] < 0.02, "nonlinear population errors {errors:?}");
    assert!(comparison_counts[1] > comparison_counts[0]);
}
