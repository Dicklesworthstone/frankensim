use super::*;
use crate::diffusion::burnup::DiffusionBurnupSources;
use crate::diffusion::{DelayedGroup, PrecursorTransfer, PrecursorTransport};

fn empty_spatial(maximum_k: f64) -> Arc<DiffusionBurnupState> {
    let original = spatial(maximum_k);
    Arc::new(
        DiffusionBurnupState::new(
            DiffusionState::new(
                Arc::clone(original.neutrons().model()),
                vec![0.0, 0.0],
                vec![],
            )
            .unwrap(),
            original
                .inventories()
                .iter()
                .map(|inventory| {
                    BurnupInventory::new(Arc::clone(inventory.model()), vec![0.0, 0.0]).unwrap()
                })
                .collect(),
            original.bindings().to_vec(),
        )
        .unwrap(),
    )
}

fn external(neutrons: f64, nuclei: f64) -> DiffusionBurnupSources {
    DiffusionBurnupSources {
        neutron_rates_s: vec![neutrons, 0.0],
        precursor_rates_s: vec![],
        nuclide_rates_s: vec![vec![0.0, 0.0], vec![nuclei, 0.0]],
    }
}

fn fed_source(state: Arc<DiffusionBurnupState>) -> Vec<DiffusionCellSource> {
    let schedule = Arc::new(
        DiffusionSourceSchedule::new(
            &state,
            external(10.0, 100.0),
            vec![
                DiffusionSourceChange {
                    time_s: 0.3,
                    sources: external(0.0, 100.0),
                },
                DiffusionSourceChange {
                    time_s: 0.55,
                    sources: external(0.0, 0.0),
                },
            ],
            2,
        )
        .unwrap(),
    );
    DiffusionCellSource::bind_with_sources(
        state,
        vec![(CellId(20), MaterialId(1)), (CellId(10), MaterialId(2))],
        schedule,
    )
    .unwrap()
}

#[test]
fn g1_scheduled_spatial_sources_create_particles_targets_and_material_heat_on_one_clock() {
    let initial = two_cells(0.0);
    let source =
        DiffusionHeating::new(&initial, fed_source(empty_spatial(900.0)), spatial_budget())
            .unwrap();
    // The requested native interval spans both independent on/off boundaries.
    let (heated, accepted, _) = source
        .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    let left = accepted.source(CellId(20)).unwrap();
    let right = accepted.source(CellId(10)).unwrap();
    assert!(Arc::ptr_eq(left.state(), right.state()));
    assert!(Arc::ptr_eq(
        left.source_schedule().unwrap(),
        right.source_schedule().unwrap()
    ));
    assert_eq!(accepted.audit().substeps, 3);
    let state = right.state();
    let captures = state.inventories()[1].cumulative_reaction_counts()[0];
    assert!(captures > 0.0 && captures < 3.0);
    close(
        state.inventories()[1].populations()[0],
        55.0 - captures,
        1e-8,
    );
    close(state.inventories()[1].populations()[1], captures, 1e-8);
    close(
        state.neutrons().neutron_counts().iter().sum::<f64>(),
        3.0 - captures,
        1e-8,
    );
    close(accepted.audit().produced_nuclei, 55.0, 1e-9);
    close(accepted.audit().captures, captures, 1e-8);
    close(accepted.audit().thermal_credit_j, 8.0 * captures, 1e-8);
    close(accepted.audit().released.radiation_j, 2.0 * captures, 1e-8);
    close(accepted.audit().released.escaped_j, captures, 1e-8);
    close(
        heated.temperature_k(CellId(10)).unwrap(),
        300.0 + 8.0 * captures / 6.0,
        1e-8,
    );
    assert_eq!(heated.temperature_k(CellId(20)).unwrap(), 300.0);
    assert_eq!(heated.cells()[&CellId(10)].mass_kg, 3.0);
    assert_eq!(accepted.time_s(), heated.time_s());
    assert_eq!(source.time_s(), 0.0);
    assert_eq!(
        source.source(CellId(10)).unwrap().state().inventories()[1].populations(),
        &[0.0, 0.0]
    );
    // At a boundary the replacement vector is active; no old-rate replay.
    assert_eq!(
        right
            .source_schedule()
            .unwrap()
            .sources_at(0.3)
            .unwrap()
            .neutron_rates_s,
        [0.0, 0.0]
    );
    assert_eq!(right.source_schedule().unwrap().next_after(0.55), None);
}

#[test]
fn g1_incoming_precursors_move_before_emission_and_drive_downstream_capture_heat() {
    let original = spatial(900.0);
    let mut data = original.neutrons().model().data().clone();
    data.delayed_groups = vec![DelayedGroup {
        data_id: "manufactured-mobile-precursor".into(),
        decay_constant_s_inv: 1.0,
        spectrum: vec![1.0],
    }];
    // No neutron face: any downstream neutron and heat must come from a
    // precursor transported there before it emits, not neutron diffusion.
    data.links.clear();
    let transport = PrecursorTransport {
        transfers: vec![PrecursorTransfer {
            from_cell: 0,
            to_cell: 1,
            volumetric_flow_m3_s: 0.01,
        }],
        outflows: vec![],
    };
    let model = Arc::new(
        DiffusionModel::new_with_transport(
            data,
            transport.clone(),
            spatial_budget().diffusion.limits,
        )
        .unwrap(),
    );
    let state = Arc::new(
        DiffusionBurnupState::new(
            DiffusionState::new(model, vec![0.0, 0.0], vec![0.0, 0.0]).unwrap(),
            original.inventories().to_vec(),
            original.bindings().to_vec(),
        )
        .unwrap(),
    );
    let forcing = DiffusionBurnupSources {
        neutron_rates_s: vec![0.0, 0.0],
        precursor_rates_s: vec![4.0, 0.0],
        nuclide_rates_s: vec![vec![0.0, 0.0], vec![0.0, 0.0]],
    };
    let schedule = Arc::new(DiffusionSourceSchedule::new(&state, forcing, vec![], 0).unwrap());
    let sources = DiffusionCellSource::bind_with_sources(
        state,
        vec![(CellId(20), MaterialId(1)), (CellId(10), MaterialId(2))],
        schedule,
    )
    .unwrap();
    let initial = two_cells(0.0);
    let source = DiffusionHeating::new(&initial, sources, spatial_budget()).unwrap();
    let (heated, accepted, _) = source
        .advance_to(&initial, 1.0, controls(0.1), || Ok(()))
        .unwrap();
    let final_state = accepted.source(CellId(10)).unwrap().state();
    let captures = 100.0 - final_state.inventories()[1].populations()[0];
    assert!(captures > 0.05);
    assert!(final_state.neutrons().neutron_counts()[1] > 0.0);
    assert_eq!(
        final_state.neutrons().model().precursor_transport(),
        &transport
    );
    let remaining = final_state
        .neutrons()
        .neutron_counts()
        .iter()
        .chain(final_state.neutrons().precursor_counts())
        .sum::<f64>();
    close(remaining + captures, 4.0, 1e-8);
    // Independently integrated dC0/dt=4-2 C0, dC1/dt=C0-C1.
    let e1 = fs_math::det::exp(-1.0);
    let e2 = fs_math::det::exp(-2.0);
    close(
        final_state.neutrons().precursor_counts()[0],
        2.0 * (1.0 - e2),
        1e-9,
    );
    close(
        final_state.neutrons().precursor_counts()[1],
        2.0 * (1.0 - 2.0 * e1 + e2),
        1e-9,
    );
    close(accepted.audit().thermal_credit_j, 8.0 * captures, 1e-8);
    close(
        heated.temperature_k(CellId(10)).unwrap(),
        300.0 + 8.0 * captures / 6.0,
        1e-8,
    );
    assert_eq!(accepted.audit().produced_nuclei, 0.0);
}

#[test]
fn g4_source_cancellation_and_material_refusal_publish_no_particle_or_heat_history() {
    let initial = two_cells(0.0);
    let source = DiffusionHeating::new(
        &initial,
        fed_source(empty_spatial(300.000001)),
        spatial_budget(),
    )
    .unwrap();
    assert!(
        source
            .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        source.advance_to(&initial, 1.0, controls(0.05), || {
            polls += 1;
            if polls == 20 {
                Err(DepositionError::MaterialDomain)
            } else {
                Ok(())
            }
        }),
        Err(FissionThermalError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(source.time_s(), 0.0);
    assert_eq!(source.audit().produced_nuclei, 0.0);
    assert_eq!(source.audit().thermal_credit_j, 0.0);
    assert_eq!(
        source.source(CellId(10)).unwrap().state().inventories()[1].populations(),
        &[0.0, 0.0]
    );
    assert_eq!(initial.temperature_k(CellId(10)).unwrap(), 300.0);
}

#[test]
fn g0_spatial_schedule_admits_all_rates_and_rejects_mixed_source_histories() {
    let state = empty_spatial(900.0);
    let mut malformed = external(1.0, 1.0);
    malformed.nuclide_rates_s[1][0] = f64::NAN;
    assert!(
        DiffusionSourceSchedule::new(
            &state,
            external(0.0, 0.0),
            vec![DiffusionSourceChange {
                time_s: 2.0,
                sources: malformed
            }],
            1
        )
        .is_err()
    );
    assert!(
        DiffusionSourceSchedule::new(
            &state,
            external(0.0, 0.0),
            vec![DiffusionSourceChange {
                time_s: 0.0,
                sources: external(1.0, 1.0)
            }],
            1
        )
        .is_err()
    );
    assert!(
        DiffusionSourceSchedule::new(
            &state,
            external(0.0, 0.0),
            vec![DiffusionSourceChange {
                time_s: 2.0,
                sources: external(1.0, 1.0)
            }],
            0
        )
        .is_err()
    );
    assert!(matches!(
        DiffusionSourceSchedule::new_with_checkpoint(
            &state,
            external(0.0, 0.0),
            vec![],
            0,
            || Err(FissionError::Cancelled)
        ),
        Err(DiffusionBurnupError::Fission(FissionError::Cancelled))
    ));
    let schedule =
        Arc::new(DiffusionSourceSchedule::new(&state, external(1.0, 2.0), vec![], 0).unwrap());
    assert!(schedule.sources_at(-1.0).is_err());
    assert!(
        DiffusionCellSource::bind_with_sources(
            empty_spatial(900.0),
            vec![(CellId(20), MaterialId(1)), (CellId(10), MaterialId(2))],
            Arc::clone(&schedule)
        )
        .is_err()
    );
    let first = fed_source(Arc::clone(&state));
    let second = DiffusionCellSource::bind_with_sources(
        state,
        vec![(CellId(20), MaterialId(1)), (CellId(10), MaterialId(2))],
        schedule,
    )
    .unwrap();
    let initial = two_cells(0.0);
    let mixed = DiffusionHeating::new(
        &initial,
        vec![first[0].clone(), second[1].clone()],
        spatial_budget(),
    )
    .unwrap();
    assert!(matches!(
        mixed.advance_to(&initial, 0.1, controls(0.1), || Ok(())),
        Err(FissionThermalError::Invalid(_))
    ));
}
