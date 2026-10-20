//! G1/G3: actual shared gas transport enters finite thermal/electrical/motion
//! owners. Source data is manufactured; this is not engine validation.
use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::network::snapshot::{self as shared, SnapshotNozzleError, SnapshotNozzleLink},
    snapshot::{GasSnapshotError, GasSnapshotSolve},
};

const DT: f64 = 0.00005;

fn snapshot_controls() -> GasSnapshotSolve {
    GasSnapshotSolve {
        maximum_withdrawn_fraction: 0.05,
        mass_tolerance_kg: 1e-13,
        energy_tolerance_j: 1e-7,
    }
}
fn coupled_controls() -> NetworkSolve {
    NetworkSolve {
        maximum_chambers: 4,
        maximum_boundary_solves: 1200,
        ..controls()
    }
}
fn manifold(molar_mass: f64) -> BTreeMap<u64, GasPiston> {
    let mut map = gases(molar_mass);
    map.insert(30, source_gas("C", 0.044, 60000.0));
    map.insert(40, source_gas("D", 0.020, 150000.0));
    map
}
fn restrictions() -> Vec<SnapshotNozzleLink> {
    [(901, [20, 10]), (502, [40, 10]), (17, [10, 30])]
        .into_iter()
        .map(|(id, parcels)| SnapshotNozzleLink {
            id,
            parcels,
            nozzle: GasNozzle {
                area_m2: 2e-6,
                discharge_coefficient: 0.8,
            },
            nozzle_solve: chemical().gas,
        })
        .collect()
}
fn gas_storage(state: &MaterialExchangeState<'_>) -> f64 {
    state
        .chambers()
        .values()
        .map(|g| {
            g.gas().mass_kg() * g.gas().state().specific_internal_energy_j_kg()
                + 0.5 * g.mass_kg() * g.velocity_m_s().powi(2)
        })
        .sum()
}
fn total_storage(state: &MaterialExchangeState<'_>) -> f64 {
    // Independently reconstruct the fixture's sole capacitor, not its report.
    state.thermal().energy_j()
        + state.circuit().charges_c()[0].powi(2) / (2.0 * 0.01)
        + gas_storage(state)
}
fn has_species(gas: &GasPiston, name: &str) -> bool {
    let id = SpeciesId::new(name).unwrap();
    gas.gas()
        .model()
        .thermodynamics()
        .components()
        .iter()
        .any(|c| c.species() == &id)
}

#[test]
fn g3_manifold_ids_do_not_change_mixing_or_the_coupled_physical_endpoint() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(thermal, circuit, manifold(0.004));
            let drives = drives(state.chambers());
            let links = restrictions();
            let (a, report) = state
                .snapshot_flow_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &links,
                    snapshot_controls(),
                    coupled_controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let mut relabelled = links.clone();
            relabelled.reverse();
            for (i, link) in relabelled.iter_mut().enumerate() {
                link.id = 100 + i as u64;
            }
            let (b, other) = state
                .snapshot_flow_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &relabelled,
                    snapshot_controls(),
                    coupled_controls(),
                    &[0.0; 2],
                )
                .unwrap();
            assert_eq!(a.chambers(), b.chambers());
            assert_eq!(a.thermal().cells(), b.thermal().cells());
            assert_eq!(a.circuit(), b.circuit());
            assert_eq!(report.flow.exchange, other.flow.exchange);
            assert!(has_species(&a.chambers()[&10], "B"));
            assert!(has_species(&a.chambers()[&10], "D"));
            assert!(has_species(&a.chambers()[&30], "A"));
            // Newly arrived B/D cannot leave the mixing node until the next epoch.
            assert!(!has_species(&a.chambers()[&30], "B"));
            assert!(!has_species(&a.chambers()[&30], "D"));
            let work: f64 = report
                .exchange
                .report
                .source
                .chambers
                .iter()
                .map(|c| c.boundary.piston.exterior_pressure_work_j)
                .sum();
            assert!((total_storage(&a) + work - total_storage(&state)).abs() < 4e-6);
            assert!(report.total_mass_residual_kg.abs() <= snapshot_controls().mass_tolerance_kg);
            assert!(
                report.total_energy_residual_j.abs() <= coupled_controls().total_energy_tolerance_j
            );
            for gas in a.chambers().values() {
                assert_eq!(gas.last_window_duration_s(), Some(DT));
                assert!(matches!(
                    gas.audit_last_window().unwrap(),
                    ModuleWindowAudit::Unavailable { .. }
                ));
            }
            // Zero connections retain the original coupled solver exactly.
            let (plain, _) = state
                .step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    coupled_controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let (closed, _) = state
                .snapshot_flow_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &[],
                    snapshot_controls(),
                    coupled_controls(),
                    &[0.0; 2],
                )
                .unwrap();
            assert_eq!(plain.chambers(), closed.chambers());
            assert_eq!(plain.thermal().cells(), closed.thermal().cells());
            assert_eq!(plain.circuit(), closed.circuit());
            assert_eq!(state.time_s(), 0.0);
        })
    });
}

#[test]
fn g1_g3_species_data_reaches_wall_heat_charge_motion_and_sound_speed_over_time() {
    let run = |molar_mass| {
        fixture([0.8, 0.8], 150.0, |assembly, thermal, circuit| {
            with_cx(|cx, _| {
                let prepared = assembly.prepare(cx, &thermal).unwrap();
                let mut state =
                    MaterialExchangeState::from_parts(thermal, circuit, manifold(molar_mass));
                let drives = drives(state.chambers());
                let initial_energy = total_storage(&state);
                let initial_mass: f64 = state.chambers().values().map(|g| g.gas().mass_kg()).sum();
                let mut work = 0.0;
                for _ in 0..24 {
                    let target = state.time_s() + DT;
                    let report = state
                        .advance_with_snapshot_flow(
                            &prepared,
                            cx,
                            target,
                            &[],
                            &drives,
                            &[],
                            &gas_fixture::patches(),
                            &restrictions(),
                            snapshot_controls(),
                            coupled_controls(),
                            &[0.0; 2],
                        )
                        .unwrap();
                    work += report
                        .exchange
                        .report
                        .source
                        .chambers
                        .iter()
                        .map(|c| c.boundary.piston.exterior_pressure_work_j)
                        .sum::<f64>();
                    assert!((total_storage(&state) + work - initial_energy).abs() < 4e-5);
                    assert!(report.flow.exchange.maximum_species_mass_residual_kg < 1e-12);
                    assert!(report.exchange.report.electrical.contacts[0].heat_j > 0.0);
                }
                let mass: f64 = state.chambers().values().map(|g| g.gas().mass_kg()).sum();
                assert!((mass - initial_mass).abs() < 1e-12);
                assert!(has_species(&state.chambers()[&30], "B"));
                (
                    state.thermal().temperature_k(CellId(30)).unwrap(),
                    state.circuit().charges_c()[0],
                    state.chambers()[&10].position_m(),
                    state.chambers()[&10].gas().state().state().sound_speed,
                )
            })
        })
    };
    let light = run(0.004);
    let heavy = run(0.044);
    assert!((light.0 - heavy.0).abs() > 1e-8);
    assert!((light.1 - heavy.1).abs() > 1e-12);
    assert!((light.2 - heavy.2).abs() > 1e-10);
    assert!((light.3 - heavy.3).abs() > 1e-3);
}

#[test]
fn g3_reaction_consumes_the_joint_transferred_inventory_once() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let law = reaction(50.0);
            let product = law
                .catalog()
                .with_composition(
                    Composition::new(CompositionBasis::MoleFraction, vec![0.0, 1.0]).unwrap(),
                )
                .unwrap();
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                BTreeMap::from([(10, piston(&law)), (20, parcel_piston(product, 180000.0))]),
            );
            let drives = drives(state.chambers());
            let links = [restrictions()[0]];
            let parcels = state
                .chambers()
                .iter()
                .map(|(&id, g)| (id, g.gas().clone()))
                .collect();
            let endpoints = drives
                .iter()
                .map(|d| (d.chamber, d.piston_solve.gas))
                .collect();
            let transferred = shared::step_candidate(
                &parcels,
                &links,
                DT,
                &endpoints,
                snapshot_controls(),
                || false,
            )
            .unwrap();
            let expected = law
                .parcel_candidate(&transferred.parcels()[&10], DT, chemical(), || false)
                .unwrap();
            let (_, report) = state
                .snapshot_flow_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[ChamberReaction {
                        chamber: 10,
                        reaction: &law,
                        controls: chemical(),
                    }],
                    &gas_fixture::patches(),
                    &links,
                    snapshot_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            assert_eq!(report.exchange.report.source.reactions.len(), 1);
            assert_eq!(
                report.exchange.report.source.reactions[0].reaction,
                *expected.step()
            );
            assert!(report.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j);
            assert_eq!(state.time_s(), 0.0);
        })
    });
}

#[test]
fn g0_g4_failed_transfer_or_coupled_solve_keeps_all_live_histories_unchanged() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, gate| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut state = MaterialExchangeState::from_parts(thermal, circuit, manifold(0.004));
            let drives = drives(state.chambers());
            let original_cells = state.thermal().cells().clone();
            let original_circuit = state.circuit().clone();
            let original_gases = state.chambers().clone();
            let mut withdrawal = snapshot_controls();
            withdrawal.maximum_withdrawn_fraction = 1e-12;
            assert!(matches!(
                state.advance_with_snapshot_flow(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &restrictions(),
                    withdrawal,
                    coupled_controls(),
                    &[0.0; 2]
                ),
                Err(MaterialFlowError::SnapshotFlow(
                    SnapshotNozzleError::Transfer(GasSnapshotError::Withdrawal { .. })
                ))
            ));
            let mut limited = coupled_controls();
            limited.network.maximum_circuit_trials = 1;
            assert!(
                state
                    .advance_with_snapshot_flow(
                        &prepared,
                        cx,
                        DT,
                        &[],
                        &drives,
                        &[],
                        &gas_fixture::patches(),
                        &restrictions(),
                        snapshot_controls(),
                        limited,
                        &[0.0; 2]
                    )
                    .is_err()
            );
            assert!(
                state
                    .advance_with_snapshot_flow(
                        &prepared,
                        cx,
                        DT,
                        &[],
                        &[drives[0]; 4],
                        &[],
                        &gas_fixture::patches(),
                        &restrictions(),
                        snapshot_controls(),
                        coupled_controls(),
                        &[0.0; 2]
                    )
                    .is_err()
            );
            assert_eq!(state.thermal().cells(), &original_cells);
            assert_eq!(state.circuit(), &original_circuit);
            assert_eq!(state.chambers(), &original_gases);
            assert_eq!(state.time_s(), 0.0);
            let (preview, _) = state
                .snapshot_flow_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &restrictions(),
                    snapshot_controls(),
                    coupled_controls(),
                    &[0.0; 2],
                )
                .unwrap();
            state
                .advance_with_snapshot_flow(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &restrictions(),
                    snapshot_controls(),
                    coupled_controls(),
                    &[0.0; 2],
                )
                .unwrap();
            assert_eq!(state.chambers(), preview.chambers());
            assert_eq!(state.circuit(), preview.circuit());
            assert_eq!(state.thermal().cells(), preview.thermal().cells());
            gate.request();
            assert!(matches!(
                state.advance_with_snapshot_flow(
                    &prepared,
                    cx,
                    2.0 * DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &restrictions(),
                    snapshot_controls(),
                    coupled_controls(),
                    &[0.0; 2]
                ),
                Err(MaterialFlowError::Cancelled)
            ));
            assert_eq!(state.chambers(), preview.chambers());
            assert_eq!(state.thermal().cells(), preview.thermal().cells());
            assert_eq!(state.circuit(), preview.circuit());
            assert_eq!(state.time_s(), DT);
        })
    });
}
