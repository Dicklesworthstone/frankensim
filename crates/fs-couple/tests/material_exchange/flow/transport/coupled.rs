//! Laminate transport through the original material, heat and motion owners.
#[path = "coupled_gas_law.rs"]
mod gas_law;

use super::*;
use fs_couple::cavity_phs::piston::thermal::ConductanceSample;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::state::flow::transport::coupled::*;
use fs_material::porous::permeability::{
    DarcyGasBarrier, IntrinsicPermeability, PorousGasLayer, PorousGasSupport,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::{
    diffusion::dusty_gas::{DustyGasFaceSolve, DustyGasPores, layered::LayeredDustyGasMedium},
    network::transport::{
        self as material, GasTransportPaths,
        dusty_gas::LayeredDustyGasNetworkLink,
    },
};

fn laminate(length_scale: f64) -> LayeredDustyGasMedium {
    let layers = [(0.00003, 2e-14), (0.00007, 4e-14)]
        .into_iter()
        .map(|(length, permeability)| {
            PorousGasLayer::try_new(
                IntrinsicPermeability::isotropic(permeability).unwrap(),
                length * length_scale,
                0.0,
                "synthetic laminate intrinsic permeability",
            )
            .unwrap()
        })
        .collect();
    let barrier = DarcyGasBarrier::try_new(
        layers,
        [1.0, 0.0, 0.0],
        PorousGasSupport {
            temperature_k: (200.0, 1800.0),
            pressure_pa: (1000.0, 1e7),
            maximum_superficial_speed_m_s: 100.0,
        },
    )
    .unwrap();
    LayeredDustyGasMedium::try_new(
        barrier,
        vec![
            DustyGasPores::try_new(1e-7, 0.4, 2.0, "synthetic first pores").unwrap(),
            DustyGasPores::try_new(2e-7, 0.3, 1.5, "synthetic second pores").unwrap(),
        ],
    )
    .unwrap()
}

fn laminate_link<'a>(
    medium: &'a LayeredDustyGasMedium,
    pairs: &'a [BinaryDiffusionMaterial<'a>],
    context: &'a QueryPoint,
) -> LayeredDustyGasNetworkLink<'a> {
    LayeredDustyGasNetworkLink {
        id: 141,
        parcels: [10, 20],
        medium,
        area_m2: 0.1,
        pairs,
        gas_context: context,
        solve: DustyGasFaceSolve {
            maximum_temperature_difference_k: 25.0,
            maximum_relative_pressure_difference: 0.25,
            gradient_tolerance_per_m: 1e-8,
        },
    }
}

fn loose_time() -> CoupledTransportControls {
    CoupledTransportControls {
        absolute: CoupledTransportScales {
            species_mass_kg: 1.0,
            gas_temperature_k: 1000.0,
            gas_pressure_pa: 1e6,
            gas_volume_m3: 1.0,
            piston_velocity_m_s: 100.0,
            cell_energy_j: 1e6,
            charge_c: 1.0,
            flux_wb: 1.0,
        },
        relative_tolerance: 0.0,
        maximum_step_s: DT,
        minimum_step_s: DT / 4096.0,
        maximum_evaluations: 1000,
    }
}

#[test]
fn mixed_laminate_adaptation_accepts_exactly_the_two_original_fine_assemblies() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let card = pair_card(2e-5);
            let point = QueryPoint::new();
            let pairs = binding(&card, &point);
            let medium = laminate(1.0);
            let layers = [laminate_link(&medium, &pairs, &point)];
            let nozzles = [restriction()];
            let paths = GasTransportPaths {
                nozzles: &nozzles,
                layered_dusty_gas: &layers,
                ..Default::default()
            };
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(thermal, circuit, chambers());
            let mut drives = drives(state.chambers());
            drives[0].heat_j = 0.01;
            drives[1].heat_j = 0.02;
            let half = drives
                .iter()
                .map(|d| ChamberDrive {
                    heat_j: d.heat_j * 0.5,
                    ..*d
                })
                .collect::<Vec<_>>();
            let (middle, first) = state
                .transport_paths_step_candidate(
                    &prepared,
                    cx,
                    DT / 2.0,
                    &[],
                    &half,
                    &[],
                    &gas_fixture::patches(),
                    paths,
                    admission(),
                    coupled(),
                    &[0.0; 2],
                )
                .unwrap();
            let (expected, second) = middle
                .transport_paths_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &half,
                    &[],
                    &gas_fixture::patches(),
                    paths,
                    admission(),
                    coupled(),
                    &[0.0; 2],
                )
                .unwrap();
            let (actual, report) = state
                .coupled_transport_paths_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    paths,
                    admission(),
                    coupled(),
                    &[0.0; 2],
                    CoupledTransportControls {
                        maximum_evaluations: 3,
                        ..loose_time()
                    },
                )
                .unwrap();
            assert_eq!(actual.chambers(), expected.chambers());
            assert_eq!(actual.circuit(), expected.circuit());
            assert_eq!(actual.thermal().cells(), expected.thermal().cells());
            assert_eq!(report.steps[0].advance.flow, first.flow);
            assert_eq!(report.steps[1].advance.flow, second.flow);
            assert!(!report.steps[0].advance.flow.layered_dusty_gas.is_empty());
            assert!(!report.steps[0].advance.flow.nozzles.is_empty());
            assert_eq!(report.evaluations, 3);
            assert_eq!(report.accepted_intervals, 1);
            assert_eq!(report.rejected_intervals, 0);
            assert!(report.estimated_error.maximum() <= 1.0);
            assert!((report.balance.prescribed_gas_heat_j - 0.03).abs() < 1e-16);
            assert!(
                report.balance.total_mass_residual_kg.abs()
                    <= admission().exchange.mass_tolerance_kg
            );
            assert!(
                report.balance.total_energy_residual_j.abs() <= coupled().total_energy_tolerance_j
            );
            for gas in actual.chambers().values() {
                assert_eq!(gas.last_window_duration_s(), Some(DT));
                assert!(matches!(
                    gas.audit_last_window().unwrap(),
                    ModuleWindowAudit::Unavailable { .. }
                ));
            }
            assert_eq!(state.time_s(), 0.0);
        })
    });
}

#[test]
fn mixed_path_withdrawal_refines_without_relaxing_inventory_or_publishing_failed_work() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, gate| {
            let card = pair_card(2e-5);
            let point = QueryPoint::new();
            let pairs = binding(&card, &point);
            let medium = laminate(1.0);
            let layers = [laminate_link(&medium, &pairs, &point)];
            let nozzles = [restriction()];
            let paths = GasTransportPaths {
                nozzles: &nozzles,
                layered_dusty_gas: &layers,
                ..Default::default()
            };
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut state = MaterialExchangeState::from_parts(thermal, circuit, chambers());
            let drives = drives(state.chambers());
            let parcels = state
                .chambers()
                .iter()
                .map(|(&id, g)| (id, g.gas().clone()))
                .collect();
            let endpoints = drives
                .iter()
                .map(|d| (d.chamber, d.piston_solve.gas))
                .collect();
            let native = material::step_candidate_with_paths(
                &parcels,
                paths,
                DT,
                &endpoints,
                admission(),
                || false,
            )
            .unwrap()
            .into_parts()
            .1;
            let mut gross = BTreeMap::<u64, f64>::new();
            for flux in &native.exchange.fluxes {
                *gross.entry(flux.transfer.donor).or_default() += flux.transfer.mass_kg;
            }
            // Every initial chamber is pure, so this is the ORIGINAL species
            // allowance. Later arriving species cannot replenish it in a step.
            let maximum = gross
                .iter()
                .map(|(id, mass)| mass / state.chambers()[id].gas().mass_kg())
                .fold(0.0_f64, f64::max);
            assert!(maximum > 0.0);
            let physical = GasTransportSolve {
                exchange: SpeciesNetworkSolve {
                    maximum_species_fraction: maximum / 4.0,
                    ..admission().exchange
                },
                ..admission()
            };
            let cells = state.thermal().cells().clone();
            let circuit = state.circuit().clone();
            let gases = state.chambers().clone();
            let refused = state.advance_with_coupled_transport_paths(
                &prepared,
                cx,
                DT,
                &[],
                &drives,
                &[],
                &gas_fixture::patches(),
                paths,
                physical,
                coupled(),
                &[0.0; 2],
                CoupledTransportControls {
                    maximum_evaluations: 3,
                    ..loose_time()
                },
            );
            assert!(matches!(
                refused,
                Err(CoupledTransportError::WorkLimit { .. })
            ));
            assert_eq!(state.chambers(), &gases);
            assert_eq!(state.circuit(), &circuit);
            assert_eq!(state.thermal().cells(), &cells);
            let (_, report) = state
                .coupled_transport_paths_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    paths,
                    physical,
                    coupled(),
                    &[0.0; 2],
                    loose_time(),
                )
                .unwrap();
            assert!(report.rejected_intervals > 0);
            assert!(report.accepted_intervals > 1);
            assert!(report.evaluations > 3 * report.accepted_intervals);
            assert_eq!(report.steps.len(), 2 * report.accepted_intervals);
            assert!(
                report.balance.maximum_species_mass_residual_kg
                    <= physical.exchange.mass_tolerance_kg
            );
            assert_eq!(report.steps.first().unwrap().start_s, 0.0);
            assert_eq!(report.steps.last().unwrap().end_s, DT);
            for adjacent in report.steps.windows(2) {
                assert_eq!(adjacent[0].end_s, adjacent[1].start_s);
            }
            gate.request();
            let cancelled = state.advance_with_coupled_transport_paths(
                &prepared,
                cx,
                DT,
                &[],
                &drives,
                &[],
                &gas_fixture::patches(),
                paths,
                physical,
                coupled(),
                &[0.0; 2],
                loose_time(),
            );
            assert!(matches!(cancelled, Err(CoupledTransportError::Cancelled)));
            assert_eq!(state.chambers(), &gases);
            assert_eq!(state.circuit(), &circuit);
            assert_eq!(state.thermal().cells(), &cells);
            assert_eq!(state.time_s(), 0.0);
        })
    });
}

#[test]
fn laminate_and_other_paths_feed_the_original_complete_assembly_once() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let d = pair_card(2e-5);
            let q = QueryPoint::new();
            let pairs = binding(&d, &q);
            let medium = laminate(1.0);
            let layers = [laminate_link(&medium, &pairs, &q)];
            let nozzles = [restriction()];
            let paths = GasTransportPaths {
                nozzles: &nozzles,
                layered_dusty_gas: &layers,
                ..Default::default()
            };
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(thermal, circuit, chambers());
            let drives = drives(state.chambers());
            let parcels = state
                .chambers()
                .iter()
                .map(|(&id, g)| (id, g.gas().clone()))
                .collect();
            let endpoints = drives
                .iter()
                .map(|d| (d.chamber, d.piston_solve.gas))
                .collect();
            let (transported, expected_transport) = material::step_candidate_with_paths(
                &parcels,
                paths,
                DT,
                &endpoints,
                admission(),
                || false,
            )
            .unwrap()
            .into_parts();
            let staged = MaterialExchangeState::from_parts(
                state.thermal().clone(),
                state.circuit().clone(),
                transported
                    .into_iter()
                    .map(|(id, gas)| (id, GasPiston::new(gas, 1.0, 0.01, 0.0).unwrap()))
                    .collect(),
            );
            let (expected, _) = staged
                .step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    coupled(),
                    &[0.0; 2],
                )
                .unwrap();
            let (actual, report) = state
                .transport_paths_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    paths,
                    admission(),
                    coupled(),
                    &[0.0; 2],
                )
                .unwrap();
            assert_eq!(report.flow, expected_transport);
            assert_eq!(actual.thermal().cells(), expected.thermal().cells());
            assert_eq!(actual.circuit(), expected.circuit());
            for (&id, gas) in actual.chambers() {
                let wanted = &expected.chambers()[&id];
                assert_eq!(gas.gas(), wanted.gas());
                assert_eq!(gas.position_m(), wanted.position_m());
                assert_eq!(gas.velocity_m_s(), wanted.velocity_m_s());
                assert!(matches!(
                    gas.audit_last_window().unwrap(),
                    ModuleWindowAudit::Unavailable { .. }
                ));
            }
            let work: f64 = report
                .exchange
                .report
                .source
                .chambers
                .iter()
                .map(|c| c.boundary.piston.exterior_pressure_work_j)
                .sum();
            assert!((storage(&actual) + work - storage(&state)).abs() < 4e-6);
            assert!(report.total_mass_residual_kg.abs() <= admission().exchange.mass_tolerance_kg);

            let patches = gas_fixture::patches();
            let surfaces = patches
                .iter()
                .map(|p| GasSurface {
                    contact: p.contact,
                    area_m2: p.geometry.area_m2,
                })
                .collect::<Vec<_>>();
            let (selected, selected_report) = state
                .transport_paths_step_with_gas_law_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &surfaces,
                    paths,
                    admission(),
                    coupled(),
                    &[0.0; 2],
                    |contact, thermal, source, gas| {
                        let patch = *patches.iter().find(|p| p.contact == contact).unwrap();
                        prepared.resolve_conductive_gas_contact(patch, thermal, source, gas)
                    },
                )
                .unwrap();
            assert_eq!(selected.chambers(), actual.chambers());
            assert_eq!(selected.circuit(), actual.circuit());
            assert_eq!(selected.thermal().cells(), actual.thermal().cells());
            assert_eq!(selected_report.flow, report.flow);

            assert_eq!(state.time_s(), 0.0);
        })
    });
}

#[test]
fn physical_layer_resistance_changes_live_heat_charge_motion_and_sound() {
    let run = |length_scale| {
        fixture([0.8, 0.8], 150.0, |assembly, thermal, circuit| {
            with_cx(|cx, _| {
                let d = pair_card(2e-5);
                let q = QueryPoint::new();
                let pairs = binding(&d, &q);
                let medium = laminate(length_scale);
                let layers = [laminate_link(&medium, &pairs, &q)];
                let prepared = assembly.prepare(cx, &thermal).unwrap();
                let mut state = MaterialExchangeState::from_parts(thermal, circuit, chambers());
                let drives = drives(state.chambers());
                let before = storage(&state);
                let mut work = 0.0;
                for _ in 0..8 {
                    let target = state.time_s() + DT;
                    let report = state
                        .advance_with_transport_paths(
                            &prepared,
                            cx,
                            target,
                            &[],
                            &drives,
                            &[],
                            &gas_fixture::patches(),
                            GasTransportPaths {
                                layered_dusty_gas: &layers,
                                ..Default::default()
                            },
                            admission(),
                            coupled(),
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
                    assert!(
                        report.total_energy_residual_j.abs() <= coupled().total_energy_tolerance_j
                    );
                }
                assert!((storage(&state) + work - before).abs() < 4e-5);
                let gas = &state.chambers()[&10];
                assert_eq!(gas.gas().model().thermodynamics().components().len(), 2);
                (
                    gas.gas().mass_kg(),
                    state.thermal().temperature_k(CellId(30)).unwrap(),
                    state.circuit().charges_c()[0],
                    gas.position_m(),
                    gas.gas().state().state().sound_speed,
                )
            })
        })
    };
    let thin = run(1.0);
    let thick = run(4.0);
    assert_ne!(thin.0, thick.0);
    assert_ne!(thin.1, thick.1);
    assert_ne!(thin.2, thick.2);
    assert_ne!(thin.3, thick.3);
    assert_ne!(thin.4, thick.4);
}

#[test]
fn late_selected_source_failure_and_cancellation_publish_no_laminate_exchange() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, gate| {
            let d = pair_card(2e-5);
            let q = QueryPoint::new();
            let pairs = binding(&d, &q);
            let medium = laminate(1.0);
            let layers = [laminate_link(&medium, &pairs, &q)];
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut state = MaterialExchangeState::from_parts(thermal, circuit, chambers());
            let drives = drives(state.chambers());
            let cells = state.thermal().cells().clone();
            let circuit = state.circuit().clone();
            let gases = state.chambers().clone();
            let surfaces = gas_fixture::patches()
                .iter()
                .map(|p| GasSurface {
                    contact: p.contact,
                    area_m2: p.geometry.area_m2,
                })
                .collect::<Vec<_>>();
            let calls = std::cell::Cell::new(0);
            let result = state.advance_with_transport_paths_and_gas_law(
                &prepared,
                cx,
                DT,
                &[],
                &drives,
                &[],
                &surfaces,
                GasTransportPaths {
                    layered_dusty_gas: &layers,
                    ..Default::default()
                },
                admission(),
                coupled(),
                &[0.0; 2],
                |_, _, _, _| -> Result<ConductanceSample<()>, &'static str> {
                    calls.set(calls.get() + 1);
                    Err("synthetic late material heat-source refusal")
                },
            );
            assert!(result.is_err());
            assert!(calls.get() > 0);
            assert_eq!(state.chambers(), &gases);
            assert_eq!(state.circuit(), &circuit);
            assert_eq!(state.thermal().cells(), &cells);
            assert_eq!(state.time_s(), 0.0);
            gate.request();
            assert!(matches!(
                state.advance_with_transport_paths(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    GasTransportPaths {
                        layered_dusty_gas: &layers,
                        ..Default::default()
                    },
                    admission(),
                    coupled(),
                    &[0.0; 2],
                ),
                Err(MaterialFlowError::Cancelled)
            ));
            assert_eq!(state.chambers(), &gases);
            assert_eq!(state.time_s(), 0.0);
        })
    });
}
