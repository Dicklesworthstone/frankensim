//! Complete adaptive bath trajectories, independent decay and atomic histories.
use super::*;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::state::flow::{
    diffusion::coupled::{CoupledDiffusionControls, CoupledDiffusionScales},
    molecular::reservoir::coupled::*,
};

fn loose_time(maximum_step_s: f64) -> CoupledMolecularReservoirControls {
    CoupledMolecularReservoirControls {
        state: CoupledDiffusionControls {
            absolute: CoupledDiffusionScales {
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
            maximum_step_s,
            minimum_step_s: maximum_step_s / 65536.0,
            maximum_evaluations: 512,
        },
        streams: MolecularReservoirStreamScales {
            species_mass_kg: 1.0,
            heat_j: 1e6,
        },
    }
}

fn fixed_budgets(balance: &MolecularReservoirWindowBalance) {
    let physical = molecular_controls().physical;
    let coupled = controls();
    // These are the ORIGINAL fixed tolerances, never multiplied by the number
    // of accepted substeps or eased when the time-accuracy controls get tighter.
    for (defect, limit) in [
        (
            balance.transport_energy_defect_bound_j,
            physical.exchange.energy_tolerance_j,
        ),
        (
            balance.feedback_energy_defect_bound_j,
            coupled.total_energy_tolerance_j,
        ),
        (
            balance.total_energy_defect_bound_j,
            coupled.total_energy_tolerance_j,
        ),
        (
            balance.maximum_cell_residual_bound_j,
            coupled.network.cell_energy_tolerance_j,
        ),
        (
            balance.thermal_energy_defect_bound_j,
            coupled.network.thermal_energy_tolerance_j,
        ),
        (
            balance.circuit_energy_defect_bound_j,
            coupled.network.circuit.energy_tolerance_j,
        ),
        (
            balance.network_energy_defect_bound_j,
            coupled.network.total_energy_tolerance_j,
        ),
        (
            balance.maximum_chamber_energy_residual_j,
            coupled.total_energy_tolerance_j,
        ),
        (
            balance.transport_entropy_deficit_bound_j_k,
            physical.entropy_tolerance_j_k,
        ),
        (
            balance.maximum_species_mass_residual_kg,
            physical.exchange.mass_tolerance_kg,
        ),
    ] {
        assert!(
            defect.is_finite() && defect >= 0.0 && defect <= limit,
            "whole-window defect {defect:e} exceeds unchanged budget {limit:e}"
        );
    }
    assert!(balance.total_energy_residual_j.abs() <= coupled.total_energy_tolerance_j);
    assert!(balance.transport_energy_residual_j.abs() <= physical.exchange.energy_tolerance_j);
    assert!(balance.network_energy_residual_j.abs() <= coupled.network.total_energy_tolerance_j);
    assert!(balance.total_mass_residual_kg.abs() <= physical.exchange.mass_tolerance_kg);
    assert!(balance.transport_entropy_production_j_k >= -physical.entropy_tolerance_j_k);
    for row in &balance.species {
        assert!(row.residual_kg.abs() <= physical.exchange.mass_tolerance_kg);
    }
    for row in &balance.chambers {
        assert!(row.energy_defect_bound_j <= coupled.total_energy_tolerance_j);
        assert!(row.energy_residual_j.abs() <= coupled.total_energy_tolerance_j);
    }
}

#[test]
fn adaptive_reservoir_accepts_exactly_two_original_complete_fine_steps() {
    fixture([0.2, 0.8], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, mut original) = molecular_fixture(MASS, 1.0, -50000.0);
            let (_, mut molecular) = molecular_fixture(MASS, 1.0, -50000.0);
            let held = [
                bath(&bank, 0.9, 430.0, 1.2e5),
                bath(&bank, 0.1, 370.0, 0.9e5),
            ];
            let retained = held.clone();
            let reservoirs = [
                ReservoirMolecularLink {
                    area_m2: 0.08,
                    ..bath_link(&held[0])
                },
                ReservoirMolecularLink {
                    id: 202,
                    parcel: 20,
                    area_m2: 0.06,
                    ..bath_link(&held[1])
                },
            ];
            let internal = [molecular_link()];
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(0.7, 390.0, 1.1e5), (0.3, 410.0, 0.95e5)], false),
            );
            let mut loads = drives(state.chambers());
            loads[0].heat_j = 0.015625;
            loads[1].heat_j = 0.03125;
            loads[0].external_force_n = 2.0;
            loads[1].external_force_n = -1.5;
            let half: Vec<_> = loads
                .iter()
                .map(|drive| ChamberDrive {
                    heat_j: 0.5 * drive.heat_j,
                    ..*drive
                })
                .collect();
            let law = reaction(50.0);
            let reactions = [ChamberReaction {
                chamber: 10,
                reaction: &law,
                controls: chemical(),
            }];
            let patches = gas_fixture::patches();
            let (middle, first) = state
                .molecular_reservoir_step_candidate(
                    &prepared,
                    cx,
                    DT / 2.0,
                    &[],
                    &half,
                    &reactions,
                    &patches,
                    &mut original,
                    &internal,
                    &reservoirs,
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let (expected, second) = middle
                .molecular_reservoir_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &half,
                    &reactions,
                    &patches,
                    &mut original,
                    &internal,
                    &reservoirs,
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let mut time = loose_time(DT);
            time.state.maximum_evaluations = 3;
            let (actual, report) = state
                .coupled_molecular_reservoir_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &reactions,
                    &patches,
                    &mut molecular,
                    &internal,
                    &reservoirs,
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                    time,
                )
                .unwrap();
            assert_eq!(actual.circuit(), expected.circuit());
            assert_eq!(actual.thermal().cells(), expected.thermal().cells());
            assert_eq!(report.steps[0].advance.flow, first.flow);
            assert_eq!(report.steps[1].advance.flow, second.flow);
            assert_eq!(report.steps[0].advance.exchange.step, first.exchange.step);
            assert_eq!(report.steps[1].advance.exchange.step, second.exchange.step);
            assert_eq!(report.evaluations, 3);
            assert_eq!(report.accepted_intervals, 1);
            assert_eq!(report.rejected_intervals, 0);
            assert_eq!(report.steps.len(), 2);
            assert_eq!(
                (report.steps[0].start_s, report.steps[0].end_s),
                (0.0, DT / 2.0)
            );
            assert_eq!(
                (report.steps[1].start_s, report.steps[1].end_s),
                (DT / 2.0, DT)
            );
            assert!(report.estimated_error.maximum() <= 1.0);
            close(report.balance.prescribed_gas_heat_j, 0.046875, 1e-16);
            close(
                report.balance.reaction_extents[0].extent_mol,
                first.exchange.report.source.reactions[0]
                    .reaction
                    .extent_mol
                    + second.exchange.report.source.reactions[0]
                        .reaction
                        .extent_mol,
                1e-16,
            );
            for (&id, chamber) in actual.chambers() {
                let oracle = &expected.chambers()[&id];
                assert_eq!(chamber.gas(), oracle.gas());
                assert_eq!(chamber.position_m(), oracle.position_m());
                assert_eq!(chamber.velocity_m_s(), oracle.velocity_m_s());
                assert_eq!(chamber.last_window_duration_s(), Some(DT));
                assert!(matches!(
                    chamber.audit_last_window().unwrap(),
                    ModuleWindowAudit::Unavailable { .. }
                ));
                let account = report
                    .balance
                    .chambers
                    .iter()
                    .find(|row| row.chamber == id)
                    .unwrap();
                let transport: f64 = [&first, &second]
                    .into_iter()
                    .map(|step| {
                        step.balance
                            .chambers
                            .iter()
                            .find(|row| row.chamber == id)
                            .unwrap()
                            .transport_energy_change_j
                    })
                    .sum();
                let feedback: f64 = [&first, &second]
                    .into_iter()
                    .map(|step| {
                        step.balance
                            .chambers
                            .iter()
                            .find(|row| row.chamber == id)
                            .unwrap()
                            .feedback_energy_change_j
                    })
                    .sum();
                close(account.transport_energy_change_j, transport, 2e-10);
                close(account.feedback_energy_change_j, feedback, 2e-10);
                assert!(account.transported_species && account.transported_heat);
            }
            for boundary in &report.balance.reservoir_heat {
                let net: f64 = [&first, &second]
                    .into_iter()
                    .map(|step| {
                        step.flow
                            .links
                            .iter()
                            .find(|link| link.id == boundary.id)
                            .unwrap()
                            .reduced_heat_into_gas_j
                    })
                    .sum();
                close(boundary.net_heat_into_gas_j, net, 2e-12);
                close(
                    boundary.heat_into_gas_j - boundary.heat_out_of_gas_j,
                    net,
                    2e-12,
                );
            }
            fixed_budgets(&report.balance);
            assert_eq!(held, retained);
            assert_eq!(state.time_s(), 0.0);
        });
    });
}

const CAPACITY: f64 = 0.625;
const WINDOW: f64 = 20.0 * DT;

fn thermal_reference(time_s: f64) -> (f64, f64) {
    // Pure A, fixed original volume, no wall patches: C=2.5*p*V/T=0.625 J/K.
    // The authored constant conductivity and exposed geometry give G=30 W/K.
    // This is the continuum ODE, independent of every native integration path.
    let rate = (CONDUCTIVITY * 0.1 / 0.0001) / CAPACITY;
    let heat = CAPACITY * 40.0 * -(-rate * time_s).exp_m1();
    (400.0 + heat / CAPACITY, heat)
}

fn pure_heat_account(
    before: &MaterialExchangeState<'_>,
    after: &MaterialExchangeState<'_>,
    report: &CoupledMolecularReservoirReport,
) -> f64 {
    let old = before.chambers()[&10].gas();
    let gas = after.chambers()[&10].gas();
    close(
        old.mass_kg() * 2.5 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / MASS,
        CAPACITY,
        1e-14,
    );
    close(gas.mass_kg(), old.mass_kg(), 1e-16);
    assert_eq!(gas.model().thermodynamics().components().len(), 1);
    assert_eq!(species_mass(gas, "B"), 0.0);
    assert!(report.balance.reservoir_species.is_empty());
    assert_eq!(report.balance.reservoir_heat.len(), 1);
    let boundary = &report.balance.reservoir_heat[0];
    assert_eq!((boundary.id, boundary.parcel), (101, 10));
    assert!(boundary.heat_into_gas_j > 0.0);
    assert_eq!(boundary.heat_out_of_gas_j, 0.0);
    close(
        boundary.net_heat_into_gas_j,
        boundary.heat_into_gas_j,
        1e-13,
    );
    let storage = CAPACITY * (gas.state().state().temperature - 400.0);
    close(boundary.net_heat_into_gas_j, storage, 2e-8);
    assert!(report.estimated_error.maximum() <= 1.0);
    assert_eq!(report.steps.len(), 2 * report.accepted_intervals);
    assert!(report.evaluations >= 3 * report.accepted_intervals);
    assert_eq!(after.time_s(), WINDOW);
    assert_ne!(after.circuit(), before.circuit());
    assert_eq!(before.time_s(), 0.0);
    fixed_budgets(&report.balance);
    boundary.net_heat_into_gas_j
}

#[test]
fn complete_bath_temperature_refinement_converges_to_independent_exponential_decay() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, _) = molecular_fixture(MASS, 1.0, 0.0);
            let held = bath(&bank, 1.0, 440.0, 1e5);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(1.0, 400.0, 1e5); 2], true),
            );
            let loads = drives(state.chambers());
            let exact = thermal_reference(WINDOW);
            let mut errors = Vec::new();
            let mut work = Vec::new();
            for tolerance in [0.03, 0.0003] {
                let (_, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
                let mut time = loose_time(WINDOW);
                time.state.absolute.gas_temperature_k = tolerance;
                let (next, report) = state
                    .coupled_molecular_reservoir_candidate(
                        &prepared,
                        cx,
                        WINDOW,
                        &[],
                        &loads,
                        &[],
                        &[],
                        &mut molecular,
                        &[],
                        &[bath_link(&held)],
                        molecular_controls(),
                        controls(),
                        &[0.0; 2],
                        time,
                    )
                    .unwrap();
                pure_heat_account(&state, &next, &report);
                errors
                    .push((next.chambers()[&10].gas().state().state().temperature - exact.0).abs());
                work.push(report.evaluations);
                if tolerance < 0.001 {
                    assert!(report.steps.len() > 2);
                }
            }
            assert!(
                errors[1] < 0.2 * errors[0],
                "temperature errors: {errors:?}"
            );
            assert!(errors[1] < 0.005, "temperature errors: {errors:?}");
            assert!(work[1] > work[0]);
        });
    });
}

#[test]
fn strict_bath_stream_accuracy_refines_even_when_all_state_scales_are_loose() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, _) = molecular_fixture(MASS, 1.0, 0.0);
            let held = bath(&bank, 1.0, 440.0, 1e5);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(1.0, 400.0, 1e5); 2], true),
            );
            let loads = drives(state.chambers());
            let exact = thermal_reference(WINDOW);
            let mut errors = Vec::new();
            let mut work = Vec::new();
            for tolerance in [1e6, 1e-4] {
                let (_, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
                let mut time = loose_time(WINDOW);
                time.streams.heat_j = tolerance;
                let (next, report) = state
                    .coupled_molecular_reservoir_candidate(
                        &prepared,
                        cx,
                        WINDOW,
                        &[],
                        &loads,
                        &[],
                        &[],
                        &mut molecular,
                        &[],
                        &[bath_link(&held)],
                        molecular_controls(),
                        controls(),
                        &[0.0; 2],
                        time,
                    )
                    .unwrap();
                let heat = pure_heat_account(&state, &next, &report);
                errors.push((heat - exact.1).abs());
                work.push(report.evaluations);
                assert!(report.estimated_error.state.maximum() < 1.0);
                assert!(report.estimated_error.streams.heat <= 1.0);
                if tolerance < 1.0 {
                    assert!(report.steps.len() > 2);
                }
            }
            assert!(
                errors[1] < 0.3 * errors[0],
                "bath-integral errors: {errors:?}"
            );
            assert!(errors[1] < 0.003, "bath-integral errors: {errors:?}");
            assert!(work[1] > work[0]);
        });
    });
}

#[test]
fn late_adaptive_work_limit_and_cancellation_publish_no_prefix_or_port_history() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, gate| {
            let (bank, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
            let held = bath(&bank, 1.0, 440.0, 1e5);
            let link = bath_link(&held);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(1.0, 400.0, 1e5); 2], true),
            );
            let loads = drives(state.chambers());
            let (preview, _) = state
                .coupled_molecular_reservoir_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &[],
                    &[],
                    &mut molecular,
                    &[],
                    &[link],
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                    loose_time(DT),
                )
                .unwrap();
            state
                .advance_with_coupled_molecular_reservoir(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &[],
                    &[],
                    &mut molecular,
                    &[],
                    &[link],
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                    loose_time(DT),
                )
                .unwrap();
            same_owners(&state, &preview);
            assert_eq!(state.chambers()[&10].last_window_duration_s(), Some(DT));
            let mut time = loose_time(DT);
            time.state.maximum_evaluations = 3;
            match state.advance_with_coupled_molecular_reservoir(
                &prepared,
                cx,
                3.0 * DT,
                &[],
                &loads,
                &[],
                &[],
                &mut molecular,
                &[],
                &[link],
                molecular_controls(),
                controls(),
                &[0.0; 2],
                time,
            ) {
                Err(CoupledMolecularReservoirError::WorkLimit {
                    time_s,
                    evaluations,
                    accepted_substeps,
                }) => {
                    assert_eq!(time_s, 2.0 * DT);
                    assert_eq!(evaluations, 3);
                    assert_eq!(accepted_substeps, 2);
                }
                other => panic!(
                    "expected a typed late adaptive work limit, got {}",
                    match other {
                        Ok(_) => "an accepted window".into(),
                        Err(error) => format!("{error:?}"),
                    }
                ),
            }
            same_owners(&state, &preview);
            gate.request();
            assert!(matches!(
                state.advance_with_coupled_molecular_reservoir(
                    &prepared,
                    cx,
                    2.0 * DT,
                    &[],
                    &loads,
                    &[],
                    &[],
                    &mut molecular,
                    &[],
                    &[link],
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                    loose_time(DT),
                ),
                Err(CoupledMolecularReservoirError::Cancelled)
            ));
            same_owners(&state, &preview);
        });
    });
}
