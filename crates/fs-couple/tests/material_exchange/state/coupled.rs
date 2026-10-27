//! Persistent advancement moves the original mixed-boundary candidates once.
use super::*;
use crate::coupled_core::{binding, implicit_reaction};
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::exchange::{
    gas_law::GasSurface, reactive::CoupledMaterialExchangeError,
};
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::reactive::CoupledChamberStep;
use std::cell::Cell;

fn same_owned_state(a: &MaterialExchangeState<'_>, b: &MaterialExchangeState<'_>) {
    assert_eq!(a.time_s(), b.time_s());
    assert_eq!(a.thermal().cells(), b.thermal().cells());
    assert_eq!(a.circuit(), b.circuit());
    assert_eq!(a.chambers(), b.chambers());
    for (&id, gas) in a.chambers() {
        assert_eq!(
            gas.last_window_duration_s(),
            b.chambers()[&id].last_window_duration_s()
        );
        assert_eq!(
            gas.last_window_transfers(),
            b.chambers()[&id].last_window_transfers()
        );
    }
}

#[test]
fn continued_live_state_matches_direct_candidates_without_repeating_reaction_or_ports() {
    fixture([0.8, 0.8], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let law = implicit_reaction(10000.0);
            let mut state = MaterialExchangeState::from_parts(thermal, circuit, gas_map(&law));
            let drives = [
                loads(20, &state.chambers()[&20]),
                loads(10, &state.chambers()[&10]),
            ];
            let reactions = [binding(&law)];
            for _ in 0..2 {
                let target = state.time_s() + 0.00005;
                let dt = target - state.time_s();
                let cells = state.thermal().cells().clone();
                let electrical = state.circuit().clone();
                let gases = state.chambers().clone();
                let time = state.time_s();
                let chambers = [
                    gas_fixture::chamber(20, &state.chambers()[&20], dt),
                    gas_fixture::chamber(10, &state.chambers()[&10], dt),
                ];
                let direct = prepared
                    .coupled_reactive_step_candidate(
                        cx,
                        state.thermal(),
                        target,
                        state.circuit(),
                        &[],
                        &chambers,
                        &reactions,
                        &gas_fixture::patches(),
                        controls(),
                        &[0.0; 2],
                    )
                    .unwrap();
                let (preview, preview_report) = state
                    .coupled_reactive_candidate(
                        &prepared,
                        cx,
                        target,
                        &[],
                        &drives,
                        &reactions,
                        &gas_fixture::patches(),
                        controls(),
                        &[0.0; 2],
                    )
                    .unwrap();
                assert_eq!(state.time_s(), time);
                assert_eq!(state.thermal().cells(), &cells);
                assert_eq!(state.circuit(), &electrical);
                assert_eq!(state.chambers(), &gases);
                let accepted = state
                    .advance_coupled_reactive(
                        &prepared,
                        cx,
                        target,
                        &[],
                        &drives,
                        &reactions,
                        &gas_fixture::patches(),
                        controls(),
                        &[0.0; 2],
                    )
                    .unwrap();
                same_owned_state(&state, &preview);
                assert_eq!(state.time_s(), target);
                assert_eq!(state.thermal().cells(), direct.thermal().cells());
                assert_eq!(state.circuit(), direct.circuit());
                assert_eq!(accepted.step, *direct.step());
                assert_eq!(accepted.step, preview_report.step);
                for (account, candidate) in accepted
                    .report
                    .source
                    .chambers
                    .iter()
                    .zip(direct.chambers())
                {
                    assert_eq!(account.id, candidate.id);
                    assert_eq!(&state.chambers()[&account.id], candidate.boundary.piston());
                    assert_eq!(account.surfaces.as_slice(), candidate.boundary.source());
                    match &account.boundary {
                        CoupledChamberStep::Reactive(step) => {
                            assert_eq!(Some(step), candidate.boundary.reactive_step());
                            assert!(step.piston.reaction.reaction.extent_mol > 0.0);
                        }
                        CoupledChamberStep::Nonreactive(step) => {
                            assert_eq!(Some(step), candidate.boundary.nonreactive_step());
                        }
                    }
                    assert_eq!(
                        state.chambers()[&account.id].last_window_duration_s(),
                        Some(dt)
                    );
                    assert_eq!(
                        state.chambers()[&account.id].last_window_transfers(),
                        candidate.boundary.piston().last_window_transfers()
                    );
                }
                assert_eq!(
                    accepted.report.electrical.contacts,
                    preview_report.report.electrical.contacts
                );
                assert!(
                    accepted.step.total_energy_residual_j.abs()
                        <= controls().total_energy_tolerance_j
                );
                let (thermal, circuit, chambers) = state.into_parts();
                state = MaterialExchangeState::from_parts(thermal, circuit, chambers);
                same_owned_state(&state, &preview);
            }
        })
    });
}

#[test]
fn empty_reactions_reduce_to_native_and_failed_live_updates_retry_without_history_loss() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, gate| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let law = implicit_reaction(10000.0);
            let mut state = MaterialExchangeState::from_parts(thermal, circuit, gas_map(&law));
            let drives = [
                loads(10, &state.chambers()[&10]),
                loads(20, &state.chambers()[&20]),
            ];
            let reactions = [binding(&law)];
            let chambers = [
                gas_fixture::chamber(10, &state.chambers()[&10], 0.00005),
                gas_fixture::chamber(20, &state.chambers()[&20], 0.00005),
            ];
            let direct = prepared
                .step_candidate(
                    cx,
                    state.thermal(),
                    0.00005,
                    state.circuit(),
                    &[],
                    &chambers,
                    &gas_fixture::patches(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let ordinary = state
                .advance_coupled_reactive(
                    &prepared,
                    cx,
                    0.00005,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            assert_eq!(state.thermal().cells(), direct.thermal().cells());
            assert_eq!(state.circuit(), direct.circuit());
            assert_eq!(ordinary.step, *direct.step());
            for (account, candidate) in ordinary
                .report
                .source
                .chambers
                .iter()
                .zip(direct.chambers())
            {
                assert_eq!(&state.chambers()[&account.id], candidate.boundary.piston());
                let CoupledChamberStep::Nonreactive(step) = &account.boundary else {
                    panic!("an omitted reaction cannot change the native integration law");
                };
                assert_eq!(step, candidate.boundary.step());
            }

            let target = state.time_s() + 0.00005;
            let time = state.time_s();
            let cells = state.thermal().cells().clone();
            let circuit = state.circuit().clone();
            let gases = state.chambers().clone();
            let patches = gas_fixture::patches();
            let surfaces = patches.map(|patch| GasSurface {
                contact: patch.contact,
                area_m2: patch.geometry.area_m2,
            });
            let calls = Cell::new(0usize);
            let (reference, _) = state
                .coupled_reactive_with_gas_law_candidate(
                    &prepared,
                    cx,
                    target,
                    &[],
                    &drives,
                    &reactions,
                    &surfaces,
                    controls(),
                    &[0.0; 2],
                    |contact, thermal, source, gas| {
                        calls.set(calls.get() + 1);
                        prepared.resolve_conductive_gas_contact(
                            *patches.iter().find(|p| p.contact.id == contact.id).unwrap(),
                            thermal,
                            source,
                            gas,
                        )
                    },
                )
                .unwrap();
            assert!(calls.get() > 1);
            let mut limited = controls();
            limited.network.maximum_circuit_trials = 1;
            assert!(
                state
                    .advance_coupled_reactive(
                        &prepared,
                        cx,
                        target,
                        &[],
                        &drives,
                        &reactions,
                        &patches,
                        limited,
                        &[0.0; 2]
                    )
                    .is_err()
            );
            assert_eq!(state.time_s(), time);
            assert_eq!(state.thermal().cells(), &cells);
            assert_eq!(state.circuit(), &circuit);
            assert_eq!(state.chambers(), &gases);
            let observed = Cell::new(0usize);
            assert!(matches!(
                state.advance_coupled_reactive_with_gas_law(
                    &prepared,
                    cx,
                    target,
                    &[],
                    &drives,
                    &reactions,
                    &surfaces,
                    controls(),
                    &[0.0; 2],
                    |contact, thermal, source, gas| {
                        observed.set(observed.get() + 1);
                        if observed.get() == calls.get() {
                            gate.request();
                        }
                        prepared.resolve_conductive_gas_contact(
                            *patches.iter().find(|p| p.contact.id == contact.id).unwrap(),
                            thermal,
                            source,
                            gas,
                        )
                    }
                ),
                Err(CoupledMaterialExchangeError::Cancelled)
            ));
            assert_eq!(observed.get(), calls.get());
            assert_eq!(state.time_s(), time);
            assert_eq!(state.thermal().cells(), &cells);
            assert_eq!(state.circuit(), &circuit);
            assert_eq!(state.chambers(), &gases);
            for (&id, gas) in state.chambers() {
                assert_eq!(
                    gas.last_window_duration_s(),
                    gases[&id].last_window_duration_s()
                );
                assert_eq!(
                    gas.last_window_transfers(),
                    gases[&id].last_window_transfers()
                );
            }
            // A fresh execution scope retries the same still-accepted physical
            // interval; cancellation does not become a new chemistry/port clock.
            with_cx(|retry_cx, _| {
                state
                    .advance_coupled_reactive(
                        &prepared,
                        retry_cx,
                        target,
                        &[],
                        &drives,
                        &reactions,
                        &patches,
                        controls(),
                        &[0.0; 2],
                    )
                    .unwrap();
            });
            same_owned_state(&state, &reference);
        })
    });
}
