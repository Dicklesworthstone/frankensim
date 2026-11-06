//! G1/G4: stiff binary mixing reaches the complete original live assembly.
use super::*;
use diffusion::binary::{self as binary_transport, BinaryDiffusionNetworkAdvanceControls};

fn binary_time_controls() -> BinaryDiffusionNetworkAdvanceControls {
    let mut temporal = adaptive_controls();
    temporal.maximum_evaluations = 3;
    BinaryDiffusionNetworkAdvanceControls {
        temporal,
        maximum_pair_evaluations: 3,
    }
}

#[test]
fn g1_adaptive_binary_diffusion_matches_the_exponential_and_original_assembly_owners() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let card = pair_card(0.2);
            let context = QueryPoint::new();
            let pairs = [binding(&card, &context)];
            let links = [face(&pairs)];
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(thermal, circuit, mixtures());
            let drives = drives(state.chambers());
            let mut physical = diffusion_controls();
            physical.exchange.maximum_species_fraction = 0.9;
            // k*DT=10 per finite volume. The explicit transfer overdraws its
            // donor even with this allowance; exact pair relaxation stays finite.
            assert!(matches!(
                state.diffusion_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &links,
                    physical,
                    controls(),
                    &[0.0; 2],
                ),
                Err(MaterialFlowError::Diffusion(_))
            ));
            let parcels = state
                .chambers()
                .iter()
                .map(|(&id, gas)| (id, gas.gas().clone()))
                .collect();
            let endpoints = drives
                .iter()
                .map(|drive| (drive.chamber, drive.piston_solve.gas))
                .collect();
            let (diffused, native) = binary_transport::advance_adaptive_candidate(
                &parcels,
                &links,
                DT,
                &endpoints,
                physical,
                binary_time_controls(),
                || false,
            )
            .unwrap()
            .into_parts();
            let fraction = |id| {
                diffused[&id]
                    .model()
                    .thermodynamics()
                    .mole_composition()
                    .fractions()[0]
            };
            // Independent two-cell finite-volume mode at equal T, P and volume:
            // gap(t)=gap(0)*exp[-D*A*(1/Vl+1/Vr)*t/distance].
            let gap = 0.5 * (-20.0_f64).exp();
            assert!(((fraction(10) - fraction(20)) - gap).abs() < 2e-12);
            let staged = MaterialExchangeState::from_parts(
                state.thermal().clone(),
                state.circuit().clone(),
                diffused
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
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let (actual, report) = state
                .adaptive_binary_diffusion_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &links,
                    physical,
                    binary_time_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            assert_eq!(report.flow, native);
            assert_eq!(actual.chambers(), expected.chambers());
            assert_eq!(actual.circuit(), expected.circuit());
            assert_eq!(actual.thermal().cells(), expected.thermal().cells());
            assert_eq!(actual.time_s(), DT);
            assert_eq!(report.flow.duration_s, DT);
            assert_eq!(report.flow.accepted_intervals, 1);
            assert_eq!(report.flow.accepted_substeps, 2);
            assert_eq!(report.flow.rejected_intervals, 0);
            assert_eq!(report.flow.evaluations, 3);
            assert_eq!(report.flow.pair_evaluations, 3);
            assert_eq!(report.flow.last_substep.len(), 1);
            assert_eq!(report.flow.last_substep[0].parcels, [10, 20]);
            assert!(report.flow.maximum_species_error_ratio <= 1.0);
            assert!(report.flow.maximum_temperature_error_ratio <= 1.0);
            assert!(
                report.flow.maximum_species_mass_residual_kg <= physical.exchange.mass_tolerance_kg
            );
            assert!(
                report.flow.endpoint_energy_change_j.abs() <= physical.exchange.energy_tolerance_j
            );
            assert!(report.flow.entropy_change_j_k > 0.0);
            assert!(report.total_mass_residual_kg.abs() <= physical.exchange.mass_tolerance_kg);
            assert!(report.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j);
            let before_mass: f64 = state
                .chambers()
                .values()
                .map(|gas| gas.gas().mass_kg())
                .sum();
            let after_mass: f64 = actual
                .chambers()
                .values()
                .map(|gas| gas.gas().mass_kg())
                .sum();
            assert!((after_mass - before_mass).abs() <= physical.exchange.mass_tolerance_kg);
            let work: f64 = report
                .exchange
                .report
                .source
                .chambers
                .iter()
                .map(|chamber| chamber.boundary.piston.exterior_pressure_work_j)
                .sum();
            assert!((storage(&actual) + work - storage(&state)).abs() < 4e-6);
            for (&id, gas) in actual.chambers() {
                assert!(
                    report.flow.transfers.iter().any(|row| {
                        row.mass_kg > 0.0 && (row.donor == id || row.receiver == id)
                    })
                );
                assert_eq!(gas.last_window_duration_s(), Some(DT));
                assert!(matches!(
                    gas.audit_last_window().unwrap(),
                    ModuleWindowAudit::Unavailable { .. }
                ));
            }
            assert_eq!(state.time_s(), 0.0);
        });
    });
}

#[test]
fn g1_adaptive_binary_reaction_uses_actual_diffused_inventory_once_at_zero_net_mass() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let law = reaction(50.0);
            let card = pair_card(0.2);
            let context = QueryPoint::new();
            let pairs = [binding(&card, &context)];
            let links = [face(&pairs)];
            let gases = [0.75, 0.25]
                .into_iter()
                .enumerate()
                .map(|(i, x)| {
                    let model = law
                        .catalog()
                        .with_composition(
                            Composition::new(CompositionBasis::MoleFraction, vec![x, 1.0 - x])
                                .unwrap(),
                        )
                        .unwrap();
                    (10 + 10 * i as u64, parcel_piston(model, 100000.0))
                })
                .collect();
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(thermal, circuit, gases);
            let drives = drives(state.chambers());
            let mut physical = diffusion_controls();
            physical.exchange.maximum_species_fraction = 0.9;
            let parcels = state
                .chambers()
                .iter()
                .map(|(&id, gas)| (id, gas.gas().clone()))
                .collect();
            let endpoints = drives
                .iter()
                .map(|drive| (drive.chamber, drive.piston_solve.gas))
                .collect();
            let diffused = binary_transport::advance_adaptive_candidate(
                &parcels,
                &links,
                DT,
                &endpoints,
                physical,
                binary_time_controls(),
                || false,
            )
            .unwrap();
            let expected = law
                .parcel_candidate(&diffused.parcels()[&10], DT, chemical(), || false)
                .unwrap();
            let unmixed = law
                .parcel_candidate(state.chambers()[&10].gas(), DT, chemical(), || false)
                .unwrap();
            assert_ne!(expected.step().extent_mol, unmixed.step().extent_mol);
            let (next, report) = state
                .adaptive_binary_diffusion_step_candidate(
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
                    physical,
                    binary_time_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            assert_eq!(report.flow, *diffused.report());
            assert_eq!(report.exchange.report.source.reactions.len(), 1);
            assert_eq!(
                report.exchange.report.source.reactions[0].reaction,
                *expected.step()
            );
            assert!(report.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j);
            let work: f64 = report
                .exchange
                .report
                .source
                .chambers
                .iter()
                .map(|chamber| chamber.boundary.piston.exterior_pressure_work_j)
                .sum();
            assert!((storage(&next) + work - storage(&state)).abs() < 4e-6);
            // Equal-molar-mass A/B counterdiffusion has zero net chamber mass,
            // but each real gross transfer must still open that chamber's audit.
            for (&id, before) in state.chambers() {
                let incoming: f64 = report
                    .flow
                    .transfers
                    .iter()
                    .filter(|row| row.receiver == id)
                    .map(|row| row.mass_kg)
                    .sum();
                let outgoing: f64 = report
                    .flow
                    .transfers
                    .iter()
                    .filter(|row| row.donor == id)
                    .map(|row| row.mass_kg)
                    .sum();
                assert!(incoming > 0.0 && outgoing > 0.0);
                assert!((incoming - outgoing).abs() < 1e-15);
                assert!(
                    (next.chambers()[&id].gas().mass_kg() - before.gas().mass_kg()).abs()
                        <= physical.exchange.mass_tolerance_kg
                );
                assert_eq!(next.chambers()[&id].last_window_duration_s(), Some(DT));
                assert!(matches!(
                    next.chambers()[&id].audit_last_window().unwrap(),
                    ModuleWindowAudit::Unavailable { .. }
                ));
            }
            assert_eq!(next.time_s(), DT);
            assert_eq!(state.time_s(), 0.0);
        });
    });
}

#[test]
fn g4_adaptive_binary_work_source_feedback_and_cancellation_publish_nothing() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, gate| {
            let card = pair_card(0.2);
            let context = QueryPoint::new();
            let pairs = [binding(&card, &context)];
            let links = [face(&pairs)];
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut state = MaterialExchangeState::from_parts(thermal, circuit, mixtures());
            let drives = drives(state.chambers());
            let saved_gas = state.chambers().clone();
            let saved_cells = state.thermal().cells().clone();
            let saved_circuit = state.circuit().clone();
            let mut physical = diffusion_controls();
            physical.exchange.maximum_species_fraction = 0.9;
            let limited = BinaryDiffusionNetworkAdvanceControls {
                maximum_pair_evaluations: 2,
                ..binary_time_controls()
            };
            assert!(matches!(
                state.advance_with_adaptive_binary_diffusion(
                    &prepared, cx, DT, &[], &drives, &[], &gas_fixture::patches(),
                    &links, physical, limited, controls(), &[0.0; 2],
                ),
                Err(MaterialFlowError::AdaptiveDiffusion(
                    DiffusionNetworkAdvanceError::PairWorkLimit {
                        elapsed_s,
                        accepted_intervals: 0,
                        pair_evaluations: 2,
                        maximum_pair_evaluations: 2,
                    }
                )) if elapsed_s == 0.0
            ));
            let missing = DiffusionNetworkLink {
                pairs: &[],
                ..links[0].clone()
            };
            assert!(matches!(
                state.advance_with_adaptive_binary_diffusion(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &[missing],
                    physical,
                    binary_time_controls(),
                    controls(),
                    &[0.0; 2],
                ),
                Err(MaterialFlowError::AdaptiveDiffusion(
                    DiffusionNetworkAdvanceError::BinaryStep { .. }
                ))
            ));
            // Transport itself succeeds; refusal occurs in the later existing
            // electrical/thermal feedback solve.
            let mut late = controls();
            late.network.maximum_circuit_trials = 1;
            assert!(matches!(
                state.advance_with_adaptive_binary_diffusion(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &links,
                    physical,
                    binary_time_controls(),
                    late,
                    &[0.0; 2],
                ),
                Err(MaterialFlowError::Exchange(_))
            ));
            assert_eq!(state.chambers(), &saved_gas);
            assert_eq!(state.circuit(), &saved_circuit);
            assert_eq!(state.thermal().cells(), &saved_cells);
            assert_eq!(state.time_s(), 0.0);
            let (preview, _) = state
                .adaptive_binary_diffusion_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &links,
                    physical,
                    binary_time_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            state
                .advance_with_adaptive_binary_diffusion(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &links,
                    physical,
                    binary_time_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            assert_eq!(state.chambers(), preview.chambers());
            assert_eq!(state.circuit(), preview.circuit());
            assert_eq!(state.thermal().cells(), preview.thermal().cells());
            assert_eq!(state.time_s(), DT);
            gate.request();
            assert!(matches!(
                state.advance_with_adaptive_binary_diffusion(
                    &prepared,
                    cx,
                    2.0 * DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &links,
                    physical,
                    binary_time_controls(),
                    controls(),
                    &[0.0; 2],
                ),
                Err(MaterialFlowError::Cancelled)
            ));
            assert_eq!(state.chambers(), preview.chambers());
            assert_eq!(state.circuit(), preview.circuit());
            assert_eq!(state.thermal().cells(), preview.thermal().cells());
            assert_eq!(state.time_s(), DT);
        });
    });
}
