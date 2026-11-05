//! Complete physical feedback and time refinement through the shared controller.
use super::*;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::state::flow::molecular::coupled::*;

fn loose_time() -> CoupledMolecularControls {
    CoupledMolecularControls {
        absolute: CoupledMolecularScales {
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
fn coupled_cross_diffusion_accepts_the_native_fine_trajectory_and_heat_once() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, mut original) = molecular_fixture(MASS, 1.0, 0.0);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(0.5, 390.0, 1e5), (0.5, 410.0, 1e5)], false),
            );
            let mut drives = drives(state.chambers());
            drives[0].heat_j = 0.01;
            drives[1].heat_j = 0.02;
            drives[0].external_force_n = 2.0;
            let half: Vec<_> = drives
                .iter()
                .map(|drive| ChamberDrive {
                    heat_j: 0.5 * drive.heat_j,
                    ..*drive
                })
                .collect();
            let links = [molecular_link()];
            let paths = GasTransportPaths::default();
            let patches = gas_fixture::patches();
            let (middle, first) = state
                .molecular_transport_step_candidate(
                    &prepared,
                    cx,
                    DT / 2.0,
                    &[],
                    &half,
                    &[],
                    &patches,
                    &mut original,
                    &links,
                    paths,
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let (expected, second) = middle
                .molecular_transport_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &half,
                    &[],
                    &patches,
                    &mut original,
                    &links,
                    paths,
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let (_, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
            let (actual, report) = state
                .coupled_molecular_transport_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &patches,
                    &mut molecular,
                    &links,
                    paths,
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                    CoupledMolecularControls {
                        maximum_evaluations: 3,
                        ..loose_time()
                    },
                )
                .unwrap();
            assert_eq!(actual.circuit(), expected.circuit());
            assert_eq!(actual.thermal().cells(), expected.thermal().cells());
            for (&id, chamber) in actual.chambers() {
                assert_eq!(chamber.gas(), expected.chambers()[&id].gas());
                assert_eq!(chamber.position_m(), expected.chambers()[&id].position_m());
                assert_eq!(
                    chamber.velocity_m_s(),
                    expected.chambers()[&id].velocity_m_s()
                );
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
                let accepted_heat: f64 = [&first, &second]
                    .into_iter()
                    .map(|step| {
                        step.flow
                            .transport
                            .exchange
                            .parcels
                            .iter()
                            .find(|row| row.id == id)
                            .unwrap()
                            .net_additional_heat_j
                    })
                    .sum();
                close(account.net_diffusive_heat_j, accepted_heat, 1e-12);
                assert!(account.transported_species && account.transported_heat);
                assert!(account.energy_residual_j.abs() <= controls().total_energy_tolerance_j);
            }
            assert_eq!(report.steps[0].advance.flow, first.flow);
            assert_eq!(report.steps[1].advance.flow, second.flow);
            assert_ne!(
                first.flow.links[0].face.enthalpy_temperature_k,
                second.flow.links[0].face.enthalpy_temperature_k
            );
            assert_eq!(report.evaluations, 3);
            assert_eq!(report.accepted_intervals, 1);
            assert_eq!(report.steps.len(), 2);
            assert!(report.estimated_error.maximum() <= 1.0);
            close(report.balance.prescribed_gas_heat_j, 0.03, 1e-16);
            assert!(
                report.balance.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j
            );
            assert!(
                report.balance.maximum_species_mass_residual_kg
                    <= molecular_controls().physical.exchange.mass_tolerance_kg
            );
            assert_eq!(state.time_s(), 0.0);
        });
    });
}

#[test]
fn coupled_molecular_refinement_converges_to_independent_binary_diffusion_decay() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, _) = molecular_fixture(MASS, 0.0, 0.0);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(0.8, 400.0, 1e5), (0.2, 400.0, 1e5)], true),
            );
            let drives = drives(state.chambers());
            let link = molecular_link();
            // Identical calorics, equal masses, fixed volumes and no wall heat:
            // x_left - x_right = 0.6 exp(-2 D A t / (V L)).
            let exact =
                0.6 * (-2.0 * DIFFUSIVITY * link.area_m2 * DT / (0.001 * link.distance_m)).exp();
            let mut errors = Vec::new();
            let mut counts = Vec::new();
            for tolerance in [1e-9, 1e-10, 1e-11] {
                let (_, mut molecular) = molecular_fixture(MASS, 0.0, 0.0);
                let mut time = loose_time();
                time.absolute.species_mass_kg = tolerance;
                let (next, report) = state
                    .coupled_molecular_transport_candidate(
                        &prepared,
                        cx,
                        DT,
                        &[],
                        &drives,
                        &[],
                        &[],
                        &mut molecular,
                        &[link],
                        GasTransportPaths::default(),
                        molecular_controls(),
                        controls(),
                        &[0.0; 2],
                        time,
                    )
                    .unwrap();
                let x = |id| {
                    next.chambers()[&id]
                        .gas()
                        .model()
                        .thermodynamics()
                        .mole_composition()
                        .fractions()[0]
                };
                errors.push((x(10) - x(20) - exact).abs());
                counts.push(report.accepted_intervals);
                assert!(report.estimated_error.maximum() <= 1.0);
                assert!(
                    report.balance.total_energy_residual_j.abs()
                        <= controls().total_energy_tolerance_j
                );
                assert!(
                    report.balance.maximum_species_mass_residual_kg
                        <= molecular_controls().physical.exchange.mass_tolerance_kg
                );
                assert_eq!(next.time_s(), DT);
            }
            for i in 1..errors.len() {
                assert!(
                    errors[i] < 0.7 * errors[i - 1],
                    "diffusion errors: {errors:?}"
                );
                assert!(counts[i] > counts[i - 1], "accepted intervals: {counts:?}");
            }
            assert!(errors[2] < 2e-7, "diffusion errors: {errors:?}");
        });
    });
}

#[test]
fn coupled_molecular_late_work_limit_and_cancellation_publish_no_physical_prefix() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, gate| {
            let (bank, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(0.5, 390.0, 1e5), (0.5, 410.0, 1e5)], false),
            );
            let drives = drives(state.chambers());
            let old_chambers = state.chambers().clone();
            let old_cells = state.thermal().cells().clone();
            let old_circuit = state.circuit().clone();
            let time = CoupledMolecularControls {
                maximum_step_s: DT / 2.0,
                maximum_evaluations: 4,
                ..loose_time()
            };
            let result = state.advance_with_coupled_molecular_transport(
                &prepared,
                cx,
                DT,
                &[],
                &drives,
                &[],
                &gas_fixture::patches(),
                &mut molecular,
                &[molecular_link()],
                GasTransportPaths::default(),
                molecular_controls(),
                controls(),
                &[0.0; 2],
                time,
            );
            assert!(matches!(
                result,
                Err(CoupledMolecularError::WorkLimit {
                    evaluations: 4,
                    accepted_substeps: 2,
                    ..
                })
            ));
            assert_eq!(state.chambers(), &old_chambers);
            assert_eq!(state.thermal().cells(), &old_cells);
            assert_eq!(state.circuit(), &old_circuit);
            assert_eq!(state.time_s(), 0.0);
            for (&id, chamber) in state.chambers() {
                assert_eq!(
                    chamber.last_window_transfers(),
                    old_chambers[&id].last_window_transfers()
                );
                assert_eq!(
                    chamber.last_window_duration_s(),
                    old_chambers[&id].last_window_duration_s()
                );
            }
            gate.request();
            let result = state.advance_with_coupled_molecular_transport(
                &prepared,
                cx,
                DT,
                &[],
                &drives,
                &[],
                &gas_fixture::patches(),
                &mut molecular,
                &[molecular_link()],
                GasTransportPaths::default(),
                molecular_controls(),
                controls(),
                &[0.0; 2],
                loose_time(),
            );
            assert!(matches!(result, Err(CoupledMolecularError::Cancelled)));
            assert_eq!(state.chambers(), &old_chambers);
            assert_eq!(state.thermal().cells(), &old_cells);
            assert_eq!(state.circuit(), &old_circuit);
            assert_eq!(state.time_s(), 0.0);
            for (&id, chamber) in state.chambers() {
                assert_eq!(
                    chamber.last_window_transfers(),
                    old_chambers[&id].last_window_transfers()
                );
                assert_eq!(
                    chamber.last_window_duration_s(),
                    old_chambers[&id].last_window_duration_s()
                );
            }
        });
    });
}
