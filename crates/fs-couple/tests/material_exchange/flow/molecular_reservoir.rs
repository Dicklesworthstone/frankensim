//! Original held baths, molecular cross effects and complete live feedback.
#[path = "molecular_reservoir_gas_law.rs"]
mod gas_law_tests;

use super::*;
use fs_gas::material::transport::parcel::reservoir::{MolecularReservoir, ReservoirMolecularLink};

fn bath(bank: &ThermochemicalGasModel, x: f64, t: f64, p: f64) -> MolecularReservoir {
    MolecularReservoir::new(
        bank.with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![x, 1.0 - x]).unwrap(),
        )
        .unwrap(),
        gas_point(t, p),
        || false,
    )
    .unwrap()
}

fn bath_link(bath: &MolecularReservoir) -> ReservoirMolecularLink<'_> {
    ReservoirMolecularLink {
        id: 101,
        parcel: 10,
        reservoir: bath,
        area_m2: 0.1,
        distance_m: 0.0001,
    }
}

fn same_owners(actual: &MaterialExchangeState<'_>, expected: &MaterialExchangeState<'_>) {
    assert_eq!(actual.time_s(), expected.time_s());
    assert_eq!(actual.chambers(), expected.chambers());
    assert_eq!(actual.circuit(), expected.circuit());
    assert_eq!(actual.thermal().cells(), expected.thermal().cells());
    for (&id, chamber) in actual.chambers() {
        assert_eq!(
            chamber.last_window_transfers(),
            expected.chambers()[&id].last_window_transfers()
        );
        assert_eq!(
            chamber.last_window_duration_s(),
            expected.chambers()[&id].last_window_duration_s()
        );
    }
}

#[test]
fn cooling_bath_matches_binary_soret_dufour_temperature_and_open_entropy() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
            let held = bath(&bank, 0.5, 390.0, 1e5);
            let retained = held.clone();
            let link = bath_link(&held);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(0.5, 410.0, 1e5), (0.5, 400.0, 1e5)], true),
            );
            let (next, report) = state
                .molecular_reservoir_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives(state.chambers()),
                    &[],
                    &[],
                    &mut molecular,
                    &[],
                    &[link],
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            // The normal points from the held 390 K bath into the 410 K gas.
            // Equal molar masses remove pressure diffusion; this original
            // binary law has Cv=2.5 R and h*_A-h*_B=100000 J/kg at 400 K.
            let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
            let gradient_t = 20.0 / link.distance_m;
            let driving = 0.25 * MASS * 100000.0 / (r * 400.0) * gradient_t / 400.0;
            let j = -1e5 / (r * 400.0) * DIFFUSIVITY * MASS * driving;
            let q = -CONDUCTIVITY * gradient_t + 100000.0 * j;
            let area_time = link.area_m2 * DT;
            assert!(j < 0.0 && q < 0.0);
            let face = &report.flow.links[0];
            close(face.mass_into_gas_kg[0], j * area_time, 1e-16);
            close(face.mass_into_gas_kg[1], -j * area_time, 1e-16);
            close(face.reduced_heat_into_gas_j, q * area_time, 1e-12);
            close(face.face.enthalpy_temperature_k, 390.0, 1e-12);
            assert_eq!(&face.reservoir, held.state());
            assert_eq!(held, retained);
            let old = state.chambers()[&10].gas();
            let gas = next.chambers()[&10].gas();
            let amount = old.mass_kg() / MASS;
            let expected_t = 410.0 + q * area_time / (amount * 2.5 * r);
            close(gas.state().state().temperature, expected_t, 2e-8);
            close(gas.mass_kg(), old.mass_kg(), 1e-16);
            close(
                species_mass(gas, "A") - species_mass(old, "A"),
                j * area_time,
                1e-16,
            );
            // The equal-Cp, equal-composition bath species entropy streams
            // cancel. The gas may lose entropy while gas plus bath gains it.
            let x = species_mass(gas, "A") / gas.mass_kg();
            let mixing = x * x.ln() + (1.0 - x) * (1.0 - x).ln() + 2.0_f64.ln();
            let gas_entropy = amount * r * (2.5 * (expected_t / 410.0).ln() - mixing);
            assert!(gas_entropy < 0.0);
            let production = gas_entropy - q * area_time / 390.0;
            assert!(production > 0.0);
            close(report.flow.audit.entropy_production_j_k, production, 2e-11);
            assert!(
                report.balance.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j
            );
            assert!(
                report.balance.maximum_species_mass_residual_kg
                    <= molecular_controls().physical.exchange.mass_tolerance_kg
            );
            assert!(matches!(
                next.chambers()[&10].audit_last_window().unwrap(),
                ModuleWindowAudit::Unavailable { .. }
            ));
            assert_ne!(next.circuit(), state.circuit());
            assert_eq!(state.time_s(), 0.0);
        });
    });
}

#[test]
fn unequal_mass_bath_and_prescribed_heat_are_invariant_to_large_formation_datums() {
    let run = |formation_b| {
        fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
            with_cx(|cx, _| {
                let (bank, mut molecular) = molecular_fixture(2.0 * MASS, 1.0, formation_b);
                let held = bath(&bank, 0.7, 420.0, 1e5);
                let link = bath_link(&held);
                let prepared = assembly.prepare(cx, &thermal).unwrap();
                let state = MaterialExchangeState::from_parts(
                    thermal,
                    circuit,
                    molecular_gases(&bank, [(0.3, 380.0, 1e5), (0.5, 400.0, 1e5)], true),
                );
                let mut loads = drives(state.chambers());
                loads[0].heat_j = 0.015;
                let (next, report) = state
                    .molecular_reservoir_step_candidate(
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
                    )
                    .unwrap();
                let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
                let gradient_t = -40.0 / link.distance_m;
                let gradient_x = -0.4 / link.distance_m;
                let centered_a = (2.0 / 3.0) * 100000.0;
                let driving =
                    gradient_x + 0.5 * MASS * centered_a / (r * 400.0 * 400.0) * gradient_t;
                let j =
                    -1e5 / (r * 400.0) * DIFFUSIVITY * (MASS * 2.0 * MASS / (1.5 * MASS)) * driving;
                let q = -CONDUCTIVITY * gradient_t + 100000.0 * j;
                let area_time = link.area_m2 * DT;
                let dm = j * area_time;
                let dn = dm / MASS - dm / (2.0 * MASS);
                let n0 = 1e5 * 0.001 / (r * 380.0);
                // Both signs of species transfer carry the 420 K bath Cp T.
                // Each species' formation contribution cancels analytically
                // against its actual changed inventory, before numeric heat.
                let expected_t =
                    (2.5 * r * n0 * 380.0 + 3.5 * r * 420.0 * dn + q * area_time + loads[0].heat_j)
                        / (2.5 * r * (n0 + dn));
                let gas = next.chambers()[&10].gas();
                close(gas.state().state().temperature, expected_t, 2e-8);
                let face = &report.flow.links[0];
                close(face.mass_into_gas_kg[0], dm, 1e-16);
                close(face.mass_into_gas_kg[1], -dm, 1e-16);
                close(face.reduced_heat_into_gas_j, q * area_time, 1e-12);
                close(face.face.enthalpy_temperature_k, 420.0, 1e-12);
                assert!(report.flow.audit.entropy_production_j_k > 0.0);
                assert!(
                    report.balance.total_energy_residual_j.abs()
                        <= controls().total_energy_tolerance_j
                );
                assert!(
                    report.balance.maximum_chamber_energy_residual_j
                        <= controls().total_energy_tolerance_j
                );
                close(
                    report.balance.total_energy_residual_j,
                    report.flow.audit.energy_residual_j
                        + report.exchange.step.total_energy_residual_j,
                    1e-14,
                );
                let account = report
                    .balance
                    .chambers
                    .iter()
                    .find(|account| account.chamber == 10)
                    .unwrap();
                close(account.feedback_energy_change_j, loads[0].heat_j, 2e-8);
                assert_eq!(
                    account.transport_energy_change_j,
                    report
                        .flow
                        .exchange
                        .parcels
                        .iter()
                        .find(|parcel| parcel.id == 10)
                        .unwrap()
                        .energy_change_j
                );
                if formation_b > 1e10 {
                    assert!(account.transport_energy_change_j.abs() > 1e6);
                }
                (
                    gas.state().state().temperature,
                    face.mass_into_gas_kg.clone(),
                    face.reduced_heat_into_gas_j,
                    report.balance.total_energy_residual_j,
                )
            })
        })
    };
    let original = run(0.0);
    let shifted = run(1e14);
    close(original.0, shifted.0, 2e-8);
    for (a, b) in original.1.iter().zip(&shifted.1) {
        close(*a, *b, 1e-16);
    }
    close(original.2, shifted.2, 1e-12);
    close(original.3, shifted.3, 2e-8);
}

#[test]
fn internal_and_two_bath_faces_match_native_transport_then_reacting_moving_wall_feedback() {
    fixture([0.2, 0.8], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let formation_b = -50000.0;
            let (bank, mut native) = molecular_fixture(MASS, 1.0, formation_b);
            let (_, mut molecular) = molecular_fixture(MASS, 1.0, formation_b);
            let baths = [
                bath(&bank, 0.9, 430.0, 1.2e5),
                bath(&bank, 0.1, 370.0, 0.9e5),
            ];
            let retained = baths.clone();
            let links = [
                ReservoirMolecularLink {
                    area_m2: 0.08,
                    ..bath_link(&baths[0])
                },
                ReservoirMolecularLink {
                    id: 202,
                    parcel: 20,
                    area_m2: 0.06,
                    ..bath_link(&baths[1])
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
            loads[0].heat_j = 0.01;
            loads[1].heat_j = -0.005;
            loads[0].external_force_n = -3.0;
            loads[1].external_force_n = 2.0;
            let law = reaction(50.0);
            let reactions = [ChamberReaction {
                chamber: 10,
                reaction: &law,
                controls: chemical(),
            }];
            let parcels = state
                .chambers()
                .iter()
                .map(|(&id, chamber)| (id, chamber.gas().clone()))
                .collect();
            let endpoints = loads
                .iter()
                .map(|drive| (drive.chamber, drive.piston_solve.gas))
                .collect();
            let (transported, native_flow) = native
                .network_candidate(
                    &parcels,
                    &internal,
                    &links,
                    DT,
                    &endpoints,
                    molecular_controls(),
                    || false,
                )
                .unwrap()
                .into_parts();
            let staged = MaterialExchangeState::from_parts(
                state.thermal().clone(),
                state.circuit().clone(),
                transported
                    .into_iter()
                    .map(|(id, parcel)| (id, GasPiston::new(parcel, 1.0, 0.01, 0.0).unwrap()))
                    .collect(),
            );
            let patches = gas_fixture::patches();
            let (expected, original_feedback) = staged
                .step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &reactions,
                    &patches,
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let (actual, report) = state
                .molecular_reservoir_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &reactions,
                    &patches,
                    &mut molecular,
                    &internal,
                    &links,
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            assert_eq!(report.flow, native_flow);
            assert_eq!(report.exchange.step, original_feedback.step);
            assert_eq!(
                report.exchange.report.source.reactions,
                original_feedback.report.source.reactions
            );
            let extent = report.exchange.report.source.reactions[0]
                .reaction
                .extent_mol;
            assert!(extent > 0.0);
            assert_eq!(actual.circuit(), expected.circuit());
            assert_eq!(actual.thermal().cells(), expected.thermal().cells());
            assert_eq!(baths, retained);
            for (&id, chamber) in actual.chambers() {
                let oracle = &expected.chambers()[&id];
                assert_eq!(chamber.gas(), oracle.gas());
                assert_eq!(chamber.position_m(), oracle.position_m());
                assert_eq!(chamber.velocity_m_s(), oracle.velocity_m_s());
                assert!(chamber.velocity_m_s().abs() > 1e-8);
                assert_eq!(chamber.last_window_duration_s(), Some(DT));
                assert!(matches!(
                    chamber.audit_last_window().unwrap(),
                    ModuleWindowAudit::Unavailable { .. }
                ));
            }
            assert!(
                report
                    .exchange
                    .report
                    .source
                    .chambers
                    .iter()
                    .any(|chamber| chamber.boundary.piston.gas.supplied_heat_j.abs() > 1e-6)
            );
            // Independent absolute-energy cross-check only at a moderate
            // formation reference, where direct endpoint subtraction is safe.
            // The public admission keeps its reference-cancelled native audit.
            let mut residual = actual.thermal().energy_j() - state.thermal().energy_j()
                + report.exchange.report.electrical.stored_after_j
                - report.exchange.report.electrical.stored_before_j
                - report.exchange.report.electrical.native.supplied
                - report.exchange.step.external_network_heat_j;
            for boundary in &report.exchange.report.source.chambers {
                let id = boundary.id;
                let before = &state.chambers()[&id];
                let after = &actual.chambers()[&id];
                residual += after.gas().mass_kg()
                    * after.gas().state().specific_internal_energy_j_kg()
                    - before.gas().mass_kg() * before.gas().state().specific_internal_energy_j_kg()
                    + 0.5
                        * before.mass_kg()
                        * (after.velocity_m_s().powi(2) - before.velocity_m_s().powi(2))
                    + boundary.boundary.piston.exterior_pressure_work_j
                    - boundary.boundary.piston.external_force_work_j
                    - loads.iter().find(|load| load.chamber == id).unwrap().heat_j;
            }
            let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
            let supplied: f64 = report
                .flow
                .links
                .iter()
                .map(|face| {
                    let t = face.reservoir.state().temperature;
                    face.mass_into_gas_kg[0] / MASS * (3.5 * r * t)
                        + face.mass_into_gas_kg[1] / MASS * (3.5 * r * t + formation_b)
                        + face.reduced_heat_into_gas_j
                })
                .sum();
            residual -= supplied;
            close(residual, 0.0, controls().total_energy_tolerance_j);
            close(residual, report.balance.total_energy_residual_j, 2e-8);
            assert!(
                report.balance.maximum_species_mass_residual_kg
                    <= molecular_controls().physical.exchange.mass_tolerance_kg
            );
            assert!(
                report.balance.maximum_chamber_energy_residual_j
                    <= controls().total_energy_tolerance_j
            );
            for row in &report.balance.species {
                let index = report
                    .flow
                    .species
                    .iter()
                    .position(|species| species == &row.species)
                    .unwrap();
                let internal_change: f64 = report
                    .flow
                    .internal_links
                    .iter()
                    .map(|link| {
                        if link.parcels[0] == row.chamber {
                            -link.mass_left_to_right_kg[index]
                        } else if link.parcels[1] == row.chamber {
                            link.mass_left_to_right_kg[index]
                        } else {
                            0.0
                        }
                    })
                    .sum();
                let reservoir_change: f64 = report
                    .flow
                    .links
                    .iter()
                    .filter(|link| link.parcel == row.chamber)
                    .map(|link| link.mass_into_gas_kg[index])
                    .sum();
                let reaction_change = if row.chamber == 10 {
                    if row.species == species_id("A") {
                        -extent * MASS
                    } else {
                        extent * MASS
                    }
                } else {
                    0.0
                };
                let name = if row.species == species_id("A") {
                    "A"
                } else {
                    "B"
                };
                close(
                    row.initial_mass_kg,
                    species_mass(state.chambers()[&row.chamber].gas(), name),
                    1e-16,
                );
                close(
                    row.final_mass_kg,
                    species_mass(actual.chambers()[&row.chamber].gas(), name),
                    1e-16,
                );
                close(row.internal_change_kg, internal_change, 1e-16);
                close(row.reservoir_change_kg, reservoir_change, 1e-16);
                close(row.reaction_change_kg, reaction_change, 1e-16);
                assert!(
                    row.residual_kg.abs()
                        <= molecular_controls().physical.exchange.mass_tolerance_kg
                );
            }
            let bath_mass: f64 = report
                .flow
                .links
                .iter()
                .flat_map(|link| &link.mass_into_gas_kg)
                .sum();
            close(report.balance.net_reservoir_mass_kg, bath_mass, 1e-16);
            let gas_mass_change: f64 = actual
                .chambers()
                .iter()
                .map(|(&id, gas)| gas.gas().mass_kg() - state.chambers()[&id].gas().mass_kg())
                .sum();
            close(
                report.balance.total_mass_residual_kg,
                gas_mass_change - bath_mass,
                1e-16,
            );
            assert_eq!(state.time_s(), 0.0);
        });
    });
}

#[test]
fn pure_bath_heat_and_late_feedback_refusal_preserve_complete_state_and_ports() {
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
            let capacity =
                state.chambers()[&10].gas().mass_kg() * 2.5 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
                    / MASS;
            let expected_heat = CONDUCTIVITY * 40.0 / link.distance_m * link.area_m2 * DT;
            let (preview, report) = state
                .molecular_reservoir_step_candidate(
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
                )
                .unwrap();
            assert!(report.flow.exchange.boundaries.is_empty());
            assert_eq!(report.flow.links[0].mass_into_gas_kg, vec![0.0; 2]);
            close(
                report.flow.links[0].reduced_heat_into_gas_j,
                expected_heat,
                1e-12,
            );
            let gas = preview.chambers()[&10].gas();
            close(
                gas.state().state().temperature,
                400.0 + expected_heat / capacity,
                2e-8,
            );
            assert_eq!(gas.model().thermodynamics().components().len(), 1);
            assert_eq!(species_mass(gas, "B"), 0.0);
            assert!(matches!(
                preview.chambers()[&10].audit_last_window().unwrap(),
                ModuleWindowAudit::Unavailable { .. }
            ));
            state
                .advance_with_molecular_reservoir(
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
                )
                .unwrap();
            same_owners(&state, &preview);
            let mut limited = controls();
            limited.network.maximum_circuit_trials = 1;
            assert!(matches!(
                state.advance_with_molecular_reservoir(
                    &prepared,
                    cx,
                    2.0 * DT,
                    &[],
                    &loads,
                    &[],
                    &gas_fixture::patches(),
                    &mut molecular,
                    &[],
                    &[link],
                    molecular_controls(),
                    limited,
                    &[0.0; 2],
                ),
                Err(MaterialFlowError::Exchange(
                    ReactiveExchangeError::Exchange(CircuitEnthalpyError::CircuitBudget)
                ))
            ));
            same_owners(&state, &preview);
            gate.request();
            assert!(matches!(
                state.advance_with_molecular_reservoir(
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
                ),
                Err(MaterialFlowError::Cancelled)
            ));
            same_owners(&state, &preview);
        });
    });
}
