//! Finite circuit/cells and sealed reacting chambers share one refined interval.
use super::*;
use fs_couple::PortedModule;
use fs_couple::cavity_phs::piston::reactive::material_rate::{
    ACTIVITY_RATE_SCALE_DIMS, ACTIVITY_RATE_SCALE_PROPERTY, MaterialRateScale,
};
use fs_couple::cavity_phs::piston::reactive::thermal::subcycling::ReactiveThermalSubcycling;
use fs_couple::cavity_phs::piston::thermal::ThermalPistonError;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::NetworkError;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::exchange::gas_law::GasSurface;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::state::{
    ChamberDrive, MaterialExchangeState, subcycling::*,
};
use std::cell::Cell;
use std::collections::BTreeMap;

fn two_steps() -> MaterialExchangeSubcycling {
    MaterialExchangeSubcycling {
        coupled: ReactiveThermalSubcycling {
            absolute_amount_mol: 1.0,
            absolute_temperature_k: 1000.0,
            absolute_pressure_pa: 1e6,
            absolute_volume_m3: 1.0,
            absolute_velocity_m_s: 100.0,
            absolute_body_enthalpy_j_kg: 1e6,
            relative_tolerance: 0.0,
            maximum_substeps: 2,
            maximum_evaluations: 3,
        },
        absolute_charge_c: 1.0,
        absolute_flux_wb: 1.0,
    }
}
fn drive(chamber: u64, gas: &GasPiston) -> ChamberDrive {
    let c = gas_fixture::chamber(chamber, gas, 0.00005);
    ChamberDrive {
        chamber,
        heat_j: 0.0,
        exterior_pressure_pa: c.input.exterior_pressure_pa,
        external_force_n: 0.0,
        piston_solve: c.piston_solve,
        exchange: c.exchange,
    }
}
fn gases(law: &GasReaction, fixed: bool) -> BTreeMap<u64, GasPiston> {
    [(10, piston(law)), (20, gas_fixture::pistons()[1].clone())]
        .into_iter()
        .map(|(id, gas)| {
            (
                id,
                if fixed {
                    gas.with_fixed_support().unwrap()
                } else {
                    gas
                },
            )
        })
        .collect()
}
fn sourced_reaction() -> GasReaction {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(ACTIVITY_RATE_SCALE_PROPERTY, ACTIVITY_RATE_SCALE_DIMS),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, 25.0), (400.0, 50.0), (1800.0, 400.0)],
                dims: ACTIVITY_RATE_SCALE_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "temperature",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                200.0,
                1800.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic coupled kinetic scale".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let material = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-rate".into(),
            phase: "gas".into(),
            process: "manufactured-temperature-law".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    reaction(50.0)
        .with_material_rate(
            MaterialRateScale::new(
                material,
                MaterialPropertySelection::SingleClaimOnly,
                fs_qty::Pressure::new(100000.0),
            )
            .unwrap(),
        )
        .unwrap()
}
fn tight(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b} within {tolerance}");
}
fn physical_energy(state: &MaterialExchangeState<'_>) -> f64 {
    // Fixture's one physical 0.01-F capacitor, independently of solver reports.
    state.thermal().energy_j()
        + 0.5 * state.circuit().charges_c()[0].powi(2) / 0.01
        + state
            .chambers()
            .values()
            .map(|g| {
                g.gas().mass_kg() * g.gas().state().specific_internal_energy_j_kg()
                    + 0.5 * g.mass_kg() * g.velocity_m_s().powi(2)
            })
            .sum::<f64>()
}

#[test]
fn material_radiation_changes_reaction_within_the_same_finite_network_interval() {
    let run = |emissivity| {
        fixture(
            [emissivity, emissivity],
            150.0,
            |assembly, thermal, circuit| {
                with_cx(|cx, _| {
                    let prepared = assembly.prepare(cx, &thermal).unwrap();
                    let law = sourced_reaction();
                    let state =
                        MaterialExchangeState::from_parts(thermal, circuit, gases(&law, true));
                    let initial = physical_energy(&state);
                    let drives = state
                        .chambers()
                        .iter()
                        .map(|(&id, gas)| drive(id, gas))
                        .collect::<Vec<_>>();
                    let bindings = [ChamberReaction {
                        chamber: 10,
                        reaction: &law,
                        controls: chemical(),
                    }];
                    let (coarse, original) = state
                        .step_candidate(
                            &prepared,
                            cx,
                            0.002,
                            &[],
                            &drives,
                            &bindings,
                            &gas_fixture::patches(),
                            controls(),
                            &[0.0; 2],
                        )
                        .unwrap();
                    let (next, report) = state
                        .step_subcycled_candidate(
                            &prepared,
                            cx,
                            0.002,
                            &[],
                            &drives,
                            &bindings,
                            &gas_fixture::patches(),
                            controls(),
                            &[0.0; 2],
                            two_steps(),
                        )
                        .unwrap();
                    assert_eq!(report.steps.len(), 2);
                    assert_eq!(report.evaluations, 3);
                    assert_eq!(report.chambers.len(), 2);
                    assert!(
                        report.maximum_cell_residual_bound_j
                            <= controls().network.cell_energy_tolerance_j
                    );
                    assert!(
                        report.network_energy_residual_j.abs()
                            <= controls().network.total_energy_tolerance_j
                    );
                    tight(
                        physical_energy(&next),
                        initial,
                        controls().total_energy_tolerance_j,
                    );
                    for interval in &report.steps {
                        assert!(
                            interval
                                .advance
                                .report
                                .source
                                .material
                                .material
                                .radiation
                                .is_some()
                        );
                        assert_eq!(
                            interval
                                .advance
                                .report
                                .source
                                .material
                                .material
                                .contacts
                                .len(),
                            1
                        );
                        assert!(interval.advance.report.electrical.contacts[0].heat_j > 0.0);
                    }
                    let scale = |i: usize| {
                        report.steps[i].advance.report.source.reactions[0]
                            .reaction
                            .rate_material
                            .as_ref()
                            .unwrap()
                            .property(ACTIVITY_RATE_SCALE_PROPERTY)
                            .unwrap()
                            .value_si()
                    };
                    assert!(scale(1) > scale(0));
                    for account in &report.chambers {
                        let gas = &next.chambers()[&account.chamber];
                        assert_eq!(gas.last_window_duration_s(), Some(0.002));
                        tight(
                            gas.last_window_transfers().unwrap()[1]
                                .energy_into()
                                .value(),
                            account.gas_boundary_heat_j,
                            1e-12,
                        );
                        assert_eq!(
                            gas.gas().mass_kg(),
                            state.chambers()[&account.chamber].gas().mass_kg()
                        );
                    }
                    assert_eq!(state.time_s(), 0.0);
                    let _ = coarse;
                    (
                        original.report.source.reactions[0].reaction.extent_mol,
                        report.chambers[0].reaction_extents_mol[0].1,
                        next.circuit().charges_c()[0],
                        next.thermal().temperature_k(CellId(30)).unwrap(),
                    )
                })
            },
        )
    };
    let low = run(0.2);
    let high = run(0.8);
    assert_eq!(low.0, high.0); // One original frozen reaction cannot see emissivity.
    assert!((low.1 - high.1).abs() > 1e-12);
    assert!((low.2 - high.2).abs() > 1e-9);
    assert!((low.3 - high.3).abs() > 1e-3);
}

#[test]
fn integrated_gas_heat_is_split_but_current_and_force_are_held_on_real_steps() {
    fixture([0.8, 0.8], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut graph = circuit.graph().clone();
            graph
                .branches
                .push((1, 0, Branch::CurrentSource { port: 0 }));
            let input_a = 0.25;
            let circuit = MaterialCircuit::new_with_contacts(
                graph,
                circuit.resistors().to_vec(),
                circuit.contacts().to_vec(),
                circuit.fluxes_wb().to_vec(),
                circuit.charges_c().to_vec(),
                &[input_a],
            )
            .unwrap();
            let law = sourced_reaction();
            let mut state = MaterialExchangeState::from_parts(thermal, circuit, gases(&law, false));
            let mut drives = state
                .chambers()
                .iter()
                .map(|(&id, gas)| drive(id, gas))
                .collect::<Vec<_>>();
            drives[0].heat_j = 0.03125;
            drives[1].heat_j = -0.015625;
            drives[0].external_force_n = -100.0;
            drives[1].external_force_n = 50.0;
            let reactions = [ChamberReaction {
                chamber: 10,
                reaction: &law,
                controls: chemical(),
            }];
            for _ in 0..2 {
                let energy = physical_energy(&state);
                let before = state.chambers().clone();
                let start = state.time_s();
                let target = start + 0.00005;
                let report = state
                    .advance_subcycled(
                        &prepared,
                        cx,
                        target,
                        &[input_a],
                        &drives,
                        &reactions,
                        &gas_fixture::patches(),
                        controls(),
                        &[0.0; 2],
                        two_steps(),
                    )
                    .unwrap();
                let mut supply = 0.0;
                for step in &report.steps {
                    let dt = step.end_s - step.start_s;
                    let circuit = &step.advance.report.electrical.native;
                    tight(circuit.supplied, dt * input_a * circuit.y[0], 1e-10);
                    supply += circuit.supplied;
                    for chamber in &step.advance.report.source.chambers {
                        let drive = drives.iter().find(|d| d.chamber == chamber.id).unwrap();
                        tight(
                            chamber.boundary.piston.gas.supplied_heat_j
                                - chamber.boundary.exchanged_heat_j,
                            drive.heat_j * dt / (target - start),
                            1e-12,
                        );
                    }
                }
                let mut mechanical = 0.0;
                for drive in &drives {
                    let old = &before[&drive.chamber];
                    let next = &state.chambers()[&drive.chamber];
                    let dv = next.gas().volume_m3() - old.gas().volume_m3();
                    mechanical += drive.external_force_n * dv / old.area_m2()
                        - drive.exterior_pressure_pa * dv;
                    let account = report
                        .chambers
                        .iter()
                        .find(|c| c.chamber == drive.chamber)
                        .unwrap();
                    let impulse: f64 = report
                        .steps
                        .iter()
                        .map(|s| {
                            let p = &s
                                .advance
                                .report
                                .source
                                .chambers
                                .iter()
                                .find(|c| c.id == drive.chamber)
                                .unwrap()
                                .boundary
                                .piston;
                            (s.end_s - s.start_s)
                                * (old.area_m2()
                                    * (p.midpoint.state().pressure - drive.exterior_pressure_pa)
                                    + drive.external_force_n)
                                + p.support_impulse_n_s
                        })
                        .sum();
                    let residual =
                        old.mass_kg() * (next.velocity_m_s() - old.velocity_m_s()) - impulse;
                    tight(residual, account.momentum_residual_n_s, 1e-12);
                    assert!(residual.abs() <= drive.piston_solve.momentum_tolerance_n_s);
                    assert_eq!(next.last_window_duration_s(), Some(target - start));
                }
                tight(report.prescribed_gas_heat_j, 0.015625, 1e-15);
                tight(
                    physical_energy(&state) - energy,
                    supply + mechanical + 0.015625,
                    controls().total_energy_tolerance_j,
                );
                assert_eq!(state.time_s(), target);
            }
        })
    });
}

#[test]
fn charge_tolerance_refines_the_real_circuit_without_mixed_unit_state_norms() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(thermal, circuit, BTreeMap::new());
            let (coarse, _) = state
                .step_candidate(
                    &prepared,
                    cx,
                    0.002,
                    &[],
                    &[],
                    &[],
                    &[],
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let (ordinary, report) = state
                .step_subcycled_candidate(
                    &prepared,
                    cx,
                    0.002,
                    &[],
                    &[],
                    &[],
                    &[],
                    controls(),
                    &[0.0; 2],
                    two_steps(),
                )
                .unwrap();
            let discrepancy =
                (coarse.circuit().charges_c()[0] - ordinary.circuit().charges_c()[0]).abs();
            assert!(discrepancy > 1e-10);
            let mut policy = two_steps();
            policy.absolute_charge_c = discrepancy / 4.0;
            policy.coupled.maximum_substeps = 64;
            policy.coupled.maximum_evaluations = 127;
            let (next, refined) = state
                .step_subcycled_candidate(
                    &prepared,
                    cx,
                    0.002,
                    &[],
                    &[],
                    &[],
                    &[],
                    controls(),
                    &[0.0; 2],
                    policy,
                )
                .unwrap();
            assert!(refined.steps.len() > report.steps.len());
            assert!(refined.estimated_error_ratio <= 1.0);
            assert!(refined.evaluations <= policy.coupled.maximum_evaluations);
            tight(
                physical_energy(&next),
                physical_energy(&state),
                controls().total_energy_tolerance_j,
            );
            assert_eq!(state.time_s(), 0.0);
        })
    });
}

#[test]
fn late_material_refusal_work_limit_and_cancellation_preserve_the_complete_assembly() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, gate| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let law = sourced_reaction();
            let mut state = MaterialExchangeState::from_parts(thermal, circuit, gases(&law, true));
            let cells = state.thermal().cells().clone();
            let circuit = state.circuit().clone();
            let gases = state.chambers().clone();
            let drives = state
                .chambers()
                .iter()
                .map(|(&id, gas)| drive(id, gas))
                .collect::<Vec<_>>();
            let reactions = [ChamberReaction {
                chamber: 10,
                reaction: &law,
                controls: chemical(),
            }];
            let patches = gas_fixture::patches();
            let surfaces = patches
                .iter()
                .map(|p| GasSurface {
                    contact: p.contact,
                    area_m2: p.geometry.area_m2,
                })
                .collect::<Vec<_>>();
            let count = Cell::new(0usize);
            let (_, report) = state
                .step_subcycled_with_gas_law_candidate(
                    &prepared,
                    cx,
                    0.00005,
                    &[],
                    &drives,
                    &reactions,
                    &surfaces,
                    controls(),
                    &[0.0; 2],
                    two_steps(),
                    |contact, state, source, gas| {
                        count.set(count.get() + 1);
                        let patch = *patches.iter().find(|p| p.contact == contact).unwrap();
                        prepared.resolve_conductive_gas_contact(patch, state, source, gas)
                    },
                )
                .unwrap();
            assert_eq!(report.steps.len(), 2);
            let last = count.get();
            assert!(last > 2);
            count.set(0);
            let failed = state.advance_subcycled_with_gas_law(
                &prepared,
                cx,
                0.00005,
                &[],
                &drives,
                &reactions,
                &surfaces,
                controls(),
                &[0.0; 2],
                two_steps(),
                |contact, state, source, gas| {
                    count.set(count.get() + 1);
                    if count.get() == last {
                        return Err("late material refusal");
                    }
                    let patch = *patches.iter().find(|p| p.contact == contact).unwrap();
                    prepared
                        .resolve_conductive_gas_contact(patch, state, source, gas)
                        .map_err(|_| "unexpected source refusal")
                },
            );
            assert!(matches!(
                failed,
                Err(MaterialExchangeSubcyclingError::Step(
                    ReactiveExchangeError::Exchange(CircuitEnthalpyError::Source(
                        NetworkError::Boundary {
                            error: ThermalPistonError::Interface("late material refusal"),
                            ..
                        }
                    ))
                ))
            ));
            let mut limited = two_steps();
            limited.coupled.maximum_evaluations = 2;
            assert!(matches!(
                state.advance_subcycled(
                    &prepared,
                    cx,
                    0.00005,
                    &[],
                    &drives,
                    &reactions,
                    &patches,
                    controls(),
                    &[0.0; 2],
                    limited
                ),
                Err(MaterialExchangeSubcyclingError::Limit { evaluations: 2, .. })
            ));
            count.set(0);
            let cancelled = state.advance_subcycled_with_gas_law(
                &prepared,
                cx,
                0.00005,
                &[],
                &drives,
                &reactions,
                &surfaces,
                controls(),
                &[0.0; 2],
                two_steps(),
                |contact, state, source, gas| {
                    count.set(count.get() + 1);
                    if count.get() == last {
                        gate.request();
                    }
                    let patch = *patches.iter().find(|p| p.contact == contact).unwrap();
                    prepared.resolve_conductive_gas_contact(patch, state, source, gas)
                },
            );
            assert!(matches!(
                cancelled,
                Err(MaterialExchangeSubcyclingError::Cancelled)
            ));
            assert_eq!(state.thermal().cells(), &cells);
            assert_eq!(state.circuit(), &circuit);
            assert_eq!(state.chambers(), &gases);
            assert_eq!(state.time_s(), 0.0);
            for (id, gas) in state.chambers() {
                assert_eq!(
                    gas.last_window_transfers(),
                    gases[id].last_window_transfers()
                );
                assert_eq!(
                    gas.last_window_duration_s(),
                    gases[id].last_window_duration_s()
                );
            }
        })
    });
}
