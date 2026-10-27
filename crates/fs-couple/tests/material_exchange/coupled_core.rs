//! Direct mixed-boundary tests of the shared circuit/enthalpy residual.
//! The existing synthetic cards are implementation fixtures, not validation data.
use super::*;
use fs_couple::cavity_phs::piston::reactive::implicit::coupled::ImplicitMechanismSolve;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::{
    self, Chamber, DrivenNetworkCandidate, DrivenNetworkError,
    reactive::{self, CoupledBoundaryError, CoupledChamberBoundary, CoupledChamberReaction},
};
use fs_couple::cavity_phs::piston::{GasPiston, reactive::GasReaction};
use fs_couple::circuit_enthalpy::CircuitEnthalpyInput;
use fs_couple::circuit_enthalpy::contacts::ContactThermalSample;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

pub(super) fn implicit_reaction(rate: f64) -> GasReaction {
    let law = reactive_tests::reaction(rate);
    let policy = ImplicitMechanismSolve {
        extent_brackets_mol: law
            .channels()
            .map(|(kinetics, _)| {
                (
                    kinetics
                        .reverse_rate_closure()
                        .equilibrium()
                        .reaction()
                        .clone(),
                    [-1.0, 1.0],
                )
            })
            .collect(),
        extent_tolerance_mol: 1e-12,
        interior_fraction: 0.05,
        jacobian_relative_step: 1e-5,
        maximum_iterations: 32,
        maximum_evaluations: 1024,
        newton: fs_solver::NewtonKrylovConfig {
            absolute_tolerance: 1e-13,
            relative_tolerance: 1e-12,
            linear_restart: 8,
            max_linear_cycles: 8,
            forcing_maximum: 1e-6,
            ..Default::default()
        },
    };
    law.with_implicit_mechanism(policy).unwrap()
}

pub(super) fn binding(law: &GasReaction) -> CoupledChamberReaction<'_> {
    CoupledChamberReaction {
        chamber: 10,
        reaction: law,
        chemistry: reactive_tests::chemical(),
        pressure_interval_pa: [50000.0, 300000.0],
        pressure_tolerance_pa: 1e-6,
    }
}

type Candidate<'t> = DrivenNetworkCandidate<
    't,
    ContactThermalSample<MaterialExchangeSample>,
    CoupledChamberBoundary<ResolvedMaterialStatePoint>,
>;
type Error = CircuitEnthalpyError<
    DrivenNetworkError<MaterialExchangeError, CoupledBoundaryError<MaterialNetworkError>>,
>;

#[allow(clippy::too_many_arguments)]
fn solve<'t>(
    cx: &Cx<'_>,
    assembly: &MaterialExchangeNetwork<'_>,
    prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
    thermal: &DepositionSimulation<'t>,
    target: f64,
    circuit: &MaterialCircuit,
    chambers: &[Chamber<'_>],
    reactions: &[CoupledChamberReaction<'_>],
    cancelled: impl FnMut() -> bool,
) -> Result<Candidate<'t>, Error> {
    let destinations: Vec<_> = assembly
        .electrical
        .bulk
        .resistors
        .iter()
        .map(|resistor| resistor.destination)
        .collect();
    let routes = assembly.electrical.heat_destinations(thermal).unwrap();
    let patches = gas_fixture::patches();
    let contacts = patches.map(|patch| patch.contact);
    reactive::step_with_contacts_candidate(
        thermal,
        target,
        CircuitEnthalpyInput {
            circuit,
            ports: &[],
            destinations: &destinations,
        },
        &routes,
        chambers,
        reactions,
        &contacts,
        controls(),
        &[0.0; 2],
        |state| prepared.resolve(cx, state),
        |contact, state, source, gas| {
            prepared.resolve_conductive_gas_contact(
                *patches
                    .iter()
                    .find(|patch| patch.contact.id == contact.id)
                    .unwrap(),
                state,
                source,
                gas,
            )
        },
        cancelled,
    )
}

fn amounts(piston: &GasPiston) -> [f64; 2] {
    let total = piston.gas().mass_kg() / 0.028;
    let model = piston.gas().model().thermodynamics();
    let mut result = [0.0; 2];
    for (species, fraction) in model
        .components()
        .iter()
        .zip(model.mole_composition().fractions())
    {
        result[match species.species().as_str() {
            "A" => 0,
            "B" => 1,
            _ => panic!("unexpected fixture species"),
        }] = total * fraction;
    }
    result
}

#[test]
fn mixed_reaction_and_ordinary_gas_share_each_finite_heat_and_electrical_balance() {
    fixture([0.8, 0.8], 150.0, |assembly, mut thermal, mut circuit| {
        with_cx(|cx, _| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let rate = 10000.0;
            let law = implicit_reaction(rate);
            let mut gases = [
                reactive_tests::piston(&law),
                gas_fixture::pistons()[1].clone(),
            ];
            for _ in 0..2 {
                let target = thermal.time_s() + 0.00005;
                let dt = target - thermal.time_s();
                // Deliberately reverse input order; physical IDs own the bindings.
                let chambers = [
                    gas_fixture::chamber(20, &gases[1], dt),
                    gas_fixture::chamber(10, &gases[0], dt),
                ];
                let result = solve(
                    cx,
                    assembly,
                    &prepared,
                    &thermal,
                    target,
                    &circuit,
                    &chambers,
                    &[binding(&law)],
                    || false,
                )
                .unwrap();
                assert_eq!(
                    result.chambers().iter().map(|c| c.id).collect::<Vec<_>>(),
                    [10, 20]
                );
                let reactive = &result.chambers()[0].boundary;
                let ordinary = &result.chambers()[1].boundary;
                let reaction = &reactive.reactive_step().unwrap().piston.reaction.reaction;
                assert!(reactive.nonreactive_step().is_none());
                assert!(ordinary.reactive_step().is_none());
                assert!(ordinary.nonreactive_step().is_some());
                let n0 = amounts(&gases[0]);
                let n1 = amounts(reactive.piston());
                let t1 = reactive.piston().gas().state().state().temperature;
                let xi = n1[1] - n0[1];
                let extent_from_rate =
                    dt * rate * R * t1 / 100000.0 * (n1[0] - n1[1] * (-50000.0 / (R * t1)).exp());
                assert!((xi - extent_from_rate).abs() < 2e-12);
                assert!((reaction.extent_mol - xi).abs() < 1e-14);
                assert!((n1[0] - n0[0] + xi).abs() < 1e-14);
                let chemical_report = reaction.implicit_mechanism.as_ref().unwrap();
                assert!(chemical_report.maximum_residual_mol <= 1e-12);
                assert!(
                    reactive.reactive_step().unwrap().chemical_evaluations
                        >= chemical_report.evaluations
                );

                let mut gas_heat = std::collections::BTreeMap::<CellId, f64>::new();
                let mut exchanged = 0.0;
                let mut gas_and_work = 0.0;
                let endpoint = &result.material().material;
                for (index, chamber) in result.chambers().iter().enumerate() {
                    let gas = chamber.boundary.piston();
                    let gas_temperature = if let Some(step) = chamber.boundary.reactive_step() {
                        let n = amounts(gas);
                        near(
                            gas.gas().mass_kg() * gas.gas().state().specific_internal_energy_j_kg(),
                            2.5 * R * (n[0] + n[1]) * gas.gas().state().state().temperature
                                - 50000.0 * n[1],
                        );
                        assert!(step.total_energy_residual_j.abs() <= 1e-7);
                        gas.gas().state().state().temperature
                    } else {
                        chamber
                            .boundary
                            .nonreactive_step()
                            .unwrap()
                            .piston
                            .midpoint
                            .state()
                            .temperature
                    };
                    for sample in chamber.boundary.source() {
                        let cell = sample.source.0.cell;
                        let temperature = result.thermal().temperature_k(cell).unwrap();
                        near(
                            sample.heat_j,
                            dt * sample.conductance_w_k * (temperature - gas_temperature),
                        );
                        let bulk = endpoint.bulk.cells.iter().find(|s| s.cell == cell).unwrap();
                        assert_eq!(bulk.conductivity.as_ref(), Some(&sample.source.1));
                        *gas_heat.entry(cell).or_default() += sample.heat_j;
                        exchanged += sample.heat_j;
                    }
                    assert_eq!(gas.gas().mass_kg(), gases[index].gas().mass_kg());
                    gas_and_work += gas.gas().mass_kg()
                        * (gas.gas().state().specific_internal_energy_j_kg()
                            - gases[index].gas().state().specific_internal_energy_j_kg())
                        + 0.5
                            * gas.mass_kg()
                            * (gas.velocity_m_s().powi(2) - gases[index].velocity_m_s().powi(2))
                        + 100000.0 * (gas.gas().volume_m3() - gases[index].gas().volume_m3());
                }
                near(exchanged, result.step().transferred_heat_j);
                let storage_change = result.thermal().energy_j() - thermal.energy_j()
                    + (result.circuit().charges_c()[0].powi(2) - circuit.charges_c()[0].powi(2))
                        / 0.02;
                assert!(
                    (storage_change + gas_and_work).abs() < controls().total_energy_tolerance_j
                );
                assert!(
                    result.step().total_energy_residual_j.abs()
                        < controls().total_energy_tolerance_j
                );
                near(result.step().external_network_heat_j, 0.0);
                let (next_thermal, next_circuit, report, _) = result.into_parts();
                // Independently reconstruct EACH receiving cell: bulk Joule + its
                // declared contact fraction + contact conduction + radiation - gas.
                let material = &report.source.material.material;
                let contact = &material.contacts[0];
                let link_heat = dt
                    * contact.conductance_w_k
                    * (next_thermal.temperature_k(contact.cells[1]).unwrap()
                        - next_thermal.temperature_k(contact.cells[0]).unwrap());
                for cell in [CellId(30), CellId(10)] {
                    let bulk_heat: f64 = report
                        .electrical
                        .resistors
                        .iter()
                        .filter(|loss| {
                            assembly.electrical.bulk.resistors.iter().any(|binding| {
                                binding.destination.cell == cell
                                    && binding.destination.branch == loss.material.branch
                            })
                        })
                        .map(|loss| loss.heat_j)
                        .sum();
                    let route = &report.source.material.destinations[0];
                    let fraction = if cell == route.cells[0] {
                        route.heat_fraction_to_first
                    } else {
                        1.0 - route.heat_fraction_to_first
                    };
                    let interface_heat = fraction * report.electrical.contacts[0].heat_j;
                    assert!(interface_heat > 0.0 && bulk_heat > 0.0);
                    let radiative = dt * material.radiation.as_ref().unwrap().sources_w()[&cell];
                    let conduction = if cell == contact.cells[0] {
                        link_heat
                    } else {
                        -link_heat
                    };
                    let actual = next_thermal.cells()[&cell].mass_kg
                        * (next_thermal.cells()[&cell].specific_h_j_kg
                            - thermal.cells()[&cell].specific_h_j_kg);
                    near(
                        actual,
                        bulk_heat + interface_heat + radiative + conduction - gas_heat[&cell],
                    );
                }
                for (index, chamber) in report.source.chambers.into_iter().enumerate() {
                    gases[index] = chamber.boundary.into_parts().0;
                }
                thermal = next_thermal;
                circuit = next_circuit;
            }
        })
    });
}

#[test]
fn absent_reactions_preserve_the_original_network_driver_exactly() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let gases = gas_fixture::pistons();
            let chambers = [
                gas_fixture::chamber(10, &gases[0], 0.00005),
                gas_fixture::chamber(20, &gases[1], 0.00005),
            ];
            let candidate = solve(
                cx,
                assembly,
                &prepared,
                &thermal,
                0.00005,
                &circuit,
                &chambers,
                &[],
                || false,
            )
            .unwrap();
            let destinations: Vec<_> = assembly
                .electrical
                .bulk
                .resistors
                .iter()
                .map(|r| r.destination)
                .collect();
            let routes = assembly.electrical.heat_destinations(&thermal).unwrap();
            let patches = gas_fixture::patches();
            let contacts = patches.map(|patch| patch.contact);
            let old = network::step_with_contacts_candidate(
                &thermal,
                0.00005,
                CircuitEnthalpyInput {
                    circuit: &circuit,
                    ports: &[],
                    destinations: &destinations,
                },
                &routes,
                &chambers,
                &contacts,
                controls(),
                &[0.0; 2],
                |state| prepared.resolve(cx, state),
                |contact, state, source, gas| {
                    prepared.resolve_conductive_gas_contact(
                        *patches
                            .iter()
                            .find(|patch| patch.contact.id == contact.id)
                            .unwrap(),
                        state,
                        source,
                        gas,
                    )
                },
                || false,
            )
            .unwrap();
            assert_eq!(candidate.thermal().cells(), old.thermal().cells());
            assert_eq!(candidate.circuit(), old.circuit());
            assert_eq!(candidate.step(), old.step());
            for (current, previous) in candidate.chambers().iter().zip(old.chambers()) {
                assert!(current.boundary.reactive_step().is_none());
                assert_eq!(current.boundary.piston(), previous.boundary.piston());
                assert_eq!(current.boundary.source(), previous.boundary.source());
                assert_eq!(
                    current.boundary.nonreactive_step().unwrap(),
                    previous.boundary.step()
                );
            }
        })
    });
}

#[test]
fn reaction_mapping_and_final_checkpoint_failure_preserve_every_accepted_history() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let law = implicit_reaction(10000.0);
            let gases = [
                reactive_tests::piston(&law),
                gas_fixture::pistons()[1].clone(),
            ];
            let original_gases = gases.clone();
            let original_cells = thermal.cells().clone();
            let original_circuit = circuit.clone();
            let chambers = [
                gas_fixture::chamber(10, &gases[0], 0.00005),
                gas_fixture::chamber(20, &gases[1], 0.00005),
            ];
            let source = binding(&law);
            let mut calls = 0;
            let accepted = solve(
                cx,
                assembly,
                &prepared,
                &thermal,
                0.00005,
                &circuit,
                &chambers,
                &[source],
                || {
                    calls += 1;
                    false
                },
            )
            .unwrap();
            assert!(
                accepted.chambers()[0]
                    .boundary
                    .reactive_step()
                    .unwrap()
                    .piston
                    .reaction
                    .reaction
                    .extent_mol
                    > 0.0
            );
            assert!(
                solve(
                    cx,
                    assembly,
                    &prepared,
                    &thermal,
                    0.00005,
                    &circuit,
                    &chambers,
                    &[source, source],
                    || false
                )
                .is_err()
            );
            assert!(
                solve(
                    cx,
                    assembly,
                    &prepared,
                    &thermal,
                    0.00005,
                    &circuit,
                    &chambers,
                    &[CoupledChamberReaction {
                        chamber: 99,
                        ..source
                    }],
                    || false
                )
                .is_err()
            );
            let mut observed = 0;
            assert!(matches!(
                solve(
                    cx,
                    assembly,
                    &prepared,
                    &thermal,
                    0.00005,
                    &circuit,
                    &chambers,
                    &[source],
                    || {
                        observed += 1;
                        observed == calls
                    }
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            assert_eq!(gases, original_gases);
            assert_eq!(thermal.cells(), &original_cells);
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(circuit, original_circuit);
        })
    });
}
