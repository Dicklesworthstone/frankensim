//! Actual electrode/RLC heat, original transport sources and finite receivers.
//! The independent two-temperature reference uses only the integrated Whitney
//! basis and explicit Fourier/gray-surface formulas, never the native solvers.

use super::super::nonlinear::heat::with_heat_material;
use super::*;
use fs_couple::circuit_enthalpy::{
    CircuitThermalCoefficients,
    contacts::{
        ContactThermalCoefficients,
        heat::{
            HeatAssembly, HeatNetworkError, HeatNetworkSample, MaterialHeatNetwork,
            PreparedHeatNetwork,
        },
        material::ElectricalContactNetwork,
    },
    material::{CellResistor, MaterialNetwork},
    spatial_induction::circuit::{
        SpatialCircuitThermalCandidate, SpatialTerminalCircuitThermalError, heat_assembly_candidate,
    },
};
use fs_exec::Cx;
use fs_matdb::ClaimSelection;
use fs_material::{
    field_storage::permanent::FLUX_DENSITY_DIMS,
    state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS},
    thermal::conduction::{
        LayerConductionConfig,
        network::{LayerLinkBinding, PreparedLayerConduction},
        source::ResolvedConductionLayer,
    },
};
use std::collections::BTreeMap;

const DT: f64 = 0.05;
const RESISTANCE: f64 = 0.5;
const MU_SLOPE: f64 = 0.0004;
const BASIS: ElectricalPropertyBasis = ElectricalPropertyBasis::Resistivity;

fn radiated(ta: f64, tb: f64, epsilon: f64) -> f64 {
    5.670_374_419e-8 * (tb.powi(4) - ta.powi(4))
        / (1.0 / (epsilon * (1.0 + 0.0001 * (ta - 300.0)))
            + 1.0 / (0.6 * (1.0 + 0.0001 * (tb - 300.0)))
            - 1.0)
}

struct Reference {
    temperature_k: [f64; 2],
    a: f64,
    voltage: f64,
    current: f64,
    field_heat: f64,
    graph_heat: f64,
    reversible_heat: f64,
}

fn reference(
    old_t: [f64; 2],
    old_a: f64,
    source_v: f64,
    transfer: impl Fn(f64, f64) -> f64,
) -> Reference {
    // int W_x=1/6, int W²=1/5, K=4/mu. Eliminate the actual electrode
    // current, circuit voltage and single free Whitney coordinate first.
    let mu_old = 1.0 + MU_SLOPE * (old_t[0] - 300.0);
    let at = |t: f64| {
        let mu = 1.0 + MU_SLOPE * (t - 300.0);
        let s = sigma(t, 0.2, 0.009, BASIS);
        let f = s / 6.0;
        let d = s / 5.0 + 2.0 * DT / mu;
        let magnetic = 2.0 * (mu_old.recip() + mu.recip()) * old_a;
        let y = s - f * f / d;
        let i0 = f * magnetic / d;
        let voltage = (source_v - RESISTANCE * i0) / (1.0 + RESISTANCE * y);
        let current = i0 + y * voltage;
        let rate = (f * voltage - magnetic) / d;
        let a = old_a + DT * rate;
        let joule = DT * s * (voltage * voltage - voltage * rate / 3.0 + rate * rate / 5.0);
        let parameter = 2.0 * old_a * a * (mu.recip() - mu_old.recip());
        let ts =
            2.0 * MU_SLOPE * (t * a * a / (mu * mu) - old_t[0] * old_a * old_a / (mu_old * mu_old));
        (a, voltage, current, joule, -parameter - ts)
    };
    let sink = |ta: f64| {
        let current = at(ta).2;
        bracket(280.0, 350.0, |tb| {
            4.0 * (tb - old_t[1]) - DT * RESISTANCE * current * current + DT * transfer(ta, tb)
        })
    };
    let ta = bracket(280.0, 350.0, |ta| {
        let (_, _, _, joule, reversible) = at(ta);
        4.0 * (ta - old_t[0]) - joule - reversible - DT * transfer(ta, sink(ta))
    });
    let (a, voltage, current, field_heat, reversible_heat) = at(ta);
    Reference {
        temperature_k: [ta, sink(ta)],
        a,
        voltage,
        current,
        field_heat,
        graph_heat: DT * RESISTANCE * current * current,
        reversible_heat,
    }
}

fn advance<'t, 'm, A: HeatAssembly + ?Sized>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    circuit: &SpatialTerminalCircuit<'m>,
    source_v: f64,
    laws: &[SpatialMagneticThermalLaw],
) -> Result<
    SpatialCircuitThermalCandidate<'t, 'm, ResolvedConductorState, A::Sample>,
    CircuitEnthalpyError<SpatialTerminalCircuitThermalError<A::Error>>,
> {
    heat_assembly_candidate(
        provider,
        cx,
        thermal,
        thermal.time_s() + DT,
        circuit,
        &[source_v],
        &[[0.0; 3]; 6],
        &[[0.0; 3]; 8],
        &regions(),
        &[ResistorCell {
            branch: 1,
            cell: CellId(8),
        }],
        1e-12,
        laws,
        thermal_controls(),
        circuit_controls(),
        assembly(),
        &[0.0; 2],
    )
}

fn check_reference<S: AsRef<HeatNetworkSample>>(
    next: &SpatialCircuitThermalCandidate<'_, '_, ResolvedConductorState, S>,
    expected: &Reference,
) {
    let report = next.report();
    let step = &report.coupled.electrical;
    for (cell, t) in [CellId(7), CellId(8)]
        .into_iter()
        .zip(expected.temperature_k)
    {
        near(next.thermal().temperature_k(cell).unwrap(), t, 1e-7);
    }
    near(
        next.circuit().field().potential_wb()[diagonal(next.circuit().field().problem())],
        expected.a,
        1e-7,
    );
    near(step.terminal_voltage_v[0], expected.voltage, 1e-7);
    near(step.terminal_current_a[0], expected.current, 2e-7);
    near(step.field.heat_j, expected.field_heat, 2e-7);
    near(step.resistor_heat_j, expected.graph_heat, 2e-7);
    near(
        report.coupled.source.reversible_heat_j,
        expected.reversible_heat,
        2e-9,
    );
    near(
        next.deposited_heat_j(),
        expected.field_heat + expected.graph_heat,
        3e-7,
    );
    near(report.external_heat_j, 0.0, 2e-8);
    near(report.total_energy_residual_j, 0.0, 2e-8);
    near(report.coupled.total_energy_residual_j, 0.0, 2e-8);
    assert!(step.cancelled_reference_heat_j > 0.0);
    for surface in report
        .coupled
        .source
        .network
        .as_ref()
        .radiation
        .as_ref()
        .unwrap()
        .surfaces()
    {
        let sample = &report.coupled.source.network.as_ref().material;
        assert_eq!(
            surface.phase,
            sample
                .cells
                .iter()
                .find(|cell| cell.cell == surface.cell)
                .unwrap()
                .phase
        );
    }
}

#[test]
fn contacts_radiation_and_integrated_layers_reach_the_actual_electrode_circuit() {
    let mut endpoints = Vec::new();
    for (epsilon, contact, layers, sourced_graph) in [
        (0.2, 0.2, false, false),
        (0.8, 0.2, false, true),
        (0.2, 2.0, false, false),
        (0.2, 0.2, true, true),
    ] {
        endpoints.push(with_heat_material(
            false,
            epsilon,
            contact,
            500.0,
            None,
            |card| card,
            |cx, _, original, mut thermal, card| {
                let problem = problem(card, BASIS);
                let mut circuit = connect(
                    initial(&problem, BASIS, 0.0, 1.0),
                    series_graph(RESISTANCE),
                    vec![],
                    vec![0.0],
                    0.7,
                );
                let laws = laws(card, &problem, 3.0);
                let resistor = [CellResistor {
                    destination: ResistorCell {
                        branch: 1,
                        cell: CellId(8),
                    },
                    selection: MaterialPropertySelection::SingleClaimOnly,
                    length_m: 1.0,
                    area_m2: 1.0,
                }];
                let bases = [(1, BASIS)];
                let bulk = MaterialNetwork {
                    cells: original.electrical.bulk.cells,
                    resistors: if sourced_graph { &resistor } else { &[] },
                    links: original.electrical.bulk.links,
                };
                let electrical = ElectricalContactNetwork {
                    bulk: &bulk,
                    contacts: &[],
                    electrical_bases: if sourced_graph { &bases } else { &[] },
                };
                let heat = MaterialHeatNetwork {
                    electrical: &electrical,
                    thermal_contacts: if layers {
                        &[]
                    } else {
                        original.thermal_contacts
                    },
                    radiation: original.radiation,
                };
                let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
                let sources = [0, 1].map(|i| {
                    ResolvedConductionLayer::new(
                        bulk.cells[i].curve,
                        bulk.cells[i].point,
                        &key,
                        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
                        SolidLiquidPhase::Solid,
                        [200.0, 500.0],
                        || Ok(()),
                    )
                    .unwrap()
                });
                let layer = PreparedLayerConduction::new(
                    &thermal,
                    vec![LayerLinkBinding {
                        link: 17,
                        layers: [&sources[0], &sources[1]],
                        distances_m: [0.5; 2],
                        contact_area_resistance_m2_k_w: contact,
                    }],
                    1,
                    || Ok(()),
                )
                .unwrap();
                let mut prepared = heat.prepare(cx, &thermal).unwrap();
                if layers {
                    prepared = prepared
                        .with_integrated_conduction(
                            cx,
                            &thermal,
                            &layer,
                            LayerConductionConfig {
                                max_iterations: 80,
                                absolute_heat_tolerance_w: 1e-10,
                                relative_heat_tolerance: 1e-12,
                            },
                        )
                        .unwrap();
                }
                let mut old_a = 0.0;
                for source_v in [12.0, -8.0] {
                    let old_t =
                        [CellId(7), CellId(8)].map(|cell| thermal.temperature_k(cell).unwrap());
                    let expected = reference(old_t, old_a, source_v, |ta, tb| {
                        let rt = if layers {
                            contact
                        } else {
                            contact * (1.0 + 0.001 * (ta - 300.0))
                        };
                        (tb - ta) / (1.0 + rt) + radiated(ta, tb, epsilon)
                    });
                    let next = advance(&prepared, cx, &thermal, &circuit, source_v, &laws).unwrap();
                    check_reference(&next, &expected);
                    let source = &next.report().coupled.source.network;
                    assert_eq!(source.thermal_contacts.len(), usize::from(!layers));
                    assert_eq!(source.integrated_conduction.len(), usize::from(layers));
                    assert_eq!(
                        next.report().coupled.electrical.graph_materials.len(),
                        usize::from(sourced_graph)
                    );
                    assert_eq!(
                        [CellId(7), CellId(8)].map(|cell| thermal.temperature_k(cell).unwrap()),
                        old_t
                    );
                    near(
                        circuit.field().potential_wb()[diagonal(&problem)],
                        old_a,
                        1e-7,
                    );
                    let parts = next.into_parts();
                    thermal = parts.0;
                    circuit = parts.1;
                    old_a = expected.a;
                }
                (thermal.temperature_k(CellId(7)).unwrap(), old_a)
            },
        ));
    }
    // Changing an actual source changes both the finite temperature and the
    // retained electromagnetic history, rather than only a diagnostic receipt.
    for endpoint in &endpoints[1..] {
        assert!((endpoint.0 - endpoints[0].0).abs() > 1e-5);
        assert!((endpoint.1 - endpoints[0].1).abs() > 1e-9);
    }
}

#[derive(Debug)]
enum ProviderError {
    Original(HeatNetworkError),
    TemperatureLimit,
}
impl From<HeatNetworkError> for ProviderError {
    fn from(error: HeatNetworkError) -> Self {
        Self::Original(error)
    }
}

struct FiniteSample {
    heat: HeatNetworkSample,
    temperatures_k: [f64; 2],
    conductance_w_k: f64,
    into_field_w: f64,
}
impl AsRef<HeatNetworkSample> for FiniteSample {
    fn as_ref(&self) -> &HeatNetworkSample {
        &self.heat
    }
}
struct FiniteLink<'p, 'n, 'a> {
    prepared: &'p PreparedHeatNetwork<'n, 'a>,
    maximum_field_temperature_k: f64,
}
fn finite_conductance(ta: f64, tb: f64) -> f64 {
    0.7 + 0.002 * (ta + tb - 600.0)
}
impl HeatAssembly for FiniteLink<'_, '_, '_> {
    type Sample = FiniteSample;
    type Error = ProviderError;
    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> {
        self.prepared
    }
    fn resolve_heat(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<FiniteSample>, CircuitEnthalpyError<ProviderError>> {
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let temperatures_k = [
            thermal.temperature_k(CellId(7))?,
            thermal.temperature_k(CellId(8))?,
        ];
        if temperatures_k[0] > self.maximum_field_temperature_k {
            return Err(CircuitEnthalpyError::Source(
                ProviderError::TemperatureLimit,
            ));
        }
        let [ta, tb] = temperatures_k;
        let conductance_w_k = finite_conductance(ta, tb);
        let coefficients = self
            .prepared
            .resolve_with_link_conductances(cx, thermal, &BTreeMap::from([(17, conductance_w_k)]))
            .map_err(|error| match error {
                HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
                other => CircuitEnthalpyError::Source(ProviderError::Original(other)),
            })?;
        let bulk = coefficients.bulk;
        Ok(ContactThermalCoefficients {
            contacts: coefficients.contacts,
            bulk: CircuitThermalCoefficients {
                resistors: bulk.resistors,
                sources_w: bulk.sources_w,
                conductances_w_k: bulk.conductances_w_k,
                source: FiniteSample {
                    heat: bulk.source,
                    temperatures_k,
                    conductance_w_k,
                    into_field_w: conductance_w_k * (tb - ta),
                },
            },
        })
    }
}

#[test]
fn finite_link_provider_retains_payloads_typed_errors_and_one_internal_heat_path() {
    // This manufactured finite-link provider checks the same associated-source
    // seam used by fs-airflow's actual coolant owners, without inventing a second
    // coolant temperature or importing a downstream crate into fs-couple.
    with_heat_material(
        false,
        0.2,
        0.2,
        500.0,
        None,
        |card| card,
        |cx, gate, original, initial_thermal, card| {
            let mut events = initial_thermal
                .cells()
                .values()
                .map(|cell| TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(cell.clone()),
                })
                .collect::<Vec<_>>();
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    conductance_w_k: 0.0,
                    ..*initial_thermal.links().get(&17).unwrap()
                }),
            });
            let laws_table = original
                .electrical
                .bulk
                .cells
                .iter()
                .map(|cell| (cell.material_id, cell.curve as &dyn ThermalLaw))
                .collect::<Vec<_>>();
            let thermal = DepositionSimulation::new(events, &laws_table, initial_thermal.budget())
                .unwrap()
                .advance_to(
                    0.0,
                    StepConfig {
                        max_dt_s: DT,
                        energy_tolerance_j: 1e-9,
                    },
                    || Ok(()),
                )
                .unwrap()
                .0;
            let bulk = MaterialNetwork {
                cells: original.electrical.bulk.cells,
                resistors: &[],
                links: &[],
            };
            let electrical = ElectricalContactNetwork {
                bulk: &bulk,
                contacts: &[],
                electrical_bases: &[],
            };
            let heat = MaterialHeatNetwork {
                electrical: &electrical,
                thermal_contacts: &[],
                radiation: original.radiation,
            };
            let prepared = heat
                .prepare_with_link_conductances(
                    cx,
                    &thermal,
                    &BTreeMap::from([(17, finite_conductance(300.0, 340.0))]),
                )
                .unwrap();
            let provider = FiniteLink {
                prepared: &prepared,
                maximum_field_temperature_k: 350.0,
            };
            let problem = problem(card, BASIS);
            let circuit = connect(
                initial(&problem, BASIS, 0.0, 1.0),
                series_graph(RESISTANCE),
                vec![],
                vec![0.0],
                0.7,
            );
            let laws = laws(card, &problem, 3.0);
            let expected = reference([300.0, 340.0], 0.0, 12.0, |ta, tb| {
                finite_conductance(ta, tb) * (tb - ta) + radiated(ta, tb, 0.2)
            });
            let next = advance(&provider, cx, &thermal, &circuit, 12.0, &laws).unwrap();
            check_reference(&next, &expected);
            let sample = &next.report().coupled.source.network;
            for (actual, expected) in sample
                .temperatures_k
                .into_iter()
                .zip(expected.temperature_k)
            {
                near(actual, expected, 1e-7);
            }
            let [ta, tb] = sample.temperatures_k;
            near(sample.conductance_w_k, finite_conductance(ta, tb), 1e-14);
            near(
                sample.into_field_w,
                sample.conductance_w_k * (tb - ta),
                1e-14,
            );
            assert!((sample.conductance_w_k - finite_conductance(300.0, 340.0)).abs() > 1e-4);
            assert!(sample.heat.thermal_contacts.is_empty());
            let failed = advance(
                &FiniteLink {
                    prepared: &prepared,
                    maximum_field_temperature_k: 299.0,
                },
                cx,
                &thermal,
                &circuit,
                12.0,
                &laws,
            );
            assert!(matches!(
                failed,
                Err(CircuitEnthalpyError::Source(
                    SpatialTerminalCircuitThermalError::Heat(ProviderError::TemperatureLimit)
                ))
            ));
            gate.request();
            assert!(matches!(
                advance(&provider, cx, &thermal, &circuit, 12.0, &laws),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(circuit.time_s(), 0.0);
            near(thermal.temperature_k(CellId(7)).unwrap(), 300.0, 1e-12);
            near(thermal.temperature_k(CellId(8)).unwrap(), 340.0, 1e-12);
            near(
                circuit.field().potential_wb()[diagonal(&problem)],
                0.0,
                1e-12,
            );
        },
    );
}

#[test]
fn field_conditioned_radiation_refuses_without_an_owned_interface_field() {
    with_heat_material(
        false,
        0.2,
        0.2,
        500.0,
        Some(FLUX_DENSITY_DIMS),
        |card| card,
        |cx, _, heat, thermal, card| {
            let problem = problem(card, BASIS);
            let circuit = connect(
                initial(&problem, BASIS, 0.0, 1.0),
                series_graph(RESISTANCE),
                vec![],
                vec![0.0],
                0.7,
            );
            let prepared = heat.prepare(cx, &thermal).unwrap();
            let result = advance(
                &prepared,
                cx,
                &thermal,
                &circuit,
                12.0,
                &laws(card, &problem, 3.0),
            );
            assert!(matches!(
                result,
                Err(CircuitEnthalpyError::Source(
                    SpatialTerminalCircuitThermalError::Heat(HeatNetworkError::Invalid(_))
                ))
            ));
            assert_eq!(circuit.time_s(), 0.0);
            assert_eq!(thermal.time_s(), 0.0);
        },
    );
}
