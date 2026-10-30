//! A distinct finite native link reaches nonlinear current and field storage.
//! Manufactured heat-law data; actual fluid implementations use this SAME trait.
use super::*;
use crate::REGIONS;
use fs_couple::circuit_enthalpy::{
    CircuitThermalCoefficients,
    conduction::{
        contacts::ConductionContactThermal,
        heat::{ConductionHeatError, step_response_field_circuit_candidate as advance},
    },
    contacts::{
        ContactThermalCoefficients,
        heat::{
            HeatAssembly, HeatNetworkError, HeatNetworkSample, MaterialHeatNetwork,
            PreparedHeatNetwork,
        },
        material::ElectricalContactNetwork,
    },
};
use fs_exec::{CancelGate, Cx};
use std::collections::BTreeMap;

const LINK: u64 = 9;
const DT: f64 = 0.1;
const FIELD_SLOPE: f64 = 0.001;
const LINK_SLOPE: f64 = 0.002;

#[derive(Debug)]
enum AddedHeatError {
    Heat(HeatNetworkError),
    Refused,
}
impl From<HeatNetworkError> for AddedHeatError {
    fn from(error: HeatNetworkError) -> Self {
        Self::Heat(error)
    }
}
impl core::fmt::Display for AddedHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Heat(error) => write!(f, "original heat source: {error}"),
            Self::Refused => f.write_str("authored additional heat source refused"),
        }
    }
}
struct AddedHeatSample {
    heat: HeatNetworkSample,
    conductance_w_k: f64,
    power_to_field_w: f64,
}
impl AsRef<HeatNetworkSample> for AddedHeatSample {
    fn as_ref(&self) -> &HeatNetworkSample {
        &self.heat
    }
}
struct AddedHeat<'p, 'n, 'a> {
    prepared: &'p PreparedHeatNetwork<'n, 'a>,
    scale: f64,
    refuse: bool,
}
impl HeatAssembly for AddedHeat<'_, '_, '_> {
    type Sample = AddedHeatSample;
    type Error = AddedHeatError;
    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> {
        self.prepared
    }
    fn resolve_heat(
        &self,
        cx: &Cx<'_>,
        state: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        if self.refuse {
            return Err(CircuitEnthalpyError::Source(AddedHeatError::Refused));
        }
        let tb = state.temperature_k(CellId(10))?;
        let tf = state.temperature_k(FIELD_CELL)?;
        let conductance_w_k = self.scale * (1.0 + LINK_SLOPE * (0.5 * (tb + tf) - 300.0));
        let resolved = self
            .prepared
            .resolve_with_link_conductances(cx, state, &BTreeMap::from([(LINK, conductance_w_k)]))
            .map_err(|error| match error {
                HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
                other => CircuitEnthalpyError::Source(AddedHeatError::Heat(other)),
            })?;
        let bulk = resolved.bulk;
        Ok(ContactThermalCoefficients {
            contacts: resolved.contacts,
            bulk: CircuitThermalCoefficients {
                resistors: bulk.resistors,
                sources_w: bulk.sources_w,
                conductances_w_k: bulk.conductances_w_k,
                source: AddedHeatSample {
                    heat: bulk.source,
                    conductance_w_k,
                    power_to_field_w: conductance_w_k * (tb - tf),
                },
            },
        })
    }
}

fn fixture<R>(
    magnetic: bool,
    scale: f64,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &mut AddedHeat<'_, '_, '_>,
        DepositionSimulation<'_>,
        ConductionResponseCircuit,
        &[ConductionContactThermal<'_>],
        &[FieldThermalBinding],
    ) -> R,
) -> R {
    let coordinate = ConductorMagnitude::ElectricField;
    let (a, ca) = material(1.0, CP[0], ALPHA[0], BETA[0], coordinate, false);
    let (b, cb) = material(2.0, CP[1], ALPHA[1], BETA[1], coordinate, false);
    let (cards, curves) = ([a, b], [ca, cb]);
    let (field_card, field_curve) = field_card(magnetic, FIELD_SLOPE);
    let interface = interface(&cards, 0.5, CONTACT_SLOPE, 1000.0);
    let context = point(300.0, coordinate);
    let contact_context = interface_context(300.0);
    let field_context = field_point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let mut cells = cells(&cards, &curves, &context, &key);
    cells.push(MaterialCell {
        cell: FIELD_CELL,
        material_id: MaterialId(2),
        curve: field_curve.curve(),
        point: ThermalMaterialPoint {
            material: &field_card,
            context: &field_context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    });
    let resistors = [resistor()];
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &crate::LINKS,
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let contacts = [binding(&cards, &interface, &contact_context, FRACTION)];
    // Circuit initialization queries the same three caloric cells. Additional
    // heat topology is then supplied by the shared assembly before any advance.
    let initial = history(&curves, &field_curve, 1.0);
    let (circuit, field) = initial_circuit(
        specimen(&cards, &interface, coordinate),
        &bulk,
        &initial,
        &field_card,
        magnetic,
    );
    let thermal = history_with_links(
        &curves,
        &field_curve,
        1.0,
        &[ThermalLink {
            id: LINK,
            a: CellId(10),
            b: FIELD_CELL,
            area_m2: 1.0,
            conductance_w_k: 0.0,
        }],
    );
    assert_eq!(initial.cells(), thermal.cells());
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 191,
                kernel_id: 57,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let prepared = heat
            .prepare_with_link_conductances(&cx, &thermal, &BTreeMap::from([(LINK, scale)]))
            .unwrap();
        let mut complete = AddedHeat {
            prepared: &prepared,
            scale,
            refuse: false,
        };
        run(
            &cx,
            &gate,
            &mut complete,
            thermal,
            circuit,
            &contacts,
            &[field],
        )
    })
}

#[test]
fn added_finite_exchange_changes_continuing_spatial_current_and_source_field_storage() {
    for magnetic in [false, true] {
        let mut outcomes = Vec::new();
        for scale in [0.1, 5.0] {
            let outcome = fixture(
                magnetic,
                scale,
                |cx, gate, heat, mut thermal, mut circuit, contacts, fields| {
                    let mut old = [300.0; 3];
                    let mut z0 = 3.0;
                    let initial = thermal.energy_j() + internal(FIELD_SLOPE, old[2], z0);
                    let original_links = thermal.links().clone();
                    let mut current = 0.0;
                    for interval in 1..=4 {
                        let target = interval as f64 * DT;
                        let dt = target - thermal.time_s();
                        let saved = thermal.cells().clone();
                        let next = advance(
                            heat,
                            cx,
                            &thermal,
                            target,
                            &circuit,
                            &[],
                            contacts,
                            fields,
                            &REGIONS,
                            1e-12,
                            limits(),
                            field_controls(),
                            controls(),
                            config(),
                            &[0.0; 3],
                        )
                        .unwrap();
                        let t = [CellId(30), CellId(10), FIELD_CELL]
                            .map(|id| next.thermal().temperature_k(id).unwrap());
                        let theta = t.map(|value| value - 300.0);
                        let z = if magnetic {
                            next.circuit().fluxes_wb()[0]
                        } else {
                            next.circuit().charges_c()[0]
                        };
                        let report = next.report();
                        let electrical = &report.coupled.electrical.circuit;
                        let source = &report.coupled.source.network;
                        let j = electrical.spatial_current_a;
                        current = j;
                        // Independent local material and circuit equations.
                        let e = [0, 1].map(|i| {
                            j * ((i + 1) as f64 + ALPHA[i] * theta[i]) / (1.0 - BETA[i] * j)
                        });
                        let rc = 0.5 * (1.0 + CONTACT_SLOPE * 0.5 * (theta[0] + theta[1]));
                        let v = e[0] + e[1] + j * rc;
                        let load = LOAD_R + LOAD_SLOPE * theta[2];
                        let total_current = j + v / load;
                        let k = 1.0 + FIELD_SLOPE * theta[2];
                        near(electrical.spatial_voltage_v, v);
                        near(z, z0 - dt * if magnetic { v } else { total_current });
                        near(0.5 * (z0 + z), k * if magnetic { total_current } else { v });
                        near(
                            next.circuit().material_resistors()[0]
                                .conductor
                                .resistance_ohm(),
                            load,
                        );
                        near(report.coupled.source.fields[0].material.coefficient(), k);
                        let g = scale * (1.0 + LINK_SLOPE * 0.5 * (theta[1] + theta[2]));
                        let power = g * (t[1] - t[2]);
                        near(source.conductance_w_k, g);
                        near(source.power_to_field_w, power);
                        let transfer = t[0] - t[1]; // native source k=1, unit link length/area
                        near(
                            CP[0] * (t[0] - old[0]),
                            dt * (j * e[0] + FRACTION * j * j * rc - transfer),
                        );
                        near(
                            CP[1] * (t[1] - old[1]),
                            dt * (j * e[1] + (1.0 - FRACTION) * j * j * rc + transfer - power),
                        );
                        let du = internal(FIELD_SLOPE, t[2], z) - internal(FIELD_SLOPE, old[2], z0);
                        near(
                            FIELD_CP * (t[2] - old[2]) + du,
                            dt * (-v * total_current + v * v / load + power + 1.0),
                        );
                        near(
                            FIELD_CP * (t[2] - old[2]),
                            dt * (v * v / load + power + 1.0)
                                + report.coupled.source.reversible_heat_j,
                        );
                        near(
                            report.coupled.deposited_joule_j,
                            dt * (j * v + v * v / load),
                        );
                        near(report.external_heat_j, dt);
                        near(report.total_energy_residual_j, 0.0);
                        near(
                            next.thermal().energy_j() + internal(FIELD_SLOPE, t[2], z),
                            initial + target,
                        );
                        assert!(report.coupled.electrical.field_solves > 1);
                        assert_eq!(source.heat.material.links.len(), 1);
                        assert_eq!(next.thermal().links(), &original_links);
                        assert_eq!(thermal.cells(), &saved);
                        assert_eq!(circuit.time_s(), thermal.time_s());
                        old = t;
                        z0 = z;
                        (thermal, circuit, _) = next.into_parts();
                    }
                    let saved = thermal.cells().clone();
                    let charges = circuit.charges_c().to_vec();
                    let fluxes = circuit.fluxes_wb().to_vec();
                    heat.refuse = true;
                    assert!(matches!(
                        advance(
                            heat,
                            cx,
                            &thermal,
                            5.0 * DT,
                            &circuit,
                            &[],
                            contacts,
                            fields,
                            &REGIONS,
                            1e-12,
                            limits(),
                            field_controls(),
                            controls(),
                            config(),
                            &[0.0; 3]
                        ),
                        Err(CircuitEnthalpyError::Source(ConductionHeatError::Heat(
                            AddedHeatError::Refused
                        )))
                    ));
                    heat.refuse = false;
                    gate.request();
                    assert!(matches!(
                        advance(
                            heat,
                            cx,
                            &thermal,
                            5.0 * DT,
                            &circuit,
                            &[],
                            contacts,
                            fields,
                            &REGIONS,
                            1e-12,
                            limits(),
                            field_controls(),
                            controls(),
                            config(),
                            &[0.0; 3]
                        ),
                        Err(CircuitEnthalpyError::Cancelled)
                    ));
                    assert_eq!(thermal.cells(), &saved);
                    assert_eq!(circuit.charges_c(), charges);
                    assert_eq!(circuit.fluxes_wb(), fluxes);
                    [old[2], z0, current]
                },
            );
            outcomes.push(outcome);
        }
        for (changed, original) in outcomes[1].iter().zip(outcomes[0]) {
            assert!((changed - original).abs() > 1e-8);
        }
    }
}
