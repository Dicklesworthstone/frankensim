//! Nonlinear spatial/contact response, sourced R/C/L and physical field U.
//! Coefficients are manufactured numerical fixtures, not measured materials.
#[path = "field_storage/assembly.rs"]
mod assembly;
#[path = "field_storage/heat.rs"]
mod heat;
use super::*;
use fs_couple::circuit_enthalpy::{
    ResistorCell,
    field_storage::FieldThermalBinding,
    material::{CellResistor, MaterialCell},
};
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_material::{
    conductor::ELECTRICAL_RESISTIVITY_DIMS,
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS,
        MAGNETIC_PERMEABILITY_PROPERTY, ResolvedCapacitor, ResolvedInductor,
        thermal::ThermalFieldStorage,
    },
    thermal::source::ThermalMaterialPoint,
};
use fs_thermal::enthalpy::{
    CellGeometry, DepositedCell, DepositionBudget, DepositionSimulation, Event, MaterialId,
    StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};

const FIELD_CELL: CellId = CellId(70);
const FIELD_CP: f64 = 2.0;
const LOAD_R: f64 = 4.0;
const LOAD_SLOPE: f64 = 0.02;

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn field_point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_quantity(), t)
        .unwrap()
}
fn field_card(magnetic: bool, alpha: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (name, dims) = if magnetic {
        (MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS)
    } else {
        (ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS)
    };
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let mut claims = ClaimSet::new();
    for (key, low, high) in [
        (
            PropertyKey::new(name, dims),
            1.0 - 50.0 * alpha,
            1.0 + 200.0 * alpha,
        ),
        (
            PropertyKey::new("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS),
            LOAD_R - 50.0 * LOAD_SLOPE,
            LOAD_R + 200.0 * LOAD_SLOPE,
        ),
        (heat.clone(), FIELD_CP, FIELD_CP),
        (density.clone(), 1.0, 1.0),
    ] {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(250.0, low), (500.0, high)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature_quantity(),
                    250.0,
                    500.0,
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured lossy field material".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic field material".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &field_point(250.0),
            upper: &field_point(500.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-9,
    )
    .unwrap();
    (card, curve)
}

fn history<'a>(
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    field: &'a ResolvedHeatCapacityCurve,
    source_w: f64,
) -> DepositionSimulation<'a> {
    history_with_links(curves, field, source_w, &[])
}

fn history_with_links<'a>(
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    field: &'a ResolvedHeatCapacityCurve,
    source_w: f64,
    additional_links: &[ThermalLink],
) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for (i, id) in [CellId(30), CellId(10), FIELD_CELL].into_iter().enumerate() {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: 1.0,
                specific_h_j_kg: if i < 2 {
                    100.0 * CP[i]
                } else {
                    50.0 * FIELD_CP
                },
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: fs_thermal::enthalpy::Boundary {
                    source_w: if i == 2 { source_w } else { 0.0 },
                    ..crate::common::boundary()
                },
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: CellId(30),
            b: CellId(10),
            area_m2: 1.0,
            conductance_w_k: 999.0,
        }),
    });
    events.extend(additional_links.iter().map(|link| TimedEvent {
        time_s: 0.0,
        event: Event::Contact(link.clone()),
    }));
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
            (MaterialId(2), field.curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 3,
            max_links: 1 + additional_links.len(),
            max_events: 4 + additional_links.len(),
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn resistor() -> CellResistor {
    CellResistor {
        destination: ResistorCell {
            branch: 1,
            cell: FIELD_CELL,
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    }
}
fn initial_circuit(
    problem: PlanarConductionResponse,
    material: &MaterialNetwork<'_>,
    thermal: &DepositionSimulation<'_>,
    card: &MaterialCard,
    magnetic: bool,
) -> (ConductionResponseCircuit, FieldThermalBinding) {
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let (branch, element, law) = if magnetic {
        let value =
            ResolvedInductor::resolve(card, &field_point(300.0), selection, 1.0, 1.0, 1).unwrap();
        let law = ThermalFieldStorage::inductor(card, &value, "T", [250.0, 500.0]).unwrap();
        (
            (0, 1, Branch::Inductor { henries: f64::NAN }),
            MaterialStorageElement::Inductor(value),
            law,
        )
    } else {
        let value =
            ResolvedCapacitor::resolve(card, &field_point(300.0), selection, 1.0, 1.0).unwrap();
        let law = ThermalFieldStorage::capacitor(card, &value, "T", [250.0, 500.0]).unwrap();
        (
            (1, 0, Branch::Capacitor { farads: f64::NAN }),
            MaterialStorageElement::Capacitor(value),
            law,
        )
    };
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::Resistor { ohms: 999.0 }),
            (1, 0, Branch::Resistor { ohms: 999.0 }),
            branch,
        ],
        transformers: vec![],
    };
    let circuit = ConductionResponseCircuit::new_with_material_storage(
        problem,
        graph,
        ConductionCircuitTerminal {
            branch: 0,
            positive: vec![5, 6],
            negative: vec![0, 3],
        },
        material.resolve(thermal, || false).unwrap().resistors,
        vec![MaterialStorage { branch: 2, element }],
        0.0,
        if magnetic { vec![3.0] } else { vec![] },
        if magnetic { vec![] } else { vec![3.0] },
        &[],
        limits(),
        field_controls(),
        controls(),
        || false,
    )
    .unwrap();
    (
        circuit,
        FieldThermalBinding {
            branch: 2,
            cell: FIELD_CELL,
            law,
            maximum_abs_coordinate: 3.0,
        },
    )
}
fn bisect(f: impl Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> f64 {
    assert!(f(lo) <= 0.0 && f(hi) >= 0.0);
    for _ in 0..90 {
        let x = 0.5 * (lo + hi);
        if f(x) > 0.0 {
            hi = x;
        } else {
            lo = x;
        }
    }
    0.5 * (lo + hi)
}
fn internal(alpha: f64, t: f64, z: f64) -> f64 {
    let k = 1.0 + alpha * (t - 300.0);
    z * z / (2.0 * k) * (1.0 + t * alpha / k)
}
// Independent reduction: spatial two-body nonlinear/contact/thermal balance
// uses the parent's analytic 2x2 elimination. Add a parallel sourced resistor
// and finite C(T) or L(T). Close the midpoint coordinate and physical field U
// plus all three caloric inventories, not the implementation's psi audit.
fn oracle(
    magnetic: bool,
    alpha: f64,
    dt: f64,
    old: [f64; 2],
    tf0: f64,
    z0: f64,
    external: f64,
) -> ([f64; 2], f64, f64, f64, f64) {
    let at = |tf: f64| {
        let k = 1.0 + alpha * (tf - 300.0);
        let load = LOAD_R + LOAD_SLOPE * (tf - 300.0);
        let quantities = |j: f64| {
            let (theta, e, rc) = endpoint(j, dt, old, [1.0, 2.0], 0.5);
            let v = e[0] + e[1] + j * rc;
            (theta, v, j + v / load)
        };
        let j = bisect(
            |j| {
                let (_, v, i) = quantities(j);
                if magnetic {
                    k * i + 0.5 * dt * v - z0
                } else {
                    k * v + 0.5 * dt * i - z0
                }
            },
            0.0,
            3.0,
        );
        let (theta, v, i) = quantities(j);
        let z = if magnetic { z0 - dt * v } else { z0 - dt * i };
        (theta, z, j, v)
    };
    let tf = bisect(
        |tf| {
            let (theta, z, _, _) = at(tf);
            CP[0] * (theta[0] - old[0])
                + CP[1] * (theta[1] - old[1])
                + FIELD_CP * (tf - tf0)
                + internal(alpha, tf, z)
                - internal(alpha, tf0, z0)
                - external
        },
        280.0,
        340.0,
    );
    let (theta, z, j, v) = at(tf);
    (theta, tf, z, j, v)
}

#[test]
fn nonlinear_contacts_material_resistance_and_thermal_field_storage_share_one_endpoint() {
    for magnetic in [false, true] {
        let mut outcomes = Vec::new();
        for alpha in [-0.001, 0.0, 0.001] {
            let coordinate = ConductorMagnitude::ElectricField;
            let (a, ca) = material(1.0, CP[0], ALPHA[0], BETA[0], coordinate, false);
            let (b, cb) = material(2.0, CP[1], ALPHA[1], BETA[1], coordinate, false);
            let (cards, curves) = ([a, b], [ca, cb]);
            let (field_card, field_curve) = field_card(magnetic, alpha);
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
            let material = MaterialNetwork {
                cells: &cells,
                resistors: &resistors,
                links: &crate::LINKS,
            };
            let contacts = [binding(&cards, &interface, &contact_context, FRACTION)];
            let mut thermal = history(&curves, &field_curve, 1.0);
            let (mut circuit, field) = initial_circuit(
                specimen(&cards, &interface, coordinate),
                &material,
                &thermal,
                &field_card,
                magnetic,
            );
            let fields = [field];
            let layout = network(&cells);
            let (dt, mut old, mut tf0, mut z0) = (0.1, [0.0; 2], 300.0, 3.0);
            let initial_energy = thermal.energy_j() + internal(alpha, tf0, z0);
            for step in 1..=3 {
                let (theta, tf, z, j, v) = oracle(magnetic, alpha, dt, old, tf0, z0, dt);
                let next = layout
                    .step_response_field_circuit_candidate(
                        &thermal,
                        step as f64 * dt,
                        &circuit,
                        &[],
                        &resistors,
                        &contacts,
                        &fields,
                        limits(),
                        field_controls(),
                        controls(),
                        config(),
                        &[0.0; 3],
                        || false,
                    )
                    .unwrap();
                let report = next.report();
                let electrical = &report.coupled.electrical;
                let field = &report.coupled.source.fields[0];
                for (i, id) in [CellId(30), CellId(10)].into_iter().enumerate() {
                    near(next.thermal().temperature_k(id).unwrap(), 300.0 + theta[i]);
                }
                near(next.thermal().temperature_k(FIELD_CELL).unwrap(), tf);
                near(
                    if magnetic {
                        next.circuit().fluxes_wb()[0]
                    } else {
                        next.circuit().charges_c()[0]
                    },
                    z,
                );
                near(electrical.circuit.spatial_current_a, j);
                near(electrical.circuit.spatial_voltage_v, v);
                let rc = 0.5 * (1.0 + CONTACT_SLOPE * (theta[0] + theta[1]) / 2.0);
                near(
                    electrical.circuit.spatial.contacts[0].response.power_w,
                    j * j * rc,
                );
                let load = LOAD_R + LOAD_SLOPE * (tf - 300.0);
                let load_heat = dt * v * v / load;
                near(electrical.circuit.graph_resistors[0].heat_j, load_heat);
                near(
                    next.circuit().material_resistors()[0]
                        .conductor
                        .resistance_ohm(),
                    load,
                );
                near(electrical.circuit.heat_j, dt * j * v + load_heat);
                near(
                    field.reversible_heat_j,
                    FIELD_CP * (tf - tf0) - load_heat - dt,
                );
                near(field.state[0].internal_energy_j, internal(alpha, tf0, z0));
                near(field.state[1].internal_energy_j, internal(alpha, tf, z));
                near(field.material.coefficient(), 1.0 + alpha * (tf - 300.0));
                near(report.external_heat_j, dt);
                near(report.total_energy_residual_j, 0.0);
                near(
                    next.thermal().energy_j() + internal(alpha, tf, z),
                    initial_energy + step as f64 * dt,
                );
                assert!(electrical.field_solves > 1);
                assert_eq!(
                    field.material.material().card_identity(),
                    field_card.content_hash()
                );
                assert_eq!(next.circuit().storage().len(), 1);
                if alpha == 0.0 {
                    assert_eq!(field.reversible_heat_j, 0.0);
                }
                assert_eq!(circuit.time_s(), (step - 1) as f64 * dt);
                old = theta;
                tf0 = tf;
                z0 = z;
                (thermal, circuit, _) = next.into_parts();
            }
            outcomes.push((tf0, z0));
        }
        assert!((outcomes[2].0 - outcomes[0].0).abs() > 1e-4);
        assert!((outcomes[2].1 - outcomes[0].1).abs() > 1e-7);
    }
}

#[test]
fn nonlinear_field_trials_refuse_missing_owners_and_preserve_all_accepted_histories() {
    let coordinate = ConductorMagnitude::ElectricField;
    let (a, ca) = material(1.0, CP[0], ALPHA[0], BETA[0], coordinate, false);
    let (b, cb) = material(2.0, CP[1], ALPHA[1], BETA[1], coordinate, false);
    let (cards, curves) = ([a, b], [ca, cb]);
    let (field_card, field_curve) = field_card(false, 0.001);
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
    let material = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &crate::LINKS,
    };
    let thermal = history(&curves, &field_curve, 0.0);
    let contacts = [binding(&cards, &interface, &contact_context, FRACTION)];
    let (circuit, field) = initial_circuit(
        specimen(&cards, &interface, coordinate),
        &material,
        &thermal,
        &field_card,
        false,
    );
    let fields = [field];
    let layout = network(&cells);
    let original = (
        thermal.energy_j(),
        circuit.graph().clone(),
        circuit.conductors().to_vec(),
        circuit.contacts().to_vec(),
        circuit.storage().to_vec(),
    );
    let mut budget = config();
    budget.maximum_circuit_trials = 1;
    assert!(matches!(
        layout.step_response_field_circuit_candidate(
            &thermal,
            0.1,
            &circuit,
            &[],
            &resistors,
            &contacts,
            &fields,
            limits(),
            field_controls(),
            controls(),
            budget,
            &[0.0; 3],
            || false
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let mut polls = 0;
    assert!(matches!(
        layout.step_response_field_circuit_candidate(
            &thermal,
            0.1,
            &circuit,
            &[],
            &resistors,
            &contacts,
            &fields,
            limits(),
            field_controls(),
            controls(),
            config(),
            &[0.0; 3],
            || {
                polls += 1;
                polls > 100
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert!(
        layout
            .step_response_field_circuit_candidate(
                &thermal,
                0.1,
                &circuit,
                &[],
                &resistors,
                &contacts,
                &[],
                limits(),
                field_controls(),
                controls(),
                config(),
                &[0.0; 3],
                || false
            )
            .is_err()
    );
    assert!(
        layout
            .step_response_field_circuit_candidate(
                &thermal,
                0.1,
                &circuit,
                &[],
                &resistors,
                &[],
                &fields,
                limits(),
                field_controls(),
                controls(),
                config(),
                &[0.0; 3],
                || false
            )
            .is_err()
    );
    assert!(
        layout
            .step_response_contact_circuit_candidate(
                &thermal,
                0.1,
                &circuit,
                &[],
                &[resistors[0].destination],
                &contacts,
                limits(),
                field_controls(),
                controls(),
                config(),
                &[0.0; 3],
                || false
            )
            .is_err()
    );
    assert_eq!(thermal.energy_j(), original.0);
    assert_eq!(circuit.graph(), &original.1);
    assert_eq!(circuit.conductors(), original.2);
    assert_eq!(circuit.contacts(), original.3);
    assert_eq!(circuit.storage(), original.4);
    assert_eq!(circuit.charges_c(), &[3.0]);
    assert_eq!(circuit.time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    layout
        .step_response_field_circuit_candidate(
            &thermal,
            0.1,
            &circuit,
            &[],
            &resistors,
            &contacts,
            &fields,
            limits(),
            field_controls(),
            controls(),
            config(),
            &[0.0; 3],
            || false,
        )
        .unwrap();
}
