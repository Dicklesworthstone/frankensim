//! G1/G3/G4: native spatial current, source R(T), field entropy and finite heat.
//! All coefficients are manufactured; no empirical material validation is claimed.
use super::*;
use fs_couple::circuit_enthalpy::{
    conduction::ConductionThermalError,
    field_storage::{FieldThermalBinding, FieldThermalError},
    material::{CellResistor, MaterialNetwork},
};
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue, Provenance,
    SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    conductor::ELECTRICAL_RESISTIVITY_DIMS,
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS,
        MAGNETIC_PERMEABILITY_PROPERTY, ResolvedCapacitor, ResolvedInductor,
        thermal::ThermalFieldStorage,
    },
    phase::{HeatCapacitySourceQuery, SolidLiquidPhase},
    state_point::SPECIFIC_HEAT_CAPACITY_DIMS,
};
use fs_qty::{
    Density, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::{
    CellGeometry, DepositedCell, DepositionBudget, DepositionSimulation, Event, StepConfig,
    ThermalLaw, TimedEvent,
};

const BETA: f64 = 0.002;
const IDS: [CellId; 3] = [CellId(30), CellId(50), CellId(70)];
const REGIONS: [ConductionRegion; 4] = [
    ConductionRegion {
        triangle: 0,
        cell: IDS[0],
    },
    ConductionRegion {
        triangle: 1,
        cell: IDS[0],
    },
    ConductionRegion {
        triangle: 2,
        cell: IDS[0],
    },
    ConductionRegion {
        triangle: 3,
        cell: IDS[0],
    },
];

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            t,
        )
        .unwrap()
}
fn card(magnetic: bool, alpha: f64, cp: f64) -> MaterialCard {
    let (name, dims) = if magnetic {
        (MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS)
    } else {
        (ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS)
    };
    let mut claims = ClaimSet::new();
    for (name, dims, low, high) in [
        (name, dims, 1.0 - 50.0 * alpha, 1.0 + 200.0 * alpha),
        (
            "electrical_resistivity",
            ELECTRICAL_RESISTIVITY_DIMS,
            1.0 - 50.0 * BETA,
            1.0 + 200.0 * BETA,
        ),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            cp,
            cp,
        ),
        ("density", Density::DIMS, 1.0, 1.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(250.0, low), (500.0, high)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    250.0,
                    500.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured spatial field/caloric source".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "numerical electrothermal fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(250.0),
            upper: &point(500.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-9,
    )
    .unwrap()
}
fn history(
    curves: &[ResolvedHeatCapacityCurve; 3],
    field_source_w: f64,
) -> DepositionSimulation<'_> {
    let events = IDS
        .iter()
        .enumerate()
        .map(|(i, &id)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: if i == 0 { 2.0 } else { 1.0 },
                specific_h_j_kg: if i == 2 { 100.0 } else { 50.0 },
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [if i == 0 { 2.0 } else { 1.0 }, 0.0, 0.0],
                    volume_m3: if i == 0 { 2.0 } else { 1.0 },
                    surface_area_m2: 6.0,
                },
                boundary: fs_thermal::enthalpy::Boundary {
                    source_w: if i == 2 { field_source_w } else { 0.0 },
                    ..common::boundary()
                },
            }),
        })
        .collect();
    let laws = curves
        .iter()
        .enumerate()
        .map(|(i, c)| (MaterialId(i as u64), c.curve() as &dyn ThermalLaw))
        .collect::<Vec<_>>();
    DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: 3,
            max_links: 0,
            max_events: 3,
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
fn cells<'a>(
    cards: &'a [MaterialCard; 3],
    curves: &'a [ResolvedHeatCapacityCurve; 3],
    context: &'a QueryPoint,
) -> Vec<MaterialCell<'a>> {
    IDS.iter()
        .enumerate()
        .map(|(i, &cell)| MaterialCell {
            cell,
            material_id: MaterialId(i as u64),
            curve: curves[i].curve(),
            point: ThermalMaterialPoint {
                material: &cards[i],
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        })
        .collect()
}
fn network<'a>(cells: &'a [MaterialCell<'a>]) -> ConductionThermalNetwork<'a> {
    ConductionThermalNetwork {
        cells,
        links: &[],
        regions: &REGIONS,
        volume_tolerance_m3: 1e-12,
    }
}
fn resistor() -> CellResistor {
    CellResistor {
        destination: ResistorCell {
            branch: 1,
            cell: IDS[1],
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    }
}
fn initial_circuit(
    material: &MaterialNetwork<'_>,
    thermal: &DepositionSimulation<'_>,
    cards: &[MaterialCard; 3],
    magnetic: bool,
) -> (ConductionCircuit, FieldThermalBinding) {
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let (kind, element, law) = if magnetic {
        let value =
            ResolvedInductor::resolve(&cards[2], &point(300.0), selection.clone(), 1.0, 1.0, 1)
                .unwrap();
        let law = ThermalFieldStorage::inductor(&cards[2], &value, "T", [250.0, 500.0]).unwrap();
        (
            Branch::Inductor { henries: f64::NAN },
            MaterialStorageElement::Inductor(value),
            law,
        )
    } else {
        let value =
            ResolvedCapacitor::resolve(&cards[2], &point(300.0), selection.clone(), 1.0, 1.0)
                .unwrap();
        let law = ThermalFieldStorage::capacitor(&cards[2], &value, "T", [250.0, 500.0]).unwrap();
        (
            Branch::Capacitor { farads: f64::NAN },
            MaterialStorageElement::Capacitor(value),
            law,
        )
    };
    let source = resolve_conductor_state_with_basis(
        &cards[0],
        &point(300.0),
        selection,
        ElectricalPropertyBasis::Resistivity,
    )
    .unwrap();
    let problem =
        PlanarConduction::new(mesh(1.0), vec![source; 4], limits(), &mut || false).unwrap();
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (2, 0, Branch::Resistor { ohms: 99.0 }),
            (1, 2, Branch::Resistor { ohms: 99.0 }),
            if magnetic { (0, 1, kind) } else { (1, 0, kind) },
        ],
        transformers: vec![],
    };
    let circuit = ConductionCircuit::new_with_material_storage(
        problem,
        graph,
        terminal(),
        material.resolve(thermal, || false).unwrap().resistors,
        vec![MaterialStorage { branch: 2, element }],
        0.0,
        if magnetic { vec![4.0] } else { vec![] },
        if magnetic { vec![] } else { vec![4.0] },
        &[],
        field_controls(),
        || false,
    )
    .unwrap();
    (
        circuit,
        FieldThermalBinding {
            branch: 2,
            cell: IDS[2],
            law,
            maximum_abs_coordinate: 4.0,
        },
    )
}

fn bisect(f: impl Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> f64 {
    assert!(f(lo) <= 0.0 && f(hi) >= 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}
fn potential(alpha: f64, temperature: f64, z: f64) -> f64 {
    let k = 1.0 + alpha * (temperature - 300.0);
    z * z / (2.0 * k) * (1.0 + temperature * alpha / k)
}
// Independent elimination: Rspatial=2*rho, Rload=rho; their capacities 2 and 1
// give equal temperature rise. Eliminate both BE heat equations, then solve
// midpoint charge/flux and the enclosing total-INTERNAL-energy equation.
fn oracle(
    magnetic: bool,
    alpha: f64,
    dt: f64,
    tr0: f64,
    tf0: f64,
    z0: f64,
    external_heat: f64,
) -> (f64, f64, f64, f64) {
    oracle_with_motion(magnetic, alpha, dt, tr0, tf0, z0, external_heat, 0.0)
}

fn oracle_with_motion(
    magnetic: bool,
    alpha: f64,
    dt: f64,
    tr0: f64,
    tf0: f64,
    z0: f64,
    external_heat: f64,
    motional_voltage: f64,
) -> (f64, f64, f64, f64) {
    let at_field_temperature = |tf: f64| {
        let k = 1.0 + alpha * (tf - 300.0);
        let resistance = |i: f64| 3.0 * (1.0 + BETA * (tr0 - 300.0)) / (1.0 - BETA * dt * i * i);
        let current = bisect(
            |i| {
                if magnetic {
                    i * (k + 0.5 * dt * resistance(i)) - z0 + 0.5 * dt * motional_voltage
                } else {
                    i * (k * resistance(i) + 0.5 * dt) - z0 + k * motional_voltage
                }
            },
            0.0,
            10.0,
        );
        let heat = dt * current * current * resistance(current);
        let z = if magnetic {
            z0 - dt * (resistance(current) * current + motional_voltage)
        } else {
            z0 - dt * current
        };
        (tr0 + heat / 3.0, z, current)
    };
    let tf = bisect(
        |tf| {
            let (tr, z, current) = at_field_temperature(tf);
            3.0 * (tr - tr0) + 2.0 * (tf - tf0) + potential(alpha, tf, z)
                - potential(alpha, tf0, z0)
                - external_heat
                + dt * motional_voltage * current
        },
        280.0,
        330.0,
    );
    let (tr, z, current) = at_field_temperature(tf);
    (tr, tf, z, current)
}

#[test]
fn dormant_source_field_cannot_hide_active_physical_energy_or_joule_heat() {
    for alpha in [-0.001, 0.001] {
        for sign in [-1.0, 1.0] {
            let cards = [
                card(false, 0.0, 1.0),
                card(false, 0.0, 1.0),
                card(false, alpha, 2.0),
            ];
            let curves = cards.each_ref().map(caloric);
            let mut thermal = history(&curves, 0.0);
            let context = point(300.0);
            let cells = cells(&cards, &curves, &context);
            let resistors = [resistor()];
            let material = MaterialNetwork {
                cells: &cells,
                resistors: &resistors,
                links: &[],
            };
            let (original, active) = initial_circuit(&material, &thermal, &cards, false);
            let dormant = ResolvedCapacitor::resolve(
                &cards[1],
                &point(300.0),
                MaterialPropertySelection::SingleClaimOnly,
                1.0,
                1.0,
            )
            .unwrap();
            let dormant_law =
                ThermalFieldStorage::capacitor(&cards[1], &dormant, "T", [250.0, 500.0]).unwrap();
            let mut graph = original.graph().clone();
            graph.node_count = 4;
            graph
                .branches
                .push((3, 0, Branch::Capacitor { farads: f64::NAN }));
            let mut storage = original.storage().to_vec();
            storage.push(MaterialStorage {
                branch: 3,
                element: MaterialStorageElement::Capacitor(dormant),
            });
            let mut circuit = ConductionCircuit::new_with_material_storage(
                original.problem().clone(),
                graph,
                original.terminal().clone(),
                original.material_resistors().to_vec(),
                storage,
                0.0,
                vec![],
                vec![sign * 4.0, 1e12],
                &[],
                field_controls(),
                || false,
            )
            .unwrap();
            // The existing owner permits independent field coordinates in a
            // shared whole caloric volume. This dormant constant-C field has
            // no entropy or capacity contribution as that load cell heats.
            let bindings = [
                active,
                FieldThermalBinding {
                    branch: 3,
                    cell: IDS[1],
                    law: dormant_law,
                    maximum_abs_coordinate: 1e12,
                },
            ];
            let layout = network(&cells);
            let (mut tr0, mut tf0, mut z0) = (300.0, 300.0, 4.0);
            let dt = 0.2;
            for interval in 1..=2 {
                let next = layout
                    .step_field_circuit_candidate(
                        &thermal,
                        interval as f64 * dt,
                        &circuit,
                        &[],
                        &resistors,
                        &bindings,
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0; 3],
                        || false,
                    )
                    .unwrap();
                let (tr, tf, z, current) = oracle(false, alpha, dt, tr0, tf0, z0, 0.0);
                let report = next.report();
                let source = &report.coupled.source;
                let electrical = &report.coupled.electrical;
                let physical_change = potential(alpha, tf, z) - potential(alpha, tf0, z0);
                // Endpoint totals really do erase the active U change. They
                // remain diagnostics, never the enclosing conservation input.
                assert_eq!(
                    source.field_internal_energy_j[1] - source.field_internal_energy_j[0],
                    0.0
                );
                near(source.field_internal_energy_change_j, physical_change);
                near(
                    source.fields[1].state[0].internal_energy_j,
                    source.fields[1].state[1].internal_energy_j,
                );
                assert_eq!(source.fields[1].reversible_heat_j, 0.0);
                near(electrical.spatial_current_a, sign * current);
                near(electrical.heat_j, 3.0 * (tr - tr0));
                near(report.coupled.deposited_joule_j, electrical.heat_j);
                near(source.fields[0].reversible_heat_j, 2.0 * (tf - tf0));
                near(next.thermal().temperature_k(IDS[0]).unwrap(), tr);
                near(next.thermal().temperature_k(IDS[2]).unwrap(), tf);
                near(
                    source.field_internal_energy_change_j
                        + report.coupled.thermal.stored_energy_change_j,
                    0.0,
                );
                near(report.external_heat_j, 0.0);
                near(report.total_energy_residual_j, 0.0);
                // Physical U includes temperature/entropy work; the frozen
                // native coordinate increment of psi is a different quantity.
                assert!(
                    (source.field_internal_energy_change_j - electrical.stored_change_j).abs()
                        > 1e-5
                );
                near(next.circuit().charges_c()[0], sign * z);
                assert_eq!(next.circuit().charges_c()[1], 1e12);
                let (accepted_thermal, accepted_circuit, _) = next.into_parts();
                thermal = accepted_thermal;
                circuit = accepted_circuit;
                (tr0, tf0, z0) = (tr, tf, z);
            }
        }
    }
}

#[test]
fn g1_spatial_rc_and_rl_follow_source_temperature_and_reversible_field_energy() {
    for magnetic in [false, true] {
        let mut outcomes = Vec::new();
        for alpha in [-0.001, 0.0, 0.001] {
            let cards = [
                card(magnetic, 0.0, 1.0),
                card(magnetic, 0.0, 1.0),
                card(magnetic, alpha, 2.0),
            ];
            let curves = cards.each_ref().map(caloric);
            let thermal = history(&curves, 0.0);
            let context = point(300.0);
            let cells = cells(&cards, &curves, &context);
            let resistors = [resistor()];
            let material = MaterialNetwork {
                cells: &cells,
                resistors: &resistors,
                links: &[],
            };
            let (circuit, binding) = initial_circuit(&material, &thermal, &cards, magnetic);
            let bindings = [binding];
            let layout = network(&cells);
            let dt = 0.2;
            let next = layout
                .step_field_circuit_candidate(
                    &thermal,
                    dt,
                    &circuit,
                    &[],
                    &resistors,
                    &bindings,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 3],
                    || false,
                )
                .unwrap();
            let (tr, tf, z, current) = oracle(magnetic, alpha, dt, 300.0, 300.0, 4.0, 0.0);
            let report = next.report();
            let step = &report.coupled.electrical;
            let field = &report.coupled.source.fields[0];
            near(next.thermal().temperature_k(IDS[0]).unwrap(), tr);
            near(next.thermal().temperature_k(IDS[1]).unwrap(), tr);
            near(next.thermal().temperature_k(IDS[2]).unwrap(), tf);
            near(
                if magnetic {
                    next.circuit().fluxes_wb()[0]
                } else {
                    next.circuit().charges_c()[0]
                },
                z,
            );
            near(step.spatial_current_a, current);
            let rho = 1.0 + BETA * (tr - 300.0);
            near(step.spatial_voltage_v, 2.0 * rho * current);
            for i in 0..4 {
                near(step.spatial.current_density_a_m2[i][0], -current);
                near(step.spatial.cell_power_w[i], 0.5 * rho * current * current);
            }
            near(step.graph_resistors[0].power_w, rho * current * current);
            near(
                step.graph_resistors[0]
                    .conductor
                    .as_ref()
                    .unwrap()
                    .resistance_ohm(),
                rho,
            );
            near(step.heat_j, 3.0 * (tr - 300.0));
            near(report.coupled.deposited_joule_j, step.heat_j);
            near(field.reversible_heat_j, 2.0 * (tf - 300.0));
            near(
                field.state[0].internal_energy_j,
                potential(alpha, 300.0, 4.0),
            );
            near(field.state[1].internal_energy_j, potential(alpha, tf, z));
            near(report.external_heat_j, 0.0);
            near(report.total_energy_residual_j, 0.0);
            near(step.native.supplied, 0.0);
            near(
                next.thermal().energy_j() + field.state[1].internal_energy_j,
                thermal.energy_j() + field.state[0].internal_energy_j,
            );
            near(field.material.coefficient(), 1.0 + alpha * (tf - 300.0));
            assert_eq!(
                field.material.material().card_identity(),
                bindings[0].law.material().card_identity()
            );
            if alpha != 0.0 {
                assert_eq!(field.reversible_heat_j.signum(), alpha.signum());
            } else {
                let fixed = layout
                    .step_material_circuit_candidate(
                        &thermal,
                        dt,
                        &circuit,
                        &[],
                        &resistors,
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0; 3],
                        || false,
                    )
                    .unwrap();
                near(fixed.report().electrical.heat_j, step.heat_j);
                near(fixed.thermal().energy_j(), next.thermal().energy_j());
                assert_eq!(field.reversible_heat_j, 0.0);
            }
            let continued = layout
                .step_field_circuit_candidate(
                    next.thermal(),
                    2.0 * dt,
                    next.circuit(),
                    &[],
                    &resistors,
                    &bindings,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 3],
                    || false,
                )
                .unwrap();
            let (tr2, tf2, z2, _) = oracle(magnetic, alpha, dt, tr, tf, z, 0.0);
            near(continued.thermal().temperature_k(IDS[0]).unwrap(), tr2);
            near(continued.thermal().temperature_k(IDS[2]).unwrap(), tf2);
            near(
                if magnetic {
                    continued.circuit().fluxes_wb()[0]
                } else {
                    continued.circuit().charges_c()[0]
                },
                z2,
            );
            let f2 = &continued.report().coupled.source.fields[0];
            near(f2.coordinate[0], z);
            near(
                f2.state[0].internal_energy_j,
                field.state[1].internal_energy_j,
            );
            near(
                continued.thermal().energy_j() + f2.state[1].internal_energy_j,
                thermal.energy_j() + field.state[0].internal_energy_j,
            );
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(circuit.time_s(), 0.0);
            outcomes.push((tf, z));
        }
        assert!(outcomes[0].0 < outcomes[1].0 && outcomes[1].0 < outcomes[2].0);
        assert!((outcomes[2].1 - outcomes[1].1).abs() > 1e-7);
    }
}

#[test]
fn g3_external_heat_and_reversible_heat_remain_distinct() {
    let alpha = -0.001;
    let cards = [
        card(false, 0.0, 1.0),
        card(false, 0.0, 1.0),
        card(false, alpha, 2.0),
    ];
    let curves = cards.each_ref().map(caloric);
    let thermal = history(&curves, 5.0);
    let context = point(300.0);
    let cells = cells(&cards, &curves, &context);
    let resistors = [resistor()];
    let material = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let (circuit, binding) = initial_circuit(&material, &thermal, &cards, false);
    let next = network(&cells)
        .step_field_circuit_candidate(
            &thermal,
            0.2,
            &circuit,
            &[],
            &resistors,
            &[binding],
            limits(),
            field_controls(),
            config(),
            &[0.0; 3],
            || false,
        )
        .unwrap();
    let (tr, tf, z, _) = oracle(false, alpha, 0.2, 300.0, 300.0, 4.0, 1.0);
    near(next.thermal().temperature_k(IDS[0]).unwrap(), tr);
    near(next.thermal().temperature_k(IDS[2]).unwrap(), tf);
    near(next.circuit().charges_c()[0], z);
    let report = next.report();
    near(report.external_heat_j, 1.0);
    near(
        report.coupled.source.reversible_heat_j,
        2.0 * (tf - 300.0) - 1.0,
    );
    near(
        report.coupled.non_electrical_heat_j,
        1.0 + report.coupled.source.reversible_heat_j,
    );
    near(report.total_energy_residual_j, 0.0);
}

#[test]
fn g4_field_admission_budget_and_late_cancellation_preserve_all_histories() {
    let cards = [
        card(false, 0.0, 1.0),
        card(false, 0.0, 1.0),
        card(false, 0.001, 2.0),
    ];
    let curves = cards.each_ref().map(caloric);
    let thermal = history(&curves, 0.0);
    let original_cells = thermal.cells().clone();
    let context = point(300.0);
    let cells = cells(&cards, &curves, &context);
    let resistors = [resistor()];
    let material = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let (circuit, binding) = initial_circuit(&material, &thermal, &cards, false);
    let sources = circuit.storage().to_vec();
    let layout = network(&cells);
    let solve = |fields: &[FieldThermalBinding], controls| {
        layout.step_field_circuit_candidate(
            &thermal,
            0.2,
            &circuit,
            &[],
            &resistors,
            fields,
            limits(),
            field_controls(),
            controls,
            &[0.0; 3],
            || false,
        )
    };
    assert!(solve(&[], config()).is_err());
    for (branch, cell, bound) in [
        (0, IDS[2], 4.0),
        (2, IDS[0], 4.0),
        (2, IDS[2], 3.9),
        (2, IDS[2], 500.0),
    ] {
        let fields = [FieldThermalBinding {
            branch,
            cell,
            law: binding.law.clone(),
            maximum_abs_coordinate: bound,
        }];
        assert!(matches!(
            solve(&fields, config()),
            Err(CircuitEnthalpyError::Source(
                ConductionThermalError::Storage(FieldThermalError::Invalid(_))
            ))
        ));
    }
    let foreign_card = card(false, -0.001, 2.0);
    let foreign = ResolvedCapacitor::resolve(
        &foreign_card,
        &context,
        MaterialPropertySelection::SingleClaimOnly,
        1.0,
        1.0,
    )
    .unwrap();
    let foreign = FieldThermalBinding {
        branch: 2,
        cell: IDS[2],
        maximum_abs_coordinate: 4.0,
        law: ThermalFieldStorage::capacitor(&foreign_card, &foreign, "T", [250.0, 500.0]).unwrap(),
    };
    assert!(solve(&[foreign], config()).is_err());
    let wrong_volume = ResolvedCapacitor::resolve(
        &cards[2],
        &context,
        MaterialPropertySelection::SingleClaimOnly,
        2.0,
        2.0,
    )
    .unwrap();
    let wrong_volume = FieldThermalBinding {
        branch: 2,
        cell: IDS[2],
        maximum_abs_coordinate: 4.0,
        law: ThermalFieldStorage::capacitor(&cards[2], &wrong_volume, "T", [250.0, 500.0]).unwrap(),
    };
    assert!(solve(&[wrong_volume], config()).is_err());
    let bindings = [binding];
    let mut exhausted = config();
    exhausted.maximum_circuit_trials = 1;
    assert!(matches!(
        solve(&bindings, exhausted),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let mut polls = 0;
    let accepted = layout
        .step_field_circuit_candidate(
            &thermal,
            0.2,
            &circuit,
            &[],
            &resistors,
            &bindings,
            limits(),
            field_controls(),
            config(),
            &[0.0; 3],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut visited = 0;
    assert!(matches!(
        layout.step_field_circuit_candidate(
            &thermal,
            0.2,
            &circuit,
            &[],
            &resistors,
            &bindings,
            limits(),
            field_controls(),
            config(),
            &[0.0; 3],
            || {
                visited += 1;
                visited == polls
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(thermal.cells(), &original_cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(circuit.charges_c(), &[4.0]);
    assert_eq!(circuit.storage(), sources.as_slice());
    let replay = solve(&bindings, config()).unwrap();
    assert_eq!(replay.circuit().charges_c(), accepted.circuit().charges_c());
    assert_eq!(replay.thermal().cells(), accepted.thermal().cells());
}

#[test]
fn moving_nonlinear_circuit_keeps_source_field_entropy_and_lorentz_work_distinct() {
    use fs_couple::circuit_enthalpy::{
        conduction::circuit::response::{
            ConductionResponseCircuit, thermal::step_moving_response_field_circuit_candidate,
        },
        contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
    };
    use fs_em::conduction::{
        motion::ConductionMotion,
        response::{ConductionResponseControls, PlanarConductionResponse},
    };
    use fs_material::conductor::response::{ConductorMagnitude, ScalarConductorResponse};

    for (alpha, velocity_y) in [(-0.001, -0.5), (0.001, 0.5)] {
        let cards = [
            card(false, 0.0, 1.0),
            card(false, 0.0, 1.0),
            card(false, alpha, 2.0),
        ];
        let curves = cards.each_ref().map(caloric);
        let mut thermal = history(&curves, 0.0);
        let context = point(300.0)
            .with_quantity(
                "operating",
                ConductorMagnitude::ElectricField.quantity(),
                0.0,
            )
            .unwrap();
        let cells = cells(&cards, &curves, &context);
        let resistors = [resistor()];
        let material = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &[],
        };
        let selection = MaterialPropertySelection::SingleClaimOnly;
        let conductor = resolve_conductor_state_with_basis(
            &cards[0],
            &context,
            selection.clone(),
            ElectricalPropertyBasis::Resistivity,
        )
        .unwrap();
        let source = ScalarConductorResponse::new(
            &cards[0],
            conductor,
            &[("operating", ConductorMagnitude::ElectricField)],
        )
        .unwrap();
        let problem =
            PlanarConductionResponse::new(mesh(1.0), vec![source; 4], limits(), &mut || false)
                .unwrap();
        let capacitor =
            ResolvedCapacitor::resolve(&cards[2], &context, selection, 1.0, 1.0).unwrap();
        let fields = [FieldThermalBinding {
            branch: 2,
            cell: IDS[2],
            law: ThermalFieldStorage::capacitor(&cards[2], &capacitor, "T", [250.0, 500.0])
                .unwrap(),
            maximum_abs_coordinate: 4.0,
        }];
        let velocity = [[0.0, velocity_y]; 4];
        let motion = ConductionMotion {
            velocity_m_s: &velocity,
            normal_flux_density_t: &[1.0; 4],
        };
        let response = ConductionResponseControls {
            maximum_iterations: 100,
            absolute_conductivity_tolerance_s_m: 1e-12,
            relative_conductivity_tolerance: 1e-12,
            current_density_tolerance_a_m2: 1e-10,
        };
        let mut circuit = ConductionResponseCircuit::new_with_motion_and_material_storage(
            problem,
            CircuitGraph {
                node_count: 3,
                branches: vec![
                    (2, 0, Branch::Resistor { ohms: 99.0 }),
                    (1, 2, Branch::Resistor { ohms: 99.0 }),
                    (1, 0, Branch::Capacitor { farads: f64::NAN }),
                ],
                transformers: vec![],
            },
            terminal(),
            material.resolve(&thermal, || false).unwrap().resistors,
            vec![MaterialStorage {
                branch: 2,
                element: MaterialStorageElement::Capacitor(capacitor),
            }],
            0.0,
            vec![],
            vec![4.0],
            &[],
            motion,
            limits(),
            field_controls(),
            response,
            1000,
            || false,
        )
        .unwrap();
        let electrical = ElectricalContactNetwork {
            bulk: &material,
            contacts: &[],
            electrical_bases: &[],
        };
        let heat = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let gate = fs_exec::CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = fs_exec::Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 79,
                    kernel_id: 66,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            let prepared = heat.prepare(&cx, &thermal).unwrap();
            let initial_total = thermal.energy_j() + potential(alpha, 300.0, 4.0);
            let mut mechanical_total = 0.0;
            let mut expected = (300.0, 300.0, 4.0);
            for interval in 1..=2 {
                let target = interval as f64 * 0.2;
                let next = step_moving_response_field_circuit_candidate(
                    &prepared,
                    &cx,
                    &thermal,
                    target,
                    &circuit,
                    &[],
                    motion,
                    &[],
                    &fields,
                    &REGIONS,
                    1e-12,
                    limits(),
                    field_controls(),
                    response,
                    config(),
                    &[0.0; 3],
                )
                .unwrap();
                let ordinary = network(&cells)
                    .step_moving_response_field_circuit_candidate(
                        &thermal,
                        target,
                        &circuit,
                        &[],
                        motion,
                        &resistors,
                        &[],
                        &fields,
                        limits(),
                        field_controls(),
                        response,
                        config(),
                        &[0.0; 3],
                        || false,
                    )
                    .unwrap();
                let report = next.report();
                let step = &report.coupled.electrical;
                let storage = &report.coupled.source.fields[0];
                let t = next.thermal().temperature_k(IDS[2]).unwrap();
                let q = next.circuit().charges_c()[0];
                let j = step.circuit.spatial_current_a;
                let (tr, tf, z, current) = oracle_with_motion(
                    false,
                    alpha,
                    0.2,
                    expected.0,
                    expected.1,
                    expected.2,
                    0.0,
                    2.0 * velocity_y,
                );
                near(next.thermal().temperature_k(IDS[0]).unwrap(), tr);
                near(next.thermal().temperature_k(IDS[1]).unwrap(), tr);
                near(t, tf);
                near(q, z);
                near(j, current);
                expected = (tr, tf, z);
                near(step.mechanical_work_j, 0.2 * 2.0 * velocity_y * j);
                near(
                    step.motion.as_ref().unwrap().mechanical_power_w,
                    2.0 * velocity_y * j,
                );
                near(storage.state[1].internal_energy_j, potential(alpha, t, q));
                near(report.total_energy_residual_j, 0.0);
                near(report.external_heat_j, 0.0);
                near(step.circuit.supplied_work_j, 0.0);
                near(report.coupled.deposited_joule_j, step.circuit.heat_j);
                assert!(storage.reversible_heat_j.abs() > 1e-8);
                for id in IDS {
                    near(
                        next.thermal().temperature_k(id).unwrap(),
                        ordinary.thermal().temperature_k(id).unwrap(),
                    );
                }
                near(q, ordinary.circuit().charges_c()[0]);
                mechanical_total += step.mechanical_work_j;
                near(
                    next.thermal().energy_j()
                        + storage.state[1].internal_energy_j
                        + mechanical_total,
                    initial_total,
                );
                (thermal, circuit, _) = next.into_parts();
            }
        });
    }
}
