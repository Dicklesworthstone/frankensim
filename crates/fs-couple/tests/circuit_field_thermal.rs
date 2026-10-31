//! Actual material C/L -> native circuit -> original finite thermal residual.
//! Coefficients are explicitly authored numerical fixtures, not measured data.
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, ResistorCell,
    field_storage::{self, FieldThermalBinding},
    material::{CellResistor, MaterialCell, MaterialNetwork},
};
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    conductor::{ELECTRICAL_RESISTIVITY_DIMS, ElectricalPropertyBasis},
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS,
        MAGNETIC_PERMEABILITY_PROPERTY, ResolvedCapacitor, ResolvedInductor,
        thermal::ThermalFieldStorage,
    },
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS},
    thermal::source::ThermalMaterialPoint,
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{
    Density, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

const CELL: CellId = CellId(7);
const MATERIAL: MaterialId = MaterialId(11);
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
    let mut claims = ClaimSet::new();
    let (name, dims) = if magnetic {
        (MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS)
    } else {
        (ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS)
    };
    for (name, dims, low, high) in [
        (name, dims, 1.0 - 50.0 * alpha, 1.0 + 200.0 * alpha),
        (
            "electrical_resistivity",
            ELECTRICAL_RESISTIVITY_DIMS,
            1.0,
            1.0,
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
                    source: "authored affine field/caloric fixture".into(),
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
            process: "numerical fixture".into(),
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
fn thermal(curve: &ResolvedHeatCapacityCurve, source_w: f64) -> DepositionSimulation<'_> {
    let h0 = curve.curve().knots()[0].specific_enthalpy_j_kg
        + 50.0
            * (curve.curve().knots()[1].specific_enthalpy_j_kg
                - curve.curve().knots()[0].specific_enthalpy_j_kg)
            / 250.0;
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELL,
                segment_id: 0,
                body_id: 0,
                material_id: MATERIAL,
                mass_kg: 1.0,
                specific_h_j_kg: h0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w,
                },
            }),
        }],
        &[(MATERIAL, curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
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
fn cell<'a>(
    card: &'a MaterialCard,
    caloric: &'a ResolvedHeatCapacityCurve,
    context: &'a QueryPoint,
) -> MaterialCell<'a> {
    MaterialCell {
        cell: CELL,
        material_id: MATERIAL,
        curve: caloric.curve(),
        point: ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }
}
fn resistor() -> CellResistor {
    CellResistor {
        destination: ResistorCell {
            branch: 1,
            cell: CELL,
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    }
}
fn initial(
    material: &MaterialNetwork<'_>,
    thermal: &DepositionSimulation<'_>,
    card: &MaterialCard,
    magnetic: bool,
) -> (MaterialCircuit, FieldThermalBinding) {
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let (branch, element, law) = if magnetic {
        let value = ResolvedInductor::resolve(card, &point(300.0), selection, 1.0, 1.0, 1).unwrap();
        let law = ThermalFieldStorage::inductor(card, &value, "T", [250.0, 500.0]).unwrap();
        (
            Branch::Inductor { henries: f64::NAN },
            MaterialStorageElement::Inductor(value),
            law,
        )
    } else {
        let value = ResolvedCapacitor::resolve(card, &point(300.0), selection, 1.0, 1.0).unwrap();
        let law = ThermalFieldStorage::capacitor(card, &value, "T", [250.0, 500.0]).unwrap();
        (
            Branch::Capacitor { farads: f64::NAN },
            MaterialStorageElement::Capacitor(value),
            law,
        )
    };
    let circuit = MaterialCircuit::new_with_storage(
        CircuitGraph {
            node_count: 2,
            branches: vec![(1, 0, branch), (1, 0, Branch::Resistor { ohms: f64::NAN })],
            transformers: vec![],
        },
        material.resolve(thermal, || false).unwrap().resistors,
        vec![MaterialStorage { branch: 0, element }],
        if magnetic { vec![10.0] } else { vec![] },
        if magnetic { vec![] } else { vec![10.0] },
        &[],
    )
    .unwrap();
    (
        circuit,
        FieldThermalBinding {
            branch: 0,
            cell: CELL,
            law,
            maximum_abs_coordinate: 10.0,
        },
    )
}
fn controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 1,
        maximum_circuit_trials: 500,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-4,
        cell_energy_tolerance_j: 1e-7,
        thermal_energy_tolerance_j: 1e-7,
        total_energy_tolerance_j: 2e-7,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-9,
            power_tolerance_w: 1e-9,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-14,
            forcing_minimum: 1e-4,
            linear_restart: 4,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-7, "{a} != {b}");
}

#[test]
fn shared_magnetic_storage_cannot_use_independent_branch_field_entropy() {
    use fs_couple::circuit_enthalpy::field_storage::FieldThermalError;
    use fs_phs::circuit::MutualInductance;
    let card = card(true, 0.001, 2.0);
    let caloric = caloric(&card);
    let thermal = thermal(&caloric, 0.0);
    let context = point(300.0);
    let cells = [cell(&card, &caloric, &context)];
    let resistors = [resistor()];
    let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &[] };
    let (base, first) = initial(&material, &thermal, &card, true);
    let mut graph = base.graph().clone();
    graph.branches.push((1, 0, Branch::Inductor { henries: f64::NAN }));
    let mut storage = base.storage().to_vec();
    let mut second = storage[0].clone();
    second.branch = 2;
    let MaterialStorageElement::Inductor(ref coil) = second.element else { unreachable!() };
    let second_binding = FieldThermalBinding { branch: 2, cell: CELL,
        law: ThermalFieldStorage::inductor(&card, coil, "T", [250.0, 500.0]).unwrap(),
        maximum_abs_coordinate: 10.0 };
    storage.push(second);
    let circuit = MaterialCircuit::new_with_storage_and_mutual_inductance(graph,
        base.resistors().to_vec(), storage,
        vec![MutualInductance { branches: (0, 2), henries: 0.5 }],
        vec![1.0, 1.0], vec![], &[]).unwrap();
    let saved = circuit.clone();
    let thermal_before = thermal.cells().clone();
    let result = field_storage::step_candidate(&material, &thermal, 0.2, &circuit,
        &[], &[first, second_binding], &[(1, ElectricalPropertyBasis::Resistivity)],
        controls(), &[0.0], || false);
    assert!(matches!(result, Err(CircuitEnthalpyError::Source(FieldThermalError::Invalid(
        "mutual magnetic thermal feedback requires a shared inductance entropy and capacity law"
    )))));
    assert_eq!(circuit, saved);
    assert_eq!(thermal.cells(), &thermal_before);
}

// Independent closed-form native RC/RL midpoint decay and total-U equation.
// R=C0=L0=1, and Cp=2; no production callback or derivative is called here.
fn oracle(alpha: f64, dt: f64, external_heat: f64) -> (f64, f64) {
    let potential = |t: f64, z: f64| {
        let k = 1.0 + alpha * (t - 300.0);
        0.5 * z * z / k * (1.0 + t * alpha / k)
    };
    let z = |t: f64| {
        let k = 1.0 + alpha * (t - 300.0);
        10.0 * (2.0 * k - dt) / (2.0 * k + dt)
    };
    let f = |t| 2.0 * (t - 300.0) + potential(t, z(t)) - potential(300.0, 10.0) - external_heat;
    let (mut lo, mut hi) = (300.0, 400.0);
    assert!(f(lo) <= 0.0 && f(hi) >= 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let t = 0.5 * (lo + hi);
    (t, z(t))
}

#[test]
fn actual_capacitor_and_inductor_temperature_feedback_matches_total_internal_energy() {
    for magnetic in [false, true] {
        for alpha in [-0.001, 0.0, 0.001] {
            let card = card(magnetic, alpha, 2.0);
            let caloric = caloric(&card);
            let thermal = thermal(&caloric, 0.0);
            let context = point(300.0);
            let cells = [cell(&card, &caloric, &context)];
            let resistors = [resistor()];
            let material = MaterialNetwork {
                cells: &cells,
                resistors: &resistors,
                links: &[],
            };
            let (circuit, field) = initial(&material, &thermal, &card, magnetic);
            let bindings = [field];
            let bases = [(1, ElectricalPropertyBasis::Resistivity)];
            let next = field_storage::step_candidate(
                &material,
                &thermal,
                0.2,
                &circuit,
                &[],
                &bindings,
                &bases,
                controls(),
                &[0.0],
                || false,
            )
            .unwrap();
            let (expected_t, expected_z) = oracle(alpha, 0.2, 0.0);
            near(next.thermal().temperature_k(CELL).unwrap(), expected_t);
            near(
                if magnetic {
                    next.circuit().fluxes_wb()[0]
                } else {
                    next.circuit().charges_c()[0]
                },
                expected_z,
            );
            let report = next.report();
            near(report.external_heat_j, 0.0);
            near(report.total_energy_residual_j, 0.0);
            let source = &report.coupled.source;
            near(
                source.field_internal_energy_j[1] - source.field_internal_energy_j[0]
                    + next.thermal().energy_j()
                    - thermal.energy_j(),
                0.0,
            );
            assert!(report.coupled.deposited_joule_j > 0.0);
            if alpha == 0.0 {
                assert_eq!(source.reversible_heat_j, 0.0);
            } else {
                assert!(source.reversible_heat_j.abs() > 0.1);
            }
            // Closed entropy must increase: independent constant-Cp chart plus
            // the field entropy. Joule heating is its irreversible source.
            let ds = 2.0 * (expected_t / 300.0).ln() + source.fields[0].state[1].entropy_j_k
                - source.fields[0].state[0].entropy_j_k;
            assert!(ds > 0.0);
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(
                if magnetic {
                    circuit.fluxes_wb()
                } else {
                    circuit.charges_c()
                },
                &[10.0]
            );
        }
    }
}

#[test]
fn authored_heat_and_repeated_steps_retain_field_caloric_and_source_identity() {
    let card = card(false, 0.001, 2.0);
    let caloric = caloric(&card);
    let mut thermal = thermal(&caloric, 2.0);
    let context = point(300.0);
    let cells = [cell(&card, &caloric, &context)];
    let resistors = [resistor()];
    let material = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let (mut circuit, binding) = initial(&material, &thermal, &card, false);
    let bindings = [binding];
    let bases = [(1, ElectricalPropertyBasis::Resistivity)];
    let initial_u = thermal.energy_j() + 65.0;
    for n in 1..=10 {
        let solve = || {
            field_storage::step_candidate(
                &material,
                &thermal,
                n as f64 * 0.1,
                &circuit,
                &[],
                &bindings,
                &bases,
                controls(),
                &[0.0],
                || false,
            )
            .unwrap()
        };
        let next = solve();
        let replay = solve();
        assert_eq!(next.circuit(), replay.circuit());
        assert_eq!(next.thermal().cells(), replay.thermal().cells());
        near(next.report().external_heat_j, 0.2);
        near(
            next.thermal().energy_j() + next.report().coupled.source.field_internal_energy_j[1],
            initial_u + 0.2 * n as f64,
        );
        let material = match &next.circuit().storage()[0].element {
            MaterialStorageElement::Capacitor(c) => c.material(),
            _ => panic!("capacitor expected"),
        };
        assert_eq!(material.card_identity(), card.content_hash());
        near(
            material.query().axes()["T"],
            next.thermal().temperature_k(CELL).unwrap(),
        );
        let (t, c, _) = next.into_parts();
        thermal = t;
        circuit = c;
    }
}

#[test]
fn invalid_capacity_mapping_source_context_and_late_failure_leave_both_histories() {
    let card = card(false, 0.001, 2.0);
    let caloric = caloric(&card);
    let thermal = thermal(&caloric, 0.0);
    let context = point(300.0);
    let cells = [cell(&card, &caloric, &context)];
    let resistors = [resistor()];
    let material = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let (circuit, binding) = initial(&material, &thermal, &card, false);
    let mut bindings = [binding];
    let bases = [(1, ElectricalPropertyBasis::Resistivity)];
    let saved = thermal.cells().clone();
    let old = circuit.clone();
    let solve = |fields: &[FieldThermalBinding], config: CircuitEnthalpySolve| {
        field_storage::step_candidate(
            &material,
            &thermal,
            0.2,
            &circuit,
            &[],
            fields,
            &bases,
            config,
            &[0.0],
            || false,
        )
    };
    assert!(solve(&[], controls()).is_err());
    bindings[0].maximum_abs_coordinate = 9.0;
    assert!(solve(&bindings, controls()).is_err());
    bindings[0].maximum_abs_coordinate = 1000.0; // negative total heat capacity envelope
    assert!(solve(&bindings, controls()).is_err());
    bindings[0].maximum_abs_coordinate = 10.0;
    bindings[0].cell = CellId(99);
    assert!(solve(&bindings, controls()).is_err());
    bindings[0].cell = CELL;
    // Equal field coefficients and the same material label do not authorize
    // replacing the immutable caloric/field source or its represented volume.
    let original_law = bindings[0].law.clone();
    let foreign_card = self::card(false, 0.001, 3.0);
    let foreign_capacitor = ResolvedCapacitor::resolve(
        &foreign_card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        1.0,
        1.0,
    )
    .unwrap();
    bindings[0].law =
        ThermalFieldStorage::capacitor(&foreign_card, &foreign_capacitor, "T", [250.0, 500.0])
            .unwrap();
    assert!(solve(&bindings, controls()).is_err());
    // Scaling both facing area and gap preserves C, but quadruples volume.
    let different_volume = ResolvedCapacitor::resolve(
        &card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        2.0,
        2.0,
    )
    .unwrap();
    bindings[0].law =
        ThermalFieldStorage::capacitor(&card, &different_volume, "T", [250.0, 500.0]).unwrap();
    assert!(solve(&bindings, controls()).is_err());
    bindings[0].law = original_law;
    let mut polls = 0;
    let next = field_storage::step_candidate(
        &material,
        &thermal,
        0.2,
        &circuit,
        &[],
        &bindings,
        &bases,
        controls(),
        &[0.0],
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut limited = controls();
    limited.maximum_circuit_trials = next.report().coupled.circuit_trials - 1;
    assert!(matches!(
        solve(&bindings, limited),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let mut count = 0;
    assert!(matches!(
        field_storage::step_candidate(
            &material,
            &thermal,
            0.2,
            &circuit,
            &[],
            &bindings,
            &bases,
            controls(),
            &[0.0],
            || {
                count += 1;
                count == polls
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(thermal.cells(), &saved);
    assert_eq!(circuit, old);
}


#[test]
fn multiple_field_branches_share_one_material_cell_without_duplicating_caloric_storage() {
    let card = card(false, 0.001, 2.0);
    let caloric = caloric(&card);
    let thermal = thermal(&caloric, 0.0);
    let context = point(300.0);
    let cells = [cell(&card, &caloric, &context)];
    let resistor = CellResistor {
        destination: ResistorCell { branch: 2, cell: CELL },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    };
    let resistors = [resistor];
    let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &[] };
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let first =
        ResolvedCapacitor::resolve(&card, &point(300.0), selection.clone(), 1.0, 1.0).unwrap();
    let second = ResolvedCapacitor::resolve(&card, &point(300.0), selection, 1.0, 1.0).unwrap();
    let law = ThermalFieldStorage::capacitor(&card, &first, "T", [250.0, 500.0]).unwrap();
    let circuit = MaterialCircuit::new_with_storage(
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: f64::NAN }),
                (1, 0, Branch::Capacitor { farads: f64::NAN }),
                (1, 0, Branch::Resistor { ohms: f64::NAN }),
            ],
            transformers: vec![],
        },
        material.resolve(&thermal, || false).unwrap().resistors,
        vec![
            MaterialStorage { branch: 0, element: MaterialStorageElement::Capacitor(first) },
            MaterialStorage { branch: 1, element: MaterialStorageElement::Capacitor(second) },
        ],
        vec![],
        vec![10.0, 6.0],
        &[],
    ).unwrap();
    let fields = [
        FieldThermalBinding { branch: 0, cell: CELL, law: law.clone(), maximum_abs_coordinate: 10.0 },
        FieldThermalBinding { branch: 1, cell: CELL, law, maximum_abs_coordinate: 6.0 },
    ];
    let bases = [(2, ElectricalPropertyBasis::Resistivity)];
    let next = field_storage::step_candidate(
        &material, &thermal, 0.2, &circuit, &[], &fields, &bases,
        controls(), &[0.0], || false,
    ).unwrap();

    assert_eq!(next.report().coupled.source.fields.len(), 2);
    assert!(next.report().coupled.source.fields.iter().all(|field| field.cell == CELL));
    near(next.report().total_energy_residual_j, 0.0);
    near(next.report().external_heat_j, 0.0);
    near(
        next.thermal().energy_j() + next.report().coupled.source.field_internal_energy_j[1],
        thermal.energy_j() + next.report().coupled.source.field_internal_energy_j[0],
    );
    assert_ne!(next.circuit().charges_c(), circuit.charges_c());
    near(
        next.report().coupled.source.reversible_heat_j,
        next.report().coupled.source.fields.iter().map(|field| field.reversible_heat_j).sum(),
    );

    // Branch identity remains unique even though cell ownership may be shared.
    let duplicate_branch = [
        FieldThermalBinding { branch: 0, cell: CELL, law: fields[0].law.clone(), maximum_abs_coordinate: 10.0 },
        FieldThermalBinding { branch: 0, cell: CELL, law: fields[1].law.clone(), maximum_abs_coordinate: 6.0 },
    ];
    assert!(field_storage::step_candidate(
        &material, &thermal, 0.2, &circuit, &[], &duplicate_branch, &bases,
        controls(), &[0.0], || false,
    ).is_err());
}
