//! Native nonlinear dielectric memory coupled to finite caloric stores.
//! Authored sources isolate physics and numerics; no measured-data claim.
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, ResistorCell,
    dielectric::{self, DebyeThermalBinding},
    material::{CellResistor, MaterialCell, MaterialNetwork},
};
use fs_couple::material_circuit::{
    MaterialCircuitSolve, capacitance::dielectric::MaterialDielectricCapacitor,
    dielectric::DebyeCircuit,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    conductor::ELECTRICAL_RESISTIVITY_DIMS,
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS,
        dielectric::{
            DebyePoleProperties, ResolvedDebyeCapacitor,
            reversible::{
                ELECTRIC_DISPLACEMENT_DIMS, ELECTRIC_FIELD_DIMS,
                REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY, ReversibleDielectricMaterial,
            },
            thermal::ThermalDebyeCapacitor,
        },
    },
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS},
    thermal::source::ThermalMaterialPoint,
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{
    Density, QuantitySpec, Temperature, Time,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

const DIELECTRIC: CellId = CellId(7);
const RESISTOR: CellId = CellId(8);

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_quantity(), t)
        .unwrap()
}
fn validity(range: [f64; 2]) -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity("T", temperature_quantity(), range[0], range[1])
}
fn fixture(cp: f64, tau_feedback: bool, rho_slope: f64, field_range: [f64; 2]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut add = |name: &str, dims, value, interpolation, validity| {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored nonlinear dielectric/caloric fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    for (name, dims, value) in [
        ("epsilon_inf", ELECTRIC_PERMITTIVITY_DIMS, 2.0),
        ("chi", ELECTRIC_PERMITTIVITY_DIMS, 1.0),
    ] {
        add(
            name,
            dims,
            PropertyValue::Scalar { value, dims },
            InterpolationPolicy::ConstantWithinValidity,
            validity([250.0, 500.0]),
        );
    }
    add(
        REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
        ELECTRIC_DISPLACEMENT_DIMS,
        PropertyValue::Curve {
            abscissa: "E".into(),
            abscissa_dims: ELECTRIC_FIELD_DIMS,
            knots: vec![(0.0, 0.0), (1.0, 2.0), (3.0, 3.0), (4.0, 7.0)],
            dims: ELECTRIC_DISPLACEMENT_DIMS,
        },
        InterpolationPolicy::LinearInside,
        validity(field_range).with_quantity(
            "E",
            QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS),
            0.0,
            4.0,
        ),
    );
    for (name, dims, knots) in [
        (
            "tau",
            Time::DIMS,
            if tau_feedback {
                vec![(250.0, 0.2), (300.0, 0.2), (301.0, 2.2), (500.0, 2.2)]
            } else {
                vec![(250.0, 0.2), (500.0, 0.2)]
            },
        ),
        (
            "electrical_resistivity",
            ELECTRICAL_RESISTIVITY_DIMS,
            vec![
                (250.0, 1.0 - 50.0 * rho_slope),
                (500.0, 1.0 + 200.0 * rho_slope),
            ],
        ),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            vec![(250.0, cp), (500.0, cp)],
        ),
        ("density", Density::DIMS, vec![(250.0, 1.0), (500.0, 1.0)]),
    ] {
        add(
            name,
            dims,
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots,
                dims,
            },
            InterpolationPolicy::LinearInside,
            validity([250.0, 500.0]),
        );
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic nonlinear dielectric".into(),
            phase: "solid".into(),
            process: "finite thermal fixture".into(),
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
fn thermal<'a>(
    curves: &[&'a ResolvedHeatCapacityCurve],
    sources: &[f64],
) -> DepositionSimulation<'a> {
    let events = curves
        .iter()
        .zip(sources)
        .enumerate()
        .map(|(i, (curve, &source_w))| {
            let knots = curve.curve().knots();
            let initial_h = knots[0].specific_enthalpy_j_kg
                + (knots[1].specific_enthalpy_j_kg - knots[0].specific_enthalpy_j_kg) / 5.0;
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(7 + i as u64),
                    segment_id: i as u64,
                    body_id: i as u64,
                    material_id: MaterialId(11 + i as u64),
                    mass_kg: 1.0,
                    specific_h_j_kg: initial_h,
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
            }
        })
        .collect();
    let laws: Vec<_> = curves
        .iter()
        .enumerate()
        .map(|(i, curve)| (MaterialId(11 + i as u64), curve.curve() as &dyn ThermalLaw))
        .collect();
    DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: curves.len(),
            max_links: 0,
            max_events: curves.len(),
            max_substeps: 1000,
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
    index: u64,
    card: &'a MaterialCard,
    curve: &'a ResolvedHeatCapacityCurve,
    context: &'a QueryPoint,
) -> MaterialCell<'a> {
    MaterialCell {
        cell: CellId(7 + index),
        material_id: MaterialId(11 + index),
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }
}
fn initial(
    card: &MaterialCard,
    q: f64,
    p: Option<f64>,
    resistor: bool,
) -> (DebyeCircuit, ThermalDebyeCapacitor) {
    let poles = if p.is_some() {
        vec![DebyePoleProperties {
            permittivity_increment: "chi".into(),
            relaxation_time: "tau".into(),
        }]
    } else {
        vec![]
    };
    let source = ResolvedDebyeCapacitor::resolve(
        card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &poles,
        1.0,
        1.0,
    )
    .unwrap();
    let law = ThermalDebyeCapacitor::new(card, &source, "T", [250.0, 500.0]).unwrap();
    let field_point = point(300.0)
        .with_quantity("E", QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS), 0.0)
        .unwrap();
    let field = ReversibleDielectricMaterial::resolve(
        card,
        &field_point,
        "E",
        4.0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let instantaneous = MaterialDielectricCapacitor::new(0, field, 1.0, 1.0).unwrap();
    let mut branches = vec![(1, 0, Branch::Capacitor { farads: f64::NAN })];
    if resistor {
        branches.push((1, 0, Branch::Resistor { ohms: 1.0 }));
    }
    let circuit = DebyeCircuit::new_with_reversible_instantaneous(
        CircuitGraph {
            node_count: 2,
            branches,
            transformers: vec![],
        },
        source,
        instantaneous,
        vec![],
        vec![q],
        p.into_iter().collect(),
        &[],
        || false,
    )
    .unwrap();
    (circuit, law)
}
fn binding<'a>(law: &'a ThermalDebyeCapacitor, poles: &'a [f64]) -> DebyeThermalBinding<'a> {
    DebyeThermalBinding {
        cell: DIELECTRIC,
        law,
        maximum_instantaneous_charge_c: 7.0,
        maximum_polarization_charges_c: poles,
    }
}
fn controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 2,
        maximum_circuit_trials: 500,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-5,
        cell_energy_tolerance_j: 2e-10,
        thermal_energy_tolerance_j: 2e-10,
        total_energy_tolerance_j: 4e-10,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-12,
            power_tolerance_w: 1e-12,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-13,
            forcing_minimum: 1e-4,
            linear_restart: 4,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}
// Independent inverse and integral of the authored D_inf(E), A=d=1.
fn density(q: f64) -> f64 {
    let q = q.abs();
    if q <= 2.0 {
        q * q / 4.0
    } else if q <= 3.0 {
        let x = q - 2.0;
        1.0 + x + x * x
    } else {
        let x = q - 3.0;
        3.0 + 3.0 * x + x * x / 8.0
    }
}
fn voltage(q: f64) -> f64 {
    let a = q.abs();
    q.signum()
        * if a <= 2.0 {
            a / 2.0
        } else if a <= 3.0 {
            1.0 + 2.0 * (a - 2.0)
        } else {
            3.0 + (a - 3.0) / 4.0
        }
}
fn storage(q: f64, p: f64) -> f64 {
    density(q) + 0.5 * p * p
}
fn tau(t: f64, feedback: bool) -> f64 {
    if feedback {
        0.2 + 2.0 * (t - 300.0)
    } else {
        0.2
    }
}

// One-dimensional solve from conserved total charge, the exact potential
// difference, original pole resistance and the finite caloric energy equation.
fn one_step_oracle(dt: f64, feedback: bool) -> (f64, f64) {
    let initial = storage(1.8, 2.2);
    let temperature = |p| 300.0 + (initial - storage(4.0 - p, p)) / 4.0;
    let residual = |p| {
        let q = 4.0 - p;
        let effort = if q == 1.8 {
            voltage(q)
        } else {
            (density(q) - density(1.8)) / (q - 1.8)
        };
        p - 2.2 - dt * (effort - 0.5 * (p + 2.2)) / tau(temperature(p), feedback)
    };
    let (mut lo, mut hi) = (5.0 / 3.0, 2.2);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..64 {
        let p = 0.5 * (lo + hi);
        if residual(p) > 0.0 {
            hi = p;
        } else {
            lo = p;
        }
    }
    let p = 0.5 * (lo + hi);
    (p, temperature(p))
}

#[test]
fn finite_absorption_heat_refreshes_tau_and_uses_actual_nonlinear_internal_energy() {
    let mut endpoints = Vec::new();
    for feedback in [false, true] {
        for sign in [-1.0, 1.0] {
            let card = fixture(4.0, feedback, 0.0, [250.0, 500.0]);
            let curve = caloric(&card);
            let heat = thermal(&[&curve], &[0.0]);
            let context = point(300.0);
            let cells = [cell(0, &card, &curve, &context)];
            let network = MaterialNetwork {
                cells: &cells,
                resistors: &[],
                links: &[],
            };
            let (circuit, law) = initial(&card, sign * 1.8, Some(sign * 2.2), false);
            let saved = circuit.state().to_vec();
            let next = dielectric::step_candidate(
                &network,
                &heat,
                0.08,
                &circuit,
                &[],
                binding(&law, &[4.0]),
                controls(),
                &[0.0],
                || false,
            )
            .unwrap();
            let (p, t) = one_step_oracle(0.08, feedback);
            let q = 4.0 - p;
            near(next.thermal().temperature_k(DIELECTRIC).unwrap(), t, 2e-10);
            near(next.circuit().polarization_charges_c()[0], sign * p, 2e-11);
            near(next.circuit().instantaneous_charge_c(), sign * q, 2e-11);
            near(
                next.circuit().terminal_charge_c().unwrap(),
                4.0 * sign,
                2e-11,
            );
            let source = &next.report().coupled.source;
            near(
                source.dielectric_state[1].internal_energy_j,
                storage(q, p),
                2e-11,
            );
            near(
                source.dielectric_state[0].free_energy_j,
                storage(1.8, 2.2),
                2e-12,
            );
            assert!(
                (storage(q, p) - (q * q / 4.0 + 0.5 * p * p)).abs() > 1e-4,
                "linear low-field proxy must be falsified"
            );
            for state in source.dielectric_state {
                assert_eq!(state.free_energy_j, state.internal_energy_j);
                assert_eq!(state.entropy_j_k, 0.0);
                assert_eq!(state.heat_capacity_j_k, 0.0);
            }
            assert_eq!(source.parameter_work_j, 0.0);
            assert_eq!(source.thermal_offset_change_j, 0.0);
            assert_eq!(source.reversible_heat_j, 0.0);
            near(
                next.report().coupled.electrical.polarization_heat_j,
                4.0 * (t - 300.0),
                2e-10,
            );
            near(
                next.thermal().energy_j() - heat.energy_j()
                    + source.dielectric_internal_energy_change_j,
                0.0,
                3e-10,
            );
            near(
                next.circuit().source().poles()[0].relaxation_time_s,
                tau(t, feedback),
                2e-11,
            );
            let field = next.circuit().reversible_instantaneous().unwrap();
            near(
                field.material().material().state_coordinate("T").unwrap(),
                t,
                2e-10,
            );
            assert!(std::ptr::eq(
                field.curve().knots().as_ptr(),
                circuit
                    .reversible_instantaneous()
                    .unwrap()
                    .curve()
                    .knots()
                    .as_ptr()
            ));
            near(
                next.report().total_energy_residual_j,
                0.0,
                controls().total_energy_tolerance_j,
            );
            assert_eq!(circuit.state(), saved);
            assert_eq!(heat.time_s(), 0.0);
            if sign > 0.0 {
                endpoints.push(p);
            }
        }
    }
    assert!(
        (endpoints[0] - endpoints[1]).abs() > 1e-2,
        "temperature-dependent relaxation must alter actual charge history"
    );
}

#[test]
fn nonlinear_self_heating_refines_to_independent_continuous_memory_dynamics() {
    let card = fixture(4.0, true, 0.0, [250.0, 500.0]);
    let curve = caloric(&card);
    let context = point(300.0);
    let cells = [cell(0, &card, &curve, &context)];
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &[],
    };
    let energy0 = storage(1.8, 2.2);
    let rhs = |p| {
        let q = 4.0 - p;
        let t = 300.0 + (energy0 - storage(q, p)) / 4.0;
        (voltage(q) - p) / tau(t, true)
    };
    let mut reference = 2.2;
    let h = 0.12 / 16000.0;
    for _ in 0..16000 {
        let a = rhs(reference);
        let b = rhs(reference + 0.5 * h * a);
        let c = rhs(reference + 0.5 * h * b);
        let d = rhs(reference + h * c);
        reference += h * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    let run = |steps: u32| {
        let mut heat = thermal(&[&curve], &[0.0]);
        let (mut circuit, law) = initial(&card, 1.8, Some(2.2), false);
        for step in 1..=steps {
            let next = dielectric::step_candidate(
                &network,
                &heat,
                0.12 * f64::from(step) / f64::from(steps),
                &circuit,
                &[],
                binding(&law, &[4.0]),
                controls(),
                &[0.0],
                || false,
            )
            .unwrap();
            near(next.report().total_energy_residual_j, 0.0, 4e-10);
            let (t, c, _) = next.into_parts();
            heat = t;
            circuit = c;
        }
        let p = circuit.polarization_charges_c()[0];
        near(
            4.0 * (heat.temperature_k(DIELECTRIC).unwrap() - 300.0)
                + storage(circuit.instantaneous_charge_c(), p),
            energy0,
            2e-9,
        );
        assert!(circuit.instantaneous_charge_c() > 2.0);
        (p - reference).abs()
    };
    let errors = [run(20), run(40), run(80)];
    for pair in errors.windows(2) {
        assert!(pair[1] < 0.7 * pair[0], "refinement errors {errors:?}");
    }
    assert!(
        errors[2] < 1e-3,
        "final independent charge error {}",
        errors[2]
    );
}

#[test]
fn zero_pole_external_heating_keeps_nonlinear_charge_energy_and_source_support() {
    let card = fixture(4.0, true, 0.0, [250.0, 500.0]);
    let curve = caloric(&card);
    let heat = thermal(&[&curve], &[2.0]);
    let context = point(300.0);
    let cells = [cell(0, &card, &curve, &context)];
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &[],
    };
    let (circuit, law) = initial(&card, 3.2, None, false);
    let next = dielectric::step_candidate(
        &network,
        &heat,
        0.2,
        &circuit,
        &[],
        binding(&law, &[]),
        controls(),
        &[0.0],
        || false,
    )
    .unwrap();
    let source = &next.report().coupled.source;
    near(
        next.thermal().temperature_k(DIELECTRIC).unwrap(),
        300.1,
        2e-11,
    );
    for state in source.dielectric_state {
        near(state.internal_energy_j, density(3.2), 2e-14);
        assert_eq!(state.entropy_j_k, 0.0);
        assert_eq!(state.heat_capacity_j_k, 0.0);
    }
    assert!((density(3.2) - 3.2 * 3.2 / 4.0).abs() > 1.0);
    assert_eq!(
        next.circuit().instantaneous_charge_c().to_bits(),
        3.2_f64.to_bits()
    );
    assert_eq!(source.dielectric_internal_energy_change_j, 0.0);
    assert_eq!(source.reversible_heat_j, 0.0);
    assert_eq!(next.report().coupled.electrical.polarization_heat_j, 0.0);
    near(next.report().external_heat_j, 0.4, 2e-12);
    near(next.report().total_energy_residual_j, 0.0, 4e-10);
    assert_eq!(
        next.circuit()
            .reversible_instantaneous()
            .unwrap()
            .material()
            .material()
            .card_identity(),
        card.content_hash()
    );
}

#[test]
fn polarization_and_temperature_dependent_conductor_heat_reach_their_original_finite_receivers() {
    let dielectric_card = fixture(4.0, true, 0.0, [250.0, 500.0]);
    let conductor_card = fixture(2.0, false, 0.01, [250.0, 500.0]);
    let a = caloric(&dielectric_card);
    let b = caloric(&conductor_card);
    let heat = thermal(&[&a, &b], &[0.0, 0.0]);
    let context = point(300.0);
    let cells = [
        cell(0, &dielectric_card, &a, &context),
        cell(1, &conductor_card, &b, &context),
    ];
    let resistors = [CellResistor {
        destination: ResistorCell {
            branch: 1,
            cell: RESISTOR,
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    }];
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let (circuit, law) = initial(&dielectric_card, 2.5, Some(3.0), true);
    let next = dielectric::step_candidate(
        &network,
        &heat,
        0.08,
        &circuit,
        &[],
        binding(&law, &[4.0]),
        controls(),
        &[0.0, 0.0],
        || false,
    )
    .unwrap();
    let report = next.report();
    let electrical = &report.coupled.electrical;
    let td = next.thermal().temperature_k(DIELECTRIC).unwrap();
    let tr = next.thermal().temperature_k(RESISTOR).unwrap();
    assert!(td > 300.0 && tr > 300.0 && (td - tr).abs() > 1e-4);
    near(4.0 * (td - 300.0), electrical.polarization_heat_j, 2e-10);
    near(2.0 * (tr - 300.0), electrical.other_resistor_heat_j, 2e-10);
    near(
        next.circuit().source().poles()[0].relaxation_time_s,
        tau(td, true),
        2e-11,
    );
    let conductor = &next.circuit().material_resistors()[0].conductor;
    near(conductor.resistance_ohm(), 1.0 + 0.01 * (tr - 300.0), 2e-13);
    assert_eq!(
        conductor.material().card_identity(),
        conductor_card.content_hash()
    );
    near(
        next.thermal().energy_j() - heat.energy_j()
            + report.coupled.source.dielectric_internal_energy_change_j,
        0.0,
        4e-10,
    );
    near(
        report.coupled.source.dielectric_state[1].internal_energy_j,
        storage(
            next.circuit().instantaneous_charge_c(),
            next.circuit().polarization_charges_c()[0],
        ),
        2e-12,
    );
    near(report.total_energy_residual_j, 0.0, 4e-10);
}

#[test]
fn unsupported_thermal_source_range_charge_envelope_and_cancellation_publish_neither_history() {
    for field_range in [[250.0, 500.0], [290.0, 310.0]] {
        let card = fixture(4.0, true, 0.0, field_range);
        let curve = caloric(&card);
        let heat = thermal(&[&curve], &[0.0]);
        let context = point(300.0);
        let cells = [cell(0, &card, &curve, &context)];
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &[],
        };
        let (circuit, law) = initial(&card, 1.8, Some(2.2), false);
        let saved = circuit.state().to_vec();
        let saved_heat = heat.energy_j();
        if field_range[0] > 250.0 {
            assert!(
                dielectric::step_candidate(
                    &network,
                    &heat,
                    0.01,
                    &circuit,
                    &[],
                    binding(&law, &[4.0]),
                    controls(),
                    &[0.0],
                    || false
                )
                .is_err(),
                "tau/caloric support cannot stand in for the complete original D(E) temperature support"
            );
        } else {
            let mut invalid = binding(&law, &[4.0]);
            invalid.maximum_instantaneous_charge_c = 8.0;
            assert!(
                dielectric::step_candidate(
                    &network,
                    &heat,
                    0.01,
                    &circuit,
                    &[],
                    invalid,
                    controls(),
                    &[0.0],
                    || false
                )
                .is_err()
            );
            let mut calls = 0;
            assert!(matches!(
                dielectric::step_candidate(
                    &network,
                    &heat,
                    0.01,
                    &circuit,
                    &[],
                    binding(&law, &[4.0]),
                    controls(),
                    &[0.0],
                    || {
                        calls += 1;
                        calls > 20
                    }
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
        }
        assert_eq!(circuit.state(), saved);
        assert_eq!(heat.energy_j(), saved_heat);
        assert_eq!(heat.time_s(), 0.0);
        assert_eq!(
            circuit
                .reversible_instantaneous()
                .unwrap()
                .material()
                .material()
                .state_coordinate("T"),
            Some(300.0)
        );
    }
}
