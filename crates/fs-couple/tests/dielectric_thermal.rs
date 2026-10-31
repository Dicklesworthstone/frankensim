//! Original native Debye history coupled to finite caloric stores.
//! Authored numerical fixtures, not measured material properties.
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, ResistorCell,
    dielectric::{self, DebyeThermalBinding},
    material::{CellResistor, MaterialCell, MaterialNetwork},
};
use fs_couple::material_circuit::{MaterialCircuitSolve, dielectric::DebyeCircuit};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    conductor::ELECTRICAL_RESISTIVITY_DIMS,
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS,
        dielectric::{DebyePoleProperties, ResolvedDebyeCapacitor, thermal::ThermalDebyeCapacitor},
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
fn card(alpha: f64, cp: f64, tau_feedback: bool, rho_slope: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut add = |name: &str, dims, knots| {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots,
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
                    source: "authored dielectric/caloric fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    add(
        "epsilon_inf",
        ELECTRIC_PERMITTIVITY_DIMS,
        vec![(250.0, 1.0 - 50.0 * alpha), (500.0, 1.0 + 200.0 * alpha)],
    );
    add(
        "delta_epsilon",
        ELECTRIC_PERMITTIVITY_DIMS,
        vec![(250.0, 0.5 - 50.0 * alpha), (500.0, 0.5 + 200.0 * alpha)],
    );
    add(
        "tau",
        Time::DIMS,
        if tau_feedback {
            vec![(250.0, 0.2), (300.0, 0.2), (300.1, 0.4), (500.0, 0.5)]
        } else {
            vec![(250.0, 0.2), (500.0, 0.2)]
        },
    );
    add(
        "electrical_resistivity",
        ELECTRICAL_RESISTIVITY_DIMS,
        vec![
            (250.0, 1.0 - 50.0 * rho_slope),
            (500.0, 1.0 + 200.0 * rho_slope),
        ],
    );
    add(
        "specific_heat_capacity",
        SPECIFIC_HEAT_CAPACITY_DIMS,
        vec![(250.0, cp), (500.0, cp)],
    );
    add("density", Density::DIMS, vec![(250.0, 1.0), (500.0, 1.0)]);
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
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
fn thermal<'a>(curve: &'a ResolvedHeatCapacityCurve, sources: &[f64]) -> DepositionSimulation<'a> {
    let knots = curve.curve().knots();
    let h0 = knots[0].specific_enthalpy_j_kg
        + (knots[1].specific_enthalpy_j_kg - knots[0].specific_enthalpy_j_kg) / 5.0;
    let events = sources
        .iter()
        .enumerate()
        .map(|(i, &source_w)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(7 + i as u64),
                segment_id: i as u64,
                body_id: i as u64,
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
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MATERIAL, curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: sources.len(),
            max_links: 0,
            max_events: sources.len(),
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
    id: CellId,
    card: &'a MaterialCard,
    curve: &'a ResolvedHeatCapacityCurve,
    context: &'a QueryPoint,
) -> MaterialCell<'a> {
    MaterialCell {
        cell: id,
        material_id: MATERIAL,
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
    instantaneous: f64,
    polarization: Option<f64>,
    resistor: bool,
) -> (DebyeCircuit, ThermalDebyeCapacitor) {
    let poles = if polarization.is_some() {
        vec![DebyePoleProperties {
            permittivity_increment: "delta_epsilon".into(),
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
    let mut branches = vec![(1, 0, Branch::Capacitor { farads: f64::NAN })];
    if resistor {
        branches.push((1, 0, Branch::Resistor { ohms: 1.0 }));
    }
    let circuit = DebyeCircuit::new(
        CircuitGraph {
            node_count: 2,
            branches,
            transformers: vec![],
        },
        0,
        source,
        vec![],
        vec![instantaneous],
        polarization.into_iter().collect(),
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
        maximum_instantaneous_charge_c: 1.0,
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
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

// Independent scalar solution of midpoint polarization and total internal
// energy. Neither production source evaluation nor thermal callbacks are used.
fn internal_energy(alpha: f64, t: f64, qi: f64, qp: f64) -> f64 {
    [(1.0, qi), (0.5, qp)]
        .into_iter()
        .map(|(c0, q)| {
            let c = c0 + alpha * (t - 300.0);
            0.5 * q * q / c * (1.0 + t * alpha / c)
        })
        .sum()
}
fn adiabatic_oracle(alpha: f64, cp: f64, dt: f64) -> (f64, f64) {
    let q = |t| {
        let ci = 1.0 + alpha * (t - 300.0);
        let ck = 0.5 + alpha * (t - 300.0);
        let equilibrium = 0.5 * ck / (ci + ck);
        let z = dt * (ci + ck) / (0.4 * ci);
        equilibrium + (0.5 - equilibrium) * (1.0 - z) / (1.0 + z)
    };
    let energy0 = internal_energy(alpha, 300.0, 0.0, 0.5);
    let residual = |t| cp * (t - 300.0) + internal_energy(alpha, t, 0.5 - q(t), q(t)) - energy0;
    let (mut lo, mut hi) = (299.0, 301.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..60 {
        let t = 0.5 * (lo + hi);
        if residual(t) > 0.0 {
            hi = t;
        } else {
            lo = t;
        }
    }
    let t = 0.5 * (lo + hi);
    (t, q(t))
}

#[test]
fn absorption_heat_and_reversible_transfer_match_independent_total_energy() {
    for alpha in [-0.001, 0.0, 0.001] {
        let card = card(alpha, 4.0, false, 0.0);
        let curve = caloric(&card);
        let thermal = thermal(&curve, &[0.0]);
        let context = point(300.0);
        let cells = [cell(DIELECTRIC, &card, &curve, &context)];
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &[],
        };
        let (circuit, law) = initial(&card, 0.0, Some(0.5), false);
        let saved_state = circuit.state().to_vec();
        let next = dielectric::step_candidate(
            &network,
            &thermal,
            0.08,
            &circuit,
            &[],
            binding(&law, &[1.0]),
            controls(),
            &[0.0],
            || false,
        )
        .unwrap();
        let (expected_t, expected_q) = adiabatic_oracle(alpha, 4.0, 0.08);
        let t = next.thermal().temperature_k(DIELECTRIC).unwrap();
        near(t, expected_t, 1e-9);
        near(
            next.circuit().polarization_charges_c()[0],
            expected_q,
            1e-11,
        );
        near(
            next.circuit().instantaneous_charge_c() + expected_q,
            0.5,
            1e-11,
        );
        near(
            next.circuit().source().instantaneous_capacitance_f(),
            1.0 + alpha * (t - 300.0),
            1e-14,
        );
        let report = next.report();
        let source = &report.coupled.source;
        assert!(t > 300.0 && report.coupled.electrical.polarization_heat_j > 0.0);
        near(
            source.dielectric_state[1].internal_energy_j,
            internal_energy(alpha, t, 0.5 - expected_q, expected_q),
            1e-11,
        );
        near(
            next.thermal().energy_j() - thermal.energy_j()
                + source.dielectric_internal_energy_change_j,
            0.0,
            2e-10,
        );
        near(report.external_heat_j, 0.0, 1e-13);
        near(report.total_energy_residual_j, 0.0, 2e-10);
        near(
            report.coupled.deposited_joule_j,
            report.coupled.electrical.polarization_heat_j,
            1e-12,
        );
        if alpha == 0.0 {
            assert_eq!(source.reversible_heat_j, 0.0);
        } else {
            assert!(source.reversible_heat_j.abs() > 1e-4 && source.parameter_work_j != 0.0);
        }
        let entropy = 4.0 * (t / 300.0).ln() + source.dielectric_state[1].entropy_j_k
            - source.dielectric_state[0].entropy_j_k;
        assert!(entropy > 0.0);
        assert_eq!(circuit.state(), saved_state);
        assert_eq!(thermal.time_s(), 0.0);
    }
}

#[test]
fn zero_pole_external_heat_preserves_charge_and_has_no_polarization_heat() {
    let card = card(0.0, 4.0, false, 0.0);
    let curve = caloric(&card);
    let thermal = thermal(&curve, &[4.0]);
    let context = point(300.0);
    let cells = [cell(DIELECTRIC, &card, &curve, &context)];
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &[],
    };
    let (circuit, law) = initial(&card, 0.8, None, false);
    let next = dielectric::step_candidate(
        &network,
        &thermal,
        0.2,
        &circuit,
        &[],
        binding(&law, &[]),
        controls(),
        &[0.0],
        || false,
    )
    .unwrap();
    near(
        next.thermal().temperature_k(DIELECTRIC).unwrap(),
        300.2,
        1e-10,
    );
    assert_eq!(
        next.circuit().instantaneous_charge_c().to_bits(),
        0.8_f64.to_bits()
    );
    assert_eq!(next.report().coupled.electrical.polarization_heat_j, 0.0);
    assert_eq!(next.report().coupled.source.reversible_heat_j, 0.0);
    near(next.report().external_heat_j, 0.8, 1e-13);
    near(next.report().total_energy_residual_j, 0.0, 1e-10);
}

#[test]
fn temperature_dependent_relaxation_refines_toward_independent_evolution() {
    let card = card(0.0, 4.0, true, 0.0);
    let curve = caloric(&card);
    let context = point(300.0);
    let cells = [cell(DIELECTRIC, &card, &curve, &context)];
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &[],
    };
    // Test-only RK4 integrates the independent scalar ODE with its conserved
    // energy-derived temperature. This is not a second production integrator.
    let rhs = |qp: f64| {
        let qi = 0.5 - qp;
        let t = 300.0 + (0.25 - 0.5 * qi * qi - qp * qp) / 4.0;
        (0.25 - 1.5 * qp) / (0.2 + 2.0 * (t - 300.0))
    };
    let mut reference = 0.5;
    let h = 0.2 / 10000.0;
    for _ in 0..10000 {
        let a = rhs(reference);
        let b = rhs(reference + 0.5 * h * a);
        let c = rhs(reference + 0.5 * h * b);
        let d = rhs(reference + h * c);
        reference += h * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    let run = |steps: usize| {
        let mut thermal = thermal(&curve, &[0.0]);
        let (mut circuit, law) = initial(&card, 0.0, Some(0.5), false);
        for step in 1..=steps {
            let next = dielectric::step_candidate(
                &network,
                &thermal,
                0.2 * step as f64 / steps as f64,
                &circuit,
                &[],
                binding(&law, &[1.0]),
                controls(),
                &[0.0],
                || false,
            )
            .unwrap();
            assert_eq!(next.report().coupled.source.reversible_heat_j, 0.0);
            near(next.report().total_energy_residual_j, 0.0, 2e-10);
            let (t, c, _) = next.into_parts();
            thermal = t;
            circuit = c;
        }
        let t = thermal.temperature_k(DIELECTRIC).unwrap();
        assert!(t > 300.0);
        near(
            circuit.source().poles()[0].relaxation_time_s,
            0.2 + 2.0 * (t - 300.0),
            2e-12,
        );
        near(
            circuit.instantaneous_charge_c() + circuit.polarization_charges_c()[0],
            0.5,
            1e-10,
        );
        (circuit.polarization_charges_c()[0] - reference).abs()
    };
    let errors = [run(20), run(40), run(80)];
    for pair in errors.windows(2) {
        assert!(pair[1] < 0.7 * pair[0], "refinement errors={errors:?}");
    }
    assert!(errors[2] < 3e-4, "fine error={}", errors[2]);
}

#[test]
fn original_resistor_and_polarization_heat_reach_distinct_finite_bodies() {
    let card = card(0.0, 4.0, false, 0.01);
    let curve = caloric(&card);
    let thermal = thermal(&curve, &[0.0, 0.0]);
    let context = point(300.0);
    let cells = [
        cell(DIELECTRIC, &card, &curve, &context),
        cell(RESISTOR, &card, &curve, &context),
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
    let (circuit, law) = initial(&card, 0.6, Some(0.0), true);
    let next = dielectric::step_candidate(
        &network,
        &thermal,
        0.08,
        &circuit,
        &[],
        binding(&law, &[1.0]),
        controls(),
        &[0.0, 0.0],
        || false,
    )
    .unwrap();
    let report = next.report();
    let td = next.thermal().temperature_k(DIELECTRIC).unwrap();
    let tr = next.thermal().temperature_k(RESISTOR).unwrap();
    let electric = &report.coupled.electrical;
    assert!(electric.polarization_heat_j > 0.0 && electric.other_resistor_heat_j > 0.0);
    near(4.0 * (td - 300.0), electric.polarization_heat_j, 2e-10);
    near(4.0 * (tr - 300.0), electric.other_resistor_heat_j, 2e-10);
    assert_eq!(electric.other_resistors.len(), 1);
    assert_eq!(electric.other_resistors[0].branch, 1);
    near(
        electric.other_resistors[0].heat_j,
        electric.other_resistor_heat_j,
        1e-14,
    );
    near(
        next.circuit().material_resistors()[0]
            .conductor
            .resistance_ohm(),
        1.0 + 0.01 * (tr - 300.0),
        1e-13,
    );
    near(
        next.thermal().energy_j() - thermal.energy_j()
            + next.circuit().dielectric_energy_j().unwrap()
            - circuit.dielectric_energy_j().unwrap(),
        0.0,
        2e-10,
    );
    near(report.total_energy_residual_j, 0.0, 2e-10);
}

#[test]
fn capacity_work_budget_and_cancellation_refuse_without_publishing_histories() {
    let card = card(0.001, 4.0, false, 0.0);
    let curve = caloric(&card);
    let thermal = thermal(&curve, &[0.0]);
    let context = point(300.0);
    let cells = [cell(DIELECTRIC, &card, &curve, &context)];
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &[],
    };
    let (circuit, law) = initial(&card, 0.0, Some(0.5), false);
    let saved_state = circuit.state().to_vec();
    let saved_thermal = thermal.cells().clone();
    let mut excess = binding(&law, &[1e6]);
    excess.maximum_instantaneous_charge_c = 1e6;
    assert!(
        dielectric::step_candidate(
            &network,
            &thermal,
            0.08,
            &circuit,
            &[],
            excess,
            controls(),
            &[0.0],
            || false
        )
        .is_err()
    );
    let mut limited = controls();
    limited.maximum_circuit_trials = 1;
    assert!(matches!(
        dielectric::step_candidate(
            &network,
            &thermal,
            0.08,
            &circuit,
            &[],
            binding(&law, &[1.0]),
            limited,
            &[0.0],
            || false
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let mut calls = 0;
    let cancelled = dielectric::step_candidate(
        &network,
        &thermal,
        0.08,
        &circuit,
        &[],
        binding(&law, &[1.0]),
        controls(),
        &[0.0],
        || {
            calls += 1;
            calls > 30
        },
    );
    assert!(matches!(cancelled, Err(CircuitEnthalpyError::Cancelled)));
    assert_eq!(circuit.state(), saved_state);
    assert_eq!(thermal.cells(), &saved_thermal);
    assert_eq!(thermal.time_s(), 0.0);
}
