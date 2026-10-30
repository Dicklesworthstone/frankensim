//! A source-only coolant substitution changes finite current and stored heat.
use super::*;

#[path = "heat_network/finite_links.rs"]
mod finite_links;
use fs_airflow::graph::thermal::medium::convection::deposition::source::heat::{
    ConvectionHeatError, SourceConvectionHeatNetwork,
};
use fs_couple::circuit_enthalpy::contacts::heat::MaterialHeatNetwork;
use fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork;
use fs_couple::circuit_enthalpy::material::{CellResistor, MaterialCell, MaterialNetwork};
use fs_couple::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpySolve, ResistorCell};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve};
use fs_matdb::SelectionPolicy;
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, ElectricalPropertyBasis};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_solver::NewtonKrylovConfig;

const BODY: CellId = CellId(71);
const SURFACE: f64 = 0.002;
const DT: f64 = 0.05;
const CAPACITANCE: f64 = 0.1;
const INITIAL_CHARGE: f64 = 0.5;

fn body_point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_quantity(), t)
        .unwrap()
        .with("lot", 7.0)
        .unwrap()
}

// The SAME synthetic card supplies caloric storage, inventory and R(T).
fn body_material() -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let cp = PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new(DENSITY_PROPERTY, Density::DIMS);
    let mut claims = ClaimSet::new();
    for (key, values) in [
        (cp.clone(), [CP, CP]),
        (density.clone(), [1000.0, 1000.0]),
        (
            PropertyKey::new("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS),
            [0.6e-4, 1.4e-4],
        ),
    ] {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, values[0]), (320.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature_quantity(), 280.0, 320.0)
                    .with("lot", 7.0, 7.0),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored coolant/electrothermal fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic conductor".into(),
            phase: "solid".into(),
            process: "source substitution test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &cp,
            density: &density,
            temperature_axis: "T",
            lower: &body_point(280.0),
            upper: &body_point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}

fn body_initial(curve: &ResolvedHeatCapacityCurve, htc: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: BODY,
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(0),
                mass_kg: MASS,
                specific_h_j_kg: CP * (310.0 - 280.0),
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [0.1, 0.0, 0.0],
                    volume_m3: 1e-6,
                    surface_area_m2: SURFACE,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: htc,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        }],
        &[(MaterialId(0), curve.curve() as &dyn ThermalLaw)],
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
            max_dt_s: DT,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 1,
        maximum_circuit_trials: 300,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-4,
        cell_energy_tolerance_j: 1e-8,
        thermal_energy_tolerance_j: 1e-8,
        total_energy_tolerance_j: 1e-8,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-9,
            power_tolerance_w: 1e-9,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 1,
            ..NewtonKrylovConfig::default()
        },
    }
}

// Independent flat-plate algebra and exact midpoint RC equation, using only
// this authored fixture's scalar data, never the production correlation/solver.
fn reference(slope: f64) -> (f64, f64, f64, f64) {
    let evaluate = |t: f64| {
        let film = 0.5 * (t + 300.0);
        let d = film - 300.0;
        let rho = 1000.0 - 0.5 * d;
        let cp = 4000.0 + 2.0 * d;
        let k = 0.5 + 0.002 * d;
        let mu = 0.001 * (1.0 - slope * d);
        let h = 0.664 * (rho * 0.1 / mu).sqrt() * (mu * cp / k).cbrt() * k / 0.1;
        let r = 1.0 + 0.02 * (t - 300.0);
        let ratio = DT / (2.0 * CAPACITANCE * r);
        let q = INITIAL_CHARGE * (1.0 - ratio) / (1.0 + ratio);
        let joule = (INITIAL_CHARGE.powi(2) - q.powi(2)) / (2.0 * CAPACITANCE);
        let external = DT * SURFACE * h * (300.0 - t);
        let current = (INITIAL_CHARGE + q) / (2.0 * CAPACITANCE * r);
        (
            MASS * CP * (t - 310.0) - joule - external,
            q,
            current,
            external,
        )
    };
    let (mut lower, mut upper) = (300.0, 312.0);
    for _ in 0..80 {
        let t = 0.5 * (lower + upper);
        if evaluate(t).0 > 0.0 {
            upper = t;
        } else {
            lower = t;
        }
    }
    let t = 0.5 * (lower + upper);
    let (residual, charge, current, external) = evaluate(t);
    assert!(residual.abs() < 1e-10);
    (t, charge, current, external)
}

#[test]
fn coolant_data_reaches_caloric_storage_resistance_current_and_finite_charge() {
    let (card, curve) = body_material();
    let context = body_point(310.0);
    let thermal = body_initial(&curve, 0.0);
    let cells = [MaterialCell {
        cell: BODY,
        material_id: MaterialId(0),
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }];
    let resistors = [CellResistor {
        destination: ResistorCell {
            branch: 1,
            cell: BODY,
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 0.1,
        area_m2: 1e-5,
    }];
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &[],
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[(1, ElectricalPropertyBasis::Resistivity)],
    };
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let circuit = MaterialCircuit::new(
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (
                    1,
                    0,
                    Branch::Capacitor {
                        farads: CAPACITANCE,
                    },
                ),
                (1, 0, Branch::Resistor { ohms: 0.0 }),
            ],
            transformers: vec![],
        },
        bulk.resolve(&thermal, || false).unwrap().resistors,
        vec![],
        vec![INITIAL_CHARGE],
        &[],
    )
    .unwrap();
    let saved = thermal.cells().clone();
    let bindings = [SourceConvectiveCell {
        law: SourceCellConvectionLaw::Forced {
            correlation: CorrelationId::FlatPlateLaminarAverage,
            speed_m_s: 1.0,
            settings: ForcedConvectionSettings::default(),
        },
        ..film(BODY)
    }];
    let gate = fs_exec::CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 171,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let prepared = heat.prepare(&cx, &thermal).unwrap();
        let mut outcomes = Vec::new();
        for slope in [0.01, 0.03] {
            let coolant = liquid_card(0.5, slope, 320.0);
            let law = liquid_law(&coolant);
            let network = SourceConvectionHeatNetwork {
                heat: &prepared,
                cells: &bindings,
                maximum_cells: 1,
                resolve_fluid: |_: &SourceConvectiveCell,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    law.at_temperature(t)
                },
            };
            let candidate = network
                .step_candidate(&cx, &thermal, DT, &circuit, &[], controls(), &[0.0])
                .unwrap();
            let report = candidate.report();
            let source = &report.source.material;
            let temperature = candidate.thermal().temperature_k(BODY).unwrap();
            let expected = reference(slope);
            near(temperature, expected.0, 1e-8);
            near(candidate.circuit().charges_c()[0], expected.1, 1e-9);
            near(report.electrical.resistors[0].current_a, expected.2, 1e-8);
            near(report.non_electrical_heat_j, expected.3, 1e-8);
            let conductor = &report.electrical.resistors[0].material.conductor;
            near(
                conductor.resistance_ohm(),
                1.0 + 0.02 * (temperature - 300.0),
                1e-10,
            );
            assert_eq!(conductor.material().card_identity(), card.content_hash());
            assert_eq!(conductor.material().query().axes()["lot"], 7.0);
            near(
                source.heat.material.cells[0].phase.temperature_k(),
                temperature,
                1e-12,
            );
            let fluid = &source.convection.materials()[0];
            assert_eq!(
                fluid.material.material().card_identity(),
                coolant.content_hash()
            );
            near(
                fluid.property_temperature_k,
                0.5 * (temperature + 300.0),
                1e-12,
            );
            assert_eq!(fluid.material.material().query().axes()["p"], 100000.0);
            assert_eq!(
                fluid.material.material().query().axes()["solute_mass_fraction"],
                0.2
            );
            for property in fluid.material.material().properties() {
                let original = law
                    .reference_state()
                    .material()
                    .property(property.requirement().name())
                    .unwrap();
                assert_eq!(
                    property.answer().receipt.selected,
                    original.answer().receipt.selected
                );
            }
            near(
                report.non_electrical_heat_j,
                DT * source.convection.convection().total_heat_into_cells_w(),
                1e-10,
            );
            near(
                candidate.thermal().energy_j() - thermal.energy_j()
                    + report.electrical.stored_after_j
                    - report.electrical.stored_before_j,
                report.non_electrical_heat_j,
                1e-8,
            );
            near(report.deposited_joule_j, report.electrical.heat_j, 1e-10);
            assert_eq!(candidate.thermal().cells()[&BODY].boundary.htc_w_m2_k, 0.0);
            assert_eq!(candidate.thermal().cells()[&BODY].boundary.source_w, 0.0);
            assert_eq!(thermal.cells(), &saved);
            assert_eq!(circuit.charges_c(), &[INITIAL_CHARGE]);
            outcomes.push((temperature, candidate.circuit().charges_c()[0]));

            // Authored Robin cooling cannot be charged a second time.
            let robin = body_initial(&curve, 2.0);
            assert!(matches!(
                network.resolve(&cx, &robin),
                Err(ConvectionHeatError::Convection(
                    SourceCellConvectionError::Convection(ConvectionSourceError::Invalid(_))
                ))
            ));
            assert_eq!(robin.time_s(), 0.0);
        }
        assert!(outcomes[1].0 < outcomes[0].0);
        assert!(outcomes[1].1 < outcomes[0].1);
        assert!((outcomes[1].1 - outcomes[0].1).abs() > 1e-6);

        let limited_card = liquid_card(0.5, 0.01, 302.0);
        let limited = liquid_law(&limited_card);
        let refused = SourceConvectionHeatNetwork {
            heat: &prepared,
            cells: &bindings,
            maximum_cells: 1,
            resolve_fluid: |_: &SourceConvectiveCell,
                            _: &DepositionSimulation<'_>,
                            t: Temperature| limited.at_temperature(t),
        };
        assert!(matches!(
            refused.step_candidate(&cx, &thermal, DT, &circuit, &[], controls(), &[0.0]),
            Err(CircuitEnthalpyError::Source(
                ConvectionHeatError::Convection(SourceCellConvectionError::Source {
                    cell: BODY,
                    ..
                })
            ))
        ));
        assert_eq!(thermal.cells(), &saved);
        assert_eq!(circuit.charges_c(), &[INITIAL_CHARGE]);
    });
}
