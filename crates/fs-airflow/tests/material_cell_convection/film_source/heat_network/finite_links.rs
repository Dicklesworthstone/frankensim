//! One source supplies a finite coolant's heat capacity and transport.
//! Its native thermal cell and the original RC circuit close one energy balance.
use super::*;
#[path = "finite_links/observer.rs"]
mod observer;
#[path = "finite_links/field_storage.rs"]
mod field_storage;
#[path = "finite_links/piezoelectric.rs"]
mod piezoelectric;
#[path = "finite_links/piezoelectric_thermal.rs"]
mod piezoelectric_thermal;
#[path = "finite_links/thermoelectric.rs"]
mod thermoelectric;
#[path = "finite_links/graph.rs"]
mod graph;
#[path = "finite_links/hysteretic.rs"]
mod hysteretic;
#[path = "finite_links/spatial_induction.rs"]
mod spatial_induction;
#[path = "finite_links/magnetostrictive.rs"]
mod magnetostrictive;
use fs_airflow::graph::thermal::medium::convection::deposition::source::links::{
    SourceConvectiveLink, evaluate_source_link_convection,
};
use fs_couple::circuit_enthalpy::material::ConductionLink;

const FLUID: CellId = CellId(72);
const LINK: u64 = 17;

fn fluid_point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_quantity(), t)
        .unwrap()
        .with("p", 100000.0)
        .unwrap()
        .with("solute_mass_fraction", 0.2)
        .unwrap()
}

fn fluid_material(
    viscosity_scale: f64,
    cp_value: f64,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let original = liquid_card(0.0, 0.01, 320.0);
    let mut claims = ClaimSet::new();
    for (_, original) in original.claims().claims_ordered() {
        let mut claim = original.clone();
        if claim.key.name() == SPECIFIC_HEAT_CAPACITY_PROPERTY {
            claim.value = PropertyValue::Scalar {
                value: cp_value,
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            };
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        if claim.key.name() == DYNAMIC_VISCOSITY_PROPERTY {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                for (_, value) in knots {
                    *value *= viscosity_scale;
                }
            }
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &fluid_point(280.0),
            upper: &fluid_point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}

fn initial_pair<'a>(
    body: &'a ResolvedHeatCapacityCurve,
    fluid: &'a ResolvedHeatCapacityCurve,
    fluid_cp: f64,
    reversed: bool,
    authored_g: f64,
) -> DepositionSimulation<'a> {
    let mut events = [(BODY, 0, CP * 30.0), (FLUID, 1, fluid_cp * 10.0)]
        .map(|(id, material, specific_h_j_kg)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: material,
                body_id: material,
                material_id: MaterialId(material),
                mass_kg: MASS,
                specific_h_j_kg,
                geometry: CellGeometry {
                    start_m: [0.0, material as f64, 0.0],
                    end_m: [0.1, material as f64, 0.0],
                    volume_m3: 1e-6,
                    surface_area_m2: SURFACE,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .to_vec();
    let [a, b] = if reversed {
        [FLUID, BODY]
    } else {
        [BODY, FLUID]
    };
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: LINK,
            a,
            b,
            area_m2: SURFACE,
            conductance_w_k: authored_g,
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), body.curve() as &dyn ThermalLaw),
            (MaterialId(1), fluid.curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
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

fn binding() -> SourceConvectiveLink {
    SourceConvectiveLink {
        link: LINK,
        fluid_cell: FLUID,
        characteristic_length: Length::new(0.1),
        property_temperature: ConvectionPropertyTemperature::FilmMean,
        law: SourceCellConvectionLaw::Forced {
            correlation: CorrelationId::FlatPlateLaminarAverage,
            speed_m_s: 1.0,
            settings: ForcedConvectionSettings::default(),
        },
    }
}

// Independent algebra: midpoint capacitor loss is the only heat supply.
// Conservation eliminates the fluid temperature, leaving one scalar equation.
fn reference_pair(viscosity_scale: f64, cp: f64) -> [f64; 5] {
    let evaluate = |wall: f64| {
        let resistance = 1.0 + 0.02 * (wall - 300.0);
        let ratio = DT / (2.0 * CAPACITANCE * resistance);
        let charge = INITIAL_CHARGE * (1.0 - ratio) / (1.0 + ratio);
        let joule = (INITIAL_CHARGE.powi(2) - charge.powi(2)) / (2.0 * CAPACITANCE);
        let fluid = 290.0 + (joule - MASS * CP * (wall - 310.0)) / (MASS * cp);
        let film = 0.5 * (wall + fluid);
        let k = 0.5 + 0.002 * (film - 300.0);
        let mu = viscosity_scale * 0.001 * (1.0 - 0.01 * (film - 300.0));
        let h = 0.664 * (1000.0 * 0.1 / mu).sqrt() * (mu * cp / k).cbrt() * k / 0.1;
        let power = SURFACE * h * (wall - fluid);
        let residual = MASS * CP * (wall - 310.0) - joule + DT * power;
        let current = (INITIAL_CHARGE + charge) / (2.0 * CAPACITANCE * resistance);
        (residual, [wall, fluid, charge, current, power])
    };
    let (mut lower, mut upper) = (290.0, 312.0);
    for _ in 0..80 {
        let wall = 0.5 * (lower + upper);
        if evaluate(wall).0 > 0.0 {
            upper = wall;
        } else {
            lower = wall;
        }
    }
    let (residual, result) = evaluate(0.5 * (lower + upper));
    assert!(residual.abs() < 1e-10);
    result
}

#[test]
fn g1_g3_finite_coolant_storage_and_transport_reach_the_original_circuit() {
    let run = |viscosity_scale, fluid_cp, reversed| {
        let (body_card, body_curve) = body_material();
        let (fluid_card, fluid_curve) = fluid_material(viscosity_scale, fluid_cp);
        let body_context = body_point(310.0);
        let fluid_context = fluid_point(290.0);
        let thermal = initial_pair(&body_curve, &fluid_curve, fluid_cp, reversed, 0.0);
        let saved = thermal.cells().clone();
        let saved_links = thermal.links().clone();
        let fluid = liquid_law(&fluid_card);
        let cells = [
            MaterialCell {
                cell: BODY,
                material_id: MaterialId(0),
                curve: body_curve.curve(),
                point: ThermalMaterialPoint {
                    material: &body_card,
                    context: &body_context,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                conductivity: None,
            },
            MaterialCell {
                cell: FLUID,
                material_id: MaterialId(1),
                curve: fluid_curve.curve(),
                point: ThermalMaterialPoint {
                    material: &fluid_card,
                    context: &fluid_context,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                conductivity: None,
            },
        ];
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
        let links = [binding()];
        let initial = evaluate_source_link_convection(
            &thermal,
            &links,
            1,
            |_, _, t| fluid.at_temperature(t),
            || false,
        )
        .unwrap();
        let material = bulk
            .prepare_with_link_conductances(
                &thermal,
                electrical.electrical_bases,
                initial.conductances_w_k(),
                || false,
            )
            .unwrap();
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
            material
                .resolve_with_link_conductances(&thermal, initial.conductances_w_k(), || false)
                .unwrap()
                .resistors,
            vec![],
            vec![INITIAL_CHARGE],
            &[],
        )
        .unwrap();
        let gate = fs_exec::CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let mut accepted = [0.0; 2];
        pool.scope(|arena| {
            let cx = fs_exec::Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 173,
                    kernel_id: 29,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            // Missing a physical link owner must never silently mean insulation.
            assert!(heat.prepare(&cx, &thermal).is_err());
            let prepared = heat
                .prepare_with_link_conductances(&cx, &thermal, initial.conductances_w_k())
                .unwrap();
            assert!(prepared.resolve(&cx, &thermal).is_err());
            assert!(
                prepared
                    .resolve_with_link_conductances(&cx, &thermal, &BTreeMap::new())
                    .is_err()
            );
            for map in [
                BTreeMap::from([(LINK, f64::NAN)]),
                BTreeMap::from([(LINK + 1, 1.0)]),
            ] {
                assert!(
                    prepared
                        .resolve_with_link_conductances(&cx, &thermal, &map)
                        .is_err()
                );
            }
            let conflicting = MaterialNetwork {
                cells: &cells,
                resistors: &resistors,
                links: &[ConductionLink {
                    link: LINK,
                    distance_a_m: 0.01,
                    distance_b_m: 0.01,
                }],
            };
            assert!(
                conflicting
                    .prepare_with_link_conductances(
                        &thermal,
                        electrical.electrical_bases,
                        initial.conductances_w_k(),
                        || false,
                    )
                    .is_err()
            );

            let network = SourceConvectionHeatNetwork {
                heat: &prepared,
                cells: &[],
                maximum_cells: 0,
                resolve_fluid: |_: &SourceConvectiveCell,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    fluid.at_temperature(t)
                },
            };
            let mut solve = controls();
            solve.maximum_cells = 2;
            solve.newton.linear_restart = 2;
            let next = network
                .step_candidate_with_links(
                    &cx,
                    &thermal,
                    DT,
                    &circuit,
                    &[],
                    solve,
                    &[0.0, 0.0],
                    &links,
                    1,
                    |_, _, t| fluid.at_temperature(t),
                )
                .unwrap();
            let expected = reference_pair(viscosity_scale, fluid_cp);
            let report = next.report();
            near(
                next.thermal().temperature_k(BODY).unwrap(),
                expected[0],
                1e-8,
            );
            near(
                next.thermal().temperature_k(FLUID).unwrap(),
                expected[1],
                1e-8,
            );
            near(next.circuit().charges_c()[0], expected[2], 1e-9);
            near(report.electrical.resistors[0].current_a, expected[3], 1e-8);
            let source = &report.source.material;
            let finite = source.finite_convection.as_ref().unwrap();
            let sample = &finite.samples()[0];
            assert_eq!(sample.wall_cell, BODY);
            assert_eq!(sample.fluid_cell, FLUID);
            near(sample.heat_into_fluid_w, expected[4], 1e-7);
            near(
                sample.property_temperature_k,
                0.5 * (next.thermal().temperature_k(BODY).unwrap()
                    + next.thermal().temperature_k(FLUID).unwrap()),
                1e-12,
            );
            assert_eq!(
                sample.material.material().card_identity(),
                fluid_card.content_hash()
            );
            assert_eq!(sample.material.material().query().axes()["p"], 100000.0);
            assert_eq!(
                sample.material.material().query().axes()["solute_mass_fraction"],
                0.2
            );
            for property in sample.material.material().properties() {
                let original = fluid
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
                MASS * fluid_cp * (expected[1] - 290.0),
                DT * sample.heat_into_fluid_w,
                1e-8,
            );
            near(
                MASS * CP * (expected[0] - 310.0),
                report.deposited_joule_j - DT * sample.heat_into_fluid_w,
                1e-8,
            );
            near(report.deposited_joule_j, report.electrical.heat_j, 1e-10);
            near(report.non_electrical_heat_j, 0.0, 1e-12);
            near(
                source.convection.convection().total_heat_into_cells_w(),
                0.0,
                1e-12,
            );
            near(
                next.thermal().energy_j() - thermal.energy_j() + report.electrical.stored_after_j
                    - report.electrical.stored_before_j,
                0.0,
                1e-8,
            );
            assert!(source.heat.material.links.is_empty());
            for sample in &source.heat.material.cells {
                near(
                    sample.phase.temperature_k(),
                    next.thermal().temperature_k(sample.cell).unwrap(),
                    1e-12,
                );
            }
            assert_eq!(next.thermal().links(), &saved_links);
            assert_eq!(thermal.cells(), &saved);
            assert_eq!(circuit.charges_c(), &[INITIAL_CHARGE]);
            accepted = [
                next.thermal().temperature_k(BODY).unwrap(),
                next.circuit().charges_c()[0],
            ];

            let authored = initial_pair(&body_curve, &fluid_curve, fluid_cp, reversed, 1.0);
            assert!(
                network
                    .resolve_with_links(&cx, &authored, &links, 1, |_, _, t| fluid
                        .at_temperature(t),)
                    .is_err()
            );
            let refused = SourceConvectionHeatNetwork {
                heat: &prepared,
                cells: &[],
                maximum_cells: 0,
                resolve_fluid: |_: &SourceConvectiveCell,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    Ok::<_, &'static str>(fluid.at_temperature(t).unwrap())
                },
            };
            assert!(matches!(
                refused.step_candidate_with_links(
                    &cx,
                    &thermal,
                    DT,
                    &circuit,
                    &[],
                    solve,
                    &[0.0, 0.0],
                    &links,
                    1,
                    |_, _, _| Err::<fs_material::liquid::ResolvedLiquidState, _>(
                        "fluid source refused"
                    ),
                ),
                Err(CircuitEnthalpyError::Source(
                    ConvectionHeatError::Convection(SourceCellConvectionError::Source {
                        cell: FLUID,
                        source: "fluid source refused",
                        ..
                    })
                ))
            ));
            assert_eq!(thermal.cells(), &saved);
            assert_eq!(thermal.links(), &saved_links);
            assert_eq!(circuit.charges_c(), &[INITIAL_CHARGE]);
        });
        accepted
    };
    let reference = run(1.0, 4000.0, false);
    let reversed = run(1.0, 4000.0, true);
    for (actual, expected) in reversed.into_iter().zip(reference) {
        near(actual, expected, 1e-9);
    }
    let viscous = run(4.0, 4000.0, false);
    let higher_capacity = run(1.0, 8000.0, false);
    assert!(viscous[0] > reference[0] && viscous[1] > reference[1]);
    assert!(higher_capacity[0] < reference[0] && higher_capacity[1] < reference[1]);
}

#[path = "finite_links/solid.rs"]
mod solid;
