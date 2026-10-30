//! Finite material coolant reaches Kelvin conversion and continued charge.
use super::*;
use fs_couple::circuit_enthalpy::contacts::heat::thermoelectric::step_thermoelectric_candidate;
use fs_couple::circuit_enthalpy::thermoelectric::{ThermoelectricLeg, ThermoelectricNetworkError};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;
use fs_material::thermoelectric::{ABSOLUTE_SEEBECK_PROPERTY, SEEBECK_DIMS};

const HOT: CellId = CellId(69);
const COLD: CellId = CellId(70);
const SOLID_K: f64 = 100.0;
const SEEBECK: f64 = 0.02;
const SLOPE: f64 = 4e-4;
const INITIAL: [f64; 4] = [310.0, 290.0, 300.0, 290.0];

fn kelvin_material() -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (base, _) = body_material();
    let mut claims = base.claims().clone();
    for (key, value, interpolation) in [
        (
            PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            PropertyValue::Scalar {
                value: SOLID_K,
                dims: THERMAL_CONDUCTIVITY_DIMS,
            },
            InterpolationPolicy::ConstantWithinValidity,
        ),
        (
            PropertyKey::new(ABSOLUTE_SEEBECK_PROPERTY, SEEBECK_DIMS),
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (280.0, SEEBECK - 20.0 * SLOPE),
                    (320.0, SEEBECK + 20.0 * SLOPE),
                ],
                dims: SEEBECK_DIMS,
            },
            InterpolationPolicy::LinearInside,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key,
                value,
                interpolation,
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature_quantity(), 280.0, 320.0)
                    .with("lot", 7.0, 7.0),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored finite-coolant Kelvin fixture, not measured alloy data"
                        .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
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

fn initial<'a>(
    solid: &'a ResolvedHeatCapacityCurve,
    fluid: &'a ResolvedHeatCapacityCurve,
    fluid_cp: f64,
) -> DepositionSimulation<'a> {
    let mut events = [HOT, COLD, BODY, FLUID]
        .into_iter()
        .enumerate()
        .map(|(i, id)| {
            let liquid = id == FLUID;
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: i as u64,
                    body_id: i as u64,
                    material_id: MaterialId(if liquid { 1 } else { 0 }),
                    mass_kg: MASS,
                    specific_h_j_kg: if liquid { fluid_cp } else { CP } * (INITIAL[i] - 280.0),
                    geometry: CellGeometry {
                        start_m: [i as f64, 0.0, 0.0],
                        end_m: [i as f64 + 0.1, 0.0, 0.0],
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
            }
        })
        .collect::<Vec<_>>();
    for (id, a, b, area) in [
        (0, HOT, BODY, 1e-5),
        (1, BODY, COLD, 1e-5),
        (LINK, HOT, FLUID, SURFACE),
    ] {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id,
                a,
                b,
                area_m2: area,
                conductance_w_k: 0.0,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), solid.curve() as &dyn ThermalLaw),
            (MaterialId(1), fluid.curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 4,
            max_links: 3,
            max_events: 7,
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

// Independent endpoint heat equations, affine-S integrals and midpoint charge
// balance. No production correlation, constitutive query or residual is reused.
fn reference(
    before: [f64; 4],
    after: [f64; 4],
    charge: f64,
    viscosity_scale: f64,
    fluid_cp: f64,
) -> (f64, f64, [f64; 3], f64) {
    let [hot, cold, bulk, fluid] = after;
    let resistance = 1.0 + 0.02 * (bulk - 300.0);
    let emf = (SEEBECK + SLOPE * (0.5 * (hot + cold) - 300.0)) * (hot - cold);
    let current = (charge / CAPACITANCE + emf) / (resistance + DT / (2.0 * CAPACITANCE));
    let next_charge = charge - DT * current;
    let left = SEEBECK + SLOPE * (0.5 * (hot + bulk) - 300.0);
    let right = SEEBECK + SLOPE * (0.5 * (bulk + cold) - 300.0);
    let kelvin = [
        -current * left * hot,
        current * right * cold,
        current * (left - right) * bulk,
    ];
    let film = 0.5 * (hot + fluid);
    let k = 0.5 + 0.002 * (film - 300.0);
    let mu = viscosity_scale * 0.001 * (1.0 - 0.01 * (film - 300.0));
    let h = 0.664 * (1000.0 * 0.1 / mu).sqrt() * (mu * fluid_cp / k).cbrt() * k / 0.1;
    let convection = SURFACE * h * (hot - fluid);
    let g = SOLID_K * 1e-5 / 0.1;
    for (i, power) in [
        g * (bulk - hot) + kelvin[0] - convection,
        g * (bulk - cold) + kelvin[1],
        g * (hot + cold - 2.0 * bulk) + kelvin[2] + current * current * resistance,
        convection,
    ]
    .into_iter()
    .enumerate()
    {
        let cp = if i == 3 { fluid_cp } else { CP };
        near(MASS * cp * (after[i] - before[i]), DT * power, 2e-8);
    }
    near(
        kelvin[0] / hot + kelvin[1] / cold + kelvin[2] / bulk,
        0.0,
        1e-13,
    );
    near(kelvin.iter().sum::<f64>() + emf * current, 0.0, 1e-12);
    (current, next_charge, kelvin, convection)
}

#[test]
fn g1_g3_coolant_material_changes_kelvin_conversion_and_continued_charge() {
    let run = |viscosity_scale, fluid_cp| {
        let (solid_card, solid_curve) = kelvin_material();
        let (fluid_card, fluid_curve) = fluid_material(viscosity_scale, fluid_cp);
        let solid_context = body_point(300.0);
        let fluid_context = fluid_point(290.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let thermal = initial(&solid_curve, &fluid_curve, fluid_cp);
        let saved_cells = thermal.cells().clone();
        let saved_links = thermal.links().clone();
        let fluid = liquid_law(&fluid_card);
        let cells = [HOT, COLD, BODY, FLUID].map(|cell| {
            let liquid = cell == FLUID;
            MaterialCell {
                cell,
                material_id: MaterialId(if liquid { 1 } else { 0 }),
                curve: if liquid {
                    fluid_curve.curve()
                } else {
                    solid_curve.curve()
                },
                point: ThermalMaterialPoint {
                    material: if liquid { &fluid_card } else { &solid_card },
                    context: if liquid {
                        &fluid_context
                    } else {
                        &solid_context
                    },
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                conductivity: if liquid {
                    None
                } else {
                    Some((&key, MaterialPropertySelection::SingleClaimOnly))
                },
            }
        });
        let resistors = [CellResistor {
            destination: ResistorCell {
                branch: 1,
                cell: BODY,
            },
            selection: MaterialPropertySelection::SingleClaimOnly,
            length_m: 0.1,
            area_m2: 1e-5,
        }];
        let conductive = [0, 1].map(|link| ConductionLink {
            link,
            distance_a_m: 0.05,
            distance_b_m: 0.05,
        });
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &conductive,
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
        let admitted = evaluate_source_link_convection(
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
                admitted.conductances_w_k(),
                || false,
            )
            .unwrap();
        let circuit = MaterialCircuit::new(
            CircuitGraph {
                node_count: 3,
                branches: vec![
                    (1, 2, Branch::VoltageSource { port: 0 }),
                    (2, 0, Branch::Resistor { ohms: 0.0 }),
                    (
                        1,
                        0,
                        Branch::Capacitor {
                            farads: CAPACITANCE,
                        },
                    ),
                ],
                transformers: vec![],
            },
            material
                .resolve_with_link_conductances(&thermal, admitted.conductances_w_k(), || false)
                .unwrap()
                .resistors,
            vec![],
            vec![0.0],
            &[0.0],
        )
        .unwrap();
        let original_circuit = circuit.clone();
        let legs = [ThermoelectricLeg {
            source_branch: 0,
            resistor_branch: 1,
            terminals: [HOT, COLD],
            thermal_links: [0, 1],
            selection: MaterialPropertySelection::SingleClaimOnly,
        }];
        let gate = fs_exec::CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let mut result = [0.0; 3];
        pool.scope(|arena| {
            let cx = fs_exec::Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 181,
                    kernel_id: 29,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            assert!(heat.prepare(&cx, &thermal).is_err());
            let prepared = heat
                .prepare_with_link_conductances(&cx, &thermal, admitted.conductances_w_k())
                .unwrap();
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
            let assembly = network.with_links(&links, 1, |_, _, t| fluid.at_temperature(t));
            let mut solve = controls();
            solve.maximum_cells = 4;
            solve.newton.linear_restart = 4;
            let mut accepted = thermal.clone();
            let mut electrical_state = circuit.clone();
            let mut before = INITIAL;
            for step in 1..=2 {
                let old_cells = accepted.cells().clone();
                let old_circuit = electrical_state.clone();
                let next = step_thermoelectric_candidate(
                    &assembly,
                    &cx,
                    &accepted,
                    step as f64 * DT,
                    &electrical_state,
                    &[0.0],
                    &legs,
                    solve,
                    &[0.0; 4],
                )
                .unwrap();
                let temperatures = [HOT, COLD, BODY, FLUID]
                    .map(|cell| next.thermal().temperature_k(cell).unwrap());
                let (current, charge, kelvin, convected) = reference(
                    before,
                    temperatures,
                    electrical_state.charges_c()[0],
                    viscosity_scale,
                    fluid_cp,
                );
                let report = next.report();
                let coupled = &report.coupled;
                let leg = &coupled.source.legs[0];
                near(leg.transfer.current_a, current, 1e-9);
                near(coupled.electrical.resistors[0].current_a, current, 1e-9);
                near(next.circuit().charges_c()[0], charge, 1e-10);
                near(leg.transfer.terminal_heat_w[0], kelvin[0], 1e-9);
                near(leg.transfer.terminal_heat_w[1], kelvin[1], 1e-9);
                near(leg.transfer.bulk_heat_w, kelvin[2], 1e-9);
                assert!(leg.transfer.bulk_heat_w.abs() > 1e-4);
                assert_eq!(
                    leg.material.bulk_material().card_identity(),
                    solid_card.content_hash()
                );
                let sample = &coupled.source.network;
                let coolant = &sample.finite_convection.as_ref().unwrap().samples()[0];
                assert_eq!((coolant.wall_cell, coolant.fluid_cell), (HOT, FLUID));
                near(coolant.heat_into_fluid_w, convected, 1e-7);
                assert_eq!(
                    coolant.material.material().card_identity(),
                    fluid_card.content_hash()
                );
                assert_eq!(coolant.material.material().query().axes()["p"], 100000.0);
                assert_eq!(
                    coolant.material.material().query().axes()["solute_mass_fraction"],
                    0.2
                );
                let conductor = &coupled.electrical.resistors[0].material.conductor;
                assert_eq!(
                    conductor.material().card_identity(),
                    solid_card.content_hash()
                );
                near(
                    conductor.material().query().axes()["T"],
                    temperatures[2],
                    1e-12,
                );
                near(
                    coupled.deposited_joule_j,
                    DT * current * current * conductor.resistance_ohm(),
                    1e-10,
                );
                near(
                    coupled.source.reversible_heat_j,
                    DT * kelvin.iter().sum::<f64>(),
                    1e-10,
                );
                near(
                    coupled.source.reversible_heat_j + coupled.source.reversible_electrical_work_j,
                    0.0,
                    1e-10,
                );
                near(report.external_heat_j, 0.0, 1e-10);
                near(report.external_electrical_work_j, 0.0, 1e-10);
                near(
                    next.thermal().energy_j() - accepted.energy_j()
                        + coupled.electrical.stored_after_j
                        - coupled.electrical.stored_before_j,
                    0.0,
                    2e-8,
                );
                near(report.total_energy_residual_j, 0.0, 2e-8);
                assert_eq!(accepted.cells(), &old_cells);
                assert_eq!(electrical_state, old_circuit);
                assert_eq!(next.thermal().links(), &saved_links);
                before = temperatures;
                result = [temperatures[0], current, charge];
                let (next_thermal, next_circuit, _) = next.into_parts();
                accepted = next_thermal;
                electrical_state = next_circuit;
            }
            let queries = std::cell::Cell::new(0);
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
            let failing = refused.with_links(&links, 1, |_, _, t| {
                queries.set(queries.get() + 1);
                if queries.get() > 1 {
                    Err("finite coolant source refused")
                } else {
                    Ok(fluid.at_temperature(t).unwrap())
                }
            });
            assert!(matches!(
                step_thermoelectric_candidate(
                    &failing,
                    &cx,
                    &thermal,
                    DT,
                    &circuit,
                    &[0.0],
                    &legs,
                    solve,
                    &[0.0; 4],
                ),
                Err(CircuitEnthalpyError::Source(
                    ThermoelectricNetworkError::Heat(ConvectionHeatError::Convection(
                        SourceCellConvectionError::Source {
                            cell: FLUID,
                            source: "finite coolant source refused",
                            ..
                        }
                    ))
                ))
            ));
            assert!(queries.get() > 1);
            gate.request();
            assert!(matches!(
                step_thermoelectric_candidate(
                    &assembly,
                    &cx,
                    &thermal,
                    DT,
                    &circuit,
                    &[0.0],
                    &legs,
                    solve,
                    &[0.0; 4],
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            assert_eq!(thermal.cells(), &saved_cells);
            assert_eq!(thermal.links(), &saved_links);
            assert_eq!(circuit, original_circuit);
        });
        result
    };
    let baseline = run(1.0, 4000.0);
    let viscous = run(4.0, 4000.0);
    let higher_capacity = run(1.0, 8000.0);
    assert!(viscous[0] > baseline[0] && viscous[1] > baseline[1] && viscous[2] < baseline[2]);
    assert!(
        higher_capacity[0] < baseline[0]
            && higher_capacity[1] < baseline[1]
            && higher_capacity[2] > baseline[2]
    );
}
