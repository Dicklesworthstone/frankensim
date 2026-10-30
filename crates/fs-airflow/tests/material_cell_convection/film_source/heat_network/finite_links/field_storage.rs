//! Source-only coolant and field-law changes reach original charge/flux history.
use super::*;
use fs_couple::circuit_enthalpy::field_storage::{
    FieldThermalBinding,
    heat::{HeatFieldError, field_storage_candidate},
};
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_material::field_storage::{
    ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS,
    MAGNETIC_PERMEABILITY_PROPERTY, ResolvedCapacitor, ResolvedInductor,
    thermal::ThermalFieldStorage,
};

// Rebuild calorics after adding the field claim: every owner must retain the
// same complete card, not merely equal numeric Cp/R values from an older card.
fn field_body(magnetic: bool, alpha: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = body_material();
    let mut claims = original.claims().clone();
    let (name, dims) = if magnetic {
        (MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS)
    } else {
        (ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS)
    };
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (280.0, CAPACITANCE * (1.0 - 30.0 * alpha)),
                    (320.0, CAPACITANCE * (1.0 + 10.0 * alpha)),
                ],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", temperature_quantity(), 280.0, 320.0)
                .with("lot", 7.0, 7.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "authored affine field/finite-coolant fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
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

fn coordinate(circuit: &MaterialCircuit, magnetic: bool) -> f64 {
    if magnetic {
        circuit.fluxes_wb()[0]
    } else {
        circuit.charges_c()[0]
    }
}

fn field_u(alpha: f64, temperature: f64, z: f64) -> f64 {
    let k = CAPACITANCE * (1.0 + alpha * (temperature - 310.0));
    let slope = CAPACITANCE * alpha;
    z * z / (2.0 * k) * (1.0 + temperature * slope / k)
}

// Independently eliminate fluid temperature by total caloric + field U.
// Original midpoint RC/RL decay and endpoint finite convection leave one scalar
// equation. No production field, convection or coupled solver is called here.
fn reference_step(magnetic: bool, alpha: f64, viscosity: f64, cp: f64, old: [f64; 3]) -> [f64; 5] {
    let evaluate = |wall: f64| {
        let k = CAPACITANCE * (1.0 + alpha * (wall - 310.0));
        let r = 1.0 + 0.02 * (wall - 300.0);
        let ratio = if magnetic {
            DT * r / (2.0 * k)
        } else {
            DT / (2.0 * k * r)
        };
        let z = old[2] * (1.0 - ratio) / (1.0 + ratio);
        let du = field_u(alpha, wall, z) - field_u(alpha, old[0], old[2]);
        let fluid = old[1] - (MASS * CP * (wall - old[0]) + du) / (MASS * cp);
        let film = 0.5 * (wall + fluid);
        let conductivity = 0.5 + 0.002 * (film - 300.0);
        let mu = viscosity * 0.001 * (1.0 - 0.01 * (film - 300.0));
        let h = 0.664 * (1000.0 * 0.1 / mu).sqrt() * (mu * cp / conductivity).cbrt() * conductivity
            / 0.1;
        let power = SURFACE * h * (wall - fluid);
        let residual = MASS * CP * (wall - old[0]) + du + DT * power;
        let current = (old[2] + z) / (2.0 * k) / if magnetic { 1.0 } else { r };
        (residual, [wall, fluid, z, current, power])
    };
    let (mut lower, mut upper) = ((old[1] - 2.0).max(280.01), (old[0] + 5.0).min(319.99));
    assert!(evaluate(lower).0 < 0.0 && evaluate(upper).0 > 0.0);
    for _ in 0..80 {
        let middle = 0.5 * (lower + upper);
        if evaluate(middle).0 > 0.0 {
            upper = middle;
        } else {
            lower = middle;
        }
    }
    let (residual, result) = evaluate(0.5 * (lower + upper));
    assert!(residual.abs() < 1e-10);
    result
}

#[test]
fn finite_coolant_and_field_data_change_native_charge_flux_and_total_internal_energy() {
    let run = |magnetic, alpha, viscosity, fluid_cp| {
        let (body_card, body_curve) = field_body(magnetic, alpha);
        let (fluid_card, fluid_curve) = fluid_material(viscosity, fluid_cp);
        let body_context = body_point(310.0);
        let fluid_context = fluid_point(290.0);
        let mut thermal = initial_pair(&body_curve, &fluid_curve, fluid_cp, false, 0.0);
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
        let material = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &[],
        };
        let electrical = ElectricalContactNetwork {
            bulk: &material,
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
        let selection = MaterialPropertySelection::SingleClaimOnly;
        // Both original reductions occupy precisely the body's 1e-6 m^3.
        let (branch, element, law) = if magnetic {
            let value =
                ResolvedInductor::resolve(&body_card, &body_context, selection, 1e-3, 1e-3, 1)
                    .unwrap();
            let law =
                ThermalFieldStorage::inductor(&body_card, &value, "T", [280.0, 320.0]).unwrap();
            (
                Branch::Inductor { henries: f64::NAN },
                MaterialStorageElement::Inductor(value),
                law,
            )
        } else {
            let value =
                ResolvedCapacitor::resolve(&body_card, &body_context, selection, 1e-3, 1e-3)
                    .unwrap();
            let law =
                ThermalFieldStorage::capacitor(&body_card, &value, "T", [280.0, 320.0]).unwrap();
            (
                Branch::Capacitor { farads: f64::NAN },
                MaterialStorageElement::Capacitor(value),
                law,
            )
        };
        let fields = [FieldThermalBinding {
            branch: 0,
            cell: BODY,
            law,
            maximum_abs_coordinate: INITIAL_CHARGE,
        }];
        let gate = fs_exec::CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = fs_exec::Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 174,
                    kernel_id: 29,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            let prepared = heat
                .prepare_with_link_conductances(&cx, &thermal, initial.conductances_w_k())
                .unwrap();
            let mut circuit = MaterialCircuit::new_with_storage(
                CircuitGraph {
                    node_count: 2,
                    branches: vec![(1, 0, branch), (1, 0, Branch::Resistor { ohms: 0.0 })],
                    transformers: vec![],
                },
                prepared
                    .resolve_with_link_conductances(&cx, &thermal, initial.conductances_w_k())
                    .unwrap()
                    .bulk
                    .resistors,
                vec![MaterialStorage { branch: 0, element }],
                if magnetic {
                    vec![INITIAL_CHARGE]
                } else {
                    vec![]
                },
                if magnetic {
                    vec![]
                } else {
                    vec![INITIAL_CHARGE]
                },
                &[],
            )
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
            solve.maximum_cells = 2;
            solve.newton.linear_restart = 2;
            let total = thermal.energy_j() + field_u(alpha, 310.0, INITIAL_CHARGE);
            let mut old = [310.0, 290.0, INITIAL_CHARGE];
            let mut endpoints = [[0.0; 3]; 2];
            for step in 1..=2 {
                let expected = reference_step(magnetic, alpha, viscosity, fluid_cp, old);
                let next = field_storage_candidate(
                    &assembly,
                    &cx,
                    &thermal,
                    step as f64 * DT,
                    &circuit,
                    &[],
                    &fields,
                    solve,
                    &[0.0; 2],
                )
                .unwrap();
                let actual = [
                    next.thermal().temperature_k(BODY).unwrap(),
                    next.thermal().temperature_k(FLUID).unwrap(),
                    coordinate(next.circuit(), magnetic),
                ];
                for (value, reference) in actual.into_iter().zip(expected) {
                    near(value, reference, 1e-8);
                }
                let report = next.report();
                near(
                    report.coupled.electrical.resistors[0].current_a.abs(),
                    expected[3],
                    1e-8,
                );
                let source = &report.coupled.source;
                let path = &source.network.finite_convection.as_ref().unwrap().samples()[0];
                near(path.heat_into_fluid_w, expected[4], 1e-7);
                near(
                    path.property_temperature_k,
                    0.5 * (actual[0] + actual[1]),
                    1e-12,
                );
                assert_eq!(
                    path.material.material().card_identity(),
                    fluid_card.content_hash()
                );
                let field = &source.fields[0];
                assert_eq!(
                    field.material.material().card_identity(),
                    body_card.content_hash()
                );
                near(
                    field.material.material().query().axes()["T"],
                    actual[0],
                    1e-12,
                );
                near(
                    field.material.coefficient(),
                    CAPACITANCE * (1.0 + alpha * (actual[0] - 310.0)),
                    1e-12,
                );
                near(
                    source.field_internal_energy_j[1],
                    field_u(alpha, expected[0], expected[2]),
                    1e-8,
                );
                near(
                    next.thermal().energy_j() + source.field_internal_energy_j[1],
                    total,
                    1e-8,
                );
                near(report.external_heat_j, 0.0, 1e-12);
                near(report.total_energy_residual_j, 0.0, 1e-8);
                near(
                    MASS * fluid_cp * (actual[1] - old[1]),
                    DT * path.heat_into_fluid_w,
                    1e-8,
                );
                assert!(source.network.heat.material.links.is_empty());
                if alpha != 0.0 {
                    assert!(source.reversible_heat_j.abs() > 1e-3);
                }
                let (t, c, _) = next.into_parts();
                thermal = t;
                circuit = c;
                endpoints[step - 1] = actual;
                old = [expected[0], expected[1], expected[2]];
            }

            if !magnetic && alpha == 0.004 && viscosity == 1.0 && fluid_cp == 4000.0 {
                let saved = thermal.cells().clone();
                let saved_links = thermal.links().clone();
                let saved_circuit = circuit.clone();
                // Omitting the finite link law cannot turn it into insulation.
                assert!(
                    field_storage_candidate(
                        &network,
                        &cx,
                        &thermal,
                        3.0 * DT,
                        &circuit,
                        &[],
                        &fields,
                        solve,
                        &[0.0; 2],
                    )
                    .is_err()
                );
                let mut exhausted = solve;
                exhausted.maximum_circuit_trials = 1;
                assert!(matches!(
                    field_storage_candidate(
                        &assembly,
                        &cx,
                        &thermal,
                        3.0 * DT,
                        &circuit,
                        &[],
                        &fields,
                        exhausted,
                        &[0.0; 2],
                    ),
                    Err(CircuitEnthalpyError::CircuitBudget)
                ));
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
                let refused = refused.with_links(&links, 1, |_, _, _| {
                    Err::<fs_material::liquid::ResolvedLiquidState, _>(
                        "finite fluid source refused",
                    )
                });
                assert!(matches!(
                    field_storage_candidate(
                        &refused,
                        &cx,
                        &thermal,
                        3.0 * DT,
                        &circuit,
                        &[],
                        &fields,
                        solve,
                        &[0.0; 2],
                    ),
                    Err(CircuitEnthalpyError::Source(HeatFieldError::Heat(
                        ConvectionHeatError::Convection(SourceCellConvectionError::Source {
                            cell: FLUID,
                            source: "finite fluid source refused",
                            ..
                        }),
                    )))
                ));
                let retry = field_storage_candidate(
                    &assembly,
                    &cx,
                    &thermal,
                    3.0 * DT,
                    &circuit,
                    &[],
                    &fields,
                    solve,
                    &[0.0; 2],
                )
                .unwrap();
                let replay = field_storage_candidate(
                    &assembly,
                    &cx,
                    &thermal,
                    3.0 * DT,
                    &circuit,
                    &[],
                    &fields,
                    solve,
                    &[0.0; 2],
                )
                .unwrap();
                assert_eq!(retry.circuit(), replay.circuit());
                assert_eq!(retry.thermal().cells(), replay.thermal().cells());
                gate.request();
                assert!(matches!(
                    field_storage_candidate(
                        &assembly,
                        &cx,
                        &thermal,
                        3.0 * DT,
                        &circuit,
                        &[],
                        &fields,
                        solve,
                        &[0.0; 2],
                    ),
                    Err(CircuitEnthalpyError::Cancelled)
                ));
                assert_eq!(thermal.cells(), &saved);
                assert_eq!(thermal.links(), &saved_links);
                assert_eq!(circuit, saved_circuit);
                assert_eq!(thermal.time_s(), 2.0 * DT);
            }
            endpoints
        })
    };
    for magnetic in [false, true] {
        let baseline = run(magnetic, 0.004, 1.0, 4000.0);
        for changed in [
            run(magnetic, 0.004, 4.0, 4000.0),
            run(magnetic, 0.004, 1.0, 8000.0),
            run(magnetic, 0.0, 1.0, 4000.0),
        ] {
            for (old, new) in baseline[1].iter().zip(changed[1]) {
                assert!((old - new).abs() > 1e-5);
            }
        }
    }
}
