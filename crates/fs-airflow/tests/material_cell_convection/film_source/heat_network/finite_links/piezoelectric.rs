//! Source-only coolant substitution changes the original reciprocal piezo rod
//! through load temperature and current, with one finite thermal inventory.
use super::*;
use fs_couple::circuit_enthalpy::piezoelectric::heat as piezo_heat;
use fs_couple::material_circuit::piezoelectric::PiezoelectricCircuit;
use fs_couple::piezoelectric::{PiezoelectricRod, SmallSignalLimits};
use fs_material::field_storage::ELECTRIC_PERMITTIVITY_DIMS;
use fs_material::piezoelectric::{
    COMPLIANCE_33_E_PROPERTY, COMPLIANCE_DIMS, LongitudinalPiezoelectric,
    PERMITTIVITY_33_T_PROPERTY, STRAIN_CHARGE_33_PROPERTY, STRAIN_CHARGE_DIMS,
};

const EXTENSION: f64 = 0.02;
const VELOCITY: f64 = 0.01;

// Authored reciprocal material, held independently of the finite load/coolant.
// L=2, A=.5 gives m=1, k=2.5, C=.4 and theta=.5 in SI units.
fn rod() -> PiezoelectricRod {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (COMPLIANCE_33_E_PROPERTY, COMPLIANCE_DIMS, 0.1),
        (STRAIN_CHARGE_33_PROPERTY, STRAIN_CHARGE_DIMS, 0.2),
        (PERMITTIVITY_33_T_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS, 2.0),
        (DENSITY_PROPERTY, Density::DIMS, 3.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic fixed-ceramic coolant coupling fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored piezoelectric solid".into(),
            phase: "solid".into(),
            process: "synthetic reciprocal longitudinal reduction".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let material = LongitudinalPiezoelectric::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    PiezoelectricRod::new(
        material,
        2.0,
        0.5,
        SmallSignalLimits {
            max_abs_strain: 0.05,
            max_abs_electric_field_v_m: 10.0,
        },
    )
    .unwrap()
}

// Independent midpoint Kirchhoff/mechanical equations. Thermal conservation
// eliminates the coolant temperature; the remaining wall balance is scalar.
// No production convection, circuit or thermal solver supplies this reference.
fn reference(viscosity_scale: f64, cp: f64, start: [f64; 5]) -> [f64; 8] {
    let [wall0, fluid0, x0, v0, q0] = start;
    let evaluate = |wall: f64| {
        let resistance = 1.0 + 0.02 * (wall - 300.0);
        let (mass, stiffness, capacitance, theta) = (1.0, 2.5, 0.4, 0.5);
        let denominator = 1.0 + DT * DT * stiffness / (4.0 * mass);
        let a = (x0 + 0.5 * DT * v0) / denominator;
        let b = DT * DT * theta / (4.0 * mass * denominator);
        let voltage = (q0 - theta * a) / (capacitance + DT / (2.0 * resistance) + theta * b);
        let midpoint_x = a + b * voltage;
        let x = 2.0 * midpoint_x - x0;
        let v = 4.0 * (midpoint_x - x0) / DT - v0;
        let current = voltage / resistance;
        let charge = q0 - DT * current;
        let joule = DT * voltage * current;
        let fluid = fluid0 + (joule - MASS * CP * (wall - wall0)) / (MASS * cp);
        let film = 0.5 * (wall + fluid);
        let k = 0.5 + 0.002 * (film - 300.0);
        let mu = viscosity_scale * 0.001 * (1.0 - 0.01 * (film - 300.0));
        let h = 0.664 * (1000.0 * 0.1 / mu).sqrt() * (mu * cp / k).cbrt() * k / 0.1;
        let power = SURFACE * h * (wall - fluid);
        let residual = MASS * CP * (wall - wall0) - joule + DT * power;
        (residual, [wall, fluid, x, v, charge, current, joule, power])
    };
    let (mut low, mut high) = (280.0, 320.0);
    assert!(evaluate(low).0 < 0.0 && evaluate(high).0 > 0.0);
    for _ in 0..100 {
        let midpoint = 0.5 * (low + high);
        if evaluate(midpoint).0 > 0.0 {
            high = midpoint;
        } else {
            low = midpoint;
        }
    }
    evaluate(0.5 * (low + high)).1
}

#[test]
fn coolant_material_changes_current_and_rod_motion_on_the_original_joint_owner() {
    let run = |viscosity_scale, fluid_cp| {
        let (body_card, body_curve) = body_material();
        let (fluid_card, fluid_curve) = fluid_material(viscosity_scale, fluid_cp);
        let body_context = body_point(310.0);
        let fluid_context = fluid_point(290.0);
        let thermal = initial_pair(&body_curve, &fluid_curve, fluid_cp, false, 0.0);
        let saved_cells = thermal.cells().clone();
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
        let model = PiezoelectricCircuit::new(
            CircuitGraph {
                node_count: 2,
                branches: vec![
                    (1, 0, Branch::Capacitor { farads: f64::NAN }),
                    (1, 0, Branch::Resistor { ohms: f64::NAN }),
                ],
                transformers: vec![],
            },
            material
                .resolve_with_link_conductances(&thermal, initial.conductances_w_k(), || false)
                .unwrap()
                .resistors,
            vec![],
            vec![INITIAL_CHARGE],
            0,
            rod(),
            EXTENSION,
            VELOCITY,
            &[],
            || false,
        )
        .unwrap();
        let saved_model = model.clone();
        let gate = fs_exec::CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let mut result = [0.0; 3];
        pool.scope(|arena| {
            let cx = fs_exec::Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 177,
                    kernel_id: 31,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            let prepared = heat
                .prepare_with_link_conductances(&cx, &thermal, initial.conductances_w_k())
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
            let mut accepted_thermal = thermal.clone();
            let mut accepted_model = model.clone();
            let mut state = [310.0, 290.0, EXTENSION, VELOCITY, INITIAL_CHARGE];
            for step in 1..=2 {
                let before_cells = accepted_thermal.cells().clone();
                let before_model = accepted_model.clone();
                let expected = reference(viscosity_scale, fluid_cp, state);
                let next = piezo_heat::step_candidate(
                    &assembly,
                    &cx,
                    &accepted_thermal,
                    step as f64 * DT,
                    &accepted_model,
                    &[],
                    0.0,
                    solve,
                    &[0.0, 0.0],
                )
                .unwrap();
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
                near(next.model().extension_m(), expected[2], 1e-10);
                near(next.model().velocity_m_s(), expected[3], 1e-10);
                near(next.model().charge_c(), expected[4], 1e-10);
                near(report.electrical.resistors[0].current_a, expected[5], 1e-9);
                let conductor = &report.electrical.resistors[0].material.conductor;
                assert_eq!(
                    conductor.material().card_identity(),
                    body_card.content_hash()
                );
                assert_eq!(conductor.material().query().axes()["lot"], 7.0);
                near(
                    conductor.material().query().axes()["T"],
                    next.thermal().temperature_k(BODY).unwrap(),
                    1e-12,
                );
                near(report.deposited_joule_j, expected[6], 1e-9);
                near(report.deposited_joule_j, report.electrical.heat_j, 1e-10);
                near(report.electrical.parameter_work_j, 0.0, 1e-14);
                near(report.non_electrical_heat_j, 0.0, 1e-10);
                near(
                    next.thermal().energy_j() - accepted_thermal.energy_j()
                        + report.electrical.stored_after_j
                        - report.electrical.stored_before_j,
                    0.0,
                    1e-8,
                );
                near(
                    MASS * fluid_cp * (expected[1] - state[1]),
                    DT * expected[7],
                    1e-8,
                );
                near(
                    MASS * CP * (expected[0] - state[0]),
                    report.deposited_joule_j - DT * expected[7],
                    1e-8,
                );
                let source = &report.source;
                let finite = source.finite_convection.as_ref().unwrap();
                let sample = &finite.samples()[0];
                assert_eq!((sample.wall_cell, sample.fluid_cell), (BODY, FLUID));
                near(sample.heat_into_fluid_w, expected[7], 1e-7);
                near(
                    sample.property_temperature_k,
                    0.5 * (expected[0] + expected[1]),
                    1e-8,
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
                assert!(source.heat.material.links.is_empty());
                for cell in &source.heat.material.cells {
                    near(
                        cell.phase.temperature_k(),
                        next.thermal().temperature_k(cell.cell).unwrap(),
                        1e-12,
                    );
                }
                near(
                    source.convection.convection().total_heat_into_cells_w(),
                    0.0,
                    1e-14,
                );
                assert_eq!(next.model().rod(), model.rod());
                assert_eq!(next.thermal().links(), &saved_links);
                assert_eq!(accepted_thermal.cells(), &before_cells);
                assert_eq!(accepted_model, before_model);
                result = [
                    next.thermal().temperature_k(BODY).unwrap(),
                    report.electrical.resistors[0].current_a,
                    next.model().extension_m(),
                ];
                state = [
                    expected[0],
                    expected[1],
                    expected[2],
                    expected[3],
                    expected[4],
                ];
                let (new_thermal, new_model, _) = next.into_parts();
                accepted_thermal = new_thermal;
                accepted_model = new_model;
            }

            // Native circuit budget refusal remains typed and publishes neither
            // caloric nor mechanical/electrical history.
            let mut insufficient = solve;
            insufficient.maximum_circuit_trials = 1;
            assert!(matches!(
                piezo_heat::step_candidate(
                    &assembly,
                    &cx,
                    &thermal,
                    DT,
                    &model,
                    &[],
                    0.0,
                    insufficient,
                    &[0.0, 0.0]
                ),
                Err(CircuitEnthalpyError::CircuitBudget)
            ));

            // Refuse after accepted source admission, during a thermal trial.
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
                    Err("coolant source refused during trial")
                } else {
                    Ok(fluid.at_temperature(t).unwrap())
                }
            });
            assert!(matches!(
                piezo_heat::step_candidate(
                    &failing,
                    &cx,
                    &thermal,
                    DT,
                    &model,
                    &[],
                    0.0,
                    solve,
                    &[0.0, 0.0]
                ),
                Err(CircuitEnthalpyError::Source(
                    ConvectionHeatError::Convection(SourceCellConvectionError::Source {
                        cell: FLUID,
                        source: "coolant source refused during trial",
                        ..
                    })
                ))
            ));
            assert!(queries.get() > 1);
            gate.request();
            assert!(matches!(
                piezo_heat::step_candidate(
                    &assembly,
                    &cx,
                    &thermal,
                    DT,
                    &model,
                    &[],
                    0.0,
                    solve,
                    &[0.0, 0.0]
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            assert_eq!(thermal.cells(), &saved_cells);
            assert_eq!(thermal.links(), &saved_links);
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(model, saved_model);
        });
        result
    };
    let baseline = run(1.0, 4000.0);
    let viscous = run(4.0, 4000.0);
    let higher_capacity = run(1.0, 8000.0);
    assert!(viscous[0] > baseline[0] && viscous[1] < baseline[1] && viscous[2] > baseline[2]);
    assert!(
        higher_capacity[0] < baseline[0]
            && higher_capacity[1] > baseline[1]
            && higher_capacity[2] < baseline[2]
    );
}
