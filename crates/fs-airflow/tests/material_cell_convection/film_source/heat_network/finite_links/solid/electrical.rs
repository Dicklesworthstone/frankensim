//! Finite coolant, sourced R(T), native J2 and charge/velocity share one endpoint.
use super::*;
use fs_couple::circuit_enthalpy::contacts::heat::HeatAssembly;
use fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics;
use fs_couple::solid_enthalpy::plastic::electrical::{
    ElectricalSolidError, ElectricalSolidSolve, ElectricalSolidSourceError,
};

fn powered_controls() -> ElectricalSolidSolve {
    ElectricalSolidSolve {
        coupled: CircuitEnthalpySolve {
            maximum_cells: 2,
            maximum_circuit_trials: 1000,
            maximum_newton_iterations: 32,
            newton: NewtonKrylovConfig {
                linear_restart: 2,
                ..controls().newton
            },
            ..controls()
        },
        maximum_mechanical_trials: 1000,
        yield_tolerance_pa: 1e-3,
    }
}

// The limited case changes only original viscosity support. Caloric and
// transport bindings still refer to this exact same immutable material card.
fn coolant(cp: f64, limited: bool) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, original_curve) = fluid_material(1.0, cp);
    if !limited {
        return (original, original_curve);
    }
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == DYNAMIC_VISCOSITY_PROPERTY {
            claim.validity = ValidityDomain::unconstrained()
                .with_quantity("T", temperature_quantity(), 300.0, 320.0)
                .with("p", 100000.0, 100000.0)
                .with("solute_mass_fraction", 0.2, 0.2);
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

// Independent radial-return, midpoint capacitor and backward-Euler finite-link
// equations. No production thermal, electrical, convection or solid solver.
fn reference(cp: f64) -> [f64; 5] {
    let evaluate = |wall: f64| {
        let (displacement, alpha) = apex_reference(wall);
        let plastic = VOLUME * (3e6 - 1e5 * (wall - 300.0)) * alpha;
        let resistance = 1.0 + 0.02 * (wall - 300.0);
        let ratio = DT / (2.0 * CAPACITANCE * resistance);
        let charge = INITIAL_CHARGE * (1.0 - ratio) / (1.0 + ratio);
        let joule = (INITIAL_CHARGE.powi(2) - charge.powi(2)) / (2.0 * CAPACITANCE);
        let fluid = 290.0 + (joule + plastic - MASS * CP * (wall - 310.0)) / (MASS * cp);
        let film = 0.5 * (wall + fluid);
        let conductivity = 0.5 + 0.002 * (film - 300.0);
        let viscosity = 0.001 * (1.0 - 0.01 * (film - 300.0));
        let h = 0.664
            * (1000.0 * 0.1 / viscosity).sqrt()
            * (viscosity * cp / conductivity).cbrt()
            * conductivity
            / 0.1;
        (
            MASS * CP * (wall - 310.0) - joule - plastic + DT * SURFACE * h * (wall - fluid),
            [wall, fluid, charge, displacement, alpha],
        )
    };
    let (mut low, mut high) = (290.0, 312.0);
    assert!(evaluate(low).0 < 0.0 && evaluate(high).0 > 0.0);
    for _ in 0..80 {
        let t = 0.5 * (low + high);
        if evaluate(t).0 > 0.0 {
            high = t;
        } else {
            low = t;
        }
    }
    let (residual, expected) = evaluate(0.5 * (low + high));
    assert!(residual.abs() < 1e-10);
    expected
}

#[test]
fn finite_coolant_closes_powered_plasticity_and_motion_without_losing_source_history() {
    let mut observations = Vec::new();
    for (fluid_cp, limited) in [(2000.0, false), (8000.0, false), (2000.0, true)] {
        let gate = fs_exec::CancelGate::new();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
            let cx = fs_exec::Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 174,
                    kernel_id: 0x534f4c45,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            let (body_card, body_curve) = solid_material();
            let (fluid_card, fluid_curve) = coolant(fluid_cp, limited);
            let thermal = initial_pair(&body_curve, &fluid_curve, fluid_cp, false, 0.0);
            let solid = ThermoplasticSolid::try_new(
                &thermal,
                &NODES,
                &TETS,
                vec![ThermoplasticCell {
                    cell: BODY,
                    query: SourceQuery {
                        card: &body_card,
                        caloric: &body_curve,
                        context: body_point(310.0),
                        temperature_axis: "T".into(),
                        selection: MaterialPropertySelection::SingleClaimOnly,
                        hardening_basis: BASIS,
                        thermal_expansion: None,
                        strain_limit: 0.04,
                    },
                }],
                Default::default(),
                &cx,
            )
            .unwrap();
            let original_cells = thermal.cells().clone();
            let original_solid = solid.solution().clone();
            let fluid_context = fluid_point(290.0);
            let fluid = liquid_law(&fluid_card);
            let roles = [None];
            let mut cells = solid.thermal_material_cells(&roles).unwrap();
            cells.push(MaterialCell {
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
            let initial = assembly.resolve_heat(&cx, &thermal).unwrap();
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
                initial.bulk.resistors,
                vec![],
                vec![INITIAL_CHARGE],
                &[],
            )
            .unwrap();
            let load = [(9, TRACTION * VOLUME / NODES[3][2])];
            let result = solid.advance_electrical_candidate(
                &thermal,
                &circuit,
                &assembly,
                &[],
                DT,
                &supports(SHEAR),
                &load,
                Default::default(),
                powered_controls(),
                &[0.0; 2],
                &cx,
            );
            if limited {
                assert!(matches!(
                    result,
                    Err(ElectricalSolidError::Coupled(CircuitEnthalpyError::Source(
                        ElectricalSolidSourceError::Heat(ConvectionHeatError::Convection(
                            SourceCellConvectionError::Source { .. }
                        ))
                    )))
                ));
                assert_eq!(thermal.cells(), &original_cells);
                assert_eq!(solid.solution(), &original_solid);
                assert_eq!(circuit.charges_c(), &[INITIAL_CHARGE]);
                return;
            }
            let loaded = result.unwrap();
            let expected = reference(fluid_cp);
            let wall = loaded.thermal().temperature_k(BODY).unwrap();
            let fluid_t = loaded.thermal().temperature_k(FLUID).unwrap();
            near(wall, expected[0], 2e-8);
            near(fluid_t, expected[1], 2e-8);
            near(loaded.circuit().charges_c()[0], expected[2], 1e-9);
            near(
                loaded.solid().solution().displacement_m[3][0],
                expected[3],
                2e-11,
            );
            near(
                loaded.solid().solution().element_states[0].alpha,
                expected[4],
                2e-10,
            );
            let report = loaded.report();
            assert!(report.physical_plastic_heat_j > 0.0 && report.numerical_loss_j > 0.0);
            assert!(report.circuit_thermal.deposited_joule_j > 0.0);
            near(report.external_heat_j, 0.0, 1e-12);
            near(report.combined_residual_j, 0.0, 1e-8);
            near(report.thermal_residual_j, 0.0, 1e-8);
            near(
                loaded.thermal().energy_j() - thermal.energy_j(),
                report.physical_plastic_heat_j + report.circuit_thermal.deposited_joule_j,
                1e-8,
            );
            let source = &report.circuit_thermal.source.material;
            let finite = source.finite_convection.as_ref().unwrap();
            let sample = &finite.samples()[0];
            near(sample.wall_temperature_k, wall, 1e-12);
            near(sample.fluid_temperature_k, fluid_t, 1e-12);
            near(sample.property_temperature_k, 0.5 * (wall + fluid_t), 1e-12);
            assert_eq!(
                sample.material.material().card_identity(),
                fluid_card.content_hash()
            );
            for property in sample.material.material().properties() {
                assert_eq!(
                    property.answer().receipt.selected,
                    fluid
                        .reference_state()
                        .material()
                        .property(property.requirement().name())
                        .unwrap()
                        .answer()
                        .receipt
                        .selected
                );
            }
            let conductor = &report.circuit_thermal.electrical.resistors[0]
                .material
                .conductor;
            near(
                conductor.resistance_ohm(),
                1.0 + 0.02 * (wall - 300.0),
                1e-10,
            );
            assert_eq!(
                conductor.material().card_identity(),
                body_card.content_hash()
            );
            assert_eq!(
                loaded.solid().materials()[0].material().card_identity(),
                body_card.content_hash()
            );
            let endpoint = assembly.resolve_heat(&cx, loaded.thermal()).unwrap();
            let native = circuit
                .step_candidate(
                    endpoint.bulk.resistors,
                    &[],
                    DT,
                    powered_controls().coupled.circuit,
                    || false,
                )
                .unwrap();
            assert_eq!(native.circuit(), loaded.circuit());
            // Omitting the original finite-link owner must remain a refusal.
            assert!(
                solid
                    .advance_electrical_candidate(
                        &thermal,
                        &circuit,
                        &prepared,
                        &[],
                        DT,
                        &supports(SHEAR),
                        &load,
                        Default::default(),
                        powered_controls(),
                        &[0.0; 2],
                        &cx
                    )
                    .is_err()
            );
            let history = loaded.solid().solution().element_states.clone();
            let continued = loaded
                .solid()
                .advance_electrical_candidate(
                    loaded.thermal(),
                    loaded.circuit(),
                    &assembly,
                    &[],
                    2.0 * DT,
                    &supports(SHEAR - 0.0005),
                    &[],
                    Default::default(),
                    powered_controls(),
                    &[0.0; 2],
                    &cx,
                )
                .unwrap();
            assert_eq!(continued.solid().solution().element_states, history);
            assert!(continued.circuit().charges_c()[0] < loaded.circuit().charges_c()[0]);
            near(continued.report().combined_residual_j, 0.0, 1e-8);
            near(
                continued.thermal().energy_j() - loaded.thermal().energy_j(),
                continued.report().circuit_thermal.deposited_joule_j,
                1e-8,
            );
            assert!(
                continued
                    .report()
                    .circuit_thermal
                    .source
                    .material
                    .finite_convection
                    .is_some()
            );

            let mut velocity = vec![[0.0; 3]; NODES.len()];
            velocity[3][0] = 2.0;
            let moving =
                ThermoplasticDynamics::try_new(solid.clone(), velocity.clone(), &cx).unwrap();
            let dynamic = moving
                .advance_electrical_candidate(
                    &thermal,
                    &circuit,
                    &assembly,
                    &[],
                    DT,
                    &supports(0.0),
                    &[],
                    Default::default(),
                    powered_controls(),
                    &[0.0; 2],
                    &cx,
                )
                .unwrap();
            let mass = MASS / 10.0;
            let stiffness = 5e8 * VOLUME / NODES[3][2].powi(2);
            near(
                dynamic.dynamics().solid().solution().displacement_m[3][0],
                mass * DT * 2.0 / (mass + stiffness * DT * DT),
                1e-12,
            );
            let dynamic_report = dynamic.report();
            assert!(dynamic_report.inertia.numerical_loss_j > 1e-5);
            near(dynamic_report.material.physical_plastic_heat_j, 0.0, 1e-10);
            near(
                dynamic.thermal().energy_j() - thermal.energy_j(),
                dynamic_report.material.circuit_thermal.deposited_joule_j,
                1e-8,
            );
            near(dynamic_report.material.external_heat_j, 0.0, 1e-12);
            near(dynamic_report.material.combined_residual_j, 0.0, 1e-8);
            assert!(
                dynamic_report
                    .material
                    .circuit_thermal
                    .source
                    .material
                    .finite_convection
                    .is_some()
            );
            gate.request();
            assert!(matches!(
                moving.advance_electrical_candidate(
                    &thermal,
                    &circuit,
                    &assembly,
                    &[],
                    DT,
                    &supports(0.0),
                    &[],
                    Default::default(),
                    powered_controls(),
                    &[0.0; 2],
                    &cx
                ),
                Err(ElectricalSolidError::Cancelled)
            ));
            assert_eq!(moving.velocity_mps(), velocity);
            assert_eq!(moving.solid().solution(), &original_solid);
            assert_eq!(thermal.cells(), &original_cells);
            assert_eq!(solid.solution(), &original_solid);
            assert_eq!(circuit.charges_c(), &[INITIAL_CHARGE]);
            observations.push([wall, fluid_t, loaded.circuit().charges_c()[0], expected[3]]);
        });
    }
    assert_eq!(observations.len(), 2);
    assert!(observations[1][0] < observations[0][0] - 0.1);
    assert!(observations[1][1] < observations[0][1]);
    assert!(observations[1][2] < observations[0][2]);
    assert!(observations[1][3] < observations[0][3]);
}
