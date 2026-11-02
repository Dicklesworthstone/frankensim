//! Native tetrahedra and contact traces share one finite source-fluid endpoint.
//! Manufactured data check the equations; these are not empirical materials.
use super::*;
use fs_couple::circuit_enthalpy::{
    conduction::{
        contacts::ConductionContactThermal,
        heat::{ConductionHeatError, step_spatial_conduction_candidate},
        spatial::{SpatialConductionEnthalpyCandidate, SpatialConductionRegion},
    },
    contacts::heat::HeatAssembly,
};
use fs_em::{conduction::ConductionControls, spatial::conduction::*};
use fs_exec::Cx;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::{
    conductor::{ResolvedConductorState, resolve_conductor_state_with_basis},
    electrical_interface::{
        ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
        resolve_electrical_interface_state,
    },
    thermal::interface::ThermalInterfacePoint,
};
use fs_solver::CheckedCgConfig;

const CROSS_SECTION: f64 = 1e-5;
const VOLTAGE: f64 = 3.0;

fn solid_point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_quantity(), t)
        .unwrap()
}

fn solid_material(resistivity_scale: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    // This field fixture has only owned temperature coordinates. The parent
    // graph fixture's untyped lot tag is not a fixed physical source coordinate.
    let (original, _) = body_material();
    let mut claims = ClaimSet::new();
    for (_, original) in original.claims().claims_ordered() {
        let mut claim = original.clone();
        claim.validity = ValidityDomain::unconstrained().with_quantity(
            "T",
            temperature_quantity(),
            280.0,
            320.0,
        );
        if claim.key.name() == "electrical_resistivity" {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                for (_, value) in knots {
                    *value *= resistivity_scale;
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
            lower: &solid_point(280.0),
            upper: &solid_point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}

fn interface_point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("Ta", temperature_quantity(), t)
        .unwrap()
        .with_quantity("Tb", temperature_quantity(), t)
        .unwrap()
}

fn interface(card: &MaterialCard, scale: f64) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
                ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Surface2D {
                axes: ["Ta".into(), "Tb".into()],
                axis_dims: [Temperature::DIMS; 2],
                samples: [
                    [280.0, 280.0],
                    [320.0, 280.0],
                    [280.0, 320.0],
                    [320.0, 320.0],
                ]
                .map(|p| {
                    (
                        p,
                        scale
                            * 0.5
                            * CROSS_SECTION
                            * (1.0 + 0.003 * (p[0] - 300.0) + 0.007 * (p[1] - 300.0)),
                    )
                })
                .to_vec(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("Ta", temperature_quantity(), 280.0, 320.0)
                .with_quantity("Tb", temperature_quantity(), 280.0, 320.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured two-sided finite-coolant interface, not empirical data"
                    .into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: card.id().clone(),
            texture_frame: "first original face".into(),
        },
        SurfaceSpec {
            material: card.id().clone(),
            texture_frame: "second original face".into(),
        },
        SystemContext {
            medium: "manufactured fixed contact".into(),
            third_body: None,
            environment: "fixed external conditions".into(),
            history: "fixed preparation".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 20,
        max_dofs: 100,
        max_triplets: 1000,
    }
}

fn specimen(card: &MaterialCard, interface: &InterfaceSystemCard) -> SpatialConduction {
    let mut positions = Vec::new();
    let mut tetrahedra = Vec::new();
    for (part, xs) in [[0.0, 0.05], [0.05, 0.1]].into_iter().enumerate() {
        for x in xs {
            for z in [0.0, 0.001] {
                for y in [0.0, 0.01] {
                    positions.push([x, y, z]);
                }
            }
        }
        tetrahedra.extend(
            [
                [0, 4, 5, 7],
                [0, 5, 1, 7],
                [0, 1, 3, 7],
                [0, 3, 2, 7],
                [0, 2, 6, 7],
                [0, 6, 4, 7],
            ]
            .map(|tet| tet.map(|vertex| vertex + 8 * part as u32)),
        );
    }
    let bulk = resolve_conductor_state_with_basis(
        card,
        &solid_point(310.0),
        MaterialPropertySelection::SingleClaimOnly,
        ElectricalPropertyBasis::Resistivity,
    )
    .unwrap();
    let surface = resolve_electrical_interface_state(
        interface,
        &interface_point(310.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    SpatialConduction::new_with_contacts(
        positions,
        tetrahedra,
        vec![bulk; 12],
        [([4, 5, 7], [8, 9, 11]), ([4, 6, 7], [8, 10, 11])]
            .map(
                |(first_vertices, second_vertices)| SpatialConductionContact {
                    first_vertices,
                    second_vertices,
                    source: surface.clone(),
                },
            )
            .to_vec(),
        limits(),
        &mut || false,
    )
    .unwrap()
}

fn boundary(problem: &SpatialConduction) -> SpatialConductionBoundary {
    let mut boundary = SpatialConductionBoundary::default();
    for face in problem.exterior_faces() {
        let xs = problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0]);
        if xs == [0.0; 3] || xs == [0.1; 3] {
            boundary.voltage.push(ConductionVoltageFace {
                face,
                voltage_v: [if xs[0] == 0.0 { 0.0 } else { VOLTAGE }; 3],
            });
        }
    }
    boundary
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn advance<'t, A: HeatAssembly + ?Sized>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    problem: &SpatialConduction,
    boundary: &SpatialConductionBoundary,
    contacts: &[ConductionContactThermal<'_>],
    regions: &[SpatialConductionRegion],
) -> Result<
    SpatialConductionEnthalpyCandidate<'t, ResolvedConductorState, A::Sample>,
    CircuitEnthalpyError<ConductionHeatError<A::Error>>,
> {
    let mut solve = controls();
    solve.maximum_cells = 2;
    solve.newton.linear_restart = 2;
    step_spatial_conduction_candidate(
        heat,
        cx,
        thermal,
        thermal.time_s() + DT,
        problem,
        boundary,
        contacts,
        regions,
        1e-14,
        limits(),
        ConductionControls {
            linear: CheckedCgConfig {
                tolerance: 1e-12,
                max_iterations: 200,
                max_corrections: 2,
            },
            current_tolerance_a: 1e-9,
            power_tolerance_w: 1e-9,
        },
        solve,
        &[0.0, 0.0],
    )
}

// Independent series-slab solution plus two backward-Euler caloric balances.
// No production field, convection evaluator or coupled residual is called.
fn reference(
    old: [f64; 2],
    viscosity: f64,
    fluid_cp: f64,
    bulk_scale: f64,
    contact_scale: f64,
) -> [f64; 6] {
    let evaluate = |wall: f64| {
        let bulk_r = bulk_scale * (1.0 + 0.02 * (wall - 300.0));
        let contact_r = contact_scale * 0.5 * (1.0 + 0.01 * (wall - 300.0));
        let current = VOLTAGE / (bulk_r + contact_r);
        let joule = DT * VOLTAGE * current;
        let fluid = old[1] + (joule - MASS * CP * (wall - old[0])) / (MASS * fluid_cp);
        let film = 0.5 * (wall + fluid);
        let k = 0.5 + 0.002 * (film - 300.0);
        let mu = viscosity * 0.001 * (1.0 - 0.01 * (film - 300.0));
        let h = 0.664 * (1000.0 * 0.1 / mu).sqrt() * (mu * fluid_cp / k).cbrt() * k / 0.1;
        let convection = SURFACE * h * (wall - fluid);
        (
            MASS * CP * (wall - old[0]) - joule + DT * convection,
            [wall, fluid, current, convection, bulk_r, contact_r],
        )
    };
    let (mut lower, mut upper) = (290.0, 312.0);
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
fn g1_g3_native_current_and_contacts_share_finite_source_coolant() {
    let run = |viscosity, fluid_cp, bulk_scale, contact_scale, reversed| {
        let (body_card, body_curve) = solid_material(bulk_scale);
        let (fluid_card, fluid_curve) = fluid_material(viscosity, fluid_cp);
        let body_context = solid_point(310.0);
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
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &[],
        };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &[],
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
        let interface = interface(&body_card, contact_scale);
        let contact_context = interface_point(310.0);
        let contacts = [0, 1].map(|contact| ConductionContactThermal {
            contact,
            point: ThermalInterfacePoint {
                interface: &interface,
                materials: [&body_card, &body_card],
                context: &contact_context,
                temperature_axes: ["Ta", "Tb"],
                liquid_fraction_axes: [None, None],
                mean_temperature_axis: None,
            },
            heat_fraction_to_first: 0.25,
        });
        let problem = specimen(&body_card, &interface);
        let boundary = boundary(&problem);
        // Both disjoint half-slabs share one explicit uniform caloric body;
        // the finite coolant has no electrical tetrahedra and retains its mass.
        let regions: Vec<_> = (0..12)
            .rev()
            .map(|tetrahedron| SpatialConductionRegion {
                tetrahedron,
                cell: BODY,
            })
            .collect();
        let gate = fs_exec::CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let mut accepted = [0.0; 3];
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 181,
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
            let next = advance(
                &assembly, &cx, &thermal, &problem, &boundary, &contacts, &regions,
            )
            .unwrap();
            let expected = reference(
                [310.0, 290.0],
                viscosity,
                fluid_cp,
                bulk_scale,
                contact_scale,
            );
            let report = next.report();
            let field = &report.electrical.field;
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
            for density in &field.current_density_a_m2 {
                near(density[0] * CROSS_SECTION, -expected[2], 1e-8);
                near(density[1] * CROSS_SECTION, 0.0, 1e-8);
                near(density[2] * CROSS_SECTION, 0.0, 1e-8);
            }
            near(
                field.cell_power_w.iter().sum(),
                expected[2].powi(2) * expected[4],
                1e-8,
            );
            for (contact, source) in report
                .electrical
                .contacts
                .iter()
                .zip(next.problem().contacts())
            {
                assert_eq!(contact.cells, [BODY, BODY]);
                near(contact.response.current_a, -0.5 * expected[2], 1e-8);
                near(
                    contact.response.power_w,
                    0.5 * expected[2].powi(2) * expected[5],
                    1e-8,
                );
                near(
                    source.source.resistance_ohm_m2(),
                    expected[5] * CROSS_SECTION,
                    1e-14,
                );
                for axis in ["Ta", "Tb"] {
                    near(
                        source.source.resolved().state_coordinate(axis).unwrap(),
                        expected[0],
                        1e-8,
                    );
                }
            }
            let finite = report.source.finite_convection.as_ref().unwrap();
            let sample = &finite.samples()[0];
            assert_eq!(sample.wall_cell, BODY);
            assert_eq!(sample.fluid_cell, FLUID);
            near(sample.heat_into_fluid_w, expected[3], 1e-7);
            near(
                sample.property_temperature_k,
                0.5 * (expected[0] + expected[1]),
                1e-8,
            );
            assert_eq!(
                sample.material.material().card_identity(),
                fluid_card.content_hash()
            );
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
            near(report.deposited_joule_j, DT * VOLTAGE * expected[2], 1e-8);
            near(report.non_electrical_heat_j, 0.0, 1e-12);
            near(
                next.thermal().energy_j() - thermal.energy_j(),
                report.electrical.supplied_work_j,
                1e-8,
            );
            near(report.total_energy_residual_j, 0.0, 1e-8);
            assert_eq!(report.source.heat.material.cells.len(), 2);
            assert!(report.source.heat.material.links.is_empty());
            assert_eq!(next.thermal().links(), &saved_links);
            accepted = [expected[0], expected[1], field.power_w / VOLTAGE];

            // Continue with exactly the same fluid/source pins and the jointly
            // returned native field plus caloric history, not a fresh field.
            let continued = advance(
                &assembly,
                &cx,
                next.thermal(),
                next.problem(),
                &boundary,
                &contacts,
                &regions,
            )
            .unwrap();
            let second = reference(
                [expected[0], expected[1]],
                viscosity,
                fluid_cp,
                bulk_scale,
                contact_scale,
            );
            near(
                continued.thermal().temperature_k(BODY).unwrap(),
                second[0],
                1e-8,
            );
            near(
                continued.thermal().temperature_k(FLUID).unwrap(),
                second[1],
                1e-8,
            );
            near(
                continued.report().electrical.field.power_w,
                VOLTAGE * second[2],
                1e-8,
            );
            near(
                continued.thermal().energy_j() - next.thermal().energy_j(),
                continued.report().electrical.supplied_work_j,
                1e-8,
            );

            // Omitting a finite link's owner cannot silently remove coolant.
            assert!(
                advance(
                    &prepared, &cx, &thermal, &problem, &boundary, &contacts, &regions
                )
                .is_err()
            );
            let refused_network = SourceConvectionHeatNetwork {
                heat: &prepared,
                cells: &[],
                maximum_cells: 0,
                resolve_fluid: |_: &SourceConvectiveCell,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    Ok::<_, &'static str>(fluid.at_temperature(t).unwrap())
                },
            };
            let reached_trial = std::cell::Cell::new(false);
            let refused =
                refused_network.with_links(&links, 1, |_, trial: &DepositionSimulation<'_>, t| {
                    if trial.temperature_k(BODY).unwrap() < 309.99 {
                        reached_trial.set(true);
                        Err("original coolant source refused")
                    } else {
                        Ok(fluid.at_temperature(t).unwrap())
                    }
                });
            assert!(matches!(
                advance(
                    &refused, &cx, &thermal, &problem, &boundary, &contacts, &regions
                ),
                Err(CircuitEnthalpyError::Source(ConductionHeatError::Heat(
                    ConvectionHeatError::Convection(SourceCellConvectionError::Source {
                        cell: FLUID,
                        source: "original coolant source refused",
                        ..
                    })
                )))
            ));
            assert!(reached_trial.get());
            assert_eq!(thermal.cells(), &saved);
            assert_eq!(thermal.links(), &saved_links);
            for contact in problem.contacts() {
                near(
                    contact.source.resolved().state_coordinate("Ta").unwrap(),
                    310.0,
                    0.0,
                );
            }
        });
        accepted
    };
    let original = run(1.0, 4000.0, 1.0, 1.0, false);
    let reversed = run(1.0, 4000.0, 1.0, 1.0, true);
    for (actual, expected) in reversed.into_iter().zip(original) {
        near(actual, expected, 1e-9);
    }
    let viscous = run(4.0, 4000.0, 1.0, 1.0, false);
    let larger_capacity = run(1.0, 8000.0, 1.0, 1.0, false);
    let other_solid = run(1.0, 4000.0, 2.0, 1.0, false);
    let other_interface = run(1.0, 4000.0, 1.0, 2.0, false);
    assert!(viscous[0] > original[0] && viscous[2] < original[2]);
    assert!(larger_capacity[0] < original[0] && larger_capacity[2] > original[2]);
    assert!(other_solid[2] < original[2] && other_interface[2] < original[2]);
}
