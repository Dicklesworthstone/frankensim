//! Native sphere convection uses bulk transport and the same fluid at the wall.
//! Synthetic property laws verify implementation, not experimental calibration.
use super::*;
use fs_airflow::graph::thermal::transport::transient::material::convection::{
    GraphConvectionError,
    reservoir::{GraphConvectiveReservoir, GraphReservoirConvectionError},
};
use fs_conduction::transient::deposition::{
    Reservoir,
    backward_euler::{BoundaryConductances, BoundaryStepCoefficients},
};

const DIAMETER: f64 = 0.04;
const PATCH_AREA: f64 = 0.002;

fn fixture(caloric: &ResolvedHeatCapacityCurve, wall_k: f64) -> DepositionSimulation<'_> {
    let original = thermal(caloric, [300.0; 6], Some(wall_k), false);
    let mut schedule: Vec<_> = original
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            if cell.id == CellId(3) {
                cell.boundary.reservoir = Some(Reservoir {
                    area_m2: PATCH_AREA,
                    conductance_w_k: 0.0,
                    temperature_k: wall_k,
                });
            } else if cell.id == CellId(6) {
                let volume = core::f64::consts::PI * DIAMETER.powi(3) / 6.0;
                cell.mass_kg = 1000.0 * volume;
                cell.geometry.volume_m3 = volume;
                cell.geometry.surface_area_m2 = core::f64::consts::PI * DIAMETER * DIAMETER;
                cell.geometry.end_m = [DIAMETER, 0.0, 0.0];
            }
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    schedule.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(original.links()[&9].clone()),
    });
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 16,
            max_links: 16,
            max_events: 16,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn bindings() -> (GraphConvectiveLink, GraphConvectiveReservoir) {
    (
        GraphConvectiveLink {
            link: 9,
            fluid_cell: CellId(2),
            characteristic_length: Length::new(DIAMETER),
            property_temperature: ConvectionPropertyTemperature::Reservoir,
            law: GraphConvectionLaw::ForcedSphere { branch: 0 },
        },
        GraphConvectiveReservoir {
            fluid_cell: CellId(3),
            characteristic_length: Length::new(DIAMETER),
            property_temperature: ConvectionPropertyTemperature::Reservoir,
            law: GraphConvectionLaw::ForcedSphere { branch: 0 },
        },
    )
}

// Independent Whitaker equation for the fixture's affine liquid viscosity.
// Only mu_surface changes at the wall; Re, Pr and k retain bulk properties.
fn coefficient(bulk_k: f64, wall_k: f64, speed: f64) -> f64 {
    let mu = 0.005 - 1e-5 * bulk_k;
    let surface_mu = 0.005 - 1e-5 * wall_k;
    let re = 1000.0 * speed * DIAMETER / mu;
    let pr = mu * 2000.0 / 0.5;
    let nu = 2.0
        + (0.4 * re.sqrt() + 0.06 * re.powf(2.0 / 3.0))
            * pr.powf(0.4)
            * (mu / surface_mu).powf(0.25);
    nu * 0.5 / DIAMETER
}

#[test]
fn g1_g3_native_sphere_uses_two_original_property_states_for_both_wall_owners() {
    context(|cx| {
        let card = card(true, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = fixture(&caloric, 400.0);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(2300.0);
        let owner = LiquidGraphSources::try_new(
            cx,
            &law,
            &caloric,
            &initial,
            SourceGraphBinding {
                graph: &graph,
                inventory: &inventory,
                laws: &laws,
                acceleration_m_s2: [0.0; 3],
            },
            caps(),
        )
        .unwrap();
        let (link, reservoir) = bindings();
        for sign in [-1.0, 1.0] {
            let q = sign * 1e-5;
            let flows = [q, 0.0, q];
            let resolved = owner.resolve(cx, &initial, &flows).unwrap();
            let report = native_report(
                cx,
                &resolved.graph,
                &inventory,
                flows,
                [sign, 0.0, sign],
                flows,
                0.1,
            );
            let links = owner
                .evaluate_links(cx, &initial, &resolved.sources, &report, &[link], 1)
                .unwrap();
            let reservoirs = owner
                .evaluate_reservoirs(cx, &initial, &resolved.sources, &report, &[reservoir], 1)
                .unwrap();
            let finite = &links.convection.samples()[0];
            let prescribed = &reservoirs.convection.samples()[0];
            for (material, thermal, surface, evaluation) in [
                (
                    &finite.material,
                    &finite.thermal,
                    finite.surface.as_ref().unwrap(),
                    finite.evaluation.as_ref().unwrap(),
                ),
                (
                    &prescribed.material,
                    &prescribed.thermal,
                    prescribed.surface.as_ref().unwrap(),
                    prescribed.evaluation.as_ref().unwrap(),
                ),
            ] {
                near(thermal.temperature_k(), 300.0);
                near(surface.thermal.temperature_k(), 400.0);
                near(thermal.transport().dynamic_viscosity_pa_s(), 0.002);
                near(surface.thermal.transport().dynamic_viscosity_pa_s(), 0.001);
                assert_eq!(material.material().card_identity(), card.content_hash());
                assert_eq!(
                    surface.material.material().card_identity(),
                    card.content_hash()
                );
                near(evaluation.nusselt.groups()["Re"], 200.0);
                near(evaluation.nusselt.groups()["Pr"], 8.0);
                near(evaluation.nusselt.groups()["mu_bulk_over_surface"], 2.0);
                near(
                    evaluation.coefficient.value.value(),
                    coefficient(300.0, 400.0, 0.01),
                );
            }
            near(finite.property_temperature_k, 300.0);
            near(prescribed.property_temperature_k, 300.0);
            near(
                finite.conductance_w_k,
                PATCH_AREA * coefficient(300.0, 400.0, 0.01),
            );
            near(finite.conductance_w_k, prescribed.conductance_w_k);
            near(finite.heat_into_fluid_w, 100.0 * finite.conductance_w_k);
            near(prescribed.heat_into_fluid_w, finite.heat_into_fluid_w);
        }
    });
}

#[test]
fn g4_sphere_interval_surface_and_geometry_refusals_preserve_native_history() {
    context(|cx| {
        let card = card(true, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = fixture(&caloric, 400.0);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(2300.0);
        let owner = LiquidGraphSources::try_new(
            cx,
            &law,
            &caloric,
            &initial,
            SourceGraphBinding {
                graph: &graph,
                inventory: &inventory,
                laws: &laws,
                acceleration_m_s2: [0.0; 3],
            },
            caps(),
        )
        .unwrap();
        let (link, reservoir) = bindings();
        let q = 1e-5;
        let resolved = owner.resolve(cx, &initial, &[q, 0.0, q]).unwrap();
        for old in [1e-7, -q] {
            // First Re_old=2 < 3.5. Second crosses zero with admissible old,
            // average and new speed magnitudes. The native drag admits both.
            let new = 2.0 * q - old;
            let drive = 1e6 * (new - old) / 0.1 + 1e5 * q;
            let report = native_report(
                cx,
                &resolved.graph,
                &inventory,
                [old, 0.0, old],
                [drive, 0.0, drive],
                [q, 0.0, q],
                0.1,
            );
            owner
                .admit_actual(cx, &initial, &resolved.sources, &report)
                .unwrap();
            assert!(matches!(
                owner.evaluate_links(cx, &initial, &resolved.sources, &report, &[link], 1),
                Err(GraphConvectionError::Envelope { link: 9, .. })
            ));
            assert!(matches!(
                owner.evaluate_reservoirs(
                    cx,
                    &initial,
                    &resolved.sources,
                    &report,
                    &[reservoir],
                    1
                ),
                Err(GraphReservoirConvectionError::Envelope {
                    cell: CellId(3),
                    ..
                })
            ));
        }
        let steady = native_report(
            cx,
            &resolved.graph,
            &inventory,
            [q, 0.0, q],
            [1.0, 0.0, 1.0],
            [q, 0.0, q],
            0.1,
        );
        let mut film = reservoir;
        film.property_temperature = ConvectionPropertyTemperature::FilmMean;
        assert!(matches!(
            owner.evaluate_reservoirs(cx, &initial, &resolved.sources, &steady, &[film], 1),
            Err(GraphReservoirConvectionError::Convection(_))
        ));
        let mut small = reservoir;
        small.characteristic_length = Length::new(0.01); // patch exceeds pi D²
        assert!(matches!(
            owner.evaluate_reservoirs(cx, &initial, &resolved.sources, &steady, &[small], 1),
            Err(GraphReservoirConvectionError::Convection(_))
        ));
        let hot = fixture(&caloric, 440.0); // mu_bulk/mu_wall = 10/3 > 3.2
        let hot_sources = owner.resolve(cx, &hot, &[q, 0.0, q]).unwrap();
        assert!(matches!(
            owner.evaluate_reservoirs(cx, &hot, &hot_sources.sources, &steady, &[reservoir], 1),
            Err(GraphReservoirConvectionError::Convection(_))
        ));
        // The bulk is supported while the required wall-property query is not.
        let limited = card_with_transport_ceiling(true, 1.0, 360.0);
        let (limited_law, limited_caloric) = source_pair(&limited);
        let limited_initial = fixture(&limited_caloric, 400.0);
        let limited_inventory = bind_inventory(cx, &graph, &limited_initial);
        let limited_owner = LiquidGraphSources::try_new(
            cx,
            &limited_law,
            &limited_caloric,
            &limited_initial,
            SourceGraphBinding {
                graph: &graph,
                inventory: &limited_inventory,
                laws: &laws,
                acceleration_m_s2: [0.0; 3],
            },
            caps(),
        )
        .unwrap();
        let supported_bulk = limited_owner
            .resolve(cx, &limited_initial, &[q, 0.0, q])
            .unwrap();
        let limited_report = native_report(
            cx,
            &supported_bulk.graph,
            &limited_inventory,
            [q, 0.0, q],
            [1.0, 0.0, 1.0],
            [q, 0.0, q],
            0.1,
        );
        assert!(matches!(
            limited_owner.evaluate_reservoirs(
                cx,
                &limited_initial,
                &supported_bulk.sources,
                &limited_report,
                &[reservoir],
                1
            ),
            Err(GraphReservoirConvectionError::Convection(_))
        ));
        near(initial.time_s(), 0.0);
        near(initial.temperature_k(CellId(2)).unwrap(), 300.0);
        near(initial.temperature_k(CellId(6)).unwrap(), 400.0);
        near(initial.links()[&9].conductance_w_k, 0.0);
    });
}

#[test]
fn g1_coupled_spheres_refresh_bulk_and_wall_viscosity_and_close_native_energy() {
    use fs_airflow::graph::thermal::transport::transient::inertial::{
        InertialGraphHeat, InertialThermalControls, InertialThermalProblem, ResolvedInertialGraph,
        advance_inertial_graph_distributed,
    };
    use fs_solver::NewtonKrylovConfig;
    context(|cx| {
        let card = card(true, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = fixture(&caloric, 400.0);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(2300.0);
        let owner = LiquidGraphSources::try_new(
            cx,
            &law,
            &caloric,
            &initial,
            SourceGraphBinding {
                graph: &graph,
                inventory: &inventory,
                laws: &laws,
                acceleration_m_s2: [0.0; 3],
            },
            caps(),
        )
        .unwrap();
        let initial_flows = [1e-5, 0.0, 1e-5];
        let first = owner.resolve(cx, &initial, &initial_flows).unwrap();
        let branches: Vec<_> = inventory
            .inertial_branches()
            .iter()
            .map(|b| b.unwrap())
            .collect();
        let boundaries = vec![FixedPressure {
            node: 0,
            pressure: Pressure::new(0.0),
        }];
        let momentum = InertialGraph::try_new(
            cx,
            &first.graph,
            Density::new(1000.0),
            &branches,
            &boundaries,
            &initial_flows,
            0.0,
            32,
        )
        .unwrap();
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: inventory.node_cells(),
            maximum_connections: 32,
        };
        let controls = InertialThermalControls {
            solver: NewtonKrylovConfig {
                absolute_tolerance: 1e-9,
                relative_tolerance: 0.0,
                ..NewtonKrylovConfig::default()
            },
            maximum_cells: 16,
            maximum_newton_iterations: 32,
            maximum_endpoint_evaluations: 2048,
            energy_scale_j: 1.0,
            flow_scale_m3_s: 1e-5,
            finite_difference_step: 1e-4,
            cell_tolerance_j: 1e-6,
            total_thermal_tolerance_j: 1e-6,
            impulse_tolerance_pa_s: 1e-6,
            mechanical_tolerance_j: 1e-8,
            average_flow_tolerance_m3_s: 1e-12,
            momentum: InertialGraphControls {
                maximum_drag_step: 2.0,
                maximum_neglected_power_w: 100.0,
                maximum_neglected_energy_j: 100.0,
            },
        };
        let (link, reservoir) = bindings();
        let dt = 1.0;
        let (next_momentum, next, report) = advance_inertial_graph_distributed(
            cx,
            &problem,
            &inventory,
            &initial,
            dt,
            controls,
            &[0.0; 7],
            &initial_flows,
            |cx, candidate, flows| {
                let resolved = owner
                    .resolve(cx, candidate, flows)
                    .map_err(|e| e.to_string())?;
                Ok::<_, String>(ResolvedInertialGraph {
                    graph: resolved.graph,
                    pressure_boundaries: boundaries.clone(),
                    driving_pressure_pa: vec![1.0, 0.0, 1.0],
                    payload: resolved.sources,
                })
            },
            |cx, candidate, sources, hydraulics| {
                let links = owner
                    .evaluate_links(cx, candidate, sources, hydraulics, &[link], 1)
                    .map_err(|e| e.to_string())?;
                let reservoirs = owner
                    .evaluate_reservoirs(cx, candidate, sources, hydraulics, &[reservoir], 1)
                    .map_err(|e| e.to_string())?;
                let mut conductances = BTreeMap::new();
                let mut boundaries = BoundaryConductances::default();
                links
                    .convection
                    .apply_to(&mut conductances, || cx.checkpoint().is_err())
                    .map_err(|e| e.to_string())?;
                reservoirs
                    .convection
                    .apply_to(&mut boundaries, || cx.checkpoint().is_err())
                    .map_err(|e| e.to_string())?;
                Ok::<_, String>(InertialGraphHeat {
                    mass_boundaries: vec![],
                    heat: BoundaryStepCoefficients {
                        bulk: StepCoefficients {
                            sources_w: BTreeMap::new(),
                            conductances_w_k: conductances,
                            payload: (links, reservoirs),
                        },
                        boundaries,
                    },
                })
            },
        )
        .unwrap();
        let (links, reservoirs) = &report.heat;
        let finite = &links.convection.samples()[0];
        let prescribed = &reservoirs.convection.samples()[0];
        let q = report.hydraulics.average_volume_flows_m3_s[0];
        assert!(q > initial_flows[0]);
        assert!(next.temperature_k(CellId(6)).unwrap() < 400.0);
        for (cell, bulk, surface, actual_h) in [
            (
                finite.fluid_cell,
                &finite.material,
                finite.surface.as_ref().unwrap(),
                finite
                    .evaluation
                    .as_ref()
                    .unwrap()
                    .coefficient
                    .value
                    .value(),
            ),
            (
                prescribed.fluid_cell,
                &prescribed.material,
                prescribed.surface.as_ref().unwrap(),
                prescribed
                    .evaluation
                    .as_ref()
                    .unwrap()
                    .coefficient
                    .value
                    .value(),
            ),
        ] {
            let bulk_k = next.temperature_k(cell).unwrap();
            let wall_k = if cell == CellId(2) {
                next.temperature_k(CellId(6)).unwrap()
            } else {
                400.0
            };
            assert!(bulk_k > 300.0);
            near(bulk.dynamic_viscosity_pa_s(), 0.005 - 1e-5 * bulk_k);
            near(surface.thermal.temperature_k(), wall_k);
            near(
                surface.thermal.transport().dynamic_viscosity_pa_s(),
                0.005 - 1e-5 * wall_k,
            );
            near(actual_h, coefficient(bulk_k, wall_k, q / PIPE_AREA));
            assert_eq!(
                surface.material.material().card_identity(),
                card.content_hash()
            );
        }
        let frozen = coefficient(300.0, 400.0, initial_flows[0] / PIPE_AREA);
        assert!(
            (finite
                .evaluation
                .as_ref()
                .unwrap()
                .coefficient
                .value
                .value()
                - frozen)
                .abs()
                > 1e-5
        );
        let stored: f64 = next
            .cells()
            .iter()
            .map(|(id, cell)| {
                cell.mass_kg * (cell.specific_h_j_kg - initial.cells()[id].specific_h_j_kg)
            })
            .sum();
        let wall_change = next.cells()[&CellId(6)].mass_kg
            * (next.cells()[&CellId(6)].specific_h_j_kg
                - initial.cells()[&CellId(6)].specific_h_j_kg);
        let external = dt * prescribed.heat_into_fluid_w;
        assert!((stored - external).abs() < 2e-6);
        assert!((report.thermal.external_heat_j - external).abs() < 2e-6);
        assert!((wall_change + dt * finite.heat_into_fluid_w).abs() < 1e-6);
        assert!(report.maximum_cell_residual_j <= controls.cell_tolerance_j);
        assert!(
            report.hydraulics.mechanical_balance_residual_j.abs()
                <= controls.mechanical_tolerance_j
        );
        near(next.time_s(), dt);
        near(next_momentum.accepted_time_s(), dt);
        near(initial.time_s(), 0.0);
        near(momentum.accepted_time_s(), 0.0);
    });
}
