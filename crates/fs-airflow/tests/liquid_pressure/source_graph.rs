//! Original source/caloric identity and physical graph closure; synthetic data
//! verify implementation and do not claim experimentally calibrated accuracy.
use super::*;
use fs_airflow::graph::inertial::{InertialGraph, InertialGraphControls};
use fs_airflow::graph::thermal::medium::convection::{
    ForcedConvectionSettings, deposition::source::ConvectionPropertyTemperature,
};
use fs_airflow::graph::thermal::natural_circulation::material_response::MaterialBranchLaw;
use fs_airflow::graph::thermal::transport::transient::{
    distributed::DistributedGraphInventory,
    material::{
        LiquidGraphSources, SourceGraphBinding, SourceGraphError, SourceGraphLimits,
        convection::{GraphConvectionLaw, GraphConvectiveLink},
    },
};
use fs_airflow::graph::{FixedPressure, GraphBranch, GraphSolveConfig, LossGraph};
use fs_airflow::{LossElement, LossResistance, SourceProvenance, ToleranceBasis};
use fs_conduction::transient::deposition::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
    backward_euler::{BackwardEulerStep, StepCoefficients},
};
use fs_qty::{Density, VolumetricFlowRate};
use std::collections::BTreeMap;

const PIPE_AREA: f64 = 1e-3;
const A: f64 = 5e7;

// Authored developed-flow support. Subdividing a passage does not reset its
// upstream development length; the original correlation still admits L/D.
fn developed_settings() -> ForcedConvectionSettings {
    ForcedConvectionSettings {
        length_over_hydraulic_diameter: Some(100.0),
        ..ForcedConvectionSettings::default()
    }
}

fn context(run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 104,
                kernel_id: 756,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}
fn source_pair(card: &MaterialCard) -> (LiquidTemperatureLaw<'_>, ResolvedHeatCapacityCurve) {
    let law = LiquidTemperatureLaw::try_new(
        card,
        &point(300.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let cp = PropertyKey::with_quantity(
        "specific_heat_capacity",
        QuantitySpec::dimensional(Dims([2, 0, -2, -1, 0, 0])),
    );
    let rho = PropertyKey::with_quantity(
        "density",
        QuantitySpec::dimensional(Dims([-3, 1, 0, 0, 0, 0])),
    );
    let caloric = ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &cp,
            density: &rho,
            temperature_axis: "T",
            lower: &point(250.0),
            upper: &point(450.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        0.0,
        0.5,
    )
    .unwrap();
    (law, caloric)
}
fn graph() -> LossGraph {
    let branches = [(0, 1, "split-a"), (0, 1, "split-b"), (1, 0, "return")]
        .into_iter()
        .map(|(from, to, name)| GraphBranch {
            from,
            to,
            loss: LossElement::new(
                name,
                LossResistance::new(1.0),
                0.0,
                SourceProvenance::new("synthetic native graph", name),
                ToleranceBasis::Analytic,
            )
            .unwrap(),
        })
        .collect();
    LossGraph::new(2, branches).unwrap()
}
fn paths() -> Vec<Vec<CellId>> {
    vec![vec![CellId(2), CellId(3)], vec![CellId(4)], vec![CellId(5)]]
}
fn laws(re_max: f64) -> Vec<Vec<MaterialBranchLaw>> {
    [vec![0.5, 0.5], vec![1.0], vec![1.0]]
        .into_iter()
        .map(|weights| {
            weights
                .into_iter()
                .map(|weight| MaterialBranchLaw {
                    resistance: MaterialFlowResistance::try_new(
                        A * weight,
                        0.0,
                        Area::new(PIPE_AREA),
                        Length::new(0.03),
                        (0.0, re_max),
                    )
                    .unwrap(),
                    source: SourceProvenance::new("authored straight section", "native-section"),
                    uncertainty_rel: 0.0,
                    tolerance_basis: ToleranceBasis::Analytic,
                    exchange: vec![],
                })
                .collect()
        })
        .collect()
}
fn caps() -> SourceGraphLimits {
    SourceGraphLimits {
        maximum_cells: 16,
        maximum_relative_density_departure: 0.15,
        maximum_neglected_section_power_w: 100.0,
        maximum_neglected_section_energy_j: 100.0,
    }
}
fn thermal<'a>(
    caloric: &'a ResolvedHeatCapacityCurve,
    temperatures: [f64; 6],
    wall: Option<f64>,
    bent: bool,
) -> DepositionSimulation<'a> {
    let middle = if bent {
        [0.2, 0.0, 0.5]
    } else {
        [0.0, 0.0, 0.5]
    };
    let geometry = [
        ([0.0, 0.0, 0.0], [0.1, 0.0, 0.0], 0.1),
        ([0.0, 0.0, 1.0], [0.1, 0.0, 1.0], 0.1),
        ([0.0, 0.0, 0.0], middle, 0.5),
        (middle, [0.0, 0.0, 1.0], 0.5),
        ([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0),
        ([0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 1.0),
    ];
    let mut schedule: Vec<_> = geometry
        .into_iter()
        .zip(temperatures)
        .enumerate()
        .map(|(i, ((start, end, length), t))| {
            let volume = PIPE_AREA * length;
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(i as u64),
                    segment_id: i as u64,
                    body_id: 0,
                    material_id: MaterialId(1),
                    mass_kg: 1000.0 * volume,
                    specific_h_j_kg: 2000.0 * (t - 250.0),
                    geometry: CellGeometry {
                        start_m: start,
                        end_m: end,
                        volume_m3: volume,
                        surface_area_m2: 0.1,
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
        .collect();
    if let Some(t) = wall {
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(6),
                segment_id: 6,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: 1.0,
                specific_h_j_kg: 2000.0 * (t - 250.0),
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [0.0, 0.1, 0.0],
                    volume_m3: 0.001,
                    surface_area_m2: 0.1,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 9,
                a: CellId(2),
                b: CellId(6),
                area_m2: 0.002,
                conductance_w_k: 0.0,
            }),
        });
    }
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
fn bind_inventory(
    cx: &Cx<'_>,
    graph: &LossGraph,
    thermal: &DepositionSimulation<'_>,
) -> DistributedGraphInventory {
    DistributedGraphInventory::new(
        cx,
        graph,
        Density::new(1000.0),
        &[CellId(0), CellId(1)],
        &paths(),
        thermal,
        32,
    )
    .unwrap()
}
fn native_report(
    cx: &Cx<'_>,
    graph: &LossGraph,
    inventory: &DistributedGraphInventory,
    initial: [f64; 3],
    drive: [f64; 3],
    candidate: [f64; 3],
    dt: f64,
) -> fs_airflow::graph::inertial::InertialGraphReport {
    let branches: Vec<_> = inventory
        .inertial_branches()
        .iter()
        .map(|b| b.unwrap())
        .collect();
    let boundaries = [FixedPressure {
        node: 0,
        pressure: Pressure::new(0.0),
    }];
    let momentum = InertialGraph::try_new(
        cx,
        graph,
        Density::new(1000.0),
        &branches,
        &boundaries,
        &initial,
        0.0,
        32,
    )
    .unwrap();
    momentum
        .propose(
            cx,
            graph,
            &boundaries,
            &drive,
            &candidate,
            dt,
            InertialGraphControls {
                maximum_drag_step: 2.0,
                maximum_neglected_power_w: 1e6,
                maximum_neglected_energy_j: 1e6,
            },
        )
        .unwrap()
        .accept(cx, 1e-6, 1e-6)
        .unwrap()
        .1
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-9 * a.abs().max(b.abs()).max(1e-9),
        "{a} != {b}"
    );
}

#[test]
fn original_spatial_viscosity_controls_branched_flow_and_reversal() {
    context(|cx| {
        let card = card(true, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = thermal(
            &caloric,
            [300.0, 300.0, 300.0, 400.0, 300.0, 300.0],
            None,
            false,
        );
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
        for sign in [-1.0, 1.0] {
            let resolved = owner
                .resolve(cx, &initial, &[sign * 1e-5, sign * 7.5e-6, -sign * 7.5e-6])
                .unwrap();
            let result = resolved
                .graph
                .solve(
                    &[
                        FixedPressure {
                            node: 0,
                            pressure: Pressure::new(sign * 0.75),
                        },
                        FixedPressure {
                            node: 1,
                            pressure: Pressure::new(0.0),
                        },
                    ],
                    GraphSolveConfig {
                        max_sweeps: 128,
                        max_node_iterations: 80,
                        absolute_flow_tolerance: VolumetricFlowRate::new(1e-12),
                        relative_flow_tolerance: 1e-10,
                    },
                    cx,
                )
                .unwrap();
            // mu(300)=.002, mu(400)=.001; both half-length sections share Q.
            near(result.branches[0].flow.value(), sign * 1e-5);
            near(result.branches[1].flow.value(), sign * 7.5e-6);
            near(result.branches[2].flow.value(), -sign * 7.5e-6);
            let sections = &resolved.sources.sections()[0];
            near(sections[0].momentum.medium.dynamic_viscosity_pa_s(), 0.002);
            near(sections[1].momentum.medium.dynamic_viscosity_pa_s(), 0.001);
            assert_eq!(
                sections[1].momentum.medium.material().card_identity(),
                law.reference_state().material().card_identity()
            );
            assert!(
                initial
                    .cells()
                    .values()
                    .all(|c| initial.uses_material_law(c.material_id, caloric.curve()))
            );
        }
    });
}

#[test]
fn cancelling_section_buoyancy_cannot_hide_local_omitted_work() {
    context(|cx| {
        let card = buoyant_card();
        let (law, caloric) = source_pair(&card);
        let initial = thermal(
            &caloric,
            [300.0, 300.0, 280.0, 320.0, 300.0, 300.0],
            None,
            false,
        );
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(2300.0);
        let mut limits = caps();
        limits.maximum_neglected_section_power_w = 0.005;
        let owner = LiquidGraphSources::try_new(
            cx,
            &law,
            &caloric,
            &initial,
            SourceGraphBinding {
                graph: &graph,
                inventory: &inventory,
                laws: &laws,
                acceleration_m_s2: [0.0, 0.0, -9.81],
            },
            limits,
        )
        .unwrap();
        let resolved = owner.resolve(cx, &initial, &[1e-4, 0.0, 1e-4]).unwrap();
        near(resolved.sources.branch_densities_kg_m3()[0], 1000.0);
        let heads = &resolved.sources.section_anomalous_heads_pa()[0];
        near(heads[0], -49.05);
        near(heads[1], 49.05);
        let report = native_report(
            cx,
            &resolved.graph,
            &inventory,
            [1e-4, 0.0, 1e-4],
            [10.0, 0.0, 10.0],
            [1e-4, 0.0, 1e-4],
            0.1,
        );
        // Net branch buoyancy cancels. Native and source ledgers both retain
        // .00981 W of local work, plus .002 W of dissipative power.
        near(report.neglected_buoyancy_energy_bound_j, 0.000981);
        match owner.admit_actual(cx, &initial, &resolved.sources, &report) {
            Err(SourceGraphError::NeglectedWork { power_w, energy_j }) => {
                near(power_w, 0.01181);
                near(energy_j, 0.001181);
            }
            other => panic!("local source-work cap did not refuse: {other:?}"),
        }
    });
}

#[test]
fn actual_endpoint_envelope_refuses_reynolds_hidden_by_reversal() {
    context(|cx| {
        let card = card(false, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = thermal(&caloric, [300.0; 6], None, false);
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
        let resolved = owner.resolve(cx, &initial, &[0.0; 3]).unwrap();
        let report = native_report(
            cx,
            &resolved.graph,
            &inventory,
            [1e-3, 0.0, 1e-3],
            [-20000.0, 0.0, -20000.0],
            [0.0; 3],
            0.1,
        );
        assert!(report.average_volume_flows_m3_s[0].abs() < 1e-12);
        assert!(matches!(
            owner.admit_actual(cx, &initial, &resolved.sources, &report),
            Err(SourceGraphError::Resistance {
                branch: 0,
                section: 0,
                ..
            })
        ));
    });
}

#[test]
fn exact_caloric_owner_source_domain_and_stale_candidates_are_not_substituted() {
    context(|cx| {
        let card = card(true, 1.0);
        let (law, caloric) = source_pair(&card);
        let (_, other_caloric) = source_pair(&card);
        let initial = thermal(&caloric, [300.0; 6], None, false);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(2300.0);
        let binding = || SourceGraphBinding {
            graph: &graph,
            inventory: &inventory,
            laws: &laws,
            acceleration_m_s2: [0.0; 3],
        };
        assert!(matches!(
            LiquidGraphSources::try_new(cx, &law, &other_caloric, &initial, binding(), caps()),
            Err(SourceGraphError::Cell { .. })
        ));
        let owner =
            LiquidGraphSources::try_new(cx, &law, &caloric, &initial, binding(), caps()).unwrap();
        let resolved = owner.resolve(cx, &initial, &[0.0; 3]).unwrap();
        let report = native_report(
            cx,
            &resolved.graph,
            &inventory,
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
            0.1,
        );
        let changed = thermal(
            &caloric,
            [300.0, 300.0, 310.0, 300.0, 300.0, 300.0],
            None,
            false,
        );
        assert!(matches!(
            owner.admit_actual(cx, &changed, &resolved.sources, &report),
            Err(SourceGraphError::Invalid(_))
        ));
        let foreign =
            LiquidGraphSources::try_new(cx, &law, &caloric, &initial, binding(), caps()).unwrap();
        assert!(matches!(
            foreign.admit_actual(cx, &initial, &resolved.sources, &report),
            Err(SourceGraphError::Invalid(_))
        ));

        let narrow = card_with_transport_ceiling(true, 1.0, 315.0);
        let (narrow_law, narrow_caloric) = source_pair(&narrow);
        let hot = thermal(
            &narrow_caloric,
            [300.0, 300.0, 350.0, 300.0, 300.0, 300.0],
            None,
            false,
        );
        let hot_inventory = bind_inventory(cx, &graph, &hot);
        let narrow_owner = LiquidGraphSources::try_new(
            cx,
            &narrow_law,
            &narrow_caloric,
            &hot,
            SourceGraphBinding {
                graph: &graph,
                inventory: &hot_inventory,
                laws: &laws,
                acceleration_m_s2: [0.0; 3],
            },
            caps(),
        )
        .unwrap();
        assert!(matches!(
            narrow_owner.resolve(cx, &hot, &[0.0; 3]),
            Err(SourceGraphError::Liquid(_))
        ));
    });
}

#[test]
fn bent_original_source_path_drives_circulation_that_an_endpoint_mean_would_erase() {
    context(|cx| {
        let card = buoyant_card();
        let (law, caloric) = source_pair(&card);
        let initial = thermal(
            &caloric,
            [300.0, 300.0, 320.0, 280.0, 300.0, 300.0],
            None,
            true,
        );
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let mut laws = laws(2300.0);
        let cells = paths();
        for (section, id) in laws[0].iter_mut().zip(&cells[0]) {
            let cell = &initial.cells()[id];
            let delta: [f64; 3] = core::array::from_fn(|axis| {
                cell.geometry.end_m[axis] - cell.geometry.start_m[axis]
            });
            let length = fs_math::det::sqrt(delta.iter().map(|d| d * d).sum());
            near(cell.mass_kg, 1000.0 * cell.geometry.volume_m3);
            section.resistance = MaterialFlowResistance::try_new(
                A * 0.5,
                0.0,
                Area::new(cell.geometry.volume_m3 / length),
                Length::new(0.03),
                (0.0, 2300.0),
            )
            .unwrap();
        }
        let owner = LiquidGraphSources::try_new(
            cx,
            &law,
            &caloric,
            &initial,
            SourceGraphBinding {
                graph: &graph,
                inventory: &inventory,
                laws: &laws,
                acceleration_m_s2: [2.0, 0.0, -9.81],
            },
            caps(),
        )
        .unwrap();
        let resolved = owner.resolve(cx, &initial, &[0.0; 3]).unwrap();
        let gravity = resolved.graph.gravity().unwrap();
        let sections = gravity.branch_sections().unwrap();
        for (section, id) in sections[0].iter().zip(&cells[0]) {
            assert_eq!(section.start_m, initial.cells()[id].geometry.start_m);
            assert_eq!(section.end_m, initial.cells()[id].geometry.end_m);
        }
        // Equal section lengths give a diagnostic mean of rho0. The actual
        // path has heads (-10)*(-4.505) + (+10)*(-5.305) = -8 Pa.
        near(resolved.sources.branch_densities_kg_m3()[0], 1000.0);
        near(resolved.sources.section_anomalous_heads_pa()[0][0], 45.05);
        near(resolved.sources.section_anomalous_heads_pa()[0][1], -53.05);
        near(gravity.departure_heads_pa()[0], -8.0);
        near(gravity.absolute_departure_heads_pa()[0], 98.1);
        let solved = resolved
            .graph
            .solve(
                &[FixedPressure {
                    node: 0,
                    pressure: Pressure::new(0.0),
                }],
                GraphSolveConfig {
                    max_sweeps: 128,
                    max_node_iterations: 96,
                    absolute_flow_tolerance: VolumetricFlowRate::new(1e-15),
                    relative_flow_tolerance: 1e-13,
                },
                cx,
            )
            .unwrap();
        // Three equal native linear resistances close the original graph;
        // the sole fixed pressure is a datum, with no imposed pressure drop.
        for (branch, expected) in
            solved
                .branches
                .iter()
                .zip([-16.0 / 300000.0, 8.0 / 300000.0, -8.0 / 300000.0])
        {
            near(branch.flow.value(), expected);
        }
        let admitted = owner
            .admit_quasi_steady(cx, &initial, &resolved.sources, &solved, 0.1)
            .unwrap();
        assert!(admitted.neglected_buoyancy_power_w > 0.0);
        assert_eq!(initial.time_s(), 0.0);
    });
}

#[path = "source_graph/mixed.rs"]
mod mixed;
#[path = "source_graph/quasi_steady.rs"]
mod quasi_steady;

#[test]
fn actual_flow_film_source_and_native_link_account_for_wall_heat_once() {
    context(|cx| {
        let card = card(true, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = thermal(&caloric, [300.0; 6], Some(400.0), false);
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
        let resolved = owner.resolve(cx, &initial, &[1e-5, 0.0, 1e-5]).unwrap();
        let report = native_report(
            cx,
            &resolved.graph,
            &inventory,
            [1e-5, 0.0, 1e-5],
            [1.0, 0.0, 1.0],
            [1e-5, 0.0, 1e-5],
            0.1,
        );
        let heat = owner
            .evaluate_links(
                cx,
                &initial,
                &resolved.sources,
                &report,
                &[GraphConvectiveLink {
                    link: 9,
                    fluid_cell: CellId(2),
                    characteristic_length: Length::new(0.03),
                    property_temperature: ConvectionPropertyTemperature::FilmMean,
                    law: GraphConvectionLaw::Forced {
                        branch: 0,
                        correlation: CorrelationId::CircularDuctLaminarCwt,
                        settings: developed_settings(),
                    },
                }],
                16,
            )
            .unwrap();
        let sample = &heat.convection.samples()[0];
        near(sample.property_temperature_k, 350.0);
        near(sample.material.dynamic_viscosity_pa_s(), 0.0015);
        let conductance = 0.002 * 3.66 * 0.5 / 0.03;
        near(sample.conductance_w_k, conductance);
        near(sample.heat_into_fluid_w, 100.0 * conductance);
        let mut coefficients = StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        };
        heat.convection
            .apply_to(&mut coefficients.conductances_w_k, || false)
            .unwrap();
        let problem = BackwardEulerStep::new(&initial, 0.1, || Ok(())).unwrap();
        let trial = problem
            .evaluate(
                &[0.0; 7],
                |_| Ok::<_, DepositionError>(coefficients),
                || Ok(()),
            )
            .unwrap();
        near(trial.residual_j()[2], -0.1 * sample.heat_into_fluid_w);
        near(trial.residual_j()[6], 0.1 * sample.heat_into_fluid_w);
        near(trial.residual_j().iter().sum(), 0.0);
    });
}

// Regressions for interval support and native link identity. The original
// resistance/correlation owners decide admissibility; these tests supply only
// explicit native trajectories and the authored original validity intervals.
#[test]
fn lower_reynolds_support_checks_old_endpoint_and_signed_crossing() {
    context(|cx| {
        let card = card(false, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = thermal(&caloric, [300.0; 6], None, false);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let mut laws = laws(1e5);
        for section in &mut laws[0] {
            section.resistance = MaterialFlowResistance::try_new(
                A * 0.5,
                0.0,
                Area::new(PIPE_AREA),
                Length::new(0.03),
                (3000.0, 20000.0),
            )
            .unwrap();
        }
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
        // First: Re_old=1500, Re_avg=4500, Re_new=7500.
        // Second: both endpoints and the average are supported, but the native
        // signed trajectory crosses zero, outside the original Re>=3000 band.
        for (old, new, refused_re) in [(1e-4, 5e-4, 1500.0), (-6e-4, 1.2e-3, 0.0)] {
            let average = 0.5 * (old + new);
            let resolved = owner
                .resolve(cx, &initial, &[average, 0.0, average])
                .unwrap();
            let drive = 1e6 * (new - old) / 0.1 + 1e5 * average;
            let report = native_report(
                cx,
                &resolved.graph,
                &inventory,
                [old, 0.0, old],
                [drive, 0.0, drive],
                [average, 0.0, average],
                0.1,
            );
            near(report.average_volume_flows_m3_s[0], average);
            near(report.final_flows_m3_s[0], new);
            match owner.admit_actual(cx,&initial,&resolved.sources,&report) {
                Err(SourceGraphError::Resistance { branch:0,section:0,
                    source:fs_airflow::graph::thermal::medium::series::hydraulic::ResistanceError::Reynolds { value,admitted } }) => {
                    near(value,refused_re);
                    assert_eq!(admitted,(3000.0,20000.0));
                }
                other => panic!("original lower Reynolds support was not enforced: {other:?}"),
            }
        }
    });
}

#[test]
fn forced_card_checks_old_endpoint_and_zero_crossing_independently_of_drag() {
    use fs_airflow::graph::thermal::transport::transient::material::convection::GraphConvectionError;
    context(|cx| {
        let card = card(false, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = thermal(&caloric, [300.0; 6], Some(400.0), false);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(1e5); // Hydraulic support includes zero for this test.
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
        let links = [GraphConvectiveLink {
            link: 9,
            fluid_cell: CellId(2),
            characteristic_length: Length::new(0.03),
            property_temperature: ConvectionPropertyTemperature::Reservoir,
            law: GraphConvectionLaw::Forced {
                branch: 0,
                correlation: CorrelationId::Gnielinski,
                settings: developed_settings(),
            },
        }];
        for (old, new) in [(1e-4, 5e-4), (-6e-4, 1.2e-3)] {
            let average = 0.5 * (old + new);
            let resolved = owner
                .resolve(cx, &initial, &[average, 0.0, average])
                .unwrap();
            let drive = 1e6 * (new - old) / 0.1 + 1e5 * average;
            let report = native_report(
                cx,
                &resolved.graph,
                &inventory,
                [old, 0.0, old],
                [drive, 0.0, drive],
                [average, 0.0, average],
                0.1,
            );
            // The original drag law allows the full interval; the original
            // Gnielinski card must independently refuse Re_old=1500 or Re=0.
            owner
                .admit_actual(cx, &initial, &resolved.sources, &report)
                .unwrap();
            assert!(matches!(
                owner.evaluate_links(cx, &initial, &resolved.sources, &report, &links, 16),
                Err(GraphConvectionError::Envelope { link: 9, .. })
            ));
        }
    });
}

#[test]
fn reversed_link_ids_preserve_each_original_film_and_correlation() {
    context(|cx| {
        let card = card(true, 1.0);
        let (law, caloric) = source_pair(&card);
        let base = thermal(
            &caloric,
            [300.0, 300.0, 300.0, 300.0, 400.0, 300.0],
            Some(350.0),
            false,
        );
        let mut schedule: Vec<_> = base
            .cells()
            .values()
            .cloned()
            .map(|cell| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            })
            .collect();
        let mut second_wall = base.cells()[&CellId(6)].clone();
        second_wall.id = CellId(7);
        second_wall.segment_id = 7;
        second_wall.body_id = 2;
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(second_wall),
        });
        for (id, fluid, wall) in [(9, CellId(2), CellId(6)), (20, CellId(4), CellId(7))] {
            schedule.push(TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id,
                    a: fluid,
                    b: wall,
                    area_m2: 0.002,
                    conductance_w_k: 0.0,
                }),
            });
        }
        let initial = DepositionSimulation::new(
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
        .0;
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(1e5);
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
        let flows = [1e-4, 1.5e-4, 2.5e-4];
        let resolved = owner.resolve(cx, &initial, &flows).unwrap();
        let report = native_report(
            cx,
            &resolved.graph,
            &inventory,
            flows,
            [10.0, 7.5, 25.0],
            flows,
            0.1,
        );
        // Link 20 has 375 K film and turbulent Re=3600. Link 9 has 325 K
        // film and laminar Re~1714. Swapping their materials makes BOTH cards
        // inadmissible, so source identity matters beyond the heat debit.
        let high = GraphConvectiveLink {
            link: 20,
            fluid_cell: CellId(4),
            characteristic_length: Length::new(0.03),
            property_temperature: ConvectionPropertyTemperature::FilmMean,
            law: GraphConvectionLaw::Forced {
                branch: 1,
                correlation: CorrelationId::Gnielinski,
                settings: developed_settings(),
            },
        };
        let low = GraphConvectiveLink {
            link: 9,
            fluid_cell: CellId(2),
            characteristic_length: Length::new(0.03),
            property_temperature: ConvectionPropertyTemperature::FilmMean,
            law: GraphConvectionLaw::Forced {
                branch: 0,
                correlation: CorrelationId::CircularDuctLaminarCwt,
                settings: developed_settings(),
            },
        };
        let reversed = owner
            .evaluate_links(cx, &initial, &resolved.sources, &report, &[high, low], 16)
            .unwrap();
        let ordered = owner
            .evaluate_links(cx, &initial, &resolved.sources, &report, &[low, high], 16)
            .unwrap();
        let samples = reversed.convection.samples();
        assert_eq!(
            samples.iter().map(|sample| sample.link).collect::<Vec<_>>(),
            vec![9, 20]
        );
        near(samples[0].property_temperature_k, 325.0);
        near(samples[1].property_temperature_k, 375.0);
        near(samples[0].material.dynamic_viscosity_pa_s(), 0.00175);
        near(samples[1].material.dynamic_viscosity_pa_s(), 0.00125);
        assert!(samples[0].heat_into_fluid_w > 0.0);
        assert!(samples[1].heat_into_fluid_w < 0.0);
        for (a, b) in samples.iter().zip(ordered.convection.samples()) {
            assert_eq!(a.link, b.link);
            near(a.conductance_w_k, b.conductance_w_k);
            near(a.heat_into_fluid_w, b.heat_into_fluid_w);
        }
    });
}

#[test]
fn original_liquid_source_advances_distributed_momentum_and_reciprocal_wall_heat_together() {
    use fs_airflow::graph::thermal::medium::convection::deposition::source::SourceCellConvectionError;
    use fs_airflow::graph::thermal::transport::transient::{
        inertial::{
            InertialGraphHeat, InertialThermalControls, InertialThermalProblem,
            ResolvedInertialGraph, advance_inertial_graph_distributed,
        },
        material::convection::GraphConvectionError,
    };
    use fs_solver::NewtonKrylovConfig;
    context(|cx| {
        let card = card(true, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = thermal(&caloric, [300.0; 6], Some(400.0), false);
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
            .map(|branch| branch.unwrap())
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
        let links = [GraphConvectiveLink {
            link: 9,
            fluid_cell: CellId(2),
            characteristic_length: Length::new(0.03),
            property_temperature: ConvectionPropertyTemperature::FilmMean,
            law: GraphConvectionLaw::Forced {
                branch: 0,
                correlation: CorrelationId::CircularDuctLaminarCwt,
                settings: developed_settings(),
            },
        }];
        let (next_momentum, next, report) = advance_inertial_graph_distributed(
            cx,
            &problem,
            &inventory,
            &initial,
            0.1,
            controls,
            &[0.0; 7],
            &initial_flows,
            |cx, candidate, flows| {
                let resolved = owner
                    .resolve(cx, candidate, flows)
                    .map_err(GraphConvectionError::Source)?;
                Ok::<_, GraphConvectionError>(ResolvedInertialGraph {
                    graph: resolved.graph,
                    pressure_boundaries: boundaries.clone(),
                    driving_pressure_pa: vec![1.0, 0.0, 1.0],
                    payload: resolved.sources,
                })
            },
            |cx, candidate, sources, hydraulics| {
                let heat = owner.evaluate_links(cx, candidate, sources, hydraulics, &links, 16)?;
                let mut conductances = BTreeMap::new();
                heat.convection
                    .apply_to(&mut conductances, || cx.checkpoint().is_err())
                    .map_err(|source| {
                        GraphConvectionError::Convection(SourceCellConvectionError::Convection(
                            source,
                        ))
                    })?;
                Ok::<_, GraphConvectionError>(InertialGraphHeat {
                    mass_boundaries: vec![],
                    heat: StepCoefficients {
                        sources_w: BTreeMap::new(),
                        conductances_w_k: conductances,
                        payload: heat,
                    }
                    .into(),
                })
            },
        )
        .unwrap();
        near(next.time_s(), 0.1);
        near(next_momentum.accepted_time_s(), 0.1);
        assert!(next.temperature_k(CellId(2)).unwrap() > 300.0);
        assert!(next.temperature_k(CellId(6)).unwrap() < 400.0);
        assert!(report.heat.convection.samples()[0].heat_into_fluid_w > 0.0);
        assert!(report.maximum_cell_residual_j <= controls.cell_tolerance_j);
        assert!(report.hydraulics.maximum_impulse_residual_pa_s <= controls.impulse_tolerance_pa_s);
        assert!(
            report.hydraulics.mechanical_balance_residual_j.abs()
                <= controls.mechanical_tolerance_j
        );
        assert!(report.thermal.residual_j.abs() <= controls.total_thermal_tolerance_j);
        assert!(report.thermal.stored_energy_change_j.abs() < 1e-6);
        near(report.thermal.external_heat_j, 0.0);
        let total_change: f64 = next
            .cells()
            .iter()
            .map(|(id, cell)| {
                cell.mass_kg * (cell.specific_h_j_kg - initial.cells()[id].specific_h_j_kg)
            })
            .sum();
        assert!(total_change.abs() < 1e-6);
        let actual_t = next.temperature_k(CellId(2)).unwrap();
        let actual_source = &report.sources.sections()[0][0].momentum;
        near(
            actual_source.medium.dynamic_viscosity_pa_s(),
            0.005 - 1e-5 * actual_t,
        );
        assert!(actual_source.medium.dynamic_viscosity_pa_s() < 0.002);
        near(
            actual_source.evaluation.volume_flow_m3_s,
            report.hydraulics.average_volume_flows_m3_s[0],
        );
        // Proposals and acceptance return new owners; no old history changed.
        near(initial.time_s(), 0.0);
        near(momentum.accepted_time_s(), 0.0);
        near(initial.temperature_k(CellId(2)).unwrap(), 300.0);
        near(initial.temperature_k(CellId(6)).unwrap(), 400.0);
    });
}

#[path = "source_graph/axial_conduction.rs"]
mod axial_conduction;

#[path = "source_graph/dissipation.rs"]
mod dissipation;

#[path = "source_graph/reservoir.rs"]
mod reservoir;

#[path = "source_graph/developing.rs"]
mod developing;
