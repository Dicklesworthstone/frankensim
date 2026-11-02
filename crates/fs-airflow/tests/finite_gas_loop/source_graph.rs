//! The original NASA/source fixture now drives a branched finite graph.
use super::*;
use fs_airflow::graph::inertial::{InertialGraph, InertialGraphControls, InertialGraphReport};
use fs_airflow::graph::thermal::medium::convection::{
    NaturalConvectionGeometry, deposition::source::ConvectionPropertyTemperature,
};
use fs_airflow::graph::thermal::natural_circulation::material_response::MaterialBranchLaw;
use fs_airflow::graph::thermal::transport::transient::{
    distributed::DistributedGraphInventory,
    material::{
        GasGraphSources, SourceGraphBinding, SourceGraphError, SourceGraphLimits,
        convection::{GraphConvectionLaw, GraphConvectiveLink},
    },
};
use fs_airflow::graph::{FixedPressure, GraphBranch, LossGraph};
use fs_airflow::{LossElement, LossResistance, SourceProvenance, ToleranceBasis};
use fs_conduction::transient::deposition::ThermalLink;
use fs_material::fluid_transport::FluidTransport;
use fs_qty::Density;

fn graph() -> LossGraph {
    LossGraph::new(
        2,
        [(0, 1, "gas-a"), (0, 1, "gas-b"), (1, 0, "gas-return")]
            .into_iter()
            .map(|(from, to, name)| GraphBranch {
                from,
                to,
                loss: LossElement::new(
                    name,
                    LossResistance::new(1.0),
                    0.0,
                    SourceProvenance::new("synthetic fixed-composition gas graph", name),
                    ToleranceBasis::Analytic,
                )
                .unwrap(),
            })
            .collect(),
    )
    .unwrap()
}
fn paths() -> Vec<Vec<CellId>> {
    vec![vec![CellId(2), CellId(3)], vec![CellId(4)], vec![CellId(5)]]
}
fn laws() -> Vec<Vec<MaterialBranchLaw>> {
    [vec![0.5, 0.5], vec![1.0], vec![1.0]]
        .into_iter()
        .map(|weights| {
            weights
                .into_iter()
                .map(|weight| MaterialBranchLaw {
                    resistance: MaterialFlowResistance::try_new(
                        VISCOUS_GEOMETRY * weight,
                        0.0,
                        Area::new(AREA),
                        Length::new(0.03),
                        (0.0, 1e10),
                    )
                    .unwrap(),
                    source: SourceProvenance::new(
                        "authored gas pipe geometry",
                        "original-gas-sections",
                    ),
                    uncertainty_rel: 0.0,
                    tolerance_basis: ToleranceBasis::Analytic,
                    exchange: vec![],
                })
                .collect()
        })
        .collect()
}
fn source_limits() -> SourceGraphLimits {
    SourceGraphLimits {
        maximum_cells: 16,
        maximum_relative_density_departure: 0.2,
        maximum_neglected_section_power_w: 100.0,
        maximum_neglected_section_energy_j: 100.0,
    }
}
fn thermal(
    law: &GasEnthalpyLaw,
    temperatures: [f64; 6],
    wall: Option<f64>,
) -> DepositionSimulation<'_> {
    let rho0 = law.model().evaluate(law.context()).unwrap().state().density;
    let geometries = [
        ([0.0, 0.0, 0.0], [0.0, 0.1, 0.0], 0.1),
        ([1.0, 0.0, 0.0], [1.0, 0.1, 0.0], 0.1),
        ([0.0, 0.0, 0.0], [0.5, 0.0, 0.0], 0.5),
        ([0.5, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5),
        ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0),
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0),
    ];
    let mut schedule: Vec<_> = geometries
        .into_iter()
        .zip(temperatures)
        .enumerate()
        .map(|(i, ((start, end, length), t))| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: 0,
                material_id: MaterialId(1),
                mass_kg: rho0 * AREA * length,
                specific_h_j_kg: law.enthalpy_at_temperature(t).unwrap(),
                geometry: CellGeometry {
                    start_m: start,
                    end_m: end,
                    volume_m3: AREA * length,
                    surface_area_m2: 0.1,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
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
                mass_kg: rho0 * AREA,
                specific_h_j_kg: law.enthalpy_at_temperature(t).unwrap(),
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [0.1, 0.0, 0.0],
                    volume_m3: AREA,
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
                id: 11,
                a: CellId(2),
                b: CellId(6),
                area_m2: 0.002,
                conductance_w_k: 0.0,
            }),
        });
    }
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), law as &dyn ThermalLaw)],
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
    law: &GasEnthalpyLaw,
    thermal: &DepositionSimulation<'_>,
) -> DistributedGraphInventory {
    let rho = law.model().evaluate(law.context()).unwrap().state().density;
    DistributedGraphInventory::new(
        cx,
        graph,
        Density::new(rho),
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
) -> InertialGraphReport {
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
        inventory.reference_density(),
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
        (a - b).abs() <= 2e-9 * a.abs().max(b.abs()).max(1e-8),
        "{a} != {b}"
    );
}

#[test]
fn fixed_pressure_composition_nasa_enthalpy_and_local_eos_survive_graph_binding() {
    context(|cx| {
        let law = gas_law();
        let initial = thermal(&law, [300.0, 300.0, 280.0, 340.0, 320.0, 300.0], None);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &law, &initial);
        let laws = laws();
        let owner = GasGraphSources::try_new(
            cx,
            &law,
            &initial,
            SourceGraphBinding {
                graph: &graph,
                inventory: &inventory,
                laws: &laws,
                acceleration_m_s2: [-9.81, 0.0, 0.0],
            },
            source_limits(),
            0.1,
        )
        .unwrap();
        let result = owner.resolve(cx, &initial, &[0.0; 3]).unwrap();
        let rho0 = owner.reference().state().density;
        for (branch, path) in paths().iter().enumerate() {
            for (section, &id) in path.iter().enumerate() {
                let actual = &result.sources.sections()[branch][section].momentum.medium;
                let t = initial.temperature_k(id).unwrap();
                let direct = law.model().evaluate(&point(t)).unwrap();
                assert_eq!(actual, &direct);
                near(actual.state().pressure, PRESSURE);
                near(actual.state().density, PRESSURE * MOLAR_MASS / (R * t));
                near(actual.specific_enthalpy_j_kg(), h(t));
                near(
                    initial.cells()[&id].mass_kg,
                    rho0 * initial.cells()[&id].geometry.volume_m3,
                );
            }
        }
        assert_eq!(owner.source().mole_composition().fractions(), &[0.75, 0.25]);
        let expected = rho0 * 0.5 * (300.0 / 280.0 + 300.0 / 340.0);
        near(result.sources.branch_densities_kg_m3()[0], expected);
        assert!((expected - rho0 * 300.0 / 310.0).abs() > 1e-3);
        let heads = result.sources.section_anomalous_heads_pa()[0]
            .iter()
            .sum::<f64>();
        near(heads, (expected - rho0) * (-9.81));
    });
}

#[test]
fn gas_endpoint_mach_cannot_hide_behind_zero_average_flow() {
    context(|cx| {
        let law = gas_law();
        let initial = thermal(&law, [300.0; 6], None);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &law, &initial);
        let laws = laws();
        let owner = GasGraphSources::try_new(
            cx,
            &law,
            &initial,
            SourceGraphBinding {
                graph: &graph,
                inventory: &inventory,
                laws: &laws,
                acceleration_m_s2: [0.0; 3],
            },
            source_limits(),
            1e-4,
        )
        .unwrap();
        let result = owner.resolve(cx, &initial, &[0.0; 3]).unwrap();
        let inertia = inventory.reference_density().value() / AREA;
        let drive = -2.0 * inertia * 1e-4 / 0.1;
        let report = native_report(
            cx,
            &result.graph,
            &inventory,
            [1e-4, 0.0, 1e-4],
            [drive, 0.0, drive],
            [0.0; 3],
            0.1,
        );
        assert!(report.average_volume_flows_m3_s[0].abs() < 1e-12);
        assert!(matches!(
            owner.admit_actual(cx, &initial, &result.sources, &report),
            Err(SourceGraphError::GasMach {
                branch: 0,
                section: 0,
                ..
            })
        ));
        assert!(matches!(
            owner.resolve(cx, &initial, &[1e-4, 0.0, 1e-4]),
            Err(SourceGraphError::GasMach { .. })
        ));
    });
}

#[test]
fn exact_registered_gas_law_is_required_even_for_identical_source_values() {
    context(|cx| {
        let law = gas_law();
        let other = gas_law();
        let initial = thermal(&law, [300.0; 6], None);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &law, &initial);
        let laws = laws();
        assert!(matches!(
            GasGraphSources::try_new(
                cx,
                &other,
                &initial,
                SourceGraphBinding {
                    graph: &graph,
                    inventory: &inventory,
                    laws: &laws,
                    acceleration_m_s2: [0.0; 3]
                },
                source_limits(),
                0.1
            ),
            Err(SourceGraphError::Cell { .. })
        ));
    });
}

#[test]
fn stagnant_gas_has_sourced_natural_wall_exchange_without_a_forced_speed_floor() {
    context(|cx| {
        let law = gas_law();
        let initial = thermal(&law, [300.0; 6], Some(330.0));
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &law, &initial);
        let laws = laws();
        let owner = GasGraphSources::try_new(
            cx,
            &law,
            &initial,
            SourceGraphBinding {
                graph: &graph,
                inventory: &inventory,
                laws: &laws,
                acceleration_m_s2: [0.0, 0.0, -9.81],
            },
            source_limits(),
            0.1,
        )
        .unwrap();
        let result = owner.resolve(cx, &initial, &[0.0; 3]).unwrap();
        let report = native_report(
            cx,
            &result.graph,
            &inventory,
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
            0.1,
        );
        let heat = owner
            .evaluate_links(
                cx,
                &initial,
                &result.sources,
                &report,
                &[GraphConvectiveLink {
                    link: 11,
                    fluid_cell: CellId(2),
                    characteristic_length: Length::new(0.03),
                    property_temperature: ConvectionPropertyTemperature::FilmMean,
                    law: GraphConvectionLaw::Natural {
                        geometry: NaturalConvectionGeometry::VerticalPlate,
                    },
                }],
                16,
            )
            .unwrap();
        let sample = &heat.convection.samples()[0];
        near(sample.property_temperature_k, 315.0);
        near(sample.material.state().pressure, PRESSURE);
        let transport = sample.material.transport_state().unwrap();
        let ra = 9.81 / 315.0 * 30.0 * 0.03_f64.powi(3)
            / (transport.kinematic_viscosity_m2_s() * transport.thermal_diffusivity_m2_s());
        let pr = transport.prandtl();
        let nu = (0.825
            + 0.387 * ra.powf(1.0 / 6.0) / (1.0 + (0.492 / pr).powf(9.0 / 16.0)).powf(8.0 / 27.0))
        .powi(2);
        let conductance = 0.002 * nu * transport.thermal_conductivity_w_m_k() / 0.03;
        near(sample.conductance_w_k, conductance);
        near(sample.heat_into_fluid_w, 30.0 * conductance);
        assert!(sample.heat_into_fluid_w > 0.0);
        assert!(
            heat.admission
                .actual_losses
                .iter()
                .flatten()
                .all(|l| l.dissipated_power_w.abs() < 1e-20)
        );
    });
}
#[path = "source_graph/quasi_steady.rs"]
mod quasi_steady;

#[path = "source_graph/axial_conduction.rs"]
mod axial_conduction;

