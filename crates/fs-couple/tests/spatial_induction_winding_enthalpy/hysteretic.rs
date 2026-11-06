//! Native spatial play, eddy currents, strands and graph bodies share finite heat.
#[allow(dead_code)]
#[path = "../../../fs-em/tests/spatial_induction_winding/hysteretic_common.rs"]
mod data;

use data::*;
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, ResistorCell,
    material::{CellResistor, MaterialCell},
    spatial_induction::{
        SpatialMagneticRegion, SpatialMagneticThermalNetwork,
        hysteretic::winding::SpatialWindingThermalNetwork,
    },
};
use fs_em::winding::nonlinear::resistors::GraphResistorMaterial;
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint, SelectionPolicy};
use fs_material::{
    conductor::ElectricalPropertyBasis,
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS},
    thermal::source::ThermalMaterialPoint,
};
use fs_qty::Density;
use fs_thermal::enthalpy::*;

const CELL: CellId = CellId(7);
const GRAPH_CELL: CellId = CellId(8);
const MAT: MaterialId = MaterialId(9);
const BASIS: ElectricalPropertyBasis = ElectricalPropertyBasis::Resistivity;

fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &material::point(200.0),
            upper: &material::point(500.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-9,
    )
    .unwrap()
}
fn thermal(curve: &ResolvedHeatCapacityCurve, separate_graph: bool) -> DepositionSimulation<'_> {
    let ids = if separate_graph {
        vec![CELL, GRAPH_CELL]
    } else {
        vec![CELL]
    };
    DepositionSimulation::new(
        ids.into_iter()
            .map(|id| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: 0,
                    body_id: id.0,
                    material_id: MAT,
                    mass_kg: 1.0,
                    specific_h_j_kg: CAPACITY * (300.0 - 200.0),
                    geometry: CellGeometry {
                        start_m: [0.0; 3],
                        end_m: [1.0, 0.0, 0.0],
                        volume_m3: 1.0,
                        surface_area_m2: 1.0,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            })
            .collect(),
        &[(MAT, curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 2,
            max_substeps: 1000,
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
fn cells<'a>(
    card: &'a MaterialCard,
    curve: &'a ResolvedHeatCapacityCurve,
    query: &'a QueryPoint,
    separate_graph: bool,
) -> Vec<MaterialCell<'a>> {
    let ids = if separate_graph {
        vec![CELL, GRAPH_CELL]
    } else {
        vec![CELL]
    };
    ids.into_iter()
        .map(|cell| MaterialCell {
            cell,
            material_id: MAT,
            curve: curve.curve(),
            point: ThermalMaterialPoint {
                material: card,
                context: query,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        })
        .collect()
}
fn regions() -> Vec<SpatialMagneticRegion> {
    (0..48)
        .map(|tetrahedron| SpatialMagneticRegion::Cell {
            tetrahedron,
            cell: CELL,
        })
        .collect()
}

#[test]
fn finite_temperature_feedback_matches_independent_voltage_play_and_four_physical_losses() {
    let card = material::card_with_heat_capacity(500.0, CAPACITY);
    let curve = caloric(&card);
    let problem = cube(&card);
    let query = material::point(300.0);
    let cells = cells(&card, &curve, &query, false);
    let regions = regions();
    let routes = [ResistorCell {
        branch: 1,
        cell: CELL,
    }];
    let network = SpatialWindingThermalNetwork {
        magnetic: SpatialMagneticThermalNetwork {
            cells: &cells,
            links: &[],
            magnetic_regions: &regions,
            volume_tolerance_m3: 1e-12,
        },
        graph_resistors: &routes,
        fill_tolerance: 1e-12,
    };
    let mut heat = thermal(&curve, false);
    let mut accepted = circuit(&problem, &material::conductor(&card, 300.0));
    let mut old = Mode::default();
    let mut old_t = 300.0;
    for (index, input) in [180.0, 0.0].into_iter().enumerate() {
        // Scalar temperature bisection and exact mode integrals are independent
        // of both the native vector return map and the thermal Newton solve.
        let physical = |t| {
            let next = voltage_endpoint(old, t, input, GRAPH_R);
            core_heat(old, next, t)
                + pinning_heat(old, next)
                + DT * (8.0 * material::resistivity(t) + GRAPH_R) * next.i.powi(2)
        };
        let t = material::bracket(old_t, 480.0, |t| CAPACITY * (t - old_t) - physical(t));
        let expected = voltage_endpoint(old, t, input, GRAPH_R);
        let before = (
            heat.cells().clone(),
            accepted.field().cell_states().to_vec(),
        );
        let candidate = network
            .step_hysteretic_candidate(
                &heat,
                (index + 1) as f64 * DT,
                &accepted,
                &[input],
                super::thermal_controls(),
                controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        let report = candidate.report();
        let step = &report.coupled.electrical;
        near(candidate.thermal().temperature_k(CELL).unwrap(), t, 3e-7);
        near(step.currents_a[0], expected.i, 4e-7);
        assert_field(candidate.circuit().field(), expected);
        near(step.field.joule_heat_j, core_heat(old, expected, t), 3e-8);
        near(
            step.field.hysteresis_heat_j,
            pinning_heat(old, expected),
            3e-8,
        );
        near(
            step.cells.iter().map(|cell| cell.heat_j).sum(),
            DT * 8.0 * material::resistivity(t) * expected.i.powi(2),
            4e-7,
        );
        near(
            step.graph_resistors[0].heat_j,
            DT * GRAPH_R * expected.i.powi(2),
            4e-7,
        );
        near(step.heat_j, physical(t), 5e-7);
        near(candidate.deposited_heat_j(), physical(t), 5e-7);
        near(
            report.coupled.thermal.stored_energy_change_j,
            physical(t),
            5e-7,
        );
        near(step.supplied_work_j, DT * input * expected.i, 7e-7);
        near(report.external_heat_j, 0.0, 1e-12);
        near(
            report.total_energy_residual_j,
            0.0,
            super::thermal_controls().total_energy_tolerance_j,
        );
        assert!(step.numerical_dissipation_j > 1.0);
        assert_eq!(heat.cells(), &before.0);
        assert_eq!(accepted.field().cell_states(), before.1);
        for source in candidate.circuit().conductors() {
            near(
                source.conductor.resistivity_ohm_m(),
                material::resistivity(t),
                1e-12,
            );
            near(source.conductor.material().query().axes()["T"], t, 3e-7);
        }
        for source in candidate.circuit().field().conductors().iter().flatten() {
            near(source.resistivity_ohm_m(), material::resistivity(t), 1e-12);
        }
        if index == 0 {
            assert!(t > 340.0 && t < 430.0);
            assert!(expected.i < voltage_endpoint(old, 300.0, input, GRAPH_R).i);
        } else {
            assert!(step.stored_change_j < 0.0 && expected.i > 0.0);
            assert!(expected.z > material::play(0.0, expected.h) + 0.1);
        }
        old = expected;
        old_t = t;
        let (next_heat, next_circuit, _) = candidate.into_parts();
        heat = next_heat;
        accepted = next_circuit;
    }
}

#[test]
fn material_graph_heat_has_its_own_receiver_and_cancelled_or_mismatched_sources_cannot_publish() {
    let card = material::card_with_heat_capacity(500.0, CAPACITY);
    let source = material::conductor(&card, 300.0);
    let curve = caloric(&card);
    let problem = cube(&card);
    let query = material::point(300.0);
    let cells = cells(&card, &curve, &query, true);
    let regions = regions();
    let routes = [ResistorCell {
        branch: 1,
        cell: GRAPH_CELL,
    }];
    let length = (GRAPH_R / material::resistivity(300.0)).sqrt();
    let area = 1.0 / length;
    let resistors = [CellResistor {
        destination: routes[0],
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: length,
        area_m2: area,
    }];
    let bases = [(1, BASIS)];
    let network = SpatialWindingThermalNetwork {
        magnetic: SpatialMagneticThermalNetwork {
            cells: &cells,
            links: &[],
            magnetic_regions: &regions,
            volume_tolerance_m3: 1e-12,
        },
        graph_resistors: &routes,
        fill_tolerance: 1e-12,
    };
    let heat = thermal(&curve, true);
    let original = circuit(&problem, &source)
        .with_graph_resistors(
            vec![GraphResistorMaterial {
                branch: 1,
                conductor: source.into_uniform_conductor(length, area).unwrap(),
            }],
            &mut || false,
        )
        .unwrap();
    let before = (
        heat.cells().clone(),
        original.state().to_vec(),
        original.field().cell_states().to_vec(),
    );
    let mut polls = 0;
    let candidate = network
        .step_hysteretic_material_candidate(
            &heat,
            DT,
            &original,
            &[120.0],
            &resistors,
            &bases,
            super::thermal_controls(),
            controls(),
            &[0.0; 2],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let next = candidate.circuit();
    let step = &candidate.report().coupled.electrical;
    let t = candidate.thermal().temperature_k(CELL).unwrap();
    let tr = candidate.thermal().temperature_k(GRAPH_CELL).unwrap();
    let rho = material::resistivity(t);
    let r = material::resistivity(tr) * length / area;
    let expected = voltage_endpoint(Mode::default(), t, 120.0, r);
    assert_field(next.field(), expected);
    near(next.currents_a()[0], expected.i, 4e-7);
    near(
        step.cells.iter().map(|cell| cell.heat_j).sum(),
        DT * 8.0 * rho * expected.i.powi(2),
        4e-7,
    );
    near(
        CAPACITY * (t - 300.0),
        step.field.heat_j + step.cells.iter().map(|cell| cell.heat_j).sum::<f64>(),
        5e-7,
    );
    near(CAPACITY * (tr - 300.0), DT * r * expected.i.powi(2), 5e-7);
    near(
        step.graph_resistors[0].heat_j,
        DT * r * expected.i.powi(2),
        4e-7,
    );
    assert!(tr > t && t > 300.0);
    assert_eq!(step.graph_materials.len(), 1);
    near(step.graph_materials[0].conductor.resistance_ohm(), r, 1e-12);
    near(
        step.graph_materials[0].conductor.material().query().axes()["T"],
        tr,
        3e-7,
    );
    near(
        candidate.deposited_heat_j(),
        CAPACITY * (t + tr - 600.0),
        5e-7,
    );
    near(candidate.report().total_energy_residual_j, 0.0, 5e-7);
    assert!(polls > 20);
    let mut seen = 0;
    assert!(matches!(
        network.step_hysteretic_material_candidate(
            &heat,
            DT,
            &original,
            &[120.0],
            &resistors,
            &bases,
            super::thermal_controls(),
            controls(),
            &[0.0; 2],
            || {
                seen += 1;
                seen >= polls - 1
            },
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let wrong = original
        .with_graph_resistors(
            vec![GraphResistorMaterial {
                branch: 1,
                conductor: material::conductor(&card, 301.0)
                    .into_uniform_conductor(length, area)
                    .unwrap(),
            }],
            &mut || false,
        )
        .unwrap();
    assert!(
        network
            .step_hysteretic_material_candidate(
                &heat,
                DT,
                &wrong,
                &[120.0],
                &resistors,
                &bases,
                super::thermal_controls(),
                controls(),
                &[0.0; 2],
                || false,
            )
            .is_err()
    );
    assert_eq!(heat.cells(), &before.0);
    assert_eq!(original.state(), before.1);
    assert_eq!(original.field().cell_states(), before.2);
    assert_eq!(original.time_s(), 0.0);
}
