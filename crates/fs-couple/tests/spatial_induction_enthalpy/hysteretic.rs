//! Native three-dimensional memory and physical loss reach finite thermal cells.
#[path = "../../../fs-em/tests/spatial_induction/hysteretic/common.rs"]
mod data;

use super::{near, regions, thermal_controls};
use data::{
    CP, DIRECTIONS, DT, MU, PIN, bracket, conductor, endpoint, field_controls, point, resistivity,
    stored,
};
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError,
    magnetic::MagneticThermalError,
    material::MaterialCell,
    spatial_induction::{SpatialMagneticRegion, SpatialMagneticThermalNetwork},
};
use fs_em::spatial::vector_potential::induction::hysteretic::SpatialHystereticProblem;
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint, SelectionPolicy};
use fs_material::{
    conductor::ElectricalPropertyBasis,
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::SPECIFIC_HEAT_CAPACITY_DIMS,
    thermal::source::ThermalMaterialPoint,
};
use fs_qty::Density;
use fs_thermal::enthalpy::*;
use std::sync::Arc;

fn cube(card: &Arc<MaterialCard>) -> SpatialHystereticProblem {
    // Reuse the original spatial induction mesh rather than introducing a
    // second geometric convention for the independent one-mode balance.
    let basis = ElectricalPropertyBasis::Resistivity;
    let geometry = super::problem(
        &super::source::absolute_affine_magnetic_card(0.2, 0.0, 1.0, 0.0, basis),
        basis,
    );
    data::problem(
        card,
        geometry.positions_m().to_vec(),
        geometry.mesh().tets.clone(),
    )
}
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(500.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-9,
    )
    .unwrap()
}
fn thermal(curve: &ResolvedHeatCapacityCurve, volume: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(7),
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(9),
                mass_kg: volume,
                specific_h_j_kg: CP * (300.0 - 200.0),
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: volume,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        }],
        &[(MaterialId(9), curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 1,
            max_events: 1,
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
) -> [MaterialCell<'a>; 1] {
    [MaterialCell {
        cell: CellId(7),
        material_id: MaterialId(9),
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: card,
            context: query,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }]
}

#[test]
fn spatial_play_and_temperature_feedback_share_one_finite_energy_balance_and_history() {
    let card = data::card(500.0);
    let curve = caloric(&card);
    let problem = cube(&card);
    let query = point(300.0);
    let cells = cells(&card, &curve, &query);
    let regions = regions();
    let network = SpatialMagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let mut heat = thermal(&curve, 1.0);
    let mut state = data::initial(&problem, Some(conductor(&card, 300.0)));
    let (mut old_h, mut old_z, mut old_b, mut old_t) = (0.0, 0.0, 0.0, 300.0);
    for (index, current) in [63.12, 0.0, 0.0].into_iter().enumerate() {
        // Nested scalar bisection uses exact cube integrals and elementary
        // scalar play, never production matrices or the vector return map.
        let physical_heat = |t: f64| {
            let (_, z, b) = endpoint(old_b, old_z, resistivity(t).recip(), DT, current);
            3.0 / (160.0 * resistivity(t) * DT) * (b - old_b).powi(2) + PIN * (z - old_z).abs()
        };
        let t = bracket(old_t, 480.0, |t| CP * (t - old_t) - physical_heat(t));
        let (h, z, b) = endpoint(old_b, old_z, resistivity(t).recip(), DT, current);
        let joule = 3.0 / (160.0 * resistivity(t) * DT) * (b - old_b).powi(2);
        let pinning = PIN * (z - old_z).abs();
        let change = stored(h, z) - stored(old_h, old_z);
        let numerical = h * (b - old_b) - change - pinning;
        let target = (index + 1) as f64 * DT;
        let before = (
            heat.cells().clone(),
            state.potential_wb().to_vec(),
            state.cell_states().to_vec(),
        );
        let mut checkpoints = 0;
        let candidate = network
            .step_hysteretic_candidate(
                &heat,
                target,
                &state,
                &[[current; 3]; 6],
                thermal_controls(),
                field_controls(),
                &[0.0],
                || {
                    checkpoints += 1;
                    false
                },
            )
            .unwrap();
        near(
            candidate.thermal().temperature_k(CellId(7)).unwrap(),
            t,
            3e-8,
        );
        let edge = problem.mesh().edges.binary_search(&[0, 7]).unwrap();
        near(candidate.field().potential_wb()[edge], 0.5 * b, 3e-9);
        for (memory, direction) in candidate.field().cell_states().iter().zip(DIRECTIONS) {
            for axis in 0..3 {
                near(memory.field_a_m()[axis], h * direction[axis], 3e-8);
                near(memory.flux_density_t()[axis], b * direction[axis], 3e-9);
                near(memory.polarizations_t()[0][axis], z * direction[axis], 3e-9);
            }
        }
        let report = candidate.report();
        let native = &report.coupled.electrical;
        near(native.joule_heat_j, joule, 3e-8);
        near(native.hysteresis_heat_j, pinning, 3e-8);
        near(native.numerical_dissipation_j, numerical, 3e-8);
        near(native.stored_change_j, change, 3e-8);
        near(candidate.deposited_heat_j(), joule + pinning, 3e-8);
        near(
            report.coupled.thermal.stored_energy_change_j,
            joule + pinning,
            3e-8,
        );
        near(report.total_energy_residual_j, 0.0, 3e-8);
        near(report.external_heat_j, 0.0, 3e-8);
        for source in candidate.field().conductors().iter().flatten() {
            near(source.resistivity_ohm_m(), resistivity(t), 1e-12);
            near(source.material().query().axes()["T"], t, 3e-8);
        }
        if index == 0 {
            assert!(joule > 1.0 && pinning > 1.0 && numerical > 1.0);
            let frozen = endpoint(old_b, old_z, resistivity(old_t).recip(), DT, current);
            assert!((h - frozen.0).abs() > 0.1);
            assert!(
                (candidate.thermal().energy_j() - heat.energy_j() - joule - pinning - numerical)
                    .abs()
                    > 1.0
            );
            assert!(checkpoints > 10);
            let mut reached = 0;
            assert!(matches!(
                network.step_hysteretic_candidate(
                    &heat,
                    target,
                    &state,
                    &[[current; 3]; 6],
                    thermal_controls(),
                    field_controls(),
                    &[0.0],
                    || {
                        reached += 1;
                        reached == checkpoints - 1
                    },
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
        }
        assert_eq!(heat.cells(), &before.0);
        assert_eq!(state.potential_wb(), before.1);
        assert_eq!(state.cell_states(), before.2);
        assert_eq!(state.time_s(), heat.time_s());
        let (next_heat, next_state, _) = candidate.into_parts();
        heat = next_heat;
        state = next_state;
        (old_h, old_z, old_b, old_t) = (h, z, b, t);
    }
    assert!(old_z > 0.4 && old_z < 1.0);
}

#[test]
fn insulating_spatial_magnets_receive_pinning_heat_and_refuse_missing_or_invalid_receivers() {
    let card = data::card(500.0);
    let curve = caloric(&card);
    let problem = cube(&card);
    let state = data::initial(&problem, None);
    let heat = thermal(&curve, 1.0);
    let query = point(300.0);
    let cells = cells(&card, &curve, &query);
    let mut regions = regions();
    let network = SpatialMagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    // With no conductors, exact H=J_component/4=12 and z=1 throughout the cube.
    let candidate = network
        .step_hysteretic_candidate(
            &heat,
            DT,
            &state,
            &[[48.0; 3]; 6],
            thermal_controls(),
            field_controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    let native = &candidate.report().coupled.electrical;
    assert!(native.cells.is_empty());
    assert_eq!(native.magnetic_cells.len(), 6);
    near(native.joule_heat_j, 0.0, 1e-12);
    near(native.hysteresis_heat_j, PIN, 3e-8);
    near(candidate.deposited_heat_j(), PIN, 3e-8);
    near(
        candidate.report().coupled.thermal.stored_energy_change_j,
        PIN,
        3e-8,
    );
    near(
        candidate.thermal().temperature_k(CellId(7)).unwrap(),
        300.0 + PIN / CP,
        3e-8,
    );
    near(native.stored_after_j, stored(12.0, 1.0), 3e-8);
    near(
        native.numerical_dissipation_j,
        12.0 * (12.0 * MU + 1.0) - stored(12.0, 1.0) - PIN,
        3e-8,
    );
    near(candidate.report().total_energy_residual_j, 0.0, 3e-8);
    for cell in &native.magnetic_cells {
        near(cell.irreversible_heat_j, PIN / 6.0, 3e-8);
    }
    regions[0] = SpatialMagneticRegion::Fixed { tetrahedron: 0 };
    // Match the remaining mapped volume: rejection must reach the insulating
    // cell's actual pinning heat, rather than a volume-accounting mismatch.
    let incomplete_heat = thermal(&curve, 5.0 / 6.0);
    let missing = SpatialMagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    assert!(matches!(
        missing.step_hysteretic_candidate(
            &incomplete_heat,
            DT,
            &state,
            &[[48.0; 3]; 6],
            thermal_controls(),
            field_controls(),
            &[0.0],
            || false,
        ),
        Err(CircuitEnthalpyError::Source(MagneticThermalError::Invalid(
            "every dissipating magnetic tetrahedron needs an actual thermal heat destination"
        )))
    ));
    assert_eq!(state.time_s(), 0.0);
    assert!(
        state
            .cell_states()
            .iter()
            .all(|s| s.polarizations_t() == [[0.0; 3]])
    );
    near(
        incomplete_heat.temperature_k(CellId(7)).unwrap(),
        300.0,
        1e-12,
    );

    let unsupported = data::card(310.0);
    let unsupported_problem = cube(&unsupported);
    let unsupported_state = data::initial(&unsupported_problem, None);
    let unsupported_curve = caloric(&unsupported);
    let unsupported_heat = thermal(&unsupported_curve, 1.0);
    let unsupported_cells = self::cells(&unsupported, &unsupported_curve, &query);
    let complete = super::regions();
    let invalid = SpatialMagneticThermalNetwork {
        cells: &unsupported_cells,
        links: &[],
        magnetic_regions: &complete,
        volume_tolerance_m3: 1e-12,
    };
    assert!(matches!(
        invalid.step_hysteretic_candidate(
            &unsupported_heat,
            DT,
            &unsupported_state,
            &[[48.0; 3]; 6],
            thermal_controls(),
            field_controls(),
            &[0.0],
            || false,
        ),
        Err(CircuitEnthalpyError::Source(
            MagneticThermalError::Magnetic(_)
        ))
    ));
    assert_eq!(unsupported_state.time_s(), 0.0);
    near(
        unsupported_heat.temperature_k(CellId(7)).unwrap(),
        300.0,
        1e-12,
    );
}
