//! Native heat, original caloric volume and actual outlet flow share one state.

#[path = "thermal_dilation/isobaric.rs"]
mod isobaric;

use super::{DEPTH, binding, initial_thermal, near, thermal_config, zero};
use std::collections::BTreeMap;

use fs_blake3::ContentHash;
use fs_couple::stokes_motion::thermal_dilation::{
    StokesThermalDilation, StokesThermalDilationBinding, StokesThermalDilationControls,
    StokesThermalDilationError, solve_at_thermal_state,
};
use fs_flux::ns::boundary::thermal_dilation::CaloricDilationError;
use fs_flux::{FluxParams, FluxSystem, StokesBoundary, StokesBoundaryError, TriMesh};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_solid::Mesh2;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionSimulation, MaterialId, ThermalLaw, ThermalStorageBasis,
};

// Explicit manufactured thermodynamic charts, not named-material predictions.
fn curve(final_density: f64, phase: SolidLiquidPhase) -> EquilibriumEnthalpyPhaseCurve {
    let liquid_mass_fraction = if phase == SolidLiquidPhase::Liquid {
        1.0
    } else {
        0.0
    };
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0xa3; 32]),
        phase,
        vec![
            EnthalpyPhaseKnot {
                temperature_k: 300.0,
                specific_enthalpy_j_kg: 0.0,
                liquid_mass_fraction,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                temperature_k: 400.0,
                specific_enthalpy_j_kg: 100.0,
                liquid_mass_fraction,
                bulk_density_kg_m3: final_density,
            },
        ],
    )
    .unwrap()
}

fn zero_traction(_: [f64; 2], _: [f64; 2]) -> [f64; 2] {
    [0.0; 2]
}

fn open_boundaries(mesh: &TriMesh) -> Vec<StokesBoundary<'static>> {
    mesh.edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, data)| {
            if data.mid[0] == 1.0 {
                StokesBoundary::TotalTraction {
                    edge,
                    value: &zero_traction,
                }
            } else if data.mid[0] == 0.0 {
                StokesBoundary::Velocity { edge, value: &zero }
            } else {
                StokesBoundary::NavierSlip {
                    edge,
                    value: &zero,
                    tangential_resistance_m_s: 0.0,
                }
            }
        })
        .collect()
}

fn solve(
    system: &FluxSystem<'_>,
    thermal: &DepositionSimulation<'_>,
    ids: &[CellId],
    caloric: &EquilibriumEnthalpyPhaseCurve,
    boundaries: &[StokesBoundary<'_>],
) -> Result<StokesThermalDilation, StokesThermalDilationError> {
    let charts = BTreeMap::from([(MaterialId(9), caloric)]);
    solve_at_thermal_state(
        system,
        thermal,
        StokesThermalDilationBinding {
            cells: binding(ids),
            depth_m: DEPTH,
            caloric_by_material: &charts,
        },
        FluxParams::default(),
        &|_, _| 2.0,
        &zero,
        boundaries,
        StokesThermalDilationControls {
            maximum_dilation_defect_per_s: 2e-11,
            maximum_volume_rate_defect_m3_s: 2e-11,
        },
        || false,
    )
}

#[test]
fn native_heat_drives_expansion_contraction_and_zero_flow_without_advancing() {
    for n in [1, 2] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, n, n));
        let system = FluxSystem::new(&mesh);
        let boundaries = open_boundaries(&mesh);
        for (final_density, h, rate, expected_dilation) in [
            (800.0, 0.0, 40.0, 0.1),
            (800.0, 50.0, -45.0, -0.1),
            (4000.0 / 3.0, 0.0, 40.0, -0.1),
            (800.0, 50.0, 0.0, 0.0),
        ] {
            let caloric = curve(final_density, SolidLiquidPhase::Liquid);
            let rho = caloric
                .state_at_specific_enthalpy(h)
                .unwrap()
                .bulk_density_kg_m3();
            let (thermal, ids) = initial_thermal(&mesh, &caloric, false, |_, cell| {
                cell.specific_h_j_kg = h;
                cell.mass_kg = rho * cell.geometry.volume_m3;
                cell.boundary.source_w = cell.mass_kg * rate;
            });
            let before = thermal.cells().clone();
            let answer = solve(&system, &thermal, &ids, &caloric, &boundaries).unwrap();
            near(
                answer.expected_volume_rate_m3_s,
                DEPTH * expected_dilation,
                1e-13,
            );
            near(
                answer.boundary_volume_rate_m3_s,
                DEPTH * expected_dilation,
                2e-11,
            );
            assert!(answer.volume_rate_defect_m3_s < 2e-11);
            for point in [[0.2, 0.3], [0.7, 0.4]] {
                let velocity = system.velocity_at(&answer.flow.flow.field.x, point);
                near(velocity[0], expected_dilation * point[0], 2e-11);
                near(velocity[1], 0.0, 2e-11);
            }
            for row in answer.cells {
                assert_eq!(row.caloric_identity, *caloric.identity().as_bytes());
                near(row.heat_power_w, before[&row.cell].boundary.source_w, 1e-12);
                near(row.material_specific_enthalpy_rate_w_kg, rate, 1e-12);
                near(
                    row.volume_rate_m3_s,
                    before[&row.cell].geometry.volume_m3 * expected_dilation,
                    1e-13,
                );
            }
            assert_eq!(answer.time_s, 0.0);
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(*thermal.cells(), before);
        }
    }
}

#[test]
fn conservative_native_conduction_can_change_volume_and_preserves_the_enthalpy_datum() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let boundaries = open_boundaries(&mesh);
    let mut first_field = None;
    for datum in [0.0, 1e6] {
        let caloric = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0xa4; 32]),
            SolidLiquidPhase::Liquid,
            [
                (300.0, 0.0, 1000.0),
                (350.0, 50.0, 800.0),
                (400.0, 100.0, 500.0),
            ]
            .into_iter()
            .map(|(t, h, rho)| EnthalpyPhaseKnot {
                temperature_k: t,
                specific_enthalpy_j_kg: datum + h,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: rho,
            })
            .collect(),
        )
        .unwrap();
        let (thermal, ids) = initial_thermal(&mesh, &caloric, true, |index, cell| {
            let h = datum + if index == 0 { 20.0 } else { 80.0 };
            cell.specific_h_j_kg = h;
            cell.mass_kg = caloric
                .state_at_specific_enthalpy(h)
                .unwrap()
                .bulk_density_kg_m3()
                * cell.geometry.volume_m3;
        });
        let answer = solve(&system, &thermal, &ids, &caloric, &boundaries).unwrap();
        // The native 0.2 W/K link transfers 12 W between 320 K and 380 K.
        // Original specific-volume slopes are 5e-6 and 15e-6 m³/J.
        near(answer.heat.net_power_w[&ids[0]], 12.0, 1e-12);
        near(answer.heat.net_power_w[&ids[1]], -12.0, 1e-12);
        near(answer.heat.external_power_w, 0.0, 1e-12);
        near(answer.cells[0].volume_rate_m3_s, 6e-5, 1e-15);
        near(answer.cells[1].volume_rate_m3_s, -18e-5, 1e-15);
        near(answer.expected_volume_rate_m3_s, -12e-5, 1e-15);
        near(answer.boundary_volume_rate_m3_s, -12e-5, 2e-11);
        if let Some(first) = &first_field {
            assert_eq!(&answer.flow.flow.field.x, first);
        } else {
            first_field = Some(answer.flow.flow.field.x);
        }
    }
}

#[test]
fn sealed_heating_refuses_without_inventing_compressibility_or_advancing_owners() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let caloric = curve(800.0, SolidLiquidPhase::Liquid);
    let (thermal, ids) = initial_thermal(&mesh, &caloric, false, |_, cell| {
        cell.mass_kg = 1000.0 * cell.geometry.volume_m3;
        cell.boundary.source_w = 40.0 * cell.mass_kg;
    });
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| StokesBoundary::Velocity { edge, value: &zero })
        .collect();
    let before = thermal.cells().clone();
    assert!(matches!(
        solve(&system, &thermal, &ids, &caloric, &boundaries),
        Err(StokesThermalDilationError::Flow(
            CaloricDilationError::Flow(StokesBoundaryError::IncompatibleFlux { .. })
        ))
    ));
    assert_eq!(*thermal.cells(), before);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn original_chart_current_inventory_storage_and_liquid_phase_are_required() {
    struct InternalEnergy<'a>(&'a EquilibriumEnthalpyPhaseCurve);
    impl ThermalLaw for InternalEnergy<'_> {
        fn storage_basis(&self) -> ThermalStorageBasis {
            ThermalStorageBasis::InternalEnergy
        }
        fn caloric_identity(&self) -> Option<[u8; 32]> {
            ThermalLaw::caloric_identity(self.0)
        }
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            ThermalLaw::temperature_k(self.0, h)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            ThermalLaw::maximum_temperature_slope(self.0)
        }
    }
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let boundaries = open_boundaries(&mesh);
    let caloric = curve(800.0, SolidLiquidPhase::Liquid);
    let wrong_chart = curve(700.0, SolidLiquidPhase::Liquid);
    let (thermal, ids) = initial_thermal(&mesh, &caloric, false, |_, cell| {
        cell.mass_kg = 1000.0 * cell.geometry.volume_m3;
        cell.boundary.source_w = 40.0 * cell.mass_kg;
    });
    assert!(matches!(
        solve(&system, &thermal, &ids, &wrong_chart, &boundaries),
        Err(StokesThermalDilationError::Binding(_))
    ));
    // Same material card and temperature are insufficient: the entire caloric
    // chart, including its volume law and energy datum, must match.
    assert_ne!(caloric.identity(), wrong_chart.identity());
    let duplicate = vec![ids[0]; ids.len()];
    assert!(matches!(
        solve(&system, &thermal, &duplicate, &caloric, &boundaries),
        Err(StokesThermalDilationError::Binding(_))
    ));
    let warmed = thermal
        .advance_to(0.1, thermal_config(), || Ok(()))
        .unwrap()
        .0;
    // Native thermal reference geometry cannot masquerade as the updated
    // parcel volume after heating at fixed mass.
    assert!(matches!(
        solve(&system, &warmed, &ids, &caloric, &boundaries),
        Err(StokesThermalDilationError::Binding(_))
    ));
    let internal_law = InternalEnergy(&caloric);
    let (internal, internal_ids) = initial_thermal(&mesh, &internal_law, false, |_, cell| {
        cell.mass_kg = 1000.0 * cell.geometry.volume_m3;
    });
    assert!(matches!(
        solve(&system, &internal, &internal_ids, &caloric, &boundaries),
        Err(StokesThermalDilationError::Binding(_))
    ));
    let solid = curve(800.0, SolidLiquidPhase::Solid);
    let (solid_thermal, solid_ids) = initial_thermal(&mesh, &solid, false, |_, cell| {
        cell.mass_kg = 1000.0 * cell.geometry.volume_m3;
    });
    assert!(matches!(
        solve(&system, &solid_thermal, &solid_ids, &solid, &boundaries),
        Err(StokesThermalDilationError::Flow(
            CaloricDilationError::NonLiquid { .. }
        ))
    ));
    let charts = BTreeMap::from([(MaterialId(9), &caloric)]);
    assert!(matches!(
        solve_at_thermal_state(
            &system,
            &thermal,
            StokesThermalDilationBinding {
                cells: binding(&ids),
                depth_m: DEPTH,
                caloric_by_material: &charts,
            },
            FluxParams::default(),
            &|_, _| 2.0,
            &zero,
            &boundaries,
            StokesThermalDilationControls {
                maximum_dilation_defect_per_s: 2e-11,
                maximum_volume_rate_defect_m3_s: 2e-11,
            },
            || true,
        ),
        Err(StokesThermalDilationError::Cancelled)
    ));
}
