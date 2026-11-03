//! Actual wall populations share thermal transport, sources and acceptance.

use super::*;
use fs_lbm::d3q19::freesurface3::{PlanarWall3, PlanarWallFieldControls3};
use fs_lbm::material::free_surface::MaterialPlanarWallVelocityField3;
use fs_lbm::material::free_surface::enthalpy::MaterialFreeSurfaceEnthalpy3;
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;

fn controls() -> PlanarWallFieldControls3 {
    PlanarWallFieldControls3 {
        maximum_domain_cells: 512,
        maximum_wall_faces: 64,
        maximum_relative_mass_defect: 1e-11,
        maximum_relative_face_mass_defect: 1e-12,
    }
}

fn motion(velocities: &[[f64; 2]]) -> MaterialPlanarWallVelocityField3<'_> {
    MaterialPlanarWallVelocityField3 {
        wall: PlanarWall3::ZMin,
        moment_origin_m: [0.004, 0.004, 0.0],
        tangential_velocity_m_s: velocities,
    }
}

fn field() -> Vec<[f64; 2]> {
    (0..64)
        .map(|i| [0.02 + 0.001 * (i % 8) as f64, -0.001 * (i / 8) as f64])
        .collect()
}

fn with_liquid<R>(
    height: usize,
    buoyant: bool,
    exercise: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let card = card(0.02, if buoyant { 980.0 } else { 1000.0 });
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point(300.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        &card,
        &point(300.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = caloric(&card);
    let gate = CancelGate::new_clock_free();
    let gravity = if buoyant {
        [0.003, 0.0, -0.02]
    } else {
        [0.0; 3]
    };
    let mut cells = initial(true);
    cells.retain(|cell| cell.coordinates[2] < height);
    for cell in &mut cells {
        cell.velocity_m_s = [0.003 * (cell.coordinates[2] as f64 - 3.5), 0.0, 0.0];
    }
    with_cx(&gate, |cx| {
        let flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells(
                MaterialFreeSurfaceDomain3 {
                    gravity_m_s2: gravity,
                    ..domain()
                },
                &cells,
                cx,
            )
            .unwrap();
        let h = |p: [usize; 3]| if p[0] < 4 { 100000.0 } else { 300000.0 };
        let mut liquid = if buoyant {
            flow.with_boussinesq_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                BoussinesqConfig {
                    gravity_m_s2: gravity,
                    max_relative_density_departure: 0.03,
                    max_relative_specific_heat_departure: 0.0,
                },
                h,
                &gate,
            )
        } else {
            flow.with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                h,
                &gate,
            )
        }
        .unwrap();
        exercise(&mut liquid, &gate)
    })
}

#[test]
fn zero_wall_retains_original_thermal_and_buoyant_trajectory() {
    for height in [4, 8] {
        for buoyant in [false, true] {
            with_liquid(height, buoyant, |liquid, gate| {
                let mut expected = liquid.clone();
                let heat = [FreeSurfaceCellHeat {
                    coordinates: [3, 3, 2],
                    energy_j: 0.001,
                }];
                let original = expected.step_with_heat(&heat, gate).unwrap();
                let report = liquid
                    .step_with_planar_wall_velocity_field_si(
                        motion(&[[0.0; 2]; 64]),
                        controls(),
                        &heat,
                        gate,
                    )
                    .unwrap();
                assert_eq!(report.flow, original);
                assert_eq!(
                    liquid.flow().lattice().grid().state_bits(),
                    expected.flow().lattice().grid().state_bits()
                );
                assert_eq!(
                    liquid.total_enthalpy_j().unwrap(),
                    expected.total_enthalpy_j().unwrap()
                );
                assert_eq!(
                    liquid.material_at([3, 3, 2]),
                    expected.material_at([3, 3, 2])
                );
                assert_eq!(report.wall.exchange.fluid_work_j, 0.0);
                assert_eq!(liquid.accepted_steps(), 1);
            });
        }
    }
}

#[test]
fn driven_faces_transport_enthalpy_and_retain_work_without_heating_it() {
    with_liquid(4, false, |liquid, gate| {
        let mut stationary = liquid.clone();
        let velocities = field();
        let initial_h = liquid.total_enthalpy_j().unwrap();
        let initial_mass = liquid.flow().ledger_mass_kg().unwrap();
        for _ in 0..4 {
            stationary.step(gate).unwrap();
            let report = liquid
                .step_with_planar_wall_velocity_field_si(motion(&velocities), controls(), &[], gate)
                .unwrap();
            near(report.flow.enthalpy_after_j, initial_h, 2e-11);
            near(report.flow.mass_after_kg, initial_mass, 1e-11);
            near(report.wall.mass_after_kg, report.flow.mass_after_kg, 1e-12);
            assert!(report.flow.mass_routes > 0);
            assert_eq!(report.flow.external_heat_j, 0.0);
            assert_eq!(report.flow.viscous_heat_j, 0.0);
            assert!(report.wall.exchange.fluid_work_j > 0.0);
            let mut face_work = 0.0;
            for face in &report.wall.faces {
                let expected = face
                    .velocity_m_s
                    .iter()
                    .zip(face.exchange.fluid_impulse_kg_m_s)
                    .map(|(u, impulse)| u * impulse)
                    .sum::<f64>();
                near(face.exchange.fluid_work_j, expected, 1e-12);
                assert_eq!(face.exchange.wall_work_j, -face.exchange.fluid_work_j);
                face_work += face.exchange.fluid_work_j;
            }
            near(report.wall.exchange.fluid_work_j, face_work, 1e-12);
        }
        assert_ne!(
            liquid.flow().lattice().grid().state_bits(),
            stationary.flow().lattice().grid().state_bits()
        );
        // The advected thermal source is live, not the reference 300 K viscosity.
        for p in [[3, 3, 1], [4, 3, 1]] {
            let h = liquid.specific_enthalpy_at(p).unwrap();
            let temperature = 300.0 + h / 2000.0;
            near(
                liquid.material_at(p).unwrap().dynamic_viscosity_pa_s(),
                0.1 - 0.0004 * (temperature - 300.0),
                1e-12,
            );
        }
    });
}

#[test]
fn molecular_wall_heat_matches_original_source_and_is_applied_once() {
    with_liquid(8, true, |liquid, gate| {
        let velocities = field();
        liquid
            .step_with_planar_wall_velocity_field_si(motion(&velocities), controls(), &[], gate)
            .unwrap();
        let mut source_probe = liquid.clone();
        let source = source_probe.step_with_viscous_heating(&[], gate).unwrap();
        let heat: Vec<_> = source
            .applied
            .iter()
            .map(|sample| FreeSurfaceCellHeat {
                coordinates: sample.coordinates,
                energy_j: sample.heat_j,
            })
            .collect();
        let mut authored = liquid.clone();
        let expected = authored
            .step_with_planar_wall_velocity_field_si(motion(&velocities), controls(), &heat, gate)
            .unwrap();
        let before = liquid.total_enthalpy_j().unwrap();
        let actual = liquid
            .step_with_planar_wall_velocity_field_and_viscous_heating_si(
                motion(&velocities),
                controls(),
                &[],
                gate,
            )
            .unwrap();
        assert_eq!(actual.applied, source.applied);
        assert_eq!(actual.flow.wall, expected.wall);
        assert_eq!(
            liquid.flow().lattice().grid().state_bits(),
            authored.flow().lattice().grid().state_bits()
        );
        assert_eq!(
            liquid.total_enthalpy_j().unwrap(),
            authored.total_enthalpy_j().unwrap()
        );
        assert!(actual.flow.flow.viscous_heat_j > 0.0);
        assert_eq!(actual.flow.flow.external_heat_j, 0.0);
        near(
            actual.flow.flow.viscous_heat_j,
            source.flow.viscous_heat_j,
            1e-13,
        );
        assert!(
            (liquid.total_enthalpy_j().unwrap() - before - actual.flow.flow.viscous_heat_j).abs()
                < 1e-8
        );
    });
}
