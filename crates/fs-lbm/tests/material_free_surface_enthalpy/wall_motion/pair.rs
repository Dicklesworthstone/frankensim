//! G1/G3/G4: both wall exchanges share actual thermal source and mass routes.

use super::*;
use fs_lbm::material::free_surface::MaterialPlanarRigidWallMotion3;

fn motions(axis: usize) -> [MaterialPlanarRigidWallMotion3; 2] {
    let planes = [
        [PlanarWall3::XMin, PlanarWall3::XMax],
        [PlanarWall3::YMin, PlanarWall3::YMax],
        [PlanarWall3::ZMin, PlanarWall3::ZMax],
    ][axis];
    core::array::from_fn(|end| {
        let mut origin_m = [0.004; 3];
        origin_m[axis] = 0.008 * end as f64;
        let mut velocity_at_origin_m_s = [0.0; 3];
        velocity_at_origin_m_s[(axis + 1) % 3] = if end == 0 { 0.025 } else { -0.015 };
        velocity_at_origin_m_s[(axis + 2) % 3] = if end == 0 { -0.006 } else { 0.004 };
        MaterialPlanarRigidWallMotion3 {
            wall: planes[end],
            origin_m,
            velocity_at_origin_m_s,
            angular_speed_rad_s: if end == 0 { 1.2 } else { -0.8 },
        }
    })
}

fn close(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_zero_pair_is_one_original_thermal_tick_on_all_axes_with_heat_and_buoyancy() {
    for axis in 0..3 {
        for buoyant in [false, true] {
            with_liquid(8, buoyant, |liquid, gate| {
                let original_bits = liquid.flow().lattice().grid().state_bits();
                let original_h = liquid.total_enthalpy_j().unwrap();
                let mut ordinary = liquid.clone();
                let heat = [FreeSurfaceCellHeat {
                    coordinates: [3, 3, 2],
                    energy_j: 0.001,
                }];
                let expected = ordinary.step_with_heat(&heat, gate).unwrap();
                let stationary = motions(axis).map(|motion| MaterialPlanarRigidWallMotion3 {
                    velocity_at_origin_m_s: [0.0; 3],
                    angular_speed_rad_s: 0.0,
                    ..motion
                });
                let (candidate, report) = liquid
                    .planar_rigid_wall_pair_candidate_si(stationary, 1e-11, &heat, gate)
                    .unwrap();
                assert_eq!(report.flow, expected);
                assert_eq!(
                    candidate.flow().lattice().grid().state_bits(),
                    ordinary.flow().lattice().grid().state_bits()
                );
                assert_eq!(
                    candidate.total_enthalpy_j().unwrap(),
                    ordinary.total_enthalpy_j().unwrap()
                );
                assert_eq!(candidate.accepted_steps(), 1);
                assert_eq!(
                    candidate.physical_time_s(),
                    candidate.flow().scaling().timestep_s()
                );
                assert_eq!(report.walls.fluid_work_j, 0.0);
                near(report.walls.mass_after_kg, report.flow.mass_after_kg, 1e-12);
                for p in [[3, 3, 2], [4, 3, 2]] {
                    assert_eq!(
                        candidate.specific_enthalpy_at(p),
                        ordinary.specific_enthalpy_at(p)
                    );
                    assert_eq!(candidate.material_at(p), ordinary.material_at(p));
                }
                assert_eq!(liquid.flow().lattice().grid().state_bits(), original_bits);
                assert_eq!(liquid.total_enthalpy_j().unwrap(), original_h);
                assert_eq!(liquid.accepted_steps(), 0);
            });
        }
    }
}

#[test]
fn g1_heated_pair_uses_advected_source_viscosity_and_buoyancy_without_heating_wall_work() {
    for buoyant in [false, true] {
        with_liquid(8, buoyant, |liquid, gate| {
            let motion = motions(2);
            // Establish native non-equilibrium stress before changing viscosity.
            for _ in 0..2 {
                *liquid = liquid
                    .planar_rigid_wall_pair_candidate_si(motion, 1e-11, &[], gate)
                    .unwrap()
                    .0;
            }
            let before_h = liquid.total_enthalpy_j().unwrap();
            let before_mass = liquid.flow().ledger_mass_kg().unwrap();
            let before_steps = liquid.accepted_steps();
            let heat = [FreeSurfaceCellHeat {
                coordinates: [3, 3, 0],
                energy_j: 0.01,
            }];
            let (unheated, baseline) = liquid
                .planar_rigid_wall_pair_candidate_si(motion, 1e-11, &[], gate)
                .unwrap();
            let (heated, actual) = liquid
                .planar_rigid_wall_pair_candidate_si(motion, 1e-11, &heat, gate)
                .unwrap();
            assert_ne!(
                heated.flow().lattice().grid().state_bits(),
                unheated.flow().lattice().grid().state_bits()
            );
            assert_ne!(
                actual.walls.walls[0].fluid_impulse_kg_m_s,
                baseline.walls.walls[0].fluid_impulse_kg_m_s
            );
            assert_eq!(actual.flow.external_heat_j, heat[0].energy_j);
            assert_eq!(actual.flow.viscous_heat_j, 0.0);
            assert_eq!(actual.flow.body_heat_j, 0.0);
            near(
                actual.flow.enthalpy_after_j,
                before_h + heat[0].energy_j,
                2e-11,
            );
            near(actual.flow.mass_after_kg, before_mass, 2e-12);
            assert!(actual.flow.mass_routes > 0);
            assert_eq!(heated.accepted_steps(), before_steps + 1);
            assert_eq!(liquid.accepted_steps(), before_steps);
            for p in [[3, 3, 0], [4, 3, 0], [3, 3, 7]] {
                let h = heated.specific_enthalpy_at(p).unwrap();
                let temperature = 300.0 + h / 2000.0;
                near(
                    heated.material_at(p).unwrap().dynamic_viscosity_pa_s(),
                    0.1 - 0.0004 * (temperature - 300.0),
                    1e-12,
                );
            }
            let mut total_work = 0.0;
            for wall in actual.walls.walls {
                let power_pair = (0..3)
                    .map(|d| wall.motion.velocity_at_origin_m_s[d] * wall.fluid_impulse_kg_m_s[d])
                    .sum::<f64>()
                    + wall.motion.angular_speed_rad_s * wall.fluid_angular_impulse_kg_m2_s[2];
                close(wall.fluid_work_j, power_pair, 2e-18);
                assert_eq!(wall.wall_work_j, -wall.fluid_work_j);
                assert_eq!(
                    wall.wall_impulse_kg_m_s,
                    wall.fluid_impulse_kg_m_s.map(|v| -v)
                );
                assert_eq!(
                    wall.wall_angular_impulse_kg_m2_s,
                    wall.fluid_angular_impulse_kg_m2_s.map(|v| -v)
                );
                total_work += wall.fluid_work_j;
            }
            close(actual.walls.fluid_work_j, total_work, 2e-18);
            assert_eq!(actual.walls.wall_work_j, -actual.walls.fluid_work_j);
            // The wall exchange can change while extensive heat still closes
            // solely to authored heat. Both work and heat have their own owners.
            assert_ne!(actual.walls.fluid_work_j, 0.0);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), before_h);
        });
    }
}

#[test]
fn g3_pair_molecular_heat_matches_the_original_source_and_is_applied_once() {
    with_liquid(8, true, |liquid, gate| {
        let motion = motions(2);
        *liquid = liquid
            .planar_rigid_wall_pair_candidate_si(motion, 1e-11, &[], gate)
            .unwrap()
            .0;
        let mut source_probe = liquid.clone();
        let source = source_probe.step_with_viscous_heating(&[], gate).unwrap();
        let authored = FreeSurfaceCellHeat {
            coordinates: [2, 3, 2],
            energy_j: 0.001,
        };
        let mut all_heat: Vec<_> = source
            .applied
            .iter()
            .map(|sample| FreeSurfaceCellHeat {
                coordinates: sample.coordinates,
                energy_j: sample.heat_j,
            })
            .collect();
        // Authored cell heat has unique coordinates. Compose both physical
        // sources at the one shared cell before the equivalent source-only
        // replay; a duplicate entry is an input refusal, not added heat.
        let shared = all_heat
            .iter_mut()
            .find(|sample| sample.coordinates == authored.coordinates)
            .unwrap();
        shared.energy_j += authored.energy_j;
        let (expected, reference) = liquid
            .planar_rigid_wall_pair_candidate_si(motion, 1e-11, &all_heat, gate)
            .unwrap();
        let before_h = liquid.total_enthalpy_j().unwrap();
        let before_steps = liquid.accepted_steps();
        let (actual, report) = liquid
            .planar_rigid_wall_pair_candidate_with_viscous_heating_si(
                motion,
                1e-11,
                &[authored],
                gate,
            )
            .unwrap();
        assert_eq!(report.applied, source.applied);
        assert_eq!(report.flow.walls, reference.walls);
        assert_eq!(
            actual.flow().lattice().grid().state_bits(),
            expected.flow().lattice().grid().state_bits()
        );
        assert_eq!(
            actual.total_enthalpy_j().unwrap(),
            expected.total_enthalpy_j().unwrap()
        );
        assert_eq!(report.flow.flow.external_heat_j, authored.energy_j);
        assert!(report.flow.flow.viscous_heat_j > 0.0);
        assert_eq!(report.flow.flow.viscous_heat_j, source.flow.viscous_heat_j);
        near(
            actual.total_enthalpy_j().unwrap(),
            before_h + authored.energy_j + report.flow.flow.viscous_heat_j,
            2e-11,
        );
        assert_eq!(actual.accepted_steps(), before_steps + 1);
        assert_eq!(liquid.accepted_steps(), before_steps);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), before_h);
    });
}

#[test]
fn g4_second_wall_source_failure_and_cancellation_leave_the_thermal_owner_retryable() {
    with_liquid(8, true, |liquid, gate| {
        let motion = motions(2);
        let bits = liquid.flow().lattice().grid().state_bits();
        let h = liquid.total_enthalpy_j().unwrap();
        let state = liquid.material_at([3, 3, 2]).cloned();
        let mut bad = motion;
        bad[1].velocity_at_origin_m_s[2] = 0.001;
        assert!(
            liquid
                .planar_rigid_wall_pair_candidate_si(bad, 1e-11, &[], gate)
                .is_err()
        );
        let excessive_heat = [FreeSurfaceCellHeat {
            coordinates: [3, 3, 2],
            energy_j: 6.0,
        }];
        assert!(
            liquid
                .planar_rigid_wall_pair_candidate_si(motion, 1e-11, &excessive_heat, gate)
                .is_err()
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            liquid
                .planar_rigid_wall_pair_candidate_si(motion, 1e-11, &[], &cancelled)
                .is_err()
        );
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), h);
        assert_eq!(liquid.material_at([3, 3, 2]).cloned(), state);
        assert_eq!(liquid.accepted_steps(), 0);
        let heat = [FreeSurfaceCellHeat {
            coordinates: [3, 3, 2],
            energy_j: 0.001,
        }];
        let (retry, report) = liquid
            .planar_rigid_wall_pair_candidate_si(motion, 1e-11, &heat, gate)
            .unwrap();
        let (reference, repeated) = liquid
            .planar_rigid_wall_pair_candidate_si(motion, 1e-11, &heat, gate)
            .unwrap();
        assert_eq!(report, repeated);
        assert_eq!(
            retry.flow().lattice().grid().state_bits(),
            reference.flow().lattice().grid().state_bits()
        );
        assert_eq!(
            retry.total_enthalpy_j().unwrap(),
            reference.total_enthalpy_j().unwrap()
        );
        assert_eq!(retry.accepted_steps(), 1);
    });
}
