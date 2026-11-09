//! Prescribed gas loading reaches finite bodies through one thermal trajectory.

use super::*;
use fs_couple::free_surface_wall::planar::pair::thermal::FreeSurfaceThermalPlanarWallPair;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceWallGasPressure, FreeSurfaceWallSources,
};
use std::cell::Cell;

const REFERENCE: f64 = 100000.0;

fn with_pressure<R>(
    pocket: bool,
    exercise: impl FnOnce(MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    with_thermal_cells(
        0.02,
        false,
        Some(REFERENCE),
        |x, y, z| {
            if pocket {
                !((2..6).contains(&x) && (2..6).contains(&y) && (2..6).contains(&z))
            } else {
                z < 2
            }
        },
        exercise,
    )
}

fn pressure(position: [f64; 3]) -> f64 {
    REFERENCE + 200.0 * (position[0] - 4.0 * DX) + 100.0 * (position[1] - 3.0 * DX)
}

fn sources<'a>(
    heat: &'a [FreeSurfaceCellHeat],
    pressure_pa: &'a dyn Fn([f64; 3]) -> f64,
) -> FreeSurfaceWallSources<'a> {
    FreeSurfaceWallSources {
        heat,
        molecular_heating: false,
        gas_pressure: Some(FreeSurfaceWallGasPressure::Cell {
            maximum_density_excursion: 0.05,
            pressure_pa,
        }),
    }
}

fn pair_motion<'s, 'm>(
    liquid: MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gate: &CancelGate,
    moving: bool,
) -> (
    FreeSurfaceThermalPlanarWallPair<'s, 'm>,
    [MassProperties; 2],
) {
    let planes = [PlanarWall3::ZMin, PlanarWall3::ZMax];
    let (p0, b0) = body_parameters(planes[0], 0.7);
    let (p1, b1) = body_parameters(planes[1], 1.8);
    let properties = [p0, p1];
    let mut bodies = [b0, b1];
    if !moving {
        for i in 0..2 {
            bodies[i] = state(
                properties[i].mass(),
                properties[i].principal_inertia_body().z,
                [0.0; 3],
                bodies[i].pose().position_world(),
                2,
            );
        }
    }
    (
        FreeSurfaceThermalPlanarWallPair::new(liquid, properties, bodies, planes, CELLS, gate)
            .unwrap(),
        properties,
    )
}

#[test]
fn g0_reference_pressure_preserves_original_single_and_pair_thermal_ticks() {
    let heat = [FreeSurfaceCellHeat {
        coordinates: [3, 3, 1],
        energy_j: 0.002,
    }];
    with_pressure(false, |liquid, gate| {
        let (properties, body) = body_parameters(PlanarWall3::ZMin, 0.7);
        let motion = FreeSurfaceThermalPlanarWall::new(
            liquid,
            properties,
            body,
            PlanarWall3::ZMin,
            CELLS,
            gate,
        )
        .unwrap();
        let original = motion
            .step_candidate(Vec3::ZERO, 0.0, &heat, controls(), gate)
            .unwrap();
        let calls = Cell::new(0);
        let reference = |_| {
            calls.set(calls.get() + 1);
            REFERENCE
        };
        let actual = motion
            .step_candidate_with_sources(
                Vec3::ZERO,
                0.0,
                sources(&heat, &reference),
                controls(),
                gate,
            )
            .unwrap();
        assert!(calls.get() > 0, "the slab must contain actual gas links");
        assert_eq!(actual.motion.body(), original.motion.body());
        assert_eq!(actual.wall, original.wall);
        assert_eq!(actual.liquid_report, original.liquid_report);
        assert_eq!(
            actual.motion.liquid().flow().lattice().grid().state_bits(),
            original
                .motion
                .liquid()
                .flow()
                .lattice()
                .grid()
                .state_bits()
        );
    });
    with_pressure(true, |liquid, gate| {
        let (motion, _) = pair_motion(liquid, gate, true);
        let original = motion
            .step_candidate([Vec3::ZERO; 2], [0.0; 2], &heat, controls(), gate)
            .unwrap();
        let calls = Cell::new(0);
        let reference = |_| {
            calls.set(calls.get() + 1);
            REFERENCE
        };
        let actual = motion
            .step_candidate_with_sources(
                [Vec3::ZERO; 2],
                [0.0; 2],
                sources(&heat, &reference),
                controls(),
                gate,
            )
            .unwrap();
        assert!(
            calls.get() > 0,
            "the enclosed pocket must contain gas links"
        );
        let link = motion
            .step_candidate_with_sources(
                [Vec3::ZERO; 2],
                [0.0; 2],
                FreeSurfaceWallSources {
                    heat: &heat,
                    gas_pressure: Some(FreeSurfaceWallGasPressure::Link {
                        maximum_density_excursion: 0.05,
                        pressure_pa: &|_| REFERENCE,
                    }),
                    ..FreeSurfaceWallSources::default()
                },
                controls(),
                gate,
            )
            .unwrap();
        for result in [actual, link] {
            assert_eq!(result.motion.bodies(), original.motion.bodies());
            assert_eq!(result.fluid, original.fluid);
            assert_eq!(result.liquid_report, original.liquid_report);
            assert_eq!(
                result.motion.liquid().flow().lattice().grid().state_bits(),
                original
                    .motion
                    .liquid()
                    .flow()
                    .lattice()
                    .grid()
                    .state_bits()
            );
        }
    });
}

#[test]
fn g1_pressure_drives_delayed_pair_recoil_with_original_momentum_and_heat_accounts() {
    with_pressure(true, |liquid, gate| {
        let (mut loaded, properties) = pair_motion(liquid, gate, false);
        let mut unloaded = loaded.clone();
        let heat = [FreeSurfaceCellHeat {
            coordinates: [3, 3, 1],
            energy_j: 0.002,
        }];
        for index in 1..=5 {
            let before = loaded.bodies();
            let before_h = loaded.liquid().total_enthalpy_j().unwrap();
            let actual = loaded
                .step_candidate_with_sources(
                    [Vec3::ZERO; 2],
                    [0.0; 2],
                    sources(&heat, &pressure),
                    controls(),
                    gate,
                )
                .unwrap();
            let reference = unloaded
                .step_candidate([Vec3::ZERO; 2], [0.0; 2], &heat, controls(), gate)
                .unwrap();
            let mut body_linear = Vec3::ZERO;
            let mut body_angular = Vec3::ZERO;
            let mut guide_linear = Vec3::ZERO;
            let mut guide_angular = Vec3::ZERO;
            let mut kinetic = 0.0;
            for i in 0..2 {
                let after = actual.motion.bodies()[i];
                let account = &actual.walls[i];
                body_linear = body_linear.add(
                    after
                        .linear_momentum_world()
                        .sub(before[i].linear_momentum_world()),
                );
                body_angular = body_angular.add(angular(after).sub(angular(before[i])));
                guide_linear = guide_linear.add(account.raw_guide_impulse_world);
                guide_angular = guide_angular
                    .add(
                        account
                            .coupling_origin_world_m
                            .cross(account.raw_guide_impulse_world),
                    )
                    .add(account.raw_guide_angular_impulse_world);
                kinetic += after.kinetic_energy(properties[i]).unwrap()
                    - before[i].kinetic_energy(properties[i]).unwrap();
            }
            close_vector(
                body_linear.add(vector(actual.fluid.fluid_impulse_kg_m_s)),
                guide_linear,
                3e-17,
            );
            close_vector(
                body_angular.add(vector(
                    actual.fluid.fluid_angular_impulse_at_si([0.0; 3]).unwrap(),
                )),
                guide_angular,
                3e-19,
            );
            near(
                kinetic + actual.fluid.fluid_work_j - actual.raw_guide_work_j,
                0.0,
                3e-18,
            );
            near(
                actual.motion.liquid().total_enthalpy_j().unwrap() - before_h,
                0.002,
                2e-10,
            );
            assert_eq!(actual.liquid_report.enthalpy.external_heat_j, 0.002);
            assert_eq!(actual.liquid_report.enthalpy.viscous_heat_j, 0.0);
            assert!(actual.liquid_report.molecular_heat.is_none());
            assert_eq!(actual.motion.accepted_steps(), index);
            assert_eq!(actual.motion.liquid().accepted_steps(), index);
            assert_eq!(actual.native_trials, 2);
            assert_eq!(actual.trial_cells, 2 * CELLS);
            near(
                actual.motion.liquid().physical_time_s(),
                37.0 + index as f64 * DT,
                1e-14,
            );
            assert_ne!(
                actual.motion.liquid().flow().lattice().grid().state_bits(),
                reference
                    .motion
                    .liquid()
                    .flow()
                    .lattice()
                    .grid()
                    .state_bits()
            );
            if index == 5 {
                assert!(actual.liquid_report.enthalpy.mass_routes > 0);
                let mut received_x = 0.0;
                for i in 0..2 {
                    let difference = actual.motion.bodies()[i]
                        .linear_momentum_world()
                        .sub(reference.motion.bodies()[i].linear_momentum_world());
                    assert!(
                        difference.norm_squared() > 1e-36,
                        "pressure must reach each finite body"
                    );
                    received_x += difference.x;
                }
                assert!(
                    received_x > 1e-18,
                    "higher pocket pressure on +x drives positive recoil"
                );
                let (replayed, report) = loaded
                    .liquid()
                    .planar_rigid_wall_pair_candidate_with_sources_si(
                        actual.fluid.walls.map(|wall| wall.motion),
                        controls().maximum_relative_mass_defect,
                        sources(&heat, &pressure),
                        gate,
                    )
                    .unwrap();
                assert_eq!(report.flow.walls, actual.fluid);
                assert_eq!(report.flow.flow, actual.liquid_report.enthalpy);
                assert_eq!(
                    replayed.flow().lattice().grid().state_bits(),
                    actual.motion.liquid().flow().lattice().grid().state_bits()
                );
                assert_eq!(
                    replayed.total_enthalpy_j().unwrap(),
                    actual.motion.liquid().total_enthalpy_j().unwrap()
                );
            }
            loaded = actual.motion;
            unloaded = reference.motion;
        }
    });
}

#[test]
fn g3_pressure_and_single_body_replay_apply_original_molecular_heat_once() {
    with_pressure(false, |liquid, gate| {
        let (properties, body) = body_parameters(PlanarWall3::ZMin, 0.7);
        let motion = FreeSurfaceThermalPlanarWall::new(
            liquid,
            properties,
            body,
            PlanarWall3::ZMin,
            CELLS,
            gate,
        )
        .unwrap()
        .step_candidate_with_sources(Vec3::ZERO, 0.0, sources(&[], &pressure), controls(), gate)
        .unwrap()
        .motion;
        let mut original = motion.liquid().clone();
        let sampled = original.step_with_viscous_heating(&[], gate).unwrap();
        let molecular: f64 = sampled.applied.iter().map(|cell| cell.heat_j).sum();
        assert!(molecular > 0.0);
        let combined: Vec<_> = sampled
            .applied
            .iter()
            .map(|cell| FreeSurfaceCellHeat {
                coordinates: cell.coordinates,
                energy_j: cell.heat_j,
            })
            .collect();
        let implicit = motion
            .step_candidate_with_sources(
                Vec3::ZERO,
                0.0,
                FreeSurfaceWallSources {
                    molecular_heating: true,
                    ..sources(&[], &pressure)
                },
                controls(),
                gate,
            )
            .unwrap();
        let explicit = motion
            .step_candidate_with_sources(
                Vec3::ZERO,
                0.0,
                sources(&combined, &pressure),
                controls(),
                gate,
            )
            .unwrap();
        assert_eq!(
            implicit.liquid_report.molecular_heat.as_ref().unwrap(),
            &sampled.applied
        );
        assert_eq!(implicit.liquid_report.enthalpy.viscous_heat_j, molecular);
        assert_eq!(implicit.liquid_report.enthalpy.external_heat_j, 0.0);
        near(
            implicit.motion.liquid().total_enthalpy_j().unwrap()
                - motion.liquid().total_enthalpy_j().unwrap(),
            molecular,
            2e-10,
        );
        assert_eq!(implicit.motion.body(), explicit.motion.body());
        assert_eq!(implicit.wall, explicit.wall);
        assert_eq!(
            implicit
                .motion
                .liquid()
                .flow()
                .lattice()
                .grid()
                .state_bits(),
            explicit
                .motion
                .liquid()
                .flow()
                .lattice()
                .grid()
                .state_bits()
        );
        assert_eq!(
            implicit.motion.liquid().total_enthalpy_j().unwrap(),
            explicit.motion.liquid().total_enthalpy_j().unwrap()
        );
        assert_eq!(
            implicit.motion.liquid().accepted_steps(),
            motion.liquid().accepted_steps() + 1
        );
        assert_eq!(
            implicit.motion.accepted_steps(),
            motion.accepted_steps() + 1
        );
    });
}

#[test]
fn g4_late_pressure_and_cancel_refusals_preserve_both_bodies_and_thermal_retry() {
    with_pressure(true, |liquid, gate| {
        let (motion, _) = pair_motion(liquid, gate, true);
        let before = (
            motion.bodies(),
            motion.liquid().flow().lattice().grid().state_bits(),
            motion.liquid().total_enthalpy_j().unwrap(),
            motion.liquid().material_at([3, 3, 1]).cloned(),
            motion.liquid().physical_time_s(),
            motion.liquid().accepted_steps(),
            motion.elapsed_time_s(),
            motion.accepted_steps(),
        );
        let heat = [FreeSurfaceCellHeat {
            coordinates: [3, 3, 1],
            energy_j: 0.002,
        }];
        let count = Cell::new(0);
        let counted = |position| {
            count.set(count.get() + 1);
            pressure(position)
        };
        let reference = motion
            .step_candidate_with_sources(
                [Vec3::ZERO; 2],
                [0.0; 2],
                sources(&heat, &counted),
                controls(),
                gate,
            )
            .unwrap();
        let late = count.get() - 2;
        assert!(late > 4);
        count.set(0);
        let invalid = |position| {
            count.set(count.get() + 1);
            if count.get() == late {
                f64::NAN
            } else {
                pressure(position)
            }
        };
        assert!(
            motion
                .step_candidate_with_sources(
                    [Vec3::ZERO; 2],
                    [0.0; 2],
                    sources(&heat, &invalid),
                    controls(),
                    gate,
                )
                .is_err()
        );
        assert_eq!(count.get(), late, "refusal must reach the final replay");
        let cancelled = CancelGate::new_clock_free();
        count.set(0);
        let interrupt = |position| {
            count.set(count.get() + 1);
            if count.get() == late {
                cancelled.request();
            }
            pressure(position)
        };
        assert!(
            motion
                .step_candidate_with_sources(
                    [Vec3::ZERO; 2],
                    [0.0; 2],
                    sources(&heat, &interrupt),
                    controls(),
                    &cancelled,
                )
                .is_err()
        );
        assert!(count.get() >= late);
        assert_eq!(
            (
                motion.bodies(),
                motion.liquid().flow().lattice().grid().state_bits(),
                motion.liquid().total_enthalpy_j().unwrap(),
                motion.liquid().material_at([3, 3, 1]).cloned(),
                motion.liquid().physical_time_s(),
                motion.liquid().accepted_steps(),
                motion.elapsed_time_s(),
                motion.accepted_steps(),
            ),
            before
        );
        let retry = motion
            .step_candidate_with_sources(
                [Vec3::ZERO; 2],
                [0.0; 2],
                sources(&heat, &pressure),
                controls(),
                gate,
            )
            .unwrap();
        assert_eq!(retry.motion.bodies(), reference.motion.bodies());
        assert_eq!(retry.fluid, reference.fluid);
        assert_eq!(retry.liquid_report, reference.liquid_report);
        assert_eq!(
            retry.motion.liquid().flow().lattice().grid().state_bits(),
            reference
                .motion
                .liquid()
                .flow()
                .lattice()
                .grid()
                .state_bits()
        );
    });
}
