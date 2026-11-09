//! Two finite walls retain one original thermal inventory and physical tick.

use super::*;
use fs_couple::free_surface_wall::planar::pair::thermal::FreeSurfaceThermalPlanarWallPair;
use fs_couple::free_surface_wall::planar::pair::{FreeSurfacePlanarWallPair, PlanarWallPairError};

fn pair_parameters(normal: usize) -> ([MassProperties; 2], [RigidBodyState; 2], [PlanarWall3; 2]) {
    let planes = [
        [PlanarWall3::XMin, PlanarWall3::XMax],
        [PlanarWall3::YMin, PlanarWall3::YMax],
        [PlanarWall3::ZMin, PlanarWall3::ZMax],
    ][normal];
    let (first_properties, first) = body_parameters(planes[0], 0.7);
    let (second_properties, second) = body_parameters(planes[1], 2.2);
    let mut center = components(second.pose().position_world());
    center[(normal + 1) % 3] = 6.0 * DX;
    center[(normal + 2) % 3] = 5.3 * DX;
    let second = state(
        second_properties.mass(),
        components(second_properties.principal_inertia_body())[normal],
        [-0.009, 0.004, -0.8],
        vector(center),
        normal,
    );
    (
        [first_properties, second_properties],
        [first, second],
        planes,
    )
}

// G1/G3: constant material laws reduce exactly to the original isothermal pair,
// including a thermal history advanced before either body is attached.
#[test]
fn g1_constant_thermal_pair_recovers_original_mechanics_on_every_axis() {
    for normal in 0..3 {
        let (properties, bodies, planes) = pair_parameters(normal);
        with_thermal_domain(planes[0], 0.1, false, true, |mut liquid, gate| {
            liquid.step(gate).unwrap();
            let initial_steps = liquid.accepted_steps();
            let initial_time = liquid.physical_time_s();
            let initial_h = liquid.total_enthalpy_j().unwrap();
            let mut plain = FreeSurfacePlanarWallPair::new(
                liquid.flow().clone(),
                properties,
                bodies,
                planes,
                CELLS,
                gate,
            )
            .unwrap();
            let mut thermal = FreeSurfaceThermalPlanarWallPair::new(
                liquid.clone(),
                properties,
                bodies,
                planes,
                CELLS,
                gate,
            )
            .unwrap();
            let forces = [
                tangent([0.0003, -0.0002, 0.0], normal),
                tangent([-0.0001, 0.0004, 0.0], normal),
            ];
            for index in 1..=2 {
                let reference = plain
                    .step_candidate(forces, [1e-7, -2e-7], controls(), gate)
                    .unwrap();
                let actual = thermal
                    .step_candidate(forces, [1e-7, -2e-7], &[], controls(), gate)
                    .unwrap();
                assert_eq!(actual.motion.bodies(), reference.motion.bodies());
                assert_eq!(actual.fluid, reference.fluid);
                assert_eq!(
                    actual.motion.liquid().flow().lattice().grid().state_bits(),
                    reference.motion.liquid().lattice().grid().state_bits()
                );
                assert_eq!(actual.native_trials, 2);
                assert_eq!(actual.trial_cells, 2 * CELLS);
                assert_eq!(actual.motion.accepted_steps(), index);
                assert_eq!(
                    actual.motion.liquid().accepted_steps(),
                    initial_steps + index
                );
                near(actual.motion.elapsed_time_s(), index as f64 * DT, 1e-15);
                near(
                    actual.motion.liquid().physical_time_s(),
                    initial_time + index as f64 * DT,
                    1e-14,
                );
                near(
                    actual.motion.liquid().total_enthalpy_j().unwrap(),
                    initial_h,
                    2e-10,
                );
                assert_eq!(actual.liquid_report.enthalpy.external_heat_j, 0.0);
                assert_eq!(actual.liquid_report.enthalpy.viscous_heat_j, 0.0);
                assert!(actual.liquid_report.molecular_heat.is_none());
                plain = reference.motion;
                thermal = actual.motion;
            }
        });
    }
}

// G1: actual endpoint momentum/kinetic changes, original local material readback,
// and one original thermal replay prove the joint heated physical trajectory.
#[test]
fn g1_heated_buoyant_pair_closes_actual_body_exchange_and_one_heat_inventory() {
    for buoyant in [false, true] {
        with_thermal_domain(PlanarWall3::ZMin, 0.02, buoyant, true, |liquid, gate| {
            let (properties, bodies, planes) = pair_parameters(2);
            let mut motion = FreeSurfaceThermalPlanarWallPair::new(
                liquid, properties, bodies, planes, CELLS, gate,
            )
            .unwrap();
            let forces = [
                Vec3::new(0.0003, -0.0002, 0.0),
                Vec3::new(-0.0001, 0.0004, 0.0),
            ];
            let torques = [1e-7, -2e-7];
            let heat = [
                FreeSurfaceCellHeat {
                    coordinates: [3, 3, 1],
                    energy_j: 0.01,
                },
                FreeSurfaceCellHeat {
                    coordinates: [4, 2, 6],
                    energy_j: 0.007,
                },
            ];
            for index in 1..=3 {
                let before = motion.bodies();
                let before_h = motion.liquid().total_enthalpy_j().unwrap();
                let before_mass = motion.liquid().flow().ledger_mass_kg().unwrap();
                let result = motion
                    .step_candidate(forces, torques, &heat, controls(), gate)
                    .unwrap();
                let mut body_linear = Vec3::ZERO;
                let mut body_angular = Vec3::ZERO;
                let mut source_linear = Vec3::ZERO;
                let mut source_angular = Vec3::ZERO;
                let mut dk = 0.0;
                let mut work = 0.0;
                for i in 0..2 {
                    let after = result.motion.bodies()[i];
                    let account = &result.walls[i];
                    body_linear = body_linear.add(
                        after
                            .linear_momentum_world()
                            .sub(before[i].linear_momentum_world()),
                    );
                    body_angular = body_angular.add(angular(after).sub(angular(before[i])));
                    let linear = forces[i].scale(DT).add(account.raw_guide_impulse_world);
                    let couple = Vec3::new(0.0, 0.0, torques[i] * DT)
                        .add(account.raw_guide_angular_impulse_world);
                    source_linear = source_linear.add(linear);
                    source_angular = source_angular
                        .add(account.coupling_origin_world_m.cross(linear))
                        .add(couple);
                    dk += after.kinetic_energy(properties[i]).unwrap()
                        - before[i].kinetic_energy(properties[i]).unwrap();
                    let midpoint_v = before[i]
                        .linear_momentum_world()
                        .add(after.linear_momentum_world())
                        .scale(0.5 / properties[i].mass());
                    let spin0 = before[i]
                        .pose()
                        .orientation()
                        .rotate_body_to_world(before[i].angular_momentum_body())
                        .z;
                    let spin1 = after
                        .pose()
                        .orientation()
                        .rotate_body_to_world(after.angular_momentum_body())
                        .z;
                    let midpoint_w =
                        0.5 * (spin0 + spin1) / properties[i].principal_inertia_body().z;
                    work += linear.dot(midpoint_v) + couple.z * midpoint_w;
                }
                close_vector(
                    body_linear.add(vector(result.fluid.fluid_impulse_kg_m_s)),
                    source_linear,
                    3e-17,
                );
                close_vector(
                    body_angular.add(vector(
                        result.fluid.fluid_angular_impulse_at_si([0.0; 3]).unwrap(),
                    )),
                    source_angular,
                    3e-19,
                );
                near(dk + result.fluid.fluid_work_j, work, 3e-18);
                near(
                    result.motion.liquid().total_enthalpy_j().unwrap() - before_h,
                    0.017,
                    2e-10,
                );
                near(
                    result.motion.liquid().flow().ledger_mass_kg().unwrap(),
                    before_mass,
                    2e-16,
                );
                near(result.liquid_report.enthalpy.external_heat_j, 0.017, 1e-17);
                assert_eq!(result.liquid_report.enthalpy.viscous_heat_j, 0.0);
                // The first reflection may only accelerate the rest fluid;
                // every later accepted tick must transport enthalpy.
                if index > 1 {
                    assert!(result.liquid_report.enthalpy.mass_routes > 0);
                }
                assert_eq!(result.motion.liquid().accepted_steps(), index);
                for cell in &heat {
                    let h = result
                        .motion
                        .liquid()
                        .specific_enthalpy_at(cell.coordinates)
                        .unwrap();
                    let mu = result
                        .motion
                        .liquid()
                        .material_at(cell.coordinates)
                        .unwrap()
                        .dynamic_viscosity_pa_s();
                    near(mu, 0.1 - 0.08 * h / 400000.0, 1e-14);
                }
                let (replayed, replay) = motion
                    .liquid()
                    .planar_rigid_wall_pair_candidate_si(
                        result.fluid.walls.map(|w| w.motion),
                        controls().maximum_relative_mass_defect,
                        &heat,
                        gate,
                    )
                    .unwrap();
                assert_eq!(replay.walls, result.fluid);
                assert_eq!(replay.flow, result.liquid_report.enthalpy);
                assert_eq!(
                    replayed.flow().lattice().grid().state_bits(),
                    result.motion.liquid().flow().lattice().grid().state_bits()
                );
                assert_eq!(
                    replayed.total_enthalpy_j().unwrap(),
                    result.motion.liquid().total_enthalpy_j().unwrap()
                );
                for cell in &heat {
                    assert_eq!(
                        replayed.material_at(cell.coordinates),
                        result.motion.liquid().material_at(cell.coordinates)
                    );
                }
                motion = result.motion;
            }
        });
    }
}

// G3: the original whole-liquid molecular sampler supplies one source. Its
// heat is not multiplied by two walls or by the discarded common pilot.
#[test]
fn g3_molecular_heat_matches_original_source_once_for_both_bodies() {
    with_thermal_domain(PlanarWall3::ZMin, 0.02, false, true, |liquid, gate| {
        let (properties, bodies, planes) = pair_parameters(2);
        let motion =
            FreeSurfaceThermalPlanarWallPair::new(liquid, properties, bodies, planes, CELLS, gate)
                .unwrap()
                .step_candidate([Vec3::ZERO; 2], [0.0; 2], &[], controls(), gate)
                .unwrap()
                .motion;
        let mut original = motion.liquid().clone();
        let reference = original.step_with_viscous_heating(&[], gate).unwrap();
        let molecular: f64 = reference.applied.iter().map(|cell| cell.heat_j).sum();
        assert!(molecular > 0.0);
        let authored: Vec<_> = reference
            .applied
            .iter()
            .map(|cell| FreeSurfaceCellHeat {
                coordinates: cell.coordinates,
                energy_j: cell.heat_j,
            })
            .collect();
        let implicit = motion
            .step_candidate_with_viscous_heating([Vec3::ZERO; 2], [0.0; 2], &[], controls(), gate)
            .unwrap();
        let explicit = motion
            .step_candidate([Vec3::ZERO; 2], [0.0; 2], &authored, controls(), gate)
            .unwrap();
        assert_eq!(
            implicit.liquid_report.molecular_heat.as_ref().unwrap(),
            &reference.applied
        );
        assert_eq!(implicit.liquid_report.enthalpy.viscous_heat_j, molecular);
        assert_eq!(implicit.liquid_report.enthalpy.external_heat_j, 0.0);
        near(
            implicit.motion.liquid().total_enthalpy_j().unwrap()
                - motion.liquid().total_enthalpy_j().unwrap(),
            molecular,
            5e-12,
        );
        assert_eq!(implicit.motion.bodies(), explicit.motion.bodies());
        assert_eq!(implicit.fluid, explicit.fluid);
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

// G4: late source or second-body rejection cannot publish the other body's
// staged endpoint, the common thermal history, or either accepted clock.
#[test]
fn g4_late_second_body_and_thermal_refusals_leave_all_owners_retryable() {
    with_thermal_domain(PlanarWall3::ZMin, 0.02, false, true, |liquid, gate| {
        let (properties, mut bodies, planes) = pair_parameters(2);
        bodies[0] = state(
            properties[0].mass(),
            properties[0].principal_inertia_body().z,
            [0.0; 3],
            Vec3::new(4.0 * DX, 4.0 * DX, -2.0 * DX),
            2,
        );
        let motion =
            FreeSurfaceThermalPlanarWallPair::new(liquid, properties, bodies, planes, CELLS, gate)
                .unwrap();
        let reference = motion.clone();
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
            energy_j: 0.01,
        }];
        let mut limited = controls();
        limited.maximum_trial_cells = 2 * CELLS - 1;
        assert!(matches!(
            motion.step_candidate([Vec3::ZERO; 2], [0.0; 2], &heat, limited, gate),
            Err(PlanarWallPairError::Wall {
                error: PlanarWallError::WorkBudget { .. },
                ..
            })
        ));
        limited = controls();
        limited.maximum_rotation_rad = 1e-6;
        assert!(matches!(
            motion.step_candidate([Vec3::ZERO; 2], [0.0; 2], &heat, limited, gate),
            Err(PlanarWallPairError::Wall {
                index: 1,
                error: PlanarWallError::Invalid(_)
            })
        ));
        limited = controls();
        limited.maximum_absolute_wall_work_j = 0.0;
        assert!(matches!(
            motion.step_candidate([Vec3::ZERO; 2], [0.0; 2], &heat, limited, gate),
            Err(PlanarWallPairError::Invalid(_))
        ));
        assert!(matches!(
            motion.step_candidate(
                [Vec3::ZERO; 2],
                [0.0; 2],
                &[FreeSurfaceCellHeat {
                    coordinates: [3, 3, 1],
                    energy_j: 10.0
                }],
                controls(),
                gate
            ),
            Err(PlanarWallPairError::Thermal(_))
        ));
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(matches!(
            motion.step_candidate([Vec3::ZERO; 2], [0.0; 2], &heat, controls(), &cancelled),
            Err(PlanarWallPairError::Cancelled)
        ));
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
            .step_candidate([Vec3::ZERO; 2], [0.0; 2], &heat, controls(), gate)
            .unwrap();
        let fresh = reference
            .step_candidate([Vec3::ZERO; 2], [0.0; 2], &heat, controls(), gate)
            .unwrap();
        assert_eq!(retry.motion.bodies(), fresh.motion.bodies());
        assert_eq!(retry.fluid, fresh.fluid);
        assert_eq!(retry.liquid_report, fresh.liquid_report);
        assert_eq!(
            retry.motion.liquid().flow().lattice().grid().state_bits(),
            fresh.motion.liquid().flow().lattice().grid().state_bits()
        );
    });
}
