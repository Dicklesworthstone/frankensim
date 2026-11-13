//! Manufactured native rigid/wall oracles, not measured rotor calibration.

use super::{CELL_MASS, CELLS, DT, DX, liquid, near, source};
use fs_couple::free_surface_wall::rotor::{
    FreeSurfacePlanarRotor, PlanarRotorControls, PlanarRotorError,
};
use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::PlanarWall3;
use fs_lbm::material::free_surface::MaterialFreeSurface3;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

// Sum (x²+y²)/3 over the 8×8 receiving wall-cell centers relative to (4,4).
const ROTATIONAL_IMPEDANCE: f64 = 224.0 * CELL_MASS * DX * DX;

fn controls() -> PlanarRotorControls {
    PlanarRotorControls {
        maximum_trial_cells: 2 * CELLS,
        maximum_relative_mass_defect: 1e-12,
        maximum_rotation_rad: 0.2,
        maximum_absolute_wall_work_j: 1e-6,
        angular_momentum_absolute_kg_m2_s: 1e-20,
        midpoint_absolute_rad_s: 1e-10,
        work_absolute_j: 1e-19,
        relative_tolerance: 1e-10,
    }
}

fn properties(inertia: f64) -> MassProperties {
    MassProperties::new(
        1e-3,
        Vec3::ZERO,
        Vec3::new(0.5 * inertia, 0.5 * inertia, inertia),
    )
    .unwrap()
}

fn body(inertia: f64, omega: f64, origin: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(origin, UnitQuaternion::IDENTITY).unwrap(),
        Vec3::ZERO,
        Vec3::new(0.0, 0.0, inertia * omega),
    )
    .unwrap()
}

fn center() -> Vec3 {
    Vec3::new(4.0 * DX, 4.0 * DX, 0.0)
}

fn liquid_spin(liquid: &MaterialFreeSurface3<'_>, origin: Vec3) -> Vec3 {
    let mut result = Vec3::ZERO;
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if let Some(sample) = liquid.sample_at_si(x, y, z).unwrap() {
                    let mass = sample.fill_fraction * sample.density_kg_m3 * DX * DX * DX;
                    let arm = Vec3::new(
                        (x as f64 + 0.5) * DX,
                        (y as f64 + 0.5) * DX,
                        (z as f64 + 0.5) * DX,
                    )
                    .sub(origin);
                    result = result.add(
                        arm.cross(
                            Vec3::new(
                                sample.velocity_m_s[0],
                                sample.velocity_m_s[1],
                                sample.velocity_m_s[2],
                            )
                            .scale(mass),
                        ),
                    );
                }
            }
        }
    }
    result
}

#[test]
fn finite_rotor_matches_independent_impedance_spin_work_and_native_axial_attitude() {
    let material = source(0.1);
    let gate = CancelGate::new_clock_free();
    for inertia in [0.25, 0.5, 2.0, 1000.0].map(|ratio| ratio * ROTATIONAL_IMPEDANCE) {
        let fluid = liquid(&material, None);
        let before_spin = liquid_spin(&fluid, center());
        let bits = fluid.lattice().grid().state_bits();
        let initial_body = body(inertia, 2.0, center());
        let rotor = FreeSurfacePlanarRotor::new(
            fluid,
            properties(inertia),
            initial_body,
            PlanarWall3::ZMin,
            CELLS,
            &gate,
        )
        .unwrap();
        let step = rotor.step_candidate(0.0, controls(), &gate).unwrap();
        let midpoint = 2.0 * inertia * 2.0 / (2.0 * inertia + ROTATIONAL_IMPEDANCE);
        let endpoint = 2.0 * midpoint - 2.0;
        near(step.midpoint_angular_speed_rad_s, midpoint, 2e-11);
        near(step.motion.angular_speed_rad_s().unwrap(), endpoint, 2e-11);
        near(step.wall.fluid_response[2][2], ROTATIONAL_IMPEDANCE, 1e-23);
        near(
            step.wall.fluid_angular_impulse_kg_m2_s[2],
            ROTATIONAL_IMPEDANCE * midpoint,
            1e-20,
        );
        near(
            liquid_spin(step.motion.liquid(), center()).z - before_spin.z,
            step.wall.fluid_angular_impulse_kg_m2_s[2],
            2e-18,
        );
        near(
            step.wall.fluid_work_j,
            ROTATIONAL_IMPEDANCE * midpoint * midpoint,
            2e-19,
        );
        near(
            step.kinetic_energy_change_j,
            0.5 * inertia * (endpoint * endpoint - 4.0),
            2e-19,
        );
        near(
            step.kinetic_energy_change_j + step.wall.fluid_work_j,
            0.0,
            2e-19,
        );
        near(step.mechanical_work_residual_j, 0.0, 2e-19);
        let expected =
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), midpoint * DT).unwrap();
        for (actual, expected) in step
            .motion
            .body()
            .pose()
            .orientation()
            .components()
            .into_iter()
            .zip(expected.components())
        {
            near(actual, expected, 1e-14);
        }
        assert_eq!(step.motion.body().pose().position_world(), center());
        assert_eq!(step.motion.body().linear_momentum_world(), Vec3::ZERO);
        assert_eq!(step.native_trials, 2);
        assert_eq!(step.trial_cells, 2 * CELLS);
        assert_eq!(step.motion.liquid().lattice().grid().steps(), 1);
        assert_eq!(step.motion.accepted_steps(), 1);
        assert_eq!(step.motion.elapsed_time_s(), DT);
        assert_eq!(step.motion.liquid().scaling().material(), &material);
        assert_eq!(rotor.body(), initial_body);
        assert_eq!(rotor.liquid().lattice().grid().state_bits(), bits);
    }
}

#[test]
fn accepted_rotor_replays_native_liquid_and_continues_with_actuator_work() {
    let material = source(0.1);
    let gate = CancelGate::new_clock_free();
    let inertia = 3.0 * ROTATIONAL_IMPEDANCE;
    let mut rotor = FreeSurfacePlanarRotor::new(
        liquid(&material, None),
        properties(inertia),
        body(inertia, 1.5, center()),
        PlanarWall3::ZMin,
        CELLS,
        &gate,
    )
    .unwrap();
    for (i, torque) in [1e-7, -2e-7, 0.0, 1e-7].into_iter().enumerate() {
        let mut replay = rotor.liquid().clone();
        let before = rotor.body();
        let step = rotor.step_candidate(torque, controls(), &gate).unwrap();
        let report = replay
            .step_with_planar_rigid_wall_motion_si(
                step.wall.motion,
                controls().maximum_relative_mass_defect,
                &gate,
            )
            .unwrap();
        assert_eq!(report, step.wall);
        assert_eq!(
            replay.lattice().grid().state_bits(),
            step.motion.liquid().lattice().grid().state_bits()
        );
        let omega_before = properties(inertia)
            .angular_velocity_body_checked(before.angular_momentum_body())
            .unwrap()
            .z;
        let omega_after = step.motion.angular_speed_rad_s().unwrap();
        near(
            step.external_work_j,
            torque * DT * (0.5 * omega_before + 0.5 * omega_after),
            1e-22,
        );
        near(
            step.kinetic_energy_change_j + step.wall.fluid_work_j
                - step.external_work_j
                - step.raw_guide_work_j,
            0.0,
            2e-19,
        );
        assert_eq!(step.mechanics.state_after, step.motion.body());
        assert_eq!(step.mechanics.state_before, before);
        assert_eq!(step.motion.accepted_steps(), (i + 1) as u64);
        assert_eq!(
            step.motion.liquid().lattice().grid().steps(),
            (i + 1) as u64
        );
        near(step.motion.elapsed_time_s(), (i + 1) as f64 * DT, 1e-18);
        rotor = step.motion;
    }
}

#[test]
fn pressure_reference_changes_physical_support_force_and_moment_without_changing_dynamics() {
    let material = source(0.1);
    let gate = CancelGate::new_clock_free();
    let inertia = 20.0 * ROTATIONAL_IMPEDANCE;
    let origin = Vec3::new(DX, 2.0 * DX, -DX);
    let make = |pressure| {
        FreeSurfacePlanarRotor::new(
            liquid(&material, pressure),
            properties(inertia),
            body(inertia, 1.0, origin),
            PlanarWall3::ZMin,
            CELLS,
            &gate,
        )
        .unwrap()
        .step_candidate(0.0, controls(), &gate)
        .unwrap()
    };
    let no_reference = make(None);
    let first = make(Some(100_000.0));
    let shifted = make(Some(101_000.0));
    assert!(no_reference.absolute_guide_angular_impulse_world.is_none());
    assert_eq!(first.motion.body(), shifted.motion.body());
    assert_eq!(
        first.motion.liquid().lattice().grid().state_bits(),
        shifted.motion.liquid().lattice().grid().state_bits()
    );
    assert_eq!(
        first.raw_guide_impulse_world,
        shifted.raw_guide_impulse_world
    );
    assert_eq!(
        first.raw_guide_angular_impulse_world,
        shifted.raw_guide_angular_impulse_world
    );
    let force_change = shifted
        .absolute_guide_impulse_world
        .unwrap()
        .sub(first.absolute_guide_impulse_world.unwrap());
    let expected_force = Vec3::new(0.0, 0.0, 1000.0 * 64.0 * DX * DX * DT);
    near(force_change.z, expected_force.z, 1e-17);
    let expected_moment = center().sub(origin).cross(expected_force);
    let moment_change = shifted
        .absolute_guide_angular_impulse_world
        .unwrap()
        .sub(first.absolute_guide_angular_impulse_world.unwrap());
    near(moment_change.x, expected_moment.x, 1e-19);
    near(moment_change.y, expected_moment.y, 1e-19);
    near(moment_change.z, 0.0, 0.0);
}

#[test]
fn unsupported_states_budgets_and_late_endpoint_refusal_are_atomic() {
    let material = source(0.1);
    let gate = CancelGate::new_clock_free();
    let inertia = 2.0 * ROTATIONAL_IMPEDANCE;
    let moving = RigidBodyState::new(
        Pose::new(center(), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(1e-5, 0.0, 0.0),
        Vec3::ZERO,
    )
    .unwrap();
    assert!(
        FreeSurfacePlanarRotor::new(
            liquid(&material, None),
            properties(inertia),
            moving,
            PlanarWall3::ZMin,
            CELLS,
            &gate
        )
        .is_err()
    );
    let tilted = RigidBodyState::new(
        Pose::new(
            center(),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), 0.2).unwrap(),
        )
        .unwrap(),
        Vec3::ZERO,
        Vec3::ZERO,
    )
    .unwrap();
    assert!(
        FreeSurfacePlanarRotor::new(
            liquid(&material, None),
            properties(inertia),
            tilted,
            PlanarWall3::ZMin,
            CELLS,
            &gate
        )
        .is_err()
    );
    assert!(
        FreeSurfacePlanarRotor::new(
            liquid(&material, None),
            properties(inertia),
            body(inertia, 40.0, center()),
            PlanarWall3::ZMin,
            CELLS,
            &gate
        )
        .is_err()
    );
    let rotor = FreeSurfacePlanarRotor::new(
        liquid(&material, None),
        properties(inertia),
        body(inertia, 0.0, center()),
        PlanarWall3::ZMin,
        CELLS,
        &gate,
    )
    .unwrap();
    let before = rotor.body();
    let bits = rotor.liquid().lattice().grid().state_bits();
    let mut limited = controls();
    limited.maximum_trial_cells = CELLS;
    assert!(
        matches!(rotor.step_candidate(0.0, limited, &gate), Err(PlanarRotorError::WorkBudget { required_trial_cells, maximum }) if required_trial_cells == 2 * CELLS && maximum == CELLS)
    );
    // A 12 rad/s midpoint is locally admissible; its 24 rad/s endpoint exceeds
    // the original 0.1 m/s source envelope. This fails AFTER final fluid replay.
    let torque = 12.0 * (2.0 * inertia + ROTATIONAL_IMPEDANCE) / DT;
    assert!(
        matches!(rotor.step_candidate(torque, controls(), &gate), Err(PlanarRotorError::Invalid(reason)) if reason.contains("endpoint"))
    );
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(matches!(
        rotor.step_candidate(0.0, controls(), &stopped),
        Err(PlanarRotorError::Cancelled)
    ));
    limited = controls();
    limited.maximum_absolute_wall_work_j = 0.0;
    assert!(rotor.step_candidate(1e-7, limited, &gate).is_err());
    limited = controls();
    limited.maximum_rotation_rad = 1e-12;
    assert!(rotor.step_candidate(1e-7, limited, &gate).is_err());
    assert_eq!(rotor.body(), before);
    assert_eq!(rotor.liquid().lattice().grid().state_bits(), bits);
    assert_eq!(rotor.accepted_steps(), 0);
    let first = rotor.step_candidate(1e-7, controls(), &gate).unwrap();
    let retry = rotor.step_candidate(1e-7, controls(), &gate).unwrap();
    assert_eq!(first.motion.body(), retry.motion.body());
    assert_eq!(
        first.motion.liquid().lattice().grid().state_bits(),
        retry.motion.liquid().lattice().grid().state_bits()
    );
}
