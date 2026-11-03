//! Independent D3Q19 impedance and native rigid momentum/work checks.

use super::{AREA, CELL_MASS, CELLS, DT, DX, domain, liquid, near, source};
use fs_couple::free_surface_wall::planar::{
    FreeSurfacePlanarWall, PlanarWallControls, PlanarWallError,
};
use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::PlanarWall3;
use fs_lbm::material::free_surface::{MaterialFreeSurface3, plan_free_surface_from_material};
use fs_material::liquid::surface::ResolvedLiquidSurfaceState;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const IMPEDANCE: f64 = CELL_MASS * 64.0 / 3.0;
const SPIN_IMPEDANCE: f64 = 224.0 * CELL_MASS * DX * DX;

fn controls() -> PlanarWallControls {
    PlanarWallControls {
        maximum_trial_cells: 2 * CELLS,
        maximum_relative_mass_defect: 1e-12,
        maximum_translation_m: 2.0 * DX,
        maximum_rotation_rad: 0.2,
        maximum_absolute_wall_work_j: 1e-6,
        momentum_absolute_kg_m_s: 1e-17,
        angular_momentum_absolute_kg_m2_s: 1e-20,
        midpoint_absolute_m_s: 1e-12,
        midpoint_absolute_rad_s: 1e-10,
        work_absolute_j: 1e-18,
        relative_tolerance: 1e-10,
    }
}

fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn axis(wall: PlanarWall3) -> usize {
    match wall {
        PlanarWall3::XMin | PlanarWall3::XMax => 0,
        PlanarWall3::YMin | PlanarWall3::YMax => 1,
        PlanarWall3::ZMin | PlanarWall3::ZMax => 2,
    }
}
fn upper(wall: PlanarWall3) -> bool {
    matches!(
        wall,
        PlanarWall3::XMax | PlanarWall3::YMax | PlanarWall3::ZMax
    )
}
fn tangent(q: [f64; 3], normal: usize) -> Vec3 {
    let mut v = [0.0; 3];
    v[(normal + 1) % 3] = q[0];
    v[(normal + 2) % 3] = q[1];
    vector(v)
}
fn inertia(mass: f64, spin: f64, normal: usize) -> MassProperties {
    let mut diagonal = [0.75 * spin; 3];
    diagonal[normal] = spin;
    MassProperties::new(mass, Vec3::ZERO, vector(diagonal)).unwrap()
}
fn state(mass: f64, spin: f64, q: [f64; 3], center: Vec3, normal: usize) -> RigidBodyState {
    let mut angular = [0.0; 3];
    angular[normal] = spin * q[2];
    RigidBodyState::new(
        Pose::new(center, UnitQuaternion::IDENTITY).unwrap(),
        tangent(q, normal).scale(mass),
        vector(angular),
    )
    .unwrap()
}
fn wall_liquid(source: &ResolvedLiquidSurfaceState, wall: PlanarWall3) -> MaterialFreeSurface3<'_> {
    plan_free_surface_from_material(source, DX, DT, 0.1)
        .unwrap()
        .free_surface3(domain(), |x, y, z| {
            let coordinate = [x, y, z][axis(wall)];
            if upper(wall) {
                coordinate >= 4
            } else {
                coordinate < 4
            }
        })
        .unwrap()
}
fn close_vector(actual: Vec3, expected: Vec3, tolerance: f64) {
    for (a, e) in components(actual).into_iter().zip(components(expected)) {
        near(a, e, tolerance);
    }
}
fn angular(body: RigidBodyState) -> Vec3 {
    body.pose()
        .position_world()
        .cross(body.linear_momentum_world())
        .add(
            body.pose()
                .orientation()
                .rotate_body_to_world(body.angular_momentum_body()),
        )
}

#[test]
fn all_planes_close_joint_translation_spin_against_independent_reflection_oracle() {
    let material = source(0.1);
    let gate = CancelGate::new_clock_free();
    for wall in [
        PlanarWall3::XMin,
        PlanarWall3::XMax,
        PlanarWall3::YMin,
        PlanarWall3::YMax,
        PlanarWall3::ZMin,
        PlanarWall3::ZMax,
    ] {
        let normal = axis(wall);
        let a = (normal + 1) % 3;
        let b = (normal + 2) % 3;
        for ratio in [0.3, 3.0, 100.0] {
            let mass = ratio * IMPEDANCE;
            let spin = 0.6 * ratio * SPIN_IMPEDANCE;
            let q0 = [0.015, -0.007, 1.3];
            let mut center = [0.0; 3];
            center[a] = 1.5 * DX;
            center[b] = 2.3 * DX;
            center[normal] = if upper(wall) { 10.0 * DX } else { -2.0 * DX };
            let center = vector(center);
            let before = state(mass, spin, q0, center, normal);
            let props = inertia(mass, spin, normal);
            let motion = FreeSurfacePlanarWall::new(
                wall_liquid(&material, wall),
                props,
                before,
                wall,
                CELLS,
                &gate,
            )
            .unwrap();
            let force = tangent([0.0003, -0.0002, 0.0], normal);
            let torque = 1e-7;
            let step = motion
                .step_candidate(force, torque, controls(), &gate)
                .unwrap();
            let midpoint_center = center.add(tangent(q0, normal).scale(0.5 * DT));
            close_vector(step.coupling_origin_world_m, midpoint_center, 1e-17);

            // Sum original D3Q19 wall directions/weights on an 8-by-8 plane.
            // At rest, C_tt = rho dx^3 N/3 and C_spin(center) is the polar
            // second moment of the receiving cell centers times rho dx^3/3.
            let arm_a = 4.0 * DX - components(midpoint_center)[a];
            let arm_b = 4.0 * DX - components(midpoint_center)[b];
            let c = [
                [IMPEDANCE, 0.0, -IMPEDANCE * arm_b],
                [0.0, IMPEDANCE, IMPEDANCE * arm_a],
                [
                    -IMPEDANCE * arm_b,
                    IMPEDANCE * arm_a,
                    SPIN_IMPEDANCE + IMPEDANCE * (arm_a * arm_a + arm_b * arm_b),
                ],
            ];
            for i in 0..3 {
                for j in 0..3 {
                    near(step.wall.fluid_response[i][j], c[i][j], 1e-18);
                }
            }
            let rhs = [
                mass * q0[0] + 0.5 * components(force)[a] * DT,
                mass * q0[1] + 0.5 * components(force)[b] * DT,
                spin * q0[2] + 0.5 * torque * DT,
            ];
            let translation = mass + 0.5 * IMPEDANCE;
            let cross_a = 0.5 * c[0][2];
            let cross_b = 0.5 * c[1][2];
            let angular = spin + 0.5 * c[2][2];
            // Independently eliminate both translations; no production matrix solve.
            let omega = (rhs[2] - (cross_a * rhs[0] + cross_b * rhs[1]) / translation)
                / (angular - (cross_a * cross_a + cross_b * cross_b) / translation);
            let qmid = [
                (rhs[0] - cross_a * omega) / translation,
                (rhs[1] - cross_b * omega) / translation,
                omega,
            ];
            close_vector(step.midpoint_velocity_m_s, tangent(qmid, normal), 2e-12);
            near(step.midpoint_angular_speed_rad_s, omega, 2e-10);
            assert!(
                (omega - rhs[2] / angular).abs() > 1e-5,
                "fixture must exercise translation/spin coupling"
            );
            let (velocity, accepted_omega) = step.motion.twist().unwrap();
            let q1 = core::array::from_fn(|i| 2.0 * qmid[i] - q0[i]);
            close_vector(velocity, tangent(q1, normal), 4e-12);
            near(accepted_omega, q1[2], 4e-10);
            close_vector(
                step.motion.body().pose().position_world(),
                center.add(tangent(qmid, normal).scale(DT)),
                1e-16,
            );
            let expected_rotation = UnitQuaternion::from_axis_angle(
                vector(core::array::from_fn(
                    |d| if d == normal { 1.0 } else { 0.0 },
                )),
                omega * DT,
            )
            .unwrap();
            for (actual, expected) in step
                .motion
                .body()
                .pose()
                .orientation()
                .components()
                .into_iter()
                .zip(expected_rotation.components())
            {
                near(actual, expected, 1e-13);
            }
            let expected_work = (0..3)
                .map(|i| qmid[i] * (0..3).map(|j| c[i][j] * qmid[j]).sum::<f64>())
                .sum::<f64>();
            near(step.wall.fluid_work_j, expected_work, 1e-18);
            near(
                step.kinetic_energy_change_j,
                step.motion.body().kinetic_energy(props).unwrap()
                    - before.kinetic_energy(props).unwrap(),
                2e-18,
            );
            near(step.mechanical_work_residual_j, 0.0, 2e-18);
            assert_eq!(step.raw_guide_work_j, 0.0);
            assert_eq!(step.native_trials, 2);
            assert_eq!(step.trial_cells, 2 * CELLS);
            assert_eq!(step.motion.accepted_steps(), 1);
            assert_eq!(step.motion.elapsed_time_s(), DT);
            assert_eq!(step.motion.liquid().scaling().material(), &material);
            assert_eq!(motion.body(), before);
        }
    }
}

#[test]
fn translated_com_continuation_replays_fluid_and_balances_full_orbital_momentum() {
    let material = source(0.1);
    let gate = CancelGate::new_clock_free();
    let mass = 2.0 * IMPEDANCE;
    let spin = 3.0 * SPIN_IMPEDANCE;
    let props = inertia(mass, spin, 2);
    let initial = state(
        mass,
        spin,
        [0.012, 0.006, 1.7],
        Vec3::new(1.5 * DX, 6.0 * DX, -DX),
        2,
    );
    let mut motion = FreeSurfacePlanarWall::new(
        liquid(&material, None),
        props,
        initial,
        PlanarWall3::ZMin,
        CELLS,
        &gate,
    )
    .unwrap();
    for (index, (force, torque)) in [
        (Vec3::new(-0.002, 0.003, 0.0), 2e-7),
        (Vec3::new(0.001, -0.002, 0.0), -1e-7),
        (Vec3::ZERO, 0.0),
    ]
    .into_iter()
    .enumerate()
    {
        let before = motion.body();
        let mut replay = motion.liquid().clone();
        let step = motion
            .step_candidate(force, torque, controls(), &gate)
            .unwrap();
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
        let wall_impulse = vector(step.wall.wall_impulse_kg_m_s);
        let wall_couple = vector(step.wall.wall_angular_impulse_kg_m2_s);
        let total_impulse = wall_impulse
            .add(step.raw_guide_impulse_world)
            .add(step.external_linear_impulse_world);
        let total_couple = wall_couple
            .add(step.raw_guide_angular_impulse_world)
            .add(step.external_couple_impulse_world);
        close_vector(
            step.motion
                .body()
                .linear_momentum_world()
                .sub(before.linear_momentum_world()),
            total_impulse,
            1e-17,
        );
        close_vector(
            angular(step.motion.body()).sub(angular(before)),
            total_couple.add(step.coupling_origin_world_m.cross(total_impulse)),
            1e-19,
        );
        let initial_arm_error = step
            .coupling_origin_world_m
            .sub(before.pose().position_world())
            .cross(total_impulse);
        if index == 0 {
            assert!(
                initial_arm_error.norm_squared() > 1e-36,
                "noncollinear motion must make the free half-drift moment arm observable"
            );
        }
        let v0 = before.linear_momentum_world().scale(1.0 / mass);
        let v1 = step.motion.body().linear_momentum_world().scale(1.0 / mass);
        let w0 = props
            .angular_velocity_body_checked(before.angular_momentum_body())
            .unwrap()
            .z;
        let w1 = step.motion.twist().unwrap().1;
        near(
            step.external_work_j,
            force.dot(v0.add(v1).scale(0.5)) * DT + torque * 0.5 * (w0 + w1) * DT,
            2e-19,
        );
        near(
            step.kinetic_energy_change_j + step.wall.fluid_work_j
                - step.external_work_j
                - step.raw_guide_work_j,
            0.0,
            2e-18,
        );
        assert_eq!(step.motion.accepted_steps(), index as u64 + 1);
        assert_eq!(
            step.motion.liquid().lattice().grid().steps(),
            index as u64 + 1
        );
        motion = step.motion;
    }
}

#[test]
fn physical_guide_pressure_changes_neither_the_joint_twist_nor_native_populations() {
    let material = source(0.1);
    let gate = CancelGate::new_clock_free();
    let mass = 3.0 * IMPEDANCE;
    let spin = 2.0 * SPIN_IMPEDANCE;
    let initial = state(
        mass,
        spin,
        [0.01, -0.005, 1.5],
        Vec3::new(DX, 2.0 * DX, -DX),
        2,
    );
    let make = |pressure| {
        FreeSurfacePlanarWall::new(
            liquid(&material, Some(pressure)),
            inertia(mass, spin, 2),
            initial,
            PlanarWall3::ZMin,
            CELLS,
            &gate,
        )
        .unwrap()
        .step_candidate(Vec3::ZERO, 0.0, controls(), &gate)
        .unwrap()
    };
    let low = make(95_000.0);
    let high = make(105_000.0);
    assert_eq!(low.motion.body(), high.motion.body());
    assert_eq!(
        low.motion.liquid().lattice().grid().state_bits(),
        high.motion.liquid().lattice().grid().state_bits()
    );
    assert_eq!(low.raw_guide_impulse_world, high.raw_guide_impulse_world);
    assert_eq!(
        low.pressure_perturbation_guide_impulse_world,
        high.pressure_perturbation_guide_impulse_world
    );
    let change = Vec3::new(0.0, 0.0, 10_000.0 * AREA * DT);
    close_vector(
        high.absolute_guide_impulse_world
            .unwrap()
            .sub(low.absolute_guide_impulse_world.unwrap()),
        change,
        1e-17,
    );
    let arm = vector(low.wall.wetted_centroid_m).sub(low.coupling_origin_world_m);
    close_vector(
        high.absolute_guide_angular_impulse_world
            .unwrap()
            .sub(low.absolute_guide_angular_impulse_world.unwrap()),
        arm.cross(change),
        1e-19,
    );
    assert_eq!(low.raw_guide_work_j, 0.0);
    assert_eq!(high.raw_guide_work_j, 0.0);
}

#[test]
fn unsupported_guides_and_late_replay_refusals_preserve_all_accepted_state() {
    let material = source(0.1);
    let gate = CancelGate::new_clock_free();
    let props = inertia(1.0, 1e-5, 2);
    let body = state(1.0, 1e-5, [0.0; 3], Vec3::new(4.0 * DX, 4.0 * DX, 0.0), 2);
    let motion = FreeSurfacePlanarWall::new(
        liquid(&material, None),
        props,
        body,
        PlanarWall3::ZMin,
        CELLS,
        &gate,
    )
    .unwrap();
    let bits = motion.liquid().lattice().grid().state_bits();
    let mut short = controls();
    short.maximum_trial_cells = CELLS;
    assert!(matches!(
        motion.step_candidate(Vec3::ZERO, 0.0, short, &gate),
        Err(PlanarWallError::WorkBudget { .. })
    ));
    assert!(matches!(
        motion.step_candidate(Vec3::new(0.0, 0.0, 1.0), 0.0, controls(), &gate),
        Err(PlanarWallError::Invalid(_))
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        motion
            .step_candidate(Vec3::ZERO, 0.0, controls(), &cancelled)
            .err(),
        Some(PlanarWallError::Cancelled)
    );
    // The native replay admits midpoint speed .075 m/s, then the actual native
    // wall endpoint .15 m/s exceeds the material's .1 m/s scenario bound.
    let force = Vec3::new((2.0 + IMPEDANCE) * 0.075 / DT, 0.0, 0.0);
    assert!(
        matches!(motion.step_candidate(force, 0.0, controls(), &gate),
        Err(PlanarWallError::Invalid(reason)) if reason.contains("endpoint"))
    );
    assert_eq!(motion.body(), body);
    assert_eq!(motion.liquid().lattice().grid().state_bits(), bits);
    assert_eq!(motion.elapsed_time_s(), 0.0);
    assert_eq!(motion.accepted_steps(), 0);
    let retry = motion
        .step_candidate(Vec3::new(0.001, -0.002, 0.0), 1e-7, controls(), &gate)
        .unwrap();
    let reference = motion
        .clone()
        .step_candidate(Vec3::new(0.001, -0.002, 0.0), 1e-7, controls(), &gate)
        .unwrap();
    assert_eq!(retry.motion.body(), reference.motion.body());
    assert_eq!(
        retry.motion.liquid().lattice().grid().state_bits(),
        reference.motion.liquid().lattice().grid().state_bits()
    );
    let invalid_body =
        RigidBodyState::new(body.pose(), Vec3::new(0.0, 0.0, 1e-3), Vec3::ZERO).unwrap();
    assert!(matches!(
        FreeSurfacePlanarWall::new(
            liquid(&material, None),
            props,
            invalid_body,
            PlanarWall3::ZMin,
            CELLS,
            &gate
        ),
        Err(PlanarWallError::Invalid(_))
    ));
}
