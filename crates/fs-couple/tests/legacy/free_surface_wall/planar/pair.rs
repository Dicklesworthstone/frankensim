//! Six finite wall modes, one native liquid history, and independent reflection algebra.

use super::{
    AREA, CELLS, DT, DX, IMPEDANCE, PlanarWall3, PlanarWallError, ResolvedLiquidSurfaceState,
    SPIN_IMPEDANCE, angular, axis, close_vector, components, controls, domain, inertia, near,
    source, state, tangent, vector,
};
use fs_couple::free_surface_wall::planar::pair::{FreeSurfacePlanarWallPair, PlanarWallPairError};
use fs_exec::CancelGate;
use fs_lbm::material::free_surface::{MaterialFreeSurface3, plan_free_surface_from_material};
use fs_mbd::{MassProperties, RigidBodyState, Vec3};

fn flooded(source: &ResolvedLiquidSurfaceState, pressure: Option<f64>) -> MaterialFreeSurface3<'_> {
    let plan = plan_free_surface_from_material(source, DX, DT, 0.1).unwrap();
    let plan = match pressure {
        Some(p) => plan.with_pressure_reference(p).unwrap(),
        None => plan,
    };
    plan.free_surface3(domain(), |_, _, _| true).unwrap()
}

fn bodies(normal: usize) -> ([MassProperties; 2], [RigidBodyState; 2]) {
    let masses = [0.7 * IMPEDANCE, 3.0 * IMPEDANCE];
    let spins = [1.1 * SPIN_IMPEDANCE, 2.7 * SPIN_IMPEDANCE];
    let properties = core::array::from_fn(|i| inertia(masses[i], spins[i], normal));
    let states = core::array::from_fn(|i| {
        let mut center = [0.0; 3];
        center[(normal + 1) % 3] = [1.6, 6.0][i] * DX;
        center[(normal + 2) % 3] = [2.1, 5.3][i] * DX;
        center[normal] = [-2.0, 10.0][i] * DX;
        state(
            masses[i],
            spins[i],
            [[0.012, -0.006, 1.3], [-0.008, 0.004, -0.7]][i],
            vector(center),
            normal,
        )
    });
    (properties, states)
}

// G0/G3: independent D3Q19 moments and scalar elimination determine all six
// midpoint modes. Shared native replay proves one liquid tick for both bodies.
#[test]
fn all_opposing_axes_match_six_mode_reflection_and_single_tick_oracles() {
    let material = source(0.1);
    let gate = CancelGate::new_clock_free();
    for planes in [
        [PlanarWall3::XMin, PlanarWall3::XMax],
        [PlanarWall3::YMin, PlanarWall3::YMax],
        [PlanarWall3::ZMin, PlanarWall3::ZMax],
    ] {
        let normal = axis(planes[0]);
        let a = (normal + 1) % 3;
        let b = (normal + 2) % 3;
        let (properties, before) = bodies(normal);
        let motion = FreeSurfacePlanarWallPair::new(
            flooded(&material, None),
            properties,
            before,
            planes,
            CELLS,
            &gate,
        )
        .unwrap();
        let q0 = motion
            .twists()
            .unwrap()
            .map(|(v, w)| [components(v)[a], components(v)[b], w]);
        let forces = [
            tangent([0.0003, -0.0002, 0.0], normal),
            tangent([-0.0001, 0.0004, 0.0], normal),
        ];
        let torques = [1e-7, -1.5e-7];
        let result = motion
            .step_candidate(forces, torques, controls(), &gate)
            .unwrap();
        let mut body_linear_change = Vec3::ZERO;
        let mut body_angular_change = Vec3::ZERO;
        let mut source_and_guide_linear = Vec3::ZERO;
        let mut source_and_guide_angular = Vec3::ZERO;
        let mut body_kinetic_change = 0.0;
        let mut source_and_guide_work = 0.0;
        for i in 0..2 {
            let account = &result.walls[i];
            let origin = before[i]
                .pose()
                .position_world()
                .add(tangent(q0[i], normal).scale(0.5 * DT));
            close_vector(account.coupling_origin_world_m, origin, 1e-17);
            let ra = 4.0 * DX - components(origin)[a];
            let rb = 4.0 * DX - components(origin)[b];
            let c = [
                [IMPEDANCE, 0.0, -IMPEDANCE * rb],
                [0.0, IMPEDANCE, IMPEDANCE * ra],
                [
                    -IMPEDANCE * rb,
                    IMPEDANCE * ra,
                    SPIN_IMPEDANCE + IMPEDANCE * (ra * ra + rb * rb),
                ],
            ];
            let m = properties[i].mass();
            let spin = components(properties[i].principal_inertia_body())[normal];
            let rhs = [
                m * q0[i][0] + 0.5 * components(forces[i])[a] * DT,
                m * q0[i][1] + 0.5 * components(forces[i])[b] * DT,
                spin * q0[i][2] + 0.5 * torques[i] * DT,
            ];
            let translation = m + 0.5 * IMPEDANCE;
            let ca = 0.5 * c[0][2];
            let cb = 0.5 * c[1][2];
            let omega = (rhs[2] - (ca * rhs[0] + cb * rhs[1]) / translation)
                / (spin + 0.5 * c[2][2] - (ca * ca + cb * cb) / translation);
            let qmid = [
                (rhs[0] - ca * omega) / translation,
                (rhs[1] - cb * omega) / translation,
                omega,
            ];
            for row in 0..3 {
                for col in 0..3 {
                    near(account.wall.fluid_response[row][col], c[row][col], 1e-18);
                }
            }
            close_vector(account.midpoint_velocity_m_s, tangent(qmid, normal), 2e-12);
            near(account.midpoint_angular_speed_rad_s, omega, 2e-10);
            let expected_work: f64 = (0..3)
                .map(|r| qmid[r] * (0..3).map(|k| c[r][k] * qmid[k]).sum::<f64>())
                .sum();
            near(account.wall.fluid_work_j, expected_work, 1e-18);
            let endpoint = result.motion.twists().unwrap()[i];
            let q1 = core::array::from_fn(|k| 2.0 * qmid[k] - q0[i][k]);
            close_vector(endpoint.0, tangent(q1, normal), 2e-12);
            near(endpoint.1, q1[2], 2e-10);
            let impulse = vector(account.wall.wall_impulse_kg_m_s)
                .add(account.raw_guide_impulse_world)
                .add(account.external_linear_impulse_world);
            let couple = vector(account.wall.wall_angular_impulse_kg_m2_s)
                .add(account.raw_guide_angular_impulse_world)
                .add(account.external_couple_impulse_world);
            close_vector(
                angular(result.motion.bodies()[i]).sub(angular(before[i])),
                origin.cross(impulse).add(couple),
                2e-19,
            );
            body_linear_change = body_linear_change.add(
                result.motion.bodies()[i]
                    .linear_momentum_world()
                    .sub(before[i].linear_momentum_world()),
            );
            body_angular_change =
                body_angular_change.add(angular(result.motion.bodies()[i]).sub(angular(before[i])));
            let source_and_guide_impulse = account.raw_guide_impulse_world.add(forces[i].scale(DT));
            let mut applied_couple = [0.0; 3];
            applied_couple[normal] = torques[i] * DT;
            let source_and_guide_couple = account
                .raw_guide_angular_impulse_world
                .add(vector(applied_couple));
            source_and_guide_linear = source_and_guide_linear.add(source_and_guide_impulse);
            source_and_guide_angular = source_and_guide_angular
                .add(origin.cross(source_and_guide_impulse))
                .add(source_and_guide_couple);
            body_kinetic_change += result.motion.bodies()[i]
                .kinetic_energy(properties[i])
                .unwrap()
                - before[i].kinetic_energy(properties[i]).unwrap();
            let actual_midpoint_velocity = result.motion.bodies()[i]
                .linear_momentum_world()
                .add(before[i].linear_momentum_world())
                .scale(0.5 / m);
            let actual_midpoint_spin = 0.5 * (endpoint.1 + q0[i][2]);
            source_and_guide_work += source_and_guide_impulse.dot(actual_midpoint_velocity)
                + components(source_and_guide_couple)[normal] * actual_midpoint_spin;
        }
        // These joint identities are formed independently from actual native
        // endpoint states and held physical sources, not receipt residuals.
        close_vector(
            body_linear_change.add(vector(result.fluid.fluid_impulse_kg_m_s)),
            source_and_guide_linear,
            3e-17,
        );
        close_vector(
            body_angular_change.add(vector(
                result.fluid.fluid_angular_impulse_at_si([0.0; 3]).unwrap(),
            )),
            source_and_guide_angular,
            3e-19,
        );
        near(
            body_kinetic_change + result.fluid.fluid_work_j,
            source_and_guide_work,
            3e-18,
        );
        let (replay, report) = motion
            .liquid()
            .planar_rigid_wall_pair_candidate_si(
                result.fluid.walls.map(|w| w.motion),
                controls().maximum_relative_mass_defect,
                &gate,
            )
            .unwrap();
        assert_eq!(report, result.fluid);
        assert_eq!(
            replay.lattice().grid().state_bits(),
            result.motion.liquid().lattice().grid().state_bits()
        );
        assert_eq!(motion.bodies(), before);
        assert_eq!(result.motion.elapsed_time_s(), DT);
        assert_eq!(result.motion.accepted_steps(), 1);
        assert_eq!(result.native_trials, 2);
        assert_eq!(result.trial_cells, 2 * CELLS);
        close_vector(result.momentum_residual_kg_m_s, Vec3::ZERO, 2e-17);
        close_vector(result.angular_momentum_residual_kg_m2_s, Vec3::ZERO, 2e-19);
        near(result.mechanical_work_residual_j, 0.0, 2e-18);
    }
}

// G3: both bodies retain their own source ordering and physical pressure datum;
// changing that datum alters only the corresponding guide-force accounts.
#[test]
fn wall_order_and_pressure_reference_preserve_the_joint_mechanics() {
    let material = source(0.1);
    let gate = CancelGate::new_clock_free();
    let (properties, states) = bodies(2);
    let planes = [PlanarWall3::ZMin, PlanarWall3::ZMax];
    let forces = [
        Vec3::new(0.0001, -0.0002, 0.0),
        Vec3::new(-0.0003, 0.0001, 0.0),
    ];
    let torques = [1e-7, -2e-7];
    let make = |pressure| {
        FreeSurfacePlanarWallPair::new(
            flooded(&material, Some(pressure)),
            properties,
            states,
            planes,
            CELLS,
            &gate,
        )
        .unwrap()
    };
    let low = make(100_000.0)
        .step_candidate(forces, torques, controls(), &gate)
        .unwrap();
    let high = make(200_000.0)
        .step_candidate(forces, torques, controls(), &gate)
        .unwrap();
    assert_eq!(low.motion.bodies(), high.motion.bodies());
    assert_eq!(
        low.motion.liquid().lattice().grid().state_bits(),
        high.motion.liquid().lattice().grid().state_bits()
    );
    for i in 0..2 {
        assert_eq!(
            low.walls[i].raw_guide_impulse_world,
            high.walls[i].raw_guide_impulse_world
        );
        let difference = high.walls[i]
            .absolute_guide_impulse_world
            .unwrap()
            .sub(low.walls[i].absolute_guide_impulse_world.unwrap());
        near(difference.z.abs(), 100_000.0 * AREA * DT, 1e-15);
    }
    let reverse = FreeSurfacePlanarWallPair::new(
        flooded(&material, Some(100_000.0)),
        [properties[1], properties[0]],
        [states[1], states[0]],
        [planes[1], planes[0]],
        CELLS,
        &gate,
    )
    .unwrap()
    .step_candidate(
        [forces[1], forces[0]],
        [torques[1], torques[0]],
        controls(),
        &gate,
    )
    .unwrap();
    assert_eq!(
        reverse.motion.bodies(),
        [low.motion.bodies()[1], low.motion.bodies()[0]]
    );
    assert_eq!(
        reverse.motion.liquid().lattice().grid().state_bits(),
        low.motion.liquid().lattice().grid().state_bits()
    );
}

// G1: one initially moving finite wall drives the other through actual evolving
// populations; each accepted physical interval advances the shared clock once.
#[test]
fn continuing_liquid_transmits_momentum_between_both_finite_walls() {
    let material = source(1.0);
    let gate = CancelGate::new_clock_free();
    let masses = [50.0 * IMPEDANCE, 7.0 * IMPEDANCE];
    let spins = [4.0 * SPIN_IMPEDANCE; 2];
    let properties = core::array::from_fn(|i| inertia(masses[i], spins[i], 2));
    let states = [
        state(
            masses[0],
            spins[0],
            [0.025, 0.0, 0.0],
            Vec3::new(4.0 * DX, 4.0 * DX, -DX),
            2,
        ),
        state(
            masses[1],
            spins[1],
            [0.0; 3],
            Vec3::new(4.0 * DX, 4.0 * DX, 9.0 * DX),
            2,
        ),
    ];
    let mut motion = FreeSurfacePlanarWallPair::new(
        flooded(&material, None),
        properties,
        states,
        [PlanarWall3::ZMin, PlanarWall3::ZMax],
        CELLS,
        &gate,
    )
    .unwrap();
    for i in 0..80 {
        let result = motion
            .step_candidate([Vec3::ZERO; 2], [0.0; 2], controls(), &gate)
            .unwrap();
        assert_eq!(result.motion.accepted_steps(), i + 1);
        near(result.motion.elapsed_time_s(), (i + 1) as f64 * DT, 1e-15);
        near(result.mechanical_work_residual_j, 0.0, 1e-18);
        if i == 0 {
            close_vector(result.motion.twists().unwrap()[1].0, Vec3::ZERO, 1e-14);
        }
        motion = result.motion;
    }
    let velocities = motion.twists().unwrap();
    assert!(
        velocities[1].0.x.abs() > 1e-6,
        "finite upper wall must receive transported momentum"
    );
    assert!(
        velocities[0].0.x < 0.025,
        "finite lower wall must recoil against the liquid"
    );
}

// G0: invalid second-wall data and late second-body admission cannot publish
// the already staged first body, shared liquid, or a partial physical interval.
#[test]
fn second_wall_failures_and_cancellation_leave_every_owner_retryable() {
    let material = source(0.1);
    let gate = CancelGate::new_clock_free();
    let properties = [inertia(2.0 * IMPEDANCE, 2.0 * SPIN_IMPEDANCE, 2); 2];
    let states = [
        state(
            properties[0].mass(),
            2.0 * SPIN_IMPEDANCE,
            [0.0; 3],
            Vec3::new(4.0 * DX, 4.0 * DX, -DX),
            2,
        ),
        state(
            properties[1].mass(),
            2.0 * SPIN_IMPEDANCE,
            [0.0; 3],
            Vec3::new(4.0 * DX, 4.0 * DX, 9.0 * DX),
            2,
        ),
    ];
    let planes = [PlanarWall3::ZMin, PlanarWall3::ZMax];
    let motion = FreeSurfacePlanarWallPair::new(
        flooded(&material, None),
        properties,
        states,
        planes,
        CELLS,
        &gate,
    )
    .unwrap();
    let bits = motion.liquid().lattice().grid().state_bits();
    assert!(matches!(
        motion.step_candidate(
            [Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0)],
            [0.0; 2],
            controls(),
            &gate
        ),
        Err(PlanarWallPairError::Wall {
            index: 1,
            error: PlanarWallError::Invalid(_)
        })
    ));
    let mut late = controls();
    late.maximum_rotation_rad = 1e-6;
    assert!(matches!(
        motion.step_candidate([Vec3::ZERO; 2], [0.0, 1e-7], late, &gate),
        Err(PlanarWallPairError::Wall {
            index: 1,
            error: PlanarWallError::Invalid(_)
        })
    ));
    let mut short = controls();
    short.maximum_trial_cells = 2 * CELLS - 1;
    assert!(matches!(
        motion.step_candidate([Vec3::ZERO; 2], [0.0; 2], short, &gate),
        Err(PlanarWallPairError::Wall {
            error: PlanarWallError::WorkBudget { .. },
            ..
        })
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        motion.step_candidate([Vec3::ZERO; 2], [0.0; 2], controls(), &cancelled),
        Err(PlanarWallPairError::Cancelled)
    ));
    assert_eq!(motion.bodies(), states);
    assert_eq!(motion.liquid().lattice().grid().state_bits(), bits);
    assert_eq!(motion.elapsed_time_s(), 0.0);
    assert_eq!(motion.accepted_steps(), 0);
    let retry = motion
        .step_candidate([Vec3::ZERO; 2], [0.0, 1e-7], controls(), &gate)
        .unwrap();
    let reference = motion
        .step_candidate([Vec3::ZERO; 2], [0.0, 1e-7], controls(), &gate)
        .unwrap();
    assert_eq!(retry.motion.bodies(), reference.motion.bodies());
    assert_eq!(
        retry.motion.liquid().lattice().grid().state_bits(),
        reference.motion.liquid().lattice().grid().state_bits()
    );
}
