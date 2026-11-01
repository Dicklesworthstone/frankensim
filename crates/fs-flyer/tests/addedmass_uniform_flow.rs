//! Uniform-current inertia, affine momentum and prescribed-flow source work.
use fs_flyer::addedmass::articulated::{
    LinkedStrip, UniformFlow, UniformFlowAssembly, assemble_articulated_strips,
    assemble_articulated_strips_in_uniform_flow,
};
use fs_flyer::addedmass::{AeroGeneralizedLoads, Strip, solve_articulated_joint};
use fs_flyer::partitioned::articulated::{
    ArticulatedLoads, ArticulatedState, ArticulatedStepConfig, step,
};
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::generalized::free_floating_terms;
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
    forward_kinematics,
};

fn vector(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}
fn close_slice(actual: &[f64], expected: &[f64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (&a, &b) in actual.iter().zip(expected) {
        close(a, b, tolerance);
    }
}
fn pose(rotation: [f64; 3], translation: [f64; 3]) -> Se3 {
    Se3::from_parts(
        So3::exp(So3Tangent::new(vector(rotation))).unwrap(),
        vector(translation),
    )
    .unwrap()
}
fn inertia(mass: f64, center: [f64; 3]) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        vector(center),
        Mat3 {
            m: [
                0.3 * mass,
                0.0,
                0.0,
                0.0,
                0.4 * mass,
                0.0,
                0.0,
                0.0,
                0.5 * mass,
            ],
        },
    )
    .unwrap()
}
fn plate(center: [f64; 3]) -> Strip {
    Strip {
        name: "finite plate in uniform current",
        chord_m: 0.8,
        span_m: 1.4,
        position_m: center,
        normal: [0.0, 0.0, 1.0],
        span_axis: [0.0, 1.0, 0.0],
        control_coord: None,
        control_gain: 0.0,
        control_angular_gain: [0.0; 3],
    }
}
fn rigid_model(mass: f64) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(mass, [0.0; 3]),
    )])
    .unwrap()
}
fn flow_assembly(
    model: &ArticulatedModel,
    base: FreeFloatingBaseState,
    q: &[f64],
    rates: &[f64],
    strips: &[LinkedStrip],
    flow: UniformFlow,
) -> UniformFlowAssembly {
    let dry = free_floating_terms(model, base, q, rates, || false).unwrap();
    assemble_articulated_strips_in_uniform_flow(1.18, strips, &dry.kinematics, flow).unwrap()
}
fn generalized_velocity(base: FreeFloatingBaseState, rates: &[f64]) -> Vec<f64> {
    let mut velocity = array(base.twist_body.linear).to_vec();
    velocity.extend_from_slice(&array(base.twist_body.angular));
    velocity.extend_from_slice(rates);
    velocity
}
fn same_mass(a: &AeroGeneralizedLoads, b: &AeroGeneralizedLoads) {
    for (x, y) in a
        .m_added_rr
        .iter()
        .flatten()
        .chain(bounded_cross(a))
        .chain(a.m_added_cc.iter().flatten())
        .zip(
            b.m_added_rr
                .iter()
                .flatten()
                .chain(bounded_cross(b))
                .chain(b.m_added_cc.iter().flatten()),
        )
    {
        close(*x, *y, 4e-13);
    }
}
fn bounded_cross(loads: &AeroGeneralizedLoads) -> impl Iterator<Item = &f64> {
    loads.m_added_rc.iter().flatten()
}

#[test]
fn zero_flow_recovers_original_blocks_and_comoving_fluid_has_no_disturbance() {
    let model = rigid_model(2.0);
    let strips = [LinkedStrip::new(0, plate([0.0; 3])).unwrap()];
    let base = FreeFloatingBaseState::new(
        Se3::identity(),
        Twist::new(vector([0.3, -0.4, 0.2]), vector([1.2, -0.4, 0.7])),
    );
    let dry = free_floating_terms(&model, base, &[], &[], || false).unwrap();
    let still = assemble_articulated_strips(1.18, &strips, &dry.kinematics).unwrap();
    let zero = assemble_articulated_strips_in_uniform_flow(
        1.18,
        &strips,
        &dry.kinematics,
        UniformFlow::default(),
    )
    .unwrap();
    assert_eq!(zero.loads(), &still);
    let acceleration = [0.2, -0.1, 0.4, -0.3, 0.2, 0.1];
    let reaction = zero.reaction(&acceleration).unwrap();
    let u = generalized_velocity(base, &[]);
    close(
        reaction.body_power_w,
        reaction
            .generalized_force
            .iter()
            .zip(u)
            .map(|(q, v)| q * v)
            .sum(),
        2e-13,
    );
    close(reaction.storage_rate_w + reaction.body_power_w, 0.0, 2e-13);
    assert_eq!(reaction.ambient_power_w, 0.0);

    let translating = FreeFloatingBaseState::new(
        Se3::identity(),
        Twist::new(Vec3::default(), vector([1.2, -0.4, 0.7])),
    );
    let dry = free_floating_terms(&model, translating, &[], &[], || false).unwrap();
    let old = assemble_articulated_strips(1.18, &strips, &dry.kinematics).unwrap();
    assert!(
        old.q_added_bias[4].abs() > 0.1,
        "still-water Munk torque distinguishes the models"
    );
    let comoving = assemble_articulated_strips_in_uniform_flow(
        1.18,
        &strips,
        &dry.kinematics,
        UniformFlow {
            velocity_world_mps: translating.twist_body.linear,
            acceleration_world_mps2: Vec3::default(),
        },
    )
    .unwrap();
    same_mass(&old, comoving.loads());
    assert_eq!(comoving.disturbance_energy_j(), 0.0);
    assert_eq!(comoving.generalized_fluid_momentum(), &[0.0; 6]);
    assert_eq!(comoving.loads().q_added_bias, vec![0.0; 6]);

    // Load callbacks may clone and combine blocks. Diagnostics remain tied to
    // their original state rather than accepting an accidentally edited block.
    let original = comoving.reaction(&acceleration).unwrap();
    let mut edited = comoving.loads().clone();
    edited.q_added_bias[2] = 1e6;
    assert_ne!(&edited, comoving.loads());
    assert_eq!(comoving.reaction(&acceleration).unwrap(), original);
    assert_eq!(
        comoving.reaction(&[0.0; 5]).unwrap_err().code,
        "uniform-flow-acceleration-invalid"
    );
    assert_eq!(
        comoving.reaction(&[f64::NAN; 6]).unwrap_err().code,
        "uniform-flow-acceleration-invalid"
    );
    assert_eq!(
        assemble_articulated_strips_in_uniform_flow(
            1.18,
            &strips,
            &dry.kinematics,
            UniformFlow {
                velocity_world_mps: vector([f64::NAN, 0.0, 0.0]),
                acceleration_world_mps2: Vec3::default(),
            },
        )
        .unwrap_err()
        .code,
        "uniform-flow-invalid"
    );
}

#[test]
fn pitching_plate_feels_tangential_current_before_it_has_normal_relative_speed() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(2.0, [0.0; 3]),
        ),
        Link::new(
            "pitching plate",
            Some(0),
            Se3::identity(),
            JointModel::revolute(vector([0.0, 1.0, 0.0]), None).unwrap(),
            inertia(1.0, [0.0; 3]),
        ),
    ])
    .unwrap();
    let geometry = plate([0.0; 3]);
    let ma = 1.18 * std::f64::consts::PI * geometry.chord_m.powi(2) * geometry.span_m / 4.0;
    let pitch_inertia = ma * geometry.chord_m.powi(2) / 32.0;
    let strips = [LinkedStrip::new(1, geometry).unwrap()];
    let pitch_rate = 1.6;
    let current = 2.3;
    let flow = flow_assembly(
        &model,
        FreeFloatingBaseState::stationary(Se3::identity()),
        &[0.0],
        &[pitch_rate],
        &strips,
        UniformFlow {
            velocity_world_mps: vector([current, 0.0, 0.0]),
            acceleration_world_mps2: Vec3::default(),
        },
    );
    let mut expected_bias = [0.0; 7];
    expected_bias[2] = ma * pitch_rate * current;
    close_slice(&flow.loads().q_added_bias, &expected_bias, 3e-13);
    close(
        flow.disturbance_energy_j(),
        0.5 * pitch_inertia * pitch_rate.powi(2),
        2e-13,
    );
    let mut expected_momentum = [0.0; 7];
    expected_momentum[4] = pitch_inertia * pitch_rate;
    expected_momentum[6] = pitch_inertia * pitch_rate;
    close_slice(flow.generalized_fluid_momentum(), &expected_momentum, 2e-13);

    let acceleration = [0.0, 0.0, 0.2, 0.0, 0.3, 0.0, -0.2];
    let reaction = flow.reaction(&acceleration).unwrap();
    expected_bias[2] -= ma * acceleration[2];
    expected_bias[4] = -pitch_inertia * (acceleration[4] + acceleration[6]);
    expected_bias[6] = expected_bias[4];
    close_slice(&reaction.generalized_force, &expected_bias, 3e-13);
    close(reaction.body_power_w, expected_bias[6] * pitch_rate, 2e-13);
    close(
        reaction.storage_rate_w,
        -expected_bias[6] * pitch_rate,
        2e-13,
    );
    close(reaction.ambient_power_w, 0.0, 2e-13);
}

fn serial_model() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        // The mounting frame is deliberately distinct from the floating base:
        // world reaction force must not mistake root-link axes for base axes.
        Link::new(
            "mounted root",
            None,
            pose([0.2, -0.1, 0.15], [0.1, -0.3, 0.2]),
            JointModel::FIXED,
            inertia(4.0, [0.1, -0.05, 0.02]),
        ),
        Link::new(
            "oblique hinge",
            Some(0),
            pose([0.1, -0.2, 0.15], [0.5, 0.2, -0.1]),
            JointModel::revolute(vector([0.3, 0.8, 0.4]), None).unwrap(),
            inertia(1.8, [0.25, 0.05, -0.1]),
        ),
        Link::new(
            "moving screw",
            Some(1),
            pose([-0.15, 0.1, 0.2], [0.4, -0.2, 0.3]),
            JointModel::helical(vector([-0.4, 0.3, 0.8]), 0.12, None).unwrap(),
            inertia(1.1, [0.2, -0.1, 0.15]),
        ),
    ])
    .unwrap()
}

// Independent span quadrature and section-pitch energy, never an assembled
// matrix or the production rank modes. Native FK supplies only physical poses
// and twists; the oracle applies the plate law directly at material points.
fn physical_energy(
    model: &ArticulatedModel,
    base: FreeFloatingBaseState,
    q: &[f64],
    rates: &[f64],
    geometry: &[Strip],
    current_world: Vec3,
) -> f64 {
    let poses = forward_kinematics(
        model,
        BaseState::prescribed(base.world_from_base, base.twist_body, Twist::zero()),
        q,
        rates,
    )
    .unwrap();
    let mut energy = 0.0;
    for (i, strip) in geometry.iter().enumerate() {
        let velocity = poses.body_twist[i];
        let current = poses.world_from_link[i]
            .rotation()
            .inverse()
            .rotate(current_world)
            .unwrap();
        let ma = 1.18 * std::f64::consts::PI * strip.chord_m.powi(2) * strip.span_m / 4.0;
        for sign in [-1.0, 1.0] {
            let point = vector(strip.position_m)
                + vector(strip.span_axis).scale(sign * strip.span_m / (2.0 * 3.0_f64.sqrt()));
            let normal_speed = (velocity.linear - current + velocity.angular.cross(point))
                .dot(vector(strip.normal));
            energy += 0.25 * ma * normal_speed.powi(2);
        }
        energy += 0.5 * ma * strip.chord_m.powi(2) / 32.0
            * velocity.angular.dot(vector(strip.span_axis)).powi(2);
    }
    energy
}

#[test]
fn rotating_serial_flow_preserves_galilean_and_world_frame_physics_and_storage_rate() {
    let model = serial_model();
    let base = FreeFloatingBaseState::new(
        pose([-0.2, 0.3, 0.1], [0.4, -0.2, 0.8]),
        Twist::new(vector([0.3, -0.25, 0.2]), vector([0.8, -0.3, 0.2])),
    );
    let q = [0.35, -0.25];
    let rates = [0.7, -0.45];
    let mut geometry = [
        plate([0.2, -0.3, 0.1]),
        plate([0.3, 0.08, -0.11]),
        plate([0.45, -0.1, 0.2]),
    ];
    geometry[1].normal = [0.6, 0.0, 0.8];
    let strips: Vec<_> = geometry
        .iter()
        .enumerate()
        .map(|(i, s)| LinkedStrip::new(i, s.clone()).unwrap())
        .collect();
    let flow = UniformFlow {
        velocity_world_mps: vector([1.1, -0.6, 0.9]),
        acceleration_world_mps2: vector([0.2, 0.1, -0.15]),
    };
    let dry = free_floating_terms(&model, base, &q, &rates, || false).unwrap();
    let assembly =
        assemble_articulated_strips_in_uniform_flow(1.18, &strips, &dry.kinematics, flow).unwrap();
    let (mut acceleration, joint) = solve_articulated_joint(&dry, assembly.loads()).unwrap();
    acceleration.extend(joint);
    let reaction = assembly.reaction(&acceleration).unwrap();
    let u = generalized_velocity(base, &rates);
    close(
        reaction.body_power_w,
        reaction
            .generalized_force
            .iter()
            .zip(&u)
            .map(|(f, v)| f * v)
            .sum(),
        3e-12,
    );
    let expected_force = base
        .world_from_base
        .rotation()
        .rotate(vector([
            reaction.generalized_force[0],
            reaction.generalized_force[1],
            reaction.generalized_force[2],
        ]))
        .unwrap();
    close_slice(&array(reaction.force_world), &array(expected_force), 3e-12);
    close(
        assembly.disturbance_energy_j(),
        physical_energy(&model, base, &q, &rates, &geometry, flow.velocity_world_mps),
        4e-13,
    );

    // Differentiate independently evaluated material-point energy along the
    // actual velocity/acceleration, including the current's explicit time law.
    let epsilon = 2e-6;
    let mut energies = [0.0; 2];
    for (k, dt) in [-epsilon, epsilon].into_iter().enumerate() {
        let shifted_base = FreeFloatingBaseState::new(
            base.world_from_base
                .compose(Se3::exp(base.twist_body.scale(dt)).unwrap())
                .unwrap(),
            base.twist_body.plus(
                Twist::new(
                    vector([acceleration[3], acceleration[4], acceleration[5]]),
                    vector([acceleration[0], acceleration[1], acceleration[2]]),
                )
                .scale(dt),
            ),
        );
        let shifted_q = [q[0] + dt * rates[0], q[1] + dt * rates[1]];
        let shifted_rates = [
            rates[0] + dt * acceleration[6],
            rates[1] + dt * acceleration[7],
        ];
        energies[k] = physical_energy(
            &model,
            shifted_base,
            &shifted_q,
            &shifted_rates,
            &geometry,
            flow.velocity_world_mps + flow.acceleration_world_mps2.scale(dt),
        );
    }
    close(
        reaction.storage_rate_w,
        (energies[1] - energies[0]) / (2.0 * epsilon),
        2e-8,
    );

    let boost = vector([1.4, -0.8, 0.6]);
    let boost_body = base
        .world_from_base
        .rotation()
        .inverse()
        .rotate(boost)
        .unwrap();
    let boosted_base = FreeFloatingBaseState::new(
        base.world_from_base,
        Twist::new(base.twist_body.angular, base.twist_body.linear + boost_body),
    );
    let boosted_dry = free_floating_terms(&model, boosted_base, &q, &rates, || false).unwrap();
    let boosted = assemble_articulated_strips_in_uniform_flow(
        1.18,
        &strips,
        &boosted_dry.kinematics,
        UniformFlow {
            velocity_world_mps: flow.velocity_world_mps + boost,
            ..flow
        },
    )
    .unwrap();
    same_mass(assembly.loads(), boosted.loads());
    close(
        boosted.disturbance_energy_j(),
        assembly.disturbance_energy_j(),
        5e-13,
    );
    close_slice(
        boosted.generalized_fluid_momentum(),
        assembly.generalized_fluid_momentum(),
        8e-13,
    );
    let (mut boosted_acceleration, joints) =
        solve_articulated_joint(&boosted_dry, boosted.loads()).unwrap();
    boosted_acceleration.extend(joints);
    let mut expected_acceleration = acceleration.clone();
    let correction = array(base.twist_body.angular.cross(boost_body));
    for i in 0..3 {
        expected_acceleration[i] -= correction[i];
    }
    close_slice(&boosted_acceleration, &expected_acceleration, 2e-11);
    let boosted_reaction = boosted.reaction(&boosted_acceleration).unwrap();
    close_slice(
        &boosted_reaction.generalized_force,
        &reaction.generalized_force,
        2e-11,
    );
    close(
        boosted_reaction.storage_rate_w,
        reaction.storage_rate_w,
        2e-11,
    );
    close(
        boosted_reaction.body_power_w,
        reaction.body_power_w + boost.dot(reaction.force_world),
        2e-11,
    );
    close(
        boosted_reaction.ambient_power_w,
        reaction.ambient_power_w + boost.dot(reaction.force_world),
        2e-11,
    );

    let scene = pose([0.4, -0.3, 0.2], [1.2, -0.7, 0.5]);
    let transformed_base = FreeFloatingBaseState::new(
        scene.compose(base.world_from_base).unwrap(),
        base.twist_body,
    );
    let transformed = flow_assembly(
        &model,
        transformed_base,
        &q,
        &rates,
        &strips,
        UniformFlow {
            velocity_world_mps: scene.rotation().rotate(flow.velocity_world_mps).unwrap(),
            acceleration_world_mps2: scene
                .rotation()
                .rotate(flow.acceleration_world_mps2)
                .unwrap(),
        },
    );
    same_mass(assembly.loads(), transformed.loads());
    close_slice(
        &transformed.loads().q_added_bias,
        &assembly.loads().q_added_bias,
        3e-12,
    );
    close(
        transformed.disturbance_energy_j(),
        assembly.disturbance_energy_j(),
        3e-12,
    );
    close_slice(
        transformed.generalized_fluid_momentum(),
        assembly.generalized_fluid_momentum(),
        3e-12,
    );
    let transformed_reaction = transformed.reaction(&acceleration).unwrap();
    close_slice(
        &array(transformed_reaction.force_world),
        &array(scene.rotation().rotate(reaction.force_world).unwrap()),
        3e-12,
    );
    close(
        transformed_reaction.body_power_w,
        reaction.body_power_w,
        3e-12,
    );
    close(
        transformed_reaction.ambient_power_w,
        reaction.ambient_power_w,
        3e-12,
    );
    close(
        transformed_reaction.storage_rate_w,
        reaction.storage_rate_w,
        3e-12,
    );
}

#[test]
fn accelerating_current_drives_actual_midpoint_with_exact_momentum_and_source_work() {
    let mass = 2.3;
    let model = rigid_model(mass);
    let rho = 1.2;
    let geometry = plate([0.0; 3]);
    let ma = rho * std::f64::consts::PI * geometry.chord_m.powi(2) * geometry.span_m / 4.0;
    let strips = [LinkedStrip::new(0, geometry).unwrap()];
    let fluid_acceleration = 1.7;
    // This is added mass alone: no displaced-volume pressure-gradient force.
    let body_acceleration = ma / (mass + ma) * fluid_acceleration;
    let mut state = ArticulatedState {
        motion: FreeFloatingMotionState {
            base: FreeFloatingBaseState::stationary(Se3::identity()),
            joint_position: vec![],
            joint_velocity: vec![],
        },
        memory: vec![],
    };
    let dt = 0.04;
    let mut source_work = 0.0;
    for tick in 0..12 {
        let t = tick as f64 * dt;
        let result = step(
            &model,
            &[],
            &state,
            &[],
            t,
            dt,
            ArticulatedStepConfig {
                max_iterations: 16,
                linear_tolerance_mps: 1e-12,
                angular_tolerance_rad_s: 1e-12,
                joint_velocity_tolerances: &[],
            },
            |time, _, map| {
                Ok(ArticulatedLoads {
                    fluid: assemble_articulated_strips_in_uniform_flow(
                        rho,
                        &strips,
                        map,
                        UniformFlow {
                            velocity_world_mps: vector([0.0, 0.0, fluid_acceleration * time]),
                            acceleration_world_mps2: vector([0.0, 0.0, fluid_acceleration]),
                        },
                    )?
                    .into_loads(),
                    world_wrenches: vec![Wrench::default()],
                    gravity_world: Vec3::default(),
                    joint_effort: vec![],
                })
            },
            || false,
        )
        .unwrap();
        state = result.state;
        let endpoint_time = t + dt;
        let current = fluid_acceleration * endpoint_time;
        let velocity = state.motion.base.twist_body.linear.z;
        close(velocity, body_acceleration * endpoint_time, 2e-12);
        close(
            state.motion.base.world_from_base.translation().z,
            0.5 * body_acceleration * endpoint_time.powi(2),
            2e-12,
        );
        let map = free_floating_terms(&model, state.motion.base, &[], &[], || false)
            .unwrap()
            .kinematics;
        let assembly = assemble_articulated_strips_in_uniform_flow(
            rho,
            &strips,
            &map,
            UniformFlow {
                velocity_world_mps: vector([0.0, 0.0, current]),
                acceleration_world_mps2: vector([0.0, 0.0, fluid_acceleration]),
            },
        )
        .unwrap();
        close(
            assembly.generalized_fluid_momentum()[2],
            ma * (velocity - current),
            2e-12,
        );
        close(
            mass * velocity + assembly.generalized_fluid_momentum()[2],
            0.0,
            2e-12,
        );
        let reaction = assembly
            .reaction(&[0.0, 0.0, body_acceleration, 0.0, 0.0, 0.0])
            .unwrap();
        close(reaction.force_world.z, mass * body_acceleration, 2e-12);
        close(
            reaction.body_power_w,
            mass * body_acceleration * velocity,
            2e-12,
        );
        close(
            reaction.ambient_power_w,
            current * mass * body_acceleration,
            2e-12,
        );
        close(
            reaction.storage_rate_w,
            ma * (velocity - current) * (body_acceleration - fluid_acceleration),
            2e-12,
        );
        source_work += dt * (fluid_acceleration * (t + 0.5 * dt)) * mass * body_acceleration;
        let total_energy = 0.5 * mass * velocity.powi(2) + assembly.disturbance_energy_j();
        close(
            total_energy,
            0.5 * mass * ma / (mass + ma) * current.powi(2),
            2e-12,
        );
        close(total_energy, source_work, 2e-12);
    }
}
