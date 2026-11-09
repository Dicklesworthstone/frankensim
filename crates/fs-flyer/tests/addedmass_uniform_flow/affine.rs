//! G0/G1/G3: affine strip flow, variational forces and actual midpoint motion.
use super::*;
use fs_flyer::addedmass::articulated::{
    AffineFlow, AffineFlowAssembly, assemble_articulated_strips_in_affine_flow as assemble_affine,
};

fn sample(
    model: &ArticulatedModel,
    base: FreeFloatingBaseState,
    q: &[f64],
    rates: &[f64],
    strips: &[LinkedStrip],
    flow: AffineFlow,
) -> AffineFlowAssembly {
    let dry = free_floating_terms(model, base, q, rates, || false).unwrap();
    assemble_affine(1.18, strips, &dry.kinematics, flow).unwrap()
}

fn shifted_flow(flow: AffineFlow, dt: f64) -> AffineFlow {
    AffineFlow {
        velocity_world_mps: flow.velocity_world_mps
            + flow.partial_time_velocity_world_mps2.scale(dt),
        gradient_world_per_s: Mat3 {
            m: std::array::from_fn(|i| {
                flow.gradient_world_per_s.m[i] + dt * flow.partial_time_gradient_world_per_s2.m[i]
            }),
        },
        ..flow
    }
}

// Independent material-point normal velocities at center and plate edges.
// No production mode, bias, or gradient-variation implementation is consumed.
fn energy(
    model: &ArticulatedModel,
    base: FreeFloatingBaseState,
    q: &[f64],
    rates: &[f64],
    geometry: &[Strip],
    flow: AffineFlow,
) -> f64 {
    let fk = forward_kinematics(
        model,
        BaseState::prescribed(base.world_from_base, base.twist_body, Twist::zero()),
        q,
        rates,
    )
    .unwrap();
    let mut total = 0.0;
    for (i, strip) in geometry.iter().enumerate() {
        let pose = fk.world_from_link[i];
        let v = fk.body_twist[i];
        let normal = pose.rotation().rotate(vector(strip.normal)).unwrap();
        let normal_speed = |point: Vec3| {
            let world = pose.transform_point(point).unwrap();
            let point_velocity = pose
                .rotation()
                .rotate(v.linear + v.angular.cross(point))
                .unwrap();
            let air = flow.velocity_world_mps
                + flow
                    .gradient_world_per_s
                    .apply(world - flow.reference_position_world_m);
            normal.dot(point_velocity - air)
        };
        let center = vector(strip.position_m);
        let span = vector(strip.span_axis);
        let chord = span.cross(vector(strip.normal));
        let chord_half = chord.scale(0.5 * strip.chord_m);
        let span_half = span.scale(0.5 * strip.span_m);
        let slope_chord =
            (normal_speed(center + chord_half) - normal_speed(center - chord_half)) / strip.chord_m;
        let slope_span =
            (normal_speed(center + span_half) - normal_speed(center - span_half)) / strip.span_m;
        let ma = 1.18 * std::f64::consts::PI * strip.chord_m.powi(2) * strip.span_m / 4.0;
        total += 0.5
            * ma
            * (normal_speed(center).powi(2)
                + strip.chord_m.powi(2) / 32.0 * slope_chord.powi(2)
                + strip.span_m.powi(2) / 12.0 * slope_span.powi(2));
    }
    total
}

#[test]
fn uniform_and_still_limits_retain_exact_original_blocks_and_receipts() {
    let model = rigid_model(2.0);
    let strips = [LinkedStrip::new(0, plate([0.2, -0.1, 0.3])).unwrap()];
    let base = FreeFloatingBaseState::new(
        pose([0.2, -0.3, 0.1], [1.0, 0.3, -0.5]),
        Twist::new(vector([0.3, 0.2, -0.4]), vector([1.2, -0.7, 0.8])),
    );
    let dry = free_floating_terms(&model, base, &[], &[], || false).unwrap();
    for uniform in [
        UniformFlow::default(),
        UniformFlow {
            velocity_world_mps: vector([1.1, -0.4, 0.9]),
            acceleration_world_mps2: vector([0.2, 0.3, -0.1]),
        },
    ] {
        let old =
            assemble_articulated_strips_in_uniform_flow(1.18, &strips, &dry.kinematics, uniform)
                .unwrap();
        let new = assemble_affine(
            1.18,
            &strips,
            &dry.kinematics,
            AffineFlow {
                reference_position_world_m: vector([3.0, -1.0, 0.5]),
                velocity_world_mps: uniform.velocity_world_mps,
                partial_time_velocity_world_mps2: uniform.acceleration_world_mps2,
                ..AffineFlow::default()
            },
        )
        .unwrap();
        assert_eq!(old.loads(), new.loads());
        assert_eq!(old.disturbance_energy_j(), new.disturbance_energy_j());
        assert_eq!(
            old.generalized_fluid_momentum(),
            new.generalized_fluid_momentum()
        );
        let acceleration = [0.3, -0.2, 0.7, 0.1, -0.5, 0.4];
        assert_eq!(
            old.reaction(&acceleration).unwrap(),
            new.reaction(&acceleration).unwrap()
        );
    }
}

fn potential_flow(k: f64) -> AffineFlow {
    AffineFlow {
        gradient_world_per_s: Mat3 {
            m: [-0.5 * k, 0.0, 0.0, 0.0, -0.5 * k, 0.0, 0.0, 0.0, k],
        },
        ..AffineFlow::default()
    }
}

#[test]
fn potential_strain_uses_fluid_material_acceleration_not_body_advection_alone() {
    let mass = 2.1;
    let model = rigid_model(mass);
    let geometry = plate([0.0; 3]);
    let ma = 1.18 * std::f64::consts::PI * geometry.chord_m.powi(2) * geometry.span_m / 4.0;
    let strips = [LinkedStrip::new(0, geometry).unwrap()];
    let flow = AffineFlow {
        velocity_world_mps: vector([0.0, 0.0, 0.7]),
        partial_time_velocity_world_mps2: vector([0.0, 0.0, 0.4]),
        ..potential_flow(1.5)
    };
    let fluid_acceleration = 0.4 + 1.5 * (0.7 + 1.5 * 2.0);
    close(
        flow.material_acceleration_world_at(vector([0.0, 0.0, 2.0]))
            .unwrap()
            .z,
        fluid_acceleration,
        1e-14,
    );
    for speed in [0.0, 4.0, -2.0] {
        let base = FreeFloatingBaseState::new(
            pose([0.0; 3], [0.0, 0.0, 2.0]),
            Twist::new(Vec3::default(), vector([0.0, 0.0, speed])),
        );
        let dry = free_floating_terms(&model, base, &[], &[], || false).unwrap();
        let assembly = assemble_affine(1.18, &strips, &dry.kinematics, flow).unwrap();
        close(
            assembly.loads().q_added_bias[2],
            ma * fluid_acceleration,
            3e-13,
        );
        let (acceleration, _) = solve_articulated_joint(&dry, assembly.loads()).unwrap();
        close(
            acceleration[2],
            ma / (mass + ma) * fluid_acceleration,
            4e-13,
        );
        let reaction = assembly.reaction(&acceleration).unwrap();
        close(
            reaction.generalized_force[2],
            ma * (fluid_acceleration - acceleration[2]),
            4e-13,
        );
    }
}

#[test]
fn rigidly_comoving_rotating_air_has_no_finite_strip_disturbance() {
    let omega = vector([0.3, -0.5, 0.7]);
    let gradient = Mat3 {
        m: [
            0.0, -omega.z, omega.y, omega.z, 0.0, -omega.x, -omega.y, omega.x, 0.0,
        ],
    };
    let position = vector([1.0, -0.7, 0.4]);
    let world_pose = pose([0.2, -0.4, 0.3], array(position));
    let inverse = world_pose.rotation().inverse();
    let base = FreeFloatingBaseState::new(
        world_pose,
        Twist::new(
            inverse.rotate(omega).unwrap(),
            inverse.rotate(gradient.apply(position)).unwrap(),
        ),
    );
    let model = rigid_model(2.0);
    let strips = [LinkedStrip::new(0, plate([0.4, -0.2, 0.3])).unwrap()];
    let assembly = sample(
        &model,
        base,
        &[],
        &[],
        &strips,
        AffineFlow {
            gradient_world_per_s: gradient,
            ..AffineFlow::default()
        },
    );
    assert!(assembly.disturbance_energy_j() < 1e-27);
    close_slice(assembly.generalized_fluid_momentum(), &[0.0; 6], 2e-13);
    close_slice(
        &assembly.reaction(&[0.0; 6]).unwrap().generalized_force,
        &[0.0; 6],
        4e-13,
    );
}

#[test]
fn serial_links_match_independent_energy_variations_in_fixed_and_changing_gradients() {
    let model = serial_model();
    let base = FreeFloatingBaseState::new(
        pose([0.2, -0.3, 0.1], [0.8, -0.4, 0.6]),
        Twist::new(vector([0.3, -0.2, 0.4]), vector([0.7, -0.5, 1.1])),
    );
    let q = [0.37, -0.21];
    let rates = [-0.6, 0.9];
    let geometry = [
        plate([0.1, -0.2, 0.3]),
        plate([0.3, 0.1, -0.2]),
        plate([-0.2, 0.3, 0.1]),
    ];
    let strips: Vec<_> = geometry
        .iter()
        .enumerate()
        .map(|(i, p)| LinkedStrip::new(i, p.clone()).unwrap())
        .collect();
    for changing in [false, true] {
        let flow = AffineFlow {
            reference_position_world_m: vector([0.2, -0.3, 0.1]),
            velocity_world_mps: vector([1.2, -0.5, 0.7]),
            partial_time_velocity_world_mps2: if changing {
                vector([0.3, 0.2, -0.15])
            } else {
                Vec3::default()
            },
            gradient_world_per_s: Mat3 {
                m: [0.2, 0.3, -0.1, -0.15, -0.35, 0.22, 0.4, -0.25, 0.15],
            },
            partial_time_gradient_world_per_s2: if changing {
                Mat3 {
                    m: [0.05, -0.03, 0.02, 0.01, -0.02, -0.04, 0.03, 0.01, -0.03],
                }
            } else {
                Mat3::zero()
            },
        };
        let dry = free_floating_terms(&model, base, &q, &rates, || false).unwrap();
        let assembly = assemble_affine(1.18, &strips, &dry.kinematics, flow).unwrap();
        let (mut acceleration, joints) = solve_articulated_joint(&dry, assembly.loads()).unwrap();
        acceleration.extend(joints);
        let reaction = assembly.reaction(&acceleration).unwrap();
        close(
            assembly.disturbance_energy_j(),
            energy(&model, base, &q, &rates, &geometry, flow),
            2e-12,
        );
        let velocity = generalized_velocity(base, &rates);
        close(
            reaction.body_power_w,
            reaction
                .generalized_force
                .iter()
                .zip(&velocity)
                .map(|(f, v)| f * v)
                .sum(),
            3e-12,
        );
        let epsilon = 2e-6;
        let mut energies = [0.0; 2];
        let mut momenta = [Vec::new(), Vec::new()];
        for (side, dt) in [-epsilon, epsilon].into_iter().enumerate() {
            let shifted = FreeFloatingBaseState::new(
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
            let next_q = [q[0] + dt * rates[0], q[1] + dt * rates[1]];
            let next_rates = [
                rates[0] + dt * acceleration[6],
                rates[1] + dt * acceleration[7],
            ];
            energies[side] = energy(
                &model,
                shifted,
                &next_q,
                &next_rates,
                &geometry,
                shifted_flow(flow, dt),
            );
            momenta[side] = sample(
                &model,
                shifted,
                &next_q,
                &next_rates,
                &strips,
                shifted_flow(flow, dt),
            )
            .generalized_fluid_momentum()
            .to_vec();
        }
        let storage = (energies[1] - energies[0]) / (2.0 * epsilon);
        close(reaction.storage_rate_w, storage, 2e-7);
        close(
            storage + reaction.body_power_w,
            reaction.ambient_power_w,
            2e-7,
        );
        let p = assembly.generalized_fluid_momentum();
        let base_momentum = Wrench::new(vector([p[3], p[4], p[5]]), vector([p[0], p[1], p[2]]));
        for column in 0..8 {
            let mut eta = [0.0; 6];
            if column < 6 {
                eta[column] = 1.0;
            }
            let eta = Twist::new(
                vector([eta[3], eta[4], eta[5]]),
                vector([eta[0], eta[1], eta[2]]),
            );
            let mut variation = [0.0; 2];
            for (side, amount) in [-epsilon, epsilon].into_iter().enumerate() {
                let mut perturbed_base = base;
                let mut perturbed_q = q;
                if column < 6 {
                    perturbed_base.world_from_base = base
                        .world_from_base
                        .compose(Se3::exp(eta.scale(amount)).unwrap())
                        .unwrap();
                } else {
                    perturbed_q[column - 6] += amount;
                }
                variation[side] = energy(
                    &model,
                    perturbed_base,
                    &perturbed_q,
                    &rates,
                    &geometry,
                    flow,
                );
            }
            let expected = -(momenta[1][column] - momenta[0][column]) / (2.0 * epsilon)
                + base_momentum.pairing(base.twist_body.bracket(eta))
                + (variation[1] - variation[0]) / (2.0 * epsilon);
            close(reaction.generalized_force[column], expected, 3e-7);
        }
    }
}

#[test]
fn potential_gradient_drives_native_midpoint_with_second_order_motion_and_energy() {
    let mass = 2.0;
    let model = rigid_model(mass);
    let geometry = plate([0.0; 3]);
    let ma = 1.18 * std::f64::consts::PI * geometry.chord_m.powi(2) * geometry.span_m / 4.0;
    let strips = [LinkedStrip::new(0, geometry).unwrap()];
    let k: f64 = 1.2;
    let flow = potential_flow(k);
    let rate = (ma / (mass + ma)).sqrt() * k;
    let time = 0.8;
    let exact = 1.1 * (rate * time).cosh() + 0.4 / rate * (rate * time).sinh();
    let mut errors = Vec::new();
    for steps in [16, 32, 64] {
        let dt = time / f64::from(steps);
        let mut state = ArticulatedState {
            motion: FreeFloatingMotionState {
                base: FreeFloatingBaseState::new(
                    pose([0.0; 3], [0.0, 0.0, 1.1]),
                    Twist::new(Vec3::default(), vector([0.0, 0.0, 0.4])),
                ),
                joint_position: vec![],
                joint_velocity: vec![],
            },
            memory: vec![],
        };
        let invariant = 0.5 * (mass + ma) * 0.4_f64.powi(2) - 0.5 * ma * k * k * 1.1_f64.powi(2);
        for tick in 0..steps {
            state = step(
                &model,
                &[],
                &state,
                &[],
                f64::from(tick) * dt,
                dt,
                ArticulatedStepConfig {
                    max_iterations: 24,
                    linear_tolerance_mps: 1e-12,
                    angular_tolerance_rad_s: 1e-12,
                    joint_velocity_tolerances: &[],
                },
                |_, _, map| {
                    Ok(ArticulatedLoads {
                        fluid: assemble_affine(1.18, &strips, map, flow)?.into_loads(),
                        world_wrenches: vec![Wrench::default()],
                        gravity_world: Vec3::default(),
                        joint_effort: vec![],
                    })
                },
                || false,
            )
            .unwrap()
            .state;
        }
        let z = state.motion.base.world_from_base.translation().z;
        let v = state.motion.base.twist_body.linear.z;
        close(
            0.5 * (mass + ma) * v * v - 0.5 * ma * k * k * z * z,
            invariant,
            2e-10,
        );
        errors.push((z - exact).abs());
    }
    assert!(
        errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8,
        "{errors:?}"
    );
}

#[test]
fn malformed_affine_flow_refuses_without_mutating_native_state() {
    let model = rigid_model(2.0);
    let base = FreeFloatingBaseState::stationary(Se3::identity());
    let dry = free_floating_terms(&model, base, &[], &[], || false).unwrap();
    let before = dry.clone();
    let strips = [LinkedStrip::new(0, plate([0.0; 3])).unwrap()];
    let mut flow = potential_flow(0.2);
    flow.partial_time_gradient_world_per_s2.m[7] = f64::NAN;
    assert_eq!(
        assemble_affine(1.18, &strips, &dry.kinematics, flow)
            .unwrap_err()
            .code,
        "affine-flow-invalid"
    );
    assert_eq!(dry, before);
    assert!(
        potential_flow(1.0)
            .material_acceleration_world_at(vector([0.0, f64::NAN, 0.0]))
            .is_err()
    );
}
