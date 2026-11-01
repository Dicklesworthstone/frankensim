//! Simultaneous sustained contact between a prismatic linkage and a free body.
//! The independent oracle eliminates the one physical contact constraint from
//! a three-coordinate linear midpoint KKT system. It never reads a production
//! mass/Delassus matrix or calls a production linear solver.

use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_flyer::addedmass::Strip;
use fs_flyer::addedmass::articulated::{LinkedStrip, assemble_articulated_strips};
use fs_flyer::partitioned::articulated::contact::systems::{ContactSystem, step_with_contacts};
use fs_flyer::partitioned::articulated::impact::friction::ArticulatedFrictionalContact;
use fs_flyer::partitioned::articulated::impact::{ContactCounterpart, LinkContactPoint};
use fs_flyer::partitioned::articulated::{
    ArticulatedLoads, ArticulatedState, ArticulatedStepConfig, JointImpedance,
};
use fs_ga::{Mat3, Se3, So3, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};

fn z(value: f64) -> Vec3 {
    Vec3::new(0., 0., value)
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * (1. + expected.abs()),
        "{actual:.16e} != {expected:.16e}",
    );
}
fn axial(actual: Vec3, expected: f64) {
    close(actual.x, 0.);
    close(actual.y, 0.);
    close(actual.z, expected);
}
fn inertia(mass: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [1., 0., 0., 0., 1., 0., 0., 0., 1.],
        },
    )
    .unwrap()
}
fn initial(height: f64, velocity: f64, joint: Option<(f64, f64)>) -> ArticulatedState {
    ArticulatedState {
        motion: FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(
                Se3::from_parts(So3::identity(), z(height)).unwrap(),
                Twist::new(Vec3::default(), z(velocity)),
            ),
            joint_position: joint.into_iter().map(|(q, _)| q).collect(),
            joint_velocity: joint.into_iter().map(|(_, rate)| rate).collect(),
        },
        memory: Vec::new(),
    }
}
fn strip(link: usize, span: f64, name: &'static str) -> LinkedStrip {
    LinkedStrip::new(
        link,
        Strip {
            name,
            chord_m: 1.,
            span_m: span,
            position_m: [0.; 3],
            normal: [0., 0., 1.],
            span_axis: [0., 1., 0.],
            control_coord: None,
            control_gain: 0.,
            control_angular_gain: [0.; 3],
        },
    )
    .unwrap()
}

#[derive(Clone, Copy)]
struct Case {
    // Carrier, child and independent other body's physical dry masses/spans.
    dry: [f64; 3],
    span: [f64; 3],
    density: f64,
    impedance: JointImpedance,
    drive: f64,
    gravity: f64,
    q0: f64,
    // Carrier velocity, relative joint rate and other body's velocity.
    u0: [f64; 3],
    dt: f64,
}
impl Case {
    fn masses(self) -> [f64; 3] {
        // Independent finite-strip normal inertia, rho*pi*c^2*span/4.
        core::array::from_fn(|i| {
            self.dry[i] + self.density * std::f64::consts::PI * self.span[i] / 4.
        })
    }
    fn predict(self) -> Prediction {
        let [a, l, b] = self.masses();
        let [mb, ml, mo] = self.dry;
        let [vb, qdot, vo] = self.u0;
        let h = 0.5 * self.dt;
        let d = h * self.impedance.damping + h * h * self.impedance.stiffness;
        let s = l + b;
        let spring = self.impedance.stiffness * (self.q0 - self.impedance.rest_position);
        // Physical generalized mass in [v_carrier, qdot, v_other] is
        // [[A+L,L,0],[L,L,0],[0,0,B]]. The midpoint effective mass adds
        // h*c+h^2*k only to the joint diagonal. Contact J=[1,1,-1].
        // Eliminating J*w=0 leaves the independent two-coordinate system
        // [[A+S,S],[S,S+d]]*w=[r0,r1], where S=L+B.
        let r0 = (a + l) * vb + l * qdot + b * vo - h * (mb + ml + mo) * self.gravity;
        let r1 = l * (vb + qdot) + b * vo + h * (self.drive - spring - (ml + mo) * self.gravity);
        let determinant = a * s + (a + s) * d;
        let wm_base = ((s + d) * r0 - s * r1) / determinant;
        let wm_joint = ((a + s) * r1 - s * r0) / determinant;
        let wm_other = wm_base + wm_joint;
        let midpoint = [wm_base, wm_joint, wm_other];
        let endpoint = core::array::from_fn(|i| 2. * midpoint[i] - self.u0[i]);
        // The other body's momentum balance fixes the FULL reaction impulse;
        // its gravity acts on dry mass, while acceleration uses dry+fluid mass.
        let impulse = 2. * b * (vo - wm_other) - self.dt * mo * self.gravity;
        Prediction {
            midpoint,
            endpoint,
            impulse,
            q1: self.q0 + self.dt * wm_joint,
        }
    }
}
struct Prediction {
    midpoint: [f64; 3],
    endpoint: [f64; 3],
    q1: f64,
    impulse: f64,
}
#[derive(Debug)]
struct AffinePathBound {
    initial_gap_m: f64,
    endpoint_gap_m: f64,
    minimum_gap_m: f64,
}

#[test]
fn moving_joint_and_other_tree_share_one_midpoint_reaction_and_complete_physical_work() {
    let dry = [2., 1., 3.];
    let carrier = ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(dry[0]),
        ),
        Link::new(
            "contacting slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(z(1.), None).unwrap(),
            inertia(dry[1]),
        ),
    ])
    .unwrap();
    let other = ArticulatedModel::new(vec![Link::new(
        "independent body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(dry[2]),
    )])
    .unwrap();
    let span = [2., 0.8, 4.];
    let carrier_strips = [
        strip(0, span[0], "carrier fluid"),
        strip(1, span[1], "slider fluid"),
    ];
    let other_strips = [strip(0, span[2], "other body fluid")];
    let gap = 1. / 16_777_216.; // Explicit positive binary gap; never projected.
    let joint_tolerances = [1e-13];
    let smooth_a = ArticulatedStepConfig {
        max_iterations: 32,
        linear_tolerance_mps: 1e-13,
        angular_tolerance_rad_s: 1e-13,
        joint_velocity_tolerances: &joint_tolerances,
    };
    let smooth_b = ArticulatedStepConfig {
        joint_velocity_tolerances: &[],
        ..smooth_a
    };
    let mut speeds_by_density = Vec::new();
    for density in [0., 1.2] {
        let case = Case {
            dry,
            span,
            density,
            impedance: JointImpedance {
                stiffness: 256.,
                damping: 8.,
                rest_position: 1. / 64.,
            },
            drive: 3.,
            gravity: 9.81,
            q0: 1. / 16.,
            u0: [-1. / 8., 3. / 16., 1. / 16.],
            dt: 1. / 32.,
        };
        let expected = case.predict();
        assert!(expected.impulse > 0.1, "fixture must remain compressive");
        let initial_a = initial(-case.q0, case.u0[0], Some((case.q0, case.u0[1])));
        let initial_b = initial(-gap, case.u0[2], None);
        let retained = [initial_a.clone(), initial_b.clone()];
        let impedance = [case.impedance];
        let mut order_results = Vec::new();
        for reverse in [false, true] {
            let mut systems = vec![
                ContactSystem {
                    model: &carrier,
                    state: &initial_a,
                    impedance: &impedance,
                    memory_inputs: &[],
                    smooth: smooth_a,
                },
                ContactSystem {
                    model: &other,
                    state: &initial_b,
                    impedance: &[],
                    memory_inputs: &[],
                    smooth: smooth_b,
                },
            ];
            if reverse {
                systems.reverse();
            }
            let (ia, ib) = if reverse { (1, 0) } else { (0, 1) };
            let contact = ArticulatedFrictionalContact {
                id: 29,
                point: LinkContactPoint {
                    system_index: ia,
                    link_index: 1,
                    point_body: Vec3::default(),
                },
                other: ContactCounterpart::Body(LinkContactPoint {
                    system_index: ib,
                    link_index: 0,
                    point_body: Vec3::default(),
                }),
                normal_world: z(1.),
                tangent_world: Vec3::new(1., 0., 0.),
                activation_tolerance: 2. * gap,
                friction_coefficient: 0.,
                normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0. },
            };
            let result = step_with_contacts(
                &systems,
                0.7,
                case.dt,
                &[contact],
                FrictionalContactConfig::default(),
                |_, stages| {
                    stages
                        .iter()
                        .enumerate()
                        .map(|(i, stage)| {
                            let strips: &[LinkedStrip] = if density == 0. {
                                &[]
                            } else if i == ia {
                                &carrier_strips
                            } else {
                                &other_strips
                            };
                            Ok(ArticulatedLoads {
                                fluid: assemble_articulated_strips(
                                    if density == 0. { 1. } else { density },
                                    strips,
                                    stage.kinematics,
                                )?,
                                world_wrenches: vec![
                                    Wrench::new(Vec3::default(), Vec3::default());
                                    if i == ia { 2 } else { 1 }
                                ],
                                gravity_world: z(-case.gravity),
                                joint_effort: if i == ia {
                                    vec![case.drive]
                                } else {
                                    Vec::new()
                                },
                            })
                        })
                        .collect()
                },
                |path| {
                    // This is a proof for this declared affine motion class,
                    // not a finite-sample CCD claim. Every base has exactly zero
                    // angular midpoint twist, and the sole joint is a fixed-axis
                    // prismatic translation. Thus each material point and their
                    // separating gap are affine over the entire accepted path.
                    for i in 0..path.systems.len() {
                        assert_eq!(
                            path.systems[i].state.motion.base.twist_body.angular,
                            Vec3::default()
                        );
                        assert_eq!(
                            path.midpoints[i].motion.base.twist_body.angular,
                            Vec3::default()
                        );
                        assert_eq!(
                            path.endpoints[i].motion.base.twist_body.angular,
                            Vec3::default()
                        );
                        assert_eq!(path.midpoints[i].motion.base.twist_body.linear.x, 0.);
                        assert_eq!(path.midpoints[i].motion.base.twist_body.linear.y, 0.);
                    }
                    let height_a = |s: &ArticulatedState| {
                        s.motion.base.world_from_base.translation().z + s.motion.joint_position[0]
                    };
                    let g0 = height_a(path.systems[ia].state)
                        - path.systems[ib]
                            .state
                            .motion
                            .base
                            .world_from_base
                            .translation()
                            .z;
                    let g1 = height_a(path.endpoints[ia])
                        - path.endpoints[ib]
                            .motion
                            .base
                            .world_from_base
                            .translation()
                            .z;
                    assert!(g0 >= 0. && g1 >= 0., "actual endpoint gaps: {g0:e}, {g1:e}");
                    // This sample only checks that the public path query follows
                    // the affine interpolant used by the endpoint bound above.
                    let fraction = 0.375;
                    let sampled = path.sample(fraction, || false)?;
                    let actual = sampled[ia].world_from_link[1].translation().z
                        - sampled[ib].world_from_link[0].translation().z;
                    close(actual, (1. - fraction) * g0 + fraction * g1);
                    Ok(AffinePathBound {
                        initial_gap_m: g0,
                        endpoint_gap_m: g1,
                        minimum_gap_m: g0.min(g1),
                    })
                },
                || false,
            )
            .unwrap();
            let end_a = &result.steps[ia].state;
            let end_b = &result.steps[ib].state;
            let actual_u = [
                end_a.motion.base.twist_body.linear.z,
                end_a.motion.joint_velocity[0],
                end_b.motion.base.twist_body.linear.z,
            ];
            for i in 0..3 {
                close(actual_u[i], expected.endpoint[i]);
            }
            axial(end_a.motion.base.twist_body.linear, expected.endpoint[0]);
            axial(end_b.motion.base.twist_body.linear, expected.endpoint[2]);
            assert_eq!(end_a.motion.base.twist_body.angular, Vec3::default());
            assert_eq!(end_b.motion.base.twist_body.angular, Vec3::default());
            close(end_a.motion.joint_position[0], expected.q1);
            close(
                end_a.motion.base.world_from_base.translation().z,
                -case.q0 + case.dt * expected.midpoint[0],
            );
            close(
                end_b.motion.base.world_from_base.translation().z,
                -gap + case.dt * expected.midpoint[2],
            );
            close(result.reaction_solve.impulses[0][0], expected.impulse);
            close(result.reaction_solve.impulses[0][1], 0.);
            close(result.reaction_solve.impulses[0][2], 0.);
            assert_eq!(result.contact_ids, vec![29]);
            close(result.initial_normal_speeds[0], 0.);
            close(result.midpoint_relative_velocities[0][0], 0.);
            close(result.endpoint_relative_velocities[0][0], 0.);
            close(result.midpoint_velocity_residual_mps, 0.);
            close(result.reaction_work_j, 0.);
            close(result.normal_work_j, 0.);
            close(result.friction_work_j, 0.);
            close(result.support_work_j, 0.);
            assert_eq!(result.external_linear_impulse_world, Vec3::default());
            assert_eq!(result.external_angular_impulse_world, Vec3::default());
            assert!(result.initial_gaps[0] >= 0. && result.endpoint_gaps[0] >= 0.);
            assert!(result.path_admission.minimum_gap_m > 0.);
            close(result.initial_gaps[0], result.path_admission.initial_gap_m);
            close(
                result.endpoint_gaps[0],
                result.path_admission.endpoint_gap_m,
            );
            close(result.endpoint_gaps[0], gap);

            let [a, l, b] = case.masses();
            let midpoint: [f64; 3] = core::array::from_fn(|i| 0.5 * (case.u0[i] + actual_u[i]));
            let child0 = case.u0[0] + case.u0[1];
            let child1 = actual_u[0] + actual_u[1];
            let q1 = end_a.motion.joint_position[0];
            let qm = 0.5 * (case.q0 + q1);
            let joint_force = case.drive
                - case.impedance.stiffness * (qm - case.impedance.rest_position)
                - case.impedance.damping * midpoint[1];
            let impulse = result.reaction_solve.impulses[0][0];
            // Three separate Newton balances check the carrier's equal/opposite
            // actuator force and transmission of the other body's fluid inertia.
            close(
                a * (actual_u[0] - case.u0[0]),
                case.dt * (-case.dry[0] * case.gravity - joint_force),
            );
            close(
                l * (child1 - child0),
                case.dt * (-case.dry[1] * case.gravity + joint_force) + impulse,
            );
            close(
                b * (actual_u[2] - case.u0[2]),
                -case.dt * case.dry[2] * case.gravity - impulse,
            );
            let momentum_change = a * (actual_u[0] - case.u0[0])
                + l * (child1 - child0)
                + b * (actual_u[2] - case.u0[2]);
            close(
                momentum_change,
                -case.dt * case.dry.iter().sum::<f64>() * case.gravity,
            );

            // Independent physical storage/work: added mass stores kinetic
            // energy but receives no invented weight. Internal contact work is
            // zero; spring storage, dry gravity, damping and actuator work close.
            let kinetic_change = 0.5 * a * (actual_u[0].powi(2) - case.u0[0].powi(2))
                + 0.5 * l * (child1.powi(2) - child0.powi(2))
                + 0.5 * b * (actual_u[2].powi(2) - case.u0[2].powi(2));
            let spring_change = 0.5
                * case.impedance.stiffness
                * ((q1 - case.impedance.rest_position).powi(2)
                    - (case.q0 - case.impedance.rest_position).powi(2));
            let dz_a = end_a.motion.base.world_from_base.translation().z
                - initial_a.motion.base.world_from_base.translation().z;
            let dz_b = end_b.motion.base.world_from_base.translation().z
                - initial_b.motion.base.world_from_base.translation().z;
            let gravity_change = case.gravity
                * (case.dry[0] * dz_a + case.dry[1] * (dz_a + q1 - case.q0) + case.dry[2] * dz_b);
            let damping_loss = case.dt * case.impedance.damping * midpoint[1].powi(2);
            let actuator_work = case.drive * (q1 - case.q0);
            assert!(damping_loss > 0. && actuator_work > 0.);
            close(
                kinetic_change + spring_change + gravity_change + damping_loss - actuator_work,
                result.reaction_work_j,
            );
            assert_eq!(initial_a, retained[0]);
            assert_eq!(initial_b, retained[1]);
            order_results.push((actual_u, impulse));
        }
        for i in 0..3 {
            close(order_results[0].0[i], order_results[1].0[i]);
        }
        close(order_results[0].1, order_results[1].1);
        speeds_by_density.push(order_results[0].0);
    }
    assert!(
        speeds_by_density[1][2] > speeds_by_density[0][2] + 0.1,
        "the other tree's complete fluid inertia must materially change shared recoil"
    );
}
