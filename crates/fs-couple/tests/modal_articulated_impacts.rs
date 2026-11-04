//! Original mixed articulated inertia and support sources coupled to modal points.
use fs_contact::frictional_multicontact::{FrictionalContactDeterminacy, NormalVelocityLaw};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::plate::impact::radiation::{Model, Pole};
use fs_couple::render::schedule::force::coupled::impulse::impact::articulated::{
    ArticulatedModalImpactConfig, ArticulatedModalImpactContact, ArticulatedModalImpactError,
    ModalImpactPoint, PreparedArticulatedModalImpacts,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalCouplingConfig, ModalCouplingError,
};
use fs_exec::CancelGate;
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_math::c64::C64;
use fs_mbd::Vec3 as WorldVec3;
use fs_mbd::articulated::contact_geometry::LinkContactPoint;
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::motion::prescribed::PrescribedBaseMotionState;
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia, apply_free_floating_impulses, prescribed_base_impulse_delta,
};
use fs_mbd::articulated_contacts::mixed::MixedArticulatedImpulseResponse;
use fs_mbd::articulated_contacts::trajectory::mixed::MixedArticulatedTrajectoryState as State;

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 32,
        max_connections: 16,
        max_setup_terms: 1_000_000,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e5,
        maximum_abs_connection_force_n: 1e5,
        solve_relative_tolerance: 1e-10,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-9,
    }
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 3e-9 + 3e-9 * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn vector(a: Vec3, b: Vec3) {
    close(a.x, b.x);
    close(a.y, b.y);
    close(a.z, b.z);
}
fn world(value: Vec3) -> WorldVec3 {
    WorldVec3::new(value.x, value.y, value.z)
}
fn model(rate: u32, velocity: [f64; 3], budget: ModalAcousticTimeBudget) -> ModalAcousticTimeModel {
    let mut result = ModalAcousticTimeModel::try_new_with_free_coordinates(
        rate,
        vec![
            ModalAcousticMode {
                angular_frequency_rad_s: 0.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO
            };
            3
        ],
        budget,
    )
    .unwrap();
    result
        .restore_states(&velocity.map(|velocity| ModalAcousticState {
            displacement_m_sqrt_kg: 0.0,
            velocity_m_sqrt_kg_per_s: velocity,
        }))
        .unwrap();
    result
}
fn network(mass: f64, velocity: [f64; 3]) -> CoupledModalSystem {
    CoupledModalSystem::new(
        vec![model(
            2000,
            velocity.map(|v| v * mass.sqrt()),
            ModalAcousticTimeBudget::audible_reference(),
        )],
        vec![],
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn modal_states(system: &CoupledModalSystem) -> Vec<ModalAcousticState> {
    system
        .components()
        .iter()
        .flat_map(|model| model.states().iter().copied())
        .collect()
}
fn modal_point(mass: f64, reference: Vec3) -> ModalImpactPoint {
    let a = mass.sqrt().recip();
    ModalImpactPoint {
        component: 0,
        reference_position_world_m: world(reference),
        shapes_world: vec![
            WorldVec3::new(a, 0., 0.),
            WorldVec3::new(0., a, 0.),
            WorldVec3::new(0., 0., a),
        ],
    }
}
fn inertia(mass: f64, center: Vec3, polar: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        center,
        Mat3 {
            m: [polar, 0., 0., 0., polar, 0., 0., 0., polar],
        },
    )
    .unwrap()
}
fn single(joint: JointModel, mass: f64, center: Vec3, polar: f64) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        joint,
        inertia(mass, center, polar),
    )])
    .unwrap()
}
fn free(base: FreeFloatingBaseState, joint_position: Vec<f64>, joint_velocity: Vec<f64>) -> State {
    State::Free(FreeFloatingMotionState {
        base,
        joint_position,
        joint_velocity,
    })
}
fn prescribed(base: BaseState, joint_position: Vec<f64>, joint_velocity: Vec<f64>) -> State {
    State::Prescribed(PrescribedBaseMotionState {
        time_s: 3.5,
        base,
        joint_position,
        joint_velocity,
    })
}
fn contact(
    id: u64,
    system: usize,
    link: usize,
    point_body: Vec3,
    modal_point: ModalImpactPoint,
    mu: f64,
    restitution: f64,
) -> ArticulatedModalImpactContact {
    ArticulatedModalImpactContact {
        id,
        point: LinkContactPoint {
            system_index: system,
            link_index: link,
            point_body,
        },
        modal_point,
        normal_world: WorldVec3::new(1., 0., 0.),
        tangent_world: WorldVec3::new(0., 1., 0.),
        activation_tolerance_m: 0.0,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution {
            coefficient: restitution,
        },
    }
}
fn prepare(
    system: &CoupledModalSystem,
    contacts: Vec<ArticulatedModalImpactContact>,
    gate: &CancelGate,
) -> PreparedArticulatedModalImpacts {
    PreparedArticulatedModalImpacts::new(
        system,
        contacts,
        ArticulatedModalImpactConfig::default(),
        gate,
    )
    .unwrap()
}

#[test]
fn free_chain_retains_joint_base_and_cross_axis_coupling_in_a_rotated_world() {
    let gate = CancelGate::new_clock_free();
    let mechanism = ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(4., Vec3::default(), 1.),
        ),
        Link::new(
            "slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(0.6, 0.8, 0.), None).unwrap(),
            inertia(1., Vec3::default(), 1.),
        ),
    ])
    .unwrap();
    for (angle, origin) in [(0.0, Vec3::default()), (0.37, Vec3::new(3., -2., 1.))] {
        let rotation = So3::exp(So3Tangent::new(Vec3::new(0., 0., angle))).unwrap();
        let pose = Se3::from_parts(rotation, origin).unwrap();
        let incoming = Vec3::new(-0.892, -0.081, 0.);
        let mut states = [free(
            FreeFloatingBaseState::new(pose, Twist::new(Vec3::default(), incoming)),
            vec![0.],
            vec![0.],
        )];
        let mut system = network(2., [0.; 3]);
        let mut description = contact(41, 0, 1, Vec3::default(), modal_point(2., origin), 0.5, 0.);
        description.normal_world = world(rotation.rotate(Vec3::new(1., 0., 0.)).unwrap());
        description.tangent_world = world(rotation.rotate(Vec3::new(0., 1., 0.)).unwrap());
        for shape in &mut description.modal_point.shapes_world {
            *shape = world(
                rotation
                    .rotate(Vec3::new(shape.x, shape.y, shape.z))
                    .unwrap(),
            );
        }
        let mut prepared = prepare(&system, vec![description], &gate);
        let receipt = system
            .apply_articulated_frictional_impacts(&mut states, &[&mechanism], &mut prepared, &gate)
            .unwrap();
        // Original child mobility 0.2 I + 0.8 a a^T; modal mass adds 0.5 I.
        for (actual, expected) in receipt
            .contact_mobility_kg_inv
            .iter()
            .zip([0.988, 0.384, 0., 0.384, 1.212, 0., 0., 0., 0.7])
        {
            close(*actual, expected);
        }
        close(receipt.contacts.impulses[0][0], 1.);
        close(receipt.contacts.impulses[0][1], -0.25);
        close(receipt.contact_work_j, -3487. / 8000.);
        close(receipt.kinetic_energy_change_j, -3487. / 8000.);
        close(receipt.support_work_j, 0.);
        close(receipt.energy_residual_j, 0.);
        vector(receipt.linear_momentum_residual_world, Vec3::default());
        vector(receipt.angular_momentum_residual_world, Vec3::default());
        let State::Free(after) = &states[0] else {
            panic!("free boundary changed")
        };
        assert_eq!(after.base.world_from_base, pose);
        assert_eq!(after.joint_position, vec![0.]);
        vector(
            after.base.twist_body.linear,
            incoming + Vec3::new(0.152, -0.114, 0.),
        );
        close(after.joint_velocity[0], 0.4);
        close(
            modal_states(&system)[0].velocity_m_sqrt_kg_per_s / 2f64.sqrt(),
            -0.5,
        );
        close(
            modal_states(&system)[1].velocity_m_sqrt_kg_per_s / 2f64.sqrt(),
            0.125,
        );
        assert_eq!(system.samples_rendered(), 0);
    }
}

#[test]
fn driven_offset_com_hinge_transmits_native_support_impulse_and_positive_source_work() {
    let gate = CancelGate::new_clock_free();
    let mechanism = single(
        JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
        2.,
        Vec3::new(1., 0., 0.),
        0.5,
    );
    for support_speed in [0., 2.] {
        let base = BaseState::prescribed(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0., support_speed, 0.)),
            Twist::zero(),
        );
        let mut states = [prescribed(base, vec![0.], vec![0.])];
        let mut system = network(3., [0., support_speed + 2., 0.]);
        let arm = Vec3::new(2., 0., 0.);
        let mut description = contact(7, 0, 0, arm, modal_point(3., arm), 0., 0.5);
        description.normal_world = WorldVec3::new(0., 1., 0.);
        description.tangent_world = WorldVec3::new(-1., 0., 0.);
        let mut prepared = prepare(&system, vec![description], &gate);
        let receipt = system
            .apply_articulated_frictional_impacts(&mut states, &[&mechanism], &mut prepared, &gate)
            .unwrap();
        let impulse = 45. / 29.;
        close(receipt.contact_mobility_kg_inv[0], 4. / 2.5 + 1. / 3.);
        close(receipt.contacts.impulses[0][0], impulse);
        close(receipt.contact_work_j, -45. / 58.);
        close(receipt.support_work_j, support_speed * 0.6 * impulse);
        close(
            receipt.kinetic_energy_change_j,
            -45. / 58. + support_speed * 0.6 * impulse,
        );
        if support_speed > 0. {
            assert!(receipt.kinetic_energy_change_j > 0.);
        }
        close(receipt.energy_residual_j, 0.);
        vector(
            receipt.support_impulse_world.force,
            Vec3::new(0., 0.6 * impulse, 0.),
        );
        vector(receipt.support_impulse_world.torque, Vec3::default());
        let MixedArticulatedImpulseResponse::Prescribed(response) =
            &receipt.articulated_impulses[0]
        else {
            panic!("support reaction missing")
        };
        close(response.joint_velocity_after[0], 36. / 29.);
        vector(
            response.support_impulse_body.force,
            receipt.support_impulse_world.force,
        );
        let State::Prescribed(after) = &states[0] else {
            panic!("prescribed boundary changed")
        };
        assert_eq!(after.base, base);
        assert_eq!(after.time_s, 3.5);
        assert_eq!(after.joint_position, vec![0.]);
        close(after.joint_velocity[0], 36. / 29.);
    }
}

#[test]
fn simultaneous_free_and_prescribed_trees_share_modal_recoil_independent_of_contact_order() {
    let gate = CancelGate::new_clock_free();
    let free_model = single(JointModel::FIXED, 1., Vec3::default(), 1.);
    let slider = single(
        JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
        3.,
        Vec3::default(),
        1.,
    );
    let mut outcomes = Vec::new();
    for reverse in [false, true] {
        let mut states = [
            free(
                FreeFloatingBaseState::new(
                    Se3::identity(),
                    Twist::new(Vec3::default(), Vec3::new(-2.5, 0., 0.)),
                ),
                vec![],
                vec![],
            ),
            prescribed(
                BaseState::stationary(Se3::identity()),
                vec![0.],
                vec![-13. / 6.],
            ),
        ];
        let mut system = network(2., [0.; 3]);
        let mut contacts = vec![
            contact(
                31,
                0,
                0,
                Vec3::default(),
                modal_point(2., Vec3::default()),
                0.,
                0.,
            ),
            contact(
                7,
                1,
                0,
                Vec3::default(),
                modal_point(2., Vec3::default()),
                0.,
                0.,
            ),
        ];
        if reverse {
            contacts.reverse();
        }
        let mut prepared = prepare(&system, contacts, &gate);
        let receipt = system
            .apply_articulated_frictional_impacts(
                &mut states,
                &[&free_model, &slider],
                &mut prepared,
                &gate,
            )
            .unwrap();
        assert_eq!(receipt.contact_ids, vec![7, 31]);
        close(receipt.contact_mobility_kg_inv[3], 0.5);
        close(receipt.contacts.impulses[0][0], 2.);
        close(receipt.contacts.impulses[1][0], 1.);
        close(receipt.contact_work_j, -41. / 12.);
        close(receipt.support_work_j, 0.);
        close(receipt.energy_residual_j, 0.);
        let State::Free(free) = &states[0] else {
            panic!()
        };
        let State::Prescribed(slider) = &states[1] else {
            panic!()
        };
        close(free.base.twist_body.linear.x, -1.5);
        close(slider.joint_velocity[0], -1.5);
        close(
            modal_states(&system)[0].velocity_m_sqrt_kg_per_s / 2f64.sqrt(),
            -1.5,
        );
        outcomes.push((states, modal_states(&system), receipt));
    }
    assert_eq!(outcomes[0], outcomes[1]);
}

#[test]
fn inverse_mass_probe_ignores_fictitious_unit_speed_but_actual_limit_refusal_is_atomic() {
    let gate = CancelGate::new_clock_free();
    let limits = JointLimits::new(-1., 1., 1., 1.).unwrap();
    let mechanism = single(
        JointModel::prismatic(Vec3::new(1., 0., 0.), Some(limits)).unwrap(),
        0.01,
        Vec3::default(),
        0.01,
    );
    let unit = prescribed_base_impulse_delta(
        &mechanism,
        &[0.],
        &[0.],
        &[Wrench::new(Vec3::default(), Vec3::new(1., 0., 0.))],
    )
    .unwrap();
    close(unit.velocity.joint[0], 100.);
    let mut states = [prescribed(
        BaseState::stationary(Se3::identity()),
        vec![0.],
        vec![0.],
    )];
    let original_states = states.clone();
    let mut system = network(1., [2., 0., 0.]);
    let original_modal = modal_states(&system);
    let description = contact(
        9,
        0,
        0,
        Vec3::default(),
        modal_point(1., Vec3::default()),
        0.,
        0.,
    );
    let mut prepared = prepare(&system, vec![description], &gate);
    assert!(matches!(
        system.apply_articulated_frictional_impacts(
            &mut states,
            &[&mechanism],
            &mut prepared,
            &gate
        ),
        Err(ArticulatedModalImpactError::Dynamics { system: 0, .. })
    ));
    assert_eq!(states, original_states);
    assert_eq!(modal_states(&system), original_modal);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        system.apply_articulated_frictional_impacts(
            &mut states,
            &[&mechanism],
            &mut prepared,
            &cancelled
        ),
        Err(ArticulatedModalImpactError::Coupling(
            ModalCouplingError::Cancelled
        ))
    ));
    assert_eq!(states, original_states);
    assert_eq!(modal_states(&system), original_modal);
    // The same prepared rows and modal scratch remain reusable after refusal.
    system.apply_modal_impulse(&[-1.9, 0., 0.], &gate).unwrap();
    let receipt = system
        .apply_articulated_frictional_impacts(&mut states, &[&mechanism], &mut prepared, &gate)
        .unwrap();
    close(receipt.contacts.impulses[0][0], 0.1 / 101.);
    let State::Prescribed(after) = &states[0] else {
        panic!()
    };
    close(after.joint_velocity[0], 10. / 101.);
    assert!(after.joint_velocity[0].abs() < limits.velocity());
    assert_eq!(after.joint_position, vec![0.]);
    assert_eq!(after.time_s, 3.5);
}

#[test]
fn redundant_rows_preserve_native_determinacy_and_late_modal_budget_refusal() {
    let gate = CancelGate::new_clock_free();
    let mechanism = single(JointModel::FIXED, 1., Vec3::default(), 1.);
    let original = free(
        FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(-2., 0., 0.)),
        ),
        vec![],
        vec![],
    );
    let mut states = [original.clone()];
    let mut system = network(1., [0.; 3]);
    let a = contact(
        1,
        0,
        0,
        Vec3::default(),
        modal_point(1., Vec3::default()),
        0.,
        0.,
    );
    let mut b = a.clone();
    b.id = 2;
    let mut prepared = prepare(&system, vec![a.clone(), b], &gate);
    let receipt = system
        .apply_articulated_frictional_impacts(&mut states, &[&mechanism], &mut prepared, &gate)
        .unwrap();
    let FrictionalContactDeterminacy::UniqueVelocity(witness) = &receipt.contacts.determinacy
    else {
        panic!("redundant impulses became unique")
    };
    assert_eq!((witness.dimension, witness.rank), (6, 3));
    assert!(receipt.contacts.impulse_work_witness.is_some());
    close(receipt.contacts.impulses.iter().map(|p| p[0]).sum(), 1.);
    close(receipt.contact_work_j, -1.);
    assert!(
        PreparedArticulatedModalImpacts::new(
            &system,
            vec![a.clone(), a.clone()],
            Default::default(),
            &gate
        )
        .is_err()
    );
    let mut incomplete = a.clone();
    incomplete.modal_point.shapes_world.pop();
    assert!(
        PreparedArticulatedModalImpacts::new(&system, vec![incomplete], Default::default(), &gate)
            .is_err()
    );
    let mut budget = ModalAcousticTimeBudget::audible_reference();
    budget.maximum_abs_velocity_m_sqrt_kg_per_s = 0.6;
    let mut system =
        CoupledModalSystem::new(vec![model(2000, [0.; 3], budget)], vec![], config(), &gate)
            .unwrap();
    let before = modal_states(&system);
    let mut states = [original];
    let before_states = states.clone();
    let mut prepared = prepare(&system, vec![a], &gate);
    assert!(matches!(
        system.apply_articulated_frictional_impacts(
            &mut states,
            &[&mechanism],
            &mut prepared,
            &gate
        ),
        Err(ArticulatedModalImpactError::Coupling(
            ModalCouplingError::Component { component: 0, .. }
        ))
    ));
    assert_eq!(states, before_states);
    assert_eq!(modal_states(&system), before);
    let mut empty = prepare(&system, vec![], &gate);
    let noop = system
        .apply_articulated_frictional_impacts(&mut states, &[&mechanism], &mut empty, &gate)
        .unwrap();
    assert!(noop.contacts.impulses.is_empty());
    assert_eq!(states, before_states);
    assert_eq!(modal_states(&system), before);
}

fn live_network() -> CoupledModalSystem {
    let gate = CancelGate::new_clock_free();
    let mut receiver = ModalAcousticTimeModel::try_new(
        2000,
        vec![
            ModalAcousticMode {
                angular_frequency_rad_s: 40.,
                damping_ratio: 0.,
                pressure_per_modal_velocity: C64::ZERO
            };
            3
        ],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    receiver
        .restore_states(&[
            ModalAcousticState {
                displacement_m_sqrt_kg: 0.01,
                velocity_m_sqrt_kg_per_s: 0.,
            },
            ModalAcousticState::default(),
            ModalAcousticState::default(),
        ])
        .unwrap();
    CoupledModalSystem::new(vec![receiver], vec![], config(), &gate)
        .unwrap()
        .with_material_relaxation(
            vec![fs_phs::RelaxationBranch {
                projection: vec![1., 0., 0., 0., 0., 0.],
                stiffness: 200.,
                relaxation_time_s: 0.05,
            }],
            Some(vec![0.02]),
            &gate,
        )
        .unwrap()
        .with_radiation_load(
            &Model {
                ports: 1,
                poles: vec![Pole {
                    omega: 60.,
                    zeta: 0.1,
                    coupling: vec![10.],
                }],
            },
            &[0],
            1,
            &gate,
        )
        .unwrap()
}

#[test]
fn articulated_impulse_keeps_material_radiation_and_forced_modal_continuation() {
    let gate = CancelGate::new_clock_free();
    let mechanism = single(JointModel::FIXED, 4., Vec3::default(), 1.);
    let mut system = live_network();
    let mut oracle = live_network();
    let held = [1e-3, -2e-4, 0.];
    for _ in 0..3 {
        system.step_under_gate(&held, &gate).unwrap();
        oracle.step_under_gate(&held, &gate).unwrap();
    }
    let before = modal_states(&system);
    let material = system.material_memory().to_vec();
    let radiation = system.radiation_state().to_vec();
    assert!(material.iter().any(|&v| v != 0.));
    assert!(radiation.iter().any(|&v| v != 0.));
    let frame = system.last_frame().unwrap().clone();
    let position = Vec3::new(
        before[0].displacement_m_sqrt_kg,
        before[1].displacement_m_sqrt_kg,
        before[2].displacement_m_sqrt_kg,
    );
    let base = FreeFloatingBaseState::new(
        Se3::from_parts(So3::identity(), position).unwrap(),
        Twist::new(Vec3::default(), Vec3::new(-0.2, 0.4, 0.)),
    );
    let mut states = [free(base, vec![], vec![])];
    let mut prepared = prepare(
        &system,
        vec![contact(
            9,
            0,
            0,
            Vec3::default(),
            modal_point(1., Vec3::default()),
            0.3,
            0.25,
        )],
        &gate,
    );
    let receipt = system
        .apply_articulated_frictional_impacts(&mut states, &[&mechanism], &mut prepared, &gate)
        .unwrap();
    let normal = 0.2 + before[0].velocity_m_sqrt_kg_per_s;
    let tangent = -0.3 * normal;
    close(receipt.contacts.impulses[0][0], normal);
    close(receipt.contacts.impulses[0][1], tangent);
    assert!(receipt.relative_velocities_after_m_s[0][1] > 0.);
    let expected = apply_free_floating_impulses(
        &mechanism,
        base,
        &[],
        &[],
        &[],
        &[Wrench::new(Vec3::default(), Vec3::new(normal, tangent, 0.))],
    )
    .unwrap();
    let State::Free(after) = &states[0] else {
        panic!()
    };
    vector(
        after.base.twist_body.linear,
        expected.base_after.twist_body.linear,
    );
    oracle
        .apply_modal_impulse(&[-normal, -tangent, 0.], &gate)
        .unwrap();
    assert_eq!(system.material_memory(), material);
    assert_eq!(system.radiation_state(), radiation);
    assert_eq!(system.last_frame(), Some(&frame));
    assert_eq!(system.samples_rendered(), 3);
    for (old, new) in before.iter().zip(modal_states(&system)) {
        assert_eq!(old.displacement_m_sqrt_kg, new.displacement_m_sqrt_kg);
    }
    for _ in 0..16 {
        let actual = system.step_under_gate(&held, &gate).unwrap();
        let expected = oracle.step_under_gate(&held, &gate).unwrap();
        close(actual.total_energy_j(), expected.total_energy_j());
        close(actual.external_work_j, expected.external_work_j);
        for (actual, expected) in modal_states(&system).iter().zip(modal_states(&oracle)) {
            close(
                actual.displacement_m_sqrt_kg,
                expected.displacement_m_sqrt_kg,
            );
            close(
                actual.velocity_m_sqrt_kg_per_s,
                expected.velocity_m_sqrt_kg_per_s,
            );
        }
    }
}
