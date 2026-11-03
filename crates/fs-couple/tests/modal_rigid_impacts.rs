//! Original finite rigid inertia exchanging momentum with continuing modal states.
use fs_contact::frictional_multicontact::{
    FrictionalContactDeterminacy, FrictionalContactError, NormalVelocityLaw,
};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::plate::impact::radiation::{Model, Pole};
use fs_couple::render::schedule::force::coupled::impulse::impact::rigid::{
    ModalImpactPoint, PreparedRigidModalImpacts, RigidModalImpactConfig, RigidModalImpactContact,
    RigidModalImpactError,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalCouplingConfig, ModalCouplingError,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

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
fn close_vector(a: Vec3, b: Vec3) {
    close(a.x, b.x);
    close(a.y, b.y);
    close(a.z, b.z);
}
fn model(
    rate: u32,
    velocity: [f64; 3],
    displacement: [f64; 3],
    budget: ModalAcousticTimeBudget,
) -> ModalAcousticTimeModel {
    let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
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
    model
        .restore_states(
            &(0..3)
                .map(|axis| ModalAcousticState {
                    displacement_m_sqrt_kg: displacement[axis],
                    velocity_m_sqrt_kg_per_s: velocity[axis],
                })
                .collect::<Vec<_>>(),
        )
        .unwrap();
    model
}
fn network(rate: u32, mass: f64, velocity: [f64; 3]) -> CoupledModalSystem {
    CoupledModalSystem::new(
        vec![model(
            rate,
            velocity.map(|v| mass.sqrt() * v),
            [0.0; 3],
            ModalAcousticTimeBudget::audible_reference(),
        )],
        vec![],
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn properties(mass: f64) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap()
}
fn rigid(mass: f64, position: Vec3, velocity: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        velocity.scale(mass),
        Vec3::ZERO,
    )
    .unwrap()
}
fn states(system: &CoupledModalSystem) -> Vec<ModalAcousticState> {
    system
        .components()
        .iter()
        .flat_map(|model| model.states().iter().copied())
        .collect()
}
fn point(mass: f64, reference: Vec3) -> ModalImpactPoint {
    let a = mass.sqrt().recip();
    ModalImpactPoint {
        component: 0,
        reference_position_world_m: reference,
        shapes_world: vec![
            Vec3::new(a, 0.0, 0.0),
            Vec3::new(0.0, a, 0.0),
            Vec3::new(0.0, 0.0, a),
        ],
    }
}
fn contact(
    id: u64,
    point: ModalImpactPoint,
    arm: Vec3,
    mu: f64,
    e: f64,
) -> RigidModalImpactContact {
    RigidModalImpactContact {
        id,
        rigid_body: 0,
        arm_body_m: arm,
        modal_point: point,
        normal_world: Vec3::new(1.0, 0.0, 0.0),
        tangent_world: Vec3::new(0.0, 1.0, 0.0),
        activation_tolerance_m: 0.0,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: e },
    }
}
fn prepare(
    system: &CoupledModalSystem,
    contacts: Vec<RigidModalImpactContact>,
    gate: &CancelGate,
) -> PreparedRigidModalImpacts {
    PreparedRigidModalImpacts::new(system, contacts, RigidModalImpactConfig::default(), gate)
        .unwrap()
}

#[test]
fn finite_rigid_and_modal_masses_share_stick_slip_and_preserve_continuation() {
    let gate = CancelGate::new_clock_free();
    for rate in [1024, 2048] {
        for (mu, tangent, loss) in [(0.25, -0.6, 3.375), (2.0, -3.2, 7.6)] {
            for angle in [0.0_f64, 0.61] {
                let mut system = network(rate, 1.0, [0.0; 3]);
                system.step_under_gate(&[0.0; 3], &gate).unwrap();
                let frame = system.last_frame().unwrap().clone();
                let mut bodies = [rigid(4.0, Vec3::ZERO, Vec3::new(-2.0, 4.0, 0.0))];
                let original_pose = bodies[0].pose();
                let mut row = contact(17, point(1.0, Vec3::ZERO), Vec3::ZERO, mu, 0.5);
                row.tangent_world = Vec3::new(0.0, angle.cos(), angle.sin());
                let mut prepared = prepare(&system, vec![row], &gate);
                let report = system
                    .apply_rigid_frictional_impacts(
                        &mut bodies,
                        &[properties(4.0)],
                        &mut prepared,
                        &gate,
                    )
                    .unwrap();
                assert_eq!(report.contact_ids, vec![17]);
                close(report.contact_mobility_kg_inv[0], 1.25);
                close(report.contacts.impulses[0][0], 2.4);
                close(report.contacts.impulses[0][1], tangent * angle.cos());
                close(report.contacts.impulses[0][2], -tangent * angle.sin());
                close_vector(
                    bodies[0].linear_momentum_world(),
                    Vec3::new(-5.6, 16.0 + tangent, 0.0),
                );
                let after = states(&system);
                for (state, expected) in after.iter().zip([-2.4, -tangent, 0.0]) {
                    close(state.velocity_m_sqrt_kg_per_s, expected);
                    assert_eq!(state.displacement_m_sqrt_kg, 0.0);
                }
                close(
                    bodies[0].linear_momentum_world().x + after[0].velocity_m_sqrt_kg_per_s,
                    -8.0,
                );
                close(
                    bodies[0].linear_momentum_world().y + after[1].velocity_m_sqrt_kg_per_s,
                    16.0,
                );
                close(report.contacts.dissipated_energy, loss);
                close(report.kinetic_energy_change_j, -loss);
                close(report.contact_work_j, -loss);
                close(report.energy_residual_j, 0.0);
                assert_eq!(bodies[0].pose(), original_pose);
                assert_eq!(system.last_frame(), Some(&frame));
                assert_eq!(system.samples_rendered(), 1);
                system.step_under_gate(&[0.0; 3], &gate).unwrap();
                for (next, reset) in states(&system).iter().zip(after) {
                    close(
                        next.displacement_m_sqrt_kg,
                        reset.velocity_m_sqrt_kg_per_s / rate as f64,
                    );
                    close(
                        next.velocity_m_sqrt_kg_per_s,
                        reset.velocity_m_sqrt_kg_per_s,
                    );
                }
            }
        }
    }
    // A small resolved loss beside a large unchanged common kinetic baseline.
    let mut budget = ModalAcousticTimeBudget::audible_reference();
    budget.maximum_abs_velocity_m_sqrt_kg_per_s = 1e9;
    budget.maximum_total_energy_j = 1e18;
    let mut cfg = config();
    cfg.maximum_total_energy_j = 1e18;
    let mut system = CoupledModalSystem::new(
        vec![model(2000, [0.0, 1e8, 0.0], [0.0; 3], budget)],
        vec![],
        cfg,
        &gate,
    )
    .unwrap();
    let mut bodies = [rigid(1.0, Vec3::ZERO, Vec3::new(-2.0, 1e8, 0.0))];
    let mut prepared = prepare(
        &system,
        vec![contact(1, point(1.0, Vec3::ZERO), Vec3::ZERO, 0.0, 0.0)],
        &gate,
    );
    let report = system
        .apply_rigid_frictional_impacts(&mut bodies, &[properties(1.0)], &mut prepared, &gate)
        .unwrap();
    close(report.kinetic_energy_change_j, -1.0);
    close(report.contact_work_j, -1.0);
}

#[test]
fn off_center_full_tensor_response_keeps_local_geometry_at_a_distant_origin() {
    let gate = CancelGate::new_clock_free();
    let arm = Vec3::new(1.0, -1.0, 1.0);
    let property = MassProperties::from_inertia_tensor(
        1.0,
        Vec3::ZERO,
        fs_ga::Mat3 {
            m: [2.0, 0.5, 0.0, 0.5, 2.0, 0.0, 0.0, 0.0, 3.0],
        },
    )
    .unwrap();
    // Analytic full relative mobility for unit rigid/modal masses and this arm.
    let expected = [
        43.0 / 15.0,
        7.0 / 15.0,
        -0.4,
        7.0 / 15.0,
        43.0 / 15.0,
        0.4,
        -0.4,
        0.4,
        2.8,
    ];
    let impulse = Vec3::new(1.0, -0.3, 0.2);
    let velocity = Vec3::new(-397.0 / 150.0, 47.0 / 150.0, -1.0 / 25.0);
    for origin in [Vec3::ZERO, Vec3::new(1e16, -1e16, 1e16)] {
        let mut system = CoupledModalSystem::new(
            vec![model(
                2000,
                [0.0; 3],
                [1.0, -1.0, 1.0],
                ModalAcousticTimeBudget::audible_reference(),
            )],
            vec![],
            config(),
            &gate,
        )
        .unwrap();
        let mut bodies = [rigid(1.0, origin, velocity)];
        let mut prepared = prepare(
            &system,
            vec![contact(8, point(1.0, origin), arm, 2.0, 0.0)],
            &gate,
        );
        let report = system
            .apply_rigid_frictional_impacts(&mut bodies, &[property], &mut prepared, &gate)
            .unwrap();
        for (&actual, expected) in report.contact_mobility_kg_inv.iter().zip(expected) {
            close(actual, expected);
        }
        for (&actual, expected) in report.contacts.impulses[0].iter().zip([1.0, -0.3, 0.2]) {
            close(actual, expected);
        }
        for value in report.relative_velocities_after_m_s[0] {
            close(value, 0.0);
        }
        close_vector(bodies[0].angular_momentum_body(), Vec3::new(0.1, 0.8, 0.7));
        close_vector(bodies[0].linear_momentum_world(), velocity.add(impulse));
        let modal = states(&system);
        let modal_momentum = Vec3::new(
            modal[0].velocity_m_sqrt_kg_per_s,
            modal[1].velocity_m_sqrt_kg_per_s,
            modal[2].velocity_m_sqrt_kg_per_s,
        );
        close_vector(modal_momentum, impulse.scale(-1.0));
        // World-origin terms cancel through the independently balanced linear
        // impulses; the surviving local orbital impulse cancels central spin.
        close_vector(
            bodies[0]
                .angular_momentum_body()
                .add(arm.cross(modal_momentum)),
            Vec3::ZERO,
        );
        close(report.kinetic_energy_change_j, -4123.0 / 3000.0);
        close(report.energy_residual_j, 0.0);
        assert_eq!(report.gaps_m, [0.0]);
    }
}

#[test]
fn shared_rigid_body_and_modal_modes_enter_one_simultaneous_matrix() {
    let gate = CancelGate::new_clock_free();
    let first = contact(
        9,
        point(1.0, Vec3::new(0.0, 1.0, 0.0)),
        Vec3::new(0.0, 1.0, 0.0),
        0.0,
        0.0,
    );
    let second = contact(
        5,
        point(1.0, Vec3::new(0.0, -1.0, 0.0)),
        Vec3::new(0.0, -1.0, 0.0),
        0.0,
        0.0,
    );
    let mut outcomes = Vec::new();
    for rows in [vec![first.clone(), second.clone()], vec![second, first]] {
        let mut system = network(2000, 1.0, [0.0; 3]);
        let mut bodies = [rigid(1.0, Vec3::ZERO, Vec3::new(-2.0, 0.0, 0.0))];
        let mut prepared = prepare(&system, rows, &gate);
        let report = system
            .apply_rigid_frictional_impacts(&mut bodies, &[properties(1.0)], &mut prepared, &gate)
            .unwrap();
        assert_eq!(report.contact_ids, vec![5, 9]);
        close(report.contact_mobility_kg_inv[0], 3.0);
        close(report.contact_mobility_kg_inv[3], 1.0);
        close(report.contact_mobility_kg_inv[3 * 6 + 3], 3.0);
        close(report.contacts.impulses[0][0], 0.5);
        close(report.contacts.impulses[1][0], 0.5);
        close_vector(bodies[0].linear_momentum_world(), Vec3::new(-1.0, 0.0, 0.0));
        close_vector(bodies[0].angular_momentum_body(), Vec3::ZERO);
        close(states(&system)[0].velocity_m_sqrt_kg_per_s, -1.0);
        close(report.contacts.dissipated_energy, 1.0);
        outcomes.push((bodies, states(&system), report));
    }
    assert_eq!(outcomes[0], outcomes[1]);
}

#[test]
fn redundant_rows_keep_native_velocity_determinacy_and_current_rigid_inertia() {
    let gate = CancelGate::new_clock_free();
    let mut system = network(2000, 1.0, [0.0; 3]);
    let description = contact(1, point(1.0, Vec3::ZERO), Vec3::ZERO, 0.0, 0.0);
    let mut duplicate = description.clone();
    duplicate.id = 2;
    let mut prepared = prepare(&system, vec![description.clone(), duplicate], &gate);
    let mut bodies = [rigid(1.0, Vec3::ZERO, Vec3::new(-2.0, 0.0, 0.0))];
    let report = system
        .apply_rigid_frictional_impacts(&mut bodies, &[properties(1.0)], &mut prepared, &gate)
        .unwrap();
    let FrictionalContactDeterminacy::UniqueVelocity(witness) = &report.contacts.determinacy else {
        panic!("duplicate contact loads must remain a gauge");
    };
    assert_eq!((witness.dimension, witness.rank), (6, 3));
    assert!(report.contacts.impulse_work_witness.is_some());
    close(report.contacts.impulses.iter().map(|p| p[0]).sum(), 1.0);
    close(report.contact_work_j, -1.0);
    assert!(
        PreparedRigidModalImpacts::new(
            &system,
            vec![description.clone(); 2],
            Default::default(),
            &gate
        )
        .is_err()
    );
    let mut incomplete = description;
    incomplete.modal_point.shapes_world.pop();
    assert!(
        PreparedRigidModalImpacts::new(&system, vec![incomplete], Default::default(), &gate)
            .is_err()
    );
    // Reuse prepared rows with another actual body mass, pose and inertia.
    // Nothing at preparation is allowed to freeze the rigid response.
    let mut system = network(2000, 1.0, [0.0; 3]);
    let arm = Vec3::new(0.0, 1.0, 0.0);
    let mut prepared = prepare(
        &system,
        vec![contact(4, point(1.0, arm), arm, 0.0, 0.0)],
        &gate,
    );
    for (mass, polar, angle) in [(1.0, 1.0, 0.0), (4.0, 2.0, std::f64::consts::FRAC_PI_2)] {
        let rotation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), angle).unwrap();
        let world_arm = rotation.rotate_body_to_world(arm);
        let position = arm.sub(world_arm);
        let property =
            MassProperties::new(mass, Vec3::ZERO, Vec3::new(polar, polar, polar)).unwrap();
        let mut bodies = [RigidBodyState::new(
            Pose::new(position, rotation).unwrap(),
            Vec3::new(-2.0 * mass, 0.0, 0.0),
            Vec3::ZERO,
        )
        .unwrap()];
        let before_modal = states(&system)[0].velocity_m_sqrt_kg_per_s;
        let expected = 1.0 + mass.recip() + angle.cos().powi(2) / polar;
        let report = system
            .apply_rigid_frictional_impacts(&mut bodies, &[property], &mut prepared, &gate)
            .unwrap();
        close(report.contact_mobility_kg_inv[0], expected);
        close(
            report.contacts.impulses[0][0],
            (2.0 + before_modal) / expected,
        );
    }
}

#[test]
fn every_receiver_remains_unchanged_after_late_state_momentum_or_geometry_refusal() {
    let gate = CancelGate::new_clock_free();
    let mut budget = ModalAcousticTimeBudget::audible_reference();
    budget.maximum_abs_velocity_m_sqrt_kg_per_s = 0.6;
    let mut system = CoupledModalSystem::new(
        vec![model(2000, [0.0; 3], [0.0; 3], budget)],
        vec![],
        config(),
        &gate,
    )
    .unwrap();
    let mut bodies = [rigid(1.0, Vec3::ZERO, Vec3::new(-1.0, 4.0, 0.0))];
    let before_bodies = bodies;
    let before_modal = states(&system);
    let description = contact(1, point(1.0, Vec3::ZERO), Vec3::ZERO, 2.0, 0.0);
    let mut prepared = prepare(&system, vec![description.clone()], &gate);
    assert!(matches!(
        system.apply_rigid_frictional_impacts(
            &mut bodies,
            &[properties(1.0)],
            &mut prepared,
            &gate
        ),
        Err(RigidModalImpactError::Coupling(
            ModalCouplingError::Component { component: 0, .. }
        ))
    ));
    assert_eq!(bodies, before_bodies);
    assert_eq!(states(&system), before_modal);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        system.apply_rigid_frictional_impacts(
            &mut bodies,
            &[properties(1.0)],
            &mut prepared,
            &cancelled
        ),
        Err(RigidModalImpactError::Coupling(
            ModalCouplingError::Cancelled
        ))
    ));
    let mut source = description.clone();
    source.friction_coefficient = 0.0;
    source.normal_law = NormalVelocityLaw::Target { speed: 2.0 };
    let mut source = prepare(&system, vec![source], &gate);
    assert!(matches!(
        system.apply_rigid_frictional_impacts(&mut bodies, &[properties(1.0)], &mut source, &gate),
        Err(RigidModalImpactError::Contact(
            FrictionalContactError::EnergyIncrease { .. }
        ))
    ));
    let mut gap = description.clone();
    gap.friction_coefficient = 0.5;
    gap.activation_tolerance_m = 1e-3;
    let mut gap = prepare(&system, vec![gap], &gate);
    let mut separated = [rigid(
        1.0,
        Vec3::new(1e-4, 0.0, 0.0),
        Vec3::new(-1.0, 4.0, 0.0),
    )];
    let original_separated = separated;
    assert!(matches!(
        system.apply_rigid_frictional_impacts(&mut separated, &[properties(1.0)], &mut gap, &gate),
        Err(RigidModalImpactError::Coupling(
            ModalCouplingError::Invalid { .. }
        ))
    ));
    assert_eq!(separated, original_separated);
    assert_eq!(states(&system), before_modal);
    let mut penetrating = [rigid(
        1.0,
        Vec3::new(-1e-10, 0.0, 0.0),
        Vec3::new(-1.0, 4.0, 0.0),
    )];
    assert!(matches!(
        system.apply_rigid_frictional_impacts(
            &mut penetrating,
            &[properties(1.0)],
            &mut gap,
            &gate
        ),
        Err(RigidModalImpactError::Contact(
            FrictionalContactError::PenetratingContact { .. }
        ))
    ));
    let mut empty = prepare(&system, vec![], &gate);
    assert!(
        system
            .apply_rigid_frictional_impacts(&mut bodies, &[properties(1.0)], &mut empty, &gate)
            .unwrap()
            .contacts
            .impulses
            .is_empty()
    );
    assert_eq!(bodies, before_bodies);
    assert_eq!(states(&system), before_modal);
    let mut retry = description;
    retry.friction_coefficient = 0.5;
    let mut retry = prepare(&system, vec![retry], &gate);
    system
        .apply_rigid_frictional_impacts(&mut bodies, &[properties(1.0)], &mut retry, &gate)
        .unwrap();
    close(states(&system)[0].velocity_m_sqrt_kg_per_s, -0.5);
    close(states(&system)[1].velocity_m_sqrt_kg_per_s, 0.25);

    // Modal staging succeeds, then the rigid body's 1 N s is swallowed by
    // its 1e16 N s background momentum. The modal reset must not leak through.
    let mut system = network(2000, 1.0, [2.0, 0.0, 0.0]);
    let before = states(&system);
    let mut bodies = [rigid(1e16, Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0))];
    let original = bodies;
    let mut prepared = prepare(
        &system,
        vec![contact(1, point(1.0, Vec3::ZERO), Vec3::ZERO, 0.0, 0.0)],
        &gate,
    );
    assert!(matches!(
        system.apply_rigid_frictional_impacts(
            &mut bodies,
            &[properties(1e16)],
            &mut prepared,
            &gate
        ),
        Err(RigidModalImpactError::Momentum {
            body: 0,
            coordinate: 0,
            ..
        })
    ));
    assert_eq!(bodies, original);
    assert_eq!(states(&system), before);
    bodies[0] = rigid(1.0, Vec3::ZERO, Vec3::new(-1.0, 0.0, 0.0));
    system
        .apply_rigid_frictional_impacts(&mut bodies, &[properties(1.0)], &mut prepared, &gate)
        .unwrap();
    close(states(&system)[0].velocity_m_sqrt_kg_per_s, 0.5);
}

fn live_network() -> CoupledModalSystem {
    let gate = CancelGate::new_clock_free();
    let mut receiver = ModalAcousticTimeModel::try_new(
        2000,
        vec![
            ModalAcousticMode {
                angular_frequency_rad_s: 40.0,
                damping_ratio: 0.0,
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
                velocity_m_sqrt_kg_per_s: 0.0,
            },
            ModalAcousticState::default(),
            ModalAcousticState::default(),
        ])
        .unwrap();
    CoupledModalSystem::new(vec![receiver], vec![], config(), &gate)
        .unwrap()
        .with_material_relaxation(
            vec![fs_phs::RelaxationBranch {
                projection: vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                stiffness: 200.0,
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
                    omega: 60.0,
                    zeta: 0.1,
                    coupling: vec![10.0],
                }],
            },
            &[0],
            1,
            &gate,
        )
        .unwrap()
}

#[test]
fn actual_rigid_recoil_preserves_material_radiation_and_forced_modal_continuation() {
    let gate = CancelGate::new_clock_free();
    let mut system = live_network();
    let mut oracle = live_network();
    let held = [1e-3, -2e-4, 0.0];
    for _ in 0..3 {
        system.step_under_gate(&held, &gate).unwrap();
        oracle.step_under_gate(&held, &gate).unwrap();
    }
    let before = states(&system);
    let material = system.material_memory().to_vec();
    let radiation = system.radiation_state().to_vec();
    assert!(material.iter().any(|&v| v != 0.0));
    assert!(radiation.iter().any(|&v| v != 0.0));
    let frame = system.last_frame().unwrap().clone();
    let position = Vec3::new(
        before[0].displacement_m_sqrt_kg,
        before[1].displacement_m_sqrt_kg,
        before[2].displacement_m_sqrt_kg,
    );
    let original_body = rigid(4.0, position, Vec3::new(-0.2, 0.4, 0.0));
    let mut bodies = [original_body];
    let mut prepared = prepare(
        &system,
        vec![contact(9, point(1.0, Vec3::ZERO), Vec3::ZERO, 0.3, 0.25)],
        &gate,
    );
    let report = system
        .apply_rigid_frictional_impacts(&mut bodies, &[properties(4.0)], &mut prepared, &gate)
        .unwrap();
    // Diagonal mobility 5/4 and e=1/4: Jn=-(v_r-v_m); sliding Jt=-0.3Jn.
    let normal = 0.2 + before[0].velocity_m_sqrt_kg_per_s;
    let tangent = -0.3 * normal;
    close(report.contacts.impulses[0][0], normal);
    close(report.contacts.impulses[0][1], tangent);
    assert!(report.relative_velocities_after_m_s[0][1] > 0.0);
    oracle
        .apply_modal_impulse(&[-normal, -tangent, 0.0], &gate)
        .unwrap();
    let expected_body = original_body
        .apply_impulse_at_body_point(properties(4.0), Vec3::ZERO, Vec3::new(normal, tangent, 0.0))
        .unwrap()
        .state_after;
    close_vector(
        bodies[0].linear_momentum_world(),
        expected_body.linear_momentum_world(),
    );
    assert_eq!(system.material_memory(), material);
    assert_eq!(system.radiation_state(), radiation);
    assert_eq!(system.last_frame(), Some(&frame));
    assert_eq!(system.samples_rendered(), 3);
    for (old, new) in before.iter().zip(states(&system)) {
        assert_eq!(old.displacement_m_sqrt_kg, new.displacement_m_sqrt_kg);
    }
    for _ in 0..16 {
        let actual = system.step_under_gate(&held, &gate).unwrap();
        let expected = oracle.step_under_gate(&held, &gate).unwrap();
        close(actual.total_energy_j(), expected.total_energy_j());
        close(actual.external_work_j, expected.external_work_j);
        for (actual, expected) in states(&system).iter().zip(states(&oracle)) {
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
