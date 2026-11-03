//! G1/G3/G4: original source coefficients through joint and floating-base recoil.
use super::sphere_plane::{original_claims, source_fixture};
use super::*;
use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_contact::interface_binding::thermal_friction::elliptic::SourceFrictionConditions;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::SourceImpulseCoefficientRole;
use fs_couple::friction_heat::directional::articulated::*;
use fs_couple::friction_heat::directional::impact::{
    SourceRigidImpactBinding, SourceRigidImpactError,
};
use fs_ga::{Mat3, Se3, So3, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath,
};
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};
use fs_mbd::articulated_contacts::friction::elliptic::{
    ArticulatedEllipticSourceError, ArticulatedImpactGeometry,
};
use fs_mbd::articulated_contacts::mixed::{
    ArticulatedContactBoundary as Boundary, MixedArticulatedContactSystem as System,
    MixedArticulatedImpulseResponse as Response,
};
use fs_mbd::articulated_contacts::step::mixed::MixedArticulatedContactPath;
use fs_mbd::articulated_contacts::step::mixed::elliptic::MixedArticulatedEllipticSourceStepError;
use fs_mbd::articulated_contacts::step::{
    ArticulatedContactStepConfig, ArticulatedContactStepError,
};
use fs_mbd::articulated_contacts::trajectory::mixed::{
    MixedArticulatedTrajectoryState as StepState, MixedArticulatedTrajectorySystem as StepSystem,
};
use fs_mbd::articulated_contacts::{
    ArticulatedContactError, ArticulatedContactSystem, ContactCounterpart, LinkContactPoint,
};

fn slider() -> ArticulatedModel {
    let inertia = |mass| {
        SpatialInertia::new(
            mass,
            Vec3::default(),
            Mat3 {
                m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            },
        )
        .unwrap()
    };
    ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(2.0),
        ),
        Link::new(
            "slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(3.0),
        ),
    ])
    .unwrap()
}

fn base(z: f64, falling_speed: f64) -> FreeFloatingBaseState {
    FreeFloatingBaseState::new(
        Se3::from_parts(So3::identity(), Vec3::new(0.0, 0.0, z)).unwrap(),
        Twist::new(Vec3::default(), Vec3::new(0.0, 0.0, -falling_speed)),
    )
}

fn contact(id: u64, system: usize) -> ArticulatedImpactGeometry {
    ArticulatedImpactGeometry {
        id,
        point: LinkContactPoint {
            system_index: system,
            link_index: 1,
            point_body: Vec3::default(),
        },
        other: ContactCounterpart::PrescribedPlane {
            // A remote anchor represents the SAME plane, not the source point.
            point_world: Vec3::new(2.0, 3.0, 0.0),
            velocity_world: Vec3::default(),
        },
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance: 0.0,
        friction_coefficient: (),
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    }
}

fn conditions(query: &QueryPoint) -> SourceFrictionConditions {
    SourceFrictionConditions {
        query: query.clone(),
        contact_temperature_k: Some(query.axes()["Ti"]),
        nominal_pressure_pa: None,
    }
}

fn binding(law: &SourceThermalEllipticCoulomb) -> SourceRigidImpactBinding<'_> {
    SourceRigidImpactBinding {
        law,
        coefficient_role: SourceImpulseCoefficientRole::Kinetic,
    }
}

fn balance(report: &SourceArticulatedImpactCandidate) {
    let m = &report.mechanical;
    let change: f64 = m
        .responses
        .iter()
        .map(|r| r.momentum_after.kinetic_energy - r.momentum_before.kinetic_energy)
        .sum();
    close(change - m.plane_work + m.contact.dissipated_energy, 0.0);
    close(m.energy_balance_residual, 0.0);
    close(m.endpoint_velocity_residual, 0.0);
    for v in [m.linear_momentum_residual, m.angular_momentum_residual] {
        for x in [v.x, v.y, v.z] {
            close(x, 0.0);
        }
    }
}

#[test]
fn source_friction_changes_joint_speed_and_normal_base_recoil_from_original_inertia() {
    let (card, law, query) = source_fixture(curve("speed", SPEED, [(0.0, 0.1), (12.0, 0.4)]));
    let model = slider();
    let q = [0.0];
    let qd = [4.0];
    let initial_base = base(0.0, 1.0);
    let systems = [ArticulatedContactSystem {
        model: &model,
        base: initial_base,
        joint_position: &q,
        joint_velocity: &qd,
    }];
    let original = law.initial_state();
    for restitution in [0.0, 0.5] {
        let mut c = contact(17, 0);
        c.normal_law = NormalVelocityLaw::NewtonRestitution {
            coefficient: restitution,
        };
        let mut calls = 0;
        let report = resolve_source_articulated_elliptic_contacts(
            &systems,
            &[c],
            core::slice::from_ref(&original),
            &[binding(&law)],
            FrictionalContactConfig::default(),
            |accepted, k| {
                calls += 1;
                assert_eq!(accepted, &original);
                assert_eq!(k.contact_id, 17);
                assert_eq!(k.original_index, 0);
                assert_eq!(k.point_world, Vec3::default());
                assert_eq!(k.counterpart_point_world, Vec3::default());
                assert_eq!(k.relative_velocity_m_s, [-1.0, 4.0, 0.0]);
                Ok::<_, &'static str>(conditions(&query))
            },
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(calls, 1);
        // A free prismatic x coordinate decouples horizontal translation:
        // Wnn=1/(2+3), Wxx=1/3, Wyy=1/(2+3). Both masses share z recoil.
        let pn = 5.0 * (1.0 + restitution);
        let pt = -0.2 * pn;
        let joint_after = 4.0 + pt / 3.0;
        close(report.mechanical.delassus[0], 0.2);
        close(report.mechanical.delassus[4], 1.0 / 3.0);
        close(report.mechanical.delassus[8], 0.2);
        close(report.mechanical.contact.impulses[0][0], pn);
        close(report.mechanical.contact.impulses[0][1], pt);
        close(report.mechanical.contact.impulses[0][2], 0.0);
        let response = &report.mechanical.responses[0];
        close(response.base_after.twist_body.linear.x, 0.0);
        close(response.base_after.twist_body.linear.z, restitution);
        close(response.joint_velocity_after[0], joint_after);
        close(response.momentum_after.linear_momentum_world.x, 12.0 + pt);
        close(
            response.momentum_after.linear_momentum_world.z,
            5.0 * restitution,
        );
        assert_eq!(
            response.base_after.world_from_base,
            initial_base.world_from_base
        );
        let expected_energy = 0.5 * 5.0 * restitution.powi(2) + 0.5 * 3.0 * joint_after.powi(2);
        close(response.momentum_after.kinetic_energy, expected_energy);
        close(
            report.mechanical.contact.dissipated_energy,
            26.5 - expected_energy,
        );
        let sample = &report.sources[0].sample;
        close(sample.queried_speed_m_s, 4.0);
        close(sample.coefficients_mu[0], 0.2);
        close(sample.coefficients_mu[1], 0.4);
        assert_eq!(sample.state.constitutive_state().slip_m(), [0.0; 2]);
        for source in sample.state.friction_sources() {
            close(source.state_coordinate("path").unwrap(), 0.25);
            close(source.state_coordinate("speed").unwrap(), 4.0);
        }
        original_claims(&card, &original, &sample.state);
        assert_eq!(systems[0].base, initial_base);
        assert_eq!(systems[0].joint_position, [0.0]);
        assert_eq!(systems[0].joint_velocity, [4.0]);
        balance(&report);
    }
}

#[test]
fn distinct_source_bindings_follow_canonical_contacts_and_skip_inactive_history() {
    let (card_a, law_a, query_a) = source_fixture(curve("speed", SPEED, [(0.0, 0.1), (12.0, 0.4)]));
    let (card_b, law_b, query_b) = source_fixture(scalar(0.3));
    let model = slider();
    let q = [0.0];
    let first_velocity = [4.0];
    let second_velocity = [-3.0];
    let systems = [
        ArticulatedContactSystem {
            model: &model,
            base: base(0.0, 1.0),
            joint_position: &q,
            joint_velocity: &first_velocity,
        },
        ArticulatedContactSystem {
            model: &model,
            base: base(0.0, 2.0),
            joint_position: &q,
            joint_velocity: &second_velocity,
        },
    ];
    let mut inactive = contact(5, 0);
    inactive.point.link_index = 0;
    inactive.other = ContactCounterpart::PrescribedPlane {
        point_world: Vec3::new(0.0, 0.0, -10.0),
        velocity_world: Vec3::default(),
    };
    let geometry = [contact(20, 0), contact(30, 1), inactive];
    let original = [
        law_a.initial_state(),
        law_b.initial_state(),
        law_a.initial_state(),
    ];
    let bindings = [binding(&law_a), binding(&law_b), binding(&law_a)];
    let mut snapshot = None;
    for order in [[0, 1, 2], [2, 0, 1], [1, 2, 0]] {
        let contacts = order.map(|i| geometry[i]);
        let histories = order.map(|i| original[i].clone());
        let before = histories.clone();
        let selected = order.map(|i| bindings[i]);
        let mut ids = Vec::new();
        let report = resolve_source_articulated_elliptic_contacts(
            &systems,
            &contacts,
            &histories,
            &selected,
            FrictionalContactConfig::default(),
            |accepted, k| {
                ids.push(k.contact_id);
                assert_eq!(accepted, &histories[k.original_index]);
                match k.contact_id {
                    20 => Ok::<_, &'static str>(conditions(&query_a)),
                    30 => Ok(conditions(&query_b)),
                    _ => panic!("inactive source was queried"),
                }
            },
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(ids, [20, 30]);
        assert_eq!(histories, before);
        assert_eq!(report.mechanical.contact_ids, [5, 20, 30]);
        assert_eq!(report.mechanical.contact.impulses[0], [0.0; 3]);
        // The inactive base point still receives the original system's z recoil.
        close(report.mechanical.relative_velocities_after[0][0], 0.0);
        for (row, pn, pt) in [(1, 5.0, -1.0), (2, 10.0, 3.0)] {
            close(report.mechanical.contact.impulses[row][0], pn);
            close(report.mechanical.contact.impulses[row][1], pt);
        }
        close(
            report.mechanical.responses[0].joint_velocity_after[0],
            11.0 / 3.0,
        );
        close(report.mechanical.responses[1].joint_velocity_after[0], -2.0);
        for source in &report.sources {
            assert_eq!(
                contacts[source.kinematics.original_index].id,
                source.kinematics.contact_id
            );
            let (card, old) = if source.kinematics.contact_id == 20 {
                (&card_a, &original[0])
            } else {
                (&card_b, &original[1])
            };
            original_claims(card, old, &source.sample.state);
        }
        let value = (
            report.mechanical.responses.clone(),
            report.mechanical.contact.impulses.clone(),
        );
        if let Some(previous) = &snapshot {
            assert_eq!(&value, previous);
        } else {
            snapshot = Some(value);
        }
        balance(&report);
    }
}

#[test]
fn later_source_refusal_and_cancellation_return_no_partial_articulated_candidate() {
    let (_, law, query) = source_fixture(scalar(0.2));
    let model = slider();
    let q = [0.0];
    let qd = [4.0];
    let systems = [0, 1].map(|_| ArticulatedContactSystem {
        model: &model,
        base: base(0.0, 1.0),
        joint_position: &q,
        joint_velocity: &qd,
    });
    let contacts = [contact(30, 1), contact(20, 0)];
    let histories = [law.initial_state(), law.initial_state()];
    let before = histories.clone();
    for cancel in [false, true] {
        let gate = CancelGate::new();
        let mut calls = 0;
        let result = resolve_source_articulated_elliptic_contacts(
            &systems,
            &contacts,
            &histories,
            &[binding(&law); 2],
            FrictionalContactConfig::default(),
            |_, k| {
                calls += 1;
                let mut boundary = conditions(&query);
                if k.contact_id == 30 {
                    if cancel {
                        gate.request();
                    } else {
                        boundary.query = boundary.query.with("Ti", 450.0).unwrap();
                        boundary.contact_temperature_k = Some(450.0);
                    }
                }
                Ok::<_, &'static str>(boundary)
            },
            &gate,
        );
        assert!(matches!(
            result,
            Err(ArticulatedEllipticSourceError::Source { contact_id: 30, .. })
        ));
        assert_eq!(calls, 2);
        assert_eq!(histories, before);
        assert_eq!(systems[0].joint_velocity, [4.0]);
        assert_eq!(systems[1].joint_velocity, [4.0]);
    }
}

#[test]
fn late_native_momentum_refusal_discards_sources_and_touching_retry_is_deterministic() {
    let (_, law, query) = source_fixture(curve("speed", SPEED, [(0.0, 0.1), (12.0, 0.4)]));
    let model = slider();
    let q = [0.0];
    let upper_velocity = [4.0];
    let lower_velocity = [0.0];
    let histories = [law.initial_state()];
    let before = histories.clone();
    let mut c = contact(17, 0);
    c.other = ContactCounterpart::Body(LinkContactPoint {
        system_index: 1,
        link_index: 1,
        point_body: Vec3::default(),
    });
    c.activation_tolerance = 1e-3;
    let mut accepted = None;
    for gap in [1e-4, 0.0, 0.0] {
        let systems = [
            ArticulatedContactSystem {
                model: &model,
                base: base(gap, 1.0),
                joint_position: &q,
                joint_velocity: &upper_velocity,
            },
            ArticulatedContactSystem {
                model: &model,
                base: base(0.0, 0.0),
                joint_position: &q,
                joint_velocity: &lower_velocity,
            },
        ];
        let mut calls = 0;
        let result = resolve_source_articulated_elliptic_contacts(
            &systems,
            &[c],
            &histories,
            &[binding(&law)],
            FrictionalContactConfig::default(),
            |_, k| {
                calls += 1;
                close(k.gap_m, gap);
                Ok::<_, &'static str>(conditions(&query))
            },
            &CancelGate::new(),
        );
        assert_eq!(calls, 1);
        assert_eq!(histories, before);
        if gap > 0.0 {
            // Opposite tangent impulses at separated points would introduce an
            // artificial internal couple. The original late audit must refuse.
            assert!(matches!(
                result,
                Err(ArticulatedEllipticSourceError::Mechanics(
                    ArticulatedContactError::Unresolved(
                        "actual global impact momentum failed admission"
                    )
                ))
            ));
        } else {
            let report = result.unwrap();
            balance(&report);
            close(report.mechanical.contact.impulses[0][0], 2.5);
            close(report.mechanical.contact.impulses[0][1], -0.5);
            assert_ne!(report.sources[0].sample.state, histories[0]);
            let value = (
                report.mechanical.responses,
                report.sources[0].sample.state.clone(),
            );
            if let Some(previous) = &accepted {
                assert_eq!(&value, previous);
            } else {
                accepted = Some(value);
            }
        }
    }
}

#[test]
fn mixed_source_friction_recoils_both_joints_and_retains_driven_support_work() {
    let (card, law, query) = source_fixture(curve("speed", SPEED, [(0.0, 0.1), (12.0, 0.4)]));
    let model = slider();
    let q = [0.0];
    let upper_rate = [4.0];
    let lower_rate = [0.0];
    let support_velocity = Vec3::new(0.5, -0.2, 0.6);
    let support = BaseState::prescribed(
        Se3::identity(),
        Twist::new(Vec3::default(), support_velocity),
        Twist::zero(),
    );
    let upper = FreeFloatingBaseState::new(
        Se3::identity(),
        Twist::new(
            Vec3::default(),
            support_velocity + Vec3::new(0.0, 0.0, -1.0),
        ),
    );
    let systems = [
        System {
            model: &model,
            boundary: Boundary::Free(upper),
            joint_position: &q,
            joint_velocity: &upper_rate,
        },
        System {
            model: &model,
            boundary: Boundary::Prescribed(support),
            joint_position: &q,
            joint_velocity: &lower_rate,
        },
    ];
    let histories = [law.initial_state()];
    let before = histories.clone();
    for restitution in [0.0, 0.5] {
        let mut c = contact(17, 0);
        c.other = ContactCounterpart::Body(LinkContactPoint {
            system_index: 1,
            link_index: 1,
            point_body: Vec3::default(),
        });
        c.normal_law = NormalVelocityLaw::NewtonRestitution {
            coefficient: restitution,
        };
        let mut calls = 0;
        let report = resolve_source_mixed_articulated_elliptic_contacts(
            &systems,
            &[c],
            &histories,
            &[binding(&law)],
            FrictionalContactConfig::default(),
            |accepted, k| {
                calls += 1;
                assert_eq!(accepted, &histories[0]);
                assert_eq!(k.contact_id, 17);
                assert_eq!(k.original_index, 0);
                assert_eq!(k.point_world, Vec3::default());
                assert_eq!(k.counterpart_point_world, Vec3::default());
                for (actual, expected) in k.relative_velocity_m_s.into_iter().zip([-1.0, 4.0, 0.0])
                {
                    close(actual, expected);
                }
                Ok::<_, &'static str>(conditions(&query))
            },
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(calls, 1);
        assert_eq!(report.sources.len(), 1);
        let m = &report.mechanical;
        // Only the free tree recoils normally; both 3 kg sliders recoil along
        // the material x axis. The original speed query gives mu_x = 0.2.
        let pn = 5.0 * (1.0 + restitution);
        let pt = -0.2 * pn;
        let slip_after = 4.0 + 2.0 * pt / 3.0;
        let loss = 0.5 * pn * (1.0 - restitution) - 0.5 * pt * (4.0 + slip_after);
        close(m.delassus[0], 0.2);
        close(m.delassus[4], 2.0 / 3.0);
        close(m.delassus[8], 0.2);
        for (actual, expected) in m.contact.impulses[0].into_iter().zip([pn, pt, 0.0]) {
            close(actual, expected);
        }
        let Response::Free(free_after) = &m.responses[0] else {
            panic!("wrong boundary")
        };
        let Response::Prescribed(supported_after) = &m.responses[1] else {
            panic!("wrong boundary")
        };
        close(free_after.joint_velocity_after[0], 4.0 + pt / 3.0);
        close(supported_after.joint_velocity_after[0], -pt / 3.0);
        close(
            free_after.base_after.twist_body.linear.z,
            support_velocity.z + restitution,
        );
        assert_eq!(supported_after.base_after, support);
        for (actual, expected) in [
            (supported_after.support_impulse_world.force.x, 0.0),
            (supported_after.support_impulse_world.force.y, 0.0),
            (supported_after.support_impulse_world.force.z, pn),
            (supported_after.support_impulse_world.torque.x, 0.0),
            (supported_after.support_impulse_world.torque.y, 0.0),
            (supported_after.support_impulse_world.torque.z, 0.0),
        ] {
            close(actual, expected);
        }
        close(m.plane_work, 0.0);
        close(m.support_work, pn * support_velocity.z);
        close(m.boundary_work, pn * support_velocity.z);
        close(m.contact.dissipated_energy, loss);
        close(
            free_after.kinetic_energy_change + supported_after.kinetic_energy_change,
            m.boundary_work - loss,
        );
        close(m.energy_balance_residual, 0.0);
        close(m.endpoint_velocity_residual, 0.0);
        for residual in [m.linear_momentum_residual, m.angular_momentum_residual] {
            for value in [residual.x, residual.y, residual.z] {
                close(value, 0.0);
            }
        }
        let sample = &report.sources[0].sample;
        close(sample.queried_speed_m_s, 4.0);
        close(sample.coefficients_mu[0], 0.2);
        close(sample.coefficients_mu[1], 0.4);
        assert_eq!(sample.state.constitutive_state().slip_m(), [0.0; 2]);
        for source in sample.state.friction_sources() {
            close(source.state_coordinate("path").unwrap(), 0.25);
            close(source.state_coordinate("speed").unwrap(), 4.0);
        }
        original_claims(&card, &histories[0], &sample.state);
        assert_eq!(histories, before);
        assert_eq!(systems[0].boundary, Boundary::Free(upper));
        assert_eq!(systems[1].boundary, Boundary::Prescribed(support));
        assert_eq!(upper_rate, [4.0]);
        assert_eq!(lower_rate, [0.0]);
    }
}

#[test]
fn mixed_late_momentum_refusal_discards_sampled_history_before_touching_retry() {
    let (_, law, query) = source_fixture(curve("speed", SPEED, [(0.0, 0.1), (12.0, 0.4)]));
    let model = slider();
    let histories = [law.initial_state()];
    let before = histories.clone();
    let q = [0.0];
    let upper_rate = [4.0];
    let lower_rate = [0.0];
    let support = BaseState::stationary(Se3::identity());
    let mut c = contact(17, 0);
    c.other = ContactCounterpart::Body(LinkContactPoint {
        system_index: 1,
        link_index: 1,
        point_body: Vec3::default(),
    });
    c.activation_tolerance = 1e-3;
    let mut snapshot = None;
    for gap in [1e-4, 0.0, 0.0] {
        let upper = base(gap, 1.0);
        let systems = [
            System {
                model: &model,
                boundary: Boundary::Free(upper),
                joint_position: &q,
                joint_velocity: &upper_rate,
            },
            System {
                model: &model,
                boundary: Boundary::Prescribed(support),
                joint_position: &q,
                joint_velocity: &lower_rate,
            },
        ];
        let mut calls = 0;
        let result = resolve_source_mixed_articulated_elliptic_contacts(
            &systems,
            &[c],
            &histories,
            &[binding(&law)],
            FrictionalContactConfig::default(),
            |_, k| {
                calls += 1;
                close(k.gap_m, gap);
                Ok::<_, &'static str>(conditions(&query))
            },
            &CancelGate::new(),
        );
        assert_eq!(calls, 1);
        assert_eq!(histories, before);
        assert_eq!(systems[0].boundary, Boundary::Free(upper));
        assert_eq!(systems[1].boundary, Boundary::Prescribed(support));
        assert_eq!(upper_rate, [4.0]);
        assert_eq!(lower_rate, [0.0]);
        if gap > 0.0 {
            // Opposite tangent impulses at separated material points introduce
            // an internal couple even when a real support carries the recoil.
            assert!(matches!(
                result,
                Err(ArticulatedEllipticSourceError::Mechanics(
                    ArticulatedContactError::Unresolved(
                        "actual global impact momentum failed admission"
                    )
                ))
            ));
        } else {
            let report = result.unwrap();
            close(report.mechanical.contact.impulses[0][0], 5.0);
            close(report.mechanical.contact.impulses[0][1], -1.0);
            close(report.mechanical.boundary_work, 0.0);
            close(report.mechanical.energy_balance_residual, 0.0);
            assert_ne!(report.sources[0].sample.state, histories[0]);
            let value = (
                report.mechanical.responses,
                report.sources[0].sample.state.clone(),
            );
            if let Some(previous) = &snapshot {
                assert_eq!(&value, previous);
            } else {
                snapshot = Some(value);
            }
        }
    }
}

fn source_step_config() -> ArticulatedContactStepConfig {
    ArticulatedContactStepConfig {
        duration: 0.125,
        angular_origin_world: Vec3::new(3.0, -2.0, 1.0),
        contact: FrictionalContactConfig::default(),
        drift_energy_tolerance: 1e-8,
        energy_tolerance: 1e-8,
        linear_momentum_tolerance: 1e-8,
        angular_momentum_tolerance: 1e-8,
    }
}

fn source_step_states(velocity: Vec3) -> [StepState; 2] {
    let twist = Twist::new(Vec3::default(), velocity);
    [
        StepState::Free(FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(Se3::identity(), twist),
            joint_position: vec![0.0],
            joint_velocity: vec![4.0],
        }),
        StepState::Prescribed(PrescribedBaseMotionState {
            time_s: 0.0,
            base: BaseState::prescribed(Se3::identity(), twist, Twist::zero()),
            joint_position: vec![0.0],
            joint_velocity: vec![0.0],
        }),
    ]
}

fn source_step_contact() -> ArticulatedImpactGeometry {
    let mut c = contact(17, 0);
    c.other = ContactCounterpart::Body(LinkContactPoint {
        system_index: 1,
        link_index: 1,
        point_body: Vec3::default(),
    });
    c.normal_law = NormalVelocityLaw::Target { speed: 0.0 };
    c
}

fn source_step_support(
    index: usize,
    state: &PrescribedBaseMotionState,
    duration: f64,
) -> Result<PrescribedBasePath, PrescribedBaseMotionError> {
    assert_eq!(index, 1);
    assert_eq!(state.base.world_from_base.rotation(), So3::identity());
    assert_eq!(state.base.twist_body.angular, Vec3::default());
    // Authored uniform translation, not interpolation of unknown drive data.
    let sample = |fraction| {
        BaseState::prescribed(
            Se3::from_parts(
                So3::identity(),
                state.base.world_from_base.translation()
                    + state.base.twist_body.linear.scale(fraction * duration),
            )
            .unwrap(),
            state.base.twist_body,
            Twist::zero(),
        )
    };
    Ok(PrescribedBasePath {
        midpoint: sample(0.5),
        endpoint: sample(1.0),
    })
}

fn source_step_path(path: &MixedArticulatedContactPath<'_, ()>) -> Result<(), &'static str> {
    // Both original models have only a fixed root and a Cartesian x slider.
    // With zero angular speed and this uniform support, the unloaded native
    // drift has constant joint rates and affine witnesses. Nonnegative endpoint
    // gaps therefore admit the complete path; this fixture has no other obstacles.
    for (before, after) in path.post_contact.iter().zip(path.endpoint) {
        let (a, b, qa, qb, qda, qdb) = match (before, after) {
            (StepState::Free(a), StepState::Free(b)) => (
                a.base.twist_body,
                b.base.twist_body,
                &a.joint_position,
                &b.joint_position,
                &a.joint_velocity,
                &b.joint_velocity,
            ),
            (StepState::Prescribed(a), StepState::Prescribed(b)) => {
                assert_eq!(a.base.acceleration_body, Twist::zero());
                assert_eq!(b.base.acceleration_body, Twist::zero());
                (
                    a.base.twist_body,
                    b.base.twist_body,
                    &a.joint_position,
                    &b.joint_position,
                    &a.joint_velocity,
                    &b.joint_velocity,
                )
            }
            _ => return Err("changed boundary kind"),
        };
        close(a.angular.norm(), 0.0);
        close(b.angular.norm(), 0.0);
        close((a.linear - b.linear).norm(), 0.0);
        close(qda[0], qdb[0]);
        close(qb[0] - qa[0], path.duration * qda[0]);
    }
    Ok(())
}

#[test]
fn source_force_contact_drift_uses_postkick_speed_and_advances_both_original_trees() {
    let (card, law, query) = source_fixture(curve("speed", SPEED, [(0.0, 0.1), (12.0, 0.4)]));
    let model = slider();
    let wrenches = [Wrench::default(); 2];
    let effort = [12.0];
    let zero = [0.0];
    for drive in [Vec3::new(0.5, -0.25, 0.5), Vec3::new(16.5, -8.25, 4.5)] {
        let states = source_step_states(drive);
        let before = states.clone();
        let histories = [law.initial_state(), law.initial_state()];
        let original_histories = histories.clone();
        let systems = [
            StepSystem {
                model: &model,
                state: &states[0],
                loads: FreeFloatingMotionLoads {
                    generalized_force: &effort,
                    external_body_wrench: &wrenches,
                    gravity_world: Vec3::new(0.0, 0.0, -8.0),
                },
            },
            StepSystem {
                model: &model,
                state: &states[1],
                loads: FreeFloatingMotionLoads {
                    generalized_force: &zero,
                    external_body_wrench: &wrenches,
                    gravity_world: Vec3::default(),
                },
            },
        ];
        let mut inactive = contact(3, 0);
        inactive.point.link_index = 0;
        inactive.other = ContactCounterpart::PrescribedPlane {
            point_world: Vec3::new(0.0, 0.0, -10.0),
            velocity_world: Vec3::default(),
        };
        let mut calls = 0;
        let candidate = step_source_mixed_articulated_elliptic_contacts(
            &systems,
            &[source_step_contact(), inactive],
            &histories,
            &[binding(&law); 2],
            0.0,
            source_step_config(),
            |accepted, k| {
                calls += 1;
                assert_eq!(accepted, &histories[0]);
                assert_eq!((k.contact_id, k.original_index), (17, 0));
                for (actual, expected) in k.relative_velocity_m_s.into_iter().zip([-1.0, 4.5, 0.0])
                {
                    close(actual, expected);
                }
                Ok::<_, &'static str>(conditions(&query))
            },
            source_step_support,
            source_step_path,
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(calls, 1);
        assert_eq!(candidate.sources.len(), 1);
        let sample = &candidate.sources[0].sample;
        close(sample.queried_speed_m_s, 4.5); // Original pre-load slip was 4.
        close(sample.coefficients_mu[0], 0.2125);
        close(sample.coefficients_mu[1], 0.4);
        assert_eq!(sample.state.constitutive_state().slip_m(), [0.0; 2]);
        for source in sample.state.friction_sources() {
            close(source.state_coordinate("path").unwrap(), 0.25);
            close(source.state_coordinate("speed").unwrap(), 4.5);
        }
        original_claims(&card, &histories[0], &sample.state);
        let step = &candidate.mechanical;
        assert_eq!(step.contact.contact_ids, [3, 17]);
        assert_eq!(step.contact.contact.impulses[0], [0.0; 3]);
        // Native load kick: J_joint=1.5; free root recoils by -1.5/2,
        // child gains 1.5/3, and the original relative joint rate gains 1.25.
        // Normal recoil moves 5 kg; tangent contact moves each 3 kg slider.
        let pn = 5.0;
        let pt = -0.2125 * pn;
        for (actual, expected) in step.contact.contact.impulses[1]
            .into_iter()
            .zip([pn, pt, 0.0])
        {
            close(actual, expected);
        }
        let StepState::Free(upper) = &step.states_after[0] else {
            panic!("free owner lost")
        };
        let StepState::Prescribed(lower) = &step.states_after[1] else {
            panic!("support owner lost")
        };
        close(upper.base.twist_body.linear.x, drive.x - 0.75);
        close(upper.base.twist_body.linear.z, drive.z);
        close(upper.joint_velocity[0], 5.25 + pt / 3.0);
        close(lower.joint_velocity[0], -pt / 3.0);
        close(upper.joint_position[0], 0.125 * (5.25 + pt / 3.0));
        close(lower.joint_position[0], -0.125 * pt / 3.0);
        close(
            (lower.base.world_from_base.translation() - drive.scale(0.125)).norm(),
            0.0,
        );
        close(lower.time_s, 0.125);
        close(step.time_before_s, 0.0);
        close(step.time_after_s, 0.125);
        let Response::Prescribed(reaction) = &step.contact.responses[1] else {
            panic!("support response lost")
        };
        close(
            (reaction.support_impulse_world.force - Vec3::new(0.0, 0.0, pn)).norm(),
            0.0,
        );
        close(step.contact.plane_work, 0.0);
        close(step.contact.boundary_work, pn * drive.z);
        close(step.contact.support_work, pn * drive.z);
        close(step.actuator_impulse_work, 1.5 * (4.0 + 5.25) / 2.0);
        close(step.gravity_impulse_work, -pn * (drive.z - 0.5));
        close(step.load_support_work, 0.0);
        close(step.drift_support_work, 0.0);
        close(
            step.contact.contact.dissipated_energy,
            0.5 * pn - 0.5 * pt * (9.0 + 2.0 * pt / 3.0),
        );
        close(step.kinetic_work_residual, 0.0);
        close(step.drift_energy_defect, 0.0);
        close(step.linear_momentum_residual.norm(), 0.0);
        close(step.angular_momentum_residual.norm(), 0.0);
        close(step.external_linear_impulse_world.norm(), 0.0);
        assert!(
            step.initial_gaps
                .iter()
                .chain(&step.endpoint_gaps)
                .all(|g| g.signed_gap >= 0.0)
        );
        assert_eq!(states, before);
        assert_eq!(histories, original_histories);
    }
}

#[test]
fn source_contact_step_keeps_histories_private_through_path_refusal_and_final_cancellation() {
    let (_, law, query) = source_fixture(curve("speed", SPEED, [(0.0, 0.1), (12.0, 0.4)]));
    let model = slider();
    let states = source_step_states(Vec3::new(0.5, -0.25, 0.5));
    let original_states = states.clone();
    let histories = [law.initial_state()];
    let original_histories = histories.clone();
    let wrenches = [Wrench::default(); 2];
    let effort = [12.0];
    let zero = [0.0];
    let systems = [
        StepSystem {
            model: &model,
            state: &states[0],
            loads: FreeFloatingMotionLoads {
                generalized_force: &effort,
                external_body_wrench: &wrenches,
                gravity_world: Vec3::new(0.0, 0.0, -8.0),
            },
        },
        StepSystem {
            model: &model,
            state: &states[1],
            loads: FreeFloatingMotionLoads {
                generalized_force: &zero,
                external_body_wrench: &wrenches,
                gravity_world: Vec3::default(),
            },
        },
    ];
    let mut retry = None;
    for mode in 0..6 {
        let gate = CancelGate::new();
        if mode == 0 {
            gate.request();
        }
        let mut source_calls = 0;
        let mut path_calls = 0;
        let result = step_source_mixed_articulated_elliptic_contacts(
            &systems,
            &[source_step_contact()],
            &histories,
            &[binding(&law)],
            0.0,
            source_step_config(),
            |accepted, _| {
                source_calls += 1;
                assert_eq!(accepted, &histories[0]);
                if mode == 1 {
                    return Err("source unavailable");
                }
                Ok(conditions(&query))
            },
            source_step_support,
            |path| {
                path_calls += 1;
                source_step_path(path)?;
                if mode == 2 {
                    return Err("whole path refused");
                }
                if mode == 3 {
                    gate.request();
                }
                Ok(())
            },
            &gate,
        );
        assert_eq!(states, original_states);
        assert_eq!(histories, original_histories);
        assert_eq!(source_calls, usize::from(mode != 0));
        assert_eq!(path_calls, usize::from(mode >= 2));
        match mode {
            0 | 3 => assert!(matches!(
                result,
                Err(MixedArticulatedEllipticSourceStepError::Step(
                    ArticulatedContactStepError::Cancelled
                ))
            )),
            1 => assert!(matches!(
                result,
                Err(MixedArticulatedEllipticSourceStepError::Source {
                    contact_id: 17,
                    original_index: 0,
                    source: SourceRigidImpactError::Conditions("source unavailable"),
                })
            )),
            2 => assert!(matches!(
                result,
                Err(MixedArticulatedEllipticSourceStepError::Step(
                    ArticulatedContactStepError::Path("whole path refused")
                ))
            )),
            _ => {
                let report = result.unwrap();
                assert_eq!(report.sources.len(), 1);
                assert_ne!(report.sources[0].sample.state, histories[0]);
                let snapshot = (
                    report.mechanical.states_after,
                    report.mechanical.contact.contact.impulses,
                    report.sources[0].sample.state.clone(),
                );
                if let Some(previous) = &retry {
                    assert_eq!(&snapshot, previous);
                } else {
                    retry = Some(snapshot);
                }
            }
        }
    }
}
