//! G1/G3/G4: original source coefficients through joint and floating-base recoil.
use super::sphere_plane::{original_claims, source_fixture};
use super::*;
use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_contact::interface_binding::thermal_friction::elliptic::SourceFrictionConditions;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::SourceImpulseCoefficientRole;
use fs_couple::friction_heat::directional::articulated::*;
use fs_couple::friction_heat::directional::impact::SourceRigidImpactBinding;
use fs_ga::{Mat3, Se3, So3, Twist, Vec3};
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
