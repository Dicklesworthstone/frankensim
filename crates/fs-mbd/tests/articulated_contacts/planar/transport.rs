//! G1/G3/G4: objective history transport through original articulated mechanics.

use super::*;
use fs_mbd::articulated_contacts::tangential::planar::transport::{
    ArticulatedFramedPlanarTangentialResult, ArticulatedFramedTangentialSite,
    advance_articulated_framed_planar_tangential, advance_articulated_framed_tangential_contacts,
};
use fs_tribo::tangential::transport::{
    FramedTangentialState, ObjectiveTangentialLaw, TangentialTransportConfig,
    TangentialTransportError,
};

fn transport_controls() -> TangentialTransportConfig {
    TangentialTransportConfig {
        absolute_energy_tolerance_j: 1e-13,
        relative_energy_tolerance: 1e-12,
    }
}

fn site<'a, L: ObjectiveTangentialLaw>(
    contact: ArticulatedPlanarTangentialContact,
    law: &'a L,
    history: &'a FramedTangentialState<L::State>,
    dt: f64,
    spin: f64,
) -> ArticulatedFramedTangentialSite<'a, L> {
    ArticulatedFramedTangentialSite {
        contact,
        law,
        history,
        drive: drive(dt),
        common_spin_angle_rad: spin,
    }
}

fn rotate2(value: [f64; 2], angle: f64) -> [f64; 2] {
    let (s, c) = angle.sin_cos();
    [c * value[0] - s * value[1], s * value[0] + c * value[1]]
}

fn world(frame: TangentFrame, value: [f64; 2]) -> Vec3 {
    let a = frame.longitudinal();
    let b = frame.lateral();
    Vec3::new(a[0], a[1], a[2]).scale(value[0]) + Vec3::new(b[0], b[1], b[2]).scale(value[1])
}

fn rotate_frame(frame: TangentFrame, rotation: So3) -> TangentFrame {
    let n = frame.normal();
    let t = frame.longitudinal();
    let n = rotation.rotate(Vec3::new(n[0], n[1], n[2])).unwrap();
    let t = rotation.rotate(Vec3::new(t[0], t[1], t[2])).unwrap();
    TangentFrame::new([n.x, n.y, n.z], [t.x, t.y, t.z]).unwrap()
}

fn audit_framed<S>(result: &ArticulatedFramedPlanarTangentialResult<S>) {
    audit(&result.mechanics);
    close(result.transport_energy_residual_j, 0.0);
    close(result.accepted_contact_work_residual_j, 0.0);
    close(result.accepted_energy_balance_residual_j, 0.0);
}

#[test]
fn tangent_gauge_rotates_elastic_and_plastic_coordinates_without_changing_world_recoil() {
    let model = body(2.0, 1.0);
    let initial = system(
        &model,
        Vec3::default(),
        Vec3::new(0.6, -0.2, 0.1),
        Vec3::new(0.2, -0.3, 0.4),
    );
    let law = layer(0.5);
    let history = FramedTangentialState::new(&law, preload(&law, [0.01, -0.02]), frame()).unwrap();
    let saved = history.clone();
    assert!(history.state().plastic_slip_m()[0].abs() > 0.001);
    let contact = port(Vec3::new(0.6, -0.8, 0.2), Vec3::new(0.1, 0.2, 0.0));
    let baseline = advance_articulated_framed_planar_tangential(
        &[initial],
        site(contact, &law, &history, 0.02, 0.0),
        config(),
        transport_controls(),
        || false,
    )
    .unwrap();
    let angle = 0.71;
    let gauged = ArticulatedPlanarTangentialContact {
        frame_world: frame().rotated(angle).unwrap(),
        ..contact
    };
    let changed = advance_articulated_framed_planar_tangential(
        &[initial],
        site(gauged, &law, &history, 0.02, 0.0),
        config(),
        transport_controls(),
        || false,
    )
    .unwrap();
    let a = &baseline.mechanics;
    let b = &changed.mechanics;
    close_vector(
        world(contact.frame_world, a.contact.traction_n),
        world(gauged.frame_world, b.contact.traction_n),
    );
    close2(b.contact.traction_n, rotate2(a.contact.traction_n, -angle));
    close2(
        b.contact.state.state().elastic_slip_m(),
        rotate2(a.contact.state.state().elastic_slip_m(), -angle),
    );
    close2(
        b.contact.state.state().plastic_slip_m(),
        rotate2(a.contact.state.state().plastic_slip_m(), -angle),
    );
    close_vector(
        b.responses[0].base_after.twist_body.linear,
        a.responses[0].base_after.twist_body.linear,
    );
    close_vector(
        b.responses[0].base_after.twist_body.angular,
        a.responses[0].base_after.twist_body.angular,
    );
    close(b.contact.dissipated_j, a.contact.dissipated_j);
    close(b.contact_energy_after_j, a.contact_energy_after_j);
    assert_eq!(b.contact.state.frame(), gauged.frame_world);
    assert_eq!(history, saved);
    audit_framed(&baseline);
    audit_framed(&changed);
}

#[test]
fn normal_tilt_keeps_preload_energy_and_releases_the_rotated_physical_wrench() {
    let model = body(2.0, 1.0);
    let initial = system(&model, Vec3::default(), Vec3::default(), Vec3::default());
    let law = layer(100.0);
    let elastic = [0.01, -0.007];
    let history = FramedTangentialState::new(&law, preload(&law, elastic), frame()).unwrap();
    // The minimum rotation sends old y to (0,.8,-.6), without shortening it.
    let contact = ArticulatedPlanarTangentialContact {
        frame_world: TangentFrame::new([0.0, 0.6, 0.8], [1.0, 0.0, 0.0]).unwrap(),
        ..port(Vec3::new(0.4, -0.3, 0.2), Vec3::default())
    };
    let dt = 0.03;
    let result = advance_articulated_framed_planar_tangential(
        &[initial],
        site(contact, &law, &history, dt, 0.0),
        config(),
        transport_controls(),
        || false,
    )
    .unwrap();
    let r = &result.mechanics;
    let gram = mobility(&[contact], 0.5, 1.0);
    let expected = elastic_force(
        [[gram[0], gram[1]], [gram[2], gram[3]]],
        dt,
        elastic,
        [0.0; 2],
    );
    close2(r.contact.traction_n, expected);
    let impulse = world(contact.frame_world, expected).scale(dt);
    close_vector(
        r.responses[0].base_after.twist_body.linear,
        impulse.scale(0.5),
    );
    close_vector(
        r.responses[0].base_after.twist_body.angular,
        contact.point.point_body.cross(impulse),
    );
    close(
        r.contact_energy_before_j,
        200.0 * (elastic[0] * elastic[0] + elastic[1] * elastic[1]),
    );
    for i in 0..2 {
        close(
            r.contact.state.state().elastic_slip_m()[i],
            elastic[i] + r.slip_increment_m[i],
        );
    }
    close2(r.contact.state.state().plastic_slip_m(), [0.0; 2]);
    close(r.contact.dissipated_j, 0.0);
    assert!(r.kinetic_energy_change_j > 0.0);
    close(
        r.kinetic_energy_change_j + law.stored_energy_j(r.contact.state.state()),
        law.stored_energy_j(history.state()),
    );
    audit_framed(&result);
}

#[test]
fn rotating_the_whole_world_preserves_tilt_spin_history_and_actual_body_response() {
    let model = body(2.0, 1.0);
    let law = layer(0.8);
    let state = preload(&law, [0.014, -0.009]);
    let arm = Vec3::new(0.4, -0.3, 0.2);
    let support = Vec3::new(-0.1, 0.2, 0.3);
    let old_frame = frame().rotated(-0.2).unwrap();
    let endpoint = TangentFrame::new([0.0, 0.6, 0.8], [1.0, 0.0, 0.0])
        .unwrap()
        .rotated(0.3)
        .unwrap();
    let mut reference = None;
    for rotation in [
        So3::identity(),
        So3::try_from_quat(Quat {
            w: 0.5,
            x: -0.5,
            y: 0.5,
            z: 0.5,
        })
        .unwrap(),
    ] {
        let initial = ArticulatedContactSystem {
            model: &model,
            base: FreeFloatingBaseState::new(
                Se3::from_parts(rotation, Vec3::default()).unwrap(),
                Twist::new(Vec3::new(0.2, -0.3, 0.4), Vec3::new(0.6, -0.2, 0.1)),
            ),
            joint_position: &[],
            joint_velocity: &[],
        };
        let history =
            FramedTangentialState::new(&law, state.clone(), rotate_frame(old_frame, rotation))
                .unwrap();
        let contact = ArticulatedPlanarTangentialContact {
            other: ContactCounterpart::PrescribedPlane {
                point_world: rotation.rotate(arm).unwrap(),
                velocity_world: rotation.rotate(support).unwrap(),
            },
            frame_world: rotate_frame(endpoint, rotation),
            ..port(arm, Vec3::default())
        };
        let result = advance_articulated_framed_planar_tangential(
            &[initial],
            site(contact, &law, &history, 0.02, 0.43),
            config(),
            transport_controls(),
            || false,
        )
        .unwrap();
        audit_framed(&result);
        let r = result.mechanics;
        if let Some((twist, force, elastic, plastic, heat, support_work)) = reference {
            let twist: Twist = twist;
            close_vector(r.responses[0].base_after.twist_body.linear, twist.linear);
            close_vector(r.responses[0].base_after.twist_body.angular, twist.angular);
            close2(r.contact.traction_n, force);
            close2(r.contact.state.state().elastic_slip_m(), elastic);
            close2(r.contact.state.state().plastic_slip_m(), plastic);
            close(r.contact.dissipated_j, heat);
            close(r.support_work_j, support_work);
        } else {
            reference = Some((
                r.responses[0].base_after.twist_body,
                r.contact.traction_n,
                r.contact.state.state().elastic_slip_m(),
                r.contact.state.state().plastic_slip_m(),
                r.contact.dissipated_j,
                r.support_work_j,
            ));
        }
    }
}

#[test]
fn common_spin_rotates_retained_vectors_without_inventing_slip_or_heat() {
    let model = body(2.0, 1.0);
    let initial = system(&model, Vec3::default(), Vec3::default(), Vec3::default());
    let law = FrictionLaw::Coulomb {
        static_mu: 0.6,
        kinetic_mu: 0.4,
    };
    let accumulated = [0.02, -0.01];
    let history = FramedTangentialState::new(&law, preload(&law, accumulated), frame()).unwrap();
    let angle = 0.67;
    for gauge in [0.0, angle] {
        let contact = ArticulatedPlanarTangentialContact {
            frame_world: frame().rotated(gauge).unwrap(),
            ..port(Vec3::new(0.4, -0.3, 0.2), Vec3::default())
        };
        let result = advance_articulated_framed_planar_tangential(
            &[initial],
            site(contact, &law, &history, 0.02, angle),
            config(),
            transport_controls(),
            || false,
        )
        .unwrap();
        close2(
            result.mechanics.contact.state.state().slip_m(),
            rotate2(accumulated, angle - gauge),
        );
        close2(result.mechanics.contact.slip_increment_m, [0.0; 2]);
        close2(result.mechanics.contact.traction_n, [0.0; 2]);
        close(result.mechanics.contact.dissipated_j, 0.0);
        assert_eq!(result.mechanics.responses[0].base_after, initial.base);
        audit_framed(&result);
    }
    let layer = layer(100.0);
    let elastic = [0.01, -0.007];
    let memory = FramedTangentialState::new(&layer, preload(&layer, elastic), frame()).unwrap();
    let result = advance_articulated_framed_planar_tangential(
        &[initial],
        site(
            port(Vec3::default(), Vec3::default()),
            &layer,
            &memory,
            0.02,
            angle,
        ),
        config(),
        transport_controls(),
        || false,
    )
    .unwrap();
    close2(
        result.mechanics.contact.traction_n,
        elastic_force(
            [[0.5, 0.0], [0.0, 0.5]],
            0.02,
            rotate2(elastic, angle),
            [0.0; 2],
        ),
    );
    close(result.mechanics.contact.dissipated_j, 0.0);
    audit_framed(&result);
}

#[test]
fn different_site_transports_feed_one_full_cross_contact_elastic_solution() {
    let model = body(2.0, 1.0);
    let initial = system(
        &model,
        Vec3::default(),
        Vec3::new(0.2, -0.1, 0.05),
        Vec3::new(0.1, 0.2, -0.1),
    );
    let laws = layers();
    let mut contacts = opposing();
    let elastic = [[0.004, -0.002], [-0.003, 0.005]];
    let histories = [
        FramedTangentialState::new(
            &laws[0],
            preload(&laws[0], elastic[0]),
            contacts[0].frame_world,
        )
        .unwrap(),
        FramedTangentialState::new(
            &laws[1],
            preload(&laws[1], elastic[1]),
            contacts[1].frame_world,
        )
        .unwrap(),
    ];
    let gauge = [0.3, -0.2];
    let spin = [-0.5, 0.4];
    contacts[0].frame_world = TangentFrame::new([0.0, 0.6, 0.8], [1.0, 0.0, 0.0])
        .unwrap()
        .rotated(gauge[0])
        .unwrap();
    contacts[1].frame_world = TangentFrame::new([0.6, 0.0, -0.8], [0.8, 0.0, 0.6])
        .unwrap()
        .rotated(gauge[1])
        .unwrap();
    let dt = 0.02;
    let inputs = [
        site(contacts[0], &laws[0], &histories[0], dt, spin[0]),
        site(contacts[1], &laws[1], &histories[1], dt, spin[1]),
    ];
    let result = advance_articulated_framed_tangential_contacts(
        &[initial],
        &inputs,
        controls(),
        transport_controls(),
        || false,
    )
    .unwrap();
    let r = &result.mechanics;
    let w = mobility(&contacts, 0.5, 1.0);
    assert!(w[2].abs() > 0.01 || w[3].abs() > 0.01);
    let transported = [
        rotate2(elastic[0], spin[0] - gauge[0]),
        rotate2(elastic[1], spin[1] - gauge[1]),
    ];
    let force = elastic_reference(
        &w,
        &released(&contacts, initial.base.twist_body, dt),
        &[400.0, 700.0],
        &transported,
        dt,
    );
    for i in 0..2 {
        close2(r.contacts[i].traction_n, [force[2 * i], force[2 * i + 1]]);
        close(
            r.contact_energy_before_j[i],
            laws[i].stored_energy_j(histories[i].state()),
        );
        for j in 0..2 {
            close(
                r.contacts[i].state.state().elastic_slip_m()[j],
                transported[i][j] + r.slip_increment_m[i][j],
            );
        }
        assert_eq!(r.contacts[i].state.frame(), contacts[i].frame_world);
        close(result.transport_energy_residuals_j[i], 0.0);
        close(result.accepted_contact_work_residuals_j[i], 0.0);
    }
    let (impulse, torque, work) = resultant(&contacts, &force, dt);
    close_vector(
        r.responses[0].base_after.twist_body.linear,
        initial.base.twist_body.linear + impulse.scale(0.5),
    );
    close_vector(
        r.responses[0].base_after.twist_body.angular,
        initial.base.twist_body.angular + torque,
    );
    close(r.support_work_j, work);
    close(result.accepted_energy_balance_residual_j, 0.0);
    audit_all(r);
}

#[test]
fn frame_source_and_late_mechanics_refusals_leave_all_accepted_histories_retryable() {
    let model = body(2.0, 1.0);
    let initial = system(
        &model,
        Vec3::default(),
        Vec3::new(0.7, -0.4, 0.1),
        Vec3::new(0.3, -0.2, 0.15),
    );
    let laws = layers();
    let mut contacts = opposing();
    let histories = [
        FramedTangentialState::new(
            &laws[0],
            preload(&laws[0], [0.004, -0.002]),
            contacts[0].frame_world,
        )
        .unwrap(),
        FramedTangentialState::new(
            &laws[1],
            preload(&laws[1], [-0.003, 0.005]),
            contacts[1].frame_world,
        )
        .unwrap(),
    ];
    let saved = histories.clone();
    contacts[0].frame_world = contacts[0].frame_world.rotated(0.3).unwrap();
    contacts[1].frame_world = contacts[1].frame_world.rotated(-0.4).unwrap();
    let mut inputs = [
        site(contacts[0], &laws[0], &histories[0], 0.02, -0.2),
        site(contacts[1], &laws[1], &histories[1], 0.02, 0.5),
    ];
    let mut polls = 0;
    let baseline = advance_articulated_framed_tangential_contacts(
        &[initial],
        &inputs,
        controls(),
        transport_controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        let refused = advance_articulated_framed_tangential_contacts(
            &[initial],
            &inputs,
            controls(),
            transport_controls(),
            || {
                seen += 1;
                seen == stop
            },
        );
        assert!(matches!(refused, Err(ArticulatedContactError::Cancelled)));
        assert_eq!(histories, saved);
    }
    let mut tight = controls();
    tight.mechanical.maximum_slip_speed_m_s = 1e-6;
    assert!(matches!(
        advance_articulated_framed_tangential_contacts(
            &[initial],
            &inputs,
            tight,
            transport_controls(),
            || false
        ),
        Err(ArticulatedContactError::Unresolved(_))
    ));
    let foreign = layer(100.0);
    inputs[1].law = &foreign;
    assert!(matches!(
        advance_articulated_framed_tangential_contacts(
            &[initial],
            &inputs,
            controls(),
            transport_controls(),
            || false
        ),
        Err(ArticulatedContactError::TangentialTransport(
            TangentialTransportError::Law(_)
        ))
    ));
    inputs[1].law = &laws[1];
    inputs[1].contact.frame_world = frame(); // antipodal to its accepted normal
    assert!(matches!(
        advance_articulated_framed_tangential_contacts(
            &[initial],
            &inputs,
            controls(),
            transport_controls(),
            || false
        ),
        Err(ArticulatedContactError::TangentialTransport(
            TangentialTransportError::Frame(_)
        ))
    ));
    inputs[1].contact = contacts[1];
    let mut small = controls();
    small.coupling.maximum_contacts = 1;
    assert!(matches!(
        advance_articulated_framed_tangential_contacts(
            &[initial],
            &inputs,
            small,
            transport_controls(),
            || false
        ),
        Err(ArticulatedContactError::Invalid(_))
    ));
    assert_eq!(histories, saved);
    let retry = advance_articulated_framed_tangential_contacts(
        &[initial],
        &inputs,
        controls(),
        transport_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(retry.mechanics.contacts, baseline.mechanics.contacts);
    assert_eq!(
        retry.mechanics.responses[0].base_after,
        baseline.mechanics.responses[0].base_after
    );
    assert_eq!(
        retry.transport_energy_residuals_j,
        baseline.transport_energy_residuals_j
    );
    assert_eq!(histories, saved);
}
