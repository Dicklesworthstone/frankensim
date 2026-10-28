//! G1/G3: original rigid bodies under simultaneous normal/friction impulses.
use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_exec::CancelGate;
use fs_mbd::frictional_impact::*;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

#[path = "frictional_impact/redundant.rs"]
mod redundant;

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-8 * a.abs().max(b.abs()).max(1.0),
        "{a:e} != {b:e}"
    );
}
fn close_vector(a: Vec3, b: Vec3) {
    close(a.x, b.x);
    close(a.y, b.y);
    close(a.z, b.z);
}
fn controls() -> RigidImpactControls {
    RigidImpactControls {
        maximum_bodies: 8,
        maximum_contacts: 8,
        geometry_tolerance_m: 1e-10,
        velocity_tolerance_mps: 2e-8,
        linear_momentum_tolerance_ns: 1e-9,
        angular_momentum_tolerance_nms: 1e-9,
        energy_tolerance_j: 2e-8,
        contact: FrictionalContactConfig::default(),
    }
}
fn property(mass: f64, inertia: Vec3) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, inertia).unwrap()
}
fn body(position: Vec3, orientation: UnitQuaternion, momentum: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, orientation).unwrap(),
        momentum,
        Vec3::ZERO,
    )
    .unwrap()
}
fn plane(point: RigidImpactPoint, friction: f64, restitution: f64) -> RigidFrictionalContact {
    RigidFrictionalContact {
        id: 1,
        point,
        other: RigidImpactCounterpart::PrescribedPlane {
            point_world: Vec3::ZERO,
            velocity_world: Vec3::ZERO,
        },
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance_m: 1e-10,
        friction_coefficient: friction,
        normal_law: NormalVelocityLaw::NewtonRestitution {
            coefficient: restitution,
        },
    }
}

#[test]
fn g1_off_centre_impact_retains_normal_tangent_coupling_and_rotated_inertia() {
    let properties = [property(2.0, Vec3::new(1.0, 2.0, 2.0))];
    for angle in [0.0, core::f64::consts::FRAC_PI_2] {
        let q = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), angle).unwrap();
        let mut bodies = [body(
            Vec3::new(0.0, 0.0, 1.0),
            q,
            q.rotate_body_to_world(Vec3::new(6.0, 0.0, -2.0)),
        )];
        let initial = bodies;
        let point = RigidImpactPoint {
            body: 0,
            arm_body: Vec3::new(1.0, 0.0, -1.0),
        };
        let mut contact = plane(point, 0.3, 0.0);
        contact.tangent_world = q.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
        let report = apply_rigid_frictional_impacts(
            &mut bodies,
            &properties,
            &[contact],
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        // Independent point mobility: Wnn=Wtt=1, Wnt=Wtn=1/2.
        // Sliding pt=-mu*pn and vn+=-1+pn+pt/2=0.
        let pn = 1.0 / (1.0 - 0.5 * 0.3);
        let pt = -0.3 * pn;
        close(report.contact.impulses[0][0], pn);
        close(report.contact.impulses[0][1], pt);
        close(report.contact.impulses[0][2], 0.0);
        close(report.contact.relative_velocities_after[0][0], 0.0);
        close(
            report.contact.relative_velocities_after[0][1],
            3.0 + 0.5 * pn + pt,
        );
        close_vector(
            bodies[0].linear_momentum_world(),
            initial[0]
                .linear_momentum_world()
                .add(q.rotate_body_to_world(Vec3::new(pt, 0.0, pn))),
        );
        close_vector(
            bodies[0].angular_momentum_body(),
            Vec3::new(0.0, -pt - pn, 0.0),
        );
        assert_eq!(bodies[0].pose(), initial[0].pose());
        close(
            report.kinetic_before_j - report.kinetic_after_j,
            report.contact.dissipated_energy,
        );
        close(report.work_energy_residual_j, 0.0);
        close_vector(report.angular_momentum_residual_nms, Vec3::ZERO);
        // Sequential normal then Coulomb clipping would leave vn+=-0.15.
        assert!((pn - 1.0).abs() > 0.1);
    }
}

#[test]
fn g1_sticking_restitution_and_moving_support_work_close_native_energy() {
    let properties = [property(2.0, Vec3::new(1.0, 1.0, 1.0))];
    let mut bodies = [body(
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(0.4, -0.2, -2.0),
    )];
    let point = RigidImpactPoint {
        body: 0,
        arm_body: Vec3::ZERO,
    };
    let report = apply_rigid_frictional_impacts(
        &mut bodies,
        &properties,
        &[plane(point, 0.8, 0.5)],
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    close_vector(bodies[0].linear_momentum_world(), Vec3::new(0.0, 0.0, 1.0));
    close(report.contact.impulses[0][0], 3.0);
    close(report.contact.impulses[0][1], -0.4);
    close(report.contact.impulses[0][2], 0.2);
    close(report.contact.dissipated_energy, 0.8);

    let properties = [property(1.0, Vec3::new(1.0, 1.0, 1.0))];
    let mut moving = [body(
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(0.0, 0.0, -1.0),
    )];
    let mut contact = plane(point, 2.0, 0.0);
    contact.other = RigidImpactCounterpart::PrescribedPlane {
        point_world: Vec3::ZERO,
        velocity_world: Vec3::new(3.0, 0.0, 0.0),
    };
    let report = apply_rigid_frictional_impacts(
        &mut moving,
        &properties,
        &[contact],
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    close_vector(moving[0].linear_momentum_world(), Vec3::new(2.0, 0.0, 0.0));
    close(report.support_work_j, 6.0);
    close(report.contact.dissipated_energy, 4.5);
    close(report.kinetic_after_j - report.kinetic_before_j, 1.5);
    close(report.total_energy_residual_j, 0.0);
}

#[test]
fn g1_three_body_chain_is_simultaneous_and_exactly_permutation_invariant() {
    let properties = [property(1.0, Vec3::new(1.0, 1.0, 1.0)); 3];
    let initial = [
        body(
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::new(2.0, 1.0, 0.0),
        ),
        body(
            Vec3::new(2.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
        ),
        body(
            Vec3::new(4.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
        ),
    ];
    let contacts = std::array::from_fn::<_, 2, _>(|i| RigidFrictionalContact {
        id: 10 + 10 * i as u64,
        point: RigidImpactPoint {
            body: i,
            arm_body: Vec3::new(1.0, 0.0, 0.0),
        },
        other: RigidImpactCounterpart::Body(RigidImpactPoint {
            body: i + 1,
            arm_body: Vec3::new(-1.0, 0.0, 0.0),
        }),
        normal_world: Vec3::new(-1.0, 0.0, 0.0),
        tangent_world: Vec3::new(0.0, 1.0, 0.0),
        activation_tolerance_m: 1e-10,
        friction_coefficient: 0.5,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    });
    let mut a = initial;
    let mut b = initial;
    let first = apply_rigid_frictional_impacts(
        &mut a,
        &properties,
        &contacts,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let second = apply_rigid_frictional_impacts(
        &mut b,
        &properties,
        &[contacts[1], contacts[0]],
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(a, b);
    assert_eq!(first.contact_ids, second.contact_ids);
    assert_eq!(first.contact, second.contact);
    // Normal block [[2,-1],[-1,2]] and initially [v1,v2]=[-2,0].
    close(first.contact.impulses[0][0], 4.0 / 3.0);
    close(first.contact.impulses[1][0], 2.0 / 3.0);
    // Tangent block is diag(4,4): middle-body translational and rotational
    // cross terms cancel for the two opposite unit arms.
    close(first.contact.impulses[0][1], -0.25);
    close(first.contact.impulses[1][1], 0.0);
    for body in a {
        close(body.linear_momentum_world().x, 2.0 / 3.0);
    }
    close(a[0].linear_momentum_world().y, 0.75);
    close(a[1].linear_momentum_world().y, 0.25);
    close(a[2].linear_momentum_world().y, 0.0);
    close_vector(first.linear_momentum_residual_ns, Vec3::ZERO);
    close_vector(first.angular_momentum_residual_nms, Vec3::ZERO);
    close(first.contact.dissipated_energy, 35.0 / 24.0);
}

#[test]
fn g3_penetration_redundancy_and_late_angular_refusal_preserve_bodies() {
    let properties = [property(1.0, Vec3::new(1.0, 1.0, 1.0)); 2];
    let initial = [
        body(
            Vec3::new(0.0, 0.0, 1e-5),
            UnitQuaternion::IDENTITY,
            Vec3::new(3.0, 0.0, -1.0),
        ),
        body(Vec3::ZERO, UnitQuaternion::IDENTITY, Vec3::ZERO),
    ];
    let mut bodies = initial;
    let mut c = plane(
        RigidImpactPoint {
            body: 0,
            arm_body: Vec3::ZERO,
        },
        0.5,
        0.0,
    );
    c.other = RigidImpactCounterpart::Body(RigidImpactPoint {
        body: 1,
        arm_body: Vec3::ZERO,
    });
    c.activation_tolerance_m = 1e-4;
    let mut config = controls();
    config.geometry_tolerance_m = 1e-4;
    config.angular_momentum_tolerance_nms = 1e-12;
    // An active finite gap and tangential impulse create a real couple; the
    // late whole-body check must refuse even though the local law was solved.
    assert!(matches!(
        apply_rigid_frictional_impacts(&mut bodies, &properties, &[c], config, &CancelGate::new()),
        Err(RigidFrictionalImpactError::Audit(_))
    ));
    assert_eq!(bodies, initial);

    c.normal_world = Vec3::new(0.0, 0.0, -1.0);
    assert!(
        apply_rigid_frictional_impacts(
            &mut bodies,
            &properties,
            &[c],
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(bodies, initial);
    c.normal_world = Vec3::new(0.0, 0.0, 1.0);
    let mut duplicate = c;
    duplicate.id = 2;
    assert!(
        apply_rigid_frictional_impacts(
            &mut bodies,
            &properties,
            &[c, duplicate],
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(bodies, initial);
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        apply_rigid_frictional_impacts(&mut bodies, &properties, &[c], controls(), &gate),
        Err(RigidFrictionalImpactError::Cancelled)
    ));
    assert_eq!(bodies, initial);
}
