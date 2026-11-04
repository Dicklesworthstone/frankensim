//! G1/G3: original rigid bodies under simultaneous normal/friction impulses.
use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_exec::CancelGate;
use fs_mbd::frictional_impact::*;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

#[path = "frictional_impact/redundant.rs"]
mod redundant;

#[path = "frictional_impact/normal_gauge.rs"]
mod normal_gauge;

#[path = "frictional_impact/elliptic.rs"]
mod elliptic;

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
fn g3_material_gap_survives_a_distant_common_origin() {
    let properties = [property(1.0, Vec3::new(1.0, 1.0, 1.0)); 2];
    for origin in [Vec3::ZERO, Vec3::new(0.0, 0.0, 2.0_f64.powi(54))] {
        for gap in [-1.0, 1.0] {
            for dynamic in [false, true] {
                let initial = [
                    body(origin, UnitQuaternion::IDENTITY, Vec3::new(0.0, 0.0, -1.0)),
                    body(origin, UnitQuaternion::IDENTITY, Vec3::ZERO),
                ];
                let mut bodies = initial;
                let mut c = plane(
                    RigidImpactPoint {
                        body: 0,
                        arm_body: Vec3::new(0.0, 0.0, gap),
                    },
                    0.0,
                    0.0,
                );
                c.other = if dynamic {
                    RigidImpactCounterpart::Body(RigidImpactPoint {
                        body: 1,
                        arm_body: Vec3::ZERO,
                    })
                } else {
                    RigidImpactCounterpart::PrescribedPlane {
                        point_world: origin,
                        velocity_world: Vec3::ZERO,
                    }
                };
                let result = apply_rigid_frictional_impacts(
                    &mut bodies,
                    &properties,
                    &[c],
                    controls(),
                    &CancelGate::new(),
                );
                if gap < 0.0 {
                    assert!(matches!(
                        result,
                        Err(RigidFrictionalImpactError::Contact(
                            fs_contact::frictional_multicontact::FrictionalContactError::PenetratingContact { .. }
                        ))
                    ), "origin={origin:?}, dynamic={dynamic}: {result:?}");
                } else {
                    let report = result.unwrap();
                    assert_eq!(report.contact.impulses, [[0.0; 3]]);
                    assert_eq!(report.kinetic_change_j, 0.0);
                }
                // A real gap cannot become an impact, and genuine penetration
                // must refuse, even when COM + arm rounds to the COM itself.
                assert_eq!(bodies, initial);
            }
        }
    }
}

#[test]
fn g3_distant_origin_cannot_hide_a_tangential_witness_mismatch() {
    let properties = [property(1.0, Vec3::new(1.0, 1.0, 1.0)); 2];
    for origin in [Vec3::ZERO, Vec3::new(2.0_f64.powi(54), 0.0, 0.0)] {
        let initial = [
            body(origin, UnitQuaternion::IDENTITY, Vec3::new(0.0, 0.0, -1.0)),
            body(origin, UnitQuaternion::IDENTITY, Vec3::ZERO),
        ];
        let mut bodies = initial;
        let mut c = plane(
            RigidImpactPoint {
                body: 0,
                arm_body: Vec3::new(1.0, 0.0, 0.0),
            },
            0.0,
            0.0,
        );
        c.other = RigidImpactCounterpart::Body(RigidImpactPoint {
            body: 1,
            arm_body: Vec3::ZERO,
        });
        assert!(matches!(
            apply_rigid_frictional_impacts(
                &mut bodies,
                &properties,
                &[c],
                controls(),
                &CancelGate::new(),
            ),
            Err(RigidFrictionalImpactError::Invalid(
                "dynamic contact witnesses do not share a normal line"
            ))
        ));
        assert_eq!(bodies, initial);
    }
}

#[test]
fn g3_distant_origin_retains_the_original_support_lever_arm() {
    let properties = [property(1.0, Vec3::new(1.0, 1.0, 1.0))];
    for origin in [Vec3::ZERO, Vec3::new(2.0_f64.powi(54), 0.0, 0.0)] {
        let initial = body(origin, UnitQuaternion::IDENTITY, Vec3::new(0.0, 0.0, -1.0));
        let mut bodies = [initial];
        let mut c = plane(
            RigidImpactPoint {
                body: 0,
                arm_body: Vec3::new(1.0, 0.0, 0.0),
            },
            0.0,
            0.0,
        );
        c.other = RigidImpactCounterpart::PrescribedPlane {
            point_world: origin,
            velocity_world: Vec3::ZERO,
        };
        let report = apply_rigid_frictional_impacts(
            &mut bodies,
            &properties,
            &[c],
            controls(),
            &CancelGate::new(),
        )
        .unwrap_or_else(|error| panic!("origin={origin:?}: {error:?}"));
        // m = I_y = arm_x = 1 gives effective mass 1/2. Plastic impact
        // transfers equal translational and rotational kinetic energy.
        close(report.contact.impulses[0][0], 0.5);
        close_vector(bodies[0].linear_momentum_world(), Vec3::new(0.0, 0.0, -0.5));
        close_vector(bodies[0].angular_momentum_body(), Vec3::new(0.0, -0.5, 0.0));
        close(report.contact.dissipated_energy, 0.25);
        close_vector(report.angular_momentum_residual_nms, Vec3::ZERO);
        assert_eq!(bodies[0].pose(), initial.pose());
    }
}

#[test]
fn g3_origin_shift_cannot_hide_an_unbalanced_friction_couple() {
    let properties = [property(1.0, Vec3::new(1.0, 1.0, 1.0)); 2];
    for origin_x in [0.0, 2.0_f64.powi(54)] {
        let initial = [
            body(
                Vec3::new(origin_x, 0.0, 1.0),
                UnitQuaternion::IDENTITY,
                Vec3::new(1.0, 0.0, -1.0),
            ),
            body(
                Vec3::new(origin_x, 0.0, 0.0),
                UnitQuaternion::IDENTITY,
                Vec3::ZERO,
            ),
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
        c.activation_tolerance_m = 2.0;
        let result = apply_rigid_frictional_impacts(
            &mut bodies,
            &properties,
            &[c],
            controls(),
            &CancelGate::new(),
        );
        // Equal/opposite friction impulses at distinct normal-separated
        // witnesses create a real -0.25 N m s couple. Translating their common
        // origin cannot authorize this violation of angular momentum.
        let Err(RigidFrictionalImpactError::Audit(report)) = result else {
            panic!("origin_x={origin_x}: {result:?}");
        };
        close(report.angular_momentum_residual_nms.y, -0.25);
        assert_eq!(bodies, initial);
    }
}

#[test]
fn g3_contact_balance_retains_small_normal_loss_under_boost_and_origin_shift() {
    let properties = [property(1.0, Vec3::new(1.0, 1.0, 1.0))];
    for boost in [0.0, 1e6] {
        for origin in [Vec3::ZERO, Vec3::new(0.7, -0.4, 1e9)] {
            for restitution in [0.0, 0.5, 1.0] {
                let initial = body(
                    origin,
                    UnitQuaternion::IDENTITY,
                    Vec3::new(boost, 0.0, -0.3),
                );
                let mut bodies = [initial];
                let mut contact = plane(
                    RigidImpactPoint {
                        body: 0,
                        arm_body: Vec3::ZERO,
                    },
                    0.0,
                    restitution,
                );
                contact.other = RigidImpactCounterpart::PrescribedPlane {
                    point_world: origin,
                    velocity_world: Vec3::new(boost, 0.0, 0.0),
                };
                let report = apply_rigid_frictional_impacts(
                    &mut bodies,
                    &properties,
                    &[contact],
                    controls(),
                    &CancelGate::new(),
                )
                .unwrap_or_else(|error| panic!("boost={boost}, origin={origin:?}: {error:?}"));
                let loss = 0.045 * (1.0 - restitution * restitution);
                close(bodies[0].linear_momentum_world().z, 0.3 * restitution);
                assert_eq!(bodies[0].linear_momentum_world().x, boost);
                assert_eq!(bodies[0].angular_momentum_body(), Vec3::ZERO);
                assert_eq!(bodies[0].pose(), initial.pose());
                close(report.contact.dissipated_energy, loss);
                close(report.kinetic_change_j, -loss);
                close(report.work_energy_residual_j, 0.0);
                close(report.total_energy_residual_j, 0.0);
                close_vector(report.linear_momentum_residual_ns, Vec3::ZERO);
                close_vector(report.angular_momentum_residual_nms, Vec3::ZERO);
            }
        }
    }
}

#[test]
fn g0_simultaneous_opposing_impulses_retain_existing_central_momentum() {
    // The centre's original small momentum is representable. Updating it once
    // per contact loses it alongside either O(mass) impulse; their RESULTANT
    // must instead be applied to the original state once.
    let mass = 2.0_f64.powi(20);
    let small = 2.0_f64.powi(-34);
    let properties = [property(mass, Vec3::new(mass, mass, mass)); 3];
    let initial = [
        body(
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::new(mass, 0.0, 0.0),
        ),
        body(
            Vec3::new(2.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(small, 0.0, 0.0),
        ),
        body(
            Vec3::new(4.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(-mass, 0.0, 0.0),
        ),
    ];
    let contacts = std::array::from_fn::<_, 2, _>(|i| RigidFrictionalContact {
        id: 10 + i as u64,
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
        friction_coefficient: 0.0,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    });
    let mut bodies = initial;
    let report = apply_rigid_frictional_impacts(
        &mut bodies,
        &properties,
        &contacts,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let left = report.contact.impulses[0][0];
    let right = report.contact.impulses[1][0];
    let expected = small + (left - right);
    assert_eq!(bodies[1].linear_momentum_world().x, expected);
    assert_ne!(
        (small + left) - right,
        expected,
        "fixture must expose sequential rounding"
    );
    assert_eq!(report.linear_momentum_residual_ns, Vec3::ZERO);
    assert_eq!(report.angular_momentum_residual_nms, Vec3::ZERO);

    let mut reordered = initial;
    let repeated = apply_rigid_frictional_impacts(
        &mut reordered,
        &properties,
        &[contacts[1], contacts[0]],
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, reordered);
    assert_eq!(report.contact, repeated.contact);
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
fn g1_strongly_coupled_friction_resolves_native_momenta_beyond_block_contraction() {
    let properties = [property(5.0, Vec3::new(1.25, 1.25, 1.25))];
    let mut bodies = [body(
        Vec3::new(0.0, 0.0, 1.0),
        UnitQuaternion::IDENTITY,
        Vec3::new(10.0, 0.0, -5.0),
    )];
    let before = bodies;
    let point = RigidImpactPoint {
        body: 0,
        arm_body: Vec3::new(1.0, 0.0, -1.0),
    };
    let report = apply_rigid_frictional_impacts(
        &mut bodies,
        &properties,
        &[plane(point, 0.8, 0.0)],
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    // Wnn=Wtt=1/m+1/I=1; Wnt=1/I=.8, all from actual body inertia.
    // Simultaneous sliding: pn+.8pt=1 and pt=-.8pn. The old block
    // contraction check refused this unique response despite positive energy.
    close(report.contact.impulses[0][0], 25.0 / 9.0);
    close(report.contact.impulses[0][1], -20.0 / 9.0);
    close_vector(
        bodies[0].linear_momentum_world(),
        Vec3::new(70.0 / 9.0, 0.0, -20.0 / 9.0),
    );
    close_vector(
        bodies[0].angular_momentum_body(),
        Vec3::new(0.0, -5.0 / 9.0, 0.0),
    );
    assert_eq!(bodies[0].pose(), before[0].pose());
    close(report.contact.relative_velocities_after[0][0], 0.0);
    close(report.contact.relative_velocities_after[0][1], 2.0);
    close(report.kinetic_before_j, 12.5);
    close(report.kinetic_after_j, 20.0 / 3.0);
    close(report.contact.dissipated_energy, 35.0 / 6.0);
    close(report.total_energy_residual_j, 0.0);
    close_vector(report.linear_momentum_residual_ns, Vec3::ZERO);
    close_vector(report.angular_momentum_residual_nms, Vec3::ZERO);
    assert!(report.contact.impulse_energy_witness.is_some());
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
