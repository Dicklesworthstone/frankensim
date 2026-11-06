//! Ordered source cards, physical slip geometry and whole-path applicability.
//! All constitutive values are synthetic fixtures, not measured material data.
use fs_contact::tangential::impact::multiple::CoupledFrictionImpulseControls;
use fs_contact::tangential::impact::multiple::relaxation::CoupledFrictionFlowControls;
use fs_mbd::rigid_rattle::distance::DistanceStepControls;
use fs_mbd::rigid_rattle::distance::graph::joint_friction::relaxation::interface::{
    InterfaceJointFrictionError, JointDryInterface, JointDrySlipGeometry,
    advance_interface_joint_friction,
};
use fs_mbd::rigid_rattle::distance::graph::joint_friction::relaxation::{
    JointFrictionFlowControls, JointFrictionFlowMotion, JointFrictionFlowPort,
};
use fs_mbd::rigid_rattle::distance::graph::{
    ArticulatedJoint, DistanceGraphControls, PrismaticJoint, RevoluteJoint,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::{
    ApplicabilityRange, DryFrictionApplicability, DryInterfaceSystemCard, FrictionLaw,
    InputAuthority, InterfaceMedium, InterfaceSystemRef, TriboError,
};

const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-10 * a.abs().max(b.abs()).max(1.0),
        "{a:e} != {b:e}"
    );
}
fn properties(inertia: f64) -> MassProperties {
    MassProperties::new(2.0, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap()
}
fn state(p: MassProperties, velocity: f64, spin: f64) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).unwrap(),
        Z.scale(velocity * p.mass()),
        p.angular_momentum_body_checked(Z.scale(spin)).unwrap(),
    )
    .unwrap()
}
fn controls() -> JointFrictionFlowControls {
    JointFrictionFlowControls {
        graph: DistanceGraphControls {
            step: DistanceStepControls {
                maximum_impulse_ns: 1000.0,
                maximum_rotation_increment_rad: 0.25,
                maximum_position_iterations: 30,
                maximum_backtracks: 16,
            },
            position_tolerance_m: 1e-10,
            velocity_tolerance_mps: 1e-10,
            minimum_scaled_pivot: 1e-10,
            maximum_bodies: 32,
            maximum_constraints: 32,
            maximum_trial_evaluations: 160,
            maximum_live_bytes: 2_000_000,
        },
        flow: CoupledFrictionFlowControls {
            projection: CoupledFrictionImpulseControls {
                maximum_iterations: 256,
                impulse_tolerance_ns: 1e-10,
                velocity_tolerance_m_s: 1e-10,
                energy_tolerance_j: 1e-9,
                relative_matrix_tolerance: 1e-12,
                maximum_impulse_ns: 1000.0,
            },
            maximum_segments: 32,
            acceleration_tolerance_m_s2: 1e-10,
            maximum_force_n: 1000.0,
        },
        work_tolerance_j: 1e-9,
        linear_impulse_tolerance_ns: 1e-9,
        angular_impulse_tolerance_nms: 1e-9,
    }
}
fn card(
    id: &str,
    authority: InputAuthority,
    mu: [f64; 2],
    speed: [f64; 2],
) -> DryInterfaceSystemCard {
    DryInterfaceSystemCard::new(
        InterfaceSystemRef::new(
            format!("{id}:A->B"),
            format!("{id}:history"),
            format!("{id}:synthetic-source"),
            authority,
            InterfaceMedium::Dry,
        )
        .unwrap(),
        FrictionLaw::Coulomb {
            static_mu: mu[0],
            kinetic_mu: mu[1],
        },
        DryFrictionApplicability::new(
            ApplicabilityRange::new(250.0, 450.0).unwrap(),
            ApplicabilityRange::new(0.0, 1e6).unwrap(),
            ApplicabilityRange::new(speed[0], speed[1]).unwrap(),
        )
        .unwrap(),
    )
    .unwrap()
}
fn slider(
    v: f64,
) -> (
    [RigidBodyState; 2],
    [MassProperties; 2],
    [ArticulatedJoint; 1],
) {
    let p = [properties(2.0); 2];
    (
        [state(p[0], v, 0.0), state(p[1], 0.0, 0.0)],
        p,
        [ArticulatedJoint::Prismatic(
            PrismaticJoint::new([0, 1], [Vec3::ZERO; 2], [UnitQuaternion::IDENTITY; 2], 1.0)
                .unwrap(),
        )],
    )
}
fn sliding<'a>(card: &'a DryInterfaceSystemCard, drive_n: f64) -> JointDryInterface<'a> {
    JointDryInterface {
        joint_index: 0,
        card,
        geometry: JointDrySlipGeometry::Prismatic { drive_n },
        normal_force_n: 2.0,
        nominal_area_m2: 1e-4,
        temperature_k: 300.0,
    }
}
fn hinge(bodies: [usize; 2], length: f64) -> ArticulatedJoint {
    ArticulatedJoint::Revolute(RevoluteJoint::new(bodies, [Vec3::ZERO; 2], [Z; 2], length).unwrap())
}
fn rotors() -> (
    [RigidBodyState; 3],
    [MassProperties; 3],
    [ArticulatedJoint; 2],
) {
    let p = [properties(10.0), properties(1.0 / 0.9), properties(10.0)];
    (
        [
            state(p[0], 0.0, 0.0),
            state(p[1], 0.0, 0.1),
            state(p[2], 0.0, 1.01),
        ],
        p,
        [hinge([0, 1], 1.0), hinge([1, 2], 1.0)],
    )
}
fn rotor_bindings<'a>(cards: [&'a DryInterfaceSystemCard; 2]) -> [JointDryInterface<'a>; 2] {
    core::array::from_fn(|i| JointDryInterface {
        joint_index: i,
        card: cards[i],
        geometry: JointDrySlipGeometry::CylindricalSleeve {
            radius_m: 1.0,
            drive_nm: 0.0,
        },
        normal_force_n: 1.0,
        nominal_area_m2: 0.01,
        temperature_k: 300.0,
    })
}

#[test]
fn same_ordered_card_holds_at_rest_and_uses_lower_kinetic_capacity_in_motion() {
    let card = card("rail", InputAuthority::Estimated, [0.8, 0.3], [0.0, 3.0]);
    for (speed, end, travel, loss) in [(0.0, 0.0, 0.0, 0.0), (1.0, 1.4, 1.2, 0.72)] {
        let (s, p, j) = slider(speed);
        let out = advance_interface_joint_friction(
            &s,
            &p,
            &j,
            &[sliding(&card, 1.0)],
            1.0,
            controls(),
            || false,
        )
        .unwrap();
        let JointFrictionFlowMotion::Prismatic {
            velocity_after_m_s,
            displacement_m,
            ..
        } = out.mechanical.ports[0].motion
        else {
            panic!("slider")
        };
        near(velocity_after_m_s, end);
        near(displacement_m, travel);
        near(out.mechanical.dissipated_j, loss);
        near(out.mechanical.drive_work_j, travel);
        near(out.mechanical.kinetic_energy_change_j, travel - loss);
        let receipt = &out.interfaces[0];
        assert!(core::ptr::eq(receipt.binding.card, &card));
        assert_eq!(receipt.ordered_bodies, [0, 1]);
        near(receipt.nominal_pressure_pa, 20_000.0);
        assert_eq!(
            receipt.binding.card.interface().provenance().authority(),
            InputAuthority::Estimated
        );
        assert_eq!(
            receipt.binding.card.interface().history_id(),
            "rail:history"
        );
        let JointFrictionFlowPort::PrismaticStiction {
            static_capacity_n,
            kinetic_capacity_n,
            ..
        } = receipt.source
        else {
            panic!("stiction source")
        };
        near(static_capacity_n, 1.6);
        near(kinetic_capacity_n, 0.6);
        near(receipt.admitted_state_range[0].slip_speed_mps(), speed);
        near(receipt.admitted_state_range[1].slip_speed_mps(), end);
    }
}

#[test]
fn sleeve_radius_controls_torque_and_speed_independently_of_native_hinge_scaling() {
    let card = card(
        "sleeve",
        InputAuthority::CallerDeclared,
        [0.8, 0.3],
        [0.0, 1.0],
    );
    let p = [properties(2.0); 2];
    let s = [state(p[0], 0.0, 0.0), state(p[1], 0.0, 1.0)];
    let binding = JointDryInterface {
        joint_index: 0,
        card: &card,
        geometry: JointDrySlipGeometry::CylindricalSleeve {
            radius_m: 0.2,
            drive_nm: 0.0,
        },
        normal_force_n: 10.0,
        nominal_area_m2: 2e-5,
        temperature_k: 300.0,
    };
    let mut endpoint = None;
    for length in [0.125, 1.0, 8.0] {
        let out = advance_interface_joint_friction(
            &s,
            &p,
            &[hinge([0, 1], length)],
            &[binding],
            0.5,
            controls(),
            || false,
        )
        .unwrap();
        let JointFrictionFlowMotion::Revolute {
            rate_after_rad_s,
            rotation_rad,
            ..
        } = out.mechanical.ports[0].motion
        else {
            panic!("sleeve")
        };
        near(rate_after_rad_s, 0.7);
        near(rotation_rad, 0.425);
        near(out.mechanical.dissipated_j, 0.255);
        near(out.interfaces[0].nominal_pressure_pa, 500_000.0);
        near(
            out.interfaces[0].admitted_state_range[0].slip_speed_mps(),
            0.14,
        );
        near(
            out.interfaces[0].admitted_state_range[1].slip_speed_mps(),
            0.2,
        );
        let JointFrictionFlowPort::RevoluteStiction {
            static_capacity_nm,
            kinetic_capacity_nm,
            ..
        } = out.interfaces[0].source
        else {
            panic!("sleeve source")
        };
        near(static_capacity_nm, 1.6);
        near(kinetic_capacity_nm, 0.6);
        if let Some(previous) = endpoint {
            assert_eq!(out.mechanical.states_after, previous);
        }
        endpoint = Some(out.mechanical.states_after);
    }
}

#[test]
fn original_pressure_temperature_and_lower_speed_domains_refuse_without_extrapolation() {
    let card = card(
        "bounded",
        InputAuthority::SyntheticFixture,
        [0.8, 0.3],
        [0.25, 2.0],
    );
    let (s, p, j) = slider(1.0);
    let original = s;
    let binding = sliding(&card, 0.0);
    let mut hot = binding;
    hot.temperature_k = 451.0;
    let mut high_pressure = binding;
    high_pressure.nominal_area_m2 = 1e-7;
    for (binding, duration, field) in [
        (hot, 0.1, "temperature_kelvin"),
        (high_pressure, 0.1, "contact_pressure_pa"),
        (binding, 2.0, "slip_speed_mps"),
    ] {
        let error =
            advance_interface_joint_friction(&s, &p, &j, &[binding], duration, controls(), || {
                false
            })
            .unwrap_err();
        assert!(
            matches!(error,InterfaceJointFrictionError::Constitutive {joint_index:0,
            error:TriboError::OutsideApplicability {field:actual,..}} if actual==field),
            "{error:?}"
        );
        assert_eq!(s, original);
    }
}

#[test]
fn coupled_interior_speed_peak_must_fit_the_original_source_even_when_outer_endpoints_do() {
    let wide = card(
        "accelerated",
        InputAuthority::Estimated,
        [0.2, 0.1],
        [0.0, 2.0],
    );
    let narrow = card(
        "accelerated",
        InputAuthority::Estimated,
        [0.2, 0.1],
        [0.0, 0.5],
    );
    let neighbour = card(
        "neighbour",
        InputAuthority::SyntheticFixture,
        [1.5, 1.0],
        [0.0, 2.0],
    );
    let (s, p, j) = rotors();
    let original = s;
    let out = advance_interface_joint_friction(
        &s,
        &p,
        &j,
        &rotor_bindings([&wide, &neighbour]),
        60.0,
        controls(),
        || false,
    )
    .unwrap();
    // W=[[1,-.9],[-.9,1]], v=[.1,.91]. Before neighbour arrest at t=1,
    // dv0/dt=.8, so v0 peaks at .9. It then decelerates to zero.
    let JointFrictionFlowMotion::Revolute {
        rate_before_rad_s,
        rate_after_rad_s,
        ..
    } = out.mechanical.ports[0].motion
    else {
        panic!("hinge")
    };
    near(rate_before_rad_s, 0.1);
    near(rate_after_rad_s, 0.0);
    near(
        out.interfaces[0].admitted_state_range[0].slip_speed_mps(),
        0.0,
    );
    near(
        out.interfaces[0].admitted_state_range[1].slip_speed_mps(),
        0.9,
    );
    assert!(rate_before_rad_s.abs() < 0.5 && rate_after_rad_s.abs() < 0.5);
    let error = advance_interface_joint_friction(
        &s,
        &p,
        &j,
        &rotor_bindings([&narrow, &neighbour]),
        60.0,
        controls(),
        || false,
    )
    .unwrap_err();
    assert!(
        matches!(error,InterfaceJointFrictionError::Constitutive {joint_index:0,
        error:TriboError::OutsideApplicability {field:"slip_speed_mps",value,maximum:0.5,..}} if value>0.8),
        "{error:?}"
    );
    assert_eq!(s, original);
}

#[test]
fn binding_permutation_retains_ordered_cards_history_authority_and_native_motion() {
    let cards = [
        card(
            "first",
            InputAuthority::CallerDeclared,
            [0.2, 0.1],
            [0.0, 2.0],
        ),
        card("second", InputAuthority::Estimated, [1.5, 1.0], [0.0, 2.0]),
    ];
    let (s, p, j) = rotors();
    let bindings = rotor_bindings([&cards[0], &cards[1]]);
    let a =
        advance_interface_joint_friction(&s, &p, &j, &bindings, 2.0, controls(), || false).unwrap();
    let b = advance_interface_joint_friction(
        &s,
        &p,
        &j,
        &[bindings[1], bindings[0]],
        2.0,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(a.mechanical, b.mechanical);
    for i in 0..2 {
        let x = &a.interfaces[i];
        let y = &b.interfaces[i];
        assert_eq!(x.binding.joint_index, i);
        assert_eq!(y.binding.joint_index, i);
        assert!(
            core::ptr::eq(x.binding.card, &cards[i]) && core::ptr::eq(y.binding.card, &cards[i])
        );
        assert_eq!(x.ordered_bodies, [i, i + 1]);
        assert_eq!(x.ordered_bodies, y.ordered_bodies);
        assert_eq!(x.source, y.source);
        assert_eq!(x.admitted_state_range, y.admitted_state_range);
        assert_eq!(
            x.binding.card.interface().provenance(),
            cards[i].interface().provenance()
        );
        assert_eq!(
            x.binding.card.interface().ordered_system_id(),
            cards[i].interface().ordered_system_id()
        );
        assert_eq!(
            x.binding.card.interface().history_id(),
            cards[i].interface().history_id()
        );
    }
}

#[test]
fn unsupported_laws_geometry_and_mid_solve_cancellation_publish_no_candidate() {
    let card = card(
        "law",
        InputAuthority::SyntheticFixture,
        [0.8, 0.3],
        [0.0, 3.0],
    );
    let stribeck = DryInterfaceSystemCard::new(
        card.interface().clone(),
        FrictionLaw::Stribeck {
            static_mu: 0.8,
            kinetic_mu: 0.3,
            characteristic_speed: 0.1,
            viscous_per_speed: 0.0,
        },
        card.applicability(),
    )
    .unwrap();
    let (s, p, j) = slider(1.0);
    let original = s;
    let error = advance_interface_joint_friction(
        &s,
        &p,
        &j,
        &[sliding(&stribeck, 0.0)],
        0.1,
        controls(),
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        InterfaceJointFrictionError::UnsupportedLaw { joint_index: 0 }
    ));
    let binding = sliding(&card, 0.0);
    let mut wrong = binding;
    wrong.geometry = JointDrySlipGeometry::CylindricalSleeve {
        radius_m: 1.0,
        drive_nm: 0.0,
    };
    let mut zero_area = binding;
    zero_area.nominal_area_m2 = 0.0;
    for invalid in [vec![wrong], vec![zero_area], vec![binding, binding], vec![]] {
        assert!(matches!(
            advance_interface_joint_friction(&s, &p, &j, &invalid, 0.1, controls(), || false),
            Err(InterfaceJointFrictionError::Invalid(_))
        ));
    }
    let mut polls = 0;
    let stopped = advance_interface_joint_friction(&s, &p, &j, &[binding], 1.0, controls(), || {
        polls += 1;
        polls > 20
    });
    assert!(stopped.is_err());
    assert!(polls > 20);
    assert_eq!(s, original);
}
