use super::*;

fn identity(component: &impl Component) -> Vec<u64> {
    let mut words = Vec::new();
    component.identity_words(&mut words);
    words
}

fn state_bits(states: &[RigidBodyState]) -> Vec<u64> {
    states
        .iter()
        .flat_map(|state| {
            let position = state.pose().position_world();
            let q = state.pose().orientation().components();
            let linear = state.linear_momentum_world();
            let angular = state.angular_momentum_body();
            [
                position.x, position.y, position.z, q[0], q[1], q[2], q[3], linear.x, linear.y,
                linear.z, angular.x, angular.y, angular.z,
            ]
            .map(f64::to_bits)
        })
        .collect()
}

fn make(
    properties: Vec<MassProperties>,
    links: Vec<ViscousNetworkLink>,
    gravity: Gravity,
    config: ViscousNetworkMotionControls,
) -> ViscousNetwork {
    ViscousNetwork::new(properties, links, gravity, config).unwrap()
}

#[derive(Debug)]
struct IdentifiedLaw(f64);

impl DampingLaw for IdentifiedLaw {
    fn checkpoint_identity(&self) -> Option<[u8; 32]> {
        let bytes: Vec<_> = [self.0, 250.0, 400.0, 300.0]
            .into_iter()
            .flat_map(f64::to_le_bytes)
            .collect();
        Some(
            *fs_blake3::hash_domain("org.frankensim.tests.viscous-network.thinning.v1", &bytes)
                .as_bytes(),
        )
    }

    fn at_temperature(
        &self,
        connector: usize,
        reference: &PassiveBushingMatrix,
        temperature_k: f64,
    ) -> Result<PassiveBushingMatrix, String> {
        Thinning { beta_per_k: self.0 }.at_temperature(connector, reference, temperature_k)
    }
}

#[test]
fn rigid_body_and_network_codecs_keep_canonical_quaternion_and_signed_bits() {
    let (properties, links, mut states) = fixture();
    // Canonical transport admits this roundoff-sized norm drift. Re-running
    // the normalizing constructor would silently change its first component.
    let components = [1.0 + 8.0 * f64::EPSILON, -0.0, 0.0, -0.0];
    let orientation = UnitQuaternion::from_canonical_components(components).unwrap();
    assert_ne!(
        UnitQuaternion::new(components[0], components[1], components[2], components[3])
            .unwrap()
            .components()
            .map(f64::to_bits),
        components.map(f64::to_bits)
    );
    states[0] = RigidBodyState::new(
        Pose::new(Vec3::new(-0.0, 2.0, -3.0), orientation).unwrap(),
        Vec3::new(0.4, -0.0, -0.2),
        Vec3::new(-0.0, 0.3, 0.2),
    )
    .unwrap();
    let network = make(properties.clone(), links, Gravity::ZERO, controls());
    let body = RigidBody::new(RigidBodyIntegrator::new(Gravity::ZERO), properties[0]).unwrap();
    let mut body_words = Vec::new();
    assert!(body.encode_state(&states[0], &mut body_words));
    assert_eq!(body_words, state_bits(&states[..1]));
    let decoded = body.decode_state(&body_words).unwrap();
    assert_eq!(state_bits(&[decoded]), body_words);
    let mut words = Vec::new();
    assert!(network.encode_state(&states, &mut words));
    assert_eq!(words, state_bits(&states));
    assert_eq!(state_bits(&network.decode_state(&words).unwrap()), words);

    for (index, value) in [
        (0, f64::NAN),
        (3, 2.0),
        (3, -1.0),
        (7, f64::INFINITY),
        (10, f64::MAX),
        (38, f64::NAN),
    ] {
        let mut invalid = words.clone();
        invalid[index] = value.to_bits();
        assert!(network.decode_state(&invalid).is_none(), "word {index}");
        if index < 13 {
            assert!(body.decode_state(&invalid[..13]).is_none());
        }
    }
    let mut zero_rotation = words.clone();
    zero_rotation[3..7].fill(0);
    assert!(network.decode_state(&zero_rotation).is_none());
    assert!(network.decode_state(&words[..words.len() - 1]).is_none());
    let mut extra = words.clone();
    extra.push(0);
    assert!(network.decode_state(&extra).is_none());
    assert!(body.decode_state(&body_words[..12]).is_none());
    assert!(body.decode_state(&words).is_none());

    let mut untouched = vec![17, 29];
    assert!(!network.encode_state(&states[..2].to_vec(), &mut untouched));
    assert_eq!(untouched, [17, 29]);
    let overflowing =
        RigidBodyState::new(Pose::identity(), Vec3::new(f64::MAX, 0.0, 0.0), Vec3::ZERO).unwrap();
    assert!(!body.encode_state(&overflowing, &mut untouched));
    states[2] = overflowing;
    assert!(!network.encode_state(&states, &mut untouched));
    assert_eq!(untouched, [17, 29]);
}

#[test]
fn body_bank_codecs_refuse_overflow_of_the_total_energy() {
    let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let state =
        RigidBodyState::new(Pose::identity(), Vec3::new(1.2e154, 0.0, 0.0), Vec3::ZERO).unwrap();
    assert!(
        RigidBodyIntegrator::new(Gravity::ZERO)
            .diagnostics(state, properties)
            .unwrap()
            .mechanical_energy
            .is_finite()
    );
    let network = make(vec![properties; 3], Vec::new(), Gravity::ZERO, controls());
    let mut words = vec![71];
    assert!(!network.encode_state(&vec![state; 3], &mut words));
    assert_eq!(words, [71]);
    assert!(network.decode_state(&state_bits(&[state; 3])).is_none());
}

#[test]
fn network_identity_binds_inertia_geometry_full_forms_controls_and_port_order() {
    let (properties, links, _) = fixture();
    let base = identity(&make(
        properties.clone(),
        links.clone(),
        Gravity::ZERO,
        controls(),
    ));
    let changed =
        |properties: Vec<MassProperties>, links: Vec<ViscousNetworkLink>, gravity, controls| {
            assert_ne!(identity(&make(properties, links, gravity, controls)), base);
        };
    for field in 0..6 {
        let mut p = properties.clone();
        let old = p[0];
        let mut moments = old.principal_inertia_body();
        match field {
            1 => moments.x += 0.1,
            2 => moments.y += 0.1,
            3 => moments.z += 0.1,
            _ => {}
        }
        p[0] = MassProperties::from_principal_frame(
            if field == 0 { 2.1 } else { old.mass() },
            if field == 4 {
                Vec3::new(-0.0, 0.0, 0.0)
            } else {
                Vec3::ZERO
            },
            moments,
            if field == 5 {
                UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.3).unwrap()
            } else {
                UnitQuaternion::IDENTITY
            },
        )
        .unwrap();
        changed(p, links.clone(), Gravity::ZERO, controls());
    }
    for axis in 0..3 {
        let mut v = [0.0; 3];
        v[axis] = 0.1;
        changed(
            properties.clone(),
            links.clone(),
            Gravity::new(Vec3::new(v[0], v[1], v[2])).unwrap(),
            controls(),
        );
    }
    for arm in 0..2 {
        for axis in 0..3 {
            let mut l = links.clone();
            let v = &mut l[0].arms_body_m[arm];
            match axis {
                0 => v.x += 0.1,
                1 => v.y += 0.1,
                _ => v.z += 0.1,
            }
            changed(properties.clone(), l, Gravity::ZERO, controls());
        }
    }
    let mut reordered = links.clone();
    reordered.reverse();
    changed(properties.clone(), reordered, Gravity::ZERO, controls());
    let mut moved = links.clone();
    moved[0].bodies = [0, 2];
    changed(properties.clone(), moved, Gravity::ZERO, controls());
    for variation in 0..4 {
        let mut l = links.clone();
        let mut native = *l[0].damping.native_coefficients();
        if variation == 0 {
            native[0][0] += 0.2;
        }
        if variation == 1 {
            native[0][1] = 0.01;
            native[1][0] = 0.01;
        }
        l[0].damping = PassiveBushingMatrix::try_new(
            native,
            BushingCoordinateScales {
                translation_m: if variation == 2 { 2.0 } else { 1.0 },
                rotation_rad: if variation == 3 { 3.0 } else { 1.0 },
            },
        )
        .unwrap();
        changed(properties.clone(), l, Gravity::ZERO, controls());
    }
    for field in 0..13 {
        let mut c = controls();
        match field {
            0 => c.source.impulse.angular_origin_world_m.x += 0.1,
            1 => c.source.impulse.angular_origin_world_m.y += 0.1,
            2 => c.source.impulse.angular_origin_world_m.z += 0.1,
            3 => c.source.impulse.linear_tolerance_ns *= 2.0,
            4 => c.source.impulse.angular_tolerance_nms *= 2.0,
            5 => c.source.impulse.work_tolerance_j *= 2.0,
            6 => c.source.linear_speed_tolerance_m_s *= 2.0,
            7 => c.source.angular_speed_tolerance_rad_s *= 2.0,
            8 => c.source_substeps += 1,
            9 => c.limits.maximum_bodies += 1,
            10 => c.limits.maximum_connectors += 1,
            11 => c.limits.maximum_pair_relaxations += 1,
            12 => c.limits.maximum_live_bytes += 1,
            _ => unreachable!(),
        }
        changed(properties.clone(), links.clone(), Gravity::ZERO, c);
    }
    let network = make(properties.clone(), links.clone(), Gravity::ZERO, controls());
    assert_ne!(
        identity(&network.clone().torqued(2).unwrap().torqued(0).unwrap()),
        identity(&network.torqued(0).unwrap().torqued(2).unwrap())
    );

    // Equal rounded matrices need not mean equal original source factors.
    let scales = BushingCoordinateScales {
        translation_m: 1.0,
        rotation_rad: 1.0,
    };
    let mut rows = [[0.0; 6]; 6];
    rows[0][0] = 1.0;
    rows[0][1] = 0.5;
    let original = PassiveBushingMatrix::from_gram_rows(rows, scales).unwrap();
    rows[0] = rows[0].map(|value| -value);
    let reversed = PassiveBushingMatrix::from_gram_rows(rows, scales).unwrap();
    assert_eq!(
        original.native_coefficients(),
        reversed.native_coefficients()
    );
    let mut left = links.clone();
    left[0].damping = original;
    let mut right = links.clone();
    right[0].damping = reversed;
    assert_ne!(
        identity(&make(properties.clone(), left, Gravity::ZERO, controls())),
        identity(&make(properties, right, Gravity::ZERO, controls()))
    );
}

#[test]
fn opaque_damping_laws_refuse_checkpoints_and_named_laws_replay_native_steps() {
    let (properties, links, mut states) = fixture();
    let gravity = Gravity::new(Vec3::new(0.0, -1.0, -9.81)).unwrap();
    let opaque = ViscousNetwork::temperature_dependent(
        properties.clone(),
        links.clone(),
        gravity,
        controls(),
        Rc::new(Thinning { beta_per_k: -0.002 }),
    )
    .unwrap();
    let mut untouched = vec![101];
    assert!(!opaque.encode_state(&states, &mut untouched));
    assert_eq!(untouched, [101]);
    assert!(opaque.decode_state(&state_bits(&states)).is_none());
    let identified = |beta| {
        ViscousNetwork::temperature_dependent(
            properties.clone(),
            links.clone(),
            gravity,
            controls(),
            Rc::new(IdentifiedLaw(beta)),
        )
        .unwrap()
        .torqued(2)
        .unwrap()
        .torqued(0)
        .unwrap()
    };
    let component = identified(-0.002);
    assert_ne!(identity(&component), identity(&identified(-0.003)));
    let peer = identified(-0.002);
    let inputs = [
        0.3, 0.0, -0.1, 0.0, 0.0, 0.0, -0.2, 0.4, 0.0, 305.0, 288.0, 310.0, 330.0, 0.02, -0.01,
        0.01, -0.01, 0.01, 0.0,
    ];
    let output_count: usize = component
        .ports()
        .iter()
        .map(|port| port.output_width())
        .sum();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for _ in 0..17 {
            let mut next = Vec::new();
            let proposal = component
                .propose(
                    &states,
                    0.001,
                    &inputs,
                    &mut next,
                    &mut vec![0.0; output_count],
                    cx,
                )
                .unwrap();
            assert!(matches!(proposal, Proposal::Step(_)));
            states = next;
        }
        let mut words = Vec::new();
        assert!(component.encode_state(&states, &mut words));
        let mut resumed = peer.decode_state(&words).unwrap();
        assert_eq!(state_bits(&resumed), state_bits(&states));
        for _ in 0..40 {
            let (mut next, mut next_resumed) = (Vec::new(), Vec::new());
            let (mut outputs, mut resumed_outputs) =
                (vec![0.0; output_count], vec![0.0; output_count]);
            let result = component
                .propose(&states, 0.001, &inputs, &mut next, &mut outputs, cx)
                .unwrap();
            let replay = peer
                .propose(
                    &resumed,
                    0.001,
                    &inputs,
                    &mut next_resumed,
                    &mut resumed_outputs,
                    cx,
                )
                .unwrap();
            assert!(matches!(result, Proposal::Step(_)));
            assert_eq!(state_bits(&next), state_bits(&next_resumed));
            assert_eq!(
                outputs
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>(),
                resumed_outputs
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>()
            );
            assert_eq!(format!("{result:?}"), format!("{replay:?}"));
            states = next;
            resumed = next_resumed;
        }
    });
}
