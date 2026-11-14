use super::super::SystemScene;
use super::*;
use fs_mbd::{Pose, UnitQuaternion, Vec3};

const PAIR: &str = include_str!("../../../../../examples/system/damped_pair.json");
const MOUNT: &str = include_str!("../../../../../examples/system/viscoelastic_mount.json");

fn close(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "{actual:?} != {expected:?}"
    );
}

fn rotation(angle: f64) -> [f64; 4] {
    let (sine, cosine) = fs_math::det::sin_cos(0.5 * angle);
    let quaternion = [cosine, sine, 0.0, 0.0];
    if cosine < 0.0 {
        quaternion.map(|x| -x)
    } else {
        quaternion
    }
}

#[test]
fn full_vector_errors_detect_equal_magnitudes_and_every_signed_branch_amplitude() {
    let (error, scale) = vector_error([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]);
    close(error, hypot(1.0, 1.0), 0.0);
    assert_eq!(scale, 1.0);

    // Each original factor row matters, including a row with no stiffness.
    // A sign reversal has equal energy but the opposite constitutive state.
    for row in 0..6 {
        let mut amplitude = [0.0; 6];
        amplitude[row] = 0.75;
        let (error, scale) = vector_error(amplitude, amplitude.map(|x| -x), amplitude);
        assert_eq!(error, 1.5);
        assert_eq!(scale, 0.75);
        let mut turned = [0.0; 6];
        turned[(row + 1) % 6] = 0.75;
        close(
            vector_error(amplitude, turned, amplitude).0,
            hypot(0.75, 0.75),
            0.0,
        );
    }
}

#[test]
fn vector_metrics_keep_representable_extremes_and_refuse_every_invalid_input() {
    let large = [1.0e200, -1.0e200, 1.0e200];
    let (error, scale) = vector_error(large, large, large);
    assert_eq!(error, 0.0);
    assert!(scale.is_finite() && scale > 1.0e200);
    let tiny = vector_error([0.0; 3], [1.0e-200, 0.0, 0.0], [0.0; 3]);
    assert_eq!(tiny, (1.0e-200, 0.0));

    for checkpoint in 0..3 {
        for row in 0..6 {
            for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let mut values = [[0.0; 6]; 3];
                values[checkpoint][row] = bad;
                let (error, scale) = vector_error(values[0], values[1], values[2]);
                assert!(error.is_nan() && scale.is_nan());
            }
        }
    }
    assert!(vector_error([0.0], [f64::MAX], [-f64::MAX]).0.is_nan());
    assert!(vector_error([f64::MAX; 3], [0.0; 3], [0.0; 3]).1.is_nan());
}

#[test]
fn rotation_metric_resolves_small_angles_norm_drift_and_canonical_half_turns() {
    let identity = [1.0, 0.0, 0.0, 0.0];
    let arbitrary = [0.5, 0.5, -0.5, 0.5];
    assert_eq!(rotation_distance(arbitrary, arbitrary), 0.0);
    assert_eq!(rotation_distance(arbitrary, arbitrary.map(|x| -x)), 0.0);
    let drift = [1.0 + 8.0 * f64::EPSILON, -0.0, 0.0, 0.0];
    let bits = drift.map(f64::to_bits);
    assert_eq!(rotation_distance(identity, drift), 0.0);
    assert_eq!(drift.map(f64::to_bits), bits);
    close(
        rotation_distance(identity, rotation(1.0e-12)),
        1.0e-12,
        1.0e-26,
    );
    close(
        rotation_distance(identity, [0.0, 1.0, 0.0, 0.0]),
        core::f64::consts::PI,
        2.0e-15,
    );

    let delta = 1.0e-7;
    let before = rotation(core::f64::consts::PI - delta);
    let after = rotation(core::f64::consts::PI + delta);
    assert!(before[1] > 0.0 && after[1] < 0.0);
    close(rotation_distance(before, after), 2.0 * delta, 2.0e-15);

    for invalid in [
        [0.0; 4],
        [f64::NAN, 0.0, 0.0, 0.0],
        [1.0, f64::INFINITY, 0.0, 0.0],
    ] {
        assert!(rotation_distance(invalid, identity).is_nan());
        assert!(rotation_distance(identity, invalid).is_nan());
    }
}

fn product(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    [
        a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
        a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
        a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
        a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
    ]
}

#[test]
fn rotation_error_and_interval_scale_ignore_quaternion_sign_and_world_frame() {
    let start = rotation(0.1);
    let coarse = rotation(0.17);
    let refined = rotation(0.16);
    let world = [0.5, -0.5, 0.5, 0.5];
    let expected = (
        rotation_distance(coarse, refined),
        rotation_distance(start, refined),
    );
    close(expected.0, 0.01, 2.0e-15);
    close(expected.1, 0.06, 2.0e-15);
    close(
        rotation_distance(product(world, coarse), product(world, refined)),
        expected.0,
        2.0e-15,
    );
    close(
        rotation_distance(product(world, start), product(world, refined)),
        expected.1,
        2.0e-15,
    );
    assert_eq!(rotation_distance(coarse.map(|x| -x), refined), expected.0);
    assert_eq!(rotation_distance(start, refined.map(|x| -x)), expected.1);
}

#[test]
fn body_vectors_use_mass_and_the_full_principal_frame_tensor() {
    let frame = UnitQuaternion::new(0.8, 0.0, 0.0, 0.6).unwrap();
    let properties =
        MassProperties::from_principal_frame(4.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0), frame)
            .unwrap();
    let angular_velocity = Vec3::new(0.25, -0.5, 0.75);
    let body = RigidBodyState::new(
        Pose::new(Vec3::new(1.0, -2.0, 3.0), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(8.0, -12.0, 16.0),
        properties
            .angular_momentum_body_checked(angular_velocity)
            .unwrap(),
    )
    .unwrap();
    assert_eq!(
        BodyQuantity::Position.vector(&body, properties),
        Some([1.0, -2.0, 3.0])
    );
    assert_eq!(
        BodyQuantity::Velocity.vector(&body, properties),
        Some([2.0, -3.0, 4.0])
    );
    let spin = BodyQuantity::AngularVelocity
        .vector(&body, properties)
        .unwrap();
    for (actual, expected) in spin.into_iter().zip([0.25, -0.5, 0.75]) {
        close(actual, expected, 2.0e-15);
    }
    let tiny = MassProperties::new(
        1.0e-300,
        Vec3::ZERO,
        Vec3::new(1.0e-300, 1.0e-300, 1.0e-300),
    )
    .unwrap();
    let extreme = RigidBodyState::new(
        Pose::identity(),
        Vec3::new(1.0e100, 0.0, 0.0),
        Vec3::new(1.0e100, 0.0, 0.0),
    )
    .unwrap();
    assert!(BodyQuantity::Velocity.vector(&extreme, tiny).is_none());
    assert!(
        BodyQuantity::AngularVelocity
            .vector(&extreme, tiny)
            .is_none()
    );
}

fn captured(source: &str, name: &str) -> (Checkpoint, NodeId, MassProperties) {
    let scene = SystemScene::from_text(source, "adaptive_metric.json").unwrap();
    let element = scene
        .device
        .elements
        .iter()
        .find(|element| element.name == name)
        .unwrap();
    (
        scene.device.system.checkpoint().unwrap(),
        element.node,
        element.bodies[0],
    )
}

#[test]
fn native_checkpoint_dispatch_keeps_vector_direction_memory_sign_and_missing_accounts_visible() {
    let (initial, node, properties) = captured(PAIR, "pair");
    let (turned, _, _) = captured(
        &PAIR.replace(
            "\"velocity_m_s\": [1.0, 0.0, 0.0]",
            "\"velocity_m_s\": [0.0, 1.0, 0.0]",
        ),
        "pair",
    );
    let coordinate = NetworkCoordinate::Body {
        body: 0,
        properties,
        quantity: BodyQuantity::Velocity,
    };
    close(
        coordinate
            .estimate(node, &initial, &turned, &initial)
            .unwrap()
            .0,
        hypot(1.0, 1.0),
        0.0,
    );

    let (memory, mount, _) = captured(MOUNT, "mount");
    let (reversed, _, _) = captured(
        &MOUNT.replace(
            "\"extension_m\": [0.5, 0.0, 0.0]",
            "\"extension_m\": [-0.5, 0.0, 0.0]",
        ),
        "mount",
    );
    let branch = NetworkCoordinate::Branch {
        connector: 0,
        branch: 0,
    };
    let (error, scale) = branch.estimate(mount, &memory, &reversed, &memory).unwrap();
    assert!(scale > 0.0);
    close(error, 2.0 * scale, 1.0e-15);

    for missing in [
        NetworkCoordinate::Body {
            body: 2,
            properties,
            quantity: BodyQuantity::Position,
        },
        NetworkCoordinate::Branch {
            connector: 1,
            branch: 0,
        },
        NetworkCoordinate::Branch {
            connector: 0,
            branch: 1,
        },
    ] {
        let mut norm = ErrorNorm::new();
        missing.include(
            &mut norm,
            mount,
            &memory,
            &memory,
            &memory,
            ErrorTolerance {
                absolute: 1.0,
                relative: 0.0,
            },
        );
        assert!(norm.value().is_nan());
    }
    assert!(
        branch
            .estimate(node, &initial, &initial, &initial)
            .is_none()
    );
}
