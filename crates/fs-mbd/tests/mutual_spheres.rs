//! G1/G3/G4: exact finite-sphere gravity, independent volume integration and refusal.
use fs_mbd::gravity::GravityError;
use fs_mbd::gravity::mutual::MutualGravity;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const FRAME: u64 = 127;
const G: f64 = 0.7;
const MASSES: [f64; 2] = [2.0, 3.0];

fn norm(value: Vec3) -> f64 {
    value.x.hypot(value.y).hypot(value.z)
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance {tolerance:e}"
    );
}

fn vector_close(actual: Vec3, expected: Vec3, tolerance: f64) {
    close(norm(actual.sub(expected)), 0.0, tolerance);
}

fn sphere_mass(mass: f64, radius: f64) -> MassProperties {
    let moment = 0.4 * mass * radius * radius;
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(moment, moment, moment)).unwrap()
}

fn state(position: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        Vec3::ZERO,
        Vec3::ZERO,
    )
    .unwrap()
}

fn bodies(distance: f64) -> [RigidBodyState; 2] {
    [state(Vec3::ZERO), state(Vec3::new(distance, 0.0, 0.0))]
}

// Four-point Gaussian quadrature integrates degree <= 7 polynomials exactly.
// The shell-volume integrands below have degree <= 5 after the source surface
// and the target's full/cap shell regions are split. This integrates the
// original source potential over target volume, not the pair kernel formula.
fn integrate(lower: f64, upper: f64, f: impl Fn(f64) -> f64) -> f64 {
    if upper <= lower {
        return 0.0;
    }
    let root = (6.0_f64 / 5.0).sqrt();
    let nodes = [
        ((3.0 + 2.0 * root) / 7.0).sqrt(),
        ((3.0 - 2.0 * root) / 7.0).sqrt(),
    ];
    let weights = [
        (18.0 - 30.0_f64.sqrt()) / 36.0,
        (18.0 + 30.0_f64.sqrt()) / 36.0,
    ];
    let midpoint = 0.5 * (lower + upper);
    let half = 0.5 * (upper - lower);
    half * nodes
        .into_iter()
        .zip(weights)
        .map(|(node, weight)| weight * (f(midpoint - half * node) + f(midpoint + half * node)))
        .sum::<f64>()
}

fn volume_potential_oracle(radius_a: f64, radius_b: f64, distance: f64) -> f64 {
    let a = radius_a.max(radius_b);
    let b = radius_a.min(radius_b);
    let potential = |s: f64| {
        if s <= a {
            -(1.5 - 0.5 * (s / a).powi(2)) / a
        } else {
            -1.0 / s
        }
    };
    let full = if distance < b {
        integrate(0.0, b - distance, |s| {
            3.0 * s * s / b.powi(3) * potential(s)
        })
    } else {
        0.0
    };
    if distance == 0.0 {
        return G * MASSES[0] * MASSES[1] * full;
    }
    let lower = (distance - b).abs();
    let upper = distance + b;
    let cap = |s: f64| {
        let area_weight = 3.0 * s * (b * b - (s - distance).powi(2)) / (4.0 * b.powi(3) * distance);
        area_weight * potential(s)
    };
    let cap_integral = integrate(lower, a.min(upper).max(lower), cap)
        + integrate(a.max(lower).min(upper), upper, cap);
    G * MASSES[0] * MASSES[1] * (full + cap_integral)
}

#[test]
fn exterior_containment_and_coincident_centers_have_exact_physical_limits() {
    let radii = [1.5, 0.4];
    let properties = [
        sphere_mass(MASSES[0], radii[0]),
        sphere_mass(MASSES[1], radii[1]),
    ];
    let point = MutualGravity::new(FRAME, G, 0.0, 1).unwrap();
    let finite = point.uniform_spheres(&radii).unwrap();
    for distance in [radii[0] + radii[1], 3.0, 8.0] {
        let states = bodies(distance);
        assert_eq!(
            finite.load(FRAME, &states, &properties, || false).unwrap(),
            point.load(FRAME, &states, &properties, || false).unwrap()
        );
    }
    for distance in [0.0, 0.2, 0.8, radii[0] - radii[1]] {
        let load = finite
            .load(FRAME, &bodies(distance), &properties, || false)
            .unwrap();
        let product = G * MASSES[0] * MASSES[1];
        let expected_u = -product / (2.0 * radii[0])
            * (3.0 - (distance / radii[0]).powi(2) - 0.6 * (radii[1] / radii[0]).powi(2));
        close(load.potential_energy_j, expected_u, 3e-14);
        vector_close(
            load.forces_world[0],
            Vec3::new(product * distance / radii[0].powi(3), 0.0, 0.0),
            3e-14,
        );
        assert_eq!(load.forces_world[1], load.forces_world[0].scale(-1.0));
        assert_eq!(load.pair_count, 1);
    }
}

#[test]
fn partial_overlap_matches_independent_volume_integration_and_force_derivatives() {
    for (radii, distance) in [
        ([1.0, 1.0], 0.2),
        ([1.0, 1.0], 1.0),
        ([1.0, 1.0], 1.9),
        ([1.6, 0.4], 1.3),
        ([1.6, 0.4], 1.95),
        ([1.0, 0.01], 0.995),
        ([1.0, 0.01], 1.005),
    ] {
        let properties = [
            sphere_mass(MASSES[0], radii[0]),
            sphere_mass(MASSES[1], radii[1]),
        ];
        let model = MutualGravity::new(FRAME, G, 0.0, 1)
            .unwrap()
            .uniform_spheres(&radii)
            .unwrap();
        let load = model
            .load(FRAME, &bodies(distance), &properties, || false)
            .unwrap();
        close(
            load.potential_energy_j,
            volume_potential_oracle(radii[0], radii[1], distance),
            3e-12,
        );
        let h = 5e-6;
        let derivative = (volume_potential_oracle(radii[0], radii[1], distance + h)
            - volume_potential_oracle(radii[0], radii[1], distance - h))
            / (2.0 * h);
        close(load.forces_world[0].x, derivative, 4e-8);
        assert!(load.forces_world[0].x > 0.0);
        assert_eq!(load.forces_world[1], load.forces_world[0].scale(-1.0));
    }
}

#[test]
fn equal_spheres_remain_finite_and_restore_linearly_near_coincidence() {
    let radii = [1.0, 1.0];
    let properties = [sphere_mass(MASSES[0], 1.0), sphere_mass(MASSES[1], 1.0)];
    let model = MutualGravity::new(FRAME, G, 0.0, 1)
        .unwrap()
        .uniform_spheres(&radii)
        .unwrap();
    let product = G * MASSES[0] * MASSES[1];
    for x in [0.0_f64, 1e-12, 1e-6, 0.25, 1.0, 1.99, 2.0] {
        let load = model
            .load(FRAME, &bodies(x), &properties, || false)
            .unwrap();
        let potential = product * (-1.2 + 0.5 * x * x - 3.0 * x.powi(3) / 16.0 + x.powi(5) / 160.0);
        let force = product * (x - 9.0 * x * x / 16.0 + x.powi(4) / 32.0);
        close(load.potential_energy_j, potential, 3e-14);
        close(
            load.forces_world[0].x,
            force,
            3e-13 * force.abs().max(1e-20),
        );
    }
}

#[test]
fn source_order_translation_and_orientation_preserve_the_same_pair_law() {
    let radii = [1.4, 0.6];
    let properties = [
        sphere_mass(MASSES[0], radii[0]),
        sphere_mass(MASSES[1], radii[1]),
    ];
    let base = MutualGravity::new(FRAME, G, 0.0, 1).unwrap();
    let states = bodies(1.2);
    let original = base
        .uniform_spheres(&radii)
        .unwrap()
        .load(FRAME, &states, &properties, || false)
        .unwrap();
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 3.0), 0.7).unwrap();
    let shift = Vec3::new(7.0, -11.0, 5.0);
    let transformed = states.map(|body| {
        RigidBodyState::new(
            Pose::new(
                rotation
                    .rotate_body_to_world(body.pose().position_world())
                    .add(shift),
                rotation,
            )
            .unwrap(),
            Vec3::ZERO,
            Vec3::ZERO,
        )
        .unwrap()
    });
    let rotated = base
        .uniform_spheres(&radii)
        .unwrap()
        .load(FRAME, &transformed, &properties, || false)
        .unwrap();
    close(
        rotated.potential_energy_j,
        original.potential_energy_j,
        3e-14,
    );
    vector_close(
        rotated.forces_world[0],
        rotation.rotate_body_to_world(original.forces_world[0]),
        3e-14,
    );
    let reverse_radii = [radii[1], radii[0]];
    let reversed = base
        .uniform_spheres(&reverse_radii)
        .unwrap()
        .load(
            FRAME,
            &[states[1], states[0]],
            &[properties[1], properties[0]],
            || false,
        )
        .unwrap();
    close(
        reversed.potential_energy_j,
        original.potential_energy_j,
        3e-14,
    );
    vector_close(reversed.forces_world[0], original.forces_world[1], 3e-14);
}

#[test]
fn source_mass_inertia_frame_budget_cutoff_and_cancellation_remain_binding() {
    let radii = [1.0, 0.5];
    let properties = [
        sphere_mass(MASSES[0], radii[0]),
        sphere_mass(MASSES[1], radii[1]),
    ];
    let initial = bodies(0.2);
    let base = MutualGravity::new(FRAME, G, 0.0, 1).unwrap();
    for invalid in [[0.0, 0.5], [-1.0, 0.5], [1.0, f64::NAN]] {
        assert!(base.uniform_spheres(&invalid).is_err());
    }
    let finite = base.uniform_spheres(&radii).unwrap();
    assert!(matches!(
        finite.load(FRAME + 1, &initial, &properties, || false),
        Err(GravityError::FrameMismatch)
    ));
    assert!(
        finite
            .load(FRAME, &initial[..1], &properties, || false)
            .is_err()
    );
    let mut inconsistent = properties;
    inconsistent[0] = MassProperties::new(MASSES[0], Vec3::ZERO, Vec3::new(0.8, 0.9, 1.0)).unwrap();
    assert!(
        finite
            .load(FRAME, &initial, &inconsistent, || false)
            .is_err()
    );
    assert!(
        MutualGravity::new(FRAME, G, 0.0, 0)
            .unwrap()
            .uniform_spheres(&radii)
            .and_then(|model| model.load(FRAME, &initial, &properties, || false))
            .is_err()
    );
    assert!(matches!(
        MutualGravity::new(FRAME, G, 0.1, 1)
            .unwrap()
            .uniform_spheres(&radii)
            .unwrap()
            .load(FRAME, &bodies(0.0), &properties, || false),
        Err(GravityError::SeparationCutoff { .. })
    ));
    assert!(matches!(
        finite.load(FRAME, &initial, &properties, || true),
        Err(GravityError::Cancelled)
    ));
    let mut polls = 0;
    finite
        .load(FRAME, &initial, &properties, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        finite.load(FRAME, &initial, &properties, || {
            seen += 1;
            seen == polls
        }),
        Err(GravityError::Cancelled)
    ));
}
