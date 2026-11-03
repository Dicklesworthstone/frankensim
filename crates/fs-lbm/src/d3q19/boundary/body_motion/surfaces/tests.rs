use super::*;
use fs_ga::Mat3;

#[test]
fn full_turn_with_clear_endpoints_still_refuses_the_sphere_sweep() {
    let body = MovingPartialSaturationBody3 {
        center: [0.0; 3],
        orientation: So3::identity(),
        mass: 1.0,
        inertia: InertiaTensor3::new(Mat3 {
            m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        })
        .unwrap(),
        velocity: [0.0; 3],
        angular_momentum_body: Vec3::new(0.0, 0.0, 4.0 * core::f64::consts::PI),
    };
    let start =
        MovingFiniteLubricationPlane3::new(body, [0.0; 3], [1.0, 0.0, 0.0], 0.1, 0.001).unwrap();
    let end = start.at_body(drift(body, 0.5).unwrap()).unwrap();
    let center = [[1.2, 0.0, 0.0]];
    for plane in [start, end] {
        assert!(dot(plane.normal(), center[0]) - 1.0 > 0.19);
    }
    assert!(
        admit_spheres(
            &[MovingLubricationSurface3::Finite(start)],
            &[MovingLubricationSurface3::Finite(end)],
            center,
            center,
            [1.0],
            &CancelGate::new(),
        )
        .is_err()
    );
}

#[test]
fn separated_parallel_endpoints_cannot_hide_a_duplicate_mid_sweep() {
    let plane = |x, v| {
        MovingLubricationSurface3::Prescribed(
            LubricationPlane3::new([x, 0.0, 0.0], [1.0, 0.0, 0.0], [v, 0.0, 0.0], 0.1, 0.001)
                .unwrap(),
        )
    };
    let surfaces = [plane(0.0, 0.04), plane(0.01, -0.04)];
    // They coincide inside the first half-drift, although neither endpoint
    // duplicates a plane. Repeated film resistance must not be admitted.
    assert!(PreparedSurfaces::new(&surfaces, &[[0.0; 6]; 2], &CancelGate::new()).is_err());
    let separated = [plane(0.0, 0.001), plane(0.01, 0.001)];
    assert!(PreparedSurfaces::new(&separated, &[[0.0; 6]; 2], &CancelGate::new()).is_ok());
}

#[test]
fn small_rotation_bounds_the_full_translated_endpoint_arm() {
    let body = MovingPartialSaturationBody3 {
        center: [0.0; 3],
        orientation: So3::identity(),
        mass: 1.0,
        inertia: InertiaTensor3::new(Mat3 {
            m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        })
        .unwrap(),
        velocity: [0.0, -0.08, 0.0],
        angular_momentum_body: Vec3::new(0.0, 0.0, 0.02),
    };
    let moving =
        MovingFiniteLubricationPlane3::new(body, [0.0; 3], [1.0, 0.0, 0.0], 0.1, 0.001).unwrap();
    let quarter = moving.at_body(drift(body, 0.25).unwrap()).unwrap();
    let halfway = moving.at_body(drift(body, 0.5).unwrap()).unwrap();
    let prescribed =
        LubricationPlane3::new(quarter.point(), quarter.normal(), [0.0; 3], 0.1, 0.001).unwrap();
    // The two planes coincide at t=1/4, but their endpoint normals differ.
    // Relative to the moving carrier's initial center, its initial arm is
    // zero and its final arm is 0.04. The rotation bound must include
    // travel*0.04, never travel^2*0.04 when travel is below one.
    assert_ne!(moving.normal(), prescribed.normal());
    assert_ne!(halfway.normal(), prescribed.normal());
    assert!(
        MovingLubricationSurface3::Finite(moving)
            .normal_travel()
            .unwrap()
            < 1.0
    );
    let surfaces = [
        MovingLubricationSurface3::Prescribed(prescribed),
        MovingLubricationSurface3::Finite(moving),
    ];
    assert!(PreparedSurfaces::new(&surfaces, &[[0.0; 6]; 2], &CancelGate::new()).is_err());
}
