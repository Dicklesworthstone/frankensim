use super::*;
use fs_ga::{Mat3, Twist};

const POSITIONS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
// f(x,y,z) = (y, 2z, 3x), with genuinely varying force directions.
const DENSITY: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 3.0],
    [1.0, 0.0, 0.0],
    [0.0, 2.0, 0.0],
];
const LIMITS: SpatialLorentzLoadLimits = SpatialLorentzLoadLimits {
    max_vertices: 5,
    max_cells: 2,
};

fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 2e-13, "{actual} != {expected}");
}
fn fixture(positions: &[[f64; 3]]) -> SpatialLorentzLoads {
    integrate(
        positions,
        [CellLoad {
            vertices: [0, 1, 2, 3],
            volume_m3: 1.0 / 6.0,
            force_density_n_m3: &DENSITY,
        }]
        .into_iter(),
        LIMITS,
        &mut || false,
    )
    .unwrap()
}
fn velocity(p: [f64; 3]) -> [f64; 3] {
    [
        1.0 + 2.0 * p[0] - p[1],
        -3.0 + p[1] + p[2],
        2.0 - 2.0 * p[0] + p[2],
    ]
}

#[test]
fn nonuniform_force_preserves_analytic_torque_and_p1_virtual_work() {
    let loads = fixture(&POSITIONS);
    // Independently integrate monomials on the unit tetrahedron:
    // integral x_i = 1/24, x_i^2 = 1/60, x_i*x_j = 1/120.
    let wrench = loads.wrench(&Se3::identity(), [0.0; 3]).unwrap();
    for (actual, expected) in wrench.to_array().into_iter().zip([
        -1.0 / 120.0,
        -1.0 / 24.0,
        0.0,
        1.0 / 24.0,
        1.0 / 12.0,
        1.0 / 8.0,
    ]) {
        near(actual, expected);
    }
    // v.f integrates to 1/24 + 7/120 - 5/60 = 1/60.
    near(
        loads
            .power_w(&POSITIONS.map(velocity), &mut || false)
            .unwrap(),
        1.0 / 60.0,
    );
    // Density sampling with lumped V/4 weights would give torque -1/12 in x.
    assert!((wrench.torque.x + 1.0 / 12.0).abs() > 0.05);
}

#[test]
fn arbitrary_frame_and_origin_preserve_rigid_virtual_work() {
    let loads = fixture(&POSITIONS);
    let placement = Se3::exp(Twist::new(
        Vec3::new(0.2, -0.4, 0.3),
        Vec3::new(2.0, -1.0, 3.0),
    ))
    .unwrap();
    let origin = [4.0, -2.0, 1.0];
    let angular = Vec3::new(0.7, -0.2, 0.4);
    let linear = Vec3::new(-0.5, 0.3, 1.1);
    let receiver_twist = Twist::new(angular, linear);
    let inverse = placement.inverse().unwrap();
    let original_velocity = POSITIONS.map(|point| {
        let arm = placement.transform_point(vector(point)).unwrap() - vector(origin);
        let v = inverse
            .transform_vector(linear + angular.cross(arm))
            .unwrap();
        [v.x, v.y, v.z]
    });
    let wrench = loads.wrench(&placement, origin).unwrap();
    near(
        wrench.pairing(receiver_twist),
        loads.power_w(&original_velocity, &mut || false).unwrap(),
    );
    let other_origin = [-1.0, 3.0, 2.0];
    let shifted = loads.wrench(&placement, other_origin).unwrap();
    let expected = wrench.torque + (vector(origin) - vector(other_origin)).cross(wrench.force);
    for (a, b) in [shifted.torque.x, shifted.torque.y, shifted.torque.z]
        .into_iter()
        .zip([expected.x, expected.y, expected.z])
    {
        near(a, b);
    }
}

#[test]
fn original_nodal_loads_drive_native_linear_and_angular_momentum() {
    use fs_mbd::articulated::{
        ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
        apply_free_floating_impulses,
    };
    let load = fixture(&POSITIONS)
        .wrench(&Se3::identity(), [0.0; 3])
        .unwrap();
    let dt = 0.07;
    let impulse = Wrench::new(load.torque.scale(dt), load.force.scale(dt));
    let inertia = SpatialInertia::new(
        2.0,
        Vec3::new(0.0, 0.0, 0.0),
        Mat3 {
            m: [3.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 5.0],
        },
    )
    .unwrap();
    let model = ArticulatedModel::new(vec![Link::new(
        "conductor",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia,
    )])
    .unwrap();
    let response = apply_free_floating_impulses(
        &model,
        FreeFloatingBaseState::stationary(Se3::identity()),
        &[],
        &[],
        &[],
        &[impulse],
    )
    .unwrap();
    let actual = Wrench::new(
        response.angular_momentum_change_world,
        response.linear_momentum_change_world,
    );
    for (a, b) in actual.to_array().into_iter().zip(impulse.to_array()) {
        near(a, b);
    }
    // Nonzero spin and translation reach the actual mechanical state.
    near(
        response.base_after.twist_body.angular.x,
        -dt / (120.0 * 3.0),
    );
    near(response.base_after.twist_body.angular.y, -dt / (24.0 * 4.0));
    near(response.base_after.twist_body.linear.z, dt / (8.0 * 2.0));
    assert!(response.kinetic_energy_change > 0.0);
    near(response.kinetic_work_residual, 0.0);
}

#[test]
fn refinement_and_shared_vertex_scatter_preserve_force_torque_and_power() {
    let positions = [
        POSITIONS[0],
        POSITIONS[1],
        POSITIONS[2],
        POSITIONS[3],
        [0.5, 0.0, 0.0],
    ];
    let left = [DENSITY[0], [0.0, 0.0, 1.5], DENSITY[2], DENSITY[3]];
    let right = [[0.0, 0.0, 1.5], DENSITY[1], DENSITY[2], DENSITY[3]];
    let loads = integrate(
        &positions,
        [
            CellLoad {
                vertices: [0, 4, 2, 3],
                volume_m3: 1.0 / 12.0,
                force_density_n_m3: &left,
            },
            CellLoad {
                vertices: [4, 1, 2, 3],
                volume_m3: 1.0 / 12.0,
                force_density_n_m3: &right,
            },
        ]
        .into_iter(),
        LIMITS,
        &mut || false,
    )
    .unwrap();
    let original = fixture(&POSITIONS);
    for (a, b) in loads
        .wrench(&Se3::identity(), [0.0; 3])
        .unwrap()
        .to_array()
        .into_iter()
        .zip(
            original
                .wrench(&Se3::identity(), [0.0; 3])
                .unwrap()
                .to_array(),
        )
    {
        near(a, b);
    }
    near(
        loads
            .power_w(&positions.map(velocity), &mut || false)
            .unwrap(),
        1.0 / 60.0,
    );
}

#[test]
fn nearby_moment_reference_retains_small_torque_under_large_translation() {
    let offset = [1e12, -2e12, 3e12];
    let translated = POSITIONS.map(|p| [p[0] + offset[0], p[1] + offset[1], p[2] + offset[2]]);
    for (a, b) in fixture(&translated)
        .wrench(&Se3::identity(), offset)
        .unwrap()
        .to_array()
        .into_iter()
        .zip(
            fixture(&POSITIONS)
                .wrench(&Se3::identity(), [0.0; 3])
                .unwrap()
                .to_array(),
        )
    {
        near(a, b);
    }
}

#[test]
fn extraction_limits_cancellation_and_invalid_values_refuse_atomically() {
    let cell = || {
        [CellLoad {
            vertices: [0, 1, 2, 3],
            volume_m3: 1.0 / 6.0,
            force_density_n_m3: &DENSITY,
        }]
        .into_iter()
    };
    assert!(matches!(
        integrate(&POSITIONS, cell(), LIMITS, &mut || true),
        Err(SpatialLorentzLoadError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        integrate(&POSITIONS, cell(), LIMITS, &mut || {
            polls += 1;
            polls == 4
        }),
        Err(SpatialLorentzLoadError::Cancelled)
    ));
    assert!(
        integrate(
            &POSITIONS,
            cell(),
            SpatialLorentzLoadLimits {
                max_vertices: 3,
                ..LIMITS
            },
            &mut || false
        )
        .is_err()
    );
    let mut bad = DENSITY;
    bad[1][2] = f64::NAN;
    assert!(
        integrate(
            &POSITIONS,
            [CellLoad {
                vertices: [0, 1, 2, 3],
                volume_m3: 1.0 / 6.0,
                force_density_n_m3: &bad
            }]
            .into_iter(),
            LIMITS,
            &mut || false
        )
        .is_err()
    );
    let loads = fixture(&POSITIONS);
    assert!(
        loads
            .wrench(&Se3::identity(), [f64::NAN, 0.0, 0.0])
            .is_err()
    );
    assert!(loads.power_w(&[[0.0; 3]; 3], &mut || false).is_err());
    assert!(matches!(
        loads.power_w(&[[0.0; 3]; 4], &mut || true),
        Err(SpatialLorentzLoadError::Cancelled)
    ));
}
