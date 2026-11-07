//! G1/G3/G4: cached hydrostatics retain analytic moments and bounded work.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3};
use fs_material::fluid_transport::FluidTransportState;
use fs_mbd::Gravity;
use fs_mbd_material::hydrostatic::{
    DisplacementMesh, HydrostaticError, HydrostaticMedium, HydrostaticReservoir,
};
use fs_mbd_material::principal::PrincipalFrame;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

// The same [-1,1]^3 cube, partitioned into actual nonoverlapping smaller cells.
fn cube(parts: usize) -> DisplacementMesh {
    let size = 2.0 / parts as f64;
    let axes = [
        Vec3::new(size, 0.0, 0.0),
        Vec3::new(0.0, size, 0.0),
        Vec3::new(0.0, 0.0, size),
    ];
    let mut cells = Vec::new();
    for z in 0..parts {
        for y in 0..parts {
            for x in 0..parts {
                let start = Vec3::new(
                    -1.0 + x as f64 * size,
                    -1.0 + y as f64 * size,
                    -1.0 + z as f64 * size,
                );
                for p in [
                    [0, 1, 2],
                    [0, 2, 1],
                    [1, 0, 2],
                    [1, 2, 0],
                    [2, 0, 1],
                    [2, 1, 0],
                ] {
                    cells.push([
                        start,
                        start + axes[p[0]],
                        start + axes[p[0]] + axes[p[1]],
                        start + Vec3::new(size, size, size),
                    ]);
                }
            }
        }
    }
    let count = cells.len();
    DisplacementMesh::new(cells, count, || false).unwrap()
}

fn medium<'a>(
    fluid: &'a FluidTransportState,
    datum: Vec3,
) -> HydrostaticMedium<'a, FluidTransportState> {
    HydrostaticMedium::new(
        fluid,
        Gravity::new(fs_mbd::Vec3::new(-3.0, -4.0, -12.0)).unwrap(),
        datum,
    )
    .unwrap()
}

#[test]
fn dry_and_submerged_query_work_is_independent_of_mesh_refinement_and_attitude() {
    let fluid = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
    let coarse = cube(1);
    let refined = cube(8);
    let up = Vec3::new(3.0 / 13.0, 4.0 / 13.0, 12.0 / 13.0);
    for rotation in [
        So3::identity(),
        So3::exp(So3Tangent::new(Vec3::new(0.4, -0.8, 0.2))).unwrap(),
    ] {
        for height in [-5.0, 5.0] {
            let pose = Se3::from_parts(rotation, up.scale(height)).unwrap();
            let mut counts = Vec::new();
            for mesh in [&coarse, &refined] {
                let mut calls = 0;
                let displaced = medium(&fluid, Vec3::new(0.0, 0.0, 0.0))
                    .displacement(mesh, pose, || {
                        calls += 1;
                        false
                    })
                    .unwrap();
                counts.push(calls);
                if height > 0.0 {
                    assert_eq!(displaced.volume_m3(), 0.0);
                    assert_eq!(displaced.centroid_body_m(), None);
                } else {
                    close(displaced.volume_m3(), 8.0);
                    let center = displaced.centroid_body_m().unwrap();
                    close(center.x, 0.0);
                    close(center.y, 0.0);
                    close(center.z, 0.0);
                }
            }
            // This checks eliminated geometry work without a timing threshold.
            assert_eq!(counts[0], counts[1]);
            assert!(
                counts[1] <= 3,
                "unnecessary whole-mesh traversal: {counts:?}"
            );
        }
    }
}

#[test]
fn large_world_datum_does_not_change_body_relative_cut_classification() {
    let fluid = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
    let mesh = cube(2);
    let gravity = Gravity::new(fs_mbd::Vec3::new(0.0, 0.0, -10.0)).unwrap();
    for height in [-3.0, 0.0, 3.0] {
        let mut results = Vec::new();
        for offset in [0.0, 1e12] {
            let medium =
                HydrostaticMedium::new(&fluid, gravity, Vec3::new(offset, -offset, offset))
                    .unwrap();
            let pose =
                Se3::from_parts(So3::identity(), Vec3::new(offset, -offset, offset + height))
                    .unwrap();
            results.push(medium.displacement(&mesh, pose, || false).unwrap());
        }
        assert_eq!(results[0], results[1]);
    }
}

#[test]
fn waterline_contacts_and_positive_thin_caps_keep_their_analytic_volume() {
    let fluid = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
    let medium = HydrostaticMedium::new(
        &fluid,
        Gravity::new(fs_mbd::Vec3::new(0.0, 0.0, -10.0)).unwrap(),
        Vec3::new(0.0, 0.0, 0.0),
    )
    .unwrap();
    let mesh = cube(1);
    for height in [-1.0, 0.0, 1.0] {
        let pose = Se3::from_parts(So3::identity(), Vec3::new(0.0, 0.0, height)).unwrap();
        let displaced = medium.displacement(&mesh, pose, || false).unwrap();
        close(displaced.volume_m3(), 4.0 * (1.0 - height));
    }
    let cap = 2.0_f64.powi(-30);
    let pose = Se3::from_parts(So3::identity(), Vec3::new(0.0, 0.0, 1.0 - cap)).unwrap();
    let displaced = medium.displacement(&mesh, pose, || false).unwrap();
    assert!(displaced.volume_m3() > 0.0);
    assert!((displaced.volume_m3() / (4.0 * cap) - 1.0).abs() < 1e-7);
    close(displaced.centroid_body_m().unwrap().z, -1.0 + 0.5 * cap);
}

#[test]
fn cached_queries_still_observe_initial_and_final_cancellation() {
    let fluid = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
    let medium = medium(&fluid, Vec3::new(0.0, 0.0, 0.0));
    let mesh = cube(8);
    for height in [-5.0, 5.0] {
        let pose = Se3::from_parts(So3::identity(), Vec3::new(0.0, 0.0, height)).unwrap();
        let mut count = 0;
        let accepted = medium
            .displacement(&mesh, pose, || {
                count += 1;
                false
            })
            .unwrap();
        for stop in 1..=count {
            let mut calls = 0;
            assert_eq!(
                medium.displacement(&mesh, pose, || {
                    calls += 1;
                    calls == stop
                }),
                Err(HydrostaticError::Cancelled)
            );
        }
        assert_eq!(
            medium.displacement(&mesh, pose, || false).unwrap(),
            accepted
        );
    }
}

fn check_large_body_origin_depth(z_offsets: [f64; 4], expected_mean_depth: f64) {
    let b = 2.0_f64.powi(52);
    let cells = vec![[
        Vec3::new(0.0, 0.0, b + z_offsets[0]),
        Vec3::new(1.0, 0.0, b + z_offsets[1]),
        Vec3::new(0.0, 1.0, b + z_offsets[2]),
        Vec3::new(0.0, 0.0, b + z_offsets[3]),
    ]];
    let mesh = DisplacementMesh::new(cells, 1, || false).unwrap();
    let fluid = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
    let medium = HydrostaticMedium::new(
        &fluid,
        Gravity::new(fs_mbd::Vec3::new(0.0, 0.0, -10.0)).unwrap(),
        Vec3::new(0.0, 0.0, b + 1.0),
    )
    .unwrap();
    // The inertial frame need not equal the displacement envelope. An origin-
    // centred frame retains the authored identity pose for this geometry test.
    let inertia = fs_mbd::articulated::SpatialInertia::new(
        1.0,
        Vec3::new(0.0, 0.0, 0.0),
        Mat3 {
            m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        },
    )
    .unwrap();
    let frame = PrincipalFrame::from_inertia(inertia, 1e-12, || false).unwrap();
    let state = frame
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let load = HydrostaticReservoir::load_frame(&medium, &mesh, &frame, state, || false).unwrap();
    close(load.displaced().volume_m3(), 1.0 / 6.0);
    close(load.wrench().force_world.z, 10000.0 / 6.0);
    // The exact mean of an affine depth over a tetrahedron is the arithmetic
    // mean of its four vertex depths; no absolute centroid subtraction is needed.
    close(
        load.pressure_potential_j(),
        (10000.0 / 6.0) * expected_mean_depth,
    );
}

#[test]
fn large_body_origin_does_not_erase_positive_full_immersion_depth() {
    // Absolute centroid b+3/4 rounds onto the surface b+1, despite mean depth 1/4.
    check_large_body_origin_depth([0.0, 1.0, 1.0, 1.0], 0.25);
}

#[test]
fn large_body_origin_does_not_distort_a_positive_full_immersion_depth() {
    // Absolute centroid b+1/4 rounds to b, giving depth 1 instead of the true 3/4.
    check_large_body_origin_depth([0.0, 0.0, 0.0, 1.0], 0.75);
}

#[test]
fn oblique_large_origin_retains_the_original_vertex_depth_owner() {
    let b = 2.0_f64.powi(52);
    let vertices = [
        Vec3::new(0.0, 0.0, b),
        Vec3::new(1.0, 0.0, b),
        Vec3::new(0.0, 1.0, b),
        Vec3::new(0.0, 0.0, b + 1.0),
    ];
    let down = Vec3::new(0.6, 0.0, 0.8);
    let origin_depth = down.dot(Vec3::new(0.0, 0.0, -b));
    let depths = vertices.map(|point| origin_depth + down.dot(point));
    assert_eq!(depths, [0.0, 0.5, 0.0, 1.0]);
    let mean_depth: f64 = depths.into_iter().map(|d| 0.25 * d).sum();
    assert_eq!(mean_depth, 0.375);
    // The anchored centroid gives 0.35 instead. Although it is closer to the
    // ideal affine geometry, changing only the fully wet branch would leave
    // the retained clipping/gradient owner with a different pressure potential.
    let mesh = DisplacementMesh::new(vec![vertices], 1, || false).unwrap();
    let fluid = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
    let medium = HydrostaticMedium::new(
        &fluid,
        Gravity::new(fs_mbd::Vec3::new(6.0, 0.0, 8.0)).unwrap(),
        Vec3::new(0.0, 0.0, b),
    )
    .unwrap();
    let inertia = fs_mbd::articulated::SpatialInertia::new(
        1.0,
        Vec3::new(0.0, 0.0, 0.0),
        Mat3 {
            m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        },
    )
    .unwrap();
    let frame = PrincipalFrame::from_inertia(inertia, 1e-12, || false).unwrap();
    let state = frame
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let load = HydrostaticReservoir::load_frame(&medium, &mesh, &frame, state, || false).unwrap();
    close(load.displaced().volume_m3(), 1.0 / 6.0);
    close(load.pressure_potential_j(), (10000.0 / 6.0) * mean_depth);
}

#[test]
fn full_cache_preserves_the_native_underflowed_volume_weight_refusal() {
    let tetrahedron = |origin: Vec3, edge: f64| {
        [
            origin,
            origin + Vec3::new(edge, 0.0, 0.0),
            origin + Vec3::new(0.0, edge, 0.0),
            origin + Vec3::new(0.0, 0.0, edge),
        ]
    };
    // Both disjoint cells have finite positive volumes and inertias, but the
    // second one's volume divided by their running total underflows to zero.
    // Mesh admission remains valid; the original query refusal must survive.
    let mesh = DisplacementMesh::new(
        vec![
            tetrahedron(Vec3::new(0.0, 0.0, 0.0), 1e60),
            tetrahedron(Vec3::new(-2e-60, -2e-60, -2e-60), 1e-60),
        ],
        2,
        || false,
    )
    .unwrap();
    assert!(mesh.volume_m3().is_finite() && mesh.volume_m3() > 0.0);
    let fluid = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
    let medium = HydrostaticMedium::new(
        &fluid,
        Gravity::new(fs_mbd::Vec3::new(0.0, 0.0, -1.0)).unwrap(),
        Vec3::new(0.0, 0.0, 0.0),
    )
    .unwrap();
    let pose = Se3::from_parts(So3::identity(), Vec3::new(0.0, 0.0, -3e60)).unwrap();
    assert_eq!(
        medium.displacement(&mesh, pose, || false),
        Err(HydrostaticError::Invalid(
            "clipped volume moments are unrepresentable",
        )),
    );
}
