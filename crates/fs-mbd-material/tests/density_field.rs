//! Sourced continuous density -> exact moments -> actual native mechanics.
mod common;

#[path = "density_field/motion.rs"]
mod motion;

#[path = "density_field/partition.rs"]
mod partition;

use common::{close, density_card, point};
use fs_ga::{Mat3, Se3, So3, Twist, Vec3};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd::Gravity;
use fs_mbd::articulated::SpatialInertia;
use fs_mbd_material::density_field::{DensityVertexBinding, TetrahedralDensityBody};
use fs_mbd_material::hydrostatic::{
    DisplacementMesh, HydrostaticMedium, HydrostaticReservoir, HydrostaticStepControls,
};
use fs_mbd_material::{MaterialInertiaError, RegionBinding, VolumeMoments, resolve_body};

fn vertices() -> [Vec3; 4] {
    [
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ]
}
fn resolve(
    positions: &[Vec3],
    densities: &[f64],
    tetrahedra: &[[usize; 4]],
) -> TetrahedralDensityBody {
    let cards: Vec<_> = densities
        .iter()
        .enumerate()
        .map(|(i, &rho)| density_card(&format!("synthetic-vertex-{i}"), rho, rho))
        .collect();
    let point = point(300.0);
    let bindings: Vec<_> = positions
        .iter()
        .zip(&cards)
        .map(|(&position_m, material)| DensityVertexBinding {
            position_m,
            material,
            point: &point,
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
        .collect();
    TetrahedralDensityBody::resolve(&bindings, tetrahedra, positions.len(), || false).unwrap()
}
fn vec_close(a: Vec3, b: Vec3) {
    close(a.x, b.x);
    close(a.y, b.y);
    close(a.z, b.z);
}
fn tensor_close(a: Mat3, b: Mat3) {
    for (a, b) in a.m.into_iter().zip(b.m) {
        close(a, b);
    }
}
fn inertia_close(a: SpatialInertia, b: SpatialInertia) {
    close(a.mass(), b.mass());
    vec_close(a.center_of_mass(), b.center_of_mass());
    tensor_close(a.inertia_com(), b.inertia_com());
}
fn ga(v: fs_mbd::Vec3) -> Vec3 {
    Vec3::new(v.x, v.y, v.z)
}

// Independent Cartesian integrals over x,y,z >= 0, x+y+z <= 1:
// integral x^a y^b z^c dV = a! b! c! / (a+b+c+3)!.
fn analytic_gradient(rho: f64, slope: f64) -> (f64, Vec3, Mat3) {
    let mass = rho / 6.0 + slope / 24.0;
    let center = Vec3::new(
        (rho / 24.0 + slope / 60.0) / mass,
        (rho / 24.0 + slope / 120.0) / mass,
        (rho / 24.0 + slope / 120.0) / mass,
    );
    let xx = rho / 60.0 + slope / 120.0;
    let yy = rho / 60.0 + slope / 360.0;
    let xy = rho / 120.0 + slope / 360.0;
    let yz = rho / 120.0 + slope / 720.0;
    let c = [center.x, center.y, center.z];
    let mut covariance = [xx, xy, xy, xy, yy, yz, xy, yz, yy];
    for row in 0..3 {
        for col in 0..3 {
            covariance[3 * row + col] -= mass * c[row] * c[col];
        }
    }
    (mass, center, Mat3 { m: covariance })
}
fn inertia_from_covariance(covariance: Mat3) -> Mat3 {
    let trace = covariance.m[0] + covariance.m[4] + covariance.m[8];
    let mut inertia = covariance.m.map(|entry| -entry);
    for i in [0, 4, 8] {
        inertia[i] += trace;
    }
    Mat3 { m: inertia }
}

#[test]
fn graded_density_matches_independent_cubic_cartesian_integrals() {
    for (rho, slope) in [(1000.0, 4000.0), (5000.0, -4000.0)] {
        let body = resolve(&vertices(), &[rho, rho + slope, rho, rho], &[[0, 1, 2, 3]]);
        let (mass, center, covariance) = analytic_gradient(rho, slope);
        close(body.inertia().mass(), mass);
        vec_close(body.inertia().center_of_mass(), center);
        tensor_close(
            body.inertia().inertia_com(),
            inertia_from_covariance(covariance),
        );
        close(body.cells()[0].geometry().volume_m3(), 1.0 / 6.0);
        for (index, vertex) in body.vertices().iter().enumerate() {
            close(
                vertex.density_kg_m3(),
                if index == 1 { rho + slope } else { rho },
            );
            assert_eq!(vertex.material().properties().len(), 1);
            assert_eq!(vertex.material().state_coordinate("T"), Some(300.0));
        }
        // A uniform inertia using the correct average density still has the
        // wrong COM and cross moments for this actual material distribution.
        let uniform = resolve(&vertices(), &[rho + slope / 4.0; 4], &[[0, 1, 2, 3]]);
        close(uniform.inertia().mass(), mass);
        assert!((uniform.inertia().center_of_mass().x - center.x).abs() > 0.02);
        assert!(
            (uniform.inertia().inertia_com().m[1] - body.inertia().inertia_com().m[1]).abs() > 0.1
        );
    }
}

#[test]
fn uniform_sources_recover_the_original_region_mass_owner() {
    let positions =
        vertices().map(|v| Vec3::new(4.0 + 2.0 * v.x + 0.3 * v.y, -2.0 + v.y, 3.0 + 0.7 * v.z));
    let graded = resolve(&positions, &[2400.0; 4], &[[0, 1, 2, 3]]);
    let card = density_card("synthetic-uniform", 2400.0, 2400.0);
    let point = point(300.0);
    let uniform = resolve_body(
        &[RegionBinding {
            key: "tetrahedron",
            geometry: VolumeMoments::tetrahedron(positions).unwrap(),
            material: &card,
            point: &point,
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap();
    inertia_close(graded.inertia(), uniform.inertia());
    assert_eq!(
        graded.deform(&positions, || false).unwrap().inertia(),
        graded.inertia()
    );
}

#[test]
fn exact_density_moments_survive_partition_and_vertex_order_changes() {
    let original = resolve(
        &vertices(),
        &[1000.0, 5000.0, 1000.0, 1000.0],
        &[[0, 1, 2, 3]],
    );
    let mut refined_vertices = vertices().to_vec();
    refined_vertices.push(Vec3::new(0.25, 0.25, 0.25));
    let cells = [[4, 1, 2, 3], [0, 4, 2, 3], [0, 1, 4, 3], [0, 1, 2, 4]];
    let refined = resolve(
        &refined_vertices,
        &[1000.0, 5000.0, 1000.0, 1000.0, 2000.0],
        &cells,
    );
    inertia_close(original.inertia(), refined.inertia());
    let reversed = cells
        .into_iter()
        .rev()
        .map(|[a, b, c, d]| [d, c, b, a])
        .collect::<Vec<_>>();
    let reordered = resolve(
        &refined_vertices,
        &[1000.0, 5000.0, 1000.0, 1000.0, 2000.0],
        &reversed,
    );
    assert_eq!(refined, reordered);
    let positions = vertices();
    let weights = [1000.0, 5000.0, 1000.0, 1000.0];
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    if a == b || a == c || a == d || b == c || b == d || c == d {
                        continue;
                    }
                    let order = [a, b, c, d];
                    let inertia = SpatialInertia::from_weighted_tetrahedron(
                        original.inertia().mass(),
                        order.map(|i| positions[i]),
                        order.map(|i| weights[i]),
                    )
                    .unwrap();
                    inertia_close(inertia, original.inertia());
                }
            }
        }
    }
}

#[test]
fn affine_expansion_and_shear_transport_the_same_graded_mass() {
    let body = resolve(
        &vertices(),
        &[1000.0, 5000.0, 1000.0, 1000.0],
        &[[0, 1, 2, 3]],
    );
    let transform = Mat3 {
        m: [2.0, 0.4, 0.0, 0.0, 3.0, -0.2, 0.0, 0.0, 0.5],
    };
    let translation = Vec3::new(10.0, -3.0, 4.0);
    let current = vertices().map(|v| transform.apply(v) + translation);
    let deformed = body.deform(&current, || false).unwrap();
    let (mass, center, covariance) = analytic_gradient(1000.0, 4000.0);
    assert_eq!(deformed.inertia().mass(), body.inertia().mass());
    close(deformed.inertia().mass(), mass);
    vec_close(
        deformed.inertia().center_of_mass(),
        transform.apply(center) + translation,
    );
    tensor_close(
        deformed.inertia().inertia_com(),
        inertia_from_covariance(transform.compose(covariance).compose(transform.transpose())),
    );
    close(deformed.cells()[0].geometry().volume_m3(), 0.5);
    for (current, reference) in deformed.cells()[0]
        .vertex_density_kg_m3()
        .into_iter()
        .zip(body.cells()[0].vertex_density_kg_m3())
    {
        close(current, reference / 3.0);
    }
    assert_eq!(deformed, body.deform(&current, || false).unwrap());
    assert_eq!(deformed.reference(), &body);
    assert_eq!(
        body.deform(&vertices(), || false).unwrap().inertia(),
        body.inertia()
    );
    deformed.principal_frame(1e-12, || false).unwrap();
}

#[test]
fn different_cell_jacobians_preserve_local_mass_without_averaging_shared_density() {
    let mut positions = vertices().to_vec();
    positions.push(Vec3::new(0.0, 0.0, -1.0));
    let body = resolve(
        &positions,
        &[1000.0, 5000.0, 1000.0, 1000.0, 1000.0],
        &[[0, 1, 2, 3], [0, 1, 2, 4]],
    );
    positions[3].z = 2.0;
    positions[4].z = -3.0;
    let deformed = body.deform(&positions, || false).unwrap();
    assert_eq!(deformed.inertia().mass(), body.inertia().mass());
    for ((current, reference), jacobian) in
        deformed.cells().iter().zip(body.cells()).zip([2.0, 3.0])
    {
        assert_eq!(current.inertia().mass(), reference.inertia().mass());
        for (rho, original) in current
            .vertex_density_kg_m3()
            .into_iter()
            .zip(reference.vertex_density_kg_m3())
        {
            close(rho, original / jacobian);
        }
    }
    assert_ne!(
        deformed.cells()[0].vertex_density_kg_m3()[0],
        deformed.cells()[1].vertex_density_kg_m3()[0]
    );
}

#[test]
fn density_gradient_changes_native_impulse_response_and_buoyant_torque() {
    let mesh = DisplacementMesh::new(vec![vertices()], 1, || false).unwrap();
    let fluid = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
    let gravity = Gravity::new(fs_mbd::Vec3::new(0.0, 0.0, -10.0)).unwrap();
    let medium = HydrostaticMedium::new(&fluid, gravity, Vec3::new(0.0, 0.0, 0.0)).unwrap();
    for scale in [1.0, 2.0] {
        let body = resolve(
            &vertices(),
            &[1000.0, 5000.0, 1000.0, 1000.0].map(|rho| rho * scale),
            &[[0, 1, 2, 3]],
        );
        let frame = body.principal_frame(1e-12, || false).unwrap();
        let pose = Se3::from_parts(So3::identity(), Vec3::new(0.0, 0.0, -2.0)).unwrap();
        let state = frame.state_from_kinematics(pose, Twist::zero()).unwrap();
        let arm = frame.arm_for_body_point(Vec3::new(0.0, 0.0, 0.0)).unwrap();
        let impulse = state
            .apply_impulse_at_body_point(frame.properties(), arm, fs_mbd::Vec3::new(0.0, 0.0, 1.0))
            .unwrap();
        let (mass, center, covariance) = analytic_gradient(1000.0 * scale, 4000.0 * scale);
        close(
            impulse
                .state_after
                .center_of_mass_velocity_world(frame.properties())
                .unwrap()
                .z,
            1.0 / mass,
        );
        let omega_principal = frame
            .properties()
            .angular_velocity_body_checked(impulse.state_after.angular_momentum_body())
            .unwrap();
        let omega = frame
            .body_from_principal()
            .rotation()
            .rotate(ga(omega_principal))
            .unwrap();
        vec_close(
            inertia_from_covariance(covariance).apply(omega),
            center.scale(-1.0).cross(Vec3::new(0.0, 0.0, 1.0)),
        );
        assert!(omega.norm() > 0.01);
        let load =
            HydrostaticReservoir::load_frame(&medium, &mesh, &frame, state, || false).unwrap();
        close(load.wrench().force_world.z, 10000.0 / 6.0);
        let torque = frame
            .body_from_principal()
            .rotation()
            .rotate(ga(load.wrench().torque_body))
            .unwrap();
        vec_close(
            torque,
            (Vec3::new(0.25, 0.25, 0.25) - center).cross(Vec3::new(0.0, 0.0, 10000.0 / 6.0)),
        );
        let step = HydrostaticReservoir::step_frame(
            &medium,
            &mesh,
            &frame,
            state,
            HydrostaticStepControls {
                duration_s: 1e-4,
                maximum_energy_defect_j: 1e-3,
            },
            || false,
        )
        .unwrap();
        close(
            step.dynamics
                .state_after
                .center_of_mass_velocity_world(frame.properties())
                .unwrap()
                .z
                / 1e-4,
            -10.0 + (10000.0 / 6.0) / mass,
        );
    }
}

#[test]
fn invalid_geometry_or_source_and_cancelled_deformation_publish_no_body() {
    let original = vertices();
    let body = resolve(
        &original,
        &[1000.0, 5000.0, 1000.0, 1000.0],
        &[[0, 1, 2, 3]],
    );
    let mut collapsed = original;
    collapsed[3] = Vec3::new(0.2, 0.2, 0.0);
    assert!(SpatialInertia::from_weighted_tetrahedron(1.0, collapsed, [1.0; 4]).is_err());
    assert!(
        SpatialInertia::from_weighted_tetrahedron(1.0, original, [1.0, 0.0, 1.0, 1.0]).is_err()
    );
    assert!(
        SpatialInertia::from_weighted_tetrahedron(1.0, original, [1.0, f64::NAN, 1.0, 1.0])
            .is_err()
    );
    let mut inverted = original;
    inverted[3].z = -1.0;
    assert!(body.deform(&inverted, || false).is_err());
    let expanded = original.map(|v| v.scale(1.5));
    let mut polls = 0;
    body.deform(&expanded, || {
        polls += 1;
        false
    })
    .unwrap();
    for checkpoint in 1..=polls {
        let mut call = 0;
        assert_eq!(
            body.deform(&expanded, || {
                call += 1;
                call == checkpoint
            }),
            Err(MaterialInertiaError::Cancelled)
        );
    }
    let card = density_card("unsupported-node-state", 1000.0, 2000.0);
    let unsupported = point(500.0);
    let bindings = original.map(|position_m| DensityVertexBinding {
        position_m,
        material: &card,
        point: &unsupported,
        selection: MaterialPropertySelection::SingleClaimOnly,
    });
    assert!(matches!(
        TetrahedralDensityBody::resolve(&bindings, &[[0, 1, 2, 3]], 4, || false),
        Err(MaterialInertiaError::Material(_))
    ));
    assert_eq!(
        body.deform(&original, || false).unwrap().inertia(),
        body.inertia()
    );
}
