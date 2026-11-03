//! G1/G3/G4: retain original graded matter through an actual topology split.

use super::*;

fn refined() -> TetrahedralDensityBody {
    let mut points = vertices().to_vec();
    points.push(Vec3::new(0.25, 0.25, 0.25));
    resolve(
        &points,
        &[1000.0, 5000.0, 1000.0, 1000.0, 2000.0],
        &[[4, 1, 2, 3], [0, 4, 2, 3], [0, 1, 4, 3], [0, 1, 2, 4]],
    )
}

fn check_source_maps(
    original: &TetrahedralDensityBody,
    part: &fs_mbd_material::density_field::TetrahedralDensityPart,
) {
    assert_eq!(
        part.body.vertices().len(),
        part.original_vertex_indices.len()
    );
    assert_eq!(part.body.cells().len(), part.original_cell_indices.len());
    assert!(part.original_vertex_indices.windows(2).all(|p| p[0] < p[1]));
    assert!(part.original_cell_indices.windows(2).all(|p| p[0] < p[1]));
    for (local, &source) in part.original_vertex_indices.iter().enumerate() {
        assert_eq!(part.body.vertices()[local], original.vertices()[source]);
    }
    for (local, &source) in part.original_cell_indices.iter().enumerate() {
        assert_eq!(part.body.cells()[local], original.cells()[source]);
        assert_eq!(
            part.body.tetrahedra()[local].map(|v| part.original_vertex_indices[v]),
            original.tetrahedra()[source],
        );
    }
}

#[test]
fn g1_partition_preserves_graded_mass_full_inertia_and_every_source_receipt() {
    let mut points = vertices().to_vec();
    points.push(Vec3::new(0.0, 0.0, -2.0));
    points.push(Vec3::new(8.0, 7.0, 6.0)); // Unused material vertex carries no cell mass.
    let body = resolve(
        &points,
        &[1000.0, 5000.0, 1000.0, 1000.0, 1000.0, 1777.0],
        &[[0, 1, 2, 3], [0, 1, 2, 4]],
    );
    let original = body.clone();
    let parts = body.partition_cells(&[1], 2, 8, || false).unwrap();
    assert_eq!(parts.selected.original_vertex_indices, [0, 1, 2, 4]);
    assert_eq!(parts.retained.original_vertex_indices, [0, 1, 2, 3]);
    assert_eq!(parts.selected.original_cell_indices, [1]);
    assert_eq!(parts.retained.original_cell_indices, [0]);
    check_source_maps(&body, &parts.selected);
    check_source_maps(&body, &parts.retained);

    // Independent Cartesian moments for rho=1000+4000*x on a unit tet.
    // The selected tet is its z-reflection/expansion with |det A|=2.
    let (mass, center, covariance) = analytic_gradient(1000.0, 4000.0);
    close(parts.retained.body.inertia().mass(), mass);
    vec_close(parts.retained.body.inertia().center_of_mass(), center);
    tensor_close(
        parts.retained.body.inertia().inertia_com(),
        inertia_from_covariance(covariance),
    );
    let a = Mat3 {
        m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -2.0],
    };
    let transformed = a.compose(covariance).compose(a.transpose());
    close(parts.selected.body.inertia().mass(), 2.0 * mass);
    vec_close(
        parts.selected.body.inertia().center_of_mass(),
        a.apply(center),
    );
    tensor_close(
        parts.selected.body.inertia().inertia_com(),
        inertia_from_covariance(Mat3 {
            m: transformed.m.map(|x| 2.0 * x),
        }),
    );
    assert_eq!(
        parts.selected.body.cells()[0].inertia().mass().to_bits(),
        body.cells()[1].inertia().mass().to_bits()
    );
    assert_eq!(body, original);
}

#[test]
fn g3_partition_order_and_original_native_twist_keep_total_momenta_and_kinetic_energy() {
    let body = refined();
    let split = body.partition_cells(&[3, 1], 4, 16, || false).unwrap();
    assert_eq!(
        split,
        body.partition_cells(&[1, 3], 4, 16, || false).unwrap()
    );
    let rotation = So3::exp(fs_ga::So3Tangent::new(Vec3::new(0.3, -0.5, 0.2))).unwrap();
    let pose = Se3::from_parts(rotation, Vec3::new(2.0, -3.0, 4.0)).unwrap();
    let twist = Twist::new(Vec3::new(0.8, -0.4, 1.2), Vec3::new(0.3, 0.2, -0.1));
    let frame = body.principal_frame(1e-12, || false).unwrap();
    let state = frame.state_from_kinematics(pose, twist).unwrap();
    let original_p = ga(state.linear_momentum_world());
    let original_l = ga(state.pose().position_world()).cross(original_p)
        + ga(state
            .pose()
            .orientation()
            .rotate_body_to_world(state.angular_momentum_body()));
    let mut p = Vec3::new(0.0, 0.0, 0.0);
    let mut l = p;
    let (mut kinetic, mut mass) = (0.0, 0.0);
    for part in [&split.selected, &split.retained] {
        check_source_maps(&body, part);
        let frame = part.body.principal_frame(1e-12, || false).unwrap();
        let state = frame.state_from_kinematics(pose, twist).unwrap();
        let momentum = ga(state.linear_momentum_world());
        p = p + momentum;
        l = l
            + ga(state.pose().position_world()).cross(momentum)
            + ga(state
                .pose()
                .orientation()
                .rotate_body_to_world(state.angular_momentum_body()));
        kinetic += state.kinetic_energy(frame.properties()).unwrap();
        mass += frame.properties().mass();
    }
    close(mass, body.inertia().mass());
    vec_close(p, original_p);
    vec_close(l, original_l);
    close(kinetic, state.kinetic_energy(frame.properties()).unwrap());
}

#[test]
fn g1_compact_parts_continue_the_original_density_and_orientation_deformation() {
    let body = refined();
    let split = body.partition_cells(&[0, 2], 4, 16, || false).unwrap();
    let a = Mat3 {
        m: [1.2, 0.2, 0.0, 0.0, 0.8, -0.1, 0.0, 0.0, 1.4],
    };
    let shift = Vec3::new(0.7, -0.2, 0.5);
    let coordinates: Vec<_> = body
        .vertices()
        .iter()
        .map(|v| a.apply(v.position_m()) + shift)
        .collect();
    let all = body.deform(&coordinates, || false).unwrap();
    for part in [&split.selected, &split.retained] {
        let reference: Vec<_> = part
            .body
            .vertices()
            .iter()
            .map(|v| v.position_m())
            .collect();
        assert_eq!(
            part.body.deform(&reference, || false).unwrap().inertia(),
            part.body.inertia()
        );
        let local: Vec<_> = part
            .original_vertex_indices
            .iter()
            .map(|&v| coordinates[v])
            .collect();
        let deformed = part.body.deform(&local, || false).unwrap();
        for (local, &source) in part.original_cell_indices.iter().enumerate() {
            assert_eq!(deformed.cells()[local], all.cells()[source]);
            assert_eq!(
                deformed.cells()[local].inertia().mass().to_bits(),
                body.cells()[source].inertia().mass().to_bits()
            );
        }
    }
}

#[test]
fn g4_duplicate_empty_budget_and_late_cancellation_refuse_without_changing_matter() {
    let body = refined();
    let original = body.clone();
    for selection in [vec![], vec![0, 0], vec![4], vec![0, 1, 2, 3]] {
        assert!(body.partition_cells(&selection, 4, 16, || false).is_err());
        assert_eq!(body, original);
    }
    assert!(body.partition_cells(&[1, 3], 3, 16, || false).is_err());
    // Five original vertices become two five-vertex meshes. Source count alone
    // cannot admit their combined retained payload.
    assert!(body.partition_cells(&[1, 3], 4, 9, || false).is_err());
    assert!(matches!(
        body.partition_cells(&[1, 3], 4, 16, || true),
        Err(MaterialInertiaError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        body.partition_cells(&[1, 3], 4, 16, || {
            polls += 1;
            polls > 24
        }),
        Err(MaterialInertiaError::Cancelled)
    ));
    assert!(polls > 24);
    assert_eq!(body, original);
    assert_eq!(
        body.partition_cells(&[1, 3], 4, 16, || false).unwrap(),
        body.partition_cells(&[3, 1], 4, 16, || false).unwrap()
    );
}
