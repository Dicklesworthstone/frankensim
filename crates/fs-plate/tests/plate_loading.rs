//! G0/G3: physical force/moment ports, interpolation, and reduced support maps.
use fs_plate::{
    AssemblyOptions, EdgeSupport, PlateMesh, PlateModel, PlateSection, assemble,
    loading::{PlateLoadBudget, PlatePointStencil},
};

fn mesh() -> PlateMesh {
    PlateMesh::from_unstructured(
        vec![(1.0, -1.0), (3.0, -1.0), (3.0, 2.0), (1.0, 2.0)],
        vec![[0, 1, 2], [0, 2, 3]],
    )
    .unwrap()
}

fn model(mesh: &PlateMesh, fixed: &[usize]) -> PlateModel {
    assemble(
        mesh,
        &PlateSection::isotropic(2e9, 0.3, 0.01, 800.0).unwrap(),
        fixed,
        &[],
        &AssemblyOptions {
            pretension: 0.0,
            support: EdgeSupport::Clamped,
        },
    )
    .unwrap()
}

fn budget() -> PlateLoadBudget {
    PlateLoadBudget {
        max_nodes: 4,
        max_triangles: 2,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-13 * (1.0 + expected.abs()),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

#[test]
fn rigid_translation_and_rotation_reconstruct_physical_motion() {
    let mesh = mesh();
    let model = model(&mesh, &[]);
    let point = [2.4, 0.2];
    let stencil = PlatePointStencil::locate(&mesh, &model, point, budget()).unwrap();
    assert_eq!(stencil.triangle(), 0);
    assert_eq!(stencil.point(), point);
    for (actual, expected) in stencil.weights().into_iter().zip([0.3, 0.3, 0.4]) {
        close(actual, expected);
    }
    let mut q = vec![0.0; model.free];
    let mut v = q.clone();
    for (node, &(x, y)) in mesh.nodes.iter().enumerate() {
        // w = translation + theta_x*y - theta_y*x; slopes are not rotations.
        let values = [0.7 + 0.3 * y - 0.2 * x, -0.2, 0.3];
        let rates = [-0.4 + 0.6 * y + 0.5 * x, 0.5, 0.6];
        for component in 0..3 {
            let index = model.dof_map[3 * node + component].unwrap();
            q[index] = values[component];
            v[index] = rates[component];
        }
    }
    let sample = stencil.sample(&q, &v).unwrap();
    close(sample.displacement, 0.7 + 0.3 * point[1] - 0.2 * point[0]);
    close(sample.slopes[0], -0.2);
    close(sample.slopes[1], 0.3);
    close(sample.velocity, -0.4 + 0.6 * point[1] + 0.5 * point[0]);
    for (actual, expected) in sample.angular_velocity.into_iter().zip([0.6, -0.5, 0.0]) {
        close(actual, expected);
    }
}

#[test]
fn physical_load_and_transpose_are_work_conjugate_for_arbitrary_nodal_motion() {
    let mesh = mesh();
    let model = model(&mesh, &[]);
    let stencil = PlatePointStencil::locate(&mesh, &model, [2.4, 0.2], budget()).unwrap();
    let q: Vec<f64> = (0..model.free).map(|i| (i as f64 - 4.0) * 0.13).collect();
    let v: Vec<f64> = (0..model.free)
        .map(|i| ((7 * i + 3) % 11) as f64 * 0.17 - 0.6)
        .collect();
    let load = [3.2, -1.7, 2.1];
    let mut nodal = vec![0.0; model.free];
    stencil.add_load(load, &mut nodal).unwrap();
    let sample = stencil.sample(&q, &v).unwrap();
    close(
        dot(&nodal, &v),
        load[0] * sample.velocity
            + load[1] * sample.angular_velocity[0]
            + load[2] * sample.angular_velocity[1],
    );
    close(
        dot(&nodal, &q),
        load[0] * sample.displacement + load[1] * sample.slopes[1] - load[2] * sample.slopes[0],
    );
    let transpose = stencil.load_vjp(&v).unwrap();
    close(dot(&nodal, &v), dot(&load, &transpose));
    close(transpose[0], sample.velocity);
    close(transpose[1], sample.angular_velocity[0]);
    close(transpose[2], sample.angular_velocity[1]);
}

#[test]
fn distributed_nodal_forces_and_couples_preserve_the_world_resultant() {
    let mesh = mesh();
    let model = model(&mesh, &[]);
    let point = [2.4, 0.2];
    let stencil = PlatePointStencil::locate(&mesh, &model, point, budget()).unwrap();
    let load = [4.0, 6.0, -2.0];
    let mut nodal = vec![0.0; model.free];
    let applied = stencil.add_load(load, &mut nodal).unwrap();
    assert_eq!(applied.fixed_dof_loads, Vec::new());
    let mut force = 0.0;
    let mut moment = [0.0; 3];
    for (node, &(x, y)) in mesh.nodes.iter().enumerate() {
        let f = nodal[model.dof_map[3 * node].unwrap()];
        force += f;
        moment[0] += y * f + nodal[model.dof_map[3 * node + 2].unwrap()];
        moment[1] -= x * f + nodal[model.dof_map[3 * node + 1].unwrap()];
    }
    close(force, load[0]);
    for (actual, expected) in applied.force.into_iter().zip([0.0, 0.0, load[0]]) {
        close(actual, expected);
    }
    let expected = [
        load[1] + point[1] * load[0],
        load[2] - point[0] * load[0],
        0.0,
    ];
    for i in 0..3 {
        close(moment[i], expected[i]);
        close(applied.moment[i], expected[i]);
    }
}

#[test]
fn a_moving_point_is_continuous_across_the_edge_and_selects_the_lowest_triangle() {
    let mesh = mesh();
    let model = model(&mesh, &[]);
    let points = [[2.0, 0.5 - 1e-8], [2.0, 0.5], [2.0, 0.5 + 1e-8]];
    let mut vectors = Vec::new();
    for (point, triangle) in points.into_iter().zip([0, 0, 1]) {
        let stencil = PlatePointStencil::locate(&mesh, &model, point, budget()).unwrap();
        assert_eq!(stencil.triangle(), triangle);
        close(stencil.weights().iter().sum(), 1.0);
        assert!(stencil.weights().iter().all(|weight| *weight >= 0.0));
        let mut output = vec![0.0; model.free];
        stencil.add_load([1.2, -0.3, 0.7], &mut output).unwrap();
        vectors.push(output);
    }
    for side in [0, 2] {
        assert!(
            vectors[side]
                .iter()
                .zip(&vectors[1])
                .all(|(a, b)| (a - b).abs() < 1e-7)
        );
    }
    assert_eq!(
        PlatePointStencil::locate(&mesh, &model, [1.0, -1.0], budget())
            .unwrap()
            .triangle(),
        0
    );
}

#[test]
fn supported_dofs_receive_separate_loads_and_sample_as_zero() {
    let mesh = mesh();
    let model = model(&mesh, &[0]);
    let stencil = PlatePointStencil::locate(&mesh, &model, [2.4, 0.2], budget()).unwrap();
    let mut output = vec![0.0; model.free];
    let applied = stencil.add_load([4.0, 6.0, -2.0], &mut output).unwrap();
    assert_eq!(
        applied
            .fixed_dof_loads
            .iter()
            .map(|(dof, _)| *dof)
            .collect::<Vec<_>>(),
        [0, 1, 2]
    );
    for ((_, value), expected) in applied.fixed_dof_loads.iter().zip([1.2, 0.6, 1.8]) {
        close(*value, expected);
    }
    for (node, weight) in [(1, 0.3), (2, 0.4), (3, 0.0)] {
        for (component, load) in [4.0, 2.0, 6.0].into_iter().enumerate() {
            close(
                output[model.dof_map[3 * node + component].unwrap()],
                weight * load,
            );
        }
    }
    let mut q = vec![0.0; model.free];
    let mut v = q.clone();
    for (full, reduced) in model.dof_map.iter().enumerate() {
        if let Some(index) = reduced {
            q[*index] = (full % 3 + 1) as f64;
            v[*index] = (full % 3 + 4) as f64;
        }
    }
    let sample = stencil.sample(&q, &v).unwrap();
    close(sample.displacement, 0.7);
    close(sample.slopes[0], 1.4);
    close(sample.slopes[1], 2.1);
    close(sample.velocity, 2.8);
    close(sample.angular_velocity[0], 4.2);
    close(sample.angular_velocity[1], -3.5);
    close(
        dot(&output, &v),
        dot(&[4.0, 6.0, -2.0], &stencil.load_vjp(&v).unwrap()),
    );
}

#[test]
fn admission_rejects_bad_points_meshes_maps_and_insufficient_budgets() {
    let mesh = mesh();
    let model = model(&mesh, &[]);
    let point = [2.4, 0.2];
    for point in [
        [f64::NAN, 0.2],
        [2.4, f64::INFINITY],
        [4.0, 0.5],
        [0.0, -1.0],
    ] {
        assert!(PlatePointStencil::locate(&mesh, &model, point, budget()).is_err());
    }
    for budget in [
        PlateLoadBudget {
            max_nodes: 3,
            ..budget()
        },
        PlateLoadBudget {
            max_triangles: 1,
            ..budget()
        },
    ] {
        assert!(PlatePointStencil::locate(&mesh, &model, point, budget).is_err());
    }
    let mut malformed = mesh.clone();
    malformed.tris[1] = [0, 0, 3];
    assert!(PlatePointStencil::locate(&malformed, &model, point, budget()).is_err());
    malformed = mesh.clone();
    malformed.tris[1][2] = mesh.nodes.len();
    assert!(PlatePointStencil::locate(&malformed, &model, point, budget()).is_err());
    malformed = mesh.clone();
    malformed.nodes[3].0 = f64::NAN;
    assert!(PlatePointStencil::locate(&malformed, &model, point, budget()).is_err());
    let mut malformed = model.clone();
    malformed.dof_map.pop();
    assert!(PlatePointStencil::locate(&mesh, &malformed, point, budget()).is_err());
    malformed = model.clone();
    malformed.dof_map[0] = Some(model.free);
    assert!(PlatePointStencil::locate(&mesh, &malformed, point, budget()).is_err());
    malformed = model.clone();
    malformed.dof_map[1] = malformed.dof_map[0];
    assert!(PlatePointStencil::locate(&mesh, &malformed, point, budget()).is_err());
}

#[test]
fn additive_loading_and_failures_never_corrupt_existing_loads() {
    let mesh = mesh();
    let model = model(&mesh, &[]);
    let stencil = PlatePointStencil::locate(&mesh, &model, [2.4, 0.2], budget()).unwrap();
    let mut combined = vec![0.0; model.free];
    let mut separate = combined.clone();
    stencil.add_load([2.0, 0.3, -0.7], &mut separate).unwrap();
    stencil.add_load([0.4, -0.2, 1.1], &mut separate).unwrap();
    stencil.add_load([2.4, 0.1, 0.4], &mut combined).unwrap();
    for (a, b) in combined.iter().zip(&separate) {
        close(*a, *b);
    }
    for load in [[f64::NAN, 0.0, 0.0], [0.0, f64::INFINITY, 0.0]] {
        let before = separate.clone();
        assert!(stencil.add_load(load, &mut separate).is_err());
        assert_eq!(separate, before);
    }
    let mut short = vec![0.25; model.free - 1];
    let before = short.clone();
    assert!(stencil.add_load([1.0; 3], &mut short).is_err());
    assert_eq!(short, before);
    separate[model.dof_map[8].unwrap()] = f64::MAX;
    let before = separate.clone();
    assert!(
        stencil
            .add_load([0.0, f64::MAX, 0.0], &mut separate)
            .is_err()
    );
    assert_eq!(separate, before);
    let mut bad = vec![0.0; model.free];
    bad[0] = f64::NAN;
    assert!(stencil.sample(&bad, &combined).is_err());
    assert!(stencil.sample(&combined, &bad).is_err());
    assert!(stencil.load_vjp(&bad).is_err());
    assert!(stencil.sample(&short, &combined).is_err());
    assert!(stencil.load_vjp(&short).is_err());
}
