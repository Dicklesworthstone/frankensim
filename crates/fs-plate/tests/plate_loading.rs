//! G0/G3: physical force/moment ports, interpolation, and reduced support maps.
use fs_plate::{
    AssemblyOptions, EdgeSupport, PlateMesh, PlateModel, PlateSection, assemble,
    loading::{PlateLoadBudget, PlatePointStencil, PlateTraceStencil},
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

#[test]
fn g0_original_trace_preserves_affine_motion_and_transpose_point_work() {
    let mesh = mesh();
    let model = model(&mesh, &[]);
    let velocity = [0.13, -0.27, 0.41, 0.09];
    let mut nodal_force = [0.0; 4];
    let mut physical_power = 0.0;
    for (point, force) in [([2.4, 0.2], 1.7), ([2.0, 0.5], -0.8), ([1.0, 2.0], 0.3)] {
        let trace = PlateTraceStencil::locate(&mesh, point, budget()).unwrap();
        let original = PlatePointStencil::locate(&mesh, &model, point, budget()).unwrap();
        assert_eq!(trace.triangle(), original.triangle());
        assert_eq!(trace.weights(), original.weights());
        let mut motion = 0.0;
        let mut point_velocity = 0.0;
        for (node, weight) in trace.nodes().into_iter().zip(trace.weights()) {
            let (x, y) = mesh.nodes[node];
            motion += weight * (0.7 + 0.3 * y - 0.2 * x);
            point_velocity += weight * velocity[node];
            nodal_force[node] += weight * force;
        }
        close(motion, 0.7 + 0.3 * point[1] - 0.2 * point[0]);
        physical_power += force * point_velocity;
    }
    close(dot(&nodal_force, &velocity), physical_power);
}

use fs_plate::loading::{PlatePatchBudget, PlatePatchPoint, PlatePatchStencil};

fn patch_points() -> [PlatePatchPoint; 3] {
    [
        PlatePatchPoint {
            point: [2.5, -0.5],
            area_m2: 1.0,
        },
        PlatePatchPoint {
            point: [1.5, 1.5],
            area_m2: 2.0,
        },
        PlatePatchPoint {
            point: [1.0, -1.0],
            area_m2: 1.0,
        },
    ]
}

fn patch_budget() -> PlatePatchBudget {
    PlatePatchBudget {
        max_points: 3,
        location: budget(),
    }
}

#[test]
fn finite_patch_crosses_triangles_and_supports_with_conservative_resultants_and_power() {
    let mesh = mesh();
    let model = model(&mesh, &[0]);
    let points = patch_points();
    let patch = PlatePatchStencil::locate(&mesh, &model, &points, patch_budget()).unwrap();
    assert_eq!(patch.area_m2(), 4.0);
    assert_eq!(patch.point_count(), 3);
    let load = [8.0, 4.0, -12.0];
    let mut reduced = vec![0.0; model.free];
    let applied = patch.add_load(load, &mut reduced).unwrap();
    let mut full = vec![0.0; model.dof_map.len()];
    for (i, index) in model.dof_map.iter().enumerate() {
        if let Some(index) = index {
            full[i] = reduced[*index];
        }
    }
    assert_eq!(
        applied
            .fixed_dof_loads
            .iter()
            .map(|entry| entry.0)
            .collect::<Vec<_>>(),
        [0, 1, 2]
    );
    for (index, value) in &applied.fixed_dof_loads {
        full[*index] = *value;
    }
    let mut resultant = [0.0; 3];
    for (node, &(x, y)) in mesh.nodes.iter().enumerate() {
        resultant[0] += full[3 * node];
        resultant[1] += y * full[3 * node] + full[3 * node + 2];
        resultant[2] += -x * full[3 * node] - full[3 * node + 1];
    }
    let centroid = [1.625, 0.375];
    close(resultant[0], load[0]);
    close(resultant[1], load[1] + centroid[1] * load[0]);
    close(resultant[2], load[2] - centroid[0] * load[0]);
    close(applied.force[2], resultant[0]);
    close(applied.moment[0], resultant[1]);
    close(applied.moment[1], resultant[2]);
    let q: Vec<f64> = (0..model.free).map(|i| 0.13 * i as f64 - 0.2).collect();
    let v: Vec<f64> = (0..model.free)
        .map(|i| ((i * 5) % 7) as f64 * 0.17)
        .collect();
    let sample = patch.sample(&q, &v).unwrap();
    close(
        dot(&reduced, &q),
        load[0] * sample.displacement + load[1] * sample.slopes[1] - load[2] * sample.slopes[0],
    );
    close(dot(&reduced, &v), patch.power(load, &v).unwrap());
    close(dot(&reduced, &q), dot(&load, &patch.load_vjp(&q).unwrap()));
    let mut pressure_load = vec![0.0; model.free];
    let pressure_applied = patch
        .add_tractions(&[[2.0, 1.0, -3.0]; 3], &mut pressure_load)
        .unwrap();
    assert_eq!(reduced, pressure_load);
    assert_eq!(applied, pressure_applied);
}

#[test]
fn nonuniform_patch_tractions_preserve_physical_power_and_quadrature_subdivision() {
    let mesh = mesh();
    let model = model(&mesh, &[]);
    let points = patch_points();
    let patch = PlatePatchStencil::locate(&mesh, &model, &points, patch_budget()).unwrap();
    let tractions = [[2.0, 1.0, -0.5], [-0.7, 0.4, 1.2], [3.1, -0.3, 0.9]];
    let v: Vec<f64> = (0..model.free)
        .map(|i| (i * 7 % 13) as f64 * 0.13 - 0.5)
        .collect();
    let mut nodal = vec![0.0; model.free];
    let applied = patch.add_tractions(&tractions, &mut nodal).unwrap();
    let mut force = 0.0;
    let mut moment = [0.0; 2];
    let mut power = 0.0;
    for (point, traction) in points.iter().zip(tractions) {
        let load = traction.map(|value| value * point.area_m2);
        force += load[0];
        moment[0] += load[1] + point.point[1] * load[0];
        moment[1] += load[2] - point.point[0] * load[0];
        let sample = PlatePointStencil::locate(&mesh, &model, point.point, budget())
            .unwrap()
            .sample(&v, &v)
            .unwrap();
        power += load[0] * sample.velocity
            + load[1] * sample.angular_velocity[0]
            + load[2] * sample.angular_velocity[1];
    }
    close(applied.force[2], force);
    close(applied.moment[0], moment[0]);
    close(applied.moment[1], moment[1]);
    close(dot(&nodal, &v), power);
    let split: Vec<_> = points
        .iter()
        .flat_map(|point| {
            [PlatePatchPoint {
                point: point.point,
                area_m2: point.area_m2 / 2.0,
            }; 2]
        })
        .collect();
    let divided = PlatePatchStencil::locate(
        &mesh,
        &model,
        &split,
        PlatePatchBudget {
            max_points: 6,
            ..patch_budget()
        },
    )
    .unwrap();
    let load = [2.3, -1.2, 0.4];
    let mut original = vec![0.0; model.free];
    let mut subdivided = original.clone();
    patch.add_load(load, &mut original).unwrap();
    divided.add_load(load, &mut subdivided).unwrap();
    for (&a, &b) in original.iter().zip(&subdivided) {
        close(a, b);
    }
    close(
        patch.power(load, &v).unwrap(),
        divided.power(load, &v).unwrap(),
    );
}

#[test]
fn patch_admission_and_late_failure_preserve_the_entire_target() {
    let mesh = mesh();
    let model = model(&mesh, &[0]);
    let points = patch_points();
    assert!(PlatePatchStencil::locate(&mesh, &model, &[], patch_budget()).is_err());
    assert!(
        PlatePatchStencil::locate(
            &mesh,
            &model,
            &points,
            PlatePatchBudget {
                max_points: 2,
                ..patch_budget()
            }
        )
        .is_err()
    );
    for area in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut bad = points;
        bad[2].area_m2 = area;
        assert!(PlatePatchStencil::locate(&mesh, &model, &bad, patch_budget()).is_err());
    }
    let mut outside = points;
    outside[2].point = [20.0, 0.0];
    assert!(PlatePatchStencil::locate(&mesh, &model, &outside, patch_budget()).is_err());
    let patch = PlatePatchStencil::locate(&mesh, &model, &points, patch_budget()).unwrap();
    let mut reduced = vec![0.3; model.free];
    let before = reduced.clone();
    // The first two points have already modified private scratch before the
    // final (fully supported) point refuses its nonfinite moment.
    assert!(
        patch
            .add_tractions(
                &[[2.0; 3], [3.0; 3], [0.0, f64::INFINITY, 0.0]],
                &mut reduced
            )
            .is_err()
    );
    assert_eq!(reduced, before);
    assert!(patch.add_tractions(&[[2.0; 3]], &mut reduced).is_err());
    assert_eq!(reduced, before);
    assert!(patch.add_load([f64::NAN, 0.0, 0.0], &mut reduced).is_err());
    assert_eq!(reduced, before);
    assert!(patch.power([f64::INFINITY; 3], &reduced).is_err());
    assert!(patch.sample(&reduced[..1], &reduced).is_err());
}
