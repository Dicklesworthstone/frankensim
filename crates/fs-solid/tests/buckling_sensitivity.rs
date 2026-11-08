//! Total sizing derivatives must follow the re-solved elastic stress field.
use fs_solid::SolidError;
use fs_solid::linear::{Formulation, LinearProblem, PlaneKind};
use fs_solid::mesh2::{Mesh2, Patch};
use fs_solid::stability::sensitivity::{
    element_scale_sensitivities, element_scale_sensitivities_cancellable,
};
use fs_solid::stability::{
    buckling_loads, eigenvalue_derivative, group_stiffness, reduced_pencil,
    reduced_pencil_with_element_scales,
};

fn problem(mesh: &Mesh2, formulation: Formulation, plane: PlaneKind) -> LinearProblem<'_> {
    LinearProblem {
        mesh,
        youngs: 1000.0,
        poisson: if plane == PlaneKind::Strain {
            0.499
        } else {
            0.28
        },
        plane,
        formulation,
        body_force: Some(&|_, _| [0.02, -0.03]),
        dirichlet: vec![(Patch::Left, &|_, _| [0.0, 0.0])],
        traction: vec![(Patch::Right, &|_, y| [-1.2 - 0.6 * y, 0.08])],
        symmetry: vec![],
    }
}

fn loads(problem: &LinearProblem<'_>, scales: &[f64], count: usize) -> Vec<f64> {
    let (k, g, map, _) = reduced_pencil_with_element_scales(problem, scales).unwrap();
    buckling_loads(&k, &g, &map, count + 1, 600).unwrap().loads
}

fn check_finite_differences(problem: &LinearProblem<'_>, scales: &[f64], count: usize) -> f64 {
    let result = element_scale_sensitivities(problem, scales, count, 600).unwrap();
    assert_eq!(result.element_scales(), scales);
    assert_eq!(result.buckling().loads.len(), count);
    assert!(
        result
            .adjoint_relative_residuals()
            .iter()
            .all(|r| *r <= 1e-8)
    );
    let (_, g, map, _) = reduced_pencil_with_element_scales(problem, scales).unwrap();
    let mut frozen_error = 0.0_f64;
    for element in 0..scales.len() {
        let step = 5e-5 * scales[element];
        let mut plus = scales.to_vec();
        let mut minus = scales.to_vec();
        plus[element] += step;
        minus[element] -= step;
        // Both endpoints independently re-solve the reference displacement,
        // recover the changed stress and solve the changed physical pencil.
        let positive = loads(problem, &plus, count);
        let negative = loads(problem, &minus, count);
        let dk = group_stiffness(problem, &map, &|e| e == element);
        for mode in 0..count {
            let fd = (positive[mode] - negative[mode]) / (2.0 * step);
            let gradient = result.gradients()[mode][element];
            let scale = fd.abs().max(1e-3 * result.buckling().loads[mode]);
            assert!(
                (gradient - fd).abs() <= 3e-5 * scale,
                "{form:?}/{plane:?} element {element}, mode {mode}: total={gradient:e}, FD={fd:e}, relative={:e}",
                (gradient - fd).abs() / scale,
                form = problem.formulation,
                plane = problem.plane
            );
            let frozen = eigenvalue_derivative(&result.buckling().modes[mode], &dk, &g);
            frozen_error = frozen_error.max((frozen - fd).abs() / result.buckling().loads[mode]);
        }
    }
    frozen_error
}

#[test]
fn heterogeneous_p1_q1_total_gradients_match_recomputed_equilibrium() {
    let quads = Mesh2::mapped_quads(2, 2, &|s, t| {
        [3.0 * s + 0.15 * s * t, 0.8 * t + 0.06 * s * t]
    });
    let triangles = Mesh2::triangles(3.0, 0.8, 2, 1);
    let scales = [0.65, 0.9, 1.25, 1.7];
    let mut frozen_error = 0.0_f64;
    for mesh in [&quads, &triangles] {
        for form in [Formulation::Standard, Formulation::BBar] {
            for plane in [PlaneKind::Stress, PlaneKind::Strain] {
                frozen_error = frozen_error.max(check_finite_differences(
                    &problem(mesh, form, plane),
                    &scales,
                    1,
                ));
            }
        }
    }
    assert!(
        frozen_error > 1e-3,
        "fixture must distinguish the missing equilibrium-adjoint chain: {frozen_error:e}"
    );
}

#[test]
fn uniform_scaling_respects_force_and_displacement_loading_homogeneity() {
    let mesh = Mesh2::quads(3.0, 0.7, 3, 1);
    let scales = [0.65, 1.1, 1.6];
    let scaled: Vec<_> = scales.iter().map(|s| 1.7 * s).collect();
    for formulation in [Formulation::Standard, Formulation::BBar] {
        let mut p = problem(&mesh, formulation, PlaneKind::Stress);
        let force = element_scale_sensitivities(&p, &scales, 1, 600).unwrap();
        let derivative: f64 = force.gradients()[0]
            .iter()
            .zip(scales)
            .map(|(g, s)| g * s)
            .sum();
        let load = force.buckling().loads[0];
        assert!((derivative - load).abs() < 2e-7 * load);
        assert!((loads(&p, &scaled, 1)[0] / load - 1.7).abs() < 2e-7);
        let u = p.solve_with_element_scales(&scales).unwrap();
        let stiff = p.solve_with_element_scales(&scaled).unwrap();
        for (a, b) in u.iter().flatten().zip(stiff.iter().flatten()) {
            assert!((a - 1.7 * b).abs() < 1e-10);
        }

        p.body_force = None;
        p.traction.clear();
        p.dirichlet
            .push((Patch::Right, &|_, y| [-0.015, 0.001 * y]));
        let imposed = element_scale_sensitivities(&p, &scales, 1, 600).unwrap();
        let derivative: f64 = imposed.gradients()[0]
            .iter()
            .zip(scales)
            .map(|(g, s)| g * s)
            .sum();
        let load = imposed.buckling().loads[0];
        assert!(
            derivative.abs() < 2e-7 * load,
            "uniform modulus scaling cannot change prescribed-load buckling: {derivative:e}"
        );
        assert!((loads(&p, &scaled, 1)[0] / load - 1.0).abs() < 2e-7);
        let u = p.solve_with_element_scales(&scales).unwrap();
        let stiff = p.solve_with_element_scales(&scaled).unwrap();
        for (a, b) in u.iter().flatten().zip(stiff.iter().flatten()) {
            assert!((a - b).abs() < 1e-10);
        }
    }
}

#[test]
fn mixed_prescribed_and_force_loading_retains_the_fixed_column_derivative() {
    let mesh = Mesh2::mapped_quads(3, 2, &|s, t| [3.0 * s + 0.1 * s * t, 0.8 * t]);
    let scales = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6];
    for formulation in [Formulation::Standard, Formulation::BBar] {
        let mut p = problem(&mesh, formulation, PlaneKind::Stress);
        p.dirichlet
            .push((Patch::Right, &|_, y| [-0.01 * (1.0 + 0.1 * y), 0.002 * y]));
        p.traction = vec![(Patch::Top, &|_, _| [-0.3, -0.2])];
        assert!(check_finite_differences(&p, &scales, 2) > 1e-3);
        let u = p.solve_with_element_scales(&scales).unwrap();
        for (patch, field) in &p.dirichlet {
            for node in mesh.patch_nodes(*patch) {
                assert_eq!(u[node], field(mesh.nodes[node][0], mesh.nodes[node][1]));
            }
        }
    }
}

#[test]
fn scaled_forward_path_preserves_legacy_equilibrium_and_rejects_invalid_designs() {
    let mesh = Mesh2::quads(2.0, 0.8, 2, 1);
    let p = problem(&mesh, Formulation::BBar, PlaneKind::Strain);
    let (k, g, map, u) = reduced_pencil(&p).unwrap();
    let (ks, gs, maps, us) = reduced_pencil_with_element_scales(&p, &[1.0, 1.0]).unwrap();
    assert_eq!(k.to_dense(), ks.to_dense());
    assert_eq!(g.to_dense(), gs.to_dense());
    assert_eq!(map, maps);
    assert_eq!(u, us);
    for scales in [
        vec![1.0],
        vec![1.0, 0.0],
        vec![-1.0, 1.0],
        vec![1.0, f64::NAN],
        vec![f64::INFINITY, 1.0],
        vec![f64::MAX, 1.0],
    ] {
        assert!(matches!(
            p.solve_with_element_scales(&scales),
            Err(SolidError::InvalidInput { .. })
        ));
    }
    let mut invalid_support = problem(&mesh, Formulation::BBar, PlaneKind::Strain);
    invalid_support.symmetry = vec![(Patch::Bottom, 2)];
    assert!(matches!(
        invalid_support.solve_with_element_scales(&[1.0, 1.0]),
        Err(SolidError::InvalidInput { .. })
    ));
    assert!(matches!(
        element_scale_sensitivities(&p, &[1.0, 1.0], 0, 10),
        Err(SolidError::InvalidInput { .. })
    ));
    assert!(matches!(
        element_scale_sensitivities_cancellable(&p, &[1.0, 1.0], 1, 100, &|| true),
        Err(SolidError::Cancelled)
    ));
}

#[test]
fn repeated_modes_refuse_an_individual_branch_gradient() {
    // Two physically disconnected, identically loaded columns have an exact
    // double first eigenvalue. Requesting only one must not hide its neighbor.
    let mut mesh = Mesh2::quads(2.0, 0.7, 1, 1);
    let duplicate = mesh.clone();
    let offset = mesh.nodes.len();
    mesh.nodes
        .extend(duplicate.nodes.iter().map(|p| [p[0], p[1] + 2.0]));
    mesh.elems.extend(
        duplicate
            .elems
            .iter()
            .map(|e| e.iter().map(|n| n + offset).collect()),
    );
    for ((_, edges), (_, additional)) in mesh.patches.iter_mut().zip(&duplicate.patches) {
        edges.extend(additional.iter().map(|(a, b)| (a + offset, b + offset)));
    }
    let mut p = problem(&mesh, Formulation::Standard, PlaneKind::Stress);
    p.poisson = 0.0;
    p.body_force = None;
    p.traction = vec![(Patch::Right, &|_, _| [-1.0, 0.0])];
    assert!(
        matches!(element_scale_sensitivities(&p, &[1.0, 1.0], 1, 600),
        Err(SolidError::InvalidInput { what }) if what.contains("separated eigenvalues"))
    );
}
