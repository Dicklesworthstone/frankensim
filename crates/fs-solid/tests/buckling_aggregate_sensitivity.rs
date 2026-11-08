//! Complete spectral packets have a smooth total derivative at mode crossings.
use fs_solid::SolidError;
use fs_solid::linear::{Formulation, LinearProblem, PlaneKind};
use fs_solid::mesh2::{Mesh2, Patch};
use fs_solid::stability::sensitivity::aggregate::{
    ElementScaleAggregateSensitivity, element_scale_aggregate_sensitivities,
    element_scale_aggregate_sensitivities_cancellable,
};
use fs_solid::stability::sensitivity::element_scale_sensitivities;
use fs_solid::stability::{buckling_loads, reduced_pencil_with_element_scales};

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

fn columns(count: usize) -> Mesh2 {
    let mut mesh = Mesh2::quads(2.0, 0.7, 1, 1);
    let unit = mesh.clone();
    for column in 1..count {
        let offset = mesh.nodes.len();
        mesh.nodes.extend(
            unit.nodes
                .iter()
                .map(|p| [p[0], p[1] + 2.0 * column as f64]),
        );
        mesh.elems.extend(
            unit.elems
                .iter()
                .map(|e| e.iter().map(|n| n + offset).collect()),
        );
        for ((_, edges), (_, additional)) in mesh.patches.iter_mut().zip(&unit.patches) {
            edges.extend(additional.iter().map(|(a, b)| (a + offset, b + offset)));
        }
    }
    mesh
}

fn column_problem(mesh: &Mesh2) -> LinearProblem<'_> {
    let mut p = problem(mesh, Formulation::Standard, PlaneKind::Stress);
    p.poisson = 0.0;
    p.body_force = None;
    p.traction = vec![(Patch::Right, &|_, _| [-1.0, 0.0])];
    p
}

fn aggregate(
    problem: &LinearProblem<'_>,
    scales: &[f64],
    count: usize,
    sharpness: f64,
) -> ElementScaleAggregateSensitivity {
    let result =
        element_scale_aggregate_sensitivities(problem, scales, count, 600, sharpness).unwrap();
    assert_eq!(result.element_scales(), scales);
    assert_eq!(result.mode_count(), count);
    assert_eq!(result.buckling().loads.len(), count);
    assert_eq!(result.buckling().modes.len(), count);
    assert_eq!(result.sharpness(), sharpness);
    assert!((result.weights().iter().sum::<f64>() - 1.0).abs() < 1e-14);
    assert!(result.adjoint_relative_residual() <= 1e-8);
    assert!(
        result
            .buckling()
            .relative_residuals
            .iter()
            .all(|r| *r <= 1e-8)
    );
    result
}

fn ks(loads: &[f64], sharpness: f64) -> f64 {
    let minimum = loads.iter().copied().fold(f64::INFINITY, f64::min);
    minimum
        - loads
            .iter()
            .map(|load| (-sharpness * (load - minimum)).exp())
            .sum::<f64>()
            .ln()
            / sharpness
}

fn recomputed_value(
    problem: &LinearProblem<'_>,
    scales: &[f64],
    count: usize,
    sharpness: f64,
) -> f64 {
    // Independent endpoints re-solve equilibrium, recover stresses and solve
    // the physical buckling pencil. No accepted modes or adjoints are reused.
    let (k, g, map, _) = reduced_pencil_with_element_scales(problem, scales).unwrap();
    let physical = buckling_loads(&k, &g, &map, count + 1, 600).unwrap();
    ks(&physical.loads[..count], sharpness)
}

fn check_finite_differences(
    problem: &LinearProblem<'_>,
    scales: &[f64],
    count: usize,
    sharpness: f64,
) -> ElementScaleAggregateSensitivity {
    let result = aggregate(problem, scales, count, sharpness);
    for element in 0..scales.len() {
        let step = 5e-5 * scales[element];
        let mut plus = scales.to_vec();
        let mut minus = scales.to_vec();
        plus[element] += step;
        minus[element] -= step;
        let fd = (recomputed_value(problem, &plus, count, sharpness)
            - recomputed_value(problem, &minus, count, sharpness))
            / (2.0 * step);
        let scale = fd.abs().max(1e-3 * result.buckling().loads[0]);
        assert!(
            (result.gradient()[element] - fd).abs() <= 4e-5 * scale,
            "{form:?}/{plane:?}, element {element}: aggregate={:e}, FD={fd:e}",
            result.gradient()[element],
            form = problem.formulation,
            plane = problem.plane,
        );
    }
    result
}

#[test]
fn complete_double_mode_packet_matches_analytic_trace_and_crossing_derivatives() {
    let mesh = columns(2);
    let p = column_problem(&mesh);
    let initial = aggregate(&p, &[1.0, 1.0], 2, 1.0);
    let first = initial.buckling().loads[0];
    assert!((initial.buckling().loads[1] / first - 1.0).abs() < 1e-10);
    let sharpness = 8.0 / first;
    for scales in [[1.0, 1.0], [0.93, 1.07], [1.07, 0.93]] {
        let result = check_finite_differences(&p, &scales, 2, sharpness);
        // Each disconnected force-loaded column's physical load is exactly
        // first*scale. These analytic branches give an independent aggregate
        // value and derivative, including the exact two-dimensional cluster.
        let exact_loads = [first * scales[0], first * scales[1]];
        let minimum = exact_loads[0].min(exact_loads[1]);
        let terms = exact_loads.map(|load| (-sharpness * (load - minimum)).exp());
        let sum = terms.iter().sum::<f64>();
        assert!((result.value() - ks(&exact_loads, sharpness)).abs() < 1e-9 * first);
        for element in 0..2 {
            let exact = first * terms[element] / sum;
            assert!((result.gradient()[element] - exact).abs() < 1e-8 * first);
        }
    }
    let at_crossing = aggregate(&p, &[1.0, 1.0], 2, sharpness);
    assert!((at_crossing.value() - (first - 2.0_f64.ln() / sharpness)).abs() < 1e-9 * first);
    assert!((at_crossing.gradient()[0] - at_crossing.gradient()[1]).abs() < 1e-9 * first);

    // Reversing physical element numbering permutes the gradient, even though
    // the native eigensolver may choose a different degenerate eigenbasis.
    let mut reordered_mesh = mesh.clone();
    reordered_mesh.elems.reverse();
    let reordered = aggregate(
        &column_problem(&reordered_mesh),
        &[1.07, 0.93],
        2,
        sharpness,
    );
    let original = aggregate(&p, &[0.93, 1.07], 2, sharpness);
    assert!((reordered.value() - original.value()).abs() < 1e-9 * first);
    assert!((reordered.gradient()[0] - original.gradient()[1]).abs() < 1e-8 * first);
    assert!((reordered.gradient()[1] - original.gradient()[0]).abs() < 1e-8 * first);
}

#[test]
fn unequal_p1_q1_packets_follow_total_equilibrium_and_native_bbar_transpose() {
    let quads = Mesh2::mapped_quads(2, 2, &|s, t| {
        [3.0 * s + 0.15 * s * t, 0.8 * t + 0.06 * s * t]
    });
    let triangles = Mesh2::triangles(3.0, 0.8, 2, 1);
    let scales = [0.65, 0.9, 1.25, 1.7];
    for mesh in [&quads, &triangles] {
        for form in [Formulation::Standard, Formulation::BBar] {
            for plane in [PlaneKind::Stress, PlaneKind::Strain] {
                let p = problem(mesh, form, plane);
                let separated = element_scale_sensitivities(&p, &scales, 2, 600).unwrap();
                let sharpness = 2.0 / separated.buckling().loads[1];
                let result = check_finite_differences(&p, &scales, 2, sharpness);
                for element in 0..scales.len() {
                    let weighted: f64 = result
                        .weights()
                        .iter()
                        .zip(separated.gradients())
                        .map(|(weight, gradient)| weight * gradient[element])
                        .sum();
                    assert!(
                        (result.gradient()[element] - weighted).abs()
                            <= 2e-7 * result.buckling().loads[0]
                    );
                }
                // Under uniform modulus scaling with dead forces, every
                // physical eigenvalue scales linearly. The derivative is the
                // weighted physical load, not the biased aggregate value.
                let dilation: f64 = result
                    .gradient()
                    .iter()
                    .zip(scales)
                    .map(|(g, s)| g * s)
                    .sum();
                let expected: f64 = result
                    .weights()
                    .iter()
                    .zip(&result.buckling().loads)
                    .map(|(w, load)| w * load)
                    .sum();
                assert!((dilation - expected).abs() < 2e-7 * expected);
                assert!((expected - result.value()).abs() > 1e-3 * expected);
            }
        }
    }
}

#[test]
fn prescribed_columns_and_mixed_forces_survive_the_single_aggregate_adjoint() {
    let mesh = Mesh2::mapped_quads(3, 2, &|s, t| [3.0 * s + 0.1 * s * t, 0.8 * t]);
    let scales = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6];
    for formulation in [Formulation::Standard, Formulation::BBar] {
        let mut p = problem(&mesh, formulation, PlaneKind::Stress);
        p.dirichlet
            .push((Patch::Right, &|_, y| [-0.01 * (1.0 + 0.1 * y), 0.002 * y]));
        p.traction = vec![(Patch::Top, &|_, _| [-0.3, -0.2])];
        let result = check_finite_differences(&p, &scales, 2, 0.03);
        assert!(result.adjoint_iterations() > 0);

        p.body_force = None;
        p.traction.clear();
        let imposed = check_finite_differences(&p, &scales, 2, 0.03);
        let dilation: f64 = imposed
            .gradient()
            .iter()
            .zip(scales)
            .map(|(g, s)| g * s)
            .sum();
        assert!(dilation.abs() < 2e-7 * imposed.buckling().loads[0]);
        let scaled: Vec<_> = scales.iter().map(|scale| 1.7 * scale).collect();
        let changed = aggregate(&p, &scaled, 2, 0.03);
        assert!((changed.value() - imposed.value()).abs() < 2e-7 * imposed.buckling().loads[0]);
    }
}

#[test]
fn packet_boundaries_cannot_cut_repeated_modes_and_shortfall_is_not_silently_truncated() {
    let mesh = columns(2);
    let p = column_problem(&mesh);
    assert!(matches!(
        element_scale_aggregate_sensitivities(&p, &[1.0, 1.0], 1, 600, 0.1),
        Err(SolidError::InvalidInput { what }) if what.contains("packet boundary")
    ));
    assert!(matches!(
        element_scale_sensitivities(&p, &[1.0, 1.0], 1, 600),
        Err(SolidError::InvalidInput { what }) if what.contains("separated eigenvalues")
    ));
    aggregate(&p, &[1.0, 1.0], 2, 0.1);
    assert!(matches!(
        element_scale_aggregate_sensitivities(&p, &[1.0, 1.0], 9, 600, 0.1),
        Err(SolidError::InvalidInput { what }) if what.contains("complete requested")
    ));
    let mut unloaded = column_problem(&mesh);
    unloaded.traction.clear();
    assert!(matches!(
        element_scale_aggregate_sensitivities(&unloaded, &[1.0, 1.0], 2, 600, 0.1),
        Err(SolidError::InvalidInput { what }) if what.contains("complete requested")
    ));
}

#[test]
fn shifted_ks_extreme_sharpness_and_reference_units_remain_finite() {
    let mesh = columns(2);
    let p = column_problem(&mesh);
    let scales = [0.9, 1.1];
    let baseline = aggregate(&p, &scales, 2, 0.01);
    let sharp = aggregate(&p, &scales, 2, f64::MAX);
    assert_eq!(sharp.value(), sharp.buckling().loads[0]);
    assert_eq!(sharp.weights(), &[1.0, 0.0]);
    let broad = aggregate(&p, &scales, 2, 1e-10);
    assert!(broad.value() < 0.0);

    // Changing reference-force units must scale physical eigenvalues and total
    // gradients reciprocally. The tiny-load case would overflow an explicitly
    // formed lambda_i*q_i, which is avoided by the scaled single adjoint.
    for force_scale in [1e-170, 1e170] {
        let traction = |_: f64, _: f64| [-force_scale, 0.0];
        let mut scaled = column_problem(&mesh);
        scaled.traction = vec![(Patch::Right, &traction)];
        let changed = aggregate(&scaled, &scales, 2, baseline.sharpness() * force_scale);
        let load_scale = baseline.buckling().loads[0];
        assert!((changed.value() * force_scale - baseline.value()).abs() < 2e-7 * load_scale);
        for element in 0..2 {
            assert!(
                (changed.gradient()[element] * force_scale - baseline.gradient()[element]).abs()
                    < 2e-7 * load_scale
            );
        }
    }
}

#[test]
fn repeated_packets_refuse_sharpness_that_resolves_native_roundoff() {
    let mesh = columns(2);
    let p = column_problem(&mesh);
    let accepted = aggregate(&p, &[1.0, 1.0], 2, 0.01);
    assert!((accepted.weights()[0] - 0.5).abs() < 1e-10);
    assert!(matches!(
        element_scale_aggregate_sensitivities(&p, &[1.0, 1.0], 2, 600, f64::MAX),
        Err(SolidError::InvalidInput { what }) if what.contains("native residual/roundoff resolution")
    ));
    // A repeated higher cluster carries no aggregate weight at this sharpness
    // and must not prevent a well-separated lowest branch from being accepted.
    let three = columns(3);
    let higher = aggregate(&column_problem(&three), &[0.9, 1.1, 1.1], 3, f64::MAX);
    assert_eq!(higher.weights(), &[1.0, 0.0, 0.0]);
    assert_eq!(higher.value(), higher.buckling().loads[0]);
}

#[test]
fn invalid_policies_and_cancellation_publish_no_partial_aggregate() {
    let mesh = columns(2);
    let p = column_problem(&mesh);
    for sharpness in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e-310] {
        assert!(matches!(
            element_scale_aggregate_sensitivities(&p, &[1.0, 1.0], 2, 600, sharpness),
            Err(SolidError::InvalidInput { .. })
        ));
    }
    for (count, iterations) in [(0, 600), (2, 0)] {
        assert!(matches!(
            element_scale_aggregate_sensitivities(&p, &[1.0, 1.0], count, iterations, 0.1),
            Err(SolidError::InvalidInput { .. })
        ));
    }
    assert!(matches!(
        element_scale_aggregate_sensitivities_cancellable(&p, &[1.0, 1.0], 2, 600, 0.1, &|| true),
        Err(SolidError::Cancelled)
    ));

    let checkpoints = std::cell::Cell::new(0usize);
    let accepted =
        element_scale_aggregate_sensitivities_cancellable(&p, &[1.0, 1.0], 2, 600, 0.1, &|| {
            checkpoints.set(checkpoints.get() + 1);
            false
        })
        .unwrap();
    let total = checkpoints.get();
    assert!(total > 5);
    let expected = accepted.gradient().to_vec();
    for limit in [total / 2, total - 2] {
        let calls = std::cell::Cell::new(0usize);
        assert!(matches!(
            element_scale_aggregate_sensitivities_cancellable(
                &p,
                &[1.0, 1.0],
                2,
                600,
                0.1,
                &|| {
                    calls.set(calls.get() + 1);
                    calls.get() >= limit
                },
            ),
            Err(SolidError::Cancelled)
        ));
        assert_eq!(accepted.gradient(), expected);
    }
}
