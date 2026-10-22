//! M2.5-W5 plain CBIE. The iterative owner keeps the dense owner's triangle
//! quadrature, so below the dense cap it must equal the dense plain-CBIE solve
//! to the FMM's accuracy. Run release on ts1 or trj.
use fs_bem::helmholtz::iterative::{
    CBIE_MIN_LEAF_RADII, IterativeOptions, IterativeRadiation, solve_radiation_iterative,
};
use fs_bem::helmholtz::{Formulation, Medium, solve_radiation, with_dense_workers};
use fs_bem::panel3d::SpherePanels;
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};

fn velocity(surface: &SpherePanels) -> Vec<C64> {
    surface
        .normals()
        .iter()
        .enumerate()
        .map(|(i, n)| C64::new(1.0 + 0.2 * n[0] + 0.03 * det::sin(i as f64), 0.1 * n[1]))
        .collect()
}
fn relative(a: &[C64], b: &[C64]) -> f64 {
    det::sqrt(
        a.iter()
            .zip(b)
            .map(|(&a, &b)| (a - b).norm_sq())
            .sum::<f64>()
            / b.iter().map(|b| b.norm_sq()).sum::<f64>(),
    )
}
fn bits(p: &[C64]) -> Vec<(u64, u64)> {
    p.iter().map(|z| (z.re.to_bits(), z.im.to_bits())).collect()
}
fn workers() -> usize {
    std::thread::available_parallelism().map_or(1, |n| n.get().min(16))
}
/// Icosphere-4 stretched to semi-axes 1, 0.7 and 0.5 m: no spherical symmetry.
fn ellipsoid() -> SpherePanels {
    let sphere = SpherePanels::icosphere(1.0, 4).unwrap();
    let triangles = sphere
        .triangles()
        .unwrap()
        .iter()
        .map(|t| t.map(|v| [v[0], 0.7 * v[1], 0.5 * v[2]]))
        .collect();
    SpherePanels::from_triangles(triangles).unwrap()
}

/// G3: iterative plain CBIE against the dense plain-CBIE solve, at 5120
/// panels, on two tree depths with order 7: leaf capacity 128 panels gives
/// depth 2 (0.5 m leaves), 32 gives depth 3 (0.25 m, 5.2 panel radii).
/// The bound is the Burton-Miller owner's 10x its 1e-6 GMRES tolerance
/// (bem_iterative_matches_dense); GMRES here converges to 1e-10 so the FMM
/// and the far panels' seven-point rule set the difference.
#[test]
fn iterative_cbie_matches_dense_cbie() {
    let gate = CancelGate::new_clock_free();
    let medium = Medium::air();
    let sphere = SpherePanels::icosphere(1.0, 4).unwrap();
    let stretched = ellipsoid();
    // k < 0.8 * pi * sqrt(3) / 2 = 2.18, below the sphere's plain-CBIE limit.
    let cases = [
        ("sphere", &sphere, 0.05, 128),
        ("sphere", &sphere, 1.0, 32),
        ("sphere", &sphere, 2.0, 32),
        ("ellipsoid", &stretched, 1.0, 32),
    ];
    for (name, surface, k, leaf_capacity) in cases {
        let v = velocity(surface);
        let dense = with_dense_workers(workers(), || {
            solve_radiation(surface, k, medium, &v, Formulation::PlainCbie)
        })
        .unwrap();
        let o = IterativeOptions {
            tol: 1e-10,
            fmm_order: 7,
            leaf_capacity,
            max_near_entries: 1 << 26,
            workers: workers(),
            ..IterativeOptions::default()
        };
        let (actual, report, receipt) =
            solve_radiation_iterative(surface, k, medium, &v, Formulation::PlainCbie, &o, &gate)
                .unwrap();
        let error = relative(&actual.pressure, &dense.pressure);
        let power = (actual.radiated_power - dense.radiated_power).abs() / dense.radiated_power;
        println!(
            "{{\"suite\":\"fs-bem\",\"case\":\"iterative-cbie-vs-dense\",\"inputs\":\"analytic {name}; SI; e^-iwt\",\"panels\":{},\"k\":{k},\"leaf_capacity\":{leaf_capacity},\"fmm_order\":7,\"iters\":{},\"rel_residual\":{},\"err\":{error},\"power_err\":{power},\"tol\":1e-5,\"verdict\":\"{}\",\"wall_s\":{},\"near_field_bytes\":{},\"fmm_stats\":[{},{}]}}",
            surface.areas().len(),
            report.iters,
            report.rel_residual,
            if error <= 1e-5 { "pass" } else { "fail" },
            receipt.wall_s,
            receipt.near_field_bytes,
            receipt.fmm_stats[0],
            receipt.fmm_stats[1],
        );
        assert!(report.converged_euclidean());
        assert!(error <= 1e-5, "{name} k={k} leaf={leaf_capacity}: {error}");
        assert!(power <= 1e-5, "{name} k={k} leaf={leaf_capacity}: {power}");
    }
}

/// G1: the FMM order is live, and the near/far split is not hiding a
/// centroid-rule floor: the difference from dense falls with the order.
#[test]
fn iterative_cbie_fmm_order_ladder_matches_dense() {
    let gate = CancelGate::new_clock_free();
    let medium = Medium::air();
    let surface = SpherePanels::icosphere(1.0, 4).unwrap();
    let (k, v) = (1.0, velocity(&surface));
    let dense = with_dense_workers(workers(), || {
        solve_radiation(&surface, k, medium, &v, Formulation::PlainCbie)
    })
    .unwrap();
    let mut errors = Vec::new();
    for p in [3, 5, 7] {
        let o = IterativeOptions {
            tol: 1e-10,
            fmm_order: p,
            leaf_capacity: 32,
            max_near_entries: 1 << 26,
            workers: workers(),
            ..IterativeOptions::default()
        };
        let (actual, report, _) =
            solve_radiation_iterative(&surface, k, medium, &v, Formulation::PlainCbie, &o, &gate)
                .unwrap();
        assert!(report.converged_euclidean());
        let error = relative(&actual.pressure, &dense.pressure);
        println!(
            "{{\"suite\":\"fs-bem\",\"case\":\"iterative-cbie-order-ladder\",\"inputs\":\"analytic sphere; SI\",\"axis\":\"fmm-order\",\"rung\":{p},\"panels\":5120,\"k\":{k},\"err\":{error},\"tol\":1e-5}}"
        );
        errors.push(error);
    }
    assert!(
        errors[0] > errors[1] && errors[1] > errors[2],
        "FMM order must be live: {errors:?}"
    );
    assert!(errors[2] <= 1e-5, "{errors:?}");
}

/// G5: declared workers and the M2L row cache change no bit of the near
/// table, the far passes or GMRES; a repeat is bitwise.
#[test]
fn iterative_cbie_workers_and_cache_are_bitwise() {
    let gate = CancelGate::new_clock_free();
    let medium = Medium::air();
    let surface = SpherePanels::icosphere(1.0, 4).unwrap();
    let v = velocity(&surface);
    let solve = |workers: usize, m2l_cache_bytes: usize| {
        let o = IterativeOptions {
            tol: 1e-8,
            fmm_order: 5,
            leaf_capacity: 32,
            max_near_entries: 1 << 26,
            workers,
            m2l_cache_bytes,
            ..IterativeOptions::default()
        };
        let (solution, report, _) =
            solve_radiation_iterative(&surface, 1.0, medium, &v, Formulation::PlainCbie, &o, &gate)
                .unwrap();
        (bits(&solution.pressure), report.iters)
    };
    let serial = solve(1, 0);
    assert_eq!(solve(1, 0), serial);
    assert_eq!(solve(5, 0), serial);
    assert_eq!(solve(16, 1 << 30), serial);
}

/// G0: the leaf-width admission is a named refusal at the measured boundary
/// of the seven-point far rule, never a silent loss of accuracy.
#[test]
fn iterative_cbie_refuses_leaves_narrower_than_the_far_rule() {
    let gate = CancelGate::new_clock_free();
    let surface = SpherePanels::icosphere(1.0, 3).unwrap();
    // 1280 panels: capacity 32 gives depth 2 (leaf 0.5 m, 5.7 panel radii);
    // capacity 16 gives depth 3 (0.25 m, 2.9 radii).
    let at = |leaf_capacity| {
        IterativeRadiation::prepare(
            &surface,
            1.0,
            Medium::air(),
            Formulation::PlainCbie,
            &IterativeOptions {
                leaf_capacity,
                ..IterativeOptions::default()
            },
            &gate,
        )
        .map(|_| ())
    };
    assert!(at(32).is_ok());
    let refused = at(16).unwrap_err().to_string();
    assert!(refused.contains("CBIE_MIN_LEAF_RADII"), "{refused}");
    assert!(CBIE_MIN_LEAF_RADII >= 4.0);
}
