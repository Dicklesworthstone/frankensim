//! G0/G1/G2/G3/G4/G5 qualification of the actual FMM Burton-Miller solve.
//! Large tests are intentional: run release, --test-threads=1, on ts1/trj.
use fs_bem::helmholtz::iterative::{
    IterativeError, IterativeOptions, IterativeRadiation, Preconditioner, ResourceReceipt,
    solve_radiation_iterative,
};
use fs_bem::helmholtz::{Formulation, MAX_DENSE_PANELS, Medium, solve_radiation};
use fs_bem::panel3d::{MAX_SURFACE_PANELS, SpherePanels};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};
use fs_solver::krylov::SolveReport;
use std::sync::Once;

fn uniform(n: usize) -> Vec<C64> {
    vec![C64::ONE; n]
}
fn velocity(surface: &SpherePanels) -> Vec<C64> {
    surface
        .normals()
        .iter()
        .enumerate()
        .map(|(i, n)| C64::new(1.0 + 0.2 * n[0] + 0.03 * det::sin(i as f64), 0.1 * n[1]))
        .collect()
}
fn opts() -> IterativeOptions {
    IterativeOptions {
        tol: 1e-6,
        ..IterativeOptions::default()
    }
}
fn oracle(ka: f64) -> C64 {
    let m = Medium::air();
    C64::new(ka * ka, -ka).scale(m.density * m.sound_speed / (1.0 + ka * ka))
}
fn mean(surface: &SpherePanels, p: &[C64]) -> C64 {
    p.iter()
        .zip(surface.areas())
        .fold(C64::ZERO, |s, (&p, &a)| s + p.scale(a))
        .scale(1.0 / surface.areas().iter().sum::<f64>())
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
fn quoted(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}
fn emit(
    case: &str,
    surface: &SpherePanels,
    ka: f64,
    options: &IterativeOptions,
    report: &SolveReport,
    receipt: &ResourceReceipt,
    error: f64,
    tolerance: f64,
    axis: &str,
    rung: usize,
) {
    static RUN: Once = Once::new();
    RUN.call_once(|| {
        let git=std::process::Command::new("git").args(["rev-parse","HEAD"]).output().unwrap();
        let dirty=std::process::Command::new("git").args(["status","--porcelain","--untracked-files=no"]).output().unwrap();
        assert!(git.status.success() && dirty.status.success());
        println!("{{\"suite\":\"fs-bem\",\"case\":\"iterative-run\",\"git_sha\":{},\"dirty\":{},\"host\":{},\"command\":{},\"profile\":{}}}",
            quoted(String::from_utf8(git.stdout).unwrap().trim()), !dirty.stdout.is_empty(), quoted(&receipt.host),
            quoted(&std::env::args().collect::<Vec<_>>().join(" ")),
            quoted(if cfg!(debug_assertions) {"debug"} else {"release"}));
    });
    println!(
        "{{\"suite\":\"fs-bem\",\"case\":{},\"inputs\":\"analytic sphere or explicitly authored lathed profile; SI; e^-iwt\",\"axis\":{},\"rung\":{},\"panels\":{},\"ka\":{},\"ppw\":{},\"precond\":{},\"iters\":{},\"rel_residual\":{},\"observed\":{},\"err\":{},\"tol\":{},\"verdict\":{},\"wall_s\":{},\"elapsed_s\":{},\"peak_rss_bytes\":{},\"matvecs\":{},\"fmm_order\":{},\"near_field_bytes\":{},\"fmm_stats\":[{},{}]}}",
        quoted(case),
        quoted(axis),
        rung,
        surface.areas().len(),
        ka,
        std::f64::consts::TAU / (ka * surface.areas().iter().copied().fold(0.0, f64::max).sqrt()),
        quoted(&format!("{:?}", options.precond)),
        report.iters,
        report.rel_residual,
        error,
        error,
        tolerance,
        quoted(if report.converged_euclidean() && error <= tolerance {
            "pass"
        } else {
            "fail"
        }),
        receipt.wall_s,
        receipt.wall_s,
        receipt
            .peak_rss_bytes
            .map_or_else(|| "null".into(), |n| n.to_string()),
        receipt.matvecs,
        receipt.fmm_order,
        receipt.near_field_bytes,
        receipt.fmm_stats[0],
        receipt.fmm_stats[1]
    );
}

#[test]
fn bem_iterative_matches_dense() {
    let gate = CancelGate::new_clock_free();
    // The top sphere rung has 5120 panels, immediately below the dense cap.
    // The larger rungs exercise M2L as well as exact adjacent-leaf P2P.
    for s in [1, 2, 3, 4] {
        let surface = SpherePanels::icosphere(1.0, s).unwrap();
        let v = velocity(&surface);
        let o = IterativeOptions {
            fmm_order: 7,
            leaf_capacity: 128,
            ..opts()
        };
        // At k=pi, s=1 resolves only 4.92 panels/wavelength (<6 required).
        // Use k=2.5 there (6.18 panels/wavelength); finer meshes retain pi
        // to exercise dense/iterative agreement at the fictitious frequency.
        let high_k = if s == 1 { 2.5 } else { std::f64::consts::PI };
        for k in [1.0, high_k] {
            let dense =
                solve_radiation(&surface, k, Medium::air(), &v, Formulation::BurtonMiller).unwrap();
            let (actual, report, receipt) = solve_radiation_iterative(
                &surface,
                k,
                Medium::air(),
                &v,
                Formulation::BurtonMiller,
                &o,
                &gate,
            )
            .unwrap();
            let error = relative(&actual.pressure, &dense.pressure);
            emit(
                "bem-iterative-dense-sphere",
                &surface,
                k,
                &o,
                &report,
                &receipt,
                error,
                10.0 * o.tol,
                "panels",
                s as usize,
            );
            assert!(error <= 10.0 * o.tol, "sphere s={s} k={k}: {error}");
        }
    }
    let bell =
        fs_bem::radiation_bake::lathe_profile(&[(0.0, 0.1), (0.4, 0.14), (0.8, 0.3)], 32).unwrap();
    let driven = bell.driven_mask();
    let surface = SpherePanels::from_triangles(bell.triangles).unwrap();
    let v: Vec<_> = driven
        .iter()
        .map(|&d| if d { C64::ONE } else { C64::ZERO })
        .collect();
    let o = IterativeOptions {
        fmm_order: 7,
        leaf_capacity: 128,
        ..opts()
    };
    let dense =
        solve_radiation(&surface, 1.0, Medium::air(), &v, Formulation::BurtonMiller).unwrap();
    let (actual, report, receipt) = solve_radiation_iterative(
        &surface,
        1.0,
        Medium::air(),
        &v,
        Formulation::BurtonMiller,
        &o,
        &gate,
    )
    .unwrap();
    let error = relative(&actual.pressure, &dense.pressure);
    emit(
        "bem-iterative-dense-lathed-bell",
        &surface,
        1.0,
        &o,
        &report,
        &receipt,
        error,
        10.0 * o.tol,
        "geometry",
        0,
    );
    assert!(error <= 10.0 * o.tol);
}

#[test]
fn bem_fmm_order_ladder_matches_dense() {
    let gate = CancelGate::new_clock_free();
    let surface = SpherePanels::icosphere(1.0, 2).unwrap();
    let v = velocity(&surface);
    let k = 2.5;
    let dense = solve_radiation(&surface, k, Medium::air(), &v, Formulation::BurtonMiller).unwrap();
    let mut errors = Vec::new();
    for p in [4, 6, 7] {
        let o = IterativeOptions {
            fmm_order: p,
            leaf_capacity: 16,
            ..opts()
        };
        let (actual, report, receipt) = solve_radiation_iterative(
            &surface,
            k,
            Medium::air(),
            &v,
            Formulation::BurtonMiller,
            &o,
            &gate,
        )
        .unwrap();
        let error = relative(&actual.pressure, &dense.pressure);
        emit(
            "bem-iterative-order-ladder",
            &surface,
            k,
            &o,
            &report,
            &receipt,
            error,
            10.0 * o.tol,
            "fmm-order",
            p,
        );
        errors.push(error);
    }
    assert!(errors[2] < errors[0], "FMM order must be live: {errors:?}");
    assert!(
        errors[2] <= 10.0 * opts().tol,
        "finest-order dense error: {errors:?}"
    );
}

#[test]
fn bem_64k_pulsating_sphere() {
    let surface = SpherePanels::icosphere(1.0, 6).unwrap();
    assert_eq!(surface.areas().len(), 81_920);
    let gate = CancelGate::new_clock_free();
    let k = 5.0 / (1.0 + 1e-9); // root-box padding stays inside W5a's kD<=10.
    let v = uniform(surface.areas().len());
    let reference = oracle(k);
    let mut errors = Vec::new();
    // All three explicit orders fit the owner's original work/storage caps.
    // Only the finest is the qualification; coarse rows are measurements.
    for p in [2, 3, 4] {
        let o = IterativeOptions {
            fmm_order: p,
            leaf_capacity: 128,
            ..opts()
        };
        let (solution, report, receipt) = solve_radiation_iterative(
            &surface,
            k,
            Medium::air(),
            &v,
            Formulation::BurtonMiller,
            &o,
            &gate,
        )
        .unwrap();
        let error = (mean(&surface, &solution.pressure) - reference).abs() / reference.abs();
        emit(
            "bem-64k-pulsating-sphere",
            &surface,
            k,
            &o,
            &report,
            &receipt,
            error,
            1e-3,
            "fmm-order",
            p,
        );
        assert!(report.converged_euclidean() && report.rel_residual <= 1e-6);
        assert!(receipt.matvecs > 0 && receipt.wall_s > 0.0 && receipt.near_field_bytes > 0);
        #[cfg(target_os = "linux")]
        assert!(receipt.peak_rss_bytes.is_some_and(|b| b > 0));
        errors.push(error);
    }
    assert!(
        errors[2] < errors[0],
        "large-rung order must be live: {errors:?}"
    );
    assert!(
        errors[2] <= 1e-3,
        "81920-panel analytic mean pressure: {errors:?}"
    );
}

#[test]
fn bem_preconditioner_iteration_counts() {
    let gate = CancelGate::new_clock_free();
    let mut counts = Vec::new();
    for (s, k) in [(1, 0.625), (2, 1.25), (3, 2.5), (4, 5.0)] {
        let surface = SpherePanels::icosphere(1.0, s).unwrap();
        let v = velocity(&surface);
        let mut row = Vec::new();
        for precond in [
            Preconditioner::None,
            Preconditioner::NearFieldBlockJacobi {
                max_block_panels: 32,
            },
        ] {
            let o = IterativeOptions { precond, ..opts() };
            let (_, report, receipt) = solve_radiation_iterative(
                &surface,
                k,
                Medium::air(),
                &v,
                Formulation::BurtonMiller,
                &o,
                &gate,
            )
            .unwrap();
            emit(
                "bem-preconditioner-iteration-counts",
                &surface,
                k,
                &o,
                &report,
                &receipt,
                report.rel_residual,
                o.tol,
                "panels-fixed-ppw",
                s as usize,
            );
            row.push(report.iters);
        }
        counts.push((surface.areas().len(), row[1]));
    }
    let (n0, i0) = counts[0];
    // Fixed before observing native data: one initial iteration count per
    // e-fold of N. Unpreconditioned counts remain the contrast, not a fit.
    for &(n, iters) in &counts {
        let bound = i0 as f64 * (1.0 + det::ln(n as f64 / n0 as f64));
        assert!(
            iters as f64 <= bound.ceil(),
            "iteration growth exceeds logarithmic envelope: {counts:?}"
        );
    }
}

#[test]
fn bem_gmres_resume_bitwise() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let gate = CancelGate::new_clock_free();
    let o = IterativeOptions {
        restart: 3,
        tol: 1e-9,
        max_iters: 600,
        ..opts()
    };
    let op = IterativeRadiation::prepare(
        &surface,
        2.5,
        Medium::air(),
        Formulation::BurtonMiller,
        &o,
        &gate,
    )
    .unwrap();
    let initial = op.start(&velocity(&surface), &gate).unwrap();
    let mut uninterrupted = initial.clone();
    op.advance(&mut uninterrupted, 600, &gate).unwrap();
    let mut first = initial;
    op.advance(&mut first, 1, &gate).unwrap();
    assert!(first.iterations() > 0 && !first.report(o.tol).converged_euclidean());
    let mut resumed = first.clone();
    op.advance(&mut resumed, 599, &gate).unwrap();
    let (a, ar, _) = op.finish(&uninterrupted).unwrap();
    let (b, br, _) = op.finish(&resumed).unwrap();
    assert_eq!(bits(&a.pressure), bits(&b.pressure));
    assert_eq!(ar.iters, br.iters);
    assert_eq!(
        ar.history.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        br.history.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
    );
}

#[test]
fn bem_near_field_cache_is_bitwise() {
    // The default caches both matvec operators' adjacent-leaf pairs; zero and
    // a budget below the cache disable it. All three solves are bit-identical.
    let surface = SpherePanels::icosphere(1.0, 2).unwrap();
    let gate = CancelGate::new_clock_free();
    let solve = |near_field_cache_bytes| {
        let o = IterativeOptions {
            leaf_capacity: 16,
            near_field_cache_bytes,
            ..opts()
        };
        solve_radiation_iterative(
            &surface,
            2.0,
            Medium::air(),
            &velocity(&surface),
            Formulation::BurtonMiller,
            &o,
            &gate,
        )
        .unwrap()
    };
    let (cached, report, receipt) = solve(IterativeOptions::default().near_field_cache_bytes);
    assert!(receipt.near_cache_bytes > 0);
    for budget in [0, receipt.near_cache_bytes / 2 - 1] {
        let (plain, plain_report, plain_receipt) = solve(budget);
        assert_eq!(plain_receipt.near_cache_bytes, 0);
        assert_eq!(bits(&cached.pressure), bits(&plain.pressure));
        assert_eq!(report.iters, plain_report.iters);
    }
}

#[test]
fn bem_iterative_refuses_unconverged() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let gate = CancelGate::new_clock_free();
    let o = IterativeOptions {
        max_iters: 1,
        restart: 12,
        tol: 1e-12,
        ..opts()
    };
    let error = solve_radiation_iterative(
        &surface,
        1.0,
        Medium::air(),
        &velocity(&surface),
        Formulation::BurtonMiller,
        &o,
        &gate,
    )
    .unwrap_err();
    let IterativeError::NotConverged {
        pressure,
        report,
        receipt,
    } = error
    else {
        panic!("{error:?}")
    };
    assert_eq!(report.iters, 1);
    assert!(!report.converged_euclidean());
    assert_eq!(pressure.len(), surface.areas().len());
    assert!(report.rel_residual > o.tol);
    println!(
        "{{\"suite\":\"fs-bem\",\"case\":\"bem-iterative-refuses-unconverged\",\"inputs\":\"analytic s1 sphere, nonuniform velocity\",\"axis\":\"max-iters\",\"rung\":1,\"observed\":{},\"err\":{},\"tol\":{},\"verdict\":\"refused\",\"elapsed_s\":{}}}",
        report.iters, report.rel_residual, o.tol, receipt.wall_s
    );
}

#[test]
fn bem_iterative_cancellation_bounded() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let gate = CancelGate::new_clock_free();
    let o = IterativeOptions {
        restart: 8,
        ..opts()
    };
    let op = IterativeRadiation::prepare(
        &surface,
        1.0,
        Medium::air(),
        Formulation::BurtonMiller,
        &o,
        &gate,
    )
    .unwrap();
    let mut state = op.start(&velocity(&surface), &gate).unwrap();
    let before = state.clone();
    let mut requested_at = None;
    let error = op
        .advance_controlled(&mut state, 20, &gate, &mut |index| {
            if index == 2 {
                requested_at = Some(index);
                gate.request();
            }
            gate.is_requested()
        })
        .unwrap_err();
    let IterativeError::Cancelled { matvecs } = error else {
        panic!("{error:?}")
    };
    assert_eq!(requested_at, Some(2));
    assert_eq!(matvecs, 2);
    assert_eq!(state.iterations(), before.iterations());
    assert!(matches!(
        op.finish(&state),
        Err(IterativeError::NotConverged { .. })
    ));
    // A new, uncancelled gate resumes the retained checkpoint exactly.
    let fresh = CancelGate::new_clock_free();
    op.advance(&mut state, 600, &fresh).unwrap();
    let mut reference = before;
    op.advance(&mut reference, 600, &fresh).unwrap();
    assert_eq!(
        bits(&op.finish(&state).unwrap().0.pressure),
        bits(&op.finish(&reference).unwrap().0.pressure)
    );
}

fn at_dense_cap(count: usize) -> SpherePanels {
    let surface = SpherePanels::icosphere(1.0, 4).unwrap();
    let extra = count - surface.areas().len();
    let mut triangles = Vec::with_capacity(count);
    for (i, &t) in surface.triangles().unwrap().iter().enumerate() {
        if i < extra / 2 {
            let c = std::array::from_fn(|a| (t[0][a] + t[1][a] + t[2][a]) / 3.0);
            triangles.extend([[t[0], t[1], c], [t[1], t[2], c], [t[2], t[0], c]]);
        } else if i == extra / 2 && extra % 2 == 1 {
            let m = std::array::from_fn(|a| f64::midpoint(t[0][a], t[1][a]));
            triangles.extend([[t[0], m, t[2]], [m, t[1], t[2]]]);
        } else {
            triangles.push(t);
        }
    }
    SpherePanels::from_triangles(triangles).unwrap()
}

#[test]
fn bem_iterative_caps() {
    let gate = CancelGate::new_clock_free();
    assert_eq!(MAX_DENSE_PANELS, 8192);
    assert_eq!(MAX_SURFACE_PANELS, 81920);
    let a = at_dense_cap(8192);
    let b = at_dense_cap(8193);
    assert!(
        solve_radiation(
            &a,
            1.0,
            Medium::air(),
            &vec![C64::ZERO; 8192],
            Formulation::BurtonMiller
        )
        .is_ok()
    );
    assert!(matches!(
        solve_radiation(
            &b,
            1.0,
            Medium::air(),
            &vec![C64::ZERO; 8193],
            Formulation::BurtonMiller
        ),
        Err(fs_bem::helmholtz::HelmholtzError::WorkCap { panels: 8193 })
    ));
    // Nonzero load proves the iterative route really solves above the cap.
    let o = opts();
    let (_, report, _) = solve_radiation_iterative(
        &b,
        1.0,
        Medium::air(),
        &uniform(8193),
        Formulation::BurtonMiller,
        &o,
        &gate,
    )
    .unwrap();
    assert!(report.converged_euclidean());
    let largest = SpherePanels::icosphere(1.0, 6).unwrap();
    let mut too_many = largest.triangles().unwrap().to_vec();
    too_many.push(too_many[0]);
    assert!(matches!(
        SpherePanels::from_triangles(too_many),
        Err(fs_bem::BemError::InvalidPanelCount { count: 81921, .. })
    ));

    let small = SpherePanels::icosphere(1.0, 1).unwrap();
    let v = velocity(&small);
    let (_, r, receipt) = solve_radiation_iterative(
        &small,
        1.0,
        Medium::air(),
        &v,
        Formulation::BurtonMiller,
        &o,
        &gate,
    )
    .unwrap();
    println!(
        "{{\"suite\":\"fs-bem\",\"case\":\"bem-iteration-cap-measured\",\"inputs\":\"analytic s1 sphere, nonuniform velocity\",\"observed\":{},\"oracle\":\"native unconstrained iteration count\",\"tol\":{},\"verdict\":\"pass\",\"elapsed_s\":{}}}",
        r.iters, o.tol, receipt.wall_s
    );
    assert!(r.iters > 0);
    let at = IterativeOptions {
        max_iters: r.iters,
        ..o.clone()
    };
    assert!(
        solve_radiation_iterative(
            &small,
            1.0,
            Medium::air(),
            &v,
            Formulation::BurtonMiller,
            &at,
            &gate
        )
        .is_ok()
    );
    let below = IterativeOptions {
        max_iters: r.iters - 1,
        ..o
    };
    assert!(matches!(
        solve_radiation_iterative(
            &small,
            1.0,
            Medium::air(),
            &v,
            Formulation::BurtonMiller,
            &below,
            &gate
        ),
        Err(IterativeError::NotConverged { .. })
    ));
}

#[test]
fn bem_iterative_input_and_near_storage_boundaries() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let gate = CancelGate::new_clock_free();
    let prepare = |o: &IterativeOptions| {
        IterativeRadiation::prepare(
            &surface,
            1.0,
            Medium::air(),
            Formulation::BurtonMiller,
            o,
            &gate,
        )
    };
    let mut o = opts();
    o.max_near_entries = 0;
    let required = match prepare(&o) {
        Err(IterativeError::NearFieldBudget {
            required_entries,
            maximum_entries: 0,
        }) => required_entries,
        _ => panic!("expected a near-field budget refusal"),
    };
    assert!(required > 0);
    o.max_near_entries = required - 1;
    assert!(matches!(
        prepare(&o),
        Err(IterativeError::NearFieldBudget { .. })
    ));
    o.max_near_entries = required;
    let op = prepare(&o).unwrap();
    assert!(matches!(
        op.start(&[], &gate),
        Err(IterativeError::Invalid { .. })
    ));
    let mut bad = uniform(surface.areas().len());
    bad[0] = C64::new(f64::NAN, 0.0);
    assert!(matches!(
        op.start(&bad, &gate),
        Err(IterativeError::Invalid { .. })
    ));
    assert!(matches!(
        op.apply(&[], &gate),
        Err(IterativeError::Invalid { .. })
    ));
    let mut zero = op
        .start(&vec![C64::ZERO; surface.areas().len()], &gate)
        .unwrap();
    assert!(
        op.advance(&mut zero, 1, &gate)
            .unwrap()
            .converged_euclidean()
    );
    assert!(
        op.finish(&zero)
            .unwrap()
            .0
            .pressure
            .iter()
            .all(|p| *p == C64::ZERO)
    );
    let mut wrong = opts();
    wrong.restart += 1;
    let other = prepare(&wrong).unwrap();
    assert!(matches!(
        other.advance(&mut zero, 1, &gate),
        Err(IterativeError::Invalid { .. })
    ));
    assert!(matches!(
        other.finish(&zero),
        Err(IterativeError::Invalid { .. })
    ));
    for block in [0, 129] {
        wrong.precond = Preconditioner::NearFieldBlockJacobi {
            max_block_panels: block,
        };
        assert!(matches!(
            prepare(&wrong),
            Err(IterativeError::Invalid { .. })
        ));
    }
    for restart in [0, 257] {
        wrong = opts();
        wrong.restart = restart;
        assert!(matches!(
            prepare(&wrong),
            Err(IterativeError::Invalid { .. })
        ));
    }
    assert!(matches!(
        IterativeRadiation::prepare(
            &surface,
            1.0,
            Medium::air(),
            Formulation::PlainCbie,
            &opts(),
            &gate
        ),
        Err(IterativeError::Invalid { .. })
    ));
    let fine = SpherePanels::icosphere(1.0, 3).unwrap();
    assert!(matches!(
        IterativeRadiation::prepare(
            &fine,
            6.0,
            Medium::air(),
            Formulation::BurtonMiller,
            &opts(),
            &gate
        ),
        Err(IterativeError::Fmm(
            fs_fmm::FmmError::HelmholtzRegime { .. }
        ))
    ));
    gate.request();
    assert!(matches!(
        prepare(&opts()),
        Err(IterativeError::Cancelled { matvecs: 0 })
    ));
}

#[test]
fn bem_directional_matches_dense_above_w5a_regime() {
    let profile: Vec<_> = (0..10)
        .map(|i| {
            let t = i as f64 / 9.0;
            (2.0 * t, 0.02 + 0.02 * t * t)
        })
        .collect();
    let bell = fs_bem::radiation_bake::lathe_profile(&profile, 32).unwrap();
    let mask = bell.driven_mask();
    let surface = SpherePanels::from_triangles(bell.triangles).unwrap();
    let v: Vec<_> = mask
        .iter()
        .map(|&m| if m { C64::ONE } else { C64::ZERO })
        .collect();
    let k = 8.0;
    let gate = CancelGate::new_clock_free();
    assert!(matches!(
        solve_radiation_iterative(
            &surface,
            k,
            Medium::air(),
            &v,
            Formulation::BurtonMiller,
            &opts(),
            &gate
        ),
        Err(IterativeError::Fmm(
            fs_fmm::FmmError::HelmholtzRegime { .. }
        ))
    ));
    let dense = solve_radiation(&surface, k, Medium::air(), &v, Formulation::BurtonMiller).unwrap();
    let mut errors = Vec::new();
    for p in [5, 6, 7] {
        let mut directional = fs_fmm::DirectionalOptions::for_tolerance(1e-3).unwrap();
        directional.order = p;
        directional.low_freq_threshold = 0.5;
        directional.leaf_capacity = 128;
        let o = IterativeOptions {
            directional: Some(directional),
            ..opts()
        };
        let (actual, report, receipt) = solve_radiation_iterative(
            &surface,
            k,
            Medium::air(),
            &v,
            Formulation::BurtonMiller,
            &o,
            &gate,
        )
        .unwrap();
        let error = relative(&actual.pressure, &dense.pressure);
        emit(
            "bem-directional-dense",
            &surface,
            k,
            &o,
            &report,
            &receipt,
            error,
            10.0 * o.tol,
            "directional-order",
            p,
        );
        assert_eq!(receipt.fmm_order, p);
        assert!(
            receipt
                .fmm_stats
                .iter()
                .all(|s| s.contains("\"directional\":true"))
        );
        assert!(report.converged_euclidean());
        errors.push(error);
    }
    assert!(
        errors[2] < errors[0] && errors[2] <= 10.0 * opts().tol,
        "directional BM dense error: {errors:?}"
    );
}
