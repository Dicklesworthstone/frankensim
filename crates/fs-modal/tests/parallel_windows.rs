//! IRP .67: the serial slice owner, actual executor workers, real 3-D bodies.
use fs_exec::{CancelGate, Cx};
use fs_math::det;
use fs_modal::{ModalError, SliceOptions, WindowPlan, slice_window, slice_windows};
use fs_sparse::{Coo, Csr, DirectOrdering, LdltOptions, NdOptions};
use std::time::Instant;

const MIB: usize = 1024 * 1024;
fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            fs_exec::StreamKey {
                seed: 67,
                kernel_id: 67,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}
fn chain(n: usize) -> Csr {
    let mut c = Coo::new(n, n);
    for i in 0..n {
        c.push(i, i, 2.0);
        if i > 0 {
            c.push(i, i - 1, -1.0);
            c.push(i - 1, i, -1.0);
        }
    }
    c.assemble()
}
fn eigen(n: usize, i: usize) -> f64 {
    2.0 - 2.0 * det::cos(std::f64::consts::PI * i as f64 / (n + 1) as f64)
}
fn windows(n: usize) -> Vec<(f64, f64)> {
    let cuts: Vec<_> = (0..=8)
        .map(|i| {
            if i == 0 {
                -0.1
            } else {
                (eigen(n, 2 * i) + eigen(n, 2 * i + 1)) * 0.5
            }
        })
        .collect();
    cuts.windows(2).map(|w| (w[0], w[1])).collect()
}
fn options() -> SliceOptions {
    SliceOptions {
        ritz_tol: 1e-12,
        ldlt: LdltOptions {
            max_factor_bytes: Some(64 * MIB),
            ..LdltOptions::default()
        },
        ..SliceOptions::default()
    }
}
fn plan(n: usize, threads: usize) -> WindowPlan {
    WindowPlan {
        windows: windows(n),
        threads,
        maximum_windows: 8,
        maximum_bytes: 2048 * MIB,
        maximum_relative_halfwidth: 1e-8,
        coordinates: None,
    }
}
fn bits(r: &fs_modal::MultiSliceReport) -> Vec<u64> {
    let mut v = Vec::new();
    for r in &r.reports {
        v.extend([
            r.window.0.to_bits(),
            r.window.1.to_bits(),
            r.below_low as u64,
            r.below_high as u64,
            r.expected as u64,
            r.stats.shift.to_bits(),
            r.stats.lanczos_iters as u64,
            r.stats.restarts as u64,
        ]);
        for m in &r.modes {
            v.extend([
                m.lambda.to_bits(),
                m.residual.to_bits(),
                m.interval.0.to_bits(),
                m.interval.1.to_bits(),
            ]);
            v.extend(m.phi.iter().map(|x| x.to_bits()));
        }
    }
    v
}
#[test]
fn g1_parallel_windows_certified_and_thread_invariant() {
    let n = 80;
    let k = chain(n);
    let m = Csr::identity(n);
    let gate = CancelGate::new_clock_free();
    for ordering in [
        DirectOrdering::Amd,
        DirectOrdering::NestedDissection(NdOptions {
            leaf_size: 16,
            refinement_passes: 4,
        }),
    ] {
        let opts = SliceOptions {
            ordering,
            ..options()
        };
        let first = with_cx(&gate, |cx| slice_windows(&k, &m, &plan(n, 1), &opts, cx)).unwrap();
        let count: usize = first.reports.iter().map(|r| r.expected).sum();
        assert_eq!(count, 16);
        assert_eq!(
            count,
            first.reports.last().unwrap().below_high - first.reports[0].below_low
        );
        for pair in first.reports.windows(2) {
            assert_eq!(pair[0].below_high, pair[1].below_low);
        }
        for report in &first.reports {
            let serial = slice_window(&k, &m, report.window, &opts).unwrap();
            assert_eq!(report.expected, serial.expected);
            for (a, b) in report.modes.iter().zip(&serial.modes) {
                assert_eq!(a.lambda.to_bits(), b.lambda.to_bits());
                assert_eq!(a.residual.to_bits(), b.residual.to_bits());
                assert_eq!(a.interval.0.to_bits(), b.interval.0.to_bits());
                assert_eq!(a.interval.1.to_bits(), b.interval.1.to_bits());
                assert_eq!(
                    a.phi.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                    b.phi.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
                );
                assert!(a.residual / a.lambda.abs() <= 1e-8);
            }
            assert_eq!(report.modes.len(), report.expected);
            let truth: Vec<_> = (1..=n)
                .map(|i| eigen(n, i))
                .filter(|&x| x > report.window.0 && x <= report.window.1)
                .collect();
            assert_eq!(truth.len(), report.expected);
            for (&truth, m) in truth.iter().zip(&report.modes) {
                assert!((truth - m.lambda).abs() <= m.residual + 32.0 * f64::EPSILON);
            }
        }
        for threads in [2, 7] {
            let r = with_cx(&gate, |cx| {
                slice_windows(&k, &m, &plan(n, threads), &opts, cx)
            })
            .unwrap();
            assert_eq!(r.receipt.threads, threads);
            assert_eq!(bits(&first), bits(&r));
        }
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.5-W1\",\"test\":\"g1_parallel_windows_certified_and_thread_invariant\",\"inputs\":{{\"fixture\":\"chain80\",\"n\":80,\"ordering\":\"{ordering:?}\",\"threads\":[1,2,7]}},\"observed\":{{\"counts\":{:?},\"nnz_l\":{},\"peak_bytes\":{},\"max_rel_halfwidth\":{}}},\"receipt\":{{\"wall_s_1_thread\":{}}},\"tolerance\":{{\"relative_halfwidth_max\":1e-8,\"thread_result_bits\":\"exact\"}},\"verdict\":\"pass\"}}",
            first.reports.iter().map(|r| r.expected).collect::<Vec<_>>(),
            first.symbolic_nnz_l,
            first.receipt.admitted_peak_bytes,
            first.receipt.maximum_relative_halfwidth,
            first.receipt.wall_s
        );
    }
}
#[test]
fn g1_one_symbolic_analysis_serves_all_shifts() {
    let n = 80;
    let k = chain(n);
    let m = Csr::identity(n);
    let gate = CancelGate::new_clock_free();
    let r = with_cx(&gate, |cx| {
        slice_windows(&k, &m, &plan(n, 7), &options(), cx)
    })
    .unwrap();
    assert_eq!(r.receipt.pencil_symbolic_analyses, 1);
    assert_eq!(r.receipt.mass_symbolic_analyses, 1);
    assert_eq!(
        r.receipt.endpoint_factorizations, 9,
        "shared endpoints factor exactly once"
    );
    assert_eq!(
        r.receipt.interior_factorizations,
        r.reports
            .iter()
            .map(|r| r.stats.factorizations)
            .sum::<usize>()
    );
    assert!(r.receipt.interior_factorizations >= 8);
    assert!(r.symbolic_nnz_l >= n);
}
#[test]
fn g0_window_count_cap_and_cap_plus_one() {
    let n = 80;
    let k = chain(n);
    let m = Csr::identity(n);
    let gate = CancelGate::new_clock_free();
    let opts = options();
    let mut p = plan(n, 1);
    let r = with_cx(&gate, |cx| slice_windows(&k, &m, &p, &opts, cx)).unwrap();
    p.maximum_bytes = r.receipt.admitted_peak_bytes;
    with_cx(&gate, |cx| slice_windows(&k, &m, &p, &opts, cx)).unwrap();
    p.maximum_bytes -= 1;
    assert!(matches!(
        with_cx(&gate, |cx| slice_windows(&k, &m, &p, &opts, cx)),
        Err(ModalError::MemoryBudget { .. })
    ));
    p = plan(n, 1);
    p.maximum_windows = 7;
    assert!(matches!(
        with_cx(&gate, |cx| slice_windows(&k, &m, &p, &opts, cx)),
        Err(ModalError::WindowBudget {
            requested: 8,
            cap: 7
        })
    ));
    p = plan(n, 1);
    p.windows[1].0 = p.windows[0].0;
    assert!(matches!(
        with_cx(&gate, |cx| slice_windows(&k, &m, &p, &opts, cx)),
        Err(ModalError::WindowPlan { .. })
    ));
    p = plan(n, 0);
    assert!(with_cx(&gate, |cx| slice_windows(&k, &m, &p, &opts, cx)).is_err());
    assert!(
        matches!(
            with_cx(&gate, |cx| slice_windows(
                &k,
                &m,
                &plan(n, 1),
                &SliceOptions::default(),
                cx
            )),
            Err(ModalError::WindowPlan { .. })
        ),
        "parallel work needs an explicit factor cap"
    );
}
#[test]
fn g4_cancellation_drains_and_returns_cancelled() {
    let n = 1000;
    let k = chain(n);
    let m = Csr::identity(n);
    let gate = CancelGate::new();
    // Cancel a real scoped call. There is no partial-report success path;
    // returning from the executor has joined its workers before this gate
    // and all borrowed matrices can leave scope. Then a fresh call succeeds.
    std::thread::scope(|scope| {
        scope.spawn(|| {
            std::thread::sleep(std::time::Duration::from_millis(2));
            gate.request();
        });
        assert!(matches!(
            with_cx(&gate, |cx| slice_windows(
                &k,
                &m,
                &plan(n, 7),
                &options(),
                cx
            )),
            Err(ModalError::Cancelled)
        ));
    });
    let fresh = CancelGate::new_clock_free();
    let r = with_cx(&fresh, |cx| {
        slice_windows(&chain(80), &Csr::identity(80), &plan(80, 2), &options(), cx)
    })
    .unwrap();
    assert_eq!(r.reports.len(), 8);
}

#[cfg(not(debug_assertions))]
fn laplacian3(s: usize) -> (Csr, Vec<[f64; 3]>) {
    let n = s * s * s;
    let mut rp = Vec::with_capacity(n + 1);
    let mut ci = Vec::with_capacity(7 * n);
    let mut va = Vec::with_capacity(7 * n);
    let mut xyz = Vec::with_capacity(n);
    rp.push(0);
    for z in 0..s {
        for y in 0..s {
            for x in 0..s {
                let v = (z * s + y) * s + x;
                let mut row = vec![(v, 6.0)];
                xyz.push([x as f64, y as f64, z as f64]);
                for (i, stride) in [(x, 1), (y, s), (z, s * s)] {
                    if i > 0 {
                        row.push((v - stride, -1.0));
                    }
                    if i + 1 < s {
                        row.push((v + stride, -1.0));
                    }
                }
                row.sort_by_key(|v| v.0);
                for (c, a) in row {
                    ci.push(c);
                    va.push(a);
                }
                rp.push(ci.len());
            }
        }
    }
    (Csr::from_parts(n, n, rp, ci, va), xyz)
}
#[cfg(not(debug_assertions))]
fn elastic_block(s: usize, quadratic: bool, cx: &Cx<'_>) -> (Csr, Csr, Vec<[f64; 3]>) {
    use fs_solid::linear3::{
        TetAssemblyBudget, TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
    };
    let index = |x, y, z| (z * (s + 1) + y) * (s + 1) + x;
    let mut nodes = Vec::new();
    let mut tets = Vec::new();
    let mut fixed = Vec::new();
    for z in 0..=s {
        for y in 0..=s {
            for x in 0..=s {
                let v = nodes.len();
                nodes.push([
                    x as f64 / s as f64,
                    0.5 * y as f64 / s as f64,
                    0.25 * z as f64 / s as f64,
                ]);
                if x == 0 {
                    fixed.extend([3 * v, 3 * v + 1, 3 * v + 2]);
                }
            }
        }
    }
    for z in 0..s {
        for y in 0..s {
            for x in 0..s {
                let v = [
                    index(x, y, z),
                    index(x + 1, y, z),
                    index(x, y + 1, z),
                    index(x + 1, y + 1, z),
                    index(x, y, z + 1),
                    index(x + 1, y, z + 1),
                    index(x, y + 1, z + 1),
                    index(x + 1, y + 1, z + 1),
                ];
                tets.extend([
                    [v[0], v[1], v[3], v[7]],
                    [v[0], v[3], v[2], v[7]],
                    [v[0], v[2], v[6], v[7]],
                    [v[0], v[6], v[4], v[7]],
                    [v[0], v[4], v[5], v[7]],
                    [v[0], v[5], v[1], v[7]],
                ]);
            }
        }
    }
    // Declared synthetic SI homogeneous elastic block; not a real specimen or
    // an instrument/material qualification. Consistent mass comes from the
    // existing 3-D P1 elasticity owner, never three decoupled scalar grids.
    let material = TetElasticMaterial::try_new(1000.0, 1e6, 0.3, [67; 32]).unwrap();
    let a = if quadratic {
        use fs_solid::linear3::tet10::{QuadraticTetProblem, TetQuadrature, promote_p1_to_p2};
        let (promoted, elements) = promote_p1_to_p2(&nodes, &tets, None, cx).unwrap();
        nodes = promoted;
        fixed = nodes
            .iter()
            .enumerate()
            .filter(|(_, p)| p[0] == 0.0)
            .flat_map(|(i, _)| [3 * i, 3 * i + 1, 3 * i + 2])
            .collect();
        QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget::standard(),
            quadrature: TetQuadrature {
                stiffness_points_per_dir: 3,
                mass_points_per_dir: 4,
            },
        }
        .assemble(cx)
        .unwrap()
    } else {
        TetLinearElasticProblem {
            nodes_m: &nodes,
            tetrahedra: &tets,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget::standard(),
        }
        .assemble(cx)
        .unwrap()
    };
    let xyz = a.free_dofs.iter().map(|&d| nodes[d / 3]).collect();
    (a.stiffness, a.mass, xyz)
}
/// Expensive ts1-only release lane. It contains the required full fixtures,
/// has no #[ignore], and never substitutes a smaller debug-scale fixture.
#[test]
#[cfg(not(debug_assertions))]
fn g2_million_dof_window_timing_receipt() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for fixture in ["elastic-block-1e5", "grid3d-1e6"] {
            let start = Instant::now();
            let (k, m, coordinates, window) = if fixture == "elastic-block-1e5" {
                let (k, m, x) = elastic_block(32, false, cx);
                assert!(k.nrows() >= 100_000);
                (k, m, x, (-1.0, 180.0))
            } else {
                let (k, x) = laplacian3(100);
                let m = Csr::identity(k.nrows());
                let low = 3.0 * eigen(100, 1);
                let next = 2.0 * eigen(100, 1) + eigen(100, 2);
                (k, m, x, (-1e-5, 0.5 * (low + next)))
            };
            let p = WindowPlan {
                windows: vec![window],
                threads: 1,
                maximum_windows: 1,
                maximum_bytes: 64 * 1024 * MIB,
                maximum_relative_halfwidth: 1e-8,
                coordinates: Some(coordinates),
            };
            let opts = SliceOptions {
                ordering: DirectOrdering::NestedDissection(NdOptions {
                    leaf_size: 64,
                    refinement_passes: 4,
                }),
                ldlt: LdltOptions {
                    max_factor_bytes: Some(24 * 1024 * MIB),
                    ..LdltOptions::default()
                },
                ritz_tol: 1e-13,
                max_lanczos: 240,
                max_restarts: 8,
                mass_diagonal_equilibration: fixture == "elastic-block-1e5",
            };
            let report = slice_windows(&k, &m, &p, &opts, cx).unwrap();
            let count = report.reports[0].expected;
            assert!(count > 0);
            assert_eq!(report.reports[0].modes.len(), count);
            if fixture == "grid3d-1e6" {
                assert_eq!(count, 1);
            }
            assert!(report.receipt.maximum_relative_halfwidth <= 1e-8);
            assert!(report.receipt.admitted_peak_bytes <= p.maximum_bytes);
            println!(
                "{{\"suite\":\"irp-s1\",\"bead\":\"M1.5-W1\",\"test\":\"g2_million_dof_window_timing_receipt\",\"inputs\":{{\"fixture\":\"{fixture}\",\"n\":{},\"window\":[{},{}],\"ordering\":\"geometric-nd\",\"threads\":1,\"memory_budget_bytes\":{}}},\"observed\":{{\"nnz_l\":{},\"peak_bytes\":{},\"peak_kind\":\"admitted-numeric-payload-not-RSS\",\"counts\":[{count}],\"max_rel_halfwidth\":{},\"factorizations\":{}}},\"receipt\":{{\"host_fingerprint\":\"{}-{}-{}\",\"wall_s\":{},\"pencil_symbolic_analyses\":{},\"endpoint_factorizations\":{}}},\"tolerance\":{{\"relative_halfwidth_max\":1e-8}},\"verdict\":\"pass\",\"elapsed_s\":{}}}",
                k.nrows(),
                window.0,
                window.1,
                p.maximum_bytes,
                report.symbolic_nnz_l,
                report.receipt.admitted_peak_bytes,
                report.receipt.maximum_relative_halfwidth,
                1 + report.receipt.endpoint_factorizations + report.receipt.interior_factorizations,
                std::env::var("HOSTNAME").unwrap_or_else(|_| "unknown".into()),
                std::env::consts::ARCH,
                std::env::consts::OS,
                report.receipt.wall_s,
                report.receipt.pencil_symbolic_analyses,
                report.receipt.endpoint_factorizations,
                start.elapsed().as_secs_f64()
            );
        }
    });
}

/// Elastic P1/TET10 ordering results are reports, never substitutes for the
/// gated seven-point Laplacian comparison or the 100k coupled elastic solve.
#[test]
#[cfg(not(debug_assertions))]
fn report_elastic_p1_tet10_fill() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for quadratic in [false, true] {
            let (k, _, coordinates) = elastic_block(4, quadratic, cx);
            for ordering in [
                DirectOrdering::Amd,
                DirectOrdering::NestedDissection(NdOptions {
                    leaf_size: 64,
                    refinement_passes: 4,
                }),
            ] {
                let start = Instant::now();
                let s = fs_sparse::SymbolicLdlt::analyze_with_coordinates(
                    &k,
                    ordering,
                    Some(&coordinates),
                )
                .unwrap();
                println!(
                    "{{\"suite\":\"irp-s1\",\"bead\":\"M1.5-W1\",\"test\":\"report_elastic_p1_tet10_fill\",\"inputs\":{{\"element\":\"{}\",\"n\":{},\"ordering\":\"{ordering:?}\"}},\"observed\":{{\"nnz_l\":{},\"estimated_flops\":{},\"estimated_numeric_payload_bytes\":{}}},\"verdict\":\"reported\",\"elapsed_s\":{}}}",
                    if quadratic { "TET10" } else { "P1" },
                    k.nrows(),
                    s.nnz_l(),
                    s.estimated_flops(),
                    s.estimated_factor_bytes(),
                    start.elapsed().as_secs_f64()
                );
            }
        }
    });
}
