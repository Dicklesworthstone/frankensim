//! fs-cutfem conformance battery (bead tfz.8).
//!
//! - cut-001 G0: certified cut classification — zero misclassification
//!   on adversarial tangent/near-tangent interfaces, enclosure
//!   containment law.
//! - cut-002 G0: cut-quadrature exactness on polynomial fixtures
//!   (linear interface: quadratic moments exact) and depth convergence
//!   on curved interfaces (error control documented).
//! - cut-003 G0: ghost-penalty conditioning independence of cut
//!   fraction (eigenvalue-verified curves; blowup without ghost).
//! - cut-004 G1: MMS optimal orders across RANDOMIZED cut
//!   configurations, slivers deliberately included.
//! - cut-005 G1: embedded-BC accuracy vs a body-fitted reference on a
//!   shared fixture.
//! - cut-006 G1: moving-interface sequence without re-meshing (the
//!   topology-optimization rehearsal).
//! - cut-007: aggregated-element fallback — conditioning restored and
//!   accuracy held with ghost penalty OFF (policy logged).
//! - cut-008: hanging-node patch test (linear reproduction on a graded
//!   tree) and estimate-driven h-refinement efficiency.

use fs_cutfem::{
    AggPolicy, CellClass, Circle, CutSdf, FemParams, HalfPlane, Quadtree, Space,
    condition_estimate, cut_cell_rules,
};
use fs_solver::krylov::CgState;
use fs_solver::op::CsrOp;
use fs_sparse::Coo;
use fs_sparse::precond::IdentityPrecond;
use std::f64::consts::PI;
use std::fmt::Write as _;

fn verdict(name: &str, pass: bool, details: &str) {
    println!(
        "{{\"test\":\"{name}\",\"verdict\":\"{}\",{details}}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "{name} failed: {details}");
}

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }
    #[allow(clippy::cast_precision_loss)]
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

// ---------------------------------------------------------------- MMS fixtures

fn sin_u(x: f64, y: f64) -> f64 {
    (PI * x).sin() * (PI * y).sin()
}
fn sin_f(x: f64, y: f64) -> f64 {
    2.0 * PI * PI * sin_u(x, y)
}
fn sin_grad(x: f64, y: f64) -> [f64; 2] {
    [
        PI * (PI * x).cos() * (PI * y).sin(),
        PI * (PI * x).sin() * (PI * y).cos(),
    ]
}

// ------------------------------------------------------------------- cut-001

#[test]
#[allow(clippy::too_many_lines)] // six fixtures declared inline: the adversarial list IS the test
fn cut_001_certified_classification() {
    let h16 = 1.0 / 16.0;
    let fixtures: Vec<(&str, Box<dyn CutSdf>)> = vec![
        (
            "circle-generic",
            Box::new(Circle {
                center: [0.5, 0.5],
                radius: 0.25,
            }),
        ),
        (
            "halfplane-on-gridline",
            Box::new(HalfPlane {
                normal: [1.0, 0.0],
                offset: 0.5,
            }),
        ),
        (
            "halfplane-through-nodes",
            Box::new(HalfPlane {
                normal: [
                    std::f64::consts::FRAC_1_SQRT_2,
                    std::f64::consts::FRAC_1_SQRT_2,
                ],
                offset: std::f64::consts::FRAC_1_SQRT_2,
            }),
        ),
        (
            "circle-tangent-to-gridline",
            Box::new(Circle {
                center: [0.5, 5.0 * h16],
                radius: 4.0 * h16,
            }),
        ),
        (
            "circle-near-tangent-1e13",
            Box::new(Circle {
                center: [0.5, 5.0 * h16],
                radius: 4.0 * h16 + 1e-13,
            }),
        ),
        (
            "tiny-circle-inside-one-cell",
            Box::new(Circle {
                center: [0.5 + 0.3 * h16, 0.5 + 0.4 * h16],
                radius: 0.05 * h16,
            }),
        ),
    ];
    let grid = Quadtree::uniform(4);
    let mut total_violations = 0usize;
    let mut total_conservative = 0usize;
    let mut rows = String::new();
    for (name, sdf) in &fixtures {
        let mut violations = 0usize;
        let mut conservative = 0usize;
        for c in grid.leaves() {
            let (lo, hi) = grid.rect(c);
            let iv = sdf.enclose(lo, hi);
            let class = if iv.hi() < 0.0 {
                CellClass::Inside
            } else if iv.lo() > 0.0 {
                CellClass::Outside
            } else {
                CellClass::Cut
            };
            let mut any_pos = false;
            let mut any_nonpos = false;
            for si in 0..=16 {
                for sj in 0..=16 {
                    let p = [
                        lo[0] + (hi[0] - lo[0]) * f64::from(si) / 16.0,
                        lo[1] + (hi[1] - lo[1]) * f64::from(sj) / 16.0,
                    ];
                    let v = sdf.value(p);
                    if !iv.contains(v) {
                        violations += 1; // containment law broken
                    }
                    if v > 0.0 {
                        any_pos = true;
                    } else {
                        any_nonpos = true;
                    }
                }
            }
            match class {
                CellClass::Inside if any_pos => violations += 1,
                CellClass::Outside if any_nonpos => violations += 1,
                CellClass::Cut if !(any_pos && any_nonpos) => conservative += 1,
                _ => {}
            }
        }
        let _ = write!(
            rows,
            "\"{name}\":{{\"violations\":{violations},\"conservative_cuts\":{conservative}}},"
        );
        total_violations += violations;
        total_conservative += conservative;
    }
    verdict(
        "cut-001",
        total_violations == 0,
        &format!(
            "\"detail\":\"certified classification, 6 adversarial fixtures\",\
             {rows}\"total_violations\":{total_violations},\
             \"total_conservative_cuts\":{total_conservative}"
        ),
    );
}

// ------------------------------------------------------------------- cut-002

/// Exact ∫ over a polygon of quadratic monomials via the degree-2
/// midpoint triangle rule on an ANALYTICALLY constructed polygon.
fn polygon_moments_exact(poly: &[[f64; 2]]) -> [f64; 6] {
    let mono =
        |p: [f64; 2]| -> [f64; 6] { [1.0, p[0], p[1], p[0] * p[0], p[0] * p[1], p[1] * p[1]] };
    let mut acc = [0.0f64; 6];
    for k in 1..poly.len() - 1 {
        let (p, q, r) = (poly[0], poly[k], poly[k + 1]);
        let area = 0.5 * ((q[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (q[1] - p[1]));
        for m in [
            [f64::midpoint(p[0], q[0]), f64::midpoint(p[1], q[1])],
            [f64::midpoint(q[0], r[0]), f64::midpoint(q[1], r[1])],
            [f64::midpoint(r[0], p[0]), f64::midpoint(r[1], p[1])],
        ] {
            let vals = mono(m);
            for (a, v) in acc.iter_mut().zip(vals) {
                *a += area / 3.0 * v;
            }
        }
    }
    acc
}

#[test]
fn cut_002_quadrature_exactness_and_depth_control() {
    // (a) Linear interface: quadratic moments EXACT (crossings by
    // bisection, polygon degree-2 rule).
    let hp = HalfPlane {
        normal: [0.6, 0.8],
        offset: 0.37,
    };
    // Analytic inside polygon on [0,1]²: walk corners + exact crossings.
    let corners = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let mut poly: Vec<[f64; 2]> = Vec::new();
    for e in 0..4 {
        let a = corners[e];
        let b = corners[(e + 1) % 4];
        let (fa, fb) = (hp.value(a), hp.value(b));
        if fa <= 0.0 {
            poly.push(a);
        }
        if (fa <= 0.0) != (fb <= 0.0) {
            let t = fa / (fa - fb); // exact for linear φ
            poly.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
        }
    }
    let exact = polygon_moments_exact(&poly);
    let mut max_rel = 0.0f64;
    for depth in [0u32, 3u32] {
        let rules = cut_cell_rules(&hp, [0.0, 0.0], [1.0, 1.0], depth);
        let mono =
            |p: [f64; 2]| -> [f64; 6] { [1.0, p[0], p[1], p[0] * p[0], p[0] * p[1], p[1] * p[1]] };
        let mut got = [0.0f64; 6];
        for &(p, w) in &rules.bulk {
            let vals = mono(p);
            for (g, v) in got.iter_mut().zip(vals) {
                *g += w * v;
            }
        }
        for (g, e) in got.iter().zip(exact) {
            max_rel = max_rel.max((g - e).abs() / e.abs().max(1e-30));
        }
    }
    // (b) Curved interface: area/perimeter converge quadratically in
    // depth (the error-control documentation).
    let circ = Circle {
        center: [0.5, 0.5],
        radius: 0.3,
    };
    let true_area = PI * 0.3 * 0.3;
    let true_len = 2.0 * PI * 0.3;
    let mut area_err = Vec::new();
    let mut len_err = Vec::new();
    for depth in [1u32, 3, 5] {
        let rules = cut_cell_rules(&circ, [0.0, 0.0], [1.0, 1.0], depth);
        let area: f64 = rules.bulk.iter().map(|&(_, w)| w).sum();
        let len: f64 = rules.iface.iter().map(|&(_, w, _)| w).sum();
        area_err.push((area - true_area).abs() / true_area);
        len_err.push((len - true_len).abs() / true_len);
    }
    let area_order = (area_err[1] / area_err[2]).log2() / 2.0;
    let len_order = (len_err[1] / len_err[2]).log2() / 2.0;
    // Depth-5 magnitudes match the chord-sagitta model (area ≈ h²/6r²,
    // length ≈ (h/r)²/24 relative): the control knob is the measured
    // second-order convergence in depth.
    let pass = max_rel < 1e-9
        && area_err[2] < 3e-3
        && len_err[2] < 1e-3
        && area_order > 1.7
        && len_order > 1.7;
    verdict(
        "cut-002",
        pass,
        &format!(
            "\"detail\":\"linear-interface quadratic moments exact; curved depth control\",\
             \"linear_max_rel\":{max_rel:.3e},\
             \"area_err\":[{:.3e},{:.3e},{:.3e}],\"area_order\":{area_order:.2},\
             \"len_err\":[{:.3e},{:.3e},{:.3e}],\"len_order\":{len_order:.2}",
            area_err[0], area_err[1], area_err[2], len_err[0], len_err[1], len_err[2]
        ),
    );
}

// ------------------------------------------------------------------- cut-003

fn strip_space_cond(eps: f64, ghost_gamma: f64, agg: Option<AggPolicy>) -> (f64, usize) {
    let grid = Quadtree::uniform(4);
    let h = 1.0 / 16.0;
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 0.5 + eps * h,
    };
    let params = FemParams {
        ghost_gamma,
        agg,
        strong_outer: true,
        ..FemParams::default()
    };
    let space = Space::build(&grid, &sdf, params).expect("strip builds");
    let (a, _b) = space.assemble(&|_, _| 1.0, &|_, _| 0.0);
    (condition_estimate(&a).cond, space.dof_count())
}

#[test]
fn cut_003_ghost_penalty_conditioning_independence() {
    let eps_sweep = [0.5, 1e-2, 1e-4, 1e-6, 1e-8];
    let mut with_ghost = Vec::new();
    let mut without = Vec::new();
    for &eps in &eps_sweep {
        with_ghost.push(strip_space_cond(eps, 0.5, None).0);
        without.push(strip_space_cond(eps, 0.0, None).0);
    }
    let ghost_ratio = with_ghost.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        / with_ghost.iter().copied().fold(f64::INFINITY, f64::min);
    let blowup = without[4] / without[0];
    let mut curves = String::new();
    for (i, &eps) in eps_sweep.iter().enumerate() {
        let _ = write!(
            curves,
            "{{\"eps\":{eps:.1e},\"cond_ghost\":{:.3e},\"cond_no_ghost\":{:.3e}}},",
            with_ghost[i], without[i]
        );
    }
    let pass = ghost_ratio < 30.0 && blowup > 1e3;
    verdict(
        "cut-003",
        pass,
        &format!(
            "\"detail\":\"eigenvalue-verified conditioning vs cut fraction\",\
             \"curves\":[{}],\"ghost_max_over_min\":{ghost_ratio:.2},\
             \"no_ghost_blowup\":{blowup:.3e}",
            curves.trim_end_matches(',')
        ),
    );
}

// ------------------------------------------------------------------- cut-004

fn disk_mms_l2_h1(circle: &Circle, level: u32) -> (f64, f64) {
    let grid = Quadtree::uniform(level);
    let space = Space::build(&grid, circle, FemParams::default()).expect("disk builds");
    let sol = space.solve(&sin_f, &sin_u).expect("disk solves");
    space.l2_h1_error(circle, &sin_u, &sin_grad, &sol.nodal)
}

#[test]
fn cut_004_mms_orders_across_randomized_cuts() {
    let mut lcg = Lcg(0x1001_2026_0707_0051);
    let mut configs: Vec<Circle> = (0..6)
        .map(|_| Circle {
            center: [
                0.5 + 0.12 * (lcg.unit() - 0.5),
                0.5 + 0.12 * (lcg.unit() - 0.5),
            ],
            radius: 0.26 + 0.07 * lcg.unit(),
        })
        .collect();
    // The deliberate sliver: tangency to four grid lines at every level.
    configs.push(Circle {
        center: [0.5, 0.5],
        radius: 0.25 + 1e-9,
    });
    let mut l2_orders = Vec::new();
    let mut h1_orders = Vec::new();
    let mut rows = String::new();
    let mut worst_fine_l2 = 0.0f64;
    for (k, c) in configs.iter().enumerate() {
        let (l2a, h1a) = disk_mms_l2_h1(c, 4);
        let (l2b, h1b) = disk_mms_l2_h1(c, 5);
        let (l2c, h1c) = disk_mms_l2_h1(c, 6);
        let ol2 = (l2b / l2c).log2();
        let oh1 = (h1b / h1c).log2();
        l2_orders.push(ol2);
        h1_orders.push(oh1);
        worst_fine_l2 = worst_fine_l2.max(l2c);
        let _ = write!(
            rows,
            "{{\"config\":{k},\"center\":[{:.4},{:.4}],\"r\":{:.6},\
             \"l2\":[{l2a:.3e},{l2b:.3e},{l2c:.3e}],\"h1\":[{h1a:.3e},{h1b:.3e},{h1c:.3e}],\
             \"l2_order\":{ol2:.2},\"h1_order\":{oh1:.2}}},",
            c.center[0], c.center[1], c.radius
        );
    }
    let median = |v: &mut Vec<f64>| -> f64 {
        v.sort_by(|a, b| a.partial_cmp(b).expect("finite orders"));
        v[v.len() / 2]
    };
    let med_l2 = median(&mut l2_orders);
    let med_h1 = median(&mut h1_orders);
    // The median alone can mask a SINGLE non-converging config — precisely the
    // pathological sliver this test deliberately includes. The "(sliver
    // included)" claim is only honest if EVERY config converges, so gate the
    // WORST order too (median() sorts the order vecs ascending in place, so
    // index 0 is the min). Floors are set to catch genuine non-convergence
    // (order → 0..1) without being brittle to the 2-level estimate's noise;
    // observed worst is ~2.16 L2 / ~1.04 H1 across all 7 configs.
    let worst_l2_order = l2_orders[0];
    let worst_h1_order = h1_orders[0];
    let pass = med_l2 >= 1.8
        && med_h1 >= 0.85
        && worst_fine_l2 < 2e-3
        && worst_l2_order >= 1.5
        && worst_h1_order >= 0.8;
    verdict(
        "cut-004",
        pass,
        &format!(
            "\"detail\":\"MMS orders across 7 randomized cut configs (sliver included)\",\
             \"configs\":[{}],\"median_l2_order\":{med_l2:.2},\
             \"median_h1_order\":{med_h1:.2},\"worst_l2_order\":{worst_l2_order:.2},\
             \"worst_h1_order\":{worst_h1_order:.2},\"worst_fine_l2\":{worst_fine_l2:.3e}",
            rows.trim_end_matches(',')
        ),
    );
}

// ------------------------------------------------------------------- cut-005

/// Body-fitted Q1 reference on (0, a) × (0, 1) with homogeneous
/// Dirichlet data (the shared fixture's exact solution vanishes on
/// every side): returns the L2 error at n×n cells.
#[allow(clippy::too_many_lines)]
fn body_fitted_l2(
    a_len: f64,
    n: usize,
    u: &dyn Fn(f64, f64) -> f64,
    f: &dyn Fn(f64, f64) -> f64,
) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let hx = a_len / n as f64;
    #[allow(clippy::cast_precision_loss)]
    let hy = 1.0 / n as f64;
    let idx = |i: usize, j: usize| -> Option<usize> {
        if i == 0 || i == n || j == 0 || j == n {
            None
        } else {
            Some((i - 1) + (j - 1) * (n - 1))
        }
    };
    let g3: [(f64, f64); 3] = [
        (-0.774_596_669_241_483_4, 0.555_555_555_555_555_6),
        (0.0, 0.888_888_888_888_889),
        (0.774_596_669_241_483_4, 0.555_555_555_555_555_6),
    ];
    let q1 = |lo: [f64; 2], hi: [f64; 2], p: [f64; 2]| -> ([f64; 4], [[f64; 2]; 4]) {
        let hx = hi[0] - lo[0];
        let hy = hi[1] - lo[1];
        let xi = (p[0] - lo[0]) / hx;
        let et = (p[1] - lo[1]) / hy;
        (
            [
                (1.0 - xi) * (1.0 - et),
                xi * (1.0 - et),
                xi * et,
                (1.0 - xi) * et,
            ],
            [
                [-(1.0 - et) / hx, -(1.0 - xi) / hy],
                [(1.0 - et) / hx, -xi / hy],
                [et / hx, xi / hy],
                [-et / hx, (1.0 - xi) / hy],
            ],
        )
    };
    let nf = (n - 1) * (n - 1);
    let mut coo = Coo::new(nf, nf);
    let mut rhs = vec![0.0f64; nf];
    for ci in 0..n {
        for cj in 0..n {
            #[allow(clippy::cast_precision_loss)]
            let lo = [ci as f64 * hx, cj as f64 * hy];
            let hi = [lo[0] + hx, lo[1] + hy];
            let local = [(ci, cj), (ci + 1, cj), (ci + 1, cj + 1), (ci, cj + 1)];
            let mut k = [[0.0f64; 4]; 4];
            let mut fl = [0.0f64; 4];
            for &(gx, wx) in &g3 {
                for &(gy, wy) in &g3 {
                    let p = [
                        f64::midpoint(lo[0], hi[0]) + 0.5 * hx * gx,
                        f64::midpoint(lo[1], hi[1]) + 0.5 * hy * gy,
                    ];
                    let w = wx * wy * 0.25 * hx * hy;
                    let (nv, gr) = q1(lo, hi, p);
                    let fv = f(p[0], p[1]);
                    for a in 0..4 {
                        for b in 0..4 {
                            k[a][b] += w * (gr[a][0] * gr[b][0] + gr[a][1] * gr[b][1]);
                        }
                        fl[a] += w * fv * nv[a];
                    }
                }
            }
            for a in 0..4 {
                let Some(ia) = idx(local[a].0, local[a].1) else {
                    continue;
                };
                rhs[ia] += fl[a];
                for b in 0..4 {
                    if let Some(ib) = idx(local[b].0, local[b].1) {
                        coo.push(ia, ib, k[a][b]);
                    }
                }
            }
        }
    }
    let op = CsrOp::symmetric(coo.assemble());
    let m = IdentityPrecond;
    let mut st = CgState::new(&op, &m, &rhs);
    let _ = st.run(&op, &m, 1e-12, 60_000);
    assert!(st.rel_residual() < 1e-8, "body-fitted CG converged");
    let mut l2 = 0.0f64;
    for ci in 0..n {
        for cj in 0..n {
            #[allow(clippy::cast_precision_loss)]
            let lo = [ci as f64 * hx, cj as f64 * hy];
            let hi = [lo[0] + hx, lo[1] + hy];
            let local = [(ci, cj), (ci + 1, cj), (ci + 1, cj + 1), (ci, cj + 1)];
            let vals: Vec<f64> = local
                .iter()
                .map(|&(i, j)| idx(i, j).map_or(0.0, |id| st.x[id]))
                .collect();
            for &(gx, wx) in &g3 {
                for &(gy, wy) in &g3 {
                    let p = [
                        f64::midpoint(lo[0], hi[0]) + 0.5 * hx * gx,
                        f64::midpoint(lo[1], hi[1]) + 0.5 * hy * gy,
                    ];
                    let w = wx * wy * 0.25 * hx * hy;
                    let (nv, _) = q1(lo, hi, p);
                    let uh: f64 = (0..4).map(|a| nv[a] * vals[a]).sum();
                    let e = u(p[0], p[1]) - uh;
                    l2 += w * e * e;
                }
            }
        }
    }
    l2.sqrt()
}

#[test]
fn cut_005_embedded_bc_matches_body_fitted() {
    let a_len = 0.55;
    let u = move |x: f64, y: f64| (PI * x / a_len).sin() * (PI * y).sin();
    let f = move |x: f64, y: f64| (PI * PI / (a_len * a_len) + PI * PI) * u(x, y);
    let grad = move |x: f64, y: f64| {
        [
            PI / a_len * (PI * x / a_len).cos() * (PI * y).sin(),
            PI * (PI * x / a_len).sin() * (PI * y).cos(),
        ]
    };
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: a_len,
    };
    let mut ratios = Vec::new();
    let mut cut_l2s = Vec::new();
    let mut rows = String::new();
    for (level, n) in [(4u32, 16usize), (5, 32)] {
        let grid = Quadtree::uniform(level);
        let params = FemParams {
            strong_outer: true,
            ..FemParams::default()
        };
        let space = Space::build(&grid, &sdf, params).expect("strip builds");
        let sol = space.solve(&f, &u).expect("strip solves");
        let (l2, _h1) = space.l2_h1_error(&sdf, &u, &grad, &sol.nodal);
        let bf = body_fitted_l2(a_len, n, &u, &f);
        ratios.push(l2 / bf);
        cut_l2s.push(l2);
        let _ = write!(
            rows,
            "{{\"level\":{level},\"cutfem_l2\":{l2:.3e},\"body_fitted_l2\":{bf:.3e},\
             \"ratio\":{:.2}}},",
            l2 / bf
        );
    }
    let order = (cut_l2s[0] / cut_l2s[1]).log2();
    let pass = ratios.iter().all(|&r| r < 3.0) && order > 1.7;
    verdict(
        "cut-005",
        pass,
        &format!(
            "\"detail\":\"embedded Nitsche vs body-fitted Q1 on the shared strip fixture\",\
             \"rows\":[{}],\"cutfem_order\":{order:.2}",
            rows.trim_end_matches(',')
        ),
    );
}

// ------------------------------------------------------------------- cut-006

#[test]
fn cut_006_moving_interface_without_remeshing() {
    let grid = Quadtree::uniform(5); // built ONCE; never touched again
    let leaves_before = grid.leaf_count();
    let mut errs = Vec::new();
    let mut rows = String::new();
    for k in 0..11 {
        let circle = Circle {
            center: [0.35 + 0.03 * f64::from(k), 0.5],
            radius: 0.3,
        };
        let space = Space::build(&grid, &circle, FemParams::default()).expect("step builds");
        let sol = space.solve(&sin_f, &sin_u).expect("step solves");
        let (l2, _) = space.l2_h1_error(&circle, &sin_u, &sin_grad, &sol.nodal);
        let _ = write!(
            rows,
            "{{\"step\":{k},\"cx\":{:.2},\"l2\":{l2:.3e},\"iters\":{}}},",
            0.35 + 0.03 * f64::from(k),
            sol.iters
        );
        errs.push(l2);
    }
    let emax = errs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let emin = errs.iter().copied().fold(f64::INFINITY, f64::min);
    let pass = leaves_before == grid.leaf_count() && emax < 8e-3 && emax / emin < 2.5;
    verdict(
        "cut-006",
        pass,
        &format!(
            "\"detail\":\"11-step level-set translation on one fixed background grid\",\
             \"steps\":[{}],\"err_max\":{emax:.3e},\"err_stability_ratio\":{:.2}",
            rows.trim_end_matches(','),
            emax / emin
        ),
    );
}

// ------------------------------------------------------------------- cut-007

#[test]
fn cut_007_aggregation_fallback() {
    // The pathological sliver from cut-003, ghost penalty OFF.
    let (cond_bare, _) = strip_space_cond(1e-8, 0.0, None);
    let (cond_agg, _) = strip_space_cond(1e-8, 0.0, Some(AggPolicy::default()));
    // Accuracy: the shared strip MMS solved with aggregation-only vs
    // ghost-only stabilization.
    let a_len = 0.5 + 1e-8 / 16.0;
    let u = move |x: f64, y: f64| (PI * x / a_len).sin() * (PI * y).sin();
    let f = move |x: f64, y: f64| (PI * PI / (a_len * a_len) + PI * PI) * u(x, y);
    let grad = move |x: f64, y: f64| {
        [
            PI / a_len * (PI * x / a_len).cos() * (PI * y).sin(),
            PI * (PI * x / a_len).sin() * (PI * y).cos(),
        ]
    };
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: a_len,
    };
    let grid = Quadtree::uniform(4);
    let solve_with = |ghost: f64, agg: Option<AggPolicy>| -> (f64, usize, Vec<String>) {
        let params = FemParams {
            ghost_gamma: ghost,
            agg,
            strong_outer: true,
            ..FemParams::default()
        };
        let space = Space::build(&grid, &sdf, params).expect("sliver builds");
        let sol = space.solve(&f, &u).expect("sliver solves");
        let (l2, _) = space.l2_h1_error(&sdf, &u, &grad, &sol.nodal);
        (l2, space.stats().aggregated, space.agg_log().to_vec())
    };
    let (l2_ghost, _, _) = solve_with(0.5, None);
    let (l2_agg, aggregated, log) = solve_with(0.0, Some(AggPolicy::default()));
    let pass = aggregated > 0 && cond_agg < cond_bare / 1e2 && l2_agg < 2.0 * l2_ghost;
    verdict(
        "cut-007",
        pass,
        &format!(
            "\"detail\":\"aggregation restores conditioning and accuracy with ghost OFF\",\
             \"cond_bare\":{cond_bare:.3e},\"cond_agg\":{cond_agg:.3e},\
             \"aggregated_nodes\":{aggregated},\"l2_ghost\":{l2_ghost:.3e},\
             \"l2_agg\":{l2_agg:.3e},\"policy_rows\":{}",
            log.len()
        ),
    );
    for row in log.iter().take(3) {
        println!("{{\"test\":\"cut-007\",\"policy\":{row}}}");
    }
}

// ------------------------------------------------------------------- cut-008

#[test]
fn cut_008_hanging_nodes_and_adaptivity() {
    // (a) Patch test: a linear exact solution is reproduced EXACTLY on
    // a graded tree with hanging nodes (straight interface, so chords
    // are exact too).
    let lin = |x: f64, y: f64| 1.0 + 2.0 * x + 3.0 * y;
    let sdf = HalfPlane {
        normal: [0.6, 0.8],
        offset: 0.55,
    };
    let mut grid = Quadtree::with_room(3, 5);
    grid.refine_toward_interface(&sdf, 5);
    let params = FemParams {
        strong_outer: true,
        solver_tol: 1e-13,
        ..FemParams::default()
    };
    let space = Space::build(&grid, &sdf, params).expect("graded builds");
    let hanging = space.stats().hanging;
    let sol = space.solve(&|_, _| 0.0, &lin).expect("patch solves");
    let mut max_nodal = 0.0f64;
    for (&n, &v) in &sol.nodal {
        let p = space_node_pos(&grid, n);
        max_nodal = max_nodal.max((v - lin(p[0], p[1])).abs());
    }
    // (b) Estimate-driven h-refinement efficiency: a bump on the
    // interface; graded tree beats uniform-coarse on error at a
    // fraction of uniform-fine's DOFs.
    let circle = Circle {
        center: [0.5, 0.5],
        radius: 0.3,
    };
    let beta = 200.0;
    let px = 0.8;
    let py = 0.5;
    let u = move |x: f64, y: f64| {
        let r2 = (x - px) * (x - px) + (y - py) * (y - py);
        sin_u(x, y) + (-beta * r2).exp()
    };
    let f = move |x: f64, y: f64| {
        let r2 = (x - px) * (x - px) + (y - py) * (y - py);
        sin_f(x, y) + (4.0 * beta - 4.0 * beta * beta * r2) * (-beta * r2).exp()
    };
    let grad = move |x: f64, y: f64| {
        let r2 = (x - px) * (x - px) + (y - py) * (y - py);
        let e = (-beta * r2).exp();
        let s = sin_grad(x, y);
        [
            s[0] - 2.0 * beta * (x - px) * e,
            s[1] - 2.0 * beta * (y - py) * e,
        ]
    };
    let solve_on = |grid: &Quadtree| -> (f64, usize) {
        let space = Space::build(grid, &circle, FemParams::default()).expect("builds");
        let sol = space.solve(&f, &u).expect("solves");
        let (l2, _) = space.l2_h1_error(&circle, &u, &grad, &sol.nodal);
        (l2, space.dof_count())
    };
    let (err_u4, _dofs_u4) = solve_on(&Quadtree::uniform(4));
    let (err_u6, dofs_u6) = solve_on(&Quadtree::uniform(6));
    let mut graded = Quadtree::with_room(4, 6);
    // The "estimate": refine where the data f is large (the bump) —
    // the dwr-adaptivity hook drives the same refine_where surface.
    graded.refine_where(6, &|lo: [f64; 2], hi: [f64; 2]| {
        let cx = f64::midpoint(lo[0], hi[0]);
        let cy = f64::midpoint(lo[1], hi[1]);
        let d2 = (cx - px) * (cx - px) + (cy - py) * (cy - py);
        d2 < 0.15 * 0.15
    });
    graded.refine_toward_interface(&circle, 6);
    let (err_graded, dofs_graded) = solve_on(&graded);
    // Efficiency claim, all measured: the graded tree lands within 2×
    // of uniform-fine ERROR at ~half its DOFs, and far below
    // uniform-coarse error.
    #[allow(clippy::cast_precision_loss)]
    let dof_ratio = dofs_graded as f64 / dofs_u6 as f64;
    let pass = hanging > 0
        && max_nodal < 1e-8
        && err_graded < 0.5 * err_u4
        && err_graded < 2.0 * err_u6
        && dof_ratio < 0.55;
    verdict(
        "cut-008",
        pass,
        &format!(
            "\"detail\":\"hanging-node patch test + estimate-driven refinement\",\
             \"hanging\":{hanging},\"patch_max_nodal_err\":{max_nodal:.3e},\
             \"err_uniform4\":{err_u4:.3e},\"err_uniform6\":{err_u6:.3e},\
             \"err_graded\":{err_graded:.3e},\"dofs_uniform6\":{dofs_u6},\
             \"dofs_graded\":{dofs_graded}"
        ),
    );
}

fn space_node_pos(grid: &Quadtree, n: (u32, u32)) -> [f64; 2] {
    grid.node_pos(n)
}
