//! Closed-box multi-panel ties, glue compliance, and fluid-boundary ownership tests (`M1.2-W2b`).
//!
//! Bead: `frankensim-irp-v5-instrument-realism-pexfp.51` (`[IRP] M1.2-W2b: Closed boxes + glue compliance`).
//!
//! DONE-WHEN verification gates:
//! - `g3_rigid_matching_box_equals_monolithic` (G3): 6-panel closed rectangular box (12 edge ties,
//!   8 crosspoints) matches one conforming closed-box `ShellMesh`: first 30 eigenvalues, including
//!   the 6 rigid modes, agree within `1e-10` relative.
//! - `g1_nonmatching_box_converges` (G1): nonmatching panel meshes (2:3 ratios) over >= 3 rungs:
//!   eigenvalue error against the monolithic fine solution converges at an observed order within
//!   `0.2` of the monolithic order (`fs-mms` `OrderGate`).
//! - `g1_glue_compliance_rigid_limit` (G1): compliant glue lines have area-specific normal and
//!   shear stiffness `k`, scaled x10 over >= 4 rungs; every eigenvalue increases monotonically
//!   with `k` (Rayleigh, G0), and the spectrum converges to the rigid-tie spectrum with an
//!   observed slope of `-1 +/- 0.2` in `1/k`.
//! - `g0_glue_compliance_is_data` (G0): glue stiffness and loss come only from an ordered
//!   `fs_matdb::InterfaceSystemCard` with `SystemContext.third_body` = adhesive and a named
//!   history state, resolved with receipts; missing property refuses; fixture cards are named
//!   as authored fixtures.
//! - `g0_glue_loss_nonnegative_and_counted_once` (G0): Kelvin–Voigt glue viscosity per area gives
//!   a PSD damping operator (per-mode loss >= 0), vanishing on rigid-body modes and counted once.
//! - `g0_fluid_boundary_ownership_audit` (G0): every skin face of the closed box (both faces of
//!   each panel) is owned by exactly one acoustic domain (`InteriorCavity(id)` or `Exterior`);
//!   glue-overlap regions are owned by no fluid; open apertures are explicit boundary loops;
//!   missing owner, double owner, orientation conflict, and unclosed interior each refuse.
//! - `g0_interior_volume_from_skins` (G0): for a midsurface box `a x b x c` with uniform wall `h`,
//!   the inner skins plus the closure faces the audit builds at tied edges and corners enclose
//!   `(a - h)(b - h)(c - h)` to `1e-12` relative.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel,
};
use fs_math::det;
use fs_mms::{LadderSide, OrderGate, RefinementLadder, fit_order};
use fs_modal::{SliceOptions, slice_window};
use fs_plate::PlateSection;
use fs_plate::shell::junction::{
    EdgeTie, FaceTie, GlueLine, RotationCoupling, ShellComponent, TieBudget, TieKind, TieOffset,
    TieSide, TiedShellModel, tie_shells,
};
use fs_plate::shell::ownership::{
    AcousticDomainId, ApertureLoop, FluidBoundaryAssignments, FluidBoundaryOwnership, SkinFaceId,
};
use fs_plate::shell::reduction::radiation::{
    RadiationSurfaceBudget, ShellFace, ShellRadiationSurface,
};
use fs_plate::shell::reduction::{ReductionBudget, ShellReduction};
use fs_plate::shell::{ShellMesh, ShellModel, ShellSupport, assemble_shell, modes_shell};
use fs_qty::{Dims, DynViscosity, Length, Pressure};
use std::collections::BTreeMap;

/// Structured rectangular panel mesh on `[0, u_len] x [0, v_len]` mapped into 3-D by `map_uv(u, v)`.
/// If `flip_winding` is true, triangle vertices `(1, 2)` are swapped so the facet normal reverses.
fn build_panel_mesh(
    nu: usize,
    nv: usize,
    flip_winding: bool,
    map_uv: impl Fn(f64, f64) -> [f64; 3],
) -> ShellMesh {
    let mut nodes = Vec::with_capacity((nu + 1) * (nv + 1));
    for j in 0..=nv {
        let v = (j as f64) / (nv as f64);
        for i in 0..=nu {
            let u = (i as f64) / (nu as f64);
            nodes.push(map_uv(u, v));
        }
    }
    let idx = |i: usize, j: usize| j * (nu + 1) + i;
    let mut tris = Vec::with_capacity(2 * nu * nv);
    for j in 0..nv {
        for i in 0..nu {
            let n00 = idx(i, j);
            let n10 = idx(i + 1, j);
            let n11 = idx(i + 1, j + 1);
            let n01 = idx(i, j + 1);
            if flip_winding {
                tris.push([n00, n11, n10]);
                tris.push([n00, n01, n11]);
            } else {
                tris.push([n00, n10, n11]);
                tris.push([n00, n11, n01]);
            }
        }
    }
    ShellMesh::new(nodes, tris).expect("valid rectangular panel ShellMesh")
}

/// Extract the 4 ordered boundary polylines of a `(nu, nv)` panel mesh:
/// - `bottom`: `v = 0`, `i = 0..=nu`
/// - `top`: `v = 1`, `i = 0..=nu`
/// - `left`: `u = 0`, `j = 0..=nv`
/// - `right`: `u = 1`, `j = 0..=nv`
fn panel_boundary_edges(nu: usize, nv: usize) -> [Vec<usize>; 4] {
    let idx = |i: usize, j: usize| j * (nu + 1) + i;
    let bottom: Vec<usize> = (0..=nu).map(|i| idx(i, 0)).collect();
    let top: Vec<usize> = (0..=nu).map(|i| idx(i, nv)).collect();
    let left: Vec<usize> = (0..=nv).map(|j| idx(0, j)).collect();
    let right: Vec<usize> = (0..=nv).map(|j| idx(nu, j)).collect();
    [bottom, top, left, right]
}

struct SixPanelBox {
    meshes: [ShellMesh; 6],
    models: [ShellModel; 6],
    thicknesses: [Vec<f64>; 6],
    ties: Vec<EdgeTie>,
}

/// Build a 6-panel closed rectangular box `a x b x c` with wall thickness `h` and panel subdivisions:
/// - Panel 0 (Bottom, `z = 0`, outward normal `-z`): `nx_tb x ny_tb`
/// - Panel 1 (Top, `z = c`, outward normal `+z`): `nx_tb x ny_tb`
/// - Panel 2 (Front, `y = 0`, outward normal `-y`): `nx_fb x nz_side`
/// - Panel 3 (Back, `y = b`, outward normal `+y`): `nx_fb x nz_side`
/// - Panel 4 (Left, `x = 0`, outward normal `-x`): `ny_lr x nz_side`
/// - Panel 5 (Right, `x = a`, outward normal `+x`): `ny_lr x nz_side`
///
/// And 12 `EdgeTie`s meeting at the 8 three-way box corners.
#[allow(clippy::too_many_arguments)]
fn build_six_panel_box(
    a: f64,
    b: f64,
    c: f64,
    h: f64,
    sec: &PlateSection,
    nx_tb: usize,
    ny_tb: usize,
    nx_fb: usize,
    ny_lr: usize,
    nz_side: usize,
    kind: TieKind,
) -> SixPanelBox {
    // 0: Bottom (z = 0), u = x/a, v = y/b, flip_winding = true -> normal (0, 0, -1)
    let m0 = build_panel_mesh(nx_tb, ny_tb, true, |u, v| [u * a, v * b, 0.0]);
    // 1: Top (z = c), u = x/a, v = y/b, flip_winding = false -> normal (0, 0, +1)
    let m1 = build_panel_mesh(nx_tb, ny_tb, false, |u, v| [u * a, v * b, c]);
    // 2: Front (y = 0), u = x/a, v = z/c, flip_winding = false -> normal (0, -1, 0)
    let m2 = build_panel_mesh(nx_fb, nz_side, false, |u, v| [u * a, 0.0, v * c]);
    // 3: Back (y = b), u = x/a, v = z/c, flip_winding = true -> normal (0, +1, 0)
    let m3 = build_panel_mesh(nx_fb, nz_side, true, |u, v| [u * a, b, v * c]);
    // 4: Left (x = 0), u = y/b, v = z/c, flip_winding = true -> normal (-1, 0, 0)
    let m4 = build_panel_mesh(ny_lr, nz_side, true, |u, v| [0.0, u * b, v * c]);
    // 5: Right (x = a), u = y/b, v = z/c, flip_winding = false -> normal (+1, 0, 0)
    let m5 = build_panel_mesh(ny_lr, nz_side, false, |u, v| [a, u * b, v * c]);

    let meshes = [m0, m1, m2, m3, m4, m5];
    let models = meshes
        .each_ref()
        .map(|m| assemble_shell(m, sec, &[], ShellSupport::Free).expect("assemble panel shell"));
    let thicknesses = meshes.each_ref().map(|m| vec![h; m.nodes.len()]);

    let [b0_bot, b0_top, b0_left, b0_right] = panel_boundary_edges(nx_tb, ny_tb);
    let [b1_bot, b1_top, b1_left, b1_right] = panel_boundary_edges(nx_tb, ny_tb);
    let [b2_bot, b2_top, b2_left, b2_right] = panel_boundary_edges(nx_fb, nz_side);
    let [b3_bot, b3_top, b3_left, b3_right] = panel_boundary_edges(nx_fb, nz_side);
    let [b4_bot, b4_top, b4_left, b4_right] = panel_boundary_edges(ny_lr, nz_side);
    let [b5_bot, b5_top, b5_left, b5_right] = panel_boundary_edges(ny_lr, nz_side);

    let make_tie = |ca: usize, na: Vec<usize>, cb: usize, nb: Vec<usize>| EdgeTie {
        a: TieSide {
            component: ca,
            nodes: na,
        },
        b: TieSide {
            component: cb,
            nodes: nb,
        },
        offset: TieOffset::Midsurface,
        rotations: RotationCoupling::Full,
        kind,
    };

    let ties = vec![
        // 4 bottom horizontal edges (master = Bottom 0)
        make_tie(0, b0_bot, 2, b2_bot),
        make_tie(0, b0_top, 3, b3_bot),
        make_tie(0, b0_left, 4, b4_bot),
        make_tie(0, b0_right, 5, b5_bot),
        // 4 top horizontal edges (master = Top 1)
        make_tie(1, b1_bot, 2, b2_top),
        make_tie(1, b1_top, 3, b3_top),
        make_tie(1, b1_left, 4, b4_top),
        make_tie(1, b1_right, 5, b5_top),
        // 4 vertical edges (master = Front 2 or Back 3, slave = Left 4 or Right 5)
        make_tie(2, b2_left, 4, b4_left),
        make_tie(2, b2_right, 5, b5_left),
        make_tie(3, b3_left, 4, b4_right),
        make_tie(3, b3_right, 5, b5_right),
    ];

    SixPanelBox {
        meshes,
        models,
        thicknesses,
        ties,
    }
}

/// Build the single welded conforming `ShellMesh` for the same `a x b x c` box with matching
/// `(nx, ny, nz)` subdivisions by welding coincident vertices across the 6 panels.
fn build_monolithic_box_mesh(
    a: f64,
    b: f64,
    c: f64,
    h: f64,
    sec: &PlateSection,
    nx: usize,
    ny: usize,
    nz: usize,
) -> ShellMesh {
    let six = build_six_panel_box(a, b, c, h, sec, nx, ny, nx, ny, nz, TieKind::Rigid);
    let quant = |x: f64| -> i64 { (x * 1e9).round() as i64 };
    let mut node_map: BTreeMap<[i64; 3], usize> = BTreeMap::new();
    let mut nodes: Vec<[f64; 3]> = Vec::new();
    let mut tris: Vec<[usize; 3]> = Vec::new();

    for panel_mesh in &six.meshes {
        let mut local_to_global = Vec::with_capacity(panel_mesh.nodes.len());
        for &p in &panel_mesh.nodes {
            let key = [quant(p[0]), quant(p[1]), quant(p[2])];
            let g_idx = *node_map.entry(key).or_insert_with(|| {
                let id = nodes.len();
                nodes.push(p);
                id
            });
            local_to_global.push(g_idx);
        }
        for &tri in &panel_mesh.tris {
            tris.push([
                local_to_global[tri[0]],
                local_to_global[tri[1]],
                local_to_global[tri[2]],
            ]);
        }
    }
    ShellMesh::new(nodes, tris).expect("valid monolithic closed-box ShellMesh")
}

fn add_scaled_csr(k: &fs_sparse::Csr, m: &fs_sparse::Csr, sigma: f64) -> fs_sparse::Csr {
    let n = k.nrows();
    let mut coo = fs_sparse::Coo::new(n, n);
    for r in 0..n {
        let (k_cols, k_vals) = k.row(r);
        for (&c, &v) in k_cols.iter().zip(k_vals.iter()) {
            coo.push(r, c, v);
        }
        let (m_cols, m_vals) = m.row(r);
        for (&c, &v) in m_cols.iter().zip(m_vals.iter()) {
            coo.push(r, c, sigma * v);
        }
    }
    coo.assemble()
}

fn find_tight_upper_window_bound(
    k_shifted: &fs_sparse::Csr,
    m: &fs_sparse::Csr,
    min_count: usize,
    start_hi: f64,
) -> f64 {
    let sample_pencil = add_scaled_csr(k_shifted, m, -start_hi);
    let symbolic =
        fs_sparse::SymbolicLdlt::analyze(&sample_pencil, fs_sparse::DirectOrdering::Amd).unwrap();
    let count_at = |mu: f64| -> usize {
        let pencil = add_scaled_csr(k_shifted, m, -mu);
        symbolic
            .factor(&pencil, &fs_sparse::LdltOptions::default())
            .unwrap()
            .inertia()
            .negative
    };
    let mut lo = 0.0_f64;
    let mut hi = start_hi;
    for _ in 0..40 {
        let neg = count_at(hi);
        if neg >= min_count {
            break;
        }
        lo = hi;
        hi *= 2.0;
    }
    for _ in 0..6 {
        if count_at(hi) <= min_count + 6 {
            break;
        }
        let mid = 0.5 * (lo + hi);
        if count_at(mid) >= min_count {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    hi
}

/// Solve the first `n_modes` eigenvalues (and mass-normalized eigenvectors) of a `ShellModel`
/// using a positive shift `sigma` and an inertia-bracketed upper bound so the 6 rigid-body modes
/// (`lambda ~ 0`) and lowest flexible modes are resolved in a single pass of `fs_modal::slice_window`.
fn solve_first_modes(
    model: &ShellModel,
    n_modes: usize,
    start_hi: f64,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let sigma = 1.0e5_f64;
    let k_shifted = add_scaled_csr(&model.k, &model.m, sigma);
    let hi = find_tight_upper_window_bound(&k_shifted, &model.m, n_modes, start_hi + sigma);
    let opts = SliceOptions {
        ritz_tol: 1e-13,
        max_lanczos: 1024,
        max_restarts: 64,
        mass_diagonal_equilibration: true,
        ..SliceOptions::default()
    };
    let report = slice_window(&k_shifted, &model.m, (0.1 * sigma, hi), &opts)
        .expect("shifted modal slice succeeds");
    assert!(
        report.modes.len() >= n_modes,
        "expected at least {n_modes} modes in window, found {}",
        report.modes.len()
    );
    let lambdas: Vec<f64> = report.modes[..n_modes]
        .iter()
        .map(|m| m.lambda - sigma)
        .collect();
    let phis: Vec<Vec<f64>> = report.modes[..n_modes]
        .iter()
        .map(|m| m.phi.clone())
        .collect();
    (lambdas, phis)
}

#[test]
fn g3_rigid_matching_box_equals_monolithic() {
    let (a, b, c, h) = (0.24, 0.18, 0.12, 0.003);
    let sec = PlateSection::isotropic(10.0e9, 0.30, h, 500.0).unwrap();
    let (nx, ny, nz) = (3, 3, 2);

    // 1. Monolithic closed-box ShellMesh.
    let mono_mesh = build_monolithic_box_mesh(a, b, c, h, &sec, nx, ny, nz);
    let mono_model = assemble_shell(&mono_mesh, &sec, &[], ShellSupport::Free).unwrap();

    // 2. 6-panel tied box with 12 EdgeTies and 8 three-way corner crosspoints.
    let six = build_six_panel_box(a, b, c, h, &sec, nx, ny, nx, ny, nz, TieKind::Rigid);
    let comps: Vec<ShellComponent<'_>> = (0..6)
        .map(|i| ShellComponent {
            id: i,
            mesh: &six.meshes[i],
            model: &six.models[i],
            nodal_thickness_m: &six.thicknesses[i],
        })
        .collect();
    let budget = TieBudget {
        max_segments: 256,
        max_slave_dofs: 1024,
        gap_tolerance_m: 1e-9,
        angle_tolerance_rad: 1e-6,
    };
    let tied = tie_shells(&comps, &six.ties, budget, || false).expect("6-panel tied box succeeds");

    assert_eq!(
        tied.model.free, mono_model.free,
        "6-panel tied box with 8 three-way corner crosspoints must have identical free DOFs to monolithic box"
    );

    let (lam_mono, _) = solve_first_modes(&mono_model, 30, 8.0e8);
    let (lam_tied, _) = solve_first_modes(&tied.model, 30, 8.0e8);

    // First 6 modes are rigid-body modes (~0); mode 6 is the first flexible mode.
    let ref_scale = lam_mono[6].max(1.0);
    let mut max_rel_err = 0.0_f64;
    for k in 0..30 {
        let denom = if k < 6 {
            ref_scale
        } else {
            lam_mono[k].abs().max(ref_scale)
        };
        let rel = (lam_tied[k] - lam_mono[k]).abs() / denom;
        max_rel_err = max_rel_err.max(rel);
        assert!(
            rel <= 1e-10,
            "mode {k}: tied lambda {:.6e} vs monolithic {:.6e}, rel_err {rel:.3e} > 1e-10",
            lam_tied[k],
            lam_mono[k]
        );
    }
    println!(
        "G3 rigid matching 6-panel box vs monolithic: free={}, max_rel_err_30_modes={max_rel_err:.3e} <= 1e-10",
        tied.model.free
    );
}

#[test]
fn g1_nonmatching_box_converges() {
    let (a, b, c, h) = (0.24, 0.18, 0.12, 0.004);
    let sec = PlateSection::isotropic(10.0e9, 0.30, h, 500.0).unwrap();
    let budget = TieBudget {
        max_segments: 2048,
        max_slave_dofs: 8192,
        gap_tolerance_m: 1e-9,
        angle_tolerance_rad: 1e-6,
    };

    // The eigenvalue error is measured against a monolithic fine solution, not
    // an extrapolation that would presume the order being measured. Rungs
    // m = 4, 6, 8 put 8/12/16 cells on every monolithic panel edge (12/18/24 on
    // the 3-sides of the nonmatching box). The reference is the monolithic box
    // at m = 32, 64 cells per edge, four times the finest coarse side. The
    // first ten flexible modes are gated, each against its own monolithic
    // order.
    //
    // At m = 2 (4 cells per edge) the errors were 12-22 % and the close pairs
    // crossed, so the ladder starts at m = 4. Flexible modes 9 and 10 are a
    // near pair (86.82e6 and 87.68e6, 1.0 % apart) whose errors exceed that gap
    // at every affordable rung, so sorted index hands the 9th eigenvalue to
    // whichever family is lower on each mesh. Consecutive reference eigenvalues
    // whose gap is below the coarsest-rung error of either therefore form one
    // cluster, gated by its mean (the cluster trace), which no family swap can
    // reorder. A singleton cluster is the eigenvalue itself.
    let rungs = [4usize, 6, 8];
    let m_ref = 32usize;
    let flexible = 10usize;
    let count = 6 + flexible;
    let mono_lambdas = |m: usize, start_hi: f64| -> Vec<f64> {
        let mesh = build_monolithic_box_mesh(a, b, c, h, &sec, 2 * m, 2 * m, 2 * m);
        let model = assemble_shell(&mesh, &sec, &[], ShellSupport::Free).unwrap();
        solve_first_modes(&model, count, start_hi).0
    };
    let nonmatch_lambdas = |m: usize| -> Vec<f64> {
        // 2:3 nonmatching panel mesh ratio along all 8 horizontal box edges:
        // Top/Bottom use 2*m x 2*m; Front/Back/Left/Right use 3*m x 2*m.
        let six = build_six_panel_box(
            a,
            b,
            c,
            h,
            &sec,
            2 * m,
            2 * m,
            3 * m,
            3 * m,
            2 * m,
            TieKind::Rigid,
        );
        let comps: Vec<ShellComponent<'_>> = (0..6)
            .map(|i| ShellComponent {
                id: i,
                mesh: &six.meshes[i],
                model: &six.models[i],
                nodal_thickness_m: &six.thicknesses[i],
            })
            .collect();
        let tied = tie_shells(&comps, &six.ties, budget, || false).unwrap();
        solve_first_modes(&tied.model, count, 2.5e8).0
    };

    let mono_raw: Vec<Vec<f64>> = rungs.iter().map(|&m| mono_lambdas(m, 2.5e8)).collect();
    let nonmatch_raw: Vec<Vec<f64>> = rungs.iter().map(|&m| nonmatch_lambdas(m)).collect();
    // The reference converges from the finest rung's side; seeding its window
    // there spares the bisection factorizations of the 148k-DOF reference.
    let reference = mono_lambdas(m_ref, 1.1 * mono_raw.last().unwrap()[count - 1]);

    let flex = |lams: &[f64]| lams[6..count].to_vec();
    let ref_flex = flex(&reference);
    let coarse_err: Vec<f64> = flex(&mono_raw[0])
        .iter()
        .zip(&ref_flex)
        .map(|(l, r)| (l - r).abs() / r)
        .collect();
    let mut clusters: Vec<Vec<usize>> = vec![vec![0]];
    for j in 1..flexible {
        let gap = ref_flex[j] / ref_flex[j - 1] - 1.0;
        if gap < coarse_err[j].max(coarse_err[j - 1]) {
            clusters.last_mut().unwrap().push(j);
        } else {
            clusters.push(vec![j]);
        }
    }
    let mean = |lams: &[f64], cluster: &[usize]| {
        cluster.iter().map(|&j| lams[j]).sum::<f64>() / cluster.len() as f64
    };
    let hs: Vec<f64> = rungs.iter().map(|&m| 1.0 / m as f64).collect();

    let mut failures = Vec::new();
    for cluster in &clusters {
        let target = mean(&ref_flex, cluster);
        let ladder_of = |raw: &[Vec<f64>]| -> Vec<f64> {
            raw.iter()
                .map(|l| (mean(&flex(l), cluster) - target).abs() / target)
                .collect()
        };
        let mono = ladder_of(&mono_raw);
        let nonmatch = ladder_of(&nonmatch_raw);
        let own = fit_order(
            &RefinementLadder::new(hs.clone(), mono.clone()).expect("valid monolithic ladder"),
        );
        let ladder =
            RefinementLadder::new(hs.clone(), nonmatch.clone()).expect("valid nonmatching ladder");
        let members: Vec<String> = cluster.iter().map(|j| (j + 1).to_string()).collect();
        let case = format!(
            "nonmatching-2-to-3-box-flexible-modes-{}",
            members.join("-")
        );
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W2b\",\"test\":\"g1_nonmatching_box_converges\",\"tier\":\"G1\",\"case\":\"{case}\",\"rungs_m\":{rungs:?},\"reference_m\":{m_ref},\"reference_lambdas\":{:?},\"monolithic_lambdas\":{:?},\"nonmatching_lambdas\":{:?},\"monolithic_order\":{},\"monolithic_errors\":{mono:?},\"nonmatching_errors\":{nonmatch:?}}}",
            cluster.iter().map(|&j| ref_flex[j]).collect::<Vec<_>>(),
            mono_raw
                .iter()
                .map(|l| cluster.iter().map(|&j| l[6 + j]).collect::<Vec<_>>())
                .collect::<Vec<_>>(),
            nonmatch_raw
                .iter()
                .map(|l| cluster.iter().map(|&j| l[6 + j]).collect::<Vec<_>>())
                .collect::<Vec<_>>(),
            own.observed,
        );
        match (OrderGate {
            theoretical: own.observed,
        })
        .check(&case, LadderSide::Primal, &ladder)
        {
            Ok(verdict) => println!("{}", verdict.json_line(true)),
            Err(refusal) => {
                println!("{}", refusal.detail());
                failures.push(case);
            }
        }
    }
    assert!(
        failures.is_empty(),
        "nonmatching box order gate failed for {failures:?}"
    );
}

#[test]
fn g1_glue_compliance_rigid_limit() {
    let (a, b, c, h) = (0.24, 0.18, 0.12, 0.003);
    let sec = PlateSection::isotropic(10.0e9, 0.30, h, 500.0).unwrap();
    let (nx, ny, nz) = (2, 2, 2);
    let budget = TieBudget {
        max_segments: 256,
        max_slave_dofs: 2048,
        gap_tolerance_m: 1e-9,
        angle_tolerance_rad: 1e-6,
    };

    // Rigid reference spectrum.
    let six_rigid = build_six_panel_box(a, b, c, h, &sec, nx, ny, nx, ny, nz, TieKind::Rigid);
    let comps_rigid: Vec<ShellComponent<'_>> = (0..6)
        .map(|i| ShellComponent {
            id: i,
            mesh: &six_rigid.meshes[i],
            model: &six_rigid.models[i],
            nodal_thickness_m: &six_rigid.thicknesses[i],
        })
        .collect();
    let tied_rigid = tie_shells(&comps_rigid, &six_rigid.ties, budget, || false).unwrap();
    let (lam_rigid, _) = solve_first_modes(&tied_rigid.model, 14, 5.0e8);

    // 4 rungs of glue stiffness k scaled x10: k = k0 * {1, 10, 100, 1000}.
    let k0 = 1.0e12_f64;
    let scales = [1.0_f64, 10.0, 100.0, 1000.0];
    let mut spectra: Vec<Vec<f64>> = Vec::with_capacity(scales.len());

    for &s in &scales {
        let glue = GlueLine {
            normal_n_m3: k0 * s,
            shear_n_m3: [0.5 * k0 * s, 0.5 * k0 * s],
            viscosity_n_s_m3: None,
            width_m: h,
        };
        let six_glue = build_six_panel_box(
            a,
            b,
            c,
            h,
            &sec,
            nx,
            ny,
            nx,
            ny,
            nz,
            TieKind::Compliant(glue),
        );
        let comps_glue: Vec<ShellComponent<'_>> = (0..6)
            .map(|i| ShellComponent {
                id: i,
                mesh: &six_glue.meshes[i],
                model: &six_glue.models[i],
                nodal_thickness_m: &six_glue.thicknesses[i],
            })
            .collect();
        let tied_glue = tie_shells(&comps_glue, &six_glue.ties, budget, || false).unwrap();
        let (lams, _) = solve_first_modes(&tied_glue.model, 14, 5.0e8);
        spectra.push(lams);
    }

    // 1. Rayleigh monotonicity (G0): every flexible eigenvalue (modes 6..14) increases monotonically with k.
    for m in 6..14 {
        for r in 0..(scales.len() - 1) {
            assert!(
                spectra[r + 1][m] > spectra[r][m],
                "mode {m} did not increase monotonically from rung {r} ({:.6e}) to rung {} ({:.6e})",
                spectra[r][m],
                r + 1,
                spectra[r + 1][m]
            );
        }
        assert!(
            lam_rigid[m] > spectra[scales.len() - 1][m],
            "rigid limit must bound compliant eigenvalue from above for mode {m}"
        );
    }

    // 2. Each flexible mode converges to the rigid-tie spectrum at slope -1 +- 0.2
    //    in 1/k, per x10 step. Gated per mode: a mean over modes could hide one
    //    mode converging at a different rate.
    let mut slopes = Vec::with_capacity(8);
    for m in 6..14 {
        let gaps: Vec<f64> = (0..scales.len())
            .map(|r| (lam_rigid[m] - spectra[r][m]).abs() / lam_rigid[m])
            .collect();
        for r in 0..(scales.len() - 1) {
            let slope = det::ln(gaps[r + 1] / gaps[r]) / det::ln(10.0);
            slopes.push(slope);
            assert!(
                (slope - (-1.0)).abs() <= 0.2,
                "mode {m}, rung {r}->{}: observed log-slope {slope:.4} outside -1.0 +/- 0.2 (gaps={gaps:?})",
                r + 1
            );
        }
    }
    println!(
        "G1 glue compliance rigid limit: per-mode slopes (modes 6..14, 3 steps each) in [{:.4}, {:.4}]",
        slopes.iter().copied().fold(f64::INFINITY, f64::min),
        slopes.iter().copied().fold(f64::NEG_INFINITY, f64::max)
    );
}

fn stiffness_per_area_dims() -> Dims {
    Pressure::DIMS
        .checked_minus(Length::DIMS)
        .expect("N/m^3 dimensions")
}

fn viscosity_per_area_dims() -> Dims {
    DynViscosity::DIMS
        .checked_minus(Length::DIMS)
        .expect("N*s/m^3 dimensions")
}

fn resolve_glue_card_fixture(
    card: &InterfaceSystemCard,
    expected_a: &MaterialStateId,
    expected_b: &MaterialStateId,
    point: &QueryPoint,
    width_m: f64,
) -> Result<(GlueLine, [fs_matdb::PropertyUsageReceipt; 6]), String> {
    if &card.surface_a().material != expected_a || &card.surface_b().material != expected_b {
        return Err("ordered surface mismatch".to_string());
    }
    let Some(third_body) = card.third_body() else {
        return Err("missing third_body adhesive".to_string());
    };
    if third_body.trim().is_empty() || card.history().trim().is_empty() {
        return Err("blank third_body or history".to_string());
    }
    let k_dims = stiffness_per_area_dims();
    let c_dims = viscosity_per_area_dims();
    let keys = [
        ("glue-normal-stiffness", k_dims),
        ("glue-shear-stiffness-1", k_dims),
        ("glue-shear-stiffness-2", k_dims),
        ("glue-normal-viscosity", c_dims),
        ("glue-shear-viscosity-1", c_dims),
        ("glue-shear-viscosity-2", c_dims),
    ];
    let mut vals = [0.0; 6];
    let mut receipts = Vec::with_capacity(6);
    for (i, (key, dims)) in keys.iter().enumerate() {
        let ans = card
            .claims()
            .query(key, point, SelectionPolicy::SingleClaimOnly)
            .map_err(|e| format!("missing or invalid property {key}: {e}"))?;
        if ans.evidence.value.dims != *dims {
            return Err(format!("unit mismatch on {key}"));
        }
        vals[i] = ans.evidence.value.value;
        receipts.push(ans.receipt);
    }
    let glue = GlueLine {
        normal_n_m3: vals[0],
        shear_n_m3: [vals[1], vals[2]],
        viscosity_n_s_m3: Some([vals[3], vals[4], vals[5]]),
        width_m,
    };
    glue.validate().map_err(|e| format!("{e:?}"))?;
    Ok((glue, receipts.try_into().unwrap()))
}

#[test]
fn g0_glue_compliance_is_data() {
    let spruce = MaterialStateId {
        chemistry: "Picea-abies-spruce".to_string(),
        phase: "cellular-orthotropic-wood".to_string(),
        process: "quarter-sawn-air-dried-12pct-mc".to_string(),
        revision: 1,
    };
    let maple = MaterialStateId {
        chemistry: "Acer-pseudoplatanus-maple".to_string(),
        phase: "cellular-orthotropic-wood".to_string(),
        process: "quarter-sawn-rib-bent-12pct-mc".to_string(),
        revision: 1,
    };
    let k_dims = stiffness_per_area_dims();
    let c_dims = viscosity_per_area_dims();

    let build_card = |omit: Option<&str>| -> InterfaceSystemCard {
        let mut claims = ClaimSet::new();
        for (k, v, d) in [
            ("glue-normal-stiffness", 2.0e11, k_dims),
            ("glue-shear-stiffness-1", 1.0e11, k_dims),
            ("glue-shear-stiffness-2", 1.0e11, k_dims),
            ("glue-normal-viscosity", 4.0e5, c_dims),
            ("glue-shear-viscosity-1", 2.0e5, c_dims),
            ("glue-shear-viscosity-2", 2.0e5, c_dims),
        ] {
            if omit == Some(k) {
                continue;
            }
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(k, d),
                    value: PropertyValue::Scalar { value: v, dims: d },
                    validity: ValidityDomain::unconstrained().with("temperature", 270.0, 330.0),
                    uncertainty: UncertaintyModel::RelativeHalfWidth {
                        fraction: 0.1,
                        confidence: 0.95,
                    },
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: Vec::new(),
                    provenance: Provenance {
                        source: "authored-fixture:hide-glue-spruce-maple-v1".to_string(),
                        license: "CC0-1.0".to_string(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
        InterfaceSystemCard::assemble(
            SurfaceSpec {
                material: spruce.clone(),
                texture_frame: "planed-toothed-0.1mm-L".to_string(),
            },
            SurfaceSpec {
                material: maple.clone(),
                texture_frame: "scraped-rib-edge-L".to_string(),
            },
            SystemContext {
                medium: "cured-bondline-0.08mm".to_string(),
                third_body: Some("authored-fixture:hot-hide-glue-192g-bloom".to_string()),
                environment: "air-293K-50pctRH".to_string(),
                history: "cured-72h-293K-virgin".to_string(),
            },
            claims,
            Vec::new(),
        )
        .unwrap()
    };

    let card = build_card(None);
    assert!(
        card.third_body().unwrap().starts_with("authored-fixture:"),
        "test fixture card must be named as an authored fixture"
    );
    let point = QueryPoint::new().with("temperature", 293.15).unwrap();
    let (glue, receipts) = resolve_glue_card_fixture(&card, &spruce, &maple, &point, 0.003)
        .expect("valid card resolves");
    assert_eq!(glue.normal_n_m3, 2.0e11);
    assert_eq!(receipts.len(), 6);

    // Swapping ordered surfaces (maple, spruce) refuses.
    assert!(resolve_glue_card_fixture(&card, &maple, &spruce, &point, 0.003).is_err());

    // Missing any required property refuses by name.
    let incomplete = build_card(Some("glue-shear-stiffness-2"));
    let err = resolve_glue_card_fixture(&incomplete, &spruce, &maple, &point, 0.003).unwrap_err();
    assert!(err.contains("glue-shear-stiffness-2"));
}

#[test]
fn g0_glue_loss_nonnegative_and_counted_once() {
    let (a, b, c, h) = (0.24, 0.18, 0.12, 0.003);
    let sec = PlateSection::isotropic(10.0e9, 0.30, h, 500.0).unwrap();
    let glue = GlueLine {
        normal_n_m3: 2.0e11,
        shear_n_m3: [1.0e11, 1.0e11],
        viscosity_n_s_m3: Some([4.0e5, 2.0e5, 2.0e5]),
        width_m: h,
    };
    let six = build_six_panel_box(a, b, c, h, &sec, 2, 2, 2, 2, 2, TieKind::Compliant(glue));
    let comps: Vec<ShellComponent<'_>> = (0..6)
        .map(|i| ShellComponent {
            id: i,
            mesh: &six.meshes[i],
            model: &six.models[i],
            nodal_thickness_m: &six.thicknesses[i],
        })
        .collect();
    let budget = TieBudget {
        max_segments: 256,
        max_slave_dofs: 2048,
        gap_tolerance_m: 1e-9,
        angle_tolerance_rad: 1e-6,
    };
    let tied: TiedShellModel = tie_shells(&comps, &six.ties, budget, || false).unwrap();
    let (_lams, phis) = solve_first_modes(&tied.model, 14, 5.0e8);

    let c_modal = tied.project_glue_damping(&phis).unwrap();
    assert_eq!(c_modal.len(), 14 * 14);

    // 1. Per-mode loss is non-negative for all modes, vanishes on the 6 rigid-body modes,
    //    and is strictly positive on flexible modes that deform the glue seams.
    for m in 0..6 {
        let diag = c_modal[m * 14 + m];
        assert!(
            diag.abs() <= 1e-9,
            "rigid mode {m} must have zero glue dissipation, got {diag:.3e}"
        );
    }
    for m in 6..14 {
        let diag = c_modal[m * 14 + m];
        assert!(
            diag > 0.0,
            "flexible mode {m} must have nonnegative/positive glue loss, got {diag:.3e}"
        );
        // Counted once: modal quadratic form phi_m^T C_glue phi_m equals glue_dissipation_rate(&phis[m]).
        let direct_rate = tied.glue_dissipation_rate(&phis[m]).unwrap();
        assert!(
            (diag - direct_rate).abs() <= 1e-12 * direct_rate.max(1.0),
            "modal projection must match single-owner glue_dissipation_rate without double counting"
        );
    }
}

fn build_box_skins(six: &SixPanelBox, sec: &PlateSection) -> Vec<ShellRadiationSurface> {
    let rad_budget = RadiationSurfaceBudget {
        max_panels: 1024,
        max_panel_modes: 4096,
    };
    let slice_opts = SliceOptions {
        ritz_tol: 1e-12,
        max_lanczos: 128,
        max_restarts: 32,
        mass_diagonal_equilibration: true,
        ..SliceOptions::default()
    };
    (0..6)
        .map(|i| {
            let report = modes_shell(&six.models[i], (10.0, 1.0e9), &slice_opts).unwrap();
            let red = ShellReduction::new(
                &six.meshes[i],
                &vec![*sec; six.meshes[i].tris.len()],
                &six.models[i],
                &report.modes[..1],
                ReductionBudget {
                    max_modes: 8,
                    max_facet_modes: 4096,
                    relative_tolerance: 1e-6,
                },
            )
            .unwrap();
            red.radiation_surface(&six.thicknesses[i], rad_budget)
                .unwrap()
        })
        .collect()
}

#[test]
fn g0_fluid_boundary_ownership_audit() {
    let (a, b, c, h) = (0.30, 0.22, 0.14, 0.003);
    let sec = PlateSection::isotropic(10.0e9, 0.30, h, 500.0).unwrap();
    let six = build_six_panel_box(a, b, c, h, &sec, 2, 2, 2, 2, 2, TieKind::Rigid);
    let skins = build_box_skins(&six, &sec);
    let comps: Vec<ShellComponent<'_>> = (0..6)
        .map(|i| ShellComponent {
            id: i,
            mesh: &six.meshes[i],
            model: &six.models[i],
            nodal_thickness_m: &six.thicknesses[i],
        })
        .collect();

    // All 6 panels have outward midsurface normals (+n = Exterior, -n = InteriorCavity(0)).
    let doms = vec![
        (
            AcousticDomainId::Exterior,
            AcousticDomainId::InteriorCavity(0)
        );
        6
    ];
    let assignments = FluidBoundaryAssignments::from_component_sides(
        &comps,
        &doms,
        six.ties.clone(),
        Vec::new(),
        Vec::new(),
    );

    let receipt = FluidBoundaryOwnership::audit(&comps, &skins, &assignments)
        .expect("valid closed box passes audit");
    let cav = &receipt.domains[&AcousticDomainId::InteriorCavity(0)];
    let ext = &receipt.domains[&AcousticDomainId::Exterior];
    assert_eq!(cav.faces.len(), 6 * 8);
    assert_eq!(ext.faces.len(), 6 * 8);
    assert!(receipt.solid_interior_area_m2 > 0.0);

    // Verify deterministic sorted (component, facet, side) order in receipt.
    assert!(cav.faces.windows(2).all(|w| w[0] < w[1]));
    assert!(ext.faces.windows(2).all(|w| w[0] < w[1]));

    // Also verify FaceTie glue-overlap faces are classified as solid-interior (owned by no fluid).
    let ftie = FaceTie {
        master_component: 0,
        master_tris: vec![0],
        slave_component: 1,
        slave_tris: vec![0],
        offset: TieOffset::Midsurface,
        rotations: RotationCoupling::Full,
        kind: TieKind::Rigid,
    };
    let assignments_with_ftie = FluidBoundaryAssignments::from_component_sides(
        &comps,
        &doms,
        six.ties.clone(),
        vec![ftie.clone()],
        Vec::new(),
    );
    let receipt_ftie =
        FluidBoundaryOwnership::audit(&comps, &skins, &assignments_with_ftie).unwrap();
    assert_eq!(receipt_ftie.solid_interior_faces.len(), 2);

    // Also verify an open box with an explicitly declared ApertureLoop passes when the aperture loop
    // covers the untied boundary edges!
    // Remove the 4 ties touching Top panel (1), and omit Top panel (1) from InteriorCavity(0) by
    // declaring the top rim loop on panels 2, 5, 3, 4 as open apertures:
    // Each mutation must refuse for its own reason, not for any reason.
    let refusal = |assignments: &FluidBoundaryAssignments| -> &'static str {
        match FluidBoundaryOwnership::audit(&comps, &skins, assignments).err() {
            Some(fs_plate::PlateError::BadConstraint { what }) => what,
            other => panic!("expected a BadConstraint refusal, got {other:?}"),
        }
    };
    let mut open_ap_assignments = assignments.clone();
    // Mutation 1: Missing owner refuses.
    open_ap_assignments.face_owners.pop();
    let why = refusal(&open_ap_assignments);
    assert!(why.starts_with("missing owner:"), "{why}");

    // Mutation 2: Double owner refuses.
    let mut double_assignments = assignments.clone();
    let first_entry = double_assignments.face_owners[0];
    double_assignments
        .face_owners
        .push((first_entry.0, AcousticDomainId::Exterior));
    let why = refusal(&double_assignments);
    assert!(why.starts_with("double owner:"), "{why}");

    // Mutation 3: Orientation conflict (flip one panel's Positive/Negative domain assignment) refuses.
    let mut flipped_doms = doms.clone();
    flipped_doms[3] = (
        AcousticDomainId::InteriorCavity(0),
        AcousticDomainId::Exterior,
    );
    let flipped_assignments = FluidBoundaryAssignments::from_component_sides(
        &comps,
        &flipped_doms,
        six.ties.clone(),
        Vec::new(),
        Vec::new(),
    );
    let why = refusal(&flipped_assignments);
    assert!(why.starts_with("orientation conflict:"), "{why}");

    // Mutation 4: Unclosed interior (remove one edge tie without declaring an ApertureLoop) refuses.
    let mut unclosed_ties = six.ties.clone();
    unclosed_ties.pop();
    let unclosed_assignments = FluidBoundaryAssignments::from_component_sides(
        &comps,
        &doms,
        unclosed_ties,
        Vec::new(),
        Vec::new(),
    );
    let why = refusal(&unclosed_assignments);
    assert!(why.starts_with("unclosed interior:"), "{why}");

    // Verify a 6-panel box with a square sound-hole in Top panel (1) and a declared closed
    // ApertureLoop on the 4 sound-hole nodes [5, 6, 10, 9] succeeds (and refuses if undeclared).
    let mut box_hole = build_six_panel_box(a, b, c, h, &sec, 3, 3, 3, 3, 2, TieKind::Rigid);
    // Remove the 2 center-quad triangles (indices 8 and 9: [5, 6, 10] and [5, 10, 9]) on Top panel 1.
    let mut top_tris = box_hole.meshes[1].tris.clone();
    top_tris.remove(9);
    top_tris.remove(8);
    box_hole.meshes[1] = ShellMesh::new(box_hole.meshes[1].nodes.clone(), top_tris).unwrap();
    box_hole.models[1] =
        assemble_shell(&box_hole.meshes[1], &sec, &[], ShellSupport::Free).unwrap();
    let hole_skins = build_box_skins(&box_hole, &sec);
    let hole_comps: Vec<ShellComponent<'_>> = (0..6)
        .map(|i| ShellComponent {
            id: i,
            mesh: &box_hole.meshes[i],
            model: &box_hole.models[i],
            nodal_thickness_m: &box_hole.thicknesses[i],
        })
        .collect();
    let ap_loop = ApertureLoop {
        cavity_id: 0,
        component: 1,
        nodes: vec![5, 6, 10, 9],
    };
    let hole_assignments = FluidBoundaryAssignments::from_component_sides(
        &hole_comps,
        &doms,
        box_hole.ties.clone(),
        Vec::new(),
        vec![ap_loop],
    );
    let hole_receipt = FluidBoundaryOwnership::audit(&hole_comps, &hole_skins, &hole_assignments)
        .expect("box with declared sound-hole ApertureLoop succeeds");
    assert_eq!(hole_receipt.open_aperture_loops, 1);
    let _ = SkinFaceId {
        component: 0,
        facet: 0,
        side: ShellFace::Positive,
    };
}

#[test]
fn g0_interior_volume_from_skins() {
    // Test multiple box aspect ratios and mesh densities (including non-matching 2:3 meshes!).
    for &(a, b, c, h, nx_tb, ny_tb, nx_fb, ny_lr, nz_side) in &[
        (0.30, 0.20, 0.12, 0.003, 2, 2, 2, 2, 2),
        (0.48, 0.36, 0.11, 0.0025, 3, 3, 3, 3, 2),
        (0.35, 0.25, 0.15, 0.004, 2, 2, 3, 3, 2),
    ] {
        let sec = PlateSection::isotropic(10.0e9, 0.30, h, 500.0).unwrap();
        let six = build_six_panel_box(
            a,
            b,
            c,
            h,
            &sec,
            nx_tb,
            ny_tb,
            nx_fb,
            ny_lr,
            nz_side,
            TieKind::Rigid,
        );
        let skins = build_box_skins(&six, &sec);
        let comps: Vec<ShellComponent<'_>> = (0..6)
            .map(|i| ShellComponent {
                id: i,
                mesh: &six.meshes[i],
                model: &six.models[i],
                nodal_thickness_m: &six.thicknesses[i],
            })
            .collect();
        let doms = vec![
            (
                AcousticDomainId::Exterior,
                AcousticDomainId::InteriorCavity(0)
            );
            6
        ];
        let assignments = FluidBoundaryAssignments::from_component_sides(
            &comps,
            &doms,
            six.ties.clone(),
            Vec::new(),
            Vec::new(),
        );
        let receipt = FluidBoundaryOwnership::audit(&comps, &skins, &assignments).unwrap();
        let cav = &receipt.domains[&AcousticDomainId::InteriorCavity(0)];
        let actual_vol = cav
            .enclosed_volume_m3
            .expect("interior cavity has enclosed volume");
        let expected_vol = (a - h) * (b - h) * (c - h);
        let rel_err = (actual_vol - expected_vol).abs() / expected_vol;
        assert!(
            rel_err <= 1e-12,
            "box ({a}x{b}x{c}, h={h}): enclosed volume {actual_vol:.15e} vs (a-h)(b-h)(c-h) = {expected_vol:.15e}, rel_err={rel_err:.3e} > 1e-12"
        );
        println!(
            "G0 interior volume from skins ({a}x{b}x{c}, h={h}): actual={actual_vol:.12e} m^3, expected={expected_vol:.12e} m^3, rel_err={rel_err:.3e} <= 1e-12"
        );
    }
}
