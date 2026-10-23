//! Integration and Gauntlet verification suite for `fs_plate::shell::junction` (`M1.2-W2a`).
//!
//! Covers both the CST/DKT shell path (`assemble_shell` / `assemble_shell_with_drilling`)
//! and the shear-deformable MITC3+ path (`assemble_mitc3p`) across:
//! - `g0_t_junction_rigid_modes_and_constant_strain_patch`
//! - `g0_tie_virtual_work_identity`
//! - `g3_matching_tie_reproduces_monolithic_pencil`
//! - `g1_nonmatching_mortar_converges_to_monolithic`
//! - `g1_skin_offset_t_junction`
//! - `g0_drilling_alpha_insensitivity_at_junction`

use fs_math::det;
use fs_mms::{LadderSide, ORDER_GATE_TOLERANCE, OrderGate, RefinementLadder, fit_order};
use fs_modal::{SliceOptions, slice_window};
use fs_plate::shell::junction::{
    EdgeTie, RotationCoupling, ShellComponent, TieBudget, TieKind, TieOffset, TieSide,
    TiedShellModel, tie_shells,
};
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, NodalDirectors, ShearShellSection, ShellConstraint, ShellMass,
    assemble_mitc3p,
};
use fs_plate::shell::{
    DRILLING_ALPHA, ShellMesh, ShellModel, ShellSupport, assemble_shell_with_drilling, modes_shell,
};
use fs_plate::{PlateError, PlateSection};
use fs_rand::StreamKey;
use fs_sparse::{Coo, Csr, DirectOrdering, LdltOptions, SymbolicLdlt};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ElementKind {
    CstDkt,
    Mitc3p,
}

impl ElementKind {
    const fn name(self) -> &'static str {
        match self {
            Self::CstDkt => "cst_dkt",
            Self::Mitc3p => "mitc3+",
        }
    }
}

fn default_tie_budget() -> TieBudget {
    TieBudget {
        max_segments: 4096,
        max_slave_dofs: 16384,
        gap_tolerance_m: 1e-6,
        angle_tolerance_rad: 1e-4,
    }
}

/// Build a structured rectangular `ShellMesh` in the `x-y` plane at `z = z0`.
fn rect_xy_mesh(x0: f64, x1: f64, y0: f64, y1: f64, z0: f64, nx: usize, ny: usize) -> ShellMesh {
    let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1));
    for j in 0..=ny {
        let y = y0 + (y1 - y0) * (j as f64 / ny as f64);
        for i in 0..=nx {
            let x = x0 + (x1 - x0) * (i as f64 / nx as f64);
            nodes.push([x, y, z0]);
        }
    }
    let idx = |i: usize, j: usize| j * (nx + 1) + i;
    let mut tris = Vec::with_capacity(2 * nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            tris.push([idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)]);
            tris.push([idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)]);
        }
    }
    ShellMesh::new(nodes, tris).expect("valid xy rectangular shell mesh")
}

/// Build a structured rectangular `ShellMesh` in the `x-z` plane at `y = y0`.
fn rect_xz_mesh(x0: f64, x1: f64, y0: f64, z0: f64, z1: f64, nx: usize, nz: usize) -> ShellMesh {
    let mut nodes = Vec::with_capacity((nx + 1) * (nz + 1));
    for k in 0..=nz {
        let z = z0 + (z1 - z0) * (k as f64 / nz as f64);
        for i in 0..=nx {
            let x = x0 + (x1 - x0) * (i as f64 / nx as f64);
            nodes.push([x, y0, z]);
        }
    }
    let idx = |i: usize, k: usize| k * (nx + 1) + i;
    let mut tris = Vec::with_capacity(2 * nx * nz);
    for k in 0..nz {
        for i in 0..nx {
            tris.push([idx(i, k), idx(i + 1, k), idx(i + 1, k + 1)]);
            tris.push([idx(i, k), idx(i + 1, k + 1), idx(i, k + 1)]);
        }
    }
    ShellMesh::new(nodes, tris).expect("valid xz rectangular shell mesh")
}

/// Assemble a single component using either CST/DKT or MITC3+ with uniform isotropic properties.
fn assemble_component_model(
    kind: ElementKind,
    mesh: &ShellMesh,
    e_pa: f64,
    nu: f64,
    h_m: f64,
    rho: f64,
    clamped_nodes: &[usize],
    drilling_alpha: f64,
) -> ShellModel {
    match kind {
        ElementKind::CstDkt => {
            let sec = PlateSection::isotropic(e_pa, nu, h_m, rho).expect("valid plate section");
            let support = if clamped_nodes.is_empty() {
                ShellSupport::Free
            } else {
                ShellSupport::Clamped
            };
            assemble_shell_with_drilling(mesh, &sec, clamped_nodes, support, drilling_alpha)
                .expect("valid CST/DKT assembly")
        }
        ElementKind::Mitc3p => {
            let g_pa = e_pa / (2.0 * (1.0 + nu));
            let sec = ShearShellSection::from_orthotropic_solid(
                e_pa,
                e_pa,
                nu,
                g_pa,
                g_pa,
                g_pa,
                h_m,
                rho,
                0.0,
                5.0 / 6.0,
            )
            .expect("valid MITC3+ section");
            let sections = vec![sec; mesh.tris.len()];
            // Compute uniform normal from first facet (all components here are planar facets).
            let n0 = mesh.facet(0).expect("valid facet").frame[2];
            let directors = NodalDirectors::new(mesh, vec![n0; mesh.nodes.len()], 0.9)
                .expect("valid directors");
            let constraints: Vec<ShellConstraint> = clamped_nodes
                .iter()
                .map(|&node| ShellConstraint {
                    node,
                    mask: [true; 6],
                    frame: ConstraintFrame::Global,
                })
                .collect();
            assemble_mitc3p(
                mesh,
                &directors,
                &sections,
                &constraints,
                ShellMass::Lumped,
                MitcBudget {
                    max_nodes: mesh.nodes.len(),
                    max_triangles: mesh.tris.len(),
                },
                drilling_alpha,
                || false,
            )
            .expect("valid MITC3+ assembly")
            .pencil()
            .clone()
        }
    }
}

/// Shift `K_sigma = K + sigma * M` so rigid modes sit at `+sigma` and can be extracted via `slice_window`.
fn add_scaled_csr(k: &Csr, m: &Csr, sigma: f64) -> Csr {
    let n = k.nrows();
    let mut coo = Coo::new(n, n);
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

fn find_upper_window_bound(k: &Csr, m: &Csr, min_count: usize, start_hi: f64) -> f64 {
    find_upper_window_bound_by(k, m, min_count, start_hi, 2.0)
}

/// Grow `start_hi` by `growth` until Sylvester inertia counts `min_count`
/// eigenvalues below it.
fn find_upper_window_bound_by(
    k: &Csr,
    m: &Csr,
    min_count: usize,
    start_hi: f64,
    growth: f64,
) -> f64 {
    let mut hi = start_hi;
    for _ in 0..40 {
        let pencil = add_scaled_csr(k, m, -hi);
        let symbolic = SymbolicLdlt::analyze(&pencil, DirectOrdering::Amd).unwrap();
        let neg = symbolic
            .factor(&pencil, &LdltOptions::default())
            .unwrap()
            .inertia()
            .negative;
        if neg >= min_count {
            return hi;
        }
        hi *= growth;
    }
    panic!("spectral window upper bound search exhausted for count {min_count}");
}

/// Collect the lowest `count` eigenvalues of `(model.k, model.m)` using a positive shift `sigma`.
fn lowest_eigenvalues_shifted(model: &ShellModel, count: usize, sigma: f64) -> Vec<f64> {
    let k_shifted = add_scaled_csr(&model.k, &model.m, sigma);
    let hi = find_upper_window_bound(&k_shifted, &model.m, count, 2.0 * sigma);
    let opts = SliceOptions {
        ritz_tol: 1e-14,
        max_lanczos: 1024,
        max_restarts: 64,
        mass_diagonal_equilibration: true,
        ..SliceOptions::default()
    };
    let report = slice_window(&k_shifted, &model.m, (0.25 * sigma, hi), &opts)
        .expect("shifted modal slice succeeds");
    assert!(
        report.modes.len() >= count,
        "expected at least {count} eigenvalues in window, got {}",
        report.modes.len()
    );
    report
        .modes
        .iter()
        .take(count)
        .map(|mode| mode.lambda - sigma)
        .collect()
}

/// Collect the lowest `count` strictly positive eigenvalues of a constrained model `(model.k, model.m)`.
fn lowest_positive_eigenvalues(model: &ShellModel, count: usize) -> Vec<f64> {
    positive_eigenvalues_from(model, count, 100.0, 2.0)
}

/// As [`lowest_positive_eigenvalues`], with the inertia search for the window
/// starting at `start_hi` and growing by `growth`.
fn positive_eigenvalues_from(
    model: &ShellModel,
    count: usize,
    start_hi: f64,
    growth: f64,
) -> Vec<f64> {
    let hi = find_upper_window_bound_by(&model.k, &model.m, count, start_hi, growth);
    let opts = SliceOptions {
        ritz_tol: 1e-14,
        max_lanczos: 1024,
        max_restarts: 64,
        mass_diagonal_equilibration: true,
        ..SliceOptions::default()
    };
    let report = modes_shell(model, (0.0, hi), &opts).expect("positive modal slice succeeds");
    assert!(
        report.modes.len() >= count,
        "expected at least {count} positive eigenvalues, got {}",
        report.modes.len()
    );
    report
        .modes
        .into_iter()
        .take(count)
        .map(|mode| mode.lambda)
        .collect()
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Reduced coordinates of a nodal field given on every component.
fn reduced_from_field(
    tied: &TiedShellModel,
    meshes: &[&ShellMesh],
    field: impl Fn([f64; 3]) -> [f64; 6],
) -> Vec<f64> {
    let mut u_r = vec![0.0; tied.model.free];
    for (c, mesh) in meshes.iter().enumerate() {
        for (n, &p) in mesh.nodes.iter().enumerate() {
            let value = field(p);
            for d in 0..6 {
                if let Some(r) = tied.map[c][n][d] {
                    u_r[r] = value[d];
                }
            }
        }
    }
    u_r
}

/// Largest deviation of the expanded components from `field`.
fn expansion_error(
    tied: &TiedShellModel,
    meshes: &[&ShellMesh],
    u_r: &[f64],
    field: impl Fn([f64; 3]) -> [f64; 6],
) -> f64 {
    let full = tied.expand_components(u_r).expect("expansion succeeds");
    let mut err = 0.0_f64;
    for (c, mesh) in meshes.iter().enumerate() {
        for (n, &p) in mesh.nodes.iter().enumerate() {
            let value = field(p);
            for d in 0..6 {
                err = err.max((full[c][6 * n + d] - value[d]).abs());
            }
        }
    }
    err
}

#[test]
fn g0_t_junction_rigid_modes_and_constant_strain_patch() {
    let started = std::time::Instant::now();
    let e_pa = 10e9;
    let nu = 0.25;
    let h_m = 0.012;
    let rho = 450.0;
    let lx = 0.30;
    let ly = 0.20;
    let hz = 0.12;

    // Three components meet at the seam y = z = 0: two flange halves and a web.
    // The left half is the master with 6 seam nodes; the right half and the web
    // are slaves with 4, a 2:3 node ratio on both mortar edges.
    let (nx_master, nx_slave, nz) = (5usize, 3usize, 2usize);
    let flange_left = rect_xy_mesh(0.0, lx, -0.5 * ly, 0.0, 0.0, nx_master, 1);
    let flange_right = rect_xy_mesh(0.0, lx, 0.0, 0.5 * ly, 0.0, nx_slave, 1);
    let web = rect_xz_mesh(0.0, lx, 0.0, -hz, 0.0, nx_slave, nz);
    let meshes = [&flange_left, &flange_right, &web];
    let master_seam: Vec<usize> = (0..=nx_master).map(|i| (nx_master + 1) + i).collect();
    let right_seam: Vec<usize> = (0..=nx_slave).collect();
    let web_seam: Vec<usize> = (0..=nx_slave).map(|i| nz * (nx_slave + 1) + i).collect();

    for (kind, rotations) in [ElementKind::CstDkt, ElementKind::Mitc3p]
        .into_iter()
        .flat_map(|kind| {
            [RotationCoupling::Full, RotationCoupling::SeamTangent]
                .map(|rotations| (kind, rotations))
        })
    {
        let models: Vec<ShellModel> = meshes
            .iter()
            .map(|mesh| {
                assemble_component_model(kind, mesh, e_pa, nu, h_m, rho, &[], DRILLING_ALPHA)
            })
            .collect();
        let thickness: Vec<Vec<f64>> = meshes.iter().map(|m| vec![h_m; m.nodes.len()]).collect();
        let comps: Vec<ShellComponent> = (0..3)
            .map(|c| ShellComponent {
                id: c,
                mesh: meshes[c],
                model: &models[c],
                nodal_thickness_m: &thickness[c],
            })
            .collect();
        let tie = |component: usize, nodes: &[usize], rotations: RotationCoupling| EdgeTie {
            a: TieSide {
                component: 0,
                nodes: master_seam.clone(),
            },
            b: TieSide {
                component,
                nodes: nodes.to_vec(),
            },
            offset: TieOffset::Midsurface,
            rotations,
            kind: TieKind::Rigid,
        };
        // The coplanar half shares all six rotations; the web is the fold.
        let ties = [
            tie(1, &right_seam, RotationCoupling::Full),
            tie(2, &web_seam, rotations),
        ];
        let tied = tie_shells(&comps, &ties, default_tie_budget(), || false)
            .expect("T-junction tie succeeds");

        // 1. Exactly six rigid-body modes (|lambda_i| / lambda_7 < 1e-12).
        let eigs = lowest_eigenvalues_shifted(&tied.model, 8, 1.0);
        let lambda_7 = eigs[6];
        assert!(
            lambda_7 > 100.0,
            "7th mode must be positive elastic: {lambda_7}"
        );
        let mut max_rigid_ratio = 0.0_f64;
        for (i, &lam) in eigs.iter().take(6).enumerate() {
            let ratio = lam.abs() / lambda_7;
            max_rigid_ratio = max_rigid_ratio.max(ratio);
            assert!(
                ratio < 1e-12,
                "{} {rotations:?}: rigid mode {i} ratio |lambda_{i}|/lambda_7 = {ratio:e} >= 1e-12 (lambda_{i}={lam:e}, lambda_7={lambda_7:e})",
                kind.name()
            );
        }

        // 2. Tie forces under each of the six rigid motions: every component's
        // internal force K_c q_c, row by row relative to its absolute sum, and
        // the expanded slave field against the rigid field.
        let mut max_force_ratio = 0.0_f64;
        let mut max_rigid_field_err = 0.0_f64;
        for motion in 0..6 {
            let field = |p: [f64; 3]| -> [f64; 6] {
                let mut t = [0.0; 3];
                let mut w = [0.0; 3];
                if motion < 3 {
                    t[motion] = 1.0;
                } else {
                    w[motion - 3] = 1.0 / lx;
                }
                let r = cross3(w, p);
                [t[0] + r[0], t[1] + r[1], t[2] + r[2], w[0], w[1], w[2]]
            };
            let u_r = reduced_from_field(&tied, &meshes, field);
            max_rigid_field_err =
                max_rigid_field_err.max(expansion_error(&tied, &meshes, &u_r, field));
            let free = tied
                .expand_component_free(&u_r)
                .expect("expansion succeeds");
            for (c, model) in models.iter().enumerate() {
                let mut force = vec![0.0; model.free];
                model.k.spmv(&free[c], &mut force);
                for (i, &f) in force.iter().enumerate() {
                    let (cols, vals) = model.k.row(i);
                    let scale: f64 = cols
                        .iter()
                        .zip(vals)
                        .map(|(&j, &v)| (v * free[c][j]).abs())
                        .sum();
                    if scale > 0.0 {
                        max_force_ratio = max_force_ratio.max(f.abs() / scale);
                    }
                }
            }
        }
        assert!(
            max_force_ratio <= 1e-12 && max_rigid_field_err <= 1e-12,
            "{} {rotations:?}: rigid-motion tie force ratio {max_force_ratio:e}, field error {max_rigid_field_err:e}",
            kind.name()
        );

        // 3. Constant membrane strain across both 2:3 mortar edges.
        let eps_0 = 1e-4;
        let patch_err =
            verify_t_junction_constant_strain_patch(&meshes, &tied, e_pa, nu, h_m, lx, eps_0);
        assert!(
            patch_err < 1e-12,
            "{} {rotations:?}: T-junction constant-strain patch relative error {patch_err:e} >= 1e-12",
            kind.name()
        );

        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W2a\",\"test\":\"g0_t_junction_rigid_modes_and_constant_strain_patch\",\"tier\":\"G0\",\"element\":\"{}\",\"rotations\":\"{rotations:?}\",\"components\":3,\"seam_nodes\":[{},{},{}],\"max_rigid_ratio\":{:e},\"lambda_7\":{},\"rigid_tie_force_ratio\":{:e},\"rigid_field_err\":{:e},\"constant_strain_patch_err\":{:e},\"elapsed_s\":{}}}",
            kind.name(),
            master_seam.len(),
            right_seam.len(),
            web_seam.len(),
            max_rigid_ratio,
            lambda_7,
            max_force_ratio,
            max_rigid_field_err,
            patch_err,
            started.elapsed().as_secs_f64()
        );
    }
}

/// Uniaxial stress sigma_xx = E eps_0 on every component: u = eps_0 (x, -nu y,
/// -nu z), no rotation, end tractions on x = 0 and x = lx. The exact field must
/// balance the projected loads at every reduced DOF, seam DOFs included, and
/// the expanded slaves must reproduce it.
#[allow(clippy::too_many_arguments)]
fn verify_t_junction_constant_strain_patch(
    meshes: &[&ShellMesh],
    tied: &TiedShellModel,
    e_pa: f64,
    nu: f64,
    h_m: f64,
    lx: f64,
    eps_0: f64,
) -> f64 {
    let field = |p: [f64; 3]| {
        [
            eps_0 * p[0],
            -nu * eps_0 * p[1],
            -nu * eps_0 * p[2],
            0.0,
            0.0,
            0.0,
        ]
    };
    let u_r_exact = reduced_from_field(tied, meshes, field);

    // Consistent end loads: sigma_xx h over each end node's tributary length.
    // Every end edge here runs along y (flange, z = 0) or z (web, y = 0), so
    // y + z is the coordinate along it.
    let sigma_xx = e_pa * eps_0;
    let mut edge_length = 0.0;
    let loads: Vec<Vec<f64>> = meshes
        .iter()
        .map(|mesh| {
            let mut f = vec![0.0; 6 * mesh.nodes.len()];
            for (end, sign) in [(0.0, -1.0), (lx, 1.0)] {
                let mut edge: Vec<(f64, usize)> = mesh
                    .nodes
                    .iter()
                    .enumerate()
                    .filter(|(_, p)| p[0] == end)
                    .map(|(n, p)| (p[1] + p[2], n))
                    .collect();
                edge.sort_by(|a, b| a.0.total_cmp(&b.0));
                if sign > 0.0 {
                    edge_length += edge.last().unwrap().0 - edge[0].0;
                }
                for i in 0..edge.len() {
                    let below = if i > 0 {
                        edge[i].0 - edge[i - 1].0
                    } else {
                        0.0
                    };
                    let above = if i + 1 < edge.len() {
                        edge[i + 1].0 - edge[i].0
                    } else {
                        0.0
                    };
                    f[6 * edge[i].1] += sign * sigma_xx * h_m * 0.5 * (below + above);
                }
            }
            f
        })
        .collect();
    let load_refs: Vec<&[f64]> = loads.iter().map(Vec::as_slice).collect();
    let f_r = tied
        .project_loads(&load_refs)
        .expect("load projection succeeds");

    let mut r_vec = vec![0.0; tied.model.free];
    tied.model.k.spmv(&u_r_exact, &mut r_vec);
    let force_scale = sigma_xx * h_m * edge_length;
    let mut max_rel_res = 0.0_f64;
    for i in 0..tied.model.free {
        max_rel_res = max_rel_res.max((r_vec[i] - f_r[i]).abs() / force_scale);
    }
    max_rel_res.max(expansion_error(tied, meshes, &u_r_exact, field) / (eps_0 * lx))
}

#[test]
fn g0_tie_virtual_work_identity() {
    let started = std::time::Instant::now();
    let e_pa = 12e9;
    let nu = 0.3;
    let h_m = 0.004;
    let rho = 500.0;

    // Non-matching 3:2 T-junction with skin offset so both mortar weights and eccentric cross-product
    // coupling S(d) are exercised in v_s = T_{sm} v_m.
    let flange_mesh = rect_xy_mesh(0.0, 0.30, -0.10, 0.10, 0.0, 3, 2);
    let web_mesh = rect_xz_mesh(0.0, 0.30, 0.0, -0.12, -0.5 * h_m, 2, 2);

    let master_seam: Vec<usize> = (0..=3).map(|i| 1 * 4 + i).collect();
    let slave_seam: Vec<usize> = (0..=2).map(|i| 2 * 3 + i).collect();

    for (kind, rotations) in [ElementKind::CstDkt, ElementKind::Mitc3p]
        .into_iter()
        .flat_map(|kind| {
            [RotationCoupling::Full, RotationCoupling::SeamTangent]
                .map(|rotations| (kind, rotations))
        })
    {
        let mod0 =
            assemble_component_model(kind, &flange_mesh, e_pa, nu, h_m, rho, &[], DRILLING_ALPHA);
        let mod1 =
            assemble_component_model(kind, &web_mesh, e_pa, nu, h_m, rho, &[], DRILLING_ALPHA);
        let h0 = vec![h_m; flange_mesh.nodes.len()];
        let h1 = vec![h_m; web_mesh.nodes.len()];
        let comps = [
            ShellComponent {
                id: 0,
                mesh: &flange_mesh,
                model: &mod0,
                nodal_thickness_m: &h0,
            },
            ShellComponent {
                id: 1,
                mesh: &web_mesh,
                model: &mod1,
                nodal_thickness_m: &h1,
            },
        ];
        let ties = [EdgeTie {
            a: TieSide {
                component: 0,
                nodes: master_seam.clone(),
            },
            b: TieSide {
                component: 1,
                nodes: slave_seam.clone(),
            },
            offset: TieOffset::Skin { side: -1 },
            rotations,
            kind: TieKind::Rigid,
        }];
        let tied = tie_shells(&comps, &ties, default_tie_budget(), || false)
            .expect("non-matching skin-offset tie succeeds");

        let mut rng = StreamKey {
            seed: 0x2026_1020_0050,
            kernel: 1,
            tile: 0,
        }
        .stream();
        let mut max_rel_power_err = 0.0_f64;
        let mut min_cancellation = f64::INFINITY;
        for _trial in 0..1000 {
            let v_r: Vec<f64> = (0..tied.model.free)
                .map(|_| 2.0 * rng.next_f64() - 1.0)
                .collect();
            let lambda: Vec<f64> = (0..tied.receipt.slave_dofs)
                .map(|_| 2.0 * rng.next_f64() - 1.0)
                .collect();
            let (p_m, p_s) = tied
                .interface_virtual_power(&v_r, &lambda)
                .expect("virtual power evaluation succeeds");
            // Relative to the absolute work exchanged through the couplings,
            // sum_k |lambda_k| sum_r |w_kr v_r|, which bounds |P_m| and |P_s|.
            // A random draw whose P cancels to a small sum has its roundoff
            // set by this scale, not by |P|.
            let scale: f64 = tied
                .slave_couplings()
                .iter()
                .zip(&lambda)
                .map(|(coupling, &lam)| {
                    lam.abs()
                        * coupling
                            .master_terms
                            .iter()
                            .map(|&(r, w)| (w * v_r[r]).abs())
                            .sum::<f64>()
                })
                .sum();
            assert!(scale > 0.0);
            let rel_err = (p_m + p_s).abs() / scale;
            max_rel_power_err = max_rel_power_err.max(rel_err);
            min_cancellation = min_cancellation.min(p_m.abs() / scale);
            assert!(
                rel_err <= 1e-13,
                "{} {rotations:?}: virtual work identity violated: P_m={p_m}, P_s={p_s}, rel_err={rel_err:e}",
                kind.name()
            );
        }

        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W2a\",\"test\":\"g0_tie_virtual_work_identity\",\"tier\":\"G0\",\"element\":\"{}\",\"rotations\":\"{rotations:?}\",\"trials\":1000,\"max_rel_power_err\":{:e},\"min_abs_power_over_scale\":{min_cancellation:e},\"elapsed_s\":{}}}",
            kind.name(),
            max_rel_power_err,
            started.elapsed().as_secs_f64()
        );
    }
}

/// One conforming mesh from parts whose shared nodes coincide bit for bit,
/// with each part's node map into it.
fn weld(parts: &[&ShellMesh]) -> (ShellMesh, Vec<Vec<usize>>) {
    let mut index = BTreeMap::new();
    let mut nodes = Vec::new();
    let mut tris = Vec::new();
    let mut maps = Vec::with_capacity(parts.len());
    for part in parts {
        let map: Vec<usize> = part
            .nodes
            .iter()
            .map(|p| {
                *index.entry(p.map(f64::to_bits)).or_insert_with(|| {
                    nodes.push(*p);
                    nodes.len() - 1
                })
            })
            .collect();
        tris.extend(part.tris.iter().map(|t| t.map(|n| map[n])));
        maps.push(map);
    }
    (ShellMesh::new(nodes, tris).expect("welded mesh"), maps)
}

/// Reduced coordinate -> welded free index, through each component node's
/// place in the weld. Refuses unless it is a bijection.
fn reduced_to_welded(
    tied: &TiedShellModel,
    mono: &ShellModel,
    node_maps: &[Vec<usize>],
) -> Vec<usize> {
    let n = mono.free;
    assert_eq!(tied.model.free, n, "tied and welded free DOF counts");
    let mut to_welded: Vec<Option<usize>> = vec![None; n];
    for (c, map) in node_maps.iter().enumerate() {
        for (node, &w) in map.iter().enumerate() {
            for d in 0..6 {
                if let Some(r) = tied.map[c][node][d] {
                    let q = mono.dof_map[6 * w + d]
                        .expect("a retained tied DOF is free in the welded mesh");
                    if let Some(p) = to_welded[r] {
                        assert_eq!(p, q, "reduced DOF {r} reaches two welded DOFs");
                    }
                    to_welded[r] = Some(q);
                }
            }
        }
    }
    let to_welded: Vec<usize> = to_welded
        .into_iter()
        .map(|q| q.expect("every reduced DOF reaches the weld"))
        .collect();
    let mut reached = vec![false; n];
    for &q in &to_welded {
        assert!(!reached[q], "welded DOF {q} reached twice");
        reached[q] = true;
    }
    to_welded
}

/// Entry gaps between a tied matrix and the welded one under `to_welded`,
/// over both sparsity patterns, keyed by welded (row, column). Each gap is
/// relative to sqrt(D_i D_j) for the welded diagonal D: for a semidefinite
/// matrix that bounds |A_ij| and the sum of the element terms' magnitudes
/// assembled into it, so a reassociated sum differs from it by a few ulp on
/// this scale.
fn entry_gaps(
    reduced: &Csr,
    welded: &Csr,
    to_welded: &[usize],
    diagonal: &[f64],
) -> BTreeMap<(usize, usize), f64> {
    let n = to_welded.len();
    let mut to_tied = vec![0; n];
    for (r, &q) in to_welded.iter().enumerate() {
        to_tied[q] = r;
    }
    let scale = |i: usize, j: usize| det::sqrt(diagonal[i] * diagonal[j]);
    let mut gaps = BTreeMap::new();
    for r in 0..n {
        let (cols, vals) = reduced.row(r);
        for (&c, &v) in cols.iter().zip(vals) {
            let (i, j) = (to_welded[r], to_welded[c]);
            gaps.insert((i, j), (v - welded.get(i, j)).abs() / scale(i, j));
        }
    }
    for i in 0..n {
        let (cols, vals) = welded.row(i);
        for (&j, &v) in cols.iter().zip(vals) {
            gaps.entry((i, j))
                .or_insert_with(|| (reduced.get(to_tied[i], to_tied[j]) - v).abs() / scale(i, j));
        }
    }
    gaps
}

fn largest(gaps: impl Iterator<Item = f64>) -> f64 {
    gaps.fold(0.0_f64, f64::max)
}

/// `wa * a + wb * b` on the union pattern.
fn combine(a: &Csr, wa: f64, b: &Csr, wb: f64) -> Csr {
    let n = a.nrows();
    let mut coo = Coo::new(n, n);
    for r in 0..n {
        for (matrix, weight) in [(a, wa), (b, wb)] {
            let (cols, vals) = matrix.row(r);
            for (&c, &v) in cols.iter().zip(vals) {
                coo.push(r, c, weight * v);
            }
        }
    }
    coo.assemble()
}

/// Positive diagonal of a welded matrix, the scale of its entry gaps.
fn positive_diagonal(a: &Csr) -> Vec<f64> {
    let diagonal: Vec<f64> = (0..a.nrows()).map(|i| a.get(i, i)).collect();
    assert!(
        diagonal.iter().all(|&d| d > 0.0),
        "welded diagonal is positive"
    );
    diagonal
}

/// Bracket lambda_max of (K, M) for a diagonal M, on the symmetric form
/// A = M^-1/2 K M^-1/2. Power steps from a fixed fs-rand start give a
/// Rayleigh quotient, which for a semidefinite A rises toward lambda_max from
/// below; Gershgorin's absolute row sums of A bound it from above.
fn lambda_max_bracket(model: &ShellModel) -> (f64, f64) {
    let n = model.free;
    let inv_sqrt_m: Vec<f64> = (0..n)
        .map(|r| {
            let (cols, vals) = model.m.row(r);
            for (&c, &v) in cols.iter().zip(vals) {
                assert!(c == r || v == 0.0, "lumped mass is diagonal");
            }
            let m_rr = model.m.get(r, r);
            assert!(m_rr > 0.0, "lumped mass is positive");
            1.0 / det::sqrt(m_rr)
        })
        .collect();
    let gershgorin = (0..n)
        .map(|r| {
            let (cols, vals) = model.k.row(r);
            cols.iter()
                .zip(vals)
                .map(|(&c, &v)| (v * inv_sqrt_m[r] * inv_sqrt_m[c]).abs())
                .sum::<f64>()
        })
        .fold(0.0_f64, f64::max);
    let mut rng = StreamKey {
        seed: 0x2026_1023_0050,
        kernel: 3,
        tile: 0,
    }
    .stream();
    let mut x: Vec<f64> = (0..n).map(|_| 2.0 * rng.next_f64() - 1.0).collect();
    let mut y = vec![0.0; n];
    let mut ky = vec![0.0; n];
    let mut rayleigh = 0.0;
    for _ in 0..500 {
        let norm = det::sqrt(x.iter().map(|v| v * v).sum::<f64>());
        for i in 0..n {
            y[i] = x[i] / norm * inv_sqrt_m[i];
        }
        model.k.spmv(&y, &mut ky);
        // With x normalised, x^T A x = y^T K y, and A x = M^-1/2 K y.
        rayleigh = y.iter().zip(&ky).map(|(a, b)| a * b).sum::<f64>();
        for i in 0..n {
            x[i] = ky[i] * inv_sqrt_m[i];
        }
    }
    (rayleigh, gershgorin)
}

/// Entry gap gate on the tied pencil, in units of f64::EPSILON on the
/// sqrt(D_i D_j) scale. The welded and tied entries sum the same element
/// terms (at most nine at a T seam node) in a different association, through
/// tie weights that are 1 to roundoff.
const G3_ENTRY_ULPS: f64 = 8.0;
/// The eigenvalue gate is C * eps * lambda_max / lambda_1 with this C, the
/// floor of two independent backward-stable solves of the same pencil. On
/// these coarse fixtures lambda_max is set by in-plane and drilling stiffness,
/// so the floor sits orders above the measured gaps; the entrywise identity
/// is the binding check.
const G3_FLOOR_C: f64 = 10.0;
/// Elastic eigenvalues compared, after the six rigid modes.
const G3_ELASTIC_MODES: usize = 20;

/// The tied pencil against the welded monolithic one.
struct PencilComparison {
    free: usize,
    seam_fan_dofs: usize,
    stiffness_entry_gap: f64,
    physical_entry_gap: f64,
    gauge_entry_gap: f64,
    gauge_gap_off_seam_fans: f64,
    mass_entry_gap: f64,
    lambda_1: f64,
    lambda_max: f64,
    lambda_max_gershgorin: f64,
    eigenvalue_floor: f64,
    elastic_gaps: Vec<f64>,
}

/// Tie `parts` rigidly (`Full`, from part 0's seam to each other part's seam),
/// weld the same parts into one conforming mesh, and compare the two pencils.
fn tied_vs_welded(
    kind: ElementKind,
    parts: &[&ShellMesh],
    seams: &[Vec<usize>],
    (e_pa, nu, h_m, rho): (f64, f64, f64, f64),
) -> PencilComparison {
    let (welded, node_maps) = weld(parts);
    let pencils = |alpha: f64| -> (TiedShellModel, ShellModel) {
        let models: Vec<ShellModel> = parts
            .iter()
            .map(|mesh| assemble_component_model(kind, mesh, e_pa, nu, h_m, rho, &[], alpha))
            .collect();
        let thickness: Vec<Vec<f64>> = parts.iter().map(|m| vec![h_m; m.nodes.len()]).collect();
        let comps: Vec<ShellComponent> = (0..parts.len())
            .map(|c| ShellComponent {
                id: c,
                mesh: parts[c],
                model: &models[c],
                nodal_thickness_m: &thickness[c],
            })
            .collect();
        let ties: Vec<EdgeTie> = (1..parts.len())
            .map(|c| EdgeTie {
                a: TieSide {
                    component: 0,
                    nodes: seams[0].clone(),
                },
                b: TieSide {
                    component: c,
                    nodes: seams[c].clone(),
                },
                offset: TieOffset::Midsurface,
                rotations: RotationCoupling::Full,
                kind: TieKind::Rigid,
            })
            .collect();
        let tied = tie_shells(&comps, &ties, default_tie_budget(), || false)
            .expect("matching tie succeeds");
        let mono = assemble_component_model(kind, &welded, e_pa, nu, h_m, rho, &[], alpha);
        (tied, mono)
    };
    let (tied, mono) = pencils(DRILLING_ALPHA);
    let (tied_2, mono_2) = pencils(2.0 * DRILLING_ALPHA);
    let to_welded = reduced_to_welded(&tied, &mono, &node_maps);
    assert_eq!(to_welded, reduced_to_welded(&tied_2, &mono_2, &node_maps));
    let k_diagonal = positive_diagonal(&mono.k);
    let m_diagonal = positive_diagonal(&mono.m);

    // Both elements' K is linear in the drilling strength alpha, K = P +
    // alpha G, so 2 K(alpha) - K(2 alpha) is the physical stiffness P and
    // K(2 alpha) - K(alpha) the drilling gauge alpha G, each to a few ulp of K.
    // CST/DKT adds its drilling penalty element by element, so its gauge sums
    // across a seam like any element term. MITC3+ fixes the gauge with one
    // relative-spin equation per node, area averaged over the node's fan
    // (assemble_mitc3p): a welded seam node averages over both sides, while
    // each tied component averages its own side and the tie sums the two
    // penalties. That gauge difference, alpha sum_i S_i (r_i - r)(r_i - r)^T,
    // lives on the seam fans only.
    let physical = |a: &ShellModel, a_2: &ShellModel| combine(&a.k, 2.0, &a_2.k, -1.0);
    let gauge = |a: &ShellModel, a_2: &ShellModel| combine(&a_2.k, 1.0, &a.k, -1.0);
    let stiffness_gaps = entry_gaps(&tied.model.k, &mono.k, &to_welded, &k_diagonal);
    let physical_gaps = entry_gaps(
        &physical(&tied.model, &tied_2.model),
        &physical(&mono, &mono_2),
        &to_welded,
        &k_diagonal,
    );
    let gauge_gaps = entry_gaps(
        &gauge(&tied.model, &tied_2.model),
        &gauge(&mono, &mono_2),
        &to_welded,
        &k_diagonal,
    );
    let mass_gaps = entry_gaps(&tied.model.m, &mono.m, &to_welded, &m_diagonal);
    // Welded DOFs of the seam fans: every node of a facet that touches a node
    // reached from two or more parts.
    let mut reached = vec![0usize; welded.nodes.len()];
    for map in &node_maps {
        for &w in map {
            reached[w] += 1;
        }
    }
    let dof_map = &mono.dof_map;
    let seam_fan_dofs: BTreeSet<usize> = welded
        .tris
        .iter()
        .filter(|tri| tri.iter().any(|&v| reached[v] > 1))
        .flat_map(|tri| tri.iter().copied())
        .flat_map(|w| (0..6).filter_map(move |d| dof_map[6 * w + d]))
        .collect();
    let gauge_gap_off_seam_fans = largest(
        gauge_gaps
            .iter()
            .filter(|((i, j), _)| !(seam_fan_dofs.contains(i) && seam_fan_dofs.contains(j)))
            .map(|(_, &g)| g),
    );

    let (lambda_max, lambda_max_gershgorin) = lambda_max_bracket(&mono);
    // Six rigid modes, then the elastic ones. The slice holds 32 so that none
    // of the compared modes sits at the window's top edge.
    let eigs_tied = lowest_eigenvalues_shifted(&tied.model, 32, 500.0);
    let eigs_mono = lowest_eigenvalues_shifted(&mono, 32, 500.0);
    for eigs in [&eigs_tied, &eigs_mono] {
        assert!(
            eigs[5].abs() < 1e-8 * eigs[6],
            "{}: six rigid modes, then elastic: {eigs:?}",
            kind.name()
        );
    }
    let lambda_1 = eigs_mono[6];
    PencilComparison {
        free: mono.free,
        seam_fan_dofs: seam_fan_dofs.len(),
        stiffness_entry_gap: largest(stiffness_gaps.into_values()),
        physical_entry_gap: largest(physical_gaps.into_values()),
        gauge_entry_gap: largest(gauge_gaps.into_values()),
        gauge_gap_off_seam_fans,
        mass_entry_gap: largest(mass_gaps.into_values()),
        lambda_1,
        lambda_max,
        lambda_max_gershgorin,
        eigenvalue_floor: G3_FLOOR_C * f64::EPSILON * lambda_max / lambda_1,
        elastic_gaps: (6..6 + G3_ELASTIC_MODES)
            .map(|m| (eigs_tied[m] - eigs_mono[m]).abs() / eigs_mono[m].abs())
            .collect(),
    }
}

#[test]
fn g3_matching_tie_reproduces_monolithic_pencil() {
    let started = std::time::Instant::now();
    let material = (11e9, 0.28, 0.0035, 480.0);
    let lx = 0.40;
    let ly = 0.30;
    let hz = 0.12;
    let (nx, ny, nz) = (6usize, 4usize, 3usize);

    // Flat split: two coplanar halves, seam x = lx/2.
    let left = rect_xy_mesh(0.0, 0.5 * lx, 0.0, ly, 0.0, nx / 2, ny);
    let right = rect_xy_mesh(0.5 * lx, lx, 0.0, ly, 0.0, nx / 2, ny);
    let flat_seams: [Vec<usize>; 2] = [
        (0..=ny).map(|j| j * (nx / 2 + 1) + nx / 2).collect(),
        (0..=ny).map(|j| j * (nx / 2 + 1)).collect(),
    ];
    // Folded L: flange y in [0, ly] and web z in [-hz, 0], seam y = z = 0.
    let flange = rect_xy_mesh(0.0, lx, 0.0, ly, 0.0, nx, ny);
    let web = rect_xz_mesh(0.0, lx, 0.0, -hz, 0.0, nx, nz);
    let flange_seam: Vec<usize> = (0..=nx).collect();
    let web_seam: Vec<usize> = (0..=nx).map(|i| nz * (nx + 1) + i).collect();
    // Folded T: a second flange half y in [-ly, 0] on the same seam, so each
    // seam edge is shared by three facets in the welded mesh.
    let flange_minus = rect_xy_mesh(0.0, lx, -ly, 0.0, 0.0, nx, ny);
    let minus_seam: Vec<usize> = (0..=nx).map(|i| ny * (nx + 1) + i).collect();

    // MITC3+ runs on the flat split only (owner ruling on the bead,
    // 2026-10-23). NodalDirectors requires every pair of facet normals at a
    // node to have a dot product of at least minimum_facet_dot > 0, and a 90
    // degree fold has exactly 0, so no conforming folded MITC3+ mesh exists to
    // compare with. MITC3+ folds are covered by the tied T-junction tests in
    // this file. CST/DKT carries global rotations on flat facets, so its
    // welded fold is an admissible monolithic mesh.
    // Averaged facet normals clear the director-to-facet test at a fold node
    // (dot 1/sqrt(2) or more), so the refusal at a seam node is the pairwise
    // facet-normal rule.
    for parts in [vec![&flange, &web], vec![&flange, &web, &flange_minus]] {
        let (welded, _) = weld(&parts);
        let mut normals = vec![[0.0; 3]; welded.nodes.len()];
        for (e, tri) in welded.tris.iter().enumerate() {
            let n = welded.facet(e).expect("valid facet").frame[2];
            for &v in tri {
                for c in 0..3 {
                    normals[v][c] += n[c];
                }
            }
        }
        let averaged: Vec<[f64; 3]> = normals
            .into_iter()
            .map(|n| n.map(|x| x / det::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2])))
            .collect();
        match NodalDirectors::new(&welded, averaged, f64::MIN_POSITIVE) {
            Err(PlateError::BadDirector { node }) => {
                let p = welded.nodes[node];
                assert!(
                    p[1] == 0.0 && p[2] == 0.0,
                    "the folded MITC3+ refusal is at a seam node, got node {node} at {p:?}"
                );
            }
            other => panic!("a folded MITC3+ mesh must refuse at the fold, got {other:?}"),
        }
    }

    let cases = [
        (
            "flat",
            ElementKind::CstDkt,
            tied_vs_welded(ElementKind::CstDkt, &[&left, &right], &flat_seams, material),
        ),
        (
            "flat",
            ElementKind::Mitc3p,
            tied_vs_welded(ElementKind::Mitc3p, &[&left, &right], &flat_seams, material),
        ),
        (
            "folded-L",
            ElementKind::CstDkt,
            tied_vs_welded(
                ElementKind::CstDkt,
                &[&flange, &web],
                &[flange_seam.clone(), web_seam.clone()],
                material,
            ),
        ),
        (
            "folded-T",
            ElementKind::CstDkt,
            tied_vs_welded(
                ElementKind::CstDkt,
                &[&flange, &web, &flange_minus],
                &[flange_seam, web_seam, minus_seam],
                material,
            ),
        ),
    ];

    let entry_gate = G3_ENTRY_ULPS * f64::EPSILON;
    for (geometry, kind, c) in &cases {
        let mut worst_mode = 0;
        let mut max_gap = 0.0_f64;
        for (m, &g) in c.elastic_gaps.iter().enumerate() {
            if g > max_gap {
                worst_mode = m + 1;
                max_gap = g;
            }
        }
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W2a\",\"test\":\"g3_matching_tie_reproduces_monolithic_pencil\",\"tier\":\"G3\",\"geometry\":\"{geometry}\",\"element\":\"{}\",\"free_dofs\":{},\"seam_fan_dofs\":{},\"entry_gate\":{entry_gate:e},\"stiffness_entry_gap\":{:e},\"physical_stiffness_entry_gap\":{:e},\"drilling_gauge_entry_gap\":{:e},\"drilling_gauge_gap_off_seam_fans\":{:e},\"mass_entry_gap\":{:e},\"lambda_1\":{:e},\"lambda_max_power\":{:e},\"lambda_max_gershgorin\":{:e},\"floor_c\":{G3_FLOOR_C},\"eigenvalue_floor\":{:e},\"elastic_modes_checked\":{G3_ELASTIC_MODES},\"max_rel_eig_gap\":{max_gap:e},\"worst_elastic_mode\":{worst_mode},\"elastic_gaps\":{:?},\"elapsed_s\":{}}}",
            kind.name(),
            c.free,
            c.seam_fan_dofs,
            c.stiffness_entry_gap,
            c.physical_entry_gap,
            c.gauge_entry_gap,
            c.gauge_gap_off_seam_fans,
            c.mass_entry_gap,
            c.lambda_1,
            c.lambda_max,
            c.lambda_max_gershgorin,
            c.eigenvalue_floor,
            c.elastic_gaps,
            started.elapsed().as_secs_f64()
        );
    }
    for (geometry, kind, c) in &cases {
        // Mass and physical stiffness equal the welded pencil entry by entry;
        // a drilling-gauge difference may sit on the seam fans only, and only
        // for MITC3+'s per-node gauge.
        assert!(
            c.mass_entry_gap <= entry_gate
                && c.physical_entry_gap <= entry_gate
                && c.gauge_gap_off_seam_fans <= entry_gate,
            "{} {geometry}: tied pencil entries differ from the welded ones by M {:e}, physical K {:e}, drilling gauge off the seam fans {:e} > {entry_gate:e}",
            kind.name(),
            c.mass_entry_gap,
            c.physical_entry_gap,
            c.gauge_gap_off_seam_fans
        );
        if *kind == ElementKind::CstDkt {
            assert!(
                c.stiffness_entry_gap <= entry_gate,
                "{geometry}: CST/DKT adds its drilling penalty per element, so its whole K must match: gap {:e} > {entry_gate:e}",
                c.stiffness_entry_gap
            );
        }
        for (m, &g) in c.elastic_gaps.iter().enumerate() {
            assert!(
                g <= c.eigenvalue_floor,
                "{} {geometry}: elastic mode {} differs by {g:e}, above the two-solve floor {:e}",
                kind.name(),
                m + 1,
                c.eigenvalue_floor
            );
        }
    }
}

/// CST/DKT (ratio, 1-based mode) pairs whose mortar ladder converges FASTER
/// than the monolithic one on these rungs, so the two-sided band cannot hold
/// for them. They are gated one-sided (owner ruling on the bead, 2026-10-23).
///
/// Measured without a reference, from signed successive eigenvalue changes at
/// k = 8, 12, 16, 24, 32 (ts1 log
/// 20261022T195820Z-d6dd965dc0f4-MagentaPanther-20261022T195820Z-b.log), as
/// local orders over consecutive rung triples:
/// - mode 3: monolithic 2.00, 1.99, 1.99; 2:3 mortar 2.28, 2.13, 2.09; 1:2
///   mortar 2.49, 2.27.
/// - mode 8: monolithic 2.67, 2.23, 2.08; 2:3 mortar 2.97, 2.66, 2.53; 1:2
///   mortar 3.25, 3.14.
///
/// Every change has the monolithic sign (DKT converges from below), and the
/// mortar orders fall toward 2 from above: a smaller leading error with
/// larger higher-order terms, not a loss of consistency. One finer rung does
/// not bring mode 8 into the band (2.53 against 2.08 on the finest 2:3
/// triple), so the ruling's one-sided fallback applies to these four pairs.
/// MITC3+ holds the two-sided band on all twenty.
const G1_SUPERCONVERGENT_CST_DKT: [(&str, usize); 4] =
    [("2:3", 3), ("2:3", 8), ("1:2", 3), ("1:2", 8)];

#[test]
fn g1_nonmatching_mortar_converges_to_monolithic() {
    let started = std::time::Instant::now();
    let e_pa = 10e9;
    let nu = 0.30;
    let h_m = 0.004;
    let rho = 500.0;
    let lx = 0.36;
    let ly = 0.24;

    // Fixture: a plate clamped on all four edges, with a coplanar seam at
    // x = lx / 2. Free edges carry a Reissner-Mindlin boundary layer about
    // h_m wide. With free y-edges (clamped at x = 0 and lx only), the MITC3+
    // monolithic ladder did not converge at 8..128 cells per half: mode 2
    // ran 1.6392e6, 1.6434e6, 1.6430e6, 1.6410e6, 1.6392e6 as the cells shrank
    // through h_m. All-clamped, the successive changes shrink 4.0-4.4x from
    // 16 cells per half.
    //
    // On rung k the coarse (slave) half has 2k x 2k cells and the fine
    // (master) half 3k (2:3) or 4k (1:2); the monolithic family has 2k per
    // half on the same rungs. The ladder starts at 16 cells per half. At 8,
    // the near pair at 38.29e6 and 38.74e6 (modes 3 and 4) crosses: its modes
    // mix (shape MAC 0.79). Eigenvalues are compared in sorted order, as
    // min-max defines "the first ten". The reference has 128 cells per half,
    // four times the finest coarse side, so its own error is about 1/16 of
    // the finest rung's.
    let rungs = [8usize, 12, 16];
    let n_ref = 128usize;
    let hs: Vec<f64> = rungs.iter().map(|&k| 1.0 / k as f64).collect();
    let edge = |p: [f64; 3], x_edge: f64| p[0] == x_edge || p[1] == 0.0 || p[1] == ly;
    let local_orders = |errors: &[f64]| -> Vec<f64> {
        errors
            .windows(2)
            .zip(hs.windows(2))
            .map(|(e, h)| det::ln(e[0] / e[1]) / det::ln(h[0] / h[1]))
            .collect()
    };
    let mut failures = Vec::new();

    for kind in [ElementKind::CstDkt, ElementKind::Mitc3p] {
        // Monolithic plate with n x n cells per half.
        let plate = |n: usize| -> ShellModel {
            let mesh = rect_xy_mesh(0.0, lx, 0.0, ly, 0.0, 2 * n, n);
            let clamped: Vec<usize> = (0..mesh.nodes.len())
                .filter(|&i| edge(mesh.nodes[i], 0.0) || mesh.nodes[i][0] == lx)
                .collect();
            assemble_component_model(kind, &mesh, e_pa, nu, h_m, rho, &clamped, DRILLING_ALPHA)
        };
        // Each window search starts at 1.1x the nearest coarser solve's tenth
        // eigenvalue and grows by 10 % until inertia certifies ten modes.
        // Doubling from scratch took about 20 factorizations of the reference.
        let mut mono_eigs: Vec<Vec<f64>> = Vec::new();
        for &k in &rungs {
            let model = plate(2 * k);
            let eigs = match mono_eigs.last() {
                None => lowest_positive_eigenvalues(&model, 10),
                Some(previous) => positive_eigenvalues_from(&model, 10, 1.1 * previous[9], 1.1),
            };
            mono_eigs.push(eigs);
        }
        let reference =
            positive_eigenvalues_from(&plate(n_ref), 10, 1.1 * mono_eigs.last().unwrap()[9], 1.1);
        let errors = |eigs: &Vec<f64>| -> Vec<f64> {
            eigs.iter()
                .zip(&reference)
                .map(|(lam, r)| (lam - r).abs() / r)
                .collect()
        };
        let mono: Vec<Vec<f64>> = mono_eigs.iter().map(errors).collect();

        for (ratio, fine) in [("2:3", 3usize), ("1:2", 4)] {
            let mortar: Vec<Vec<f64>> = rungs
                .iter()
                .zip(&mono_eigs)
                .map(|(&k, same_rung)| {
                    let n_a = fine * k;
                    let n_b = 2 * k;
                    let left_mesh = rect_xy_mesh(0.0, 0.5 * lx, 0.0, ly, 0.0, n_a, n_a);
                    let right_mesh = rect_xy_mesh(0.5 * lx, lx, 0.0, ly, 0.0, n_b, n_b);
                    let left_clamped: Vec<usize> = (0..left_mesh.nodes.len())
                        .filter(|&i| edge(left_mesh.nodes[i], 0.0))
                        .collect();
                    let right_clamped: Vec<usize> = (0..right_mesh.nodes.len())
                        .filter(|&i| edge(right_mesh.nodes[i], lx))
                        .collect();
                    let mod_left = assemble_component_model(
                        kind,
                        &left_mesh,
                        e_pa,
                        nu,
                        h_m,
                        rho,
                        &left_clamped,
                        DRILLING_ALPHA,
                    );
                    let mod_right = assemble_component_model(
                        kind,
                        &right_mesh,
                        e_pa,
                        nu,
                        h_m,
                        rho,
                        &right_clamped,
                        DRILLING_ALPHA,
                    );
                    let h_left = vec![h_m; left_mesh.nodes.len()];
                    let h_right = vec![h_m; right_mesh.nodes.len()];
                    let comps = [
                        ShellComponent {
                            id: 0,
                            mesh: &left_mesh,
                            model: &mod_left,
                            nodal_thickness_m: &h_left,
                        },
                        ShellComponent {
                            id: 1,
                            mesh: &right_mesh,
                            model: &mod_right,
                            nodal_thickness_m: &h_right,
                        },
                    ];
                    // Finer side as master, coarser side as slave. A coplanar seam
                    // ties all six DOFs, as the monolithic mesh shares them.
                    let ties = [EdgeTie {
                        a: TieSide {
                            component: 0,
                            nodes: (0..=n_a).map(|j| j * (n_a + 1) + n_a).collect(),
                        },
                        b: TieSide {
                            component: 1,
                            nodes: (0..=n_b).map(|j| j * (n_b + 1)).collect(),
                        },
                        offset: TieOffset::Midsurface,
                        rotations: RotationCoupling::Full,
                        kind: TieKind::Rigid,
                    }];
                    let tied = tie_shells(&comps, &ties, default_tie_budget(), || false)
                        .expect("nonmatching tie succeeds");
                    errors(&positive_eigenvalues_from(
                        &tied.model,
                        10,
                        1.1 * same_rung[9],
                        1.1,
                    ))
                })
                .collect();

            for m in 0..10 {
                let column = |rows: &[Vec<f64>]| rows.iter().map(|r| r[m]).collect::<Vec<_>>();
                let (mono_errors, mortar_errors) = (column(&mono), column(&mortar));
                let own = fit_order(
                    &RefinementLadder::new(hs.clone(), mono_errors.clone())
                        .expect("valid monolithic ladder"),
                );
                let ladder = RefinementLadder::new(hs.clone(), mortar_errors.clone())
                    .expect("valid mortar ladder");
                let observed = fit_order(&ladder).observed;
                let superconvergent = kind == ElementKind::CstDkt
                    && G1_SUPERCONVERGENT_CST_DKT.contains(&(ratio, m + 1));
                let case = format!("{}-{ratio}-mode-{}", kind.name(), m + 1);
                println!(
                    "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W2a\",\"test\":\"g1_nonmatching_mortar_converges_to_monolithic\",\"tier\":\"G1\",\"case\":\"{case}\",\"fixture\":\"clamped-all-edges\",\"rungs_k\":{rungs:?},\"reference_cells_per_half\":{n_ref},\"monolithic_order\":{},\"monolithic_rms_residual\":{:e},\"mortar_order\":{observed},\"gate\":\"{}\",\"monolithic_local_orders\":{:?},\"mortar_local_orders\":{:?},\"monolithic_errors\":{mono_errors:?},\"mortar_errors\":{mortar_errors:?}}}",
                    own.observed,
                    own.rms_residual,
                    if superconvergent {
                        "one-sided"
                    } else {
                        "two-sided"
                    },
                    local_orders(&mono_errors),
                    local_orders(&mortar_errors)
                );
                // The mortar error falls on every rung.
                if !mortar_errors.windows(2).all(|w| w[1] < w[0]) {
                    println!("{case}: mortar error does not fall on every rung");
                    failures.push(case);
                    continue;
                }
                if superconvergent {
                    // No order loss, and the declaration is still live: a pair
                    // that returns to the band must leave the list.
                    if observed < own.observed - ORDER_GATE_TOLERANCE {
                        println!(
                            "{case}: mortar order {observed} below monolithic {} - {ORDER_GATE_TOLERANCE}",
                            own.observed
                        );
                        failures.push(case);
                    } else if observed <= own.observed + ORDER_GATE_TOLERANCE {
                        println!(
                            "{case}: declared superconvergent but inside the two-sided band; drop it from G1_SUPERCONVERGENT_CST_DKT"
                        );
                        failures.push(case);
                    }
                    continue;
                }
                // The mortar converges at the monolithic mesh's own order.
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
        }
    }
    println!(
        "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W2a\",\"test\":\"g1_nonmatching_mortar_converges_to_monolithic\",\"failures\":{failures:?},\"elapsed_s\":{}}}",
        started.elapsed().as_secs_f64()
    );
    assert!(
        failures.is_empty(),
        "mortar order gate failed for {failures:?}"
    );
}

#[test]
fn g1_skin_offset_t_junction() {
    let started = std::time::Instant::now();
    let e_pa = 12e9;
    let nu = 0.25;
    let h_plate = 0.008;
    let h_rib = 0.006;
    let rho = 450.0;
    let lx = 0.36;
    let ly = 0.04;
    let h_web = 0.012;
    let nx = 8usize;
    let ny = 2usize;
    let nz = 2usize;

    // Horizontal plate in z = 0 plane, clamped at x = 0 and x = lx.
    let plate_mesh = rect_xy_mesh(0.0, lx, -0.5 * ly, 0.5 * ly, 0.0, nx, ny);
    let plate_clamped: Vec<usize> = (0..=ny)
        .flat_map(|j| [j * (nx + 1), j * (nx + 1) + nx])
        .collect();
    let master_seam: Vec<usize> = (0..=nx).map(|i| 1 * (nx + 1) + i).collect();

    // Downward rib in y = 0 plane, z in [-h_web, 0] (seam at top edge z = 0 -> row k = nz).
    let rib_below_mesh = rect_xz_mesh(0.0, lx, 0.0, -h_web, 0.0, nx, nz);
    let rib_below_clamped: Vec<usize> = (0..=nz)
        .flat_map(|k| [k * (nx + 1), k * (nx + 1) + nx])
        .collect();
    let slave_seam_below: Vec<usize> = (0..=nx).map(|i| nz * (nx + 1) + i).collect();

    // Upward rib in y = 0 plane, z in [0, +h_web] (seam at bottom edge z = 0 -> row k = 0).
    let rib_above_mesh = rect_xz_mesh(0.0, lx, 0.0, 0.0, h_web, nx, nz);
    let rib_above_clamped: Vec<usize> = (0..=nz)
        .flat_map(|k| [k * (nx + 1), k * (nx + 1) + nx])
        .collect();
    let slave_seam_above: Vec<usize> = (0..=nx).collect();

    for kind in [ElementKind::CstDkt, ElementKind::Mitc3p] {
        let mod_plate = assemble_component_model(
            kind,
            &plate_mesh,
            e_pa,
            nu,
            h_plate,
            rho,
            &plate_clamped,
            DRILLING_ALPHA,
        );
        let mod_rib_below = assemble_component_model(
            kind,
            &rib_below_mesh,
            e_pa,
            nu,
            h_rib,
            rho,
            &rib_below_clamped,
            DRILLING_ALPHA,
        );
        let mod_rib_above = assemble_component_model(
            kind,
            &rib_above_mesh,
            e_pa,
            nu,
            h_rib,
            rho,
            &rib_above_clamped,
            DRILLING_ALPHA,
        );

        let h0 = vec![h_plate; plate_mesh.nodes.len()];
        let h1 = vec![h_rib; rib_below_mesh.nodes.len()];

        let eval_tied = |rib_mesh: &ShellMesh,
                         mod_rib: &ShellModel,
                         slave_seam: &[usize],
                         offset: TieOffset|
         -> TiedShellModel {
            let comps = [
                ShellComponent {
                    id: 0,
                    mesh: &plate_mesh,
                    model: &mod_plate,
                    nodal_thickness_m: &h0,
                },
                ShellComponent {
                    id: 1,
                    mesh: rib_mesh,
                    model: mod_rib,
                    nodal_thickness_m: &h1,
                },
            ];
            let ties = [EdgeTie {
                a: TieSide {
                    component: 0,
                    nodes: master_seam.clone(),
                },
                b: TieSide {
                    component: 1,
                    nodes: slave_seam.to_vec(),
                },
                offset,
                rotations: RotationCoupling::Full,
                kind: TieKind::Rigid,
            }];
            tie_shells(&comps, &ties, default_tie_budget(), || false)
                .expect("skin-offset T-junction tie succeeds")
        };

        // 1. Downward rib: Skin { side: -1 } increases effective lever arm by h_plate / 2
        // (parallel-axis theorem: I_eff increases by A_rib * (h_plate/2)^2 + ...),
        // increasing the fundamental longitudinal bending frequency relative to Midsurface.
        let tied_mid_below = eval_tied(
            &rib_below_mesh,
            &mod_rib_below,
            &slave_seam_below,
            TieOffset::Midsurface,
        );
        let tied_skin_minus = eval_tied(
            &rib_below_mesh,
            &mod_rib_below,
            &slave_seam_below,
            TieOffset::Skin { side: -1 },
        );
        let tied_skin_plus_below = eval_tied(
            &rib_below_mesh,
            &mod_rib_below,
            &slave_seam_below,
            TieOffset::Skin { side: 1 },
        );

        let f_mid = det::sqrt(lowest_positive_eigenvalues(&tied_mid_below.model, 1)[0])
            / (2.0 * core::f64::consts::PI);
        let f_minus = det::sqrt(lowest_positive_eigenvalues(&tied_skin_minus.model, 1)[0])
            / (2.0 * core::f64::consts::PI);
        let f_plus_below =
            det::sqrt(lowest_positive_eigenvalues(&tied_skin_plus_below.model, 1)[0])
                / (2.0 * core::f64::consts::PI);

        assert!(
            f_minus > f_mid * 1.01,
            "{}: Skin {{ side: -1 }} on downward rib must increase f1 over Midsurface: f_minus={f_minus}, f_mid={f_mid}",
            kind.name()
        );
        assert!(
            f_plus_below < f_mid * 0.99,
            "{}: flipping side to +1 on downward rib must decrease f1 below Midsurface: f_plus={f_plus_below}, f_mid={f_mid}",
            kind.name()
        );

        // 2. Also verify asymmetric bending-membrane coupling sign flip when rib sits above (side = +1)
        // vs below (side = -1) under a uniform downward transverse load F_z on the plate centerline:
        // Neutral axis shifts toward the rib, so plate midsurface axial displacement u_x at x = lx/4
        // has opposite sign for rib-below (side = -1) vs rib-above (side = +1)!
        let tied_skin_plus_above = eval_tied(
            &rib_above_mesh,
            &mod_rib_above,
            &slave_seam_above,
            TieOffset::Skin { side: 1 },
        );

        let ux_below = midspan_bending_membrane_ux(
            &plate_mesh,
            &rib_below_mesh,
            &tied_skin_minus,
            &master_seam,
        );
        let ux_above = midspan_bending_membrane_ux(
            &plate_mesh,
            &rib_above_mesh,
            &tied_skin_plus_above,
            &master_seam,
        );
        assert!(
            ux_below < 0.0 && ux_above > 0.0 && (ux_below + ux_above).abs() < 0.05 * ux_below.abs(),
            "{}: bending-membrane coupling u_x must flip sign symmetrically between rib-below (side=-1) and rib-above (side=+1): ux_below={ux_below:e}, ux_above={ux_above:e}",
            kind.name()
        );

        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W2a\",\"test\":\"g1_skin_offset_t_junction\",\"tier\":\"G1\",\"element\":\"{}\",\"f1_midsurface_hz\":{},\"f1_skin_minus_hz\":{},\"f1_skin_plus_below_hz\":{},\"ux_below\":{:e},\"ux_above\":{:e},\"elapsed_s\":{}}}",
            kind.name(),
            f_mid,
            f_minus,
            f_plus_below,
            ux_below,
            ux_above,
            started.elapsed().as_secs_f64()
        );
    }
}

fn midspan_bending_membrane_ux(
    plate_mesh: &ShellMesh,
    rib_mesh: &ShellMesh,
    tied: &TiedShellModel,
    master_seam: &[usize],
) -> f64 {
    let mut f0 = vec![0.0; 6 * plate_mesh.nodes.len()];
    let f1 = vec![0.0; 6 * rib_mesh.nodes.len()];
    // Apply +z load at midspan seam node.
    let mid_node = master_seam[master_seam.len() / 2];
    f0[6 * mid_node + 2] = 100.0;
    let f_r = tied.project_loads(&[&f0, &f1]).unwrap();
    let sym = SymbolicLdlt::analyze(&tied.model.k, DirectOrdering::Amd).unwrap();
    let num = sym.factor(&tied.model.k, &LdltOptions::default()).unwrap();
    let u_r = num.solve(&f_r);
    let full = tied.expand_components(&u_r).unwrap();
    let quarter_node = master_seam[1];
    full[0][6 * quarter_node]
}

#[test]
fn g0_drilling_alpha_insensitivity_at_junction() {
    let started = std::time::Instant::now();
    let e_pa = 10e9;
    let nu = 0.28;
    let h_m = 0.006;
    let rho = 460.0;
    let lx = 0.30;
    let ly = 0.20;
    let hz = 0.12;
    let nx = 8usize;
    let ny = 4usize;
    let nz = 4usize;

    let flange_mesh = rect_xy_mesh(0.0, lx, -0.5 * ly, 0.5 * ly, 0.0, nx, ny);
    let web_mesh = rect_xz_mesh(0.0, lx, 0.0, -hz, 0.0, nx, nz);
    let master_seam: Vec<usize> = (0..=nx).map(|i| (ny / 2) * (nx + 1) + i).collect();
    let slave_seam: Vec<usize> = (0..=nx).map(|i| nz * (nx + 1) + i).collect();

    // A fold ties only the rotation about the seam. `Full` would also tie the
    // flange's bending slope to the web's drilling spin (and the reverse), so
    // the drilling penalty would enter the junction. Measured at origin
    // 8378586d9, across alpha 1e-4..1e-2: Full moves the first ten
    // frequencies by 6.6e-4 (CST/DKT) and 2.4e-3 (MITC3+), SeamTangent by
    // 3.2e-6 and 1.4e-6, and the flange alone by 1.6e-13 and 1.2e-13.
    for kind in [ElementKind::CstDkt, ElementKind::Mitc3p] {
        let solve_for_alpha = |alpha: f64| -> Vec<f64> {
            let mod0 = assemble_component_model(kind, &flange_mesh, e_pa, nu, h_m, rho, &[], alpha);
            let mod1 = assemble_component_model(kind, &web_mesh, e_pa, nu, h_m, rho, &[], alpha);
            let h0 = vec![h_m; flange_mesh.nodes.len()];
            let h1 = vec![h_m; web_mesh.nodes.len()];
            let comps = [
                ShellComponent {
                    id: 0,
                    mesh: &flange_mesh,
                    model: &mod0,
                    nodal_thickness_m: &h0,
                },
                ShellComponent {
                    id: 1,
                    mesh: &web_mesh,
                    model: &mod1,
                    nodal_thickness_m: &h1,
                },
            ];
            let ties = [EdgeTie {
                a: TieSide {
                    component: 0,
                    nodes: master_seam.clone(),
                },
                b: TieSide {
                    component: 1,
                    nodes: slave_seam.clone(),
                },
                offset: TieOffset::Midsurface,
                rotations: RotationCoupling::SeamTangent,
                kind: TieKind::Rigid,
            }];
            let tied = tie_shells(&comps, &ties, default_tie_budget(), || false).unwrap();
            let eigs = lowest_eigenvalues_shifted(&tied.model, 16, 500.0);
            eigs.into_iter()
                .skip(6)
                .take(10)
                .map(|lam| det::sqrt(lam) / (2.0 * core::f64::consts::PI))
                .collect()
        };

        let freqs_lo = solve_for_alpha(1e-4);
        let freqs_mid = solve_for_alpha(1e-3);
        let freqs_hi = solve_for_alpha(1e-2);

        let mut max_rel_change = 0.0_f64;
        for m in 0..10 {
            let rel_lo_hi = (freqs_hi[m] - freqs_lo[m]).abs() / freqs_mid[m];
            max_rel_change = max_rel_change.max(rel_lo_hi);
            assert!(
                rel_lo_hi < 1e-4,
                "{}: mode {} frequency changed by {rel_lo_hi:e} >= 1e-4 across alpha in [1e-4, 1e-2] (f_lo={}, f_mid={}, f_hi={})",
                kind.name(),
                m + 1,
                freqs_lo[m],
                freqs_mid[m],
                freqs_hi[m]
            );
        }

        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W2a\",\"test\":\"g0_drilling_alpha_insensitivity_at_junction\",\"tier\":\"G0\",\"element\":\"{}\",\"max_rel_freq_change\":{:e},\"elapsed_s\":{}}}",
            kind.name(),
            max_rel_change,
            started.elapsed().as_secs_f64()
        );
    }
}
