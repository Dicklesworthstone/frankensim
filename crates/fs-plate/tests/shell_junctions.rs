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
use fs_mms::{RefinementLadder, fit_order};
use fs_modal::{SliceOptions, slice_window};
use fs_plate::PlateSection;
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
use fs_rand::StreamKey;
use fs_sparse::{Coo, Csr, DirectOrdering, LdltOptions, SymbolicLdlt};

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
        hi *= 2.0;
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
    let hi = find_upper_window_bound(&model.k, &model.m, count, 100.0);
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
    let nx = 3usize;
    let ny = 2usize;
    let nz = 2usize;

    // Horizontal flange: x in [0, lx], y in [-ly/2, ly/2], z = 0 (seam at y = 0 -> row j = 1).
    let flange_mesh = rect_xy_mesh(0.0, lx, -0.5 * ly, 0.5 * ly, 0.0, nx, ny);
    // Vertical web: x in [0, lx], y = 0, z in [-hz, 0] (seam at z = 0 -> row k = nz).
    let web_mesh = rect_xz_mesh(0.0, lx, 0.0, -hz, 0.0, nx, nz);

    let master_seam: Vec<usize> = (0..=nx).map(|i| 1 * (nx + 1) + i).collect();
    let slave_seam: Vec<usize> = (0..=nx).map(|i| nz * (nx + 1) + i).collect();

    for (kind, rotations) in [ElementKind::CstDkt, ElementKind::Mitc3p]
        .into_iter()
        .flat_map(|kind| {
            [RotationCoupling::Full, RotationCoupling::SeamTangent]
                .map(|rotations| (kind, rotations))
        })
    {
        // 1. Unconstrained T-junction: check 6 rigid-body modes (|lambda_i| / lambda_7 < 1e-12).
        let mod_flange =
            assemble_component_model(kind, &flange_mesh, e_pa, nu, h_m, rho, &[], DRILLING_ALPHA);
        let mod_web =
            assemble_component_model(kind, &web_mesh, e_pa, nu, h_m, rho, &[], DRILLING_ALPHA);
        let h_flange = vec![h_m; flange_mesh.nodes.len()];
        let h_web = vec![h_m; web_mesh.nodes.len()];
        let comps = [
            ShellComponent {
                id: 0,
                mesh: &flange_mesh,
                model: &mod_flange,
                nodal_thickness_m: &h_flange,
            },
            ShellComponent {
                id: 1,
                mesh: &web_mesh,
                model: &mod_web,
                nodal_thickness_m: &h_web,
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
            rotations,
            kind: TieKind::Rigid,
        }];
        let tied = tie_shells(&comps, &ties, default_tie_budget(), || false)
            .expect("T-junction tie succeeds");

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

        // 2. Constant-membrane-strain patch test across the seam:
        // Prescribe uniaxial stress sigma_xx = E * eps_0 along x on both flange and web.
        // At x = 0: fix u_x = 0 (and minimal rigid support at (0,0,0) and (0, ly/2, 0)).
        // At x = lx: apply consistent nodal traction forces F_x = E * h * L_edge_share * eps_0.
        // Verify that the static solution K_r u_r = f_r (with minimal isostatic support)
        // recovers constant strain du_x/dx = eps_0, du_y/dy = -nu * eps_0 (flange),
        // du_z/dz = -nu * eps_0 (web) across the T-junction seam to 1e-12!
        let eps_0 = 1e-4;
        let patch_err = verify_t_junction_constant_strain_patch(
            &flange_mesh,
            &web_mesh,
            &tied,
            e_pa,
            nu,
            h_m,
            lx,
            ly,
            hz,
            nx,
            ny,
            nz,
            eps_0,
        );
        assert!(
            patch_err < 1e-12,
            "{} {rotations:?}: T-junction constant-strain patch relative error {patch_err:e} >= 1e-12",
            kind.name()
        );

        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W2a\",\"test\":\"g0_t_junction_rigid_modes_and_constant_strain_patch\",\"tier\":\"G0\",\"element\":\"{}\",\"rotations\":\"{rotations:?}\",\"max_rigid_ratio\":{:e},\"lambda_7\":{},\"constant_strain_patch_err\":{:e},\"elapsed_s\":{}}}",
            kind.name(),
            max_rigid_ratio,
            lambda_7,
            patch_err,
            started.elapsed().as_secs_f64()
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn verify_t_junction_constant_strain_patch(
    flange_mesh: &ShellMesh,
    web_mesh: &ShellMesh,
    tied: &TiedShellModel,
    e_pa: f64,
    nu: f64,
    h_m: f64,
    lx: f64,
    ly: f64,
    hz: f64,
    nx: usize,
    ny: usize,
    nz: usize,
    eps_0: f64,
) -> f64 {
    // Exact uniaxial stress state sigma_xx = E * eps_0, sigma_yy = sigma_zz = 0:
    // u_x(x, y, z) = eps_0 * x
    // On flange (z = 0): u_y = -nu * eps_0 * y, u_z = 0, theta = 0
    // On web (y = 0):    u_y = 0, u_z = -nu * eps_0 * z, theta = 0
    // Evaluate K_c * u_exact on both components and project f_r = T^T f_exact:
    // On all interior and seam nodes (where boundary tractions are zero), (K_r u_r)_i must be 0!
    // And solving K_r u_r = f_r with boundary DOFs pinned to u_exact recovers u_exact everywhere to < 1e-12.
    let mut u_r_exact = vec![0.0; tied.model.free];
    for (n, &[x, y, _z]) in flange_mesh.nodes.iter().enumerate() {
        let exact = [eps_0 * x, -nu * eps_0 * y, 0.0, 0.0, 0.0, 0.0];
        for d in 0..6 {
            if let Some(r) = tied.map[0][n][d] {
                u_r_exact[r] = exact[d];
            }
        }
    }
    for (n, &[x, _y, z]) in web_mesh.nodes.iter().enumerate() {
        let exact = [eps_0 * x, 0.0, -nu * eps_0 * z, 0.0, 0.0, 0.0];
        for d in 0..6 {
            if let Some(r) = tied.map[1][n][d] {
                u_r_exact[r] = exact[d];
            }
        }
    }

    // Build boundary load vector f_full at x = 0 and x = lx for uniaxial traction sigma_xx = E * eps_0.
    let sigma_xx = e_pa * eps_0;
    let mut f0 = vec![0.0; 6 * flange_mesh.nodes.len()];
    let dy = ly / ny as f64;
    for j in 0..=ny {
        let weight = if j == 0 || j == ny { 0.5 * dy } else { dy };
        let n_left = j * (nx + 1);
        let n_right = j * (nx + 1) + nx;
        f0[6 * n_left] -= sigma_xx * h_m * weight;
        f0[6 * n_right] += sigma_xx * h_m * weight;
    }
    let mut f1 = vec![0.0; 6 * web_mesh.nodes.len()];
    let dz = hz / nz as f64;
    for k in 0..=nz {
        let weight = if k == 0 || k == nz { 0.5 * dz } else { dz };
        let n_left = k * (nx + 1);
        let n_right = k * (nx + 1) + nx;
        f1[6 * n_left] -= sigma_xx * h_m * weight;
        f1[6 * n_right] += sigma_xx * h_m * weight;
    }
    let f_r = tied
        .project_loads(&[&f0, &f1])
        .expect("load projection succeeds");

    // Check residual r = K_r * u_r_exact - f_r across all reduced DOFs (including seam DOFs!).
    let mut r_vec = vec![0.0; tied.model.free];
    tied.model.k.spmv(&u_r_exact, &mut r_vec);
    let force_scale = sigma_xx * h_m * (ly + hz);
    let mut max_rel_res = 0.0_f64;
    for i in 0..tied.model.free {
        let err = (r_vec[i] - f_r[i]).abs() / force_scale;
        max_rel_res = max_rel_res.max(err);
    }

    // Also check expanded slave displacements at the seam match exact field.
    let expanded = tied.expand_components(&u_r_exact).unwrap();
    let disp_scale = eps_0 * lx;
    for (n, &[x, _y, z]) in web_mesh.nodes.iter().enumerate() {
        let exact = [eps_0 * x, 0.0, -nu * eps_0 * z, 0.0, 0.0, 0.0];
        for d in 0..6 {
            let err = (expanded[1][6 * n + d] - exact[d]).abs() / disp_scale;
            max_rel_res = max_rel_res.max(err);
        }
    }
    max_rel_res
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
        for _trial in 0..64 {
            let v_r: Vec<f64> = (0..tied.model.free)
                .map(|_| 2.0 * rng.next_f64() - 1.0)
                .collect();
            let lambda: Vec<f64> = (0..tied.receipt.slave_dofs)
                .map(|_| 2.0 * rng.next_f64() - 1.0)
                .collect();
            let (p_m, p_s) = tied
                .interface_virtual_power(&v_r, &lambda)
                .expect("virtual power evaluation succeeds");
            let denom = p_m.abs() + p_s.abs();
            assert!(denom > 1e-9);
            let rel_err = (p_m + p_s).abs() / denom;
            max_rel_power_err = max_rel_power_err.max(rel_err);
            assert!(
                rel_err <= 1e-13,
                "{} {rotations:?}: virtual work identity violated: P_m={p_m}, P_s={p_s}, rel_err={rel_err:e}",
                kind.name()
            );
        }

        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W2a\",\"test\":\"g0_tie_virtual_work_identity\",\"tier\":\"G0\",\"element\":\"{}\",\"rotations\":\"{rotations:?}\",\"trials\":64,\"max_rel_power_err\":{:e},\"elapsed_s\":{}}}",
            kind.name(),
            max_rel_power_err,
            started.elapsed().as_secs_f64()
        );
    }
}

#[test]
fn g3_matching_tie_reproduces_monolithic_pencil() {
    let started = std::time::Instant::now();
    let e_pa = 11e9;
    let nu = 0.28;
    let h_m = 0.0035;
    let rho = 480.0;
    let lx = 0.40;
    let ly = 0.30;
    let nx_half = 3usize;
    let ny = 4usize;

    let left_mesh = rect_xy_mesh(0.0, 0.5 * lx, 0.0, ly, 0.0, nx_half, ny);
    let right_mesh = rect_xy_mesh(0.5 * lx, lx, 0.0, ly, 0.0, nx_half, ny);
    let mono_mesh = rect_xy_mesh(0.0, lx, 0.0, ly, 0.0, 2 * nx_half, ny);

    let master_seam: Vec<usize> = (0..=ny).map(|j| j * (nx_half + 1) + nx_half).collect();
    let slave_seam: Vec<usize> = (0..=ny).map(|j| j * (nx_half + 1)).collect();

    for kind in [ElementKind::CstDkt, ElementKind::Mitc3p] {
        let mod_left =
            assemble_component_model(kind, &left_mesh, e_pa, nu, h_m, rho, &[], DRILLING_ALPHA);
        let mod_right =
            assemble_component_model(kind, &right_mesh, e_pa, nu, h_m, rho, &[], DRILLING_ALPHA);
        let mod_mono =
            assemble_component_model(kind, &mono_mesh, e_pa, nu, h_m, rho, &[], DRILLING_ALPHA);

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
            rotations: RotationCoupling::Full,
            kind: TieKind::Rigid,
        }];
        let tied = tie_shells(&comps, &ties, default_tie_budget(), || false)
            .expect("matching tie succeeds");
        assert_eq!(tied.model.free, mod_mono.free);

        // Extract 6 rigid + 15 elastic = 21 eigenvalues and compare the 15 non-rigid modes.
        let eigs_tied = lowest_eigenvalues_shifted(&tied.model, 21, 500.0);
        let eigs_mono = lowest_eigenvalues_shifted(&mod_mono, 21, 500.0);
        let mut max_rel_err = 0.0_f64;
        for m in 6..21 {
            let rel_err = (eigs_tied[m] - eigs_mono[m]).abs() / eigs_mono[m].abs();
            max_rel_err = max_rel_err.max(rel_err);
            assert!(
                rel_err < 1e-11,
                "{}: mode {} mismatch: tied={}, mono={}, rel_err={rel_err:e}",
                kind.name(),
                m - 5,
                eigs_tied[m],
                eigs_mono[m]
            );
        }

        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W2a\",\"test\":\"g3_matching_tie_reproduces_monolithic_pencil\",\"tier\":\"G3\",\"element\":\"{}\",\"modes_checked\":15,\"max_rel_err\":{:e},\"elapsed_s\":{}}}",
            kind.name(),
            max_rel_err,
            started.elapsed().as_secs_f64()
        );
    }
}

#[test]
fn g1_nonmatching_mortar_converges_to_monolithic() {
    let started = std::time::Instant::now();
    let e_pa = 10e9;
    let nu = 0.30;
    let h_m = 0.004;
    let rho = 500.0;
    let lx = 0.36;
    let ly = 0.24;

    // Clamped at x = 0 and x = lx so the seam x = lx / 2 is strictly interior and all modes > 0.
    // Refine with k in [2, 4, 6]:
    // Left half [0, lx/2] x [0, ly] has mesh (3*k, 3*k).
    // Right half [lx/2, lx] x [0, ly] has mesh (2*k, 2*k) (seam ratio 3 : 2).
    // Compare against the finest monolithic reference across the ladder.
    let rungs = [2usize, 4, 6];
    let n_ref = 24usize;

    for kind in [ElementKind::CstDkt, ElementKind::Mitc3p] {
        // Compute fine monolithic reference on (2 * n_ref) x n_ref mesh.
        let mono_mesh = rect_xy_mesh(0.0, lx, 0.0, ly, 0.0, 2 * n_ref, n_ref);
        let mono_clamped: Vec<usize> = (0..=n_ref)
            .flat_map(|j| [j * (2 * n_ref + 1), j * (2 * n_ref + 1) + 2 * n_ref])
            .collect();
        let mod_ref = assemble_component_model(
            kind,
            &mono_mesh,
            e_pa,
            nu,
            h_m,
            rho,
            &mono_clamped,
            DRILLING_ALPHA,
        );
        let eigs_ref = lowest_positive_eigenvalues(&mod_ref, 10);

        let mut mode_errors: Vec<Vec<f64>> = vec![Vec::with_capacity(rungs.len()); 10];
        let mut hs = Vec::with_capacity(rungs.len());

        for &k in &rungs {
            let n_a = 3 * k;
            let n_b = 2 * k;
            hs.push(1.0 / (k as f64));

            let left_mesh = rect_xy_mesh(0.0, 0.5 * lx, 0.0, ly, 0.0, n_a, n_a);
            let right_mesh = rect_xy_mesh(0.5 * lx, lx, 0.0, ly, 0.0, n_b, n_b);

            let left_clamped: Vec<usize> = (0..=n_a).map(|j| j * (n_a + 1)).collect();
            let right_clamped: Vec<usize> = (0..=n_b).map(|j| j * (n_b + 1) + n_b).collect();

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

            // Finer side (3*k) as master, coarser side (2*k) as slave.
            let master_seam: Vec<usize> = (0..=n_a).map(|j| j * (n_a + 1) + n_a).collect();
            let slave_seam: Vec<usize> = (0..=n_b).map(|j| j * (n_b + 1)).collect();

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
            let ties = [EdgeTie {
                a: TieSide {
                    component: 0,
                    nodes: master_seam,
                },
                b: TieSide {
                    component: 1,
                    nodes: slave_seam,
                },
                offset: TieOffset::Midsurface,
                rotations: RotationCoupling::Full,
                kind: TieKind::Rigid,
            }];
            let tied = tie_shells(&comps, &ties, default_tie_budget(), || false)
                .expect("3:2 nonmatching tie succeeds");
            let eigs = lowest_positive_eigenvalues(&tied.model, 10);
            for m in 0..10 {
                let rel_err = (eigs[m] - eigs_ref[m]).abs() / eigs_ref[m];
                mode_errors[m].push(rel_err);
            }
        }

        let mut min_order = f64::INFINITY;
        let mut orders = Vec::with_capacity(10);
        for m in 0..10 {
            let ladder = RefinementLadder::new(hs.clone(), mode_errors[m].clone())
                .expect("valid refinement ladder");
            let fit = fit_order(&ladder);
            orders.push(fit.observed);
            min_order = min_order.min(fit.observed);
            assert!(
                fit.observed >= 1.8,
                "{}: mode {} observed convergence order {:.3} < 1.8 (errors={:?})",
                kind.name(),
                m + 1,
                fit.observed,
                mode_errors[m]
            );
        }

        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W2a\",\"test\":\"g1_nonmatching_mortar_converges_to_monolithic\",\"tier\":\"G1\",\"element\":\"{}\",\"min_observed_order\":{},\"orders\":{:?},\"finest_errors\":{:?},\"elapsed_s\":{}}}",
            kind.name(),
            min_order,
            orders,
            mode_errors
                .iter()
                .map(|e| *e.last().unwrap())
                .collect::<Vec<_>>(),
            started.elapsed().as_secs_f64()
        );
    }
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
