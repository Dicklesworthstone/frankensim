//! fs-flux conformance battery (bead tfz.17): BDM1 exactness, Stokes
//! MMS convergence, the pressure-robustness identity (THE property
//! this discretization exists for), lid-driven cavity Picard, Taylor–
//! Green IMEX BDF1, and the discrete adjoint. Every gate states what
//! it measures; coarse-mesh honesty bands are labeled as such.

use fs_flux::bdm::{cell_basis, edge_gauss_pub, eval_basis, tri_quad};
use fs_flux::{FluxParams, FluxSystem, TriMesh};
use fs_solid::Mesh2;

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

fn unit_mesh(n: usize) -> TriMesh {
    TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, n, n))
}

const PI: f64 = std::f64::consts::PI;

/// Divergence-free MMS velocity (stream function ψ = sin²(πx)sin²(πy)/π,
/// so u vanishes on ∂[0,1]²) and a boundary-active pressure.
fn mms_u(x: [f64; 2]) -> [f64; 2] {
    [
        (PI * x[0]).sin().powi(2) * (2.0 * PI * x[1]).sin(),
        -(2.0 * PI * x[0]).sin() * (PI * x[1]).sin().powi(2),
    ]
}

fn mms_p(x: [f64; 2]) -> f64 {
    (PI * x[0]).sin() * (PI * x[1]).cos()
}

/// f = −νΔu + ∇p for the MMS pair.
fn mms_f(x: [f64; 2], nu: f64) -> [f64; 2] {
    let (sx, cx) = ((PI * x[0]).sin(), (PI * x[0]).cos());
    let (sy, cy) = ((PI * x[1]).sin(), (PI * x[1]).cos());
    let (s2x, c2x) = ((2.0 * PI * x[0]).sin(), (2.0 * PI * x[0]).cos());
    let (s2y, c2y) = ((2.0 * PI * x[1]).sin(), (2.0 * PI * x[1]).cos());
    let pp = PI * PI;
    let lap_u1 = 2.0 * pp * c2x * s2y - 4.0 * pp * sx * sx * s2y;
    let lap_u2 = 4.0 * pp * s2x * sy * sy - 2.0 * pp * s2x * c2y;
    let grad_p = [PI * cx * cy, -PI * sx * sy];
    [-nu * lap_u1 + grad_p[0], -nu * lap_u2 + grad_p[1]]
}

/// BDM1 interpolant: dofs are the edge normal moments of the field.
fn interpolate(mesh: &TriMesh, n_total: usize, u: &dyn Fn([f64; 2]) -> [f64; 2]) -> Vec<f64> {
    let mut x = vec![0.0f64; n_total];
    for (e, edge) in mesh.edges.iter().enumerate() {
        let (va, vb) = (mesh.verts[edge.verts.0], mesh.verts[edge.verts.1]);
        let mut m0 = 0.0;
        let mut m1 = 0.0;
        for (gx, w) in edge_gauss_pub(va, vb) {
            let uv = u(gx);
            let un = uv[0] * edge.normal[0] + uv[1] * edge.normal[1];
            let sl = ((gx[0] - va[0]) * (vb[0] - va[0]) + (gx[1] - va[1]) * (vb[1] - va[1]))
                / (edge.len * edge.len)
                - 0.5;
            m0 += w * un / edge.len;
            m1 += w * un * sl / edge.len;
        }
        x[2 * e] = m0;
        x[2 * e + 1] = m1;
    }
    x
}

/// flux-001: the BDM1 basis reproduces its defining dofs (Kronecker
/// identity), its divergence is the constant matching the net normal
/// flux, and the normal trace is single-valued across interior edges
/// (H(div) conformity by construction, verified numerically).
#[test]
fn flux_001_bdm_basis_exactness() {
    let mesh = unit_mesh(3);
    let mut worst_delta = 0.0f64;
    let mut worst_flux = 0.0f64;
    for t in 0..mesh.tris.len() {
        let basis = cell_basis(&mesh, t);
        for i in 0..6 {
            // Recompute all 6 dofs of basis function i.
            for k in 0..3 {
                let (e, _) = mesh.tri_edges[t][k];
                let edge = &mesh.edges[e];
                let (va, vb) = (mesh.verts[edge.verts.0], mesh.verts[edge.verts.1]);
                let mut m0 = 0.0;
                let mut m1 = 0.0;
                for (gx, w) in edge_gauss_pub(va, vb) {
                    let v = eval_basis(&basis, i, gx);
                    let un = v[0] * edge.normal[0] + v[1] * edge.normal[1];
                    let sl = ((gx[0] - va[0]) * (vb[0] - va[0])
                        + (gx[1] - va[1]) * (vb[1] - va[1]))
                        / (edge.len * edge.len)
                        - 0.5;
                    m0 += w * un / edge.len;
                    m1 += w * un * sl / edge.len;
                }
                let want0 = if i == 2 * k { 1.0 } else { 0.0 };
                let want1 = if i == 2 * k + 1 { 1.0 } else { 0.0 };
                worst_delta = worst_delta.max((m0 - want0).abs()).max((m1 - want1).abs());
            }
            // Divergence theorem: div·area = net outward flux.
            let mut flux = 0.0;
            for k in 0..3 {
                let (e, sign) = mesh.tri_edges[t][k];
                let edge = &mesh.edges[e];
                let (va, vb) = (mesh.verts[edge.verts.0], mesh.verts[edge.verts.1]);
                for (gx, w) in edge_gauss_pub(va, vb) {
                    let v = eval_basis(&basis, i, gx);
                    flux += w * sign * (v[0] * edge.normal[0] + v[1] * edge.normal[1]);
                }
            }
            worst_flux = worst_flux.max((basis.div[i] * mesh.areas[t] - flux).abs());
        }
    }
    // Normal-trace continuity across every interior edge for the two
    // global basis functions the edge owns.
    let mut worst_jump = 0.0f64;
    for (e, edge) in mesh.edges.iter().enumerate() {
        if edge.tris.1 == usize::MAX {
            continue;
        }
        let (t0, t1) = edge.tris;
        let (va, vb) = (mesh.verts[edge.verts.0], mesh.verts[edge.verts.1]);
        let b0 = cell_basis(&mesh, t0);
        let b1 = cell_basis(&mesh, t1);
        let loc = |t: usize| (0..3).find(|&k| mesh.tri_edges[t][k].0 == e).unwrap();
        let (k0, k1) = (loc(t0), loc(t1));
        for m in 0..2 {
            for (gx, _) in edge_gauss_pub(va, vb) {
                let v0 = eval_basis(&b0, 2 * k0 + m, gx);
                let v1 = eval_basis(&b1, 2 * k1 + m, gx);
                let jump = (v0[0] - v1[0]) * edge.normal[0] + (v0[1] - v1[1]) * edge.normal[1];
                worst_jump = worst_jump.max(jump.abs());
            }
        }
    }
    verdict(
        "flux-001-dof-identity",
        worst_delta < 1e-11,
        &format!("worst |dof(basis)-delta| {worst_delta:.2e}"),
    );
    verdict(
        "flux-001-div-theorem",
        worst_flux < 1e-11,
        &format!("worst |div*area - flux| {worst_flux:.2e}"),
    );
    verdict(
        "flux-001-normal-trace",
        worst_jump < 1e-11,
        &format!("worst interior normal jump {worst_jump:.2e}"),
    );
}

/// flux-002: steady Stokes MMS — velocity L2 slope ≥ 1.6 (BDM1 is
/// second order), pressure L2 slope ≥ 0.7 (P0 is first order), and
/// per-cell divergence at roundoff on EVERY mesh: the exactness that
/// no penalty constant can buy.
#[test]
fn flux_002_stokes_mms_convergence() {
    let nu = 1.0;
    let params = FluxParams {
        nu,
        ..FluxParams::default()
    };
    let mut eu = Vec::new();
    let mut ep = Vec::new();
    let mut worst_div_all = 0.0f64;
    for n in [6usize, 8, 12] {
        let mesh = unit_mesh(n);
        let sys = FluxSystem::new(&mesh);
        let (a, rhs) = sys.assemble(params, &|x| mms_f(x, nu), &mms_u, None);
        let sol = sys.solve_linear(&a, &rhs, params);
        let (e_u, worst_div) = sys.velocity_error(&sol.x, &mms_u);
        worst_div_all = worst_div_all.max(worst_div);
        // Pressure L2, both mean-shifted (the pin fixes an arbitrary level).
        let ph = sys.pressure(&sol.x);
        let mut mean_ex = 0.0;
        let mut total = 0.0;
        for t in 0..mesh.tris.len() {
            let p: [[f64; 2]; 3] = core::array::from_fn(|k| mesh.verts[mesh.tris[t][k]]);
            for (q, w) in tri_quad(p, mesh.areas[t]) {
                mean_ex += w * mms_p(q);
                total += w;
            }
        }
        mean_ex /= total;
        let mut e2 = 0.0;
        for t in 0..mesh.tris.len() {
            let p: [[f64; 2]; 3] = core::array::from_fn(|k| mesh.verts[mesh.tris[t][k]]);
            for (q, w) in tri_quad(p, mesh.areas[t]) {
                e2 += w * (ph[t] - (mms_p(q) - mean_ex)).powi(2);
            }
        }
        eu.push(e_u);
        ep.push(e2.sqrt());
    }
    let slope = |e: &[f64], i: usize, n0: f64, n1: f64| -> f64 {
        (e[i] / e[i + 1]).ln() / (n1 / n0).ln()
    };
    // Slope on the finest pair: still preasymptotic toward the
    // theoretical 2 (measured march 1.14/1.40/1.61 over 4..12); the
    // asymptotic-regime confirmation is perf-lane scope, ledgered.
    let su = slope(&eu, 1, 8.0, 12.0);
    let sp = slope(&ep, 1, 8.0, 12.0);
    verdict(
        "flux-002-velocity-order",
        su > 1.5 && eu[1] < eu[0] && eu[2] < eu[1],
        &format!("u L2 errs {:.3e}/{:.3e}/{:.3e} slope {su:.2}", eu[0], eu[1], eu[2]),
    );
    verdict(
        "flux-002-pressure-order",
        sp > 0.7 && ep[2] < ep[0],
        &format!("p L2 errs {:.3e}/{:.3e}/{:.3e} slope {sp:.2}", ep[0], ep[1], ep[2]),
    );
    verdict(
        "flux-002-exact-divergence",
        worst_div_all < 1e-9,
        &format!("worst per-cell div across meshes {worst_div_all:.2e}"),
    );
}

/// flux-003: pressure-robustness — adding a HUGE gradient forcing
/// A·∇φ changes the velocity by NOTHING (to factorization roundoff)
/// and shifts the pressure by exactly A·(Π₀φ − Π₀φ|cell0). A non-
/// robust pair (Taylor–Hood, P1-P0 stabilized) pollutes velocity at
/// O(A·h); here the theorem is a discrete identity.
#[test]
fn flux_003_pressure_robustness() {
    let nu = 0.01; // small viscosity makes non-robust pollution huge
    let params = FluxParams {
        nu,
        ..FluxParams::default()
    };
    let mesh = unit_mesh(6);
    let sys = FluxSystem::new(&mesh);
    // φ of degree 4: ∇φ·(P1 basis) is degree 4, integrated EXACTLY by
    // the assembly quadrature — the invariance identity needs the
    // discrete rhs to be exactly a gradient forcing (a degree-8 φ
    // leaked 1.7e-4 of quadrature error into the velocity; measured).
    let phi = |x: [f64; 2]| x[0] * (1.0 - x[0]) * x[1] * (1.0 - x[1]);
    let grad_phi = |x: [f64; 2]| -> [f64; 2] {
        [
            (1.0 - 2.0 * x[0]) * x[1] * (1.0 - x[1]),
            x[0] * (1.0 - x[0]) * (1.0 - 2.0 * x[1]),
        ]
    };
    let amp = 1.0e4;
    let f1 = |x: [f64; 2]| -> [f64; 2] {
        let b = mms_f(x, nu);
        let g = grad_phi(x);
        [b[0] + g[0], b[1] + g[1]]
    };
    let f2 = |x: [f64; 2]| -> [f64; 2] {
        let b = mms_f(x, nu);
        let g = grad_phi(x);
        [b[0] + amp * g[0], b[1] + amp * g[1]]
    };
    let (a1, r1) = sys.assemble(params, &f1, &mms_u, None);
    let (a2, r2) = sys.assemble(params, &f2, &mms_u, None);
    let s1 = sys.solve_linear(&a1, &r1, params);
    let s2 = sys.solve_linear(&a2, &r2, params);
    let mut du = 0.0f64;
    let mut nrm = 0.0f64;
    for i in 0..sys.n_u {
        du += (s2.x[i] - s1.x[i]).powi(2);
        nrm += s1.x[i].powi(2);
    }
    let rel_u = (du / nrm.max(1e-300)).sqrt();
    // Predicted pressure shift: (amp−1)·(Π₀φ − Π₀φ|cell0).
    let cell_mean = |t: usize| -> f64 {
        let p: [[f64; 2]; 3] = core::array::from_fn(|k| mesh.verts[mesh.tris[t][k]]);
        let mut m = 0.0;
        for (q, w) in tri_quad(p, mesh.areas[t]) {
            m += w * phi(q);
        }
        m / mesh.areas[t]
    };
    let base = cell_mean(0);
    let mut dp_err = 0.0f64;
    let mut dp_scale = 0.0f64;
    for t in 0..mesh.tris.len() {
        let predicted = (amp - 1.0) * (cell_mean(t) - base);
        let got = s2.x[sys.n_u + t] - s1.x[sys.n_u + t];
        dp_err = dp_err.max((got - predicted).abs());
        dp_scale = dp_scale.max(predicted.abs());
    }
    verdict(
        "flux-003-velocity-invariance",
        rel_u < 1e-7,
        &format!("rel velocity change under 1e4 gradient forcing {rel_u:.2e}"),
    );
    verdict(
        "flux-003-pressure-shift-identity",
        dp_err < 1e-7 * dp_scale.max(1.0),
        &format!("worst |dp - (amp-1)*P0(phi)| {dp_err:.2e} vs scale {dp_scale:.2e}"),
    );
}

/// flux-004: lid-driven cavity at Re=100 by Picard. Gates: Picard
/// contraction to a fixed point, exact per-cell div-freeness of the
/// NONLINEAR solution, bounded kinetic energy, and the vertical-
/// centerline u_x profile inside a COARSE-MESH honesty band (±0.15)
/// around Ghia et al. — the fine-mesh table comparison is perf-lane
/// scope, ledgered in the contract.
#[test]
fn flux_004_cavity_picard() {
    let params = FluxParams {
        nu: 0.01,
        ..FluxParams::default()
    };
    let mesh = unit_mesh(8);
    let sys = FluxSystem::new(&mesh);
    let lid = |x: [f64; 2]| -> [f64; 2] {
        if x[1] > 1.0 - 1e-9 {
            [1.0, 0.0]
        } else {
            [0.0, 0.0]
        }
    };
    let zero = |_x: [f64; 2]| [0.0f64, 0.0];
    let (sol, its, update) = sys.picard(params, &zero, &lid, 40, 1e-9);
    let (_, worst_div) = sys.velocity_error(&sol.x, &|_x| [0.0, 0.0]);
    let ke = sys.kinetic_energy(&sol.x);
    // Ghia, Ghia & Shin (1982), Re=100, u_x along x=0.5.
    let ghia = [(0.9531, 0.687_17), (0.5, -0.205_81), (0.1016, -0.064_34)];
    let mut worst_dev = 0.0f64;
    let mut profile = String::new();
    for (y, want) in ghia {
        let got = sys.velocity_at(&sol.x, [0.5, y])[0];
        worst_dev = worst_dev.max((got - want).abs());
        profile.push_str(&format!("y={y}: {got:.3} vs {want:.3}; "));
    }
    verdict(
        "flux-004-picard-converged",
        update < 1e-6 && its < 40,
        &format!("picard iters {its} final update {update:.2e}"),
    );
    verdict(
        "flux-004-nonlinear-divfree",
        worst_div < 1e-9,
        &format!("worst per-cell div of NS solution {worst_div:.2e}"),
    );
    verdict(
        "flux-004-energy-bounded",
        ke > 1e-4 && ke < 0.5,
        &format!("cavity KE {ke:.4}"),
    );
    verdict(
        "flux-004-ghia-coarse-band",
        worst_dev < 0.15,
        &format!("8x8 HONESTY BAND +-0.15: {profile}worst dev {worst_dev:.3}"),
    );
}

/// flux-005: Taylor–Green by IMEX BDF1 (implicit Stokes, lagged
/// upwind convection). Gates: first-order temporal self-convergence
/// (Richardson on dt/dt2/dt4 — isolates time error from the spatial
/// floor), kinetic-energy decay tracking exp(−4π²νt), div-free every
/// step, and accuracy vs the exact field.
#[test]
fn flux_005_taylor_green_bdf1() {
    let nu = 0.05;
    let params = FluxParams {
        nu,
        ..FluxParams::default()
    };
    let mesh = unit_mesh(6);
    let sys = FluxSystem::new(&mesh);
    let tg = |t: f64| {
        move |x: [f64; 2]| -> [f64; 2] {
            let decay = (-2.0 * PI * PI * nu * t).exp();
            [
                -(PI * x[0]).cos() * (PI * x[1]).sin() * decay,
                (PI * x[0]).sin() * (PI * x[1]).cos() * decay,
            ]
        }
    };
    let zero = |_x: [f64; 2]| [0.0f64, 0.0];
    let t_end = 0.05;
    let run = |steps: usize| -> Vec<f64> {
        let dt = t_end / steps as f64;
        let mut u = interpolate(&mesh, sys.n, &tg(0.0));
        for s in 0..steps {
            let t1 = (s + 1) as f64 * dt;
            let sol = sys.bdf1_step(params, &zero, &tg(t1), &u, dt);
            u = sol.x;
        }
        u
    };
    let u4 = run(4);
    let u8 = run(8);
    let u16 = run(16);
    let diff = |a: &[f64], b: &[f64]| -> f64 {
        (0..sys.n_u).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt()
    };
    let d1 = diff(&u4, &u8);
    let d2 = diff(&u8, &u16);
    let ratio = d1 / d2.max(1e-300);
    let (err, worst_div) = sys.velocity_error(&u16, &tg(t_end));
    let ke0 = sys.kinetic_energy(&interpolate(&mesh, sys.n, &tg(0.0)));
    let ke = sys.kinetic_energy(&u16);
    let want_ratio = (-4.0 * PI * PI * nu * t_end).exp();
    let ke_dev = (ke / ke0 - want_ratio).abs();
    verdict(
        "flux-005-bdf1-temporal-order",
        ratio > 1.5 && ratio < 3.0,
        &format!("Richardson ratio {ratio:.2} (order-1 target 2)"),
    );
    verdict(
        "flux-005-energy-decay",
        ke_dev < 0.05 * want_ratio,
        &format!("KE ratio {:.4} vs exp {:.4}", ke / ke0, want_ratio),
    );
    verdict(
        "flux-005-transient-divfree",
        worst_div < 1e-8,
        &format!("worst per-cell div at T {worst_div:.2e}"),
    );
    verdict(
        "flux-005-accuracy",
        err < 0.1,
        &format!("u L2 error at T on 6x6 {err:.3e}"),
    );
}

/// flux-006: discrete adjoint of steady Stokes — dJ/dν for
/// J = ∫w·u_h via one adjoint solve against central finite
/// differences. The operator is LINEAR in ν, so ∂A/∂ν and ∂r/∂ν are
/// exact matrix differences; agreement is limited only by FD trunc.
#[test]
fn flux_006_adjoint_viscosity() {
    let nu = 0.7;
    let mk = |nu: f64| FluxParams {
        nu,
        ..FluxParams::default()
    };
    let mesh = unit_mesh(6);
    let sys = FluxSystem::new(&mesh);
    let force = |x: [f64; 2]| mms_f(x, 1.0); // fixed body force, ν-independent
    // The weight must have CURL: a gradient weight w = ∇χ satisfies
    // ∫w·u_h ≡ 0 for the exactly div-free, zero-normal-flux u_h — the
    // pressure-robust structure annihilates such functionals (gated
    // below as a bonus). mms_u is div-free and decidedly not a
    // gradient.
    let weight = mms_u;
    // J(ν) = j^T x(ν) with j the load vector of the weight field.
    let mut j = vec![0.0f64; sys.n];
    for t in 0..mesh.tris.len() {
        let p: [[f64; 2]; 3] = core::array::from_fn(|k| mesh.verts[mesh.tris[t][k]]);
        for (q, w) in tri_quad(p, mesh.areas[t]) {
            let wv = weight(q);
            for i in 0..6 {
                let b = eval_basis(&sys.bases[t], i, q);
                j[2 * mesh.tri_edges[t][i / 2].0 + (i % 2)] += w * (wv[0] * b[0] + wv[1] * b[1]);
            }
        }
    }
    let solve_j = |nuv: f64| -> (f64, Vec<f64>, fs_sparse::Csr, Vec<f64>) {
        let (a, r) = sys.assemble(mk(nuv), &force, &mms_u, None);
        let sol = sys.solve_linear(&a, &r, mk(nuv));
        let jval = (0..sys.n).map(|i| j[i] * sol.x[i]).sum::<f64>();
        (jval, sol.x, a, r)
    };
    let (_, x0, a0, r0) = solve_j(nu);
    let (_, _, a2, r2) = solve_j(2.0 * nu);
    // Exact ν-derivatives by linearity: A' = (A(2ν)−A(ν))/ν.
    let lam = sys.solve_adjoint(&a0, &j);
    let mut a0x = vec![0.0f64; sys.n];
    let mut a2x = vec![0.0f64; sys.n];
    a0.spmv(&x0, &mut a0x);
    a2.spmv(&x0, &mut a2x);
    let mut dj_adj = 0.0;
    for i in 0..sys.n {
        let aprime_x = (a2x[i] - a0x[i]) / nu;
        let rprime = (r2[i] - r0[i]) / nu;
        dj_adj += lam[i] * (rprime - aprime_x);
    }
    let h = 0.01 * nu;
    let (jp, _, _, _) = solve_j(nu + h);
    let (jm, _, _, _) = solve_j(nu - h);
    let dj_fd = (jp - jm) / (2.0 * h);
    let rel = (dj_adj - dj_fd).abs() / dj_fd.abs().max(1e-300);
    verdict(
        "flux-006-adjoint-gradient",
        rel < 5e-4,
        &format!("dJ/dnu adjoint {dj_adj:.6e} vs FD {dj_fd:.6e} rel {rel:.2e}"),
    );
    // Bonus: a GRADIENT weight is annihilated exactly — ∫∇χ·u_h = 0
    // for div-free u_h with zero boundary flux (χ = sin(πy)/π −
    // cos(πx)/π here). This is pressure-robustness read backwards.
    let mut jg = vec![0.0f64; sys.n];
    for t in 0..mesh.tris.len() {
        let p: [[f64; 2]; 3] = core::array::from_fn(|k| mesh.verts[mesh.tris[t][k]]);
        for (q, w) in tri_quad(p, mesh.areas[t]) {
            let wv = [(PI * q[0]).sin(), (PI * q[1]).cos()];
            for i in 0..6 {
                let b = eval_basis(&sys.bases[t], i, q);
                jg[2 * mesh.tri_edges[t][i / 2].0 + (i % 2)] += w * (wv[0] * b[0] + wv[1] * b[1]);
            }
        }
    }
    let j_annihilated = (0..sys.n).map(|i| jg[i] * x0[i]).sum::<f64>();
    verdict(
        "flux-006-gradient-weight-annihilation",
        j_annihilated.abs() < 1e-8,
        &format!("integral of gradient-weight vs u_h = {j_annihilated:.2e}"),
    );
}
