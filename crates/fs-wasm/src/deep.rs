//! fs-wasm · DEEP tier — the newly-unlocked upper-stack kernels in the browser.
//!
//! Every function here runs the *real* compiled FrankenSim kernel — the same
//! code the native workspace builds — targeted at `wasm32-unknown-unknown`
//! through asupersync's `wasm-browser-prod` surface. No mocks, no re-derived
//! math: the exterior-calculus Hodge decomposition, the FEEC-native transient
//! Navier–Stokes solver, the exact Gaussian-process Cholesky inference, the
//! information-geometric CMA-ES, log-domain Sinkhorn OT, DFT-block-diagonalized
//! circulant solves, the resumable Krylov stack (CG/MINRES/GMRES), the
//! interval-certified CutFEM quadtree, and Bernstein free-form deformation are
//! all invoked directly.
//!
//! SAFETY CONTRACT (identical to the rest of the crate): `unsafe_code` is
//! forbidden, every input is clamped to a safe range, every fallible kernel
//! result is folded to `NaN` / an empty vector, and every loop is capped.
//! Nothing here can trap — a wasm trap would kill the whole page.

use fs_bo::acq::expected_improvement;
use fs_bo::gp::{Gp, Kernel, Matern};
use fs_cutfem::{Circle, CutSdf, Quadtree, cut_cell_rules};
use fs_dfo::cma::{CmaParams, cmaes};
use fs_dfo::ot::{cost_sq_1d, sinkhorn};
use fs_feec::{
    betti_numbers, deram1, element_geometry, harmonic_basis, hodge_decompose,
    hodge_diagonal_barycentric, masked_cube_grid,
};
use fs_flux::{FluxParams, FluxSystem, TriMesh};
use fs_rep_mesh::TetComplex;
use fs_solid::Mesh2;
use fs_solver::krylov::{CgState, GmresState, MinresState};
use fs_solver::op::CsrOp;
use fs_sparse::precond::IdentityPrecond;
use fs_symmetry::{CyclicGroup, solve_circulant};
use fs_xform::{FfdLattice, Parameterization, Point3, Vec3, det3};

/* ======================================================================= */
/*  1 · Hodge decomposition (fs-feec)                                       */
/* ======================================================================= */

/// One of three deterministic multiply-connected slab fixtures, built with
/// the real `masked_cube_grid` Kuhn/Freudenthal tet subdivision:
/// `shape` 0 = solid 3×3 disk (b₁=0), 1 = 3×3 annulus / one hole (b₁=1),
/// 2 = 5×3 domain with two holes (b₁=2). One cell layer thick in z.
fn build_shape(shape: u32) -> (TetComplex, Vec<[f64; 3]>) {
    match shape {
        0 => masked_cube_grid(3, 3, 1, &|_i, _j, _k| true),
        1 => masked_cube_grid(3, 3, 1, &|i, j, _k| !(i == 1 && j == 1)),
        _ => masked_cube_grid(5, 3, 1, &|i, j, _k| !((i == 1 || i == 3) && j == 1)),
    }
}

/// The center of the (first) hole for the seeding vector field — the axis the
/// harmonic circulation wraps.
fn shape_hole_center(shape: u32) -> (f64, f64) {
    match shape {
        0 | 1 => (1.5, 1.5),
        _ => (1.5, 1.5),
    }
}

/// Discrete **Hodge decomposition** of a real vector field on a multiply-
/// connected 2D slab, via `fs-feec` exterior calculus: the field is de-Rham-
/// mapped to a 1-cochain (`deram1`, edge line integrals), then split by
/// `hodge_decompose` into EXACT (gradient, `im d`), COEXACT (curl, `im δ`) and
/// HARMONIC (cohomology) parts against the diagonal Hodge stars. The harmonic
/// dimension equals the first Betti number — geometry (integer-rank Betti) and
/// physics (the projection) agreeing is the internal consistency proof.
///
/// `shape`: 0 = disk, 1 = annulus (one hole), 2 = two holes.
///
/// Output layout (length `12 + 6·E`, `E` = edge count):
/// - `[0]`        — `shape` echoed.
/// - `[1..5]`     — Betti numbers `b₀, b₁, b₂, b₃` (exact integer rank).
/// - `[5..8]`     — M-weighted orthogonality residuals of the split,
///   `(exact·coexact, exact·harmonic, coexact·harmonic)` relative to `‖x‖²_M`
///   (all ≈ 1e-12 — the decomposition is orthogonal to machine precision).
/// - `[8..11]`    — component energies `[‖exact‖²_M, ‖coexact‖²_M, ‖harmonic‖²_M]`
///   (the Hodge energy split; sums to the field's total M-energy).
/// - `[11]`       — `E`, the edge count.
/// - `[12 .. 12+6E]` — per edge `e` (canonical sorted order), six values:
///   `[mx, my, mz, exact_e, coexact_e, harmonic_e]`, where `(mx,my,mz)` is the
///   edge midpoint and the three scalars are that edge's component of the
///   decomposed 1-cochain (units: circulation along the edge).
pub fn hodge_decomposition(shape: u32) -> Vec<f64> {
    let s = shape.min(2);
    let (complex, positions) = build_shape(s);
    let geo = element_geometry(&complex, &positions);
    let (cx, cy) = shape_hole_center(s);
    // A field with all three ingredients: a 1/r circulation around the hole
    // (the harmonic seed), plus a radial gradient part.
    let field = |p: [f64; 3]| -> [f64; 3] {
        let dx = p[0] - cx;
        let dy = p[1] - cy;
        let r2 = (dx * dx + dy * dy).max(0.05);
        [-dy / r2 + 0.3 * dx, dx / r2 + 0.3 * dy, 0.0]
    };
    let x = deram1(&complex, &positions, &field);
    let parts = hodge_decompose(&complex, &positions, &geo, 1, &x);
    let m = hodge_diagonal_barycentric(&complex, &positions, &geo, 1);
    let betti = betti_numbers(&complex);
    let energy = |v: &[f64]| -> f64 { m.iter().zip(v).map(|(w, a)| w * a * a).sum() };

    let e = complex.edges.len();
    let mut out = Vec::with_capacity(12 + 6 * e);
    out.push(s as f64);
    for b in betti {
        out.push(b as f64);
    }
    out.extend_from_slice(&parts.ortho_residuals);
    out.push(energy(&parts.exact));
    out.push(energy(&parts.coexact));
    out.push(energy(&parts.harmonic));
    out.push(e as f64);
    for (idx, &[u, v]) in complex.edges.iter().enumerate() {
        let pu = positions[u as usize];
        let pv = positions[v as usize];
        out.push(f64::midpoint(pu[0], pv[0]));
        out.push(f64::midpoint(pu[1], pv[1]));
        out.push(f64::midpoint(pu[2], pv[2]));
        out.push(parts.exact.get(idx).copied().unwrap_or(f64::NAN));
        out.push(parts.coexact.get(idx).copied().unwrap_or(f64::NAN));
        out.push(parts.harmonic.get(idx).copied().unwrap_or(f64::NAN));
    }
    out
}

/* ======================================================================= */
/*  10 · Betti numbers + harmonic representatives (fs-feec)                 */
/* ======================================================================= */

/// **Betti numbers + harmonic representative fields** for a chosen shape, via
/// `fs-feec`: `betti_numbers` (exact fraction-free integer rank of the incidence
/// operators — no tolerance knobs) and `harmonic_basis` (an M-orthonormal basis
/// of the degree-1 harmonic space, one representative 1-cochain per hole).
///
/// `shape`: 0 = disk (b₁=0), 1 = annulus (b₁=1), 2 = two holes (b₁=2).
///
/// Output layout (length `6 + 3E + H·E`):
/// - `[0..4]`     — Betti numbers `b₀, b₁, b₂, b₃`.
/// - `[4]`        — `E`, the edge count.
/// - `[5 .. 5+3E]`— edge midpoints `[mx, my, mz]` per edge (canonical order).
/// - `[5+3E]`     — `H`, the number of harmonic representatives found (= `b₁`).
/// - then `H` blocks of `E` values (starting at `[6+3E]`): harmonic
///   representative `h`'s value on each
///   edge (an M-orthonormal harmonic 1-cochain — the circulation loop around
///   hole `h`).
pub fn betti_shapes(shape: u32) -> Vec<f64> {
    let s = shape.min(2);
    let (complex, positions) = build_shape(s);
    let geo = element_geometry(&complex, &positions);
    let betti = betti_numbers(&complex);
    let b1 = betti[1];
    let basis = harmonic_basis(&complex, &positions, &geo, 1, b1 + 2);
    let e = complex.edges.len();

    let mut out = Vec::with_capacity(6 + 3 * e + 1 + basis.len() * e);
    for b in betti {
        out.push(b as f64);
    }
    out.push(e as f64);
    for &[u, v] in &complex.edges {
        let pu = positions[u as usize];
        let pv = positions[v as usize];
        out.push(f64::midpoint(pu[0], pv[0]));
        out.push(f64::midpoint(pu[1], pv[1]));
        out.push(f64::midpoint(pu[2], pv[2]));
    }
    out.push(basis.len() as f64);
    for h in &basis {
        for k in 0..e {
            out.push(h.get(k).copied().unwrap_or(0.0));
        }
    }
    out
}

/* ======================================================================= */
/*  2 · Navier–Stokes transient — lid-driven cavity (fs-flux)               */
/* ======================================================================= */

/// **Transient incompressible Navier–Stokes** on a lid-driven unit-square
/// cavity, via `fs-flux`: FEEC-native BDM1–P0 (exactly divergence-free,
/// pressure-robust) on a triangulated mesh (`Mesh2::triangles` → `TriMesh`),
/// time-stepped with the real IMEX `bdf1_step` (implicit Stokes + upwinded DG
/// convection lagged on the single-valued H(div) face flux), solved through the
/// exact dense saddle factorization. The top lid moves at u=(1,0); the other
/// three walls are no-slip.
///
/// `cells` (per side, clamped 3..=8), `frames` (1..=30), `re` (Reynolds,
/// 10..=500 → viscosity ν = 1/Re), `steps_per_frame` (1..=6, BDF1 steps of
/// dt=0.03 between recorded frames). Grid sample resolution is fixed at G=20.
///
/// Output layout (length `2 + F·2·G·G`):
/// - `[0]` — `G` (the sample-grid side; 20).
/// - `[1]` — `F` (frame count).
/// - then `F` frame blocks, each `2·G·G` values, row-major over the G×G grid
///   covering `(0,1)²` (index `iy·G + ix`, x fastest):
///   - first `G·G` — velocity magnitude `|u|` at each grid point.
///   - next `G·G`  — vorticity `ω = ∂v/∂x − ∂u/∂y` (central differences).
pub fn navier_stokes_cavity(
    cells_in: usize,
    frames_in: usize,
    re: f64,
    steps_per_frame_in: usize,
) -> Vec<f64> {
    let cells = cells_in.clamp(3, 8);
    let frames = frames_in.clamp(1, 30);
    let spf = steps_per_frame_in.clamp(1, 6);
    let reynolds = re.clamp(10.0, 500.0);
    let grid = 20usize;
    let dt = 0.03;

    let mesh2 = Mesh2::triangles(1.0, 1.0, cells, cells);
    let mesh = TriMesh::from_mesh2(&mesh2);
    let sys = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 1.0 / reynolds,
        sigma: 10.0,
        tol: 1e-10,
        max_iters: 40_000,
    };
    // Lid-driven cavity: u=(1,0) on the top wall (y≈1), no-slip elsewhere.
    let g = |x: [f64; 2]| -> [f64; 2] {
        if x[1] > 1.0 - 1e-9 {
            [1.0, 0.0]
        } else {
            [0.0, 0.0]
        }
    };
    let f = |_x: [f64; 2]| -> [f64; 2] { [0.0, 0.0] };

    let mut uprev = vec![0.0f64; sys.n];
    let mut out = Vec::with_capacity(2 + frames * 2 * grid * grid);
    out.push(grid as f64);
    out.push(frames as f64);

    let h = 1.0 / (grid as f64 - 1.0);
    for _ in 0..frames {
        for _ in 0..spf {
            let sol = sys.bdf1_step(params, &f, &g, &uprev, dt);
            uprev = sol.x;
        }
        // Sample u and v on the grid.
        let mut uu = vec![0.0f64; grid * grid];
        let mut vv = vec![0.0f64; grid * grid];
        for iy in 0..grid {
            for ix in 0..grid {
                let px = ix as f64 * h;
                let py = iy as f64 * h;
                let vel = sys.velocity_at(&uprev, [px, py]);
                uu[iy * grid + ix] = vel[0];
                vv[iy * grid + ix] = vel[1];
            }
        }
        // Speed, then vorticity by central differences.
        for iy in 0..grid {
            for ix in 0..grid {
                let vx = uu[iy * grid + ix];
                let vy = vv[iy * grid + ix];
                out.push((vx * vx + vy * vy).sqrt());
            }
        }
        for iy in 0..grid {
            for ix in 0..grid {
                let ixp = (ix + 1).min(grid - 1);
                let ixm = ix.saturating_sub(1);
                let iyp = (iy + 1).min(grid - 1);
                let iym = iy.saturating_sub(1);
                let dvdx =
                    (vv[iy * grid + ixp] - vv[iy * grid + ixm]) / (((ixp - ixm).max(1)) as f64 * h);
                let dudy =
                    (uu[iyp * grid + ix] - uu[iym * grid + ix]) / (((iyp - iym).max(1)) as f64 * h);
                out.push(dvdx - dudy);
            }
        }
    }
    out
}

/* ======================================================================= */
/*  3 · Gaussian-process regression + Expected Improvement (fs-bo)          */
/* ======================================================================= */

/// The 1D objective the GP demo regresses / minimizes (a smooth, mildly
/// multimodal target on `[0, 6]`).
fn gp_objective(x: f64) -> f64 {
    (1.8 * x).sin() * (-0.15 * x).exp()
}

/// **Gaussian-process regression + Bayesian-optimization acquisition** via
/// `fs-bo`: an exact Matérn-5/2 GP (Cholesky inference, `Gp::try_fit`) is fitted
/// to a handful of samples of a 1D function, then its posterior mean/variance
/// and the closed-form Expected-Improvement acquisition are evaluated densely.
/// The EI argmax is the next point a BO loop would query.
///
/// `n_train` (training points, clamped 3..=16), `samples` (dense evaluation
/// points, 16..=512).
///
/// Output layout (length `1 + 2·Ntr + 1 + 4·S + 2`):
/// - `[0]` — `Ntr`, the training-point count.
/// - `[1 .. 1+2Ntr]` — training points interleaved `[x, y]`.
/// - `[1+2Ntr]` — `S`, the dense sample count.
/// - then `S` blocks of 4: `[x, mean, variance, ei]` — posterior mean, posterior
///   variance (σ², latent, ≥0), and Expected Improvement at `x` (minimization,
///   margin ξ=0.01; ≥0). Draw the band as `mean ± k·√variance`.
/// - last two values: `[x_next, ei_max]` — the EI argmax (next BO query) and its
///   EI value. On a fit failure the whole tail is `NaN`.
pub fn gp_regression(n_train_in: usize, samples_in: usize) -> Vec<f64> {
    let ntr = n_train_in.clamp(3, 16);
    let samples = samples_in.clamp(16, 512);
    let (lo, hi) = (0.0f64, 6.0f64);

    // Training points spread across the domain (deterministic).
    let xs: Vec<Vec<f64>> = (0..ntr)
        .map(|i| {
            let t = (i as f64 + 0.5) / ntr as f64;
            vec![lo + (hi - lo) * t]
        })
        .collect();
    let ys: Vec<f64> = xs.iter().map(|p| gp_objective(p[0])).collect();
    let f_best = ys.iter().copied().fold(f64::INFINITY, f64::min);

    let mut out = Vec::with_capacity(1 + 2 * ntr + 1 + 4 * samples + 2);
    out.push(ntr as f64);
    for (p, &y) in xs.iter().zip(&ys) {
        out.push(p[0]);
        out.push(y);
    }
    out.push(samples as f64);

    let kernel = Kernel {
        family: Matern::FiveHalves,
        signal: 1.0,
        lengthscales: vec![0.6],
    };
    let Some(gp) = Gp::try_fit(&xs, &ys, kernel, 1e-4) else {
        for _ in 0..(4 * samples + 2) {
            out.push(f64::NAN);
        }
        return out;
    };

    let mut best_x = f64::NAN;
    let mut best_ei = f64::NEG_INFINITY;
    for i in 0..samples {
        let t = i as f64 / (samples as f64 - 1.0);
        let x = lo + (hi - lo) * t;
        let (mean, var) = gp.predict(&[x]);
        let ei = expected_improvement(&gp, &[x], f_best, 0.01);
        out.push(x);
        out.push(mean);
        out.push(var);
        out.push(ei);
        if ei > best_ei {
            best_ei = ei;
            best_x = x;
        }
    }
    out.push(best_x);
    out.push(best_ei.max(0.0));
    out
}

/* ======================================================================= */
/*  4 · CMA-ES search trajectory (fs-dfo)                                   */
/* ======================================================================= */

/// The Himmelblau test function (four equal global minima at f=0).
fn himmelblau(p: &[f64]) -> f64 {
    let x = p[0];
    let y = p[1];
    let a = x * x + y - 11.0;
    let b = x + y * y - 7.0;
    a * a + b * b
}

/// **CMA-ES search trajectory** minimizing the 2D Himmelblau function via the
/// real information-geometric `fs-dfo::cmaes` (natural-gradient covariance +
/// step-size adaptation). Because the public `CmaReport` exposes the best point,
/// best value and step size σ — but NOT the internal covariance matrix — the
/// "search ellipse" is reported as the isotropic step-size radius σ. Each
/// checkpoint reruns `cmaes` (a pure function of the seed) with a one-generation-
/// larger budget, so the checkpoints are consistent prefixes of one true run.
///
/// `seed` (RNG seed), `gens` (number of generation checkpoints, 2..=60).
///
/// Output layout (length `1 + 5·G + 8`):
/// - `[0]` — `G`, the checkpoint count.
/// - then `G` blocks of 5: `[generation, f_best, x_best_x, x_best_y, sigma]` —
///   the best point/value after that many generations and the current step-size
///   radius σ (draw a circle of radius σ centered at the mean/best point).
/// - last 8 values — the four known Himmelblau minima `[x, y]` (reference
///   targets the trajectory converges onto).
pub fn cmaes_trace(seed_in: u32, gens_in: usize) -> Vec<f64> {
    let gens = gens_in.clamp(2, 60);
    let seed = seed_in as u64;
    let lambda = 6usize; // 4 + floor(3·ln 2) for n=2
    let sigma0 = 2.0;
    let x0 = [0.0f64, 0.0f64];

    let mut out = Vec::with_capacity(1 + 5 * gens + 8);
    out.push(gens as f64);
    for gcp in 1..=gens {
        let budget = 1 + gcp * lambda;
        let params = CmaParams::standard(2, sigma0, budget, 1e-12);
        let mut obj = |p: &[f64]| himmelblau(p);
        let rep = cmaes(&mut obj, &x0, &params, seed);
        out.push(rep.generations as f64);
        out.push(rep.f_best);
        out.push(rep.x_best.first().copied().unwrap_or(f64::NAN));
        out.push(rep.x_best.get(1).copied().unwrap_or(f64::NAN));
        out.push(rep.sigma);
    }
    // Known minima.
    for &(x, y) in &[
        (3.0, 2.0),
        (-2.805118, 3.131312),
        (-3.779310, -3.283186),
        (3.584428, -1.848126),
    ] {
        out.push(x);
        out.push(y);
    }
    out
}

/* ======================================================================= */
/*  5 · Optimal transport — log-domain Sinkhorn (fs-dfo)                    */
/* ======================================================================= */

/// **Entropic optimal transport** between two 1D distributions via the real
/// log-domain, logsumexp-stabilized `fs-dfo::sinkhorn`. Source `a` is a single
/// Gaussian bump; target `b` is a two-bump mixture; the squared-distance cost
/// matrix comes from `cost_sq_1d`. Returns the full coupling (transport plan) so
/// the viz can animate mass flowing from `a` to `b`.
///
/// `n` (bins, clamped 8..=48), `epsilon` (entropic regularization, 0.005..=1.0;
/// smaller ⇒ sharper, more transport-like plan).
///
/// Output layout (length `2 + 3·n + n·n + 1`):
/// - `[0]` — `n`.
/// - `[1]` — the achieved marginal residual (max ‖P·1−a‖∞, ‖Pᵀ1−b‖∞; ≈0).
/// - `[2 .. 2+n]`     — bin positions `x_i = i/(n−1)` on `[0,1]`.
/// - `[2+n .. 2+2n]`  — source marginal `a` (sums to 1).
/// - `[2+2n .. 2+3n]` — target marginal `b` (sums to 1).
/// - `[2+3n .. 2+3n+n·n]` — transport plan `P`, row-major `n×n` (`P[i·n+j]` =
///   mass moved from bin `i` to bin `j`; rows sum to `a`, columns to `b`).
/// - last value — transport cost `⟨P, C⟩`.
pub fn optimal_transport(n_in: usize, epsilon: f64) -> Vec<f64> {
    let n = n_in.clamp(8, 48);
    let eps = epsilon.clamp(0.005, 1.0);
    let floor = 1e-6;
    let pos: Vec<f64> = (0..n).map(|i| i as f64 / (n as f64 - 1.0)).collect();
    let bump = |x: f64, c: f64, w: f64| (-((x - c) * (x - c)) / (2.0 * w * w)).exp();

    let raw_a: Vec<f64> = pos.iter().map(|&x| bump(x, 0.30, 0.09) + floor).collect();
    let raw_b: Vec<f64> = pos
        .iter()
        .map(|&x| bump(x, 0.65, 0.07) + 0.8 * bump(x, 0.85, 0.05) + floor)
        .collect();
    let sa: f64 = raw_a.iter().sum();
    let sb: f64 = raw_b.iter().sum();
    let a: Vec<f64> = raw_a.iter().map(|v| v / sa).collect();
    let b: Vec<f64> = raw_b.iter().map(|v| v / sb).collect();
    let c = cost_sq_1d(&pos, &pos);
    let report = sinkhorn(&a, &b, &c, eps, 400);

    let mut out = Vec::with_capacity(2 + 3 * n + n * n + 1);
    out.push(n as f64);
    out.push(report.marginal_residual);
    out.extend_from_slice(&pos);
    out.extend_from_slice(&a);
    out.extend_from_slice(&b);
    out.extend_from_slice(&report.plan);
    out.push(report.cost);
    out
}

/* ======================================================================= */
/*  6 · Cyclic symmetry — circulant block-diagonalization (fs-symmetry)     */
/* ======================================================================= */

/// **Cyclic-symmetry solve** of an N-fold symmetric ring structure via
/// `fs-symmetry`: the ring stiffness is a circulant operator, `solve_circulant`
/// diagonalizes it by the DFT (the isotypic projection for the cyclic group
/// Cₙ), turning one n×n solve into n scalar solves. A point load is applied at
/// sector 0; the response spreads around the ring. The per-harmonic (irrep)
/// content is the DFT of the response, built from the exact `CyclicGroup`
/// character table (roots of unity).
///
/// `n` (sectors, clamped 4..=64), `stiffness` (ring coupling κ added to the
/// diagonal to keep the operator nonsingular, 0.05..=4.0).
///
/// Output layout (length `1 + 4·n`):
/// - `[0]` — `n`.
/// - `[1 .. 1+n]`     — the circulant operator's first row (the ring stencil).
/// - `[1+n .. 1+2n]`  — the right-hand side (point load at sector 0).
/// - `[1+2n .. 1+3n]` — the solution `x` around the ring (displacement per
///   sector). All `NaN` if the operator is reported singular.
/// - `[1+3n .. 1+4n]` — per-harmonic magnitude `|X̂_k|` for irreps `k = 0..n−1`
///   (the irrep/mode content of the response).
pub fn cyclic_symmetry(n_in: usize, stiffness: f64) -> Vec<f64> {
    let n = n_in.clamp(4, 64);
    let kappa = stiffness.clamp(0.05, 4.0);
    // Ring stencil: (2+κ) on the diagonal, −1 to each neighbor (wraps).
    let mut first_row = vec![0.0f64; n];
    first_row[0] = 2.0 + kappa;
    first_row[1] = -1.0;
    first_row[n - 1] = -1.0;
    // Point load at sector 0.
    let mut rhs = vec![0.0f64; n];
    rhs[0] = 1.0;

    let sol = solve_circulant(&first_row, &rhs).unwrap_or_else(|_| vec![f64::NAN; n]);

    // Per-harmonic content via the exact character table.
    let grp = CyclicGroup::new(n);
    let harmonics: Vec<f64> = (0..n)
        .map(|k| {
            let mut re = 0.0f64;
            let mut im = 0.0f64;
            for (j, &xj) in sol.iter().enumerate() {
                let (c, s) = grp.character(k, j); // e^{+2πi kj/n}
                re += xj * c; // projection onto conjugate: real part
                im -= xj * s;
            }
            re.hypot(im) / n as f64
        })
        .collect();

    let mut out = Vec::with_capacity(1 + 4 * n);
    out.push(n as f64);
    out.extend_from_slice(&first_row);
    out.extend_from_slice(&rhs);
    out.extend_from_slice(&sol);
    out.extend_from_slice(&harmonics);
    out
}

/* ======================================================================= */
/*  7 · Krylov solver race — CG vs MINRES vs GMRES (fs-solver)              */
/* ======================================================================= */

/// **Krylov convergence race** on a real SPD 2D Laplacian (the assembled
/// 5-point stencil), solved three ways through `fs-solver`'s resumable
/// matrix-free stack: preconditioned CG (`CgState`, identity preconditioner),
/// MINRES (`MinresState`), and restarted GMRES(30) (`GmresState`). Each solver
/// returns its full relative-residual history so the viz can race the curves.
///
/// `n` (interior grid side, clamped 4..=24 → system size m = n²), `maxit`
/// (iteration/inner cap, 20..=1500).
///
/// Output layout (a flat array the viz slices by length prefixes):
/// - `[0]` — `m`, the system size (n²).
/// - `[1]` — `Lcg`; then `Lcg` values — CG relative residual `‖r‖/‖b‖` per
///   iteration.
/// - next `[.]` — `Lminres`; then `Lminres` values — MINRES per iteration.
/// - next `[.]` — `Lgmres`; then `Lgmres` values — GMRES relative residual per
///   completed restart cycle (coarser-grained than the other two).
pub fn krylov_convergence(n_in: usize, maxit_in: usize) -> Vec<f64> {
    let n = n_in.clamp(4, 24);
    let maxit = maxit_in.clamp(20, 1500);
    let m = n * n;
    let a = crate::laplacian_5pt(n);
    let op = CsrOp::symmetric(a);
    let precond = IdentityPrecond;
    let b = vec![1.0f64; m];
    let tol = 1e-10;

    let mut out = Vec::new();
    out.push(m as f64);

    // CG.
    let mut cg = CgState::new(&op, &precond, &b);
    let _ = cg.run(&op, &precond, tol, maxit);
    out.push(cg.history.len() as f64);
    out.extend_from_slice(&cg.history);

    // MINRES.
    let mut mr = MinresState::new(&op, &b);
    let _ = mr.run(&op, tol, maxit);
    out.push(mr.history.len() as f64);
    out.extend_from_slice(&mr.history);

    // GMRES(30).
    let restart = 30usize;
    let mut gm = GmresState::new(&b, restart);
    let cycles = (maxit / restart).max(1);
    let _ = gm.run(&op, &b, tol, cycles, false);
    out.push(gm.history.len() as f64);
    out.extend_from_slice(&gm.history);

    out
}

/* ======================================================================= */
/*  8 · CutFEM quadtree — adaptive mesh + cut quadrature (fs-cutfem)         */
/* ======================================================================= */

/// **CutFEM adaptive quadtree + cut-cell quadrature** for a circular SDF
/// boundary via `fs-cutfem`: a 2:1-balanced `Quadtree` is refined toward the
/// interface (`refine_toward_interface`), each leaf is classified by the
/// interval-certified SDF enclosure (Inside / Outside / Cut — no misclassified
/// cell), and every cut cell gets its real `cut_cell_rules` bulk + interface
/// quadrature (marching-squares polygonization + bisected crossings). The
/// domain is the unit square; the boundary is a circle at (0.5, 0.5).
///
/// `base` (starting uniform level, 2..=4), `target` (max refinement level,
/// base..=7), `radius` (circle radius, 0.15..=0.45).
///
/// Output layout (a flat array the viz slices):
/// - `[0]` — `L`, the leaf count.
/// - `[1 .. 1+4L]` — per leaf: `[cx, cy, size, class]` — cell center, side
///   length, and class (0 = outside Ω, 1 = inside Ω, 2 = cut).
/// - `[1+4L]` — `Qi`, the interface quadrature-point count; then `Qi` blocks of
///   4: `[x, y, nx, ny]` — point position and OUTWARD unit normal.
/// - next `[.]` — `Qb`, the bulk (inside-region) quadrature-point count; then
///   `Qb` blocks of 3: `[x, y, weight]` — point position and integration weight
///   (weights sum to the inside area).
pub fn cutfem_quadtree(base_in: usize, target_in: usize, radius: f64) -> Vec<f64> {
    let base = base_in.clamp(2, 4) as u32;
    let target = (target_in as u32).clamp(base, 7);
    let r = radius.clamp(0.15, 0.45);
    let circle = Circle {
        center: [0.5, 0.5],
        radius: r,
    };
    let mut qt = Quadtree::with_room(base, target);
    qt.refine_toward_interface(&circle, target);

    let depth = 3u32; // cut-cell subdivision depth for quadrature accuracy
    let mut leaves: Vec<f64> = Vec::new();
    let mut iface: Vec<f64> = Vec::new();
    let mut bulk: Vec<f64> = Vec::new();
    let mut leaf_count = 0usize;
    let mut iface_count = 0usize;
    let mut bulk_count = 0usize;

    for c in qt.leaves() {
        let (lo, hi) = qt.rect(c);
        let size = qt.cell_h(c);
        let iv = circle.enclose(lo, hi);
        let class = if iv.hi() < 0.0 {
            1.0 // inside Ω
        } else if iv.lo() > 0.0 {
            0.0 // outside Ω
        } else {
            2.0 // cut
        };
        leaves.push(f64::midpoint(lo[0], hi[0]));
        leaves.push(f64::midpoint(lo[1], hi[1]));
        leaves.push(size);
        leaves.push(class);
        leaf_count += 1;

        if class == 2.0 {
            let rules = cut_cell_rules(&circle, lo, hi, depth);
            for (p, _w, nrm) in &rules.iface {
                iface.push(p[0]);
                iface.push(p[1]);
                iface.push(nrm[0]);
                iface.push(nrm[1]);
                iface_count += 1;
            }
            for (p, w) in &rules.bulk {
                bulk.push(p[0]);
                bulk.push(p[1]);
                bulk.push(*w);
                bulk_count += 1;
            }
        }
    }

    let mut out = Vec::with_capacity(1 + leaves.len() + 1 + iface.len() + 1 + bulk.len());
    out.push(leaf_count as f64);
    out.extend_from_slice(&leaves);
    out.push(iface_count as f64);
    out.extend_from_slice(&iface);
    out.push(bulk_count as f64);
    out.extend_from_slice(&bulk);
    out
}

/* ======================================================================= */
/*  9 · Free-form deformation (fs-xform)                                    */
/* ======================================================================= */

/// **Bernstein free-form deformation** of a point cloud via `fs-xform`: a grid
/// of shape sample points in the z=0.5 mid-plane is warped through a trivariate
/// Bernstein FFD control lattice (`FfdLattice::apply`, exact `bernstein` /
/// `bernstein_derivative` basis), with control-point offsets from a chosen mode.
/// The spatial Jacobian determinant (`spatial_jacobian` + `det3`) is checked at
/// every sample for fold-over (det ≤ 0 — the warp becoming non-invertible).
///
/// `grid` (sample points per side, clamped 4..=20 → P = grid²), `controls`
/// (control lattice count per in-plane axis, 2..=5), `amp` (offset amplitude,
/// clamped ±1.5), `mode`: 0 = shear, 1 = bulge, 2 = twist, 3 = pinch (large amp
/// folds space).
///
/// Output layout (length `1 + 4P + 1 + 4C + 2`):
/// - `[0]` — `P`, the sample-point count.
/// - `[1 .. 1+4P]` — per sample: `[ox, oy, dx, dy]` — original then deformed
///   in-plane position.
/// - `[1+4P]` — `C`, the control-lattice node count (controls²); then `C` blocks
///   of 4: `[ox, oy, dx, dy]` — original then displaced control node (the cage).
/// - next value — fold-over flag (1.0 if any sampled det(∂T/∂x) ≤ 0, else 0.0).
/// - last value — the minimum sampled Jacobian determinant (≤0 ⇒ folded).
pub fn ffd_deform(grid_in: usize, controls_in: usize, amp: f64, mode: u32) -> Vec<f64> {
    let grid = grid_in.clamp(4, 20);
    let nc = controls_in.clamp(2, 5);
    let a = amp.clamp(-1.5, 1.5);
    let ffd = FfdLattice {
        origin: Point3::new(0.0, 0.0, 0.0),
        size: Vec3::new(1.0, 1.0, 1.0),
        counts: [nc, nc, 2],
    };
    let dof = ffd.dof();
    let mut theta = vec![0.0f64; dof];

    // Control-point offsets per mode (in-plane only).
    let offset = |u: f64, v: f64| -> (f64, f64) {
        match mode {
            0 => (a * v, 0.0),                   // shear: top slides in x
            1 => (a * (u - 0.5), a * (v - 0.5)), // bulge: push out from center
            2 => {
                // twist: rotate control node about the center by angle ∝ v
                let ang = a * (v - 0.5) * 2.0;
                let (cx, cy) = (u - 0.5, v - 0.5);
                let nx = cx * ang.cos() - cy * ang.sin();
                let ny = cx * ang.sin() + cy * ang.cos();
                (nx - cx, ny - cy)
            }
            _ => (-a * (u - 0.5) * 2.0, 0.0), // pinch: strong inward x → can fold
        }
    };
    for i in 0..nc {
        for j in 0..nc {
            for k in 0..2 {
                let node = (i * nc + j) * 2 + k;
                let u = i as f64 / (nc as f64 - 1.0);
                let v = j as f64 / (nc as f64 - 1.0);
                let (dx, dy) = offset(u, v);
                theta[3 * node] = dx;
                theta[3 * node + 1] = dy;
                theta[3 * node + 2] = 0.0;
            }
        }
    }

    let p = grid * grid;
    let mut out = Vec::with_capacity(1 + 4 * p + 1 + 4 * nc * nc + 2);
    out.push(p as f64);
    let mut foldover = false;
    let mut min_det = f64::INFINITY;
    for iy in 0..grid {
        for ix in 0..grid {
            let ox = 0.05 + 0.9 * ix as f64 / (grid as f64 - 1.0);
            let oy = 0.05 + 0.9 * iy as f64 / (grid as f64 - 1.0);
            let orig = Point3::new(ox, oy, 0.5);
            let (dx, dy) = match ffd.apply(&theta, orig) {
                Ok(q) => (q.x, q.y),
                Err(_) => (ox, oy),
            };
            let det = match ffd.spatial_jacobian(&theta, orig) {
                Ok(j) => det3(&j),
                Err(_) => 1.0,
            };
            if det <= 0.0 {
                foldover = true;
            }
            min_det = min_det.min(det);
            out.push(ox);
            out.push(oy);
            out.push(dx);
            out.push(dy);
        }
    }
    // Cage nodes (k=0 layer).
    out.push((nc * nc) as f64);
    for i in 0..nc {
        for j in 0..nc {
            let u = i as f64 / (nc as f64 - 1.0);
            let v = j as f64 / (nc as f64 - 1.0);
            let ox = u; // origin (0,0), size (1,1)
            let oy = v;
            let node = (i * nc + j) * 2;
            let dx = ox + theta[3 * node];
            let dy = oy + theta[3 * node + 1];
            out.push(ox);
            out.push(oy);
            out.push(dx);
            out.push(dy);
        }
    }
    out.push(if foldover { 1.0 } else { 0.0 });
    out.push(if min_det.is_finite() {
        min_det
    } else {
        f64::NAN
    });
    out
}

/* ======================================================================= */
/*  The JavaScript boundary (wasm32 only)                                   */
/* ======================================================================= */

#[cfg(target_arch = "wasm32")]
mod wasm {
    use wasm_bindgen::prelude::*;

    #[wasm_bindgen]
    pub fn hodge_decomposition(shape: u32) -> Vec<f64> {
        super::hodge_decomposition(shape)
    }

    #[wasm_bindgen]
    pub fn betti_shapes(shape: u32) -> Vec<f64> {
        super::betti_shapes(shape)
    }

    #[wasm_bindgen]
    pub fn navier_stokes_cavity(
        cells: usize,
        frames: usize,
        re: f64,
        steps_per_frame: usize,
    ) -> Vec<f64> {
        super::navier_stokes_cavity(cells, frames, re, steps_per_frame)
    }

    #[wasm_bindgen]
    pub fn gp_regression(n_train: usize, samples: usize) -> Vec<f64> {
        super::gp_regression(n_train, samples)
    }

    #[wasm_bindgen]
    pub fn cmaes_trace(seed: u32, gens: usize) -> Vec<f64> {
        super::cmaes_trace(seed, gens)
    }

    #[wasm_bindgen]
    pub fn optimal_transport(n: usize, epsilon: f64) -> Vec<f64> {
        super::optimal_transport(n, epsilon)
    }

    #[wasm_bindgen]
    pub fn cyclic_symmetry(n: usize, stiffness: f64) -> Vec<f64> {
        super::cyclic_symmetry(n, stiffness)
    }

    #[wasm_bindgen]
    pub fn krylov_convergence(n: usize, maxit: usize) -> Vec<f64> {
        super::krylov_convergence(n, maxit)
    }

    #[wasm_bindgen]
    pub fn cutfem_quadtree(base: usize, target: usize, radius: f64) -> Vec<f64> {
        super::cutfem_quadtree(base, target, radius)
    }

    #[wasm_bindgen]
    pub fn ffd_deform(grid: usize, controls: usize, amp: f64, mode: u32) -> Vec<f64> {
        super::ffd_deform(grid, controls, amp, mode)
    }
}
