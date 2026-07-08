//! Stability analysis (plan §8.2, bead tfz.15): geometric-stiffness
//! buckling pencils and eigenvalue derivatives — because the
//! steel-minimizing frame that ignores buckling IS A PAPERCLIP.
//!
//! The pencil `(K + λ K_G(σ₀)) φ = 0` is solved MATRIX-FREE with
//! fs-la's LOBPCG (the named consumer) on the Cholesky-reduced
//! symmetric operator Ã = L⁻¹(−K_G)L⁻ᵀ (K = LLᵀ): the largest Ritz
//! value θ gives the critical load λ = 1/θ. The dense factorization is
//! fixture-gated (the cond-probe discipline); production-scale shifted
//! iterations are a recorded no-claim.
//!
//! Eigenvalue derivatives ship with the pencil: the DIRECT derivative
//! `dλ/ds = φᵀ(∂K/∂s)φ / (φᵀ(−K_G)φ)` at frozen prebuckling stress,
//! plus the documented clustered-eigenvalue handling — individual
//! eigenvalues are NOT differentiable where branches cross, so
//! optimization consumes the smooth Kreisselmeier–Steinhauser
//! aggregate `λ_ks = −ln(Σ e^{−ρλᵢ})/ρ` whose derivative is the
//! softmax-weighted mode sum (the classic trap, measured in the
//! battery). The prebuckling-adjoint chain (σ₀ depending on s) is the
//! fs-ad/ASCENT integration successor.

use crate::SolidError;
use crate::linear::LinearProblem;
use crate::mesh2::{quad_points, shapes_at};
use fs_la::eigen::{LobpcgState, lobpcg_run};
use fs_la::factor::cholesky;
use fs_sparse::{Coo, Csr};

/// Buckling analysis output.
#[derive(Debug, Clone)]
pub struct BucklingResult {
    /// Critical load factors, ascending (multipliers on the reference
    /// load producing σ₀).
    pub loads: Vec<f64>,
    /// Modes (free-DOF vectors, `−K_G`-normalized), one per load.
    pub modes: Vec<Vec<f64>>,
    /// Free-DOF map: `dof_map[full_dof] = Some(free_index)`.
    pub dof_map: Vec<Option<usize>>,
    /// LOBPCG iterations spent.
    pub iters: usize,
}

/// The reduced (free-DOF) stiffness, geometric stiffness, and DOF map
/// for a linear problem at its reference load.
///
/// # Errors
/// Propagates the prebuckling solve's [`SolidError`].
pub fn reduced_pencil(
    problem: &LinearProblem<'_>,
) -> Result<(Csr, Csr, Vec<Option<usize>>, Vec<[f64; 2]>), SolidError> {
    let u0 = problem.solve()?;
    let (lambda, mu) = crate::linear::lame(problem.youngs, problem.poisson, problem.plane);
    let mesh = problem.mesh;
    let n = mesh.node_count();
    // Fixed DOFs: every component on a Dirichlet patch.
    let mut is_fixed = vec![false; 2 * n];
    for (patch, _) in &problem.dirichlet {
        for node in mesh.patch_nodes(*patch) {
            is_fixed[2 * node] = true;
            is_fixed[2 * node + 1] = true;
        }
    }
    let mut dof_map: Vec<Option<usize>> = Vec::with_capacity(2 * n);
    let mut nf = 0usize;
    for &f in &is_fixed {
        if f {
            dof_map.push(None);
        } else {
            dof_map.push(Some(nf));
            nf += 1;
        }
    }
    let mut kc = Coo::new(nf, nf);
    let mut gc = Coo::new(nf, nf);
    for conn in &mesh.elems {
        let nn = conn.len();
        let (ke, _) = problem.element(conn, lambda, mu);
        // Geometric stiffness from the prebuckling stress state σ₀(u₀):
        // K_G[ai][bj] = δ_ij ∫ ∇N_aᵀ σ₀ ∇N_b.
        let mut ge = vec![vec![0.0f64; nn]; nn];
        for &(xi, eta, w) in &quad_points(nn) {
            let (_, grads, det) = shapes_at(&mesh.nodes, conn, xi, eta);
            let wq = w * det;
            // Small-strain stress at this point from u₀.
            let mut eps = [0.0f64; 3]; // [εxx, εyy, γxy]
            for (a, &node) in conn.iter().enumerate() {
                eps[0] += grads[a][0] * u0[node][0];
                eps[1] += grads[a][1] * u0[node][1];
                eps[2] += grads[a][1] * u0[node][0] + grads[a][0] * u0[node][1];
            }
            let sxx = (lambda + 2.0 * mu) * eps[0] + lambda * eps[1];
            let syy = lambda * eps[0] + (lambda + 2.0 * mu) * eps[1];
            let sxy = mu * eps[2];
            for a in 0..nn {
                let ga = grads[a];
                let sg = [sxx * ga[0] + sxy * ga[1], sxy * ga[0] + syy * ga[1]];
                for b in 0..nn {
                    ge[a][b] += wq * (sg[0] * grads[b][0] + sg[1] * grads[b][1]);
                }
            }
        }
        for a in 0..nn {
            for ca in 0..2 {
                let Some(ia) = dof_map[2 * conn[a] + ca] else {
                    continue;
                };
                for b in 0..nn {
                    for cb in 0..2 {
                        let Some(ib) = dof_map[2 * conn[b] + cb] else {
                            continue;
                        };
                        let kv = ke[2 * a + ca][2 * b + cb];
                        if kv != 0.0 {
                            kc.push(ia, ib, kv);
                        }
                        if ca == cb {
                            let gv = ge[a][b];
                            if gv != 0.0 {
                                gc.push(ia, ib, gv);
                            }
                        }
                    }
                }
            }
        }
    }
    Ok((kc.assemble(), gc.assemble(), dof_map, u0))
}

/// Solve the buckling pencil `(K + λ K_G) φ = 0` for the `nev` lowest
/// positive critical loads, matrix-free LOBPCG on the reduced
/// operator.
///
/// # Errors
/// [`SolidError::SolveFailed`] if the stiffness is not SPD (fixture
/// misconfiguration) — carried as a structured refusal.
///
/// # Panics
/// On fixture-gate violation (free DOFs > 4096; the dense Cholesky
/// reduction is for conformance-scale systems).
pub fn buckling_loads(
    k: &Csr,
    kg: &Csr,
    dof_map: &[Option<usize>],
    nev: usize,
    steps: usize,
) -> Result<BucklingResult, SolidError> {
    let n = k.nrows();
    assert!(n <= 4096, "dense pencil reduction is fixture-gated (n = {n})");
    let chol = cholesky(&k.to_dense(), n).map_err(|_| SolidError::SolveFailed {
        iters: 0,
        rel_residual: f64::INFINITY,
    })?;
    // Lower-triangular access for the two half-solves.
    let l_at = |i: usize, j: usize| chol.l(i, j);
    let fwd = |b: &mut [f64]| {
        for i in 0..n {
            let mut v = b[i];
            for (kk, bv) in b.iter().enumerate().take(i) {
                v -= l_at(i, kk) * bv;
            }
            b[i] = v / l_at(i, i);
        }
    };
    let bwd = |b: &mut [f64]| {
        for i in (0..n).rev() {
            let mut v = b[i];
            for (kk, bv) in b.iter().enumerate().skip(i + 1) {
                v -= l_at(kk, i) * bv;
            }
            b[i] = v / l_at(i, i);
        }
    };
    // Ã v = L⁻¹ (−K_G) L⁻ᵀ v, symmetric by construction.
    let op = |v: &[f64], out: &mut [f64]| {
        let mut t = v.to_vec();
        bwd(&mut t);
        let mut s = vec![0.0f64; n];
        kg.spmv(&t, &mut s);
        for x in &mut s {
            *x = -*x;
        }
        fwd(&mut s);
        out.copy_from_slice(&s);
    };
    let ident = |r: &[f64], out: &mut [f64]| out.copy_from_slice(r);
    let mut state = LobpcgState::new(n, nev.min(n));
    let pairs = lobpcg_run(&op, &mut state, steps, true, &ident);
    let mut loads = Vec::new();
    let mut modes = Vec::new();
    for p in &pairs {
        if p.value <= 1e-12 {
            continue; // not a compressive (positive-load) mode
        }
        let lam = 1.0 / p.value;
        // φ = L⁻ᵀ x, normalized so φᵀ(−K_G)φ = 1.
        let mut phi = p.vector.clone();
        bwd(&mut phi);
        let mut kg_phi = vec![0.0f64; n];
        kg.spmv(&phi, &mut kg_phi);
        let m: f64 = -phi.iter().zip(&kg_phi).map(|(a, b)| a * b).sum::<f64>();
        if m > 0.0 {
            let s = 1.0 / m.sqrt();
            for x in &mut phi {
                *x *= s;
            }
        }
        loads.push(lam);
        modes.push(phi);
    }
    let mut order: Vec<usize> = (0..loads.len()).collect();
    order.sort_by(|&a, &b| loads[a].partial_cmp(&loads[b]).expect("finite loads"));
    Ok(BucklingResult {
        loads: order.iter().map(|&i| loads[i]).collect(),
        modes: order.iter().map(|&i| modes[i].clone()).collect(),
        dof_map: dof_map.to_vec(),
        iters: state.iters,
    })
}

/// DIRECT pencil eigenvalue derivative at frozen prebuckling stress:
/// `dλ/ds = φᵀ(∂K/∂s)φ / (φᵀ(−K_G)φ)` with the mode −K_G-normalized
/// (denominator 1). `dk_ds` is the derivative of the reduced stiffness
/// (e.g. the group-restricted stiffness for a per-group scale
/// parameter).
#[must_use]
pub fn eigenvalue_derivative(mode: &[f64], dk_ds: &Csr, kg: &Csr) -> f64 {
    let n = mode.len();
    let mut t = vec![0.0f64; n];
    dk_ds.spmv(mode, &mut t);
    let num: f64 = mode.iter().zip(&t).map(|(a, b)| a * b).sum();
    let mut g = vec![0.0f64; n];
    kg.spmv(mode, &mut g);
    let den: f64 = -mode.iter().zip(&g).map(|(a, b)| a * b).sum::<f64>();
    num / den
}

/// The smooth Kreisselmeier–Steinhauser aggregate of a (possibly
/// clustered) set of critical loads: `λ_ks = −ln(Σ e^{−ρ(λᵢ−λ₁)})/ρ +
/// λ₁` (shifted for overflow safety), a CONSERVATIVE lower envelope
/// (λ_ks ≤ min λᵢ) that stays differentiable where branches cross.
#[must_use]
pub fn ks_aggregate(loads: &[f64], rho: f64) -> f64 {
    let lmin = loads.iter().copied().fold(f64::INFINITY, f64::min);
    let sum: f64 = loads.iter().map(|&l| (-rho * (l - lmin)).exp()).sum();
    lmin - sum.ln() / rho
}

/// Derivative of the KS aggregate from per-eigenvalue derivatives:
/// the softmax-weighted sum — smooth through clusters where the
/// individual `dλᵢ/ds` are not even well-defined.
#[must_use]
pub fn ks_aggregate_derivative(loads: &[f64], dloads: &[f64], rho: f64) -> f64 {
    let lmin = loads.iter().copied().fold(f64::INFINITY, f64::min);
    let weights: Vec<f64> = loads.iter().map(|&l| (-rho * (l - lmin)).exp()).collect();
    let total: f64 = weights.iter().sum();
    weights
        .iter()
        .zip(dloads)
        .map(|(w, d)| w / total * d)
        .sum()
}

/// Restrict the derivative-carrying stiffness to an element group: the
/// reduced `∂K/∂s` for `K(s) = K_rest + s·K_group` (per-group Young's
/// scale — the topo-sizing lever).
#[must_use]
pub fn group_stiffness(
    problem: &LinearProblem<'_>,
    dof_map: &[Option<usize>],
    group: &dyn Fn(usize) -> bool,
) -> Csr {
    let (lambda, mu) = crate::linear::lame(problem.youngs, problem.poisson, problem.plane);
    let nf = dof_map.iter().flatten().count();
    let mut coo = Coo::new(nf, nf);
    for (e, conn) in problem.mesh.elems.iter().enumerate() {
        if !group(e) {
            continue;
        }
        let nn = conn.len();
        let (ke, _) = problem.element(conn, lambda, mu);
        for a in 0..nn {
            for ca in 0..2 {
                let Some(ia) = dof_map[2 * conn[a] + ca] else {
                    continue;
                };
                for b in 0..nn {
                    for cb in 0..2 {
                        let Some(ib) = dof_map[2 * conn[b] + cb] else {
                            continue;
                        };
                        let v = ke[2 * a + ca][2 * b + cb];
                        if v != 0.0 {
                            coo.push(ia, ib, v);
                        }
                    }
                }
            }
        }
    }
    coo.assemble()
}

/// Richardson discretization-error indicator on a critical load from
/// two mesh levels (coarse, fine, refinement ratio 2, order 2):
/// returns (extrapolated λ, |indicator| on the fine value).
#[must_use]
pub fn lambda_indicator(coarse: f64, fine: f64) -> (f64, f64) {
    let extrap = fine + (fine - coarse) / 3.0;
    (extrap, (extrap - fine).abs())
}

/// Evidence row for a buckling analysis (ledger-style JSON).
#[must_use]
pub fn evidence_row(loads: &[f64], indicator: f64) -> String {
    use std::fmt::Write as _;
    let mut s = String::from("{\"kind\":\"buckling\",\"loads\":[");
    for (i, l) in loads.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        let _ = write!(s, "{l:.6e}");
    }
    let _ = write!(s, "],\"lambda_indicator\":{indicator:.3e}}}");
    s
}

/// Map a free-DOF mode back to full nodal vectors (fixed DOFs zero).
#[must_use]
pub fn expand_mode(mode: &[f64], dof_map: &[Option<usize>]) -> Vec<[f64; 2]> {
    let n = dof_map.len() / 2;
    let mut out = vec![[0.0f64; 2]; n];
    for (full, slot) in dof_map.iter().enumerate() {
        if let Some(free) = slot {
            out[full / 2][full % 2] = mode[*free];
        }
    }
    out
}
