//! Stability analysis (plan §8.2, bead tfz.15): geometric-stiffness
//! buckling pencils and eigenvalue derivatives — because the
//! steel-minimizing frame that ignores buckling IS A PAPERCLIP.
//!
//! The pencil `(K + λ K_G(σ₀)) φ = 0` is solved MATRIX-FREE with
//! fs-la's LOBPCG (the named consumer) on the Cholesky-reduced
//! symmetric operator Ã = L⁻¹(−K_G)L⁻ᵀ (K = LLᵀ): the largest Ritz
//! value θ gives the critical load λ = 1/θ. The dense factorization is
//! fixture-gated (the cond-probe discipline); production-scale shifted
//! iterations are a recorded no-claim. Independently scaled matrices preserve
//! physical load factors; transformed and original-pencil residuals gate every
//! result. Small pencils that cannot fit a LOBPCG block use fs-la's native dense
//! symmetric eigensolver within an explicit 64-DOF bound.
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
use crate::linear::{Formulation, LinearProblem};
use crate::mesh2::{quad_points, shapes_at};
use fs_la::eigen::{EigenPair, LobpcgState, householder_eigh, lobpcg_run};
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
    /// Recomputed relative original-pencil residuals, each <= 1e-8.
    pub relative_residuals: Vec<f64>,
    /// Free-DOF map: `dof_map[full_dof] = Some(free_index)`.
    pub dof_map: Vec<Option<usize>>,
    /// LOBPCG iterations spent.
    pub iters: usize,
}

/// The reduced pencil: (stiffness, geometric stiffness, DOF map,
/// prebuckling displacement).
pub type ReducedPencil = (Csr, Csr, Vec<Option<usize>>, Vec<[f64; 2]>);

/// The reduced (free-DOF) stiffness, geometric stiffness, and DOF map
/// for a linear problem at its reference load.
///
/// # Errors
/// Propagates the prebuckling solve's [`SolidError`].
pub fn reduced_pencil(problem: &LinearProblem<'_>) -> Result<ReducedPencil, SolidError> {
    let u0 = problem.solve()?;
    let (lambda, mu) = crate::linear::lame(problem.youngs, problem.poisson, problem.plane);
    let mesh = problem.mesh;
    let n = mesh.node_count();
    // Buckling perturbations obey the same homogeneous essential constraints
    // as the reference solve, including component-wise symmetry supports.
    let mut is_fixed = vec![false; 2 * n];
    for (patch, _) in &problem.dirichlet {
        for node in mesh.patch_nodes(*patch) {
            is_fixed[2 * node] = true;
            is_fixed[2 * node + 1] = true;
        }
    }
    for (patch, component) in &problem.symmetry {
        for node in mesh.patch_nodes(*patch) {
            is_fixed[2 * node + component] = true;
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
        let pts = quad_points(nn);
        // The B-bar equilibrium solve uses an element-average dilatation.
        // Recover that same stress here: restoring the pointwise dilatation
        // would multiply its unresolved oscillations by the bulk modulus and
        // corrupt the geometric stiffness near incompressibility.
        let mut mean_dilatation = 0.0;
        if problem.formulation == Formulation::BBar {
            let mut volume = 0.0;
            for &(xi, eta, w) in &pts {
                let (_, grads, det) = shapes_at(&mesh.nodes, conn, xi, eta);
                let eps = reference_strain(conn, &grads, &u0);
                mean_dilatation += w * det * (eps[0] + eps[1]);
                volume += w * det;
            }
            mean_dilatation /= volume;
        }
        // Geometric stiffness from the prebuckling stress state σ₀(u₀):
        // K_G[ai][bj] = δ_ij ∫ ∇N_aᵀ σ₀ ∇N_b.
        let mut ge = vec![vec![0.0f64; nn]; nn];
        for &(xi, eta, w) in &pts {
            let (_, grads, det) = shapes_at(&mesh.nodes, conn, xi, eta);
            let wq = w * det;
            let eps = reference_strain(conn, &grads, &u0);
            let dilatation = match problem.formulation {
                Formulation::Standard => eps[0] + eps[1],
                Formulation::BBar => mean_dilatation,
            };
            // This split also avoids cancelling large Lamé terms when the
            // resolved deformation is nearly isochoric.
            let mean_stress = (lambda + mu) * dilatation;
            let deviator = mu * (eps[0] - eps[1]);
            let sxx = mean_stress + deviator;
            let syy = mean_stress - deviator;
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

/// Small strain [εxx, εyy, γxy], removing rigid displacement before contraction.
fn reference_strain(conn: &[usize], grads: &[[f64; 2]], u: &[[f64; 2]]) -> [f64; 3] {
    let origin = u[conn[0]];
    let mut eps = [0.0; 3];
    for (a, &node) in conn.iter().enumerate().skip(1) {
        let du = [u[node][0] - origin[0], u[node][1] - origin[1]];
        eps[0] += grads[a][0] * du[0];
        eps[1] += grads[a][1] * du[1];
        eps[2] += grads[a][1] * du[0] + grads[a][0] * du[1];
    }
    eps
}

/// Solve `(K + lambda K_G) phi = 0` for the lowest positive critical loads.
///
/// Independent matrix scaling removes physical units before Cholesky/LOBPCG.
/// Both transformed and original-pencil relative residuals must be <= 1e-8;
/// exhausted iterations refuse instead of publishing an unconverged prediction.
/// No absolute cutoff on the inverse load can erase a small reference load.
///
/// If the requested block cannot fit LOBPCG's search space, pencils with at
/// most 64 free DOFs use the native dense symmetric eigensolver. The Cholesky
/// reduction is bounded to 4096 free DOFs. This is not a scalable sparse
/// buckling algorithm, a nonlinear stability limit or a convergence certificate.
///
/// # Errors
/// Invalid, nonfinite, nonsymmetric or oversized pencils and invalid DOF maps
/// return `InvalidInput`. Non-SPD stiffness or unconverged eigenpairs return
/// `SolveFailed`. Near-zero unresolved eigenvalues and unrepresentable loads
/// or work-normalized modes refuse rather than being silently discarded.
#[allow(clippy::too_many_lines)]
pub fn buckling_loads(
    k: &Csr,
    kg: &Csr,
    dof_map: &[Option<usize>],
    nev: usize,
    steps: usize,
) -> Result<BucklingResult, SolidError> {
    let invalid = |what: &str| SolidError::InvalidInput { what: what.into() };
    let n = k.nrows();
    if n == 0 || n > 4096 || k.ncols() != n || kg.nrows() != n || kg.ncols() != n {
        return Err(invalid(
            "buckling requires matching nonempty square pencils with at most 4096 DOFs",
        ));
    }
    if nev == 0 || steps == 0 || !dof_map.len().is_multiple_of(2) {
        return Err(invalid(
            "buckling requires positive mode/iteration allowances and a planar DOF map",
        ));
    }
    let mut seen = vec![false; n];
    for &free in dof_map.iter().flatten() {
        if free >= n || seen[free] {
            return Err(invalid(
                "buckling DOF map must cover each free DOF exactly once",
            ));
        }
        seen[free] = true;
    }
    if seen.contains(&false) {
        return Err(invalid("buckling DOF map is missing a free DOF"));
    }
    let count = nev.min(n);
    let mut iterative = LobpcgState::try_new(n, count);
    if iterative.is_none() && n > 64 {
        return Err(invalid(
            "requested buckling block exceeds the bounded eigensolver workspace",
        ));
    }
    let (k_scale, scaled_k) = scaled_symmetric(k)?;
    let (g_scale, scaled_g) = scaled_symmetric(kg)?;
    if k_scale == 0.0 {
        return Err(SolidError::SolveFailed {
            iters: 0,
            rel_residual: f64::INFINITY,
        });
    }
    let chol = cholesky(&scaled_k.to_dense(), n).map_err(|_| SolidError::SolveFailed {
        iters: 0,
        rel_residual: f64::INFINITY,
    })?;
    if g_scale == 0.0 {
        return Ok(BucklingResult {
            loads: Vec::new(),
            modes: Vec::new(),
            relative_residuals: Vec::new(),
            dof_map: dof_map.to_vec(),
            iters: 0,
        });
    }
    let fwd = |b: &mut [f64]| {
        for i in 0..n {
            let mut v = b[i];
            for (j, bv) in b.iter().enumerate().take(i) {
                v -= chol.l(i, j) * bv;
            }
            b[i] = v / chol.l(i, i);
        }
    };
    let bwd = |b: &mut [f64]| {
        for i in (0..n).rev() {
            let mut v = b[i];
            for (j, bv) in b.iter().enumerate().skip(i + 1) {
                v -= chol.l(j, i) * bv;
            }
            b[i] = v / chol.l(i, i);
        }
    };
    // Ahat = Lhat^-1 (-Ghat) Lhat^-T. Reuse the operator's two work vectors.
    let scratch = std::cell::RefCell::new((vec![0.0; n], vec![0.0; n]));
    let op = |x: &[f64], out: &mut [f64]| {
        let mut scratch = scratch.borrow_mut();
        let (t, s) = &mut *scratch;
        t.copy_from_slice(x);
        bwd(t);
        scaled_g.spmv(t, s);
        for value in s.iter_mut() {
            *value = -*value;
        }
        fwd(s);
        out.copy_from_slice(s);
    };
    let (pairs, iters) = if let Some(state) = &mut iterative {
        let pairs = lobpcg_run(&op, state, steps, true, &|r, out| out.copy_from_slice(r));
        (pairs, state.iters)
    } else {
        (small_buckling_pairs(n, count, &op)?, 0)
    };
    let mut loads = Vec::new();
    let mut modes = Vec::new();
    let mut relative_residuals = Vec::new();
    let mut k_phi = vec![0.0; n];
    let mut g_phi = vec![0.0; n];
    for pair in &pairs {
        let transformed = if pair.value == 0.0 && pair.residual == 0.0 {
            0.0
        } else {
            pair.residual / pair.value.abs()
        };
        if !pair.value.is_finite()
            || !transformed.is_finite()
            || transformed > 1e-8
            || pair.vector.iter().any(|v| !v.is_finite())
        {
            return Err(SolidError::SolveFailed {
                iters,
                rel_residual: transformed,
            });
        }
        if pair.value <= 0.0 {
            continue;
        }
        let mut load = (k_scale / g_scale) / pair.value;
        if !load.is_finite() || load <= 0.0 {
            load = (k_scale / pair.value) / g_scale;
        }
        if !load.is_finite() || load <= 0.0 {
            return Err(invalid(
                "critical buckling load factor is not representable",
            ));
        }
        let mut phi = pair.vector.clone();
        bwd(&mut phi);
        // Remove arbitrary amplitude before evaluating work and residuals.
        let amplitude = phi.iter().fold(0.0_f64, |scale, v| scale.max(v.abs()));
        if !amplitude.is_finite() || amplitude <= 0.0 {
            return Err(invalid("buckling mode is not representable"));
        }
        for value in &mut phi {
            *value /= amplitude;
        }
        scaled_k.spmv(&phi, &mut k_phi);
        scaled_g.spmv(&phi, &mut g_phi);
        // ||theta Khat phi + Ghat phi|| / (|theta| ||Khat phi|| + ||Ghat phi||),
        // dividing out the larger scalar before arithmetic, without dimensional floors.
        let (ks, gs) = if pair.value <= 1.0 {
            (pair.value, 1.0)
        } else {
            (1.0, 1.0 / pair.value)
        };
        let (mut residual, mut knorm, mut gnorm) = (0.0_f64, 0.0_f64, 0.0_f64);
        for i in 0..n {
            residual = residual.hypot(ks.mul_add(k_phi[i], gs * g_phi[i]));
            knorm = knorm.hypot(k_phi[i]);
            gnorm = gnorm.hypot(g_phi[i]);
        }
        let relative = residual / (ks * knorm + gs * gnorm);
        if !relative.is_finite() || relative > 1e-8 {
            return Err(SolidError::SolveFailed {
                iters,
                rel_residual: relative,
            });
        }
        let work = -phi.iter().zip(&g_phi).map(|(a, b)| a * b).sum::<f64>();
        // Avoid forming work*g_scale: that can overflow even for finite modes.
        let normalization = work.sqrt() * g_scale.sqrt();
        if !normalization.is_finite() || normalization <= 0.0 {
            return Err(invalid(
                "buckling mode has no representable positive geometric work",
            ));
        }
        for value in &mut phi {
            *value /= normalization;
        }
        if phi.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "work-normalized buckling mode is not representable",
            ));
        }
        loads.push(load);
        modes.push(phi);
        relative_residuals.push(relative);
    }
    let mut order: Vec<usize> = (0..loads.len()).collect();
    order.sort_by(|&a, &b| loads[a].total_cmp(&loads[b]));
    Ok(BucklingResult {
        loads: order.iter().map(|&i| loads[i]).collect(),
        modes: order.iter().map(|&i| modes[i].clone()).collect(),
        relative_residuals: order.iter().map(|&i| relative_residuals[i]).collect(),
        dof_map: dof_map.to_vec(),
        iters,
    })
}

/// Remove physical units without changing sparsity or hiding asymmetry.
fn scaled_symmetric(matrix: &Csr) -> Result<(f64, Csr), SolidError> {
    let mut scale = 0.0_f64;
    for row in 0..matrix.nrows() {
        for &value in matrix.row(row).1 {
            if !value.is_finite() {
                return Err(SolidError::InvalidInput {
                    what: "buckling matrices must be finite".into(),
                });
            }
            scale = scale.max(value.abs());
        }
    }
    let mut scaled = matrix.clone();
    if scale > 0.0 {
        for row in 0..matrix.nrows() {
            let (columns, values) = matrix.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                if (value / scale - matrix.get(column, row) / scale).abs() > 64.0 * f64::EPSILON {
                    return Err(SolidError::InvalidInput {
                        what: "buckling matrices must be symmetric".into(),
                    });
                }
            }
            for value in scaled.row_values_mut(row) {
                *value /= scale;
            }
        }
    }
    Ok((scale, scaled))
}

/// Exact native dense reduction for small pencils that cannot fit a 3b search space.
fn small_buckling_pairs(
    n: usize,
    count: usize,
    op: &impl Fn(&[f64], &mut [f64]),
) -> Result<Vec<EigenPair>, SolidError> {
    let mut matrix = vec![0.0; n * n];
    let mut unit = vec![0.0; n];
    let mut column = vec![0.0; n];
    for j in 0..n {
        unit[j] = 1.0;
        op(&unit, &mut column);
        unit[j] = 0.0;
        for i in 0..n {
            matrix[i * n + j] = column[i];
        }
    }
    if matrix.iter().any(|v| !v.is_finite()) {
        return Err(SolidError::InvalidInput {
            what: "buckling reduction is not representable".into(),
        });
    }
    for i in 0..n {
        for j in 0..i {
            let mean = f64::midpoint(matrix[i * n + j], matrix[j * n + i]);
            matrix[i * n + j] = mean;
            matrix[j * n + i] = mean;
        }
    }
    let (values, vectors) = householder_eigh(&matrix, n);
    Ok((0..n)
        .rev()
        .take(count)
        .map(|j| {
            let vector: Vec<_> = (0..n).map(|i| vectors[i * n + j]).collect();
            op(&vector, &mut column);
            let residual = (0..n).fold(0.0_f64, |norm, i| {
                norm.hypot((-values[j]).mul_add(vector[i], column[i]))
            });
            EigenPair {
                value: values[j],
                vector,
                residual,
            }
        })
        .collect())
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
    weights.iter().zip(dloads).map(|(w, d)| w / total * d).sum()
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
