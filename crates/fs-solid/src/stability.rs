//! Stability analysis (plan §8.2, bead tfz.15): geometric-stiffness
//! buckling pencils and eigenvalue derivatives — because the
//! steel-minimizing frame that ignores buckling IS A PAPERCLIP.
//!
//! The pencil `(K + λ K_G(σ₀)) φ = 0` uses a symmetric elastic-factor
//! reduction. Small pencils retain the dense Cholesky/LOBPCG route; above
//! 4096 free DOFs, bounded native sparse LDLᵀ whitening delegates to the
//! generic fs-spectral service. Largest algebraic inverse-load roots select
//! the lowest positive critical loads. Independently scaled matrices preserve
//! physical units; transformed and original-pencil residuals qualify the
//! actual rounded load returned to the caller. Small pencils that cannot
//! fit a LOBPCG block use fs-la's native dense symmetric eigensolver within
//! an explicit 64-DOF bound. Sparse fill/work limits remain explicit and
//! numerical convergence does not certify spectrum completeness.
//!
//! Eigenvalue derivatives ship with the pencil: the DIRECT derivative
//! `dλ/ds = φᵀ(∂K/∂s)φ / (φᵀ(−K_G)φ)` at frozen prebuckling stress,
//! plus the documented clustered-eigenvalue handling — individual
//! eigenvalues are NOT differentiable where branches cross, so
//! optimization consumes the smooth Kreisselmeier–Steinhauser
//! aggregate `λ_ks = −ln(Σ e^{−ρλᵢ})/ρ` whose derivative is the
//! softmax-weighted mode sum (the classic trap, measured in the
//! battery). [`sensitivity`] supplies the complete element-modulus derivative:
//! direct elastic and geometric terms plus one reference-equilibrium adjoint
//! per mode, with the original Standard/B-bar strain operator and supports.
//! Individual derivatives refuse unresolved clusters; the complete-packet
//! aggregate uses one equilibrium adjoint and admits repeated interior modes.

use crate::SolidError;
use crate::linear::{Formulation, LinearProblem};
use crate::mesh2::{quad_points, shapes_at};
use fs_la::eigen::{EigenPair, LobpcgState, householder_eigh, lobpcg_run};
use fs_la::factor::{FactorRunError, cholesky_cancellable};
use fs_sparse::{Coo, Csr};

/// Total buckling derivatives through the reference elastic equilibrium.
pub mod sensitivity;
/// Proportional buckling of source-bound three-dimensional prestress assemblies.
pub mod solid3;
/// Buckling with bounded sparse factorization and the generic spectral service.
pub mod sparse;

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
    reduced_pencil_scaled(problem, None, &|| false)
}

/// Reference equilibrium and buckling pencil with one positive Young's-modulus
/// multiplier per original element. The same multipliers enter stiffness,
/// prescribed-displacement lifting and recovered geometric prestress.
/// Dead loads, geometry, Poisson ratio and supports remain fixed.
///
/// # Errors
/// Propagates scale admission and the reference equilibrium's `SolidError`.
pub fn reduced_pencil_with_element_scales(
    problem: &LinearProblem<'_>,
    scales: &[f64],
) -> Result<ReducedPencil, SolidError> {
    reduced_pencil_scaled(problem, Some(scales), &|| false)
}

fn reduced_pencil_scaled(
    problem: &LinearProblem<'_>,
    scales: Option<&[f64]>,
    cancelled: &impl Fn() -> bool,
) -> Result<ReducedPencil, SolidError> {
    let u0 = problem.solve_scaled(scales, cancelled)?;
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
    for (element, conn) in mesh.elems.iter().enumerate() {
        check_cancelled(cancelled)?;
        let scale = scales.map_or(1.0, |values| values[element]);
        let (lambda, mu) = (lambda * scale, mu * scale);
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
    check_cancelled(cancelled)?;
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
/// most 64 free DOFs use the native dense symmetric eigensolver. Up to 4096
/// free DOFs retain the Cholesky route. Larger pencils use bounded sparse
/// positive-factor whitening with default [`sparse::SparseBucklingLimits`].
/// The explicit sparse entry accepts caller-selected resource allowances.
/// Neither route certifies missed-mode exclusion, nonlinear collapse or mesh
/// convergence; exhausted fill and work allowances remain explicit refusals.
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
    buckling_loads_cancellable(k, kg, dof_map, nev, steps, &|| false)
}

/// The same physical pencil with cancellation during factorization and between
/// bounded LOBPCG advances. The map is component-independent: it can represent
/// planar, three-dimensional, or other displacement DOFs without padding.
///
/// The callback returns true to cancel. Cancellation publishes no partial modes.
/// Small pencils retain the dense factor/workspace caps; larger pencils use
/// [`sparse::SparseBucklingLimits`] and the generic spectral service. Both
/// qualify the original physical equation and preserve cancellation semantics.
#[allow(clippy::too_many_lines)]
pub fn buckling_loads_cancellable(
    k: &Csr,
    kg: &Csr,
    dof_map: &[Option<usize>],
    nev: usize,
    steps: usize,
    cancelled: &impl Fn() -> bool,
) -> Result<BucklingResult, SolidError> {
    check_cancelled(cancelled)?;
    let invalid = |what: &str| SolidError::InvalidInput { what: what.into() };
    let n = k.nrows();
    if n > 4096 {
        return sparse::buckling_loads_cancellable(
            k,
            kg,
            dof_map,
            nev,
            steps,
            &sparse::SparseBucklingLimits::default(),
            cancelled,
        );
    }
    if n == 0 || n > 4096 || k.ncols() != n || kg.nrows() != n || kg.ncols() != n {
        return Err(invalid(
            "buckling requires matching nonempty square pencils with at most 4096 DOFs",
        ));
    }
    if nev == 0 || steps == 0 {
        return Err(invalid(
            "buckling requires positive mode and iteration allowances",
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
    let (k_scale, scaled_k) = scaled_symmetric(k, cancelled)?;
    let (g_scale, scaled_g) = scaled_symmetric(kg, cancelled)?;
    if k_scale == 0.0 {
        return Err(SolidError::SolveFailed {
            iters: 0,
            rel_residual: f64::INFINITY,
        });
    }
    check_cancelled(cancelled)?;
    let chol =
        cholesky_cancellable(&scaled_k.to_dense(), n, cancelled).map_err(|error| match error {
            FactorRunError::Cancelled => SolidError::Cancelled,
            FactorRunError::Factor(_) => SolidError::SolveFailed {
                iters: 0,
                rel_residual: f64::INFINITY,
            },
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
        let mut pairs = Vec::new();
        for _ in 0..steps {
            check_cancelled(cancelled)?;
            let before = state.iters;
            pairs = lobpcg_run(&op, state, 1, true, &|r, out| out.copy_from_slice(r));
            check_cancelled(cancelled)?;
            // The native state is resumable. Stop once its genuine residuals
            // resolve the requested block, leaving headroom for the original
            // physical-pencil residual gate after the coordinate transform.
            if state.iters == before || pairs.iter().all(|pair| transformed_residual(pair) <= 1e-10)
            {
                break;
            }
        }
        (pairs, state.iters)
    } else {
        (small_buckling_pairs(n, count, &op, cancelled)?, 0)
    };
    let mut loads = Vec::new();
    let mut modes = Vec::new();
    let mut relative_residuals = Vec::new();
    let mut k_phi = vec![0.0; n];
    let mut g_phi = vec![0.0; n];
    for pair in &pairs {
        check_cancelled(cancelled)?;
        let transformed = transformed_residual(pair);
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
        let (load, (ks, gs)) = restored_load(k_scale, g_scale, pair.value)?;
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
        // Verify the actual rounded physical load against the original pencil.
        // A precise transformed root cannot justify a coarsely restored load.
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
    check_cancelled(cancelled)?;
    Ok(BucklingResult {
        loads: order.iter().map(|&i| loads[i]).collect(),
        modes: order.iter().map(|&i| modes[i].clone()).collect(),
        relative_residuals: order.iter().map(|&i| relative_residuals[i]).collect(),
        dof_map: dof_map.to_vec(),
        iters,
    })
}

fn check_cancelled(cancelled: &impl Fn() -> bool) -> Result<(), SolidError> {
    if cancelled() {
        Err(SolidError::Cancelled)
    } else {
        Ok(())
    }
}

/// Restore a physical load and coefficients for its actual rounded pencil.
///
/// Splitting binary exponents prevents a subnormal or overflowing intermediate
/// k_scale/g_scale from corrupting an otherwise representable load. The final
/// load is rounded once into its physical exponent range. Crucially, the
/// residual coefficients are then reconstructed from that returned value, not
/// from the unrounded transformed Ritz root.
///
/// The returned coefficients are proportional to (k_scale, load*g_scale), with
/// a common power of two removed. Their relative rounding must be at most
/// 1e-12, leaving headroom for the physical 1e-8 residual gate. A coarse
/// subnormal load may still be returned here: that residual gate decides
/// whether its actual physical equation is accurate enough to publish.
#[allow(clippy::too_many_lines)]
fn restored_load(k_scale: f64, g_scale: f64, theta: f64) -> Result<(f64, (f64, f64)), SolidError> {
    let invalid = |what: &str| SolidError::InvalidInput { what: what.into() };
    if [k_scale, g_scale, theta]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(invalid(
            "buckling load restoration requires finite positive scales and root",
        ));
    }

    // Exact positive x = mantissa * 2^exponent, 1 <= mantissa < 2.
    // Lifting a subnormal by 2^54 is exact and makes its exponent explicit.
    fn parts(x: f64) -> (f64, i32) {
        let (lifted, correction) = if x < f64::MIN_POSITIVE {
            (x * f64::from_bits(1077_u64 << 52), -54)
        } else {
            (x, 0)
        };
        let bits = lifted.to_bits();
        let exponent = ((bits >> 52) & 0x7ff) as i32 - 1023 + correction;
        let mantissa = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | (1023_u64 << 52));
        (mantissa, exponent)
    }

    // A subnormal has only the final multiplication's rounding, not a series
    // of halvings or an already-underflowed power-of-two intermediate.
    fn scale(mantissa: f64, exponent: i32) -> f64 {
        let (mantissa, own) = parts(mantissa);
        let exponent = exponent + own;
        if exponent > 1023 {
            f64::INFINITY
        } else if exponent < -1075 {
            0.0
        } else if exponent >= -1022 {
            mantissa * f64::from_bits(((exponent + 1023) as u64) << 52)
        } else if exponent == -1075 {
            (mantissa * 0.5) * f64::from_bits(1)
        } else {
            mantissa * f64::from_bits(1_u64 << (exponent + 1074))
        }
    }

    fn coefficient(mantissa: f64, exponent: i32) -> Option<f64> {
        let value = scale(mantissa, exponent);
        if !value.is_finite() || value <= 0.0 {
            return None;
        }
        // Compare in the normalized range, even when value is subnormal.
        // Positive rounding can only retain or raise the leading exponent by
        // one here; reconstructing the mantissa never forms an extreme ratio.
        let (rounded, rounded_exponent) = parts(value);
        let delta = rounded_exponent - exponent;
        if !(0..=1).contains(&delta) {
            return None;
        }
        let recovered = rounded * f64::from_bits(((delta + 1023) as u64) << 52);
        ((recovered - mantissa).abs() <= 1e-12 * mantissa).then_some(value)
    }

    let (km, ke) = parts(k_scale);
    let (gm, ge) = parts(g_scale);
    let (tm, te) = parts(theta);
    let load = scale(km / gm / tm, ke - ge - te);
    if !load.is_finite() || load <= 0.0 {
        return Err(invalid(
            "critical buckling load factor is not representable",
        ));
    }

    let (lm, le) = parts(load);
    let (right_mantissa, right_adjustment) = parts(lm * gm);
    let right_exponent = le + ge + right_adjustment;
    let common = ke.max(right_exponent);
    let ks = coefficient(km, ke - common);
    let gs = coefficient(right_mantissa, right_exponent - common);
    match (ks, gs) {
        (Some(ks), Some(gs)) => Ok((load, (ks, gs))),
        _ => Err(invalid(
            "returned buckling load has no faithfully representable residual coefficients",
        )),
    }
}

fn transformed_residual(pair: &EigenPair) -> f64 {
    if pair.value == 0.0 && pair.residual == 0.0 {
        0.0
    } else {
        pair.residual / pair.value.abs()
    }
}

/// Remove physical units without changing sparsity or hiding asymmetry.
fn scaled_symmetric(matrix: &Csr, cancelled: &impl Fn() -> bool) -> Result<(f64, Csr), SolidError> {
    let mut scale = 0.0_f64;
    for row in 0..matrix.nrows() {
        if row % 256 == 0 {
            check_cancelled(cancelled)?;
        }
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
            if row % 256 == 0 {
                check_cancelled(cancelled)?;
            }
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
    cancelled: &impl Fn() -> bool,
) -> Result<Vec<EigenPair>, SolidError> {
    let mut matrix = vec![0.0; n * n];
    let mut unit = vec![0.0; n];
    let mut column = vec![0.0; n];
    for j in 0..n {
        check_cancelled(cancelled)?;
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
    check_cancelled(cancelled)?;
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
