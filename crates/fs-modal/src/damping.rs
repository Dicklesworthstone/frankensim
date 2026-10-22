//! Certified modal damping (IRP M0.2-W5): when a lightly damped body may
//! keep one decay rate per mode, and which modes need a full block.
//!
//! **Setting.** Work in mass-normalized modal coordinates. The undamped
//! frequencies are `Ω = diag(ω_i)` and the symmetric modal damping matrix is
//! `C` [1/s]. The free-decay poles are the 2n roots of
//! `det Q(λ) = 0`, where `Q(λ) = λ²I + λC + Ω²`.
//!
//! Split `C = D + F` into its diagonal and its off-diagonals. The diagonal
//! (energy-quotient) estimate of mode i solves `λ² + λD_ii + ω_i² = 0`:
//! `λ̂_i = −D_ii/2 + iν_i`, with `ν_i = √((ω_i − D_ii/2)(ω_i + D_ii/2))` and
//! `ζ_i = D_ii/(2ω_i)`.
//!
//! **Certificate (Schur complement and Rouché).** Partition the modes into
//! a block B and its complement K. Let `Λ_k(λ) = λ² + λD_kk + ω_k²`; it
//! keeps each mode's own damping exact. Then `Q_KK = Λ_K + λF_KK` and
//! `det Q = det Q_KK · det S_B`, where
//! `S_B(λ) = Q_BB(λ) − λ² F_BK Q_KK(λ)⁻¹ F_KB`.
//!
//! On a disc `𝒟 = {|λ − c| ≤ ρ}` in the open upper half-plane:
//! 1. `|Λ_k(λ)| ≥ m_k(ρ) = (|c − λ̂_k| − ρ)(|c − λ̂_k*| − ρ)`, and both
//!    factors must be positive, so no complement pole lies in 𝒟.
//! 2. `|λ| ≤ a = |c| + ρ`.
//! 3. `Q_KK = Λ^{1/2}(I + λE)Λ^{1/2}` with
//!    `|E_jk| ≤ Ē_jk = |F_jk|/√(m_j m_k)`. Since Ē is symmetric and
//!    nonnegative, `‖E‖₂ ≤ ρ(Ē)`, and a Collatz–Wielandt quotient bounds
//!    that from above. If `a·e < 1`, Q_KK is invertible on all of 𝒟.
//! 4. `‖λ² F_BK Q_KK⁻¹ F_KB‖₂ ≤ P(ρ) = a² S_B(ρ)/(1 − a·e(ρ))`, where
//!    `S_B(ρ) = Σ_{j∈B, k∈K} F_jk²/m_k(ρ)`.
//!
//! **Singletons.** For `B = {i}` and `c = λ̂_i`, on the circle
//! `|q_i(λ)| = ρ|λ − λ̂_i*| ≥ ρ(2ν_i − ρ)`. If `P(ρ) < ρ(2ν_i − ρ)`, Rouché
//! puts exactly one eigenvalue of Q in 𝒟. P is quadratic in F, so the
//! smallest such ρ is second order: `ρ ≈ ω_i² Σ_k F_ik²/(m_k · 2ν_i)`,
//! with `m_k ≈ 2ω_i · gap_ik`. The smallest ρ is the least fixed point of the
//! increasing map `ρ ↦ P(ρ)/(2ν_i − ρ)`. It is found by iteration from 0,
//! inflated by 1e-6, and then *verified*: the verification is the
//! certificate.
//!
//! **Blocks.** Modes whose radius breaks the budget merge, union-find style,
//! with their strongest partner (largest `|C_ij|/|ω_i − ω_j|`; the nearest
//! frequency when no coupling is known). The block's 2b roots `μ_l` come
//! from a dense complex solve of its energy-coordinate matrix
//! `[[0, Ω_B], [−Ω_B, −C_BB]]`. That matrix has the pencil's poles and a
//! norm of order `ω`, not `ω²`.
//!
//! Take a circle around a group G of those roots. Then
//! `‖Q_BB⁻¹‖ ≤ X^{b−1}/Δ`, where `X ≥ ‖Q_BB‖_F` on 𝒟 and
//! `Δ = Π_{l∈G}(ρ − |c − μ_l|) · Π_{l∉G}(|c − μ_l| − ρ)` bounds `|det Q_BB|`.
//! If `X^{b−1} P < Δ`, then det Q_BB and det S_B have the same zero count in
//! 𝒟 (Gohberg–Sigal): exactly |G|.
//!
//! Every certified disc lies in the open upper half-plane, and the discs are
//! pairwise disjoint (overlapping discs merge their blocks). The disc counts
//! therefore sum to n, and every pole is enclosed, with its conjugate in the
//! mirrored disc. A `Full` block is realized by the owner's existing
//! channels restricted to it; this module adds no projection.
//!
//! **Scalable path.** [`certify_modal_damping_windowed`] never forms C.
//! - It reads `C_jk` only for `|ω_j − ω_k| < window`
//!   ([`damping_window_neighbours`]).
//! - It bounds the rest of each row by Bessel's inequality:
//!   `Σ_far C_jk² ≤ ‖C e_j‖² − Σ_near C_jk²`, with the caller's column bound
//!   `‖C e_j‖ ≤ √(xᵀM⁻¹x)`, `x = C_FE φ̂_j`. A difference of squares keeps
//!   only half the digits, so this bound never falls below about `√ε‖C e_j‖`.
//! - [`certify_modal_damping_windowed_far_rows`] instead takes the far norm
//!   itself. Its natural source is the mass norm of the residual
//!   `r_j = M⁻¹x_j − Σ_{k∈near∪{j}} φ̂_k C_kj`: for M-orthonormal modes,
//!   `C_kj = ⟨φ̂_k, r_j⟩_M` for every unread k, so `Σ_far C_kj² ≤ ‖r_j‖²_M`.
//!   The residual is formed by subtracting vectors, not squares, and its
//!   floor is about `ε‖C e_j‖`. [`residual_far_row_bounds`] forms it for a
//!   lumped mass, charging the measured defect of ΦᵀMΦ = I and the rounding.
//! - Splitting the complement into modes within `window/2` of the disc and the
//!   rest gives `e ≤ max(‖Ē_NN‖, ‖Ē_FF‖) + ‖Ē_NF‖`. Only the first term needs
//!   matrix entries.
//!
//! **Floating point.** Comparisons carry a relative safety margin
//! ([`ROUCHE_SAFETY`]); this is not interval arithmetic.
//! - Singleton centres are closed-form roots within `δ = 2εν` of the exact
//!   ones. That δ enters both `|q_i|` on the circle and every `m_k`, so
//!   radii never fall below the centres' own rounding.
//! - Block roots are seeded by the dense solve. They are then refined by
//!   inverse iteration and the two-sided Rayleigh functional, and carry 8×
//!   their first-order forward-error estimate into Δ. That is an
//!   a-posteriori estimate, not an enclosure, so an ill-conditioned block can
//!   break it. The dense-truth tests are the check.
//!
//! **Counted once.** Only material damping is certified here. Radiation
//! and thermoelastic losses remain separate owners.

use crate::ModalError;
use fs_la::eigen_complex::{eig as complex_eig, lu_complex};
use fs_math::c64::C64;
use fs_math::det;
use std::collections::BTreeMap;

/// Largest certified cluster: the 64-slot modal connection ceiling.
pub const MAX_DAMPING_CLUSTER: usize = 64;

/// Relative margin by which every Rouché comparison must hold.
pub const ROUCHE_SAFETY: f64 = 1e-9;

/// Fixed-point iterations allowed per disc before declaring it uncertified.
const MAX_RADIUS_ITERATIONS: usize = 200;
/// Collatz–Wielandt power steps for the Perron bound.
const PERRON_STEPS: usize = 6;

/// Required accuracy of the certified eigenvalues. There is no default:
/// every field must be finite and positive, and the cluster cap must lie in
/// `1..=MAX_DAMPING_CLUSTER`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModalDampingBudget {
    /// Largest certified `|Δ Re λ| / |Re λ|`: the relative decay-rate error.
    pub max_rel_decay_error: f64,
    /// Largest certified frequency error in cents, from `|Δ Im λ|/Im λ`.
    pub max_freq_error_cents: f64,
    /// Largest block the fallback may form.
    pub max_cluster_size: usize,
}

/// One certified piece of the modal damping.
#[derive(Debug, Clone, PartialEq)]
pub enum ModalDampingBlock {
    /// One decay rate per mode is certified: the exact pole lies within
    /// `radius` of `eigenvalue = −ζω + iν`.
    Diagonal {
        /// Mode index.
        mode: usize,
        /// Energy-quotient damping ratio `D_ii/(2ω_i)`.
        zeta: f64,
        /// Diagonal estimate `λ̂_i` (upper half-plane).
        eigenvalue: C64,
        /// Certified radius [1/s].
        radius: f64,
    },
    /// The modes must be integrated as one block. `eigenvalues` are the
    /// block's upper-half-plane poles, each within `radii[l]` of a true pole.
    Full {
        /// Mode indices, ascending.
        modes: Vec<usize>,
        /// Block poles, upper half-plane, ascending imaginary part.
        eigenvalues: Vec<C64>,
        /// Certified radius per pole [1/s].
        radii: Vec<f64>,
        /// Largest of `radii`.
        radius: f64,
    },
}

/// Far coupling used by the certificate. The dense path has no window and
/// zero far norms.
#[derive(Debug, Clone, PartialEq)]
pub struct FarCouplingBound {
    /// Coupling window [rad/s], or `None` for the dense path.
    pub window_rad_s: Option<f64>,
    /// Upper bound on the norm of each row's off-diagonal entries outside the
    /// window (Bessel).
    pub far_row_norm: Vec<f64>,
}

/// The certified decomposition of a modal damping matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct CertifiedModalDamping {
    /// Blocks in ascending order of their first mode.
    pub blocks: Vec<ModalDampingBlock>,
    /// Largest certified relative decay error over all poles.
    pub worst_rel_decay: f64,
    /// Largest certified frequency error over all poles [cents].
    pub worst_freq_cents: f64,
    /// Far-coupling bound used.
    pub far_coupling: FarCouplingBound,
}

/// Certify the modal damping of a dense symmetric `C` (row-major, n×n) on
/// mass-normalized modes with undamped frequencies `omega` [rad/s]. This is
/// the reference path.
///
/// # Errors
/// - [`ModalError::DampingInput`] for a bad budget, frequency, shape,
///   non-finite or asymmetric entry, or negative diagonal.
/// - [`ModalError::DampingNotUnderdamped`] when a diagonal estimate or a
///   block has real poles.
/// - [`ModalError::DampingClusterTooLarge`] when the minimal certified
///   cluster exceeds the cap.
/// - [`ModalError::QuadraticEig`] when a block solve fails.
pub fn certify_modal_damping(
    omega: &[f64],
    c: &[f64],
    budget: &ModalDampingBudget,
) -> Result<CertifiedModalDamping, ModalError> {
    let n = omega.len();
    if n == 0 || c.len() != n * n {
        return Err(input("damping matrix shape", c.len()));
    }
    validate(omega, budget)?;
    let mut diag = Vec::with_capacity(n);
    let mut near = Vec::with_capacity(n);
    for i in 0..n {
        let mut row = Vec::new();
        for j in 0..n {
            let value = c[n * i + j];
            if !value.is_finite() {
                return Err(input("non-finite damping entry", n * i + j));
            }
            // Exact symmetry; finite values differ iff their difference does.
            if value - c[n * j + i] != 0.0 {
                return Err(input("asymmetric damping entry", n * i + j));
            }
            if j != i && value != 0.0 {
                row.push((j, value));
            }
        }
        diag.push(c[n * i + i]);
        near.push(row);
    }
    let coupling = Coupling::new(diag, near, vec![0.0; n], None);
    Certifier::new(omega, coupling, budget)?.run()
}

/// Certify without forming C. `near(i, j)` returns `C_ij`. It is called only
/// for `|ω_i − ω_j| < window`, including `i == j`. `column_norm_bound[j]`
/// must bound `‖C e_j‖₂` over every mode of the underlying space, not only
/// the retained ones; Bessel's inequality then bounds the unread entries.
///
/// # Errors
/// - As [`certify_modal_damping`].
/// - [`ModalError::DampingInput`] for a bad window, a column bound below
///   the entries actually read, or an asymmetric `near`.
/// - [`ModalError::DampingBlockSpansWindow`] when a required block needs
///   entries outside the window.
pub fn certify_modal_damping_windowed(
    omega: &[f64],
    near: &dyn Fn(usize, usize) -> f64,
    column_norm_bound: &[f64],
    window: f64,
    budget: &ModalDampingBudget,
) -> Result<CertifiedModalDamping, ModalError> {
    windowed(
        omega,
        near,
        FarInput::ColumnNorm(column_norm_bound),
        window,
        budget,
    )
}

/// Certify without forming C, with the caller's own bound on what is not
/// read. `near` is as in [`certify_modal_damping_windowed`].
/// `far_row_norm_bound[j]` must bound the norm of row j's off-diagonal
/// entries outside [`damping_window_neighbours`]`(omega, window)[j]`. The
/// mass norm of the residual `M⁻¹x_j − Σ_{k∈near∪{j}} φ̂_k C_kj` is such a
/// bound for M-orthonormal modes, over every mode of the underlying space,
/// without the √ε floor of a difference of squares.
///
/// # Errors
/// - As [`certify_modal_damping`].
/// - [`ModalError::DampingInput`] for a bad window, a negative or
///   non-finite far bound, or an asymmetric `near`.
/// - [`ModalError::DampingBlockSpansWindow`] when a required block needs
///   entries outside the window.
pub fn certify_modal_damping_windowed_far_rows(
    omega: &[f64],
    near: &dyn Fn(usize, usize) -> f64,
    far_row_norm_bound: &[f64],
    window: f64,
    budget: &ModalDampingBudget,
) -> Result<CertifiedModalDamping, ModalError> {
    windowed(
        omega,
        near,
        FarInput::FarRows(far_row_norm_bound),
        window,
        budget,
    )
}

/// The off-diagonal entries a windowed certificate reads: for each mode i,
/// every other mode k with `|ω_i − ω_k| < window`, ascending. The relation
/// is symmetric. A caller that bounds the unread entries itself must exclude
/// exactly these.
///
/// # Errors
/// [`ModalError::DampingInput`] for a non-positive or non-finite window or
/// frequency.
pub fn damping_window_neighbours(
    omega: &[f64],
    window: f64,
) -> Result<Vec<Vec<usize>>, ModalError> {
    if !window.is_finite() || window <= 0.0 {
        return Err(input("coupling window", 0));
    }
    for (i, &w) in omega.iter().enumerate() {
        if !(w.is_finite() && w > 0.0) {
            return Err(input("undamped frequency", i));
        }
    }
    let n = omega.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| omega[a].total_cmp(&omega[b]).then(a.cmp(&b)));
    let mut neighbours = vec![Vec::new(); n];
    for (position, &i) in order.iter().enumerate() {
        // Sorted sweep: only modes inside the window are ever visited.
        let low = order[..position]
            .iter()
            .rev()
            .take_while(|&&j| omega[i] - omega[j] < window);
        let high = order[position + 1..]
            .iter()
            .take_while(|&&j| omega[j] - omega[i] < window);
        let mut row: Vec<usize> = low.chain(high).copied().collect();
        row.sort_unstable();
        neighbours[i] = row;
    }
    Ok(neighbours)
}

/// Largest radius around `mu` that still meets the budget.
fn budget_radius(budget: &ModalDampingBudget, mu: C64) -> f64 {
    let k = budget.max_freq_error_cents * core::f64::consts::LN_2 / 1200.0;
    (budget.max_rel_decay_error * mu.re.abs()).min(mu.im * k / (1.0 + k))
}

/// [`certify_modal_damping`] with the diagonal decay rates owned by the
/// caller. `zeta_owner[i]` is the caller's own evaluator for mode i (for a
/// plate, its energy-quotient image), which can differ from `C_ii/(2ω_i)` in
/// its last bits.
///
/// The certificate is computed on C exactly as in [`certify_modal_damping`].
/// Each `Diagonal` block then reports `zeta_owner[i]` and its closed-form pole
/// `λ̂ᵒ = −ζω + i√((ω − ζω)(ω + ζω))`. Its radius grows by `|λ̂ᵒ − λ̂|` plus
/// the new centre's rounding, so by the triangle inequality the disc still
/// encloses the exact pole. `Full` blocks are unchanged.
///
/// # Errors
/// - As [`certify_modal_damping`].
/// - [`ModalError::DampingInput`] for a wrong length or a negative or
///   non-finite owner ζ, for an owner ζ whose grown disc breaks the budget,
///   or for one that reaches another certified disc.
/// - [`ModalError::DampingNotUnderdamped`] when an owner ζ is not underdamped.
pub fn certify_modal_damping_owned_diagonal(
    omega: &[f64],
    c: &[f64],
    zeta_owner: &[f64],
    budget: &ModalDampingBudget,
) -> Result<CertifiedModalDamping, ModalError> {
    if zeta_owner.len() != omega.len() {
        return Err(input("owner zeta length", zeta_owner.len()));
    }
    if let Some(i) = zeta_owner
        .iter()
        .position(|z| !(z.is_finite() && *z >= 0.0))
    {
        return Err(input("owner zeta", i));
    }
    let mut result = certify_modal_damping(omega, c, budget)?;
    for block in &mut result.blocks {
        if let ModalDampingBlock::Diagonal {
            mode,
            zeta,
            eigenvalue,
            radius,
        } = block
        {
            let (w, z) = (omega[*mode], zeta_owner[*mode]);
            let half = z * w;
            if half >= w {
                return Err(ModalError::DampingNotUnderdamped { modes: vec![*mode] });
            }
            let nu = det::sqrt((w - half) * (w + half));
            let centre = C64::new(-half, nu);
            // The owner's closed form: one product for the real part, and
            // as in `Certifier::new`, 2εν for the imaginary part.
            let rounding = f64::EPSILON * half + 2.0 * f64::EPSILON * nu;
            let grown = (*radius + (centre - *eigenvalue).abs() + rounding) * (1.0 + ROUCHE_SAFETY);
            if !(grown <= budget_radius(budget, centre)) {
                return Err(input(
                    "owner zeta disagrees with the damping diagonal beyond the budget",
                    *mode,
                ));
            }
            (*zeta, *eigenvalue, *radius) = (z, centre, grown);
        }
    }
    // Every pole's disc, as (centre, radius, block). A `Full` pole's disc
    // contains its group's certified disc, so only pairs with a moved
    // `Diagonal` disc are checked: the certifier made the rest disjoint.
    let mut discs: Vec<(C64, f64, usize)> = Vec::new();
    for (b, block) in result.blocks.iter().enumerate() {
        match block {
            ModalDampingBlock::Diagonal {
                eigenvalue, radius, ..
            } => discs.push((*eigenvalue, *radius, b)),
            ModalDampingBlock::Full {
                eigenvalues, radii, ..
            } => discs.extend(eigenvalues.iter().zip(radii).map(|(&e, &r)| (e, r, b))),
        }
    }
    for (x, &(cx, rx, bx)) in discs.iter().enumerate() {
        for &(cy, ry, by) in &discs[x + 1..] {
            let moved = |b: usize| matches!(result.blocks[b], ModalDampingBlock::Diagonal { .. });
            if bx != by && (moved(bx) || moved(by)) && (cx - cy).abs() <= rx + ry {
                let mode = match &result.blocks[bx] {
                    ModalDampingBlock::Diagonal { mode, .. } => *mode,
                    ModalDampingBlock::Full { modes, .. } => modes[0],
                };
                return Err(input(
                    "owner zeta moves a certified disc into another",
                    mode,
                ));
            }
        }
    }
    let (mut worst_rel_decay, mut worst_freq_cents) = (0.0_f64, 0.0_f64);
    for &(centre, radius, _) in &discs {
        if radius > 0.0 {
            worst_rel_decay = worst_rel_decay.max(radius / centre.re.abs());
            worst_freq_cents = worst_freq_cents.max(cents_bound(radius, centre.im));
        }
    }
    result.worst_rel_decay = worst_rel_decay;
    result.worst_freq_cents = worst_freq_cents;
    Ok(result)
}

/// Window entries and Bessel bounds on the rest of each row, from
/// [`residual_far_row_bounds`].
#[derive(Debug, Clone, PartialEq)]
pub struct ResidualFarBounds {
    /// Row j holds `(k, C_jk)` for k = j and every window neighbour,
    /// ascending. Each pair is formed once, so `C_jk` and `C_kj` are one value.
    pub near: Vec<Vec<(usize, f64)>>,
    /// Bound on the norm of each row's unformed off-diagonal entries.
    pub far_row_norm: Vec<f64>,
    /// Bound on each column norm `‖C e_j‖`.
    pub column_norm: Vec<f64>,
    /// Bound on `‖ΦᵀMΦ − I‖₂` of the given modes, charged to both bounds.
    pub basis_defect: f64,
}

/// Form `C_kj = φ_kᵀx_j` inside the window and bound the rest of each row
/// from the residual `r_j = M⁻¹x_j − Σ_{k∈near∪{j}} φ_k C_kj`.
///
/// - `mass` is the diagonal of a lumped positive mass matrix.
/// - `modes` are the retained shapes φ_j on the same DOFs, mass-normalized.
///   They need not be exactly M-orthonormal: the measured defect
///   `E = ΦᵀMΦ − I` is charged.
/// - `applied[j]` is `x_j = C_FE φ_j` as the caller computed it, and
///   `applied_error[j]` bounds its distance from the exact product in the
///   M⁻¹ norm.
///
/// For every unformed retained k, and for ANY near coefficients,
/// `C_kj = ⟨φ_k, r_j⟩_M + Σ_l E_kl C_lj − φ_kᵀδx_j`. Hence
/// `‖far_j‖ ≤ √(1+‖E‖)(‖r_j‖_M + ‖δx_j‖_{M⁻¹}) + ‖E‖‖C_near,j‖`, and Bessel
/// covers every mode of the underlying space, retained or not. The rounding
/// of r_j, of the Gram sums and of the norms enters as first-order `γ_k`
/// allowances; this is not interval arithmetic. Cost: O(n²·dofs) for the
/// defect plus O(Σ|near|·dofs).
///
/// # Errors
/// [`ModalError::DampingInput`] for mismatched shapes, a bad frequency,
/// window or mass, a non-finite mode or applied vector, a negative or
/// non-finite error bound, or a bound that overflows.
pub fn residual_far_row_bounds(
    omega: &[f64],
    window: f64,
    mass: &[f64],
    modes: &[Vec<f64>],
    applied: &[Vec<f64>],
    applied_error: &[f64],
) -> Result<ResidualFarBounds, ModalError> {
    let n = omega.len();
    let dofs = mass.len();
    if n == 0 || modes.len() != n || applied.len() != n || applied_error.len() != n {
        return Err(input("residual bound shape", n));
    }
    let neighbours = damping_window_neighbours(omega, window)?;
    for (d, &m) in mass.iter().enumerate() {
        if !(m.is_finite() && m > 0.0) {
            return Err(input("lumped mass", d));
        }
    }
    for j in 0..n {
        if modes[j].len() != dofs || applied[j].len() != dofs {
            return Err(input("residual bound shape", j));
        }
        if modes[j].iter().chain(&applied[j]).any(|v| !v.is_finite()) {
            return Err(input("non-finite mode or applied vector", j));
        }
        if !(applied_error[j].is_finite() && applied_error[j] >= 0.0) {
            return Err(input("applied error bound", j));
        }
    }
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(p, q)| p * q).sum::<f64>();
    let mass_norm = |v: &[f64]| det::sqrt(v.iter().zip(mass).map(|(p, m)| m * p * p).sum::<f64>());
    let inverse_norm =
        |v: &[f64]| det::sqrt(v.iter().zip(mass).map(|(p, m)| p * p / m).sum::<f64>());
    let mut near: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for j in 0..n {
        near[j].push((j, dot(&modes[j], &applied[j])));
        for &k in neighbours[j].iter().filter(|&&k| k > j) {
            let value = dot(&modes[k], &applied[j]);
            near[j].push((k, value));
            near[k].push((j, value));
        }
    }
    for row in &mut near {
        row.sort_by_key(|&(k, _)| k);
    }
    // ‖E‖_F ≥ ‖E‖₂. Each Gram sum is within γ‖φ_k‖_M‖φ_l‖_M ≤ γ(1 + ‖E‖)
    // of exact, so ‖E‖ ≤ (g + nγ)/(1 − nγ) for the computed Frobenius norm g.
    let mut defect_sq = 0.0_f64;
    for k in 0..n {
        for l in k..n {
            let gram: f64 = modes[k]
                .iter()
                .zip(&modes[l])
                .zip(mass)
                .map(|((p, q), m)| p * m * q)
                .sum();
            let e = gram - if k == l { 1.0 } else { 0.0 };
            defect_sq += if k == l { e * e } else { 2.0 * e * e };
        }
    }
    let spread = n as f64 * gamma(dofs + 2);
    if !(spread < 0.5) {
        return Err(input("residual bound size", n));
    }
    let g = det::sqrt(defect_sq) * (1.0 + gamma(n * n + 2));
    let basis_defect = (g + spread) / (1.0 - spread);
    let inflate = det::sqrt(1.0 + basis_defect);
    // A computed norm is within γ_{dofs+2} of exact (positive terms, sqrt).
    let norms = 1.0 + gamma(dofs + 2);
    let mut far = Vec::with_capacity(n);
    let mut column = Vec::with_capacity(n);
    for j in 0..n {
        let mut r: Vec<f64> = applied[j].iter().zip(mass).map(|(x, m)| x / m).collect();
        let mut magnitude: Vec<f64> = r.iter().map(|v| v.abs()).collect();
        let mut near_sq = 0.0_f64;
        for &(k, value) in &near[j] {
            near_sq = value.mul_add(value, near_sq);
            for ((r, size), &p) in r.iter_mut().zip(&mut magnitude).zip(&modes[k]) {
                let t = p * value;
                *r -= t;
                *size += t.abs();
            }
        }
        // Each residual entry took one division and |near| products and
        // subtractions: within γ_{2|near|+2} of its magnitude sum.
        let residual_rounding = gamma(2 * near[j].len() + 2) * mass_norm(&magnitude) * norms;
        column.push(inflate * (inverse_norm(&applied[j]) * norms + applied_error[j]));
        far.push(
            inflate * (mass_norm(&r) * norms + residual_rounding + applied_error[j])
                + basis_defect * det::sqrt(near_sq) * norms,
        );
    }
    if let Some(j) = (0..n).find(|&j| {
        !(far[j].is_finite() && column[j].is_finite() && near[j].iter().all(|(_, v)| v.is_finite()))
    }) {
        return Err(input("residual bound overflow", j));
    }
    Ok(ResidualFarBounds {
        near,
        far_row_norm: far,
        column_norm: column,
        basis_defect,
    })
}

/// `γ_k = kε/(1 − kε)`: the relative error bound of k rounded operations.
fn gamma(k: usize) -> f64 {
    let ke = k as f64 * f64::EPSILON;
    ke / (1.0 - ke)
}

/// How a windowed certificate bounds the entries it does not read.
#[derive(Clone, Copy)]
enum FarInput<'a> {
    /// `‖C e_j‖₂` over the complete space; Bessel by a difference of squares.
    ColumnNorm(&'a [f64]),
    /// The caller's bound on each row's unread off-diagonal norm.
    FarRows(&'a [f64]),
}

fn windowed(
    omega: &[f64],
    near: &dyn Fn(usize, usize) -> f64,
    far_input: FarInput<'_>,
    window: f64,
    budget: &ModalDampingBudget,
) -> Result<CertifiedModalDamping, ModalError> {
    let n = omega.len();
    let (supplied, length_field) = match far_input {
        FarInput::ColumnNorm(bounds) => (bounds, "column bound length"),
        FarInput::FarRows(bounds) => (bounds, "far row bound length"),
    };
    if n == 0 || supplied.len() != n {
        return Err(input(length_field, supplied.len()));
    }
    let neighbours = damping_window_neighbours(omega, window)?;
    validate(omega, budget)?;
    let mut diag = vec![0.0; n];
    let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut far = vec![0.0; n];
    for i in 0..n {
        let value = near(i, i);
        if !value.is_finite() {
            return Err(input("non-finite damping entry", n * i + i));
        }
        diag[i] = value;
        let mut sum_sq = value * value;
        for &j in &neighbours[i] {
            let value = near(i, j);
            if !value.is_finite() {
                return Err(input("non-finite damping entry", n * i + j));
            }
            if value - near(j, i) != 0.0 {
                return Err(input("asymmetric damping entry", n * i + j));
            }
            if value != 0.0 {
                rows[i].push((j, value));
            }
            sum_sq = value.mul_add(value, sum_sq);
        }
        let bound = supplied[i];
        far[i] = match far_input {
            FarInput::ColumnNorm(_) => {
                if !bound.is_finite() || bound < 0.0 || bound * bound < sum_sq * (1.0 - 1e-12) {
                    return Err(input("column norm bound below the entries read", i));
                }
                // Bessel's difference, inflated for the rounding of the sums
                // above (a difference of squares loses half the digits: the
                // floor is about √ε·bound, never a vanishing bound).
                let terms = (rows[i].len() + 2) as f64;
                let rounding = 4.0 * terms * f64::EPSILON * bound * bound;
                det::sqrt((bound * bound - sum_sq).max(0.0) + rounding)
            }
            FarInput::FarRows(_) => {
                if !bound.is_finite() || bound < 0.0 {
                    return Err(input("far row norm bound", i));
                }
                bound
            }
        };
    }
    let coupling = Coupling::new(diag, rows, far, Some(window));
    Certifier::new(omega, coupling, budget)?.run()
}

fn input(field: &'static str, index: usize) -> ModalError {
    ModalError::DampingInput { field, index }
}

fn validate(omega: &[f64], budget: &ModalDampingBudget) -> Result<(), ModalError> {
    if !(budget.max_rel_decay_error.is_finite() && budget.max_rel_decay_error > 0.0) {
        return Err(input("budget max_rel_decay_error", 0));
    }
    if !(budget.max_freq_error_cents.is_finite() && budget.max_freq_error_cents > 0.0) {
        return Err(input("budget max_freq_error_cents", 0));
    }
    if !(1..=MAX_DAMPING_CLUSTER).contains(&budget.max_cluster_size) {
        return Err(input("budget max_cluster_size", budget.max_cluster_size));
    }
    for (i, &w) in omega.iter().enumerate() {
        if !(w.is_finite() && w > 0.0) {
            return Err(input("undamped frequency", i));
        }
    }
    Ok(())
}

/// Known off-diagonal entries and bounds on the unknown ones.
struct Coupling {
    diag: Vec<f64>,
    /// Known nonzero off-diagonal entries of each row, ascending column.
    near: Vec<Vec<(usize, f64)>>,
    /// Upper bound on the norm of each row's unknown off-diagonal entries.
    far: Vec<f64>,
    /// Upper bound on each row's off-diagonal 1-norm.
    l1: Vec<f64>,
    window: Option<f64>,
}

impl Coupling {
    fn new(
        diag: Vec<f64>,
        near: Vec<Vec<(usize, f64)>>,
        far: Vec<f64>,
        window: Option<f64>,
    ) -> Self {
        let n = diag.len();
        let l1 = near
            .iter()
            .zip(&far)
            .map(|(row, &far)| {
                let known: f64 = row.iter().map(|&(_, v)| v.abs()).sum();
                let unknown = n.saturating_sub(1 + row.len()) as f64;
                // Cauchy–Schwarz over the unread entries.
                known + far * det::sqrt(unknown)
            })
            .collect();
        Self {
            diag,
            near,
            far,
            l1,
            window,
        }
    }

    fn entry(&self, i: usize, j: usize) -> f64 {
        if i == j {
            return self.diag[i];
        }
        self.near[i]
            .binary_search_by_key(&j, |&(k, _)| k)
            .map_or(0.0, |at| self.near[i][at].1)
    }
}

/// A certified disc: `count` poles of the full pencil lie inside it.
#[derive(Debug, Clone)]
struct Disc {
    center: C64,
    radius: f64,
    /// Block poles inside, each with its own radius `|μ − c| + ρ`.
    members: Vec<(C64, f64)>,
}

#[derive(Clone)]
enum BlockOutcome {
    Certified(Vec<Disc>),
    /// No disc within budget: merge with a partner.
    Uncertified,
}

struct Certifier<'a> {
    omega: &'a [f64],
    coupling: Coupling,
    budget: ModalDampingBudget,
    lam: Vec<C64>,
    nu: Vec<f64>,
    /// Bound on `|λ̂_k − exact root of Λ_k|`: the closed form's rounding.
    delta: Vec<f64>,
}

impl<'a> Certifier<'a> {
    fn new(
        omega: &'a [f64],
        coupling: Coupling,
        budget: &ModalDampingBudget,
    ) -> Result<Self, ModalError> {
        let n = omega.len();
        let mut lam = Vec::with_capacity(n);
        let mut nu = Vec::with_capacity(n);
        let mut delta = Vec::with_capacity(n);
        let mut overdamped = Vec::new();
        for (i, (&d, &w)) in coupling.diag.iter().zip(omega).enumerate() {
            if d < 0.0 {
                return Err(input("negative modal damping diagonal", i));
            }
            let half = 0.5 * d;
            if half >= w {
                overdamped.push(i);
                lam.push(C64::ZERO);
                nu.push(0.0);
                delta.push(0.0);
                continue;
            }
            // Undamped: iω is exact. Otherwise two roundings, a product and
            // a square root put ν within 1.25ε·ν of the exact root.
            let (v, err) = if d == 0.0 {
                (w, 0.0)
            } else {
                let v = det::sqrt((w - half) * (w + half));
                (v, 2.0 * f64::EPSILON * v)
            };
            lam.push(C64::new(-half, v));
            nu.push(v);
            delta.push(err);
        }
        if !overdamped.is_empty() {
            return Err(ModalError::DampingNotUnderdamped { modes: overdamped });
        }
        Ok(Self {
            omega,
            coupling,
            budget: *budget,
            lam,
            nu,
            delta,
        })
    }

    fn n(&self) -> usize {
        self.omega.len()
    }

    /// Largest radius around `mu` that still meets the budget.
    fn budget_radius(&self, mu: C64) -> f64 {
        budget_radius(&self.budget, mu)
    }

    fn run(&self) -> Result<CertifiedModalDamping, ModalError> {
        let n = self.n();
        let mut parent: Vec<usize> = (0..n).collect();
        // A block's complement is every other mode however it is grouped, so
        // an unchanged block keeps its outcome across rounds.
        let mut cache: BTreeMap<Vec<usize>, BlockOutcome> = BTreeMap::new();
        loop {
            let blocks = groups(&mut parent);
            let mut outcomes = Vec::with_capacity(blocks.len());
            for block in &blocks {
                let outcome = if let Some(outcome) = cache.get(block) {
                    outcome.clone()
                } else {
                    let outcome = self.certify_block(block)?;
                    cache.insert(block.clone(), outcome.clone());
                    outcome
                };
                outcomes.push(outcome);
            }
            let mut merges: Vec<(usize, usize)> = Vec::new();
            for (block, outcome) in blocks.iter().zip(&outcomes) {
                if matches!(outcome, BlockOutcome::Uncertified) {
                    merges.push((block[0], self.partner(block)));
                }
            }
            // Discs of different blocks must be disjoint for the counts to
            // add up.
            let discs: Vec<(usize, &Disc)> = blocks
                .iter()
                .zip(&outcomes)
                .filter_map(|(block, outcome)| match outcome {
                    BlockOutcome::Certified(discs) => {
                        Some(discs.iter().map(move |d| (block[0], d)))
                    }
                    BlockOutcome::Uncertified => None,
                })
                .flatten()
                .collect();
            for (a, &(owner_a, disc_a)) in discs.iter().enumerate() {
                for &(owner_b, disc_b) in &discs[a + 1..] {
                    if owner_a != owner_b
                        && (disc_a.center - disc_b.center).abs() <= disc_a.radius + disc_b.radius
                    {
                        merges.push((owner_a, owner_b));
                    }
                }
            }
            if merges.is_empty() {
                return Ok(self.assemble(&blocks, outcomes));
            }
            for (a, b) in merges {
                let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                if ra != rb {
                    let (low, high) = (ra.min(rb), ra.max(rb));
                    parent[high] = low;
                }
            }
            for block in groups(&mut parent) {
                if block.len() > self.budget.max_cluster_size {
                    return Err(ModalError::DampingClusterTooLarge {
                        modes: block,
                        cap: self.budget.max_cluster_size,
                    });
                }
            }
        }
    }

    /// Strongest coupling partner outside `block`, or else the nearest
    /// frequency. Ties go to the lowest index.
    fn partner(&self, block: &[usize]) -> usize {
        let inside = |k: usize| block.binary_search(&k).is_ok();
        let mut best: Option<(f64, usize)> = None;
        for &j in block {
            for &(k, value) in &self.coupling.near[j] {
                if inside(k) {
                    continue;
                }
                let gap = (self.omega[j] - self.omega[k]).abs();
                let score = if gap == 0.0 {
                    f64::INFINITY
                } else {
                    value.abs() / gap
                };
                if best.is_none_or(|(s, b)| match score.total_cmp(&s) {
                    core::cmp::Ordering::Greater => true,
                    core::cmp::Ordering::Equal => k < b,
                    core::cmp::Ordering::Less => false,
                }) {
                    best = Some((score, k));
                }
            }
        }
        if let Some((_, k)) = best {
            return k;
        }
        let mut nearest = (f64::INFINITY, usize::MAX);
        for &j in block {
            for k in (0..self.n()).filter(|&k| !inside(k)) {
                let gap = (self.omega[j] - self.omega[k]).abs();
                if gap.total_cmp(&nearest.0).then(k.cmp(&nearest.1)).is_lt() {
                    nearest = (gap, k);
                }
            }
        }
        nearest.1
    }

    fn certify_block(&self, block: &[usize]) -> Result<BlockOutcome, ModalError> {
        if let Some(window) = self.coupling.window {
            for &j in block {
                for &k in block {
                    if (self.omega[j] - self.omega[k]).abs() >= window {
                        return Err(ModalError::DampingBlockSpansWindow {
                            modes: block.to_vec(),
                            window,
                        });
                    }
                }
            }
        }
        if block.len() == self.n() {
            // No complement: the block poles are the poles, to the accuracy
            // of their refinement.
            return self.block_poles(block).map(|poles| {
                BlockOutcome::Certified(
                    poles
                        .into_iter()
                        .map(|(mu, slack)| Disc {
                            center: mu,
                            radius: slack,
                            members: vec![(mu, slack)],
                        })
                        .collect(),
                )
            });
        }
        if let [i] = *block {
            return Ok(self.certify_singleton(i));
        }
        let poles = self.block_poles(block)?;
        let mut all_roots: Vec<(C64, f64)> = poles.clone();
        all_roots.extend(poles.iter().map(|&(p, s)| (p.conj(), s)));
        // Groups of upper poles that share one disc.
        let mut groups: Vec<Vec<usize>> = (0..poles.len()).map(|l| vec![l]).collect();
        loop {
            let mut discs: Vec<Disc> = Vec::with_capacity(groups.len());
            let mut merge = None;
            for (g, members) in groups.iter().enumerate() {
                match self.block_disc(block, &poles, &all_roots, members) {
                    GroupDisc::Certified(disc) => discs.push(disc),
                    GroupDisc::Uncertified => return Ok(BlockOutcome::Uncertified),
                    GroupDisc::Hits(other) => {
                        match groups.iter().position(|m| m.contains(&other)) {
                            Some(h) if h != g => {
                                merge = Some((g, h));
                                break;
                            }
                            _ => return Ok(BlockOutcome::Uncertified),
                        }
                    }
                }
            }
            if merge.is_none() {
                'overlap: for a in 0..discs.len() {
                    for b in a + 1..discs.len() {
                        if (discs[a].center - discs[b].center).abs()
                            <= discs[a].radius + discs[b].radius
                        {
                            merge = Some((a, b));
                            break 'overlap;
                        }
                    }
                }
            }
            match merge {
                Some((a, b)) => merge_groups(&mut groups, a, b),
                None => return Ok(BlockOutcome::Certified(discs)),
            }
        }
    }

    /// Upper-half-plane poles of the block's own pencil, ascending in
    /// imaginary part, each with an error allowance.
    ///
    /// The dense solve only seeds them: fs-la's complex QR can sit ~1e-12
    /// relative off. Each seed is refined by inverse iteration on the complex
    /// symmetric `Q_BB(λ)` and the root of `xᵀQ_BB(λ)x` (x is also the left
    /// eigenvector). The allowance is 8× the first-order forward error
    /// `‖Q_BB(μ)x‖‖x‖/|xᵀQ_BB′(μ)x|` plus 8ε|μ|. That is an a-posteriori
    /// estimate, not an enclosure. A refinement that wanders or collides
    /// keeps its seed with a 1e-6 relative allowance.
    fn block_poles(&self, block: &[usize]) -> Result<Vec<(C64, f64)>, ModalError> {
        let b = block.len();
        if b == 1 {
            let i = block[0];
            return Ok(vec![(self.lam[i], self.delta[i])]);
        }
        let d = 2 * b;
        let mut a = vec![C64::ZERO; d * d];
        for (r, &j) in block.iter().enumerate() {
            a[r * d + b + r] = C64::from_re(self.omega[j]);
            a[(b + r) * d + r] = C64::from_re(-self.omega[j]);
            for (s, &k) in block.iter().enumerate() {
                a[(b + r) * d + b + s] = C64::from_re(-self.coupling.entry(j, k));
            }
        }
        let roots = complex_eig(&a, d).map_err(|e| ModalError::QuadraticEig {
            detail: format!("certified damping block {block:?}: {e:?}"),
        })?;
        let mut seeds: Vec<C64> = roots.into_iter().filter(|r| r.im > 0.0).collect();
        if seeds.len() != b {
            return Err(ModalError::DampingNotUnderdamped {
                modes: block.to_vec(),
            });
        }
        seeds.sort_by(|x, y| x.im.total_cmp(&y.im).then(x.re.total_cmp(&y.re)));
        let fallback = |seed: C64| (seed, 1e-6 * seed.abs());
        let mut poles: Vec<(C64, f64)> = seeds
            .iter()
            .map(|&seed| {
                let (mu, estimate) = self.refine_block_pole(block, seed);
                if (mu - seed).abs() <= 1e-6 * seed.abs() && estimate.is_finite() {
                    (mu, 8.0 * estimate + 8.0 * f64::EPSILON * mu.abs())
                } else {
                    fallback(seed)
                }
            })
            .collect();
        // Two seeds refined onto one root would lose a root of det Q_BB.
        for x in 0..b {
            for y in x + 1..b {
                let (px, py) = (poles[x], poles[y]);
                if (px.0 - py.0).abs() <= px.1 + py.1 {
                    poles[x] = fallback(seeds[x]);
                    poles[y] = fallback(seeds[y]);
                }
            }
        }
        poles.sort_by(|x, y| x.0.im.total_cmp(&y.0.im).then(x.0.re.total_cmp(&y.0.re)));
        Ok(poles)
    }

    /// `Q_BB(μ) = μ²I + μC_BB + Ω_B²`, row-major.
    fn block_pencil(&self, block: &[usize], mu: C64) -> Vec<C64> {
        let b = block.len();
        let mut q = vec![C64::ZERO; b * b];
        for (r, &j) in block.iter().enumerate() {
            for (s, &k) in block.iter().enumerate() {
                q[r * b + s] = mu * C64::from_re(self.coupling.entry(j, k));
            }
            q[r * b + r] = q[r * b + r] + mu * mu + C64::from_re(self.omega[j] * self.omega[j]);
        }
        q
    }

    /// Inverse iteration plus the two-sided Rayleigh functional; returns the
    /// refined pole and its first-order forward-error estimate.
    fn refine_block_pole(&self, block: &[usize], seed: C64) -> (C64, f64) {
        let b = block.len();
        let mut mu = seed;
        let mut x = vec![C64::ONE; b];
        let normalize = |x: &mut Vec<C64>| {
            let scale = x.iter().map(|v| v.abs()).fold(0.0, f64::max);
            if scale > 0.0 && scale.is_finite() {
                for v in x.iter_mut() {
                    *v = v.scale(1.0 / scale);
                }
                true
            } else {
                false
            }
        };
        for _ in 0..8 {
            let Ok(lu) = lu_complex(&self.block_pencil(block, mu), b) else {
                // Exactly singular: μ is a root to working precision.
                return (mu, 0.0);
            };
            lu.solve(&mut x);
            if !normalize(&mut x) {
                return (mu, f64::INFINITY);
            }
            let (mut a2, mut a1, mut a0) = (C64::ZERO, C64::ZERO, C64::ZERO);
            for (r, &j) in block.iter().enumerate() {
                a2 = a2 + x[r] * x[r];
                a0 = a0 + x[r] * x[r] * C64::from_re(self.omega[j] * self.omega[j]);
                for (s, &k) in block.iter().enumerate() {
                    a1 = a1 + x[r] * C64::from_re(self.coupling.entry(j, k)) * x[s];
                }
            }
            let root = (a1 * a1 - C64::from_re(4.0) * a2 * a0).sqrt();
            let two_a2 = C64::from_re(2.0) * a2;
            let (r1, r2) = (
                (C64::ZERO - a1 + root) / two_a2,
                (C64::ZERO - a1 - root) / two_a2,
            );
            let next = if (r1 - mu).abs() <= (r2 - mu).abs() {
                r1
            } else {
                r2
            };
            let step = (next - mu).abs();
            mu = next;
            if !(step > 2.0 * f64::EPSILON * mu.abs()) {
                break;
            }
        }
        // One more inverse step aligns x with the final μ.
        if let Ok(lu) = lu_complex(&self.block_pencil(block, mu), b) {
            lu.solve(&mut x);
            if !normalize(&mut x) {
                return (mu, f64::INFINITY);
            }
        } else {
            return (mu, 0.0);
        }
        let q = self.block_pencil(block, mu);
        let (mut residual, mut x_norm, mut slope) = (0.0_f64, 0.0_f64, C64::ZERO);
        for r in 0..b {
            let mut row = C64::ZERO;
            for s in 0..b {
                row = row + q[r * b + s] * x[s];
            }
            residual += row.norm_sq();
            x_norm += x[r].norm_sq();
        }
        for (r, &j) in block.iter().enumerate() {
            slope = slope + C64::from_re(2.0) * mu * x[r] * x[r];
            for (s, &k) in block.iter().enumerate() {
                slope = slope + x[r] * C64::from_re(self.coupling.entry(j, k)) * x[s];
            }
        }
        let denominator = slope.abs();
        if denominator == 0.0 {
            return (mu, f64::INFINITY);
        }
        (mu, det::sqrt(residual) * det::sqrt(x_norm) / denominator)
    }

    fn certify_singleton(&self, i: usize) -> BlockOutcome {
        let center = self.lam[i];
        let nu = self.nu[i];
        let disc = |radius| Disc {
            center,
            radius,
            members: vec![(center, radius)],
        };
        // The exact roots of q_i lie within δ of λ̂_i and its conjugate.
        let delta = self.delta[i];
        if self.coupling.near[i].is_empty() && self.coupling.far[i] == 0.0 {
            // Row i of F vanishes: the pole is the root of q_i itself.
            return BlockOutcome::Certified(vec![disc(delta)]);
        }
        let block = [i];
        let cap = self.budget_radius(center).min(0.5 * nu);
        // |q_i| on the circle: distance to the root inside and to its mirror.
        let unperturbed = |rho: f64| (rho - delta) * (2.0 * nu - rho - delta);
        let mut rho = delta;
        for _ in 0..MAX_RADIUS_ITERATIONS {
            let far_factor = 2.0 * nu - rho - delta;
            let limit = (cap - delta) * far_factor / (1.0 + ROUCHE_SAFETY);
            let Some(p) = self.perturbation(&block, center, rho, limit) else {
                return BlockOutcome::Uncertified;
            };
            let next = delta + p * (1.0 + ROUCHE_SAFETY) / far_factor;
            if !(next <= cap) {
                return BlockOutcome::Uncertified;
            }
            if next <= rho * (1.0 + 1e-12) {
                break;
            }
            rho = next;
        }
        let rho = rho * (1.0 + 1e-6) + f64::MIN_POSITIVE;
        match self.perturbation(&block, center, rho, f64::INFINITY) {
            Some(p)
                if rho <= cap
                    && p * (1.0 + ROUCHE_SAFETY) < unperturbed(rho) * (1.0 - ROUCHE_SAFETY) =>
            {
                BlockOutcome::Certified(vec![disc(rho)])
            }
            _ => BlockOutcome::Uncertified,
        }
    }

    fn block_disc(
        &self,
        block: &[usize],
        poles: &[(C64, f64)],
        roots: &[(C64, f64)],
        members: &[usize],
    ) -> GroupDisc {
        let b = block.len();
        let center = members
            .iter()
            .fold(C64::ZERO, |acc, &l| acc + poles[l].0)
            .scale(1.0 / members.len() as f64);
        let inner: Vec<f64> = members
            .iter()
            .map(|&l| (poles[l].0 - center).abs() + poles[l].1)
            .collect();
        let outer: Vec<(usize, f64)> = roots
            .iter()
            .enumerate()
            .filter(|(l, _)| !members.contains(l))
            .map(|(l, &(r, slack))| (l, (r - center).abs() - slack))
            .collect();
        let reach = inner.iter().copied().fold(0.0_f64, f64::max);
        let cap = members
            .iter()
            .map(|&l| self.budget_radius(poles[l].0) - (poles[l].0 - center).abs())
            .fold(f64::INFINITY, f64::min)
            .min(0.5 * center.im);
        let inside = |rho: f64| inner.iter().map(|&d| rho - d).product::<f64>();
        // Everything except the inner product of Δ, as a function of ρ.
        // With `bounded`, a value whose solve would pass `cap` is cut short.
        let rest = |rho: f64, bounded: bool| -> Result<f64, Option<usize>> {
            let mut outside = 1.0_f64;
            for &(l, d) in &outer {
                if d - rho <= 0.0 {
                    return Err((l < poles.len()).then_some(l));
                }
                outside *= d - rho;
            }
            let x = det::powi(
                self.block_operator_bound(block, center, rho),
                (b - 1) as i32,
            );
            let limit = if bounded && cap > reach {
                inside(cap) * outside / (x * (1.0 + ROUCHE_SAFETY))
            } else {
                f64::INFINITY
            };
            let p = self.perturbation(block, center, rho, limit).ok_or(None)?;
            Ok(x * p * (1.0 + ROUCHE_SAFETY) / outside)
        };
        let mut rho = reach * (1.0 + 1e-12) + f64::MIN_POSITIVE;
        for _ in 0..MAX_RADIUS_ITERATIONS {
            let target = match rest(rho, true) {
                Ok(t) => t,
                Err(Some(l)) => return GroupDisc::Hits(l),
                Err(None) => return GroupDisc::Uncertified,
            };
            // Smallest ρ > reach with Π(ρ − d) = target; Π is increasing.
            let next = solve_increasing(&inside, reach, target);
            if !(next <= cap) {
                return GroupDisc::Uncertified;
            }
            if next <= rho * (1.0 + 1e-12) {
                break;
            }
            rho = next;
        }
        let rho = rho * (1.0 + 1e-6);
        match rest(rho, false) {
            Ok(t) if rho <= cap && t < inside(rho) * (1.0 - ROUCHE_SAFETY) => {
                GroupDisc::Certified(Disc {
                    center,
                    radius: rho,
                    members: members
                        .iter()
                        .map(|&l| (poles[l].0, (poles[l].0 - center).abs() + rho))
                        .collect(),
                })
            }
            Err(Some(l)) => GroupDisc::Hits(l),
            _ => GroupDisc::Uncertified,
        }
    }

    /// `X ≥ ‖Q_BB(λ)‖_F` on the disc.
    fn block_operator_bound(&self, block: &[usize], center: C64, rho: f64) -> f64 {
        let a = center.abs() + rho;
        let mut sum = 0.0_f64;
        for &j in block {
            let (near, mirror) = (
                (center - self.lam[j]).abs() + rho + self.delta[j],
                (center - self.lam[j].conj()).abs() + rho + self.delta[j],
            );
            let diag = near * mirror;
            sum = diag.mul_add(diag, sum);
            for &k in block {
                if k != j {
                    let v = a * self.coupling.entry(j, k);
                    sum = v.mul_add(v, sum);
                }
            }
        }
        det::sqrt(sum)
    }

    /// `P(ρ)`: a bound on `‖λ² F_BK Q_KK⁻¹ F_KB‖₂` over the disc, or `None`
    /// when a complement pole enters the disc or `a·e ≥ 1`.
    ///
    /// `limit` short-circuits: `P ≥ a²S_B`, so once `a²S_B > limit` the
    /// caller's acceptance test cannot pass and the Perron bound is skipped.
    fn perturbation(&self, block: &[usize], center: C64, rho: f64, limit: f64) -> Option<f64> {
        let n = self.n();
        let in_block = |k: usize| block.binary_search(&k).is_ok();
        let a = center.abs() + rho;
        let y = center.im;
        let mut m = vec![0.0_f64; n];
        let mut m_all = f64::INFINITY;
        let mut m_far = f64::INFINITY;
        let half_window = self.coupling.window.map(|w| 0.5 * w);
        let is_near = |k: usize| half_window.is_none_or(|h| (self.omega[k] - y).abs() < h);
        for k in (0..n).filter(|&k| !in_block(k)) {
            let f1 = (center - self.lam[k]).abs() - rho - self.delta[k];
            let f2 = (center - self.lam[k].conj()).abs() - rho - self.delta[k];
            if !(f1 > 0.0 && f2 > 0.0) {
                return None;
            }
            m[k] = f1 * f2;
            m_all = m_all.min(m[k]);
            if !is_near(k) {
                m_far = m_far.min(m[k]);
            }
        }
        // S_B: block rows against the complement.
        let mut s = 0.0_f64;
        for &j in block {
            for &(k, f) in &self.coupling.near[j] {
                if !in_block(k) {
                    s += f * f / m[k];
                }
            }
            let far = self.coupling.far[j];
            if far > 0.0 {
                s += far * far / m_all;
            }
        }
        if a * a * s > limit {
            return None;
        }
        let near_set: Vec<usize> = (0..n).filter(|&k| !in_block(k) && is_near(k)).collect();
        let e_nn = self.perron_bound(&near_set, &m, &in_block, &is_near);
        let (mut e_ff, mut e_nf_sq) = (0.0_f64, 0.0_f64);
        if m_far.is_finite() {
            for k in (0..n).filter(|&k| !in_block(k) && !is_near(k)) {
                e_ff = e_ff.max(self.coupling.l1[k] / det::sqrt(m[k] * m_far));
            }
            for &j in &near_set {
                let mut row = 0.0_f64;
                for &(k, f) in &self.coupling.near[j] {
                    if !in_block(k) && !is_near(k) {
                        row += f * f / m[k];
                    }
                }
                let far = self.coupling.far[j];
                row += far * far / m_far;
                e_nf_sq += row / m[j];
            }
        }
        let e = e_nn.max(e_ff) + det::sqrt(e_nf_sq);
        let ae = a * e * (1.0 + ROUCHE_SAFETY);
        if !(ae < 1.0) {
            return None;
        }
        Some(a * a * s / (1.0 - ae))
    }

    /// Collatz–Wielandt upper bound on the Perron root of `Ē` restricted to
    /// `set`, using the known entries.
    fn perron_bound(
        &self,
        set: &[usize],
        m: &[f64],
        in_block: &dyn Fn(usize) -> bool,
        is_near: &dyn Fn(usize) -> bool,
    ) -> f64 {
        if set.is_empty() {
            return 0.0;
        }
        let n = self.n();
        let scale: Vec<f64> = m
            .iter()
            .map(|&v| if v > 0.0 { 1.0 / det::sqrt(v) } else { 0.0 })
            .collect();
        let mut x = vec![0.0_f64; n];
        for &k in set {
            x[k] = 1.0;
        }
        let mut best = f64::INFINITY;
        for _ in 0..PERRON_STEPS {
            let mut y = vec![0.0_f64; n];
            let mut quotient = 0.0_f64;
            let mut largest = 0.0_f64;
            for &j in set {
                let mut acc = 0.0_f64;
                for &(k, f) in &self.coupling.near[j] {
                    if !in_block(k) && is_near(k) {
                        acc += f.abs() * scale[k] * x[k];
                    }
                }
                y[j] = acc * scale[j];
                quotient = quotient.max(y[j] / x[j]);
                largest = largest.max(y[j]);
            }
            best = best.min(quotient);
            if largest == 0.0 {
                return 0.0;
            }
            // Keep x positive on the set; any positive x gives a valid bound.
            for &j in set {
                x[j] = y[j] + 1e-3 * largest;
            }
        }
        best
    }

    fn assemble(
        &self,
        blocks: &[Vec<usize>],
        outcomes: Vec<BlockOutcome>,
    ) -> CertifiedModalDamping {
        let mut result = Vec::with_capacity(blocks.len());
        let (mut worst_rel_decay, mut worst_freq_cents) = (0.0_f64, 0.0_f64);
        for (block, outcome) in blocks.iter().zip(outcomes) {
            let BlockOutcome::Certified(discs) = outcome else {
                unreachable!("assembled only when every block is certified");
            };
            let mut members: Vec<(C64, f64)> = discs.into_iter().flat_map(|d| d.members).collect();
            members.sort_by(|x, y| x.0.im.total_cmp(&y.0.im).then(x.0.re.total_cmp(&y.0.re)));
            for &(mu, r) in &members {
                if r > 0.0 {
                    worst_rel_decay = worst_rel_decay.max(r / mu.re.abs());
                    worst_freq_cents = worst_freq_cents.max(cents_bound(r, mu.im));
                }
            }
            if let [i] = **block {
                let (eigenvalue, radius) = members[0];
                result.push(ModalDampingBlock::Diagonal {
                    mode: i,
                    zeta: self.coupling.diag[i] / (2.0 * self.omega[i]),
                    eigenvalue,
                    radius,
                });
            } else {
                let radius = members.iter().map(|m| m.1).fold(0.0, f64::max);
                result.push(ModalDampingBlock::Full {
                    modes: block.clone(),
                    eigenvalues: members.iter().map(|m| m.0).collect(),
                    radii: members.iter().map(|m| m.1).collect(),
                    radius,
                });
            }
        }
        CertifiedModalDamping {
            blocks: result,
            worst_rel_decay,
            worst_freq_cents,
            far_coupling: FarCouplingBound {
                window_rad_s: self.coupling.window,
                far_row_norm: self.coupling.far.clone(),
            },
        }
    }
}

enum GroupDisc {
    Certified(Disc),
    Uncertified,
    /// The disc would reach this other root (an index into the upper poles,
    /// or a conjugate when out of range).
    Hits(usize),
}

/// Upper bound on the frequency error in cents of a pole within `r` of
/// one with imaginary part `im`: `−1200 log2(1 − r/im) ≤ 1200/ln2 · x/(1 − x)`.
fn cents_bound(r: f64, im: f64) -> f64 {
    let x = r / im;
    if x >= 1.0 {
        return f64::INFINITY;
    }
    1200.0 / core::f64::consts::LN_2 * x / (1.0 - x)
}

/// Smallest `ρ > reach` with `f(ρ) = target` for increasing `f`, `f(reach) = 0`.
fn solve_increasing(f: &dyn Fn(f64) -> f64, reach: f64, target: f64) -> f64 {
    let mut lo = reach;
    let mut step = reach.max(f64::MIN_POSITIVE);
    let mut hi = reach + step;
    while f(hi) < target {
        step *= 2.0;
        hi = reach + step;
        if !hi.is_finite() {
            return f64::INFINITY;
        }
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if !(lo < mid && mid < hi) {
            break;
        }
        if f(mid) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    hi
}

fn find(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

/// Current blocks, each ascending, ordered by first mode.
fn groups(parent: &mut [usize]) -> Vec<Vec<usize>> {
    let n = parent.len();
    let mut by_root: Vec<Vec<usize>> = vec![Vec::new(); n];
    for x in 0..n {
        let r = find(parent, x);
        by_root[r].push(x);
    }
    let mut blocks: Vec<Vec<usize>> = by_root.into_iter().filter(|b| !b.is_empty()).collect();
    blocks.sort_by_key(|b| b[0]);
    blocks
}

fn merge_groups(groups: &mut Vec<Vec<usize>>, a: usize, b: usize) {
    let (low, high) = (a.min(b), a.max(b));
    let moved = groups.remove(high);
    groups[low].extend(moved);
    groups[low].sort_unstable();
}
