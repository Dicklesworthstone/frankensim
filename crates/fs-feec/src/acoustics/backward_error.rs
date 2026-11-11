//! The backward error of one row of a refined solve: componentwise, with the rounding floor of the returned iterate.
//!
//! The solvers that use this return an `f64` solution after at most three residual evaluations in double-double and accept it when,
//! row by row, `|r_i| <= limit * d_i`, where `r_i = b_i - sum_j a_ij x_j` and `d_i = |b_i| + sum_j |a_ij| |x_j|` is the sum of the
//! row's absolute terms. That componentwise bound certifies a row only if the solution is accurate componentwise in that row, and
//! the returned `x` is not. Refinement converges in norm: every component carries the rounding of its last `f64` addition, at most
//! `eps |x_j|`, and a component whose exact value is 0 carries the dust that the contraction of the last corrections has not yet
//! removed. What such an `x` can reach on row `i` is therefore not `eps d_i` but
//!
//! ```text
//! |r_i| <= |sum_j a_ij dx_j| <= eps * N_i,        N_i = (sum_j |a_ij|) * max_j |x_j| + |b_i|,
//! ```
//!
//! the row's normwise scale `N_i` (the refinement is assumed to have contracted its error below that rounding, which is what it is
//! for: a solve that has not, with `kappa * eps` near 1, leaves a residual far above `eps N_i` and fails). For a row whose data are
//! exactly 0 and whose exact solution is 0 on its columns, `d_i` is itself dust and `|r_i| / d_i` divides roundoff by roundoff. It
//! was 0.294, 0.0419 and 2.9e-6 for `lossless_wall_pole_retains_motion_storage_and_is_crossed_by_the_full_survey`,
//! `pressure_release_wall_spring_balance_and_accelerating_fluid_have_exact_work` and
//! `g1_original_memory_has_correct_reactive_impedance_storage_and_heat` (rows with `d_i` of 1.6e-33, 7.4e-22 and 2.0e-48 beside
//! largest row scales of 0.447, 1.1e-5 and 2.0), while `|r_i|` divided by the largest row scale of its solve was 1e-33, 3e-18 and 3e-54.
//!
//! Arioli, Demmel and Duff, "Solving sparse linear systems with sparse backward error", SIAM J. Matrix Anal. Appl. 10 (1989), bound
//! such rows normwise and the others componentwise (their mixed criterion). The same is one expression here: the ratio is
//! `|r_i| / (d_i + (eps / limit) N_i)`, so the gate `ratio <= limit` is `|r_i| <= limit d_i + eps N_i`, the componentwise tolerance
//! plus the rounding floor. A row with substantial terms (`d_i` far above `(eps / limit) N_i`, i.e. above 4e-4 N_i at the limit
//! 5e-13) is unchanged to the first order; a wrong value in a zero row is still refused as soon as it exceeds `eps N_i`.
//!
//! The mixed ratio decides ACCEPTANCE only. A refinement's early exit ("exact to double-double", `< 1e-18`) keeps the strict
//! componentwise ratio: the mixed one is smaller, and stopping on it ends the refinement at the first iterate (the power balance of
//! `curved_moving_wall_retains_full_vector_no_slip_and_reciprocal_work` was 1.2e-7 against its 2e-8 when it did).

/// The denominator of the ratio `|r_i| / denominator` that a solver compares with `limit`: the row's absolute `terms` plus the
/// normwise floor `(eps / limit) * norm`, `norm` being the row's normwise scale `N_i`.
#[must_use]
pub(crate) fn row_scale(terms: f64, norm: f64, limit: f64) -> f64 {
    terms + norm * (f64::EPSILON / limit)
}

#[cfg(test)]
mod tests;
