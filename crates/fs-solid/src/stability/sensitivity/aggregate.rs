//! Smooth total buckling gradients for complete packets of lowest modes.
//!
//! The conservative KS value is `min(lambda) - log(sum(exp(-rho *
//! (lambda_i - min(lambda)))))/rho`. Its positive weights sum to one.
//! Complete repeated eigenspaces are admitted: equal eigenvalue weights make
//! the aggregate derivative a basis-independent trace. Cutting a cluster at
//! the retained packet boundary is refused. No individual repeated-mode
//! derivatives are published. If sharpness would amplify native residual or
//! roundoff resolution inside an active cluster into a log-weight uncertainty
//! above 1e-6, analysis refuses. This is a numerical admission policy, not a
//! certified eigenvalue forward-error bound.
//!
//! Only one static adjoint is solved, regardless of the retained mode count.
//! The right-hand side sums each mode's displacement pullback with its own
//! physical load, rather than the aggregate load. Element contractions reuse
//! the original P1/Q1, Standard/B-bar quadrature and its exact transpose.

use super::{
    BucklingResult, CheckedCgConfig, CheckedCgError, CsrOp, Jacobi, LinearProblem, SolidError,
    buckling_loads_cancellable, check_cancelled, checked_cg, contractions, expand_mode, invalid,
    reduced_pencil_scaled,
};

/// An accepted complete spectral packet, its KS policy and its total gradient.
///
/// The physical modes, element design and smoothing policy cannot be replaced
/// independently after analysis. This owner exposes only an aggregate gradient;
/// individual branch gradients need not exist at repeated eigenvalues.
#[derive(Debug, Clone)]
pub struct ElementScaleAggregateSensitivity {
    buckling: BucklingResult,
    scales: Vec<f64>,
    value: f64,
    gradient: Vec<f64>,
    weights: Vec<f64>,
    sharpness: f64,
    adjoint_relative_residual: f64,
    adjoint_iterations: usize,
}

impl ElementScaleAggregateSensitivity {
    /// Exactly the requested lowest positive modes and their pencil residuals.
    #[must_use]
    pub fn buckling(&self) -> &BucklingResult {
        &self.buckling
    }

    /// Original dimensionless element Young's-modulus multipliers.
    #[must_use]
    pub fn element_scales(&self) -> &[f64] {
        &self.scales
    }

    /// Conservative smooth lower bound on the packet's smallest physical load.
    ///
    /// Finite negative values are valid when the smoothing bias exceeds the
    /// smallest load; a caller must not interpret them as physical eigenvalues.
    #[must_use]
    pub fn value(&self) -> f64 {
        self.value
    }

    /// Total derivative of the KS value with respect to each element multiplier.
    #[must_use]
    pub fn gradient(&self) -> &[f64] {
        &self.gradient
    }

    /// `d KS / d lambda_i`, in the same order as the retained physical modes.
    #[must_use]
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    /// Fixed packet size used by this analysis.
    #[must_use]
    pub fn mode_count(&self) -> usize {
        self.weights.len()
    }

    /// Positive smoothing sharpness, in inverse physical load-factor units.
    #[must_use]
    pub fn sharpness(&self) -> f64 {
        self.sharpness
    }

    /// Recomputed true relative residual of the single scaled static adjoint.
    #[must_use]
    pub fn adjoint_relative_residual(&self) -> f64 {
        self.adjoint_relative_residual
    }

    /// Total native CG iterations, including any residual defect corrections.
    #[must_use]
    pub fn adjoint_iterations(&self) -> usize {
        self.adjoint_iterations
    }
}

/// Compute a smooth total sizing derivative through a complete buckling packet.
///
/// Exactly `modes` lowest positive modes must be available. One additional pair
/// checks the packet boundary; a relative boundary gap at most `1e-6` refuses
/// because an omitted member of a repeated cluster destroys trace invariance.
/// Repeated or crossing modes strictly inside the packet are allowed. The
/// original native buckling solver's factor-size and iteration limits apply.
///
/// Geometry, Poisson ratio, dead loads and prescribed values are design
/// independent. The derivative includes the re-equilibrated stress field and
/// the full prescribed-column lift. Apply any density/material parameterization
/// chain rule outside this owner. This is infinitesimal initial-stress buckling,
/// not a derivative of a nonlinear or load-following equilibrium.
///
/// # Errors
/// Invalid scales, nonpositive/nonfinite sharpness, an incomplete packet,
/// unresolved boundary separation, excessive cluster sharpness, unrepresentable results or failed native
/// equilibrium, eigenpair and adjoint residual gates return `SolidError`.
pub fn element_scale_aggregate_sensitivities(
    problem: &LinearProblem<'_>,
    scales: &[f64],
    modes: usize,
    iterations: usize,
    sharpness: f64,
) -> Result<ElementScaleAggregateSensitivity, SolidError> {
    element_scale_aggregate_sensitivities_cancellable(
        problem,
        scales,
        modes,
        iterations,
        sharpness,
        &|| false,
    )
}

/// The same complete-packet analysis with the original cancellation callback.
///
/// Cancellation is checked during assembly, eigenvalue and adjoint iterations,
/// and element contractions. No partial aggregate result is published.
///
/// # Errors
/// As [`element_scale_aggregate_sensitivities`], plus `SolidError::Cancelled`.
pub fn element_scale_aggregate_sensitivities_cancellable(
    problem: &LinearProblem<'_>,
    scales: &[f64],
    modes: usize,
    iterations: usize,
    sharpness: f64,
    cancelled: &impl Fn() -> bool,
) -> Result<ElementScaleAggregateSensitivity, SolidError> {
    check_cancelled(cancelled)?;
    if modes == 0 || iterations == 0 || !sharpness.is_finite() || sharpness <= 0.0 {
        return Err(invalid(
            "aggregate buckling sensitivity requires positive mode and iteration allowances and finite positive sharpness",
        ));
    }
    let (k, kg, map, reference) = reduced_pencil_scaled(problem, Some(scales), cancelled)?;
    let mut buckling = buckling_loads_cancellable(
        &k,
        &kg,
        &map,
        modes.saturating_add(1),
        iterations,
        cancelled,
    )?;
    if buckling.loads.len() < modes {
        return Err(invalid(
            "aggregate buckling sensitivity requires the complete requested positive mode packet",
        ));
    }
    if let Some(&next) = buckling.loads.get(modes) {
        if (next - buckling.loads[modes - 1]) / next <= 1e-6 {
            return Err(invalid(
                "aggregate buckling sensitivity packet boundary splits an unresolved eigenvalue cluster",
            ));
        }
    }
    buckling.loads.truncate(modes);
    buckling.modes.truncate(modes);
    buckling.relative_residuals.truncate(modes);
    let (value, weights) = ks_value_and_weights(&buckling.loads, sharpness)?;
    admit_cluster_weight_resolution(&buckling, &weights, sharpness)?;

    // Native modes have unit geometric work, but use their measured work in
    // every contraction rather than silently assuming perfect normalization.
    let mut denominators = Vec::with_capacity(modes);
    let mut coefficients = Vec::with_capacity(modes);
    for (mode, (&load, &weight)) in buckling
        .modes
        .iter()
        .zip(buckling.loads.iter().zip(&weights))
    {
        check_cancelled(cancelled)?;
        let mut geometric_mode = vec![0.0; mode.len()];
        kg.spmv(mode, &mut geometric_mode);
        let denominator = -mode
            .iter()
            .zip(&geometric_mode)
            .map(|(a, b)| a * b)
            .sum::<f64>();
        let coefficient = (weight * load) / denominator;
        if !denominator.is_finite() || denominator <= 0.0 || !coefficient.is_finite() {
            return Err(invalid(
                "aggregate buckling sensitivity requires finite positive geometric modal work and finite adjoint weights",
            ));
        }
        denominators.push(denominator);
        coefficients.push(coefficient);
    }
    let adjoint_scale = coefficients.iter().copied().fold(0.0_f64, f64::max);
    if !adjoint_scale.is_finite() || adjoint_scale <= 0.0 {
        return Err(invalid(
            "aggregate buckling adjoint scaling is not representable",
        ));
    }

    // q = sum_i w_i lambda_i q_i / D_i. Solve for z / adjoint_scale
    // instead: lambda_i*q_i can overflow unnecessarily under tiny reference
    // loads, even when the physical total derivative is representable.
    let mut rhs = vec![0.0; k.nrows()];
    let mut direct = vec![0.0; scales.len()];
    let mut reference_forces = None;
    for i in 0..modes {
        check_cancelled(cancelled)?;
        if weights[i] == 0.0 {
            continue;
        }
        let full_mode = expand_mode(&buckling.modes[i], &map);
        let (mode_rhs, mode_direct, forces) = contractions(
            problem,
            scales,
            &map,
            &reference,
            &full_mode,
            buckling.loads[i],
            cancelled,
        )?;
        let rhs_weight = coefficients[i] / adjoint_scale;
        let direct_weight = weights[i] / denominators[i];
        for (sum, contribution) in rhs.iter_mut().zip(mode_rhs) {
            *sum = rhs_weight.mul_add(contribution, *sum);
        }
        for (sum, contribution) in direct.iter_mut().zip(mode_direct) {
            *sum = direct_weight.mul_add(contribution, *sum);
        }
        if reference_forces.is_none() {
            // K'_e u is mode independent; retain only one element-force array.
            reference_forces = Some(forces);
        }
    }
    if rhs.iter().chain(&direct).any(|value| !value.is_finite()) {
        return Err(invalid(
            "aggregate buckling contractions are not representable",
        ));
    }
    let reference_forces = reference_forces.ok_or_else(|| {
        invalid("aggregate buckling sensitivity has no representable spectral weight")
    })?;
    let preconditioner = Jacobi::new(&k);
    let op = CsrOp::symmetric(k);
    let adjoint = checked_cg(
        &op,
        &preconditioner,
        &rhs,
        CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 40_000,
            max_corrections: 2,
        },
        |_| check_cancelled(cancelled),
    )
    .map_err(|error| match error {
        CheckedCgError::InvalidInput(what) => invalid(what),
        CheckedCgError::Interrupted(error) => error,
    })?;
    let residual = adjoint.report.rel_residual;
    if !residual.is_finite() || residual > 1e-8 {
        return Err(SolidError::SolveFailed {
            iters: adjoint.report.iters,
            rel_residual: residual,
        });
    }
    let mut gradient = Vec::with_capacity(scales.len());
    for (element, conn) in problem.mesh.elems.iter().enumerate() {
        check_cancelled(cancelled)?;
        let mut redistribution = 0.0;
        for (a, &node) in conn.iter().enumerate() {
            for c in 0..2 {
                if let Some(free) = map[2 * node + c] {
                    redistribution = adjoint.x[free]
                        .mul_add(reference_forces[element][2 * a + c], redistribution);
                }
            }
        }
        let derivative = (-adjoint_scale).mul_add(redistribution, direct[element]);
        if !derivative.is_finite() {
            return Err(invalid(
                "aggregate total buckling derivative is not representable",
            ));
        }
        gradient.push(derivative);
    }
    check_cancelled(cancelled)?;
    Ok(ElementScaleAggregateSensitivity {
        buckling,
        scales: scales.to_vec(),
        value,
        gradient,
        weights,
        sharpness,
        adjoint_relative_residual: residual,
        adjoint_iterations: adjoint.report.iters,
    })
}

/// An almost repeated spectrum cannot support arbitrary discrimination between
/// its branches. Native pencil residuals are backward-error measurements, not
/// eigenvalue error bounds; combine them with a 64-epsilon roundoff floor as a
/// conservative *admission indicator*. The exponent difference controls each
/// softmax weight ratio, so limit its resolution amplification to 1e-6.
///
/// Only interior unresolved clusters with a nonzero aggregate weight need this
/// gate. Separated modes and fully underflowed higher clusters may still use a
/// sharpness that saturates the conservative minimum. Do not average physical
/// eigenvalues or silently replace the requested smoothing policy.
fn admit_cluster_weight_resolution(
    buckling: &BucklingResult,
    weights: &[f64],
    sharpness: f64,
) -> Result<(), SolidError> {
    for i in 0..buckling.loads.len().saturating_sub(1) {
        let left = buckling.loads[i];
        let right = buckling.loads[i + 1];
        if (right - left) / right > 1e-6 || (weights[i] == 0.0 && weights[i + 1] == 0.0) {
            continue;
        }
        let roundoff = 64.0 * f64::EPSILON;
        let resolution = left * buckling.relative_residuals[i].max(roundoff)
            + right * buckling.relative_residuals[i + 1].max(roundoff);
        if sharpness * resolution > 1e-6 {
            return Err(invalid(
                "aggregate buckling sharpness exceeds the native residual/roundoff resolution of an active eigenvalue cluster",
            ));
        }
    }
    Ok(())
}

/// Shift before exponentiation: at least one term is exactly one, while large
/// gaps may safely underflow to zero. An overflowing negative exponent is also
/// benign; overflow in the resulting finite-load aggregate is not.
fn ks_value_and_weights(loads: &[f64], sharpness: f64) -> Result<(f64, Vec<f64>), SolidError> {
    let minimum = loads.iter().copied().fold(f64::INFINITY, f64::min);
    let mut weights: Vec<f64> = loads
        .iter()
        .map(|&load| (-sharpness * (load - minimum)).exp())
        .collect();
    let sum: f64 = weights.iter().sum();
    let value = minimum - sum.ln() / sharpness;
    if !value.is_finite() || !sum.is_finite() || sum <= 0.0 {
        return Err(invalid("aggregate buckling KS value is not representable"));
    }
    for weight in &mut weights {
        *weight /= sum;
    }
    if weights
        .iter()
        .any(|weight| !weight.is_finite() || *weight < 0.0)
    {
        return Err(invalid(
            "aggregate buckling KS weights are not representable",
        ));
    }
    Ok((value, weights))
}
