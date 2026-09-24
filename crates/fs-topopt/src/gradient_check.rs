//! Budgeted finite-difference audits of the complete multi-load design chain.
//!
//! Real filter, projection and elasticity evaluations are used at every probe.
//! Two spatially varying, bound-aware directions collectively touch every density;
//! second-order differences are compared with BOTH the compliance pullback and
//! the physical-volume pullback. These are numerical directional checks, not
//! componentwise proofs, error enclosures, or a stationarity certificate.

use crate::control::{EvaluationStop, SolveControl, SolveWork};
use crate::elasticity::DensityElasticity;
use crate::pipeline::{DesignPipeline, LoadCase, SimpParams, assert_valid_load_cases};

/// Explicit finite-difference experiment; no implicit parameter rescaling.
#[derive(Debug, Clone, Copy)]
pub struct GradientCheckOptions {
    /// Literal raw-density step. Every active probe uses rho+h*d and rho+2h*d.
    pub step: f64,
    /// Maximum relative discrepancy for each objective and volume direction.
    pub relative_tolerance: f64,
}

impl Default for GradientCheckOptions {
    fn default() -> Self {
        Self { step: 1e-4, relative_tolerance: 5e-4 }
    }
}

impl GradientCheckOptions {
    pub(crate) fn assert_valid(self) {
        assert!(self.step.is_finite() && self.step > 0.0 && self.step <= 0.125,
            "gradient-check step must be finite and lie in (0, 0.125]");
        assert!(self.relative_tolerance.is_finite() && self.relative_tolerance > 0.0
            && self.relative_tolerance < 1.0,
            "gradient-check relative tolerance must be finite and lie in (0, 1)");
    }
}

/// Deterministic direction reconstructed from the baseline densities and step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GradientDirection {
    /// Positive alternating weights (1 on even indices, 1/2 on odd indices),
    /// masked to zero where rho+2h*d would exceed one.
    Increasing,
    /// Negative complementary weights (-1/2 on even indices, -1 on odd
    /// indices), masked to zero where rho+2h*d would fall below zero.
    Decreasing,
}

/// One complete second-order one-sided directional experiment.
#[derive(Debug, Clone)]
pub struct GradientProbe {
    /// Which bound-aware mask was used.
    pub direction: GradientDirection,
    /// Number of nonzero direction entries; never zero in a retained probe.
    pub active_densities: usize,
    /// Weighted independent-load gradient dotted with the probe direction.
    pub compliance_analytic: f64,
    /// The independently re-solved compliance finite difference.
    pub compliance_difference: f64,
    /// Scale-free discrepancy, with no fixed absolute floor hiding tiny loads.
    pub compliance_relative_error: f64,
    /// Physical-volume gradient dotted with the same direction.
    pub volume_analytic: f64,
    /// The independently evaluated physical-volume finite difference.
    pub volume_difference: f64,
    /// Scale-free physical-volume discrepancy.
    pub volume_relative_error: f64,
}

/// Completed numerical evidence at ONE baseline and ONE model configuration.
#[derive(Debug, Clone)]
pub struct MultiLoadGradientCheck {
    /// Parameters of the tested design map.
    pub params: SimpParams,
    /// Exact experiment settings.
    pub options: GradientCheckOptions,
    /// Baseline weighted compliance; not a perturbed or partially solved value.
    pub baseline_compliance: f64,
    /// Baseline physical material volume fraction.
    pub baseline_volume_fraction: f64,
    /// Nonempty directional experiments. Empty masks are skipped, never counted
    /// as successful checks. Their union touches every design coordinate.
    pub probes: Vec<GradientProbe>,
    /// Cumulative control work, including earlier work using the same control.
    pub work: SolveWork,
}

impl MultiLoadGradientCheck {
    /// Whether every retained numerical discrepancy meets the requested gate.
    /// This does not promote the result to a mathematical certificate.
    #[must_use]
    pub fn passed(&self) -> bool {
        !self.probes.is_empty() && self.probes.iter().all(|probe| {
            probe.active_densities > 0
                && probe.compliance_relative_error.is_finite()
                && probe.volume_relative_error.is_finite()
                && probe.compliance_relative_error <= self.options.relative_tolerance
                && probe.volume_relative_error <= self.options.relative_tolerance
        })
    }
}

fn breakdown(stage: &'static str) -> EvaluationStop {
    EvaluationStop::Breakdown { stage }
}

fn difference(base: f64, one: f64, two: f64, step: f64) -> Result<f64, EvaluationStop> {
    // Difference first: forming -3*base+4*one-two would unnecessarily overflow
    // at large but otherwise representable compliance scales.
    let value = (2.0 * (one - base) - 0.5 * (two - base)) / step;
    if value.is_finite() { Ok(value) }
    else { Err(breakdown("gradient-check-difference")) }
}

fn relative_error(analytic: f64, measured: f64) -> Result<f64, EvaluationStop> {
    if !analytic.is_finite() || !measured.is_finite() {
        return Err(breakdown("gradient-check-comparison"));
    }
    let scale = analytic.abs().max(measured.abs());
    if scale == 0.0 { return Ok(0.0); }
    // Divide before subtraction: opposite near-MAX values still produce the
    // finite discrepancy 2, not an infinity that a caller might discard.
    Ok((analytic / scale - measured / scale).abs())
}

fn volume(
    pipeline: &DesignPipeline, rho: &[f64], volumes: &[f64], total: f64,
    control: &mut SolveControl<'_>,
) -> Result<f64, EvaluationStop> {
    let (_, projected, _) = pipeline.try_forward(rho, control)?;
    let value: f64 = projected.iter().zip(volumes).map(|(r, v)| r * (v / total)).sum();
    if value.is_finite() { Ok(value) }
    else { Err(breakdown("gradient-check-volume")) }
}

/// Shared numerical experiment for the tetrahedral and cut-cell design maps.
/// The owner supplies actual equilibria and retains its own rollback boundary.
pub(crate) struct GradientSample {
    pub compliance: f64,
    pub compliance_gradient: Vec<f64>,
    pub volume: f64,
    pub volume_gradient: Vec<f64>,
}

pub(crate) fn audit_design_map(
    params: SimpParams,
    rho: &[f64],
    options: GradientCheckOptions,
    control: &mut SolveControl<'_>,
    mut evaluate: impl FnMut(&[f64], bool, &mut SolveControl<'_>) -> Result<GradientSample, EvaluationStop>,
) -> Result<MultiLoadGradientCheck, EvaluationStop> {
    control.checkpoint("gradient-check")?;
    let baseline = evaluate(rho, true, control)?;
    let mut probes = Vec::with_capacity(2);
    for direction in [GradientDirection::Increasing, GradientDirection::Decreasing] {
        control.checkpoint("gradient-probe")?;
        let h = options.step;
        // Complementary heterogeneous weights exercise spatial pullbacks;
        // uniform material scalings could conceal a transpose error.
        let mask: Vec<f64> = rho.iter().enumerate().map(|(index, &r)| {
            let even = index % 2 == 0;
            let d = match direction {
                GradientDirection::Increasing => if even { 1.0 } else { 0.5 },
                GradientDirection::Decreasing => if even { -0.5 } else { -1.0 },
            };
            if (0.0..=1.0).contains(&(r + (2.0 * h) * d)) { d } else { 0.0 }
        }).collect();
        let active = mask.iter().filter(|&&d| d != 0.0).count();
        if active == 0 { continue; }
        let one: Vec<f64> = rho.iter().zip(&mask).map(|(r, d)| r + h * d).collect();
        let two: Vec<f64> = rho.iter().zip(&mask).map(|(r, d)| r + (2.0 * h) * d).collect();
        for (((&r, &d), &a), &b) in rho.iter().zip(&mask).zip(&one).zip(&two) {
            if !a.is_finite() || !b.is_finite() || !(0.0..=1.0).contains(&a)
                || !(0.0..=1.0).contains(&b) || (d != 0.0 && (a == r || a == b)) {
                return Err(breakdown("gradient-check-stencil"));
            }
        }
        let ca: f64 = baseline.compliance_gradient.iter().zip(&mask).map(|(g, d)| g * d).sum();
        let va: f64 = baseline.volume_gradient.iter().zip(&mask).map(|(g, d)| g * d).sum();
        let first = evaluate(&one, false, control)?;
        let second = evaluate(&two, false, control)?;
        let cd = difference(baseline.compliance, first.compliance, second.compliance, h)?;
        let vd = difference(baseline.volume, first.volume, second.volume, h)?;
        probes.push(GradientProbe {
            direction, active_densities: active,
            compliance_analytic: ca, compliance_difference: cd,
            compliance_relative_error: relative_error(ca, cd)?,
            volume_analytic: va, volume_difference: vd,
            volume_relative_error: relative_error(va, vd)?,
        });
    }
    if probes.is_empty() { return Err(breakdown("gradient-check-empty")); }
    control.checkpoint("gradient-check-publish")?;
    Ok(MultiLoadGradientCheck {
        params, options, baseline_compliance: baseline.compliance,
        baseline_volume_fraction: baseline.volume, probes, work: control.work(),
    })
}

/// Audit the full design-to-physics chain without changing the caller's operator.
///
/// Uses at most five independent-load evaluations: baseline plus two second-
/// order one-sided probes. Their complementary alternating weights are not
/// uniform material scalings. Probe directions do not leave the density box.
/// Each probe uses (-3J(rho)+4J(rho+h*d)-J(rho+2h*d))/(2h), not differentiation
/// through Krylov iterations. All real solves share `control`, including work
/// discarded on failure. Cancellation or unrepresentable arithmetic returns no
/// partial audit. Moduli are restored on BOTH success and any returned stop.
/// Invalid model inputs panic before operator mutation.
///
/// A passed audit tests two directions at this baseline only. It does not prove
/// every gradient component, later optimization iterates, or continuum accuracy.
#[allow(clippy::too_many_arguments)]
pub fn controlled_multi_load_gradient_check(
    pipeline: &DesignPipeline,
    elasticity: &mut DensityElasticity,
    loads: &[LoadCase<'_>],
    rho: &[f64],
    cell_vol: &[f64],
    options: GradientCheckOptions,
    control: &mut SolveControl<'_>,
) -> Result<MultiLoadGradientCheck, EvaluationStop> {
    options.assert_valid();
    pipeline.params.assert_valid();
    assert_valid_load_cases(elasticity, loads);
    assert_eq!(rho.len(), elasticity.cells(), "one density per elasticity cell is required");
    assert_eq!(rho.len(), cell_vol.len(), "one volume per design cell is required");
    assert!(!rho.is_empty() && rho.iter().all(|r| r.is_finite() && (0.0..=1.0).contains(r)),
        "gradient-check densities must be nonempty, finite, and lie in [0, 1]");
    assert!(cell_vol.iter().all(|v| v.is_finite() && *v > 0.0),
        "gradient-check cell volumes must be finite and positive");
    let total: f64 = cell_vol.iter().sum();
    assert!(total.is_finite(), "total cell volume must be finite");
    let previous_moduli = elasticity.moduli.clone();
    let outcome = audit_design_map(pipeline.params, rho, options, control, |point, baseline, control| {
        let before = if baseline {
            Some(pipeline.try_volume_and_gradient(point, cell_vol, control)?)
        } else { None };
        let solved = pipeline.try_multi_load_compliance_and_gradient(elasticity, point, loads, control)?;
        let (volume, volume_gradient) = match before {
            Some(values) => values,
            None => (volume(pipeline, point, cell_vol, total, control)?, Vec::new()),
        };
        Ok(GradientSample {
            compliance: solved.compliance, compliance_gradient: solved.gradient,
            volume, volume_gradient,
        })
    });
    elasticity.moduli = previous_moduli;
    outcome
}

#[cfg(test)]
mod tests {
    use super::{difference, relative_error};

    #[test]
    fn g0_second_order_quadratic_difference_respects_extreme_load_scales() {
        for scale in [1e-200, 1.0, 1e200] {
            let f = |x: f64| scale * (x * x + 2.0 * x + 1.0);
            let (x, h) = (0.25, 1e-3);
            let derivative = difference(f(x), f(x + h), f(x + 2.0 * h), h).unwrap();
            assert!(relative_error(2.5 * scale, derivative).unwrap() < 1e-10);
        }
    }

    #[test]
    fn g0_wrong_or_nonfinite_gradients_cannot_pass_via_an_absolute_floor() {
        for scale in [1e-200, 1.0, 1e200] {
            assert_eq!(relative_error(0.0, scale).unwrap(), 1.0);
            assert_eq!(relative_error(-scale, scale).unwrap(), 2.0);
        }
        assert_eq!(relative_error(-f64::MAX, f64::MAX).unwrap(), 2.0);
        assert!(relative_error(f64::NAN, 1.0).is_err());
        assert!(difference(0.0, f64::MAX, -f64::MAX, 1.0).is_err());
    }
}
