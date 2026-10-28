//! Minimum dilated material under one authored cap for ALL f=Fq, ||q||<=1.
//! This is not a list of independent per-column compliance requirements.
pub mod start;
pub mod binary;
pub mod refinement;
pub mod refit;
pub mod campaign;
use super::*;
pub use crate::sdf3::robust::envelope::limits::ComplianceLimitedOptions3 as LoadSpaceLimitedOptions3;

/// A physical upper compliance requirement for an entire zero-centered load
/// ellipsoid. Unlike LoadSpace3's reference, this scale IS a hard requirement.
/// Basis columns carry amplitudes/correlations and stay fixed throughout a run.
#[derive(Debug, Clone, PartialEq)]
pub struct LoadSpaceCap3 { policy: LoadSpace3 }
impl LoadSpaceCap3 {
    pub fn new(dimension: usize, cap: f64, temperature: f64) -> Result<Self, &'static str> {
        let policy = LoadSpace3::new(dimension, cap, temperature)?;
        if policy.bias_bound() >= 1.0 {
            return Err("load-space smoothing consumes the entire compliance margin");
        }
        Ok(Self { policy })
    }
    #[must_use] pub fn cap(&self) -> f64 { self.policy.reference() }
    #[must_use] pub fn policy(&self) -> &LoadSpace3 { &self.policy }

    /// Numerical admission of a COMPLETE production load-space evaluation.
    /// Check physical worst-force work directly, as division can round a
    /// violating ratio to one. Neither AL tolerance nor smoothing waives it.
    /// The eigensolver's admitted maximum must also be <=1, not just the basis
    /// diagonals. This remains numerical evidence, not a rigorous enclosure.
    #[must_use]
    pub fn admits(&self, e: &LoadSpaceEvaluation3, smooth_tolerance: f64) -> bool {
        self.admits_spectrum(&e.policy, &e.spectrum, e.worst_compliance, smooth_tolerance)
    }

    // Both gray and binary evaluations use the same physical requirement.
    // This is numerical result admission, not a certificate for caller-edited data.
    fn admits_spectrum(&self, policy: &LoadSpace3, spectrum: &LoadSpaceSpectrum3,
        worst_compliance: f64, smooth_tolerance: f64) -> bool {
        policy == &self.policy && smooth_tolerance.is_finite()
            && (0.0..1.0).contains(&smooth_tolerance)
            && worst_compliance.is_finite() && worst_compliance >= 0.0
            && worst_compliance <= self.cap()
            && spectrum.maximum.is_finite() && (0.0..=1.0).contains(&spectrum.maximum)
            && spectrum.smooth_maximum.is_finite()
            && spectrum.smooth_maximum - 1.0 <= smooth_tolerance
    }
}

// Only shared AL machinery consumes this internal option shape. There is NO
// additional volume constraint in the minimum-material problem: full-domain
// volume here is merely the common initializer's admission ceiling.
fn engine_options(options: LoadSpaceLimitedOptions3) -> EnvelopeDesignOptions3 {
    EnvelopeDesignOptions3 { robustness: options.robustness, volume_limit: 1.0,
        density_floor: options.density_floor, gradient_check: options.gradient_check,
        optimizer: options.optimizer }
}

/// Minimize dilated volume while the smooth upper spectral maximum is <=1
/// and retained feasible results satisfy the actual combined-load cap.
/// Uses the existing load-space evaluator, gradient gate and projected AL loop.
/// No basis-direction sampling, cached worst-direction derivative, new optimizer
/// or weakening of physical requirements. Initial controls may be infeasible.
///
/// `best_feasible` means least DILATED volume, not least compliance. KKT values
/// describe the final accepted AL iterate, which may differ from that incumbent.
/// Geometry, supports, load basis, material model, cap and temperature stay fixed.
/// All physical work is cumulative; errors leave incoming nominal scales intact.
/// Fixed-grid linear elasticity only, with no stress/buckling/continuum guarantee.
pub struct LoadSpaceLimitedDesign3<'a, 'callback, O: Sdf3Elasticity> {
    inner: LoadSpaceDesignStudy3<'a, 'callback, O>,
}
impl<'a, 'callback, O: Sdf3Elasticity> LoadSpaceLimitedDesign3<'a, 'callback, O> {
    pub fn new(study: &'a mut CutDensityStudy3<O>, basis: &'a [&'a [f64]], rho: &[f64],
        cap: LoadSpaceCap3, options: LoadSpaceLimitedOptions3, control: &'a mut SolveControl<'callback>)
        -> Result<Self, ProjectedAlError<EvaluationStop>> {
        let policy = cap.policy.clone();
        let inner = LoadSpaceDesignStudy3::initialize(study, basis, rho, policy,
            engine_options(options), Some(cap), control)?;
        Ok(Self { inner })
    }
    /// Restore accepted AL state under the same force ellipsoid and physical
    /// cap. Fresh accepted/incumbent solves spend the original allowances;
    /// neither earlier search steps nor the original gradient gate are replayed.
    pub fn restore(study: &'a mut CutDensityStudy3<O>, basis: &'a [&'a [f64]],
        checkpoint: LoadSpaceDesignCheckpoint3, cap: LoadSpaceCap3,
        options: LoadSpaceLimitedOptions3, control: &'a mut SolveControl<'callback>)
        -> Result<Self, ProjectedAlError<EvaluationStop>> {
        let inner = LoadSpaceDesignStudy3::restore_family(study, basis, checkpoint,
            cap.policy.clone(), engine_options(options), Some(cap), control)?;
        Ok(Self { inner })
    }
    #[must_use]
    pub fn checkpoint(&self) -> LoadSpaceDesignCheckpoint3 { self.inner.checkpoint() }
    #[must_use] pub fn cap(&self) -> &LoadSpaceCap3 { self.inner.cap.as_ref().expect("cap-bound engine") }
    #[must_use] pub fn accepted(&self) -> &LoadSpaceDesignPoint3 { self.inner.accepted() }
    #[must_use] pub fn best_feasible(&self) -> Option<&LoadSpaceDesignPoint3> { self.inner.best_feasible() }
    /// Feasibility here means the combined-load requirement, not a volume cap.
    #[must_use] pub fn history(&self) -> &[LoadSpaceDesignIteration3] { self.inner.history() }
    #[must_use] pub fn gradient_check(&self) -> &LoadSpaceGradientCheck3 { self.inner.gradient_check() }
    #[must_use] pub fn work(&self) -> SolveWork { self.inner.work() }
    #[must_use] pub fn optimizer_work(&self) -> ProjectedAlWork { self.inner.optimizer_work() }
    /// Read the fixed geometry and physical field layout of accepted results.
    /// This cannot mutate the optimizer's model or its cumulative work budget.
    #[must_use] pub fn study(&self) -> &CutDensityStudy3<O> { self.inner.study }
    pub fn run(&mut self, additional_steps: usize) -> Result<ProjectedAlReport, ProjectedAlError<EvaluationStop>> {
        self.inner.run(additional_steps)
    }

    /// Fresh full-basis solve, cross-compliance spectrum and worst-load witness.
    /// This is not a check of a previously selected direction alone. Charge the
    /// original physical budget without changing AL state/history or its counters.
    pub fn recheck_best(&mut self) -> Result<Option<LoadSpaceDesignPoint3>, EvaluationStop> {
        let inner = &mut self.inner;
        let Some(best) = &inner.best else { return Ok(None); };
        inner.control.checkpoint("sdf3-load-space-cap-recheck")?;
        let evaluation = inner.study.evaluate_robust_load_space(&best.rho, inner.basis,
            inner.options.robustness, &inner.policy, inner.control)?;
        if !inner.cap.as_ref().expect("cap-bound engine").admits(&evaluation, inner.options.optimizer.tolerance) {
            return Err(failure("sdf3-load-space-cap-recheck"));
        }
        Ok(Some(LoadSpaceDesignPoint3 { rho: best.rho.clone(), evaluation }))
    }
}

#[cfg(test)]
#[path = "limits_tests.rs"]
mod tests;
