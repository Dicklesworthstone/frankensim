//! Authored finite-counterflow production onset for a seeded Vinen tangle.
//!
//! The steady imposed-flow fit `sqrt(L)=gamma*(W-W0)` is discussed in
//! Baggaley and Laizet, Phys. Fluids 25, 115101 (2013), equation (8),
//! <https://arxiv.org/html/1310.5890v1>. It does not uniquely determine a
//! transient law. HERE the explicit phenomenological choice is
//! `dL/dt=a*max(W-W0,0)*L^(3/2)-b*L^2`, with unchanged mutual friction.
//! W0 is measured/authored for the application, not predicted from helium
//! temperature, geometry, or a universal Landau critical velocity.

use super::{MaterialError, VinenTangle, decay, invalid, nonnegative, positive};

impl VinenTangle {
    /// Return this closure with a finite nonnegative counterflow onset [m/s].
    ///
    /// Production uses `max(|v_n-v_s|-onset,0)`; existing lines continue to
    /// decay and exert friction below onset. A zero seed remains zero even
    /// above onset. This is neither spontaneous nucleation nor a prediction
    /// of the TI/TII transition. Zero retains the original arithmetic.
    ///
    /// Reducing production leaves the original global heat-admission
    /// inequality sufficient. Only the immutable closure changes; no channel
    /// history or line density is changed by constructing this value.
    pub fn with_counterflow_onset_m_s(mut self, onset_m_s: f64) -> Result<Self, MaterialError> {
        if !nonnegative(onset_m_s) {
            return Err(invalid(
                "Vinen production onset must be finite and nonnegative",
            ));
        }
        self.production_onset_m_s = onset_m_s;
        Ok(self)
    }

    pub const fn counterflow_onset_m_s(&self) -> f64 {
        self.production_onset_m_s
    }

    /// Nonnegative source equilibrium at an EXTERNALLY HELD relative speed.
    ///
    /// This is `(a/b)^2 * max(W-W0,0)^2`, not an equilibrium of the freely
    /// decaying coupled velocity/density system. It is useful for comparing
    /// an authored coefficient set to steady counterflow measurements.
    pub fn imposed_counterflow_equilibrium_per_m2(
        &self,
        relative_speed_m_s: f64,
    ) -> Result<f64, MaterialError> {
        if !nonnegative(relative_speed_m_s) {
            return Err(invalid(
                "imposed Vinen counterflow speed must be finite and nonnegative",
            ));
        }
        let excess = (relative_speed_m_s - self.production_onset_m_s).max(0.0);
        if excess == 0.0 || self.parameters.production_coefficient == 0.0 {
            return Ok(0.0);
        }
        let root = self.parameters.production_coefficient
            * (excess / self.parameters.decay_coefficient_m2_s);
        let density = root * root;
        if !nonnegative(density)
            || (self.parameters.production_coefficient > 0.0 && excess > 0.0 && density == 0.0)
        {
            return Err(invalid("imposed Vinen equilibrium is unrepresentable"));
        }
        Ok(density)
    }

    /// Positive, thermodynamically admitted exponential-midpoint source.
    /// Each local map uses a constant a_eff in [0,a] and therefore retains
    /// the original nonnegative heat proof. Sampling a_eff and the physical
    /// clock at the predicted midpoint gives second order in smooth regions.
    /// The authored rate fraction controls refinement at the onset kink.
    pub(super) fn exposure_with_onset(
        &self,
        initial_w: f64,
        initial_root: f64,
        duration: f64,
    ) -> Result<(f64, usize, f64, f64), MaterialError> {
        let mut elapsed = 0.0;
        let mut exposure = 0.0;
        let mut exposure_correction = 0.0;
        let mut root = initial_root;
        let mut root_change = 0.0;
        let mut root_correction = 0.0;
        for substep in 1..=self.controls.maximum_substeps {
            let w = decay::attenuate(initial_w, self.drag_per_line_m2_s * exposure);
            let remaining = duration - elapsed;
            let mut local = *self;
            let (h, increment, next_root, change) = if w <= self.production_onset_m_s
                || self.parameters.production_coefficient == 0.0
            {
                // W can only decrease in a source step, so once production
                // stops it stays stopped. Finish the entire tail exactly.
                local.parameters.production_coefficient = 0.0;
                let (increment, _) = local.exposure(w, root, remaining)?;
                let (next_root, change) = local.root_at_exposure(w, root, increment)?;
                (remaining, increment, next_root, change)
            } else {
                let density = root * root;
                let b = self.parameters.decay_coefficient_m2_s;
                let a = self.parameters.production_coefficient;
                let rate = (self.drag_per_line_m2_s * density)
                    .max(0.5 * b * density)
                    .max(0.5 * a * (w - self.production_onset_m_s) * root);
                if !positive(rate) || !positive(density) {
                    return Err(invalid(
                        "Vinen onset source rate or density is unrepresentable",
                    ));
                }
                let h = remaining.min(self.controls.maximum_fractional_step / rate);
                if !positive(h) || elapsed + h <= elapsed {
                    return Err(invalid(
                        "Vinen onset source cannot advance its physical clock",
                    ));
                }
                local.parameters.production_coefficient = a * ((w - self.production_onset_m_s) / w);
                let half_exposure = 0.5 * h * density;
                let (mid_root, _) = local.root_at_exposure(w, root, half_exposure)?;
                let mid_w = decay::attenuate(w, self.drag_per_line_m2_s * half_exposure);
                local.parameters.production_coefficient = if mid_w > self.production_onset_m_s {
                    a * ((mid_w - self.production_onset_m_s) / mid_w)
                } else {
                    0.0
                };
                let increment = h * mid_root * mid_root;
                let (next_root, change) = local.root_at_exposure(w, root, increment)?;
                (h, increment, next_root, change)
            };
            let next_elapsed = if h == remaining {
                duration
            } else {
                elapsed + h
            };
            if !positive(increment) || next_elapsed <= elapsed {
                return Err(invalid(
                    "Vinen onset source clock or exposure cannot advance",
                ));
            }
            let corrected = increment - exposure_correction;
            let next_exposure = exposure + corrected;
            exposure_correction = (next_exposure - exposure) - corrected;
            if !positive(next_exposure) || next_exposure <= exposure {
                return Err(invalid("Vinen onset integrated density is unrepresentable"));
            }
            let corrected_root = change - root_correction;
            let next_change = root_change + corrected_root;
            root_correction = (next_change - root_change) - corrected_root;
            if !next_change.is_finite() {
                return Err(invalid("Vinen onset line increment is unrepresentable"));
            }
            exposure = next_exposure;
            root = next_root;
            root_change = next_change;
            elapsed = next_elapsed;
            if elapsed == duration {
                return Ok((exposure, substep, root, root_change));
            }
        }
        Err(invalid(
            "Vinen onset source exhausted maximum_substeps; refine duration or increase the explicit budget",
        ))
    }
}
