//! Frozen linear support and the combined adhesive potential barrier.
use super::{
    CancelGate, JkrNormalMotion, JkrNormalMotionError, JkrNormalMotionState,
    SourceJkrSphereContact, finite, poll,
};

pub(super) struct Support {
    pub stiffness: f64,
    pub anchor: f64,
    minimum_indentation: f64,
    minimum_reaction: f64,
    fold_reaction: f64,
}
impl Support {
    pub fn new(
        normal: &SourceJkrSphereContact,
        stiffness: f64,
        anchor: f64,
        iterations: usize,
    ) -> Result<Self, JkrNormalMotionError> {
        if !stiffness.is_finite() || stiffness < 0.0 || !anchor.is_finite() {
            return Err(JkrNormalMotionError::Invalid(
                "finite nonnegative elastic support",
            ));
        }
        let fold = normal.minimum_attached_indentation_m();
        let x = 3.0_f64.cbrt();
        let original_minimum = normal.indentation_scale_m() * (x * x * x * x - 4.0 * x);
        let minimum = if stiffness == 0.0 {
            original_minimum
        } else {
            // The original JKR tangent strictly increases from -infinity at
            // the fold to zero at its force minimum. Adding a positive spring
            // leaves one combined force minimum. Resolve it once, cold.
            let mut lo = fold;
            let mut hi = original_minimum;
            let tolerance = resolution(normal, fold);
            let mut result = None;
            for _ in 0..iterations {
                let mid = lo.midpoint(hi);
                if mid == fold {
                    return Err(JkrNormalMotionError::Unresolved(
                        "support minimum at source resolution",
                    ));
                }
                let tangent = normal
                    .evaluate_attached(mid)?
                    .tangent_n_per_m
                    .ok_or(JkrNormalMotionError::Unresolved("support minimum tangent"))?;
                let slope = finite(tangent + stiffness, "supported contact tangent")?;
                if slope < 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
                if 0.5 * (hi - lo) <= tolerance || mid == lo.midpoint(hi) {
                    result = Some(lo.midpoint(hi));
                    break;
                }
            }
            result.ok_or(JkrNormalMotionError::Unresolved(
                "support minimum iterations",
            ))?
        };
        if minimum <= fold {
            return Err(JkrNormalMotionError::Unresolved(
                "unrepresentable support minimum",
            ));
        }
        let mut support = Self {
            stiffness,
            anchor,
            minimum_indentation: minimum,
            minimum_reaction: 0.0,
            fold_reaction: 0.0,
        };
        support.minimum_reaction = finite(
            normal.evaluate_attached(minimum)?.normal_force_n + support.force(minimum)?,
            "minimum supported reaction",
        )?;
        support.fold_reaction = finite(
            normal.evaluate_attached(fold)?.normal_force_n + support.force(fold)?,
            "fold supported reaction",
        )?;
        Ok(support)
    }

    pub fn force(&self, indentation: f64) -> Result<f64, JkrNormalMotionError> {
        if self.stiffness == 0.0 {
            return Ok(0.0);
        }
        finite(
            self.stiffness * (indentation - self.anchor),
            "support force",
        )
    }
    pub fn mean(&self, from: f64, to: f64) -> Result<f64, JkrNormalMotionError> {
        self.force(from.midpoint(to))
    }
    pub fn energy(&self, indentation: f64) -> Result<f64, JkrNormalMotionError> {
        if self.stiffness == 0.0 {
            return Ok(0.0);
        }
        finite(
            0.5 * self.force(indentation)? * (indentation - self.anchor),
            "support energy",
        )
    }
    pub fn work(&self, from: f64, to: f64) -> Result<f64, JkrNormalMotionError> {
        finite(self.mean(from, to)? * (to - from), "support energy change")
    }

    /// Locate at most one unstable equilibrium per held-load advance; substeps
    /// reuse this value. It is the maximum of U_JKR + U_support - P*delta.
    pub fn barrier(
        &self,
        normal: &SourceJkrSphereContact,
        load: f64,
        iterations: usize,
        gate: &CancelGate,
    ) -> Result<Option<f64>, JkrNormalMotionError> {
        if self.stiffness == 0.0 {
            let scaled = load / normal.force_scale_n;
            if !(-9.0..=-5.0).contains(&scaled) {
                return Ok(None);
            }
            let x = (3.0 - (9.0 + scaled).sqrt()).cbrt();
            return Ok(Some(
                normal.indentation_scale_m() * (x * x * x * x - 4.0 * x),
            ));
        }
        if load <= self.minimum_reaction || load >= self.fold_reaction {
            return Ok(None);
        }
        let mut lo = normal.minimum_attached_indentation_m();
        let mut hi = self.minimum_indentation;
        let mut point = lo.midpoint(hi);
        for _ in 0..iterations {
            poll(gate)?;
            let response = normal.evaluate_attached(point)?;
            let residual = finite(
                response.normal_force_n + self.force(point)? - load,
                "supported barrier residual",
            )?;
            if residual == 0.0 {
                return Ok(Some(point));
            }
            if residual > 0.0 {
                lo = point;
            } else {
                hi = point;
            }
            if 0.5 * (hi - lo) <= resolution(normal, point) {
                return Ok(Some(lo.midpoint(hi)));
            }
            let slope = finite(
                response
                    .tangent_n_per_m
                    .ok_or(JkrNormalMotionError::Unresolved(
                        "supported barrier tangent",
                    ))?
                    + self.stiffness,
                "supported barrier slope",
            )?;
            let newton = point - residual / slope;
            point = if slope < 0.0 && newton > lo && newton < hi {
                newton
            } else {
                lo.midpoint(hi)
            };
        }
        Err(JkrNormalMotionError::Unresolved(
            "supported barrier iterations",
        ))
    }
}

fn resolution(normal: &SourceJkrSphereContact, point: f64) -> f64 {
    8.0 * f64::EPSILON * point.abs().max(normal.indentation_scale_m())
}

pub(super) fn check_barrier(
    law: &JkrNormalMotion<'_>,
    state: JkrNormalMotionState,
    load: f64,
    low: f64,
    high: f64,
    barrier: Option<f64>,
) -> Result<(), JkrNormalMotionError> {
    if let Some(critical) = barrier.filter(|d| *d > low && *d < high) {
        let from = state.contact.indentation_m();
        let needed = finite(
            law.normal.reversible_work_between(from, critical)?
                + law.support.work(from, critical)?
                - load * (critical - from),
            "supported adhesive barrier work",
        )?;
        let kinetic = finite(
            0.5 * law.mass_kg * state.velocity_m_s * state.velocity_m_s,
            "normal barrier kinetic energy",
        )?;
        let tolerance = law.budget.energy_absolute_tolerance_j
            + law.budget.energy_relative_tolerance * needed.abs().max(kinetic);
        if needed > kinetic + tolerance {
            return Err(JkrNormalMotionError::Unresolved(
                "unresolved adhesive turning barrier",
            ));
        }
    }
    Ok(())
}
