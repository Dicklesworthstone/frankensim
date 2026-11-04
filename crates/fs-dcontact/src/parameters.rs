//! Parameter derivatives of the original total unilateral resisting force.
use super::{DContactError, Obstacle, det};

impl Obstacle {
    /// Contract a generalized force seed with derivatives of the TOTAL positive
    /// resisting force, including the original Hunt--Crossley loss increment.
    /// Output order is `[K, chi, gap_0, ..., gap_(n_points-1)]`. Exponent,
    /// collocation and physical quadrature measures remain fixed.
    ///
    /// `x` has the original interleaved `[q, momentum, ...]` layout and may have
    /// additional finite internal coordinates. Velocity and seed each contain
    /// `n_modes` generalized entries; no inner-storage force is differentiated.
    /// With `delta = [b.q-gap]_+`, `rate = b.velocity`, and `s = b.seed`, each
    /// uncapped site contributes `s*w*delta^alpha*(1+chi*rate)` to the K
    /// derivative, `s*w*K*delta^alpha*rate` to chi, and
    /// `-s*w*K*alpha*delta^(alpha-1)*(1+chi*rate)` to its gap.
    ///
    /// Separation (`delta <= 0`) and strictly capped unloading (`chi*rate < -1`)
    /// select zero derivatives. At the cap equality the unsaturated one-sided
    /// derivative is selected, matching the existing state tangent. The K=0 and
    /// chi=0 boundaries retain their admitted-side derivatives without division
    /// by a coefficient or measure.
    ///
    /// Allocation-free. Shape, malformed-owner and nonfinite-input refusals
    /// leave output untouched. Arithmetic refusal invalidates all output scratch;
    /// callers must discard it on any error. Success overwrites every entry.
    pub fn resistance_parameter_vjp_into(
        &self,
        n_modes: usize,
        x: &[f64],
        velocity: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), DContactError> {
        if self.n_points.checked_mul(n_modes) != Some(self.collocation.len())
            || self.gaps.len() != self.n_points
            || self.weights.len() != self.n_points
            || n_modes.checked_mul(2).is_none_or(|prefix| x.len() < prefix)
            || velocity.len() != n_modes
            || seed.len() != n_modes
            || self.n_points.checked_add(2) != Some(output.len())
        {
            return Err(DContactError::Shape {
                what: "contact parameter VJP state/collocation/velocity/seed/output length",
            });
        }
        if !self.stiffness.is_finite()
            || self.stiffness < 0.0
            || !self.alpha.is_finite()
            || self.alpha < 1.0
            || !self.internal_loss.is_finite()
            || self.internal_loss < 0.0
            || self.weights.iter().any(|w| *w < 0.0)
            || self
                .collocation
                .iter()
                .chain(&self.gaps)
                .chain(&self.weights)
                .chain(x)
                .chain(velocity)
                .chain(seed)
                .any(|v| !v.is_finite())
        {
            return Err(DContactError::Parameter {
                what: "invalid contact parameter VJP owner or nonfinite input",
            });
        }

        output.fill(0.0);
        for i in 0..self.n_points {
            let delta = self.penetration_at(n_modes, x, i);
            if !delta.is_finite() {
                return Err(derivative_overflow());
            }
            if delta <= 0.0 {
                continue;
            }
            let row = &self.collocation[i * n_modes..(i + 1) * n_modes];
            let mut rate = 0.0;
            let mut projected_seed = 0.0;
            for (k, b) in row.iter().enumerate() {
                rate += b * velocity[k];
                projected_seed += b * seed[k];
            }
            let loss_factor = self.internal_loss * rate;
            if ![rate, projected_seed, loss_factor]
                .iter()
                .all(|v| v.is_finite())
            {
                return Err(derivative_overflow());
            }
            if loss_factor < -1.0 {
                continue;
            }
            let force_factor = 1.0 + loss_factor;
            let weighted_power = self.weights[i] * det::pow(delta, self.alpha);
            let elastic = self.stiffness * weighted_power;
            let slope =
                self.weights[i] * self.stiffness * self.alpha * det::pow(delta, self.alpha - 1.0);
            output[0] += projected_seed * weighted_power * force_factor;
            output[1] += projected_seed * elastic * rate;
            output[2 + i] = -projected_seed * slope * force_factor;
        }
        if output.iter().any(|v| !v.is_finite()) {
            return Err(derivative_overflow());
        }
        Ok(())
    }
}

fn derivative_overflow() -> DContactError {
    DContactError::Parameter {
        what: "contact parameter derivative overflowed",
    }
}
