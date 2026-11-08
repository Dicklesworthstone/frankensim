//! Restitution calibration and derivatives of the unilateral contact law.
use super::{DContactError, Obstacle, det};

/// Calibrate Hunt–Crossley `chi` [s/m] from normal restitution at an approach
/// speed [m/s]. The returned coefficient remains fixed during subsequent impacts.
///
/// For an isolated, initially uncompressed normal impact with positive effective
/// mass and stiffness, `F = K delta^alpha (1 + chi * delta_dot)` gives the exact
/// relation `x - ln(1+x) = -e*x - ln(1-e*x)`, where `x = chi * approach_speed`
/// and `e` is restitution. This solves its unique positive root; the mass,
/// stiffness, and exponent cancel. It does not use a weak-damping approximation.
///
/// Restitution is matched at the supplied speed only. External loading, coupled
/// modes and simultaneous contacts can change the rebound. This is parameter
/// calibration, not an impact impulse or a replacement for the dissipative port.
///
/// `e = 1` returns exactly zero loss. Perfectly inelastic `e = 0` would require
/// an infinite coefficient and is refused. Allocation-free, deterministic, and
/// intended for setup rather than the timestep loop.
///
/// # Errors
/// [`DContactError::Parameter`] for nonfinite values, `e` outside `(0, 1]`, a
/// nonpositive speed, or a coefficient outside the representable positive range.
pub fn hunt_crossley_loss_from_restitution(
    restitution: f64,
    approach_speed_m_s: f64,
) -> Result<f64, DContactError> {
    if !restitution.is_finite() || restitution <= 0.0 || restitution > 1.0 {
        return Err(DContactError::Parameter {
            what: "Hunt-Crossley restitution must be finite and in (0, 1]",
        });
    }
    if !approach_speed_m_s.is_finite() || approach_speed_m_s <= 0.0 {
        return Err(DContactError::Parameter {
            what: "Hunt-Crossley reference approach speed must be finite and positive",
        });
    }
    if restitution == 1.0 {
        return Ok(0.0);
    }

    let (mut lo, mut hi) = (0.0, 1.0 / restitution);
    if !hi.is_finite() {
        return Err(DContactError::Parameter {
            what: "Hunt-Crossley dimensionless loss is not representable",
        });
    }
    // Dividing out the double root at zero keeps nearly elastic impacts from
    // converging spuriously to chi=0. The residual is positive at zero and
    // negative at the singular unloading limit x=1/e.
    for _ in 0..128 {
        let x = f64::midpoint(lo, hi);
        if x == lo || x == hi {
            break;
        }
        if restitution_residual(restitution, x) > 0.0 {
            lo = x;
        } else {
            hi = x;
        }
    }
    let chi = f64::midpoint(lo, hi) / approach_speed_m_s;
    if !chi.is_finite() || chi <= 0.0 {
        return Err(DContactError::Parameter {
            what: "Hunt-Crossley calibrated chi is not representable",
        });
    }
    Ok(chi)
}

fn restitution_residual(e: f64, x: f64) -> f64 {
    if x <= 0.25 {
        // F(x)/x^2 = sum_n [(-1)^n - e^(n+2)] x^n/(n+2).
        // Factoring 1-e^2 preserves measured loss close to e=1. At x<=1/4,
        // the omitted tail after n=32 is below 2e-21 in absolute magnitude.
        let mut residual = 0.5 * (1.0 - e) * (1.0 + e);
        let (mut x_power, mut e_power) = (1.0, e * e);
        for n in 1..=32 {
            x_power *= x;
            e_power *= e;
            let coefficient = if n % 2 == 0 {
                1.0 - e_power
            } else {
                -1.0 - e_power
            };
            residual += coefficient * x_power / f64::from(n + 2);
        }
        residual
    } else {
        let unloading = 1.0 - e * x;
        if unloading <= 0.0 {
            return f64::NEG_INFINITY;
        }
        (1.0 + e) * x + det::ln(unloading) - det::ln(1.0 + x)
    }
}

impl Obstacle {
    /// Set a fixed Hunt–Crossley loss coefficient from restitution measured at
    /// the given normal approach speed. See [`hunt_crossley_loss_from_restitution`]
    /// for the isolated-impact assumptions and speed dependence. The obstacle's
    /// caller-supplied provenance is retained; this does not certify a source.
    ///
    /// # Errors
    /// [`DContactError::Parameter`] for invalid or unrepresentable calibration.
    pub fn with_restitution_loss(
        self,
        restitution: f64,
        approach_speed_m_s: f64,
    ) -> Result<Self, DContactError> {
        self.with_internal_loss(hunt_crossley_loss_from_restitution(
            restitution,
            approach_speed_m_s,
        )?)
    }

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
