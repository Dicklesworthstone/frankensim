//! Normalized primary-creep kinetics on the original Norton history owner.
//!
//! `rate(q,T,alpha) = rate_Norton(q,T) h(alpha)` with
//! `h(alpha) = f + (1-f) (1+alpha/a)^(-m)`. The positive reference strain `a`
//! removes the virgin-state singularity; `m >= 0` slows primary creep and the
//! explicitly selected `f` sets the asymptotic secondary/initial rate ratio.
//! This is an empirical strain-hardening power law, not energetic plastic
//! hardening: no extra recoverable energy is inferred from the rate reduction.
//! See COMSOL's Creep and Viscoplasticity theory, hardening-function section:
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.sme/sme_ug_theory.06.039.html>.
//! The normalized shift and optional residual-rate term are explicit model
//! choices; they require their own calibration. No cyclic recovery, tertiary
//! damage, rupture, or material-specific predictive accuracy is implied.

use super::{CreepError, NortonCreep, Relaxation, Result, logistic_inverse, positive, require};
use fs_math::det;

fn log_one_plus(x: f64) -> f64 {
    if x < 1e-4 {
        // Omitted term <= 2e-21. Keep tiny shifts before a large exponent
        // multiplies them instead of rounding 1+x to exactly one.
        x * (1.0 + x * (-0.5 + x * (1.0 / 3.0 - x * 0.25)))
    } else {
        det::ln(1.0 + x)
    }
}

/// Dimensionless primary-creep parameters. The original Norton reference rate
/// is the rate at zero accumulated creep, not the eventual secondary rate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StrainHardeningCreep {
    /// Positive regularizing equivalent-strain scale.
    pub reference_strain: f64,
    /// Nonnegative exponent of the decreasing strain-hardening multiplier.
    pub exponent: f64,
    /// Asymptotic/initial rate ratio in [0,1]. Zero gives pure primary creep;
    /// one, or exponent zero, reduces exactly to the original Norton law.
    pub residual_rate_fraction: f64,
}

impl StrainHardeningCreep {
    pub(crate) fn validate(&self) -> Result<()> {
        require(
            positive(self.reference_strain)
                && self.exponent.is_finite()
                && self.exponent >= 0.0
                && self.residual_rate_fraction.is_finite()
                && (0.0..=1.0).contains(&self.residual_rate_fraction),
            "primary creep needs positive reference strain, exponent >= 0 and residual fraction in [0,1]",
        )
    }

    pub(crate) fn active(&self) -> bool {
        self.exponent > 0.0 && self.residual_rate_fraction < 1.0
    }

    // Return log(h) and the fraction of h supplied by its decaying term.
    pub(crate) fn value(&self, alpha: f64) -> Result<(f64, f64)> {
        if !self.active() {
            return Ok((0.0, 0.0));
        }
        if alpha == 0.0 {
            return Ok((0.0, 1.0 - self.residual_rate_fraction));
        }
        let ratio = alpha / self.reference_strain;
        let log_shift = if ratio.is_finite() {
            log_one_plus(ratio)
        } else {
            det::ln(alpha) - det::ln(self.reference_strain)
        };
        let log_primary = -self.exponent * log_shift;
        let f = self.residual_rate_fraction;
        let (log_h, weight) = if f == 0.0 {
            (log_primary, 1.0)
        } else {
            let a = det::ln(1.0 - f) + log_primary;
            let b = det::ln(f);
            let shift = a.max(b);
            let primary = det::exp(a - shift);
            let sum = primary + det::exp(b - shift);
            (shift + det::ln(sum), primary / sum)
        };
        require(
            log_h.is_finite(),
            "primary creep multiplier is not representable",
        )?;
        Ok((log_h, weight))
    }

    // -increment * d(log h)/d alpha. Scaling the ratio first avoids both an
    // overflowing alpha+a and a derivative overflow when a is very small.
    pub(crate) fn sensitivity(&self, alpha: f64, increment: f64, weight: f64) -> f64 {
        let scale = alpha.max(self.reference_strain);
        let fraction = (increment / scale) / (alpha / scale + self.reference_strain / scale);
        (self.exponent * fraction) * weight
    }
}

impl NortonCreep {
    /// Select an explicit primary/secondary strain-hardening model. Retain this
    /// immutable selection when preparing another physical interval. Sourced
    /// laws must use the corresponding source constructor to retain all claims.
    pub fn with_strain_hardening(mut self, parameters: StrainHardeningCreep) -> Result<Self> {
        parameters.validate()?;
        self.strain_hardening = Some(parameters);
        Ok(self)
    }

    /// The immutable optional primary-creep selection, including exact Norton
    /// reductions. History continuation must preserve these parameters.
    pub const fn strain_hardening(&self) -> Option<&StrainHardeningCreep> {
        self.strain_hardening.as_ref()
    }

    pub(super) fn hardening_value(&self, alpha: f64) -> Result<(f64, f64)> {
        match &self.strain_hardening {
            Some(parameters) => parameters.value(alpha),
            None => Ok((0.0, 0.0)),
        }
    }

    /// Shared tensor/axial/joint return. With monotone decreasing h, the
    /// equation q + modulus*dt*rate(q,alpha_old+delta) = q_trial has one root.
    /// Solve for z=log(modulus*delta/q). Logistic coordinates retain BOTH
    /// tiny relaxed stress and tiny creep history, even when their companion
    /// rounds to its trial value. No accepted history is mutated.
    pub(super) fn relaxation_from(
        &self,
        q_trial: f64,
        log_step_scale: f64,
        log_modulus: f64,
        old_alpha: f64,
    ) -> Result<Relaxation> {
        let Some(parameters) = self.strain_hardening.filter(StrainHardeningCreep::active) else {
            return self.relaxation(q_trial, log_step_scale);
        };
        let returned = parameters.relaxation(
            q_trial,
            log_step_scale,
            log_modulus,
            old_alpha,
            self.parameters.stress_exponent,
            self.log_reference_stress,
        )?;
        Ok(Relaxation {
            q: returned.q,
            delta_alpha: returned.delta_alpha,
            ratio: returned.ratio,
            radial: returned.radial,
        })
    }
}

/// Scalar kinetic return shared by additive and multiplicative creep owners.
pub(crate) struct PrimaryCreepReturn {
    pub(crate) q: f64,
    pub(crate) delta_alpha: f64,
    pub(crate) ratio: f64,
    pub(crate) radial: f64,
    /// Stable dq/dq_trial - q/q_trial, including positive primary resistance.
    pub(crate) radial_correction: f64,
}

impl StrainHardeningCreep {
    /// Solve the same dimensionless endpoint equation for an explicitly supplied
    /// elastic modulus and stress normalization. Stress measure and kinematics
    /// remain owned by the caller; Kirchhoff and Cauchy fits are not interchanged.
    pub(crate) fn relaxation(
        &self,
        q_trial: f64,
        log_step_scale: f64,
        log_modulus: f64,
        old_alpha: f64,
        n: f64,
        log_reference_stress: f64,
    ) -> Result<PrimaryCreepReturn> {
        let parameters = self;
        require(
            log_modulus.is_finite(),
            "primary creep stiffness is not representable",
        )?;
        let old_log_h = parameters.value(old_alpha)?.0;
        if q_trial == 0.0 {
            let ratio = if n == 1.0 {
                logistic_inverse(log_step_scale + old_log_h)
            } else {
                1.0
            };
            return Ok(PrimaryCreepReturn {
                q: 0.0,
                delta_alpha: 0.0,
                ratio,
                radial: ratio,
                radial_correction: 0.0,
            });
        }
        let log_q_trial = det::ln(q_trial);
        let log_compliance = log_q_trial - log_modulus;
        let log_k = log_step_scale + (n - 1.0) * (log_q_trial - log_reference_stress);
        let log_initial_k = log_k + old_log_h;
        require(
            log_initial_k.is_finite(),
            "primary creep step is not representable",
        )?;
        let evaluate = |z: f64| -> Result<(f64, f64, f64, f64, f64, f64, f64)> {
            let (log_ratio, log_fraction) = if z >= 0.0 {
                let correction = log_one_plus(det::exp(-z));
                (-z - correction, -correction)
            } else {
                let correction = log_one_plus(det::exp(z));
                (-correction, z - correction)
            };
            let delta = det::exp(log_compliance + log_fraction);
            let alpha = old_alpha + delta;
            require(
                alpha.is_finite(),
                "primary creep history is not representable",
            )?;
            let (log_h, weight) = parameters.value(alpha)?;
            // Enforce the HISTORY equation. A sub-ulp stress change can still
            // harden a material with a correspondingly small reference strain.
            let residual = log_fraction - n * log_ratio - log_k - log_h;
            let roundoff = 8.0
                * f64::EPSILON
                * (1.0 + log_fraction.abs() + (n * log_ratio).abs() + log_k.abs() + log_h.abs());
            require(
                residual.is_finite() && roundoff.is_finite(),
                "primary creep residual is not representable",
            )?;
            let beta = parameters.sensitivity(alpha, delta, weight);
            let ratio = det::exp(log_ratio);
            let fraction = det::exp(log_fraction);
            let slope = 1.0 + (n - 1.0) * fraction + beta * ratio;
            Ok((residual, roundoff, ratio, delta, beta, slope, fraction))
        };
        // h(alpha_old) bounds the rate from above. Bounded bracket expansion
        // accounts for the further hardening inside this physical interval.
        let mut upper = if log_initial_k > 0.0 {
            log_initial_k / n
        } else {
            log_initial_k
        };
        let mut width = 1.0;
        let mut lower = upper - width;
        for _ in 0..32 {
            if evaluate(lower)?.0 <= 0.0 {
                break;
            }
            width *= 2.0;
            lower = upper - width;
        }
        require(
            evaluate(lower)?.0 <= 0.0,
            "primary creep root could not be bracketed",
        )?;
        let mut z = upper;
        for _ in 0..96 {
            let (residual, roundoff, ratio, delta_alpha, beta, slope, fraction) = evaluate(z)?;
            if residual.abs() <= roundoff {
                let q = q_trial * ratio;
                require(
                    positive(q) && delta_alpha.is_finite() && delta_alpha >= 0.0,
                    "primary creep stress or increment is not representable",
                )?;
                let radial = logistic_inverse(det::ln(n) + z - log_one_plus(beta));
                require(
                    radial.is_finite(),
                    "primary creep tangent is not representable",
                )?;
                let radial_correction = radial * fraction * ((beta - (n - 1.0)) / (1.0 + beta));
                require(
                    radial_correction.is_finite(),
                    "primary creep radial correction is not representable",
                )?;
                return Ok(PrimaryCreepReturn {
                    q,
                    delta_alpha,
                    ratio,
                    radial,
                    radial_correction,
                });
            }
            if residual > 0.0 {
                upper = z;
            } else {
                lower = z;
            }
            let candidate = z - residual / slope;
            z = if candidate.is_finite() && candidate > lower && candidate < upper && candidate != z
            {
                candidate
            } else {
                0.5 * lower + 0.5 * upper
            };
        }
        Err(CreepError::Invalid(
            "bounded primary creep root did not converge",
        ))
    }
}

impl NortonCreep {
    /// Perfect-plastic joint branch: q is fixed at yield, but its creep increment
    /// must still solve the implicit hardening equation. The branch test supplies
    /// a finite strict upper bound before any plastic increment can be accepted.
    pub(super) fn fixed_stress_increment(
        &self,
        q: f64,
        old_alpha: f64,
        maximum: f64,
    ) -> Result<f64> {
        let Some(parameters) = self.strain_hardening.filter(StrainHardeningCreep::active) else {
            return Ok(det::exp(self.log_duration + self.log_rate(q)));
        };
        require(
            positive(maximum),
            "primary creep increment bound must be positive",
        )?;
        let log_rate_step = self.log_duration + self.log_rate(q);
        let initial = log_rate_step + parameters.value(old_alpha)?.0;
        if det::exp(initial) == 0.0 {
            return Ok(0.0);
        }
        let mut lower = det::ln(f64::from_bits(1));
        let mut upper = det::ln(maximum);
        let mut y = initial.min(upper);
        for _ in 0..96 {
            let delta = det::exp(y);
            let alpha = old_alpha + delta;
            require(
                alpha.is_finite(),
                "primary creep history is not representable",
            )?;
            let (log_h, weight) = parameters.value(alpha)?;
            let residual = y - log_rate_step - log_h;
            let roundoff = 8.0 * f64::EPSILON * (1.0 + y.abs() + log_rate_step.abs() + log_h.abs());
            if residual.abs() <= roundoff {
                return Ok(delta);
            }
            if residual > 0.0 {
                upper = y;
            } else {
                lower = y;
            }
            let slope = 1.0 + parameters.sensitivity(alpha, delta, weight);
            let candidate = y - residual / slope;
            y = if candidate.is_finite() && candidate > lower && candidate < upper && candidate != y
            {
                candidate
            } else {
                0.5 * lower + 0.5 * upper
            };
        }
        Err(CreepError::Invalid(
            "bounded primary creep rate root did not converge",
        ))
    }
}
