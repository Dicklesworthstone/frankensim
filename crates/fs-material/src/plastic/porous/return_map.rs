//! Shared pressure-sensitive return in an elastic strain / conjugate stress pair.
//! Kinematics and admissible history belong to each caller: infinitesimal strain
//! for PorousPlasticity, logarithmic elastic strain and Kirchhoff stress for
//! finite::porous. No strain-domain test is bypassed by either public law.

use super::{hypot, invalid, norm, positive};
use crate::elastic::bulk_shear_stiffness;
use crate::tensor::{deviator, trace};
use crate::{MaterialError, Tangent6, Voigt};
use fs_math::det;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(in crate::plastic) struct ReturnParameters {
    pub bulk: f64,
    pub shear: f64,
    pub c: f64,
    a2: f64,
    dev_return: f64,
    vol_return: f64,
    yield_stress: f64,
    hardening: f64,
    sensitivity: f64,
}

/// Derivative is with respect to the supplied elastic predictor strain.
pub(in crate::plastic) struct PressureReturn {
    pub stress: Voigt,
    pub tangent: Tangent6,
    pub elastic_strain: Voigt,
    pub plastic_increment: Voigt,
    pub alpha: f64,
    pub multiplier: f64,
    pub yield_residual: f64,
}

impl ReturnParameters {
    pub fn new(
        bulk: f64,
        shear: f64,
        yield_stress: f64,
        hardening: f64,
        sensitivity: f64,
    ) -> Result<Self, MaterialError> {
        if !positive(bulk)
            || !positive(shear)
            || !positive(yield_stress)
            || !(hardening.is_finite() && hardening >= 0.0)
            || !(sensitivity.is_finite() && sensitivity >= 0.0)
        {
            return Err(MaterialError::Parameters {
                what:
                    "porous bulk/shear elasticity, yield, hardening or sensitivity is inadmissible"
                        .into(),
            });
        }
        let a2 = sensitivity * sensitivity;
        let c = 1.0 + a2 / 9.0;
        let dev_return = 3.0 * (shear / c);
        let vol_return = bulk * (a2 / c);
        if !a2.is_finite()
            || !positive(c)
            || !positive(dev_return)
            || !vol_return.is_finite()
            || (sensitivity > 0.0 && !positive(vol_return))
            || !(hardening + dev_return).is_finite()
            || !(hardening + vol_return).is_finite()
            || bulk_shear_stiffness(bulk, shear)
                .iter()
                .flatten()
                .any(|v| !v.is_finite())
        {
            return Err(MaterialError::Parameters {
                what: "porous elastic or return moduli are unrepresentable".into(),
            });
        }
        Ok(Self {
            bulk,
            shear,
            c,
            a2,
            dev_return,
            vol_return,
            yield_stress,
            hardening,
            sensitivity,
        })
    }

    pub fn response(
        &self,
        elastic: &Voigt,
        alpha: f64,
        yield_roundoff: f64,
    ) -> Result<PressureReturn, MaterialError> {
        if elastic.iter().any(|x| !x.is_finite())
            || !(alpha.is_finite() && alpha >= 0.0)
            || !(yield_roundoff.is_finite() && yield_roundoff >= 0.0)
        {
            return Err(invalid(
                "return predictor, accumulation or numerical zero band is inadmissible",
            ));
        }
        let s_trial = deviator(elastic).map(|v| 2.0 * self.shear * v);
        let p_trial = self.bulk * trace(elastic);
        let q_trial = det::sqrt(1.5) * norm(&s_trial);
        let root_c = det::sqrt(self.c);
        let qn = q_trial / root_c;
        let pn = (self.sensitivity / root_c) * p_trial;
        let phi_trial = hypot(qn, pn);
        let y0 = self.yield_stress + self.hardening * alpha;
        if !phi_trial.is_finite() || !positive(y0) || s_trial.iter().any(|x| !x.is_finite()) {
            return Err(invalid("trial stress or yield strength is unrepresentable"));
        }
        if phi_trial - y0 <= yield_roundoff {
            return Ok(PressureReturn {
                stress: core::array::from_fn(|i| s_trial[i] + if i < 3 { p_trial } else { 0.0 }),
                tangent: bulk_shear_stiffness(self.bulk, self.shear),
                elastic_strain: *elastic,
                plastic_increment: [0.0; 6],
                alpha,
                multiplier: 0.0,
                yield_residual: phi_trial - y0,
            });
        }
        let (minimum, maximum) = match (qn != 0.0, pn != 0.0) {
            (true, true) => (
                self.dev_return.min(self.vol_return),
                self.dev_return.max(self.vol_return),
            ),
            (true, false) => (self.dev_return, self.dev_return),
            (false, true) => (self.vol_return, self.vol_return),
            _ => return Err(invalid("nonzero yielding requires a represented stress")),
        };
        let mut lower = (phi_trial - y0) / (self.hardening + maximum);
        let mut upper = (phi_trial - y0) / (self.hardening + minimum);
        if !lower.is_finite() || !positive(upper) {
            return Err(invalid("return bracket is unrepresentable"));
        }
        // Starting at the lower bound resolves the tiny-sensitivity J2 limit.
        let mut dl = lower;
        let mut solved = false;
        for _ in 0..96 {
            let dd = y0 + (self.hardening + self.dev_return) * dl;
            let dp = y0 + (self.hardening + self.vol_return) * dl;
            let (u, v) = (qn / dd, pn / dp);
            let value = hypot(u, v);
            if !(positive(dd) && positive(dp) && value.is_finite()) {
                return Err(invalid("scalar return left the finite set"));
            }
            let residual = value - 1.0;
            if residual.abs() <= 32.0 * f64::EPSILON {
                solved = true;
                break;
            }
            if residual > 0.0 {
                lower = dl;
            } else {
                upper = dl;
            }
            let derivative = -(u * u * ((self.hardening + self.dev_return) / dd)
                + v * v * ((self.hardening + self.vol_return) / dp))
                / value;
            let newton = dl - residual / derivative;
            dl = if newton.is_finite() && newton > lower && newton < upper {
                newton
            } else {
                0.5 * lower + 0.5 * upper
            };
        }
        if !solved {
            return Err(invalid("bounded scalar return did not converge"));
        }
        let y = y0 + self.hardening * dl;
        let dd = y0 + (self.hardening + self.dev_return) * dl;
        let dp = y0 + (self.hardening + self.vol_return) * dl;
        let theta_d = y / dd;
        let theta_p = y / dp;
        let sd = s_trial.map(|v| v / dd);
        let pd = p_trial / dp;
        let plastic_increment: Voigt = core::array::from_fn(|i| {
            dl * (1.5 / self.c * sd[i]
                + if i < 3 {
                    (self.a2 / self.c) * pd / 3.0
                } else {
                    0.0
                })
        });
        let next_alpha = alpha + dl;
        if !(next_alpha.is_finite() && next_alpha > alpha) {
            return Err(invalid(
                "plastic increment cannot advance the represented history",
            ));
        }
        let stress = core::array::from_fn(|i| {
            theta_d * s_trial[i] + if i < 3 { theta_p * p_trial } else { 0.0 }
        });
        // h=.5[(qn/Dd)^2+(pn/Dp)^2-1]; h_lambda=-denominator,
        // h_strain=v, stress_lambda=-Y0*v. Tensor-Voigt shear columns weigh two.
        let v: Voigt = core::array::from_fn(|i| {
            (self.dev_return / dd) * sd[i]
                + if i < 3 {
                    (self.vol_return / dp) * pd
                } else {
                    0.0
                }
        });
        let denominator = (qn / dd).powi(2) * ((self.hardening + self.dev_return) / dd)
            + (pn / dp).powi(2) * ((self.hardening + self.vol_return) / dp);
        if !positive(denominator) {
            return Err(invalid("consistent return derivative is singular"));
        }
        let mut tangent = bulk_shear_stiffness(self.bulk * theta_p, self.shear * theta_d);
        for (i, row) in tangent.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                *value -= (y0 / denominator) * v[i] * v[j] * if j < 3 { 1.0 } else { 2.0 };
            }
        }
        let elastic_deviator = deviator(elastic);
        let elastic_mean = trace(elastic) / 3.0;
        let elastic_strain = core::array::from_fn(|i| {
            theta_d * elastic_deviator[i] + if i < 3 { theta_p * elastic_mean } else { 0.0 }
        });
        let yield_residual =
            hypot(qn * theta_d, pn * theta_p) - (self.yield_stress + self.hardening * next_alpha);
        if stress
            .iter()
            .chain(tangent.iter().flatten())
            .chain(plastic_increment.iter())
            .any(|v| !v.is_finite())
            || !yield_residual.is_finite()
        {
            return Err(invalid(
                "corrected stress, tangent or flow is unrepresentable",
            ));
        }
        Ok(PressureReturn {
            stress,
            tangent,
            elastic_strain,
            plastic_increment,
            alpha: next_alpha,
            multiplier: dl,
            yield_residual,
        })
    }
}
