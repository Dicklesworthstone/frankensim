//! Energetic isotropic saturation in the original multiplicative Hencky law.
//!
//! The Kirchhoff yield radius is `Y + H alpha + Q (1 - exp(-b alpha))`.
//! The [Voce hardening family](https://neml.readthedocs.io/en/stable/hardening/simple/iso_voce.html)
//! uses one accumulated history; the existing plastic metric still owns every
//! tensor component. Q is an added saturation strength, not the final yield
//! stress. No small-strain material fit is silently promoted to this law.
//!
//! The additional recoverable energy is `Q [alpha + expm1(-b alpha)/b]`.
//! Consequently yield heat remains `Y delta_alpha`; neither hardening storage
//! nor the endpoint integration defect is heat. The same bracketed scalar
//! solve includes explicitly selected Perzyna overstress. Its endpoint
//! derivative drives the existing full dP/dF tensor tangent. The zero-Q case
//! retains the original linear law, including its arithmetic and fast path.

use super::rate::FiniteOverstressJ2;
use super::{ElasticEndpoint, FiniteJ2, FiniteJ2Error, FiniteJ2State, RadialReturn, ReturnMap};
use fs_math::det;

/// Explicit, immutable Voce coefficients at one temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteVoceHardening {
    saturation_pa: f64,
    rate_per_strain: f64,
}

impl FiniteVoceHardening {
    /// Added asymptotic isotropic strength Q [Pa].
    pub const fn saturation_pa(self) -> f64 {
        self.saturation_pa
    }

    /// Saturation rate b per unit equivalent plastic strain.
    pub const fn rate_per_strain(self) -> f64 {
        self.rate_per_strain
    }

    pub(super) fn strength(self, alpha: f64) -> f64 {
        self.saturation_pa * -det::expm1(-self.rate_per_strain * alpha)
    }

    pub(super) fn tangent(self, alpha: f64) -> f64 {
        (self.saturation_pa * self.rate_per_strain) * det::exp(-self.rate_per_strain * alpha)
    }

    pub(super) fn strength_increment(self, alpha: f64, increment: f64) -> f64 {
        (self.saturation_pa * det::exp(-self.rate_per_strain * alpha))
            * -det::expm1(-self.rate_per_strain * increment)
    }

    pub(super) fn stored_energy(self, alpha: f64) -> f64 {
        self.saturation_pa * saturation_primitive(alpha, self.rate_per_strain)
    }

    // Delta R - R'(alpha + delta) delta, a nonnegative concavity remainder.
    // This stable expression avoids a near-cancellation in the radial tangent.
    pub(super) fn curvature_remainder(self, alpha: f64, increment: f64) -> f64 {
        let x = self.rate_per_strain * increment;
        let shape = if x < 0.1 {
            let polynomial = -1.0 / 45360.0 + x / 403200.0;
            let polynomial = 1.0 / 5760.0 + x * polynomial;
            let polynomial = -1.0 / 840.0 + x * polynomial;
            let polynomial = 1.0 / 144.0 + x * polynomial;
            let polynomial = -1.0 / 30.0 + x * polynomial;
            let polynomial = 1.0 / 8.0 + x * polynomial;
            let polynomial = -1.0 / 3.0 + x * polynomial;
            x * x * (0.5 + x * polynomial)
        } else if x.is_infinite() {
            1.0
        } else {
            -det::expm1(-x) - x * det::exp(-x)
        };
        (self.saturation_pa * det::exp(-self.rate_per_strain * alpha)) * shape
    }
}

/// Integral of `1-exp(-b alpha)`, without loss of its quadratic small-alpha
/// limit. Shared with the affine-temperature entropy and coefficient accounts.
pub(super) fn saturation_primitive(alpha: f64, rate: f64) -> f64 {
    let x = rate * alpha;
    let shape = if x < 0.1 {
        let polynomial = -1.0 / 362880.0 + x / 3628800.0;
        let polynomial = 1.0 / 40320.0 + x * polynomial;
        let polynomial = -1.0 / 5040.0 + x * polynomial;
        let polynomial = 1.0 / 720.0 + x * polynomial;
        let polynomial = -1.0 / 120.0 + x * polynomial;
        let polynomial = 1.0 / 24.0 + x * polynomial;
        let polynomial = -1.0 / 6.0 + x * polynomial;
        x * (0.5 + x * polynomial)
    } else {
        1.0 + det::expm1(-x) / x
    };
    alpha * shape
}

impl FiniteJ2 {
    /// Add explicit energetic Voce saturation to the linear isotropic law.
    /// Q >= 0, b > 0, all coefficients and the full admitted strength, energy
    /// and return modulus must be representable. This does not add kinematic
    /// hardening, softening, damage, anisotropy or an empirical calibration.
    /// Construct the complete law before creating a retained material history.
    pub fn with_voce_hardening(
        mut self,
        saturation_pa: f64,
        rate_per_strain: f64,
    ) -> Result<Self, FiniteJ2Error> {
        let voce = FiniteVoceHardening {
            saturation_pa,
            rate_per_strain,
        };
        let slope = saturation_pa * rate_per_strain;
        if !saturation_pa.is_finite()
            || saturation_pa < 0.0
            || !rate_per_strain.is_finite()
            || rate_per_strain <= 0.0
            || !slope.is_finite()
            || (saturation_pa > 0.0 && slope == 0.0)
            || self.parameters.yield_stress_pa <= 0.0
            || !(3.0 * self.parameters.shear_modulus_pa
                + self.parameters.hardening_modulus_pa
                + slope)
                .is_finite()
        {
            return Err(FiniteJ2Error::Parameters(
                "finite Voce requires positive yield, finite Q >= 0 and b > 0, and representable slope/return modulus",
            ));
        }
        let selection = (saturation_pa != 0.0).then_some(voce);
        if self.voce.is_some() && self.voce != selection {
            return Err(FiniteJ2Error::Parameters(
                "finite Voce coefficients cannot be silently replaced",
            ));
        }
        self.voce = selection;
        let maximum = self.limits.maximum_accumulated_plastic_strain;
        if !self.yield_radius(maximum).is_finite() || !self.hardening_energy(maximum).is_finite() {
            return Err(FiniteJ2Error::Parameters(
                "finite Voce strength and storage must be representable throughout the history domain",
            ));
        }
        Ok(self)
    }

    /// Explicit saturation selection; `None` is exactly the original linear law.
    pub const fn voce_hardening(&self) -> Option<FiniteVoceHardening> {
        self.voce
    }

    pub(super) fn hardening_tangent(&self, alpha: f64) -> f64 {
        self.parameters.hardening_modulus_pa + self.voce.map_or(0.0, |v| v.tangent(alpha))
    }

    pub(super) fn return_map_voce(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
        trial: ElasticEndpoint,
        q_trial: f64,
        kinetics: Option<&FiniteOverstressJ2>,
    ) -> Result<ReturnMap, FiniteJ2Error> {
        let voce = self.voce.expect("saturated return has admitted coefficients");
        let radius_old = self.yield_radius(old.alpha);
        let elastic = 3.0 * self.parameters.shear_modulus_pa;
        let excess = q_trial - radius_old;
        let roundoff = 64.0
            * f64::EPSILON
            * (2.0 * self.parameters.shear_modulus_pa)
                .max(q_trial)
                .max(radius_old);
        let mut radial = RadialReturn {
            plastic_increment: 0.0,
            creep_increment: 0.0,
            scale: 1.0,
            tangent_correction: 0.0,
            plastic_dissipation: 0.0,
            creep_dissipation: 0.0,
            primary_potential: 0.0,
            rate_potential: 0.0,
        };
        if excess > roundoff {
            let increment = self.saturated_increment(excess, old.alpha, kinetics)?;
            let alpha = old.alpha + increment;
            let (overstress, viscous_tangent, transverse, potential_fraction) = match kinetics {
                Some(law) => {
                    let stress = law.endpoint_overstress_pa(increment)?;
                    let inverse_n = 1.0 / law.parameters().exponent;
                    (
                        stress,
                        (stress / increment) * inverse_n,
                        (1.0 - inverse_n) * stress,
                        1.0 / (1.0 + inverse_n),
                    )
                }
                None => (0.0, 0.0, 0.0, 0.0),
            };
            let denominator = elastic + self.hardening_tangent(alpha) + viscous_tangent;
            if !denominator.is_finite() || !viscous_tangent.is_finite() {
                return Err(FiniteJ2Error::Numerical(
                    "finite Voce derivative is not representable",
                ));
            }
            let viscous_work = overstress * increment;
            radial.plastic_increment = increment;
            radial.scale = (self.yield_radius(alpha) + overstress) / q_trial;
            radial.tangent_correction = -(elastic / denominator)
                * ((radius_old + voce.curvature_remainder(old.alpha, increment) + transverse)
                    / q_trial);
            radial.primary_potential = self.parameters.yield_stress_pa * increment;
            radial.plastic_dissipation = radial.primary_potential + viscous_work;
            radial.rate_potential = potential_fraction * viscous_work;
        }
        self.finish_return(f, old, trial, q_trial, radial)
    }

    // The hardening increment is evaluated directly, not by subtracting two
    // large endpoint yield radii. Bracket with the exact linear-rate owner;
    // Newton then differentiates saturation and kinetics together.
    pub(super) fn saturated_increment(
        &self,
        excess: f64,
        alpha: f64,
        kinetics: Option<&FiniteOverstressJ2>,
    ) -> Result<f64, FiniteJ2Error> {
        let voce = self.voce.expect("saturated return has admitted coefficients");
        let a = 3.0 * self.parameters.shear_modulus_pa + self.parameters.hardening_modulus_pa;
        let solve_linear = |modulus| match kinetics {
            Some(law) => law.solve(excess, modulus).map(|flow| flow.increment),
            None => Ok(excess / modulus),
        };
        let mut lower = 0.0;
        let mut upper = solve_linear(a)?;
        let mut increment = solve_linear(a + voce.tangent(alpha))?;
        for _ in 0..96 {
            if !increment.is_finite() || increment <= 0.0 {
                break;
            }
            let viscous = kinetics
                .map(|law| law.endpoint_overstress_pa(increment))
                .transpose()?
                .unwrap_or(0.0);
            let viscous_tangent =
                kinetics.map_or(0.0, |law| (viscous / increment) / law.parameters().exponent);
            let residual =
                a.mul_add(increment, voce.strength_increment(alpha, increment)) + viscous - excess;
            if residual.is_finite() && residual.abs() <= 64.0 * f64::EPSILON * excess {
                return Ok(increment);
            }
            if !residual.is_finite() || !viscous_tangent.is_finite() {
                break;
            }
            if residual > 0.0 {
                upper = increment;
            } else {
                lower = increment;
            }
            let derivative = a + voce.tangent(alpha + increment) + viscous_tangent;
            let newton = increment - residual / derivative;
            let middle = 0.5 * lower + 0.5 * upper;
            if middle == lower || middle == upper {
                break;
            }
            increment = if newton.is_finite() && newton > lower && newton < upper {
                newton
            } else {
                middle
            };
        }
        Err(FiniteJ2Error::Numerical(
            "bounded finite Voce return did not converge",
        ))
    }
}
