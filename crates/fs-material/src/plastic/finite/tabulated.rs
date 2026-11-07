//! Tabulated energetic hardening on the original multiplicative J2 model.
//!
//! The table contains absolute **Kirchhoff** yield stress versus accumulated
//! equivalent plastic strain. A Cauchy or engineering tensile curve is not
//! silently relabelled as this finite constitutive law. Every history remains
//! inside both the supplied curve and the original finite-strain limits.
//!
//! Piecewise-linear isotropic hardening is the standard radial-return extension
//! described by the [MOOSE isotropic plasticity model](https://mooseframework.inl.gov/moose/source/materials/IsotropicPlasticityStressUpdate.html).
//! Here the existing finite tensor return, plastic metric, matrix-log derivative
//! and work convention remain unchanged. Only the scalar consistency relation
//! changes. At a table kink the right-hand slope is selected (the final segment
//! at the last knot); a classical two-sided tangent does not exist there.
//!
//! Hardening is explicitly energetic: the integral of `Y(alpha)-Y(0)` is stored,
//! while `Y(0)*delta_alpha` is physical yield loss. An empirical flow curve alone
//! does not identify that partition for a real specimen. Perzyna kinetics adds
//! its original viscous work and potential without changing the table's domain.

/// Original material-card Kirchhoff curves on the native finite owner.
pub mod source;

use super::rate::FiniteOverstressJ2;
use super::{ElasticEndpoint, FiniteJ2, FiniteJ2Error, FiniteJ2State, RadialReturn, ReturnMap};
use crate::plastic::tabulated::TabulatedHardening;

impl FiniteJ2 {
    /// Select an explicit bounded Kirchhoff-yield curve. The original initial
    /// yield must equal the first ordinate and the original linear modulus must
    /// be zero: the table is the complete isotropic resistance, not an extra
    /// term. Its support must cover the entire declared accumulated-strain
    /// domain. Construct this complete law before retaining any history.
    pub fn with_tabulated_hardening(
        mut self,
        table: TabulatedHardening,
    ) -> Result<Self, FiniteJ2Error> {
        if self.voce.is_some()
            || self.parameters.hardening_modulus_pa != 0.0
            || self.parameters.yield_stress_pa != table.initial_yield_stress_pa()
            || self.limits.maximum_accumulated_plastic_strain > table.maximum_plastic_strain()
            || self.table.as_ref().is_some_and(|old| old != &table)
        {
            return Err(FiniteJ2Error::Parameters(
                "finite table requires matching initial yield, zero linear hardening, full history support and no replacement or Voce selection",
            ));
        }
        let elastic = 3.0 * self.parameters.shear_modulus_pa;
        for point in table.points() {
            let slope = table
                .tangent_pa(point.plastic_strain)
                .map_err(FiniteJ2Error::Tabulated)?;
            if !(elastic + slope).is_finite() {
                return Err(FiniteJ2Error::Parameters(
                    "tabulated finite return modulus must be representable on every segment",
                ));
            }
        }
        self.table = Some(table);
        Ok(self)
    }

    /// The exact immutable table, if explicitly selected. None retains the
    /// original linear/Voce law and its original arithmetic.
    pub fn tabulated_hardening(&self) -> Option<&TabulatedHardening> {
        self.table.as_ref()
    }

    pub(super) fn return_map_tabulated(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
        trial: ElasticEndpoint,
        q_trial: f64,
        kinetics: Option<&FiniteOverstressJ2>,
    ) -> Result<ReturnMap, FiniteJ2Error> {
        let table = self.table.as_ref().expect("selected tabulated return");
        let radius_old = table
            .yield_stress_pa(old.alpha)
            .map_err(FiniteJ2Error::Tabulated)?;
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
            let increment = match kinetics {
                Some(law) => self.tabulated_rate_increment(excess, old.alpha, law)?,
                None => table
                    .solve_increment(q_trial, old.alpha, elastic)
                    .map_err(FiniteJ2Error::Tabulated)?,
            };
            let alpha = old.alpha + increment;
            let radius = table
                .yield_stress_pa(alpha)
                .map_err(FiniteJ2Error::Tabulated)?;
            let hardening = table.tangent_pa(alpha).map_err(FiniteJ2Error::Tabulated)?;
            let strength_increment = table
                .strength_increment_pa(old.alpha, increment)
                .map_err(FiniteJ2Error::Tabulated)?;
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
            let denominator = elastic + hardening + viscous_tangent;
            if !denominator.is_finite() {
                return Err(FiniteJ2Error::Numerical(
                    "tabulated finite derivative is not representable",
                ));
            }
            let viscous_work = overstress * increment;
            radial.plastic_increment = increment;
            radial.scale = (radius + overstress) / q_trial;
            // dq_return/dq_trial - q_return/q_trial. The increment form avoids
            // subtracting two nearly equal radial moduli on an affine segment.
            radial.tangent_correction = -(elastic / denominator)
                * ((radius_old + (strength_increment - hardening * increment) + transverse)
                    / q_trial);
            radial.primary_potential = self.parameters.yield_stress_pa * increment;
            radial.plastic_dissipation = radial.primary_potential + viscous_work;
            radial.rate_potential = potential_fraction * viscous_work;
        }
        self.finish_return(f, old, trial, q_trial, radial)
    }

    /// Safeguarded scalar solve including the original Perzyna inverse rate.
    /// All curve evaluations are inside the declared domain, including bracket
    /// trials. A root beyond the last supplied point refuses, never extrapolates.
    fn tabulated_rate_increment(
        &self,
        excess: f64,
        alpha: f64,
        kinetics: &FiniteOverstressJ2,
    ) -> Result<f64, FiniteJ2Error> {
        let table = self.table.as_ref().expect("selected tabulated return");
        let elastic = 3.0 * self.parameters.shear_modulus_pa;
        let maximum = self.limits.maximum_accumulated_plastic_strain;
        let remaining = maximum - alpha;
        let unconstrained_upper = kinetics.solve(excess, elastic)?.increment;
        if remaining <= 0.0 {
            return Err(FiniteJ2Error::Domain {
                quantity: "tabulated accumulated plastic strain",
                value: alpha + unconstrained_upper,
                limit: maximum,
            });
        }
        let residual = |delta: f64| -> Result<(f64, f64), FiniteJ2Error> {
            let increment = table
                .strength_increment_pa(alpha, delta)
                .map_err(FiniteJ2Error::Tabulated)?;
            let slope = table
                .tangent_pa(alpha + delta)
                .map_err(FiniteJ2Error::Tabulated)?;
            let viscous = kinetics.endpoint_overstress_pa(delta)?;
            let derivative = elastic + slope + (viscous / delta) / kinetics.parameters().exponent;
            Ok((
                elastic.mul_add(delta, increment) + viscous - excess,
                derivative,
            ))
        };
        let tolerance = 64.0 * f64::EPSILON * excess;
        let mut upper = unconstrained_upper.min(remaining);
        let (upper_residual, _) = residual(upper)?;
        if upper_residual < -tolerance {
            return Err(FiniteJ2Error::Domain {
                quantity: "tabulated accumulated plastic strain",
                value: alpha + unconstrained_upper,
                limit: maximum,
            });
        }
        if upper_residual.abs() <= tolerance {
            return Ok(upper);
        }
        let slope = table.tangent_pa(alpha).map_err(FiniteJ2Error::Tabulated)?;
        let mut delta = kinetics
            .solve(excess, elastic + slope)?
            .increment
            .min(upper);
        let mut lower = 0.0;
        for _ in 0..96 {
            let (value, derivative) = residual(delta)?;
            if value.is_finite() && value.abs() <= tolerance {
                return Ok(delta);
            }
            if !value.is_finite() || !derivative.is_finite() || derivative <= 0.0 {
                break;
            }
            if value > 0.0 {
                upper = delta;
            } else {
                lower = delta;
            }
            let middle = 0.5 * lower + 0.5 * upper;
            if middle == lower || middle == upper {
                break;
            }
            let newton = delta - value / derivative;
            delta = if newton.is_finite() && newton > lower && newton < upper {
                newton
            } else {
                middle
            };
        }
        Err(FiniteJ2Error::Numerical(
            "bounded tabulated finite rate return did not converge",
        ))
    }
}
