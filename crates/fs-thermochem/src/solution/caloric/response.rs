//! Local isobaric storage response from the original liquid caloric owner.

use super::{BinaryCaloricErrorV1, BinaryCaloricStateV1, R, Sum, checked, positive};
use crate::solution::BinaryMixingStabilityV1;
use fs_math::det::{exp, expm1, sqrt};

/// Whether a local storage response is a physical derivative or a declared
/// homogeneous-side tangent for a nonsmooth or unresolved phase boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BinaryCaloricHeatCapacityQualificationV1 {
    /// Fixed composition in one retained homogeneous phase.
    Homogeneous,
    /// Equilibrium redistribution between two distinct liquid phases.
    Coexisting,
    /// Exact binodal/critical state: select the homogeneous one-sided tangent.
    HomogeneousPhaseBoundaryTangent,
    /// Feed lies within the retained numerical binodal bracket. The selected
    /// homogeneous tangent is an approximate Newton derivative, not a claim
    /// that the physical state lies on the homogeneous side of that bracket.
    HomogeneousBinodalBracketTangent,
}

/// Positive total dH/dT [J/K], qualified for the retained phase state.
/// Boundary qualifications explicitly select a tangent suitable for an
/// implicit thermal solve; they do not assert a unique physical heat capacity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryCaloricHeatCapacityResponseV1 {
    heat_capacity_j_per_k: f64,
    qualification: BinaryCaloricHeatCapacityQualificationV1,
}

impl BinaryCaloricHeatCapacityResponseV1 {
    /// Actual smooth heat capacity or the explicitly qualified boundary tangent.
    pub const fn heat_capacity_j_per_k(&self) -> f64 {
        self.heat_capacity_j_per_k
    }

    pub const fn qualification(&self) -> BinaryCaloricHeatCapacityQualificationV1 {
        self.qualification
    }
}

impl BinaryCaloricStateV1 {
    /// Total isobaric heat capacity of this conserved inventory, dH/dT [J/K].
    ///
    /// An equilibrated state includes the change of coexistence compositions,
    /// in addition to each original pure-component Cp(T). An explicitly
    /// homogeneous state retains its fixed-composition constraint; it does not
    /// acquire an equilibrium response merely by requesting this derivative.
    /// This is an enthalpy derivative, not a constant-volume heat capacity.
    ///
    /// With temperature-dependent interactions, equilibrium phase-onset and
    /// critical points have different one-sided derivatives. Those points,
    /// and feeds whose branch is unresolved by the retained numerical binodal
    /// bracket, refuse instead of choosing a side. When the interaction's
    /// enthalpy coefficient vanishes at this temperature, its local phase
    /// redistribution contributes no heat capacity and needs no such refusal.
    /// The inherited phase-root accuracy limits the response; no new phase
    /// solve, numerical temperature differencing or allocation is performed.
    pub fn isobaric_heat_capacity_j_per_k(&self) -> Result<f64, BinaryCaloricErrorV1> {
        let response = self.isobaric_heat_capacity_response()?;
        match response.qualification {
            BinaryCaloricHeatCapacityQualificationV1::Homogeneous
            | BinaryCaloricHeatCapacityQualificationV1::Coexisting => {
                Ok(response.heat_capacity_j_per_k)
            }
            BinaryCaloricHeatCapacityQualificationV1::HomogeneousPhaseBoundaryTangent => {
                Err(BinaryCaloricErrorV1::Invalid(
                    "equilibrium heat capacity has distinct phase-onset or critical limits",
                ))
            }
            BinaryCaloricHeatCapacityQualificationV1::HomogeneousBinodalBracketTangent => {
                Err(BinaryCaloricErrorV1::Invalid(
                    "heat-capacity branch unresolved by the retained binodal bracket",
                ))
            }
        }
    }

    /// Qualified local isobaric response for an implicit thermal solver.
    ///
    /// Smooth homogeneous/coexisting states return their physical heat capacity.
    /// At a phase boundary this selects the homogeneous one-sided tangent;
    /// inside the retained binodal bracket it selects that same positive
    /// approximate Newton tangent. The qualification always identifies either
    /// selection. Only these two boundary cases use a one-sided response:
    /// unsupported sources and unrepresentable arithmetic still refuse.
    /// The actual energy residual must continue to use the original enthalpy.
    pub fn isobaric_heat_capacity_response(
        &self,
    ) -> Result<BinaryCaloricHeatCapacityResponseV1, BinaryCaloricErrorV1> {
        let mut capacity = Sum::default();
        for (moles, source) in self.component_moles.iter().zip(&self.owner.heat_capacity) {
            if *moles > 0.0 {
                let contribution = checked(
                    moles * source.heat_capacity_at(self.temperature)?.value(),
                    "pure-liquid heat capacity contribution",
                )?;
                if contribution <= 0.0 {
                    return Err(BinaryCaloricErrorV1::Arithmetic(
                        "positive pure-liquid heat capacity contribution underflows",
                    ));
                }
                capacity.add(contribution)?;
            }
        }
        let (redistribution, qualification) = self.equilibrium_excess_heat_capacity_j_per_k()?;
        let coefficient = self
            .owner
            .thermodynamics
            .excess_heat_capacity_coefficient_j_per_mol_k();
        if coefficient != 0.0 {
            let fractions = match qualification {
                BinaryCaloricHeatCapacityQualificationV1::HomogeneousPhaseBoundaryTangent
                | BinaryCaloricHeatCapacityQualificationV1::HomogeneousBinodalBracketTangent => {
                    self.equilibrium
                        .as_ref()
                        .ok_or(BinaryCaloricErrorV1::Invalid(
                            "qualified heat-capacity boundary requires an equilibrium state",
                        ))?
                        .feed()
                        .mole_fractions()
                }
                _ => self.mixing.mole_fractions(),
            };
            // Retain the actual phase composition (or the explicitly selected
            // homogeneous-side tangent). This term may have either sign;
            // model admission proves positive total fixed-composition Cp.
            let excess = checked(
                (coefficient * fractions[0].max(fractions[1]))
                    * (self.total_moles * fractions[0].min(fractions[1])),
                "fixed-composition excess heat capacity",
            )?;
            capacity.add(excess)?;
        }
        capacity.add(redistribution)?;
        let value = capacity.value()?;
        if !positive(value) {
            return Err(BinaryCaloricErrorV1::Arithmetic(
                "positive total liquid heat capacity is unrepresentable",
            ));
        }
        Ok(BinaryCaloricHeatCapacityResponseV1 {
            heat_capacity_j_per_k: value,
            qualification,
        })
    }

    fn equilibrium_excess_heat_capacity_j_per_k(
        &self,
    ) -> Result<(f64, BinaryCaloricHeatCapacityQualificationV1), BinaryCaloricErrorV1> {
        use BinaryCaloricHeatCapacityQualificationV1 as Q;
        let Some(equilibrium) = self.equilibrium.as_ref() else {
            return Ok((0.0, Q::Homogeneous));
        };
        if equilibrium.coexistence().is_none()
            && equilibrium.feed_stability() != BinaryMixingStabilityV1::Critical
        {
            return Ok((0.0, Q::Homogeneous));
        }
        let pair = &self.owner.thermodynamics;
        let legacy = pair.excess_heat_capacity_coefficient_j_per_mol_k() == 0.0;
        // -T² dchi/dT = H_interaction(T)/R. Keep the original B path for
        // zero excess Cp, including its established scaling/roundoff behavior.
        let interaction = if legacy {
            pair.chi_inverse_temperature_k
        } else {
            pair.interaction_enthalpy_j_per_mol(self.temperature)?
        };
        let split = equilibrium.phases().len() == 2;
        if interaction == 0.0 {
            return Ok((0.0, if split { Q::Coexisting } else { Q::Homogeneous }));
        }
        if matches!(
            equilibrium.feed_stability(),
            BinaryMixingStabilityV1::Critical | BinaryMixingStabilityV1::Binodal
        ) {
            return Ok((0.0, Q::HomogeneousPhaseBoundaryTangent));
        }
        let Some(coexistence) = equilibrium.coexistence() else {
            return Ok((0.0, Q::Homogeneous));
        };
        let [low, high] = coexistence.log_odds_bracket();
        // x(u)=exp(-u)/(1+exp(-u)) decreases monotonically. Retain the
        // original root's complete numerical bracket, expanded for endpoint
        // rounding, instead of treating its midpoint as an exact phase edge.
        // This is the same numerical qualification as the phase owner, not
        // an interval proof about the caller's source uncertainty.
        let lower_dilute = dilute_fraction(high).next_down();
        let upper_dilute = dilute_fraction(low).next_up();
        let feed = equilibrium.feed().mole_fractions();
        let dilute_feed = feed[0].min(feed[1]);
        let branch_resolved = if split {
            dilute_feed > upper_dilute
        } else {
            dilute_feed < lower_dilute
        };
        if !branch_resolved {
            return Ok((0.0, Q::HomogeneousBinodalBracketTangent));
        }
        if !split {
            return Ok((0.0, Q::Homogeneous));
        }

        // The original binodal satisfies ln(x/(1-x))+chi(1-2x)=0.
        // Differentiating it at fixed total inventory gives
        // Cp_redistribution = n (H_interaction/T)^2/R *
        //                    (1-2x)^2 x(1-x) / Gamma,
        // Gamma=1-2 chi x(1-x)>0. Both numerator and Gamma vanish at the
        // critical point; evaluate their finite ratio from the same accepted
        // log-odds root rather than dividing rounded, cancelling curvatures.
        let root = 0.5 * low + 0.5 * high;
        let geometry = coexistence_heat_capacity_factor(root);
        // Scale before squaring: the intensive factor can underflow when total
        // capacity of the supplied inventory remains representable.
        let amplitude = (sqrt(self.total_moles) * (interaction.abs() / self.temperature.value()))
            * sqrt(if legacy { R * geometry } else { geometry / R });
        let excess = checked(amplitude * amplitude, "equilibrium excess heat capacity")?;
        if !positive(excess) {
            return Err(BinaryCaloricErrorV1::Arithmetic(
                "positive equilibrium excess heat capacity is unrepresentable",
            ));
        }
        Ok((excess, Q::Coexisting))
    }
}

fn dilute_fraction(log_odds: f64) -> f64 {
    let e = exp(-log_odds);
    e / (1.0 + e)
}

// For u=ln((1-x)/x), the binodal relation gives
// (1-2x)^2 x(1-x)/Gamma = tanh(u/2)^3 / [2(sinh(u)-u)].
// This ratio tends to 3/8 on the split side. Factoring u^3 avoids loss
// of significance without replacing any finite source state by that limit.
fn coexistence_heat_capacity_factor(u: f64) -> f64 {
    let one_minus_e = -expm1(-u);
    if u <= 1.0 {
        let u2 = u * u;
        // (sinh(u)-u)/u^3. All terms are positive; at u<=1 the
        // omitted 1/25! term is below double-precision roundoff.
        let sinh_remainder = [
            1.0 / 6.0,
            1.0 / 120.0,
            1.0 / 5040.0,
            1.0 / 362880.0,
            1.0 / 39916800.0,
            1.0 / 6227020800.0,
            1.0 / 1307674368000.0,
            1.0 / 355687428096000.0,
            1.0 / 121645100408832000.0,
            1.0 / 51090942171709440000.0,
            1.0 / 25852016738884976640000.0,
        ]
        .into_iter()
        .rev()
        .fold(0.0, |sum, coefficient| coefficient + u2 * sum);
        let scaled_tanh = (one_minus_e / u) / (2.0 - one_minus_e);
        (scaled_tanh * scaled_tanh * scaled_tanh) / (2.0 * sinh_remainder)
    } else {
        let e = exp(-u);
        let one_plus_e = 1.0 + e;
        let difference = one_minus_e / one_plus_e;
        let product = (e / one_plus_e) / one_plus_e;
        let curvature = 1.0 - (2.0 * u * e) / (one_minus_e * one_plus_e);
        (difference * difference) * (product / curvature)
    }
}

#[cfg(test)]
mod tests;
