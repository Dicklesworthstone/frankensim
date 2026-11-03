//! Source-bound, signed reversible cell potentials for an explicit reaction.
//!
//! IUPAC defines E° = -Delta_r G°/(n F):
//! <https://goldbook.iupac.org/terms/view/S05913>.
//! The reversible potential is the zero-current thermodynamic limit:
//! <https://goldbook.iupac.org/terms/view/E01934>.
//!
//! The sign belongs to the written reaction. A negative value for copper
//! electrowinning describes its externally driven productive direction. The
//! opposing reversible voltage magnitude does not predict the operating cell
//! voltage, overpotential, resistance, kinetics or current efficiency. No
//! efficiency factor is applied to this equilibrium thermodynamic quantity.

use fs_qty::{Dimensionless, Voltage};

use super::{ElectrolysisReaction, Result};
use crate::{
    MolarEnergyQuantityV1, SpeciesId,
    stoichiometry::{
        FARADAY_C_PER_MOL, ReactionBalanceError, product, representable,
        thermodynamics::{FormationGibbsSource, ReactionActivityEvaluation, StandardReactionGibbs},
    },
};

/// Reversible standard potential retaining its original equation, electron
/// association and complete source-backed Gibbs result.
#[derive(Clone, Debug)]
pub struct StandardCellPotential {
    model: ElectrolysisReaction,
    gibbs: StandardReactionGibbs,
    molar_transferred_charge: f64,
    potential: Voltage,
}

impl ElectrolysisReaction {
    /// Evaluate the original reaction's formation-Gibbs sources, then divide
    /// its signed standard Gibbs energy by the original equation's n F.
    ///
    /// A Gibbs result for another equation cannot be supplied. The declared
    /// electron association remains a caller input, just as in Faraday yield;
    /// atom/charge balance does not authenticate an electrode mechanism.
    /// Electron counts that cannot convert exactly to binary64 refuse.
    pub fn standard_potential(
        &self,
        sources: &[FormationGibbsSource<'_>],
    ) -> Result<StandardCellPotential> {
        let charge = molar_transferred_charge(self.electrons_per_equation())?;
        let gibbs = self.reaction().standard_gibbs_change(sources)?;
        let potential = signed_potential(gibbs.molar(), charge)?;
        Ok(StandardCellPotential {
            model: self.clone(),
            gibbs,
            molar_transferred_charge: charge,
            potential,
        })
    }
}

impl StandardCellPotential {
    /// Original conserved equation and declared Faraday product association.
    #[must_use]
    pub const fn model(&self) -> &ElectrolysisReaction {
        &self.model
    }

    /// Original standard Gibbs result, including each selected source receipt.
    #[must_use]
    pub const fn gibbs(&self) -> &StandardReactionGibbs {
        &self.gibbs
    }

    /// Signed E° = -Delta_r G°/(n F) [V], for the equation as written.
    #[must_use]
    pub const fn potential(&self) -> Voltage {
        self.potential
    }

    /// Evaluate E = -Delta_r G/(n F) using the retained source conventions
    /// and one finite positive dimensionless activity per active species.
    ///
    /// This delegates activity admission and ln Q to the original Gibbs
    /// consumer. It retains the local dG/dxi result and does not integrate a
    /// finite extent or infer an operating voltage at nonzero current.
    pub fn at_activities(
        &self,
        activities: &[(SpeciesId, Dimensionless)],
    ) -> Result<ActivityCellPotential<'_>> {
        let gibbs = self.gibbs.at_activities(activities)?;
        let potential = signed_potential(gibbs.molar, self.molar_transferred_charge)?;
        Ok(ActivityCellPotential {
            standard: self,
            gibbs,
            potential,
        })
    }
}

/// Reversible potential at supplied activities, with the original activity
/// evaluation and source/electron basis retained without alteration.
#[derive(Debug)]
pub struct ActivityCellPotential<'a> {
    standard: &'a StandardCellPotential,
    gibbs: ReactionActivityEvaluation<'a>,
    potential: Voltage,
}

impl<'a> ActivityCellPotential<'a> {
    /// Source-backed standard state and original Faraday association.
    #[must_use]
    pub const fn standard(&self) -> &'a StandardCellPotential {
        self.standard
    }

    /// Original local Gibbs evaluation and all explicitly supplied activities.
    #[must_use]
    pub const fn gibbs(&self) -> &ReactionActivityEvaluation<'a> {
        &self.gibbs
    }

    /// Signed E = -Delta_r G/(n F) [V] at the supplied activity state.
    #[must_use]
    pub const fn potential(&self) -> Voltage {
        self.potential
    }
}

fn molar_transferred_charge(electrons: u128) -> Result<f64> {
    if electrons == 0 {
        return Err(ReactionBalanceError::Invalid(
            "electron count per equation must be positive",
        ));
    }
    // Removing powers of two leaves the exact significand width. This also
    // handles integers above 2^53 without a saturating float-to-int round trip.
    let significand_bits = u128::BITS - electrons.leading_zeros() - electrons.trailing_zeros();
    if significand_bits > f64::MANTISSA_DIGITS {
        return Err(ReactionBalanceError::Invalid(
            "electron count per equation is not exactly representable in binary64",
        ));
    }
    product(
        electrons as f64,
        FARADAY_C_PER_MOL,
        "molar transferred charge",
    )
}

fn signed_potential(gibbs: MolarEnergyQuantityV1, charge: f64) -> Result<Voltage> {
    let value = -gibbs.value() / charge;
    representable(value, gibbs.value() != 0.0, "reversible cell potential")?;
    Ok(Voltage::new(value))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_electron_conversion_accepts_large_powers_and_refuses_rounded_counts() {
        for n in [1, 4, 1_u128 << 53, (1_u128 << 53) + 2, 1_u128 << 84] {
            assert!(molar_transferred_charge(n).unwrap().is_finite());
        }
        for n in [0, (1_u128 << 53) + 1, (1_u128 << 54) + 2, u128::MAX] {
            assert!(matches!(
                molar_transferred_charge(n),
                Err(ReactionBalanceError::Invalid(_))
            ));
        }
    }

    #[test]
    fn potential_preserves_reference_zero_and_refuses_nonzero_underflow() {
        let charge = molar_transferred_charge(4).unwrap();
        assert_eq!(
            signed_potential(MolarEnergyQuantityV1::new(0.0), charge)
                .unwrap()
                .value(),
            0.0
        );
        assert!(matches!(
            signed_potential(MolarEnergyQuantityV1::new(f64::from_bits(1)), charge),
            Err(ReactionBalanceError::Arithmetic(
                "reversible cell potential"
            ))
        ));
    }
}
