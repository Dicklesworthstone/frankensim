//! Relative entropy for the same parent/crystal calorics used by kinetics.
//!
//! Caloric branches alone do not determine their relative entropy datum. This
//! explicit model supplies one: the caller admits coexistence at Tm, so
//! `s_parent(Tm)-s_crystal(Tm) = [h_parent(Tm)-h_crystal(Tm)]/Tm`.
//! Both phases then retain their ORIGINAL integrals of `dh/T`. At fixed
//! pressure/composition the mixture uses crystalline MASS fraction X:
//! `s(T,X) = (1-X)*s_parent(T) + X*s_crystal(T)`. Separate macroscopic phases
//! acquire no invented ideal-mixing term, apparent Cp or extra latent heater.
//! Consequently `T ds = dh + (g_parent-g_crystal) dX`: an adiabatic crystal
//! increment produces entropy only where the original phase driving gap permits it.
//!
//! Tm is an admitted thermodynamic input, not a value inferred from a kinetic
//! fit or a material name. This model is unsuitable when the supplied branches
//! do not describe phases that can coexist at that temperature. It does not
//! establish stability over the rest of the domain; the reported Gibbs gap
//! exposes that question without altering the kinetic rate or accepted state.
//! Elastic, interfacial, chemical and variable-pressure entropy remain outside
//! this caloric model. A signed material entropy change is not by itself an
//! entropy-production certificate: external heat needs its actual boundary
//! entropy, and all omitted work/internal variables need their own accounting.

use super::{CrystallizationCaloricBranch, CrystallizationCaloricPair, invalid};
use crate::phase::PhaseStateError;
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::crystallization::{
    CrystallizationError, CrystallizationLaw, CrystallizationSolve, CrystallizationState,
    CrystallizationStep,
};

#[derive(Clone, Debug, PartialEq)]
pub enum CrystallizationEntropyError {
    Phase(PhaseStateError),
    Kinetics(CrystallizationError),
    /// Original branch identity, source domain or coexistence temperature differs.
    LawMismatch,
}
impl From<PhaseStateError> for CrystallizationEntropyError {
    fn from(error: PhaseStateError) -> Self {
        Self::Phase(error)
    }
}
impl From<CrystallizationError> for CrystallizationEntropyError {
    fn from(error: CrystallizationError) -> Self {
        Self::Kinetics(error)
    }
}
impl core::fmt::Display for CrystallizationEntropyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "crystallization entropy: {self:?}")
    }
}
impl std::error::Error for CrystallizationEntropyError {}

/// Explicit coexistence reference for two immutable original material charts.
///
/// Preparation caches the coexistence calorics; queries allocate nothing and
/// reuse the chart owner's exact segment integrals. Fixed-fraction changes
/// visit only crossed temperature segments, without revisiting Tm. Fraction
/// changes additionally integrate the original branches between Tm and the
/// warmer endpoint. No numerical quadrature or replacement Cp law is added.
#[derive(Clone, Copy, Debug)]
pub struct CrystallizationEntropy<'a> {
    pair: CrystallizationCaloricPair<'a>,
    coexistence_temperature_k: f64,
    parent_h_at_coexistence_j_kg: f64,
    crystal_h_at_coexistence_j_kg: f64,
    entropy_gap_at_coexistence_j_kg_k: f64,
    identity: ContentHash,
}

impl<'a> CrystallizationEntropy<'a> {
    /// Declare coexistence at Tm inside the original pair's source domain.
    /// The caller owns the applicability of this phase-equilibrium assumption.
    pub fn new(
        pair: &CrystallizationCaloricPair<'a>,
        coexistence_temperature_k: f64,
    ) -> Result<Self, CrystallizationEntropyError> {
        pair.admit_temperature(coexistence_temperature_k)?;
        let parent_h = pair
            .parent()
            .enthalpy_at_temperature(coexistence_temperature_k)?;
        let crystal_h = pair
            .crystal()
            .enthalpy_at_temperature(coexistence_temperature_k)?;
        let entropy_gap = finite((parent_h - crystal_h) / coexistence_temperature_k)?;
        if entropy_gap <= 0.0 {
            return Err(
                invalid("coexistence entropy gap must be positive and representable").into(),
            );
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-material.crystallization-entropy.v1");
        hash.update(pair.identity().as_bytes());
        hash.update(&coexistence_temperature_k.to_bits().to_le_bytes());
        Ok(Self {
            pair: *pair,
            coexistence_temperature_k,
            parent_h_at_coexistence_j_kg: parent_h,
            crystal_h_at_coexistence_j_kg: crystal_h,
            entropy_gap_at_coexistence_j_kg_k: entropy_gap,
            identity: hash.finalize(),
        })
    }

    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn caloric_pair(&self) -> &CrystallizationCaloricPair<'a> {
        &self.pair
    }
    pub const fn coexistence_temperature_k(&self) -> f64 {
        self.coexistence_temperature_k
    }

    /// `s_parent(T)-s_crystal(T)`, J/(kg K), under the declared reference.
    pub fn phase_entropy_gap_j_kg_k(
        &self,
        temperature_k: f64,
    ) -> Result<f64, CrystallizationEntropyError> {
        let (parent, crystal) = self.entropy_from_coexistence(temperature_k)?;
        sum([self.entropy_gap_at_coexistence_j_kg_k, parent, -crystal])
    }

    /// `g_parent(T)-g_crystal(T)`, J/kg, at fixed pressure/composition.
    /// Positive values favor conversion to crystal; negative values favor the
    /// parent. This is the declared two-branch potential, not a global stability
    /// proof or a newly calibrated kinetic law. At Tm the result is exactly zero.
    pub fn phase_gibbs_gap_j_kg(
        &self,
        temperature_k: f64,
    ) -> Result<f64, CrystallizationEntropyError> {
        let (parent_s, crystal_s) = self.entropy_from_coexistence(temperature_k)?;
        let parent_h = self.pair.parent().enthalpy_at_temperature(temperature_k)?;
        let crystal_h = self.pair.crystal().enthalpy_at_temperature(temperature_k)?;
        // Integrate from coexistence instead of subtracting two large absolute
        // Gibbs energies. The first term also preserves small T-Tm changes.
        sum([
            (self.coexistence_temperature_k - temperature_k)
                * self.entropy_gap_at_coexistence_j_kg_k,
            parent_h - self.parent_h_at_coexistence_j_kg,
            -(crystal_h - self.crystal_h_at_coexistence_j_kg),
            -temperature_k * (parent_s - crystal_s),
        ])
    }

    /// Signed mixture entropy change, J/(kg K), for explicit T/X endpoints.
    ///
    /// X is crystalline mass fraction, not equilibrium liquid or volume fraction.
    /// No phase equilibration or kinetic update is performed. The original
    /// charts must explicitly admit both temperatures, including any metastable
    /// state. Reversing the endpoints reverses the sign exactly on the same ISA.
    pub fn specific_entropy_change_j_kg_k(
        &self,
        initial_temperature_k: f64,
        initial_crystalline_mass_fraction: f64,
        final_temperature_k: f64,
        final_crystalline_mass_fraction: f64,
    ) -> Result<f64, CrystallizationEntropyError> {
        for (temperature, fraction) in [
            (initial_temperature_k, initial_crystalline_mass_fraction),
            (final_temperature_k, final_crystalline_mass_fraction),
        ] {
            self.pair.admit_temperature(temperature)?;
            if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
                return Err(invalid("crystalline mass fraction must be in [0, 1]").into());
            }
        }
        let initial = (initial_temperature_k, initial_crystalline_mass_fraction);
        let final_state = (final_temperature_k, final_crystalline_mass_fraction);
        if initial == final_state {
            return Ok(0.0);
        }
        // Canonical endpoint order makes reversal use the same arithmetic.
        let reverse = final_state < initial;
        let ((ta, xa), (tb, xb)) = if reverse {
            (final_state, initial)
        } else {
            (initial, final_state)
        };
        let parent = if xa == 1.0 {
            0.0
        } else {
            (1.0 - xa) * entropy_between(self.pair.parent(), ta, tb)?
        };
        let crystal = if xa == 0.0 {
            0.0
        } else {
            xa * entropy_between(self.pair.crystal(), ta, tb)?
        };
        // Heat at the old fraction, then change fraction at the final T.
        // This state-function identity avoids subtracting arbitrary absolute
        // phase entropies and keeps the phase-change term explicit.
        let conversion = if xa == xb {
            0.0
        } else {
            -(xb - xa) * self.phase_entropy_gap_j_kg_k(tb)?
        };
        let change = sum([parent, crystal, conversion])?;
        Ok(if reverse { -change } else { change })
    }

    /// Validate the original sources and both accepted kinetic states before
    /// evaluating their entropy difference. The caller retains the caloric
    /// reconstruction tolerance used for those states; it is not an entropy
    /// production tolerance or an enlargement of either source domain.
    pub fn state_entropy_change_j_kg_k(
        &self,
        law: &CrystallizationLaw<'_>,
        initial: CrystallizationState,
        final_state: CrystallizationState,
        solve: CrystallizationSolve,
    ) -> Result<f64, CrystallizationEntropyError> {
        let config = law.config();
        let range = self.pair.temperature_range_k();
        if law.parent_caloric_identity() != Some(*self.pair.parent().identity().as_bytes())
            || law.crystalline_caloric_identity()
                != Some(*self.pair.crystal().identity().as_bytes())
            || config.melting_temperature_k != self.coexistence_temperature_k
            || config.temperature_range_k[0] < range[0]
            || config.temperature_range_k[1] > range[1]
        {
            return Err(CrystallizationEntropyError::LawMismatch);
        }
        law.validate_state(initial, solve)?;
        law.validate_state(final_state, solve)?;
        self.specific_entropy_change_j_kg_k(
            initial.temperature_k(),
            initial.crystalline_mass_fraction(),
            final_state.temperature_k(),
            final_state.crystalline_mass_fraction(),
        )
    }

    /// Material entropy change of the two states retained by a kinetic step.
    /// This revalidates the endpoints, not the reported heat or time accuracy.
    /// Closed adiabatic use additionally requires every other energy and
    /// entropy exchange to be absent; externally supplied heat needs its own
    /// reservoir entropy account at the actual boundary temperature.
    pub fn step_entropy_change_j_kg_k(
        &self,
        law: &CrystallizationLaw<'_>,
        step: &CrystallizationStep,
        solve: CrystallizationSolve,
    ) -> Result<f64, CrystallizationEntropyError> {
        self.state_entropy_change_j_kg_k(law, step.initial, step.state, solve)
    }

    fn entropy_from_coexistence(
        &self,
        temperature_k: f64,
    ) -> Result<(f64, f64), CrystallizationEntropyError> {
        self.pair.admit_temperature(temperature_k)?;
        let parent_h = self.pair.parent().enthalpy_at_temperature(temperature_k)?;
        let crystal_h = self.pair.crystal().enthalpy_at_temperature(temperature_k)?;
        Ok((
            self.pair
                .parent()
                .curve()
                .specific_entropy_change_j_kg_k(self.parent_h_at_coexistence_j_kg, parent_h)?,
            self.pair
                .crystal()
                .curve()
                .specific_entropy_change_j_kg_k(self.crystal_h_at_coexistence_j_kg, crystal_h)?,
        ))
    }
}

fn entropy_between(
    branch: &CrystallizationCaloricBranch<'_>,
    initial_temperature_k: f64,
    final_temperature_k: f64,
) -> Result<f64, CrystallizationEntropyError> {
    let initial_h = branch.enthalpy_at_temperature(initial_temperature_k)?;
    let final_h = branch.enthalpy_at_temperature(final_temperature_k)?;
    Ok(branch
        .curve()
        .specific_entropy_change_j_kg_k(initial_h, final_h)?)
}

fn finite(value: f64) -> Result<f64, CrystallizationEntropyError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("crystallization entropy is not finite and representable").into())
    }
}

fn sum<const N: usize>(terms: [f64; N]) -> Result<f64, CrystallizationEntropyError> {
    let mut total = 0.0;
    let mut correction = 0.0;
    for term in terms {
        let increment = finite(term)? - correction;
        let next = finite(total + increment)?;
        correction = (next - total) - increment;
        total = next;
    }
    Ok(total)
}

#[cfg(test)]
mod tests;
