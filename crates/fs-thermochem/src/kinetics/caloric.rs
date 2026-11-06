//! Shared source-path admission for kinetic reactors and compressible flow.

use super::{BoundReaction, IdealGasKineticsV1, KineticsErrorV1, R, poll, positive};
pub(super) use crate::caloric::ReactionCaloricBasis;
use crate::caloric::ReactionCaloricSum;
use crate::equilibrium::composition::IdealGasEquilibriumV1;
use crate::{MolarEnergyQuantityV1, MolarThermalQuantityV1, Nasa9StandardStateModelV1};
use fs_qty::Temperature;

/// Original standard-state caloric changes per mole of reaction progress.
/// Chemical formation differences remain; only conserved references cancel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactionCaloricChangeV1 {
    /// Delta h = sum(nu_i h_i), including ideal-gas pressure-flow work.
    pub molar_enthalpy: MolarEnergyQuantityV1,
    /// Delta u = sum(nu_i u_i); Delta h - Delta u = Delta nu R T.
    pub molar_internal_energy: MolarEnergyQuantityV1,
    /// Delta g under the original standard-state pressure and entropy reference.
    pub molar_gibbs_energy: MolarEnergyQuantityV1,
    /// Delta s from the original entropy coefficients, not (Delta h - Delta g)/T.
    pub molar_entropy: MolarThermalQuantityV1,
}

impl IdealGasKineticsV1 {
    /// Source-consistent caloric changes for one canonical reaction index.
    ///
    /// The same coefficient reduction used by full kinetics evaluates Delta h
    /// and Delta u independently, retaining the pressure-work difference even
    /// when either change is small. No absolute formation-loaded species
    /// energies are subtracted. Distinct source formation energies, raw region
    /// jumps and explicit caloric corrections remain physically active.
    ///
    /// Every original source, including spectators, must admit the temperature
    /// and positive Cv. Successful calls allocate no workspace. This state
    /// function is not an additional heat source for a total-chemical-energy
    /// formulation; consumers must retain their original energy convention.
    pub fn reaction_caloric_change_at(
        &self,
        reaction: usize,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactionCaloricChangeV1, KineticsErrorV1> {
        poll(&mut cancelled)?;
        let reaction = self
            .bound
            .get(reaction)
            .ok_or(KineticsErrorV1::Invalid("reaction caloric index"))?;
        let basis = self.admitted_caloric_basis(temperature, &mut cancelled)?;
        let change = reaction_calorics(&self.thermodynamics, reaction, &basis)?;
        poll(&mut cancelled)?;
        Ok(change)
    }

    /// Original molar internal-energy change of each reaction, J/mol of progress.
    ///
    /// Combine the original NASA-9 coefficients with the exact stoichiometric
    /// column before evaluating the polynomial. Conserved formation references
    /// cancel before rounding; distinct formation energies, original region
    /// jumps and explicit caloric corrections remain. This is a reaction state
    /// function, not a heat source to add to a total-chemical-energy equation.
    /// Every original source, including spectators, must support the temperature
    /// with finite properties and positive Cv, as in `evaluate_thermodynamics`.
    /// Output follows `reactions()`. No successful-call allocation; output is
    /// scratch after any refusal or cancellation.
    pub fn reaction_internal_energies_at(
        &self,
        temperature: Temperature,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), KineticsErrorV1> {
        poll(&mut cancelled)?;
        if !positive(temperature.value()) || out.len() != self.bound.len() {
            return Err(KineticsErrorV1::Invalid(
                "reaction caloric temperature or output shape",
            ));
        }
        let basis = self.admitted_caloric_basis(temperature, &mut cancelled)?;
        for (reaction, value) in self.bound.iter().zip(out) {
            poll(&mut cancelled)?;
            *value = reaction_calorics(&self.thermodynamics, reaction, &basis)?
                .molar_internal_energy
                .value();
        }
        poll(&mut cancelled)
    }

    /// Original reaction enthalpy changes in canonical reaction order.
    /// Like `reaction_internal_energies_at`, this combines source coefficients
    /// before evaluation, admits every original species and allocates nothing.
    /// It is the fixed-pressure caloric derivative, not additional reaction heat.
    pub fn reaction_enthalpies_at(
        &self,
        temperature: Temperature,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), KineticsErrorV1> {
        poll(&mut cancelled)?;
        if !positive(temperature.value()) || out.len() != self.bound.len() {
            return Err(KineticsErrorV1::Invalid(
                "reaction caloric temperature or output shape",
            ));
        }
        let basis = self.admitted_caloric_basis(temperature, &mut cancelled)?;
        for (reaction, value) in self.bound.iter().zip(out) {
            poll(&mut cancelled)?;
            *value = reaction_calorics(&self.thermodynamics, reaction, &basis)?
                .molar_enthalpy
                .value();
        }
        poll(&mut cancelled)
    }

    /// Original standard-state reaction entropy changes, J/(mol K) of progress.
    ///
    /// Source coefficients are reduced before evaluation, cancelling conserved
    /// elemental entropy references while retaining every reaction entropy
    /// difference and explicit source correction. This does not subtract Gibbs
    /// energy from enthalpy, which can lose entropy at large formation energies.
    /// Ideal-gas concentration/mixing entropy remains the consumer's responsibility.
    /// As with the energy queries, every original source is admitted and successful
    /// calls allocate nothing; output follows `reactions()` and is scratch on error.
    pub fn reaction_entropies_at(
        &self,
        temperature: Temperature,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), KineticsErrorV1> {
        poll(&mut cancelled)?;
        if !positive(temperature.value()) || out.len() != self.bound.len() {
            return Err(KineticsErrorV1::Invalid(
                "reaction caloric temperature or output shape",
            ));
        }
        let basis = self.admitted_caloric_basis(temperature, &mut cancelled)?;
        for (reaction, value) in self.bound.iter().zip(out) {
            poll(&mut cancelled)?;
            *value = reaction_calorics(&self.thermodynamics, reaction, &basis)?
                .molar_entropy
                .value();
        }
        poll(&mut cancelled)
    }

    // Full rate evaluation has already admitted these properties. Standalone
    // reaction queries enforce that same source semantics once per call.
    fn admitted_caloric_basis(
        &self,
        temperature: Temperature,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ReactionCaloricBasis, KineticsErrorV1> {
        if !positive(temperature.value()) {
            return Err(KineticsErrorV1::Invalid("reaction caloric temperature"));
        }
        for model in self.thermodynamics.models() {
            poll(cancelled)?;
            let properties = model.evaluate_properties(temperature).map_err(|source| {
                KineticsErrorV1::SpeciesEvaluation {
                    species: model.species().clone(),
                    source,
                }
            })?;
            if !positive(properties.cp().value() - R) {
                return Err(KineticsErrorV1::Invalid(
                    "kinetic source requires positive species Cv",
                ));
            }
        }
        Ok(ReactionCaloricBasis::new(temperature))
    }
}

// Equilibrium constants and kinetic reverse rates share the original source
// coefficient reduction, including explicit h/s corrections.
pub(super) fn reaction_calorics(
    source: &IdealGasEquilibriumV1,
    reaction: &BoundReaction,
    basis: &ReactionCaloricBasis,
) -> Result<ReactionCaloricChangeV1, KineticsErrorV1> {
    let mut calorics = ReactionCaloricSum::default();
    for &(i, coefficient) in &reaction.net {
        let model = &source.models()[i];
        calorics
            .add_model(model, coefficient, basis.temperature)
            .map_err(|source| KineticsErrorV1::SpeciesEvaluation {
                species: model.species().clone(),
                source,
            })?;
    }
    let (h, u, g, s) = calorics.evaluate(basis, reaction.delta_molecules);
    if !h.is_finite() || !u.is_finite() || !g.is_finite() || !s.is_finite() {
        return Err(KineticsErrorV1::Arithmetic("reaction caloric reduction"));
    }
    Ok(ReactionCaloricChangeV1 {
        molar_enthalpy: MolarEnergyQuantityV1::new(h),
        molar_internal_energy: MolarEnergyQuantityV1::new(u),
        molar_gibbs_energy: MolarEnergyQuantityV1::new(g),
        molar_entropy: MolarThermalQuantityV1::new(s),
    })
}

/// Admit source support, explicit caloric continuity and strictly positive Cv
/// over the entire original NASA-9 temperature interval.
///
/// Outward-rounded original monomial bounds are subdivided at most12levels and
///4095pieces per source region. An unresolved bound refuses; sampling positive
/// values is never treated as proof that an intervening polynomial stays positive.
/// No source coefficients or formation/entropy references are altered.
pub fn admit_caloric_interval(
    model: &Nasa9StandardStateModelV1,
    lower: Temperature,
    upper: Temperature,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(), KineticsErrorV1> {
    poll(&mut cancelled)?;
    let lo = lower.value();
    let hi = upper.value();
    if !positive(lo) || !positive(hi) || lo > hi {
        return Err(KineticsErrorV1::Invalid("invalid caloric source interval"));
    }
    let first = model
        .region_for_temperature(lo)
        .ok_or(KineticsErrorV1::Invalid(
            "lower temperature lacks source support",
        ))?
        .0;
    let last = model
        .region_for_temperature(hi)
        .ok_or(KineticsErrorV1::Invalid(
            "upper temperature lacks source support",
        ))?
        .0;
    if first != last && model.caloric_reference().is_none() {
        return Err(KineticsErrorV1::Invalid(
            "cross-region path requires explicit caloric continuity",
        ));
    }
    for region in &model.regions()[first..=last] {
        positive_cv(
            *region.coefficients(),
            lo.max(region.temperature_min().value()),
            hi.min(region.temperature_max().value()),
            &mut cancelled,
        )?;
    }
    Ok(())
}

// Each original Cv/R monomial is monotone on positive T. Directed rounding
// bounds every multiply and sum, including the cp/R constant minus one.
fn cv_lower_bound(a: [f64; 9], lo: f64, hi: f64) -> f64 {
    let mut low = (a[2] - 1.0).next_down();
    for (i, (coefficient, l, h)) in [
        (a[0], (1.0 / hi).next_down(), (1.0 / lo).next_up()),
        (a[1], (1.0 / hi).next_down(), (1.0 / lo).next_up()),
        (a[3], lo, hi),
        (a[4], lo, hi),
        (a[5], lo, hi),
        (a[6], lo, hi),
    ]
    .into_iter()
    .enumerate()
    {
        if coefficient == 0.0 {
            continue;
        }
        let power = match i {
            0 => 2,
            1 | 2 => 1,
            3 => 2,
            4 => 3,
            _ => 4,
        };
        let (mut lower, mut upper) = (1.0_f64, 1.0_f64);
        for _ in 0..power {
            lower = (lower * l).next_down();
            upper = (upper * h).next_up();
        }
        let term = coefficient * if coefficient > 0.0 { lower } else { upper };
        low = (low + term.next_down()).next_down();
    }
    low
}
fn positive_cv(
    a: [f64; 9],
    lo: f64,
    hi: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), KineticsErrorV1> {
    let mut stack = [(0.0, 0.0, 0_u8); 14];
    stack[0] = (lo, hi, 0);
    let mut n = 1;
    let mut work = 0;
    while n > 0 {
        poll(cancelled)?;
        if work >= 4095 {
            return Err(KineticsErrorV1::Invalid(
                "positive Cv interval subdivision budget exhausted",
            ));
        }
        work += 1;
        n -= 1;
        let (l, h, depth) = stack[n];
        let bound = cv_lower_bound(a, l, h);
        if bound.is_finite() && bound > 0.0 {
            continue;
        }
        let mid = l + 0.5 * (h - l);
        if depth == 12 || mid <= l || mid >= h {
            return Err(KineticsErrorV1::Invalid(
                "positive component Cv could not be bounded on the source path",
            ));
        }
        stack[n] = (mid, h, depth + 1);
        stack[n + 1] = (l, mid, depth + 1);
        n += 2;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g0_interval_cv_refuses_an_interior_negative_pocket_and_preserves_cancellation() {
        // Cv/R=(T-500)^2/10000-1: both endpoints are positive, midpoint negative.
        let mut a = [0.0; 9];
        a[2] = 25.0;
        a[3] = -0.1;
        a[4] = 0.0001;
        assert!(cv_lower_bound(a, 200.0, 800.0) < 0.0);
        assert!(positive_cv(a, 200.0, 800.0, &mut || false).is_err());
        a[2] = 27.0;
        assert!(positive_cv(a, 200.0, 800.0, &mut || false).is_ok());
        assert_eq!(
            positive_cv(a, 200.0, 800.0, &mut || true),
            Err(KineticsErrorV1::Cancelled)
        );
    }
}
