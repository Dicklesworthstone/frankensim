//! Original material calorics for finite-rate crystallization.
//!
//! A kinetic crystal fraction is an independent internal variable. It cannot
//! be recovered from the equilibrium melting chart or substituted for a liquid
//! fraction: an untransformed polymer can itself be an amorphous solid. This
//! adapter borrows two explicitly supplied single-phase charts over their
//! common temperature domain. Their enthalpies must use one physical datum;
//! the difference `h_parent(T) - h_crystal(T)` owns latent energy exactly once.
//! No supercooled branch, phase datum, melting point or density is extrapolated.
//!
//! Source-backed charts retain their original Cp/density receipts. Authored
//! charts remain authored data. Neither this preparation nor a kinetic fit
//! proves nucleation, growth morphology, thermodynamic phase stability or
//! experimental validity. The kinetic owner's explicit melting/domain limit
//! remains mandatory. Density supplies free volume change, not mechanical
//! locking, residual stress, a moving boundary or pressure work.

/// Prepared source-bound temperature-dependent Nakamura rates.
pub mod kinetics;

use super::{EquilibriumEnthalpyPhaseCurve, PhaseStateError};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::PropertyUsageReceipt;
use fs_thermal::crystallization::{
    CrystallizationCaloric, CrystallizationConfig, CrystallizationError, CrystallizationLaw,
    CrystallizationRate,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// Explicit stress-free locking and original-density shrinkage in native mechanics.
pub mod mechanics;

/// Borrowed, strictly sensible branch with its original source identity.
/// Preparation is linear in knots; temperature queries allocate nothing and
/// use binary search. Every generated Cp breakpoint remains in the chart.
#[derive(Clone, Copy, Debug)]
pub struct CrystallizationCaloricBranch<'a> {
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    minimum_heat_capacity_j_kg_k: f64,
}

impl<'a> CrystallizationCaloricBranch<'a> {
    /// Admit a complete pure solid or pure liquid chart. A coexistence plateau
    /// cannot describe the metastable parent at independently specified X.
    pub fn new(curve: &'a EquilibriumEnthalpyPhaseCurve) -> Result<Self, PhaseStateError> {
        let knots = curve.knots();
        let fraction = knots[0].liquid_mass_fraction;
        if (fraction != 0.0 && fraction != 1.0)
            || knots
                .iter()
                .any(|knot| knot.liquid_mass_fraction != fraction)
        {
            return Err(invalid(
                "crystallization requires separate pure-phase caloric branches",
            ));
        }
        let mut minimum_heat_capacity_j_kg_k = f64::INFINITY;
        for pair in knots.windows(2) {
            let dt = pair[1].temperature_k - pair[0].temperature_k;
            let dh = pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg;
            let cp = dh / dt;
            if !dt.is_finite()
                || dt <= 0.0
                || !cp.is_finite()
                || cp <= 0.0
                || !cp.recip().is_finite()
                || cp.recip() <= 0.0
            {
                return Err(invalid(
                    "crystallization branch must have finite positive sensible Cp",
                ));
            }
            minimum_heat_capacity_j_kg_k = minimum_heat_capacity_j_kg_k.min(cp);
        }
        Ok(Self {
            curve,
            minimum_heat_capacity_j_kg_k,
        })
    }

    #[must_use]
    pub const fn curve(&self) -> &'a EquilibriumEnthalpyPhaseCurve {
        self.curve
    }

    #[must_use]
    pub fn source_receipts(&self) -> &'a [PropertyUsageReceipt] {
        self.curve.source_receipts()
    }

    #[must_use]
    pub fn identity(&self) -> ContentHash {
        self.curve.identity()
    }

    #[must_use]
    pub fn temperature_range_k(&self) -> [f64; 2] {
        let knots = self.curve.knots();
        [knots[0].temperature_k, knots[knots.len() - 1].temperature_k]
    }

    /// Lower bound on dh/dT for the actual piecewise-linear h(T) chart,
    /// including its stated sensible interpolation error; not a source bound.
    #[must_use]
    pub const fn minimum_heat_capacity_j_kg_k(&self) -> f64 {
        self.minimum_heat_capacity_j_kg_k
    }

    /// Invert the original sensible chart without replacing its interpolation.
    pub fn enthalpy_at_temperature(&self, temperature_k: f64) -> Result<f64, PhaseStateError> {
        let [low, high] = self.temperature_range_k();
        if !temperature_k.is_finite() || temperature_k < low || temperature_k > high {
            return Err(invalid(
                "temperature is outside the original crystallization branch",
            ));
        }
        let knots = self.curve.knots();
        let upper = knots.partition_point(|knot| knot.temperature_k < temperature_k);
        if knots[upper].temperature_k == temperature_k {
            return Ok(knots[upper].specific_enthalpy_j_kg);
        }
        let a = knots[upper - 1];
        let b = knots[upper];
        let weight = (temperature_k - a.temperature_k) / (b.temperature_k - a.temperature_k);
        let dh = b.specific_enthalpy_j_kg - a.specific_enthalpy_j_kg;
        let h = if weight <= 0.5 {
            a.specific_enthalpy_j_kg + weight * dh
        } else {
            b.specific_enthalpy_j_kg - (1.0 - weight) * dh
        };
        if !h.is_finite() {
            return Err(invalid(
                "crystallization branch enthalpy is not representable",
            ));
        }
        Ok(h)
    }

    /// Original source density at the same caloric coordinate. No volume or
    /// density interpolation is invented by the crystallization adapter.
    pub fn density_at_temperature(&self, temperature_k: f64) -> Result<f64, PhaseStateError> {
        let h = self.enthalpy_at_temperature(temperature_k)?;
        self.curve
            .values_at_specific_enthalpy(h)
            .map(|state| state.bulk_density_kg_m3)
    }
}

impl CrystallizationCaloric for CrystallizationCaloricBranch<'_> {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some(*self.identity().as_bytes())
    }
    fn minimum_heat_capacity_j_kg_k(&self) -> Option<f64> {
        Some(self.minimum_heat_capacity_j_kg_k)
    }
    fn specific_enthalpy_j_kg(&self, temperature_k: f64) -> Result<f64, CrystallizationError> {
        self.enthalpy_at_temperature(temperature_k).map_err(|_| {
            CrystallizationError::Invalid(
                "temperature is outside the original material caloric branch",
            )
        })
    }
}

impl ThermalLaw for CrystallizationCaloricBranch<'_> {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some(*self.identity().as_bytes())
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        ThermalLaw::temperature_k(self.curve, h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        self.minimum_heat_capacity_j_kg_k.recip()
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h: f64,
        final_h: f64,
    ) -> Result<Option<f64>, DepositionError> {
        ThermalLaw::specific_entropy_change_j_kg_k(self.curve, initial_h, final_h)
    }
}

/// Two immutable material branches with a positive crystallization enthalpy
/// gap over a declared common domain. The parent can be solid or liquid; the
/// product must be the supplied pure crystalline solid. Pressure, composition,
/// phase identity and a compatible absolute energy reference are caller inputs,
/// retained through the original chart identities and receipts.
#[derive(Clone, Copy, Debug)]
pub struct CrystallizationCaloricPair<'a> {
    parent: CrystallizationCaloricBranch<'a>,
    crystal: CrystallizationCaloricBranch<'a>,
    temperature_range_k: [f64; 2],
    identity: ContentHash,
}

impl<'a> CrystallizationCaloricPair<'a> {
    pub fn new(
        parent: &'a EquilibriumEnthalpyPhaseCurve,
        crystal: &'a EquilibriumEnthalpyPhaseCurve,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, PhaseStateError> {
        if parent.material_card_identity() != crystal.material_card_identity() {
            return Err(invalid(
                "parent and crystal branches must retain the same material card",
            ));
        }
        let parent = CrystallizationCaloricBranch::new(parent)?;
        let crystal = CrystallizationCaloricBranch::new(crystal)?;
        if crystal.curve.knots()[0].liquid_mass_fraction != 0.0 {
            return Err(invalid(
                "the crystalline product must be an explicitly solid branch",
            ));
        }
        let [low, high] = temperature_range_k;
        if !low.is_finite() || !high.is_finite() || low <= 0.0 || high <= low {
            return Err(invalid(
                "crystallization needs a finite positive temperature interval",
            ));
        }
        for branch in [parent, crystal] {
            let span = branch.temperature_range_k();
            if low < span[0] || high > span[1] {
                return Err(invalid(
                    "both original branches must cover the entire crystallization interval",
                ));
            }
        }
        // h_parent-h_crystal is affine between the UNION of both knot sets.
        // Positive endpoints therefore establish a positive gap throughout the
        // represented domain, including unequal and temperature-dependent Cp.
        for temperature_k in [low, high].into_iter().chain(
            parent
                .curve
                .knots()
                .iter()
                .chain(crystal.curve.knots())
                .map(|knot| knot.temperature_k)
                .filter(|t| *t > low && *t < high),
        ) {
            let gap = parent.enthalpy_at_temperature(temperature_k)?
                - crystal.enthalpy_at_temperature(temperature_k)?;
            if !gap.is_finite() || gap <= 0.0 {
                return Err(invalid(
                    "crystallization requires a positive original-phase enthalpy gap throughout its domain",
                ));
            }
        }
        let mut hash =
            DomainHasher::new("org.frankensim.fs-material.crystallization-caloric-pair.v1");
        hash.update(parent.identity().as_bytes());
        hash.update(crystal.identity().as_bytes());
        for value in temperature_range_k {
            hash.update(&value.to_bits().to_le_bytes());
        }
        Ok(Self {
            parent,
            crystal,
            temperature_range_k,
            identity: hash.finalize(),
        })
    }

    #[must_use]
    pub const fn parent(&self) -> &CrystallizationCaloricBranch<'a> {
        &self.parent
    }
    #[must_use]
    pub const fn crystal(&self) -> &CrystallizationCaloricBranch<'a> {
        &self.crystal
    }
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    #[must_use]
    pub const fn temperature_range_k(&self) -> [f64; 2] {
        self.temperature_range_k
    }

    /// Bind this exact caloric pair and explicitly identified rate law to the
    /// shared kinetic owner. `kinetic_source_identity` must identify the supplied
    /// K(T) data, process and applicability; a digest is not calibration evidence.
    /// All parameters and both original chart identities enter the state identity.
    /// The parent/crystal source domains must explicitly support supercooling.
    pub fn law<'b>(
        &'b self,
        rate: &'b dyn CrystallizationRate,
        kinetic_source_identity: ContentHash,
        avrami_exponent: f64,
        maximum_crystalline_mass_fraction: f64,
        melting_temperature_k: f64,
    ) -> Result<CrystallizationLaw<'b>, CrystallizationError> {
        let mut hash = DomainHasher::new("org.frankensim.fs-material.crystallization-law.v1");
        hash.update(self.identity.as_bytes());
        hash.update(kinetic_source_identity.as_bytes());
        for value in [
            avrami_exponent,
            maximum_crystalline_mass_fraction,
            melting_temperature_k,
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        CrystallizationLaw::new(
            &self.parent,
            &self.crystal,
            rate,
            CrystallizationConfig {
                identity: *hash.finalize().as_bytes(),
                avrami_exponent,
                maximum_crystalline_mass_fraction,
                temperature_range_k: self.temperature_range_k,
                melting_temperature_k,
            },
        )
    }

    /// Global thermal slope bound at any FIXED crystalline mass fraction.
    /// A kinetics integrator must additionally resolve changing crystallinity;
    /// this conduction bound alone does not bound kinetic integration error.
    #[must_use]
    pub fn maximum_temperature_slope(&self) -> f64 {
        self.parent
            .minimum_heat_capacity_j_kg_k
            .min(self.crystal.minimum_heat_capacity_j_kg_k)
            .recip()
    }

    pub fn latent_enthalpy_j_kg(&self, temperature_k: f64) -> Result<f64, PhaseStateError> {
        self.admit_temperature(temperature_k)?;
        Ok(self.parent.enthalpy_at_temperature(temperature_k)?
            - self.crystal.enthalpy_at_temperature(temperature_k)?)
    }

    /// Shared additive-specific-volume state at the ACTUAL kinetic fraction.
    /// This does not re-equilibrate X from T. Mechanical and moving-volume
    /// consumers must account separately for constraints and pressure work.
    pub fn volume_state(
        &self,
        temperature_k: f64,
        crystalline_mass_fraction: f64,
    ) -> Result<CrystallizationVolumeState, PhaseStateError> {
        self.admit_temperature(temperature_k)?;
        let x = crystalline_mass_fraction;
        if !x.is_finite() || !(0.0..=1.0).contains(&x) {
            return Err(invalid("crystalline mass fraction must be in [0, 1]"));
        }
        let parent_density_kg_m3 = self.parent.density_at_temperature(temperature_k)?;
        let crystal_density_kg_m3 = self.crystal.density_at_temperature(temperature_k)?;
        let vc = x / crystal_density_kg_m3;
        let vp = (1.0 - x) / parent_density_kg_m3;
        let specific_volume_m3_kg = vc + vp;
        let bulk_density_kg_m3 = specific_volume_m3_kg.recip();
        if !specific_volume_m3_kg.is_finite()
            || specific_volume_m3_kg <= 0.0
            || !bulk_density_kg_m3.is_finite()
            || bulk_density_kg_m3 <= 0.0
            || (x > 0.0 && vc == 0.0)
            || (x < 1.0 && vp == 0.0)
        {
            return Err(invalid(
                "crystallization phase volumes are not representable",
            ));
        }
        let crystalline_volume_fraction = if x == 0.0 || x == 1.0 {
            x
        } else if vc <= vp {
            vc / specific_volume_m3_kg
        } else {
            1.0 - vp / specific_volume_m3_kg
        };
        Ok(CrystallizationVolumeState {
            caloric_pair_identity: self.identity,
            temperature_k,
            crystalline_mass_fraction: x,
            crystalline_volume_fraction,
            parent_density_kg_m3,
            crystal_density_kg_m3,
            specific_volume_m3_kg,
            bulk_density_kg_m3,
        })
    }

    fn admit_temperature(&self, temperature_k: f64) -> Result<(), PhaseStateError> {
        if !temperature_k.is_finite()
            || temperature_k < self.temperature_range_k[0]
            || temperature_k > self.temperature_range_k[1]
        {
            Err(invalid(
                "temperature is outside the admitted crystallization pair",
            ))
        } else {
            Ok(())
        }
    }
}

/// One fixed-mass mixture volume, with mass and volume fractions distinguished.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrystallizationVolumeState {
    pub caloric_pair_identity: ContentHash,
    pub temperature_k: f64,
    pub crystalline_mass_fraction: f64,
    pub crystalline_volume_fraction: f64,
    pub parent_density_kg_m3: f64,
    pub crystal_density_kg_m3: f64,
    pub specific_volume_m3_kg: f64,
    pub bulk_density_kg_m3: f64,
}

fn invalid(what: &'static str) -> PhaseStateError {
    PhaseStateError::InvalidCurve { what }
}

#[cfg(test)]
mod tests;
