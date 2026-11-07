//! Verified same-feed coexistence between two original local phase ensembles.
//!
//! A genuine latent interval requires equal temperature and every present
//! species' chemical potential, not merely overlapping temperatures or a
//! crossing of total Gibbs energies. Both endpoints retain their original
//! physical caloric sources. Their actual positive DeltaH and DeltaS must
//! close DeltaH/T = DeltaS within the explicitly funded numerical tolerance.
//! The measured residual is retained; no entropy value is fitted to the heat.
//!
//! This admits one offered common tangent. It does not discover its temperature,
//! establish global Gibbs stability, or introduce nucleation/interphase kinetics.

use super::super::{MAX_REGULAR_SOLUTION_SPECIES_V1, Sum, positive};
use super::transactions::MulticomponentCaloricInventoryRefV1;
use super::{
    MulticomponentEquilibriumCaloricErrorV1, MulticomponentEquilibriumCaloricStateV1,
    MulticomponentLiquidCaloricModelV1, Result, Temperature, checked, poll,
};
use crate::solution::multicomponent::equilibrium::MulticomponentEquilibriumPhaseV1;

/// Accuracy and retained-phase work limits for an offered common tangent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentCaloricCoexistenceControlV1 {
    /// Absolute allowance for every original dimensionless mu_i/(R T).
    pub chemical_potential_tolerance: f64,
    /// Absolute allowance for original DeltaS - DeltaH/T [J/K]. A 64-epsilon
    /// arithmetic allowance is added using the two actual latent quantities.
    pub entropy_tolerance_j_per_k: f64,
    /// Bound on the sum of both original ensemble phase counts. No duplicate
    /// phase identification or unreported morphology reduction is performed.
    pub maximum_retained_phases: usize,
}
impl Default for MulticomponentCaloricCoexistenceControlV1 {
    fn default() -> Self {
        Self {
            chemical_potential_tolerance: 1e-10,
            entropy_tolerance_j_per_k: 1e-8,
            maximum_retained_phases: 2 * MAX_REGULAR_SOLUTION_SPECIES_V1,
        }
    }
}

/// Two immutable same-feed ensembles and their original latent storage.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentCaloricCoexistenceV1 {
    endpoints: [MulticomponentEquilibriumCaloricStateV1; 2],
    latent_enthalpy_j: f64,
    latent_entropy_j_per_k: f64,
    latent_clausius_residual_j_per_k: f64,
    maximum_chemical_potential_residual: f64,
    controls: MulticomponentCaloricCoexistenceControlV1,
}

impl MulticomponentLiquidCaloricModelV1 {
    /// Admit an offered lower-/higher-enthalpy coexistence pair at exactly the
    /// same temperature and total species inventory. Every original endpoint
    /// must belong to this complete caloric owner and carry a resolved local
    /// phase response. In particular, finite discovery alone cannot mint this
    /// common-tangent relation or fill a first-order enthalpy gap.
    pub fn coexistence_plateau(
        &self,
        lower: &MulticomponentEquilibriumCaloricStateV1,
        upper: &MulticomponentEquilibriumCaloricStateV1,
        controls: MulticomponentCaloricCoexistenceControlV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentCaloricCoexistenceV1> {
        poll(&mut cancelled)?;
        self.admit_equilibrium_caloric(lower)?;
        self.admit_equilibrium_caloric(upper)?;
        if !controls.chemical_potential_tolerance.is_finite()
            || !(4.0 * f64::EPSILON..=1e-3).contains(&controls.chemical_potential_tolerance)
            || !controls.entropy_tolerance_j_per_k.is_finite()
            || controls.entropy_tolerance_j_per_k < 0.0
            || !(1..=2 * MAX_REGULAR_SOLUTION_SPECIES_V1)
                .contains(&controls.maximum_retained_phases)
            || lower.temperature() != upper.temperature()
            || lower.component_moles() != upper.component_moles()
        {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(
                "coexistence needs bounded controls, one exact temperature and the same feed",
            ));
        }
        if lower.candidate().phases().len() + upper.candidate().phases().len()
            > controls.maximum_retained_phases
        {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(
                "coexistence retained-phase budget is exhausted",
            ));
        }
        let reference = lower.candidate().phases()[0]
            .thermodynamics()
            .reduced_chemical_potentials();
        let mut maximum_chemical_potential_residual: f64 = 0.0;
        for phase in lower
            .candidate()
            .phases()
            .iter()
            .chain(upper.candidate().phases())
        {
            poll(&mut cancelled)?;
            let potential = phase.thermodynamics().reduced_chemical_potentials();
            for (i, &amount) in lower.component_moles().iter().enumerate() {
                if amount == 0.0 {
                    continue;
                }
                let a = reference[i].expect("accepted present component");
                let b = potential[i].expect("accepted present component");
                maximum_chemical_potential_residual = maximum_chemical_potential_residual.max(
                    checked((b - a).abs(), "coexistence chemical-potential difference")?,
                );
            }
        }
        if maximum_chemical_potential_residual > controls.chemical_potential_tolerance {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(
                "coexisting ensembles do not share every original chemical potential",
            ));
        }
        // The same conserved homogeneous feed cancels before arithmetic, so
        // neither pure reference enthalpies nor a phase-count pairing enter.
        let feed = lower.candidate().feed().mole_fractions();
        let low = self.mixing_changes(
            MulticomponentCaloricInventoryRefV1::from(lower).phases(),
            feed,
            &mut cancelled,
        )?;
        let high = self.mixing_changes(
            MulticomponentCaloricInventoryRefV1::from(upper).phases(),
            feed,
            &mut cancelled,
        )?;
        let difference = |index: usize| -> Result<f64> {
            let mut sum = Sum::default();
            sum.add(low[index]);
            sum.add(-high[index]);
            Ok(checked(sum.value(), "original coexistence latent storage")?)
        };
        let latent_enthalpy_j = difference(0)?;
        let latent_entropy_j_per_k = difference(1)?;
        if !positive(latent_enthalpy_j) || !positive(latent_entropy_j_per_k) {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(
                "coexistence requires ordered positive original latent heat and entropy",
            ));
        }
        let heat_entropy = checked(
            latent_enthalpy_j / lower.temperature().value(),
            "latent Clausius entropy",
        )?;
        let latent_clausius_residual_j_per_k = checked(
            latent_entropy_j_per_k - heat_entropy,
            "original latent Clausius residual",
        )?;
        let allowance = checked(
            controls.entropy_tolerance_j_per_k
                + 64.0 * f64::EPSILON * latent_entropy_j_per_k.max(heat_entropy),
            "latent Clausius allowance",
        )?;
        if latent_clausius_residual_j_per_k.abs() > allowance {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(
                "original latent enthalpy and entropy do not close their common-temperature balance",
            ));
        }
        let coexistence = MulticomponentCaloricCoexistenceV1 {
            endpoints: [lower.clone(), upper.clone()],
            latent_enthalpy_j,
            latent_entropy_j_per_k,
            latent_clausius_residual_j_per_k,
            maximum_chemical_potential_residual,
            controls,
        };
        poll(&mut cancelled)?;
        Ok(coexistence)
    }
}

impl MulticomponentCaloricCoexistenceV1 {
    pub const fn lower(&self) -> &MulticomponentEquilibriumCaloricStateV1 {
        &self.endpoints[0]
    }
    pub const fn upper(&self) -> &MulticomponentEquilibriumCaloricStateV1 {
        &self.endpoints[1]
    }
    pub fn temperature(&self) -> Temperature {
        self.lower().temperature()
    }
    pub fn component_moles(&self) -> &[f64] {
        self.lower().component_moles()
    }
    pub const fn latent_enthalpy_j(&self) -> f64 {
        self.latent_enthalpy_j
    }
    pub const fn latent_entropy_j_per_k(&self) -> f64 {
        self.latent_entropy_j_per_k
    }
    /// Actual source DeltaS minus DeltaH/T; it is never fitted away.
    pub const fn latent_clausius_residual_j_per_k(&self) -> f64 {
        self.latent_clausius_residual_j_per_k
    }
    pub const fn maximum_chemical_potential_residual(&self) -> f64 {
        self.maximum_chemical_potential_residual
    }
    pub const fn controls(&self) -> MulticomponentCaloricCoexistenceControlV1 {
        self.controls
    }

    /// Recover the actual convex phase inventory at a heat coordinate in
    /// [0, DeltaH]. The material mass fraction equals this ensemble fraction
    /// because both endpoints have exactly the same conserved species feed.
    pub fn state_at_enthalpy_increment(
        &self,
        heat_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentCoexistenceStateV1<'_>> {
        poll(&mut cancelled)?;
        if !heat_j.is_finite() || heat_j < 0.0 || heat_j > self.latent_enthalpy_j {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(
                "heat outside the original coexistence interval",
            ));
        }
        let fraction = heat_j / self.latent_enthalpy_j;
        if heat_j > 0.0 && !positive(fraction) {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Unrepresentable(
                "positive coexistence fraction",
            ));
        }
        let state = MulticomponentCoexistenceStateV1 {
            coexistence: self,
            upper_fraction: fraction,
            heat_j,
        };
        for phase in state.phases() {
            poll(&mut cancelled)?;
            if !positive(phase.total_moles())
                || phase
                    .phase
                    .component_moles()
                    .iter()
                    .any(|&amount| amount > 0.0 && !positive(amount * phase.ensemble_fraction))
            {
                return Err(MulticomponentEquilibriumCaloricErrorV1::Unrepresentable(
                    "positive coexistence phase or component inventory",
                ));
            }
        }
        poll(&mut cancelled)?;
        Ok(state)
    }
}

/// Borrowed physical phase inventory on an admitted isothermal latent interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentCoexistenceStateV1<'a> {
    coexistence: &'a MulticomponentCaloricCoexistenceV1,
    upper_fraction: f64,
    heat_j: f64,
}
impl<'a> MulticomponentCoexistenceStateV1<'a> {
    pub const fn coexistence(&self) -> &'a MulticomponentCaloricCoexistenceV1 {
        self.coexistence
    }
    pub fn temperature(&self) -> Temperature {
        self.coexistence.temperature()
    }
    pub const fn upper_ensemble_fraction(&self) -> f64 {
        self.upper_fraction
    }
    pub const fn enthalpy_increment_j(&self) -> f64 {
        self.heat_j
    }
    /// Actual original source entropy carried by the interpolated inventory.
    pub fn entropy_increment_j_per_k(&self) -> f64 {
        self.upper_fraction * self.coexistence.latent_entropy_j_per_k
    }
    /// Rounding of the represented convex inventory relative to supplied heat.
    pub fn enthalpy_residual_j(&self) -> f64 {
        self.upper_fraction * self.coexistence.latent_enthalpy_j - self.heat_j
    }
    pub fn phases(&self) -> impl Iterator<Item = MulticomponentCoexistingPhaseV1<'a>> {
        let low = 1.0 - self.upper_fraction;
        let high = self.upper_fraction;
        self.coexistence
            .lower()
            .candidate()
            .phases()
            .iter()
            .map(move |phase| MulticomponentCoexistingPhaseV1 {
                phase,
                ensemble_fraction: low,
            })
            .chain(
                self.coexistence
                    .upper()
                    .candidate()
                    .phases()
                    .iter()
                    .map(move |phase| MulticomponentCoexistingPhaseV1 {
                        phase,
                        ensemble_fraction: high,
                    }),
            )
            .filter(|phase| phase.ensemble_fraction > 0.0)
    }
}

/// One original phase, with its amount scaled by the actual ensemble fraction.
/// Its original composition, chemical potentials and source identity survive.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentCoexistingPhaseV1<'a> {
    phase: &'a MulticomponentEquilibriumPhaseV1,
    ensemble_fraction: f64,
}
impl<'a> MulticomponentCoexistingPhaseV1<'a> {
    pub const fn original_phase(&self) -> &'a MulticomponentEquilibriumPhaseV1 {
        self.phase
    }
    pub const fn ensemble_fraction(&self) -> f64 {
        self.ensemble_fraction
    }
    pub fn total_moles(&self) -> f64 {
        self.ensemble_fraction * self.phase.total_moles()
    }
    pub fn component_moles(&self) -> impl Iterator<Item = f64> + '_ {
        self.phase
            .component_moles()
            .iter()
            .map(|&n| self.ensemble_fraction * n)
    }
}
