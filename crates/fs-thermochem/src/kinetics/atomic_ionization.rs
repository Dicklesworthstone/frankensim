//! One ground-state atomic ionization channel at separate Te and Th.
//!
//! The retained reaction is `e + A <-> I + 2e`, with a neutral monatomic A,
//! its singly charged ground-state ion I, Maxwellian electrons, and stationary
//! massive heavy targets. All three active sources must have Cp/R=5/2, fixed
//! formation constants and fixed ground-state statistical weights. Molecular
//! internal modes, excited-state partition functions, radiation and Coulomb
//! corrections are not inferred. This is a reversible collisional CHANNEL,
//! not a complete collisional-radiative or total recombination model.
//!
//! Maxwellian microscopic reciprocity gives `k_3(Te)=k_ion(Te)/S(Te)`;
//! see NIST FLYCHK manual, pp.29–30, including the differential Fowler relation:
//! <https://nlte.nist.gov/FLY/Doc/Manual_FLYCHK_Nov08.pdf>.
//! A Saha relation at Te describes this specified kinetic partial balance,
//! not arbitrary two-temperature plasma equilibrium: van der Mullen et al.,
//! Phys. Rev. E 50,3925 (1994), <https://doi.org/10.1103/PhysRevE.50.3925>.
//!
//! The molar concentration factor is derived from ORIGINAL monatomic NASA9
//! constants: `ln Kc = delta_a8 - 5/2 + ln(p0/R) + 3/2 ln Te - chi/(R Te)`.
//! Its statistical value is independently checked against the declared level
//! weights and actual source masses, including `(m_I/m_A)^(3/2)`. No generic
//! heavy-species Gibbs function is evaluated at the electron temperature.
//!
//! With ground-state chemical energy chi [J/mol], electron sensible energy
//! changes by `-chi*q`, while heavy sensible energy is unchanged in the
//! stationary-target approximation. Changing electron population also changes
//! its heat capacity: `Cv_e*dTe/dt=-(chi+3RT_e/2)*q`. The population term is
//! not an additional heat sink in the conserved energy equation. Under these
//! assumptions the complete two-pool ideal entropy production is exactly
//! `R*q*ln(Kc*c_A/(c_I*c_e)) >= 0`, independently of Th. Heavy translational
//! temperature factors cancel because one heavy atom becomes one heavy ion.

use super::{
    ElectronImpactRateV1, ElementaryReactionV1, ForwardRateV1, IdealGasKineticsV1, KineticsErrorV1,
    net_rate, rate_from_log,
};
use crate::electron_collision::{
    AVOGADRO_PER_MOL, BOLTZMANN_J_PER_K, ELECTRON_MASS_KG, MaxwellianRateWeightingV1,
};
use crate::equilibrium::composition::IdealGasEquilibriumV1;
use crate::{ReactionId, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_math::det::{ln, sqrt};
use fs_qty::Temperature;

const PLANCK_J_S: f64 = 6.626_070_15e-34;
/// Maximum dimensionless discrepancy between the source reaction entropy
/// constant and declared ground-state statistics. This is a numerical source
/// admission tolerance, not an empirical uncertainty or a fitted correction.
pub const MAX_STATISTICAL_ENTROPY_DISCREPANCY_V1: f64 = 1e-10;

/// Explicit ground-state level degeneracies. Electron spin contributes exactly
/// two states; it is never counted again in either heavy statistical weight.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AtomicStatisticalWeightsV1 {
    /// Neutral atomic ground-state statistical weight, greater than zero.
    pub neutral: u32,
    /// Singly charged ionic ground-state statistical weight, greater than zero.
    pub ion: u32,
}

/// Original source admission, trial evaluation or arithmetic refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum AtomicIonizationErrorV1 {
    /// Invalid state, unsupported atomic assumptions or inconsistent sources.
    Invalid(&'static str),
    /// Original exact reaction, mass, collision or species source refusal.
    Kinetics(KineticsErrorV1),
    /// Source entropy constants disagree with explicit ground-state statistics.
    StatisticalEntropyMismatch {
        /// Original delta-a8 minus the statistical value, dimensionless.
        discrepancy: f64,
    },
    /// Cancellation before a complete result exists.
    Cancelled,
}
impl core::fmt::Display for AtomicIonizationErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ground-state atomic ionization: {self:?}")
    }
}
impl std::error::Error for AtomicIonizationErrorV1 {}
impl From<KineticsErrorV1> for AtomicIonizationErrorV1 {
    fn from(error: KineticsErrorV1) -> Self {
        match error {
            KineticsErrorV1::Cancelled => Self::Cancelled,
            error => Self::Kinetics(error),
        }
    }
}

/// Source-bound atomic channel; immutable, allocation-free trial evaluation.
#[derive(Clone, Debug, PartialEq)]
pub struct AtomicIonizationV1 {
    mechanism: IdealGasKineticsV1,
    indices: [usize; 3],
    statistical_weights: AtomicStatisticalWeightsV1,
    formation_energies_j_mol: [f64; 3],
    ionization_energy_j_mol: f64,
    log_concentration_factor_at_one_kelvin: f64,
    statistical_entropy_discrepancy: f64,
}
impl AtomicIonizationV1 {
    /// Bind the original species bank and electron-impact source. Exact element,
    /// charge, molecular mass, electron identity and binary-order admission are
    /// delegated to the existing kinetics owner. No source values are adjusted.
    ///
    /// Active calorics must be ground-state monatomic on every retained region.
    /// Region-varying a7/a8 and opt-in caloric-reference transformations refuse:
    /// this first atomic closure does not silently reinterpret raw coefficients.
    /// The cross-section threshold may exceed the source energy gap (an explicit
    /// approximate cross section); it may not lie below it except roundoff.
    pub fn new(
        thermodynamics: IdealGasEquilibriumV1,
        reaction: ReactionId,
        ion: SpeciesId,
        forward: ElectronImpactRateV1,
        statistical_weights: AtomicStatisticalWeightsV1,
        relative_mass_tolerance: f64,
    ) -> Result<Self, AtomicIonizationErrorV1> {
        if statistical_weights.neutral == 0 || statistical_weights.ion == 0 {
            return Err(AtomicIonizationErrorV1::Invalid(
                "positive ground-state statistical weights required",
            ));
        }
        let neutral = forward.cross_section().target_species().clone();
        let electron = forward.electron_species().clone();
        let channel = ElementaryReactionV1::new(
            reaction,
            vec![(electron.clone(), 1), (neutral.clone(), 1)],
            vec![(ion.clone(), 1), (electron.clone(), 2)],
            forward,
        )?;
        let mechanism = IdealGasKineticsV1::new(
            thermodynamics,
            "atomic-ground-state-ionization",
            1,
            vec![channel],
            relative_mass_tolerance,
        )?;
        let bank = mechanism.thermodynamics();
        let species = bank.elemental_matrix().species();
        let indices = [neutral, ion, electron].map(|id| {
            species
                .binary_search(&id)
                .expect("admitted reaction species")
        });
        let mut nuclei = 0;
        for row in 0..bank.elemental_matrix().row_count() {
            match bank.elemental_matrix().get(row, indices[0]) {
                Some(0) => {}
                Some(1) => nuclei += 1,
                _ => {
                    return Err(AtomicIonizationErrorV1::Invalid(
                        "ground-state channel requires one monatomic nucleus",
                    ));
                }
            }
        }
        if nuclei != 1 {
            return Err(AtomicIonizationErrorV1::Invalid(
                "ground-state channel requires one monatomic nucleus",
            ));
        }
        let mut formation = [0.0; 3];
        let mut entropy = [0.0; 3];
        for (slot, &index) in indices.iter().enumerate() {
            let model = &bank.models()[index];
            if model.caloric_reference().is_some() {
                return Err(AtomicIonizationErrorV1::Invalid(
                    "atomic channel requires untransformed fixed level references",
                ));
            }
            let first = model.regions()[0].coefficients();
            for region in model.regions() {
                let c = region.coefficients();
                if c[2] != 2.5
                    || [0, 1, 3, 4, 5, 6].iter().any(|&i| c[i] != 0.0)
                    || c[7] != first[7]
                    || c[8] != first[8]
                {
                    return Err(AtomicIonizationErrorV1::Invalid(
                        "fixed monatomic Cp/R=5/2 and ground-state a7/a8 required",
                    ));
                }
            }
            formation[slot] = R * first[7];
            entropy[slot] = first[8];
        }
        let energy_gap = (formation[1] - formation[0]) + formation[2];
        if !super::positive(energy_gap) {
            return Err(AtomicIonizationErrorV1::Invalid(
                "positive finite source ionization energy gap required",
            ));
        }
        let impact = match mechanism.reactions()[0].forward_rate() {
            ForwardRateV1::ElectronImpact(impact) => impact,
            _ => unreachable!("constructed source impact channel"),
        };
        let threshold = impact.cross_section().threshold().value() * AVOGADRO_PER_MOL;
        if threshold < energy_gap && (energy_gap - threshold) / energy_gap > 64.0 * f64::EPSILON {
            return Err(AtomicIonizationErrorV1::Invalid(
                "cross section supports collisions below source ionization energy",
            ));
        }
        let masses = indices.map(|i| bank.models()[i].molar_mass().value() / AVOGADRO_PER_MOL);
        if sqrt(ELECTRON_MASS_KG / masses[1]) > impact.maximum_target_thermal_speed_ratio() {
            return Err(AtomicIonizationErrorV1::Invalid(
                "ionic recoil exceeds massive-target domain",
            ));
        }
        let entropy_difference = (entropy[1] - entropy[0]) + entropy[2];
        let pressure = bank.convention().reference_pressure().value();
        let log_factor = entropy_difference - 2.5 + ln(pressure) - ln(R);
        let statistical_factor = ln(2.0) + ln(f64::from(statistical_weights.ion))
            - ln(f64::from(statistical_weights.neutral))
            + 1.5
                * (ln(2.0 * core::f64::consts::PI * BOLTZMANN_J_PER_K / (PLANCK_J_S * PLANCK_J_S))
                    + ln(masses[1])
                    + ln(masses[2])
                    - ln(masses[0]))
            - ln(AVOGADRO_PER_MOL);
        let discrepancy = log_factor - statistical_factor;
        if !discrepancy.is_finite() || discrepancy.abs() > MAX_STATISTICAL_ENTROPY_DISCREPANCY_V1 {
            return Err(AtomicIonizationErrorV1::StatisticalEntropyMismatch { discrepancy });
        }
        Ok(Self {
            mechanism,
            indices,
            statistical_weights,
            formation_energies_j_mol: formation,
            ionization_energy_j_mol: energy_gap,
            log_concentration_factor_at_one_kelvin: log_factor,
            statistical_entropy_discrepancy: discrepancy,
        })
    }

    /// Original source bank, including every inert species and source card.
    pub const fn thermodynamics(&self) -> &IdealGasEquilibriumV1 {
        self.mechanism.thermodynamics()
    }
    /// Existing exact-conservation one-channel owner, for common-T comparison.
    pub const fn common_temperature_mechanism(&self) -> &IdealGasKineticsV1 {
        &self.mechanism
    }
    /// Canonical source-bank index of neutral A.
    pub const fn neutral_index(&self) -> usize {
        self.indices[0]
    }
    /// Canonical source-bank index of singly charged I.
    pub const fn ion_index(&self) -> usize {
        self.indices[1]
    }
    /// Canonical source-bank index of the structural electron.
    pub const fn electron_index(&self) -> usize {
        self.indices[2]
    }
    /// Original electron-impact cross section and numerical controls.
    pub fn forward_rate(&self) -> &ElectronImpactRateV1 {
        match self.mechanism.reactions()[0].forward_rate() {
            ForwardRateV1::ElectronImpact(impact) => impact,
            _ => unreachable!("immutable source impact channel"),
        }
    }
    /// Ground-state statistical weights used to check source entropy, unmodified.
    pub const fn statistical_weights(&self) -> AtomicStatisticalWeightsV1 {
        self.statistical_weights
    }
    /// Source statistical entropy discrepancy, dimensionless; never corrected.
    pub const fn statistical_entropy_discrepancy(&self) -> f64 {
        self.statistical_entropy_discrepancy
    }
    /// Source ionization energy chi=epsilon_I+epsilon_e-epsilon_A [J/mol].
    pub const fn ionization_energy_j_mol(&self) -> f64 {
        self.ionization_energy_j_mol
    }
    /// Heavy source formation-energy gain epsilon_I-epsilon_A [J/mol].
    /// Unlike chi this retains any authored electron/charge reference offset.
    pub fn heavy_energy_gap_j_mol(&self) -> f64 {
        self.formation_energies_j_mol[1] - self.formation_energies_j_mol[0]
    }
    /// Original [neutral, ion, electron] ground-state energies [J/mol].
    pub const fn formation_energies_j_mol(&self) -> [f64; 3] {
        self.formation_energies_j_mol
    }

    /// ln(Kc) for concentrations normalized by1mol/m³. Source support is
    /// checked at Te for the electron only. Heavy species are at Th and are
    /// not extrapolated to Te; their fixed level constants are already admitted.
    pub fn log_equilibrium_constant(
        &self,
        electron_temperature: Temperature,
    ) -> Result<f64, AtomicIonizationErrorV1> {
        self.admit_source_temperature(self.indices[2], electron_temperature)?;
        let t = electron_temperature.value();
        let log = self.log_concentration_factor_at_one_kelvin + 1.5 * ln(t)
            - (self.ionization_energy_j_mol / R) / t;
        if !log.is_finite() {
            return Err(AtomicIonizationErrorV1::Invalid(
                "unrepresentable atomic concentration factor",
            ));
        }
        Ok(log)
    }

    /// Allocation-free rates and analytic derivatives at actual Th,Te.
    /// Input order is `[neutral, ion, electron]`, in mol/m³, independently of
    /// the bank's canonical order. Exact absence is allowed; no seed is inserted.
    pub fn evaluate(
        &self,
        heavy_temperature: Temperature,
        electron_temperature: Temperature,
        concentrations: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AtomicIonizationEvaluationV1, AtomicIonizationErrorV1> {
        super::poll(&mut cancelled)?;
        if concentrations.iter().any(|c| !c.is_finite() || *c < 0.0) {
            return Err(AtomicIonizationErrorV1::Invalid(
                "nonnegative finite atomic concentrations required",
            ));
        }
        self.admit_source_temperature(self.indices[0], heavy_temperature)?;
        self.admit_source_temperature(self.indices[1], heavy_temperature)?;
        let log_kc = self.log_equilibrium_constant(electron_temperature)?;
        let te = electron_temperature.value();
        let th = heavy_temperature.value();
        let impact = self.forward_rate();
        for &index in &self.indices[..2] {
            let mass =
                self.thermodynamics().models()[index].molar_mass().value() / AVOGADRO_PER_MOL;
            let ratio = (sqrt(th) / sqrt(te)) * (sqrt(ELECTRON_MASS_KG) / sqrt(mass));
            if !ratio.is_finite() || ratio > impact.maximum_target_thermal_speed_ratio() {
                return Err(AtomicIonizationErrorV1::Invalid(
                    "heavy thermal motion exceeds electron-impact source domain",
                ));
            }
        }
        let collision = impact
            .cross_section()
            .evaluate_properties(
                electron_temperature,
                MaxwellianRateWeightingV1::Event,
                impact.controls(),
                &mut cancelled,
            )
            .map_err(|error| match error {
                crate::electron_collision::ElectronCollisionErrorV1::Cancelled => {
                    KineticsErrorV1::Cancelled
                }
                error => KineticsErrorV1::ElectronCollision(error),
            })?;
        let log_kf = collision.log_coefficient_m3_s() + ln(AVOGADRO_PER_MOL);
        let log_kr = log_kf - log_kc;
        let logs = concentrations.map(|c| if c == 0.0 { f64::NEG_INFINITY } else { ln(c) });
        let forward = rate_from_log(log_kf + logs[0] + logs[2])?;
        let reverse = rate_from_log(log_kr + logs[1] + 2.0 * logs[2])?;
        let affinity = match (
            concentrations[0] == 0.0,
            concentrations[1] == 0.0 || concentrations[2] == 0.0,
        ) {
            (true, true) => None,
            (true, false) => Some(f64::NEG_INFINITY),
            (false, true) => Some(f64::INFINITY),
            (false, false) => Some(log_kc + logs[0] - logs[1] - logs[2]),
        };
        let net = net_rate(forward, reverse, affinity.unwrap_or(0.0))?;
        let concentration_derivatives = [
            rate_from_log(log_kf + logs[2])?,
            -rate_from_log(log_kr + 2.0 * logs[2])?,
            if concentrations[2] > 0.0 {
                (net - reverse) / concentrations[2]
            } else {
                rate_from_log(log_kf + logs[0])?
            },
        ];
        let dlog_kc = (1.5 + (self.ionization_energy_j_mol / R) / te) / te;
        let temperature_derivative =
            net * collision.logarithmic_temperature_derivative_per_k() + reverse * dlog_kc;
        let force = affinity.unwrap_or(0.0);
        let entropy_production = if net == 0.0 { 0.0 } else { R * (net * force) };
        if concentration_derivatives.iter().any(|d| !d.is_finite())
            || !temperature_derivative.is_finite()
            || entropy_production.is_nan()
            || entropy_production < 0.0
            || (force.is_finite() && !entropy_production.is_finite())
        {
            return Err(KineticsErrorV1::Arithmetic(
                "atomic rate derivative or entropy production",
            )
            .into());
        }
        super::poll(&mut cancelled)?;
        Ok(AtomicIonizationEvaluationV1 {
            progress_rate_mol_m3_s: net,
            forward_rate_mol_m3_s: forward,
            reverse_rate_mol_m3_s: reverse,
            concentration_derivatives_per_s: concentration_derivatives,
            electron_temperature_derivative_mol_m3_s_k: temperature_derivative,
            log_equilibrium_constant: log_kc,
            dimensionless_affinity: affinity,
            entropy_production_w_m3_k: entropy_production,
        })
    }

    fn admit_source_temperature(
        &self,
        index: usize,
        temperature: Temperature,
    ) -> Result<(), AtomicIonizationErrorV1> {
        let model = &self.thermodynamics().models()[index];
        model.evaluate_properties(temperature).map_err(|source| {
            KineticsErrorV1::SpeciesEvaluation {
                species: model.species().clone(),
                source,
            }
        })?;
        Ok(())
    }
}

/// Complete trial rates in coherent SI. Derivatives hold all other supplied
/// concentrations/temperatures fixed; the L3 energy owner adds population dilution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AtomicIonizationEvaluationV1 {
    /// Signed forward-minus-reverse progress [mol/(m³ s)].
    pub progress_rate_mol_m3_s: f64,
    /// Forward ionization progress [mol/(m³ s)].
    pub forward_rate_mol_m3_s: f64,
    /// Reverse three-body ground-state recombination progress [mol/(m³ s)].
    pub reverse_rate_mol_m3_s: f64,
    /// dq/d[c_A,c_I,c_e], each in 1/s, in the explicit input order.
    pub concentration_derivatives_per_s: [f64; 3],
    /// dq/dTe at fixed concentrations and Th [mol/(m³ s K)].
    pub electron_temperature_derivative_mol_m3_s_k: f64,
    /// Natural logarithm of the source monatomic Kc, normalized by1mol/m³.
    pub log_equilibrium_constant: f64,
    /// Dimensionless two-pool chemical entropy force. One-sided absence gives
    /// signed infinity; competing absent products/reactants give None.
    pub dimensionless_affinity: Option<f64>,
    /// R*q*affinity [W/(m³ K)], nonnegative. Positive infinity denotes an exact
    /// ideal-mixture boundary derivative, not a finite-step entropy increment.
    pub entropy_production_w_m3_k: f64,
}

#[cfg(test)]
mod tests;
