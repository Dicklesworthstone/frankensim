//! Primary creep with a separate resistance memory that recovers during dwell.
//!
//! The permanent path still obeys `alpha_dot = rate_Norton(q,T) h(r)`, while
//! `r_dot = alpha_dot - k(T) r`. The explicitly selected primary multiplier is
//! `h(r) = f + (1-f) (1+r/a)^(-m)`. Recovery changes subsequent creep resistance;
//! it never reverses signed permanent strain or decreases accumulated alpha.
//!
//! With `c = delta_alpha` and `g = 1 + k(T) dt`, backward Euler gives
//! `r_new = (r_old+c)/g`. Thus the original primary return is used exactly with
//! reference strain `a*g` and old kinetic memory `r_old`. Its consistent tangent
//! includes this coupled resistance evolution without an additional root solve.
//! Only the private prepared law carries the transformed parameters.
//!
//! This is an authored phenomenological resistance model requiring independent
//! recovery calibration. Resistance carries no stored energy in this model;
//! the original native work account and physical creep heat `q*c` are retained.
//! No viscoelastic strain recovery, tertiary damage, rupture or finite-strain
//! behavior is implied.

/// Simultaneous resistance recovery and native energetic isotropic J2 plasticity.
pub mod plastic;
/// Immutable source binding for independently selected recovery coefficients.
pub mod source;
/// Exact free-sided axial restriction with the original resistance history.
pub mod uniaxial;

use super::thermal::expansion::CreepHistoryStrain;
use super::{
    CreepIncrement, CreepLaw, CreepState, CreepTolerances, GAS_CONSTANT_J_MOL_K, NortonCreep,
    Result, Update, positive, require,
};
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard, ValidityDomain};
use fs_math::det;

/// Independently calibrated first-order recovery kinetics. The original Norton
/// temperature interval also bounds this model and its recovery reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CreepRecoveryParameters {
    /// Nonnegative resistance recovery rate at the reference temperature [1/s].
    /// Zero disables recovery exactly; it does not disable permanent creep.
    pub reference_rate_s: f64,
    /// Nonnegative recovery activation energy [J/mol], independent of creep Q.
    pub activation_energy_j_mol: f64,
    /// Absolute temperature of the recovery reference rate [K].
    pub reference_temperature_k: f64,
}

/// Permanent creep and kinetic resistance are distinct histories. For states
/// born from the virgin material, resistance is nonnegative and no larger than
/// `native.equivalent_creep_strain`.
#[derive(Clone, Debug, PartialEq)]
pub struct RecoveringCreepState {
    pub native: CreepState,
    /// Recoverable primary-resistance memory, not a signed strain or path length.
    pub resistance_strain: f64,
}

impl CreepHistoryStrain for RecoveringCreepState {
    const STRAIN_FIELDS: usize = 1;

    fn inelastic_strain(&self) -> Voigt {
        self.native.creep_strain
    }
}

// One preparation and resistance evolution for pure creep and simultaneous J2.
// The transformed Norton law is returned to the original mechanical owner;
// recovery never supplies another constitutive return or stored-energy model.
#[derive(Clone, Copy, Debug, PartialEq)]
struct RecoveryStep {
    parameters: CreepRecoveryParameters,
    rate_s: f64,
    divisor: f64,
}

impl RecoveryStep {
    fn prepare(
        native: &NortonCreep,
        recovery: CreepRecoveryParameters,
    ) -> Result<(NortonCreep, Self)> {
        let mut primary = *native.strain_hardening().ok_or(super::CreepError::Invalid(
            "resistance recovery requires explicitly selected primary creep coefficients",
        ))?;
        let [lo, hi] = native.parameters().temperature_range_k;
        require(
            recovery.reference_rate_s.is_finite()
                && recovery.reference_rate_s >= 0.0
                && recovery.activation_energy_j_mol.is_finite()
                && recovery.activation_energy_j_mol >= 0.0
                && recovery.reference_temperature_k.is_finite()
                && recovery.reference_temperature_k >= lo
                && recovery.reference_temperature_k <= hi,
            "recovery needs nonnegative finite rate/activation energy and a reference temperature in the native domain",
        )?;
        let recovery_rate_s = if recovery.reference_rate_s == 0.0 {
            0.0
        } else if native.temperature_k() == recovery.reference_temperature_k
            || recovery.activation_energy_j_mol == 0.0
        {
            recovery.reference_rate_s
        } else {
            let shift = (recovery.activation_energy_j_mol / GAS_CONSTANT_J_MOL_K)
                * ((native.temperature_k() - recovery.reference_temperature_k)
                    / recovery.reference_temperature_k
                    / native.temperature_k());
            let log_rate = det::ln(recovery.reference_rate_s) + shift;
            require(
                log_rate.is_finite(),
                "recovery Arrhenius shift is not representable",
            )?;
            let rate = det::exp(log_rate);
            require(positive(rate), "recovery rate is not representable")?;
            rate
        };
        let recovery_step = recovery_rate_s * native.duration_s();
        let recovery_divisor = 1.0 + recovery_step;
        require(
            recovery_step.is_finite() && positive(recovery_divisor),
            "recovery interval is not representable",
        )?;
        if primary.active() {
            primary.reference_strain *= recovery_divisor;
            require(
                positive(primary.reference_strain),
                "recovered primary strain scale is not representable",
            )?;
        }
        let prepared = native.clone().with_strain_hardening(primary)?;
        Ok((
            prepared,
            Self {
                parameters: recovery,
                rate_s: recovery_rate_s,
                divisor: recovery_divisor,
            },
        ))
    }

    fn next_resistance(&self, old: f64, delta: f64, alpha: f64) -> Result<f64> {
        let numerator = old + delta;
        let resistance = numerator / self.divisor;
        require(
            numerator.is_finite() && (numerator == 0.0 || positive(resistance)),
            "recovered resistance history is not representable",
        )?;
        validate_resistance(resistance, alpha)?;
        Ok(resistance)
    }
}

/// A fixed-temperature, fixed-duration native Norton interval with recovering
/// primary resistance. Rebuild from the original parameters for each interval.
#[derive(Clone, Debug, PartialEq)]
pub struct RecoveringNortonCreep {
    native: NortonCreep,
    prepared: NortonCreep,
    recovery: RecoveryStep,
}

impl RecoveringNortonCreep {
    pub fn new(native: NortonCreep, recovery: CreepRecoveryParameters) -> Result<Self> {
        let (prepared, recovery) = RecoveryStep::prepare(&native, recovery)?;
        Ok(Self {
            native,
            prepared,
            recovery,
        })
    }

    /// Original immutable Norton law, including the authored primary selection.
    /// Timestep-transformed coefficients are never exposed as material inputs.
    pub const fn norton(&self) -> &NortonCreep {
        &self.native
    }

    pub const fn recovery_parameters(&self) -> &CreepRecoveryParameters {
        &self.recovery.parameters
    }

    /// Independently shifted recovery rate at this interval's temperature [1/s].
    pub const fn recovery_rate_s(&self) -> f64 {
        self.recovery.rate_s
    }

    pub(super) const fn prepared_norton(&self) -> &NortonCreep {
        &self.prepared
    }

    /// Validate the original permanent history and its independent resistance.
    pub fn validate_history(&self, strain: &Voigt, state: &RecoveringCreepState) -> Result<()> {
        self.native.validate_history(strain, &state.native)?;
        validate_resistance(
            state.resistance_strain,
            state.native.equivalent_creep_strain,
        )
    }

    pub(super) fn next_resistance(&self, old: f64, delta: f64, alpha: f64) -> Result<f64> {
        self.recovery.next_resistance(old, delta, alpha)
    }

    fn update(&self, strain: &Voigt, old: &RecoveringCreepState) -> Result<(Update, f64)> {
        validate_resistance(old.resistance_strain, old.native.equivalent_creep_strain)?;
        let next = self
            .prepared
            .update_with_memory(strain, &old.native, old.resistance_strain)?;
        let resistance = self.next_resistance(
            old.resistance_strain,
            next.delta_alpha,
            next.state.equivalent_creep_strain,
        )?;
        Ok((next, resistance))
    }

    /// Return one candidate using the original native energy account and actual
    /// solved creep increment. Recovery of kinetic resistance creates no heat.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old: &RecoveringCreepState,
        strain: &Voigt,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<RecoveringCreepState>> {
        validate_resistance(old.resistance_strain, old.native.equivalent_creep_strain)?;
        let (next, delta) = self.prepared.advance_with_memory(
            old_strain,
            &old.native,
            strain,
            tolerances,
            old.resistance_strain,
        )?;
        let resistance_strain = self.next_resistance(
            old.resistance_strain,
            delta,
            next.state.equivalent_creep_strain,
        )?;
        Ok(CreepIncrement {
            strain: next.strain,
            stress_pa: next.stress_pa,
            tangent_pa: next.tangent_pa,
            state: RecoveringCreepState {
                native: next.state,
                resistance_strain,
            },
            energy: next.energy,
        })
    }
}

fn validate_resistance(resistance: f64, alpha: f64) -> Result<()> {
    require(
        resistance.is_finite() && resistance >= 0.0 && alpha.is_finite() && alpha >= 0.0,
        "resistance and accumulated creep must be finite and nonnegative",
    )?;
    // Match the native path-consistency allowance for independently accumulated
    // floating histories. A relative bound has no absolute floor at virgin alpha;
    // compare differences so alpha + allowance cannot overflow.
    let allowance = 256.0 * f64::EPSILON * alpha.max(resistance);
    require(
        resistance <= alpha || resistance - alpha <= allowance,
        "recovered resistance cannot exceed the accumulated permanent creep path",
    )
}

impl CreepLaw for RecoveringNortonCreep {
    fn kinetics(&self) -> &NortonCreep {
        self.norton()
    }

    fn validate_creep_history(&self, strain: &Voigt, state: &RecoveringCreepState) -> Result<()> {
        self.validate_history(strain, state)
    }

    fn creep_increment(
        &self,
        old_strain: &Voigt,
        old: &RecoveringCreepState,
        strain: &Voigt,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<RecoveringCreepState>> {
        self.advance_with_energy(old_strain, old, strain, tolerances)
    }
}

impl SmallStrainLaw for RecoveringNortonCreep {
    type State = RecoveringCreepState;

    fn initial_state(&self) -> Self::State {
        RecoveringCreepState {
            native: self.native.initial_state(),
            resistance_strain: 0.0,
        }
    }

    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt {
        self.update(strain, state)
            .map(|(v, _)| v.stress)
            .unwrap_or([f64::NAN; 6])
    }

    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6 {
        self.update(strain, state)
            .map(|(v, _)| v.tangent)
            .unwrap_or([[f64::NAN; 6]; 6])
    }

    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State {
        self.update(strain, state)
            .map(|(v, resistance_strain)| RecoveringCreepState {
                native: v.state,
                resistance_strain,
            })
            .unwrap_or(RecoveringCreepState {
                native: CreepState {
                    creep_strain: [f64::NAN; 6],
                    equivalent_creep_strain: f64::NAN,
                },
                resistance_strain: f64::NAN,
            })
    }

    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true,
            polyconvex: None,
            tangent_symmetric: true,
            failure_envelope: "bounded small-strain primary/secondary creep with kinetic resistance recovery; no damage or rupture law",
        }
    }

    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.recovering-primary-norton-creep",
            "0.1.0",
            Ambition::Solid,
            vec![
                "isotropic small-strain Norton creep with independent primary resistance memory".into(),
                "first-order resistance recovery with independently selected Arrhenius kinetics".into(),
                "backward Euler; permanent accumulated creep never decreases".into(),
                "native elastic storage and creep heat; no stored energy assigned to resistance".into(),
            ],
            ValidityDomain::unconstrained()
                .with("strain-magnitude", 0.0, self.native.elastic().strain_limit)
                .with(
                    "temperature-k",
                    self.native.parameters().temperature_range_k[0],
                    self.native.parameters().temperature_range_k[1],
                ),
            vec![
                "no viscoelastic strain recovery, plasticity, finite deformation, tertiary damage or rupture".into(),
                "independent recovery calibration required; no measured material qualification without source binding".into(),
            ],
            0.05,
        )
    }
}
