//! Finite creep with independent resistance recovery and the original thermal law.
//!
//! The base family owns affine Kirchhoff coefficients, multiplicative thermal
//! stretch, Norton viscosity, free energy, entropy and thermal-capacity bounds.
//! This extension selects independent recovery kinetics and retains both the
//! permanent finite history and resistance. It supports standalone Norton creep
//! and simultaneous Norton/J2, including the base's energetic hardening.
//!
//! At an endpoint temperature, the original finite recovery owner solves
//! r_new=(r_old+delta_creep)/(1+k(T)*dt), with its original inverse-rate potential
//! and analytic tangent. Resistance has no stored energy or entropy. Recovery
//! never erases the permanent metric or accumulated path and creates no heat.
//! Physical creep work, plastic yield heat, reversible thermal exchange and the
//! signed endpoint-work defect remain distinct in the shared thermal account.

use super::{
    AffineThermalFiniteJ2, FiniteJ2Response, FiniteJ2ThermalError, FiniteJ2ThermalIncrement,
    FiniteJ2Thermodynamics, FiniteThermalInelasticMaterial, FiniteThermalInelasticPoint,
    FrozenNative, FrozenThermalFiniteJ2,
};
use crate::plastic::finite::creep::FiniteNortonCreep;
use crate::plastic::finite::creep::recovery::{
    FiniteRecoveringCreepState, FiniteRecoveringNortonCreep, FiniteRecoveringNortonJ2,
    validate_resistance,
};
use fs_math::det;

/// Independent resistance-recovery activation, sharing the base's reference T.
/// These coefficients are kinetic; neither contributes stored energy or entropy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteCreepThermalRecovery {
    /// Nonnegative recovery rate at the base's reference temperature [1/s].
    /// Zero disables recovery exactly, while permanent creep remains active.
    pub reference_recovery_rate_s: f64,
    /// Nonnegative Q_recovery/R in kelvin, independently calibrated from creep.
    /// k(T)=k_ref*exp[-(Q_recovery/R)*(1/T-1/T_ref)].
    pub activation_temperature_k: f64,
}

/// Complete thermal family with independently retained recovering resistance.
#[derive(Clone, Debug, PartialEq)]
pub struct AffineThermalRecoveringCreep {
    base: AffineThermalFiniteJ2,
    recovery: FiniteCreepThermalRecovery,
}

#[derive(Clone, Debug, PartialEq)]
enum FrozenRecovery {
    Standalone(FiniteRecoveringNortonCreep),
    Combined(FiniteRecoveringNortonJ2),
}

/// One owned temperature view. Inspection has no constitutive interval and
/// cannot advance a history; a physical step retains its complete recovery law.
#[derive(Clone, Debug, PartialEq)]
pub struct FrozenThermalRecoveringCreep {
    thermal: FrozenThermalFiniteJ2,
    native: Option<FrozenRecovery>,
    combined: bool,
}

impl AffineThermalRecoveringCreep {
    /// Select an explicit primary Norton family and independent recovery.
    /// Rate-independent plasticity and Perzyna-only bases refuse. A disabled
    /// primary multiplier remains admissible and preserves secondary mechanics.
    pub fn new(
        base: AffineThermalFiniteJ2,
        recovery: FiniteCreepThermalRecovery,
    ) -> Result<Self, FiniteJ2ThermalError> {
        if (!base.is_norton_creep() && base.maximum_accumulated_creep_strain().is_none())
            || base.norton_strain_hardening().is_none()
        {
            return Err(FiniteJ2ThermalError::Parameters(
                "finite thermal recovery requires explicitly selected standalone or combined primary Norton creep",
            ));
        }
        if !recovery.reference_recovery_rate_s.is_finite()
            || recovery.reference_recovery_rate_s < 0.0
            || !recovery.activation_temperature_k.is_finite()
            || recovery.activation_temperature_k < 0.0
        {
            return Err(FiniteJ2ThermalError::Parameters(
                "finite thermal recovery requires nonnegative finite rate and activation temperature",
            ));
        }
        let family = Self { base, recovery };
        // Arrhenius recovery is monotone on the admitted positive-T interval.
        // Endpoints bound its representability without inventing any duration.
        for temperature in family.temperature_range_k() {
            family.recovery_rate_s(temperature)?;
        }
        Ok(family)
    }

    /// Original authored thermal coefficients, kinetics and applicability bounds.
    pub const fn base(&self) -> &AffineThermalFiniteJ2 {
        &self.base
    }

    pub const fn recovery(&self) -> &FiniteCreepThermalRecovery {
        &self.recovery
    }

    pub fn reference_temperature_k(&self) -> f64 {
        self.base.parameters().reference_temperature_k
    }

    pub fn temperature_range_k(&self) -> [f64; 2] {
        self.base.parameters().temperature_range_k
    }

    /// Recovery is purely kinetic, so the original thermal-capacity bound holds.
    pub fn minimum_mechanical_heat_capacity_j_m3_k(&self) -> f64 {
        self.base.minimum_mechanical_heat_capacity_j_m3_k()
    }

    pub fn initial_state(&self) -> FiniteRecoveringCreepState {
        FiniteRecoveringCreepState {
            native: self.base.initial_state(),
            resistance_strain: 0.0,
        }
    }

    /// Recovery rate at one admitted temperature, without a constitutive step.
    /// The logarithmic product retains representable rates even when a bare
    /// Arrhenius multiplier would overflow or underflow.
    pub fn recovery_rate_s(&self, temperature_k: f64) -> Result<f64, FiniteJ2ThermalError> {
        self.base.admit_temperature(temperature_k)?;
        let reference = self.reference_temperature_k();
        let parameters = self.recovery;
        if parameters.reference_recovery_rate_s == 0.0
            || temperature_k == reference
            || parameters.activation_temperature_k == 0.0
        {
            return Ok(parameters.reference_recovery_rate_s);
        }
        // Avoid subtracting nearby reciprocals or forming T*T_ref.
        let shift = ((temperature_k - reference) / temperature_k)
            * (parameters.activation_temperature_k / reference);
        let log_rate = det::ln(parameters.reference_recovery_rate_s) + shift;
        let rate = det::exp(log_rate);
        if !shift.is_finite() || !log_rate.is_finite() || !rate.is_finite() || rate <= 0.0 {
            return Err(FiniteJ2ThermalError::Numerical(
                "Arrhenius recovery rate must remain representable and positive over its temperature domain",
            ));
        }
        Ok(rate)
    }

    /// Inspect complete history and thermodynamics without flow or recovery.
    pub fn at_temperature_for_inspection(
        &self,
        temperature_k: f64,
    ) -> Result<FrozenThermalRecoveringCreep, FiniteJ2ThermalError> {
        Ok(FrozenThermalRecoveringCreep {
            thermal: self.base.at_temperature_for_inspection(temperature_k)?,
            native: None,
            combined: self.base.maximum_accumulated_creep_strain().is_some(),
        })
    }

    /// Freeze both physical kinetics at the same endpoint temperature and actual
    /// interval. Temperature trials own no specimen history and cannot age it.
    pub fn at_temperature_for_step(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<FrozenThermalRecoveringCreep, FiniteJ2ThermalError> {
        let thermal = self
            .base
            .at_temperature_for_step(temperature_k, duration_s)?;
        let rate = self.recovery_rate_s(temperature_k)?;
        let native = match &thermal.native {
            FrozenNative::Overstress(law) if self.base.is_norton_creep() => {
                // Reuse the frozen original coefficients and primary selection.
                // The native recovery owner retains the one finite return map.
                let p = law.backbone().parameters();
                let primary =
                    *law.norton_strain_hardening()
                        .ok_or(FiniteJ2ThermalError::Parameters(
                            "finite thermal recovery lost its primary law",
                        ))?;
                let norton = FiniteNortonCreep::new(
                    p.shear_modulus_pa,
                    p.bulk_modulus_pa,
                    *law.backbone().limits(),
                    *law.parameters(),
                    duration_s,
                )?
                .with_strain_hardening(primary)?;
                FrozenRecovery::Standalone(FiniteRecoveringNortonCreep::new(norton, rate)?)
            }
            FrozenNative::Combined(law) => {
                FrozenRecovery::Combined(FiniteRecoveringNortonJ2::new(law.clone(), rate)?)
            }
            _ => {
                return Err(FiniteJ2ThermalError::Parameters(
                    "finite thermal recovery requires its complete selected Norton owner",
                ));
            }
        };
        Ok(FrozenThermalRecoveringCreep {
            thermal,
            native: Some(native),
            combined: self.base.maximum_accumulated_creep_strain().is_some(),
        })
    }

    pub fn thermodynamics(
        &self,
        f: &[f64; 9],
        state: &FiniteRecoveringCreepState,
        temperature_k: f64,
    ) -> Result<FiniteJ2Thermodynamics, FiniteJ2ThermalError> {
        self.at_temperature_for_inspection(temperature_k)?
            .thermodynamics(f, state)
    }

    pub fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FiniteRecoveringCreepState,
        temperature_k: f64,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2ThermalError> {
        self.at_temperature_for_inspection(temperature_k)?
            .validate_committed(f, state, yield_tolerance_pa)
    }

    /// One finite return with the original two-temperature thermodynamic account.
    /// Accepted native history is checked at OLD T before it can be advanced;
    /// coefficient changes at NEW T may legitimately activate plastic yielding.
    #[allow(clippy::too_many_arguments)]
    pub fn checked_increment_for_step(
        &self,
        old_f: &[f64; 9],
        old: &FiniteRecoveringCreepState,
        old_temperature_k: f64,
        new_f: &[f64; 9],
        new_temperature_k: f64,
        yield_tolerance_pa: f64,
        duration_s: f64,
    ) -> Result<FiniteJ2ThermalIncrement<FiniteRecoveringCreepState>, FiniteJ2ThermalError> {
        let old_law = self.at_temperature_for_inspection(old_temperature_k)?;
        let new_law = self.at_temperature_for_step(new_temperature_k, duration_s)?;
        old_law.validate_history(old)?;
        AffineThermalFiniteJ2::increment_with_response(
            old_f,
            &old.native,
            new_f,
            yield_tolerance_pa,
            &old_law.thermal,
            &new_law.thermal,
            |fm| new_law.native_response(fm, old),
            |state| &state.native,
        )
    }
}

impl FrozenThermalRecoveringCreep {
    pub fn temperature_k(&self) -> f64 {
        self.thermal.temperature_k()
    }

    pub fn thermal_stretch(&self) -> f64 {
        self.thermal.thermal_stretch()
    }

    pub fn log_thermal_stretch(&self) -> f64 {
        self.thermal.log_thermal_stretch()
    }

    /// None denotes inspection; it never means a zero or unit physical step.
    pub fn duration_s(&self) -> Option<f64> {
        self.native.as_ref().map(|native| match native {
            FrozenRecovery::Standalone(law) => law.norton().duration_s(),
            FrozenRecovery::Combined(law) => law.norton_j2().duration_s(),
        })
    }

    pub fn validate_history(
        &self,
        state: &FiniteRecoveringCreepState,
    ) -> Result<(), FiniteJ2ThermalError> {
        self.thermal.validate_history(&state.native)?;
        let accumulated = if self.combined {
            state.native.creep_alpha
        } else {
            state.native.alpha
        };
        Ok(validate_resistance(state.resistance_strain, accumulated)?)
    }

    fn native_response(
        &self,
        fm: &[f64; 9],
        old: &FiniteRecoveringCreepState,
    ) -> Result<FiniteJ2Response<FiniteRecoveringCreepState>, FiniteJ2ThermalError> {
        Ok(
            match self
                .native
                .as_ref()
                .ok_or(FiniteJ2ThermalError::DurationRequired)?
            {
                FrozenRecovery::Standalone(law) => law.response(fm, old)?,
                FrozenRecovery::Combined(law) => law.response(fm, old)?,
            },
        )
    }

    pub fn response(
        &self,
        f: &[f64; 9],
        old: &FiniteRecoveringCreepState,
    ) -> Result<FiniteJ2Response<FiniteRecoveringCreepState>, FiniteJ2ThermalError> {
        let (fm, j) = self.thermal.mechanical_deformation(f)?;
        self.thermal
            .map_response(self.native_response(&fm, old)?, j)
    }

    pub fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &FiniteRecoveringCreepState,
    ) -> Result<([f64; 9], f64), FiniteJ2ThermalError> {
        let (fm, _) = self.thermal.mechanical_deformation(f)?;
        let (piola, potential) = match self
            .native
            .as_ref()
            .ok_or(FiniteJ2ThermalError::DurationRequired)?
        {
            FrozenRecovery::Standalone(law) => law.stress_and_potential(&fm, old)?,
            FrozenRecovery::Combined(law) => law.stress_and_potential(&fm, old)?,
        };
        let piola = piola.map(|value| value * self.thermal.inverse_thermal_stretch);
        if piola.iter().any(|value| !value.is_finite()) {
            return Err(FiniteJ2ThermalError::Numerical(
                "thermal recovering Piola chain is not representable",
            ));
        }
        Ok((piola, potential))
    }

    pub fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FiniteRecoveringCreepState,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2ThermalError> {
        self.validate_history(state)?;
        self.thermal
            .validate_committed(f, &state.native, yield_tolerance_pa)
    }

    pub fn thermodynamics(
        &self,
        f: &[f64; 9],
        state: &FiniteRecoveringCreepState,
    ) -> Result<FiniteJ2Thermodynamics, FiniteJ2ThermalError> {
        self.validate_history(state)?;
        self.thermal.thermodynamics(f, &state.native)
    }
}

impl FiniteThermalInelasticMaterial for AffineThermalRecoveringCreep {
    type State = FiniteRecoveringCreepState;
    type Frozen = FrozenThermalRecoveringCreep;

    fn initial_state(&self) -> Self::State {
        Self::initial_state(self)
    }
    fn reference_temperature_k(&self) -> f64 {
        Self::reference_temperature_k(self)
    }
    fn temperature_range_k(&self) -> [f64; 2] {
        Self::temperature_range_k(self)
    }
    fn minimum_mechanical_heat_capacity_j_m3_k(&self) -> f64 {
        Self::minimum_mechanical_heat_capacity_j_m3_k(self)
    }
    fn at_temperature_for_inspection(
        &self,
        temperature_k: f64,
    ) -> Result<Self::Frozen, FiniteJ2ThermalError> {
        Self::at_temperature_for_inspection(self, temperature_k)
    }
    fn at_temperature_for_step(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<Self::Frozen, FiniteJ2ThermalError> {
        Self::at_temperature_for_step(self, temperature_k, duration_s)
    }
    fn checked_increment_for_step(
        &self,
        old_f: &[f64; 9],
        old: &Self::State,
        old_temperature_k: f64,
        new_f: &[f64; 9],
        new_temperature_k: f64,
        yield_tolerance_pa: f64,
        duration_s: f64,
    ) -> Result<FiniteJ2ThermalIncrement<Self::State>, FiniteJ2ThermalError> {
        Self::checked_increment_for_step(
            self,
            old_f,
            old,
            old_temperature_k,
            new_f,
            new_temperature_k,
            yield_tolerance_pa,
            duration_s,
        )
    }
}

impl FiniteThermalInelasticPoint for FrozenThermalRecoveringCreep {
    type State = FiniteRecoveringCreepState;

    fn response(
        &self,
        f: &[f64; 9],
        old: &Self::State,
    ) -> Result<FiniteJ2Response<Self::State>, FiniteJ2ThermalError> {
        Self::response(self, f, old)
    }
    fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &Self::State,
    ) -> Result<([f64; 9], f64), FiniteJ2ThermalError> {
        Self::stress_and_potential(self, f, old)
    }
    fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &Self::State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2ThermalError> {
        Self::validate_committed(self, f, state, yield_tolerance_pa)
    }
    fn thermodynamics(
        &self,
        f: &[f64; 9],
        state: &Self::State,
    ) -> Result<FiniteJ2Thermodynamics, FiniteJ2ThermalError> {
        Self::thermodynamics(self, f, state)
    }
}
