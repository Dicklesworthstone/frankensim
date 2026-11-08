//! Independent primary resistance recovery on the native finite Norton owners.
//!
//! At this explicitly fixed temperature, `alpha_dot=rate_Norton(q) h(r)` and
//! `r_dot=alpha_dot-k*r`, with equivalent Kirchhoff stress q and authored
//! `k >= 0` in 1/s. Backward Euler gives `r_new=(r_old+delta_alpha)/(1+k*dt)`.
//! The original primary return therefore uses `a*(1+k*dt)` and old resistance
//! as its kinetic memory. The scalar return, inverse-rate potential, finite
//! metric reconstruction and analytic tangent stay with their original owners.
//!
//! Resistance is independent of the irreversible path and carries no stored
//! energy. Recovery never restores the creep metric, reduces accumulated creep
//! or generates heat; physical creep heat remains `q*delta_alpha`. As in
//! [`super::FiniteNortonCreep`], the native `alpha`, `plastic_increment` and
//! `plastic_dissipation_j_m3` names refer to this sole creep mechanism.
//! For simultaneous creep and J2 yielding, [`FiniteRecoveringNortonJ2`] retains
//! separate permanent plastic `alpha` and creep `creep_alpha`; recovery affects
//! only creep resistance, leaving energetic plastic hardening unchanged.
//!
//! This is the same authored phenomenological resistance model as
//! [`crate::creep::recovery`], now with finite kinematics. Its coefficients
//! require independent calibration. This owner has no automatic temperature
//! shift, source-card binding, damage or rupture law.
//! Consumers must retain the complete typed history, including resistance.

use super::FiniteNortonCreep;
use crate::creep::StrainHardeningCreep;
use crate::plastic::finite::combined::FiniteNortonJ2;
use crate::plastic::finite::{FiniteJ2Error, FiniteJ2Response, FiniteJ2State, ReturnMap};

/// Permanent finite creep metric/path and separately recovering resistance.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteRecoveringCreepState {
    /// Original native history. Sole creep uses alpha for permanent creep;
    /// combined J2 uses alpha for plastic and creep_alpha for permanent creep.
    /// Neither counter is recovering resistance.
    pub native: FiniteJ2State,
    /// Nonnegative primary resistance, at most the accumulated permanent path.
    pub resistance_strain: f64,
}

/// One authored, fixed-temperature and fixed-duration recovery/creep interval.
/// Rebuild from the original coefficients for another physical interval.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteRecoveringNortonCreep {
    native: FiniteNortonCreep,
    prepared: FiniteNortonCreep,
    recovery: RecoveryStep,
}

impl FiniteRecoveringNortonCreep {
    /// Admit an explicitly selected primary law and nonnegative finite recovery
    /// rate [1/s] at the same fixed temperature as the native Norton kinetics.
    /// Zero recovery preserves the original primary law exactly when r=alpha.
    pub fn new(native: FiniteNortonCreep, recovery_rate_s: f64) -> Result<Self, FiniteJ2Error> {
        let (primary, recovery) = RecoveryStep::prepare(
            native.strain_hardening(),
            native.duration_s(),
            recovery_rate_s,
        )?;
        let prepared = native.clone().with_strain_hardening(primary)?;
        Ok(Self {
            native,
            prepared,
            recovery,
        })
    }

    /// Original authored Norton coefficients, limits, primary law and duration.
    /// The private duration-transformed primary scale is not a material input.
    pub const fn norton(&self) -> &FiniteNortonCreep {
        &self.native
    }

    /// Authored recovery rate at this interval's fixed temperature [1/s].
    pub const fn recovery_rate_s(&self) -> f64 {
        self.recovery.rate_s
    }

    /// Virgin metric, zero permanent path and zero primary resistance.
    pub fn initial_state(&self) -> FiniteRecoveringCreepState {
        FiniteRecoveringCreepState {
            native: self.native.initial_state(),
            resistance_strain: 0.0,
        }
    }

    /// Check all original finite-history bounds and independent resistance.
    pub fn validate_history(
        &self,
        state: &FiniteRecoveringCreepState,
    ) -> Result<(), FiniteJ2Error> {
        self.native.validate_history(&state.native)?;
        validate_resistance(state.resistance_strain, state.native.alpha)
    }

    /// Inspect original elastic storage at fixed complete history, without
    /// advancing resistance or consuming another physical creep interval.
    pub fn stored_energy(
        &self,
        f: &[f64; 9],
        state: &FiniteRecoveringCreepState,
    ) -> Result<f64, FiniteJ2Error> {
        self.validate_history(state)?;
        self.native.stored_energy(f, &state.native)
    }

    /// Inspect an accepted state without a return, yield projection or recovery.
    pub fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FiniteRecoveringCreepState,
    ) -> Result<(), FiniteJ2Error> {
        self.validate_history(state)?;
        self.native.validate_committed(f, &state.native)
    }

    /// Original finite return and full analytic dP/dF, now retaining both
    /// irreversible creep and recovering resistance in one candidate state.
    pub fn response(
        &self,
        f: &[f64; 9],
        old: &FiniteRecoveringCreepState,
    ) -> Result<FiniteJ2Response<FiniteRecoveringCreepState>, FiniteJ2Error> {
        let (returned, resistance_strain) = self.return_map(f, old)?;
        self.prepared
            .native
            .backbone()
            .response_from_return(f, &old.native, returned)
            .map(|response| {
                response.map_state(|native| FiniteRecoveringCreepState {
                    native,
                    resistance_strain,
                })
            })
    }

    /// Identical return and complete-history checks, without tangent assembly
    /// for Newton line search. The borrowed accepted state is never changed.
    pub fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &FiniteRecoveringCreepState,
    ) -> Result<([f64; 9], f64), FiniteJ2Error> {
        let (returned, _) = self.return_map(f, old)?;
        Ok((returned.first_piola, returned.potential))
    }

    /// Inverse endpoint kinetic relation at an explicit creep increment, using
    /// accepted resistance rather than the independently accumulated path.
    pub fn endpoint_equivalent_stress_from(
        &self,
        increment: f64,
        old: &FiniteRecoveringCreepState,
    ) -> Result<f64, FiniteJ2Error> {
        self.validate_history(old)?;
        self.prepared
            .endpoint_equivalent_stress_from(increment, old.resistance_strain)
    }

    fn return_map(
        &self,
        f: &[f64; 9],
        old: &FiniteRecoveringCreepState,
    ) -> Result<(ReturnMap, f64), FiniteJ2Error> {
        self.validate_history(old)?;
        let law = &self.prepared.native;
        let returned = if let Some(primary) = law.active_primary() {
            super::primary::standalone_with_memory(
                law,
                primary,
                f,
                &old.native,
                old.resistance_strain,
            )?
        } else {
            // Disabled primary resistance preserves the exact secondary path.
            law.backbone()
                .return_map_with_rate(f, &old.native, Some(law))?
        };
        let resistance = self.recovery.next_resistance(
            old.resistance_strain,
            returned.dgamma,
            returned.state.alpha,
        )?;
        Ok((returned, resistance))
    }
}

/// Simultaneous finite primary creep, energetic J2 plasticity and independent
/// resistance recovery, at one authored temperature and physical duration.
///
/// `native.alpha` remains permanent plastic accumulation and stores the native
/// linear/Voce/tabulated hardening energy. `native.creep_alpha` remains permanent
/// creep accumulation. Only `resistance_strain` recovers during an unloaded hold.
/// Both mechanisms retain the original joint return, potential and full tangent.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteRecoveringNortonJ2 {
    native: FiniteNortonJ2,
    prepared: FiniteNortonJ2,
    recovery: RecoveryStep,
}

impl FiniteRecoveringNortonJ2 {
    /// Select recovery of explicitly authored primary creep resistance. The
    /// rate is finite and nonnegative [1/s]; zero recovery with r=creep_alpha
    /// preserves the original combined primary law exactly.
    pub fn new(native: FiniteNortonJ2, recovery_rate_s: f64) -> Result<Self, FiniteJ2Error> {
        let (primary, recovery) = RecoveryStep::prepare(
            native.strain_hardening(),
            native.duration_s(),
            recovery_rate_s,
        )?;
        let prepared = native.clone().with_strain_hardening(primary)?;
        Ok(Self {
            native,
            prepared,
            recovery,
        })
    }

    /// Original authored elastic, energetic hardening and creep coefficients.
    /// The private duration-transformed primary scale is not a material input.
    pub const fn norton_j2(&self) -> &FiniteNortonJ2 {
        &self.native
    }

    /// Authored recovery rate at this interval's fixed temperature [1/s].
    pub const fn recovery_rate_s(&self) -> f64 {
        self.recovery.rate_s
    }

    /// Frozen physical interval for both creep and resistance recovery [s].
    pub fn duration_s(&self) -> f64 {
        self.native.duration_s()
    }

    pub fn initial_state(&self) -> FiniteRecoveringCreepState {
        FiniteRecoveringCreepState {
            native: self.native.initial_state(),
            resistance_strain: 0.0,
        }
    }

    /// Check both permanent paths and independent resistance. Plastic alpha
    /// neither bounds nor supplies primary creep resistance.
    pub fn validate_history(
        &self,
        state: &FiniteRecoveringCreepState,
    ) -> Result<(), FiniteJ2Error> {
        self.native.validate_history(&state.native)?;
        validate_resistance(state.resistance_strain, state.native.creep_alpha)
    }

    /// Inspect original elastic and energetic plastic storage at fixed history.
    pub fn stored_energy(
        &self,
        f: &[f64; 9],
        state: &FiniteRecoveringCreepState,
    ) -> Result<f64, FiniteJ2Error> {
        self.validate_history(state)?;
        self.native.stored_energy(f, &state.native)
    }

    /// Inspect an accepted state against its plastic yield surface without
    /// advancing either mechanism or recovering resistance a second time.
    pub fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FiniteRecoveringCreepState,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2Error> {
        self.validate_history(state)?;
        self.native
            .validate_committed(f, &state.native, yield_tolerance_pa)
    }

    /// Native joint return, physical mechanism heats and full analytic dP/dF.
    pub fn response(
        &self,
        f: &[f64; 9],
        old: &FiniteRecoveringCreepState,
    ) -> Result<FiniteJ2Response<FiniteRecoveringCreepState>, FiniteJ2Error> {
        let (returned, resistance_strain) = self.return_map(f, old)?;
        self.prepared
            .backbone()
            .response_from_return(f, &old.native, returned)
            .map(|response| {
                response.map_state(|native| FiniteRecoveringCreepState {
                    native,
                    resistance_strain,
                })
            })
    }

    /// Identical joint return without tangent assembly for Newton line search.
    pub fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &FiniteRecoveringCreepState,
    ) -> Result<([f64; 9], f64), FiniteJ2Error> {
        let (returned, _) = self.return_map(f, old)?;
        Ok((returned.first_piola, returned.potential))
    }

    fn return_map(
        &self,
        f: &[f64; 9],
        old: &FiniteRecoveringCreepState,
    ) -> Result<(ReturnMap, f64), FiniteJ2Error> {
        self.validate_history(old)?;
        let returned =
            self.prepared
                .return_map_with_creep_memory(f, &old.native, old.resistance_strain)?;
        // Use the original solved creep increment, not a rounded subtraction
        // of accumulated counters and never the accompanying plastic increment.
        let resistance = self.recovery.next_resistance(
            old.resistance_strain,
            returned.dcreep,
            returned.state.creep_alpha,
        )?;
        Ok((returned, resistance))
    }
}

/// The same backward-Euler resistance equation serves pure and combined creep.
/// Its transformed scale enters only the existing primary kinetic owner.
#[derive(Clone, Copy, Debug, PartialEq)]
struct RecoveryStep {
    rate_s: f64,
    divisor: f64,
}

impl RecoveryStep {
    fn prepare(
        primary: Option<&StrainHardeningCreep>,
        duration_s: f64,
        rate_s: f64,
    ) -> Result<(StrainHardeningCreep, Self), FiniteJ2Error> {
        if !rate_s.is_finite() || rate_s < 0.0 {
            return Err(FiniteJ2Error::Parameters(
                "finite primary recovery requires a finite nonnegative recovery rate",
            ));
        }
        let mut primary = *primary.ok_or(FiniteJ2Error::Parameters(
            "finite resistance recovery requires explicitly selected primary creep coefficients",
        ))?;
        let interval = rate_s * duration_s;
        let divisor = 1.0 + interval;
        if !interval.is_finite() || !divisor.is_finite() {
            return Err(FiniteJ2Error::Parameters(
                "finite resistance recovery interval is not representable",
            ));
        }
        if primary.active() {
            primary.reference_strain *= divisor;
            if !primary.reference_strain.is_finite() || primary.reference_strain <= 0.0 {
                return Err(FiniteJ2Error::Parameters(
                    "finite recovered primary strain scale is not representable",
                ));
            }
        }
        Ok((primary, Self { rate_s, divisor }))
    }

    fn next_resistance(
        &self,
        old: f64,
        increment: f64,
        accumulated: f64,
    ) -> Result<f64, FiniteJ2Error> {
        let numerator = old + increment;
        let resistance = numerator / self.divisor;
        if !numerator.is_finite() || (numerator > 0.0 && resistance <= 0.0) {
            return Err(FiniteJ2Error::Numerical(
                "finite recovered resistance history is not representable",
            ));
        }
        validate_resistance(resistance, accumulated)?;
        Ok(resistance)
    }
}

pub(in crate::plastic::finite) fn validate_resistance(
    resistance: f64,
    accumulated: f64,
) -> Result<(), FiniteJ2Error> {
    if !resistance.is_finite() || resistance < 0.0 {
        return Err(FiniteJ2Error::State(
            "finite primary resistance must be finite and nonnegative",
        ));
    }
    // Relative roundoff only: a virgin path cannot hide a positive resistance
    // behind an absolute tolerance, and accumulated+allowance need not be formed.
    let allowance = 256.0 * f64::EPSILON * accumulated.max(resistance);
    if resistance > accumulated && resistance - accumulated > allowance {
        return Err(FiniteJ2Error::State(
            "finite recovered resistance cannot exceed the accumulated permanent creep path",
        ));
    }
    Ok(())
}
