//! Independent primary resistance recovery on the native finite Norton owner.
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
//!
//! This is the same authored phenomenological resistance model as
//! [`crate::creep::recovery`], now with finite kinematics. Its coefficients
//! require independent calibration. This owner has no automatic temperature
//! shift, source-card binding, simultaneous yielding, damage or rupture law.
//! Consumers must retain the complete typed history, including resistance.

use super::FiniteNortonCreep;
use crate::plastic::finite::{FiniteJ2Error, FiniteJ2Response, FiniteJ2State, ReturnMap};

/// Permanent finite creep metric/path and separately recovering resistance.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteRecoveringCreepState {
    /// Original native history: alpha is cumulative creep, never resistance.
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
    recovery_rate_s: f64,
    recovery_divisor: f64,
}

impl FiniteRecoveringNortonCreep {
    /// Admit an explicitly selected primary law and nonnegative finite recovery
    /// rate [1/s] at the same fixed temperature as the native Norton kinetics.
    /// Zero recovery preserves the original primary law exactly when r=alpha.
    pub fn new(native: FiniteNortonCreep, recovery_rate_s: f64) -> Result<Self, FiniteJ2Error> {
        if !recovery_rate_s.is_finite() || recovery_rate_s < 0.0 {
            return Err(FiniteJ2Error::Parameters(
                "finite primary recovery requires a finite nonnegative recovery rate",
            ));
        }
        let mut primary = *native.strain_hardening().ok_or(FiniteJ2Error::Parameters(
            "finite resistance recovery requires explicitly selected primary creep coefficients",
        ))?;
        let interval = recovery_rate_s * native.duration_s();
        let recovery_divisor = 1.0 + interval;
        if !interval.is_finite() || !recovery_divisor.is_finite() {
            return Err(FiniteJ2Error::Parameters(
                "finite resistance recovery interval is not representable",
            ));
        }
        if primary.active() {
            primary.reference_strain *= recovery_divisor;
            if !primary.reference_strain.is_finite() || primary.reference_strain <= 0.0 {
                return Err(FiniteJ2Error::Parameters(
                    "finite recovered primary strain scale is not representable",
                ));
            }
        }
        let prepared = native.clone().with_strain_hardening(primary)?;
        Ok(Self {
            native,
            prepared,
            recovery_rate_s,
            recovery_divisor,
        })
    }

    /// Original authored Norton coefficients, limits, primary law and duration.
    /// The private duration-transformed primary scale is not a material input.
    pub const fn norton(&self) -> &FiniteNortonCreep {
        &self.native
    }

    /// Authored recovery rate at this interval's fixed temperature [1/s].
    pub const fn recovery_rate_s(&self) -> f64 {
        self.recovery_rate_s
    }

    /// Virgin metric, zero permanent path and zero primary resistance.
    pub fn initial_state(&self) -> FiniteRecoveringCreepState {
        FiniteRecoveringCreepState {
            native: self.native.initial_state(),
            resistance_strain: 0.0,
        }
    }

    /// Check all original finite-history bounds and independent resistance.
    pub fn validate_history(&self, state: &FiniteRecoveringCreepState) -> Result<(), FiniteJ2Error> {
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
            law.backbone().return_map_with_rate(f, &old.native, Some(law))?
        };
        let numerator = old.resistance_strain + returned.dgamma;
        let resistance = numerator / self.recovery_divisor;
        if !numerator.is_finite() || (numerator > 0.0 && resistance <= 0.0) {
            return Err(FiniteJ2Error::Numerical(
                "finite recovered resistance history is not representable",
            ));
        }
        validate_resistance(resistance, returned.state.alpha)?;
        Ok((returned, resistance))
    }
}

fn validate_resistance(resistance: f64, alpha: f64) -> Result<(), FiniteJ2Error> {
    if !resistance.is_finite() || resistance < 0.0 {
        return Err(FiniteJ2Error::State(
            "finite primary resistance must be finite and nonnegative",
        ));
    }
    // Relative roundoff only: a virgin path cannot hide a positive resistance
    // behind an absolute tolerance, and alpha+allowance need not be formed.
    let allowance = 256.0 * f64::EPSILON * alpha.max(resistance);
    if resistance > alpha && resistance - alpha > allowance {
        return Err(FiniteJ2Error::State(
            "finite recovered resistance cannot exceed the accumulated permanent creep path",
        ));
    }
    Ok(())
}
