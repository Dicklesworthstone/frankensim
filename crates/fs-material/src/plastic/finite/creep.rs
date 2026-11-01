//! Threshold-free Norton creep with the native multiplicative finite history.
//!
//! The explicit law is `alpha_dot=(s0/eta)*(q/s0)^n`, where `q` is equivalent
//! KIRCHHOFF stress, `eta` is in Pa s, `s0` is in Pa and `n>=1`. The linear
//! case is `alpha_dot=q/eta`; there is no yield radius or hardening storage.
//! The same implicit finite return, metric reconstruction, domain checks and
//! full analytic tangent as [`super::rate`] solve
//! `q + 3mu*dgamma = q_trial` and `dgamma=dt*alpha_dot(q)` simultaneously.
//!
//! `F=Fe*Fc`, with retained `Cc^-1` and accumulated equivalent creep. For
//! compatibility these use [`FiniteJ2State`]: its `inverse_plastic_cauchy_green`
//! and `alpha` fields mean the creep metric and accumulated creep here. The
//! `plastic_*` response and limit fields likewise refer to this sole creep
//! mechanism. This is not a simultaneous creep-plus-plasticity model; source
//! selection must retain this distinction rather than relabel old history.
//!
//! Physical heat is `q*dgamma`. The incremental potential adds
//! `n/(n+1)*q*dgamma` to elastic storage; its gradient is Piola stress. There
//! is no artificial energetic hardening, equilibrium yield projection or
//! inference of creep parameters from melting temperature. Every positive
//! resolved stress drives flow; an increment below retained-history precision
//! refuses explicitly. The exact linear relaxed tangent holds at zero stress.
//! Finite strain and accumulated-history bounds remain explicit. Primary or
//! tertiary creep, rupture, damage and empirical calibration are outside this
//! bounded secondary-creep law.

use super::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use super::{
    FiniteJ2, FiniteJ2Error, FiniteJ2Limits, FiniteJ2Parameters, FiniteJ2Response, FiniteJ2State,
};

/// One explicitly authored finite secondary-creep law and physical interval.
/// It deliberately has no rate-independent `FiniteJ2Material` adapter.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteNortonCreep {
    native: FiniteOverstressJ2,
}

impl FiniteNortonCreep {
    /// Admit independent elastic moduli, finite-strain bounds, physical Norton
    /// kinetics and a strictly positive finite physical duration.
    pub fn new(
        shear_modulus_pa: f64,
        bulk_modulus_pa: f64,
        limits: FiniteJ2Limits,
        kinetics: FiniteOverstressParameters,
        duration_s: f64,
    ) -> Result<Self, FiniteJ2Error> {
        let backbone = FiniteJ2::new_creep_backbone(
            FiniteJ2Parameters {
                shear_modulus_pa,
                bulk_modulus_pa,
                yield_stress_pa: 0.0,
                hardening_modulus_pa: 0.0,
            },
            limits,
        )?;
        Ok(Self {
            native: FiniteOverstressJ2::new(backbone, kinetics, duration_s)?,
        })
    }

    /// Original-reference-volume Hencky shear modulus.
    pub fn shear_modulus_pa(&self) -> f64 {
        self.native.backbone().parameters().shear_modulus_pa
    }
    /// Original-reference-volume Hencky bulk modulus.
    pub fn bulk_modulus_pa(&self) -> f64 {
        self.native.backbone().parameters().bulk_modulus_pa
    }
    /// Native applicability bounds; the plastic-named fields limit creep here.
    pub fn limits(&self) -> &FiniteJ2Limits {
        self.native.backbone().limits()
    }
    /// Explicit Norton viscosity, stress scale and exponent at this temperature.
    pub fn kinetics(&self) -> &FiniteOverstressParameters {
        self.native.parameters()
    }
    /// Actual physical interval used by the implicit return.
    pub fn duration_s(&self) -> f64 {
        self.native.duration_s()
    }
    /// Virgin creep metric and zero accumulated equivalent creep.
    pub fn initial_state(&self) -> FiniteJ2State {
        self.native.initial_state()
    }
    /// Validate the entire isochoric metric, accumulated path and finite bounds.
    pub fn validate_history(&self, state: &FiniteJ2State) -> Result<(), FiniteJ2Error> {
        self.native.validate_history(state)
    }
    /// Inspect elastic storage at fixed creep history; this performs no return.
    pub fn stored_energy(&self, f: &[f64; 9], state: &FiniteJ2State) -> Result<f64, FiniteJ2Error> {
        self.native.stored_energy(f, state)
    }
    /// Inspect finite accepted creep stress and history without projecting it
    /// to zero stress or advancing it through another physical interval.
    pub fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
    ) -> Result<(), FiniteJ2Error> {
        self.native.validate_committed(f, state, 0.0)
    }
    /// Native finite implicit return and exact derivative of the resulting P.
    pub fn response(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
    ) -> Result<FiniteJ2Response, FiniteJ2Error> {
        self.native.response(f, old)
    }
    /// Identical stress and reduced potential without allocating a full tangent.
    pub fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
    ) -> Result<([f64; 9], f64), FiniteJ2Error> {
        self.native.stress_and_potential(f, old)
    }
    /// Endpoint equivalent Kirchhoff stress implied by a nonnegative increment.
    pub fn endpoint_equivalent_stress_pa(&self, increment: f64) -> Result<f64, FiniteJ2Error> {
        self.native.endpoint_overstress_pa(increment)
    }
}
