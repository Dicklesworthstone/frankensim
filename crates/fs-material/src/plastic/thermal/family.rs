//! Complete native histories on a shared thermoplastic feedback path.
//!
//! A family resolves one physical interval and audits its original constitutive
//! update. Spatial consumers retain the associated native state: J2 history for
//! the original source, and exactly N tensor backstresses for thermal Chaboche.
//! No erased or padded per-point history, replacement return map, or additional
//! thermal integrator is introduced. All trial intervals reuse the same accepted
//! history, old enthalpy and duration.

use super::chaboche::{ThermalChaboche, ThermalChabocheError, ThermalChabocheSource};
use super::{FeedbackError, ThermalJ2, ThermoplasticIncrement, ThermoplasticSource};
use crate::plastic::J2State;
use crate::plastic::chaboche::{ChabocheError, ChabocheState};
use crate::plastic::energy::{EnergyError, EnergyTolerances};
use crate::plastic::source::ResolvedJ2Plasticity;
use crate::{SmallStrainLaw, Voigt};
use core::fmt;

/// One immutable endpoint law, retaining the complete original source bundle.
#[derive(Clone, Debug)]
pub struct ThermoplasticStep<L> {
    pub temperature_k: f64,
    pub material: ResolvedJ2Plasticity,
    pub law: L,
}

/// Complete native candidate and its shared caloric/mechanical accounts.
#[derive(Clone, Debug)]
pub struct ThermoplasticFamilyIncrement<S> {
    /// The full native history; never reconstruct it from the J2 projection.
    pub state: S,
    /// Includes every branch's storage, entropy and physical recovery heat.
    /// `mechanical.state` retains the original J2 projection; `state` above
    /// supplies the complete state that must agree with native equilibrium.
    pub thermal: ThermoplasticIncrement,
}

/// Original source, history, return-map or physical-account refusal.
#[derive(Debug)]
pub enum ThermoplasticFamilyError {
    Thermal(FeedbackError),
    Energy(EnergyError),
    Chaboche(ChabocheError),
}
impl fmt::Display for ThermoplasticFamilyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Thermal(error) => error.fmt(f),
            Self::Energy(error) => error.fmt(f),
            Self::Chaboche(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for ThermoplasticFamilyError {}
impl From<FeedbackError> for ThermoplasticFamilyError {
    fn from(error: FeedbackError) -> Self {
        Self::Thermal(error)
    }
}
impl From<EnergyError> for ThermoplasticFamilyError {
    fn from(error: EnergyError) -> Self {
        Self::Energy(error)
    }
}
impl From<ChabocheError> for ThermoplasticFamilyError {
    fn from(error: ChabocheError) -> Self {
        Self::Chaboche(error)
    }
}
impl From<ThermalChabocheError> for ThermoplasticFamilyError {
    fn from(error: ThermalChabocheError) -> Self {
        match error {
            ThermalChabocheError::Thermal(error) => Self::Thermal(error),
            ThermalChabocheError::Mechanical(error) => Self::Chaboche(error),
        }
    }
}

/// Source-pinned native histories admitted by the shared thermal owner.
///
/// `resolve_step` and `checked_increment` represent the same physical interval.
/// Their old enthalpy/history and duration stay frozen during all equilibrium
/// and thermal trials. `validate_accepted` only inspects committed history; it
/// must neither return-map that history again nor advance constitutive time.
pub trait ThermoplasticSourceFamily<'a>: Clone + fmt::Debug {
    type Material: SmallStrainLaw + Clone + fmt::Debug;

    /// Original caloric, reference-mass and material/source identity owner.
    fn thermal_source(&self) -> &ThermoplasticSource<'a>;

    fn resolve_step(
        &self,
        old_h: f64,
        new_h: f64,
        duration_s: f64,
    ) -> Result<ThermoplasticStep<Self::Material>, ThermoplasticFamilyError>;

    #[allow(clippy::too_many_arguments)]
    fn checked_increment(
        &self,
        old_h: f64,
        new_h: f64,
        old_strain: &Voigt,
        old: &<Self::Material as SmallStrainLaw>::State,
        new_strain: &Voigt,
        tolerances: EnergyTolerances,
        duration_s: f64,
    ) -> Result<
        ThermoplasticFamilyIncrement<<Self::Material as SmallStrainLaw>::State>,
        ThermoplasticFamilyError,
    >;

    fn validate_accepted(
        &self,
        h: f64,
        total_strain: &Voigt,
        state: &<Self::Material as SmallStrainLaw>::State,
        yield_tolerance_pa: f64,
    ) -> Result<(), ThermoplasticFamilyError>;
}

impl<'a> ThermoplasticSourceFamily<'a> for ThermoplasticSource<'a> {
    type Material = ThermalJ2;

    fn thermal_source(&self) -> &ThermoplasticSource<'a> {
        self
    }

    fn resolve_step(
        &self,
        _old_h: f64,
        new_h: f64,
        duration_s: f64,
    ) -> Result<ThermoplasticStep<Self::Material>, ThermoplasticFamilyError> {
        let (phase, law) = self.resolve_expanded_enthalpy_for_step(new_h, duration_s)?;
        Ok(ThermoplasticStep {
            temperature_k: phase.temperature_k(),
            material: law.resolved().clone(),
            law,
        })
    }

    fn checked_increment(
        &self,
        old_h: f64,
        new_h: f64,
        old_strain: &Voigt,
        old: &J2State,
        new_strain: &Voigt,
        tolerances: EnergyTolerances,
        duration_s: f64,
    ) -> Result<ThermoplasticFamilyIncrement<J2State>, ThermoplasticFamilyError> {
        let thermal = self.advance_with_energy_for_step(
            old_h, new_h, old_strain, old, new_strain, tolerances, duration_s,
        )?;
        Ok(ThermoplasticFamilyIncrement {
            state: thermal.mechanical.state.clone(),
            thermal,
        })
    }

    fn validate_accepted(
        &self,
        h: f64,
        total_strain: &Voigt,
        state: &J2State,
        yield_tolerance_pa: f64,
    ) -> Result<(), ThermoplasticFamilyError> {
        valid_yield_tolerance(yield_tolerance_pa)?;
        let (phase, material) = self.resolve_enthalpy(h)?;
        let mechanical = self.mechanical_strain(phase.temperature_k(), total_strain)?;
        self.validate_committed_history(&material, &mechanical, state, yield_tolerance_pa)?;
        Ok(())
    }
}

impl<'a, const N: usize> ThermoplasticSourceFamily<'a> for ThermalChabocheSource<'a, N> {
    type Material = ThermalChaboche<N>;

    fn thermal_source(&self) -> &ThermoplasticSource<'a> {
        self.base()
    }

    fn resolve_step(
        &self,
        old_h: f64,
        new_h: f64,
        duration_s: f64,
    ) -> Result<ThermoplasticStep<Self::Material>, ThermoplasticFamilyError> {
        let (phase, law) = self.resolve_enthalpy_for_step(old_h, new_h, duration_s)?;
        Ok(ThermoplasticStep {
            temperature_k: phase.temperature_k(),
            material: law.resolved().clone(),
            law,
        })
    }

    fn checked_increment(
        &self,
        old_h: f64,
        new_h: f64,
        old_strain: &Voigt,
        old: &ChabocheState<N>,
        new_strain: &Voigt,
        tolerances: EnergyTolerances,
        duration_s: f64,
    ) -> Result<ThermoplasticFamilyIncrement<ChabocheState<N>>, ThermoplasticFamilyError> {
        let increment = self.advance_with_energy_for_step(
            old_h, new_h, old_strain, old, new_strain, tolerances, duration_s,
        )?;
        Ok(ThermoplasticFamilyIncrement {
            state: increment.state,
            thermal: increment.thermal,
        })
    }

    fn validate_accepted(
        &self,
        h: f64,
        total_strain: &Voigt,
        state: &ChabocheState<N>,
        yield_tolerance_pa: f64,
    ) -> Result<(), ThermoplasticFamilyError> {
        valid_yield_tolerance(yield_tolerance_pa)?;
        let (phase, material) = self.base().resolve_enthalpy(h)?;
        let mechanical = self
            .base()
            .mechanical_strain(phase.temperature_k(), total_strain)?;
        let law = self.law_at(&material)?;
        law.validate_history(&mechanical, state)?;

        // Recover already accepted stress elastically. Another constitutive
        // response at the same strain would relax a viscous state for a second
        // interval, and could conceal an invalid rate-independent endpoint.
        let elastic = core::array::from_fn(|i| mechanical[i] - state.plastic.plastic_strain[i]);
        let stress = law.backbone().backbone().elastic.stress(&elastic, &());
        let residual = law.yield_function_with_state(&stress, state);
        if !residual.is_finite() {
            return Err(
                EnergyError::Invalid("accepted cyclic yield residual is unrepresentable").into(),
            );
        }
        if !self.base().has_overstress() && residual > yield_tolerance_pa {
            return Err(EnergyError::YieldResidual {
                residual_pa: residual,
                tolerance_pa: yield_tolerance_pa,
            }
            .into());
        }
        Ok(())
    }
}

fn valid_yield_tolerance(value: f64) -> Result<(), EnergyError> {
    if !value.is_finite() || value <= 0.0 {
        Err(EnergyError::Invalid(
            "yield allowance must be finite and positive",
        ))
    } else {
        Ok(())
    }
}
