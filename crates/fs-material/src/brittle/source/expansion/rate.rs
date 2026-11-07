//! Physical-time fracture-plane viscosity with source-bound thermal eigenstrain.
//!
//! The original source first transfers coefficients at held mechanical strain.
//! A straight TOTAL-strain path and an explicitly affine FREE-strain path then
//! give a straight mechanical path under those fixed new coefficients:
//! `delta_epsilon_m = delta_epsilon_total - delta_free * I`. The native timed
//! brittle/cohesive owner supplies the full response, history, derivatives and
//! separately integrated viscous loss; no second constitutive state is stored.
//!
//! Its total path-mean stress is conjugate to both ports. Displacement work is
//! `mean_stress : delta_epsilon_total`; eigenstrain work is
//! `-trace(mean_stress) * delta_free`. Their sum is native mechanical-strain
//! work, including intrinsic fracture and viscous loss. The motion owner retains
//! the OLD source to evaluate this split and adds coefficient work separately.
//! Neither port nor either dissipation is automatically converted into heat.
//!
//! Interval endpoint stress uses the mechanical secant rate. To observe an
//! actual endpoint velocity, subtract the declared free-strain rate from the
//! total strain rate before evaluating the surviving-bond dashpots. Affine
//! free strain supplies `(free_next-free_old)/dt`, not alpha at either endpoint
//! times a guessed temperature rate. A separately supplied actual temperature
//! trajectory can use the source helper in [`source`].

/// Timed coefficient and eigenstrain evolution at explicitly held TOTAL strain.
pub mod restrained;
/// Original fracture, viscosity, expansion and optional reference-mass sources.
pub mod source;
pub use source::{
    resolve_thermal_viscous_brittle, resolve_thermal_viscous_brittle_dynamic_material,
    thermal_viscous_brittle_dynamic_requirements, thermal_viscous_brittle_requirements,
};

use super::{BrittleSourceError, ResolvedThermalBrittle};
use crate::brittle::rate::{BrittleRatePathStep, BrittleRateStep};
use crate::brittle::{BrittleState, finite, finite_voigt, positive};
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::ModelCard;

/// TOTAL-strain endpoint-force adapter at one fixed positive physical duration.
/// Its sole history is the original native MECHANICAL-strain state.
#[derive(Clone, Copy, Debug)]
pub struct ThermalBrittleRateSmallStrain<'a> {
    owner: &'a ResolvedThermalBrittle,
    dt_s: f64,
}

/// TOTAL-strain work-conjugate mean-force adapter at fixed new coefficients.
/// A missing native classical path tangent remains an explicit refusal.
#[derive(Clone, Copy, Debug)]
pub struct ThermalBrittleRateMeanSmallStrain<'a> {
    owner: &'a ResolvedThermalBrittle,
    dt_s: f64,
}

impl ResolvedThermalBrittle {
    /// Timed native endpoint at a TOTAL tensor-strain endpoint.
    ///
    /// First transfer original coefficients/history into this source with
    /// `transfer_coefficients`. Its retained OLD mechanical strain is then the
    /// start of the physical interval; subtracting the NEW free strain here
    /// includes the whole eigenstrain increment in the native secant rate.
    /// Returned work is mechanical-strain work, not displacement work alone.
    pub fn trial_total_over_time<'a>(
        &'a self,
        total: &Voigt,
        state: &BrittleState<'a>,
        dt_s: f64,
    ) -> Result<BrittleRateStep<'a>, BrittleSourceError> {
        let mechanical = self.mechanical_strain(total)?;
        Ok(self
            .band()
            .law()
            .trial_over_time(&mechanical, state, dt_s)?)
    }

    /// Timed endpoint plus the exact TOTAL mean stress on the same mechanical
    /// path. Both physical work ports use this stress, including its viscosity.
    /// This force-only query avoids evaluating an unused path derivative.
    pub fn trial_total_over_time_with_mean_stress<'a>(
        &'a self,
        total: &Voigt,
        state: &BrittleState<'a>,
        dt_s: f64,
    ) -> Result<BrittleRatePathStep<'a>, BrittleSourceError> {
        let mechanical = self.mechanical_strain(total)?;
        Ok(self
            .band()
            .law()
            .trial_over_time_with_mean_stress(&mechanical, state, dt_s)?)
    }

    /// Timed path and d(mean stress)/d(NEW TOTAL strain), at fixed original
    /// history, new source temperature and dt. Additive eigenstrain has no
    /// displacement derivative, so this is exactly the native mechanical
    /// derivative. It is not a derivative with respect to temperature.
    pub fn trial_total_over_time_with_mean_stress_tangent<'a>(
        &'a self,
        total: &Voigt,
        state: &BrittleState<'a>,
        dt_s: f64,
    ) -> Result<(BrittleRatePathStep<'a>, Option<Tangent6>), BrittleSourceError> {
        let mechanical = self.mechanical_strain(total)?;
        Ok(self
            .band()
            .law()
            .trial_over_time_with_mean_stress_tangent(&mechanical, state, dt_s)?)
    }

    /// Borrow the timed endpoint force while retaining the native history.
    pub fn small_strain_law_over_time(
        &self,
        dt_s: f64,
    ) -> Result<ThermalBrittleRateSmallStrain<'_>, BrittleSourceError> {
        positive(dt_s)?;
        Ok(ThermalBrittleRateSmallStrain { owner: self, dt_s })
    }

    /// Borrow the timed work-conjugate mean force at one fixed physical dt.
    pub fn mean_small_strain_law_over_time(
        &self,
        dt_s: f64,
    ) -> Result<ThermalBrittleRateMeanSmallStrain<'_>, BrittleSourceError> {
        positive(dt_s)?;
        Ok(ThermalBrittleRateMeanSmallStrain { owner: self, dt_s })
    }

    /// Convert an actual TOTAL tensor strain rate to mechanical rate [1/s].
    /// `free_rate_per_s` is the explicitly supplied derivative of the original
    /// isotropic free linear strain, not a temperature rate. Shear is unchanged.
    /// No duration or thermal trajectory is inferred from a source snapshot.
    #[allow(clippy::unused_self)] // Keep the rate map with the source-bound strain map.
    pub fn mechanical_strain_rate(
        &self,
        total_rate: &Voigt,
        free_rate_per_s: f64,
    ) -> Result<Voigt, BrittleSourceError> {
        finite_voigt(total_rate)?;
        finite(free_rate_per_s)?;
        let mechanical =
            core::array::from_fn(|i| total_rate[i] - if i < 3 { free_rate_per_s } else { 0.0 });
        finite_voigt(&mechanical)?;
        Ok(mechanical)
    }

    /// Observe total stress at an accepted endpoint and its actual velocity.
    ///
    /// The caller supplies the free-strain rate of its declared thermal path;
    /// an affine free-strain increment uses `delta_free/dt` even when the
    /// instantaneous displacement velocity differs from its secant. Original
    /// owner, strain and mechanical jump-rate admission apply. This read-only
    /// observation does not advance damage or add energy to any ledger.
    /// Without surviving dashpots, finite supplied rates remain mandatory but
    /// their mechanically irrelevant subtraction/projection is not evaluated.
    pub fn stress_at_rate(
        &self,
        state: &BrittleState<'_>,
        total_rate: &Voigt,
        free_rate_per_s: f64,
    ) -> Result<Voigt, BrittleSourceError> {
        self.total_strain(state)?;
        let mechanical_rate = if self.band().law().viscosity().is_some() && !state.failed() {
            self.mechanical_strain_rate(total_rate, free_rate_per_s)?
        } else {
            finite_voigt(total_rate)?;
            finite(free_rate_per_s)?;
            [0.0; 6]
        };
        Ok(self.band().law().stress_at_rate(state, &mechanical_rate)?)
    }
}

impl<'a> ThermalBrittleRateSmallStrain<'a> {
    /// Explicit physical interval shared by all constitutive trial iterates.
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.dt_s
    }

    /// Checked native endpoint; the returned history has not been published.
    pub fn trial(
        &self,
        total: &Voigt,
        state: &BrittleState<'a>,
    ) -> Result<BrittleRateStep<'a>, BrittleSourceError> {
        self.owner.trial_total_over_time(total, state, self.dt_s)
    }
}

impl<'a> ThermalBrittleRateMeanSmallStrain<'a> {
    /// Explicit physical interval shared by all constitutive trial iterates.
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.dt_s
    }

    /// Checked native path and optional classical TOTAL-strain derivative.
    pub fn trial(
        &self,
        total: &Voigt,
        state: &BrittleState<'a>,
    ) -> Result<(BrittleRatePathStep<'a>, Option<Tangent6>), BrittleSourceError> {
        self.owner
            .trial_total_over_time_with_mean_stress_tangent(total, state, self.dt_s)
    }
}

impl<'a> SmallStrainLaw for ThermalBrittleRateSmallStrain<'a> {
    type State = BrittleState<'a>;

    fn initial_state(&self) -> Self::State {
        self.owner.band().law().initial_state()
    }
    fn stress(&self, total: &Voigt, state: &Self::State) -> Voigt {
        self.trial(total, state)
            .map(|step| step.stress_pa)
            .unwrap_or([f64::NAN; 6])
    }
    fn tangent(&self, total: &Voigt, state: &Self::State) -> Tangent6 {
        self.trial(total, state)
            .map(|step| step.branch_tangent_pa)
            .unwrap_or([[f64::NAN; 6]; 6])
    }
    fn update_state(&self, total: &Voigt, state: &Self::State) -> Self::State {
        self.trial(total, state)
            .map(|step| step.equilibrium.state)
            .unwrap_or_else(|_| invalid_state(self.owner))
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.owner.small_strain_law().admissibility()
    }
    fn card(&self) -> ModelCard {
        rate_card(self.owner, false)
    }
}

impl<'a> SmallStrainLaw for ThermalBrittleRateMeanSmallStrain<'a> {
    type State = BrittleState<'a>;

    fn initial_state(&self) -> Self::State {
        self.owner.band().law().initial_state()
    }
    fn stress(&self, total: &Voigt, state: &Self::State) -> Voigt {
        self.owner
            .trial_total_over_time_with_mean_stress(total, state, self.dt_s)
            .map(|path| path.mean_stress_pa)
            .unwrap_or([f64::NAN; 6])
    }
    fn tangent(&self, total: &Voigt, state: &Self::State) -> Tangent6 {
        self.trial(total, state)
            .ok()
            .and_then(|(_, tangent)| tangent)
            .unwrap_or([[f64::NAN; 6]; 6])
    }
    fn update_state(&self, total: &Voigt, state: &Self::State) -> Self::State {
        self.owner
            .trial_total_over_time(total, state, self.dt_s)
            .map(|step| step.equilibrium.state)
            .unwrap_or_else(|_| invalid_state(self.owner))
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.owner.small_strain_law().admissibility()
    }
    fn card(&self) -> ModelCard {
        rate_card(self.owner, true)
    }
}

fn invalid_state(owner: &ResolvedThermalBrittle) -> BrittleState<'_> {
    BrittleState {
        owner: owner.band().law(),
        strain: [f64::NAN; 6],
        point: owner.band().law().cohesive.bonded_point(),
    }
}

fn rate_card(owner: &ResolvedThermalBrittle, mean: bool) -> ModelCard {
    let mut card = owner.small_strain_law().card();
    card.name = if mean {
        "material.brittle-source-thermal-rate-path-mean"
    } else {
        "material.brittle-source-thermal-rate-endpoint"
    }
    .into();
    card.assumptions.push("fixed new coefficients after original coefficient-only history transfer; total strain and free isotropic strain follow one affine path over explicit positive physical dt".into());
    card.assumptions.push("original fracture-plane viscosity acts on mechanical jump rate in Pa s/m; intrinsic fracture, physical viscous loss, coefficient work and eigenstrain work remain separate".into());
    card.known_failures.push("interval secant-rate stress differs from actual endpoint-velocity stress; the latter requires the explicit free-strain rate of the declared thermal path".into());
    card.known_failures.push("no arbitrary-width objectivity of additional viscous work, rate-dependent intrinsic strength, thermal Jacobian or unresolved thermal-trajectory accuracy is inferred".into());
    if mean {
        card.known_failures.push("native classical mean tangent can be absent at nonsmooth contact/damage paths; only the consumer may admit a fallback for this same mean force".into());
    }
    card.assumptions.sort_unstable();
    card.known_failures.sort_unstable();
    card
}
