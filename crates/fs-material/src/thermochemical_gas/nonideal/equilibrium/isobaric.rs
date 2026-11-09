//! Pure-fluid sensible/latent storage on the original finite thermal timeline.
//!
//! The original PR coexistence pressure is held fixed. Conserved enthalpy
//! selects temperature and phase amount; latent heat and pressure-volume work
//! already belong to that coordinate. Geometry is a reference mesh, not a
//! rigid-volume constraint. No heat integrator, gas/liquid transport average,
//! nucleation rate or multicomponent phase allocation is inferred.

/// Explicit bulk conductivity on this same sensible/latent caloric path.
pub mod conduction;

use super::{PureFluidEquilibriumError, PureFluidEquilibriumPhase, PureFluidEquilibriumState,
    ThermochemicalPengRobinsonModel};
use crate::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_math::dd::Dd;
use fs_qty::{Pressure, Temperature};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw, ThermalStorageBasis};
use fs_thermochem::eos::peng_robinson::saturation::{
    PengRobinsonSaturationErrorV1, PengRobinsonSaturationOptionsV1, PengRobinsonSaturationV1,
    isobaric::{PengRobinsonIsobaricCaloricV1, PureFluidIsobaricErrorV1, PureFluidIsobaricPhaseV1,
        PureFluidIsobaricSolutionV1, PureFluidIsobaricSolveV1, PureFluidIsobaricStateV1},
};

/// Original source or phase refusal leaves every caller-owned state unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum PureFluidIsobaricStorageError {
    Caloric(PureFluidIsobaricErrorV1),
    Material(PureFluidEquilibriumError),
    Capacity(GasEnergyError),
    Phase(PengRobinsonSaturationErrorV1),
    Invalid(&'static str),
    Cancelled,
}
impl core::fmt::Display for PureFluidIsobaricStorageError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pure-fluid isobaric thermal storage: {self:?}")
    }
}
impl std::error::Error for PureFluidIsobaricStorageError {}
impl From<PureFluidIsobaricErrorV1> for PureFluidIsobaricStorageError {
    fn from(error: PureFluidIsobaricErrorV1) -> Self {
        if matches!(error, PureFluidIsobaricErrorV1::Cancelled) { Self::Cancelled }
        else { Self::Caloric(error) }
    }
}
type E = PureFluidIsobaricStorageError;

/// The original complete phase material state and the independently retained
/// source-relative enthalpy residual/phase derivatives. No effective transport
/// or automatic frozen/equilibrium acoustic policy is attached.
#[derive(Clone, Debug, PartialEq)]
pub struct PureFluidIsobaricThermalState {
    material: PureFluidEquilibriumState,
    caloric: PureFluidIsobaricSolutionV1,
}
impl PureFluidIsobaricThermalState {
    pub const fn material(&self) -> &PureFluidEquilibriumState { &self.material }
    pub const fn caloric(&self) -> &PureFluidIsobaricSolutionV1 { &self.caloric }
}

/// Immutable original-source enthalpy law with a genuine latent interval.
///
/// h=0 is saturated liquid, so compressed-liquid sensible h is negative.
/// Thermal entropy is the defining integral dh/T. The separate original-source
/// entropy query retains the small accepted coexistence H-T*Delta S defect.
/// These coordinates carry no implicit open-flow/advection datum binding.
pub struct PureFluidIsobaricThermalStorage<'a> {
    model: &'a ThermochemicalPengRobinsonModel,
    caloric: PengRobinsonIsobaricCaloricV1<'a>,
    maximum_slope: f64,
}

impl ThermochemicalPengRobinsonModel {
    /// Prepare original coexistence at a declared transition temperature, then
    /// keep its accepted pressure fixed over the sensible/latent storage path.
    /// The native saturation owner retains all work/source/refusal semantics;
    /// the returned pressure is not replaced by a fitted boiling-point curve.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_isobaric_thermal_storage(
        &self,
        transition_temperature: Temperature,
        pressure_bracket: [Pressure; 2],
        temperature_domain: [Temperature; 2],
        controls: PureFluidIsobaricSolveV1,
        saturation_options: PengRobinsonSaturationOptionsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PureFluidIsobaricThermalStorage<'_>, E> {
        poll(&mut cancelled)?;
        controls.validate()?;
        let saturation = self.mechanics().pure_saturation_at_temperature_cancellable(
            transition_temperature, pressure_bracket, saturation_options, &mut cancelled)
            .map_err(|error| if matches!(error, PengRobinsonSaturationErrorV1::Cancelled) {
                E::Cancelled
            } else { E::Phase(error) })?;
        PureFluidIsobaricThermalStorage::new(self, saturation, temperature_domain, controls, cancelled)
    }
}

impl<'a> PureFluidIsobaricThermalStorage<'a> {
    /// Bind an original accepted saturation state and prove a conservative
    /// global inverse slope from original source coefficients, not Cp samples.
    /// Both endpoint material states must admit the same unchanged source bank.
    pub fn new(
        model: &'a ThermochemicalPengRobinsonModel,
        saturation: PengRobinsonSaturationV1,
        temperature_domain: [Temperature; 2],
        controls: PureFluidIsobaricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, E> {
        poll(&mut cancelled)?;
        if model.source().thermodynamics().components().len() != 1 {
            return Err(E::Invalid("pure-fluid thermal storage requires exactly one original species"));
        }
        let source = &model.source().thermodynamics().components()[0];
        let caloric = PengRobinsonIsobaricCaloricV1::new(model.mechanics(), source,
            saturation, temperature_domain, controls, &mut cancelled)?;
        // The native owner proves a''(T)>=0 on this original pure subcritical
        // PR alpha branch. Its stable sensible roots therefore have
        // Cp = Cv - T*p_T^2/p_v >= Cv >= Cv_ideal. On coexistence dT/dh=0.
        // Reuse the original outward coefficient bound of the complete NASA
        // interval; no frozen mixture Cp sample is a global bound.
        let minimum_cv = model.source().prepare_caloric_branch(GasEnergySolve {
            minimum_temperature_k: temperature_domain[0].value(),
            maximum_temperature_k: temperature_domain[1].value(),
            energy_tolerance_j_kg: controls.enthalpy_tolerance_j_kg,
            maximum_iterations: controls.maximum_evaluations.min(128),
        }, &mut cancelled).map_err(|error| if matches!(error, GasEnergyError::Cancelled) {
            E::Cancelled
        } else { E::Capacity(error) })?.minimum_cv_j_kg_k();
        let maximum_slope = (1.0 / minimum_cv).next_up();
        if !positive(minimum_cv) || !positive(maximum_slope) {
            return Err(E::Invalid("no represented global source Cv bound for the pure-fluid storage path"));
        }
        let law = Self { model, caloric, maximum_slope };
        for h in law.specific_enthalpy_domain_j_kg() {
            law.resolve_specific_enthalpy(h, &mut cancelled)?;
        }
        // Admit both exact latent endpoints through the same original material
        // publication owner, without selecting or averaging gas transport data.
        for h in [0.0, law.caloric.latent_enthalpy_j_kg()] {
            law.resolve_specific_enthalpy(h, &mut cancelled)?;
        }
        poll(&mut cancelled)?;
        Ok(law)
    }
    pub const fn model(&self) -> &'a ThermochemicalPengRobinsonModel { self.model }
    pub const fn original_calorics(&self) -> &PengRobinsonIsobaricCaloricV1<'a> { &self.caloric }
    pub fn pressure(&self) -> Pressure { self.caloric.pressure() }
    pub fn temperature_domain(&self) -> [Temperature; 2] { self.caloric.temperature_domain() }
    pub fn specific_enthalpy_domain_j_kg(&self) -> [f64; 2] { self.caloric.specific_enthalpy_domain_j_kg() }
    pub fn coexistence_specific_enthalpy_domain_j_kg(&self) -> [f64; 2] {
        [0.0, self.caloric.latent_enthalpy_j_kg()]
    }
    /// Original measured latent H-T*Delta S defect, not fitted out of the source.
    pub fn latent_closure_residual_j_kg(&self) -> f64 { self.caloric.latent_closure_residual_j_kg() }

    /// The original source temperature, phase amount, density and derivatives.
    pub fn caloric_state_at_specific_enthalpy(&self, h: f64,
        mut cancelled: impl FnMut() -> bool) -> Result<PureFluidIsobaricSolutionV1, E> {
        let solved = self.caloric.sample_specific_enthalpy(h, &mut cancelled)?;
        self.bound(&solved.state)?;
        poll(&mut cancelled)?;
        Ok(solved)
    }
    /// Temperature alone cannot determine enthalpy or phase amount exactly at
    /// coexistence. The original set-valued refusal remains explicit here.
    pub fn caloric_state_at_temperature(&self, temperature: Temperature,
        mut cancelled: impl FnMut() -> bool) -> Result<PureFluidIsobaricStateV1, E> {
        let state = self.caloric.sample_temperature(temperature, &mut cancelled)?;
        self.bound(&state)?;
        poll(&mut cancelled)?;
        Ok(state)
    }
    /// Publish the accepted native phase through the EXISTING complete material
    /// owner. Original NASA receipts, latent calorics and Cv remain unchanged.
    pub fn resolve_specific_enthalpy(&self, h: f64,
        mut cancelled: impl FnMut() -> bool) -> Result<PureFluidIsobaricThermalState, E> {
        let caloric = self.caloric_state_at_specific_enthalpy(h, &mut cancelled)?;
        let phase = match &caloric.state.phase {
            PureFluidIsobaricPhaseV1::Liquid(root) => PureFluidEquilibriumPhase::Liquid(root.clone()),
            PureFluidIsobaricPhaseV1::Vapor(root) => PureFluidEquilibriumPhase::Vapor(root.clone()),
            PureFluidIsobaricPhaseV1::Coexistence(state) => PureFluidEquilibriumPhase::Coexistence(state.clone()),
        };
        let material = self.model.publish_pure_equilibrium_phase(caloric.state.temperature,
            self.pressure(), caloric.state.density_kg_m3, phase).map_err(E::Material)?;
        poll(&mut cancelled)?;
        Ok(PureFluidIsobaricThermalState { material, caloric })
    }
    /// Actual direct original-source entropy, including the accepted coexistence
    /// chemical-potential defect. This is distinct from the thermal integral.
    pub fn specific_source_entropy_change_j_kg_k(&self, a: f64, b: f64,
        cancelled: impl FnMut() -> bool) -> Result<f64, E> {
        Ok(self.caloric.change_between_specific_enthalpies(a,b,cancelled)?.source_entropy_change_j_kg_k)
    }
    /// Defining caloric entropy integral for the represented thermal coordinates.
    /// Original local H must close the endpoint inverse residuals. The source
    /// sensible/latent inverse-T mean must lie in the actual endpoint interval;
    /// unresolved intervals refuse instead of clamping entropy into acceptance.
    pub fn thermal_entropy_change_j_kg_k(&self, a: f64, b: f64,
        mut cancelled: impl FnMut() -> bool) -> Result<f64, E> {
        let change = self.caloric.change_between_specific_enthalpies(a,b,&mut cancelled)?;
        self.bound(&change.initial.state)?;
        self.bound(&change.final_state.state)?;
        if a == b { poll(&mut cancelled)?; return Ok(0.0); }
        let residual = (Dd::from_f64(change.source_enthalpy_change_j_kg)
            - Dd::from_f64(change.requested_enthalpy_change_j_kg)
            - Dd::from_f64(change.final_state.residual_j_kg)
            + Dd::from_f64(change.initial.residual_j_kg)).to_f64();
        if !residual.is_finite() || residual.abs() > self.caloric.controls().enthalpy_tolerance_j_kg {
            return Err(E::Invalid("original sensible/latent H does not close the stored thermal interval"));
        }
        let dh = b-a;
        let t0=change.initial.state.temperature.value();
        let t1=change.final_state.state.temperature.value();
        if (dh>0.0 && t1<t0) || (dh<0.0 && t1>t0) {
            return Err(E::Invalid("nonmonotone pure-fluid enthalpy inverse"));
        }
        let mean = if t0==t1 { 1.0/t0 } else {
            if change.source_enthalpy_change_j_kg==0.0 || change.caloric_entropy_change_j_kg_k==0.0
                || change.source_enthalpy_change_j_kg.is_sign_positive()!=dh.is_sign_positive()
                || change.caloric_entropy_change_j_kg_k.is_sign_positive()!=dh.is_sign_positive() {
                return Err(E::Invalid("unresolved sensible/latent source entropy interval"));
            }
            let mean=change.caloric_entropy_change_j_kg_k/change.source_enthalpy_change_j_kg;
            let lower=(1.0/t0.max(t1)).next_down();
            let upper=(1.0/t0.min(t1)).next_up();
            let allowance=(64.0*f64::EPSILON*upper).next_up();
            if !positive(mean) || mean<lower-allowance || mean>upper+allowance {
                return Err(E::Invalid("source entropy mean is outside the accepted thermal interval"));
            }
            mean
        };
        let entropy=dh*mean;
        if !entropy.is_finite() || entropy==0.0 { return Err(E::Invalid("unrepresentable caloric entropy increment")); }
        poll(&mut cancelled)?;
        Ok(entropy)
    }
    fn bound(&self, state: &PureFluidIsobaricStateV1) -> Result<(), E> {
        if !state.temperature_slope_kg_k_j.is_finite() || state.temperature_slope_kg_k_j<0.0
            || state.temperature_slope_kg_k_j>self.maximum_slope || !positive(state.density_kg_m3) {
            return Err(E::Invalid("original pure-fluid response violates its global source bound"));
        }
        Ok(())
    }
}
impl ThermalLaw for PureFluidIsobaricThermalStorage<'_> {
    fn storage_basis(&self) -> ThermalStorageBasis { ThermalStorageBasis::Enthalpy }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.caloric_state_at_specific_enthalpy(h,||false)
            .map(|s|s.state.temperature.value()).map_err(|_|DepositionError::MaterialDomain)
    }
    fn bulk_density_kg_m3(&self, h: f64) -> Result<Option<f64>, DepositionError> {
        self.caloric_state_at_specific_enthalpy(h,||false)
            .map(|s|Some(s.state.density_kg_m3)).map_err(|_|DepositionError::MaterialDomain)
    }
    fn maximum_temperature_slope(&self) -> f64 { self.maximum_slope }
    fn specific_entropy_change_j_kg_k(&self, a: f64, b: f64) -> Result<Option<f64>,DepositionError> {
        self.thermal_entropy_change_j_kg_k(a,b,||false).map(Some).map_err(|_|DepositionError::MaterialDomain)
    }
}
fn positive(value:f64)->bool {value.is_finite()&&value>0.0}
fn poll(cancelled:&mut impl FnMut()->bool)->Result<(),E> {
    if cancelled(){Err(E::Cancelled)}else{Ok(())}
}

#[cfg(test)]
mod tests;
