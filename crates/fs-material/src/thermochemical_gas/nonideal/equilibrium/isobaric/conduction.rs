//! Original sensible/latent calorics on the existing implicit heat operator.
//!
//! The authored bulk k(T) curve explicitly assumes quality-independent
//! conductivity on the isothermal interval. Neither phase morphology nor an
//! average of gas and liquid transport cards is inferred here.

use super::PureFluidIsobaricThermalStorage;
use crate::thermal::conduction::ConductivityCurve;
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw, ImplicitConductionObservables, ImplicitConductionState,
};

/// One exact borrowed storage owner plus an explicit bulk conductivity law.
pub struct PureFluidIsobaricConductionLaw<'storage, 'source> {
    caloric: &'storage PureFluidIsobaricThermalStorage<'source>,
    conductivity: ConductivityCurve,
}
impl<'storage, 'source> PureFluidIsobaricConductionLaw<'storage, 'source> {
    /// Require complete numerical support of the original temperature span.
    /// Applicability of k(T), including its plateau assumption, is caller data.
    pub fn try_new(caloric: &'storage PureFluidIsobaricThermalStorage<'source>,
        conductivity: ConductivityCurve) -> Result<Self, DepositionError> {
        let [lo,hi]=caloric.temperature_domain();
        let [kl,kh]=conductivity.temperature_range_k();
        if kl>lo.value() || kh<hi.value() { return Err(DepositionError::MaterialDomain); }
        Ok(Self {caloric,conductivity})
    }
    pub const fn caloric(&self)->&'storage PureFluidIsobaricThermalStorage<'source> {self.caloric}
    pub const fn conductivity(&self)->&ConductivityCurve {&self.conductivity}
    /// Original h inverse and chain rule, with exact zero dT/dh and dk/dh on
    /// the inclusive latent interval. Those endpoint tangents are one-sided.
    pub fn state_at_specific_enthalpy(&self,h:f64)->Result<ImplicitConductionState,DepositionError> {
        let state=self.caloric.caloric_state_at_specific_enthalpy(h,||false)
            .map_err(|_|DepositionError::MaterialDomain)?.state;
        let temperature_k=state.temperature.value();
        let temperature_slope_kg_k_j=state.temperature_slope_kg_k_j;
        let (conductivity_w_m_k,dkdt)=self.conductivity
            .conductivity_and_temperature_derivative(temperature_k)?;
        let conductivity_slope_kg_w_m_k_j=dkdt*temperature_slope_kg_k_j;
        if !conductivity_slope_kg_w_m_k_j.is_finite()
            || (dkdt!=0.0 && temperature_slope_kg_k_j!=0.0 && conductivity_slope_kg_w_m_k_j==0.0) {
            return Err(DepositionError::Arithmetic);
        }
        Ok(ImplicitConductionState {temperature_k,temperature_slope_kg_k_j,
            conductivity_w_m_k,conductivity_slope_kg_w_m_k_j})
    }
}
impl ConductivityLaw for PureFluidIsobaricConductionLaw<'_, '_> {
    fn conductivity_w_m_k(&self,h:f64)->Result<f64,DepositionError> {
        self.conductivity.conductivity_w_m_k(self.caloric.temperature_k(h)?)
    }
}
impl ImplicitConductionLaw for PureFluidIsobaricConductionLaw<'_, '_> {
    fn caloric_law(&self)->&dyn ThermalLaw {self.caloric}
    fn conductivity_law(&self)->&dyn ConductivityLaw {self}
    fn enthalpy_bounds_j_kg(&self)->[f64;2] {self.caloric.specific_enthalpy_domain_j_kg()}
    /// Finite stepping uses the latent one-sided tangent at h=0 and h=L.
    /// A harmonic perturbation crosses into a sensible branch there, so its
    /// two-sided derivative is unavailable. Interior latent states remain valid.
    fn admit_harmonic_operating_point(&self,h:f64)->Result<(),DepositionError> {
        let [liquid,vapor]=self.caloric.coexistence_specific_enthalpy_domain_j_kg();
        if h==liquid || h==vapor { return Err(DepositionError::MaterialDomain); }
        self.state_at_specific_enthalpy(h).map(|_|())
    }
    fn state_at_specific_enthalpy(&self,h:f64)->Result<ImplicitConductionState,DepositionError> {
        Self::state_at_specific_enthalpy(self,h)
    }
    fn specific_volume_derivative_at_specific_enthalpy(&self,h:f64)->Result<Option<f64>,DepositionError> {
        self.caloric.caloric_state_at_specific_enthalpy(h,||false)
            .map(|s|Some(s.state.specific_volume_slope_m3_j)).map_err(|_|DepositionError::MaterialDomain)
    }
    fn time_observables(&self,h:f64)->Result<ImplicitConductionObservables,DepositionError> {
        let state=self.caloric.caloric_state_at_specific_enthalpy(h,||false)
            .map_err(|_|DepositionError::MaterialDomain)?.state;
        Ok(ImplicitConductionObservables {temperature_k:state.temperature.value(),
            phase_mass_fraction:Some(state.vapor_fraction())})
    }
}
