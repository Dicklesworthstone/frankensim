//! Explicit bulk heat conduction on the same sealed two-phase internal energy.
//!
//! k(T) is authored for this declared fixed-density equilibrium path. Neither
//! gas transport nor phase morphology supplies an implicit mixture rule here.

use super::PureFluidIsochoricThermalStorage;
use crate::thermal::conduction::ConductivityCurve;
use fs_thermal::enthalpy::{DepositionError,ThermalLaw};
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw,ImplicitConductionObservables,ImplicitConductionState,
};

pub struct PureFluidIsochoricConductionLaw<'storage,'source> {
    caloric:&'storage PureFluidIsochoricThermalStorage<'source>,
    conductivity:ConductivityCurve,
}
impl<'storage,'source> PureFluidIsochoricConductionLaw<'storage,'source> {
    pub fn try_new(caloric:&'storage PureFluidIsochoricThermalStorage<'source>,
        conductivity:ConductivityCurve)->Result<Self,DepositionError> {
        let [lo,hi]=caloric.temperature_domain();
        let [kl,kh]=conductivity.temperature_range_k();
        if kl>lo.value()||kh<hi.value() {return Err(DepositionError::MaterialDomain);}
        Ok(Self{caloric,conductivity})
    }
    pub const fn caloric(&self)->&'storage PureFluidIsochoricThermalStorage<'source> {self.caloric}
    pub const fn conductivity(&self)->&ConductivityCurve {&self.conductivity}
    /// Trait h-labelled coordinates carry INTERNAL ENERGY for this law.
    pub fn state_at_specific_internal_energy(&self,u:f64)
        ->Result<ImplicitConductionState,DepositionError> {
        let state=self.caloric.resolve_specific_internal_energy(u,||false)
            .map_err(|_|DepositionError::MaterialDomain)?;
        let temperature_k=state.material().temperature_k();
        let temperature_slope_kg_k_j=state.temperature_slope_kg_k_j();
        let (conductivity_w_m_k,dkdt)=self.conductivity
            .conductivity_and_temperature_derivative(temperature_k)?;
        let conductivity_slope_kg_w_m_k_j=dkdt*temperature_slope_kg_k_j;
        if !conductivity_slope_kg_w_m_k_j.is_finite()
            || (dkdt!=0.0&&conductivity_slope_kg_w_m_k_j==0.0) {
            return Err(DepositionError::Arithmetic);
        }
        Ok(ImplicitConductionState{temperature_k,temperature_slope_kg_k_j,
            conductivity_w_m_k,conductivity_slope_kg_w_m_k_j})
    }
}
impl ConductivityLaw for PureFluidIsochoricConductionLaw<'_, '_> {
    fn conductivity_w_m_k(&self,u:f64)->Result<f64,DepositionError> {
        self.conductivity.conductivity_w_m_k(self.caloric.temperature_k(u)?)
    }
}
impl ImplicitConductionLaw for PureFluidIsochoricConductionLaw<'_, '_> {
    fn caloric_law(&self)->&dyn ThermalLaw {self.caloric}
    fn conductivity_law(&self)->&dyn ConductivityLaw {self}
    fn enthalpy_bounds_j_kg(&self)->[f64;2] {self.caloric.specific_internal_energy_domain_j_kg()}
    fn state_at_specific_enthalpy(&self,u:f64)->Result<ImplicitConductionState,DepositionError> {
        self.state_at_specific_internal_energy(u)
    }
    fn admit_harmonic_operating_point(&self,u:f64)->Result<(),DepositionError> {
        let t=self.state_at_specific_internal_energy(u)?.temperature_k;
        // The fluid lane is strictly inside coexistence. A declared k(T) kink
        // still lacks a two-sided harmonic derivative; one-sided stepping stays valid.
        if self.conductivity.knots().windows(3).any(|k|k[1][0]==t
            &&(k[1][1]-k[0][1])/(k[1][0]-k[0][0])!=(k[2][1]-k[1][1])/(k[2][0]-k[1][0])) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(())
    }
    fn specific_volume_derivative_at_specific_enthalpy(&self,u:f64)
        ->Result<Option<f64>,DepositionError> {
        self.caloric.admit_specific_internal_energy(u)
            .map(|_|Some(0.0)).map_err(|_|DepositionError::MaterialDomain)
    }
    fn time_observables(&self,u:f64)->Result<ImplicitConductionObservables,DepositionError> {
        let state=self.caloric.resolve_specific_internal_energy(u,||false)
            .map_err(|_|DepositionError::MaterialDomain)?;
        Ok(ImplicitConductionObservables{temperature_k:state.material().temperature_k(),
            phase_mass_fraction:Some(state.material().vapor_mass_fraction())})
    }
}
