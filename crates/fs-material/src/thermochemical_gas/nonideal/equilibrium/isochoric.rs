//! Sealed rigid pure-fluid coexistence on the original finite thermal timeline.
//!
//! Conserved mass and volume fix density. Internal energy selects T, saturation
//! pressure and phase amount through the EXISTING UV/phase owner. Unlike an
//! isobaric latent interval, this equilibrium path has positive dT/du and rising
//! saturation pressure. Latent energy is already inside U; no pressure work or
//! second heat store is added. This lane stays inside a declared two-phase
//! subcritical interval and refuses phase departure, kinetics and interfaces.

pub mod conduction;

use super::{
    E, PureFluidEquilibriumControls, PureFluidEquilibriumEnergySolution,
    PureFluidEquilibriumPhase, PureFluidEquilibriumState, Target,
    ThermochemicalPengRobinsonModel, Work, failure, poll, positive,
};
use crate::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_math::dd::Dd;
use fs_qty::{Density, Temperature};
use fs_thermal::enthalpy::{DepositionError, DepositedCell, ThermalLaw, ThermalStorageBasis};
use fs_thermochem::eos::peng_robinson::saturation::{
    PengRobinsonEquilibriumResponseV1, PengRobinsonSaturatedMixtureV1,
    isochoric::{PengRobinsonIsochoricCaloricReferenceV1, PureFluidIsochoricErrorV1},
};

pub(super) fn caloric_error(error:PureFluidIsochoricErrorV1)->E {
    if matches!(error,PureFluidIsochoricErrorV1::Cancelled) {E::Cancelled}
    else {E::Material(failure("original fixed-density phase calorics",error))}
}

/// Complete original material/UV solution and same-source moving-phase response.
#[derive(Clone,Debug,PartialEq)]
pub struct PureFluidIsochoricThermalState {
    solution:PureFluidEquilibriumEnergySolution,
    response:PengRobinsonEquilibriumResponseV1,
    temperature_slope_kg_k_j:f64,
    pressure_slope_pa_kg_j:f64,
}
impl PureFluidIsochoricThermalState {
    pub const fn material(&self)->&PureFluidEquilibriumState {self.solution.state()}
    pub const fn solution(&self)->&PureFluidEquilibriumEnergySolution {&self.solution}
    pub const fn response(&self)->&PengRobinsonEquilibriumResponseV1 {&self.response}
    pub const fn temperature_slope_kg_k_j(&self)->f64 {self.temperature_slope_kg_k_j}
    pub const fn pressure_slope_pa_kg_j(&self)->f64 {self.pressure_slope_pa_kg_j}
}

/// Immutable, repeatable internal-energy law at one exact declared density.
///
/// u=0 is the prepared cold endpoint; the original absolute formation datum
/// remains in full material states. A thermal cell must retain m=rho*V within
/// the original floating-point volume allowance. Geometry cannot silently
/// implement an isobaric expansion or a second boundary-work debit.
pub struct PureFluidIsochoricThermalStorage<'a> {
    model:&'a ThermochemicalPengRobinsonModel,
    reference:PengRobinsonIsochoricCaloricReferenceV1<'a>,
    controls:PureFluidEquilibriumControls,
    endpoints:[PureFluidIsochoricThermalState;2],
    maximum_u:f64,
    maximum_slope:f64,
    preparation_temperature_evaluations:usize,
    preparation_eos_evaluations:usize,
}

impl ThermochemicalPengRobinsonModel {
    pub fn prepare_isochoric_thermal_storage(
        &self,density:Density,controls:PureFluidEquilibriumControls,
        mut cancelled:impl FnMut()->bool,
    )->Result<PureFluidIsochoricThermalStorage<'_>,E> {
        poll(&mut cancelled)?;
        self.admit_pure_equilibrium(density.value(),controls)?;
        let mut work=Work::default();
        let cold=self.equilibrium_sample(density.value(),controls.temperature_bracket[0],
            controls,&mut work,&mut cancelled)?;
        let reference=PengRobinsonIsochoricCaloricReferenceV1::new(self.mechanics(),
            &self.source().thermodynamics().components()[0],phase(&cold)?.clone(),
            density,controls.temperature_bracket,&mut cancelled).map_err(caloric_error)?;
        let branch=self.source().prepare_caloric_branch(GasEnergySolve {
            minimum_temperature_k:controls.temperature_bracket[0].value(),
            maximum_temperature_k:controls.temperature_bracket[1].value(),
            energy_tolerance_j_kg:controls.energy_tolerance_j_kg,
            maximum_iterations:controls.maximum_temperature_evaluations.min(128),
        },&mut cancelled).map_err(|error| {
            if matches!(error,GasEnergyError::Cancelled) {E::Cancelled}
            else {E::Material(failure("original isochoric capacity bound",error))}
        })?;
        let minimum_cv=branch.minimum_cv_j_kg_k();
        let maximum_slope=(1.0/minimum_cv).next_up();
        if !positive(minimum_cv)||!positive(maximum_slope) {
            return Err(E::Invalid("original source has no represented global Cv floor"));
        }
        let hot=self.equilibrium_sample(density.value(),controls.temperature_bracket[1],
            controls,&mut work,&mut cancelled)?;
        admit_nested_dome(phase(&cold)?,phase(&hot)?)?;
        if controls.maximum_eos_evaluations.saturating_sub(work.eos)<4 {
            return Err(work.exhausted());
        }
        work.eos+=4;
        let change=reference.change_to(phase(&hot)?,&mut cancelled).map_err(caloric_error)?;
        let maximum_u=change.specific_internal_energy_change_j_kg;
        if !positive(maximum_u) {
            return Err(E::Invalid("represented two-phase internal energy must rise over the interval"));
        }
        let zero=PureFluidEquilibriumEnergySolution {
            state:cold,residual_j_kg:0.0,accepted_increment_j_kg:Some(0.0),
            temperature_evaluations:0,eos_evaluations:0,
        };
        let upper=PureFluidEquilibriumEnergySolution {
            state:hot,residual_j_kg:0.0,accepted_increment_j_kg:Some(maximum_u),
            temperature_evaluations:0,eos_evaluations:0,
        };
        let cold=resolved(&reference,maximum_slope,zero,&mut cancelled)?;
        let hot=resolved(&reference,maximum_slope,upper,&mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(PureFluidIsochoricThermalStorage{
            model:self,reference,controls,endpoints:[cold,hot],maximum_u,maximum_slope,
            preparation_temperature_evaluations:work.temperatures,
            preparation_eos_evaluations:work.eos,
        })
    }
}
impl<'a> PureFluidIsochoricThermalStorage<'a> {
    pub const fn model(&self)->&'a ThermochemicalPengRobinsonModel {self.model}
    pub const fn original_calorics(&self)->&PengRobinsonIsochoricCaloricReferenceV1<'a> {&self.reference}
    pub const fn controls(&self)->PureFluidEquilibriumControls {self.controls}
    pub fn density(&self)->Density {self.reference.density()}
    pub const fn temperature_domain(&self)->[Temperature;2] {self.controls.temperature_bracket}
    pub fn specific_internal_energy_domain_j_kg(&self)->[f64;2] {[0.0,self.maximum_u]}
    pub const fn preparation_temperature_evaluations(&self)->usize {self.preparation_temperature_evaluations}
    pub const fn preparation_eos_evaluations(&self)->usize {self.preparation_eos_evaluations}

    /// Original UV inversion with the cold source reference fixed across every
    /// query. All root trials keep the original coexistence/domain/work checks.
    /// Cached endpoints are exact accepted original states, with zero new work.
    pub fn resolve_specific_internal_energy(&self,u:f64,
        mut cancelled:impl FnMut()->bool)->Result<PureFluidIsochoricThermalState,E> {
        poll(&mut cancelled)?;
        self.admit_specific_internal_energy(u)?;
        if u==0.0||u==self.maximum_u {
            let result=self.endpoints[usize::from(u!=0.0)].clone();
            poll(&mut cancelled)?;
            return Ok(result);
        }
        let solution=self.model.solve_pure_energy(self.density().value(),
            Target::Isochoric{reference:&self.reference,requested:u},self.controls,
            Work::default(),&mut cancelled)?;
        let result=resolved(&self.reference,self.maximum_slope,solution,&mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Original source-integrated entropy, including its accepted numerical
    /// coexistence defect. No ideal mixing entropy is invented between phases.
    pub fn specific_source_entropy_change_j_kg_k(&self,a:f64,b:f64,
        mut cancelled:impl FnMut()->bool)->Result<f64,E> {
        let initial=self.resolve_specific_internal_energy(a,&mut cancelled)?;
        let final_state=self.resolve_specific_internal_energy(b,&mut cancelled)?;
        let local=self.local_reference(&initial,&mut cancelled)?;
        let result=local.change_to(phase(final_state.material())?,&mut cancelled)
            .map_err(caloric_error)?.source_entropy_change_j_kg_k;
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Stored-coordinate entropy uses the actual source inverse-temperature
    /// mean only after direct local U closes both original UV residuals.
    /// Same-T/nonzero-U or out-of-interval source means refuse; they are never
    /// clamped or reinterpreted as an isothermal latent plateau.
    pub fn thermal_entropy_change_j_kg_k(&self,a:f64,b:f64,
        mut cancelled:impl FnMut()->bool)->Result<f64,E> {
        let initial=self.resolve_specific_internal_energy(a,&mut cancelled)?;
        if a==b {poll(&mut cancelled)?;return Ok(0.0);}
        let final_state=self.resolve_specific_internal_energy(b,&mut cancelled)?;
        let local=self.local_reference(&initial,&mut cancelled)?;
        let change=local.change_to(phase(final_state.material())?,&mut cancelled)
            .map_err(caloric_error)?;
        let requested=Dd::from_f64(b)-Dd::from_f64(a);
        let closure=(Dd::from_f64(change.specific_internal_energy_change_j_kg)-requested
            -Dd::from_f64(final_state.solution.residual_j_kg())
            +Dd::from_f64(initial.solution.residual_j_kg())).to_f64();
        if !closure.is_finite()||closure.abs()>self.controls.energy_tolerance_j_kg {
            return Err(E::Invalid("local source U does not close the stored rigid-cell interval"));
        }
        let du=requested.to_f64();
        let t0=initial.material().temperature_k();let t1=final_state.material().temperature_k();
        let source_u=change.specific_internal_energy_change_j_kg;
        let source_s=change.source_entropy_change_j_kg_k;
        if t0==t1 || (du>0.0&&t1<t0)||(du<0.0&&t1>t0)
            || source_u==0.0||source_s==0.0
            || source_u.is_sign_positive()!=du.is_sign_positive()
            || source_s.is_sign_positive()!=du.is_sign_positive() {
            return Err(E::Invalid("unresolved original isochoric temperature/entropy interval"));
        }
        let mean=source_s/source_u;
        let lower=(1.0/t0.max(t1)).next_down();
        let upper=(1.0/t0.min(t1)).next_up();
        let allowance=(64.0*f64::EPSILON*upper).next_up();
        if !positive(mean)||mean<lower-allowance||mean>upper+allowance {
            return Err(E::Invalid("original isochoric entropy mean is outside its endpoint temperatures"));
        }
        let result=du*mean;
        if !result.is_finite()||result==0.0 {
            return Err(E::Invalid("rigid-cell entropy increment is not representable"));
        }
        poll(&mut cancelled)?;
        Ok(result)
    }
    // Density and rigid-volume derivatives are fixed by the conservation
    // constraint. Their domain admission needs no repeated UV inversion.
    fn admit_specific_internal_energy(&self,u:f64)->Result<(),E> {
        if !u.is_finite()||u<0.0||u>self.maximum_u {
            return Err(E::Invalid("specific internal energy outside the sealed two-phase interval"));
        }
        Ok(())
    }
    fn local_reference(&self,state:&PureFluidIsochoricThermalState,
        cancelled:&mut impl FnMut()->bool)->Result<PengRobinsonIsochoricCaloricReferenceV1<'a>,E> {
        PengRobinsonIsochoricCaloricReferenceV1::new(self.model.mechanics(),
            self.reference.source(),phase(state.material())?.clone(),self.density(),
            self.controls.temperature_bracket,cancelled).map_err(caloric_error)
    }
}
fn phase(state:&PureFluidEquilibriumState)->Result<&PengRobinsonSaturatedMixtureV1,E> {
    let PureFluidEquilibriumPhase::Coexistence(phase)=state.phase() else {
        return Err(E::Invalid("rigid thermal storage requires original interior phase coexistence"));
    };
    if !(0.0<phase.vapor_fraction()&&phase.vapor_fraction()<1.0) {
        return Err(E::Invalid("pure-phase endpoints need a different caloric domain"));
    }
    Ok(phase)
}
fn resolved(reference:&PengRobinsonIsochoricCaloricReferenceV1<'_>,maximum_slope:f64,
    solution:PureFluidEquilibriumEnergySolution,cancelled:&mut impl FnMut()->bool)
    ->Result<PureFluidIsochoricThermalState,E> {
    let response=reference.response_at(phase(solution.state())?,cancelled).map_err(caloric_error)?;
    let temperature_slope_kg_k_j=reference.source().molar_mass().value()/response.molar_cv_j_mol_k;
    let pressure_slope_pa_kg_j=response.dp_dt_v_pa_k*temperature_slope_kg_k_j;
    if !positive(temperature_slope_kg_k_j)||temperature_slope_kg_k_j>maximum_slope
        || !positive(pressure_slope_pa_kg_j) {
        return Err(E::Invalid("original moving-phase response violates the global source Cv bound"));
    }
    Ok(PureFluidIsochoricThermalState {
        solution,response,temperature_slope_kg_k_j,pressure_slope_pa_kg_j,
    })
}
fn admit_nested_dome(cold:&PengRobinsonSaturatedMixtureV1,
    hot:&PengRobinsonSaturatedMixtureV1)->Result<(),E> {
    // For this ORIGINAL pure PR alpha below Tc, a'<=0 and
    // d(p_T)/dv=-R/(v-b)^2 + a'*(2v+2b)/(v^2+2bv-b^2)^2 <0.
    // Exact coexistence makes p_sat' the mean of p_T across the two volumes:
    // p_T,v < p_sat' < p_T,l. Both p_v are negative, so saturated liquid
    // volume rises and vapor volume falls with T. The hot interior volume is
    // therefore inside every colder dome in this admitted branch.
    // Observe the same signs/order on the accepted original numerical roots;
    // every UV trial STILL runs native phase and rigid-volume admission.
    let a=cold.saturation();let b=hot.saturation();
    for s in [a,b] {
        let slope=s.saturation_pressure_slope_pa_k();
        let (ptl,pvl)=s.liquid().pressure_derivatives();
        let (ptv,pvv)=s.vapor().pressure_derivatives();
        if !(ptv<slope&&slope<ptl&&pvl<0.0&&pvv<0.0) {
            return Err(E::Invalid("original saturation volume nesting is numerically unresolved"));
        }
    }
    if !(a.liquid().molar_volume_m3_per_mol()<b.liquid().molar_volume_m3_per_mol()
        && a.vapor().molar_volume_m3_per_mol()>b.vapor().molar_volume_m3_per_mol()
        && a.pressure().value()<b.pressure().value()) {
        return Err(E::Invalid("original subcritical coexistence interval is not nested"));
    }
    Ok(())
}
impl ThermalLaw for PureFluidIsochoricThermalStorage<'_> {
    fn storage_basis(&self)->ThermalStorageBasis {ThermalStorageBasis::InternalEnergy}
    fn validate_cell(&self,cell:&DepositedCell)->Result<(),DepositionError> {
        if !positive(cell.mass_kg)||!positive(cell.geometry.volume_m3) {
            return Err(DepositionError::MaterialDomain);
        }
        let relative=((Dd::from_f64(self.density().value())*Dd::from_f64(cell.geometry.volume_m3)
            -Dd::from_f64(cell.mass_kg))/Dd::from_f64(cell.mass_kg)).to_f64();
        if !relative.is_finite()||relative.abs()>128.0*f64::EPSILON {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(())
    }
    fn temperature_k(&self,u:f64)->Result<f64,DepositionError> {
        self.resolve_specific_internal_energy(u,||false).map(|s|s.material().temperature_k())
            .map_err(|_|DepositionError::MaterialDomain)
    }
    fn bulk_density_kg_m3(&self,u:f64)->Result<Option<f64>,DepositionError> {
        self.admit_specific_internal_energy(u).map(|_|Some(self.density().value()))
            .map_err(|_|DepositionError::MaterialDomain)
    }
    fn maximum_temperature_slope(&self)->f64 {self.maximum_slope}
    fn specific_entropy_change_j_kg_k(&self,a:f64,b:f64)->Result<Option<f64>,DepositionError> {
        self.thermal_entropy_change_j_kg_k(a,b,||false).map(Some).map_err(|_|DepositionError::MaterialDomain)
    }
}

#[cfg(test)]
mod tests;
