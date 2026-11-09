//! Source-relative calorics between original pure-fluid coexistence states.
//!
//! The original saturation owner selects phases and the volume lever rule.
//! This owner binds a fixed density and integrates changes of the original
//! NASA/PR phase laws before summing changing phase amounts. No UV root finder,
//! phase flash, heat integrator or apparent heat capacity is introduced.
//!
//! For any phase property x, the symmetric exact product identity is
//! Delta[(1-q)x_l+q*x_v] = (1-q_bar)Delta x_l + q_bar Delta x_v
//!                              + Delta q * mean(x_v-x_l).
//! The common ideal source is added once; identical species have no mixing
//! entropy. Source entropy retains the accepted coexistence Gibbs defect.

use super::super::{MixingRecord, PengRobinsonErrorV1, PengRobinsonMixtureV1};
use super::{
    PengRobinsonEquilibriumResponseV1, PengRobinsonSaturatedMixtureV1,
    PengRobinsonSaturationErrorV1,
};
use crate::{
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, StandardStatePhaseV1,
    ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    kinetics::{KineticsErrorV1, admit_caloric_interval},
};
use fs_math::dd::Dd;
use fs_qty::{Density, Temperature};

/// Original source, phase or numerical refusal; no query mutates its owner.
#[derive(Clone, Debug, PartialEq)]
pub enum PureFluidIsochoricErrorV1 {
    Invalid(&'static str),
    Eos(PengRobinsonErrorV1),
    Source(ThermochemErrorV1),
    Caloric(KineticsErrorV1),
    Phase(PengRobinsonSaturationErrorV1),
    Cancelled,
}
impl core::fmt::Display for PureFluidIsochoricErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pure PR isochoric calorics: {self:?}")
    }
}
impl std::error::Error for PureFluidIsochoricErrorV1 {}
impl From<PengRobinsonErrorV1> for PureFluidIsochoricErrorV1 {
    fn from(value: PengRobinsonErrorV1) -> Self { Self::Eos(value) }
}
impl From<ThermochemErrorV1> for PureFluidIsochoricErrorV1 {
    fn from(value: ThermochemErrorV1) -> Self { Self::Source(value) }
}
impl From<PengRobinsonSaturationErrorV1> for PureFluidIsochoricErrorV1 {
    fn from(value: PengRobinsonSaturationErrorV1) -> Self {
        if matches!(value, PengRobinsonSaturationErrorV1::Cancelled) { Self::Cancelled }
        else { Self::Phase(value) }
    }
}
type E = PureFluidIsochoricErrorV1;

/// Direct final-minus-reference calorics and the ORIGINAL final response.
/// No pressure work belongs in this fixed-volume internal-energy increment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PureFluidIsochoricChangeV1 {
    pub specific_internal_energy_change_j_kg: f64,
    pub source_entropy_change_j_kg_k: f64,
    /// Original moving-phase Cv, pressure, entropy and acoustic response.
    pub response: PengRobinsonEquilibriumResponseV1,
    /// Mixture-volume minus declared fixed-volume, normalized by fixed volume.
    pub relative_volume_residual: f64,
    /// Original phase caloric evaluations: zero for an identical state, else four.
    /// These are explicit algebra evaluations, not additional saturation solves.
    pub caloric_evaluations: usize,
}

/// Immutable reference to one pure NASA/PR source and one accepted two-phase
/// state at declared fixed density. Every queried state must retain this same
/// density, source component and subcritical temperature interval.
///
/// Only resolved interior phase amounts are admitted. The existing material UV
/// owner can select broader phases; this prepared caloric lane deliberately
/// refuses departure from coexistence rather than silently changing its model.
pub struct PengRobinsonIsochoricCaloricReferenceV1<'a> {
    model: &'a PengRobinsonMixtureV1,
    source: &'a Nasa9StandardStateModelV1,
    initial: PengRobinsonSaturatedMixtureV1,
    density: Density,
    molar_volume: f64,
    domain: [Temperature; 2],
}
impl<'a> PengRobinsonIsochoricCaloricReferenceV1<'a> {
    pub fn new(
        model: &'a PengRobinsonMixtureV1,
        source: &'a Nasa9StandardStateModelV1,
        initial: PengRobinsonSaturatedMixtureV1,
        density: Density,
        domain: [Temperature; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, E> {
        poll(&mut cancelled)?;
        if model.components().len()!=1 || model.mole_fractions()!=[1.0]
            || model.components()[0]!=*initial.saturation().component()
            || model.molar_mass()!=source.molar_mass()
            || source.convention().phase()!=StandardStatePhaseV1::Gas
            || source.convention().eos()!=ReferenceEquationOfStateV1::IdealGas
            || !positive(density.value()) || !positive(domain[0].value())
            || !domain[1].value().is_finite() || domain[0].value()>=domain[1].value()
            || domain[1].value()>=model.components()[0].critical_temperature.value()
        {
            return Err(E::Invalid("same pure source, positive fixed density and subcritical interval required"));
        }
        let molar_volume=source.molar_mass().value()/density.value();
        if !positive(molar_volume) { return Err(E::Invalid("fixed molar volume is not representable")); }
        admit_caloric_interval(source,domain[0],domain[1],&mut cancelled).map_err(|error| {
            if matches!(error,KineticsErrorV1::Cancelled) { E::Cancelled } else { E::Caloric(error) }
        })?;
        // This is the original PR alpha coefficient, not an estimated Cv floor.
        // kappa>=0 below Tc implies a''>=0 and hence each residual Cv>=0.
        // The original moving-phase term -T*(p_sat'-p_T)^2/p_v is also >=0
        // on the two mechanically stable roots. Material owners can therefore
        // use their original analytic ideal-Cv lower bound throughout the lane.
        let mut nonnegative=false;
        model.mixing_caloric(domain[0].value(),|record| {
            if let MixingRecord::Component{kappa,..}=record { nonnegative=kappa>=0.0; }
        })?;
        if !nonnegative { return Err(E::Invalid("original alpha curvature has no nonnegative Cv proof")); }
        let result=Self{model,source,initial,density,molar_volume,domain};
        result.admit(&result.initial)?;
        poll(&mut cancelled)?;
        Ok(result)
    }
    pub const fn source(&self)->&'a Nasa9StandardStateModelV1 {self.source}
    pub const fn constitutive_law(&self)->&'a PengRobinsonMixtureV1 {self.model}
    pub const fn initial_state(&self)->&PengRobinsonSaturatedMixtureV1 {&self.initial}
    pub const fn density(&self)->Density {self.density}
    pub const fn temperature_domain(&self)->[Temperature;2] {self.domain}

    /// Complete original response at an independently accepted same-density
    /// coexistence state. No frozen-quality or ideal Cp-Cv shortcut is used.
    pub fn response_at(&self, state:&PengRobinsonSaturatedMixtureV1,
        mut cancelled:impl FnMut()->bool)->Result<PengRobinsonEquilibriumResponseV1,E> {
        poll(&mut cancelled)?;
        let (_,response)=self.admit(state)?;
        poll(&mut cancelled)?;
        Ok(response)
    }

    /// Direct local U/S increment, retaining formation/entropy datum invariance.
    /// This uses two prepared original density-caloric changes (four bounded
    /// phase evaluations) and polls between them and before publication.
    pub fn change_to(&self, final_state:&PengRobinsonSaturatedMixtureV1,
        mut cancelled:impl FnMut()->bool)->Result<PureFluidIsochoricChangeV1,E> {
        poll(&mut cancelled)?;
        let (relative_volume_residual,response)=self.admit(final_state)?;
        if final_state==&self.initial {
            poll(&mut cancelled)?;
            return Ok(PureFluidIsochoricChangeV1{
                specific_internal_energy_change_j_kg:0.0,source_entropy_change_j_kg_k:0.0,
                response,relative_volume_residual,caloric_evaluations:0,
            });
        }
        let a=self.initial.saturation();
        let b=final_state.saturation();
        let t0=a.temperature();let t1=b.temperature();
        let mut du=[Dd::ZERO;2];let mut ds=[Dd::ZERO;2];
        for (i,(before,after)) in [(a.liquid(),b.liquid()),(a.vapor(),b.vapor())].into_iter().enumerate() {
            poll(&mut cancelled)?;
            let change=self.model.prepare_density_caloric_change(t0,
                Density::new(before.mass_density_kg_per_m3()),
                Density::new(after.mass_density_kg_per_m3()))?.evaluate(t1)?;
            du[i]=Dd::from_f64(change.residual_internal_energy_change_j_mol());
            ds[i]=Dd::from_f64(change.residual_entropy_at_density_change_j_mol_k())
                +Dd::from_f64(change.ideal_entropy_density_correction_j_mol_k());
        }
        let q0=Dd::from_f64(self.initial.vapor_fraction());
        let q1=Dd::from_f64(final_state.vapor_fraction());
        let mean=(q0+q1)*Dd::from_f64(0.5);
        let dq=q1-q0;
        let latent_u=|s:&super::PengRobinsonSaturationV1| {
            Dd::from_f64(s.vapor().residual_internal_energy_j_mol())
                -Dd::from_f64(s.liquid().residual_internal_energy_j_mol())
        };
        let latent_s=|s:&super::PengRobinsonSaturationV1| {
            Dd::from_f64(s.vapor().residual_entropy_j_mol_k())
                -Dd::from_f64(s.liquid().residual_entropy_j_mol_k())
        };
        let mass=Dd::from_f64(self.source.molar_mass().value());
        let energy=(Dd::from_f64(self.source.molar_internal_energy_change(t0,t1)?.value())
            +(Dd::from_f64(1.0)-mean)*du[0]+mean*du[1]
            +dq*(latent_u(a)+latent_u(b))*Dd::from_f64(0.5))/mass;
        let entropy=(Dd::from_f64(self.source.molar_entropy_change(t0,t1)?.value())
            +(Dd::from_f64(1.0)-mean)*ds[0]+mean*ds[1]
            +dq*(latent_s(a)+latent_s(b))*Dd::from_f64(0.5))/mass;
        let specific_internal_energy_change_j_kg=represented(energy)?;
        let source_entropy_change_j_kg_k=represented(entropy)?;
        poll(&mut cancelled)?;
        Ok(PureFluidIsochoricChangeV1{
            specific_internal_energy_change_j_kg,source_entropy_change_j_kg_k,
            response,relative_volume_residual,caloric_evaluations:4,
        })
    }

    fn admit(&self,state:&PengRobinsonSaturatedMixtureV1)
        ->Result<(f64,PengRobinsonEquilibriumResponseV1),E> {
        let saturation=state.saturation();let t=saturation.temperature().value();
        if saturation.component()!=&self.model.components()[0]
            || t<self.domain[0].value() || t>self.domain[1].value()
            || !(0.0<state.vapor_fraction() && state.vapor_fraction()<1.0)
        {
            return Err(E::Invalid("same source, supported temperature and interior two-phase state required"));
        }
        let relative=((Dd::from_f64(state.molar_volume().value())
            -Dd::from_f64(self.molar_volume))/Dd::from_f64(self.molar_volume)).to_f64();
        // Same geometric allowance as the original saturation volume lever rule.
        if !relative.is_finite() || relative.abs()>128.0*f64::EPSILON {
            return Err(E::Invalid("coexistence state does not retain the declared rigid density"));
        }
        let ideal=self.source.evaluate_properties(saturation.temperature())?.cp().value()-R;
        let response=state.equilibrium_response(ideal)?;
        if !positive(ideal) || response.molar_cv_j_mol_k<ideal {
            return Err(E::Invalid("original moving-phase response violates nonnegative residual Cv"));
        }
        Ok((relative,response))
    }
}
fn positive(x:f64)->bool {x.is_finite()&&x>0.0}
fn represented(x:Dd)->Result<f64,E> {
    let y=x.to_f64();
    if !y.is_finite() || (x!=Dd::ZERO&&y==0.0) {
        Err(E::Invalid("original source increment is not representable"))
    } else {Ok(y)}
}
fn poll(cancelled:&mut impl FnMut()->bool)->Result<(),E> {
    if cancelled(){Err(E::Cancelled)}else{Ok(())}
}

#[cfg(test)]
mod tests;
