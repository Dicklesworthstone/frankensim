//! Isothermal finite gas parcels coupled to the original dense-membrane owner.
//!
//! Both reservoirs and sorbed cells share one implicit amount solve. Actual gas
//! inventories, EOS, caloric law, transport and source support remain owned by
//! ClosedGasParcel. Its source-relative boundary operation receives the exact
//! isothermal thermostat term -dm*(R*T/M), canceling pressure-flow work before
//! formation references are evaluated. No rounded absolute energy difference
//! is used to recover that small thermal balance.
//!
//! These are GAS-SIDE energy accounts. D and Henry S do not supply sorbed-phase
//! internal energy or heat of sorption, so no whole-system energy closure or
//! membrane heat prediction is claimed. Temperature, geometry and pure-species
//! identity stay fixed; the two gas masses and pressures evolve.

use std::sync::Arc;
use super::{Bindings,ResolvedMembraneDiffusion,ResolvedMembraneDiffusionError,
    MembraneDiffusionSourceSupport,MembraneDiffusionState,MembraneDiffusionSolve,
    MembraneDiffusionError,GasEnergyError,source_error};
use super::slab::finite::{MolarPressureCapacity,FiniteMembraneDiffusionState,
    FiniteMembraneDiffusionWorkspace,FiniteMembraneDiffusionStep};
use crate::thermochemical_gas::ThermochemicalGasModel;
use crate::thermochemical_gas::energy::parcel::ClosedGasParcel;
use crate::thermochemical_gas::energy::parcel::transfer::species::{catalog,Sum};
use crate::thermochemical_gas::energy::parcel::transfer::species::boundary::{
    GasSpeciesBoundarySolve,GasSpeciesBoundaryStep,SpeciesBoundaryRelativeTransfer};
use fs_qty::Temperature;
use fs_thermochem::{SpeciesId,UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

#[derive(Clone,Copy,Debug,PartialEq)]
pub struct IsothermalGasMembraneSolve {
    pub membrane:MembraneDiffusionSolve,
    pub gas:[GasSpeciesBoundarySolve;2],
    /// Combined actual gas plus sorbed-species mass budget, after EOS reseeding.
    pub mass_tolerance_kg:f64,
}
impl IsothermalGasMembraneSolve {
    pub(super) fn validate(self,temperature_k:f64)->Result<(),IsothermalGasMembraneError> {
        self.membrane.validate()?;
        if !positive(self.mass_tolerance_kg) {return Err(Error::Invalid("positive combined membrane/gas mass tolerance required"));}
        for (i,gas) in self.gas.iter().enumerate() {
            gas.gas.validate().map_err(|error|gas_error(i,error))?;
            if !positive(gas.mass_tolerance_kg)||!positive(gas.energy_tolerance_j)
                ||!(gas.gas.minimum_temperature_k..=gas.gas.maximum_temperature_k).contains(&temperature_k) {
                return Err(Error::Invalid("original gas controls need positive budgets and a bracket containing the held temperature"));
            }
        }
        Ok(())
    }
}

#[derive(Clone,Debug,PartialEq)]
pub enum IsothermalGasMembraneError {
    Invalid(&'static str),
    Source(ResolvedMembraneDiffusionError),
    Gas {reservoir:usize,error:GasEnergyError},
    PressureClosure {reservoir:usize,mass_error_kg:f64,tolerance_kg:f64},
    MassBalance {residual_kg:f64,tolerance_kg:f64},
    Cancelled,
}
impl core::fmt::Display for IsothermalGasMembraneError {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {write!(f,"isothermal gas/membrane: {self:?}")}
}
impl std::error::Error for IsothermalGasMembraneError {}
impl From<ResolvedMembraneDiffusionError> for IsothermalGasMembraneError {
    fn from(error:ResolvedMembraneDiffusionError)->Self {
        if error==ResolvedMembraneDiffusionError::Cancelled {Self::Cancelled}else{Self::Source(error)}
    }
}
impl From<MembraneDiffusionError> for IsothermalGasMembraneError {
    fn from(error:MembraneDiffusionError)->Self {ResolvedMembraneDiffusionError::from(error).into()}
}
impl From<GasEnergyError> for IsothermalGasMembraneError {
    fn from(error:GasEnergyError)->Self {ResolvedMembraneDiffusionError::from(error).into()}
}
use IsothermalGasMembraneError as Error;
pub(super) fn positive(x:f64)->bool {x.is_finite()&&x>0.0}
pub(super) fn finite(x:f64)->Result<f64,Error> {
    if x.is_finite() {Ok(x)}else{Err(Error::Invalid("gas/membrane conversion is unrepresentable"))}
}
pub(super) fn poll(cancelled:&mut impl FnMut()->bool)->Result<(),Error> {
    if cancelled() {Err(Error::Cancelled)}else{Ok(())}
}
pub(super) fn gas_error(reservoir:usize,error:GasEnergyError)->Error {
    if error==GasEnergyError::Cancelled {Error::Cancelled}else{Error::Gas {reservoir,error}}
}


/// One original gas enthalpy stream and its exact isothermal flow-work offset.
pub(super) fn isothermal_stream(species:&SpeciesId,amount_mol:f64,molar_mass:f64,
    molar_enthalpy:f64,temperature:f64)->Result<(SpeciesBoundaryRelativeTransfer,f64),Error>
{
    let dm=finite(amount_mol*molar_mass)?;
    if amount_mol!=0.0&&dm==0.0 {return Err(Error::Invalid("gas/membrane species mass transfer underflowed"));}
    // Match boundary.rs grouping so the source-relative thermal target is
    // exactly zero before any arbitrarily large formation datum enters.
    let heat=finite(-dm*(R*temperature/molar_mass))?;
    let enthalpy=finite(dm*molar_enthalpy/molar_mass)?;
    Ok((SpeciesBoundaryRelativeTransfer {species:species.clone(),mass_into_gas_kg:dm,
        relative_energy_into_gas_j:heat},enthalpy))
}

#[derive(Clone,Debug)]
pub struct IsothermalGasMembrane<'a> {
    bindings:Arc<Bindings<'a>>,
    finite:FiniteMembraneDiffusionState,
    gases:[ClosedGasParcel;2],
    catalog:Arc<ThermochemicalGasModel>,
    support:MembraneDiffusionSourceSupport,
    controls:IsothermalGasMembraneSolve,
    molar_mass_kg_mol:f64,
    molar_enthalpy_j_mol:f64,
}

#[derive(Clone,Debug,PartialEq)]
pub struct IsothermalGasMembraneStep {
    pub diffusion:FiniteMembraneDiffusionStep,
    pub gas:[GasSpeciesBoundaryStep;2],
    /// Original-temperature gas partial enthalpy times signed gas mass input.
    /// These streams do not assign an enthalpy to the sorbed phase.
    pub stream_enthalpy_into_gas_j:[f64;2],
    pub thermostat_heat_into_gas_j:[f64;2],
    pub gas_energy_residual_j:f64,
    pub actual_reservoir_pressures_pa:[f64;2],
    pub pressure_mismatch_pa:[f64;2],
    pub pressure_mismatch_mass_kg:[f64;2],
    pub actual_gas_mass_changes_kg:[f64;2],
    pub sorbed_mass_change_kg:f64,
    pub total_mass_before_kg:f64,
    pub total_mass_after_kg:f64,
    /// Sum of actual gas changes and cellwise sorbed changes, retaining small increments.
    pub mass_residual_kg:f64,
    /// Independent direct audit of accepted actual totals after pressure reseeding.
    pub total_mass_residual_kg:f64,
}
#[derive(Clone,Debug)]
pub struct IsothermalGasMembraneCandidate<'a> {
    state:IsothermalGasMembrane<'a>,
    step:IsothermalGasMembraneStep,
}
impl<'a> IsothermalGasMembraneCandidate<'a> {
    pub const fn state(&self)->&IsothermalGasMembrane<'a> {&self.state}
    pub const fn step(&self)->&IsothermalGasMembraneStep {&self.step}
    pub fn into_parts(self)->(IsothermalGasMembrane<'a>,IsothermalGasMembraneStep) {(self.state,self.step)}
}

impl<'a> ResolvedMembraneDiffusion<'a> {
    pub fn with_finite_gas_reservoirs(self,gases:[ClosedGasParcel;2],controls:IsothermalGasMembraneSolve,
        cancelled:impl FnMut()->bool)->Result<IsothermalGasMembrane<'a>,Error>
    {IsothermalGasMembrane::new(self,gases,controls,cancelled)}
}

impl<'a> IsothermalGasMembrane<'a> {
    /// Consume one source-qualified membrane and two original pure ideal gases.
    /// The catalog join checks original caloric/transport/pack identity; no gas
    /// source or mixture is inferred from a material or species name.
    pub fn new(membrane:ResolvedMembraneDiffusion<'a>,gases:[ClosedGasParcel;2],
        controls:IsothermalGasMembraneSolve,mut cancelled:impl FnMut()->bool)->Result<Self,Error>
    {
        poll(&mut cancelled)?;controls.validate(membrane.temperature_k())?;
        for gas in &gases {
            let components=gas.model().thermodynamics().components();
            if components.len()!=1||components[0].species()!=membrane.species()
                ||gas.model().thermodynamics().mole_composition().fractions()!=[1.0]
                ||gas.state().state().temperature!=membrane.temperature_k() {
                return Err(Error::Invalid("finite membrane reservoirs require the same pure species and exact common held temperature"));
            }
        }
        let catalog=Arc::new(catalog(&gases[0],&gases[1])?);
        let source=&catalog.thermodynamics().components()[0];
        let molar_mass_kg_mol=source.molar_mass().value();
        let molar_enthalpy_j_mol=source.evaluate_properties(Temperature::new(membrane.temperature_k()))
            .map_err(source_error)?.h().value();
        let rt=finite(R*membrane.temperature_k())?;
        let capacities=[MolarPressureCapacity::try_new(gases[0].volume_m3()/rt)?,
            MolarPressureCapacity::try_new(gases[1].volume_m3()/rt)?];
        let pressures=gases.each_ref().map(|gas|gas.state().state().pressure);
        let ResolvedMembraneDiffusion {bindings,state,..}=membrane;
        let flux=state.sample_fluxes(pressures,&mut cancelled)?;
        let support=bindings.admit([(&state,flux.face_partial_pressures_pa.as_slice())],&mut cancelled)?;
        let finite_state=FiniteMembraneDiffusionState::new(state,capacities,pressures,&mut cancelled)?;
        let result=Self {bindings,finite:finite_state,gases,catalog,support,controls,molar_mass_kg_mol,molar_enthalpy_j_mol};
        result.admit_capacity_mass()?;result.total_mass_kg()?;
        poll(&mut cancelled)?;Ok(result)
    }
    pub fn species(&self)->&SpeciesId {&self.bindings.species}
    pub fn temperature_k(&self)->f64 {self.bindings.temperature_k}
    pub fn membrane(&self)->&MembraneDiffusionState {self.finite.membrane()}
    pub const fn gases(&self)->&[ClosedGasParcel;2] {&self.gases}
    pub fn catalog(&self)->&ThermochemicalGasModel {&self.catalog}
    pub const fn finite_state(&self)->&FiniteMembraneDiffusionState {&self.finite}
    pub const fn source_support(&self)->&MembraneDiffusionSourceSupport {&self.support}
    pub const fn controls(&self)->IsothermalGasMembraneSolve {self.controls}
    pub fn total_mass_kg(&self)->Result<f64,Error> {
        let mut mass=Sum::default();
        mass.add(self.gases[0].mass_kg())?;
        mass.add(finite(self.membrane().stored_amount_mol()?*self.molar_mass_kg_mol)?)?;
        mass.add(self.gases[1].mass_kg())?;Ok(mass.value)
    }
    fn pressure_budget(&self,i:usize)->f64 {
        self.controls.mass_tolerance_kg.min(self.controls.gas[i].mass_tolerance_kg)
    }
    fn admit_capacity_mass(&self)->Result<(),Error> {
        let amounts=self.finite.reservoir_amounts_mol()?;let mut residual=Sum::default();
        for i in 0..2 {
            let difference=finite(self.gases[i].mass_kg()-finite(amounts[i]*self.molar_mass_kg_mol)?)?;
            if difference.abs()>self.pressure_budget(i) {
                return Err(Error::PressureClosure {reservoir:i,mass_error_kg:difference,tolerance_kg:self.pressure_budget(i)});
            }
            residual.add(difference)?;
        }
        self.admit_mass(residual.value)
    }
    fn admit_mass(&self,residual_kg:f64)->Result<(),Error> {
        if !residual_kg.is_finite()||residual_kg.abs()>self.controls.mass_tolerance_kg {
            Err(Error::MassBalance {residual_kg,tolerance_kg:self.controls.mass_tolerance_kg})
        }else{Ok(())}
    }

    /// Advance the finite inventories, then reconstruct both gases through the
    /// original source-relative boundary owner. All physical/source and gas-side
    /// energy checks precede publication; no accepted state is mutated.
    pub fn step_candidate(&self,duration_s:f64,work:&mut FiniteMembraneDiffusionWorkspace,
        mut cancelled:impl FnMut()->bool)->Result<IsothermalGasMembraneCandidate<'a>,Error>
    {
        poll(&mut cancelled)?;
        let candidate=self.finite.step_candidate(duration_s,self.controls.membrane,work,&mut cancelled)?;
        let (numerical,diffusion)=candidate.into_parts();
        let amounts=[-diffusion.membrane.left_reservoir_to_membrane_mol,diffusion.membrane.membrane_to_right_reservoir_mol];
        let mut gas_rows=Vec::with_capacity(2);let mut gases=Vec::with_capacity(2);
        let mut stream=[0.0;2];let mut heat=[0.0;2];let mut pressure_error=[0.0;2];let mut mass_error=[0.0;2];
        for i in 0..2 {
            poll(&mut cancelled)?;
            let (transfer,enthalpy)=isothermal_stream(self.species(),amounts[i],self.molar_mass_kg_mol,
                self.molar_enthalpy_j_mol,self.temperature_k())?;
            heat[i]=transfer.relative_energy_into_gas_j;stream[i]=enthalpy;
            let proposed=self.gases[i].exchange_species_boundary_relative(&self.catalog,
                &[transfer],self.controls.gas[i],&mut cancelled)
                .map_err(|error|gas_error(i,error))?;
            let (gas,_,step)=proposed.into_parts();
            if gas.state().state().temperature!=self.temperature_k() {
                return Err(Error::Invalid("original gas reconstruction changed the held membrane temperature"));
            }
            pressure_error[i]=finite(gas.state().state().pressure-numerical.reservoir_pressures_pa()[i])?;
            mass_error[i]=finite((numerical.capacities()[i].value_si()*pressure_error[i])*self.molar_mass_kg_mol)?;
            if mass_error[i].abs()>self.pressure_budget(i) {
                return Err(Error::PressureClosure {reservoir:i,mass_error_kg:mass_error[i],tolerance_kg:self.pressure_budget(i)});
            }
            gas_rows.push(step);gases.push(gas);
        }
        let gases:[ClosedGasParcel;2]=gases.try_into().map_err(|_|Error::Invalid("finite membrane needs exactly two accepted gases"))?;
        let gas:[GasSpeciesBoundaryStep;2]=gas_rows.try_into().map_err(|_|Error::Invalid("finite membrane needs exactly two gas reports"))?;
        let pressures=gases.each_ref().map(|gas|gas.state().state().pressure);
        let old=self.membrane().sample_fluxes(self.finite.reservoir_pressures_pa(),&mut cancelled)?;
        let actual=numerical.membrane().sample_fluxes(pressures,&mut cancelled)?;
        let support=self.bindings.admit([
            (self.membrane(),old.face_partial_pressures_pa.as_slice()),
            (numerical.membrane(),diffusion.membrane.face_partial_pressures_pa.as_slice()),
            (numerical.membrane(),actual.face_partial_pressures_pa.as_slice())],&mut cancelled)?;
        let capacities=numerical.capacities();
        let finite_state=FiniteMembraneDiffusionState::new(numerical.into_membrane(),capacities,pressures,&mut cancelled)?;
        let state=Self {bindings:self.bindings.clone(),finite:finite_state,gases,catalog:self.catalog.clone(),support,
            controls:self.controls,molar_mass_kg_mol:self.molar_mass_kg_mol,molar_enthalpy_j_mol:self.molar_enthalpy_j_mol};
        // Audit the ACTUAL original gas inventories after the coordinate reseed.
        state.admit_capacity_mass()?;
        let gas_changes=[state.gases[0].mass_kg()-self.gases[0].mass_kg(),state.gases[1].mass_kg()-self.gases[1].mass_kg()];
        let sorbed_change=finite(diffusion.membrane.stored_change_mol*self.molar_mass_kg_mol)?;
        let mut balance=Sum::default();balance.add(gas_changes[0])?;balance.add(gas_changes[1])?;balance.add(sorbed_change)?;
        state.admit_mass(balance.value)?;
        let before=self.total_mass_kg()?;let after=state.total_mass_kg()?;
        state.admit_mass(after-before)?;
        let gas_energy_residual_j=finite(gas[0].energy_residual_j+gas[1].energy_residual_j)?;
        let step=IsothermalGasMembraneStep {diffusion,gas,stream_enthalpy_into_gas_j:stream,
            thermostat_heat_into_gas_j:heat,gas_energy_residual_j,actual_reservoir_pressures_pa:pressures,
            pressure_mismatch_pa:pressure_error,pressure_mismatch_mass_kg:mass_error,
            actual_gas_mass_changes_kg:gas_changes,sorbed_mass_change_kg:sorbed_change,
            total_mass_before_kg:before,total_mass_after_kg:after,mass_residual_kg:balance.value,
            total_mass_residual_kg:after-before};
        poll(&mut cancelled)?;Ok(IsothermalGasMembraneCandidate {state,step})
    }
}

