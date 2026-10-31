//! Finite gas-mixture separation by explicitly independent Henry/Fick channels.
//!
//! Each catalog species has its own sourced D, S and optional interfacial
//! permeance on ONE shared rigid laminate. Flux is relative to that stationary
//! polymer and is driven by species partial pressure; it is not a zero-sum
//! barycentric gas-diffusion law. The original finite-reservoir solve advances
//! each independent channel, then ONE original gas reconstruction per endpoint
//! updates all species, EOS and transport together.
//!
//! Independence is an explicit constitutive assumption. The caller must supply
//! original claims applicable to this mixture; pure-gas measurements are not
//! promoted to mixed-gas evidence. Composition, humidity and concentration
//! dependent/restricted claims refuse in this rung. Fixed material, specimen,
//! process and protocol identities remain the original bindings. Competitive
//! sorption, cross diffusion, swelling and heat of sorption are not inferred.
//!
//! Source applicability uses summed gas-equivalent Henry pressures at actual
//! cells and ordered faces, including the actual gas EOS total pressure at each
//! exterior boundary. Partial pressures remain the transport coordinates.
//! Every catalog species needs an explicit channel; D=0 preserves blocked gas
//! and sorbed inventories without inventing an absent species law. Temperature
//! is externally held. Energy reports cover the two original gases and their
//! explicit thermostats, not an unknown sorbed-phase internal-energy store.

use std::sync::Arc;
use super::{Bindings,ResolvedMembraneDiffusion,MembraneDiffusionLayerSource,
    MembraneInterfaceMaterial,MembraneDiffusionSourceSupport};
use super::slab::finite::{MolarPressureCapacity,FiniteMembraneDiffusionState,
    FiniteMembraneDiffusionWorkspace,FiniteMembraneDiffusionStep};
use super::reservoir::{IsothermalGasMembraneSolve,IsothermalGasMembraneError as Error,
    finite,positive,poll,gas_error,isothermal_stream};
use crate::state_point::ResolvedScalarProperty;
use crate::thermochemical_gas::{ThermochemicalGasModel,
    energy::parcel::{ClosedGasParcel,transfer::species::{catalog,inventory,Sum,
        boundary::{GasSpeciesBoundaryStep,SpeciesBoundaryRelativeTransfer}}}};
use fs_qty::{Temperature,semantic::QuantityKind};
use fs_thermochem::{SpeciesId,UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    mixture::MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1};

/// Existing original source bindings plus an explicit sorbed loading per cell.
/// All channels must share thicknesses and cell counts in the same layer order.
#[derive(Clone,Debug)]
pub struct MembraneMixtureChannelSource<'a> {
    pub layers:Vec<MembraneDiffusionLayerSource<'a>>,
    pub interfaces:Vec<Option<MembraneInterfaceMaterial<'a>>>,
    pub initial_partial_pressures_pa:Vec<f64>,
}

#[derive(Clone,Debug)]
struct Channel<'a> {
    bindings:Arc<Bindings<'a>>,
    finite:FiniteMembraneDiffusionState,
    molar_mass:f64,
    molar_enthalpy:f64,
}
impl Channel<'_> {
    fn species(&self)->&SpeciesId {&self.bindings.species}
}

#[derive(Clone,Debug)]
pub struct IsothermalGasMembraneMixture<'a> {
    channels:Vec<Channel<'a>>,
    gases:[ClosedGasParcel;2],
    catalog:Arc<ThermochemicalGasModel>,
    support:Vec<MembraneDiffusionSourceSupport>,
    controls:IsothermalGasMembraneSolve,
}

/// Each channel retains the original factor cache keyed by its exact model,
/// both reservoir capacities and duration. No species matrix is refactored
/// merely because the shared gas composition changed.
#[derive(Debug)]
pub struct IsothermalGasMembraneMixtureWorkspace {
    channels:Vec<FiniteMembraneDiffusionWorkspace>,
}

#[derive(Clone,Debug,PartialEq)]
pub struct MembraneMixtureSpeciesStep {
    pub species:SpeciesId,
    pub diffusion:FiniteMembraneDiffusionStep,
    pub actual_partial_pressures_pa:[f64;2],
    pub pressure_mismatch_pa:[f64;2],
    pub pressure_mismatch_mass_kg:[f64;2],
    pub stream_enthalpy_into_gas_j:[f64;2],
    pub thermostat_heat_into_gas_j:[f64;2],
    pub actual_gas_mass_changes_kg:[f64;2],
    pub sorbed_mass_change_kg:f64,
    pub total_mass_before_kg:f64,
    pub total_mass_after_kg:f64,
    pub mass_residual_kg:f64,
    pub total_mass_residual_kg:f64,
}
#[derive(Clone,Debug,PartialEq)]
pub struct IsothermalGasMembraneMixtureStep {
    /// Canonical SpeciesId order, independent of authored channel ordering.
    pub species:Vec<MembraneMixtureSpeciesStep>,
    pub gas:[GasSpeciesBoundaryStep;2],
    pub actual_reservoir_pressures_pa:[f64;2],
    pub stream_enthalpy_into_gas_j:[f64;2],
    pub thermostat_heat_into_gas_j:[f64;2],
    pub gas_energy_residual_j:f64,
    pub actual_gas_mass_changes_kg:[f64;2],
    pub sorbed_mass_change_kg:f64,
    pub total_mass_before_kg:f64,
    pub total_mass_after_kg:f64,
    pub mass_residual_kg:f64,
    pub total_mass_residual_kg:f64,
}
#[derive(Clone,Debug)]
pub struct IsothermalGasMembraneMixtureCandidate<'a> {
    state:IsothermalGasMembraneMixture<'a>,
    step:IsothermalGasMembraneMixtureStep,
}
impl<'a> IsothermalGasMembraneMixtureCandidate<'a> {
    pub const fn state(&self)->&IsothermalGasMembraneMixture<'a> {&self.state}
    pub const fn step(&self)->&IsothermalGasMembraneMixtureStep {&self.step}
    pub fn into_parts(self)->(IsothermalGasMembraneMixture<'a>,IsothermalGasMembraneMixtureStep) {(self.state,self.step)}
}

fn admit_independent(property:&ResolvedScalarProperty)->Result<(),Error> {
    // The shared Henry admission already requires static fixed conditions and
    // canonical T/p. These additional gas conditions no longer stay fixed when
    // a mixture enriches, so a constant scalar with such validity also refuses.
    let receipt=&property.answer().receipt;
    if property.answer().evidence.model.validity.bounds().keys().any(|axis|
        receipt.axis_quantities.get(axis).and_then(|q|q.semantic_type()).is_some_and(|q|
            matches!(q.kind(),QuantityKind::Composition(_)|QuantityKind::Moisture(_)
                |QuantityKind::MassConcentration|QuantityKind::AmountConcentration))) {
        return Err(Error::Invalid("independent membrane mixture cannot freeze composition, humidity or concentration source applicability"));
    }
    Ok(())
}
fn admit_bindings(bindings:&Bindings<'_>)->Result<(),Error> {
    for layer in &bindings.layers {
        for key in [super::MEMBRANE_DIFFUSIVITY_PROPERTY,super::MOLAR_HENRY_SOLUBILITY_PROPERTY] {
            admit_independent(layer.original.property(key)
                .ok_or(Error::Invalid("mixture channel lost its original D or S claim"))?)?;
        }
    }
    for interface in bindings.interfaces.iter().flatten() {
        admit_independent(interface.original.property(super::MOLAR_INTERFACE_PERMEANCE_PROPERTY)
            .ok_or(Error::Invalid("mixture channel lost its original interface claim"))?)?;
    }
    Ok(())
}

struct PressureSample {
    cells:Vec<f64>,
    faces:Vec<[Option<f64>;2]>,
}
impl PressureSample {
    fn collect(states:&[&FiniteMembraneDiffusionState],boundary_totals:[f64;2],
        cancelled:&mut impl FnMut()->bool)->Result<Self,Error>
    {
        let count=states[0].membrane().model().cell_count();
        let mut cells:Vec<Sum>=(0..count).map(|_|Sum::default()).collect();
        let mut faces:Vec<[Sum;2]>=(0..=count).map(|_|[Sum::default(),Sum::default()]).collect();
        let mut known=vec![[true;2];count+1];
        for state in states {
            poll(cancelled)?;
            let flux=state.membrane().sample_fluxes(state.reservoir_pressures_pa(),&mut *cancelled)?;
            for (i,&q) in state.membrane().partial_pressures_pa().iter().enumerate() {
                if i%128==0 {poll(cancelled)?;}cells[i].add(q)?;
            }
            for (i,trace) in flux.face_partial_pressures_pa.iter().enumerate() {
                if i%128==0 {poll(cancelled)?;}
                for side in 0..2 {
                    if let Some(q)=trace[side] {faces[i][side].add(q)?;}else{known[i][side]=false;}
                }
            }
        }
        let mut faces:Vec<_>=faces.into_iter().enumerate().map(|(i,row)| {
            let [left,right]=row;
            [known[i][0].then_some(left.value),known[i][1].then_some(right.value)]
        }).collect();
        // The exterior half-cell contact is Dirichlet on both traces, including
        // a sealed D=0 channel. Keep the original EOS total exactly at a source
        // boundary instead of perturbing it by summing rounded mole fractions.
        faces[0]=[Some(boundary_totals[0]);2];
        faces[count]=[Some(boundary_totals[1]);2];
        Ok(Self {cells:cells.into_iter().map(|sum|sum.value).collect(),faces})
    }
}

impl<'a> IsothermalGasMembraneMixture<'a> {
    pub fn new(temperature_k:f64,area_m2:f64,sources:&[MembraneMixtureChannelSource<'a>],
        gases:[ClosedGasParcel;2],controls:IsothermalGasMembraneSolve,
        mut cancelled:impl FnMut()->bool)->Result<Self,Error>
    {
        poll(&mut cancelled)?;controls.validate(temperature_k)?;
        if !positive(temperature_k)||sources.is_empty()||sources.len()>MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1
            ||gases.iter().any(|gas|gas.state().state().temperature!=temperature_k) {
            return Err(Error::Invalid("bounded mixture channels and both original gases must share one exact positive held temperature"));
        }
        let catalog=Arc::new(catalog(&gases[0],&gases[1])?);
        if sources.len()!=catalog.thermodynamics().components().len() {
            return Err(Error::Invalid("every catalog species requires one explicit membrane channel, including blocked species"));
        }
        let rt=finite(R*temperature_k)?;
        let capacities=[MolarPressureCapacity::try_new(gases[0].volume_m3()/rt)?,
            MolarPressureCapacity::try_new(gases[1].volume_m3()/rt)?];
        let mut channels=Vec::with_capacity(sources.len());
        for source in sources {
            poll(&mut cancelled)?;
            let (bindings,state)=ResolvedMembraneDiffusion::resolve_parts(temperature_k,area_m2,
                &source.layers,&source.interfaces,&source.initial_partial_pressures_pa,&mut cancelled)?;
            admit_bindings(&bindings)?;
            let gas_source=catalog.thermodynamics().components().iter()
                .find(|law|law.species()==&bindings.species)
                .ok_or(Error::Invalid("membrane species is absent from the original gas catalog"))?;
            let molar_mass=gas_source.molar_mass().value();
            let molar_enthalpy=gas_source.evaluate_properties(Temperature::new(temperature_k))
                .map_err(super::source_error)?.h().value();
            let partial=[gases[0].species_partial_pressure_pa(&bindings.species)?,
                gases[1].species_partial_pressure_pa(&bindings.species)?];
            let finite=FiniteMembraneDiffusionState::new(state,capacities,partial,&mut cancelled)?;
            channels.push(Channel {bindings,finite,molar_mass,molar_enthalpy});
        }
        channels.sort_by(|a,b|a.species().cmp(b.species()));
        let first=channels[0].finite.membrane();
        for (i,channel) in channels.iter().enumerate() {
            let state=channel.finite.membrane();
            if i>0&&channel.species()==channels[i-1].species()
                ||state.elapsed_s()!=first.elapsed_s()
                ||state.model().layers().len()!=first.model().layers().len()
                ||state.model().layers().iter().zip(first.model().layers()).any(|(a,b)|
                    a.cells!=b.cells||a.thickness_m!=b.thickness_m) {
                return Err(Error::Invalid("unique mixture species must share the exact layer geometry, grid and time"));
            }
        }
        let samples:Vec<_>=channels.iter().map(|channel|&channel.finite).collect();
        let total=PressureSample::collect(&samples,gases.each_ref().map(|gas|gas.state().state().pressure),&mut cancelled)?;
        let mut support=Vec::with_capacity(channels.len());
        for channel in &channels {
            poll(&mut cancelled)?;
            support.push(channel.bindings.admit_pressures([
                (channel.finite.membrane(),total.cells.as_slice(),total.faces.as_slice())],&mut cancelled)?);
        }
        let state=Self {channels,gases,catalog,support,controls};
        state.admit_capacity_mass()?;state.total_mass_kg()?;poll(&mut cancelled)?;Ok(state)
    }
    pub fn temperature_k(&self)->f64 {self.channels[0].bindings.temperature_k}
    pub fn elapsed_s(&self)->f64 {self.channels[0].finite.membrane().elapsed_s()}
    pub const fn gases(&self)->&[ClosedGasParcel;2] {&self.gases}
    pub fn catalog(&self)->&ThermochemicalGasModel {&self.catalog}
    pub const fn controls(&self)->IsothermalGasMembraneSolve {self.controls}
    pub fn species(&self)->impl ExactSizeIterator<Item=&SpeciesId> {self.channels.iter().map(Channel::species)}
    fn index(&self,species:&SpeciesId)->Option<usize> {self.channels.binary_search_by(|row|row.species().cmp(species)).ok()}
    pub fn finite_state(&self,species:&SpeciesId)->Option<&FiniteMembraneDiffusionState> {
        self.index(species).map(|i|&self.channels[i].finite)
    }
    pub fn source_support(&self,species:&SpeciesId)->Option<&MembraneDiffusionSourceSupport> {
        self.index(species).map(|i|&self.support[i])
    }
    pub fn workspace(&self)->Result<IsothermalGasMembraneMixtureWorkspace,Error> {
        Ok(IsothermalGasMembraneMixtureWorkspace {channels:self.channels.iter()
            .map(|row|FiniteMembraneDiffusionWorkspace::new(row.finite.membrane().model().cell_count()))
            .collect::<Result<_,_>>()?})
    }
    pub fn total_mass_kg(&self)->Result<f64,Error> {
        let mut mass=Sum::default();mass.add(self.gases[0].mass_kg())?;mass.add(self.gases[1].mass_kg())?;
        for channel in &self.channels {mass.add(finite(channel.finite.membrane().stored_amount_mol()?*channel.molar_mass)?)?;}
        Ok(mass.value)
    }
    fn admit_mass(&self,residual_kg:f64)->Result<(),Error> {
        if !residual_kg.is_finite()||residual_kg.abs()>self.controls.mass_tolerance_kg {
            Err(Error::MassBalance {residual_kg,tolerance_kg:self.controls.mass_tolerance_kg})
        }else{Ok(())}
    }
    fn pressure_budget(&self,side:usize)->f64 {self.controls.mass_tolerance_kg.min(self.controls.gas[side].mass_tolerance_kg)}
    fn admit_capacity_mass(&self)->Result<(),Error> {
        let inventories=[inventory(self.gases[0].model(),self.gases[0].mass_kg())?,
            inventory(self.gases[1].model(),self.gases[1].mass_kg())?];
        let mut total=Sum::default();
        for channel in &self.channels {
            let amounts=channel.finite.reservoir_amounts_mol()?;
            for side in 0..2 {
                let actual=inventories[side].get(channel.species()).copied().unwrap_or(0.0);
                let error=finite(actual-finite(amounts[side]*channel.molar_mass)?)?;
                if error.abs()>self.pressure_budget(side) {
                    return Err(Error::PressureClosure {reservoir:side,mass_error_kg:error,tolerance_kg:self.pressure_budget(side)});
                }
                total.add(error.abs())?;
            }
        }
        self.admit_mass(total.value)
    }

    /// All channels start from the same original gases. Exactly one original
    /// source-relative reconstruction per gas accepts their combined streams.
    /// Failure at any stage leaves this owner and all its inventories unchanged.
    pub fn step_candidate(&self,duration_s:f64,work:&mut IsothermalGasMembraneMixtureWorkspace,
        mut cancelled:impl FnMut()->bool)->Result<IsothermalGasMembraneMixtureCandidate<'a>,Error>
    {
        poll(&mut cancelled)?;
        if work.channels.len()!=self.channels.len() {
            return Err(Error::Invalid("mixture workspace must cover every explicit species channel"));
        }
        let mut numerical=Vec::with_capacity(self.channels.len());
        let mut diffusion=Vec::with_capacity(self.channels.len());
        let mut transfers:[Vec<SpeciesBoundaryRelativeTransfer>;2]=[Vec::new(),Vec::new()];
        let mut streams=Vec::with_capacity(self.channels.len());let mut heats=Vec::with_capacity(self.channels.len());
        for (channel,workspace) in self.channels.iter().zip(&mut work.channels) {
            poll(&mut cancelled)?;
            let (state,report)=channel.finite.step_candidate(duration_s,self.controls.membrane,workspace,&mut cancelled)?.into_parts();
            let amounts=[-report.membrane.left_reservoir_to_membrane_mol,report.membrane.membrane_to_right_reservoir_mol];
            let mut stream=[0.0;2];let mut heat=[0.0;2];
            for side in 0..2 {
                let (transfer,enthalpy)=isothermal_stream(channel.species(),amounts[side],channel.molar_mass,
                    channel.molar_enthalpy,self.temperature_k())?;
                stream[side]=enthalpy;heat[side]=transfer.relative_energy_into_gas_j;transfers[side].push(transfer);
            }
            numerical.push(state);diffusion.push(report);streams.push(stream);heats.push(heat);
        }
        let mut gases=Vec::with_capacity(2);let mut gas_reports=Vec::with_capacity(2);
        for side in 0..2 {
            poll(&mut cancelled)?;
            let candidate=self.gases[side].exchange_species_boundary_relative(&self.catalog,&transfers[side],
                self.controls.gas[side],&mut cancelled).map_err(|error|gas_error(side,error))?;
            let (gas,_,report)=candidate.into_parts();
            if gas.state().state().temperature!=self.temperature_k() {
                return Err(Error::Invalid("original mixture reconstruction changed the held membrane temperature"));
            }
            gases.push(gas);gas_reports.push(report);
        }
        let gases:[ClosedGasParcel;2]=gases.try_into().map_err(|_|Error::Invalid("mixture needs two gas endpoints"))?;
        let gas:[GasSpeciesBoundaryStep;2]=gas_reports.try_into().map_err(|_|Error::Invalid("mixture needs two gas reports"))?;
        let before_inventory=[inventory(self.gases[0].model(),self.gases[0].mass_kg())?,
            inventory(self.gases[1].model(),self.gases[1].mass_kg())?];
        let after_inventory=[inventory(gases[0].model(),gases[0].mass_kg())?,inventory(gases[1].model(),gases[1].mass_kg())?];
        let old_refs:Vec<_>=self.channels.iter().map(|row|&row.finite).collect();
        let old=PressureSample::collect(&old_refs,self.gases.each_ref().map(|gas|gas.state().state().pressure),&mut cancelled)?;
        let mut total_pressure=[Sum::default(),Sum::default()];
        for state in &numerical {for side in 0..2 {total_pressure[side].add(state.reservoir_pressures_pa()[side])?;}}
        let numerical_refs:Vec<_>=numerical.iter().collect();
        let solved=PressureSample::collect(&numerical_refs,total_pressure.map(|sum|sum.value),&mut cancelled)?;
        let mut channels=Vec::with_capacity(self.channels.len());let mut rows=Vec::with_capacity(self.channels.len());
        let mut storage=Sum::default();let mut stream_total=[Sum::default(),Sum::default()];let mut heat_total=[Sum::default(),Sum::default()];
        for (i,(numerical,diffusion)) in numerical.into_iter().zip(diffusion).enumerate() {
            poll(&mut cancelled)?;let channel=&self.channels[i];let id=channel.species();
            let partial=[gases[0].species_partial_pressure_pa(id)?,gases[1].species_partial_pressure_pa(id)?];
            let mut pressure_error=[0.0;2];let mut mass_error=[0.0;2];let mut changes=[0.0;2];
            let mut mass_before=Sum::default();let mut mass_after=Sum::default();let mut balance=Sum::default();
            for side in 0..2 {
                pressure_error[side]=finite(partial[side]-numerical.reservoir_pressures_pa()[side])?;
                mass_error[side]=finite((numerical.capacities()[side].value_si()*pressure_error[side])*channel.molar_mass)?;
                if mass_error[side].abs()>self.pressure_budget(side) {
                    return Err(Error::PressureClosure {reservoir:side,mass_error_kg:mass_error[side],tolerance_kg:self.pressure_budget(side)});
                }
                let before=before_inventory[side].get(id).copied().unwrap_or(0.0);
                let after=after_inventory[side].get(id).copied().unwrap_or(0.0);
                changes[side]=finite(after-before)?;balance.add(changes[side])?;
                mass_before.add(before)?;mass_after.add(after)?;
                stream_total[side].add(streams[i][side])?;heat_total[side].add(heats[i][side])?;
            }
            let stored_change=finite(diffusion.membrane.stored_change_mol*channel.molar_mass)?;
            balance.add(stored_change)?;storage.add(stored_change)?;self.admit_mass(balance.value)?;
            mass_before.add(finite(channel.finite.membrane().stored_amount_mol()?*channel.molar_mass)?)?;
            mass_after.add(finite(numerical.membrane().stored_amount_mol()?*channel.molar_mass)?)?;
            self.admit_mass(mass_after.value-mass_before.value)?;
            let capacities=numerical.capacities();
            let finite_state=FiniteMembraneDiffusionState::new(numerical.into_membrane(),capacities,partial,&mut cancelled)?;
            channels.push(Channel {bindings:channel.bindings.clone(),finite:finite_state,
                molar_mass:channel.molar_mass,molar_enthalpy:channel.molar_enthalpy});
            rows.push(MembraneMixtureSpeciesStep {species:id.clone(),diffusion,actual_partial_pressures_pa:partial,
                pressure_mismatch_pa:pressure_error,pressure_mismatch_mass_kg:mass_error,
                stream_enthalpy_into_gas_j:streams[i],thermostat_heat_into_gas_j:heats[i],
                actual_gas_mass_changes_kg:changes,sorbed_mass_change_kg:stored_change,
                total_mass_before_kg:mass_before.value,total_mass_after_kg:mass_after.value,
                mass_residual_kg:balance.value,total_mass_residual_kg:mass_after.value-mass_before.value});
        }
        let actual_refs:Vec<_>=channels.iter().map(|row|&row.finite).collect();
        let pressures=gases.each_ref().map(|gas|gas.state().state().pressure);
        let actual=PressureSample::collect(&actual_refs,pressures,&mut cancelled)?;
        let mut support=Vec::with_capacity(channels.len());
        for channel in &channels {
            poll(&mut cancelled)?;let state=channel.finite.membrane();
            support.push(channel.bindings.admit_pressures([
                (state,old.cells.as_slice(),old.faces.as_slice()),
                (state,solved.cells.as_slice(),solved.faces.as_slice()),
                (state,actual.cells.as_slice(),actual.faces.as_slice())],&mut cancelled)?);
        }
        let state=Self {channels,gases,catalog:self.catalog.clone(),support,controls:self.controls};
        state.admit_capacity_mass()?;
        let changes=[state.gases[0].mass_kg()-self.gases[0].mass_kg(),state.gases[1].mass_kg()-self.gases[1].mass_kg()];
        let mut balance=Sum::default();balance.add(changes[0])?;balance.add(changes[1])?;balance.add(storage.value)?;
        state.admit_mass(balance.value)?;
        let before=self.total_mass_kg()?;let after=state.total_mass_kg()?;state.admit_mass(after-before)?;
        let gas_energy_residual_j=finite(gas[0].energy_residual_j+gas[1].energy_residual_j)?;
        let step=IsothermalGasMembraneMixtureStep {species:rows,gas,actual_reservoir_pressures_pa:pressures,
            stream_enthalpy_into_gas_j:stream_total.map(|sum|sum.value),
            thermostat_heat_into_gas_j:heat_total.map(|sum|sum.value),gas_energy_residual_j,
            actual_gas_mass_changes_kg:changes,sorbed_mass_change_kg:storage.value,
            total_mass_before_kg:before,total_mass_after_kg:after,mass_residual_kg:balance.value,total_mass_residual_kg:after-before};
        poll(&mut cancelled)?;Ok(IsothermalGasMembraneMixtureCandidate {state,step})
    }
}
