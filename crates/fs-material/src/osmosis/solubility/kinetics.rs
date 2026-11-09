//! Actual dissolved and solid inventories under finite-rate phase exchange.
//!
//! The original thermochemical law advances each ordered neutral component;
//! only the actual dissolved amount enters shared aqueous thermodynamics. A
//! change of temperature or solvent volume preserves both phase stocks until a
//! kinetic step is taken. It does not silently impose instantaneous saturation.
//!
//! Coefficients and surface applicability remain explicit caller inputs. An
//! original law can retain a fixed surface or a self-similar particle population
//! whose area changes with actual solid stock. Dilution and transported dissolved
//! amounts never reset that law's reference particle stock.
//! `advance` is isothermal during each step. [`caloric`] instead couples these
//! actual phases to the original finite thermal bath and prescribed heat input.
//! Explicitly binding [`KineticSolubilitySolution::with_nonideal_thermodynamics`] adds
//! activity-consistent conditional phase exchange on the original virial model.
//! Neither operation supplies particle breakup/agglomeration or a nucleation law.

/// Coupled finite-rate phase exchange and the original shared-bath heat balance.
pub mod caloric;

mod nonideal;

use super::{FiniteSolubilityError, FiniteSolubilitySolution};
use crate::osmosis::{
    AqueousSolution, IdealDiluteSolution, MAX_OSMOTIC_SOLUTES, OsmosisError,
    water_activity_from_osmotic_pressure,
};
use fs_qty::{Amount, Temperature, Time, Volume, VolumetricFlowRate};
use fs_thermochem::solution::solubility::{
    SolubilityErrorV1,
    kinetics::{
        PrecipitationPolicyV1, SolubilityKineticStateV1, SolubilityKineticStepV1,
        SolubilityKineticsV1,
    },
};
use std::sync::Arc;

/// Original finite-phase/aqueous refusal or transactional cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum KineticSolubilityError {
    Material(FiniteSolubilityError),
    Cancelled,
}
impl core::fmt::Display for KineticSolubilityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-rate soluble material: {self:?}")
    }
}
impl std::error::Error for KineticSolubilityError {}
impl From<FiniteSolubilityError> for KineticSolubilityError {
    fn from(error: FiniteSolubilityError) -> Self {
        Self::Material(error)
    }
}
impl From<SolubilityErrorV1> for KineticSolubilityError {
    fn from(error: SolubilityErrorV1) -> Self {
        Self::Material(error.into())
    }
}
impl From<OsmosisError> for KineticSolubilityError {
    fn from(error: OsmosisError) -> Self {
        Self::Material(error.into())
    }
}
type Result<T> = core::result::Result<T, KineticSolubilityError>;

/// Prepared runtime stocks and original ordered kinetic/solubility sources.
/// Clones share immutable laws and own independent writable phase buffers.
#[derive(Clone, Debug, PartialEq)]
pub struct KineticSolubilitySolution {
    solution: AqueousSolution,
    models: Arc<[SolubilityKineticsV1]>,
    steps: Box<[SolubilityKineticStepV1]>,
    calibration: Option<Arc<nonideal::Calibration>>,
}
impl KineticSolubilitySolution {
    /// Prepare 1..=64 distinct neutral solutes sharing one original solvent.
    /// Initial phase stocks may be out of equilibrium; initial receipts have
    /// zero elapsed time and zero transfer.
    pub fn new(
        temperature: Temperature,
        solvent_molar_volume_m3_mol: f64,
        solvent_volume: Volume,
        components: Vec<(SolubilityKineticsV1, SolubilityKineticStateV1)>,
    ) -> Result<Self> {
        if !(1..=MAX_OSMOTIC_SOLUTES).contains(&components.len()) {
            return Err(invalid("one through 64 kinetic components required"));
        }
        for (index, (model, _)) in components.iter().enumerate() {
            if model.law().solvent() != components[0].0.law().solvent()
                || components[..index]
                    .iter()
                    .any(|(other, _)| other.law().solute() == model.law().solute())
            {
                return Err(invalid(
                    "common solvent and unique ordered kinetic solutes required",
                ));
            }
        }
        let (models, phases): (Vec<_>, Vec<_>) = components.into_iter().unzip();
        Self::prepare(
            temperature,
            solvent_molar_volume_m3_mol,
            solvent_volume,
            models.into(),
            phases.into_iter(),
        )
    }

    /// Start finite-rate evolution from the ORIGINAL equilibrium owner's stocks.
    /// Fixed-area conductances and growth policies follow its exact ordered solute
    /// layout. Use `new` with the original laws' `with_particle_surface` option
    /// when a retained reference particle stock should control evolving area.
    pub fn from_equilibrium(
        original: &FiniteSolubilitySolution,
        transfer: &[(VolumetricFlowRate, PrecipitationPolicyV1)],
    ) -> Result<Self> {
        if transfer.len() != original.components().len() {
            return Err(invalid(
                "kinetic coefficients must match the original species layout",
            ));
        }
        let components = original
            .components()
            .iter()
            .zip(original.phase_partitions())
            .zip(transfer)
            .map(|((component, phase), &(conductance, policy))| {
                Ok((
                    SolubilityKineticsV1::new(component.law().clone(), conductance, policy)?,
                    SolubilityKineticStateV1::new(phase.total_solute, phase.dissolved)?,
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        Self::new(
            Temperature::new(original.solution().temperature_k()),
            original.solution().solvent_molar_volume_m3_mol(),
            Volume::new(original.solution().solvent_volume_m3()),
            components,
        )
    }

    fn prepare(
        temperature: Temperature,
        solvent_molar_volume_m3_mol: f64,
        solvent_volume: Volume,
        models: Arc<[SolubilityKineticsV1]>,
        phases: impl ExactSizeIterator<Item = SolubilityKineticStateV1>,
    ) -> Result<Self> {
        if phases.len() != models.len() {
            return Err(invalid("original kinetic phase layout"));
        }
        let mut steps = Vec::with_capacity(models.len());
        let mut dissolved = Vec::with_capacity(models.len());
        for (model, phase) in models.iter().zip(phases) {
            let step = model.advance(temperature, solvent_volume, &phase, Time::new(0.0))?;
            dissolved.push(step.state().dissolved().value());
            steps.push(step);
        }
        let solution: AqueousSolution = IdealDiluteSolution::new(
            temperature.value(),
            solvent_molar_volume_m3_mol,
            solvent_volume.value(),
            dissolved,
        )?
        .into();
        solution.water_activity()?;
        Ok(Self {
            solution,
            models,
            steps: steps.into_boxed_slice(),
            calibration: None,
        })
    }

    /// Actual dissolved-particle view; the equilibrium target is not substituted.
    pub const fn solution(&self) -> &AqueousSolution {
        &self.solution
    }
    pub fn models(&self) -> &[SolubilityKineticsV1] {
        &self.models
    }
    /// Original common solvent identity; the admitted component list is nonempty.
    pub fn solvent(&self) -> &fs_thermochem::SpeciesId {
        self.models[0].law().solvent()
    }

    /// Owned prepared storage, conservatively counting shared source payloads.
    /// Advancing or copying a compatible prepared state allocates no memory.
    pub fn heap_bytes(&self) -> usize {
        self.models.iter().fold(
            self.solution
                .heap_bytes()
                .saturating_add(core::mem::size_of_val(self.models.as_ref()))
                .saturating_add(core::mem::size_of_val(self.steps.as_ref()))
                .saturating_add(2 * core::mem::size_of::<usize>())
                .saturating_add(core::mem::align_of::<SolubilityKineticsV1>())
                .saturating_add(
                    self.calibration
                        .as_ref()
                        .map_or(0, |data| data.heap_bytes()),
                ),
            |bytes, model| bytes.saturating_add(model.law().heap_bytes()),
        )
    }

    /// Replace transported solvent and dissolved molecules while retaining every
    /// actual solid stock. A bare membrane cannot transport or equilibrate the
    /// crystals. New component totals are derived from the two explicit phases;
    /// a trace solid is retained even when it is below one ulp of the solution.
    /// All components and the aqueous endpoint are checked before publication.
    pub fn replace_dissolved_inventories(
        &mut self,
        volume_m3: f64,
        dissolved_moles: &[f64],
    ) -> Result<()> {
        if dissolved_moles.len() != self.models.len() {
            return Err(invalid("transported dissolved inventory layout"));
        }
        let mut candidate = [None; MAX_OSMOTIC_SOLUTES];
        let temperature = Temperature::new(self.solution.temperature_k());
        self.validate_aqueous(volume_m3, dissolved_moles)?;
        let concentrations = self.phase_concentrations(volume_m3, dissolved_moles)?;
        for (index, ((model, old), &dissolved)) in self
            .models
            .iter()
            .zip(self.steps.iter())
            .zip(dissolved_moles)
            .enumerate()
        {
            let phase =
                SolubilityKineticStateV1::from_phases(Amount::new(dissolved), old.state().solid())?;
            candidate[index] = Some(if self.calibration.is_some() {
                model.advance_at_saturation(
                    temperature,
                    Volume::new(volume_m3),
                    &phase,
                    Time::new(0.0),
                    self.conditional_saturation_unchecked(index, &concentrations)?,
                )?
            } else {
                model.advance(temperature, Volume::new(volume_m3), &phase, Time::new(0.0))?
            });
        }
        self.publish_prepared(volume_m3, dissolved_moles, &candidate)
    }

    /// Copy an accepted compatible state into prepared native transport buffers.
    /// The exact phase receipts, including sub-ulp donor tails, are copied without
    /// repartitioning, re-equilibrating or allocating. Sources, temperature and
    /// solvent convention must match before any accepted storage can change.
    pub fn copy_state_from(&mut self, source: &Self) -> Result<()> {
        if !(Arc::ptr_eq(&self.models, &source.models) || self.models == source.models)
            || self.solution.temperature_k() != source.solution.temperature_k()
            || self.solution.solvent_molar_volume_m3_mol()
                != source.solution.solvent_molar_volume_m3_mol()
            || self.solution.nonideal_model() != source.solution.nonideal_model()
            || self.calibration != source.calibration
        {
            return Err(invalid("incompatible kinetic source or thermal state"));
        }
        self.solution.replace_inventories(
            source.solution.solvent_volume_m3(),
            source.solution.solute_moles(),
        )?;
        self.steps.copy_from_slice(&source.steps);
        Ok(())
    }
    /// Original species order; each receipt retains current total/dissolved/solid
    /// stocks, signed last transfer, target and any exhaustion event. With
    /// interacting solutes the target holds the other concentrations at that
    /// coordinate stage; later stages can move it. It is not a simultaneous
    /// coupled equilibrium prediction.
    pub fn last_steps(&self) -> &[SolubilityKineticStepV1] {
        &self.steps
    }

    /// Advance every component at the original common T and solvent volume.
    /// Uses only prepared buffers and bounded stack storage; no heap allocation
    /// occurs. Particle surfaces use the original bounded scalar phase solve.
    /// Cancellation is polled before work, between components and before
    /// publishing any accepted state. A bound nonideal model uses a first-order
    /// sequential phase split in the declared species order. External adaptive
    /// full/two-half stepping must control its cross-component time error.
    pub fn advance(
        &mut self,
        duration: Time,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<&[SolubilityKineticStepV1]> {
        poll(&mut cancelled)?;
        let mut candidate = [None; MAX_OSMOTIC_SOLUTES];
        let mut dissolved = [0.0; MAX_OSMOTIC_SOLUTES];
        dissolved[..self.models.len()].copy_from_slice(self.solution.solute_moles());
        let temperature = Temperature::new(self.solution.temperature_k());
        let volume = Volume::new(self.solution.solvent_volume_m3());
        let mut concentrations =
            self.phase_concentrations(volume.value(), self.solution.solute_moles())?;
        for (index, (model, old)) in self.models.iter().zip(self.steps.iter()).enumerate() {
            poll(&mut cancelled)?;
            let step = if self.calibration.is_some() {
                model.advance_at_saturation(
                    temperature,
                    volume,
                    old.state(),
                    duration,
                    self.conditional_saturation_unchecked(index, &concentrations)?,
                )?
            } else {
                model.advance(temperature, volume, old.state(), duration)?
            };
            dissolved[index] = step.state().dissolved().value();
            concentrations[index] = dissolved[index] / volume.value();
            candidate[index] = Some(step);
        }
        let count = self.models.len();
        self.validate_aqueous(volume.value(), &dissolved[..count])?;
        poll(&mut cancelled)?;
        self.publish_prepared(volume.value(), &dissolved[..count], &candidate)?;
        Ok(&self.steps)
    }

    fn validate_aqueous(&self, volume_m3: f64, dissolved: &[f64]) -> Result<()> {
        let pressure = self.solution.osmotic_pressure_at(volume_m3, dissolved)?;
        water_activity_from_osmotic_pressure(
            self.solution.temperature_k(),
            pressure,
            self.solution.solvent_molar_volume_m3_mol(),
        )?;
        Ok(())
    }

    fn publish_prepared(
        &mut self,
        volume_m3: f64,
        dissolved: &[f64],
        candidate: &[Option<SolubilityKineticStepV1>; MAX_OSMOTIC_SOLUTES],
    ) -> Result<()> {
        self.solution.replace_inventories(volume_m3, dissolved)?;
        for (accepted, prepared) in self.steps.iter_mut().zip(candidate) {
            *accepted = prepared.expect("every original kinetic component was prepared");
        }
        Ok(())
    }

    /// Prescribe pure-solvent dilution/concentration without re-equilibrating.
    /// Both actual phase stocks are unchanged; target capacities are refreshed.
    pub fn with_solvent_volume_m3(&self, volume_m3: f64) -> Result<Self> {
        if self.calibration.is_some() {
            let mut next = self.clone();
            next.solution = self.solution.with_solvent_volume_m3(volume_m3)?;
            next.advance(Time::new(0.0), || false)?;
            return Ok(next);
        }
        Self::prepare(
            Temperature::new(self.solution.temperature_k()),
            self.solution.solvent_molar_volume_m3_mol(),
            Volume::new(volume_m3),
            Arc::clone(&self.models),
            self.steps.iter().map(|step| *step.state()),
        )
    }

    /// Externally prescribe a temperature inside every original pair domain.
    /// Actual stocks remain unchanged. The caller owns heat/thermostat accounting
    /// and must ensure each declared film coefficient applies at the new
    /// temperature. Particle reference stocks and their area response are retained.
    /// An explicitly bound virial model admits only its original temperature;
    /// changing that fixed source requires new thermodynamic data.
    pub fn with_temperature(&self, temperature: Temperature) -> Result<Self> {
        if self.calibration.is_some() {
            if temperature.value().to_bits() != self.solution.temperature_k().to_bits() {
                return Err(invalid(
                    "a fixed-temperature virial source cannot be extrapolated to a new temperature",
                ));
            }
            return Ok(self.clone());
        }
        Self::prepare(
            temperature,
            self.solution.solvent_molar_volume_m3_mol(),
            Volume::new(self.solution.solvent_volume_m3()),
            Arc::clone(&self.models),
            self.steps.iter().map(|step| *step.state()),
        )
    }
}

fn invalid(message: &'static str) -> KineticSolubilityError {
    FiniteSolubilityError::Invalid(message).into()
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(KineticSolubilityError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
