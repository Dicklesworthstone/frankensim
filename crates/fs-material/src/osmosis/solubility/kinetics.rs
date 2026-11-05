//! Actual dissolved and solid inventories under finite-rate phase exchange.
//!
//! The original thermochemical law advances each independent neutral solute;
//! only the actual dissolved amount enters shared aqueous thermodynamics. A
//! change of temperature or solvent volume preserves both phase stocks until a
//! kinetic step is taken. It does not silently impose instantaneous saturation.
//!
//! Coefficients and fixed surface applicability remain explicit caller inputs.
//! This owner is isothermal during each step, supplies no simultaneous membrane
//! or moisture transport, and infers no thermostat heat, nonideal chemistry,
//! shrinking particle geometry, or nucleation law.

use super::{FiniteSolubilityError, FiniteSolubilitySolution};
use crate::osmosis::{
    AqueousSolution, IdealDiluteSolution, MAX_OSMOTIC_SOLUTES, OsmosisError,
    water_activity_from_osmotic_pressure,
};
use fs_qty::{Temperature, Time, Volume, VolumetricFlowRate};
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
    /// Conductances and growth policies follow its exact ordered solute layout.
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
        })
    }

    /// Actual dissolved-particle view; the equilibrium target is not substituted.
    pub const fn solution(&self) -> &AqueousSolution {
        &self.solution
    }
    pub fn models(&self) -> &[SolubilityKineticsV1] {
        &self.models
    }
    /// Original species order; each receipt retains current total/dissolved/solid
    /// stocks, signed last transfer, target and any exhaustion event.
    pub fn last_steps(&self) -> &[SolubilityKineticStepV1] {
        &self.steps
    }

    /// Advance every component at the original common T and solvent volume.
    /// Uses only prepared buffers and bounded stack storage; no heap allocation
    /// or iterative phase solve occurs during this step. Cancellation is polled
    /// before work, between components and before publishing any accepted state.
    pub fn advance(
        &mut self,
        duration: Time,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<&[SolubilityKineticStepV1]> {
        poll(&mut cancelled)?;
        let mut candidate = [None; MAX_OSMOTIC_SOLUTES];
        let mut dissolved = [0.0; MAX_OSMOTIC_SOLUTES];
        let temperature = Temperature::new(self.solution.temperature_k());
        let volume = Volume::new(self.solution.solvent_volume_m3());
        for (index, (model, old)) in self.models.iter().zip(self.steps.iter()).enumerate() {
            poll(&mut cancelled)?;
            let step = model.advance(temperature, volume, old.state(), duration)?;
            dissolved[index] = step.state().dissolved().value();
            candidate[index] = Some(step);
        }
        let count = self.models.len();
        let pressure = self
            .solution
            .osmotic_pressure_at(volume.value(), &dissolved[..count])?;
        water_activity_from_osmotic_pressure(
            temperature.value(),
            pressure,
            self.solution.solvent_molar_volume_m3_mol(),
        )?;
        poll(&mut cancelled)?;
        self.solution
            .replace_inventories(volume.value(), &dissolved[..count])?;
        for (accepted, prepared) in self.steps.iter_mut().zip(candidate) {
            *accepted = prepared.expect("every original kinetic component was prepared");
        }
        Ok(&self.steps)
    }

    /// Prescribe pure-solvent dilution/concentration without re-equilibrating.
    /// Both actual phase stocks are unchanged; target capacities are refreshed.
    pub fn with_solvent_volume_m3(&self, volume_m3: f64) -> Result<Self> {
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
    /// and must ensure each fixed conductance applies at the new temperature.
    pub fn with_temperature(&self, temperature: Temperature) -> Result<Self> {
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
