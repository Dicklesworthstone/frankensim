//! Composition-dependent radioactive spectra on the original heat transaction.
//!
//! Prepared per-parent energy partitions preserve the original nuclear energy
//! account. Only the radioactive part passes through the supplied group response;
//! other outward energy, such as prompt capture/fission radiation, escapes and is
//! reported explicitly. No attenuation coefficient or material response is inferred.

use super::*;
use crate::decay::thermal::DecayCellSource;
use crate::endf::EndfError;
use crate::radiation::decay::{BurnupDecayRadiation, DecayRadiationStep, PreparedDecayRadiation};
use crate::radiation::grouped::GroupedRadiationResponse;

/// A prepared evaluated energy partition for one original physical source cell.
#[derive(Clone, Debug)]
pub struct DecayRadiationEmitter {
    pub cell_id: CellId,
    pub partition: Arc<PreparedDecayRadiation>,
}

/// Original source laws whose accepted reports expose evaluated decay events.
/// Sealed through NuclearHeatSource; no source inventory is copied or re-evolved.
pub trait DecayRadiationHeatSource: NuclearHeatSource + 'static {
    /// Immutable, original-model-bound projection shared by trial clones.
    type RadiationBinding: core::fmt::Debug + Send + Sync + 'static;

    fn bind_decay_radiation(
        &self,
        partition: Arc<PreparedDecayRadiation>,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Self::RadiationBinding, EndfError>;

    fn decay_radiation_step(
        &self,
        binding: &Self::RadiationBinding,
        step: &Self::Step,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError>;
}

trait Projection<S: NuclearHeatSource>: core::fmt::Debug + Send + Sync {
    fn project(
        &self,
        source: &S,
        step: &S::Step,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError>;
}

#[derive(Debug)]
struct PreparedProjection<S: DecayRadiationHeatSource> {
    binding: S::RadiationBinding,
}

impl<S: DecayRadiationHeatSource> Projection<S> for PreparedProjection<S> {
    fn project(
        &self,
        source: &S,
        step: &S::Step,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        source.decay_radiation_step(&self.binding, step, checkpoint)
    }
}

/// This immutable binding contains no mutable nuclear or thermal history.
#[derive(Debug)]
pub(super) struct BoundSpectralRadiation<S: NuclearHeatSource> {
    response: Arc<GroupedRadiationResponse>,
    projections: BTreeMap<CellId, Arc<dyn Projection<S>>>,
    receivers: BTreeMap<CellId, ThermalBinding>,
}

impl<S: DecayRadiationHeatSource> NuclearHeating<S> {
    /// Apply an explicit energy-dependent response to evaluated radioactive energy.
    ///
    /// Emitters must cover the original cohort exactly and their immutable models,
    /// radiation kind and energy grid must agree. Each prepared partition already
    /// checks spectral completeness and explicitly reconciles evaluated mean-energy
    /// rounding to the original radioactive energy account. Other outward energy
    /// remains ungrouped escape; it is never assigned a guessed decay spectrum.
    ///
    /// This replaces any previous aggregate or grouped response at an accepted
    /// boundary. The material/geometry bindings, source trials, heat credits and
    /// cancellation all use the existing NuclearHeating transaction. Response
    /// coefficients must be independently supplied for their declared applicability.
    pub fn with_grouped_decay_radiation(
        mut self,
        thermal: &DepositionSimulation<'_>,
        response: Arc<GroupedRadiationResponse>,
        emitters: Vec<DecayRadiationEmitter>,
        receivers: Vec<RadiationReceiver>,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, FissionThermalError> {
        checkpoint()?;
        self.validate(thermal, self.time_s, &mut || Ok(()))?;
        if self.pending.is_some() {
            return Err(FissionThermalError::Invalid(
                "cannot change grouped radiation during a pending interval",
            ));
        }
        if emitters.len() != self.sources.len()
            || !response
                .source_ids()
                .eq(self.sources.keys().map(|cell| cell.0))
        {
            return Err(FissionThermalError::Invalid(
                "grouped radiation must cover the original source cohort exactly",
            ));
        }
        let mut projections: BTreeMap<CellId, Arc<dyn Projection<S>>> = BTreeMap::new();
        for emitter in emitters {
            checkpoint()?;
            let source = self
                .sources
                .get(&emitter.cell_id)
                .ok_or(FissionThermalError::Invalid(
                    "unknown grouped radiation source cell",
                ))?;
            if emitter.partition.radiation() != response.radiation()
                || emitter.partition.group_edges_j() != response.group_edges_j()
            {
                return Err(FissionThermalError::Invalid(
                    "radiation kind and energy grid must match the evaluated partition",
                ));
            }
            let binding = source
                .source
                .bind_decay_radiation(emitter.partition, &mut checkpoint)?;
            if projections
                .insert(
                    emitter.cell_id,
                    Arc::new(PreparedProjection::<S> { binding }),
                )
                .is_some()
            {
                return Err(FissionThermalError::Invalid(
                    "duplicate grouped radiation emitter",
                ));
            }
        }
        let receivers =
            radiation::bind_receivers(thermal, receivers, response.receiver_ids(), &mut || {
                checkpoint().map_err(Into::into)
            })?;
        checkpoint()?;
        self.radiation = None;
        self.photon_feedback = None;
        self.spectral_radiation = Some(Arc::new(BoundSpectralRadiation {
            response,
            projections,
            receivers,
        }));
        Ok(self)
    }
}

impl<S: NuclearHeatSource> NuclearHeating<S> {
    /// The grouped response used for future accepted source intervals.
    pub fn grouped_radiation_response(&self) -> Option<&Arc<GroupedRadiationResponse>> {
        self.spectral_radiation
            .as_ref()
            .map(|binding| &binding.response)
    }
}

impl<S: NuclearHeatSource> BoundSpectralRadiation<S> {
    // The shared temperature model retains identical photon/source/receiver
    // grids. Replacing only its response shares all original decay projections
    // and keeps the native receiver correspondence unchanged.
    pub(super) fn with_response(&self, response: Arc<GroupedRadiationResponse>) -> Self {
        Self {
            response,
            projections: self.projections.clone(),
            receivers: self.receivers.clone(),
        }
    }

    pub(super) fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), FissionThermalError> {
        for (&id, binding) in &self.receivers {
            checkpoint()?;
            binding.validate(thermal, id)?;
            thermal.temperature_k(id)?;
        }
        Ok(())
    }

    pub(super) fn project(
        &self,
        id: CellId,
        source: &S,
        report: &S::Step,
        account: &mut Account,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<f64>, FissionThermalError> {
        let mut checkpoint_error = None;
        let emission = self.projections[&id].project(source, report, &mut || {
            checkpoint().map_err(|error| {
                checkpoint_error = Some(error);
                FissionError::Cancelled
            })
        });
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        let emission = emission?;
        if emission.radiation != self.response.radiation()
            || emission.group_edges_j.as_ref() != self.response.group_edges_j()
            || !nonnegative(emission.original_radiation_j)
            || !nonnegative(emission.reconciled_spectral_energy_j)
        {
            return Err(FissionThermalError::Invalid(
                "grouped radioactive source changed its binding",
            ));
        }
        // Compare the projected radioactive energy with the original interval
        // account. The remaining prompt/non-decay energy has no supplied spectrum.
        let mut remainder = Sum::default();
        remainder.add(S::released(report).radiation_j)?;
        remainder.add(-emission.original_radiation_j)?;
        let remainder_j = remainder.value()?;
        account
            .radiation_roundoff
            .add(remainder.export_roundoff()?)?;
        if remainder_j < 0.0 {
            // A small numerical overshoot is charged, never hidden by clipping;
            // the shared advance-wide audit rejects it if the allowance is spent.
            account.radiation_roundoff.add(-remainder_j)?;
        } else {
            account.radiation_escaped.add(remainder_j)?;
            account.radiation_ungrouped.add(remainder_j)?;
        }
        let mut closure = Sum::default();
        closure.add(-emission.original_radiation_j)?;
        for &energy in &emission.emitted_energy_j {
            checkpoint()?;
            if !nonnegative(energy) {
                return Err(FissionThermalError::Invalid(
                    "invalid grouped radioactive energy",
                ));
            }
            closure.add(energy)?;
        }
        account.radiation_roundoff.add(closure.value()?.abs())?;
        account.radiation_roundoff.add(closure.export_roundoff()?)?;
        account
            .radiation_reconciled
            .add(emission.reconciled_spectral_energy_j)?;
        Ok(emission.emitted_energy_j)
    }

    pub(super) fn deposit(
        &self,
        emitted: &[(u64, Vec<f64>)],
        tolerance_j: f64,
        heat_j: &mut BTreeMap<CellId, Sum>,
        account: &mut Account,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<RadiationDeposition, FissionThermalError> {
        let borrowed: Vec<_> = emitted
            .iter()
            .map(|(id, values)| (*id, values.as_slice()))
            .collect();
        let mut checkpoint_error = None;
        let deposition = self.response.apply(&borrowed, tolerance_j, || {
            checkpoint().map_err(|error| {
                checkpoint_error = Some(error);
                RadiationError::Cancelled
            })
        });
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        radiation::credit_deposition(deposition?, heat_j, account, checkpoint)
    }
}

fn bind_decay(
    model: &Arc<crate::decay::DecayModel>,
    partition: Arc<PreparedDecayRadiation>,
) -> Result<Arc<PreparedDecayRadiation>, EndfError> {
    if !Arc::ptr_eq(model, partition.source().assembly().model()) {
        return Err(FissionError::InvalidInput(
            "radiation partition belongs to another decay model",
        )
        .into());
    }
    Ok(partition)
}

impl DecayRadiationHeatSource for DecayCellSource {
    type RadiationBinding = Arc<PreparedDecayRadiation>;
    fn bind_decay_radiation(
        &self,
        partition: Arc<PreparedDecayRadiation>,
        _: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Self::RadiationBinding, EndfError> {
        bind_decay(self.inventory.model(), partition)
    }
    fn decay_radiation_step(
        &self,
        binding: &Self::RadiationBinding,
        step: &Self::Step,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        binding.from_decay_step(self.inventory.model(), step, checkpoint)
    }
}

impl DecayRadiationHeatSource for DecayProductionCellSource {
    type RadiationBinding = Arc<PreparedDecayRadiation>;
    fn bind_decay_radiation(
        &self,
        partition: Arc<PreparedDecayRadiation>,
        _: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Self::RadiationBinding, EndfError> {
        bind_decay(self.inventory().model(), partition)
    }
    fn decay_radiation_step(
        &self,
        binding: &Self::RadiationBinding,
        step: &Self::Step,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        binding.from_decay_step(self.inventory().model(), &step.decay, checkpoint)
    }
}

impl DecayRadiationHeatSource for BurnupCellSource {
    type RadiationBinding = BurnupDecayRadiation;
    fn bind_decay_radiation(
        &self,
        partition: Arc<PreparedDecayRadiation>,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Self::RadiationBinding, EndfError> {
        partition.bind_burnup(Arc::clone(self.inventory.model()), checkpoint)
    }
    fn decay_radiation_step(
        &self,
        binding: &Self::RadiationBinding,
        step: &Self::Step,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        binding.from_burnup_step(self.inventory.model(), step, checkpoint)
    }
}

impl DecayRadiationHeatSource for BurnupProductionCellSource {
    type RadiationBinding = BurnupDecayRadiation;
    fn bind_decay_radiation(
        &self,
        partition: Arc<PreparedDecayRadiation>,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Self::RadiationBinding, EndfError> {
        partition.bind_burnup(Arc::clone(self.inventory().model()), checkpoint)
    }
    fn decay_radiation_step(
        &self,
        binding: &Self::RadiationBinding,
        step: &Self::Step,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        binding.from_burnup_step(self.inventory().model(), &step.burnup, checkpoint)
    }
}

impl DecayRadiationHeatSource for BurnupTransportCellSource {
    type RadiationBinding = BurnupDecayRadiation;
    fn bind_decay_radiation(
        &self,
        partition: Arc<PreparedDecayRadiation>,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Self::RadiationBinding, EndfError> {
        partition.bind_burnup(
            Arc::clone(self.state().inventories()[self.compartment_index()].model()),
            checkpoint,
        )
    }
    fn decay_radiation_step(
        &self,
        binding: &Self::RadiationBinding,
        step: &Self::Step,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        binding.from_burnup_step(
            self.state().inventories()[self.compartment_index()].model(),
            step.burnup(),
            checkpoint,
        )
    }
}

impl DecayRadiationHeatSource for CaptureActivationCellSource {
    type RadiationBinding = Arc<PreparedDecayRadiation>;

    fn bind_decay_radiation(
        &self,
        partition: Arc<PreparedDecayRadiation>,
        _: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Self::RadiationBinding, EndfError> {
        bind_decay(self.state.products().model(), partition)
    }

    fn decay_radiation_step(
        &self,
        binding: &Self::RadiationBinding,
        step: &Self::Step,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        binding.from_decay_step(
            self.state.products().model(),
            &step.products.decay,
            checkpoint,
        )
    }
}

impl DecayRadiationHeatSource for DiffusionCellSource {
    type RadiationBinding = BurnupDecayRadiation;

    fn bind_decay_radiation(
        &self,
        partition: Arc<PreparedDecayRadiation>,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Self::RadiationBinding, EndfError> {
        partition.bind_burnup(
            Arc::clone(self.state().inventories()[self.cell_index()].model()),
            checkpoint,
        )
    }

    fn decay_radiation_step(
        &self,
        binding: &Self::RadiationBinding,
        step: &Self::Step,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        binding.from_burnup_step(
            self.state().inventories()[self.cell_index()].model(),
            step.burnup(),
            checkpoint,
        )
    }
}

impl DecayRadiationHeatSource for FusionNeutronActivationCellSource {
    type RadiationBinding = BurnupDecayRadiation;

    fn bind_decay_radiation(
        &self,
        partition: Arc<PreparedDecayRadiation>,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Self::RadiationBinding, EndfError> {
        partition.bind_burnup(
            Arc::clone(self.state().inventories()[self.cell_index()].model()),
            checkpoint,
        )
    }

    fn decay_radiation_step(
        &self,
        binding: &Self::RadiationBinding,
        step: &Self::Step,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        binding.from_burnup_step(
            self.state().inventories()[self.cell_index()].model(),
            step.burnup(),
            checkpoint,
        )
    }
}

impl<S: FluxDrivenHeatSource + DecayRadiationHeatSource> DecayRadiationHeatSource
    for KineticCellSource<S>
{
    type RadiationBinding = S::RadiationBinding;

    fn bind_decay_radiation(
        &self,
        partition: Arc<PreparedDecayRadiation>,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Self::RadiationBinding, EndfError> {
        self.irradiation.bind_decay_radiation(partition, checkpoint)
    }

    fn decay_radiation_step(
        &self,
        binding: &Self::RadiationBinding,
        step: &Self::Step,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        self.irradiation
            .decay_radiation_step(binding, &step.irradiation, checkpoint)
    }
}

impl DecayRadiationHeatSource for StochasticDecayCellSource {
    type RadiationBinding = Arc<PreparedDecayRadiation>;

    fn bind_decay_radiation(
        &self,
        partition: Arc<PreparedDecayRadiation>,
        _: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Self::RadiationBinding, EndfError> {
        bind_decay(self.inventory.model(), partition)
    }

    fn decay_radiation_step(
        &self,
        binding: &Self::RadiationBinding,
        step: &Self::Step,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        binding.from_stochastic_step(self.inventory.model(), step, checkpoint)
    }
}

impl DecayRadiationHeatSource for StochasticDecayProductionCellSource {
    type RadiationBinding = Arc<PreparedDecayRadiation>;

    fn bind_decay_radiation(
        &self,
        partition: Arc<PreparedDecayRadiation>,
        _: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Self::RadiationBinding, EndfError> {
        bind_decay(self.production.inventory().model(), partition)
    }

    fn decay_radiation_step(
        &self,
        binding: &Self::RadiationBinding,
        step: &Self::Step,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        binding.from_stochastic_step(self.production.inventory().model(), &step.decay, checkpoint)
    }
}

impl<C: StochasticDecayCohort + 'static> DecayRadiationHeatSource
    for StochasticDecayCohortCellSource<C>
{
    type RadiationBinding = Arc<PreparedDecayRadiation>;

    fn bind_decay_radiation(
        &self,
        partition: Arc<PreparedDecayRadiation>,
        _: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Self::RadiationBinding, EndfError> {
        bind_decay(
            self.state().inventories()[self.compartment_index()].model(),
            partition,
        )
    }

    fn decay_radiation_step(
        &self,
        binding: &Self::RadiationBinding,
        step: &Self::Step,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        binding.from_stochastic_step(
            self.state().inventories()[self.compartment_index()].model(),
            step.decay(),
            checkpoint,
        )
    }
}
