//! Native temperature feedback for the existing prescribed photon source.
//!
//! Every source interval uses its original start-state caloric temperatures.
//! Mass, volume and layer geometry remain fixed; only explicitly supplied mass
//! energy coefficients change. The existing transport solve and optional
//! spatial refinement are reused. This is first-order frozen-coefficient time
//! coupling; reducing StepConfig::max_dt_s is a separate accuracy obligation.

use super::*;
use crate::radiation::transport::admit_slab_dimensions;
use crate::radiation::transport::material::temperature::TemperaturePhotonMassMaterial;

/// One original material cell in physical order along the fixed planar stack.
#[derive(Clone, Debug)]
pub struct TemperaturePhotonSlabLayer {
    /// Original native thermal cell; its stored mass and geometry stay fixed.
    pub cell_id: CellId,
    /// Original material identity and complete supported photon energy kernels.
    pub material: Arc<TemperaturePhotonMassMaterial>,
}

/// Storage and repeated-transport resources for native photon feedback.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonFeedbackBudget {
    /// Sum of retained table scalars over layer assignments, including repeated
    /// assignments to one shared table. This bounds work before interpolation.
    pub max_material_values: usize,
    /// Initial construction, or response rebuilds during one whole advance.
    /// Unchanged temperatures reuse the prepared response without a rebuild.
    pub max_transport_preparations: usize,
    /// Sum of the original transport owner's charged work over those rebuilds.
    /// An optional spatial-refinement solve includes all its candidate grids.
    pub max_transport_work: usize,
}

impl Default for PhotonFeedbackBudget {
    fn default() -> Self {
        Self {
            max_material_values: 1_048_576,
            max_transport_preparations: 1024,
            max_transport_work: 200_000_000,
        }
    }
}

/// Actual accepted response preparation in the last native advance.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PhotonFeedbackWork {
    /// Responses actually rebuilt; unchanged-temperature cache hits cost none.
    pub preparations: usize,
    /// Charged work from the original slab or spatial-refinement owner.
    pub transport_work: usize,
}

#[derive(Debug)]
struct FeedbackModel {
    data_id: String,
    geometry: PhotonSlabGeometry,
    layers: Arc<[TemperaturePhotonSlabLayer]>,
    sources: Arc<[PhotonSlabSource]>,
    ordinates: Arc<[PhotonOrdinate]>,
    limits: PhotonSlabLimits,
    refinement: Option<PhotonSlabRefinement>,
    budget: PhotonFeedbackBudget,
}

#[derive(Clone, Debug)]
pub(super) struct PhotonTemperatureFeedback {
    model: Arc<FeedbackModel>,
    temperatures_k: Arc<[(CellId, f64)]>,
    work: PhotonFeedbackWork,
    // Additional resources of one whole-trajectory refinement candidate. The
    // original model budget stays immutable, including on the returned owner.
    trial_work_limit: Option<PhotonFeedbackWork>,
}

#[derive(Clone, Debug)]
pub(super) struct PhotonFeedbackCandidate {
    pub(super) response: Arc<GroupedRadiationResponse>,
    pub(super) feedback: PhotonTemperatureFeedback,
}

impl ExternalPhotonHeating {
    /// Bind explicitly supplied temperature-dependent photon interactions to
    /// the original native heat cells and external illumination histories.
    ///
    /// Preparation derives density from original mass/volume and layer depth
    /// from volume/area, using the same binding as `from_slab`. Every changed
    /// temperature vector rebuilds the original energy-conserving multigroup
    /// response, optionally through its existing spatial refinement. Heat and
    /// downscattered energy change together; no attenuation-to-heat conversion
    /// or implicit density/phase law is introduced.
    ///
    /// Coefficients are held at the start of each explicit native interval.
    /// Refine the thermal timestep for this feedback; the transport's spatial
    /// agreement and accounting tolerance do not bound that temporal error.
    /// Endpoints and events must remain inside every original table's domain.
    #[allow(clippy::too_many_arguments)]
    pub fn from_temperature_slab(
        thermal: &DepositionSimulation<'_>,
        data_id: String,
        geometry: PhotonSlabGeometry,
        layers: Vec<TemperaturePhotonSlabLayer>,
        sources: Vec<PhotonIllumination>,
        ordinates: Vec<PhotonOrdinate>,
        limits: PhotonSlabLimits,
        refinement: Option<PhotonSlabRefinement>,
        budget: PhotonHeatingBudget,
        feedback_budget: PhotonFeedbackBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, FissionThermalError> {
        checkpoint()?;
        let groups = layers
            .first()
            .map_or(0, |layer| layer.material.group_edges_j().len() - 1);
        admit_slab_dimensions(layers.len(), groups, sources.len(), ordinates.len(), limits)?;
        if feedback_budget.max_transport_preparations == 0 {
            return Err(RadiationError::Budget("photon feedback preparation count").into());
        }
        let mut values = 0_usize;
        for layer in &layers {
            checkpoint()?;
            values = values
                .checked_add(layer.material.value_count())
                .filter(|&n| n <= feedback_budget.max_material_values)
                .ok_or(RadiationError::Budget("photon feedback material storage"))?;
        }
        let (geometries, powers) = split_sources(sources, budget, &mut checkpoint)?;
        admit_spectra(
            &powers,
            layers[0].material.radiation(),
            layers[0].material.group_edges_j(),
            &mut checkpoint,
        )?;
        let model = Arc::new(FeedbackModel {
            data_id,
            geometry,
            layers: layers.into(),
            sources: geometries.into(),
            ordinates: ordinates.into(),
            limits,
            refinement,
            budget: feedback_budget,
        });
        let temperatures = model.temperatures(thermal, &mut checkpoint)?;
        let (prepared, receivers) = model.prepare(
            thermal,
            &temperatures,
            feedback_budget.max_transport_work,
            &mut checkpoint,
        )?;
        let mut heating = Self::new(
            thermal,
            Arc::clone(&prepared.response),
            powers,
            receivers,
            budget,
            &mut checkpoint,
        )?;
        heating.feedback = Some(PhotonTemperatureFeedback {
            model,
            temperatures_k: temperatures.into(),
            work: PhotonFeedbackWork {
                preparations: 1,
                transport_work: prepared.work,
            },
            trial_work_limit: None,
        });
        checkpoint()?;
        Ok(heating)
    }

    /// Native temperatures at which `response()` was last prepared, in the
    /// original geometric order. They are source state, not inferred density.
    pub fn prepared_temperatures_k(&self) -> Option<&[(CellId, f64)]> {
        self.feedback
            .as_ref()
            .map(|feedback| feedback.temperatures_k.as_ref())
    }

    /// Construction work before the first advance; accepted work in the last
    /// whole advance afterward. Trial failures never change this owner.
    pub fn feedback_work(&self) -> Option<PhotonFeedbackWork> {
        self.feedback.as_ref().map(|feedback| feedback.work)
    }
}

impl FeedbackModel {
    fn temperatures(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<(CellId, f64)>, FissionThermalError> {
        let mut temperatures = Vec::with_capacity(self.layers.len());
        for layer in self.layers.iter() {
            checkpoint()?;
            let temperature = thermal.temperature_k(layer.cell_id)?;
            admit_temperature(layer, temperature)?;
            temperatures.push((layer.cell_id, temperature));
        }
        Ok(temperatures)
    }

    fn prepare(
        &self,
        thermal: &DepositionSimulation<'_>,
        temperatures: &[(CellId, f64)],
        remaining_work: usize,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(PhotonSlabResponse, Vec<RadiationReceiver>), FissionThermalError> {
        if remaining_work == 0 {
            return Err(RadiationError::Budget("photon feedback transport work").into());
        }
        let limits = PhotonSlabLimits {
            max_work: self.limits.max_work.min(remaining_work),
            ..self.limits
        };
        // Resource refusal precedes material sampling and native grid expansion.
        admit_slab_dimensions(
            self.layers.len(),
            self.layers[0].material.group_edges_j().len() - 1,
            self.sources.len(),
            self.ordinates.len(),
            limits,
        )?;
        let layers = photon_call(checkpoint, |poll| {
            self.layers
                .iter()
                .zip(temperatures)
                .map(|(layer, &(_, temperature))| {
                    Ok(PhotonSlabLayer {
                        cell_id: layer.cell_id,
                        material: Arc::new(layer.material.sample(temperature, &mut *poll)?),
                    })
                })
                .collect::<Result<Vec<_>, RadiationError>>()
        })?;
        let input = photon_material_input(
            thermal,
            self.geometry,
            layers,
            self.sources.iter().map(|source| source.source_id).collect(),
            false,
            self.ordinates.len(),
            limits,
            checkpoint,
        )?;
        let transport = photon_call(checkpoint, |poll| {
            if let Some(refinement) = self.refinement {
                Ok(refine_slab_response_with_sources(
                    self.data_id.clone(),
                    input.radiation,
                    input.edges_j,
                    input.cells,
                    self.sources.to_vec(),
                    self.ordinates.to_vec(),
                    limits,
                    refinement,
                    poll,
                )?
                .transport)
            } else {
                prepare_slab_response_with_sources(
                    self.data_id.clone(),
                    input.radiation,
                    input.edges_j,
                    input.cells,
                    self.sources.to_vec(),
                    self.ordinates.to_vec(),
                    limits,
                    poll,
                )
            }
        })?;
        Ok((transport, input.receivers))
    }
}

impl PhotonTemperatureFeedback {
    pub(super) fn limit_trial_work(&mut self, limit: Option<PhotonFeedbackWork>) {
        self.trial_work_limit = limit;
    }

    pub(super) fn begin(&mut self) {
        self.work = PhotonFeedbackWork::default();
    }

    // Support is checked at every native endpoint/event, even when that instant
    // is the final requested state and there will be no later source query.
    pub(super) fn validate_temperature_support(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), FissionThermalError> {
        for layer in self.model.layers.iter() {
            checkpoint()?;
            admit_temperature(layer, thermal.temperature_k(layer.cell_id)?)?;
        }
        Ok(())
    }

    pub(super) fn candidate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<PhotonFeedbackCandidate>, FissionThermalError> {
        let mut changed = false;
        for &(id, temperature) in self.temperatures_k.iter() {
            checkpoint()?;
            changed |= thermal.temperature_k(id)? != temperature;
        }
        if !changed {
            return Ok(None);
        }
        let maximum_preparations =
            self.trial_work_limit
                .map_or(self.model.budget.max_transport_preparations, |limit| {
                    limit
                        .preparations
                        .min(self.model.budget.max_transport_preparations)
                });
        let maximum_work =
            self.trial_work_limit
                .map_or(self.model.budget.max_transport_work, |limit| {
                    limit
                        .transport_work
                        .min(self.model.budget.max_transport_work)
                });
        if self.work.preparations >= maximum_preparations {
            return Err(RadiationError::Budget("photon feedback preparation count").into());
        }
        let remaining = maximum_work
            .checked_sub(self.work.transport_work)
            .ok_or(RadiationError::Budget("photon feedback transport work"))?;
        // Keep an arbitrary original thermal callback error exact through the
        // existing material/transport API's FissionError cancellation adapter.
        let mut checkpoint_error = None;
        let result = (|| {
            let mut poll = || {
                checkpoint().map_err(|error| {
                    checkpoint_error = Some(error);
                    FissionError::Cancelled
                })
            };
            let temperatures = self.model.temperatures(thermal, &mut poll)?;
            let (prepared, _) = self
                .model
                .prepare(thermal, &temperatures, remaining, &mut poll)?;
            let mut feedback = self.clone();
            feedback.temperatures_k = temperatures.into();
            feedback.work.preparations = feedback
                .work
                .preparations
                .checked_add(1)
                .ok_or(FissionThermalError::Arithmetic)?;
            feedback.work.transport_work = feedback
                .work
                .transport_work
                .checked_add(prepared.work)
                .filter(|&work| work <= maximum_work)
                .ok_or(RadiationError::Budget("photon feedback transport work"))?;
            Ok(Some(PhotonFeedbackCandidate {
                response: prepared.response,
                feedback,
            }))
        })();
        if let Some(error) = checkpoint_error {
            Err(error.into())
        } else {
            result
        }
    }
}

fn admit_temperature(
    layer: &TemperaturePhotonSlabLayer,
    temperature_k: f64,
) -> Result<(), FissionThermalError> {
    let (minimum_k, maximum_k) = layer.material.temperature_bounds_k();
    if !temperature_k.is_finite() || temperature_k < minimum_k || temperature_k > maximum_k {
        Err(FissionThermalError::PhotonTemperatureDomain {
            cell_id: layer.cell_id,
            temperature_k,
            minimum_k,
            maximum_k,
        })
    } else {
        Ok(())
    }
}
