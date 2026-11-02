//! Bind state-dependent diffuse emissivity to the existing surface heat owner.

use super::diffuse_feedback::SurfaceDiffuseLaw;
use super::*;
use crate::thermal_surface::SpectralEmissivity;
use std::collections::BTreeSet;

#[derive(Clone, Debug)]
pub struct StatefulDiffuseSurfaceIlluminationSource<const N: usize> {
    pub source_id: u64,
    pub thermal_cell: CellId,
    pub material_id: MaterialId,
    pub temperature_bounds_k: (f64, f64),
    pub wavelengths_nm: [f64; N],
    pub optics: Arc<SurfaceDiffuseLaw>,
    pub power: SpectralPowerSchedule<N>,
}

#[derive(Debug)]
struct BoundLaw {
    law: Arc<SurfaceDiffuseLaw>,
    cell: CellId,
    material: MaterialId,
    source_indices: Vec<usize>,
}

#[derive(Clone, Copy, Debug)]
struct BoundResponse<const N: usize> {
    law: usize,
    wavelengths_nm: [f64; N],
}

#[derive(Debug)]
struct FeedbackModel<const N: usize> {
    laws: Vec<BoundLaw>,
    sources: Vec<Option<BoundResponse<N>>>,
}

#[derive(Clone, Debug)]
struct FeedbackCache<const N: usize> {
    emissivities: Vec<SpectralEmissivity>,
    responses: Vec<Option<SurfacePowerResponse<N>>>,
}

#[derive(Clone, Debug)]
pub(super) struct SurfaceDiffuseFeedbackState<const N: usize> {
    model: Arc<FeedbackModel<N>>,
    cache: Arc<FeedbackCache<N>>,
}

impl<const N: usize> SurfaceOpticalHeating<N> {
    /// Admit fixed and state-dependent diffuse sources into one native heat owner.
    ///
    /// The state law supplies spectral emissivity and therefore reciprocal
    /// absorptivity. Optical state is evaluated at each accepted interval start,
    /// exactly like conductor feedback; nonlinear time accuracy still requires
    /// refinement. Constant rows reuse their compiled response bit-for-bit.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_stateful_diffuse(
        fixed: &[SurfaceIlluminationSource<N>],
        stateful: &[StatefulDiffuseSurfaceIlluminationSource<N>],
        thermal: &DepositionSimulation<'_>,
        budget: SurfaceHeatingBudget,
        max_total_optical_knots: usize,
        tolerance: SurfaceHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SurfaceHeatingError> {
        checkpoint()?;
        let count = fixed
            .len()
            .checked_add(stateful.len())
            .ok_or(SurfaceHeatingError::Budget)?;
        if count > budget.max_sources {
            return Err(SurfaceHeatingError::Budget);
        }
        if stateful.is_empty() {
            return Self::try_new(fixed, thermal, budget, tolerance, checkpoint);
        }
        let mut knots = 0usize;
        let mut law_ids = BTreeSet::new();
        let mut segments = 0usize;
        for source in fixed {
            checkpoint()?;
            segments = segments
                .checked_add(source.power.segment_count())
                .ok_or(SurfaceHeatingError::Budget)?;
        }
        for source in stateful {
            checkpoint()?;
            if law_ids.insert(*source.optics.identity().as_bytes()) {
                knots = knots
                    .checked_add(source.optics.row_count())
                    .ok_or(SurfaceHeatingError::Budget)?;
                if knots > max_total_optical_knots {
                    return Err(SurfaceHeatingError::Budget);
                }
            }
            segments = segments
                .checked_add(source.power.segment_count())
                .ok_or(SurfaceHeatingError::Budget)?;
        }
        if segments > budget.max_total_power_segments {
            return Err(SurfaceHeatingError::Budget);
        }

        let mut prepared = Vec::new();
        prepared
            .try_reserve_exact(count)
            .map_err(|_| SurfaceHeatingError::Budget)?;
        prepared.extend_from_slice(fixed);
        let mut bindings = Vec::<BoundLaw>::new();
        let mut emissivities = Vec::new();
        let mut slots = BTreeMap::new();
        let mut by_source = BTreeMap::new();
        for source in stateful {
            checkpoint()?;
            let key = (
                source.thermal_cell,
                source.material_id,
                *source.optics.identity().as_bytes(),
            );
            let slot = match slots.get(&key) {
                Some(&slot) => slot,
                None => {
                    let slot = bindings.len();
                    emissivities.push(source.optics.resolve(
                        thermal,
                        source.thermal_cell,
                        source.material_id,
                    )?);
                    bindings.push(BoundLaw {
                        law: Arc::clone(&source.optics),
                        cell: source.thermal_cell,
                        material: source.material_id,
                        source_indices: Vec::new(),
                    });
                    slots.insert(key, slot);
                    slot
                }
            };
            if by_source
                .insert(
                    source.source_id,
                    BoundResponse {
                        law: slot,
                        wavelengths_nm: source.wavelengths_nm,
                    },
                )
                .is_some()
            {
                return Err(invalid("surface illumination source IDs must be unique"));
            }
            prepared.push(SurfaceIlluminationSource {
                source_id: source.source_id,
                thermal_cell: source.thermal_cell,
                material_id: source.material_id,
                temperature_bounds_k: source.temperature_bounds_k,
                response: SurfacePowerResponse::from_spectral_diffuse(
                    source.wavelengths_nm,
                    &emissivities[slot],
                    source.optics.identity(),
                )?,
                power: source.power.clone(),
            });
        }

        let mut result = Self::try_new(&prepared, thermal, budget, tolerance, &mut checkpoint)?;
        let mut sources = Vec::new();
        let mut responses = Vec::new();
        sources
            .try_reserve_exact(count)
            .map_err(|_| SurfaceHeatingError::Budget)?;
        responses
            .try_reserve_exact(count)
            .map_err(|_| SurfaceHeatingError::Budget)?;
        let mut hash =
            DomainHasher::new("org.frankensim.render.stateful-diffuse-surface-heating.v1");
        hash.update(&(N as u64).to_le_bytes());
        hash.update(&(count as u64).to_le_bytes());
        for (index, source) in result.model.sources.iter().enumerate() {
            checkpoint()?;
            for id in [
                source.source_id,
                source.thermal_cell.0,
                source.material_id.0,
            ] {
                hash.update(&id.to_le_bytes());
            }
            for t in [source.temperature_bounds_k.0, source.temperature_bounds_k.1] {
                hash.update(&t.to_bits().to_le_bytes());
            }
            hash.update(source.power.identity().as_bytes());
            let bound = by_source.get(&source.source_id).copied();
            match bound {
                Some(bound) => {
                    hash.update(&[1]);
                    hash.update(bindings[bound.law].law.identity().as_bytes());
                    for wavelength in bound.wavelengths_nm {
                        hash.update(&wavelength.to_bits().to_le_bytes());
                    }
                    bindings[bound.law].source_indices.push(index);
                    responses.push(Some(source.response));
                }
                None => {
                    hash.update(&[0]);
                    hash.update(source.response.identity().as_bytes());
                    responses.push(None);
                }
            }
            sources.push(bound);
        }
        Arc::get_mut(&mut result.model)
            .expect("new surface model is uniquely owned")
            .identity = hash.finalize();
        result.diffuse_feedback = Some(SurfaceDiffuseFeedbackState {
            model: Arc::new(FeedbackModel {
                laws: bindings,
                sources,
            }),
            cache: Arc::new(FeedbackCache {
                emissivities,
                responses,
            }),
        });
        checkpoint()?;
        Ok(result)
    }
}

impl<const N: usize> SurfaceDiffuseFeedbackState<N> {
    pub(super) fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), SurfaceHeatingError> {
        for bound in &self.model.laws {
            checkpoint()?;
            bound
                .law
                .coordinate_at(thermal, bound.cell, bound.material)?;
        }
        Ok(())
    }

    pub(super) fn refresh(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), SurfaceHeatingError> {
        for (index, bound) in self.model.laws.iter().enumerate() {
            checkpoint()?;
            let coordinate = bound
                .law
                .coordinate_at(thermal, bound.cell, bound.material)?;
            let emissivity = bound.law.emissivity_at_coordinate(coordinate)?;
            if emissivity != self.cache.emissivities[index] {
                let cache = Arc::make_mut(&mut self.cache);
                cache.emissivities[index] = emissivity;
                for &source in &bound.source_indices {
                    cache.responses[source] = None;
                }
            }
        }
        Ok(())
    }

    pub(super) fn response(
        &mut self,
        source_index: usize,
        illuminated: bool,
        evaluations: &mut usize,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<SurfacePowerResponse<N>>, SurfaceHeatingError> {
        let Some(bound) = self.model.sources[source_index] else {
            return Ok(None);
        };
        if let Some(response) = self.cache.responses[source_index] {
            return Ok(Some(response));
        }
        if !illuminated {
            return Ok(None);
        }
        checkpoint()?;
        *evaluations = evaluations
            .checked_add(1)
            .ok_or(SurfaceHeatingError::Budget)?;
        let response = SurfacePowerResponse::from_spectral_diffuse(
            bound.wavelengths_nm,
            &self.cache.emissivities[bound.law],
            self.model.laws[bound.law].law.identity(),
        )?;
        Arc::make_mut(&mut self.cache).responses[source_index] = Some(response);
        Ok(Some(response))
    }
}
