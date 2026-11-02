//! Bind explicit optical constitutive laws to the original native state.
//!
//! Immutable laws and resolved packets are shared across candidate clones.
//! A distinct receiver/law pair resolves once per source interval; only changed
//! optical constants invalidate Fresnel packets, and darkness needs no solve.

use super::feedback::SurfaceConductorLaw;
use super::*;
use crate::conductor::ConductorOptics;
use std::collections::BTreeSet;

/// Prescribed geometry and polarization for one intercepted-power source.
/// Wavelengths are vacuum nm; indices and s fractions are dimensionless.
/// Area and cosine projection must already be included in supplied watts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceIncidence<const N: usize> {
    pub wavelengths_nm: [f64; N],
    pub incident_eta: [f64; N],
    pub incidence_cosine: f64,
    pub s_power_fraction: [f64; N],
}

impl<const N: usize> SurfaceIncidence<N> {
    fn response(
        &self,
        optics: &ConductorOptics,
    ) -> Result<SurfacePowerResponse<N>, SurfaceHeatingError> {
        SurfacePowerResponse::from_conductor(
            self.wavelengths_nm,
            optics,
            self.incident_eta,
            self.incidence_cosine,
            self.s_power_fraction,
        )
    }

    fn hash(&self, hash: &mut DomainHasher) {
        hash.update(&self.incidence_cosine.to_bits().to_le_bytes());
        for values in [
            &self.wavelengths_nm,
            &self.incident_eta,
            &self.s_power_fraction,
        ] {
            for value in values {
                hash.update(&value.to_bits().to_le_bytes());
            }
        }
    }
}

/// An incident source whose substrate n/k follow an explicitly authored law.
/// The original native caloric owner supplies T and h; no phase mixing or
/// material-card coordinate mapping is inferred. The additional temperature
/// bounds intersect the law's own support, including while illumination is off.
#[derive(Clone, Debug)]
pub struct StatefulSurfaceIlluminationSource<const N: usize> {
    pub source_id: u64,
    pub thermal_cell: CellId,
    pub material_id: MaterialId,
    pub temperature_bounds_k: (f64, f64),
    pub incidence: SurfaceIncidence<N>,
    pub optics: Arc<SurfaceConductorLaw>,
    pub power: SpectralPowerSchedule<N>,
}

#[derive(Debug)]
struct BoundLaw {
    law: Arc<SurfaceConductorLaw>,
    cell: CellId,
    material: MaterialId,
    source_indices: Vec<usize>,
}

#[derive(Clone, Copy, Debug)]
struct BoundResponse<const N: usize> {
    law: usize,
    incidence: SurfaceIncidence<N>,
}

#[derive(Debug)]
struct FeedbackModel<const N: usize> {
    laws: Vec<BoundLaw>,
    sources: Vec<Option<BoundResponse<N>>>,
}

#[derive(Clone, Debug)]
struct FeedbackCache<const N: usize> {
    optics: Vec<ConductorOptics>,
    responses: Vec<Option<SurfacePowerResponse<N>>>,
}

#[derive(Clone, Debug)]
pub(super) struct SurfaceFeedbackState<const N: usize> {
    model: Arc<FeedbackModel<N>>,
    cache: Arc<FeedbackCache<N>>,
}

impl<const N: usize> SurfaceOpticalHeating<N> {
    /// Admit fixed and state-dependent sources into the SAME surface heat owner.
    ///
    /// `max_total_optical_knots` bounds state-coordinate rows across distinct
    /// law identities, counting each shared law once. Each row has the original
    /// fixed nine wavelength knots. The existing source/segment/evaluation
    /// budgets cover both source kinds and bound all per-interval optical work.
    ///
    /// Optical constants are evaluated at the accepted interval START; native
    /// thermal advancement supplies its end state. Time refinement is required
    /// for nonlinear feedback. Exact pulse fluence does not eliminate this
    /// temporal approximation or establish an endpoint-implicit optical solve.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_stateful(
        fixed: &[SurfaceIlluminationSource<N>],
        stateful: &[StatefulSurfaceIlluminationSource<N>],
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
        bindings
            .try_reserve_exact(stateful.len())
            .map_err(|_| SurfaceHeatingError::Budget)?;
        let mut optics = Vec::new();
        optics
            .try_reserve_exact(stateful.len())
            .map_err(|_| SurfaceHeatingError::Budget)?;
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
                    optics.push(source.optics.resolve(
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
                        incidence: source.incidence,
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
                response: source.incidence.response(&optics[slot])?,
                power: source.power.clone(),
            });
        }
        // The existing admission owns canonical source order, receiver domains,
        // wavelength consistency, schedules and all original fixed-source rules.
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
            DomainHasher::new("org.frankensim.render.stateful-surface-optical-heating.v1");
        hash.update(&(N as u64).to_le_bytes());
        hash.update(&(count as u64).to_le_bytes());
        for (index, source) in result.model.sources.iter().enumerate() {
            checkpoint()?;
            let receiver = source
                .receiver
                .expect("single-surface source retains its receiver");
            for id in [
                source.source_id,
                receiver.thermal_cell.0,
                receiver.material_id.0,
            ] {
                hash.update(&id.to_le_bytes());
            }
            for t in [
                receiver.temperature_bounds_k.0,
                receiver.temperature_bounds_k.1,
            ] {
                hash.update(&t.to_bits().to_le_bytes());
            }
            hash.update(source.power.identity().as_bytes());
            let bound = by_source.get(&source.source_id).copied();
            match bound {
                Some(bound) => {
                    hash.update(&[1]);
                    hash.update(bindings[bound.law].law.identity().as_bytes());
                    bound.incidence.hash(&mut hash);
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
        // Identity describes the constitutive program, not its admission state.
        // Two starts within the same supported law retain the same model identity.
        Arc::get_mut(&mut result.model)
            .expect("new surface model is uniquely owned")
            .identity = hash.finalize();
        result.feedback = Some(SurfaceFeedbackState {
            model: Arc::new(FeedbackModel {
                laws: bindings,
                sources,
            }),
            cache: Arc::new(FeedbackCache { optics, responses }),
        });
        checkpoint()?;
        Ok(result)
    }
}

impl<const N: usize> SurfaceFeedbackState<N> {
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
            let optics = bound.law.optics_at_coordinate(coordinate)?;
            if optics != self.cache.optics[index] {
                let cache = Arc::make_mut(&mut self.cache);
                cache.optics[index] = optics;
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
        // Domain validation still ran, but zero intercepted watts cannot heat.
        // Retain the invalidated packet until nonzero light needs its response.
        if !illuminated {
            return Ok(None);
        }
        checkpoint()?;
        *evaluations = evaluations
            .checked_add(1)
            .ok_or(SurfaceHeatingError::Budget)?;
        let response = bound.incidence.response(&self.cache.optics[bound.law])?;
        Arc::make_mut(&mut self.cache).responses[source_index] = Some(response);
        Ok(Some(response))
    }
}
