//! Finite-inventory emitted-photon feedback on the accepted thermal clock.
//!
//! Complete ballistic rays or elastic walks are sampled once from the original cells. Their
//! geometry is reused while molecular absorption follows the same evolving
//! ground inventory as all external pumps. Spatial/angular sampling, kinetic
//! refinement and thermal feedback are distinct approximations and budgets.
//! Photon flight is quasistatic. Elastic walks use the original scattering
//! spectra and phase laws; interfaces and coherent light are not inferred.

use super::cross_sections::{CrossSections, EmitterAbsorptionCrossSections};
use super::{
    CheckedSum, ContentHash, DepositionError, DepositionSimulation, DomainHasher,
    LuminescenceSource, OpticalHeatingError, SaturableEmitterCell, SaturableLuminescentHeating,
    SingletTripletSource, finite, invalid, preserve_checkpoint, resolve_model,
};
use crate::volumes::transport::luminescence::LuminescenceGrid;
use crate::volumes::transport::scattering::ScatteringGrid;
use crate::volumes::transport::{AbsorbingGrid, VolumeTransportError};
use fs_exec::Cx;
use fs_thermal::luminescence::transfer::pump_beams::recycling::{
    EmissionTerminal, SaturableEmissionPath, SaturablePhotonRecycling,
};
use fs_thermal::luminescence::transfer::pump_beams::{SaturablePumpPath, SaturablePumpSegment};
use std::{collections::BTreeMap, sync::Arc};

pub use fs_thermal::luminescence::transfer::pump_beams::recycling::SaturableRecyclingLimits;

mod evolution;
mod geometry;
pub(super) use evolution::step_sources;

/// Source-bound molecular absorption at the common emission nodes. Both
/// possible destination states compete for the same emitted photon.
pub type RecyclingEmitterCrossSections<const M: usize> = EmitterAbsorptionCrossSections<M>;

/// One reusable spatial/angular quadrature over complete open-grid rays.
/// Increasing ray count or changing the seed assesses sampling sensitivity;
/// the nonlinear result is not an unbiased estimator or a confidence bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonRecyclingGeometry {
    pub histories_per_emitter: u32,
    /// Includes external paths and every emitted spectral-pair/spin path.
    /// Scattered walks use one physical wavelength per path and two spin copies.
    pub maximum_paths: usize,
    /// Includes every retained segment over all those paths, not per emitter.
    pub maximum_total_segments: usize,
    pub seed: u64,
    pub source_id: u64,
    /// Explicit validity domain of the immutable passive host optical snapshot,
    /// including its scattering coefficients and phase laws when present.
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    /// Arithmetic allowance for complete emission-probability rows.
    pub fraction_roundoff: f64,
}

#[derive(Clone, Debug)]
struct RayCell<const M: usize> {
    cell: usize,
    length_m: f64,
    volume_m3: f64,
    passive_depth: [f64; M],
}

#[derive(Clone, Debug)]
struct EmissionRay<const M: usize> {
    source_cell: usize,
    /// Ballistic rays share all bins; an elastic walk samples one wavelength.
    lane: Option<usize>,
    cells: Vec<RayCell<M>>,
}

#[derive(Debug)]
struct Geometry<const M: usize> {
    config: PhotonRecyclingGeometry,
    cross_sections: Vec<RecyclingEmitterCrossSections<M>>,
    volumes_m3: Vec<f64>,
    rays: Vec<EmissionRay<M>>,
    limits: SaturableRecyclingLimits,
}

#[derive(Debug)]
pub(super) struct PreparedRecycling {
    identity: ContentHash,
    temperatures: Vec<u64>,
    pub(super) solver: SaturablePhotonRecycling<2>,
    /// Original common emission-bin indices, in exactly the solver path order.
    pub(super) path_bins: Vec<[Option<usize>; 2]>,
}

#[derive(Clone, Debug)]
pub(super) struct PhotonRecycling<const M: usize> {
    geometry: Arc<Geometry<M>>,
    prepared: Option<Arc<PreparedRecycling>>,
}

impl<const M: usize> SaturableLuminescentHeating<M> {
    /// Construct the finite owner and its complete emitted-photon transport
    /// together. This admits passive host absorption without first building an
    /// unused fixed-probability operator. Every positive-density host cell needs its
    /// original thermal binding; a zero emitter inventory is permitted.
    /// Maintained prescribed excitation and external depleting beams address
    /// the same inventories; beams may be bound before or after construction.
    /// Dark precharged stores also work without an external pump. Existing pair
    /// laws remain independently bound.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_photon_recycling(
        cx: &Cx<'_>,
        medium: &AbsorbingGrid<'_>,
        singlet: LuminescenceSource<M>,
        triplet: SingletTripletSource<M>,
        bindings: &[SaturableEmitterCell],
        thermal: &DepositionSimulation<'_>,
        max_grid_cells: usize,
        energy_tolerance_j: f64,
        cross_sections: &[RecyclingEmitterCrossSections<M>],
        geometry: PhotonRecyclingGeometry,
        limits: SaturableRecyclingLimits,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        Self::try_new_impl(
            medium,
            singlet,
            triplet,
            bindings,
            thermal,
            max_grid_cells,
            energy_tolerance_j,
            None,
            true,
            &mut checkpoint,
        )?
        .with_photon_recycling(
            cx,
            medium,
            cross_sections,
            thermal,
            geometry,
            limits,
            checkpoint,
        )
    }

    /// Bind active photon return to every original finite emitter inventory.
    /// The host medium excludes the molecular cross sections declared here.
    /// Existing passive transport is superseded by these complete rays, while
    /// its explicit thermal-domain admission remains checked if it was bound.
    /// No excitation probability or missing coefficient is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn with_photon_recycling(
        self,
        cx: &Cx<'_>,
        medium: &AbsorbingGrid<'_>,
        cross_sections: &[RecyclingEmitterCrossSections<M>],
        thermal: &DepositionSimulation<'_>,
        config: PhotonRecyclingGeometry,
        limits: SaturableRecyclingLimits,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        self.with_photon_recycling_impl(
            cx, medium, None, cross_sections, thermal, config, limits, checkpoint,
        )
    }

    /// Construct finite populations with complete, wavelength-resolved elastic
    /// walks through the original scattering medium. Passive host absorption
    /// excludes the molecular cross sections supplied here. Every occupied host
    /// cell needs its original thermal binding; zero emitter counts are allowed.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_scattered_photon_recycling(
        cx: &Cx<'_>,
        medium: &ScatteringGrid<'_>,
        singlet: LuminescenceSource<M>,
        triplet: SingletTripletSource<M>,
        bindings: &[SaturableEmitterCell],
        thermal: &DepositionSimulation<'_>,
        max_grid_cells: usize,
        energy_tolerance_j: f64,
        cross_sections: &[RecyclingEmitterCrossSections<M>],
        geometry: PhotonRecyclingGeometry,
        maximum_scattering_events_per_history: u32,
        limits: SaturableRecyclingLimits,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        Self::try_new_impl(
            medium.absorption_grid(),
            singlet,
            triplet,
            bindings,
            thermal,
            max_grid_cells,
            energy_tolerance_j,
            None,
            true,
            &mut checkpoint,
        )?
        .with_scattered_photon_recycling(
            cx,
            medium,
            cross_sections,
            thermal,
            geometry,
            maximum_scattering_events_per_history,
            limits,
            checkpoint,
        )
    }

    /// Bind original elastic scattering to the same finite ground inventory.
    /// Each emission bin samples independent complete walks with its own
    /// original phase law; no spectral mixture weights or fixed molecular
    /// absorption probabilities are introduced. Absorption is integrated along
    /// every ordered visit, including returns to a previously visited cell.
    ///
    /// A collision or traversal cap refuses the entire binding without
    /// classifying unfinished photons as escape or heat. Zero scattering takes
    /// the exact ballistic route, preserving its samples and source identity.
    #[allow(clippy::too_many_arguments)]
    pub fn with_scattered_photon_recycling(
        self,
        cx: &Cx<'_>,
        medium: &ScatteringGrid<'_>,
        cross_sections: &[RecyclingEmitterCrossSections<M>],
        thermal: &DepositionSimulation<'_>,
        config: PhotonRecyclingGeometry,
        maximum_scattering_events_per_history: u32,
        limits: SaturableRecyclingLimits,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        let (lo, hi) = medium.wavelength_bounds_nm();
        if self
            .singlet
            .fluorescence()
            .wavelengths_nm()
            .iter()
            .any(|wavelength| !(lo..=hi).contains(wavelength))
        {
            return Err(invalid("emission exceeds original scattering snapshot support"));
        }
        if !medium.has_scattering() {
            return self.with_photon_recycling(
                cx,
                medium.absorption_grid(),
                cross_sections,
                thermal,
                config,
                limits,
                checkpoint,
            );
        }
        self.with_photon_recycling_impl(
            cx,
            medium.absorption_grid(),
            Some((medium, maximum_scattering_events_per_history)),
            cross_sections,
            thermal,
            config,
            limits,
            checkpoint,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn with_photon_recycling_impl(
        mut self,
        cx: &Cx<'_>,
        medium: &AbsorbingGrid<'_>,
        scattering: Option<(&ScatteringGrid<'_>, u32)>,
        cross_sections: &[RecyclingEmitterCrossSections<M>],
        thermal: &DepositionSimulation<'_>,
        config: PhotonRecyclingGeometry,
        limits: SaturableRecyclingLimits,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalHeatingError> {
        self.check_state(thermal, true, &mut checkpoint)?;
        if !(1..=16).contains(&M)
            || self.recycling.is_some()
            || self.pending.is_some()
            || self.energy.duration_s != 0.0
            || self.emitters.is_empty()
            || self.medium_identity != medium.content_identity()
        {
            return Err(invalid(
                "finite photon recycling requires original unadvanced inventories and matching medium",
            ));
        }
        let mut hash = DomainHasher::new(if scattering.is_some() {
            "org.frankensim.render.finite-scattered-photon-recycling.v1"
        } else {
            "org.frankensim.render.finite-photon-recycling.v1"
        });
        hash.update(self.source_identity.as_bytes());
        if let Some((medium, maximum_events)) = scattering {
            hash.update(medium.content_identity().as_bytes());
            hash.update(&maximum_events.to_le_bytes());
        }
        let geometry = geometry::compile(
            &self,
            cx,
            medium,
            scattering,
            cross_sections,
            config,
            limits,
            &mut hash,
            &mut checkpoint,
        )?;
        let feedback = PhotonRecycling {
            geometry: Arc::new(geometry),
            prepared: None,
        };
        // Include the complete geometry/source declaration before preparing its
        // cache so subsequent beam/pulse bindings naturally invalidate it.
        self.source_identity = hash.finalize();
        let prepared = feedback.prepare(&self, thermal, &mut checkpoint)?;
        checkpoint()?;
        self.recycling = Some(PhotonRecycling {
            prepared: Some(prepared),
            ..feedback
        });
        Ok(self)
    }

    pub fn has_photon_recycling(&self) -> bool {
        self.recycling.is_some()
    }

    /// Observe accepted endpoint emission with its accepted finite-population
    /// molecular opacity and original passive host. The camera's exact cell
    /// integral remains distinct from the sampled recycling quadrature.
    /// Interval-mean powers cannot generally be combined with endpoint opacity
    /// to obtain a changing-opacity exposure; resolve exposure times explicitly.
    pub fn recycling_emission_grid<'m, 'g>(
        &self,
        cx: &Cx<'_>,
        medium: &'m AbsorbingGrid<'g>,
        max_cells: usize,
    ) -> Result<LuminescenceGrid<'m, 'g, M>, OpticalHeatingError> {
        if medium.content_identity() != self.medium_identity {
            return Err(invalid(
                "recycling observation differs from original optical medium",
            ));
        }
        let recycling = self
            .recycling
            .as_ref()
            .ok_or_else(|| invalid("active absorption observation requires photon recycling"))?;
        let prepared = recycling
            .prepared
            .as_ref()
            .ok_or_else(|| invalid("accepted molecular absorption is not prepared"))?;
        let wavelengths = *self.singlet.fluorescence().wavelengths_nm();
        let mut absorption = Vec::with_capacity(self.emitters.len());
        for (index, emitter) in self.emitters.iter().enumerate() {
            cx.checkpoint()
                .map_err(|_| OpticalHeatingError::Thermal(DepositionError::Cancelled))?;
            let sections = recycling.geometry.cross_sections[index]
                .evaluate(wavelengths, f64::from_bits(prepared.temperatures[index]))?;
            let mut values = [0.0; M];
            for lane in 0..M {
                let mut total = CheckedSum::default();
                total.add(sections[lane][0])?;
                total.add(sections[lane][1])?;
                values[lane] = super::super::rate_power(
                    total.value()?,
                    emitter.populations[0],
                    recycling.geometry.volumes_m3[index],
                )?;
                if total.value()? > 0.0 && emitter.populations[0] > 0.0 && values[lane] == 0.0 {
                    return Err(invalid(
                        "positive accepted molecular absorption is below representation",
                    ));
                }
            }
            absorption.push((emitter.binding.grid_cell, values));
        }
        LuminescenceGrid::try_new(
            cx,
            medium,
            self.source_identity,
            wavelengths,
            &self.emission_w,
            max_cells,
        )
        .and_then(|field| {
            field.with_additional_absorption(cx, self.source_identity, &absorption, max_cells)
        })
        .map_err(transport_error)
    }
}

impl<const M: usize> PhotonRecycling<M> {
    pub(super) fn prepare(
        &self,
        owner: &SaturableLuminescentHeating<M>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Arc<PreparedRecycling>, OpticalHeatingError> {
        let domain = self.geometry.config;
        let mut temperatures = Vec::with_capacity(owner.emitters.len());
        for emitter in &owner.emitters {
            checkpoint()?;
            let temperature = thermal.temperature_k(emitter.binding.thermal_cell)?;
            if !temperature.is_finite()
                || !(domain.minimum_temperature_k..=domain.maximum_temperature_k)
                    .contains(&temperature)
            {
                return Err(invalid(
                    "native temperature exceeds passive recycling snapshot support",
                ));
            }
            temperatures.push(temperature.to_bits());
        }
        if let Some(prepared) = &self.prepared
            && prepared.identity == owner.source_identity
            && prepared.temperatures == temperatures
        {
            return Ok(Arc::clone(prepared));
        }
        let wavelengths = *owner.singlet.fluorescence().wavelengths_nm();
        let mut models = Vec::with_capacity(owner.emitters.len());
        let mut sections = Vec::with_capacity(owner.emitters.len());
        let mut cache = None;
        for (index, &temperature) in temperatures.iter().enumerate() {
            checkpoint()?;
            let temperature = f64::from_bits(temperature);
            models.push(resolve_model(
                &owner.singlet,
                &owner.triplet,
                temperature,
                &mut cache,
            )?);
            let values = self.geometry.cross_sections[index].evaluate(wavelengths, temperature)?;
            for lane in 0..M {
                for target in 0..2 {
                    let source = if target == 0 {
                        &owner.singlet
                    } else {
                        owner.triplet.triplet()
                    };
                    let (lo, hi) = source.fluorescence().excitation_bounds_nm();
                    if values[lane][target] > 0.0 && !(lo..=hi).contains(&wavelengths[lane]) {
                        return Err(invalid(
                            "recycled absorption exceeds original excitation wavelength support",
                        ));
                    }
                }
            }
            sections.push(values);
        }
        let external = match &owner.beam {
            Some(beam) => beam.recycling_paths(owner, thermal, checkpoint)?,
            None => Vec::new(),
        };
        let mut paths = Vec::new();
        let mut path_bins = Vec::new();
        let ray_weight = 1.0 / f64::from(self.geometry.config.histories_per_emitter);
        for ray in &self.geometry.rays {
            for state in 0..2 {
                let spectrum = if state == 0 {
                    owner.singlet.fluorescence()
                } else {
                    owner.triplet.triplet().fluorescence()
                };
                let pairs = if ray.lane.is_some() { 1 } else { M.div_ceil(2) };
                for pair in 0..pairs {
                    checkpoint()?;
                    let first = ray.lane.unwrap_or(2 * pair);
                    let bins = [
                        Some(first),
                        (ray.lane.is_none() && first + 1 < M).then_some(first + 1),
                    ];
                    let mut photon_energies_j = [0.0; 2];
                    let mut photon_fractions = [0.0; 2];
                    for lane in 0..2 {
                        let bin = bins[lane].unwrap_or(first);
                        photon_energies_j[lane] = finite(super::HC_J_NM / wavelengths[bin])?;
                        if bins[lane].is_some() {
                            photon_fractions[lane] =
                                finite(spectrum.photon_fractions()[bin] * ray_weight)?;
                            if spectrum.photon_fractions()[bin] > 0.0
                                && photon_fractions[lane] == 0.0
                            {
                                return Err(invalid(
                                    "positive sampled emission weight is below representation",
                                ));
                            }
                        }
                    }
                    let mut segments = Vec::with_capacity(ray.cells.len());
                    for cell in &ray.cells {
                        checkpoint()?;
                        let mut optical_depth_per_ground = [[0.0; 2]; 2];
                        let mut passive_optical_depth = [0.0; 2];
                        for lane in 0..2 {
                            if let Some(bin) = bins[lane] {
                                passive_optical_depth[lane] = cell.passive_depth[bin];
                                for target in 0..2 {
                                    optical_depth_per_ground[lane][target] =
                                        super::super::rate_power(
                                            sections[cell.cell][bin][target],
                                            cell.length_m,
                                            cell.volume_m3,
                                        )?;
                                    if sections[cell.cell][bin][target] > 0.0
                                        && optical_depth_per_ground[lane][target] == 0.0
                                    {
                                        return Err(invalid(
                                            "positive recycled molecular opacity is below representation",
                                        ));
                                    }
                                }
                            }
                        }
                        segments.push(SaturablePumpSegment {
                            cell_index: cell.cell,
                            optical_depth_per_ground,
                            passive_optical_depth,
                        });
                    }
                    paths.push(SaturableEmissionPath {
                        source_cell: ray.source_cell,
                        source_state: state,
                        photon_fractions,
                        path: SaturablePumpPath {
                            photon_energies_j,
                            segments,
                        },
                        terminal: EmissionTerminal::Escape,
                    });
                    path_bins.push(bins);
                }
            }
        }
        let solver = preserve_checkpoint(checkpoint, |cancelled| {
            SaturablePhotonRecycling::try_new(
                &models,
                &external,
                &paths,
                domain.fraction_roundoff,
                models.len(),
                domain.maximum_paths,
                domain.maximum_total_segments,
                cancelled,
            )
        })?;
        checkpoint()?;
        Ok(Arc::new(PreparedRecycling {
            identity: owner.source_identity,
            temperatures,
            solver,
            path_bins,
        }))
    }

    pub(super) fn accept(&mut self, prepared: Arc<PreparedRecycling>) {
        self.prepared = Some(prepared);
    }
}

fn transport_error(error: VolumeTransportError) -> OpticalHeatingError {
    match error {
        VolumeTransportError::Cancelled => OpticalHeatingError::Thermal(DepositionError::Cancelled),
        _ => invalid("finite recycling geometry, observation or traversal budget"),
    }
}
