//! Full polarized mirror histories driven by the original native thermal state.
//!
//! Smooth isotropic reflection changes amplitude and phase, but not the specular
//! direction. A native mirror-only itinerary can therefore retain its geometry
//! while propagating a fresh Stokes beam through every hit. Updating only each
//! local absorptance would leave downstream incident power and polarization stale.

use super::mirror_path::{
    MirrorPathIllumination, checked_absorbed_fraction, hash_beam, packet_lane, packet_schedule,
    ratio,
};
use super::response::{response_hash, validate_packet};
use super::*;
use crate::conductor::ConductorOptics;
use crate::polarization::dielectric::scene::{MirrorPathResult, PathTermination};
use crate::polarization::dielectric::transport::{BeamMedium, PolarizedBeam};
use fs_geom::Vec3;
use std::collections::BTreeSet;

/// One native object's heat receiver and fixed or state-dependent metal optics.
///
/// A law reads this exact receiver's temperature/enthalpy and material binding.
/// The full conductor owner retains its attached coherent coating. Film and
/// substrate heat remain aggregated in this cell; no layer partition is inferred.
#[derive(Clone, Debug)]
pub struct StatefulMirrorReceiver {
    pub receiver: SurfaceHeatingReceiver,
    pub optics: SubstrateOpticsModel,
}

#[derive(Debug)]
struct BoundOptics {
    cell: CellId,
    material: MaterialId,
    optics: SubstrateOpticsModel,
    paths: Vec<usize>,
    wavelengths: Vec<f64>,
}

#[derive(Debug)]
struct MirrorStep {
    normal: Vec3,
    distance_m: f64,
    optics: usize,
    receiver: SurfaceHeatingReceiver,
}

#[derive(Debug)]
struct MirrorItinerary {
    launch: PolarizedBeam,
    steps: Vec<MirrorStep>,
    lane: usize,
    identity: ContentHash,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct MirrorPacket<const N: usize> {
    pub(super) incident_scale: f64,
    pub(super) response: SurfacePowerResponse<N>,
}

#[derive(Debug)]
struct Model<const N: usize> {
    wavelengths_nm: [f64; N],
    optics: Vec<BoundOptics>,
    paths: Vec<MirrorItinerary>,
    sources: Vec<(usize, usize)>,
}

#[derive(Clone, Debug)]
struct Cache<const N: usize> {
    optics: Vec<ConductorOptics>,
    paths: Vec<Option<Arc<[MirrorPacket<N>]>>>,
}

#[derive(Clone, Debug)]
pub(super) struct MirrorFeedbackState<const N: usize> {
    model: Arc<Model<N>>,
    cache: Arc<Cache<N>>,
}

impl<const N: usize> SurfaceOpticalHeating<N> {
    /// Heat native cells with full mirror-path feedback from their T or h.
    ///
    /// Geometry, homogeneous host media, wavelengths and launch polarization
    /// stay fixed. Every retained surface event must be a smooth mirror, apart
    /// from an optional final detector/opaque termination. Dielectric branching,
    /// rough/diffuse scattering and extinguished or collision-limited histories
    /// refuse. A surface/distance cap remains an explicit surviving-light limit.
    ///
    /// At each accepted interval START, existing conductor laws resolve the
    /// native state. Changed optics invalidate the complete affected history;
    /// native Beer transport and complex polarized reflection then recompute
    /// every downstream incident watt and heat fraction. Unchanged histories
    /// reuse packets; darkness still checks all law domains without replay.
    /// Exact pulse fluence does not remove this explicit temporal approximation.
    ///
    /// Launch schedules carry physical spectral/ensemble-weighted watts, as in
    /// `try_new_from_mirror_paths`. IDs follow path/hit order. The usual source,
    /// segment and interval budgets apply; the extra knot limit counts distinct
    /// immutable laws, each shared law once. Native conduction, caloric storage,
    /// event ordering and transactional rollback remain the existing owner's.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_from_stateful_mirror_paths(
        wavelengths_nm: [f64; N],
        paths: &[MirrorPathIllumination<'_>],
        first_source_id: u64,
        receivers: &BTreeMap<u64, StatefulMirrorReceiver>,
        thermal: &DepositionSimulation<'_>,
        budget: SurfaceHeatingBudget,
        max_total_optical_knots: usize,
        tolerance: SurfaceHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SurfaceHeatingError> {
        checkpoint()?;
        validate_packet(&wavelengths_nm, &[0.5; N])?;
        if paths.is_empty() {
            return Err(invalid("at least one native mirror itinerary is required"));
        }
        let mut count = 0usize;
        let mut segments = 0usize;
        for source in paths {
            checkpoint()?;
            admit_itinerary(source.path)?;
            packet_lane(&wavelengths_nm, source.path.source_beam().wavelength_nm())?;
            count = count
                .checked_add(source.path.mirrors.len())
                .ok_or(SurfaceHeatingError::Budget)?;
            segments = source
                .path
                .mirrors
                .len()
                .checked_mul(source.power.segment_count())
                .and_then(|n| segments.checked_add(n))
                .ok_or(SurfaceHeatingError::Budget)?;
            if count > budget.max_sources || segments > budget.max_total_power_segments {
                return Err(SurfaceHeatingError::Budget);
            }
        }
        first_source_id
            .checked_add(u64::try_from(count - 1).map_err(|_| SurfaceHeatingError::Budget)?)
            .ok_or_else(|| invalid("mirror thermal source ID range overflows"))?;
        let mut prepared = Vec::new();
        prepared
            .try_reserve_exact(count)
            .map_err(|_| SurfaceHeatingError::Budget)?;
        let mut model = Model {
            wavelengths_nm,
            optics: Vec::new(),
            paths: Vec::new(),
            sources: Vec::new(),
        };
        model
            .optics
            .try_reserve_exact(count)
            .map_err(|_| SurfaceHeatingError::Budget)?;
        model
            .paths
            .try_reserve_exact(paths.len())
            .map_err(|_| SurfaceHeatingError::Budget)?;
        model
            .sources
            .try_reserve_exact(count)
            .map_err(|_| SurfaceHeatingError::Budget)?;
        let mut cache = Cache {
            optics: Vec::new(),
            paths: Vec::new(),
        };
        cache
            .optics
            .try_reserve_exact(count)
            .map_err(|_| SurfaceHeatingError::Budget)?;
        cache
            .paths
            .try_reserve_exact(paths.len())
            .map_err(|_| SurfaceHeatingError::Budget)?;
        let mut slots = BTreeMap::new();
        let mut laws = BTreeSet::new();
        let mut knots = 0usize;
        let mut identity = DomainHasher::new("org.frankensim.render.stateful-mirror-heating.v1");
        identity.update(&(N as u64).to_le_bytes());
        for wavelength in wavelengths_nm {
            identity.update(&wavelength.to_bits().to_le_bytes());
        }
        for (path_index, source) in paths.iter().enumerate() {
            checkpoint()?;
            let launch = *source.path.source_beam();
            let lane = packet_lane(&wavelengths_nm, launch.wavelength_nm())?;
            let mut hash = DomainHasher::new("org.frankensim.render.native-mirror-itinerary.v1");
            hash_beam(&mut hash, &launch);
            hash.update(&source.path.seed().to_le_bytes());
            hash.update(&source.path.path_id().to_le_bytes());
            let mut steps = Vec::new();
            steps
                .try_reserve_exact(source.path.mirrors.len())
                .map_err(|_| SurfaceHeatingError::Budget)?;
            for hit in &source.path.mirrors {
                checkpoint()?;
                let bound = receivers.get(&hit.hit.object_id).ok_or_else(|| {
                    invalid("each replayed mirror needs an explicit native receiver and optics")
                })?;
                let (kind, optical_identity) = optics_identity(&bound.optics);
                let key = (
                    bound.receiver.thermal_cell,
                    bound.receiver.material_id,
                    kind,
                    optical_identity,
                );
                let slot = match slots.get(&key) {
                    Some(&slot) => slot,
                    None => {
                        if let SubstrateOpticsModel::Law(law) = &bound.optics {
                            if laws.insert(optical_identity) {
                                knots = knots
                                    .checked_add(law.budget_size())
                                    .ok_or(SurfaceHeatingError::Budget)?;
                                if knots > max_total_optical_knots {
                                    return Err(SurfaceHeatingError::Budget);
                                }
                            }
                        }
                        let resolved = resolve(
                            &bound.optics,
                            thermal,
                            bound.receiver.thermal_cell,
                            bound.receiver.material_id,
                        )?;
                        let slot = model.optics.len();
                        cache.optics.push(resolved);
                        model.optics.push(BoundOptics {
                            cell: bound.receiver.thermal_cell,
                            material: bound.receiver.material_id,
                            optics: bound.optics.clone(),
                            paths: Vec::new(),
                            wavelengths: Vec::new(),
                        });
                        slots.insert(key, slot);
                        slot
                    }
                };
                let optical = &mut model.optics[slot];
                if optical.paths.last() != Some(&path_index) {
                    optical.paths.push(path_index);
                }
                if !optical.wavelengths.contains(&launch.wavelength_nm()) {
                    optical.wavelengths.push(launch.wavelength_nm());
                }
                let normal = hit
                    .hit
                    .hit
                    .normal
                    .ok_or_else(|| invalid("native mirror normal is missing"))?;
                for value in [hit.propagation_distance_m(), normal.x, normal.y, normal.z] {
                    hash.update(&value.to_bits().to_le_bytes());
                }
                hash.update(&hit.hit.object_id.to_le_bytes());
                hash.update(hit.hit.geometry_identity.as_bytes());
                hash.update(hit.hit.frame_identity.as_bytes());
                let source_id = first_source_id
                    .checked_add(
                        u64::try_from(prepared.len() + steps.len())
                            .map_err(|_| SurfaceHeatingError::Budget)?,
                    )
                    .ok_or_else(|| invalid("mirror thermal source ID range overflows"))?;
                for id in [
                    source_id,
                    bound.receiver.thermal_cell.0,
                    bound.receiver.material_id.0,
                ] {
                    identity.update(&id.to_le_bytes());
                }
                identity.update(&[kind]);
                identity.update(&optical_identity);
                for t in [
                    bound.receiver.temperature_bounds_k.0,
                    bound.receiver.temperature_bounds_k.1,
                ] {
                    identity.update(&t.to_bits().to_le_bytes());
                }
                steps.push(MirrorStep {
                    normal,
                    distance_m: hit.propagation_distance_m(),
                    optics: slot,
                    receiver: bound.receiver,
                });
            }
            let itinerary = MirrorItinerary {
                launch,
                steps,
                lane,
                identity: hash.finalize(),
            };
            identity.update(itinerary.identity.as_bytes());
            identity.update(source.power.identity().as_bytes());
            let packets = replay(&itinerary, &cache.optics, wavelengths_nm, &mut checkpoint)?;
            let launch_power = packet_schedule(source.power, lane, &mut checkpoint)?;
            for (hit_index, (step, packet)) in itinerary.steps.iter().zip(&packets).enumerate() {
                let source_id = first_source_id
                    .checked_add(
                        u64::try_from(prepared.len()).map_err(|_| SurfaceHeatingError::Budget)?,
                    )
                    .ok_or_else(|| invalid("mirror thermal source ID range overflows"))?;
                prepared.push(SurfaceIlluminationSource {
                    source_id,
                    thermal_cell: step.receiver.thermal_cell,
                    material_id: step.receiver.material_id,
                    temperature_bounds_k: step.receiver.temperature_bounds_k,
                    response: packet.response,
                    power: launch_power.clone(),
                });
                model.sources.push((path_index, hit_index));
            }
            model.paths.push(itinerary);
            cache.paths.push(Some(packets.into()));
        }
        let mut result = Self::try_new(&prepared, thermal, budget, tolerance, &mut checkpoint)?;
        Arc::get_mut(&mut result.model)
            .expect("new surface model is uniquely owned")
            .identity = identity.finalize();
        result.mirror_feedback = Some(MirrorFeedbackState {
            model: Arc::new(model),
            cache: Arc::new(cache),
        });
        checkpoint()?;
        Ok(result)
    }
}

impl<const N: usize> MirrorFeedbackState<N> {
    pub(super) fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), SurfaceHeatingError> {
        for bound in &self.model.optics {
            checkpoint()?;
            if let SubstrateOpticsModel::Law(law) = &bound.optics {
                if let Some(optics) = law.validate_native(thermal, bound.cell, bound.material)? {
                    for &wavelength in &bound.wavelengths {
                        optics.complex_ior(wavelength)?;
                    }
                }
            }
        }
        Ok(())
    }

    pub(super) fn refresh(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), SurfaceHeatingError> {
        for (index, bound) in self.model.optics.iter().enumerate() {
            checkpoint()?;
            let resolved = resolve(&bound.optics, thermal, bound.cell, bound.material)?;
            if resolved != self.cache.optics[index] {
                let cache = Arc::make_mut(&mut self.cache);
                cache.optics[index] = resolved;
                for &path in &bound.paths {
                    cache.paths[path] = None;
                }
            }
        }
        Ok(())
    }

    pub(super) fn response(
        &mut self,
        source: usize,
        illuminated: bool,
        evaluations: &mut usize,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<MirrorPacket<N>>, SurfaceHeatingError> {
        let (path, hit) = self.model.sources[source];
        if self.cache.paths[path].is_none() {
            if !illuminated {
                return Ok(None);
            }
            let itinerary = &self.model.paths[path];
            *evaluations = evaluations
                .checked_add(itinerary.steps.len())
                .ok_or(SurfaceHeatingError::Budget)?;
            let packets = replay(
                itinerary,
                &self.cache.optics,
                self.model.wavelengths_nm,
                checkpoint,
            )?;
            Arc::make_mut(&mut self.cache).paths[path] = Some(packets.into());
        }
        Ok(Some(
            self.cache.paths[path]
                .as_ref()
                .expect("illuminated itinerary was resolved")[hit],
        ))
    }
}

fn optics_identity(optics: &SubstrateOpticsModel) -> (u8, [u8; 32]) {
    match optics {
        SubstrateOpticsModel::Fixed(value) => (0, *value.content_identity().as_bytes()),
        SubstrateOpticsModel::Law(law) => (1, *law.identity().as_bytes()),
    }
}

fn resolve(
    optics: &SubstrateOpticsModel,
    thermal: &DepositionSimulation<'_>,
    cell: CellId,
    material: MaterialId,
) -> Result<ConductorOptics, SurfaceHeatingError> {
    match optics {
        SubstrateOpticsModel::Fixed(value) => Ok(*value),
        SubstrateOpticsModel::Law(law) => law.resolve(thermal, cell, material),
    }
}

fn admit_itinerary(path: &MirrorPathResult) -> Result<(), SurfaceHeatingError> {
    let launch = path.source_beam();
    let count = u32::try_from(path.mirrors.len()).map_err(|_| SurfaceHeatingError::Budget)?;
    if count == 0
        || !launch.power().i.is_finite()
        || launch.power().i <= 0.0
        || !matches!(launch.medium(), BeamMedium::Homogeneous(_))
    {
        return Err(invalid(
            "mirror replay requires a positive homogeneous native launch and mirror hits",
        ));
    }
    let terminal = match path.path.termination() {
        PathTermination::Detected | PathTermination::Opaque => 1,
        PathTermination::SurfaceLimit
        | PathTermination::DistanceLimit
        | PathTermination::Escaped => 0,
        _ => {
            return Err(invalid(
                "this native path termination cannot be replayed as a mirror itinerary",
            ));
        }
    };
    if count.checked_add(terminal) != Some(path.path.surface_events())
        || path.mirrors.iter().enumerate().any(|(index, hit)| {
            hit.surface_event_index() as usize != index
                || !hit.propagation_distance_m().is_finite()
                || hit.propagation_distance_m() <= 0.0
        })
    {
        return Err(invalid(
            "stateful mirror replay refuses intervening non-mirror surface events",
        ));
    }
    Ok(())
}

fn replay<const N: usize>(
    path: &MirrorItinerary,
    optics: &[ConductorOptics],
    wavelengths_nm: [f64; N],
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<Vec<MirrorPacket<N>>, SurfaceHeatingError> {
    let mut beam = path.launch;
    let reference = beam.power().i;
    let mut packets = Vec::new();
    packets
        .try_reserve_exact(path.steps.len())
        .map_err(|_| SurfaceHeatingError::Budget)?;
    for (index, step) in path.steps.iter().enumerate() {
        checkpoint()?;
        beam = beam.advance(step.distance_m)?;
        let balance = beam.reflect_conductor_with_balance(step.normal, &optics[step.optics])?;
        let incoming = balance.incident_power_w();
        if !incoming.is_finite() || incoming <= 0.0 {
            return Err(invalid(
                "replayed mirror reference power is unresolved; rescale the launch beam or shorten its itinerary",
            ));
        }
        let mut reflected = [1.0; N];
        let mut absorbed = [0.0; N];
        reflected[path.lane] = ratio(balance.reflection().beam.power().i, incoming)?;
        absorbed[path.lane] = checked_absorbed_fraction(&balance)?;
        let mut hash = response_hash(b"stateful-native-mirror-hit", balance.optics_identity())?;
        hash.update(path.identity.as_bytes());
        hash.update(&(index as u64).to_le_bytes());
        hash_beam(&mut hash, &beam);
        packets.push(MirrorPacket {
            incident_scale: ratio(incoming, reference)?,
            response: SurfacePowerResponse::finish(
                hash,
                wavelengths_nm,
                reflected,
                [0.0; N],
                absorbed,
            )?,
        });
        beam = balance.reflection().beam;
    }
    Ok(packets)
}
