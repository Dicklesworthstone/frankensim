//! Traced mirror absorption on the original native thermal-cell timeline.
//!
//! The scene has already applied visibility, incidence, polarization, previous
//! reflections and bulk attenuation. Each hit becomes an ordinary native
//! surface heating source, with its original positive electromagnetic loss.
//! Repeated hits and distinct mirrors assigned to one cell add through the
//! existing compensated receiver reduction. Geometry and optics remain frozen
//! over the receivers' explicitly declared temperature domains.

use super::response::{response_hash, validate_packet};
use super::{
    DepositionError, DepositionSimulation, DomainHasher, SpectralPowerSchedule,
    SurfaceHeatingBudget, SurfaceHeatingError, SurfaceHeatingReceiver, SurfaceHeatingTolerance,
    SurfaceIlluminationSource, SurfaceOpticalHeating, SurfacePowerResponse, invalid,
};
use crate::polarization::dielectric::mirrors::ConductorBeamBalance;
use crate::polarization::dielectric::scene::MirrorPathResult;
use crate::polarization::dielectric::transport::PolarizedBeam;
use crate::volumes::heating::illumination::schedule::SpectralPowerSegment;
use std::collections::BTreeMap;

/// One native monochromatic history and its physical launch-power schedule.
///
/// Power already includes that history's spatial, spectral and ensemble
/// quadrature weight. Splitting one source into several histories therefore
/// requires splitting its scheduled watts; no implicit ray-count average is
/// applied. All histories in one batch use the same receiver-object namespace.
#[derive(Clone, Copy, Debug)]
pub struct MirrorPathIllumination<'a> {
    pub path: &'a MirrorPathResult,
    pub power: &'a SpectralPowerSchedule<1>,
}

impl SurfaceIlluminationSource<1> {
    /// Compile one monochromatic native mirror path into real thermal sources.
    ///
    /// `power` is physical launch power [W] versus time. The retained launch
    /// beam supplies its positive reference normalization; every hit schedule
    /// is scaled by that hit's original incident/reference power. No additional
    /// cosine, area, attenuation, Fresnel factor or ray trace is applied.
    /// `receivers` maps native geometry object IDs to existing thermal cells.
    /// Coating and opaque-substrate absorption heat that one assigned cell;
    /// this query does not infer a per-film or depth partition.
    ///
    /// IDs are consecutive from `first_source_id`, in original hit order.
    /// Several paths may be compiled and passed together to the ordinary
    /// `SurfaceOpticalHeating::try_new`, with disjoint ID ranges. Sampled paths
    /// retain their original estimator weights; the caller supplies ensemble
    /// normalization in the launch schedules before thermal integration.
    ///
    /// Optical energy totals are sums of LOCAL mirror incidences/reflections:
    /// recirculating light can cross several mirrors. Only absorbed power is
    /// deposited, exactly once per hit. Bulk/diffuse/coating-dielectric losses,
    /// detectors and unresolved surviving rays receive no extra thermal credit.
    /// There must be at least one mirror hit and resolvable positive reference
    /// powers. Any failure returns no partially prepared source set.
    pub fn from_mirror_path(
        path: &MirrorPathResult,
        first_source_id: u64,
        power: &SpectralPowerSchedule<1>,
        receivers: &BTreeMap<u64, SurfaceHeatingReceiver>,
        budget: SurfaceHeatingBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<Self>, SurfaceHeatingError> {
        checkpoint()?;
        let launch = path.source_beam();
        let reference = launch.power().i;
        if !reference.is_finite() || reference <= 0.0 || path.mirrors.is_empty() {
            return Err(invalid(
                "a positive native launch beam and at least one mirror hit are required",
            ));
        }
        let count = path.mirrors.len();
        if count > budget.max_sources
            || count
                .checked_mul(power.segment_count())
                .is_none_or(|n| n > budget.max_total_power_segments)
        {
            return Err(SurfaceHeatingError::Budget);
        }
        let mut path_hash = DomainHasher::new("org.frankensim.render.mirror-path-heating.v1");
        hash_beam(&mut path_hash, launch);
        path_hash.update(&path.seed().to_le_bytes());
        path_hash.update(&path.path_id().to_le_bytes());
        let path_identity = path_hash.finalize();
        let wavelength = launch.wavelength_nm();
        let mut sources = Vec::new();
        sources
            .try_reserve_exact(count)
            .map_err(|_| SurfaceHeatingError::Budget)?;
        for (index, hit) in path.mirrors.iter().enumerate() {
            checkpoint()?;
            let receiver = receivers.get(&hit.hit.object_id).ok_or_else(|| {
                invalid("each encountered mirror needs an explicit native thermal receiver")
            })?;
            let incoming = hit.balance.incident_power_w();
            let reflected = hit.balance.reflection().beam;
            if reflected.wavelength_nm().to_bits() != wavelength.to_bits()
                || !incoming.is_finite()
                || incoming <= 0.0
            {
                return Err(invalid(
                    "mirror hits must retain the launch wavelength and positive reference power",
                ));
            }
            let offset = u64::try_from(index).map_err(|_| SurfaceHeatingError::Budget)?;
            let source_id = first_source_id
                .checked_add(offset)
                .ok_or_else(|| invalid("mirror thermal source ID range overflows"))?;
            let mut hash = response_hash(b"native-mirror-hit", hit.balance.optics_identity())?;
            hash.update(path_identity.as_bytes());
            hash.update(&offset.to_le_bytes());
            hash.update(&hit.hit.object_id.to_le_bytes());
            hash.update(hit.hit.geometry_identity.as_bytes());
            hash.update(hit.hit.frame_identity.as_bytes());
            hash_beam(&mut hash, &reflected);
            hash.update(&incoming.to_bits().to_le_bytes());
            let absorbed = checked_absorbed_fraction(&hit.balance)?;
            let response = SurfacePowerResponse::finish(
                hash,
                [wavelength],
                [ratio(reflected.power().i, incoming)?],
                [0.0],
                [absorbed],
            )?;
            sources.push(Self {
                source_id,
                thermal_cell: receiver.thermal_cell,
                material_id: receiver.material_id,
                temperature_bounds_k: receiver.temperature_bounds_k,
                response,
                power: scaled_schedule(power, ratio(incoming, reference)?, &mut checkpoint)?,
            });
        }
        checkpoint()?;
        Ok(sources)
    }
}

impl<const N: usize> SurfaceIlluminationSource<N> {
    /// Compile multiple wavelengths and histories into one native heat owner.
    ///
    /// Each launch wavelength must exactly match a lane in the strictly ordered
    /// shared packet. Histories may have different hit counts, polarization,
    /// reference powers and pulse schedules. Their physical heat adds on one
    /// thermal timeline, including repeated visits and different colors sharing
    /// a receiver. IDs follow input-path order, then hit order within each path.
    ///
    /// Only the traced lane carries incident watts. Other response lanes are
    /// lossless padding with identically zero power, not evaluations of that
    /// material at untraced wavelengths. There is no spectral interpolation,
    /// coherent interference between histories, extra Fresnel solve or implicit
    /// averaging. Geometry and optical state remain frozen as in the single
    /// path adapter. The aggregate source/segment budget applies to the batch.
    pub fn from_mirror_paths(
        wavelengths_nm: [f64; N],
        paths: &[MirrorPathIllumination<'_>],
        first_source_id: u64,
        receivers: &BTreeMap<u64, SurfaceHeatingReceiver>,
        budget: SurfaceHeatingBudget,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<Self>, SurfaceHeatingError> {
        checkpoint()?;
        validate_packet(&wavelengths_nm, &[0.5; N])?;
        if paths.is_empty() {
            return Err(invalid("at least one native mirror history is required"));
        }
        let mut count = 0usize;
        let mut segments = 0usize;
        for input in paths {
            checkpoint()?;
            let hits = input.path.mirrors.len();
            if hits == 0 {
                return Err(invalid("each native mirror history needs at least one hit"));
            }
            packet_lane(&wavelengths_nm, input.path.source_beam().wavelength_nm())?;
            count = count.checked_add(hits).ok_or(SurfaceHeatingError::Budget)?;
            segments = hits
                .checked_mul(input.power.segment_count())
                .and_then(|n| segments.checked_add(n))
                .ok_or(SurfaceHeatingError::Budget)?;
            if count > budget.max_sources || segments > budget.max_total_power_segments {
                return Err(SurfaceHeatingError::Budget);
            }
        }
        let last_offset = u64::try_from(count - 1).map_err(|_| SurfaceHeatingError::Budget)?;
        first_source_id
            .checked_add(last_offset)
            .ok_or_else(|| invalid("mirror thermal source ID range overflows"))?;
        let mut sources = Vec::new();
        sources
            .try_reserve_exact(count)
            .map_err(|_| SurfaceHeatingError::Budget)?;
        for input in paths {
            checkpoint()?;
            let lane = packet_lane(&wavelengths_nm, input.path.source_beam().wavelength_nm())?;
            let offset = u64::try_from(sources.len()).map_err(|_| SurfaceHeatingError::Budget)?;
            let start_id = first_source_id
                .checked_add(offset)
                .ok_or_else(|| invalid("mirror thermal source ID range overflows"))?;
            // Reuse the original per-hit admission and normalization, including
            // refusal of unresolved reference loss. This performs no optics.
            for source in SurfaceIlluminationSource::from_mirror_path(
                input.path,
                start_id,
                input.power,
                receivers,
                budget,
                &mut checkpoint,
            )? {
                checkpoint()?;
                let mut reflected = [1.0; N];
                let mut transmitted = [0.0; N];
                let mut absorbed = [0.0; N];
                reflected[lane] = source.response.reflected_fraction()[0];
                transmitted[lane] = source.response.transmitted_fraction()[0];
                absorbed[lane] = source.response.absorbed_fraction()[0];
                let response = SurfacePowerResponse::finish(
                    response_hash(b"native-mirror-packet", source.response.identity())?,
                    wavelengths_nm,
                    reflected,
                    transmitted,
                    absorbed,
                )?;
                sources.push(Self {
                    source_id: source.source_id,
                    thermal_cell: source.thermal_cell,
                    material_id: source.material_id,
                    temperature_bounds_k: source.temperature_bounds_k,
                    response,
                    power: packet_schedule(&source.power, lane, &mut checkpoint)?,
                });
            }
        }
        checkpoint()?;
        Ok(sources)
    }
}

impl<const N: usize> SurfaceOpticalHeating<N> {
    /// Bind simultaneous spectral mirror histories to one native thermal clock.
    ///
    /// All schedules share the original physical time axis. The existing
    /// observer cuts steps at every source transition and sums their absorbed
    /// watts into the same native enthalpy cells. No history is retraced.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_from_mirror_paths(
        wavelengths_nm: [f64; N],
        paths: &[MirrorPathIllumination<'_>],
        first_source_id: u64,
        receivers: &BTreeMap<u64, SurfaceHeatingReceiver>,
        thermal: &DepositionSimulation<'_>,
        budget: SurfaceHeatingBudget,
        tolerance: SurfaceHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SurfaceHeatingError> {
        let sources = SurfaceIlluminationSource::from_mirror_paths(
            wavelengths_nm,
            paths,
            first_source_id,
            receivers,
            budget,
            &mut checkpoint,
        )?;
        Self::try_new(&sources, thermal, budget, tolerance, checkpoint)
    }
}

impl SurfaceOpticalHeating<1> {
    /// Bind a traced mirror path directly to the existing native thermal clock.
    ///
    /// The returned observer is consumed by `DepositionSimulation::advance_observed`.
    /// Its normal interval integration, material/temperature admission, enthalpy
    /// storage, conduction, schedule cuts and transactional rollback apply.
    /// Optical responses and visibility are frozen; no path is retraced during
    /// thermal steps. Rebuild from the optical owner when that scene changes.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_from_mirror_path(
        path: &MirrorPathResult,
        first_source_id: u64,
        power: &SpectralPowerSchedule<1>,
        receivers: &BTreeMap<u64, SurfaceHeatingReceiver>,
        thermal: &DepositionSimulation<'_>,
        budget: SurfaceHeatingBudget,
        tolerance: SurfaceHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SurfaceHeatingError> {
        let sources = SurfaceIlluminationSource::from_mirror_path(
            path,
            first_source_id,
            power,
            receivers,
            budget,
            &mut checkpoint,
        )?;
        Self::try_new(&sources, thermal, budget, tolerance, checkpoint)
    }
}

pub(super) fn packet_lane<const N: usize>(
    wavelengths_nm: &[f64; N],
    wavelength_nm: f64,
) -> Result<usize, SurfaceHeatingError> {
    wavelengths_nm
        .binary_search_by(|value| value.total_cmp(&wavelength_nm))
        .map_err(|_| invalid("each native launch wavelength must match an explicit packet lane"))
}

pub(super) fn packet_schedule<const N: usize>(
    original: &SpectralPowerSchedule<1>,
    lane: usize,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<SpectralPowerSchedule<N>, SurfaceHeatingError> {
    let mut segments = Vec::new();
    segments
        .try_reserve_exact(original.segment_count())
        .map_err(|_| SurfaceHeatingError::Budget)?;
    for segment in original.segments() {
        checkpoint()?;
        let mut start_power_w = [0.0; N];
        let mut end_power_w = [0.0; N];
        start_power_w[lane] = segment.start_power_w[0];
        end_power_w[lane] = segment.end_power_w[0];
        segments.push(SpectralPowerSegment {
            start_s: segment.start_s,
            end_s: segment.end_s,
            start_power_w,
            end_power_w,
        });
    }
    Ok(SpectralPowerSchedule::try_new(
        &segments,
        original.segment_count(),
    )?)
}

pub(super) fn checked_absorbed_fraction(
    balance: &ConductorBeamBalance,
) -> Result<f64, SurfaceHeatingError> {
    let absorbed = ratio(balance.absorbed_power_w(), balance.incident_power_w())?;
    let native = balance
        .absorbed_fraction()
        .ok_or_else(|| invalid("a dark mirror has no physical power normalization"))?;
    if (absorbed - native).abs() > 128.0 * f64::EPSILON * native {
        return Err(invalid(
            "reference mirror loss is unresolved; rescale the launch beam and retrace",
        ));
    }
    Ok(absorbed)
}

pub(super) fn ratio(numerator: f64, denominator: f64) -> Result<f64, SurfaceHeatingError> {
    let value = numerator / denominator;
    if !value.is_finite() || value < 0.0 || (numerator > 0.0 && value == 0.0) {
        return Err(invalid(
            "native mirror power normalization is not representable",
        ));
    }
    Ok(value)
}

fn scaled_schedule(
    original: &SpectralPowerSchedule<1>,
    scale: f64,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<SpectralPowerSchedule<1>, SurfaceHeatingError> {
    if scale == 1.0 {
        return Ok(original.clone());
    }
    let mut segments = Vec::new();
    segments
        .try_reserve_exact(original.segment_count())
        .map_err(|_| SurfaceHeatingError::Budget)?;
    for segment in original.segments() {
        checkpoint()?;
        let scaled = |value: f64| {
            let result = value * scale;
            if !result.is_finite() || (value > 0.0 && result == 0.0) {
                Err(invalid(
                    "native mirror incident schedule is not representable",
                ))
            } else {
                Ok([result])
            }
        };
        segments.push(SpectralPowerSegment {
            start_s: segment.start_s,
            end_s: segment.end_s,
            start_power_w: scaled(segment.start_power_w[0])?,
            end_power_w: scaled(segment.end_power_w[0])?,
        });
    }
    Ok(SpectralPowerSchedule::try_new(
        &segments,
        original.segment_count(),
    )?)
}

pub(super) fn hash_beam(hash: &mut DomainHasher, beam: &PolarizedBeam) {
    let origin = beam.ray().origin;
    let k = beam.frame().direction();
    let h = beam.frame().horizontal();
    let s = beam.power();
    for value in [
        beam.wavelength_nm(),
        origin.x,
        origin.y,
        origin.z,
        k.x,
        k.y,
        k.z,
        h.x,
        h.y,
        h.z,
        s.i,
        s.q,
        s.u,
        s.v,
    ] {
        hash.update(&value.to_bits().to_le_bytes());
    }
}

#[cfg(test)]
mod tests;
