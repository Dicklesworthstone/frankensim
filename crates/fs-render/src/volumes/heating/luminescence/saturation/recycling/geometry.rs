//! Reusable complete emission paths through the original DDA geometry.

use super::*;
use crate::volumes::transport::positive_product;
use crate::volumes::transport::scattering::exchange::emitted_ray;
use crate::volumes::transport::scattering::population::PopulationPathSampler;
use fs_rand::StreamKey;

#[allow(clippy::too_many_arguments)]
pub(super) fn compile<const M: usize>(
    owner: &SaturableLuminescentHeating<M>,
    cx: &Cx<'_>,
    medium: &AbsorbingGrid<'_>,
    scattering: Option<(&ScatteringGrid<'_>, u32)>,
    sources: &[RecyclingEmitterCrossSections<M>],
    config: PhotonRecyclingGeometry,
    limits: SaturableRecyclingLimits,
    hash: &mut DomainHasher,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<Geometry<M>, OpticalHeatingError> {
    if config.histories_per_emitter == 0
        || config.maximum_paths == 0
        || config.maximum_total_segments == 0
        || !config.minimum_temperature_k.is_finite()
        || config.minimum_temperature_k <= 0.0
        || !config.maximum_temperature_k.is_finite()
        || config.maximum_temperature_k < config.minimum_temperature_k
        || !config.fraction_roundoff.is_finite()
        || !(0.0..1.0).contains(&config.fraction_roundoff)
        || sources.len() != owner.emitters.len()
    {
        return Err(invalid(
            "invalid finite recycling sampling, domain or source bindings",
        ));
    }
    let lane_count = if scattering.is_some() { M } else { 1 };
    let copies = 2usize
        .checked_mul(if scattering.is_some() { 1 } else { M.div_ceil(2) })
        .ok_or_else(|| invalid("emission spectrum path count overflow"))?;
    let ray_count = owner
        .emitters
        .len()
        .checked_mul(config.histories_per_emitter as usize)
        .and_then(|count| count.checked_mul(lane_count))
        .ok_or_else(|| invalid("emission ray count overflow"))?;
    let path_count = ray_count
        .checked_mul(copies)
        .ok_or_else(|| invalid("emission path count overflow"))?;
    if path_count > config.maximum_paths {
        return Err(invalid("emitted rays exceed aggregate path allowance"));
    }
    let mut bindings = BTreeMap::new();
    for source in sources {
        checkpoint()?;
        if bindings
            .insert(source.grid_cell(), source.clone())
            .is_some()
        {
            return Err(OpticalHeatingError::DuplicateCell(source.grid_cell()));
        }
    }
    let mut cross_sections = Vec::with_capacity(owner.emitters.len());
    for emitter in &owner.emitters {
        let source = bindings
            .remove(&emitter.binding.grid_cell)
            .ok_or(OpticalHeatingError::UnmappedCell(emitter.binding.grid_cell))?;
        hash.update(&(emitter.binding.grid_cell as u64).to_le_bytes());
        hash.update(source.identity().as_bytes());
        match &source.source {
            CrossSections::Authored { values, .. } => {
                for value in values.iter().flatten() {
                    hash.update(&value.to_bits().to_le_bytes());
                }
            }
            #[cfg(feature = "material-state")]
            CrossSections::Material { .. } => {}
        }
        cross_sections.push(source);
    }
    for cell in 0..medium.cell_count() {
        checkpoint()?;
        // Native thermal admission is required even for currently unseen
        // passive absorbers; another heat source can change their validity.
        medium
            .coefficients_at_cell(cell, owner.singlet.fluorescence().wavelengths_nm())
            .map_err(transport_error)?;
        if medium.cell_density_kg_m3(cell).is_some_and(|rho| rho > 0.0)
            && owner
                .emitters
                .binary_search_by_key(&cell, |e| e.binding.grid_cell)
                .is_err()
        {
            return Err(OpticalHeatingError::UnmappedCell(cell));
        }
    }
    for value in [
        config.histories_per_emitter as u64,
        config.maximum_paths as u64,
        config.maximum_total_segments as u64,
        config.seed,
        config.source_id,
    ] {
        hash.update(&value.to_le_bytes());
    }
    for value in [
        config.minimum_temperature_k,
        config.maximum_temperature_k,
        config.fraction_roundoff,
    ] {
        hash.update(&value.to_bits().to_le_bytes());
    }
    let sampler = scattering
        .map(|(medium, _)| {
            PopulationPathSampler::try_new(cx, medium, *owner.singlet.fluorescence().wavelengths_nm())
        })
        .transpose()
        .map_err(transport_error)?;
    let mut rays = Vec::new();
    rays.try_reserve_exact(ray_count)
        .map_err(|_| invalid("emission ray allocation budget"))?;
    let mut volumes_m3 = vec![0.0; owner.emitters.len()];
    let mut retained = 0usize;
    for (source_cell, emitter) in owner.emitters.iter().enumerate() {
        for lane in 0..lane_count {
            let selected_lane = scattering.map(|_| lane);
            let mut key = DomainHasher::new(if scattering.is_some() {
                "org.frankensim.render.finite-scattered-emission-rays.v1"
            } else {
                "org.frankensim.render.finite-emission-rays.v1"
            });
            for value in [
                config.seed,
                config.source_id,
                emitter.binding.grid_cell as u64,
            ] {
                key.update(&value.to_le_bytes());
            }
            if let Some(lane) = selected_lane {
                key.update(&(lane as u64).to_le_bytes());
            }
            let seed = u64::from_le_bytes(
                key.finalize().as_bytes()[..8]
                    .try_into()
                    .expect("eight hash bytes"),
            );
            let mut stream = StreamKey {
                seed,
                kernel: 0x6672_6379,
                tile: 0,
            }
            .stream();
            for _ in 0..config.histories_per_emitter {
                checkpoint()?;
                let (origin, direction) =
                    emitted_ray(medium, emitter.binding.grid_cell, &mut || stream.next_f64())
                        .map_err(transport_error)?;
                let remaining = (config.maximum_total_segments - retained) / copies;
                let path = if let Some(sampler) = &sampler {
                    sampler.sample_path(
                        cx,
                        lane,
                        origin,
                        direction,
                        scattering.expect("prepared scattering snapshot").1,
                        remaining,
                        &mut || stream.next_f64(),
                    )
                } else {
                    medium.population_beam_cell_path(cx, origin, direction, f64::INFINITY, remaining)
                }
                .map_err(transport_error)?;
                retained = retained
                    .checked_add(
                        path.len()
                            .checked_mul(copies)
                            .ok_or_else(|| invalid("emission segment count overflow"))?,
                    )
                    .ok_or_else(|| invalid("emission segment count overflow"))?;
                let volume = path
                    .first()
                    .filter(|s| s.cell_index == emitter.binding.grid_cell)
                    .ok_or_else(|| invalid("sampled source ray lacks its original cell"))?
                    .cell_volume_m3;
                if volumes_m3[source_cell] != 0.0 && volumes_m3[source_cell] != volume {
                    return Err(invalid("inconsistent original emitter volume"));
                }
                volumes_m3[source_cell] = volume;
                let mut cells = Vec::with_capacity(path.len());
                if let Some(lane) = selected_lane {
                    hash.update(&(lane as u64).to_le_bytes());
                }
                for segment in path {
                    checkpoint()?;
                    let coefficients = if let Some(sampler) = &sampler {
                        let mut coefficients = [0.0; M];
                        coefficients[lane] = sampler
                            .passive_absorption(segment.cell_index, lane)
                            .map_err(transport_error)?;
                        coefficients
                    } else {
                        medium
                            .coefficients_at_cell(
                                segment.cell_index,
                                owner.singlet.fluorescence().wavelengths_nm(),
                            )
                            .map_err(transport_error)?
                    };
                    let density = medium
                        .cell_density_kg_m3(segment.cell_index)
                        .ok_or(OpticalHeatingError::UnmappedCell(segment.cell_index))?;
                    let mut passive_depth = [0.0; M];
                    for lane in 0..M {
                        if density > 0.0 && coefficients[lane] > 0.0 {
                            passive_depth[lane] = finite(positive_product(
                                density,
                                coefficients[lane],
                                segment.length_m,
                            ))?;
                            if passive_depth[lane] == 0.0 {
                                return Err(invalid(
                                    "positive passive emission opacity is below representation",
                                ));
                            }
                        }
                    }
                    let cell = match owner
                        .emitters
                        .binary_search_by_key(&segment.cell_index, |e| e.binding.grid_cell)
                    {
                        Ok(cell) => cell,
                        Err(_) if passive_depth == [0.0; M] => continue,
                        Err(_) => return Err(OpticalHeatingError::UnmappedCell(segment.cell_index)),
                    };
                    hash.update(&(segment.cell_index as u64).to_le_bytes());
                    for value in [segment.length_m, segment.cell_volume_m3]
                        .into_iter()
                        .chain(passive_depth)
                    {
                        hash.update(&value.to_bits().to_le_bytes());
                    }
                    cells.push(RayCell {
                        cell,
                        length_m: segment.length_m,
                        volume_m3: segment.cell_volume_m3,
                        passive_depth,
                    });
                }
                hash.update(&(cells.len() as u64).to_le_bytes());
                rays.push(EmissionRay {
                    source_cell,
                    lane: selected_lane,
                    cells,
                });
            }
        }
    }
    checkpoint()?;
    Ok(Geometry {
        config,
        cross_sections,
        volumes_m3,
        rays,
        limits,
    })
}
