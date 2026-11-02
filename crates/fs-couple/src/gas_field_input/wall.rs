//! Authored stationary wall cells for both original acoustic performance paths.
//! Geometry is taken from the same bore; material and transport laws are not
//! reimplemented here. A conical or profiled cell requires its own wall law.

use super::{SteadyGasInputError, SteadyGasReadError, bad, poll, row};
use crate::material_source::thermal::{MaterialWallSource, ResolvedMaterialWallSource};
use crate::material_source::{MAX_RECORD_AXES, MaterialRecord};
use fs_duct::Segment;
use fs_matdb::NormalizedMaterialCardPack;
use fs_material::gas_field::steady::wall::MaterialGasWall;
use fs_material::gas_field::steady::{
    SteadyAxialGas, SteadyGasError, SteadyGasSolution, SteadyGasSolve,
};
use fs_material::phase::SolidLiquidPhase;
use fs_material::thermal::conduction::{LayerConductionConfig, ReservoirLayerGeometry};
use fs_thermal::enthalpy::DepositionError;
use std::{collections::BTreeMap, str::SplitAsciiWhitespace, sync::Arc};

/// The actual source and geometry used to replace one cell's numeric wall law.
#[derive(Clone, Debug)]
pub struct SteadyGasMaterialWall {
    pub cell: usize,
    pub source: ResolvedMaterialWallSource,
    pub inner_radius_m: f64,
    pub thickness_m: f64,
    pub length_m: f64,
    pub film_resistance_m2_k_w: f64,
    line: usize,
}

#[derive(Default)]
pub(super) struct Sources {
    packs: BTreeMap<String, Arc<NormalizedMaterialCardPack>>,
    bytes: usize,
}

fn source_error(line: usize, error: impl core::fmt::Display) -> SteadyGasInputError {
    SteadyGasInputError::MaterialSource {
        line,
        source: error.to_string(),
    }
}

impl SteadyGasMaterialWall {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn read<'a, E>(
        next: &mut impl FnMut() -> Result<(usize, SplitAsciiWhitespace<'a>), E>,
        cell: usize,
        segments: &[Segment],
        start: f64,
        end: f64,
        sources: &mut Sources,
        load: &mut dyn FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, SteadyGasReadError<E>> {
        let mut material = row(next, "gas_wall_material")?;
        let line = material.line;
        let tokens: Vec<_> = material.words.by_ref().take(MAX_RECORD_AXES + 3).collect();
        material.finish()?;
        let record = MaterialRecord::parse(&tokens).map_err(|error| source_error(line, error))?;
        let mut row = row(next, "gas_wall_source")?;
        let temperature_axis = row.word()?.to_owned();
        let phase = match row.word()? {
            "solid" => SolidLiquidPhase::Solid,
            "liquid" => SolidLiquidPhase::Liquid,
            _ => {
                return Err(bad(
                    row.line,
                    "wall source requires an explicit solid or liquid phase",
                )
                .into());
            }
        };
        let property_names = [
            row.word()?.to_owned(),
            row.word()?.to_owned(),
            row.word()?.to_owned(),
        ];
        let temperature_span_k = [row.positive()?, row.positive()?];
        let caloric_interpolation_tolerance_j_kg = row.positive()?;
        row.finish()?;
        let source = MaterialWallSource {
            material: record,
            temperature_axis,
            phase,
            property_names,
            temperature_span_k,
            caloric_interpolation_tolerance_j_kg,
        };
        let mut geometry = super::row(next, "gas_wall_geometry")?;
        if geometry.word()? != "cylindrical" {
            return Err(bad(geometry.line, "material wall geometry must be cylindrical").into());
        }
        let thickness_m = geometry.positive()?;
        let film_resistance_m2_k_w = geometry.scalar()?;
        if film_resistance_m2_k_w < 0.0 {
            return Err(bad(
                geometry.line,
                "gas film resistance must be nonnegative in m2 K/W",
            )
            .into());
        }
        geometry.finish()?;
        let mut offset = 0.0;
        let mut inner_radius = None;
        for segment in segments {
            poll(line, cancelled)?;
            let stop = offset + segment.axial_length();
            if stop > start && offset < end {
                let Segment::Cylinder { radius, .. } = segment else {
                    return Err(bad(
                        geometry.line,
                        "material wall cells must lie in constant-radius cylindrical bore geometry",
                    )
                    .into());
                };
                if inner_radius.is_some_and(|prior| prior != *radius) {
                    return Err(bad(
                        geometry.line,
                        "material wall cells cannot span a change in cylindrical bore radius",
                    )
                    .into());
                }
                inner_radius = Some(*radius);
            }
            offset = stop;
        }
        let inner_radius_m =
            inner_radius.ok_or_else(|| bad(geometry.line, "wall cell has no bore geometry"))?;
        // Admit the physical shell before loading any material artifact.
        ReservoirLayerGeometry::cylindrical(
            inner_radius_m,
            inner_radius_m + thickness_m,
            end - start,
        )
        .map_err(|error| source_error(geometry.line, error))?;
        poll(line, cancelled)?;
        if !sources.packs.contains_key(&source.material.pack) {
            let pack = load(&source.material.pack).map_err(|error| source_error(line, error))?;
            poll(line, cancelled)?;
            let bytes = pack.to_bytes().len();
            if bytes > 8 * 1024 * 1024
                || bytes > (64 * 1024 * 1024_usize).saturating_sub(sources.bytes)
            {
                return Err(bad(
                    line,
                    "wall sources exceed 8 MiB per pack or 64 MiB combined",
                )
                .into());
            }
            sources.bytes += bytes;
            sources
                .packs
                .insert(source.material.pack.clone(), Arc::new(pack));
        }
        let pack = Arc::clone(&sources.packs[&source.material.pack]);
        let source = ResolvedMaterialWallSource::new(source, pack)
            .map_err(|error| source_error(line, error))?;
        poll(line, cancelled)?;
        Ok(Self {
            cell,
            source,
            inner_radius_m,
            thickness_m,
            length_m: end - start,
            film_resistance_m2_k_w,
            line,
        })
    }
}

pub(super) fn solve(
    problem: &SteadyAxialGas,
    walls: &[SteadyGasMaterialWall],
    wall_controls: LayerConductionConfig,
    controls: SteadyGasSolve,
    line: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SteadyGasSolution, SteadyGasInputError> {
    let layers = walls
        .iter()
        .map(|wall| {
            wall.source
                .layer(|| {
                    if cancelled() {
                        Err(DepositionError::Cancelled)
                    } else {
                        Ok(())
                    }
                })
                .map_err(|error| source_error(wall.line, error))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let bindings = walls
        .iter()
        .zip(&layers)
        .map(|(wall, layer)| {
            let geometry = ReservoirLayerGeometry::cylindrical(
                wall.inner_radius_m,
                wall.inner_radius_m + wall.thickness_m,
                wall.length_m,
            )
            .map_err(|error| source_error(wall.line, error))?;
            MaterialGasWall::new(
                wall.cell,
                layer,
                geometry,
                wall.film_resistance_m2_k_w / geometry.surface_area_m2(),
            )
            .map_err(|source| SteadyGasInputError::SteadyGas {
                line: wall.line,
                source,
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    problem
        .solve_with_material_walls(&bindings, wall_controls, controls, cancelled)
        .map_err(|source: SteadyGasError| SteadyGasInputError::SteadyGas { line, source })
}

#[cfg(test)]
mod tests;
