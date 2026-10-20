//! Flue adapter over the shared material-boundary grammar and physical builder.
use super::material::Sources;
use super::{FlueInputError, Rows, bad};
use crate::bernoulli_aperture::material_boundary::{
    BoundaryContext, MaterialBoundaryError, ParsedBoundaries as SharedBoundaries,
    Rows as BoundaryRows,
};
use crate::flue_loop::{
    FlueError, FlueVoice,
    network::{FlueBoundaryAperture, FlueSourcePaths},
};
use fs_exec::CancelGate;
use fs_matdb::NormalizedMaterialCardPack;
use fs_scenario::entity::EntityDeclaration;

pub use crate::bernoulli_aperture::material_boundary::{
    MaterialBoundaryRecord as FlueBoundaryRecord, boundary_closing_force_path,
};

/// Typed pressure/angle paths matching `control_source ASSEMBLY PART`.
#[must_use]
pub fn flue_source_control_paths(assembly: &str, part: &str) -> FlueSourcePaths {
    let assembly = EntityDeclaration::assembly(assembly).identity();
    FlueSourcePaths::for_entity(EntityDeclaration::part(assembly, part).identity())
}

pub(super) fn read_control_source(
    rows: &mut Rows<'_>,
) -> Result<Option<FlueSourcePaths>, FlueInputError> {
    if !rows.next_is("control_source") {
        return Ok(None);
    }
    let mut row = rows.row("control_source")?;
    let assembly = row
        .fields
        .next()
        .ok_or_else(|| bad(row.line, "missing source assembly"))?;
    let part = row
        .fields
        .next()
        .ok_or_else(|| bad(row.line, "missing source part"))?;
    if assembly.len() > 128 || part.len() > 128 {
        return Err(bad(row.line, "source declaration name exceeds 128 bytes"));
    }
    row.end()?;
    Ok(Some(flue_source_control_paths(assembly, part)))
}

fn error(source: MaterialBoundaryError) -> FlueInputError {
    match source {
        MaterialBoundaryError::Input { line, what } => FlueInputError::Input { line, what },
        MaterialBoundaryError::MaterialSource { line, source } => {
            FlueInputError::MaterialSource { line, source }
        }
        MaterialBoundaryError::Material { line, source } => {
            FlueInputError::Material { line, source }
        }
        MaterialBoundaryError::Physics(source) => {
            FlueInputError::Physics(FlueError::Physics(source))
        }
    }
}
fn source_error(source: FlueInputError) -> MaterialBoundaryError {
    match source {
        FlueInputError::Input { line, what } => MaterialBoundaryError::Input { line, what },
        FlueInputError::MaterialSource { line, source } => {
            MaterialBoundaryError::MaterialSource { line, source }
        }
        FlueInputError::Material { line, source } => {
            MaterialBoundaryError::Material { line, source }
        }
        other => MaterialBoundaryError::MaterialSource {
            line: 0,
            source: other.to_string(),
        },
    }
}
pub(super) struct ParsedBoundaries(SharedBoundaries);
impl ParsedBoundaries {
    pub fn read(rows: &mut Rows<'_>) -> Result<Self, FlueInputError> {
        let mut shared = BoundaryRows {
            lines: rows.lines.clone(),
            line: rows.line,
        };
        let parsed = SharedBoundaries::read(&mut shared).map_err(error)?;
        rows.lines = shared.lines;
        rows.line = shared.line;
        Ok(Self(parsed))
    }
    pub fn bind(
        self,
        voice: FlueVoice,
        max_steps: u64,
        sources: &mut Sources,
        loader: &mut impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
        gate: &CancelGate,
    ) -> Result<(FlueVoice, Vec<FlueBoundaryRecord>), FlueInputError> {
        let network = voice.network().expect("file construction selected a graph");
        let built = self
            .0
            .build(
                BoundaryContext {
                    graph: network.graph(),
                    gases: network.section_gases(),
                    represented: network.represented_sections(),
                    sample_period_s: voice.sample_period_s(),
                    max_steps,
                },
                &mut |path, line| sources.load(path, line, loader, gate).map_err(source_error),
                gate,
            )
            .map_err(error)?;
        if built.mechanics.is_empty() {
            return Ok((voice, built.records));
        }
        let mechanics: Vec<_> = built
            .mechanics
            .into_iter()
            .map(|key| FlueBoundaryAperture {
                boundary: key.boundary,
                closing_force: key.closing_force,
                initial_closing_force_n: key.initial_closing_force_n,
                exterior_pressure_pa: key.exterior_pressure_pa,
            })
            .collect();
        let count = mechanics.len();
        let voice = voice
            .with_boundary_apertures(mechanics, count, gate)
            .map_err(FlueInputError::Physics)?;
        Ok((voice, built.records))
    }
}
