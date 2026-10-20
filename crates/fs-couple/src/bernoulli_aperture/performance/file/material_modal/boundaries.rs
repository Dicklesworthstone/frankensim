//! Shared material key declarations on either original acoustic side.
use super::*;
use crate::bernoulli_aperture::material_boundary::{
    BoundaryContext, MaterialBoundaryError, MaterialBoundaryRecord, ParsedBoundaries,
    Rows as BoundaryRows,
};
use crate::bernoulli_aperture::network::boundaries::BoundaryDrive;
use crate::bernoulli_aperture::performance::boundaries::BoundaryProgram;
use crate::bernoulli_aperture::two_sided::observer::schedule::boundaries::{
    BoundaryForceBinding, BoundaryGestureProgram, source_pressure_path,
};
use fs_scenario::assembly::ParameterPath;
use fs_scenario::gesture::GestureSchedule;

pub(super) struct Inputs {
    sides: [ParsedBoundaries; 2],
    pressure: Option<ParameterPath>,
}
fn next_is(rows: &Reader<'_>, key: &str) -> bool {
    rows.lines
        .clone()
        .next()
        .and_then(|line| line.split_ascii_whitespace().next())
        == Some(key)
}
impl Inputs {
    pub(super) fn read(rows: &mut Reader<'_>) -> Result<Self, PlateValveInputError> {
        let mut sides = [ParsedBoundaries::none(), ParsedBoundaries::none()];
        for (index, key) in ["upstream_boundaries", "downstream_boundaries"]
            .into_iter()
            .enumerate()
        {
            if !next_is(rows, key) {
                continue;
            }
            rows.row(key)?.finish()?;
            if !next_is(rows, "boundaries") {
                return Err(bad(
                    rows.line,
                    "a selected boundary side requires its shared boundaries declaration",
                ));
            }
            let mut shared = BoundaryRows {
                lines: rows.lines.clone(),
                line: rows.line,
            };
            sides[index] = ParsedBoundaries::read(&mut shared).map_err(error)?;
            rows.lines = shared.lines;
            rows.line = shared.line;
        }
        let pressure = if next_is(rows, "control_source") {
            let mut row = rows.row("control_source")?;
            let assembly = row.word()?;
            let part = row.word()?;
            if assembly.len() > 128 || part.len() > 128 {
                return Err(bad(row.line, "control source names exceed 128 bytes"));
            }
            row.finish()?;
            Some(source_pressure_path(assembly, part))
        } else {
            None
        };
        Ok(Self { sides, pressure })
    }

    pub(super) fn decode_schedule(
        &self,
        text: &str,
    ) -> Result<GestureSchedule, PlateValveInputError> {
        super::super::decode_schedule_with_track_limit(
            text,
            if self.pressure.is_some() { 129 } else { 1 },
        )
    }

    /// Build original material keys at the actual lowered impedance/clock, then
    /// hand all controls to the retained mechanical/receiver scheduler. The
    /// existing blade cache is the only material loader and aggregate budget.
    pub(super) fn compile(
        self,
        mut source: TwoSidedPressure,
        schedule: GestureSchedule,
        config: AperturePerformanceConfig,
        cache: &mut MaterialCache,
        loader: &mut impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
        gate: &CancelGate,
    ) -> Result<
        (
            ScheduledTwoSidedPressure,
            Vec<(ValveSide, MaterialBoundaryRecord)>,
        ),
        PlateValveInputError,
    > {
        let mut programs: [Option<BoundaryGestureProgram>; 2] = [None, None];
        let mut records = Vec::new();
        for (side, inputs) in [ValveSide::Upstream, ValveSide::Downstream]
            .into_iter()
            .zip(self.sides)
        {
            let system = source.system();
            let built = inputs
                .build(
                    BoundaryContext {
                        graph: system.graph(side),
                        gases: system.section_gases(side),
                        represented: system.represented_sections(side),
                        sample_period_s: system.valve().time_step_s(),
                        max_steps: system.valve().max_steps(),
                    },
                    &mut |path, line| cache.load(path, line, loader, gate).map_err(source_error),
                    gate,
                )
                .map_err(error)?;
            if built.mechanics.is_empty() {
                continue;
            }
            let mut initial = Vec::with_capacity(built.mechanics.len());
            let mut bindings = Vec::with_capacity(built.mechanics.len());
            let mut mechanics = Vec::with_capacity(built.mechanics.len());
            for key in built.mechanics {
                initial.push(BoundaryDrive {
                    node: key.boundary.node,
                    exterior_pressure_pa: key.exterior_pressure_pa,
                    opening_force_n: -key.initial_closing_force_n,
                });
                bindings.push(BoundaryForceBinding {
                    node: key.boundary.node,
                    closing_force: key.closing_force,
                });
                mechanics.push(key.boundary);
            }
            let count = mechanics.len();
            source = source
                .with_boundary_apertures(side, mechanics, count)
                .map_err(|e| PlateValveInputError::Gesture(GestureCompileError::Render(e)))?;
            let index = match side {
                ValveSide::Upstream => 0,
                ValveSide::Downstream => 1,
            };
            programs[index] = Some(BoundaryGestureProgram { initial, bindings });
            records.extend(built.records.into_iter().map(|record| (side, record)));
        }
        let renderer = if let Some(pressure) = self.pressure {
            ScheduledTwoSidedPressure::with_typed_boundary_schedule(
                source, schedule, config, pressure, programs,
            )
            .map_err(PlateValveInputError::Gesture)?
        } else {
            let mut renderer = ScheduledTwoSidedPressure::new(source, schedule, config)
                .map_err(PlateValveInputError::Gesture)?;
            for (side, program) in [ValveSide::Upstream, ValveSide::Downstream]
                .into_iter()
                .zip(programs)
            {
                if let Some(program) = program {
                    renderer = renderer
                        .with_boundary_controls(
                            side,
                            BoundaryProgram {
                                initial: program.initial,
                                controls: Vec::new(),
                            },
                        )
                        .map_err(PlateValveInputError::Gesture)?;
                }
            }
            renderer
        };
        Ok((renderer, records))
    }
}

fn error(source: MaterialBoundaryError) -> PlateValveInputError {
    match source {
        MaterialBoundaryError::Input { line, what } => PlateValveInputError::Input { line, what },
        MaterialBoundaryError::MaterialSource { line, source } => {
            PlateValveInputError::MaterialSource { line, source }
        }
        MaterialBoundaryError::Material { line, source } => {
            PlateValveInputError::Material { line, source }
        }
        MaterialBoundaryError::Physics(source) => PlateValveInputError::Physics(source),
    }
}
fn source_error(source: PlateValveInputError) -> MaterialBoundaryError {
    match source {
        PlateValveInputError::Input { line, what } => MaterialBoundaryError::Input { line, what },
        PlateValveInputError::MaterialSource { line, source } => {
            MaterialBoundaryError::MaterialSource { line, source }
        }
        PlateValveInputError::Material { line, source } => {
            MaterialBoundaryError::Material { line, source }
        }
        PlateValveInputError::Physics(source) => MaterialBoundaryError::Physics(source),
        other => MaterialBoundaryError::MaterialSource {
            line: 0,
            source: other.to_string(),
        },
    }
}
