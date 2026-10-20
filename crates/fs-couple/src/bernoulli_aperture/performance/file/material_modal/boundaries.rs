//! Shared material key declarations on either original acoustic side.
use super::*;
use crate::bernoulli_aperture::material_boundary::{
    BoundaryContext, MaterialBoundaryError, MaterialBoundaryRecord, ParsedBoundaries,
    Rows as BoundaryRows,
};
use crate::bernoulli_aperture::network::boundaries::BoundaryDrive;
use crate::bernoulli_aperture::performance::boundaries::{BoundaryControl, BoundaryProgram};
use crate::bernoulli_aperture::two_sided::TwoSidedAperture;
use crate::bernoulli_aperture::two_sided::observer::schedule::boundaries::{
    BoundaryForceBinding, BoundaryGestureProgram, source_pressure_path,
};
use fs_scenario::assembly::ParameterPath;
use fs_scenario::gesture::GestureSchedule;

pub(super) struct Inputs {
    sides: [ParsedBoundaries; 2],
    pressure: Option<ParameterPath>,
    sampled: [Vec<BoundaryControl>; 2],
}
pub(super) struct Programs {
    sides: [Option<BoundaryGestureProgram>; 2],
    pressure: Option<ParameterPath>,
    sampled: [Vec<BoundaryControl>; 2],
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
        let mut sampled = [Vec::new(), Vec::new()];
        let explicit_samples = next_is(rows, "boundary_groups");
        if explicit_samples {
            let count = rows.count("boundary_groups", 2)?;
            if count == 0 {
                return Err(bad(rows.line, "boundary_groups needs one or two explicit keyed sides"));
            }
            let mut selected = [false; 2];
            let mut assignments = 0usize;
            for _ in 0..count {
                let mut row = rows.row("boundary_side")?;
                let index = match row.word()? {
                    "upstream" => 0,
                    "downstream" => 1,
                    _ => return Err(bad(row.line, "boundary_side must name upstream or downstream")),
                };
                row.finish()?;
                if selected[index] {
                    return Err(bad(rows.line, "each keyed side must be declared exactly once"));
                }
                selected[index] = true;
                if rows.lines.clone().next().and_then(|line| line.split_ascii_whitespace().nth(1))
                    == Some("none") {
                    return Err(bad(rows.line + 1, "each sampled boundary group needs material moving slits"));
                }
                sides[index] = read_side(rows)?;
                let count = rows.count("boundary_controls", MAX_CONTROLS)?;
                assignments = assignments.checked_add(count)
                    .ok_or_else(|| bad(rows.line, "boundary assignment count overflow"))?;
                if assignments > MAX_CONTROLS {
                    return Err(bad(rows.line, "both sides share bounded boundary control storage"));
                }
                sampled[index].reserve(count);
                for _ in 0..count {
                    let mut row = rows.row("boundary_control")?;
                    let sample = row.parse()?;
                    let node = row.parse()?;
                    let closing_force = row.scalar()?;
                    let exterior_pressure_pa = row.scalar()?;
                    if closing_force < 0.0 {
                        return Err(bad(row.line, "boundary closing force must be nonnegative newtons"));
                    }
                    row.finish()?;
                    sampled[index].push(BoundaryControl { sample, drive: BoundaryDrive {
                        node, exterior_pressure_pa, opening_force_n: -closing_force,
                    }});
                }
            }
            if next_is(rows, "control_source") || next_is(rows, "upstream_boundaries")
                || next_is(rows, "downstream_boundaries") {
                return Err(bad(rows.line + 1,
                    "sampled boundary_groups cannot mix with typed or held boundary declarations"));
            }
        } else {
            for (index, key) in ["upstream_boundaries", "downstream_boundaries"]
            .into_iter()
            .enumerate()
            {
                if !next_is(rows, key) {
                    continue;
                }
                rows.row(key)?.finish()?;
                sides[index] = read_side(rows)?;
            }
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
        if next_is(rows, "boundary_groups") {
            return Err(bad(rows.line + 1, "one native boundary control selector is allowed"));
        }
        Ok(Self { sides, pressure, sampled })
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

    /// Attach original material keys to the cold physical core BEFORE receiver
    /// admission, at its actual lowered gas/load/clock. The existing blade cache
    /// remains the only material loader and aggregate source-byte allowance.
    #[allow(clippy::type_complexity)]
    pub(super) fn bind(
        self,
        mut system: TwoSidedAperture,
        cache: &mut MaterialCache,
        loader: &mut impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
        gate: &CancelGate,
    ) -> Result<
        (
            TwoSidedAperture,
            Programs,
            Vec<(ValveSide, MaterialBoundaryRecord)>,
        ),
        PlateValveInputError,
    > {
        let mut programs: [Option<BoundaryGestureProgram>; 2] = [None, None];
        let mut records: Vec<(ValveSide, MaterialBoundaryRecord)> = Vec::new();
        for (side, inputs) in [ValveSide::Upstream, ValveSide::Downstream]
            .into_iter()
            .zip(self.sides)
        {
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
            for record in &built.records {
                if records.iter().any(|(_, original)| original.closing_force() == record.closing_force()) {
                    return Err(bad(0, "material keys on both sides need distinct physical declaration paths"));
                }
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
            system = system.with_boundary_apertures(side, mechanics, count, gate)
                .map_err(PlateValveInputError::Physics)?;
            let index = match side {
                ValveSide::Upstream => 0,
                ValveSide::Downstream => 1,
            };
            programs[index] = Some(BoundaryGestureProgram { initial, bindings });
            records.extend(built.records.into_iter().map(|record| (side, record)));
        }
        Ok((system, Programs { sides: programs, pressure: self.pressure, sampled: self.sampled }, records))
    }
}

impl Programs {
    /// Preserve the original all-typed gesture owner or compile one explicit
    /// complete drive program per side. Never attach two competing schedules.
    pub(super) fn compile(self, source: TwoSidedPressure, schedule: GestureSchedule,
        config: AperturePerformanceConfig) -> Result<ScheduledTwoSidedPressure, PlateValveInputError>
    {
        let renderer = if let Some(pressure) = self.pressure {
            ScheduledTwoSidedPressure::with_typed_boundary_schedule(
                source, schedule, config, pressure, self.sides,
            )
            .map_err(PlateValveInputError::Gesture)?
        } else {
            let mut renderer = ScheduledTwoSidedPressure::new(source, schedule, config)
                .map_err(PlateValveInputError::Gesture)?;
            for ((side, program), controls) in [ValveSide::Upstream, ValveSide::Downstream]
                .into_iter()
                .zip(self.sides).zip(self.sampled)
            {
                if let Some(program) = program {
                    renderer = renderer
                        .with_boundary_controls(
                            side,
                            BoundaryProgram {
                                initial: program.initial,
                                controls,
                            },
                        )
                        .map_err(PlateValveInputError::Gesture)?;
                }
            }
            renderer
        };
        Ok(renderer)
    }
}

fn read_side(rows: &mut Reader<'_>) -> Result<ParsedBoundaries, PlateValveInputError> {
    if !next_is(rows, "boundaries") {
        return Err(bad(rows.line + 1, "a selected side requires its shared boundaries declaration"));
    }
    let mut shared = BoundaryRows { lines: rows.lines.clone(), line: rows.line };
    let parsed = ParsedBoundaries::read(&mut shared).map_err(error)?;
    rows.lines = shared.lines;
    rows.line = shared.line;
    Ok(parsed)
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
