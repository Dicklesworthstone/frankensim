//! Reed acoustic-axis adapter for the shared stationary gas input.

use super::{
    GasFieldResolution, GasRecord, Reader, ReedPerformanceError, ThermochemicalSources,
    read_resolution,
};
use fs_duct::Segment;
use crate::acoustic_realize::AcousticRealizeError;
use crate::gas_field_input::{SteadyGasInput, SteadyGasInputError};
pub use crate::gas_field_input::SteadyGasPerformance;

impl From<SteadyGasInputError> for ReedPerformanceError {
    fn from(error: SteadyGasInputError) -> Self {
        match error {
            SteadyGasInputError::Input { line, what } => Self::Input { line, what },
            SteadyGasInputError::Gas { line, source } => Self::Gas { line, source },
            SteadyGasInputError::Geometry(source) => Self::Model(AcousticRealizeError::Duct(source)),
            SteadyGasInputError::SteadyGas { line, source } => Self::SteadyGas { line, source },
        }
    }
}

/// gas_steady MODE RESOLUTION CELLS INJECTIONS mole SPECIES_COUNT ID...
/// gas_inlet MASS_FLOW_KG_S T_K FRACTION...
/// gas_cell END_M WALL_CONDUCTANCE_W_K WALL_T_K INITIAL_T_K (per cell)
/// gas_injection CELL MASS_FLOW_KG_S T_K FRACTION... (zero-based cell)
/// gas_solve MIN_T MAX_T H_TOL MAX_ITER DT MAX_RELAX MAX_SUBSTEPS P_TOL E_TOL MAX_MACH
///
/// The wall conductance is an explicit total W/K, never an inferred transport
/// coefficient. This one-dimensional generator has no side-branch mass balance;
/// a tone-hole network must supply its own prescribed field instead.
pub(in super::super) fn read(
    reader: &mut Reader<'_>,
    sources: Option<&ThermochemicalSources<'_>>,
    pressure_pa: f64,
    segments: &[Segment],
    records: &mut Vec<GasRecord>,
) -> Result<(SteadyGasPerformance, GasFieldResolution), ReedPerformanceError> {
    let mut row = reader.row("gas_steady")?;
    let resolution = read_resolution(&mut row)?;
    let input = SteadyGasInput::header(row.line, &mut row.words, sources)?;
    let generated = input.read(
        &mut || {
            let row = reader.next()?;
            Ok::<_, ReedPerformanceError>((row.line, row.words))
        },
        pressure_pa, segments, records, &mut || false,
    ).map_err(|error| error.into_error(Into::into))?;
    Ok((generated, resolution))
}

// Keep the existing geometry regression seam on the one shared implementation.
#[cfg(test)]
fn geometry(
    segments: &[Segment],
    start: f64,
    end: f64,
    line: usize,
    geometry_steps: &mut usize,
) -> Result<(f64, f64), ReedPerformanceError> {
    crate::gas_field_input::geometry(segments, start, end, line, geometry_steps, &mut || false)
        .map_err(Into::into)
}

#[cfg(test)]
mod tests;
