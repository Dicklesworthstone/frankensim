//! Shared authored DG bore input for moving apertures and flue sources.
//!
//! Geometry, primitive material gas, acoustic mesh, wall quadrature and pole
//! resolution are data. All propagation and material evaluation remain with
//! their original owners. This parser opens no files and creates no oscillator.

use super::DgPortLoad;
use crate::acoustic_realize::AcousticRealizeError;
use fs_duct::dg::{
    DgBoundary, DgDuct, DgMesh, DgWallResolution, MAX_DG_ELEMENTS, MAX_DG_ORDER,
    MAX_DG_PROPAGATOR_DOFS, PortId,
};
use fs_duct::profile::gas_field::GasFieldDuctError;
use fs_duct::{Duct, Segment};
use fs_exec::CancelGate;
use fs_matdb::QueryPoint;
use fs_material::gas::GasState;
use fs_material::gas_field::{AxialGasField, MAX_FIELD_STATIONS};
use fs_material::gas_field::input::PrimitiveGasColumns;
use fs_material::gas_species::input::{GasInputError, GasRecord, ThermochemicalSources};
use std::str::SplitAsciiWhitespace;

/// Keep source-line errors distinct from the original material and solver error.
#[derive(Debug)]
pub enum DgSourceError {
    Input { line: usize, what: &'static str },
    Gas { line: usize, source: fs_material::MaterialError },
    Physics(AcousticRealizeError),
}
impl core::fmt::Display for DgSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Input { line, what } => write!(f, "DG source line {line}: {what}"),
            Self::Gas { line, source } => write!(f, "DG material line {line}: {source}"),
            Self::Physics(source) => write!(f, "DG source physics: {source}"),
        }
    }
}
impl std::error::Error for DgSourceError {}
fn bad(line: usize, what: &'static str) -> DgSourceError {
    DgSourceError::Input { line, what }
}
fn checkpoint(gate: &CancelGate, line: usize) -> Result<(), DgSourceError> {
    if gate.is_requested() { Err(bad(line, "DG source preparation cancelled")) } else { Ok(()) }
}
fn gas_error(line: usize, error: GasInputError) -> DgSourceError {
    match error {
        GasInputError::Input(what) => bad(line, what),
        GasInputError::Material(source) => DgSourceError::Gas { line, source },
    }
}
fn field_error(line: usize, error: GasFieldDuctError) -> DgSourceError {
    match error {
        GasFieldDuctError::Material(source) => DgSourceError::Gas { line, source },
        GasFieldDuctError::Duct(source) =>
            DgSourceError::Physics(AcousticRealizeError::Duct(source)),
    }
}

struct Row<'a, 'b> { line: usize, fields: &'a mut SplitAsciiWhitespace<'b> }
impl<'a, 'b> Row<'a, 'b> {
    fn word(&mut self) -> Result<&'b str, DgSourceError> {
        self.fields.next().ok_or_else(|| bad(self.line, "missing DG source field"))
    }
    fn scalar(&mut self) -> Result<f64, DgSourceError> {
        let value: f64 = self.word()?.parse()
            .map_err(|_| bad(self.line, "invalid DG numeric field"))?;
        if !value.is_finite() { return Err(bad(self.line, "DG values must be finite")); }
        Ok(value)
    }
    fn positive(&mut self) -> Result<f64, DgSourceError> {
        let value = self.scalar()?;
        if value <= 0.0 { return Err(bad(self.line, "DG dimensions and tolerances must be positive")); }
        Ok(value)
    }
    fn count(&mut self, maximum: usize) -> Result<usize, DgSourceError> {
        let value: usize = self.word()?.parse()
            .map_err(|_| bad(self.line, "invalid DG count"))?;
        if value == 0 || value > maximum {
            return Err(bad(self.line, "DG count is outside its positive source budget"));
        }
        Ok(value)
    }
    fn finish(&mut self) -> Result<(), DgSourceError> {
        if self.fields.next().is_some() {
            return Err(bad(self.line, "unexpected extra DG source field"));
        }
        Ok(())
    }
}

/// Complete source description retained until one original DG load is prepared.
pub struct DgSourceInput {
    geometry: Duct,
    field: Option<AxialGasField>,
    outlet: DgBoundary,
    mesh: DgMesh,
    poles: usize,
    band_hz: [f64; 2],
    maximum_relative_error: f64,
    exponential_tolerance: f64,
    maximum_dofs: usize,
    wall_quadrature: Option<usize>,
    line: usize,
}

impl DgSourceInput {
    /// Read the shared bore grammar after its record kind has been consumed.
    ///
    /// dg_tube L R OUTLET K P POLES FMIN FMAX WALL_ERROR EXP_ERROR MAX_DOFS
    ///
    /// dg_bore SEGMENTS OUTLET K P POLES FMIN FMAX WALL_ERROR EXP_ERROR MAX_DOFS Q
    /// followed by exactly SEGMENTS rows:
    /// dg_segment cylinder R L
    /// dg_segment cone R_IN R_OUT L
    ///
    /// The dg_bore header may append:
    /// gas_field STATIONS mole SPECIES_COUNT ID...
    /// Its station rows follow the geometry: gas_station X T FRACTION...
    /// The already declared ambient supplies common static pressure. Every
    /// local acoustic/slit coefficient comes from the bound inlet/field.
    ///
    /// Error type E lets callers retain their existing line-aware row reader.
    #[allow(clippy::too_many_arguments)]
    pub fn read<'a, E: From<DgSourceError>>(
        kind: &str,
        line: usize,
        fields: &mut SplitAsciiWhitespace<'a>,
        next_row: &mut impl FnMut() -> Result<(usize, SplitAsciiWhitespace<'a>), E>,
        ambient: &GasState,
        sources: Option<&ThermochemicalSources<'_>>,
        records: &mut Vec<GasRecord>,
        gate: &CancelGate,
    ) -> Result<Self, E> {
        checkpoint(gate, line)?;
        let mut row = Row { line, fields };
        let (segment_count, mut segments) = match kind {
            "dg_tube" => {
                let length = row.positive()?;
                let radius = row.positive()?;
                (1, vec![Segment::Cylinder { radius, length }])
            }
            "dg_bore" => (row.count(MAX_DG_ELEMENTS)?, Vec::new()),
            _ => return Err(bad(line, "expected dg_tube or dg_bore").into()),
        };
        let outlet = match row.word()? {
            "open" => DgBoundary::IdealOpen,
            "closed" => DgBoundary::Closed,
            "matched" => DgBoundary::Matched,
            _ => return Err(bad(line, "DG outlet must be open, closed or matched").into()),
        };
        let mesh = DgMesh {
            elements: row.count(MAX_DG_ELEMENTS)?,
            order: row.count(MAX_DG_ORDER)?,
        };
        let poles = row.count(fs_phs::MAX_ZK_DIFFUSIVE_POLES)?;
        let band_hz = [row.positive()?, row.positive()?];
        let maximum_relative_error = row.positive()?;
        let exponential_tolerance = row.positive()?;
        let maximum_dofs = row.count(MAX_DG_PROPAGATOR_DOFS)?;
        let wall_quadrature = if kind == "dg_bore" {
            Some(row.count(64)?)
        } else { None };
        if mesh.elements < segment_count || band_hz[0] >= band_hz[1]
            || maximum_relative_error >= 1.0 || exponential_tolerance >= 1.0
            || wall_quadrature.is_some_and(|q| q < mesh.order + 1)
        {
            return Err(bad(line, "unresolved DG mesh, band, wall rule or tolerance").into());
        }
        let columns = if row.fields.clone().next() == Some("gas_field") {
            if wall_quadrature.is_none() {
                return Err(bad(line, "a gas field needs explicit dg_bore wall quadrature").into());
            }
            row.word()?;
            let stations = row.count(MAX_FIELD_STATIONS)?;
            if stations < 2 { return Err(bad(line, "a gas field needs at least two stations").into()); }
            let columns = PrimitiveGasColumns::read(row.fields, sources, stations)
                .map_err(|error| gas_error(line, error))?;
            Some((stations, columns))
        } else { None };
        row.finish()?;
        if kind == "dg_bore" {
            segments.reserve(segment_count);
            for _ in 0..segment_count {
                checkpoint(gate, line)?;
                let (line, mut fields) = next_row()?;
                let mut row = Row { line, fields: &mut fields };
                if row.word()? != "dg_segment" {
                    return Err(bad(line, "expected original dg_segment geometry").into());
                }
                let segment = match row.word()? {
                    "cylinder" => Segment::Cylinder {
                        radius: row.positive()?, length: row.positive()?,
                    },
                    "cone" => Segment::Cone {
                        inlet_radius: row.positive()?, outlet_radius: row.positive()?,
                        length: row.positive()?,
                    },
                    _ => return Err(bad(line, "DG source segment must be cylinder or cone").into()),
                };
                row.finish()?;
                segments.push(segment);
            }
        }
        let mut source = Self {
            geometry: Duct { segments }, field: None, outlet, mesh, poles, band_hz,
            maximum_relative_error, exponential_tolerance, maximum_dofs, wall_quadrature, line,
        };
        // Admission before collecting/evaluating a potentially large field.
        source.check_state_budget(ambient)?;
        if let Some((stations, columns)) = columns {
            let mut positions = Vec::with_capacity(stations);
            let mut states = Vec::with_capacity(stations);
            let mut lines = Vec::with_capacity(stations);
            for _ in 0..stations {
                checkpoint(gate, line)?;
                let (station_line, mut fields) = next_row()?;
                let mut row = Row { line: station_line, fields: &mut fields };
                if row.word()? != "gas_station" {
                    return Err(bad(station_line, "expected primitive gas_station").into());
                }
                positions.push(row.scalar()?);
                states.push(columns.read_state(row.fields)
                    .map_err(|error| gas_error(station_line, error))?);
                row.finish()?;
                lines.push(station_line);
            }
            let (model, phase) = columns.into_parts();
            let field = AxialGasField::new(
                positions, ambient.pressure, states, model, QueryPoint::new(), phase,
            ).map_err(|source| DgSourceError::Gas { line, source })?;
            if field.length_m() != source.geometry.segments.iter().map(Segment::axial_length).sum::<f64>() {
                return Err(bad(line, "gas field must cover the exact complete DG bore length").into());
            }
            // The shared observation surface retains each original knot receipt.
            for (&position, &station_line) in field.stations_m().iter().zip(&lines) {
                let state = field.state_at(position)
                    .map_err(|source| DgSourceError::Gas { line: station_line, source })?;
                records.push(GasRecord::Thermochemical(state));
            }
            source.field = Some(field);
        }
        Ok(source)
    }

    fn check_state_budget(&self, ambient: &GasState) -> Result<(), DgSourceError> {
        // Spatial attachment reserves both sets of poles, even if a particular
        // sample has exactly zero thermal loss. It never drops a physical state
        // to satisfy this allowance. Uniform input retains its original count.
        let nodes = self.mesh.elements * (self.mesh.order + 1);
        let memory = if let Some(points) = self.wall_quadrature {
            self.mesh.elements * points * 2 * self.poles
        } else {
            nodes * (self.poles + if ambient.gamma == 1.0 { 0 } else { self.poles })
        };
        if 2 * nodes + memory > self.maximum_dofs {
            return Err(bad(self.line, "DG acoustic and wall-memory extent exceeds the declared DOF allowance"));
        }
        Ok(())
    }

    /// Build once from the original geometry and field, then bind its exact
    /// mean-port load. Tapered walls use the explicit local-cylinder Webster
    /// rung; no conical three-dimensional thermoviscous accuracy is implied.
    pub fn prepare(
        self, ambient: &GasState, dt: f64, gate: &CancelGate,
    ) -> Result<DgPortLoad, DgSourceError> {
        checkpoint(gate, self.line)?;
        self.check_state_budget(ambient)?;
        let owner = if let Some(field) = &self.field {
            DgDuct::new_with_gas_field(
                &self.geometry, field, self.mesh, DgBoundary::Port(PortId(0)), self.outlet,
            ).map_err(|error| field_error(self.line, error))?
        } else {
            DgDuct::new(
                &self.geometry, ambient, self.mesh, DgBoundary::Port(PortId(0)), self.outlet,
            ).map_err(|error| DgSourceError::Physics(AcousticRealizeError::Duct(error)))?
        };
        checkpoint(gate, self.line)?;
        let owner = if let Some(points) = self.wall_quadrature {
            owner.with_spatial_viscothermal(DgWallResolution {
                poles_per_field: self.poles,
                quadrature_points: points,
                band_rad_s: self.band_hz.map(|f| core::f64::consts::TAU * f),
                maximum_relative_error: self.maximum_relative_error,
            }).map_err(|error| field_error(self.line, error))?
        } else {
            owner.with_viscothermal(self.poles, self.band_hz, self.maximum_relative_error)
                .map_err(|error| DgSourceError::Physics(AcousticRealizeError::Duct(error)))?
        };
        checkpoint(gate, self.line)?;
        let load = DgPortLoad::new(owner, dt, self.exponential_tolerance, self.maximum_dofs)
            .map_err(DgSourceError::Physics)?;
        checkpoint(gate, self.line)?;
        Ok(load)
    }
}

#[cfg(test)]
mod tests;
