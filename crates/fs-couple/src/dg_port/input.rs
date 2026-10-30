//! Shared authored DG bore input for moving apertures and flue sources.
//!
//! Geometry, primitive material gas, acoustic mesh, wall quadrature and pole
//! resolution are data. All propagation and material evaluation remain with
//! their original owners. This parser opens no files and creates no oscillator.

use super::{
    DgBranch, DgPortLoad, DgTerminalChannelInput, DgTerminalChannelKind,
    MAX_DG_TERMINAL_TRANSVERSE_CELLS,
};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::radiation::environment::multiport::bake::MAX_BAFFLED_MIMO_PORTS;
use crate::bernoulli_aperture::radiation::environment::multiport::input::MimoRadiationInputError;
use crate::bernoulli_aperture::radiation::environment::multiport::source::ComputedBaffledInput;
use crate::pcm_wav::baffled::openings::CircularOpening;
use fs_duct::dg::{
    DgBoundary, DgDuct, DgMesh, DgSidePort, DgWallResolution, MAX_DG_ELEMENTS, MAX_DG_ORDER,
    MAX_DG_PROPAGATOR_DOFS, MAX_DG_SIDE_PORTS, PortId,
};
use fs_duct::profile::gas_field::GasFieldDuctError;
use fs_duct::{Duct, Segment};
use fs_exec::CancelGate;
use fs_matdb::QueryPoint;
use fs_material::gas::GasState;
use fs_material::gas_field::input::PrimitiveGasColumns;
use fs_material::gas_field::{AxialGasField, MAX_FIELD_STATIONS};
use fs_material::gas_species::input::{
    GasInputError, GasRecord, ThermochemicalSources, read_gas_record,
};
use crate::gas_field_input::{SteadyGasInput, SteadyGasInputError};
pub use crate::gas_field_input::SteadyGasPerformance;
use std::str::SplitAsciiWhitespace;

/// Shared source-backed solid wall declaration for every authored DG wind path.
pub mod wall;

/// Keep source-line errors distinct from the original material and solver error.
#[derive(Debug)]
pub enum DgSourceError {
    Input { line: usize, what: &'static str },
    Gas { line: usize, source: fs_material::MaterialError },
    SteadyGas { line: usize, source: fs_material::gas_field::steady::SteadyGasError },
    Material { line: usize, source: fs_material::state_point::MaterialStatePointError },
    Physics(AcousticRealizeError),
    Radiation { line: usize, source: MimoRadiationInputError },
}
impl core::fmt::Display for DgSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Input { line, what } => write!(f, "DG source line {line}: {what}"),
            Self::Gas { line, source } => write!(f, "DG material line {line}: {source}"),
            Self::SteadyGas { line, source } => write!(f, "DG steady gas line {line}: {source}"),
            Self::Material { line, source } => write!(f, "DG solid material line {line}: {source}"),
            Self::Physics(source) => write!(f, "DG source physics: {source}"),
            Self::Radiation { line, source } => write!(f, "DG radiation line {line}: {source}"),
        }
    }
}
impl std::error::Error for DgSourceError {}
impl From<SteadyGasInputError> for DgSourceError {
    fn from(error: SteadyGasInputError) -> Self {
        match error {
            SteadyGasInputError::Input { line, what } => Self::Input { line, what },
            SteadyGasInputError::Gas { line, source } => Self::Gas { line, source },
            SteadyGasInputError::Geometry(source) => Self::Physics(AcousticRealizeError::Duct(source)),
            SteadyGasInputError::SteadyGas { line, source } => Self::SteadyGas { line, source },
        }
    }
}
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
    steady_gas: Option<SteadyGasPerformance>,
    outlet: DgBoundary,
    radiation: Option<(ComputedBaffledInput, GasRecord)>,
    side_openings: Vec<(DgSidePort, CircularOpening)>,
    chimneys: Vec<DgChimneyInput>,
    terminal_channels: Vec<DgTerminalChannelInput>,
    maximum_terminal_transverse_cells: usize,
    mesh: DgMesh,
    poles: usize,
    band_hz: [f64; 2],
    maximum_relative_error: f64,
    exponential_tolerance: f64,
    maximum_dofs: usize,
    wall_quadrature: Option<usize>,
    molecular: Option<fs_material::relaxation::RelaxationState>,
    molecular_law: Option<fs_material::relaxation::RelaxationAbsorption>,
    line: usize,
}

/// A finite, uniformly filled cylindrical side branch. Its gas is an original
/// shared material record; its wave and wall states belong to `fs_duct`.
struct DgChimneyInput {
    port: PortId,
    radius_m: f64,
    length_m: f64,
    mesh: DgMesh,
    gas: GasRecord,
}

impl DgChimneyInput {
    fn state_dofs(&self, poles: usize) -> usize {
        let nodes = self.mesh.elements * (self.mesh.order + 1);
        nodes * (2 + poles + if self.gas.state().gamma == 1.0 { 0 } else { poles })
    }
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
    /// Alternatively append gas_steady CELLS INJECTIONS mole SPECIES_COUNT ID...
    /// and place the shared gas_inlet/gas_cell/gas_injection/gas_solve records
    /// after geometry. The original stationary heat/species owner generates
    /// the SAME field; DG mesh/order/wall quadrature remain its acoustic axes.
    /// This route requires an unbranched bore without finite terminal passages;
    /// those mean-flow/thermal volumes are not part of the axial generator.
    ///
    /// OUTLET may also be baffled, followed after geometry and any gas stations
    /// by dg_radiation computed-baffled plus the original fourteen bake fields,
    /// then dg_exterior plus one complete shared gas record. The actual outlet
    /// area supplies its circle; the exterior record retains its own material
    /// sources. Bore/exterior equilibrium pressures must agree. The original
    /// bake covers mechanical Nyquist, independent of the microphone band.
    ///
    /// Before dg_radiation, optional dg_side_openings COUNT declares up to the
    /// original DG side-port budget (32), followed by exactly COUNT rows.
    /// Up to 33 actual exterior openings, including the axial outlet, can share
    /// the radiation owner within its original state/memory budgets. A flush circular opening uses:
    /// dg_side_opening PORT_ID DOWNSTREAM_ELEMENT X_M Y_M RADIUS_M
    /// Each selects an existing interior mesh face, a unique ID above 1, and
    /// independent physical opening geometry in the common rigid baffle. Its
    /// radius is not inferred from bore area. A flush row has zero internal
    /// length. A finite cylindrical chimney instead uses:
    /// dg_side_chimney PORT_ID DOWNSTREAM_ELEMENT X_M Y_M RADIUS_M LENGTH_M K P
    /// dg_chimney_gas followed by one complete shared material gas record.
    /// It retains its own original DG propagation and viscothermal wall states,
    /// using the parent wall poles/band/tolerance and common physical clock.
    /// The total DOF allowance includes every branch. The exterior circle is
    /// the actual chimney outlet; all terminal frames use its outward flow.
    /// Internal 3D junction corrections and mean flow are absent.
    /// The axial outlet is matrix channel zero, followed by the authored order.
    ///
    /// After any side openings and before dg_radiation, an optional block
    /// declares finite fixed slit/circular channels or exact closures in series with
    /// existing terminals:
    /// dg_terminal_channels COUNT MAX_TOTAL_TRANSVERSE_CELLS
    /// dg_terminal_closed PORT_ID
    /// dg_terminal_slit PORT_ID WIDTH_M HEIGHT_M LENGTH_M CELLS MAX_MACH
    /// dg_terminal_circular PORT_ID RADIUS_M LENGTH_M CELLS MAX_MACH
    /// A closed finite chimney retains its complete DG volume and wall memory
    /// with a rigid distal boundary, and consumes no exterior matrix row.
    /// Omitted ports remain directly connected to the full exterior matrix.
    /// Each passage uses the actual terminal gas and the original unsteady transverse
    /// momentum owner. Its fluid memory has a separate explicit cell budget;
    /// it is not a resolved moving pad, annular gap or turbulent jet model.
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
            "baffled" => DgBoundary::Port(PortId(1)),
            _ => return Err(bad(line, "DG outlet must be open, closed, matched or baffled").into()),
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
        let mut steady_input = None;
        let columns = if matches!(row.fields.clone().next(), Some("gas_field" | "gas_steady")) {
            if wall_quadrature.is_none() {
                return Err(bad(line, "a gas field needs explicit dg_bore wall quadrature").into());
            }
            if row.word()? == "gas_steady" {
                steady_input = Some(SteadyGasInput::header(line, row.fields, sources)
                    .map_err(DgSourceError::from)?);
                None
            } else {
                let stations = row.count(MAX_FIELD_STATIONS)?;
                if stations < 2 { return Err(bad(line, "a gas field needs at least two stations").into()); }
                let columns = PrimitiveGasColumns::read(row.fields, sources, stations)
                    .map_err(|error| gas_error(line, error))?;
                Some((stations, columns))
            }
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
            geometry: Duct { segments }, field: None, steady_gas: None, outlet, radiation: None,
            side_openings: Vec::new(), chimneys: Vec::new(), mesh, poles, band_hz,
            terminal_channels: Vec::new(), maximum_terminal_transverse_cells: 0,
            maximum_relative_error, exponential_tolerance, maximum_dofs, wall_quadrature,
            molecular: None, molecular_law: None, line,
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
        if let Some(input) = steady_input {
            let generated = input.read(next_row, ambient.pressure, &source.geometry.segments,
                records, &mut || gate.is_requested())
                .map_err(|error| error.into_error(|source| DgSourceError::from(source).into()))?;
            source.field = Some(generated.solution().field().clone());
            source.steady_gas = Some(generated);
        }
        if matches!(outlet, DgBoundary::Port(_)) {
            checkpoint(gate, line)?;
            let (mut radiation_line, mut fields) = next_row()?;
            if fields.clone().next() == Some("dg_side_openings") {
                if source.steady_gas.is_some() {
                    return Err(bad(radiation_line,
                        "gas_steady needs an unbranched DG bore; side-branch mean flow is not specified").into());
                }
                let mut row = Row { line: radiation_line, fields: &mut fields };
                row.word()?;
                let count = row.count(MAX_DG_SIDE_PORTS)?;
                row.finish()?;
                source.side_openings.reserve(count);
                for _ in 0..count {
                    checkpoint(gate, radiation_line)?;
                    let (line, mut fields) = next_row()?;
                    let mut row = Row { line, fields: &mut fields };
                    let opening_kind = row.word()?;
                    if !matches!(opening_kind, "dg_side_opening" | "dg_side_chimney") {
                        return Err(bad(line, "expected a DG flush opening or finite chimney").into());
                    }
                    let id = PortId(row.word()?.parse::<u32>()
                        .map_err(|_| bad(line, "DG side opening requires an unsigned port ID"))?);
                    let downstream_element = row.count(mesh.elements - 1)?;
                    let aperture = CircularOpening {
                        center_m: [row.scalar()?, row.scalar()?], radius_m: row.positive()?,
                    };
                    let chimney = if opening_kind == "dg_side_chimney" {
                        Some((row.positive()?, DgMesh {
                            elements: row.count(MAX_DG_ELEMENTS)?,
                            order: row.count(MAX_DG_ORDER)?,
                        }))
                    } else { None };
                    row.finish()?;
                    if id.0 <= 1 || source.side_openings.iter().any(|(port, _)|
                        port.id == id || port.downstream_element == downstream_element)
                    {
                        return Err(bad(line, "DG side openings require unique interior faces and unique port IDs above 1").into());
                    }
                    source.side_openings.push((DgSidePort { id, downstream_element }, aperture));
                    if let Some((length_m, mesh)) = chimney {
                        checkpoint(gate, line)?;
                        let (gas_line, mut fields) = next_row()?;
                        let mut row = Row { line: gas_line, fields: &mut fields };
                        if row.word()? != "dg_chimney_gas" {
                            return Err(bad(gas_line, "finite DG chimney requires its material gas record").into());
                        }
                        let gas = read_gas_record(row.fields, None, sources)
                            .map_err(|error| gas_error(gas_line, error))?;
                        row.finish()?;
                        if gas.state().pressure.to_bits() != ambient.pressure.to_bits() {
                            return Err(bad(gas_line, "DG chimney and bore require equal equilibrium pressures").into());
                        }
                        records.push(gas.clone());
                        source.chimneys.push(DgChimneyInput {
                            port: id, radius_m: aperture.radius_m, length_m, mesh, gas,
                        });
                    }
                }
                (radiation_line, fields) = next_row()?;
            }
            if fields.clone().next() == Some("dg_terminal_channels") {
                if source.steady_gas.is_some() {
                    return Err(bad(radiation_line,
                        "gas_steady does not include finite DG terminal passages or closures").into());
                }
                let mut row = Row { line: radiation_line, fields: &mut fields };
                row.word()?;
                let count = row.count(source.side_openings.len() + 1)?;
                let allowance = row.word()?.parse::<usize>()
                    .map_err(|_| bad(radiation_line, "invalid terminal transverse-cell allowance"))?;
                if allowance > MAX_DG_TERMINAL_TRANSVERSE_CELLS {
                    return Err(bad(radiation_line, "terminal channel cells exceed the shared hard ceiling").into());
                }
                row.finish()?;
                source.maximum_terminal_transverse_cells = allowance;
                source.terminal_channels.reserve(count);
                let mut total_cells = 0usize;
                for _ in 0..count {
                    checkpoint(gate, radiation_line)?;
                    let (line, mut fields) = next_row()?;
                    let mut row = Row { line, fields: &mut fields };
                    let channel_kind = row.word()?;
                    let port = PortId(row.word()?.parse::<u32>()
                        .map_err(|_| bad(line, "terminal channel requires an unsigned port ID"))?);
                    if (port != PortId(1) && !source.side_openings.iter().any(|(side, _)| side.id == port))
                        || source.terminal_channels.iter().any(|channel| channel.port == port)
                    {
                        return Err(bad(line, "terminal channels require distinct declared exterior ports").into());
                    }
                    let kind = match channel_kind {
                        "dg_terminal_closed" => DgTerminalChannelKind::Closed,
                        "dg_terminal_slit" => {
                            let width_m = row.positive()?;
                            let height_m = row.positive()?;
                            let length_m = row.positive()?;
                            let transverse_cells = row.count(MAX_DG_TERMINAL_TRANSVERSE_CELLS)?;
                            let maximum_mach = row.positive()?;
                            if maximum_mach > 0.3 {
                                return Err(bad(line, "terminal slit Mach ceiling must be in (0, 0.3]").into());
                            }
                            DgTerminalChannelKind::Slit {
                                width_m, height_m, length_m, transverse_cells, maximum_mach,
                            }
                        },
                        "dg_terminal_circular" => DgTerminalChannelKind::Circular {
                            radius_m: row.positive()?,
                            length_m: row.positive()?,
                            radial_cells: row.count(MAX_DG_TERMINAL_TRANSVERSE_CELLS)?,
                            maximum_mach: row.positive()?,
                        },
                        _ => return Err(bad(line, "expected a closed terminal or finite slit/circular channel").into()),
                    };
                    row.finish()?;
                    let channel = DgTerminalChannelInput { port, kind };
                    channel.validate().map_err(DgSourceError::Physics)?;
                    if let Some((_, cells)) = kind.section() {
                        total_cells = total_cells.checked_add(cells)
                            .ok_or_else(|| bad(line, "complete terminal cell count overflow"))?;
                        if total_cells > allowance {
                            return Err(bad(line, "complete channel state exceeds its transverse-cell allowance").into());
                        }
                    }
                    source.terminal_channels.push(channel);
                }
                (radiation_line, fields) = next_row()?;
            }
            let radiating = source.side_openings.len() + 1 - source.chimneys.iter()
                .filter(|chimney| source.closed_chimney(chimney.port)).count();
            if radiating > MAX_BAFFLED_MIMO_PORTS {
                return Err(bad(radiation_line,
                    "unsealed DG openings exceed the exterior radiation aperture budget").into());
            }
            let mut row = Row { line: radiation_line, fields: &mut fields };
            if row.word()? != "dg_radiation" || row.word()? != "computed-baffled" {
                return Err(bad(radiation_line, "baffled DG outlet requires its computed radiation declaration").into());
            }
            let declaration = row.fields.by_ref().take(ComputedBaffledInput::FIELDS + 1)
                .collect::<Vec<_>>().join(" ");
            let controls = ComputedBaffledInput::parse(&declaration)
                .map_err(|source| DgSourceError::Radiation { line: radiation_line, source })?;
            row.finish()?;
            let (exterior_line, mut fields) = next_row()?;
            let mut row = Row { line: exterior_line, fields: &mut fields };
            if row.word()? != "dg_exterior" {
                return Err(bad(exterior_line, "computed DG radiation requires an explicit exterior gas record").into());
            }
            let exterior = read_gas_record(row.fields, None, sources)
                .map_err(|error| gas_error(exterior_line, error))?;
            row.finish()?;
            if exterior.state().pressure.to_bits() != ambient.pressure.to_bits() {
                return Err(bad(exterior_line, "quiescent DG bore and exterior require equal static pressures").into());
            }
            records.push(exterior.clone());
            source.radiation = Some((controls, exterior));
        }
        source.check_state_budget(ambient)?;
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
        let active = |state: &fs_material::relaxation::RelaxationState|
            state.processes().iter().filter(|p| p.heat_capacity_j_kg_k > 0.0).count();
        // Molecular reservoirs live at the original thermal coordinates:
        // exact-mass nodal arms for dg_tube, axial quadrature arms for dg_bore.
        let thermal_points = self.wall_quadrature
            .map_or(nodes, |points| self.mesh.elements * points);
        let molecular_dofs = self.molecular.as_ref()
            .map_or(0, |state| thermal_points * self.poles * active(state));
        self.check_complete_state_budget(2 * nodes + memory + molecular_dofs)
    }

    // A field's reservoir count belongs to the original DG quadrature owner.
    // Before construction, check_state_budget reserves its acoustic/wall state;
    // prepare checks the exact augmented count here before allocating it.
    fn check_complete_state_budget(&self, bore_dofs: usize) -> Result<(), DgSourceError> {
        let mut total = bore_dofs + self.chimneys.iter()
            .map(|chimney| chimney.state_dofs(self.poles)).sum::<usize>();
        if let Some(law) = &self.molecular_law {
            for chimney in &self.chimneys {
                let medium = fs_material::relaxation::input::medium_from_record(&chimney.gas)
                    .map_err(|source| DgSourceError::Gas { line: self.line, source })?;
                let state = law.evaluate(&medium)
                    .map_err(|source| DgSourceError::Gas { line: self.line, source })?;
                let active = state.processes().iter()
                    .filter(|process| process.heat_capacity_j_kg_k > 0.0).count();
                total += chimney.mesh.elements * (chimney.mesh.order + 1)
                    * self.poles * active;
            }
        }
        if total > self.maximum_dofs {
            return Err(bad(self.line, "complete DG bore and chimney state exceeds the declared DOF allowance"));
        }
        Ok(())
    }

    /// Attach the shared species law to the original internal gas. The same
    /// declaration is evaluated independently for every finite chimney; no
    /// atmosphere, species composition or material coefficient is copied from
    /// another region. Arbitrary admitted bore geometry uses its original
    /// local-cylinder spatial wall quadrature. An axial field retains its
    /// original thermochemical sources and evaluates the law pointwise.
    pub fn with_molecular_relaxation(
        mut self,
        law: fs_material::relaxation::RelaxationAbsorption,
        ambient: &GasRecord,
    ) -> Result<Self, DgSourceError> {
        if self.molecular_law.is_some() {
            return Err(bad(self.line,
                "DG molecular relaxation requires one constitutive declaration"));
        }
        if let Some(field) = &self.field {
            // Validate original source knots now; the core additionally checks
            // every mass and wall quadrature point before making a propagator.
            // Ambient is not a proxy for a varying interior material.
            for &position in field.stations_m() {
                let state = field.state_at(position)
                    .map_err(|source| DgSourceError::Gas { line: self.line, source })?;
                let medium = fs_material::relaxation::RelaxingMedium::from_thermochemical(&state)
                    .map_err(|source| DgSourceError::Gas { line: self.line, source })?;
                law.evaluate(&medium)
                    .map_err(|source| DgSourceError::Gas { line: self.line, source })?;
            }
        } else {
            let medium = fs_material::relaxation::input::medium_from_record(ambient)
                .map_err(|source| DgSourceError::Gas { line: self.line, source })?;
            self.molecular = Some(law.evaluate(&medium)
                .map_err(|source| DgSourceError::Gas { line: self.line, source })?);
        }
        // Validate every region now. A missing partner or unsupported local
        // temperature must refuse before a propagator or radiation bake starts.
        for chimney in &self.chimneys {
            let medium = fs_material::relaxation::input::medium_from_record(&chimney.gas)
                .map_err(|source| DgSourceError::Gas { line: self.line, source })?;
            law.evaluate(&medium)
                .map_err(|source| DgSourceError::Gas { line: self.line, source })?;
        }
        self.molecular_law = Some(law);
        self.check_state_budget(ambient.state())?;
        Ok(self)
    }

    /// Original exterior source, required by composition-sensitive observers.
    pub fn exterior_gas_record(&self) -> Option<&GasRecord> {
        self.radiation.as_ref().map(|(_, record)| record)
    }

    /// Original stationary problem and solution used to bind this bore.
    /// Prescribed fields and homogeneous sources return None.
    pub fn steady_gas(&self) -> Option<&SteadyGasPerformance> {
        self.steady_gas.as_ref()
    }

    fn closed_chimney(&self, port: PortId) -> bool {
        self.chimneys.iter().any(|chimney| chimney.port == port)
            && self.terminal_channels.iter().any(|channel| {
                channel.port == port && matches!(channel.kind, DgTerminalChannelKind::Closed)
            })
    }

    fn prepare_chimneys(&self, gate: &CancelGate) -> Result<Vec<DgBranch>, DgSourceError> {
        self.chimneys.iter().map(|chimney| {
            checkpoint(gate, self.line)?;
            let geometry = Duct { segments: vec![Segment::Cylinder {
                radius: chimney.radius_m, length: chimney.length_m,
            }] };
            let outlet = if self.closed_chimney(chimney.port) {
                DgBoundary::Closed
            } else { DgBoundary::Port(PortId(1)) };
            let duct = DgDuct::new(
                &geometry, chimney.gas.state(), chimney.mesh,
                DgBoundary::Port(PortId(0)), outlet,
            ).and_then(|duct| duct.with_viscothermal(
                self.poles, self.band_hz, self.maximum_relative_error,
            )).map_err(|error| DgSourceError::Physics(AcousticRealizeError::Duct(error)))?;
            let duct = if let Some(law) = &self.molecular_law {
                let medium = fs_material::relaxation::input::medium_from_record(&chimney.gas)
                    .map_err(|source| DgSourceError::Gas { line: self.line, source })?;
                let state = law.evaluate(&medium)
                    .map_err(|source| DgSourceError::Gas { line: self.line, source })?;
                duct.with_molecular_relaxation_controlled(state, &mut || gate.is_requested())
                    .map_err(|error| DgSourceError::Physics(AcousticRealizeError::Duct(error)))?
            } else { duct };
            Ok(DgBranch { port: chimney.port, duct })
        }).collect()
    }

    /// Build once from the original geometry and field, then bind its exact
    /// mean-port load. Tapered walls use the explicit local-cylinder Webster
    /// rung; no conical three-dimensional thermoviscous accuracy is implied.
    pub fn prepare(
        self, ambient: &GasState, dt: f64, gate: &CancelGate,
    ) -> Result<DgPortLoad, DgSourceError> {
        self.prepare_with_material_wall(ambient, dt, gate, None)
    }

    /// Prepare the same DG gas, ports and radiation with an optional original
    /// source-backed wall. All wall inertia, elasticity and positive Maxwell
    /// arms enter the ONE native DG exponential, alongside gas wall memory.
    /// The full bore must satisfy the original uniform-cylinder/thin-wall
    /// admission. Side chimneys retain their explicitly rigid wall model.
    pub fn prepare_with_material_wall(
        self, ambient: &GasState, dt: f64, gate: &CancelGate,
        wall: Option<&crate::bernoulli_aperture::wall::material::LocalMaterialWall>,
    ) -> Result<DgPortLoad, DgSourceError> {
        self.prepare_with_solid_wall(ambient, dt, gate, wall.map(|wall| wall.band_hz()), |owner| {
            if let Some(wall) = wall {
                wall.clone().bind_dg_controlled(owner, &mut || gate.is_requested())
                    .map(|bound| bound.into_duct()).map_err(DgSourceError::Physics)
            } else { Ok(owner) }
        }, Ok)
    }

    /// Bind physical main-bore and optional finite-chimney wall assignments to
    /// the same native gas/port/radiation path. Regions never follow mesh indices;
    /// all original branch gas, wall and hereditary states enter the shared DOF budget.
    pub fn prepare_with_material_wall_layout(
        self, ambient: &GasState, dt: f64, gate: &CancelGate,
        wall: &wall::regions::DgMaterialWallLayout,
    ) -> Result<DgPortLoad, DgSourceError> {
        self.prepare_with_solid_wall(ambient, dt, gate, wall.complete_band_hz(),
            |owner| wall.bind(owner, gate), |branches| wall.bind_chimneys(branches, gate))
    }

    fn prepare_with_solid_wall(
        self,
        ambient: &GasState,
        dt: f64,
        gate: &CancelGate,
        wall_band_hz: Option<[f64; 2]>,
        bind_wall: impl FnOnce(DgDuct) -> Result<DgDuct, DgSourceError>,
        bind_branches: impl FnOnce(Vec<DgBranch>) -> Result<Vec<DgBranch>, DgSourceError>,
    ) -> Result<DgPortLoad, DgSourceError> {
        checkpoint(gate, self.line)?;
        if wall_band_hz.is_some_and(|band|
            band[0] > self.band_hz[0] || band[1] < self.band_hz[1]
        ) {
            return Err(bad(self.line,
                "DG solid wall material band must cover the complete propagation band"));
        }
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
        let ports: Vec<_> = self.side_openings.iter().map(|(port, _)| *port).collect();
        let owner = owner.with_side_ports(&ports)
            .map_err(|error| DgSourceError::Physics(AcousticRealizeError::Duct(error)))?;
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
        let owner = if let (Some(_), Some(law)) = (&self.field, &self.molecular_law) {
            let dofs = owner.molecular_relaxation_field_dofs_controlled(law, &mut || gate.is_requested())
                .map_err(|error| field_error(self.line, error))?;
            self.check_complete_state_budget(dofs)?;
            checkpoint(gate, self.line)?;
            owner.with_molecular_relaxation_field_controlled(law.clone(), &mut || gate.is_requested())
                .map_err(|error| field_error(self.line, error))?
        } else if let Some(molecular) = &self.molecular {
            owner.with_molecular_relaxation_controlled(molecular.clone(), &mut || gate.is_requested())
                .map_err(|error| DgSourceError::Physics(AcousticRealizeError::Duct(error)))?
        } else { owner };
        let owner = bind_wall(owner)?;
        checkpoint(gate, self.line)?;
        let branches = bind_branches(self.prepare_chimneys(gate)?)?;
        let sealed: Vec<_> = branches
            .iter()
            .filter(|branch| branch.duct.outlet() == DgBoundary::Closed)
            .map(|branch| branch.port)
            .collect();
        let load = if let Some((controls, exterior)) = self.radiation {
            let mut openings = Vec::with_capacity(self.side_openings.len() + 1);
            openings.push((PortId(1), CircularOpening {
                center_m: [0.0; 2],
                radius_m: fs_math::det::sqrt(owner.outlet_area_m2() / core::f64::consts::PI),
            }));
            openings.extend(self.side_openings.iter()
                .filter(|(port, _)| !sealed.contains(&port.id))
                .map(|(port, aperture)| (port.id, *aperture)));
            DgPortLoad::new_with_computed_baffled_branches_at_dt(
                owner, dt, self.exponential_tolerance, self.maximum_dofs,
                &exterior, controls, &openings, branches, gate,
            ).map_err(|source| DgSourceError::Radiation { line: self.line, source })?
        } else {
            DgPortLoad::new(owner, dt, self.exponential_tolerance, self.maximum_dofs)
                .map_err(DgSourceError::Physics)?
        };
        checkpoint(gate, self.line)?;
        let terminal_channels = self
            .terminal_channels
            .into_iter()
            .filter(|channel| !sealed.contains(&channel.port))
            .collect();
        load.with_terminal_channels(terminal_channels, self.maximum_terminal_transverse_cells)
            .map_err(DgSourceError::Physics)
    }
}

#[cfg(test)]
mod tests;
