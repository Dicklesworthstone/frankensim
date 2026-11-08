//! Series electrical bindings on an original nonuniform thermal filament.
//!
//! This assembles native material equations only. The caller's prepared heat
//! network retains axial conduction, contacts, surface radiation and storage;
//! its graph/thermal step owns the nonlinear solve and physical clock.

use crate::constitutive_graph::enthalpy::{
    GraphHeatPort, GraphLiquidMassFraction, GraphTemperature, GraphThermalInput, ThermalCoordinate,
};
use crate::constitutive_graph::{MaterialEquation, MaterialUnknown, MaterialValue, PortKey};
use fs_exec::Cx;
use fs_material::conductor::graph::MaterialOhmicNode;
use fs_material::conductor::{ConductorError, OhmicDrive};
use fs_material::graph::{ConstitutiveGraph, GraphError};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, PhaseStateError};
use fs_qty::Dims;
use fs_thermal::enthalpy::{CellGeometry, CellId, DepositionSimulation, MaterialId};
use std::collections::{BTreeMap, BTreeSet};

const CURRENT: Dims = Dims([0, 0, 0, 0, 1, 0]);
const VOLTAGE: Dims = Dims([2, 1, -3, 0, -1, 0]);

/// One original finite cell and its fixed-geometry material conductor.
pub struct FilamentSegment<'a> {
    pub cell: CellId,
    /// The exact caloric object registered with this thermal cell.
    pub curve: &'a EquilibriumEnthalpyPhaseCurve,
    pub source: MaterialOhmicNode,
}

struct CellBinding<'a> {
    cell: CellId,
    material: MaterialId,
    body: u64,
    segment: u64,
    mass_kg: f64,
    geometry: CellGeometry,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    node: String,
}

/// Owned native Ohmic graph and its series equations, with borrowed caloric laws.
///
/// Each segment receives its own endpoint T/f and deposits its own Joule watts
/// exactly once. Different original material cards and caloric charts are valid.
/// An ideal current source fixes the same signed current everywhere; voltage
/// control solves equal currents and the sum of native voltage drops. No current
/// density, conductivity, emissivity, expansion or mixing rule is invented.
pub struct PreparedSeriesFilament<'a> {
    graph: ConstitutiveGraph,
    cells: Vec<CellBinding<'a>>,
    fixed: BTreeMap<PortKey, f64>,
    unknowns: Vec<MaterialUnknown>,
    equations: Vec<MaterialEquation>,
    initial_coordinates: Vec<f64>,
    temperatures: Vec<GraphTemperature>,
    fractions: Vec<GraphLiquidMassFraction<'a>>,
    heat: Vec<GraphHeatPort>,
}

/// Original owner refusals remain distinct from a bad filament correspondence.
#[derive(Debug)]
pub enum FilamentError {
    Invalid(&'static str),
    Electrical(ConductorError),
    Caloric(PhaseStateError),
    Graph(GraphError),
    Cancelled,
}

impl core::fmt::Display for FilamentError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "series filament: {what}"),
            Self::Electrical(error) => write!(f, "series filament electrical source: {error}"),
            Self::Caloric(error) => write!(f, "series filament caloric source: {error}"),
            Self::Graph(error) => write!(f, "series filament material graph: {error}"),
            Self::Cancelled => write!(f, "series filament preparation cancelled"),
        }
    }
}
impl std::error::Error for FilamentError {}

impl<'a> PreparedSeriesFilament<'a> {
    /// Assemble transactionally against existing cells, with a caller work cap.
    ///
    /// Segment length and L*A must match the original centerline and volume.
    /// Complete source support follows every original caloric interval, including
    /// latent motion and crossed material knots. Sourced caloric conditions must
    /// agree with the electrical source; authored charts keep their authority.
    /// Electrical resistance at the initial state supplies only a Newton seed.
    /// Actual resistance and Joule heat are reevaluated at each thermal trial.
    pub fn try_new(
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        segments: Vec<FilamentSegment<'a>>,
        drive: OhmicDrive,
        maximum_segments: usize,
    ) -> Result<Self, FilamentError> {
        cx.checkpoint().map_err(|_| FilamentError::Cancelled)?;
        let drive_value = match drive {
            OhmicDrive::Current(value) | OhmicDrive::Voltage(value) => value,
        };
        if segments.is_empty() || segments.len() > maximum_segments || !drive_value.is_finite() {
            return Err(FilamentError::Invalid(
                "invalid segment count, work cap or drive",
            ));
        }
        let mut prepared = Self {
            graph: ConstitutiveGraph::new(),
            cells: Vec::with_capacity(segments.len()),
            fixed: BTreeMap::new(),
            unknowns: Vec::new(),
            equations: Vec::new(),
            initial_coordinates: Vec::new(),
            temperatures: Vec::with_capacity(segments.len()),
            fractions: Vec::new(),
            heat: Vec::with_capacity(segments.len()),
        };
        let mut seen = BTreeSet::new();
        let mut resistances = Vec::new();
        for segment in segments {
            cx.checkpoint().map_err(|_| FilamentError::Cancelled)?;
            let cell = thermal
                .cells()
                .get(&segment.cell)
                .ok_or(FilamentError::Invalid(
                    "segment cell is absent from the original thermal state",
                ))?;
            let source = segment.source;
            if !seen.insert(segment.cell)
                || !thermal.uses_material_law(cell.material_id, segment.curve)
            {
                return Err(FilamentError::Invalid(
                    "duplicate cell or foreign caloric object",
                ));
            }
            let delta = std::array::from_fn::<_, 3, _>(|i| {
                cell.geometry.end_m[i] - cell.geometry.start_m[i]
            });
            let magnitude = delta.iter().copied().map(f64::abs).fold(0.0_f64, f64::max);
            let length = magnitude
                * delta
                    .iter()
                    .map(|x| (x / magnitude) * (x / magnitude))
                    .sum::<f64>()
                    .sqrt();
            if !same_positive(length, source.length_m())
                || !same_positive(
                    cell.geometry.volume_m3,
                    source.length_m() * source.area_m2(),
                )
            {
                return Err(FilamentError::Invalid(
                    "conductor and original cell geometry differ",
                ));
            }
            // Zero drive admits the full original source path without assigning
            // the terminal voltage independently to every series segment.
            source
                .maximum_positive_joule_specific_enthalpy_slope(
                    OhmicDrive::Current(0.0),
                    segment.curve,
                )
                .map_err(FilamentError::Electrical)?;
            let original = source.source().original().material();
            for receipt in segment.curve.source_receipts() {
                cx.checkpoint().map_err(|_| FilamentError::Cancelled)?;
                for (axis, value) in &receipt.query_point {
                    if receipt.axis_quantities.get(axis) != original.axis_quantities().get(axis)
                        || (axis != source.temperature_axis()
                            && Some(axis.as_str()) != source.liquid_mass_fraction_axis()
                            && original.state_coordinate(axis) != Some(*value))
                    {
                        return Err(FilamentError::Invalid(
                            "caloric and electrical source conditions differ",
                        ));
                    }
                }
            }
            let state = segment
                .curve
                .values_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(FilamentError::Caloric)?;
            let conductor = if source.liquid_mass_fraction_axis().is_some() {
                source.resolve_at_temperature_and_liquid_mass_fraction(
                    state.temperature_k,
                    state.liquid_mass_fraction,
                )
            } else {
                source.resolve_at_temperature_k(state.temperature_k)
            }
            .map_err(FilamentError::Electrical)?;
            resistances.push(conductor.resistance_ohm());
            let node = format!("filament-{}", cell.id.0);
            prepared.temperatures.push(GraphTemperature {
                input: key(&node, "temperature_change"),
                cell: cell.id,
                coordinate: ThermalCoordinate::ChangeFrom(source.reference_temperature_k()),
            });
            if let Some(reference) = source.reference_liquid_mass_fraction() {
                prepared.fractions.push(GraphLiquidMassFraction {
                    input: key(&node, "liquid_mass_fraction_change"),
                    cell: cell.id,
                    curve: segment.curve,
                    reference,
                });
            }
            prepared.heat.push(GraphHeatPort {
                output: key(&node, "heat"),
                cell: cell.id,
            });
            prepared
                .graph
                .add_node(&node, Box::new(source))
                .map_err(FilamentError::Graph)?;
            prepared.cells.push(CellBinding {
                cell: cell.id,
                material: cell.material_id,
                body: cell.body_id,
                segment: cell.segment_id,
                mass_kg: cell.mass_kg,
                geometry: cell.geometry,
                curve: segment.curve,
                node,
            });
        }
        match drive {
            OhmicDrive::Current(current) => {
                for cell in &prepared.cells {
                    prepared.fixed.insert(key(&cell.node, "current"), current);
                }
            }
            OhmicDrive::Voltage(voltage) => {
                let maximum = resistances.iter().copied().fold(0.0_f64, f64::max);
                let scaled_sum = resistances.iter().map(|r| r / maximum).sum::<f64>();
                let mut current = (voltage / maximum) / scaled_sum;
                if !current.is_finite() {
                    current = (voltage / scaled_sum) / maximum;
                }
                if !current.is_finite() || (voltage != 0.0 && current == 0.0) {
                    return Err(FilamentError::Invalid(
                        "unrepresentable initial series current",
                    ));
                }
                let current_scale = if current == 0.0 { 1.0 } else { current.abs() };
                let first = key(&prepared.cells[0].node, "current");
                let mut drops = Vec::with_capacity(prepared.cells.len());
                for (index, cell) in prepared.cells.iter().enumerate() {
                    cx.checkpoint().map_err(|_| FilamentError::Cancelled)?;
                    let input = key(&cell.node, "current");
                    prepared.unknowns.push(MaterialUnknown {
                        input: input.clone(),
                        dims: CURRENT,
                        scale_si: current_scale,
                    });
                    prepared.initial_coordinates.push(current / current_scale);
                    if index != 0 {
                        prepared.equations.push(MaterialEquation {
                            terms: vec![
                                (1.0, MaterialValue::Input(input)),
                                (-1.0, MaterialValue::Input(first.clone())),
                            ],
                            dims: CURRENT,
                            target_si: 0.0,
                            scale_si: current_scale,
                        });
                    }
                    drops.push((1.0, MaterialValue::Output(key(&cell.node, "voltage"))));
                }
                prepared.equations.push(MaterialEquation {
                    terms: drops,
                    dims: VOLTAGE,
                    target_si: voltage,
                    scale_si: if voltage == 0.0 { 1.0 } else { voltage.abs() },
                });
            }
        }
        cx.checkpoint().map_err(|_| FilamentError::Cancelled)?;
        Ok(prepared)
    }

    /// Native graph input for a step on these same original finite cells.
    ///
    /// Revalidating the caloric object, mass and geometry prevents reusing this
    /// binding on a changed specimen. Enthalpy and thermal boundaries remain
    /// owned by the caller's evolving thermal state and prepared heat network.
    pub fn thermal_input<'s>(
        &'s self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<GraphThermalInput<'s>, FilamentError> {
        for binding in &self.cells {
            cx.checkpoint().map_err(|_| FilamentError::Cancelled)?;
            let cell = thermal
                .cells()
                .get(&binding.cell)
                .ok_or(FilamentError::Invalid("prepared filament cell is absent"))?;
            if cell.material_id != binding.material
                || cell.body_id != binding.body
                || cell.segment_id != binding.segment
                || cell.mass_kg != binding.mass_kg
                || cell.geometry != binding.geometry
                || !thermal.uses_material_law(cell.material_id, binding.curve)
            {
                return Err(FilamentError::Invalid(
                    "prepared filament cell binding changed",
                ));
            }
        }
        Ok(GraphThermalInput {
            graph: &self.graph,
            version: self.graph.state_schema().version(),
            state: &[],
            fixed: &self.fixed,
            unknowns: &self.unknowns,
            equations: &self.equations,
            initial_coordinates: &self.initial_coordinates,
            temperatures: &self.temperatures,
            liquid_mass_fractions: &self.fractions,
            durations: &[],
            heat: &self.heat,
        })
    }

    /// Stable native node identifier for inspecting per-cell currents and heat.
    #[must_use]
    pub fn node_id(&self, cell: CellId) -> Option<&str> {
        self.cells
            .iter()
            .find(|binding| binding.cell == cell)
            .map(|binding| binding.node.as_str())
    }
}

fn key(node: &str, port: &str) -> PortKey {
    (node.to_owned(), port.to_owned())
}

fn same_positive(a: f64, b: f64) -> bool {
    a.is_finite()
        && b.is_finite()
        && a > 0.0
        && b > 0.0
        && (a - b).abs() <= 64.0 * f64::EPSILON * a.max(b)
}
