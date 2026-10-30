//! Material-driven convection on the shared finite-cell enthalpy state.
//!
//! Prescribed-reservoir boundaries and reciprocal finite-fluid links use the
//! same material/correlation owner. One original transport snapshot supplies
//! Re or Ra, Pr and h per binding.
//! The original thermal owner supplies wall temperature and currently exposed
//! area. Implicit solves evaluate every thermal trial; explicit observers
//! re-evaluate at each native substep. Neither freezes the initial heat source.
//! Borrowed-snapshot boundaries freeze fluid properties. The source-backed
//! adapter re-evaluates transport and expansion at declared reservoir or film
//! temperature inside the same residual. Boundary and finite-link APIs retain
//! their own storage ownership. Geometry stays prescribed; neither resolves flow.

pub mod source;

use std::collections::{BTreeMap, BTreeSet};
use fs_conduction::transient::deposition::{CellId, DepositionError, DepositionSimulation};
use super::*;

/// Explicit reservoir and forced-flow geometry for one whole exposed cell.
/// Several cells may borrow the SAME gas mixture or liquid source. The wall's
/// temperature is solved; the fluid source's property state is authored.
#[derive(Clone, Copy)]
pub struct ConvectiveCell<'a> {
    pub cell: CellId,
    pub fluid: &'a dyn FluidTransport,
    pub reference_temperature: Temperature,
    pub characteristic_length: Length,
    pub correlation: CorrelationId,
    pub speed_m_s: f64,
    /// Geometry settings; Dittus-Boelter heating/cooling direction is derived
    /// from the actual wall/reservoir temperatures at each trial.
    pub settings: ForcedConvectionSettings,
}

/// A boundary law evaluated from the same material and current thermal cell.
/// Each exposed cell has exactly one law; forced and natural coefficients are
/// not summed into an unsupported mixed-convection correlation.
#[derive(Clone, Copy, Debug)]
pub enum CellConvectionLaw {
    /// Prescribed relative flow through an existing forced-convection card.
    Forced {
        correlation: CorrelationId,
        speed_m_s: f64,
        settings: ForcedConvectionSettings,
    },
    /// Churchill-Chu isothermal vertical plate. Gravity is a nonnegative
    /// magnitude [m/s²]; volumetric expansion [1/K] is explicitly supplied.
    /// Negative expansion is admitted by the original vertical-plate owner.
    /// No liquid expansion is inferred from a gas temperature or material name.
    VerticalPlate {
        gravity_m_s2: f64,
        expansion_per_k: f64,
    },
}

/// One material source, reservoir and law on the whole currently exposed cell.
/// The property state is explicitly frozen. The caller must retain its source
/// and justify the selected card's property-evaluation temperature convention.
#[derive(Clone, Copy)]
pub struct MaterialConvectiveCell<'a> {
    pub cell: CellId,
    pub fluid: &'a dyn FluidTransport,
    pub reference_temperature: Temperature,
    pub characteristic_length: Length,
    pub law: CellConvectionLaw,
}

impl<'a> From<ConvectiveCell<'a>> for MaterialConvectiveCell<'a> {
    fn from(binding: ConvectiveCell<'a>) -> Self {
        Self {
            cell: binding.cell,
            fluid: binding.fluid,
            reference_temperature: binding.reference_temperature,
            characteristic_length: binding.characteristic_length,
            law: CellConvectionLaw::Forced {
                correlation: binding.correlation,
                speed_m_s: binding.speed_m_s,
                settings: binding.settings,
            },
        }
    }
}

/// One supported correlation evaluation, without an invented equilibrium limit.
#[derive(Debug)]
pub struct CellConvectionEvaluation {
    pub nusselt: NusseltEvaluation,
    pub coefficient: Evidence<HeatTransferCoefficient>,
}

/// Original dimensional/correlation result and the exact thermal destination.
#[derive(Debug)]
pub struct CellConvectionSample {
    pub cell: CellId,
    pub wall_temperature_k: f64,
    pub reference_temperature_k: f64,
    pub exposed_area_m2: f64,
    pub transport: FluidTransportState,
    /// Absent only at exact natural-convection thermal equilibrium: the power
    /// is zero, but a correlation at unsupported Ra=0 has not been evaluated.
    pub evaluation: Option<CellConvectionEvaluation>,
    /// Signed external power into the finite body. Cooling is negative.
    pub heat_into_cell_w: f64,
}

#[derive(Debug)]
pub enum ConvectionSourceError {
    Invalid(&'static str),
    Material(MaterialConvectionError),
    Thermal(DepositionError),
    Cancelled,
}
impl From<DepositionError> for ConvectionSourceError {
    fn from(error: DepositionError) -> Self {
        if error == DepositionError::Cancelled { Self::Cancelled } else { Self::Thermal(error) }
    }
}
impl fmt::Display for ConvectionSourceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "material convection source: {self:?}")
    }
}
impl std::error::Error for ConvectionSourceError {}
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), ConvectionSourceError> {
    if cancelled() { Err(ConvectionSourceError::Cancelled) } else { Ok(()) }
}
fn finite(value: f64) -> Result<f64, ConvectionSourceError> {
    if value.is_finite() { Ok(value) } else {
        Err(ConvectionSourceError::Invalid("convective power is unrepresentable"))
    }
}

/// Endpoint powers and original source evaluations, with no state advanced.
#[derive(Debug)]
pub struct CellConvectionSources {
    sources_w: BTreeMap<CellId, f64>,
    samples: Vec<CellConvectionSample>,
    total_w: f64,
}
impl CellConvectionSources {
    pub fn sources_w(&self) -> &BTreeMap<CellId, f64> { &self.sources_w }
    pub fn samples(&self) -> &[CellConvectionSample] { &self.samples }
    /// External heat into the simulated bodies. The prescribed reservoir receives
    /// the opposite power; it has no invented finite mass or temperature history.
    pub const fn total_heat_into_cells_w(&self) -> f64 { self.total_w }

    /// Add once to a trial's additional source map, retaining independent heaters,
    /// contact/Joule heat and radiation. Validate every sum before publication.
    /// The shared thermal owner separately retains authored boundary sources.
    pub fn add_to(&self, destination: &mut BTreeMap<CellId, f64>,
        mut cancelled: impl FnMut() -> bool) -> Result<(), ConvectionSourceError> {
        check(&mut cancelled)?;
        let mut staged = Vec::with_capacity(self.sources_w.len());
        for (&cell, &power) in &self.sources_w {
            check(&mut cancelled)?;
            let old = destination.get(&cell).copied().unwrap_or(0.0);
            let next = finite(finite(old)? + power)?;
            if power != 0.0 && next == old {
                return Err(ConvectionSourceError::Invalid("convective power disappears in source addition"));
            }
            staged.push((cell, next));
        }
        check(&mut cancelled)?;
        for (cell, power) in staged { destination.insert(cell, power); }
        Ok(())
    }
}

/// Resolve forced convection at the ACTUAL proposed caloric endpoint.
///
/// The thermal owner's exposed area already excludes active contact/reservoir
/// coverage. Never add a second ambient Robin term: a bound cell with nonzero
/// authored `htc_w_m2_k` refuses. Independent authored heaters remain valid.
/// Missing/duplicate cells, invalid reservoirs, correlation-domain failure or
/// cancellation return no partial map. `maximum_cells` bounds material queries.
/// Source order is stable CellId order, not input order.
///
/// A caller must retain the borrowed fluid source alongside the report to retain
/// its material receipts. A transport tuple is not independent material evidence.
/// No molecular identity, fluid expansion, correlation or film temperature is
/// inferred from the material name. Unsupported zero-flow groups still refuse.
pub fn evaluate_cell_convection(
    thermal: &DepositionSimulation<'_>, bindings: &[ConvectiveCell<'_>],
    maximum_cells: usize, cancelled: impl FnMut() -> bool,
) -> Result<CellConvectionSources, ConvectionSourceError> {
    evaluate_bound_cells(thermal, bindings.iter().copied().map(MaterialConvectiveCell::from),
        maximum_cells, cancelled)
}

/// Evaluate forced and natural boundaries through one thermal ownership loop.
///
/// Use this inside the original thermal residual, where `thermal` is the actual
/// proposed enthalpy endpoint. A natural boundary's Rayleigh number and h then
/// follow its current wall/reservoir temperature difference. Material properties
/// and expansion remain frozen at the explicitly borrowed source state.
/// Exact natural-convection thermal equilibrium returns zero power with no
/// correlation evaluation, after admitting the original fluid and physical
/// inputs. Nonzero temperature differences still require an admitted Rayleigh
/// number; zero gravity/expansion do not invent an out-of-domain coefficient.
///
/// Exposed-area/contact accounting, independent source addition, stable ordering,
/// work limits and transactional refusals are identical to
/// `evaluate_cell_convection`. A repeated cell across different laws also
/// refuses: this API does not compose a mixed-convection law on one surface.
pub fn evaluate_material_cell_convection(
    thermal: &DepositionSimulation<'_>, bindings: &[MaterialConvectiveCell<'_>],
    maximum_cells: usize, cancelled: impl FnMut() -> bool,
) -> Result<CellConvectionSources, ConvectionSourceError> {
    evaluate_bound_cells(thermal, bindings.iter().copied(), maximum_cells, cancelled)
}

fn evaluate_bound_cells<'a>(
    thermal: &DepositionSimulation<'_>,
    bindings: impl ExactSizeIterator<Item = MaterialConvectiveCell<'a>>,
    maximum_cells: usize, mut cancelled: impl FnMut() -> bool,
) -> Result<CellConvectionSources, ConvectionSourceError> {
    check(&mut cancelled)?;
    if bindings.len() > maximum_cells || bindings.len() > thermal.cells().len() {
        return Err(ConvectionSourceError::Invalid("convective cells exceed their explicit work bound"));
    }
    let mut ids = BTreeSet::new();
    let mut ordered = Vec::with_capacity(bindings.len());
    for binding in bindings {
        check(&mut cancelled)?;
        if !ids.insert(binding.cell) {
            return Err(ConvectionSourceError::Invalid("duplicate convective cell"));
        }
        let cell = thermal.cells().get(&binding.cell)
            .ok_or(ConvectionSourceError::Invalid("inactive convective cell"))?;
        if cell.boundary.htc_w_m2_k != 0.0 {
            return Err(ConvectionSourceError::Invalid("material convection would duplicate the authored ambient coefficient"));
        }
        let reference = binding.reference_temperature.value();
        if !reference.is_finite() || reference <= 0.0 {
            return Err(ConvectionSourceError::Invalid("reservoir temperature must be positive finite kelvin"));
        }
        ordered.push(binding);
    }
    ordered.sort_by_key(|binding| binding.cell);
    let mut samples = Vec::with_capacity(ordered.len());
    let mut sources_w = BTreeMap::new();
    let mut total_w = 0.0;
    for binding in ordered {
        check(&mut cancelled)?;
        let wall = thermal.temperature_k(binding.cell)?;
        let area = thermal.exposed_area_m2(binding.cell)?;
        // Query the original source exactly once. Both group evaluation and
        // dimensional conversion below consume this same admitted snapshot.
        let transport = binding.fluid.transport_state()
            .map_err(|error| ConvectionSourceError::Material(error.into()))?;
        let evaluation = evaluate_convection_state(&transport, wall,
            binding.reference_temperature, binding.characteristic_length, binding.law)
            .map_err(ConvectionSourceError::Material)?;
        let delta = finite(binding.reference_temperature.value() - wall)?;
        let power = match &evaluation {
            None => 0.0,
            Some(value) => {
                let conductance = finite(value.coefficient.value.value() * area)?;
                let power = finite(conductance * delta)?;
                if (area > 0.0 && conductance == 0.0)
                    || (conductance != 0.0 && delta != 0.0 && power == 0.0) {
                    return Err(ConvectionSourceError::Invalid("convective transfer underflowed to zero"));
                }
                power
            }
        };
        total_w = finite(total_w + power)?;
        sources_w.insert(binding.cell, power);
        samples.push(CellConvectionSample { cell: binding.cell, wall_temperature_k: wall,
            reference_temperature_k: binding.reference_temperature.value(), exposed_area_m2: area,
            transport, evaluation, heat_into_cell_w: power });
    }
    check(&mut cancelled)?;
    Ok(CellConvectionSources { sources_w, samples, total_w })
}

// All finite-cell boundaries and internal fluid links use this same correlation
// admission, actual thermal direction and exact-equilibrium convention.
fn evaluate_convection_state(
    transport: &FluidTransportState, wall: f64, reference: Temperature,
    length: Length, law: CellConvectionLaw,
) -> Result<Option<CellConvectionEvaluation>, MaterialConvectionError> {
    let result = match law {
        CellConvectionLaw::Forced { correlation, speed_m_s, mut settings } => {
            if wall > reference.value() {
                settings.direction = ThermalDirection::HeatingFluid;
            } else if wall < reference.value() {
                settings.direction = ThermalDirection::CoolingFluid;
            }
            evaluate_forced_from_material(correlation, transport,
                speed_m_s, length, settings).map(Some)
        }
        CellConvectionLaw::VerticalPlate { gravity_m_s2, expansion_per_k } => {
            if wall == reference.value() {
                transport.rayleigh(gravity_m_s2, expansion_per_k, 0.0, length.value())
                    .map_err(MaterialConvectionError::from).map(|_| None)
            } else {
                evaluate_vertical_plate_from_material(transport, gravity_m_s2,
                    expansion_per_k, Temperature::new(wall), reference, length).map(Some)
            }
        }
    }?;
    Ok(result.map(|value| {
        let (_, nusselt, coefficient) = value.into_parts();
        CellConvectionEvaluation { nusselt, coefficient }
    }))
}
