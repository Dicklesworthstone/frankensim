//! Material thermocouples and Peltier legs on the existing circuit/enthalpy solve.
//!
//! A leg is an original material resistor in series with one original voltage
//! source branch. Its two junctions and original bulk cell receive reciprocal
//! Kelvin heat from source-integrated thermal edges; the bulk receives the
//! variable-Seebeck interior transfer as well as Joule heat exactly once. Original thermal
//! links carry Fourier conduction using the SAME material/geometry bindings.
//! All material, electrical and thermal candidates use the same proposed state.

/// Reciprocal thermoelectric, field-motion and finite-thermal evolution.
pub mod motion;
/// Exact caloric and reversible Kelvin entropy on the same accepted interval.
pub mod entropy;

use super::{CircuitEnthalpyError, CircuitEnthalpyInput,
    CircuitEnthalpyReport, CircuitEnthalpySolve, CircuitThermalCoefficients,
    DrivenCoefficients, ThermalDriver, coupled_driver_candidate_entropy_option};
use super::entropy::CircuitEntropyAudit;
use super::material::{MaterialNetwork, MaterialNetworkError, MaterialNetworkSample};
use crate::material_circuit::{MaterialCircuit, MaterialCircuitCandidate,
    MaterialCircuitError, MaterialResistor};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermoelectric::{ABSOLUTE_SEEBECK_PROPERTY, ResolvedSeebeckSpan,
    ThermoelectricError, ThermoelectricTransfer};
use fs_phs::circuit::Branch;
use fs_phs::StepRecord;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::cell::RefCell;
use std::collections::BTreeSet;

/// A source-bound leg; its material and geometry come from `MaterialNetwork`.
pub struct ThermoelectricLeg {
    /// Original source branch (a, middle, VoltageSource), not its port index.
    pub source_branch: usize,
    /// Original material resistor (middle, b, Resistor), in matching orientation.
    pub resistor_branch: usize,
    /// Thermal junction cells a and b. Both are distinct from the bulk heat cell.
    pub terminals: [CellId; 2],
    /// Existing bulk-to-junction thermal links in the SAME a,b order. Their
    /// areas must equal the conductor area and bulk distances must sum to L.
    pub thermal_links: [u64; 2],
    /// Absolute Seebeck source choice on the resistor's original card.
    /// Constant scalars and supported piecewise-linear temperature paths,
    /// including fixed-context slices through two-coordinate sources, are admitted.
    pub selection: MaterialPropertySelection,
}

/// Original constitutive refusals survive all numerical trial callbacks.
#[derive(Debug)]
pub enum ThermoelectricNetworkError {
    /// Original common caloric/electrical/conductive material owner.
    Network(MaterialNetworkError),
    /// Original reversible thermoelectric source owner.
    Seebeck(ThermoelectricError),
    /// Unsupported leg topology, phase or geometry.
    Invalid(&'static str),
}
impl core::fmt::Display for ThermoelectricNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermoelectric network: {self:?}")
    }
}
impl std::error::Error for ThermoelectricNetworkError {}

/// Exact source, orientation and solved current for one accepted leg.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermoelectricLegSample {
    /// Original graph source branch.
    pub source_branch: usize,
    /// Original graph resistor branch.
    pub resistor_branch: usize,
    /// Junction destination order used for signed reversible heat.
    pub terminals: [CellId; 2],
    /// Original resistor/caloric cell receiving the interior reversible heat.
    pub bulk: CellId,
    /// Original endpoint and bulk receipts and integrated thermal-edge law.
    pub material: ResolvedSeebeckSpan,
    /// Transfer at the ACTUAL Kirchhoff midpoint source current.
    pub transfer: ThermoelectricTransfer,
}
/// Common material payload plus reversible transfers. Joule heat remains in
/// the native circuit report, not in a second thermoelectric heat model.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermoelectricSample {
    /// Original density/caloric/electrical/conduction material evaluations.
    pub network: MaterialNetworkSample,
    /// Leg records in original input order.
    pub legs: Vec<ThermoelectricLegSample>,
    /// Internal Seebeck-port work supplied to electrical storage/loss [J].
    pub reversible_electrical_work_j: f64,
    /// Internal signed reversible heat in junction AND bulk cells [J].
    pub reversible_heat_j: f64,
}

/// Both native histories, with internal cross-domain transfer removed from
/// external supply accounting. Numerical convergence is not a time-error bound.
pub struct ThermoelectricReport {
    /// Existing circuit and thermal reports; its `non_electrical_heat_j` includes
    /// reversible junction/bulk transfer and its native port supply includes Seebeck work.
    pub coupled: CircuitEnthalpyReport<ThermoelectricSample>,
    /// Work supplied through the remaining authored electrical ports [J].
    pub external_electrical_work_j: f64,
    /// Authored thermal boundary/source heat, excluding Joule AND Kelvin transfer [J].
    pub external_heat_j: f64,
    /// Total electrical + thermal storage minus external supply [J].
    pub total_energy_residual_j: f64,
}
/// Atomic electrical/thermal endpoint; neither accepted input is changed.
pub struct ThermoelectricCandidate<'m> {
    thermal: DepositionSimulation<'m>,
    circuit: MaterialCircuit,
    report: ThermoelectricReport,
}
impl<'m> ThermoelectricCandidate<'m> {
    /// Proposed thermal history with original event semantics.
    pub const fn thermal(&self) -> &DepositionSimulation<'m> { &self.thermal }
    /// Proposed conserved charges and flux linkages.
    pub const fn circuit(&self) -> &MaterialCircuit { &self.circuit }
    /// Complete sources and independent energy account.
    pub const fn report(&self) -> &ThermoelectricReport { &self.report }
    /// Publish both histories together with their original native evidence.
    pub fn into_parts(self) -> (DepositionSimulation<'m>, MaterialCircuit, ThermoelectricReport) {
        (self.thermal, self.circuit, self.report)
    }
}

type Error = CircuitEnthalpyError<ThermoelectricNetworkError>;
fn invalid(message: &'static str) -> Error { CircuitEnthalpyError::Invalid(message) }
fn network(error: MaterialNetworkError) -> Error {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ThermoelectricNetworkError::Network(other)),
    }
}
fn finite(value: f64) -> Result<f64, Error> {
    if value.is_finite() { Ok(value) } else { Err(invalid("unrepresentable thermoelectric work")) }
}
fn same_geometry(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && a > 0.0 && b > 0.0
        && (a - b).abs() <= 64.0 * f64::EPSILON * a.max(b)
}

// Return the original source port and material-cell index for every leg.
fn admit(
    thermal: &DepositionSimulation<'_>, material: &MaterialNetwork<'_>,
    circuit: &MaterialCircuit, ports: &[f64], legs: &[ThermoelectricLeg],
) -> Result<Vec<(usize, usize)>, Error> {
    let graph = circuit.graph();
    if legs.is_empty() || legs.len() > material.resistors.len()
        || material.resistors.len() != circuit.resistors().len()
        || !circuit.contacts().is_empty() {
        return Err(invalid("bounded thermoelectric legs need bulk-only material resistor branches"));
    }
    let mut sources = BTreeSet::new();
    let mut resistors = BTreeSet::new();
    let mut routing = Vec::with_capacity(legs.len());
    for leg in legs {
        let &(a, middle, source) = graph.branches.get(leg.source_branch)
            .ok_or(invalid("thermoelectric source branch is absent"))?;
        let &(r_middle, _, resistor) = graph.branches.get(leg.resistor_branch)
            .ok_or(invalid("thermoelectric resistor branch is absent"))?;
        let Branch::VoltageSource { port } = source else {
            return Err(invalid("Seebeck branch must be an original voltage source"));
        };
        if !matches!(resistor, Branch::Resistor { .. }) || middle != r_middle || a == middle
            || !sources.insert(leg.source_branch) || !resistors.insert(leg.resistor_branch)
            || graph.branches.iter().filter(|(x, y, _)| *x == middle || *y == middle).count() != 2
            || graph.transformers.iter().any(|t| [t.primary.0, t.primary.1,
                t.secondary.0, t.secondary.1].contains(&middle)) {
            return Err(invalid("each leg needs one oriented resistor in exclusive series with its Seebeck source"));
        }
        if ports.get(port) != Some(&0.0) || graph.branches.iter().filter(|(_, _, b)|
            matches!(b, Branch::VoltageSource { port: p } | Branch::CurrentSource { port: p } if *p == port)).count() != 1 {
            return Err(invalid("each internal Seebeck port must be unique with a zero external-input placeholder"));
        }
        let bulk = material.resistors.iter().find(|r| r.destination.branch == leg.resistor_branch)
            .ok_or(invalid("thermoelectric resistor has no shared material binding"))?;
        let core = bulk.destination.cell;
        let cell_index = material.cells.iter().position(|c| c.cell == core)
            .ok_or(invalid("thermoelectric bulk cell is absent"))?;
        if leg.terminals[0] == leg.terminals[1] || leg.terminals.contains(&core)
            || leg.thermal_links[0] == leg.thermal_links[1] {
            return Err(invalid("distinct bulk and junction cells with two conductive links are required"));
        }
        let mut bulk_length = 0.0;
        for (&terminal, &id) in leg.terminals.iter().zip(&leg.thermal_links) {
            let link = thermal.links().get(&id).ok_or(invalid("thermoelectric thermal link is absent"))?;
            let binding = material.links.iter().find(|l| l.link == id)
                .ok_or(invalid("thermoelectric thermal link needs the shared conductivity binding"))?;
            let distance = if link.a == core && link.b == terminal { binding.distance_a_m }
                else if link.b == core && link.a == terminal { binding.distance_b_m }
                else { return Err(invalid("leg conductive links must connect its bulk and ordered junctions")); };
            if !same_geometry(link.area_m2, bulk.area_m2) {
                return Err(invalid("electrical and thermal leg cross-sections disagree"));
            }
            bulk_length += distance;
        }
        if !same_geometry(bulk_length, bulk.length_m) {
            return Err(invalid("electrical length and both thermal bulk distances disagree"));
        }
        routing.push((port, cell_index));
    }
    Ok(routing)
}

fn span(
    material: &MaterialNetwork<'_>, state: &DepositionSimulation<'_>,
    coefficients: &CircuitThermalCoefficients<MaterialNetworkSample>,
    leg: &ThermoelectricLeg, cell_index: usize, selection: MaterialPropertySelection,
) -> Result<ResolvedSeebeckSpan, Error> {
    let binding = &material.cells[cell_index];
    let phase = coefficients.source.cells.iter().find(|s| s.cell == binding.cell)
        .ok_or(invalid("thermoelectric material phase sample is absent"))?.phase;
    let temperatures = [state.temperature_k(leg.terminals[0])?, state.temperature_k(leg.terminals[1])?];
    let knots = binding.curve.knots();
    if knots.iter().any(|k| k.liquid_mass_fraction != 0.0)
        || temperatures.iter().any(|t| *t < knots[0].temperature_k
            || *t > knots[knots.len() - 1].temperature_k) {
        return Err(invalid("thermoelectric solid leg requires one supported single-phase caloric span"));
    }
    let context = binding.point.at(phase).map_err(|e| network(MaterialNetworkError::Material(e)))?;
    ResolvedSeebeckSpan::resolve(binding.point.material, &context, binding.point.temperature_axis,
        temperatures, selection).map_err(|e| CircuitEnthalpyError::Source(ThermoelectricNetworkError::Seebeck(e)))
}

/// Advance a thermocouple/Peltier network through the original coupled solver.
///
/// Material cards supply rho_e or sigma, k, caloric storage and ABSOLUTE S(T).
/// Transport and Seebeck claim choices are pinned before nonlinear trials.
/// S may be a validity-wide scalar or a supported piecewise-linear temperature
/// path through curve or surface data at fixed external conditions. Its
/// integrals on junction-to-bulk thermal edges drive the SAME native
/// circuit current and reversible heat into all three original caloric cells.
/// Each edge conserves electrical/thermal energy and transfers zero net entropy.
/// For variable S, the bulk receives an interior Thomson contribution; the edge
/// averages also include portions in the terminal heat. These are finite-volume
/// nodal transfers, not exact continuum pointwise Peltier fluxes or a whole-leg
/// Thomson source lumped into the bulk. Refine the physical thermal network to
/// resolve that spatial distribution. Magnetic fields, electrical contact jumps
/// and phase transitions remain outside this solid-leg model.
/// Internal source ports require zero input placeholders.
///
/// Multiple legs can meet at shared junction cells: their signed Peltier terms
/// then produce the physical coefficient DIFFERENCE without a junction preset.
/// Fourier conduction, lumped bulk Joule deposition, source ports and latent
/// storage in OTHER thermal cells remain with the original owners. Spatial
/// lumping and backward-Euler endpoint coefficients require refinement; this is
/// not a resolved semiconductor, electronic heat-capacity or AC/tensor model.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'m>(
    material: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
    circuit: &MaterialCircuit, ports: &[f64], legs: &[ThermoelectricLeg],
    electrical_bases: &[(usize, ElectricalPropertyBasis)], controls: CircuitEnthalpySolve,
    guess_j: &[f64], mut cancelled: impl FnMut() -> bool,
) -> Result<ThermoelectricCandidate<'m>, Error> {
    let (thermal, circuit, coupled) = coupled_candidate(material, thermal, target_s,
        circuit, ports, legs, electrical_bases, controls, guess_j,
        |sources, inputs, dt, cancelled| circuit.step_candidate(sources, inputs, dt,
            controls.circuit, cancelled), &mut cancelled)?;
    finish_candidate(thermal, circuit, coupled, controls, cancelled)
}

fn finish_candidate<'m>(
    thermal: DepositionSimulation<'m>, circuit: MaterialCircuit,
    coupled: CircuitEnthalpyReport<ThermoelectricSample>, controls: CircuitEnthalpySolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ThermoelectricCandidate<'m>, Error> {
    let external_electrical_work_j = finite(coupled.electrical.native.supplied
        - coupled.source.reversible_electrical_work_j)?;
    let external_heat_j = finite(coupled.non_electrical_heat_j - coupled.source.reversible_heat_j)?;
    let total_energy_residual_j = finite(coupled.electrical.stored_after_j
        - coupled.electrical.stored_before_j + coupled.thermal.stored_energy_change_j
        - external_electrical_work_j - external_heat_j)?;
    if total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
        return Err(invalid("combined thermoelectric and thermal storage do not close external supply"));
    }
    if cancelled() { return Err(CircuitEnthalpyError::Cancelled); }
    Ok(ThermoelectricCandidate { thermal, circuit, report: ThermoelectricReport {
        coupled, external_electrical_work_j, external_heat_j, total_energy_residual_j } })
}

// Both native candidates expose the SAME Kirchhoff source-port ordering. The
// moving driver merely appends its mechanical port; it cannot author Seebeck
// current, loss, or a second thermal/electrical update.
trait ThermoelectricDriver: ThermalDriver {
    fn native(&self) -> &StepRecord;
}
impl ThermoelectricDriver for MaterialCircuitCandidate {
    fn native(&self) -> &StepRecord { &self.step().native }
}

#[allow(clippy::too_many_arguments)]
fn coupled_candidate<'m, D: ThermoelectricDriver>(
    material: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
    circuit: &MaterialCircuit, ports: &[f64], legs: &[ThermoelectricLeg],
    electrical_bases: &[(usize, ElectricalPropertyBasis)], controls: CircuitEnthalpySolve,
    guess_j: &[f64],
    advance: impl Fn(Vec<MaterialResistor>, &[f64], f64, &mut dyn FnMut() -> bool)
        -> Result<D, MaterialCircuitError>,
    cancelled: impl FnMut() -> bool,
) -> Result<(DepositionSimulation<'m>, D::Model,
    CircuitEnthalpyReport<ThermoelectricSample, D::Step>), Error> {
    coupled_candidate_entropy_option(material, thermal, target_s, circuit, ports,
        legs, electrical_bases, controls, guess_j, advance, false, cancelled)
        .map(|(thermal, model, report, _)| (thermal, model, report))
}

#[allow(clippy::too_many_arguments)]
fn coupled_candidate_entropy_option<'m, D: ThermoelectricDriver>(
    material: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
    circuit: &MaterialCircuit, ports: &[f64], legs: &[ThermoelectricLeg],
    electrical_bases: &[(usize, ElectricalPropertyBasis)], controls: CircuitEnthalpySolve,
    guess_j: &[f64],
    advance: impl Fn(Vec<MaterialResistor>, &[f64], f64, &mut dyn FnMut() -> bool)
        -> Result<D, MaterialCircuitError>,
    with_entropy: bool,
    cancelled: impl FnMut() -> bool,
) -> Result<(DepositionSimulation<'m>, D::Model,
    CircuitEnthalpyReport<ThermoelectricSample, D::Step>, Option<CircuitEntropyAudit>), Error> {
    let cancel = RefCell::new(cancelled);
    let check = || -> Result<(), Error> {
        if (*cancel.borrow_mut())() { Err(CircuitEnthalpyError::Cancelled) } else { Ok(()) }
    };
    check()?;
    let routing = admit(thermal, material, circuit, ports, legs)?;
    let prepared = material.prepare_with_bases(thermal, electrical_bases,
        || (*cancel.borrow_mut())()).map_err(network)?;
    let initial = prepared.resolve(thermal, || (*cancel.borrow_mut())()).map_err(network)?;
    let mut pins = Vec::with_capacity(legs.len());
    for (leg, &(_, cell_index)) in legs.iter().zip(&routing) {
        check()?;
        let source = span(material, thermal, &initial, leg, cell_index, leg.selection.clone())?;
        let selected = source.material()[0].property(ABSOLUTE_SEEBECK_PROPERTY)
            .ok_or(invalid("admitted Seebeck claim is absent"))?.answer().receipt.selected;
        pins.push(MaterialPropertySelection::PinnedByProperty(vec![
            (ABSOLUTE_SEEBECK_PROPERTY.to_owned(), selected)]));
    }
    let dt = target_s - thermal.time_s();
    let destinations = material.resistors.iter().map(|r| r.destination).collect::<Vec<_>>();
    coupled_driver_candidate_entropy_option(
        thermal, target_s, CircuitEnthalpyInput { circuit, ports, destinations: &destinations },
        &[], controls, guess_j, |state, _remaining| {
            let mut coefficients = prepared.resolve(state, || (*cancel.borrow_mut())()).map_err(network)?;
            // This material owner supplies only transport; authored independent
            // powers remain on original cell boundaries. Future source additions
            // must provide their classification rather than silently disappear.
            if with_entropy && !coefficients.sources_w.is_empty() {
                return Err(invalid("thermoelectric entropy requires explicit independent source semantics"));
            }
            let mut inputs = ports.to_vec();
            let mut spans = Vec::with_capacity(legs.len());
            for ((leg, &(port, cell_index)), pin) in legs.iter().zip(&routing).zip(&pins) {
                check()?;
                let source = span(material, state, &coefficients, leg, cell_index, pin.clone())?;
                inputs[port] = source.transfer(0.0).map_err(|e|
                    CircuitEnthalpyError::Source(ThermoelectricNetworkError::Seebeck(e)))?.voltage_drop_v;
                spans.push(source);
            }
            let electrical = advance(coefficients.resistors, &inputs, dt,
                &mut || (*cancel.borrow_mut())())?;
            let mut samples = Vec::with_capacity(legs.len());
            let mut heat_w = 0.0;
            let mut electrical_work_j = 0.0;
            for ((leg, &(port, cell_index)), source) in legs.iter().zip(&routing).zip(spans) {
                check()?;
                // Native y is minus the a-to-b voltage-source branch current.
                let current_a = -electrical.native().y[port];
                let transfer = source.transfer(current_a).map_err(|e|
                    CircuitEnthalpyError::Source(ThermoelectricNetworkError::Seebeck(e)))?;
                electrical_work_j = finite(electrical_work_j + dt * inputs[port] * electrical.native().y[port])?;
                let bulk = material.cells[cell_index].cell;
                for (cell, power) in [(leg.terminals[0], transfer.terminal_heat_w[0]),
                    (bulk, transfer.bulk_heat_w), (leg.terminals[1], transfer.terminal_heat_w[1])] {
                    let total = coefficients.sources_w.entry(cell).or_insert(0.0);
                    *total = finite(*total + power)?;
                    heat_w = finite(heat_w + power)?;
                }
                samples.push(ThermoelectricLegSample { source_branch: leg.source_branch,
                    resistor_branch: leg.resistor_branch, terminals: leg.terminals, bulk, material: source, transfer });
            }
            let heat_j = finite(dt * heat_w)?;
            if finite(heat_j + electrical_work_j)?.abs() > controls.total_energy_tolerance_j {
                return Err(invalid("Seebeck port work and reciprocal junction/bulk heat do not close"));
            }
            Ok(DrivenCoefficients { circuit: electrical, sources_w: coefficients.sources_w,
                conductances_w_k: coefficients.conductances_w_k, circuit_trials: 1,
                source: ThermoelectricSample { network: coefficients.source, legs: samples,
                    reversible_electrical_work_j: electrical_work_j, reversible_heat_j: heat_j } })
        }, with_entropy, || (*cancel.borrow_mut())())
}
