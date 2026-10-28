//! Original frame viscosity heats its own finite source-bound material layers.
//!
//! The native coupled propagator integrates each physical port's work. Common
//! axial strain partitions it by eta*A, and common curvature by eta*I, using
//! the original section integrals. Every whole-layer receiver has the SAME
//! material card, caloric law, reference volume and mass as its frame source.
//! The existing enthalpy/Newton owner advances conduction and reservoirs, and
//! both histories publish only after independent work, heat and storage gates.
//!
//! Geometry, density and E remain fixed; both original viscosity laws may vary
//! with temperature. Each interval uses its accepted starting resistances. Its
//! actual work heats the finite cells, then their accepted endpoint sources
//! rebuild the SAME mechanical propagator for the next interval. This is a
//! first-order thermal/coefficient split; refine the shared clock for accuracy.
//! A dashpot stores no energy, so changing eta adds no parameter-work impulse.
//! Elastic thermal feedback and expansion use [`thermoelastic`]. Melting needs
//! a different material/geometry owner.

pub mod thermoelastic;

use super::{FrameVibrationBasis, FrameVibrationComponent, FrameVibrationError};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, DrivenCoefficients, DriverAccount, ThermalDriver,
    cell_driver_candidate_entropy_option,
    entropy::CircuitEntropyAudit,
    material::{
        MaterialNetwork, MaterialNetworkError, MaterialNetworkSample, PreparedMaterialNetwork,
    },
};
use crate::material_circuit::MaterialCircuitSolve;
use crate::render::schedule::force::coupled::{
    CoupledModalFrame, CoupledModalSystem, ModalCouplingConfig, ModalCouplingError,
    PreparedModalBudget, damping::ConnectionDamping,
};
use fs_exec::CancelGate;
use fs_matdb::EvaluationDecision;
use fs_material::section::ConcentricElasticSection;
use fs_material::section_temperature::ConcentricSectionTemperatureLaw;
use fs_material::state_point::{DENSITY_PROPERTY, MaterialStatePointError, YOUNG_MODULUS_PROPERTY};
use fs_solid::frame2d::FrameViscousSection;
use fs_solver::{NewtonKrylovConfig, NewtonReport};
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionSimulation, EnergyAudit,
};
use std::collections::{BTreeMap, BTreeSet};

/// Spatial correspondence only; heat fractions and thermal mass are not inputs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FrameLayerCell {
    pub element: usize,
    pub layer: usize,
    pub cell: CellId,
}

/// Original thermal solve controls, with no unrelated circuit parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameThermalSolve {
    pub maximum_cells: usize,
    pub maximum_material_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    pub cell_energy_tolerance_j: f64,
    pub thermal_energy_tolerance_j: f64,
    pub total_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}
impl FrameThermalSolve {
    fn driver(self) -> CircuitEnthalpySolve {
        CircuitEnthalpySolve {
            maximum_cells: self.maximum_cells,
            maximum_circuit_trials: self.maximum_material_trials,
            maximum_newton_iterations: self.maximum_newton_iterations,
            difference_step_j: self.difference_step_j,
            cell_energy_tolerance_j: self.cell_energy_tolerance_j,
            thermal_energy_tolerance_j: self.thermal_energy_tolerance_j,
            total_energy_tolerance_j: self.total_energy_tolerance_j,
            newton: self.newton,
            // The shared cell-driver residual does not call a circuit solver.
            circuit: MaterialCircuitSolve {
                energy_tolerance_j: 0.0,
                power_tolerance_w: 0.0,
            },
        }
    }
}

#[derive(Debug)]
pub enum FrameThermalSourceError {
    Frame(FrameVibrationError),
    Material(MaterialStatePointError),
    Network(MaterialNetworkError),
    Port(fs_ports::CoupleError),
}
impl core::fmt::Display for FrameThermalSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "frame material heat: {self:?}")
    }
}
impl std::error::Error for FrameThermalSourceError {}
pub type FrameThermalError = CircuitEnthalpyError<FrameThermalSourceError>;
type Result<T> = core::result::Result<T, FrameThermalError>;

impl From<ModalCouplingError> for FrameThermalError {
    fn from(error: ModalCouplingError) -> Self {
        match error {
            ModalCouplingError::Cancelled => Self::Cancelled,
            other => Self::Source(FrameThermalSourceError::Frame(other.into())),
        }
    }
}
impl From<fs_ports::CoupleError> for FrameThermalError {
    fn from(error: fs_ports::CoupleError) -> Self {
        Self::Source(FrameThermalSourceError::Port(error))
    }
}
fn invalid(message: &'static str) -> FrameThermalError {
    CircuitEnthalpyError::Invalid(message)
}
fn material(error: MaterialStatePointError) -> FrameThermalError {
    CircuitEnthalpyError::Source(FrameThermalSourceError::Material(error))
}
fn network(error: MaterialNetworkError) -> FrameThermalError {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        other => CircuitEnthalpyError::Source(FrameThermalSourceError::Network(other)),
    }
}
fn check(gate: &CancelGate) -> Result<()> {
    if gate.is_requested() {
        Err(CircuitEnthalpyError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable frame heat arithmetic"))
    }
}
fn same_positive(actual: f64, expected: f64) -> bool {
    actual.is_finite()
        && expected.is_finite()
        && actual > 0.0
        && expected > 0.0
        && (actual - expected).abs() <= 1e-10 * expected
}
fn share(part: f64, whole: f64) -> Result<f64> {
    if whole == 0.0 {
        return Ok(0.0);
    }
    let fraction = finite(part / whole)?;
    if part > 0.0 && fraction == 0.0 {
        return Err(invalid("layer heat fraction underflow"));
    }
    Ok(fraction)
}

struct Receiver {
    binding: FrameLayerCell,
    original: DepositedCell,
    temperature_axis: String,
    axial_share: f64,
    bending_share: f64,
}
struct Region<'a> {
    element: usize,
    length_m: f64,
    law: ConcentricSectionTemperatureLaw<'a>,
    receivers: Vec<Receiver>,
}

/// One original modal trajectory and retained source/finite-cell correspondence.
/// The thermal timeline stays with its existing owner and is passed to `step`.
pub struct FrameThermalSystem<'n, 'a> {
    basis: FrameVibrationBasis,
    mechanical: CoupledModalSystem,
    prepared: PreparedMaterialNetwork<'n, 'a>,
    regions: Vec<Region<'a>>,
    viscosity: Vec<(usize, FrameViscousSection)>,
    origin_s: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameLayerHeat {
    pub destination: FrameLayerCell,
    pub axial_heat_j: f64,
    pub bending_heat_j: f64,
}
impl FrameLayerHeat {
    pub fn total_j(self) -> f64 {
        self.axial_heat_j + self.bending_heat_j
    }
}

/// Original entropy integrals and endpoint quadrature, in J/K. The signed
/// temporal/algebraic defects are numerical, never additional physical heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameHeatEntropy {
    pub thermal_storage_change_j_k: f64,
    pub reservoir_entropy_input_j_k: f64,
    pub conduction_generation_j_k: f64,
    /// Sum of actual layer heat divided by its receiving endpoint Kelvin.
    pub material_loss_generation_j_k: f64,
    pub unclassified_source_at_cell_j_k: f64,
    pub has_unclassified_sources: bool,
    pub temporal_defect_j_k: f64,
    pub algebraic_residual_j_k: f64,
    pub accounting_residual_j_k: f64,
}
impl From<CircuitEntropyAudit> for FrameHeatEntropy {
    fn from(a: CircuitEntropyAudit) -> Self {
        Self {
            thermal_storage_change_j_k: a.thermal_storage_change_j_k,
            reservoir_entropy_input_j_k: a.reservoir_entropy_input_j_k,
            conduction_generation_j_k: a.conduction_generation_j_k,
            material_loss_generation_j_k: a.electrical_loss_generation_j_k,
            unclassified_source_at_cell_j_k: a.unclassified_source_at_cell_j_k,
            has_unclassified_sources: a.has_unclassified_sources,
            temporal_defect_j_k: a.temporal_defect_j_k,
            algebraic_residual_j_k: a.algebraic_residual_j_k,
            accounting_residual_j_k: a.accounting_residual_j_k,
        }
    }
}

/// One jointly accepted interval. Mechanical material heat is internal and
/// cancels once; additional authored modal loss leaves through an explicit sink.
pub struct FrameThermalReport {
    pub mechanical: CoupledModalFrame,
    pub thermal: EnergyAudit,
    pub layer_heat: Vec<FrameLayerHeat>,
    pub deposited_material_heat_j: f64,
    pub external_thermal_heat_j: f64,
    pub external_mechanical_loss_j: f64,
    pub combined_energy_residual_j: f64,
    pub entropy: FrameHeatEntropy,
    /// Actual endpoint conductivity/phase evaluation of the original network.
    pub source: MaterialNetworkSample,
    /// Every original section re-evaluated at its actual layer temperatures.
    pub sections: Vec<(usize, ConcentricElasticSection)>,
    /// Source viscosity used during this interval, before its heat was received.
    pub interval_viscosity: Vec<(usize, FrameViscousSection)>,
    /// Accepted endpoint viscosity, installed for the NEXT mechanical interval.
    pub next_viscosity: Vec<(usize, FrameViscousSection)>,
    pub nonlinear: NewtonReport,
    pub maximum_cell_residual_j: f64,
    pub material_trials: usize,
}

impl FrameVibrationComponent {
    /// Bind every viscous member layer to its whole physical caloric cell.
    /// Additional passive network cells and sourced conductive links are valid.
    /// No resistor branch is allowed: its electrical owner must close that work.
    /// Restored nonzero modal state survives; the thermal clock sets the origin.
    #[allow(clippy::too_many_arguments)]
    pub fn into_thermal_system<'n, 'a>(
        self,
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        destinations: &[FrameLayerCell],
        config: ModalCouplingConfig,
        prepared: Option<PreparedModalBudget>,
        gate: &CancelGate,
    ) -> Result<FrameThermalSystem<'n, 'a>> {
        let (binding, regions, viscosity) =
            bind_regions(&self.basis, thermal, source, destinations, true, gate)?;
        let (basis, mechanical) =
            self.into_system(config, prepared, gate)
                .map_err(|e| match e {
                    FrameVibrationError::Cancelled => CircuitEnthalpyError::Cancelled,
                    other => CircuitEnthalpyError::Source(FrameThermalSourceError::Frame(other)),
                })?;
        check(gate)?;
        Ok(FrameThermalSystem {
            basis,
            mechanical,
            prepared: binding,
            regions,
            viscosity,
            origin_s: thermal.time_s(),
        })
    }
}

// Share the original source/card/caloric/inventory admission across the two
// physical thermal consumers. The thermoelastic child supplies the original
// Helmholtz law instead of the frozen-E restriction.
fn bind_regions<'n, 'a>(
    basis: &FrameVibrationBasis,
    thermal: &DepositionSimulation<'_>,
    source: &'n MaterialNetwork<'a>,
    destinations: &[FrameLayerCell],
    frozen_elasticity: bool,
    gate: &CancelGate,
) -> Result<(
    PreparedMaterialNetwork<'n, 'a>,
    Vec<Region<'a>>,
    Vec<(usize, FrameViscousSection)>,
)> {
    check(gate)?;
    if !source.resistors.is_empty() {
        return Err(invalid("frame heat needs a thermal-only material network"));
    }
    let binding = source
        .prepare(thermal, || gate.is_requested())
        .map_err(network)?;
    let operator = basis
        .material_damping()
        .ok_or_else(|| invalid("frame heat requires its original material viscosity"))?;
    let count: usize = operator
        .sections()
        .iter()
        .map(|(_, s)| s.elastic().layers().len())
        .sum();
    if destinations.len() != count {
        return Err(invalid(
            "exactly one destination for every original member layer required",
        ));
    }
    let routes: BTreeMap<_, _> = destinations
        .iter()
        .map(|d| ((d.element, d.layer), d.cell))
        .collect();
    if routes.len() != count {
        return Err(invalid(
            "exactly one destination for every original member layer required",
        ));
    }
    let mut used = BTreeSet::new();
    let mut regions = Vec::with_capacity(operator.sections().len());
    for (index, (element, section)) in operator.sections().iter().enumerate() {
        check(gate)?;
        // The native operator's documented symmetric-bending factor has
        // rotational weights L: use that SAME reference geometry, not an
        // independently authored heat length or a modal wavelength.
        let length = operator.ports()[3 * index + 1].displacement_weights[2];
        let mut cards = Vec::new();
        let mut axes = Vec::new();
        let mut receivers = Vec::new();
        for (layer, (solid, geometry)) in section
            .elastic()
            .layers()
            .iter()
            .zip(section.elastic().layer_integrals())
            .enumerate()
        {
            check(gate)?;
            let cell = *routes
                .get(&(*element, layer))
                .ok_or_else(|| invalid("foreign or missing frame-layer heat destination"))?;
            if !used.insert(cell) {
                return Err(invalid(
                    "distinct physical layers require distinct whole-cell inventories",
                ));
            }
            let caloric = source
                .cells
                .iter()
                .find(|b| b.cell == cell)
                .ok_or_else(|| {
                    invalid("frame layer lacks its original caloric material binding")
                })?;
            let original = thermal
                .cells()
                .get(&cell)
                .ok_or_else(|| invalid("inactive frame heat destination"))?;
            let density = geometry.line_density_kg_m() / geometry.area_m2();
            if !same_positive(original.geometry.volume_m3, geometry.area_m2() * length)
                || !same_positive(original.mass_kg, geometry.line_density_kg_m() * length)
                || caloric.curve.knots().iter().any(|k| {
                    k.liquid_mass_fraction != 0.0 || !same_positive(k.bulk_density_kg_m3, density)
                })
                || solid.material.query() != *caloric.point.context
                || thermal.temperature_k(cell)?.to_bits()
                    != caloric.point.context.axes()[caloric.point.temperature_axis].to_bits()
            {
                return Err(invalid(
                    "frame layer must retain its whole fixed-solid mass, volume and original physical source context",
                ));
            }
            for key in [DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY]
                .into_iter()
                .filter(|key| frozen_elasticity || *key == DENSITY_PROPERTY)
            {
                let coefficient = solid
                    .material
                    .property(key)
                    .ok_or_else(|| invalid("frame heat lost an original mechanical source"))?;
                if !coefficient.source_axes().is_empty()
                    || coefficient.answer().receipt.decision
                        != EvaluationDecision::ConstantWithinValidity
                {
                    return Err(invalid(
                        "frame self-heating requires source-declared temperature-independent density and elasticity",
                    ));
                }
            }
            receivers.push(Receiver {
                binding: FrameLayerCell {
                    element: *element,
                    layer,
                    cell,
                },
                original: original.clone(),
                temperature_axis: caloric.point.temperature_axis.to_owned(),
                axial_share: share(
                    section.axial().layer_viscous_rigidity_n_s()[layer],
                    section.axial().viscous_rigidity_n_s(),
                )?,
                bending_share: share(
                    section.bending().layer_viscous_stiffness_n_m2_s()[layer],
                    section.bending().viscous_stiffness_n_m2_s(),
                )?,
            });
            cards.push(caloric.point.material);
            axes.push(caloric.point.temperature_axis);
        }
        let law =
            ConcentricSectionTemperatureLaw::try_new(section.elastic().clone(), &cards, &axes)
                .map_err(material)?;
        regions.push(Region {
            element: *element,
            length_m: length,
            law,
            receivers,
        });
    }
    let viscosity = operator.sections().to_vec();
    Ok((binding, regions, viscosity))
}

impl FrameThermalSystem<'_, '_> {
    /// Fixed mass/elastic reference and physical maps. Its source damping
    /// records describe construction; `viscous_sections()` is the live state.
    pub const fn basis(&self) -> &FrameVibrationBasis {
        &self.basis
    }
    pub const fn mechanical(&self) -> &CoupledModalSystem {
        &self.mechanical
    }
    /// Current source-bound resistances at each accepted layer temperature.
    pub fn viscous_sections(&self) -> &[(usize, FrameViscousSection)] {
        &self.viscosity
    }

    /// Advance one original held-force modal sample and the same finite thermal
    /// interval. Heat uses independently integrated port work, never a stored-
    /// energy defect. The native implicit thermal residual evaluates sourced
    /// conduction at the endpoint. Viscosity is frozen over this interval, then
    /// refreshed at its accepted temperature without changing stored mechanical
    /// energy or histories. Refine the clock for this first-order feedback split.
    ///
    /// Start/inside-interval thermal events retain the original residual's
    /// refusal rules. On failure, BOTH accepted histories, modal sample/ledger,
    /// and thermal event cursor remain unchanged. Retry uses the same inputs.
    pub fn step(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        external_modal_forces: &[f64],
        controls: FrameThermalSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
    ) -> Result<FrameThermalReport> {
        check(gate)?;
        let sample = self.mechanical.samples_rendered();
        let next = sample
            .checked_add(1)
            .ok_or_else(|| invalid("frame heat clock overflow"))?;
        let dt = self.mechanical.sample_period_s();
        let start = finite(self.origin_s + sample as f64 * dt)?;
        let target = finite(self.origin_s + next as f64 * dt)?;
        if thermal.time_s().to_bits() != start.to_bits() || target <= start {
            return Err(invalid(
                "mechanical and thermal frame histories must share one advancing clock",
            ));
        }
        let before = self.mechanical.total_energy_j()?;
        let regions = &self.regions;
        let viscosity = &self.viscosity;
        for (region, (_, section)) in regions.iter().zip(viscosity) {
            for (receiver, layer) in region.receivers.iter().zip(section.elastic().layers()) {
                let temperature = thermal.temperature_k(receiver.binding.cell)?;
                let source = layer
                    .material
                    .query_point()
                    .iter()
                    .find(|(axis, _)| axis == &receiver.temperature_axis)
                    .map(|(_, temperature)| *temperature)
                    .ok_or_else(|| {
                        invalid("frame viscosity lost its original temperature coordinate")
                    })?;
                if temperature.to_bits() != source.to_bits() {
                    return Err(invalid(
                        "frame viscosity and finite cells must retain their jointly accepted thermal state",
                    ));
                }
            }
        }
        let prepared = &self.prepared;
        let (next_thermal, report, damping, next_shares, next_viscosity) = self.mechanical.step_with_admission(
            external_modal_forces, gate, |mechanical, frame, transfers| {
                check(gate)?;
                if transfers.len() != 3 * regions.len() { return Err(invalid("frame heat port layout changed")); }
                let mut layer_heat = Vec::new();
                let mut powers = BTreeMap::new();
                for (index, region) in regions.iter().enumerate() {
                    check(gate)?;
                    let mut loss = [0.0; 3];
                    for (k, value) in loss.iter_mut().enumerate() {
                        *value = finite(-transfers[3 * index + k].energy_into().value())?;
                        if *value < 0.0 { return Err(invalid("material dashpot port supplied negative heat")); }
                    }
                    for receiver in &region.receivers {
                        let partition = |share: f64, joules: f64| -> Result<f64> {
                            let heat = finite(share * joules)?;
                            if share > 0.0 && joules > 0.0 && heat == 0.0 {
                                return Err(invalid("physical layer heat underflow"));
                            }
                            Ok(heat)
                        };
                        let heat = FrameLayerHeat {
                            destination: receiver.binding,
                            axial_heat_j: partition(receiver.axial_share, loss[0])?,
                            bending_heat_j: partition(receiver.bending_share, finite(loss[1] + loss[2])?)?,
                        };
                        let power = finite(heat.total_j() / (target - start))?;
                        if heat.total_j() != 0.0 && power == 0.0 { return Err(invalid("frame heat power underflow")); }
                        powers.insert(receiver.binding.cell, power);
                        layer_heat.push(heat);
                    }
                }
                let outside = finite(frame.component_dissipation_j + frame.material_dissipation_j + frame.radiation_dissipation_j)?;
                let (candidate, (), audit, entropy) = cell_driver_candidate_entropy_option(
                    thermal, target, controls.driver(), energy_guess_j,
                    |state, _remaining| {
                        check(gate)?;
                        let source = prepared.resolve(state, || gate.is_requested()).map_err(network)?;
                        let mut sections = Vec::with_capacity(regions.len());
                        for (region, (_, previous)) in regions.iter().zip(viscosity) {
                            check(gate)?;
                            let mut temperatures = Vec::with_capacity(region.receivers.len());
                            for receiver in &region.receivers {
                                let original = &receiver.original;
                                let cell = state.cells().get(&receiver.binding.cell)
                                    .ok_or_else(|| invalid("frame heat receiver disappeared"))?;
                                if cell.id != original.id || cell.material_id != original.material_id
                                    || cell.segment_id != original.segment_id || cell.body_id != original.body_id
                                    || cell.mass_kg.to_bits() != original.mass_kg.to_bits() || cell.geometry != original.geometry {
                                    return Err(invalid("frame heat receiver changed its original material inventory or geometry"));
                                }
                                temperatures.push(state.temperature_k(cell.id)?);
                            }
                            let section = region.law.at_temperatures_k(&temperatures).map_err(material)?;
                            let section = FrameViscousSection::try_new(section,
                                previous.axial().viscosity_property(), previous.bending().viscosity_property())
                                .map_err(|error| CircuitEnthalpyError::Source(FrameThermalSourceError::Frame(error.into())))?;
                            let band = section.omega_band_rad_s();
                            if band.0 != 0.0 || band.1 < core::f64::consts::PI / dt {
                                return Err(invalid("refreshed frame viscosity must support DC through the actual Nyquist band"));
                            }
                            sections.push((region.element, section));
                        }
                        Ok(DrivenCoefficients {
                            circuit: MechanicalCandidate { frame: frame.clone(), stored_before_j: before, outside_loss_j: outside },
                            sources_w: source.sources_w, conductances_w_k: source.conductances_w_k,
                            source: (source.source, sections), circuit_trials: 1,
                        })
                    },
                    |_driver, _dt, checkpoint| { checkpoint()?; Ok(powers.clone()) },
                    true, || gate.is_requested(),
                )?;
                let entropy = entropy.ok_or_else(|| invalid("frame heat requires its original entropy account"))?;
                let mut updates = Vec::with_capacity(3 * regions.len());
                let mut next_shares = Vec::new();
                for (index, (region, (_, section))) in regions.iter().zip(&audit.source.1).enumerate() {
                    check(gate)?;
                    let axial = finite(section.axial().viscous_rigidity_n_s() / region.length_m)?;
                    let bending = finite(section.bending().viscous_stiffness_n_m2_s()
                        / region.length_m / region.length_m / region.length_m)?;
                    if (section.axial().viscous_rigidity_n_s() > 0.0 && axial == 0.0)
                        || (section.bending().viscous_stiffness_n_m2_s() > 0.0 && bending == 0.0) {
                        return Err(invalid("refreshed frame viscosity is unrepresentable on the original geometry"));
                    }
                    for (k, damping_n_s_m) in [axial, finite(3.0 * bending)?, bending].into_iter().enumerate() {
                        updates.push(ConnectionDamping { connection: 3 * index + k, damping_n_s_m });
                    }
                    for layer in 0..region.receivers.len() {
                        next_shares.push((
                            share(section.axial().layer_viscous_rigidity_n_s()[layer], section.axial().viscous_rigidity_n_s())?,
                            share(section.bending().layer_viscous_stiffness_n_m2_s()[layer], section.bending().viscous_stiffness_n_m2_s())?,
                        ));
                    }
                }
                let damping = mechanical.prepare_staged_connection_damping(&updates, gate)?;
                let sections = audit.source.1.iter().map(|(element, section)| (*element, section.elastic().clone())).collect();
                let report = FrameThermalReport {
                    mechanical: audit.electrical, thermal: audit.thermal, layer_heat,
                    deposited_material_heat_j: audit.deposited_joule_j,
                    external_thermal_heat_j: audit.non_electrical_heat_j,
                    external_mechanical_loss_j: outside,
                    combined_energy_residual_j: audit.total_energy_residual_j,
                    entropy: entropy.into(), source: audit.source.0, sections,
                    interval_viscosity: viscosity.to_vec(), next_viscosity: audit.source.1,
                    nonlinear: audit.nonlinear, maximum_cell_residual_j: audit.maximum_cell_residual_j,
                    material_trials: audit.circuit_trials,
                };
                let next_viscosity = report.next_viscosity.clone();
                Ok((candidate, report, damping, next_shares, next_viscosity))
            },
        )?;
        // The native stage/admit/publish boundary made the final cancellation
        // check. These infallible moves publish the matching finite-cell state.
        self.mechanical.install_connection_damping(damping);
        for (receiver, (axial, bending)) in self
            .regions
            .iter_mut()
            .flat_map(|region| &mut region.receivers)
            .zip(next_shares)
        {
            receiver.axial_share = axial;
            receiver.bending_share = bending;
        }
        self.viscosity = next_viscosity;
        *thermal = next_thermal;
        Ok(report)
    }
}

struct MechanicalCandidate {
    frame: CoupledModalFrame,
    stored_before_j: f64,
    outside_loss_j: f64,
}
impl ThermalDriver for MechanicalCandidate {
    type Model = ();
    type Step = CoupledModalFrame;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.frame.connection_dissipation_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.frame.total_energy_j() - self.stored_before_j,
            supplied_j: self.frame.external_work_j - self.outside_loss_j,
        }
    }
    fn into_parts(self) -> ((), CoupledModalFrame) {
        ((), self.frame)
    }
}
