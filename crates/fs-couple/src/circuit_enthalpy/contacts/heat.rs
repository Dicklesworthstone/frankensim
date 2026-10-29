//! Material radiation and thermal contacts on the existing electrical/enthalpy solve.
//! All mechanisms see the same trial cells. Bulk conduction, contact resistance,
//! radiosity, electrical loss and time integration keep their original owners.

use super::material::{ElectricalContactNetwork, ElectricalContactNetworkError,
    PreparedElectricalContactNetwork};
use super::{CircuitEnthalpyCandidate, CircuitEnthalpyError, CircuitEnthalpyInput,
    CircuitEnthalpySolve, CircuitThermalCoefficients, ContactThermalCoefficients,
    ContactThermalSample, DepositionError, DepositionSimulation, MaterialCircuit};
use crate::circuit_enthalpy::material::{MaterialNetworkError, MaterialNetworkSample};
use fs_conduction::ConductionError;
use fs_conduction::transient::deposition::contact::{
    ContactLinkSample, MaterialContactLink, apply_material_contacts};
use fs_conduction::transient::deposition::radiation::{
    MaterialRadiation, RadiatingCell, RadiationCellError, RadiationCellSources};
use fs_exec::Cx;
use fs_matdb::ClaimSelection;
use fs_material::state_point::MaterialPropertySelection;

mod extension;
mod observer;
mod layers;
pub use observer::MaterialHeatObserver;
pub mod solid;

/// A single declarative assembly, with optional physical heat-transfer paths.
/// Empty electrical contacts are valid; bulk resistor bindings remain required.
/// None for radiation means no authored enclosure, not missing emissivity data.
/// Thermal and electrical contact resistance are independent source properties.
pub struct MaterialHeatNetwork<'a> {
    /// Existing bulk/electrical material bindings and explicit contact heat splits.
    pub electrical: &'a ElectricalContactNetwork<'a>,
    /// Existing thermal links whose bulk conductance needs series contact resistance.
    pub thermal_contacts: &'a [MaterialContactLink<'a>],
    /// Closed gray-diffuse enclosure on the SAME cell cards and physical conditions.
    pub radiation: Option<&'a MaterialRadiation<'a>>,
}

/// Original source failures, without replacing them by zero flux or a fallback.
#[derive(Debug)]
pub enum HeatNetworkError {
    /// Incompatible material bindings or unrepresentable combined work.
    Invalid(&'static str),
    /// Existing electrical/bulk material query refusal.
    Electrical(ElectricalContactNetworkError),
    /// Existing thermal-contact or series-resistance refusal.
    Contact(ConductionError),
    /// Existing radiation/material/geometry refusal.
    Radiation(RadiationCellError),
    /// Original source-integrated layer conduction refusal.
    Layer(fs_material::thermal::conduction::source::LayerSourceError),
    /// Original shared-timeline checkpoint or thermal-state failure.
    Thermal(DepositionError),
    /// Cancellation before returning a source or any physical endpoint.
    Cancelled,
}
impl From<DepositionError> for HeatNetworkError {
    fn from(e: DepositionError) -> Self {
        match e { DepositionError::Cancelled => Self::Cancelled, other => Self::Thermal(other) }
    }
}
impl From<ElectricalContactNetworkError> for HeatNetworkError {
    fn from(e: ElectricalContactNetworkError) -> Self {
        match e {
            ElectricalContactNetworkError::Thermal(DepositionError::Cancelled)
            | ElectricalContactNetworkError::Bulk(MaterialNetworkError::Thermal(DepositionError::Cancelled)) => Self::Cancelled,
            other => Self::Electrical(other),
        }
    }
}
impl From<ConductionError> for HeatNetworkError {
    fn from(e: ConductionError) -> Self {
        match e { ConductionError::Cancelled { .. } => Self::Cancelled, other => Self::Contact(other) }
    }
}
impl From<RadiationCellError> for HeatNetworkError {
    fn from(e: RadiationCellError) -> Self {
        match e { RadiationCellError::Cancelled => Self::Cancelled, other => Self::Radiation(other) }
    }
}
impl From<fs_material::thermal::conduction::source::LayerSourceError> for HeatNetworkError {
    fn from(e: fs_material::thermal::conduction::source::LayerSourceError) -> Self {
        match e {
            fs_material::thermal::conduction::source::LayerSourceError::Thermal(
                DepositionError::Cancelled,
            ) => Self::Cancelled,
            other => Self::Layer(other),
        }
    }
}
impl core::fmt::Display for HeatNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material heat network: {self:?}")
    }
}
impl std::error::Error for HeatNetworkError {}
fn check(cx: &Cx<'_>) -> Result<(), HeatNetworkError> {
    cx.checkpoint().map_err(|_| HeatNetworkError::Cancelled)
}

/// Original results from EVERY mechanism at the same proposed thermal endpoint.
#[derive(Debug)]
pub struct HeatNetworkSample {
    /// Complete bulk material snapshots, reused by electrical and gas consumers.
    pub material: MaterialNetworkSample,
    /// Original bulk-plus-contact resistance and both endpoint phase states.
    pub thermal_contacts: Vec<ContactLinkSample>,
    /// Source-integrated nonlinear layer fluxes replacing named bulk links.
    /// These transfers already enter the native conductance map, not sources.
    pub integrated_conduction: Vec<fs_material::thermal::conduction::network::LayerLinkSample>,
    /// Original radiosity, source receipts and per-cell signed radiation powers.
    pub radiation: Option<RadiationCellSources>,
}
impl HeatNetworkSample {
    /// Roundoff-net watts into this closed enclosure's cells. Physically this is
    /// internal exchange, NOT external supply; preserve rather than redistribute it.
    pub fn radiation_net_power_w(&self) -> f64 {
        self.radiation.as_ref().map_or(0.0, RadiationCellSources::net_heat_w)
    }
}

/// Selected source laws, not cached physical states or coefficient values.
/// The immutable assembly supplies topology; original owners re-query all values.
/// Keep this across steps to retain sources, or prepare again to select explicitly.
pub struct PreparedHeatNetwork<'n, 'a> {
    network: &'n MaterialHeatNetwork<'a>,
    electrical: PreparedElectricalContactNetwork<'n, 'a>,
    thermal_contacts: Vec<MaterialContactLink<'a>>,
    radiating_cells: Vec<RadiatingCell<'a>>,
    integrated_conduction: Option<layers::IntegratedConduction<'n, 'a>>,
}

impl<'a> MaterialHeatNetwork<'a> {
    /// Select every bulk/interface/emissivity claim once at the accepted state.
    /// Selection and physical validity are delegated to the original evaluators.
    /// Later queries remain pinned, including through overlapping source domains.
    /// Failed preparation returns no partially usable binding.
    pub fn prepare<'n>(&'n self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>)
        -> Result<PreparedHeatNetwork<'n, 'a>, HeatNetworkError>
    {
        check(cx)?;
        let electrical = self.electrical.prepare(thermal, || cx.checkpoint().is_err())?;
        self.prepare_bound(cx, thermal, electrical, self.thermal_contacts,
            self.radiation.map_or(&[], |radiation| radiation.cells))
    }

    // Fresh preparation and topology continuation share the original evaluators.
    // Supplied bindings may already contain retained pins; never select over them.
    fn prepare_bound<'n>(&'n self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>,
        electrical: PreparedElectricalContactNetwork<'n, 'a>,
        contact_bindings: &[MaterialContactLink<'a>], surface_bindings: &[RadiatingCell<'a>])
        -> Result<PreparedHeatNetwork<'n, 'a>, HeatNetworkError>
    {
        check(cx)?;
        let mut initial = electrical.resolve(thermal, || cx.checkpoint().is_err())?;
        let contacts = apply_material_contacts(thermal, contact_bindings,
            &mut initial.bulk.conductances_w_k, || cx.checkpoint().is_err())?;
        let mut thermal_contacts = Vec::with_capacity(contacts.len());
        for sample in &contacts {
            check(cx)?;
            let mut binding = *contact_bindings.iter().find(|c| c.link == sample.link)
                .ok_or(HeatNetworkError::Invalid("thermal contact binding disappeared"))?;
            binding.selection = ClaimSelection::Pinned(sample.interface.receipt().selected);
            thermal_contacts.push(binding);
        }
        let mut radiating_cells = Vec::new();
        if let Some(original) = self.radiation {
            let radiation = MaterialRadiation { enclosure: original.enclosure, cells: surface_bindings,
                power_tolerance_w: original.power_tolerance_w,
                radiosity_tolerance_w_m2: original.radiosity_tolerance_w_m2 };
            let evaluated = radiation.resolve(cx, thermal)?;
            for binding in surface_bindings {
                check(cx)?;
                let cell = self.electrical.bulk.cells.iter().find(|c| c.cell == binding.cell)
                    .ok_or(HeatNetworkError::Invalid("radiation cell lacks its bulk material binding"))?;
                let sample = &evaluated.surfaces()[binding.surface];
                // Different source properties may have different receipts, but
                // they cannot silently describe different states of the same body.
                let bulk_point = cell.point.at(sample.phase)
                    .map_err(ElectricalContactNetworkError::State)?;
                let radiation_point = binding.material.at(sample.phase)
                    .map_err(ElectricalContactNetworkError::State)?;
                if cell.point.material.content_hash() != binding.material.material.content_hash()
                    || cell.point.temperature_axis != binding.material.temperature_axis
                    || cell.point.liquid_fraction_axis != binding.material.liquid_fraction_axis
                    || bulk_point != radiation_point
                {
                    return Err(HeatNetworkError::Invalid("radiation and bulk must share the same complete material conditions"));
                }
                let pins = sample.material.properties().iter()
                    .map(|p| (p.requirement().name().to_owned(), p.answer().receipt.selected)).collect();
                radiating_cells.push(RadiatingCell { surface: binding.surface, cell: binding.cell,
                    material_id: binding.material_id, curve: binding.curve, material: binding.material,
                    selection: MaterialPropertySelection::PinnedByProperty(pins) });
            }
        }
        check(cx)?;
        Ok(PreparedHeatNetwork { network: self, electrical, thermal_contacts, radiating_cells,
            integrated_conduction: None })
    }
}

impl PreparedHeatNetwork<'_, '_> {
    // Native field/circuit consumers may resolve ordinary graph resistors but
    // have no series-contact owner. Never silently drop electrical interfaces.
    pub(crate) fn field_material_network(
        &self,
    ) -> Result<&crate::circuit_enthalpy::material::MaterialNetwork<'_>, HeatNetworkError> {
        if !self.network.electrical.contacts.is_empty() {
            return Err(HeatNetworkError::Invalid(
                "field circuits require their own electrical-contact work and heat owner"));
        }
        Ok(self.network.electrical.bulk)
    }

    // Spatial electrical consumers retain this exact original material assembly.
    // A graph resistor/contact must not disappear when no graph is being stepped.
    pub(crate) fn thermal_material_network(
        &self,
    ) -> Result<&crate::circuit_enthalpy::material::MaterialNetwork<'_>, HeatNetworkError> {
        let electrical = self.network.electrical;
        if !electrical.bulk.resistors.is_empty() || !electrical.contacts.is_empty()
            || !electrical.electrical_bases.is_empty() {
            return Err(HeatNetworkError::Invalid(
                "spatial conduction cannot discard graph resistors, electrical contacts or circuit source bindings"));
        }
        Ok(electrical.bulk)
    }

    /// Original bulk heat routes for consumers composing additional physics.
    pub fn resistor_destinations(&self) -> Vec<super::ResistorCell> {
        self.electrical.resistor_destinations()
    }

    /// Original A/B contact heat partitions; no independent caller heat map.
    pub fn heat_destinations(&self, thermal: &DepositionSimulation<'_>)
        -> Result<Vec<super::ContactHeatCells>, HeatNetworkError> {
        self.electrical.heat_destinations(thermal).map_err(HeatNetworkError::from)
    }

    /// Re-evaluate all mechanisms at ONE trial. Contact resistance is applied to
    /// fresh bulk-only conductances once; radiation is added to the original
    /// source map once. Electrical loss is left to the existing circuit owner.
    /// All intermediate maps are owned here; late refusal cannot modify callers.
    pub fn resolve(&self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>)
        -> Result<ContactThermalCoefficients<HeatNetworkSample>, HeatNetworkError>
    {
        check(cx)?;
        let mut coefficients = self.electrical.resolve(thermal, || cx.checkpoint().is_err())?;
        let integrated_conduction = self.resolve_integrated_conduction(cx, thermal,
            &mut coefficients.bulk.conductances_w_k, &mut coefficients.bulk.source)?;
        let thermal_contacts = apply_material_contacts(thermal, &self.thermal_contacts,
            &mut coefficients.bulk.conductances_w_k, || cx.checkpoint().is_err())?;
        let radiation = if let Some(original) = self.network.radiation {
            let source = MaterialRadiation { enclosure: original.enclosure, cells: &self.radiating_cells,
                power_tolerance_w: original.power_tolerance_w,
                radiosity_tolerance_w_m2: original.radiosity_tolerance_w_m2 }.resolve(cx, thermal)?;
            source.add_to(&mut coefficients.bulk.sources_w)?;
            Some(source)
        } else { None };
        check(cx)?;
        let bulk = coefficients.bulk;
        Ok(ContactThermalCoefficients { contacts: coefficients.contacts,
            bulk: CircuitThermalCoefficients { resistors: bulk.resistors, sources_w: bulk.sources_w,
                conductances_w_k: bulk.conductances_w_k,
                source: HeatNetworkSample { material: bulk.source, thermal_contacts,
                    integrated_conduction, radiation } } })
    }

    /// Advance the original coupled circuit/enthalpy solve with all authored
    /// heat-transfer mechanisms. Same unknowns, source budgets and local gates.
    /// Bulk/interface losses and radiative exchange are never external supply.
    /// The returned standard report retains roundoff-net radiation in its thermal
    /// source account; an ADDITIONAL combined gate excludes that internal power.
    /// True external heat is non_electrical_heat_j - dt*radiation_net_power_w().
    /// Geometry/view factors, caloric laws and C/L remain fixed over the step;
    /// the endpoint/midpoint composition remains generally first order.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_candidate<'m>(&self, cx: &Cx<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
        circuit: &MaterialCircuit, ports: &[f64], controls: CircuitEnthalpySolve, guess_j: &[f64])
        -> Result<CircuitEnthalpyCandidate<'m, ContactThermalSample<HeatNetworkSample>>,
            CircuitEnthalpyError<HeatNetworkError>>
    {
        check(cx).map_err(as_circuit_error)?;
        let routes = self.network.electrical.heat_destinations(thermal)
            .map_err(HeatNetworkError::from).map_err(as_circuit_error)?;
        let destinations = self.network.electrical.bulk.resistors.iter()
            .map(|r| r.destination).collect::<Vec<_>>();
        let proposal = super::step_candidate(thermal, target_s,
            CircuitEnthalpyInput { circuit, ports, destinations: &destinations }, &routes,
            controls, guess_j, |state| self.resolve(cx, state), || cx.checkpoint().is_err());
        let proposal = match proposal {
            Err(CircuitEnthalpyError::Source(HeatNetworkError::Cancelled)) => return Err(CircuitEnthalpyError::Cancelled),
            other => other?,
        };
        let report = proposal.report();
        closed_balance(target_s - thermal.time_s(), report.total_energy_residual_j,
            &report.source.material, controls.total_energy_tolerance_j).map_err(as_circuit_error)?;
        check(cx).map_err(as_circuit_error)?;
        Ok(proposal)
    }
}

fn as_circuit_error(e: HeatNetworkError) -> CircuitEnthalpyError<HeatNetworkError> {
    match e { HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled, other => CircuitEnthalpyError::Source(other) }
}
pub(crate) fn closed_balance(dt: f64, residual: f64, sample: &HeatNetworkSample, tolerance: f64)
    -> Result<(), HeatNetworkError>
{
    let corrected = residual + dt * sample.radiation_net_power_w();
    if !corrected.is_finite() || corrected.abs() > tolerance {
        return Err(HeatNetworkError::Invalid("closed-enclosure power cannot masquerade as external energy"));
    }
    Ok(())
}
