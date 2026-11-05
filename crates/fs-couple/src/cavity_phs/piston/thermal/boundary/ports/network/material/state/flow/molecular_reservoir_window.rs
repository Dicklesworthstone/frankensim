//! Accepted open streams, fixed physical budgets and complete-window histories.
//!
//! Gas transport retains its original source-reference cancellation. Feedback
//! retains the existing reaction/caloric owner. Actual endpoint species,
//! circuit, finite-cell and mechanical storage independently close the window;
//! no sum of changing local enthalpy charts is treated as held-bath energy.

use super::*;
use crate::cavity_phs::piston::thermal::boundary::ports::network::material::state::flow::diffusion::coupled::{
    CoupledDiffusionError, CoupledReactionExtent, PhysicalSubstep,
};
use crate::port_adoption::{PortTransfer, StepAccount};
use fs_gas::material::transport::parcel::adaptive::{
    ParcelMolecularIntegratedHeat, ParcelMolecularIntegratedTransfer,
};
use fs_gas::material::transport::parcel::reservoir::adaptive::{
    ReservoirMolecularIntegratedHeat, ReservoirMolecularIntegratedSpecies,
};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::MAX_SNAPSHOT_TRANSFERS;
use fs_thermochem::{ReactionId, SpeciesId};

/// Whole-window gas and wall accounts, retaining each original caloric leg.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularReservoirWindowChamberBalance {
    pub chamber: u64,
    pub mass_change_kg: f64,
    /// Sum of native intended source-reference species increments plus actual
    /// sensible storage. Inventory rounding retains its separate mass budget.
    /// This is not a held bath's integrated boundary enthalpy.
    pub transport_energy_change_j: f64,
    /// Sum of original post-transport reaction/heat/motion gas-energy increments.
    pub feedback_energy_change_j: f64,
    /// Independently measured original-to-final wall kinetic-energy change.
    pub kinetic_energy_change_j: f64,
    /// Sum of accepted gas entropy-storage observations; not entropy production.
    pub entropy_change_j_k: f64,
    pub gas_boundary_heat_j: f64,
    /// Sum of actual accepted Q/T_midpoint, preserving changes of heat direction.
    pub gas_boundary_entropy_j_k: f64,
    pub external_force_work_j: f64,
    pub exterior_pressure_work_j: f64,
    pub support_work_j: f64,
    pub support_impulse_n_s: f64,
    pub momentum_residual_n_s: f64,
    pub transport_energy_residual_j: f64,
    /// Original caloric increments against actual endpoint mechanical work and
    /// the independently retained surface heat and prescribed whole-window heat.
    pub feedback_energy_residual_j: f64,
    pub energy_residual_j: f64,
    /// Absolute accepted defects plus the independent feedback account drift.
    pub energy_defect_bound_j: f64,
    /// Any accepted gross species stream, including a cancelling throughflow.
    pub transported_species: bool,
    /// Any original extra heat path, even if its net whole-window energy is zero.
    pub transported_heat: bool,
}

/// Complete-window physical budgets, independent of subdivision and rejection.
///
/// The energy account composes reference-cancelled native transport with the
/// original feedback calorics, actual finite-body/circuit storage and endpoint
/// mechanics. It never subtracts gross absolute gas U and bath enthalpy. Native
/// stream receipts preserve directions; a zero net exchange does not erase them.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularReservoirWindowBalance {
    pub thermal_energy_change_j: f64,
    pub electrical_energy_change_j: f64,
    pub electrical_supply_work_j: f64,
    pub deposited_joule_j: f64,
    pub external_network_heat_j: f64,
    pub prescribed_gas_heat_j: f64,
    pub non_electrical_heat_j: f64,
    /// Sum of native per-step maximum cell defects, bounding every fixed cell.
    pub maximum_cell_residual_bound_j: f64,
    pub thermal_energy_defect_bound_j: f64,
    pub circuit_energy_defect_bound_j: f64,
    pub network_energy_defect_bound_j: f64,
    pub network_energy_residual_j: f64,
    pub transport_energy_residual_j: f64,
    /// Independent window feedback closure, not just summed substep residuals.
    pub feedback_energy_residual_j: f64,
    pub accumulated_feedback_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    /// Sum of per-step largest absolute native global or chamber transport
    /// defects, bounded by the original single native transport allowance.
    pub transport_energy_defect_bound_j: f64,
    /// Absolute accepted feedback defects plus independent endpoint-account drift.
    pub feedback_energy_defect_bound_j: f64,
    pub total_energy_defect_bound_j: f64,
    pub maximum_chamber_energy_residual_j: f64,
    pub net_reservoir_mass_kg: f64,
    pub total_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    /// Sum of the original bath-inclusive TRANSPORT-stage entropy productions.
    /// Chemistry and heat/work feedback are excluded; this is not a whole-
    /// assembly or original-to-final bath entropy-production certificate.
    pub transport_entropy_production_j_k: f64,
    /// Sum of negative native transport productions, using one original allowance.
    pub transport_entropy_deficit_bound_j_k: f64,
    pub species: Vec<MolecularReservoirSpeciesBalance>,
    pub transfers: Vec<ParcelMolecularIntegratedTransfer>,
    pub heat: Vec<ParcelMolecularIntegratedHeat>,
    pub reservoir_species: Vec<ReservoirMolecularIntegratedSpecies>,
    pub reservoir_heat: Vec<ReservoirMolecularIntegratedHeat>,
    pub reaction_extents: Vec<CoupledReactionExtent>,
    pub chambers: Vec<MolecularReservoirWindowChamberBalance>,
}

type Error = CoupledDiffusionError;
fn check(cx: &Cx<'_>) -> Result<(), Error> {
    cx.checkpoint().map_err(|_| Error::Cancelled)
}
fn finite(value: f64) -> Result<f64, Error> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::Balance(
            "open molecular window arithmetic is unrepresentable",
        ))
    }
}
#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), Error> {
        let value = finite(value)?;
        let next = finite(self.value + value)?;
        self.correction = finite(
            self.correction
                + if self.value.abs() >= value.abs() {
                    (self.value - next) + value
                } else {
                    (value - next) + self.value
                },
        )?;
        self.value = next;
        Ok(())
    }
    fn total(&self) -> Result<f64, Error> {
        finite(self.value + self.correction)
    }
}
fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64, Error> {
    let mut result = Sum::default();
    for value in values {
        result.add(value)?;
    }
    result.total()
}
#[derive(Default)]
struct Directional {
    positive: Sum,
    negative: Sum,
    net: Sum,
}
impl Directional {
    fn add(&mut self, value: f64) -> Result<(), Error> {
        self.net.add(value)?;
        if value < 0.0 {
            self.negative.add(-value)
        } else {
            self.positive.add(value)
        }
    }
}
#[derive(Default)]
struct SpeciesLedger {
    internal: Sum,
    reservoir: Sum,
    reaction: Sum,
}
#[derive(Default)]
struct ChamberLedger {
    species: BTreeMap<SpeciesId, SpeciesLedger>,
    mass: Sum,
    transport_energy: Sum,
    feedback_energy: Sum,
    transport_residual: Sum,
    transport_defect: Sum,
    feedback_residual: Sum,
    feedback_defect: Sum,
    kinetic: Sum,
    entropy: Sum,
    heat: Sum,
    heat_entropy: Sum,
    evaluated_surface_heat: Sum,
    force_work: Sum,
    pressure_work: Sum,
    support_work: Sum,
    gas_impulse: Sum,
    support_impulse: Sum,
    moved: bool,
    transported_heat: bool,
}
#[derive(Default)]
struct Accepted {
    chambers: BTreeMap<u64, ChamberLedger>,
    transfers: BTreeMap<(u64, u64, SpeciesId), (Sum, Sum)>,
    heat: BTreeMap<[u64; 2], Directional>,
    reservoir_species: BTreeMap<(u64, u64, SpeciesId), Directional>,
    reservoir_heat: BTreeMap<(u64, u64), Directional>,
    extents: BTreeMap<(u64, ReactionId), Sum>,
    supply: Sum,
    joule: Sum,
    circuit_heat: Sum,
    external: Sum,
    non_electrical: Sum,
    thermal_input: Sum,
    cell_defect: Sum,
    thermal_defect: Sum,
    circuit_defect: Sum,
    network_defect: Sum,
    transport_residual: Sum,
    feedback_residual: Sum,
    transport_defect: Sum,
    feedback_defect: Sum,
    transport_entropy: Sum,
    transport_entropy_deficit: Sum,
    bath_mass: Sum,
}
impl Accepted {
    fn chamber(&mut self, id: u64) -> Result<&mut ChamberLedger, Error> {
        self.chambers
            .get_mut(&id)
            .ok_or(Error::Balance("unknown accepted open-window chamber"))
    }
    fn stream_entries(&self) -> Option<usize> {
        let counts = [
            self.transfers.len(),
            self.heat.len(),
            self.reservoir_species.len(),
            self.reservoir_heat.len(),
            self.extents.len(),
        ];
        counts
            .into_iter()
            .try_fold(0usize, |count, more| count.checked_add(more))
    }
    fn bound(&self) -> Result<(), Error> {
        let entries = self.stream_entries().and_then(|count| {
            self.chambers.values().try_fold(count, |count, chamber| {
                count.checked_add(chamber.species.len())
            })
        });
        if entries.is_none_or(|count| count > MAX_SNAPSHOT_TRANSFERS) {
            Err(Error::Balance(
                "combined accepted open-window ledgers exceed native transfer work bound",
            ))
        } else {
            Ok(())
        }
    }

    fn record(
        &mut self,
        initial: &MaterialExchangeState<'_>,
        step: &PhysicalSubstep<MolecularReservoirAdvance>,
        reactions: &[ChamberReaction<'_>],
        cx: &Cx<'_>,
    ) -> Result<(), Error> {
        check(cx)?;
        let advance = &step.advance;
        let flow = &advance.flow;
        let report = &advance.exchange.report;
        let dt = finite(step.end_s - step.start_s)?;
        if dt <= 0.0 || flow.duration_s != dt {
            return Err(Error::Balance(
                "accepted open receipt has a different physical interval",
            ));
        }
        self.supply.add(report.electrical.native.supplied)?;
        self.joule.add(report.deposited_joule_j)?;
        self.circuit_heat.add(report.electrical.heat_j)?;
        self.external
            .add(advance.exchange.step.external_network_heat_j)?;
        self.non_electrical.add(report.non_electrical_heat_j)?;
        self.thermal_input.add(report.thermal.external_heat_j)?;
        self.thermal_input.add(report.thermal.incoming_energy_j)?;
        self.cell_defect.add(report.maximum_cell_residual_j.abs())?;
        self.thermal_defect.add(report.thermal.residual_j.abs())?;
        self.circuit_defect
            .add(report.electrical.energy_residual_j.abs())?;
        self.network_defect
            .add(report.total_energy_residual_j.abs())?;
        self.transport_residual.add(flow.audit.energy_residual_j)?;
        self.feedback_residual
            .add(advance.exchange.step.total_energy_residual_j)?;
        self.transport_entropy
            .add(flow.audit.entropy_production_j_k)?;
        self.transport_entropy_deficit
            .add((-finite(flow.audit.entropy_production_j_k)?).max(0.0))?;
        let mut transport_defect = finite(flow.audit.energy_residual_j)?
            .abs()
            .max(finite(flow.exchange.total_energy_residual_j)?.abs());
        let mut feedback_defect = finite(advance.exchange.step.total_energy_residual_j)?.abs();
        for flux in &flow.exchange.fluxes {
            check(cx)?;
            let transfer = &flux.transfer;
            let gross = self
                .transfers
                .entry((transfer.donor, transfer.receiver, transfer.species.clone()))
                .or_default();
            gross.0.add(transfer.mass_kg)?;
            gross.1.add(flux.carried_enthalpy_j)?;
            for (id, sign) in [(transfer.donor, -1.0), (transfer.receiver, 1.0)] {
                let chamber = self.chamber(id)?;
                chamber
                    .species
                    .entry(transfer.species.clone())
                    .or_default()
                    .internal
                    .add(sign * transfer.mass_kg)?;
                chamber.moved |= transfer.mass_kg > 0.0;
            }
        }
        for flux in &flow.exchange.boundaries {
            check(cx)?;
            let transfer = &flux.transfer;
            let chamber = self.chamber(flux.parcel)?;
            chamber
                .species
                .entry(transfer.species.clone())
                .or_default()
                .reservoir
                .add(transfer.net_mass_into_gas_kg)?;
            chamber.moved |= transfer.incoming_mass_kg > 0.0 || transfer.outgoing_mass_kg > 0.0;
            self.bath_mass.add(transfer.net_mass_into_gas_kg)?;
        }
        // Keep directional held-bath streams separately from the native local
        // sensible correction, whose reference temperature changes each step.
        for link in &flow.links {
            check(cx)?;
            if link.mass_into_gas_kg.len() != flow.species.len() {
                return Err(Error::Balance(
                    "reservoir face lost its canonical species report",
                ));
            }
            for (species, &mass) in flow.species.iter().zip(&link.mass_into_gas_kg) {
                check(cx)?;
                if mass != 0.0 {
                    self.reservoir_species
                        .entry((link.id, link.parcel, species.clone()))
                        .or_default()
                        .add(mass)?;
                }
            }
            self.reservoir_heat
                .entry((link.id, link.parcel))
                .or_default()
                .add(link.reduced_heat_into_gas_j)?;
            if link.reduced_heat_into_gas_j != 0.0 {
                self.chamber(link.parcel)?.transported_heat = true;
            }
        }
        for heat in &flow.exchange.heat {
            check(cx)?;
            self.heat
                .entry(heat.parcels)
                .or_default()
                .add(heat.heat_left_to_right_j)?;
            if heat.heat_left_to_right_j != 0.0 {
                for id in heat.parcels {
                    self.chamber(id)?.transported_heat = true;
                }
            }
        }
        // Combined native heat may cancel between parallel original faces.
        for link in &flow.internal_links {
            check(cx)?;
            if link.heat_correction_left_to_right_j != 0.0 {
                for id in link.parcels {
                    self.chamber(id)?.transported_heat = true;
                }
            }
        }
        for chemical in &report.source.reactions {
            check(cx)?;
            let binding = reactions
                .iter()
                .find(|reaction| reaction.chamber == chemical.chamber)
                .ok_or(Error::Balance("unbound accepted open-window reaction"))?;
            let sources = binding.reaction.catalog().thermodynamics().components();
            for species in &chemical.reaction.species {
                check(cx)?;
                let source = sources
                    .iter()
                    .find(|source| source.species() == &species.species)
                    .ok_or(Error::Balance(
                        "open-window reaction omitted its original species source",
                    ))?;
                self.chamber(chemical.chamber)?
                    .species
                    .entry(species.species.clone())
                    .or_default()
                    .reaction
                    .add(species.stoichiometric_change_mol * source.molar_mass().value())?;
            }
            for (reaction, extent) in chemical.reaction.reaction_extents() {
                self.extents
                    .entry((chemical.chamber, reaction.clone()))
                    .or_default()
                    .add(extent)?;
            }
        }
        let native: BTreeMap<_, _> = flow
            .exchange
            .parcels
            .iter()
            .map(|row| (row.id, row))
            .collect();
        let feedback: BTreeMap<_, _> = report
            .source
            .chambers
            .iter()
            .map(|row| (row.id, row))
            .collect();
        let accounts: BTreeMap<_, _> = advance
            .balance
            .chambers
            .iter()
            .map(|row| (row.chamber, row))
            .collect();
        if !initial.chambers.keys().eq(native.keys())
            || !initial.chambers.keys().eq(feedback.keys())
            || !initial.chambers.keys().eq(accounts.keys())
        {
            return Err(Error::Balance(
                "open-window receipt omitted an original chamber",
            ));
        }
        for (&id, original) in &initial.chambers {
            check(cx)?;
            let account = accounts[&id];
            let boundary = &feedback[&id].boundary;
            let piston = &boundary.piston;
            let chemical = report
                .source
                .reactions
                .iter()
                .find(|reaction| reaction.chamber == id);
            let feedback_residual = finite(
                chemical.map_or(boundary.total_energy_residual_j, |reaction| {
                    reaction.boundary_energy_residual_j
                }),
            )?;
            let transport_residual = finite(native[&id].energy_residual_j)?;
            transport_defect = transport_defect.max(transport_residual.abs());
            feedback_defect = feedback_defect.max(feedback_residual.abs());
            let chamber = self.chamber(id)?;
            chamber.mass.add(account.mass_change_kg)?;
            chamber
                .transport_energy
                .add(account.transport_energy_change_j)?;
            chamber
                .feedback_energy
                .add(account.feedback_energy_change_j)?;
            chamber.transport_residual.add(transport_residual)?;
            chamber.transport_defect.add(transport_residual.abs())?;
            chamber.feedback_residual.add(feedback_residual)?;
            chamber.feedback_defect.add(feedback_residual.abs())?;
            chamber.kinetic.add(piston.kinetic_energy_change_j)?;
            chamber.entropy.add(account.entropy_change_j_k)?;
            chamber.heat.add(piston.gas.supplied_heat_j)?;
            chamber
                .heat_entropy
                .add(piston.gas.supplied_heat_j / piston.midpoint.state().temperature)?;
            chamber
                .evaluated_surface_heat
                .add(boundary.evaluated_heat_j)?;
            chamber.force_work.add(piston.external_force_work_j)?;
            chamber.pressure_work.add(piston.exterior_pressure_work_j)?;
            chamber.support_work.add(piston.support_work_j)?;
            chamber
                .gas_impulse
                .add(dt * original.area_m2() * piston.midpoint.state().pressure)?;
            chamber.support_impulse.add(piston.support_impulse_n_s)?;
        }
        self.transport_defect.add(transport_defect)?;
        self.feedback_defect.add(feedback_defect)?;
        self.bound()?;
        check(cx)
    }
}

fn inventory(gas: &ClosedGasParcel) -> Result<BTreeMap<SpeciesId, f64>, Error> {
    let amount = finite(
        gas.mass_kg()
            / gas
                .state()
                .thermodynamics()
                .properties()
                .molar_mass()
                .value(),
    )?;
    gas.model()
        .thermodynamics()
        .components()
        .iter()
        .zip(gas.model().thermodynamics().mole_composition().fractions())
        .map(|(source, &fraction)| {
            let mass = finite((amount * fraction) * source.molar_mass().value())?;
            if mass < 0.0 || (fraction > 0.0 && mass == 0.0) {
                return Err(Error::Balance(
                    "actual open-window species inventory is unrepresentable",
                ));
            }
            Ok((source.species().clone(), mass))
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
pub(super) fn summarize(
    initial: &MaterialExchangeState<'_>,
    final_state: &MaterialExchangeState<'_>,
    cx: &Cx<'_>,
    ports: &[f64],
    drives: &[ChamberDrive],
    reactions: &[ChamberReaction<'_>],
    molecular_controls: ParcelMolecularSolve,
    controls: NetworkSolve,
    steps: &[PhysicalSubstep<MolecularReservoirAdvance>],
) -> Result<MolecularReservoirWindowBalance, Error> {
    check(cx)?;
    let duration = finite(final_state.time_s() - initial.time_s())?;
    if steps.is_empty()
        || duration <= 0.0
        || !initial.chambers.keys().eq(final_state.chambers.keys())
        || !initial
            .thermal
            .cells()
            .keys()
            .eq(final_state.thermal.cells().keys())
        || initial.thermal.links() != final_state.thermal.links()
        || initial.thermal.cells().iter().any(|(id, before)| {
            let after = &final_state.thermal.cells()[id];
            before.mass_kg != after.mass_kg
                || before.material_id != after.material_id
                || before.geometry != after.geometry
                || before.boundary != after.boundary
                || before.body_id != after.body_id
                || before.segment_id != after.segment_id
        })
        || initial.chambers.iter().any(|(id, before)| {
            let after = &final_state.chambers[id];
            before.mass_kg() != after.mass_kg()
                || before.area_m2() != after.area_m2()
                || before.support() != after.support()
        })
    {
        return Err(Error::Balance(
            "open window must retain fixed physical owners and one forward clock",
        ));
    }
    let mut loads = BTreeMap::new();
    for drive in drives {
        check(cx)?;
        if !initial.chambers.contains_key(&drive.chamber)
            || loads.insert(drive.chamber, drive).is_some()
        {
            return Err(Error::Balance(
                "one held load per original open-window chamber required",
            ));
        }
    }
    if !initial.chambers.keys().eq(loads.keys()) {
        return Err(Error::Balance(
            "open window omitted an original chamber drive",
        ));
    }
    let mut accepted = Accepted {
        chambers: initial
            .chambers
            .keys()
            .map(|&id| (id, ChamberLedger::default()))
            .collect(),
        ..Accepted::default()
    };
    let mut cursor = initial.time_s();
    for step in steps {
        check(cx)?;
        if step.start_s != cursor
            || !step.end_s.is_finite()
            || step.end_s <= step.start_s
            || step.end_s > final_state.time_s()
        {
            return Err(Error::Balance(
                "accepted open-window receipts must cover one contiguous interval",
            ));
        }
        accepted.record(initial, step, reactions, cx)?;
        cursor = step.end_s;
    }
    if cursor != final_state.time_s() {
        return Err(Error::Balance(
            "accepted open-window receipts do not reach the final owner clock",
        ));
    }
    let mut thermal = Sum::default();
    for (id, before) in initial.thermal.cells() {
        check(cx)?;
        thermal.add(
            before.mass_kg
                * (final_state.thermal.cells()[id].specific_h_j_kg - before.specific_h_j_kg),
        )?;
    }
    let dh = thermal.total()?;
    let de = initial
        .circuit
        .fixed_storage_change_j(&final_state.circuit, ports)
        .map_err(Error::Circuit)?;
    let supplied = accepted.supply.total()?;
    let deposited = accepted.joule.total()?;
    let circuit_heat = accepted.circuit_heat.total()?;
    let external = accepted.external.total()?;
    let non_electrical = accepted.non_electrical.total()?;
    let cell_defect = accepted.cell_defect.total()?;
    let thermal_defect = accepted.thermal_defect.total()?;
    let circuit_defect = accepted.circuit_defect.total()?;
    let network_defect = accepted.network_defect.total()?;
    let network_residual = sum([de, dh, -supplied, -non_electrical])?;
    if cell_defect > controls.network.cell_energy_tolerance_j
        || thermal_defect > controls.network.thermal_energy_tolerance_j
        || circuit_defect > controls.network.circuit.energy_tolerance_j
        || network_defect > controls.network.total_energy_tolerance_j
        || sum([dh, -accepted.thermal_input.total()?])?.abs()
            > controls.network.thermal_energy_tolerance_j
        || sum([de, circuit_heat, -supplied])?.abs() > controls.network.circuit.energy_tolerance_j
        || sum([deposited, -circuit_heat])?.abs() > controls.network.total_energy_tolerance_j
        || network_residual.abs() > controls.network.total_energy_tolerance_j
    {
        return Err(Error::Balance(
            "open-window actual cell/circuit storage exceeds original fixed energy budgets",
        ));
    }
    let mut species_rows = Vec::new();
    let stream_entries = accepted.stream_entries().ok_or(Error::Balance(
        "open-window ledger entry count is unrepresentable",
    ))?;
    let mut global_species: BTreeMap<SpeciesId, Sum> = BTreeMap::new();
    let mut mass = Sum::default();
    let mut feedback_window = Sum::default();
    for value in [dh, de, -supplied, -external] {
        feedback_window.add(value)?;
    }
    let mut prescribed = Sum::default();
    let mut maximum_species = 0.0_f64;
    let mut maximum_energy = 0.0_f64;
    let mut chambers = Vec::with_capacity(initial.chambers.len());
    for (&id, before) in &initial.chambers {
        check(cx)?;
        let after = &final_state.chambers[&id];
        let drive = loads[&id];
        let ledger = &accepted.chambers[&id];
        let old = inventory(before.gas())?;
        let new = inventory(after.gas())?;
        let species: BTreeSet<_> = old
            .keys()
            .chain(new.keys())
            .chain(ledger.species.keys())
            .cloned()
            .collect();
        for species in species {
            check(cx)?;
            let initial_mass_kg = old.get(&species).copied().unwrap_or(0.0);
            let final_mass_kg = new.get(&species).copied().unwrap_or(0.0);
            let (internal_change_kg, reservoir_change_kg, reaction_change_kg) =
                if let Some(row) = ledger.species.get(&species) {
                    (
                        row.internal.total()?,
                        row.reservoir.total()?,
                        row.reaction.total()?,
                    )
                } else {
                    (0.0, 0.0, 0.0)
                };
            let change = finite(final_mass_kg - initial_mass_kg)?;
            let residual_kg = sum([
                change,
                -internal_change_kg,
                -reservoir_change_kg,
                -reaction_change_kg,
            ])?;
            maximum_species = maximum_species.max(residual_kg.abs());
            let global = global_species.entry(species.clone()).or_default();
            for value in [change, -reservoir_change_kg, -reaction_change_kg] {
                global.add(value)?;
            }
            if stream_entries
                .checked_add(species_rows.len())
                .is_none_or(|count| count >= MAX_SNAPSHOT_TRANSFERS)
            {
                return Err(Error::Balance(
                    "accepted open-window streams and species exceed native work bound",
                ));
            }
            species_rows.push(MolecularReservoirSpeciesBalance {
                chamber: id,
                species,
                initial_mass_kg,
                final_mass_kg,
                internal_change_kg,
                reservoir_change_kg,
                reaction_change_kg,
                residual_kg,
            });
        }
        let dm = finite(after.gas().mass_kg() - before.gas().mass_kg())?;
        mass.add(dm)?;
        let dk = finite(
            (0.5 * before.mass_kg() * (after.velocity_m_s() - before.velocity_m_s()))
                * (after.velocity_m_s() + before.velocity_m_s()),
        )?;
        let dv = finite(after.gas().volume_m3() - before.gas().volume_m3())?;
        let force_work = finite(drive.external_force_n * (dv / before.area_m2()))?;
        let pressure_work = finite(drive.exterior_pressure_pa * dv)?;
        let support_work = ledger.support_work.total()?;
        let support_impulse = ledger.support_impulse.total()?;
        let feedback_energy = ledger.feedback_energy.total()?;
        let heat = ledger.heat.total()?;
        let transport_residual = ledger.transport_residual.total()?;
        let feedback_residual = sum([
            feedback_energy,
            dk,
            pressure_work,
            -force_work,
            -support_work,
            -drive.heat_j,
            -ledger.evaluated_surface_heat.total()?,
        ])?;
        let energy_residual = sum([transport_residual, feedback_residual])?;
        let feedback_account_drift =
            sum([feedback_residual, -ledger.feedback_residual.total()?])?.abs();
        let transport_defect = ledger.transport_defect.total()?;
        let feedback_defect = sum([ledger.feedback_defect.total()?, feedback_account_drift])?;
        let energy_defect = sum([transport_defect, feedback_defect])?;
        let impulse = sum([
            before.mass_kg() * (after.velocity_m_s() - before.velocity_m_s()),
            -ledger.gas_impulse.total()?,
            duration * before.area_m2() * drive.exterior_pressure_pa,
            -duration * drive.external_force_n,
            -support_impulse,
        ])?;
        if sum([dm, -ledger.mass.total()?])?.abs()
            > molecular_controls.physical.exchange.mass_tolerance_kg
            || transport_defect > molecular_controls.physical.exchange.energy_tolerance_j
            || feedback_defect > drive.exchange.energy_tolerance_j
            || energy_defect > controls.total_energy_tolerance_j
            || energy_residual.abs() > controls.total_energy_tolerance_j
            || sum([
                feedback_energy,
                dk,
                pressure_work,
                -force_work,
                -support_work,
                -heat,
            ])?
            .abs()
                > drive.piston_solve.energy_tolerance_j
            || impulse.abs() > drive.piston_solve.momentum_tolerance_n_s
            || sum([dk, -ledger.kinetic.total()?])?.abs() > drive.piston_solve.energy_tolerance_j
            || sum([force_work, -ledger.force_work.total()?])?.abs()
                > drive.piston_solve.energy_tolerance_j
            || sum([pressure_work, -ledger.pressure_work.total()?])?.abs()
                > drive.piston_solve.energy_tolerance_j
        {
            return Err(Error::Balance(
                "open-window original caloric legs and actual wall storage/work exceed fixed budgets",
            ));
        }
        maximum_energy = maximum_energy.max(energy_residual.abs());
        for value in [
            feedback_energy,
            dk,
            pressure_work,
            -force_work,
            -support_work,
            -drive.heat_j,
        ] {
            feedback_window.add(value)?;
        }
        prescribed.add(drive.heat_j)?;
        chambers.push(MolecularReservoirWindowChamberBalance {
            chamber: id,
            mass_change_kg: dm,
            transport_energy_change_j: ledger.transport_energy.total()?,
            feedback_energy_change_j: feedback_energy,
            kinetic_energy_change_j: dk,
            entropy_change_j_k: ledger.entropy.total()?,
            gas_boundary_heat_j: heat,
            gas_boundary_entropy_j_k: ledger.heat_entropy.total()?,
            external_force_work_j: force_work,
            exterior_pressure_work_j: pressure_work,
            support_work_j: support_work,
            support_impulse_n_s: support_impulse,
            momentum_residual_n_s: impulse,
            transport_energy_residual_j: transport_residual,
            feedback_energy_residual_j: feedback_residual,
            energy_residual_j: energy_residual,
            energy_defect_bound_j: energy_defect,
            transported_species: ledger.moved,
            transported_heat: ledger.transported_heat,
        });
    }
    for row in global_species.values() {
        check(cx)?;
        maximum_species = maximum_species.max(row.total()?.abs());
    }
    let bath_mass = accepted.bath_mass.total()?;
    mass.add(-bath_mass)?;
    let mass_residual = mass.total()?;
    let transport_residual = accepted.transport_residual.total()?;
    let accumulated_feedback_residual = accepted.feedback_residual.total()?;
    let feedback_residual = feedback_window.total()?;
    let total_residual = sum([transport_residual, feedback_residual])?;
    let transport_defect = accepted.transport_defect.total()?;
    let feedback_defect = sum([
        accepted.feedback_defect.total()?,
        sum([feedback_residual, -accumulated_feedback_residual])?.abs(),
    ])?;
    let total_defect = sum([transport_defect, feedback_defect])?;
    let transport_entropy = accepted.transport_entropy.total()?;
    let entropy_deficit = accepted.transport_entropy_deficit.total()?;
    if mass_residual.abs() > molecular_controls.physical.exchange.mass_tolerance_kg
        || maximum_species > molecular_controls.physical.exchange.mass_tolerance_kg
        || transport_defect > molecular_controls.physical.exchange.energy_tolerance_j
        || total_defect > controls.total_energy_tolerance_j
        || total_residual.abs() > controls.total_energy_tolerance_j
        || entropy_deficit > molecular_controls.physical.entropy_tolerance_j_k
    {
        return Err(Error::Balance(
            "open-window accepted defects or endpoint species exceed one fixed physical budget",
        ));
    }
    accepted.bound()?;
    let transfers = accepted
        .transfers
        .into_iter()
        .map(|((donor, receiver, species), (mass, energy))| {
            check(cx)?;
            Ok(ParcelMolecularIntegratedTransfer {
                donor,
                receiver,
                species,
                mass_kg: mass.total()?,
                carried_enthalpy_j: energy.total()?,
            })
        })
        .collect::<Result<_, Error>>()?;
    let heat = accepted
        .heat
        .into_iter()
        .map(|(parcels, values)| {
            check(cx)?;
            Ok(ParcelMolecularIntegratedHeat {
                parcels,
                left_to_right_j: values.positive.total()?,
                right_to_left_j: values.negative.total()?,
            })
        })
        .collect::<Result<_, Error>>()?;
    let reservoir_species = accepted
        .reservoir_species
        .into_iter()
        .map(|((id, parcel, species), values)| {
            check(cx)?;
            Ok(ReservoirMolecularIntegratedSpecies {
                id,
                parcel,
                species,
                mass_into_gas_kg: values.positive.total()?,
                mass_out_of_gas_kg: values.negative.total()?,
                net_mass_into_gas_kg: values.net.total()?,
            })
        })
        .collect::<Result<_, Error>>()?;
    let reservoir_heat = accepted
        .reservoir_heat
        .into_iter()
        .map(|((id, parcel), values)| {
            check(cx)?;
            Ok(ReservoirMolecularIntegratedHeat {
                id,
                parcel,
                heat_into_gas_j: values.positive.total()?,
                heat_out_of_gas_j: values.negative.total()?,
                net_heat_into_gas_j: values.net.total()?,
            })
        })
        .collect::<Result<_, Error>>()?;
    let reaction_extents = accepted
        .extents
        .into_iter()
        .map(|((chamber, reaction), extent)| {
            check(cx)?;
            Ok(CoupledReactionExtent {
                chamber,
                reaction,
                extent_mol: extent.total()?,
            })
        })
        .collect::<Result<_, Error>>()?;
    check(cx)?;
    Ok(MolecularReservoirWindowBalance {
        thermal_energy_change_j: dh,
        electrical_energy_change_j: de,
        electrical_supply_work_j: supplied,
        deposited_joule_j: deposited,
        external_network_heat_j: external,
        prescribed_gas_heat_j: prescribed.total()?,
        non_electrical_heat_j: non_electrical,
        maximum_cell_residual_bound_j: cell_defect,
        thermal_energy_defect_bound_j: thermal_defect,
        circuit_energy_defect_bound_j: circuit_defect,
        network_energy_defect_bound_j: network_defect,
        network_energy_residual_j: network_residual,
        transport_energy_residual_j: transport_residual,
        feedback_energy_residual_j: feedback_residual,
        accumulated_feedback_energy_residual_j: accumulated_feedback_residual,
        total_energy_residual_j: total_residual,
        transport_energy_defect_bound_j: transport_defect,
        feedback_energy_defect_bound_j: feedback_defect,
        total_energy_defect_bound_j: total_defect,
        maximum_chamber_energy_residual_j: maximum_energy,
        net_reservoir_mass_kg: bath_mass,
        total_mass_residual_kg: mass_residual,
        maximum_species_mass_residual_kg: maximum_species,
        transport_entropy_production_j_k: transport_entropy,
        transport_entropy_deficit_bound_j_k: entropy_deficit,
        species: species_rows,
        transfers,
        heat,
        reservoir_species,
        reservoir_heat,
        reaction_extents,
        chambers,
    })
}

pub(super) fn macro_ports(
    initial: &MaterialExchangeState<'_>,
    next: &mut MaterialExchangeState<'_>,
    cx: &Cx<'_>,
    controls: NetworkSolve,
    balance: &MolecularReservoirWindowBalance,
) -> Result<(), Error> {
    let duration = finite(next.time_s() - initial.time_s())?;
    for account in &balance.chambers {
        check(cx)?;
        let old = initial
            .chambers
            .get(&account.chamber)
            .ok_or(Error::Balance(
                "open-window history names an unknown initial chamber",
            ))?;
        let after = next
            .chambers
            .get_mut(&account.chamber)
            .ok_or(Error::Balance(
                "open-window history names an unknown final chamber",
            ))?;
        after.ports = old.ports.clone();
        after.ports.accept(StepAccount {
            tick: after.ports.next_tick(), duration_s: duration,
            storage_change_j: Some(sum([account.transport_energy_change_j,
                account.feedback_energy_change_j, account.kinetic_energy_change_j])?),
            loss_j: Some(0.0), included_work_j: 0.0,
            entropy_change_j_k: account.entropy_change_j_k,
            tolerance_j: Some(controls.total_energy_tolerance_j),
            unavailable: (account.transported_species || account.transported_heat).then_some(
                "held-bath species/heat streams require the whole-window open molecular balance"),
        }, [
            PortTransfer::work(sum([account.external_force_work_j, -account.exterior_pressure_work_j,
                account.support_work_j])?, true),
            PortTransfer::heat_profile(account.gas_boundary_heat_j, account.gas_boundary_entropy_j_k),
        ]);
    }
    check(cx)
}
