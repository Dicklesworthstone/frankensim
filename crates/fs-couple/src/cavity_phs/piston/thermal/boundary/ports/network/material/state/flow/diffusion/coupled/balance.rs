//! Accepted native receipts plus independent original-to-final storage.
use super::*;
use super::PhysicalControlError as Error;
use fs_thermochem::ReactionId;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::diffusion::network::adaptive::DiffusionNetworkIntegratedTransfer;
use crate::port_adoption::{PortTransfer, StepAccount};

#[derive(Clone, Debug, PartialEq)]
pub struct CoupledDiffusionSpeciesBalance {
    pub chamber: u64,
    pub species: SpeciesId,
    pub initial_mass_kg: f64,
    pub final_mass_kg: f64,
    /// Accepted signed molecular transfer into this chamber [kg].
    pub diffusive_change_kg: f64,
    /// Original stoichiometric changes converted by original species molar masses.
    pub reaction_change_kg: f64,
    pub residual_kg: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledReactionExtent {
    pub chamber: u64,
    pub reaction: ReactionId,
    pub extent_mol: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledDiffusionChamberBalance {
    pub chamber: u64,
    pub mass_change_kg: f64,
    pub gas_energy_change_j: f64,
    pub kinetic_energy_change_j: f64,
    pub entropy_change_j_k: f64,
    pub net_diffusive_enthalpy_j: f64,
    /// Additional paired gas/gas heat, separate from transported partial enthalpy.
    pub net_diffusive_heat_j: f64,
    pub gas_boundary_heat_j: f64,
    pub gas_boundary_entropy_j_k: f64,
    pub external_force_work_j: f64,
    pub exterior_pressure_work_j: f64,
    pub support_impulse_n_s: f64,
    pub momentum_residual_n_s: f64,
    pub energy_residual_j: f64,
    /// Gross accepted stream, even if final net mass or species change is zero.
    pub transported_species: bool,
    /// A separate gas/gas heat path was active, even if its net integral is zero.
    pub transported_heat: bool,
}
/// All budgets retain their original absolute physical allowances across the
/// complete window. Neither a timestep count nor a rejected probe enlarges them.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledDiffusionBalance {
    pub thermal_energy_change_j: f64,
    pub electrical_energy_change_j: f64,
    pub electrical_supply_work_j: f64,
    pub deposited_joule_j: f64,
    pub external_network_heat_j: f64,
    pub prescribed_gas_heat_j: f64,
    pub non_electrical_heat_j: f64,
    /// Sum of accepted native per-step maximum cell residuals. This bounds
    /// every fixed cell's accumulated defect even when cell errors cancel.
    pub maximum_cell_residual_bound_j: f64,
    pub network_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub total_mass_residual_kg: f64,
    /// Maximum of per-chamber and whole-network species residuals after chemistry.
    pub maximum_species_mass_residual_kg: f64,
    pub species: Vec<CoupledDiffusionSpeciesBalance>,
    /// Canonical donor/receiver/species rows. Reverse transport is a separate row.
    pub transfers: Vec<DiffusionNetworkIntegratedTransfer>,
    pub reaction_extents: Vec<CoupledReactionExtent>,
    pub chambers: Vec<CoupledDiffusionChamberBalance>,
}

#[derive(Clone, Copy, Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), Error> {
        let value = finite(value)?;
        let next = finite(self.value + value)?;
        self.correction = finite(self.correction + if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else { (value - next) + self.value })?;
        self.value = next;
        Ok(())
    }
    fn total(&self) -> Result<f64, Error> { finite(self.value + self.correction) }
}
#[derive(Default)]
struct SpeciesLedger { transport: Sum, chemical: Sum }
#[derive(Default)]
struct ChamberLedger {
    enthalpy: Sum, transport_heat: Sum, heat: Sum, heat_entropy: Sum,
    force_work: Sum, pressure_work: Sum, gas_impulse: Sum, support_impulse: Sum,
    moved: bool, transported_heat: bool,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn summarize<F: TransportLedger, I>(
    initial: &MaterialExchangeState<'_>, final_state: &MaterialExchangeState<'_>, cx: &Cx<'_>,
    ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
    species_controls: SpeciesNetworkSolve, controls: NetworkSolve, steps: &[CoupledDiffusionSubstep<F, I>],
) -> Result<CoupledDiffusionBalance, Error> {
    check(cx)?;
    let mut ledgers: BTreeMap<(u64, SpeciesId), SpeciesLedger> = BTreeMap::new();
    let mut transfers: BTreeMap<(u64, u64, SpeciesId), (Sum, Sum)> = BTreeMap::new();
    let mut extents: BTreeMap<(u64, ReactionId), Sum> = BTreeMap::new();
    let mut chamber_ledgers: BTreeMap<u64, ChamberLedger> = initial.chambers.keys()
        .map(|&id| (id, ChamberLedger::default())).collect();
    let (mut supply, mut joule, mut circuit_heat, mut external, mut non_electrical, mut thermal_input) =
        (Sum::default(), Sum::default(), Sum::default(), Sum::default(), Sum::default(), Sum::default());
    let mut cell_residual_bound = Sum::default();
    for step in steps {
        check(cx)?;
        let r = &step.advance.exchange.report;
        cell_residual_bound.add(r.maximum_cell_residual_j)?;
        let exchange = step.advance.flow.species_exchange();
        if !exchange.boundaries.is_empty() || !exchange.external_heat.is_empty() {
            return Err(Error::Balance("external transport streams need their additional finite owner in the full window balance"));
        }
        supply.add(r.electrical.native.supplied)?;
        joule.add(r.deposited_joule_j)?;
        circuit_heat.add(r.electrical.heat_j)?;
        external.add(step.advance.exchange.step.external_network_heat_j)?;
        non_electrical.add(r.non_electrical_heat_j)?;
        thermal_input.add(r.thermal.external_heat_j)?;
        thermal_input.add(r.thermal.incoming_energy_j)?;
        for flux in &exchange.fluxes {
            check(cx)?;
            let t = &flux.transfer;
            let gross = transfers.entry((t.donor, t.receiver, t.species.clone())).or_default();
            gross.0.add(t.mass_kg)?; gross.1.add(flux.carried_enthalpy_j)?;
            for (id, sign) in [(t.donor, -1.0), (t.receiver, 1.0)] {
                ledgers.entry((id, t.species.clone())).or_default().transport.add(sign * t.mass_kg)?;
                chamber_ledgers.get_mut(&id).ok_or(Error::Balance("unknown accepted diffusion endpoint"))?
                    .moved |= t.mass_kg > 0.0;
            }
        }
        for parcel in &exchange.parcels {
            let ledger = chamber_ledgers.get_mut(&parcel.id).ok_or(Error::Balance("missing diffusion chamber"))?;
            ledger.enthalpy.add(parcel.net_enthalpy_j)?;
            ledger.transport_heat.add(parcel.net_additional_heat_j)?;
        }
        for heat in &exchange.heat {
            check(cx)?;
            if heat.heat_left_to_right_j != 0.0 {
                for id in heat.parcels {
                    chamber_ledgers.get_mut(&id).ok_or(Error::Balance("unknown paired transport-heat endpoint"))?
                        .transported_heat = true;
                }
            }
        }
        for chemical in &r.source.reactions {
            check(cx)?;
            let binding = reactions.iter().find(|r| r.chamber == chemical.chamber)
                .ok_or(Error::Balance("unbound accepted reaction"))?;
            let catalog = binding.reaction.catalog().thermodynamics().components();
            for species in &chemical.reaction.species {
                let original = catalog.iter().find(|s| s.species() == &species.species)
                    .ok_or(Error::Balance("accepted reaction omitted its original species source"))?;
                ledgers.entry((chemical.chamber, species.species.clone())).or_default().chemical
                    .add(species.stoichiometric_change_mol * original.molar_mass().value())?;
            }
            for (id, extent) in chemical.reaction.reaction_extents() {
                extents.entry((chemical.chamber, id.clone())).or_default().add(extent)?;
            }
        }
        for chamber in &r.source.chambers {
            check(cx)?;
            let sums = chamber_ledgers.get_mut(&chamber.id).ok_or(Error::Balance("unknown boundary chamber"))?;
            let p = &chamber.boundary.piston;
            let dt = step.end_s - step.start_s;
            sums.heat.add(p.gas.supplied_heat_j)?;
            sums.heat_entropy.add(p.gas.supplied_heat_j / p.midpoint.state().temperature)?;
            sums.force_work.add(p.external_force_work_j)?;
            sums.pressure_work.add(p.exterior_pressure_work_j)?;
            sums.gas_impulse.add(dt * initial.chambers[&chamber.id].area_m2() * p.midpoint.state().pressure)?;
            sums.support_impulse.add(p.support_impulse_n_s)?;
        }
    }
    let mut thermal = Sum::default();
    for (id, before) in initial.thermal.cells() {
        check(cx)?;
        thermal.add(before.mass_kg * (final_state.thermal.cells()[id].specific_h_j_kg - before.specific_h_j_kg))?;
    }
    let dh = thermal.total()?;
    let de = initial.circuit.fixed_storage_change_j(&final_state.circuit, ports).map_err(Error::Circuit)?;
    let supplied = supply.total()?; let deposited = joule.total()?;
    let outside_heat = external.total()?; let non_electrical = non_electrical.total()?;
    let cell_residual_bound = cell_residual_bound.total()?;
    let network_residual = finite((de + dh) - supplied - non_electrical)?;
    if cell_residual_bound > controls.network.cell_energy_tolerance_j
        || (dh - thermal_input.total()?).abs() > controls.network.thermal_energy_tolerance_j
        || (de + circuit_heat.total()? - supplied).abs() > controls.network.circuit.energy_tolerance_j
        || (deposited - circuit_heat.total()?).abs() > controls.network.total_energy_tolerance_j
        || network_residual.abs() > controls.network.total_energy_tolerance_j
    {
        return Err(Error::Balance("whole-window native cell/circuit storage exceeds original energy allowances"));
    }
    let mut mass = Sum::default();
    let mut total = Sum::default();
    for value in [dh, de, -supplied, -outside_heat] { total.add(value)?; }
    let mut prescribed_heat = Sum::default();
    let mut species_rows = Vec::new();
    let mut global_species: BTreeMap<SpeciesId, Sum> = BTreeMap::new();
    let mut maximum_species_residual = 0.0_f64;
    let mut chambers = Vec::new();
    let duration = final_state.time_s() - initial.time_s();
    for (&id, before) in &initial.chambers {
        check(cx)?;
        let after = &final_state.chambers[&id];
        let drive = drives.iter().find(|d| d.chamber == id)
            .ok_or(Error::Balance("missing admitted chamber drive"))?;
        let ledger = &chamber_ledgers[&id];
        let old = inventory(before.gas())?; let new = inventory(after.gas())?;
        let species: BTreeSet<_> = old.keys().chain(new.keys()).cloned()
            .chain(ledgers.keys().filter(|(chamber, _)| *chamber == id).map(|(_, species)| species.clone())).collect();
        for species in species {
            check(cx)?;
            let entry = ledgers.entry((id, species.clone())).or_default();
            let initial_mass_kg = old.get(&species).copied().unwrap_or(0.0);
            let final_mass_kg = new.get(&species).copied().unwrap_or(0.0);
            let diffusive_change_kg = entry.transport.total()?;
            let reaction_change_kg = entry.chemical.total()?;
            let change = final_mass_kg - initial_mass_kg;
            let residual_kg = finite((change - diffusive_change_kg) - reaction_change_kg)?;
            maximum_species_residual = maximum_species_residual.max(residual_kg.abs());
            global_species.entry(species.clone()).or_default().add(change - reaction_change_kg)?;
            species_rows.push(CoupledDiffusionSpeciesBalance { chamber: id, species,
                initial_mass_kg, final_mass_kg, diffusive_change_kg, reaction_change_kg, residual_kg });
        }
        let dm = finite(after.gas().mass_kg() - before.gas().mass_kg())?;
        mass.add(dm)?;
        // Same extensive source account as the original full assembly, with
        // before/after measured over this entire requested window.
        let u0 = before.gas().state().specific_internal_energy_j_kg();
        let du = finite(after.gas().mass_kg() * (after.gas().state().specific_internal_energy_j_kg() - u0) + dm * u0)?;
        let s0 = before.gas().state().specific_entropy_j_kg_k().map_err(|error| Error::Gas { chamber: id, error })?;
        let s1 = after.gas().state().specific_entropy_j_kg_k().map_err(|error| Error::Gas { chamber: id, error })?;
        let ds = finite(after.gas().mass_kg() * (s1 - s0) + dm * s0)?;
        let dk = finite((0.5 * before.mass_kg() * (after.velocity_m_s() - before.velocity_m_s()))
            * (after.velocity_m_s() + before.velocity_m_s()))?;
        let dv = finite(after.gas().volume_m3() - before.gas().volume_m3())?;
        let force_work = finite(drive.external_force_n * (dv / before.area_m2()))?;
        let pressure_work = finite(drive.exterior_pressure_pa * dv)?;
        let heat = ledger.heat.total()?; let enthalpy = ledger.enthalpy.total()?;
        let transport_heat = ledger.transport_heat.total()?;
        let energy = finite((du + dk) + pressure_work - force_work - heat - enthalpy - transport_heat)?;
        let impulse = finite(before.mass_kg() * (after.velocity_m_s() - before.velocity_m_s())
            - ledger.gas_impulse.total()? + duration * before.area_m2() * drive.exterior_pressure_pa
            - duration * drive.external_force_n - ledger.support_impulse.total()?)?;
        if energy.abs() > controls.total_energy_tolerance_j
            || impulse.abs() > drive.piston_solve.momentum_tolerance_n_s
            || (force_work - ledger.force_work.total()?).abs() > drive.piston_solve.energy_tolerance_j
            || (pressure_work - ledger.pressure_work.total()?).abs() > drive.piston_solve.energy_tolerance_j
        {
            return Err(Error::Balance("whole-window chamber streams, energy or momentum exceed original allowances"));
        }
        for value in [du, dk, pressure_work, -force_work, -drive.heat_j] { total.add(value)?; }
        prescribed_heat.add(drive.heat_j)?;
        chambers.push(CoupledDiffusionChamberBalance { chamber: id, mass_change_kg: dm,
            gas_energy_change_j: du, kinetic_energy_change_j: dk, entropy_change_j_k: ds,
            net_diffusive_enthalpy_j: enthalpy, net_diffusive_heat_j: transport_heat, gas_boundary_heat_j: heat,
            gas_boundary_entropy_j_k: ledger.heat_entropy.total()?, external_force_work_j: force_work,
            exterior_pressure_work_j: pressure_work, support_impulse_n_s: ledger.support_impulse.total()?,
            momentum_residual_n_s: impulse, energy_residual_j: energy,
            transported_species: ledger.moved, transported_heat: ledger.transported_heat });
    }
    for value in global_species.values() {
        maximum_species_residual = maximum_species_residual.max(value.total()?.abs());
    }
    let mass = mass.total()?; let energy = total.total()?;
    if mass.abs() > species_controls.mass_tolerance_kg
        || maximum_species_residual > species_controls.mass_tolerance_kg
        || energy.abs() > controls.total_energy_tolerance_j
    {
        return Err(Error::Balance("whole-window actual species, mass or energy exceeds original absolute budgets"));
    }
    let transfers = transfers.into_iter().map(|((donor, receiver, species), (mass, heat))| {
        Ok(DiffusionNetworkIntegratedTransfer { donor, receiver, species,
            mass_kg: mass.total()?, carried_enthalpy_j: heat.total()? })
    }).collect::<Result<_, Error>>()?;
    let reaction_extents = extents.into_iter().map(|((chamber, reaction), extent)| {
        Ok(CoupledReactionExtent { chamber, reaction, extent_mol: extent.total()? })
    }).collect::<Result<_, Error>>()?;
    check(cx)?;
    Ok(CoupledDiffusionBalance { thermal_energy_change_j: dh, electrical_energy_change_j: de,
        electrical_supply_work_j: supplied, deposited_joule_j: deposited,
        external_network_heat_j: outside_heat, prescribed_gas_heat_j: prescribed_heat.total()?,
        non_electrical_heat_j: non_electrical, network_energy_residual_j: network_residual,
        maximum_cell_residual_bound_j: cell_residual_bound,
        total_energy_residual_j: energy, total_mass_residual_kg: mass,
        maximum_species_mass_residual_kg: maximum_species_residual,
        species: species_rows, transfers, reaction_extents, chambers })
}

pub(super) fn macro_ports(initial: &MaterialExchangeState<'_>, next: &mut MaterialExchangeState<'_>,
    cx: &Cx<'_>, controls: NetworkSolve, balance: &CoupledDiffusionBalance,
) -> Result<(), Error> {
    let duration = next.time_s() - initial.time_s();
    for account in &balance.chambers {
        check(cx)?;
        let old = &initial.chambers[&account.chamber];
        let after = next.chambers.get_mut(&account.chamber).ok_or(Error::Balance("missing accepted chamber"))?;
        after.ports = old.ports.clone();
        after.ports.accept(StepAccount {
            tick: after.ports.next_tick(), duration_s: duration,
            storage_change_j: Some(finite(account.gas_energy_change_j + account.kinetic_energy_change_j)?),
            loss_j: Some(0.0), included_work_j: 0.0, entropy_change_j_k: account.entropy_change_j_k,
            tolerance_j: Some(controls.total_energy_tolerance_j),
            unavailable: (account.transported_species || account.transported_heat).then_some(
                "gas transport needs stream/heat ports; use the whole-window coupled transport balance"),
        }, [PortTransfer::work(account.external_force_work_j - account.exterior_pressure_work_j, true),
            PortTransfer::heat_profile(account.gas_boundary_heat_j, account.gas_boundary_entropy_j_k)]);
    }
    Ok(())
}
