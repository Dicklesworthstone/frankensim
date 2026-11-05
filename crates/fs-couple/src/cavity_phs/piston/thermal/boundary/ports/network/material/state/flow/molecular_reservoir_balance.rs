//! Open streams retain source-relative energy and independent endpoint species.
use super::*;
use crate::port_adoption::{PortTransfer, StepAccount};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_thermochem::SpeciesId;

/// Actual endpoint inventory minus accepted internal, bath and chemical changes.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularReservoirSpeciesBalance {
    pub chamber: u64,
    pub species: SpeciesId,
    pub initial_mass_kg: f64,
    pub final_mass_kg: f64,
    pub internal_change_kg: f64,
    pub reservoir_change_kg: f64,
    pub reaction_change_kg: f64,
    pub residual_kg: f64,
}

/// Two original caloric legs retained separately so large transported formation
/// energy does not erase a small subsequent thermal or mechanical increment.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularReservoirChamberBalance {
    pub chamber: u64,
    pub mass_change_kg: f64,
    /// Intended source-reference species change plus actual native sensible
    /// storage; inventory rounding has its separate species/mass budget.
    pub transport_energy_change_j: f64,
    /// Original reaction/heat/motion gas-energy increment and caloric conventions.
    pub feedback_energy_change_j: f64,
    pub kinetic_energy_change_j: f64,
    pub entropy_change_j_k: f64,
    /// Sum of native reference-cancelled transport and coupled boundary defects.
    pub energy_residual_j: f64,
}

/// One fixed absolute budget over both independently checked physical stages.
/// The energy residual composes actual native transport/feedback legs; it does
/// not subtract gross absolute U and reservoir enthalpy with large source datums.
/// The native flow.audit retains bath-inclusive transport entropy production;
/// no whole-assembly entropy-production certificate is inferred from that stage.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularReservoirBalance {
    pub net_reservoir_mass_kg: f64,
    pub transport_energy_residual_j: f64,
    pub feedback_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub maximum_chamber_energy_residual_j: f64,
    pub total_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub species: Vec<MolecularReservoirSpeciesBalance>,
    pub chambers: Vec<MolecularReservoirChamberBalance>,
}

type Error = FlowBalanceError;
fn checkpoint(cx: &Cx<'_>) -> Result<(), Error> {
    cx.checkpoint().map_err(|_| Error::Cancelled)
}
fn finite(value: f64) -> Result<f64, Error> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::Invalid(
            "open molecular assembly arithmetic is unrepresentable",
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
#[derive(Default)]
struct SpeciesLedger {
    internal: Sum,
    reservoir: Sum,
    reaction: Sum,
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
        .map(|(species, &fraction)| {
            let mass = finite((amount * fraction) * species.molar_mass().value())?;
            if mass < 0.0 {
                return Err(Error::Invalid("negative actual species inventory"));
            }
            Ok((species.species().clone(), mass))
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
pub(super) fn finish<I>(
    initial: &MaterialExchangeState<'_>,
    staged: &BTreeMap<u64, GasPiston>,
    next: &mut MaterialExchangeState<'_>,
    cx: &Cx<'_>,
    flow: &ReservoirMolecularStep,
    feedback: &MaterialExchangeAdvance<I>,
    reactions: &[ChamberReaction<'_>],
    molecular_controls: ParcelMolecularSolve,
    controls: NetworkSolve,
) -> Result<MolecularReservoirBalance, Error> {
    checkpoint(cx)?;
    if !initial.chambers.keys().eq(staged.keys())
        || !initial.chambers.keys().eq(next.chambers.keys())
    {
        return Err(Error::Invalid(
            "open transport and feedback must retain every chamber identity",
        ));
    }
    let native: BTreeMap<_, _> = flow
        .exchange
        .parcels
        .iter()
        .map(|row| (row.id, row))
        .collect();
    let boundaries: BTreeMap<_, _> = feedback
        .report
        .source
        .chambers
        .iter()
        .map(|row| (row.id, row))
        .collect();
    if !initial.chambers.keys().eq(native.keys()) || !initial.chambers.keys().eq(boundaries.keys())
    {
        return Err(Error::Invalid(
            "missing native open-assembly chamber account",
        ));
    }
    let mut ledgers: BTreeMap<(u64, SpeciesId), SpeciesLedger> = BTreeMap::new();
    let mut moved = BTreeSet::new();
    let mut bath_mass = Sum::default();
    for row in &flow.exchange.fluxes {
        checkpoint(cx)?;
        let t = &row.transfer;
        for (id, sign) in [(t.donor, -1.0), (t.receiver, 1.0)] {
            ledgers
                .entry((id, t.species.clone()))
                .or_default()
                .internal
                .add(sign * t.mass_kg)?;
            if t.mass_kg > 0.0 {
                moved.insert(id);
            }
        }
    }
    for row in &flow.exchange.boundaries {
        checkpoint(cx)?;
        ledgers
            .entry((row.parcel, row.transfer.species.clone()))
            .or_default()
            .reservoir
            .add(row.transfer.net_mass_into_gas_kg)?;
        bath_mass.add(row.transfer.net_mass_into_gas_kg)?;
        if row.transfer.incoming_mass_kg > 0.0 || row.transfer.outgoing_mass_kg > 0.0 {
            moved.insert(row.parcel);
        }
    }
    for row in &flow.internal_links {
        checkpoint(cx)?;
        if row.heat_correction_left_to_right_j != 0.0 {
            moved.extend(row.parcels);
        }
    }
    for row in &flow.exchange.heat {
        checkpoint(cx)?;
        if row.heat_left_to_right_j != 0.0 {
            moved.extend(row.parcels);
        }
    }
    for row in &flow.exchange.external_heat {
        checkpoint(cx)?;
        if row.heat_into_gas_j != 0.0 {
            moved.insert(row.parcel);
        }
    }
    for row in &feedback.report.source.reactions {
        checkpoint(cx)?;
        let binding = reactions
            .iter()
            .find(|r| r.chamber == row.chamber)
            .ok_or(Error::Invalid("unbound accepted open-assembly reaction"))?;
        let sources = binding.reaction.catalog().thermodynamics().components();
        for species in &row.reaction.species {
            checkpoint(cx)?;
            let source = sources
                .iter()
                .find(|s| s.species() == &species.species)
                .ok_or(Error::Invalid(
                    "reaction omitted its original caloric species",
                ))?;
            ledgers
                .entry((row.chamber, species.species.clone()))
                .or_default()
                .reaction
                .add(species.stoichiometric_change_mol * source.molar_mass().value())?;
        }
    }
    let mut species_rows = Vec::new();
    let mut global_species: BTreeMap<SpeciesId, Sum> = BTreeMap::new();
    let mut mass = Sum::default();
    let mut maximum_species = 0.0_f64;
    let mut maximum_energy = 0.0_f64;
    let mut chambers = Vec::with_capacity(initial.chambers.len());
    let duration_s = next.time_s() - initial.time_s();
    for (&id, before) in &initial.chambers {
        checkpoint(cx)?;
        let after = &next.chambers[&id];
        let old = inventory(before.gas())?;
        let new = inventory(after.gas())?;
        let species: BTreeSet<_> = old
            .keys()
            .chain(new.keys())
            .cloned()
            .chain(
                ledgers
                    .keys()
                    .filter(|(chamber, _)| *chamber == id)
                    .map(|(_, species)| species.clone()),
            )
            .collect();
        for species in species {
            checkpoint(cx)?;
            let row = ledgers.entry((id, species.clone())).or_default();
            let initial_mass_kg = old.get(&species).copied().unwrap_or(0.0);
            let final_mass_kg = new.get(&species).copied().unwrap_or(0.0);
            let internal_change_kg = row.internal.total()?;
            let reservoir_change_kg = row.reservoir.total()?;
            let reaction_change_kg = row.reaction.total()?;
            let mut residual = Sum::default();
            for value in [
                final_mass_kg - initial_mass_kg,
                -internal_change_kg,
                -reservoir_change_kg,
                -reaction_change_kg,
            ] {
                residual.add(value)?;
            }
            let residual_kg = residual.total()?;
            maximum_species = maximum_species.max(residual_kg.abs());
            let global = global_species.entry(species.clone()).or_default();
            for value in [
                final_mass_kg - initial_mass_kg,
                -reservoir_change_kg,
                -reaction_change_kg,
            ] {
                global.add(value)?;
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
        let boundary = &boundaries[&id].boundary;
        let chemical = feedback
            .report
            .source
            .reactions
            .iter()
            .find(|r| r.chamber == id);
        let feedback_energy = chemical.map_or(boundary.piston.gas.stored_energy_change_j, |r| {
            r.gas_energy_change_j
        });
        let feedback_residual = chemical.map_or(boundary.total_energy_residual_j, |r| {
            r.boundary_energy_residual_j
        });
        let native = native[&id];
        let energy_residual_j = finite(native.energy_residual_j + feedback_residual)?;
        maximum_energy = maximum_energy.max(energy_residual_j.abs());
        let transport_energy_change_j = finite(native.energy_change_j)?;
        let feedback_energy_change_j = finite(feedback_energy)?;
        let kinetic_energy_change_j = finite(boundary.piston.kinetic_energy_change_j)?;
        // Entropy is retained as a storage observation. The open transport's
        // source-relative entropy production remains in the original audit.
        let transported = staged[&id].gas();
        let old_s = before
            .gas()
            .state()
            .specific_entropy_j_kg_k()
            .map_err(|error| Error::Gas { chamber: id, error })?;
        let staged_s = transported
            .state()
            .specific_entropy_j_kg_k()
            .map_err(|error| Error::Gas { chamber: id, error })?;
        let transport_entropy = transported.mass_kg() * (staged_s - old_s)
            + (transported.mass_kg() - before.gas().mass_kg()) * old_s;
        let entropy_change_j_k = finite(
            transport_entropy
                + boundary.piston.gas.entropy_change_j_k
                + chemical.map_or(0.0, |r| r.reaction.entropy_change_j_k),
        )?;
        if moved.contains(&id) {
            let after = next
                .chambers
                .get_mut(&id)
                .ok_or(Error::Invalid("missing final open chamber"))?;
            after.ports = before.ports.clone();
            after.ports.accept(StepAccount {
                tick: after.ports.next_tick(), duration_s,
                storage_change_j: Some(finite(transport_energy_change_j
                    + feedback_energy_change_j + kinetic_energy_change_j)?),
                loss_j: Some(0.0), included_work_j: 0.0,
                entropy_change_j_k, tolerance_j: Some(controls.total_energy_tolerance_j),
                unavailable: Some("held-bath species/heat streams require the enclosing open molecular balance"),
            }, [
                PortTransfer::work(boundary.piston.external_force_work_j
                    - boundary.piston.exterior_pressure_work_j, true),
                PortTransfer::heat(boundary.piston.gas.supplied_heat_j,
                    boundary.piston.midpoint.state().temperature),
            ]);
        }
        chambers.push(MolecularReservoirChamberBalance {
            chamber: id,
            mass_change_kg: dm,
            transport_energy_change_j,
            feedback_energy_change_j,
            kinetic_energy_change_j,
            entropy_change_j_k,
            energy_residual_j,
        });
    }
    for row in global_species.values() {
        maximum_species = maximum_species.max(row.total()?.abs());
    }
    let net_reservoir_mass_kg = bath_mass.total()?;
    mass.add(-net_reservoir_mass_kg)?;
    let total_mass_residual_kg = mass.total()?;
    let transport_energy_residual_j = finite(flow.audit.energy_residual_j)?;
    let feedback_energy_residual_j = finite(feedback.step.total_energy_residual_j)?;
    let total_energy_residual_j = finite(transport_energy_residual_j + feedback_energy_residual_j)?;
    if total_mass_residual_kg.abs() > molecular_controls.physical.exchange.mass_tolerance_kg
        || maximum_species > molecular_controls.physical.exchange.mass_tolerance_kg
        || maximum_energy > controls.total_energy_tolerance_j
        || total_energy_residual_j.abs() > controls.total_energy_tolerance_j
    {
        return Err(Error::Invalid(
            "open molecular stages exceed one fixed species, mass or energy budget",
        ));
    }
    checkpoint(cx)?;
    Ok(MolecularReservoirBalance {
        net_reservoir_mass_kg,
        transport_energy_residual_j,
        feedback_energy_residual_j,
        total_energy_residual_j,
        maximum_chamber_energy_residual_j: maximum_energy,
        total_mass_residual_kg,
        maximum_species_mass_residual_kg: maximum_species,
        species: species_rows,
        chambers,
    })
}
