//! Reserved future-neutron energy on the original delayed-family history.
//!
//! A precursor is a count in DiffusionState, never a second state or solver.
//! The reserve transfers internally at its native decay and travels with its
//! original material flow. Exterior precursor energy is not neutron leakage.

use super::{
    Accumulator, DiffusionModel, DiffusionReactionKind, DiffusionState, DiffusionStep,
    FissionError, InducedFissionPrecursorFunding, NeutronKineticError, NeutronKineticModel,
    checked, nonnegative, positive, sum, weighted_sum,
};

pub(super) fn energies(
    diffusion: &DiffusionModel,
    neutron_energies: &[f64],
) -> Result<Vec<f64>, NeutronKineticError> {
    diffusion
        .delayed_groups()
        .iter()
        .map(|family| {
            let energy = weighted_sum(&family.spectrum, neutron_energies)?;
            if !positive(energy) {
                return Err(NeutronKineticError::Arithmetic(
                    "precursor reserve must be representable and positive",
                ));
            }
            Ok(energy)
        })
        .collect()
}

/// Retain the actual future-neutron coefficient; a declared cap is an admission
/// bound, not a source value. Every reaction is checked even when inactive.
pub(super) fn funding(
    diffusion: &DiffusionModel,
    energies: &[f64],
    declarations: &[Option<InducedFissionPrecursorFunding>],
) -> Result<Vec<f64>, NeutronKineticError> {
    let mut reserved = Vec::with_capacity(diffusion.reactions().len());
    for (reaction, declaration) in diffusion.reactions().iter().zip(declarations) {
        let DiffusionReactionKind::Fission { delayed_yields, .. } = &reaction.kind else {
            if declaration.is_some() {
                return Err(NeutronKineticError::Invalid(
                    "induced precursor funding must follow original fission reactions only",
                ));
            }
            reserved.push(0.0);
            continue;
        };
        let mut total = Accumulator::new(0.0);
        for (&yield_, &energy) in delayed_yields.iter().zip(energies) {
            let value = checked(yield_ * energy)?;
            if yield_ > 0.0 && value == 0.0 {
                return Err(NeutronKineticError::Arithmetic(
                    "induced-fission precursor energy underflow",
                ));
            }
            total.add(value)?;
        }
        match declaration {
            Some(InducedFissionPrecursorFunding::AdditionalDelayedConversion {
                maximum_energy_j_per_fission,
            }) if nonnegative(*maximum_energy_j_per_fission)
                && total.value() <= *maximum_energy_j_per_fission => {}
            None if delayed_yields.iter().all(|&yield_| yield_ == 0.0) => {}
            _ => {
                return Err(NeutronKineticError::Invalid(
                    "induced-fission future-neutron energy requires its own finite sufficient funding cap",
                ));
            }
        }
        reserved.push(total.value());
    }
    Ok(reserved)
}

pub(super) struct Accounts {
    pub initial_j: f64,
    pub final_j: f64,
    pub injected_j: f64,
    pub converted_j: f64,
    pub cell_induced_j: Vec<f64>,
    pub induced_j: f64,
    pub cell_exported_j: Vec<f64>,
    pub exported_j: f64,
}

/// Audit family storage separately, then remove its actual delayed-neutron
/// births from each neutron balance. This prevents opposite storage defects
/// from cancelling while retaining one common absolute energy allowance.
pub(super) fn account(
    model: &NeutronKineticModel,
    before: &DiffusionState,
    after: &DiffusionState,
    step: &DiffusionStep,
    neutron_balance: &mut [Accumulator],
    defect: &mut Accumulator,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<Accounts, NeutronKineticError> {
    let families = model.precursor_energies_j.len();
    let cells = neutron_balance.len();
    let nodes = cells * families;
    let transport = after.model().precursor_transport();
    if after.model().delayed_groups() != model.diffusion.delayed_groups()
        || before.model().delayed_groups() != model.diffusion.delayed_groups()
        || before.precursor_counts().len() != nodes
        || after.precursor_counts().len() != nodes
        || step.delayed_decays.len() != nodes
        || step.reaction_counts.len() != model.diffusion.reactions().len()
        || step.injected_precursor_counts.len() != nodes
        || step.precursor_transfer_counts.len() != transport.transfers.len() * families
        || step.precursor_outflow_counts.len() != transport.outflows.len() * families
    {
        return Err(NeutronKineticError::Invalid(
            "unchanged original delayed families and complete native histories required",
        ));
    }
    let mut initial = Accumulator::new(0.0);
    let mut final_energy = Accumulator::new(0.0);
    let mut injected = Accumulator::new(0.0);
    let mut converted = Accumulator::new(0.0);
    let mut induced = Accumulator::new(0.0);
    let mut induced_counts = Accumulator::new(0.0);
    let mut cell_induced = vec![Accumulator::new(0.0); cells];
    let mut balance = vec![Accumulator::new(0.0); nodes];
    for cell in 0..cells {
        for (family, &energy) in model.precursor_energies_j.iter().enumerate() {
            checkpoint()?;
            let node = cell * families + family;
            let start = product(before.precursor_counts()[node], energy, defect)?;
            let end = product(after.precursor_counts()[node], energy, defect)?;
            let input = product(step.injected_precursor_counts[node], energy, defect)?;
            let decays = step.delayed_decays[node];
            let conversion = product(decays, energy, defect)?;
            add(&mut initial, start, defect)?;
            add(&mut final_energy, end, defect)?;
            add(&mut injected, input, defect)?;
            add(&mut converted, conversion, defect)?;
            for value in [end, -start, -input, conversion] {
                add(&mut balance[node], value, defect)?;
            }
            // Native production rounds count * spectrum for each group. Retain
            // those actual products and charge their difference from the reserve;
            // changing or normalizing a family spectrum would change its owner.
            let mut emitted = Accumulator::new(0.0);
            for (&weight, &neutron_energy) in model.diffusion.delayed_groups()[family]
                .spectrum
                .iter()
                .zip(&model.energies_j)
            {
                let born = checked(decays * weight)?;
                if decays > 0.0 && weight > 0.0 && born == 0.0 {
                    return Err(NeutronKineticError::Arithmetic(
                        "delayed neutron birth underflow",
                    ));
                }
                defect.add(checked(decays.mul_add(weight, -born) * neutron_energy)?.abs())?;
                let energy = product(born, neutron_energy, defect)?;
                add(&mut emitted, energy, defect)?;
            }
            defect.add(checked(sum(&[emitted.value(), -conversion]))?.abs())?;
            add(&mut neutron_balance[cell], -emitted.value(), defect)?;
        }
    }
    // Native induced births are not injection. Reconstruct the exact native
    // count * yield products, retaining their family identity and rounding.
    for &index in &model.induced_fission_reactions {
        checkpoint()?;
        let reaction = &model.diffusion.reactions()[index];
        let DiffusionReactionKind::Fission { delayed_yields, .. } = &reaction.kind else {
            continue;
        };
        let count = step.reaction_counts[index];
        for (family, (&yield_, &energy)) in delayed_yields
            .iter()
            .zip(&model.precursor_energies_j)
            .enumerate()
        {
            checkpoint()?;
            let born = checked(count * yield_)?;
            if count > 0.0 && yield_ > 0.0 && born == 0.0 {
                return Err(NeutronKineticError::Arithmetic(
                    "induced-fission precursor birth underflow",
                ));
            }
            defect.add(checked(count.mul_add(yield_, -born) * energy)?.abs())?;
            let reserve = product(born, energy, defect)?;
            let count_roundoff = induced_counts.add(born)?.abs();
            defect.add(checked(count_roundoff * energy)?.abs())?;
            add(
                &mut balance[reaction.cell * families + family],
                -reserve,
                defect,
            )?;
            add(&mut cell_induced[reaction.cell], reserve, defect)?;
            add(&mut induced, reserve, defect)?;
        }
    }
    let count_residual = checked(sum(&[
        induced_counts.value(),
        -step.delayed_precursors_produced,
    ]))?;
    let largest_reserve = model
        .precursor_energies_j
        .iter()
        .copied()
        .fold(0.0_f64, f64::max);
    defect.add(checked(count_residual.abs() * largest_reserve)?)?;
    for (face, transfer) in transport.transfers.iter().enumerate() {
        checkpoint()?;
        for (family, &energy) in model.precursor_energies_j.iter().enumerate() {
            let moved = product(
                step.precursor_transfer_counts[face * families + family],
                energy,
                defect,
            )?;
            add(
                &mut balance[transfer.from_cell * families + family],
                moved,
                defect,
            )?;
            add(
                &mut balance[transfer.to_cell * families + family],
                -moved,
                defect,
            )?;
        }
    }
    let mut exported = Accumulator::new(0.0);
    let mut cell_exported = vec![Accumulator::new(0.0); cells];
    for (face, outflow) in transport.outflows.iter().enumerate() {
        checkpoint()?;
        for (family, &energy) in model.precursor_energies_j.iter().enumerate() {
            let lost = product(
                step.precursor_outflow_counts[face * families + family],
                energy,
                defect,
            )?;
            add(&mut balance[outflow.cell * families + family], lost, defect)?;
            add(&mut cell_exported[outflow.cell], lost, defect)?;
            add(&mut exported, lost, defect)?;
        }
    }
    for entry in balance {
        checkpoint()?;
        defect.add(entry.value().abs())?;
    }
    Ok(Accounts {
        initial_j: initial.value(),
        final_j: final_energy.value(),
        injected_j: injected.value(),
        converted_j: converted.value(),
        cell_induced_j: cell_induced.into_iter().map(Accumulator::value).collect(),
        induced_j: induced.value(),
        cell_exported_j: cell_exported.into_iter().map(Accumulator::value).collect(),
        exported_j: exported.value(),
    })
}

fn add(
    account: &mut Accumulator,
    value: f64,
    defect: &mut Accumulator,
) -> Result<(), NeutronKineticError> {
    defect.add(account.add(value)?.abs())?;
    Ok(())
}

fn product(count: f64, energy: f64, defect: &mut Accumulator) -> Result<f64, NeutronKineticError> {
    let value = checked(count * energy)?;
    if count > 0.0 && energy > 0.0 && value == 0.0 {
        return Err(NeutronKineticError::Arithmetic(
            "precursor energy underflow",
        ));
    }
    defect.add(count.mul_add(energy, -value).abs())?;
    Ok(value)
}
