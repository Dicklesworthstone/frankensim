//! Internal neutron-energy transfer from the original finite radioactive stock.
//! Parent-wide spectra use the admitted decay-branch rates, not capture energy
//! or renormalized nuclear data. This module accounts an accepted trajectory;
//! it never prepares or advances a second neutron or inventory propagator.

use super::*;

#[derive(Clone, Copy, Debug)]
struct Partition {
    neutrons_per_decay: f64,
    precursors_per_decay: f64,
    neutron_j: f64,
    precursor_j: f64,
    radiation_j: f64,
    roundoff_j: f64,
}

#[derive(Clone, Debug)]
pub(super) struct DecayPartitions(Vec<Vec<Option<Partition>>>);

pub(super) struct DecayTransfer {
    pub remaining_radiation_j: Vec<f64>,
    pub decay_neutrons: f64,
    pub decay_precursors: f64,
    pub transferred_energy_j: f64,
    pub reserved_precursor_energy_j: f64,
    pub external_neutrons: f64,
    pub external_precursors: f64,
    pub external_energy_j: f64,
    pub external_precursor_energy_j: f64,
    pub number_roundoff: f64,
    pub energy_roundoff_j: f64,
}

impl DecayPartitions {
    pub fn new(
        model: &NeutronKineticModel,
        spatial: &DiffusionBurnupState,
    ) -> Result<Option<Arc<Self>>, FusionNeutronActivationError> {
        let Some(emissions) = spatial.decay_neutrons() else {
            return Ok(None);
        };
        let mut cells = Vec::with_capacity(emissions.len());
        for emission in emissions {
            let mut parents = Vec::with_capacity(emission.neutron_yields().len());
            for (parent, yields) in emission.neutron_yields().iter().enumerate() {
                let mut error = Accumulator::new(0.0);
                let (neutrons_per_decay, neutron_j) =
                    moments(yields, model.group_energies_j(), &mut error)?;
                let (precursors_per_decay, precursor_j) = moments(
                    &emission.precursor_yields()[parent],
                    model.precursor_energies_j(),
                    &mut error,
                )?;
                let reserved = checked(sum(&[neutron_j, precursor_j]))?;
                error.add(sum(&[neutron_j, precursor_j, -reserved]).abs())?;
                if reserved == 0.0 {
                    parents.push(None);
                    continue;
                }
                let lambda = emission.model().nuclides()[parent]
                    .law
                    .decay_constant()
                    .value();
                let branches: Vec<_> = emission
                    .model()
                    .reactions()
                    .iter()
                    .filter_map(|reaction| {
                        if reaction.parent == parent {
                            if let BurnupRate::Decay { branch_fraction } = reaction.rate {
                                return Some((
                                    lambda * branch_fraction,
                                    reaction.energy_per_event.radiation_j,
                                ));
                            }
                        }
                        None
                    })
                    .collect();
                let scale = branches
                    .iter()
                    .map(|branch| branch.0)
                    .fold(0.0_f64, f64::max);
                if !scale.is_finite() || scale <= 0.0 {
                    return Err(FusionNeutronActivationError::Invalid(
                        "emitted kinetic energy requires positive original decay rates",
                    ));
                }
                // A rate-weighted conditional mean is not a change to the
                // supplied branches. Scaling the rates avoids overflow and
                // preserves equality when every branch spends its full share.
                let mut weights = Accumulator::new(0.0);
                let mut weight_error = Accumulator::new(0.0);
                let mut remaining = Accumulator::new(0.0);
                let mut partition_error = Accumulator::new(0.0);
                for (rate, radiation) in branches {
                    let weight = checked(rate / scale)?;
                    if rate > 0.0 && weight == 0.0 {
                        return Err(FissionError::Arithmetic.into());
                    }
                    let difference = checked(sum(&[radiation, -reserved]))?;
                    let division_error = checked(weight.mul_add(scale, -rate).abs() / scale)?;
                    partition_error.add(checked(division_error * difference.abs())?)?;
                    partition_error.add(checked(
                        sum(&[radiation, -reserved, -difference]).abs() * weight,
                    )?)?;
                    let value = product(weight, difference, &mut partition_error)?;
                    add(&mut remaining, value, &mut partition_error)?;
                    add(&mut weights, weight, &mut weight_error)?;
                }
                let denominator = export(weights, &mut weight_error)?;
                let numerator = export(remaining, &mut partition_error)?;
                if numerator < 0.0 {
                    return Err(FusionNeutronActivationError::Invalid(
                        "immediate and reserved delayed-neutron energy exceeds original decay radiation",
                    ));
                }
                let radiation_j = checked(numerator / denominator)?;
                if numerator > 0.0 && radiation_j == 0.0 {
                    return Err(FissionError::Arithmetic.into());
                }
                partition_error.add(radiation_j.mul_add(denominator, -numerator).abs())?;
                partition_error.add(checked(radiation_j * upper(weight_error)?)?)?;
                error.add(checked(upper(partition_error)? / denominator)?)?;
                parents.push(Some(Partition {
                    neutrons_per_decay,
                    precursors_per_decay,
                    neutron_j,
                    precursor_j,
                    radiation_j,
                    roundoff_j: upper(error)?,
                }));
            }
            cells.push(parents);
        }
        Ok(Some(Arc::new(Self(cells))))
    }

    pub fn account(
        &self,
        model: &NeutronKineticModel,
        inventories: &[BurnupInventory],
        step: &DiffusionBurnupStep,
        duration_s: f64,
        external_rates_s: &[f64],
        external_precursor_rates_s: &[f64],
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayTransfer, FusionNeutronActivationError> {
        let births =
            step.decay_neutron_counts
                .as_ref()
                .ok_or(FusionNeutronActivationError::Invalid(
                    "accepted decay emission account is missing",
                ))?;
        let groups = model.group_energies_j().len();
        let families = model.precursor_energies_j().len();
        let precursor_births = step.decay_precursor_counts.as_ref().ok_or(
            FusionNeutronActivationError::Invalid("accepted decay precursor account is missing"),
        )?;
        if births.len() != self.0.len()
            || births.iter().any(|row| row.len() != groups)
            || precursor_births.len() != self.0.len()
            || precursor_births.iter().any(|row| row.len() != families)
            || external_rates_s.len() != self.0.len() * groups
            || external_precursor_rates_s.len() != self.0.len() * families
        {
            return Err(FusionNeutronActivationError::Invalid(
                "accepted decay emission layout changed",
            ));
        }
        let mut energy_error = Accumulator::new(0.0);
        let mut number_error = Accumulator::new(0.0);
        let mut emitted = Accumulator::new(0.0);
        let mut born_precursors = Accumulator::new(0.0);
        let mut transferred = Accumulator::new(0.0);
        let mut reserved_precursors = Accumulator::new(0.0);
        let mut external = Accumulator::new(0.0);
        let mut external_precursors = Accumulator::new(0.0);
        let mut external_energy = Accumulator::new(0.0);
        let mut external_precursor_energy = Accumulator::new(0.0);
        let mut remaining_radiation_j = Vec::with_capacity(self.0.len());
        for (cell, ((partitions, inventory), report)) in
            self.0.iter().zip(inventories).zip(&step.burnup).enumerate()
        {
            checkpoint()?;
            let mut parent_decays = vec![Accumulator::new(0.0); partitions.len()];
            let mut parent_count_errors = vec![Accumulator::new(0.0); partitions.len()];
            let mut remaining = Accumulator::new(0.0);
            let mut cell_transfer = Accumulator::new(0.0);
            let mut cell_reservation = Accumulator::new(0.0);
            let mut cell_born_energy = Accumulator::new(0.0);
            let mut cell_precursor_energy = Accumulator::new(0.0);
            for (reaction, &count) in inventory
                .model()
                .reactions()
                .iter()
                .zip(&report.reaction_counts)
            {
                checkpoint()?;
                if matches!(reaction.rate, BurnupRate::Decay { .. })
                    && partitions[reaction.parent].is_some()
                {
                    add(
                        &mut parent_decays[reaction.parent],
                        count,
                        &mut parent_count_errors[reaction.parent],
                    )?;
                } else {
                    let value = product(
                        count,
                        reaction.energy_per_event.radiation_j,
                        &mut energy_error,
                    )?;
                    add(&mut remaining, value, &mut energy_error)?;
                }
            }
            for ((counts, mut count_error), partition) in parent_decays
                .into_iter()
                .zip(parent_count_errors)
                .zip(partitions)
            {
                if let Some(partition) = partition {
                    let count = export(counts, &mut count_error)?;
                    let count_error = upper(count_error)?;
                    let neutron_error =
                        product(count_error, partition.neutrons_per_decay, &mut number_error)?;
                    number_error.add(neutron_error)?;
                    let precursor_error =
                        product(count_error, partition.precursors_per_decay, &mut number_error)?;
                    number_error.add(precursor_error)?;
                    let weighted_error = product(
                        count_error,
                        checked(sum(&[
                            partition.neutron_j,
                            partition.precursor_j,
                            partition.radiation_j,
                        ]))?,
                        &mut energy_error,
                    )?;
                    energy_error.add(weighted_error)?;
                    energy_error.add(checked(count * partition.roundoff_j)?)?;
                    let radiation = product(count, partition.radiation_j, &mut energy_error)?;
                    let transfer = product(count, partition.neutron_j, &mut energy_error)?;
                    let reservation = product(count, partition.precursor_j, &mut energy_error)?;
                    add(&mut remaining, radiation, &mut energy_error)?;
                    add(&mut cell_transfer, transfer, &mut energy_error)?;
                    add(&mut cell_reservation, reservation, &mut energy_error)?;
                }
            }
            let radiation = export(remaining, &mut energy_error)?;
            let transfer = export(cell_transfer, &mut energy_error)?;
            let reservation = export(cell_reservation, &mut energy_error)?;
            energy_error
                .add(checked(sum(&[radiation, transfer, reservation, -report.released.radiation_j]))?.abs())?;
            add(&mut transferred, transfer, &mut energy_error)?;
            add(&mut reserved_precursors, reservation, &mut energy_error)?;
            remaining_radiation_j.push(radiation);
            for (group, &internal) in births[cell].iter().enumerate() {
                checkpoint()?;
                let index = cell * groups + group;
                let supplied = product(external_rates_s[index], duration_s, &mut number_error)?;
                number_error.add(
                    checked(sum(&[
                        step.neutrons.injected_neutron_counts[index],
                        -internal,
                        -supplied,
                    ]))?
                    .abs(),
                )?;
                add(&mut emitted, internal, &mut number_error)?;
                add(&mut external, supplied, &mut number_error)?;
                let energy = model.group_energies_j()[group];
                let internal_energy = product(internal, energy, &mut energy_error)?;
                let supplied_energy = product(supplied, energy, &mut energy_error)?;
                add(&mut cell_born_energy, internal_energy, &mut energy_error)?;
                add(&mut external_energy, supplied_energy, &mut energy_error)?;
            }
            let born_energy_j = export(cell_born_energy, &mut energy_error)?;
            energy_error.add(checked(sum(&[transfer, -born_energy_j]))?.abs())?;
            for (family, &internal) in precursor_births[cell].iter().enumerate() {
                checkpoint()?;
                let index = cell * families + family;
                let supplied = product(external_precursor_rates_s[index], duration_s, &mut number_error)?;
                number_error.add(checked(sum(&[
                    step.neutrons.injected_precursor_counts[index],
                    -internal,
                    -supplied,
                ]))?.abs())?;
                add(&mut born_precursors, internal, &mut number_error)?;
                add(&mut external_precursors, supplied, &mut number_error)?;
                let energy = model.precursor_energies_j()[family];
                let internal_energy = product(internal, energy, &mut energy_error)?;
                let supplied_energy = product(supplied, energy, &mut energy_error)?;
                add(&mut cell_precursor_energy, internal_energy, &mut energy_error)?;
                add(&mut external_precursor_energy, supplied_energy, &mut energy_error)?;
            }
            let born_precursor_energy_j = export(cell_precursor_energy, &mut energy_error)?;
            energy_error.add(checked(sum(&[reservation, -born_precursor_energy_j]))?.abs())?;
        }
        let decay_neutrons = export(emitted, &mut number_error)?;
        let decay_precursors = export(born_precursors, &mut number_error)?;
        let external_neutrons = export(external, &mut number_error)?;
        let external_precursors = export(external_precursors, &mut number_error)?;
        let transferred_energy_j = export(transferred, &mut energy_error)?;
        let reserved_precursor_energy_j = export(reserved_precursors, &mut energy_error)?;
        let external_energy_j = export(external_energy, &mut energy_error)?;
        let external_precursor_energy_j = export(external_precursor_energy, &mut energy_error)?;
        let number_roundoff = upper(number_error)?;
        energy_error.add(checked(
            number_roundoff
                * model
                    .group_energies_j()
                    .iter()
                    .chain(model.precursor_energies_j())
                    .copied()
                    .fold(0.0_f64, f64::max),
        )?)?;
        Ok(DecayTransfer {
            remaining_radiation_j,
            decay_neutrons,
            decay_precursors,
            transferred_energy_j,
            reserved_precursor_energy_j,
            external_neutrons,
            external_precursors,
            external_energy_j,
            external_precursor_energy_j,
            number_roundoff,
            energy_roundoff_j: upper(energy_error)?,
        })
    }
}

fn moments(
    yields: &[f64],
    energies: &[f64],
    error: &mut Accumulator,
) -> Result<(f64, f64), FusionNeutronActivationError> {
    let mut energy = Accumulator::new(0.0);
    let mut multiplicity = Accumulator::new(0.0);
    let mut count_error = Accumulator::new(0.0);
    for (&yield_, &coefficient) in yields.iter().zip(energies) {
        add(&mut multiplicity, yield_, &mut count_error)?;
        let value = product(yield_, coefficient, error)?;
        add(&mut energy, value, error)?;
    }
    multiplicity.add(count_error.hi)?;
    multiplicity.add(count_error.lo)?;
    Ok((upper(multiplicity)?, export(energy, error)?))
}

fn product(a: f64, b: f64, error: &mut Accumulator) -> Result<f64, FusionNeutronActivationError> {
    let value = checked(a * b)?;
    if a != 0.0 && b != 0.0 && value == 0.0 {
        return Err(FissionError::Arithmetic.into());
    }
    error.add(a.mul_add(b, -value).abs())?;
    Ok(value)
}
fn add(
    account: &mut Accumulator,
    value: f64,
    error: &mut Accumulator,
) -> Result<(), FusionNeutronActivationError> {
    error.add(account.add(value)?.abs())?;
    Ok(())
}
fn export(
    account: Accumulator,
    error: &mut Accumulator,
) -> Result<f64, FusionNeutronActivationError> {
    let value = checked(account.value())?;
    error.add(checked(sum(&[account.hi, account.lo, -value]))?.abs())?;
    Ok(value)
}
fn upper(account: Accumulator) -> Result<f64, FusionNeutronActivationError> {
    let value = checked(account.value())?;
    Ok(if sum(&[account.hi, account.lo, -value]) > 0.0 {
        checked(value.next_up())?
    } else {
        value
    })
}
