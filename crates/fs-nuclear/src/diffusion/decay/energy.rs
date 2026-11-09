//! Radioactive neutron kinetic energy as an internal, disjoint energy transfer.
//!
//! Each parent's emitted group energy is debited from its original outward
//! radiation account. The original neutron energy owner attaches absorption,
//! scattering, prompt Q, precursor storage and leakage to the SAME coupled
//! diffusion history. Separately funded induced-fission families retain their
//! future-neutron conversion and outflow accounts on that original history.
//! Only actual local deposition is heat; neutrons still in flight retain their
//! kinetic energy, and leakage remains an outgoing export. Neither emitted
//! neutron energy nor its later absorbed energy is counted twice.
//!
//! These are expected parent-wide spectra. A branch's sampled energy or daughter
//! correlation is not inferred. The caller supplies compatible group energies,
//! transport coefficients and outward routing; the evaluated spectral adapter
//! can bind those inputs without an interval-average source approximation.

use std::sync::Arc;

use super::{DecayNeutronBudget, DecayNeutronError, DecayNeutronState, DecayNeutronStep};
use crate::fusion_neutrons::{
    NeutronKineticError, NeutronKineticModel, NeutronKineticState, NeutronKineticStep,
};
use crate::{Accumulator, EnergyChannels, FissionError, channel_values, nonnegative, sum};

/// One joint count budget and one noncancelling combined energy allowance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecayNeutronEnergyBudget {
    pub spatial: DecayNeutronBudget,
    pub energy_tolerance_j: f64,
}

#[derive(Clone, Copy, Debug)]
struct Partition {
    neutron_j: f64,
    radiation_j: f64,
    roundoff_j: f64,
}

/// Complete accepted radioactive, neutron and disjoint energy histories.
/// Earlier source emission counters remain owned by the supplied inventories;
/// this adapter's combined cell exports start at zero at binding.
#[derive(Clone, Debug)]
pub struct DecayNeutronEnergyState {
    spatial: DecayNeutronState,
    kinetic: NeutronKineticState,
    partitions: Arc<Vec<Vec<Partition>>>,
    maximum_source_energy_j: f64,
    maximum_neutron_energy_j: f64,
    released: Vec<[Accumulator; 3]>,
}

#[derive(Clone, Debug)]
pub struct DecayNeutronEnergyStep {
    pub spatial: DecayNeutronStep,
    /// Energy of the same neutron history; no second propagation occurs.
    pub kinetic: NeutronKineticStep,
    /// Original raw source exports, including energy sent into the neutron field.
    /// Do not add this diagnostic to the disjoint combined exports below.
    pub decay_released: EnergyChannels,
    /// Outward decay radiation internally transferred into neutron storage [J].
    pub transferred_neutron_energy_j: f64,
    /// Disjoint local deposition, remaining outward radiation and escape by cell.
    pub cell_released: Vec<EnergyChannels>,
    pub released: EnergyChannels,
    /// Final neutron and precursor storage - initial + disjoint exports
    /// - radioactive release - prompt Q - additional induced conversion [J].
    pub energy_balance_residual_j: f64,
    pub energy_roundoff_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DecayNeutronEnergyError {
    Invalid(&'static str),
    Budget(&'static str),
    Source(DecayNeutronError),
    Kinetic(NeutronKineticError),
    Fission(FissionError),
}
impl From<DecayNeutronError> for DecayNeutronEnergyError {
    fn from(error: DecayNeutronError) -> Self {
        Self::Source(error)
    }
}
impl From<NeutronKineticError> for DecayNeutronEnergyError {
    fn from(error: NeutronKineticError) -> Self {
        Self::Kinetic(error)
    }
}
impl From<FissionError> for DecayNeutronEnergyError {
    fn from(error: FissionError) -> Self {
        Self::Fission(error)
    }
}
impl core::fmt::Display for DecayNeutronEnergyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "radioactive neutron energy refused: {self:?}")
    }
}
impl std::error::Error for DecayNeutronEnergyError {}

impl DecayNeutronEnergyState {
    /// Bind the exact original neutron energy model. Each parent's assigned
    /// kinetic energy must fit within its existing outward radiation account.
    /// Locally deposited or escaped energy cannot also become a neutron source.
    pub fn new(
        model: Arc<NeutronKineticModel>,
        spatial: DecayNeutronState,
    ) -> Result<Self, DecayNeutronEnergyError> {
        Self::prepare(model, spatial, false)
    }

    /// Admit separately funded induced families in the prescribed background.
    ///
    /// The caller declares that these background fission products are outside
    /// every represented radioactive source cohort. The group-only background
    /// has no fragment isotope mapping from which this can be inferred. Use
    /// finite activation with explicit product bindings for represented targets.
    /// This declaration does not reset original source or precursor populations.
    pub fn new_with_unrepresented_fission_precursors(
        model: Arc<NeutronKineticModel>,
        spatial: DecayNeutronState,
    ) -> Result<Self, DecayNeutronEnergyError> {
        Self::prepare(model, spatial, true)
    }

    fn prepare(
        model: Arc<NeutronKineticModel>,
        spatial: DecayNeutronState,
        unrepresented_fission_products: bool,
    ) -> Result<Self, DecayNeutronEnergyError> {
        if !unrepresented_fission_products
            && model
                .fission_precursor_energy_j()
                .iter()
                .any(|&energy| energy > 0.0)
        {
            return Err(DecayNeutronEnergyError::Invalid(
                "funded induced families require an explicit unrepresented background-fragment cohort",
            ));
        }
        if !Arc::ptr_eq(model.diffusion(), spatial.neutrons().model()) {
            return Err(DecayNeutronEnergyError::Invalid(
                "radioactive source requires its original neutron energy model",
            ));
        }
        let kinetic = NeutronKineticState::new_with_precursors(
            Arc::clone(&model),
            spatial.neutrons().neutron_counts().to_vec(),
            spatial.neutrons().precursor_counts().to_vec(),
        )?;
        let mut partitions = Vec::with_capacity(spatial.sources().len());
        let mut maximum_source_energy_j = 0.0_f64;
        for source in spatial.sources() {
            let nodes = source.inventory.model().nuclides().len();
            let mut outward = vec![Accumulator::new(0.0); nodes];
            let mut defects = vec![Accumulator::new(0.0); nodes];
            for branch in source.inventory.model().branches() {
                let energy = multiply(
                    branch.probability,
                    branch.energy_per_decay.radiation_j,
                    &mut defects[branch.parent],
                )?;
                add(
                    &mut outward[branch.parent],
                    energy,
                    &mut defects[branch.parent],
                )?;
                maximum_source_energy_j =
                    maximum_source_energy_j.max(branch.energy_per_decay.total_j());
            }
            let mut source_partitions = Vec::with_capacity(nodes);
            for (parent, row) in source.emission.neutron_yields().iter().enumerate() {
                let mut neutrons = Accumulator::new(0.0);
                for (&yield_, &energy) in row.iter().zip(model.group_energies_j()) {
                    let value = multiply(yield_, energy, &mut defects[parent])?;
                    add(&mut neutrons, value, &mut defects[parent])?;
                }
                let neutron_j = finite(neutrons.value())?;
                let radiation_j = finite(sum(&[
                    outward[parent].hi,
                    outward[parent].lo,
                    -neutrons.hi,
                    -neutrons.lo,
                ]))?;
                if radiation_j < 0.0 {
                    return Err(DecayNeutronEnergyError::Invalid(
                        "neutron kinetic emission exceeds the parent's outward radiation",
                    ));
                }
                defects[parent].add(finite(sum(&[neutrons.hi, neutrons.lo, -neutron_j]))?.abs())?;
                defects[parent].add(
                    finite(sum(&[
                        radiation_j,
                        neutron_j,
                        -outward[parent].hi,
                        -outward[parent].lo,
                    ]))?
                    .abs(),
                )?;
                source_partitions.push(Partition {
                    neutron_j,
                    radiation_j,
                    roundoff_j: bound(defects[parent])?,
                });
            }
            partitions.push(source_partitions);
        }
        let maximum_neutron_energy_j = finite(sum(&[
            model
                .group_energies_j()
                .iter()
                .chain(model.precursor_energies_j())
                .copied()
                .fold(0.0_f64, f64::max),
            model
                .reaction_energies()
                .iter()
                .map(|energy| energy.nuclear_q_j.abs())
                .fold(0.0_f64, f64::max),
            model
                .fission_precursor_energy_j()
                .iter()
                .copied()
                .fold(0.0_f64, f64::max),
        ]))?;
        let released = vec![[Accumulator::new(0.0); 3]; model.diffusion().cell_count()];
        Ok(Self {
            spatial,
            kinetic,
            partitions: Arc::new(partitions),
            maximum_source_energy_j,
            maximum_neutron_energy_j,
            released,
        })
    }

    pub fn model(&self) -> &Arc<NeutronKineticModel> {
        self.kinetic.model()
    }
    pub fn spatial(&self) -> &DecayNeutronState {
        &self.spatial
    }
    pub fn kinetic(&self) -> &NeutronKineticState {
        &self.kinetic
    }
    pub fn cell_released(&self, cell: usize) -> Option<EnergyChannels> {
        self.released.get(cell).copied().map(channel_values)
    }
    pub fn cumulative_cell_released(&self) -> Vec<EnergyChannels> {
        self.released.iter().copied().map(channel_values).collect()
    }

    pub fn trial(
        &self,
        duration_s: f64,
        budget: DecayNeutronEnergyBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, DecayNeutronEnergyStep), DecayNeutronEnergyError> {
        checkpoint()?;
        if !nonnegative(budget.energy_tolerance_j) {
            return Err(DecayNeutronEnergyError::Invalid(
                "finite nonnegative combined energy allowance required",
            ));
        }
        let (spatial, step) = self
            .spatial
            .trial(duration_s, budget.spatial, &mut checkpoint)?;
        let (kinetic, kinetic_step) = self.kinetic.account_accepted_diffusion(
            spatial.neutrons().clone(),
            step.neutrons.clone(),
            budget.energy_tolerance_j,
            &mut checkpoint,
        )?;
        let mut error = Accumulator::new(0.0);
        for defect in [
            kinetic_step.energy_roundoff_j,
            step.decay_energy_roundoff_j,
            finite(step.decay_population_roundoff * self.maximum_source_energy_j)?,
            finite(step.neutrons.population_roundoff * self.maximum_neutron_energy_j)?,
        ] {
            error.add(defect)?;
        }
        let mut cells = vec![[Accumulator::new(0.0); 3]; self.released.len()];
        for (cell, energy) in cells.iter_mut().zip(&kinetic_step.cell_released) {
            for (total, value) in cell.iter_mut().zip(values(*energy)) {
                add(total, value, &mut error)?;
            }
        }
        let mut raw = [Accumulator::new(0.0); 3];
        let mut transferred = Accumulator::new(0.0);
        for ((source, report), partitions) in self
            .spatial
            .sources()
            .iter()
            .zip(&step.decay)
            .zip(self.partitions.iter())
        {
            checkpoint()?;
            let mut remaining = Accumulator::new(0.0);
            let mut emitted = Accumulator::new(0.0);
            for (&count, partition) in report.decays.iter().zip(partitions) {
                error.add(finite(count * partition.roundoff_j)?)?;
                let radiated = multiply(count, partition.radiation_j, &mut error)?;
                let born = multiply(count, partition.neutron_j, &mut error)?;
                add(&mut remaining, radiated, &mut error)?;
                add(&mut emitted, born, &mut error)?;
            }
            error.add(
                finite(sum(&[
                    remaining.hi,
                    remaining.lo,
                    emitted.hi,
                    emitted.lo,
                    -report.released.radiation_j,
                ]))?
                .abs(),
            )?;
            add(&mut transferred, emitted.value(), &mut error)?;
            for (total, value) in raw.iter_mut().zip(values(report.released)) {
                add(total, value, &mut error)?;
            }
            for (total, value) in cells[source.cell].iter_mut().zip([
                report.released.deposited_j,
                remaining.value(),
                report.released.escaped_j,
            ]) {
                add(total, value, &mut error)?;
            }
            // Charge retained low parts lost on the independent source exports.
            for account in [remaining, emitted] {
                error.add(finite(sum(&[account.hi, account.lo, -account.value()]))?.abs())?;
            }
        }
        let transferred_neutron_energy_j = finite(transferred.value())?;
        error.add(
            finite(sum(&[
                transferred.hi,
                transferred.lo,
                -transferred_neutron_energy_j,
            ]))?
            .abs(),
        )?;
        error.add(
            finite(sum(&[
                transferred_neutron_energy_j,
                -kinetic_step.injected_kinetic_energy_j,
            ]))?
            .abs(),
        )?;
        let decay_released = export(raw, &mut error)?;
        let mut released = self.released.clone();
        let mut total = [Accumulator::new(0.0); 3];
        let mut cell_released = Vec::with_capacity(cells.len());
        for (cell, accumulated) in cells.into_iter().zip(&mut released) {
            checkpoint()?;
            let channels = export(cell, &mut error)?;
            for ((cumulative, total), value) in
                accumulated.iter_mut().zip(&mut total).zip(values(channels))
            {
                add(cumulative, value, &mut error)?;
                add(total, value, &mut error)?;
            }
            finite(channel_values(*accumulated).total_j())?;
            cell_released.push(channels);
        }
        let combined = export(total, &mut error)?;
        let energy_balance_residual_j = finite(sum(&[
            kinetic_step.final_kinetic_energy_j,
            -kinetic_step.initial_kinetic_energy_j,
            kinetic_step.final_precursor_energy_j,
            -kinetic_step.initial_precursor_energy_j,
            kinetic_step.exported_precursor_energy_j,
            combined.deposited_j,
            combined.radiation_j,
            combined.escaped_j,
            -decay_released.deposited_j,
            -decay_released.radiation_j,
            -decay_released.escaped_j,
            -kinetic_step.nuclear_reaction_energy_j,
            -kinetic_step.induced_precursor_energy_j,
        ]))?;
        error.add(energy_balance_residual_j.abs())?;
        let energy_roundoff_j = bound(error)?;
        if energy_roundoff_j > budget.energy_tolerance_j {
            return Err(DecayNeutronEnergyError::Budget(
                "combined radioactive/neutron energy accounts",
            ));
        }
        checkpoint()?;
        Ok((
            Self {
                spatial,
                kinetic,
                partitions: Arc::clone(&self.partitions),
                maximum_source_energy_j: self.maximum_source_energy_j,
                maximum_neutron_energy_j: self.maximum_neutron_energy_j,
                released,
            },
            DecayNeutronEnergyStep {
                spatial: step,
                kinetic: kinetic_step,
                decay_released,
                transferred_neutron_energy_j,
                cell_released,
                released: combined,
                energy_balance_residual_j,
                energy_roundoff_j,
            },
        ))
    }

    pub fn advance(
        &mut self,
        duration_s: f64,
        budget: DecayNeutronEnergyBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayNeutronEnergyStep, DecayNeutronEnergyError> {
        let (next, step) = self.trial(duration_s, budget, checkpoint)?;
        *self = next;
        Ok(step)
    }
}

fn values(energy: EnergyChannels) -> [f64; 3] {
    [energy.deposited_j, energy.radiation_j, energy.escaped_j]
}
fn finite(value: f64) -> Result<f64, DecayNeutronEnergyError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DecayNeutronEnergyError::Invalid(
            "unrepresentable radioactive neutron energy",
        ))
    }
}
fn multiply(a: f64, b: f64, error: &mut Accumulator) -> Result<f64, DecayNeutronEnergyError> {
    let result = finite(a * b)?;
    if a > 0.0 && b > 0.0 && result == 0.0 {
        return Err(DecayNeutronEnergyError::Invalid(
            "positive radioactive neutron energy underflow",
        ));
    }
    error.add(a.mul_add(b, -result).abs())?;
    Ok(result)
}
fn add(
    total: &mut Accumulator,
    value: f64,
    error: &mut Accumulator,
) -> Result<(), DecayNeutronEnergyError> {
    error.add(total.add(value)?.abs())?;
    Ok(())
}
fn export(
    accounts: [Accumulator; 3],
    error: &mut Accumulator,
) -> Result<EnergyChannels, DecayNeutronEnergyError> {
    let energy = channel_values(accounts);
    finite(energy.total_j())?;
    for (account, value) in accounts.iter().zip(values(energy)) {
        error.add(finite(sum(&[account.hi, account.lo, -value]))?.abs())?;
    }
    Ok(energy)
}
fn bound(account: Accumulator) -> Result<f64, DecayNeutronEnergyError> {
    let value = finite(account.value())?;
    let remainder = finite(sum(&[account.hi, account.lo, -value]))?;
    finite(if remainder > 0.0 {
        value.next_up()
    } else {
        value
    })
}

#[cfg(test)]
mod tests;

#[cfg(test)]
pub(crate) use tests::material_induced_decay_fixture;
