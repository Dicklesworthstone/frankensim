//! Source-owned coordinates for the shared time-refinement comparison.
//! No populations are re-evolved here and no lifetime counter is subtracted to
//! reconstruct a small interval release. The original three interval reports
//! supply reaction/emission comparisons even after a long previous history.

use super::*;
use crate::burnup::{BurnupInventory, BurnupStep};
use crate::decay::thermal::DecayCellSource;
use crate::irradiation::IrradiationState;

mod diffusion;

/// Physical state and interval comparisons used by the shared adaptive owner.
/// Sealed through `NuclearHeatSource`; every implementation retains the original
/// nuclide/group order and compares immutable-model-bound candidate states.
pub trait AdaptiveNuclearHeatSource: NuclearHeatSource {
    /// Compare this accepted start state with coarse and fine endpoints.
    fn time_state_error(
        &self,
        coarse: &Self,
        fine: &Self,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError>;

    /// Compare one coarse interval with the sum of its two physical halves.
    /// Reports are current-interval values, never lifetime counter differences.
    fn time_step_error(
        coarse: &Self::Step,
        first: &Self::Step,
        second: &Self::Step,
        _duration_s: f64,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        let mut metric = Comparison::new(tolerances, checkpoint);
        common_steps::<Self>(&mut metric, coarse, first, second)?;
        Ok(metric.error)
    }
}

pub(super) struct Comparison<'a> {
    tolerances: NuclearTimeTolerances,
    pub(super) error: NuclearTimeError,
    checkpoint: &'a mut dyn FnMut() -> Result<(), DepositionError>,
    samples: usize,
}

impl<'a> Comparison<'a> {
    pub(super) fn new(
        tolerances: NuclearTimeTolerances,
        checkpoint: &'a mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Self {
        Self {
            tolerances,
            error: NuclearTimeError::default(),
            checkpoint,
            samples: 0,
        }
    }

    fn coordinate(
        &mut self,
        initial: f64,
        coarse: f64,
        fine: f64,
        absolute: f64,
    ) -> Result<f64, FissionThermalError> {
        if self.samples % 256 == 0 {
            (self.checkpoint)()?;
        }
        self.samples = self
            .samples
            .checked_add(1)
            .ok_or(FissionThermalError::Arithmetic)?;
        let difference = finite(coarse - fine)?.abs();
        let scale = finite(coarse - initial)?
            .abs()
            .max(finite(fine - initial)?.abs());
        let allowance = absolute + self.tolerances.relative * scale;
        if !absolute.is_finite()
            || absolute <= 0.0
            || !allowance.is_finite()
            || !self.tolerances.relative.is_finite()
            || !(0.0..1.0).contains(&self.tolerances.relative)
        {
            return Err(FissionThermalError::Invalid(
                "invalid temporal coordinate allowance",
            ));
        }
        // Overflow of a ratio for a very strict scale is a rejected estimate,
        // not permission to discard a finite physical difference.
        self.error.normalized = self.error.normalized.max(difference / allowance);
        Ok(difference)
    }
    pub(super) fn specific_energy(
        &mut self,
        start: f64,
        coarse: f64,
        fine: f64,
    ) -> Result<(), FissionThermalError> {
        let difference =
            self.coordinate(start, coarse, fine, self.tolerances.specific_energy_j_kg)?;
        self.error.specific_energy_difference_j_kg =
            self.error.specific_energy_difference_j_kg.max(difference);
        Ok(())
    }
    pub(super) fn temperature(
        &mut self,
        start: f64,
        coarse: f64,
        fine: f64,
    ) -> Result<(), FissionThermalError> {
        let difference = self.coordinate(start, coarse, fine, self.tolerances.temperature_k)?;
        self.error.temperature_difference_k = self.error.temperature_difference_k.max(difference);
        Ok(())
    }
    fn nuclei(&mut self, start: f64, coarse: f64, fine: f64) -> Result<(), FissionThermalError> {
        let difference = self.coordinate(start, coarse, fine, self.tolerances.nuclei)?;
        self.error.nuclei_difference = self.error.nuclei_difference.max(difference);
        Ok(())
    }
    fn energy(&mut self, start: f64, coarse: f64, fine: f64) -> Result<(), FissionThermalError> {
        let difference = self.coordinate(start, coarse, fine, self.tolerances.nuclear_energy_j)?;
        self.error.nuclear_energy_difference_j =
            self.error.nuclear_energy_difference_j.max(difference);
        Ok(())
    }
    fn neutron(&mut self, start: f64, coarse: f64, fine: f64) -> Result<(), FissionThermalError> {
        let difference =
            self.coordinate(start, coarse, fine, self.tolerances.neutron_population)?;
        self.error.neutron_population_difference =
            self.error.neutron_population_difference.max(difference);
        Ok(())
    }
    fn neutron_count(
        &mut self,
        start: f64,
        coarse: f64,
        fine: f64,
    ) -> Result<(), FissionThermalError> {
        let difference = self.coordinate(start, coarse, fine, self.tolerances.neutron_count)?;
        self.error.neutron_count_difference = self.error.neutron_count_difference.max(difference);
        Ok(())
    }
    fn populations(
        &mut self,
        start: &[f64],
        coarse: &[f64],
        fine: &[f64],
    ) -> Result<(), FissionThermalError> {
        if start.len() != coarse.len() || start.len() != fine.len() {
            return Err(FissionThermalError::Invalid(
                "refinement changed population ordering",
            ));
        }
        for ((&a, &b), &c) in start.iter().zip(coarse).zip(fine) {
            self.nuclei(a, b, c)?;
        }
        Ok(())
    }
    fn reservoirs(
        &mut self,
        start: &[f64],
        coarse: &[f64],
        fine: &[f64],
    ) -> Result<(), FissionThermalError> {
        if start.len() != coarse.len() || start.len() != fine.len() {
            return Err(FissionThermalError::Invalid(
                "refinement changed delayed-energy ordering",
            ));
        }
        for ((&a, &b), &c) in start.iter().zip(coarse).zip(fine) {
            self.energy(a, b, c)?;
        }
        Ok(())
    }
    fn reaction_counts(
        &mut self,
        coarse: &[f64],
        first: &[f64],
        second: &[f64],
    ) -> Result<(), FissionThermalError> {
        if coarse.len() != first.len() || coarse.len() != second.len() {
            return Err(FissionThermalError::Invalid(
                "refinement changed reaction ordering",
            ));
        }
        for ((&a, &b), &c) in coarse.iter().zip(first).zip(second) {
            self.nuclei(0.0, a, finite(b + c)?)?;
        }
        Ok(())
    }
    fn particle_counts(
        &mut self,
        coarse: &[f64],
        first: &[f64],
        second: &[f64],
    ) -> Result<(), FissionThermalError> {
        if coarse.len() != first.len() || coarse.len() != second.len() {
            return Err(FissionThermalError::Invalid(
                "refinement changed spatial transfer ordering",
            ));
        }
        for ((&a, &b), &c) in coarse.iter().zip(first).zip(second) {
            self.neutron_count(0.0, a, finite(b + c)?)?;
        }
        Ok(())
    }
}

fn same_binding<S: NuclearHeatSource>(
    initial: &S,
    coarse: &S,
    fine: &S,
) -> Result<(), FissionThermalError> {
    if initial.cell_id() != coarse.cell_id()
        || initial.cell_id() != fine.cell_id()
        || initial.material_id() != coarse.material_id()
        || initial.material_id() != fine.material_id()
    {
        return Err(FissionThermalError::Invalid(
            "refinement changed physical source binding",
        ));
    }
    Ok(())
}

fn same_model<T: ?Sized>(initial: &T, coarse: &T, fine: &T) -> Result<(), FissionThermalError> {
    if !std::ptr::eq(initial, coarse) || !std::ptr::eq(initial, fine) {
        return Err(FissionThermalError::Invalid(
            "refinement requires the same immutable nuclear model",
        ));
    }
    Ok(())
}

fn common_steps<S: NuclearHeatSource>(
    metric: &mut Comparison<'_>,
    coarse: &S::Step,
    first: &S::Step,
    second: &S::Step,
) -> Result<(), FissionThermalError> {
    let a = channels(S::released(coarse));
    let b = channels(S::released(first));
    let c = channels(S::released(second));
    for i in 0..3 {
        metric.energy(0.0, a[i], finite(b[i] + c[i])?)?;
    }
    for counts in [S::fissions, S::captures, S::decays, S::produced_nuclei] {
        metric.nuclei(0.0, counts(coarse), finite(counts(first) + counts(second))?)?;
    }
    Ok(())
}

fn channels(value: EnergyChannels) -> [f64; 3] {
    [value.deposited_j, value.radiation_j, value.escaped_j]
}

fn irradiation(
    metric: &mut Comparison<'_>,
    initial: &IrradiationState,
    coarse: &IrradiationState,
    fine: &IrradiationState,
) -> Result<(), FissionThermalError> {
    same_model(initial.model(), coarse.model(), fine.model())?;
    metric.nuclei(
        initial.remaining_parent_nuclei(),
        coarse.remaining_parent_nuclei(),
        fine.remaining_parent_nuclei(),
    )?;
    metric.reservoirs(
        initial.delayed_energy_j(),
        coarse.delayed_energy_j(),
        fine.delayed_energy_j(),
    )
}

fn burnup(
    metric: &mut Comparison<'_>,
    initial: &BurnupInventory,
    coarse: &BurnupInventory,
    fine: &BurnupInventory,
) -> Result<(), FissionThermalError> {
    same_model(
        initial.model().as_ref(),
        coarse.model().as_ref(),
        fine.model().as_ref(),
    )?;
    metric.populations(
        initial.populations(),
        coarse.populations(),
        fine.populations(),
    )
}

impl AdaptiveNuclearHeatSource for FissionCellSource {
    fn time_state_error(
        &self,
        coarse: &Self,
        fine: &Self,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        same_binding(self, coarse, fine)?;
        same_model(self.state.model(), coarse.state.model(), fine.state.model())?;
        let mut metric = Comparison::new(tolerances, checkpoint);
        metric.nuclei(
            self.state.remaining_fuel_nuclei(),
            coarse.state.remaining_fuel_nuclei(),
            fine.state.remaining_fuel_nuclei(),
        )?;
        metric.reservoirs(
            self.state.delayed_energy_j(),
            coarse.state.delayed_energy_j(),
            fine.state.delayed_energy_j(),
        )?;
        Ok(metric.error)
    }
}

impl AdaptiveNuclearHeatSource for IrradiationCellSource {
    fn time_state_error(
        &self,
        coarse: &Self,
        fine: &Self,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        same_binding(self, coarse, fine)?;
        if self.inventories.len() != coarse.inventories.len()
            || self.inventories.len() != fine.inventories.len()
        {
            return Err(FissionThermalError::Invalid(
                "refinement changed isotope ordering",
            ));
        }
        let mut metric = Comparison::new(tolerances, checkpoint);
        for ((initial, coarse), fine) in self
            .inventories
            .iter()
            .zip(&coarse.inventories)
            .zip(&fine.inventories)
        {
            irradiation(&mut metric, initial, coarse, fine)?;
        }
        Ok(metric.error)
    }
    fn time_step_error(
        coarse: &Self::Step,
        first: &Self::Step,
        second: &Self::Step,
        _: f64,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        if coarse.isotopes.len() != first.isotopes.len()
            || coarse.isotopes.len() != second.isotopes.len()
        {
            return Err(FissionThermalError::Invalid(
                "refinement changed isotope reports",
            ));
        }
        let mut metric = Comparison::new(tolerances, checkpoint);
        common_steps::<Self>(&mut metric, coarse, first, second)?;
        for ((a, b), c) in coarse
            .isotopes
            .iter()
            .zip(&first.isotopes)
            .zip(&second.isotopes)
        {
            metric.nuclei(0.0, a.fissions, finite(b.fissions + c.fissions)?)?;
            metric.nuclei(0.0, a.captures, finite(b.captures + c.captures)?)?;
            for ((a, b), c) in channels(a.released)
                .into_iter()
                .zip(channels(b.released))
                .zip(channels(c.released))
            {
                metric.energy(0.0, a, finite(b + c)?)?;
            }
        }
        Ok(metric.error)
    }
}

impl AdaptiveNuclearHeatSource for BurnupCellSource {
    fn time_state_error(
        &self,
        coarse: &Self,
        fine: &Self,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        same_binding(self, coarse, fine)?;
        let mut metric = Comparison::new(tolerances, checkpoint);
        burnup(
            &mut metric,
            &self.inventory,
            &coarse.inventory,
            &fine.inventory,
        )?;
        Ok(metric.error)
    }
    fn time_step_error(
        coarse: &Self::Step,
        first: &Self::Step,
        second: &Self::Step,
        _: f64,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        let mut metric = Comparison::new(tolerances, checkpoint);
        common_steps::<Self>(&mut metric, coarse, first, second)?;
        burnup_steps(&mut metric, coarse, first, second)?;
        Ok(metric.error)
    }
}

fn burnup_steps(
    metric: &mut Comparison<'_>,
    coarse: &BurnupStep,
    first: &BurnupStep,
    second: &BurnupStep,
) -> Result<(), FissionThermalError> {
    metric.reaction_counts(
        &coarse.reaction_counts,
        &first.reaction_counts,
        &second.reaction_counts,
    )
}

impl AdaptiveNuclearHeatSource for BurnupProductionCellSource {
    fn time_state_error(
        &self,
        coarse: &Self,
        fine: &Self,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        self.source()
            .time_state_error(coarse.source(), fine.source(), tolerances, checkpoint)
    }
    fn time_step_error(
        coarse: &Self::Step,
        first: &Self::Step,
        second: &Self::Step,
        _: f64,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        let mut metric = Comparison::new(tolerances, checkpoint);
        common_steps::<Self>(&mut metric, coarse, first, second)?;
        burnup_steps(&mut metric, &coarse.burnup, &first.burnup, &second.burnup)?;
        metric.reaction_counts(
            &coarse.produced_nuclei,
            &first.produced_nuclei,
            &second.produced_nuclei,
        )?;
        Ok(metric.error)
    }
}

impl AdaptiveNuclearHeatSource for DecayProductionCellSource {
    fn time_state_error(
        &self,
        coarse: &Self,
        fine: &Self,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        same_binding(self, coarse, fine)?;
        same_model(
            self.inventory().model().as_ref(),
            coarse.inventory().model().as_ref(),
            fine.inventory().model().as_ref(),
        )?;
        let mut metric = Comparison::new(tolerances, checkpoint);
        metric.populations(
            self.inventory().populations(),
            coarse.inventory().populations(),
            fine.inventory().populations(),
        )?;
        Ok(metric.error)
    }
    fn time_step_error(
        coarse: &Self::Step,
        first: &Self::Step,
        second: &Self::Step,
        _: f64,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        let mut metric = Comparison::new(tolerances, checkpoint);
        common_steps::<Self>(&mut metric, coarse, first, second)?;
        metric.reaction_counts(
            &coarse.decay.decays,
            &first.decay.decays,
            &second.decay.decays,
        )?;
        metric.reaction_counts(
            &coarse.produced_nuclei,
            &first.produced_nuclei,
            &second.produced_nuclei,
        )?;
        Ok(metric.error)
    }
}

impl AdaptiveNuclearHeatSource for DecayCellSource {
    fn time_state_error(
        &self,
        coarse: &Self,
        fine: &Self,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        same_binding(self, coarse, fine)?;
        same_model(
            self.inventory.model().as_ref(),
            coarse.inventory.model().as_ref(),
            fine.inventory.model().as_ref(),
        )?;
        let mut metric = Comparison::new(tolerances, checkpoint);
        metric.populations(
            self.inventory.populations(),
            coarse.inventory.populations(),
            fine.inventory.populations(),
        )?;
        Ok(metric.error)
    }
    fn time_step_error(
        coarse: &Self::Step,
        first: &Self::Step,
        second: &Self::Step,
        _: f64,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        let mut metric = Comparison::new(tolerances, checkpoint);
        common_steps::<Self>(&mut metric, coarse, first, second)?;
        metric.reaction_counts(&coarse.decays, &first.decays, &second.decays)?;
        metric.reaction_counts(
            &coarse.branch_decays,
            &first.branch_decays,
            &second.branch_decays,
        )?;
        Ok(metric.error)
    }
}

impl AdaptiveNuclearHeatSource for CaptureActivationCellSource {
    fn time_state_error(
        &self,
        coarse: &Self,
        fine: &Self,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        same_binding(self, coarse, fine)?;
        same_model(
            self.state.products().model().as_ref(),
            coarse.state.products().model().as_ref(),
            fine.state.products().model().as_ref(),
        )?;
        let mut metric = Comparison::new(tolerances, checkpoint);
        irradiation(
            &mut metric,
            self.state.parent(),
            coarse.state.parent(),
            fine.state.parent(),
        )?;
        metric.populations(
            self.state.products().populations(),
            coarse.state.products().populations(),
            fine.state.products().populations(),
        )?;
        Ok(metric.error)
    }
    fn time_step_error(
        coarse: &Self::Step,
        first: &Self::Step,
        second: &Self::Step,
        _: f64,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        let mut metric = Comparison::new(tolerances, checkpoint);
        common_steps::<Self>(&mut metric, coarse, first, second)?;
        metric.reaction_counts(
            &coarse.products.decay.decays,
            &first.products.decay.decays,
            &second.products.decay.decays,
        )?;
        metric.reaction_counts(
            &coarse.products.produced_nuclei,
            &first.products.produced_nuclei,
            &second.products.produced_nuclei,
        )?;
        Ok(metric.error)
    }
}

impl<S: AdaptiveNuclearHeatSource + FluxDrivenHeatSource> AdaptiveNuclearHeatSource
    for KineticCellSource<S>
{
    fn time_state_error(
        &self,
        coarse: &Self,
        fine: &Self,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        same_model(
            self.kinetics.model(),
            coarse.kinetics.model(),
            fine.kinetics.model(),
        )?;
        let mut error = self.irradiation.time_state_error(
            &coarse.irradiation,
            &fine.irradiation,
            tolerances,
            checkpoint,
        )?;
        let mut metric = Comparison::new(tolerances, checkpoint);
        metric.neutron(
            self.kinetics.neutron_population(),
            coarse.kinetics.neutron_population(),
            fine.kinetics.neutron_population(),
        )?;
        let a = self.kinetics.precursor_populations();
        let b = coarse.kinetics.precursor_populations();
        let c = fine.kinetics.precursor_populations();
        if a.len() != b.len() || a.len() != c.len() {
            return Err(FissionThermalError::Invalid(
                "refinement changed precursor ordering",
            ));
        }
        for ((&a, &b), &c) in a.iter().zip(b).zip(c) {
            metric.neutron(a, b, c)?;
        }
        error.merge(metric.error);
        Ok(error)
    }
    fn time_step_error(
        coarse: &Self::Step,
        first: &Self::Step,
        second: &Self::Step,
        duration_s: f64,
        tolerances: NuclearTimeTolerances,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<NuclearTimeError, FissionThermalError> {
        let mut error = S::time_step_error(
            &coarse.irradiation,
            &first.irradiation,
            &second.irradiation,
            duration_s,
            tolerances,
            checkpoint,
        )?;
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(FissionThermalError::Invalid(
                "positive refinement duration required",
            ));
        }
        // Compare fluence in mean-amplitude units using the existing neutron
        // integral, independently of whether any isotope currently absorbs it.
        let fine_integral = finite(
            first.kinetics.integrated_population_s + second.kinetics.integrated_population_s,
        )?;
        let mut metric = Comparison::new(tolerances, checkpoint);
        metric.neutron(
            0.0,
            finite(coarse.kinetics.integrated_population_s / duration_s)?,
            finite(fine_integral / duration_s)?,
        )?;
        error.merge(metric.error);
        Ok(error)
    }
}
