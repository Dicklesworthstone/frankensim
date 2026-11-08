//! Repeated conserved U,V increments from one original phase-equilibrium state.
//!
//! The native energy/volume owner retains both conservation residuals and its
//! pressure/PH bracketing policy. Preparation only avoids selecting the same
//! reference again; it never substitutes a caloric table or chains trial states.

use super::{EnergyTarget, PhaseEnergyVolumeControls, PhaseEnergyVolumeSolution, Sample, Work};
use crate::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use crate::thermochemical_gas::nonideal::flash::{
    CaloricSample, E, enthalpy::spacing, poll, positive,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::eos::peng_robinson::phase_screen::PengRobinsonPhaseScreenV1;

/// Immutable original reference for repeated internal-energy/volume queries.
///
/// The original mixture, phase policy, caloric interval, pressure bracket and
/// tolerances remain fixed. Every query starts from this reference, even if a
/// previous query used another volume. Charge preparation once to an enclosing
/// operation; query counters include only newly attempted work.
#[derive(Clone, Debug)]
pub struct PreparedPhaseEnergyVolumePath<'m> {
    model: &'m ThermochemicalPengRobinsonModel,
    initial: Sample,
    controls: PhaseEnergyVolumeControls,
    preparation_pt_evaluations: usize,
    preparation_eos_evaluation_upper_bound: usize,
}

impl ThermochemicalPengRobinsonModel {
    /// Select one original reference for repeated conserved U,V increments.
    ///
    /// This performs one native PT selection and admits the original caloric
    /// interval once. No process heat, boundary work or entropy law is inferred.
    /// The prepared pressure and PT caps apply separately to each query, which
    /// can also supply smaller remaining allowances for an enclosing operation.
    pub fn prepare_phase_energy_volume_path(
        &self,
        initial_temperature: Temperature,
        initial_pressure: Pressure,
        controls: PhaseEnergyVolumeControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedPhaseEnergyVolumePath<'_>, E> {
        controls.validate()?;
        if !positive(initial_pressure.value()) {
            return Err(E::Invalid(
                "positive prepared U,V reference pressure required",
            ));
        }
        self.admit_enthalpy_flash(
            controls.pressure_bracket[0],
            controls.phase.enthalpy,
            Some(initial_temperature),
            &mut cancelled,
        )?;
        let mut work = Work::default();
        let initial = self.energy_volume_sample(
            initial_temperature,
            initial_pressure,
            controls,
            &mut work,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(PreparedPhaseEnergyVolumePath {
            model: self,
            initial,
            controls,
            preparation_pt_evaluations: work.pt,
            preparation_eos_evaluation_upper_bound: work.eos_bound,
        })
    }
}

impl PreparedPhaseEnergyVolumePath<'_> {
    /// Original PT selections performed during preparation.
    #[must_use]
    pub const fn preparation_pt_evaluations(&self) -> usize {
        self.preparation_pt_evaluations
    }

    /// Preparation's original EOS work ceiling, not a measured EOS count.
    #[must_use]
    pub const fn preparation_eos_evaluation_upper_bound(&self) -> usize {
        self.preparation_eos_evaluation_upper_bound
    }

    /// Original mass-specific reference volume [m³/kg], without another flash.
    #[must_use]
    pub fn reference_specific_volume_m3_kg(&self) -> f64 {
        self.reference_sample().molar_volume / self.model.mechanics().molar_mass().value()
    }

    // Composition-changing native owners cancel species datums against this
    // original caloric sample. It cannot be forged through the public API.
    pub(in crate::thermochemical_gas::nonideal::flash) fn reference_sample(
        &self,
    ) -> &CaloricSample<PengRobinsonPhaseScreenV1> {
        &self.initial
    }

    /// Resolve a new volume and original-reference internal-energy increment.
    ///
    /// Both allowances bound new work and are capped by the prepared controls.
    /// Exact zero energy at unchanged volume needs no new work when the reference
    /// lies in both brackets, so it remains available with zero allowances.
    /// A volume change with zero energy still requires the native U,V solve.
    ///
    /// Failure or cancellation leaves this path reusable. The original
    /// source-integrated energy and entropy increments retain their species
    /// datum cancellation, including during phase conversion.
    pub fn solve_increment(
        &self,
        final_specific_volume_m3_kg: f64,
        change_j_kg: f64,
        remaining_pressure_evaluations: usize,
        remaining_pt_evaluations: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseEnergyVolumeSolution, E> {
        poll(&mut cancelled)?;
        if !positive(final_specific_volume_m3_kg)
            || spacing(final_specific_volume_m3_kg) > self.controls.specific_volume_tolerance_m3_kg
            || !change_j_kg.is_finite()
        {
            return Err(E::Invalid(
                "representable positive prepared U,V volume and finite energy increment required",
            ));
        }
        let temperature = self.initial.temperature.value();
        let pressure = self.initial.pressure.value();
        let temperature_bracket = self.controls.phase.enthalpy.temperature_bracket;
        let pressure_bracket = self.controls.pressure_bracket;
        let exact_rest = change_j_kg == 0.0
            && final_specific_volume_m3_kg == self.reference_specific_volume_m3_kg()
            && temperature >= temperature_bracket[0].value()
            && temperature <= temperature_bracket[1].value()
            && pressure >= pressure_bracket[0].value()
            && pressure <= pressure_bracket[1].value();
        if remaining_pressure_evaluations == 0 && !exact_rest {
            return Err(E::PressureWorkLimit {
                pressure_evaluations: 0,
            });
        }
        if remaining_pt_evaluations == 0 && !exact_rest {
            return Err(E::WorkLimit {
                flash_evaluations: 0,
            });
        }
        let mut controls = self.controls;
        // Retain valid controls for exact rest, which makes no pressure or PT
        // query. Every other solve is bounded by the actual remaining allowance.
        controls.maximum_pressure_evaluations = controls
            .maximum_pressure_evaluations
            .min(remaining_pressure_evaluations.max(1));
        controls.maximum_pt_evaluations = controls
            .maximum_pt_evaluations
            .min(remaining_pt_evaluations.max(1));
        self.model.solve_energy_volume(
            final_specific_volume_m3_kg,
            EnergyTarget::Increment {
                initial: self.initial.clone(),
                requested: change_j_kg,
            },
            controls,
            Work::default(),
            &mut cancelled,
        )
    }
}
