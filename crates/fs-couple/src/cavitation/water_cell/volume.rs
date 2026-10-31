//! Sealed equilibrium water driven by prescribed boundary volume and heat.
//!
//! The existing finite-water owner retains M and conserved U; the original
//! IF97 source recovers the endpoint at v=V/M. A single native scalar solve
//! enforces U1-U0 = Q - (p0+p1)(V1-V0)/2. The same mean ABSOLUTE pressure is
//! exposed to the mechanical owner, with positive work done by the water.
//! This trapezoidal pressure work is second order on smooth resolved paths;
//! phase crossings and source-domain changes require interval refinement.
//!
//! No wall momentum, compliance, imposed pressure law, kinetic energy,
//! evaporation rate, extra latent heat or liquid equation of state is invented.
//! This is a sealed, uniform, quasistatic equilibrium-volume step. The caller
//! supplies the actual mechanical volume and must use its returned work in
//! that same mechanical solve. A proposal can be inspected before publication.
//! Rigid-volume stream transport remains on the original `advance` method.
//!
//! The source image and its inversion defects remain distinct from conserved
//! U. Raw IF97 entropy changes include source seams and temporal error; without
//! a heat-boundary temperature they are not certified entropy production.
//! Source: IAPWS R7-97(2012), sections 5, 6, 8 and 10.2,
//! <https://iapws.org/technical-guidance/release/IF97-Rev.download>.

use super::{
    EquilibriumWaterCell, WaterCellError, WaterCellStep, budget, finite, image, poll, positive,
};
use crate::cavity_phs::piston::thermal::feedback;
use fs_material::phase::water::isochoric::IsochoricWater;
use fs_qty::Temperature;
use fs_thermochem::water::MassSpecificVolume;
use fs_thermochem::water::if97::saturation_pressure_derivative;

/// Physical interval bounds and the original bounded scalar-solve budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterVolumeControls {
    pub maximum_step_s: f64,
    /// In (0,0.1]; endpoint volume change divided by initial physical volume.
    pub maximum_relative_volume_change: f64,
    pub maximum_temperature_increment_k: f64,
    /// Source first-law root residual. Original cell energy/volume budgets
    /// additionally admit the final original-source image.
    pub energy_tolerance_j: f64,
    pub maximum_iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum WaterVolumeError {
    Water(WaterCellError),
    Invalid(&'static str),
    OutsideBracket { residual_j: [f64; 2] },
    NotConverged { iterations: usize, residual_j: f64 },
}
impl From<WaterCellError> for WaterVolumeError {
    fn from(value: WaterCellError) -> Self {
        Self::Water(value)
    }
}
impl core::fmt::Display for WaterVolumeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite water volume-work step: {self:?}")
    }
}
impl std::error::Error for WaterVolumeError {}

/// Original finite-water accounts with reciprocal mechanical work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterVolumeStep {
    pub water: WaterCellStep,
    /// (p_initial+p_endpoint)/2 in Pa. Its product with delta V is exactly
    /// the boundary_work_j used in the conserved water update.
    pub mean_boundary_pressure_pa: f64,
    pub volume_change_m3: f64,
    /// Source endpoint U minus conserved endpoint U, as evaluated by the
    /// scalar work solve. Kept separately from the actual stored-source delta.
    pub source_energy_root_residual_j: f64,
    /// Complete endpoint proposals, bounded by maximum_iterations + 2.
    pub source_trials: usize,
}

/// Complete candidate, borrowing the existing finite-water owner.
/// Dropping this value preserves its inventory, volume and accepted clock.
#[must_use = "inspect and commit to publish the complete volume-work interval"]
pub struct WaterVolumeCandidate<'a> {
    owner: &'a mut EquilibriumWaterCell,
    step: WaterVolumeStep,
}
impl WaterVolumeCandidate<'_> {
    pub const fn step(&self) -> &WaterVolumeStep {
        &self.step
    }
    pub const fn state(&self) -> &super::WaterCellState {
        &self.step.water.state
    }
    /// Infallible publication after downstream mechanical work has been admitted.
    pub fn commit(self) -> WaterVolumeStep {
        self.owner.state = self.step.water.state;
        self.step
    }
}

impl EquilibriumWaterCell {
    /// Solve and publish a sealed prescribed-volume interval atomically.
    pub fn advance_volume(
        &mut self,
        volume_m3: f64,
        supplied_heat_j: f64,
        end_time_s: f64,
        controls: WaterVolumeControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WaterVolumeStep, WaterVolumeError> {
        Ok(self
            .propose_volume(
                volume_m3,
                supplied_heat_j,
                end_time_s,
                controls,
                &mut cancelled,
            )?
            .commit())
    }

    /// Prepare a sealed interval on the same native M,U/source state.
    ///
    /// Volume/temperature/time ceilings are admission limits, not estimates of
    /// truncation error. Resolve the prescribed path and phase crossings by
    /// interval refinement. The caller's mechanics must account for the returned
    /// signed work before committing. Callback side effects are not rolled back.
    pub fn propose_volume(
        &mut self,
        volume_m3: f64,
        supplied_heat_j: f64,
        end_time_s: f64,
        controls: WaterVolumeControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WaterVolumeCandidate<'_>, WaterVolumeError> {
        poll(&mut cancelled)?;
        let initial = self.state;
        let duration = end_time_s - initial.time_s();
        let dv = volume_m3 - initial.volume_m3();
        if !positive(volume_m3)
            || !supplied_heat_j.is_finite()
            || !end_time_s.is_finite()
            || !positive(duration)
            || !positive(controls.maximum_step_s)
            || duration > controls.maximum_step_s
            || !positive(controls.maximum_relative_volume_change)
            || controls.maximum_relative_volume_change > 0.1
            || !positive(controls.maximum_temperature_increment_k)
            || !positive(controls.energy_tolerance_j)
            || !(1..=128).contains(&controls.maximum_iterations)
            || !dv.is_finite()
            || dv.abs() / initial.volume_m3() > controls.maximum_relative_volume_change
        {
            return Err(WaterVolumeError::Invalid(
                "finite forward interval and bounded physical volume/temperature/energy controls required",
            ));
        }
        if dv == 0.0 {
            let water = self.propose(&[], 0.0, supplied_heat_j, end_time_s, &mut cancelled)?;
            if (water.state.source_state().temperature_k() - initial.source_state().temperature_k())
                .abs()
                > controls.maximum_temperature_increment_k
            {
                return Err(WaterVolumeError::Invalid(
                    "water temperature increment exceeded",
                ));
            }
            budget(
                "water volume-work energy",
                water.energy_residual_j,
                controls.energy_tolerance_j,
            )?;
            let pressure = 0.5 * initial.source_state().pressure_pa()
                + 0.5 * water.state.source_state().pressure_pa();
            poll(&mut cancelled)?;
            return Ok(WaterVolumeCandidate {
                owner: self,
                step: WaterVolumeStep {
                    mean_boundary_pressure_pa: pressure,
                    volume_change_m3: 0.0,
                    source_trials: 0,
                    source_energy_root_residual_j: water.state.source_energy_residual_j(),
                    water,
                },
            });
        }
        let mass = initial.mass_kg();
        let specific_volume = volume_m3 / mass;
        if !positive(specific_volume) {
            return Err(WaterVolumeError::Invalid(
                "specific volume is unrepresentable",
            ));
        }
        let law = IsochoricWater::new(
            MassSpecificVolume::new(specific_volume),
            self.controls.source_controls(mass),
            &mut cancelled,
        )
        .map_err(WaterCellError::from)?;
        let t0 = initial.source_state().temperature_k();
        let u0 = initial.source_state().specific_internal_energy_j_kg();
        let p0 = initial.source_state().pressure_pa();
        let lo = (law.controls().minimum_temperature_k - t0)
            .max(-controls.maximum_temperature_increment_k);
        let hi = (law.controls().maximum_temperature_k - t0)
            .min(controls.maximum_temperature_increment_k);
        // Solve in a temperature increment so the original zero-containing
        // bracket owner is reused, including its bounded work and cancellation.
        let tolerance = controls
            .energy_tolerance_j
            .min(self.controls.energy_tolerance_j)
            * 0.25;
        if !positive(tolerance) {
            return Err(WaterVolumeError::Invalid(
                "root energy tolerance underflows",
            ));
        }
        let mut pressure_iterations = 0;
        let (source, mean_pressure, work, source_residual, trials) = feedback::solve_scalar(
            [lo, hi],
            tolerance,
            controls.maximum_iterations,
            |dt, trials, cancelled| {
                let temperature = t0 + dt;
                if dt != 0.0 && temperature == t0 {
                    return Err(WaterVolumeError::Invalid(
                        "temperature increment is unrepresentable",
                    ));
                }
                let source = law
                    .state_at_temperature(Temperature::new(temperature), &mut *cancelled)
                    .map_err(WaterCellError::from)?;
                pressure_iterations += source.pressure_iterations();
                let mean_pressure = finite(
                    0.5 * p0 + 0.5 * source.pressure_pa(),
                    "mean boundary pressure",
                )?;
                let work = finite(mean_pressure * dv, "water pressure-volume work")?;
                if work == 0.0 {
                    return Err(WaterVolumeError::Invalid(
                        "nonzero boundary work underflows",
                    ));
                }
                // Original du/dT at fixed v plus the endpoint pressure-work
                // derivative. An unstable sampled root branch is refused;
                // this local check is not a global monotonicity enclosure.
                let pressure_slope = if let Some(single) = source.single_phase_source() {
                    -single.dv_dt_at_pressure().value() / single.dv_dp_at_temperature().value()
                } else {
                    saturation_pressure_derivative(Temperature::new(temperature))
                        .map_err(|error| WaterCellError::Source(error.into()))?
                        .value()
                };
                let slope =
                    mass * source.heat_capacity_at_volume_j_kg_k() + 0.5 * dv * pressure_slope;
                if !positive(slope) {
                    return Err(WaterVolumeError::Invalid(
                        "nonpositive local water/work caloric slope",
                    ));
                }
                let change = finite(
                    mass * (source.specific_internal_energy_j_kg() - u0),
                    "source water energy change",
                )?;
                let residual = finite(
                    (change + initial.source_energy_residual_j()) + work - supplied_heat_j,
                    "source volume-work balance",
                )?;
                Ok(((source, mean_pressure, work, residual, trials), residual))
            },
            &mut cancelled,
        )
        .map_err(|error| match error {
            feedback::ScalarError::Evaluation(e) => e,
            feedback::ScalarError::Invalid(s) => WaterVolumeError::Invalid(s),
            feedback::ScalarError::Cancelled => WaterCellError::Cancelled.into(),
            feedback::ScalarError::OutsideBracket { lower, upper } => {
                WaterVolumeError::OutsideBracket {
                    residual_j: [lower, upper],
                }
            }
            feedback::ScalarError::NotConverged {
                iterations,
                residual,
            } => WaterVolumeError::NotConverged {
                iterations,
                residual_j: residual,
            },
        })?;
        let increment = finite(supplied_heat_j - work, "conserved water energy increment")?;
        let energy = finite(
            initial.internal_energy_j() + increment,
            "conserved water energy",
        )?;
        if increment != 0.0 && energy == initial.internal_energy_j() {
            return Err(WaterVolumeError::Invalid(
                "conserved water energy increment is unrepresentable",
            ));
        }
        let state = image(volume_m3, mass, energy, end_time_s, source, self.controls)?;
        let energy_change = state.internal_energy_j() - initial.internal_energy_j();
        let energy_residual = finite(
            energy_change + work - supplied_heat_j,
            "conserved water work balance",
        )?;
        budget(
            "conserved water volume-work",
            energy_residual,
            controls.energy_tolerance_j,
        )?;
        let source_change = state.source_internal_energy_j() - initial.source_internal_energy_j();
        let actual_source_residual = finite(
            source_change + work - supplied_heat_j,
            "source water work balance",
        )?;
        let entropy_change = finite(
            mass * (source.specific_entropy_j_kg_k()
                - initial.source_state().specific_entropy_j_kg_k()),
            "source water entropy change",
        )?;
        // This is the actual stored-source delta balance. It differs from the
        // root's balance by the initial source inversion defect and roundoff.
        let water = WaterCellStep {
            initial,
            state,
            incoming_mass_kg: 0.0,
            incoming_enthalpy_j: 0.0,
            outgoing: None,
            supplied_heat_j,
            boundary_work_j: work,
            mass_change_kg: 0.0,
            internal_energy_change_j: energy_change,
            mass_residual_kg: 0.0,
            energy_residual_j: energy_residual,
            source_first_law_residual_j: actual_source_residual,
            source_entropy_change_j_k: entropy_change,
            advected_entropy_j_k: 0.0,
            unaccounted_entropy_change_j_k: entropy_change,
            // Counts queried endpoint inversions; original chart preparation
            // has its own bounded pair of pressure inversions.
            temperature_iterations: trials,
            pressure_iterations,
        };
        poll(&mut cancelled)?;
        Ok(WaterVolumeCandidate {
            owner: self,
            step: WaterVolumeStep {
                water,
                mean_boundary_pressure_pa: mean_pressure,
                volume_change_m3: dv,
                source_trials: trials,
                source_energy_root_residual_j: source_residual,
            },
        })
    }
}

#[cfg(test)]
#[path = "volume/tests.rs"]
mod tests;
