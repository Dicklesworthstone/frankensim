//! Scheduled source-bound contact heating through the original pressure owner.
//! Thermal storage and constitutive feedback belong to the existing finite
//! enthalpy/contact transaction. This wrapper changes only which accepted step
//! the shared bow/force/observer loop calls; it never creates another voice.

use super::{MultiBowPerformance, render_invalid};
use crate::acoustic_realize::AcousticRealizeError;
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::render::schedule::force::coupled::ModalCouplingError;
use crate::render::schedule::force::coupled::contact::friction::multiple::layer::thermal::{
    LayerThermalCell, LayerThermalFrame, LayerThermalSolve, SeriesDriveHeatDestination,
};
use fs_exec::CancelGate;
use fs_thermal::enthalpy::{DepositionSimulation, backward_euler::StepCoefficients};
use std::marker::PhantomData;

/// A finite canonical performance and its original thermal network on one clock.
/// Source material/caloric cards remain borrowed, and the accepted network can
/// be recovered for continuation. No bath, partition fraction or thermal mass
/// is inferred. Each source-temperature change reaches the actual traction and
/// the original body/cavity or selected shell-pressure observer.
///
/// This inherits the fixed compressed-film, affine-modulus and solid-state
/// limits of `step_with_contact_thermal_network`. Moving stops/footprints refuse
/// because that owner does not account their thermal transport/work. Initial
/// material matching and all per-step domains remain checked by that owner.
pub struct ContactThermalPerformance<'p, 't, 'm, 'g, F, S, E> {
    performance: &'p mut MultiBowPerformance,
    thermal: DepositionSimulation<'t>,
    bindings: Vec<LayerThermalCell<'m>>,
    series_heat: Vec<SeriesDriveHeatDestination>,
    controls: LayerThermalSolve,
    gate: &'g CancelGate,
    resolve: F,
    energy_guess_j: Vec<f64>,
    last_frame: Option<LayerThermalFrame<S>>,
    error: PhantomData<fn() -> E>,
}

impl MultiBowPerformance {
    /// Render canonical bow/force controls with the existing contact/heat solve.
    /// `resolve` supplies OTHER heat and original-link conductances only; routed
    /// plastic/reversible/series-drive heat is already owned by the transaction.
    /// Thermal time must match the bank's absolute clock, including prior motion.
    /// Bindings and fixed thermal cell count are bounded by `controls`.
    #[allow(clippy::too_many_arguments)]
    pub fn with_contact_thermal_network<'p, 't, 'm, 'g, F, S, E>(
        &'p mut self,
        thermal: DepositionSimulation<'t>,
        bindings: Vec<LayerThermalCell<'m>>,
        series_heat: Vec<SeriesDriveHeatDestination>,
        controls: LayerThermalSolve,
        gate: &'g CancelGate,
        resolve: F,
    ) -> Result<ContactThermalPerformance<'p, 't, 'm, 'g, F, S, E>, RenderError>
    where
        F: Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        E: core::fmt::Display,
    {
        self.validate_sample_count(0)?;
        if gate.is_requested() {
            return Err(RenderError::Coupled(ModalCouplingError::Cancelled));
        }
        let system = &self.bank.system;
        if self.stop_motion.is_some()
            || self.stop_station_motion.is_some()
            || self.bank.bow_station_motion_limits().is_some()
            || thermal.time_s() != system.samples_rendered() as f64 * system.sample_period_s()
            || thermal.cells().is_empty()
            || thermal.cells().len() > controls.maximum_cells
            || (bindings.is_empty() && series_heat.is_empty())
        {
            return Err(render_invalid(
                "contact heat needs matching physical clocks, finite cells, explicit routes and fixed contact geometry",
            ));
        }
        let energy_guess_j = vec![0.0; thermal.cells().len()];
        Ok(ContactThermalPerformance {
            performance: self,
            thermal,
            bindings,
            series_heat,
            controls,
            gate,
            resolve,
            energy_guess_j,
            last_frame: None,
            error: PhantomData,
        })
    }
}

impl<F, S, E> ContactThermalPerformance<'_, '_, '_, '_, F, S, E> {
    /// Original physical performance with the accepted material and bow history.
    #[must_use]
    pub fn performance(&self) -> &MultiBowPerformance {
        self.performance
    }
    /// Original finite thermal state, never an observer-only temperature effect.
    #[must_use]
    pub fn thermal(&self) -> &DepositionSimulation<'_> {
        &self.thermal
    }
    /// Most recent accepted combined heat/work/pressure report.
    #[must_use]
    pub fn last_frame(&self) -> Option<&LayerThermalFrame<S>> {
        self.last_frame.as_ref()
    }
}
impl<'t, F, S, E> ContactThermalPerformance<'_, 't, '_, '_, F, S, E> {
    /// Release the performance borrow and recover all accepted heat history.
    pub fn into_thermal(self) -> DepositionSimulation<'t> {
        self.thermal
    }
}

impl<'t, F, S, E> PressureRenderer for ContactThermalPerformance<'_, 't, '_, '_, F, S, E>
where
    F: Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
    E: core::fmt::Display,
{
    fn samples_rendered(&self) -> u64 {
        self.performance.samples_rendered()
    }
    fn max_block_len(&self) -> usize {
        self.performance.max_block_len()
    }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.performance.validate_sample_rate(rate)
    }
    fn validate_sample_count(&self, count: u64) -> Result<(), RenderError> {
        self.performance.validate_sample_count(count)
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        let Self {
            performance,
            thermal,
            bindings,
            series_heat,
            controls,
            gate,
            resolve,
            energy_guess_j,
            last_frame,
            ..
        } = self;
        performance.block_with_step(output, |bank, forces, drives| {
            if thermal.cells().len() != energy_guess_j.len() {
                return Err(render_invalid(
                    "thermal topology changed during a fixed contact performance",
                ));
            }
            for (guess, cell) in energy_guess_j.iter_mut().zip(thermal.cells().values()) {
                *guess = cell.mass_kg * cell.specific_h_j_kg;
            }
            let accepted = bank
                .step_with_contact_thermal_network(
                    thermal,
                    bindings,
                    series_heat,
                    *controls,
                    energy_guess_j,
                    forces,
                    drives,
                    gate,
                    &*resolve,
                )
                .map_err(|error| {
                    RenderError::Voice(AcousticRealizeError::Nonlinear(error.to_string()))
                })?;
            let pressure = accepted.mechanics.observer_pressure_pa;
            *last_frame = Some(accepted);
            Ok(pressure)
        })
    }
}
