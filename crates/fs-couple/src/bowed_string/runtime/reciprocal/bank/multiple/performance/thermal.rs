//! Scheduled source-bound contact heating through the original pressure owner.
//! Thermal storage and constitutive feedback belong to the existing finite
//! enthalpy/contact transaction. This wrapper changes only which accepted step
//! the shared bow/force/observer loop calls; it never creates another voice.

use super::{MultiBowPerformance, render_invalid};
use super::super::shear_layer::dry::BowDryContactHeat;
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::render::schedule::force::coupled::ModalCouplingError;
use crate::render::schedule::force::coupled::contact::friction::multiple::FrictionPortFrame;
use crate::render::schedule::force::coupled::contact::friction::multiple::layer::thermal::{
    LayerThermalCell, LayerThermalFrame, LayerThermalSolve, SeriesDriveHeatDestination,
};
use fs_exec::CancelGate;
use fs_thermal::enthalpy::{DepositionSimulation, backward_euler::StepCoefficients};
use std::marker::PhantomData;

type TrialHeatSource<'p, 't, S, E> = dyn Fn(&DepositionSimulation<'t>, &[FrictionPortFrame])
    -> Result<StepCoefficients<S>, E> + Send + Sync + 'p;

/// A finite canonical performance and its original thermal network on one clock.
/// Source material/caloric cards remain borrowed, and the accepted network can
/// be recovered for continuation. No bath, partition fraction or thermal mass
/// is inferred. Each source-temperature change reaches the actual traction and
/// the original body/cavity or selected shell-pressure observer.
///
/// This inherits the fixed compressed-film, affine-modulus and solid-state
/// limits of the shared contact/enthalpy owner. Independent finger/fret supports
/// may press, release and migrate: their actuator work joins the SAME mechanical
/// transaction at every thermal trial. Film footprints and thermal cells stay
/// fixed; bow/film migration still needs transported material history. Initial
/// material matching and all per-step domains remain checked by that owner.
pub struct ContactThermalPerformance<'p, 't, 'm, 'g, F, S, E> {
    performance: &'p mut MultiBowPerformance,
    thermal: DepositionSimulation<'t>,
    bindings: Vec<LayerThermalCell<'m>>,
    series_heat: Vec<SeriesDriveHeatDestination>,
    dry: Vec<BowDryContactHeat<'m>>,
    controls: LayerThermalSolve,
    gate: &'g CancelGate,
    resolve: F,
    trial_resolve: Option<Box<TrialHeatSource<'p, 't, S, E>>>,
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
        if self.bank.bow_station_motion_limits().is_some()
            || thermal.time_s() != system.samples_rendered() as f64 * system.sample_period_s()
            || thermal.cells().is_empty()
            || thermal.cells().len() > controls.maximum_cells
            || (bindings.is_empty() && series_heat.is_empty())
        {
            return Err(render_invalid(
                "contact heat needs matching physical clocks, finite cells, explicit routes and fixed bow/film geometry",
            ));
        }
        let energy_guess_j = vec![0.0; thermal.cells().len()];
        Ok(ContactThermalPerformance {
            performance: self,
            thermal,
            bindings,
            series_heat,
            dry: Vec::new(),
            controls,
            gate,
            resolve,
            trial_resolve: None,
            energy_guess_j,
            last_frame: None,
            error: PhantomData,
        })
    }
}

impl MultiBowPerformance {
    /// Use the SAME canonical performance loop with moving material-coordinate
    /// dry-contact heat. Bow and stop programs, fiber support acceleration and
    /// pressure observations pass through the existing joint owner. This does
    /// not replace the source dry law with a fitted thermal friction curve.
    #[allow(clippy::too_many_arguments)]
    pub fn with_dry_contact_thermal_network<'p, 't, 'm, 'g, F, S, E>(
        &'p mut self, thermal: DepositionSimulation<'t>, dry: Vec<BowDryContactHeat<'m>>,
        series_heat: Vec<SeriesDriveHeatDestination>, controls: LayerThermalSolve,
        gate: &'g CancelGate, resolve: F,
    ) -> Result<ContactThermalPerformance<'p, 't, 'm, 'g, F, S, E>, RenderError>
    where F: Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>, E: core::fmt::Display {
        self.validate_sample_count(0)?;
        if gate.is_requested() { return Err(RenderError::Coupled(ModalCouplingError::Cancelled)); }
        let system = &self.bank.system;
        if dry.is_empty() || thermal.cells().is_empty() || thermal.cells().len() > controls.maximum_cells
            || thermal.time_s() != system.samples_rendered() as f64 * system.sample_period_s() {
            return Err(render_invalid("dry contact heat needs original finite cells, explicit routes and one physical clock"));
        }
        let energy_guess_j = vec![0.0; thermal.cells().len()];
        Ok(ContactThermalPerformance { performance: self, thermal, bindings: Vec::new(),
            series_heat, dry, controls, gate, resolve, trial_resolve: None,
            energy_guess_j, last_frame: None, error: PhantomData })
    }
}

impl<'p, 't, F, S, E> ContactThermalPerformance<'p, 't, '_, '_, F, S, E> {
    /// Replace the state-only source callback with a source that observes the
    /// same candidate's actual slip. Retained-film conduction can then move over
    /// finite counterface material without a commanded or lagged stroke cursor.
    /// Original film mass/history and the fixed mechanical patch are retained.
    pub fn with_trial_heat_source(
        mut self,
        resolve: impl Fn(&DepositionSimulation<'t>, &[FrictionPortFrame])
            -> Result<StepCoefficients<S>, E> + Send + Sync + 'p,
    ) -> Result<Self, RenderError> {
        if !self.dry.is_empty() || self.bindings.is_empty() {
            return Err(render_invalid("trial film heat requires the original retained-film thermal owner"));
        }
        self.trial_resolve = Some(Box::new(resolve));
        Ok(self)
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
            dry,
            controls,
            gate,
            resolve,
            trial_resolve,
            energy_guess_j,
            last_frame,
            ..
        } = self;
        performance.block_with_step(output, |bank, forces, drives, motion| {
            if thermal.cells().len() != energy_guess_j.len() {
                return Err(render_invalid(
                    "thermal topology changed during a fixed contact performance",
                ));
            }
            for (guess, cell) in energy_guess_j.iter_mut().zip(thermal.cells().values()) {
                *guess = cell.mass_kg * cell.specific_h_j_kg;
            }
            let accepted = if let Some(source) = trial_resolve.as_ref() {
                bank.step_with_contact_thermal_feedback_and_stops(thermal, bindings, series_heat,
                    *controls, energy_guess_j, forces, drives, motion.stop_stations,
                    motion.stop_offsets_m, gate, &**source)?
            } else if dry.is_empty() {
                bank.step_with_contact_thermal_network_and_stops(
                    thermal,
                    bindings,
                    series_heat,
                    *controls,
                    energy_guess_j,
                    forces,
                    drives,
                    motion.stop_stations,
                    motion.stop_offsets_m,
                    gate,
                    &*resolve,
                )?
            } else {
                bank.step_with_dry_contact_thermal_motion(thermal, dry, series_heat,
                    *controls, energy_guess_j, forces, drives, motion.bow_stations,
                    motion.stop_stations, motion.stop_offsets_m, gate, &*resolve)?
            };
            let pressure = accepted.mechanics.observer_pressure_pa;
            *last_frame = Some(accepted);
            Ok(pressure)
        })
    }
}
