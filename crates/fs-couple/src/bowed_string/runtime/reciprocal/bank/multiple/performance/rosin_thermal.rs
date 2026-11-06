//! Smith & Woodhouse (2000) thermal rosin on the canonical scheduled bow loop.
//! Every friction port of the bank is one cell of the same rosin contact
//! (`contact::friction::multiple::rosin_thermal`). Before each joint step every
//! cell takes its step-start Coulomb coefficient; the ordinary step owner then
//! solves the sample, and only the step it accepts advances the heat. Schedule,
//! forces, motion, observers and the PCM route are the original performance's.
//! A later refusal on the same sample (an observer or coating check) poisons
//! the performance as it does without rosin; it is not a rollback. The
//! performance is borrowed, so its owner keeps it, with the last coefficients
//! on its ports, for the receipts after the render.

use super::{MultiBowPerformance, ordinary_step, render_invalid};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::render::schedule::force::coupled::contact::friction::multiple::rosin_thermal::RosinThermalPorts;
use fs_tribo::rosin_thermal::SmithWoodhouseParams;

/// A scheduled bow performance whose contacts carry thermal rosin.
pub struct RosinThermalPerformance<'p> {
    performance: &'p mut MultiBowPerformance,
    rosin: RosinThermalPorts,
}

impl MultiBowPerformance {
    /// Bind every bow contact to one thermal rosin contact before the first
    /// sample. `params.contact` is the WHOLE contact: its radius at the whole
    /// bow load. Every contact must be a plain Coulomb port; the law's step
    /// must be the bank's sample period. A refusal leaves the performance as it was.
    pub fn with_rosin_thermal(&mut self, params: SmithWoodhouseParams)
        -> Result<RosinThermalPerformance<'_>, RenderError>
    {
        self.validate_sample_count(0)?;
        if self.samples_rendered() != 0 {
            return Err(render_invalid("thermal rosin binds before the performance's first sample"));
        }
        let ports = (0..self.bank.system.port_count()).collect();
        let rosin = RosinThermalPorts::new(params, ports, &mut self.bank.system).map_err(RenderError::Friction)?;
        Ok(RosinThermalPerformance { performance: self, rosin })
    }
}

impl RosinThermalPerformance<'_> {
    /// The borrowed performance, with its accepted mechanics.
    #[must_use] pub fn performance(&self) -> &MultiBowPerformance { self.performance }
    /// Per-cell temperatures, coefficients and the heat-step receipt.
    #[must_use] pub const fn rosin(&self) -> &RosinThermalPorts { &self.rosin }
    /// End the binding and keep the rosin receipt.
    #[must_use] pub fn into_rosin(self) -> RosinThermalPorts { self.rosin }
}

impl PressureRenderer for RosinThermalPerformance<'_> {
    fn samples_rendered(&self) -> u64 { self.performance.samples_rendered() }
    fn max_block_len(&self) -> usize { self.performance.max_block_len() }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.performance.validate_sample_rate(rate)
    }
    fn validate_sample_count(&self, count: u64) -> Result<(), RenderError> {
        self.performance.validate_sample_count(count)
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        let Self { performance, rosin } = self;
        performance.block_with_step(output, |bank, forces, drives, motion| {
            rosin.apply(&mut bank.system).map_err(RenderError::Friction)?;
            let frame = ordinary_step(bank, forces, drives, motion)?;
            rosin.accept(frame).map_err(RenderError::Friction)?;
            Ok(frame.observer_pressure_pa)
        })
    }
}
