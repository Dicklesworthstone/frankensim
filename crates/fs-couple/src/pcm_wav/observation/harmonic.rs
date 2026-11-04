//! Physical harmonic pressure through the existing PCM streaming owner.
//!
//! This observes an already propagated, settled acoustic field. It adds no
//! source dynamics, attack envelope, normalization, or second WAV encoder.

use super::PressureRenderer;
use crate::render::RenderError;
use fs_aeroac::observation::{HarmonicSampleBudget, HarmonicSampleClock, PreparedHarmonicField};

mod channels;
pub use channels::{HarmonicOutputChannel, HarmonicOutputSource, HarmonicPressureFrameRenderer};

/// A fixed receiver's coherent harmonic pressure in physical pascals.
///
/// The original field is borrowed, retaining its source and propagation
/// assumptions. Construction allocates one bounded callback buffer. Successful
/// callbacks allocate nothing and their bits do not depend on block partition.
/// A failed callback leaves both the caller's output and sample counter intact.
pub struct HarmonicPressureRenderer<'a, const D: usize> {
    field: &'a PreparedHarmonicField<D>,
    clock: HarmonicSampleClock,
    budget: HarmonicSampleBudget,
    scratch: Vec<f64>,
    completed: u64,
}

impl<'a, const D: usize> HarmonicPressureRenderer<'a, D> {
    pub fn new(
        field: &'a PreparedHarmonicField<D>,
        clock: HarmonicSampleClock,
        budget: HarmonicSampleBudget,
    ) -> Result<Self, RenderError> {
        if budget.max_samples == 0
            || budget
                .max_samples
                .checked_mul(field.frequency_count())
                .is_none_or(|work| work > budget.max_tone_samples)
        {
            return Err(RenderError::Sizing {
                what: "harmonic callback exceeds its sample/work budget",
            });
        }
        field
            .validate_sample_range(clock, 0, 0)
            .map_err(RenderError::Acoustic)?;
        let mut scratch = Vec::new();
        scratch
            .try_reserve_exact(budget.max_samples)
            .map_err(|_| RenderError::Sizing {
                what: "harmonic callback allocation failed",
            })?;
        scratch.resize(budget.max_samples, 0.0);
        Ok(Self {
            field,
            clock,
            budget,
            scratch,
            completed: 0,
        })
    }

    pub const fn clock(&self) -> HarmonicSampleClock {
        self.clock
    }
    pub const fn field(&self) -> &PreparedHarmonicField<D> {
        self.field
    }
}

impl<const D: usize> PressureRenderer for HarmonicPressureRenderer<'_, D> {
    fn samples_rendered(&self) -> u64 {
        self.completed
    }
    fn max_block_len(&self) -> usize {
        self.budget.max_samples
    }

    fn validate_sample_rate(&self, sample_rate_hz: u32) -> Result<(), RenderError> {
        if sample_rate_hz != self.clock.sample_rate_hz() {
            return Err(RenderError::Sizing {
                what: "harmonic field and output sample clocks differ",
            });
        }
        Ok(())
    }

    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        self.field
            .validate_sample_range(self.clock, self.completed, samples)
            .map_err(RenderError::Acoustic)
    }

    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        if output.is_empty() {
            return Err(RenderError::EmptyBlock);
        }
        if output.len() > self.budget.max_samples {
            return Err(RenderError::Sizing {
                what: "harmonic output exceeds callback capacity",
            });
        }
        let count = u64::try_from(output.len()).map_err(|_| RenderError::Sizing {
            what: "harmonic output count exceeds u64",
        })?;
        self.validate_sample_count(count)?;
        let scratch = &mut self.scratch[..output.len()];
        self.field
            .sample_pressure_into(self.clock, self.completed, scratch, self.budget)
            .map_err(RenderError::Acoustic)?;
        output.copy_from_slice(scratch);
        self.completed += count;
        Ok(())
    }
}
