//! Source-derived mechanics on the same finite scheduled reed/duct renderer.
//!
//! Structural reduction and its admission stay with `BeatingReedSource`.
//! This bridge adds no mechanical law, coefficient fitting or material lookup.
//! It preserves the source beside the renderer, so geometry, material identity,
//! retained modes and reduction evidence remain available to its caller.

use super::{GasFieldResolution, ReedPerformance, ReedPerformanceError};
use crate::pcm_wav::observation::{ObservationPolicy, PressureRenderer};
use crate::render::RenderError;
use fs_material::gas_field::AxialGasField;
use fs_material::gas_species::input::ThermochemicalSources;
use fs_scenario::{BeatingReed, reed_source::BeatingReedSource};

/// A structural admission and an acoustic/performance refusal keep their
/// original error types. Failed reduction never selects authored coefficients.
#[derive(Debug)]
pub enum SourceReedPerformanceError<E> {
    /// The physical source refused its requested reduction.
    Source(E),
    /// The shared performance owner refused the input, model or schedule.
    Performance(ReedPerformanceError),
}

impl<E: core::fmt::Display> core::fmt::Display for SourceReedPerformanceError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(error) => write!(f, "reed source: {error}"),
            Self::Performance(error) => error.fmt(f),
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for SourceReedPerformanceError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Source(error) => Some(error),
            Self::Performance(error) => Some(error),
        }
    }
}

/// A source-retaining finite performance using the existing massive reed, TMM
/// feedback, pressure gestures and observation owners.
///
/// The native input must explicitly select `reed source` in place of its six
/// numerical mechanics fields. Optional measured lay fields retain the same
/// syntax. The source supplies opening, width, closing pressure, mass,
/// stiffness and damping through its checked reduction; the gesture schedule
/// exclusively supplies blowing pressure and attack history.
/// `tongue_contact` refuses: its authored cantilever projection cannot represent
/// an arbitrary source mode. A source-projected contact map is required first.
///
/// This is the source's admitted single-mode approximation. It does not turn
/// that approximation into a resolved double reed or a measured instrument.
/// The inner performance's input hash identifies only its authored bytes, not
/// this separately retained structural source or its admission.
pub struct SourceReedPerformance<S: BeatingReedSource> {
    performance: ReedPerformance,
    source: S,
    admission: S::Admission,
    reed: BeatingReed,
}

impl<S: BeatingReedSource> SourceReedPerformance<S> {
    /// Compile the physical source and the complete native performance before
    /// advancing any sample. Gas declarations use the same source catalogue,
    /// local media and physical observation path as `ReedPerformance`.
    pub fn from_bytes(
        bytes: &[u8],
        max_block: usize,
        policy: ObservationPolicy,
        source: S,
        admission: S::Admission,
        gas_sources: Option<&ThermochemicalSources<'_>>,
    ) -> Result<Self, SourceReedPerformanceError<S::Error>> {
        Self::from_bytes_with_environment(
            bytes,
            max_block,
            policy,
            source,
            admission,
            gas_sources,
            None,
        )
    }

    /// Bind the source-derived mechanics to the same prescribed or converged
    /// axial gas field. Exact inlet gas drives the reed, and the shared field
    /// drives both duct feedback and physical terminal observation.
    #[allow(clippy::too_many_arguments)]
    pub fn from_bytes_with_field_resolution(
        bytes: &[u8],
        max_block: usize,
        policy: ObservationPolicy,
        source: S,
        admission: S::Admission,
        gas_sources: Option<&ThermochemicalSources<'_>>,
        field: AxialGasField,
        resolution: GasFieldResolution,
    ) -> Result<Self, SourceReedPerformanceError<S::Error>> {
        Self::from_bytes_with_environment(
            bytes,
            max_block,
            policy,
            source,
            admission,
            gas_sources,
            Some((field, resolution)),
        )
    }

    fn from_bytes_with_environment(
        bytes: &[u8],
        max_block: usize,
        policy: ObservationPolicy,
        source: S,
        admission: S::Admission,
        gas_sources: Option<&ThermochemicalSources<'_>>,
        gas_field: Option<(AxialGasField, GasFieldResolution)>,
    ) -> Result<Self, SourceReedPerformanceError<S::Error>> {
        let reed = source
            .compile_beating_reed(0.0, 0.0, &admission)
            .map_err(SourceReedPerformanceError::Source)?;
        let performance = ReedPerformance::from_bytes_with_environment_and_reed(
            bytes,
            max_block,
            policy,
            gas_sources,
            gas_field,
            Some(reed),
        )
        .map_err(SourceReedPerformanceError::Performance)?;
        Ok(Self {
            performance,
            source,
            admission,
            reed,
        })
    }

    /// The original source, retaining its physical evidence and identity.
    #[must_use]
    pub const fn source(&self) -> &S {
        &self.source
    }

    /// The explicit admission used by the original structural owner.
    #[must_use]
    pub const fn admission(&self) -> &S::Admission {
        &self.admission
    }

    /// Exactly the mechanics supplied to the shared massive-reed runtime.
    #[must_use]
    pub const fn reed(&self) -> BeatingReed {
        self.reed
    }

    /// Shared gas, schedule, observer, horizon and error receipts.
    #[must_use]
    pub const fn performance(&self) -> &ReedPerformance {
        &self.performance
    }

    /// Emit pressure on the original finite observation clock. Late physical
    /// failures preserve the shared renderer's poison-on-failure behavior.
    pub fn pressure_block(&mut self, output: &mut [f64]) -> Result<(), ReedPerformanceError> {
        self.performance.pressure_block(output)
    }
}

impl<S: BeatingReedSource> PressureRenderer for SourceReedPerformance<S> {
    fn samples_rendered(&self) -> u64 {
        self.performance.samples_rendered()
    }
    fn max_block_len(&self) -> usize {
        self.performance.max_block_len()
    }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.performance.validate_sample_rate(rate)
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        self.performance.validate_sample_count(samples)
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.performance.block(output)
    }
}
