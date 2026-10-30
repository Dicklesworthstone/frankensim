//! Accepted physical medium snapshots over an exposure interval.
//!
//! The timeline holds each immutable state over its declared interval. It does
//! not interpolate density, temperature, phase, coefficients or grid topology.
//! The renderer selects one state at the camera's absolute shutter time and
//! retains it for the complete path, including secondary and shadow rays.
//! This is a quasi-static exposure model, not finite-light-travel-time transport.

use super::AbsorbingGrid;
use super::scattering::ScatteringGrid;
use crate::motion::ShutterInterval;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;

/// An accepted optical snapshot whose simulation buffers remain borrowed.
#[derive(Clone, Copy)]
pub enum PhysicalMedium<'view, 'field> {
    /// Deterministic absorption and local-equilibrium emission.
    Absorbing(&'view AbsorbingGrid<'field>),
    /// Spectral absorption, emission and multiple scattering.
    Scattering(&'view ScatteringGrid<'field>),
}

impl PhysicalMedium<'_, '_> {
    /// Exact state/model identity, independent of the borrowed address.
    #[must_use]
    pub fn content_identity(self) -> ContentHash {
        match self {
            Self::Absorbing(medium) => medium.content_identity(),
            Self::Scattering(medium) => medium.content_identity(),
        }
    }

    /// Whether the accepted snapshot has zero extinction and emission.
    #[must_use]
    pub fn is_vacuum(self) -> bool {
        match self {
            Self::Absorbing(medium) => medium.is_vacuum(),
            Self::Scattering(medium) => medium.is_vacuum(),
        }
    }

    /// Whether the state has an absorption-linked positive-temperature source.
    #[must_use]
    pub fn has_thermal_emission(self) -> bool {
        match self {
            Self::Absorbing(medium) => medium.has_thermal_emission(),
            Self::Scattering(medium) => medium.has_thermal_emission(),
        }
    }
}

impl core::fmt::Debug for PhysicalMedium<'_, '_> {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let name = match self {
            Self::Absorbing(_) => "Absorbing",
            Self::Scattering(_) => "Scattering",
        };
        formatter
            .debug_tuple(name)
            .field(&self.content_identity())
            .finish()
    }
}

/// One accepted state held over a finite, nonzero physical interval.
///
/// Internal intervals are `[start_s, end_s)`: at a discontinuity the later
/// snapshot wins. Only the final history endpoint is inclusive, admitting a
/// zero-duration shutter at that endpoint without extrapolation.
#[derive(Clone, Copy, Debug)]
pub struct MediumFrame<'view, 'field> {
    /// Inclusive state start in absolute physical seconds.
    pub start_s: f64,
    /// Exclusive state end, except for the final history endpoint.
    pub end_s: f64,
    /// Immutable accepted optical state, borrowed without copying fields.
    pub medium: PhysicalMedium<'view, 'field>,
}

/// Caller-owned admission work limit; each physical state's own ray limits
/// remain binding when that state is rendered.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MediumTimelineBudget {
    /// Maximum number of accepted intervals inspected and hashed.
    pub max_frames: usize,
}

/// No timeline or exposure is admitted with missing or ambiguous time support.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TimelineError {
    /// The operation was cancelled.
    Cancelled,
    /// At least one finite interval is required.
    Empty,
    /// The interval scan exceeds the caller's budget.
    Budget,
    /// A state interval has nonfinite, equal or reversed endpoints.
    InvalidInterval {
        /// Invalid interval index.
        index: usize,
    },
    /// Adjacent intervals leave physical time without an accepted state.
    Gap {
        /// Index of the later interval.
        index: usize,
    },
    /// Adjacent intervals overlap or are out of order.
    Overlap {
        /// Index of the later interval.
        index: usize,
    },
    /// The requested sample time is nonfinite or outside the accepted history.
    TimeOutsideHistory,
    /// The complete shutter does not belong to the accepted history.
    ShutterOutsideHistory,
}

impl core::fmt::Display for TimelineError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "physical medium history refused: {self:?}")
    }
}

impl core::error::Error for TimelineError {}

/// Immutable, contiguous accepted history. The frame list and every underlying
/// density/temperature array remain borrowed; admission allocates nothing.
/// Content identity binds exact intervals, model kinds and accepted states.
#[derive(Debug)]
pub struct MediumTimeline<'frames, 'field> {
    frames: &'frames [MediumFrame<'frames, 'field>],
    identity: ContentHash,
    vacuum: bool,
    thermal_emission: bool,
}

impl<'frames, 'field> MediumTimeline<'frames, 'field> {
    /// Admit exact contiguous intervals once, checking cancellation per frame.
    ///
    /// A missing interval is never filled by clamping or an invented state.
    /// Successive endpoints must compare exactly equal; the unique physical
    /// zero is canonicalized for hashing. Different grid resolutions and
    /// optical models may meet at a declared discontinuity.
    pub fn try_new(
        cx: &Cx<'_>,
        frames: &'frames [MediumFrame<'frames, 'field>],
        budget: MediumTimelineBudget,
    ) -> Result<Self, TimelineError> {
        checkpoint(cx)?;
        if frames.is_empty() {
            return Err(TimelineError::Empty);
        }
        if frames.len() > budget.max_frames {
            return Err(TimelineError::Budget);
        }
        let mut hash = DomainHasher::new("org.frankensim.render.medium-timeline.v1");
        hash.update(&(frames.len() as u64).to_le_bytes());
        let mut vacuum = true;
        let mut thermal_emission = false;
        for (index, frame) in frames.iter().enumerate() {
            checkpoint(cx)?;
            if !frame.start_s.is_finite()
                || !frame.end_s.is_finite()
                || frame.start_s >= frame.end_s
            {
                return Err(TimelineError::InvalidInterval { index });
            }
            if index > 0 {
                let previous_end = frames[index - 1].end_s;
                if frame.start_s > previous_end {
                    return Err(TimelineError::Gap { index });
                }
                if frame.start_s < previous_end {
                    return Err(TimelineError::Overlap { index });
                }
            }
            hash.update(&canonical_zero(frame.start_s).to_bits().to_le_bytes());
            hash.update(&canonical_zero(frame.end_s).to_bits().to_le_bytes());
            hash.update(&[match frame.medium {
                PhysicalMedium::Absorbing(_) => 0,
                PhysicalMedium::Scattering(_) => 1,
            }]);
            hash.update(frame.medium.content_identity().as_bytes());
            vacuum &= frame.medium.is_vacuum();
            thermal_emission |= frame.medium.has_thermal_emission();
        }
        checkpoint(cx)?;
        Ok(Self {
            frames,
            identity: hash.finalize(),
            vacuum,
            thermal_emission,
        })
    }

    /// Exact history identity for progressive exposure binding.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Number of held accepted states.
    #[must_use]
    pub const fn frame_count(&self) -> usize {
        self.frames.len()
    }

    /// Inclusive first accepted physical time.
    #[must_use]
    pub fn start_s(&self) -> f64 {
        canonical_zero(self.frames[0].start_s)
    }

    /// Inclusive final history endpoint.
    #[must_use]
    pub fn end_s(&self) -> f64 {
        canonical_zero(self.frames[self.frames.len() - 1].end_s)
    }

    /// Every accepted state is a vacuum.
    #[must_use]
    pub const fn is_vacuum(&self) -> bool {
        self.vacuum
    }

    /// At least one accepted state has an absorption-linked thermal source.
    #[must_use]
    pub const fn has_thermal_emission(&self) -> bool {
        self.thermal_emission
    }

    /// Admit complete shutter coverage, including zero-duration endpoint shots.
    pub fn admit_shutter(
        &self,
        cx: &Cx<'_>,
        shutter: ShutterInterval,
    ) -> Result<(), TimelineError> {
        checkpoint(cx)?;
        if shutter.open_s() < self.start_s() || shutter.close_s() > self.end_s() {
            return Err(TimelineError::ShutterOutsideHistory);
        }
        Ok(())
    }

    /// Select the held accepted state in logarithmic time, without allocation.
    /// The returned borrow must remain fixed for the complete optical path.
    pub fn at_time(&self, time_s: f64) -> Result<PhysicalMedium<'frames, 'field>, TimelineError> {
        if !time_s.is_finite() || time_s < self.start_s() || time_s > self.end_s() {
            return Err(TimelineError::TimeOutsideHistory);
        }
        let index = self.frames.partition_point(|frame| frame.start_s <= time_s) - 1;
        Ok(self.frames[index].medium)
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), TimelineError> {
    cx.checkpoint().map_err(|_| TimelineError::Cancelled)
}

fn canonical_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

#[cfg(test)]
mod tests;
