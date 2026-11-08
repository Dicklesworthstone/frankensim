//! Bounded first-order ghost paths and their joint path/pupil proposal.

use super::{CameraError, LensSample, SequentialLens};
use crate::polarization::{Mueller, Stokes};

#[cfg(test)]
mod tests;

/// A reciprocal lens path with exactly two interface reflections.
///
/// Surfaces are indexed sensor to scene. The first reflection must be at a
/// larger index than the second: the ray returns toward the sensor, reflects
/// again, then exits toward the scene. All other encounters transmit. These
/// are incoherent first-order ghost images, not a coherent bulk etalon or an
/// approximation to the sum of arbitrarily many internal round trips.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LensGhostPath {
    first: usize,
    second: usize,
}

impl LensGhostPath {
    /// Declare the first and second reflected faces in reciprocal traversal
    /// order. The consuming prescription also checks that both faces exist.
    pub fn try_new(first: usize, second: usize) -> Result<Self, CameraError> {
        if first >= SequentialLens::MAX_SURFACES || second >= first {
            return Err(CameraError::InvalidSequentialLens);
        }
        Ok(Self { first, second })
    }

    /// First reflected face, on the initial sensor-to-scene leg.
    #[must_use]
    pub const fn first_surface_index(self) -> usize {
        self.first
    }

    /// Second reflected face, on the return leg toward the sensor.
    #[must_use]
    pub const fn second_surface_index(self) -> usize {
        self.second
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct GhostSampling {
    direct_probability: f64,
    pair_count: usize,
}

impl GhostSampling {
    pub(super) fn try_new(surfaces: usize, direct_probability: f64) -> Result<Self, CameraError> {
        if !(2..=SequentialLens::MAX_SURFACES).contains(&surfaces)
            || !direct_probability.is_finite()
            || direct_probability <= 0.0
            || direct_probability >= 1.0
            || !(1.0 / direct_probability).is_finite()
        {
            return Err(CameraError::InvalidSequentialLens);
        }
        Ok(Self {
            direct_probability,
            pair_count: surfaces * (surfaces - 1) / 2,
        })
    }

    pub(super) const fn direct_probability(self) -> f64 {
        self.direct_probability
    }

    /// Partition u into path intervals and remap within the selected interval.
    /// The conditional pupil coordinate remains uniform on [0,1), independent
    /// of the path choice; v is untouched. No independent random dimension or
    /// path-dependent resampling is introduced. The returned probability is
    /// consumed exactly once by the camera's measure weight.
    pub(super) fn sample(
        self,
        sample: LensSample,
    ) -> Result<(Option<LensGhostPath>, LensSample, f64), CameraError> {
        let (path, local_u, probability) = if sample.u() < self.direct_probability {
            (
                None,
                sample.u() / self.direct_probability,
                self.direct_probability,
            )
        } else {
            let residual_probability = 1.0 - self.direct_probability;
            let count = self.pair_count as f64;
            let coordinate = ((sample.u() - self.direct_probability) / residual_probability
                * count)
                .min(count.next_down());
            let pair = coordinate as usize;
            let mut second = pair;
            let mut first = 1;
            while second >= first {
                second -= first;
                first += 1;
            }
            (
                Some(LensGhostPath { first, second }),
                coordinate - pair as f64,
                residual_probability / count,
            )
        };
        let pupil = LensSample::try_new(local_u.min(1.0_f64.next_down()), sample.v())?;
        Ok((path, pupil, probability))
    }
}

pub(super) struct PathStep {
    pub(super) surface: usize,
    pub(super) reverse: bool,
    pub(super) reflect: bool,
}

/// A direct path visits N faces; a pair (first,second) visits
/// N+2*(first-second), at most 3N-2. The path is structural, with no queue,
/// heap allocation, Russian roulette, or intensity-based truncation.
pub(super) struct PathSteps {
    surfaces: usize,
    next: usize,
    phase: u8,
    ghost: Option<LensGhostPath>,
}

impl PathSteps {
    pub(super) fn new(surfaces: usize, ghost: Option<LensGhostPath>) -> Self {
        Self {
            surfaces,
            next: 0,
            phase: 0,
            ghost,
        }
    }
}

impl Iterator for PathSteps {
    type Item = PathStep;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next >= self.surfaces {
            return None;
        }
        let surface = self.next;
        let reverse = self.phase == 1;
        let reflect = self.ghost.is_some_and(|ghost| {
            (self.phase == 0 && surface == ghost.first)
                || (self.phase == 1 && surface == ghost.second)
        });
        if reflect {
            self.phase += 1;
        }
        self.next = if self.phase == 1 {
            surface - 1
        } else {
            surface + 1
        };
        Some(PathStep {
            surface,
            reverse,
            reflect,
        })
    }
}

// Use exact zero support of the native diagonal s/p reflection operator.
// Checking input channels without halving avoids mistaking subnormal input
// underflow for a physically extinguished polarization.
pub(super) fn reflection_has_zero_support(matrix: Mueller, input: Stokes) -> bool {
    // The native Jones apply reconstructs a weak channel from the U/V block
    // even when mean +/- difference loses it near Brewster incidence. A
    // nonzero cross term therefore proves that neither channel is absent.
    if matrix.0[2][2] != 0.0 || matrix.0[3][2] != 0.0 {
        return false;
    }
    let difference = matrix.0[0][1];
    if matrix.0[0][0] + difference.abs() == 0.0 {
        return true;
    }
    if difference >= 0.0 {
        input.i == -input.q
    } else {
        input.i == input.q
    }
}
