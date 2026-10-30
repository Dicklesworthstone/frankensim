//! Bounded passive replication of an ideal paraxial projector light field.

use std::sync::Arc;

use crate::photometry::source::SourceError;

/// Maximum number of outgoing beams in one ideal projector splitter.
pub const MAX_PROJECTOR_BEAMS: usize = 32;

/// One achromatic outgoing beam of an ideal paraxial splitter.
///
/// The dimensionless offset adds to the outgoing slopes `w.x/w.z` and
/// `w.y/w.z` in the projector's optical frame. It is not an angle in radians.
/// At a real focus distance `q`, the projected image shifts by `q * offset`
/// metres. Throughput is a fraction of the incident ray power, not a relative
/// weight that will be normalized. Every beam retains the complete pupil.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProjectorBeam {
    slope_offset: [f64; 2],
    throughput: f64,
}

impl ProjectorBeam {
    /// Admit finite slope offsets and a passive throughput in `[0,1]`.
    pub fn try_new(slope_offset: [f64; 2], throughput: f64) -> Result<Self, SourceError> {
        if slope_offset.iter().any(|value| !value.is_finite()) {
            return Err(SourceError::InvalidDistribution);
        }
        if !throughput.is_finite() || !(0.0..=1.0).contains(&throughput) {
            return Err(SourceError::InvalidFilter);
        }
        Ok(Self {
            slope_offset,
            throughput,
        })
    }

    /// Dimensionless outgoing slope offsets in optical X and Y.
    #[must_use]
    pub const fn slope_offset(&self) -> [f64; 2] {
        self.slope_offset
    }

    /// Absolute fraction of incident ray power assigned to this beam.
    #[must_use]
    pub const fn throughput(&self) -> f64 {
        self.throughput
    }
}

/// Incoherent, achromatic replication of a projector's exit-pupil light field.
///
/// The sum of the supplied beam throughputs must not exceed one. Missing power
/// is absorbed; overlapping outgoing beams add without renormalization. This
/// is a declared ideal paraxial source model, not a solved refractive prism:
/// facet-dependent pupil partitions, dispersion, diffraction, internal
/// reflections and aberrations are not represented.
#[derive(Debug, Clone)]
pub struct ParaxialBeamSplitter {
    beams: Arc<[ProjectorBeam]>,
    transmitted_fraction: f64,
    maximum_abs_offset: [f64; 2],
}

impl ParaxialBeamSplitter {
    /// Admit one to [`MAX_PROJECTOR_BEAMS`] beams in the supplied stable order.
    /// All-zero throughputs explicitly represent an opaque splitter. Zero-power
    /// beams do not enlarge the source's angular support or radiance bound.
    pub fn try_new(beams: Vec<ProjectorBeam>) -> Result<Self, SourceError> {
        if beams.is_empty() {
            return Err(SourceError::InvalidDistribution);
        }
        if beams.len() > MAX_PROJECTOR_BEAMS {
            return Err(SourceError::DistributionTooLarge);
        }
        let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
        let mut maximum_abs_offset = [0.0_f64; 2];
        for beam in &beams {
            let value = beam.throughput;
            let next = sum + value;
            correction += if sum >= value {
                (sum - next) + value
            } else {
                (value - next) + sum
            };
            sum = next;
            if value > 0.0 {
                for (maximum, offset) in maximum_abs_offset.iter_mut().zip(beam.slope_offset) {
                    *maximum = maximum.max(offset.abs());
                }
            }
        }
        // Compare the retained residual too: rounding 1 + a tiny positive
        // throughput back to 1 must not silently admit an amplifying splitter.
        if correction > 1.0 - sum {
            return Err(SourceError::InvalidFilter);
        }
        let transmitted_fraction = sum + correction;
        if !transmitted_fraction.is_finite() || !(0.0..=1.0).contains(&transmitted_fraction) {
            return Err(SourceError::NumericalRange);
        }
        Ok(Self {
            beams: beams.into(),
            transmitted_fraction,
            maximum_abs_offset,
        })
    }

    /// Immutable beam data in the supplied deterministic summation order.
    #[must_use]
    pub fn beams(&self) -> &[ProjectorBeam] {
        &self.beams
    }

    /// Total surviving power fraction, independent of the projector's gobo.
    #[must_use]
    pub const fn transmitted_fraction(&self) -> f64 {
        self.transmitted_fraction
    }

    pub(super) const fn maximum_abs_offset(&self) -> [f64; 2] {
        self.maximum_abs_offset
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g0_splitter_preserves_absolute_power_and_admits_explicit_absorption() {
        let splitter = ParaxialBeamSplitter::try_new(vec![
            ProjectorBeam::try_new([-0.1, 0.0], 0.25).unwrap(),
            ProjectorBeam::try_new([0.2, 0.05], 0.5).unwrap(),
            ProjectorBeam::try_new([f64::MAX, -f64::MAX], 0.0).unwrap(),
        ])
        .unwrap();
        assert_eq!(splitter.transmitted_fraction(), 0.75);
        assert_eq!(splitter.maximum_abs_offset(), [0.2, 0.05]);
        assert_eq!(splitter.beams()[0].slope_offset(), [-0.1, 0.0]);
        assert_eq!(splitter.beams()[1].throughput(), 0.5);
        let opaque = ParaxialBeamSplitter::try_new(vec![
            ProjectorBeam::try_new([f64::MAX; 2], 0.0).unwrap(),
        ])
        .unwrap();
        assert_eq!(opaque.transmitted_fraction(), 0.0);
        assert_eq!(opaque.maximum_abs_offset(), [0.0; 2]);
    }

    #[test]
    fn g0_splitter_rejects_amplification_invalid_beams_and_unbounded_work() {
        assert!(ProjectorBeam::try_new([f64::NAN, 0.0], 1.0).is_err());
        for throughput in [-0.1, 1.01, f64::INFINITY] {
            assert!(ProjectorBeam::try_new([0.0; 2], throughput).is_err());
        }
        assert!(ParaxialBeamSplitter::try_new(Vec::new()).is_err());
        assert!(
            ParaxialBeamSplitter::try_new(vec![
                ProjectorBeam::try_new([0.0; 2], 0.0).unwrap();
                MAX_PROJECTOR_BEAMS + 1
            ])
            .is_err()
        );
        for weights in [[0.6, 0.6], [1.0, f64::MIN_POSITIVE]] {
            assert!(
                ParaxialBeamSplitter::try_new(
                    weights
                        .map(|weight| ProjectorBeam::try_new([0.0; 2], weight).unwrap())
                        .to_vec()
                )
                .is_err()
            );
        }
    }
}
