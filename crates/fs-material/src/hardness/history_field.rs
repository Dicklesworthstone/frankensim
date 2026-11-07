//! Compact source-bound hardness transferred from accepted native histories.
//!
//! One fixed native material/prepared hardness owner supplies an ordered field.
//! Actual law, strain and complete history are checked at every point by the
//! original J2/Voce/Chaboche accepted observer. Each point retains only alpha
//! and H; original mechanical and indentation receipts are shared once. There
//! is no per-point database query, return map or source-receipt allocation.
//!
//! This transfers an independently sourced property of an externally evolved
//! or preconditioned bulk history. It does not transfer residual stress, plastic
//! displacements, geometry or elastic compliance. A contact consumer retains
//! its original geometry, homogeneous half-space assumptions and explicit
//! cell-to-point map; contact depth and wear recession never stand in for alpha.
//! Bulk indentation plasticity and a surface indentation approximation must
//! not own the same deformation or its work twice. The ordered inputs are
//! caller-declared endpoints, not an attestation of a trajectory or mesh origin.

use std::sync::Arc;

use fs_blake3::{ContentHash, DomainHasher};

use super::IndentationHardnessStatePoint;
use super::cyclic_history::{CyclicHardnessError, ResolvedCyclicHardness, SourceBoundHardnessLaw};
use super::plastic_history::ResolvedJ2Hardness;
use crate::Voigt;
use crate::plastic::J2State;
use crate::plastic::source::ResolvedJ2Plasticity;
use crate::state_point::ResolvedMaterialStatePoint;

/// Numeric observation admitted from one complete native accepted endpoint.
/// Private values cannot be changed independently of the original source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcceptedHardnessSample {
    alpha: f64,
    hardness_pa: f64,
}

impl AcceptedHardnessSample {
    /// Native accumulated equivalent plastic strain, never signed total strain.
    #[must_use]
    pub const fn alpha(&self) -> f64 {
        self.alpha
    }

    /// Positive original force/area indentation hardness [Pa].
    #[must_use]
    pub const fn hardness_pa(&self) -> f64 {
        self.hardness_pa
    }
}

#[derive(Clone, Debug, PartialEq)]
struct Sources {
    identity: ContentHash,
    mechanical: ResolvedMaterialStatePoint,
    hardness: IndentationHardnessStatePoint,
    maximum_alpha: f64,
}

/// Complete ordered numeric field and one shared original source descriptor.
/// All points use one fixed native source, with independently retained alpha.
#[derive(Clone, Debug, PartialEq)]
pub struct AcceptedHardnessField {
    sources: Arc<Sources>,
    samples: Vec<AcceptedHardnessSample>,
    identity: ContentHash,
}

/// A refusal returns no partial field and leaves all input histories untouched.
#[derive(Clone, Debug, PartialEq)]
pub enum HardnessFieldError {
    Invalid(&'static str),
    Point {
        point: usize,
        error: CyclicHardnessError,
    },
    /// The caller's checkpoint refused continued construction.
    Interrupted,
}

impl core::fmt::Display for HardnessFieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "accepted hardness field: {reason}"),
            Self::Point { point, error } => write!(f, "accepted hardness point {point}: {error}"),
            Self::Interrupted => f.write_str("accepted hardness field interrupted"),
        }
    }
}
impl std::error::Error for HardnessFieldError {}

impl AcceptedHardnessField {
    /// Observe actual native J2 law/strain/state tuples in the retained point
    /// order. Every actual law must match the one prepared source exactly.
    pub fn from_j2<'a>(
        source: &ResolvedJ2Hardness,
        points: impl ExactSizeIterator<Item = (&'a ResolvedJ2Plasticity, &'a Voigt, &'a J2State)>,
        yield_tolerance_pa: f64,
    ) -> Result<Self, HardnessFieldError> {
        Self::from_j2_checked(source, points, yield_tolerance_pa, || Ok(()))
    }

    /// J2 observation with a caller-owned checkpoint before allocation, at
    /// least every 256 points and before publication. The material layer owns
    /// no execution runtime; a native spatial owner supplies its original Cx.
    pub fn from_j2_checked<'a>(
        source: &ResolvedJ2Hardness,
        points: impl ExactSizeIterator<Item = (&'a ResolvedJ2Plasticity, &'a Voigt, &'a J2State)>,
        yield_tolerance_pa: f64,
        checkpoint: impl FnMut() -> Result<(), HardnessFieldError>,
    ) -> Result<Self, HardnessFieldError> {
        Self::build(
            source.identity(),
            source.plastic().material(),
            source.hardness_source(),
            source.maximum_accumulated_plastic_strain(),
            points.map(|(native, strain, state)| {
                let hardness_pa =
                    source.hardness_pa_at_accepted(native, strain, state, yield_tolerance_pa)?;
                Ok(AcceptedHardnessSample {
                    alpha: state.alpha,
                    hardness_pa,
                })
            }),
            checkpoint,
        )
    }

    /// Observe complete actual Voce or Chaboche tuples in their retained point
    /// order, including all branches and the native nonlinear yield surface.
    pub fn from_cyclic<'a, L: SourceBoundHardnessLaw + 'a>(
        source: &ResolvedCyclicHardness<L>,
        points: impl ExactSizeIterator<Item = (&'a L, &'a Voigt, &'a L::State)>,
        yield_tolerance_pa: f64,
    ) -> Result<Self, HardnessFieldError>
    where
        L::State: 'a,
    {
        Self::from_cyclic_checked(source, points, yield_tolerance_pa, || Ok(()))
    }

    /// Cyclic observation with the same caller-owned checkpoint protocol.
    pub fn from_cyclic_checked<'a, L: SourceBoundHardnessLaw + 'a>(
        source: &ResolvedCyclicHardness<L>,
        points: impl ExactSizeIterator<Item = (&'a L, &'a Voigt, &'a L::State)>,
        yield_tolerance_pa: f64,
        checkpoint: impl FnMut() -> Result<(), HardnessFieldError>,
    ) -> Result<Self, HardnessFieldError>
    where
        L::State: 'a,
    {
        Self::build(
            source.identity(),
            source.plastic().hardness_backbone().material(),
            source.hardness_source(),
            source.maximum_accumulated_plastic_strain(),
            points.map(|(native, strain, state)| {
                let (alpha, hardness_pa) =
                    source.accepted_values(native, strain, state, yield_tolerance_pa)?;
                Ok(AcceptedHardnessSample { alpha, hardness_pa })
            }),
            checkpoint,
        )
    }

    fn build(
        source_identity: ContentHash,
        mechanical: &ResolvedMaterialStatePoint,
        hardness: &IndentationHardnessStatePoint,
        maximum_alpha: f64,
        points: impl ExactSizeIterator<Item = Result<AcceptedHardnessSample, CyclicHardnessError>>,
        mut checkpoint: impl FnMut() -> Result<(), HardnessFieldError>,
    ) -> Result<Self, HardnessFieldError> {
        checkpoint()?;
        let count = points.len();
        if count == 0 {
            return Err(HardnessFieldError::Invalid(
                "a nonempty ordered native point field is required",
            ));
        }
        let mut samples = Vec::new();
        samples.try_reserve_exact(count).map_err(|_| {
            HardnessFieldError::Invalid(
                "native hardness point field exceeds representable allocation",
            )
        })?;
        let mut identity =
            DomainHasher::new("org.frankensim.fs-material.accepted-hardness-field.v1");
        identity.update(source_identity.as_bytes());
        identity.update(&(count as u64).to_le_bytes());
        for (point, sample) in points.enumerate() {
            if point % 256 == 0 {
                checkpoint()?;
            }
            if point >= count {
                return Err(HardnessFieldError::Invalid(
                    "native point count changed during observation",
                ));
            }
            let sample = sample.map_err(|error| HardnessFieldError::Point { point, error })?;
            identity.update(&sample.alpha.to_bits().to_le_bytes());
            identity.update(&sample.hardness_pa.to_bits().to_le_bytes());
            samples.push(sample);
        }
        if samples.len() != count {
            return Err(HardnessFieldError::Invalid(
                "native point count changed during observation",
            ));
        }
        checkpoint()?;
        Ok(Self {
            sources: Arc::new(Sources {
                identity: source_identity,
                mechanical: mechanical.clone(),
                hardness: hardness.clone(),
                maximum_alpha,
            }),
            samples,
            identity: identity.finalize(),
        })
    }

    /// Original mechanical law/family/branch order, fixed conditions, exact
    /// indentation test and prepared alpha support, independent of point values.
    #[must_use]
    pub fn source_identity(&self) -> ContentHash {
        self.sources.identity
    }

    /// Source identity plus every ordered numeric observation. This does not
    /// attest a mesh, time, displacement field or unobserved loading trajectory.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Complete original fixed-condition native mechanical property bundle.
    #[must_use]
    pub fn mechanical_source(&self) -> &ResolvedMaterialStatePoint {
        &self.sources.mechanical
    }

    /// Original prepared alpha=0 hardness and test receipts. Local accepted H
    /// and alpha are in `samples`; this reference value is not their replacement.
    #[must_use]
    pub fn hardness_source(&self) -> &IndentationHardnessStatePoint {
        &self.sources.hardness
    }

    /// Original continuous support ceiling used to admit every point.
    #[must_use]
    pub fn maximum_accumulated_plastic_strain(&self) -> f64 {
        self.sources.maximum_alpha
    }

    /// Number of retained native points in their original declared order.
    #[must_use]
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// Admitted fields are nonempty; retained for conventional collection use.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }

    /// Every original point without interpolating or averaging its response.
    #[must_use]
    pub fn samples(&self) -> &[AcceptedHardnessSample] {
        &self.samples
    }

    /// Checked lookup for an explicit consumer cell-to-point mapping.
    #[must_use]
    pub fn sample(&self, point: usize) -> Option<&AcceptedHardnessSample> {
        self.samples.get(point)
    }

    /// Local original H [Pa], or None for an unmapped/out-of-range point.
    #[must_use]
    pub fn hardness_pa(&self, point: usize) -> Option<f64> {
        self.sample(point).map(AcceptedHardnessSample::hardness_pa)
    }
}
