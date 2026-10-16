//! As-built FIELD extraction: signed deviation, thickness/gap, warpage form,
//! and profile roughness from registered scan samples.
//!
//! Pose is the first-order as-built correction; FIELDS are the second one.
//! [`crate::rigid3`] answers "where is the part?"; this module answers "what
//! shape is it, and how well do we know that?" — the bond-line thickness that
//! actually got dispensed, the warpage that changed the contact pattern, the
//! gap under a cover, the roughness that governs contact conductance.
//!
//! # Correspondence is SUPPLIED, never discovered
//!
//! Every entry point consumes [`SurfaceSample`]s: a nominal surface point, the
//! outward unit normal declared there, and the measured point asserted to
//! correspond to it. This is the same discipline [`crate::rigid3`] applies to
//! fiducials, and for the same reason — correspondence-free matching (ICP,
//! closest-point projection onto a chart) is an optimization with its own
//! failure modes, and hiding it inside a field extractor would let a bad match
//! masquerade as a measured deviation. Chart-driven correspondence search is
//! an explicit follow-on, not a silent default. It is why this module needs no
//! geometry-query dependency at all.
//!
//! # What the uncertainty is composed from
//!
//! Each admitted point's deviation carries three declared, independently
//! reported contributions:
//!
//! - SCAN RANGE noise, the along-probe component, projected onto the normal.
//!   A range error `e` along the probe direction `v` displaces the point by
//!   `e v`, which moves the signed normal deviation by `e (n . v)`. It
//!   therefore SHRINKS at oblique incidence rather than growing.
//! - REGISTRATION pose uncertainty, through the calibrated 6-dof covariance,
//!   using [`crate::rigid3::CalibratedRigid3Registration::normal_deviation_sensitivity`]
//!   so the pivot convention cannot drift.
//! - CORRESPONDENCE AMBIGUITY, the lateral term. The measurement footprint
//!   elongates as `1 / (n . v)` at oblique incidence, and a lateral offset `d`
//!   projects onto the normal with worst-case factor `sqrt(1 - (n.v)^2)`, so
//!   the contribution is `lateral_sigma * sqrt(1 - c^2) / c` with `c = |n . v|`.
//!   This diverges as the surface turns edge-on to the probe, which is exactly
//!   the near-vertical-wall failure that must be flagged rather than silently
//!   reported. [`ProbeModel::grazing_floor`] refuses those samples; they are
//!   COUNTED in [`DeviationField::grazing`], never dropped.
//!
//! The three are combined in quadrature. That composition declares the field
//! points' measurement errors independent of the fiducial errors that fixed
//! the pose — which is false when one instrument measured both. See the
//! crate CONTRACT's no-claim boundaries.
//!
//! # Determinism and cancellation
//!
//! Every scan is a fixed-order loop over supplied samples with no
//! scheduling-dependent reduction, polling at
//! [`crate::AS_BUILT_POLL_STRIDE_POINTS`] and publishing no partial result on
//! cancellation. Replay on one platform is bitwise.

use crate::AS_BUILT_POLL_STRIDE_POINTS;
use crate::canonical_zero;
use crate::propagate::CoveragePolicy;
use crate::rigid3::{
    CalibratedRigid3Registration, Point3, add3, admit_unit_normal, dot3, mat3_vec, sub3,
};
use fs_evidence::uncertainty::{
    EngineeringUncertaintyKind, EngineeringUncertaintyTerm, TermValue, UncertaintyArtifactRef,
};

/// Schema version of the retained field records.
pub const FIELD_SCHEMA_VERSION: u32 = 1;

/// Maximum samples accepted by one field extraction.
pub const MAX_FIELD_SAMPLES: usize = crate::MAX_AS_BUILT_POINTS;

/// Maximum admitted polynomial total degree for a warpage form fit.
///
/// Warpage is a LOW-ORDER shape statement. Raising the order until the
/// residual vanishes would fit the roughness and call it form.
pub const MAX_FORM_ORDER: u8 = 3;

const FIELD_IDENTITY_DOMAIN: &str = "org.frankensim.fs-asbuilt.field.v1";
const FIELD_IDENTITY_KIND: &[u8] = b"fs-asbuilt-field-v1";
const THICKNESS_IDENTITY_KIND: &[u8] = b"fs-asbuilt-thickness-v1";

/// A structured as-built field-extraction failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FieldError {
    /// No samples were supplied.
    NoSamples,
    /// More samples than one extraction accepts.
    TooManySamples {
        /// Supplied sample count.
        have: usize,
        /// Admitted maximum.
        max: usize,
    },
    /// Two collections that must be index-paired had different lengths.
    LengthMismatch {
        /// Which collection disagreed.
        field: &'static str,
        /// Expected length.
        expected: usize,
        /// Supplied length.
        found: usize,
    },
    /// A supplied scalar failed its admission grammar.
    InvalidScalar {
        /// Offending field.
        field: &'static str,
        /// Violated requirement.
        requirement: &'static str,
    },
    /// The registration refused a pose sensitivity for a sample.
    Registration {
        /// Sample ordinal in supplied order.
        index: usize,
    },
    /// Every supplied sample was refused, so no field exists.
    NoAdmittedSamples {
        /// Count refused for grazing incidence.
        grazing: usize,
    },
    /// A form fit had fewer admitted points than its basis has coefficients.
    UnderdeterminedFit {
        /// Admitted point count.
        points: usize,
        /// Basis size for the declared order.
        coefficients: usize,
    },
    /// The form-fit design matrix was rank-deficient at the declared order.
    ///
    /// Collinear or coincident in-plane coordinates cannot determine a tilt or
    /// curvature term; the fit refuses rather than returning an arbitrary
    /// member of the solution family.
    RankDeficientFit {
        /// Declared total degree.
        order: u8,
    },
    /// A declared polynomial order exceeded [`MAX_FORM_ORDER`].
    FormOrderTooHigh {
        /// Declared order.
        order: u8,
        /// Admitted maximum.
        max: u8,
    },
    /// A roughness profile had too few points for its declared segmentation.
    ProfileTooShort {
        /// Supplied height count.
        points: usize,
        /// Minimum required.
        need: usize,
    },
    /// An arithmetic step left the finite range.
    ArithmeticOverflow {
        /// Offending quantity.
        field: &'static str,
    },
    /// A capacity reservation failed.
    AllocationFailed,
    /// Cancellation was observed; no partial field is published.
    Cancelled {
        /// Phase at which cancellation was observed.
        phase: &'static str,
    },
    /// An uncertainty term could not be admitted by `fs-evidence`.
    TermRefused,
}

impl core::fmt::Display for FieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoSamples => write!(f, "no as-built surface samples were supplied"),
            Self::TooManySamples { have, max } => {
                write!(f, "{have} samples exceeds the admitted maximum {max}")
            }
            Self::LengthMismatch {
                field,
                expected,
                found,
            } => write!(f, "{field} has {found} entries, expected {expected}"),
            Self::InvalidScalar { field, requirement } => {
                write!(f, "{field} must be {requirement}")
            }
            Self::Registration { index } => {
                write!(
                    f,
                    "registration refused a pose sensitivity at sample {index}"
                )
            }
            Self::NoAdmittedSamples { grazing } => write!(
                f,
                "every sample was refused ({grazing} for grazing incidence); no field exists"
            ),
            Self::UnderdeterminedFit {
                points,
                coefficients,
            } => write!(
                f,
                "a {coefficients}-coefficient form fit needs at least that many points, got {points}"
            ),
            Self::RankDeficientFit { order } => write!(
                f,
                "the in-plane coordinates cannot determine an order-{order} form fit"
            ),
            Self::FormOrderTooHigh { order, max } => {
                write!(f, "form order {order} exceeds the admitted maximum {max}")
            }
            Self::ProfileTooShort { points, need } => {
                write!(f, "a profile of {points} points needs at least {need}")
            }
            Self::ArithmeticOverflow { field } => write!(f, "{field} left the finite range"),
            Self::AllocationFailed => write!(f, "a field allocation failed"),
            Self::Cancelled { phase } => write!(f, "cancelled during {phase}"),
            Self::TermRefused => write!(f, "the composed uncertainty term was refused"),
        }
    }
}

impl std::error::Error for FieldError {}

/// One supplied surface correspondence.
///
/// `design` and `normal` are DESIGN-frame; `measured` is MEASURED-frame. The
/// normal is the outward unit normal of the nominal surface at `design`, where
/// "outward" means away from the solid — the sign convention that makes a
/// positive deviation mean material displaced outward.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceSample {
    design: Point3,
    normal: [f64; 3],
    measured: Point3,
}

impl SurfaceSample {
    /// Declare one correspondence.
    ///
    /// # Errors
    /// A non-finite or non-unit normal.
    pub fn new(
        design: Point3,
        normal: [f64; 3],
        measured: Point3,
    ) -> Result<Self, crate::rigid3::Rigid3Error> {
        admit_unit_normal(normal)?;
        Ok(Self {
            design,
            normal,
            measured,
        })
    }

    /// Nominal design-frame surface point.
    #[must_use]
    pub const fn design(&self) -> Point3 {
        self.design
    }

    /// Declared design-frame outward unit normal.
    #[must_use]
    pub const fn normal(&self) -> [f64; 3] {
        self.normal
    }

    /// Measured-frame scan point.
    #[must_use]
    pub const fn measured(&self) -> Point3 {
        self.measured
    }
}

/// The scan instrument's declared first-order error model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProbeModel {
    direction: [f64; 3],
    grazing_floor: f64,
    lateral_sigma: f64,
    range_sigma: f64,
}

impl ProbeModel {
    /// Declare the probe model.
    ///
    /// `direction` is the MEASURED-frame unit probe/line-of-sight direction.
    /// `grazing_floor` is the smallest admitted `|n . v|`; samples below it
    /// are flagged, not measured. `lateral_sigma` and `range_sigma` are
    /// one-sigma lateral and along-probe point-position uncertainties in the
    /// deviation's length unit.
    ///
    /// # Errors
    /// A non-unit direction, a floor outside `(0, 1]`, or a negative or
    /// non-finite sigma.
    pub fn new(
        direction: [f64; 3],
        grazing_floor: f64,
        lateral_sigma: f64,
        range_sigma: f64,
    ) -> Result<Self, FieldError> {
        admit_unit_normal(direction).map_err(|_| FieldError::InvalidScalar {
            field: "probe.direction",
            requirement: "finite and unit length",
        })?;
        if !grazing_floor.is_finite() || grazing_floor <= 0.0 || grazing_floor > 1.0 {
            return Err(FieldError::InvalidScalar {
                field: "probe.grazing_floor",
                requirement: "inside (0, 1]",
            });
        }
        for (value, field) in [
            (lateral_sigma, "probe.lateral_sigma"),
            (range_sigma, "probe.range_sigma"),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(FieldError::InvalidScalar {
                    field,
                    requirement: "finite and non-negative",
                });
            }
        }
        Ok(Self {
            direction,
            grazing_floor,
            lateral_sigma,
            range_sigma,
        })
    }

    /// Measured-frame unit probe direction.
    #[must_use]
    pub const fn direction(&self) -> [f64; 3] {
        self.direction
    }

    /// Smallest admitted incidence cosine.
    #[must_use]
    pub const fn grazing_floor(&self) -> f64 {
        self.grazing_floor
    }

    /// Declared one-sigma lateral position uncertainty.
    #[must_use]
    pub const fn lateral_sigma(&self) -> f64 {
        self.lateral_sigma
    }

    /// Declared one-sigma along-probe range uncertainty.
    #[must_use]
    pub const fn range_sigma(&self) -> f64 {
        self.range_sigma
    }
}

/// Why a supplied sample did not become a field point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SampleDisposition {
    /// Incidence cosine fell at or below the declared grazing floor.
    Grazing,
}

/// One admitted field point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FieldPoint {
    deviation: f64,
    incidence_cosine: f64,
    scan_sigma: f64,
    registration_sigma: f64,
    ambiguity_sigma: f64,
    half_width: f64,
}

impl FieldPoint {
    /// Signed deviation along the registered outward normal. Positive means
    /// the measured surface sits outboard of nominal.
    #[must_use]
    pub const fn deviation(&self) -> f64 {
        self.deviation
    }

    /// `|n . v|` at this point.
    #[must_use]
    pub const fn incidence_cosine(&self) -> f64 {
        self.incidence_cosine
    }

    /// Range-noise contribution, one sigma.
    #[must_use]
    pub const fn scan_sigma(&self) -> f64 {
        self.scan_sigma
    }

    /// Registration-pose contribution, one sigma.
    #[must_use]
    pub const fn registration_sigma(&self) -> f64 {
        self.registration_sigma
    }

    /// Correspondence-ambiguity contribution, one sigma.
    #[must_use]
    pub const fn ambiguity_sigma(&self) -> f64 {
        self.ambiguity_sigma
    }

    /// Coverage-scaled half-width of the composed deviation uncertainty.
    #[must_use]
    pub const fn half_width(&self) -> f64 {
        self.half_width
    }

    /// Whether the deviation exceeds its own coverage-scaled half-width.
    ///
    /// This is a resolution screen, not a tolerance verdict: it says the
    /// measurement can tell this deviation from zero under the declared model,
    /// not that the part conforms.
    #[must_use]
    pub fn resolved(&self) -> bool {
        self.deviation.abs() > self.half_width
    }
}

/// An extracted signed-deviation field over supplied correspondences.
#[derive(Debug, Clone, PartialEq)]
pub struct DeviationField {
    points: Vec<FieldPoint>,
    admitted_indices: Vec<usize>,
    grazing: Vec<usize>,
    registration_ref: fs_blake3::ContentHash,
    coverage: CoveragePolicy,
    identity: fs_blake3::ContentHash,
}

impl DeviationField {
    /// Extract signed normal deviations with composed per-point uncertainty.
    ///
    /// Samples whose incidence cosine falls at or below
    /// [`ProbeModel::grazing_floor`] are refused and counted; the remainder
    /// become field points in supplied order.
    ///
    /// # Errors
    /// An empty or oversized sample set, a registration that refuses a
    /// sensitivity, an all-refused sample set, arithmetic overflow, or
    /// cancellation.
    pub fn extract(
        samples: &[SurfaceSample],
        registration: &CalibratedRigid3Registration,
        probe: &ProbeModel,
        coverage: CoveragePolicy,
        cx: &fs_exec::Cx<'_>,
    ) -> Result<Self, FieldError> {
        admit_sample_count(samples.len())?;
        let rotation = registration.registration().rotation();
        let translation = registration.registration().translation();
        let covariance = registration.covariance();

        let mut points = Vec::new();
        let mut admitted_indices = Vec::new();
        let mut grazing = Vec::new();
        points
            .try_reserve_exact(samples.len())
            .map_err(|_| FieldError::AllocationFailed)?;
        admitted_indices
            .try_reserve_exact(samples.len())
            .map_err(|_| FieldError::AllocationFailed)?;

        for (index, sample) in samples.iter().enumerate() {
            checkpoint(cx, index, "field.deviation")?;
            let rotated_normal = mat3_vec(rotation, sample.normal);
            let incidence = dot3(rotated_normal, probe.direction).abs();
            if !incidence.is_finite() {
                return Err(FieldError::ArithmeticOverflow {
                    field: "incidence cosine",
                });
            }
            if incidence <= probe.grazing_floor {
                grazing
                    .try_reserve(1)
                    .map_err(|_| FieldError::AllocationFailed)?;
                grazing.push(index);
                continue;
            }

            let registered = add3(mat3_vec(rotation, sample.design.coords()), translation);
            let offset = sub3(sample.measured.coords(), registered);
            let deviation = dot3(offset, rotated_normal);

            let gradient = registration
                .normal_deviation_sensitivity(sample.design, sample.normal)
                .map_err(|_| FieldError::Registration { index })?;
            let registration_sigma = quadratic_form_sigma(covariance, &gradient)?;

            let scan_sigma = probe.range_sigma * incidence;
            let tangential = (1.0 - incidence * incidence).max(0.0).sqrt();
            let ambiguity_sigma = probe.lateral_sigma * tangential / incidence;

            let combined = quadrature3(scan_sigma, registration_sigma, ambiguity_sigma)?;
            let half_width = coverage.factor() * combined;
            for (value, field) in [
                (deviation, "deviation"),
                (registration_sigma, "registration sigma"),
                (ambiguity_sigma, "ambiguity sigma"),
                (half_width, "half width"),
            ] {
                if !value.is_finite() {
                    return Err(FieldError::ArithmeticOverflow { field });
                }
            }

            points.push(FieldPoint {
                deviation: canonical_zero(deviation),
                incidence_cosine: canonical_zero(incidence),
                scan_sigma: canonical_zero(scan_sigma),
                registration_sigma: canonical_zero(registration_sigma),
                ambiguity_sigma: canonical_zero(ambiguity_sigma),
                half_width: canonical_zero(half_width),
            });
            admitted_indices.push(index);
        }

        if points.is_empty() {
            return Err(FieldError::NoAdmittedSamples {
                grazing: grazing.len(),
            });
        }

        cx.checkpoint().map_err(|_| FieldError::Cancelled {
            phase: "field.deviation-publish",
        })?;

        let registration_ref = registration.model_identity();
        let identity = field_identity(
            FIELD_IDENTITY_KIND,
            registration_ref,
            &coverage,
            probe,
            &points,
            &grazing,
        );
        Ok(Self {
            points,
            admitted_indices,
            grazing,
            registration_ref,
            coverage,
            identity,
        })
    }

    /// Admitted field points, in supplied order.
    #[must_use]
    pub fn points(&self) -> &[FieldPoint] {
        &self.points
    }

    /// Supplied-order index of each admitted point.
    #[must_use]
    pub fn admitted_indices(&self) -> &[usize] {
        &self.admitted_indices
    }

    /// Supplied-order indices refused for grazing incidence, counted not
    /// dropped.
    #[must_use]
    pub fn grazing(&self) -> &[usize] {
        &self.grazing
    }

    /// Disposition of one refused sample.
    #[must_use]
    pub fn disposition(&self, index: usize) -> Option<SampleDisposition> {
        self.grazing
            .binary_search(&index)
            .ok()
            .map(|_| SampleDisposition::Grazing)
    }

    /// Total supplied samples, admitted plus refused.
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.points.len() + self.grazing.len()
    }

    /// Cited calibrated-registration identity.
    #[must_use]
    pub const fn registration_ref(&self) -> fs_blake3::ContentHash {
        self.registration_ref
    }

    /// Declared coverage policy.
    #[must_use]
    pub const fn coverage(&self) -> &CoveragePolicy {
        &self.coverage
    }

    /// Domain-separated identity of this field record. An integrity address,
    /// not authentication.
    #[must_use]
    pub const fn identity(&self) -> fs_blake3::ContentHash {
        self.identity
    }

    /// Peak-to-valley span of the admitted deviations.
    #[must_use]
    pub fn peak_to_valley(&self) -> f64 {
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        for point in &self.points {
            low = low.min(point.deviation);
            high = high.max(point.deviation);
        }
        high - low
    }

    /// Largest coverage-scaled half-width over the admitted points.
    #[must_use]
    pub fn max_half_width(&self) -> f64 {
        self.points
            .iter()
            .fold(0.0f64, |acc, point| acc.max(point.half_width))
    }

    /// Smallest coverage-scaled half-width over the admitted points.
    #[must_use]
    pub fn min_half_width(&self) -> f64 {
        self.points
            .iter()
            .fold(f64::INFINITY, |acc, point| acc.min(point.half_width))
    }

    /// Project the field's measurement uncertainty into one eight-term budget
    /// slot.
    ///
    /// The bracket is `[min, max]` of the per-point coverage-scaled
    /// half-widths, which is what [`TermValue::IntervalBound`] means: it
    /// brackets the HALF-WIDTH, and aggregation reads the upper endpoint. The
    /// term is [`EngineeringUncertaintyKind::Measurement`], not `Geometry`:
    /// the pose contribution reaches a budget through
    /// [`crate::propagate`], and double-counting it here would inflate the
    /// budget while looking more rigorous.
    ///
    /// # Errors
    /// A role identifier or bracket that `fs-evidence` refuses.
    pub fn measurement_term(&self, role: &str) -> Result<EngineeringUncertaintyTerm, FieldError> {
        let artifact = UncertaintyArtifactRef::new(role, self.identity)
            .map_err(|_| FieldError::TermRefused)?;
        let value = TermValue::interval(self.min_half_width(), self.max_half_width())
            .map_err(|_| FieldError::TermRefused)?;
        EngineeringUncertaintyTerm::try_new(
            EngineeringUncertaintyKind::Measurement,
            value,
            artifact,
        )
        .map_err(|_| FieldError::TermRefused)
    }
}

/// One paired local-thickness point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThicknessPoint {
    nominal: f64,
    thickness: f64,
    registration_sigma: f64,
    measurement_sigma: f64,
    half_width: f64,
}

impl ThicknessPoint {
    /// Declared nominal thickness at this station.
    #[must_use]
    pub const fn nominal(&self) -> f64 {
        self.nominal
    }

    /// As-built local thickness.
    #[must_use]
    pub const fn thickness(&self) -> f64 {
        self.thickness
    }

    /// Signed departure from nominal.
    #[must_use]
    pub fn change(&self) -> f64 {
        self.thickness - self.nominal
    }

    /// Registration contribution to the thickness uncertainty, one sigma.
    ///
    /// Computed from the SUMMED pose sensitivity of both faces, not by adding
    /// the two faces' sigmas in quadrature. Both deviations ride the same
    /// pose, so their pose errors are perfectly correlated: on opposing faces
    /// the translation sensitivities are near-opposite and cancel. Treating
    /// them as independent would manufacture an uncertainty that the rigid
    /// motion cannot produce.
    #[must_use]
    pub const fn registration_sigma(&self) -> f64 {
        self.registration_sigma
    }

    /// Combined per-face measurement contribution, one sigma.
    #[must_use]
    pub const fn measurement_sigma(&self) -> f64 {
        self.measurement_sigma
    }

    /// Coverage-scaled half-width of the thickness uncertainty.
    #[must_use]
    pub const fn half_width(&self) -> f64 {
        self.half_width
    }

    /// Whether the departure from nominal exceeds its own half-width.
    #[must_use]
    pub fn resolved(&self) -> bool {
        self.change().abs() > self.half_width
    }
}

/// An extracted local-thickness (or gap) field over paired faces.
#[derive(Debug, Clone, PartialEq)]
pub struct ThicknessField {
    points: Vec<ThicknessPoint>,
    grazing: Vec<usize>,
    registration_ref: fs_blake3::ContentHash,
    coverage: CoveragePolicy,
    identity: fs_blake3::ContentHash,
}

impl ThicknessField {
    /// Extract local thickness from index-paired opposing-face samples.
    ///
    /// `side_a[i]` and `side_b[i]` are the two faces at station `i`, each with
    /// its own OUTWARD normal (pointing away from the material between them),
    /// and `nominal[i]` is the design thickness there. As-built thickness is
    /// `nominal + d_a + d_b`: an outward deviation on either face adds
    /// material.
    ///
    /// For a gap between two parts, supply the two facing surfaces with their
    /// outward normals and the nominal gap; a positive deviation on either
    /// face then means material intruding into the gap, so the same sum
    /// reports the as-built clearance with the sign already handled by the
    /// declared normals.
    ///
    /// A station is refused if EITHER face is grazing, and its index is
    /// counted in [`ThicknessField::grazing`].
    ///
    /// # Errors
    /// Mismatched lengths, an empty or oversized set, a non-finite nominal, a
    /// refusing registration, an all-refused set, or cancellation.
    #[allow(clippy::too_many_lines)] // One linear paired-station scan.
    pub fn extract(
        side_a: &[SurfaceSample],
        side_b: &[SurfaceSample],
        nominal: &[f64],
        registration: &CalibratedRigid3Registration,
        probe: &ProbeModel,
        coverage: CoveragePolicy,
        cx: &fs_exec::Cx<'_>,
    ) -> Result<Self, FieldError> {
        admit_sample_count(side_a.len())?;
        if side_b.len() != side_a.len() {
            return Err(FieldError::LengthMismatch {
                field: "side_b",
                expected: side_a.len(),
                found: side_b.len(),
            });
        }
        if nominal.len() != side_a.len() {
            return Err(FieldError::LengthMismatch {
                field: "nominal",
                expected: side_a.len(),
                found: nominal.len(),
            });
        }

        let rotation = registration.registration().rotation();
        let translation = registration.registration().translation();
        let covariance = registration.covariance();

        let mut points = Vec::new();
        let mut grazing = Vec::new();
        points
            .try_reserve_exact(side_a.len())
            .map_err(|_| FieldError::AllocationFailed)?;

        for index in 0..side_a.len() {
            checkpoint(cx, index, "field.thickness")?;
            let station_nominal = nominal[index];
            if !station_nominal.is_finite() {
                return Err(FieldError::InvalidScalar {
                    field: "nominal thickness",
                    requirement: "finite",
                });
            }

            let mut deviations = [0.0f64; 2];
            let mut gradients = [[0.0f64; 6]; 2];
            let mut per_face_sigma = [0.0f64; 2];
            let mut refused = false;
            for (slot, sample) in [side_a[index], side_b[index]].into_iter().enumerate() {
                let rotated_normal = mat3_vec(rotation, sample.normal);
                let incidence = dot3(rotated_normal, probe.direction).abs();
                if !incidence.is_finite() {
                    return Err(FieldError::ArithmeticOverflow {
                        field: "incidence cosine",
                    });
                }
                if incidence <= probe.grazing_floor {
                    refused = true;
                    break;
                }
                let registered = add3(mat3_vec(rotation, sample.design.coords()), translation);
                let offset = sub3(sample.measured.coords(), registered);
                deviations[slot] = dot3(offset, rotated_normal);
                gradients[slot] = registration
                    .normal_deviation_sensitivity(sample.design, sample.normal)
                    .map_err(|_| FieldError::Registration { index })?;
                let tangential = (1.0 - incidence * incidence).max(0.0).sqrt();
                per_face_sigma[slot] = quadrature2(
                    probe.range_sigma * incidence,
                    probe.lateral_sigma * tangential / incidence,
                )?;
            }
            if refused {
                grazing
                    .try_reserve(1)
                    .map_err(|_| FieldError::AllocationFailed)?;
                grazing.push(index);
                continue;
            }

            // The two faces ride ONE pose: sum the sensitivities, then form
            // the quadratic form once. Quadrature over per-face sigmas would
            // discard the cancellation that makes thickness robust to rigid
            // motion.
            let mut summed = [0.0f64; 6];
            for axis in 0..6 {
                summed[axis] = gradients[0][axis] + gradients[1][axis];
            }
            let registration_sigma = quadratic_form_sigma(covariance, &summed)?;
            let measurement_sigma = quadrature2(per_face_sigma[0], per_face_sigma[1])?;
            let combined = quadrature2(registration_sigma, measurement_sigma)?;
            let half_width = coverage.factor() * combined;
            let thickness = station_nominal + deviations[0] + deviations[1];
            for (value, field) in [
                (thickness, "thickness"),
                (registration_sigma, "thickness registration sigma"),
                (half_width, "thickness half width"),
            ] {
                if !value.is_finite() {
                    return Err(FieldError::ArithmeticOverflow { field });
                }
            }

            points.push(ThicknessPoint {
                nominal: canonical_zero(station_nominal),
                thickness: canonical_zero(thickness),
                registration_sigma: canonical_zero(registration_sigma),
                measurement_sigma: canonical_zero(measurement_sigma),
                half_width: canonical_zero(half_width),
            });
        }

        if points.is_empty() {
            return Err(FieldError::NoAdmittedSamples {
                grazing: grazing.len(),
            });
        }

        cx.checkpoint().map_err(|_| FieldError::Cancelled {
            phase: "field.thickness-publish",
        })?;

        let registration_ref = registration.model_identity();
        let identity = thickness_identity(registration_ref, &coverage, probe, &points, &grazing);
        Ok(Self {
            points,
            grazing,
            registration_ref,
            coverage,
            identity,
        })
    }

    /// Admitted thickness stations, in supplied order.
    #[must_use]
    pub fn points(&self) -> &[ThicknessPoint] {
        &self.points
    }

    /// Station indices refused because a face was grazing.
    #[must_use]
    pub fn grazing(&self) -> &[usize] {
        &self.grazing
    }

    /// Cited calibrated-registration identity.
    #[must_use]
    pub const fn registration_ref(&self) -> fs_blake3::ContentHash {
        self.registration_ref
    }

    /// Declared coverage policy.
    #[must_use]
    pub const fn coverage(&self) -> &CoveragePolicy {
        &self.coverage
    }

    /// Domain-separated identity of this thickness record.
    #[must_use]
    pub const fn identity(&self) -> fs_blake3::ContentHash {
        self.identity
    }

    /// Smallest admitted as-built thickness.
    #[must_use]
    pub fn min_thickness(&self) -> f64 {
        self.points
            .iter()
            .fold(f64::INFINITY, |acc, point| acc.min(point.thickness))
    }

    /// Largest admitted as-built thickness.
    #[must_use]
    pub fn max_thickness(&self) -> f64 {
        self.points
            .iter()
            .fold(f64::NEG_INFINITY, |acc, point| acc.max(point.thickness))
    }
}

/// A low-order polynomial form (warpage) fit of a deviation field.
#[derive(Debug, Clone, PartialEq)]
pub struct FormFit {
    order: u8,
    coefficients: Vec<f64>,
    center: [f64; 2],
    half_extent: [f64; 2],
    rms_residual: f64,
    max_abs_residual: f64,
    peak_to_valley: f64,
    point_count: usize,
}

impl FormFit {
    /// Declared total polynomial degree.
    #[must_use]
    pub const fn order(&self) -> u8 {
        self.order
    }

    /// Fitted coefficients in the NORMALIZED coordinate frame, in the
    /// monomial order produced by [`FormFit::basis_labels`].
    #[must_use]
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    /// Monomial labels matching [`FormFit::coefficients`], in the same order.
    #[must_use]
    pub fn basis_labels(&self) -> Vec<String> {
        basis_exponents(self.order)
            .into_iter()
            .map(|(i, j)| match (i, j) {
                (0, 0) => "1".to_string(),
                (i, 0) => format!("u^{i}"),
                (0, j) => format!("v^{j}"),
                (i, j) => format!("u^{i} v^{j}"),
            })
            .collect()
    }

    /// Normalization center of the in-plane coordinates.
    #[must_use]
    pub const fn center(&self) -> [f64; 2] {
        self.center
    }

    /// Normalization half-extent of the in-plane coordinates.
    #[must_use]
    pub const fn half_extent(&self) -> [f64; 2] {
        self.half_extent
    }

    /// Root-mean-square fit residual.
    ///
    /// This is what the declared order did NOT capture. A warpage number
    /// quoted without it is a claim about the fit, not about the part.
    #[must_use]
    pub const fn rms_residual(&self) -> f64 {
        self.rms_residual
    }

    /// Largest absolute fit residual.
    #[must_use]
    pub const fn max_abs_residual(&self) -> f64 {
        self.max_abs_residual
    }

    /// Peak-to-valley span of the FITTED surface over the sample stations.
    ///
    /// This is the warpage statement: the low-order shape's span, evaluated
    /// where there is data. It is not the field's own peak-to-valley, which
    /// includes everything the fit rejected.
    #[must_use]
    pub const fn peak_to_valley(&self) -> f64 {
        self.peak_to_valley
    }

    /// Points that entered the fit.
    #[must_use]
    pub const fn point_count(&self) -> usize {
        self.point_count
    }

    /// Evaluate the fitted form at raw in-plane coordinates.
    #[must_use]
    pub fn evaluate(&self, coordinate: [f64; 2]) -> f64 {
        let u = normalize_axis(coordinate[0], self.center[0], self.half_extent[0]);
        let v = normalize_axis(coordinate[1], self.center[1], self.half_extent[1]);
        let mut value = 0.0f64;
        for (coefficient, (i, j)) in self.coefficients.iter().zip(basis_exponents(self.order)) {
            value = coefficient.mul_add(powi_exact(u, i) * powi_exact(v, j), value);
        }
        value
    }
}

/// Fit a low-order polynomial form to a deviation field.
///
/// `coordinates` are the in-plane stations of the ADMITTED field points, in
/// [`DeviationField::points`] order — one per admitted point, not one per
/// supplied sample. Coordinates are normalized to `[-1, 1]` internally so the
/// least-squares problem stays well conditioned; the returned coefficients are
/// in that normalized frame and [`FormFit::evaluate`] applies the mapping.
///
/// # Errors
/// A length mismatch, an order above [`MAX_FORM_ORDER`], fewer points than
/// coefficients, degenerate or non-finite coordinates, or a rank-deficient
/// design matrix.
#[allow(clippy::too_many_lines)] // One linear admit-normalize-solve-residual pass.
pub fn fit_form(
    field: &DeviationField,
    coordinates: &[[f64; 2]],
    order: u8,
) -> Result<FormFit, FieldError> {
    if order > MAX_FORM_ORDER {
        return Err(FieldError::FormOrderTooHigh {
            order,
            max: MAX_FORM_ORDER,
        });
    }
    let points = field.points();
    if coordinates.len() != points.len() {
        return Err(FieldError::LengthMismatch {
            field: "coordinates",
            expected: points.len(),
            found: coordinates.len(),
        });
    }
    let exponents = basis_exponents(order);
    let width = exponents.len();
    if points.len() < width {
        return Err(FieldError::UnderdeterminedFit {
            points: points.len(),
            coefficients: width,
        });
    }

    let mut low = [f64::INFINITY; 2];
    let mut high = [f64::NEG_INFINITY; 2];
    for coordinate in coordinates {
        for axis in 0..2 {
            if !coordinate[axis].is_finite() {
                return Err(FieldError::InvalidScalar {
                    field: "in-plane coordinate",
                    requirement: "finite",
                });
            }
            low[axis] = low[axis].min(coordinate[axis]);
            high[axis] = high[axis].max(coordinate[axis]);
        }
    }
    let mut center = [0.0f64; 2];
    let mut half_extent = [0.0f64; 2];
    for axis in 0..2 {
        center[axis] = f64::midpoint(low[axis], high[axis]);
        let extent = 0.5 * (high[axis] - low[axis]);
        // A degenerate axis (all stations on one line) is legal for orders
        // that do not use it; the rank check below is what refuses a fit that
        // actually needs the missing direction.
        half_extent[axis] = if extent > 0.0 { extent } else { 1.0 };
    }

    let rows = points.len();
    let mut design = Vec::new();
    design
        .try_reserve_exact(rows * width)
        .map_err(|_| FieldError::AllocationFailed)?;
    for coordinate in coordinates {
        let u = normalize_axis(coordinate[0], center[0], half_extent[0]);
        let v = normalize_axis(coordinate[1], center[1], half_extent[1]);
        for &(i, j) in &exponents {
            design.push(powi_exact(u, i) * powi_exact(v, j));
        }
    }
    let mut rhs = Vec::new();
    rhs.try_reserve_exact(rows)
        .map_err(|_| FieldError::AllocationFailed)?;
    for point in points {
        rhs.push(point.deviation());
    }

    let factorization = fs_la::factor::qr(&design, rows, width);
    // A rank-deficient design matrix makes solve_ls divide by a zero pivot and
    // return non-finite coefficients. Refuse on the factor itself rather than
    // discovering it downstream in an evaluated surface.
    let mut smallest = f64::INFINITY;
    let mut largest = 0.0f64;
    for k in 0..width {
        let pivot = factorization.r(k, k).abs();
        smallest = smallest.min(pivot);
        largest = largest.max(pivot);
    }
    if !(smallest.is_finite() && largest > 0.0) || smallest <= largest * 1e-12 {
        return Err(FieldError::RankDeficientFit { order });
    }
    let coefficients = factorization.solve_ls(&rhs);
    if coefficients.len() != width || !coefficients.iter().all(|value| value.is_finite()) {
        return Err(FieldError::RankDeficientFit { order });
    }

    let mut sum_squares = 0.0f64;
    let mut max_abs_residual = 0.0f64;
    let mut fitted_low = f64::INFINITY;
    let mut fitted_high = f64::NEG_INFINITY;
    for (row, point) in points.iter().enumerate() {
        let mut fitted = 0.0f64;
        for column in 0..width {
            fitted = coefficients[column].mul_add(design[row * width + column], fitted);
        }
        let residual = point.deviation() - fitted;
        sum_squares = residual.mul_add(residual, sum_squares);
        max_abs_residual = max_abs_residual.max(residual.abs());
        fitted_low = fitted_low.min(fitted);
        fitted_high = fitted_high.max(fitted);
    }
    let rms_residual = (sum_squares / rows as f64).sqrt();
    let peak_to_valley = fitted_high - fitted_low;
    for (value, field) in [
        (rms_residual, "form rms residual"),
        (peak_to_valley, "form peak-to-valley"),
    ] {
        if !value.is_finite() {
            return Err(FieldError::ArithmeticOverflow { field });
        }
    }

    Ok(FormFit {
        order,
        coefficients: coefficients.into_iter().map(canonical_zero).collect(),
        center,
        half_extent,
        rms_residual: canonical_zero(rms_residual),
        max_abs_residual: canonical_zero(max_abs_residual),
        peak_to_valley: canonical_zero(peak_to_valley),
        point_count: rows,
    })
}

/// Declared form removed from a profile before roughness statistics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProfileForm {
    /// Subtract the arithmetic mean only.
    Mean,
    /// Subtract the least-squares straight mean line.
    Line,
    /// Subtract the least-squares parabolic mean line.
    Parabola,
}

impl ProfileForm {
    const fn coefficients(self) -> usize {
        match self {
            Self::Mean => 1,
            Self::Line => 2,
            Self::Parabola => 3,
        }
    }

    /// Stable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Mean => "mean",
            Self::Line => "line",
            Self::Parabola => "parabola",
        }
    }
}

/// Unfiltered profile-height statistics after declared form removal.
///
/// # These are NOT ISO 4287 parameters
///
/// ISO 4287 `Ra`/`Rq`/`Rz` are defined on a profile separated into roughness
/// and waviness by a phase-correct Gaussian filter at a declared cutoff
/// `lambda_c` (ISO 16610-21), evaluated over a stated number of sampling
/// lengths. NO SUCH FILTER IS APPLIED HERE. These are the same arithmetic
/// applied to a merely form-removed profile, so any waviness inside the
/// supplied trace is counted as roughness and the values will generally read
/// HIGH against an instrument that filters. They are usable as a
/// self-consistent relative statistic and as an input declaration; they are
/// not a conformance value and must not be reported as `Ra` without the
/// unfiltered qualifier.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProfileStatistics {
    form: ProfileForm,
    ra: f64,
    rq: f64,
    rt: f64,
    rz: f64,
    segments: usize,
    point_count: usize,
}

impl ProfileStatistics {
    /// Declared form removal.
    #[must_use]
    pub const fn form(self) -> ProfileForm {
        self.form
    }

    /// Arithmetic mean deviation of the form-removed profile (unfiltered).
    #[must_use]
    pub const fn ra(self) -> f64 {
        self.ra
    }

    /// Root-mean-square deviation of the form-removed profile (unfiltered).
    #[must_use]
    pub const fn rq(self) -> f64 {
        self.rq
    }

    /// Total height: max peak to min valley over the whole trace (unfiltered).
    #[must_use]
    pub const fn rt(self) -> f64 {
        self.rt
    }

    /// Mean per-segment peak-to-valley over the declared segments
    /// (unfiltered).
    #[must_use]
    pub const fn rz(self) -> f64 {
        self.rz
    }

    /// Segments the trace was divided into for `Rz`.
    #[must_use]
    pub const fn segments(self) -> usize {
        self.segments
    }

    /// Profile points consumed.
    #[must_use]
    pub const fn point_count(self) -> usize {
        self.point_count
    }
}

/// Compute unfiltered profile statistics after removing a declared form.
///
/// `heights` are ordered profile heights at uniform spacing. `segments` is the
/// number of equal sampling segments used for `Rz`; the trace is split by
/// index, and a trailing partial segment is merged into the last full one so
/// every point is counted exactly once.
///
/// # Errors
/// A non-finite height, zero segments, or a trace too short for the declared
/// form and segmentation.
pub fn profile_statistics(
    heights: &[f64],
    form: ProfileForm,
    segments: usize,
) -> Result<ProfileStatistics, FieldError> {
    if segments == 0 {
        return Err(FieldError::InvalidScalar {
            field: "segments",
            requirement: "at least one",
        });
    }
    let need = form.coefficients().max(segments);
    if heights.len() < need {
        return Err(FieldError::ProfileTooShort {
            points: heights.len(),
            need,
        });
    }
    if !heights.iter().all(|value| value.is_finite()) {
        return Err(FieldError::InvalidScalar {
            field: "profile height",
            requirement: "finite",
        });
    }

    let count = heights.len();
    let width = form.coefficients();
    // Normalized abscissa keeps the mean-line solve well conditioned for long
    // traces; the fitted line is subtracted, so the frame is immaterial.
    let denominator = if count > 1 { (count - 1) as f64 } else { 1.0 };
    let mut design = Vec::new();
    design
        .try_reserve_exact(count * width)
        .map_err(|_| FieldError::AllocationFailed)?;
    for index in 0..count {
        let x = 2.0 * (index as f64) / denominator - 1.0;
        for power in 0..width {
            design.push(powi_exact(x, power));
        }
    }
    let factorization = fs_la::factor::qr(&design, count, width);
    let coefficients = factorization.solve_ls(heights);
    if coefficients.len() != width || !coefficients.iter().all(|value| value.is_finite()) {
        return Err(FieldError::ArithmeticOverflow {
            field: "profile mean line",
        });
    }

    let mut residuals = Vec::new();
    residuals
        .try_reserve_exact(count)
        .map_err(|_| FieldError::AllocationFailed)?;
    for (index, height) in heights.iter().enumerate() {
        let mut fitted = 0.0f64;
        for column in 0..width {
            fitted = coefficients[column].mul_add(design[index * width + column], fitted);
        }
        residuals.push(height - fitted);
    }

    let mut sum_abs = 0.0f64;
    let mut sum_squares = 0.0f64;
    let mut low = f64::INFINITY;
    let mut high = f64::NEG_INFINITY;
    for residual in &residuals {
        sum_abs += residual.abs();
        sum_squares = residual.mul_add(*residual, sum_squares);
        low = low.min(*residual);
        high = high.max(*residual);
    }
    let points = count as f64;
    let ra = sum_abs / points;
    let rq = (sum_squares / points).sqrt();
    let rt = high - low;

    // Rz: equal index segments, trailing remainder folded into the last one.
    let span = count / segments;
    let mut sum_segment_span = 0.0f64;
    for segment in 0..segments {
        let start = segment * span;
        let end = if segment + 1 == segments {
            count
        } else {
            (segment + 1) * span
        };
        let mut segment_low = f64::INFINITY;
        let mut segment_high = f64::NEG_INFINITY;
        for residual in &residuals[start..end] {
            segment_low = segment_low.min(*residual);
            segment_high = segment_high.max(*residual);
        }
        sum_segment_span += segment_high - segment_low;
    }
    let rz = sum_segment_span / segments as f64;

    for (value, field) in [(ra, "Ra"), (rq, "Rq"), (rt, "Rt"), (rz, "Rz")] {
        if !value.is_finite() {
            return Err(FieldError::ArithmeticOverflow { field });
        }
    }

    Ok(ProfileStatistics {
        form,
        ra: canonical_zero(ra),
        rq: canonical_zero(rq),
        rt: canonical_zero(rt),
        rz: canonical_zero(rz),
        segments,
        point_count: count,
    })
}

// ---------------------------------------------------------------------------
// Internal helpers.
// ---------------------------------------------------------------------------

fn admit_sample_count(count: usize) -> Result<(), FieldError> {
    if count == 0 {
        return Err(FieldError::NoSamples);
    }
    if count > MAX_FIELD_SAMPLES {
        return Err(FieldError::TooManySamples {
            have: count,
            max: MAX_FIELD_SAMPLES,
        });
    }
    Ok(())
}

fn checkpoint(cx: &fs_exec::Cx<'_>, ordinal: usize, phase: &'static str) -> Result<(), FieldError> {
    if ordinal.is_multiple_of(AS_BUILT_POLL_STRIDE_POINTS) {
        cx.checkpoint()
            .map_err(|_| FieldError::Cancelled { phase })?;
    }
    Ok(())
}

/// `sqrt(g^T C g)`, the first-order standard deviation of a linear functional.
fn quadratic_form_sigma(
    covariance: &[[f64; 6]; 6],
    gradient: &[f64; 6],
) -> Result<f64, FieldError> {
    let mut variance = 0.0f64;
    for row in 0..6 {
        let mut inner = 0.0f64;
        for column in 0..6 {
            inner = covariance[row][column].mul_add(gradient[column], inner);
        }
        variance = gradient[row].mul_add(inner, variance);
    }
    if !variance.is_finite() {
        return Err(FieldError::ArithmeticOverflow {
            field: "pose quadratic form",
        });
    }
    // A symmetric positive-semidefinite form cannot be negative in exact
    // arithmetic; rounding can still produce a tiny negative near zero.
    Ok(variance.max(0.0).sqrt())
}

fn quadrature2(a: f64, b: f64) -> Result<f64, FieldError> {
    let value = a.mul_add(a, b * b).sqrt();
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FieldError::ArithmeticOverflow {
            field: "quadrature sum",
        })
    }
}

fn quadrature3(a: f64, b: f64, c: f64) -> Result<f64, FieldError> {
    let value = a.mul_add(a, b.mul_add(b, c * c)).sqrt();
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FieldError::ArithmeticOverflow {
            field: "quadrature sum",
        })
    }
}

fn normalize_axis(value: f64, center: f64, half_extent: f64) -> f64 {
    (value - center) / half_extent
}

/// Exact small-integer power by repeated multiplication.
///
/// `f64::powi` is not pinned across build modes in this workspace, and a
/// monomial basis entry that changes with the optimization level would move
/// every fitted coefficient.
fn powi_exact(value: f64, power: usize) -> f64 {
    let mut result = 1.0f64;
    for _ in 0..power {
        result *= value;
    }
    result
}

/// Total-degree monomial exponents `(i, j)` in canonical order: by total
/// degree, then by the first exponent.
fn basis_exponents(order: u8) -> Vec<(usize, usize)> {
    let order = order as usize;
    let mut exponents = Vec::new();
    for degree in 0..=order {
        for i in 0..=degree {
            exponents.push((i, degree - i));
        }
    }
    exponents
}

struct FieldEncoder {
    hasher: fs_blake3::Blake3,
}

impl FieldEncoder {
    fn new(kind: &[u8]) -> Self {
        let mut encoder = Self {
            hasher: fs_blake3::Blake3::new(),
        };
        encoder.bytes(FIELD_IDENTITY_DOMAIN.as_bytes());
        encoder.bytes(kind);
        encoder.u32(FIELD_SCHEMA_VERSION);
        encoder
    }

    fn bytes(&mut self, bytes: &[u8]) {
        self.hasher.update(&(bytes.len() as u64).to_le_bytes());
        self.hasher.update(bytes);
    }

    fn u32(&mut self, value: u32) {
        self.bytes(&value.to_le_bytes());
    }

    fn u64(&mut self, value: u64) {
        self.bytes(&value.to_le_bytes());
    }

    fn f64(&mut self, value: f64) {
        self.bytes(&canonical_zero(value).to_bits().to_le_bytes());
    }

    fn probe(&mut self, probe: &ProbeModel) {
        for axis in probe.direction {
            self.f64(axis);
        }
        self.f64(probe.grazing_floor);
        self.f64(probe.lateral_sigma);
        self.f64(probe.range_sigma);
    }

    fn coverage(&mut self, coverage: &CoveragePolicy) {
        self.f64(coverage.level());
        self.f64(coverage.factor());
    }

    fn refusals(&mut self, grazing: &[usize]) {
        self.u64(grazing.len() as u64);
        for index in grazing {
            self.u64(*index as u64);
        }
    }

    fn finish(self) -> fs_blake3::ContentHash {
        let preimage = self.hasher.finalize();
        fs_blake3::hash_domain(FIELD_IDENTITY_DOMAIN, preimage.as_bytes())
    }
}

fn field_identity(
    kind: &[u8],
    registration_ref: fs_blake3::ContentHash,
    coverage: &CoveragePolicy,
    probe: &ProbeModel,
    points: &[FieldPoint],
    grazing: &[usize],
) -> fs_blake3::ContentHash {
    let mut encoder = FieldEncoder::new(kind);
    encoder.bytes(registration_ref.as_bytes());
    encoder.coverage(coverage);
    encoder.probe(probe);
    encoder.u64(points.len() as u64);
    for point in points {
        encoder.f64(point.deviation);
        encoder.f64(point.incidence_cosine);
        encoder.f64(point.scan_sigma);
        encoder.f64(point.registration_sigma);
        encoder.f64(point.ambiguity_sigma);
        encoder.f64(point.half_width);
    }
    encoder.refusals(grazing);
    encoder.finish()
}

fn thickness_identity(
    registration_ref: fs_blake3::ContentHash,
    coverage: &CoveragePolicy,
    probe: &ProbeModel,
    points: &[ThicknessPoint],
    grazing: &[usize],
) -> fs_blake3::ContentHash {
    let mut encoder = FieldEncoder::new(THICKNESS_IDENTITY_KIND);
    encoder.bytes(registration_ref.as_bytes());
    encoder.coverage(coverage);
    encoder.probe(probe);
    encoder.u64(points.len() as u64);
    for point in points {
        encoder.f64(point.nominal);
        encoder.f64(point.thickness);
        encoder.f64(point.registration_sigma);
        encoder.f64(point.measurement_sigma);
        encoder.f64(point.half_width);
    }
    encoder.refusals(grazing);
    encoder.finish()
}
