//! Discrepancy models and model bracketing (patch Rev B mechanisms 2–3).
//!
//! A [`DiscrepancyModel`] is fit from paired two-fidelity evaluations (the
//! ledger's accumulating corpus) and answers "how wrong is the cheap model
//! HERE" — refusing with a teaching [`OutOfDomain`] when asked outside the
//! region it has data for or with a parameter key set that differs from the
//! exact training schema (the surrogate out-of-distribution guard). v1 is
//! deliberately statistics-light: an observed parameter box plus
//! mean/max relative discrepancy — honest bookkeeping, not learning
//! (learned discrepancy models arrive with FrankenTorch — CONTRACT
//! no-claims).
//!
//! A [`ModelBracket`] handles weakly-understood physics (the vessel
//! flagship's contact-line mitigation): run EVERY plausible model, report
//! the QoI spread as an enclosure plus a model-form band — sensitivity to
//! the modeling choice, not pretended certainty.

use crate::{
    Evidence, ModelEvidence, NumericalCertificate, ProvenanceHash, SensitivitySummary,
    StatisticalCertificate, ValidityDomain, color_identity_reason, color_leaf_identity_reason,
};
use core::fmt;
use std::collections::BTreeMap;

const MAX_TRAINING_PAIRS: usize = 65_536;
const MAX_TRAINING_PARAMETERS: usize = 1_024;
const MAX_TRAINING_COORDINATES: usize = 1_048_576;
pub(crate) const MIN_BRACKET_MEMBERS: usize = 2;
pub(crate) const MAX_BRACKET_MEMBERS: usize = 1_024;

/// One paired two-fidelity evaluation at a parameter point.
#[derive(Debug, Clone, PartialEq)]
pub struct FidelityPair {
    /// Where in parameter space the pair was evaluated.
    pub params: BTreeMap<String, f64>,
    /// Low-fidelity QoI.
    pub lo_fi: f64,
    /// High-fidelity QoI (the reference).
    pub hi_fi: f64,
}

/// The OBSERVED discrepancy statistics of the training corpus.
///
/// Both fields are sample statistics over the training pairs. Neither is a
/// bound on the discrepancy anywhere else: nothing in v1 establishes
/// smoothness, sample density, or a Lipschitz constant between training
/// points, so the observed maximum is not a conservative band at an unvisited
/// parameter point (bead frankensim-extreal-program-f85xj.2.7). The field
/// names say `observed` so a reader cannot infer otherwise.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiscrepancyBand {
    /// Mean relative discrepancy ACROSS THE TRAINING PAIRS.
    pub mean_observed_rel: f64,
    /// Worst relative discrepancy OBSERVED AT A TRAINING PAIR.
    pub max_observed_rel: f64,
}

/// Structured discrepancy-model or model-bracket refusal.
///
/// The historical name is retained for API compatibility; bracket construction
/// uses the same error so every public discrepancy-evidence path fails with a
/// nameable, deterministic reason.
#[derive(Debug, Clone, PartialEq)]
pub enum FitError {
    /// No paired observations were supplied.
    EmptyTrainingSet,
    /// The paired corpus exceeds the bounded v1 bookkeeping budget.
    TooManyTrainingPairs {
        /// Supplied pair count.
        count: usize,
        /// Maximum admitted pair count.
        maximum: usize,
    },
    /// A pair did not declare any parameter coordinates.
    EmptyParameterSchema,
    /// The parameter schema exceeds the bounded validity-domain budget.
    TooManyTrainingParameters {
        /// Supplied parameter count.
        count: usize,
        /// Maximum admitted parameter count.
        maximum: usize,
    },
    /// Pair count times parameter count exceeds the bounded fit-work budget.
    TooManyTrainingCoordinates {
        /// Supplied pair count.
        pairs: usize,
        /// Parameters per pair.
        parameters: usize,
        /// Maximum admitted coordinate count.
        maximum: usize,
    },
    /// A pair's parameter names differ from the first pair's schema.
    InconsistentParameterSchema {
        /// Zero-based pair index.
        pair_index: usize,
    },
    /// A parameter name cannot become a validity-domain identity.
    InvalidParameterIdentity {
        /// Zero-based pair index.
        pair_index: usize,
        /// Shared identity-grammar rejection reason.
        reason: &'static str,
    },
    /// A training pair contains a NaN or infinite QoI.
    NonFiniteTrainingQoi {
        /// Zero-based pair index.
        pair_index: usize,
    },
    /// A training coordinate is NaN or infinite.
    NonFiniteTrainingParameter {
        /// Zero-based pair index.
        pair_index: usize,
        /// Bounded, already-validated parameter identity.
        param: String,
    },
    /// `evidence_at` was given an unusable model-card identity.
    InvalidCardIdentity {
        /// Shared leaf-identity rejection reason.
        reason: &'static str,
    },
    /// An otherwise valid query is outside the observed parameter box.
    QueryOutOfDomain(OutOfDomain),
    /// A bracket cannot measure model-form spread with fewer than two models.
    TooFewBracketMembers {
        /// Number of admitted unique members.
        count: usize,
        /// Required minimum.
        minimum: usize,
    },
    /// A bracket exceeded its bounded member budget.
    TooManyBracketMembers {
        /// Maximum admitted member count.
        maximum: usize,
    },
    /// The exact member name could not reserve its bounded storage.
    BracketMemberNameAllocationFailed {
        /// Exact UTF-8 byte capacity requested.
        requested_bytes: usize,
    },
    /// The canonical member list could not reserve one more slot.
    BracketMemberListAllocationFailed {
        /// Exact resulting member capacity requested.
        requested_members: usize,
    },
    /// A member name cannot become a model-card identity.
    InvalidBracketMemberIdentity {
        /// Shared leaf-identity rejection reason.
        reason: &'static str,
    },
    /// Two member QoIs claimed the same model identity.
    DuplicateBracketMember {
        /// The duplicate, bounded identity.
        name: String,
    },
    /// A member QoI was NaN or infinite.
    NonFiniteBracketQoi {
        /// The member's bounded identity.
        name: String,
    },
    /// A declared fill distance is not a usable normalized radius.
    InvalidFillDistance {
        /// The offending declaration.
        fill_distance: f64,
    },
    /// The query lies inside the observed bounding box but farther from every
    /// training point than the declared fill distance, so the observed band is
    /// extrapolation rather than interpolation.
    QueryBeyondFillDistance {
        /// Normalized (per-axis box-extent) distance to the nearest pair.
        distance: f64,
        /// The declared normalized fill distance.
        fill_distance: f64,
    },
}

impl fmt::Display for FitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FitError::EmptyTrainingSet => {
                write!(f, "discrepancy fit refused: the training set is empty")
            }
            FitError::TooManyTrainingPairs { count, maximum } => write!(
                f,
                "discrepancy fit refused: {count} training pairs exceed the bounded maximum {maximum}"
            ),
            FitError::EmptyParameterSchema => write!(
                f,
                "discrepancy fit refused: training pairs must share a non-empty parameter schema"
            ),
            FitError::TooManyTrainingParameters { count, maximum } => write!(
                f,
                "discrepancy fit refused: {count} parameters exceed the bounded maximum {maximum}"
            ),
            FitError::TooManyTrainingCoordinates {
                pairs,
                parameters,
                maximum,
            } => write!(
                f,
                "discrepancy fit refused: {pairs} pairs x {parameters} parameters exceed the bounded {maximum}-coordinate work budget"
            ),
            FitError::InconsistentParameterSchema { pair_index } => write!(
                f,
                "discrepancy fit refused: training pair {pair_index} does not have exactly the first pair's parameter schema"
            ),
            FitError::InvalidParameterIdentity { pair_index, reason } => write!(
                f,
                "discrepancy fit refused: training pair {pair_index} has an invalid parameter identity ({reason})"
            ),
            FitError::NonFiniteTrainingQoi { pair_index } => write!(
                f,
                "discrepancy fit refused: training pair {pair_index} has a non-finite QoI"
            ),
            FitError::NonFiniteTrainingParameter { pair_index, param } => write!(
                f,
                "discrepancy fit refused: training pair {pair_index} parameter `{param}` is non-finite"
            ),
            FitError::InvalidCardIdentity { reason } => write!(
                f,
                "discrepancy evidence refused: the model-card identity is invalid ({reason})"
            ),
            FitError::QueryOutOfDomain(error) => error.fmt(f),
            FitError::TooFewBracketMembers { count, minimum } => write!(
                f,
                "model bracket refused: {count} unique model member(s) cannot measure model-choice spread; supply at least {minimum}"
            ),
            FitError::TooManyBracketMembers { maximum } => write!(
                f,
                "model bracket refused: member count exceeds the bounded maximum {maximum}"
            ),
            FitError::BracketMemberNameAllocationFailed { requested_bytes } => write!(
                f,
                "model bracket refused: could not reserve {requested_bytes} bytes for the member name"
            ),
            FitError::BracketMemberListAllocationFailed { requested_members } => write!(
                f,
                "model bracket refused: could not reserve storage for {requested_members} members"
            ),
            FitError::InvalidBracketMemberIdentity { reason } => write!(
                f,
                "model bracket refused: a member identity is invalid ({reason})"
            ),
            FitError::DuplicateBracketMember { name } => write!(
                f,
                "model bracket refused: duplicate model identity `{name}` is ambiguous"
            ),
            FitError::NonFiniteBracketQoi { name } => write!(
                f,
                "model bracket refused: member `{name}` has a non-finite QoI"
            ),
            FitError::InvalidFillDistance { fill_distance } => write!(
                f,
                "discrepancy fill distance refused: {fill_distance} is not a finite radius in (0, 1]"
            ),
            FitError::QueryBeyondFillDistance {
                distance,
                fill_distance,
            } => write!(
                f,
                "discrepancy evidence refused: the query is {distance} (normalized) from the \
                 nearest training pair, beyond the declared fill distance {fill_distance} — the \
                 observed band there would be extrapolation, not interpolation; gather pairs \
                 nearer the query or declare the wider assumption explicitly"
            ),
        }
    }
}

impl core::error::Error for FitError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            FitError::QueryOutOfDomain(error) => Some(error),
            _ => None,
        }
    }
}

/// The out-of-distribution refusal, naming the violated parameter (the
/// diagnosis an agent needs to decide between gathering data and
/// escalating fidelity).
#[derive(Debug, Clone, PartialEq)]
pub struct OutOfDomain {
    /// The parameter outside the trained box, missing from the query, or absent
    /// from the training schema.
    pub param: String,
    /// The queried value (`None` = the query omitted the parameter).
    pub value: Option<f64>,
    /// The trained box for that parameter. `None` means the query supplied an
    /// unexpected parameter that was never part of the training schema.
    pub trained: Option<(f64, f64)>,
}

impl fmt::Display for OutOfDomain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match (self.value, self.trained) {
            (Some(v), Some(trained)) => write!(
                f,
                "discrepancy query out of domain: `{}` = {v} lies outside the trained box \
                 [{}, {}] — the band would be extrapolation; gather pairs there or escalate \
                 fidelity",
                self.param, trained.0, trained.1
            ),
            (None, Some(trained)) => write!(
                f,
                "discrepancy query out of domain: `{}` was not supplied but the model was \
                 trained on it (trained box [{}, {}])",
                self.param, trained.0, trained.1
            ),
            (Some(v), None) => write!(
                f,
                "discrepancy query out of domain: unexpected parameter `{}` = {v} was not \
                 present in the exact training schema; remove it or refit the discrepancy model \
                 with that dimension",
                self.param
            ),
            (None, None) => write!(
                f,
                "discrepancy query out of domain: unexpected parameter `{}` was not present in \
                 the exact training schema",
                self.param
            ),
        }
    }
}

impl core::error::Error for OutOfDomain {}

/// A two-fidelity discrepancy model: observed box + mean/max relative
/// discrepancy (v1 bookkeeping — see module docs).
#[derive(Debug, Clone, PartialEq)]
pub struct DiscrepancyModel {
    observed: ValidityDomain,
    trained_bounds: BTreeMap<String, (f64, f64)>,
    /// Training coordinates in `trained_bounds` key order — the support set the
    /// fill-distance guard measures against.
    trained_points: Vec<Vec<f64>>,
    band: DiscrepancyBand,
    pairs: usize,
    /// Caller-DECLARED normalized fill distance. `0.0` (the default) admits
    /// only exact training points: nothing about a bounding box establishes
    /// that a point between samples is supported.
    fill_distance: f64,
}

fn validate_training_shape(pair_count: usize, parameter_count: usize) -> Result<(), FitError> {
    if pair_count == 0 {
        return Err(FitError::EmptyTrainingSet);
    }
    if pair_count > MAX_TRAINING_PAIRS {
        return Err(FitError::TooManyTrainingPairs {
            count: pair_count,
            maximum: MAX_TRAINING_PAIRS,
        });
    }
    if parameter_count == 0 {
        return Err(FitError::EmptyParameterSchema);
    }
    if parameter_count > MAX_TRAINING_PARAMETERS {
        return Err(FitError::TooManyTrainingParameters {
            count: parameter_count,
            maximum: MAX_TRAINING_PARAMETERS,
        });
    }
    if pair_count
        .checked_mul(parameter_count)
        .is_none_or(|coordinates| coordinates > MAX_TRAINING_COORDINATES)
    {
        return Err(FitError::TooManyTrainingCoordinates {
            pairs: pair_count,
            parameters: parameter_count,
            maximum: MAX_TRAINING_COORDINATES,
        });
    }
    Ok(())
}

impl DiscrepancyModel {
    /// Fit from paired evaluations. The observed box is the per-parameter
    /// min/max over training points; the band is mean/max of
    /// `|hi - lo| / max(|hi|, tiny)`.
    ///
    /// # Errors
    /// [`FitError`] when the corpus is empty, oversized, non-finite, uses an
    /// invalid parameter identity, or does not have one exact shared non-empty
    /// parameter schema.
    pub fn fit(pairs: &[FidelityPair]) -> Result<Self, FitError> {
        let Some(first_pair) = pairs.first() else {
            return Err(FitError::EmptyTrainingSet);
        };
        let parameter_count = first_pair.params.len();
        validate_training_shape(pairs.len(), parameter_count)?;

        let mut bounds: BTreeMap<String, (f64, f64)> = BTreeMap::new();
        for (param, &value) in &first_pair.params {
            if let Some(reason) = color_identity_reason(param) {
                return Err(FitError::InvalidParameterIdentity {
                    pair_index: 0,
                    reason,
                });
            }
            bounds.insert(param.clone(), (value, value));
        }

        let mut mean_observed_rel = 0.0_f64;
        let mut max_observed_rel = 0.0_f64;
        for (pair_index, pair) in pairs.iter().enumerate() {
            if pair.params.len() != bounds.len() || !pair.params.keys().eq(bounds.keys()) {
                return Err(FitError::InconsistentParameterSchema { pair_index });
            }
            if !pair.lo_fi.is_finite() || !pair.hi_fi.is_finite() {
                return Err(FitError::NonFiniteTrainingQoi { pair_index });
            }
            for (param, &value) in &pair.params {
                if !value.is_finite() {
                    return Err(FitError::NonFiniteTrainingParameter {
                        pair_index,
                        param: param.clone(),
                    });
                }
                let Some((lo, hi)) = bounds.get_mut(param) else {
                    return Err(FitError::InconsistentParameterSchema { pair_index });
                };
                *lo = lo.min(value);
                *hi = hi.max(value);
            }

            let rel = (pair.hi_fi - pair.lo_fi).abs() / pair.hi_fi.abs().max(f64::MIN_POSITIVE);
            max_observed_rel = max_observed_rel.max(rel);
            if rel.is_infinite() {
                mean_observed_rel = f64::INFINITY;
            } else if mean_observed_rel.is_finite() {
                let count = (pair_index + 1) as f64;
                mean_observed_rel += (rel - mean_observed_rel) / count;
                mean_observed_rel = mean_observed_rel.min(max_observed_rel);
            }
        }
        let mut observed = ValidityDomain::unconstrained();
        for (k, &(lo, hi)) in &bounds {
            observed = observed.with(k.clone(), lo, hi);
        }
        let trained_points = pairs
            .iter()
            .map(|pair| {
                bounds
                    .keys()
                    .map(|param| pair.params[param.as_str()])
                    .collect()
            })
            .collect();
        Ok(DiscrepancyModel {
            observed,
            trained_bounds: bounds,
            trained_points,
            band: DiscrepancyBand {
                mean_observed_rel,
                max_observed_rel,
            },
            pairs: pairs.len(),
            fill_distance: 0.0,
        })
    }

    /// DECLARE the fill distance under which the observed band may be read as
    /// interpolation: a normalized (per-axis box-extent) radius in `(0, 1]`
    /// around each training point.
    ///
    /// This is an ASSUMPTION the caller takes responsibility for, not a fact
    /// the corpus establishes, and it travels in the assumptions of every
    /// [`ModelEvidence`] the model emits. Without it a model admits only
    /// exact training points, because a bounding box says nothing about the
    /// unvisited interior (bead frankensim-extreal-program-f85xj.2.7).
    ///
    /// # Errors
    /// [`FitError::InvalidFillDistance`] when the radius is not finite in
    /// `(0, 1]`.
    pub fn with_declared_fill_distance(mut self, fill_distance: f64) -> Result<Self, FitError> {
        if !fill_distance.is_finite() || fill_distance <= 0.0 || fill_distance > 1.0 {
            return Err(FitError::InvalidFillDistance { fill_distance });
        }
        self.fill_distance = fill_distance;
        Ok(self)
    }

    /// The declared normalized fill distance (`0.0` = none declared).
    #[must_use]
    pub fn declared_fill_distance(&self) -> f64 {
        self.fill_distance
    }

    /// Normalized distance from `point` to the NEAREST training pair, using
    /// the per-axis box extent as the scale and the max-norm across axes. A
    /// zero-extent axis contributes `0.0` (the box check already pinned it).
    ///
    /// The point is assumed to have passed [`DiscrepancyModel::query`], so
    /// every trained axis is present and finite.
    fn nearest_training_distance(&self, point: &BTreeMap<String, f64>) -> f64 {
        let scales: Vec<f64> = self
            .trained_bounds
            .values()
            .map(|&(lo, hi)| hi - lo)
            .collect();
        let query: Vec<f64> = self
            .trained_bounds
            .keys()
            .map(|param| point[param.as_str()])
            .collect();
        self.trained_points
            .iter()
            .map(|training| {
                training
                    .iter()
                    .zip(&query)
                    .zip(&scales)
                    .map(|((&t, &q), &scale)| {
                        if scale > 0.0 && scale.is_finite() {
                            ((q - t) / scale).abs()
                        } else {
                            0.0
                        }
                    })
                    .fold(0.0_f64, f64::max)
            })
            .fold(f64::INFINITY, f64::min)
    }

    /// Number of training pairs.
    #[must_use]
    pub fn pairs(&self) -> usize {
        self.pairs
    }

    /// The observed (trained) parameter box.
    #[must_use]
    pub fn trained_domain(&self) -> &ValidityDomain {
        &self.observed
    }

    /// The OBSERVED training statistics, for a `point` inside the observed
    /// bounding box.
    ///
    /// The box test is bounding-box membership, not an in-distribution test:
    /// the axis-aligned hull of the training coordinates contains corners and
    /// interiors where no pair was ever evaluated. Reading the returned
    /// statistics as a band AT `point` is an interpolation assumption; the
    /// guard for that lives on [`DiscrepancyModel::evidence_at`], which
    /// requires a declared fill distance.
    ///
    /// # Errors
    /// [`OutOfDomain`] naming the first missing, unexpected, non-finite, or
    /// out-of-box parameter (BTreeMap order — deterministic diagnosis). The
    /// query key set must equal the training schema exactly; silently ignoring
    /// an extra physical dimension would make even the box claim unsound.
    pub fn query(&self, point: &BTreeMap<String, f64>) -> Result<DiscrepancyBand, OutOfDomain> {
        for (param, &(lo, hi)) in &self.trained_bounds {
            match point.get(param) {
                None => {
                    return Err(OutOfDomain {
                        param: param.clone(),
                        value: None,
                        trained: Some((lo, hi)),
                    });
                }
                Some(&v) if !v.is_finite() || v < lo || v > hi => {
                    return Err(OutOfDomain {
                        param: param.clone(),
                        value: Some(v),
                        trained: Some((lo, hi)),
                    });
                }
                Some(_) => {}
            }
        }
        if let Some((param, &value)) = point
            .iter()
            .find(|(param, _)| !self.trained_bounds.contains_key(param.as_str()))
        {
            return Err(OutOfDomain {
                param: param.clone(),
                value: Some(value),
                trained: None,
            });
        }
        Ok(self.band)
    }

    /// Model evidence for the LOW-fidelity model at a SUPPORTED point,
    /// carrying the observed maximum as its band.
    ///
    /// `in_domain: true` here means exactly two things: `point` is inside the
    /// observed bounding box, AND it is within the caller-declared fill
    /// distance of an actual training pair. Both are recorded in the returned
    /// assumptions. Without a declared fill distance only exact training
    /// points qualify — the bounding box alone never established that the
    /// unvisited interior is supported, and admitting it published a sampled
    /// maximum as a model-form band that could then certify (bead
    /// frankensim-extreal-program-f85xj.2.7).
    ///
    /// # Errors
    /// [`FitError::InvalidCardIdentity`] for an unusable evidence identity,
    /// [`FitError::QueryOutOfDomain`] as in [`DiscrepancyModel::query`], or
    /// [`FitError::QueryBeyondFillDistance`] when the query is inside the box
    /// but not supported by any training pair.
    pub fn evidence_at(
        &self,
        card_name: &str,
        point: &BTreeMap<String, f64>,
    ) -> Result<ModelEvidence, FitError> {
        if let Some(reason) = color_leaf_identity_reason(card_name) {
            return Err(FitError::InvalidCardIdentity { reason });
        }
        let band = self.query(point).map_err(FitError::QueryOutOfDomain)?;
        let distance = self.nearest_training_distance(point);
        if distance.is_nan() || distance > self.fill_distance {
            return Err(FitError::QueryBeyondFillDistance {
                distance,
                fill_distance: self.fill_distance,
            });
        }
        Ok(ModelEvidence {
            cards: vec![card_name.to_string()],
            assumptions: vec![
                format!(
                    "lo-fi accuracy from a {}-pair two-fidelity discrepancy model",
                    self.pairs
                ),
                format!(
                    "band is the OBSERVED training maximum read at a point {distance} \
                     (normalized) from the nearest pair, under a declared fill distance of {}",
                    self.fill_distance
                ),
            ],
            validity: self.observed.clone(),
            discrepancy_rel: band.max_observed_rel,
            in_domain: true,
        })
    }
}

/// Model bracketing: N plausible models of weakly-understood physics, one
/// QoI each; the evidence is the SPREAD, not a pretended point value.
///
/// Admitted members are stored in exact model-name order so construction order
/// is presentation-only and every equivalent bracket has one structural form.
#[derive(Debug, Clone, PartialEq)]
pub struct ModelBracket {
    members: Vec<(String, f64)>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct BracketReservationError;

impl ModelBracket {
    /// Start a bracket.
    #[must_use]
    pub fn new() -> Self {
        ModelBracket {
            members: Vec::new(),
        }
    }

    fn push_member_with_reservations<RN, RM>(
        &mut self,
        name: &str,
        qoi: f64,
        reserve_name: RN,
        reserve_members: RM,
    ) -> Result<(), FitError>
    where
        RN: FnOnce(&mut String, usize) -> Result<(), BracketReservationError>,
        RM: FnOnce(&mut Vec<(String, f64)>, usize) -> Result<(), BracketReservationError>,
    {
        if self.members.len() >= MAX_BRACKET_MEMBERS {
            return Err(FitError::TooManyBracketMembers {
                maximum: MAX_BRACKET_MEMBERS,
            });
        }
        if let Some(reason) = color_leaf_identity_reason(name) {
            return Err(FitError::InvalidBracketMemberIdentity { reason });
        }
        let mut owned_name = String::new();
        reserve_name(&mut owned_name, name.len()).map_err(|BracketReservationError| {
            FitError::BracketMemberNameAllocationFailed {
                requested_bytes: name.len(),
            }
        })?;
        owned_name.push_str(name);
        if !qoi.is_finite() {
            return Err(FitError::NonFiniteBracketQoi { name: owned_name });
        }
        let Err(position) = self
            .members
            .binary_search_by(|(member, _)| member.as_str().cmp(owned_name.as_str()))
        else {
            return Err(FitError::DuplicateBracketMember { name: owned_name });
        };
        let requested_members = self.members.len() + 1;
        reserve_members(&mut self.members, 1).map_err(|BracketReservationError| {
            FitError::BracketMemberListAllocationFailed { requested_members }
        })?;
        self.members.insert(position, (owned_name, qoi));
        Ok(())
    }

    fn push_member(&mut self, name: &str, qoi: f64) -> Result<(), FitError> {
        self.push_member_with_reservations(
            name,
            qoi,
            |value, additional| {
                value
                    .try_reserve_exact(additional)
                    .map_err(|_| BracketReservationError)
            },
            |values, additional| {
                values
                    .try_reserve_exact(additional)
                    .map_err(|_| BracketReservationError)
            },
        )
    }

    /// Exact canonical member-name/QoI rows for the strong-identity helper.
    pub(crate) fn identity_members(&self) -> &[(String, f64)] {
        &self.members
    }

    /// Add a member model's QoI and report admission failure immediately.
    ///
    /// # Errors
    /// [`FitError`] when the member identity or QoI is invalid, duplicated, or
    /// exceeds the bounded bracket-member budget, or when bounded name/member
    /// storage cannot be reserved.
    pub fn try_with_member(mut self, name: impl AsRef<str>, qoi: f64) -> Result<Self, FitError> {
        self.push_member(name.as_ref(), qoi)?;
        Ok(self)
    }

    /// Collapse the bracket into evidence: the numerical slice encloses
    /// every member's QoI (outward-rounded); the model slice records the
    /// bracketing and carries the relative spread as its band; the
    /// representative value is the member MIDRANGE (deterministic). At
    /// least two uniquely named finite members are required.
    ///
    /// The model slice is `in_domain: false`. A bracket holds only
    /// `(name, qoi)` rows: no member validity domain, evaluation point, or
    /// model card is ever supplied, so nothing here can establish that any
    /// member was used inside its stated validity range. Bracket spread
    /// quantifies model-CHOICE sensitivity, which is a different claim — and
    /// a hard-coded `in_domain: true` was exactly what let this evidence
    /// certify (bead frankensim-extreal-program-f85xj.2.10). Use
    /// [`ModelBracket::evidence_at`] to supply the domain and evaluation point
    /// the verdict actually needs.
    ///
    /// # Errors
    /// [`FitError::TooFewBracketMembers`] when fewer than two valid models were
    /// supplied.
    pub fn evidence(&self, provenance: ProvenanceHash) -> Result<Evidence<f64>, FitError> {
        self.evidence_with_domain(provenance, ValidityDomain::unconstrained(), false)
    }

    /// Bracket evidence whose `in_domain` verdict is COMPUTED from a declared
    /// member validity domain and the evaluation point, instead of asserted.
    ///
    /// # Errors
    /// [`FitError::TooFewBracketMembers`] when fewer than two valid models were
    /// supplied.
    pub fn evidence_at(
        &self,
        provenance: ProvenanceHash,
        validity: ValidityDomain,
        point: &BTreeMap<String, f64>,
    ) -> Result<Evidence<f64>, FitError> {
        let in_domain = !validity.bounds().is_empty() && validity.contains(point);
        self.evidence_with_domain(provenance, validity, in_domain)
    }

    fn evidence_with_domain(
        &self,
        provenance: ProvenanceHash,
        validity: ValidityDomain,
        in_domain: bool,
    ) -> Result<Evidence<f64>, FitError> {
        if self.members.len() < MIN_BRACKET_MEMBERS {
            return Err(FitError::TooFewBracketMembers {
                count: self.members.len(),
                minimum: MIN_BRACKET_MEMBERS,
            });
        }
        let mut qois = self.members.iter().map(|(_, qoi)| *qoi);
        let Some(first) = qois.next() else {
            return Err(FitError::TooFewBracketMembers {
                count: 0,
                minimum: MIN_BRACKET_MEMBERS,
            });
        };
        let (lo, hi) = qois.fold((first, first), |(lo, hi), qoi| (lo.min(qoi), hi.max(qoi)));
        let mid = f64::midpoint(lo, hi);
        let spread_rel = (hi - lo) / mid.abs().max(f64::MIN_POSITIVE);
        let names: Vec<String> = self.members.iter().map(|(n, _)| n.clone()).collect();
        let mut sensitivity = SensitivitySummary::default();
        sensitivity
            .d_qoi
            .insert("model-choice(bracket-spread)".to_string(), hi - lo);
        let enclosure_lo = lo.next_down();
        let enclosure_hi = hi.next_up();
        Ok(Evidence {
            value: mid,
            qoi: mid,
            numerical: NumericalCertificate::enclosure(
                if enclosure_lo.is_finite() {
                    enclosure_lo
                } else {
                    lo
                },
                if enclosure_hi.is_finite() {
                    enclosure_hi
                } else {
                    hi
                },
            ),
            statistical: StatisticalCertificate::None,
            model: ModelEvidence {
                assumptions: vec![
                    format!("model-bracketed over: {}", names.join(", ")),
                    if in_domain {
                        "every bracket member was declared valid at the supplied evaluation point"
                            .to_string()
                    } else {
                        "bracket spread measures model-CHOICE sensitivity only; no member \
                         validity domain was supplied, so no member is claimed to have been \
                         used inside its stated validity range"
                            .to_string()
                    },
                ],
                cards: names,
                validity,
                discrepancy_rel: spread_rel,
                in_domain,
            },
            sensitivity,
            provenance,
            adjoint_ref: None,
        })
    }
}

impl Default for ModelBracket {
    fn default() -> Self {
        ModelBracket::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pt(pairs: &[(&str, f64)]) -> BTreeMap<String, f64> {
        pairs.iter().map(|&(k, v)| (k.to_string(), v)).collect()
    }

    fn bracket(members: &[(&str, f64)]) -> ModelBracket {
        members
            .iter()
            .try_fold(ModelBracket::new(), |bracket, &(name, qoi)| {
                bracket.try_with_member(name, qoi)
            })
            .expect("valid bracket fixture")
    }

    #[test]
    fn fit_refuses_empty_and_non_finite_training_sets() {
        let err = DiscrepancyModel::fit(&[]).expect_err("empty");
        assert!(matches!(&err, FitError::EmptyTrainingSet));
        assert!(err.to_string().contains("empty"), "{err}");
        let err = DiscrepancyModel::fit(&[FidelityPair {
            params: pt(&[("Re", 1e4)]),
            lo_fi: f64::NAN,
            hi_fi: 1.0,
        }])
        .expect_err("nan");
        assert!(err.to_string().contains("non-finite"), "{err}");
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = DiscrepancyModel::fit(&[FidelityPair {
                params: pt(&[("Re", value)]),
                lo_fi: 1.0,
                hi_fi: 1.0,
            }])
            .expect_err("non-finite parameter");
            assert!(err.to_string().contains("parameter"), "{err}");
        }
    }

    #[test]
    fn fit_requires_one_exact_bounded_parameter_schema() {
        let empty_schema = FidelityPair {
            params: BTreeMap::new(),
            lo_fi: 1.0,
            hi_fi: 1.0,
        };
        assert!(matches!(
            DiscrepancyModel::fit(&[empty_schema]),
            Err(FitError::EmptyParameterSchema)
        ));

        let inconsistent = [
            FidelityPair {
                params: pt(&[("Re", 1.0), ("Ma", 0.1)]),
                lo_fi: 1.0,
                hi_fi: 1.0,
            },
            FidelityPair {
                params: pt(&[("Re", 2.0)]),
                lo_fi: 1.0,
                hi_fi: 1.0,
            },
        ];
        assert!(matches!(
            DiscrepancyModel::fit(&inconsistent),
            Err(FitError::InconsistentParameterSchema { pair_index: 1 })
        ));

        for param in ["", " Re", "pending", "control\naxis"] {
            assert!(matches!(
                DiscrepancyModel::fit(&[FidelityPair {
                    params: pt(&[(param, 1.0)]),
                    lo_fi: 1.0,
                    hi_fi: 1.0,
                }]),
                Err(FitError::InvalidParameterIdentity { pair_index: 0, .. })
            ));
        }
        let oversized = "x".repeat(crate::MAX_COLOR_IDENTITY_BYTES + 1);
        assert!(matches!(
            DiscrepancyModel::fit(&[FidelityPair {
                params: pt(&[(oversized.as_str(), 1.0)]),
                lo_fi: 1.0,
                hi_fi: 1.0,
            }]),
            Err(FitError::InvalidParameterIdentity {
                pair_index: 0,
                reason: "too-long"
            })
        ));

        assert!(matches!(
            validate_training_shape(MAX_TRAINING_PAIRS + 1, 1),
            Err(FitError::TooManyTrainingPairs { .. })
        ));
        assert!(matches!(
            validate_training_shape(1, MAX_TRAINING_PARAMETERS + 1),
            Err(FitError::TooManyTrainingParameters { .. })
        ));
        assert!(matches!(
            validate_training_shape(
                MAX_TRAINING_PAIRS,
                MAX_TRAINING_COORDINATES / MAX_TRAINING_PAIRS + 1,
            ),
            Err(FitError::TooManyTrainingCoordinates { .. })
        ));
    }

    #[test]
    fn fit_mean_is_bounded_by_max_even_when_naive_sum_would_overflow() {
        let pairs = [
            FidelityPair {
                params: pt(&[("x", 0.0)]),
                lo_fi: -f64::MAX,
                hi_fi: 1.0,
            },
            FidelityPair {
                params: pt(&[("x", 1.0)]),
                lo_fi: -f64::MAX,
                hi_fi: 1.0,
            },
        ];
        let model = DiscrepancyModel::fit(&pairs).expect("finite extreme fit");
        let band = model.query(&pt(&[("x", 0.5)])).expect("in domain");
        assert!(band.mean_observed_rel.is_finite());
        assert!(band.mean_observed_rel <= band.max_observed_rel);

        let unbounded = DiscrepancyModel::fit(&[
            FidelityPair {
                params: pt(&[("x", 0.0)]),
                lo_fi: -f64::MAX,
                hi_fi: f64::MAX,
            },
            FidelityPair {
                params: pt(&[("x", 1.0)]),
                lo_fi: 1.0,
                hi_fi: 1.0,
            },
        ])
        .expect("derived overflow is an honest unbounded band");
        let band = unbounded.query(&pt(&[("x", 0.5)])).expect("in domain");
        assert!(band.mean_observed_rel.is_infinite() && band.max_observed_rel.is_infinite());
    }

    #[test]
    fn in_domain_queries_report_the_band_and_out_of_domain_refuses() {
        let pairs: Vec<FidelityPair> = (0..10)
            .map(|i| {
                let re = 1e4 + f64::from(i) * 1e4;
                FidelityPair {
                    params: pt(&[("Re", re)]),
                    lo_fi: 1.0 + 0.05 * f64::from(i % 3),
                    hi_fi: 1.0,
                }
            })
            .collect();
        let model = DiscrepancyModel::fit(&pairs).expect("fit");
        let band = model.query(&pt(&[("Re", 5e4)])).expect("in domain");
        assert!(band.max_observed_rel >= band.mean_observed_rel && band.max_observed_rel <= 0.2);
        let err = model.query(&pt(&[("Re", 1e6)])).expect_err("extrapolation");
        assert_eq!(err.param, "Re");
        assert!(err.to_string().contains("extrapolation"), "{err}");
        let err = model.query(&pt(&[("Ma", 0.1)])).expect_err("missing param");
        assert!(err.to_string().contains("not supplied"), "{err}");
        let err = model
            .query(&pt(&[("Re", 5e4), ("Mach", 0.1)]))
            .expect_err("an untrained query dimension must not be ignored");
        assert_eq!(err.param, "Mach");
        assert_eq!(err.trained, None);
        assert!(err.to_string().contains("exact training schema"), "{err}");
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = model
                .query(&pt(&[("Re", value)]))
                .expect_err("non-finite query");
            assert_eq!(err.param, "Re");
            assert_eq!(err.value.map(f64::to_bits), Some(value.to_bits()));
        }
        assert!(matches!(
            model.evidence_at("pending", &pt(&[("Re", 5e4)])),
            Err(FitError::InvalidCardIdentity {
                reason: "placeholder"
            })
        ));
        assert!(matches!(
            model.evidence_at("panel-vs-les", &pt(&[("Re", 1e6)])),
            Err(FitError::QueryOutOfDomain(OutOfDomain { ref param, .. }))
                if param == "Re"
        ));
    }

    #[test]
    fn brackets_enclose_every_member_and_report_the_spread() {
        let bracket = bracket(&[
            ("contact-angle-60", 0.90),
            ("contact-angle-90", 1.00),
            ("contact-angle-120", 1.16),
        ]);
        let ev = bracket
            .evidence(ProvenanceHash::of_bytes(b"vessel-lip"))
            .expect("nonempty bracket");
        assert!(ev.numerical.lo <= 0.90 && ev.numerical.hi >= 1.16);
        assert!(
            ev.model.discrepancy_rel > 0.2,
            "{}",
            ev.model.discrepancy_rel
        );
        assert!(ev.model.assumptions[0].contains("model-bracketed"));
        // A bracket supplies no member validity domain, so the colour carries
        // no defensible spread claim (bead f85xj.2.10). It used to publish the
        // bracket spread as the dispersion off a hard-coded `in_domain: true`.
        assert!(matches!(
            crate::color_of(&ev.numerical, &ev.model),
            crate::Color::Estimated { dispersion, .. } if dispersion.is_infinite()
        ));
        assert!(matches!(
            ModelBracket::new().evidence(ProvenanceHash(0)),
            Err(FitError::TooFewBracketMembers {
                count: 0,
                minimum: MIN_BRACKET_MEMBERS
            })
        ));
        let single = ModelBracket::new()
            .try_with_member("only-model", 1.0)
            .expect("valid member");
        assert!(matches!(
            single.evidence(ProvenanceHash(0)),
            Err(FitError::TooFewBracketMembers { count: 1, .. })
        ));
    }

    #[test]
    fn bracket_refuses_non_finite_duplicate_and_invalid_members() {
        for qoi in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                ModelBracket::new().try_with_member("bad-model", qoi),
                Err(FitError::NonFiniteBracketQoi { ref name }) if name == "bad-model"
            ));
        }
        let one_member = ModelBracket::new()
            .try_with_member("same-model", 1.0)
            .expect("first identity");
        assert!(matches!(
            one_member.try_with_member("same-model", 2.0),
            Err(FitError::DuplicateBracketMember { ref name }) if name == "same-model"
        ));
        for name in ["", " pending", "pending", "derived:v2:forged"] {
            assert!(matches!(
                ModelBracket::new().try_with_member(name, 1.0),
                Err(FitError::InvalidBracketMemberIdentity { .. })
            ));
        }
        let two_members = bracket(&[("a", 1.0), ("b", 2.0)]);
        assert!(matches!(
            two_members.try_with_member("b", 3.0),
            Err(FitError::DuplicateBracketMember { .. })
        ));
        let full = ModelBracket {
            members: (0..MAX_BRACKET_MEMBERS)
                .map(|index| (format!("model-{index}"), index as f64))
                .collect(),
        };
        assert!(matches!(
            full.try_with_member("one-too-many", 0.0),
            Err(FitError::TooManyBracketMembers {
                maximum: MAX_BRACKET_MEMBERS
            })
        ));
    }

    #[test]
    fn bracket_allocation_refusals_are_typed_and_atomic() {
        let original = ModelBracket::new()
            .try_with_member("model-a", 1.0)
            .expect("baseline member");

        let mut name_failure = original.clone();
        let error = name_failure
            .push_member_with_reservations(
                "model-b",
                2.0,
                |_, _| Err(BracketReservationError),
                |_, _| panic!("member-list reservation must not follow name refusal"),
            )
            .expect_err("injected name reservation failure");
        assert_eq!(
            error,
            FitError::BracketMemberNameAllocationFailed {
                requested_bytes: "model-b".len(),
            }
        );
        assert_eq!(name_failure, original);

        let mut list_failure = original.clone();
        let error = list_failure
            .push_member_with_reservations(
                "model-b",
                2.0,
                |value, additional| {
                    value
                        .try_reserve_exact(additional)
                        .map_err(|_| BracketReservationError)
                },
                |_, _| Err(BracketReservationError),
            )
            .expect_err("injected member-list reservation failure");
        assert_eq!(
            error,
            FitError::BracketMemberListAllocationFailed {
                requested_members: 2,
            }
        );
        assert_eq!(list_failure, original);
    }

    #[test]
    fn finite_extreme_bracket_keeps_finite_numerical_evidence() {
        let evidence = ModelBracket::new()
            .try_with_member("negative-extreme", -f64::MAX)
            .expect("first member")
            .try_with_member("positive-extreme", f64::MAX)
            .expect("second member")
            .evidence(ProvenanceHash(0))
            .expect("finite bracket");
        assert!(evidence.numerical.lo.is_finite());
        assert!(evidence.numerical.hi.is_finite());
        assert!(evidence.model.discrepancy_rel.is_infinite());
        // Bracket evidence is no longer certifiable off a hard-coded in-domain
        // literal (bead f85xj.2.10); the numerical slice stays finite either
        // way, which is what this test is about.
        assert!(matches!(
            evidence.clone().certified(),
            Err(crate::CertifyError::OutOfDomain)
        ));
        assert!(evidence.breakdown().model_rel.is_infinite());
    }

    /// Regression, bead frankensim-extreal-program-f85xj.2.10 — a bracket holds
    /// only `(name, qoi)` rows, so `evidence()` may not assert that any member
    /// was evaluated inside its validity range. Before the fix `in_domain` was
    /// the literal `true` and that literal was exactly what let the evidence
    /// certify.
    #[test]
    fn bracket_evidence_does_not_forge_an_in_domain_verdict_f85xj_2_10() {
        let ev = bracket(&[("contact-angle-60", 0.90), ("contact-angle-120", 1.16)])
            .evidence(ProvenanceHash::of_bytes(b"vessel-lip"))
            .expect("valid bracket");
        assert!(
            !ev.model.in_domain,
            "a bracket cannot certify domain membership it never checked"
        );
        assert!(ev.model.validity.bounds().is_empty());
        assert!(
            ev.model
                .assumptions
                .iter()
                .any(|a| a.contains("model-CHOICE sensitivity")),
            "{:?}",
            ev.model.assumptions
        );
        assert!(matches!(
            ev.clone().certified(),
            Err(crate::CertifyError::OutOfDomain)
        ));
        assert!(ev.breakdown().model_rel.is_infinite());

        // The honest door: supply the domain and the evaluation point and the
        // verdict is COMPUTED, not asserted.
        let validity = ValidityDomain::unconstrained().with("Re", 1e4, 1e6);
        let inside = bracket(&[("contact-angle-60", 0.90), ("contact-angle-120", 1.16)])
            .evidence_at(
                ProvenanceHash::of_bytes(b"vessel-lip"),
                validity.clone(),
                &pt(&[("Re", 2e5)]),
            )
            .expect("valid bracket");
        assert!(inside.model.in_domain);
        assert!(inside.clone().certified().is_ok());
        let outside = bracket(&[("contact-angle-60", 0.90), ("contact-angle-120", 1.16)])
            .evidence_at(
                ProvenanceHash::of_bytes(b"vessel-lip"),
                validity,
                &pt(&[("Re", 2e7)]),
            )
            .expect("valid bracket");
        assert!(!outside.model.in_domain);
        assert!(matches!(
            outside.certified(),
            Err(crate::CertifyError::OutOfDomain)
        ));
    }

    /// Regression, bead frankensim-extreal-program-f85xj.2.7 — an axis-aligned
    /// bounding box is not an in-distribution test. Before the fix a box corner
    /// where no pair was ever evaluated returned
    /// `ModelEvidence{discrepancy_rel: 0.0, in_domain: true}`, which certified.
    #[test]
    fn discrepancy_evidence_refuses_unsupported_box_corners_f85xj_2_7() {
        let model = DiscrepancyModel::fit(&[
            FidelityPair {
                params: pt(&[("Re", 1e4), ("Ma", 0.1)]),
                lo_fi: 1.0,
                hi_fi: 1.0,
            },
            FidelityPair {
                params: pt(&[("Re", 1e6), ("Ma", 0.9)]),
                lo_fi: 1.0,
                hi_fi: 1.0,
            },
        ])
        .expect("two-corner fit");
        let corner = pt(&[("Re", 1e6), ("Ma", 0.1)]);

        // The box test still admits the corner: it IS inside the hull.
        let band = model.query(&corner).expect("inside the observed box");
        assert_eq!(band.max_observed_rel.to_bits(), 0.0_f64.to_bits());

        // Minting model evidence there does not.
        let refused = model
            .evidence_at("panel-vs-les", &corner)
            .expect_err("an unvisited box corner is not supported evidence");
        assert!(matches!(
            refused,
            FitError::QueryBeyondFillDistance {
                distance,
                fill_distance: 0.0,
            } if (distance - 1.0).abs() < 1e-12
        ));
        assert!(refused.to_string().contains("extrapolation"), "{refused}");

        // A declared fill distance is an explicit, recorded assumption — and it
        // still has to cover the query.
        let narrow = model
            .clone()
            .with_declared_fill_distance(0.25)
            .expect("valid radius");
        assert!(matches!(
            narrow.evidence_at("panel-vs-les", &corner),
            Err(FitError::QueryBeyondFillDistance { .. })
        ));
        let wide = model
            .clone()
            .with_declared_fill_distance(1.0)
            .expect("valid radius");
        let evidence = wide
            .evidence_at("panel-vs-les", &corner)
            .expect("declared support covers the corner");
        assert!(evidence.in_domain);
        assert!(
            evidence
                .assumptions
                .iter()
                .any(|a| a.contains("declared fill distance")),
            "{:?}",
            evidence.assumptions
        );
        for bad in [0.0, -0.5, 1.5, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                model.clone().with_declared_fill_distance(bad),
                Err(FitError::InvalidFillDistance { .. })
            ));
        }

        // An exact training point needs no declaration at all.
        assert!(
            model
                .evidence_at("panel-vs-les", &pt(&[("Re", 1e4), ("Ma", 0.1)]))
                .is_ok()
        );
    }
}
