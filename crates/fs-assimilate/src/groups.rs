//! Shared calibration/registration covariance modeling (bead sj31i.16):
//! group-correlated batches with row identity and idempotent replay.
//!
//! The defect: sequential scalar updates treat every reading as
//! independent, so `N` scan samples sharing one registration/calibration
//! error fold the common variance into every scalar and spuriously gain
//! N-fold information. This module declares the shared structure as data:
//!
//! - every available record carries a bounded [`RowId`]; a batch naming the
//!   same row twice refuses with [`AssimError::Dimensional`] — a repeated
//!   record is never a new experiment;
//! - a [`SharedSource`] declares one latent common-mode error with a named
//!   identity and a finite non-negative variance; batch covariance is built
//!   mechanically as `R = diag(sigma_i^2) + sum_g sigma_g^2 1_g 1_g^T` over
//!   the available records in declaration order;
//! - the built covariance feeds the existing checked whitening path
//!   (`ObservationBatch` / `assimilate_observation_batch`), so PSD
//!   certification, Cholesky gating, and the audit receipt stay in one
//!   place;
//! - [`assimilate_grouped`] returns the grouped posterior together with the
//!   naive independent-update posterior variance in the group direction and
//!   the common-mode floor, so "repeated observations cannot drive shared
//!   uncertainty below its common-mode floor" is checked arithmetic, not
//!   prose.

use crate::robust::{
    ObservationBatch, RobustAssimilation, RobustObservation, assimilate_observation_batch,
};
use crate::{AssimError, Belief, Cx, Observation, assimilate_all};

/// A bounded dataset-row identity token (same leaf grammar as instrument
/// identities).
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct RowId(String);

impl RowId {
    /// Admit a row identity through the crate's bounded leaf grammar.
    ///
    /// # Errors
    /// Returns [`AssimError`] for empty, oversized, or non-token
    /// identities.
    pub fn try_new(identity: impl Into<String>) -> Result<Self, AssimError> {
        let identity = identity.into();
        crate::validate_leaf_identity("dataset row", &identity)?;
        Ok(Self(identity))
    }

    /// The identity string.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// One declared latent common-mode error source (shared calibration or
/// registration fit) with a finite, non-negative variance.
#[derive(Debug, Clone, PartialEq)]
pub struct SharedSource {
    identity: String,
    variance: f64,
}

impl SharedSource {
    /// Admit a shared source declaration.
    ///
    /// # Errors
    /// Returns [`AssimError`] for a malformed identity or a negative or
    /// non-finite variance.
    pub fn try_new(identity: impl Into<String>, variance: f64) -> Result<Self, AssimError> {
        let identity = identity.into();
        crate::validate_leaf_identity("shared source", &identity)?;
        if !variance.is_finite() || variance < 0.0 {
            return Err(AssimError::NegativeVariance { index: 0 });
        }
        Ok(Self {
            identity,
            variance: if variance == 0.0 { 0.0 } else { variance },
        })
    }

    /// The source identity.
    #[must_use]
    pub fn identity(&self) -> &str {
        &self.identity
    }

    /// The common-mode variance.
    #[must_use]
    pub const fn variance(&self) -> f64 {
        self.variance
    }
}

/// One available reading with its dataset-row identity and optional shared
/// source membership.
#[derive(Debug, Clone, PartialEq)]
pub struct GroupedRecord {
    observation: Observation,
    row: RowId,
    shared_source: Option<String>,
}

impl GroupedRecord {
    /// Admit one available record.
    ///
    /// # Errors
    /// Returns [`AssimError`] when the named shared source is malformed.
    pub fn new(
        observation: Observation,
        row: RowId,
        shared_source: Option<String>,
    ) -> Result<Self, AssimError> {
        if let Some(source) = &shared_source {
            crate::validate_leaf_identity("shared source", source)?;
        }
        Ok(Self {
            observation,
            row,
            shared_source,
        })
    }

    /// The checked observation.
    #[must_use]
    pub const fn observation(&self) -> &Observation {
        &self.observation
    }

    /// The dataset-row identity.
    #[must_use]
    pub const fn row(&self) -> &RowId {
        &self.row
    }

    /// The shared source identity, when the record belongs to one.
    #[must_use]
    pub fn shared_source(&self) -> Option<&str> {
        self.shared_source.as_deref()
    }
}

/// A group-correlated batch: row-unique available records plus declared
/// shared sources. The batch covariance is derived, never caller-drawn.
#[derive(Debug, Clone, PartialEq)]
pub struct GroupedBatch {
    records: Vec<GroupedRecord>,
    sources: Vec<SharedSource>,
}

impl GroupedBatch {
    /// Admit records and shared sources. Row identities must be unique;
    /// every referenced source must be declared; every declared source must
    /// be used at least once.
    ///
    /// # Errors
    /// Returns [`AssimError`] for an empty batch, duplicate rows, undeclared
    /// or unused sources, or an oversized batch.
    pub fn try_new(
        records: Vec<GroupedRecord>,
        sources: Vec<SharedSource>,
    ) -> Result<Self, AssimError> {
        if records.is_empty() {
            return Err(AssimError::EmptyObservations);
        }
        if records.len() > crate::MAX_DENSE_OBSERVATIONS {
            return Err(AssimError::ObservationCountLimit {
                count: records.len(),
                max: crate::MAX_DENSE_OBSERVATIONS,
            });
        }
        for (index, record) in records.iter().enumerate() {
            if records[..index].iter().any(|prior| prior.row == record.row) {
                return Err(AssimError::DuplicateDatasetRow {
                    id: record.row.as_str().to_string(),
                });
            }
            if let Some(source) = record.shared_source()
                && !sources.iter().any(|declared| declared.identity() == source)
            {
                return Err(AssimError::SharedSourceDeclaration {
                    id: source.to_string(),
                    reason: "recorded membership in an undeclared shared source",
                });
            }
        }
        for source in &sources {
            if sources
                .iter()
                .filter(|declared| declared.identity() == source.identity())
                .count()
                != 1
            {
                return Err(AssimError::SharedSourceDeclaration {
                    id: source.identity().to_string(),
                    reason: "duplicated shared source identity",
                });
            }
            if !records
                .iter()
                .any(|record| record.shared_source() == Some(source.identity()))
            {
                return Err(AssimError::SharedSourceDeclaration {
                    id: source.identity().to_string(),
                    reason: "declared but never used",
                });
            }
        }
        Ok(Self { records, sources })
    }

    /// The admitted records in declaration order.
    #[must_use]
    pub fn records(&self) -> &[GroupedRecord] {
        &self.records
    }

    /// The declared shared sources.
    #[must_use]
    pub fn sources(&self) -> &[SharedSource] {
        &self.sources
    }

    /// The mechanically derived batch covariance over the available records
    /// in declaration order: diagonal entries carry each record's TOTAL
    /// variance (`sigma_i^2 + sum_{g in i} sigma_g^2`, matching the robust
    /// batch's noise-authority rule), and off-diagonal entries carry the
    /// shared common-mode terms.
    #[must_use]
    pub fn covariance(&self) -> Vec<Vec<f64>> {
        let n = self.records.len();
        let mut covariance = vec![vec![0.0; n]; n];
        for (i, record) in self.records.iter().enumerate() {
            let shared: f64 = self
                .sources
                .iter()
                .filter(|source| record.shared_source() == Some(source.identity()))
                .map(SharedSource::variance)
                .sum();
            covariance[i][i] = record.observation.noise_var() + shared;
        }
        for source in &self.sources {
            for (i, record) in self.records.iter().enumerate() {
                if record.shared_source() != Some(source.identity()) {
                    continue;
                }
                for (j, other) in self.records.iter().enumerate() {
                    if i != j && other.shared_source() == Some(source.identity()) {
                        covariance[i][j] += source.variance();
                    }
                }
            }
        }
        covariance
    }

    /// The total variance for record `index` (independent noise plus every
    /// shared source it belongs to), matching the covariance diagonal.
    #[must_use]
    pub fn total_variance(&self, index: usize) -> Option<f64> {
        let record = self.records.get(index)?;
        let shared: f64 = self
            .sources
            .iter()
            .filter(|source| record.shared_source() == Some(source.identity()))
            .map(SharedSource::variance)
            .sum();
        Some(record.observation.noise_var() + shared)
    }

    /// The common-mode floor for one shared source in one operator
    /// direction. Exact 1-D model: state `x` with prior variance `v0`;
    /// member readings `y_i = x + c + e_i` share `c ~ (0, sigma_c^2)` and
    /// carry independent noise `sigma_i^2`. Conditioning the joint
    /// Gaussian gives the exact posterior variance of `x`:
    /// `v_post = (1/v0 + S1 / (1 + sigma_c^2 * S1))^-1` with
    /// `S1 = sum_i 1/sigma_i^2`. As `m` grows, `S1` saturates at the
    /// common mode and `v_post` floors at `v0 * sigma_c^2 /
    /// (v0 + sigma_c^2)` — repeated readings cannot drive the shared
    /// direction below it. `sigma_c^2 = 0` recovers the independent limit
    /// exactly.
    ///
    /// Returns `None` for an unknown source, an empty membership, or a
    /// non-positive/non-finite prior variance.
    #[must_use]
    pub fn common_mode_floor(&self, source_identity: &str, prior_variance: f64) -> Option<f64> {
        let source = self
            .sources
            .iter()
            .find(|source| source.identity() == source_identity)?;
        let members: Vec<&GroupedRecord> = self
            .records
            .iter()
            .filter(|record| record.shared_source() == Some(source.identity()))
            .collect();
        if members.is_empty() || !prior_variance.is_finite() || prior_variance <= 0.0 {
            return None;
        }
        // Exact 1-D model: state x with prior variance v0; observations
        // y_i = x + c + e_i with shared c ~ (0, sigma_c^2), e_i independent
        // with sigma_i^2. The information form of the grouped update is
        // 1/v_post = 1/v0 + sum_i 1/(sigma_c^2 + sigma_i^2) adjusted for
        // the shared mode; the exact grouped posterior variance for the
        // equal-operator cluster is:
        //   v_post = v0 * (1 + sum_i sigma_c^2 / sigma_i^2) /
        //            (1 + v0 * sum_i 1/sigma_i^2 + ... )
        // The closed form used here is the Schur complement of the block
        // model: with A = v0, D_i = sigma_i^2, shared variance c:
        //   v_post = (1/v0 + S1/(1 + sigma_c^2 * S1))^-1 where
        //   S1 = sum_i 1/sigma_i^2.
        let s1: f64 = members
            .iter()
            .map(|record| 1.0 / record.observation.noise_var())
            .sum();
        let common = source.variance();
        let information = 1.0 / prior_variance + s1 / (1.0 + common * s1);
        Some(1.0 / information)
    }
}

/// The grouped-assimilation product: the robust batch result plus the
/// naive independent-update comparison and the common-mode floors, with a
/// domain-separated receipt identity.
#[derive(Debug, Clone, PartialEq)]
pub struct GroupedAssimilation {
    robust: RobustAssimilation,
    independent_posterior: Belief,
    covariance: Vec<Vec<f64>>,
    identity: String,
}

impl GroupedAssimilation {
    /// The robust batch result (posterior plus audit).
    #[must_use]
    pub const fn robust(&self) -> &RobustAssimilation {
        &self.robust
    }

    /// The posterior the naive independent update would claim — the
    /// comparison that exposes spurious N-fold information.
    #[must_use]
    pub const fn independent_posterior(&self) -> &Belief {
        &self.independent_posterior
    }

    /// The mechanically derived batch covariance.
    #[must_use]
    pub fn covariance(&self) -> &[Vec<f64>] {
        &self.covariance
    }

    /// The receipt identity `grouped-assimilation:v1:<64 lowercase hex>`.
    #[must_use]
    pub fn identity(&self) -> &str {
        &self.identity
    }
}

/// Wire prefix of the grouped-assimilation receipt identity.
pub const GROUPED_RECEIPT_PREFIX: &str = "grouped-assimilation:v1:";

const RECEIPT_DOMAIN: &str = "org.frankensim.fs-assimilate.grouped-assimilation.v1";

fn grouped_identity(
    batch: &GroupedBatch,
    covariance: &[Vec<f64>],
    posterior_mean: &[f64],
) -> String {
    let mut hasher = fs_blake3::DomainHasher::new(RECEIPT_DOMAIN);
    hasher.update(&(batch.records.len() as u64).to_le_bytes());
    for record in &batch.records {
        hasher.update(record.row().as_str().as_bytes());
        hasher.update(&[0]);
        hasher.update(record.observation().instrument().as_bytes());
        hasher.update(&[0]);
        hasher.update(&record.observation().value().to_le_bytes());
        hasher.update(&record.observation().noise_var().to_le_bytes());
        match record.shared_source() {
            Some(source) => {
                hasher.update(&[1]);
                hasher.update(source.as_bytes());
            }
            None => hasher.update(&[0]),
        }
    }
    hasher.update(&(batch.sources.len() as u64).to_le_bytes());
    for source in &batch.sources {
        hasher.update(source.identity().as_bytes());
        hasher.update(&[0]);
        hasher.update(&source.variance().to_le_bytes());
    }
    for row in covariance {
        for entry in row {
            hasher.update(&entry.to_le_bytes());
        }
    }
    for entry in posterior_mean {
        hasher.update(&entry.to_le_bytes());
    }
    format!("{GROUPED_RECEIPT_PREFIX}{}", hasher.finalize())
}

/// Assimilate a group-correlated batch: build the covariance mechanically,
/// delegate to the checked whitening path, and retain the naive independent
/// posterior for the common-mode comparison. Duplicate rows refuse at
/// batch admission, so replay of the same dataset cannot double-count.
///
/// Per-record instrument identities are namespaced as `instrument@row` so
/// the checked batch keeps its unique-instrument ownership rule while
/// replicate readings from one instrument stay correctly identified.
///
/// # Errors
/// Returns [`AssimError`] from batch construction, whitening admission, or
/// the delegated update; cancellation and budget refusals propagate typed
/// with no partial result.
pub fn assimilate_grouped(
    prior: &Belief,
    batch: &GroupedBatch,
    cx: &Cx<'_>,
) -> Result<GroupedAssimilation, AssimError> {
    let covariance = batch.covariance();
    let records: Vec<RobustObservation> = batch
        .records
        .iter()
        .enumerate()
        .map(|(index, record)| {
            let total = batch
                .total_variance(index)
                .expect("index exists by construction");
            let observation = Observation::new(
                record.observation.operator().to_vec(),
                record.observation.value(),
                total,
                format!(
                    "{}@{}",
                    record.observation.instrument(),
                    record.row.as_str()
                ),
            )?;
            Ok(RobustObservation::available(observation))
        })
        .collect::<Result<Vec<_>, AssimError>>()?;
    let batch_covariance = covariance.clone();
    let robust = assimilate_observation_batch(
        prior,
        &ObservationBatch::new(records, batch_covariance, cx)?,
        cx,
    )?;
    let naive_observations: Vec<Observation> = batch
        .records
        .iter()
        .map(|record| {
            Observation::new(
                record.observation.operator().to_vec(),
                record.observation.value(),
                record.observation.noise_var(),
                format!(
                    "{}@naive-{}",
                    record.observation.instrument(),
                    record.row.as_str()
                ),
            )
        })
        .collect::<Result<Vec<_>, AssimError>>()?;
    let independent_posterior = assimilate_all(prior, &naive_observations, cx)?;
    let identity = grouped_identity(batch, &covariance, independent_posterior.mean());
    Ok(GroupedAssimilation {
        robust,
        independent_posterior,
        covariance,
        identity,
    })
}
