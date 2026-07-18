//! Closed-form tolerance allocation for admitted perfect-dependence groups.
//!
//! Members of one group share a perfectly positively correlated standardized
//! manufacturing shock. Different groups are mutually uncorrelated at second
//! order. Sensitivities are strictly positive magnitudes, so there is no
//! within-group cancellation. This bounded model admits a unique convex
//! optimum without an iterative solver or arbitrary callbacks.

use std::{collections::BTreeMap, num::NonZeroU64};

use fs_math::det;

use super::{Action, Feature, ScalarIssue, action_for};

/// Fixed receipt and evaluation schema for grouped allocation.
pub const GROUPED_ALLOCATION_SCHEMA_V1: u32 = 1;

/// Maximum number of perfect-dependence groups in one model.
pub const MAX_DEPENDENCY_GROUPS_V1: usize = 128;

/// Maximum number of grouped features in one model.
pub const MAX_GROUPED_FEATURES_V1: usize = 128;

/// Maximum byte length of one group or feature name.
pub const MAX_DEPENDENCY_NAME_BYTES_V1: usize = 256;

/// Maximum byte length of one external model namespace.
pub const MAX_DEPENDENCY_NAMESPACE_BYTES_V1: usize = 256;

/// Stable ordinal identity of one perfect-dependence group.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DependencyGroupId(
    /// Zero-based group ordinal in the retained model.
    pub u16,
);

impl DependencyGroupId {
    /// Zero-based group index.
    #[must_use]
    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

/// One declared perfect-positive-dependence group.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DependencyGroup {
    /// Stable lowercase ASCII group key.
    pub key: String,
}

/// One legacy feature bound to an explicit dependency group.
#[derive(Debug, Clone, PartialEq)]
pub struct GroupedFeature {
    /// Existing feature semantics: sensitivity, color, cost, and baseline.
    pub feature: Feature,
    /// Perfect-dependence group used by this feature.
    pub group: DependencyGroupId,
}

/// Caller-supplied external identity retained with grouped allocation.
///
/// Admission checks grammar and a nonzero digest. It neither computes the
/// digest from model bytes nor authenticates the external owner.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GroupedDependenceIdentity {
    namespace: Box<str>,
    schema_version: NonZeroU64,
    semantic_digest: [u8; 32],
}

impl GroupedDependenceIdentity {
    /// Admit one bounded slash-separated external identity.
    ///
    /// # Errors
    ///
    /// Refuses a noncanonical namespace or all-zero digest.
    pub fn try_new(
        namespace: impl Into<String>,
        schema_version: NonZeroU64,
        semantic_digest: [u8; 32],
    ) -> Result<Self, GroupedAllocationError> {
        let namespace = namespace.into();
        validate_namespace(&namespace)?;
        if semantic_digest == [0; 32] {
            return Err(GroupedAllocationError::ZeroSemanticDigest);
        }
        Ok(Self {
            namespace: namespace.into_boxed_str(),
            schema_version,
            semantic_digest,
        })
    }

    /// External model namespace.
    #[must_use]
    pub fn namespace(&self) -> &str {
        &self.namespace
    }

    /// Explicit external schema version.
    #[must_use]
    pub const fn schema_version(&self) -> NonZeroU64 {
        self.schema_version
    }

    /// Exact caller-supplied semantic digest.
    #[must_use]
    pub const fn semantic_digest(&self) -> [u8; 32] {
        self.semantic_digest
    }
}

/// Complete grouped-dependence allocation input.
#[derive(Debug, Clone, PartialEq)]
pub struct GroupedDependenceModel {
    /// Caller-supplied identity retained alongside the model.
    pub identity: GroupedDependenceIdentity,
    /// Mutually uncorrelated groups in stable ordinal order.
    pub groups: Vec<DependencyGroup>,
    /// Features in caller order, each assigned to exactly one group.
    pub features: Vec<GroupedFeature>,
}

/// Bounded resource named by a grouped-allocation refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GroupedAllocationResource {
    /// Perfect-dependence groups.
    Groups,
    /// Grouped features.
    Features,
}

/// Derived grouped-allocation quantity that failed safe representation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GroupedDerivedQuantity {
    /// One max-shifted log-sum-exp contribution.
    LogSumExpContribution,
    /// One group shape log J_g.
    GroupLogShape,
    /// Global shape log D.
    GlobalLogShape,
    /// Common optimum normalization log alpha.
    AllocationNormalization,
    /// One allocated tolerance.
    Tolerance,
    /// One feature manufacturing-cost contribution.
    CostContribution,
    /// One feature standard-deviation loading s_i t_i / k.
    StandardDeviationLoading,
    /// One feature counterfactual independent variance contribution.
    IndependentVarianceContribution,
    /// One group standard deviation.
    GroupStandardDeviation,
    /// One group coherent variance.
    GroupVariance,
    /// One group counterfactual independent variance.
    GroupIndependentVariance,
    /// One group manufacturing cost.
    GroupCost,
    /// Accumulated manufacturing cost.
    TotalCost,
    /// Closed-form optimum cost oracle.
    ClosedFormCost,
    /// Accumulated coherent variance.
    AchievedVariance,
    /// Accumulated counterfactual independent variance.
    IndependentVariance,
    /// Signed coherent-minus-independent variance delta.
    DependencyVarianceDelta,
    /// Maximum log-domain KKT stationarity residual.
    StationarityResidual,
}

/// Refusal from grouped-model admission or closed-form allocation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GroupedAllocationError {
    /// The external namespace is not a bounded canonical slash-separated key.
    InvalidNamespace {
        /// Bounded prefix of the rejected namespace.
        namespace: String,
        /// Stable grammar explanation.
        reason: &'static str,
    },
    /// An all-zero digest cannot identify an external model.
    ZeroSemanticDigest,
    /// A bounded resource is empty or exceeds its schema cap.
    ResourceLimit {
        /// Resource class.
        resource: GroupedAllocationResource,
        /// Supplied count.
        actual: usize,
        /// Inclusive schema maximum.
        max: usize,
    },
    /// One group key is empty, oversized, or noncanonical.
    InvalidGroupKey {
        /// Group ordinal.
        index: usize,
        /// Bounded prefix of the rejected key.
        key: String,
        /// Stable grammar explanation.
        reason: &'static str,
    },
    /// Two group keys collide exactly in the canonical lowercase grammar.
    DuplicateGroupKey {
        /// First group ordinal.
        first_index: usize,
        /// Colliding group ordinal.
        duplicate_index: usize,
        /// Exact canonical key.
        key: String,
    },
    /// One feature name is empty, oversized, or unstable.
    InvalidFeatureName {
        /// Feature ordinal.
        index: usize,
        /// Bounded prefix of the rejected name.
        name: String,
        /// Stable explanation.
        reason: &'static str,
    },
    /// Two feature names collide under deterministic lowercase comparison.
    DuplicateFeatureName {
        /// First feature ordinal.
        first_index: usize,
        /// Colliding feature ordinal.
        duplicate_index: usize,
        /// Canonical lowercase key.
        canonical_name: String,
    },
    /// One feature references a group outside the supplied group table.
    InvalidGroupReference {
        /// Feature ordinal.
        feature_index: usize,
        /// Supplied group ordinal.
        group: DependencyGroupId,
        /// Number of available groups.
        available: usize,
    },
    /// A declared group has no feature and therefore no allocation semantics.
    EmptyGroup {
        /// Empty group ordinal.
        group: DependencyGroupId,
        /// Stable group key.
        key: String,
    },
    /// One feature scalar is outside its strictly positive finite domain.
    InvalidFeatureField {
        /// Feature ordinal.
        index: usize,
        /// Feature name.
        feature: String,
        /// Rejected field.
        field: &'static str,
        /// Stable numeric issue.
        issue: ScalarIssue,
    },
    /// A public scalar argument is not strictly positive and finite.
    InvalidArgument {
        /// Argument name.
        argument: &'static str,
        /// Stable numeric issue.
        issue: ScalarIssue,
    },
    /// Finite admitted inputs produced an unsafe derived result.
    InvalidDerived {
        /// Failed quantity.
        quantity: GroupedDerivedQuantity,
        /// Feature ordinal when feature-local.
        feature_index: Option<usize>,
        /// Group ordinal when group-local.
        group: Option<DependencyGroupId>,
        /// Stable numeric issue.
        issue: ScalarIssue,
    },
}

/// Privately constructed allocation for one grouped feature.
#[derive(Debug, Clone, PartialEq)]
pub struct GroupedAllocationItem {
    feature_index: usize,
    group: DependencyGroupId,
    tolerance: f64,
    action: Action,
    cost_contribution: f64,
    standard_deviation_loading: f64,
    independent_variance_contribution: f64,
}

impl GroupedAllocationItem {
    /// Feature ordinal in the retained model.
    #[must_use]
    pub const fn feature_index(&self) -> usize {
        self.feature_index
    }

    /// Declared dependency group.
    #[must_use]
    pub const fn group(&self) -> DependencyGroupId {
        self.group
    }

    /// Cost-optimal allocated tolerance.
    #[must_use]
    pub const fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Tighten, loosen, or unchanged relative to the retained baseline.
    #[must_use]
    pub const fn action(&self) -> Action {
        self.action
    }

    /// Manufacturing-cost contribution c_i / t_i.
    #[must_use]
    pub const fn cost_contribution(&self) -> f64 {
        self.cost_contribution
    }

    /// Positive coherent loading s_i t_i / k.
    #[must_use]
    pub const fn standard_deviation_loading(&self) -> f64 {
        self.standard_deviation_loading
    }

    /// Counterfactual variance if this feature shock were independent.
    #[must_use]
    pub const fn independent_variance_contribution(&self) -> f64 {
        self.independent_variance_contribution
    }
}

/// Retained moments and cost for one perfect-dependence group.
#[derive(Debug, Clone, PartialEq)]
pub struct DependencyGroupReceipt {
    group: DependencyGroupId,
    feature_count: usize,
    log_shape: f64,
    standard_deviation: f64,
    variance: f64,
    independent_variance: f64,
    dependency_variance_delta: f64,
    total_cost: f64,
}

impl DependencyGroupReceipt {
    /// Group identity.
    #[must_use]
    pub const fn group(&self) -> DependencyGroupId {
        self.group
    }

    /// Number of features assigned to this group.
    #[must_use]
    pub const fn feature_count(&self) -> usize {
        self.feature_count
    }

    /// Scale-safe log J_g, where J_g is the sum of square-root cost-load terms.
    #[must_use]
    pub const fn log_shape(&self) -> f64 {
        self.log_shape
    }

    /// Coherent group standard deviation, the sum of member loadings.
    #[must_use]
    pub const fn standard_deviation(&self) -> f64 {
        self.standard_deviation
    }

    /// Coherent group variance.
    #[must_use]
    pub const fn variance(&self) -> f64 {
        self.variance
    }

    /// Counterfactual variance under independent member shocks.
    #[must_use]
    pub const fn independent_variance(&self) -> f64 {
        self.independent_variance
    }

    /// Signed coherent variance minus counterfactual independent variance.
    #[must_use]
    pub const fn dependency_variance_delta(&self) -> f64 {
        self.dependency_variance_delta
    }

    /// Manufacturing cost accumulated for this group's features.
    #[must_use]
    pub const fn total_cost(&self) -> f64 {
        self.total_cost
    }
}

/// Sealed result of one admitted perfect-dependence allocation.
#[derive(Debug, Clone, PartialEq)]
pub struct GroupedAllocationReceipt {
    schema_version: u32,
    model: GroupedDependenceModel,
    variance_budget: f64,
    tolerance_to_sigma: f64,
    items: Box<[GroupedAllocationItem]>,
    groups: Box<[DependencyGroupReceipt]>,
    total_cost: f64,
    closed_form_cost: f64,
    cost_residual: f64,
    log_scale_correction: f64,
    achieved_variance: f64,
    budget_residual: f64,
    independent_variance: f64,
    dependency_variance_delta: f64,
    max_stationarity_log_residual: f64,
}

impl GroupedAllocationReceipt {
    /// Fixed grouped-allocation receipt/evaluation schema.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Exact retained input model.
    #[must_use]
    pub const fn model(&self) -> &GroupedDependenceModel {
        &self.model
    }

    /// Caller-supplied coherent variance budget.
    #[must_use]
    pub const fn variance_budget(&self) -> f64 {
        self.variance_budget
    }

    /// Caller-supplied tolerance-to-standard-deviation factor k.
    #[must_use]
    pub const fn tolerance_to_sigma(&self) -> f64 {
        self.tolerance_to_sigma
    }

    /// Allocated items in caller feature order.
    #[must_use]
    pub fn items(&self) -> &[GroupedAllocationItem] {
        &self.items
    }

    /// Per-group receipts in declared group order.
    #[must_use]
    pub fn groups(&self) -> &[DependencyGroupReceipt] {
        &self.groups
    }

    /// Cost accumulated from published feature tolerances.
    #[must_use]
    pub const fn total_cost(&self) -> f64 {
        self.total_cost
    }

    /// Separately evaluated grouped closed-form optimum cost oracle.
    #[must_use]
    pub const fn closed_form_cost(&self) -> f64 {
        self.closed_form_cost
    }

    /// Accumulated cost minus the closed-form oracle.
    #[must_use]
    pub const fn cost_residual(&self) -> f64 {
        self.cost_residual
    }

    /// Common log-tolerance correction applied after binary64 normalization.
    #[must_use]
    pub const fn log_scale_correction(&self) -> f64 {
        self.log_scale_correction
    }

    /// Coherent grouped variance accumulated from group receipts.
    #[must_use]
    pub const fn achieved_variance(&self) -> f64 {
        self.achieved_variance
    }

    /// Achieved coherent variance minus the requested budget.
    #[must_use]
    pub const fn budget_residual(&self) -> f64 {
        self.budget_residual
    }

    /// Counterfactual variance if every feature shock were independent.
    #[must_use]
    pub const fn independent_variance(&self) -> f64 {
        self.independent_variance
    }

    /// Signed coherent variance minus counterfactual independent variance.
    #[must_use]
    pub const fn dependency_variance_delta(&self) -> f64 {
        self.dependency_variance_delta
    }

    /// Largest absolute feature-to-feature KKT log-stationarity mismatch.
    #[must_use]
    pub const fn max_stationarity_log_residual(&self) -> f64 {
        self.max_stationarity_log_residual
    }
}

fn bounded_utf8_prefix(value: &str, max_bytes: usize) -> String {
    let mut end = value.len().min(max_bytes);
    while !value.is_char_boundary(end) {
        end -= 1;
    }
    value[..end].to_string()
}

fn is_canonical_key_bytes(key: &[u8]) -> bool {
    key.first().is_some_and(u8::is_ascii_lowercase)
        && key.last() != Some(&b'-')
        && !key.windows(2).any(|pair| pair == b"--")
        && key
            .iter()
            .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || *byte == b'-')
}

fn validate_namespace(namespace: &str) -> Result<(), GroupedAllocationError> {
    let reason = if namespace.is_empty() {
        Some("namespace must not be empty")
    } else if namespace.len() > MAX_DEPENDENCY_NAMESPACE_BYTES_V1 {
        Some("namespace exceeds the versioned byte cap")
    } else if namespace
        .as_bytes()
        .split(|byte| *byte == b'/')
        .any(|segment| segment.is_empty())
    {
        Some("namespace segments must not be empty")
    } else if namespace
        .as_bytes()
        .split(|byte| *byte == b'/')
        .any(|segment| !is_canonical_key_bytes(segment))
    {
        Some(
            "namespace segments must start lowercase and use lowercase ASCII letters, digits, or single interior hyphens",
        )
    } else {
        None
    };
    if let Some(reason) = reason {
        Err(GroupedAllocationError::InvalidNamespace {
            namespace: bounded_utf8_prefix(namespace, MAX_DEPENDENCY_NAMESPACE_BYTES_V1),
            reason,
        })
    } else {
        Ok(())
    }
}

fn validate_group_key(index: usize, key: &str) -> Result<(), GroupedAllocationError> {
    let reason = if key.is_empty() {
        Some("key must not be empty")
    } else if key.len() > MAX_DEPENDENCY_NAME_BYTES_V1 {
        Some("key exceeds the versioned byte cap")
    } else if !key.is_ascii() || !is_canonical_key_bytes(key.as_bytes()) {
        Some(
            "key must start lowercase and use lowercase ASCII letters, digits, or single interior hyphens",
        )
    } else {
        None
    };
    if let Some(reason) = reason {
        Err(GroupedAllocationError::InvalidGroupKey {
            index,
            key: bounded_utf8_prefix(key, MAX_DEPENDENCY_NAME_BYTES_V1),
            reason,
        })
    } else {
        Ok(())
    }
}

fn canonical_feature_name(index: usize, name: &str) -> Result<String, GroupedAllocationError> {
    let reason = if name.is_empty() {
        Some("name must not be empty")
    } else if name.len() > MAX_DEPENDENCY_NAME_BYTES_V1 {
        Some("name exceeds the versioned byte cap")
    } else if name.trim() != name {
        Some("name must not have leading or trailing whitespace")
    } else if name.chars().any(char::is_control) {
        Some("name must not contain control characters")
    } else {
        None
    };
    if let Some(reason) = reason {
        return Err(GroupedAllocationError::InvalidFeatureName {
            index,
            name: bounded_utf8_prefix(name, MAX_DEPENDENCY_NAME_BYTES_V1),
            reason,
        });
    }
    let canonical = name.to_lowercase();
    if canonical.len() > MAX_DEPENDENCY_NAME_BYTES_V1 {
        Err(GroupedAllocationError::InvalidFeatureName {
            index,
            name: name.to_string(),
            reason: "lowercase comparison key exceeds the versioned byte cap",
        })
    } else {
        Ok(canonical)
    }
}

fn validate_positive_argument(
    argument: &'static str,
    value: f64,
) -> Result<(), GroupedAllocationError> {
    let issue = if !value.is_finite() {
        Some(ScalarIssue::NonFinite)
    } else if value <= 0.0 {
        Some(ScalarIssue::NonPositive)
    } else {
        None
    };
    if let Some(issue) = issue {
        Err(GroupedAllocationError::InvalidArgument { argument, issue })
    } else {
        Ok(())
    }
}

fn validate_positive_feature(
    index: usize,
    feature: &Feature,
    field: &'static str,
    value: f64,
) -> Result<(), GroupedAllocationError> {
    let issue = if !value.is_finite() {
        Some(ScalarIssue::NonFinite)
    } else if value <= 0.0 {
        Some(ScalarIssue::NonPositive)
    } else {
        None
    };
    if let Some(issue) = issue {
        Err(GroupedAllocationError::InvalidFeatureField {
            index,
            feature: feature.name.clone(),
            field,
            issue,
        })
    } else {
        Ok(())
    }
}

fn invalid_derived(
    quantity: GroupedDerivedQuantity,
    feature_index: Option<usize>,
    group: Option<DependencyGroupId>,
    issue: ScalarIssue,
) -> GroupedAllocationError {
    GroupedAllocationError::InvalidDerived {
        quantity,
        feature_index,
        group,
        issue,
    }
}

fn validate_positive_derived(
    quantity: GroupedDerivedQuantity,
    feature_index: Option<usize>,
    group: Option<DependencyGroupId>,
    value: f64,
) -> Result<(), GroupedAllocationError> {
    if !value.is_finite() {
        Err(invalid_derived(
            quantity,
            feature_index,
            group,
            ScalarIssue::NonFinite,
        ))
    } else if value <= 0.0 {
        Err(invalid_derived(
            quantity,
            feature_index,
            group,
            if value == 0.0 {
                ScalarIssue::Underflow
            } else {
                ScalarIssue::NonPositive
            },
        ))
    } else {
        Ok(())
    }
}

fn validate_nonnegative_derived(
    quantity: GroupedDerivedQuantity,
    feature_index: Option<usize>,
    group: Option<DependencyGroupId>,
    value: f64,
) -> Result<(), GroupedAllocationError> {
    if !value.is_finite() {
        Err(invalid_derived(
            quantity,
            feature_index,
            group,
            ScalarIssue::NonFinite,
        ))
    } else if value < 0.0 {
        Err(invalid_derived(
            quantity,
            feature_index,
            group,
            ScalarIssue::Negative,
        ))
    } else {
        Ok(())
    }
}

fn validate_dependency_variance_delta(
    group: Option<DependencyGroupId>,
    has_coherent_cross_terms: bool,
    value: f64,
) -> Result<(), GroupedAllocationError> {
    validate_nonnegative_derived(
        GroupedDerivedQuantity::DependencyVarianceDelta,
        None,
        group,
        value,
    )?;
    // Every pair of strictly positive coherent loadings contributes the
    // strictly positive cross term 2 a_i a_j. A multi-member group, or a
    // grouped model containing one, can therefore never have an exact zero
    // coherent-minus-independent delta.
    if has_coherent_cross_terms && value == 0.0 {
        return Err(invalid_derived(
            GroupedDerivedQuantity::DependencyVarianceDelta,
            None,
            group,
            ScalarIssue::Underflow,
        ));
    }
    Ok(())
}

fn compensated_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut sum = 0.0_f64;
    let mut correction = 0.0_f64;
    for value in values {
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}

fn log_sum_exp(
    terms: &[f64],
    quantity: GroupedDerivedQuantity,
    group: Option<DependencyGroupId>,
) -> Result<f64, GroupedAllocationError> {
    if terms.is_empty() || terms.iter().any(|term| !term.is_finite()) {
        return Err(invalid_derived(
            quantity,
            None,
            group,
            ScalarIssue::NonFinite,
        ));
    }
    let maximum = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut scaled_terms = Vec::with_capacity(terms.len());
    for term in terms {
        let scaled = det::exp(*term - maximum);
        if !scaled.is_finite() {
            return Err(invalid_derived(
                GroupedDerivedQuantity::LogSumExpContribution,
                None,
                group,
                ScalarIssue::NonFinite,
            ));
        }
        if scaled == 0.0 {
            return Err(invalid_derived(
                GroupedDerivedQuantity::LogSumExpContribution,
                None,
                group,
                ScalarIssue::Underflow,
            ));
        }
        scaled_terms.push(scaled);
    }
    let scaled_sum = compensated_sum(scaled_terms.iter().copied());
    if !scaled_sum.is_finite() || scaled_sum <= 0.0 {
        return Err(invalid_derived(
            quantity,
            None,
            group,
            ScalarIssue::NonFinite,
        ));
    }
    // Max shifting makes at least one contribution exactly 1.0. With more
    // than one term, the exact scaled sum is therefore strictly greater than
    // 1.0. Refuse if binary64 accumulation erases that entire positive tail:
    // accepting 1.0 would silently reduce a multi-member/group expression to
    // its single dominant term.
    if terms.len() > 1 && scaled_sum == 1.0 {
        return Err(invalid_derived(
            GroupedDerivedQuantity::LogSumExpContribution,
            None,
            group,
            ScalarIssue::Underflow,
        ));
    }
    let result = maximum + det::ln(scaled_sum);
    if !result.is_finite() {
        return Err(invalid_derived(
            quantity,
            None,
            group,
            ScalarIssue::NonFinite,
        ));
    }
    // Even when the scaled tail survives accumulation, adding its logarithm
    // back to a much larger maximum can round to that maximum. The exact LSE
    // of multiple finite terms is strictly greater than its largest term.
    if terms.len() > 1 && result <= maximum {
        return Err(invalid_derived(
            GroupedDerivedQuantity::LogSumExpContribution,
            None,
            group,
            ScalarIssue::Underflow,
        ));
    }
    // Aggregate monotonicity is insufficient: one visible tail term can move
    // the sum while a smaller retained term vanishes behind it. Audit every
    // contribution against the deterministic leave-one-out counterfactual.
    // The full scaled sum and reconstructed LSE must each be strictly greater
    // than their counterpart without that positive term. The quadratic audit
    // is bounded by the schema's 128-term maximum.
    if scaled_terms.len() > 1 {
        for omitted_index in 0..scaled_terms.len() {
            let leave_one_out_sum = compensated_sum(
                scaled_terms
                    .iter()
                    .enumerate()
                    .filter(|(index, _)| *index != omitted_index)
                    .map(|(_, scaled)| *scaled),
            );
            if !leave_one_out_sum.is_finite() {
                return Err(invalid_derived(
                    GroupedDerivedQuantity::LogSumExpContribution,
                    None,
                    group,
                    ScalarIssue::NonFinite,
                ));
            }
            if leave_one_out_sum <= 0.0 || scaled_sum <= leave_one_out_sum {
                return Err(invalid_derived(
                    GroupedDerivedQuantity::LogSumExpContribution,
                    None,
                    group,
                    ScalarIssue::Underflow,
                ));
            }
            // Re-center the reconstructed counterfactual on its own maximum.
            // The common-scale sum above audits accumulation influence; this
            // independent max shift audits the actual leave-one-out LSE
            // without inheriting avoidable precision loss when the omitted
            // term is the unique maximum.
            let leave_one_out_maximum = terms
                .iter()
                .enumerate()
                .filter(|(index, _)| *index != omitted_index)
                .map(|(_, term)| *term)
                .fold(f64::NEG_INFINITY, f64::max);
            let leave_one_out_recentered_sum = compensated_sum(
                terms
                    .iter()
                    .enumerate()
                    .filter(|(index, _)| *index != omitted_index)
                    .map(|(_, term)| det::exp(*term - leave_one_out_maximum)),
            );
            if !leave_one_out_recentered_sum.is_finite() {
                return Err(invalid_derived(
                    GroupedDerivedQuantity::LogSumExpContribution,
                    None,
                    group,
                    ScalarIssue::NonFinite,
                ));
            }
            if leave_one_out_recentered_sum <= 0.0 {
                return Err(invalid_derived(
                    GroupedDerivedQuantity::LogSumExpContribution,
                    None,
                    group,
                    ScalarIssue::Underflow,
                ));
            }
            let leave_one_out_result =
                leave_one_out_maximum + det::ln(leave_one_out_recentered_sum);
            if !leave_one_out_result.is_finite() {
                return Err(invalid_derived(
                    GroupedDerivedQuantity::LogSumExpContribution,
                    None,
                    group,
                    ScalarIssue::NonFinite,
                ));
            }
            if result <= leave_one_out_result {
                return Err(invalid_derived(
                    GroupedDerivedQuantity::LogSumExpContribution,
                    None,
                    group,
                    ScalarIssue::Underflow,
                ));
            }
        }
    }
    Ok(result)
}

fn validate_model(
    model: &GroupedDependenceModel,
) -> Result<Vec<Vec<usize>>, GroupedAllocationError> {
    if model.groups.is_empty() || model.groups.len() > MAX_DEPENDENCY_GROUPS_V1 {
        return Err(GroupedAllocationError::ResourceLimit {
            resource: GroupedAllocationResource::Groups,
            actual: model.groups.len(),
            max: MAX_DEPENDENCY_GROUPS_V1,
        });
    }
    if model.features.is_empty() || model.features.len() > MAX_GROUPED_FEATURES_V1 {
        return Err(GroupedAllocationError::ResourceLimit {
            resource: GroupedAllocationResource::Features,
            actual: model.features.len(),
            max: MAX_GROUPED_FEATURES_V1,
        });
    }

    let mut group_keys = BTreeMap::new();
    for (index, group) in model.groups.iter().enumerate() {
        validate_group_key(index, &group.key)?;
        if let Some(&first_index) = group_keys.get(&group.key) {
            return Err(GroupedAllocationError::DuplicateGroupKey {
                first_index,
                duplicate_index: index,
                key: group.key.clone(),
            });
        }
        group_keys.insert(group.key.clone(), index);
    }

    let mut feature_names = BTreeMap::new();
    let mut group_features = vec![Vec::new(); model.groups.len()];
    for (index, grouped) in model.features.iter().enumerate() {
        let canonical_name = canonical_feature_name(index, &grouped.feature.name)?;
        if let Some(&first_index) = feature_names.get(&canonical_name) {
            return Err(GroupedAllocationError::DuplicateFeatureName {
                first_index,
                duplicate_index: index,
                canonical_name,
            });
        }
        feature_names.insert(canonical_name, index);
        if grouped.group.index() >= model.groups.len() {
            return Err(GroupedAllocationError::InvalidGroupReference {
                feature_index: index,
                group: grouped.group,
                available: model.groups.len(),
            });
        }
        validate_positive_feature(
            index,
            &grouped.feature,
            "sensitivity",
            grouped.feature.sensitivity,
        )?;
        validate_positive_feature(
            index,
            &grouped.feature,
            "cost_coeff",
            grouped.feature.cost_coeff,
        )?;
        validate_positive_feature(
            index,
            &grouped.feature,
            "baseline_tolerance",
            grouped.feature.baseline_tolerance,
        )?;
        group_features[grouped.group.index()].push(index);
    }
    for (index, members) in group_features.iter().enumerate() {
        if members.is_empty() {
            return Err(GroupedAllocationError::EmptyGroup {
                group: DependencyGroupId(index as u16),
                key: model.groups[index].key.clone(),
            });
        }
    }
    Ok(group_features)
}

fn compute_group_log_standard_deviations(
    model: &GroupedDependenceModel,
    group_features: &[Vec<usize>],
    log_tolerances: &[f64],
    log_k: f64,
) -> Result<Vec<f64>, GroupedAllocationError> {
    let mut result = Vec::with_capacity(group_features.len());
    for (group_index, members) in group_features.iter().enumerate() {
        let id = DependencyGroupId(group_index as u16);
        let terms = members
            .iter()
            .map(|index| {
                det::ln(model.features[*index].feature.sensitivity) + log_tolerances[*index] - log_k
            })
            .collect::<Vec<_>>();
        result.push(log_sum_exp(
            &terms,
            GroupedDerivedQuantity::GroupStandardDeviation,
            Some(id),
        )?);
    }
    Ok(result)
}

fn achieved_variance_from_logs(
    group_log_standard_deviations: &[f64],
) -> Result<f64, GroupedAllocationError> {
    let contributions = group_log_standard_deviations
        .iter()
        .map(|log_standard_deviation| det::exp(2.0 * *log_standard_deviation))
        .collect::<Vec<_>>();
    for contribution in &contributions {
        validate_positive_derived(
            GroupedDerivedQuantity::GroupVariance,
            None,
            None,
            *contribution,
        )?;
    }
    let total = compensated_sum(contributions);
    validate_positive_derived(GroupedDerivedQuantity::AchievedVariance, None, None, total)?;
    Ok(total)
}

/// Allocate the unique minimum-cost tolerances for perfect-positive groups.
///
/// For group g, J_g is the sum over its features of sqrt(c_i s_i), D is the
/// sum over groups of J_g^(4/3), and alpha is sqrt(B / D). The unique optimum
/// is t_i = k alpha sqrt(c_i / s_i) J_g^(-1/3). Evaluation uses deterministic
/// max-shifted log-sum-exp and one global scale correction before publishing
/// checked binary64 results.
///
/// # Errors
///
/// Refuses malformed identity, empty or oversized tables, unstable names,
/// unused or invalid groups, nonpositive feature/argument scalars, and any
/// required positive result that overflows, underflows to zero, or otherwise
/// becomes non-finite. No partial receipt is returned.
pub fn allocate_grouped(
    model: &GroupedDependenceModel,
    variance_budget: f64,
    k: f64,
) -> Result<GroupedAllocationReceipt, GroupedAllocationError> {
    validate_positive_argument("variance_budget", variance_budget)?;
    validate_positive_argument("k", k)?;
    let group_features = validate_model(model)?;

    let log_budget = det::ln(variance_budget);
    let log_k = det::ln(k);
    if !log_budget.is_finite() || !log_k.is_finite() {
        return Err(invalid_derived(
            GroupedDerivedQuantity::AllocationNormalization,
            None,
            None,
            ScalarIssue::NonFinite,
        ));
    }

    let mut feature_log_sensitivity = Vec::with_capacity(model.features.len());
    let mut feature_log_cost = Vec::with_capacity(model.features.len());
    let mut half_shape_terms = Vec::with_capacity(model.features.len());
    for grouped in &model.features {
        let log_sensitivity = det::ln(grouped.feature.sensitivity);
        let log_cost = det::ln(grouped.feature.cost_coeff);
        let half_shape = 0.5 * (log_sensitivity + log_cost);
        if !log_sensitivity.is_finite() || !log_cost.is_finite() || !half_shape.is_finite() {
            return Err(invalid_derived(
                GroupedDerivedQuantity::GroupLogShape,
                None,
                Some(grouped.group),
                ScalarIssue::NonFinite,
            ));
        }
        feature_log_sensitivity.push(log_sensitivity);
        feature_log_cost.push(log_cost);
        half_shape_terms.push(half_shape);
    }

    let mut group_log_shapes = Vec::with_capacity(model.groups.len());
    for (group_index, members) in group_features.iter().enumerate() {
        let id = DependencyGroupId(group_index as u16);
        let terms = members
            .iter()
            .map(|index| half_shape_terms[*index])
            .collect::<Vec<_>>();
        group_log_shapes.push(log_sum_exp(
            &terms,
            GroupedDerivedQuantity::GroupLogShape,
            Some(id),
        )?);
    }
    let global_terms = group_log_shapes
        .iter()
        .map(|log_shape| (4.0 / 3.0) * *log_shape)
        .collect::<Vec<_>>();
    let global_log_shape =
        log_sum_exp(&global_terms, GroupedDerivedQuantity::GlobalLogShape, None)?;
    let log_alpha = 0.5 * (log_budget - global_log_shape);
    if !log_alpha.is_finite() {
        return Err(invalid_derived(
            GroupedDerivedQuantity::AllocationNormalization,
            None,
            None,
            ScalarIssue::NonFinite,
        ));
    }

    let mut log_tolerances = Vec::with_capacity(model.features.len());
    for (index, grouped) in model.features.iter().enumerate() {
        let log_tolerance =
            log_k + log_alpha + 0.5 * (feature_log_cost[index] - feature_log_sensitivity[index])
                - group_log_shapes[grouped.group.index()] / 3.0;
        if !log_tolerance.is_finite() {
            return Err(invalid_derived(
                GroupedDerivedQuantity::Tolerance,
                Some(index),
                Some(grouped.group),
                ScalarIssue::NonFinite,
            ));
        }
        log_tolerances.push(log_tolerance);
    }

    // Correct the one common scale from the actual binary64 grouped variance.
    // This preserves all KKT ratios and tightens the published budget residual.
    let initial_group_logs =
        compute_group_log_standard_deviations(model, &group_features, &log_tolerances, log_k)?;
    let initial_variance = achieved_variance_from_logs(&initial_group_logs)?;
    let correction = 0.5 * (log_budget - det::ln(initial_variance));
    if !correction.is_finite() {
        return Err(invalid_derived(
            GroupedDerivedQuantity::AllocationNormalization,
            None,
            None,
            ScalarIssue::NonFinite,
        ));
    }
    for log_tolerance in &mut log_tolerances {
        *log_tolerance += correction;
        if !log_tolerance.is_finite() {
            return Err(invalid_derived(
                GroupedDerivedQuantity::Tolerance,
                None,
                None,
                ScalarIssue::NonFinite,
            ));
        }
    }

    let mut items = Vec::with_capacity(model.features.len());
    for (index, grouped) in model.features.iter().enumerate() {
        let id = grouped.group;
        let tolerance = det::exp(log_tolerances[index]);
        validate_positive_derived(
            GroupedDerivedQuantity::Tolerance,
            Some(index),
            Some(id),
            tolerance,
        )?;
        let published_log_tolerance = det::ln(tolerance);
        let cost_contribution = det::exp(feature_log_cost[index] - published_log_tolerance);
        validate_positive_derived(
            GroupedDerivedQuantity::CostContribution,
            Some(index),
            Some(id),
            cost_contribution,
        )?;
        let log_loading = feature_log_sensitivity[index] + published_log_tolerance - log_k;
        let standard_deviation_loading = det::exp(log_loading);
        validate_positive_derived(
            GroupedDerivedQuantity::StandardDeviationLoading,
            Some(index),
            Some(id),
            standard_deviation_loading,
        )?;
        let independent_variance_contribution =
            standard_deviation_loading * standard_deviation_loading;
        validate_positive_derived(
            GroupedDerivedQuantity::IndependentVarianceContribution,
            Some(index),
            Some(id),
            independent_variance_contribution,
        )?;
        items.push(GroupedAllocationItem {
            feature_index: index,
            group: id,
            tolerance,
            action: action_for(tolerance, grouped.feature.baseline_tolerance),
            cost_contribution,
            standard_deviation_loading,
            independent_variance_contribution,
        });
    }

    let mut group_receipts = Vec::with_capacity(model.groups.len());
    for (group_index, members) in group_features.iter().enumerate() {
        let id = DependencyGroupId(group_index as u16);
        let standard_deviation = compensated_sum(
            members
                .iter()
                .map(|index| items[*index].standard_deviation_loading),
        );
        validate_positive_derived(
            GroupedDerivedQuantity::GroupStandardDeviation,
            None,
            Some(id),
            standard_deviation,
        )?;
        let variance = standard_deviation * standard_deviation;
        validate_positive_derived(
            GroupedDerivedQuantity::GroupVariance,
            None,
            Some(id),
            variance,
        )?;
        let independent_variance = compensated_sum(
            members
                .iter()
                .map(|index| items[*index].independent_variance_contribution),
        );
        validate_positive_derived(
            GroupedDerivedQuantity::GroupIndependentVariance,
            None,
            Some(id),
            independent_variance,
        )?;
        let dependency_variance_delta = variance - independent_variance;
        validate_dependency_variance_delta(Some(id), members.len() > 1, dependency_variance_delta)?;
        let total_cost =
            compensated_sum(members.iter().map(|index| items[*index].cost_contribution));
        validate_positive_derived(
            GroupedDerivedQuantity::GroupCost,
            None,
            Some(id),
            total_cost,
        )?;
        group_receipts.push(DependencyGroupReceipt {
            group: id,
            feature_count: members.len(),
            log_shape: group_log_shapes[group_index],
            standard_deviation,
            variance,
            independent_variance,
            dependency_variance_delta,
            total_cost,
        });
    }

    let total_cost = compensated_sum(items.iter().map(GroupedAllocationItem::cost_contribution));
    validate_positive_derived(GroupedDerivedQuantity::TotalCost, None, None, total_cost)?;
    let closed_form_cost = det::exp(1.5 * global_log_shape - log_k - 0.5 * log_budget);
    validate_positive_derived(
        GroupedDerivedQuantity::ClosedFormCost,
        None,
        None,
        closed_form_cost,
    )?;
    let cost_residual = total_cost - closed_form_cost;
    if !cost_residual.is_finite() {
        return Err(invalid_derived(
            GroupedDerivedQuantity::TotalCost,
            None,
            None,
            ScalarIssue::NonFinite,
        ));
    }

    let achieved_variance =
        compensated_sum(group_receipts.iter().map(DependencyGroupReceipt::variance));
    validate_positive_derived(
        GroupedDerivedQuantity::AchievedVariance,
        None,
        None,
        achieved_variance,
    )?;
    let independent_variance = compensated_sum(
        group_receipts
            .iter()
            .map(DependencyGroupReceipt::independent_variance),
    );
    validate_positive_derived(
        GroupedDerivedQuantity::IndependentVariance,
        None,
        None,
        independent_variance,
    )?;
    let dependency_variance_delta = achieved_variance - independent_variance;
    validate_dependency_variance_delta(
        None,
        group_features.iter().any(|members| members.len() > 1),
        dependency_variance_delta,
    )?;

    let mut stationarity_terms = Vec::with_capacity(model.features.len());
    for (index, grouped) in model.features.iter().enumerate() {
        let group_log_standard_deviation =
            det::ln(group_receipts[grouped.group.index()].standard_deviation);
        let published_log_tolerance = det::ln(items[index].tolerance);
        let term = feature_log_cost[index]
            - 2.0 * published_log_tolerance
            - feature_log_sensitivity[index]
            - group_log_standard_deviation;
        if !term.is_finite() {
            return Err(invalid_derived(
                GroupedDerivedQuantity::StationarityResidual,
                Some(index),
                Some(grouped.group),
                ScalarIssue::NonFinite,
            ));
        }
        stationarity_terms.push(term);
    }
    let minimum_stationarity = stationarity_terms
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let maximum_stationarity = stationarity_terms
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let max_stationarity_log_residual = maximum_stationarity - minimum_stationarity;
    validate_nonnegative_derived(
        GroupedDerivedQuantity::StationarityResidual,
        None,
        None,
        max_stationarity_log_residual,
    )?;

    Ok(GroupedAllocationReceipt {
        schema_version: GROUPED_ALLOCATION_SCHEMA_V1,
        model: model.clone(),
        variance_budget,
        tolerance_to_sigma: k,
        items: items.into_boxed_slice(),
        groups: group_receipts.into_boxed_slice(),
        total_cost,
        closed_form_cost,
        cost_residual,
        log_scale_correction: correction,
        achieved_variance,
        budget_residual: achieved_variance - variance_budget,
        independent_variance,
        dependency_variance_delta,
        max_stationarity_log_residual,
    })
}

#[cfg(test)]
mod tests {
    use super::{
        DependencyGroupId, GroupedAllocationError, GroupedDerivedQuantity, ScalarIssue,
        validate_dependency_variance_delta,
    };

    #[test]
    fn coherent_dependency_delta_guard_refuses_zero_at_group_and_model_scope() {
        for group in [Some(DependencyGroupId(3)), None] {
            assert_eq!(
                validate_dependency_variance_delta(group, true, 0.0),
                Err(GroupedAllocationError::InvalidDerived {
                    quantity: GroupedDerivedQuantity::DependencyVarianceDelta,
                    feature_index: None,
                    group,
                    issue: ScalarIssue::Underflow,
                })
            );
        }
        assert_eq!(validate_dependency_variance_delta(None, false, 0.0), Ok(()));
    }
}
