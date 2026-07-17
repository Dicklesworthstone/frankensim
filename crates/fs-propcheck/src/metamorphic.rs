//! Gauntlet G3 metamorphic relations over the G0 generator, replay, and
//! shrink engine.
//!
//! This module is dependency-free on purpose. Consumer crates provide their
//! domain transforms (motors, unit systems, refinement requests, derivative
//! directions, or conversion routes) while this layer owns deterministic
//! relation identity, joint input/transform shrinking, tolerance metadata,
//! and structured failure receipts.

use crate::{Shrink, Stream, StructuredFailure, StructuredVerdict, check_structured};

/// The canonical G3 relation families from the FrankenSim Gauntlet.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CanonicalRelation {
    /// Apply one rigid motion to all geometric inputs.
    RigidMotion,
    /// Express one physical problem in a different coherent unit system.
    UnitRescaling,
    /// Refine an estimator or discretization without worsening its claim.
    RefinementMonotonicity,
    /// Compare an adjoint directional derivative with finite differences.
    AdjointFiniteDifference,
    /// Compare two admissible representation-conversion paths.
    ConversionPathIndependence,
    /// Map dimensional and nondimensionalized formulations of one problem
    /// onto each other through the declared scaling (6nb.4 patch Rev A).
    RegimeScalingCoherence,
}

impl CanonicalRelation {
    /// Stable receipt label.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::RigidMotion => "rigid-motion",
            Self::UnitRescaling => "unit-rescaling",
            Self::RefinementMonotonicity => "refinement-monotonicity",
            Self::AdjointFiniteDifference => "adjoint-finite-difference",
            Self::ConversionPathIndependence => "conversion-path-independence",
            Self::RegimeScalingCoherence => "regime-scaling-coherence",
        }
    }
}

/// Declared scalar tolerance semantics for a relation.
///
/// Custom output comparators receive this value explicitly. The built-in
/// [`Tolerance::evaluate_scalar`] helper applies it fail-closed to finite
/// `f64` observations and returns a signed admission margin.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Tolerance {
    /// Finite IEEE-754 values must have identical bits.
    Exact,
    /// `|candidate - reference| <= max_abs`.
    Absolute {
        /// Maximum admitted absolute error.
        max_abs: f64,
    },
    /// Absolute-or-relative tolerance using the larger finite allowance.
    AbsoluteRelative {
        /// Absolute error floor.
        max_abs: f64,
        /// Relative error multiplier applied to the larger magnitude.
        max_relative: f64,
    },
    /// The candidate may increase by at most the declared slack.
    NonIncreasing {
        /// Maximum admitted increase.
        max_increase: f64,
    },
}

/// Invalid tolerance metadata.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToleranceError {
    /// A tolerance limit is NaN or infinite.
    NonFinite(&'static str),
    /// A tolerance limit is negative.
    Negative(&'static str),
}

impl core::fmt::Display for ToleranceError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonFinite(name) => write!(formatter, "relation tolerance `{name}` is not finite"),
            Self::Negative(name) => write!(formatter, "relation tolerance `{name}` is negative"),
        }
    }
}

impl std::error::Error for ToleranceError {}

impl Tolerance {
    /// Validate all numeric limits before case generation begins.
    ///
    /// # Errors
    /// [`ToleranceError`] for a non-finite or negative limit.
    pub fn validate(self) -> Result<(), ToleranceError> {
        fn limit(name: &'static str, value: f64) -> Result<(), ToleranceError> {
            if !value.is_finite() {
                Err(ToleranceError::NonFinite(name))
            } else if value < 0.0 {
                Err(ToleranceError::Negative(name))
            } else {
                Ok(())
            }
        }

        match self {
            Self::Exact => Ok(()),
            Self::Absolute { max_abs } => limit("max_abs", max_abs),
            Self::AbsoluteRelative {
                max_abs,
                max_relative,
            } => {
                limit("max_abs", max_abs)?;
                limit("max_relative", max_relative)
            }
            Self::NonIncreasing { max_increase } => limit("max_increase", max_increase),
        }
    }

    /// Compare two finite scalars and return a signed admission margin.
    /// Non-finite observations and overflowed error calculations fail closed.
    #[must_use]
    pub fn evaluate_scalar(self, reference: f64, candidate: f64) -> RelationObservation {
        if self.validate().is_err() {
            return RelationObservation::new(
                -f64::MAX,
                "invalid tolerance metadata reached scalar comparison",
            );
        }
        if !reference.is_finite() || !candidate.is_finite() {
            return RelationObservation::new(
                -f64::MAX,
                "relation scalar observations must be finite",
            );
        }

        match self {
            Self::Exact => {
                if reference.to_bits() == candidate.to_bits() {
                    RelationObservation::new(0.0, "finite scalar bits match exactly")
                } else {
                    RelationObservation::new(-1.0, "finite scalar bits differ")
                }
            }
            Self::Absolute { max_abs } => {
                let error = (candidate - reference).abs();
                finite_margin(max_abs, error, "absolute scalar tolerance")
            }
            Self::AbsoluteRelative {
                max_abs,
                max_relative,
            } => {
                let relative = max_relative * reference.abs().max(candidate.abs());
                if !relative.is_finite() {
                    return RelationObservation::new(
                        -f64::MAX,
                        "relative scalar allowance overflowed",
                    );
                }
                let error = (candidate - reference).abs();
                finite_margin(
                    max_abs.max(relative),
                    error,
                    "absolute-relative scalar tolerance",
                )
            }
            Self::NonIncreasing { max_increase } => {
                let increase = if candidate <= reference {
                    0.0
                } else {
                    candidate - reference
                };
                finite_margin(max_increase, increase, "non-increasing scalar tolerance")
            }
        }
    }
}

fn finite_margin(allowance: f64, error: f64, detail: &'static str) -> RelationObservation {
    if !allowance.is_finite() || !error.is_finite() {
        RelationObservation::new(-f64::MAX, format!("{detail} arithmetic overflowed"))
    } else {
        RelationObservation::new(allowance - error, detail)
    }
}

/// One relation evaluation. A nonnegative finite margin is admitted.
#[derive(Debug, Clone, PartialEq)]
pub struct RelationObservation {
    margin: f64,
    detail: String,
}

impl RelationObservation {
    /// Construct an observation from a finite signed margin.
    ///
    /// # Panics
    /// Panics when `margin` is NaN or infinite; relation evidence must remain
    /// numerically decidable.
    #[must_use]
    pub fn new(margin: f64, detail: impl Into<String>) -> Self {
        assert!(
            margin.is_finite(),
            "relation admission margin must be finite"
        );
        Self {
            margin,
            detail: detail.into(),
        }
    }

    /// Whether the relation admitted this observation.
    #[must_use]
    pub fn admitted(&self) -> bool {
        self.margin >= 0.0
    }

    /// Signed admission margin; nonnegative values pass.
    #[must_use]
    pub fn margin(&self) -> f64 {
        self.margin
    }

    /// Relation-specific diagnostic.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

/// A generated base input paired with the transform to apply.
#[derive(Debug, Clone, PartialEq)]
pub struct RelationCase<Input, Transform> {
    /// Base operator input.
    pub input: Input,
    /// Relation transform parameters.
    pub transform: Transform,
}

impl<Input, Transform> RelationCase<Input, Transform> {
    /// Pair an input with its relation transform.
    #[must_use]
    pub const fn new(input: Input, transform: Transform) -> Self {
        Self { input, transform }
    }
}

impl<Input: Shrink, Transform: Shrink> Shrink for RelationCase<Input, Transform> {
    fn shrink_candidates(&self) -> Vec<Self> {
        (self.input.clone(), self.transform.clone())
            .shrink_candidates()
            .into_iter()
            .map(|(input, transform)| Self { input, transform })
            .collect()
    }
}

/// An operator that can be exercised by the metamorphic harness.
pub trait MetamorphicOperator<Input, Output> {
    /// Evaluate one admitted input.
    fn evaluate(&self, input: &Input) -> Output;
}

impl<Input, Output, Function> MetamorphicOperator<Input, Output> for Function
where
    Function: Fn(&Input) -> Output,
{
    fn evaluate(&self, input: &Input) -> Output {
        self(input)
    }
}

/// A declared transform and output expectation for one canonical G3 family.
pub trait MetamorphicRelation<Input, Output> {
    /// Shrinkable parameters selecting the input transformation.
    type Transform: Shrink + core::fmt::Debug;

    /// Stable relation identity within the operator battery.
    fn id(&self) -> &'static str;

    /// Canonical G3 family.
    fn kind(&self) -> CanonicalRelation;

    /// Declared tolerance semantics.
    fn tolerance(&self) -> Tolerance;

    /// Produce the transformed input.
    fn transform(&self, input: &Input, transform: &Self::Transform) -> Input;

    /// Compare base and transformed outputs.
    fn compare(
        &self,
        base: &Output,
        transformed: &Output,
        transform: &Self::Transform,
    ) -> RelationObservation;
}

/// Closure-backed declaration used by the five canonical constructors.
pub struct DeclaredRelation<Transform, TransformInput, CompareOutput> {
    id: &'static str,
    kind: CanonicalRelation,
    tolerance: Tolerance,
    transform_input: TransformInput,
    compare_output: CompareOutput,
    transform_marker: core::marker::PhantomData<fn() -> Transform>,
}

impl<Transform, TransformInput, CompareOutput>
    DeclaredRelation<Transform, TransformInput, CompareOutput>
{
    fn new(
        id: &'static str,
        kind: CanonicalRelation,
        tolerance: Tolerance,
        transform_input: TransformInput,
        compare_output: CompareOutput,
    ) -> Self {
        assert!(
            !id.trim().is_empty(),
            "metamorphic relation id must not be empty"
        );
        tolerance
            .validate()
            .unwrap_or_else(|error| panic!("invalid metamorphic relation `{id}`: {error}"));
        Self {
            id,
            kind,
            tolerance,
            transform_input,
            compare_output,
            transform_marker: core::marker::PhantomData,
        }
    }

    /// Stable relation identity.
    #[must_use]
    pub const fn id(&self) -> &'static str {
        self.id
    }

    /// Canonical relation family.
    #[must_use]
    pub const fn kind(&self) -> CanonicalRelation {
        self.kind
    }

    /// Declared tolerance semantics.
    #[must_use]
    pub const fn tolerance(&self) -> Tolerance {
        self.tolerance
    }
}

impl<Input, Output, Transform, TransformInput, CompareOutput> MetamorphicRelation<Input, Output>
    for DeclaredRelation<Transform, TransformInput, CompareOutput>
where
    Transform: Shrink + core::fmt::Debug,
    TransformInput: Fn(&Input, &Transform) -> Input,
    CompareOutput: Fn(&Output, &Output, &Transform, Tolerance) -> RelationObservation,
{
    type Transform = Transform;

    fn id(&self) -> &'static str {
        self.id
    }

    fn kind(&self) -> CanonicalRelation {
        self.kind
    }

    fn tolerance(&self) -> Tolerance {
        self.tolerance
    }

    fn transform(&self, input: &Input, transform: &Transform) -> Input {
        (self.transform_input)(input, transform)
    }

    fn compare(
        &self,
        base: &Output,
        transformed: &Output,
        transform: &Transform,
    ) -> RelationObservation {
        (self.compare_output)(base, transformed, transform, self.tolerance)
    }
}

macro_rules! canonical_constructor {
    ($function:ident, $kind:ident, $docs:literal) => {
        #[doc = $docs]
        #[must_use]
        pub fn $function<Input, Output, Transform, TransformInput, CompareOutput>(
            id: &'static str,
            tolerance: Tolerance,
            transform_input: TransformInput,
            compare_output: CompareOutput,
        ) -> DeclaredRelation<Transform, TransformInput, CompareOutput>
        where
            Transform: Shrink + core::fmt::Debug,
            TransformInput: Fn(&Input, &Transform) -> Input,
            CompareOutput: Fn(&Output, &Output, &Transform, Tolerance) -> RelationObservation,
        {
            DeclaredRelation::new(
                id,
                CanonicalRelation::$kind,
                tolerance,
                transform_input,
                compare_output,
            )
        }
    };
}

canonical_constructor!(
    rigid_motion,
    RigidMotion,
    "Declare a rigid-motion relation. Domain motor semantics stay in the consumer."
);
canonical_constructor!(
    unit_rescaling,
    UnitRescaling,
    "Declare a coherent unit-rescaling relation. Unit transforms stay in the consumer."
);
canonical_constructor!(
    refinement_monotonicity,
    RefinementMonotonicity,
    "Declare a refinement-monotonicity relation. The output comparator owns estimator meaning."
);
canonical_constructor!(
    adjoint_finite_difference,
    AdjointFiniteDifference,
    "Declare an adjoint-versus-finite-difference relation. Direction and step are transform data."
);
canonical_constructor!(
    conversion_path_independence,
    ConversionPathIndependence,
    "Declare a representation path-independence relation. Route certificates stay in the consumer."
);
canonical_constructor!(
    regime_scaling_coherence,
    RegimeScalingCoherence,
    "Declare a regime-scaling coherence relation. The declared scaling map stays in the consumer."
);

/// Run one declared relation through the standard deterministic case stream,
/// replay selector, joint shrink descent, and JSONL failure artifact.
///
/// A returned violation's failure context records the final shrunk relation
/// id/family, transform, both outputs, tolerance, and signed margin. Transform
/// or operator panics are caught by the shared runner and remain generic panic
/// failures rather than skipped cases; outputs that were never produced cannot
/// appear in that panic receipt.
///
/// # Panics
/// On invalid harness metadata, malformed replay configuration, or the first
/// shrunk relation violation or transform/operator/comparator panic.
pub fn check_relation<Input, Output, Operator, Relation>(
    operator_name: &str,
    suite_seed: u64,
    cases: u64,
    generate: impl Fn(
        &mut Stream,
    ) -> RelationCase<
        Input,
        <Relation as MetamorphicRelation<Input, Output>>::Transform,
    >,
    operator: &Operator,
    relation: &Relation,
) where
    Input: Shrink + core::fmt::Debug,
    Output: core::fmt::Debug,
    Operator: MetamorphicOperator<Input, Output>,
    Relation: MetamorphicRelation<Input, Output>,
{
    assert!(
        !operator_name.trim().is_empty(),
        "metamorphic operator name must not be empty"
    );
    assert!(
        !relation.id().trim().is_empty(),
        "metamorphic relation id must not be empty"
    );
    assert!(cases > 0, "metamorphic relation requires at least one case");
    relation.tolerance().validate().unwrap_or_else(|error| {
        panic!("invalid metamorphic relation `{}`: {error}", relation.id())
    });

    let property_name = format!("{operator_name}::{}", relation.id());
    check_structured(&property_name, suite_seed, cases, generate, |case| {
        let base_output = operator.evaluate(&case.input);
        let transformed_input = relation.transform(&case.input, &case.transform);
        let transformed_output = operator.evaluate(&transformed_input);
        let observation = relation.compare(&base_output, &transformed_output, &case.transform);
        if observation.admitted() {
            StructuredVerdict::Pass
        } else {
            StructuredVerdict::Fail(StructuredFailure::new(
                "relation-violation",
                observation.detail(),
                vec![
                    ("relation_id", relation.id().to_string()),
                    ("relation_kind", relation.kind().slug().to_string()),
                    ("transform", format!("{:?}", case.transform)),
                    ("base_output", format!("{base_output:?}")),
                    ("transformed_output", format!("{transformed_output:?}")),
                    ("tolerance", format!("{:?}", relation.tolerance())),
                    ("margin", format!("{:.17e}", observation.margin())),
                ],
            ))
        }
    });
}
