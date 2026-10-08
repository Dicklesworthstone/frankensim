//! Typed Maslov--Krein--Evans theorem statements and hypothesis lattices.
//!
//! This module is deliberately a statement boundary, not a proof engine. It
//! keeps four signed-count objects nonfungible, records every convention and
//! theorem hypothesis, preregisters falsifiers, and gives the resulting
//! lattice a deterministic identity. Opaque witness and proof-artifact IDs are
//! references only: successful validation never upgrades them into scientific
//! authority.

use core::fmt;
use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, NeverCancel, ProblemSemanticId, WireType,
};
use fs_exec::Cx;

/// Current semantic and wire version of the bridge-lattice schema.
pub const BRIDGE_LATTICE_SCHEMA_VERSION_V1: u32 = 1;
/// Hard cap on theorem nodes, checked before sorting or graph work.
pub const MAX_BRIDGE_NODES_V1: usize = 128;
/// Hard cap on implications, checked before graph construction.
pub const MAX_BRIDGE_IMPLICATIONS_V1: usize = 512;
/// Hard cap on hypotheses retained by one theorem node.
pub const MAX_BRIDGE_HYPOTHESES_PER_NODE_V1: usize = 64;
/// Hard cap on falsifiers retained by one theorem node.
pub const MAX_BRIDGE_FALSIFIERS_PER_NODE_V1: usize = 32;

const BRIDGE_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(1 << 20, 1 << 19, 16, 8192, 8192);

trait DigestBytes {
    fn digest_bytes(&self) -> &[u8; 32];
}

macro_rules! opaque_bridge_id {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name([u8; 32]);

        impl $name {
            /// Construct a typed content reference from exact digest bytes.
            /// The bytes alone confer no theorem authority.
            #[must_use]
            pub const fn from_bytes(bytes: [u8; 32]) -> Self {
                Self(bytes)
            }

            /// Exact digest bytes.
            #[must_use]
            pub const fn as_bytes(&self) -> &[u8; 32] {
                &self.0
            }
        }

        impl DigestBytes for $name {
            fn digest_bytes(&self) -> &[u8; 32] {
                self.as_bytes()
            }
        }
    };
}

opaque_bridge_id!(
    /// Exact theorem-domain identity.
    BridgeDomainIdV1
);
opaque_bridge_id!(
    /// Versioned Hamiltonian, monodromy, or spatial-dynamics operator family.
    BridgeOperatorFamilyIdV1
);
opaque_bridge_id!(
    /// Versioned parameter path and direction convention.
    BridgeParameterizationIdV1
);
opaque_bridge_id!(
    /// Exact Lagrangian path used by the Maslov index.
    LagrangianPathIdV1
);
opaque_bridge_id!(
    /// Reference Lagrangian subspace or boundary plane.
    LagrangianReferenceIdV1
);
opaque_bridge_id!(
    /// Symplectic-form artifact.
    BridgeSymplecticFormIdV1
);
opaque_bridge_id!(
    /// Crossing-form family and its sign convention.
    CrossingFormFamilyIdV1
);
opaque_bridge_id!(
    /// Pontryagin-space or indefinite-metric path.
    PontryaginPathIdV1
);
opaque_bridge_id!(
    /// Nondegenerate Krein-form artifact.
    BridgeKreinFormIdV1
);
opaque_bridge_id!(
    /// Analytic Fredholm or Evans operator family.
    AnalyticFredholmFamilyIdV1
);
opaque_bridge_id!(
    /// Analytic domain on which an Evans family is stated.
    EvansAnalyticDomainIdV1
);
opaque_bridge_id!(
    /// Oriented Evans/argument-principle contour.
    EvansContourIdV1
);
opaque_bridge_id!(
    /// Evans-function normalization and branch convention.
    EvansNormalizationIdV1
);
opaque_bridge_id!(
    /// Essential-spectrum exclusion or exponential-dichotomy artifact.
    EssentialSpectrumExclusionIdV1
);
opaque_bridge_id!(
    /// Multiplicity convention shared by a claimed count relation.
    BridgeMultiplicityRuleIdV1
);
opaque_bridge_id!(
    /// Physical machine-instability definition/model artifact.
    MachineInstabilityModelIdV1
);
opaque_bridge_id!(
    /// Exact map between two stated count conventions or domains.
    BridgeCorrespondenceMapIdV1
);
opaque_bridge_id!(
    /// Stable identity chosen for one theorem node.
    BridgeTheoremNodeIdV1
);
opaque_bridge_id!(
    /// Retained hypothesis witness reference.
    BridgeHypothesisWitnessIdV1
);
opaque_bridge_id!(
    /// Retained proof artifact reference.
    BridgeProofArtifactIdV1
);
opaque_bridge_id!(
    /// Proof checker/verifier implementation identity.
    BridgeVerifierIdV1
);
opaque_bridge_id!(
    /// Policy under which an external proof artifact was reviewed.
    BridgeProofPolicyIdV1
);
opaque_bridge_id!(
    /// Formal-system or native checker version.
    BridgeFormalSystemVersionIdV1
);
opaque_bridge_id!(
    /// Version identity of the theorem statement vocabulary.
    BridgeStatementVersionIdV1
);
opaque_bridge_id!(
    /// Trusted-computing-base component-set identity.
    BridgeTcbIdV1
);
opaque_bridge_id!(
    /// Explicit no-claim statement.
    BridgeNoClaimIdV1
);
opaque_bridge_id!(
    /// Retained falsifier fixture or oracle identity.
    BridgeFalsifierArtifactIdV1
);
opaque_bridge_id!(
    /// Retained mathematical or numerical counterexample identity.
    BridgeCounterexampleIdV1
);

/// Domain-separated identity schema for one canonical bridge lattice.
pub enum BridgeLatticeIdentitySchemaV1 {}

impl CanonicalSchema for BridgeLatticeIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-spectral.maslov-krein-evans-bridge.v1";
    const NAME: &'static str = "maslov-krein-evans-hypothesis-lattice";
    const VERSION: u32 = BRIDGE_LATTICE_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "typed count objects, exact theorem domains, hypotheses, implications, falsifiers, proof state, TCB, and validation budget";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("statement-version", WireType::Bytes),
        FieldSpec::required("count-objects", WireType::Bytes),
        FieldSpec::required("theorem-nodes", WireType::CanonicalSet),
        FieldSpec::required("implications", WireType::CanonicalSet),
        FieldSpec::required("tcb", WireType::Bytes),
        FieldSpec::required("budget", WireType::Bytes),
    ];
}

/// Typed deterministic identity of one validated statement lattice.
pub type BridgeLatticeIdV1 = ProblemSemanticId<BridgeLatticeIdentitySchemaV1>;

/// Positive versus reversed orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CountOrientationV1 {
    /// Canonical positive direction.
    Positive,
    /// Reversed direction.
    Negative,
}

impl CountOrientationV1 {
    const fn sign(self) -> i64 {
        match self {
            Self::Positive => 1,
            Self::Negative => -1,
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Positive => 0,
            Self::Negative => 1,
        }
    }
}

/// Endpoint weight expressed in doubled-index units.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EndpointWeightV1 {
    /// Do not count the endpoint crossing.
    Excluded,
    /// Count half the crossing signature.
    Half,
    /// Count the full crossing signature.
    Full,
}

impl EndpointWeightV1 {
    const fn doubled_weight(self) -> i64 {
        match self {
            Self::Excluded => 0,
            Self::Half => 1,
            Self::Full => 2,
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Excluded => 0,
            Self::Half => 1,
            Self::Full => 2,
        }
    }
}

/// Explicit left/right endpoint convention.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EndpointConventionV1 {
    /// Weight assigned to the initial endpoint in the chosen orientation.
    pub left: EndpointWeightV1,
    /// Weight assigned to the final endpoint in the chosen orientation.
    pub right: EndpointWeightV1,
}

/// Complete signed-count convention relevant to executable transformations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SignedCountConventionV1 {
    /// Parameter, contour, or geometric orientation.
    pub orientation: CountOrientationV1,
    /// Endpoint convention.
    pub endpoints: EndpointConventionV1,
}

/// Endpoint signatures in canonical positive orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EndpointSignatureTraceV1 {
    /// Left signature, or `None` when it has not been resolved.
    pub left: Option<i64>,
    /// Right signature, or `None` when it has not been resolved.
    pub right: Option<i64>,
}

/// A signed count multiplied by two so half-signature endpoint conventions are
/// represented exactly without floating point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DoubledSignedCountV1(
    /// Exact doubled-index value.
    pub i64,
);

/// Executable affine transformation between two signed-count conventions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SignedCountTransformV1 {
    /// Exact theorem/correspondence map that justifies the transformation.
    pub map: BridgeCorrespondenceMapIdV1,
    /// Multiply the doubled source count by this sign (`1` or `-1`).
    sign: i8,
    /// Add this exact doubled-count endpoint correction.
    doubled_offset: i64,
}

impl SignedCountTransformV1 {
    /// Construct an explicit transform. Only signs `1` and `-1` are admitted.
    #[must_use]
    pub const fn new(
        map: BridgeCorrespondenceMapIdV1,
        sign: i8,
        doubled_offset: i64,
    ) -> Option<Self> {
        if sign == 1 || sign == -1 {
            Some(Self {
                map,
                sign,
                doubled_offset,
            })
        } else {
            None
        }
    }

    /// Multiplicative sign.
    #[must_use]
    pub const fn sign(self) -> i8 {
        self.sign
    }

    /// Exact doubled-count offset.
    #[must_use]
    pub const fn doubled_offset(self) -> i64 {
        self.doubled_offset
    }

    /// Apply the transformation, refusing integer overflow.
    #[must_use]
    pub fn apply(self, count: DoubledSignedCountV1) -> Option<DoubledSignedCountV1> {
        count
            .0
            .checked_mul(i64::from(self.sign))
            .and_then(|value| value.checked_add(self.doubled_offset))
            .map(DoubledSignedCountV1)
    }
}

/// Failure to derive a convention change from explicit endpoint data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConventionTransformErrorV1 {
    /// A changed left endpoint weight requires a missing signature.
    LeftEndpointUnresolved,
    /// A changed right endpoint weight requires a missing signature.
    RightEndpointUnresolved,
    /// Exact integer arithmetic overflowed.
    Overflow,
}

/// Derive the executable orientation/endpoint transform between conventions.
///
/// Endpoint signatures are expressed in canonical positive orientation. The
/// result therefore handles orientation reversal and endpoint-rule mutation in
/// one exact doubled-integer calculation.
///
/// # Errors
/// Returns [`ConventionTransformErrorV1`] when a changed endpoint convention
/// lacks its signature or exact arithmetic overflows.
pub fn derive_convention_transform_v1(
    map: BridgeCorrespondenceMapIdV1,
    source: SignedCountConventionV1,
    target: SignedCountConventionV1,
    endpoints: EndpointSignatureTraceV1,
) -> Result<SignedCountTransformV1, ConventionTransformErrorV1> {
    let sign_i64 = target.orientation.sign() * source.orientation.sign();
    let (source_left, source_right) = canonical_endpoint_weights(source);
    let (target_left, target_right) = canonical_endpoint_weights(target);
    let mut doubled_offset = 0_i64;
    for (source_weight, target_weight, signature, unresolved) in [
        (
            source_left,
            target_left,
            endpoints.left,
            ConventionTransformErrorV1::LeftEndpointUnresolved,
        ),
        (
            source_right,
            target_right,
            endpoints.right,
            ConventionTransformErrorV1::RightEndpointUnresolved,
        ),
    ] {
        let weight_delta = target_weight.doubled_weight() - source_weight.doubled_weight();
        if weight_delta == 0 {
            continue;
        }
        let signature = signature.ok_or(unresolved)?;
        let contribution = weight_delta
            .checked_mul(signature)
            .and_then(|value| value.checked_mul(target.orientation.sign()))
            .ok_or(ConventionTransformErrorV1::Overflow)?;
        doubled_offset = doubled_offset
            .checked_add(contribution)
            .ok_or(ConventionTransformErrorV1::Overflow)?;
    }
    SignedCountTransformV1::new(map, sign_i64 as i8, doubled_offset)
        .ok_or(ConventionTransformErrorV1::Overflow)
}

fn canonical_endpoint_weights(
    convention: SignedCountConventionV1,
) -> (EndpointWeightV1, EndpointWeightV1) {
    match convention.orientation {
        CountOrientationV1::Positive => (convention.endpoints.left, convention.endpoints.right),
        CountOrientationV1::Negative => (convention.endpoints.right, convention.endpoints.left),
    }
}

/// Neutral-signature handling for the Krein count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NeutralKreinPolicyV1 {
    /// Neutral directions are excluded by a retained hypothesis witness.
    Excluded {
        /// Exact exclusion witness.
        witness: BridgeHypothesisWitnessIdV1,
    },
    /// A named perturbation/splitting map resolves neutral directions.
    Resolved {
        /// Resolution map.
        map: BridgeCorrespondenceMapIdV1,
        /// Boundary on what the resolution does not establish.
        no_claim: BridgeNoClaimIdV1,
    },
    /// Neutral signature remains unresolved; equality nodes must stay
    /// conjectural or refuted.
    Unresolved {
        /// Explicit no-claim artifact.
        no_claim: BridgeNoClaimIdV1,
    },
}

/// Typed Maslov-index object. It is not interchangeable with spectral flow.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MaslovIndexObjectV1 {
    /// Lagrangian path.
    pub path: LagrangianPathIdV1,
    /// Reference Lagrangian.
    pub reference: LagrangianReferenceIdV1,
    /// Supporting symplectic form.
    pub symplectic_form: BridgeSymplecticFormIdV1,
    /// Crossing-form family and sign convention.
    pub crossing_forms: CrossingFormFamilyIdV1,
    /// Signed-count convention.
    pub convention: SignedCountConventionV1,
    /// Multiplicity convention.
    pub multiplicity: BridgeMultiplicityRuleIdV1,
}

/// Typed Krein spectral-flow object. It is not a Maslov index by label.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KreinSpectralFlowObjectV1 {
    /// Pontryagin/Krein operator path.
    pub path: PontryaginPathIdV1,
    /// Supporting nondegenerate indefinite form.
    pub krein_form: BridgeKreinFormIdV1,
    /// Neutral-signature policy.
    pub neutral_policy: NeutralKreinPolicyV1,
    /// Signed-count convention.
    pub convention: SignedCountConventionV1,
    /// Multiplicity convention.
    pub multiplicity: BridgeMultiplicityRuleIdV1,
}

/// Typed Evans/argument-principle winding object.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EvansWindingObjectV1 {
    /// Analytic Fredholm/Evans operator family.
    pub family: AnalyticFredholmFamilyIdV1,
    /// Analytic domain.
    pub analytic_domain: EvansAnalyticDomainIdV1,
    /// Oriented contour.
    pub contour: EvansContourIdV1,
    /// Essential-spectrum exclusion artifact.
    pub essential_spectrum_exclusion: EssentialSpectrumExclusionIdV1,
    /// Evans normalization/branch convention.
    pub normalization: EvansNormalizationIdV1,
    /// Signed-count convention.
    pub convention: SignedCountConventionV1,
    /// Multiplicity convention.
    pub multiplicity: BridgeMultiplicityRuleIdV1,
}

/// Typed machine-instability count. This remains distinct from all spectral
/// counts until a theorem node supplies an explicit interpretation map.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MachineInstabilityCountObjectV1 {
    /// Physical instability definition and admitted model family.
    pub model: MachineInstabilityModelIdV1,
    /// Physical parameter path.
    pub parameterization: BridgeParameterizationIdV1,
    /// Signed-count convention.
    pub convention: SignedCountConventionV1,
    /// Multiplicity convention.
    pub multiplicity: BridgeMultiplicityRuleIdV1,
}

/// Four nonfungible signed-count kinds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BridgeCountKindV1 {
    /// Maslov index.
    Maslov,
    /// Krein spectral flow.
    Krein,
    /// Evans winding/argument-principle count.
    Evans,
    /// Physical machine-instability count.
    MachineInstability,
}

impl BridgeCountKindV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::Maslov => 0,
            Self::Krein => 1,
            Self::Evans => 2,
            Self::MachineInstability => 3,
        }
    }
}

/// Exact theorem-domain class. Dimensions and extension maps are statement
/// data, not inferred from a label.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BridgeTheoremScopeV1 {
    /// Classical finite-dimensional Hamiltonian path.
    ClassicalFiniteHamiltonian {
        /// Even phase-space dimension.
        phase_dimension: u32,
    },
    /// Periodic monodromy/Floquet extension.
    PeriodicMonodromy {
        /// State dimension of the return map.
        state_dimension: u32,
        /// Exact flow-to-monodromy correspondence.
        monodromy_map: BridgeCorrespondenceMapIdV1,
    },
    /// Spatial-dynamics analytic Fredholm/Evans extension.
    SpatialDynamicsEvans {
        /// Stable-bundle dimension on the admitted contour/domain.
        stable_dimension: u32,
        /// Unstable-bundle dimension on the admitted contour/domain.
        unstable_dimension: u32,
        /// Exact spatial-dynamics correspondence map.
        spatial_map: BridgeCorrespondenceMapIdV1,
    },
    /// Bold maximal family connecting all three count objects.
    MaximalMaslovKreinEvans {
        /// Finite-to-periodic extension map.
        finite_to_periodic: BridgeCorrespondenceMapIdV1,
        /// Periodic-to-spatial/Evans extension map.
        periodic_to_evans: BridgeCorrespondenceMapIdV1,
    },
    /// Explicit machine-instability corollary.
    MachineInstabilityCorollary {
        /// Spectral-to-physical interpretation map.
        interpretation: BridgeCorrespondenceMapIdV1,
    },
}

impl BridgeTheoremScopeV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::ClassicalFiniteHamiltonian { .. } => 0,
            Self::PeriodicMonodromy { .. } => 1,
            Self::SpatialDynamicsEvans { .. } => 2,
            Self::MaximalMaslovKreinEvans { .. } => 3,
            Self::MachineInstabilityCorollary { .. } => 4,
        }
    }

    const fn extension_rank(self) -> u8 {
        match self {
            Self::MachineInstabilityCorollary { .. } => 0,
            Self::ClassicalFiniteHamiltonian { .. } => 1,
            Self::PeriodicMonodromy { .. } | Self::SpatialDynamicsEvans { .. } => 2,
            Self::MaximalMaslovKreinEvans { .. } => 3,
        }
    }
}

/// One exact conclusion in the theorem family.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BridgeConclusionV1 {
    /// Equality of two distinct typed count objects after exact transforms to
    /// a common convention.
    PairwiseEquality {
        /// Left count kind.
        left: BridgeCountKindV1,
        /// Right count kind.
        right: BridgeCountKindV1,
        /// Left-to-common transform.
        left_to_common: SignedCountTransformV1,
        /// Right-to-common transform.
        right_to_common: SignedCountTransformV1,
    },
    /// Maximal equality of Maslov, Krein, and Evans counts.
    TripleEquality {
        /// Maslov-to-common transform.
        maslov_to_common: SignedCountTransformV1,
        /// Krein-to-common transform.
        krein_to_common: SignedCountTransformV1,
        /// Evans-to-common transform.
        evans_to_common: SignedCountTransformV1,
    },
    /// Explicit bridge from one spectral count to physical instability.
    SpectralToMachine {
        /// Spectral count being interpreted.
        spectral: BridgeCountKindV1,
        /// Spectral-to-machine count transform.
        spectral_to_machine: SignedCountTransformV1,
        /// Exact interpretation map, kept separate from convention changes.
        interpretation: BridgeCorrespondenceMapIdV1,
    },
}

impl BridgeConclusionV1 {
    /// Typed count objects participating in this conclusion.
    #[must_use]
    pub fn count_kinds(self) -> BTreeSet<BridgeCountKindV1> {
        match self {
            Self::PairwiseEquality { left, right, .. } => [left, right].into_iter().collect(),
            Self::TripleEquality { .. } => [
                BridgeCountKindV1::Maslov,
                BridgeCountKindV1::Krein,
                BridgeCountKindV1::Evans,
            ]
            .into_iter()
            .collect(),
            Self::SpectralToMachine { spectral, .. } => {
                [spectral, BridgeCountKindV1::MachineInstability]
                    .into_iter()
                    .collect()
            }
        }
    }
}

/// Atomic theorem hypotheses. Nodes list their exact product set; validation
/// never infers a missing atom from prose.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BridgeHypothesisKindV1 {
    /// Continuous symplectic path on the stated form.
    SymplecticPath,
    /// Fredholm Lagrangian pair with the stated reference.
    LagrangianFredholmPair,
    /// Regular crossings or an exact multiplicity-preserving resolution.
    CrossingFormsControlled,
    /// Endpoint convention and endpoint crossings are resolved.
    EndpointConventionResolved,
    /// Continuous Pontryagin/Krein operator path.
    PontryaginPathContinuous,
    /// Supporting Krein form is nondegenerate on the stated space.
    KreinFormNondegenerate,
    /// Neutral Krein directions are excluded or explicitly resolved.
    NeutralSignaturePolicyClosed,
    /// Analytic Fredholm/Evans family on the stated domain.
    AnalyticFredholmFamily,
    /// Essential spectrum is excluded from the theorem domain and contour.
    EssentialSpectrumExcluded,
    /// The contour is closed, oriented, and avoids point spectrum on its edge.
    ContourAdmissible,
    /// Evans normalization and analytic branch are fixed.
    EvansNormalizationFixed,
    /// Periodic orbit, phase, section, and monodromy correspondence are exact.
    PeriodicMonodromyCorrespondence,
    /// Stable/unstable spatial dichotomies and bundles are admitted.
    SpatialDichotomy,
    /// Parameter direction is exact and shared by all compared counts.
    ParameterDirectionAligned,
    /// Algebraic/geometric multiplicity interpretation is preserved.
    MultiplicityPreserved,
    /// Every pairwise correspondence map in a maximal node is exact.
    CorrespondenceMapsExact,
    /// Spectral crossing has an exact physical-instability interpretation.
    MachineInterpretationExact,
}

impl BridgeHypothesisKindV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::SymplecticPath => 0,
            Self::LagrangianFredholmPair => 1,
            Self::CrossingFormsControlled => 2,
            Self::EndpointConventionResolved => 3,
            Self::PontryaginPathContinuous => 4,
            Self::KreinFormNondegenerate => 5,
            Self::NeutralSignaturePolicyClosed => 6,
            Self::AnalyticFredholmFamily => 7,
            Self::EssentialSpectrumExcluded => 8,
            Self::ContourAdmissible => 9,
            Self::EvansNormalizationFixed => 10,
            Self::PeriodicMonodromyCorrespondence => 11,
            Self::SpatialDichotomy => 12,
            Self::ParameterDirectionAligned => 13,
            Self::MultiplicityPreserved => 14,
            Self::CorrespondenceMapsExact => 15,
            Self::MachineInterpretationExact => 16,
        }
    }
}

/// Evidence state attached to one required hypothesis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BridgeHypothesisStateV1 {
    /// Opaque witness reference. Validation checks binding/canonical shape, not
    /// scientific correctness.
    WitnessReferenced {
        /// Retained witness.
        witness: BridgeHypothesisWitnessIdV1,
    },
    /// Requirement is intentionally retained but unresolved.
    Unresolved {
        /// Explicit no-claim artifact.
        no_claim: BridgeNoClaimIdV1,
    },
    /// Requirement is falsified by a retained counterexample.
    Refuted {
        /// Counterexample artifact.
        counterexample: BridgeCounterexampleIdV1,
    },
}

/// One typed hypothesis and its current evidence state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BridgeHypothesisV1 {
    /// Hypothesis kind.
    pub kind: BridgeHypothesisKindV1,
    /// Referenced, unresolved, or refuted state.
    pub state: BridgeHypothesisStateV1,
}

/// Required adversarial falsifier categories.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BridgeFalsifierKindV1 {
    /// Multiple or degenerate crossing form.
    DegenerateCrossing,
    /// Tangential crossing under parameter mutation.
    TangentialCrossing,
    /// Crossing placed exactly at a path endpoint.
    EndpointCrossing,
    /// Neutral Krein signature.
    NeutralKreinSignature,
    /// Failure of Fredholmness.
    NonFredholmFamily,
    /// Evans contour contact/pinch.
    ContourContact,
    /// Essential-spectrum contact with the admitted domain.
    EssentialSpectrumContact,
    /// Reversed orientation or parameter direction.
    OrientationMutation,
    /// Algebraic/geometric multiplicity mutation.
    MultiplicityMutation,
    /// Reviewer mutation of one statement field or convention.
    StatementMutation,
}

impl BridgeFalsifierKindV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::DegenerateCrossing => 0,
            Self::TangentialCrossing => 1,
            Self::EndpointCrossing => 2,
            Self::NeutralKreinSignature => 3,
            Self::NonFredholmFamily => 4,
            Self::ContourContact => 5,
            Self::EssentialSpectrumContact => 6,
            Self::OrientationMutation => 7,
            Self::MultiplicityMutation => 8,
            Self::StatementMutation => 9,
        }
    }
}

/// Expected response when a preregistered falsifier fires.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BridgeFalsifierResponseV1 {
    /// The affected theorem node must refuse admission.
    RefuseNode,
    /// The node remains visible but must demote to conjecture/no-claim.
    DemoteToNoClaim,
    /// The node must become explicitly refuted.
    MarkRefuted,
}

/// One preregistered falsifier fixture/oracle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BridgeFalsifierV1 {
    /// Falsifier category.
    pub kind: BridgeFalsifierKindV1,
    /// Retained fixture/oracle artifact.
    pub artifact: BridgeFalsifierArtifactIdV1,
    /// Required response.
    pub response: BridgeFalsifierResponseV1,
}

/// Proof state of a theorem node. Even the referenced-proof variant remains
/// epistemically unauthoritative until a later proof-admission layer exists.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BridgeProofStateV1 {
    /// Bold statement retained without a proof artifact.
    StatementOnly {
        /// Exact no-proof/no-claim artifact.
        no_claim: BridgeNoClaimIdV1,
    },
    /// External or native proof artifact is referenced but not verified here.
    ProofArtifactReferenced {
        /// Proof artifact.
        artifact: BridgeProofArtifactIdV1,
        /// Checker/verifier identity.
        verifier: BridgeVerifierIdV1,
        /// Formal-system/checker version.
        formal_system: BridgeFormalSystemVersionIdV1,
        /// Explicit boundary that reference admission is not theorem truth.
        no_claim: BridgeNoClaimIdV1,
    },
    /// Node is refuted by a retained counterexample.
    Refuted {
        /// Counterexample artifact.
        counterexample: BridgeCounterexampleIdV1,
    },
}

/// One theorem node in the implication lattice.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BridgeTheoremNodeV1 {
    /// Stable node identity chosen by the statement author.
    pub id: BridgeTheoremNodeIdV1,
    /// Finite, periodic, spatial/Evans, maximal, or machine scope.
    pub scope: BridgeTheoremScopeV1,
    /// Exact mathematical domain.
    pub domain: BridgeDomainIdV1,
    /// Exact versioned operator family.
    pub operator_family: BridgeOperatorFamilyIdV1,
    /// Exact parameter path and direction.
    pub parameterization: BridgeParameterizationIdV1,
    /// Exact count equality/correspondence statement.
    pub conclusion: BridgeConclusionV1,
    /// Product set of typed hypotheses and their evidence states.
    pub hypotheses: Vec<BridgeHypothesisV1>,
    /// Preregistered adversarial falsifiers.
    pub falsifiers: Vec<BridgeFalsifierV1>,
    /// Current proof/reference/refutation state.
    pub proof: BridgeProofStateV1,
}

/// Evidence state for one implication edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BridgeImplicationStateV1 {
    /// Projection/implication witness is referenced but not verified here.
    WitnessReferenced {
        /// Witness artifact.
        witness: BridgeHypothesisWitnessIdV1,
    },
    /// Implication is retained as an unresolved conjecture.
    Unresolved {
        /// Explicit no-claim artifact.
        no_claim: BridgeNoClaimIdV1,
    },
    /// Implication has a retained counterexample.
    Refuted {
        /// Counterexample artifact.
        counterexample: BridgeCounterexampleIdV1,
    },
}

/// Directed edge from a stronger/bolder node to a weaker visible theorem.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BridgeImplicationV1 {
    /// Stronger theorem node.
    pub stronger: BridgeTheoremNodeIdV1,
    /// Weaker theorem node retained by the lattice.
    pub weaker: BridgeTheoremNodeIdV1,
    /// Exact domain/count projection map.
    pub projection: BridgeCorrespondenceMapIdV1,
    /// Referenced, unresolved, or refuted implication state.
    pub state: BridgeImplicationStateV1,
}

/// Trusted-computing-base declaration for statement/proof replay.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BridgeTrustedComputingBaseV1 {
    /// Exact TCB component-set identity.
    pub tcb: BridgeTcbIdV1,
    /// Checker/verifier implementation.
    pub verifier: BridgeVerifierIdV1,
    /// Proof-review/admission policy.
    pub policy: BridgeProofPolicyIdV1,
    /// Formal/native checker version.
    pub formal_system: BridgeFormalSystemVersionIdV1,
    /// Explicit boundary on what this TCB declaration proves.
    pub no_claim: BridgeNoClaimIdV1,
}

/// Caller-declared resource budget, enforced beneath hard schema caps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BridgeValidationBudgetV1 {
    /// Maximum theorem nodes admitted for this artifact.
    pub max_nodes: u32,
    /// Maximum implications admitted for this artifact.
    pub max_implications: u32,
    /// Maximum hypotheses per node.
    pub max_hypotheses_per_node: u32,
    /// Maximum falsifiers per node.
    pub max_falsifiers_per_node: u32,
}

/// Raw theorem-family and hypothesis-lattice descriptor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BridgeLatticeSpecV1 {
    schema_version: u32,
    statement_version: BridgeStatementVersionIdV1,
    maslov: MaslovIndexObjectV1,
    krein: KreinSpectralFlowObjectV1,
    evans: EvansWindingObjectV1,
    machine: Option<MachineInstabilityCountObjectV1>,
    nodes: Vec<BridgeTheoremNodeV1>,
    implications: Vec<BridgeImplicationV1>,
    tcb: BridgeTrustedComputingBaseV1,
    budget: BridgeValidationBudgetV1,
}

impl BridgeLatticeSpecV1 {
    /// Construct a current-version bridge statement.
    #[allow(clippy::too_many_arguments)] // Independent theorem axes remain explicit.
    #[must_use]
    pub fn new(
        statement_version: BridgeStatementVersionIdV1,
        maslov: MaslovIndexObjectV1,
        krein: KreinSpectralFlowObjectV1,
        evans: EvansWindingObjectV1,
        machine: Option<MachineInstabilityCountObjectV1>,
        nodes: Vec<BridgeTheoremNodeV1>,
        implications: Vec<BridgeImplicationV1>,
        tcb: BridgeTrustedComputingBaseV1,
        budget: BridgeValidationBudgetV1,
    ) -> Self {
        Self::with_schema_version(
            BRIDGE_LATTICE_SCHEMA_VERSION_V1,
            statement_version,
            maslov,
            krein,
            evans,
            machine,
            nodes,
            implications,
            tcb,
            budget,
        )
    }

    /// Construct decoded versioned input. Unsupported versions fail closed.
    #[allow(clippy::too_many_arguments)] // Decoder surface mirrors the wire product.
    #[must_use]
    pub fn with_schema_version(
        schema_version: u32,
        statement_version: BridgeStatementVersionIdV1,
        maslov: MaslovIndexObjectV1,
        krein: KreinSpectralFlowObjectV1,
        evans: EvansWindingObjectV1,
        machine: Option<MachineInstabilityCountObjectV1>,
        nodes: Vec<BridgeTheoremNodeV1>,
        implications: Vec<BridgeImplicationV1>,
        tcb: BridgeTrustedComputingBaseV1,
        budget: BridgeValidationBudgetV1,
    ) -> Self {
        Self {
            schema_version,
            statement_version,
            maslov,
            krein,
            evans,
            machine,
            nodes,
            implications,
            tcb,
            budget,
        }
    }

    /// Declared wire schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Version identity of the theorem statement vocabulary.
    #[must_use]
    pub const fn statement_version(&self) -> BridgeStatementVersionIdV1 {
        self.statement_version
    }

    /// Typed Maslov object.
    #[must_use]
    pub const fn maslov(&self) -> MaslovIndexObjectV1 {
        self.maslov
    }

    /// Typed Krein object.
    #[must_use]
    pub const fn krein(&self) -> KreinSpectralFlowObjectV1 {
        self.krein
    }

    /// Typed Evans object.
    #[must_use]
    pub const fn evans(&self) -> EvansWindingObjectV1 {
        self.evans
    }

    /// Optional, still-distinct physical instability count.
    #[must_use]
    pub const fn machine(&self) -> Option<MachineInstabilityCountObjectV1> {
        self.machine
    }

    /// Raw theorem nodes.
    #[must_use]
    pub fn nodes(&self) -> &[BridgeTheoremNodeV1] {
        &self.nodes
    }

    /// Raw implication edges.
    #[must_use]
    pub fn implications(&self) -> &[BridgeImplicationV1] {
        &self.implications
    }

    /// TCB declaration.
    #[must_use]
    pub const fn tcb(&self) -> BridgeTrustedComputingBaseV1 {
        self.tcb
    }

    /// Resource budget.
    #[must_use]
    pub const fn budget(&self) -> BridgeValidationBudgetV1 {
        self.budget
    }
}

/// Node disposition derived without assigning theorem truth.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BridgeNodeDispositionV1 {
    /// At least one hypothesis remains unresolved, or the node is
    /// intentionally statement-only.
    ConjectureOnly,
    /// Every hypothesis is referenced and a proof artifact is referenced, but
    /// neither was scientifically verified by this module.
    ReferencedNotVerified,
    /// A hypothesis or proof state carries a retained counterexample.
    Refuted,
}

/// Node disposition bound to its exact node identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BridgeNodeDispositionRecordV1 {
    /// Theorem node.
    pub node: BridgeTheoremNodeIdV1,
    /// Conservative local disposition.
    pub disposition: BridgeNodeDispositionV1,
}

/// Scientific authority exposed by this statement validator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BridgeScientificAuthorityV1 {
    /// Shape, references, and deterministic identity are validated; theorem
    /// correctness is not.
    ScientificCorrectnessNotProven,
}

/// Structured fail-closed validation issue.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BridgeValidationIssueV1 {
    /// Validation was cancelled at a bounded checkpoint.
    Cancelled,
    /// Unsupported wire schema.
    UnsupportedSchemaVersion {
        /// Supplied version.
        found: u32,
        /// Sole supported version.
        supported: u32,
    },
    /// A caller budget is zero or exceeds a hard schema cap.
    InvalidBudget,
    /// Node collection exceeds a hard cap or caller budget.
    TooManyNodes {
        /// Items supplied.
        found: usize,
        /// Effective limit.
        limit: usize,
    },
    /// Implication collection exceeds a hard cap or caller budget.
    TooManyImplications {
        /// Items supplied.
        found: usize,
        /// Effective limit.
        limit: usize,
    },
    /// One node exceeds its hypothesis cap or budget.
    TooManyHypotheses {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
        /// Items supplied.
        found: usize,
        /// Effective limit.
        limit: usize,
    },
    /// One node exceeds its falsifier cap or budget.
    TooManyFalsifiers {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
        /// Items supplied.
        found: usize,
        /// Effective limit.
        limit: usize,
    },
    /// Finite/periodic/spatial dimension is zero or otherwise inadmissible.
    InvalidScopeDimension {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
    },
    /// Required theorem-family scope is absent.
    MissingScope {
        /// Missing scope discriminant.
        scope: &'static str,
    },
    /// Duplicate theorem node identity.
    DuplicateNode {
        /// Duplicate ID.
        node: BridgeTheoremNodeIdV1,
    },
    /// Pairwise relation compares a count kind with itself.
    ReflexiveCountEquality {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
    },
    /// A machine relation names machine as its alleged spectral source.
    InvalidMachineSource {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
    },
    /// A pairwise relation mentions the physical machine count instead of
    /// using the dedicated spectral-to-machine relation.
    MachineCountNeedsExplicitRelation {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
    },
    /// A spectral-to-machine relation is stated outside the dedicated machine
    /// corollary scope.
    MachineRelationScopeMismatch {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
    },
    /// Machine scope and conclusion name different interpretation maps.
    MachineInterpretationMismatch {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
    },
    /// Machine node parameterization differs from the physical count object.
    MachineParameterizationMismatch {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
    },
    /// A machine theorem is present without a distinct machine-count object.
    MachineObjectMissing {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
    },
    /// A maximal scope does not state the triple equality.
    MaximalNodeNotTriple {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
    },
    /// A machine-corollary scope does not state a spectral-to-machine map.
    MachineNodeNotPhysical {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
    },
    /// The count object says neutral Krein type is unresolved while the node
    /// presents that hypothesis as witnessed.
    NeutralPolicyOverclaim {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
    },
    /// A referenced proof names a verifier or formal-system version outside
    /// the lattice's declared TCB.
    ProofTcbMismatch {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
    },
    /// Duplicate hypothesis kind within one product set.
    DuplicateHypothesis {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
        /// Duplicate kind.
        kind: BridgeHypothesisKindV1,
    },
    /// Required hypothesis atom is absent.
    MissingHypothesis {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
        /// Missing kind.
        kind: BridgeHypothesisKindV1,
    },
    /// Duplicate falsifier category.
    DuplicateFalsifier {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
        /// Duplicate category.
        kind: BridgeFalsifierKindV1,
    },
    /// Required falsifier category is absent.
    MissingFalsifier {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
        /// Missing category.
        kind: BridgeFalsifierKindV1,
    },
    /// Implication references an absent endpoint node.
    UnknownImplicationNode,
    /// Implication points from a node to itself.
    SelfImplication {
        /// Affected node.
        node: BridgeTheoremNodeIdV1,
    },
    /// Duplicate stronger-to-weaker edge.
    DuplicateImplication,
    /// An implication points from a lower extension scope to a higher one.
    ImplicationScopeOrderMismatch {
        /// Alleged stronger node.
        stronger: BridgeTheoremNodeIdV1,
        /// Alleged weaker node.
        weaker: BridgeTheoremNodeIdV1,
    },
    /// Implication graph contains a directed cycle.
    ImplicationCycle,
    /// A non-refuted maximal theorem exposes no non-refuted weaker projection.
    MaximalNodeHasNoWeakerProjection {
        /// Affected maximal node.
        node: BridgeTheoremNodeIdV1,
    },
    /// A maximal node has no non-refuted branch of the named scope that reaches
    /// a classical finite theorem.
    MaximalProjectionCoverageMissing {
        /// Affected maximal node.
        node: BridgeTheoremNodeIdV1,
        /// Missing periodic or spatial branch.
        scope: &'static str,
    },
    /// Canonical identity construction failed.
    Identity(CanonicalError),
}

/// Deterministic refusal report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BridgeValidationReportV1 {
    issues: Vec<BridgeValidationIssueV1>,
}

impl BridgeValidationReportV1 {
    fn one(issue: BridgeValidationIssueV1) -> Self {
        Self {
            issues: vec![issue],
        }
    }

    fn new(issues: Vec<BridgeValidationIssueV1>) -> Self {
        Self { issues }
    }

    /// Deterministically ordered issues.
    #[must_use]
    pub fn issues(&self) -> &[BridgeValidationIssueV1] {
        &self.issues
    }
}

impl fmt::Display for BridgeValidationReportV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Maslov-Krein-Evans bridge lattice refused with {} issue(s)",
            self.issues.len()
        )
    }
}

impl core::error::Error for BridgeValidationReportV1 {}

/// Canonical bridge statement. This token validates schema closure and replay,
/// not theorem truth.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValidatedBridgeLatticeV1 {
    spec: BridgeLatticeSpecV1,
    dispositions: Vec<BridgeNodeDispositionRecordV1>,
    receipt: IdentityReceipt<BridgeLatticeIdV1>,
}

impl ValidatedBridgeLatticeV1 {
    /// Canonicalized observational statement view.
    #[must_use]
    pub const fn spec(&self) -> &BridgeLatticeSpecV1 {
        &self.spec
    }

    /// Canonically ordered conservative node dispositions.
    #[must_use]
    pub fn dispositions(&self) -> &[BridgeNodeDispositionRecordV1] {
        &self.dispositions
    }

    /// Deterministic semantic identity.
    #[must_use]
    pub const fn lattice_id(&self) -> BridgeLatticeIdV1 {
        self.receipt.id()
    }

    /// Identity plus canonical-preimage receipt.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<BridgeLatticeIdV1> {
        self.receipt
    }

    /// This module's fixed no-theorem authority boundary.
    #[must_use]
    pub const fn scientific_authority(&self) -> BridgeScientificAuthorityV1 {
        BridgeScientificAuthorityV1::ScientificCorrectnessNotProven
    }
}

/// Validate and canonicalize a theorem-family hypothesis lattice.
///
/// Collection caps are checked before sorting, graph construction, or identity
/// work. Cancellation is polled initially and at every bounded node/edge loop.
/// No partial token escapes on failure.
///
/// # Errors
/// Returns [`BridgeValidationReportV1`] for malformed scopes, missing
/// hypotheses/falsifiers, graph defects, cancellation, or identity failure.
#[must_use = "the theorem statement must be validated before use"]
#[allow(clippy::too_many_lines)] // The fail-closed phases remain visible in one audit trail.
pub fn validate_bridge_lattice_v1(
    mut spec: BridgeLatticeSpecV1,
    cx: &Cx<'_>,
) -> Result<ValidatedBridgeLatticeV1, BridgeValidationReportV1> {
    checkpoint(cx)?;
    let budget = spec.budget;
    if !budget_is_valid(budget) {
        return Err(BridgeValidationReportV1::one(
            BridgeValidationIssueV1::InvalidBudget,
        ));
    }
    let node_limit = usize::min(MAX_BRIDGE_NODES_V1, budget.max_nodes as usize);
    if spec.nodes.len() > node_limit {
        return Err(BridgeValidationReportV1::one(
            BridgeValidationIssueV1::TooManyNodes {
                found: spec.nodes.len(),
                limit: node_limit,
            },
        ));
    }
    let implication_limit =
        usize::min(MAX_BRIDGE_IMPLICATIONS_V1, budget.max_implications as usize);
    if spec.implications.len() > implication_limit {
        return Err(BridgeValidationReportV1::one(
            BridgeValidationIssueV1::TooManyImplications {
                found: spec.implications.len(),
                limit: implication_limit,
            },
        ));
    }
    let hypothesis_limit = usize::min(
        MAX_BRIDGE_HYPOTHESES_PER_NODE_V1,
        budget.max_hypotheses_per_node as usize,
    );
    let falsifier_limit = usize::min(
        MAX_BRIDGE_FALSIFIERS_PER_NODE_V1,
        budget.max_falsifiers_per_node as usize,
    );
    let mut oversized_hypotheses = None;
    let mut oversized_falsifiers = None;
    for node in &spec.nodes {
        if node.hypotheses.len() > hypothesis_limit
            && oversized_hypotheses
                .as_ref()
                .is_none_or(|(id, found)| (node.id, node.hypotheses.len()) < (*id, *found))
        {
            oversized_hypotheses = Some((node.id, node.hypotheses.len()));
        }
        if node.falsifiers.len() > falsifier_limit
            && oversized_falsifiers
                .as_ref()
                .is_none_or(|(id, found)| (node.id, node.falsifiers.len()) < (*id, *found))
        {
            oversized_falsifiers = Some((node.id, node.falsifiers.len()));
        }
    }
    if let Some((node, found)) = oversized_hypotheses {
        return Err(BridgeValidationReportV1::one(
            BridgeValidationIssueV1::TooManyHypotheses {
                node,
                found,
                limit: hypothesis_limit,
            },
        ));
    }
    if let Some((node, found)) = oversized_falsifiers {
        return Err(BridgeValidationReportV1::one(
            BridgeValidationIssueV1::TooManyFalsifiers {
                node,
                found,
                limit: falsifier_limit,
            },
        ));
    }

    for node in &mut spec.nodes {
        node.hypotheses.sort_by_cached_key(|item| {
            let mut encoded = Vec::with_capacity(34);
            push_hypothesis(&mut encoded, *item);
            encoded
        });
        node.falsifiers.sort_by_cached_key(|item| {
            let mut encoded = Vec::with_capacity(34);
            push_falsifier(&mut encoded, *item);
            encoded
        });
    }
    spec.nodes
        .sort_by_cached_key(|node| (*node.id.as_bytes(), node_bytes(node)));
    spec.implications.sort_by_cached_key(|implication| {
        (
            implication_sort_key(implication),
            implication_bytes(implication),
        )
    });

    let mut issues = Vec::new();
    if spec.schema_version != BRIDGE_LATTICE_SCHEMA_VERSION_V1 {
        issues.push(BridgeValidationIssueV1::UnsupportedSchemaVersion {
            found: spec.schema_version,
            supported: BRIDGE_LATTICE_SCHEMA_VERSION_V1,
        });
    }

    let mut node_ids = BTreeSet::new();
    let mut scopes = BTreeSet::new();
    for node in &spec.nodes {
        checkpoint(cx)?;
        scopes.insert(node.scope.tag());
        if !node_ids.insert(node.id) {
            issues.push(BridgeValidationIssueV1::DuplicateNode { node: node.id });
        }
        validate_scope(node, &mut issues);
        validate_conclusion(node, spec.machine, spec.krein.neutral_policy, &mut issues);
        validate_hypotheses(node, &mut issues);
        validate_falsifiers(node, &mut issues);
        validate_proof_tcb(node, spec.tcb, &mut issues);
    }
    for (tag, name) in [
        (0, "classical-finite-Hamiltonian"),
        (1, "periodic-monodromy"),
        (2, "spatial-dynamics-Evans"),
        (3, "maximal-Maslov-Krein-Evans"),
    ] {
        if !scopes.contains(&tag) {
            issues.push(BridgeValidationIssueV1::MissingScope { scope: name });
        }
    }

    let mut edges = BTreeSet::new();
    let mut adjacency: BTreeMap<BridgeTheoremNodeIdV1, Vec<BridgeTheoremNodeIdV1>> =
        BTreeMap::new();
    let mut non_refuted_adjacency: BTreeMap<BridgeTheoremNodeIdV1, Vec<BridgeTheoremNodeIdV1>> =
        BTreeMap::new();
    let scopes_by_node: BTreeMap<_, _> = spec
        .nodes
        .iter()
        .map(|node| (node.id, node.scope))
        .collect();
    let non_refuted_scopes: BTreeMap<_, _> = spec
        .nodes
        .iter()
        .filter(|node| !node_is_refuted(node))
        .map(|node| (node.id, node.scope))
        .collect();
    for implication in &spec.implications {
        checkpoint(cx)?;
        if !node_ids.contains(&implication.stronger) || !node_ids.contains(&implication.weaker) {
            issues.push(BridgeValidationIssueV1::UnknownImplicationNode);
            continue;
        }
        if implication.stronger == implication.weaker {
            issues.push(BridgeValidationIssueV1::SelfImplication {
                node: implication.stronger,
            });
            continue;
        }
        if !edges.insert((implication.stronger, implication.weaker)) {
            issues.push(BridgeValidationIssueV1::DuplicateImplication);
            continue;
        }
        let stronger_scope = scopes_by_node[&implication.stronger];
        let weaker_scope = scopes_by_node[&implication.weaker];
        let scope_order_valid = stronger_scope.extension_rank() >= weaker_scope.extension_rank();
        if !scope_order_valid {
            issues.push(BridgeValidationIssueV1::ImplicationScopeOrderMismatch {
                stronger: implication.stronger,
                weaker: implication.weaker,
            });
        }
        adjacency
            .entry(implication.stronger)
            .or_default()
            .push(implication.weaker);
        if scope_order_valid
            && !matches!(implication.state, BridgeImplicationStateV1::Refuted { .. })
            && non_refuted_scopes.contains_key(&implication.stronger)
            && non_refuted_scopes.contains_key(&implication.weaker)
        {
            non_refuted_adjacency
                .entry(implication.stronger)
                .or_default()
                .push(implication.weaker);
        }
    }
    if implication_cycle(&node_ids, &adjacency, cx)? {
        issues.push(BridgeValidationIssueV1::ImplicationCycle);
    }
    let reaches_classical =
        nodes_reaching_classical(&non_refuted_scopes, &non_refuted_adjacency, cx)?;
    for node in &spec.nodes {
        let is_non_refuted_maximal = matches!(
            node.scope,
            BridgeTheoremScopeV1::MaximalMaslovKreinEvans { .. }
        ) && !node_is_refuted(node);
        if is_non_refuted_maximal
            && non_refuted_adjacency
                .get(&node.id)
                .is_none_or(Vec::is_empty)
        {
            issues
                .push(BridgeValidationIssueV1::MaximalNodeHasNoWeakerProjection { node: node.id });
        }
        if is_non_refuted_maximal {
            let direct_projections = non_refuted_adjacency
                .get(&node.id)
                .map_or(&[][..], Vec::as_slice);
            for (branch_tag, scope) in [(1, "periodic-monodromy"), (2, "spatial-dynamics-Evans")] {
                if !direct_projections.iter().any(|projection| {
                    reaches_classical.contains(projection)
                        && scopes_by_node
                            .get(projection)
                            .is_some_and(|scope| scope.tag() == branch_tag)
                }) {
                    issues.push(BridgeValidationIssueV1::MaximalProjectionCoverageMissing {
                        node: node.id,
                        scope,
                    });
                }
            }
        }
    }

    if !issues.is_empty() {
        return Err(BridgeValidationReportV1::new(issues));
    }
    checkpoint(cx)?;
    let receipt = bridge_lattice_receipt(&spec)
        .map_err(|error| BridgeValidationReportV1::one(BridgeValidationIssueV1::Identity(error)))?;
    let dispositions = spec
        .nodes
        .iter()
        .map(|node| BridgeNodeDispositionRecordV1 {
            node: node.id,
            disposition: node_disposition(node),
        })
        .collect();
    Ok(ValidatedBridgeLatticeV1 {
        spec,
        dispositions,
        receipt,
    })
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), BridgeValidationReportV1> {
    cx.checkpoint()
        .map_err(|_| BridgeValidationReportV1::one(BridgeValidationIssueV1::Cancelled))
}

fn budget_is_valid(budget: BridgeValidationBudgetV1) -> bool {
    budget.max_nodes > 0
        && budget.max_nodes as usize <= MAX_BRIDGE_NODES_V1
        && budget.max_implications > 0
        && budget.max_implications as usize <= MAX_BRIDGE_IMPLICATIONS_V1
        && budget.max_hypotheses_per_node > 0
        && budget.max_hypotheses_per_node as usize <= MAX_BRIDGE_HYPOTHESES_PER_NODE_V1
        && budget.max_falsifiers_per_node > 0
        && budget.max_falsifiers_per_node as usize <= MAX_BRIDGE_FALSIFIERS_PER_NODE_V1
}

fn validate_scope(node: &BridgeTheoremNodeV1, issues: &mut Vec<BridgeValidationIssueV1>) {
    let invalid = match node.scope {
        BridgeTheoremScopeV1::ClassicalFiniteHamiltonian { phase_dimension } => {
            phase_dimension == 0 || phase_dimension % 2 != 0
        }
        BridgeTheoremScopeV1::PeriodicMonodromy {
            state_dimension, ..
        } => {
            state_dimension == 0
                || (node
                    .conclusion
                    .count_kinds()
                    .contains(&BridgeCountKindV1::Maslov)
                    && state_dimension % 2 != 0)
        }
        BridgeTheoremScopeV1::SpatialDynamicsEvans {
            stable_dimension,
            unstable_dimension,
            ..
        } => stable_dimension == 0 || unstable_dimension == 0,
        BridgeTheoremScopeV1::MaximalMaslovKreinEvans { .. }
        | BridgeTheoremScopeV1::MachineInstabilityCorollary { .. } => false,
    };
    if invalid {
        issues.push(BridgeValidationIssueV1::InvalidScopeDimension { node: node.id });
    }
}

fn validate_conclusion(
    node: &BridgeTheoremNodeV1,
    machine: Option<MachineInstabilityCountObjectV1>,
    neutral_policy: NeutralKreinPolicyV1,
    issues: &mut Vec<BridgeValidationIssueV1>,
) {
    match node.conclusion {
        BridgeConclusionV1::PairwiseEquality { left, right, .. } => {
            if left == right {
                issues.push(BridgeValidationIssueV1::ReflexiveCountEquality { node: node.id });
            }
            if left == BridgeCountKindV1::MachineInstability
                || right == BridgeCountKindV1::MachineInstability
            {
                issues.push(BridgeValidationIssueV1::MachineCountNeedsExplicitRelation {
                    node: node.id,
                });
            }
        }
        BridgeConclusionV1::SpectralToMachine {
            spectral,
            interpretation,
            ..
        } => {
            if spectral == BridgeCountKindV1::MachineInstability {
                issues.push(BridgeValidationIssueV1::InvalidMachineSource { node: node.id });
            }
            if machine.is_none() {
                issues.push(BridgeValidationIssueV1::MachineObjectMissing { node: node.id });
            }
            match node.scope {
                BridgeTheoremScopeV1::MachineInstabilityCorollary {
                    interpretation: scoped,
                } => {
                    if scoped != interpretation {
                        issues.push(BridgeValidationIssueV1::MachineInterpretationMismatch {
                            node: node.id,
                        });
                    }
                    if machine
                        .is_some_and(|object| object.parameterization != node.parameterization)
                    {
                        issues.push(BridgeValidationIssueV1::MachineParameterizationMismatch {
                            node: node.id,
                        });
                    }
                }
                _ => {
                    issues.push(BridgeValidationIssueV1::MachineRelationScopeMismatch {
                        node: node.id,
                    });
                }
            }
        }
        BridgeConclusionV1::TripleEquality { .. } => {}
    }
    if matches!(
        node.scope,
        BridgeTheoremScopeV1::MaximalMaslovKreinEvans { .. }
    ) && !matches!(node.conclusion, BridgeConclusionV1::TripleEquality { .. })
    {
        issues.push(BridgeValidationIssueV1::MaximalNodeNotTriple { node: node.id });
    }
    if matches!(
        node.scope,
        BridgeTheoremScopeV1::MachineInstabilityCorollary { .. }
    ) && !matches!(
        node.conclusion,
        BridgeConclusionV1::SpectralToMachine { .. }
    ) {
        issues.push(BridgeValidationIssueV1::MachineNodeNotPhysical { node: node.id });
    }
    if node
        .conclusion
        .count_kinds()
        .contains(&BridgeCountKindV1::Krein)
        && matches!(neutral_policy, NeutralKreinPolicyV1::Unresolved { .. })
        && node.hypotheses.iter().any(|hypothesis| {
            hypothesis.kind == BridgeHypothesisKindV1::NeutralSignaturePolicyClosed
                && matches!(
                    hypothesis.state,
                    BridgeHypothesisStateV1::WitnessReferenced { .. }
                )
        })
    {
        issues.push(BridgeValidationIssueV1::NeutralPolicyOverclaim { node: node.id });
    }
}

fn validate_hypotheses(node: &BridgeTheoremNodeV1, issues: &mut Vec<BridgeValidationIssueV1>) {
    let mut present = BTreeSet::new();
    for hypothesis in &node.hypotheses {
        if !present.insert(hypothesis.kind) {
            issues.push(BridgeValidationIssueV1::DuplicateHypothesis {
                node: node.id,
                kind: hypothesis.kind,
            });
        }
    }
    for kind in required_hypotheses(node) {
        if !present.contains(&kind) {
            issues.push(BridgeValidationIssueV1::MissingHypothesis {
                node: node.id,
                kind,
            });
        }
    }
}

fn validate_proof_tcb(
    node: &BridgeTheoremNodeV1,
    tcb: BridgeTrustedComputingBaseV1,
    issues: &mut Vec<BridgeValidationIssueV1>,
) {
    if let BridgeProofStateV1::ProofArtifactReferenced {
        verifier,
        formal_system,
        ..
    } = node.proof
        && (verifier != tcb.verifier || formal_system != tcb.formal_system)
    {
        issues.push(BridgeValidationIssueV1::ProofTcbMismatch { node: node.id });
    }
}

fn required_hypotheses(node: &BridgeTheoremNodeV1) -> BTreeSet<BridgeHypothesisKindV1> {
    use BridgeCountKindV1::{Evans, Krein, MachineInstability, Maslov};
    use BridgeHypothesisKindV1 as H;

    let kinds = node.conclusion.count_kinds();
    let mut required = BTreeSet::from([H::ParameterDirectionAligned, H::MultiplicityPreserved]);
    if kinds.contains(&Maslov) {
        required.extend([
            H::SymplecticPath,
            H::LagrangianFredholmPair,
            H::CrossingFormsControlled,
            H::EndpointConventionResolved,
        ]);
    }
    if kinds.contains(&Krein) {
        required.extend([
            H::PontryaginPathContinuous,
            H::KreinFormNondegenerate,
            H::NeutralSignaturePolicyClosed,
            H::EndpointConventionResolved,
        ]);
    }
    if kinds.contains(&Evans) {
        required.extend([
            H::AnalyticFredholmFamily,
            H::EssentialSpectrumExcluded,
            H::ContourAdmissible,
            H::EvansNormalizationFixed,
        ]);
    }
    if kinds.contains(&MachineInstability) {
        required.insert(H::MachineInterpretationExact);
    }
    match node.scope {
        BridgeTheoremScopeV1::PeriodicMonodromy { .. } => {
            required.insert(H::PeriodicMonodromyCorrespondence);
        }
        BridgeTheoremScopeV1::SpatialDynamicsEvans { .. } => {
            required.insert(H::SpatialDichotomy);
        }
        BridgeTheoremScopeV1::MaximalMaslovKreinEvans { .. } => {
            required.extend([
                H::PeriodicMonodromyCorrespondence,
                H::SpatialDichotomy,
                H::CorrespondenceMapsExact,
            ]);
        }
        BridgeTheoremScopeV1::ClassicalFiniteHamiltonian { .. }
        | BridgeTheoremScopeV1::MachineInstabilityCorollary { .. } => {}
    }
    required
}

fn validate_falsifiers(node: &BridgeTheoremNodeV1, issues: &mut Vec<BridgeValidationIssueV1>) {
    let mut present = BTreeSet::new();
    for falsifier in &node.falsifiers {
        if !present.insert(falsifier.kind) {
            issues.push(BridgeValidationIssueV1::DuplicateFalsifier {
                node: node.id,
                kind: falsifier.kind,
            });
        }
    }
    for kind in required_falsifiers(node.conclusion) {
        if !present.contains(&kind) {
            issues.push(BridgeValidationIssueV1::MissingFalsifier {
                node: node.id,
                kind,
            });
        }
    }
}

fn required_falsifiers(conclusion: BridgeConclusionV1) -> BTreeSet<BridgeFalsifierKindV1> {
    use BridgeCountKindV1::{Evans, Krein, Maslov};
    use BridgeFalsifierKindV1 as F;

    let kinds = conclusion.count_kinds();
    let mut required = BTreeSet::from([
        F::OrientationMutation,
        F::MultiplicityMutation,
        F::StatementMutation,
    ]);
    if kinds.contains(&Maslov) || kinds.contains(&Krein) {
        required.extend([
            F::DegenerateCrossing,
            F::TangentialCrossing,
            F::EndpointCrossing,
        ]);
    }
    if kinds.contains(&Krein) {
        required.insert(F::NeutralKreinSignature);
    }
    if kinds.contains(&Evans) {
        required.extend([
            F::NonFredholmFamily,
            F::ContourContact,
            F::EssentialSpectrumContact,
        ]);
    }
    required
}

fn implication_cycle(
    nodes: &BTreeSet<BridgeTheoremNodeIdV1>,
    adjacency: &BTreeMap<BridgeTheoremNodeIdV1, Vec<BridgeTheoremNodeIdV1>>,
    cx: &Cx<'_>,
) -> Result<bool, BridgeValidationReportV1> {
    let mut indegree: BTreeMap<BridgeTheoremNodeIdV1, usize> =
        nodes.iter().copied().map(|node| (node, 0)).collect();
    for weaker_nodes in adjacency.values() {
        for weaker in weaker_nodes {
            *indegree.entry(*weaker).or_default() += 1;
        }
    }
    let mut ready: BTreeSet<_> = indegree
        .iter()
        .filter_map(|(node, degree)| (*degree == 0).then_some(*node))
        .collect();
    let mut visited = 0_usize;
    while let Some(node) = ready.pop_first() {
        checkpoint(cx)?;
        visited += 1;
        if let Some(weaker_nodes) = adjacency.get(&node) {
            for weaker in weaker_nodes {
                let Some(degree) = indegree.get_mut(weaker) else {
                    // The caller already records unknown endpoints. Preserve
                    // panic freedom if this internal invariant is ever
                    // changed without the graph builder being updated.
                    return Ok(true);
                };
                *degree -= 1;
                if *degree == 0 {
                    ready.insert(*weaker);
                }
            }
        }
    }
    Ok(visited != nodes.len())
}

fn nodes_reaching_classical(
    scopes: &BTreeMap<BridgeTheoremNodeIdV1, BridgeTheoremScopeV1>,
    adjacency: &BTreeMap<BridgeTheoremNodeIdV1, Vec<BridgeTheoremNodeIdV1>>,
    cx: &Cx<'_>,
) -> Result<BTreeSet<BridgeTheoremNodeIdV1>, BridgeValidationReportV1> {
    let mut reverse: BTreeMap<_, Vec<_>> = BTreeMap::new();
    for (stronger, weaker_nodes) in adjacency {
        for weaker in weaker_nodes {
            reverse.entry(*weaker).or_default().push(*stronger);
        }
    }
    let mut reaches_classical = BTreeSet::new();
    let mut pending: Vec<_> = scopes
        .iter()
        .filter_map(|(node, scope)| (scope.tag() == 0).then_some(*node))
        .collect();
    while let Some(node) = pending.pop() {
        checkpoint(cx)?;
        if !reaches_classical.insert(node) {
            continue;
        }
        if let Some(parents) = reverse.get(&node) {
            pending.extend(parents.iter().copied());
        }
    }
    Ok(reaches_classical)
}

fn node_is_refuted(node: &BridgeTheoremNodeV1) -> bool {
    matches!(node.proof, BridgeProofStateV1::Refuted { .. })
        || node
            .hypotheses
            .iter()
            .any(|item| matches!(item.state, BridgeHypothesisStateV1::Refuted { .. }))
}

fn node_disposition(node: &BridgeTheoremNodeV1) -> BridgeNodeDispositionV1 {
    if node_is_refuted(node) {
        return BridgeNodeDispositionV1::Refuted;
    }
    let all_hypotheses_referenced = node.hypotheses.iter().all(|item| {
        matches!(
            item.state,
            BridgeHypothesisStateV1::WitnessReferenced { .. }
        )
    });
    if all_hypotheses_referenced
        && matches!(
            node.proof,
            BridgeProofStateV1::ProofArtifactReferenced { .. }
        )
    {
        BridgeNodeDispositionV1::ReferencedNotVerified
    } else {
        BridgeNodeDispositionV1::ConjectureOnly
    }
}

fn implication_sort_key(implication: &BridgeImplicationV1) -> ([u8; 32], [u8; 32], [u8; 32]) {
    (
        *implication.stronger.as_bytes(),
        *implication.weaker.as_bytes(),
        *implication.projection.as_bytes(),
    )
}

fn push_id<I: DigestBytes + Copy>(out: &mut Vec<u8>, id: I) {
    out.extend_from_slice(id.digest_bytes());
}

fn push_u32(out: &mut Vec<u8>, value: u32) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn push_i64(out: &mut Vec<u8>, value: i64) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn push_convention(out: &mut Vec<u8>, convention: SignedCountConventionV1) {
    out.push(convention.orientation.tag());
    out.push(convention.endpoints.left.tag());
    out.push(convention.endpoints.right.tag());
}

fn push_transform(out: &mut Vec<u8>, transform: SignedCountTransformV1) {
    push_id(out, transform.map);
    out.push(transform.sign as u8);
    push_i64(out, transform.doubled_offset);
}

fn count_objects_bytes(spec: &BridgeLatticeSpecV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(1024);
    push_id(&mut out, spec.maslov.path);
    push_id(&mut out, spec.maslov.reference);
    push_id(&mut out, spec.maslov.symplectic_form);
    push_id(&mut out, spec.maslov.crossing_forms);
    push_convention(&mut out, spec.maslov.convention);
    push_id(&mut out, spec.maslov.multiplicity);

    push_id(&mut out, spec.krein.path);
    push_id(&mut out, spec.krein.krein_form);
    match spec.krein.neutral_policy {
        NeutralKreinPolicyV1::Excluded { witness } => {
            out.push(0);
            push_id(&mut out, witness);
        }
        NeutralKreinPolicyV1::Resolved { map, no_claim } => {
            out.push(1);
            push_id(&mut out, map);
            push_id(&mut out, no_claim);
        }
        NeutralKreinPolicyV1::Unresolved { no_claim } => {
            out.push(2);
            push_id(&mut out, no_claim);
        }
    }
    push_convention(&mut out, spec.krein.convention);
    push_id(&mut out, spec.krein.multiplicity);

    push_id(&mut out, spec.evans.family);
    push_id(&mut out, spec.evans.analytic_domain);
    push_id(&mut out, spec.evans.contour);
    push_id(&mut out, spec.evans.essential_spectrum_exclusion);
    push_id(&mut out, spec.evans.normalization);
    push_convention(&mut out, spec.evans.convention);
    push_id(&mut out, spec.evans.multiplicity);

    match spec.machine {
        Some(machine) => {
            out.push(1);
            push_id(&mut out, machine.model);
            push_id(&mut out, machine.parameterization);
            push_convention(&mut out, machine.convention);
            push_id(&mut out, machine.multiplicity);
        }
        None => out.push(0),
    }
    out
}

fn push_scope(out: &mut Vec<u8>, scope: BridgeTheoremScopeV1) {
    out.push(scope.tag());
    match scope {
        BridgeTheoremScopeV1::ClassicalFiniteHamiltonian { phase_dimension } => {
            push_u32(out, phase_dimension);
        }
        BridgeTheoremScopeV1::PeriodicMonodromy {
            state_dimension,
            monodromy_map,
        } => {
            push_u32(out, state_dimension);
            push_id(out, monodromy_map);
        }
        BridgeTheoremScopeV1::SpatialDynamicsEvans {
            stable_dimension,
            unstable_dimension,
            spatial_map,
        } => {
            push_u32(out, stable_dimension);
            push_u32(out, unstable_dimension);
            push_id(out, spatial_map);
        }
        BridgeTheoremScopeV1::MaximalMaslovKreinEvans {
            finite_to_periodic,
            periodic_to_evans,
        } => {
            push_id(out, finite_to_periodic);
            push_id(out, periodic_to_evans);
        }
        BridgeTheoremScopeV1::MachineInstabilityCorollary { interpretation } => {
            push_id(out, interpretation);
        }
    }
}

fn push_conclusion(out: &mut Vec<u8>, conclusion: BridgeConclusionV1) {
    match conclusion {
        BridgeConclusionV1::PairwiseEquality {
            left,
            right,
            left_to_common,
            right_to_common,
        } => {
            out.push(0);
            out.push(left.tag());
            out.push(right.tag());
            push_transform(out, left_to_common);
            push_transform(out, right_to_common);
        }
        BridgeConclusionV1::TripleEquality {
            maslov_to_common,
            krein_to_common,
            evans_to_common,
        } => {
            out.push(1);
            push_transform(out, maslov_to_common);
            push_transform(out, krein_to_common);
            push_transform(out, evans_to_common);
        }
        BridgeConclusionV1::SpectralToMachine {
            spectral,
            spectral_to_machine,
            interpretation,
        } => {
            out.push(2);
            out.push(spectral.tag());
            push_transform(out, spectral_to_machine);
            push_id(out, interpretation);
        }
    }
}

fn push_hypothesis(out: &mut Vec<u8>, hypothesis: BridgeHypothesisV1) {
    out.push(hypothesis.kind.tag());
    match hypothesis.state {
        BridgeHypothesisStateV1::WitnessReferenced { witness } => {
            out.push(0);
            push_id(out, witness);
        }
        BridgeHypothesisStateV1::Unresolved { no_claim } => {
            out.push(1);
            push_id(out, no_claim);
        }
        BridgeHypothesisStateV1::Refuted { counterexample } => {
            out.push(2);
            push_id(out, counterexample);
        }
    }
}

fn push_falsifier(out: &mut Vec<u8>, falsifier: BridgeFalsifierV1) {
    out.push(falsifier.kind.tag());
    push_id(out, falsifier.artifact);
    out.push(match falsifier.response {
        BridgeFalsifierResponseV1::RefuseNode => 0,
        BridgeFalsifierResponseV1::DemoteToNoClaim => 1,
        BridgeFalsifierResponseV1::MarkRefuted => 2,
    });
}

fn node_bytes(node: &BridgeTheoremNodeV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(2048);
    push_id(&mut out, node.id);
    push_scope(&mut out, node.scope);
    push_id(&mut out, node.domain);
    push_id(&mut out, node.operator_family);
    push_id(&mut out, node.parameterization);
    push_conclusion(&mut out, node.conclusion);
    push_u32(&mut out, node.hypotheses.len() as u32);
    for hypothesis in &node.hypotheses {
        push_hypothesis(&mut out, *hypothesis);
    }
    push_u32(&mut out, node.falsifiers.len() as u32);
    for falsifier in &node.falsifiers {
        push_falsifier(&mut out, *falsifier);
    }
    match node.proof {
        BridgeProofStateV1::StatementOnly { no_claim } => {
            out.push(0);
            push_id(&mut out, no_claim);
        }
        BridgeProofStateV1::ProofArtifactReferenced {
            artifact,
            verifier,
            formal_system,
            no_claim,
        } => {
            out.push(1);
            push_id(&mut out, artifact);
            push_id(&mut out, verifier);
            push_id(&mut out, formal_system);
            push_id(&mut out, no_claim);
        }
        BridgeProofStateV1::Refuted { counterexample } => {
            out.push(2);
            push_id(&mut out, counterexample);
        }
    }
    out
}

fn implication_bytes(implication: &BridgeImplicationV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(192);
    push_id(&mut out, implication.stronger);
    push_id(&mut out, implication.weaker);
    push_id(&mut out, implication.projection);
    match implication.state {
        BridgeImplicationStateV1::WitnessReferenced { witness } => {
            out.push(0);
            push_id(&mut out, witness);
        }
        BridgeImplicationStateV1::Unresolved { no_claim } => {
            out.push(1);
            push_id(&mut out, no_claim);
        }
        BridgeImplicationStateV1::Refuted { counterexample } => {
            out.push(2);
            push_id(&mut out, counterexample);
        }
    }
    out
}

fn tcb_bytes(tcb: BridgeTrustedComputingBaseV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(160);
    push_id(&mut out, tcb.tcb);
    push_id(&mut out, tcb.verifier);
    push_id(&mut out, tcb.policy);
    push_id(&mut out, tcb.formal_system);
    push_id(&mut out, tcb.no_claim);
    out
}

fn budget_bytes(budget: BridgeValidationBudgetV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(16);
    push_u32(&mut out, budget.max_nodes);
    push_u32(&mut out, budget.max_implications);
    push_u32(&mut out, budget.max_hypotheses_per_node);
    push_u32(&mut out, budget.max_falsifiers_per_node);
    out
}

fn bridge_lattice_receipt(
    spec: &BridgeLatticeSpecV1,
) -> Result<IdentityReceipt<BridgeLatticeIdV1>, CanonicalError> {
    let counts = count_objects_bytes(spec);
    let nodes: Vec<Vec<u8>> = spec.nodes.iter().map(node_bytes).collect();
    let implications: Vec<Vec<u8>> = spec.implications.iter().map(implication_bytes).collect();
    let tcb = tcb_bytes(spec.tcb);
    let budget = budget_bytes(spec.budget);
    CanonicalEncoder::<BridgeLatticeIdV1, _>::new(BRIDGE_IDENTITY_LIMITS, NeverCancel)?
        .bytes(
            Field::new(0, "statement-version"),
            spec.statement_version.as_bytes(),
        )?
        .bytes(Field::new(1, "count-objects"), &counts)?
        .canonical_set(
            Field::new(2, "theorem-nodes"),
            nodes.len() as u64,
            nodes.iter().map(Vec::as_slice),
        )?
        .canonical_set(
            Field::new(3, "implications"),
            implications.len() as u64,
            implications.iter().map(Vec::as_slice),
        )?
        .bytes(Field::new(4, "tcb"), &tcb)?
        .bytes(Field::new(5, "budget"), &budget)?
        .finish()
}
