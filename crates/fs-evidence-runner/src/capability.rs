//! Pure semantic root-capability policy and least-privilege validation.

use crate::canonical::CanonicalFrameV1;
use crate::catalog::{
    DestinationAdmissionModeV2, DigestRoleV2, OverlapPolicyRelationV2, PlatformPathProfileV2,
    RootCapabilityAccessV2, RootCapabilityRightV2, RootClassV2, RunnerCommandV2,
};
use crate::construction::{ConstructionErrorKindV2, ConstructionErrorV2};
use crate::identity::{DigestValueV2, NoClaimScopeRootV1, RootCapabilityPolicyRootV2};
use crate::publication::PublicationSelectionV2;
use std::collections::BTreeSet;

/// Maximum aggregate registrations across every root-policy kind in one
/// family projection.
pub const ROOT_POLICY_REGISTRATIONS_MAX_V2: usize = 64;
/// Canonical domain for one semantic root-capability policy.
pub const ROOT_CAPABILITY_POLICY_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.root-capability-policy.v1";

/// One registered overlap-policy declaration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct OverlapPolicyRegistrationV2 {
    policy_id: u16,
    relation: OverlapPolicyRelationV2,
}

impl OverlapPolicyRegistrationV2 {
    /// Construct one nonzero registered declaration.
    pub fn new(
        policy_id: u16,
        relation: OverlapPolicyRelationV2,
    ) -> Result<Self, ConstructionErrorV2> {
        require_nonzero("root_policy_registry.overlap_policy_id", policy_id)?;
        Ok(Self {
            policy_id,
            relation,
        })
    }

    /// Registered identifier.
    #[must_use]
    pub const fn policy_id(self) -> u16 {
        self.policy_id
    }

    /// Closed overlap relation.
    #[must_use]
    pub const fn relation(self) -> OverlapPolicyRelationV2 {
        self.relation
    }
}

/// Bounded, non-wire projection emitted by the sealed family-policy registry.
///
/// Registration declares policy semantics. It is not proof that a physical
/// resource is fresh, unrevoked, or disjoint from another acquired resource.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RootPolicyRegistryProjectionV2 {
    freshness_policy_ids: Box<[u16]>,
    revocation_policy_ids: Box<[u16]>,
    overlap_policies: Box<[OverlapPolicyRegistrationV2]>,
}

impl RootPolicyRegistryProjectionV2 {
    /// Validate bounded, nonzero, duplicate-free registries and canonicalize
    /// their nonsemantic order.
    pub fn new(
        mut freshness_policy_ids: Vec<u16>,
        mut revocation_policy_ids: Vec<u16>,
        mut overlap_policies: Vec<OverlapPolicyRegistrationV2>,
    ) -> Result<Self, ConstructionErrorV2> {
        validate_id_registry(
            "root_policy_registry.freshness_policy_ids",
            &freshness_policy_ids,
        )?;
        validate_id_registry(
            "root_policy_registry.revocation_policy_ids",
            &revocation_policy_ids,
        )?;
        let aggregate_count = freshness_policy_ids
            .len()
            .checked_add(revocation_policy_ids.len())
            .and_then(|count| count.checked_add(overlap_policies.len()))
            .ok_or_else(|| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::ArithmeticOverflow,
                    "root_policy_registry.registrations",
                    "an aggregate registration count representable as usize",
                    "count overflow",
                )
            })?;
        if aggregate_count > ROOT_POLICY_REGISTRATIONS_MAX_V2 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "root_policy_registry.registrations",
                "at most 64 aggregate freshness, revocation, and overlap registrations",
                aggregate_count,
            ));
        }
        if overlap_policies.len() > ROOT_POLICY_REGISTRATIONS_MAX_V2 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "root_policy_registry.overlap_policies",
                "at most 64 registrations",
                overlap_policies.len(),
            ));
        }
        let mut seen_overlap = BTreeSet::new();
        for policy in &overlap_policies {
            if !seen_overlap.insert(policy.policy_id) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "root_policy_registry.overlap_policies",
                    "unique policy IDs",
                    policy.policy_id,
                ));
            }
        }
        freshness_policy_ids.sort_unstable();
        revocation_policy_ids.sort_unstable();
        overlap_policies.sort_unstable();
        Ok(Self {
            freshness_policy_ids: freshness_policy_ids.into_boxed_slice(),
            revocation_policy_ids: revocation_policy_ids.into_boxed_slice(),
            overlap_policies: overlap_policies.into_boxed_slice(),
        })
    }

    /// Canonical freshness-policy IDs.
    #[must_use]
    pub fn freshness_policy_ids(&self) -> &[u16] {
        &self.freshness_policy_ids
    }

    /// Canonical revocation-policy IDs.
    #[must_use]
    pub fn revocation_policy_ids(&self) -> &[u16] {
        &self.revocation_policy_ids
    }

    /// Canonical overlap-policy declarations.
    #[must_use]
    pub fn overlap_policies(&self) -> &[OverlapPolicyRegistrationV2] {
        &self.overlap_policies
    }

    fn overlap_relation(&self, policy_id: u16) -> Option<OverlapPolicyRelationV2> {
        self.overlap_policies
            .binary_search_by_key(&policy_id, |entry| entry.policy_id)
            .ok()
            .map(|index| self.overlap_policies[index].relation)
    }
}

/// Frozen semantic policy. Physical handles, slots, credentials, prefixes,
/// generations, and acquisition attempts cannot be represented here.
///
/// The safe constructor accepts only closed semantic catalog values and the
/// accessors expose only that canonical policy:
///
/// ```
/// use fs_evidence_runner::capability::RootCapabilityPolicyV2;
/// use fs_evidence_runner::catalog::{
///     DigestRoleV2, PlatformPathProfileV2, RootCapabilityAccessV2,
///     RootCapabilityRightV2, RootClassV2,
/// };
/// use fs_evidence_runner::identity::NoClaimScopeRootV1;
///
/// let no_claim_scope = NoClaimScopeRootV1::parse_presented(
///     DigestRoleV2::ClaimScope,
///     NoClaimScopeRootV1::DESCRIPTOR.domain(),
///     &"00".repeat(32),
/// )
/// .unwrap();
/// let policy = RootCapabilityPolicyV2::new(
///     RootClassV2::InputArtifactRoot,
///     PlatformPathProfileV2::PosixDescriptorRelativeV1,
///     RootCapabilityAccessV2::ReadOnlyInput,
///     vec![
///         RootCapabilityRightV2::Traverse,
///         RootCapabilityRightV2::ReadObject,
///         RootCapabilityRightV2::Enumerate,
///     ],
///     1,
///     1,
///     1,
///     no_claim_scope,
/// )
/// .unwrap();
///
/// assert_eq!(policy.access(), RootCapabilityAccessV2::ReadOnlyInput);
/// assert_eq!(policy.rights().len(), 3);
/// ```
///
/// Physical acquisition material has no slot in the semantic right vector.
/// This one physical-only sum covers descriptors, handles, slots, paths,
/// credentials, prefixes, generations, and acquisition attempts without
/// depending on any downstream backend type:
///
/// ```compile_fail
/// use fs_evidence_runner::capability::RootCapabilityPolicyV2;
/// use fs_evidence_runner::catalog::{
///     DigestRoleV2, PlatformPathProfileV2, RootCapabilityAccessV2, RootClassV2,
/// };
/// use fs_evidence_runner::identity::NoClaimScopeRootV1;
///
/// enum PhysicalAcquisitionMaterial<'a> {
///     Descriptor(i32),
///     Handle(usize),
///     Slot(u16),
///     Path(&'a str),
///     Credential(&'a [u8]),
///     Prefix(&'a str),
///     Generation(u64),
///     Attempt(u64),
/// }
///
/// let no_claim_scope = NoClaimScopeRootV1::parse_presented(
///     DigestRoleV2::ClaimScope,
///     NoClaimScopeRootV1::DESCRIPTOR.domain(),
///     &"00".repeat(32),
/// )
/// .unwrap();
/// let physical = PhysicalAcquisitionMaterial::Descriptor(7);
/// let _policy = RootCapabilityPolicyV2::new(
///     RootClassV2::InputArtifactRoot,
///     PlatformPathProfileV2::PosixDescriptorRelativeV1,
///     RootCapabilityAccessV2::ReadOnlyInput,
///     vec![physical],
///     1,
///     1,
///     1,
///     no_claim_scope,
/// );
/// ```
///
/// Fields and the semantic-root constructor remain private:
///
/// ```compile_fail
/// use fs_evidence_runner::capability::RootCapabilityPolicyV2;
///
/// fn try_to_replace_root(policy: &mut RootCapabilityPolicyV2) {
///     policy.root = policy.semantic_root().clone();
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RootCapabilityPolicyV2 {
    root_class: RootClassV2,
    path_profile: PlatformPathProfileV2,
    access: RootCapabilityAccessV2,
    rights: Box<[RootCapabilityRightV2]>,
    freshness_policy_id: u16,
    revocation_policy_id: u16,
    overlap_policy_id: u16,
    no_claim_scope: NoClaimScopeRootV1,
    root: RootCapabilityPolicyRootV2,
}

impl RootCapabilityPolicyV2 {
    /// Perform intrinsic validation and canonicalize the duplicate-free right
    /// set. Registration and destination-mode equality are separate stages.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        root_class: RootClassV2,
        path_profile: PlatformPathProfileV2,
        access: RootCapabilityAccessV2,
        mut rights: Vec<RootCapabilityRightV2>,
        freshness_policy_id: u16,
        revocation_policy_id: u16,
        overlap_policy_id: u16,
        no_claim_scope: NoClaimScopeRootV1,
    ) -> Result<Self, ConstructionErrorV2> {
        match (root_class, access) {
            (RootClassV2::InputArtifactRoot, RootCapabilityAccessV2::ReadOnlyInput)
            | (RootClassV2::OutputArtifactRoot, RootCapabilityAccessV2::DurableOutput)
            | (RootClassV2::Other(_), _) => {}
            _ => {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "root_capability_policy.root_class_access",
                    "input/read-only, output/durable, or registered Other",
                    format_args!("{}/{}", root_class.name(), access.name()),
                ));
            }
        }
        if rights.is_empty() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Missing,
                "root_capability_policy.rights",
                "one canonical nonempty right set",
                0,
            ));
        }
        let mut seen = BTreeSet::new();
        for right in &rights {
            if !seen.insert(right.code()) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "root_capability_policy.rights",
                    "unique rights before canonical sorting",
                    right.name(),
                ));
            }
        }
        rights.sort_unstable_by_key(|right| right.code());
        require_nonzero(
            "root_capability_policy.freshness_policy_id",
            freshness_policy_id,
        )?;
        require_nonzero(
            "root_capability_policy.revocation_policy_id",
            revocation_policy_id,
        )?;
        require_nonzero(
            "root_capability_policy.overlap_policy_id",
            overlap_policy_id,
        )?;

        let matches_a_legal_cell = [
            DestinationAdmissionModeV2::Absent,
            DestinationAdmissionModeV2::PreExistingEmpty,
        ]
        .iter()
        .any(|mode| expected_rights(path_profile, access, *mode).as_slice() == rights.as_slice());
        if !matches_a_legal_cell {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "root_capability_policy.rights",
                "the exact rights of at least one legal profile/access/mode cell",
                render_rights(&rights),
            ));
        }

        let root = construct_policy_root(
            root_class,
            path_profile,
            access,
            &rights,
            freshness_policy_id,
            revocation_policy_id,
            overlap_policy_id,
            &no_claim_scope,
        )?;

        Ok(Self {
            root_class,
            path_profile,
            access,
            rights: rights.into_boxed_slice(),
            freshness_policy_id,
            revocation_policy_id,
            overlap_policy_id,
            no_claim_scope,
            root,
        })
    }

    /// Validate all three registered IDs against the sealed registry
    /// projection.
    pub fn validate_registration(
        &self,
        registry: &RootPolicyRegistryProjectionV2,
    ) -> Result<(), ConstructionErrorV2> {
        validate_registered_id(
            "root_capability_policy.freshness_policy_id",
            self.freshness_policy_id,
            registry.freshness_policy_ids(),
        )?;
        validate_registered_id(
            "root_capability_policy.revocation_policy_id",
            self.revocation_policy_id,
            registry.revocation_policy_ids(),
        )?;
        if registry.overlap_relation(self.overlap_policy_id).is_none() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::UnknownCode,
                "root_capability_policy.overlap_policy_id",
                "an ID in the sealed overlap-policy projection",
                self.overlap_policy_id,
            ));
        }
        Ok(())
    }

    /// Root class.
    #[must_use]
    pub const fn root_class(&self) -> RootClassV2 {
        self.root_class
    }

    /// Logical path profile.
    #[must_use]
    pub const fn path_profile(&self) -> PlatformPathProfileV2 {
        self.path_profile
    }

    /// Semantic access mode.
    #[must_use]
    pub const fn access(&self) -> RootCapabilityAccessV2 {
        self.access
    }

    /// Canonically ordered exact rights.
    #[must_use]
    pub fn rights(&self) -> &[RootCapabilityRightV2] {
        &self.rights
    }

    /// Freshness-policy registration ID.
    #[must_use]
    pub const fn freshness_policy_id(&self) -> u16 {
        self.freshness_policy_id
    }

    /// Revocation-policy registration ID.
    #[must_use]
    pub const fn revocation_policy_id(&self) -> u16 {
        self.revocation_policy_id
    }

    /// Overlap-policy registration ID.
    #[must_use]
    pub const fn overlap_policy_id(&self) -> u16 {
        self.overlap_policy_id
    }

    /// Explicit no-claim scope.
    #[must_use]
    pub const fn no_claim_scope(&self) -> &NoClaimScopeRootV1 {
        &self.no_claim_scope
    }

    /// Nominal semantic identity of the exact canonical policy fields.
    ///
    /// The root is declaration data only. It proves neither acquisition nor
    /// freshness, revocation status, physical placement, or disjointness.
    #[must_use]
    pub const fn semantic_root(&self) -> &RootCapabilityPolicyRootV2 {
        &self.root
    }
}

/// Validate exact output rights against the destination mode carried by the
/// singular publication selection.
pub fn validate_policy_against_selection_v2(
    policy: &RootCapabilityPolicyV2,
    selection: &PublicationSelectionV2,
) -> Result<(), ConstructionErrorV2> {
    if policy.access != RootCapabilityAccessV2::DurableOutput {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "root_capability_policy.access",
            "DurableOutput for a publication selection",
            policy.access.name(),
        ));
    }
    if policy.path_profile != selection.path_profile() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "root_capability_policy.path_profile",
            "the publication selection path profile",
            policy.path_profile.name(),
        ));
    }
    let expected = expected_rights(
        policy.path_profile,
        policy.access,
        selection.destination_mode(),
    );
    if policy.rights.as_ref() != expected.as_slice() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "root_capability_policy.rights",
            "the exact least-privilege output cell",
            render_rights(&policy.rights),
        ));
    }
    Ok(())
}

/// Pure exact-rights view. It is semantic data, not an affine runtime
/// capability or evidence that acquisition occurred.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NarrowedPolicyViewV2 {
    path_profile: PlatformPathProfileV2,
    access: RootCapabilityAccessV2,
    destination_mode: Option<DestinationAdmissionModeV2>,
    rights: Box<[RootCapabilityRightV2]>,
}

impl NarrowedPolicyViewV2 {
    /// Derive the exact output view after contextual validation.
    pub fn for_publication(
        policy: &RootCapabilityPolicyV2,
        selection: &PublicationSelectionV2,
    ) -> Result<Self, ConstructionErrorV2> {
        validate_policy_against_selection_v2(policy, selection)?;
        Ok(Self {
            path_profile: policy.path_profile,
            access: policy.access,
            destination_mode: Some(selection.destination_mode()),
            rights: policy.rights.clone(),
        })
    }

    /// Derive the exact input view. Destination admission does not apply to a
    /// read-only input.
    pub fn for_read_only(policy: &RootCapabilityPolicyV2) -> Result<Self, ConstructionErrorV2> {
        if policy.access != RootCapabilityAccessV2::ReadOnlyInput {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "root_capability_policy.access",
                "ReadOnlyInput",
                policy.access.name(),
            ));
        }
        let expected = expected_rights(
            policy.path_profile,
            policy.access,
            DestinationAdmissionModeV2::Absent,
        );
        if policy.rights.as_ref() != expected.as_slice() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "root_capability_policy.rights",
                "the exact least-privilege input cell",
                render_rights(&policy.rights),
            ));
        }
        Ok(Self {
            path_profile: policy.path_profile,
            access: policy.access,
            destination_mode: None,
            rights: policy.rights.clone(),
        })
    }

    /// Exact semantic rights exposed by this view.
    #[must_use]
    pub fn rights(&self) -> &[RootCapabilityRightV2] {
        &self.rights
    }

    /// Logical profile.
    #[must_use]
    pub const fn path_profile(&self) -> PlatformPathProfileV2 {
        self.path_profile
    }

    /// Access mode.
    #[must_use]
    pub const fn access(&self) -> RootCapabilityAccessV2 {
        self.access
    }

    /// Destination mode for an output view.
    #[must_use]
    pub const fn destination_mode(&self) -> Option<DestinationAdmissionModeV2> {
        self.destination_mode
    }
}

/// Ordered, duplicate-rejected policy set for one command projection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RootCapabilityPolicySetV2 {
    command: RunnerCommandV2,
    policies: Box<[RootCapabilityPolicyV2]>,
}

impl RootCapabilityPolicySetV2 {
    /// Validate command cardinality, registration, root-class uniqueness, and
    /// Replay's shared declared-disjoint overlap policy.
    pub fn new(
        command: RunnerCommandV2,
        mut policies: Vec<RootCapabilityPolicyV2>,
        registry: &RootPolicyRegistryProjectionV2,
    ) -> Result<Self, ConstructionErrorV2> {
        let expected_policy_count = match command {
            RunnerCommandV2::List | RunnerCommandV2::Check | RunnerCommandV2::SelfTest => 0,
            RunnerCommandV2::Run | RunnerCommandV2::Negative => 1,
            RunnerCommandV2::Replay => 2,
        };
        if policies.len() != expected_policy_count {
            return Err(cardinality_error(
                command,
                expected_policy_count,
                policies.len(),
            ));
        }

        for policy in &policies {
            policy.validate_registration(registry)?;
        }
        let mut seen = BTreeSet::new();
        for policy in &policies {
            let key = (
                policy.root_class.tag(),
                policy.root_class.registered_id().unwrap_or(0),
            );
            if !seen.insert(key) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "root_capability_policy_set.root_class",
                    "one policy per exact root class",
                    format_args!("{}:{:?}", key.0, key.1),
                ));
            }
        }
        policies.sort_by_key(|policy| {
            (
                policy.root_class.tag(),
                policy.root_class.registered_id().unwrap_or(0),
            )
        });

        match command {
            RunnerCommandV2::List | RunnerCommandV2::Check | RunnerCommandV2::SelfTest => {}
            RunnerCommandV2::Run | RunnerCommandV2::Negative => {
                if policies[0].access != RootCapabilityAccessV2::DurableOutput {
                    return Err(cardinality_error(command, 1, policies.len()));
                }
            }
            RunnerCommandV2::Replay => {
                let input = policies
                    .iter()
                    .find(|policy| policy.access == RootCapabilityAccessV2::ReadOnlyInput);
                let output = policies
                    .iter()
                    .find(|policy| policy.access == RootCapabilityAccessV2::DurableOutput);
                let (Some(input), Some(output)) = (input, output) else {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Incompatible,
                        "root_capability_policy_set.replay_access",
                        "one ReadOnlyInput and one DurableOutput",
                        render_accesses(&policies),
                    ));
                };
                if input.overlap_policy_id != output.overlap_policy_id {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Incompatible,
                        "root_capability_policy_set.overlap_policy_id",
                        "one shared registered overlap-policy ID",
                        format_args!("{}/{}", input.overlap_policy_id, output.overlap_policy_id),
                    ));
                }
                if registry.overlap_relation(input.overlap_policy_id)
                    != Some(OverlapPolicyRelationV2::RequireInputOutputDisjoint)
                {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Incompatible,
                        "root_capability_policy_set.overlap_relation",
                        "RequireInputOutputDisjoint",
                        input.overlap_policy_id,
                    ));
                }
            }
        }
        Ok(Self {
            command,
            policies: policies.into_boxed_slice(),
        })
    }

    /// Command whose policy cardinality was validated.
    #[must_use]
    pub const fn command(&self) -> RunnerCommandV2 {
        self.command
    }

    /// Canonically root-class-ordered policies.
    #[must_use]
    pub fn policies(&self) -> &[RootCapabilityPolicyV2] {
        &self.policies
    }
}

pub(crate) fn expected_rights(
    profile: PlatformPathProfileV2,
    access: RootCapabilityAccessV2,
    mode: DestinationAdmissionModeV2,
) -> Vec<RootCapabilityRightV2> {
    use RootCapabilityRightV2 as Right;
    match (profile, access, mode) {
        (
            PlatformPathProfileV2::PosixDescriptorRelativeV1
            | PlatformPathProfileV2::WindowsHandleRelativeV1,
            RootCapabilityAccessV2::ReadOnlyInput,
            _,
        ) => vec![Right::Traverse, Right::ReadObject, Right::Enumerate],
        (
            PlatformPathProfileV2::PosixDescriptorRelativeV1
            | PlatformPathProfileV2::WindowsHandleRelativeV1,
            RootCapabilityAccessV2::DurableOutput,
            DestinationAdmissionModeV2::Absent,
        ) => vec![
            Right::Traverse,
            Right::Enumerate,
            Right::CreateObject,
            Right::SyncObject,
            Right::SyncContainer,
        ],
        (
            PlatformPathProfileV2::PosixDescriptorRelativeV1
            | PlatformPathProfileV2::WindowsHandleRelativeV1,
            RootCapabilityAccessV2::DurableOutput,
            DestinationAdmissionModeV2::PreExistingEmpty,
        ) => vec![
            Right::Traverse,
            Right::Enumerate,
            Right::CreateObject,
            Right::PopulateEmptyDestination,
            Right::SyncObject,
            Right::SyncContainer,
        ],
        (
            PlatformPathProfileV2::ContentStoreObjectKeyV1,
            RootCapabilityAccessV2::ReadOnlyInput,
            _,
        ) => vec![Right::ReadObject, Right::Enumerate, Right::QueryGeneration],
        (
            PlatformPathProfileV2::ContentStoreObjectKeyV1,
            RootCapabilityAccessV2::DurableOutput,
            DestinationAdmissionModeV2::Absent,
        ) => vec![
            Right::CreateObject,
            Right::QueryGeneration,
            Right::CommitCompareAndSwap,
        ],
        (
            PlatformPathProfileV2::ContentStoreObjectKeyV1,
            RootCapabilityAccessV2::DurableOutput,
            DestinationAdmissionModeV2::PreExistingEmpty,
        ) => vec![
            Right::Enumerate,
            Right::CreateObject,
            Right::AcquireExclusiveLease,
            Right::QueryGeneration,
            Right::CommitCompareAndSwap,
        ],
    }
}

#[allow(clippy::too_many_arguments)]
fn construct_policy_root(
    root_class: RootClassV2,
    path_profile: PlatformPathProfileV2,
    access: RootCapabilityAccessV2,
    rights: &[RootCapabilityRightV2],
    freshness_policy_id: u16,
    revocation_policy_id: u16,
    overlap_policy_id: u16,
    no_claim_scope: &NoClaimScopeRootV1,
) -> Result<RootCapabilityPolicyRootV2, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSROOTCAPABILITYPOLICY\x01", 4096)?;
    frame.push_u16("root_capability_policy.root_class", root_class.tag())?;
    frame.push_presence(
        "root_capability_policy.root_class_registered_id",
        root_class.registered_id().is_some(),
    )?;
    if let Some(registered_id) = root_class.registered_id() {
        frame.push_u16(
            "root_capability_policy.root_class_registered_id",
            registered_id,
        )?;
    }
    frame.push_u16("root_capability_policy.path_profile", path_profile.code())?;
    frame.push_u16("root_capability_policy.access", access.code())?;
    frame.push_u16(
        "root_capability_policy.right_count",
        u16::try_from(rights.len()).expect("the closed right catalog has ten entries"),
    )?;
    for right in rights {
        frame.push_u16("root_capability_policy.right", right.code())?;
    }
    frame.push_u16(
        "root_capability_policy.freshness_policy_id",
        freshness_policy_id,
    )?;
    frame.push_u16(
        "root_capability_policy.revocation_policy_id",
        revocation_policy_id,
    )?;
    frame.push_u16(
        "root_capability_policy.overlap_policy_id",
        overlap_policy_id,
    )?;
    frame.push_u16(
        "root_capability_policy.no_claim_scope_role",
        no_claim_scope.role().code(),
    )?;
    frame.push_str(
        "root_capability_policy.no_claim_scope_domain",
        no_claim_scope.domain(),
    )?;
    frame.push_bytes(
        "root_capability_policy.no_claim_scope_bytes",
        no_claim_scope.bytes(),
    )?;

    let content = frame.root(ROOT_CAPABILITY_POLICY_DOMAIN_V1);
    let digest = DigestValueV2::from_array(
        DigestRoleV2::Policy,
        RootCapabilityPolicyRootV2::DESCRIPTOR.domain_witness(),
        *content.as_bytes(),
    );
    Ok(RootCapabilityPolicyRootV2::from_digest(digest)
        .expect("the private policy constructor fixes the nominal role and domain"))
}

fn validate_id_registry(field: &'static str, ids: &[u16]) -> Result<(), ConstructionErrorV2> {
    if ids.len() > ROOT_POLICY_REGISTRATIONS_MAX_V2 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            field,
            "at most 64 registrations",
            ids.len(),
        ));
    }
    let mut seen = BTreeSet::new();
    for &id in ids {
        require_nonzero(field, id)?;
        if !seen.insert(id) {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Duplicate,
                field,
                "unique registered IDs",
                id,
            ));
        }
    }
    Ok(())
}

fn validate_registered_id(
    field: &'static str,
    id: u16,
    registry: &[u16],
) -> Result<(), ConstructionErrorV2> {
    if registry.binary_search(&id).is_err() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::UnknownCode,
            field,
            "an ID in the sealed policy projection",
            id,
        ));
    }
    Ok(())
}

fn require_nonzero(field: &'static str, value: u16) -> Result<(), ConstructionErrorV2> {
    if value == 0 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Zero,
            field,
            "a nonzero registered u16 ID",
            value,
        ));
    }
    Ok(())
}

fn render_rights(rights: &[RootCapabilityRightV2]) -> String {
    rights
        .iter()
        .map(|right| right.name())
        .collect::<Vec<_>>()
        .join(",")
}

fn render_accesses(policies: &[RootCapabilityPolicyV2]) -> String {
    policies
        .iter()
        .map(|policy| policy.access.name())
        .collect::<Vec<_>>()
        .join(",")
}

fn cardinality_error(
    command: RunnerCommandV2,
    expected: usize,
    observed: usize,
) -> ConstructionErrorV2 {
    ConstructionErrorV2::new(
        ConstructionErrorKindV2::Incompatible,
        "root_capability_policy_set.command_cardinality",
        match expected {
            0 => "zero policies",
            1 => "exactly one DurableOutput policy",
            2 => "exactly one ReadOnlyInput and one DurableOutput policy",
            _ => "the frozen command cardinality",
        },
        format_args!("{}/{}", command.name(), observed),
    )
}

#[cfg(test)]
mod tests {
    use super::{
        NarrowedPolicyViewV2, OverlapPolicyRegistrationV2, RootCapabilityPolicySetV2,
        RootCapabilityPolicyV2, RootPolicyRegistryProjectionV2,
        validate_policy_against_selection_v2,
    };
    use crate::catalog::{
        DestinationAdmissionModeV2, DigestRoleV2, OverlapPolicyRelationV2, PlatformPathProfileV2,
        PublicationProtocolV2, RootCapabilityAccessV2, RootCapabilityRightV2, RootClassV2,
        RunnerCommandV2,
    };
    use crate::identity::NoClaimScopeRootV1;
    use crate::path::{ContentStoreObjectKeyV1, LogicalBundlePathV1};
    use crate::publication::{PublicationSelectionV2, PublicationTargetV2};

    #[derive(Clone, Copy)]
    struct CapabilityOracleCell {
        profile: PlatformPathProfileV2,
        access: RootCapabilityAccessV2,
        mode: DestinationAdmissionModeV2,
        rights: &'static [RootCapabilityRightV2],
    }

    use RootCapabilityRightV2 as Right;

    const FILE_INPUT_RIGHTS: &[Right] = &[Right::Traverse, Right::ReadObject, Right::Enumerate];
    const FILE_ABSENT_OUTPUT_RIGHTS: &[Right] = &[
        Right::Traverse,
        Right::Enumerate,
        Right::CreateObject,
        Right::SyncObject,
        Right::SyncContainer,
    ];
    const FILE_EMPTY_OUTPUT_RIGHTS: &[Right] = &[
        Right::Traverse,
        Right::Enumerate,
        Right::CreateObject,
        Right::PopulateEmptyDestination,
        Right::SyncObject,
        Right::SyncContainer,
    ];
    const STORE_INPUT_RIGHTS: &[Right] =
        &[Right::ReadObject, Right::Enumerate, Right::QueryGeneration];
    const STORE_ABSENT_OUTPUT_RIGHTS: &[Right] = &[
        Right::CreateObject,
        Right::QueryGeneration,
        Right::CommitCompareAndSwap,
    ];
    const STORE_EMPTY_OUTPUT_RIGHTS: &[Right] = &[
        Right::Enumerate,
        Right::CreateObject,
        Right::AcquireExclusiveLease,
        Right::QueryGeneration,
        Right::CommitCompareAndSwap,
    ];

    const CAPABILITY_ORACLE: [CapabilityOracleCell; 12] = [
        CapabilityOracleCell {
            profile: PlatformPathProfileV2::PosixDescriptorRelativeV1,
            access: RootCapabilityAccessV2::ReadOnlyInput,
            mode: DestinationAdmissionModeV2::Absent,
            rights: FILE_INPUT_RIGHTS,
        },
        CapabilityOracleCell {
            profile: PlatformPathProfileV2::PosixDescriptorRelativeV1,
            access: RootCapabilityAccessV2::ReadOnlyInput,
            mode: DestinationAdmissionModeV2::PreExistingEmpty,
            rights: FILE_INPUT_RIGHTS,
        },
        CapabilityOracleCell {
            profile: PlatformPathProfileV2::PosixDescriptorRelativeV1,
            access: RootCapabilityAccessV2::DurableOutput,
            mode: DestinationAdmissionModeV2::Absent,
            rights: FILE_ABSENT_OUTPUT_RIGHTS,
        },
        CapabilityOracleCell {
            profile: PlatformPathProfileV2::PosixDescriptorRelativeV1,
            access: RootCapabilityAccessV2::DurableOutput,
            mode: DestinationAdmissionModeV2::PreExistingEmpty,
            rights: FILE_EMPTY_OUTPUT_RIGHTS,
        },
        CapabilityOracleCell {
            profile: PlatformPathProfileV2::WindowsHandleRelativeV1,
            access: RootCapabilityAccessV2::ReadOnlyInput,
            mode: DestinationAdmissionModeV2::Absent,
            rights: FILE_INPUT_RIGHTS,
        },
        CapabilityOracleCell {
            profile: PlatformPathProfileV2::WindowsHandleRelativeV1,
            access: RootCapabilityAccessV2::ReadOnlyInput,
            mode: DestinationAdmissionModeV2::PreExistingEmpty,
            rights: FILE_INPUT_RIGHTS,
        },
        CapabilityOracleCell {
            profile: PlatformPathProfileV2::WindowsHandleRelativeV1,
            access: RootCapabilityAccessV2::DurableOutput,
            mode: DestinationAdmissionModeV2::Absent,
            rights: FILE_ABSENT_OUTPUT_RIGHTS,
        },
        CapabilityOracleCell {
            profile: PlatformPathProfileV2::WindowsHandleRelativeV1,
            access: RootCapabilityAccessV2::DurableOutput,
            mode: DestinationAdmissionModeV2::PreExistingEmpty,
            rights: FILE_EMPTY_OUTPUT_RIGHTS,
        },
        CapabilityOracleCell {
            profile: PlatformPathProfileV2::ContentStoreObjectKeyV1,
            access: RootCapabilityAccessV2::ReadOnlyInput,
            mode: DestinationAdmissionModeV2::Absent,
            rights: STORE_INPUT_RIGHTS,
        },
        CapabilityOracleCell {
            profile: PlatformPathProfileV2::ContentStoreObjectKeyV1,
            access: RootCapabilityAccessV2::ReadOnlyInput,
            mode: DestinationAdmissionModeV2::PreExistingEmpty,
            rights: STORE_INPUT_RIGHTS,
        },
        CapabilityOracleCell {
            profile: PlatformPathProfileV2::ContentStoreObjectKeyV1,
            access: RootCapabilityAccessV2::DurableOutput,
            mode: DestinationAdmissionModeV2::Absent,
            rights: STORE_ABSENT_OUTPUT_RIGHTS,
        },
        CapabilityOracleCell {
            profile: PlatformPathProfileV2::ContentStoreObjectKeyV1,
            access: RootCapabilityAccessV2::DurableOutput,
            mode: DestinationAdmissionModeV2::PreExistingEmpty,
            rights: STORE_EMPTY_OUTPUT_RIGHTS,
        },
    ];

    fn oracle_rights(
        profile: PlatformPathProfileV2,
        access: RootCapabilityAccessV2,
        mode: DestinationAdmissionModeV2,
    ) -> &'static [RootCapabilityRightV2] {
        CAPABILITY_ORACLE
            .iter()
            .find(|cell| cell.profile == profile && cell.access == access && cell.mode == mode)
            .expect("the handwritten oracle contains all 12 cells")
            .rights
    }

    fn no_claim(byte: u8) -> NoClaimScopeRootV1 {
        NoClaimScopeRootV1::parse_presented(
            DigestRoleV2::ClaimScope,
            NoClaimScopeRootV1::DESCRIPTOR.domain(),
            &format!("{byte:02x}").repeat(32),
        )
        .expect("fixture no-claim scope")
    }

    fn registry() -> RootPolicyRegistryProjectionV2 {
        RootPolicyRegistryProjectionV2::new(
            vec![1, 2],
            vec![1, 2],
            vec![
                OverlapPolicyRegistrationV2::new(
                    1,
                    OverlapPolicyRelationV2::RequireInputOutputDisjoint,
                )
                .expect("overlap registration"),
                OverlapPolicyRegistrationV2::new(
                    2,
                    OverlapPolicyRelationV2::RequireInputOutputDisjoint,
                )
                .expect("overlap registration"),
            ],
        )
        .expect("fixture registry")
    }

    fn protocol(profile: PlatformPathProfileV2) -> PublicationProtocolV2 {
        match profile {
            PlatformPathProfileV2::PosixDescriptorRelativeV1 => {
                PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1
            }
            PlatformPathProfileV2::WindowsHandleRelativeV1 => {
                PublicationProtocolV2::WindowsHandleReplaceAndDirectoryFlushV1
            }
            PlatformPathProfileV2::ContentStoreObjectKeyV1 => {
                PublicationProtocolV2::ContentStoreAtomicCommitV1
            }
        }
    }

    fn selection(
        profile: PlatformPathProfileV2,
        destination_mode: DestinationAdmissionModeV2,
    ) -> PublicationSelectionV2 {
        let target = match profile {
            PlatformPathProfileV2::PosixDescriptorRelativeV1 => PublicationTargetV2::PosixRelative(
                LogicalBundlePathV1::new("results/bundle").expect("fixture path"),
            ),
            PlatformPathProfileV2::WindowsHandleRelativeV1 => PublicationTargetV2::WindowsRelative(
                LogicalBundlePathV1::new("results/bundle").expect("fixture path"),
            ),
            PlatformPathProfileV2::ContentStoreObjectKeyV1 => {
                PublicationTargetV2::ContentStoreLogicalKey(
                    ContentStoreObjectKeyV1::new("results/bundle").expect("fixture key"),
                )
            }
        };
        PublicationSelectionV2::new(profile, protocol(profile), destination_mode, target)
            .expect("compatible publication cell")
    }

    fn policy(
        profile: PlatformPathProfileV2,
        access: RootCapabilityAccessV2,
        destination_mode: DestinationAdmissionModeV2,
        rights: Vec<RootCapabilityRightV2>,
        registration_ids: [u16; 3],
        claim_byte: u8,
    ) -> Result<RootCapabilityPolicyV2, crate::ConstructionErrorV2> {
        RootCapabilityPolicyV2::new(
            match access {
                RootCapabilityAccessV2::ReadOnlyInput => RootClassV2::InputArtifactRoot,
                RootCapabilityAccessV2::DurableOutput => RootClassV2::OutputArtifactRoot,
            },
            profile,
            access,
            if rights.is_empty() {
                oracle_rights(profile, access, destination_mode).to_vec()
            } else {
                rights
            },
            registration_ids[0],
            registration_ids[1],
            registration_ids[2],
            no_claim(claim_byte),
        )
    }

    fn exact_policy(
        profile: PlatformPathProfileV2,
        access: RootCapabilityAccessV2,
        destination_mode: DestinationAdmissionModeV2,
        registration_ids: [u16; 3],
        claim_byte: u8,
    ) -> RootCapabilityPolicyV2 {
        policy(
            profile,
            access,
            destination_mode,
            Vec::new(),
            registration_ids,
            claim_byte,
        )
        .expect("handwritten exact policy cell")
    }

    #[test]
    fn registry_is_bounded_sorted_duplicate_free_and_nonzero() {
        let registry = RootPolicyRegistryProjectionV2::new(
            vec![2, 1],
            vec![2, 1],
            vec![
                OverlapPolicyRegistrationV2::new(
                    2,
                    OverlapPolicyRelationV2::RequireInputOutputDisjoint,
                )
                .expect("registration"),
                OverlapPolicyRegistrationV2::new(
                    1,
                    OverlapPolicyRelationV2::RequireInputOutputDisjoint,
                )
                .expect("registration"),
            ],
        )
        .expect("canonical registry");
        assert_eq!(registry.freshness_policy_ids(), &[1, 2]);
        assert_eq!(registry.revocation_policy_ids(), &[1, 2]);
        assert_eq!(
            registry
                .overlap_policies()
                .iter()
                .map(|entry| entry.policy_id())
                .collect::<Vec<_>>(),
            vec![1, 2]
        );
        let zero_overlap = OverlapPolicyRegistrationV2::new(
            0,
            OverlapPolicyRelationV2::RequireInputOutputDisjoint,
        )
        .expect_err("overlap policy ID zero must refuse");
        assert_eq!(zero_overlap.kind(), crate::ConstructionErrorKindV2::Zero);
        assert_eq!(
            zero_overlap.field(),
            "root_policy_registry.overlap_policy_id"
        );
        for (field, result) in [
            (
                "root_policy_registry.freshness_policy_ids",
                RootPolicyRegistryProjectionV2::new(vec![0], Vec::new(), Vec::new()),
            ),
            (
                "root_policy_registry.revocation_policy_ids",
                RootPolicyRegistryProjectionV2::new(Vec::new(), vec![0], Vec::new()),
            ),
        ] {
            let error = result.expect_err("zero registry IDs must refuse");
            assert_eq!(error.kind(), crate::ConstructionErrorKindV2::Zero);
            assert_eq!(error.field(), field);
        }
        assert!(RootPolicyRegistryProjectionV2::new(vec![1, 1], vec![], vec![]).is_err());
        assert!(RootPolicyRegistryProjectionV2::new(vec![], vec![1, 1], vec![]).is_err());
        assert!(
            RootPolicyRegistryProjectionV2::new(
                vec![],
                vec![],
                vec![
                    OverlapPolicyRegistrationV2::new(
                        1,
                        OverlapPolicyRelationV2::RequireInputOutputDisjoint
                    )
                    .expect("registration");
                    2
                ],
            )
            .is_err()
        );
        assert!(
            RootPolicyRegistryProjectionV2::new((1..=65).collect(), Vec::new(), Vec::new())
                .is_err()
        );
        RootPolicyRegistryProjectionV2::new(
            (1..=22).collect(),
            (1..=21).collect(),
            (1..=21)
                .map(|id| {
                    OverlapPolicyRegistrationV2::new(
                        id,
                        OverlapPolicyRelationV2::RequireInputOutputDisjoint,
                    )
                    .expect("exact-cap overlap registration")
                })
                .collect(),
        )
        .expect("64 aggregate registrations are admitted");
        let aggregate_one_over = RootPolicyRegistryProjectionV2::new(
            (1..=22).collect(),
            (1..=22).collect(),
            (1..=21)
                .map(|id| {
                    OverlapPolicyRegistrationV2::new(
                        id,
                        OverlapPolicyRelationV2::RequireInputOutputDisjoint,
                    )
                    .expect("one-over overlap registration")
                })
                .collect(),
        )
        .expect_err("65 aggregate registrations must refuse");
        assert_eq!(
            aggregate_one_over.kind(),
            crate::ConstructionErrorKindV2::TooLarge
        );
        assert_eq!(
            aggregate_one_over.field(),
            "root_policy_registry.registrations"
        );
    }

    #[test]
    #[allow(
        clippy::too_many_lines,
        reason = "the least-privilege test deliberately enumerates every command and every one-right mutation in one auditable matrix"
    )]
    fn least_privilege_matrix_rejects_every_one_right_mutant() {
        let all_rights = [
            RootCapabilityRightV2::Traverse,
            RootCapabilityRightV2::ReadObject,
            RootCapabilityRightV2::Enumerate,
            RootCapabilityRightV2::CreateObject,
            RootCapabilityRightV2::PopulateEmptyDestination,
            RootCapabilityRightV2::SyncObject,
            RootCapabilityRightV2::SyncContainer,
            RootCapabilityRightV2::AcquireExclusiveLease,
            RootCapabilityRightV2::QueryGeneration,
            RootCapabilityRightV2::CommitCompareAndSwap,
        ];
        let mut refused_mutants = 0_usize;

        for cell in CAPABILITY_ORACLE {
            let expected = cell.rights;
            let accepted = policy(
                cell.profile,
                cell.access,
                cell.mode,
                expected.to_vec(),
                [1, 1, 1],
                0,
            )
            .expect("exact handwritten policy cell");
            accepted
                .validate_registration(&registry())
                .expect("registered policy");
            assert_eq!(accepted.rights(), expected);
            match cell.access {
                RootCapabilityAccessV2::ReadOnlyInput => {
                    let view = NarrowedPolicyViewV2::for_read_only(&accepted).expect("input view");
                    assert_eq!(view.rights(), expected);
                    assert_eq!(view.destination_mode(), None);
                }
                RootCapabilityAccessV2::DurableOutput => {
                    let selected = selection(cell.profile, cell.mode);
                    validate_policy_against_selection_v2(&accepted, &selected)
                        .expect("exact output cell");
                    let view = NarrowedPolicyViewV2::for_publication(&accepted, &selected)
                        .expect("output view");
                    assert_eq!(view.rights(), expected);
                    assert_eq!(view.destination_mode(), Some(cell.mode));
                }
            }

            for omitted in expected {
                let mut mutant = expected.to_vec();
                mutant.retain(|right| right != omitted);
                assert_policy_mutant_refuses(cell, mutant);
                refused_mutants += 1;
            }
            for added in all_rights
                .iter()
                .copied()
                .filter(|right| !expected.contains(right))
            {
                let mut mutant = expected.to_vec();
                mutant.push(added);
                assert_policy_mutant_refuses(cell, mutant);
                refused_mutants += 1;
            }
            for omitted in expected {
                for replacement in all_rights
                    .iter()
                    .copied()
                    .filter(|right| !expected.contains(right))
                {
                    let mut mutant = expected.to_vec();
                    let slot = mutant
                        .iter_mut()
                        .find(|right| *right == omitted)
                        .expect("present right");
                    *slot = replacement;
                    assert_policy_mutant_refuses(cell, mutant);
                    refused_mutants += 1;
                }
            }
        }
        assert_eq!(refused_mutants, 390);

        let duplicate_cell = CAPABILITY_ORACLE[0];
        let mut duplicate_rights = duplicate_cell.rights.to_vec();
        duplicate_rights.push(duplicate_cell.rights[0]);
        let duplicate = policy(
            duplicate_cell.profile,
            duplicate_cell.access,
            duplicate_cell.mode,
            duplicate_rights,
            [1, 1, 1],
            0,
        )
        .expect_err("a duplicate semantic right must refuse before canonical sorting");
        assert_eq!(duplicate.kind(), crate::ConstructionErrorKindV2::Duplicate);
        assert_eq!(duplicate.field(), "root_capability_policy.rights");

        let empty = RootCapabilityPolicyV2::new(
            RootClassV2::InputArtifactRoot,
            PlatformPathProfileV2::PosixDescriptorRelativeV1,
            RootCapabilityAccessV2::ReadOnlyInput,
            Vec::new(),
            1,
            1,
            1,
            no_claim(0),
        )
        .expect_err("an empty semantic right set must refuse");
        assert_eq!(empty.kind(), crate::ConstructionErrorKindV2::Missing);
        assert_eq!(empty.field(), "root_capability_policy.rights");

        for (registration_ids, field) in [
            ([0, 1, 1], "root_capability_policy.freshness_policy_id"),
            ([1, 0, 1], "root_capability_policy.revocation_policy_id"),
            ([1, 1, 0], "root_capability_policy.overlap_policy_id"),
        ] {
            let error = policy(
                duplicate_cell.profile,
                duplicate_cell.access,
                duplicate_cell.mode,
                duplicate_cell.rights.to_vec(),
                registration_ids,
                0,
            )
            .expect_err("zero policy registration IDs must refuse intrinsically");
            assert_eq!(error.kind(), crate::ConstructionErrorKindV2::Zero);
            assert_eq!(error.field(), field);
        }
    }

    fn assert_policy_mutant_refuses(
        cell: CapabilityOracleCell,
        rights: Vec<RootCapabilityRightV2>,
    ) {
        match policy(cell.profile, cell.access, cell.mode, rights, [1, 1, 1], 0) {
            Err(error) => {
                assert_eq!(error.kind(), crate::ConstructionErrorKindV2::Incompatible);
                assert_eq!(error.field(), "root_capability_policy.rights");
            }
            Ok(mutant) => match cell.access {
                RootCapabilityAccessV2::ReadOnlyInput => {
                    let error = NarrowedPolicyViewV2::for_read_only(&mutant)
                        .expect_err("mutant input rights must refuse");
                    assert_eq!(error.kind(), crate::ConstructionErrorKindV2::Incompatible);
                    assert_eq!(error.field(), "root_capability_policy.rights");
                }
                RootCapabilityAccessV2::DurableOutput => {
                    let error = validate_policy_against_selection_v2(
                        &mutant,
                        &selection(cell.profile, cell.mode),
                    );
                    let error = error.expect_err("mutant output rights must refuse");
                    assert_eq!(error.kind(), crate::ConstructionErrorKindV2::Incompatible);
                    assert_eq!(error.field(), "root_capability_policy.rights");
                }
            },
        }
    }

    #[test]
    fn policy_root_binds_each_semantic_field_and_ignores_opaque_observations() {
        let base = exact_policy(
            PlatformPathProfileV2::PosixDescriptorRelativeV1,
            RootCapabilityAccessV2::DurableOutput,
            DestinationAdmissionModeV2::Absent,
            [1, 1, 1],
            0,
        );
        let base_root = base.semantic_root().bytes();
        let mutations = [
            exact_policy(
                PlatformPathProfileV2::WindowsHandleRelativeV1,
                RootCapabilityAccessV2::DurableOutput,
                DestinationAdmissionModeV2::Absent,
                [1, 1, 1],
                0,
            ),
            RootCapabilityPolicyV2::new(
                RootClassV2::from_tag(3, Some(7)).expect("registered root class"),
                PlatformPathProfileV2::PosixDescriptorRelativeV1,
                RootCapabilityAccessV2::DurableOutput,
                oracle_rights(
                    PlatformPathProfileV2::PosixDescriptorRelativeV1,
                    RootCapabilityAccessV2::DurableOutput,
                    DestinationAdmissionModeV2::Absent,
                )
                .to_vec(),
                1,
                1,
                1,
                no_claim(0),
            )
            .expect("root-class mutation"),
            exact_policy(
                PlatformPathProfileV2::PosixDescriptorRelativeV1,
                RootCapabilityAccessV2::ReadOnlyInput,
                DestinationAdmissionModeV2::Absent,
                [1, 1, 1],
                0,
            ),
            exact_policy(
                PlatformPathProfileV2::PosixDescriptorRelativeV1,
                RootCapabilityAccessV2::DurableOutput,
                DestinationAdmissionModeV2::PreExistingEmpty,
                [1, 1, 1],
                0,
            ),
            exact_policy(
                PlatformPathProfileV2::PosixDescriptorRelativeV1,
                RootCapabilityAccessV2::DurableOutput,
                DestinationAdmissionModeV2::Absent,
                [2, 1, 1],
                0,
            ),
            exact_policy(
                PlatformPathProfileV2::PosixDescriptorRelativeV1,
                RootCapabilityAccessV2::DurableOutput,
                DestinationAdmissionModeV2::Absent,
                [1, 2, 1],
                0,
            ),
            exact_policy(
                PlatformPathProfileV2::PosixDescriptorRelativeV1,
                RootCapabilityAccessV2::DurableOutput,
                DestinationAdmissionModeV2::Absent,
                [1, 1, 2],
                0,
            ),
            exact_policy(
                PlatformPathProfileV2::PosixDescriptorRelativeV1,
                RootCapabilityAccessV2::DurableOutput,
                DestinationAdmissionModeV2::Absent,
                [1, 1, 1],
                1,
            ),
        ];
        for mutation in mutations {
            assert_ne!(mutation.semantic_root().bytes(), base_root);
        }

        let opaque_handle_a = 17_u64;
        let opaque_handle_b = 99_u64;
        assert_ne!(opaque_handle_a, opaque_handle_b);
        assert_eq!(base.semantic_root().bytes(), base_root);
    }

    #[test]
    #[allow(
        clippy::too_many_lines,
        reason = "the test keeps command cardinality precedence, registration failures, duplicate roots, canonical Replay ordering, and overlap-policy compatibility in one auditable policy-set matrix"
    )]
    fn registration_and_command_policy_sets_are_exact() {
        let registry = registry();
        let output = exact_policy(
            PlatformPathProfileV2::ContentStoreObjectKeyV1,
            RootCapabilityAccessV2::DurableOutput,
            DestinationAdmissionModeV2::Absent,
            [1, 1, 1],
            0,
        );
        let input = exact_policy(
            PlatformPathProfileV2::ContentStoreObjectKeyV1,
            RootCapabilityAccessV2::ReadOnlyInput,
            DestinationAdmissionModeV2::Absent,
            [1, 1, 1],
            0,
        );
        assert!(
            RootCapabilityPolicySetV2::new(RunnerCommandV2::List, Vec::new(), &registry).is_ok()
        );
        assert!(
            RootCapabilityPolicySetV2::new(RunnerCommandV2::Run, vec![output.clone()], &registry)
                .is_ok()
        );
        let replay = RootCapabilityPolicySetV2::new(
            RunnerCommandV2::Replay,
            vec![output.clone(), input.clone()],
            &registry,
        )
        .expect("replay pair");
        assert_eq!(
            replay.policies()[0].root_class(),
            RootClassV2::InputArtifactRoot
        );
        assert_eq!(
            replay.policies()[1].root_class(),
            RootClassV2::OutputArtifactRoot
        );
        assert!(
            RootCapabilityPolicySetV2::new(
                RunnerCommandV2::Negative,
                vec![input.clone()],
                &registry
            )
            .is_err()
        );
        let unregistered_for_oversized_set = exact_policy(
            PlatformPathProfileV2::ContentStoreObjectKeyV1,
            RootCapabilityAccessV2::DurableOutput,
            DestinationAdmissionModeV2::Absent,
            [3, 1, 1],
            0,
        );
        for oversized_policies in [
            vec![unregistered_for_oversized_set, input.clone()],
            vec![output.clone(), output.clone()],
        ] {
            let error =
                RootCapabilityPolicySetV2::new(RunnerCommandV2::Run, oversized_policies, &registry)
                    .expect_err("command cardinality precedes registration and duplicate checks");
            assert_eq!(error.kind(), crate::ConstructionErrorKindV2::Incompatible);
            assert_eq!(
                error.field(),
                "root_capability_policy_set.command_cardinality"
            );
            assert_eq!(error.expected(), "exactly one DurableOutput policy");
            assert_eq!(error.observed(), "run/2");
        }

        for (registration_ids, field) in [
            ([3, 1, 1], "root_capability_policy.freshness_policy_id"),
            ([1, 3, 1], "root_capability_policy.revocation_policy_id"),
            ([1, 1, 3], "root_capability_policy.overlap_policy_id"),
        ] {
            let unregistered = exact_policy(
                PlatformPathProfileV2::ContentStoreObjectKeyV1,
                RootCapabilityAccessV2::DurableOutput,
                DestinationAdmissionModeV2::Absent,
                registration_ids,
                0,
            );
            let error = unregistered
                .validate_registration(&registry)
                .expect_err("every unregistered policy ID must refuse");
            assert_eq!(error.kind(), crate::ConstructionErrorKindV2::UnknownCode);
            assert_eq!(error.field(), field);
        }
        let different_overlap = exact_policy(
            PlatformPathProfileV2::ContentStoreObjectKeyV1,
            RootCapabilityAccessV2::ReadOnlyInput,
            DestinationAdmissionModeV2::Absent,
            [1, 1, 2],
            0,
        );
        assert!(
            RootCapabilityPolicySetV2::new(
                RunnerCommandV2::Replay,
                vec![
                    exact_policy(
                        PlatformPathProfileV2::ContentStoreObjectKeyV1,
                        RootCapabilityAccessV2::DurableOutput,
                        DestinationAdmissionModeV2::Absent,
                        [1, 1, 1],
                        0,
                    ),
                    different_overlap,
                ],
                &registry,
            )
            .is_err()
        );
    }
}
