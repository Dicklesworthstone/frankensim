//! Shared, declaration-only Runner V2 extension-registry primitives.
//!
//! This module owns bounded data and pure validation only. It does not load a
//! family, execute a callback, encode or decode an artifact, inspect a frame,
//! grant a capability, establish registry membership from a bare numeric ID, or
//! mint scientific, lifecycle, durability, admission, or authority claims.

use crate::catalog::{ArtifactRoleV2, LogicalExtentAxisV2, LogicalUnitV2};
use crate::construction::{
    ConstructionClosedSemanticV2, ConstructionErrorKindV2, ConstructionErrorV2,
    ConstructionFixedObservationV2, ConstructionObservedDataClassV2, ConstructionObservedV2,
};
use crate::identity::{NoClaimScopeRootV1, RunnerLimitsRootV2};
use crate::limits::RunnerLimitsV2;
use crate::value::{RationalV2, StableTokenV2, UnitV2};
use core::num::NonZeroU16;
use fs_blake3::{ContentHash, hash_domain};
use std::collections::BTreeSet;

/// Domain for the first canonical, non-wire base extension-registry
/// projection.
pub const BASE_EXTENSION_REGISTRY_PROJECTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-extension-registry-projection.v1";

/// Domain for the first canonical, non-wire logical-extent projection.
pub const LOGICAL_EXTENT_PROJECTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.logical-extent-projection.v1";

/// Absolute pre-allocation cap for one registered axis's allowed-unit rows.
pub const LOGICAL_AXIS_ALLOWED_UNITS_MAX_V2: usize = 4_096;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ExtensionConstructionObservationV2 {
    DimensionMismatch,
}

impl ConstructionClosedSemanticV2 for ExtensionConstructionObservationV2 {
    fn construction_stable_name(&self) -> &'static str {
        match self {
            Self::DimensionMismatch => "dimension mismatch",
        }
    }
}

/// Fixed Runner V2 artifact codecs.
///
/// The catalog is declaration-only. In particular, this type provides no
/// encoder, decoder, frame validator, checksum, or artifact-inventory
/// behavior.
///
/// ```compile_fail,E0599
/// use fs_evidence_runner::ArtifactCodecIdV2;
///
/// let _encoded = ArtifactCodecIdV2::Identity.encode(b"payload");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum ArtifactCodecIdV2 {
    /// Encoded bytes are the identity representation.
    Identity = 0,
    /// A Zstandard frame under the separately owned V1 byte contract.
    ZstdFrameV1 = 1,
}

impl ArtifactCodecIdV2 {
    /// Exact catalog order.
    pub const ALL: [Self; 2] = [Self::Identity, Self::ZstdFrameV1];

    /// Resolve an exact fixed codec code.
    pub fn from_code(code: u16) -> Result<Self, ConstructionErrorV2> {
        match code {
            0 => Ok(Self::Identity),
            1 => Ok(Self::ZstdFrameV1),
            _ => Err(refusal(
                ConstructionErrorKindV2::UnknownCode,
                "artifact_codec.code",
                "fixed codec code 0 or 1",
                code,
            )),
        }
    }

    /// Exact fixed wire code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Stable codec name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Identity => "identity",
            Self::ZstdFrameV1 => "zstd-frame-v1",
        }
    }
}

/// One registered artifact-role descriptor.
///
/// Private fields prevent a syntactic `RegisteredFamilyRole(id)` value from
/// claiming semantic membership without lookup in a checked
/// [`BaseExtensionRegistryProjectionV2`].
///
/// ```compile_fail,E0616
/// use fs_evidence_runner::extension::RegisteredArtifactRoleDescriptorV2;
///
/// fn expose_unchecked_id(descriptor: &RegisteredArtifactRoleDescriptorV2) {
///     let _ = descriptor.id;
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct RegisteredArtifactRoleDescriptorV2 {
    id: NonZeroU16,
    name: StableTokenV2,
    owner: StableTokenV2,
    no_claim_scope: NoClaimScopeRootV1,
}

impl RegisteredArtifactRoleDescriptorV2 {
    /// Construct a bounded, non-executable descriptor.
    pub fn new(
        role: ArtifactRoleV2,
        name: StableTokenV2,
        owner: StableTokenV2,
        no_claim_scope: NoClaimScopeRootV1,
    ) -> Result<Self, ConstructionErrorV2> {
        Ok(Self {
            id: require_registered_variant(
                role.registered_id(),
                "extension.artifact_role.id",
                "ArtifactRoleV2::RegisteredFamilyRole",
            )?,
            name: require_globally_namespaced(name, "extension.artifact_role.name")?,
            owner,
            no_claim_scope,
        })
    }

    /// Namespace-local nonzero identifier.
    #[must_use]
    pub const fn id(&self) -> u16 {
        self.id.get()
    }

    /// Syntactic catalog value corresponding to this registered descriptor.
    #[must_use]
    pub const fn role(&self) -> ArtifactRoleV2 {
        ArtifactRoleV2::RegisteredFamilyRole(self.id)
    }

    /// Globally namespaced stable name.
    #[must_use]
    pub fn name(&self) -> &StableTokenV2 {
        &self.name
    }

    /// Sole semantic owner declaration.
    #[must_use]
    pub fn owner(&self) -> &StableTokenV2 {
        &self.owner
    }

    /// Explicit no-claim boundary.
    #[must_use]
    pub const fn no_claim_scope(&self) -> &NoClaimScopeRootV1 {
        &self.no_claim_scope
    }
}

/// One registered logical-unit descriptor.
///
/// This is registry data, not a conversion row. Conversion factors are scoped
/// to an owning logical axis.
///
/// Typed namespaces cannot be cross-substituted:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::{
///     ArtifactRoleV2, RegisteredLogicalUnitDescriptorV2,
/// };
///
/// let wrong_namespace = ArtifactRoleV2::from_tag(8, Some(7)).unwrap();
/// let _ = RegisteredLogicalUnitDescriptorV2::new(
///     wrong_namespace,
///     todo!(),
///     todo!(),
///     todo!(),
/// );
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct RegisteredLogicalUnitDescriptorV2 {
    id: NonZeroU16,
    name: StableTokenV2,
    owner: StableTokenV2,
    no_claim_scope: NoClaimScopeRootV1,
}

impl RegisteredLogicalUnitDescriptorV2 {
    /// Construct a bounded, non-executable descriptor.
    pub fn new(
        unit: LogicalUnitV2,
        name: StableTokenV2,
        owner: StableTokenV2,
        no_claim_scope: NoClaimScopeRootV1,
    ) -> Result<Self, ConstructionErrorV2> {
        Ok(Self {
            id: require_registered_variant(
                unit.registered_id(),
                "extension.logical_unit.id",
                "LogicalUnitV2::RegisteredUnit",
            )?,
            name: require_globally_namespaced(name, "extension.logical_unit.name")?,
            owner,
            no_claim_scope,
        })
    }

    /// Namespace-local nonzero identifier.
    #[must_use]
    pub const fn id(&self) -> u16 {
        self.id.get()
    }

    /// Syntactic catalog value corresponding to this registered descriptor.
    #[must_use]
    pub const fn unit(&self) -> LogicalUnitV2 {
        LogicalUnitV2::RegisteredUnit(self.id)
    }

    /// Globally namespaced stable name.
    #[must_use]
    pub fn name(&self) -> &StableTokenV2 {
        &self.name
    }

    /// Sole semantic owner declaration.
    #[must_use]
    pub fn owner(&self) -> &StableTokenV2 {
        &self.owner
    }

    /// Explicit no-claim boundary.
    #[must_use]
    pub const fn no_claim_scope(&self) -> &NoClaimScopeRootV1 {
        &self.no_claim_scope
    }
}

/// One unit admitted by an owning logical axis.
///
/// `scale_to_canonical` means that one value in `unit`, multiplied by this
/// positive exact scale, equals the value in the axis's canonical unit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LogicalUnitScaleToCanonicalV2 {
    unit: LogicalUnitV2,
    scale_to_canonical: RationalV2,
}

impl LogicalUnitScaleToCanonicalV2 {
    /// Admit raw scale parts only when they are already canonical and
    /// strictly positive.
    ///
    /// This is the presentation-boundary constructor. It deliberately refuses
    /// a reducible row such as `2/2` instead of silently normalizing registry
    /// data whose exact bytes are part of the projection identity.
    pub fn from_canonical_parts(
        unit: LogicalUnitV2,
        numerator: i128,
        denominator: u128,
    ) -> Result<Self, ConstructionErrorV2> {
        let scale_to_canonical =
            RationalV2::from_canonical_parts(numerator, denominator).map_err(|_| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "extension.axis.allowed_unit.scale_to_canonical",
                    "a positive exact rational already in lowest terms",
                    ConstructionObservedV2::signed_unsigned_pair(numerator, denominator),
                )
            })?;
        Self::new(unit, scale_to_canonical)
    }

    /// Construct one positive, already normalized scale row.
    pub fn new(
        unit: LogicalUnitV2,
        scale_to_canonical: RationalV2,
    ) -> Result<Self, ConstructionErrorV2> {
        if !scale_to_canonical.is_positive() {
            return Err(refusal(
                ConstructionErrorKindV2::OutOfRange,
                "extension.axis.allowed_unit.scale_to_canonical",
                "a positive exact reduced rational",
                ConstructionObservedV2::signed_unsigned_pair(
                    scale_to_canonical.numerator(),
                    scale_to_canonical.denominator(),
                ),
            ));
        }
        Ok(Self {
            unit,
            scale_to_canonical,
        })
    }

    /// Admitted unit.
    #[must_use]
    pub const fn unit(self) -> LogicalUnitV2 {
        self.unit
    }

    /// Exact positive scale into the owning axis's canonical unit.
    #[must_use]
    pub const fn scale_to_canonical(self) -> RationalV2 {
        self.scale_to_canonical
    }
}

/// One registered logical-extent-axis descriptor.
///
/// The allowed-unit vector is canonical, nonempty, and duplicate-free. The
/// canonical unit occurs exactly once with scale `1/1`.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct RegisteredLogicalExtentAxisDescriptorV2 {
    id: NonZeroU16,
    name: StableTokenV2,
    owner: StableTokenV2,
    no_claim_scope: NoClaimScopeRootV1,
    canonical_unit: LogicalUnitV2,
    allowed_units: Box<[LogicalUnitScaleToCanonicalV2]>,
}

impl RegisteredLogicalExtentAxisDescriptorV2 {
    /// Construct one bounded, normalized registered-axis descriptor.
    pub fn new(
        axis: LogicalExtentAxisV2,
        name: StableTokenV2,
        owner: StableTokenV2,
        no_claim_scope: NoClaimScopeRootV1,
        canonical_unit: LogicalUnitV2,
        allowed_units: &[LogicalUnitScaleToCanonicalV2],
    ) -> Result<Self, ConstructionErrorV2> {
        let id = require_registered_variant(
            axis.registered_id(),
            "extension.logical_axis.id",
            "LogicalExtentAxisV2::RegisteredAxis",
        )?;
        let name = require_globally_namespaced(name, "extension.logical_axis.name")?;
        if allowed_units.len() > LOGICAL_AXIS_ALLOWED_UNITS_MAX_V2 {
            return Err(refusal(
                ConstructionErrorKindV2::TooLarge,
                "extension.logical_axis.allowed_units",
                "at most 4096 allowed-unit rows before cloning",
                allowed_units.len(),
            ));
        }
        if allowed_units.is_empty() {
            return Err(refusal(
                ConstructionErrorKindV2::Missing,
                "extension.logical_axis.allowed_units",
                "a nonempty allowed-unit set",
                0,
            ));
        }

        let mut seen = BTreeSet::new();
        for row in allowed_units {
            if !seen.insert(row.unit()) {
                return Err(refusal(
                    ConstructionErrorKindV2::Duplicate,
                    "extension.logical_axis.allowed_units.unit",
                    "one row per typed logical unit",
                    ConstructionObservedV2::tag_and_optional_id(
                        row.unit().tag(),
                        row.unit().registered_id(),
                    ),
                ));
            }
        }

        let canonical_rows = allowed_units
            .iter()
            .filter(|row| row.unit() == canonical_unit)
            .collect::<Vec<_>>();
        if canonical_rows.len() != 1 {
            return Err(refusal(
                if canonical_rows.is_empty() {
                    ConstructionErrorKindV2::Missing
                } else {
                    ConstructionErrorKindV2::Duplicate
                },
                "extension.logical_axis.canonical_unit",
                "exactly one canonical-unit row",
                canonical_rows.len(),
            ));
        }
        if canonical_rows[0].scale_to_canonical() != rational_one() {
            let canonical_scale = canonical_rows[0].scale_to_canonical();
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "extension.logical_axis.canonical_scale",
                "the exact rational 1/1",
                ConstructionObservedV2::signed_unsigned_pair(
                    canonical_scale.numerator(),
                    canonical_scale.denominator(),
                ),
            ));
        }

        let mut canonical = allowed_units.to_vec();
        canonical.sort_by_key(|row| row.unit());
        Ok(Self {
            id,
            name,
            owner,
            no_claim_scope,
            canonical_unit,
            allowed_units: canonical.into_boxed_slice(),
        })
    }

    /// Namespace-local nonzero identifier.
    #[must_use]
    pub const fn id(&self) -> u16 {
        self.id.get()
    }

    /// Syntactic catalog value corresponding to this registered descriptor.
    #[must_use]
    pub const fn axis(&self) -> LogicalExtentAxisV2 {
        LogicalExtentAxisV2::RegisteredAxis(self.id)
    }

    /// Globally namespaced stable name.
    #[must_use]
    pub fn name(&self) -> &StableTokenV2 {
        &self.name
    }

    /// Sole semantic owner declaration.
    #[must_use]
    pub fn owner(&self) -> &StableTokenV2 {
        &self.owner
    }

    /// Explicit no-claim boundary.
    #[must_use]
    pub const fn no_claim_scope(&self) -> &NoClaimScopeRootV1 {
        &self.no_claim_scope
    }

    /// Canonical logical unit for this axis.
    #[must_use]
    pub const fn canonical_unit(&self) -> LogicalUnitV2 {
        self.canonical_unit
    }

    /// Canonically ordered, nonempty allowed-unit rows.
    #[must_use]
    pub fn allowed_units(&self) -> &[LogicalUnitScaleToCanonicalV2] {
        &self.allowed_units
    }

    fn scale_for_unit(&self, unit: LogicalUnitV2) -> Option<RationalV2> {
        self.allowed_units
            .binary_search_by_key(&unit, |row| row.unit())
            .ok()
            .map(|index| self.allowed_units[index].scale_to_canonical())
    }
}

/// Exact logical extent shared by evaluator and artifact-inventory schemas.
///
/// Raw fields are private so a registered ID cannot bypass semantic registry
/// lookup.
///
/// ```compile_fail,E0616
/// use fs_evidence_runner::extension::LogicalExtentV2;
///
/// fn expose_unchecked_axis(extent: &LogicalExtentV2) {
///     let _ = extent.axis;
/// }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LogicalExtentV2 {
    axis: LogicalExtentAxisV2,
    value: u128,
    unit: LogicalUnitV2,
    root: ContentHash,
}

/// Exact ordered field schema for [`LogicalExtentV2`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum LogicalExtentFieldV1 {
    /// Typed logical-axis tag and optional registered ID.
    Axis = 1,
    /// Exact unsigned 128-bit magnitude.
    Value = 2,
    /// Typed logical-unit tag and optional registered ID.
    Unit = 3,
}

impl LogicalExtentFieldV1 {
    /// Exact semantic field order.
    pub const ALL: [Self; 3] = [Self::Axis, Self::Value, Self::Unit];

    /// Exact non-wire schema ordinal.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Stable field name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Axis => "axis",
            Self::Value => "value",
            Self::Unit => "unit",
        }
    }
}

impl LogicalExtentV2 {
    /// Construct a base-axis extent without a family registry.
    ///
    /// Registered axes or units always refuse on this registry-free path.
    pub fn try_new_base(
        axis: LogicalExtentAxisV2,
        value: u128,
        unit: LogicalUnitV2,
    ) -> Result<Self, ConstructionErrorV2> {
        if axis.registered_id().is_some() {
            return Err(refusal(
                ConstructionErrorKindV2::UnknownCode,
                "logical_extent.axis",
                "a base axis or a registered axis resolved by a checked registry",
                ConstructionObservedV2::tag_and_optional_id(axis.tag(), axis.registered_id()),
            ));
        }
        if unit.registered_id().is_some() {
            return Err(refusal(
                ConstructionErrorKindV2::UnknownCode,
                "logical_extent.unit",
                "a base unit on the registry-free path",
                ConstructionObservedV2::tag_and_optional_id(unit.tag(), unit.registered_id()),
            ));
        }
        if base_axis_scale(axis, unit).is_none() {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "logical_extent.axis_unit",
                "a unit admitted by the selected base axis",
                logical_axis_unit_observation(axis, unit),
            ));
        }
        Ok(Self {
            axis,
            value,
            unit,
            root: logical_extent_root(axis, value, unit),
        })
    }

    /// Exact logical axis.
    #[must_use]
    pub const fn axis(self) -> LogicalExtentAxisV2 {
        self.axis
    }

    /// Exact unsigned 128-bit magnitude.
    #[must_use]
    pub const fn value(self) -> u128 {
        self.value
    }

    /// Exact logical unit.
    #[must_use]
    pub const fn unit(self) -> LogicalUnitV2 {
        self.unit
    }

    /// Domain-separated canonical identity in exact axis/value/unit order.
    #[must_use]
    pub const fn semantic_root(&self) -> &ContentHash {
        &self.root
    }
}

/// Bounded, canonical projection of the three extension-registry categories
/// owned by the base schema.
///
/// Population and family sealing remain downstream. The projection proves only
/// that the presented descriptor data is internally canonical and within its
/// explicitly named caps.
///
/// Base and registered-extension close-capability IDs are nominally distinct:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::coverage::{
///     BaseCoverageCloseCapabilityIdV1,
///     BaseCoverageCloseRegisteredExtensionCapabilityIdV1,
/// };
///
/// fn extension_capability_id_rejects_base_capability_id(
///     _id: BaseCoverageCloseRegisteredExtensionCapabilityIdV1,
/// ) {
/// }
///
/// let base_id = BaseCoverageCloseCapabilityIdV1::new(1).unwrap();
/// extension_capability_id_rejects_base_capability_id(base_id);
/// ```
///
/// The reverse substitution is rejected independently:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::coverage::{
///     BaseCoverageCloseCapabilityIdV1,
///     BaseCoverageCloseRegisteredExtensionCapabilityIdV1,
/// };
///
/// fn base_capability_id_rejects_extension_capability_id(
///     _id: BaseCoverageCloseCapabilityIdV1,
/// ) {
/// }
///
/// let extension_id =
///     BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(1).unwrap();
/// base_capability_id_rejects_extension_capability_id(extension_id);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseExtensionRegistryProjectionV2 {
    limits_root: RunnerLimitsRootV2,
    artifact_roles: Box<[RegisteredArtifactRoleDescriptorV2]>,
    logical_units: Box<[RegisteredLogicalUnitDescriptorV2]>,
    logical_axes: Box<[RegisteredLogicalExtentAxisDescriptorV2]>,
    root: ContentHash,
}

impl BaseExtensionRegistryProjectionV2 {
    /// Validate and canonicalize the base-owned extension categories.
    ///
    /// Duplicate IDs and global names are detected in caller order before
    /// canonical sorting. Numeric IDs are namespace-local: the same nonzero
    /// value may occur once in each different typed category.
    pub fn try_new(
        limits: &RunnerLimitsV2,
        artifact_roles: &[RegisteredArtifactRoleDescriptorV2],
        logical_units: &[RegisteredLogicalUnitDescriptorV2],
        logical_axes: &[RegisteredLogicalExtentAxisDescriptorV2],
    ) -> Result<Self, ConstructionErrorV2> {
        require_count_cap(
            artifact_roles.len(),
            limits.artifact_roles_per_family(),
            "extension.artifact_roles",
        )?;
        require_count_cap(
            logical_units.len(),
            limits.registered_units_per_family(),
            "extension.logical_units",
        )?;
        require_count_cap(
            logical_axes.len(),
            limits.registered_extent_axes_per_family(),
            "extension.logical_axes",
        )?;

        detect_duplicate_ids(
            artifact_roles
                .iter()
                .map(RegisteredArtifactRoleDescriptorV2::id),
            "extension.artifact_roles.id",
        )?;
        detect_duplicate_ids(
            logical_units
                .iter()
                .map(RegisteredLogicalUnitDescriptorV2::id),
            "extension.logical_units.id",
        )?;
        detect_duplicate_ids(
            logical_axes
                .iter()
                .map(RegisteredLogicalExtentAxisDescriptorV2::id),
            "extension.logical_axes.id",
        )?;

        let mut global_names = BTreeSet::new();
        for name in artifact_roles
            .iter()
            .map(RegisteredArtifactRoleDescriptorV2::name)
            .chain(
                logical_units
                    .iter()
                    .map(RegisteredLogicalUnitDescriptorV2::name),
            )
            .chain(
                logical_axes
                    .iter()
                    .map(RegisteredLogicalExtentAxisDescriptorV2::name),
            )
        {
            if !global_names.insert(name.as_str()) {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Duplicate,
                    "extension.global_name",
                    "one globally unambiguous namespaced name across all categories",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
        }

        let registered_unit_ids = logical_units
            .iter()
            .map(RegisteredLogicalUnitDescriptorV2::id)
            .collect::<BTreeSet<_>>();
        for axis in logical_axes {
            require_count_cap(
                axis.allowed_units().len(),
                limits.generic_array_items(),
                "extension.logical_axis.allowed_units",
            )?;
            for row in axis.allowed_units() {
                if let Some(id) = row.unit().registered_id()
                    && !registered_unit_ids.contains(&id)
                {
                    return Err(refusal(
                        ConstructionErrorKindV2::UnknownCode,
                        "extension.logical_axis.allowed_units.registered_unit",
                        "a registered logical unit present in this exact projection",
                        id,
                    ));
                }
            }
        }

        let mut artifact_roles = artifact_roles.to_vec();
        artifact_roles.sort_by_key(RegisteredArtifactRoleDescriptorV2::id);
        let mut logical_units = logical_units.to_vec();
        logical_units.sort_by_key(RegisteredLogicalUnitDescriptorV2::id);
        let mut logical_axes = logical_axes.to_vec();
        logical_axes.sort_by_key(RegisteredLogicalExtentAxisDescriptorV2::id);
        let limits_root = limits.semantic_root();
        let root =
            extension_registry_root(&limits_root, &artifact_roles, &logical_units, &logical_axes);

        Ok(Self {
            limits_root,
            artifact_roles: artifact_roles.into_boxed_slice(),
            logical_units: logical_units.into_boxed_slice(),
            logical_axes: logical_axes.into_boxed_slice(),
            root,
        })
    }

    /// Reconstruct an exact projection and refuse missing, extra, or mutated
    /// semantic descriptor data. Caller ordering may differ because the
    /// projection is explicitly permutation-invariant.
    pub fn reconstruct_exact(
        &self,
        limits: &RunnerLimitsV2,
        artifact_roles: &[RegisteredArtifactRoleDescriptorV2],
        logical_units: &[RegisteredLogicalUnitDescriptorV2],
        logical_axes: &[RegisteredLogicalExtentAxisDescriptorV2],
    ) -> Result<Self, ConstructionErrorV2> {
        let limits_root = limits.semantic_root();
        if limits_root != self.limits_root {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "extension.registry.limits_root",
                "the exact admitted Runner-limit semantic root",
                ConstructionObservedDataClassV2::CapabilityOrResource,
            ));
        }
        // Exact-set reconstruction owns missing/extra classification. Perform
        // the complete permutation-insensitive category checks before
        // cross-category referential validation so, for example, an axis that
        // still names an omitted or substituted registered unit cannot turn a
        // precise Missing or Incompatible result into the lower-level
        // UnknownCode used by first-time construction.
        require_exact_category(
            "extension.registry.artifact_roles",
            self.artifact_roles.len(),
            artifact_roles.len(),
            exact_permutation(self.artifact_roles.as_ref(), artifact_roles),
        )?;
        require_exact_category(
            "extension.registry.logical_units",
            self.logical_units.len(),
            logical_units.len(),
            exact_permutation(self.logical_units.as_ref(), logical_units),
        )?;
        require_exact_category(
            "extension.registry.logical_axes",
            self.logical_axes.len(),
            logical_axes.len(),
            exact_permutation(self.logical_axes.as_ref(), logical_axes),
        )?;
        let reconstructed = Self::try_new(limits, artifact_roles, logical_units, logical_axes)?;
        require_exact_category(
            "extension.registry.artifact_roles",
            self.artifact_roles.len(),
            reconstructed.artifact_roles.len(),
            self.artifact_roles.as_ref() == reconstructed.artifact_roles.as_ref(),
        )?;
        require_exact_category(
            "extension.registry.logical_units",
            self.logical_units.len(),
            reconstructed.logical_units.len(),
            self.logical_units.as_ref() == reconstructed.logical_units.as_ref(),
        )?;
        require_exact_category(
            "extension.registry.logical_axes",
            self.logical_axes.len(),
            reconstructed.logical_axes.len(),
            self.logical_axes.as_ref() == reconstructed.logical_axes.as_ref(),
        )?;
        if reconstructed != *self {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "extension.registry.exact_set",
                "the exact canonical descriptor set and semantic root",
                ConstructionObservedDataClassV2::CapabilityOrResource,
            ));
        }
        Ok(reconstructed)
    }

    /// Canonically ordered registered artifact roles.
    #[must_use]
    pub fn artifact_roles(&self) -> &[RegisteredArtifactRoleDescriptorV2] {
        &self.artifact_roles
    }

    /// Exact admitted Runner-limit vector bound into this projection.
    #[must_use]
    pub const fn limits_root(&self) -> &RunnerLimitsRootV2 {
        &self.limits_root
    }

    /// Canonically ordered registered logical units.
    #[must_use]
    pub fn logical_units(&self) -> &[RegisteredLogicalUnitDescriptorV2] {
        &self.logical_units
    }

    /// Canonically ordered registered logical axes.
    #[must_use]
    pub fn logical_axes(&self) -> &[RegisteredLogicalExtentAxisDescriptorV2] {
        &self.logical_axes
    }

    /// Domain-separated semantic projection root.
    #[must_use]
    pub const fn root(&self) -> &ContentHash {
        &self.root
    }

    /// Resolve one artifact-role ID in its own namespace.
    pub fn artifact_role(
        &self,
        id: u16,
    ) -> Result<&RegisteredArtifactRoleDescriptorV2, ConstructionErrorV2> {
        lookup_by_id(
            &self.artifact_roles,
            id,
            RegisteredArtifactRoleDescriptorV2::id,
            "extension.artifact_role.lookup",
        )
    }

    /// Resolve one logical-unit ID in its own namespace.
    pub fn logical_unit(
        &self,
        id: u16,
    ) -> Result<&RegisteredLogicalUnitDescriptorV2, ConstructionErrorV2> {
        lookup_by_id(
            &self.logical_units,
            id,
            RegisteredLogicalUnitDescriptorV2::id,
            "extension.logical_unit.lookup",
        )
    }

    /// Resolve one logical-axis ID in its own namespace.
    pub fn logical_axis(
        &self,
        id: u16,
    ) -> Result<&RegisteredLogicalExtentAxisDescriptorV2, ConstructionErrorV2> {
        lookup_by_id(
            &self.logical_axes,
            id,
            RegisteredLogicalExtentAxisDescriptorV2::id,
            "extension.logical_axis.lookup",
        )
    }

    /// Admit a base or registered logical extent against this exact registry.
    pub fn try_extent(
        &self,
        axis: LogicalExtentAxisV2,
        value: u128,
        unit: LogicalUnitV2,
    ) -> Result<LogicalExtentV2, ConstructionErrorV2> {
        let scale = self.axis_scale(axis, unit)?;
        let _ = scale;
        Ok(LogicalExtentV2 {
            axis,
            value,
            unit,
            root: logical_extent_root(axis, value, unit),
        })
    }

    /// Convert an admitted extent within its owning axis.
    ///
    /// The result remains an exact `u128`; a fractional or overflowing target
    /// magnitude refuses rather than rounding or saturating.
    pub fn convert_extent(
        &self,
        extent: LogicalExtentV2,
        target_unit: LogicalUnitV2,
    ) -> Result<LogicalExtentV2, ConstructionErrorV2> {
        let source_scale = self.axis_scale(extent.axis, extent.unit)?;
        let target_scale = self.axis_scale(extent.axis, target_unit)?;
        let value = checked_scale_u128(extent.value, source_scale, target_scale)?;
        Ok(LogicalExtentV2 {
            axis: extent.axis,
            value,
            unit: target_unit,
            root: logical_extent_root(extent.axis, value, target_unit),
        })
    }

    fn axis_scale(
        &self,
        axis: LogicalExtentAxisV2,
        unit: LogicalUnitV2,
    ) -> Result<RationalV2, ConstructionErrorV2> {
        if let Some(axis_id) = axis.registered_id() {
            let descriptor = self.logical_axis(axis_id)?;
            if let Some(unit_id) = unit.registered_id() {
                self.logical_unit(unit_id)?;
            }
            descriptor.scale_for_unit(unit).ok_or_else(|| {
                refusal(
                    ConstructionErrorKindV2::Incompatible,
                    "logical_extent.axis_unit",
                    "a unit explicitly admitted by the registered axis",
                    logical_axis_unit_observation(axis, unit),
                )
            })
        } else {
            if let Some(unit_id) = unit.registered_id() {
                return Err(refusal(
                    ConstructionErrorKindV2::Incompatible,
                    "logical_extent.base_axis_registered_unit",
                    "a frozen base unit admitted by the selected base axis",
                    unit_id,
                ));
            }
            base_axis_scale(axis, unit).ok_or_else(|| {
                refusal(
                    ConstructionErrorKindV2::Incompatible,
                    "logical_extent.axis_unit",
                    "a unit admitted by the selected base axis",
                    logical_axis_unit_observation(axis, unit),
                )
            })
        }
    }
}

fn require_exact_category(
    field: &'static str,
    expected_count: usize,
    observed_count: usize,
    exact_rows_match: bool,
) -> Result<(), ConstructionErrorV2> {
    if observed_count < expected_count {
        return Err(refusal(
            ConstructionErrorKindV2::Missing,
            field,
            "the exact category row set",
            count_pair_observation(observed_count, expected_count),
        ));
    }
    if observed_count > expected_count {
        return Err(refusal(
            ConstructionErrorKindV2::Unexpected,
            field,
            "no row beyond the exact category row set",
            count_pair_observation(observed_count, expected_count),
        ));
    }
    if !exact_rows_match {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            field,
            "the exact canonical category descriptors",
            observed_count,
        ));
    }
    Ok(())
}

fn exact_permutation<T: PartialEq>(expected: &[T], observed: &[T]) -> bool {
    expected.len() == observed.len() && expected.iter().all(|row| observed.contains(row))
}

/// Return the exact normalized scale ratio between dimension-compatible
/// [`UnitV2`] values.
///
/// Multiplying a value expressed in `source` by the returned ratio produces
/// the corresponding value in `target`.
pub fn normalized_unit_scale_ratio_v2(
    source: UnitV2,
    target: UnitV2,
) -> Result<RationalV2, ConstructionErrorV2> {
    if source.exponents() != target.exponents() {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "unit_conversion.dimensions",
            "identical seven-axis SI exponents",
            ConstructionObservedV2::closed(&ExtensionConstructionObservationV2::DimensionMismatch),
        ));
    }
    checked_positive_ratio(source.scale(), target.scale())
}

/// Convert one exact rational quantity between dimension-compatible units.
///
/// All normalization and overflow checks are exact. No floating point,
/// rounding, saturation, or pairwise conversion registry is used.
pub fn convert_rational_quantity_v2(
    value: RationalV2,
    source: UnitV2,
    target: UnitV2,
) -> Result<RationalV2, ConstructionErrorV2> {
    let ratio = normalized_unit_scale_ratio_v2(source, target)?;
    checked_rational_product(value, ratio)
}

fn base_axis_scale(axis: LogicalExtentAxisV2, unit: LogicalUnitV2) -> Option<RationalV2> {
    let admitted = match axis {
        LogicalExtentAxisV2::Payload => unit == LogicalUnitV2::LogicalBytes,
        LogicalExtentAxisV2::Records => unit == LogicalUnitV2::Records,
        LogicalExtentAxisV2::Rows => unit == LogicalUnitV2::Rows,
        LogicalExtentAxisV2::Elements => unit == LogicalUnitV2::Elements,
        LogicalExtentAxisV2::Samples => unit == LogicalUnitV2::Samples,
        LogicalExtentAxisV2::Iterations => unit == LogicalUnitV2::Iterations,
        LogicalExtentAxisV2::Operations => unit == LogicalUnitV2::Operations,
        LogicalExtentAxisV2::Cycles => unit == LogicalUnitV2::Cycles,
        LogicalExtentAxisV2::Duration => {
            unit == LogicalUnitV2::Nanoseconds || unit == LogicalUnitV2::Seconds
        }
        LogicalExtentAxisV2::RegisteredAxis(_) => false,
    };
    if !admitted {
        return None;
    }
    Some(
        if axis == LogicalExtentAxisV2::Duration && unit == LogicalUnitV2::Seconds {
            RationalV2::from_canonical_parts(1_000_000_000, 1)
                .expect("the frozen seconds-to-nanoseconds scale is canonical")
        } else {
            rational_one()
        },
    )
}

fn logical_extent_root(axis: LogicalExtentAxisV2, value: u128, unit: LogicalUnitV2) -> ContentHash {
    let mut bytes = Vec::with_capacity(26);
    bytes.extend_from_slice(&axis.tag().to_be_bytes());
    match axis.registered_id() {
        Some(id) => {
            bytes.push(1);
            bytes.extend_from_slice(&id.to_be_bytes());
        }
        None => bytes.push(0),
    }
    bytes.extend_from_slice(&value.to_be_bytes());
    bytes.extend_from_slice(&unit.tag().to_be_bytes());
    match unit.registered_id() {
        Some(id) => {
            bytes.push(1);
            bytes.extend_from_slice(&id.to_be_bytes());
        }
        None => bytes.push(0),
    }
    hash_domain(LOGICAL_EXTENT_PROJECTION_DOMAIN_V1, &bytes)
}

fn checked_scale_u128(
    value: u128,
    source_scale: RationalV2,
    target_scale: RationalV2,
) -> Result<u128, ConstructionErrorV2> {
    let source_numerator = positive_numerator(source_scale, "logical_extent.source_scale")?;
    let target_numerator = positive_numerator(target_scale, "logical_extent.target_scale")?;
    let mut numerators = [value, source_numerator, target_scale.denominator()];
    let mut denominators = [source_scale.denominator(), target_numerator];
    cancel_factors(&mut numerators, &mut denominators);
    let denominator = checked_product(&denominators, "logical_extent.conversion.denominator")?;
    if denominator != 1 {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "logical_extent.conversion",
            "an exact integral u128 target extent",
            denominator,
        ));
    }
    checked_product(&numerators, "logical_extent.conversion.value")
}

fn checked_positive_ratio(
    source: RationalV2,
    target: RationalV2,
) -> Result<RationalV2, ConstructionErrorV2> {
    let source_numerator = positive_numerator(source, "unit_conversion.source_scale")?;
    let target_numerator = positive_numerator(target, "unit_conversion.target_scale")?;
    let mut numerators = [source_numerator, target.denominator()];
    let mut denominators = [source.denominator(), target_numerator];
    cancel_factors(&mut numerators, &mut denominators);
    let numerator = checked_product(&numerators, "unit_conversion.ratio_numerator")?;
    let denominator = checked_product(&denominators, "unit_conversion.ratio_denominator")?;
    let numerator = i128::try_from(numerator).map_err(|_| {
        refusal(
            ConstructionErrorKindV2::ArithmeticOverflow,
            "unit_conversion.ratio_numerator",
            "a positive i128 numerator",
            numerator,
        )
    })?;
    RationalV2::from_canonical_parts(numerator, denominator).map_err(|_| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "unit_conversion.ratio",
            "one canonical exact rational",
            ConstructionObservedV2::signed_unsigned_pair(numerator, denominator),
        )
    })
}

fn checked_rational_product(
    left: RationalV2,
    right: RationalV2,
) -> Result<RationalV2, ConstructionErrorV2> {
    let negative = left.numerator().is_negative() ^ right.numerator().is_negative();
    let mut numerators = [
        left.numerator().unsigned_abs(),
        right.numerator().unsigned_abs(),
    ];
    let mut denominators = [left.denominator(), right.denominator()];
    cancel_factors(&mut numerators, &mut denominators);
    let magnitude = checked_product(&numerators, "unit_conversion.value_numerator")?;
    let denominator = checked_product(&denominators, "unit_conversion.value_denominator")?;
    let signed = signed_from_magnitude_checked(negative, magnitude)?;
    RationalV2::from_canonical_parts(signed, denominator).map_err(|_| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "unit_conversion.value",
            "one canonical exact rational",
            ConstructionObservedV2::signed_unsigned_pair(signed, denominator),
        )
    })
}

fn positive_numerator(value: RationalV2, field: &'static str) -> Result<u128, ConstructionErrorV2> {
    if !value.is_positive() {
        return Err(refusal(
            ConstructionErrorKindV2::OutOfRange,
            field,
            "a positive exact rational",
            value.numerator(),
        ));
    }
    Ok(value.numerator().unsigned_abs())
}

fn cancel_factors<const N: usize, const D: usize>(
    numerators: &mut [u128; N],
    denominators: &mut [u128; D],
) {
    for numerator in numerators {
        for denominator in denominators.iter_mut() {
            let divisor = gcd_u128(*numerator, *denominator);
            *numerator /= divisor;
            *denominator /= divisor;
        }
    }
}

fn checked_product(factors: &[u128], field: &'static str) -> Result<u128, ConstructionErrorV2> {
    factors.iter().try_fold(1_u128, |product, factor| {
        product.checked_mul(*factor).ok_or_else(|| {
            refusal(
                ConstructionErrorKindV2::ArithmeticOverflow,
                field,
                "a product representable as u128",
                ConstructionObservedV2::fixed(ConstructionFixedObservationV2::Overflow),
            )
        })
    })
}

fn signed_from_magnitude_checked(
    negative: bool,
    magnitude: u128,
) -> Result<i128, ConstructionErrorV2> {
    const I128_MIN_MAGNITUDE: u128 = 1_u128 << 127;
    if negative {
        if magnitude == I128_MIN_MAGNITUDE {
            return Ok(i128::MIN);
        }
        return i128::try_from(magnitude).map(|value| -value).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::ArithmeticOverflow,
                "unit_conversion.value_numerator",
                "a negative i128 numerator",
                magnitude,
            )
        });
    }
    i128::try_from(magnitude).map_err(|_| {
        refusal(
            ConstructionErrorKindV2::ArithmeticOverflow,
            "unit_conversion.value_numerator",
            "a nonnegative i128 numerator",
            magnitude,
        )
    })
}

fn extension_registry_root(
    limits_root: &RunnerLimitsRootV2,
    artifact_roles: &[RegisteredArtifactRoleDescriptorV2],
    logical_units: &[RegisteredLogicalUnitDescriptorV2],
    logical_axes: &[RegisteredLogicalExtentAxisDescriptorV2],
) -> ContentHash {
    let mut bytes = Vec::new();
    push_str(&mut bytes, limits_root.domain());
    bytes.extend_from_slice(limits_root.bytes());
    push_count(&mut bytes, artifact_roles.len());
    for row in artifact_roles {
        bytes.extend_from_slice(&row.id().to_be_bytes());
        push_str(&mut bytes, row.name().as_str());
        push_str(&mut bytes, row.owner().as_str());
        push_no_claim(&mut bytes, row.no_claim_scope());
    }
    push_count(&mut bytes, logical_units.len());
    for row in logical_units {
        bytes.extend_from_slice(&row.id().to_be_bytes());
        push_str(&mut bytes, row.name().as_str());
        push_str(&mut bytes, row.owner().as_str());
        push_no_claim(&mut bytes, row.no_claim_scope());
    }
    push_count(&mut bytes, logical_axes.len());
    for row in logical_axes {
        bytes.extend_from_slice(&row.id().to_be_bytes());
        push_str(&mut bytes, row.name().as_str());
        push_str(&mut bytes, row.owner().as_str());
        push_no_claim(&mut bytes, row.no_claim_scope());
        push_logical_unit(&mut bytes, row.canonical_unit());
        push_count(&mut bytes, row.allowed_units().len());
        for allowed in row.allowed_units() {
            push_logical_unit(&mut bytes, allowed.unit());
            bytes.extend_from_slice(&allowed.scale_to_canonical().numerator().to_be_bytes());
            bytes.extend_from_slice(&allowed.scale_to_canonical().denominator().to_be_bytes());
        }
    }
    hash_domain(BASE_EXTENSION_REGISTRY_PROJECTION_DOMAIN_V1, &bytes)
}

fn push_count(bytes: &mut Vec<u8>, count: usize) {
    bytes.extend_from_slice(
        &u32::try_from(count)
            .expect("every admitted registry count fits u32")
            .to_be_bytes(),
    );
}

fn push_str(bytes: &mut Vec<u8>, value: &str) {
    push_count(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_no_claim(bytes: &mut Vec<u8>, value: &NoClaimScopeRootV1) {
    bytes.extend_from_slice(&value.role().code().to_be_bytes());
    push_str(bytes, value.domain());
    bytes.extend_from_slice(value.bytes());
}

fn push_logical_unit(bytes: &mut Vec<u8>, value: LogicalUnitV2) {
    bytes.extend_from_slice(&value.tag().to_be_bytes());
    bytes.extend_from_slice(&value.registered_id().unwrap_or(0).to_be_bytes());
}

fn require_registered_variant(
    id: Option<u16>,
    field: &'static str,
    expected_variant: &'static str,
) -> Result<NonZeroU16, ConstructionErrorV2> {
    let Some(id) = id else {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            field,
            expected_variant,
            ConstructionObservedV2::fixed(ConstructionFixedObservationV2::Absent),
        ));
    };
    NonZeroU16::new(id).ok_or_else(|| {
        refusal(
            ConstructionErrorKindV2::Zero,
            field,
            "a nonzero namespace-local u16 identifier",
            0,
        )
    })
}

fn require_globally_namespaced(
    name: StableTokenV2,
    field: &'static str,
) -> Result<StableTokenV2, ConstructionErrorV2> {
    if name.as_str().split('.').count() < 3 {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            field,
            "a globally namespaced stable token with at least three dot-separated segments",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    Ok(name)
}

fn require_count_cap(
    observed: usize,
    maximum: u32,
    field: &'static str,
) -> Result<(), ConstructionErrorV2> {
    let Ok(maximum) = usize::try_from(maximum) else {
        // If the u32 cap does not fit usize, no collection representable on
        // this target can exceed it.
        return Ok(());
    };
    if observed > maximum {
        return Err(refusal(
            ConstructionErrorKindV2::TooLarge,
            field,
            "a collection within its semantically named RunnerLimitsV2 cap",
            count_pair_observation(observed, maximum),
        ));
    }
    Ok(())
}

fn detect_duplicate_ids(
    ids: impl IntoIterator<Item = u16>,
    field: &'static str,
) -> Result<(), ConstructionErrorV2> {
    let mut seen = BTreeSet::new();
    for id in ids {
        if !seen.insert(id) {
            return Err(refusal(
                ConstructionErrorKindV2::Duplicate,
                field,
                "one descriptor per namespace-local ID",
                id,
            ));
        }
    }
    Ok(())
}

fn lookup_by_id<'a, T>(
    rows: &'a [T],
    id: u16,
    id_of: impl Fn(&T) -> u16,
    field: &'static str,
) -> Result<&'a T, ConstructionErrorV2> {
    if id == 0 {
        return Err(refusal(
            ConstructionErrorKindV2::Zero,
            field,
            "a registered nonzero u16 ID present in this exact projection",
            0,
        ));
    }
    rows.binary_search_by_key(&id, id_of)
        .map(|index| &rows[index])
        .map_err(|_| {
            refusal(
                ConstructionErrorKindV2::UnknownCode,
                field,
                "a registered nonzero u16 ID present in this exact projection",
                id,
            )
        })
}

fn rational_one() -> RationalV2 {
    RationalV2::from_canonical_parts(1, 1).expect("1/1 is canonical")
}

fn gcd_u128(mut left: u128, mut right: u128) -> u128 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left
}

fn logical_axis_unit_observation(
    axis: LogicalExtentAxisV2,
    unit: LogicalUnitV2,
) -> ConstructionObservedV2 {
    ConstructionObservedV2::unsigned_quad(
        u64::from(axis.tag()),
        u64::from(axis.registered_id().unwrap_or(0)),
        u64::from(unit.tag()),
        u64::from(unit.registered_id().unwrap_or(0)),
    )
}

fn count_pair_observation(observed: usize, expected_or_maximum: usize) -> ConstructionObservedV2 {
    ConstructionObservedV2::unsigned_pair(
        u64::try_from(observed).expect("Runner collection counts fit u64"),
        u64::try_from(expected_or_maximum).expect("Runner collection limits fit u64"),
    )
}

fn refusal(
    kind: ConstructionErrorKindV2,
    field: &'static str,
    expected: &'static str,
    observed: impl Into<ConstructionObservedV2>,
) -> ConstructionErrorV2 {
    ConstructionErrorV2::new(kind, field, expected, observed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalog::{DigestRoleV2, RunProfileV2};

    fn token(value: impl Into<String>) -> StableTokenV2 {
        StableTokenV2::new(value).expect("valid stable token fixture")
    }

    fn no_claim(byte: u8) -> NoClaimScopeRootV1 {
        NoClaimScopeRootV1::parse_presented(
            DigestRoleV2::ClaimScope,
            NoClaimScopeRootV1::DESCRIPTOR.domain(),
            &format!("{byte:02x}").repeat(32),
        )
        .expect("valid presented no-claim fixture")
    }

    fn registered_role(id: u16) -> ArtifactRoleV2 {
        ArtifactRoleV2::from_tag(8, Some(id)).expect("valid registered role fixture")
    }

    fn registered_unit(id: u16) -> LogicalUnitV2 {
        LogicalUnitV2::from_tag(16, Some(id)).expect("valid registered unit fixture")
    }

    fn registered_axis(id: u16) -> LogicalExtentAxisV2 {
        LogicalExtentAxisV2::from_tag(10, Some(id)).expect("valid registered axis fixture")
    }

    fn role(id: u16, suffix: &str, claim_byte: u8) -> RegisteredArtifactRoleDescriptorV2 {
        RegisteredArtifactRoleDescriptorV2::new(
            registered_role(id),
            token(format!("org.example.role.{suffix}")),
            token("org.example.owner"),
            no_claim(claim_byte),
        )
        .expect("valid registered artifact-role fixture")
    }

    fn unit(id: u16, suffix: &str, claim_byte: u8) -> RegisteredLogicalUnitDescriptorV2 {
        RegisteredLogicalUnitDescriptorV2::new(
            registered_unit(id),
            token(format!("org.example.unit.{suffix}")),
            token("org.example.owner"),
            no_claim(claim_byte),
        )
        .expect("valid registered logical-unit fixture")
    }

    fn scale(
        unit: LogicalUnitV2,
        numerator: i128,
        denominator: u128,
    ) -> LogicalUnitScaleToCanonicalV2 {
        LogicalUnitScaleToCanonicalV2::new(
            unit,
            RationalV2::from_canonical_parts(numerator, denominator)
                .expect("canonical rational fixture"),
        )
        .expect("positive scale fixture")
    }

    fn axis(
        id: u16,
        suffix: &str,
        claim_byte: u8,
        canonical_unit: LogicalUnitV2,
        allowed_units: &[LogicalUnitScaleToCanonicalV2],
    ) -> RegisteredLogicalExtentAxisDescriptorV2 {
        RegisteredLogicalExtentAxisDescriptorV2::new(
            registered_axis(id),
            token(format!("org.example.axis.{suffix}")),
            token("org.example.owner"),
            no_claim(claim_byte),
            canonical_unit,
            allowed_units,
        )
        .expect("valid registered axis fixture")
    }

    #[test]
    fn artifact_codec_catalog_is_exact_closed_and_registry_independent() {
        assert_eq!(
            ArtifactCodecIdV2::ALL,
            [ArtifactCodecIdV2::Identity, ArtifactCodecIdV2::ZstdFrameV1]
        );
        assert_eq!(ArtifactCodecIdV2::Identity.code(), 0);
        assert_eq!(ArtifactCodecIdV2::Identity.name(), "identity");
        assert_eq!(ArtifactCodecIdV2::ZstdFrameV1.code(), 1);
        assert_eq!(ArtifactCodecIdV2::ZstdFrameV1.name(), "zstd-frame-v1");
        assert_eq!(
            ArtifactCodecIdV2::from_code(0),
            Ok(ArtifactCodecIdV2::Identity)
        );
        assert_eq!(
            ArtifactCodecIdV2::from_code(1),
            Ok(ArtifactCodecIdV2::ZstdFrameV1)
        );
        for code in [2, u16::MAX] {
            let error = ArtifactCodecIdV2::from_code(code).expect_err("unknown codec");
            assert_eq!(error.kind(), ConstructionErrorKindV2::UnknownCode);
            assert_eq!(error.field(), "artifact_codec.code");
        }

        let empty = BaseExtensionRegistryProjectionV2::try_new(
            &RunnerLimitsV2::base(RunProfileV2::Smoke),
            &[],
            &[],
            &[],
        )
        .expect("fixed codecs consume no family registry slots");
        assert!(empty.artifact_roles().is_empty());
        assert!(empty.logical_units().is_empty());
        assert!(empty.logical_axes().is_empty());
    }

    #[test]
    fn logical_extent_base_axis_unit_table_and_u128_extrema_are_exact() {
        let exact = [
            (LogicalExtentAxisV2::Payload, LogicalUnitV2::LogicalBytes),
            (LogicalExtentAxisV2::Records, LogicalUnitV2::Records),
            (LogicalExtentAxisV2::Rows, LogicalUnitV2::Rows),
            (LogicalExtentAxisV2::Elements, LogicalUnitV2::Elements),
            (LogicalExtentAxisV2::Samples, LogicalUnitV2::Samples),
            (LogicalExtentAxisV2::Iterations, LogicalUnitV2::Iterations),
            (LogicalExtentAxisV2::Operations, LogicalUnitV2::Operations),
            (LogicalExtentAxisV2::Cycles, LogicalUnitV2::Cycles),
            (LogicalExtentAxisV2::Duration, LogicalUnitV2::Nanoseconds),
            (LogicalExtentAxisV2::Duration, LogicalUnitV2::Seconds),
        ];
        for (axis, unit) in exact {
            for value in [0, 1, u128::MAX] {
                let extent =
                    LogicalExtentV2::try_new_base(axis, value, unit).expect("exact base cell");
                assert_eq!(extent.axis(), axis);
                assert_eq!(extent.value(), value);
                assert_eq!(extent.unit(), unit);
            }
        }

        let error =
            LogicalExtentV2::try_new_base(LogicalExtentAxisV2::Payload, 1, LogicalUnitV2::Rows)
                .expect_err("cross-axis unit");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);

        let registered_unit = LogicalUnitV2::from_tag(16, Some(7)).expect("syntactic ID");
        let registered_axis = LogicalExtentAxisV2::from_tag(10, Some(7)).expect("syntactic ID");
        assert_eq!(
            LogicalExtentV2::try_new_base(LogicalExtentAxisV2::Payload, 1, registered_unit)
                .expect_err("bare registered unit")
                .kind(),
            ConstructionErrorKindV2::UnknownCode
        );
        assert_eq!(
            LogicalExtentV2::try_new_base(registered_axis, 1, LogicalUnitV2::Count)
                .expect_err("bare registered axis")
                .kind(),
            ConstructionErrorKindV2::UnknownCode
        );
    }

    #[test]
    fn logical_extent_schema_and_root_bind_axis_value_unit_in_exact_order() {
        assert_eq!(
            LogicalExtentFieldV1::ALL,
            [
                LogicalExtentFieldV1::Axis,
                LogicalExtentFieldV1::Value,
                LogicalExtentFieldV1::Unit,
            ]
        );
        assert_eq!(
            LogicalExtentFieldV1::ALL.map(LogicalExtentFieldV1::code),
            [1, 2, 3]
        );
        assert_eq!(
            LogicalExtentFieldV1::ALL.map(LogicalExtentFieldV1::name),
            ["axis", "value", "unit"]
        );

        let base = LogicalExtentV2::try_new_base(
            LogicalExtentAxisV2::Duration,
            1,
            LogicalUnitV2::Nanoseconds,
        )
        .expect("base extent");
        let mut independent_bytes = Vec::new();
        independent_bytes.extend_from_slice(&9_u16.to_be_bytes());
        independent_bytes.push(0);
        independent_bytes.extend_from_slice(&1_u128.to_be_bytes());
        independent_bytes.extend_from_slice(&13_u16.to_be_bytes());
        independent_bytes.push(0);
        assert_eq!(
            *base.semantic_root(),
            hash_domain(LOGICAL_EXTENT_PROJECTION_DOMAIN_V1, &independent_bytes)
        );

        let value_mutated = LogicalExtentV2::try_new_base(
            LogicalExtentAxisV2::Duration,
            2,
            LogicalUnitV2::Nanoseconds,
        )
        .expect("value mutant");
        let unit_mutated =
            LogicalExtentV2::try_new_base(LogicalExtentAxisV2::Duration, 1, LogicalUnitV2::Seconds)
                .expect("unit mutant");
        assert_ne!(value_mutated.semantic_root(), base.semantic_root());
        assert_ne!(unit_mutated.semantic_root(), base.semantic_root());

        let axis_descriptor = axis(
            1,
            "duration-compatible",
            1,
            LogicalUnitV2::Nanoseconds,
            &[scale(LogicalUnitV2::Nanoseconds, 1, 1)],
        );
        let registry = BaseExtensionRegistryProjectionV2::try_new(
            &RunnerLimitsV2::base(RunProfileV2::Smoke),
            &[],
            &[],
            &[axis_descriptor],
        )
        .expect("axis-only registry");
        let axis_mutated = registry
            .try_extent(registered_axis(1), 1, LogicalUnitV2::Nanoseconds)
            .expect("axis mutant");
        assert_ne!(axis_mutated.semantic_root(), base.semantic_root());
    }

    #[test]
    fn registered_descriptors_enforce_names_ids_and_canonical_allowed_units() {
        let claim = no_claim(1);
        assert_eq!(
            RegisteredArtifactRoleDescriptorV2::new(
                ArtifactRoleV2::Observation,
                token("org.example.role.zero"),
                token("org.example.owner"),
                claim.clone(),
            )
            .expect_err("fixed role is not a registered descriptor")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            RegisteredLogicalUnitDescriptorV2::new(
                registered_unit(1),
                token("not-namespaced"),
                token("org.example.owner"),
                claim.clone(),
            )
            .expect_err("global name")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        assert_eq!(
            RegisteredLogicalExtentAxisDescriptorV2::new(
                registered_axis(1),
                token("org.example.axis.empty"),
                token("org.example.owner"),
                claim.clone(),
                LogicalUnitV2::Count,
                &[],
            )
            .expect_err("empty allowed set")
            .kind(),
            ConstructionErrorKindV2::Missing
        );
        let too_many_allowed =
            vec![scale(LogicalUnitV2::Count, 1, 1); LOGICAL_AXIS_ALLOWED_UNITS_MAX_V2 + 1];
        let too_many_error = RegisteredLogicalExtentAxisDescriptorV2::new(
            registered_axis(1),
            token("org.example.axis.too-many-units"),
            token("org.example.owner"),
            claim.clone(),
            LogicalUnitV2::Count,
            &too_many_allowed,
        )
        .expect_err("pre-allocation cap is enforced before duplicate scanning");
        assert_eq!(too_many_error.kind(), ConstructionErrorKindV2::TooLarge);
        assert_eq!(
            too_many_error.field(),
            "extension.logical_axis.allowed_units"
        );
        assert_eq!(
            RegisteredLogicalExtentAxisDescriptorV2::new(
                registered_axis(1),
                token("org.example.axis.missing-canonical"),
                token("org.example.owner"),
                claim.clone(),
                LogicalUnitV2::Count,
                &[scale(LogicalUnitV2::Rows, 1, 1)],
            )
            .expect_err("missing canonical row")
            .kind(),
            ConstructionErrorKindV2::Missing
        );
        assert_eq!(
            RegisteredLogicalExtentAxisDescriptorV2::new(
                registered_axis(1),
                token("org.example.axis.duplicate"),
                token("org.example.owner"),
                claim.clone(),
                LogicalUnitV2::Count,
                &[
                    scale(LogicalUnitV2::Count, 1, 1),
                    scale(LogicalUnitV2::Count, 2, 1),
                ],
            )
            .expect_err("duplicate unit")
            .kind(),
            ConstructionErrorKindV2::Duplicate
        );
        assert_eq!(
            RegisteredLogicalExtentAxisDescriptorV2::new(
                registered_axis(1),
                token("org.example.axis.nonidentity-canonical"),
                token("org.example.owner"),
                claim,
                LogicalUnitV2::Count,
                &[scale(LogicalUnitV2::Count, 2, 1)],
            )
            .expect_err("canonical scale must be one")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            LogicalUnitScaleToCanonicalV2::new(
                LogicalUnitV2::Count,
                RationalV2::new(0, 1).unwrap(),
            )
            .expect_err("zero scale")
            .kind(),
            ConstructionErrorKindV2::OutOfRange
        );
        assert_eq!(
            LogicalUnitScaleToCanonicalV2::from_canonical_parts(LogicalUnitV2::Count, 2, 2,)
                .expect_err("reducible presented ratio")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            LogicalUnitScaleToCanonicalV2::from_canonical_parts(LogicalUnitV2::Count, 1, 0,)
                .expect_err("zero denominator")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            LogicalUnitScaleToCanonicalV2::from_canonical_parts(LogicalUnitV2::Count, -1, 1,)
                .expect_err("negative scale")
                .kind(),
            ConstructionErrorKindV2::OutOfRange
        );
    }

    #[test]
    fn registry_is_typed_permutation_invariant_and_exact_set_reconstructible() {
        let registered_unit_value = LogicalUnitV2::from_tag(16, Some(7)).unwrap();
        let roles = [role(9, "nine", 1), role(7, "seven", 2)];
        let units = [unit(9, "nine", 3), unit(7, "seven", 4)];
        let axes = [
            axis(
                9,
                "nine",
                5,
                LogicalUnitV2::Count,
                &[scale(LogicalUnitV2::Count, 1, 1)],
            ),
            axis(
                7,
                "seven",
                6,
                registered_unit_value,
                &[scale(registered_unit_value, 1, 1)],
            ),
        ];
        let limits = RunnerLimitsV2::base(RunProfileV2::Smoke);
        let first =
            BaseExtensionRegistryProjectionV2::try_new(&limits, &roles, &units, &axes).unwrap();
        let second = BaseExtensionRegistryProjectionV2::try_new(
            &limits,
            &[roles[1].clone(), roles[0].clone()],
            &[units[1].clone(), units[0].clone()],
            &[axes[1].clone(), axes[0].clone()],
        )
        .unwrap();
        assert_eq!(first, second);
        assert_eq!(
            first
                .artifact_roles()
                .iter()
                .map(RegisteredArtifactRoleDescriptorV2::id)
                .collect::<Vec<_>>(),
            [7, 9]
        );
        assert_eq!(
            first.artifact_role(7).unwrap().role().registered_id(),
            Some(7)
        );
        assert_eq!(
            first.logical_unit(7).unwrap().unit().registered_id(),
            Some(7)
        );
        assert_eq!(
            first.logical_axis(7).unwrap().axis().registered_id(),
            Some(7)
        );
        first
            .reconstruct_exact(&limits, &roles, &units, &axes)
            .expect("permuted exact set");

        let roles_extra = [roles[0].clone(), roles[1].clone(), role(11, "eleven", 7)];
        let units_extra = [units[0].clone(), units[1].clone(), unit(11, "eleven", 8)];
        let axes_extra = [
            axes[0].clone(),
            axes[1].clone(),
            axis(
                11,
                "eleven",
                9,
                LogicalUnitV2::Count,
                &[scale(LogicalUnitV2::Count, 1, 1)],
            ),
        ];
        for (result, expected_kind, expected_field) in [
            (
                first.reconstruct_exact(&limits, &roles[..1], &units, &axes),
                ConstructionErrorKindV2::Missing,
                "extension.registry.artifact_roles",
            ),
            (
                first.reconstruct_exact(&limits, &roles_extra, &units, &axes),
                ConstructionErrorKindV2::Unexpected,
                "extension.registry.artifact_roles",
            ),
            (
                first.reconstruct_exact(&limits, &roles, &units[..1], &axes),
                ConstructionErrorKindV2::Missing,
                "extension.registry.logical_units",
            ),
            (
                first.reconstruct_exact(&limits, &roles, &units_extra, &axes),
                ConstructionErrorKindV2::Unexpected,
                "extension.registry.logical_units",
            ),
            (
                first.reconstruct_exact(&limits, &roles, &units, &axes[..1]),
                ConstructionErrorKindV2::Missing,
                "extension.registry.logical_axes",
            ),
            (
                first.reconstruct_exact(&limits, &roles, &units, &axes_extra),
                ConstructionErrorKindV2::Unexpected,
                "extension.registry.logical_axes",
            ),
        ] {
            let error = result.expect_err("category-specific exact-set mismatch");
            assert_eq!(error.kind(), expected_kind);
            assert_eq!(error.field(), expected_field);
        }

        let mutated_roles = [roles[0].clone(), role(8, "eight", 10)];
        let error = first
            .reconstruct_exact(&limits, &mutated_roles, &units, &axes)
            .expect_err("same-count registered-role substitution");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
        assert_eq!(error.field(), "extension.registry.artifact_roles");

        let mutated_units = [units[0].clone(), unit(8, "eight", 10)];
        let error = first
            .reconstruct_exact(&limits, &roles, &mutated_units, &axes)
            .expect_err("same-count registered-unit substitution");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
        assert_eq!(error.field(), "extension.registry.logical_units");

        let substituted_unit_value = LogicalUnitV2::from_tag(16, Some(8)).unwrap();
        let mutated_axes = [
            axes[0].clone(),
            axis(
                7,
                "seven",
                6,
                substituted_unit_value,
                &[scale(substituted_unit_value, 1, 1)],
            ),
        ];
        let error = first
            .reconstruct_exact(&limits, &roles, &units, &mutated_axes)
            .expect_err("same-count unresolved axis-unit substitution");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
        assert_eq!(error.field(), "extension.registry.logical_axes");

        let error = first
            .reconstruct_exact(
                &RunnerLimitsV2::base(RunProfileV2::Full),
                &roles,
                &units,
                &axes,
            )
            .expect_err("different limits root");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
        assert_eq!(error.field(), "extension.registry.limits_root");
    }

    #[test]
    fn registry_refuses_duplicate_collision_unknown_and_over_cap_data() {
        let limits = RunnerLimitsV2::base(RunProfileV2::Smoke);
        let duplicate = [role(1, "first", 1), role(1, "second", 2)];
        assert_eq!(
            BaseExtensionRegistryProjectionV2::try_new(&limits, &duplicate, &[], &[])
                .expect_err("duplicate role ID")
                .kind(),
            ConstructionErrorKindV2::Duplicate
        );

        let same_name_role = role(1, "collision", 1);
        let same_name_unit = RegisteredLogicalUnitDescriptorV2::new(
            registered_unit(1),
            same_name_role.name().clone(),
            token("org.example.owner"),
            no_claim(2),
        )
        .unwrap();
        assert_eq!(
            BaseExtensionRegistryProjectionV2::try_new(
                &limits,
                &[same_name_role],
                &[same_name_unit],
                &[],
            )
            .expect_err("global name collision")
            .kind(),
            ConstructionErrorKindV2::Duplicate
        );

        let unknown_registered_unit = LogicalUnitV2::from_tag(16, Some(77)).unwrap();
        let unknown_axis = axis(
            1,
            "unknown-unit",
            3,
            unknown_registered_unit,
            &[scale(unknown_registered_unit, 1, 1)],
        );
        assert_eq!(
            BaseExtensionRegistryProjectionV2::try_new(&limits, &[], &[], &[unknown_axis])
                .expect_err("unknown unit")
                .kind(),
            ConstructionErrorKindV2::UnknownCode
        );

        let exact = (1_u16..=64)
            .map(|id| role(id, &format!("r{id}"), (id % 251) as u8))
            .collect::<Vec<_>>();
        BaseExtensionRegistryProjectionV2::try_new(&limits, &exact, &[], &[])
            .expect("exact 64 role cap");
        let over = (1_u16..=65)
            .map(|id| role(id, &format!("over{id}"), (id % 251) as u8))
            .collect::<Vec<_>>();
        assert_eq!(
            BaseExtensionRegistryProjectionV2::try_new(&limits, &over, &[], &[])
                .expect_err("65 exceeds role cap")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
    }

    #[test]
    fn registry_category_caps_are_independent_at_zero_one_64_and_65() {
        let limits = RunnerLimitsV2::base(RunProfileV2::Smoke);
        let empty =
            BaseExtensionRegistryProjectionV2::try_new(&limits, &[], &[], &[]).expect("zero rows");
        assert!(empty.artifact_roles().is_empty());
        assert!(empty.logical_units().is_empty());
        assert!(empty.logical_axes().is_empty());

        let one_role = role(1, "one", 1);
        let one_unit = unit(1, "one", 2);
        let one_axis = axis(
            1,
            "one",
            3,
            LogicalUnitV2::Count,
            &[scale(LogicalUnitV2::Count, 1, 1)],
        );
        let one = BaseExtensionRegistryProjectionV2::try_new(
            &limits,
            &[one_role],
            &[one_unit],
            &[one_axis],
        )
        .expect("one row in every typed category");
        assert_eq!(one.artifact_roles().len(), 1);
        assert_eq!(one.logical_units().len(), 1);
        assert_eq!(one.logical_axes().len(), 1);

        let roles_64 = (1_u16..=64)
            .map(|id| role(id, &format!("cap{id}"), (id % 251) as u8))
            .collect::<Vec<_>>();
        let units_64 = (1_u16..=64)
            .map(|id| unit(id, &format!("cap{id}"), ((id + 64) % 251) as u8))
            .collect::<Vec<_>>();
        let axes_64 = (1_u16..=64)
            .map(|id| {
                axis(
                    id,
                    &format!("cap{id}"),
                    ((id + 128) % 251) as u8,
                    LogicalUnitV2::Count,
                    &[scale(LogicalUnitV2::Count, 1, 1)],
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            BaseExtensionRegistryProjectionV2::try_new(&limits, &roles_64, &[], &[])
                .expect("exact role cap")
                .artifact_roles()
                .len(),
            64
        );
        assert_eq!(
            BaseExtensionRegistryProjectionV2::try_new(&limits, &[], &units_64, &[])
                .expect("exact unit cap")
                .logical_units()
                .len(),
            64
        );
        assert_eq!(
            BaseExtensionRegistryProjectionV2::try_new(&limits, &[], &[], &axes_64)
                .expect("exact axis cap")
                .logical_axes()
                .len(),
            64
        );

        let roles_65 = (1_u16..=65)
            .map(|id| role(id, &format!("over-cap{id}"), (id % 251) as u8))
            .collect::<Vec<_>>();
        let units_65 = (1_u16..=65)
            .map(|id| unit(id, &format!("over-cap{id}"), ((id + 64) % 251) as u8))
            .collect::<Vec<_>>();
        let axes_65 = (1_u16..=65)
            .map(|id| {
                axis(
                    id,
                    &format!("over-cap{id}"),
                    ((id + 128) % 251) as u8,
                    LogicalUnitV2::Count,
                    &[scale(LogicalUnitV2::Count, 1, 1)],
                )
            })
            .collect::<Vec<_>>();
        for result in [
            BaseExtensionRegistryProjectionV2::try_new(&limits, &roles_65, &[], &[]),
            BaseExtensionRegistryProjectionV2::try_new(&limits, &[], &units_65, &[]),
            BaseExtensionRegistryProjectionV2::try_new(&limits, &[], &[], &axes_65),
        ] {
            assert_eq!(
                result.expect_err("65 exceeds its own category cap").kind(),
                ConstructionErrorKindV2::TooLarge
            );
        }
    }

    #[test]
    fn u16_max_is_unknown_unless_registered_in_the_exact_typed_namespace() {
        let limits = RunnerLimitsV2::base(RunProfileV2::Smoke);
        let empty = BaseExtensionRegistryProjectionV2::try_new(&limits, &[], &[], &[]).unwrap();
        assert_eq!(
            empty
                .logical_unit(u16::MAX)
                .expect_err("unregistered maximum")
                .kind(),
            ConstructionErrorKindV2::UnknownCode
        );

        let maximum = unit(u16::MAX, "maximum", 1);
        let registered =
            BaseExtensionRegistryProjectionV2::try_new(&limits, &[], &[maximum], &[]).unwrap();
        assert_eq!(registered.logical_unit(u16::MAX).unwrap().id(), u16::MAX);
        assert_eq!(
            registered
                .artifact_role(u16::MAX)
                .expect_err("same numeric ID in another typed namespace remains unknown")
                .kind(),
            ConstructionErrorKindV2::UnknownCode
        );
    }

    #[test]
    fn duration_and_registered_axis_conversions_are_exact_and_bounded() {
        let limits = RunnerLimitsV2::base(RunProfileV2::Smoke);
        let registered_unit_value = LogicalUnitV2::from_tag(16, Some(7)).unwrap();
        let units = [unit(7, "sevens", 1)];
        let axes = [axis(
            7,
            "work",
            2,
            LogicalUnitV2::Operations,
            &[
                scale(LogicalUnitV2::Operations, 1, 1),
                scale(registered_unit_value, 3, 2),
            ],
        )];
        let registry =
            BaseExtensionRegistryProjectionV2::try_new(&limits, &[], &units, &axes).unwrap();

        let seconds = registry
            .try_extent(LogicalExtentAxisV2::Duration, 2, LogicalUnitV2::Seconds)
            .unwrap();
        let nanoseconds = registry
            .convert_extent(seconds, LogicalUnitV2::Nanoseconds)
            .unwrap();
        assert_eq!(nanoseconds.value(), 2_000_000_000);
        assert_eq!(
            registry
                .convert_extent(nanoseconds, LogicalUnitV2::Seconds)
                .unwrap(),
            seconds
        );

        let registered_axis = LogicalExtentAxisV2::from_tag(10, Some(7)).unwrap();
        let alternate = registry
            .try_extent(registered_axis, 2, registered_unit_value)
            .unwrap();
        let canonical = registry
            .convert_extent(alternate, LogicalUnitV2::Operations)
            .unwrap();
        assert_eq!(canonical.value(), 3);
        assert_eq!(
            registry
                .convert_extent(canonical, registered_unit_value)
                .expect("three operations is exactly two alternate units"),
            alternate
        );
        assert_eq!(
            registry
                .convert_extent(
                    registry
                        .try_extent(registered_axis, 1, LogicalUnitV2::Operations)
                        .unwrap(),
                    registered_unit_value,
                )
                .expect_err("one operation is a fractional alternate unit")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        assert_eq!(
            registry
                .try_extent(registered_axis, 1, LogicalUnitV2::Rows)
                .expect_err("unavailable unit")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            registry
                .convert_extent(
                    registry
                        .try_extent(
                            LogicalExtentAxisV2::Duration,
                            u128::MAX,
                            LogicalUnitV2::Seconds,
                        )
                        .unwrap(),
                    LogicalUnitV2::Nanoseconds,
                )
                .expect_err("conversion overflow")
                .kind(),
            ConstructionErrorKindV2::ArithmeticOverflow
        );
    }

    #[test]
    fn unit_conversion_checks_dimensions_normalization_extrema_and_overflow() {
        let second = UnitV2::from_parts(1, 1, [0, 0, 1, 0, 0, 0, 0]).unwrap();
        let nanosecond = UnitV2::from_parts(1, 1_000_000_000, [0, 0, 1, 0, 0, 0, 0]).unwrap();
        assert_eq!(
            normalized_unit_scale_ratio_v2(second, nanosecond).unwrap(),
            RationalV2::from_canonical_parts(1_000_000_000, 1).unwrap()
        );
        assert_eq!(
            convert_rational_quantity_v2(
                RationalV2::from_canonical_parts(2, 1).unwrap(),
                second,
                nanosecond,
            )
            .unwrap(),
            RationalV2::from_canonical_parts(2_000_000_000, 1).unwrap()
        );
        let length = UnitV2::from_parts(1, 1, [1, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(
            normalized_unit_scale_ratio_v2(second, length)
                .expect_err("cross-dimension conversion")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let extreme_scale = UnitV2::from_parts(i128::MAX, 1, [0, 0, 1, 0, 0, 0, 0]).unwrap();
        let value = RationalV2::from_canonical_parts(i128::MAX, 1).unwrap();
        assert_eq!(
            convert_rational_quantity_v2(value, extreme_scale, second)
                .expect_err("i128 product overflow")
                .kind(),
            ConstructionErrorKindV2::ArithmeticOverflow
        );
        assert_eq!(
            convert_rational_quantity_v2(
                RationalV2::from_canonical_parts(i128::MIN, 1).unwrap(),
                second,
                second,
            )
            .unwrap(),
            RationalV2::from_canonical_parts(i128::MIN, 1).unwrap()
        );
    }

    #[test]
    fn every_descriptor_field_and_allowed_unit_mutation_moves_the_registry_root() {
        let limits = RunnerLimitsV2::base(RunProfileV2::Smoke);
        let registered_unit_value = LogicalUnitV2::from_tag(16, Some(7)).unwrap();
        let base_role = role(7, "role", 1);
        let base_unit = unit(7, "unit", 2);
        let base_axis = axis(
            7,
            "axis",
            3,
            LogicalUnitV2::Operations,
            &[
                scale(LogicalUnitV2::Operations, 1, 1),
                scale(registered_unit_value, 3, 2),
            ],
        );
        let root = |role: RegisteredArtifactRoleDescriptorV2,
                    unit: RegisteredLogicalUnitDescriptorV2,
                    axis: RegisteredLogicalExtentAxisDescriptorV2| {
            BaseExtensionRegistryProjectionV2::try_new(&limits, &[role], &[unit], &[axis])
                .map(|projection| *projection.root())
        };
        let expected = root(base_role.clone(), base_unit.clone(), base_axis.clone()).unwrap();

        let role_mutants = [
            role(8, "role", 1),
            role(7, "role-name", 1),
            RegisteredArtifactRoleDescriptorV2::new(
                registered_role(7),
                base_role.name().clone(),
                token("org.example.other-owner"),
                no_claim(1),
            )
            .unwrap(),
            role(7, "role", 9),
        ];
        for mutant in role_mutants {
            match root(mutant, base_unit.clone(), base_axis.clone()) {
                Ok(mutant_root) => assert_ne!(mutant_root, expected),
                Err(error) => assert_ne!(error.kind(), ConstructionErrorKindV2::Unsupported),
            }
        }

        let unit_mutants = [
            unit(8, "unit", 2),
            unit(7, "unit-name", 2),
            RegisteredLogicalUnitDescriptorV2::new(
                registered_unit(7),
                base_unit.name().clone(),
                token("org.example.other-owner"),
                no_claim(2),
            )
            .unwrap(),
            unit(7, "unit", 9),
        ];
        for mutant in unit_mutants {
            match root(base_role.clone(), mutant, base_axis.clone()) {
                Ok(mutant_root) => assert_ne!(mutant_root, expected),
                Err(error) => assert_ne!(error.kind(), ConstructionErrorKindV2::Unsupported),
            }
        }

        let axis_mutants = [
            axis(
                8,
                "axis",
                3,
                LogicalUnitV2::Operations,
                &[scale(LogicalUnitV2::Operations, 1, 1)],
            ),
            axis(
                7,
                "axis-name",
                3,
                LogicalUnitV2::Operations,
                &[scale(LogicalUnitV2::Operations, 1, 1)],
            ),
            RegisteredLogicalExtentAxisDescriptorV2::new(
                registered_axis(7),
                base_axis.name().clone(),
                token("org.example.other-owner"),
                no_claim(3),
                LogicalUnitV2::Operations,
                &[scale(LogicalUnitV2::Operations, 1, 1)],
            )
            .unwrap(),
            axis(
                7,
                "axis",
                9,
                LogicalUnitV2::Operations,
                &[scale(LogicalUnitV2::Operations, 1, 1)],
            ),
            axis(
                7,
                "axis",
                3,
                registered_unit_value,
                &[scale(registered_unit_value, 1, 1)],
            ),
            axis(
                7,
                "axis",
                3,
                LogicalUnitV2::Operations,
                &[
                    scale(LogicalUnitV2::Operations, 1, 1),
                    scale(registered_unit_value, 5, 2),
                ],
            ),
        ];
        for mutant in axis_mutants {
            match root(base_role.clone(), base_unit.clone(), mutant) {
                Ok(mutant_root) => assert_ne!(mutant_root, expected),
                Err(error) => assert_ne!(error.kind(), ConstructionErrorKindV2::Unsupported),
            }
        }
    }
}
