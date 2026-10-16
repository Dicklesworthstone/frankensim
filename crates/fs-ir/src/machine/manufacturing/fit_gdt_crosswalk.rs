//! Explicit applicability links between admitted fit endpoints and geometric
//! tolerance controls.
//!
//! This module closes a semantic traceability gap without inventing geometry:
//! every internal/external endpoint of an admitted fit requirement names one
//! admitted geometric-tolerance control on the same caller-declared body. The
//! receipt retains the contact feature and controlled surface patch as distinct
//! identities. It does not prove that either selector is physically contained
//! by that body, that they describe the same geometric feature, or that the
//! selected control is sufficient for fit, clearance, inspection, or assembly.

use core::fmt;

use std::collections::BTreeMap;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, NeverCancel, ProblemSemanticId, StrongIdentity, WireType,
};

use crate::IR_VERSION;

use super::super::{BodyId, ContactFeatureId, MachineGraphIdV1, SurfacePatchId};
use super::datum_system::DatumReferenceFrameIdV1;
use super::fit_clearance::{
    AdmittedMachineFitClearanceV1, FitRequirementIdV1, MachineFitClearanceIdV1,
};
use super::geometric_tolerance::{
    AdmittedMachineGeometricToleranceV1, GeometricCharacteristicV1, GeometricToleranceControlIdV1,
    GeometricToleranceLengthV1, MachineGeometricToleranceIdV1,
};

/// Identity/admission schema for the fit/GD&T applicability crosswalk.
pub const MACHINE_FIT_GDT_CROSSWALK_SCHEMA_VERSION_V1: u32 = 1;
/// Maximum endpoint links retained by version one.
pub const MAX_MACHINE_FIT_GDT_LINKS_V1: usize = 8_192;

const FIT_GDT_CROSSWALK_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(12 * 1_024 * 1_024, 8 * 1_024 * 1_024, 6, 8_192, 4_096);

/// Role of one fit endpoint in the role-complete crosswalk.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum FitGdtEndpointRoleV1 {
    /// Internal fit feature.
    Internal = 1,
    /// External fit feature.
    External = 2,
}

impl FitGdtEndpointRoleV1 {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }
}

/// Caller-declared applicability link from one fit endpoint to one control.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FitGdtEndpointLinkV1 {
    fit_requirement: FitRequirementIdV1,
    role: FitGdtEndpointRoleV1,
    geometric_control: GeometricToleranceControlIdV1,
}

impl FitGdtEndpointLinkV1 {
    /// Construct one authority-free endpoint/control association.
    #[must_use]
    pub const fn new(
        fit_requirement: FitRequirementIdV1,
        role: FitGdtEndpointRoleV1,
        geometric_control: GeometricToleranceControlIdV1,
    ) -> Self {
        Self {
            fit_requirement,
            role,
            geometric_control,
        }
    }

    /// Referenced fit requirement.
    #[must_use]
    pub const fn fit_requirement(&self) -> &FitRequirementIdV1 {
        &self.fit_requirement
    }

    /// Internal or external endpoint role.
    #[must_use]
    pub const fn role(&self) -> FitGdtEndpointRoleV1 {
        self.role
    }

    /// Referenced geometric-tolerance control.
    #[must_use]
    pub const fn geometric_control(&self) -> &GeometricToleranceControlIdV1 {
        &self.geometric_control
    }
}

/// One resolved applicability row retained by the admitted crosswalk.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FitGdtEndpointReceiptV1 {
    link: FitGdtEndpointLinkV1,
    declared_body: BodyId,
    fit_feature: ContactFeatureId,
    controlled_patch: SurfacePatchId,
    characteristic: GeometricCharacteristicV1,
    zone_width: GeometricToleranceLengthV1,
    datum_frame: Option<DatumReferenceFrameIdV1>,
}

impl FitGdtEndpointReceiptV1 {
    /// Exact submitted applicability tuple.
    #[must_use]
    pub const fn link(&self) -> &FitGdtEndpointLinkV1 {
        &self.link
    }

    /// Caller-declared body shared by the two upstream selectors.
    #[must_use]
    pub const fn declared_body(&self) -> &BodyId {
        &self.declared_body
    }

    /// Durable contact feature selected by the fit endpoint.
    #[must_use]
    pub const fn fit_feature(&self) -> &ContactFeatureId {
        &self.fit_feature
    }

    /// Durable surface patch selected by the geometric control.
    #[must_use]
    pub const fn controlled_patch(&self) -> &SurfacePatchId {
        &self.controlled_patch
    }

    /// Closed upstream geometric characteristic.
    #[must_use]
    pub const fn characteristic(&self) -> GeometricCharacteristicV1 {
        self.characteristic
    }

    /// Upstream positive zone width, retained without clearance arithmetic.
    #[must_use]
    pub const fn zone_width(&self) -> GeometricToleranceLengthV1 {
        self.zone_width
    }

    /// Upstream datum-frame reference, when the control is an orientation control.
    #[must_use]
    pub const fn datum_frame(&self) -> Option<&DatumReferenceFrameIdV1> {
        self.datum_frame.as_ref()
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(768);
        append_bytes(
            &mut row,
            self.link.fit_requirement.canonical_key().as_bytes(),
        );
        row.push(self.link.role.tag());
        append_bytes(&mut row, self.declared_body.identity().as_bytes());
        append_bytes(&mut row, self.fit_feature.identity().as_bytes());
        append_bytes(
            &mut row,
            self.link.geometric_control.canonical_key().as_bytes(),
        );
        append_bytes(&mut row, self.controlled_patch.identity().as_bytes());
        row.push(self.characteristic.tag());
        row.extend_from_slice(&self.zone_width.submitted_value().to_bits().to_le_bytes());
        row.push(self.zone_width.unit().tag());
        row.extend_from_slice(&self.zone_width.metres().to_bits().to_le_bytes());
        match &self.datum_frame {
            Some(frame) => {
                row.push(1);
                append_bytes(&mut row, frame.canonical_key().as_bytes());
            }
            None => row.push(0),
        }
        row
    }
}

/// Mutable-by-construction endpoint applicability table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineFitGdtCrosswalkDraftV1 {
    /// Endpoint/control links in non-semantic caller order.
    pub links: Vec<FitGdtEndpointLinkV1>,
}

impl MachineFitGdtCrosswalkDraftV1 {
    /// Admit a role-complete crosswalk against two already-admitted catalogs.
    ///
    /// # Errors
    /// Refuses graph mismatch, empty/oversized or aliased links, unknown IDs,
    /// body mismatch, incomplete endpoint coverage, or identity failure.
    #[allow(clippy::result_large_err)] // Preserve exact owned IDs in refusals.
    pub fn admit_against(
        self,
        fit: &AdmittedMachineFitClearanceV1,
        geometric_tolerance: &AdmittedMachineGeometricToleranceV1,
    ) -> Result<AdmittedMachineFitGdtCrosswalkV1, MachineFitGdtCrosswalkAdmissionErrorV1> {
        admit_fit_gdt_crosswalk(self, fit, geometric_tolerance)
    }
}

/// Canonical identity schema for one fit/GD&T applicability crosswalk.
pub enum MachineFitGdtCrosswalkIdentitySchemaV1 {}

impl CanonicalSchema for MachineFitGdtCrosswalkIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.fit-gdt-crosswalk.v1";
    const NAME: &'static str = "admitted-machine-fit-gdt-crosswalk";
    const VERSION: u32 = MACHINE_FIT_GDT_CROSSWALK_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "one exact fit catalog and geometric-tolerance catalog plus role-complete applicability links";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("fit-gdt-crosswalk-schema-version", WireType::U64),
        FieldSpec::required("frankenscript-ir-version", WireType::U64),
        FieldSpec::required("machine-graph", WireType::Bytes),
        FieldSpec::required("fit-catalog", WireType::Bytes),
        FieldSpec::required("geometric-tolerance-catalog", WireType::Bytes),
        FieldSpec::required("endpoint-links", WireType::OrderedBytes),
    ];
}

/// Strong semantic identity of one admitted fit/GD&T crosswalk.
pub type MachineFitGdtCrosswalkIdV1 = ProblemSemanticId<MachineFitGdtCrosswalkIdentitySchemaV1>;

/// Canonically ordered graph/catalog-bound applicability receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedMachineFitGdtCrosswalkV1 {
    graph: MachineGraphIdV1,
    fit_catalog: IdentityReceipt<MachineFitClearanceIdV1>,
    geometric_tolerance_catalog: IdentityReceipt<MachineGeometricToleranceIdV1>,
    endpoints: Vec<FitGdtEndpointReceiptV1>,
    receipt: IdentityReceipt<MachineFitGdtCrosswalkIdV1>,
}

impl AdmittedMachineFitGdtCrosswalkV1 {
    /// Exact Machine graph shared by both input catalogs.
    #[must_use]
    pub const fn graph(&self) -> MachineGraphIdV1 {
        self.graph
    }

    /// Exact admitted fit-catalog identity.
    #[must_use]
    pub const fn fit_catalog(&self) -> MachineFitClearanceIdV1 {
        self.fit_catalog.id()
    }

    /// Complete upstream fit-catalog canonical-preimage receipt.
    #[must_use]
    pub const fn fit_catalog_receipt(&self) -> IdentityReceipt<MachineFitClearanceIdV1> {
        self.fit_catalog
    }

    /// Exact admitted geometric-tolerance-catalog identity.
    #[must_use]
    pub const fn geometric_tolerance_catalog(&self) -> MachineGeometricToleranceIdV1 {
        self.geometric_tolerance_catalog.id()
    }

    /// Complete upstream geometric-tolerance canonical-preimage receipt.
    #[must_use]
    pub const fn geometric_tolerance_catalog_receipt(
        &self,
    ) -> IdentityReceipt<MachineGeometricToleranceIdV1> {
        self.geometric_tolerance_catalog
    }

    /// Resolved links in fit-requirement then endpoint-role order.
    #[must_use]
    pub fn endpoints(&self) -> &[FitGdtEndpointReceiptV1] {
        &self.endpoints
    }

    /// Domain-separated aggregate identity.
    #[must_use]
    pub const fn identity(&self) -> MachineFitGdtCrosswalkIdV1 {
        self.receipt.id()
    }

    /// Complete canonical-preimage receipt for collision adjudication.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<MachineFitGdtCrosswalkIdV1> {
        self.receipt
    }
}

/// Structured refusal from fit/GD&T crosswalk admission.
#[allow(clippy::large_enum_variant)] // Preserve exact owned IDs in rich refusals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MachineFitGdtCrosswalkAdmissionErrorV1 {
    /// At least one endpoint link is required.
    NoLinks,
    /// Raw submissions exceeded the fixed endpoint-link cap.
    LinkLimit {
        /// Submitted link count.
        actual: usize,
        /// Maximum admitted count.
        max: usize,
    },
    /// The admitted catalogs bind different Machine graphs.
    CatalogGraphMismatch {
        /// Graph bound by the fit catalog.
        fit_graph: MachineGraphIdV1,
        /// Graph bound by the geometric-tolerance catalog.
        geometric_tolerance_graph: MachineGraphIdV1,
    },
    /// More than one link targeted the same fit endpoint.
    DuplicateEndpoint {
        /// Fit requirement containing the duplicated role.
        requirement: FitRequirementIdV1,
        /// Duplicated endpoint role.
        role: FitGdtEndpointRoleV1,
        /// Control named by the first canonical row.
        first_control: GeometricToleranceControlIdV1,
        /// Control named by the duplicate row.
        duplicate_control: GeometricToleranceControlIdV1,
    },
    /// A link named no admitted fit requirement.
    UnknownFitRequirement {
        /// Missing fit requirement.
        requirement: FitRequirementIdV1,
        /// Endpoint role requested from it.
        role: FitGdtEndpointRoleV1,
    },
    /// A link named no admitted geometric-tolerance control.
    UnknownGeometricControl {
        /// Fit requirement containing the link.
        requirement: FitRequirementIdV1,
        /// Endpoint role containing the link.
        role: FitGdtEndpointRoleV1,
        /// Missing geometric-tolerance control.
        control: GeometricToleranceControlIdV1,
    },
    /// The selected fit endpoint and control declare different bodies.
    DeclaredBodyMismatch {
        /// Fit requirement containing the link.
        requirement: FitRequirementIdV1,
        /// Endpoint role containing the link.
        role: FitGdtEndpointRoleV1,
        /// Selected geometric-tolerance control.
        control: GeometricToleranceControlIdV1,
        /// Body declared by the fit endpoint.
        fit_body: BodyId,
        /// Body declared by the geometric-tolerance control.
        geometric_tolerance_body: BodyId,
    },
    /// One admitted fit endpoint had no applicability link.
    MissingEndpoint {
        /// Fit requirement missing coverage.
        requirement: FitRequirementIdV1,
        /// Missing endpoint role.
        role: FitGdtEndpointRoleV1,
    },
    /// Canonical aggregate identity publication failed.
    Identity(CanonicalError),
}

impl MachineFitGdtCrosswalkAdmissionErrorV1 {
    /// Stable machine-actionable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::NoLinks => "MachineFitGdtNoLinks",
            Self::LinkLimit { .. } => "MachineFitGdtLinkLimit",
            Self::CatalogGraphMismatch { .. } => "MachineFitGdtCatalogGraphMismatch",
            Self::DuplicateEndpoint { .. } => "MachineFitGdtDuplicateEndpoint",
            Self::UnknownFitRequirement { .. } => "MachineFitGdtUnknownFitRequirement",
            Self::UnknownGeometricControl { .. } => "MachineFitGdtUnknownGeometricControl",
            Self::DeclaredBodyMismatch { .. } => "MachineFitGdtDeclaredBodyMismatch",
            Self::MissingEndpoint { .. } => "MachineFitGdtMissingEndpoint",
            Self::Identity(_) => "MachineFitGdtIdentity",
        }
    }
}

impl fmt::Display for MachineFitGdtCrosswalkAdmissionErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoLinks => formatter.write_str("fit/GD&T crosswalk requires endpoint links"),
            Self::LinkLimit { actual, max } => write!(
                formatter,
                "fit/GD&T crosswalk has {actual} links; maximum is {max}"
            ),
            Self::CatalogGraphMismatch {
                fit_graph,
                geometric_tolerance_graph,
            } => write!(
                formatter,
                "fit catalog graph {fit_graph} differs from geometric-tolerance graph \
                 {geometric_tolerance_graph}"
            ),
            Self::DuplicateEndpoint {
                requirement,
                role,
                first_control,
                duplicate_control,
            } => write!(
                formatter,
                "fit requirement {requirement} {role:?} endpoint links both control \
                 {first_control} and {duplicate_control}"
            ),
            Self::UnknownFitRequirement { requirement, role } => write!(
                formatter,
                "fit/GD&T link names unknown fit requirement {requirement} at role {role:?}"
            ),
            Self::UnknownGeometricControl {
                requirement,
                role,
                control,
            } => write!(
                formatter,
                "fit requirement {requirement} {role:?} link names unknown control {control}"
            ),
            Self::DeclaredBodyMismatch {
                requirement,
                role,
                control,
                fit_body,
                geometric_tolerance_body,
            } => write!(
                formatter,
                "fit requirement {requirement} {role:?} declares body {fit_body}, but control \
                 {control} declares body {geometric_tolerance_body}"
            ),
            Self::MissingEndpoint { requirement, role } => write!(
                formatter,
                "fit requirement {requirement} is missing its {role:?} applicability link"
            ),
            Self::Identity(error) => write!(formatter, "fit/GD&T identity refused: {error}"),
        }
    }
}

impl std::error::Error for MachineFitGdtCrosswalkAdmissionErrorV1 {}

impl From<CanonicalError> for MachineFitGdtCrosswalkAdmissionErrorV1 {
    fn from(error: CanonicalError) -> Self {
        Self::Identity(error)
    }
}

#[allow(clippy::too_many_lines)] // One atomic cross-catalog resolution and publication pass.
#[allow(clippy::result_large_err)] // Preserve exact owned IDs in structured refusals.
fn admit_fit_gdt_crosswalk(
    draft: MachineFitGdtCrosswalkDraftV1,
    fit: &AdmittedMachineFitClearanceV1,
    geometric_tolerance: &AdmittedMachineGeometricToleranceV1,
) -> Result<AdmittedMachineFitGdtCrosswalkV1, MachineFitGdtCrosswalkAdmissionErrorV1> {
    if draft.links.is_empty() {
        return Err(MachineFitGdtCrosswalkAdmissionErrorV1::NoLinks);
    }
    if draft.links.len() > MAX_MACHINE_FIT_GDT_LINKS_V1 {
        return Err(MachineFitGdtCrosswalkAdmissionErrorV1::LinkLimit {
            actual: draft.links.len(),
            max: MAX_MACHINE_FIT_GDT_LINKS_V1,
        });
    }
    let graph = fit.graph();
    if geometric_tolerance.graph() != graph {
        return Err(
            MachineFitGdtCrosswalkAdmissionErrorV1::CatalogGraphMismatch {
                fit_graph: graph,
                geometric_tolerance_graph: geometric_tolerance.graph(),
            },
        );
    }

    let mut links = draft.links;
    links.sort_by(|left, right| {
        left.fit_requirement
            .cmp(&right.fit_requirement)
            .then_with(|| left.role.cmp(&right.role))
            .then_with(|| left.geometric_control.cmp(&right.geometric_control))
    });
    if let Some(pair) = links.windows(2).find(|pair| {
        pair[0].fit_requirement == pair[1].fit_requirement && pair[0].role == pair[1].role
    }) {
        return Err(MachineFitGdtCrosswalkAdmissionErrorV1::DuplicateEndpoint {
            requirement: pair[0].fit_requirement.clone(),
            role: pair[0].role,
            first_control: pair[0].geometric_control.clone(),
            duplicate_control: pair[1].geometric_control.clone(),
        });
    }

    let fit_requirements = fit
        .requirements()
        .iter()
        .map(|requirement| (requirement.id().clone(), requirement))
        .collect::<BTreeMap<_, _>>();
    let controls = geometric_tolerance
        .controls()
        .iter()
        .map(|control| (control.id().clone(), control))
        .collect::<BTreeMap<_, _>>();
    let mut coverage = fit_requirements
        .keys()
        .cloned()
        .map(|requirement| (requirement, 0_u8))
        .collect::<BTreeMap<_, _>>();
    let mut endpoints = Vec::with_capacity(links.len());

    for link in links {
        let Some(requirement) = fit_requirements.get(&link.fit_requirement) else {
            return Err(
                MachineFitGdtCrosswalkAdmissionErrorV1::UnknownFitRequirement {
                    requirement: link.fit_requirement,
                    role: link.role,
                },
            );
        };
        let selector = match link.role {
            FitGdtEndpointRoleV1::Internal => requirement.target().internal(),
            FitGdtEndpointRoleV1::External => requirement.target().external(),
        };
        let Some(control) = controls.get(&link.geometric_control) else {
            return Err(
                MachineFitGdtCrosswalkAdmissionErrorV1::UnknownGeometricControl {
                    requirement: link.fit_requirement,
                    role: link.role,
                    control: link.geometric_control,
                },
            );
        };
        if selector.declared_body() != control.declared_body() {
            return Err(
                MachineFitGdtCrosswalkAdmissionErrorV1::DeclaredBodyMismatch {
                    requirement: link.fit_requirement,
                    role: link.role,
                    control: link.geometric_control,
                    fit_body: selector.declared_body().clone(),
                    geometric_tolerance_body: control.declared_body().clone(),
                },
            );
        }
        let bit = match link.role {
            FitGdtEndpointRoleV1::Internal => 1,
            FitGdtEndpointRoleV1::External => 2,
        };
        let Some(mask) = coverage.get_mut(&link.fit_requirement) else {
            return Err(
                MachineFitGdtCrosswalkAdmissionErrorV1::UnknownFitRequirement {
                    requirement: link.fit_requirement,
                    role: link.role,
                },
            );
        };
        *mask |= bit;
        endpoints.push(FitGdtEndpointReceiptV1 {
            link,
            declared_body: selector.declared_body().clone(),
            fit_feature: selector.feature().clone(),
            controlled_patch: control.controlled_patch().clone(),
            characteristic: control.characteristic(),
            zone_width: control.zone_width(),
            datum_frame: control.datum_frame().cloned(),
        });
    }

    for (requirement, mask) in coverage {
        if mask & 1 == 0 {
            return Err(MachineFitGdtCrosswalkAdmissionErrorV1::MissingEndpoint {
                requirement,
                role: FitGdtEndpointRoleV1::Internal,
            });
        }
        if mask & 2 == 0 {
            return Err(MachineFitGdtCrosswalkAdmissionErrorV1::MissingEndpoint {
                requirement,
                role: FitGdtEndpointRoleV1::External,
            });
        }
    }

    let rows = endpoints
        .iter()
        .map(FitGdtEndpointReceiptV1::canonical_row)
        .collect::<Vec<_>>();
    let fit_catalog = fit.identity_receipt();
    let geometric_tolerance_catalog = geometric_tolerance.identity_receipt();
    let receipt = CanonicalEncoder::<MachineFitGdtCrosswalkIdV1, _>::new(
        FIT_GDT_CROSSWALK_IDENTITY_LIMITS,
        NeverCancel,
    )?
    .u64(
        Field::new(0, "fit-gdt-crosswalk-schema-version"),
        u64::from(MACHINE_FIT_GDT_CROSSWALK_SCHEMA_VERSION_V1),
    )?
    .u64(
        Field::new(1, "frankenscript-ir-version"),
        u64::from(IR_VERSION),
    )?
    .bytes(Field::new(2, "machine-graph"), graph.as_bytes())?
    .bytes(Field::new(3, "fit-catalog"), fit_catalog.id().as_bytes())?
    .bytes(
        Field::new(4, "geometric-tolerance-catalog"),
        geometric_tolerance_catalog.id().as_bytes(),
    )?
    .ordered_bytes(
        Field::new(5, "endpoint-links"),
        rows.len() as u64,
        rows.iter().map(Vec::as_slice),
    )?
    .finish()?;

    Ok(AdmittedMachineFitGdtCrosswalkV1 {
        graph,
        fit_catalog,
        geometric_tolerance_catalog,
        endpoints,
        receipt,
    })
}

fn append_bytes(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}
