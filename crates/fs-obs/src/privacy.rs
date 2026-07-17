//! Deterministic redaction, retention, and share-policy model (i94v.7.3.4).
//!
//! Values arrive with field-level labels; this module does not guess secrets
//! by scanning bytes. It produces a bounded share manifest whose disclosures
//! cannot contain a protected field unless the exact sensitivity, license,
//! export, audience, and privilege predicates all admit it.

use crate::{EvidenceCompleteness, EvidenceIntegrity, OperationalSupport, PromotionEffect};
use core::fmt;

/// Current share/redaction policy semantics.
pub const SHARE_POLICY_VERSION: u32 = 1;

fn valid_token(value: &str, max_len: usize) -> bool {
    !value.is_empty() && value.len() <= max_len && value.chars().all(|ch| !ch.is_control())
}

/// Field sensitivity. Labels are explicit metadata, never inferred from text.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Sensitivity {
    /// Safe for public disclosure when license/export policy also permits it.
    Public,
    /// Organization-internal information.
    Internal,
    /// Personal or identifying information.
    PersonalData,
    /// Non-credential confidential secret.
    Secret,
    /// Credential, private key, bearer token, or authentication material.
    Credential,
}

impl Sensitivity {
    /// Stable policy name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Public => "public",
            Self::Internal => "internal",
            Self::PersonalData => "personal_data",
            Self::Secret => "secret",
            Self::Credential => "credential",
        }
    }

    fn redaction_marker(self) -> String {
        format!("[REDACTED:{}]", self.name())
    }
}

/// Redistribution/license realm attached to a field or artifact reference.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum LicenseRealm {
    /// Public-domain or explicitly redistributable material.
    Redistributable,
    /// Named entitlement required outside the local process.
    Restricted(String),
    /// Never leaves the local trust boundary.
    LocalOnly,
}

/// Export-control realm attached to a field or artifact reference.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum ExportRealm {
    /// No export restriction represented by this policy.
    Unrestricted,
    /// Named authorization/capability required for disclosure.
    Controlled(String),
    /// Never leaves the local trust boundary.
    LocalOnly,
}

/// Retention requirement. Time values are caller-supplied Unix nanoseconds;
/// this pure module reads no clock.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum RetentionClass {
    /// Discard after the current session unless another policy overrides it.
    Ephemeral,
    /// Retain until the supplied timestamp, then treat as expired.
    Until(u64),
    /// Retain without a policy-level expiry.
    Permanent,
    /// Preserve locally under a named legal hold. A hold never grants sharing.
    LegalHold(String),
}

/// Complete policy label for one field.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FieldPolicy {
    /// Sensitivity label.
    pub sensitivity: Sensitivity,
    /// License/redistribution realm.
    pub license: LicenseRealm,
    /// Export realm.
    pub export: ExportRealm,
    /// Retention class.
    pub retention: RetentionClass,
}

impl FieldPolicy {
    fn validate(&self) -> Result<(), FieldError> {
        match &self.license {
            LicenseRealm::Restricted(realm) if !valid_token(realm, 256) => {
                return Err(FieldError::InvalidRealm {
                    kind: "license",
                    realm: realm.clone(),
                });
            }
            LicenseRealm::Redistributable
            | LicenseRealm::Restricted(_)
            | LicenseRealm::LocalOnly => {}
        }
        match &self.export {
            ExportRealm::Controlled(realm) if !valid_token(realm, 256) => {
                return Err(FieldError::InvalidRealm {
                    kind: "export",
                    realm: realm.clone(),
                });
            }
            ExportRealm::Unrestricted | ExportRealm::Controlled(_) | ExportRealm::LocalOnly => {}
        }
        if let RetentionClass::LegalHold(hold) = &self.retention {
            if !valid_token(hold, 256) {
                return Err(FieldError::InvalidRealm {
                    kind: "legal_hold",
                    realm: hold.clone(),
                });
            }
        }
        Ok(())
    }
}

/// One path-addressed, explicitly labeled value. Bytes are opaque, so Unicode,
/// binary, and hostile encodings cannot bypass the label.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LabeledField {
    path: String,
    value: Vec<u8>,
    policy: FieldPolicy,
    required_for_replay: bool,
}

impl LabeledField {
    /// Validate one labeled field.
    ///
    /// # Errors
    /// [`FieldError`] when its path or a named policy realm is malformed.
    pub fn new(
        path: impl Into<String>,
        value: impl Into<Vec<u8>>,
        policy: FieldPolicy,
        required_for_replay: bool,
    ) -> Result<Self, FieldError> {
        let path = path.into();
        if !valid_token(&path, 512) {
            return Err(FieldError::InvalidPath { path });
        }
        policy.validate()?;
        Ok(Self {
            path,
            value: value.into(),
            policy,
            required_for_replay,
        })
    }

    /// Stable structured path.
    #[must_use]
    pub fn path(&self) -> &str {
        &self.path
    }

    /// Explicit policy label.
    #[must_use]
    pub const fn policy(&self) -> &FieldPolicy {
        &self.policy
    }

    /// Whether omission weakens replay/evidence completeness.
    #[must_use]
    pub const fn required_for_replay(&self) -> bool {
        self.required_for_replay
    }
}

/// Refusal for malformed labeled fields.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FieldError {
    /// Field path is empty, oversized, or control-bearing.
    InvalidPath {
        /// Rejected path.
        path: String,
    },
    /// Named license/export/hold realm is malformed.
    InvalidRealm {
        /// Realm category.
        kind: &'static str,
        /// Rejected realm.
        realm: String,
    },
}

impl fmt::Display for FieldError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPath { path } => write!(
                f,
                "field path must be 1..=512 bytes without control characters; got {path:?}"
            ),
            Self::InvalidRealm { kind, realm } => write!(
                f,
                "{kind} realm must be 1..=256 bytes without control characters; got {realm:?}"
            ),
        }
    }
}

impl core::error::Error for FieldError {}

/// External correlation-token method. Tokens are supplied by an admitted
/// cryptographic owner; this std-only crate does not implement a password hash,
/// KDF, MAC, or secret-bearing salt store.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum ExternalCorrelationMethod {
    /// Salted digest. Not sufficient for PII/secrets vulnerable to dictionary
    /// attack; policy rejects those combinations.
    Salted { salt_id: String },
    /// Keyed opaque token suitable for PII/secret correlation when admitted by
    /// the caller's trust boundary.
    Keyed { key_id: String },
}

/// Already-computed external correlation token.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ExternalCorrelationToken {
    method: ExternalCorrelationMethod,
    token: String,
}

impl ExternalCorrelationToken {
    /// Validate an externally computed token and its non-secret key/salt ID.
    ///
    /// # Errors
    /// [`CorrelationTokenError`] for malformed IDs or non-hex tokens.
    pub fn new(
        method: ExternalCorrelationMethod,
        token: impl Into<String>,
    ) -> Result<Self, CorrelationTokenError> {
        let identity = match &method {
            ExternalCorrelationMethod::Salted { salt_id } => salt_id,
            ExternalCorrelationMethod::Keyed { key_id } => key_id,
        };
        if !valid_token(identity, 256) {
            return Err(CorrelationTokenError::InvalidAuthorityId {
                identity: identity.clone(),
            });
        }
        let token = token.into();
        if !(16..=256).contains(&token.len()) || !token.bytes().all(|byte| byte.is_ascii_hexdigit())
        {
            return Err(CorrelationTokenError::InvalidToken { token });
        }
        Ok(Self { method, token })
    }
}

/// Refusal for a malformed external correlation token.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CorrelationTokenError {
    /// Salt/key identity is malformed.
    InvalidAuthorityId {
        /// Rejected identity.
        identity: String,
    },
    /// Token is not bounded hexadecimal data.
    InvalidToken {
        /// Rejected token.
        token: String,
    },
}

impl fmt::Display for CorrelationTokenError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAuthorityId { identity } => {
                write!(f, "correlation authority ID is invalid: {identity:?}")
            }
            Self::InvalidToken { token } => write!(f, "correlation token is invalid: {token:?}"),
        }
    }
}

impl core::error::Error for CorrelationTokenError {}

/// Correlation disclosure requested for redacted fields.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CorrelationPolicy {
    /// Fixed marker only; safest default.
    None,
    /// Unsalted FNV correlation, allowed only for public/internal data and
    /// explicitly non-cryptographic.
    UnsaltedPublic,
    /// Caller-supplied salted/keyed opaque token.
    External(ExternalCorrelationToken),
}

/// Destination audience.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShareAudience {
    /// Current local trust boundary.
    Local,
    /// Named organization boundary.
    Organization,
    /// Public/untrusted publication.
    Public,
}

/// Validated, bounded share request.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ShareRequest {
    audience: ShareAudience,
    now_ns: u64,
    max_fields: usize,
    max_input_bytes: usize,
    privileged: bool,
    allow_personal_data: bool,
    allow_secrets: bool,
    licensed_realms: Vec<String>,
    export_authorizations: Vec<String>,
    correlation: CorrelationPolicy,
}

impl ShareRequest {
    /// Construct a bounded request with deny-by-default privilege and
    /// correlation settings.
    ///
    /// # Errors
    /// [`ShareRequestError`] when a bound is zero.
    pub fn new(
        audience: ShareAudience,
        now_ns: u64,
        max_fields: usize,
        max_input_bytes: usize,
    ) -> Result<Self, ShareRequestError> {
        if max_fields == 0 {
            return Err(ShareRequestError::ZeroBound {
                field: "max_fields",
            });
        }
        if max_input_bytes == 0 {
            return Err(ShareRequestError::ZeroBound {
                field: "max_input_bytes",
            });
        }
        Ok(Self {
            audience,
            now_ns,
            max_fields,
            max_input_bytes,
            privileged: false,
            allow_personal_data: false,
            allow_secrets: false,
            licensed_realms: Vec::new(),
            export_authorizations: Vec::new(),
            correlation: CorrelationPolicy::None,
        })
    }

    /// Set privilege explicitly. Privilege alone never reveals credentials.
    #[must_use]
    pub fn with_privilege(mut self, privileged: bool) -> Self {
        self.privileged = privileged;
        self
    }

    /// Permit personal data only inside a privileged local request.
    #[must_use]
    pub fn with_personal_data(mut self, allow: bool) -> Self {
        self.allow_personal_data = allow;
        self
    }

    /// Permit non-credential secrets only inside a privileged local request.
    #[must_use]
    pub fn with_secrets(mut self, allow: bool) -> Self {
        self.allow_secrets = allow;
        self
    }

    /// Replace named license entitlements with a sorted, duplicate-free list.
    ///
    /// # Errors
    /// [`ShareRequestError::InvalidRealm`] for malformed realm names.
    pub fn with_license_realms(
        mut self,
        realms: impl IntoIterator<Item = String>,
    ) -> Result<Self, ShareRequestError> {
        self.licensed_realms = normalize_realms("license", realms)?;
        Ok(self)
    }

    /// Replace named export authorizations with a sorted, duplicate-free list.
    ///
    /// # Errors
    /// [`ShareRequestError::InvalidRealm`] for malformed realm names.
    pub fn with_export_authorizations(
        mut self,
        realms: impl IntoIterator<Item = String>,
    ) -> Result<Self, ShareRequestError> {
        self.export_authorizations = normalize_realms("export", realms)?;
        Ok(self)
    }

    /// Select correlation policy for redacted fields.
    #[must_use]
    pub fn with_correlation(mut self, correlation: CorrelationPolicy) -> Self {
        self.correlation = correlation;
        self
    }
}

fn normalize_realms(
    kind: &'static str,
    realms: impl IntoIterator<Item = String>,
) -> Result<Vec<String>, ShareRequestError> {
    let mut realms: Vec<_> = realms.into_iter().collect();
    if let Some(realm) = realms.iter().find(|realm| !valid_token(realm, 256)) {
        return Err(ShareRequestError::InvalidRealm {
            kind,
            realm: realm.clone(),
        });
    }
    realms.sort();
    realms.dedup();
    Ok(realms)
}

/// Refusal for a malformed share request.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ShareRequestError {
    /// A resident/input bound is zero.
    ZeroBound {
        /// Invalid field.
        field: &'static str,
    },
    /// Named entitlement is malformed.
    InvalidRealm {
        /// Entitlement category.
        kind: &'static str,
        /// Rejected value.
        realm: String,
    },
}

impl fmt::Display for ShareRequestError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroBound { field } => write!(f, "share request bound {field} must be non-zero"),
            Self::InvalidRealm { kind, realm } => {
                write!(f, "invalid {kind} realm in share request: {realm:?}")
            }
        }
    }
}

impl core::error::Error for ShareRequestError {}

/// Why a field was not disclosed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShareBlock {
    /// Audience/privilege does not admit the sensitivity label.
    Sensitivity,
    /// Credentials are never emitted, including privileged local views.
    Credential,
    /// License realm does not admit the destination.
    License,
    /// Export realm does not admit the destination.
    Export,
    /// Retention window has expired.
    Expired,
}

/// Correlation token that is safe to include in the selected manifest.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CorrelationDisclosure {
    /// Non-cryptographic public/internal fingerprint.
    UnsaltedFnv1a64(String),
    /// Externally computed salted token.
    ExternalSalted {
        /// Non-secret salt identity.
        salt_id: String,
        /// Opaque hexadecimal token.
        token: String,
    },
    /// Externally computed keyed token.
    ExternalKeyed {
        /// Non-secret key identity.
        key_id: String,
        /// Opaque hexadecimal token.
        token: String,
    },
}

/// Publishable field projection. Protected raw bytes exist only in `Plain`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Disclosure {
    /// Raw bytes were admitted.
    Plain(Vec<u8>),
    /// Stable marker replaced raw bytes.
    Redacted {
        /// Deterministic label-only marker.
        marker: String,
        /// Optional safe correlation token.
        correlation: Option<CorrelationDisclosure>,
        /// Why raw bytes were withheld.
        reason: ShareBlock,
    },
    /// Field is absent from payload but retained as a manifest decision.
    Omitted {
        /// Why no payload representation is safe.
        reason: ShareBlock,
    },
}

/// Retention action retained alongside a share decision.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RetentionDisposition {
    /// Delete after session closure.
    Ephemeral,
    /// Retain until caller-supplied timestamp.
    Until(u64),
    /// Retain without policy expiry.
    Permanent,
    /// Retain locally under named legal hold.
    LegalHold(String),
    /// Already expired at evaluation time.
    Expired,
}

/// One deterministic share-manifest row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ShareEntry {
    /// Stable field path.
    pub path: String,
    /// Original sensitivity label.
    pub sensitivity: Sensitivity,
    /// Whether omission weakens required replay.
    pub required_for_replay: bool,
    /// Publishable projection.
    pub disclosure: Disclosure,
    /// Independent retention action.
    pub retention: RetentionDisposition,
}

/// Bounded, deterministic share decision. It cannot promote evidence by
/// itself; missing required material demotes completeness/promotion while an
/// intentional redaction leaves integrity intact.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ShareManifest {
    /// Policy schema version.
    pub policy_version: u32,
    /// Destination audience.
    pub audience: ShareAudience,
    /// Entries sorted by path.
    pub entries: Vec<ShareEntry>,
    /// Whether all required replay material is disclosed.
    pub completeness: EvidenceCompleteness,
    /// Redaction is policy, not corruption; unsafe input is refused instead.
    pub integrity: EvidenceIntegrity,
    /// Share admission never grants promotion and demotes incomplete replay.
    pub promotion: PromotionEffect,
    /// Missing license/export authority yields unsupported; other redactions
    /// yield degraded support.
    pub support: OperationalSupport,
}

/// Refusal while constructing a share manifest.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ShareError {
    /// Input exceeds a declared bound.
    BudgetExceeded {
        /// Bound category.
        field: &'static str,
        /// Configured bound.
        limit: usize,
        /// Observed demand.
        observed: usize,
    },
    /// Two values claim the same structured path.
    DuplicatePath {
        /// Duplicate path.
        path: String,
    },
    /// Requested correlation leaks dictionary-testable identity for the label.
    UnsafeCorrelation {
        /// Field path.
        path: String,
        /// Protected label.
        sensitivity: Sensitivity,
        /// Rejected method.
        method: &'static str,
    },
}

impl fmt::Display for ShareError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::BudgetExceeded {
                field,
                limit,
                observed,
            } => write!(f, "share {field} budget exceeded: {observed} > {limit}"),
            Self::DuplicatePath { path } => write!(f, "duplicate share field path {path:?}"),
            Self::UnsafeCorrelation {
                path,
                sensitivity,
                method,
            } => write!(
                f,
                "unsafe {method} correlation for {} field {path:?}",
                sensitivity.name()
            ),
        }
    }
}

impl core::error::Error for ShareError {}

/// Evaluate labeled fields under one bounded share request.
///
/// # Errors
/// [`ShareError`] on budget overflow, duplicate paths, or unsafe correlation.
pub fn evaluate_share(
    mut fields: Vec<LabeledField>,
    request: &ShareRequest,
) -> Result<ShareManifest, ShareError> {
    if fields.len() > request.max_fields {
        return Err(ShareError::BudgetExceeded {
            field: "field_count",
            limit: request.max_fields,
            observed: fields.len(),
        });
    }
    let mut total_bytes = 0usize;
    for field in &fields {
        total_bytes =
            total_bytes
                .checked_add(field.value.len())
                .ok_or(ShareError::BudgetExceeded {
                    field: "input_bytes",
                    limit: request.max_input_bytes,
                    observed: usize::MAX,
                })?;
    }
    if total_bytes > request.max_input_bytes {
        return Err(ShareError::BudgetExceeded {
            field: "input_bytes",
            limit: request.max_input_bytes,
            observed: total_bytes,
        });
    }

    fields.sort_by(|left, right| left.path.cmp(&right.path));
    if let Some(pair) = fields.windows(2).find(|pair| pair[0].path == pair[1].path) {
        return Err(ShareError::DuplicatePath {
            path: pair[0].path.clone(),
        });
    }

    let mut entries = Vec::with_capacity(fields.len());
    let mut required_missing = false;
    let mut support = OperationalSupport::Supported;
    for field in fields {
        let retention = retention_disposition(&field.policy.retention, request.now_ns);
        let block = realm_block(&field.policy, request, &retention);
        let disclosure = if let Some(block) = block {
            support = match block {
                ShareBlock::License | ShareBlock::Export => OperationalSupport::Unsupported,
                ShareBlock::Expired | ShareBlock::Sensitivity | ShareBlock::Credential => {
                    degrade(support)
                }
            };
            Disclosure::Omitted { reason: block }
        } else if sensitivity_admitted(field.policy.sensitivity, request) {
            Disclosure::Plain(field.value)
        } else {
            let reason = if field.policy.sensitivity == Sensitivity::Credential {
                ShareBlock::Credential
            } else {
                ShareBlock::Sensitivity
            };
            let correlation = correlation_for(
                &field.path,
                field.policy.sensitivity,
                &field.value,
                &request.correlation,
            )?;
            support = degrade(support);
            Disclosure::Redacted {
                marker: field.policy.sensitivity.redaction_marker(),
                correlation,
                reason,
            }
        };
        if field.required_for_replay && !matches!(&disclosure, Disclosure::Plain(_)) {
            required_missing = true;
        }
        entries.push(ShareEntry {
            path: field.path,
            sensitivity: field.policy.sensitivity,
            required_for_replay: field.required_for_replay,
            disclosure,
            retention,
        });
    }

    Ok(ShareManifest {
        policy_version: SHARE_POLICY_VERSION,
        audience: request.audience,
        entries,
        completeness: if required_missing {
            EvidenceCompleteness::Incomplete
        } else {
            EvidenceCompleteness::Complete
        },
        integrity: EvidenceIntegrity::Intact,
        promotion: if required_missing {
            PromotionEffect::Demoted
        } else {
            PromotionEffect::Unchanged
        },
        support,
    })
}

fn retention_disposition(class: &RetentionClass, now_ns: u64) -> RetentionDisposition {
    match class {
        RetentionClass::Ephemeral => RetentionDisposition::Ephemeral,
        RetentionClass::Until(deadline) if now_ns >= *deadline => RetentionDisposition::Expired,
        RetentionClass::Until(deadline) => RetentionDisposition::Until(*deadline),
        RetentionClass::Permanent => RetentionDisposition::Permanent,
        RetentionClass::LegalHold(hold) => RetentionDisposition::LegalHold(hold.clone()),
    }
}

fn realm_block(
    policy: &FieldPolicy,
    request: &ShareRequest,
    retention: &RetentionDisposition,
) -> Option<ShareBlock> {
    if matches!(retention, RetentionDisposition::Expired) {
        return Some(ShareBlock::Expired);
    }
    match &policy.license {
        LicenseRealm::Redistributable => {}
        LicenseRealm::Restricted(realm) => match request.audience {
            ShareAudience::Local => {}
            ShareAudience::Organization if request.licensed_realms.binary_search(realm).is_ok() => {
            }
            ShareAudience::Organization | ShareAudience::Public => {
                return Some(ShareBlock::License);
            }
        },
        LicenseRealm::LocalOnly if request.audience != ShareAudience::Local => {
            return Some(ShareBlock::License);
        }
        LicenseRealm::LocalOnly => {}
    }
    match &policy.export {
        ExportRealm::Unrestricted => {}
        ExportRealm::Controlled(realm)
            if request.export_authorizations.binary_search(realm).is_ok() => {}
        ExportRealm::Controlled(_) => return Some(ShareBlock::Export),
        ExportRealm::LocalOnly if request.audience != ShareAudience::Local => {
            return Some(ShareBlock::Export);
        }
        ExportRealm::LocalOnly => {}
    }
    None
}

fn sensitivity_admitted(sensitivity: Sensitivity, request: &ShareRequest) -> bool {
    match sensitivity {
        Sensitivity::Public => true,
        Sensitivity::Internal => !matches!(request.audience, ShareAudience::Public),
        Sensitivity::PersonalData => {
            matches!(request.audience, ShareAudience::Local)
                && request.privileged
                && request.allow_personal_data
        }
        Sensitivity::Secret => {
            matches!(request.audience, ShareAudience::Local)
                && request.privileged
                && request.allow_secrets
        }
        Sensitivity::Credential => false,
    }
}

fn correlation_for(
    path: &str,
    sensitivity: Sensitivity,
    value: &[u8],
    policy: &CorrelationPolicy,
) -> Result<Option<CorrelationDisclosure>, ShareError> {
    match policy {
        CorrelationPolicy::None => Ok(None),
        CorrelationPolicy::UnsaltedPublic => match sensitivity {
            Sensitivity::Public | Sensitivity::Internal => Ok(Some(
                CorrelationDisclosure::UnsaltedFnv1a64(format!("{:016x}", crate::fnv1a64(value))),
            )),
            Sensitivity::PersonalData | Sensitivity::Secret | Sensitivity::Credential => {
                Err(ShareError::UnsafeCorrelation {
                    path: path.to_string(),
                    sensitivity,
                    method: "unsalted",
                })
            }
        },
        CorrelationPolicy::External(token) => match (&token.method, sensitivity) {
            (ExternalCorrelationMethod::Salted { .. }, Sensitivity::PersonalData)
            | (ExternalCorrelationMethod::Salted { .. }, Sensitivity::Secret)
            | (ExternalCorrelationMethod::Salted { .. }, Sensitivity::Credential) => {
                Err(ShareError::UnsafeCorrelation {
                    path: path.to_string(),
                    sensitivity,
                    method: "salted_dictionary_testable",
                })
            }
            (ExternalCorrelationMethod::Keyed { .. }, Sensitivity::Credential) => {
                Err(ShareError::UnsafeCorrelation {
                    path: path.to_string(),
                    sensitivity,
                    method: "credential_correlation",
                })
            }
            (ExternalCorrelationMethod::Salted { salt_id }, _) => {
                Ok(Some(CorrelationDisclosure::ExternalSalted {
                    salt_id: salt_id.clone(),
                    token: token.token.clone(),
                }))
            }
            (ExternalCorrelationMethod::Keyed { key_id }, _) => {
                Ok(Some(CorrelationDisclosure::ExternalKeyed {
                    key_id: key_id.clone(),
                    token: token.token.clone(),
                }))
            }
        },
    }
}

const fn degrade(support: OperationalSupport) -> OperationalSupport {
    match support {
        OperationalSupport::Supported => OperationalSupport::Degraded,
        OperationalSupport::Degraded => OperationalSupport::Degraded,
        OperationalSupport::Unsupported => OperationalSupport::Unsupported,
        OperationalSupport::NotApplicable => OperationalSupport::Degraded,
    }
}
