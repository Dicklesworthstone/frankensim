//! Authenticated session grants (bead aeq7, increment 1).
//!
//! [`crate::CapabilityToken`] is caller-declared DATA — bounded and
//! validated, but freely constructible, so it can never be authority
//! by itself. This module adds the authority boundary: an
//! [`IssuerPolicy`] (injected, deny-all by default) evaluates a token
//! as a REQUEST and mints an opaque [`SessionGrant`] whose fields are
//! private and whose canonical digest binds issuer identity, policy
//! fingerprint, session, ledger scope, exact operator set, budgets,
//! issuance/expiry, and the policy's revocation generation.
//!
//! Dynamic enforcement starts here too: a [`CoreLeaseBook`] meters
//! concurrent cores per session and refuses ungranted verbs at lease
//! acquisition, so execution cannot exceed admitted concurrency or run
//! an operator the grant never named.
//!
//! Increment scope (code-first): the grant/verifier/lease core plus
//! unit tests. Governor consumption (`open_session` taking a grant,
//! per-operation lease acquisition on the metering path) is the next
//! increment, coordinated with the idempotency lane. Until then the
//! no-claim boundary is explicit: minting a grant does not yet gate
//! the legacy token path.

use crate::SessionError;
use crate::token::{CapabilityToken, SessionId};

/// Maximum bytes for issuer id and policy fingerprint fields.
pub const MAX_ISSUER_FIELD_BYTES: usize = 128;

/// Domain for the canonical grant digest.
const GRANT_DIGEST_DOMAIN: &str = "org.frankensim.fs-session.session-grant.v1";

/// Who signed off, and under which exact policy revision.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IssuerIdentity {
    issuer_id: String,
    policy_fingerprint: String,
}

impl IssuerIdentity {
    /// Validated construction: both fields bounded ASCII-graphic.
    ///
    /// # Errors
    /// [`SessionError::InvalidIssuerField`] otherwise.
    pub fn new(issuer_id: &str, policy_fingerprint: &str) -> Result<IssuerIdentity, SessionError> {
        for (field, value) in [
            ("issuer_id", issuer_id),
            ("policy_fingerprint", policy_fingerprint),
        ] {
            if value.is_empty()
                || value.len() > MAX_ISSUER_FIELD_BYTES
                || !value.bytes().all(|byte| byte.is_ascii_graphic())
            {
                return Err(SessionError::InvalidIssuerField {
                    field,
                    observed_bytes: value.len(),
                });
            }
        }
        Ok(IssuerIdentity {
            issuer_id: issuer_id.to_string(),
            policy_fingerprint: policy_fingerprint.to_string(),
        })
    }

    /// The stable issuer identifier.
    #[must_use]
    pub fn issuer_id(&self) -> &str {
        &self.issuer_id
    }

    /// The exact policy revision fingerprint.
    #[must_use]
    pub fn policy_fingerprint(&self) -> &str {
        &self.policy_fingerprint
    }
}

/// One policy evaluation outcome.
#[derive(Debug, Clone, PartialEq)]
pub enum PolicyDecision {
    /// Grant, valid until the given ledger-time expiry.
    Granted {
        /// Exclusive expiry bound, ledger nanoseconds.
        expiry_ns: i64,
    },
    /// Refusal with a teaching reason.
    Denied {
        /// Why the request is not entitled to its asks.
        reason: String,
    },
}

/// The injected authority that turns requests into grants. The
/// default posture is deny-all ([`NoIssuerPolicy`]); production
/// deployments inject their own registry-backed policy. Signature
/// verification over external issuer messages is FUTURE scope — this
/// trait is the boundary where it lands.
pub trait IssuerPolicy: Send + Sync {
    /// The issuer identity this policy answers for.
    fn issuer(&self) -> &IssuerIdentity;

    /// Current revocation generation: bumping it invalidates every
    /// grant minted under earlier generations.
    fn revocation_generation(&self) -> u64;

    /// Evaluate a request against the policy.
    fn evaluate(&self, request: &CapabilityToken, issuance_ns: i64) -> PolicyDecision;
}

/// The deny-all default: no caller is entitled to anything.
pub struct NoIssuerPolicy {
    identity: IssuerIdentity,
}

impl NoIssuerPolicy {
    /// The default deny-all policy.
    #[must_use]
    pub fn new() -> NoIssuerPolicy {
        NoIssuerPolicy {
            identity: IssuerIdentity {
                issuer_id: "deny-all".to_string(),
                policy_fingerprint: "deny-all".to_string(),
            },
        }
    }
}

impl Default for NoIssuerPolicy {
    fn default() -> Self {
        Self::new()
    }
}

impl IssuerPolicy for NoIssuerPolicy {
    fn issuer(&self) -> &IssuerIdentity {
        &self.identity
    }

    fn revocation_generation(&self) -> u64 {
        0
    }

    fn evaluate(&self, _request: &CapabilityToken, _issuance_ns: i64) -> PolicyDecision {
        PolicyDecision::Denied {
            reason: "deny-all default policy: inject a deployment issuer policy to mint grants"
                .to_string(),
        }
    }
}

/// An opaque admitted grant. Private fields: only [`mint_grant`] can
/// construct one, so holding a `SessionGrant` IS the proof that an
/// injected policy admitted this exact request.
#[derive(Debug, Clone, PartialEq)]
pub struct SessionGrant {
    issuer: IssuerIdentity,
    session: SessionId,
    ledger_scope: String,
    ops: Vec<String>,
    core_s: f64,
    mem_bytes: u64,
    wall_s: f64,
    cores: u64,
    issuance_ns: i64,
    expiry_ns: i64,
    revocation_generation: u64,
    digest: String,
}

/// Mint a grant from an untrusted request through an injected policy.
///
/// The request's structural validation (bounded canonical operator
/// grants, canonical ledger scope, finite non-negative budgets) runs
/// first; the policy then decides entitlement. The minted grant's
/// digest canonically binds every admitted field under
/// `org.frankensim.fs-session.session-grant.v1`.
///
/// # Errors
/// Structural refusals from the token validators;
/// [`SessionError::GrantDenied`] when the policy refuses;
/// [`SessionError::InvalidResource`] for non-finite budgets or an
/// expiry at/before issuance.
pub fn mint_grant(
    policy: &dyn IssuerPolicy,
    request: &CapabilityToken,
    issuance_ns: i64,
) -> Result<SessionGrant, SessionError> {
    request.validate_operator_grants()?;
    CapabilityToken::validate_ledger_scope(&request.ledger_scope)?;
    for (resource, value) in [
        ("grant core-seconds", request.core_s),
        ("grant wall-seconds", request.wall_s),
    ] {
        if !value.is_finite() || value < 0.0 {
            return Err(SessionError::InvalidResource {
                resource,
                value,
                requirement: "must be finite and non-negative",
            });
        }
    }
    let expiry_ns = match policy.evaluate(request, issuance_ns) {
        PolicyDecision::Granted { expiry_ns } => expiry_ns,
        PolicyDecision::Denied { reason } => return Err(SessionError::GrantDenied { reason }),
    };
    if expiry_ns <= issuance_ns {
        #[allow(clippy::cast_precision_loss)]
        return Err(SessionError::InvalidResource {
            resource: "grant expiry",
            value: expiry_ns as f64,
            requirement: "must lie strictly after issuance",
        });
    }
    let mut ops = request.ops.clone();
    ops.sort_unstable();
    let mut grant = SessionGrant {
        issuer: policy.issuer().clone(),
        session: request.session,
        ledger_scope: request.ledger_scope.clone(),
        ops,
        core_s: request.core_s,
        mem_bytes: request.mem_bytes,
        wall_s: request.wall_s,
        cores: request.cores,
        issuance_ns,
        expiry_ns,
        revocation_generation: policy.revocation_generation(),
        digest: String::new(),
    };
    grant.digest = grant.compute_digest();
    Ok(grant)
}

impl SessionGrant {
    fn compute_digest(&self) -> String {
        let mut preimage = Vec::new();
        let mut push = |bytes: &[u8]| {
            preimage.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
            preimage.extend_from_slice(bytes);
        };
        push(self.issuer.issuer_id.as_bytes());
        push(self.issuer.policy_fingerprint.as_bytes());
        push(&self.session.0.to_le_bytes());
        push(self.ledger_scope.as_bytes());
        for op in &self.ops {
            push(op.as_bytes());
        }
        push(&self.core_s.to_bits().to_le_bytes());
        push(&self.mem_bytes.to_le_bytes());
        push(&self.wall_s.to_bits().to_le_bytes());
        push(&self.cores.to_le_bytes());
        push(&self.issuance_ns.to_le_bytes());
        push(&self.expiry_ns.to_le_bytes());
        push(&self.revocation_generation.to_le_bytes());
        fs_blake3::hash_domain(GRANT_DIGEST_DOMAIN, &preimage).to_string()
    }

    /// Re-verify against the issuing policy at `now_ns`: same issuer
    /// and policy revision, unexpired, current revocation generation,
    /// and an intact digest.
    ///
    /// # Errors
    /// [`SessionError::GrantForged`] on issuer/fingerprint/digest
    /// mismatch; [`SessionError::GrantExpired`] past expiry;
    /// [`SessionError::GrantRevoked`] on generation advance.
    pub fn verify_fresh(&self, policy: &dyn IssuerPolicy, now_ns: i64) -> Result<(), SessionError> {
        if self.issuer != *policy.issuer() || self.digest != self.compute_digest() {
            return Err(SessionError::GrantForged {
                session: self.session.0,
            });
        }
        if now_ns >= self.expiry_ns {
            return Err(SessionError::GrantExpired {
                session: self.session.0,
                expiry_ns: self.expiry_ns,
                now_ns,
            });
        }
        if self.revocation_generation != policy.revocation_generation() {
            return Err(SessionError::GrantRevoked {
                session: self.session.0,
                granted_generation: self.revocation_generation,
                current_generation: policy.revocation_generation(),
            });
        }
        Ok(())
    }

    /// Whether the admitted operator set covers `verb` (exact name or
    /// `ns.*` namespace wildcard — same semantics as the request type).
    #[must_use]
    pub fn grants_op(&self, verb: &str) -> bool {
        self.ops.iter().any(|grant| {
            grant == verb
                || grant.strip_suffix(".*").is_some_and(|namespace| {
                    verb.strip_prefix(namespace)
                        .and_then(|rest| rest.strip_prefix('.'))
                        .is_some_and(|tail| !tail.is_empty())
                })
        })
    }

    /// The static-admission data view of the ADMITTED (not requested)
    /// authority.
    #[must_use]
    pub fn to_admission(&self) -> fs_ir::admission::SessionCapability {
        fs_ir::admission::SessionCapability {
            ops: self.ops.clone(),
            cores: self.cores,
            mem_bytes: self.mem_bytes,
            wall_s: self.wall_s,
        }
    }

    /// The bound session.
    #[must_use]
    pub fn session(&self) -> SessionId {
        self.session
    }

    /// The canonical digest binding every admitted field.
    #[must_use]
    pub fn digest(&self) -> &str {
        &self.digest
    }

    /// The issuing identity.
    #[must_use]
    pub fn issuer(&self) -> &IssuerIdentity {
        &self.issuer
    }

    /// Concurrent-core grant.
    #[must_use]
    pub fn cores(&self) -> u64 {
        self.cores
    }
}

/// Shared concurrent-core accounting: execution acquires a lease per
/// operation and cannot exceed the grant's concurrency or run an
/// operator the grant never named. Leases release on drop.
#[derive(Debug, Default)]
pub struct CoreLeaseBook {
    active: std::sync::Mutex<std::collections::BTreeMap<u64, u64>>,
}

/// One active core lease (drop to release).
#[derive(Debug)]
pub struct CoreLease<'a> {
    book: &'a CoreLeaseBook,
    session: u64,
    cores: u64,
}

impl CoreLeaseBook {
    /// An empty book.
    #[must_use]
    pub fn new() -> CoreLeaseBook {
        CoreLeaseBook::default()
    }

    /// Acquire `cores` for `verb` under a grant re-verified against
    /// the issuing policy at `now_ns`.
    ///
    /// # Errors
    /// Freshness refusals from [`SessionGrant::verify_fresh`];
    /// [`SessionError::UngrantedVerb`] when the admitted set does not
    /// cover `verb`; [`SessionError::CoreLeaseExceeded`] when the
    /// session's active cores plus this ask exceed the grant.
    ///
    /// # Panics
    /// Only on a poisoned internal mutex (a prior panic mid-update).
    pub fn acquire<'a>(
        &'a self,
        grant: &SessionGrant,
        policy: &dyn IssuerPolicy,
        verb: &str,
        cores: u64,
        now_ns: i64,
    ) -> Result<CoreLease<'a>, SessionError> {
        grant.verify_fresh(policy, now_ns)?;
        if !grant.grants_op(verb) {
            return Err(SessionError::UngrantedVerb {
                session: grant.session.0,
                verb: verb.to_string(),
            });
        }
        let mut active = self.active.lock().expect("core lease book poisoned");
        let current = active.get(&grant.session.0).copied().unwrap_or(0);
        let next = current.checked_add(cores);
        match next {
            Some(next) if next <= grant.cores => {
                active.insert(grant.session.0, next);
                Ok(CoreLease {
                    book: self,
                    session: grant.session.0,
                    cores,
                })
            }
            _ => Err(SessionError::CoreLeaseExceeded {
                session: grant.session.0,
                granted: grant.cores,
                active: current,
                requested: cores,
            }),
        }
    }

    /// Currently leased cores for a session.
    ///
    /// # Panics
    /// Only on a poisoned internal mutex.
    #[must_use]
    pub fn active_cores(&self, session: SessionId) -> u64 {
        self.active
            .lock()
            .expect("core lease book poisoned")
            .get(&session.0)
            .copied()
            .unwrap_or(0)
    }
}

impl Drop for CoreLease<'_> {
    fn drop(&mut self) {
        let mut active = self.book.active.lock().expect("core lease book poisoned");
        if let Some(entry) = active.get_mut(&self.session) {
            *entry = entry.saturating_sub(self.cores);
            if *entry == 0 {
                active.remove(&self.session);
            }
        }
    }
}
