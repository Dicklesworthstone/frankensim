//! Capability tokens: the EXPLICIT grant every IR program executes under
//! — operators (globs), core-seconds, resident memory, wall time, ledger
//! scope. Admission checks the token statically (fs-ir admission consumes
//! the bridge type); the governor meters against it continuously.

use fs_ir::admission::SessionCapability;
use std::collections::BTreeSet;

use crate::SessionError;

/// Maximum byte length of one canonical ledger scope.
pub const MAX_LEDGER_SCOPE_BYTES: usize = 128;
/// Maximum operator grants carried by one token.
pub const MAX_CAPABILITY_OPS: usize = 256;
/// Maximum bytes in one canonical operator grant.
pub const MAX_CAPABILITY_OP_BYTES: usize = 128;
/// Maximum aggregate operator-grant bytes carried by one token.
pub const MAX_CAPABILITY_TOTAL_OP_BYTES: usize = 8 * 1024;
const LEDGER_SCOPE_REQUIREMENT: &str =
    "must contain 1..=128 ASCII graphic bytes (0x21..=0x7e), with no whitespace or controls";
const OPERATOR_GRANT_REQUIREMENT: &str = "must be a 1..=128 byte canonical exact operator name or namespace wildcard ending in .*, using only ASCII letters, digits, '.', '_', and '-'";

fn scope_error_preview(scope: &str) -> String {
    let mut end = 0;
    for (index, ch) in scope.char_indices() {
        let next = index + ch.len_utf8();
        if next > MAX_LEDGER_SCOPE_BYTES {
            break;
        }
        end = next;
    }
    scope[..end].to_string()
}

/// A session identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SessionId(pub u64);

/// The capability token: explicit, bounded, ledger-scoped.
#[derive(Debug, Clone, PartialEq)]
pub struct CapabilityToken {
    /// The session this token is bound to.
    pub session: SessionId,
    /// Granted operator globs (`flux.*`, `ascent.optimize`, …).
    pub ops: Vec<String>,
    /// Core-seconds grant (CPU time across all cores).
    pub core_s: f64,
    /// Resident-memory grant in bytes.
    pub mem_bytes: u64,
    /// Wall-clock grant in seconds.
    pub wall_s: f64,
    /// Cores the session may occupy at once.
    pub cores: u64,
    /// Canonical ledger scope (exact branch/namespace this session may write).
    /// [`crate::Governor::open_session`] validates it before registration.
    pub ledger_scope: String,
}

impl CapabilityToken {
    /// Validate the bounded canonical operator authority carried by this token.
    pub fn validate_operator_grants(&self) -> Result<(), SessionError> {
        if self.ops.len() > MAX_CAPABILITY_OPS {
            return Err(SessionError::LimitExceeded {
                resource: "capability_operator_grants",
                limit: MAX_CAPABILITY_OPS,
                observed_at_least: self.ops.len(),
            });
        }
        let mut total_bytes = 0usize;
        let mut unique = BTreeSet::new();
        for (index, grant) in self.ops.iter().enumerate() {
            if grant.is_empty()
                || grant.len() > MAX_CAPABILITY_OP_BYTES
                || !fs_ir::admission::valid_operator_pattern(grant)
            {
                return Err(SessionError::InvalidOperatorGrant {
                    index,
                    grant_preview: scope_error_preview(grant),
                    grant_bytes: grant.len(),
                    requirement: OPERATOR_GRANT_REQUIREMENT,
                });
            }
            total_bytes =
                total_bytes
                    .checked_add(grant.len())
                    .ok_or(SessionError::LimitExceeded {
                        resource: "capability_operator_bytes",
                        limit: MAX_CAPABILITY_TOTAL_OP_BYTES,
                        observed_at_least: usize::MAX,
                    })?;
            if total_bytes > MAX_CAPABILITY_TOTAL_OP_BYTES {
                return Err(SessionError::LimitExceeded {
                    resource: "capability_operator_bytes",
                    limit: MAX_CAPABILITY_TOTAL_OP_BYTES,
                    observed_at_least: total_bytes,
                });
            }
            if !unique.insert(grant.as_str()) {
                return Err(SessionError::DuplicateOperatorGrant {
                    grant: grant.clone(),
                });
            }
        }
        Ok(())
    }

    /// Validate the exact ledger namespace carried as session authority.
    ///
    /// Restricting scopes to bounded ASCII graphic bytes makes byte equality
    /// the canonical identity: there are no whitespace, control-character, or
    /// Unicode-normalization aliases for the same apparent namespace.
    pub fn validate_ledger_scope(scope: &str) -> Result<(), SessionError> {
        if scope.is_empty()
            || scope.len() > MAX_LEDGER_SCOPE_BYTES
            || !scope.bytes().all(|byte| byte.is_ascii_graphic())
        {
            return Err(SessionError::InvalidLedgerScope {
                scope_preview: scope_error_preview(scope),
                scope_bytes: scope.len(),
                requirement: LEDGER_SCOPE_REQUIREMENT,
            });
        }
        Ok(())
    }

    /// The static-admission data view of this token. This projection does not
    /// prove that the token was registered or minted by an external issuer.
    /// Memory remains an exact `u64` through the bridge, so static admission and
    /// governor metering compare the same byte value.
    pub fn to_admission(&self) -> Result<SessionCapability, SessionError> {
        self.validate_operator_grants()?;
        Ok(SessionCapability {
            ops: self.ops.clone(),
            cores: self.cores,
            mem_bytes: self.mem_bytes,
            wall_s: self.wall_s,
        })
    }

    /// Does the token grant an operator?
    #[must_use]
    pub fn grants_op(&self, verb: &str) -> bool {
        fs_ir::admission::valid_operator_pattern(verb)
            && !verb.contains('*')
            && self.validate_operator_grants().is_ok()
            && self.ops.iter().any(|p| {
                p.strip_suffix('*')
                    .map_or(p == verb, |prefix| verb.starts_with(prefix))
            })
    }
}
