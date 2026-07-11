//! Capability tokens: the EXPLICIT grant every IR program executes under
//! — operators (globs), core-seconds, resident memory, wall time, ledger
//! scope. Admission checks the token statically (fs-ir admission consumes
//! the bridge type); the governor meters against it continuously.

use fs_ir::admission::SessionCapability;

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
    pub cores: f64,
    /// Ledger scope (branch/namespace this session may write).
    pub ledger_scope: String,
}

impl CapabilityToken {
    /// The static-admission view of this token. fs-ir's current capability
    /// vocabulary represents declared memory asks as `f64`, so this is a coarse
    /// planning projection above 2^53 bytes. The enforcing governor retains and
    /// compares the original integer byte grant exactly.
    #[must_use]
    pub fn to_admission(&self) -> SessionCapability {
        SessionCapability {
            ops: self.ops.clone(),
            cores: self.cores,
            #[allow(clippy::cast_precision_loss)] // fs-ir's admission vocabulary is f64
            mem_bytes: self.mem_bytes as f64,
            wall_s: self.wall_s,
        }
    }

    /// Does the token grant an operator?
    #[must_use]
    pub fn grants_op(&self, verb: &str) -> bool {
        self.ops.iter().any(|p| {
            p.strip_suffix('*')
                .map_or(p == verb, |prefix| verb.starts_with(prefix))
        })
    }
}
