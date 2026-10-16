//! fs-truss — ground-structure truss layout optimization (plan §9.5
//! [S/F], bead 7tv.13): candidate members at scale → plastic-design
//! LP layout with explicit convergence diagnostics → Euler/code-based
//! sizing with catalog snapping → fs-solid rod re-analysis. The
//! steel-and-concrete flagship's engine (§15.2).
//!
//! Layer: L4 (ASCENT).
//! - [`ground`]: node grids + candidate members under fabrication
//!   rules (length bounds, angle sets, neighbor caps), built on
//!   FrankenNetworkx graphs; generation reproducible and ledgered.
//! - [`lp`]: the member-force (plastic design) LP — minimize volume
//!   subject to nodal equilibrium and yield, tension/compression
//!   split — solved by an in-house PDHG (Chambolle–Pock) first-order
//!   primal-dual iteration: sparse-matvec-dominated, deterministic,
//!   warm-startable, with primal/dual objective-separation and KKT residual
//!   tracking (diagnostics, never promoted by themselves; the vetted Michell
//!   closed-form catalogue
//!   comparison is a ledgered row pending the fs-fab oracle spec's
//!   constants — stated, never silently skipped).
//! - [`certificate`]: a separate, bounded cold proof path. It repairs signed
//!   forces through a deterministic member basis, proves an exact equilibrium
//!   point exists with a Neumann interval bound, outward-checks a scaled dual,
//!   and content-binds finite optimum bounds. Rank/conditioning/overflow or
//!   cancellation failures return no finite certificate.
//! - [`sizing`]: continuous areas from yield, EULER buckling floors
//!   for compression members, joint parsimony with MANDATORY
//!   equilibrium re-verification after pruning, catalog UP-snapping
//!   (feasibility preserved by construction) and the member-by-member
//!   code-check audit table (fs-constraint `Code` rows).
//! - [`rodcheck`]: the critical compression member re-analyzed as an
//!   fs-solid Cosserat rod with a seeded imperfection to 1.3× design
//!   load — the tfz.14/tfz.15 global-buckling spot check.

pub mod certificate;
pub mod ground;
pub mod lp;
pub mod rodcheck;
pub mod sizing;

pub use certificate::{
    ArithmeticEnclosureMethod, CertifiedObjectiveBounds, HARD_MAX_CERTIFICATE_ACTIVE_ROWS,
    HARD_MAX_CERTIFICATE_DENSE_ENTRIES, HARD_MAX_CERTIFICATE_MEMBERS,
    HARD_MAX_CERTIFICATE_OPERATIONS, LayoutCertificateError, LayoutCertificateIdentity,
    LayoutCertificateLimits, LayoutCertificateProblem, LayoutCertificateRefusal,
    LayoutCertificateStatus, LayoutOptimalityCertificate, PrimalCorrectionMethod,
};
pub use ground::{
    ESTIMATED_GRAPH_BYTES_PER_MEMBER, ESTIMATED_GRAPH_BYTES_PER_NODE, GroundLimits, GroundRules,
    GroundStructure, TrussConstructionError,
};
pub use lp::{
    LayoutCase, LayoutLimits, LayoutLp, MAX_PDHG_ITERS, PdhgError, PdhgReport, PdhgSettings,
};
pub use rodcheck::rod_buckling_check;
pub use sizing::{CatalogAudit, SizedMember, size_and_snap};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
