//! fs-marquee — admission/status shell for the P2 marquee study.
//!
//! The actual marquee pipeline is a frontier integration lane: raw SDF
//! geometry, CutFEM physics, DWR certificates, ledgered evidence, and
//! LUMEN renders. This crate exists so the workspace can name and gate
//! that lane without pretending the end-to-end runner is already shipped.

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// THE STUDY RUNNER (the smoke tier of the marquee demo): raw-SDF
/// design, CutFEM state solves, self-adjoint shape gradients, the
/// composed per-iteration certificate, bit-replayable traces.
#[cfg(feature = "marquee")]
pub mod study;

/// Current availability of the marquee lane.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MarqueeStatus {
    /// The feature gate is disabled; no runner is exposed.
    Disabled,
    /// The feature gate is enabled and the SMOKE-TIER runner
    /// ([`study::run_study`]) is available; the full-resolution nightly
    /// golden lane remains a no-claim boundary until it lands in CI.
    SmokeRunnerAvailable,
}

/// Return the current status of the marquee lane.
#[must_use]
pub const fn status() -> MarqueeStatus {
    if cfg!(feature = "marquee") {
        MarqueeStatus::SmokeRunnerAvailable
    } else {
        MarqueeStatus::Disabled
    }
}

/// Human-readable scope statement for agent diagnostics and ledgers.
#[must_use]
pub const fn scope_summary() -> &'static str {
    "P2 marquee lane: raw SDF -> CutFEM/DWR/evidence/render integration; \
     smoke runner shipped, nightly golden pending"
}

#[cfg(test)]
mod tests {
    use super::{MarqueeStatus, VERSION, scope_summary, status};

    #[test]
    fn version_is_stamped() {
        assert!(!VERSION.is_empty());
    }

    #[test]
    fn status_matches_feature_gate() {
        let expected = if cfg!(feature = "marquee") {
            MarqueeStatus::SmokeRunnerAvailable
        } else {
            MarqueeStatus::Disabled
        };
        assert_eq!(status(), expected);
        assert!(scope_summary().contains("nightly golden pending"));
    }
}
