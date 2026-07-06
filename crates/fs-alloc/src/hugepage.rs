//! Hugepage policy: request 2 MiB-aligned, THP-eligible chunks where the
//! platform can honor it, and RECORD the decision either way (plan §5.1
//! consequence 6: "graceful fallback with the choice recorded").
//!
//! Honesty boundary (see CONTRACT.md no-claims): fs-alloc never issues
//! `madvise` — Decalogue P1 forbids FFI and std exposes no page-attribute
//! control — so on Linux the crate can only make chunks THP-*eligible*
//! (size + alignment) and record the kernel's configured THP mode. Actual
//! backing is the kernel's choice and is NOT claimed. On Apple platforms
//! the base page is 16 KiB and no user-space THP control exists; the
//! decision records that too.

use std::fmt::Write as _;

/// Transparent-hugepage size targeted on Linux.
pub const HUGEPAGE_BYTES: usize = 2 * 1024 * 1024;

/// Caller intent, set in `ArenaConfig`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum HugepagePolicy {
    /// Use 2 MiB-aligned chunks when the chunk size and platform allow it.
    #[default]
    Auto,
    /// Never attempt hugepage-eligible chunks.
    Never,
}

impl HugepagePolicy {
    fn name(self) -> &'static str {
        match self {
            HugepagePolicy::Auto => "auto",
            HugepagePolicy::Never => "never",
        }
    }
}

/// What the probe concluded (recorded once per `ArenaPool`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HugepageOutcome {
    /// Linux, THP mode `always`: chunks >= 2 MiB are allocated 2 MiB-aligned
    /// and are THP-eligible. Backing is still the kernel's choice.
    AlignedForThp,
    /// Linux, but THP is configured `madvise`/`never` (or unreadable): since
    /// fs-alloc issues no madvise (P1: no FFI), alignment would buy nothing;
    /// chunks use the ordinary 128-byte base alignment.
    ThpNotEnabled,
    /// Not a platform with user-space transparent-hugepage control
    /// (e.g. Apple Silicon's 16 KiB base pages).
    UnsupportedPlatform,
    /// Policy was `Never`, or the configured chunk size is below 2 MiB.
    NotRequested,
}

impl HugepageOutcome {
    fn name(self) -> &'static str {
        match self {
            HugepageOutcome::AlignedForThp => "aligned_for_thp",
            HugepageOutcome::ThpNotEnabled => "thp_not_enabled",
            HugepageOutcome::UnsupportedPlatform => "unsupported_platform",
            HugepageOutcome::NotRequested => "not_requested",
        }
    }
}

/// The recorded hugepage decision: policy in, outcome + human/agent-readable
/// reason out. Ledger-bound via [`HugepageDecision::to_json`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HugepageDecision {
    /// The policy the pool was configured with.
    pub policy: HugepagePolicy,
    /// What was decided.
    pub outcome: HugepageOutcome,
    /// Why — one teaching sentence (deterministic per machine config).
    pub detail: String,
}

impl HugepageDecision {
    /// Probe once for a pool configured with `policy` and `chunk_bytes`.
    pub(crate) fn probe(policy: HugepagePolicy, chunk_bytes: usize) -> Self {
        let (outcome, detail) = match policy {
            HugepagePolicy::Never => (
                HugepageOutcome::NotRequested,
                "policy=never; chunks use the 128-byte base alignment".to_string(),
            ),
            HugepagePolicy::Auto if chunk_bytes < HUGEPAGE_BYTES => (
                HugepageOutcome::NotRequested,
                format!(
                    "chunk_bytes {chunk_bytes} < 2 MiB threshold; raise \
                     ArenaConfig::chunk_bytes to opt in"
                ),
            ),
            HugepagePolicy::Auto => probe_platform(),
        };
        HugepageDecision {
            policy,
            outcome,
            detail,
        }
    }

    /// Chunk base alignment implied by this decision for a chunk of
    /// `bytes` bytes.
    pub(crate) fn chunk_align(&self, bytes: usize) -> usize {
        if self.outcome == HugepageOutcome::AlignedForThp && bytes >= HUGEPAGE_BYTES {
            HUGEPAGE_BYTES
        } else {
            crate::ALLOC_ALIGN
        }
    }

    /// Canonical JSON object (deterministic field order) for ledger rows and
    /// `fs_obs::EventKind::Custom` payloads.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::with_capacity(96);
        let _ = write!(
            s,
            "{{\"policy\":\"{}\",\"outcome\":\"{}\",\"detail\":\"",
            self.policy.name(),
            self.outcome.name()
        );
        for c in self.detail.chars() {
            match c {
                '"' => s.push_str("\\\""),
                '\\' => s.push_str("\\\\"),
                '\n' => s.push_str("\\n"),
                c => s.push(c),
            }
        }
        s.push_str("\"}");
        s
    }
}

#[cfg(target_os = "linux")]
fn probe_platform() -> (HugepageOutcome, String) {
    match std::fs::read_to_string("/sys/kernel/mm/transparent_hugepage/enabled") {
        Ok(mode) if mode.contains("[always]") => (
            HugepageOutcome::AlignedForThp,
            "THP mode=always; chunks >= 2 MiB are 2 MiB-aligned and THP-eligible \
             (backing is the kernel's choice; not claimed)"
                .to_string(),
        ),
        Ok(mode) => (
            HugepageOutcome::ThpNotEnabled,
            format!(
                "THP mode is {} and fs-alloc issues no madvise (Decalogue P1: no FFI); \
                 chunks use the 128-byte base alignment",
                mode.trim()
            ),
        ),
        Err(e) => (
            HugepageOutcome::ThpNotEnabled,
            format!("THP sysfs unreadable ({e}); chunks use the 128-byte base alignment"),
        ),
    }
}

#[cfg(all(target_os = "macos", target_arch = "aarch64"))]
fn probe_platform() -> (HugepageOutcome, String) {
    (
        HugepageOutcome::UnsupportedPlatform,
        "Apple Silicon uses 16 KiB base pages and exposes no user-space transparent-hugepage \
         control via safe std; chunks use the 128-byte base alignment"
            .to_string(),
    )
}

#[cfg(not(any(target_os = "linux", all(target_os = "macos", target_arch = "aarch64"))))]
fn probe_platform() -> (HugepageOutcome, String) {
    (
        HugepageOutcome::UnsupportedPlatform,
        "no transparent-hugepage control on this platform; chunks use the 128-byte \
         base alignment"
            .to_string(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn never_and_small_chunks_are_not_requested() {
        let d = HugepageDecision::probe(HugepagePolicy::Never, 64 << 20);
        assert_eq!(d.outcome, HugepageOutcome::NotRequested);
        assert_eq!(d.chunk_align(64 << 20), crate::ALLOC_ALIGN);

        let d = HugepageDecision::probe(HugepagePolicy::Auto, HUGEPAGE_BYTES - 1);
        assert_eq!(d.outcome, HugepageOutcome::NotRequested);
        assert!(d.detail.contains("threshold"), "teaching detail: {d:?}");
    }

    #[test]
    fn auto_probe_records_a_platform_verdict() {
        let d = HugepageDecision::probe(HugepagePolicy::Auto, HUGEPAGE_BYTES);
        // Platform-dependent outcome, but never NotRequested at this size,
        // and always with a non-empty recorded reason.
        assert_ne!(d.outcome, HugepageOutcome::NotRequested);
        assert!(!d.detail.is_empty());
        // Alignment follows the outcome.
        let align = d.chunk_align(HUGEPAGE_BYTES);
        if d.outcome == HugepageOutcome::AlignedForThp {
            assert_eq!(align, HUGEPAGE_BYTES);
        } else {
            assert_eq!(align, crate::ALLOC_ALIGN);
        }
    }

    #[test]
    fn json_is_canonical_and_escaped() {
        let d = HugepageDecision {
            policy: HugepagePolicy::Auto,
            outcome: HugepageOutcome::NotRequested,
            detail: "quote\" backslash\\".to_string(),
        };
        let j = d.to_json();
        assert!(j.starts_with("{\"policy\":\"auto\",\"outcome\":\"not_requested\","));
        assert!(j.contains("quote\\\" backslash\\\\"));
    }
}
