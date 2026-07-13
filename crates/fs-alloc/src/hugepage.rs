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

use core::mem::size_of;
use std::fmt::Write as _;

/// Transparent-hugepage size targeted on Linux.
pub const HUGEPAGE_BYTES: usize = 2 * 1024 * 1024;

/// Semantic version of the retained hugepage-decision transport.
pub const HUGEPAGE_DECISION_IDENTITY_VERSION: u32 = 1;

/// Exact domain carried by every retained hugepage-decision transport.
pub const HUGEPAGE_DECISION_IDENTITY_DOMAIN: &str = "org.frankensim.fs-alloc.hugepage-decision.v1";

/// Owner-local declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const HUGEPAGE_DECISION_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-alloc:hugepage-decision",
    "version_const=HUGEPAGE_DECISION_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-alloc.hugepage-decision.v1",
    "domain_const=HUGEPAGE_DECISION_IDENTITY_DOMAIN",
    "encoder=HugepageDecision::to_canonical_bytes",
    "encoder_helpers=hugepage_decision_canonical_bytes_with_schema,append_identity_u32,append_identity_u64",
    "schema_constants=HUGEPAGE_DECISION_IDENTITY_VERSION,HUGEPAGE_DECISION_IDENTITY_DOMAIN",
    "schema_functions=HugepagePolicy::identity_tag,HugepagePolicy::from_identity_tag,HugepagePolicy::from_name,HugepageOutcome::identity_tag,HugepageOutcome::from_identity_tag,HugepageOutcome::from_name,HugepageDecision::to_json,HugepageDecision::from_json,read_canonical_json_string,read_identity_u32,read_identity_u64,take_identity_bytes",
    "schema_dependencies=none",
    "digest=none-exact-canonical-transport",
    "encoding=canonical-transport-exact-bits",
    "sources=HugepageDecision",
    "source_fields=HugepageDecision.policy:semantic,HugepageDecision.outcome:semantic,HugepageDecision.detail:semantic",
    "source_bindings=HugepageDecision.policy>policy-tag,HugepageDecision.outcome>outcome-tag,HugepageDecision.detail>detail-byte-count+detail-utf8",
    "external_semantic_fields=artifact-domain,domain-byte-count,identity-version",
    "semantic_fields=artifact-domain,domain-byte-count,identity-version,policy-tag,outcome-tag,detail-byte-count,detail-utf8",
    "excluded_fields=none",
    "consumers=HugepageDecision::to_canonical_bytes,HugepageDecision::from_canonical_bytes,HugepageDecision::to_json,HugepageDecision::from_json,ArenaPool::hugepage_decision,fs-exec:tilepool-placement",
    "mutations=artifact-domain:crates/fs-alloc/src/hugepage.rs#hugepage_decision_identity_fields_move_independently,domain-byte-count:crates/fs-alloc/src/hugepage.rs#hugepage_decision_identity_fields_move_independently,identity-version:crates/fs-alloc/src/hugepage.rs#hugepage_decision_identity_fields_move_independently,policy-tag:crates/fs-alloc/src/hugepage.rs#hugepage_decision_identity_fields_move_independently,outcome-tag:crates/fs-alloc/src/hugepage.rs#hugepage_decision_identity_fields_move_independently,detail-byte-count:crates/fs-alloc/src/hugepage.rs#hugepage_decision_identity_fields_move_independently,detail-utf8:crates/fs-alloc/src/hugepage.rs#hugepage_decision_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_hugepage_decision_identity_fields",
    "transport_guard=HugepageDecision::from_canonical_bytes",
    "version_guard=crates/fs-alloc/src/hugepage.rs#hugepage_decision_identity_versions_fail_closed",
    "coupling_surface=fs-alloc:hugepage-decision",
];

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

    const fn identity_tag(self) -> u8 {
        match self {
            HugepagePolicy::Auto => 0,
            HugepagePolicy::Never => 1,
        }
    }

    const fn from_identity_tag(tag: u8) -> Option<Self> {
        match tag {
            0 => Some(Self::Auto),
            1 => Some(Self::Never),
            _ => None,
        }
    }

    fn from_name(name: &str) -> Option<Self> {
        match name {
            "auto" => Some(Self::Auto),
            "never" => Some(Self::Never),
            _ => None,
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

    const fn identity_tag(self) -> u8 {
        match self {
            HugepageOutcome::AlignedForThp => 0,
            HugepageOutcome::ThpNotEnabled => 1,
            HugepageOutcome::UnsupportedPlatform => 2,
            HugepageOutcome::NotRequested => 3,
        }
    }

    const fn from_identity_tag(tag: u8) -> Option<Self> {
        match tag {
            0 => Some(Self::AlignedForThp),
            1 => Some(Self::ThpNotEnabled),
            2 => Some(Self::UnsupportedPlatform),
            3 => Some(Self::NotRequested),
            _ => None,
        }
    }

    fn from_name(name: &str) -> Option<Self> {
        match name {
            "aligned_for_thp" => Some(Self::AlignedForThp),
            "thp_not_enabled" => Some(Self::ThpNotEnabled),
            "unsupported_platform" => Some(Self::UnsupportedPlatform),
            "not_requested" => Some(Self::NotRequested),
            _ => None,
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
    /// Current producer version for the exact retained identity transport.
    pub const IDENTITY_VERSION: u32 = HUGEPAGE_DECISION_IDENTITY_VERSION;

    /// Current exact retained identity domain.
    pub const IDENTITY_DOMAIN: &str = HUGEPAGE_DECISION_IDENTITY_DOMAIN;

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
                '\r' => s.push_str("\\r"),
                '\t' => s.push_str("\\t"),
                c if c.is_control() => {
                    let _ = write!(s, "\\u{:04x}", c as u32);
                }
                c => s.push(c),
            }
        }
        s.push_str("\"}");
        s
    }

    /// Decode and admit the one canonical JSON spelling emitted by
    /// [`Self::to_json`].
    ///
    /// Field reordering, unknown enum names, non-canonical escapes, literal
    /// control characters, trailing bytes, and duplicate fields are refused.
    /// This exact fixed-point decoder is what permits tile-pool placement v2
    /// to bind the JSON bytes without treating presentation output as an
    /// unchecked identity input.
    ///
    /// # Errors
    /// Returns a stable refusal message when `json` is not the exact canonical
    /// transport for one complete decision.
    pub fn from_json(json: &str) -> Result<Self, &'static str> {
        let mut rest = json;
        rest = rest
            .strip_prefix("{\"policy\":")
            .ok_or("hugepage decision JSON has a non-canonical prefix")?;
        let (policy_name, next) = read_canonical_json_string(rest)?;
        rest = next
            .strip_prefix(",\"outcome\":")
            .ok_or("hugepage decision JSON has a non-canonical outcome field")?;
        let (outcome_name, next) = read_canonical_json_string(rest)?;
        rest = next
            .strip_prefix(",\"detail\":")
            .ok_or("hugepage decision JSON has a non-canonical detail field")?;
        let (detail, next) = read_canonical_json_string(rest)?;
        if next != "}" {
            return Err("hugepage decision JSON has trailing or reordered content");
        }
        let decision = Self {
            policy: HugepagePolicy::from_name(&policy_name)
                .ok_or("hugepage decision JSON policy is unsupported")?,
            outcome: HugepageOutcome::from_name(&outcome_name)
                .ok_or("hugepage decision JSON outcome is unsupported")?,
            detail,
        };
        if decision.to_json() != json {
            return Err("hugepage decision JSON is not canonically encoded");
        }
        Ok(decision)
    }

    /// Exact retained transport for the decision's complete semantic state.
    ///
    /// Layout (all counts and integers little-endian): domain byte count,
    /// domain bytes, identity version, policy tag, outcome tag, detail byte
    /// count, and exact UTF-8 detail bytes. The human-readable JSON remains a
    /// separate event format and is not used as a decoder.
    #[must_use]
    pub fn to_canonical_bytes(&self) -> Vec<u8> {
        hugepage_decision_canonical_bytes_with_schema(
            self,
            HUGEPAGE_DECISION_IDENTITY_DOMAIN,
            HUGEPAGE_DECISION_IDENTITY_VERSION,
            HUGEPAGE_DECISION_IDENTITY_DOMAIN.len(),
            self.detail.len(),
        )
    }

    /// Decode and admit an exact retained decision.
    ///
    /// Stale/future versions, foreign domains, unknown tags, malformed UTF-8,
    /// count mismatches, truncation, and suffix bytes are all refused before
    /// the value becomes replay authority.
    ///
    /// # Errors
    /// Returns a stable refusal message for any non-canonical or unsupported
    /// retained transport.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, &'static str> {
        let mut cursor = 0usize;
        let domain_len = usize::try_from(read_identity_u64(bytes, &mut cursor)?)
            .map_err(|_| "hugepage identity domain length exceeds usize")?;
        if domain_len != HUGEPAGE_DECISION_IDENTITY_DOMAIN.len() {
            return Err("hugepage identity domain length is not canonical");
        }
        let domain = take_identity_bytes(bytes, &mut cursor, domain_len)?;
        if domain != HUGEPAGE_DECISION_IDENTITY_DOMAIN.as_bytes() {
            return Err("hugepage identity domain is unsupported");
        }
        let version = read_identity_u32(bytes, &mut cursor)?;
        if version != HUGEPAGE_DECISION_IDENTITY_VERSION {
            return Err("hugepage identity version is unsupported");
        }
        let policy = HugepagePolicy::from_identity_tag(
            *take_identity_bytes(bytes, &mut cursor, 1)?
                .first()
                .ok_or("hugepage identity policy tag is missing")?,
        )
        .ok_or("hugepage identity policy tag is unsupported")?;
        let outcome = HugepageOutcome::from_identity_tag(
            *take_identity_bytes(bytes, &mut cursor, 1)?
                .first()
                .ok_or("hugepage identity outcome tag is missing")?,
        )
        .ok_or("hugepage identity outcome tag is unsupported")?;
        let detail_len = usize::try_from(read_identity_u64(bytes, &mut cursor)?)
            .map_err(|_| "hugepage identity detail length exceeds usize")?;
        let detail_bytes = take_identity_bytes(bytes, &mut cursor, detail_len)?;
        if cursor != bytes.len() {
            return Err("hugepage identity has non-canonical suffix bytes");
        }
        let detail = core::str::from_utf8(detail_bytes)
            .map_err(|_| "hugepage identity detail is not UTF-8")?
            .to_string();
        Ok(Self {
            policy,
            outcome,
            detail,
        })
    }
}

fn read_canonical_json_string(input: &str) -> Result<(String, &str), &'static str> {
    let mut rest = input
        .strip_prefix('"')
        .ok_or("hugepage decision JSON string is missing its opening quote")?;
    let mut value = String::new();
    loop {
        let ch = rest
            .chars()
            .next()
            .ok_or("hugepage decision JSON string is unterminated")?;
        rest = &rest[ch.len_utf8()..];
        match ch {
            '"' => return Ok((value, rest)),
            '\\' => {
                let escaped = rest
                    .chars()
                    .next()
                    .ok_or("hugepage decision JSON escape is truncated")?;
                rest = &rest[escaped.len_utf8()..];
                let decoded = match escaped {
                    '"' => '"',
                    '\\' => '\\',
                    'n' => '\n',
                    'r' => '\r',
                    't' => '\t',
                    'u' => {
                        let hex = rest
                            .get(..4)
                            .ok_or("hugepage decision JSON unicode escape is truncated")?;
                        if !hex.bytes().all(|byte| byte.is_ascii_hexdigit()) {
                            return Err("hugepage decision JSON unicode escape is invalid");
                        }
                        rest = &rest[4..];
                        char::from_u32(
                            u32::from_str_radix(hex, 16)
                                .map_err(|_| "hugepage decision JSON unicode escape is invalid")?,
                        )
                        .ok_or("hugepage decision JSON unicode escape is invalid")?
                    }
                    _ => return Err("hugepage decision JSON escape is unsupported"),
                };
                value.push(decoded);
            }
            c if c.is_control() => {
                return Err("hugepage decision JSON contains a literal control character");
            }
            c => value.push(c),
        }
    }
}

fn hugepage_decision_canonical_bytes_with_schema(
    decision: &HugepageDecision,
    domain: &str,
    version: u32,
    declared_domain_len: usize,
    declared_detail_len: usize,
) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(
        size_of::<u64>() * 2 + size_of::<u32>() + 2 + domain.len() + decision.detail.len(),
    );
    append_identity_u64(&mut bytes, declared_domain_len);
    bytes.extend_from_slice(domain.as_bytes());
    append_identity_u32(&mut bytes, version);
    bytes.push(decision.policy.identity_tag());
    bytes.push(decision.outcome.identity_tag());
    append_identity_u64(&mut bytes, declared_detail_len);
    bytes.extend_from_slice(decision.detail.as_bytes());
    bytes
}

fn append_identity_u32(bytes: &mut Vec<u8>, value: u32) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn append_identity_u64(bytes: &mut Vec<u8>, value: usize) {
    bytes.extend_from_slice(
        &u64::try_from(value)
            .expect("hugepage identity byte count exceeds u64")
            .to_le_bytes(),
    );
}

fn read_identity_u32(bytes: &[u8], cursor: &mut usize) -> Result<u32, &'static str> {
    let raw = take_identity_bytes(bytes, cursor, size_of::<u32>())?;
    Ok(u32::from_le_bytes(
        raw.try_into()
            .map_err(|_| "hugepage identity u32 is truncated")?,
    ))
}

fn read_identity_u64(bytes: &[u8], cursor: &mut usize) -> Result<u64, &'static str> {
    let raw = take_identity_bytes(bytes, cursor, size_of::<u64>())?;
    Ok(u64::from_le_bytes(
        raw.try_into()
            .map_err(|_| "hugepage identity u64 is truncated")?,
    ))
}

fn take_identity_bytes<'a>(
    bytes: &'a [u8],
    cursor: &mut usize,
    len: usize,
) -> Result<&'a [u8], &'static str> {
    let end = cursor
        .checked_add(len)
        .ok_or("hugepage identity byte extent overflows usize")?;
    let value = bytes
        .get(*cursor..end)
        .ok_or("hugepage identity is truncated")?;
    *cursor = end;
    Ok(value)
}

#[allow(dead_code)]
fn classify_hugepage_decision_identity_fields(decision: &HugepageDecision) {
    let HugepageDecision {
        policy,
        outcome,
        detail,
    } = decision;
    let _ = (policy, outcome, detail);
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
            detail: "quote\" backslash\\ newline\n carriage\r tab\t nul\0".to_string(),
        };
        let j = d.to_json();
        assert_eq!(
            j,
            "{\"policy\":\"auto\",\"outcome\":\"not_requested\",\"detail\":\"quote\\\" backslash\\\\ newline\\n carriage\\r tab\\t nul\\u0000\"}"
        );
        assert_eq!(HugepageDecision::from_json(&j), Ok(d));

        for noncanonical in [
            "{\"outcome\":\"not_requested\",\"policy\":\"auto\",\"detail\":\"x\"}",
            "{\"policy\":\"automatic\",\"outcome\":\"not_requested\",\"detail\":\"x\"}",
            "{\"policy\":\"auto\",\"outcome\":\"not_requested\",\"detail\":\"\\u0078\"}",
            "{\"policy\":\"auto\",\"outcome\":\"not_requested\",\"detail\":\"x\"} trailing",
            "{\"policy\":\"auto\",\"outcome\":\"not_requested\",\"detail\":\"literal\nnewline\"}",
        ] {
            assert!(
                HugepageDecision::from_json(noncanonical).is_err(),
                "accepted non-canonical JSON: {noncanonical:?}"
            );
        }
    }

    fn fixture_decision() -> HugepageDecision {
        HugepageDecision {
            policy: HugepagePolicy::Auto,
            outcome: HugepageOutcome::ThpNotEnabled,
            detail: "fixture detail with \\\"quoted\\\" policy".to_string(),
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn hugepage_decision_identity_fields_move_independently() {
        let fixture = fixture_decision();
        let canonical = fixture.to_canonical_bytes();
        let encode = |decision: &HugepageDecision,
                      domain: &str,
                      version: u32,
                      domain_len: usize,
                      detail_len: usize| {
            hugepage_decision_canonical_bytes_with_schema(
                decision, domain, version, domain_len, detail_len,
            )
        };
        let assert_moves = |field: &str, changed: Vec<u8>| {
            assert_ne!(
                changed, canonical,
                "semantic field {field} did not move bytes"
            );
        };

        let foreign_domain = "org.frankensim.fs-alloc.hugepage-decision.w1";
        assert_eq!(
            foreign_domain.len(),
            HUGEPAGE_DECISION_IDENTITY_DOMAIN.len()
        );
        assert_moves(
            "artifact-domain",
            encode(
                &fixture,
                foreign_domain,
                HUGEPAGE_DECISION_IDENTITY_VERSION,
                foreign_domain.len(),
                fixture.detail.len(),
            ),
        );
        assert_moves(
            "domain-byte-count",
            encode(
                &fixture,
                HUGEPAGE_DECISION_IDENTITY_DOMAIN,
                HUGEPAGE_DECISION_IDENTITY_VERSION,
                HUGEPAGE_DECISION_IDENTITY_DOMAIN.len() + 1,
                fixture.detail.len(),
            ),
        );
        assert_moves(
            "identity-version",
            encode(
                &fixture,
                HUGEPAGE_DECISION_IDENTITY_DOMAIN,
                HUGEPAGE_DECISION_IDENTITY_VERSION + 1,
                HUGEPAGE_DECISION_IDENTITY_DOMAIN.len(),
                fixture.detail.len(),
            ),
        );

        let mut changed = fixture.clone();
        changed.policy = HugepagePolicy::Never;
        assert_moves("policy-tag", changed.to_canonical_bytes());
        let mut changed = fixture.clone();
        changed.outcome = HugepageOutcome::UnsupportedPlatform;
        assert_moves("outcome-tag", changed.to_canonical_bytes());
        assert_moves(
            "detail-byte-count",
            encode(
                &fixture,
                HUGEPAGE_DECISION_IDENTITY_DOMAIN,
                HUGEPAGE_DECISION_IDENTITY_VERSION,
                HUGEPAGE_DECISION_IDENTITY_DOMAIN.len(),
                fixture.detail.len() + 1,
            ),
        );
        let mut changed = fixture.clone();
        changed.detail.push('!');
        assert_moves("detail-utf8", changed.to_canonical_bytes());
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn hugepage_decision_identity_versions_fail_closed() {
        let fixture = fixture_decision();
        let canonical = fixture.to_canonical_bytes();
        assert_eq!(
            HugepageDecision::from_canonical_bytes(&canonical),
            Ok(fixture.clone())
        );

        for version in [
            HUGEPAGE_DECISION_IDENTITY_VERSION - 1,
            HUGEPAGE_DECISION_IDENTITY_VERSION + 1,
        ] {
            let retained = hugepage_decision_canonical_bytes_with_schema(
                &fixture,
                HUGEPAGE_DECISION_IDENTITY_DOMAIN,
                version,
                HUGEPAGE_DECISION_IDENTITY_DOMAIN.len(),
                fixture.detail.len(),
            );
            assert!(
                HugepageDecision::from_canonical_bytes(&retained).is_err(),
                "retained producer version {version} must be refused"
            );
        }

        let foreign_domain = hugepage_decision_canonical_bytes_with_schema(
            &fixture,
            "org.frankensim.fs-alloc.hugepage-decision.w1",
            HUGEPAGE_DECISION_IDENTITY_VERSION,
            HUGEPAGE_DECISION_IDENTITY_DOMAIN.len(),
            fixture.detail.len(),
        );
        assert!(HugepageDecision::from_canonical_bytes(&foreign_domain).is_err());

        let wrong_domain_count = hugepage_decision_canonical_bytes_with_schema(
            &fixture,
            HUGEPAGE_DECISION_IDENTITY_DOMAIN,
            HUGEPAGE_DECISION_IDENTITY_VERSION,
            HUGEPAGE_DECISION_IDENTITY_DOMAIN.len() + 1,
            fixture.detail.len(),
        );
        assert!(HugepageDecision::from_canonical_bytes(&wrong_domain_count).is_err());
        let wrong_detail_count = hugepage_decision_canonical_bytes_with_schema(
            &fixture,
            HUGEPAGE_DECISION_IDENTITY_DOMAIN,
            HUGEPAGE_DECISION_IDENTITY_VERSION,
            HUGEPAGE_DECISION_IDENTITY_DOMAIN.len(),
            fixture.detail.len() + 1,
        );
        assert!(HugepageDecision::from_canonical_bytes(&wrong_detail_count).is_err());

        let mut suffixed = canonical.clone();
        suffixed.push(0);
        assert!(HugepageDecision::from_canonical_bytes(&suffixed).is_err());
        for truncated_len in [0, canonical.len() - 1] {
            assert!(HugepageDecision::from_canonical_bytes(&canonical[..truncated_len]).is_err());
        }

        let domain_frame = size_of::<u64>() + HUGEPAGE_DECISION_IDENTITY_DOMAIN.len();
        let tag_offset = domain_frame + size_of::<u32>();
        let mut unknown_policy = canonical.clone();
        unknown_policy[tag_offset] = u8::MAX;
        assert!(HugepageDecision::from_canonical_bytes(&unknown_policy).is_err());
        let mut unknown_outcome = canonical;
        unknown_outcome[tag_offset + 1] = u8::MAX;
        assert!(HugepageDecision::from_canonical_bytes(&unknown_outcome).is_err());

        let mut invalid_utf8 = fixture.to_canonical_bytes();
        *invalid_utf8
            .last_mut()
            .expect("fixture transport has a nonempty detail") = u8::MAX;
        assert!(HugepageDecision::from_canonical_bytes(&invalid_utf8).is_err());
    }
}
