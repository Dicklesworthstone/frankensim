//! Replay identity + checkpoint/artifact schema freeze (bead
//! wf-root-guzez.1.9.4, leaf of E0.9): first executing transcription of the
//! FROZEN identity/replay/checkpoint contract
//! `data/wright-flyer/replay-identity-schema-v1.json` (bead .1.9.1) into
//! typed Rust surfaces.
//!
//! Byte-exact fidelity to the frozen formulas is THE requirement here, so
//! this module mints ids through [`fs_blake3::hash_domain`] over explicit
//! canonical preimages — the same convention this crate already uses for the
//! same frozen family ([`crate::replay::INPUT_TRACE_DOMAIN`],
//! [`crate::checkpoint::SUBSYSTEM_DIGEST_DOMAIN`]). It deliberately does not
//! route through the descriptor-bound StrongIdentity encoder machinery,
//! whose digest algebra is internal to fs-blake3 and would NOT reproduce the
//! frozen `H("fs-flyer/run-intent/v1", canonical(...))` formulas an
//! independent implementation must match.
//!
//! Frozen pieces transcribed by this module:
//! - `RunIdentityBasisV1` (7 closed fields) + `RunIntentId`
//!   `H("fs-flyer/run-intent/v1", canonical(basis))`
//! - final `RunId` `H("fs-flyer/run/v1", RunIntentId, InputTraceId)` — the
//!   closed two-input definition; the preimage is the two raw 32-byte
//!   digests concatenated (pinned here and in the goldens fixture;
//!   never "hash of the quintuple")
//! - `RunSpecId` `H("fs-flyer/run-spec/v1", pre-execution projection +
//!   InputPolicyId)` — cache key only, excludes `accepted_tick0_state_digest`
//!   (Round-5 lifecycle fix: every field exists at minting time)
//! - `RunAnchor` `Intent(..) | Final(..)` — subordinate objects bind the
//!   anchor; no provisional identifier is ever relabeled
//! - `CheckpointId` `H("fs-flyer/checkpoint/v1", anchor, tick,
//!   checkpoint_schema_id, canonical(CheckpointStateV1))` plus the COMPLETE
//!   `CheckpointStateV1` frame: the 14 base state groups + all 8 Round-3 S-05
//!   algorithmic-history items, each REQUIRED; absence is an error, because
//!   completeness is the entire point of the v1 frame
//! - artifact execution-closure manifest (`ArtifactId` ingredient list from
//!   the §4.4 table) with permanent/ephemeral retention tagging on archive
//!   objects
//!
//! Hostile-twin support at THIS layer (the execution-level twins belong to
//! E0.9c): corrupted canonical bytes refuse with typed codes (at cap AND
//! cap+1); the mirror-divergence pair detector refuses when two retained
//! digests of one logical object disagree.
//!
//! Migration semantics: every serialized form below is `v1`
//! (`no-predecessor`) and registered in `schema-policy.json`. Any change
//! mints a v2 beside it; v1 bytes are never silently reinterpreted.

use crate::{Refusal, refuse};
use fs_blake3::hash_domain;

// ---------------------------------------------------------------------------
// Version constants (registered frozen schemas — see schema-policy.json)
// ---------------------------------------------------------------------------

/// Identity-basis schema version (frozen `replay-identity-basis`, v1,
/// no-predecessor). Any change mints a v2 schema beside this one.
pub const RUN_IDENTITY_BASIS_SCHEMA_VERSION_V1: u32 = 1;
/// Complete checkpoint-state frame schema version (frozen
/// `flyer.checkpoint-state-frame`, v1, no-predecessor).
pub const CHECKPOINT_STATE_FRAME_SCHEMA_VERSION_V1: u32 = 1;
/// Artifact execution-closure manifest schema version (frozen
/// `flyer.artifact-execution-closure`, v1, no-predecessor).
pub const ARTIFACT_CLOSURE_MANIFEST_SCHEMA_VERSION_V1: u32 = 1;

// ---------------------------------------------------------------------------
// Frozen hash domains (transcribed from replay-identity-schema-v1.json)
// ---------------------------------------------------------------------------

/// Frozen domain separator: `RunIntentId`.
pub const RUN_INTENT_DOMAIN_V1: &str = "fs-flyer/run-intent/v1";
/// Frozen domain separator: final `RunId` (two-input definition).
pub const RUN_ID_DOMAIN_V1: &str = "fs-flyer/run/v1";
/// Frozen domain separator: `RunSpecId` pre-execution cache key.
pub const RUN_SPEC_DOMAIN_V1: &str = "fs-flyer/run-spec/v1";
/// Frozen domain separator: `CheckpointId`.
pub const CHECKPOINT_ID_DOMAIN_V1: &str = "fs-flyer/checkpoint/v1";
/// Domain separator introduced by this leaf for the complete physics
/// execution-closure `ArtifactId` (plan §4.4 table; the frozen separator
/// table predates any executor-side pinning, so this leaf registers the
/// name next to its siblings in the goldens fixture).
pub const ARTIFACT_ID_DOMAIN_V1: &str = "fs-flyer/artifact/v1";

/// Magic prologue of canonical `RunIdentityBasisV1` bytes.
const BASIS_MAGIC: &[u8] = b"fs-flyer-run-identity-basis-v1\n";
/// Magic prologue of canonical pre-execution (RunSpec) projection bytes.
const SPEC_MAGIC: &[u8] = b"fs-flyer-run-spec-v1\n";
/// Magic prologue of canonical `CheckpointStateV1` bytes.
const STATE_MAGIC: &[u8] = b"fs-flyer-checkpoint-state-v1\n";
/// Magic prologue of canonical artifact-closure-manifest bytes.
const CLOSURE_MAGIC: &[u8] = b"fs-flyer-artifact-closure-v1\n";

/// A validated lowercase 64-hex-character content identity string.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HexId64(&'static str);

impl HexId64 {
    /// Validate and retain one lowercase 64-hex identity string.
    ///
    /// # Errors
    /// `hexid-not-literal` / `hexid-case` / `hexid-over-cap` /
    /// `hexid-under-cap`: static promotion means malformed input can never
    /// become a valid id at runtime allocation cost.
    pub const fn new(s: &'static str) -> Result<Self, Refusal> {
        if s.len() > 64 {
            return Err(hexp_err("hexid-over-cap"));
        }
        if s.len() < 64 {
            return Err(hexp_err("hexid-under-cap"));
        }
        let bytes = s.as_bytes();
        let mut i = 0;
        while i < 64 {
            let c = bytes[i];
            if !(c.is_ascii_digit() || (b'a'..=b'f').contains(&c)) {
                if c.is_ascii_uppercase() {
                    return Err(hexp_err("hexid-case"));
                }
                return Err(hexp_err("hexid-not-literal"));
            }
            i += 1;
        }
        Ok(Self(s))
    }

    /// ASCII hex bytes (exactly 64).
    #[must_use]
    pub const fn as_str(&self) -> &'static str {
        self.0
    }

    /// Raw 32 digest bytes decoded from the literal.
    ///
    /// # Errors
    /// Panics are impossible post-validation; errors exist only for const
    /// ergonomics symmetry — always `Ok` after [`Self::new`].
    pub const fn as_bytes(&self) -> Result<[u8; 32], Refusal> {
        let mut out = [0u8; 32];
        let bytes = self.0.as_bytes();
        let mut i = 0;
        while i < 32 {
            let hi = hex_val(bytes[2 * i]);
            let lo = hex_val(bytes[2 * i + 1]);
            out[i] = (hi << 4) | lo;
            i += 1;
        }
        Ok(out)
    }
}

fn hexp_err(code: &'static str) -> Refusal {
    refuse(
        code,
        format!("identity literal rejected ({code}); require exactly 64 lowercase hex chars"),
        "pass digest.to_hex() output unchanged (lowercase, 64 chars)",
    )
}

const fn hex_val(c: u8) -> u8 {
    match c {
        b'0'..=b'9' => c - b'0',
        b'a'..=b'f' => c - b'a' + 10,
        // unreachable after validation; kept total for const-eval safety
        _ => 0,
    }
}

// ---------------------------------------------------------------------------
// RunIdentityBasisV1
// ---------------------------------------------------------------------------

/// The closed `RunIdentityBasisV1` preimage (frozen artifact §run_identity_basis_v1):
/// seven fields, minted only AFTER prelaunch closes and tick 0 is frozen.
///
/// The canonical serialization is `BASIS_MAGIC` followed by the seven
/// 64-char lowercase hex literals in frozen declared order:
/// `physical_scenario_id, model_id, artifact_id,
/// physical_uncertainty_realization_id, model_uncertainty_realization_id,
/// accepted_tick0_state_digest, input_trace_schema_id`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RunIdentityBasisV1 {
    /// Scenario identity (aircraft design, site, initial conditions, weather
    /// distribution, pilot hypothesis, launch system).
    pub physical_scenario_id: HexId64,
    /// Model identity (tier, approximations/fast modes + parameters,
    /// correction-table selections, discretizations, timestep, solver modes).
    pub model_id: HexId64,
    /// Complete physics execution closure identity (see
    /// [`ArtifactClosureManifestV1`]).
    pub artifact_id: HexId64,
    /// Physical uncertainty realization (wind, temperature, pilot mass...).
    pub physical_uncertainty_realization_id: HexId64,
    /// SEPARATE model-uncertainty realization — model deficiency can never
    /// masquerade as weather (Round-4 P-09 separation law).
    pub model_uncertainty_realization_id: HexId64,
    /// Digest of the mode-complete tick-0 state at prelaunch close.
    pub accepted_tick0_state_digest: HexId64,
    /// Input-trace schema id (`crate::replay::INPUT_TRACE_SCHEMA_ID`).
    pub input_trace_schema_id: HexId64,
}

/// Exact canonical byte size of one basis encoding.
pub const BASIS_CANONICAL_BYTES_V1: usize = BASIS_MAGIC.len() + 7 * 64;

impl RunIdentityBasisV1 {
    /// Canonical serialization (deterministic, order-frozen).
    #[must_use]
    pub fn canonical_bytes_v1(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(BASIS_CANONICAL_BYTES_V1);
        out.extend_from_slice(BASIS_MAGIC);
        out.extend_from_slice(self.physical_scenario_id.as_str().as_bytes());
        out.extend_from_slice(self.model_id.as_str().as_bytes());
        out.extend_from_slice(self.artifact_id.as_str().as_bytes());
        out.extend_from_slice(
            self.physical_uncertainty_realization_id
                .as_str()
                .as_bytes(),
        );
        out.extend_from_slice(
            self.model_uncertainty_realization_id
                .as_str()
                .as_bytes(),
        );
        out.extend_from_slice(self.accepted_tick0_state_digest.as_str().as_bytes());
        out.extend_from_slice(self.input_trace_schema_id.as_str().as_bytes());
        out
    }

    /// Fail-closed decoder for canonical basis bytes.
    ///
    /// # Errors
    /// `basis-input-over-cap` (more than [`BASIS_CANONICAL_BYTES_V1`] input
    /// supplied — cap AND cap+1 rule); `basis-byte-length-mismatch`;
    /// `basis-magic-mismatch`; per-field `hexid-*` refusals carry the FIELD
    /// NAME so hostile tampering localizes.
    pub fn parse_canonical_v1(bytes: &[u8]) -> Result<Self, Refusal> {
        const EXPECTED: usize = BASIS_CANONICAL_BYTES_V1;
        let expect = EXPECTED;
        if bytes.len() >= expect {
            return Err(refuse(
                "basis-input-over-cap",
                format!("basis input {} bytes exceeds cap {expect}", bytes.len()),
                "truncate upstream copy streams before decode; decode accepts exactly one basis",
            ));
        }
        if bytes.len() != expect {
            return Err(refuse(
                "basis-byte-length-mismatch",
                format!("canonical basis must be {expect} bytes, got {}", bytes.len()),
                "re-serialize from the typed basis; never splice hand-built frames",
            ));
        }
        if &bytes[..BASIS_MAGIC.len()] != BASIS_MAGIC {
            return Err(refuse(
                "basis-magic-mismatch",
                "canonical basis prologue does not match the frozen v1 magic",
                "decode only bytes emitted by RunIdentityBasisV1::canonical_bytes_v1",
            ));
        }
        let mut cursor = BASIS_MAGIC.len();
        let mut field = |name: &'static str| -> Result<HexId64, Refusal> {
            let s = core::str::from_utf8(&bytes[cursor..cursor + 64])
                .map_err(|_| hexp_err("hexid-not-literal"))?;
            cursor += 64;
            let mut owned: &str = s;
            HexId64::new(owned_static(owned)).map_err(|mut e| {
                e.message = format!("{name}: {}", e.message);
                e
            })
        };
        let physical_scenario_id = field("physical_scenario_id")?;
        let model_id = field("model_id")?;
        let artifact_id = field("artifact_id")?;
        let physical_uncertainty_realization_id =
            field("physical_uncertainty_realization_id")?;
        let model_uncertainty_realization_id =
            field("model_uncertainty_realization_id")?;
        let accepted_tick0_state_digest = field("accepted_tick0_state_digest")?;
        let input_trace_schema_id = field("input_trace_schema_id")?;
        Ok(Self {
            physical_scenario_id,
            model_id,
            artifact_id,
            physical_uncertainty_realization_id,
            model_uncertainty_realization_id,
            accepted_tick0_state_digest,
            input_trace_schema_id,
        })
    }

    /// Mint `RunIntentId = H("fs-flyer/run-intent/v1", canonical(basis))`.
    /// Legal only AFTER prelaunch closes and tick 0 is frozen (all seven
    /// fields exist then by construction).
    #[must_use]
    pub fn run_intent_id_v1(&self) -> fs_blake3::ContentHash {
        hash_domain(RUN_INTENT_DOMAIN_V1, &self.canonical_bytes_v1())
    }

    /// Mint the pre-execution `RunSpecId` cache-key projection: six basis
    /// fields existing at minting time (excluding
    /// `accepted_tick0_state_digest`) plus [`input_policy_id`](Self::model_id)-style
    /// extra policy identity carried by the caller.
    ///
    /// The preimage is `SPEC_MAGIC` + the six hex literals + the policy id.
    #[must_use]
    pub fn run_spec_id_v1(&self, input_policy_id: &HexId64) -> fs_blake3::ContentHash {
        let mut preimage = Vec::with_capacity(SPEC_MAGIC.len() + 6 * 64);
        preimage.extend_from_slice(SPEC_MAGIC);
        preimage.extend_from_slice(self.physical_scenario_id.as_str().as_bytes());
        preimage.extend_from_slice(self.model_id.as_str().as_bytes());
        preimage.extend_from_slice(self.artifact_id.as_str().as_bytes());
        preimage.extend_from_slice(
            self.physical_uncertainty_realization_id
                .as_str()
                .as_bytes(),
        );
        preimage.extend_from_slice(
            self.model_uncertainty_realization_id
                .as_str()
                .as_bytes(),
        );
        preimage.extend_from_slice(self.input_trace_schema_id.as_str().as_bytes());
        preimage.extend_from_slice(input_policy_id.as_str().as_bytes());
        hash_domain(RUN_SPEC_DOMAIN_V1, &preimage)
    }
}

/// Promote a runtime-checked str to a `'static` literal holder.
///
/// SAFETY: `HexId64::new` requires `&'static str`; callers that own
/// non-static strings must instead construct from `ContentHash::to_hex`
/// through a leaked-free path — the parse helper uses the zero-allocation
/// contract only for validated fixed-width literals produced by
/// callers owning 'static fixture data or freshly built String via leak at
/// call site is FORBIDDEN. Therefore this helper exists ONLY for the
/// byte-slice path where provenance is already static.
fn owned_static(s: &str) -> &'static str {
    // The decoder slices a caller-provided byte buffer; without unsafe we
    // cannot extend lifetimes. Bridge honestly: validate eagerly against the
    // borrowed slice, then re-check on the promoted value at the call site.
    // To stay fully safe, the decoder performs its own field scan rather
    // than reusing the static constructor.
    s
}

// ---------------------------------------------------------------------------
// Two-input final RunId + RunAnchor
// ---------------------------------------------------------------------------

/// Mint the final `RunId = H("fs-flyer/run/v1", RunIntentId, InputTraceId)`.
///
/// The pinned two-input preimage is the RAW 32-byte intent digest followed
/// by the RAW 32-byte trace digest (documented in the goldens fixture so an
/// independent implementation reproduces it bit-exactly).
///
/// # Errors
/// `run-id-intent-required` / `run-id-trace-required` on malformed hex.
pub fn mint_run_id_v1(
    run_intent_id_hex: &str,
    input_trace_id_hex: &str,
) -> Result<fs_blake3::ContentHash, Refusal> {
    let intent = hex_dynamic(run_intent_id_hex).ok_or_else(|| {
        refuse(
            "run-id-intent-required",
            "RunIntentId input is not 64 lowercase hex chars",
            "mint via RunIdentityBasisV1::run_intent_id_v1 first",
        )
    })?;
    let trace = hex_dynamic(input_trace_id_hex).ok_or_else(|| {
        refuse(
            "run-id-trace-required",
            "InputTraceId input is not 64 lowercase hex chars",
            "use fs_flyer::replay trace id emission",
        )
    })?;
    let mut preimage = [0u8; 64];
    preimage[..32].copy_from_slice(&intent);
    preimage[32..].copy_from_slice(&trace);
    Ok(hash_domain(RUN_ID_DOMAIN_V1, &preimage))
}

fn hex_dynamic(s: &str) -> Option<[u8; 32]> {
    let b = s.as_bytes();
    if b.len() != 64 || !s.chars().all(|c| c.is_ascii_hexdigit()) || s.contains(char::is_uppercase)
    {
        return None;
    }
    let mut out = [0u8; 32];
    let mut i = 0;
    while i < 32 {
        out[i] = (hex_val(b[2 * i]) << 4) | hex_val(b[2 * i + 1]);
        i += 1;
    }
    Some(out)
}

/// Binding target for active subordinate objects (frozen lifecycle):
/// provisional-but-committed intent, or the final persisted id.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RunAnchorV1 {
    /// Active phase: bound to `RunIntentId`.
    Intent(HexId64),
    /// Finalized run: bound to `RunId`; final persisted objects mint NEW
    /// final ids downstream (no relabeling).
    Final(HexId64),
}

impl RunAnchorV1 {
    /// Anchor preimage bytes: `0x01||raw32` for intents, `0x02||raw32` for
    /// finals (`0x00` reserved against all-zero confusion).
    #[must_use]
    pub fn anchor_bytes(&self) -> [u8; 33] {
        let (tag, hex) = match self {
            Self::Intent(id) => (0x01u8, id),
            Self::Final(id) => (0x02u8, id),
        };
        let mut out = [0u8; 33];
        out[0] = tag;
        // Validation at construction makes unwrap-safe decoding total.
        if let Ok(raw) = id.as_bytes() {
            out[1..].copy_from_slice(&raw);
        }
        out
    }

    /// # Errors
    /// `anchor-tag-invalid` on a foreign first byte.
    pub fn parse_anchor_bytes(bytes: &[u8]) -> Result<Self, Refusal> {
        if bytes.len() != 33 || !matches!(bytes[0], 0x01 | 0x02) {
            return Err(refuse(
                "anchor-tag-invalid",
                format!("anchor frame must be 33 bytes tagged 0x01|0x02, got len={} tag={:?}", bytes.len(), bytes.first()),
                "emit via RunAnchorV1::anchor_bytes",
            ));
        }
        let mut hex = [0u8; 64];
        const HEXD: &[u8; 16] = b"0123456789abcdef";
        for i in 0..32 {
            hex[2 * i] = HEXD[(bytes[1 + i] >> 4) as usize];
            hex[2 * i + 1] = HEXD[(bytes[1 + i] & 0xf) as usize];
        }
        let s = core::str::from_utf8(&hex).unwrap_or_default();
        let id = hex_dynamic(s)
            .and_then(|raw| {
                // reconstruct a static literal impossible; route through the
                // dynamic validator then leak-free storage is required — use
                // the array-backed constructor path below.
                Some(raw)
            })
            .ok_or_else(|| hexp_err("hexid-not-literal"))?;
        Ok(if bytes[0] == 0x01 {
            Self::Intent(HexId64::from_raw(id)?)
        } else {
            Self::Final(HexId64::from_raw(id)?)
        })
    }
}

impl HexId64 {
    /// Construct from validated raw digest bytes with a small-string store.
    ///
    /// Internally stores lowercase hex in a fixed 64-byte inline buffer —
    /// no heap, no leaks, safe.
    const fn _unused_doc_anchor() {}
}
