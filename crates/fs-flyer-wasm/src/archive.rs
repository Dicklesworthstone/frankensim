//! Archive-fixture loader (bead wf-root-guzez.1.9.3, E0.9c vendor-independent
//! slice). Implements the verifier-mediated loading discipline of the frozen
//! E0.9a contract against a LOCAL content-addressed store: every archived byte
//! source is verified (size, then BLAKE3 content identity) before parsing,
//! dual-publication read-back checks both copies independently, and a
//! backward-playback replay re-executes an archived generation and compares
//! trajectory digests. The parser is STRICT fail-closed — exact key order,
//! no tolerated deviations — because a tolerant line reader is exactly how a
//! gate goes silently dead (workspace incident on golden-couplings.json).
//!
//! Vendor swap (R2/S3, TUF metadata, real trust roots) is E0.9b; this module
//! is the loader mechanics those vendors plug into.

use crate::{Refusal, hello_digest};
use fs_blake3::hash_bytes;

/// One entry of the targets manifest: path, exact size, content identity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ArchiveTarget {
    /// Store-relative path of the byte source.
    pub path: String,
    /// Exact byte length (checked BEFORE hashing).
    pub size_bytes: u64,
    /// Lowercase BLAKE3 hex of the exact bytes.
    pub blake3_hex: String,
}

/// Verify archived bytes against a target: size first, then content identity.
///
/// # Errors
/// `archive-size-mismatch` / `archive-content-digest-mismatch`.
pub fn verify_target_bytes(target: &ArchiveTarget, bytes: &[u8]) -> Result<(), Refusal> {
    if bytes.len() as u64 != target.size_bytes {
        return Err(Refusal {
            code: "archive-size-mismatch",
            message: format!(
                "{}: {} bytes on read-back, manifest says {}",
                target.path,
                bytes.len(),
                target.size_bytes
            ),
            ranked_repairs: vec![
                "re-fetch the object; a truncated read is the common cause".into(),
                "if the store object truly changed, the archive is corrupt — restore from the mirror".into(),
            ],
        });
    }
    let got = hash_bytes(bytes).to_hex();
    if got != target.blake3_hex {
        return Err(Refusal {
            code: "archive-content-digest-mismatch",
            message: format!("{}: BLAKE3 {} != manifest {}", target.path, got, target.blake3_hex),
            ranked_repairs: vec![
                "restore the object from the WORM mirror".into(),
                "if this is a NEW generation, publish it under a NEW path — paths are content-addressed and immutable".into(),
            ],
        });
    }
    Ok(())
}

/// Dual-publication read-back verification: BOTH copies must independently
/// satisfy the manifest AND be byte-identical (the E0.9a publication rule —
/// never provider-native replication).
///
/// # Errors
/// Per-copy refusals from [`verify_target_bytes`], or
/// `archive-mirror-divergence` when both verify individually yet differ
/// (impossible under an honest manifest, kept as a fail-closed tripwire).
pub fn verify_dual_publication(
    target: &ArchiveTarget,
    primary: &[u8],
    mirror: &[u8],
) -> Result<(), Refusal> {
    verify_target_bytes(target, primary)?;
    verify_target_bytes(target, mirror)?;
    if primary != mirror {
        return Err(Refusal {
            code: "archive-mirror-divergence",
            message: format!("{}: primary and mirror verify yet differ", target.path),
            ranked_repairs: vec!["treat the store as compromised; audit both providers".into()],
        });
    }
    Ok(())
}

/// A caller-supplied immutable archive endpoint.
///
/// Implementors own provider authentication, conditional-create semantics,
/// retention controls, and any provider-specific receipt. The coordinator
/// never fabricates those facts: it only orders independent writes and checks
/// the bytes it reads back.
pub trait ArchivePublicationStore {
    /// Stable identity of the independently administered provider endpoint.
    ///
    /// This is a configuration guard, not an authentication claim. A real
    /// deployment must bind it to provider-side authority outside this crate.
    fn provider_id(&self) -> &str;

    /// Create this content-addressed object only when the endpoint admits an
    /// immutable publication.
    fn put_immutable(&mut self, target: &ArchiveTarget, bytes: &[u8]) -> Result<(), Refusal>;

    /// Read the published object back through the endpoint's normal retrieval
    /// path. The coordinator verifies the returned bytes before proceeding.
    fn read_back(&mut self, target: &ArchiveTarget) -> Result<Vec<u8>, Refusal>;
}

/// Caller-supplied authority for publishing TUF targets metadata.
///
/// Its implementation owns role keys, threshold policy, signatures, and
/// rollback/freeze protection. It is invoked only after both archive copies
/// independently pass read-back verification.
pub trait TufTargetsPublisher {
    /// Publish metadata that names this already verified target.
    fn publish_target(&mut self, target: &ArchiveTarget) -> Result<(), Refusal>;
}

fn invalid_publication_topology(detail: String) -> Refusal {
    Refusal {
        code: "archive-publication-topology-invalid",
        message: detail,
        ranked_repairs: vec![
            "configure distinct, independently administered primary and mirror endpoints".into(),
            "bind each endpoint identity to its provider-side authorization policy".into(),
        ],
    }
}

/// Publish to two independently identified endpoints, read both copies back,
/// verify them against the immutable target, then allow TUF target publication.
///
/// The input is verified before any write. Target metadata therefore cannot be
/// published after a failed upload, read-back, or byte verification. This is
/// deliberately transport- and key-agnostic: it does not create cloud
/// resources, select retention controls, sign metadata, or claim that an
/// endpoint's reported identity is authentic.
///
/// # Errors
/// Returns input/read-back verification errors, endpoint refusals, or
/// `archive-publication-topology-invalid` when endpoints are blank or not
/// independently identified. TUF metadata publication is attempted only after
/// both copies verify.
pub fn publish_dual_with_verification(
    target: &ArchiveTarget,
    bytes: &[u8],
    primary: &mut impl ArchivePublicationStore,
    mirror: &mut impl ArchivePublicationStore,
    tuf_targets: &mut impl TufTargetsPublisher,
) -> Result<(), Refusal> {
    verify_target_bytes(target, bytes)?;

    let primary_id = primary.provider_id();
    let mirror_id = mirror.provider_id();
    if primary_id.is_empty() || mirror_id.is_empty() || primary_id == mirror_id {
        return Err(invalid_publication_topology(format!(
            "primary {:?} and mirror {:?} must be distinct non-empty provider identities",
            primary_id, mirror_id
        )));
    }

    primary.put_immutable(target, bytes)?;
    let primary_bytes = primary.read_back(target)?;

    mirror.put_immutable(target, bytes)?;
    let mirror_bytes = mirror.read_back(target)?;

    verify_dual_publication(target, &primary_bytes, &mirror_bytes)?;
    tuf_targets.publish_target(target)
}

/// The archived hello generation: exact scenario parameters plus the pinned
/// trajectory digest. dt is an INTEGER RATIO (plan integer-ratio doctrine) so
/// the envelope never round-trips a float through text formatting.
#[derive(Clone, Debug, PartialEq)]
pub struct HelloArchiveEnvelope {
    /// Generation ordinal within the fixture store.
    pub generation: u32,
    /// Principal inertia triple [kg·m²].
    pub inertia_kg_m2: [f64; 3],
    /// Initial unit quaternion (w, x, y, z).
    pub q0: [f64; 4],
    /// Initial body angular velocity [rad/s].
    pub omega0: [f64; 3],
    /// Timestep numerator [s].
    pub dt_num: u32,
    /// Timestep denominator.
    pub dt_den: u32,
    /// Integration steps.
    pub steps: u32,
    /// Pinned trajectory digest (the archived generation's identity).
    pub trajectory_digest_hex: String,
}

/// Envelope schema id — line 1 of every v1 envelope.
pub const HELLO_ENVELOPE_SCHEMA: &str = "org.frankensim.wright-flyer.hello-archive-envelope.v1";

const ENVELOPE_KEYS: [&str; 9] = [
    "schema",
    "generation",
    "inertia",
    "q0",
    "omega0",
    "dt_num",
    "dt_den",
    "steps",
    "trajectory_digest",
];

fn malformed(detail: String) -> Refusal {
    Refusal {
        code: "archive-envelope-malformed",
        message: detail,
        ranked_repairs: vec![
            "regenerate the envelope with the canonical writer — hand edits are not admitted"
                .into(),
        ],
    }
}

fn parse_f64_list<const N: usize>(key: &str, raw: &str) -> Result<[f64; N], Refusal> {
    let parts: Vec<&str> = raw.split(',').collect();
    if parts.len() != N {
        return Err(malformed(format!(
            "{key}: expected {N} comma-separated values"
        )));
    }
    let mut out = [0.0; N];
    for (i, p) in parts.iter().enumerate() {
        out[i] = p
            .parse::<f64>()
            .map_err(|e| malformed(format!("{key}[{i}] = {p:?}: {e}")))?;
        if !out[i].is_finite() {
            return Err(malformed(format!("{key}[{i}] is not finite")));
        }
    }
    Ok(out)
}

/// Parse a v1 envelope. STRICT: exactly the nine keys, in canonical order,
/// `key=value` per line, nothing else. Any deviation is a typed refusal —
/// tolerance here is how integrity gates die silently.
///
/// # Errors
/// `archive-envelope-malformed` with the exact offending line.
pub fn parse_hello_envelope(text: &str) -> Result<HelloArchiveEnvelope, Refusal> {
    let lines: Vec<&str> = text.trim_end_matches('\n').split('\n').collect();
    if lines.len() != ENVELOPE_KEYS.len() {
        return Err(malformed(format!(
            "expected exactly {} lines, got {}",
            ENVELOPE_KEYS.len(),
            lines.len()
        )));
    }
    let mut values: Vec<&str> = Vec::with_capacity(ENVELOPE_KEYS.len());
    for (i, (line, key)) in lines.iter().zip(ENVELOPE_KEYS.iter()).enumerate() {
        let Some((k, v)) = line.split_once('=') else {
            return Err(malformed(format!("line {}: missing '='", i + 1)));
        };
        if k != *key {
            return Err(malformed(format!(
                "line {}: key {k:?} out of canonical order (expected {key:?})",
                i + 1
            )));
        }
        values.push(v);
    }
    if values[0] != HELLO_ENVELOPE_SCHEMA {
        return Err(malformed(format!("unknown schema {:?}", values[0])));
    }
    let parse_u32 = |key: &str, raw: &str| -> Result<u32, Refusal> {
        raw.parse::<u32>()
            .map_err(|e| malformed(format!("{key} = {raw:?}: {e}")))
    };
    let env = HelloArchiveEnvelope {
        generation: parse_u32("generation", values[1])?,
        inertia_kg_m2: parse_f64_list::<3>("inertia", values[2])?,
        q0: parse_f64_list::<4>("q0", values[3])?,
        omega0: parse_f64_list::<3>("omega0", values[4])?,
        dt_num: parse_u32("dt_num", values[5])?,
        dt_den: parse_u32("dt_den", values[6])?,
        steps: parse_u32("steps", values[7])?,
        trajectory_digest_hex: values[8].to_string(),
    };
    if env.dt_den == 0 {
        return Err(malformed("dt_den must be positive".into()));
    }
    Ok(env)
}

/// Backward playback: re-execute the archived generation on the CURRENT
/// kernel and compare trajectory digests. Equality is the "old-exact"
/// receipt; divergence is a typed refusal naming both digests.
///
/// # Errors
/// Kernel refusals pass through; `archive-replay-digest-mismatch` when the
/// replayed trajectory does not reproduce the archived identity.
pub fn replay_generation(env: &HelloArchiveEnvelope) -> Result<String, Refusal> {
    let dt_s = f64::from(env.dt_num) / f64::from(env.dt_den);
    let got = hello_digest(env.inertia_kg_m2, env.q0, env.omega0, dt_s, env.steps)?;
    if got != env.trajectory_digest_hex {
        return Err(Refusal {
            code: "archive-replay-digest-mismatch",
            message: format!(
                "generation {}: replayed digest {} != archived {}",
                env.generation, got, env.trajectory_digest_hex
            ),
            ranked_repairs: vec![
                "the current kernel diverged from the archived physics — old-exact playback requires the ARCHIVED artifact, not the current one".into(),
                "if the kernel change is intentional, mint a new generation under the golden-bump protocol".into(),
            ],
        });
    }
    Ok(got)
}
