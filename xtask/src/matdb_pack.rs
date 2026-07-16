//! Strict offline compiler for canonical fs-matdb scalar/curve packs.
//!
//! This tool boundary owns source-format and redistribution decisions. The
//! runtime fs-matdb crate only admits the normalized artifact. Schema v1 is
//! deliberately honest: it supports scalar/curve material claims and typed
//! joint statistics, while species, NASA-9, and kinetics profiles refuse until
//! their runtime codecs exist.
//!
//! The bounded manifest is tab-separated: one pack_id, redistribution,
//! citation, and license record, followed by source records containing logical
//! id, safe relative path, and profile. A material-tsv-v1 source starts with
//! SOURCE_HEADER and then uses observation, scalar, curve, uncertainty,
//! validity, frame, and joint records. Every numeric token has a separate
//! explicit unit field; confidence and correlation use the exact basis 1.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File, OpenOptions};
use std::io::{Read as _, Write as _};
use std::path::{Component, Path, PathBuf};

#[cfg(unix)]
use std::os::unix::fs::MetadataExt as _;

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimId, ClaimSet, InterpolationPolicy, JointStatistics, MATDB_PACK_TARGET_BASIS,
    NormalizationReceipt, NormalizationTarget, NormalizedPack, ObservationDataset, ObservationId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, StatisticMember, UncertaintyModel,
    ValidityBoundSide,
};
use fs_qty::Dims;
use fs_qty::parse::{ParseBudget, parse_qty_with_budget};

const COMPILER_ID: &str = "frankensim-matdb-pack-compiler-v1";
const MANIFEST_HEADER: &str = "frankensim.matdb-manifest.v1";
const SOURCE_HEADER: &str = "frankensim.matdb-source.v1";
const MATERIAL_PROFILE: &str = "material-tsv-v1";
const SOURCE_ENVELOPE_DOMAIN: &str = "org.frankensim.xtask.matdb-pack.source-envelope.v1";
const SOURCE_FILE_DOMAIN: &str = "org.frankensim.xtask.matdb-pack.source-file.v1";
const SOURCE_RECORD_DOMAIN: &str = "org.frankensim.xtask.matdb-pack.source-record.v1";
const SOURCE_LITERAL_DOMAIN: &str = "org.frankensim.xtask.matdb-pack.source-literal.v1";
const DECISION_ID_DOMAIN: &str = "org.frankensim.xtask.matdb-pack.decision.v1";
const MAX_MANIFEST_BYTES: usize = 1_048_576;
const MAX_SOURCE_BYTES: usize = 64 * 1_048_576;
const MAX_TOTAL_SOURCE_BYTES: usize = 256 * 1_048_576;
const MAX_SOURCES: usize = 256;
const MAX_RECORDS: usize = 200_000;
const MAX_LINE_BYTES: usize = 1_048_576;
const MAX_IDENTIFIER_BYTES: usize = 128;
const MAX_TEXT_BYTES: usize = 1_048_576;
const MAX_CURVE_KNOTS: usize = 1_000_000;
const MAX_JOINT_MEMBERS: usize = 256;
const MAX_NORMALIZATION_RECEIPTS: usize = 100_000;
const MAX_REPEATED_PROVENANCE_BYTES: usize = 64 * 1_048_576;
const QTY_BUDGET: ParseBudget = ParseBudget::new(4_096, 256, 64, 256);

#[derive(Debug, Clone, PartialEq, Eq)]
struct CompileError {
    code: &'static str,
    subject: String,
    detail: String,
    input_hash: Option<ContentHash>,
}

impl CompileError {
    fn new(code: &'static str, subject: impl Into<String>, detail: impl Into<String>) -> Self {
        Self {
            code,
            subject: subject.into(),
            detail: detail.into(),
            input_hash: None,
        }
    }

    fn with_input_hash(mut self, input_hash: ContentHash) -> Self {
        if self.input_hash.is_none() {
            self.input_hash = Some(input_hash);
        }
        self
    }
}

impl std::fmt::Display for CompileError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "matdb pack refused [{}] {}: {}",
            self.code, self.subject, self.detail
        )
    }
}

impl std::error::Error for CompileError {}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Decision {
    subject: String,
    verdict: &'static str,
    reason_code: &'static str,
    detail: String,
    source_hash: Option<ContentHash>,
    pack_hash: Option<ContentHash>,
}

impl Decision {
    fn admit(
        subject: impl Into<String>,
        reason_code: &'static str,
        detail: impl Into<String>,
        source_hash: Option<ContentHash>,
    ) -> Self {
        Self {
            subject: subject.into(),
            verdict: "admit",
            reason_code,
            detail: detail.into(),
            source_hash,
            pack_hash: None,
        }
    }

    fn refusal(error: &CompileError, input_hash: Option<ContentHash>) -> Self {
        Self {
            subject: error.subject.clone(),
            verdict: "refuse",
            reason_code: error.code,
            detail: error.detail.clone(),
            source_hash: input_hash,
            pack_hash: None,
        }
    }

    fn refusal_for_pack(error: &CompileError, pack_hash: ContentHash) -> Self {
        let mut decision = Self::refusal(error, None);
        decision.pack_hash = Some(pack_hash);
        decision
    }

    fn canonical_preimage(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        for part in [
            COMPILER_ID.as_bytes(),
            self.subject.as_bytes(),
            self.verdict.as_bytes(),
            self.reason_code.as_bytes(),
            self.detail.as_bytes(),
        ] {
            push_part(&mut bytes, part);
        }
        push_part(&mut bytes, b"source_hash");
        push_part(&mut bytes, &[u8::from(self.source_hash.is_some())]);
        if let Some(hash) = self.source_hash {
            push_part(&mut bytes, hash.as_bytes());
        }
        push_part(&mut bytes, b"pack_hash");
        push_part(&mut bytes, &[u8::from(self.pack_hash.is_some())]);
        if let Some(hash) = self.pack_hash {
            push_part(&mut bytes, hash.as_bytes());
        }
        bytes
    }
}

#[derive(Debug)]
struct CompileOutput {
    pack: NormalizedPack,
    bytes: Vec<u8>,
    decisions: Vec<Decision>,
}

#[derive(Debug)]
struct Manifest {
    pack_id: String,
    redistribution_terms: String,
    citation: String,
    license: String,
    sources: Vec<SourceSpec>,
}

#[derive(Debug, Clone)]
struct SourceSpec {
    id: String,
    relative_path: PathBuf,
    profile: String,
}

#[derive(Debug)]
struct LoadedSource {
    spec: SourceSpec,
    bytes: Vec<u8>,
    text: String,
    hash: ContentHash,
}

#[derive(Debug)]
struct RawObservation {
    id: String,
    specimen: String,
    method: String,
    caveats: String,
    source_id: String,
    source_hash: ContentHash,
    record_hash: ContentHash,
}

#[derive(Debug)]
enum RawClaimValue {
    Scalar {
        number: String,
        unit: String,
    },
    Curve {
        abscissa: String,
        abscissa_unit: String,
        ordinate_unit: String,
        points: Vec<(String, String)>,
    },
}

#[derive(Debug)]
struct RawClaim {
    id: String,
    observations: Vec<String>,
    property: String,
    value: RawClaimValue,
    interpolation: InterpolationPolicy,
    source_id: String,
    source_hash: ContentHash,
    record_hash: ContentHash,
}

#[derive(Debug)]
enum RawUncertaintyKind {
    Unstated,
    Absolute {
        number: String,
        unit: String,
        confidence: String,
        confidence_unit: String,
    },
    Relative {
        number: String,
        unit: String,
        confidence: String,
        confidence_unit: String,
    },
}

#[derive(Debug)]
struct RawUncertainty {
    kind: RawUncertaintyKind,
    source_hash: ContentHash,
}

#[derive(Debug)]
struct RawValidity {
    axis: String,
    lower: String,
    upper: String,
    unit: String,
    source_hash: ContentHash,
}

#[derive(Debug, Clone)]
struct RawFrame {
    source: String,
    target: String,
    source_hash: ContentHash,
}

#[derive(Debug)]
struct RawJoint {
    observation: String,
    block_id: String,
    members: Vec<String>,
    covariance: Vec<String>,
    correlation: Option<Vec<String>>,
    source_hash: ContentHash,
}

#[derive(Debug, Default)]
struct RawDatabase {
    observations: BTreeMap<String, RawObservation>,
    claims: BTreeMap<String, RawClaim>,
    uncertainties: BTreeMap<String, RawUncertainty>,
    validities: BTreeMap<String, BTreeMap<String, RawValidity>>,
    frames: BTreeMap<String, RawFrame>,
    joints: BTreeMap<(String, String), RawJoint>,
    decisions: Vec<Decision>,
    records: usize,
    curve_knots: usize,
}

#[derive(Debug, Clone)]
struct ParsedQuantity {
    value: f64,
    dims: Dims,
    scale: f64,
    offset: f64,
    unit: String,
    literal_hash: ContentHash,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum LocalComponent {
    Scalar,
    CurveAbscissa(u32),
    CurveOrdinate(u32),
}

#[derive(Debug, Clone)]
struct ComponentMeta {
    member: StatisticMember,
    dims: Dims,
    scale: f64,
    unit: String,
    observations: BTreeSet<ObservationId>,
}

/// Compile a strict source manifest twice and publish one verified new artifact.
pub(super) fn cmd_matdb_pack(root: &Path, args: &[String]) -> Result<(), String> {
    let (manifest_path, output_path) = match parse_cli(root, args) {
        Ok(paths) => paths,
        Err(error) => {
            println!(
                "{}",
                render_decision(&Decision::refusal(&error, Some(cli_input_hash(args))))
            );
            return Err(error.to_string());
        }
    };
    let mut first = match compile_manifest(&manifest_path) {
        Ok(output) => output,
        Err(error) => {
            println!(
                "{}",
                render_decision(&Decision::refusal(&error, error.input_hash,))
            );
            return Err(error.to_string());
        }
    };
    let second = match compile_manifest(&manifest_path) {
        Ok(output) => output,
        Err(error) => {
            println!(
                "{}",
                render_decision(&Decision::refusal(&error, error.input_hash,))
            );
            return Err(format!("second deterministic pass failed: {error}"));
        }
    };
    let first_source = first.pack.source_artifact();
    let second_source = second.pack.source_artifact();
    if first_source != second_source {
        let mut pair = Vec::new();
        push_part(&mut pair, b"first_source_artifact");
        push_part(&mut pair, first_source.as_bytes());
        push_part(&mut pair, b"second_source_artifact");
        push_part(&mut pair, second_source.as_bytes());
        let error = CompileError::new(
            "input_changed_between_passes",
            "sources",
            "source envelope changed between the two complete compilation passes",
        )
        .with_input_hash(hash_domain(SOURCE_ENVELOPE_DOMAIN, &pair));
        println!(
            "{}",
            render_decision(&Decision::refusal(&error, error.input_hash))
        );
        return Err(error.to_string());
    }
    if first.bytes != second.bytes || first.decisions != second.decisions {
        let error = CompileError::new(
            "nondeterministic_compilation",
            "pack",
            "two complete source-to-pack passes produced different artifacts or decisions",
        )
        .with_input_hash(first_source);
        println!(
            "{}",
            render_decision(&Decision::refusal(&error, error.input_hash,))
        );
        return Err(error.to_string());
    }

    let pack_hash = first.pack.content_hash();
    write_new_verified(&output_path, &first.bytes, pack_hash).map_err(|error| {
        println!(
            "{}",
            render_decision(&Decision::refusal_for_pack(&error, pack_hash))
        );
        error.to_string()
    })?;
    first.decisions.push(Decision {
        subject: "output".to_string(),
        verdict: "admit",
        reason_code: "published_new_verified_artifact",
        detail: "created one previously absent output and re-admitted its exact bytes".to_string(),
        source_hash: None,
        pack_hash: Some(pack_hash),
    });
    sort_decisions(&mut first.decisions);
    for decision in &first.decisions {
        println!("{}", render_decision(decision));
    }
    eprintln!(
        "matdb pack OK: {} bytes, hash {}, {} decisions",
        first.bytes.len(),
        pack_hash,
        first.decisions.len()
    );
    Ok(())
}

fn parse_cli(root: &Path, args: &[String]) -> Result<(PathBuf, PathBuf), CompileError> {
    let mut manifest = None;
    let mut output = None;
    let mut index = 0;
    while index < args.len() {
        let option = &args[index];
        let value = args.get(index + 1).ok_or_else(|| {
            CompileError::new(
                "invalid_cli",
                "arguments",
                format!("option {option:?} requires one path"),
            )
        })?;
        let slot = match option.as_str() {
            "--manifest" => &mut manifest,
            "--out" => &mut output,
            _ => {
                return Err(CompileError::new(
                    "invalid_cli",
                    "arguments",
                    format!("unknown option {option:?}; expected --manifest and --out"),
                ));
            }
        };
        if slot.replace(resolve_cli_path(root, value)).is_some() {
            return Err(CompileError::new(
                "invalid_cli",
                "arguments",
                format!("option {option:?} was supplied more than once"),
            ));
        }
        index += 2;
    }
    let manifest = manifest.ok_or_else(|| {
        CompileError::new(
            "invalid_cli",
            "arguments",
            "missing required --manifest path",
        )
    })?;
    let output = output.ok_or_else(|| {
        CompileError::new("invalid_cli", "arguments", "missing required --out path")
    })?;
    if manifest == output {
        return Err(CompileError::new(
            "invalid_cli",
            "arguments",
            "manifest and output paths must differ",
        ));
    }
    Ok((manifest, output))
}

fn resolve_cli_path(root: &Path, raw: &str) -> PathBuf {
    let path = PathBuf::from(raw);
    if path.is_absolute() {
        path
    } else {
        root.join(path)
    }
}

fn cli_input_hash(args: &[String]) -> ContentHash {
    let mut payload = Vec::new();
    for argument in args {
        push_part(&mut payload, argument.as_bytes());
    }
    hash_domain(SOURCE_ENVELOPE_DOMAIN, &payload)
}

fn parse_manifest(text: &str) -> Result<Manifest, CompileError> {
    if text.as_bytes().contains(&0) {
        return Err(CompileError::new(
            "invalid_manifest_encoding",
            "manifest",
            "NUL bytes are forbidden",
        ));
    }
    let mut lines = text.lines();
    if lines.next() != Some(MANIFEST_HEADER) {
        return Err(CompileError::new(
            "unsupported_manifest_schema",
            "manifest",
            format!("first line must be {MANIFEST_HEADER:?}"),
        ));
    }
    let mut pack_id = None;
    let mut redistribution_terms = None;
    let mut citation = None;
    let mut license = None;
    let mut sources = Vec::new();
    for (offset, line) in lines.enumerate() {
        let line_number = offset + 2;
        if line.len() > MAX_LINE_BYTES {
            return Err(CompileError::new(
                "resource_limit",
                format!("manifest:line:{line_number}"),
                "line exceeds the public byte budget",
            ));
        }
        if line.is_empty() {
            return Err(CompileError::new(
                "invalid_manifest_record",
                format!("manifest:line:{line_number}"),
                "blank records are forbidden",
            ));
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields.first().copied() {
            Some("pack_id") => {
                require_field_count(&fields, 2, "manifest", line_number)?;
                set_once(
                    &mut pack_id,
                    require_text(fields[1], "pack_id", "manifest")?,
                    "duplicate_pack_id",
                    "manifest",
                )?;
            }
            Some("redistribution") => {
                require_field_count(&fields, 3, "manifest", line_number)?;
                if fields[1] != "permitted" {
                    return Err(CompileError::new(
                        "redistribution_refused",
                        "manifest",
                        "redistribution policy must be explicitly permitted",
                    ));
                }
                set_once(
                    &mut redistribution_terms,
                    require_text(fields[2], "redistribution terms", "manifest")?,
                    "duplicate_redistribution",
                    "manifest",
                )?;
            }
            Some("citation") => {
                require_field_count(&fields, 2, "manifest", line_number)?;
                set_once(
                    &mut citation,
                    require_text(fields[1], "citation", "manifest")?,
                    "duplicate_citation",
                    "manifest",
                )?;
            }
            Some("license") => {
                require_field_count(&fields, 2, "manifest", line_number)?;
                set_once(
                    &mut license,
                    require_text(fields[1], "license", "manifest")?,
                    "duplicate_license",
                    "manifest",
                )?;
            }
            Some("source") => {
                require_field_count(&fields, 4, "manifest", line_number)?;
                if sources.len() == MAX_SOURCES {
                    return Err(CompileError::new(
                        "resource_limit",
                        "manifest:sources",
                        format!("at most {MAX_SOURCES} source files are admitted"),
                    ));
                }
                let id = require_identifier(fields[1], "source id", "manifest")?;
                let relative_path = safe_relative_path(fields[2])?;
                let profile = require_text(fields[3], "source profile", "manifest")?;
                sources.push(SourceSpec {
                    id,
                    relative_path,
                    profile,
                });
            }
            Some(other) => {
                return Err(CompileError::new(
                    "unknown_manifest_record",
                    format!("manifest:line:{line_number}"),
                    format!("unknown record type {other:?}"),
                ));
            }
            None => unreachable!("split always returns one field"),
        }
    }
    if sources.is_empty() {
        return Err(CompileError::new(
            "missing_sources",
            "manifest",
            "at least one source record is required",
        ));
    }
    sources.sort_by(|left, right| left.id.cmp(&right.id));
    let mut paths = BTreeSet::new();
    for pair in sources.windows(2) {
        if pair[0].id == pair[1].id {
            return Err(CompileError::new(
                "duplicate_source_id",
                format!("source:{}", pair[0].id),
                "source ids must be unique",
            ));
        }
    }
    for source in &sources {
        if !paths.insert(source.relative_path.clone()) {
            return Err(CompileError::new(
                "duplicate_source_path",
                format!("source:{}", source.id),
                "one physical source path may not be admitted under two ids",
            ));
        }
    }
    Ok(Manifest {
        pack_id: pack_id.ok_or_else(|| {
            CompileError::new("missing_pack_id", "manifest", "pack_id record is required")
        })?,
        redistribution_terms: redistribution_terms.ok_or_else(|| {
            CompileError::new(
                "missing_redistribution_terms",
                "manifest",
                "an explicit permitted redistribution record is required",
            )
        })?,
        citation: citation.ok_or_else(|| {
            CompileError::new(
                "missing_citation",
                "manifest",
                "citation record is required",
            )
        })?,
        license: license.ok_or_else(|| {
            CompileError::new("missing_license", "manifest", "license record is required")
        })?,
        sources,
    })
}

fn load_sources(
    manifest_path: &Path,
    manifest: &Manifest,
) -> Result<Vec<LoadedSource>, CompileError> {
    let parent = manifest_path.parent().unwrap_or_else(|| Path::new("."));
    let canonical_parent = parent.canonicalize().map_err(|error| {
        CompileError::new(
            "source_root_unavailable",
            "manifest",
            format!("cannot resolve manifest directory: {error}"),
        )
    })?;
    let mut loaded = Vec::with_capacity(manifest.sources.len());
    let mut total = 0usize;
    for spec in &manifest.sources {
        if spec.profile != MATERIAL_PROFILE {
            return Err(CompileError::new(
                "unsupported_source_profile",
                format!("source:{}", spec.id),
                format!(
                    "profile {:?} has no runtime-loadable v1 codec; supported profile is {MATERIAL_PROFILE:?}",
                    spec.profile
                ),
            ));
        }
        let path = resolve_safe_source(&canonical_parent, &spec.relative_path, &spec.id)?;
        let bytes = read_bounded_regular(&path, MAX_SOURCE_BYTES, &format!("source:{}", spec.id))?;
        let hash = hash_domain(SOURCE_FILE_DOMAIN, &bytes);
        total = total.checked_add(bytes.len()).ok_or_else(|| {
            CompileError::new(
                "resource_limit",
                "sources",
                "aggregate source byte count overflowed",
            )
            .with_input_hash(hash)
        })?;
        if total > MAX_TOTAL_SOURCE_BYTES {
            return Err(CompileError::new(
                "resource_limit",
                "sources",
                format!("aggregate source bytes exceed {MAX_TOTAL_SOURCE_BYTES}"),
            )
            .with_input_hash(hash));
        }
        let text = String::from_utf8(bytes.clone()).map_err(|_| {
            CompileError::new(
                "invalid_source_encoding",
                format!("source:{}", spec.id),
                "source must be valid UTF-8",
            )
            .with_input_hash(hash)
        })?;
        loaded.push(LoadedSource {
            spec: spec.clone(),
            bytes,
            text,
            hash,
        });
    }
    Ok(loaded)
}

fn resolve_safe_source(root: &Path, relative: &Path, id: &str) -> Result<PathBuf, CompileError> {
    let mut cursor = root.to_path_buf();
    for component in relative.components() {
        let Component::Normal(part) = component else {
            return Err(CompileError::new(
                "unsafe_source_path",
                format!("source:{id}"),
                "source path must contain only ordinary relative components",
            ));
        };
        cursor.push(part);
        let metadata = fs::symlink_metadata(&cursor).map_err(|error| {
            CompileError::new(
                "source_unavailable",
                format!("source:{id}"),
                format!("cannot inspect source path component: {error}"),
            )
        })?;
        if metadata.file_type().is_symlink() {
            return Err(CompileError::new(
                "unsafe_source_path",
                format!("source:{id}"),
                "symlinks are forbidden in source paths",
            ));
        }
    }
    let canonical = cursor.canonicalize().map_err(|error| {
        CompileError::new(
            "source_unavailable",
            format!("source:{id}"),
            format!("cannot resolve source file: {error}"),
        )
    })?;
    if !canonical.starts_with(root) {
        return Err(CompileError::new(
            "unsafe_source_path",
            format!("source:{id}"),
            "resolved source escapes the manifest directory",
        ));
    }
    Ok(canonical)
}

fn read_bounded_regular(path: &Path, limit: usize, subject: &str) -> Result<Vec<u8>, CompileError> {
    let mut file = File::open(path).map_err(|error| {
        CompileError::new(
            "input_unavailable",
            subject,
            format!("cannot open input: {error}"),
        )
    })?;
    let opened = file.metadata().map_err(|error| {
        CompileError::new(
            "input_unavailable",
            subject,
            format!("cannot inspect opened input: {error}"),
        )
    })?;
    let path_metadata = fs::symlink_metadata(path).map_err(|error| {
        CompileError::new(
            "input_unavailable",
            subject,
            format!("cannot inspect input: {error}"),
        )
    })?;
    if !opened.file_type().is_file()
        || !path_metadata.file_type().is_file()
        || path_metadata.file_type().is_symlink()
    {
        return Err(CompileError::new(
            "invalid_input_type",
            subject,
            "input must be a regular non-symlink file",
        ));
    }
    #[cfg(unix)]
    if opened.dev() != path_metadata.dev() || opened.ino() != path_metadata.ino() {
        return Err(CompileError::new(
            "input_changed_during_admission",
            subject,
            "opened descriptor no longer names the inspected source inode",
        ));
    }
    let declared = usize::try_from(opened.len()).unwrap_or(usize::MAX);
    if declared > limit {
        return Err(CompileError::new(
            "resource_limit",
            subject,
            format!("input exceeds {limit} bytes"),
        ));
    }
    let take_limit = u64::try_from(limit).unwrap_or(u64::MAX).saturating_add(1);
    let mut bytes = Vec::with_capacity(declared);
    file.take(take_limit)
        .read_to_end(&mut bytes)
        .map_err(|error| {
            CompileError::new(
                "input_unavailable",
                subject,
                format!("cannot read opened input: {error}"),
            )
        })?;
    if bytes.len() > limit {
        return Err(CompileError::new(
            "resource_limit",
            subject,
            format!("input grew beyond {limit} bytes while being read"),
        ));
    }
    Ok(bytes)
}

fn safe_relative_path(raw: &str) -> Result<PathBuf, CompileError> {
    if raw.is_empty() || raw.len() > MAX_TEXT_BYTES {
        return Err(CompileError::new(
            "unsafe_source_path",
            "manifest",
            "source path must be nonempty and bounded",
        ));
    }
    let path = PathBuf::from(raw);
    if path.is_absolute()
        || path
            .components()
            .any(|component| !matches!(component, Component::Normal(_)))
    {
        return Err(CompileError::new(
            "unsafe_source_path",
            "manifest",
            "source paths must be relative and contain no dot or parent components",
        ));
    }
    Ok(path)
}

fn set_once<T>(
    slot: &mut Option<T>,
    value: T,
    code: &'static str,
    subject: &str,
) -> Result<(), CompileError> {
    if slot.replace(value).is_some() {
        Err(CompileError::new(
            code,
            subject,
            "record may appear exactly once",
        ))
    } else {
        Ok(())
    }
}

fn require_field_count(
    fields: &[&str],
    expected: usize,
    source_id: &str,
    line: usize,
) -> Result<(), CompileError> {
    if fields.len() == expected {
        Ok(())
    } else {
        Err(CompileError::new(
            "invalid_record_shape",
            format!("source:{source_id}:line:{line}"),
            format!(
                "expected {expected} tab-separated fields, found {}",
                fields.len()
            ),
        ))
    }
}

fn require_text(raw: &str, field: &str, subject: &str) -> Result<String, CompileError> {
    if raw.trim().is_empty() {
        return Err(CompileError::new(
            "missing_required_field",
            subject,
            format!("{field} must be nonblank"),
        ));
    }
    if raw.len() > MAX_TEXT_BYTES {
        return Err(CompileError::new(
            "resource_limit",
            subject,
            format!("{field} exceeds {MAX_TEXT_BYTES} bytes"),
        ));
    }
    if raw.chars().any(char::is_control) {
        return Err(CompileError::new(
            "invalid_text_field",
            subject,
            format!("{field} contains a control character"),
        ));
    }
    Ok(raw.to_string())
}

fn require_identifier(raw: &str, field: &str, subject: &str) -> Result<String, CompileError> {
    if raw.is_empty() || raw.len() > MAX_IDENTIFIER_BYTES {
        return Err(CompileError::new(
            "invalid_identifier",
            subject,
            format!("{field} must contain 1..={MAX_IDENTIFIER_BYTES} bytes"),
        ));
    }
    if !raw
        .bytes()
        .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.'))
    {
        return Err(CompileError::new(
            "invalid_identifier",
            subject,
            format!("{field} may contain only ASCII letters, digits, dot, dash, and underscore"),
        ));
    }
    Ok(raw.to_string())
}

fn push_part(output: &mut Vec<u8>, part: &[u8]) {
    output.extend_from_slice(&u64::try_from(part.len()).unwrap_or(u64::MAX).to_le_bytes());
    output.extend_from_slice(part);
}

#[allow(clippy::too_many_lines)] // Exact record grammar and per-record refusals stay co-located.
fn parse_source(source: &LoadedSource, raw: &mut RawDatabase) -> Result<(), CompileError> {
    if source.text.as_bytes().contains(&0) {
        return Err(CompileError::new(
            "invalid_source_encoding",
            format!("source:{}", source.spec.id),
            "NUL bytes are forbidden",
        ));
    }
    let mut lines = source.text.lines();
    if lines.next() != Some(SOURCE_HEADER) {
        return Err(CompileError::new(
            "unsupported_source_schema",
            format!("source:{}", source.spec.id),
            format!("first line must be {SOURCE_HEADER:?}"),
        ));
    }
    for (offset, line) in lines.enumerate() {
        let line_number = offset + 2;
        let subject = format!("source:{}:line:{line_number}", source.spec.id);
        if line.len() > MAX_LINE_BYTES {
            return Err(CompileError::new(
                "resource_limit",
                subject,
                "line exceeds the public byte budget",
            ));
        }
        if line.is_empty() {
            return Err(CompileError::new(
                "invalid_source_record",
                subject,
                "blank records are forbidden",
            ));
        }
        raw.records = raw.records.checked_add(1).ok_or_else(|| {
            CompileError::new("resource_limit", "records", "record count overflowed")
        })?;
        if raw.records > MAX_RECORDS {
            return Err(CompileError::new(
                "resource_limit",
                "records",
                format!("at most {MAX_RECORDS} source records are admitted"),
            ));
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let record_hash = hash_domain(SOURCE_RECORD_DOMAIN, line.as_bytes());
        match fields.first().copied() {
            Some("observation") => {
                require_field_count(&fields, 5, &source.spec.id, line_number)?;
                let id = require_identifier(fields[1], "observation id", &subject)?;
                let observation = RawObservation {
                    id: id.clone(),
                    specimen: require_text(fields[2], "specimen", &subject)?,
                    method: require_text(fields[3], "method", &subject)?,
                    caveats: bounded_optional_text(fields[4], "caveats", &subject)?,
                    source_id: source.spec.id.clone(),
                    source_hash: source.hash,
                    record_hash,
                };
                if raw.observations.insert(id.clone(), observation).is_some() {
                    return Err(CompileError::new(
                        "duplicate_observation",
                        format!("observation:{id}"),
                        "observation ids are globally unique",
                    ));
                }
            }
            Some("scalar") => {
                require_field_count(&fields, 7, &source.spec.id, line_number)?;
                let id = require_identifier(fields[1], "claim id", &subject)?;
                let claim = RawClaim {
                    id: id.clone(),
                    observations: parse_identifier_list(
                        fields[2],
                        "observation references",
                        &subject,
                    )?,
                    property: require_identifier(fields[3], "property", &subject)?,
                    value: RawClaimValue::Scalar {
                        number: require_number_token(fields[4], "scalar value", &subject)?,
                        unit: require_unit(fields[5], "scalar unit", &subject)?,
                    },
                    interpolation: parse_interpolation(fields[6], &subject)?,
                    source_id: source.spec.id.clone(),
                    source_hash: source.hash,
                    record_hash,
                };
                insert_claim(raw, id, claim)?;
            }
            Some("curve") => {
                require_field_count(&fields, 9, &source.spec.id, line_number)?;
                let id = require_identifier(fields[1], "claim id", &subject)?;
                let points = parse_points(fields[7], &subject)?;
                raw.curve_knots = raw.curve_knots.checked_add(points.len()).ok_or_else(|| {
                    CompileError::new(
                        "resource_limit",
                        "curve_knots",
                        "aggregate knot count overflowed",
                    )
                })?;
                if raw.curve_knots > MAX_CURVE_KNOTS {
                    return Err(CompileError::new(
                        "resource_limit",
                        "curve_knots",
                        format!("at most {MAX_CURVE_KNOTS} curve knots are admitted"),
                    ));
                }
                let claim = RawClaim {
                    id: id.clone(),
                    observations: parse_identifier_list(
                        fields[2],
                        "observation references",
                        &subject,
                    )?,
                    property: require_identifier(fields[3], "property", &subject)?,
                    value: RawClaimValue::Curve {
                        abscissa: require_identifier(fields[4], "abscissa axis", &subject)?,
                        abscissa_unit: require_unit(fields[5], "abscissa unit", &subject)?,
                        ordinate_unit: require_unit(fields[6], "ordinate unit", &subject)?,
                        points,
                    },
                    interpolation: parse_interpolation(fields[8], &subject)?,
                    source_id: source.spec.id.clone(),
                    source_hash: source.hash,
                    record_hash,
                };
                insert_claim(raw, id, claim)?;
            }
            Some("uncertainty") => {
                require_field_count(&fields, 7, &source.spec.id, line_number)?;
                let claim = require_identifier(fields[1], "claim id", &subject)?;
                let kind = match fields[2] {
                    "unstated" => {
                        if fields[3..] != ["-", "-", "-", "-"] {
                            return Err(CompileError::new(
                                "invalid_uncertainty",
                                format!("claim:{claim}:uncertainty"),
                                "unstated uncertainty requires four '-' placeholders",
                            ));
                        }
                        RawUncertaintyKind::Unstated
                    }
                    "absolute" => RawUncertaintyKind::Absolute {
                        number: require_number_token(fields[3], "absolute uncertainty", &subject)?,
                        unit: require_unit(fields[4], "uncertainty unit", &subject)?,
                        confidence: require_number_token(fields[5], "confidence", &subject)?,
                        confidence_unit: require_unit(fields[6], "confidence unit", &subject)?,
                    },
                    "relative" => RawUncertaintyKind::Relative {
                        number: require_number_token(fields[3], "relative uncertainty", &subject)?,
                        unit: require_unit(fields[4], "uncertainty unit", &subject)?,
                        confidence: require_number_token(fields[5], "confidence", &subject)?,
                        confidence_unit: require_unit(fields[6], "confidence unit", &subject)?,
                    },
                    _ => {
                        return Err(CompileError::new(
                            "invalid_uncertainty",
                            format!("claim:{claim}:uncertainty"),
                            "kind must be unstated, absolute, or relative",
                        ));
                    }
                };
                if raw
                    .uncertainties
                    .insert(
                        claim.clone(),
                        RawUncertainty {
                            kind,
                            source_hash: record_hash,
                        },
                    )
                    .is_some()
                {
                    return Err(CompileError::new(
                        "duplicate_uncertainty",
                        format!("claim:{claim}:uncertainty"),
                        "each claim requires exactly one uncertainty record",
                    ));
                }
            }
            Some("validity") => {
                require_field_count(&fields, 6, &source.spec.id, line_number)?;
                let claim = require_identifier(fields[1], "claim id", &subject)?;
                let axis = require_identifier(fields[2], "validity axis", &subject)?;
                let validity = RawValidity {
                    axis: axis.clone(),
                    lower: require_number_token(fields[3], "validity lower bound", &subject)?,
                    upper: require_number_token(fields[4], "validity upper bound", &subject)?,
                    unit: require_unit(fields[5], "validity unit", &subject)?,
                    source_hash: record_hash,
                };
                if raw
                    .validities
                    .entry(claim.clone())
                    .or_default()
                    .insert(axis.clone(), validity)
                    .is_some()
                {
                    return Err(CompileError::new(
                        "duplicate_validity_axis",
                        format!("claim:{claim}:validity:{axis}"),
                        "one claim may name each validity axis once",
                    ));
                }
            }
            Some("frame") => {
                require_field_count(&fields, 4, &source.spec.id, line_number)?;
                let claim = require_identifier(fields[1], "claim id", &subject)?;
                let frame = RawFrame {
                    source: require_identifier(fields[2], "source frame", &subject)?,
                    target: require_identifier(fields[3], "target frame", &subject)?,
                    source_hash: record_hash,
                };
                if raw.frames.insert(claim.clone(), frame).is_some() {
                    return Err(CompileError::new(
                        "duplicate_frame",
                        format!("claim:{claim}:frame"),
                        "one claim may declare at most one paired frame mapping",
                    ));
                }
            }
            Some("joint") => {
                require_field_count(&fields, 7, &source.spec.id, line_number)?;
                let observation = require_identifier(fields[1], "observation id", &subject)?;
                let block_id = require_identifier(fields[2], "joint block id", &subject)?;
                let members = parse_member_list(fields[3], &subject)?;
                if members.len() > MAX_JOINT_MEMBERS {
                    return Err(CompileError::new(
                        "resource_limit",
                        format!("joint:{observation}:{block_id}"),
                        format!("at most {MAX_JOINT_MEMBERS} members are admitted per block"),
                    ));
                }
                let covariance = parse_numeric_list(fields[4], "covariance", &subject)?;
                let correlation = if fields[5] == "-" {
                    if fields[6] != "-" {
                        return Err(CompileError::new(
                            "invalid_correlation_unit",
                            format!("joint:{observation}:{block_id}"),
                            "absent correlation requires a '-' unit placeholder",
                        ));
                    }
                    None
                } else {
                    if fields[6] != "1" {
                        return Err(CompileError::new(
                            "invalid_correlation_unit",
                            format!("joint:{observation}:{block_id}"),
                            "correlation must explicitly declare the dimensionless unit '1'",
                        ));
                    }
                    Some(parse_numeric_list(fields[5], "correlation", &subject)?)
                };
                let joint = RawJoint {
                    observation: observation.clone(),
                    block_id: block_id.clone(),
                    members,
                    covariance,
                    correlation,
                    source_hash: record_hash,
                };
                if raw
                    .joints
                    .insert((observation.clone(), block_id.clone()), joint)
                    .is_some()
                {
                    return Err(CompileError::new(
                        "duplicate_joint_block",
                        format!("joint:{observation}:{block_id}"),
                        "joint block ids are unique within an observation",
                    ));
                }
            }
            Some("species" | "nasa9" | "kinetics") => {
                return Err(CompileError::new(
                    "unsupported_source_record",
                    subject,
                    "species, NASA-9, and kinetics records await a runtime-loadable codec",
                ));
            }
            Some(other) => {
                return Err(CompileError::new(
                    "unknown_source_record",
                    subject,
                    format!("unknown record type {other:?}"),
                ));
            }
            None => unreachable!("split always returns one field"),
        }
    }
    raw.decisions.push(Decision::admit(
        format!("source:{}", source.spec.id),
        "source_schema_admitted",
        "bounded material-tsv-v1 source parsed without inference",
        Some(source.hash),
    ));
    Ok(())
}

fn insert_claim(raw: &mut RawDatabase, id: String, claim: RawClaim) -> Result<(), CompileError> {
    if raw.claims.insert(id.clone(), claim).is_some() {
        Err(CompileError::new(
            "duplicate_claim",
            format!("claim:{id}"),
            "claim ids are globally unique",
        ))
    } else {
        Ok(())
    }
}

fn bounded_optional_text(raw: &str, field: &str, subject: &str) -> Result<String, CompileError> {
    if raw.len() > MAX_TEXT_BYTES {
        return Err(CompileError::new(
            "resource_limit",
            subject,
            format!("{field} exceeds {MAX_TEXT_BYTES} bytes"),
        ));
    }
    if raw.chars().any(char::is_control) {
        return Err(CompileError::new(
            "invalid_text_field",
            subject,
            format!("{field} contains a control character"),
        ));
    }
    Ok(raw.to_string())
}

fn require_number_token(raw: &str, field: &str, subject: &str) -> Result<String, CompileError> {
    if raw.is_empty() || raw.trim() != raw || raw.len() > 128 {
        return Err(CompileError::new(
            "invalid_number",
            subject,
            format!("{field} must be a bounded whitespace-free numeric token"),
        ));
    }
    parse_source_f64(raw, field, subject)?;
    Ok(raw.to_string())
}

fn require_unit(raw: &str, field: &str, subject: &str) -> Result<String, CompileError> {
    if raw.is_empty() {
        return Err(CompileError::new(
            "ambiguous_unit",
            subject,
            format!("{field} is required; dimensionless values must explicitly use '1'"),
        ));
    }
    if raw.trim() != raw || raw.len() > 256 || raw.chars().any(char::is_whitespace) {
        return Err(CompileError::new(
            "invalid_unit",
            subject,
            format!("{field} must be a bounded whitespace-free fs-qty expression"),
        ));
    }
    Ok(raw.to_string())
}

fn parse_identifier_list(
    raw: &str,
    field: &str,
    subject: &str,
) -> Result<Vec<String>, CompileError> {
    if raw.is_empty() {
        return Err(CompileError::new(
            "missing_observation_reference",
            subject,
            format!("{field} must name at least one observation"),
        ));
    }
    let mut values = Vec::new();
    for item in raw.split(',') {
        values.push(require_identifier(item, field, subject)?);
    }
    values.sort();
    if values.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(CompileError::new(
            "duplicate_reference",
            subject,
            format!("{field} contains a duplicate id"),
        ));
    }
    Ok(values)
}

fn parse_points(raw: &str, subject: &str) -> Result<Vec<(String, String)>, CompileError> {
    if raw.is_empty() {
        return Err(CompileError::new(
            "malformed_curve",
            subject,
            "curve point list is empty",
        ));
    }
    let mut points = Vec::new();
    for point in raw.split(',') {
        let (x, y) = point.split_once(':').ok_or_else(|| {
            CompileError::new(
                "malformed_curve",
                subject,
                "each curve point must be encoded as x:y",
            )
        })?;
        if y.contains(':') {
            return Err(CompileError::new(
                "malformed_curve",
                subject,
                "curve points contain too many ':' separators",
            ));
        }
        points.push((
            require_number_token(x, "curve abscissa", subject)?,
            require_number_token(y, "curve ordinate", subject)?,
        ));
    }
    if points.len() < 2 {
        return Err(CompileError::new(
            "malformed_curve",
            subject,
            "a curve requires at least two points",
        ));
    }
    Ok(points)
}

fn parse_member_list(raw: &str, subject: &str) -> Result<Vec<String>, CompileError> {
    if raw.is_empty() {
        return Err(CompileError::new(
            "invalid_joint_members",
            subject,
            "joint block must name at least one member",
        ));
    }
    let members: Vec<String> = raw.split(',').map(str::to_string).collect();
    if members.iter().any(String::is_empty) {
        return Err(CompileError::new(
            "invalid_joint_members",
            subject,
            "joint member references may not be blank",
        ));
    }
    Ok(members)
}

fn parse_numeric_list(raw: &str, field: &str, subject: &str) -> Result<Vec<String>, CompileError> {
    if raw.is_empty() {
        return Err(CompileError::new(
            "invalid_matrix",
            subject,
            format!("{field} matrix is empty"),
        ));
    }
    raw.split(',')
        .map(|item| require_number_token(item, field, subject))
        .collect()
}

fn parse_interpolation(raw: &str, subject: &str) -> Result<InterpolationPolicy, CompileError> {
    match raw {
        "linear" => Ok(InterpolationPolicy::LinearInside),
        "constant" => Ok(InterpolationPolicy::ConstantWithinValidity),
        "tabulated" => Ok(InterpolationPolicy::TabulatedOnly),
        _ => Err(CompileError::new(
            "invalid_interpolation",
            subject,
            "interpolation must be linear, constant, or tabulated",
        )),
    }
}

fn parse_source_f64(raw: &str, field: &str, subject: &str) -> Result<f64, CompileError> {
    let value = raw.parse::<f64>().map_err(|_| {
        CompileError::new(
            "invalid_number",
            subject,
            format!("{field} is not an exact finite f64 token"),
        )
    })?;
    if !value.is_finite() {
        return Err(CompileError::new(
            "non_finite_number",
            subject,
            format!("{field} must be finite"),
        ));
    }
    let exponent = raw.find(|character| matches!(character, 'e' | 'E'));
    let significand = exponent.map_or(raw, |index| &raw[..index]);
    if value == 0.0 && significand.bytes().any(|byte| matches!(byte, b'1'..=b'9')) {
        return Err(CompileError::new(
            "numeric_underflow",
            subject,
            format!("{field} has a nonzero significand that underflows f64"),
        ));
    }
    if value.to_bits() == (-0.0f64).to_bits() {
        return Err(CompileError::new(
            "negative_zero",
            subject,
            format!("{field} uses forbidden signed zero"),
        ));
    }
    Ok(value)
}

fn parse_confidence(raw: &str, unit: &str, subject: &str) -> Result<f64, CompileError> {
    if unit != "1" {
        return Err(CompileError::new(
            "invalid_confidence_unit",
            subject,
            "confidence must use the exact dimensionless unit '1' because schema v1 has no \
             separate confidence normalization receipt",
        ));
    }
    let quantity = parse_linear_quantity(raw, unit, subject)?;
    if quantity.dims != Dims::NONE {
        return Err(CompileError::new(
            "invalid_confidence_unit",
            subject,
            "confidence must explicitly use a dimensionless unit",
        ));
    }
    let value = quantity.value;
    if value > 0.0 && value < 1.0 {
        Ok(value)
    } else {
        Err(CompileError::new(
            "invalid_confidence",
            subject,
            "confidence must be strictly between zero and one",
        ))
    }
}

fn unit_transform(unit: &str, subject: &str) -> Result<(Dims, f64, f64), CompileError> {
    if unit == "1" {
        return Ok((Dims::NONE, 1.0, 0.0));
    }
    if !unit
        .chars()
        .next()
        .is_some_and(|first| first.is_alphabetic() || first == '%' || first == 'µ')
    {
        return Err(CompileError::new(
            "invalid_unit",
            subject,
            "unit must begin with a unit symbol; numeric continuations are forbidden",
        ));
    }
    let zero = parse_qty_with_budget(&format!("0 {unit}"), QTY_BUDGET).map_err(|_| {
        CompileError::new(
            "invalid_unit",
            subject,
            "fs-qty refused the explicit unit expression",
        )
    })?;
    let one = parse_qty_with_budget(&format!("1 {unit}"), QTY_BUDGET).map_err(|_| {
        CompileError::new(
            "invalid_unit",
            subject,
            "fs-qty refused the explicit unit expression",
        )
    })?;
    if zero.dims != one.dims {
        return Err(CompileError::new(
            "unstable_unit_dimensions",
            subject,
            "unit dimensions differ between zero and one probes",
        ));
    }
    let scale = one.value - zero.value;
    if !scale.is_finite() || scale <= 0.0 {
        return Err(CompileError::new(
            "invalid_unit_scale",
            subject,
            "unit transform must have a finite positive scale",
        ));
    }
    let offset = positive_zero(zero.value);
    Ok((zero.dims, scale, offset))
}

fn parse_quantity(number: &str, unit: &str, subject: &str) -> Result<ParsedQuantity, CompileError> {
    let source_value = parse_source_f64(number, "quantity", subject)?;
    let (dims, scale, offset) = unit_transform(unit, subject)?;
    let parsed_value = if unit == "1" {
        source_value
    } else {
        let parsed =
            parse_qty_with_budget(&format!("{number} {unit}"), QTY_BUDGET).map_err(|_| {
                CompileError::new(
                    "invalid_quantity",
                    subject,
                    "fs-qty refused the explicit number and unit",
                )
            })?;
        if parsed.dims != dims {
            return Err(CompileError::new(
                "unstable_unit_dimensions",
                subject,
                "literal dimensions differ from the unit probes",
            ));
        }
        parsed.value
    };
    if !parsed_value.is_finite() {
        return Err(CompileError::new(
            "non_finite_normalization",
            subject,
            "normalized value is not finite",
        ));
    }
    let receipt_value = source_value * scale + offset;
    let tolerance = f64::EPSILON * 8.0 * parsed_value.abs().max(receipt_value.abs()).max(1.0);
    if !receipt_value.is_finite() || (receipt_value - parsed_value).abs() > tolerance {
        return Err(CompileError::new(
            "non_reproducible_unit_transform",
            subject,
            "derived affine receipt disagrees materially with the independently parsed quantity",
        ));
    }
    if offset == 0.0 && source_value != 0.0 && receipt_value == 0.0 {
        return Err(CompileError::new(
            "normalization_underflow",
            subject,
            "nonzero source quantity underflows the normalized representation",
        ));
    }
    Ok(ParsedQuantity {
        value: positive_zero(receipt_value),
        dims,
        scale,
        offset,
        unit: unit.to_string(),
        literal_hash: source_literal_hash(number, unit),
    })
}

fn parse_linear_quantity(
    number: &str,
    unit: &str,
    subject: &str,
) -> Result<ParsedQuantity, CompileError> {
    let source_value = parse_source_f64(number, "quantity", subject)?;
    let (dims, scale, _) = unit_transform(unit, subject)?;
    let value = source_value * scale;
    if !value.is_finite() {
        return Err(CompileError::new(
            "non_finite_normalization",
            subject,
            "normalized linear magnitude is not finite",
        ));
    }
    if source_value != 0.0 && value == 0.0 {
        return Err(CompileError::new(
            "normalization_underflow",
            subject,
            "nonzero source magnitude underflows the normalized representation",
        ));
    }
    Ok(ParsedQuantity {
        value: positive_zero(value),
        dims,
        scale,
        offset: 0.0,
        unit: unit.to_string(),
        literal_hash: source_literal_hash(number, unit),
    })
}

fn source_literal_hash(number: &str, unit: &str) -> ContentHash {
    let mut bytes = Vec::new();
    push_part(&mut bytes, number.as_bytes());
    push_part(&mut bytes, unit.as_bytes());
    hash_domain(SOURCE_LITERAL_DOMAIN, &bytes)
}

fn positive_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

#[allow(clippy::too_many_lines)] // The ordered admission pipeline is easier to audit in one pass.
fn compile_manifest(manifest_path: &Path) -> Result<CompileOutput, CompileError> {
    let manifest_bytes = read_bounded_regular(manifest_path, MAX_MANIFEST_BYTES, "manifest")?;
    let manifest_snapshot = hash_domain(SOURCE_ENVELOPE_DOMAIN, &manifest_bytes);
    let manifest_text = std::str::from_utf8(&manifest_bytes).map_err(|_| {
        CompileError::new(
            "invalid_manifest_encoding",
            "manifest",
            "manifest must be valid UTF-8",
        )
        .with_input_hash(manifest_snapshot)
    })?;
    let manifest =
        parse_manifest(manifest_text).map_err(|error| error.with_input_hash(manifest_snapshot))?;
    let sources = load_sources(manifest_path, &manifest)
        .map_err(|error| error.with_input_hash(manifest_snapshot))?;
    let source_artifact = source_envelope_hash(&manifest, &sources);
    let result = (|| -> Result<CompileOutput, CompileError> {
        let mut raw = RawDatabase::default();
        for source in &sources {
            parse_source(source, &mut raw)?;
        }
        validate_raw_references(&raw)?;
        preflight_provenance_budget(&manifest, &raw)?;
        preflight_normalization_budget(&raw)?;

        let mut claims = ClaimSet::new();
        let mut observation_ids = BTreeMap::new();
        let mut components = BTreeMap::new();
        let mut normalizations = Vec::new();
        let mut decisions = raw.decisions;
        decisions.push(Decision::admit(
            "manifest:redistribution",
            "redistribution_permitted",
            "explicit permitted policy and nonblank retained terms admitted before compilation",
            Some(source_artifact),
        ));
        decisions.push(Decision::admit(
            "manifest:provenance",
            "licensed_citation_admitted",
            "nonblank source citation and license attached to every emitted datum",
            Some(source_artifact),
        ));
        let mut validity_dims: BTreeMap<String, Dims> = BTreeMap::new();

        for (local_id, observation) in &raw.observations {
            debug_assert_eq!(local_id, &observation.id);
            let id = claims
                .register_observation(ObservationDataset {
                    specimen: observation.specimen.clone(),
                    method: observation.method.clone(),
                    artifact: observation.record_hash,
                    caveats: observation.caveats.clone(),
                    provenance: provenance_for(
                        &manifest,
                        &observation.source_id,
                        observation.source_hash,
                    ),
                })
                .map_err(|error| {
                    CompileError::new(
                        "observation_admission_refused",
                        format!("observation:{local_id}"),
                        error.to_string(),
                    )
                })?;
            observation_ids.insert(local_id.clone(), id);
            decisions.push(Decision::admit(
                format!("observation:{local_id}"),
                "observation_normalized",
                "registered provenance-complete observation dataset",
                Some(observation.record_hash),
            ));
        }

        for (local_id, raw_claim) in &raw.claims {
            debug_assert_eq!(local_id, &raw_claim.id);
            let frame = raw.frames.get(local_id);
            let (value, parsed_components) = compile_value(raw_claim)?;
            let value_dims = value.dims();
            let (validity, validity_receipts) =
                compile_validity(local_id, raw.validities.get(local_id), &mut validity_dims)?;
            if let PropertyValue::Curve {
                abscissa,
                abscissa_dims,
                ..
            } = &value
            {
                let declared = raw
                    .validities
                    .get(local_id)
                    .and_then(|axes| axes.get(abscissa))
                    .ok_or_else(|| {
                        CompileError::new(
                            "missing_curve_validity",
                            format!("claim:{local_id}"),
                            "curve abscissa must have an explicit validity interval",
                        )
                    })?;
                let parsed = parse_quantity(
                    &declared.lower,
                    &declared.unit,
                    &format!("claim:{local_id}:validity:{abscissa}:lower"),
                )?;
                if parsed.dims != *abscissa_dims {
                    return Err(CompileError::new(
                        "curve_validity_dims_mismatch",
                        format!("claim:{local_id}"),
                        "curve abscissa dimensions disagree with its validity axis",
                    ));
                }
            }
            let raw_uncertainty = raw.uncertainties.get(local_id).ok_or_else(|| {
                CompileError::new(
                    "missing_uncertainty",
                    format!("claim:{local_id}:uncertainty"),
                    "each source claim must explicitly declare its uncertainty, including unstated",
                )
            })?;
            let (uncertainty, uncertainty_quantity) =
                compile_uncertainty(local_id, raw_uncertainty, value_dims)?;
            let mut observations = raw_claim
                .observations
                .iter()
                .map(|reference| {
                    observation_ids.get(reference).copied().ok_or_else(|| {
                        CompileError::new(
                            "unknown_observation",
                            format!("claim:{local_id}"),
                            format!("observation reference {reference:?} is not registered"),
                        )
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            observations.sort_unstable();
            let component_observations: BTreeSet<ObservationId> =
                observations.iter().copied().collect();
            let id = claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(raw_claim.property.clone(), value_dims),
                    value,
                    validity,
                    uncertainty,
                    interpolation: raw_claim.interpolation,
                    observations,
                    provenance: provenance_for(
                        &manifest,
                        &raw_claim.source_id,
                        raw_claim.source_hash,
                    ),
                })
                .map_err(|error| {
                    CompileError::new(
                        "claim_admission_refused",
                        format!("claim:{local_id}"),
                        error.to_string(),
                    )
                })?;
            for (component, quantity) in parsed_components {
                let member = statistic_member(id, component);
                normalizations.push(quantity_receipt(
                    NormalizationTarget::ClaimValue(member),
                    &quantity,
                    frame,
                ));
                components.insert(
                    (local_id.clone(), component),
                    ComponentMeta {
                        member,
                        dims: quantity.dims,
                        scale: quantity.scale,
                        unit: quantity.unit.clone(),
                        observations: component_observations.clone(),
                    },
                );
                decisions.push(Decision::admit(
                    normalization_subject(local_id, component),
                    "quantity_normalized",
                    "explicit source quantity converted to the six-base SI target",
                    Some(quantity.literal_hash),
                ));
            }
            if let Some(quantity) = uncertainty_quantity {
                normalizations.push(quantity_receipt(
                    NormalizationTarget::ClaimUncertainty { claim: id },
                    &quantity,
                    frame,
                ));
                decisions.push(Decision::admit(
                    format!("normalization:{local_id}:uncertainty"),
                    "uncertainty_normalized",
                    "uncertainty magnitude converted without an affine translation",
                    Some(quantity.literal_hash),
                ));
            }
            decisions.push(Decision::admit(
                format!("claim:{local_id}:uncertainty"),
                "uncertainty_policy_admitted",
                "source uncertainty kind and confidence were admitted without inference",
                Some(raw_uncertainty.source_hash),
            ));
            for (axis, lower, upper) in validity_receipts {
                let validity_source_hash = raw
                    .validities
                    .get(local_id)
                    .and_then(|values| values.get(&axis))
                    .map(|value| value.source_hash)
                    .ok_or_else(|| {
                        CompileError::new(
                            "internal_linkage_error",
                            format!("claim:{local_id}:validity:{axis}"),
                            "compiled validity receipt lost its admitted source record",
                        )
                    })?;
                normalizations.push(quantity_receipt(
                    NormalizationTarget::ValidityBound {
                        claim: id,
                        axis: axis.clone(),
                        side: ValidityBoundSide::Lower,
                    },
                    &lower,
                    None,
                ));
                normalizations.push(quantity_receipt(
                    NormalizationTarget::ValidityBound {
                        claim: id,
                        axis: axis.clone(),
                        side: ValidityBoundSide::Upper,
                    },
                    &upper,
                    None,
                ));
                decisions.push(Decision::admit(
                    format!("normalization:{local_id}:validity:{axis}"),
                    "validity_normalized",
                    "both validity endpoints converted under one explicit unit transform",
                    Some(validity_source_hash),
                ));
            }
            decisions.push(Decision::admit(
                format!("claim:{local_id}"),
                "claim_normalized",
                "inserted one immutable dims-checked scalar/curve claim",
                Some(raw_claim.record_hash),
            ));
            if let Some(frame) = frame {
                decisions.push(Decision::admit(
                    format!("claim:{local_id}:frame"),
                    "frame_pair_retained",
                    "paired source and target frames retained as normalization metadata",
                    Some(frame.source_hash),
                ));
            }
        }

        let mut joint_statistics = Vec::new();
        for ((observation, block_id), raw_joint) in &raw.joints {
            let observation_id = observation_ids.get(observation).copied().ok_or_else(|| {
                CompileError::new(
                    "unknown_observation",
                    format!("joint:{observation}:{block_id}"),
                    "joint block references an unknown observation",
                )
            })?;
            let (joint, receipts, joint_decisions) =
                compile_joint(raw_joint, observation_id, &components)?;
            joint_statistics.push(joint);
            normalizations.extend(receipts);
            decisions.extend(joint_decisions);
            decisions.push(Decision::admit(
                format!("joint:{observation}:{block_id}"),
                "joint_statistics_normalized",
                "member order and packed matrices were permuted together before SI scaling",
                Some(raw_joint.source_hash),
            ));
        }

        let pack = NormalizedPack::new(
            manifest.pack_id,
            COMPILER_ID,
            source_artifact,
            manifest.redistribution_terms,
            claims,
            joint_statistics,
            normalizations,
        )
        .map_err(|error| {
            CompileError::new(
                "normalized_pack_refused",
                "pack",
                format!("runtime pack admission refused the compiled artifact: {error}"),
            )
        })?;
        let bytes = pack.to_bytes();
        let pack_hash = pack.content_hash();
        let decoded = NormalizedPack::from_bytes_verified(pack_hash, &bytes).map_err(|error| {
            CompileError::new(
                "self_verification_failed",
                "pack",
                format!("fresh artifact failed verified runtime decode: {error}"),
            )
        })?;
        if decoded.to_bytes() != bytes {
            return Err(CompileError::new(
                "self_verification_failed",
                "pack",
                "verified runtime decode did not reproduce canonical bytes",
            ));
        }
        decisions.push(Decision::admit(
            "pack",
            "runtime_pack_self_verified",
            "canonical bytes decoded under their externally pinned content hash",
            Some(source_artifact),
        ));
        for decision in &mut decisions {
            decision.pack_hash = Some(pack_hash);
        }
        sort_decisions(&mut decisions);
        Ok(CompileOutput {
            pack,
            bytes,
            decisions,
        })
    })();
    result.map_err(|error| error.with_input_hash(source_artifact))
}

fn source_envelope_hash(manifest: &Manifest, sources: &[LoadedSource]) -> ContentHash {
    let mut payload = Vec::new();
    for part in [
        manifest.pack_id.as_bytes(),
        manifest.redistribution_terms.as_bytes(),
        manifest.citation.as_bytes(),
        manifest.license.as_bytes(),
    ] {
        push_part(&mut payload, part);
    }
    for source in sources {
        push_part(&mut payload, source.spec.id.as_bytes());
        push_part(&mut payload, source.spec.profile.as_bytes());
        push_part(&mut payload, &source.bytes);
    }
    hash_domain(SOURCE_ENVELOPE_DOMAIN, &payload)
}

fn provenance_for(manifest: &Manifest, source_id: &str, source_hash: ContentHash) -> Provenance {
    Provenance {
        source: format!("{} [source:{source_id}]", manifest.citation),
        license: manifest.license.clone(),
        artifact: Some(source_hash),
    }
}

fn validate_raw_references(raw: &RawDatabase) -> Result<(), CompileError> {
    if raw.observations.is_empty() {
        return Err(CompileError::new(
            "missing_observations",
            "sources",
            "material pack must contain at least one observation",
        ));
    }
    if raw.claims.is_empty() {
        return Err(CompileError::new(
            "missing_claims",
            "sources",
            "material pack must contain at least one scalar or curve claim",
        ));
    }
    for claim in raw.uncertainties.keys() {
        if !raw.claims.contains_key(claim) {
            return Err(CompileError::new(
                "unknown_claim",
                format!("claim:{claim}:uncertainty"),
                "uncertainty record references an unknown claim",
            ));
        }
    }
    for claim in raw.validities.keys() {
        if !raw.claims.contains_key(claim) {
            return Err(CompileError::new(
                "unknown_claim",
                format!("claim:{claim}:validity"),
                "validity record references an unknown claim",
            ));
        }
    }
    for claim in raw.frames.keys() {
        if !raw.claims.contains_key(claim) {
            return Err(CompileError::new(
                "unknown_claim",
                format!("claim:{claim}:frame"),
                "frame record references an unknown claim",
            ));
        }
    }
    for raw_claim in raw.claims.values() {
        for observation in &raw_claim.observations {
            if !raw.observations.contains_key(observation) {
                return Err(CompileError::new(
                    "unknown_observation",
                    format!("claim:{}", raw_claim.id),
                    format!("observation reference {observation:?} is not declared"),
                ));
            }
        }
    }
    for raw_joint in raw.joints.values() {
        if !raw.observations.contains_key(&raw_joint.observation) {
            return Err(CompileError::new(
                "unknown_observation",
                format!("joint:{}:{}", raw_joint.observation, raw_joint.block_id),
                "joint block references an unknown observation",
            ));
        }
    }
    Ok(())
}

fn preflight_provenance_budget(manifest: &Manifest, raw: &RawDatabase) -> Result<(), CompileError> {
    let source_ids = raw
        .observations
        .values()
        .map(|value| value.source_id.as_str())
        .chain(raw.claims.values().map(|value| value.source_id.as_str()));
    let mut repeated = 0usize;
    for source_id in source_ids {
        // Two wire string lengths, the artifact tag/hash, and the literal
        // " [source:...]" expansion are included. This is a conservative
        // lower-level budget solely for repeated global provenance; ordinary
        // source payloads remain bounded independently.
        let record = manifest
            .citation
            .len()
            .checked_add(manifest.license.len())
            .and_then(|value| value.checked_add(source_id.len()))
            .and_then(|value| value.checked_add(64))
            .ok_or_else(|| {
                CompileError::new(
                    "resource_limit",
                    "provenance",
                    "repeated provenance byte estimate overflowed",
                )
            })?;
        repeated = repeated.checked_add(record).ok_or_else(|| {
            CompileError::new(
                "resource_limit",
                "provenance",
                "aggregate repeated provenance bytes overflowed",
            )
        })?;
        if repeated > MAX_REPEATED_PROVENANCE_BYTES {
            return Err(CompileError::new(
                "resource_limit",
                "provenance",
                format!(
                    "repeated citation/license expansion exceeds \
                     {MAX_REPEATED_PROVENANCE_BYTES} bytes"
                ),
            ));
        }
    }
    Ok(())
}

fn preflight_normalization_budget(raw: &RawDatabase) -> Result<(), CompileError> {
    let mut receipts = 0usize;
    let mut add = |count: usize, subject: &str| -> Result<(), CompileError> {
        receipts = receipts.checked_add(count).ok_or_else(|| {
            CompileError::new(
                "resource_limit",
                subject,
                "normalization receipt count overflowed",
            )
        })?;
        if receipts > MAX_NORMALIZATION_RECEIPTS {
            return Err(CompileError::new(
                "resource_limit",
                "normalizations",
                format!(
                    "compiled source requires {receipts} receipts; schema v1 admits at most \
                     {MAX_NORMALIZATION_RECEIPTS}"
                ),
            ));
        }
        Ok(())
    };
    for (claim_id, claim) in &raw.claims {
        let value_receipts = match &claim.value {
            RawClaimValue::Scalar { .. } => 1,
            RawClaimValue::Curve { points, .. } => {
                points.len().checked_mul(2).ok_or_else(|| {
                    CompileError::new(
                        "resource_limit",
                        format!("claim:{claim_id}"),
                        "curve normalization count overflowed",
                    )
                })?
            }
        };
        add(value_receipts, &format!("claim:{claim_id}"))?;
        if raw
            .uncertainties
            .get(claim_id)
            .is_some_and(|value| !matches!(&value.kind, RawUncertaintyKind::Unstated))
        {
            add(1, &format!("claim:{claim_id}:uncertainty"))?;
        }
        let validity_receipts = raw
            .validities
            .get(claim_id)
            .map_or(0, |axes| axes.len().saturating_mul(2));
        add(validity_receipts, &format!("claim:{claim_id}:validity"))?;
    }
    for ((observation, block_id), joint) in &raw.joints {
        let packed = joint
            .members
            .len()
            .checked_mul(joint.members.len() + 1)
            .and_then(|value| value.checked_div(2))
            .ok_or_else(|| {
                CompileError::new(
                    "resource_limit",
                    format!("joint:{observation}:{block_id}"),
                    "joint normalization count overflowed",
                )
            })?;
        add(packed, &format!("joint:{observation}:{block_id}"))?;
    }
    Ok(())
}

fn compile_value(
    claim: &RawClaim,
) -> Result<(PropertyValue, Vec<(LocalComponent, ParsedQuantity)>), CompileError> {
    let subject = format!("claim:{}", claim.id);
    match &claim.value {
        RawClaimValue::Scalar { number, unit } => {
            if claim.interpolation == InterpolationPolicy::LinearInside {
                return Err(CompileError::new(
                    "incompatible_interpolation",
                    &subject,
                    "scalar claims admit only constant or tabulated interpolation",
                ));
            }
            let quantity = parse_quantity(number, unit, &format!("{subject}:value"))?;
            Ok((
                PropertyValue::Scalar {
                    value: quantity.value,
                    dims: quantity.dims,
                },
                vec![(LocalComponent::Scalar, quantity)],
            ))
        }
        RawClaimValue::Curve {
            abscissa,
            abscissa_unit,
            ordinate_unit,
            points,
        } => {
            if claim.interpolation == InterpolationPolicy::ConstantWithinValidity {
                return Err(CompileError::new(
                    "incompatible_interpolation",
                    &subject,
                    "curve claims admit only linear or tabulated interpolation",
                ));
            }
            let mut knots = Vec::with_capacity(points.len());
            let mut components = Vec::with_capacity(points.len() * 2);
            let mut abscissa_dims = None;
            let mut ordinate_dims = None;
            for (index, (x, y)) in points.iter().enumerate() {
                let knot = u32::try_from(index).map_err(|_| {
                    CompileError::new(
                        "resource_limit",
                        &subject,
                        "curve knot index exceeds the runtime wire format",
                    )
                })?;
                let parsed_x =
                    parse_quantity(x, abscissa_unit, &format!("{subject}:curve:x:{index}"))?;
                let parsed_y =
                    parse_quantity(y, ordinate_unit, &format!("{subject}:curve:y:{index}"))?;
                if let Some(expected) = abscissa_dims {
                    if expected != parsed_x.dims {
                        return Err(CompileError::new(
                            "curve_dims_mismatch",
                            &subject,
                            "curve abscissa dimensions changed between knots",
                        ));
                    }
                } else {
                    abscissa_dims = Some(parsed_x.dims);
                }
                if let Some(expected) = ordinate_dims {
                    if expected != parsed_y.dims {
                        return Err(CompileError::new(
                            "curve_dims_mismatch",
                            &subject,
                            "curve ordinate dimensions changed between knots",
                        ));
                    }
                } else {
                    ordinate_dims = Some(parsed_y.dims);
                }
                knots.push((parsed_x.value, parsed_y.value));
                components.push((LocalComponent::CurveAbscissa(knot), parsed_x));
                components.push((LocalComponent::CurveOrdinate(knot), parsed_y));
            }
            if !knots.windows(2).all(|pair| pair[0].0 < pair[1].0) {
                return Err(CompileError::new(
                    "malformed_curve",
                    &subject,
                    "normalized curve abscissae must be strictly increasing",
                ));
            }
            let abscissa_dims = abscissa_dims.ok_or_else(|| {
                CompileError::new(
                    "malformed_curve",
                    &subject,
                    "curve has no parsed abscissa dimensions",
                )
            })?;
            let ordinate_dims = ordinate_dims.ok_or_else(|| {
                CompileError::new(
                    "malformed_curve",
                    &subject,
                    "curve has no parsed ordinate dimensions",
                )
            })?;
            Ok((
                PropertyValue::Curve {
                    abscissa: abscissa.clone(),
                    abscissa_dims,
                    knots,
                    dims: ordinate_dims,
                },
                components,
            ))
        }
    }
}

fn compile_validity(
    claim_id: &str,
    axes: Option<&BTreeMap<String, RawValidity>>,
    global_dims: &mut BTreeMap<String, Dims>,
) -> Result<
    (
        ValidityDomain,
        Vec<(String, ParsedQuantity, ParsedQuantity)>,
    ),
    CompileError,
> {
    let mut domain = ValidityDomain::unconstrained();
    let mut receipts = Vec::new();
    for (axis, raw) in axes.into_iter().flatten() {
        debug_assert_eq!(axis, &raw.axis);
        let lower = parse_quantity(
            &raw.lower,
            &raw.unit,
            &format!("claim:{claim_id}:validity:{axis}:lower"),
        )?;
        let upper = parse_quantity(
            &raw.upper,
            &raw.unit,
            &format!("claim:{claim_id}:validity:{axis}:upper"),
        )?;
        if lower.dims != upper.dims {
            return Err(CompileError::new(
                "validity_dims_mismatch",
                format!("claim:{claim_id}:validity:{axis}"),
                "validity endpoints have different dimensions",
            ));
        }
        if lower.value > upper.value {
            return Err(CompileError::new(
                "reversed_validity",
                format!("claim:{claim_id}:validity:{axis}"),
                "normalized lower endpoint exceeds upper endpoint",
            ));
        }
        if let Some(previous) = global_dims.insert(axis.clone(), lower.dims)
            && previous != lower.dims
        {
            return Err(CompileError::new(
                "validity_axis_dims_mismatch",
                format!("validity:{axis}"),
                "one named validity axis has contradictory dimensions across claims",
            ));
        }
        domain = domain.with(axis.clone(), lower.value, upper.value);
        receipts.push((axis.clone(), lower, upper));
    }
    Ok((domain, receipts))
}

fn compile_uncertainty(
    claim_id: &str,
    raw: &RawUncertainty,
    value_dims: Dims,
) -> Result<(UncertaintyModel, Option<ParsedQuantity>), CompileError> {
    let subject = format!("claim:{claim_id}:uncertainty");
    match &raw.kind {
        RawUncertaintyKind::Unstated => Ok((UncertaintyModel::Unstated, None)),
        RawUncertaintyKind::Absolute {
            number,
            unit,
            confidence,
            confidence_unit,
        } => {
            let quantity = parse_linear_quantity(number, unit, &subject)?;
            if quantity.dims != value_dims {
                return Err(CompileError::new(
                    "uncertainty_dims_mismatch",
                    &subject,
                    "absolute uncertainty dimensions disagree with the claim value",
                ));
            }
            Ok((
                UncertaintyModel::HalfWidth {
                    half_width: quantity.value,
                    confidence: parse_confidence(confidence, confidence_unit, &subject)?,
                },
                Some(quantity),
            ))
        }
        RawUncertaintyKind::Relative {
            number,
            unit,
            confidence,
            confidence_unit,
        } => {
            let quantity = parse_linear_quantity(number, unit, &subject)?;
            if quantity.dims != Dims::NONE {
                return Err(CompileError::new(
                    "uncertainty_dims_mismatch",
                    &subject,
                    "relative uncertainty must use an explicit dimensionless unit",
                ));
            }
            Ok((
                UncertaintyModel::RelativeHalfWidth {
                    fraction: quantity.value,
                    confidence: parse_confidence(confidence, confidence_unit, &subject)?,
                },
                Some(quantity),
            ))
        }
    }
}

fn statistic_member(claim: ClaimId, component: LocalComponent) -> StatisticMember {
    match component {
        LocalComponent::Scalar => StatisticMember::scalar(claim),
        LocalComponent::CurveAbscissa(knot) => StatisticMember::curve_abscissa(claim, knot),
        LocalComponent::CurveOrdinate(knot) => StatisticMember::curve_ordinate(claim, knot),
    }
}

fn quantity_receipt(
    target: NormalizationTarget,
    quantity: &ParsedQuantity,
    frame: Option<&RawFrame>,
) -> NormalizationReceipt {
    NormalizationReceipt::new(
        target,
        quantity.literal_hash,
        quantity.dims,
        quantity.scale,
        quantity.offset,
        quantity.unit.clone(),
        MATDB_PACK_TARGET_BASIS,
        frame.map(|value| value.source.clone()),
        frame.map(|value| value.target.clone()),
    )
}

fn normalization_subject(claim: &str, component: LocalComponent) -> String {
    match component {
        LocalComponent::Scalar => format!("normalization:{claim}:scalar"),
        LocalComponent::CurveAbscissa(knot) => {
            format!("normalization:{claim}:curve:x:{knot}")
        }
        LocalComponent::CurveOrdinate(knot) => {
            format!("normalization:{claim}:curve:y:{knot}")
        }
    }
}

fn compile_joint(
    raw: &RawJoint,
    observation: ObservationId,
    components: &BTreeMap<(String, LocalComponent), ComponentMeta>,
) -> Result<(JointStatistics, Vec<NormalizationReceipt>, Vec<Decision>), CompileError> {
    let subject = format!("joint:{}:{}", raw.observation, raw.block_id);
    let size = raw.members.len();
    let packed_len = size
        .checked_mul(size + 1)
        .and_then(|value| value.checked_div(2))
        .ok_or_else(|| {
            CompileError::new(
                "resource_limit",
                &subject,
                "packed matrix length overflowed",
            )
        })?;
    if raw.covariance.len() != packed_len {
        return Err(CompileError::new(
            "invalid_covariance_shape",
            &subject,
            format!(
                "{} members require {packed_len} packed covariance entries, found {}",
                size,
                raw.covariance.len()
            ),
        ));
    }
    if let Some(correlation) = &raw.correlation
        && correlation.len() != packed_len
    {
        return Err(CompileError::new(
            "invalid_correlation_shape",
            &subject,
            format!(
                "{} members require {packed_len} packed correlation entries, found {}",
                size,
                correlation.len()
            ),
        ));
    }

    let mut resolved = Vec::with_capacity(size);
    for (old_index, reference) in raw.members.iter().enumerate() {
        let (claim, component) = parse_member_reference(reference, &subject)?;
        let meta = components
            .get(&(claim.clone(), component))
            .cloned()
            .ok_or_else(|| {
                CompileError::new(
                    "unknown_joint_member",
                    &subject,
                    format!("member {reference:?} does not name a compiled claim component"),
                )
            })?;
        if !meta.observations.contains(&observation) {
            return Err(CompileError::new(
                "joint_observation_mismatch",
                &subject,
                format!("member {reference:?} is not backed by the joint block observation"),
            ));
        }
        resolved.push((old_index, meta));
    }
    resolved.sort_by_key(|(_, meta)| meta.member);
    if resolved
        .windows(2)
        .any(|pair| pair[0].1.member == pair[1].1.member)
    {
        return Err(CompileError::new(
            "duplicate_joint_member",
            &subject,
            "two source member references resolve to the same claim component",
        ));
    }

    let raw_covariance = raw
        .covariance
        .iter()
        .map(|token| parse_source_f64(token, "covariance entry", &subject))
        .collect::<Result<Vec<_>, _>>()?;
    let raw_correlation = raw
        .correlation
        .as_ref()
        .map(|values| {
            values
                .iter()
                .map(|token| parse_source_f64(token, "correlation entry", &subject))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()?;
    let mut members = Vec::with_capacity(size);
    let mut covariance = Vec::with_capacity(packed_len);
    let mut correlation = raw_correlation
        .as_ref()
        .map(|_| Vec::with_capacity(packed_len));
    let mut receipts = Vec::with_capacity(packed_len);
    let mut decisions = Vec::with_capacity(packed_len);
    for (_, meta) in &resolved {
        members.push(meta.member);
    }
    for row in 0..size {
        for column in 0..=row {
            let (old_row, row_meta) = &resolved[row];
            let (old_column, column_meta) = &resolved[column];
            let old_index = packed_index(*old_row, *old_column);
            let raw_value = raw_covariance[old_index];
            let scale = row_meta.scale * column_meta.scale;
            let value = raw_value * scale;
            if !scale.is_finite() || scale <= 0.0 || !value.is_finite() {
                return Err(CompileError::new(
                    "non_finite_covariance_normalization",
                    &subject,
                    "covariance scaling produced a non-finite or non-positive transform",
                ));
            }
            if raw_value != 0.0 && value == 0.0 {
                return Err(CompileError::new(
                    "normalization_underflow",
                    &subject,
                    "nonzero covariance entry underflows the normalized representation",
                ));
            }
            let dims = row_meta
                .dims
                .checked_plus(column_meta.dims)
                .ok_or_else(|| {
                    CompileError::new(
                        "covariance_dimension_overflow",
                        &subject,
                        "member dimension product exceeds the six-base exponent range",
                    )
                })?;
            let basis = format!("({})*({})", row_meta.unit, column_meta.unit);
            let literal_hash = source_literal_hash(&raw.covariance[old_index], &basis);
            let row = u32::try_from(row).map_err(|_| {
                CompileError::new(
                    "resource_limit",
                    &subject,
                    "joint covariance row exceeds the runtime wire index",
                )
            })?;
            let column = u32::try_from(column).map_err(|_| {
                CompileError::new(
                    "resource_limit",
                    &subject,
                    "joint covariance column exceeds the runtime wire index",
                )
            })?;
            covariance.push(positive_zero(value));
            receipts.push(NormalizationReceipt::new(
                NormalizationTarget::JointCovariance {
                    observation,
                    block_id: raw.block_id.clone(),
                    row,
                    column,
                },
                literal_hash,
                dims,
                scale,
                0.0,
                basis,
                MATDB_PACK_TARGET_BASIS,
                None,
                None,
            ));
            decisions.push(Decision::admit(
                format!("normalization:{}:covariance:{row}:{column}", subject),
                "covariance_normalized",
                "packed covariance entry scaled by both canonically permuted member transforms",
                Some(literal_hash),
            ));
            if let (Some(source), Some(output)) = (&raw_correlation, &mut correlation) {
                output.push(source[old_index]);
            }
        }
    }
    Ok((
        JointStatistics::new(
            observation,
            raw.block_id.clone(),
            members,
            covariance,
            correlation,
        ),
        receipts,
        decisions,
    ))
}

fn parse_member_reference(
    raw: &str,
    subject: &str,
) -> Result<(String, LocalComponent), CompileError> {
    let fields: Vec<&str> = raw.split(':').collect();
    match fields.as_slice() {
        [claim, "scalar"] => Ok((
            require_identifier(claim, "joint claim id", subject)?,
            LocalComponent::Scalar,
        )),
        [claim, "x", knot] | [claim, "y", knot] => {
            let claim = require_identifier(claim, "joint claim id", subject)?;
            let knot = knot.parse::<u32>().map_err(|_| {
                CompileError::new(
                    "invalid_joint_member",
                    subject,
                    "curve member knot must be a canonical u32",
                )
            })?;
            if knot.to_string() != fields[2] {
                return Err(CompileError::new(
                    "invalid_joint_member",
                    subject,
                    "curve member knot must use canonical unsigned decimal spelling",
                ));
            }
            let component = if fields[1] == "x" {
                LocalComponent::CurveAbscissa(knot)
            } else {
                LocalComponent::CurveOrdinate(knot)
            };
            Ok((claim, component))
        }
        _ => Err(CompileError::new(
            "invalid_joint_member",
            subject,
            "member must be claim:scalar, claim:x:knot, or claim:y:knot",
        )),
    }
}

fn packed_index(row: usize, column: usize) -> usize {
    let (row, column) = if row >= column {
        (row, column)
    } else {
        (column, row)
    };
    row * (row + 1) / 2 + column
}

fn sort_decisions(decisions: &mut [Decision]) {
    decisions.sort_by(|left, right| {
        (
            left.subject.as_str(),
            left.reason_code,
            left.verdict,
            left.source_hash,
            left.detail.as_str(),
        )
            .cmp(&(
                right.subject.as_str(),
                right.reason_code,
                right.verdict,
                right.source_hash,
                right.detail.as_str(),
            ))
    });
}

fn render_decision(decision: &Decision) -> String {
    let case_id = hash_domain(DECISION_ID_DOMAIN, &decision.canonical_preimage()).to_hex();
    let source_hash = decision
        .source_hash
        .map_or_else(String::new, |hash| hash.to_hex());
    let pack_hash = decision
        .pack_hash
        .map_or_else(String::new, |hash| hash.to_hex());
    format!(
        "{{\"check\":\"matdb-pack\",\"compiler\":\"{}\",\"case_id\":\"{}\",\"subject\":\"{}\",\"verdict\":\"{}\",\
         \"reason_code\":\"{}\",\"source_hash\":\"{}\",\"pack_hash\":\"{}\",\"detail\":\"{}\"}}",
        COMPILER_ID,
        case_id,
        super::json_escape(&decision.subject),
        decision.verdict,
        decision.reason_code,
        source_hash,
        pack_hash,
        super::json_escape(&decision.detail),
    )
}

#[allow(clippy::too_many_lines)] // Publication barriers are intentionally explicit and sequential.
fn write_new_verified(
    output: &Path,
    bytes: &[u8],
    expected: ContentHash,
) -> Result<(), CompileError> {
    let parent = output.parent().unwrap_or_else(|| Path::new("."));
    let file_name = output.file_name().ok_or_else(|| {
        CompileError::new(
            "invalid_output_path",
            "output",
            "output must name a file within an existing directory",
        )
    })?;
    let metadata = fs::symlink_metadata(parent).map_err(|error| {
        CompileError::new(
            "output_unavailable",
            "output",
            format!("cannot inspect output directory: {error}"),
        )
    })?;
    if !metadata.file_type().is_dir() || metadata.file_type().is_symlink() {
        return Err(CompileError::new(
            "invalid_output_directory",
            "output",
            "output parent must be an existing non-symlink directory",
        ));
    }
    let canonical_parent = parent.canonicalize().map_err(|error| {
        CompileError::new(
            "output_unavailable",
            "output",
            format!("cannot resolve output directory: {error}"),
        )
    })?;
    let canonical_output = canonical_parent.join(file_name);
    match fs::symlink_metadata(&canonical_output) {
        Ok(_) => {
            return Err(CompileError::new(
                "output_exists",
                "output",
                "refusing to replace an existing file, directory, or symlink",
            ));
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => {
            return Err(CompileError::new(
                "output_unavailable",
                "output",
                format!("cannot inspect output path: {error}"),
            ));
        }
    }

    // The verified candidate deliberately remains as a same-inode retained
    // evidence link. This obeys the repository's no-deletion rule while
    // allowing hard_link to publish the final name atomically and no-clobber.
    let mut staged = None;
    for sequence in 0..1_024u16 {
        let candidate = canonical_parent.join(format!(
            ".frankensim-matdb-{}-{}-{sequence}.verified",
            expected,
            std::process::id()
        ));
        if candidate == canonical_output {
            continue;
        }
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&candidate)
        {
            Ok(file) => {
                staged = Some((candidate, file));
                break;
            }
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(error) => {
                return Err(CompileError::new(
                    "output_stage_failed",
                    "output",
                    format!("cannot create same-directory candidate: {error}"),
                ));
            }
        }
    }
    let (stage_path, mut file) = staged.ok_or_else(|| {
        CompileError::new(
            "output_stage_failed",
            "output",
            "exhausted the bounded same-directory candidate namespace",
        )
    })?;
    file.write_all(bytes).map_err(|error| {
        CompileError::new(
            "output_stage_failed",
            "output",
            format!("cannot write complete candidate: {error}"),
        )
    })?;
    file.sync_all().map_err(|error| {
        CompileError::new(
            "output_stage_failed",
            "output",
            format!("cannot sync complete candidate: {error}"),
        )
    })?;
    drop(file);
    let staged_bytes = read_bounded_regular(&stage_path, bytes.len(), "output-stage")?;
    if staged_bytes != bytes {
        return Err(CompileError::new(
            "output_validation_failed",
            "output",
            "candidate bytes differ from the verified compiler artifact",
        ));
    }
    NormalizedPack::from_bytes_verified(expected, &staged_bytes).map_err(|error| {
        CompileError::new(
            "output_validation_failed",
            "output",
            format!("candidate failed runtime admission: {error}"),
        )
    })?;
    fs::hard_link(&stage_path, &canonical_output).map_err(|error| {
        let (code, detail) = if error.kind() == std::io::ErrorKind::AlreadyExists {
            (
                "output_exists",
                "refusing to replace an output created during compilation".to_string(),
            )
        } else {
            (
                "output_publish_failed",
                format!("atomic no-clobber publication failed: {error}"),
            )
        };
        CompileError::new(code, "output", detail)
    })?;
    let published = read_bounded_regular(&canonical_output, bytes.len(), "output")?;
    if published != bytes {
        return Err(CompileError::new(
            "output_validation_failed",
            "output",
            "published hard link does not retain the verified candidate bytes",
        ));
    }
    NormalizedPack::from_bytes_verified(expected, &published).map_err(|error| {
        CompileError::new(
            "output_validation_failed",
            "output",
            format!("published output failed runtime admission: {error}"),
        )
    })?;
    File::open(&canonical_parent)
        .and_then(|directory| directory.sync_all())
        .map_err(|error| {
            CompileError::new(
                "output_sync_failed",
                "output",
                format!("cannot sync output directory: {error}"),
            )
        })?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicU64, Ordering};

    use super::*;

    static NEXT_FIXTURE: AtomicU64 = AtomicU64::new(0);

    const SOURCE: &str = concat!(
        "frankensim.matdb-source.v1\n",
        "observation\tcoupon\talloy-X-solution-treated\tASTM-fixture\tjoint coupon series\n",
        "scalar\tdensity\tcoupon\tdensity\t7.85\tg/cm3\tconstant\n",
        "uncertainty\tdensity\tabsolute\t0.005\tg/cm3\t0.95\t1\n",
        "validity\tdensity\ttemperature\t0\t100\tdegC\n",
        "curve\tmodulus\tcoupon\tyoung_modulus\ttemperature\tdegC\tGPa\t0:210,100:202\tlinear\n",
        "uncertainty\tmodulus\trelative\t2\t%\t0.95\t1\n",
        "validity\tmodulus\ttemperature\t0\t100\tdegC\n",
        "frame\tmodulus\tspecimen\tlab\n",
        "joint\tcoupon\tdensity-modulus\tdensity:scalar,modulus:y:0\t0.000025,0,0.000009\t1,0,1\t1\n",
    );

    fn fixture_dir() -> PathBuf {
        loop {
            let sequence = NEXT_FIXTURE.fetch_add(1, Ordering::Relaxed);
            let path = std::env::temp_dir().join(format!(
                "frankensim-matdb-pack-test-{}-{sequence}",
                std::process::id()
            ));
            match fs::create_dir(&path) {
                Ok(()) => return path,
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
                Err(error) => panic!("unique fixture directory: {error}"),
            }
        }
    }

    fn manifest(profile: &str, redistribution: bool) -> String {
        let redistribution = if redistribution {
            "redistribution\tpermitted\tCC-BY-4.0 redistribution with attribution\n"
        } else {
            ""
        };
        format!(
            "{MANIFEST_HEADER}\n\
             pack_id\tfixture-alloy-x\n\
             {redistribution}\
             citation\tfixture handbook table 7\n\
             license\tCC-BY-4.0\n\
             source\tprimary\tsource.tsv\t{profile}\n"
        )
    }

    fn write_fixture(manifest_text: &str, source_text: &str) -> (PathBuf, PathBuf) {
        let directory = fixture_dir();
        let manifest_path = directory.join("manifest.tsv");
        let source_path = directory.join("source.tsv");
        fs::write(&manifest_path, manifest_text).expect("manifest fixture");
        fs::write(&source_path, source_text).expect("source fixture");
        (manifest_path, directory)
    }

    #[test]
    fn compiles_scalar_curve_covariance_and_is_two_run_deterministic() {
        let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), SOURCE);
        let first = compile_manifest(&manifest_path).expect("fixture compiles");
        let second = compile_manifest(&manifest_path).expect("independent repeat compiles");

        assert_eq!(first.bytes, second.bytes);
        assert_eq!(first.decisions, second.decisions);
        assert_eq!(first.pack.content_hash(), second.pack.content_hash());
        let decoded = NormalizedPack::from_bytes_verified(first.pack.content_hash(), &first.bytes)
            .expect("compiler output is runtime-loadable");
        assert_eq!(decoded.to_bytes(), first.bytes);
        assert_eq!(decoded.claims().claim_count(), 2);
        assert_eq!(decoded.joint_statistics().len(), 1);
        let (modulus_id, modulus) = decoded.claims().claims_for("young_modulus")[0];
        let PropertyValue::Curve { knots, .. } = &modulus.value else {
            panic!("young_modulus remains a normalized curve");
        };
        assert_eq!(knots.as_slice(), &[(273.15, 210.0e9), (373.15, 202.0e9)]);
        assert_eq!(
            modulus.validity.bound("temperature"),
            Some((273.15, 373.15))
        );

        let temperature_value_receipts: Vec<_> = decoded
            .normalizations()
            .iter()
            .filter(|receipt| {
                matches!(
                    receipt.target(),
                    NormalizationTarget::ClaimValue(member)
                        if member.claim() == modulus_id && receipt.source_basis() == "degC"
                )
            })
            .collect();
        assert_eq!(temperature_value_receipts.len(), 2);
        for receipt in temperature_value_receipts {
            assert_eq!(receipt.scale(), 1.0);
            assert_eq!(receipt.offset(), 273.15);
            assert_eq!(receipt.source_frame(), Some("specimen"));
            assert_eq!(receipt.target_frame(), Some("lab"));
        }

        let temperature_validity_receipts: Vec<_> = decoded
            .normalizations()
            .iter()
            .filter(|receipt| {
                matches!(
                    receipt.target(),
                    NormalizationTarget::ValidityBound { claim, axis, .. }
                        if *claim == modulus_id && axis == "temperature"
                )
            })
            .collect();
        assert_eq!(temperature_validity_receipts.len(), 2);
        for receipt in temperature_validity_receipts {
            assert_eq!(receipt.source_basis(), "degC");
            assert_eq!(receipt.scale(), 1.0);
            assert_eq!(receipt.offset(), 273.15);
            assert_eq!(receipt.source_frame(), None);
            assert_eq!(receipt.target_frame(), None);
        }

        let modulus_uncertainty_receipt = decoded
            .normalizations()
            .iter()
            .find(|receipt| {
                matches!(
                    receipt.target(),
                    NormalizationTarget::ClaimUncertainty { claim } if *claim == modulus_id
                )
            })
            .expect("modulus uncertainty receipt");
        assert_eq!(modulus_uncertainty_receipt.offset(), 0.0);
        assert_eq!(modulus_uncertainty_receipt.source_frame(), Some("specimen"));
        assert_eq!(modulus_uncertainty_receipt.target_frame(), Some("lab"));
        assert!(decoded.normalizations().iter().all(|receipt| {
            !matches!(
                receipt.target(),
                NormalizationTarget::JointCovariance { .. }
            ) || (receipt.offset() == 0.0
                && receipt.source_frame().is_none()
                && receipt.target_frame().is_none())
        }));

        let joint = &decoded.joint_statistics()[0];
        assert_eq!(joint.members().len(), 2);
        assert_eq!(joint.correlation().expect("source correlation")[1], 0.0);
        for (row, member) in joint.members().iter().copied().enumerate() {
            let claim = decoded
                .claims()
                .claim(member.claim())
                .expect("joint member claim");
            let diagonal = joint.covariance()[packed_index(row, row)];
            let receipt = decoded
                .normalizations()
                .iter()
                .find(|receipt| {
                    matches!(
                        receipt.target(),
                        NormalizationTarget::JointCovariance {
                            row: target_row,
                            column,
                            ..
                        } if usize::try_from(*target_row).ok() == Some(row)
                            && *column == *target_row
                    )
                })
                .expect("diagonal covariance receipt");
            match claim.key.name() {
                "density" => {
                    assert!((diagonal - 25.0).abs() <= 25.0 * f64::EPSILON);
                    assert_eq!(receipt.dims(), Dims([-6, 2, 0, 0, 0, 0]));
                    assert!((receipt.scale() - 1.0e6).abs() <= 1.0e6 * f64::EPSILON);
                }
                "young_modulus" => {
                    assert!((diagonal - 9.0e12).abs() <= 9.0e12 * f64::EPSILON);
                    assert_eq!(receipt.dims(), Dims([-2, 2, -4, 0, 0, 0]));
                    assert!((receipt.scale() - 1.0e18).abs() <= 1.0e18 * f64::EPSILON);
                }
                other => panic!("unexpected covariance member property {other}"),
            }
        }
        assert_eq!(decoded.normalizations().len(), 14);
        assert!(first.decisions.iter().all(|decision| {
            decision.verdict == "admit" && decision.pack_hash == Some(first.pack.content_hash())
        }));
    }

    #[test]
    fn multi_observation_claims_are_canonicalized_by_content_id() {
        let source = SOURCE
            .replacen(
                "observation\tcoupon\talloy-X-solution-treated\tASTM-fixture\tjoint coupon series\n",
                concat!(
                    "observation\tcoupon\talloy-X-solution-treated\tASTM-fixture\tjoint coupon series\n",
                    "observation\tsecond\talloy-X-retest\tASTM-fixture\trepeat coupon series\n",
                ),
                1,
            )
            .replacen(
                "scalar\tdensity\tcoupon\tdensity",
                "scalar\tdensity\tcoupon,second\tdensity",
                1,
            );
        let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &source);
        let compiled = compile_manifest(&manifest_path).expect("two observations compile");
        let density = compiled.pack.claims().claims_for("density")[0].1;

        assert_eq!(density.observations.len(), 2);
        assert!(
            density
                .observations
                .windows(2)
                .all(|pair| pair[0] < pair[1])
        );
    }

    #[test]
    fn absent_redistribution_terms_refuse_before_source_loading() {
        let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, false), SOURCE);
        let error = compile_manifest(&manifest_path).expect_err("terms are load-bearing");

        assert_eq!(error.code, "missing_redistribution_terms");
        assert_eq!(error.subject, "manifest");
    }

    #[test]
    fn missing_license_refuses_before_source_loading() {
        let manifest = manifest(MATERIAL_PROFILE, true).replacen("license\tCC-BY-4.0\n", "", 1);
        let (manifest_path, _) = write_fixture(&manifest, SOURCE);
        let error = compile_manifest(&manifest_path).expect_err("license is load-bearing");

        assert_eq!(error.code, "missing_license");
        assert_eq!(error.subject, "manifest");
    }

    #[test]
    fn blank_unit_refuses_instead_of_inference() {
        let source = SOURCE.replacen(
            "scalar\tdensity\tcoupon\tdensity\t7.85\tg/cm3\tconstant",
            "scalar\tdensity\tcoupon\tdensity\t7.85\t\tconstant",
            1,
        );
        let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &source);
        let error = compile_manifest(&manifest_path).expect_err("unit inference is forbidden");

        assert_eq!(error.code, "ambiguous_unit");
        assert!(error.subject.contains("source:primary:line:"));
    }

    #[test]
    fn numeric_continuation_cannot_masquerade_as_a_unit() {
        let source = SOURCE.replacen(
            "scalar\tdensity\tcoupon\tdensity\t7.85\tg/cm3\tconstant",
            "scalar\tdensity\tcoupon\tdensity\t7.85\te3\tconstant",
            1,
        );
        let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &source);
        let error = compile_manifest(&manifest_path).expect_err("e3 is not a unit token");

        assert_eq!(error.code, "invalid_unit");
    }

    #[test]
    fn confidence_requires_the_unconverted_dimensionless_basis() {
        let source = SOURCE.replacen(
            "uncertainty\tdensity\tabsolute\t0.005\tg/cm3\t0.95\t1",
            "uncertainty\tdensity\tabsolute\t0.005\tg/cm3\t0.95\t%",
            1,
        );
        let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &source);
        let error = compile_manifest(&manifest_path).expect_err("confidence has no receipt target");

        assert_eq!(error.code, "invalid_confidence_unit");
    }

    #[test]
    fn unevaluable_value_interpolation_pairs_refuse() {
        let scalar = SOURCE.replacen(
            "scalar\tdensity\tcoupon\tdensity\t7.85\tg/cm3\tconstant",
            "scalar\tdensity\tcoupon\tdensity\t7.85\tg/cm3\tlinear",
            1,
        );
        let curve = SOURCE.replacen(
            "curve\tmodulus\tcoupon\tyoung_modulus\ttemperature\tdegC\tGPa\t0:210,100:202\tlinear",
            "curve\tmodulus\tcoupon\tyoung_modulus\ttemperature\tdegC\tGPa\t0:210,100:202\tconstant",
            1,
        );
        for source in [scalar, curve] {
            let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &source);
            let error =
                compile_manifest(&manifest_path).expect_err("claim could never be evaluated");
            assert_eq!(error.code, "incompatible_interpolation");
        }
    }

    #[test]
    fn decimal_underflow_refuses_before_signed_zero_can_be_laundered() {
        for number in ["1e-999", "-1e-999"] {
            let source = SOURCE.replacen(
                "scalar\tdensity\tcoupon\tdensity\t7.85\tg/cm3\tconstant",
                &format!("scalar\tdensity\tcoupon\tdensity\t{number}\t1\tconstant"),
                1,
            );
            let (manifest_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &source);
            let error =
                compile_manifest(&manifest_path).expect_err("nonzero significand underflows");
            assert_eq!(error.code, "numeric_underflow");
        }
    }

    #[test]
    fn unsupported_species_profile_is_a_typed_refusal() {
        let (manifest_path, _) = write_fixture(&manifest("nasa9-v1", true), SOURCE);
        let error =
            compile_manifest(&manifest_path).expect_err("no runtime species codec exists yet");

        assert_eq!(error.code, "unsupported_source_profile");
        assert_eq!(error.subject, "source:primary");
    }

    #[test]
    fn covariance_member_permutation_preserves_the_matrix_semantics() {
        let density = ClaimId(ContentHash([1; 32]));
        let modulus = ClaimId(ContentHash([2; 32]));
        let observation = ObservationId(ContentHash([3; 32]));
        let mut components = BTreeMap::new();
        components.insert(
            ("density".to_string(), LocalComponent::Scalar),
            ComponentMeta {
                member: StatisticMember::scalar(density),
                dims: Dims([-3, 1, 0, 0, 0, 0]),
                scale: 1_000.0,
                unit: "g/cm3".to_string(),
                observations: BTreeSet::from([observation]),
            },
        );
        components.insert(
            ("modulus".to_string(), LocalComponent::CurveAbscissa(0)),
            ComponentMeta {
                member: StatisticMember::curve_abscissa(modulus, 0),
                dims: Dims([0, 0, 0, 1, 0, 0]),
                scale: 1.0,
                unit: "degC".to_string(),
                observations: BTreeSet::from([observation]),
            },
        );
        components.insert(
            ("modulus".to_string(), LocalComponent::CurveOrdinate(0)),
            ComponentMeta {
                member: StatisticMember::curve_ordinate(modulus, 0),
                dims: Dims([-1, 1, -2, 0, 0, 0]),
                scale: 1.0e9,
                unit: "GPa".to_string(),
                observations: BTreeSet::from([observation]),
            },
        );
        let first = RawJoint {
            observation: "coupon".to_string(),
            block_id: "block".to_string(),
            members: vec![
                "density:scalar".to_string(),
                "modulus:x:0".to_string(),
                "modulus:y:0".to_string(),
            ],
            covariance: vec![
                "4".to_string(),
                "1".to_string(),
                "9".to_string(),
                "2".to_string(),
                "3".to_string(),
                "16".to_string(),
            ],
            correlation: Some(vec![
                "1".to_string(),
                "0.1".to_string(),
                "1".to_string(),
                "0.2".to_string(),
                "0.3".to_string(),
                "1".to_string(),
            ]),
            source_hash: ContentHash([4; 32]),
        };
        let second = RawJoint {
            observation: "coupon".to_string(),
            block_id: "block".to_string(),
            members: vec![
                "modulus:y:0".to_string(),
                "density:scalar".to_string(),
                "modulus:x:0".to_string(),
            ],
            covariance: vec![
                "16".to_string(),
                "2".to_string(),
                "4".to_string(),
                "3".to_string(),
                "1".to_string(),
                "9".to_string(),
            ],
            correlation: Some(vec![
                "1".to_string(),
                "0.2".to_string(),
                "1".to_string(),
                "0.3".to_string(),
                "0.1".to_string(),
                "1".to_string(),
            ]),
            source_hash: ContentHash([4; 32]),
        };
        let (first, _, _) =
            compile_joint(&first, observation, &components).expect("first member order");
        let (second, _, _) =
            compile_joint(&second, observation, &components).expect("permuted member order");

        assert_eq!(first.members(), second.members());
        assert_eq!(first.covariance(), second.covariance());
        assert_eq!(first.correlation(), second.correlation());
    }

    #[test]
    fn publication_never_replaces_an_existing_path() {
        let (manifest_path, directory) = write_fixture(&manifest(MATERIAL_PROFILE, true), SOURCE);
        let compiled = compile_manifest(&manifest_path).expect("fixture compiles");
        let output = directory.join("pack.fsmatpk");
        fs::write(&output, b"sentinel").expect("preexisting output");

        let error = write_new_verified(&output, &compiled.bytes, compiled.pack.content_hash())
            .expect_err("replacement is forbidden");
        assert_eq!(error.code, "output_exists");
        assert_eq!(fs::read(output).expect("sentinel retained"), b"sentinel");
    }

    #[test]
    fn publication_exposes_only_a_complete_runtime_verified_hard_link() {
        let (manifest_path, directory) = write_fixture(&manifest(MATERIAL_PROFILE, true), SOURCE);
        let compiled = compile_manifest(&manifest_path).expect("fixture compiles");
        let output = directory.join("pack.fsmatpk");

        write_new_verified(&output, &compiled.bytes, compiled.pack.content_hash())
            .expect("verified no-clobber publication");
        let retained = read_bounded_regular(&output, compiled.bytes.len(), "test-output")
            .expect("published output");
        assert_eq!(retained, compiled.bytes);
        assert_eq!(
            NormalizedPack::from_bytes_verified(compiled.pack.content_hash(), &retained)
                .expect("published runtime artifact"),
            compiled.pack
        );
    }

    #[test]
    fn caller_cannot_alias_the_hidden_candidate_namespace() {
        let (manifest_path, directory) = write_fixture(&manifest(MATERIAL_PROFILE, true), SOURCE);
        let compiled = compile_manifest(&manifest_path).expect("fixture compiles");
        let output = directory.join(format!(
            ".frankensim-matdb-{}-{}-0.verified",
            compiled.pack.content_hash(),
            std::process::id()
        ));

        write_new_verified(&output, &compiled.bytes, compiled.pack.content_hash())
            .expect("candidate alias is skipped");
        assert_eq!(
            read_bounded_regular(&output, compiled.bytes.len(), "aliased-output")
                .expect("published aliased output"),
            compiled.bytes
        );
    }

    #[test]
    fn receipt_budget_refuses_before_normalized_allocations() {
        let mut raw = RawDatabase::default();
        raw.claims.insert(
            "oversized".to_string(),
            RawClaim {
                id: "oversized".to_string(),
                observations: vec!["observation".to_string()],
                property: "oversized".to_string(),
                value: RawClaimValue::Curve {
                    abscissa: "x".to_string(),
                    abscissa_unit: "1".to_string(),
                    ordinate_unit: "1".to_string(),
                    points: vec![("0".to_string(), "0".to_string()); 50_001],
                },
                interpolation: InterpolationPolicy::TabulatedOnly,
                source_id: "source".to_string(),
                source_hash: ContentHash([1; 32]),
                record_hash: ContentHash([2; 32]),
            },
        );

        let error =
            preflight_normalization_budget(&raw).expect_err("two receipts per knot exceed cap");
        assert_eq!(error.code, "resource_limit");
        assert_eq!(error.subject, "normalizations");
    }

    #[test]
    fn repeated_provenance_budget_refuses_before_claim_set_allocations() {
        let manifest = Manifest {
            pack_id: "fixture".to_string(),
            redistribution_terms: "permitted".to_string(),
            citation: "x".repeat(1_048_000),
            license: "MIT".to_string(),
            sources: Vec::new(),
        };
        let mut raw = RawDatabase::default();
        for index in 0..64 {
            let id = format!("observation-{index}");
            raw.observations.insert(
                id.clone(),
                RawObservation {
                    id,
                    specimen: String::new(),
                    method: String::new(),
                    caveats: String::new(),
                    source_id: "source".to_string(),
                    source_hash: ContentHash([1; 32]),
                    record_hash: ContentHash([2; 32]),
                },
            );
        }
        preflight_provenance_budget(&manifest, &raw).expect("64 records remain under the cap");
        raw.observations.insert(
            "observation-64".to_string(),
            RawObservation {
                id: "observation-64".to_string(),
                specimen: String::new(),
                method: String::new(),
                caveats: String::new(),
                source_id: "source".to_string(),
                source_hash: ContentHash([1; 32]),
                record_hash: ContentHash([2; 32]),
            },
        );

        let error = preflight_provenance_budget(&manifest, &raw)
            .expect_err("65 repeated megabyte-scale citations exceed the cap");
        assert_eq!(error.code, "resource_limit");
        assert_eq!(error.subject, "provenance");
    }

    #[test]
    fn invalid_utf8_refusals_are_bound_to_the_offending_source_bytes() {
        let mut hashes = Vec::new();
        for invalid in [0xff_u8, 0xfe] {
            let directory = fixture_dir();
            let manifest_path = directory.join("manifest.tsv");
            let source_path = directory.join("source.tsv");
            fs::write(&manifest_path, manifest(MATERIAL_PROFILE, true)).expect("manifest fixture");
            let mut source = format!("{SOURCE_HEADER}\n").into_bytes();
            source.push(invalid);
            fs::write(&source_path, &source).expect("invalid UTF-8 source fixture");

            let error = compile_manifest(&manifest_path).expect_err("invalid UTF-8 refuses");
            assert_eq!(error.code, "invalid_source_encoding");
            let source_hash = error.input_hash.expect("offending byte hash");
            assert_eq!(source_hash, hash_domain(SOURCE_FILE_DOMAIN, &source));
            hashes.push(source_hash);
        }

        assert_ne!(hashes[0], hashes[1]);
    }

    #[test]
    fn refusal_evidence_is_bound_to_the_rejected_input_envelope() {
        let first_source = SOURCE.replacen("\t7.85\tg/cm3\t", "\t7.85\t\t", 1);
        let second_source = SOURCE.replacen("\t7.85\tg/cm3\t", "\t8.00\t\t", 1);
        let (first_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &first_source);
        let (second_path, _) = write_fixture(&manifest(MATERIAL_PROFILE, true), &second_source);
        let first_error = compile_manifest(&first_path).expect_err("blank unit refuses");
        let second_error = compile_manifest(&second_path).expect_err("blank unit refuses");
        let first_hash = first_error.input_hash.expect("first envelope hash");
        let second_hash = second_error.input_hash.expect("second envelope hash");

        assert_ne!(first_hash, second_hash);
        let first = Decision::refusal(&first_error, Some(first_hash));
        let second = Decision::refusal(&second_error, Some(second_hash));
        assert_ne!(first.canonical_preimage(), second.canonical_preimage());
        assert!(render_decision(&first).contains(&first_hash.to_hex()));
    }

    #[test]
    fn decision_identity_distinguishes_source_and_pack_hash_roles() {
        let error = CompileError::new("fixture_refusal", "fixture", "typed refusal");
        let hash = ContentHash([9; 32]);
        let source = Decision::refusal(&error, Some(hash));
        let pack = Decision::refusal_for_pack(&error, hash);

        assert_ne!(source.canonical_preimage(), pack.canonical_preimage());
        assert_ne!(
            hash_domain(DECISION_ID_DOMAIN, &source.canonical_preimage()),
            hash_domain(DECISION_ID_DOMAIN, &pack.canonical_preimage())
        );
    }

    #[test]
    fn decision_rows_are_stable_json_lines_without_host_paths() {
        let (manifest_path, directory) = write_fixture(&manifest(MATERIAL_PROFILE, true), SOURCE);
        let compiled = compile_manifest(&manifest_path).expect("fixture compiles");
        let rows: Vec<String> = compiled.decisions.iter().map(render_decision).collect();

        assert_eq!(
            rows,
            compiled
                .decisions
                .iter()
                .map(render_decision)
                .collect::<Vec<_>>()
        );
        assert!(
            rows.iter()
                .all(|row| row.starts_with("{\"check\":\"matdb-pack\""))
        );
        assert!(
            rows.iter()
                .all(|row| !row.contains(&directory.display().to_string()))
        );
        assert!(rows.iter().all(|row| row.contains("\"pack_hash\":\"")));
    }
}
