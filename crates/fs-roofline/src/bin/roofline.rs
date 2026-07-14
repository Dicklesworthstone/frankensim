//! Roofline harness CLI (plan §14.4 nightly lane).
//!
//! Usage:
//!   roofline [--n <elements>] [--warmup <k>] [--reps <k>] [--ledger <db>]
//!            [--baseline <jsonl>] [--firmware <identity>]
//!            [--authority-policy <tsv>] [--retained-receipts <txt>]
//!   roofline promote --store <jsonl> --firmware <identity>
//!            --operator <name> --justification <text>
//!            [--probes <k≥3>] [--age-days <d>]
//!
//! Probes the machine axes, runs the default kernel registry, prints one
//! JSON line per kernel (plus the axes line and the §14.1 coverage table),
//! and — when `--ledger` is given — records the run as ledger provenance
//! and reports staleness for every registered kernel.

use fs_roofline::authority::ConfiguredPromotionAuthority;
use fs_roofline::production::{
    ProductionProbe, ProductionRun, ProductionRunConfig, ReportOnlyProductionRun,
};
use fs_roofline::{
    AttestedAxisBaselinePolicy, AttestedBaselineStore, AxisBaselinePolicy, BaselineAxes,
    BaselineIdentity, BaselineStore, CUSTOM_REGISTRY_PROTOCOL_VERSION, ContentHash, MachineAxes,
    PRODUCTION_PROTOCOL_VERSION, SECTION_14_1_TARGETS, STALENESS_MAX_AGE_NS, days_since_epoch_now,
    staleness,
};
use std::collections::BTreeSet;
use std::ffi::OsString;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

#[cfg(unix)]
use std::os::unix::fs::MetadataExt as _;

fn json_escape(value: &str) -> String {
    use core::fmt::Write as _;

    let mut escaped = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            c if c.is_control() => {
                let _ = write!(escaped, "\\u{:04x}", u32::from(c));
            }
            c => escaped.push(c),
        }
    }
    escaped
}

fn evidence_admission_json(citation_eligible: bool, refusal: Option<&str>) -> String {
    let reason = refusal.map_or_else(
        || "\"not_recorded\"".to_string(),
        |reason| format!("\"{}\"", json_escape(reason)),
    );
    format!(
        "{{\"schema\":\"fs-roofline-evidence-admission-v2\",\"citation_eligible\":{citation_eligible},\"recorded\":false,\"citable\":false,\"reason\":{reason}}}"
    )
}

fn fail(detail: &str) -> std::process::ExitCode {
    eprintln!(
        "{{\"error\":\"Roofline\",\"detail\":\"{}\"}}",
        json_escape(detail)
    );
    std::process::ExitCode::FAILURE
}

#[derive(Debug, PartialEq, Eq)]
struct CliArgs {
    n: usize,
    warmup: usize,
    reps: usize,
    ledger_path: Option<String>,
    baseline_path: Option<String>,
    firmware: Option<String>,
    authority_policy_path: Option<String>,
    retained_receipts_path: Option<String>,
}

impl Default for CliArgs {
    fn default() -> Self {
        Self {
            n: 1 << 22, // 32 MiB per f64 buffer: streams past every L2/L3
            warmup: 2,
            reps: 9,
            ledger_path: None,
            baseline_path: None,
            firmware: None,
            authority_policy_path: None,
            retained_receipts_path: None,
        }
    }
}

fn positive_usize(flag: &str, value: &str) -> Result<usize, String> {
    value
        .parse::<usize>()
        .ok()
        .filter(|value| *value > 0)
        .ok_or_else(|| format!("{flag} must be a positive integer"))
}

const MAX_PROMOTION_PROBES: usize = 1_000;
const MAX_BASELINE_AGE_DAYS: u32 = 36_500;
const MAX_AUTHORITY_INPUT_BYTES: usize =
    fs_roofline::authority::MAX_PROMOTION_AUTHORITY_POLICY_BYTES;
const MAX_RETAINED_RECEIPTS_INPUT_BYTES: usize = fs_roofline::baseline::MAX_BASELINE_STORE_BYTES;

/// `roofline promote` — the operator bootstrap for governed baselines
/// (bead c40j): probe the machine axes N ≥ 3 times, build candidates,
/// run [`fs_roofline::promote_baseline`] (which REFUSES on a loaded
/// host — the drift bands are the point), and create-or-update the
/// JSONL store. This creates an operator-promoted candidate; attestation
/// and promotion-authority admission are separate operations.
struct PromoteArgs {
    store: String,
    firmware: String,
    operator: String,
    justification: String,
    probes: usize,
    age_days: u32,
}

fn parse_promote_args(args: &[String]) -> Result<PromoteArgs, String> {
    let (mut store, mut firmware, mut operator, mut justification) = (None, None, None, None);
    let mut probes = 3usize;
    let mut age_days = 90u32;
    let mut seen = std::collections::BTreeSet::new();
    let mut args = args.iter().skip(2);
    while let Some(flag) = args.next() {
        if !matches!(
            flag.as_str(),
            "--store" | "--firmware" | "--operator" | "--justification" | "--probes" | "--age-days"
        ) {
            return Err(format!("unknown promote argument {flag:?}"));
        }
        if !seen.insert(flag.as_str()) {
            return Err(format!("duplicate promote argument {flag:?}"));
        }
        let value = args
            .next()
            .filter(|value| !value.starts_with("--"))
            .ok_or_else(|| format!("{flag} requires a value"))?;
        if value.is_empty() {
            return Err(format!("{flag} requires a non-empty value"));
        }
        match flag.as_str() {
            "--store" => store = Some(value.clone()),
            "--firmware" => firmware = Some(value.clone()),
            "--operator" => operator = Some(value.clone()),
            "--justification" => justification = Some(value.clone()),
            "--probes" => probes = positive_usize(flag, value)?,
            "--age-days" => {
                age_days = value
                    .parse::<u32>()
                    .ok()
                    .filter(|v| *v > 0)
                    .ok_or_else(|| format!("{flag} must be a positive integer"))?;
            }
            _ => unreachable!("flag list checked above"),
        }
    }
    if probes < 3 {
        return Err(
            "--probes must be at least 3 (governed promotion needs mutual agreement)".to_string(),
        );
    }
    if probes > MAX_PROMOTION_PROBES {
        return Err(format!(
            "--probes must be at most {MAX_PROMOTION_PROBES}, got {probes}"
        ));
    }
    if age_days > MAX_BASELINE_AGE_DAYS {
        return Err(format!(
            "--age-days must be at most {MAX_BASELINE_AGE_DAYS}, got {age_days}"
        ));
    }
    Ok(PromoteArgs {
        store: store.ok_or("promote requires --store <jsonl>")?,
        firmware: firmware.ok_or("promote requires --firmware <identity>")?,
        operator: operator.ok_or("promote requires --operator <name>")?,
        justification: justification.ok_or("promote requires --justification <text>")?,
        probes,
        age_days,
    })
}

fn run_promote(args: &PromoteArgs) -> Result<(), String> {
    use fs_roofline::{BaselineCandidate, promote_baseline};
    let mut candidates = Vec::with_capacity(args.probes);
    for ordinal in 0..args.probes {
        let axes = MachineAxes::probe();
        println!("{}", axes.to_jsonl());
        let identity = BaselineIdentity::current(&axes, args.firmware.clone())
            .map_err(|error| format!("probe {ordinal}: {error}"))?;
        // A content-derived source receipt: the probe's own canonical bytes
        // under a CLI-specific domain. Promotion remains a plain candidate;
        // attestation and configured authority admission are separate steps.
        let receipt = fs_blake3::hash_domain(
            "fs-roofline.cli-baseline-source.v1",
            axes.to_jsonl().as_bytes(),
        );
        let candidate = BaselineCandidate::from_receipt(axes, identity, receipt)
            .map_err(|error| format!("probe {ordinal}: {error}"))?;
        candidates.push(candidate);
    }
    let now_day = days_since_epoch_now().map_err(|error| error.to_string())?;
    let baseline = promote_baseline(
        &candidates,
        args.operator.clone(),
        args.justification.clone(),
        now_day,
        args.age_days,
    )
    .map_err(|error| format!("promotion refused: {error}"))?;
    update_promoted_store(Path::new(&args.store), baseline.clone())?;
    println!("{}", baseline.to_jsonl());
    println!(
        "{{\"promote\":\"ok\",\"fingerprint\":\"{:016x}\",\"store\":\"{}\",\"probes\":{},\"operator\":\"{}\"}}",
        baseline.identity().fingerprint(),
        json_escape(&args.store),
        args.probes,
        json_escape(&args.operator)
    );
    Ok(())
}

fn sidecar_path(store: &Path, suffix: &str) -> Result<PathBuf, String> {
    let parent = store
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let file_name = store
        .file_name()
        .ok_or_else(|| format!("baseline store path {} has no file name", store.display()))?;
    let mut sidecar = OsString::from(".");
    sidecar.push(file_name);
    sidecar.push(suffix);
    Ok(parent.join(sidecar))
}

fn promotion_lock_path(store: &Path) -> Result<PathBuf, String> {
    let parent = store
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or(Path::new("."))
        .canonicalize()
        .map_err(|error| format!("cannot resolve baseline-store directory: {error}"))?;
    let file_name = store
        .file_name()
        .ok_or_else(|| format!("baseline store path {} has no file name", store.display()))?;
    // `Debug` is lossless for the platform OsStr and avoids two non-UTF paths
    // aliasing through a lossy display conversion.
    #[allow(clippy::unnecessary_debug_formatting)]
    let identity = format!("{:?}/{file_name:?}", parent.as_os_str());
    let digest = fs_blake3::hash_domain(
        "fs-roofline.cli-baseline-store-lock.v1",
        identity.as_bytes(),
    );
    Ok(std::env::temp_dir().join(format!("fs-roofline-baseline-{digest}.lock")))
}

#[cfg(unix)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PromotionFileIdentity {
    device: u64,
    inode: u64,
    links: u64,
    len: u64,
    mode: u32,
    modified_seconds: i64,
    modified_nanoseconds: i64,
    changed_seconds: i64,
    changed_nanoseconds: i64,
}

#[cfg(unix)]
fn promotion_file_identity(
    path: &Path,
    metadata: &std::fs::Metadata,
) -> Result<PromotionFileIdentity, String> {
    if !metadata.file_type().is_file() {
        return Err(format!("{} must be a regular file", path.display()));
    }
    Ok(PromotionFileIdentity {
        device: metadata.dev(),
        inode: metadata.ino(),
        links: metadata.nlink(),
        len: metadata.len(),
        mode: metadata.mode(),
        modified_seconds: metadata.mtime(),
        modified_nanoseconds: metadata.mtime_nsec(),
        changed_seconds: metadata.ctime(),
        changed_nanoseconds: metadata.ctime_nsec(),
    })
}

#[cfg(unix)]
fn validate_promotion_path_identity(
    path: &Path,
    expected: PromotionFileIdentity,
    purpose: &str,
) -> Result<(), String> {
    let metadata = std::fs::symlink_metadata(path)
        .map_err(|error| format!("cannot re-inspect {purpose} {}: {error}", path.display()))?;
    let observed = promotion_file_identity(path, &metadata)?;
    if observed == expected {
        Ok(())
    } else {
        Err(format!(
            "{purpose} {} changed during promotion: expected {expected:?}, observed {observed:?}",
            path.display()
        ))
    }
}

#[cfg(unix)]
struct OpenedPromotionStore {
    file: std::fs::File,
    identity: PromotionFileIdentity,
    permissions: std::fs::Permissions,
}

#[cfg(unix)]
fn open_promotion_store(path: &Path) -> Result<Option<OpenedPromotionStore>, String> {
    let file = match std::fs::File::open(path) {
        Ok(file) => file,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            return match std::fs::symlink_metadata(path) {
                Err(metadata_error) if metadata_error.kind() == std::io::ErrorKind::NotFound => {
                    Ok(None)
                }
                Ok(_) => Err(format!(
                    "baseline store {} exists but is not an openable regular file",
                    path.display()
                )),
                Err(metadata_error) => Err(format!(
                    "cannot inspect baseline store {}: {metadata_error}",
                    path.display()
                )),
            };
        }
        Err(error) => {
            return Err(format!(
                "cannot open baseline store {}: {error}",
                path.display()
            ));
        }
    };
    let handle_metadata = file
        .metadata()
        .map_err(|error| format!("cannot inspect open store {}: {error}", path.display()))?;
    let identity = promotion_file_identity(path, &handle_metadata)?;
    if identity.links != 1 {
        return Err(format!(
            "baseline store {} must have exactly one hard link, observed {}",
            path.display(),
            identity.links
        ));
    }
    validate_promotion_path_identity(path, identity, "baseline store")?;
    Ok(Some(OpenedPromotionStore {
        file,
        identity,
        permissions: handle_metadata.permissions(),
    }))
}

#[cfg(unix)]
fn promotion_staging_path(store: &Path, nonce: u128, ordinal: u64) -> Result<PathBuf, String> {
    sidecar_path(
        store,
        &format!(".fs-roofline-next-{nonce:032x}-{ordinal:016x}"),
    )
}

#[cfg(unix)]
fn create_promotion_staging_file(store: &Path) -> Result<(PathBuf, std::fs::File), String> {
    static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let time = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |duration| duration.as_nanos());
    let nonce = time ^ (u128::from(std::process::id()) << 64);
    for _ in 0..128 {
        let ordinal = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let path = promotion_staging_path(store, nonce, ordinal)?;
        match std::fs::OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&path)
        {
            Ok(file) => return Ok((path, file)),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(error) => {
                return Err(format!(
                    "cannot create baseline staging file {}: {error}",
                    path.display()
                ));
            }
        }
    }
    Err("cannot allocate a unique create-new baseline staging generation".to_string())
}

/// Serialize promotion writers, re-read under the lock, and replace the store
/// only after the complete bounded next generation is durable. The stable lock
/// file lives outside the source tree. Each same-directory staging generation
/// is opened with create-new semantics, identity-checked through its open
/// handle, and never aliases or truncates a generation left by an earlier
/// crash.
#[cfg(not(unix))]
fn update_promoted_store(
    _store_path: &Path,
    _baseline: fs_roofline::BaselineAxes,
) -> Result<(), String> {
    Err("durable atomic baseline promotion currently requires Unix file identities".to_string())
}

#[cfg(unix)]
#[allow(clippy::too_many_lines)] // One lock/read/stage/replace durability transaction.
fn update_promoted_store(
    store_path: &Path,
    baseline: fs_roofline::BaselineAxes,
) -> Result<(), String> {
    let lock_path = promotion_lock_path(store_path)?;
    let lock = std::fs::OpenOptions::new()
        .create(true)
        .truncate(false)
        .read(true)
        .write(true)
        .open(&lock_path)
        .map_err(|error| {
            format!(
                "cannot open promotion lock {}: {error}",
                lock_path.display()
            )
        })?;
    lock.try_lock().map_err(|error| {
        format!(
            "another promotion is updating baseline store {}: {error}",
            store_path.display()
        )
    })?;

    let existing = open_promotion_store(store_path)?;
    let mut store = if let Some(existing) = &existing {
        let parsed =
            parse_bounded_baseline_store(&existing.file, &store_path.display().to_string())?;
        let handle_metadata = existing.file.metadata().map_err(|error| {
            format!(
                "cannot re-inspect open store {}: {error}",
                store_path.display()
            )
        })?;
        if promotion_file_identity(store_path, &handle_metadata)? != existing.identity {
            return Err(format!(
                "baseline store {} changed while it was read",
                store_path.display()
            ));
        }
        validate_promotion_path_identity(store_path, existing.identity, "baseline store")?;
        parsed
    } else {
        BaselineStore::new()
    };
    store.admit(baseline).map_err(|error| error.to_string())?;
    let rendered = store.to_jsonl();
    if rendered.len() > fs_roofline::baseline::MAX_BASELINE_STORE_BYTES {
        return Err("promoted baseline store exceeded its canonical byte bound".to_string());
    }

    let (next_path, mut next) = create_promotion_staging_file(store_path)?;
    if let Some(existing) = &existing {
        next.set_permissions(existing.permissions.clone())
            .map_err(|error| {
                format!(
                    "cannot preserve baseline-store permissions on {}: {error}",
                    next_path.display()
                )
            })?;
    }
    next.write_all(rendered.as_bytes())
        .and_then(|()| next.sync_all())
        .map_err(|error| {
            format!(
                "cannot durably stage baseline store {}: {error}",
                next_path.display()
            )
        })?;
    let staged_metadata = next.metadata().map_err(|error| {
        format!(
            "cannot inspect staged baseline {}: {error}",
            next_path.display()
        )
    })?;
    let staged_identity = promotion_file_identity(&next_path, &staged_metadata)?;
    if staged_identity.links != 1 {
        return Err(format!(
            "baseline staging file {} unexpectedly has {} hard links",
            next_path.display(),
            staged_identity.links
        ));
    }
    validate_promotion_path_identity(&next_path, staged_identity, "baseline staging file")?;

    match &existing {
        Some(existing) => {
            let handle_metadata = existing.file.metadata().map_err(|error| {
                format!(
                    "cannot re-inspect open store {}: {error}",
                    store_path.display()
                )
            })?;
            if promotion_file_identity(store_path, &handle_metadata)? != existing.identity {
                return Err(format!(
                    "baseline store {} changed before replacement",
                    store_path.display()
                ));
            }
            validate_promotion_path_identity(store_path, existing.identity, "baseline store")?;
        }
        None => match std::fs::symlink_metadata(store_path) {
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Ok(_) => {
                return Err(format!(
                    "baseline store {} appeared during promotion; refusing to overwrite it",
                    store_path.display()
                ));
            }
            Err(error) => {
                return Err(format!(
                    "cannot re-inspect absent baseline store {}: {error}",
                    store_path.display()
                ));
            }
        },
    }
    std::fs::rename(&next_path, store_path).map_err(|error| {
        format!(
            "cannot atomically replace baseline store {} from {}: {error}",
            store_path.display(),
            next_path.display()
        )
    })?;
    validate_promotion_path_identity(store_path, staged_identity, "promoted baseline store")?;
    drop(next);
    let parent = store_path
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    std::fs::File::open(parent)
        .and_then(|directory| directory.sync_all())
        .map_err(|error| {
            format!(
                "baseline store {} was replaced but its directory durability could not be confirmed: {error}",
                store_path.display()
            )
        })?;
    Ok(())
}

fn parse_args(args: &[String]) -> Result<CliArgs, String> {
    let mut parsed = CliArgs::default();
    let mut seen = std::collections::BTreeSet::new();
    let mut args = args.iter().skip(1);
    while let Some(flag) = args.next() {
        if !matches!(
            flag.as_str(),
            "--n"
                | "--warmup"
                | "--reps"
                | "--ledger"
                | "--baseline"
                | "--firmware"
                | "--authority-policy"
                | "--retained-receipts"
        ) {
            return Err(format!("unknown roofline argument {flag:?}"));
        }
        if !seen.insert(flag.as_str()) {
            return Err(format!("duplicate roofline argument {flag:?}"));
        }
        let value = args
            .next()
            .filter(|value| !value.starts_with("--"))
            .ok_or_else(|| format!("{flag} requires a value"))?;
        if value.is_empty() {
            return Err(format!("{flag} requires a non-empty value"));
        }
        match flag.as_str() {
            "--n" => parsed.n = positive_usize(flag, value)?,
            "--warmup" => parsed.warmup = positive_usize(flag, value)?,
            "--reps" => parsed.reps = positive_usize(flag, value)?,
            "--ledger" => parsed.ledger_path = Some(value.clone()),
            "--baseline" => parsed.baseline_path = Some(value.clone()),
            "--firmware" => parsed.firmware = Some(value.clone()),
            "--authority-policy" => parsed.authority_policy_path = Some(value.clone()),
            "--retained-receipts" => parsed.retained_receipts_path = Some(value.clone()),
            _ => return Err(format!("unknown roofline argument {flag:?}")),
        }
    }
    ProductionRunConfig {
        n: parsed.n,
        warmup: parsed.warmup,
        reps: parsed.reps,
    }
    .validate()?;
    Ok(parsed)
}

enum BaselineSource {
    None,
    Plain(BaselineStore),
    Attested(AttestedBaselineStore),
    Invalid(String),
}

impl BaselineSource {
    fn candidate_for_fingerprint(&self, fingerprint: u64) -> Option<&BaselineAxes> {
        match self {
            Self::Plain(store) => store.for_fingerprint(fingerprint),
            Self::Attested(store) => store.for_fingerprint(fingerprint),
            Self::None | Self::Invalid(_) => None,
        }
    }
}

struct BaselineInputs {
    source: BaselineSource,
    identity: BaselineIdentity,
    now_day: u64,
    preliminary_refusal: Option<String>,
    authority: Result<ConfiguredPromotionAuthority, String>,
    retained_receipts: Result<BTreeSet<ContentHash>, String>,
}

enum RunBaselinePolicy<'a> {
    Attested(AttestedAxisBaselinePolicy),
    ReportOnly {
        policy: AxisBaselinePolicy<'a>,
        refusal: String,
    },
}

impl BaselineInputs {
    fn report_only(&self, fingerprint: u64, refusal: String) -> RunBaselinePolicy<'_> {
        RunBaselinePolicy::ReportOnly {
            policy: AxisBaselinePolicy::new(
                self.source.candidate_for_fingerprint(fingerprint),
                &self.identity,
                self.now_day,
            ),
            refusal,
        }
    }

    fn policy(&self, fingerprint: u64) -> RunBaselinePolicy<'_> {
        if let Some(refusal) = &self.preliminary_refusal {
            let source_context = match &self.source {
                BaselineSource::None => Some("no baseline store was supplied"),
                BaselineSource::Plain(_) => {
                    Some("plain baseline stores are candidate/report-only inputs")
                }
                BaselineSource::Invalid(error) => Some(error.as_str()),
                BaselineSource::Attested(_) => None,
            };
            return self.report_only(
                fingerprint,
                source_context.map_or_else(
                    || refusal.clone(),
                    |context| format!("{refusal}; {context}"),
                ),
            );
        }
        let BaselineSource::Attested(store) = &self.source else {
            let refusal = match &self.source {
                BaselineSource::None => "no baseline store was supplied".to_string(),
                BaselineSource::Plain(_) => {
                    "plain baseline stores are candidate/report-only inputs".to_string()
                }
                BaselineSource::Invalid(error) => error.clone(),
                BaselineSource::Attested(_) => unreachable!(),
            };
            return self.report_only(fingerprint, refusal);
        };
        let authority = match &self.authority {
            Ok(authority) => authority,
            Err(error) => return self.report_only(fingerprint, error.clone()),
        };
        let retained_receipts = match &self.retained_receipts {
            Ok(receipts) => receipts,
            Err(error) => return self.report_only(fingerprint, error.clone()),
        };
        match store.policy_for_run(&self.identity, authority, retained_receipts) {
            Ok(policy) => RunBaselinePolicy::Attested(policy),
            Err(error) => self.report_only(
                fingerprint,
                format!("attested baseline authority refused: {error}"),
            ),
        }
    }
}

fn read_bounded_utf8(
    reader: impl Read,
    source: &str,
    kind: &str,
    limit: usize,
) -> Result<String, String> {
    let bounded_bytes = limit
        .checked_add(1)
        .ok_or_else(|| format!("{kind} read bound overflows usize"))?;
    let read_limit =
        u64::try_from(bounded_bytes).map_err(|_| format!("{kind} read bound does not fit u64"))?;
    let mut bytes = Vec::with_capacity(bounded_bytes);
    reader
        .take(read_limit)
        .read_to_end(&mut bytes)
        .map_err(|error| format!("cannot read {kind} {source:?}: {error}"))?;
    if bytes.len() > limit {
        return Err(format!("{kind} {source:?} exceeds the {limit}-byte bound"));
    }
    String::from_utf8(bytes).map_err(|_| format!("{kind} {source:?} is not UTF-8"))
}

fn parse_bounded_baseline_store(reader: impl Read, source: &str) -> Result<BaselineStore, String> {
    let text = read_bounded_utf8(
        reader,
        source,
        "baseline store",
        fs_roofline::baseline::MAX_BASELINE_STORE_BYTES,
    )?;
    BaselineStore::from_jsonl(&text).map_err(|error| error.to_string())
}

fn read_bounded_path(path: &str, kind: &str, limit: usize) -> Result<String, String> {
    let file = std::fs::File::open(path)
        .map_err(|error| format!("cannot read {kind} {path:?}: {error}"))?;
    read_bounded_utf8(file, path, kind, limit)
}

fn parse_baseline_source(text: &str) -> Result<BaselineSource, String> {
    if text.starts_with("{\"record\":") {
        AttestedBaselineStore::from_jsonl(text)
            .map(BaselineSource::Attested)
            .map_err(|error| format!("invalid attested baseline store: {error}"))
    } else {
        BaselineStore::from_jsonl(text)
            .map(BaselineSource::Plain)
            .map_err(|error| format!("invalid plain baseline store: {error}"))
    }
}

fn parse_retained_receipts(
    reader: impl Read,
    source: &str,
) -> Result<BTreeSet<ContentHash>, String> {
    let text = read_bounded_utf8(
        reader,
        source,
        "retained-receipt set",
        MAX_RETAINED_RECEIPTS_INPUT_BYTES,
    )?;
    let body = text.strip_suffix('\n').ok_or_else(|| {
        "retained-receipt set must be canonical newline-terminated lowercase hex".to_string()
    })?;
    if body.is_empty() {
        return Err("retained-receipt set must contain at least one receipt".to_string());
    }
    let mut receipts = BTreeSet::new();
    let mut previous = None;
    for (index, line) in body.split('\n').enumerate() {
        if line.len() != 64
            || !line
                .bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        {
            return Err(format!(
                "retained-receipt line {} must be exactly 64 lowercase hexadecimal bytes",
                index + 1
            ));
        }
        let receipt = ContentHash::from_hex(line)
            .ok_or_else(|| format!("retained-receipt line {} is not a content hash", index + 1))?;
        if previous.is_some_and(|prior| receipt <= prior) {
            return Err(format!(
                "retained-receipt line {} is not in strict ascending order",
                index + 1
            ));
        }
        previous = Some(receipt);
        let inserted = receipts.insert(receipt);
        debug_assert!(inserted);
    }
    Ok(receipts)
}

fn load_baseline_inputs(args: &CliArgs, axes: &MachineAxes) -> Result<BaselineInputs, String> {
    let declared_firmware = args.firmware.as_deref();
    let mut preliminary_refusal = declared_firmware
        .is_none()
        .then(|| "--firmware was not supplied".to_string());
    let identity =
        match BaselineIdentity::current(axes, declared_firmware.unwrap_or("unbaselined-candidate"))
        {
            Ok(identity) => identity,
            Err(error) => {
                preliminary_refusal = Some(format!("invalid baseline identity: {error}"));
                BaselineIdentity::current(axes, "unbaselined-candidate")
                    .map_err(|fallback_error| fallback_error.to_string())?
            }
        };
    let now_day = match days_since_epoch_now() {
        Ok(day) => day,
        Err(error) => {
            let clock_refusal = format!("cannot establish baseline age: {error}");
            preliminary_refusal = Some(match preliminary_refusal {
                Some(prior) => format!("{prior}; {clock_refusal}"),
                None => clock_refusal,
            });
            0
        }
    };
    let source = match args.baseline_path.as_deref() {
        Some(path) => match read_bounded_path(
            path,
            "baseline store",
            fs_roofline::baseline::MAX_BASELINE_STORE_BYTES,
        ) {
            Ok(text) => parse_baseline_source(&text).unwrap_or_else(BaselineSource::Invalid),
            Err(error) => BaselineSource::Invalid(error),
        },
        None => BaselineSource::None,
    };
    let authority = args.authority_policy_path.as_deref().map_or_else(
        || Err("--authority-policy was not supplied".to_string()),
        |path| {
            read_bounded_path(
                path,
                "promotion-authority policy",
                MAX_AUTHORITY_INPUT_BYTES,
            )
            .and_then(|text| {
                ConfiguredPromotionAuthority::from_text(&text)
                    .map_err(|error| format!("invalid promotion-authority policy: {error}"))
            })
        },
    );
    let retained_receipts = args.retained_receipts_path.as_deref().map_or_else(
        || Err("--retained-receipts was not supplied".to_string()),
        |path| {
            let file = std::fs::File::open(path)
                .map_err(|error| format!("cannot read retained-receipt set {path:?}: {error}"))?;
            parse_retained_receipts(file, path)
        },
    );
    Ok(BaselineInputs {
        source,
        identity,
        now_day,
        preliminary_refusal,
        authority,
        retained_receipts,
    })
}

enum CliProductionRun {
    Attested(ProductionRun),
    ReportOnly(ReportOnlyProductionRun),
}

impl CliProductionRun {
    fn axes(&self) -> &MachineAxes {
        match self {
            Self::Attested(run) => run.axes(),
            Self::ReportOnly(run) => run.axes(),
        }
    }

    fn post_axes(&self) -> &MachineAxes {
        match self {
            Self::Attested(run) => run.post_axes(),
            Self::ReportOnly(run) => run.post_axes(),
        }
    }

    fn results(&self) -> &[fs_roofline::Attainment] {
        match self {
            Self::Attested(run) => run.results(),
            Self::ReportOnly(run) => run.results(),
        }
    }

    fn receipt_json(&self) -> &str {
        match self {
            Self::Attested(run) => run.receipt_json(),
            Self::ReportOnly(run) => run.receipt_json(),
        }
    }

    fn baseline_hash(&self) -> Option<ContentHash> {
        match self {
            Self::Attested(run) => run.baseline_hash(),
            Self::ReportOnly(run) => run.baseline_hash(),
        }
    }

    fn evidence_admission(&self) -> (bool, Option<String>) {
        match self {
            Self::Attested(run) => {
                let refusal = run.admission_error();
                (refusal.is_none(), refusal)
            }
            Self::ReportOnly(run) => (false, run.admission_error()),
        }
    }

    fn protocol(&self) -> &'static str {
        match self {
            Self::Attested(_) => PRODUCTION_PROTOCOL_VERSION,
            Self::ReportOnly(_) => CUSTOM_REGISTRY_PROTOCOL_VERSION,
        }
    }

    fn record(self, ledger: &fs_ledger::Ledger) -> Result<i64, fs_ledger::LedgerError> {
        match self {
            Self::Attested(run) => run.record(ledger),
            Self::ReportOnly(run) => run.record(ledger),
        }
    }
}

fn main() -> std::process::ExitCode {
    let raw_args: Vec<String> = std::env::args().collect();
    if raw_args.get(1).is_some_and(|arg| arg == "promote") {
        return match parse_promote_args(&raw_args).and_then(|args| run_promote(&args)) {
            Ok(()) => std::process::ExitCode::SUCCESS,
            Err(error) => fail(&error),
        };
    }
    let args = match parse_args(&raw_args) {
        Ok(args) => args,
        Err(error) => return fail(&error),
    };

    let tune_ledger = match args.ledger_path.as_deref() {
        Some(path) => match fs_ledger::Ledger::open(path) {
            Ok(ledger) => Some(ledger),
            Err(error) => return fail(&error.to_string()),
        },
        None => None,
    };

    // Sealed production protocol (bead fz2.5): the CLI never supplies axes,
    // kernels, or the post-probe — it observes the probe (baseline selection
    // needs the fingerprint), then hands the whole lifecycle to the protocol.
    let probe = ProductionProbe::observe();
    println!("{}", probe.axes().to_jsonl());

    let baseline_inputs = match load_baseline_inputs(&args, probe.axes()) {
        Ok(inputs) => inputs,
        Err(error) => return fail(&error),
    };
    let baseline_policy = baseline_inputs.policy(probe.axes().fingerprint);

    let config = ProductionRunConfig {
        n: args.n,
        warmup: args.warmup,
        reps: args.reps,
    };
    let run = match baseline_policy {
        RunBaselinePolicy::Attested(policy) => match probe.run(config, policy, tune_ledger) {
            Ok(run) => CliProductionRun::Attested(run),
            Err(error) => return fail(&error),
        },
        RunBaselinePolicy::ReportOnly { policy, refusal } => {
            match probe
                .run_report_only(config, policy, tune_ledger)
                .and_then(|run| run.with_configuration_refusal(refusal))
            {
                Ok(run) => CliProductionRun::ReportOnly(run),
                Err(error) => return fail(&error),
            }
        }
    };
    println!("{}", run.post_axes().to_jsonl());
    println!("{}", run.receipt_json());
    let (citation_eligible, admission_error) = run.evidence_admission();
    println!(
        "{}",
        evidence_admission_json(citation_eligible, admission_error.as_deref())
    );
    for r in run.results() {
        println!("{}", r.to_jsonl());
    }
    for row in SECTION_14_1_TARGETS {
        println!(
            "{{\"target\":\"{}\",\"statement\":\"{}\",\"landed\":{}}}",
            json_escape(row.kernel),
            json_escape(row.statement),
            row.landed
        );
    }

    if let Some(db) = args.ledger_path {
        let ledger = match fs_ledger::Ledger::open(&db) {
            Ok(l) => l,
            Err(e) => return fail(&e.to_string()),
        };
        let fingerprint = run.axes().fingerprint;
        let kernel_ids: Vec<(String, String)> = run
            .results()
            .iter()
            .map(|r| (r.kernel.clone(), r.version.clone()))
            .collect();
        let baseline_hash = run.baseline_hash();
        let protocol = run.protocol();
        let op = match run.record(&ledger) {
            Ok(op) => op,
            Err(e) => return fail(&e.to_string()),
        };
        let mut citable = citation_eligible;
        for (kernel, version) in &kernel_ids {
            match staleness(&ledger, kernel, version, fingerprint, baseline_hash) {
                Ok(s) => {
                    citable &= s == fs_roofline::Staleness::Fresh;
                    println!(
                        "{{\"kernel\":\"{}\",\"staleness\":\"{s:?}\",\"max_age_ns\":{STALENESS_MAX_AGE_NS}}}",
                        json_escape(kernel),
                    );
                }
                Err(e) => return fail(&e.to_string()),
            }
        }
        let reason = if citable {
            "null"
        } else if citation_eligible {
            "\"post_commit_validation_failed\""
        } else {
            "\"admission_refused\""
        };
        println!(
            "{{\"schema\":\"fs-roofline-recorded-evidence-v1\",\"ledgered\":true,\"citation_eligible\":{citation_eligible},\"citable\":{citable},\"reason\":{reason},\"protocol\":\"{protocol}\",\"op\":{op},\"db\":\"{}\"}}",
            json_escape(&db)
        );
    }
    std::process::ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::{
        BaselineInputs, BaselineSource, MAX_BASELINE_AGE_DAYS, MAX_PROMOTION_PROBES,
        RunBaselinePolicy, evidence_admission_json, json_escape, load_baseline_inputs, parse_args,
        parse_baseline_source, parse_bounded_baseline_store, parse_promote_args,
        parse_retained_receipts, promotion_lock_path, sidecar_path,
    };
    #[cfg(unix)]
    use super::{open_promotion_store, promotion_staging_path};
    use fs_roofline::authority::{ConfiguredPromotionAuthority, PromotionAttestation};
    use fs_roofline::{
        AttestedBaselineStore, AxisAdmissionSnapshot, BaselineCandidate, BaselineIdentity,
        MachineAxes, days_since_epoch_now, promote_baseline,
    };
    use std::io::Cursor;

    fn args(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_string()).collect()
    }

    fn stable_attested_snapshot(
        store_text: &str,
        identity: &BaselineIdentity,
        authority_text: &str,
        retained: &std::collections::BTreeSet<fs_blake3::ContentHash>,
        axes: &MachineAxes,
    ) -> AxisAdmissionSnapshot {
        for _attempt in 0..2 {
            let day_before = days_since_epoch_now().expect("test clock follows the Unix epoch");
            let inputs = BaselineInputs {
                source: BaselineSource::Attested(
                    AttestedBaselineStore::from_jsonl(store_text)
                        .expect("attested store round trip"),
                ),
                identity: identity.clone(),
                now_day: day_before,
                preliminary_refusal: None,
                authority: Ok(ConfiguredPromotionAuthority::from_text(authority_text)
                    .expect("canonical configured authority")),
                retained_receipts: Ok(retained.clone()),
            };
            let RunBaselinePolicy::Attested(policy) = inputs.policy(axes.fingerprint) else {
                panic!("the fully configured CLI path must mint an opaque attested policy");
            };
            let snapshot = policy.decide(axes, axes);
            let day_after = days_since_epoch_now().expect("test clock follows the Unix epoch");
            if day_before == day_after {
                return snapshot;
            }
        }
        panic!("test clock crossed UTC midnight twice while minting an attested snapshot");
    }

    #[test]
    fn citation_eligibility_is_never_reported_as_precommit_citation() {
        let eligible = evidence_admission_json(true, None);
        assert!(eligible.contains("\"citation_eligible\":true"));
        assert!(eligible.contains("\"recorded\":false"));
        assert!(eligible.contains("\"citable\":false"));
        assert!(!eligible.contains("\"citable\":true"));

        let refused = evidence_admission_json(false, Some("development salt"));
        assert!(refused.contains("\"citation_eligible\":false"));
        assert!(refused.contains("development salt"));
        assert!(refused.contains("\"citable\":false"));
    }

    #[test]
    fn manual_json_fields_escape_hostile_paths_and_diagnostics() {
        assert_eq!(
            json_escape("ledger\\\"row\n\t\u{0001}.db"),
            "ledger\\\\\\\"row\\n\\t\\u0001.db"
        );
    }

    #[test]
    fn missing_firmware_or_baseline_remains_report_only() {
        let axes = MachineAxes {
            fingerprint: 1,
            cpu_brand: "synthetic".to_string(),
            logical_cpus: 1,
            bandwidth_single_gbs: 10.0,
            bandwidth_all_core_gbs: 10.0,
            peak_single_gflops: 10.0,
            peak_all_core_gflops: 10.0,
        };
        let incomplete = parse_args(&args(&["roofline", "--baseline", "x"]))
            .expect("partial authority configuration remains measurable");
        let inputs = load_baseline_inputs(&incomplete, &axes)
            .expect("partial authority configuration is report-only");
        assert!(matches!(
            inputs.policy(axes.fingerprint),
            RunBaselinePolicy::ReportOnly { .. }
        ));

        let parsed = parse_args(&args(&["roofline"])).expect("default invocation parses");
        let candidate =
            load_baseline_inputs(&parsed, &axes).expect("report-only invocation remains available");
        assert!(matches!(
            candidate.policy(axes.fingerprint),
            RunBaselinePolicy::ReportOnly { .. }
        ));
    }

    #[test]
    fn configured_cli_policy_mints_exact_snapshot_and_refusals_stay_report_only() {
        let axes = MachineAxes {
            fingerprint: 0xC11,
            cpu_brand: "synthetic-cli-authority".to_string(),
            logical_cpus: 4,
            bandwidth_single_gbs: 100.0,
            bandwidth_all_core_gbs: 200.0,
            peak_single_gflops: 300.0,
            peak_all_core_gflops: 600.0,
        };
        let identity =
            BaselineIdentity::current(&axes, "cli-firmware-a").expect("synthetic CLI identity");
        let now_day = days_since_epoch_now().expect("test clock follows the Unix epoch");
        let mut retained = std::collections::BTreeSet::new();
        let candidates = (0_u64..3)
            .map(|ordinal| {
                let receipt = fs_blake3::hash_domain(
                    "fs-roofline.cli-configured-test-source.v1",
                    &ordinal.to_le_bytes(),
                );
                assert!(retained.insert(receipt));
                BaselineCandidate::from_receipt(axes.clone(), identity.clone(), receipt)
                    .expect("synthetic retained candidate")
            })
            .collect::<Vec<_>>();
        let baseline = promote_baseline(
            &candidates,
            "cli-test-operator",
            "configured CLI authority fixture",
            now_day,
            90,
        )
        .expect("synthetic CLI baseline");
        let authority_text = format!(
            "authorize\tops/cli-test\t{}\tcli-test-signature\n",
            baseline.content_hash()
        );
        let authority = ConfiguredPromotionAuthority::from_text(&authority_text)
            .expect("canonical configured authority");
        let policy_receipt = authority.policy_receipt();
        let mut store = AttestedBaselineStore::new();
        store
            .admit_verified(
                baseline.clone(),
                PromotionAttestation::new("ops/cli-test", "cli-test-signature"),
                &authority,
                &retained,
            )
            .expect("configured authority admits the exact baseline");
        let store_text = store.to_jsonl();

        let snapshot =
            stable_attested_snapshot(&store_text, &identity, &authority_text, &retained, &axes);
        assert!(snapshot.authority_admitted());
        assert!(snapshot.verdict().trusted());
        assert!(snapshot.baseline_citation_eligible());
        assert!(snapshot.receipt_json().contains("\"tier\":\"attested\""));
        assert!(
            snapshot
                .receipt_json()
                .contains(&format!("\"policy_receipt\":\"{policy_receipt}\""))
        );

        let rotated_authority_text = format!(
            "authorize\tops/cli-rotated\t{}\tcli-rotated-signature\nrevoke\tops/cli-test\n",
            baseline.content_hash()
        );
        let rotated_authority = ConfiguredPromotionAuthority::from_text(&rotated_authority_text)
            .expect("canonical rotated authority");
        let rotated_policy_receipt = rotated_authority.policy_receipt();
        let mut rotated_store = AttestedBaselineStore::new();
        rotated_store
            .admit_verified(
                baseline.clone(),
                PromotionAttestation::new("ops/cli-rotated", "cli-rotated-signature"),
                &rotated_authority,
                &retained,
            )
            .expect("the same immutable baseline can be re-attested under the rotated key");
        let rotated_store_text = rotated_store.to_jsonl();
        let rotated_snapshot = stable_attested_snapshot(
            &rotated_store_text,
            &identity,
            &rotated_authority_text,
            &retained,
            &axes,
        );
        assert_eq!(
            rotated_snapshot.baseline_hash(),
            Some(baseline.content_hash()),
            "rotation re-endorses the same immutable baseline"
        );
        assert!(rotated_snapshot.baseline_citation_eligible());
        assert!(
            rotated_snapshot
                .receipt_json()
                .contains("\"key_id\":\"ops/cli-rotated\"")
        );
        assert!(
            rotated_snapshot
                .receipt_json()
                .contains(&format!("\"policy_receipt\":\"{rotated_policy_receipt}\""))
        );

        let assert_report_only = |candidate_store: &str,
                                  candidate_identity: BaselineIdentity,
                                  policy_text: &str,
                                  available: std::collections::BTreeSet<fs_blake3::ContentHash>,
                                  expected: &str| {
            let inputs = BaselineInputs {
                source: BaselineSource::Attested(
                    AttestedBaselineStore::from_jsonl(candidate_store)
                        .expect("structural attested store"),
                ),
                identity: candidate_identity,
                now_day,
                preliminary_refusal: None,
                authority: ConfiguredPromotionAuthority::from_text(policy_text)
                    .map_err(|error| error.to_string()),
                retained_receipts: Ok(available),
            };
            match inputs.policy(axes.fingerprint) {
                RunBaselinePolicy::ReportOnly { refusal, .. } => assert!(
                    refusal.contains(expected),
                    "unexpected CLI authority refusal: {refusal}"
                ),
                RunBaselinePolicy::Attested(_) => {
                    panic!("{expected} fixture must remain report-only")
                }
            }
        };
        assert_report_only(
            &store_text,
            identity.clone(),
            "",
            retained.clone(),
            "unknown-key",
        );
        assert_report_only(
            &store_text,
            identity.clone(),
            "revoke\tops/cli-test\n",
            retained.clone(),
            "revoked-key",
        );
        assert_report_only(
            &store_text,
            identity.clone(),
            &rotated_authority_text,
            retained.clone(),
            "revoked-key",
        );
        let tampered = store_text.replace("cli-test-signature", "tampered-signature");
        assert_report_only(
            &tampered,
            identity.clone(),
            &authority_text,
            retained.clone(),
            "wrong-signature",
        );
        let mut missing = retained.clone();
        let dropped = *baseline
            .provenance()
            .source_receipts()
            .first()
            .expect("promoted source receipt");
        assert!(missing.remove(&dropped));
        assert_report_only(
            &store_text,
            identity.clone(),
            &authority_text,
            missing,
            &dropped.to_string(),
        );
        let cross_machine =
            BaselineIdentity::current(&axes, "cli-firmware-b").expect("cross-machine CLI identity");
        assert_report_only(
            &store_text,
            cross_machine,
            &authority_text,
            retained,
            "does not match",
        );
    }

    #[test]
    fn parser_rejects_unknown_duplicate_missing_and_invalid_values() {
        for (case, expected) in [
            (vec!["roofline", "--unknown", "x"], "unknown"),
            (vec!["roofline", "--n", "1", "--n", "2"], "duplicate"),
            (vec!["roofline", "--ledger"], "requires a value"),
            (vec!["roofline", "--ledger", "--n", "1"], "requires a value"),
            (vec!["roofline", "--n", "0"], "positive integer"),
            (vec!["roofline", "--reps", "nope"], "positive integer"),
        ] {
            let error = parse_args(&args(&case)).expect_err("malformed argv must fail");
            assert!(
                error.contains(expected),
                "{error:?} did not contain {expected:?}"
            );
        }
    }

    #[test]
    fn parser_accepts_every_flag_once_and_preserves_report_only_default() {
        let defaults = parse_args(&args(&["roofline"])).expect("defaults");
        assert!(defaults.baseline_path.is_none());
        assert!(defaults.ledger_path.is_none());
        assert!(defaults.authority_policy_path.is_none());
        assert!(defaults.retained_receipts_path.is_none());

        let parsed = parse_args(&args(&[
            "roofline",
            "--n",
            "8",
            "--warmup",
            "1",
            "--reps",
            "2",
            "--ledger",
            "run.db",
            "--baseline",
            "axes.jsonl",
            "--firmware",
            "os-build-1",
            "--authority-policy",
            "authority.tsv",
            "--retained-receipts",
            "receipts.txt",
        ]))
        .expect("complete argv");
        assert_eq!(parsed.n, 8);
        assert_eq!(parsed.warmup, 1);
        assert_eq!(parsed.reps, 2);
        assert_eq!(parsed.ledger_path.as_deref(), Some("run.db"));
        assert_eq!(parsed.baseline_path.as_deref(), Some("axes.jsonl"));
        assert_eq!(parsed.firmware.as_deref(), Some("os-build-1"));
        assert_eq!(
            parsed.authority_policy_path.as_deref(),
            Some("authority.tsv")
        );
        assert_eq!(
            parsed.retained_receipts_path.as_deref(),
            Some("receipts.txt")
        );
    }

    #[test]
    fn parser_refuses_resource_inputs_above_the_production_envelope() {
        let too_many_elements = fs_roofline::production::MAX_PRODUCTION_ELEMENTS.saturating_add(1);
        let too_many_warmups = fs_roofline::production::MAX_PRODUCTION_WARMUP.saturating_add(1);
        let too_many_reps = fs_roofline::production::MAX_PRODUCTION_REPS.saturating_add(1);
        for (flag, value, expected) in [
            ("--n", too_many_elements, "production n"),
            ("--warmup", too_many_warmups, "production warmup"),
            ("--reps", too_many_reps, "production reps"),
        ] {
            let error = parse_args(&args(&["roofline", flag, &value.to_string()]))
                .expect_err("out-of-envelope resource input must fail before probing");
            assert!(error.contains(expected), "unexpected diagnostic: {error}");
        }
        let max_n = fs_roofline::production::MAX_PRODUCTION_ELEMENTS.to_string();
        let max_warmup = fs_roofline::production::MAX_PRODUCTION_WARMUP.to_string();
        let error = parse_args(&args(&["roofline", "--n", &max_n, "--warmup", &max_warmup]))
            .expect_err("an oversized combined loop budget must fail before probing");
        assert!(error.contains("warmup + reps"));
    }

    #[test]
    fn promotion_parser_bounds_probe_allocation_and_age() {
        let common = [
            "roofline",
            "promote",
            "--store",
            "store.jsonl",
            "--firmware",
            "firmware",
            "--operator",
            "operator",
            "--justification",
            "calibration",
        ];
        let mut probes = common.to_vec();
        let probes_limit = MAX_PROMOTION_PROBES.saturating_add(1).to_string();
        probes.extend(["--probes", probes_limit.as_str()]);
        assert!(
            parse_promote_args(&args(&probes))
                .err()
                .expect("oversized probe count must fail")
                .contains("at most")
        );

        let mut age = common.to_vec();
        let age_limit = MAX_BASELINE_AGE_DAYS.saturating_add(1).to_string();
        age.extend(["--age-days", age_limit.as_str()]);
        assert!(
            parse_promote_args(&args(&age))
                .err()
                .expect("oversized age must fail")
                .contains("at most")
        );
    }

    #[test]
    fn baseline_reader_stops_at_the_store_bound_plus_one_byte() {
        let oversized = vec![b'x'; fs_roofline::baseline::MAX_BASELINE_STORE_BYTES + 1];
        let error = parse_bounded_baseline_store(Cursor::new(oversized), "oversized.jsonl")
            .expect_err("oversized input must fail before parsing");
        assert!(error.contains("exceeds"));
    }

    #[test]
    fn retained_receipts_require_canonical_unique_lowercase_lines() {
        let first = "00".repeat(32);
        let second = "ab".repeat(32);
        let canonical = format!("{first}\n{second}\n");
        let parsed = parse_retained_receipts(Cursor::new(canonical), "receipts.txt")
            .expect("canonical retained receipts");
        assert_eq!(parsed.len(), 2);

        for malformed in [
            first.clone(),
            format!("{}\n", "AB".repeat(32)),
            format!("{first}\n{first}\n"),
            format!("{first}\n\n"),
            format!("{second}\n{first}\n"),
        ] {
            assert!(
                parse_retained_receipts(Cursor::new(malformed), "receipts.txt").is_err(),
                "malformed retained-receipt set must be refused"
            );
        }
    }

    #[test]
    fn retained_receipt_reader_stops_at_the_input_bound_plus_one_byte() {
        let oversized = vec![b'0'; super::MAX_RETAINED_RECEIPTS_INPUT_BYTES + 1];
        let error = parse_retained_receipts(Cursor::new(oversized), "oversized.txt")
            .expect_err("oversized retained-receipt input must fail before parsing");
        assert!(error.contains("exceeds"));
    }

    #[test]
    fn baseline_envelope_prefix_selects_attested_parser_without_a_new_flag() {
        let attested = parse_baseline_source("{\"record\":malformed}\n")
            .err()
            .expect("malformed attested envelope must be refused");
        assert!(attested.contains("attested"), "{attested}");

        let plain = parse_baseline_source("malformed\n")
            .err()
            .expect("malformed plain store must be refused");
        assert!(plain.contains("plain"), "{plain}");
    }

    #[test]
    fn promotion_sidecars_do_not_alias_the_store_and_lock_identity_is_stable() {
        let store = std::env::temp_dir().join("fs-roofline-sidecar-fixture.jsonl");
        let lock_a = promotion_lock_path(&store).expect("lock path");
        let lock_b = promotion_lock_path(&store).expect("repeat lock path");
        assert_ne!(lock_a, store);
        assert_eq!(lock_a, lock_b);
        assert!(lock_a.starts_with(std::env::temp_dir()));

        let ordinary_sidecar = sidecar_path(&store, ".fixture").expect("sidecar path");
        assert_ne!(ordinary_sidecar, store);
        assert_eq!(ordinary_sidecar.parent(), store.parent());
    }

    #[cfg(unix)]
    #[test]
    fn promotion_staging_generations_are_unique_and_same_directory() {
        let store = std::env::temp_dir().join("fs-roofline-staging-fixture.jsonl");
        let first = promotion_staging_path(&store, 7, 11).expect("first staging path");
        let next_nonce = promotion_staging_path(&store, 8, 11).expect("next nonce path");
        let next_ordinal = promotion_staging_path(&store, 7, 12).expect("next ordinal path");

        assert_ne!(first, store);
        assert_ne!(first, next_nonce);
        assert_ne!(first, next_ordinal);
        assert_eq!(first.parent(), store.parent());
        assert_eq!(next_nonce.parent(), store.parent());
        assert_eq!(next_ordinal.parent(), store.parent());
    }

    #[cfg(unix)]
    #[test]
    fn promotion_store_refuses_real_symlink_and_hardlink_paths() {
        use std::os::unix::fs::symlink;

        static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("test clock follows Unix epoch")
            .as_nanos();
        let ordinal = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let prefix = format!(
            "fs-roofline-store-identity-{}-{nonce}-{ordinal}",
            std::process::id()
        );
        let original = std::env::temp_dir().join(format!("{prefix}-original"));
        let symlink_path = std::env::temp_dir().join(format!("{prefix}-symlink"));
        let hardlink_path = std::env::temp_dir().join(format!("{prefix}-hardlink"));
        std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&original)
            .expect("create unique regular-file fixture");

        symlink(&original, &symlink_path).expect("create symlink fixture");
        let Err(symlink_error) = open_promotion_store(&symlink_path) else {
            panic!("symlink store must be refused");
        };
        assert!(symlink_error.contains("regular file"), "{symlink_error}");

        std::fs::hard_link(&original, &hardlink_path).expect("create hardlink fixture");
        let Err(hardlink_error) = open_promotion_store(&hardlink_path) else {
            panic!("hardlinked store must be refused");
        };
        assert!(hardlink_error.contains("hard link"), "{hardlink_error}");
    }
}
