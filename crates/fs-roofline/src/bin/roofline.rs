//! Roofline harness CLI (plan §14.4 nightly lane).
//!
//! Usage:
//!   roofline [--n <elements>] [--warmup <k>] [--reps <k>] [--ledger <db>]
//!            [--baseline <jsonl>] [--firmware <identity>]
//!   roofline promote --store <jsonl> --firmware <identity>
//!            --operator <name> --justification <text>
//!            [--probes <k≥3>] [--age-days <d>]
//!
//! Probes the machine axes, runs the default kernel registry, prints one
//! JSON line per kernel (plus the axes line and the §14.1 coverage table),
//! and — when `--ledger` is given — records the run as ledger provenance
//! and reports staleness for every registered kernel.

use fs_roofline::production::{ProductionProbe, ProductionRunConfig};
use fs_roofline::{
    AxisBaselinePolicy, BaselineIdentity, BaselineStore, MachineAxes, PRODUCTION_PROTOCOL_VERSION,
    SECTION_14_1_TARGETS, STALENESS_MAX_AGE_NS, days_since_epoch_now, staleness,
};
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

/// `roofline promote` — the operator bootstrap for governed baselines
/// (bead c40j): probe the machine axes N ≥ 3 times, build candidates,
/// run [`fs_roofline::promote_baseline`] (which REFUSES on a loaded
/// host — the drift bands are the point), and create-or-update the
/// JSONL store. Until fz2.7 lands signatures, the store is
/// operator-trusted/tamper-evident, not independently verified.
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
        // A content-derived source receipt: the probe's own canonical
        // bytes under a CLI-specific domain (structural traceability;
        // authentication is fz2.7's layer, stated in the store README).
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
            "--n" | "--warmup" | "--reps" | "--ledger" | "--baseline" | "--firmware"
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
            _ => return Err(format!("unknown roofline argument {flag:?}")),
        }
    }
    if parsed.baseline_path.is_some() && parsed.firmware.is_none() {
        return Err("--firmware is required when --baseline is supplied".to_string());
    }
    ProductionRunConfig {
        n: parsed.n,
        warmup: parsed.warmup,
        reps: parsed.reps,
    }
    .validate()?;
    Ok(parsed)
}

struct BaselineInputs {
    store: Option<BaselineStore>,
    identity: BaselineIdentity,
    now_day: u64,
}

impl BaselineInputs {
    fn policy(&self, fingerprint: u64) -> AxisBaselinePolicy<'_> {
        AxisBaselinePolicy::new(
            self.store
                .as_ref()
                .and_then(|store| store.for_fingerprint(fingerprint)),
            &self.identity,
            self.now_day,
        )
    }
}

fn parse_bounded_baseline_store(reader: impl Read, source: &str) -> Result<BaselineStore, String> {
    let limit = fs_roofline::baseline::MAX_BASELINE_STORE_BYTES;
    let bounded_bytes = limit
        .checked_add(1)
        .ok_or_else(|| "baseline-store read bound overflows usize".to_string())?;
    let read_limit = u64::try_from(bounded_bytes)
        .map_err(|_| "baseline-store read bound does not fit u64".to_string())?;
    let mut bytes = Vec::with_capacity(bounded_bytes);
    reader
        .take(read_limit)
        .read_to_end(&mut bytes)
        .map_err(|error| format!("cannot read baseline store {source:?}: {error}"))?;
    if bytes.len() > limit {
        return Err(format!(
            "baseline store {source:?} exceeds the {limit}-byte bound"
        ));
    }
    let text =
        String::from_utf8(bytes).map_err(|_| format!("baseline store {source:?} is not UTF-8"))?;
    BaselineStore::from_jsonl(&text).map_err(|error| error.to_string())
}

fn load_baseline_inputs(args: &CliArgs, axes: &MachineAxes) -> Result<BaselineInputs, String> {
    let identity = BaselineIdentity::current(
        axes,
        args.firmware.as_deref().unwrap_or("unbaselined-candidate"),
    )
    .map_err(|error| error.to_string())?;
    let now_day = days_since_epoch_now().map_err(|error| error.to_string())?;
    let store = match args.baseline_path.as_deref() {
        Some(path) => {
            let file = std::fs::File::open(path)
                .map_err(|error| format!("cannot read baseline store {path:?}: {error}"))?;
            Some(parse_bounded_baseline_store(file, path)?)
        }
        None => None,
    };
    Ok(BaselineInputs {
        store,
        identity,
        now_day,
    })
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
    let run = match probe.run(config, baseline_policy, tune_ledger) {
        Ok(run) => run,
        Err(error) => return fail(&error),
    };
    println!("{}", run.post_axes().to_jsonl());
    println!("{}", run.receipt_json());
    let citation_eligible = run.citation_eligible();
    let admission_error = run.admission_error();
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
        let op = match run.record(&ledger) {
            Ok(op) => op,
            Err(e) => return fail(&e.to_string()),
        };
        let mut citable = citation_eligible;
        for (kernel, version) in &kernel_ids {
            match staleness(
                &ledger,
                kernel,
                version,
                fingerprint,
                baseline_policy.baseline_hash(),
            ) {
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
            "{{\"schema\":\"fs-roofline-recorded-evidence-v1\",\"ledgered\":true,\"citation_eligible\":{citation_eligible},\"citable\":{citable},\"reason\":{reason},\"protocol\":\"{PRODUCTION_PROTOCOL_VERSION}\",\"op\":{op},\"db\":\"{}\"}}",
            json_escape(&db)
        );
    }
    std::process::ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::{
        MAX_BASELINE_AGE_DAYS, MAX_PROMOTION_PROBES, evidence_admission_json, json_escape,
        load_baseline_inputs, parse_args, parse_bounded_baseline_store, parse_promote_args,
        promotion_lock_path, sidecar_path,
    };
    #[cfg(unix)]
    use super::{open_promotion_store, promotion_staging_path};
    use fs_roofline::MachineAxes;
    use std::io::Cursor;

    fn args(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_string()).collect()
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
    fn baseline_store_requires_explicit_firmware_identity() {
        let axes = MachineAxes {
            fingerprint: 1,
            cpu_brand: "synthetic".to_string(),
            logical_cpus: 1,
            bandwidth_single_gbs: 10.0,
            bandwidth_all_core_gbs: 10.0,
            peak_single_gflops: 10.0,
            peak_all_core_gflops: 10.0,
        };
        let error = parse_args(&args(&["roofline", "--baseline", "x"]))
            .expect_err("firmware omission must fail before file access");
        assert!(error.contains("--firmware"));

        let parsed = parse_args(&args(&["roofline"])).expect("default invocation parses");
        let candidate =
            load_baseline_inputs(&parsed, &axes).expect("report-only invocation remains available");
        assert!(
            !candidate
                .policy(axes.fingerprint)
                .verdict(&axes, &axes)
                .trusted()
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
        ]))
        .expect("complete argv");
        assert_eq!(parsed.n, 8);
        assert_eq!(parsed.warmup, 1);
        assert_eq!(parsed.reps, 2);
        assert_eq!(parsed.ledger_path.as_deref(), Some("run.db"));
        assert_eq!(parsed.baseline_path.as_deref(), Some("axes.jsonl"));
        assert_eq!(parsed.firmware.as_deref(), Some("os-build-1"));
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
