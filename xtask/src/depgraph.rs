//! Canonical, operator-observed dependency receipts for GEMM build identity.
//!
//! Stable Cargo does not expose the invocation-exact unit graph to a build
//! script. This command therefore supports one deliberately narrow evidence
//! surface: one explicit production root package, its normal and build edges,
//! an optional target triple, and root feature selection. Test/dev/target-kind
//! and profile selection are refused instead of being silently approximated.
//! The filesystem offers no transaction spanning Cargo graph discovery and
//! package hashing, so the command derives the entire canonical receipt twice
//! and emits nothing unless both byte strings agree exactly.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitStatus, Stdio};
use std::thread;

#[cfg(unix)]
use std::os::unix::fs::MetadataExt as _;

#[path = "../../crates/fs-la/depgraph_receipt_format.rs"]
#[allow(dead_code)] // xtask consumes only the receipt grammar; fs-la owns identity helpers.
mod receipt_format;

const CLOSURE_ROOT: &str = "fs-la";
const MAX_SELECTION_ARGS: usize = 64;
const MAX_ARG_BYTES: usize = 1_024;
// `--no-dedupe` is intentional so the fs-la subtree cannot be truncated by an
// earlier sibling visit. The current orchestration roots expand to ~20 MiB /
// ~580k rows, so bounds must cover that real graph while remaining finite.
const MAX_TREE_BYTES: usize = 32 * 1_024 * 1_024;
const MAX_TREE_LINES: usize = 1_000_000;
const MAX_METADATA_BYTES: usize = 8 * 1_024 * 1_024;
const MAX_SUBPROCESS_STDERR_BYTES: usize = 1_024 * 1_024;
const MAX_FEATURES_PER_PACKAGE: usize = 1_024;
const MAX_FEATURE_BYTES: usize = 256;
const MAX_PATH_PACKAGE_FILES: usize = 100_000;
const MAX_PATH_PACKAGE_DIRECTORIES: usize = 100_000;
const MAX_PATH_PACKAGE_DEPTH: usize = 256;
const MAX_PATH_PACKAGE_BYTES: u64 = 1_073_741_824;
const MAX_PATH_MANIFEST_BYTES: usize = 64 * 1_024 * 1_024;
// Per-package limits alone permit an 8,192-row receipt to trigger terabytes of
// reads. Preflight the complete distinct local-package set, then meter both
// coherence snapshots through one closure-wide work budget.
const MAX_PATH_CLOSURE_PACKAGES: usize = 512;
const MAX_PATH_CLOSURE_READ_BYTES: u64 = 4_294_967_296;
const MAX_PATH_CLOSURE_ENTRIES: usize = 2_000_000;
const MAX_CARGO_EXECUTABLE_BYTES: u64 = 268_435_456;
const PATH_PACKAGE_DOMAIN: &str = "org.frankensim.depgraph.path-package.v2";

#[derive(Debug, Clone, PartialEq, Eq)]
struct Selection {
    package: String,
    target: Option<String>,
    features: BTreeSet<String>,
    all_features: bool,
    default_features: bool,
}

impl Selection {
    fn cargo_args(&self) -> Vec<String> {
        let mut args = vec!["--package".to_string(), self.package.clone()];
        if let Some(target) = &self.target {
            args.extend(["--target".to_string(), target.clone()]);
        }
        if !self.default_features {
            args.push("--no-default-features".to_string());
        }
        if self.all_features {
            args.push("--all-features".to_string());
        } else if !self.features.is_empty() {
            args.push("--features".to_string());
            args.push(self.features.iter().cloned().collect::<Vec<_>>().join(","));
        }
        args
    }
}

/// Run `depgraph-receipt`.
///
/// Command options precede `--`; the typed Cargo selection follows it. The
/// supported shape is `--package <ROOT>` plus optional `--target <TRIPLE>`,
/// `--features <CSV>`, `--all-features`, and `--no-default-features`.
pub fn cmd_depgraph_receipt(root: &Path, raw_args: &[String]) -> Result<(), String> {
    let (verify, selection_args) = parse_command_args(raw_args)?;
    let selection = Selection::parse(selection_args)?;
    let supplied = if verify {
        let supplied = std::env::var("FRANKENSIM_DEPGRAPH_RECEIPT").map_err(|_| {
            "verify mode requires FRANKENSIM_DEPGRAPH_RECEIPT in the environment".to_string()
        })?;
        receipt_format::parse(&supplied)
            .map_err(|error| format!("supplied receipt is malformed/non-canonical: {error}"))?;
        Some(supplied)
    } else {
        None
    };
    let receipt = derive_receipt(root, &selection)?;
    let confirmation = derive_receipt(root, &selection)?;
    require_matching_complete_derivations(&receipt, &confirmation)?;
    if let Some(supplied) = supplied {
        if supplied == receipt {
            println!("depgraph receipt verified: {} bytes", receipt.len());
            Ok(())
        } else {
            Err(format!(
                "depgraph receipt mismatch: environment carries {} bytes, recomputation yields {} bytes; \
                 the selected dependency graph changed or the selection differs",
                supplied.len(),
                receipt.len()
            ))
        }
    } else {
        println!("{receipt}");
        Ok(())
    }
}

fn require_matching_complete_derivations(first: &str, second: &str) -> Result<(), String> {
    if first == second {
        Ok(())
    } else {
        Err(
            "dependency graph or package inputs changed between two complete receipt derivations; \
             the operator-observed filesystem view was not coherent and no receipt was emitted"
                .to_string(),
        )
    }
}

fn parse_command_args(raw_args: &[String]) -> Result<(bool, &[String]), String> {
    if raw_args.len() > MAX_SELECTION_ARGS + 2 {
        return Err(format!(
            "depgraph receipt command exceeds the {}-argument processing bound",
            MAX_SELECTION_ARGS + 2
        ));
    }
    let mut verify = false;
    for (index, arg) in raw_args.iter().enumerate() {
        if arg.len() > MAX_ARG_BYTES {
            return Err(format!(
                "argument {index} exceeds the {MAX_ARG_BYTES}-byte processing bound"
            ));
        }
        match arg.as_str() {
            "--verify" if !verify => verify = true,
            "--verify" => return Err("--verify may be supplied at most once".to_string()),
            "--" => return Ok((verify, &raw_args[index + 1..])),
            _ => {
                return Err(format!(
                    "unsupported receipt command option {arg:?}; put exactly one typed Cargo \
                     selection after `--` (for example `-- --package fs-roofline`)"
                ));
            }
        }
    }
    Err("missing `--` and explicit production root; no dependency-graph claim was made".to_string())
}

impl Selection {
    #[allow(clippy::too_many_lines)]
    fn parse(args: &[String]) -> Result<Self, String> {
        if args.is_empty() {
            return Err(
                "receipt selection requires exactly one explicit `--package <ROOT>`; no graph \
                 claim was made"
                    .to_string(),
            );
        }
        if args.len() > MAX_SELECTION_ARGS {
            return Err(format!(
                "receipt selection exceeds the {MAX_SELECTION_ARGS}-argument processing bound"
            ));
        }
        let mut package = None;
        let mut target = None;
        let mut features = BTreeSet::new();
        let mut all_features = false;
        let mut default_features = true;
        let mut saw_all_features = false;
        let mut saw_no_default_features = false;
        let mut index = 0;
        while index < args.len() {
            let arg = &args[index];
            let (flag, inline_value) = split_long_value(arg);
            match flag {
                "-p" | "--package" => {
                    let value = take_value(args, &mut index, flag, inline_value)?;
                    if package.is_some() {
                        return Err(
                            "receipt v1 supports exactly one explicit production root package; \
                             repeated/multiple --package selection is ambiguous and makes no claim"
                                .to_string(),
                        );
                    }
                    validate_machine_token("package", value, 128)?;
                    package = Some(value.to_string());
                }
                "--target" => {
                    let value = take_value(args, &mut index, flag, inline_value)?;
                    if target.is_some() {
                        return Err(
                            "receipt v1 supports exactly one target triple; repeated --target is \
                             ambiguous and makes no claim"
                                .to_string(),
                        );
                    }
                    validate_machine_token("target triple", value, 256)?;
                    target = Some(value.to_string());
                }
                "-F" | "--features" => {
                    let value = take_value(args, &mut index, flag, inline_value)?;
                    parse_features(value, &mut features)?;
                }
                "--all-features" if inline_value.is_none() => {
                    if saw_all_features {
                        return Err("--all-features may be supplied at most once".to_string());
                    }
                    saw_all_features = true;
                    all_features = true;
                }
                "--no-default-features" if inline_value.is_none() => {
                    if saw_no_default_features {
                        return Err(
                            "--no-default-features may be supplied at most once".to_string()
                        );
                    }
                    saw_no_default_features = true;
                    default_features = false;
                }
                _ => return Err(unsupported_selection(arg)),
            }
            index += 1;
        }
        if all_features && !features.is_empty() {
            return Err(
                "--all-features combined with explicit --features is redundant/ambiguous in \
                 receipt v1; no graph claim was made"
                    .to_string(),
            );
        }
        let package = package.ok_or_else(|| {
            "receipt selection requires exactly one explicit `--package <ROOT>`; workspace/default \
             selection is not a supported graph claim"
                .to_string()
        })?;
        Ok(Self {
            package,
            target,
            features,
            all_features,
            default_features,
        })
    }
}

fn split_long_value(arg: &str) -> (&str, Option<&str>) {
    if arg.starts_with("--") {
        arg.split_once('=')
            .map_or((arg, None), |(flag, value)| (flag, Some(value)))
    } else if let Some(value) = arg.strip_prefix("-p") {
        if value.is_empty() {
            (arg, None)
        } else {
            ("-p", Some(value))
        }
    } else if let Some(value) = arg.strip_prefix("-F") {
        if value.is_empty() {
            (arg, None)
        } else {
            ("-F", Some(value))
        }
    } else {
        (arg, None)
    }
}

fn take_value<'a>(
    args: &'a [String],
    index: &mut usize,
    flag: &str,
    inline: Option<&'a str>,
) -> Result<&'a str, String> {
    let value = if let Some(value) = inline {
        value
    } else {
        *index += 1;
        args.get(*index)
            .ok_or_else(|| format!("{flag} requires a value"))?
    };
    if value.is_empty() {
        return Err(format!("{flag} requires a non-empty value"));
    }
    if value.len() > MAX_ARG_BYTES {
        return Err(format!("{flag} value exceeds {MAX_ARG_BYTES} bytes"));
    }
    Ok(value)
}

fn validate_machine_token(label: &str, value: &str, max: usize) -> Result<(), String> {
    if value.is_empty()
        || value.len() > max
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'.'))
    {
        return Err(format!(
            "{label} must be 1..={max} ASCII machine characters [A-Za-z0-9_.-], got {value:?}"
        ));
    }
    Ok(())
}

fn parse_features(value: &str, features: &mut BTreeSet<String>) -> Result<(), String> {
    let mut found = false;
    for feature in value.split(|character: char| character == ',' || character.is_whitespace()) {
        if feature.is_empty() {
            continue;
        }
        found = true;
        if feature.len() > MAX_FEATURE_BYTES
            || !feature.bytes().all(|byte| {
                byte.is_ascii_alphanumeric()
                    || matches!(byte, b'_' | b'-' | b'+' | b'.' | b'/' | b'?')
            })
        {
            return Err(format!(
                "feature must be at most {MAX_FEATURE_BYTES} ASCII Cargo feature characters, got \
                 {feature:?}"
            ));
        }
        features.insert(feature.to_string());
        if features.len() > MAX_FEATURES_PER_PACKAGE {
            return Err(format!(
                "feature selection exceeds the {MAX_FEATURES_PER_PACKAGE}-feature processing bound"
            ));
        }
    }
    if !found {
        return Err("--features requires at least one feature".to_string());
    }
    Ok(())
}

fn unsupported_selection(arg: &str) -> String {
    let category =
        if matches!(arg, "--workspace" | "--all" | "--exclude") || arg.starts_with("--exclude=") {
            "workspace/multi-root selection"
        } else if matches!(
            arg,
            "--test"
                | "--tests"
                | "--bench"
                | "--benches"
                | "--example"
                | "--examples"
                | "--all-targets"
        ) || arg.starts_with("--test=")
            || arg.starts_with("--bench=")
            || arg.starts_with("--example=")
        {
            "test/dev/all-target selection"
        } else if matches!(arg, "--bin" | "--bins" | "--lib") || arg.starts_with("--bin=") {
            "target-kind selection"
        } else if matches!(arg, "--release" | "--profile") || arg.starts_with("--profile=") {
            "profile selection"
        } else {
            "unrecognized selection"
        };
    format!(
        "receipt v1 refuses {category} flag {arg:?}; it supports one production root plus only \
         graph-relevant target/features, so no dependency-graph claim was made"
    )
}

struct CollectedOutput {
    status: ExitStatus,
    stdout: Vec<u8>,
    stderr: Vec<u8>,
}

struct CappedRead {
    bytes: Vec<u8>,
    total: u64,
    overflowed: bool,
}

fn drain_capped(mut reader: impl Read, limit: usize) -> Result<CappedRead, String> {
    let mut bytes = Vec::with_capacity(limit.min(64 * 1_024));
    let mut buffer = vec![0u8; 64 * 1_024].into_boxed_slice();
    let mut total = 0u64;
    loop {
        let read = reader
            .read(&mut buffer)
            .map_err(|error| format!("cannot read subprocess pipe: {error}"))?;
        if read == 0 {
            break;
        }
        total = total
            .checked_add(read as u64)
            .ok_or_else(|| "subprocess byte count overflow".to_string())?;
        if bytes.len() < limit {
            let retained = (limit - bytes.len()).min(read);
            bytes.extend_from_slice(&buffer[..retained]);
        }
    }
    Ok(CappedRead {
        bytes,
        total,
        overflowed: total > limit as u64,
    })
}

fn run_capped(
    command: &mut Command,
    context: &str,
    stdout_limit: usize,
    stderr_limit: usize,
) -> Result<CollectedOutput, String> {
    let mut child = command
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|error| format!("cannot execute {context}: {error}"))?;
    let stdout = child
        .stdout
        .take()
        .ok_or_else(|| format!("{context} stdout pipe was not created"))?;
    let stderr = child
        .stderr
        .take()
        .ok_or_else(|| format!("{context} stderr pipe was not created"))?;
    let stdout_thread = thread::spawn(move || drain_capped(stdout, stdout_limit));
    let stderr_thread = thread::spawn(move || drain_capped(stderr, stderr_limit));
    let status = child
        .wait()
        .map_err(|error| format!("cannot wait for {context}: {error}"))?;
    let stdout = stdout_thread
        .join()
        .map_err(|_| format!("{context} stdout collector panicked"))??;
    let stderr = stderr_thread
        .join()
        .map_err(|_| format!("{context} stderr collector panicked"))??;
    if stdout.overflowed {
        return Err(format!(
            "{context} stdout emitted {} bytes, exceeding the {stdout_limit}-byte collection bound",
            stdout.total
        ));
    }
    if stderr.overflowed {
        return Err(format!(
            "{context} stderr emitted {} bytes, exceeding the {stderr_limit}-byte collection bound",
            stderr.total
        ));
    }
    Ok(CollectedOutput {
        status,
        stdout: stdout.bytes,
        stderr: stderr.bytes,
    })
}

fn resolve_executable(value: &str) -> Result<PathBuf, String> {
    let declared = PathBuf::from(value);
    let candidate = if declared.is_absolute() || declared.components().count() > 1 {
        declared
    } else {
        std::env::var_os("PATH")
            .into_iter()
            .flat_map(|path| std::env::split_paths(&path).collect::<Vec<_>>())
            .map(|directory| directory.join(&declared))
            .find(|path| path.is_file())
            .ok_or_else(|| format!("cannot resolve Cargo executable {value:?} on PATH"))?
    };
    candidate.canonicalize().map_err(|error| {
        format!(
            "cannot canonicalize Cargo executable {}: {error}",
            candidate.display()
        )
    })
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct FileIdentity {
    len: u64,
    readonly: bool,
    #[cfg(unix)]
    device: u64,
    #[cfg(unix)]
    inode: u64,
    #[cfg(unix)]
    mode: u32,
    #[cfg(unix)]
    modified_seconds: i64,
    #[cfg(unix)]
    modified_nanoseconds: i64,
    #[cfg(unix)]
    changed_seconds: i64,
    #[cfg(unix)]
    changed_nanoseconds: i64,
    #[cfg(not(unix))]
    modified: std::time::SystemTime,
}

#[allow(clippy::unnecessary_wraps)] // Non-Unix metadata modification times are fallible.
fn file_identity(metadata: &std::fs::Metadata) -> Result<FileIdentity, String> {
    Ok(FileIdentity {
        len: metadata.len(),
        readonly: metadata.permissions().readonly(),
        #[cfg(unix)]
        device: metadata.dev(),
        #[cfg(unix)]
        inode: metadata.ino(),
        #[cfg(unix)]
        mode: metadata.mode(),
        #[cfg(unix)]
        modified_seconds: metadata.mtime(),
        #[cfg(unix)]
        modified_nanoseconds: metadata.mtime_nsec(),
        #[cfg(unix)]
        changed_seconds: metadata.ctime(),
        #[cfg(unix)]
        changed_nanoseconds: metadata.ctime_nsec(),
        #[cfg(not(unix))]
        modified: metadata
            .modified()
            .map_err(|error| format!("cannot observe file modification time: {error}"))?,
    })
}

fn ensure_same_file(
    path: &Path,
    expected: &FileIdentity,
    observed: &std::fs::Metadata,
    phase: &str,
) -> Result<(), String> {
    let observed = file_identity(observed)?;
    if &observed == expected {
        Ok(())
    } else {
        Err(format!(
            "{} changed during {phase}: expected {expected:?}, observed {observed:?}",
            path.display()
        ))
    }
}

fn hash_open_file(
    path: &Path,
    file: &mut File,
    expected: &FileIdentity,
    max_bytes: u64,
) -> Result<(String, u64), String> {
    let mut hasher = fs_blake3::Blake3::new();
    let mut buffer = vec![0u8; 64 * 1_024].into_boxed_slice();
    let mut total = 0u64;
    loop {
        let read = file
            .read(&mut buffer)
            .map_err(|error| format!("cannot read {}: {error}", path.display()))?;
        if read == 0 {
            break;
        }
        total = total
            .checked_add(read as u64)
            .ok_or_else(|| format!("byte count overflow for {}", path.display()))?;
        if total > max_bytes {
            return Err(format!(
                "{} changed beyond its byte bound while hashing",
                path.display()
            ));
        }
        hasher.update(&buffer[..read]);
    }
    if total != expected.len {
        return Err(format!("{} changed size while hashing", path.display()));
    }
    let after = file
        .metadata()
        .map_err(|error| format!("cannot re-inspect open file {}: {error}", path.display()))?;
    ensure_same_file(path, expected, &after, "open-handle hashing")?;
    Ok((hasher.finalize().to_string(), total))
}

fn hash_regular_file(path: &Path, max_bytes: u64) -> Result<(String, u64), String> {
    let path_before = std::fs::symlink_metadata(path)
        .map_err(|error| format!("cannot inspect {}: {error}", path.display()))?;
    if !path_before.file_type().is_file() || path_before.len() > max_bytes {
        return Err(format!(
            "{} must be a non-symlink regular file no larger than {max_bytes} bytes",
            path.display()
        ));
    }
    let expected = file_identity(&path_before)?;
    let mut file =
        File::open(path).map_err(|error| format!("cannot open {}: {error}", path.display()))?;
    let opened = file
        .metadata()
        .map_err(|error| format!("cannot inspect open file {}: {error}", path.display()))?;
    ensure_same_file(path, &expected, &opened, "path opening")?;
    let hashed = hash_open_file(path, &mut file, &expected, max_bytes)?;
    let path_after = std::fs::symlink_metadata(path)
        .map_err(|error| format!("cannot re-inspect {}: {error}", path.display()))?;
    if !path_after.file_type().is_file() {
        return Err(format!(
            "{} changed file type while hashing",
            path.display()
        ));
    }
    ensure_same_file(path, &expected, &path_after, "path hashing")?;
    Ok(hashed)
}

fn ensure_cargo_identity_stable(
    before: &receipt_format::CargoIdentity,
    after: &receipt_format::CargoIdentity,
) -> Result<(), String> {
    if before == after {
        Ok(())
    } else {
        Err(format!(
            "Cargo executable identity moved during dependency observation: before={before:?}, after={after:?}"
        ))
    }
}

fn cargo_identity(cargo: &Path) -> Result<receipt_format::CargoIdentity, String> {
    let (executable_digest, _) = hash_regular_file(cargo, MAX_CARGO_EXECUTABLE_BYTES)?;
    let output = run_capped(
        Command::new(cargo).args(["--version", "--verbose"]),
        "cargo --version --verbose",
        64 * 1_024,
        MAX_SUBPROCESS_STDERR_BYTES,
    )?;
    if !output.status.success() {
        return Err(format!(
            "cargo --version --verbose failed: {}",
            String::from_utf8_lossy(&output.stderr)
        ));
    }
    let version = String::from_utf8(output.stdout)
        .map_err(|error| format!("Cargo version output is not UTF-8: {error}"))?
        .trim_end()
        .to_string();
    if version.is_empty() {
        return Err("Cargo version identity is empty".to_string());
    }
    let identity = receipt_format::CargoIdentity {
        executable_digest,
        version,
    };
    let (confirmed_digest, _) = hash_regular_file(cargo, MAX_CARGO_EXECUTABLE_BYTES)?;
    if confirmed_digest != identity.executable_digest {
        return Err(format!(
            "Cargo executable moved while its version was observed: before={}, after={confirmed_digest}",
            identity.executable_digest
        ));
    }
    Ok(identity)
}

#[derive(Debug)]
pub(super) enum JsonValue {
    Null,
    String(String),
    Array(Vec<JsonValue>),
    Object(BTreeMap<String, JsonValue>),
    Bool,
    Number(String),
}

pub(super) struct JsonParser<'a> {
    input: &'a str,
    cursor: usize,
    nodes: usize,
}

impl<'a> JsonParser<'a> {
    pub(super) fn new(input: &'a str) -> Self {
        Self {
            input,
            cursor: 0,
            nodes: 0,
        }
    }

    fn skip_space(&mut self) {
        while self
            .input
            .as_bytes()
            .get(self.cursor)
            .is_some_and(|byte| matches!(byte, b' ' | b'\t' | b'\r' | b'\n'))
        {
            self.cursor += 1;
        }
    }

    fn byte(&mut self) -> Result<u8, String> {
        let byte = *self
            .input
            .as_bytes()
            .get(self.cursor)
            .ok_or_else(|| "unexpected end of Cargo metadata JSON".to_string())?;
        self.cursor += 1;
        Ok(byte)
    }

    fn expect(&mut self, literal: &[u8]) -> Result<(), String> {
        if self
            .input
            .as_bytes()
            .get(self.cursor..self.cursor + literal.len())
            == Some(literal)
        {
            self.cursor += literal.len();
            Ok(())
        } else {
            Err(format!(
                "expected {:?} at Cargo metadata byte {}",
                String::from_utf8_lossy(literal),
                self.cursor
            ))
        }
    }

    fn hex4(&mut self) -> Result<u16, String> {
        let mut value = 0u16;
        for _ in 0..4 {
            let byte = self.byte()?;
            let digit = match byte {
                b'0'..=b'9' => u16::from(byte - b'0'),
                b'a'..=b'f' => u16::from(byte - b'a' + 10),
                b'A'..=b'F' => u16::from(byte - b'A' + 10),
                _ => return Err("invalid Cargo metadata Unicode escape".to_string()),
            };
            value = (value << 4) | digit;
        }
        Ok(value)
    }

    fn string(&mut self) -> Result<String, String> {
        self.expect(b"\"")?;
        let mut out = String::new();
        loop {
            let byte = self.byte()?;
            match byte {
                b'"' => break,
                b'\\' => match self.byte()? {
                    b'"' => out.push('"'),
                    b'\\' => out.push('\\'),
                    b'/' => out.push('/'),
                    b'b' => out.push('\u{08}'),
                    b'f' => out.push('\u{0c}'),
                    b'n' => out.push('\n'),
                    b'r' => out.push('\r'),
                    b't' => out.push('\t'),
                    b'u' => {
                        let high = self.hex4()?;
                        let codepoint = if (0xd800..=0xdbff).contains(&high) {
                            self.expect(b"\\u")?;
                            let low = self.hex4()?;
                            if !(0xdc00..=0xdfff).contains(&low) {
                                return Err("invalid Cargo metadata low surrogate".to_string());
                            }
                            0x1_0000
                                + ((u32::from(high) - 0xd800) << 10)
                                + (u32::from(low) - 0xdc00)
                        } else {
                            if (0xdc00..=0xdfff).contains(&high) {
                                return Err("unpaired Cargo metadata low surrogate".to_string());
                            }
                            u32::from(high)
                        };
                        out.push(
                            char::from_u32(codepoint)
                                .ok_or_else(|| "invalid Cargo metadata scalar".to_string())?,
                        );
                    }
                    _ => return Err("invalid Cargo metadata JSON escape".to_string()),
                },
                0x00..=0x1f => return Err("raw control in Cargo metadata string".to_string()),
                0x20..=0x7f => out.push(char::from(byte)),
                _ => {
                    self.cursor -= 1;
                    let character = self.input[self.cursor..]
                        .chars()
                        .next()
                        .ok_or_else(|| "invalid Cargo metadata UTF-8".to_string())?;
                    self.cursor += character.len_utf8();
                    out.push(character);
                }
            }
            if out.len() > receipt_format::MAX_STRING_BYTES * 4 {
                return Err("Cargo metadata string exceeds processing bound".to_string());
            }
        }
        Ok(out)
    }

    fn number(&mut self) -> Result<String, String> {
        let start = self.cursor;
        if self.input.as_bytes().get(self.cursor) == Some(&b'-') {
            self.cursor += 1;
        }
        match self.input.as_bytes().get(self.cursor).copied() {
            Some(b'0') => {
                self.cursor += 1;
                if self
                    .input
                    .as_bytes()
                    .get(self.cursor)
                    .is_some_and(u8::is_ascii_digit)
                {
                    return Err("invalid Cargo metadata number with a leading zero".to_string());
                }
            }
            Some(b'1'..=b'9') => {
                self.cursor += 1;
                while self
                    .input
                    .as_bytes()
                    .get(self.cursor)
                    .is_some_and(u8::is_ascii_digit)
                {
                    self.cursor += 1;
                }
            }
            _ => return Err("invalid Cargo metadata number".to_string()),
        }
        if self.input.as_bytes().get(self.cursor) == Some(&b'.') {
            self.cursor += 1;
            let fraction = self.cursor;
            while self
                .input
                .as_bytes()
                .get(self.cursor)
                .is_some_and(u8::is_ascii_digit)
            {
                self.cursor += 1;
            }
            if self.cursor == fraction {
                return Err("invalid Cargo metadata number fraction".to_string());
            }
        }
        if self
            .input
            .as_bytes()
            .get(self.cursor)
            .is_some_and(|byte| matches!(byte, b'e' | b'E'))
        {
            self.cursor += 1;
            if self
                .input
                .as_bytes()
                .get(self.cursor)
                .is_some_and(|byte| matches!(byte, b'+' | b'-'))
            {
                self.cursor += 1;
            }
            let exponent = self.cursor;
            while self
                .input
                .as_bytes()
                .get(self.cursor)
                .is_some_and(u8::is_ascii_digit)
            {
                self.cursor += 1;
            }
            if self.cursor == exponent {
                return Err("invalid Cargo metadata number exponent".to_string());
            }
        }
        Ok(self.input[start..self.cursor].to_string())
    }

    fn value(&mut self, depth: usize) -> Result<JsonValue, String> {
        if depth > 128 || self.nodes >= 2_000_000 {
            return Err("Cargo metadata JSON exceeds depth/node bounds".to_string());
        }
        self.nodes += 1;
        self.skip_space();
        match self.input.as_bytes().get(self.cursor).copied() {
            Some(b'"') => self.string().map(JsonValue::String),
            Some(b'[') => {
                self.cursor += 1;
                let mut values = Vec::new();
                self.skip_space();
                if self.input.as_bytes().get(self.cursor) == Some(&b']') {
                    self.cursor += 1;
                    return Ok(JsonValue::Array(values));
                }
                loop {
                    values.push(self.value(depth + 1)?);
                    self.skip_space();
                    match self.byte()? {
                        b',' => {}
                        b']' => break,
                        _ => return Err("invalid Cargo metadata array".to_string()),
                    }
                }
                Ok(JsonValue::Array(values))
            }
            Some(b'{') => {
                self.cursor += 1;
                let mut values = BTreeMap::new();
                self.skip_space();
                if self.input.as_bytes().get(self.cursor) == Some(&b'}') {
                    self.cursor += 1;
                    return Ok(JsonValue::Object(values));
                }
                loop {
                    self.skip_space();
                    let key = self.string()?;
                    self.skip_space();
                    self.expect(b":")?;
                    let value = self.value(depth + 1)?;
                    if values.insert(key.clone(), value).is_some() {
                        return Err(format!("duplicate Cargo metadata key {key:?}"));
                    }
                    self.skip_space();
                    match self.byte()? {
                        b',' => {}
                        b'}' => break,
                        _ => return Err("invalid Cargo metadata object".to_string()),
                    }
                }
                Ok(JsonValue::Object(values))
            }
            Some(b'n') => {
                self.expect(b"null")?;
                Ok(JsonValue::Null)
            }
            Some(b't') => {
                self.expect(b"true")?;
                Ok(JsonValue::Bool)
            }
            Some(b'f') => {
                self.expect(b"false")?;
                Ok(JsonValue::Bool)
            }
            Some(b'-' | b'0'..=b'9') => self.number().map(JsonValue::Number),
            _ => Err(format!(
                "unexpected Cargo metadata token at byte {}",
                self.cursor
            )),
        }
    }

    pub(super) fn finish(mut self) -> Result<JsonValue, String> {
        let value = self.value(0)?;
        self.skip_space();
        if self.cursor != self.input.len() {
            return Err("trailing bytes after Cargo metadata JSON".to_string());
        }
        Ok(value)
    }
}

fn json_object(value: &JsonValue) -> Result<&BTreeMap<String, JsonValue>, String> {
    if let JsonValue::Object(value) = value {
        Ok(value)
    } else {
        Err("Cargo metadata field is not an object".to_string())
    }
}

fn json_array(value: &JsonValue) -> Result<&[JsonValue], String> {
    if let JsonValue::Array(value) = value {
        Ok(value)
    } else {
        Err("Cargo metadata field is not an array".to_string())
    }
}

fn json_string(value: &JsonValue) -> Result<&str, String> {
    if let JsonValue::String(value) = value {
        Ok(value)
    } else {
        Err("Cargo metadata field is not a string".to_string())
    }
}

#[derive(Debug, Clone)]
struct MetadataPackage {
    metadata_id: String,
    name: String,
    version: String,
    source_id: Option<String>,
    manifest_dir: PathBuf,
}

struct MetadataIndex {
    packages: Vec<MetadataPackage>,
    by_manifest_dir: BTreeMap<PathBuf, usize>,
    by_name_version: BTreeMap<(String, String), Vec<usize>>,
}

fn parse_metadata(text: &str) -> Result<MetadataIndex, String> {
    let root = JsonParser::new(text).finish()?;
    let object = json_object(&root)?;
    let packages = json_array(
        object
            .get("packages")
            .ok_or_else(|| "Cargo metadata lacks packages".to_string())?,
    )?;
    if packages.len() > 100_000 {
        return Err("Cargo metadata package count exceeds 100000".to_string());
    }
    let mut parsed = Vec::with_capacity(packages.len());
    let mut by_manifest_dir = BTreeMap::new();
    let mut by_name_version: BTreeMap<(String, String), Vec<usize>> = BTreeMap::new();
    for package in packages {
        let package = json_object(package)?;
        let metadata_id = json_string(
            package
                .get("id")
                .ok_or_else(|| "Cargo metadata package lacks id".to_string())?,
        )?
        .to_string();
        let name = json_string(
            package
                .get("name")
                .ok_or_else(|| "Cargo metadata package lacks name".to_string())?,
        )?
        .to_string();
        let version = json_string(
            package
                .get("version")
                .ok_or_else(|| "Cargo metadata package lacks version".to_string())?,
        )?
        .to_string();
        let source_id = match package
            .get("source")
            .ok_or_else(|| "Cargo metadata package lacks source".to_string())?
        {
            JsonValue::Null => None,
            value => Some(json_string(value)?.to_string()),
        };
        let manifest_path =
            PathBuf::from(json_string(package.get("manifest_path").ok_or_else(
                || "Cargo metadata package lacks manifest_path".to_string(),
            )?)?);
        let manifest_dir = manifest_path
            .parent()
            .ok_or_else(|| {
                format!(
                    "metadata manifest has no parent: {}",
                    manifest_path.display()
                )
            })?
            .canonicalize()
            .map_err(|error| {
                format!(
                    "cannot canonicalize metadata package directory {}: {error}",
                    manifest_path.display()
                )
            })?;
        let index = parsed.len();
        if by_manifest_dir
            .insert(manifest_dir.clone(), index)
            .is_some()
        {
            return Err(format!(
                "Cargo metadata maps multiple packages to {}",
                manifest_dir.display()
            ));
        }
        by_name_version
            .entry((name.clone(), version.clone()))
            .or_default()
            .push(index);
        parsed.push(MetadataPackage {
            metadata_id,
            name,
            version,
            source_id,
            manifest_dir,
        });
    }
    Ok(MetadataIndex {
        packages: parsed,
        by_manifest_dir,
        by_name_version,
    })
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PackageFile {
    path: PathBuf,
    relative: String,
    symlink_target: Option<String>,
    bytes: u64,
    content_digest: String,
    identity: FileIdentity,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PackageSnapshot {
    files: Vec<PackageFile>,
    directories: Vec<String>,
    bytes: u64,
    manifest_bytes: usize,
}

#[derive(Debug)]
struct PathClosureBudget {
    read_bytes: u64,
    entries: usize,
}

impl PathClosureBudget {
    const fn new() -> Self {
        Self {
            read_bytes: 0,
            entries: 0,
        }
    }

    fn charge_entry(&mut self, root: &Path) -> Result<(), String> {
        let next = self
            .entries
            .checked_add(1)
            .ok_or_else(|| "path-closure entry-count overflow".to_string())?;
        if next > MAX_PATH_CLOSURE_ENTRIES {
            return Err(format!(
                "path closure exceeds the {MAX_PATH_CLOSURE_ENTRIES}-entry aggregate work bound while scanning {}",
                root.display()
            ));
        }
        self.entries = next;
        Ok(())
    }

    fn remaining_read_bytes(&self) -> u64 {
        MAX_PATH_CLOSURE_READ_BYTES.saturating_sub(self.read_bytes)
    }

    fn charge_read(&mut self, root: &Path, bytes: u64) -> Result<(), String> {
        let next = self
            .read_bytes
            .checked_add(bytes)
            .ok_or_else(|| "path-closure byte-count overflow".to_string())?;
        if next > MAX_PATH_CLOSURE_READ_BYTES {
            return Err(format!(
                "path closure exceeds the {MAX_PATH_CLOSURE_READ_BYTES}-byte aggregate read bound while hashing {}",
                root.display()
            ));
        }
        self.read_bytes = next;
        Ok(())
    }
}

fn excluded_build_entry(root: &Path, path: &Path) -> bool {
    path.parent() == Some(root)
        && path
            .file_name()
            .and_then(|name| name.to_str())
            .is_some_and(|name| matches!(name, ".git" | "target"))
}

fn points_into_excluded_root(relative: &Path) -> bool {
    relative
        .components()
        .next()
        .is_some_and(|component| matches!(component.as_os_str().to_str(), Some(".git" | "target")))
}

fn enqueue_directory(
    root: &Path,
    directory: &Path,
    depth: usize,
    scheduled: &mut usize,
    pending: &mut Vec<(PathBuf, usize)>,
) -> Result<(), String> {
    if depth > MAX_PATH_PACKAGE_DEPTH {
        return Err(format!(
            "path package {} exceeds the {MAX_PATH_PACKAGE_DEPTH}-level directory-depth bound at {}",
            root.display(),
            directory.display()
        ));
    }
    let next = scheduled
        .checked_add(1)
        .ok_or_else(|| "package directory-count overflow".to_string())?;
    if next > MAX_PATH_PACKAGE_DIRECTORIES {
        return Err(format!(
            "path package {} exceeds the {MAX_PATH_PACKAGE_DIRECTORIES}-directory processing bound",
            root.display()
        ));
    }
    pending
        .try_reserve(1)
        .map_err(|_| "cannot reserve bounded package-directory queue".to_string())?;
    pending.push((directory.to_path_buf(), depth));
    *scheduled = next;
    Ok(())
}

fn normalized_relative(root: &Path, path: &Path) -> Result<String, String> {
    let relative = path.strip_prefix(root).map_err(|error| {
        format!(
            "package input {} escaped root {}: {error}",
            path.display(),
            root.display()
        )
    })?;
    let relative = relative
        .to_str()
        .ok_or_else(|| format!("package input path is not Unicode: {}", relative.display()))?;
    Ok(relative.replace(std::path::MAIN_SEPARATOR, "/"))
}

fn observe_package_file(
    root: &Path,
    path: &Path,
) -> Result<(Option<String>, FileIdentity), String> {
    let metadata = std::fs::symlink_metadata(path)
        .map_err(|error| format!("cannot inspect package input {}: {error}", path.display()))?;
    if metadata.file_type().is_symlink() {
        let raw_target = std::fs::read_link(path)
            .map_err(|error| format!("cannot read package symlink {}: {error}", path.display()))?;
        let raw_target = raw_target
            .to_str()
            .ok_or_else(|| format!("package symlink target is not Unicode: {}", path.display()))?;
        let target = path.canonicalize().map_err(|error| {
            format!("cannot resolve package symlink {}: {error}", path.display())
        })?;
        if !target.starts_with(root) {
            return Err(format!(
                "package symlink {} escapes source root {}",
                path.display(),
                root.display()
            ));
        }
        let relative_target = normalized_relative(root, &target)?;
        if points_into_excluded_root(Path::new(&relative_target)) {
            return Err(format!(
                "package symlink {} points into excluded VCS/build output",
                path.display()
            ));
        }
        let target_metadata = target.metadata().map_err(|error| {
            format!(
                "cannot inspect symlink target {}: {error}",
                target.display()
            )
        })?;
        if !target_metadata.is_file() {
            return Err(format!(
                "package symlink {} does not resolve to a regular file",
                path.display()
            ));
        }
        let raw_target = raw_target.replace(std::path::MAIN_SEPARATOR, "/");
        Ok((
            Some(format!("raw={raw_target}\nresolved={relative_target}")),
            file_identity(&target_metadata)?,
        ))
    } else if metadata.is_file() {
        Ok((None, file_identity(&metadata)?))
    } else {
        Err(format!(
            "package input {} is neither a directory nor a regular file",
            path.display()
        ))
    }
}

#[allow(clippy::too_many_lines)] // One bounded, identity-checked package-tree snapshot.
fn capture_package_tree(
    root: &Path,
    closure_budget: &mut PathClosureBudget,
) -> Result<PackageSnapshot, String> {
    let mut pending = Vec::new();
    let mut scheduled_directories = 0usize;
    enqueue_directory(root, root, 0, &mut scheduled_directories, &mut pending)?;
    closure_budget.charge_entry(root)?;
    let mut files = Vec::new();
    let mut directories = Vec::new();
    let mut total_bytes = 0u64;
    while let Some((directory, depth)) = pending.pop() {
        let directory_metadata = std::fs::symlink_metadata(&directory).map_err(|error| {
            format!("cannot inspect directory {}: {error}", directory.display())
        })?;
        if !directory_metadata.file_type().is_dir() {
            return Err(format!(
                "package directory {} changed file type during enumeration",
                directory.display()
            ));
        }
        let resolved = directory.canonicalize().map_err(|error| {
            format!(
                "cannot resolve package directory {}: {error}",
                directory.display()
            )
        })?;
        if !resolved.starts_with(root) {
            return Err(format!(
                "package directory {} escaped source root {}",
                directory.display(),
                root.display()
            ));
        }
        let relative_directory = normalized_relative(root, &directory)?;
        directories
            .try_reserve(1)
            .map_err(|_| "cannot reserve bounded package-directory manifest".to_string())?;
        directories.push(relative_directory);
        for entry in std::fs::read_dir(&directory)
            .map_err(|error| format!("cannot enumerate {}: {error}", directory.display()))?
        {
            let entry = entry.map_err(|error| {
                format!("cannot read entry in {}: {error}", directory.display())
            })?;
            closure_budget.charge_entry(root)?;
            let path = entry.path();
            if excluded_build_entry(root, &path) {
                continue;
            }
            let metadata = std::fs::symlink_metadata(&path)
                .map_err(|error| format!("cannot inspect {}: {error}", path.display()))?;
            if metadata.file_type().is_dir() {
                enqueue_directory(
                    root,
                    &path,
                    depth.saturating_add(1),
                    &mut scheduled_directories,
                    &mut pending,
                )?;
                continue;
            }
            if files.len() >= MAX_PATH_PACKAGE_FILES {
                return Err(format!(
                    "path package {} changed beyond the {MAX_PATH_PACKAGE_FILES}-file collection bound",
                    root.display()
                ));
            }
            let remaining = MAX_PATH_PACKAGE_BYTES.saturating_sub(total_bytes);
            let (symlink_target, expected_identity) = observe_package_file(root, &path)?;
            if expected_identity.len > remaining {
                return Err(format!(
                    "path package {} exceeds the {MAX_PATH_PACKAGE_BYTES}-byte collection bound",
                    root.display()
                ));
            }
            let aggregate_remaining = closure_budget.remaining_read_bytes();
            if expected_identity.len > aggregate_remaining {
                return Err(format!(
                    "path closure exceeds the {MAX_PATH_CLOSURE_READ_BYTES}-byte aggregate read bound before hashing {}",
                    path.display()
                ));
            }
            let mut opened = File::open(&path).map_err(|error| {
                format!("cannot open package input {}: {error}", path.display())
            })?;
            let opened_metadata = opened.metadata().map_err(|error| {
                format!(
                    "cannot inspect open package input {}: {error}",
                    path.display()
                )
            })?;
            ensure_same_file(
                &path,
                &expected_identity,
                &opened_metadata,
                "package opening",
            )?;
            let (content_digest, bytes) = hash_open_file(
                &path,
                &mut opened,
                &expected_identity,
                remaining.min(aggregate_remaining),
            )?;
            closure_budget.charge_read(root, bytes)?;
            let (confirmed_symlink_target, confirmed_identity) = observe_package_file(root, &path)?;
            if confirmed_symlink_target != symlink_target || confirmed_identity != expected_identity
            {
                return Err(format!(
                    "package input {} changed path, symlink target, or identity while hashing",
                    path.display()
                ));
            }
            total_bytes = total_bytes
                .checked_add(bytes)
                .ok_or_else(|| "package collection byte-count overflow".to_string())?;
            let relative = normalized_relative(root, &path)?;
            files
                .try_reserve(1)
                .map_err(|_| "cannot reserve bounded package-file manifest".to_string())?;
            files.push(PackageFile {
                relative,
                path,
                symlink_target,
                bytes,
                content_digest,
                identity: expected_identity,
            });
        }
    }
    files.sort_by(|left, right| left.relative.cmp(&right.relative));
    directories.sort();
    let manifest_bytes = files
        .iter()
        .try_fold(0usize, |total, file| {
            total.checked_add(
                file.relative.len()
                    + file.symlink_target.as_ref().map_or(0, String::len)
                    + file.content_digest.len()
                    + 128,
            )
        })
        .and_then(|total| {
            directories.iter().try_fold(total, |total, directory| {
                total.checked_add(directory.len() + 32)
            })
        })
        .ok_or_else(|| "package manifest-size overflow".to_string())?;
    if manifest_bytes > MAX_PATH_MANIFEST_BYTES {
        return Err(format!(
            "path package {} exceeds the {MAX_PATH_MANIFEST_BYTES}-byte manifest bound",
            root.display()
        ));
    }
    Ok(PackageSnapshot {
        files,
        directories,
        bytes: total_bytes,
        manifest_bytes,
    })
}

fn ensure_matching_package_snapshots(
    root: &Path,
    before: &PackageSnapshot,
    after: &PackageSnapshot,
) -> Result<(), String> {
    if before == after {
        Ok(())
    } else {
        Err(format!(
            "path package {} changed between two complete bounded content snapshots",
            root.display()
        ))
    }
}

fn append_digest_field(payload: &mut Vec<u8>, name: &str, bytes: &[u8]) -> Result<(), String> {
    let name_len =
        u64::try_from(name.len()).map_err(|_| "field name length overflow".to_string())?;
    let byte_len =
        u64::try_from(bytes.len()).map_err(|_| "field byte length overflow".to_string())?;
    payload.extend_from_slice(&name_len.to_le_bytes());
    payload.extend_from_slice(name.as_bytes());
    payload.extend_from_slice(&byte_len.to_le_bytes());
    payload.extend_from_slice(bytes);
    Ok(())
}

fn hash_path_package(
    root: &Path,
    closure_budget: &mut PathClosureBudget,
) -> Result<String, String> {
    let root = root.canonicalize().map_err(|error| {
        format!(
            "cannot canonicalize path package {}: {error}",
            root.display()
        )
    })?;
    let first = capture_package_tree(&root, closure_budget)?;
    if first.files.is_empty() {
        return Err(format!(
            "path package {} contains no source inputs",
            root.display()
        ));
    }
    let second = capture_package_tree(&root, closure_budget)?;
    ensure_matching_package_snapshots(&root, &first, &second)?;
    let mut payload = Vec::new();
    payload
        .try_reserve_exact(first.manifest_bytes)
        .map_err(|_| "cannot reserve bounded path-package digest manifest".to_string())?;
    append_digest_field(&mut payload, "schema", b"path-package-source-build-v2")?;
    for directory in &first.directories {
        append_digest_field(&mut payload, "directory", directory.as_bytes())?;
    }
    for file in &first.files {
        append_digest_field(&mut payload, "path", file.relative.as_bytes())?;
        append_digest_field(
            &mut payload,
            "symlink-target",
            file.symlink_target.as_deref().unwrap_or("").as_bytes(),
        )?;
        append_digest_field(&mut payload, "bytes", &file.bytes.to_le_bytes())?;
        append_digest_field(
            &mut payload,
            "content-blake3",
            file.content_digest.as_bytes(),
        )?;
    }
    if first.bytes != second.bytes || payload.len() > MAX_PATH_MANIFEST_BYTES {
        return Err(format!(
            "path package {} changed or exceeded bounds",
            root.display()
        ));
    }
    Ok(fs_blake3::hash_domain(PATH_PACKAGE_DOMAIN, &payload).to_string())
}

fn derive_receipt(root: &Path, selection: &Selection) -> Result<String, String> {
    let cargo_value = std::env::var("CARGO").unwrap_or_else(|_| "cargo".to_string());
    let cargo = resolve_executable(&cargo_value)?;
    let cargo_before = cargo_identity(&cargo)?;
    let metadata_output = run_capped(
        Command::new(&cargo).current_dir(root).args([
            "metadata",
            "--format-version",
            "1",
            "--locked",
        ]),
        "cargo metadata --locked",
        MAX_METADATA_BYTES,
        MAX_SUBPROCESS_STDERR_BYTES,
    )?;
    if !metadata_output.status.success() {
        return Err(format!(
            "cargo metadata --locked failed: {}",
            String::from_utf8_lossy(&metadata_output.stderr)
        ));
    }
    let metadata_text = String::from_utf8(metadata_output.stdout)
        .map_err(|error| format!("Cargo metadata is not UTF-8: {error}"))?;
    let metadata = parse_metadata(&metadata_text)?;

    let mut tree_command = Command::new(&cargo);
    tree_command.current_dir(root).args([
        "tree",
        "--locked",
        "-e",
        "normal,build",
        "--no-dedupe",
        "--format",
        "{p}|{f}",
        "--prefix",
        "depth",
    ]);
    tree_command.args(selection.cargo_args());
    let tree_output = run_capped(
        &mut tree_command,
        "cargo tree --locked",
        MAX_TREE_BYTES,
        MAX_SUBPROCESS_STDERR_BYTES,
    )?;
    if !tree_output.status.success() {
        return Err(format!(
            "cargo tree --locked failed for normalized selection {:?}: {}",
            selection.cargo_args(),
            String::from_utf8_lossy(&tree_output.stderr)
        ));
    }
    let tree_text = String::from_utf8(tree_output.stdout)
        .map_err(|error| format!("cargo tree emitted non-UTF-8 output: {error}"))?;
    let cargo_after = cargo_identity(&cargo)?;
    ensure_cargo_identity_stable(&cargo_before, &cargo_after)?;
    let receipt = fs_la_closure(
        &tree_text,
        &selection.package,
        selection,
        cargo_before,
        &metadata,
    )?;
    receipt_format::emit(&receipt)
}

#[derive(Debug)]
struct TreeRow {
    depth: usize,
    name: String,
    version: String,
    source_hint: Option<String>,
    features: BTreeSet<String>,
}

fn parse_tree_row(line: &str) -> Result<TreeRow, String> {
    let digits = line.bytes().take_while(u8::is_ascii_digit).count();
    if digits == 0 {
        return Err(format!("cargo tree line has no depth prefix: {line:?}"));
    }
    let depth = line[..digits]
        .parse()
        .map_err(|error| format!("bad depth prefix in {line:?}: {error}"))?;
    let rest = &line[digits..];
    let (display, feature_csv) = rest
        .split_once('|')
        .ok_or_else(|| format!("cargo tree line missing feature separator: {line:?}"))?;
    let (name, version, source_hint) = parse_tree_package_display(display)?;
    let features = resolved_features(feature_csv)?;
    Ok(TreeRow {
        depth,
        name,
        version,
        source_hint,
        features,
    })
}

fn parse_tree_package_display(display: &str) -> Result<(String, String, Option<String>), String> {
    let mut fields = display.splitn(3, ' ');
    let name = fields
        .next()
        .filter(|value| !value.is_empty())
        .ok_or_else(|| format!("cargo tree package has no name: {display:?}"))?;
    validate_machine_token("resolved package name", name, 128)?;
    let version = fields
        .next()
        .and_then(|value| value.strip_prefix('v'))
        .filter(|value| !value.is_empty())
        .ok_or_else(|| format!("cargo tree package has no canonical version: {display:?}"))?;
    if version.len() > 256 || !version.bytes().all(|byte| byte.is_ascii_graphic()) {
        return Err(format!("invalid resolved package version in {display:?}"));
    }
    let mut source_display = fields.next().unwrap_or("").trim();
    if let Some(rest) = source_display.strip_prefix("(proc-macro)") {
        source_display = rest.trim();
    }
    let source_hint = if source_display.is_empty() {
        None
    } else {
        let source = source_display
            .strip_prefix('(')
            .and_then(|value| value.strip_suffix(')'))
            .ok_or_else(|| format!("unrecognized cargo tree package source: {display:?}"))?;
        if source.is_empty() {
            return Err(format!("empty cargo tree package source: {display:?}"));
        }
        Some(source.to_string())
    };
    Ok((name.to_string(), version.to_string(), source_hint))
}

fn resolve_tree_package(row: &TreeRow, metadata: &MetadataIndex) -> Result<usize, String> {
    if let Some(source) = &row.source_hint
        && !source.contains("://")
        && !source.starts_with("git+")
        && !source.starts_with("ssh:")
    {
        let directory = PathBuf::from(source).canonicalize().map_err(|error| {
            format!("cannot canonicalize cargo-tree package path {source:?}: {error}")
        })?;
        let index = *metadata.by_manifest_dir.get(&directory).ok_or_else(|| {
            format!(
                "cargo tree path {} has no exact cargo-metadata package identity",
                directory.display()
            )
        })?;
        let package = &metadata.packages[index];
        if package.name != row.name || package.version != row.version {
            return Err(format!(
                "cargo tree package {} {} disagrees with metadata identity {} {} at {}",
                row.name,
                row.version,
                package.name,
                package.version,
                directory.display()
            ));
        }
        return Ok(index);
    }
    let candidates = metadata
        .by_name_version
        .get(&(row.name.clone(), row.version.clone()))
        .ok_or_else(|| {
            format!(
                "cargo tree package {} {} has no cargo-metadata identity",
                row.name, row.version
            )
        })?;
    let matched: Vec<usize> = if let Some(source) = &row.source_hint {
        candidates
            .iter()
            .copied()
            .filter(|index| {
                let package = &metadata.packages[*index];
                package
                    .source_id
                    .as_ref()
                    .is_some_and(|identity| identity.contains(source))
                    || package.metadata_id.contains(source)
            })
            .collect()
    } else {
        candidates.clone()
    };
    if matched.len() != 1 {
        return Err(format!(
            "cargo tree display for {} {} source {:?} maps to {} metadata packages; refusing \
             collision-prone human identity",
            row.name,
            row.version,
            row.source_hint,
            matched.len()
        ));
    }
    Ok(matched[0])
}

fn resolved_features(csv: &str) -> Result<BTreeSet<String>, String> {
    let csv = csv.strip_suffix(" (*)").unwrap_or(csv);
    let mut features = BTreeSet::new();
    if csv.is_empty() {
        return Ok(features);
    }
    for feature in csv.split(',') {
        if feature.is_empty() || feature.len() > MAX_FEATURE_BYTES {
            return Err(format!("invalid resolved feature {feature:?}"));
        }
        if !feature.bytes().all(|byte| byte.is_ascii_graphic()) {
            return Err(format!("non-ASCII resolved feature {feature:?}"));
        }
        features.insert(feature.to_string());
        if features.len() > MAX_FEATURES_PER_PACKAGE {
            return Err(format!(
                "resolved package exceeds the {MAX_FEATURES_PER_PACKAGE}-feature processing bound"
            ));
        }
    }
    Ok(features)
}

fn package_identity(
    index: usize,
    metadata: &MetadataIndex,
    path_digests: &BTreeMap<PathBuf, String>,
) -> Result<receipt_format::PackageIdentity, String> {
    let package = metadata
        .packages
        .get(index)
        .ok_or_else(|| format!("metadata package index {index} is out of bounds"))?;
    if let Some(source_id) = &package.source_id {
        return Ok(receipt_format::PackageIdentity {
            name: package.name.clone(),
            version: package.version.clone(),
            package_id: package.metadata_id.clone(),
            source_id: Some(source_id.clone()),
            path_digest: None,
        });
    }
    let digest = path_digests.get(&package.manifest_dir).ok_or_else(|| {
        format!(
            "local package {} was not included in the preflighted path-package set",
            package.manifest_dir.display()
        )
    })?;
    Ok(receipt_format::PackageIdentity {
        name: package.name.clone(),
        version: package.version.clone(),
        package_id: format!("path+blake3:{digest}#{}@{}", package.name, package.version),
        source_id: None,
        path_digest: Some(digest.clone()),
    })
}

fn distinct_path_package_roots(
    root_index: usize,
    units: &BTreeSet<(usize, BTreeSet<String>)>,
    metadata: &MetadataIndex,
) -> Result<BTreeSet<PathBuf>, String> {
    let mut package_indices = BTreeSet::new();
    package_indices.insert(root_index);
    package_indices.extend(units.iter().map(|(index, _)| *index));

    let mut roots = BTreeSet::new();
    for index in package_indices {
        let package = metadata
            .packages
            .get(index)
            .ok_or_else(|| format!("metadata package index {index} is out of bounds"))?;
        if package.source_id.is_none() {
            roots.insert(package.manifest_dir.clone());
        }
    }
    Ok(roots)
}

fn hash_distinct_path_packages_with<F>(
    roots: &BTreeSet<PathBuf>,
    mut hash: F,
) -> Result<BTreeMap<PathBuf, String>, String>
where
    F: FnMut(&Path, &mut PathClosureBudget) -> Result<String, String>,
{
    if roots.len() > MAX_PATH_CLOSURE_PACKAGES {
        return Err(format!(
            "dependency closure contains {} distinct path packages, exceeding the {MAX_PATH_CLOSURE_PACKAGES}-package aggregate hashing bound",
            roots.len()
        ));
    }
    let mut digests = BTreeMap::new();
    let mut budget = PathClosureBudget::new();
    for root in roots {
        let digest = hash(root, &mut budget)?;
        digests.insert(root.clone(), digest);
    }
    Ok(digests)
}

fn hash_distinct_path_packages(
    roots: &BTreeSet<PathBuf>,
) -> Result<BTreeMap<PathBuf, String>, String> {
    hash_distinct_path_packages_with(roots, hash_path_package)
}

/// Extract metadata-resolved unit variants in every fs-la-rooted subtree.
#[allow(clippy::too_many_arguments)]
fn fs_la_closure(
    tree: &str,
    selected_root: &str,
    selection: &Selection,
    cargo: receipt_format::CargoIdentity,
    metadata: &MetadataIndex,
) -> Result<receipt_format::Receipt, String> {
    let mut units: BTreeSet<(usize, BTreeSet<String>)> = BTreeSet::new();
    let mut root: Option<(usize, BTreeSet<String>)> = None;
    let mut in_subtree_at: Option<usize> = None;
    let mut line_count = 0;
    for line in tree.lines() {
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }
        line_count += 1;
        if line_count > MAX_TREE_LINES {
            return Err(format!(
                "cargo tree output exceeds the {MAX_TREE_LINES}-line processing bound"
            ));
        }
        let row = parse_tree_row(line)?;
        let package_index = resolve_tree_package(&row, metadata)?;
        if row.depth == 0 {
            if root.is_some() {
                return Err(
                    "cargo tree emitted multiple depth-zero roots for a single-root selection; \
                     refusing ambiguous evidence"
                        .to_string(),
                );
            }
            if row.name != selected_root {
                return Err(format!(
                    "cargo tree root {:?} does not match explicit selection {selected_root:?}",
                    row.name
                ));
            }
            root = Some((package_index, row.features.clone()));
        }
        if let Some(root_depth) = in_subtree_at
            && row.depth <= root_depth
        {
            in_subtree_at = None;
        }
        if row.name == CLOSURE_ROOT && in_subtree_at.is_none() {
            in_subtree_at = Some(row.depth);
        }
        if in_subtree_at.is_some() {
            if units.len() >= receipt_format::MAX_PACKAGES
                && !units.contains(&(package_index, row.features.clone()))
            {
                return Err(format!(
                    "fs-la closure exceeds the {} package-unit bound",
                    receipt_format::MAX_PACKAGES
                ));
            }
            units.insert((package_index, row.features));
        }
    }
    let (root_index, root_features) =
        root.ok_or_else(|| "cargo tree emitted no depth-zero root".to_string())?;
    if units.is_empty() {
        return Err(format!(
            "production root {selected_root:?} does not reach fs-la through normal/build edges; \
             no depgraph receipt applies"
        ));
    }
    // Discover and admit the entire distinct local-package set before the first
    // source byte is read. package_identity is deliberately lookup-only.
    let path_roots = distinct_path_package_roots(root_index, &units, metadata)?;
    let path_digests = hash_distinct_path_packages(&path_roots)?;
    let root_identity = package_identity(root_index, metadata, &path_digests)?;
    let mut packages = Vec::with_capacity(units.len());
    for (index, features) in units {
        packages.push(receipt_format::PackageRow {
            identity: package_identity(index, metadata, &path_digests)?,
            features,
        });
    }
    packages.sort();
    let receipt = receipt_format::Receipt {
        cargo,
        root: receipt_format::RootRow {
            identity: root_identity,
            features: root_features,
        },
        selection: receipt_format::SelectionRow {
            target: selection.target.clone(),
            features: selection.features.clone(),
            all_features: selection.all_features,
            default_features: selection.default_features,
        },
        packages,
    };
    receipt_format::validate(&receipt)?;
    Ok(receipt)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn selection_normalizes_features_and_inline_forms() {
        let selection = Selection::parse(&[
            "--package=fs-roofline".into(),
            "--target".into(),
            "x86_64-unknown-linux-gnu".into(),
            "--features=z,a,a".into(),
            "-Fa/b".into(),
            "--no-default-features".into(),
        ])
        .expect("supported selection");
        assert_eq!(selection.package, "fs-roofline");
        assert_eq!(
            selection.features.into_iter().collect::<Vec<_>>(),
            ["a", "a/b", "z"]
        );
        assert!(!selection.default_features);
    }

    #[test]
    fn complete_derivation_guard_refuses_a_mixed_filesystem_view() {
        let stable = "canonical-receipt".to_string();
        require_matching_complete_derivations(&stable, &stable)
            .expect("identical complete derivations");
        let error = require_matching_complete_derivations(
            "metadata-and-tree-before",
            "metadata-and-tree-after",
        )
        .expect_err("moving graph must not produce a receipt");
        assert!(error.contains("two complete receipt derivations"));
        assert!(error.contains("no receipt was emitted"));
    }

    #[test]
    fn selection_refuses_ambiguous_cargo_surfaces() {
        for flag in [
            "--workspace",
            "--tests",
            "--all-targets",
            "--profile=release",
            "--release",
            "--bin=probe",
            "--lib",
        ] {
            let error = Selection::parse(&["--package".into(), "fs-roofline".into(), flag.into()])
                .expect_err(flag);
            assert!(
                error.contains("no dependency-graph claim was made"),
                "{error}"
            );
        }
    }

    #[test]
    fn selection_requires_exactly_one_explicit_root() {
        let missing = Selection::parse(&["--target=x86_64-unknown-linux-gnu".into()])
            .expect_err("missing root");
        assert!(missing.contains("exactly one explicit"));
        let duplicate =
            Selection::parse(&["-pfs-la".into(), "--package".into(), "fs-roofline".into()])
                .expect_err("multiple roots");
        assert!(duplicate.contains("exactly one explicit"));
    }

    #[test]
    fn cargo_tree_rows_preserve_distinct_unit_feature_sets() {
        let first = parse_tree_row("1semver v1.0.28|default,std").expect("first");
        let second = parse_tree_row("1semver v1.0.28|default,serde,std").expect("second");
        assert_eq!(first.name, second.name);
        assert_ne!(first.features, second.features);
    }

    #[test]
    fn metadata_parser_preserves_structured_package_and_source_ids() {
        let manifest = Path::new(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
        let text = format!(
            "{{\"packages\":[{{\"id\":\"registry+https://example.invalid/index#demo@1.2.3\",\
             \"name\":\"demo\",\"version\":\"1.2.3\",\
             \"source\":\"registry+https://example.invalid/index\",\
             \"manifest_path\":{:?}}}]}}",
            manifest.to_string_lossy()
        );
        let metadata = parse_metadata(&text).expect("metadata");
        let package = &metadata.packages[0];
        assert_eq!(
            package.metadata_id,
            "registry+https://example.invalid/index#demo@1.2.3"
        );
        assert_eq!(
            package.source_id.as_deref(),
            Some("registry+https://example.invalid/index")
        );
    }

    #[test]
    fn human_tree_identity_must_map_to_one_metadata_package() {
        let package = MetadataPackage {
            metadata_id: "registry+a#demo@1.0.0".to_string(),
            name: "demo".to_string(),
            version: "1.0.0".to_string(),
            source_id: Some("registry+a".to_string()),
            manifest_dir: PathBuf::from("/unused"),
        };
        let metadata = MetadataIndex {
            packages: vec![
                package.clone(),
                MetadataPackage {
                    metadata_id: "registry+b#demo@1.0.0".to_string(),
                    source_id: Some("registry+b".to_string()),
                    ..package
                },
            ],
            by_manifest_dir: BTreeMap::new(),
            by_name_version: BTreeMap::from([(
                ("demo".to_string(), "1.0.0".to_string()),
                vec![0, 1],
            )]),
        };
        let row = parse_tree_row("1demo v1.0.0|").expect("row");
        assert!(resolve_tree_package(&row, &metadata).is_err());
    }

    #[test]
    fn capped_collector_retains_only_its_declared_limit() {
        let read = drain_capped(std::io::Cursor::new(b"0123456789"), 4).expect("read");
        assert_eq!(read.bytes, b"0123");
        assert_eq!(read.total, 10);
        assert!(read.overflowed);
    }

    #[test]
    fn processing_bounds_refuse_before_unbounded_growth() {
        let oversized = "x".repeat(MAX_ARG_BYTES + 1);
        assert!(Selection::parse(&[oversized]).is_err());
        let root = Path::new("/workspace/package");
        let mut scheduled = MAX_PATH_PACKAGE_DIRECTORIES;
        let mut pending = Vec::new();
        assert!(
            enqueue_directory(
                root,
                &root.join("one-too-many"),
                1,
                &mut scheduled,
                &mut pending,
            )
            .expect_err("wide empty-directory tree must be bounded")
            .contains("directory processing bound")
        );
        assert!(
            pending.is_empty(),
            "refusal must happen before queue growth"
        );
        let mut scheduled = 0;
        let too_deep = root.join(
            std::iter::repeat_n("nested", MAX_PATH_PACKAGE_DEPTH + 1)
                .collect::<Vec<_>>()
                .join("/"),
        );
        assert!(
            enqueue_directory(
                root,
                &too_deep,
                MAX_PATH_PACKAGE_DEPTH + 1,
                &mut scheduled,
                &mut pending,
            )
            .expect_err("deep empty-directory chain must be bounded")
            .contains("directory-depth bound")
        );
        assert!(
            pending.is_empty(),
            "depth refusal must happen before queue growth"
        );
    }

    #[test]
    fn aggregate_path_package_admission_refuses_before_any_hashing() {
        let roots = (0..=MAX_PATH_CLOSURE_PACKAGES)
            .map(|index| PathBuf::from(format!("/individually-bounded/package-{index:04}")))
            .collect::<BTreeSet<_>>();
        let hash_calls = std::cell::Cell::new(0usize);
        let error = hash_distinct_path_packages_with(&roots, |_, _| {
            hash_calls.set(hash_calls.get() + 1);
            Ok("0".repeat(64))
        })
        .expect_err("aggregate package count must refuse before hashing");
        assert!(error.contains("aggregate hashing bound"), "{error}");
        assert_eq!(
            hash_calls.get(),
            0,
            "whole-closure preflight must precede the first package hash"
        );
    }

    #[test]
    fn aggregate_path_package_budget_is_checked_across_snapshots() {
        let root = Path::new("/workspace/package");
        let mut bytes = PathClosureBudget {
            read_bytes: MAX_PATH_CLOSURE_READ_BYTES - 1,
            entries: 0,
        };
        let error = bytes
            .charge_read(root, 2)
            .expect_err("aggregate byte overflow must refuse");
        assert!(error.contains("aggregate read bound"), "{error}");
        assert_eq!(bytes.read_bytes, MAX_PATH_CLOSURE_READ_BYTES - 1);

        let mut entries = PathClosureBudget {
            read_bytes: 0,
            entries: MAX_PATH_CLOSURE_ENTRIES,
        };
        let error = entries
            .charge_entry(root)
            .expect_err("aggregate entry overflow must refuse");
        assert!(error.contains("aggregate work bound"), "{error}");
        assert_eq!(entries.entries, MAX_PATH_CLOSURE_ENTRIES);
    }

    fn synthetic_file_identity(tag: u64) -> FileIdentity {
        FileIdentity {
            len: 4,
            readonly: false,
            #[cfg(unix)]
            device: 1,
            #[cfg(unix)]
            inode: tag,
            #[cfg(unix)]
            mode: 0o100_644,
            #[cfg(unix)]
            modified_seconds: 10,
            #[cfg(unix)]
            modified_nanoseconds: 20,
            #[cfg(unix)]
            changed_seconds: 30,
            #[cfg(unix)]
            changed_nanoseconds: i64::try_from(tag).expect("fixture tag fits i64"),
            #[cfg(not(unix))]
            modified: std::time::UNIX_EPOCH + std::time::Duration::from_secs(tag),
        }
    }

    fn synthetic_package_snapshot() -> PackageSnapshot {
        PackageSnapshot {
            files: vec![PackageFile {
                path: PathBuf::from("/workspace/package/src/lib.rs"),
                relative: "src/lib.rs".to_string(),
                symlink_target: None,
                bytes: 4,
                content_digest: "a".repeat(64),
                identity: synthetic_file_identity(7),
            }],
            directories: vec![String::new(), "src".to_string()],
            bytes: 4,
            manifest_bytes: 256,
        }
    }

    #[test]
    fn coherent_snapshot_guard_rejects_same_length_content_symlink_and_path_movement() {
        let root = Path::new("/workspace/package");
        let before = synthetic_package_snapshot();
        assert!(ensure_matching_package_snapshots(root, &before, &before).is_ok());

        let mut same_length_mutation = before.clone();
        same_length_mutation.files[0].content_digest = "b".repeat(64);
        assert!(
            ensure_matching_package_snapshots(root, &before, &same_length_mutation).is_err(),
            "equal-length content drift must not survive the second snapshot"
        );

        let mut symlink_movement = before.clone();
        symlink_movement.files[0].symlink_target =
            Some("raw=other.rs\nresolved=src/other.rs".to_string());
        assert!(ensure_matching_package_snapshots(root, &before, &symlink_movement).is_err());

        let mut path_set_movement = before.clone();
        path_set_movement.files[0].relative = "src/moved.rs".to_string();
        assert!(ensure_matching_package_snapshots(root, &before, &path_set_movement).is_err());
    }

    #[test]
    fn cargo_identity_guard_rejects_executable_or_version_movement() {
        let before = receipt_format::CargoIdentity {
            executable_digest: "a".repeat(64),
            version: "cargo 1.90.0".to_string(),
        };
        assert!(ensure_cargo_identity_stable(&before, &before).is_ok());
        let mut digest_moved = before.clone();
        digest_moved.executable_digest = "b".repeat(64);
        assert!(ensure_cargo_identity_stable(&before, &digest_moved).is_err());
        let mut version_moved = before.clone();
        version_moved.version.push_str(" drifted");
        assert!(ensure_cargo_identity_stable(&before, &version_moved).is_err());
    }

    #[test]
    fn path_receipt_excludes_only_package_root_build_and_vcs_directories() {
        let root = Path::new("/workspace/package");
        assert!(excluded_build_entry(root, &root.join("target")));
        assert!(excluded_build_entry(root, &root.join(".git")));
        assert!(points_into_excluded_root(Path::new("target/generated.rs")));
        assert!(points_into_excluded_root(Path::new(".git/config")));
        assert!(
            !excluded_build_entry(root, &root.join("src/target")),
            "a source module directory named target is a compiler input"
        );
        assert!(
            !excluded_build_entry(root, &root.join("fixtures/.git")),
            "nested data named .git remains part of the package-root bytes"
        );
        assert!(!points_into_excluded_root(Path::new("src/target/table.rs")));
    }
}
