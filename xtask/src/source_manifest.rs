//! Canonical structural source manifest for the FrankenSim trust cone.
//!
//! The tracked artifact deliberately does not serialize the Git commit that
//! will contain it: a file cannot bind its own eventual commit without an
//! impossible self-reference. Instead it binds the exact non-tracker tracked
//! Git-index source bytes, workspace and standalone crate inventory, constellation pins,
//! toolchain configuration, unsafe-capsule registry, and isolated external
//! surfaces. E13.3 owns the release envelope that must additionally bind the
//! final commit, complete tree snapshot, build host, and retained bootstrap
//! provenance receipt.

use super::depgraph::{JsonParser, JsonValue};
use super::{Violation, constellation_assessment};
use std::collections::BTreeMap;
#[cfg(test)]
use std::collections::BTreeSet;
use std::fmt::Write as _;
use std::io::{BufRead as _, BufReader, Read as _, Write as _};
use std::path::{Component, Path};
use std::process::Stdio;

pub(crate) const CHECK: &str = "source-manifest";
pub(crate) const MANIFEST_PATH: &str = "frankensim-source-manifest.json";

const SCHEMA: &str = "frankensim-source-manifest-v1";
const IDENTITY_DOMAIN: &str = "org.frankensim.xtask.source-manifest.v1";
const SOURCE_ROOT_DOMAIN: &str = "org.frankensim.xtask.source-root.v1";
const BEAD_ID: &str = "frankensim-extreal-program-f85xj.13.2";
const REPOSITORY: &str = "https://github.com/Dicklesworthstone/frankensim";
const MAX_TRACKED_FILES: usize = 10_000;
const MAX_FILE_BYTES: u64 = 64 * 1024 * 1024;
const MAX_TOTAL_BYTES: u64 = 512 * 1024 * 1024;
const EXCLUDED_TRACKED_PREFIXES: &[&str] = &[".beads/"];
const REQUIRED_NEW_SOURCE: &str = "xtask/src/source_manifest.rs";
const FIRST_GENERATION_OBJECT: &str = "untracked-during-first-generation";
const ARCHIVE_INVENTORY_OBJECT: &str = "retained-source-manifest-inventory";

#[derive(Clone, Debug, PartialEq, Eq)]
struct IndexedPath {
    path: String,
    git_mode: String,
    git_object: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct SourceFile {
    path: String,
    git_mode: String,
    bytes: u64,
    content_blake3: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Toolchain {
    channel: String,
    components: Vec<String>,
    config_blake3: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CrateRow {
    name: String,
    version: String,
    layer: String,
    manifest: String,
    workspace: &'static str,
    unsafe_capsules: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CapsuleSummary {
    schema: String,
    registry_blake3: String,
    total: usize,
    by_crate: BTreeMap<String, usize>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ExternalSurface {
    name: &'static str,
    boundary: &'static str,
    manifest: &'static str,
    manifest_blake3: String,
    lock: Option<&'static str>,
    lock_blake3: Option<String>,
    external_inputs: &'static [&'static str],
    isolation: &'static str,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ManifestModel {
    workspace_version: String,
    source_root: String,
    source_files: Vec<SourceFile>,
    lock_hash: String,
    lock_blake3: String,
    siblings: Vec<SiblingRow>,
    toolchain: Toolchain,
    crates: Vec<CrateRow>,
    capsules: CapsuleSummary,
    external_surfaces: Vec<ExternalSurface>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct SiblingRow {
    lib: String,
    version: String,
    git_head: String,
    remote: String,
    boundary: &'static str,
    runtime_consumers: Vec<String>,
    dev_consumers: Vec<String>,
    production_references: usize,
    test_references: usize,
}

fn json_string(output: &mut String, value: &str) {
    output.push('"');
    output.push_str(&super::json_escape(value));
    output.push('"');
}

fn json_strings<'a>(output: &mut String, values: impl IntoIterator<Item = &'a str>) {
    output.push('[');
    for (index, value) in values.into_iter().enumerate() {
        if index > 0 {
            output.push_str(", ");
        }
        json_string(output, value);
    }
    output.push(']');
}

fn blake3(bytes: &[u8]) -> String {
    let mut hasher = fs_blake3::Blake3::new();
    hasher.update(bytes);
    hasher.finalize().to_string()
}

fn append_identity_field(payload: &mut Vec<u8>, label: &str, value: &[u8]) {
    payload.extend_from_slice(&(label.len() as u64).to_le_bytes());
    payload.extend_from_slice(label.as_bytes());
    payload.extend_from_slice(&(value.len() as u64).to_le_bytes());
    payload.extend_from_slice(value);
}

fn validate_relative_path(path: &str) -> Result<(), String> {
    if path.is_empty()
        || path.contains('\\')
        || path.chars().any(char::is_control)
        || Path::new(path).is_absolute()
        || Path::new(path)
            .components()
            .any(|component| !matches!(component, Component::Normal(_)))
    {
        return Err(format!("source-manifest path is not canonical: {path:?}"));
    }
    Ok(())
}

fn excluded_from_structural_source(path: &str) -> bool {
    path == MANIFEST_PATH
        || EXCLUDED_TRACKED_PREFIXES
            .iter()
            .any(|prefix| path.starts_with(prefix))
}

fn git_index_paths(root: &Path) -> Result<Vec<IndexedPath>, String> {
    let output = super::constellation_cleanliness::sanitized_git_command(
        root,
        &["ls-files", "--stage", "-z"],
    )
    .output()
    .map_err(|error| format!("cannot enumerate tracked source files: {error}"))?;
    if !output.status.success() {
        return Err(format!(
            "git ls-files --stage failed while building {MANIFEST_PATH}: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        ));
    }
    let mut rows = Vec::new();
    for raw in output
        .stdout
        .split(|byte| *byte == 0)
        .filter(|row| !row.is_empty())
    {
        let row = std::str::from_utf8(raw)
            .map_err(|error| format!("git index path is not UTF-8: {error}"))?;
        let (metadata, path) = row
            .split_once('\t')
            .ok_or_else(|| format!("malformed git ls-files row: {row:?}"))?;
        let mut fields = metadata.split_whitespace();
        let mode = fields
            .next()
            .ok_or_else(|| format!("missing mode in git row: {row:?}"))?;
        let object = fields
            .next()
            .ok_or_else(|| format!("missing object in git row: {row:?}"))?;
        let stage = fields
            .next()
            .ok_or_else(|| format!("missing stage in git row: {row:?}"))?;
        if fields.next().is_some() || stage != "0" {
            return Err(format!(
                "source manifest refuses non-stage-zero or malformed index row: {row:?}"
            ));
        }
        validate_relative_path(path)?;
        if excluded_from_structural_source(path) {
            continue;
        }
        if !matches!(mode, "100644" | "100755" | "120000" | "160000") {
            return Err(format!(
                "source manifest does not support git mode {mode:?} for {path}"
            ));
        }
        rows.push(IndexedPath {
            path: path.to_string(),
            git_mode: mode.to_string(),
            git_object: object.to_string(),
        });
        if rows.len() > MAX_TRACKED_FILES {
            return Err(format!(
                "source manifest exceeds the {MAX_TRACKED_FILES}-file bound"
            ));
        }
    }
    if root.join(REQUIRED_NEW_SOURCE).is_file()
        && !rows.iter().any(|row| row.path == REQUIRED_NEW_SOURCE)
    {
        rows.push(IndexedPath {
            path: REQUIRED_NEW_SOURCE.to_string(),
            git_mode: "100644".to_string(),
            git_object: FIRST_GENERATION_OBJECT.to_string(),
        });
    }
    rows.sort_by(|left, right| left.path.cmp(&right.path));
    if rows.is_empty() {
        return Err("source manifest discovered no tracked source files".to_string());
    }
    if rows.windows(2).any(|pair| pair[0].path == pair[1].path) {
        return Err("source manifest discovered a duplicate tracked path".to_string());
    }
    Ok(rows)
}

fn json_object<'a>(
    value: &'a JsonValue,
    context: &str,
) -> Result<&'a BTreeMap<String, JsonValue>, String> {
    if let JsonValue::Object(object) = value {
        Ok(object)
    } else {
        Err(format!("{context} must be a JSON object"))
    }
}

fn json_array<'a>(value: &'a JsonValue, context: &str) -> Result<&'a [JsonValue], String> {
    if let JsonValue::Array(values) = value {
        Ok(values)
    } else {
        Err(format!("{context} must be a JSON array"))
    }
}

fn json_value_string<'a>(value: &'a JsonValue, context: &str) -> Result<&'a str, String> {
    if let JsonValue::String(value) = value {
        Ok(value)
    } else {
        Err(format!("{context} must be a JSON string"))
    }
}

fn json_value_u64(value: &JsonValue, context: &str) -> Result<u64, String> {
    if let JsonValue::Number(value) = value {
        value
            .parse()
            .map_err(|error| format!("{context} is not a canonical u64: {error}"))
    } else {
        Err(format!("{context} must be a JSON number"))
    }
}

fn required_json_field<'a>(
    object: &'a BTreeMap<String, JsonValue>,
    field: &str,
    context: &str,
) -> Result<&'a JsonValue, String> {
    object
        .get(field)
        .ok_or_else(|| format!("{context} is missing {field:?}"))
}

fn retained_index_paths_from_text(text: &str) -> Result<Vec<IndexedPath>, String> {
    let value = JsonParser::new(text)
        .finish()
        .map_err(|error| format!("retained source manifest is not valid JSON: {error}"))?;
    let root = json_object(&value, "retained source manifest")?;
    let schema = json_value_string(
        required_json_field(root, "schema", "retained source manifest")?,
        "retained source manifest schema",
    )?;
    if schema != SCHEMA {
        return Err(format!(
            "retained source manifest has schema {schema:?}, expected {SCHEMA:?}"
        ));
    }
    let frankensim = json_object(
        required_json_field(root, "frankensim", "retained source manifest")?,
        "retained source manifest frankensim section",
    )?;
    let declared_count = json_value_u64(
        required_json_field(
            frankensim,
            "tracked_file_count",
            "retained source manifest frankensim section",
        )?,
        "retained source manifest tracked_file_count",
    )?;
    let files = json_array(
        required_json_field(
            frankensim,
            "files",
            "retained source manifest frankensim section",
        )?,
        "retained source manifest files",
    )?;
    if files.len() > MAX_TRACKED_FILES {
        return Err(format!(
            "retained source manifest exceeds the {MAX_TRACKED_FILES}-file bound"
        ));
    }
    if declared_count
        != u64::try_from(files.len())
            .map_err(|_| "retained source-manifest file count does not fit u64".to_string())?
    {
        return Err(format!(
            "retained source manifest declares {declared_count} files but contains {} rows",
            files.len()
        ));
    }

    let mut total = 0_u64;
    let mut rows = Vec::with_capacity(files.len());
    for (index, value) in files.iter().enumerate() {
        let context = format!("retained source manifest file row {index}");
        let object = json_object(value, &context)?;
        if object.len() != 4
            || !["path", "git_mode", "bytes", "content_blake3"]
                .iter()
                .all(|field| object.contains_key(*field))
        {
            return Err(format!(
                "{context} must contain exactly path, git_mode, bytes, and content_blake3"
            ));
        }
        let path = json_value_string(
            required_json_field(object, "path", &context)?,
            &format!("{context} path"),
        )?;
        validate_relative_path(path)?;
        if excluded_from_structural_source(path) {
            return Err(format!("{context} contains excluded path {path:?}"));
        }
        let git_mode = json_value_string(
            required_json_field(object, "git_mode", &context)?,
            &format!("{context} git_mode"),
        )?;
        if !matches!(git_mode, "100644" | "100755" | "120000") {
            return Err(format!(
                "{context} cannot rehydrate git mode {git_mode:?} without Git metadata"
            ));
        }
        let bytes = json_value_u64(
            required_json_field(object, "bytes", &context)?,
            &format!("{context} bytes"),
        )?;
        if bytes > MAX_FILE_BYTES {
            return Err(format!(
                "{context} exceeds the {MAX_FILE_BYTES}-byte per-file bound"
            ));
        }
        total = total
            .checked_add(bytes)
            .ok_or_else(|| "retained source-manifest byte count overflow".to_string())?;
        if total > MAX_TOTAL_BYTES {
            return Err(format!(
                "retained source manifest exceeds the {MAX_TOTAL_BYTES}-byte total bound"
            ));
        }
        let digest = json_value_string(
            required_json_field(object, "content_blake3", &context)?,
            &format!("{context} content_blake3"),
        )?;
        if digest.len() != 64
            || !digest
                .bytes()
                .all(|byte| byte.is_ascii_digit() || matches!(byte, b'a'..=b'f'))
        {
            return Err(format!("{context} has a non-canonical BLAKE3 digest"));
        }
        rows.push(IndexedPath {
            path: path.to_string(),
            git_mode: git_mode.to_string(),
            git_object: ARCHIVE_INVENTORY_OBJECT.to_string(),
        });
    }
    if rows.is_empty() {
        return Err("retained source manifest contains no source files".to_string());
    }
    if rows.windows(2).any(|pair| pair[0].path >= pair[1].path) {
        return Err(
            "retained source-manifest paths must be unique and strictly sorted".to_string(),
        );
    }
    Ok(rows)
}

fn indexed_paths(root: &Path) -> Result<Vec<IndexedPath>, String> {
    match git_index_paths(root) {
        Ok(rows) => Ok(rows),
        Err(git_error) if !root.join(".git").exists() => {
            let text = std::fs::read_to_string(root.join(MANIFEST_PATH)).map_err(|error| {
                format!(
                    "{git_error}; cannot read {MANIFEST_PATH} as the archive inventory fallback: {error}"
                )
            })?;
            retained_index_paths_from_text(&text).map_err(|fallback_error| {
                format!(
                    "{git_error}; {MANIFEST_PATH} archive inventory fallback failed: {fallback_error}"
                )
            })
        }
        Err(error) => Err(error),
    }
}

fn capture_worktree_file(root: &Path, indexed: &IndexedPath) -> Result<SourceFile, String> {
    let path = root.join(&indexed.path);
    let metadata = std::fs::symlink_metadata(&path)
        .map_err(|error| format!("cannot inspect tracked source {}: {error}", indexed.path))?;
    let bytes = if indexed.git_mode == "120000" {
        if !metadata.file_type().is_symlink() {
            return Err(format!(
                "tracked symlink {} is not a symlink in the worktree",
                indexed.path
            ));
        }
        std::fs::read_link(&path)
            .map_err(|error| format!("cannot read symlink {}: {error}", indexed.path))?
            .to_str()
            .ok_or_else(|| format!("symlink target for {} is not UTF-8", indexed.path))?
            .as_bytes()
            .to_vec()
    } else {
        if !metadata.file_type().is_file() {
            return Err(format!(
                "tracked source {} is not a regular file",
                indexed.path
            ));
        }
        if metadata.len() > MAX_FILE_BYTES {
            return Err(format!(
                "tracked source {} exceeds the {MAX_FILE_BYTES}-byte per-file bound",
                indexed.path
            ));
        }
        std::fs::read(&path)
            .map_err(|error| format!("cannot read tracked source {}: {error}", indexed.path))?
    };
    let byte_count = u64::try_from(bytes.len())
        .map_err(|_| format!("byte count does not fit u64 for {}", indexed.path))?;
    Ok(SourceFile {
        path: indexed.path.clone(),
        git_mode: indexed.git_mode.clone(),
        bytes: byte_count,
        content_blake3: blake3(&bytes),
    })
}

fn is_worktree_inventory(indexed: &IndexedPath) -> bool {
    matches!(
        indexed.git_object.as_str(),
        FIRST_GENERATION_OBJECT | ARCHIVE_INVENTORY_OBJECT
    )
}

fn git_blob_queries(indexed: &[IndexedPath]) -> Vec<u8> {
    let mut queries = Vec::new();
    for row in indexed {
        if row.git_mode != "160000" && !is_worktree_inventory(row) {
            queries.extend_from_slice(row.git_object.as_bytes());
            queries.push(b'\n');
        }
    }
    queries
}

fn read_git_blob(
    reader: &mut BufReader<impl std::io::Read>,
    indexed: &IndexedPath,
) -> Result<Vec<u8>, String> {
    let mut header = String::new();
    reader.read_line(&mut header).map_err(|error| {
        format!(
            "cannot read git cat-file header for {}: {error}",
            indexed.path
        )
    })?;
    if header.len() > 256 {
        return Err(format!(
            "git cat-file returned an oversized header for {}",
            indexed.path
        ));
    }
    let mut fields = header.split_whitespace();
    let object = fields
        .next()
        .ok_or_else(|| format!("git cat-file returned no object for {}", indexed.path))?;
    let kind = fields
        .next()
        .ok_or_else(|| format!("git cat-file returned no kind for {}", indexed.path))?;
    let bytes = fields
        .next()
        .ok_or_else(|| format!("git cat-file returned no size for {}", indexed.path))?
        .parse::<u64>()
        .map_err(|error| {
            format!(
                "git cat-file returned a malformed size for {}: {error}",
                indexed.path
            )
        })?;
    if fields.next().is_some() || object != indexed.git_object || kind != "blob" {
        return Err(format!(
            "git cat-file returned an unexpected header for {}: {:?}",
            indexed.path,
            header.trim_end()
        ));
    }
    if bytes > MAX_FILE_BYTES {
        return Err(format!(
            "indexed source {} exceeds the {MAX_FILE_BYTES}-byte per-file bound",
            indexed.path
        ));
    }
    let length = usize::try_from(bytes).map_err(|_| {
        format!(
            "indexed source size does not fit usize for {}",
            indexed.path
        )
    })?;
    let mut content = vec![0_u8; length];
    reader.read_exact(&mut content).map_err(|error| {
        format!(
            "cannot read {bytes} git-index bytes for {}: {error}",
            indexed.path
        )
    })?;
    let mut delimiter = [0_u8; 1];
    reader.read_exact(&mut delimiter).map_err(|error| {
        format!(
            "cannot read git cat-file delimiter for {}: {error}",
            indexed.path
        )
    })?;
    if delimiter != [b'\n'] {
        return Err(format!(
            "git cat-file returned a malformed blob delimiter for {}",
            indexed.path
        ));
    }
    Ok(content)
}

fn capture_indexed_source(root: &Path, indexed: &[IndexedPath]) -> Result<Vec<SourceFile>, String> {
    let queries = git_blob_queries(indexed);
    let mut child = if queries.is_empty() {
        None
    } else {
        Some(
            super::constellation_cleanliness::sanitized_git_command(root, &["cat-file", "--batch"])
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|error| format!("cannot start git cat-file --batch: {error}"))?,
        )
    };
    let mut writer = if let Some(child) = child.as_mut() {
        let mut stdin = child
            .stdin
            .take()
            .ok_or_else(|| "git cat-file stdin was not piped".to_string())?;
        Some(std::thread::spawn(move || {
            stdin
                .write_all(&queries)
                .map_err(|error| format!("cannot query git-index blobs: {error}"))
        }))
    } else {
        None
    };
    let mut stdout = child
        .as_mut()
        .and_then(|child| child.stdout.take())
        .map(BufReader::new);

    let captured = (|| {
        let mut total = 0_u64;
        let mut rows = Vec::with_capacity(indexed.len());
        for row in indexed {
            let captured = if row.git_mode == "160000" {
                SourceFile {
                    path: row.path.clone(),
                    git_mode: row.git_mode.clone(),
                    bytes: u64::try_from(row.git_object.len()).map_err(|_| {
                        format!("gitlink object size does not fit u64 for {}", row.path)
                    })?,
                    content_blake3: blake3(row.git_object.as_bytes()),
                }
            } else if is_worktree_inventory(row) {
                capture_worktree_file(root, row)?
            } else {
                let content = read_git_blob(
                    stdout
                        .as_mut()
                        .ok_or_else(|| "git cat-file stdout was not piped".to_string())?,
                    row,
                )?;
                SourceFile {
                    path: row.path.clone(),
                    git_mode: row.git_mode.clone(),
                    bytes: u64::try_from(content.len())
                        .map_err(|_| format!("byte count does not fit u64 for {}", row.path))?,
                    content_blake3: blake3(&content),
                }
            };
            total = total
                .checked_add(captured.bytes)
                .ok_or_else(|| "source-manifest byte count overflow".to_string())?;
            if total > MAX_TOTAL_BYTES {
                return Err(format!(
                    "source manifest exceeds the {MAX_TOTAL_BYTES}-byte total bound"
                ));
            }
            rows.push(captured);
        }
        Ok(rows)
    })();

    if let Err(error) = captured {
        if let Some(child) = child.as_mut() {
            let _ = child.kill();
        }
        drop(stdout);
        if let Some(writer) = writer.take() {
            let _ = writer.join();
        }
        if let Some(mut child) = child {
            let _ = child.wait();
        }
        return Err(error);
    }
    if let Some(writer) = writer.take() {
        match writer.join() {
            Ok(Ok(())) => {}
            Ok(Err(error)) => {
                if let Some(child) = child.as_mut() {
                    let _ = child.kill();
                    let _ = child.wait();
                }
                return Err(error);
            }
            Err(_) => {
                if let Some(child) = child.as_mut() {
                    let _ = child.kill();
                    let _ = child.wait();
                }
                return Err("git cat-file query writer panicked".to_string());
            }
        }
    }
    drop(stdout);
    if let Some(mut child) = child {
        let status = child
            .wait()
            .map_err(|error| format!("cannot wait for git cat-file --batch: {error}"))?;
        if !status.success() {
            let mut stderr = String::new();
            if let Some(mut stream) = child.stderr.take() {
                stream
                    .read_to_string(&mut stderr)
                    .map_err(|error| format!("cannot read git cat-file stderr: {error}"))?;
            }
            return Err(format!(
                "git cat-file --batch failed while capturing the source index: {}",
                stderr.trim()
            ));
        }
    }
    captured
}

fn capture_source_files_once(root: &Path) -> Result<Vec<SourceFile>, String> {
    let indexed = indexed_paths(root)?;
    capture_indexed_source(root, &indexed)
}

fn capture_source_files(root: &Path) -> Result<Vec<SourceFile>, String> {
    let first = capture_source_files_once(root)?;
    let second = capture_source_files_once(root)?;
    if first != second {
        return Err(
            "tracked source changed between two complete source-manifest captures".to_string(),
        );
    }
    Ok(first)
}

fn source_root(files: &[SourceFile]) -> String {
    let mut payload = Vec::new();
    append_identity_field(&mut payload, "schema", b"frankensim-source-root-v1");
    for file in files {
        append_identity_field(&mut payload, "path", file.path.as_bytes());
        append_identity_field(&mut payload, "git-mode", file.git_mode.as_bytes());
        append_identity_field(&mut payload, "bytes", &file.bytes.to_le_bytes());
        append_identity_field(
            &mut payload,
            "content-blake3",
            file.content_blake3.as_bytes(),
        );
    }
    fs_blake3::hash_domain(SOURCE_ROOT_DOMAIN, &payload).to_string()
}

fn package_value(text: &str, section_name: &str, key_name: &str) -> Result<String, String> {
    let mut section = "";
    for raw in text.lines() {
        let line = raw.trim();
        if line.starts_with('[') {
            section = line;
            continue;
        }
        if section != section_name {
            continue;
        }
        let Some((key, value)) = line.split_once('=') else {
            continue;
        };
        if key.trim() == key_name {
            return super::casual_manifest_string(value.trim()).map_err(|error| {
                format!("cannot parse {section_name}.{key_name} from Cargo.toml: {error}")
            });
        }
    }
    Err(format!("Cargo.toml is missing {section_name}.{key_name}"))
}

fn package_version(text: &str) -> Result<String, String> {
    package_value(text, "[package]", "version")
}

fn workspace_version(text: &str) -> Result<String, String> {
    package_value(text, "[workspace.package]", "version")
}

fn toolchain(root: &Path) -> Result<Toolchain, String> {
    let path = root.join("rust-toolchain.toml");
    let text = std::fs::read_to_string(&path)
        .map_err(|error| format!("cannot read rust-toolchain.toml: {error}"))?;
    let mut section = "";
    let mut channel = None;
    let mut components = None;
    for raw in text.lines() {
        let line = raw.trim();
        if line.starts_with('[') {
            section = line;
            continue;
        }
        if section != "[toolchain]" {
            continue;
        }
        let Some((key, value)) = line.split_once('=') else {
            continue;
        };
        match key.trim() {
            "channel" => {
                channel = Some(
                    super::casual_manifest_string(value.trim())
                        .map_err(|error| format!("invalid toolchain channel: {error}"))?,
                );
            }
            "components" => {
                components = Some(
                    super::casual_manifest_string_array(value.trim())
                        .map_err(|error| format!("invalid toolchain components: {error}"))?,
                );
            }
            _ => {}
        }
    }
    let mut components =
        components.ok_or_else(|| "rust-toolchain.toml has no components".to_string())?;
    components.sort();
    components.dedup();
    if components.is_empty() {
        return Err("rust-toolchain.toml has an empty component set".to_string());
    }
    Ok(Toolchain {
        channel: channel.ok_or_else(|| "rust-toolchain.toml has no channel".to_string())?,
        components,
        config_blake3: blake3(text.as_bytes()),
    })
}

fn json_line_string(line: &str, field: &str) -> Option<String> {
    let marker = format!("\"{field}\": \"");
    let start = line.find(&marker)? + marker.len();
    let rest = &line[start..];
    let end = rest.find('"')?;
    Some(rest[..end].to_string())
}

fn capsule_summary(root: &Path) -> Result<CapsuleSummary, String> {
    let text = std::fs::read_to_string(root.join("unsafe-capsules.json"))
        .map_err(|error| format!("cannot read unsafe-capsules.json: {error}"))?;
    let schema = text
        .lines()
        .find_map(|line| json_line_string(line, "schema"))
        .ok_or_else(|| "unsafe-capsules.json has no schema".to_string())?;
    let mut by_crate = BTreeMap::new();
    let mut total = 0_usize;
    for line in text.lines().map(str::trim) {
        if !line.starts_with("{\"crate\"") {
            continue;
        }
        let crate_name = json_line_string(line, "crate")
            .ok_or_else(|| "unsafe capsule row has no crate".to_string())?;
        *by_crate.entry(crate_name).or_insert(0) += 1;
        total = total
            .checked_add(1)
            .ok_or_else(|| "unsafe capsule count overflow".to_string())?;
    }
    if total != super::registry_modules(&text).len() || total == 0 {
        return Err("unsafe-capsule summary disagrees with the policy registry parser".to_string());
    }
    Ok(CapsuleSummary {
        schema,
        registry_blake3: blake3(text.as_bytes()),
        total,
        by_crate,
    })
}

fn crate_inventory(
    root: &Path,
    workspace_version: &str,
    capsules: &CapsuleSummary,
) -> Result<Vec<CrateRow>, String> {
    let root_manifest = std::fs::read_to_string(root.join("Cargo.toml"))
        .map_err(|error| format!("cannot read root Cargo.toml: {error}"))?;
    let mut rows = Vec::new();
    for manifest in super::load_workspace(root)? {
        let relative = manifest
            .dir
            .join("Cargo.toml")
            .strip_prefix(root)
            .map_err(|_| format!("crate {} escaped the workspace root", manifest.name))?
            .to_string_lossy()
            .replace('\\', "/");
        let member_marker = format!("\"crates/{}\"", manifest.name);
        let workspace = if root_manifest.contains(&member_marker) {
            "native"
        } else {
            "standalone"
        };
        let version = if workspace == "native" {
            workspace_version.to_string()
        } else {
            let text = std::fs::read_to_string(root.join(&relative))
                .map_err(|error| format!("cannot read {relative}: {error}"))?;
            package_version(&text)?
        };
        rows.push(CrateRow {
            unsafe_capsules: capsules.by_crate.get(&manifest.name).copied().unwrap_or(0),
            name: manifest.name,
            version,
            layer: manifest.layer.name().to_string(),
            manifest: relative,
            workspace,
        });
    }
    rows.sort_by(|left, right| left.name.cmp(&right.name));
    Ok(rows)
}

fn external_surface(
    root: &Path,
    name: &'static str,
    boundary: &'static str,
    manifest: &'static str,
    lock: Option<&'static str>,
    external_inputs: &'static [&'static str],
    isolation: &'static str,
) -> Result<ExternalSurface, String> {
    let manifest_bytes = std::fs::read(root.join(manifest))
        .map_err(|error| format!("cannot read {manifest}: {error}"))?;
    let lock_blake3 = lock
        .map(|path| {
            std::fs::read(root.join(path))
                .map(|bytes| blake3(&bytes))
                .map_err(|error| format!("cannot read {path}: {error}"))
        })
        .transpose()?;
    Ok(ExternalSurface {
        name,
        boundary,
        manifest,
        manifest_blake3: blake3(&manifest_bytes),
        lock,
        lock_blake3,
        external_inputs,
        isolation,
    })
}

fn sibling_rows(root: &Path) -> Result<(String, String, Vec<SiblingRow>), String> {
    let lock_text = super::read_constellation_lock(&root.join("constellation.lock"))?;
    let lock_blake3 = blake3(lock_text.as_bytes());
    let (lock_hash, rows) = super::parse_lock_rows(&lock_text)?;
    let usage = constellation_assessment::source_manifest_usage(root)?
        .into_iter()
        .map(|row| (row.lib.clone(), row))
        .collect::<BTreeMap<_, _>>();
    let mut siblings = Vec::with_capacity(rows.len());
    for row in rows {
        let measured = usage
            .get(&row.lib)
            .ok_or_else(|| format!("no measured usage row for {}", row.lib))?;
        siblings.push(SiblingRow {
            lib: row.lib,
            version: row.version,
            git_head: row.git_head,
            remote: row.remote,
            boundary: measured.boundary,
            runtime_consumers: measured.runtime_consumers.clone(),
            dev_consumers: measured.dev_consumers.clone(),
            production_references: measured.production_references,
            test_references: measured.test_references,
        });
    }
    Ok((lock_hash, lock_blake3, siblings))
}

fn build_model(root: &Path) -> Result<ManifestModel, String> {
    let root_manifest = std::fs::read_to_string(root.join("Cargo.toml"))
        .map_err(|error| format!("cannot read root Cargo.toml: {error}"))?;
    if !root_manifest.contains(&format!("repository = \"{REPOSITORY}\"")) {
        return Err("workspace repository authority moved or is missing".to_string());
    }
    let workspace_version = workspace_version(&root_manifest)?;
    let source_files = capture_source_files(root)?;
    let source_root = source_root(&source_files);
    let capsules = capsule_summary(root)?;
    let crates = crate_inventory(root, &workspace_version, &capsules)?;
    let (lock_hash, lock_blake3, siblings) = sibling_rows(root)?;
    let external_surfaces = vec![
        external_surface(
            root,
            "constellation-bootstrap",
            "bootstrap-tool",
            "tools/bootstrap/Cargo.toml",
            Some("tools/bootstrap/Cargo.lock"),
            &[],
            "standalone zero-dependency pre-workspace bootstrap; not a runtime library",
        )?,
        external_surface(
            root,
            "certified-arithmetic-n-version-kernel",
            "development-reference",
            "tools/cert-kernel/Cargo.toml",
            Some("tools/cert-kernel/Cargo.lock"),
            &[],
            "nested zero-external-dependency checker; its optional fs-ivl edge is comparison-only and cannot promote production certificate authority",
        )?,
        external_surface(
            root,
            "high-precision-oracle",
            "development-reference",
            "tools/oracle/Cargo.toml",
            Some("tools/oracle/Cargo.lock"),
            &["rug = 1.30.0", "MPFR/GMP transitively through rug"],
            "nested workspace excluded from the native production graph",
        )?,
        external_surface(
            root,
            "fs-wasm",
            "standalone-distribution",
            "crates/fs-wasm/Cargo.toml",
            Some("crates/fs-wasm/Cargo.lock"),
            &["wasm-bindgen = 0.2", "getrandom = 0.4 with wasm_js"],
            "standalone nested workspace; browser-only dependencies are outside the native production graph",
        )?,
    ];
    Ok(ManifestModel {
        workspace_version,
        source_root,
        source_files,
        lock_hash,
        lock_blake3,
        siblings,
        toolchain: toolchain(root)?,
        crates,
        capsules,
        external_surfaces,
    })
}

fn render_body(model: &ManifestModel) -> String {
    let mut output = String::new();
    let _ = write!(output, "  \"bead_id\": ");
    json_string(&mut output, BEAD_ID);
    output.push_str(",\n  \"authority\": \"canonical structural inventory over exact Git-index non-tracker bytes; release observations are intentionally unbound here\",\n");
    output.push_str("  \"release_binding\": {\n");
    output.push_str("    \"status\": \"required-not-yet-attached\",\n");
    output.push_str("    \"owner\": \"frankensim-extreal-program-f85xj.13.3\",\n");
    output.push_str("    \"required_fields\": [\"frankensim_git_commit\", \"complete_tree_snapshot\", \"rustc_host\", \"retained_bootstrap_provenance_identity\"],\n");
    output.push_str("    \"no_claim\": \"This tracked artifact is not a release SBOM, self-contained bundle, build attestation, or proof that any sibling is correct.\"\n");
    output.push_str("  },\n");
    output.push_str("  \"frankensim\": {\n    \"repository\": ");
    json_string(&mut output, REPOSITORY);
    output.push_str(",\n    \"workspace_version\": ");
    json_string(&mut output, &model.workspace_version);
    output.push_str(",\n    \"source_root_domain\": ");
    json_string(&mut output, SOURCE_ROOT_DOMAIN);
    output.push_str(",\n    \"source_root\": ");
    json_string(&mut output, &model.source_root);
    let _ = write!(
        output,
        ",\n    \"tracked_file_count\": {},\n    \"excluded_paths\": [\".beads/**\", \"{}\"],\n    \"files\": [\n",
        model.source_files.len(),
        MANIFEST_PATH
    );
    for (index, file) in model.source_files.iter().enumerate() {
        output.push_str("      {\"path\": ");
        json_string(&mut output, &file.path);
        output.push_str(", \"git_mode\": ");
        json_string(&mut output, &file.git_mode);
        let _ = write!(output, ", \"bytes\": {}, \"content_blake3\": ", file.bytes);
        json_string(&mut output, &file.content_blake3);
        output.push('}');
        if index + 1 != model.source_files.len() {
            output.push(',');
        }
        output.push('\n');
    }
    output.push_str("    ]\n  },\n");
    output.push_str("  \"constellation\": {\n    \"lock_schema\": \"frankensim-constellation-lock-v2\",\n    \"lock_identity\": ");
    json_string(&mut output, &model.lock_hash);
    output.push_str(",\n    \"lock_blake3\": ");
    json_string(&mut output, &model.lock_blake3);
    output.push_str(",\n    \"bootstrap_provenance_requirement\": {\"schema\": \"frankensim-constellation-bootstrap-v2\", \"identity_domain\": \"org.frankensim.xtask.constellation-bootstrap-provenance.v3\", \"identity_version\": 3, \"retention\": \"release-envelope-input-not-retained-in-this-structural-manifest\"},\n");
    output.push_str("    \"siblings\": [\n");
    for (index, sibling) in model.siblings.iter().enumerate() {
        output.push_str("      {\"lib\": ");
        json_string(&mut output, &sibling.lib);
        output.push_str(", \"version\": ");
        json_string(&mut output, &sibling.version);
        output.push_str(", \"git_head\": ");
        json_string(&mut output, &sibling.git_head);
        output.push_str(", \"remote\": ");
        json_string(&mut output, &sibling.remote);
        output.push_str(", \"boundary\": ");
        json_string(&mut output, sibling.boundary);
        output.push_str(", \"runtime_consumers\": ");
        json_strings(
            &mut output,
            sibling.runtime_consumers.iter().map(String::as_str),
        );
        output.push_str(", \"dev_consumers\": ");
        json_strings(
            &mut output,
            sibling.dev_consumers.iter().map(String::as_str),
        );
        let _ = write!(
            output,
            ", \"production_api_references\": {}, \"test_api_references\": {}, \"verification_refs\": [",
            sibling.production_references, sibling.test_references
        );
        for (reference_index, reference) in [
            format!("constellation.lock#libraries/{}", sibling.lib),
            format!("constellation-bootstrap.json#libraries/{}", sibling.lib),
            format!(
                "constellation-trust-assessment.json#siblings/{}",
                sibling.lib
            ),
        ]
        .iter()
        .enumerate()
        {
            if reference_index > 0 {
                output.push_str(", ");
            }
            json_string(&mut output, reference);
        }
        output.push_str("]}");
        if index + 1 != model.siblings.len() {
            output.push(',');
        }
        output.push('\n');
    }
    output.push_str("    ]\n  },\n");
    output.push_str("  \"toolchain\": {\n    \"channel\": ");
    json_string(&mut output, &model.toolchain.channel);
    output.push_str(",\n    \"components\": ");
    json_strings(
        &mut output,
        model.toolchain.components.iter().map(String::as_str),
    );
    output.push_str(",\n    \"config_blake3\": ");
    json_string(&mut output, &model.toolchain.config_blake3);
    output.push_str(",\n    \"release_host\": null,\n    \"release_host_status\": \"required-in-e13.3-release-envelope\"\n  },\n");
    let native = model
        .crates
        .iter()
        .filter(|row| row.workspace == "native")
        .count();
    let standalone = model.crates.len() - native;
    let _ = write!(
        output,
        "  \"workspace\": {{\n    \"native_fs_crates\": {native},\n    \"standalone_fs_crates\": {standalone},\n    \"crates\": [\n"
    );
    for (index, row) in model.crates.iter().enumerate() {
        output.push_str("      {\"name\": ");
        json_string(&mut output, &row.name);
        output.push_str(", \"version\": ");
        json_string(&mut output, &row.version);
        output.push_str(", \"layer\": ");
        json_string(&mut output, &row.layer);
        output.push_str(", \"manifest\": ");
        json_string(&mut output, &row.manifest);
        output.push_str(", \"workspace\": ");
        json_string(&mut output, row.workspace);
        let _ = write!(output, ", \"unsafe_capsules\": {}}}", row.unsafe_capsules);
        if index + 1 != model.crates.len() {
            output.push(',');
        }
        output.push('\n');
    }
    output.push_str("    ]\n  },\n");
    output.push_str("  \"unsafe_capsules\": {\n    \"schema\": ");
    json_string(&mut output, &model.capsules.schema);
    output.push_str(",\n    \"registry_blake3\": ");
    json_string(&mut output, &model.capsules.registry_blake3);
    let _ = write!(
        output,
        ",\n    \"total\": {},\n    \"by_crate\": [",
        model.capsules.total
    );
    for (index, (crate_name, count)) in model.capsules.by_crate.iter().enumerate() {
        if index > 0 {
            output.push_str(", ");
        }
        output.push_str("{\"crate\": ");
        json_string(&mut output, crate_name);
        let _ = write!(output, ", \"count\": {count}}}");
    }
    output.push_str("]\n  },\n");
    output.push_str("  \"external_boundaries\": [\n");
    for (index, surface) in model.external_surfaces.iter().enumerate() {
        output.push_str("    {\"name\": ");
        json_string(&mut output, surface.name);
        output.push_str(", \"boundary\": ");
        json_string(&mut output, surface.boundary);
        output.push_str(", \"manifest\": ");
        json_string(&mut output, surface.manifest);
        output.push_str(", \"manifest_blake3\": ");
        json_string(&mut output, &surface.manifest_blake3);
        output.push_str(", \"lock\": ");
        if let Some(lock) = surface.lock {
            json_string(&mut output, lock);
        } else {
            output.push_str("null");
        }
        output.push_str(", \"lock_blake3\": ");
        if let Some(hash) = &surface.lock_blake3 {
            json_string(&mut output, hash);
        } else {
            output.push_str("null");
        }
        output.push_str(", \"external_inputs\": ");
        json_strings(&mut output, surface.external_inputs.iter().copied());
        output.push_str(", \"isolation\": ");
        json_string(&mut output, surface.isolation);
        output.push('}');
        if index + 1 != model.external_surfaces.len() {
            output.push(',');
        }
        output.push('\n');
    }
    output.push_str("  ],\n");
    output.push_str("  \"package_citation\": {\"status\": \"not-yet-wired\", \"current_authority\": \"EvidencePackage v8 root-binds code_version and constellation_lock only\", \"required_follow_on\": \"a versioned fs-package migration must bind this manifest identity before package-citation authority is claimed\"},\n");
    output.push_str("  \"standard_renderings\": {\"spdx\": \"staged-follow-on-after-the-in-house-content-schema-stabilizes\"}\n");
    output
}

fn render(model: &ManifestModel) -> String {
    let body = render_body(model);
    let identity = fs_blake3::hash_domain(IDENTITY_DOMAIN, body.as_bytes()).to_string();
    let mut output = String::new();
    output.push_str("{\n  \"schema\": ");
    json_string(&mut output, SCHEMA);
    output.push_str(",\n  \"identity_domain\": ");
    json_string(&mut output, IDENTITY_DOMAIN);
    output.push_str(",\n  \"identity_version\": 1,\n  \"manifest_identity\": ");
    json_string(&mut output, &identity);
    output.push_str(",\n");
    output.push_str(&body);
    output.push_str("}\n");
    output
}

fn expected_artifact(root: &Path) -> Result<String, String> {
    build_model(root).map(|model| render(&model))
}

pub(crate) fn generate(root: &Path) -> Result<(), String> {
    let manifest = expected_artifact(root)?;
    std::fs::write(root.join(MANIFEST_PATH), manifest)
        .map_err(|error| format!("cannot write {MANIFEST_PATH}: {error}"))
}

pub(crate) fn check(root: &Path) -> Vec<Violation> {
    let expected = match expected_artifact(root) {
        Ok(expected) => expected,
        Err(detail) => {
            return vec![Violation {
                check: CHECK,
                crate_name: "<repo>".to_string(),
                detail,
            }];
        }
    };
    match std::fs::read_to_string(root.join(MANIFEST_PATH)) {
        Ok(actual) if actual == expected => Vec::new(),
        Ok(_) => vec![Violation {
            check: CHECK,
            crate_name: MANIFEST_PATH.to_string(),
            detail: format!(
                "tracked source manifest is stale; run cargo run -p xtask -- generate-source-manifest"
            ),
        }],
        Err(error) => vec![Violation {
            check: CHECK,
            crate_name: MANIFEST_PATH.to_string(),
            detail: format!("cannot read retained source manifest: {error}"),
        }],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn source_file(path: &str, digest: &str) -> SourceFile {
        SourceFile {
            path: path.to_string(),
            git_mode: "100644".to_string(),
            bytes: 3,
            content_blake3: digest.to_string(),
        }
    }

    #[test]
    fn source_root_moves_with_path_mode_size_and_content() {
        let baseline = vec![source_file("crates/fs-a/src/lib.rs", "aaa")];
        let baseline_root = source_root(&baseline);
        for changed in [
            vec![source_file("crates/fs-b/src/lib.rs", "aaa")],
            vec![SourceFile {
                git_mode: "100755".to_string(),
                ..baseline[0].clone()
            }],
            vec![SourceFile {
                bytes: 4,
                ..baseline[0].clone()
            }],
            vec![source_file("crates/fs-a/src/lib.rs", "bbb")],
        ] {
            assert_ne!(baseline_root, source_root(&changed));
        }
    }

    #[test]
    fn tracker_and_manifest_paths_are_explicitly_excluded() {
        assert!(excluded_from_structural_source(".beads/issues.jsonl"));
        assert!(excluded_from_structural_source(MANIFEST_PATH));
        assert!(!excluded_from_structural_source("Cargo.toml"));
    }

    #[test]
    fn retained_inventory_is_a_bounded_archive_fallback() {
        let text = r#"{
          "schema": "frankensim-source-manifest-v1",
          "frankensim": {
            "tracked_file_count": 2,
            "files": [
              {"path": "Cargo.toml", "git_mode": "100644", "bytes": 3, "content_blake3": "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"},
              {"path": "xtask/src/main.rs", "git_mode": "100644", "bytes": 4, "content_blake3": "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"}
            ]
          }
        }"#;
        let rows = retained_index_paths_from_text(text).expect("retained inventory parses");
        assert_eq!(
            rows.iter().map(|row| row.path.as_str()).collect::<Vec<_>>(),
            ["Cargo.toml", "xtask/src/main.rs"]
        );

        let unsorted = text.replace(
            "\"Cargo.toml\", \"git_mode\": \"100644\", \"bytes\": 3",
            "\"z.toml\", \"git_mode\": \"100644\", \"bytes\": 3",
        );
        assert!(
            retained_index_paths_from_text(&unsorted)
                .expect_err("unsorted retained inventory is rejected")
                .contains("strictly sorted")
        );
    }

    #[test]
    fn git_blob_batch_uses_only_index_blob_objects() {
        let rows = [
            IndexedPath {
                path: "Cargo.toml".to_string(),
                git_mode: "100644".to_string(),
                git_object: "a".repeat(40),
            },
            IndexedPath {
                path: "submodule".to_string(),
                git_mode: "160000".to_string(),
                git_object: "b".repeat(40),
            },
            IndexedPath {
                path: REQUIRED_NEW_SOURCE.to_string(),
                git_mode: "100644".to_string(),
                git_object: FIRST_GENERATION_OBJECT.to_string(),
            },
            IndexedPath {
                path: "README.md".to_string(),
                git_mode: "100644".to_string(),
                git_object: ARCHIVE_INVENTORY_OBJECT.to_string(),
            },
        ];
        assert_eq!(
            git_blob_queries(&rows),
            format!("{}\n", "a".repeat(40)).into_bytes()
        );
    }

    #[test]
    fn live_model_retains_all_pins_and_separates_release_obligations() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("xtask has a workspace parent");
        let model = build_model(root).expect("live source manifest builds");
        assert_eq!(model.siblings.len(), 7);
        assert_eq!(
            model
                .siblings
                .iter()
                .filter(|row| row.boundary == "pinned-unused-planned")
                .map(|row| row.lib.as_str())
                .collect::<BTreeSet<_>>(),
            BTreeSet::from(["frankenpandas"])
        );
        assert_eq!(
            model
                .crates
                .iter()
                .filter(|row| row.workspace == "standalone")
                .map(|row| row.name.as_str())
                .collect::<BTreeSet<_>>(),
            BTreeSet::from(["fs-wasm"])
        );
        assert_eq!(
            model
                .external_surfaces
                .iter()
                .map(|row| row.name)
                .collect::<BTreeSet<_>>(),
            BTreeSet::from([
                "certified-arithmetic-n-version-kernel",
                "constellation-bootstrap",
                "fs-wasm",
                "high-precision-oracle",
            ])
        );
        let rendered = render(&model);
        assert!(rendered.contains("\"frankensim_git_commit\""));
        assert!(rendered.contains("\"release_host\": null"));
        assert!(rendered.contains("\"package_citation\": {\"status\": \"not-yet-wired\""));
        assert_eq!(rendered, render(&model));
    }
}
