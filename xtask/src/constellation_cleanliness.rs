//! Fail-closed cleanliness inspection shared by both constellation bootstraps.
//!
//! Git's ordinary status can inherit `submodule.*.ignore`, `.gitmodules`
//! `ignore`, repository-local excludes, and hidden index flags.  A pinned
//! constellation checkout must not let any of those local policies conceal
//! initialized nested source, so this module walks stage-0 gitlinks itself and
//! applies the same checks at every materialized repository boundary.

use std::collections::BTreeMap;
use std::fs::Metadata;
use std::io::{Read as _, Write as _};
use std::path::{Component, Path, PathBuf};
use std::process::{Command, Stdio};

const MAX_NESTED_REPOSITORY_DEPTH: usize = 32;
const MAX_RAW_HASH_PATHS_PER_BATCH: usize = 128;
const MAX_RAW_HASH_ARGUMENT_BYTES: usize = 64 * 1024;
const MAX_RAW_INDEX_BYTES: u64 = 1024 * 1024 * 1024;
const MAX_GRAFTS_BYTES: u64 = 1024 * 1024;
#[cfg(unix)]
const CORE_FILE_MODE_OBSERVATION: &str = "core.fileMode=true";
#[cfg(not(unix))]
const CORE_FILE_MODE_OBSERVATION: &str = "core.fileMode=false";
const TRACKED_STATUS_ARGS: &[&str] = &[
    "-c",
    CORE_FILE_MODE_OBSERVATION,
    "status",
    "--porcelain=v1",
    "-z",
    "--untracked-files=no",
    "--ignore-submodules=none",
    "--no-renames",
];
const STAGED_STATUS_ARGS: &[&str] = &[
    "diff",
    "--cached",
    "--name-only",
    "-z",
    "--no-renames",
    "--no-ext-diff",
    "--ignore-submodules=none",
    "--",
];

#[derive(Debug, Clone, PartialEq, Eq)]
struct StageZeroEntry {
    mode: String,
    expected_oid: String,
    path: PathBuf,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Gitlink {
    expected_head: String,
    path: PathBuf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WorktreeKind {
    Missing,
    Regular,
    Symlink,
    Directory,
    Other,
    UnsupportedIndexMode,
}

impl WorktreeKind {
    const fn name(self) -> &'static str {
        match self {
            Self::Missing => "missing",
            Self::Regular => "regular file",
            Self::Symlink => "symbolic link",
            Self::Directory => "directory",
            Self::Other => "special file",
            Self::UnsupportedIndexMode => "unsupported index mode",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RawTrackedState {
    expected_mode: String,
    expected_oid: String,
    path: PathBuf,
    actual_kind: WorktreeKind,
    actual_mode: Option<String>,
    actual_oid: Option<String>,
    link_target: Option<Vec<u8>>,
    metadata_identity: Option<MetadataIdentity>,
}

#[cfg(unix)]
#[derive(Debug, Clone, PartialEq, Eq)]
struct MetadataIdentity {
    device: u64,
    inode: u64,
    mode: u32,
    size: u64,
    modified_seconds: i64,
    modified_nanoseconds: i64,
    changed_seconds: i64,
    changed_nanoseconds: i64,
}

#[cfg(not(unix))]
#[derive(Debug, Clone, PartialEq, Eq)]
struct MetadataIdentity {
    kind: WorktreeKind,
    size: u64,
    readonly: bool,
}

#[cfg(unix)]
type RawAuthorityIdentity = MetadataIdentity;

#[cfg(windows)]
#[derive(Debug, Clone, PartialEq, Eq)]
struct RawAuthorityIdentity {
    volume_serial_number: u32,
    file_index: u64,
    attributes: u32,
    links: u32,
    size: u64,
    created: u64,
    modified: u64,
}

#[cfg(not(any(unix, windows)))]
#[derive(Debug, Clone, PartialEq, Eq)]
struct RawAuthorityIdentity;

#[derive(Debug, Clone, PartialEq, Eq)]
struct GitlinkState {
    expected_head: String,
    path: PathBuf,
    child: Option<Box<RepositorySnapshot>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct IndexMaterializationState {
    logical_path: PathBuf,
    kind: WorktreeKind,
    metadata_identity: MetadataIdentity,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RawIndexInspection {
    entry_count: usize,
    fsmonitor_extension_present: bool,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct RepositoryBoundaryObservation {
    head: Option<String>,
    index: Vec<u8>,
    raw_index: Option<Vec<u8>>,
    entries: Vec<StageZeroEntry>,
    index_materialization: Vec<IndexMaterializationState>,
    staged_status: Vec<u8>,
    forced_status: Vec<u8>,
    untracked: Vec<u8>,
    untracked_gitignores: Vec<u8>,
    assume_and_skip_flags: Vec<u8>,
    fsmonitor_extension_present: bool,
    replacement_refs: Vec<u8>,
    grafts: Option<Vec<u8>>,
    raw_tracked: Vec<RawTrackedState>,
    findings: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RepositorySnapshot {
    repository: PathBuf,
    scope: PathBuf,
    expected_head: Option<String>,
    head: Option<String>,
    index: Vec<u8>,
    raw_index: Option<Vec<u8>>,
    tracked_status: Vec<u8>,
    untracked: Vec<u8>,
    untracked_gitignores: Vec<u8>,
    assume_and_skip_flags: Vec<u8>,
    fsmonitor_extension_present: bool,
    replacement_refs: Vec<u8>,
    grafts: Option<Vec<u8>>,
    index_materialization: Vec<IndexMaterializationState>,
    raw_tracked: Vec<RawTrackedState>,
    gitlinks: Vec<GitlinkState>,
    findings: Vec<String>,
}

impl RepositorySnapshot {
    fn collect_findings(&self, findings: &mut Vec<String>) {
        findings.extend(self.findings.iter().cloned());
        for gitlink in &self.gitlinks {
            if let Some(child) = &gitlink.child {
                child.collect_findings(findings);
            }
        }
    }
}

pub(crate) fn repository_worktree_status(root: &Path) -> Result<String, String> {
    repository_worktree_status_with_expected_head(root, None)
}

pub(crate) fn pinned_repository_worktree_status(
    root: &Path,
    expected_head: &str,
) -> Result<String, String> {
    repository_worktree_status_with_expected_head(root, Some(expected_head))
}

fn repository_worktree_status_with_expected_head(
    root: &Path,
    expected_head: Option<&str>,
) -> Result<String, String> {
    if expected_head.is_some_and(|head| !valid_object_id(head)) {
        return Err("expected repository HEAD is not a canonical Git object ID".to_string());
    }
    let root_identity = ordinary_root_identity(root)?;
    let canonical_root = root
        .canonicalize()
        .map_err(|error| format!("cannot canonicalize {}: {error}", escaped_path(root)))?;
    let first = observe_repository(&canonical_root, Path::new("."), expected_head, 0)?;
    let second = observe_repository(&canonical_root, Path::new("."), expected_head, 0)?;
    require_equal_observations(&canonical_root, &first, &second)?;
    require_unchanged_root(root, &canonical_root, &root_identity)?;
    let mut findings = Vec::new();
    first.collect_findings(&mut findings);
    Ok(findings.join("\n"))
}

pub(crate) fn verify_two_complete_passes<T, F>(items: &[T], mut verify: F) -> Result<(), String>
where
    F: FnMut(&T) -> Result<(), String>,
{
    for pass in 1..=2 {
        for item in items {
            verify(item).map_err(|error| {
                format!("complete constellation verification pass {pass} failed: {error}")
            })?;
        }
    }
    Ok(())
}

fn require_equal_observations<T: PartialEq>(
    root: &Path,
    first: &T,
    second: &T,
) -> Result<(), String> {
    if first == second {
        Ok(())
    } else {
        Err(format!(
            "{} moved between two complete recursive worktree observations",
            escaped_path(root)
        ))
    }
}

fn require_equal_repository_boundaries(
    repository: &Path,
    before: &RepositoryBoundaryObservation,
    after: &RepositoryBoundaryObservation,
) -> Result<(), String> {
    if before == after {
        Ok(())
    } else {
        Err(format!(
            "{} moved during recursive child admission or complete boundary re-observation",
            escaped_path(repository)
        ))
    }
}

#[allow(clippy::too_many_lines)] // one complete non-recursive repository-boundary observation
fn observe_repository_boundary(
    repository: &Path,
    scope: &Path,
    include_forced_status: bool,
) -> Result<RepositoryBoundaryObservation, String> {
    require_supported_git_marker(repository, scope)?;
    require_no_executable_git_configuration(repository, scope)?;
    ensure_exact_repository_root(repository, scope)?;
    let head = git_head(repository)?;
    let index = git_bytes(repository, &["ls-files", "--stage", "-z", "--"])?;
    let entries = parse_stage_zero_entries(&index)?;
    let raw_index = read_raw_git_index(repository, scope)?;
    let raw_index_inspection = inspect_optional_raw_git_index(repository, raw_index.as_deref())?;
    require_raw_index_for_inventory(scope, &index, raw_index_inspection.as_ref())?;
    let fsmonitor_extension_present =
        raw_index_inspection.is_some_and(|inspection| inspection.fsmonitor_extension_present);
    let index_materialization = observe_index_materialization(repository, scope, &entries)?;
    let staged_status = git_bytes(repository, STAGED_STATUS_ARGS)?;
    let forced_status = if include_forced_status {
        git_bytes(repository, TRACKED_STATUS_ARGS)?
    } else {
        Vec::new()
    };
    let untracked = git_bytes(
        repository,
        &[
            "ls-files",
            "-z",
            "--others",
            "--exclude-per-directory=.gitignore",
            "--",
        ],
    )?;
    let untracked_gitignores = git_bytes(
        repository,
        &[
            "ls-files",
            "-z",
            "--others",
            "--",
            ".gitignore",
            ":(glob)**/.gitignore",
        ],
    )?;
    let assume_and_skip_flags = git_bytes(repository, &["ls-files", "-v", "-z", "--"])?;
    let (replacement_refs, grafts, mut findings) =
        observe_replacement_authorities(repository, scope)?;
    let (raw_tracked, mut raw_findings) = observe_raw_tracked(repository, scope, &entries)?;
    findings.append(&mut raw_findings);

    let confirmed_index_materialization =
        observe_index_materialization(repository, scope, &entries)?;
    if confirmed_index_materialization != index_materialization {
        return Err(format!(
            "{} index path materialization moved during boundary observation",
            escaped_path(scope)
        ));
    }
    ensure_exact_repository_root(repository, scope)?;
    Ok(RepositoryBoundaryObservation {
        head,
        index,
        raw_index,
        entries,
        index_materialization,
        staged_status,
        forced_status,
        untracked,
        untracked_gitignores,
        assume_and_skip_flags,
        fsmonitor_extension_present,
        replacement_refs,
        grafts,
        raw_tracked,
        findings,
    })
}

#[allow(clippy::too_many_lines)] // one complete canonical recursive Git observation
fn observe_repository(
    repository: &Path,
    scope: &Path,
    expected_head: Option<&str>,
    depth: usize,
) -> Result<RepositorySnapshot, String> {
    if depth > MAX_NESTED_REPOSITORY_DEPTH {
        return Err(format!(
            "nested repository depth exceeds {MAX_NESTED_REPOSITORY_DEPTH} at {}",
            escaped_path(scope)
        ));
    }

    let mut boundary = observe_repository_boundary(repository, scope, false)?;

    let mut gitlinks = Vec::new();
    for entry in boundary
        .entries
        .iter()
        .filter(|entry| entry.mode == "160000")
    {
        let gitlink = Gitlink {
            expected_head: entry.expected_oid.clone(),
            path: entry.path.clone(),
        };
        let child_scope = if scope == Path::new(".") {
            gitlink.path.clone()
        } else {
            scope.join(&gitlink.path)
        };
        let child = match initialized_gitlink(repository, &gitlink.path, &child_scope)? {
            Some(child) => Some(Box::new(observe_repository(
                &child,
                &child_scope,
                Some(&gitlink.expected_head),
                depth + 1,
            )?)),
            None => None,
        };
        gitlinks.push(GitlinkState {
            expected_head: gitlink.expected_head,
            path: gitlink.path,
            child,
        });
    }

    // A forced-visible parent status may inspect initialized submodule
    // worktrees. Admit every initialized descendant first so that status cannot
    // execute a child's untrusted local filter or other executable Git config.
    boundary.forced_status = git_bytes(repository, TRACKED_STATUS_ARGS)?;
    let confirmed_boundary = observe_repository_boundary(repository, scope, true)?;
    require_equal_repository_boundaries(repository, &boundary, &confirmed_boundary)?;

    let RepositoryBoundaryObservation {
        head,
        index,
        raw_index,
        entries: _,
        index_materialization,
        staged_status,
        forced_status,
        untracked,
        untracked_gitignores,
        assume_and_skip_flags,
        fsmonitor_extension_present,
        replacement_refs,
        grafts,
        raw_tracked,
        mut findings,
    } = boundary;

    if let Some(expected) = expected_head {
        if head.as_deref() != Some(expected) {
            findings.push(format!(
                "{}: initialized gitlink HEAD drift (expected {expected}, actual {})",
                escaped_path(scope),
                head.as_deref().unwrap_or("<unborn>"),
            ));
        }
    }
    for path in nul_records(&untracked) {
        findings.push(format!(
            "{}: locally visible untracked path {}",
            escaped_path(scope),
            escape_bytes(path)
        ));
    }
    for path in nul_records(&untracked_gitignores) {
        findings.push(untracked_ignore_policy_finding(scope, path));
    }
    for record in nul_records(&assume_and_skip_flags) {
        if hidden_assume_or_skip_record(record) {
            findings.push(format!(
                "{}: index flag hides worktree state ({})",
                escaped_path(scope),
                escape_bytes(record)
            ));
        }
    }
    if fsmonitor_extension_present {
        findings.push(format!(
            "{}: fsmonitor-valid index state is persisted by the raw FSMN extension and is not admissible during mandatory raw inspection",
            escaped_path(scope)
        ));
    }

    let mut tracked_status = staged_status;
    tracked_status.extend_from_slice(&forced_status);
    if !tracked_status.is_empty() {
        findings.push(format!(
            "{}: tracked, staged, or initialized-submodule dirt ({})",
            escaped_path(scope),
            escape_bytes(&tracked_status)
        ));
    }

    Ok(RepositorySnapshot {
        repository: repository.to_path_buf(),
        scope: scope.to_path_buf(),
        expected_head: expected_head.map(str::to_string),
        head,
        index,
        raw_index,
        tracked_status,
        untracked,
        untracked_gitignores,
        assume_and_skip_flags,
        fsmonitor_extension_present,
        replacement_refs,
        grafts,
        index_materialization,
        raw_tracked,
        gitlinks,
        findings,
    })
}

fn initialized_gitlink(
    repository: &Path,
    relative: &Path,
    scope: &Path,
) -> Result<Option<PathBuf>, String> {
    validate_relative_git_path(relative)?;
    let candidate = repository.join(relative);
    let metadata = match std::fs::symlink_metadata(&candidate) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(error) => {
            return Err(format!(
                "cannot inspect initialized gitlink {}: {error}",
                escaped_path(scope)
            ));
        }
    };
    if is_redirecting_entry(&metadata) || !metadata.is_dir() {
        return Err(format!(
            "gitlink {} is materialized as an unsupported non-directory entry",
            escaped_path(scope)
        ));
    }

    let marker = candidate.join(".git");
    match std::fs::symlink_metadata(&marker) {
        Ok(marker_metadata) => {
            if is_redirecting_entry(&marker_metadata)
                || !(marker_metadata.is_file() || marker_metadata.is_dir())
            {
                return Err(format!(
                    "initialized gitlink {} has an unsupported .git marker",
                    escaped_path(scope)
                ));
            }
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            let empty = std::fs::read_dir(&candidate)
                .map_err(|read_error| {
                    format!(
                        "cannot inspect uninitialized gitlink {}: {read_error}",
                        escaped_path(scope)
                    )
                })?
                .next()
                .is_none();
            if empty {
                return Ok(None);
            }
            return Err(format!(
                "gitlink {} has materialized content but no repository marker",
                escaped_path(scope)
            ));
        }
        Err(error) => {
            return Err(format!(
                "cannot inspect gitlink marker {}: {error}",
                escaped_path(scope)
            ));
        }
    }

    require_no_executable_git_configuration(&candidate, scope)?;

    let candidate = candidate.canonicalize().map_err(|error| {
        format!(
            "cannot canonicalize initialized gitlink {}: {error}",
            escaped_path(scope)
        )
    })?;
    if !candidate.starts_with(repository) || candidate == repository {
        return Err(format!(
            "initialized gitlink {} escapes its containing repository",
            escaped_path(scope)
        ));
    }
    ensure_exact_repository_root(&candidate, scope)?;
    Ok(Some(candidate))
}

fn ensure_exact_repository_root(repository: &Path, scope: &Path) -> Result<(), String> {
    let inside = git_bytes(repository, &["rev-parse", "--is-inside-work-tree"])?;
    let prefix = git_bytes(repository, &["rev-parse", "--show-prefix"])?;
    if inside.as_slice() != b"true\n" || prefix.as_slice() != b"\n" {
        return Err(format!(
            "{} resolves through an ancestor repository instead of its own .git boundary",
            escaped_path(scope)
        ));
    }
    Ok(())
}

fn require_supported_git_marker(repository: &Path, scope: &Path) -> Result<(), String> {
    let marker = repository.join(".git");
    let metadata = std::fs::symlink_metadata(&marker).map_err(|error| {
        format!(
            "cannot inspect repository marker for {}: {error}",
            escaped_path(scope)
        )
    })?;
    if is_redirecting_entry(&metadata) || !(metadata.is_file() || metadata.is_dir()) {
        return Err(format!(
            "repository marker for {} is a link, reparse point, or unsupported entry",
            escaped_path(scope)
        ));
    }
    Ok(())
}

fn ordinary_root_identity(root: &Path) -> Result<MetadataIdentity, String> {
    let metadata = std::fs::symlink_metadata(root).map_err(|error| {
        format!(
            "cannot inspect repository root {}: {error}",
            escaped_path(root)
        )
    })?;
    if is_redirecting_entry(&metadata) || !metadata.is_dir() {
        return Err(format!(
            "repository root {} must be an ordinary directory, not {}",
            escaped_path(root),
            worktree_kind(&metadata).name()
        ));
    }
    Ok(metadata_identity(&metadata))
}

fn require_unchanged_root(
    original_root: &Path,
    canonical_root: &Path,
    expected_identity: &MetadataIdentity,
) -> Result<(), String> {
    let confirmed_identity = ordinary_root_identity(original_root)?;
    let confirmed_root = original_root.canonicalize().map_err(|error| {
        format!(
            "cannot confirm repository root {}: {error}",
            escaped_path(original_root)
        )
    })?;
    if &confirmed_identity != expected_identity || confirmed_root != canonical_root {
        return Err(format!(
            "repository root {} moved during recursive inspection",
            escaped_path(original_root)
        ));
    }
    Ok(())
}

fn observe_replacement_authorities(
    repository: &Path,
    scope: &Path,
) -> Result<(Vec<u8>, Option<Vec<u8>>, Vec<String>), String> {
    let replacement_refs = git_bytes(
        repository,
        &[
            "for-each-ref",
            "--format=%(refname)%00%(objectname)",
            "refs/replace/",
        ],
    )?;
    let grafts_path = git_path_output(
        repository,
        &["rev-parse", "--git-path", "info/grafts"],
        "grafts authority path",
    )?;
    if grafts_path.as_os_str().is_empty() {
        return Err(format!(
            "{}: git returned an empty grafts authority path",
            escaped_path(scope)
        ));
    }
    let grafts_path = if grafts_path.is_absolute() {
        grafts_path
    } else {
        repository.join(grafts_path)
    };
    let grafts = match std::fs::symlink_metadata(&grafts_path) {
        Ok(metadata) => {
            if metadata.file_type().is_symlink() || !metadata.is_file() {
                return Err(format!(
                    "{}: grafts authority {} is not an ordinary file",
                    escaped_path(scope),
                    escaped_path(&grafts_path)
                ));
            }
            let expected_identity = metadata_identity(&metadata);
            let mut bytes = Vec::new();
            let file = std::fs::File::open(&grafts_path).map_err(|error| {
                format!(
                    "cannot open grafts authority {}: {error}",
                    escaped_path(&grafts_path)
                )
            })?;
            let opened = file.metadata().map_err(|error| {
                format!(
                    "cannot inspect opened grafts authority {}: {error}",
                    escaped_path(&grafts_path)
                )
            })?;
            if worktree_kind(&opened) != WorktreeKind::Regular
                || metadata_identity(&opened) != expected_identity
            {
                return Err(format!(
                    "grafts authority {} moved while it was being opened",
                    escaped_path(&grafts_path)
                ));
            }
            file.take(MAX_GRAFTS_BYTES + 1)
                .read_to_end(&mut bytes)
                .map_err(|error| {
                    format!(
                        "cannot read grafts authority {}: {error}",
                        escaped_path(&grafts_path)
                    )
                })?;
            if u64::try_from(bytes.len()).unwrap_or(u64::MAX) > MAX_GRAFTS_BYTES {
                return Err(format!(
                    "grafts authority {} exceeds the {MAX_GRAFTS_BYTES}-byte inspection bound",
                    escaped_path(&grafts_path)
                ));
            }
            let confirmed = std::fs::symlink_metadata(&grafts_path).map_err(|error| {
                format!(
                    "cannot confirm grafts authority {}: {error}",
                    escaped_path(&grafts_path)
                )
            })?;
            if worktree_kind(&confirmed) != WorktreeKind::Regular
                || metadata_identity(&confirmed) != expected_identity
            {
                return Err(format!(
                    "grafts authority {} moved while it was being read",
                    escaped_path(&grafts_path)
                ));
            }
            Some(bytes)
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => None,
        Err(error) => {
            return Err(format!(
                "cannot inspect grafts authority {}: {error}",
                escaped_path(&grafts_path)
            ));
        }
    };

    let mut findings = Vec::new();
    if !replacement_refs.is_empty() {
        findings.push(format!(
            "{}: kind=replace-ref-authority detail={}",
            escaped_path(scope),
            escape_bytes(&replacement_refs)
        ));
    }
    if grafts.as_ref().is_some_and(|bytes| !bytes.is_empty()) {
        findings.push(format!(
            "{}: kind=grafts-authority path={} detail=nonempty local grafts can replace locked commit ancestry",
            escaped_path(scope),
            escaped_path(&grafts_path)
        ));
    }
    Ok((replacement_refs, grafts, findings))
}

fn require_raw_index_for_inventory(
    scope: &Path,
    inventory: &[u8],
    raw_index: Option<&RawIndexInspection>,
) -> Result<(), String> {
    let inventory_count = nul_records(inventory).count();
    match raw_index {
        None if inventory_count == 0 => Ok(()),
        None => Err(format!(
            "{}: Git reported {inventory_count} tracked index entries but the sealed raw primary index is missing",
            escaped_path(scope)
        )),
        Some(raw_index) if raw_index.entry_count == inventory_count => Ok(()),
        Some(raw_index) => Err(format!(
            "{}: raw primary index entry count {} disagrees with Git's complete staged inventory count {inventory_count}",
            escaped_path(scope),
            raw_index.entry_count
        )),
    }
}

fn raw_git_dir_authority(
    repository: &Path,
    scope: &Path,
) -> Result<(PathBuf, RawAuthorityIdentity), String> {
    let git_dir = git_path_output(
        repository,
        &["rev-parse", "--absolute-git-dir"],
        "absolute Git directory",
    )?;
    if !git_dir.is_absolute() {
        return Err(format!(
            "{}: git returned a non-absolute administrative directory {}",
            escaped_path(scope),
            escaped_path(&git_dir)
        ));
    }
    let metadata = std::fs::symlink_metadata(&git_dir).map_err(|error| {
        format!(
            "cannot inspect Git administrative directory {}: {error}",
            escaped_path(&git_dir)
        )
    })?;
    if is_redirecting_entry(&metadata) || !metadata.is_dir() {
        return Err(format!(
            "Git administrative directory {} is a link, reparse point, or unsupported entry",
            escaped_path(&git_dir)
        ));
    }
    let canonical = git_dir.canonicalize().map_err(|error| {
        format!(
            "cannot canonicalize Git administrative directory {}: {error}",
            escaped_path(&git_dir)
        )
    })?;
    #[cfg(unix)]
    if canonical != git_dir {
        return Err(format!(
            "Git administrative directory {} is not reported in canonical form",
            escaped_path(&git_dir)
        ));
    }
    let canonical_metadata = std::fs::symlink_metadata(&canonical).map_err(|error| {
        format!(
            "cannot inspect canonical Git administrative directory {}: {error}",
            escaped_path(&canonical)
        )
    })?;
    let reported_identity = raw_authority_identity(&metadata, &git_dir)?;
    let canonical_identity = raw_authority_identity(&canonical_metadata, &canonical)?;
    if is_redirecting_entry(&canonical_metadata)
        || !canonical_metadata.is_dir()
        || !reported_identity.eq(&canonical_identity)
    {
        return Err(format!(
            "Git administrative directory {} moved during canonicalization",
            escaped_path(&git_dir)
        ));
    }
    Ok((canonical, canonical_identity))
}

fn require_unchanged_raw_git_dir(
    git_dir: &Path,
    expected_identity: &RawAuthorityIdentity,
) -> Result<(), String> {
    let metadata = std::fs::symlink_metadata(git_dir).map_err(|error| {
        format!(
            "cannot confirm Git administrative directory {}: {error}",
            escaped_path(git_dir)
        )
    })?;
    let canonical = git_dir.canonicalize().map_err(|error| {
        format!(
            "cannot recanonicalize Git administrative directory {}: {error}",
            escaped_path(git_dir)
        )
    })?;
    let canonical_metadata = std::fs::symlink_metadata(&canonical).map_err(|error| {
        format!(
            "cannot confirm canonical Git administrative directory {}: {error}",
            escaped_path(&canonical)
        )
    })?;
    let confirmed_identity = raw_authority_identity(&metadata, git_dir)?;
    let canonical_identity = raw_authority_identity(&canonical_metadata, &canonical)?;
    if is_redirecting_entry(&metadata)
        || !metadata.is_dir()
        || !confirmed_identity.eq(expected_identity)
        || is_redirecting_entry(&canonical_metadata)
        || !canonical_metadata.is_dir()
        || !canonical_identity.eq(expected_identity)
    {
        return Err(format!(
            "Git administrative directory {} moved during raw index inspection",
            escaped_path(git_dir)
        ));
    }
    #[cfg(unix)]
    if canonical != git_dir {
        return Err(format!(
            "Git administrative directory {} changed canonical spelling during raw index inspection",
            escaped_path(git_dir)
        ));
    }
    Ok(())
}

fn require_primary_index_path(git_dir: &Path, index_path: &Path) -> Result<(), String> {
    if !index_path.is_absolute() || index_path.file_name() != Some(std::ffi::OsStr::new("index")) {
        return Err(format!(
            "raw primary index path {} is not a direct child of canonical Git directory {}",
            escaped_path(index_path),
            escaped_path(git_dir)
        ));
    }
    let parent = index_path.parent().ok_or_else(|| {
        format!(
            "raw primary index path {} has no administrative parent",
            escaped_path(index_path)
        )
    })?;
    #[cfg(unix)]
    if parent != git_dir {
        return Err(format!(
            "raw primary index path {} is not a direct child of canonical Git directory {}",
            escaped_path(index_path),
            escaped_path(git_dir)
        ));
    }
    #[cfg(not(unix))]
    {
        let canonical_parent = parent.canonicalize().map_err(|error| {
            format!(
                "cannot canonicalize raw primary index parent {}: {error}",
                escaped_path(parent)
            )
        })?;
        if canonical_parent != git_dir {
            return Err(format!(
                "raw primary index path {} does not resolve directly beneath canonical Git directory {}",
                escaped_path(index_path),
                escaped_path(git_dir)
            ));
        }
    }
    Ok(())
}

fn raw_primary_index_path(
    repository: &Path,
    scope: &Path,
    git_dir: &Path,
) -> Result<PathBuf, String> {
    let index_path = git_path_output(
        repository,
        &["rev-parse", "--path-format=absolute", "--git-path", "index"],
        "absolute raw index path",
    )?;
    if index_path.as_os_str().is_empty() {
        return Err(format!(
            "{}: git returned an empty raw index path",
            escaped_path(scope)
        ));
    }
    require_primary_index_path(git_dir, &index_path)?;
    Ok(git_dir.join("index"))
}

fn reconfirm_raw_git_index_authority(
    repository: &Path,
    scope: &Path,
    expected_object_id_bytes: usize,
    expected_git_dir: &Path,
    expected_git_dir_identity: &RawAuthorityIdentity,
    expected_index_path: &Path,
) -> Result<(), String> {
    let confirmed_object_id_bytes = raw_index_object_id_bytes(repository)?;
    if confirmed_object_id_bytes != expected_object_id_bytes {
        return Err(format!(
            "{}: Git object format moved during raw index inspection",
            escaped_path(scope)
        ));
    }
    let (confirmed_git_dir, confirmed_git_dir_identity) = raw_git_dir_authority(repository, scope)?;
    if confirmed_git_dir.as_path() != expected_git_dir
        || !confirmed_git_dir_identity.eq(expected_git_dir_identity)
    {
        return Err(format!(
            "{}: Git administrative authority moved during raw index inspection",
            escaped_path(scope)
        ));
    }
    let confirmed_index_path = raw_primary_index_path(repository, scope, &confirmed_git_dir)?;
    if confirmed_index_path != expected_index_path {
        return Err(format!(
            "{}: raw primary index authority moved during inspection",
            escaped_path(scope)
        ));
    }
    require_unchanged_raw_git_dir(expected_git_dir, expected_git_dir_identity)
}

fn confirm_raw_index_remains_absent(
    repository: &Path,
    scope: &Path,
    object_id_bytes: usize,
    git_dir: &Path,
    git_dir_identity: &RawAuthorityIdentity,
    index_path: &Path,
) -> Result<(), String> {
    reconfirm_raw_git_index_authority(
        repository,
        scope,
        object_id_bytes,
        git_dir,
        git_dir_identity,
        index_path,
    )?;
    match std::fs::symlink_metadata(index_path) {
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Ok(_) => Err(format!(
            "raw Git index {} appeared while its absence was being confirmed",
            escaped_path(index_path)
        )),
        Err(error) => Err(format!(
            "cannot confirm raw Git index absence {}: {error}",
            escaped_path(index_path)
        )),
    }
}

#[allow(clippy::too_many_lines)] // path/handle seals and bounded read are one atomic observation
fn read_raw_git_index(repository: &Path, scope: &Path) -> Result<Option<Vec<u8>>, String> {
    let object_id_bytes = raw_index_object_id_bytes(repository)?;
    let (git_dir, git_dir_identity) = raw_git_dir_authority(repository, scope)?;
    let index_path = raw_primary_index_path(repository, scope, &git_dir)?;
    let path_metadata = match std::fs::symlink_metadata(&index_path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            confirm_raw_index_remains_absent(
                repository,
                scope,
                object_id_bytes,
                &git_dir,
                &git_dir_identity,
                &index_path,
            )?;
            return Ok(None);
        }
        Err(error) => {
            return Err(format!(
                "cannot inspect raw Git index {}: {error}",
                escaped_path(&index_path)
            ));
        }
    };
    if is_redirecting_entry(&path_metadata) || !path_metadata.is_file() {
        return Err(format!(
            "raw Git index {} is a link, reparse point, or unsupported entry",
            escaped_path(&index_path)
        ));
    }
    if path_metadata.len() > MAX_RAW_INDEX_BYTES {
        return Err(format!(
            "raw Git index {} exceeds the {MAX_RAW_INDEX_BYTES}-byte inspection bound",
            escaped_path(&index_path)
        ));
    }

    let expected_identity = raw_authority_identity(&path_metadata, &index_path)?;
    let mut file = match std::fs::File::open(&index_path) {
        Ok(file) => file,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            confirm_raw_index_remains_absent(
                repository,
                scope,
                object_id_bytes,
                &git_dir,
                &git_dir_identity,
                &index_path,
            )?;
            return Err(format!(
                "raw Git index {} disappeared while it was being opened",
                escaped_path(&index_path)
            ));
        }
        Err(error) => {
            return Err(format!(
                "cannot open raw Git index {}: {error}",
                escaped_path(&index_path)
            ));
        }
    };
    let opened_metadata = file.metadata().map_err(|error| {
        format!(
            "cannot inspect opened raw Git index {}: {error}",
            escaped_path(&index_path)
        )
    })?;
    if worktree_kind(&opened_metadata) != WorktreeKind::Regular
        || !raw_authority_identity(&opened_metadata, &index_path)?.eq(&expected_identity)
    {
        return Err(format!(
            "raw Git index {} moved while it was being opened",
            escaped_path(&index_path)
        ));
    }

    let mut bytes = Vec::new();
    std::io::Read::by_ref(&mut file)
        .take(MAX_RAW_INDEX_BYTES + 1)
        .read_to_end(&mut bytes)
        .map_err(|error| {
            format!(
                "cannot read raw Git index {}: {error}",
                escaped_path(&index_path)
            )
        })?;
    if u64::try_from(bytes.len()).unwrap_or(u64::MAX) > MAX_RAW_INDEX_BYTES {
        return Err(format!(
            "raw Git index {} exceeds the {MAX_RAW_INDEX_BYTES}-byte inspection bound",
            escaped_path(&index_path)
        ));
    }
    let opened_after = file.metadata().map_err(|error| {
        format!(
            "cannot confirm opened raw Git index {}: {error}",
            escaped_path(&index_path)
        )
    })?;
    let path_after = std::fs::symlink_metadata(&index_path).map_err(|error| {
        format!(
            "cannot confirm raw Git index {}: {error}",
            escaped_path(&index_path)
        )
    })?;
    let bytes_len = u64::try_from(bytes.len()).unwrap_or(u64::MAX);
    if worktree_kind(&opened_after) != WorktreeKind::Regular
        || !raw_authority_identity(&opened_after, &index_path)?.eq(&expected_identity)
        || is_redirecting_entry(&path_after)
        || !path_after.is_file()
        || !raw_authority_identity(&path_after, &index_path)?.eq(&expected_identity)
        || bytes_len != path_metadata.len()
    {
        return Err(format!(
            "raw Git index {} moved while it was being read",
            escaped_path(&index_path)
        ));
    }
    reconfirm_raw_git_index_authority(
        repository,
        scope,
        object_id_bytes,
        &git_dir,
        &git_dir_identity,
        &index_path,
    )?;
    Ok(Some(bytes))
}

fn validate_git_index_path_bytes(path: &[u8], is_symlink: bool) -> Result<(), String> {
    if path.is_empty() {
        return Err("Git index path is empty".to_string());
    }
    if path.first() == Some(&b'/') {
        return Err(format!("Git index path {} is rooted", escape_bytes(path)));
    }
    if path.last() == Some(&b'/') {
        return Err(format!(
            "Git index path {} has a trailing separator",
            escape_bytes(path)
        ));
    }
    if contains_git_ignored_unicode(path) {
        return Err(format!(
            "Git index path {} contains a Unicode scalar ignored by Git/HFS canonicalization",
            escape_bytes(path)
        ));
    }

    #[cfg(windows)]
    {
        if path.contains(&b'\\') {
            return Err(format!(
                "Git index path {} contains a Windows separator or rooted prefix",
                escape_bytes(path)
            ));
        }
        if path.contains(&b':') {
            return Err(format!(
                "Git index path {} contains a Windows drive or alternate-stream prefix",
                escape_bytes(path)
            ));
        }
    }

    for component in path.split(|byte| *byte == b'/') {
        if component.is_empty() {
            return Err(format!(
                "Git index path {} contains an empty component",
                escape_bytes(path)
            ));
        }
        if component == b"." || component == b".." {
            return Err(format!(
                "Git index path {} contains a dot traversal component",
                escape_bytes(path)
            ));
        }
        if component.eq_ignore_ascii_case(b".git") {
            return Err(format!(
                "Git index path {} contains the administrative .git component",
                escape_bytes(path)
            ));
        }
        if is_symlink && component.eq_ignore_ascii_case(b".gitmodules") {
            return Err(format!(
                "Git index symlink path {} contains the administrative .gitmodules component",
                escape_bytes(path)
            ));
        }
        #[cfg(windows)]
        if windows_git_admin_alias(component, is_symlink) {
            return Err(format!(
                "Git index path {} contains a Windows administrative alias",
                escape_bytes(path)
            ));
        }
    }
    Ok(())
}

fn contains_git_ignored_unicode(path: &[u8]) -> bool {
    path.windows(3).any(|scalar| {
        (scalar[0] == 0xe2 && scalar[1] == 0x80 && matches!(scalar[2], 0x8c..=0x8f | 0xaa..=0xae))
            || (scalar[0] == 0xe2 && scalar[1] == 0x81 && matches!(scalar[2], 0xaa..=0xaf))
            || (scalar[0] == 0xef && scalar[1] == 0xbb && scalar[2] == 0xbf)
    })
}

#[cfg(windows)]
fn windows_git_admin_alias(component: &[u8], is_symlink: bool) -> bool {
    let trimmed_length = component
        .iter()
        .rposition(|byte| !matches!(byte, b'.' | b' '))
        .map_or(0, |index| index + 1);
    let trimmed = &component[..trimmed_length];
    if trimmed.eq_ignore_ascii_case(b".git") {
        return true;
    }
    let dot_git_short_name = trimmed.len() > 4
        && trimmed[..4].eq_ignore_ascii_case(b"git~")
        && trimmed[4..].len() <= 6
        && trimmed[4..].iter().all(u8::is_ascii_digit);
    if dot_git_short_name || !is_symlink {
        return dot_git_short_name;
    }
    if trimmed.eq_ignore_ascii_case(b".gitmodules") {
        return true;
    }
    windows_gitmodules_short_alias(trimmed) || windows_ntfs_fallback_alias(trimmed, b"gi7eba")
}

#[cfg(windows)]
fn windows_gitmodules_short_alias(component: &[u8]) -> bool {
    component.len() == 8
        && component[..7].eq_ignore_ascii_case(b"gitmod~")
        && matches!(component[7], b'1'..=b'4')
}

#[cfg(windows)]
fn windows_ntfs_fallback_alias(component: &[u8], prefix: &[u8]) -> bool {
    if component.len() != 8 || prefix.len() != 6 {
        return false;
    }
    let mut saw_tilde = false;
    for (index, byte) in component.iter().copied().enumerate() {
        if saw_tilde {
            if !byte.is_ascii_digit() {
                return false;
            }
        } else if byte == b'~' {
            if index + 1 >= component.len() || !matches!(component[index + 1], b'1'..=b'9') {
                return false;
            }
            saw_tilde = true;
        } else if index >= prefix.len() || !byte.eq_ignore_ascii_case(&prefix[index]) {
            return false;
        }
    }
    saw_tilde
}

fn validate_relative_git_path(path: &Path) -> Result<(), String> {
    if path.as_os_str().is_empty()
        || path
            .components()
            .any(|component| !matches!(component, Component::Normal(_)))
    {
        return Err(format!(
            "git emitted a non-relative or traversing index path {}",
            escaped_path(path)
        ));
    }
    Ok(())
}

fn observe_index_materialization(
    repository: &Path,
    scope: &Path,
    entries: &[StageZeroEntry],
) -> Result<Vec<IndexMaterializationState>, String> {
    let mut states = BTreeMap::<PathBuf, IndexMaterializationState>::new();
    #[cfg(unix)]
    let mut identity_owners = BTreeMap::<(u64, u64), PathBuf>::new();

    for entry in entries {
        validate_relative_git_path(&entry.path)?;
        let mut components = entry.path.components().peekable();
        let mut logical_path = PathBuf::new();
        let mut materialized_path = repository.to_path_buf();
        while let Some(component) = components.next() {
            let Component::Normal(name) = component else {
                return Err(format!(
                    "git emitted a non-relative or traversing index path {}",
                    escaped_path(&entry.path)
                ));
            };
            #[cfg(not(unix))]
            let parent = materialized_path.clone();
            logical_path.push(name);
            materialized_path.push(name);
            let metadata = match std::fs::symlink_metadata(&materialized_path) {
                Ok(metadata) => metadata,
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => break,
                Err(error) => {
                    return Err(format!(
                        "cannot inspect tracked index prefix {}: {error}",
                        escaped_scoped_path(scope, &logical_path)
                    ));
                }
            };

            #[cfg(not(unix))]
            require_exact_directory_component(&parent, name, scope, &logical_path)?;

            let kind = worktree_kind(&metadata);
            let is_final = components.peek().is_none();
            let is_redirecting = is_redirecting_entry(&metadata);
            if (is_redirecting && (!is_final || kind != WorktreeKind::Symlink))
                || (!is_final && kind != WorktreeKind::Directory)
            {
                return Err(format!(
                    "tracked index prefix {} is materialized as {}; ancestor links, directory reparse points, and non-directories are not admissible",
                    escaped_scoped_path(scope, &logical_path),
                    kind.name()
                ));
            }

            let metadata_identity = metadata_identity(&metadata);
            if let Some(previous) = states.get(&logical_path) {
                if previous.kind != kind || previous.metadata_identity != metadata_identity {
                    return Err(format!(
                        "tracked index prefix {} moved during materialization inspection",
                        escaped_scoped_path(scope, &logical_path)
                    ));
                }
            } else {
                states.insert(
                    logical_path.clone(),
                    IndexMaterializationState {
                        logical_path: logical_path.clone(),
                        kind,
                        metadata_identity: metadata_identity.clone(),
                    },
                );
            }

            #[cfg(unix)]
            {
                let identity = (metadata_identity.device, metadata_identity.inode);
                if let Some(owner) = identity_owners.get(&identity) {
                    if owner != &logical_path {
                        return Err(format!(
                            "tracked index prefixes {} and {} resolve to one filesystem identity; case-folding, normalization, and hard-link aliases are not admissible",
                            escaped_scoped_path(scope, owner),
                            escaped_scoped_path(scope, &logical_path)
                        ));
                    }
                } else {
                    identity_owners.insert(identity, logical_path.clone());
                }
            }
        }
    }

    Ok(states.into_values().collect())
}

#[cfg(not(unix))]
fn require_exact_directory_component(
    parent: &Path,
    expected_name: &std::ffi::OsStr,
    scope: &Path,
    logical_path: &Path,
) -> Result<(), String> {
    let mut exact_matches = 0_u8;
    let entries = std::fs::read_dir(parent).map_err(|error| {
        format!(
            "cannot enumerate tracked index prefix parent {}: {error}",
            escaped_path(parent)
        )
    })?;
    for entry in entries {
        let entry = entry.map_err(|error| {
            format!(
                "cannot enumerate tracked index prefix parent {}: {error}",
                escaped_path(parent)
            )
        })?;
        if entry.file_name() == expected_name {
            exact_matches = exact_matches.saturating_add(1);
        }
    }
    if exact_matches == 1 {
        Ok(())
    } else {
        Err(format!(
            "tracked index prefix {} is not materialized with one exact directory-entry spelling",
            escaped_scoped_path(scope, logical_path)
        ))
    }
}

#[cfg(windows)]
pub(crate) fn is_redirecting_entry(metadata: &Metadata) -> bool {
    use std::os::windows::fs::MetadataExt as _;

    const FILE_ATTRIBUTE_REPARSE_POINT: u32 = 0x400;
    metadata.file_attributes() & FILE_ATTRIBUTE_REPARSE_POINT != 0
}

#[cfg(not(windows))]
pub(crate) fn is_redirecting_entry(metadata: &Metadata) -> bool {
    metadata.file_type().is_symlink()
}

#[allow(clippy::too_many_lines)] // raw type, mode, and object state are one atomic observation
fn observe_raw_tracked(
    repository: &Path,
    scope: &Path,
    entries: &[StageZeroEntry],
) -> Result<(Vec<RawTrackedState>, Vec<String>), String> {
    let mut states = Vec::new();
    let mut regular_indices = Vec::new();

    for entry in entries {
        let candidate = repository.join(&entry.path);
        match entry.mode.as_str() {
            "100644" | "100755" => {
                let (actual_kind, actual_mode, observed_identity) =
                    match std::fs::symlink_metadata(&candidate) {
                        Ok(metadata) => {
                            let kind = worktree_kind(&metadata);
                            let mode = if kind == WorktreeKind::Regular {
                                regular_worktree_mode(&metadata)
                            } else {
                                None
                            };
                            (kind, mode, Some(metadata_identity(&metadata)))
                        }
                        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                            (WorktreeKind::Missing, None, None)
                        }
                        Err(error) => {
                            return Err(format!(
                                "cannot inspect tracked path {}: {error}",
                                escaped_scoped_path(scope, &entry.path)
                            ));
                        }
                    };
                let index = states.len();
                states.push(RawTrackedState {
                    expected_mode: entry.mode.clone(),
                    expected_oid: entry.expected_oid.clone(),
                    path: entry.path.clone(),
                    actual_kind,
                    actual_mode,
                    actual_oid: None,
                    link_target: None,
                    metadata_identity: observed_identity,
                });
                if actual_kind == WorktreeKind::Regular {
                    regular_indices.push(index);
                }
            }
            "120000" => {
                let metadata = match std::fs::symlink_metadata(&candidate) {
                    Ok(metadata) => Some(metadata),
                    Err(error) if error.kind() == std::io::ErrorKind::NotFound => None,
                    Err(error) => {
                        return Err(format!(
                            "cannot inspect tracked symlink {}: {error}",
                            escaped_scoped_path(scope, &entry.path)
                        ));
                    }
                };
                let actual_kind = metadata
                    .as_ref()
                    .map_or(WorktreeKind::Missing, worktree_kind);
                let observed_identity = metadata.as_ref().map(metadata_identity);
                let (actual_mode, actual_oid, link_target) = if actual_kind == WorktreeKind::Symlink
                {
                    let target = std::fs::read_link(&candidate).map_err(|error| {
                        format!(
                            "cannot read tracked symlink {}: {error}",
                            escaped_scoped_path(scope, &entry.path)
                        )
                    })?;
                    let target = symlink_target_bytes(&target)?;
                    let oid = hash_symlink_target(repository, &target)?;
                    confirm_symlink_unchanged(
                        &candidate,
                        observed_identity.as_ref().ok_or_else(|| {
                            "tracked symlink identity was not captured".to_string()
                        })?,
                        &target,
                        scope,
                        &entry.path,
                    )?;
                    (Some("120000".to_string()), Some(oid), Some(target))
                } else {
                    (None, None, None)
                };
                states.push(RawTrackedState {
                    expected_mode: entry.mode.clone(),
                    expected_oid: entry.expected_oid.clone(),
                    path: entry.path.clone(),
                    actual_kind,
                    actual_mode,
                    actual_oid,
                    link_target,
                    metadata_identity: observed_identity,
                });
            }
            "160000" => {}
            _ => states.push(RawTrackedState {
                expected_mode: entry.mode.clone(),
                expected_oid: entry.expected_oid.clone(),
                path: entry.path.clone(),
                actual_kind: WorktreeKind::UnsupportedIndexMode,
                actual_mode: None,
                actual_oid: None,
                link_target: None,
                metadata_identity: None,
            }),
        }
    }

    populate_regular_object_ids(repository, scope, &mut states, &regular_indices)?;

    let findings = raw_tracked_findings(scope, &states);
    Ok((states, findings))
}

fn raw_tracked_findings(scope: &Path, states: &[RawTrackedState]) -> Vec<String> {
    let mut findings = Vec::new();
    for state in states {
        let path = escaped_scoped_path(scope, &state.path);
        if state.actual_kind == WorktreeKind::UnsupportedIndexMode {
            findings.push(format!(
                "{path}: kind=unsupported-stage-zero-mode detail=mode {}",
                state.expected_mode
            ));
            continue;
        }
        let expected_kind = if state.expected_mode == "120000" {
            WorktreeKind::Symlink
        } else {
            WorktreeKind::Regular
        };
        if state.actual_kind != expected_kind {
            findings.push(format!(
                "{path}: kind=raw-tracked-source-mismatch detail=type expected={} actual={}",
                expected_kind.name(),
                state.actual_kind.name()
            ));
            continue;
        }
        if raw_mode_mismatch(
            state.actual_kind,
            &state.expected_mode,
            state.actual_mode.as_deref(),
            cfg!(unix),
        ) {
            findings.push(format!(
                "{path}: kind=raw-tracked-source-mismatch detail=mode expected={} actual={}",
                state.expected_mode,
                state.actual_mode.as_deref().unwrap_or("<unavailable>")
            ));
        }
        if state.actual_oid.as_deref() != Some(state.expected_oid.as_str()) {
            findings.push(format!(
                "{path}: kind=raw-tracked-source-mismatch detail=object expected={} actual={}",
                state.expected_oid,
                state.actual_oid.as_deref().unwrap_or("<unavailable>")
            ));
        }
    }
    findings
}

fn raw_mode_mismatch(
    actual_kind: WorktreeKind,
    expected: &str,
    actual: Option<&str>,
    regular_mode_is_observable: bool,
) -> bool {
    match actual {
        Some(actual) => actual != expected,
        None => actual_kind != WorktreeKind::Regular || regular_mode_is_observable,
    }
}

fn populate_regular_object_ids(
    repository: &Path,
    scope: &Path,
    states: &mut [RawTrackedState],
    regular_indices: &[usize],
) -> Result<(), String> {
    let mut batch = Vec::new();
    let mut argument_bytes: usize = 0;
    for &index in regular_indices {
        let path_bytes = git_path_argument_bytes(&states[index].path)?.saturating_add(1);
        if path_bytes > MAX_RAW_HASH_ARGUMENT_BYTES {
            return Err(format!(
                "tracked path exceeds raw hash argument bound: {}",
                escaped_scoped_path(scope, &states[index].path)
            ));
        }
        if !batch.is_empty()
            && (batch.len() == MAX_RAW_HASH_PATHS_PER_BATCH
                || argument_bytes.saturating_add(path_bytes) > MAX_RAW_HASH_ARGUMENT_BYTES)
        {
            hash_regular_batch(repository, scope, states, &batch)?;
            batch.clear();
            argument_bytes = 0;
        }
        batch.push(index);
        argument_bytes += path_bytes;
    }
    if !batch.is_empty() {
        hash_regular_batch(repository, scope, states, &batch)?;
    }
    Ok(())
}

fn hash_regular_batch(
    repository: &Path,
    scope: &Path,
    states: &mut [RawTrackedState],
    batch: &[usize],
) -> Result<(), String> {
    let mut command = observation_git_command(repository, &["hash-object", "--no-filters", "--"]);
    for &index in batch {
        command.arg(&states[index].path);
    }
    let output = command.output().map_err(|error| {
        format!(
            "git hash-object raw batch in {} failed to spawn: {error}",
            escaped_path(repository)
        )
    })?;
    if !output.status.success() {
        return Err(format!(
            "git hash-object raw batch in {} failed: {}",
            escaped_path(repository),
            escaped_git_stderr(&output.stderr)
        ));
    }
    let object_ids = parse_object_ids(&output.stdout, batch.len(), "raw tracked batch")?;
    for (&index, object_id) in batch.iter().zip(object_ids) {
        let candidate = repository.join(&states[index].path);
        let confirmed = std::fs::symlink_metadata(&candidate).map_err(|error| {
            format!(
                "cannot confirm tracked path {} after raw hashing: {error}",
                escaped_scoped_path(scope, &states[index].path)
            )
        })?;
        let confirmed_identity = metadata_identity(&confirmed);
        if worktree_kind(&confirmed) != WorktreeKind::Regular
            || states[index].metadata_identity.as_ref() != Some(&confirmed_identity)
        {
            return Err(format!(
                "tracked path {} moved while its raw bytes were being hashed",
                escaped_scoped_path(scope, &states[index].path)
            ));
        }
        states[index].actual_oid = Some(object_id);
    }
    Ok(())
}

fn confirm_symlink_unchanged(
    candidate: &Path,
    expected_identity: &MetadataIdentity,
    expected_target: &[u8],
    scope: &Path,
    relative: &Path,
) -> Result<(), String> {
    let metadata = std::fs::symlink_metadata(candidate).map_err(|error| {
        format!(
            "cannot confirm tracked symlink {}: {error}",
            escaped_scoped_path(scope, relative)
        )
    })?;
    if worktree_kind(&metadata) != WorktreeKind::Symlink
        || &metadata_identity(&metadata) != expected_identity
    {
        return Err(format!(
            "tracked symlink {} moved while its target was being hashed",
            escaped_scoped_path(scope, relative)
        ));
    }
    let target = std::fs::read_link(candidate).map_err(|error| {
        format!(
            "cannot confirm tracked symlink target {}: {error}",
            escaped_scoped_path(scope, relative)
        )
    })?;
    if symlink_target_bytes(&target)? != expected_target {
        return Err(format!(
            "tracked symlink {} changed target while being hashed",
            escaped_scoped_path(scope, relative)
        ));
    }
    Ok(())
}

fn hash_symlink_target(repository: &Path, target: &[u8]) -> Result<String, String> {
    let mut child =
        observation_git_command(repository, &["hash-object", "--no-filters", "--stdin"])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|error| format!("git hash-object for symlink failed to spawn: {error}"))?;
    let write_result = child
        .stdin
        .take()
        .ok_or_else(|| "git hash-object did not provide piped stdin".to_string())?
        .write_all(target);
    let output = child
        .wait_with_output()
        .map_err(|error| format!("git hash-object for symlink failed to wait: {error}"))?;
    write_result.map_err(|error| format!("cannot stream symlink target to git: {error}"))?;
    if !output.status.success() {
        return Err(format!(
            "git hash-object for symlink in {} failed: {}",
            escaped_path(repository),
            escaped_git_stderr(&output.stderr)
        ));
    }
    parse_object_ids(&output.stdout, 1, "symlink target")?
        .into_iter()
        .next()
        .ok_or_else(|| "git hash-object omitted the symlink target object ID".to_string())
}

fn parse_object_ids(bytes: &[u8], expected: usize, context: &str) -> Result<Vec<String>, String> {
    let text = std::str::from_utf8(bytes)
        .map_err(|error| format!("git hash-object {context} returned non-UTF-8 text: {error}"))?;
    let object_ids: Vec<_> = text.lines().map(str::to_string).collect();
    if object_ids.len() != expected {
        return Err(format!(
            "git hash-object {context} returned {} object IDs; expected {expected}",
            object_ids.len()
        ));
    }
    if object_ids
        .iter()
        .any(|object_id| !valid_object_id(object_id))
    {
        return Err(format!(
            "git hash-object {context} returned a malformed object ID"
        ));
    }
    Ok(object_ids)
}

fn valid_object_id(object_id: &str) -> bool {
    matches!(object_id.len(), 40 | 64)
        && object_id
            .bytes()
            .all(|byte| byte.is_ascii_digit() || matches!(byte, b'a'..=b'f'))
}

fn worktree_kind(metadata: &Metadata) -> WorktreeKind {
    let file_type = metadata.file_type();
    if file_type.is_symlink() {
        WorktreeKind::Symlink
    } else if file_type.is_file() {
        WorktreeKind::Regular
    } else if file_type.is_dir() {
        WorktreeKind::Directory
    } else {
        WorktreeKind::Other
    }
}

#[cfg(unix)]
fn metadata_identity(metadata: &Metadata) -> MetadataIdentity {
    use std::os::unix::fs::MetadataExt as _;

    MetadataIdentity {
        device: metadata.dev(),
        inode: metadata.ino(),
        mode: metadata.mode(),
        size: metadata.size(),
        modified_seconds: metadata.mtime(),
        modified_nanoseconds: metadata.mtime_nsec(),
        changed_seconds: metadata.ctime(),
        changed_nanoseconds: metadata.ctime_nsec(),
    }
}

#[cfg(not(unix))]
fn metadata_identity(metadata: &Metadata) -> MetadataIdentity {
    MetadataIdentity {
        kind: worktree_kind(metadata),
        size: metadata.len(),
        readonly: metadata.permissions().readonly(),
    }
}

#[cfg(unix)]
fn raw_authority_identity(
    metadata: &Metadata,
    _authority: &Path,
) -> Result<RawAuthorityIdentity, String> {
    Ok(metadata_identity(metadata))
}

#[cfg(windows)]
fn raw_authority_identity(
    metadata: &Metadata,
    authority: &Path,
) -> Result<RawAuthorityIdentity, String> {
    use std::os::windows::fs::MetadataExt as _;

    let volume_serial_number = metadata.volume_serial_number().ok_or_else(|| {
        format!(
            "Windows filesystem did not expose a volume serial number for raw Git authority {}",
            escaped_path(authority)
        )
    })?;
    let file_index = metadata.file_index().ok_or_else(|| {
        format!(
            "Windows filesystem did not expose a file index for raw Git authority {}",
            escaped_path(authority)
        )
    })?;
    let links = metadata.number_of_links().ok_or_else(|| {
        format!(
            "Windows filesystem did not expose a link count for raw Git authority {}",
            escaped_path(authority)
        )
    })?;
    Ok(RawAuthorityIdentity {
        volume_serial_number,
        file_index,
        attributes: metadata.file_attributes(),
        links,
        size: metadata.file_size(),
        created: metadata.creation_time(),
        modified: metadata.last_write_time(),
    })
}

#[cfg(not(any(unix, windows)))]
fn raw_authority_identity(
    _metadata: &Metadata,
    authority: &Path,
) -> Result<RawAuthorityIdentity, String> {
    Err(format!(
        "raw Git authority sealing is unsupported on this platform for {}: std exposes no stable file-object identity",
        escaped_path(authority)
    ))
}

#[cfg(unix)]
const fn git_mode_from_unix_permissions(mode: u32) -> &'static str {
    if mode & 0o100 == 0 {
        "100644"
    } else {
        "100755"
    }
}

#[cfg(unix)]
fn regular_worktree_mode(metadata: &Metadata) -> Option<String> {
    use std::os::unix::fs::PermissionsExt as _;

    Some(git_mode_from_unix_permissions(metadata.permissions().mode()).to_string())
}

#[cfg(not(unix))]
fn regular_worktree_mode(_metadata: &Metadata) -> Option<String> {
    None
}

#[cfg(unix)]
fn symlink_target_bytes(target: &Path) -> Result<Vec<u8>, String> {
    use std::os::unix::ffi::OsStrExt as _;

    Ok(target.as_os_str().as_bytes().to_vec())
}

#[cfg(not(unix))]
fn symlink_target_bytes(target: &Path) -> Result<Vec<u8>, String> {
    target
        .to_str()
        .map(|target| target.as_bytes().to_vec())
        .ok_or_else(|| "non-UTF-8 symlink targets are unsupported and refuse closed".to_string())
}

#[cfg(unix)]
fn git_path_from_bytes(bytes: &[u8]) -> Result<PathBuf, String> {
    use std::ffi::OsString;
    use std::os::unix::ffi::OsStringExt as _;

    Ok(PathBuf::from(OsString::from_vec(bytes.to_vec())))
}

#[cfg(not(unix))]
fn git_path_from_bytes(bytes: &[u8]) -> Result<PathBuf, String> {
    std::str::from_utf8(bytes)
        .map(PathBuf::from)
        .map_err(|_| "non-UTF-8 tracked paths are unsupported and refuse closed".to_string())
}

#[cfg(unix)]
fn git_path_argument_bytes(path: &Path) -> Result<usize, String> {
    use std::os::unix::ffi::OsStrExt as _;

    Ok(path.as_os_str().as_bytes().len())
}

#[cfg(not(unix))]
fn git_path_argument_bytes(path: &Path) -> Result<usize, String> {
    path.to_str()
        .map(str::len)
        .ok_or_else(|| "non-UTF-8 tracked paths are unsupported and refuse closed".to_string())
}

fn scoped_path(scope: &Path, relative: &Path) -> PathBuf {
    if scope == Path::new(".") {
        relative.to_path_buf()
    } else {
        scope.join(relative)
    }
}

fn escaped_scoped_path(scope: &Path, relative: &Path) -> String {
    escaped_path(&scoped_path(scope, relative))
}

fn untracked_ignore_policy_finding(scope: &Path, path: &[u8]) -> String {
    format!(
        "{}: kind=untracked-ignore-policy path={} detail=only tracked .gitignore files may define project ignore semantics",
        escaped_path(scope),
        escape_bytes(path)
    )
}

fn escaped_path(path: &Path) -> String {
    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStrExt as _;

        return escape_bytes(path.as_os_str().as_bytes());
    }
    #[cfg(not(unix))]
    {
        escape_bytes(path.to_string_lossy().as_bytes())
    }
}

fn raw_index_object_id_bytes(repository: &Path) -> Result<usize, String> {
    match git_bytes(repository, &["rev-parse", "--show-object-format=storage"])?.as_slice() {
        b"sha1\n" => Ok(20),
        b"sha256\n" => Ok(32),
        bytes => Err(format!(
            "unsupported or ambiguous Git object format for raw index inspection: {}",
            escape_bytes(bytes)
        )),
    }
}

fn inspect_optional_raw_git_index(
    repository: &Path,
    raw_index: Option<&[u8]>,
) -> Result<Option<RawIndexInspection>, String> {
    let Some(raw_index) = raw_index else {
        return Ok(None);
    };
    inspect_raw_git_index(raw_index, raw_index_object_id_bytes(repository)?).map(Some)
}

fn raw_index_range<'a>(
    index: &'a [u8],
    start: usize,
    length: usize,
    context: &str,
) -> Result<&'a [u8], String> {
    let end = start
        .checked_add(length)
        .ok_or_else(|| format!("raw Git index {context} bounds overflow"))?;
    index
        .get(start..end)
        .ok_or_else(|| format!("raw Git index is truncated in {context}"))
}

fn raw_index_u16(index: &[u8], offset: usize, context: &str) -> Result<u16, String> {
    let bytes = raw_index_range(index, offset, 2, context)?;
    Ok(u16::from_be_bytes([bytes[0], bytes[1]]))
}

fn raw_index_u32(index: &[u8], offset: usize, context: &str) -> Result<u32, String> {
    let bytes = raw_index_range(index, offset, 4, context)?;
    Ok(u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
}

fn digest_bit_length(bytes: &[u8], algorithm: &str) -> Result<u64, String> {
    u64::try_from(bytes.len())
        .ok()
        .and_then(|length| length.checked_mul(8))
        .ok_or_else(|| format!("raw Git index is too large for {algorithm} length encoding"))
}

#[allow(clippy::needless_range_loop)] // the message schedule is an indexed recurrence
fn sha1_compress(state: &mut [u32; 5], block: &[u8]) {
    let mut words = [0_u32; 80];
    for (index, word) in words[..16].iter_mut().enumerate() {
        let offset = index * 4;
        *word = u32::from_be_bytes([
            block[offset],
            block[offset + 1],
            block[offset + 2],
            block[offset + 3],
        ]);
    }
    for index in 16..80 {
        words[index] =
            (words[index - 3] ^ words[index - 8] ^ words[index - 14] ^ words[index - 16])
                .rotate_left(1);
    }

    let [mut a, mut b, mut c, mut d, mut e] = *state;
    for (index, word) in words.iter().enumerate() {
        let (choice, constant) = match index {
            0..20 => ((b & c) | ((!b) & d), 0x5a82_7999),
            20..40 => (b ^ c ^ d, 0x6ed9_eba1),
            40..60 => ((b & c) | (b & d) | (c & d), 0x8f1b_bcdc),
            _ => (b ^ c ^ d, 0xca62_c1d6),
        };
        let next = a
            .rotate_left(5)
            .wrapping_add(choice)
            .wrapping_add(e)
            .wrapping_add(constant)
            .wrapping_add(*word);
        e = d;
        d = c;
        c = b.rotate_left(30);
        b = a;
        a = next;
    }
    state[0] = state[0].wrapping_add(a);
    state[1] = state[1].wrapping_add(b);
    state[2] = state[2].wrapping_add(c);
    state[3] = state[3].wrapping_add(d);
    state[4] = state[4].wrapping_add(e);
}

fn sha1_digest(bytes: &[u8]) -> Result<[u8; 20], String> {
    let bit_length = digest_bit_length(bytes, "SHA-1")?;
    let mut state = [
        0x6745_2301,
        0xefcd_ab89,
        0x98ba_dcfe,
        0x1032_5476,
        0xc3d2_e1f0,
    ];
    let mut chunks = bytes.chunks_exact(64);
    for block in chunks.by_ref() {
        sha1_compress(&mut state, block);
    }
    let remainder = chunks.remainder();
    let mut tail = [0_u8; 128];
    tail[..remainder.len()].copy_from_slice(remainder);
    tail[remainder.len()] = 0x80;
    let padded_length = if remainder.len() < 56 { 64 } else { 128 };
    tail[padded_length - 8..padded_length].copy_from_slice(&bit_length.to_be_bytes());
    for block in tail[..padded_length].chunks_exact(64) {
        sha1_compress(&mut state, block);
    }

    let mut digest = [0_u8; 20];
    for (word, output) in state.iter().zip(digest.chunks_exact_mut(4)) {
        output.copy_from_slice(&word.to_be_bytes());
    }
    Ok(digest)
}

const SHA256_ROUND_CONSTANTS: [u32; 64] = [
    0x428a_2f98,
    0x7137_4491,
    0xb5c0_fbcf,
    0xe9b5_dba5,
    0x3956_c25b,
    0x59f1_11f1,
    0x923f_82a4,
    0xab1c_5ed5,
    0xd807_aa98,
    0x1283_5b01,
    0x2431_85be,
    0x550c_7dc3,
    0x72be_5d74,
    0x80de_b1fe,
    0x9bdc_06a7,
    0xc19b_f174,
    0xe49b_69c1,
    0xefbe_4786,
    0x0fc1_9dc6,
    0x240c_a1cc,
    0x2de9_2c6f,
    0x4a74_84aa,
    0x5cb0_a9dc,
    0x76f9_88da,
    0x983e_5152,
    0xa831_c66d,
    0xb003_27c8,
    0xbf59_7fc7,
    0xc6e0_0bf3,
    0xd5a7_9147,
    0x06ca_6351,
    0x1429_2967,
    0x27b7_0a85,
    0x2e1b_2138,
    0x4d2c_6dfc,
    0x5338_0d13,
    0x650a_7354,
    0x766a_0abb,
    0x81c2_c92e,
    0x9272_2c85,
    0xa2bf_e8a1,
    0xa81a_664b,
    0xc24b_8b70,
    0xc76c_51a3,
    0xd192_e819,
    0xd699_0624,
    0xf40e_3585,
    0x106a_a070,
    0x19a4_c116,
    0x1e37_6c08,
    0x2748_774c,
    0x34b0_bcb5,
    0x391c_0cb3,
    0x4ed8_aa4a,
    0x5b9c_ca4f,
    0x682e_6ff3,
    0x748f_82ee,
    0x78a5_636f,
    0x84c8_7814,
    0x8cc7_0208,
    0x90be_fffa,
    0xa450_6ceb,
    0xbef9_a3f7,
    0xc671_78f2,
];

#[allow(clippy::needless_range_loop)] // the message schedule is an indexed recurrence
fn sha256_compress(state: &mut [u32; 8], block: &[u8]) {
    let mut words = [0_u32; 64];
    for (index, word) in words[..16].iter_mut().enumerate() {
        let offset = index * 4;
        *word = u32::from_be_bytes([
            block[offset],
            block[offset + 1],
            block[offset + 2],
            block[offset + 3],
        ]);
    }
    for index in 16..64 {
        let sigma0 = words[index - 15].rotate_right(7)
            ^ words[index - 15].rotate_right(18)
            ^ (words[index - 15] >> 3);
        let sigma1 = words[index - 2].rotate_right(17)
            ^ words[index - 2].rotate_right(19)
            ^ (words[index - 2] >> 10);
        words[index] = words[index - 16]
            .wrapping_add(sigma0)
            .wrapping_add(words[index - 7])
            .wrapping_add(sigma1);
    }

    let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = *state;
    for (word, constant) in words.iter().zip(SHA256_ROUND_CONSTANTS) {
        let sum1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
        let choice = (e & f) ^ ((!e) & g);
        let temporary1 = h
            .wrapping_add(sum1)
            .wrapping_add(choice)
            .wrapping_add(constant)
            .wrapping_add(*word);
        let sum0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
        let majority = (a & b) ^ (a & c) ^ (b & c);
        let temporary2 = sum0.wrapping_add(majority);
        h = g;
        g = f;
        f = e;
        e = d.wrapping_add(temporary1);
        d = c;
        c = b;
        b = a;
        a = temporary1.wrapping_add(temporary2);
    }
    state[0] = state[0].wrapping_add(a);
    state[1] = state[1].wrapping_add(b);
    state[2] = state[2].wrapping_add(c);
    state[3] = state[3].wrapping_add(d);
    state[4] = state[4].wrapping_add(e);
    state[5] = state[5].wrapping_add(f);
    state[6] = state[6].wrapping_add(g);
    state[7] = state[7].wrapping_add(h);
}

fn sha256_digest(bytes: &[u8]) -> Result<[u8; 32], String> {
    let bit_length = digest_bit_length(bytes, "SHA-256")?;
    let mut state = [
        0x6a09_e667,
        0xbb67_ae85,
        0x3c6e_f372,
        0xa54f_f53a,
        0x510e_527f,
        0x9b05_688c,
        0x1f83_d9ab,
        0x5be0_cd19,
    ];
    let mut chunks = bytes.chunks_exact(64);
    for block in chunks.by_ref() {
        sha256_compress(&mut state, block);
    }
    let remainder = chunks.remainder();
    let mut tail = [0_u8; 128];
    tail[..remainder.len()].copy_from_slice(remainder);
    tail[remainder.len()] = 0x80;
    let padded_length = if remainder.len() < 56 { 64 } else { 128 };
    tail[padded_length - 8..padded_length].copy_from_slice(&bit_length.to_be_bytes());
    for block in tail[..padded_length].chunks_exact(64) {
        sha256_compress(&mut state, block);
    }

    let mut digest = [0_u8; 32];
    for (word, output) in state.iter().zip(digest.chunks_exact_mut(4)) {
        output.copy_from_slice(&word.to_be_bytes());
    }
    Ok(digest)
}

fn authenticate_raw_index_checksum(
    index: &[u8],
    checksum_start: usize,
    object_id_bytes: usize,
) -> Result<(), String> {
    let payload = raw_index_range(index, 0, checksum_start, "checksummed payload")?;
    let checksum = raw_index_range(index, checksum_start, object_id_bytes, "trailing checksum")?;
    let authentic = match object_id_bytes {
        20 => sha1_digest(payload)?.as_slice() == checksum,
        32 => sha256_digest(payload)?.as_slice() == checksum,
        _ => false,
    };
    if authentic {
        Ok(())
    } else {
        Err("raw Git index trailing checksum does not authenticate its payload".to_string())
    }
}

/// Inspect the primary on-disk index without asking Git to reconstruct its
/// in-memory CE_FSMONITOR_VALID bits. Git does not serialize that per-entry
/// bit; an FSMN extension instead carries a bitmap of entries that are dirty.
/// Any FSMN authority is inadmissible here because this verifier performs a
/// complete raw worktree inspection. Split indexes refuse closed because the
/// primary index alone cannot account for extensions in its shared base.
#[allow(clippy::too_many_lines)] // entry-boundary proof and extension scan must remain one parser
fn inspect_raw_git_index(
    index: &[u8],
    object_id_bytes: usize,
) -> Result<RawIndexInspection, String> {
    if !matches!(object_id_bytes, 20 | 32) {
        return Err(format!(
            "unsupported raw Git index object ID width: {object_id_bytes} bytes"
        ));
    }
    let minimum_length = 12_usize
        .checked_add(object_id_bytes)
        .ok_or_else(|| "raw Git index minimum length overflow".to_string())?;
    if index.len() < minimum_length {
        return Err("raw Git index is truncated before its checksum".to_string());
    }
    if raw_index_range(index, 0, 4, "signature")? != b"DIRC" {
        return Err("raw Git index has an invalid signature".to_string());
    }
    let version = raw_index_u32(index, 4, "version")?;
    if !matches!(version, 2 | 3) {
        return Err(format!(
            "raw Git index version {version} is unsupported; only unambiguous v2/v3 layouts are admissible"
        ));
    }
    let entry_count = usize::try_from(raw_index_u32(index, 8, "entry count")?)
        .map_err(|_| "raw Git index entry count does not fit this platform".to_string())?;
    let checksum_start = index
        .len()
        .checked_sub(object_id_bytes)
        .ok_or_else(|| "raw Git index checksum bounds underflow".to_string())?;
    let minimum_entry_length = 40_usize
        .checked_add(object_id_bytes)
        .and_then(|length| length.checked_add(3))
        .ok_or_else(|| "raw Git index entry length overflow".to_string())?;
    let maximum_entry_count = checksum_start
        .saturating_sub(12)
        .checked_div(minimum_entry_length)
        .unwrap_or(0);
    if entry_count > maximum_entry_count {
        return Err("raw Git index entry count exceeds its bounded payload".to_string());
    }

    let mut offset = 12_usize;
    let mut entry_order = Vec::with_capacity(entry_count);
    for ordinal in 0..entry_count {
        let entry_start = offset;
        let mode_offset = entry_start
            .checked_add(24)
            .ok_or_else(|| "raw Git index entry mode offset overflow".to_string())?;
        let mode = raw_index_u32(index, mode_offset, "entry mode")?;
        let is_symlink = mode & 0o170_000 == 0o120_000;
        let flags_offset = entry_start
            .checked_add(40)
            .and_then(|value| value.checked_add(object_id_bytes))
            .ok_or_else(|| "raw Git index entry flags offset overflow".to_string())?;
        let flags = raw_index_u16(index, flags_offset, "entry flags")?;
        offset = flags_offset
            .checked_add(2)
            .ok_or_else(|| "raw Git index entry header overflow".to_string())?;
        if offset > checksum_start {
            return Err(format!(
                "raw Git index entry {ordinal} overlaps its checksum"
            ));
        }

        if flags & 0x4000 != 0 {
            if version != 3 {
                return Err(format!(
                    "raw Git index v{version} entry {ordinal} uses v3 extended flags"
                ));
            }
            let extended = raw_index_u16(index, offset, "extended entry flags")?;
            // Only CE_INTENT_TO_ADD and CE_SKIP_WORKTREE are serialized in
            // the v3 extended word. In particular, the in-memory-only
            // CE_FSMONITOR_VALID bit position must never be interpreted here.
            if extended & !0x6000 != 0 {
                return Err(format!(
                    "raw Git index entry {ordinal} has unsupported extended flags 0x{extended:04x}"
                ));
            }
            offset = offset
                .checked_add(2)
                .ok_or_else(|| "raw Git index extended flags overflow".to_string())?;
            if offset > checksum_start {
                return Err(format!(
                    "raw Git index entry {ordinal} extended flags overlap its checksum"
                ));
            }
        }

        let path_start = offset;
        let relative_nul = index
            .get(path_start..checksum_start)
            .and_then(|bytes| bytes.iter().position(|byte| *byte == 0))
            .ok_or_else(|| {
                format!("raw Git index entry {ordinal} has no bounded path terminator")
            })?;
        let encoded_path_length = usize::from(flags & 0x0fff);
        if (encoded_path_length < 0x0fff && relative_nul != encoded_path_length)
            || (encoded_path_length == 0x0fff && relative_nul < encoded_path_length)
        {
            return Err(format!(
                "raw Git index entry {ordinal} path length disagrees with its flags"
            ));
        }
        let nul_offset = path_start
            .checked_add(relative_nul)
            .ok_or_else(|| "raw Git index path terminator overflow".to_string())?;
        let consumed = nul_offset
            .checked_add(1)
            .and_then(|value| value.checked_sub(entry_start))
            .ok_or_else(|| "raw Git index entry length overflow".to_string())?;
        let padded_length = consumed
            .checked_add(7)
            .map(|value| value & !7)
            .ok_or_else(|| "raw Git index entry padding overflow".to_string())?;
        let entry_end = entry_start
            .checked_add(padded_length)
            .ok_or_else(|| "raw Git index entry end overflow".to_string())?;
        if entry_end > checksum_start {
            return Err(format!(
                "raw Git index entry {ordinal} padding overlaps its checksum"
            ));
        }
        if index[nul_offset + 1..entry_end]
            .iter()
            .any(|byte| *byte != 0)
        {
            return Err(format!("raw Git index entry {ordinal} has non-NUL padding"));
        }
        entry_order.push((
            index[path_start..nul_offset].to_vec(),
            (flags >> 12) & 0x0003,
            is_symlink,
        ));
        offset = entry_end;
    }

    let mut fsmonitor_extension_present = false;
    while offset < checksum_start {
        let remaining = checksum_start - offset;
        if remaining < 8 {
            return Err("raw Git index has a truncated extension header".to_string());
        }
        let signature = raw_index_range(index, offset, 4, "extension signature")?;
        let payload_length = usize::try_from(raw_index_u32(
            index,
            offset + 4,
            "extension payload length",
        )?)
        .map_err(|_| "raw Git index extension length does not fit this platform".to_string())?;
        let payload_start = offset
            .checked_add(8)
            .ok_or_else(|| "raw Git index extension header overflow".to_string())?;
        let payload_end = payload_start
            .checked_add(payload_length)
            .ok_or_else(|| "raw Git index extension payload overflow".to_string())?;
        if payload_end > checksum_start {
            return Err(format!(
                "raw Git index extension {} overlaps its checksum",
                escape_bytes(signature)
            ));
        }
        if signature == b"link" {
            return Err(
                "raw Git index uses a split-index link extension; the primary index is not a complete authority"
                    .to_string(),
            );
        }
        if signature == b"sdir" {
            return Err(
                "raw Git index uses the required sparse-index sdir extension; collapsed entries are not a complete authority"
                    .to_string(),
            );
        }
        if signature[0].is_ascii_lowercase() {
            return Err(format!(
                "raw Git index uses unsupported required extension {}",
                escape_bytes(signature)
            ));
        }
        if !signature[0].is_ascii_uppercase() {
            return Err(format!(
                "raw Git index has an invalid extension signature {}",
                escape_bytes(signature)
            ));
        }
        if signature == b"FSMN" {
            fsmonitor_extension_present = true;
        }
        offset = payload_end;
    }
    if offset != checksum_start {
        return Err("raw Git index extension layout is ambiguous".to_string());
    }
    for (ordinal, (path, _, is_symlink)) in entry_order.iter().enumerate() {
        validate_git_index_path_bytes(path, *is_symlink).map_err(|error| {
            format!("raw Git index entry {ordinal} has a noncanonical path: {error}")
        })?;
    }
    if let Some((ordinal, _pair)) = entry_order
        .windows(2)
        .enumerate()
        .find(|(_, pair)| (pair[0].0.as_slice(), pair[0].1) >= (pair[1].0.as_slice(), pair[1].1))
    {
        return Err(format!(
            "raw Git index entries {} and {} are not strictly sorted by path and stage",
            ordinal,
            ordinal + 1
        ));
    }
    authenticate_raw_index_checksum(index, checksum_start, object_id_bytes)?;
    Ok(RawIndexInspection {
        entry_count,
        fsmonitor_extension_present,
    })
}

#[cfg(test)]
fn raw_index_contains_fsmonitor_extension(
    index: &[u8],
    object_id_bytes: usize,
) -> Result<bool, String> {
    inspect_raw_git_index(index, object_id_bytes)
        .map(|inspection| inspection.fsmonitor_extension_present)
}

fn parse_stage_zero_entries(index: &[u8]) -> Result<Vec<StageZeroEntry>, String> {
    if !index.is_empty() && !index.ends_with(&[0]) {
        return Err("git index inventory is not NUL-terminated".to_string());
    }
    let mut entries = Vec::new();
    for record in nul_records(index) {
        let tab = record
            .iter()
            .position(|byte| *byte == b'\t')
            .ok_or_else(|| format!("malformed git index record: {}", escape_bytes(record)))?;
        let header = &record[..tab];
        let path = &record[tab + 1..];
        let mut fields = header.split(|byte| *byte == b' ');
        let mode = fields.next().unwrap_or_default();
        let object = fields.next().unwrap_or_default();
        let stage = fields.next().unwrap_or_default();
        if fields.next().is_some() || mode.is_empty() || object.is_empty() || stage.is_empty() {
            return Err(format!(
                "malformed git index record: {}",
                escape_bytes(record)
            ));
        }
        if !matches!(stage, b"0" | b"1" | b"2" | b"3") {
            return Err(format!(
                "malformed git index stage: {}",
                escape_bytes(record)
            ));
        }
        let is_symlink = mode == b"120000";
        validate_git_index_path_bytes(path, is_symlink)?;
        let mode =
            std::str::from_utf8(mode).map_err(|_| "git index mode is not UTF-8".to_string())?;
        let expected_oid = std::str::from_utf8(object)
            .map_err(|_| "git index object ID is not UTF-8".to_string())?;
        if !valid_object_id(expected_oid) {
            return Err("git index object ID is malformed".to_string());
        }
        let path = git_path_from_bytes(path)?;
        validate_relative_git_path(&path)?;
        if stage == b"0" {
            entries.push(StageZeroEntry {
                mode: mode.to_string(),
                expected_oid: expected_oid.to_string(),
                path,
            });
        }
    }
    entries.sort_by(|left, right| left.path.cmp(&right.path));
    if let Some(duplicate) = entries.windows(2).find(|pair| pair[0].path == pair[1].path) {
        return Err(format!(
            "git index inventory contains duplicate stage-0 path {}",
            escaped_path(&duplicate[0].path)
        ));
    }
    Ok(entries)
}

#[cfg(test)]
fn parse_stage_zero_gitlinks(index: &[u8]) -> Result<Vec<Gitlink>, String> {
    let mut gitlinks: Vec<_> = parse_stage_zero_entries(index)?
        .into_iter()
        .filter(|entry| entry.mode == "160000")
        .map(|entry| Gitlink {
            expected_head: entry.expected_oid,
            path: entry.path,
        })
        .collect();
    gitlinks.sort_by(|left, right| left.path.cmp(&right.path));
    Ok(gitlinks)
}

fn hidden_assume_or_skip_record(record: &[u8]) -> bool {
    record
        .first()
        .is_some_and(|tag| *tag == b'S' || tag.is_ascii_lowercase())
}

fn nul_records(bytes: &[u8]) -> impl Iterator<Item = &[u8]> {
    bytes
        .split(|byte| *byte == 0)
        .filter(|record| !record.is_empty())
}

fn git_path_output(repository: &Path, args: &[&str], context: &str) -> Result<PathBuf, String> {
    let bytes = git_bytes(repository, args)?;
    let path = bytes
        .strip_suffix(b"\n")
        .ok_or_else(|| format!("git {context} output is not newline-terminated"))?;
    git_path_from_bytes(path)
}

fn git_head(repository: &Path) -> Result<Option<String>, String> {
    let output = observation_git_command(repository, &["rev-parse", "--verify", "--quiet", "HEAD"])
        .output()
        .map_err(|error| format!("git rev-parse HEAD failed to spawn: {error}"))?;
    if output.status.success() {
        return parse_object_ids(&output.stdout, 1, "HEAD").map(|mut object_ids| object_ids.pop());
    }
    if output.status.code() == Some(1) && output.stdout.is_empty() {
        return Ok(None);
    }
    Err(format!(
        "git rev-parse HEAD in {} failed: {}",
        escaped_path(repository),
        escaped_git_stderr(&output.stderr)
    ))
}

fn git_bytes(repository: &Path, args: &[&str]) -> Result<Vec<u8>, String> {
    let output = observation_git_command(repository, args)
        .output()
        .map_err(|error| format!("git {args:?} failed to spawn: {error}"))?;
    if output.status.success() {
        Ok(output.stdout)
    } else {
        Err(format!(
            "git {args:?} in {} failed: {}",
            escaped_path(repository),
            escaped_git_stderr(&output.stderr)
        ))
    }
}

fn require_no_executable_git_configuration(repository: &Path, scope: &Path) -> Result<(), String> {
    let worktree_config_enabled = worktree_git_configuration_enabled(repository)?;
    for config_scope in ["--local", "--worktree"] {
        if config_scope == "--worktree" && !worktree_config_enabled {
            continue;
        }
        let output = sanitized_git_command(
            repository,
            &[
                "config",
                config_scope,
                "--null",
                "--name-only",
                "--no-includes",
                "--list",
            ],
        )
        .output()
        .map_err(|error| format!("git config inspection failed to spawn: {error}"))?;
        if !output.status.success() {
            return Err(format!(
                "cannot inspect {config_scope} Git configuration in {}: {}",
                escaped_path(repository),
                escaped_git_stderr(&output.stderr)
            ));
        }
        if !output.stdout.is_empty() && !output.stdout.ends_with(&[0]) {
            return Err(format!(
                "{config_scope} Git configuration inventory in {} is not NUL-terminated",
                escaped_path(repository)
            ));
        }
        for key in nul_records(&output.stdout) {
            let key = std::str::from_utf8(key).map_err(|_| {
                format!(
                    "{config_scope} Git configuration key in {} is not UTF-8",
                    escaped_path(repository)
                )
            })?;
            if executable_git_configuration_key(key) {
                return Err(format!(
                    "{}: {config_scope} Git configuration authority {key:?} may execute code or redirect history; remove it before bootstrap verification",
                    escaped_path(scope)
                ));
            }
        }
    }
    Ok(())
}

fn worktree_git_configuration_enabled(repository: &Path) -> Result<bool, String> {
    let output = sanitized_git_command(
        repository,
        &[
            "config",
            "--local",
            "--no-includes",
            "--type=bool",
            "--get-all",
            "extensions.worktreeConfig",
        ],
    )
    .output()
    .map_err(|error| format!("git worktree-config inspection failed to spawn: {error}"))?;
    parse_worktree_config_query(output.status.code(), &output.stdout, &output.stderr).map_err(
        |detail| {
            format!(
                "cannot determine worktree Git configuration authority in {}: {detail}",
                escaped_path(repository)
            )
        },
    )
}

fn parse_worktree_config_query(
    exit_code: Option<i32>,
    stdout: &[u8],
    stderr: &[u8],
) -> Result<bool, String> {
    match (exit_code, stdout, stderr) {
        (Some(0), b"true\n", b"") => Ok(true),
        (Some(0), b"false\n", b"") => Ok(false),
        (Some(0), _, _) => Err("expected exactly one canonical boolean value".to_string()),
        (Some(1), b"", b"") => Ok(false),
        (Some(code), _, _) => Err(format!(
            "git config exited with status {code}: {}",
            escaped_git_stderr(stderr)
        )),
        (None, _, _) => Err("git config terminated without an exit status".to_string()),
    }
}

fn executable_git_configuration_key(key: &str) -> bool {
    let key = key.to_ascii_lowercase();
    key.starts_with("include.")
        || key.starts_with("includeif.")
        || key.starts_with("filter.")
        || matches!(
            key.as_str(),
            "core.fsmonitor"
                | "core.hookspath"
                | "core.sshcommand"
                | "core.gitproxy"
                | "core.alternaterefscommand"
                | "core.askpass"
                | "credential.helper"
                | "extensions.refstorage"
                | "fetch.bundleuri"
                | "fetch.uriprotocols"
                | "gc.recentobjectshook"
                | "protocol.allow"
        )
        || (key.starts_with("protocol.") && key.ends_with(".allow"))
        || (key.starts_with("credential.") && key.ends_with(".helper"))
        || (key.starts_with("remote.")
            && (key.ends_with(".uploadpack")
                || key.ends_with(".receivepack")
                || key.ends_with(".vcs")))
        || (key.starts_with("submodule.") && key.ends_with(".update"))
        || (key.starts_with("url.")
            && (key.ends_with(".insteadof") || key.ends_with(".pushinsteadof")))
}

#[cfg(windows)]
const NULL_GIT_PATH: &str = "NUL";
#[cfg(not(windows))]
const NULL_GIT_PATH: &str = "/dev/null";

fn sanitized_git_prefix(repository: &Path) -> Command {
    let mut command = Command::new("git");
    command
        .arg("-c")
        .arg(format!("core.hooksPath={NULL_GIT_PATH}"))
        .arg("-c")
        .arg(format!("core.attributesFile={NULL_GIT_PATH}"))
        .arg("-c")
        .arg("credential.helper=")
        .arg("-c")
        .arg("protocol.allow=never")
        .arg("-c")
        .arg("protocol.file.allow=always")
        .arg("-c")
        .arg("protocol.https.allow=always")
        .arg("-c")
        .arg("protocol.ssh.allow=always")
        .arg("-C")
        .arg(repository)
        .env("GIT_NO_REPLACE_OBJECTS", "1")
        .env("GIT_NO_LAZY_FETCH", "1")
        .env("GIT_CONFIG_NOSYSTEM", "1")
        .env("GIT_ATTR_NOSYSTEM", "1")
        .env("GIT_CONFIG_GLOBAL", NULL_GIT_PATH)
        .env("GIT_TERMINAL_PROMPT", "0");
    for inherited in [
        "GIT_DIR",
        "GIT_WORK_TREE",
        "GIT_INDEX_FILE",
        "GIT_COMMON_DIR",
        "GIT_OBJECT_DIRECTORY",
        "GIT_ALTERNATE_OBJECT_DIRECTORIES",
        "GIT_NAMESPACE",
        "GIT_PREFIX",
        "GIT_REPLACE_REF_BASE",
        "GIT_SHALLOW_FILE",
        "GIT_CEILING_DIRECTORIES",
        "GIT_DISCOVERY_ACROSS_FILESYSTEM",
        "GIT_CONFIG",
        "GIT_CONFIG_COUNT",
        "GIT_CONFIG_PARAMETERS",
        "GIT_CONFIG_GLOBAL",
        "GIT_CONFIG_SYSTEM",
        "GIT_CONFIG_NOSYSTEM",
        "GIT_ATTR_NOSYSTEM",
        "GIT_ATTR_SOURCE",
        "GIT_TEMPLATE_DIR",
        "GIT_DEFAULT_HASH",
        "GIT_DEFAULT_REF_FORMAT",
        "GIT_REFERENCE_BACKEND",
        "GIT_EXEC_PATH",
        "GIT_EXTERNAL_DIFF",
        "GIT_ASKPASS",
        "SSH_ASKPASS",
        "GIT_SSH",
        "GIT_SSH_COMMAND",
        "GIT_PROXY_COMMAND",
        "GIT_ALLOW_PROTOCOL",
        "GIT_PROTOCOL_FROM_USER",
        "GIT_TERMINAL_PROMPT",
        "GIT_QUARANTINE_PATH",
        "GIT_OPTIONAL_LOCKS",
        "GIT_REDIRECT_STDIN",
        "GIT_REDIRECT_STDOUT",
        "GIT_REDIRECT_STDERR",
        "GIT_LITERAL_PATHSPECS",
        "GIT_GLOB_PATHSPECS",
        "GIT_NOGLOB_PATHSPECS",
        "GIT_ICASE_PATHSPECS",
        "GIT_TRACE",
        "GIT_TRACE_CURL",
        "GIT_TRACE_CURL_NO_DATA",
        "GIT_TRACE_FSMONITOR",
        "GIT_TRACE_PACK_ACCESS",
        "GIT_TRACE_PACKET",
        "GIT_TRACE_PACKFILE",
        "GIT_TRACE_PERFORMANCE",
        "GIT_TRACE_REFS",
        "GIT_TRACE_SETUP",
        "GIT_TRACE_SHALLOW",
        "GIT_TRACE2",
        "GIT_TRACE2_EVENT",
        "GIT_TRACE2_PERF",
    ] {
        command.env_remove(inherited);
    }
    command
        .env("GIT_CONFIG_NOSYSTEM", "1")
        .env("GIT_ATTR_NOSYSTEM", "1")
        .env("GIT_CONFIG_GLOBAL", NULL_GIT_PATH)
        .env("GIT_TERMINAL_PROMPT", "0");
    command
}

pub(crate) fn sanitized_git_command(repository: &Path, args: &[&str]) -> Command {
    let mut command = sanitized_git_prefix(repository);
    command.args(args);
    command
}

fn observation_git_command(repository: &Path, args: &[&str]) -> Command {
    let mut command = sanitized_git_prefix(repository);
    command
        .arg("-c")
        .arg("core.fsmonitor=false")
        .arg("-c")
        .arg("core.untrackedCache=false")
        .args(args)
        .env("GIT_OPTIONAL_LOCKS", "0");
    command
}

fn escaped_git_stderr(stderr: &[u8]) -> String {
    let stderr = stderr.strip_suffix(b"\n").unwrap_or(stderr);
    if stderr.is_empty() {
        "<no stderr>".to_string()
    } else {
        escape_bytes(stderr)
    }
}

fn escape_bytes(bytes: &[u8]) -> String {
    let mut escaped = String::new();
    for &byte in bytes {
        match byte {
            b' '..=b'~' if byte != b'\\' => escaped.push(char::from(byte)),
            b'\\' => escaped.push_str("\\\\"),
            b'\n' => escaped.push_str("\\n"),
            b'\r' => escaped.push_str("\\r"),
            b'\t' => escaped.push_str("\\t"),
            0 => escaped.push_str("\\0"),
            _ => {
                use std::fmt::Write as _;
                let _ = write!(escaped, "\\x{byte:02x}");
            }
        }
    }
    escaped
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::OsStr;

    fn command_env<'a>(command: &'a Command, key: &str) -> Option<Option<&'a OsStr>> {
        command
            .get_envs()
            .find_map(|(name, value)| (name == OsStr::new(key)).then_some(value))
    }

    fn hex(bytes: &[u8]) -> String {
        use std::fmt::Write as _;

        let mut encoded = String::with_capacity(bytes.len() * 2);
        for byte in bytes {
            write!(encoded, "{byte:02x}").expect("writing to String cannot fail");
        }
        encoded
    }

    fn raw_index_fixture(
        version: u32,
        object_id_bytes: usize,
        entries: &[(&[u8], u16, Option<u16>)],
        extensions: &[(&[u8; 4], &[u8])],
    ) -> Vec<u8> {
        let mut index = b"DIRC".to_vec();
        index.extend_from_slice(&version.to_be_bytes());
        index.extend_from_slice(
            &u32::try_from(entries.len())
                .expect("fixture entry count")
                .to_be_bytes(),
        );
        for (path, stage, extended) in entries {
            let entry_start = index.len();
            index.extend_from_slice(&[0; 40]);
            index[entry_start + 24..entry_start + 28].copy_from_slice(&0o100_644_u32.to_be_bytes());
            index.extend(std::iter::repeat_n(0x5a, object_id_bytes));
            let encoded_path_length = u16::try_from(path.len().min(0x0fff))
                .expect("fixture path length fits the index name mask");
            let mut flags = encoded_path_length | ((*stage & 0x0003) << 12);
            if extended.is_some() {
                flags |= 0x4000;
            }
            index.extend_from_slice(&flags.to_be_bytes());
            if let Some(extended) = extended {
                index.extend_from_slice(&extended.to_be_bytes());
            }
            index.extend_from_slice(path);
            index.push(0);
            while (index.len() - entry_start) % 8 != 0 {
                index.push(0);
            }
        }
        for (signature, payload) in extensions {
            index.extend_from_slice(*signature);
            index.extend_from_slice(
                &u32::try_from(payload.len())
                    .expect("fixture extension length")
                    .to_be_bytes(),
            );
            index.extend_from_slice(payload);
        }
        match object_id_bytes {
            20 => {
                let digest = sha1_digest(&index).expect("fixture SHA-1");
                index.extend_from_slice(&digest);
            }
            32 => {
                let digest = sha256_digest(&index).expect("fixture SHA-256");
                index.extend_from_slice(&digest);
            }
            _ => index.extend(std::iter::repeat_n(0xa5, object_id_bytes)),
        }
        index
    }

    fn set_first_raw_index_mode(index: &mut Vec<u8>, object_id_bytes: usize, mode: u32) {
        let checksum_start = index
            .len()
            .checked_sub(object_id_bytes)
            .expect("fixture checksum width");
        index[36..40].copy_from_slice(&mode.to_be_bytes());
        index.truncate(checksum_start);
        match object_id_bytes {
            20 => {
                let digest = sha1_digest(index).expect("resigned fixture SHA-1");
                index.extend_from_slice(&digest);
            }
            32 => {
                let digest = sha256_digest(index).expect("resigned fixture SHA-256");
                index.extend_from_slice(&digest);
            }
            _ => panic!("unsupported fixture checksum width"),
        }
    }

    #[test]
    fn internal_sha1_and_sha256_match_known_vectors() {
        let long = b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
        assert_eq!(
            hex(&sha1_digest(b"").expect("SHA-1 empty vector")),
            "da39a3ee5e6b4b0d3255bfef95601890afd80709"
        );
        assert_eq!(
            hex(&sha1_digest(b"abc").expect("SHA-1 abc vector")),
            "a9993e364706816aba3e25717850c26c9cd0d89d"
        );
        assert_eq!(
            hex(&sha1_digest(long).expect("SHA-1 two-block padding vector")),
            "84983e441c3bd26ebaae4aa1f95129e5e54670f1"
        );
        assert_eq!(
            hex(&sha256_digest(b"").expect("SHA-256 empty vector")),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
        assert_eq!(
            hex(&sha256_digest(b"abc").expect("SHA-256 abc vector")),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        assert_eq!(
            hex(&sha256_digest(long).expect("SHA-256 two-block padding vector")),
            "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
        );
    }

    #[test]
    fn stage_parser_is_nul_safe_for_control_characters_in_paths() {
        let index = b"160000 0123456789abcdef0123456789abcdef01234567 0\tnested dir/with\ttab\nand-newline\0\
100644 89abcdef0123456789abcdef0123456789abcdef 0\tordinary.rs\0";
        assert_eq!(
            parse_stage_zero_gitlinks(index),
            Ok(vec![Gitlink {
                expected_head: "0123456789abcdef0123456789abcdef01234567".to_string(),
                path: PathBuf::from("nested dir/with\ttab\nand-newline"),
            }])
        );
    }

    #[test]
    fn git_index_paths_are_validated_as_canonical_bytes() {
        for invalid in [
            b"".as_slice(),
            b"/rooted",
            b"trailing/",
            b"double//slash",
            b".",
            b"..",
            b"dir/./file",
            b"dir/../file",
            b".git/config",
            b".GIT/config",
            b"dir/.git/config",
        ] {
            assert!(
                validate_git_index_path_bytes(invalid, false).is_err(),
                "{}",
                escape_bytes(invalid)
            );
        }
        assert_eq!(
            validate_git_index_path_bytes(b"ordinary/path", false),
            Ok(())
        );
        assert_eq!(
            validate_git_index_path_bytes(b".gitlike/path", false),
            Ok(())
        );

        let raw = raw_index_fixture(2, 20, &[(b"dir/../escape", 0, None)], &[]);
        let error = inspect_raw_git_index(&raw, 20)
            .expect_err("the raw parser must apply byte-level path validation");
        assert!(error.contains("noncanonical path"), "{error}");

        let staged = b"100644 0123456789abcdef0123456789abcdef01234567 0\tdir/../escape\0";
        let error = parse_stage_zero_entries(staged)
            .expect_err("the staged inventory must reject bytes before PathBuf normalization");
        assert!(error.contains("dot traversal component"), "{error}");
    }

    #[test]
    fn git_hfs_ignored_unicode_scalars_are_always_inadmissible() {
        let ignored_scalars = [
            [0xe2, 0x80, 0x8c],
            [0xe2, 0x80, 0x8d],
            [0xe2, 0x80, 0x8e],
            [0xe2, 0x80, 0x8f],
            [0xe2, 0x80, 0xaa],
            [0xe2, 0x80, 0xab],
            [0xe2, 0x80, 0xac],
            [0xe2, 0x80, 0xad],
            [0xe2, 0x80, 0xae],
            [0xe2, 0x81, 0xaa],
            [0xe2, 0x81, 0xab],
            [0xe2, 0x81, 0xac],
            [0xe2, 0x81, 0xad],
            [0xe2, 0x81, 0xae],
            [0xe2, 0x81, 0xaf],
            [0xef, 0xbb, 0xbf],
        ];
        for scalar in ignored_scalars {
            let mut path = b"ordinary".to_vec();
            path.extend_from_slice(&scalar);
            path.extend_from_slice(b"/path");
            let error = validate_git_index_path_bytes(&path, false)
                .expect_err("Git/HFS-ignored Unicode must refuse on every platform");
            assert!(error.contains("Git/HFS canonicalization"), "{error}");
        }
    }

    #[test]
    fn gitmodules_is_forbidden_for_symlinks_but_valid_for_regular_entries() {
        assert_eq!(validate_git_index_path_bytes(b".gitmodules", false), Ok(()));
        assert!(validate_git_index_path_bytes(b".GITMODULES", true).is_err());

        let regular_staged = b"100644 0123456789abcdef0123456789abcdef01234567 0\t.gitmodules\0";
        assert!(parse_stage_zero_entries(regular_staged).is_ok());
        let symlink_staged = b"120000 0123456789abcdef0123456789abcdef01234567 0\t.GITMODULES\0";
        assert!(parse_stage_zero_entries(symlink_staged).is_err());

        let regular_raw = raw_index_fixture(2, 20, &[(b".gitmodules", 0, None)], &[]);
        assert!(inspect_raw_git_index(&regular_raw, 20).is_ok());
        let mut symlink_raw = raw_index_fixture(2, 20, &[(b".GITMODULES", 0, None)], &[]);
        set_first_raw_index_mode(&mut symlink_raw, 20, 0o120_000);
        let error = inspect_raw_git_index(&symlink_raw, 20)
            .expect_err("raw ce_mode must activate symlink .gitmodules refusal");
        assert!(error.contains(".gitmodules"), "{error}");
    }

    #[cfg(unix)]
    #[test]
    fn unix_git_index_paths_retain_legal_backslashes() {
        assert_eq!(
            validate_git_index_path_bytes(b"dir\\literal", false),
            Ok(())
        );
    }

    #[cfg(windows)]
    #[test]
    fn windows_git_index_paths_reject_platform_aliases_and_prefixes() {
        for invalid in [
            b"dir\\file".as_slice(),
            b"C:/rooted",
            b"//server/share",
            b".GIT/config",
            b".git./config",
            b".git /config",
            b"git~1/config",
            b"file:stream",
        ] {
            assert!(
                validate_git_index_path_bytes(invalid, false).is_err(),
                "{}",
                escape_bytes(invalid)
            );
        }
        for alias in [
            b"gitmod~1".as_slice(),
            b"GITMOD~2",
            b"gitmod~3",
            b"gitmod~4",
            b"~1234567",
            b"g~123456",
            b"gi~12345",
            b"gi7~1234",
            b"gi7e~123",
            b"gi7eb~12",
            b"gi7eba~1",
            b"GI7EBA~9",
        ] {
            assert_eq!(validate_git_index_path_bytes(alias, false), Ok(()));
            assert!(validate_git_index_path_bytes(alias, true).is_err());
        }
        for non_alias in [b"gi7eba~0".as_slice(), b"gi7ebx~1", b"gi7e~12x"] {
            assert_eq!(validate_git_index_path_bytes(non_alias, true), Ok(()));
        }
    }

    #[test]
    fn hidden_index_flags_are_always_dirty() {
        assert!(!hidden_assume_or_skip_record(b"H ordinary.rs"));
        assert!(hidden_assume_or_skip_record(b"S skipped.rs"));
        assert!(hidden_assume_or_skip_record(b"h assumed.rs"));
        assert!(hidden_assume_or_skip_record(b"m lowercase-modified.rs"));
    }

    #[test]
    fn raw_index_v2_v3_and_both_hash_widths_are_unambiguous() {
        assert_eq!(
            inspect_optional_raw_git_index(Path::new("does-not-exist"), None),
            Ok(None)
        );
        let v2_sha1 =
            raw_index_fixture(2, 20, &[(b"a.rs", 0, None), (b"conflict.rs", 1, None)], &[]);
        assert_eq!(
            inspect_raw_git_index(&v2_sha1, 20).map(|inspection| inspection.entry_count),
            Ok(2)
        );
        assert_eq!(
            raw_index_contains_fsmonitor_extension(&v2_sha1, 20),
            Ok(false)
        );

        let v3_sha256 = raw_index_fixture(
            3,
            32,
            &[
                (b"intent.rs", 0, Some(0x2000)),
                (b"skip.rs", 0, Some(0x4000)),
            ],
            &[(b"TREE", b"optional extension")],
        );
        assert_eq!(
            raw_index_contains_fsmonitor_extension(&v3_sha256, 32),
            Ok(false)
        );

        let empty_sha1 = raw_index_fixture(2, 20, &[], &[]);
        let empty_sha256 = raw_index_fixture(3, 32, &[], &[]);
        assert_eq!(
            raw_index_contains_fsmonitor_extension(&empty_sha1, 20),
            Ok(false)
        );
        assert_eq!(
            raw_index_contains_fsmonitor_extension(&empty_sha256, 32),
            Ok(false)
        );
    }

    #[test]
    fn raw_index_presence_and_entry_count_match_the_complete_inventory() {
        assert_eq!(
            require_raw_index_for_inventory(Path::new("."), b"", None),
            Ok(())
        );

        let missing_error = require_raw_index_for_inventory(
            Path::new("child"),
            b"100644 object 0\ttracked.rs\0",
            None,
        )
        .expect_err("a nonempty inventory requires a sealed raw index");
        assert!(
            missing_error.contains("reported 1 tracked index entries"),
            "{missing_error}"
        );

        let inspection = RawIndexInspection {
            entry_count: 2,
            fsmonitor_extension_present: false,
        };
        let mismatch_error = require_raw_index_for_inventory(
            Path::new("child"),
            b"100644 object 0\ttracked.rs\0",
            Some(&inspection),
        )
        .expect_err("raw and complete inventory counts must agree");
        assert!(
            mismatch_error.contains("entry count 2 disagrees"),
            "{mismatch_error}"
        );

        let conflicts = b"100644 object 1\tconflict.rs\0\
100644 object 2\tconflict.rs\0\
100644 object 3\tconflict.rs\0";
        let conflict_inspection = RawIndexInspection {
            entry_count: 3,
            fsmonitor_extension_present: false,
        };
        assert_eq!(
            require_raw_index_for_inventory(
                Path::new("child"),
                conflicts,
                Some(&conflict_inspection)
            ),
            Ok(())
        );
    }

    #[cfg(unix)]
    #[test]
    fn primary_index_must_be_a_direct_child_of_the_canonical_git_dir() {
        let root = std::env::current_dir().expect("current directory is absolute");
        let git_dir = root.join("fixture.git");
        assert_eq!(
            require_primary_index_path(&git_dir, &git_dir.join("index")),
            Ok(())
        );
        assert!(require_primary_index_path(&git_dir, Path::new("index")).is_err());
        assert!(require_primary_index_path(&git_dir, &git_dir.join("INDEX")).is_err());
        assert!(
            require_primary_index_path(&git_dir, &root.join("other.git").join("index")).is_err()
        );
    }

    #[cfg(not(unix))]
    #[test]
    fn primary_index_parent_is_canonicalized_on_non_unix() {
        let git_dir = std::env::current_dir()
            .expect("current directory")
            .canonicalize()
            .expect("canonical current directory");
        assert_eq!(
            require_primary_index_path(&git_dir, &git_dir.join("index")),
            Ok(())
        );
        assert!(require_primary_index_path(&git_dir, Path::new("index")).is_err());
        assert!(require_primary_index_path(&git_dir, &git_dir.join("INDEX")).is_err());
    }

    #[test]
    fn raw_index_fsmn_extension_is_always_inadmissible() {
        let fsmn_v1 = raw_index_fixture(
            2,
            20,
            &[(b"tracked.rs", 0, None)],
            &[(b"FSMN", b"\0\0\0\x01opaque-v1-ewah")],
        );
        let fsmn_v2 = raw_index_fixture(
            3,
            32,
            &[(b"tracked.rs", 0, None)],
            &[(b"FSMN", b"\0\0\0\x02opaque-token\0opaque-v2-ewah")],
        );
        assert_eq!(
            raw_index_contains_fsmonitor_extension(&fsmn_v1, 20),
            Ok(true)
        );
        assert_eq!(
            raw_index_contains_fsmonitor_extension(&fsmn_v2, 32),
            Ok(true)
        );
    }

    #[test]
    fn raw_index_never_treats_the_in_memory_fsmonitor_bit_as_serialized() {
        let reserved_extended_bit =
            raw_index_fixture(3, 20, &[(b"tracked.rs", 0, Some(0x0020))], &[]);
        let error = raw_index_contains_fsmonitor_extension(&reserved_extended_bit, 20)
            .expect_err("CE_FSMONITOR_VALID's in-memory bit position is not an on-disk flag");
        assert!(
            error.contains("unsupported extended flags 0x0020"),
            "{error}"
        );

        let v2_extended = raw_index_fixture(2, 20, &[(b"tracked.rs", 0, Some(0x2000))], &[]);
        let error = raw_index_contains_fsmonitor_extension(&v2_extended, 20)
            .expect_err("v2 entries cannot carry the v3 extended word");
        assert!(
            error.contains("v2 entry 0 uses v3 extended flags"),
            "{error}"
        );
    }

    #[test]
    fn raw_index_refuses_split_sparse_and_unknown_required_extensions() {
        for (signature, expected) in [
            (b"link", "split-index link extension"),
            (b"sdir", "sparse-index sdir extension"),
            (b"abcd", "unsupported required extension abcd"),
        ] {
            let index = raw_index_fixture(2, 20, &[], &[(signature, b"opaque")]);
            let error = raw_index_contains_fsmonitor_extension(&index, 20)
                .expect_err("required extension semantics must refuse closed");
            assert!(error.contains(expected), "{error}");
        }

        let invalid_signature = raw_index_fixture(2, 20, &[], &[(b"1bad", b"opaque")]);
        let error = raw_index_contains_fsmonitor_extension(&invalid_signature, 20)
            .expect_err("only uppercase-leading optional extensions may be skipped");
        assert!(
            error.contains("invalid extension signature 1bad"),
            "{error}"
        );
    }

    #[test]
    fn raw_index_checksum_authentication_rejects_payload_and_digest_corruption() {
        for object_id_bytes in [20, 32] {
            let mut corrupt_payload =
                raw_index_fixture(2, object_id_bytes, &[(b"tracked.rs", 0, None)], &[]);
            corrupt_payload[52] ^= 0x01;
            let error = inspect_raw_git_index(&corrupt_payload, object_id_bytes)
                .expect_err("an altered object ID byte must invalidate the index checksum");
            assert!(error.contains("checksum"), "{error}");

            let mut corrupt_digest =
                raw_index_fixture(2, object_id_bytes, &[(b"tracked.rs", 0, None)], &[]);
            let last = corrupt_digest
                .last_mut()
                .expect("fixture always has a checksum");
            *last ^= 0x01;
            let error = inspect_raw_git_index(&corrupt_digest, object_id_bytes)
                .expect_err("an altered trailing digest must invalidate the index checksum");
            assert!(error.contains("checksum"), "{error}");
        }
    }

    #[test]
    fn raw_index_refuses_truncated_or_ambiguous_bounds() {
        assert!(raw_index_contains_fsmonitor_extension(b"DIRC", 20).is_err());

        let mut bad_padding = raw_index_fixture(2, 20, &[(b"ab", 0, None)], &[]);
        bad_padding[77] = 1;
        let error = raw_index_contains_fsmonitor_extension(&bad_padding, 20)
            .expect_err("entry padding must be entirely NUL");
        assert!(error.contains("non-NUL padding"), "{error}");

        let mut oversized_extension = raw_index_fixture(2, 20, &[], &[(b"FSMN", b"short")]);
        oversized_extension[16..20].copy_from_slice(&u32::MAX.to_be_bytes());
        let error = raw_index_contains_fsmonitor_extension(&oversized_extension, 20)
            .expect_err("extension bounds cannot overlap the checksum");
        assert!(error.contains("overlaps its checksum"), "{error}");

        let mut short_checksum = raw_index_fixture(3, 32, &[], &[]);
        short_checksum.truncate(12 + 31);
        assert!(raw_index_contains_fsmonitor_extension(&short_checksum, 32).is_err());

        let sha1_layout = raw_index_fixture(2, 20, &[(b"tracked.rs", 0, None)], &[]);
        assert!(raw_index_contains_fsmonitor_extension(&sha1_layout, 32).is_err());
        assert!(raw_index_contains_fsmonitor_extension(&sha1_layout, 24).is_err());
    }

    #[test]
    fn raw_index_refuses_v4_and_noncanonical_entry_order() {
        let v4 = raw_index_fixture(4, 20, &[(b"tracked.rs", 0, None)], &[]);
        let error = raw_index_contains_fsmonitor_extension(&v4, 20)
            .expect_err("v4 path compression is outside the admitted parser");
        assert!(error.contains("version 4 is unsupported"), "{error}");

        let unsorted = raw_index_fixture(2, 20, &[(b"z.rs", 0, None), (b"a.rs", 0, None)], &[]);
        let error = raw_index_contains_fsmonitor_extension(&unsorted, 20)
            .expect_err("unsorted entries make the authority noncanonical");
        assert!(error.contains("not strictly sorted"), "{error}");

        let duplicate_stage =
            raw_index_fixture(2, 20, &[(b"same.rs", 2, None), (b"same.rs", 2, None)], &[]);
        assert!(raw_index_contains_fsmonitor_extension(&duplicate_stage, 20).is_err());
    }

    #[test]
    fn diagnostics_escape_nul_and_path_control_bytes() {
        assert_eq!(
            escape_bytes(b"a b\tline\nend\0\xff"),
            "a b\\tline\\nend\\0\\xff"
        );
        assert_eq!(
            escaped_git_stderr(b"fatal:\tbad\xff\n"),
            "fatal:\\tbad\\xff"
        );
        assert_eq!(escaped_git_stderr(b"\n"), "<no stderr>");
    }

    #[test]
    fn stable_policy_kinds_are_present_in_diagnostics() {
        assert_eq!(
            untracked_ignore_policy_finding(Path::new("nested\nrepo"), b"dir/.gitignore\t"),
            "nested\\nrepo: kind=untracked-ignore-policy path=dir/.gitignore\\t detail=only tracked .gitignore files may define project ignore semantics"
        );

        let states = [RawTrackedState {
            expected_mode: "100755".to_string(),
            expected_oid: "expected".to_string(),
            path: PathBuf::from("tool"),
            actual_kind: WorktreeKind::Regular,
            actual_mode: Some("100644".to_string()),
            actual_oid: Some("actual".to_string()),
            link_target: None,
            metadata_identity: None,
        }];
        assert_eq!(
            raw_tracked_findings(Path::new("child"), &states),
            vec![
                "child/tool: kind=raw-tracked-source-mismatch detail=mode expected=100755 actual=100644".to_string(),
                "child/tool: kind=raw-tracked-source-mismatch detail=object expected=expected actual=actual".to_string(),
            ]
        );
    }

    #[test]
    fn raw_mode_comparison_respects_platform_observability() {
        assert_eq!(
            TRACKED_STATUS_ARGS[1],
            if cfg!(unix) {
                "core.fileMode=true"
            } else {
                "core.fileMode=false"
            }
        );
        assert!(TRACKED_STATUS_ARGS.contains(&"--ignore-submodules=none"));
        assert!(STAGED_STATUS_ARGS.contains(&"--ignore-submodules=none"));
        assert!(STAGED_STATUS_ARGS.contains(&"--no-ext-diff"));
        assert!(raw_mode_mismatch(
            WorktreeKind::Regular,
            "100755",
            None,
            true
        ));
        assert!(!raw_mode_mismatch(
            WorktreeKind::Regular,
            "100755",
            None,
            false
        ));
        assert!(raw_mode_mismatch(
            WorktreeKind::Symlink,
            "120000",
            None,
            false
        ));
        assert!(raw_mode_mismatch(
            WorktreeKind::Regular,
            "100755",
            Some("100644"),
            false
        ));
        assert!(!raw_mode_mismatch(
            WorktreeKind::Regular,
            "100755",
            Some("100755"),
            true
        ));
    }

    #[test]
    fn intra_pass_boundary_comparison_rejects_forced_status_movement() {
        let before = RepositoryBoundaryObservation::default();
        assert_eq!(
            require_equal_repository_boundaries(Path::new("fixture"), &before, &before),
            Ok(())
        );

        let after = RepositoryBoundaryObservation {
            forced_status: b" M initialized-child\0".to_vec(),
            ..before.clone()
        };
        let error = require_equal_repository_boundaries(Path::new("fixture"), &before, &after)
            .expect_err("a changed forced-visible status must refuse within the same pass");
        assert!(
            error.contains("complete boundary re-observation"),
            "{error}"
        );
    }

    #[test]
    fn complete_pass_helper_finishes_pass_one_before_an_injected_pass_two_failure() {
        let mut observed = Vec::new();
        let error = verify_two_complete_passes(&[1_u8, 2, 3], |item| {
            observed.push(*item);
            if observed.len() == 5 {
                Err("injected".to_string())
            } else {
                Ok(())
            }
        })
        .expect_err("the injected second-pass failure must propagate");
        assert_eq!(observed, vec![1, 2, 3, 1, 2]);
        assert_eq!(
            error,
            "complete constellation verification pass 2 failed: injected"
        );
    }

    #[test]
    fn pinned_status_rejects_malformed_heads_before_filesystem_observation() {
        assert_eq!(
            pinned_repository_worktree_status(Path::new("does-not-exist"), "not-an-object-id"),
            Err("expected repository HEAD is not a canonical Git object ID".to_string())
        );
        assert!(valid_object_id("0123456789abcdef0123456789abcdef01234567"));
        assert!(!valid_object_id("0123456789ABCDEF0123456789ABCDEF01234567"));
    }

    #[test]
    fn ordinary_root_check_rejects_a_non_directory_before_canonicalization() {
        let executable = std::env::current_exe().expect("test executable path");
        let error = ordinary_root_identity(&executable)
            .expect_err("a regular file cannot be accepted as a repository root");
        assert!(error.contains("must be an ordinary directory"), "{error}");
    }

    #[test]
    fn sanitized_git_factory_is_mutation_safe_and_observers_disable_locks() {
        let command = sanitized_git_command(Path::new("/fixture/repository"), &["status"]);
        assert_eq!(
            command_env(&command, "GIT_NO_REPLACE_OBJECTS"),
            Some(Some(OsStr::new("1")))
        );
        assert_eq!(
            command_env(&command, "GIT_NO_LAZY_FETCH"),
            Some(Some(OsStr::new("1")))
        );
        assert_eq!(command_env(&command, "GIT_DIR"), Some(None));
        assert_eq!(command_env(&command, "GIT_WORK_TREE"), Some(None));
        assert_eq!(command_env(&command, "GIT_OPTIONAL_LOCKS"), Some(None));
        assert_eq!(command_env(&command, "GIT_CONFIG"), Some(None));
        assert_eq!(command_env(&command, "GIT_ATTR_SOURCE"), Some(None));
        assert_eq!(command_env(&command, "GIT_DEFAULT_HASH"), Some(None));
        assert_eq!(command_env(&command, "GIT_DEFAULT_REF_FORMAT"), Some(None));
        assert_eq!(command_env(&command, "GIT_REFERENCE_BACKEND"), Some(None));
        assert_eq!(command_env(&command, "GIT_REDIRECT_STDIN"), Some(None));
        assert_eq!(command_env(&command, "GIT_REDIRECT_STDOUT"), Some(None));
        assert_eq!(command_env(&command, "GIT_REDIRECT_STDERR"), Some(None));
        assert_eq!(command_env(&command, "GIT_LITERAL_PATHSPECS"), Some(None));
        assert_eq!(command_env(&command, "GIT_GLOB_PATHSPECS"), Some(None));
        assert_eq!(command_env(&command, "GIT_NOGLOB_PATHSPECS"), Some(None));
        assert_eq!(command_env(&command, "GIT_ICASE_PATHSPECS"), Some(None));
        assert_eq!(command_env(&command, "GIT_TRACE"), Some(None));
        assert_eq!(command_env(&command, "GIT_TRACE2_EVENT"), Some(None));
        assert_eq!(
            command_env(&command, "GIT_CONFIG_GLOBAL"),
            Some(Some(OsStr::new(NULL_GIT_PATH)))
        );
        assert_eq!(
            command_env(&command, "GIT_CONFIG_NOSYSTEM"),
            Some(Some(OsStr::new("1")))
        );
        assert_eq!(
            command.get_args().collect::<Vec<_>>(),
            vec![
                OsStr::new("-c"),
                OsStr::new(if cfg!(windows) {
                    "core.hooksPath=NUL"
                } else {
                    "core.hooksPath=/dev/null"
                }),
                OsStr::new("-c"),
                OsStr::new(if cfg!(windows) {
                    "core.attributesFile=NUL"
                } else {
                    "core.attributesFile=/dev/null"
                }),
                OsStr::new("-c"),
                OsStr::new("credential.helper="),
                OsStr::new("-c"),
                OsStr::new("protocol.allow=never"),
                OsStr::new("-c"),
                OsStr::new("protocol.file.allow=always"),
                OsStr::new("-c"),
                OsStr::new("protocol.https.allow=always"),
                OsStr::new("-c"),
                OsStr::new("protocol.ssh.allow=always"),
                OsStr::new("-C"),
                OsStr::new("/fixture/repository"),
                OsStr::new("status")
            ]
        );

        let observer = observation_git_command(Path::new("/fixture/repository"), &["status"]);
        assert_eq!(
            command_env(&observer, "GIT_OPTIONAL_LOCKS"),
            Some(Some(OsStr::new("0")))
        );
    }

    #[test]
    fn executable_local_git_configuration_is_fail_closed() {
        for key in [
            "include.path",
            "includeIf.onbranch:main.path",
            "filter.canonical.process",
            "core.fsmonitor",
            "core.hooksPath",
            "core.sshCommand",
            "core.askPass",
            "credential.https://example.invalid.helper",
            "extensions.refStorage",
            "fetch.bundleURI",
            "fetch.uriProtocols",
            "gc.recentObjectsHook",
            "remote.origin.uploadpack",
            "remote.origin.vcs",
            "submodule.child.update",
            "protocol.allow",
            "protocol.ext.allow",
            "protocol.https.allow",
            "url.file:///redirect/.insteadOf",
        ] {
            assert!(executable_git_configuration_key(key), "{key}");
        }
        for key in [
            "core.autocrlf",
            "frankensim.bootstrapIncomplete",
            "remote.origin.url",
        ] {
            assert!(!executable_git_configuration_key(key), "{key}");
        }
    }

    #[test]
    fn inactive_worktree_configuration_scope_is_skipped_without_ambiguity() {
        assert_eq!(parse_worktree_config_query(Some(1), b"", b""), Ok(false));
        assert_eq!(
            parse_worktree_config_query(Some(0), b"false\n", b""),
            Ok(false)
        );
        assert_eq!(
            parse_worktree_config_query(Some(0), b"true\n", b""),
            Ok(true)
        );

        for (exit_code, stdout, stderr) in [
            (Some(0), b"".as_slice(), b"".as_slice()),
            (Some(0), b"true\nfalse\n", b""),
            (Some(0), b"yes\n", b""),
            (Some(0), b"true\n", b"warning"),
            (Some(1), b"false\n", b""),
            (Some(2), b"", b"invalid value"),
            (None, b"", b""),
        ] {
            assert!(
                parse_worktree_config_query(exit_code, stdout, stderr).is_err(),
                "exit_code={exit_code:?} stdout={stdout:?} stderr={stderr:?}"
            );
        }
    }

    #[test]
    fn stage_parser_rejects_truncated_and_duplicate_stage_zero_inventory() {
        let record = b"100644 0123456789abcdef0123456789abcdef01234567 0\tduplicate";
        assert_eq!(
            parse_stage_zero_entries(record),
            Err("git index inventory is not NUL-terminated".to_string())
        );

        let mut duplicate = Vec::new();
        duplicate.extend_from_slice(record);
        duplicate.push(0);
        duplicate.extend_from_slice(record);
        duplicate.push(0);
        let duplicate_error = parse_stage_zero_entries(&duplicate)
            .expect_err("duplicate stage-zero paths must refuse closed");
        assert!(duplicate_error.contains("duplicate stage-0 path duplicate"));
    }

    #[cfg(unix)]
    #[test]
    fn unix_stage_parser_preserves_non_utf8_paths_and_owner_execute_mode() {
        use std::os::unix::ffi::OsStrExt as _;

        let mut index = b"100644 0123456789abcdef0123456789abcdef01234567 0\tnon-utf8-".to_vec();
        index.extend_from_slice(b"\xff\0");
        let entries = parse_stage_zero_entries(&index).expect("non-UTF-8 Unix paths are valid");
        assert_eq!(entries[0].path.as_os_str().as_bytes(), b"non-utf8-\xff");
        assert_eq!(escaped_path(&entries[0].path), "non-utf8-\\xff");
        assert_eq!(git_mode_from_unix_permissions(0o010), "100644");
        assert_eq!(git_mode_from_unix_permissions(0o100), "100755");
    }

    #[test]
    fn recursive_observation_has_a_hard_depth_bound() {
        let error = observe_repository(
            Path::new("unused"),
            Path::new("nested"),
            None,
            MAX_NESTED_REPOSITORY_DEPTH + 1,
        )
        .expect_err("over-depth recursion must fail before touching the filesystem");
        assert!(error.contains("nested repository depth exceeds"));
    }
}
