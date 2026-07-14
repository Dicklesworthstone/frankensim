//! Canonical constellation-bootstrap provenance shared by the workspace xtask
//! and the zero-dependency standalone bootstrap binary.

use std::fmt::Write as _;
#[cfg(unix)]
use std::fs::Metadata;
use std::fs::{File, OpenOptions};
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_TEMP_SUFFIX: AtomicU64 = AtomicU64::new(0);

pub(crate) const BOOTSTRAP_PROVENANCE_SCHEMA: &str = "frankensim-constellation-bootstrap-v2";
pub(crate) const BOOTSTRAP_PROVENANCE_IDENTITY_VERSION: u32 = 3;
pub(crate) const BOOTSTRAP_PROVENANCE_IDENTITY_DOMAIN: &str =
    "org.frankensim.xtask.constellation-bootstrap-provenance.v3";

/// The ordinary, single-link file object written and fsynced before the final
/// publication barrier. Portable `Metadata` has no cross-platform object ID,
/// so targets without a safe `std` identity surface refuse publication.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct StagingSeal {
    #[cfg(unix)]
    device: u64,
    #[cfg(unix)]
    inode: u64,
    #[cfg(unix)]
    mode: u32,
    #[cfg(unix)]
    links: u64,
    #[cfg(unix)]
    owner: u32,
    #[cfg(unix)]
    group: u32,
    #[cfg(unix)]
    size: u64,
    #[cfg(unix)]
    modified_seconds: i64,
    #[cfg(unix)]
    modified_nanoseconds: i64,
    #[cfg(unix)]
    changed_seconds: i64,
    #[cfg(unix)]
    changed_nanoseconds: i64,
    #[cfg(windows)]
    volume_serial_number: u32,
    #[cfg(windows)]
    file_index: u64,
    #[cfg(windows)]
    attributes: u32,
    #[cfg(windows)]
    links: u32,
    #[cfg(windows)]
    size: u64,
    #[cfg(windows)]
    created: u64,
    #[cfg(windows)]
    modified: u64,
}

impl StagingSeal {
    fn capture(staging: &File) -> std::io::Result<Self> {
        bootstrap_provenance_support_preflight()?;
        let metadata = staging.metadata()?;
        #[cfg(unix)]
        {
            Self::from_unix_metadata(&metadata)
        }
        #[cfg(windows)]
        {
            Self::from_windows_metadata(&metadata)
        }
        #[cfg(not(any(unix, windows)))]
        {
            let _ = metadata;
            Err(unsupported_staging_identity_error())
        }
    }

    fn verify(&self, staging: &File, staging_path: &Path) -> std::io::Result<()> {
        #[cfg(unix)]
        {
            let pinned = Self::from_unix_metadata(&staging.metadata()?)?;
            let visible = Self::from_unix_metadata(&std::fs::symlink_metadata(staging_path)?)?;
            self.require_matching_authority(pinned, visible, staging_path)
        }
        #[cfg(windows)]
        {
            use std::os::windows::fs::OpenOptionsExt as _;

            // Open the directory entry itself rather than following a reparse
            // point, so a symlink/junction cannot borrow the sealed target's ID.
            const FILE_FLAG_OPEN_REPARSE_POINT: u32 = 0x0020_0000;
            let visible_file = OpenOptions::new()
                .read(true)
                .custom_flags(FILE_FLAG_OPEN_REPARSE_POINT)
                .open(staging_path)?;
            let pinned = Self::from_windows_metadata(&staging.metadata()?)?;
            let visible = Self::from_windows_metadata(&visible_file.metadata()?)?;
            self.require_matching_authority(pinned, visible, staging_path)
        }
        #[cfg(not(any(unix, windows)))]
        {
            let _ = (self, staging, staging_path);
            Err(unsupported_staging_identity_error())
        }
    }

    #[cfg(any(unix, windows))]
    fn require_matching_authority(
        &self,
        pinned: Self,
        visible: Self,
        staging_path: &Path,
    ) -> std::io::Result<()> {
        if pinned == *self && visible == *self {
            Ok(())
        } else {
            Err(std::io::Error::other(format!(
                "the still-open staging file or visible pathname {} no longer matches the sealed ordinary single-link file",
                staging_path.display()
            )))
        }
    }

    #[cfg(unix)]
    fn from_unix_metadata(metadata: &Metadata) -> std::io::Result<Self> {
        use std::os::unix::fs::MetadataExt as _;

        if !metadata.file_type().is_file() {
            return Err(std::io::Error::other(
                "bootstrap provenance staging authority is not an ordinary file",
            ));
        }
        let seal = Self {
            device: metadata.dev(),
            inode: metadata.ino(),
            mode: metadata.mode(),
            links: metadata.nlink(),
            owner: metadata.uid(),
            group: metadata.gid(),
            size: metadata.size(),
            modified_seconds: metadata.mtime(),
            modified_nanoseconds: metadata.mtime_nsec(),
            changed_seconds: metadata.ctime(),
            changed_nanoseconds: metadata.ctime_nsec(),
        };
        if seal.links != 1 {
            return Err(std::io::Error::other(format!(
                "bootstrap provenance staging authority has {} hard links; expected exactly one",
                seal.links
            )));
        }
        Ok(seal)
    }

    #[cfg(windows)]
    fn from_windows_metadata(metadata: &std::fs::Metadata) -> std::io::Result<Self> {
        use std::os::windows::fs::MetadataExt as _;

        const FILE_ATTRIBUTE_REPARSE_POINT: u32 = 0x0000_0400;
        let attributes = metadata.file_attributes();
        if !metadata.file_type().is_file() || attributes & FILE_ATTRIBUTE_REPARSE_POINT != 0 {
            return Err(std::io::Error::other(
                "bootstrap provenance staging authority is not an ordinary reparse-free file",
            ));
        }
        let seal = Self {
            volume_serial_number: metadata.volume_serial_number().ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::Unsupported,
                    "Windows filesystem did not expose a staging volume serial number",
                )
            })?,
            file_index: metadata.file_index().ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::Unsupported,
                    "Windows filesystem did not expose a staging file index",
                )
            })?,
            attributes,
            links: metadata.number_of_links().ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::Unsupported,
                    "Windows filesystem did not expose a staging hard-link count",
                )
            })?,
            size: metadata.file_size(),
            created: metadata.creation_time(),
            modified: metadata.last_write_time(),
        };
        if seal.links != 1 {
            return Err(std::io::Error::other(format!(
                "bootstrap provenance staging authority has {} hard links; expected exactly one",
                seal.links
            )));
        }
        Ok(seal)
    }
}

#[cfg(not(any(unix, windows)))]
fn unsupported_staging_identity_error() -> std::io::Error {
    std::io::Error::new(
        std::io::ErrorKind::Unsupported,
        "bootstrap provenance publication is unsupported on this target: safe std exposes no file-object identity that can prove the visible staging pathname still names the written and fsynced file",
    )
}

pub(crate) fn bootstrap_provenance_support_preflight() -> std::io::Result<()> {
    #[cfg(any(unix, windows))]
    {
        Ok(())
    }
    #[cfg(not(any(unix, windows)))]
    {
        Err(unsupported_staging_identity_error())
    }
}

pub(crate) fn provenance_path_text(path: &Path) -> std::io::Result<&str> {
    path.to_str().ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!(
                "bootstrap provenance destination is not valid UTF-8: {}",
                path.display()
            ),
        )
    })
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct BootstrapProvenanceRow {
    pub(crate) lib: String,
    pub(crate) git_head: String,
    pub(crate) remote: String,
    pub(crate) selected_transport: String,
    pub(crate) transport_used: bool,
    pub(crate) state: String,
}

impl BootstrapProvenanceRow {
    pub(crate) fn new(
        lib: &str,
        git_head: &str,
        remote: &str,
        selected_transport: &str,
        transport_used: bool,
        state: &str,
    ) -> Self {
        Self {
            lib: lib.to_string(),
            git_head: git_head.to_string(),
            remote: remote.to_string(),
            selected_transport: selected_transport.to_string(),
            transport_used,
            state: state.to_string(),
        }
    }
}

fn json_escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            control @ '\u{0000}'..='\u{001f}' => {
                let _ = write!(out, "\\u{:04x}", u32::from(control));
            }
            other => out.push(other),
        }
    }
    out
}

pub(crate) fn render_bootstrap_provenance_row(row: &BootstrapProvenanceRow) -> String {
    format!(
        "{{\"lib\": \"{}\", \"git_head\": \"{}\", \"remote\": \"{}\", \"selected_transport\": \"{}\", \"transport_used\": {}, \"state\": \"{}\"}}",
        json_escape(&row.lib),
        json_escape(&row.git_head),
        json_escape(&row.remote),
        json_escape(&row.selected_transport),
        row.transport_used,
        json_escape(&row.state),
    )
}

pub(crate) fn render_bootstrap_provenance(
    lock_hash: &str,
    dest: &str,
    rows: &[BootstrapProvenanceRow],
) -> String {
    let libraries = rows
        .iter()
        .map(render_bootstrap_provenance_row)
        .collect::<Vec<_>>()
        .join(",\n");
    format!(
        "{{\n\"schema\": \"{BOOTSTRAP_PROVENANCE_SCHEMA}\",\n\"identity_domain\": \"{BOOTSTRAP_PROVENANCE_IDENTITY_DOMAIN}\",\n\"identity_version\": {BOOTSTRAP_PROVENANCE_IDENTITY_VERSION},\n\"lock_hash\": \"{}\",\n\"dest\": \"{}\",\n\"libraries\": [\n{libraries}\n]\n}}\n",
        json_escape(lock_hash),
        json_escape(dest),
    )
}

pub(crate) fn write_bootstrap_provenance<Validate>(
    path: &Path,
    lock_hash: &str,
    dest: &str,
    rows: &[BootstrapProvenanceRow],
    validate: Validate,
) -> std::io::Result<()>
where
    Validate: FnOnce() -> Result<(), String>,
{
    let identity_epoch = (
        BOOTSTRAP_PROVENANCE_IDENTITY_DOMAIN,
        BOOTSTRAP_PROVENANCE_IDENTITY_VERSION,
    );
    let document = render_bootstrap_provenance(lock_hash, dest, rows);
    debug_assert!(document.contains(identity_epoch.0));
    debug_assert!(document.contains(&format!("\"identity_version\": {}", identity_epoch.1)));
    let (temporary, staging_file) = reserve_same_directory_temporary(path)?;
    stage_and_replace(
        staging_file,
        &temporary,
        path,
        &document,
        |staging, document| {
            staging.write_all(document)?;
            staging.sync_all()?;
            StagingSeal::capture(staging)
        },
        validate,
        |staging, staging_path, seal| seal.verify(staging, staging_path),
        |staging, destination| std::fs::rename(staging, destination),
    )?;
    sync_parent_best_effort(path);
    Ok(())
}

fn stage_and_replace<Stage, Seal, WriteStage, Validate, VerifyStage, Rename>(
    mut staging: Stage,
    staging_path: &Path,
    destination: &Path,
    document: &str,
    write_stage: WriteStage,
    validate: Validate,
    verify_stage: VerifyStage,
    rename: Rename,
) -> std::io::Result<()>
where
    WriteStage: FnOnce(&mut Stage, &[u8]) -> std::io::Result<Seal>,
    Validate: FnOnce() -> Result<(), String>,
    VerifyStage: FnOnce(&Stage, &Path, &Seal) -> std::io::Result<()>,
    Rename: FnOnce(&Path, &Path) -> std::io::Result<()>,
{
    let seal = write_stage(&mut staging, document.as_bytes()).map_err(|error| {
        std::io::Error::new(
            error.kind(),
            format!(
                "cannot stage bootstrap provenance at {} for {}: {error}",
                staging_path.display(),
                destination.display()
            ),
        )
    })?;
    validate().map_err(|error| {
        std::io::Error::other(format!(
            "bootstrap provenance publication barrier failed after staging {} for {} (no replacement or cleanup attempted): {error}",
            staging_path.display(),
            destination.display()
        ))
    })?;
    // Keep the original handle live and place this authority check immediately
    // beside the same-directory pathname rename. Safe `std` has no
    // handle-relative rename primitive; this closes the long validation window.
    verify_stage(&staging, staging_path, &seal).map_err(|error| {
        std::io::Error::new(
            error.kind(),
            format!(
                "bootstrap provenance staging authority check failed for {} before replacing {} (no replacement or cleanup attempted): {error}",
                staging_path.display(),
                destination.display()
            ),
        )
    })?;
    rename(staging_path, destination).map_err(|error| {
        std::io::Error::new(
            error.kind(),
            format!(
                "cannot replace bootstrap provenance {} from staging path {}: {error}",
                destination.display(),
                staging_path.display()
            ),
        )
    })
}

fn sync_parent_best_effort(path: &Path) {
    // The file itself is durable before rename. Directory fsync is not portable
    // across every supported filesystem, so it remains an explicit best-effort
    // crash-durability boundary rather than turning a published receipt into a
    // reported write failure after replacement already succeeded.
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    if let Ok(directory) = File::open(parent) {
        let _ = directory.sync_all();
    }
}

fn reserve_same_directory_temporary(path: &Path) -> std::io::Result<(PathBuf, File)> {
    let file_name = path.file_name().ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "bootstrap provenance path has no file name",
        )
    })?;
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty());
    let parent = parent.unwrap_or_else(|| Path::new("."));
    for _ in 0..128 {
        let suffix = NEXT_TEMP_SUFFIX.fetch_add(1, Ordering::Relaxed);
        let mut temporary_name = file_name.to_os_string();
        temporary_name.push(format!(".tmp.{}.{suffix}", std::process::id()));
        let temporary = parent.join(temporary_name);
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
        {
            Ok(file) => return Ok((temporary, file)),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(error),
        }
    }
    Err(std::io::Error::new(
        std::io::ErrorKind::AlreadyExists,
        "could not reserve a unique bootstrap provenance staging path",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_escape_preserves_utf8_and_only_escapes_json_c0() {
        let value = "\"\\/\n\r\t\u{0000}\u{0008}\u{000c}\u{001f}\u{007f}\u{0085}é";
        assert_eq!(
            json_escape(value),
            "\\\"\\\\/\\n\\r\\t\\u0000\\u0008\\u000c\\u001f\u{007f}\u{0085}é"
        );
    }

    #[test]
    fn v3_document_shape_is_exact_and_row_order_is_retained() {
        let rows = vec![
            BootstrapProvenanceRow::new("alpha", "a", "upstream-a", "mirror-a", true, "cloned"),
            BootstrapProvenanceRow::new("beta", "b", "upstream-b", "upstream-b", false, "verified"),
        ];
        assert_eq!(
            render_bootstrap_provenance("lock", "/dest", &rows),
            "{\n\"schema\": \"frankensim-constellation-bootstrap-v2\",\n\"identity_domain\": \"org.frankensim.xtask.constellation-bootstrap-provenance.v3\",\n\"identity_version\": 3,\n\"lock_hash\": \"lock\",\n\"dest\": \"/dest\",\n\"libraries\": [\n{\"lib\": \"alpha\", \"git_head\": \"a\", \"remote\": \"upstream-a\", \"selected_transport\": \"mirror-a\", \"transport_used\": true, \"state\": \"cloned\"},\n{\"lib\": \"beta\", \"git_head\": \"b\", \"remote\": \"upstream-b\", \"selected_transport\": \"upstream-b\", \"transport_used\": false, \"state\": \"verified\"}\n]\n}\n"
        );
    }

    #[test]
    fn failed_staging_write_never_attempts_destination_replacement() {
        let mut replacement_attempted = false;
        let result = stage_and_replace(
            (),
            Path::new("receipt.tmp"),
            Path::new("receipt.json"),
            "complete document",
            |_staging, _document| Err::<(), _>(std::io::Error::other("injected staging failure")),
            || Ok(()),
            |_staging, _staging_path, _seal| Ok(()),
            |_staging, _destination| {
                replacement_attempted = true;
                Ok(())
            },
        );
        assert!(result.is_err());
        assert!(!replacement_attempted);
    }

    #[test]
    fn replacement_failure_reports_both_paths_and_retains_staging_authority() {
        let result = stage_and_replace(
            (),
            Path::new("receipt.tmp.7"),
            Path::new("receipt.json"),
            "complete document",
            |_staging, _document| Ok(()),
            || Ok(()),
            |_staging, _staging_path, _seal| Ok(()),
            |_staging, _destination| Err(std::io::Error::other("injected rename failure")),
        )
        .expect_err("replacement failure must propagate");
        let detail = result.to_string();
        assert!(detail.contains("receipt.tmp.7"), "{detail}");
        assert!(detail.contains("receipt.json"), "{detail}");
    }

    #[test]
    fn publication_barrier_runs_after_durable_staging_and_before_replacement() {
        let events = std::cell::RefCell::new(Vec::new());
        let error = stage_and_replace(
            (),
            Path::new("receipt.tmp.8"),
            Path::new("receipt.json"),
            "complete document",
            |_staging, _document| {
                events.borrow_mut().push("staged");
                Ok(())
            },
            || {
                events.borrow_mut().push("validated");
                Err("source moved".to_string())
            },
            |_staging, _staging_path, _seal| {
                events.borrow_mut().push("verified");
                Ok(())
            },
            |_staging, _destination| {
                events.borrow_mut().push("replaced");
                Ok(())
            },
        )
        .expect_err("a failed final barrier must refuse replacement");
        assert_eq!(*events.borrow(), vec!["staged", "validated"]);
        assert!(error.to_string().contains("source moved"), "{error}");
    }

    #[cfg(any(unix, windows))]
    #[test]
    fn substituted_visible_staging_path_is_refused_before_replacement() {
        let suffix = NEXT_TEMP_SUFFIX.fetch_add(1, Ordering::Relaxed);
        let timestamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("wall clock is after the Unix epoch")
            .as_nanos();
        let fixture = std::env::temp_dir().join(format!(
            "frankensim-bootstrap-stage-substitution-{}-{suffix}-{timestamp}",
            std::process::id()
        ));
        std::fs::create_dir(&fixture).expect("reserve unique fixture directory");
        let destination = fixture.join("receipt.json");
        std::fs::write(&destination, b"prior receipt").expect("seed prior destination");
        let (staging_path, staging_file) =
            reserve_same_directory_temporary(&destination).expect("reserve staging file");
        let displaced_path = fixture.join("original-staging-retained");
        let validation_staging_path = staging_path.clone();
        let validation_displaced_path = displaced_path.clone();
        let replacement_attempted = std::cell::Cell::new(false);

        let error = stage_and_replace(
            staging_file,
            &staging_path,
            &destination,
            "complete document",
            |staging, document| {
                staging.write_all(document)?;
                staging.sync_all()?;
                StagingSeal::capture(staging)
            },
            move || {
                std::fs::rename(&validation_staging_path, &validation_displaced_path)
                    .map_err(|error| error.to_string())?;
                let mut substitute = OpenOptions::new()
                    .write(true)
                    .create_new(true)
                    .open(&validation_staging_path)
                    .map_err(|error| error.to_string())?;
                substitute
                    .write_all(b"substitute")
                    .map_err(|error| error.to_string())?;
                substitute.sync_all().map_err(|error| error.to_string())?;
                Ok(())
            },
            |staging, staging_path, seal| seal.verify(staging, staging_path),
            |_staging, _destination| {
                replacement_attempted.set(true);
                Ok(())
            },
        )
        .expect_err("a substituted staging pathname must refuse publication");

        let detail = error.to_string();
        assert!(
            detail.contains("staging authority check failed"),
            "{detail}"
        );
        assert!(!replacement_attempted.get());
        assert_eq!(
            std::fs::read(&destination).expect("read prior destination"),
            b"prior receipt"
        );
        assert_eq!(
            std::fs::read(&displaced_path).expect("read displaced original"),
            b"complete document"
        );
        assert_eq!(
            std::fs::read(&staging_path).expect("read visible substitute"),
            b"substitute"
        );
    }

    #[cfg(unix)]
    #[test]
    fn non_utf8_provenance_paths_are_refused() {
        use std::ffi::OsString;
        use std::os::unix::ffi::OsStringExt as _;

        let path = PathBuf::from(OsString::from_vec(vec![b'd', 0xff, b'r']));
        let error = provenance_path_text(&path).expect_err("non-UTF-8 path must be refused");
        assert_eq!(error.kind(), std::io::ErrorKind::InvalidInput);
    }
}
