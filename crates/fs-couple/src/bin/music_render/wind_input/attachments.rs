//! Explicit input snapshots for the local wind adapter.
//! Source callbacks never fall through to the host filesystem in restricted mode.
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File, Metadata};
use std::io::Read;
use std::path::{Path, PathBuf};

pub(super) const MAX_FILES: usize = 32;
const MAX_FILE_BYTES: usize = 8 * 1024 * 1024;
const MAX_TOTAL_BYTES: usize = 32 * 1024 * 1024;

#[derive(Default)]
pub(super) struct AttachedSources {
    files: BTreeMap<String, Vec<u8>>,
}

fn admit_name(name: &str) -> Result<(), String> {
    if name.is_empty()
        || name.len() > 512
        || name.starts_with('-')
        || name
            .bytes()
            .any(|b| b <= b' ' || b == 127 || b == b'\\' || b == b':')
        || name
            .split('/')
            .any(|part| part.is_empty() || part.starts_with('.'))
    {
        return Err(
            "--attached-source needs an explicit relative path without aliases or traversal".into(),
        );
    }
    Ok(())
}

fn checked_path(root: &Path, name: &str) -> Result<(PathBuf, Metadata), String> {
    let mut path = root.to_path_buf();
    let mut parts = name.split('/').peekable();
    let mut last = None;
    while let Some(part) = parts.next() {
        path.push(part);
        let metadata = fs::symlink_metadata(&path)
            .map_err(|e| format!("attached source {name:?} metadata failed: {e}"))?;
        if metadata.file_type().is_symlink()
            || (parts.peek().is_some() && !metadata.is_dir())
            || (parts.peek().is_none() && !metadata.is_file())
        {
            return Err(format!(
                "attached source {name:?} requires regular files and directories without symlinks"
            ));
        }
        last = Some(metadata);
    }
    let canonical = fs::canonicalize(&path)
        .map_err(|e| format!("attached source {name:?} resolution failed: {e}"))?;
    if !canonical.starts_with(root) {
        return Err(format!(
            "attached source {name:?} escapes its performance directory"
        ));
    }
    Ok((canonical, last.expect("admitted nonempty path")))
}

fn same_file(before: &Metadata, after: &Metadata) -> bool {
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        before.dev() == after.dev() && before.ino() == after.ino()
    }
    #[cfg(not(unix))]
    {
        before.is_file() == after.is_file()
            && before.len() == after.len()
            && before.modified().ok() == after.modified().ok()
    }
}

impl AttachedSources {
    /// Snapshot only explicitly named files before any physical source is parsed.
    /// The web adapter supplies a fresh private directory containing regular files.
    pub(super) fn load(parent: &Path, names: &[&str]) -> Result<Self, String> {
        if names.len() > MAX_FILES {
            return Err("--attached-source allows at most 32 files".into());
        }
        let mut admitted = BTreeSet::new();
        for &name in names {
            admit_name(name)?;
            if !admitted.insert(name) {
                return Err(format!("duplicate --attached-source path {name:?}"));
            }
        }
        if names.is_empty() {
            return Ok(Self::default());
        }
        let parent = if parent.as_os_str().is_empty() {
            Path::new(".")
        } else {
            parent
        };
        let root = fs::canonicalize(parent)
            .map_err(|e| format!("attached-source directory resolution failed: {e}"))?;
        let mut files = BTreeMap::new();
        let mut canonical_paths = BTreeSet::new();
        #[cfg(unix)]
        let mut file_ids = BTreeSet::new();
        let mut total = 0usize;
        for &name in names {
            let (path, before) = checked_path(&root, name)?;
            if !canonical_paths.insert(path.clone()) {
                return Err(format!(
                    "attached source {name:?} aliases an earlier source"
                ));
            }
            #[cfg(unix)]
            {
                use std::os::unix::fs::MetadataExt;
                if !file_ids.insert((before.dev(), before.ino())) {
                    return Err(format!(
                        "attached source {name:?} aliases an earlier source"
                    ));
                }
            }
            if before.len() == 0 || before.len() > MAX_FILE_BYTES as u64 {
                return Err(format!(
                    "attached source {name:?} needs 1..=8 MiB of original bytes"
                ));
            }
            let file = File::open(&path)
                .map_err(|e| format!("attached source {name:?} open failed: {e}"))?;
            let opened = file
                .metadata()
                .map_err(|e| format!("attached source {name:?} metadata failed: {e}"))?;
            if !opened.is_file() || !same_file(&before, &opened) {
                return Err(format!("attached source {name:?} changed during admission"));
            }
            let mut bytes = Vec::new();
            file.take((MAX_FILE_BYTES + 1) as u64)
                .read_to_end(&mut bytes)
                .map_err(|e| format!("attached source {name:?} read failed: {e}"))?;
            if bytes.is_empty() || bytes.len() > MAX_FILE_BYTES {
                return Err(format!(
                    "attached source {name:?} needs 1..=8 MiB of original bytes"
                ));
            }
            let (after_path, after) = checked_path(&root, name)?;
            if path != after_path
                || !same_file(&opened, &after)
                || opened.len() != bytes.len() as u64
                || after.len() != bytes.len() as u64
                || opened.modified().ok() != after.modified().ok()
            {
                return Err(format!(
                    "attached source {name:?} changed during its snapshot"
                ));
            }
            total += bytes.len();
            if total > MAX_TOTAL_BYTES {
                return Err("attached sources exceed 32 MiB combined".into());
            }
            files.insert(name.to_owned(), bytes);
        }
        Ok(Self { files })
    }

    /// Exact source name only: never normalize a callback or open a missing path.
    pub(super) fn get(&self, name: &str) -> Option<&[u8]> {
        self.files.get(name).map(Vec::as_slice)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};
    static NEXT: AtomicU64 = AtomicU64::new(0);

    fn directory() -> PathBuf {
        let timestamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "fs-attached-sources-{}-{timestamp}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&dir).unwrap();
        dir
    }

    #[test]
    fn exact_nested_bytes_are_frozen_and_unlisted_callbacks_have_no_authority() {
        let dir = directory();
        fs::create_dir(dir.join("materials")).unwrap();
        let path = dir.join("materials/wall.pack");
        fs::write(&path, [0, 255, 128, 13, 10]).unwrap();
        fs::write(dir.join("unlisted.pack"), b"never admitted").unwrap();
        let sources = AttachedSources::load(&dir, &["materials/wall.pack"]).unwrap();
        fs::write(&path, b"later material replacement").unwrap();
        assert_eq!(
            sources.get("materials/wall.pack"),
            Some([0, 255, 128, 13, 10].as_slice())
        );
        for name in [
            "unlisted.pack",
            "./materials/wall.pack",
            "materials//wall.pack",
            "../materials/wall.pack",
            "/etc/passwd",
        ] {
            assert!(sources.get(name).is_none());
        }
        let changed = AttachedSources::load(&dir, &["materials/wall.pack"]).unwrap();
        assert_eq!(
            changed.get("materials/wall.pack"),
            Some(b"later material replacement".as_slice())
        );
    }

    #[test]
    fn names_counts_and_file_sizes_refuse_before_publication() {
        let absent = Path::new("/directory-that-need-not-exist-for-admission");
        for name in [
            "",
            ".",
            "..",
            "../wall",
            "/wall",
            "a//b",
            "a/./b",
            "a/../b",
            "a\\b",
            "a:b",
            ".hidden",
            "a/.hidden",
            "a b",
            "-flag",
        ] {
            assert!(AttachedSources::load(absent, &[name]).is_err(), "{name:?}");
        }
        assert!(AttachedSources::load(absent, &["wall", "wall"]).is_err());
        assert!(AttachedSources::load(absent, &["wall"; MAX_FILES + 1]).is_err());
        assert!(
            AttachedSources::load(absent, &[])
                .unwrap()
                .get("wall")
                .is_none()
        );
        let dir = directory();
        fs::write(dir.join("empty"), []).unwrap();
        assert!(AttachedSources::load(&dir, &["empty"]).is_err());
        File::create(dir.join("large"))
            .unwrap()
            .set_len((MAX_FILE_BYTES + 1) as u64)
            .unwrap();
        assert!(AttachedSources::load(&dir, &["large"]).is_err());
        assert!(AttachedSources::load(&dir, &["missing"]).is_err());
    }

    #[cfg(unix)]
    #[test]
    fn file_symlinks_directory_symlinks_and_hard_link_aliases_refuse() {
        use std::os::unix::fs::symlink;
        let dir = directory();
        let outside = directory();
        fs::write(outside.join("wall"), b"outside source").unwrap();
        symlink(outside.join("wall"), dir.join("file")).unwrap();
        symlink(&outside, dir.join("directory")).unwrap();
        assert!(AttachedSources::load(&dir, &["file"]).is_err());
        assert!(AttachedSources::load(&dir, &["directory/wall"]).is_err());
        fs::write(dir.join("original"), b"original source").unwrap();
        fs::hard_link(dir.join("original"), dir.join("alias")).unwrap();
        assert!(AttachedSources::load(&dir, &["original", "alias"]).is_err());
    }
}
