//! The complete string performer consumes the same source snapshots as wind rendering.
//! Ordinary CLI paths retain their existing resolution; restricted callbacks never open files.
use super::{Result, attachments::AttachedSources, read_bounded};
use std::{ffi::OsString, path::{Path, PathBuf}};

#[derive(Default)]
pub(super) struct SourceOptions {
    restricted: bool,
    names: Vec<String>,
    pub(super) decimate: bool,
}

impl SourceOptions {
    pub(super) fn extract(args: &[OsString]) -> Result<(Vec<OsString>, Self)> {
        let mut options = Self::default();
        let mut retained = args.iter().take(3).cloned().collect::<Vec<_>>();
        let mut rest = args.iter().skip(3);
        while let Some(arg) = rest.next() {
            if arg == "--self-contained" {
                if options.restricted { return Err("--self-contained cannot be repeated".into()); }
                options.restricted = true;
            } else if arg == "--decimate" {
                if options.decimate { return Err("--decimate cannot be repeated".into()); }
                options.decimate = true;
            } else if arg == "--attached-source" {
                let name = rest.next().and_then(|s| s.to_str())
                    .ok_or("--attached-source needs a UTF-8 relative path")?;
                options.names.push(name.to_owned());
                if options.names.len() > super::attachments::MAX_FILES {
                    return Err("at most 32 attached source files may be supplied".into());
                }
            } else { retained.push(arg.clone()); }
        }
        if !options.restricted && !options.names.is_empty() {
            return Err("--attached-source requires --self-contained".into());
        }
        Ok((retained, options))
    }

    pub(super) fn load(&self, root: &Path) -> Result<SourceFiles> {
        let attached = if self.restricted {
            Some(AttachedSources::load(root, &self.names.iter().map(String::as_str).collect::<Vec<_>>())?)
        } else { None };
        Ok(SourceFiles { root: root.to_owned(), attached })
    }
}

pub(super) struct SourceFiles {
    root: PathBuf,
    attached: Option<AttachedSources>,
}

impl SourceFiles {
    pub(super) fn is_restricted(&self) -> bool { self.attached.is_some() }

    /// CLI profile/gas paths are relative to the caller in ordinary mode and exact
    /// attachment names in restricted mode. Each consuming parser keeps its own cap.
    pub(super) fn read_option(&self, path: &Path, cap: usize) -> Result<Vec<u8>> {
        let Some(attached) = &self.attached else { return read_bounded(path, cap); };
        let name = path.to_str().ok_or("attached source paths must be UTF-8")?;
        let bytes = attached.get(name).ok_or_else(||
            format!("--self-contained refuses unlisted source {name:?}"))?;
        if bytes.len() > cap { return Err(format!("attached source {name:?} exceeds {cap} bytes").into()); }
        Ok(bytes.to_vec())
    }

    /// Material and interface names in a performance are relative to that file.
    pub(super) fn read_relative(&self, name: &str, cap: usize) -> Result<Vec<u8>> {
        if self.is_restricted() { self.read_option(Path::new(name), cap) }
        else { read_bounded(&self.root.join(name), cap) }
    }
}
