//! Product file adapter for the existing fs-uq checkpoint format.
//! A checksum detects corruption, not authenticity; resume only trusted files.

use super::{Failure, MAX_PRODUCT_SAMPLES, Result};
use super::mean_control::MeanControl;
use fs_blake3::{ContentHash, DomainHasher};
use fs_uq::{QmcConfig, QmcExecution, SobolExecution, UqExecution, UqPlan};
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

// Raw observations use eight bytes each. The 4 KiB framing allowance also
// covers the frozen adjoint and all 256 admitted parameter coefficients.
// Reads remain bounded; each owner checks its exact version and payload length.
const MAX_CHECKPOINT_BYTES: u64 = 4096 + 8 * MAX_PRODUCT_SAMPLES as u64;

fn failure(message: impl Into<String>) -> Failure {
    Failure { code: "cooling-network-uq-checkpoint", message: message.into() }
}

/// Bind fixed inputs and target lowering to the exact executable used for the
/// child solves. The library additionally binds EVERY UqPlan field. Per-run
/// wall time and chunk length are resource controls, not sampling semantics.
pub(super) fn model_identity(base_text: &str, bindings: &str) -> Result<ContentHash> {
    let executable = super::child::executable_path()
        .map_err(|error| failure(format!("cannot locate evaluator: {error}")))?;
    let mut file = File::open(&executable)
        .map_err(|error| failure(format!("cannot read evaluator identity: {error}")))?;
    let mut binary = DomainHasher::new("org.frankensim.cooling-uq.executable.v1");
    let mut buffer = [0_u8; 65_536];
    loop {
        let count = file.read(&mut buffer)
            .map_err(|error| failure(format!("cannot hash evaluator: {error}")))?;
        if count == 0 { break; }
        binary.update(&buffer[..count]);
    }
    let mut identity = DomainHasher::new("org.frankensim.cooling-uq.model.v1");
    identity.update(binary.finalize().as_bytes());
    for text in [base_text, bindings] {
        identity.update(&(text.len() as u64).to_le_bytes());
        identity.update(text.as_bytes());
    }
    Ok(identity.finalize())
}

pub(super) fn restore(path: &Path, plan: &UqPlan, identity: ContentHash) -> Result<UqExecution> {
    UqExecution::restore(plan, identity, &read(path)?)
        .map_err(|error| failure(format!("{}: {error}", path.display())))
}

pub(super) fn restore_controlled(
    path: &Path, plan: &UqPlan, identity: ContentHash,
) -> Result<(UqExecution, MeanControl)> {
    MeanControl::restore_bytes(plan, identity, &read(path)?)
        .map_err(|error| failure(format!("{}: {error}", path.display())))
}

pub(super) fn restore_qmc(
    path: &Path, plan: &UqPlan, layout: QmcConfig, identity: ContentHash,
) -> Result<QmcExecution> {
    QmcExecution::restore(plan, layout, identity, &read(path)?)
        .map_err(|error| failure(format!("{}: {error}", path.display())))
}

pub(super) fn restore_sobol(
    path: &Path, plan: &UqPlan, identity: ContentHash,
) -> Result<SobolExecution> {
    SobolExecution::restore(plan, identity, &read(path)?)
        .map_err(|error| failure(format!("{}: {error}", path.display())))
}

fn read(path: &Path) -> Result<Vec<u8>> {
    let file = File::open(path)
        .map_err(|error| failure(format!("cannot open {}: {error}", path.display())))?;
    let mut bytes = Vec::new();
    file.take(MAX_CHECKPOINT_BYTES + 1).read_to_end(&mut bytes)
        .map_err(|error| failure(format!("cannot read {}: {error}", path.display())))?;
    if bytes.len() as u64 > MAX_CHECKPOINT_BYTES {
        return Err(failure("checkpoint exceeds the product observation budget"));
    }
    Ok(bytes)
}

/// One fresh output owned by this invocation. Never opens an existing file
/// for replacement, so a resume input or another user's file cannot be lost.
/// Successful sample boundaries replace this output atomically; a partial
/// staging write cannot destroy the last valid prefix.
pub(super) struct Output {
    path: PathBuf,
    staging: PathBuf,
}

impl Output {
    pub(super) fn reserve(path: &Path) -> Result<Self> {
        let name = path.file_name().ok_or_else(|| failure("checkpoint needs a file name"))?;
        let mut staging_name = name.to_os_string();
        staging_name.push(".pending");
        let staging = path.with_file_name(staging_name);
        // Validate the staging path before reserving the output. An abandoned
        // staging file is left intact for diagnosis, never silently removed.
        match fs::symlink_metadata(&staging) {
            Ok(_) => return Err(failure(format!("checkpoint staging path already exists: {}", staging.display()))),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => return Err(failure(format!("cannot inspect staging path: {error}"))),
        }
        OpenOptions::new().write(true).create_new(true).open(path)
            .map_err(|error| failure(format!("checkpoint output must be new ({}): {error}", path.display())))?;
        Ok(Self { path: path.to_path_buf(), staging })
    }

    pub(super) fn save(&self, execution: &UqExecution, identity: ContentHash) -> Result<()> {
        let bytes = execution.checkpoint(identity).map_err(|error| failure(error.to_string()))?;
        self.publish(&bytes)
    }

    pub(super) fn save_mc(
        &self, execution: &UqExecution, identity: ContentHash, control: Option<&MeanControl>,
    ) -> Result<()> {
        match control {
            Some(control) => self.publish(&control.checkpoint_bytes(execution, identity)?),
            None => self.save(execution, identity),
        }
    }

    pub(super) fn save_qmc(&self, execution: &QmcExecution, identity: ContentHash) -> Result<()> {
        let bytes = execution.checkpoint(identity).map_err(|error| failure(error.to_string()))?;
        self.publish(&bytes)
    }

    pub(super) fn save_sobol(&self, execution: &SobolExecution, identity: ContentHash) -> Result<()> {
        let bytes = execution.checkpoint(identity).map_err(|error| failure(error.to_string()))?;
        self.publish(&bytes)
    }

    /// A genuine failed solve is terminal. Retain its diagnosis in a format
    /// the checkpoint decoder refuses rather than allowing sample filtering.
    pub(super) fn invalidate(&self, message: &str) -> Result<()> {
        self.publish(format!("FRANKENSIM-UQ-FAILED\n{message}\n").as_bytes())
    }

    pub(super) fn publish(&self, bytes: &[u8]) -> Result<()> {
        let mut staging = OpenOptions::new().write(true).create_new(true).open(&self.staging)
            .map_err(|error| failure(format!("cannot stage checkpoint: {error}")))?;
        staging.write_all(bytes).and_then(|()| staging.sync_all())
            .map_err(|error| failure(format!("cannot persist checkpoint: {error}")))?;
        drop(staging);
        fs::rename(&self.staging, &self.path)
            .map_err(|error| failure(format!("cannot publish checkpoint: {error}")))?;
        // This CLI's cooling child uses /dev/stdin. On its Unix platforms,
        // syncing the directory also persists the rename, not only the bytes.
        #[cfg(unix)]
        {
            let parent = self.path.parent().filter(|p| !p.as_os_str().is_empty())
                .unwrap_or_else(|| Path::new("."));
            File::open(parent).and_then(|file| file.sync_all())
                .map_err(|error| failure(format!("cannot persist checkpoint directory: {error}")))?;
        }
        Ok(())
    }
}
