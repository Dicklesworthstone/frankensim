//! Version-bump machinery for the `.fsim` envelope, proven before it is ever
//! needed: readers admit exactly [`crate::FSIM_VERSION`], and any older
//! envelope must pass through [`migrate_envelope`], which re-emits canonical
//! bytes under the current version and binds both byte strings into a
//! [`ProjectMigrationReceipt`]. There is no implicit migration path.
//!
//! Version 1 is the first real schema version, so the only registered rule is
//! the synthetic proof rule from version 0: a version-0 envelope differing
//! from version 1 solely in its declared version. It exists to prove the
//! receipt machinery end to end (bead f85xj.6.1's acceptance criterion), and
//! it says so in its name.

use fs_blake3::ContentHash;

use crate::FSIM_VERSION;
use crate::wire::{DecodedProject, ProjectError, canonical_hash, parse_sexpr};

/// The named semantic rule a migration applied.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MigrationRule {
    /// Synthetic version-0 proof rule: the payload grammar is identical and
    /// only the envelope version is rewritten. Registered to prove the
    /// machinery; no released artifact ever carried version 0.
    SyntheticV0EnvelopeRewrite,
    /// Version 1 to 2: the cooling section gains the optional
    /// `(fan-system ...)` subsection (bead frn2i.1). Every version-1
    /// document is a valid version-2 document unchanged, so the rule is the
    /// receipted envelope rewrite — there is no semantic compatibility
    /// default for documents that do carry the new section: version-1
    /// readers refuse them as unknown subsections.
    CoolingFanSystemV2,
}

impl MigrationRule {
    /// Stable machine-readable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            MigrationRule::SyntheticV0EnvelopeRewrite => "synthetic-v0-envelope-rewrite",
            MigrationRule::CoolingFanSystemV2 => "cooling-fan-system-v2",
        }
    }

    /// The source version this rule migrates from.
    #[must_use]
    pub const fn source_version(self) -> u32 {
        match self {
            MigrationRule::SyntheticV0EnvelopeRewrite => 0,
            MigrationRule::CoolingFanSystemV2 => 1,
        }
    }
}

/// Receipt binding one migration's exact input and output bytes to its rule.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProjectMigrationReceipt {
    /// Version migrated from.
    pub source_version: u32,
    /// Version migrated to (always [`FSIM_VERSION`]).
    pub target_version: u32,
    /// Hash of the exact historical bytes.
    pub old_hash: ContentHash,
    /// Hash of the exact canonical bytes after migration.
    pub new_hash: ContentHash,
    /// The named rule that was applied.
    pub rule: MigrationRule,
}

impl ProjectMigrationReceipt {
    /// Re-verify this receipt against the exact byte strings it binds.
    #[must_use]
    pub fn verifies(&self, old_bytes: &[u8], new_bytes: &[u8]) -> bool {
        self.source_version == self.rule.source_version()
            && self.target_version == FSIM_VERSION
            && canonical_hash(old_bytes) == self.old_hash
            && canonical_hash(new_bytes) == self.new_hash
    }
}

/// A migrated project: the decoded current-version document plus the receipt
/// that proves where it came from.
#[derive(Debug, Clone, PartialEq)]
pub struct MigratedProject {
    /// The decoded project at the current version.
    pub decoded: DecodedProject,
    /// The migration receipt.
    pub receipt: ProjectMigrationReceipt,
}

/// Parse a canonical `.fsim` document, migrating an older envelope through
/// its registered receipted rule when needed (the `auto-migration-receipt`
/// obligation). The returned document is always the current version; the
/// migration receipt, when one fired, is retained for the caller's ledger.
///
/// # Errors
/// Returns [`ProjectError`] from the strict parser, or from
/// [`migrate_envelope`] when the declared version has no registered rule.
pub fn parse_sexpr_migrating(source: &str) -> Result<MigratedOrNative, ProjectError> {
    match parse_sexpr(source) {
        Ok(decoded) => Ok(MigratedOrNative {
            decoded,
            migration: None,
        }),
        Err(error) if error.code == "fsim-unsupported-version" => {
            let declared = declared_envelope_version(source).ok_or_else(|| ProjectError {
                code: "fsim-migration-shape",
                detail: "the envelope version could not be read for migration".to_string(),
                hint:
                    "migrate exactly the bytes that were persisted; do not hand-edit the envelope"
                        .to_string(),
            })?;
            let migrated = migrate_envelope(source, declared)?;
            Ok(MigratedOrNative {
                decoded: migrated.decoded,
                migration: Some(migrated.receipt),
            })
        }
        Err(error) => Err(error),
    }
}

/// A parsed project plus the migration receipt when an older envelope was
/// rewritten by a registered rule.
#[derive(Debug, Clone, PartialEq)]
pub struct MigratedOrNative {
    /// The decoded project at the current version.
    pub decoded: DecodedProject,
    /// The migration receipt, present exactly when a rule fired.
    pub migration: Option<ProjectMigrationReceipt>,
}

/// Read the envelope's declared `:version` integer without full parsing.
fn declared_envelope_version(source: &str) -> Option<u32> {
    let prefix = source.strip_prefix("(fsim-project :version ")?;
    let digits: String = prefix
        .chars()
        .take_while(|character| character.is_ascii_digit())
        .collect();
    digits.parse().ok()
}

/// Migrate an older envelope to the current version under a registered rule,
/// refusing unknown versions. The output is strictly canonical.
pub fn migrate_envelope(
    source: &str,
    declared_version: u32,
) -> Result<MigratedProject, ProjectError> {
    let rule = match declared_version {
        0 => MigrationRule::SyntheticV0EnvelopeRewrite,
        1 => MigrationRule::CoolingFanSystemV2,
        v if v == FSIM_VERSION => {
            return Err(ProjectError {
                code: "fsim-migration-not-needed",
                detail: format!("the document already declares version {v}"),
                hint: "parse it directly; migration is only for older envelopes".to_string(),
            });
        }
        v => {
            return Err(ProjectError {
                code: "fsim-migration-unknown-version",
                detail: format!("no registered migration rule covers version {v}"),
                hint: "register an explicit rule with its receipt semantics before reading this document".to_string(),
            });
        }
    };

    let old_prefix = format!("(fsim-project :version {declared_version}");
    let new_prefix = format!("(fsim-project :version {FSIM_VERSION}");
    let Some(rest) = source.strip_prefix(old_prefix.as_str()) else {
        return Err(ProjectError {
            code: "fsim-migration-shape",
            detail: format!(
                "the document does not open with `{old_prefix}` as its declared version requires"
            ),
            hint: "migrate exactly the bytes that were persisted; do not hand-edit the envelope"
                .to_string(),
        });
    };
    let migrated = format!("{new_prefix}{rest}");
    let decoded = parse_sexpr(&migrated)?;
    let receipt = ProjectMigrationReceipt {
        source_version: declared_version,
        target_version: FSIM_VERSION,
        old_hash: canonical_hash(source.as_bytes()),
        new_hash: canonical_hash(decoded.canonical.as_bytes()),
        rule,
    };
    Ok(MigratedProject { decoded, receipt })
}
