//! One native importer for original, externally produced radiation operators.
//!
//! The artifact independently declares its gas and aperture basis. Import never
//! substitutes the instrument's current atmosphere, fits coefficients, or
//! attests the declared bake identities. The exact byte hash is separate from
//! those identities. Physical admission belongs to `MimoEnvironment` and the
//! existing MIMO terminal runtime when the source is attached to its duct.

use super::{MimoRadiationLimits, MimoRadiationModel};
use crate::pcm_wav::baffled::openings::CircularOpening;
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::CancelGate;
use fs_material::gas_species::input::{GasInputError, GasRecord, ThermochemicalSources, read_gas_record};
use fs_math::c64::C64;
use fs_vfit::{mimo::MimoModel, model::PoleTerm};
use std::str::{Lines, SplitAsciiWhitespace};

/// Exact native grammar selector, independent of filename or instrument name.
pub const MIMO_RADIATION_SCHEMA: &str = "frankensim-mimo-radiation-v1";
/// Bound source bytes before decoding; file-system adapters must also bound reads.
pub const MAX_MIMO_RADIATION_BYTES: usize = 8 * 1024 * 1024;
const MAX_MEMORY_BYTES: usize = 64 * 1024 * 1024;

/// Source syntax, loading, gas admission and cancellation retain their cause.
#[derive(Debug)]
pub enum MimoRadiationInputError {
    /// Invalid or over-budget source; zero denotes whole-artifact/selector input.
    Input { line: usize, what: &'static str },
    /// The original shared gas owner refused this independently declared gas.
    Gas { line: usize, source: GasInputError },
    /// Caller-owned artifact loader refused the original relative path.
    Source { path: String, source: String },
    /// Cancellation before publication of a decoded source.
    Cancelled,
}
impl core::fmt::Display for MimoRadiationInputError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Input { line, what } => write!(f, "MIMO radiation input line {line}: {what}"),
            Self::Gas { line, source } => write!(f, "MIMO radiation gas line {line}: {source}"),
            Self::Source { path, source } => write!(f, "MIMO radiation artifact {path}: {source}"),
            Self::Cancelled => f.write_str("MIMO radiation input cancelled"),
        }
    }
}
impl std::error::Error for MimoRadiationInputError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Gas { source, .. } => Some(source), _ => None }
    }
}
fn bad(line: usize, what: &'static str) -> MimoRadiationInputError {
    MimoRadiationInputError::Input { line, what }
}
fn poll(gate: &CancelGate) -> Result<(), MimoRadiationInputError> {
    if gate.is_requested() { Err(MimoRadiationInputError::Cancelled) } else { Ok(()) }
}

/// Authored `ARTIFACT MAX_STATES MAX_MEMORY_BYTES PASSIVITY_POINTS` selection.
#[derive(Clone, Debug)]
pub struct MimoRadiationInput {
    path: String,
    limits: MimoRadiationLimits,
}
impl MimoRadiationInput {
    /// Decode one exact selector without opening a file or resolving its gas.
    pub fn parse(text: &str) -> Result<Self, MimoRadiationInputError> {
        let mut row = Row { fields: text.split_ascii_whitespace(), line: 0 };
        let path = row.word()?;
        if path.len() > 512 || path.contains(['\\', ':'])
            || path.split('/').any(|s| s.is_empty() || s == "." || s == "..")
        {
            return Err(bad(0, "artifact needs a relative path without empty, dot or parent components"));
        }
        let limits = MimoRadiationLimits {
            max_states: row.count(256)?,
            max_memory_bytes: row.count(MAX_MEMORY_BYTES)?,
            sampled_passivity_points: row.count(65_537)?,
        };
        row.finish()?;
        if limits.max_memory_bytes == 0 || limits.sampled_passivity_points < 17 {
            return Err(bad(0, "positive memory and 17..=65537 sampled passivity points are required"));
        }
        Ok(Self { path: path.into(), limits })
    }
    /// Original relative source address; callers decide how to resolve it.
    pub fn path(&self) -> &str { &self.path }
    /// Original explicit resource limits, passed unchanged to the runtime.
    pub const fn limits(&self) -> MimoRadiationLimits { self.limits }

    /// Load once through the caller and retain separate bytes/source identities.
    /// The caller's registry resolves only the artifact's own `exterior` record.
    /// No caller gas is accepted here, so coefficients cannot be relabeled.
    pub fn load(
        &self,
        sources: Option<&ThermochemicalSources<'_>>,
        loader: &mut impl FnMut(&str) -> Result<Vec<u8>, String>,
        gate: &CancelGate,
    ) -> Result<(MimoRadiationModel, MimoRadiationRecord), MimoRadiationInputError> {
        poll(gate)?;
        let bytes = loader(&self.path).map_err(|source| MimoRadiationInputError::Source {
            path: self.path.clone(), source,
        })?;
        poll(gate)?;
        let artifact = MimoRadiationArtifact::from_bytes_with_sources(&bytes, sources, gate)?;
        let record = MimoRadiationRecord {
            path: self.path.clone(), artifact_hash: artifact.content_hash,
            declaration_identity: artifact.source.declaration_identity,
            realization_identity: artifact.source.realization_identity,
            gas_record: artifact.gas_record, band_hz: artifact.source.band_hz,
            apertures: artifact.source.apertures.clone(), limits: self.limits,
        };
        Ok((artifact.source, record))
    }
}

/// Inspectable native source receipts, without a second retained filter/runtime.
#[derive(Clone, Debug)]
pub struct MimoRadiationRecord {
    path: String,
    artifact_hash: ContentHash,
    declaration_identity: [u8; 32],
    realization_identity: [u8; 32],
    gas_record: GasRecord,
    band_hz: [f64; 2],
    apertures: Vec<CircularOpening>,
    limits: MimoRadiationLimits,
}
impl MimoRadiationRecord {
    /// Original source path, not a source or realization identity.
    pub fn path(&self) -> &str { &self.path }
    /// Domain-separated BLAKE3 of the exact complete imported artifact bytes.
    pub const fn artifact_hash(&self) -> ContentHash { self.artifact_hash }
    /// Original declared geometry/material identity; not an import attestation.
    pub const fn declaration_identity(&self) -> [u8; 32] { self.declaration_identity }
    /// Original declared rational-bake identity; never inferred from its path.
    pub const fn realization_identity(&self) -> [u8; 32] { self.realization_identity }
    /// Original gas declaration, including source-pack receipts when applicable.
    pub const fn gas_record(&self) -> &GasRecord { &self.gas_record }
    /// Original declared physical band [Hz].
    pub const fn band_hz(&self) -> [f64; 2] { self.band_hz }
    /// Original centers/radii in exact matrix row and column order.
    pub fn apertures(&self) -> &[CircularOpening] { &self.apertures }
    /// Authored storage and sampled-passivity resource limits.
    pub const fn limits(&self) -> MimoRadiationLimits { self.limits }
}

/// Decoded source data only. Attachment still performs all physical admissions.
#[derive(Clone, Debug)]
pub struct MimoRadiationArtifact {
    source: MimoRadiationModel,
    gas_record: GasRecord,
    content_hash: ContentHash,
}
impl MimoRadiationArtifact {
    /// Decode bounded LF records, row-major matrices and conjugate-closed poles.
    /// Counts and real state dimension use the terminal owner's hard caps before
    /// allocation. Stability, properness, reciprocity, passivity, graph/gas/clock
    /// and the selected tighter resource budget remain runtime admissions.
    pub fn from_bytes_with_sources(
        bytes: &[u8],
        sources: Option<&ThermochemicalSources<'_>>,
        gate: &CancelGate,
    ) -> Result<Self, MimoRadiationInputError> {
        poll(gate)?;
        if bytes.len() > MAX_MIMO_RADIATION_BYTES || !bytes.ends_with(b"\n")
            || bytes.contains(&b'\r') || bytes.contains(&0)
        {
            return Err(bad(0, "artifact must be at most 8 MiB, LF-delimited and end with LF"));
        }
        let text = core::str::from_utf8(bytes).map_err(|_| bad(0, "artifact must be UTF-8"))?;
        let mut rows = Rows { lines: text.lines(), line: 0 };
        rows.row(MIMO_RADIATION_SCHEMA)?.finish()?;
        let declaration_identity = rows.identity("declaration")?;
        let realization_identity = rows.identity("realization")?;
        let mut row = rows.row("exterior")?;
        let gas_record = read_gas_record(&mut row.fields, None, sources)
            .map_err(|source| MimoRadiationInputError::Gas { line: row.line, source })?;
        row.finish()?;
        rows.verify_gas_source(&gas_record, gate)?;
        let mut row = rows.row("band")?;
        let band_hz = [row.scalar()?, row.scalar()?];
        row.finish()?;
        if band_hz[0] <= 0.0 || band_hz[0] >= band_hz[1] {
            return Err(bad(rows.line, "band needs positive increasing physical frequencies in Hz"));
        }
        let mut row = rows.row("apertures")?;
        let ports = row.count(128)?;
        row.finish()?;
        if ports == 0 { return Err(bad(rows.line, "at least one aperture is required")); }
        let mut apertures = Vec::with_capacity(ports);
        for _ in 0..ports {
            poll(gate)?;
            let mut row = rows.row("aperture")?;
            let aperture = CircularOpening { center_m: [row.scalar()?, row.scalar()?], radius_m: row.scalar()? };
            row.finish()?;
            if aperture.radius_m <= 0.0 { return Err(bad(rows.line, "aperture radius must be positive metres")); }
            apertures.push(aperture);
        }
        let mut row = rows.row("reciprocal")?;
        let reciprocal = match row.word()? {
            "true" => true, "false" => false,
            _ => return Err(bad(row.line, "reciprocal must be true or false")),
        };
        row.finish()?;
        let mut row = rows.row("poles")?;
        let count = row.count(256)?;
        row.finish()?;
        let (mut poles, mut residues) = (Vec::with_capacity(count), Vec::with_capacity(count));
        let cells = ports * ports; // ports was bounded before multiplication
        let mut states = 0;
        for _ in 0..count {
            poll(gate)?;
            let mut row = rows.row("pole")?;
            let pole = match row.word()? {
                "real" => PoleTerm::Real { pole: row.scalar()?, residue: 0.0 },
                "pair" => PoleTerm::Pair { pole: C64::new(row.scalar()?, row.scalar()?), residue: C64::ZERO },
                _ => return Err(bad(row.line, "pole must be real RE or pair RE IM")),
            };
            row.finish()?;
            states += pole.state_dim() * ports;
            if states > 256 { return Err(bad(rows.line, "original real MIMO state dimension exceeds 256")); }
            let mut row = rows.row("residue")?;
            let mut matrix = Vec::with_capacity(cells);
            for _ in 0..cells { matrix.push(C64::new(row.scalar()?, row.scalar()?)); }
            row.finish()?;
            poles.push(pole);
            residues.push(matrix);
        }
        let d = rows.matrix("direct", cells)?;
        let e = rows.matrix("improper", cells)?;
        if rows.lines.next().is_some() { return Err(bad(rows.line + 1, "unexpected trailing artifact record")); }
        poll(gate)?;
        Ok(Self {
            source: MimoRadiationModel {
                declaration_identity, realization_identity, exterior: *gas_record.state(),
                apertures, band_hz, impedance: MimoModel { ports, poles, residues, d, e, reciprocal },
            },
            gas_record,
            content_hash: hash_domain("org.frankensim.fs-couple.mimo-radiation.v1", bytes),
        })
    }
    /// Original source coefficients, atmosphere and ordered aperture basis.
    pub const fn source(&self) -> &MimoRadiationModel { &self.source }
    /// Original independently resolved gas and any thermochemical pack receipts.
    pub const fn gas_record(&self) -> &GasRecord { &self.gas_record }
    /// Domain-separated identity of the complete imported bytes, including IDs.
    pub const fn content_hash(&self) -> ContentHash { self.content_hash }
}

struct Rows<'a> { lines: Lines<'a>, line: usize }
impl<'a> Rows<'a> {
    fn verify_gas_source(&mut self, record: &GasRecord, gate: &CancelGate)
        -> Result<(), MimoRadiationInputError>
    {
        let GasRecord::Thermochemical(resolved) = record else { return Ok(()); };
        // Species names alone do not bind a bake's source data: a caller may
        // supply a different pack under the same id. Compare the independently
        // encoded original receipts in the gas owner's canonical order.
        let mut row = self.row("gas_sources")?;
        let count = row.count(128)?;
        row.finish()?;
        if count != resolved.source_packs().len() {
            return Err(bad(self.line, "thermochemical radiation must bind every original source pack"));
        }
        for source in resolved.source_packs() {
            poll(gate)?;
            let mut row = self.row("gas_source")?;
            if row.word()? != source.species.as_str()
                || row.hash()? != source.species_pack
                || row.hash()? != source.nasa9_pack
                || row.hash()? != source.transport_pack
            {
                return Err(bad(self.line, "thermochemical radiation source pack/order differs from its original artifact"));
            }
            row.finish()?;
        }
        // Receipts do not alone bind a transport claim-selection or water
        // policy. The original complete numerical gas also lives in the
        // artifact and must match bit-for-bit. Keep the resolved GasRecord;
        // never overwrite it or reconstruct a replacement from these numbers.
        let snapshot = self.matrix("gas_state", 12)?;
        if snapshot.iter().zip(gas_values(resolved.state()))
            .any(|(expected, actual)| expected.to_bits() != actual.to_bits())
        {
            return Err(bad(self.line, "thermochemical radiation gas snapshot differs from its original artifact"));
        }
        Ok(())
    }
    fn row(&mut self, key: &'static str) -> Result<Row<'a>, MimoRadiationInputError> {
        self.line += 1;
        let line = self.lines.next().ok_or_else(|| bad(self.line, "missing artifact record"))?;
        let mut row = Row { fields: line.split_ascii_whitespace(), line: self.line };
        if row.word()? != key { return Err(bad(self.line, "unexpected artifact record or order")); }
        Ok(row)
    }
    fn identity(&mut self, key: &'static str) -> Result<[u8; 32], MimoRadiationInputError> {
        let mut row = self.row(key)?;
        let identity = row.hash()?;
        row.finish()?;
        Ok(identity.0)
    }
    fn matrix(&mut self, key: &'static str, cells: usize) -> Result<Vec<f64>, MimoRadiationInputError> {
        let mut row = self.row(key)?;
        let mut matrix = Vec::with_capacity(cells);
        for _ in 0..cells { matrix.push(row.scalar()?); }
        row.finish()?;
        Ok(matrix)
    }
}
struct Row<'a> { fields: SplitAsciiWhitespace<'a>, line: usize }
impl<'a> Row<'a> {
    fn word(&mut self) -> Result<&'a str, MimoRadiationInputError> {
        self.fields.next().ok_or_else(|| bad(self.line, "missing artifact field"))
    }
    fn scalar(&mut self) -> Result<f64, MimoRadiationInputError> {
        let value: f64 = self.word()?.parse().map_err(|_| bad(self.line, "invalid scalar"))?;
        if value.is_finite() { Ok(value) } else { Err(bad(self.line, "scalar must be finite")) }
    }
    fn hash(&mut self) -> Result<ContentHash, MimoRadiationInputError> {
        ContentHash::from_hex(self.word()?)
            .ok_or_else(|| bad(self.line, "identity needs exactly 64 hexadecimal characters"))
    }
    fn count(&mut self, maximum: usize) -> Result<usize, MimoRadiationInputError> {
        let count = self.word()?.parse().map_err(|_| bad(self.line, "invalid unsigned count"))?;
        if count <= maximum { Ok(count) } else { Err(bad(self.line, "count exceeds bounded source budget")) }
    }
    fn finish(mut self) -> Result<(), MimoRadiationInputError> {
        if self.fields.next().is_none() { Ok(()) } else { Err(bad(self.line, "unexpected extra field")) }
    }
}

fn gas_values(g: &fs_material::gas::GasState) -> [f64; 12] {
    [g.temperature, g.pressure, g.density, g.sound_speed, g.dynamic_viscosity,
        g.thermal_conductivity, g.gamma, g.specific_gas_constant, g.specific_heat_cp,
        g.prandtl, g.characteristic_impedance, g.water_mole_fraction]
}

#[cfg(test)]
mod tests;
