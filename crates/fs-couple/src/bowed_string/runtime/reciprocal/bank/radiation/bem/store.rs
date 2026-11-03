//! A kept two-way skin load: the fitted passive model and the numbers its
//! report states, as exact f64 bits.
//!
//! The model and every kept number are functions of the loading declaration
//! alone. Its identity binds the basis identity, every skin triangle, every
//! trace, the frequency grid, the medium, the matrix allowance and the angular
//! refinement. The key is that identity, so an entry is found before any
//! boundary work. The report's declared fields are rebuilt from the live
//! inputs and never read back. A trailing hash covers every preceding byte.
//! An entry that does not decode completely, under its own key and port count,
//! refuses, and nothing is recomputed behind it.
//!
//! The key cannot see a change to the boundary solver or the fit that keeps
//! their hash domains, so a store holds the loads of one build.

use super::ShellBemError;
use crate::render::plate::impact::radiation::{MAX_POLES, MAX_PORTS, Model, Pole};
use fs_bem::modal_loading::ReciprocalLoadingDiagnostics;
use fs_blake3::ContentHash;

const MAGIC: &[u8] = b"frankensim-skin-exterior-load-v1\n";
const ENTRY_DOMAIN: &str = "fs-couple:skin-exterior-load-entry:v1";

/// Largest entry an admitted load can produce: the key, the diagnostics, the
/// gate numbers, the complete pole and port model and the hash.
pub const MAX_SKIN_LOAD_BYTES: usize =
    MAGIC.len() + 32 + 8 * 8 + 8 * 5 + 8 + 8 + MAX_POLES * (8 + 8 + 8 * MAX_PORTS) + 32;

/// Where computed skin loads are kept between renders. The load verifies
/// every byte a store returns.
pub trait SkinLoadStore {
    /// The entry kept under `key`, if any.
    ///
    /// # Errors
    /// The store's own read failure.
    fn load(&self, key: &ContentHash) -> Result<Option<Vec<u8>>, String>;
    /// Keep a new entry under `key`. An entry already kept there stays as it is.
    ///
    /// # Errors
    /// The store's own write failure.
    fn keep(&self, key: &ContentHash, entry: &[u8]) -> Result<(), String>;
}

/// The computed numbers of one load's report, beside its model.
#[derive(Clone, Copy, Debug)]
pub(super) struct LoadNumbers {
    pub(super) reciprocal: ReciprocalLoadingDiagnostics,
    pub(super) reciprocity_relative_residual: f64,
    pub(super) peak_error: f64,
    pub(super) rms_error: f64,
    pub(super) resistance_peak_error: f64,
    pub(super) resistance_rms_error: f64,
}

fn word(out: &mut Vec<u8>, n: u64) {
    out.extend_from_slice(&n.to_le_bytes());
}
fn real(out: &mut Vec<u8>, x: f64) {
    word(out, x.to_bits());
}

/// The entry for one computed load.
pub(super) fn encode(key: &ContentHash, model: &Model, numbers: &LoadNumbers) -> Vec<u8> {
    let d = &numbers.reciprocal;
    let mut out = MAGIC.to_vec();
    out.extend_from_slice(key.as_bytes());
    for n in [
        u64::from(d.angular_refinement),
        d.unknowns as u64,
        d.cache_bytes as u64,
        d.frequency_work,
        d.sweep_work,
    ] {
        word(&mut out, n);
    }
    for x in [
        d.minimum_panels_per_wavelength,
        d.maximum_backward_residual,
        d.maximum_flux_balance,
        numbers.reciprocity_relative_residual,
        numbers.peak_error,
        numbers.rms_error,
        numbers.resistance_peak_error,
        numbers.resistance_rms_error,
    ] {
        real(&mut out, x);
    }
    word(&mut out, model.ports as u64);
    word(&mut out, model.poles.len() as u64);
    for pole in &model.poles {
        real(&mut out, pole.omega);
        real(&mut out, pole.zeta);
        for &c in &pole.coupling {
            real(&mut out, c);
        }
    }
    let hash = fs_blake3::hash_domain(ENTRY_DOMAIN, &out);
    out.extend_from_slice(hash.as_bytes());
    out
}

struct Reader<'a>(&'a [u8]);
impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> Option<&'a [u8]> {
        if n > self.0.len() {
            return None;
        }
        let (head, tail) = self.0.split_at(n);
        self.0 = tail;
        Some(head)
    }
    fn word(&mut self) -> Option<u64> {
        Some(u64::from_le_bytes(self.take(8)?.try_into().ok()?))
    }
    fn count(&mut self, cap: usize) -> Option<usize> {
        usize::try_from(self.word()?).ok().filter(|&n| n <= cap)
    }
    fn real(&mut self) -> Option<f64> {
        self.word().map(f64::from_bits)
    }
}

/// Decode the entry kept under `key` for a load on `ports` ports. The model
/// must pass its own validation; the caller applies its pole allowance.
pub(super) fn decode(
    entry: &[u8],
    key: &ContentHash,
    ports: usize,
    angular_refinement: u32,
) -> Result<(Model, LoadNumbers), ShellBemError> {
    read(entry, key, ports, angular_refinement)
        .filter(|(model, _)| model.validate().is_ok())
        .ok_or(ShellBemError::Invalid(
            "a kept skin load is damaged, incomplete or another load's",
        ))
}

fn read(
    entry: &[u8],
    key: &ContentHash,
    ports: usize,
    angular_refinement: u32,
) -> Option<(Model, LoadNumbers)> {
    if entry.len() > MAX_SKIN_LOAD_BYTES || entry.len() < MAGIC.len() + 64 {
        return None;
    }
    let (body, hash) = entry.split_at(entry.len() - 32);
    if fs_blake3::hash_domain(ENTRY_DOMAIN, body).as_bytes()[..] != *hash {
        return None;
    }
    let mut r = Reader(body);
    if r.take(MAGIC.len())? != MAGIC || r.take(32)? != key.as_bytes() {
        return None;
    }
    let refinement = u32::try_from(r.word()?).ok()?;
    let (unknowns, cache_bytes) = (r.count(usize::MAX)?, r.count(usize::MAX)?);
    let (frequency_work, sweep_work) = (r.word()?, r.word()?);
    let mut x = [0.0; 8];
    for value in &mut x {
        *value = r.real()?;
    }
    if refinement != angular_refinement || r.count(MAX_PORTS)? != ports {
        return None;
    }
    let poles = (0..r.count(MAX_POLES)?)
        .map(|_| {
            Some(Pole {
                omega: r.real()?,
                zeta: r.real()?,
                coupling: (0..ports).map(|_| r.real()).collect::<Option<_>>()?,
            })
        })
        .collect::<Option<Vec<_>>>()?;
    let numbers = LoadNumbers {
        reciprocal: ReciprocalLoadingDiagnostics {
            angular_refinement: refinement,
            unknowns,
            cache_bytes,
            frequency_work,
            sweep_work,
            minimum_panels_per_wavelength: x[0],
            maximum_backward_residual: x[1],
            maximum_flux_balance: x[2],
        },
        reciprocity_relative_residual: x[3],
        peak_error: x[4],
        rms_error: x[5],
        resistance_peak_error: x[6],
        resistance_rms_error: x[7],
    };
    r.0.is_empty().then_some((Model { ports, poles }, numbers))
}
