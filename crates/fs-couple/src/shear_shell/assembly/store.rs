//! A kept modal solve of an assembled shell: the complete [`SliceReport`] of the tied pencil's window, as exact f64
//! bits, so that a later reduction of the same pencil reads it back and does not solve again.
//!
//! The key binds the solve's own inputs: both matrices of the tied pencil (their dimensions, each row's columns and
//! the bits of its values), the window and every [`SliceOptions`] field, by their `Debug` text, which names each field
//! and writes each finite `f64` so that it reads back to the same bits. A solve is kept as one entry for each mode,
//! under a key made from the solve's key and the mode's place, and then one entry for the report's own numbers under
//! the solve's key, so that neither a keep nor a load holds more than one mode's bytes beside the report. Every entry
//! carries the solve's key, a mode's entry its place, and each ends in a hash of every preceding byte. An entry that
//! does not decode completely, under its key, its place and the pencil's order, refuses, and nothing is solved behind
//! it.
//!
//! The key cannot see a change to the eigensolver that keeps these domains, so a store holds the solves of one build.

use super::{AcousticRealizeError, invalid};
use crate::bowed_string::runtime::reciprocal::bank::radiation::bem::SkinLoadStore;
use fs_blake3::{ContentHash, DomainHasher};
use fs_modal::{ModePair, SliceOptions, SliceReport, SliceStats};
use fs_plate::shell::ShellModel;

const HEAD: &[u8] = b"frankensim-shell-modes-v1\n";
const MODE: &[u8] = b"frankensim-shell-mode-v1\n";
const ENTRY_DOMAIN: &str = "fs-couple:shell-modes-entry:v1";

/// Largest entry a store is asked to keep or to return: one mode of a pencil of about 134 million coordinates. A mode
/// that would not fit refuses its keep.
pub const MAX_SHELL_MODES_ENTRY_BYTES: usize = 1 << 30;

fn word(out: &mut Vec<u8>, n: u64) {
    out.extend_from_slice(&n.to_le_bytes());
}
fn real(out: &mut Vec<u8>, x: f64) {
    word(out, x.to_bits());
}

/// The key of the solve of `model`'s pencil in `window` under `options`.
pub(super) fn key(model: &ShellModel, window: (f64, f64), options: &SliceOptions) -> ContentHash {
    let mut hasher = DomainHasher::new("fs-couple:shell-modes-key:v1");
    let mut row = Vec::new();
    for matrix in [&model.k, &model.m] {
        for n in [matrix.nrows(), matrix.ncols(), matrix.nnz()] {
            hasher.update(&(n as u64).to_le_bytes());
        }
        for r in 0..matrix.nrows() {
            let (columns, values) = matrix.row(r);
            row.clear();
            word(&mut row, columns.len() as u64);
            for (&column, &value) in columns.iter().zip(values) {
                word(&mut row, column as u64);
                real(&mut row, value);
            }
            hasher.update(&row);
        }
    }
    for x in [window.0, window.1] {
        hasher.update(&x.to_bits().to_le_bytes());
    }
    hasher.update(format!("{options:?}").as_bytes());
    hasher.finalize()
}

fn mode_key(key: &ContentHash, place: usize) -> ContentHash {
    let mut hasher = DomainHasher::new("fs-couple:shell-mode-key:v1");
    hasher.update(key.as_bytes());
    hasher.update(&(place as u64).to_le_bytes());
    hasher.finalize()
}

fn sealed(mut out: Vec<u8>) -> Vec<u8> {
    let hash = fs_blake3::hash_domain(ENTRY_DOMAIN, &out);
    out.extend_from_slice(hash.as_bytes());
    out
}

fn failed(mut e: String) -> AcousticRealizeError {
    e.insert_str(0, "shell modes store: ");
    AcousticRealizeError::Nonlinear(e)
}

/// Keep `report` under `key`: every mode first, then the report's own numbers, so that the entry under `key` is
/// written only after every mode it counts.
pub(super) fn keep(
    store: &dyn SkinLoadStore,
    key: &ContentHash,
    report: &SliceReport,
) -> Result<(), AcousticRealizeError> {
    for (place, mode) in report.modes.iter().enumerate() {
        let mut out = Vec::with_capacity(MODE.len() + 32 + 8 * (6 + mode.phi.len()) + 32);
        out.extend_from_slice(MODE);
        out.extend_from_slice(key.as_bytes());
        word(&mut out, place as u64);
        for x in [mode.lambda, mode.residual, mode.interval.0, mode.interval.1] {
            real(&mut out, x);
        }
        word(&mut out, mode.phi.len() as u64);
        for &x in &mode.phi {
            real(&mut out, x);
        }
        let out = sealed(out);
        if out.len() > MAX_SHELL_MODES_ENTRY_BYTES {
            return Err(invalid("a shell mode exceeds the store's entry ceiling"));
        }
        store.keep(&mode_key(key, place), &out).map_err(failed)?;
    }
    let s = &report.stats;
    let mut out = HEAD.to_vec();
    out.extend_from_slice(key.as_bytes());
    for x in [report.window.0, report.window.1, s.shift] {
        real(&mut out, x);
    }
    for n in [
        report.below_low,
        report.below_high,
        report.expected,
        s.factorizations,
        s.lanczos_iters,
        s.restarts,
        s.factor_nnz_l,
        s.factor_peak_bytes,
        s.pivots_delayed,
        report.modes.len(),
    ] {
        word(&mut out, n as u64);
    }
    store.keep(key, &sealed(out)).map_err(failed)
}

/// The report kept under `key` for a pencil of `order` coordinates, if any, read one mode's entry at a time.
pub(super) fn load(
    store: &dyn SkinLoadStore,
    key: &ContentHash,
    order: usize,
) -> Result<Option<SliceReport>, AcousticRealizeError> {
    let damaged = || invalid("kept shell modes are damaged, incomplete or another solve's");
    let Some(head) = store.load(key).map_err(failed)? else {
        return Ok(None);
    };
    let (mut report, count) = read_head(&head, key, order).ok_or_else(damaged)?;
    drop(head);
    for place in 0..count {
        let entry = store.load(&mode_key(key, place)).map_err(failed)?;
        let mode = entry.and_then(|entry| read_mode(&entry, key, place, order));
        report.modes.push(mode.ok_or_else(damaged)?);
    }
    Ok(Some(report))
}

struct Reader<'a>(&'a [u8]);
impl<'a> Reader<'a> {
    /// The body of an entry that ends in its own hash and begins with `magic` and `key`.
    fn open(entry: &'a [u8], magic: &[u8], key: &ContentHash) -> Option<Self> {
        if entry.len() > MAX_SHELL_MODES_ENTRY_BYTES || entry.len() < magic.len() + 64 {
            return None;
        }
        let (body, hash) = entry.split_at(entry.len() - 32);
        let mut r = Self(body);
        let whole = fs_blake3::hash_domain(ENTRY_DOMAIN, body).as_bytes()[..] == *hash;
        (whole && r.take(magic.len())? == magic && r.take(32)? == key.as_bytes()).then_some(r)
    }
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
    fn count(&mut self) -> Option<usize> {
        usize::try_from(self.word()?).ok()
    }
    fn real(&mut self) -> Option<f64> {
        self.word().map(f64::from_bits)
    }
}

/// The report's own numbers, with room for its modes, and their count, which a pencil of `order` coordinates bounds.
fn read_head(entry: &[u8], key: &ContentHash, order: usize) -> Option<(SliceReport, usize)> {
    let mut r = Reader::open(entry, HEAD, key)?;
    let (low, high, shift) = (r.real()?, r.real()?, r.real()?);
    let mut n = [0; 10];
    for value in &mut n {
        *value = r.count()?;
    }
    let report = SliceReport {
        window: (low, high),
        below_low: n[0],
        below_high: n[1],
        expected: n[2],
        modes: Vec::with_capacity(n[9].min(order)),
        stats: SliceStats {
            shift,
            factorizations: n[3],
            lanczos_iters: n[4],
            restarts: n[5],
            factor_nnz_l: n[6],
            factor_peak_bytes: n[7],
            pivots_delayed: n[8],
        },
    };
    (r.0.is_empty() && n[9] <= order).then_some((report, n[9]))
}

fn read_mode(entry: &[u8], key: &ContentHash, place: usize, order: usize) -> Option<ModePair> {
    let mut r = Reader::open(entry, MODE, key)?;
    if r.count()? != place {
        return None;
    }
    let (lambda, residual, low, high) = (r.real()?, r.real()?, r.real()?, r.real()?);
    if r.count()? != order {
        return None;
    }
    let mut phi = Vec::with_capacity(order.min(r.0.len() / 8));
    for _ in 0..order {
        phi.push(r.real()?);
    }
    r.0.is_empty().then_some(ModePair {
        lambda,
        phi,
        residual,
        interval: (low, high),
    })
}
