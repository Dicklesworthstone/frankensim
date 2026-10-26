//! A kept one-way observer setup: the baked transfer's computed numbers and the
//! receipt's impedance rows, as exact f64 bits.
//!
//! Both are functions of the transfer declaration and the body-mode
//! frequencies alone. The declaration hash binds the skin operator (panels,
//! every source trace, medium, clearance and rule), the sources, the receiver,
//! the clock, the delays, the band and the ladder. The entry key hashes exactly
//! those inputs. The bake's declared fields are rebuilt from the live
//! declaration and never read back. A trailing hash covers every preceding
//! byte. An entry that does not decode completely, under its own key and mode
//! count, refuses, and nothing is recomputed behind it.
//!
//! The key cannot see a change to the boundary solver or to the bake that keeps
//! their hash domains, so a store holds the setups of one build.

use super::{ShellObserverError, invalid};
use fs_bem::observation_transfers::{
    LadderRung, MAX_OBSERVATION_RUNGS, ObservationTransferBake, TransferDeclaration,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_math::c64::C64;
use fs_vfit::observation::{
    FirKernel, MAX_OBSERVATION_ENTRIES, MAX_OBSERVATION_INPUTS, MAX_OBSERVATION_OUTPUTS,
};

const MAGIC: &[u8] = b"frankensim-one-way-observer-setup-v1\n";
const ENTRY_DOMAIN: &str = "fs-couple:one-way-observer-setup-entry:v1";
const MAX_KERNELS: usize = MAX_OBSERVATION_INPUTS * MAX_OBSERVATION_OUTPUTS;

/// Largest entry an admitted declaration can produce: every kernel header and
/// tap, every rung, a complete impedance matrix, the key and the hash.
pub const MAX_ONE_WAY_SETUP_BYTES: usize = MAGIC.len()
    + 32
    + 8
    + 32 * MAX_KERNELS
    + 8 * MAX_OBSERVATION_ENTRIES
    + 8
    + 33 * MAX_OBSERVATION_RUNGS
    + 8
    + 16 * MAX_OBSERVATION_INPUTS * MAX_OBSERVATION_INPUTS
    + 8
    + 32;

/// Where one-way observer setups are kept between renders. The observer
/// verifies every byte a store returns.
pub trait OneWaySetupStore {
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

/// One setup's store key: the transfer declaration and every body-mode angular
/// frequency, in order.
pub(super) fn setup_key(declaration: &TransferDeclaration, omega: &[f64]) -> ContentHash {
    let mut h = DomainHasher::new("fs-couple:one-way-observer-setup:v1");
    h.update(declaration.declaration_hash().as_bytes());
    h.update(&(omega.len() as u64).to_le_bytes());
    for w in omega {
        h.update(&w.to_bits().to_le_bytes());
    }
    h.finalize()
}

fn word(out: &mut Vec<u8>, n: usize) {
    out.extend_from_slice(&(n as u64).to_le_bytes());
}
fn real(out: &mut Vec<u8>, x: f64) {
    out.extend_from_slice(&x.to_bits().to_le_bytes());
}

/// The entry for one computed setup.
pub(super) fn encode(
    key: &ContentHash,
    transfer: &ObservationTransferBake,
    impedance: &[Vec<C64>],
    plain_cbie_limit_rad_per_m: f64,
) -> Vec<u8> {
    let mut out = MAGIC.to_vec();
    out.extend_from_slice(key.as_bytes());
    word(&mut out, transfer.kernels.len());
    for k in &transfer.kernels {
        for n in [k.input, k.output, k.integer_delay_samples, k.taps.len()] {
            word(&mut out, n);
        }
        for &tap in &k.taps {
            real(&mut out, tap);
        }
    }
    word(&mut out, transfer.ladder.len());
    for rung in &transfer.ladder {
        word(&mut out, rung.n_fft);
        real(&mut out, rung.tail_energy_fraction);
        real(&mut out, rung.wrap_energy_fraction);
        out.push(u8::from(rung.dtft_relative_change.is_some()));
        real(&mut out, rung.dtft_relative_change.unwrap_or(0.0));
    }
    word(&mut out, impedance.len());
    for z in impedance.iter().flatten() {
        real(&mut out, z.re);
        real(&mut out, z.im);
    }
    real(&mut out, plain_cbie_limit_rad_per_m);
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
        self.word().map(f64::from_bits).filter(|x| x.is_finite())
    }
}

/// The computed part of one kept setup.
#[derive(Clone, Debug)]
pub(super) struct KeptSetup {
    kernels: Vec<FirKernel>,
    ladder: Vec<LadderRung>,
    impedance: Vec<Vec<C64>>,
    plain_cbie_limit_rad_per_m: f64,
}

impl KeptSetup {
    /// Decode the entry kept under `key` for `modes` body modes.
    pub(super) fn decode(
        entry: &[u8],
        key: &ContentHash,
        modes: usize,
    ) -> Result<Self, ShellObserverError> {
        Self::read(entry, key, modes).ok_or_else(|| {
            invalid("a kept one-way setup is damaged, incomplete or another setup's")
        })
    }

    fn read(entry: &[u8], key: &ContentHash, modes: usize) -> Option<Self> {
        if entry.len() > MAX_ONE_WAY_SETUP_BYTES || entry.len() < MAGIC.len() + 64 {
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
        let mut kernels = Vec::new();
        let mut taps_total = 0_usize;
        for _ in 0..r.count(MAX_KERNELS)? {
            let (input, output, integer_delay_samples) = (
                r.count(MAX_OBSERVATION_INPUTS)?,
                r.count(MAX_OBSERVATION_OUTPUTS)?,
                r.count(usize::MAX)?,
            );
            let taps = r.count(MAX_OBSERVATION_ENTRIES)?;
            taps_total = taps_total
                .checked_add(taps)
                .filter(|&n| n <= MAX_OBSERVATION_ENTRIES)?;
            kernels.push(FirKernel {
                input,
                output,
                integer_delay_samples,
                taps: (0..taps).map(|_| r.real()).collect::<Option<_>>()?,
            });
        }
        let mut ladder = Vec::new();
        for _ in 0..r.count(MAX_OBSERVATION_RUNGS)? {
            let (n_fft, tail_energy_fraction, wrap_energy_fraction) =
                (r.count(usize::MAX)?, r.real()?, r.real()?);
            let (tag, change) = (r.take(1)?[0], r.real()?);
            ladder.push(LadderRung {
                n_fft,
                tail_energy_fraction,
                wrap_energy_fraction,
                dtft_relative_change: match tag {
                    0 if change.to_bits() == 0 => None,
                    1 => Some(change),
                    _ => return None,
                },
            });
        }
        if r.count(MAX_OBSERVATION_INPUTS)? != modes {
            return None;
        }
        let impedance = (0..modes)
            .map(|_| {
                (0..modes)
                    .map(|_| Some(C64::new(r.real()?, r.real()?)))
                    .collect::<Option<Vec<_>>>()
            })
            .collect::<Option<Vec<_>>>()?;
        let plain_cbie_limit_rad_per_m = r.real()?;
        r.0.is_empty().then_some(Self {
            kernels,
            ladder,
            impedance,
            plain_cbie_limit_rad_per_m,
        })
    }

    /// The bake for `declaration` exactly as `BemObservation::bake` returns it,
    /// the impedance rows and the plain-CBIE limit.
    pub(super) fn into_parts(
        self,
        declaration: &TransferDeclaration,
    ) -> (ObservationTransferBake, Vec<Vec<C64>>, f64) {
        let transfer = ObservationTransferBake {
            declaration_hash: declaration.declaration_hash(),
            sources: declaration.sources.clone(),
            receivers: declaration.receivers.clone(),
            fs_hz: declaration.fs_hz,
            geometric_delays_s: declaration.geometric_delays_s.clone(),
            filter_latency_samples: declaration.filter_latency_samples,
            pass_band_upper_hz: declaration.pass_band_upper_hz,
            taper_stop_hz: declaration.taper_stop_hz,
            producer_uncertainty: declaration.producer_uncertainty,
            kernels: self.kernels,
            ladder: self.ladder,
        };
        (transfer, self.impedance, self.plain_cbie_limit_rad_per_m)
    }
}
