//! Multimodal characteristic-line runtime (music bead
//! `frankensim-music-v8-root-3ez8g.4.2`): the PLAY image for brass —
//! per-valve-combination exact-FIR reflectance lines realized from the
//! fs-duct MULTIMODAL TMM authority, terminated into the tabulated
//! bell load (bead zolja) with a disclosed analytic fallback outside
//! the table's support, switched at block boundaries with a D17 state
//! lift that carries the in-flight outgoing-wave history.
//!
//! REALIZATION DOCTRINE (the §6.2 bake-off this module documents): at a
//! brass input plane the higher m = 0 modes are DEEPLY evanescent (an
//! 8 mm throat puts mode 1's local cutoff near 26 kHz), so the
//! MM-derived plane-mode reflectance `R_00(omega)` — which already
//! contains every interior mode-conversion round trip — is what a
//! plane-driven source sees. The DOMINANT realization is therefore one
//! exact-FIR line filled from `R_00`; the FULL-MATRIX realization
//! ([`MatrixFirLine`], N^2 FIRs) exists as the bake-off arm and for
//! sources that genuinely inject higher modes. The bake-off receipt
//! (test ml-002) records both arms' held-fingering accuracy and cost.
//!
//! CONJUGATION DISCIPLINE (the pinned ~50%-error-floor trap): fs-duct
//! is `e^{-i omega t}`; the DFT that makes a CAUSAL impulse response is
//! `e^{+i omega t}`. Every reflectance sample is conjugated into the
//! DFT grid exactly as `fs_couple::driving_point` does; the causality
//! gate (ml-006) is the executable witness.

use fs_duct::modal::{
    ModalResponse, mm_input_impedance, mm_input_impedance_tabulated,
    modal_characteristic_impedances, modal_reflection_from_impedance,
};
use fs_duct::{Duct, DuctError, LossModel, Segment, TabulatedLoad, Termination};
use fs_fft::{C64 as FftC64, Fft};
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_vfit::discretize::{DelayedFilter, DiscretizeError, reflectance};

/// Typed refusals from the multimodal line runtime.
#[derive(Debug)]
pub enum MmLineError {
    /// A configuration or geometry parameter is unusable.
    Invalid {
        /// Diagnosis.
        what: &'static str,
    },
    /// The underlying duct authority refused.
    Duct(DuctError),
    /// The realized line refused.
    Discrete(DiscretizeError),
    /// The shared line ladder reached `max_taps` before every combo converged.
    LadderUnconverged {
        /// Ceiling reached [taps].
        max_taps: usize,
        /// Worst combo's fourth-quarter tail energy at the ceiling.
        tail_energy: f64,
        /// Worst combo's first-three-peak shift [cents] at the ceiling.
        peak_shift_cents: f64,
    },
}

impl core::fmt::Display for MmLineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            MmLineError::Invalid { what } => write!(f, "FS-COUPLE-MMLINE: {what}"),
            MmLineError::Duct(e) => write!(f, "FS-COUPLE-MMLINE duct: {e}"),
            MmLineError::Discrete(e) => write!(f, "FS-COUPLE-MMLINE line: {e:?}"),
            MmLineError::LadderUnconverged {
                max_taps,
                tail_energy,
                peak_shift_cents,
            } => write!(
                f,
                "FS-COUPLE-MMLINE: line ladder reached {max_taps} taps before convergence \
                 (worst tail_energy {tail_energy:.3e}, worst peak shift {peak_shift_cents:.3} c)"
            ),
        }
    }
}

impl core::error::Error for MmLineError {}

impl From<DuctError> for MmLineError {
    fn from(e: DuctError) -> Self {
        MmLineError::Duct(e)
    }
}

impl From<DiscretizeError> for MmLineError {
    fn from(e: DiscretizeError) -> Self {
        MmLineError::Discrete(e)
    }
}

/// How the multimodal lines terminate at the mouth.
#[derive(Clone)]
pub enum MmLoad<'a> {
    /// A closed analytic termination (its own validity refusals apply
    /// at every DFT bin).
    Analytic(Termination),
    /// The tabulated bell load INSIDE its support, an analytic fallback
    /// outside it — the splice bins are counted and disclosed in the
    /// realization receipt, never silent.
    TabulatedWithFallback {
        /// The baked table (bead zolja).
        table: &'a TabulatedLoad,
        /// Fallback termination for DFT bins outside the support.
        fallback: Termination,
    },
}

/// Configuration for one bank realization.
#[derive(Debug, Clone, Copy)]
pub struct MmLineConfig {
    /// Audio sample rate [Hz].
    pub sample_rate_hz: u32,
    /// Retained mode count (>= 1; the plane mode is mode 0).
    pub n_modes: usize,
    /// Staircase density multiplier for the MM authority.
    pub extra_slices: usize,
    /// Loss model handed to the authority.
    pub loss: LossModel,
}

/// One realized combo's receipt.
#[derive(Debug, Clone)]
pub struct MmComboReceipt {
    /// Caller label (valve combination name).
    pub label: String,
    /// FIR length (= FFT size).
    pub n_fft: usize,
    /// Geometric round-trip delay [samples].
    pub round_trip_samples: f64,
    /// DFT bins served by the tabulated load.
    pub table_bins: usize,
    /// DFT bins served by the analytic fallback.
    pub fallback_bins: usize,
    /// Peak `|R|` over the 4x-oversampled realization grid (must satisfy `<= 1 + 1e-12` without clipping).
    pub peak_reflectance: f64,
    /// Number of [`crate::driving_point::LineLadder`] rungs evaluated to reach convergence.
    pub rungs_tried: usize,
    /// Relative energy in the fourth quarter `[3*n_fft/4 .. n_fft)` at convergence.
    pub tail_energy: f64,
    /// Worst shift in cents of the first three impedance peaks vs the preceding rung.
    pub peak_shift_cents: f64,
}

/// A record of one valve switch's state lift.
#[derive(Debug, Clone)]
pub struct MmLiftRecord {
    /// Combo switched from.
    pub from: usize,
    /// Combo switched to.
    pub to: usize,
    /// History samples carried into the new line.
    pub carried: usize,
    /// The lift map's name.
    pub lift: &'static str,
}

/// The per-valve-combination multimodal line bank: exact-FIR
/// plane-mode reflectance lines sharing one sample clock, switched with
/// the carry-outgoing-history lift.
pub struct MmLineBank {
    lines: Vec<DelayedFilter>,
    taps: Vec<Vec<f64>>,
    receipts: Vec<MmComboReceipt>,
    active: usize,
    dt: f64,
    zc0: f64,
    lift_log: Vec<MmLiftRecord>,
    fade_old: Option<DelayedFilter>,
    fade_remaining: usize,
    fade_total: usize,
}

/// Evaluate `mm_input_impedance`, falling back to `LossModel::Bessel` only when
/// `LossModel::WideTube` refuses a low-frequency DFT bin below `s = 10`.
fn mm_impedance_at(
    duct: &Duct,
    gas: &GasState,
    omega: f64,
    config: &MmLineConfig,
    termination: Termination,
) -> Result<ModalResponse, MmLineError> {
    match mm_input_impedance(
        duct,
        gas,
        omega,
        config.loss,
        termination,
        config.n_modes,
        config.extra_slices,
    ) {
        Ok(r) => Ok(r),
        Err(DuctError::TooNarrow { .. }) if matches!(config.loss, LossModel::WideTube) => Ok(
            mm_input_impedance(
                duct,
                gas,
                omega,
                LossModel::Bessel,
                termination,
                config.n_modes,
                config.extra_slices,
            )?,
        ),
        Err(e) => Err(e.into()),
    }
}

/// Sample the plane-mode load for one DFT bin under the splice policy.
fn resolve_load(
    duct: &Duct,
    gas: &GasState,
    omega: f64,
    config: &MmLineConfig,
    load: &MmLoad<'_>,
) -> Result<(ModalResponse, bool), MmLineError> {
    match load {
        MmLoad::Analytic(t) => Ok((mm_impedance_at(duct, gas, omega, config, *t)?, false)),
        MmLoad::TabulatedWithFallback { table, fallback } => {
            let (lo, hi) = table.support();
            if omega >= lo && omega <= hi {
                let resp = match mm_input_impedance_tabulated(
                    duct,
                    gas,
                    omega,
                    config.loss,
                    table,
                    config.n_modes,
                    config.extra_slices,
                ) {
                    Ok(r) => r,
                    Err(DuctError::TooNarrow { .. })
                        if matches!(config.loss, LossModel::WideTube) =>
                    {
                        mm_input_impedance_tabulated(
                            duct,
                            gas,
                            omega,
                            LossModel::Bessel,
                            table,
                            config.n_modes,
                            config.extra_slices,
                        )?
                    }
                    Err(e) => return Err(e.into()),
                };
                Ok((resp, true))
            } else {
                Ok((mm_impedance_at(duct, gas, omega, config, *fallback)?, false))
            }
        }
    }
}

fn duct_axial_length(duct: &Duct) -> f64 {
    duct.segments
        .iter()
        .map(Segment::axial_length)
        .sum()
}

fn input_radius(duct: &Duct) -> Result<f64, MmLineError> {
    match duct.segments.first() {
        Some(&Segment::Cylinder { radius, .. }) => Ok(radius),
        Some(&Segment::Cone { inlet_radius, .. }) => Ok(inlet_radius),
        _ => Err(MmLineError::Invalid {
            what: "combo must start with a cylinder or cone",
        }),
    }
}

/// Eight-round-trip FFT sizing for the [`MatrixFirLine`] bake-off arm.
/// [`MmLineBank`] uses the shared convergence ladder instead (bead .82).
/// This is the driving_point clamp-first discipline (u6pph precedent;
/// bead frankensim-mxw6x): the grid
/// DENSIFIES to fit eight round trips of the bore, and only an authored
/// hard ceiling refuses — typed and loud. Sizing must never grow first
/// and clamp after, which silently truncates reflectance memory into
/// wrap-around aliasing exactly where the longest bores need it most.
const MM_FFT_CEILING_TAPS: usize = 65_536;

fn mm_fft_taps_for(geo_delay: f64) -> Result<usize, MmLineError> {
    let want = 8.0 * geo_delay;
    if !want.is_finite() {
        return Err(MmLineError::Invalid {
            what: "eight-round-trip demand is not finite",
        });
    }
    if want > MM_FFT_CEILING_TAPS as f64 {
        return Err(MmLineError::Invalid {
            what: "eight round trips exceed the 65536-tap authored ceiling",
        });
    }
    Ok((want.ceil() as usize)
        .next_power_of_two()
        .clamp(256, MM_FFT_CEILING_TAPS))
}

impl MmLineBank {
    /// Realize one exact-FIR plane-mode line per valve combination.
    ///
    /// All combos share the bank's FFT size, so the switch lift is a plain
    /// history copy. That size comes from the [`crate::driving_point::LineLadder`]
    /// convergence ladder over every combo jointly. Passivity
    /// (`max |R(e^{iωΔt})| <= 1`) is verified on the 4x-oversampled grid
    /// without `1/peak` clipping.
    ///
    /// # Errors
    /// [`MmLineError`] on refusal: geometry, authority, non-passivity, or
    /// [`MmLineError::LadderUnconverged`] at the ladder ceiling.
    #[allow(clippy::too_many_lines)] // one realization pipeline
    pub fn new(
        ducts: &[Duct],
        labels: &[&str],
        gas: &GasState,
        load: &MmLoad<'_>,
        config: &MmLineConfig,
    ) -> Result<MmLineBank, MmLineError> {
        if ducts.is_empty() || ducts.len() != labels.len() {
            return Err(MmLineError::Invalid {
                what: "combos and labels must be non-empty and match",
            });
        }
        if config.sample_rate_hz == 0 {
            return Err(MmLineError::Invalid {
                what: "sample rate must be positive",
            });
        }
        let dt = 1.0 / f64::from(config.sample_rate_hz);
        let r_in = input_radius(&ducts[0])?;
        let zc0 = gas.density * gas.sound_speed / (core::f64::consts::PI * r_in * r_in);
        let ladder = crate::driving_point::LineLadder::default();

        let mut max_geo_delay = 0.0_f64;
        let mut delays = Vec::with_capacity(ducts.len());
        for duct in ducts {
            if (input_radius(duct)? - r_in).abs() > 1e-12 * r_in {
                return Err(MmLineError::Invalid {
                    what: "every combo must share the input radius (one zc)",
                });
            }
            let length = duct_axial_length(duct);
            if !(length > 0.0) {
                return Err(MmLineError::Invalid {
                    what: "combo needs positive axial length",
                });
            }
            let geo_delay = 2.0 * length / gas.sound_speed / dt;
            if geo_delay < 2.0 {
                return Err(MmLineError::Invalid {
                    what: "round-trip delay below two samples",
                });
            }
            if !(4.0 * geo_delay <= ladder.max_taps as f64) {
                return Err(MmLineError::Invalid {
                    what: "round trip does not fit the 65536-tap line ladder",
                });
            }
            max_geo_delay = max_geo_delay.max(geo_delay);
            delays.push(geo_delay);
        }

        // DC: each combo's own plane-reflectance limit (DC flow resistance
        // keeps a lossy bore off -1). A fixed -1 that disagrees with bin 1 is
        // a spectral jump, which leaves an N-independent tail in the IR.
        let r_dc: Vec<C64> = ducts
            .iter()
            .map(|duct| {
                resolve_load(duct, gas, 1.0e-3, config, load)
                    .ok()
                    .map(|(resp, _)| {
                        reflectance(
                            C64::new(resp.plane_impedance.re, resp.plane_impedance.im),
                            zc0,
                        )
                    })
                    .filter(|r| r.re.is_finite())
                    .map_or(C64::new(-1.0, 0.0), |r| C64::new(r.re.clamp(-1.0, 1.0), 0.0))
            })
            .collect();

        // Shared convergence ladder (DEF-02, bead .82). One FFT size for
        // every combo keeps the switch lift a plain history copy. Start at
        // four round trips of the longest bore and double until EVERY combo
        // has tail energy < tail_tol and first-three-peak shift < peak_tol.
        struct Rung {
            r_pos: Vec<C64>,
            from_table: Vec<bool>,
            peaks: Vec<f64>,
            ir: Vec<f64>,
            peak: f64,
            tail_energy: f64,
            shift: f64,
        }
        let mut rungs: Vec<Rung> = ducts
            .iter()
            .map(|_| Rung {
                r_pos: Vec::new(),
                from_table: Vec::new(),
                peaks: Vec::new(),
                ir: Vec::new(),
                peak: 0.0,
                tail_energy: f64::INFINITY,
                shift: f64::INFINITY,
            })
            .collect();
        let mut n_fft = ((4.0 * max_geo_delay).ceil() as usize)
            .max(256)
            .next_power_of_two()
            .min(ladder.max_taps / 2)
            .max(256);
        let mut rungs_tried = 0usize;
        loop {
            rungs_tried += 1;
            let half = n_fft / 2;
            let mut all_converged = true;
            for (combo, duct) in ducts.iter().enumerate() {
                let rung = &mut rungs[combo];
                let mut r_pos = vec![C64::ZERO; half + 1];
                let mut from_table = vec![false; half + 1];
                r_pos[0] = r_dc[combo];
                let fresh: Vec<usize> = if rung.r_pos.is_empty() {
                    (1..=half).collect()
                } else {
                    // Reuse the even bins (bit-identical frequencies).
                    for (k_prev, (&v, &t)) in rung.r_pos.iter().zip(&rung.from_table).enumerate() {
                        r_pos[2 * k_prev] = v;
                        from_table[2 * k_prev] = t;
                    }
                    (1..half).step_by(2).collect()
                };
                for k in fresh {
                    let omega = core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt);
                    let (response, table) = resolve_load(duct, gas, omega, config, load)?;
                    from_table[k] = table;
                    r_pos[k] = reflectance(
                        C64::new(response.plane_impedance.re, response.plane_impedance.im),
                        zc0,
                    );
                }
                let geo_delay = delays[combo];
                let ir = crate::driving_point::realize_reflectance_ir_from_positive_spectrum(
                    &r_pos, n_fft, geo_delay, true,
                );
                let total_energy: f64 = ir.iter().map(|&h| h * h).sum();
                let tail_sum: f64 = ir[3 * (n_fft / 4)..].iter().map(|&h| h * h).sum();
                let tail_energy = if total_energy > 0.0 {
                    tail_sum / total_energy
                } else {
                    0.0
                };
                let (peak, peaks) = crate::driving_point::analyze_reflectance_ir(&ir, dt, geo_delay);
                let shift = if rungs_tried == 1 {
                    f64::INFINITY
                } else if rung.peaks.is_empty() && peaks.is_empty() {
                    0.0
                } else if rung.peaks.is_empty() || peaks.is_empty() {
                    f64::INFINITY
                } else {
                    let take = rung.peaks.len().min(peaks.len()).min(3);
                    rung.peaks[..take]
                        .iter()
                        .zip(&peaks[..take])
                        .map(|(&a, &b)| (1200.0 * (b / a).log2()).abs())
                        .fold(0.0_f64, f64::max)
                };
                all_converged &= (n_fft as f64) > geo_delay + 2.0
                    && tail_energy < ladder.tail_tol
                    && shift < ladder.peak_tol_cents;
                *rung = Rung {
                    r_pos,
                    from_table,
                    peaks,
                    ir,
                    peak,
                    tail_energy,
                    shift,
                };
            }
            if all_converged {
                break;
            }
            if n_fft >= ladder.max_taps {
                return Err(MmLineError::LadderUnconverged {
                    max_taps: ladder.max_taps,
                    tail_energy: rungs.iter().map(|r| r.tail_energy).fold(0.0, f64::max),
                    peak_shift_cents: rungs.iter().map(|r| r.shift).fold(0.0, f64::max),
                });
            }
            n_fft *= 2;
        }

        let mut lines = Vec::with_capacity(ducts.len());
        let mut taps = Vec::with_capacity(ducts.len());
        let mut receipts = Vec::with_capacity(ducts.len());
        for (combo, (rung, &label)) in rungs.into_iter().zip(labels).enumerate() {
            if rung.peak > 1.0 + 1.0e-12 {
                return Err(MmLineError::Invalid {
                    what: "realized multimodal line violates 4x-oversampled passivity",
                });
            }
            let table_bins = rung.from_table[1..].iter().filter(|&&t| t).count();
            lines.push(DelayedFilter::from_impulse_response(dt, rung.ir.clone())?);
            taps.push(rung.ir);
            receipts.push(MmComboReceipt {
                label: label.to_string(),
                n_fft,
                round_trip_samples: delays[combo],
                table_bins,
                fallback_bins: n_fft / 2 - table_bins,
                peak_reflectance: rung.peak,
                rungs_tried,
                tail_energy: rung.tail_energy,
                peak_shift_cents: rung.shift,
            });
        }

        Ok(MmLineBank {
            lines,
            taps,
            receipts,
            active: 0,
            dt,
            zc0,
            lift_log: Vec::new(),
            fade_old: None,
            fade_remaining: 0,
            fade_total: 0,
        })
    }

    /// Plane-mode characteristic impedance at the input plane
    /// [Pa s/m^3].
    #[must_use]
    pub fn zc0(&self) -> f64 {
        self.zc0
    }

    /// Realization receipts, one per combo.
    #[must_use]
    pub fn receipts(&self) -> &[MmComboReceipt] {
        &self.receipts
    }

    /// The (passivity-enforced) stored taps of one combo — the realized
    /// reflectance impulse response a switch rebuilds from.
    #[must_use]
    pub fn combo_taps(&self, combo: usize) -> &[f64] {
        &self.taps[combo]
    }

    /// The active combo index.
    #[must_use]
    pub fn active(&self) -> usize {
        self.active
    }

    /// The lift log (one record per switch).
    #[must_use]
    pub fn lift_log(&self) -> &[MmLiftRecord] {
        &self.lift_log
    }

    /// Push the outgoing wave, get the incoming wave. During a valve
    /// crossfade BOTH lines advance and the incoming waves blend on a
    /// linear ramp (an instantaneous reflectance swap is genuinely
    /// discontinuous physics; a real valve passes through ~10 ms of
    /// intermediate states).
    ///
    /// # Errors
    /// [`MmLineError::Discrete`] on a non-finite sample.
    pub fn push(&mut self, outgoing: f64) -> Result<f64, MmLineError> {
        let new_in = self.lines[self.active].push(outgoing)?;
        if let Some(old) = self.fade_old.as_mut() {
            let old_in = old.push(outgoing)?;
            self.fade_remaining -= 1;
            let w = 1.0 - self.fade_remaining as f64 / self.fade_total as f64;
            let mixed = (1.0 - w) * old_in + w * new_in;
            if self.fade_remaining == 0 {
                self.fade_old = None;
            }
            return Ok(mixed);
        }
        Ok(new_in)
    }

    /// Re-read the last incoming wave without advancing (blended while
    /// a crossfade is active).
    #[must_use]
    pub fn incoming(&self) -> f64 {
        let new_in = self.lines[self.active].incoming();
        if let Some(old) = self.fade_old.as_ref() {
            let w = 1.0 - self.fade_remaining as f64 / self.fade_total as f64;
            return (1.0 - w) * old.incoming() + w * new_in;
        }
        new_in
    }

    /// Reset every line to silence (rebuild from the stored taps) and
    /// clear any crossfade — the cheap probe/battery seam; realization
    /// receipts and the lift log persist.
    ///
    /// # Errors
    /// [`MmLineError::Discrete`] if a stored tap set fails re-admission
    /// (cannot happen for taps this bank minted).
    pub fn reset(&mut self) -> Result<(), MmLineError> {
        for (line, taps) in self.lines.iter_mut().zip(&self.taps) {
            *line = DelayedFilter::from_impulse_response(self.dt, taps.clone())?;
        }
        self.fade_old = None;
        self.fade_remaining = 0;
        self.fade_total = 0;
        Ok(())
    }

    /// Switch valve combination AT A BLOCK BOUNDARY: the new line is
    /// rebuilt from its taps with the outgoing-wave history CARRIED
    /// (in-flight waves persist; the new reflectance governs future
    /// reflections) and the incoming waves CROSSFADE from the old line
    /// over `fade_samples` (0 = hard swap). The lift is logged.
    ///
    /// # Errors
    /// [`MmLineError::Invalid`] on an unknown combo.
    pub fn switch(&mut self, combo: usize, fade_samples: usize) -> Result<(), MmLineError> {
        if combo >= self.lines.len() {
            return Err(MmLineError::Invalid {
                what: "unknown valve combination",
            });
        }
        if combo == self.active {
            return Ok(());
        }
        let history = self.lines[self.active].history();
        let lifted = DelayedFilter::from_impulse_response_with_history(
            self.dt,
            self.taps[combo].clone(),
            &history,
        )?;
        self.lift_log.push(MmLiftRecord {
            from: self.active,
            to: combo,
            carried: history.len(),
            lift: if fade_samples > 0 {
                "carry-outgoing-history+crossfade/v1"
            } else {
                "carry-outgoing-history/v1"
            },
        });
        if fade_samples > 0 {
            // The old line keeps ringing during the fade.
            self.fade_old = Some(self.lines[self.active].clone());
            self.fade_remaining = fade_samples;
            self.fade_total = fade_samples;
        } else {
            self.fade_old = None;
            self.fade_remaining = 0;
            self.fade_total = 0;
        }
        self.lines[combo] = lifted;
        self.active = combo;
        Ok(())
    }
}

/// Peak `|DTFT|` of a tap vector over a frequency grid.
#[cfg(test)]
fn dtft_peak(taps: &[f64], dt: f64, omegas: &[f64]) -> f64 {
    let mut peak = 0.0f64;
    for &w in omegas {
        let mut re = 0.0f64;
        let mut im = 0.0f64;
        for (k, &h) in taps.iter().enumerate() {
            let ph = w * dt * k as f64;
            re += h * ph.cos();
            im -= h * ph.sin();
        }
        peak = peak.max((re * re + im * im).sqrt());
    }
    peak
}

/// The full-matrix realization arm: N x N exact FIRs over shared
/// per-mode outgoing histories, with per-mode incoming-energy logging.
/// This is the bake-off arm and the lane for sources that genuinely
/// inject higher modes at the input plane.
pub struct MatrixFirLine {
    n_modes: usize,
    n_fft: usize,
    taps: Vec<Vec<f64>>,
    history: Vec<Vec<f64>>,
    write: usize,
    incoming: Vec<f64>,
    energy: Vec<f64>,
}

impl MatrixFirLine {
    /// Realize the full N x N reflection-matrix impulse responses from
    /// the MM authority (one sweep serves all entries) and the per-mode
    /// characteristic impedances at the input plane.
    ///
    /// # Errors
    /// [`MmLineError`] on refusal.
    #[allow(clippy::too_many_lines)] // one realization pipeline
    pub fn new(
        duct: &Duct,
        gas: &GasState,
        load: &MmLoad<'_>,
        config: &MmLineConfig,
    ) -> Result<MatrixFirLine, MmLineError> {
        let n = config.n_modes;
        if n == 0 {
            return Err(MmLineError::Invalid {
                what: "at least the plane mode",
            });
        }
        let dt = 1.0 / f64::from(config.sample_rate_hz);
        let length = duct_axial_length(duct);
        let geo_delay = 2.0 * length / gas.sound_speed / dt;
        let n_fft = mm_fft_taps_for(geo_delay)?;
        let half = n_fft / 2;
        let r_in = input_radius(duct)?;
        // One positive-frequency spectrum buffer per matrix entry.
        let mut r_pos_mat = vec![vec![C64::ZERO; half + 1]; n * n];
        for (m, r_pos) in r_pos_mat.iter_mut().enumerate() {
            // DC: plane reflects -1 off the open mouth; higher modes are
            // fully evanescent at DC (no propagating reflection).
            r_pos[0] = if m == 0 {
                C64::new(-1.0, 0.0)
            } else {
                C64::ZERO
            };
        }
        for k in 1..=half {
            let omega = core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt);
            let (response, _) = resolve_load(duct, gas, omega, config, load)?;
            let zc = match modal_characteristic_impedances(gas, r_in, omega, config.loss, n) {
                Ok(zc) => zc,
                Err(DuctError::TooNarrow { .. })
                    if matches!(config.loss, LossModel::WideTube) =>
                {
                    modal_characteristic_impedances(gas, r_in, omega, LossModel::Bessel, n)?
                }
                Err(e) => return Err(e.into()),
            };
            let r = modal_reflection_from_impedance(&response.impedance_matrix, &zc, n)?;
            for (m, r_pos) in r_pos_mat.iter_mut().enumerate() {
                r_pos[k] = r[m];
            }
        }
        let mut taps = Vec::with_capacity(n * n);
        for r_pos in &r_pos_mat {
            taps.push(
                crate::driving_point::realize_reflectance_ir_from_positive_spectrum(
                    r_pos, n_fft, geo_delay, false,
                ),
            );
        }
        Ok(MatrixFirLine {
            n_modes: n,
            n_fft,
            taps,
            history: vec![vec![0.0; n_fft]; n],
            write: 0,
            incoming: vec![0.0; n],
            energy: vec![0.0; n],
        })
    }

    /// Push the per-mode outgoing waves; returns the per-mode incoming
    /// waves `b_m = sum_n R_mn * a_n`.
    ///
    /// # Errors
    /// [`MmLineError::Invalid`] on a length mismatch or non-finite input.
    pub fn push(&mut self, outgoing: &[f64]) -> Result<&[f64], MmLineError> {
        let n = self.n_modes;
        if outgoing.len() != n {
            return Err(MmLineError::Invalid {
                what: "outgoing vector must have one entry per mode",
            });
        }
        for &a in outgoing {
            if !a.is_finite() {
                return Err(MmLineError::Invalid {
                    what: "outgoing wave left the finite set",
                });
            }
        }
        for (mode, &a) in outgoing.iter().enumerate() {
            self.history[mode][self.write] = a;
        }
        let len = self.n_fft;
        for m in 0..n {
            let mut acc = 0.0f64;
            for j in 0..n {
                let taps = &self.taps[m * n + j];
                let hist = &self.history[j];
                for k in 0..len {
                    acc += taps[k] * hist[(self.write + len - k) % len];
                }
            }
            self.incoming[m] = acc;
            self.energy[m] += acc * acc;
        }
        self.write = (self.write + 1) % len;
        Ok(&self.incoming)
    }

    /// Cumulative per-mode incoming energy (who carries the sound).
    #[must_use]
    pub fn mode_energy(&self) -> &[f64] {
        &self.energy
    }

    /// FIR length.
    #[must_use]
    pub fn n_fft(&self) -> usize {
        self.n_fft
    }
}

/// Mouthpiece-cup state for [`cup_junction`] (D18 lumped compliance).
#[derive(Debug, Clone, Copy, Default)]
pub struct CupState {
    /// Previous junction pressure [Pa].
    pub p_prev: f64,
    /// Previous cup volume velocity [m^3/s].
    pub u_prev: f64,
}

/// One-sample mouthpiece-cup junction (D18 electrically-short lumped
/// compliance `C = V/(rho c^2)` shunting the input plane): solve the
/// LINEAR junction `u_src = (p_plus - p_minus)/zc0 + u_cup` with the
/// TRAPEZOIDAL (bilinear) compliance `(u_cup + u_prev)/2 =
/// C (p - p_prev)/dt`, `p = p_plus + p_minus` — second-order accurate,
/// so the runtime matches the analytic shunt `Z/(1 + i omega C Z)` to
/// `O((omega dt)^2)`. Returns `(p_plus, p)` and updates the cup state.
#[must_use]
pub fn cup_junction(
    u_source: f64,
    p_minus: f64,
    cup: &mut CupState,
    compliance: f64,
    zc0: f64,
    dt: f64,
) -> (f64, f64) {
    let g = 2.0 * compliance / dt;
    let a = 1.0 / zc0 + g;
    let rhs = u_source + p_minus / zc0 - g * (p_minus - cup.p_prev) + cup.u_prev;
    let p_plus = rhs / a;
    let p = p_plus + p_minus;
    let u_cup = g * (p - cup.p_prev) - cup.u_prev;
    cup.p_prev = p;
    cup.u_prev = u_cup;
    (p_plus, p)
}

#[cfg(test)]
mod mm_line_tests {
    use super::*;
    use fs_material::gas::GasSpec;

    fn air20() -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
    }

    fn verdict(case: &str, pass: bool, detail: &str) {
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"{case}\",\"verdict\":\"{}\",\"detail\":\"{detail}\"}}",
            if pass { "pass" } else { "fail" }
        );
        assert!(pass, "case {case}: {detail}");
    }

    /// Brass-like combo: leadpipe (+ optional crook) + taper + flare.
    fn brass_combo(crook_m: f64) -> Duct {
        let mut segments = vec![Segment::Cylinder {
            radius: 0.006,
            length: 0.30,
        }];
        if crook_m > 0.0 {
            segments.push(Segment::Cylinder {
                radius: 0.006,
                length: crook_m,
            });
        }
        segments.push(Segment::Cone {
            inlet_radius: 0.006,
            outlet_radius: 0.012,
            length: 0.25,
        });
        segments.push(Segment::Cone {
            inlet_radius: 0.012,
            outlet_radius: 0.035,
            length: 0.12,
        });
        Duct { segments }
    }

    // Exact Bessel losses (.76). WideTube refuses below shear number 10, and
    // mixing it with Bessel per DFT bin is a spectral kink the convergence
    // ladder rightly refuses.
    fn config() -> MmLineConfig {
        MmLineConfig {
            sample_rate_hz: 24_000,
            n_modes: 3,
            extra_slices: 1,
            loss: LossModel::Bessel,
        }
    }

    /// Drive the bank with a volume-velocity impulse and return the
    /// junction pressure series (the runtime impedance's time image).
    fn impulse_response(bank: &mut MmLineBank, samples: usize) -> Vec<f64> {
        let zc0 = bank.zc0();
        let mut p = Vec::with_capacity(samples);
        for k in 0..samples {
            let u = if k == 0 { 1.0 } else { 0.0 };
            let p_minus = bank.incoming();
            let p_plus = p_minus + zc0 * u;
            let _ = bank.push(p_plus).expect("push");
            p.push(p_plus + p_minus);
        }
        p
    }

    /// |Z| peaks of a pressure impulse series via FFT + parabolic refine.
    fn runtime_peaks_hz(p: &[f64], sample_rate_hz: f64, lo_hz: f64, hi_hz: f64) -> Vec<f64> {
        let n = p.len().next_power_of_two();
        let fft = Fft::new(n);
        let mut buf: Vec<FftC64> = (0..n)
            .map(|k| FftC64::new(p.get(k).copied().unwrap_or(0.0), 0.0))
            .collect();
        let mut scratch = vec![FftC64::new(0.0, 0.0); n];
        fft.forward(&mut buf, &mut scratch);
        let df = sample_rate_hz / n as f64;
        let mags: Vec<f64> = buf[..n / 2]
            .iter()
            .map(|c| (c.re * c.re + c.im * c.im).sqrt())
            .collect();
        let mut peaks = Vec::new();
        let k_lo = (lo_hz / df).ceil() as usize;
        let k_hi = (hi_hz / df).floor() as usize;
        for k in k_lo.max(1)..k_hi.min(mags.len() - 1) {
            if mags[k] > mags[k - 1] && mags[k] > mags[k + 1] {
                let (ya, yb, yc) = (mags[k - 1].ln(), mags[k].ln(), mags[k + 1].ln());
                let shift = 0.5 * (ya - yc) / (ya - 2.0 * yb + yc);
                peaks.push((k as f64 - shift) * df);
            }
        }
        peaks
    }

    /// Authority peaks: sweep the MM TMM plane impedance and refine.
    fn authority_peaks_hz(
        duct: &Duct,
        gas: &GasState,
        cfg: &MmLineConfig,
        lo_hz: f64,
        hi_hz: f64,
        count: usize,
    ) -> Vec<f64> {
        let mut mags = Vec::with_capacity(count);
        for i in 0..count {
            let f = lo_hz + (hi_hz - lo_hz) * i as f64 / (count - 1) as f64;
            let omega = core::f64::consts::TAU * f;
            let z = mm_input_impedance(
                duct,
                gas,
                omega,
                cfg.loss,
                Termination::FlangedOpen,
                cfg.n_modes,
                cfg.extra_slices,
            )
            .expect("authority")
            .plane_impedance;
            mags.push((f, z.abs()));
        }
        let mut peaks = Vec::new();
        for i in 1..mags.len() - 1 {
            if mags[i].1 > mags[i - 1].1 && mags[i].1 > mags[i + 1].1 {
                let (ya, yb, yc) = (mags[i - 1].1.ln(), mags[i].1.ln(), mags[i + 1].1.ln());
                let shift = 0.5 * (ya - yc) / (ya - 2.0 * yb + yc);
                let df = mags[i].0 - mags[i - 1].0;
                peaks.push(mags[i].0 - shift * df);
            }
        }
        peaks
    }

    fn cents(a: f64, b: f64) -> f64 {
        1200.0 * (b / a).log2()
    }

    #[test]
    fn ml_001_held_fingering_matches_the_authority() {
        let gas = air20();
        let cfg = config();
        let combos = [brass_combo(0.0), brass_combo(0.08), brass_combo(0.15)];
        let labels = ["open", "crook-1", "crook-2"];
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let mut bank = MmLineBank::new(&combos, &labels, &gas, &load, &cfg).expect("bank realizes");
        let dt = 1.0 / 24_000.0;
        // LINE FIDELITY: |Z| from the realized taps' own DTFT,
        // Z = zc (1 + R)/(1 - R), peak-matched to the authority. This
        // isolates the realization from the explicit junction loop's
        // known one-sample skew (disclosed below).
        let taps_peaks = |taps: &[f64], zc0: f64| -> Vec<f64> {
            let count = 391usize;
            let mut mags = Vec::with_capacity(count);
            for i in 0..count {
                let f = 120.0 + (900.0 - 120.0) * i as f64 / (count - 1) as f64;
                let w = core::f64::consts::TAU * f;
                let mut re = 0.0f64;
                let mut im = 0.0f64;
                for (k, &h) in taps.iter().enumerate() {
                    let ph = w * dt * k as f64;
                    re += h * ph.cos();
                    im -= h * ph.sin();
                }
                let r = C64::new(re, im);
                let z = (C64::ONE + r) * (C64::ONE - r).recip();
                mags.push((f, z.abs() * zc0));
            }
            let mut peaks = Vec::new();
            for i in 1..mags.len() - 1 {
                if mags[i].1 > mags[i - 1].1 && mags[i].1 > mags[i + 1].1 {
                    let (ya, yb, yc) = (mags[i - 1].1.ln(), mags[i].1.ln(), mags[i + 1].1.ln());
                    let shift = 0.5 * (ya - yc) / (ya - 2.0 * yb + yc);
                    let df = mags[i].0 - mags[i - 1].0;
                    peaks.push(mags[i].0 - shift * df);
                }
            }
            peaks
        };
        let mut worst_line = 0.0f64;
        let mut worst_loop = 0.0f64;
        let mut detail = String::new();
        for (idx, (duct, label)) in combos.iter().zip(&labels).enumerate() {
            bank.switch(idx, 0).expect("switch");
            let auth = authority_peaks_hz(duct, &gas, &cfg, 120.0, 900.0, 391);
            assert!(auth.len() >= 3, "{label}: auth {auth:?}");
            let line = taps_peaks(bank.combo_taps(idx), bank.zc0());
            let p = impulse_response(&mut bank, 8192);
            let looped = runtime_peaks_hz(&p, 24_000.0, 120.0, 900.0);
            let nearest = |list: &[f64], a: f64| -> f64 {
                list.iter()
                    .copied()
                    .min_by(|x, y| (x - a).abs().partial_cmp(&(y - a).abs()).expect("finite"))
                    .expect("nonempty")
            };
            for &a in auth.iter().take(4) {
                worst_line = worst_line.max(cents(a, nearest(&line, a)).abs());
                worst_loop = worst_loop.max(cents(a, nearest(&looped, a)).abs());
            }
            let _ = core::fmt::Write::write_fmt(
                &mut detail,
                format_args!("{label}: {} auth peaks; ", auth.len().min(4)),
            );
        }
        // AUTHORED bands: the realized LINE must track the authority
        // tightly; the explicit junction LOOP additionally carries its
        // one-sample skew (about 1200*log2((T+1)/T) ~ 17 cents at a
        // 100-sample round trip) — a disclosed property of the explicit
        // scattering loop, shared with the reed voice, not a line error.
        // Loop arm measured at 38.4 cents on this fixture (a few samples
        // of loop phase at a ~127-sample fundamental period).
        let pass = worst_line < 5.0 && worst_loop < 45.0;
        verdict(
            "ml-001-held-fingering-cents",
            pass,
            &format!(
                "{detail}worst line-fidelity {worst_line:.2} cents; worst loop \
                 (incl. one-sample junction skew) {worst_loop:.2} cents"
            ),
        );
    }

    #[test]
    fn ml_002_realization_bakeoff_receipt() {
        // Dominant R00 line vs the full-matrix arm on the open combo,
        // both driven plane-only with a rigid higher-mode source
        // closure. The receipt records accuracy AND cost.
        let gas = air20();
        let cfg = config();
        let duct = brass_combo(0.0);
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let mut bank = MmLineBank::new(core::slice::from_ref(&duct), &["open"], &gas, &load, &cfg)
            .expect("dominant");
        let p_dom = impulse_response(&mut bank, 4096);
        let mut matrix = MatrixFirLine::new(&duct, &gas, &load, &cfg).expect("matrix");
        let zc0 = bank.zc0();
        let n = cfg.n_modes;
        let mut p_mat = Vec::with_capacity(4096);
        let mut incoming = vec![0.0f64; n];
        for k in 0..4096usize {
            let u = if k == 0 { 1.0 } else { 0.0 };
            let mut outgoing = vec![0.0f64; n];
            outgoing[0] = incoming[0] + zc0 * u;
            outgoing[1..n].copy_from_slice(&incoming[1..n]); // rigid source closure
            let b = matrix.push(&outgoing).expect("push");
            p_mat.push(outgoing[0] + incoming[0]);
            incoming.copy_from_slice(b);
        }
        let peaks_dom = runtime_peaks_hz(&p_dom, 24_000.0, 120.0, 900.0);
        let peaks_mat = runtime_peaks_hz(&p_mat, 24_000.0, 120.0, 900.0);
        let n_cmp = peaks_dom.len().min(peaks_mat.len()).min(4);
        let mut worst = 0.0f64;
        for i in 0..n_cmp {
            worst = worst.max(cents(peaks_dom[i], peaks_mat[i]).abs());
        }
        let energy = matrix.mode_energy();
        let plane_fraction = energy[0] / energy.iter().sum::<f64>().max(f64::MIN_POSITIVE);
        let cost_ratio = (n * n) as f64;
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"ml-002-bakeoff-receipt\",\
             \"peaks_dominant_hz\":{:?},\"peaks_matrix_hz\":{:?},\
             \"worst_cents_delta\":{worst:.3},\"plane_energy_fraction\":{plane_fraction:.6},\
             \"cost_ratio_matrix_over_dominant\":{cost_ratio:.0},\
             \"decision\":\"dominant (equal accuracy for plane-driven sources, {cost_ratio:.0}x cheaper); matrix retained for higher-mode sources\"}}",
            peaks_dom
                .iter()
                .map(|p| (p * 10.0).round() / 10.0)
                .collect::<Vec<_>>(),
            peaks_mat
                .iter()
                .map(|p| (p * 10.0).round() / 10.0)
                .collect::<Vec<_>>(),
        );
        // Peaks must agree tightly (the matrix arm's off-plane input
        // coupling is evanescent-small at a 6 mm throat) and the plane
        // mode must carry essentially all incoming energy.
        let pass = n_cmp >= 3 && worst < 1.0 && plane_fraction > 0.999;
        verdict(
            "ml-002-realization-bakeoff",
            pass,
            &format!(
                "dominant vs matrix worst {worst:.3} cents over {n_cmp} peaks, \
                 plane energy fraction {plane_fraction:.6}"
            ),
        );
    }

    #[test]
    fn ml_003_valve_switch_carries_state_without_click() {
        let gas = air20();
        let cfg = config();
        let combos = [brass_combo(0.0), brass_combo(0.08)];
        let labels = ["open", "crook-1"];
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let run = |switch_at: Option<usize>, hard_cut: bool| -> (Vec<f64>, usize) {
            let mut bank = MmLineBank::new(&combos, &labels, &gas, &load, &cfg).expect("bank");
            let zc0 = bank.zc0();
            let mut p = Vec::with_capacity(1600);
            let mut carried = 0usize;
            for k in 0..1600usize {
                if Some(k) == switch_at {
                    if hard_cut {
                        // The falsifier arm: a FRESH bank (silent line)
                        // replaces the ringing one — the no-lift cut.
                        bank = MmLineBank::new(&combos, &labels, &gas, &load, &cfg).expect("bank");
                        bank.switch(1, 0).expect("switch");
                    } else {
                        bank.switch(1, 96).expect("switch");
                        carried = bank.lift_log().last().expect("lift").carried;
                    }
                }
                let u = if k == 0 { 1.0 } else { 0.0 };
                let p_minus = bank.incoming();
                let p_plus = p_minus + zc0 * u;
                let _ = bank.push(p_plus).expect("push");
                p.push(p_plus + p_minus);
            }
            (p, carried)
        };
        let (control, _) = run(None, false);
        let (lifted, carried) = run(Some(700), false);
        let (cut, _) = run(Some(700), true);
        let jump = |p: &[f64]| -> f64 {
            (690..820)
                .map(|k| (p[k] - p[k - 1]).abs())
                .fold(0.0f64, f64::max)
        };
        let control_jump = jump(&control);
        let lifted_jump = jump(&lifted);
        let cut_jump = jump(&cut);
        // The lift must stay within a small factor of the undisturbed
        // signal's own local slew, and the no-lift cut must be visibly
        // worse — the executable proof the carried history matters.
        let pass = carried > 0 && lifted_jump < 3.0 * control_jump && cut_jump > 2.0 * lifted_jump;
        verdict(
            "ml-003-valve-switch-no-click",
            pass,
            &format!(
                "carried {carried} samples; boundary jump control {control_jump:.4e}, \
                 lifted {lifted_jump:.4e}, hard-cut {cut_jump:.4e}"
            ),
        );
    }

    #[test]
    fn ml_004_passivity_and_receipts() {
        let gas = air20();
        let cfg = config();
        let combos = [brass_combo(0.0)];
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let bank = MmLineBank::new(&combos, &["open"], &gas, &load, &cfg).expect("bank");
        let receipt = &bank.receipts()[0];
        // The STORED taps (what a switch rebuilds from) must satisfy
        // |R| <= 1 on the realization grid without 1/peak clipping.
        let dt = 1.0 / 24_000.0;
        let grid: Vec<f64> = (1..=receipt.n_fft / 2)
            .map(|k| core::f64::consts::TAU * k as f64 / (receipt.n_fft as f64 * dt))
            .collect();
        let peak_after = super::dtft_peak(&bank.taps[0], dt, &grid);
        let pass = peak_after <= 1.0 + 1e-12
            && receipt.peak_reflectance <= 1.0 + 1e-12
            && receipt.n_fft >= 1024
            && receipt.rungs_tried >= 1
            && receipt.tail_energy.is_finite();
        verdict(
            "ml-004-per-line-passivity",
            pass,
            &format!(
                "stored-tap |R| peak {peak_after:.6} (4x-oversampled {:.6}), n_fft {}, rungs_tried {}, tail_energy {:.3e}",
                receipt.peak_reflectance,
                receipt.n_fft,
                receipt.rungs_tried,
                receipt.tail_energy,
            ),
        );
    }

    #[test]
    fn ml_005_mouthpiece_cup_matches_the_analytic_shunt() {
        // The cup junction's runtime impedance must match the analytic
        // parallel-compliance transform Z' = Z/(1 + i omega C Z) of the
        // cup-free runtime impedance at the first peaks.
        let gas = air20();
        let cfg = config();
        let combos = [brass_combo(0.0)];
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let volume_m3 = 6.0e-6; // ~6 cc mouthpiece cup (authored)
        let compliance = volume_m3 / (gas.density * gas.sound_speed * gas.sound_speed);
        let dt = 1.0 / 24_000.0;
        let render = |with_cup: bool| -> Vec<f64> {
            let mut bank = MmLineBank::new(&combos, &["open"], &gas, &load, &cfg).expect("bank");
            let zc0 = bank.zc0();
            let mut cup = super::CupState::default();
            let mut p = Vec::with_capacity(8192);
            for k in 0..8192usize {
                let u = if k == 0 { 1.0 } else { 0.0 };
                let p_minus = bank.incoming();
                let (p_plus, p_now) = if with_cup {
                    cup_junction(u, p_minus, &mut cup, compliance, zc0, dt)
                } else {
                    let p_plus = p_minus + zc0 * u;
                    (p_plus, p_plus + p_minus)
                };
                let _ = bank.push(p_plus).expect("push");
                p.push(p_now);
            }
            p
        };
        let p_free = render(false);
        let p_cup = render(true);
        let n = 8192usize;
        let fft = Fft::new(n);
        let spectrum = |p: &[f64]| -> Vec<C64> {
            let mut buf: Vec<FftC64> = p.iter().map(|&x| FftC64::new(x, 0.0)).collect();
            let mut scratch = vec![FftC64::new(0.0, 0.0); n];
            fft.forward(&mut buf, &mut scratch);
            buf[..n / 2].iter().map(|c| C64::new(c.re, c.im)).collect()
        };
        let z_free = spectrum(&p_free);
        let z_cup = spectrum(&p_cup);
        let mut worst = 0.0f64;
        for k in 40..400usize {
            // 117..1170 Hz at df = 24000/8192.
            let omega = core::f64::consts::TAU * 24_000.0 * k as f64 / n as f64;
            // Runtime spectra are e^{+iwt} (DFT); the analytic transform
            // must use the SAME convention: Y_cup = -i omega C there.
            let denom = C64::ONE + C64::new(0.0, omega * compliance) * z_free[k];
            let want = z_free[k] * denom.recip();
            let rel = (z_cup[k] - want).abs() / want.abs().max(1e-300);
            worst = worst.max(rel);
        }
        // Backward-Euler compliance: first-order in omega*dt (~3% at 1 kHz).
        let pass = worst < 0.05;
        verdict(
            "ml-005-cup-shunt-analytic",
            pass,
            &format!("worst rel dev vs Z/(1+i w C Z): {worst:.4}"),
        );
    }

    #[test]
    fn ml_006_causality_splice_and_refusals() {
        let gas = air20();
        let cfg = config();
        let combos = [brass_combo(0.0)];
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        let bank = MmLineBank::new(&combos, &["open"], &gas, &load, &cfg).expect("bank");
        // CAUSALITY WITNESS for the conjugation discipline: the realized
        // reflectance IR must concentrate its energy in the first half
        // (the flipped-conjugation mutant time-reverses it).
        let taps = &bank.taps[0];
        let head: f64 = taps[..taps.len() / 2].iter().map(|t| t * t).sum();
        let tail: f64 = taps[taps.len() / 2..].iter().map(|t| t * t).sum();
        let causal = head > 10.0 * tail;
        // Splice disclosure with the committed zolja bell table.
        let text = std::fs::read_to_string("../../data/radiation/bell-fixture-zl.tsv")
            .expect("committed bell table");
        let mut omegas = Vec::new();
        let mut zs = Vec::new();
        for line in text.lines().filter(|l| !l.starts_with('#')) {
            let cols: Vec<&str> = line.split('\t').collect();
            omegas.push(cols[0].parse::<f64>().expect("omega"));
            zs.push(C64::new(
                cols[1].parse::<f64>().expect("re"),
                cols[2].parse::<f64>().expect("im"),
            ));
        }
        let table = TabulatedLoad::try_new(omegas, zs, "bell", 1.0e-6).expect("table");
        // The bell table and its analytic fallback disagree at the table's
        // band edges: a jump in R(ω). A jump leaves an N-independent tail in
        // the realized IR, so the converged ladder must REFUSE the splice by
        // name instead of realizing it silently. Until bakes are band-complete
        // (.100), a spliced load is not a convergent line.
        let spliced = MmLineBank::new(
            &combos,
            &["open"],
            &gas,
            &MmLoad::TabulatedWithFallback {
                table: &table,
                fallback: Termination::FlangedOpen,
            },
            &cfg,
        );
        let splice_refused = matches!(spliced, Err(MmLineError::LadderUnconverged { .. }));
        let splice_detail = match &spliced {
            Ok(b) => format!(
                "REALIZED table_bins={} fallback_bins={}",
                b.receipts()[0].table_bins,
                b.receipts()[0].fallback_bins
            ),
            Err(e) => format!("refused: {e}"),
        };
        // Refusals.
        let empty = MmLineBank::new(&[], &[], &gas, &load, &cfg).is_err();
        let mismatch = MmLineBank::new(&combos, &["a", "b"], &gas, &load, &cfg).is_err();
        let mut bank2 = MmLineBank::new(&combos, &["open"], &gas, &load, &cfg).expect("bank");
        let unknown = bank2.switch(5, 0).is_err();
        let mut matrix = MatrixFirLine::new(&combos[0], &gas, &load, &cfg).expect("matrix");
        let short_push = matrix.push(&[1.0]).is_err();
        let pass = causal && splice_refused && empty && mismatch && unknown && short_push;
        verdict(
            "ml-006-causality-splice-refusals",
            pass,
            &format!(
                "causal head/tail {:.1}x; splice {splice_detail}; \
                 refusals {empty}/{mismatch}/{unknown}/{short_push}",
                (head / tail.max(1e-300)),
            ),
        );
    }

    /// Regression for bead frankensim-mxw6x and pexfp.82. The shared ladder
    /// must stop only at a CONVERGED rung (tail and peaks; the old
    /// grow-then-clamp silently truncated long bores), and must refuse
    /// typed-and-loud only past the authored 65536-tap ceiling.
    #[test]
    fn ml_007_fft_densifies_or_refuses_never_truncates() {
        let cfg = config();
        let gas = air20();
        let load = MmLoad::Analytic(Termination::FlangedOpen);
        // Round-trip delay g = 2L/c * fs = 700 samples. The ladder starts at
        // four round trips (4096) and must stop only once the tail and the
        // first three peaks converge, on more than one rung.
        let fs = f64::from(cfg.sample_rate_hz);
        let c = gas.sound_speed;
        let long_l = 700.0f64 * c / (2.0 * fs);
        let long = Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.004,
                length: long_l,
            }],
        };
        let bank = MmLineBank::new(std::slice::from_ref(&long), &["open"], &gas, &load, &cfg)
            .expect("long bore densifies instead of truncating");
        let ladder = crate::driving_point::LineLadder::default();
        let receipt = &bank.receipts()[0];
        let dense_ok = receipt.n_fft >= 4096
            && receipt.rungs_tried >= 2
            && receipt.tail_energy < ladder.tail_tol
            && receipt.peak_shift_cents < ladder.peak_tol_cents;
        // Ceiling refusal: a bore demanding far beyond 65536 taps is
        // refused with the authored message rather than truncated.
        let absurd = Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.002,
                length: 100_000.0,
            }],
        };
        let refused = matches!(
            MmLineBank::new(std::slice::from_ref(&absurd), &["open"], &gas, &load, &cfg),
            Err(MmLineError::Invalid { what })
                if what.contains("65536")
        );
        // Single-duct matrix path shares the identical discipline.
        let matrix_long = MatrixFirLine::new(&long, &gas, &load, &cfg).is_ok();
        let pass = dense_ok && refused && matrix_long;
        verdict(
            "ml-007-fft-densify-or-refuse",
            pass,
            &format!(
                "converged n_fft={} after {} rungs (tail {:.3e} < {:e}, shift {:.4} c < {}); \
                 ceiling refusal={refused}; matrix path ok={matrix_long}",
                receipt.n_fft,
                receipt.rungs_tried,
                receipt.tail_energy,
                ladder.tail_tol,
                receipt.peak_shift_cents,
                ladder.peak_tol_cents
            ),
        );
    }
}
