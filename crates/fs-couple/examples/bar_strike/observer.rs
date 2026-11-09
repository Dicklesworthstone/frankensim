//! The engine's broadband modal observer on a struck bar, as an optional card row.
//!
//! The narrow-band receiver holds each mode's transfer at the mode's own frequency. While the mallet is on
//! the bar a high mode's displacement follows the force and not its own frequency, so that hold gives a
//! pressure that is not a wave: the file does not arrive at a second receiver as a delayed copy. This module
//! samples every mode's Burton-Miller receiver transfer across a band and hands the samples to
//! `BroadbandModalObserver`, whose continuous fit integrates the actual within-step modal trajectory.
//!
//! The observer is one-way. It steps its own copy of the mechanics under the forces the bar already felt, so
//! the bar, the contact and the mallet cannot change. With the row absent none of this runs.
//!
//! The grid follows the bar_s1 test: log-spaced to 2 kHz, then a 250 Hz step, one factorization per frequency
//! with every mode as a right-hand side. Its top is the BEM's six-panels-per-wavelength floor on this surface
//! or the grader's 16 kHz, whichever is lower. Each mode's own frequency is also a training sample.
//!
//! With `observer_rigid_rows on` the rigid coordinates the strike drives become right-hand sides too. While
//! the mallet is on, the whole bar accelerates as F/m, and that is a dipole the elastic rows do not contain.
//! A rigid coordinate's velocity stays constant after contact and a true transfer vanishes at zero frequency,
//! but the fit is free below the grid's lowest frequency, so its gain there is reported and not assumed.
use crate::Structure;
use crate::bar::{TAU, hz};
use crate::card::Card;
use crate::radiation::{Transfers, surface_and_fields_at};
use fs_bem::helmholtz::{
    Formulation, HelmholtzError, MIN_PANELS_PER_WAVELENGTH, Medium, exterior_pressure_at_points,
    solve_radiation_batch,
};
use fs_bem::panel3d::SpherePanels;
use fs_couple::modal_acoustic_time::ModalAcousticTimeModel;
use fs_couple::modal_acoustic_time::broadband::{
    BroadbandModalControls, BroadbandModalObserver, ModalTransferSamples,
};
use fs_math::{c64::C64, det};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

const GRID_LOW_HZ: f64 = 50.0;
const GRID_SPLIT_HZ: f64 = 2_000.0;
const GRID_PER_OCTAVE: f64 = 6.0;
const GRID_STEP_HZ: f64 = 250.0;
/// The grader's top band edge (apps/instrument-lab/tools/instrgrade/descriptors.py, ERB_FMAX_HZ).
const RECEIVER_BAND_HZ: f64 = 16_000.0;
const FIT_ITERATIONS: usize = 32;
const DELAY_KERNEL_HALF_LEN: usize = 32;
/// The orders the cost report tries. The bar_s1 test's 16, 32 and 48, and lower ones: how many poles a row
/// will take depends on the row and the band, and a row that is too smooth refuses with `SigmaCollapse`.
const FIT_ORDERS: [usize; 5] = [8, 12, 16, 32, 48];
/// Steps timed in the cost report.
const TIMED_STEPS: usize = 100_000;

/// A batch solve at one's own frequency and the narrow-band solve must be the same number.
const WIRING_TOLERANCE: f64 = 1e-9;

pub fn controls(fit_order: usize, budget: f64) -> BroadbandModalControls {
    BroadbandModalControls {
        fit_order,
        fit_iterations: FIT_ITERATIONS,
        delay_kernel_half_len: DELAY_KERNEL_HALF_LEN,
        max_held_out_rel_error: budget,
    }
}

/// Training grid with each mode's own frequency in it, and a held-out grid from the plain grid alone.
pub fn grid(top_hz: f64, own_hz: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut plain = Vec::new();
    let mut f = GRID_LOW_HZ;
    while f < GRID_SPLIT_HZ {
        plain.push(f);
        f *= det::pow(2.0, 1.0 / GRID_PER_OCTAVE);
    }
    let mut f = GRID_SPLIT_HZ;
    while f <= top_hz {
        plain.push(f);
        f += GRID_STEP_HZ;
    }
    let held = plain
        .windows(2)
        .step_by(4)
        .map(|w| det::sqrt(w[0] * w[1]))
        .collect();
    let mut training = plain;
    training.extend_from_slice(own_hz);
    training.sort_by(f64::total_cmp);
    training.dedup();
    (training, held)
}

/// The earliest possible arrival: the receiver's distance to the bar's bounding box.
pub fn earliest_arrival_m(card: &Card, nodes: &[[f64; 3]]) -> f64 {
    let mut sum = 0.0;
    for axis in 0..3 {
        let low = nodes.iter().map(|n| n[axis]).fold(f64::INFINITY, f64::min);
        let high = nodes
            .iter()
            .map(|n| n[axis])
            .fold(f64::NEG_INFINITY, f64::max);
        let d = card.receiver_m[axis] - card.receiver_m[axis].clamp(low, high);
        sum += d * d;
    }
    det::sqrt(sum)
}

struct Solved {
    transfers: Vec<C64>,
    panels_per_wavelength: f64,
}

/// Receiver pressure of every mode at one frequency, one factorization. `None` below the BEM's floor.
fn at_frequency(
    card: &Card,
    medium: Medium,
    surface: &SpherePanels,
    shapes: &[Vec<f64>],
    f_hz: f64,
) -> Result<Option<Solved>, String> {
    let fields: Vec<Vec<C64>> = shapes
        .iter()
        .map(|s| s.iter().map(|&v| C64::from_re(v)).collect())
        .collect();
    let refs: Vec<&[C64]> = fields.iter().map(Vec::as_slice).collect();
    let k = TAU * f_hz / medium.sound_speed;
    match solve_radiation_batch(surface, k, medium, &refs, Formulation::BurtonMiller) {
        Ok(solved) => {
            let panels_per_wavelength = solved[0].panels_per_wavelength;
            let mut transfers = Vec::with_capacity(solved.len());
            for s in &solved {
                let p = exterior_pressure_at_points(surface, s, medium, &[card.receiver_m])
                    .map_err(|e| format!("receiver pressure at {f_hz:.1} Hz: {e:?}"))?;
                transfers.push(p[0]);
            }
            Ok(Some(Solved {
                transfers,
                panels_per_wavelength,
            }))
        }
        Err(HelmholtzError::TooCoarse { .. }) => Ok(None),
        Err(e) => Err(format!("band BEM solve at {f_hz:.1} Hz: {e:?}")),
    }
}

/// One frequency of the band, as solved.
pub struct BandSolve {
    pub f_hz: f64,
    pub training: bool,
    pub transfers: Vec<C64>,
    pub panels_per_wavelength: f64,
    /// Wall time of this frequency on the full surface [s].
    pub wall_s: f64,
    /// max |H - H_coarse| / max |H| over the modes, where the coarser surface is admitted.
    pub ladder_change: Option<f64>,
}

/// The whole bar's translation along the thickness against a compact dipole at the band's lowest frequency,
/// which fixes the sign and the size of the rigid rows' wiring by a route that doesn't use them.
#[derive(Clone, Copy, Debug)]
pub struct RigidCheck {
    pub f_hz: f64,
    /// The bar's own volume, from the coordinate's mass norm and the card's density [m3].
    pub volume_m3: f64,
    /// The volume a compact dipole needs to give this transfer [m3].
    pub effective_volume_m3: f64,
    /// The bar's volume plus a flat strip's added volume, pi W^2 L / 4 [m3].
    pub strip_bound_m3: f64,
    /// Phase of the transfer against the dipole [deg]; zero for a compact dipole.
    pub phase_deg: f64,
}

/// `h` is the receiver transfer of rigid coordinate 2, which is a pure translation along the thickness.
fn rigid_check(card: &Card, medium: Medium, s: &Structure, f_hz: f64, h: C64) -> RigidCheck {
    let phi = s.rigid_fields[2][2];
    let mass = 1.0 / (phi * phi);
    let volume_m3 = mass / card.density_kg_m3;
    let mut d = [0.0; 3];
    for (axis, slot) in d.iter_mut().enumerate() {
        let low = s
            .mesh
            .nodes
            .iter()
            .map(|n| n[axis])
            .fold(f64::INFINITY, f64::min);
        let high = s
            .mesh
            .nodes
            .iter()
            .map(|n| n[axis])
            .fold(f64::NEG_INFINITY, f64::max);
        *slot = card.receiver_m[axis] - 0.5 * (low + high);
    }
    let r = det::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    let w = TAU * f_hz;
    let k = w / medium.sound_speed;
    // exp(-iwt): p = rho V a cos(theta) (1 - ikr) exp(ikr) / (4 pi r^2), with a = -iw v.
    let scale = medium.density * w * phi * (d[2] / r) / (4.0 * core::f64::consts::PI * r * r);
    let dipole =
        C64::new(0.0, -scale) * C64::new(1.0, -k * r) * C64::new(det::cos(k * r), det::sin(k * r));
    let ratio = h * dipole.conj();
    let norm = dipole.abs() * dipole.abs();
    RigidCheck {
        f_hz,
        volume_m3,
        effective_volume_m3: ratio.re / norm,
        strip_bound_m3: volume_m3
            + core::f64::consts::PI * card.width_m * card.width_m * card.length_m / 4.0,
        phase_deg: det::atan2(ratio.im, ratio.re) * 180.0 / core::f64::consts::PI,
    }
}

pub struct Band {
    pub samples: ModalTransferSamples,
    pub solves: Vec<BandSolve>,
    pub top_hz: f64,
    pub panels: usize,
    pub coarse_panels: Option<usize>,
    /// Wall time of the whole band [s], and the threads that shared it.
    pub wall_s: f64,
    pub threads: usize,
    /// Worst |H_batch - H_narrow| / |H_narrow| at the modes' own frequencies.
    pub wiring_error: f64,
    /// The rigid coordinates heard, by index: empty unless the card turns them on. A coordinate the strike
    /// has no force on (projection exactly zero) stays silent.
    pub driven_rigid: Vec<usize>,
    pub rigid_projection: Vec<f64>,
    pub rigid_check: Option<RigidCheck>,
}

/// Solve the band on the card's surface, and on a surface twice as coarse along the length and the width
/// where the grid divides, which gives the producer's estimate. The thickness count stays the same on that
/// surface, so the estimate says nothing about refining through the thickness.
#[allow(clippy::too_many_lines)] // one coherent band solve: grid, threaded surface solves, wiring check
pub fn band(
    card: &Card,
    medium: Medium,
    surface: &SpherePanels,
    shapes: &[Vec<f64>],
    s: &Structure,
    narrow: &Transfers,
) -> Result<Band, String> {
    let started = Instant::now();
    let own_hz: Vec<f64> = s.modes.iter().map(hz).collect();
    let highest = own_hz.iter().copied().fold(0.0, f64::max);
    let top_hz = (highest * narrow.min_panels_per_wavelength / MIN_PANELS_PER_WAVELENGTH)
        .min(RECEIVER_BAND_HZ);
    if !(top_hz > GRID_SPLIT_HZ && top_hz > highest) {
        return Err(format!(
            "the BEM floor at {top_hz:.1} Hz leaves no band above the top mode at {highest:.1} Hz"
        ));
    }
    let (training, held) = grid(top_hz, &own_hz);
    let driven: Vec<usize> = if card.observer.is_some_and(|o| o.rigid_rows) {
        (0..s.rigid.len()).filter(|&r| s.rigid[r] != 0.0).collect()
    } else {
        Vec::new()
    };
    let rigid_fields: Vec<&[f64]> = driven
        .iter()
        .map(|&r| s.rigid_fields[r].as_slice())
        .collect();
    let mut shapes: Vec<Vec<f64>> = shapes.to_vec();
    if !rigid_fields.is_empty() {
        let (_, extra) = surface_and_fields_at(card, &s.mesh, &rigid_fields, [card.bem_stride; 3])?;
        shapes.extend(extra);
    }
    let shapes = &shapes;
    let coarse = {
        let strides = [2 * card.bem_stride, 2 * card.bem_stride, card.bem_stride];
        if (0..3).all(|a| s.mesh.counts[a].is_multiple_of(strides[a])) {
            let all: Vec<&[f64]> = s
                .modes
                .iter()
                .map(|m| m.phi.as_slice())
                .chain(rigid_fields.iter().copied())
                .collect();
            Some(surface_and_fields_at(card, &s.mesh, &all, strides)?)
        } else {
            None
        }
    };
    let all: Vec<(f64, bool)> = training
        .iter()
        .map(|&f| (f, true))
        .chain(held.iter().map(|&f| (f, false)))
        .collect();
    let next = AtomicUsize::new(0);
    let results: Vec<Mutex<Option<Result<BandSolve, String>>>> =
        all.iter().map(|_| Mutex::new(None)).collect();
    let threads = std::thread::available_parallelism()
        .map_or(4, std::num::NonZero::get)
        .min(16);
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= all.len() {
                        break;
                    }
                    let (f_hz, is_training) = all[i];
                    let one = Instant::now();
                    let fine = at_frequency(card, medium, surface, shapes, f_hz);
                    let wall_s = one.elapsed().as_secs_f64();
                    let out = fine.and_then(|solved| {
                        let solved = solved.ok_or_else(|| {
                            format!("the surface is below the BEM floor at {f_hz:.1} Hz")
                        })?;
                        let ladder_change = match &coarse {
                            None => None,
                            Some((c_surface, c_shapes)) => {
                                at_frequency(card, medium, c_surface, c_shapes, f_hz)?.map(|c| {
                                    let scale = solved
                                        .transfers
                                        .iter()
                                        .map(|p| p.abs())
                                        .fold(0.0, f64::max);
                                    solved
                                        .transfers
                                        .iter()
                                        .zip(&c.transfers)
                                        .map(|(a, b)| (*a - *b).abs())
                                        .fold(0.0, f64::max)
                                        / scale
                                })
                            }
                        };
                        Ok(BandSolve {
                            f_hz,
                            training: is_training,
                            transfers: solved.transfers,
                            panels_per_wavelength: solved.panels_per_wavelength,
                            wall_s,
                            ladder_change,
                        })
                    });
                    if let Ok(mut slot) = results[i].lock() {
                        *slot = Some(out);
                    }
                }
            });
        }
    });
    let mut solves = Vec::with_capacity(all.len());
    for slot in results {
        solves.push(
            slot.into_inner()
                .map_err(|_| "poisoned band result".to_string())?
                .ok_or_else(|| "a frequency was never solved".to_string())??,
        );
    }
    let mut wiring_error = 0.0_f64;
    for (mode, f) in own_hz.iter().enumerate() {
        let own = solves
            .iter()
            // `grid` copies each own frequency into the training grid, so its bits match exactly.
            .find(|b| b.training && b.f_hz.to_bits() == f.to_bits())
            .ok_or_else(|| "a mode's own frequency left the training grid".to_string())?;
        let error =
            (own.transfers[mode] - narrow.pressure[mode]).abs() / narrow.pressure[mode].abs();
        wiring_error = wiring_error.max(error);
    }
    if !(wiring_error <= WIRING_TOLERANCE) {
        return Err(format!(
            "the band's own-frequency transfers miss the narrow-band ones by {wiring_error:e}"
        ));
    }
    let rigid = s.rigid.len();
    let modes = own_hz.len();
    let row = |p: &[C64]| -> Vec<C64> {
        let mut out: Vec<C64> = p[..modes].to_vec();
        out.resize(modes + rigid, C64::ZERO);
        for (slot, &r) in driven.iter().enumerate() {
            out[modes + r] = p[modes + slot];
        }
        out
    };
    let rigid_check_row = driven.iter().position(|&r| r == 2).and_then(|slot| {
        solves
            .iter()
            .find(|b| b.training)
            .map(|b| rigid_check(card, medium, s, b.f_hz, b.transfers[modes + slot]))
    });
    let mut per_mode = vec![Vec::with_capacity(training.len()); own_hz.len() + rigid];
    for b in solves.iter().filter(|b| b.training) {
        for (column, value) in per_mode.iter_mut().zip(row(&b.transfers)) {
            column.push(value);
        }
    }
    let held_out = solves
        .iter()
        .filter(|b| !b.training)
        .map(|b| (TAU * b.f_hz, row(&b.transfers)))
        .collect();
    let producer = solves
        .iter()
        .filter_map(|b| b.ladder_change)
        .fold(0.0, f64::max);
    Ok(Band {
        samples: ModalTransferSamples {
            source_id: format!(
                "fs-bem:BurtonMiller:bar-strike:{}:{}-panels:receiver-{:?}",
                card.name, narrow.panels, card.receiver_m
            ),
            omega_rad_s: training.iter().map(|&f| TAU * f).collect(),
            per_mode,
            held_out,
            propagation_delay_s: earliest_arrival_m(card, &s.mesh.nodes) / medium.sound_speed,
            producer_relative_error_estimate: producer,
        },
        solves,
        top_hz,
        panels: narrow.panels,
        coarse_panels: coarse.as_ref().map(|(c, _)| c.centroids().len()),
        wall_s: started.elapsed().as_secs_f64(),
        threads,
        wiring_error,
        rigid_projection: s.rigid.clone(),
        rigid_check: rigid_check_row,
        driven_rigid: driven,
    })
}

/// A value a JSON row may lack: the value itself, or `null` (`{:?}` would print `Some(..)` and `None`).
pub fn json_opt<T: std::fmt::Display>(v: Option<T>) -> String {
    v.map_or_else(|| "null".to_string(), |x| x.to_string())
}

/// One `observer-band-solve` JSON row: one frequency of the band, its transfers and what it cost.
pub fn band_solve_row(card: &Card, b: &BandSolve) -> String {
    let h: Vec<[f64; 2]> = b.transfers.iter().map(|p| [p.re, p.im]).collect();
    format!(
        "{{\"suite\":\"bar-strike\",\"case\":\"observer-band-solve\",\"bar\":\"{}\",\"f_hz\":{},\"role\":\"{}\",\"panels_per_wavelength\":{},\"wall_s\":{},\"ladder_change\":{},\"h_re_im\":{h:?}}}",
        card.name,
        b.f_hz,
        if b.training { "training" } else { "held-out" },
        b.panels_per_wavelength,
        b.wall_s,
        json_opt(b.ladder_change)
    )
}

/// One JSON row per frequency and one for the band: what the band cost and how the surface behaved.
pub fn report_band(card: &Card, band: &Band) {
    for b in &band.solves {
        println!("{}", band_solve_row(card, b));
    }
    let thread_s: f64 = band.solves.iter().map(|b| b.wall_s).sum();
    println!(
        "{{\"suite\":\"bar-strike\",\"case\":\"observer-band\",\"bar\":\"{}\",\"panels\":{},\"coarse_panels\":{},\"top_hz\":{},\"training\":{},\"held_out\":{},\"threads\":{},\"band_wall_s\":{},\"fine_thread_s\":{thread_s},\"fine_s_per_frequency\":{},\"delay_s\":{},\"producer_estimate\":{},\"wiring_error\":{:e}}}",
        card.name,
        band.panels,
        json_opt(band.coarse_panels),
        band.top_hz,
        band.samples.omega_rad_s.len(),
        band.samples.held_out.len(),
        band.threads,
        band.wall_s,
        thread_s / band.solves.len() as f64,
        band.samples.propagation_delay_s,
        band.samples.producer_relative_error_estimate,
        band.wiring_error
    );
    println!(
        "{{\"suite\":\"bar-strike\",\"case\":\"observer-rigid-rows\",\"bar\":\"{}\",\"driven\":{:?},\"projection\":{:?},\"check\":\"{:?}\"}}",
        card.name, band.driven_rigid, band.rigid_projection, band.rigid_check
    );
}

/// Every fit order: how long the fit takes, each mode's error on the withheld frequencies and at its own
/// frequency, and the cost of a step. Nothing here renders.
pub fn report_fits(
    card: &Card,
    model: &ModalAcousticTimeModel,
    band: &Band,
    narrow: &[C64],
    own_hz: &[f64],
) {
    let rate = card.rate_hz * card.substeps;
    let steps = f64::from(rate) * card.duration_s;
    for order in FIT_ORDERS {
        let started = Instant::now();
        let built = BroadbandModalObserver::build(model, &band.samples, controls(order, 1e9), rate);
        let fit_wall_s = started.elapsed().as_secs_f64();
        let (mut observer, report) = match built {
            Ok(b) => b,
            Err(e) => {
                println!(
                    "{{\"suite\":\"bar-strike\",\"case\":\"observer-fit\",\"bar\":\"{}\",\"order\":{order},\"fit_wall_s\":{fit_wall_s},\"refused\":\"{}\"}}",
                    card.name,
                    format!("{e:?}").replace('"', "'")
                );
                continue;
            }
        };
        let modes = own_hz.len();
        // Every row the receiver hears: the elastic modes, then any rigid coordinate that has a transfer.
        let rows = band.samples.per_mode.len();
        let mut held_worst = vec![0.0_f64; rows];
        let mut held_worst_hz = vec![0.0_f64; rows];
        let mut held_all: Vec<f64> = Vec::new();
        for (w, row) in &band.samples.held_out {
            for mode in 0..rows {
                if row[mode] == C64::ZERO {
                    continue;
                }
                let error = observer
                    .transfer_at(mode, *w)
                    .map_or(f64::INFINITY, |h| (h - row[mode]).abs() / row[mode].abs());
                held_all.push(error);
                if error > held_worst[mode] {
                    held_worst[mode] = error;
                    held_worst_hz[mode] = *w / TAU;
                }
            }
        }
        held_all.sort_by(f64::total_cmp);
        let held_median = held_all
            .get(held_all.len() / 2)
            .copied()
            .unwrap_or(f64::NAN);
        let own_db: Vec<f64> = (0..modes)
            .map(|mode| {
                observer
                    .transfer_at(mode, TAU * own_hz[mode])
                    .map_or(f64::NAN, |h| 20.0 * (h.abs() / narrow[mode].abs()).log10())
            })
            .collect();
        let own_error: Vec<f64> = (0..modes)
            .map(|mode| {
                observer
                    .transfer_at(mode, TAU * own_hz[mode])
                    .map_or(f64::INFINITY, |h| {
                        (h - narrow[mode]).abs() / narrow[mode].abs()
                    })
            })
            .collect();
        // A rigid coordinate keeps its velocity after contact, so the fit's gain near zero frequency is the file's
        // offset per unit of that velocity. [row, gain near zero, fit at the grid's lowest, BEM there].
        let low = band.samples.omega_rad_s[0];
        let rigid_low: Vec<[f64; 4]> = band
            .driven_rigid
            .iter()
            .map(|&r| {
                let row = modes + r;
                [
                    r as f64,
                    observer.transfer_at(row, 1e-3).map_or(f64::NAN, C64::abs),
                    observer.transfer_at(row, low).map_or(f64::NAN, C64::abs),
                    band.samples.per_mode[row][0].abs(),
                ]
            })
            .collect();
        let forces = vec![0.0; model.modes().len()];
        let timed = Instant::now();
        let mut stepped = Ok(0.0);
        for _ in 0..TIMED_STEPS {
            stepped = observer.step(&forces);
            if stepped.is_err() {
                break;
            }
        }
        let per_step_s = timed.elapsed().as_secs_f64() / TIMED_STEPS as f64;
        println!(
            "{{\"suite\":\"bar-strike\",\"case\":\"observer-fit\",\"bar\":\"{}\",\"order\":{order},\"fit_wall_s\":{fit_wall_s},\"held_out_max\":{},\"held_out_median\":{held_median},\"held_out_worst_per_mode\":{held_worst:?},\"held_out_worst_hz\":{held_worst_hz:?},\"own_frequency_error\":{own_error:?},\"own_frequency_db\":{own_db:?},\"rigid_low_gain\":{rigid_low:?},\"delay_samples\":{},\"step_ok\":{},\"per_step_s\":{per_step_s},\"render_forecast_s\":{}}}",
            card.name,
            report.held_out_max_rel_error,
            report.delay_samples,
            stepped.is_ok(),
            per_step_s * steps
        );
    }
}
