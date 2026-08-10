//! Jet-labium base flow and dipole source extraction — the
//! production half of bead 9ok02: an fs-lbm slot jet impinging on a
//! sharp-edged splitter plate (the classic edge-tone geometry),
//! surface force recorded through the momentum-exchange machinery,
//! and the transverse-force spectrum radiated via the 2D Curle
//! dipole.
//!
//! Architecture: the domain is PERIODIC in x with a FRINGE layer
//! (per-row-profile [`fs_lbm::sponge::Sponge2`]) re-conditioning the
//! outflow toward the authored slot-jet profile, so the wrap
//! delivers fresh inflow — the spectral-DNS fringe method, chosen
//! because the crate's uniform-face open boundary cannot prescribe a
//! slot, and because the fringe doubles as the MEASURED acoustic
//! absorber (R ~ 1e-4, fs-lbm sponge battery). The labium is a
//! two-cell-thick splitter plate with a sharp leading edge — a
//! DISCLOSED simplification of a wedge (the edge, not the wedge
//! angle, drives the oscillation class).
//!
//! The rig is MULTI-STABLE (executed: at Re 144 the selected
//! attractor — St 0.0366 vs the neighboring 0.0458 state — varies
//! with seed amplitude AND lattice resolution, consistent with the
//! documented hysteresis of edge tones). Fixed-parameter runs
//! therefore cannot demonstrate stage selection; the
//! [`run_adiabatic_ramp`] protocol follows ONE attractor
//! continuously by ramping viscosity (Reynolds number) per-step at
//! fixed jet velocity, measuring a spectral rung after each
//! quasi-static transition, up the stage ladder and back down. The
//! up/down legs at equal Reynolds expose hysteresis directly.
//!
//! Every run carries diagnostics (max Mach, plate-plane vs
//! fringe-plane mass-flux ratio, Reynolds number) and the crate's
//! [`crate::SCOPE_STATEMENT`]: spectra are SHAPE/SCALING authorities,
//! never absolute SPL.

use crate::curle2d::dipole_pressure;
use crate::{AeroacError, SCOPE_STATEMENT};
use fs_lbm::core2::{Cell, Grid};
use fs_lbm::sponge::{Sponge2, SpongeSide};
use fs_math::c64::C64;
use fs_math::det;

/// Configuration of one jet-labium run (lattice units throughout).
#[derive(Debug, Clone)]
pub struct JetLabiumConfig {
    /// Domain size (x is periodic through the fringe).
    pub nx: usize,
    /// Domain height (y periodic; the jet must decay well inside).
    pub ny: usize,
    /// Slot half-height parameter of the smoothed top-hat profile.
    pub slot_half: f64,
    /// Profile edge smoothing width [cells].
    pub slot_smoothing: f64,
    /// Jet peak velocity [lu/step].
    pub u_jet: f64,
    /// Relaxation time (viscosity = (tau - 1/2)/3).
    pub tau: f64,
    /// Distance from the wrap plane (x = 0) to the plate leading
    /// edge.
    pub edge_distance: usize,
    /// Plate length downstream of the edge.
    pub plate_length: usize,
    /// Fringe width at the right side.
    pub fringe_width: usize,
    /// Fringe strength.
    pub fringe_sigma: f64,
    /// Settling steps before recording.
    pub steps_settle: usize,
    /// Recorded steps (power of two, for the caller's FFT).
    pub steps_record: usize,
    /// Transverse seed amplitude RELATIVE to `u_jet` (0 disables).
    /// THE VACUOUS-OSCILLATION TRAP (executed): the rig is mirror-
    /// symmetric, and unseeded runs preserve that symmetry to
    /// roundoff — the force spectrum then shows structured, high-
    /// prominence peaks that are AMPLIFIED MACHINE NOISE (~1e-15
    /// lattice force). Frequency-selection claims survive (the
    /// instability amplifies roundoff at the physically selected
    /// frequency) but ANY amplitude-bearing claim needs this seed to
    /// reach the saturated limit cycle.
    pub seed_amplitude: f64,
    /// Nozzle wall thickness at the jet root (0 = no nozzle). A
    /// nozzle wall with a slit of the slot height provides the
    /// RECEPTIVITY edge the classic edge-tone feedback loop closes
    /// at — without it the rig oscillates at the free jet's own
    /// most-amplified frequency instead of the Brown stage ladder
    /// (executed: St 0.46 vs stage-I 0.036 at h/delta = 10).
    pub nozzle_thickness: usize,
}

/// Per-run diagnostics (the bead's mandated honesty block).
#[derive(Debug, Clone)]
pub struct JetLabiumDiagnostics {
    /// Maximum |u| observed over sampled recording steps [lu/step].
    pub mach_max_lattice: f64,
    /// Mean x mass flux through a plane just upstream of the plate.
    pub flux_plate_plane: f64,
    /// Mean x mass flux through a plane just before the fringe.
    pub flux_fringe_plane: f64,
    /// Jet Reynolds number `u_jet * 2 slot_half / nu`.
    pub reynolds: f64,
}

/// One completed run: the recorded surface-force series and
/// diagnostics. Radiation happens separately
/// ([`dipole_spectrum_line`]) so consumers can window/average as
/// they choose.
#[derive(Debug, Clone)]
pub struct JetLabiumRun {
    /// Per-step force on the plate `(Fx, Fy)` [lattice momentum /
    /// step], length `steps_record`.
    pub force_series: Vec<[f64; 2]>,
    /// Diagnostics.
    pub diagnostics: JetLabiumDiagnostics,
    /// The honest-scope statement (embedded in every output).
    pub scope: &'static str,
}

/// The smoothed top-hat slot profile.
fn jet_profile(cfg: &JetLabiumConfig, y: usize) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let yc = cfg.ny as f64 / 2.0 - 0.5;
    #[allow(clippy::cast_precision_loss)]
    let dy = (y as f64 - yc).abs();
    cfg.u_jet * 0.5 * (1.0 + det::tanh((cfg.slot_half - dy) / cfg.slot_smoothing))
}

/// Shared parameter validation for fixed runs and ramps.
fn validate_config(cfg: &JetLabiumConfig) -> Result<(), AeroacError> {
    for (v, what) in [
        (cfg.slot_half, "slot_half"),
        (cfg.slot_smoothing, "slot_smoothing"),
        (cfg.u_jet, "u_jet"),
        (cfg.seed_amplitude, "seed_amplitude"),
        (cfg.tau, "tau"),
        (cfg.fringe_sigma, "fringe_sigma"),
    ] {
        if !v.is_finite() {
            return Err(AeroacError::NonFinite { what });
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let ny_f = cfg.ny as f64;
    if cfg.nx < 64
        || cfg.ny < 32
        || cfg.u_jet <= 0.0
        || cfg.u_jet * cfg.u_jet >= 0.03
        || cfg.tau <= 0.5
        || cfg.slot_half <= 1.0
        || 2.0 * cfg.slot_half >= 0.5 * ny_f
        || cfg.slot_smoothing <= 0.0
        || !(0.0..0.2).contains(&cfg.seed_amplitude)
    {
        return Err(AeroacError::InvalidParameter {
            what: "jet geometry (domain, slot, speed, tau) out of range",
        });
    }
    if cfg.nozzle_thickness + 2 + cfg.edge_distance + cfg.plate_length + cfg.fringe_width + 8
        > cfg.nx
    {
        return Err(AeroacError::InvalidParameter {
            what: "nozzle + plate + fringe do not fit in the domain",
        });
    }
    Ok(())
}

/// The assembled lattice rig: grid with nozzle and plate walls, the
/// fringe, and the diagnostic planes. Construction order and the
/// step/fringe/record loops are EXACTLY those of the original
/// fixed-parameter fixture so [`run_jet_labium`] stays bitwise
/// identical across the refactor.
struct Rig {
    grid: Grid,
    fringe: Sponge2,
    plate_mask: Vec<bool>,
    x_plate_plane: usize,
    x_fringe_plane: usize,
}

/// One recorded segment: the force series plus sampled diagnostics
/// (fluxes already averaged over the samples).
struct RecordedSegment {
    force_series: Vec<[f64; 2]>,
    mach_max: f64,
    flux_plate: f64,
    flux_fringe: f64,
}

impl Rig {
    fn build(cfg: &JetLabiumConfig) -> Rig {
        // --- Build the grid: jet everywhere, plate as walls. ---
        let mut grid = Grid::uniform(cfg.nx, cfg.ny, cfg.tau);
        // Fringe target: WITH a nozzle the target must be the BINARY
        // slit profile matching the wall opening exactly — a smooth
        // (tanh) target carries momentum outside the slit that slams
        // into the wall shoulder on every wrap, continuously feeding
        // the slit-lip shear layer (executed: a St_delta ~ 2-3 lip
        // mode outgrew the edge tone and blocked grid convergence).
        // Without a nozzle the smooth profile stands (no wall to
        // mismatch).
        #[allow(clippy::cast_precision_loss)]
        let yc_prof = cfg.ny as f64 / 2.0 - 0.5;
        let profile: Vec<(f64, [f64; 2])> = (0..cfg.ny)
            .map(|y| {
                let ux = if cfg.nozzle_thickness > 0 {
                    #[allow(clippy::cast_precision_loss)]
                    let open = (y as f64 - yc_prof).abs() < cfg.slot_half;
                    if open { cfg.u_jet } else { 0.0 }
                } else {
                    jet_profile(cfg, y)
                };
                (1.0, [ux, 0.0])
            })
            .collect();
        for x in 0..cfg.nx {
            #[allow(clippy::cast_precision_loss)]
            let phase = 2.0 * core::f64::consts::PI * x as f64 / cfg.nx as f64;
            let seed = cfg.seed_amplitude * cfg.u_jet * det::sin(phase);
            for (y, row) in profile.iter().enumerate() {
                let i = grid.idx(x, y);
                // Sinuous-symmetry transverse seed shaped by the jet
                // profile: breaks the mirror symmetry
                // deterministically so the oscillation saturates
                // instead of riding on roundoff.
                let vy = seed * row.1[0] / cfg.u_jet;
                grid.f[i] = fs_lbm::equilibrium(1.0, row.1[0], vy);
            }
        }
        // Nozzle wall: a slit of the slot height in a solid column at
        // the domain start; the fringe (which targets the slot
        // profile) sits at the domain END, so the wrap feeds the
        // nozzle plenum.
        #[allow(clippy::cast_precision_loss)]
        let yc = cfg.ny as f64 / 2.0 - 0.5;
        for x in 0..cfg.nozzle_thickness {
            for y in 0..cfg.ny {
                #[allow(clippy::cast_precision_loss)]
                let open = (y as f64 - yc).abs() < cfg.slot_half;
                if !open {
                    let i = grid.idx(x, y);
                    grid.flags[i] = Cell::Wall;
                }
            }
        }
        let plate_x0 = cfg.nozzle_thickness + cfg.edge_distance;
        let y_plate_lo = cfg.ny / 2 - 1;
        let mut plate_mask = vec![false; cfg.nx * cfg.ny];
        for x in plate_x0..plate_x0 + cfg.plate_length {
            for y in [y_plate_lo, y_plate_lo + 1] {
                let i = grid.idx(x, y);
                grid.flags[i] = Cell::Wall;
                plate_mask[i] = true;
            }
        }
        let fringe = Sponge2::with_profile(
            SpongeSide::RightX,
            cfg.fringe_width,
            cfg.fringe_sigma,
            &profile,
        );
        let x_plate_plane = plate_x0.saturating_sub(6).max(cfg.nozzle_thickness + 1);
        let x_fringe_plane = cfg.nx - cfg.fringe_width - 4;
        Rig {
            grid,
            fringe,
            plate_mask,
            x_plate_plane,
            x_fringe_plane,
        }
    }

    fn settle(&mut self, steps: usize) {
        let mut scratch = Vec::new();
        for _ in 0..steps {
            self.grid.step(&mut scratch);
            self.fringe.apply(&mut self.grid);
        }
    }

    /// Set the uniform relaxation time on every cell (the ramp's
    /// per-step actuator; viscosity = (tau - 1/2)/3).
    fn set_tau(&mut self, tau: f64) {
        for t in &mut self.grid.tau {
            *t = tau;
        }
    }

    fn record(&mut self, steps: usize) -> RecordedSegment {
        let mut scratch = Vec::new();
        let mut force_series = Vec::with_capacity(steps);
        let mut mach_max = 0.0f64;
        let mut flux_plate = 0.0f64;
        let mut flux_fringe = 0.0f64;
        for t in 0..steps {
            let exchange = self
                .grid
                .step_with_wall_momentum(&mut scratch, &self.plate_mask);
            self.fringe.apply(&mut self.grid);
            force_series.push(exchange.wall_impulse);
            // Sampled diagnostics (every 16 steps keeps cost
            // negligible).
            if t % 16 == 0 {
                for y in 0..self.grid.ny {
                    let m = self.grid.moments(self.grid.idx(self.x_plate_plane, y));
                    flux_plate += m.rho * m.u[0];
                    let sp = det::sqrt(m.u[0] * m.u[0] + m.u[1] * m.u[1]);
                    mach_max = mach_max.max(sp);
                    let m2 = self.grid.moments(self.grid.idx(self.x_fringe_plane, y));
                    flux_fringe += m2.rho * m2.u[0];
                    let sp2 = det::sqrt(m2.u[0] * m2.u[0] + m2.u[1] * m2.u[1]);
                    mach_max = mach_max.max(sp2);
                }
            }
        }
        #[allow(clippy::cast_precision_loss)]
        let samples = (steps / 16) as f64;
        RecordedSegment {
            force_series,
            mach_max,
            flux_plate: flux_plate / samples,
            flux_fringe: flux_fringe / samples,
        }
    }
}

/// Run the jet-labium fixture.
///
/// # Errors
/// [`AeroacError::InvalidParameter`] on inconsistent geometry (plate
/// reaching into the fringe, slot taller than the domain, non-power-
/// of-two record length, degenerate sizes);
/// [`AeroacError::NonFinite`] on bad reals.
pub fn run_jet_labium(cfg: &JetLabiumConfig) -> Result<JetLabiumRun, AeroacError> {
    validate_config(cfg)?;
    if !cfg.steps_record.is_power_of_two() || cfg.steps_record < 64 {
        return Err(AeroacError::InvalidParameter {
            what: "steps_record must be a power of two >= 64",
        });
    }
    let mut rig = Rig::build(cfg);
    // --- Settle, then record. ---
    rig.settle(cfg.steps_settle);
    let seg = rig.record(cfg.steps_record);
    let nu = (cfg.tau - 0.5) / 3.0;
    Ok(JetLabiumRun {
        force_series: seg.force_series,
        diagnostics: JetLabiumDiagnostics {
            mach_max_lattice: seg.mach_max,
            flux_plate_plane: seg.flux_plate,
            flux_fringe_plane: seg.flux_fringe,
            reynolds: cfg.u_jet * 2.0 * cfg.slot_half / nu,
        },
        scope: SCOPE_STATEMENT,
    })
}

/// One measured spectral peak of a transverse-force record.
#[derive(Debug, Clone)]
pub struct ForcePeak {
    /// Strouhal number `f * 2 slot_half / u_jet` at the peak bin.
    pub strouhal: f64,
    /// Peak bin index within the length-n record (frequency =
    /// `bin / n` cycles per step); the Strouhal quantization is
    /// `+-1/(2 bin)` relative, so REPORT the bin with any claim.
    pub bin: usize,
    /// Peak power over the median admitted-bin power (the
    /// oscillation-presence discriminant).
    pub prominence: f64,
    /// Unwindowed RMS of the mean-removed transverse force — the
    /// amplitude floor that catches the vacuous-oscillation trap
    /// (roundoff-riding spectra have structured peaks at ~1e-15
    /// force).
    pub force_rms: f64,
}

/// Hann-windowed periodogram peak of the transverse (`[1]`) force
/// component: mean removal, Hann window, power over the positive
/// half-spectrum, peak search skipping the first `skip_bins` bins
/// (fringe-transient / residual-DC leakage guard), prominence
/// against the median admitted power.
///
/// # Errors
/// [`AeroacError::InvalidParameter`] unless the record length is a
/// power of two >= 64 and `1 <= skip_bins < n/4`;
/// [`AeroacError::NonFinite`]/[`AeroacError::InvalidParameter`] on
/// bad `slot_half`/`u_jet`.
pub fn transverse_force_peak(
    force_series: &[[f64; 2]],
    slot_half: f64,
    u_jet: f64,
    skip_bins: usize,
) -> Result<ForcePeak, AeroacError> {
    let n = force_series.len();
    if !n.is_power_of_two() || n < 64 {
        return Err(AeroacError::InvalidParameter {
            what: "force record length must be a power of two >= 64",
        });
    }
    if skip_bins < 1 || skip_bins >= n / 4 {
        return Err(AeroacError::InvalidParameter {
            what: "skip_bins must be in 1..n/4",
        });
    }
    for (v, what) in [(slot_half, "slot_half"), (u_jet, "u_jet")] {
        if !v.is_finite() {
            return Err(AeroacError::NonFinite { what });
        }
    }
    if slot_half <= 0.0 || u_jet <= 0.0 {
        return Err(AeroacError::InvalidParameter {
            what: "slot_half and u_jet must be positive",
        });
    }
    #[allow(clippy::cast_precision_loss)]
    let n_f = n as f64;
    let mean = force_series.iter().map(|f| f[1]).sum::<f64>() / n_f;
    let rms = det::sqrt(
        force_series
            .iter()
            .map(|f| (f[1] - mean) * (f[1] - mean))
            .sum::<f64>()
            / n_f,
    );
    let fft = fs_fft::Fft::new(n);
    let mut buf: Vec<fs_fft::C64> = force_series
        .iter()
        .enumerate()
        .map(|(i, f)| {
            #[allow(clippy::cast_precision_loss)]
            let w = 0.5 - 0.5 * det::cos((2.0 * core::f64::consts::PI * i as f64) / (n_f - 1.0));
            fs_fft::C64::new((f[1] - mean) * w, 0.0)
        })
        .collect();
    let mut scratch = vec![fs_fft::C64::new(0.0, 0.0); n];
    fft.forward(&mut buf, &mut scratch);
    let power: Vec<f64> = buf[..n / 2].iter().map(|c| c.norm_sq()).collect();
    let (peak_bin, peak_pow) = power
        .iter()
        .enumerate()
        .skip(skip_bins)
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(i, p)| (i, *p))
        .ok_or(AeroacError::InvalidParameter {
            what: "empty admitted spectrum",
        })?;
    let mut sorted = power[skip_bins..].to_vec();
    sorted.sort_by(f64::total_cmp);
    let prominence = peak_pow / sorted[sorted.len() / 2].max(1e-300);
    #[allow(clippy::cast_precision_loss)]
    let strouhal = (peak_bin as f64 / n_f) * 2.0 * slot_half / u_jet;
    Ok(ForcePeak {
        strouhal,
        bin: peak_bin,
        prominence,
        force_rms: rms,
    })
}

/// Which leg of the ramp a rung was measured on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RampDirection {
    /// Reynolds increasing from the base value toward
    /// [`RampConfig::reynolds_end`].
    Up,
    /// Reynolds decreasing back toward the base value.
    Down,
}

/// Configuration of the adiabatic hysteresis-following ramp.
#[derive(Debug, Clone)]
pub struct RampConfig {
    /// The rig at the STARTING Reynolds number; `base.steps_settle`
    /// is the initial lock-in settle and `base.steps_record` is
    /// unused (per-rung records use [`RampConfig::steps_rung_record`]).
    pub base: JetLabiumConfig,
    /// Reynolds number at the top of the ramp. The sweep actuator is
    /// the relaxation time (viscosity) at FIXED `u_jet`, so Mach,
    /// the fringe profile, the nozzle slit, and the seed are all
    /// unchanged along the ramp.
    pub reynolds_end: f64,
    /// Number of measurement rungs per leg (>= 2), inclusive of both
    /// endpoints; rung Reynolds values are uniformly spaced.
    pub rungs: usize,
    /// Steps over which tau is ramped LINEARLY PER STEP between
    /// adjacent rungs. Adiabaticity is the caller's obligation:
    /// choose >= a few oscillation periods (stage I at the staging
    /// geometry is ~2000 steps/period) or the state will be kicked,
    /// not followed.
    pub steps_ramp: usize,
    /// Re-lock settle steps after each transition, before recording.
    pub steps_rung_settle: usize,
    /// Per-rung force record length (power of two >= 64).
    pub steps_rung_record: usize,
    /// Peak-search guard band for [`transverse_force_peak`].
    pub skip_bins: usize,
}

/// One measured rung of the ramp.
#[derive(Debug, Clone)]
pub struct RampRung {
    /// Leg.
    pub direction: RampDirection,
    /// Target Reynolds number of this rung (exact grid value).
    pub reynolds: f64,
    /// Relaxation time realized at this rung.
    pub tau: f64,
    /// Spectral peak of the rung's force record.
    pub peak: ForcePeak,
    /// Max sampled |u| during the rung record [lu/step].
    pub mach_max_lattice: f64,
    /// |plate-plane − fringe-plane| / |plate-plane| mass-flux ratio.
    pub flux_imbalance: f64,
}

/// The completed ramp: `2 * rungs - 1` rung measurements (up leg
/// inclusive of both endpoints, down leg re-measuring every rung
/// except the shared top one).
#[derive(Debug, Clone)]
pub struct RampReport {
    /// Rungs in execution order (up leg then down leg).
    pub rungs: Vec<RampRung>,
    /// The honest-scope statement (embedded in every output).
    pub scope: &'static str,
}

/// Ramp tau linearly per step from `from` to `to` while stepping the
/// flow, then pin the exact endpoint.
fn ramp_tau(rig: &mut Rig, from: f64, to: f64, steps: usize) {
    let mut scratch = Vec::new();
    for s in 1..=steps {
        #[allow(clippy::cast_precision_loss)]
        let tau = from + (to - from) * s as f64 / steps as f64;
        rig.set_tau(tau);
        rig.grid.step(&mut scratch);
        rig.fringe.apply(&mut rig.grid);
    }
    rig.set_tau(to);
}

/// Run the adiabatic hysteresis-following ramp protocol.
///
/// Motivation (executed findings on bead 9ok02): the edge-tone rig
/// is multi-stable — fixed-parameter runs land on seed- and
/// resolution-dependent attractors, so neither stage selection nor
/// attractor convergence can be demonstrated by independent runs.
/// This protocol locks the base state once, then quasi-statically
/// sweeps Reynolds (viscosity at fixed jet speed) so the flow
/// FOLLOWS one attractor branch; discrete Strouhal jumps between
/// rungs are stage transitions, and disagreement between the up and
/// down legs at equal Reynolds is the hysteresis measurement.
///
/// The returned rungs are measurements, not claims: stage
/// identification against the published ladder is the caller's
/// comparison, and [`crate::SCOPE_STATEMENT`] applies to every
/// spectral quantity.
///
/// # Errors
/// [`AeroacError::InvalidParameter`] on the base-config violations
/// of [`run_jet_labium`], fewer than two rungs, zero ramp/settle
/// steps, a non-power-of-two record, a degenerate (< 1% relative)
/// Reynolds span, or a top-of-ramp relaxation time at or below the
/// 0.5005 stability floor; [`AeroacError::NonFinite`] on bad reals.
pub fn run_adiabatic_ramp(cfg: &RampConfig) -> Result<RampReport, AeroacError> {
    validate_config(&cfg.base)?;
    if cfg.rungs < 2 {
        return Err(AeroacError::InvalidParameter {
            what: "ramp needs at least two rungs",
        });
    }
    if cfg.steps_ramp == 0 || cfg.steps_rung_settle == 0 {
        return Err(AeroacError::InvalidParameter {
            what: "steps_ramp and steps_rung_settle must be positive",
        });
    }
    if !cfg.steps_rung_record.is_power_of_two() || cfg.steps_rung_record < 64 {
        return Err(AeroacError::InvalidParameter {
            what: "steps_rung_record must be a power of two >= 64",
        });
    }
    if !cfg.reynolds_end.is_finite() {
        return Err(AeroacError::NonFinite {
            what: "reynolds_end",
        });
    }
    let delta = 2.0 * cfg.base.slot_half;
    let nu0 = (cfg.base.tau - 0.5) / 3.0;
    let re0 = cfg.base.u_jet * delta / nu0;
    if cfg.reynolds_end <= 0.0 || (cfg.reynolds_end - re0).abs() / re0 < 0.01 {
        return Err(AeroacError::InvalidParameter {
            what: "reynolds_end must differ from the base Reynolds by >= 1%",
        });
    }
    let tau_of = |re: f64| 0.5 + 3.0 * cfg.base.u_jet * delta / re;
    // The highest Reynolds on the ramp has the lowest tau; refuse
    // below the authored stability floor rather than running an
    // under-relaxed lattice.
    if tau_of(cfg.reynolds_end.max(re0)) <= 0.5005 {
        return Err(AeroacError::InvalidParameter {
            what: "top-of-ramp tau at or below the 0.5005 stability floor",
        });
    }
    let rung_re: Vec<f64> = (0..cfg.rungs)
        .map(|i| {
            #[allow(clippy::cast_precision_loss)]
            let frac = i as f64 / (cfg.rungs - 1) as f64;
            re0 + (cfg.reynolds_end - re0) * frac
        })
        .collect();
    let mut rig = Rig::build(&cfg.base);
    rig.settle(cfg.base.steps_settle);
    let mut tau_now = cfg.base.tau;
    let mut rungs = Vec::with_capacity(2 * cfg.rungs - 1);
    let measure = |rig: &mut Rig,
                   reynolds: f64,
                   tau: f64,
                   direction: RampDirection|
     -> Result<RampRung, AeroacError> {
        let seg = rig.record(cfg.steps_rung_record);
        let peak = transverse_force_peak(
            &seg.force_series,
            cfg.base.slot_half,
            cfg.base.u_jet,
            cfg.skip_bins,
        )?;
        let flux_imbalance = (seg.flux_plate - seg.flux_fringe).abs() / seg.flux_plate.abs();
        Ok(RampRung {
            direction,
            reynolds,
            tau,
            peak,
            mach_max_lattice: seg.mach_max,
            flux_imbalance,
        })
    };
    // Up leg (rung 0 is measured directly out of the base settle).
    for (i, &re) in rung_re.iter().enumerate() {
        if i > 0 {
            let tau_next = tau_of(re);
            ramp_tau(&mut rig, tau_now, tau_next, cfg.steps_ramp);
            tau_now = tau_next;
            rig.settle(cfg.steps_rung_settle);
        }
        rungs.push(measure(&mut rig, re, tau_now, RampDirection::Up)?);
    }
    // Down leg (the top rung is shared, not re-measured).
    for &re in rung_re[..cfg.rungs - 1].iter().rev() {
        let tau_next = tau_of(re);
        ramp_tau(&mut rig, tau_now, tau_next, cfg.steps_ramp);
        tau_now = tau_next;
        rig.settle(cfg.steps_rung_settle);
        rungs.push(measure(&mut rig, re, tau_now, RampDirection::Down)?);
    }
    Ok(RampReport {
        rungs,
        scope: SCOPE_STATEMENT,
    })
}

/// Radiate one spectral line of the recorded force through the 2D
/// Curle dipole: given the complex force amplitude `f_hat` at
/// wavenumber `k` (the caller FFTs `force_series` and converts
/// frequency to `k = omega / c`), the pressure at `observer` relative
/// to the plate edge at `source`.
///
/// # Errors
/// As [`crate::curle2d::dipole_pressure`].
pub fn dipole_spectrum_line(
    f_hat: [C64; 2],
    k: f64,
    observer: [f64; 2],
    source: [f64; 2],
) -> Result<C64, AeroacError> {
    Ok(dipole_pressure(f_hat, k, observer, source)?.pressure)
}
