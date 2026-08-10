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

/// Run the jet-labium fixture.
///
/// # Errors
/// [`AeroacError::InvalidParameter`] on inconsistent geometry (plate
/// reaching into the fringe, slot taller than the domain, non-power-
/// of-two record length, degenerate sizes);
/// [`AeroacError::NonFinite`] on bad reals.
#[allow(clippy::too_many_lines)]
pub fn run_jet_labium(cfg: &JetLabiumConfig) -> Result<JetLabiumRun, AeroacError> {
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
    if !cfg.steps_record.is_power_of_two() || cfg.steps_record < 64 {
        return Err(AeroacError::InvalidParameter {
            what: "steps_record must be a power of two >= 64",
        });
    }
    // --- Build the grid: jet everywhere, plate as walls. ---
    let mut grid = Grid::uniform(cfg.nx, cfg.ny, cfg.tau);
    let profile: Vec<(f64, [f64; 2])> = (0..cfg.ny)
        .map(|y| (1.0, [jet_profile(cfg, y), 0.0]))
        .collect();
    for x in 0..cfg.nx {
        #[allow(clippy::cast_precision_loss)]
        let phase = 2.0 * core::f64::consts::PI * x as f64 / cfg.nx as f64;
        let seed = cfg.seed_amplitude * cfg.u_jet * det::sin(phase);
        for (y, row) in profile.iter().enumerate() {
            let i = grid.idx(x, y);
            // Sinuous-symmetry transverse seed shaped by the jet
            // profile: breaks the mirror symmetry deterministically
            // so the oscillation saturates instead of riding on
            // roundoff.
            let vy = seed * row.1[0] / cfg.u_jet;
            grid.f[i] = fs_lbm::equilibrium(1.0, row.1[0], vy);
        }
    }
    // Nozzle wall: a slit of the slot height in a solid column at
    // the domain start; the fringe (which targets the slot profile)
    // sits at the domain END, so the wrap feeds the nozzle plenum.
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
    // --- Settle, then record. ---
    let mut scratch = Vec::new();
    for _ in 0..cfg.steps_settle {
        grid.step(&mut scratch);
        fringe.apply(&mut grid);
    }
    let mut force_series = Vec::with_capacity(cfg.steps_record);
    let mut mach_max = 0.0f64;
    let x_plate_plane = plate_x0.saturating_sub(6).max(cfg.nozzle_thickness + 1);
    let x_fringe_plane = cfg.nx - cfg.fringe_width - 4;
    let mut flux_plate = 0.0f64;
    let mut flux_fringe = 0.0f64;
    for t in 0..cfg.steps_record {
        let exchange = grid.step_with_wall_momentum(&mut scratch, &plate_mask);
        fringe.apply(&mut grid);
        force_series.push(exchange.wall_impulse);
        // Sampled diagnostics (every 16 steps keeps cost negligible).
        if t % 16 == 0 {
            for y in 0..cfg.ny {
                let m = grid.moments(grid.idx(x_plate_plane, y));
                flux_plate += m.rho * m.u[0];
                let sp = det::sqrt(m.u[0] * m.u[0] + m.u[1] * m.u[1]);
                mach_max = mach_max.max(sp);
                let m2 = grid.moments(grid.idx(x_fringe_plane, y));
                flux_fringe += m2.rho * m2.u[0];
                let sp2 = det::sqrt(m2.u[0] * m2.u[0] + m2.u[1] * m2.u[1]);
                mach_max = mach_max.max(sp2);
            }
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let samples = (cfg.steps_record / 16) as f64;
    let nu = (cfg.tau - 0.5) / 3.0;
    Ok(JetLabiumRun {
        force_series,
        diagnostics: JetLabiumDiagnostics {
            mach_max_lattice: mach_max,
            flux_plate_plane: flux_plate / samples,
            flux_fringe_plane: flux_fringe / samples,
            reynolds: cfg.u_jet * 2.0 * cfg.slot_half / nu,
        },
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
