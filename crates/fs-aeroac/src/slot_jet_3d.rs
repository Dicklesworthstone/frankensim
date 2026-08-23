//! The 3-D slot-jet rig on the D3Q19 central-moment operator — the
//! staged-but-never-executed follow-up named by the recorded 2-D
//! tonal-lock refusal ([`crate::regime`], bead frankensim-l011o) and
//! owned by bead frankensim-music-v8-root-3ez8g.10.1.
//!
//! Architecture: the 2-D fringe recipe ported to the third dimension.
//! The domain is PERIODIC in all axes; a right-side FRINGE layer
//! (per-y-row equilibrium targets, uniform in z) re-conditions the
//! outflow toward the authored slot-jet profile, so the x wrap
//! delivers fresh inflow — the fringe doubles as the measured
//! acoustic absorber (its 3-D reflection coefficient is re-measured,
//! not assumed: the 2-D fs-lbm sponge battery number does not
//! transfer automatically). The nozzle is a one-column binary slit
//! wall (the executed 2-D lesson: the smooth-shoulder mismatch fed a
//! slit-lip mode and was fixed once by going BINARY — kept binary
//! here); the splitter plate is a two-cell-thick wall strip with a
//! sharp leading edge spanning the full z span.
//!
//! DISCIPLINE LAWS on this rig family (bead 10.1):
//! - explicit symmetry-breaking seed: a sinuous y seed shaped by the
//!   profile AND a z-antisymmetric factor (a z-uniform seed would
//!   preserve the spanwise mirror symmetry and let amplitude claims
//!   ride amplified roundoff);
//! - every amplitude-bearing claim carries the measured force RMS
//!   against an absolute floor ([`FORCE_RMS_AMPLITUDE_FLOOR`]);
//! - fs-lbm's non-finite-density panic is contained at settle/record
//!   boundaries and surfaced as the typed [`AeroacError::NonFinite`]
//!   refusal;
//! - FFT bin quantization is disclosed on every receipt
//!   ([`SlotJet3dDiagnostics::strouhal_bin_width`]).
//!
//! HONEST OUTCOMES: EITHER a demonstrated broadband rung (flatness
//! above [`crate::regime::TONAL_FLATNESS_CEILING`] with its
//! Re/resolution boundary) OR a quantified refusal extending the 2-D
//! one. Both are wins; neither mints an experimental claim. The
//! [`crate::SCOPE_STATEMENT`] applies to every spectral quantity.

use crate::AeroacError;
use crate::SCOPE_STATEMENT;
use crate::jetlab::transverse_force_spectrum;
use crate::regime::{SpectrumClass, classify_spectrum, measure_spectral_flatness};
use fs_lbm::d3q19::{BoundaryGrid3, BoundarySpec3, CollisionModel3, E3, Q3, equilibrium3};
use fs_math::det;

/// Absolute transverse-force RMS floor [lattice momentum/step] below
/// which a rung's AMPLITUDE carries no claim (the vacuous-oscillation
/// trap: mirror-symmetric rigs oscillate in amplified roundoff at
/// ~1e-15 lattice force while prominence cannot detect it).
pub const FORCE_RMS_AMPLITUDE_FLOOR: f64 = 1.0e-12;

/// Configuration of one 3-D slot-jet run (lattice units throughout).
#[derive(Debug, Clone)]
pub struct SlotJet3dConfig {
    /// Domain size in x (periodic through the fringe; multiple of 4).
    pub nx: usize,
    /// Domain height in y, the cross-stream direction (multiple of 4).
    pub ny: usize,
    /// Domain span in z, the spanwise direction (multiple of 4). The
    /// box-sensitivity octave doubles this extent.
    pub nz: usize,
    /// Binary slit half-width in y [cells]; the opening is every cell
    /// whose center satisfies `|y + 1/2 - ny/2| < slot_half`.
    pub slot_half: f64,
    /// Jet peak velocity [lu/step].
    pub u_jet: f64,
    /// Collision law. The lane is pinned to the D3Q19
    /// [`CollisionModel3::CentralMoment`] operator (the unlocked 3-D
    /// path); other models are refused.
    pub collision: CollisionModel3,
    /// Nozzle wall thickness at the jet root (>= 1; the binary slit
    /// receptivity edge the feedback loop closes at).
    pub nozzle_thickness: usize,
    /// Distance from the nozzle exit to the plate leading edge.
    pub edge_distance: usize,
    /// Plate length downstream of the edge (two cells thick in y,
    /// full z span).
    pub plate_length: usize,
    /// Fringe width at the right side.
    pub fringe_width: usize,
    /// Fringe strength.
    pub fringe_sigma: f64,
    /// Sinuous seed amplitude RELATIVE to `u_jet` (0 disables; a
    /// disabled seed voids amplitude claims by the roundoff law).
    pub seed_amplitude: f64,
    /// Settling steps before recording.
    pub steps_settle: usize,
    /// Recorded steps (power of two, for the shared FFT pipeline).
    pub steps_record: usize,
}

/// Per-run diagnostics (the mandated honesty block).
#[derive(Debug, Clone)]
pub struct SlotJet3dDiagnostics {
    /// Maximum sampled plane speed over recording steps [lu/step].
    pub mach_max_lattice: f64,
    /// Mean x mass flux through a plane just upstream of the plate.
    pub flux_plate_plane: f64,
    /// Mean x mass flux through a plane just before the fringe.
    pub flux_fringe_plane: f64,
    /// Jet Reynolds number `u_jet * 2 slot_half / nu` with the
    /// central-moment viscosity.
    pub reynolds: f64,
    /// Recorded length (FFT input size).
    pub record_len: usize,
    /// FFT bin quantization disclosure: Strouhal spacing
    /// `(1/n) * 2 slot_half / u_jet`. Every Strouhal claim carries
    /// this; a peak is a BIN, not a continuous frequency.
    pub strouhal_bin_width: f64,
}

/// One completed run: the recorded plate-force series and
/// diagnostics. Classification happens separately ([`classify_rung`])
/// through the shared 2-D measurement pipeline.
#[derive(Debug, Clone)]
pub struct SlotJet3dRun {
    /// Per-step impulse on the splitter plate `(Fx, Fy)` [lattice
    /// momentum/step], length `steps_record`. The plate spans the
    /// full periodic z, so the net z impulse integrates to roundoff
    /// and is not recorded.
    pub force_series: Vec<[f64; 2]>,
    /// Diagnostics.
    pub diagnostics: SlotJet3dDiagnostics,
    /// The honest-scope statement (embedded in every output).
    pub scope: &'static str,
}

/// One classified rung of the Re ladder.
#[derive(Debug, Clone)]
pub struct SlotJet3dRung {
    /// Jet Reynolds number.
    pub reynolds: f64,
    /// Central-moment second-order rate (the Re actuator at fixed
    /// `u_jet`).
    pub second_order_rate: f64,
    /// Central-moment higher-order rate.
    pub higher_order_rate: f64,
    /// Measured spectral flatness (geometric/arithmetic mean).
    pub flatness: f64,
    /// Classification against [`crate::regime::TONAL_FLATNESS_CEILING`].
    pub tonal: bool,
    /// Peak Strouhal (a BIN — see `strouhal_bin_width`).
    pub strouhal: f64,
    /// Peak bin index.
    pub peak_bin: usize,
    /// Peak prominence over the median admitted bin.
    pub prominence: f64,
    /// Measured transverse-force RMS.
    pub force_rms: f64,
    /// Whether the force RMS clears [`FORCE_RMS_AMPLITUDE_FLOOR`] —
    /// amplitude claims require this.
    pub amplitude_qualified: bool,
    /// Strouhal bin width disclosure.
    pub strouhal_bin_width: f64,
    /// Maximum sampled plane speed.
    pub mach_max_lattice: f64,
    /// Relative plate-vs-fringe plane mass-flux imbalance.
    pub flux_imbalance: f64,
}

/// The right-side fringe layer over a [`BoundaryGrid3`]: per-y-row
/// equilibrium targets (uniform in z), quadratic sigma ramp growing
/// into the layer — the 3-D port of [`fs_lbm::sponge::Sponge2`]
/// with-profile mode. Exposed so the reflection battery can drive it
/// directly.
#[derive(Debug, Clone)]
pub struct Fringe3 {
    width: usize,
    sigma_max: f64,
    /// One equilibrium target per y-row.
    targets: Vec<[f64; Q3]>,
}

impl Fringe3 {
    /// Construct a fringe whose target varies per y-row:
    /// `profile[y] = (rho, [ux, uy, uz])`.
    ///
    /// # Panics
    /// On zero width, `sigma_max` outside (0, 1], an empty profile,
    /// non-finite/non-positive targets, or a target speed outside
    /// the low-Mach envelope (crate boundary convention).
    #[must_use]
    pub fn with_profile(width: usize, sigma_max: f64, profile: &[(f64, [f64; 3])]) -> Self {
        assert!(width > 0, "fringe width must be positive");
        assert!(
            sigma_max.is_finite() && sigma_max > 0.0 && sigma_max <= 1.0,
            "fringe sigma_max must lie in (0, 1]"
        );
        assert!(!profile.is_empty(), "fringe profile must be non-empty");
        let targets = profile
            .iter()
            .map(|&(rho, u)| {
                assert!(
                    rho.is_finite() && rho > 0.0,
                    "fringe target density must be positive and finite"
                );
                assert!(
                    u.iter().all(|c| c.is_finite()),
                    "fringe target must be finite"
                );
                assert!(
                    u[0] * u[0] + u[1] * u[1] + u[2] * u[2] < 0.03,
                    "fringe target velocity exceeds the low-Mach boundary envelope"
                );
                equilibrium3(rho, u)
            })
            .collect();
        Fringe3 {
            width,
            sigma_max,
            targets,
        }
    }

    /// Blend strength at depth `d` into the layer (0 = inner edge).
    #[must_use]
    fn sigma(&self, d: usize) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let t = (d as f64 + 1.0) / self.width as f64;
        self.sigma_max * t * t
    }

    /// Apply one blending pass (call once per LBM step). Only fluid
    /// cells are touched; solid cells keep their bounce-back state.
    ///
    /// # Panics
    /// If the layer is wider than the grid or a profile length does
    /// not match `ny`.
    pub fn apply(&self, grid: &mut BoundaryGrid3) {
        let [nx, ny, nz] = grid.dimensions();
        assert!(self.width <= nx, "fringe layer wider than the grid");
        assert!(
            self.targets.len() == ny,
            "fringe profile length must equal the grid height"
        );
        for d in 0..self.width {
            let x = nx - self.width + d;
            let s = self.sigma(d);
            for y in 0..ny {
                let target = &self.targets[y];
                for z in 0..nz {
                    if grid.is_solid(x, y, z) {
                        continue;
                    }
                    let mut f = grid.populations(x, y, z);
                    for (fq, tq) in f.iter_mut().zip(target) {
                        *fq += s * (tq - *fq);
                    }
                    grid.set_populations(x, y, z, &f);
                }
            }
        }
    }
}

/// The binary slit opening predicate at cell index `y` for height
/// `ny`: cell centers with `|y + 1/2 - ny/2| < slot_half`.
#[must_use]
fn slit_open(y: usize, ny: usize, slot_half: f64) -> bool {
    #[allow(clippy::cast_precision_loss)]
    let dy = ((y as f64 + 0.5) - ny as f64 / 2.0).abs();
    dy < slot_half
}

/// Shared parameter validation.
fn validate_config(cfg: &SlotJet3dConfig) -> Result<(), AeroacError> {
    for (v, what) in [
        (cfg.slot_half, "slot_half"),
        (cfg.u_jet, "u_jet"),
        (cfg.seed_amplitude, "seed_amplitude"),
        (cfg.fringe_sigma, "fringe_sigma"),
    ] {
        if !v.is_finite() {
            return Err(AeroacError::NonFinite { what });
        }
    }
    if !matches!(cfg.collision, CollisionModel3::CentralMoment { .. }) {
        return Err(AeroacError::InvalidParameter {
            what: "the 3-D slot-jet lane is pinned to the D3Q19 central-moment operator",
        });
    }
    cfg.collision
        .validate()
        .map_err(|_| AeroacError::InvalidParameter {
            what: "3-D collision relaxation is outside the D3Q19 physical window",
        })?;
    #[allow(clippy::cast_precision_loss)]
    let ny_f = cfg.ny as f64;
    if cfg.nx == 0
        || cfg.ny == 0
        || cfg.nz == 0
        || !cfg.nx.is_multiple_of(4)
        || !cfg.ny.is_multiple_of(4)
        || !cfg.nz.is_multiple_of(4)
        || cfg.nx < 64
        || cfg.ny < 32
        || cfg.nz < 8
        || cfg.u_jet <= 0.0
        || cfg.u_jet * cfg.u_jet >= 0.03
        || cfg.slot_half <= 1.0
        || 2.0 * cfg.slot_half >= 0.5 * ny_f
        || cfg.nozzle_thickness < 1
        || cfg.edge_distance < 8
        || cfg.plate_length < 2
        || !(0.0..0.2).contains(&cfg.seed_amplitude)
        || cfg.steps_settle == 0
    {
        return Err(AeroacError::InvalidParameter {
            what: "jet geometry (domain, slot, speed, rates, layout) out of range",
        });
    }
    if cfg.nozzle_thickness + cfg.edge_distance + cfg.plate_length + cfg.fringe_width + 8 > cfg.nx {
        return Err(AeroacError::InvalidParameter {
            what: "nozzle + plate + fringe do not fit in the domain",
        });
    }
    if !cfg.steps_record.is_power_of_two() || cfg.steps_record < 64 {
        return Err(AeroacError::InvalidParameter {
            what: "steps_record must be a power of two >= 64",
        });
    }
    Ok(())
}

/// The assembled 3-D rig: periodic grid with nozzle/plate solids, the
/// fringe, and the diagnostic planes.
struct Rig3 {
    grid: BoundaryGrid3,
    fringe: Fringe3,
    /// Plate solid cells (for the momentum-exchange force).
    plate_cells: Vec<(usize, usize, usize)>,
    x_plate_plane: usize,
    x_fringe_plane: usize,
}

impl Rig3 {
    fn build(cfg: &SlotJet3dConfig) -> Rig3 {
        let CollisionModel3::CentralMoment {
            second_order_rate,
            higher_order_rate,
        } = cfg.collision
        else {
            unreachable!("validate_config pins the central-moment operator");
        };
        let mut grid = BoundaryGrid3::with_collision_model(
            cfg.nx,
            cfg.ny,
            cfg.nz,
            CollisionModel3::CentralMoment {
                second_order_rate,
                higher_order_rate,
            },
            [0.0; 3],
            BoundarySpec3::periodic(),
        );
        // Nozzle (binary slit) + splitter plate through the sanctioned
        // geometry entry. The closure receives cell-center samples;
        // only the sign is consumed.
        let plate_x0 = cfg.nozzle_thickness + cfg.edge_distance;
        #[allow(clippy::cast_precision_loss)]
        let plate_yc_lo = (cfg.ny / 2 - 1) as f64 + 0.5;
        #[allow(clippy::cast_precision_loss)]
        let plate_yc_hi = (cfg.ny / 2) as f64 + 0.5;
        let nozzle = cfg.nozzle_thickness;
        let slot = cfg.slot_half;
        let plate_len = cfg.plate_length;
        #[allow(clippy::cast_precision_loss)]
        let ny_center = cfg.ny as f64 / 2.0;
        grid.voxelize_sdf(|[sx, sy, _]| {
            let open = (sy - ny_center).abs() < slot;
            let in_nozzle_solid = sx < nozzle as f64 && !open;
            let in_plate = sx >= plate_x0 as f64
                && sx < (plate_x0 + plate_len) as f64
                && sy >= plate_yc_lo
                && sy <= plate_yc_hi;
            if in_nozzle_solid || in_plate {
                -1.0
            } else {
                1.0
            }
        });
        // Initial field: jet profile everywhere + the sinuous seed
        // with a z-antisymmetric factor (the spanwise mirror-symmetry
        // breaker).
        for z in 0..cfg.nz {
            #[allow(clippy::cast_precision_loss)]
            let z_break =
                1.0 + 0.5 * det::sin(2.0 * core::f64::consts::PI * z as f64 / cfg.nz as f64);
            for y in 0..cfg.ny {
                let ux = if slit_open(y, cfg.ny, cfg.slot_half) {
                    cfg.u_jet
                } else {
                    0.0
                };
                for x in 0..cfg.nx {
                    if grid.is_solid(x, y, z) {
                        continue;
                    }
                    #[allow(clippy::cast_precision_loss)]
                    let phase = 2.0 * core::f64::consts::PI * x as f64 / cfg.nx as f64;
                    let vy = cfg.seed_amplitude * det::sin(phase) * ux * z_break;
                    let f = equilibrium3(1.0, [ux, vy, 0.0]);
                    grid.set_populations(x, y, z, &f);
                }
            }
        }
        // Fringe target: the BINARY slit profile matching the nozzle
        // opening exactly (the executed 2-D lesson).
        let profile: Vec<(f64, [f64; 3])> = (0..cfg.ny)
            .map(|y| {
                let ux = if slit_open(y, cfg.ny, cfg.slot_half) {
                    cfg.u_jet
                } else {
                    0.0
                };
                (1.0, [ux, 0.0, 0.0])
            })
            .collect();
        let fringe = Fringe3::with_profile(cfg.fringe_width, cfg.fringe_sigma, &profile);
        // Plate solid cells (momentum-exchange surface).
        let mut plate_cells = Vec::new();
        for x in plate_x0..plate_x0 + cfg.plate_length {
            for y in [cfg.ny / 2 - 1, cfg.ny / 2] {
                for z in 0..cfg.nz {
                    if grid.is_solid(x, y, z) {
                        plate_cells.push((x, y, z));
                    }
                }
            }
        }
        let x_plate_plane = plate_x0.saturating_sub(6).max(cfg.nozzle_thickness + 1);
        let x_fringe_plane = cfg.nx - cfg.fringe_width - 4;
        Rig3 {
            grid,
            fringe,
            plate_cells,
            x_plate_plane,
            x_fringe_plane,
        }
    }

    /// Settle without recording, with the fs-lbm density-assert panic
    /// contained and surfaced as the typed refusal.
    fn settle(&mut self, steps: usize) -> Result<(), AeroacError> {
        std::panic::catch_unwind(core::panic::AssertUnwindSafe(|| {
            for _ in 0..steps {
                self.grid.step();
                self.fringe.apply(&mut self.grid);
            }
        }))
        .map_err(|_| AeroacError::NonFinite {
            what: "3-D lattice destabilized during settle (fs-lbm density assert)",
        })
    }

    /// Typed stability guard: finiteness of every fluid population
    /// and positivity of every cell density.
    fn check_stable(&self) -> Result<(), AeroacError> {
        let [nx, ny, nz] = self.grid.dimensions();
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    if self.grid.is_solid(x, y, z) {
                        continue;
                    }
                    let f = self.grid.populations(x, y, z);
                    let mut rho = 0.0;
                    for q in f {
                        if !q.is_finite() {
                            return Err(AeroacError::NonFinite {
                                what: "3-D lattice destabilized (non-finite distribution)",
                            });
                        }
                        rho += q;
                    }
                    if rho <= 0.0 {
                        return Err(AeroacError::InvalidParameter {
                            what: "3-D lattice destabilized (non-positive cell density)",
                        });
                    }
                }
            }
        }
        Ok(())
    }

    /// Momentum-exchange impulse on the plate: every population at a
    /// fluid cell directed into a plate solid cell transfers
    /// `2 f_i e_i` through the halfway bounce-back.
    fn plate_impulse(&self) -> [f64; 2] {
        let [nx, ny, nz] = self.grid.dimensions();
        let mut impulse = [0.0f64; 2];
        for &(x, y, z) in &self.plate_cells {
            for (i, e) in E3.iter().enumerate().skip(1) {
                // Fluid cell whose e_i neighbor is this plate cell.
                // rem_euclid keeps the periodic wrap exact for
                // negative lattice directions.
                let qx = (x as isize - e.0 as isize).rem_euclid(nx as isize) as usize;
                let qy = (y as isize - e.1 as isize).rem_euclid(ny as isize) as usize;
                let qz = (z as isize - e.2 as isize).rem_euclid(nz as isize) as usize;
                if self.grid.is_solid(qx, qy, qz) {
                    continue;
                }
                let fi = self.grid.populations(qx, qy, qz)[i];
                impulse[0] += 2.0 * fi * f64::from(e.0);
                impulse[1] += 2.0 * fi * f64::from(e.1);
            }
        }
        impulse
    }

    fn record(&mut self, steps: usize) -> Result<Vec<[f64; 2]>, AeroacError> {
        self.check_stable()?;
        let this = &mut *self;
        std::panic::catch_unwind(core::panic::AssertUnwindSafe(move || {
            let mut force_series = Vec::with_capacity(steps);
            for _ in 0..steps {
                this.grid.step();
                this.fringe.apply(&mut this.grid);
                let impulse = this.plate_impulse();
                if !impulse[0].is_finite() || !impulse[1].is_finite() {
                    return Err(AeroacError::NonFinite {
                        what: "3-D lattice destabilized (non-finite plate force)",
                    });
                }
                force_series.push(impulse);
            }
            Ok(force_series)
        }))
        .unwrap_or(Err(AeroacError::NonFinite {
            what: "3-D lattice destabilized during record (fs-lbm density assert)",
        }))
    }

    /// Sampled plane diagnostics.
    fn plane_diagnostics(&self) -> (f64, f64, f64) {
        let [_, ny, nz] = self.grid.dimensions();
        let mut mach_max = 0.0f64;
        let mut flux_plate = 0.0;
        let mut flux_fringe = 0.0;
        for z in 0..nz {
            for y in 0..ny {
                if self.grid.is_solid(self.x_plate_plane, y, z) {
                    continue;
                }
                let f = self.grid.populations(self.x_plate_plane, y, z);
                let rho: f64 = f.iter().sum();
                let ux: f64 = f.iter().zip(E3).map(|(q, e)| q * f64::from(e.0)).sum();
                flux_plate += rho * ux;
                mach_max = mach_max.max(det::sqrt(ux * ux));
                if self.grid.is_solid(self.x_fringe_plane, y, z) {
                    continue;
                }
                let g = self.grid.populations(self.x_fringe_plane, y, z);
                let rho2: f64 = g.iter().sum();
                let ux2: f64 = g.iter().zip(E3).map(|(q, e)| q * f64::from(e.0)).sum();
                flux_fringe += rho2 * ux2;
                mach_max = mach_max.max(det::sqrt(ux2 * ux2));
            }
        }
        #[allow(clippy::cast_precision_loss)]
        let cells = (ny * nz) as f64;
        (mach_max, flux_plate / cells, flux_fringe / cells)
    }
}

/// Run the 3-D slot-jet fixture.
///
/// # Errors
/// [`AeroacError::InvalidParameter`] on inconsistent geometry or a
/// non-central-moment collision; [`AeroacError::NonFinite`] on bad
/// reals and (with named messages) on lattice destabilization.
pub fn run_slot_jet_3d(cfg: &SlotJet3dConfig) -> Result<SlotJet3dRun, AeroacError> {
    validate_config(cfg)?;
    let mut rig = Rig3::build(cfg);
    rig.settle(cfg.steps_settle)?;
    let force_series = rig.record(cfg.steps_record)?;
    let (mach_max, flux_plate, flux_fringe) = rig.plane_diagnostics();
    let CollisionModel3::CentralMoment {
        second_order_rate, ..
    } = cfg.collision
    else {
        unreachable!("validate_config pins the central-moment operator");
    };
    let nu = (1.0 / second_order_rate - 0.5) / 3.0;
    #[allow(clippy::cast_precision_loss)]
    let reynolds = cfg.u_jet * 2.0 * cfg.slot_half / nu;
    #[allow(clippy::cast_precision_loss)]
    let strouhal_bin_width = (1.0 / cfg.steps_record as f64) * 2.0 * cfg.slot_half / cfg.u_jet;
    Ok(SlotJet3dRun {
        force_series,
        diagnostics: SlotJet3dDiagnostics {
            mach_max_lattice: mach_max,
            flux_plate_plane: flux_plate,
            flux_fringe_plane: flux_fringe,
            reynolds,
            record_len: cfg.steps_record,
            strouhal_bin_width,
        },
        scope: SCOPE_STATEMENT,
    })
}

/// Classify one run through the shared measurement pipeline
/// ([`crate::jetlab::transverse_force_spectrum`] +
/// [`crate::regime::measure_spectral_flatness`]).
///
/// # Errors
/// Pipeline refusals (record length, degenerate spectrum).
pub fn classify_rung(
    run: &SlotJet3dRun,
    cfg: &SlotJet3dConfig,
) -> Result<SlotJet3dRung, AeroacError> {
    let CollisionModel3::CentralMoment {
        second_order_rate,
        higher_order_rate,
    } = cfg.collision
    else {
        return Err(AeroacError::InvalidParameter {
            what: "the 3-D slot-jet lane is pinned to the D3Q19 central-moment operator",
        });
    };
    // skip_bins must satisfy 1 <= skip < n/4 (pipeline refusal);
    // n/8 keeps the fringe-transient guard inside the window.
    let (power, peak) = transverse_force_spectrum(
        &run.force_series,
        cfg.slot_half,
        cfg.u_jet,
        run.diagnostics.record_len / 8,
    )?;
    let flatness = measure_spectral_flatness(&power)?;
    let tonal = matches!(classify_spectrum(&power)?, SpectrumClass::Tonal { .. });
    let force_rms = peak.force_rms;
    let imbalance = if run.diagnostics.flux_plate_plane.abs() > f64::EPSILON {
        (run.diagnostics.flux_plate_plane - run.diagnostics.flux_fringe_plane).abs()
            / run.diagnostics.flux_plate_plane.abs()
    } else {
        f64::INFINITY
    };
    Ok(SlotJet3dRung {
        reynolds: run.diagnostics.reynolds,
        second_order_rate,
        higher_order_rate,
        flatness,
        tonal,
        strouhal: peak.strouhal,
        peak_bin: peak.bin,
        prominence: peak.prominence,
        force_rms,
        amplitude_qualified: force_rms > FORCE_RMS_AMPLITUDE_FLOOR,
        strouhal_bin_width: run.diagnostics.strouhal_bin_width,
        mach_max_lattice: run.diagnostics.mach_max_lattice,
        flux_imbalance: imbalance,
    })
}

impl SlotJet3dRung {
    /// One deterministic JSONL receipt line (schema
    /// `fs-aeroac.slot-jet-3d.rung/v1`). Field order is pinned; the
    /// receipt carries the FFT bin disclosure and the amplitude
    /// qualification on every line.
    #[must_use]
    pub fn to_jsonl(&self) -> String {
        format!(
            "{{\"schema\":\"fs-aeroac.slot-jet-3d.rung/v1\",\"reynolds\":{},\"second_order_rate\":{},\
\"higher_order_rate\":{},\"flatness\":{},\"tonal\":{},\"strouhal\":{},\"peak_bin\":{},\
\"prominence\":{},\"force_rms\":{},\"amplitude_qualified\":{},\"strouhal_bin_width\":{},\
\"mach_max_lattice\":{},\"flux_imbalance\":{}}}",
            self.reynolds,
            self.second_order_rate,
            self.higher_order_rate,
            self.flatness,
            self.tonal,
            self.strouhal,
            self.peak_bin,
            self.prominence,
            self.force_rms,
            self.amplitude_qualified,
            self.strouhal_bin_width,
            self.mach_max_lattice,
            self.flux_imbalance
        )
    }
}

/// Outcome of one chunked invocation.
#[derive(Debug, Clone)]
pub enum SweepProgress {
    /// The full settle+record finished; the run is complete.
    Complete(Box<SlotJet3dRun>),
    /// The step budget expired first; resume with the same
    /// checkpoint directory and configuration to continue.
    Partial {
        /// Total lattice steps completed across all chunks.
        steps_done: u64,
    },
}

/// Canonical configuration fingerprint (FNV-1a over the pinned
/// fields) binding a checkpoint to exactly one geometry/law.
#[must_use]
pub fn config_fingerprint(cfg: &SlotJet3dConfig) -> u64 {
    let CollisionModel3::CentralMoment {
        second_order_rate,
        higher_order_rate,
    } = cfg.collision
    else {
        return 0;
    };
    let tag = format!(
        "v1|{},{},{},{},{},{},{},{},{},{},{},{},{}",
        cfg.nx,
        cfg.ny,
        cfg.nz,
        cfg.slot_half,
        cfg.u_jet,
        second_order_rate,
        higher_order_rate,
        cfg.nozzle_thickness,
        cfg.edge_distance,
        cfg.plate_length,
        cfg.fringe_width,
        cfg.fringe_sigma,
        cfg.seed_amplitude
    );
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for byte in tag.as_bytes() {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

struct ChunkState {
    settle_done: u64,
    record_done: u64,
    force_series: Vec<[f64; 2]>,
}

const CHECKPOINT_MAGIC: &[u8; 8] = b"FSJSCKPT";

fn write_checkpoint(
    dir: &std::path::Path,
    cfg: &SlotJet3dConfig,
    rig: &Rig3,
    state: &ChunkState,
) -> Result<(), AeroacError> {
    use std::io::Write;
    std::fs::create_dir_all(dir).map_err(|_| AeroacError::InvalidParameter {
        what: "cannot create checkpoint directory",
    })?;
    let path = dir.join("state.bin");
    let tmp = dir.join("state.bin.tmp");
    let file = std::fs::File::create(&tmp).map_err(|_| AeroacError::InvalidParameter {
        what: "cannot create checkpoint file",
    })?;
    let mut w = std::io::BufWriter::new(file);
    let refuse = |_e: std::io::Error| AeroacError::InvalidParameter {
        what: "checkpoint write",
    };
    w.write_all(CHECKPOINT_MAGIC).map_err(refuse)?;
    w.write_all(&1u32.to_le_bytes()).map_err(refuse)?;
    w.write_all(&config_fingerprint(cfg).to_le_bytes())
        .map_err(refuse)?;
    w.write_all(&state.settle_done.to_le_bytes())
        .map_err(refuse)?;
    w.write_all(&state.record_done.to_le_bytes())
        .map_err(refuse)?;
    #[allow(clippy::cast_possible_truncation)]
    let record_len = state.force_series.len() as u64;
    w.write_all(&record_len.to_le_bytes()).map_err(refuse)?;
    for pair in &state.force_series {
        w.write_all(&pair[0].to_le_bytes()).map_err(refuse)?;
        w.write_all(&pair[1].to_le_bytes()).map_err(refuse)?;
    }
    // Fluid-cell populations in canonical x-fastest order; solid
    // cells are skipped (reconstructed by the rig on load).
    let [nx, ny, nz] = rig.grid.dimensions();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if rig.grid.is_solid(x, y, z) {
                    continue;
                }
                for q in rig.grid.populations(x, y, z) {
                    w.write_all(&q.to_le_bytes()).map_err(refuse)?;
                }
            }
        }
    }
    w.flush().map_err(refuse)?;
    std::fs::rename(&tmp, &path).map_err(|_| AeroacError::InvalidParameter {
        what: "cannot finalize checkpoint (atomic rename)",
    })?;
    Ok(())
}

/// Run the sweep under an explicit per-invocation step budget,
/// checkpointing atomically so repeated invocations chain into one
/// bit-identical execution (the recorded checkpoint-between-rungs
/// recipe under bounded worker job walls).
///
/// # Errors
/// Config validation, checkpoint I/O refusals, fingerprint
/// mismatch, and the usual destabilization refusals.
pub fn run_slot_jet_3d_chunked(
    cfg: &SlotJet3dConfig,
    checkpoint_dir: &std::path::Path,
    step_budget: usize,
) -> Result<SweepProgress, AeroacError> {
    validate_config(cfg)?;
    if step_budget == 0 {
        return Err(AeroacError::InvalidParameter {
            what: "step budget must be positive",
        });
    }
    let ckpt_path = checkpoint_dir.join("state.bin");
    let (mut rig, mut state) = if ckpt_path.exists() {
        load_checkpoint(cfg, &ckpt_path)?
    } else {
        (
            Rig3::build(cfg),
            ChunkState {
                settle_done: 0,
                record_done: 0,
                force_series: Vec::new(),
            },
        )
    };
    step_chunk_loop(cfg, &mut rig, &mut state, step_budget)?;
    let total_settle = u64::try_from(cfg.steps_settle).unwrap_or(u64::MAX);
    let total_record = u64::try_from(cfg.steps_record).unwrap_or(u64::MAX);
    if state.settle_done >= total_settle && state.record_done >= total_record {
        let _ = std::fs::remove_file(&ckpt_path);
        let (mach_max, flux_plate, flux_fringe) = rig.plane_diagnostics();
        return Ok(SweepProgress::Complete(Box::new(finish_run(
            cfg,
            state.force_series,
            mach_max,
            flux_plate,
            flux_fringe,
        ))));
    }
    write_checkpoint(checkpoint_dir, cfg, &rig, &state)?;
    Ok(SweepProgress::Partial {
        steps_done: state.settle_done + state.record_done,
    })
}

/// Jet Reynolds number for a validated central-moment config.
#[must_use]
pub fn reynolds_of(cfg: &SlotJet3dConfig) -> f64 {
    let CollisionModel3::CentralMoment {
        second_order_rate, ..
    } = cfg.collision
    else {
        return f64::NAN;
    };
    let nu = (1.0 / second_order_rate - 0.5) / 3.0;
    // No casts here; the allow below is unnecessary — plain float math.
    cfg.u_jet * 2.0 * cfg.slot_half / nu
}

/// FFT Strouhal bin width for the record length (the quantization
/// disclosure carried on every receipt).
#[must_use]
pub fn strouhal_bin_width_of(cfg: &SlotJet3dConfig) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let inv = 1.0 / cfg.steps_record as f64;
    inv * 2.0 * cfg.slot_half / cfg.u_jet
}

fn finish_run(
    cfg: &SlotJet3dConfig,
    force_series: Vec<[f64; 2]>,
    mach_max: f64,
    flux_plate: f64,
    flux_fringe: f64,
) -> SlotJet3dRun {
    SlotJet3dRun {
        force_series,
        diagnostics: SlotJet3dDiagnostics {
            mach_max_lattice: mach_max,
            flux_plate_plane: flux_plate,
            flux_fringe_plane: flux_fringe,
            reynolds: reynolds_of(cfg),
            record_len: cfg.steps_record,
            strouhal_bin_width: strouhal_bin_width_of(cfg),
        },
        scope: SCOPE_STATEMENT,
    }
}

fn load_checkpoint(
    cfg: &SlotJet3dConfig,
    ckpt_path: &std::path::Path,
) -> Result<(Rig3, ChunkState), AeroacError> {
    let bytes = std::fs::read(ckpt_path).map_err(|_| AeroacError::InvalidParameter {
        what: "cannot read checkpoint",
    })?;
    if bytes.len() < 44 || bytes[0..8] != *CHECKPOINT_MAGIC {
        return Err(AeroacError::InvalidParameter {
            what: "checkpoint magic/version mismatch",
        });
    }
    let stored_fp = u64::from_le_bytes(bytes[12..20].try_into().expect("fixed slice"));
    if stored_fp != config_fingerprint(cfg) {
        return Err(AeroacError::InvalidParameter {
            what: "checkpoint belongs to a different configuration",
        });
    }
    let settle_done = u64::from_le_bytes(bytes[20..28].try_into().expect("fixed slice"));
    let record_done = u64::from_le_bytes(bytes[28..36].try_into().expect("fixed slice"));
    let record_len = usize::try_from(u64::from_le_bytes(
        bytes[36..44].try_into().expect("fixed slice"),
    ))
    .unwrap_or(usize::MAX);
    if bytes.len() < 44 + 16 * record_len {
        return Err(AeroacError::InvalidParameter {
            what: "checkpoint truncated (force series)",
        });
    }
    let mut force_series = Vec::with_capacity(record_len);
    let mut off = 44usize;
    for _ in 0..record_len {
        let fx = f64::from_le_bytes(bytes[off..off + 8].try_into().expect("fixed slice"));
        let fy = f64::from_le_bytes(bytes[off + 8..off + 16].try_into().expect("fixed slice"));
        force_series.push([fx, fy]);
        off += 16;
    }
    let mut rig = Rig3::build(cfg);
    for z in 0..cfg.nz {
        for y in 0..cfg.ny {
            for x in 0..cfg.nx {
                if rig.grid.is_solid(x, y, z) {
                    continue;
                }
                let mut pops = [0.0f64; Q3];
                for q in &mut pops {
                    if off + 8 > bytes.len() {
                        return Err(AeroacError::InvalidParameter {
                            what: "checkpoint truncated (lattice)",
                        });
                    }
                    *q = f64::from_le_bytes(bytes[off..off + 8].try_into().expect("fixed slice"));
                    off += 8;
                }
                if !pops.iter().all(|v| v.is_finite()) {
                    return Err(AeroacError::NonFinite {
                        what: "checkpoint carries non-finite populations",
                    });
                }
                rig.grid.set_populations(x, y, z, &pops);
            }
        }
    }
    Ok((
        rig,
        ChunkState {
            settle_done,
            record_done,
            force_series,
        },
    ))
}

fn step_chunk_loop(
    cfg: &SlotJet3dConfig,
    rig: &mut Rig3,
    state: &mut ChunkState,
    step_budget: usize,
) -> Result<(), AeroacError> {
    let total_settle = u64::try_from(cfg.steps_settle).unwrap_or(u64::MAX);
    let total_record = u64::try_from(cfg.steps_record).unwrap_or(u64::MAX);
    let mut budget = step_budget as u64;
    while budget > 0 {
        if state.settle_done < total_settle {
            let stepped = std::panic::catch_unwind(core::panic::AssertUnwindSafe(|| {
                rig.grid.step();
                rig.fringe.apply(&mut rig.grid);
            }));
            if stepped.is_err() {
                return Err(AeroacError::NonFinite {
                    what: "3-D lattice destabilized during chunked settle",
                });
            }
            state.settle_done += 1;
        } else if state.record_done < total_record {
            let stepped = std::panic::catch_unwind(core::panic::AssertUnwindSafe(|| {
                rig.grid.step();
                rig.fringe.apply(&mut rig.grid);
                rig.plate_impulse()
            }));
            let impulse = match stepped {
                Ok(impulse) => impulse,
                Err(_) => {
                    return Err(AeroacError::NonFinite {
                        what: "3-D lattice destabilized during chunked record",
                    });
                }
            };
            if !impulse[0].is_finite() || !impulse[1].is_finite() {
                return Err(AeroacError::NonFinite {
                    what: "3-D lattice destabilized (non-finite plate force)",
                });
            }
            state.force_series.push(impulse);
            state.record_done += 1;
        } else {
            break;
        }
        budget -= 1;
    }
    Ok(())
}
