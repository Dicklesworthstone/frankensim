//! Free-control branch/set semantics (bead wf-root-guzez.5.14.2,
//! E4.6b-ii). Plan §5.1.3 + Round-2 free-control revision: with the
//! pilot released (disconnected), the canard DOF's behavior is a
//! SET-VALUED object, never a single trim point:
//!
//!  - the EQUILIBRIUM SET: every zero crossing of the aero hinge moment
//!    over the travel range, each with a branch identity and a
//!    stability classification (moment slope sign),
//!  - the STICTION SET: the interval(s) where |M_aero| ≤ M_stiction —
//!    the surface can rest ANYWHERE inside them,
//!  - SLIDING segments: outside stiction, quasi-static motion direction
//!    = sign(M_aero).
//!
//! The overbalance claim carrier (V-02b, sign/tendency ceiling): with
//! the hinge in the dossier prior band near the center of pressure, the
//! center-adjacent equilibrium is UNSTABLE — Orville's 'balanced too
//! near the center ... tendency to turn itself'. Quantitative levels
//! remain Estimated (A7a is the only promotion path).

use crate::Refusal;
use fs_blake3::hash_domain;

/// Sweep sample caps (absurd-input guards).
pub const MIN_SWEEP_SAMPLES: usize = 8;
/// Upper sweep cap.
pub const MAX_SWEEP_SAMPLES: usize = 4096;
/// Fixed bisection iterations (deterministic root refinement).
pub const BISECTION_ITERATIONS: u32 = 60;

/// Sweep specification.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SweepSpec {
    /// Lower deflection bound [rad].
    pub delta_min_rad: f64,
    /// Upper deflection bound [rad].
    pub delta_max_rad: f64,
    /// Sample count (uniform grid).
    pub samples: usize,
    /// Stiction threshold [N·m] (≥ 0).
    pub stiction_nm: f64,
}

/// Branch stability from the local moment slope.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BranchStability {
    /// dM/dδ < 0: restoring.
    Stable,
    /// dM/dδ > 0: self-driving.
    Unstable,
}

/// One equilibrium branch.
#[derive(Clone, Debug, PartialEq)]
pub struct Equilibrium {
    /// Deflection [rad].
    pub delta_rad: f64,
    /// Local moment slope [N·m/rad].
    pub slope_nm_per_rad: f64,
    /// Stability class.
    pub stability: BranchStability,
    /// Stable branch identity (digest of index + class + grid cell).
    pub branch_id: String,
}

/// A stiction interval (the surface rests anywhere inside).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StictionInterval {
    /// Lower edge [rad].
    pub lo_rad: f64,
    /// Upper edge [rad].
    pub hi_rad: f64,
}

/// A sliding segment (quasi-static drift direction).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlidingSegment {
    /// Lower edge [rad].
    pub lo_rad: f64,
    /// Upper edge [rad].
    pub hi_rad: f64,
    /// +1 (toward +δ) or −1.
    pub direction: i8,
}

/// The full free-control report (set-valued by construction).
#[derive(Clone, Debug, PartialEq)]
pub struct FreeControlReport {
    /// ALL equilibria found in the range (the SET).
    pub equilibria: Vec<Equilibrium>,
    /// Stiction interval set.
    pub stiction: Vec<StictionInterval>,
    /// Sliding segments between stiction/stops.
    pub sliding: Vec<SlidingSegment>,
    /// The sampled sweep (receipt).
    pub sweep: Vec<(f64, f64)>,
    /// The overbalance tendency: is the equilibrium NEAREST the center
    /// unstable (self-driving)? None if no equilibrium exists.
    pub self_driving_near_center: Option<bool>,
}

/// Run the free-control analysis on a hinge-moment function (the caller
/// supplies the coupled-solution sweep — fs-wing hinge_load per point).
///
/// # Errors
/// `sweep-spec-invalid` (bounds/order, sample caps at cap AND cap+1,
/// negative stiction); moment-function refusals pass through;
/// `sweep-moment-nonfinite`.
pub fn free_control_analysis(
    moment: &dyn Fn(f64) -> Result<f64, Refusal>,
    spec: &SweepSpec,
) -> Result<FreeControlReport, Refusal> {
    if !(spec.delta_min_rad.is_finite()
        && spec.delta_max_rad.is_finite()
        && spec.delta_min_rad < spec.delta_max_rad
        && spec.stiction_nm.is_finite()
        && spec.stiction_nm >= 0.0
        && (MIN_SWEEP_SAMPLES..=MAX_SWEEP_SAMPLES).contains(&spec.samples))
    {
        return Err(Refusal {
            code: "sweep-spec-invalid",
            message: format!("{spec:?}"),
            ranked_repairs: vec![format!(
                "ordered finite bounds; samples in [{MIN_SWEEP_SAMPLES}, {MAX_SWEEP_SAMPLES}]; stiction >= 0"
            )],
        });
    }
    let n = spec.samples;
    let dx = (spec.delta_max_rad - spec.delta_min_rad) / (n - 1) as f64;
    let mut sweep = Vec::with_capacity(n);
    for i in 0..n {
        let d = spec.delta_min_rad + dx * i as f64;
        let m = moment(d)?;
        if !m.is_finite() {
            return Err(Refusal {
                code: "sweep-moment-nonfinite",
                message: format!("M({d}) = {m:?}"),
                ranked_repairs: vec!["the moment function must be finite over the travel".into()],
            });
        }
        sweep.push((d, m));
    }
    // Equilibria: grid-exact zeros (each sample is the left endpoint of
    // exactly one cell, plus the final sample) and sign changes refined
    // by fixed-count bisection.
    let mut equilibria = Vec::new();
    let fd_slope = |d: f64| -> Result<f64, Refusal> {
        let h = dx / 64.0;
        Ok((moment(d + h)? - moment(d - h)?) / (2.0 * h))
    };
    for w in 0..(n - 1) {
        let (d0, m0) = sweep[w];
        let (d1, m1) = sweep[w + 1];
        if m0 == 0.0 {
            // Grid-exact root (measured: uniform grids land on rational
            // roots exactly): classify by the central FD slope.
            push_equilibrium(&mut equilibria, d0, fd_slope(d0)?, w);
            continue;
        }
        if m0 * m1 < 0.0 {
            let (mut lo, mut hi, mut mlo) = (d0, d1, m0);
            for _ in 0..BISECTION_ITERATIONS {
                let mid = 0.5 * (lo + hi);
                let mm = moment(mid)?;
                if mlo * mm <= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                    mlo = mm;
                }
            }
            let root = 0.5 * (lo + hi);
            // Central-FD slope at the root (h = one grid cell / 64).
            let h = dx / 64.0;
            let slope = (moment(root + h)? - moment(root - h)?) / (2.0 * h);
            push_equilibrium(&mut equilibria, root, slope, w);
        }
    }
    if sweep[n - 1].1 == 0.0 {
        let d = sweep[n - 1].0;
        push_equilibrium(&mut equilibria, d, fd_slope(d)?, n - 1);
    }
    // Stiction set: merge contiguous |M| <= stiction samples.
    let mut stiction = Vec::new();
    let mut open: Option<f64> = None;
    for &(d, m) in &sweep {
        if m.abs() <= spec.stiction_nm {
            if open.is_none() {
                open = Some(d);
            }
        } else if let Some(lo) = open.take() {
            stiction.push(StictionInterval {
                lo_rad: lo,
                hi_rad: d - dx,
            });
        }
    }
    if let Some(lo) = open {
        stiction.push(StictionInterval {
            lo_rad: lo,
            hi_rad: spec.delta_max_rad,
        });
    }
    // Sliding segments: contiguous runs with |M| > stiction, direction
    // from the moment sign (sign flips split segments via equilibria).
    let mut sliding = Vec::new();
    let mut run: Option<(f64, i8)> = None;
    for &(d, m) in &sweep {
        let dir = if m > spec.stiction_nm {
            1i8
        } else if m < -spec.stiction_nm {
            -1i8
        } else {
            0i8
        };
        match (run, dir) {
            (None, 0) => {}
            (None, _) => run = Some((d, dir)),
            (Some((lo, rd)), 0) => {
                sliding.push(SlidingSegment {
                    lo_rad: lo,
                    hi_rad: d - dx,
                    direction: rd,
                });
                run = None;
            }
            (Some((lo, rd)), _) if dir != rd => {
                sliding.push(SlidingSegment {
                    lo_rad: lo,
                    hi_rad: d - dx,
                    direction: rd,
                });
                run = Some((d, dir));
            }
            (Some(_), _) => {}
        }
    }
    if let Some((lo, rd)) = run {
        sliding.push(SlidingSegment {
            lo_rad: lo,
            hi_rad: spec.delta_max_rad,
            direction: rd,
        });
    }
    // Overbalance tendency: the equilibrium nearest δ = 0.
    let self_driving_near_center = equilibria
        .iter()
        .min_by(|a, b| {
            a.delta_rad
                .abs()
                .partial_cmp(&b.delta_rad.abs())
                .unwrap_or(core::cmp::Ordering::Equal)
        })
        .map(|e| e.stability == BranchStability::Unstable);
    Ok(FreeControlReport {
        equilibria,
        stiction,
        sliding,
        sweep,
        self_driving_near_center,
    })
}

fn push_equilibrium(out: &mut Vec<Equilibrium>, delta: f64, slope: f64, cell: usize) {
    let stability = if slope > 0.0 {
        BranchStability::Unstable
    } else {
        BranchStability::Stable
    };
    let mut p = Vec::new();
    p.extend_from_slice(&(cell as u64).to_le_bytes());
    p.push(match stability {
        BranchStability::Stable => 0,
        BranchStability::Unstable => 1,
    });
    p.extend_from_slice(&(out.len() as u64).to_le_bytes());
    let branch_id = hash_domain("org.frankensim.fs-flyer.free-control-branch.v1", &p).to_hex();
    out.push(Equilibrium {
        delta_rad: delta,
        slope_nm_per_rad: slope,
        stability,
        branch_id,
    });
}
