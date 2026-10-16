//! E4.2c — extracted added-mass CROSS terms (bead wf-root-guzez.5.5).
//! OFFLINE extraction over a pinned (canard deflection × height ×
//! warp) grid, an interpolation artifact with ANALYTIC derivatives,
//! and the discrepancy battery vs the AnalyticStrip baseline.
//!
//! Tier honesty (the AddedMassMode ladder doctrine): fs-wing does not
//! yet expose the surface potential, so this extraction is
//! `config-resolved-strip-v1` — the rigid↔canard cross column of the
//! REAL strip assembly evaluated on the DEFORMED configuration
//! (canard normal rotated by its deflection, tip strips rotated by
//! warp, moment arms shifted by height). The artifact DECLARES its
//! tier; a potential-integrating panel tier slots into the same
//! artifact schema when fs-wing grows a Φ surface. The declared
//! approximation is carried, never hidden.
//!
//! Runtime law: consumers evaluate the artifact's trilinear
//! interpolant, whose derivatives are ANALYTIC (exact closed form of
//! the interpolant). Runtime finite-differencing is FORBIDDEN and the
//! entry point that would do it refuses.

use crate::addedmass::{Strip, assemble_analytic_strip};
use crate::{Refusal, refuse};
use fs_blake3::hash_domain;
use fs_math::det;

/// Extraction algorithm version (enters the artifact id).
pub const EXTRACTION_VERSION: &str = "config-resolved-strip-v1";

/// Grid-axis cap (per axis).
pub const MAX_AXIS: usize = 16;

/// The pinned extraction grid.
#[derive(Clone, Debug, PartialEq)]
pub struct ExtractionGrid {
    /// Canard deflection samples [rad], strictly ascending.
    pub deflection_rad: Vec<f64>,
    /// Height samples [m], strictly ascending.
    pub height_m: Vec<f64>,
    /// Warp samples [rad], strictly ascending.
    pub warp_rad: Vec<f64>,
}

impl ExtractionGrid {
    fn admit(&self) -> Result<(), Refusal> {
        for (name, axis) in [
            ("deflection", &self.deflection_rad),
            ("height", &self.height_m),
            ("warp", &self.warp_rad),
        ] {
            if axis.len() < 2 || axis.len() > MAX_AXIS {
                return Err(refuse(
                    "extraction-grid-invalid",
                    format!("{name}: {} samples outside [2, {MAX_AXIS}]", axis.len()),
                    "a pinned grid is small and declared",
                ));
            }
            if axis.windows(2).any(|w| !(w[1] > w[0])) || axis.iter().any(|v| !v.is_finite()) {
                return Err(refuse(
                    "extraction-grid-invalid",
                    format!("{name} axis not strictly ascending finite"),
                    "sort and dedupe the pins",
                ));
            }
        }
        Ok(())
    }
}

/// The registered 1903-class strip set at a deformed configuration.
/// Deflection rotates the canard plate normal; warp rotates the tip
/// strips oppositely; height shifts every vertical moment arm.
#[must_use]
pub fn strips_at(deflection_rad: f64, height_m: f64, warp_rad: f64) -> Vec<Strip> {
    let (sd, cd) = (det::sin(deflection_rad), det::cos(deflection_rad));
    let (sw, cw) = (det::sin(warp_rad), det::cos(warp_rad));
    vec![
        Strip {
            name: "canard",
            chord_m: 0.61,
            span_m: 3.66,
            position_m: [2.23, 0.0, -height_m],
            normal: [sd, 0.0, cd],
            control_coord: Some(0),
            control_gain: 0.9,
        },
        Strip {
            name: "wing-center",
            chord_m: 1.981,
            span_m: 6.0,
            position_m: [0.0, 0.0, -height_m],
            normal: [0.0, 0.0, 1.0],
            control_coord: None,
            control_gain: 0.0,
        },
        Strip {
            name: "wing-tip-left",
            chord_m: 1.981,
            span_m: 3.0,
            position_m: [0.0, -4.6, -height_m],
            normal: [0.0, -sw, cw],
            control_coord: None,
            control_gain: 0.0,
        },
        Strip {
            name: "wing-tip-right",
            chord_m: 1.981,
            span_m: 3.0,
            position_m: [0.0, 4.6, -height_m],
            normal: [0.0, sw, cw],
            control_coord: None,
            control_gain: 0.0,
        },
    ]
}

/// The extracted artifact: the rigid↔canard cross column m_rc[0..6]
/// per grid node, identity-hashed.
#[derive(Clone, Debug, PartialEq)]
pub struct PanelExtractedCrossTermsV1 {
    /// Declared extraction tier.
    pub extraction_tier: &'static str,
    /// The pinned grid.
    pub grid: ExtractionGrid,
    /// values[d][h][w] = 6-vector cross column [kg·m-class].
    pub values: Vec<Vec<Vec<[f64; 6]>>>,
    /// Worst |extracted − baseline| per component, over the grid
    /// (the recorded discrepancy vs AnalyticStrip at the NOMINAL
    /// configuration — data, never forced).
    pub baseline_discrepancy: [f64; 6],
    /// Identity hash over tier + grid + values.
    pub artifact_id: String,
}

/// Air density the extraction pins (Dec-17 class).
pub const EXTRACTION_RHO: f64 = 1.294;

fn cross_column(deflection: f64, height: f64, warp: f64) -> Result<[f64; 6], Refusal> {
    let strips = strips_at(deflection, height, warp);
    let loads = assemble_analytic_strip(EXTRACTION_RHO, &strips, 1, &[0.0; 6], &[0.0])?;
    Ok(loads.m_added_rc[0])
}

/// Run the OFFLINE extraction over the pinned grid.
///
/// # Errors
/// Grid refusals; strip-assembly refusals pass through.
pub fn extract(grid: ExtractionGrid) -> Result<PanelExtractedCrossTermsV1, Refusal> {
    grid.admit()?;
    let baseline = cross_column(0.0, 0.0, 0.0)?;
    let mut values = Vec::with_capacity(grid.deflection_rad.len());
    let mut disc = [0.0f64; 6];
    for d in &grid.deflection_rad {
        let mut plane = Vec::with_capacity(grid.height_m.len());
        for h in &grid.height_m {
            let mut row = Vec::with_capacity(grid.warp_rad.len());
            for w in &grid.warp_rad {
                let col = cross_column(*d, *h, *w)?;
                for k in 0..6 {
                    disc[k] = disc[k].max((col[k] - baseline[k]).abs());
                }
                row.push(col);
            }
            plane.push(row);
        }
        values.push(plane);
    }
    let mut b = EXTRACTION_VERSION.as_bytes().to_vec();
    for axis in [&grid.deflection_rad, &grid.height_m, &grid.warp_rad] {
        for v in axis {
            b.extend_from_slice(&v.to_bits().to_le_bytes());
        }
    }
    for plane in &values {
        for row in plane {
            for col in row {
                for v in col {
                    b.extend_from_slice(&v.to_bits().to_le_bytes());
                }
            }
        }
    }
    let artifact_id = hash_domain("org.frankensim.wf.panel-crossterms.v1", &b).to_hex();
    Ok(PanelExtractedCrossTermsV1 {
        extraction_tier: EXTRACTION_VERSION,
        grid,
        values,
        baseline_discrepancy: disc,
        artifact_id,
    })
}

/// One axis's bracketing index + normalized coordinate.
fn bracket(axis: &[f64], x: f64) -> Result<(usize, f64), Refusal> {
    let lo = *axis.first().expect("admitted axis");
    let hi = *axis.last().expect("admitted axis");
    if !(x >= lo && x <= hi) {
        return Err(refuse(
            "crossterm-out-of-domain",
            format!("{x} outside [{lo}, {hi}]"),
            "the artifact never extrapolates",
        ));
    }
    let mut i = 0;
    while i + 2 < axis.len() && x > axis[i + 1] {
        i += 1;
    }
    let t = (x - axis[i]) / (axis[i + 1] - axis[i]);
    Ok((i, t))
}

/// Trilinear evaluation with ANALYTIC derivatives: the value and its
/// exact gradient d/d(deflection, height, warp) of the interpolant.
///
/// # Errors
/// `crossterm-out-of-domain` (AT the bounds admits; beyond refuses).
pub fn eval_with_derivatives(
    art: &PanelExtractedCrossTermsV1,
    deflection_rad: f64,
    height_m: f64,
    warp_rad: f64,
) -> Result<([f64; 6], [[f64; 3]; 6]), Refusal> {
    let (di, dt) = bracket(&art.grid.deflection_rad, deflection_rad)?;
    let (hi, ht) = bracket(&art.grid.height_m, height_m)?;
    let (wi, wt) = bracket(&art.grid.warp_rad, warp_rad)?;
    let dd = art.grid.deflection_rad[di + 1] - art.grid.deflection_rad[di];
    let dh = art.grid.height_m[hi + 1] - art.grid.height_m[hi];
    let dw = art.grid.warp_rad[wi + 1] - art.grid.warp_rad[wi];
    let mut value = [0.0f64; 6];
    let mut grad = [[0.0f64; 3]; 6];
    for k in 0..6 {
        let c = |a: usize, b: usize, cc: usize| art.values[di + a][hi + b][wi + cc][k];
        // Trilinear basis and its exact partials.
        let lerp3 = |td: f64, th: f64, tw: f64| -> f64 {
            let c00 = c(0, 0, 0) * (1.0 - td) + c(1, 0, 0) * td;
            let c01 = c(0, 0, 1) * (1.0 - td) + c(1, 0, 1) * td;
            let c10 = c(0, 1, 0) * (1.0 - td) + c(1, 1, 0) * td;
            let c11 = c(0, 1, 1) * (1.0 - td) + c(1, 1, 1) * td;
            let c0 = c00 * (1.0 - th) + c10 * th;
            let c1 = c01 * (1.0 - th) + c11 * th;
            c0 * (1.0 - tw) + c1 * tw
        };
        value[k] = lerp3(dt, ht, wt);
        // Analytic partials of the trilinear form (exact, closed
        // form — NOT finite differences).
        grad[k][0] = (lerp3(1.0, ht, wt) - lerp3(0.0, ht, wt)) / dd;
        grad[k][1] = (lerp3(dt, 1.0, wt) - lerp3(dt, 0.0, wt)) / dh;
        grad[k][2] = (lerp3(dt, ht, 1.0) - lerp3(dt, ht, 0.0)) / dw;
    }
    Ok((value, grad))
}

/// The FORBIDDEN runtime path: finite-differencing the artifact at
/// runtime. It exists only to refuse — the analytic derivatives above
/// are the lawful path, and a runtime FD would silently disagree with
/// them at cell boundaries.
///
/// # Errors
/// Always `runtime-fd-forbidden`.
pub fn runtime_finite_difference(
    _art: &PanelExtractedCrossTermsV1,
    _deflection_rad: f64,
    _height_m: f64,
    _warp_rad: f64,
) -> Result<[[f64; 3]; 6], Refusal> {
    Err(refuse(
        "runtime-fd-forbidden",
        "runtime finite-differencing of the cross-term artifact".into(),
        "use eval_with_derivatives — the interpolant's derivatives are analytic",
    ))
}
