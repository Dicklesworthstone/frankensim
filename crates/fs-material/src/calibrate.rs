//! Parameter calibration from experimental stress–strain data with
//! uncertainty envelopes: segmented least squares for bilinear
//! (elastic + hardening) laws — the round-trip-tested v0 of the
//! calibration story. Fits recover (E, σ_y, H) with standard errors from
//! residual variance; the envelope feeds the model card.

use crate::MaterialError;
use fs_math::det;

/// A calibrated bilinear fit with uncertainty.
#[derive(Debug, Clone, PartialEq)]
pub struct CalibrationFit {
    /// Elastic modulus E (slope of the first segment).
    pub youngs: f64,
    /// Standard error of E.
    pub youngs_se: f64,
    /// Post-yield modulus (slope of the second segment).
    pub post_yield: f64,
    /// Standard error of the post-yield slope.
    pub post_yield_se: f64,
    /// Yield stress (segment-intersection ordinate).
    pub yield_stress: f64,
    /// Root-mean-square residual (Pa).
    pub rms_residual: f64,
    /// Index of the fitted breakpoint.
    pub break_index: usize,
}

/// Least-squares line through points (through the origin when `origin`).
fn fit_line(pts: &[(f64, f64)], through_origin: bool) -> (f64, f64, f64) {
    // Returns (slope, intercept, sse).
    let n = pts.len() as f64;
    if through_origin {
        let sxy: f64 = pts.iter().map(|(x, y)| x * y).sum();
        let sxx: f64 = pts.iter().map(|(x, _)| x * x).sum();
        let slope = sxy / sxx;
        let sse: f64 = pts.iter().map(|(x, y)| (y - slope * x).powi(2)).sum();
        return (slope, 0.0, sse);
    }
    let sx: f64 = pts.iter().map(|(x, _)| x).sum();
    let sy: f64 = pts.iter().map(|(_, y)| y).sum();
    let sxx: f64 = pts.iter().map(|(x, _)| x * x).sum();
    let sxy: f64 = pts.iter().map(|(x, y)| x * y).sum();
    let denom = n * sxx - sx * sx;
    let slope = (n * sxy - sx * sy) / denom;
    let intercept = (sy - slope * sx) / n;
    let sse: f64 = pts
        .iter()
        .map(|(x, y)| (y - slope * x - intercept).powi(2))
        .sum();
    (slope, intercept, sse)
}

fn slope_se(pts: &[(f64, f64)], sse: f64) -> f64 {
    let n = pts.len() as f64;
    if n <= 2.0 {
        return f64::INFINITY;
    }
    let mean_x: f64 = pts.iter().map(|(x, _)| x).sum::<f64>() / n;
    let sxx: f64 = pts.iter().map(|(x, _)| (x - mean_x).powi(2)).sum();
    det::sqrt(sse / (n - 2.0) / sxx.max(f64::MIN_POSITIVE))
}

/// Fit a bilinear (elastic/hardening) law to monotonic uniaxial data by
/// scanning every admissible breakpoint and keeping the least-SSE split.
///
/// # Errors
/// [`MaterialError::Calibration`] for fewer than 6 points or degenerate
/// (non-increasing-strain) data.
pub fn calibrate_bilinear(data: &[(f64, f64)]) -> Result<CalibrationFit, MaterialError> {
    if data.len() < 6 {
        return Err(MaterialError::Calibration {
            what: format!("need at least 6 points, got {}", data.len()),
        });
    }
    if data.windows(2).any(|w| w[1].0 <= w[0].0) {
        return Err(MaterialError::Calibration {
            what: "strains must be strictly increasing".to_string(),
        });
    }
    let mut best: Option<(usize, f64)> = None;
    for k in 3..data.len() - 2 {
        let (_, _, sse1) = fit_line(&data[..k], true);
        let (_, _, sse2) = fit_line(&data[k..], false);
        let total = sse1 + sse2;
        if best.is_none_or(|(_, b)| total < b) {
            best = Some((k, total));
        }
    }
    let (k, _) = best.expect("len >= 6 guarantees candidates");
    let seg1 = &data[..k];
    let seg2 = &data[k..];
    let (e, _, sse1) = fit_line(seg1, true);
    let (h, c, sse2) = fit_line(seg2, false);
    // Yield point: intersection of σ = E·ε with σ = h·ε + c.
    let eps_y = c / (e - h);
    let sigma_y = e * eps_y;
    let n = data.len() as f64;
    let rms = det::sqrt((sse1 + sse2) / n);
    Ok(CalibrationFit {
        youngs: e,
        youngs_se: slope_se(seg1, sse1),
        post_yield: h,
        post_yield_se: slope_se(seg2, sse2),
        yield_stress: sigma_y,
        rms_residual: rms,
        break_index: k,
    })
}
