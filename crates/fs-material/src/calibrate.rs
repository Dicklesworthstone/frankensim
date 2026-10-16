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

fn slope_se(pts: &[(f64, f64)], sse: f64, through_origin: bool) -> f64 {
    let n = pts.len() as f64;
    if through_origin {
        // y = β·x: Var(β̂) = σ²/Σx² with n−1 degrees of freedom (one parameter).
        if n <= 1.0 {
            return f64::INFINITY;
        }
        let sxx_raw: f64 = pts.iter().map(|(x, _)| x * x).sum();
        det::sqrt(sse / (n - 1.0) / sxx_raw.max(f64::MIN_POSITIVE))
    } else {
        // y = β·x + α: Var(β̂) = σ²/Σ(x−x̄)² with n−2 degrees of freedom.
        if n <= 2.0 {
            return f64::INFINITY;
        }
        let mean_x: f64 = pts.iter().map(|(x, _)| x).sum::<f64>() / n;
        let sxx: f64 = pts.iter().map(|(x, _)| (x - mean_x).powi(2)).sum();
        det::sqrt(sse / (n - 2.0) / sxx.max(f64::MIN_POSITIVE))
    }
}

/// Fit a bilinear (elastic/hardening) law to monotonic uniaxial data by
/// scanning every admissible breakpoint and keeping the least-SSE split.
///
/// # Errors
/// [`MaterialError::Calibration`] for fewer than 6 points, non-finite or
/// non-increasing observations, non-finite fits, or segments whose yield
/// intersection is unidentifiable.
pub fn calibrate_bilinear(data: &[(f64, f64)]) -> Result<CalibrationFit, MaterialError> {
    if data.len() < 6 {
        return Err(MaterialError::Calibration {
            what: format!("need at least 6 points, got {}", data.len()),
        });
    }
    for (index, &(strain, stress)) in data.iter().enumerate() {
        if !strain.is_finite() || !stress.is_finite() {
            return Err(MaterialError::Calibration {
                what: format!(
                    "sample {index} must have finite strain and stress, got ({strain}, {stress})"
                ),
            });
        }
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
        if total.is_finite() && best.is_none_or(|(_, b)| total < b) {
            best = Some((k, total));
        }
    }
    let (k, _) = best.ok_or_else(|| MaterialError::Calibration {
        what: "no finite bilinear fit candidate".to_string(),
    })?;
    let seg1 = &data[..k];
    let seg2 = &data[k..];
    let (e, _, sse1) = fit_line(seg1, true);
    let (h, c, sse2) = fit_line(seg2, false);
    // Yield point: intersection of σ = E·ε with σ = h·ε + c.
    let slope_gap = e - h;
    if !e.is_finite()
        || !h.is_finite()
        || !c.is_finite()
        || !sse1.is_finite()
        || !sse2.is_finite()
        || !slope_gap.is_finite()
        || slope_gap == 0.0
    {
        return Err(MaterialError::Calibration {
            what: "bilinear segments must have distinct finite slopes".to_string(),
        });
    }
    let eps_y = c / slope_gap;
    let sigma_y = e * eps_y;
    let n = data.len() as f64;
    let rms = det::sqrt((sse1 + sse2) / n);
    let fit = CalibrationFit {
        youngs: e,
        youngs_se: slope_se(seg1, sse1, true),
        post_yield: h,
        post_yield_se: slope_se(seg2, sse2, false),
        yield_stress: sigma_y,
        rms_residual: rms,
        break_index: k,
    };
    if [
        fit.youngs,
        fit.youngs_se,
        fit.post_yield,
        fit.post_yield_se,
        fit.yield_stress,
        fit.rms_residual,
    ]
    .iter()
    .any(|value| !value.is_finite())
    {
        return Err(MaterialError::Calibration {
            what: "bilinear fit produced non-finite parameters".to_string(),
        });
    }
    Ok(fit)
}
