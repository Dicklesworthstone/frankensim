//! À-trous wavelet denoiser with feature guides (plan §10.5) — CLEARLY
//! LABELED AS BIASED: the output type carries a mandatory provenance tag,
//! so a denoised image can never masquerade as a converged render in a
//! comparison (the honesty tag is part of the API, not a docstring).
//!
//! Method: iterated 5×5 B3-spline à-trous convolution with edge-stopping
//! weights from optional albedo/normal guides (SVGF-lineage, single
//! frame). Deterministic: fixed traversal order, pure f64 accumulation.

use crate::ImgError;

/// Why these pixels are what they are — the mandatory honesty tag.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PixelProvenance {
    /// Untouched estimator output (may be noisy; unbiased).
    RawEstimate,
    /// Smoothed by the à-trous denoiser: BIASED — usable for previews and
    /// perceptual metrics, NEVER as ground truth in a comparison.
    BiasedDenoised {
        /// À-trous iterations applied.
        iterations: u8,
    },
}

/// A single-channel image plane with mandatory provenance.
#[derive(Debug, Clone, PartialEq)]
pub struct LabeledPlane {
    /// Width.
    pub width: usize,
    /// Height.
    pub height: usize,
    /// Row-major samples.
    pub data: Vec<f32>,
    /// The honesty tag.
    pub provenance: PixelProvenance,
}

/// Edge-stopping parameters.
#[derive(Debug, Clone, Copy)]
pub struct DenoiseParams {
    /// À-trous iterations (hole size doubles each pass).
    pub iterations: u8,
    /// Color-difference sigma.
    pub sigma_color: f32,
    /// Albedo-difference sigma (used when an albedo guide is given).
    pub sigma_albedo: f32,
}

impl Default for DenoiseParams {
    fn default() -> Self {
        DenoiseParams {
            iterations: 3,
            sigma_color: 0.25,
            sigma_albedo: 0.1,
        }
    }
}

const B3: [f32; 5] = [1.0 / 16.0, 4.0 / 16.0, 6.0 / 16.0, 4.0 / 16.0, 1.0 / 16.0];

/// Denoise one plane, optionally guided by an albedo plane of the same
/// shape. The output is PERMANENTLY tagged `BiasedDenoised`.
///
/// # Errors
/// [`ImgError::Shape`] on plane-shape disagreement.
pub fn atrous_denoise(
    noisy: &LabeledPlane,
    albedo: Option<&LabeledPlane>,
    params: &DenoiseParams,
) -> Result<LabeledPlane, ImgError> {
    let n = noisy.width * noisy.height;
    if noisy.data.len() != n {
        return Err(ImgError::Shape {
            expected: n,
            got: noisy.data.len(),
            context: "noisy plane",
        });
    }
    if let Some(a) = albedo
        && (a.width != noisy.width || a.height != noisy.height || a.data.len() != n)
    {
        return Err(ImgError::Shape {
            expected: n,
            got: a.data.len(),
            context: "albedo guide shape",
        });
    }
    let (w, h) = (noisy.width as isize, noisy.height as isize);
    let mut current = noisy.data.clone();
    for it in 0..params.iterations {
        let step = 1isize << it;
        let mut next = vec![0.0f32; current.len()];
        for y in 0..h {
            for x in 0..w {
                let center = current[(y * w + x) as usize];
                let center_albedo = albedo.map(|a| a.data[(y * w + x) as usize]);
                let mut acc = 0.0f64;
                let mut wsum = 0.0f64;
                for (kj, &wy) in B3.iter().enumerate() {
                    for (ki, &wx) in B3.iter().enumerate() {
                        let sx = (x + (ki as isize - 2) * step).clamp(0, w - 1);
                        let sy = (y + (kj as isize - 2) * step).clamp(0, h - 1);
                        let sample = current[(sy * w + sx) as usize];
                        let mut weight = f64::from(wx * wy);
                        let dc = f64::from(sample - center) / f64::from(params.sigma_color);
                        weight *= (-dc * dc).exp();
                        if let (Some(ca), Some(a)) = (center_albedo, albedo) {
                            let da = f64::from(a.data[(sy * w + sx) as usize] - ca)
                                / f64::from(params.sigma_albedo);
                            weight *= (-da * da).exp();
                        }
                        acc += weight * f64::from(sample);
                        wsum += weight;
                    }
                }
                next[(y * w + x) as usize] = if wsum > 0.0 {
                    (acc / wsum) as f32
                } else {
                    center
                };
            }
        }
        current = next;
    }
    Ok(LabeledPlane {
        width: noisy.width,
        height: noisy.height,
        data: current,
        provenance: PixelProvenance::BiasedDenoised {
            iterations: params.iterations,
        },
    })
}

/// Mean squared error between planes (the improvement metric).
///
/// # Errors
/// [`ImgError::Shape`] on length disagreement.
pub fn mse(a: &[f32], b: &[f32]) -> Result<f64, ImgError> {
    if a.len() != b.len() {
        return Err(ImgError::Shape {
            expected: a.len(),
            got: b.len(),
            context: "mse operands",
        });
    }
    if a.is_empty() {
        return Ok(0.0);
    }
    let sum: f64 = a
        .iter()
        .zip(b)
        .map(|(&x, &y)| {
            let d = f64::from(x) - f64::from(y);
            d * d
        })
        .sum();
    Ok(sum / a.len() as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_is_permanently_labeled_biased() {
        let noisy = LabeledPlane {
            width: 4,
            height: 4,
            data: vec![0.5; 16],
            provenance: PixelProvenance::RawEstimate,
        };
        let out = atrous_denoise(&noisy, None, &DenoiseParams::default()).unwrap();
        assert_eq!(
            out.provenance,
            PixelProvenance::BiasedDenoised { iterations: 3 }
        );
        // A constant image stays constant (partition of unity).
        for &v in &out.data {
            assert!((v - 0.5).abs() < 1e-6);
        }
    }
}
