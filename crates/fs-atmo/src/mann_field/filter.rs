//! Exact uniform spatial averages of the finite Fourier representation.
//!
//! For half-edge vectors `a_j`, integrate `u(x + sum(s_j a_j))` over
//! independent `s_j` uniform on [-1,1]. A mode receives
//! `G(k) = product_j sinc(k dot a_j)`, with `sinc(0)=1`. This includes rotated
//! boxes, parallelepipeds, faces, lines and points. It is an observation
//! operator, not an anti-aliasing guarantee or a new turbulence spectrum.
//! Both polarizations receive the same gain, so the analytic divergence stays
//! zero; covariance receives G squared. No lost energy is rescaled back in.

use super::{MannParams, MannQuadrature, MannSpectralModel, Refusal, det, dot, refuse};

/// Prepared scalar gains, bound to the complete immutable source model.
/// Preparation is O(modes), sampling allocates nothing and draws no randomness.
/// Numerical gain storage is at most 8*MAX_MODES bytes (plus fixed metadata).
#[derive(Clone, Debug, PartialEq)]
pub struct MannBoxFilter {
    geometry_seed: u64,
    parameters: MannParams,
    quadrature: MannQuadrature,
    advection_mps: [f64; 3],
    half_edges_m: [[f64; 3]; 3],
    gains: Vec<f64>,
}

fn error(code: &'static str, message: &str) -> Refusal {
    refuse(
        code,
        message.into(),
        "use the same admitted Mann model and finite physical half edges within the gain-storage budget",
    )
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(error(
            "mann-filter-cancelled",
            "filter preparation cancelled",
        ))
    } else {
        Ok(())
    }
}

fn sinc(x: f64) -> f64 {
    if x.abs() < 1e-4 {
        // Even series avoids sin(x)/x cancellation at a collapsed dimension.
        let x2 = x * x;
        1.0 + x2 * (-1.0 / 6.0 + x2 * (1.0 / 120.0 - x2 / 5040.0))
    } else {
        det::sin(x) / x
    }
}

impl MannBoxFilter {
    /// Prepare a normalized uniform average in the field's physical frame.
    /// Zero half edges collapse dimensions; all-zero edges use the exact
    /// original point-sample arithmetic and require no gain storage.
    /// The budget covers allocated numerical gain capacity, not fixed metadata.
    /// Cancellation is polled before allocation and between every mode.
    pub fn prepare(
        model: &MannSpectralModel,
        half_edges_m: [[f64; 3]; 3],
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Refusal> {
        poll(&mut cancelled)?;
        if !half_edges_m.iter().flatten().all(|v| v.is_finite()) {
            return Err(error(
                "mann-filter-invalid",
                "filter half edges must be finite",
            ));
        }
        let point = half_edges_m.iter().flatten().all(|&v| v == 0.0);
        let count = if point { 0 } else { model.mode_count() };
        let required = count * core::mem::size_of::<f64>();
        if required > maximum_workspace_bytes {
            return Err(error(
                "mann-filter-budget",
                "filter gain storage exceeds its byte budget",
            ));
        }
        let mut gains = Vec::new();
        gains
            .try_reserve_exact(count)
            .map_err(|_| error("mann-filter-allocation", "unable to reserve filter gains"))?;
        if gains.capacity() * core::mem::size_of::<f64>() > maximum_workspace_bytes {
            return Err(error(
                "mann-filter-budget",
                "allocated gain capacity exceeds its byte budget",
            ));
        }
        if !point {
            for mode in &model.modes {
                poll(&mut cancelled)?;
                let mut gain = 1.0;
                for edge in half_edges_m {
                    let phase = dot(mode.k, edge);
                    if !phase.is_finite() {
                        return Err(error(
                            "mann-filter-invalid",
                            "filter phase is not representable",
                        ));
                    }
                    gain *= sinc(phase);
                }
                if !gain.is_finite() {
                    return Err(error(
                        "mann-filter-invalid",
                        "filter gain is not representable",
                    ));
                }
                gains.push(gain);
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            geometry_seed: model.geometry_seed,
            parameters: model.parameters,
            quadrature: model.quadrature,
            advection_mps: model.advection_mps,
            half_edges_m,
            gains,
        })
    }

    #[must_use]
    pub const fn half_edges_m(&self) -> [[f64; 3]; 3] {
        self.half_edges_m
    }

    #[must_use]
    pub fn workspace_bytes(&self) -> usize {
        self.gains.capacity() * core::mem::size_of::<f64>()
    }

    pub(super) fn gains(&self, model: &MannSpectralModel) -> Result<Option<&[f64]>, Refusal> {
        // These are every input of the immutable deterministic model builder;
        // no modal vector clone, hash collision or O(modes) comparison is needed.
        if self.geometry_seed != model.geometry_seed
            || self.parameters != model.parameters
            || self.quadrature != model.quadrature
            || self.advection_mps != model.advection_mps
        {
            return Err(error(
                "mann-filter-source-mismatch",
                "filter belongs to a different Mann spectrum or frame",
            ));
        }
        Ok(if self.gains.is_empty() {
            None
        } else {
            Some(&self.gains)
        })
    }
}

#[cfg(test)]
#[path = "filter/tests.rs"]
mod tests;
