//! Odd-extension PSE at one stationary flat stress-free wall.
//!
//! For physical supports i,j, r=K(x_i-x_j), m=K(x_i-x_j*) and
//! omega=Gamma/A, the image contribution is -m(A_i Gamma_j+A_j Gamma_i).
//! Combine it with the real pair as a positive exchange with coefficient
//! r-m, plus diagonal losses 2m A_j Gamma_i and 2m A_i Gamma_j.
//! Each self-image contributes -2 K(x_i-x_i*) A_i Gamma_i.
//!
//! A palindrome of exact exchange and loss maps is second order at fixed
//! sites. Every map preserves nonnegativity and contracts area-weighted
//! enstrophy. Bulk circulation can change through the wall; no image is a
//! physical quadrature site or an additional returned vortex.

use super::{PairKernel, PseDiffusion, PseKernel, exchange_fraction, scaled_distance2};
use crate::planar_wall::PlanarWall2;
use crate::{AdmittedBudget, Cx, VortexParticle, VpmError, WorkMeter};
use fs_math::det;

pub(super) fn pair_work(n: usize, steps: usize, sweeps: u64) -> Result<u64, VpmError> {
    let overflow = || VpmError::PlanOverflow {
        resource: "wall diffusion pair evaluations",
    };
    let n = u64::try_from(n).map_err(|_| overflow())?;
    let steps = u64::try_from(steps).map_err(|_| overflow())?;
    n.checked_mul(n)
        .and_then(|v| v.checked_mul(2))
        .and_then(|v| v.checked_mul(steps))
        .and_then(|v| v.checked_mul(sweeps))
        .ok_or_else(overflow)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn symmetric_sweep(
    state: &mut [VortexParticle],
    areas: &[f64],
    law: PseDiffusion,
    dt: f64,
    pair_fraction: f64,
    visits: &mut u64,
    meter: &mut WorkMeter,
    accountant: &mut AdmittedBudget<'_>,
    cx: &Cx<'_>,
    wall: PlanarWall2,
) -> Result<(), VpmError> {
    let kernel = PairKernel::new(law, dt, pair_fraction);
    let n = state.len();
    let mut visit = |i: usize, j: usize, reverse: bool| {
        let rates = pair_rates(&kernel, wall, state[i], state[j], [areas[i], areas[j]]);
        let mut apply = |loss: bool| {
            meter.unit(accountant, cx, "fs-vpm.wall-diffusion-pair", || {
                *visits = visits.checked_add(1).ok_or(VpmError::PlanOverflow {
                    resource: "executed wall diffusion pairs",
                })?;
                if loss {
                    state[i].circulation = decay(state[i].circulation, rates.loss[0], i)?;
                    if i != j {
                        state[j].circulation = decay(state[j].circulation, rates.loss[1], j)?;
                    }
                    Ok(())
                } else {
                    exchange_fraction(state, areas, rates.exchange, i, j)
                }
            })
        };
        if i == j {
            apply(true)
        } else {
            apply(reverse)?;
            apply(!reverse)
        }
    };
    for i in 0..n {
        visit(i, i, false)?;
        for j in i + 1..n {
            visit(i, j, false)?;
        }
    }
    for i in (0..n).rev() {
        for j in (i + 1..n).rev() {
            visit(i, j, true)?;
        }
        visit(i, i, true)?;
    }
    meter.boundary(accountant, cx, "fs-vpm.wall-diffusion-sweep-complete")
}

#[derive(Clone, Copy)]
pub(super) struct Rates {
    pub(super) exchange: f64,
    pub(super) loss: [f64; 2],
}

pub(super) fn pair_rates(
    kernel: &PairKernel,
    wall: PlanarWall2,
    a: VortexParticle,
    b: VortexParticle,
    areas: [f64; 2],
) -> Rates {
    let real2 = scaled_distance2(a, b, kernel.width);
    if real2.is_infinite() || (kernel.shape == PseKernel::CompactWendland && real2 >= 1.0) {
        return Rates {
            exchange: 0.0,
            loss: [0.0; 2],
        };
    }
    let ha = wall.reflection(a.pos).0;
    let hb = wall.reflection(b.pos).0;
    // |x_i-x_j*|² - |x_i-x_j|² = 4 h_i h_j. Form the difference
    // directly: subtracting two rounded squared distances loses weak
    // exchange beside the wall.
    let log_extra = || det::ln(4.0) + det::ln(ha) + det::ln(hb) - 2.0 * det::ln(kernel.width);
    let mut extra = (2.0 * (ha / kernel.width)) * (2.0 * (hb / kernel.width));
    if !extra.is_normal() {
        extra = det::exp(log_extra());
    }
    match kernel.shape {
        PseKernel::Gaussian => {
            let retained = -det::expm1(-extra);
            let image2 = real2 + extra;
            rates(
                kernel,
                areas,
                det::exp(-real2) * retained,
                || {
                    -real2
                        + if retained > 0.0 {
                            det::ln(retained)
                        } else {
                            log_extra()
                        }
                },
                det::exp(-image2),
                || -image2,
            )
        }
        PseKernel::CompactWendland => {
            let q = det::sqrt(real2);
            let edge = (1.0 - real2) / (1.0 + q);
            let square = edge * edge;
            let real_weight = square * square * (1.0 + 4.0 * q);
            let image2 = real2 + extra;
            if image2 >= 1.0 {
                return rates(
                    kernel,
                    areas,
                    real_weight,
                    || 4.0 * det::ln(edge) + det::ln(1.0 + 4.0 * q),
                    0.0,
                    || f64::NEG_INFINITY,
                );
            }
            let qi = det::sqrt(image2);
            let ei = (1.0 - image2) / (1.0 + qi);
            let sum = q + qi;
            let gap = if sum > 0.0 { extra / sum } else { 0.0 };
            // Integral of 20 q (1-q)^3 from q to qi, written only
            // with nonnegative terms. Direct subtraction of the compact
            // weights would erase a thin wall-adjacent exchange.
            let e2 = ei * ei;
            let g2 = gap * gap;
            let polynomial = q * (ei * e2 + 1.5 * e2 * gap + ei * g2 + 0.25 * gap * g2)
                + gap * (0.5 * ei * e2 + 0.5 * e2 * gap + 0.25 * ei * g2 + 0.05 * gap * g2);
            rates(
                kernel,
                areas,
                20.0 * gap * polynomial,
                || {
                    if gap > 0.0 && polynomial > 0.0 {
                        det::ln(20.0) + det::ln(gap) + det::ln(polynomial)
                    } else {
                        // In the unresolved-gap limit d(phi)/d(q²)
                        // = -10(1-q)^3; retain its representable rate.
                        det::ln(10.0) + log_extra() + 3.0 * det::ln(edge)
                    }
                },
                e2 * e2 * (1.0 + 4.0 * qi),
                || 4.0 * det::ln(ei) + det::ln(1.0 + 4.0 * qi),
            )
        }
    }
}

fn rates(
    kernel: &PairKernel,
    areas: [f64; 2],
    exchange_weight: f64,
    exchange_log_weight: impl FnOnce() -> f64,
    image_weight: f64,
    image_log_weight: impl Fn() -> f64,
) -> Rates {
    Rates {
        exchange: -det::expm1(-kernel.exponent_weight(areas, exchange_weight, exchange_log_weight)),
        loss: [
            kernel.exponent_weight([areas[1]; 2], image_weight, &image_log_weight),
            kernel.exponent_weight([areas[0]; 2], image_weight, image_log_weight),
        ],
    }
}

// Fixed geometry reuses the exponential across steps. Keep its exponent too:
// an underflowed factor can still leave representable large circulation.
#[derive(Clone, Copy)]
pub(super) struct CachedLoss {
    exponent: f64,
    factor: f64,
}

impl CachedLoss {
    pub(super) fn new(exponent: f64, particle: usize) -> Result<Self, VpmError> {
        if exponent.is_nan() || exponent < 0.0 {
            return Err(VpmError::NonFiniteResult {
                phase: "fs-vpm.wall-diffusion-loss",
                particle,
            });
        }
        Ok(Self {
            exponent,
            factor: det::exp(-exponent),
        })
    }

    pub(super) fn apply(self, value: f64, particle: usize) -> Result<f64, VpmError> {
        if self.exponent == 0.0 || value == 0.0 {
            return Ok(value);
        }
        let next = if self.factor.is_normal() {
            value * self.factor
        } else {
            det::exp(det::ln(value.abs()) - self.exponent).copysign(value)
        };
        if next.is_finite() {
            Ok(next)
        } else {
            Err(VpmError::NonFiniteResult {
                phase: "fs-vpm.wall-diffusion-loss",
                particle,
            })
        }
    }
}

fn decay(value: f64, exponent: f64, particle: usize) -> Result<f64, VpmError> {
    if exponent == 0.0 || value == 0.0 {
        return Ok(value);
    }
    CachedLoss::new(exponent, particle)?.apply(value, particle)
}
