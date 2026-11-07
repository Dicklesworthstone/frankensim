//! Specular incoherent exchange between layers with distinct real indices.
//!
//! The conserved Snell coordinate is p = n sin(theta). In every propagating
//! layer, n² 2 mu dmu = -d(p²), so one common positive measure carries both
//! exterior illumination and internal thermal radiation without angular mixing.
//! Split at every critical p=n: an unavailable propagating channel is a perfect
//! reflector, never an invented evanescent transmission through a finite gap.
//!
//! Fresnel power scattering and incoherent propagation are the ray-intensity
//! limit of multilayer optics (Byrnes, https://arxiv.org/abs/1603.02720, §6).
//! Here the thermal source is local equilibrium radiance divided by n².
//! Prefix/suffix reflection sweeps eliminate all round trips in O(N); the
//! complete symmetric thermal exchange then costs O(N²), its output size.
//! There are no exponentially growing transfer matrices or iterative bounces.

use super::refraction::power_fractions;
use super::{AbsorbingLayer, ConductionError, MAX_SLAB_LAYERS, add, checkpoint};
use fs_math::decay::{decayed_fraction, remaining_fraction};
use fs_math::det;

#[allow(clippy::too_many_arguments)]
pub(super) fn rule(
    layers: &[AbsorbingLayer],
    band: [f64; 2],
    lo: f64,
    hi: f64,
    nodes: &[(f64, f64)],
    output: &mut [f64],
    correction: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
    evaluations: &mut usize,
) -> Result<(), ConductionError> {
    output.fill(0.0);
    correction.fill(0.0);
    let n = layers.len();
    let [lower, upper] = band;
    let width = (upper - lower) * (upper + lower);
    let midpoint = 0.5 * (lo + hi);
    let half = 0.5 * (hi - lo);
    let mut indices = [1.0; MAX_SLAB_LAYERS + 2];
    for (i, layer) in layers.iter().enumerate() {
        indices[i + 1] = layer.refractive_index();
    }
    for &(node, weight) in nodes {
        for sign in [-1.0, 1.0] {
            checkpoint(cancelled, *evaluations)?;
            let ordinate = midpoint + sign * half * node;
            let (sin, cos) = det::sin_cos(core::f64::consts::FRAC_PI_2 * ordinate);
            // p² = lower² + (upper²-lower²)sin²(theta). Its vanishing
            // endpoint Jacobian removes the critical-angle square-root cusp.
            let q = width * (core::f64::consts::PI * sin * cos * half * weight);
            let mut mu = [0.0; MAX_SLAB_LAYERS + 2];
            for i in 0..n + 2 {
                let index = indices[i];
                if index > lower {
                    // All active indices are >= upper. Form index²-p² as a
                    // positive sum, retaining the grazing cosine near p=index.
                    let gap = ((index - upper) / index) * ((index + upper) / index);
                    let range = ((upper - lower) / index) * ((upper + lower) / index);
                    mu[i] = (gap + range * cos * cos).sqrt();
                }
            }
            let mut t = [0.0; MAX_SLAB_LAYERS];
            let mut a = [0.0; MAX_SLAB_LAYERS];
            for (i, layer) in layers.iter().enumerate() {
                if indices[i + 1] > lower {
                    let depth = if layer.optical_depth() == 0.0 {
                        0.0
                    } else {
                        layer.optical_depth() / mu[i + 1]
                    };
                    t[i] = remaining_fraction(depth);
                    a[i] = decayed_fraction(depth);
                }
            }
            for polarization in 0..2 {
                let mut r = [1.0; MAX_SLAB_LAYERS + 1];
                let mut e = [0.0; MAX_SLAB_LAYERS + 1];
                for i in 0..=n {
                    if indices[i] > lower && indices[i + 1] > lower {
                        (r[i], e[i]) = if indices[i] == indices[i + 1] {
                            (0.0, 1.0)
                        } else if polarization == 0 {
                            power_fractions(indices[i] * mu[i], indices[i + 1] * mu[i + 1])
                        } else {
                            power_fractions(indices[i] * mu[i + 1], indices[i + 1] * mu[i])
                        };
                    }
                }
                accumulate(
                    &t[..n],
                    &a[..n],
                    &r[..=n],
                    &e[..=n],
                    lower == 0.0,
                    true,
                    0.5 * q,
                    output,
                    correction,
                );
            }
            *evaluations += 1;
        }
    }
    Ok(())
}

// Reflection seen through one interface followed by a layer and its suffix.
// Keep the deficit independently of R: a weakly absorbing closed cavity may
// have R rounded to one while its finite thermal exchange remains resolvable.
fn joined(r: f64, e: f64, t: f64, a: f64, back: f64, deficit: f64) -> (f64, f64, f64) {
    if e == 0.0 {
        return (1.0, 0.0, 0.0);
    }
    let loss = a * (1.0 + t) + t * t * deficit;
    let denominator = e + r * loss;
    let crossing = e / denominator;
    (r + e * (t * t * back * crossing), loss * crossing, crossing)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn accumulate(
    t: &[f64],
    a: &[f64],
    r: &[f64],
    e: &[f64],
    exterior: bool,
    thermal_exchange: bool,
    q: f64,
    output: &mut [f64],
    correction: &mut [f64],
) {
    let n = a.len();
    let mut left = [0.0; MAX_SLAB_LAYERS];
    let mut right = [0.0; MAX_SLAB_LAYERS];
    let mut left_deficit = [0.0; MAX_SLAB_LAYERS];
    let mut right_deficit = [0.0; MAX_SLAB_LAYERS];
    let mut cross_left = [0.0; MAX_SLAB_LAYERS];
    let mut cross_right = [0.0; MAX_SLAB_LAYERS];
    left[0] = r[0];
    left_deficit[0] = e[0];
    for i in 1..n {
        (left[i], left_deficit[i], cross_left[i]) = joined(
            r[i],
            e[i],
            t[i - 1],
            a[i - 1],
            left[i - 1],
            left_deficit[i - 1],
        );
    }
    right[n - 1] = r[n];
    right_deficit[n - 1] = e[n];
    for i in (0..n - 1).rev() {
        (right[i], right_deficit[i], cross_right[i]) = joined(
            r[i + 1],
            e[i + 1],
            t[i + 1],
            a[i + 1],
            right[i + 1],
            right_deficit[i + 1],
        );
    }
    if exterior {
        let (front, _, mut field) = joined(r[0], e[0], t[0], a[0], right[0], right_deficit[0]);
        let front_channel = output.len() - 2;
        add(
            &mut output[front_channel],
            &mut correction[front_channel],
            q * front,
        );
        for i in 0..n {
            let channel = 1 + 2 * i;
            let absorbed = a[i] * (1.0 + t[i] * right[i]) * field;
            add(&mut output[channel], &mut correction[channel], q * absorbed);
            field *= t[i];
            if i + 1 < n {
                field *= cross_right[i];
            }
        }
        add(&mut output[0], &mut correction[0], q * field * e[n]);
        let (back, _, mut field) = joined(
            r[n],
            e[n],
            t[n - 1],
            a[n - 1],
            left[n - 1],
            left_deficit[n - 1],
        );
        let back_channel = output.len() - 1;
        add(
            &mut output[back_channel],
            &mut correction[back_channel],
            q * back,
        );
        for i in (0..n).rev() {
            let channel = 2 + 2 * i;
            let absorbed = a[i] * (1.0 + t[i] * left[i]) * field;
            add(&mut output[channel], &mut correction[channel], q * absorbed);
            field *= t[i];
            if i > 0 {
                field *= cross_left[i];
            }
        }
    }
    if !thermal_exchange {
        return;
    }
    for i in 0..n {
        if a[i] == 0.0 {
            continue;
        }
        let denominator =
            a[i] * (1.0 + t[i]) + t[i] * t[i] * (left_deficit[i] + left[i] * right_deficit[i]);
        let emitted = (a[i] / denominator) * (1.0 + t[i] * left[i]);
        let mut field = emitted * cross_right[i];
        for j in i + 1..n {
            if field == 0.0 {
                break;
            }
            let channel = 1 + 2 * n + i * n + j;
            let exchange = field * a[j] * (1.0 + t[j] * right[j]);
            add(&mut output[channel], &mut correction[channel], q * exchange);
            field *= t[j] * cross_right[j];
        }
    }
}
