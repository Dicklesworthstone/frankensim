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
    exterior_refractive_index: f64,
    band: [f64; 2],
    lo: f64,
    hi: f64,
    nodes: &[(f64, f64)],
    output: &mut [f64],
    correction: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
    evaluations: &mut usize,
    interval: Option<(usize, [f64; 2])>,
) -> Result<(), ConductionError> {
    output.fill(0.0);
    correction.fill(0.0);
    let n = layers.len();
    let [lower, upper] = band;
    let width = (upper - lower) * (upper + lower);
    let midpoint = 0.5 * (lo + hi);
    let half = 0.5 * (hi - lo);
    let mut indices = [exterior_refractive_index; MAX_SLAB_LAYERS + 2];
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
                mu[i] = band_direction_cosine(indices[i], band, cos);
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
                if let Some((target, [start, end])) = interval {
                    if indices[target + 1] > lower {
                        let material = &layers[target];
                        let path = |distance| {
                            super::beam::optical_path(
                                material.absorption_m_inv(),
                                distance,
                                mu[target + 1],
                            )
                        };
                        let loss = decayed_fraction(path(end - start));
                        let attenuation = [
                            remaining_fraction(path(start)),
                            remaining_fraction(path(material.thickness_m() - end)),
                        ];
                        accumulate_depth(
                            &t[..n],
                            &a[..n],
                            &r[..=n],
                            &e[..=n],
                            lower < exterior_refractive_index,
                            target,
                            loss,
                            attenuation,
                            0.5 * q,
                            0.5 * (q / exterior_refractive_index.powi(2)),
                            output,
                            correction,
                        );
                    }
                } else {
                    accumulate(
                        &t[..n],
                        &a[..n],
                        &r[..=n],
                        &e[..=n],
                        lower < exterior_refractive_index,
                        true,
                        0.5 * q,
                        0.5 * (q / exterior_refractive_index.powi(2)),
                        output,
                        correction,
                        None,
                    );
                }
            }
            *evaluations += 1;
        }
    }
    Ok(())
}

// Within a band split at EVERY refractive index, each active index is at
// least upper. Preserve the positive Snell discriminant near critical angles,
// where reconstructing it from a rounded exterior cosine loses accuracy.
pub(super) fn band_direction_cosine(index: f64, [lower, upper]: [f64; 2], cosine: f64) -> f64 {
    if index <= lower {
        0.0
    } else {
        let gap = ((index - upper) / index) * ((index + upper) / index);
        let range = ((upper - lower) / index) * ((upper + lower) / index);
        (gap + range * cosine * cosine).sqrt()
    }
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

struct Returns {
    left: [f64; MAX_SLAB_LAYERS],
    right: [f64; MAX_SLAB_LAYERS],
    left_deficit: [f64; MAX_SLAB_LAYERS],
    right_deficit: [f64; MAX_SLAB_LAYERS],
    cross_left: [f64; MAX_SLAB_LAYERS],
    cross_right: [f64; MAX_SLAB_LAYERS],
}

fn returns(t: &[f64], a: &[f64], r: &[f64], e: &[f64]) -> Returns {
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
    Returns {
        left,
        right,
        left_deficit,
        right_deficit,
        cross_left,
        cross_right,
    }
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
    exterior_q: f64,
    output: &mut [f64],
    correction: &mut [f64],
    mut beam_fields: Option<&mut [[[f64; 2]; 2]]>,
) {
    let n = a.len();
    let Returns {
        left,
        right,
        left_deficit,
        right_deficit,
        cross_left,
        cross_right,
    } = returns(t, a, r, e);
    if exterior {
        // Exterior fractions use 2 mu_ext dmu_ext; internal exchange uses
        // the absolute optical-mode measure d(p²), retaining n² radiance.
        let q = exterior_q;
        let (front, _, mut field) = joined(r[0], e[0], t[0], a[0], right[0], right_deficit[0]);
        let front_channel = output.len() - 2;
        add(
            &mut output[front_channel],
            &mut correction[front_channel],
            q * front,
        );
        for i in 0..n {
            if let Some(fields) = beam_fields.as_deref_mut() {
                // Direction [forward@front, backward@back], incident side.
                fields[i][0][0] += q * field;
                fields[i][1][0] += q * (t[i] * right[i] * field);
            }
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
            if let Some(fields) = beam_fields.as_deref_mut() {
                fields[i][1][1] += q * field;
                fields[i][0][1] += q * (t[i] * left[i] * field);
            }
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

// Couplings into an interval from the two exterior faces and every other
// original layer. Source-layer self-return is eliminated by detailed balance
// in the caller's difference form, not by subtracting large emission terms.
#[allow(clippy::too_many_arguments)]
fn accumulate_depth(
    t: &[f64],
    a: &[f64],
    r: &[f64],
    e: &[f64],
    exterior: bool,
    target: usize,
    loss: f64,
    attenuation: [f64; 2],
    q: f64,
    exterior_q: f64,
    output: &mut [f64],
    correction: &mut [f64],
) {
    if loss == 0.0 {
        return;
    }
    let n = a.len();
    let Returns {
        left,
        right,
        left_deficit,
        right_deficit,
        cross_left,
        cross_right,
    } = returns(t, a, r, e);
    let from_front = loss * (attenuation[0] + t[target] * right[target] * attenuation[1]);
    let from_back = loss * (attenuation[1] + t[target] * left[target] * attenuation[0]);
    if exterior {
        let q = exterior_q;
        let (_, _, mut field) = joined(r[0], e[0], t[0], a[0], right[0], right_deficit[0]);
        for i in 0..target {
            field *= t[i];
            field *= cross_right[i];
        }
        add(&mut output[0], &mut correction[0], q * from_front * field);
        let (_, _, mut field) = joined(
            r[n],
            e[n],
            t[n - 1],
            a[n - 1],
            left[n - 1],
            left_deficit[n - 1],
        );
        for i in (target + 1..n).rev() {
            field *= t[i];
            field *= cross_left[i];
        }
        add(&mut output[1], &mut correction[1], q * from_back * field);
    }
    let mut path = 1.0;
    for source in (0..target).rev() {
        path *= cross_right[source];
        if path == 0.0 {
            break;
        }
        if a[source] > 0.0 {
            let denominator = a[source] * (1.0 + t[source])
                + t[source]
                    * t[source]
                    * (left_deficit[source] + left[source] * right_deficit[source]);
            let emitted = (a[source] / denominator) * (1.0 + t[source] * left[source]);
            add(
                &mut output[2 + source],
                &mut correction[2 + source],
                q * emitted * path * from_front,
            );
        }
        path *= t[source];
    }
    let mut path = 1.0;
    for source in target + 1..n {
        path *= cross_left[source];
        if path == 0.0 {
            break;
        }
        if a[source] > 0.0 {
            let denominator = a[source] * (1.0 + t[source])
                + t[source]
                    * t[source]
                    * (left_deficit[source] + left[source] * right_deficit[source]);
            let emitted = (a[source] / denominator) * (1.0 + t[source] * right[source]);
            add(
                &mut output[2 + source],
                &mut correction[2 + source],
                q * emitted * path * from_back,
            );
        }
        path *= t[source];
    }
}
