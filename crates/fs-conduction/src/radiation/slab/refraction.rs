//! Incoherent Fresnel reflection around a common-index absorbing stack.
//!
//! Snell's law and polarized Fresnel power fractions follow PBRT 4e §9.3:
//! https://pbr-book.org/4ed/Reflection_Models/Specular_Reflection_and_Transmission
//! The thermal n² factor and trapped internal radiation are discussed by
//! Siegel and Spuckler (1992), https://ntrs.nasa.gov/citations/20150021009.
//! Here boundaries are specular; the cited paper's diffuse closure is NOT used.

use super::{AbsorbingLayer, ConductionError, MAX_SLAB_LAYERS, add, checkpoint};
use fs_math::decay::{decayed_fraction, remaining_fraction};

#[allow(clippy::too_many_arguments)]
pub(super) fn rule(
    layers: &[AbsorbingLayer],
    index: f64,
    lo: f64,
    hi: f64,
    trapped: bool,
    nodes: &[(f64, f64)],
    output: &mut [f64],
    correction: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
    evaluations: &mut usize,
) -> Result<(), ConductionError> {
    output.fill(0.0);
    correction.fill(0.0);
    let n = layers.len();
    let midpoint = 0.5 * (lo + hi);
    let half = 0.5 * (hi - lo);
    let inverse = index.recip();
    let critical_squared = (1.0 - inverse) * (1.0 + inverse);
    let total_depth: f64 = layers.iter().map(AbsorbingLayer::optical_depth).sum();
    for &(node, weight) in nodes {
        for sign in [-1.0, 1.0] {
            checkpoint(cancelled, *evaluations)?;
            let ordinate = midpoint + sign * half * node;
            let (mu, q, fresnel) = if trapped {
                // Both polarizations reflect completely outside the escape cone.
                (
                    ordinate,
                    2.0 * ordinate * half * weight * index * index,
                    [(1.0, 0.0); 2],
                )
            } else {
                let mu = (critical_squared + (ordinate * inverse).powi(2)).sqrt();
                (
                    mu,
                    2.0 * ordinate * half * weight,
                    [
                        power_fractions(ordinate * inverse, mu),
                        power_fractions(ordinate, mu * inverse),
                    ],
                )
            };
            let mut t = [0.0; MAX_SLAB_LAYERS];
            let mut a = [0.0; MAX_SLAB_LAYERS];
            let mut left = [1.0; MAX_SLAB_LAYERS];
            let mut right = [1.0; MAX_SLAB_LAYERS];
            let mut through = 1.0;
            for (i, layer) in layers.iter().enumerate() {
                left[i] = through;
                let depth = layer.optical_depth() / mu;
                t[i] = remaining_fraction(depth);
                a[i] = decayed_fraction(depth);
                through *= t[i];
            }
            let mut suffix = 1.0;
            for i in (0..n).rev() {
                right[i] = suffix;
                suffix *= t[i];
            }
            let absorbed_pass = decayed_fraction(total_depth / mu);
            for (r, escape) in fresnel {
                let q = 0.5 * q;
                // 1-r²T² = [(1-r)+r(1-T)](1+rT). This retains very weak
                // absorption in a nearly perfectly reflecting cavity.
                let denominator = (escape + r * absorbed_pass) * (1.0 + r * through);
                if denominator == 0.0 {
                    // A lossless trapped mode neither heats nor reaches a port.
                    continue;
                }
                let entered = escape / denominator;
                if !trapped {
                    let transmission = escape * through * entered;
                    let reflection = r + escape * r * through * through * entered;
                    add(&mut output[0], &mut correction[0], q * transmission);
                    let reverse = output.len() - 3;
                    add(
                        &mut output[reverse],
                        &mut correction[reverse],
                        q * transmission,
                    );
                    for last in output.len() - 2..output.len() {
                        add(&mut output[last], &mut correction[last], q * reflection);
                    }
                    for i in 0..n {
                        let front = 1 + 2 * i;
                        let back = front + 1;
                        add(
                            &mut output[front],
                            &mut correction[front],
                            q * a[i] * entered * (left[i] + r * through * right[i]),
                        );
                        add(
                            &mut output[back],
                            &mut correction[back],
                            q * a[i] * entered * (right[i] + r * through * left[i]),
                        );
                    }
                }
                for i in 0..n {
                    if a[i] == 0.0 {
                        continue;
                    }
                    let mut between = 1.0;
                    for j in i + 1..n {
                        // Direct i->j transport plus the two boundary-return
                        // fields. Every term is positive and symmetric in i,j.
                        let returns = left[i] * left[j]
                            + right[i] * right[j]
                            + r * through * (right[i] * left[j] + left[i] * right[j]);
                        let exchange =
                            a[i] * a[j] * between + a[i] * (a[j] / denominator) * r * returns;
                        let channel = 1 + 2 * n + i * n + j;
                        add(&mut output[channel], &mut correction[channel], q * exchange);
                        between *= t[j];
                    }
                }
            }
            *evaluations += 1;
        }
    }
    Ok(())
}

// R and 1-R are computed independently: a squared amplitude resolves small R,
// while 4ab/(a+b)² retains weak transmission near grazing incidence.
pub(super) fn power_fractions(a: f64, b: f64) -> (f64, f64) {
    let ratio = a.min(b) / a.max(b);
    let sum = 1.0 + ratio;
    (((1.0 - ratio) / sum).powi(2), 4.0 * ratio / (sum * sum))
}
