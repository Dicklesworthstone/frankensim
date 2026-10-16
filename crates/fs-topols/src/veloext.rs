//! Velocity extension: normal speeds are DEFINED on the interface
//! (energy densities on cut cells) and must reach the whole band
//! before advection. The extension solves ∇v·∇φ = 0 outward by
//! processing nodes in ascending |φ| (each node upwinds from
//! already-final smaller-|φ| neighbors) — the FIM-adjacent ordered
//! sweep, deterministic by construction.

#![allow(clippy::cast_possible_wrap)] // lattice indices ≤ 2^17 — i64 stencil arithmetic is exact

use crate::gridsdf::GridSdf;

/// Extend nodal values known on `seeded` nodes to the whole grid
/// along interface normals. `values` holds the seeds on entry and the
/// extended field on exit.
pub fn extend_velocity(phi: &GridSdf, values: &mut [f64], seeded: &[bool]) {
    let n = phi.n();
    let stride = n + 1;
    assert_eq!(values.len(), stride * stride, "nodal length");
    assert_eq!(seeded.len(), stride * stride, "seed mask length");
    // Order nodes by ascending |φ| (deterministic tie-break on index).
    let mut order: Vec<usize> = (0..stride * stride).collect();
    order.sort_by(|&a, &b| {
        phi.nodes()[a]
            .abs()
            .partial_cmp(&phi.nodes()[b].abs())
            .expect("finite phi")
            .then(a.cmp(&b))
    });
    let mut done: Vec<bool> = seeded.to_vec();
    for &k in &order {
        if done[k] {
            continue;
        }
        let (i, j) = (k % stride, k / stride);
        // Upwind weights: |φ| gradient direction, taking only
        // FINALIZED neighbors with smaller |φ|.
        let mut wsum = 0.0f64;
        let mut vsum = 0.0f64;
        let me = phi.nodes()[k].abs();
        for (di, dj) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)] {
            let (ni, nj) = (i as i64 + di, j as i64 + dj);
            if ni < 0 || nj < 0 || ni > n as i64 || nj > n as i64 {
                continue;
            }
            #[allow(clippy::cast_sign_loss)]
            let nk = ni as usize + nj as usize * stride;
            if !done[nk] {
                continue;
            }
            let nb = phi.nodes()[nk].abs();
            if nb <= me {
                let w = (me - nb).max(1e-12);
                wsum += w;
                vsum += w * values[nk];
            }
        }
        if wsum > 0.0 {
            values[k] = vsum / wsum;
        }
        done[k] = true;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extension_is_constant_along_normals_of_a_plane() {
        // φ = x − 0.5: normals point in x; a seed profile varying in y
        // must extend unchanged along x.
        let phi = GridSdf::from_fn(16, &|x, _| x - 0.5);
        let stride = 17usize;
        let mut vals = vec![0.0f64; stride * stride];
        let mut seeds = vec![false; stride * stride];
        for j in 0..stride {
            // Seed the two columns bracketing the interface.
            {
                let i = 8usize;
                let v = (0.3 * j as f64).sin();
                vals[i + j * stride] = v;
                seeds[i + j * stride] = true;
            }
        }
        extend_velocity(&phi, &mut vals, &seeds);
        for j in 0..stride {
            let want = (0.3 * j as f64).sin();
            for i in 0..stride {
                assert!(
                    (vals[i + j * stride] - want).abs() < 1e-9,
                    "normal-constant extension"
                );
            }
        }
    }
}
