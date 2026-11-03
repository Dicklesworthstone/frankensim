//! Exact spatial derivatives of the same finite-core straight-segment
//! Biot–Savart integral used by the connected wake. Differentiating the
//! rationalized exterior expression retains thin, distant filament fields.

/// Velocity and its spatial Jacobian for unit circulation [m²/s].
/// Multiply both fields by the signed circulation of the source segment.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SegmentKinematics {
    /// Induced velocity per unit circulation [1/m].
    pub velocity: [f64; 3],
    /// `gradient[i][j] = ∂velocity[i]/∂position[j]` [1/m²].
    pub gradient: [[f64; 3]; 3],
}

/// Analytically integrate and differentiate a Rosenhead–Moore segment a→b.
///
/// This uses the same length-squared core as
/// [`super::segment_velocity_with_core`], with no finite-difference scale or
/// per-call allocation. At a positive-core filament axis velocity vanishes,
/// but its Jacobian and vorticity do not. A zero-length segment contributes
/// zero. `None` denotes non-finite/negative inputs, unrepresentable results,
/// or a probe on an unregularized segment, where the derivative is singular.
/// Probes on the extended axis outside a zero-core segment are regular.
///
/// With tangent `t`, transverse vector `c = t × (p-a)`, and line integral
/// `F = ∫ ds / (|p-a-st|² + core²)^(3/2)`, the derivative is
/// `∂j u = ((t × e_j) F + c ∂j F)/(4π)`. Thus curl is obtained from the
/// actual regularized velocity, not from a separately prescribed core profile.
#[must_use]
pub fn segment_kinematics_with_core(
    p: [f64; 3],
    a: [f64; 3],
    b: [f64; 3],
    core2_m2: f64,
) -> Option<SegmentKinematics> {
    if !core2_m2.is_finite()
        || core2_m2 < 0.0
        || !p
            .iter()
            .chain(a.iter())
            .chain(b.iter())
            .all(|x| x.is_finite())
    {
        return None;
    }
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let length = (ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2]).sqrt();
    if !length.is_finite() {
        return None;
    }
    if length == 0.0 {
        return ab
            .iter()
            .all(|x| *x == 0.0)
            .then_some(SegmentKinematics::default());
    }
    let t = ab.map(|x| x / length);
    let r1 = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let r2 = [p[0] - b[0], p[1] - b[1], p[2] - b[2]];
    let c = cross(t, r1);
    let n = cross(c, t);
    let transverse2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2] + core2_m2;
    let z1 = t[0] * r1[0] + t[1] * r1[1] + t[2] * r1[2];
    let z2 = t[0] * r2[0] + t[1] * r2[1] + t[2] * r2[2];
    let l1 = (z1 * z1 + transverse2).sqrt();
    let l2 = (z2 * z2 + transverse2).sqrt();
    let exterior = (z1 > 0.0 && z2 > 0.0) || (z1 < 0.0 && z2 < 0.0);
    if transverse2 == 0.0 && !exterior {
        return None;
    }
    let f = if exterior {
        length / (l1 * l2) * ((z1 + z2) / (z1 * l2 + z2 * l1))
    } else {
        (z1 / l1 - z2 / l2) / transverse2
    };
    let mut df = [0.0; 3];
    if exterior {
        // Logarithmic differentiation avoids subtracting nearly equal
        // endpoint integrals outside the segment, including its extended axis.
        let numerator = z1 + z2;
        let denominator = z1 * l2 + z2 * l1;
        for j in 0..3 {
            let dl1 = (z1 * t[j] + n[j]) / l1;
            let dl2 = (z2 * t[j] + n[j]) / l2;
            df[j] = f
                * (2.0 * t[j] / numerator
                    - (t[j] * (l1 + l2) + z1 * dl2 + z2 * dl1) / denominator
                    - dl1 / l1
                    - dl2 / l2);
        }
    } else {
        let inv_l1_3 = 1.0 / (l1 * l1 * l1);
        let inv_l2_3 = 1.0 / (l2 * l2 * l2);
        for j in 0..3 {
            df[j] = t[j] * (inv_l1_3 - inv_l2_3)
                - n[j] / transverse2 * (z1 * inv_l1_3 - z2 * inv_l2_3 + 2.0 * f);
        }
    }
    let inv_four_pi = 1.0 / (4.0 * core::f64::consts::PI);
    let mut out = SegmentKinematics {
        velocity: c.map(|x| x * (f * inv_four_pi)),
        gradient: [[0.0; 3]; 3],
    };
    for j in 0..3 {
        let mut e = [0.0; 3];
        e[j] = 1.0;
        let dc = cross(t, e);
        for i in 0..3 {
            out.gradient[i][j] = (dc[i] * f + c[i] * df[j]) * inv_four_pi;
        }
    }
    out.velocity
        .iter()
        .chain(out.gradient.iter().flatten())
        .all(|x| x.is_finite())
        .then_some(out)
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
