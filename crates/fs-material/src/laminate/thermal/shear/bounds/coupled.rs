//! Analytic box enclosure of the same equilibrated cylindrical complementary
//! energy as `laminate::shear::equilibrated`. Original source interval jets feed
//! the existing absolute derivative-row arithmetic; no temperatures are sampled.
//!
//! All scales below are fixed over the box and cancel algebraically. The small
//! Cholesky solves and signed probe inverse must remain regular over the ENTIRE
//! box. The free-face identity A N = B is algebraic; the native floating residual
//! check still runs at each resolved state, independently of this capacity proof.

use super::{
    AffineLaminaThermoelastic, Envelope, Interval, Jet, LaminateShearThermoelasticError, check,
    invalid, positive, scale_jet, square, upper_add, upper_mul, zero_jet,
};

struct Ply {
    h: Interval,
    q: [Jet<Interval>; 9],
    g: [Jet<Interval>; 4],
}

pub(super) fn hessian_rows(
    laws: &[AffineLaminaThermoelastic],
    angle: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, LaminateShearThermoelasticError> {
    let n = laws.len();
    let mut plies = Vec::with_capacity(n);
    let (mut h, mut q_scale, mut g_scale) = (Interval::from(0.0), 0.0_f64, 0.0_f64);
    for law in laws {
        check(cancelled)?;
        let (q, g) = law.stiffness_derivative_bounds();
        let input = law.source().input();
        let (s, c) = input.orientation_rad.sin_cos();
        let (s, c) = (Interval::from(s), Interval::from(c));
        let (ss, cc, cs) = (square(s), square(c), s * c);
        let t = [
            [cc, ss, cs],
            [ss, cc, -cs],
            [-cs * 2.0.into(), cs * 2.0.into(), cc - ss],
        ];
        let mut rotated = [zero_jet(); 9];
        for i in 0..3 {
            for j in i..3 {
                let mut entry = zero_jet();
                for a in 0..3 {
                    for b in 0..3 {
                        entry = entry + scale_jet(q[3 * a + b], t[a][i] * t[b][j]);
                    }
                }
                rotated[3 * i + j] = entry;
                rotated[3 * j + i] = entry;
            }
        }
        let cross = scale_jet(g[0] - g[1], cs);
        let g = [
            scale_jet(g[0], cc) + scale_jet(g[1], ss),
            cross,
            cross,
            scale_jet(g[0], ss) + scale_jet(g[1], cc),
        ];
        for coefficient in rotated.iter().chain(&g) {
            if [coefficient.v, coefficient.d, coefficient.dd]
                .iter()
                .any(|v| !v.finite())
            {
                return Err(invalid(
                    "whole box cannot bound the full coupled source tensors",
                ));
            }
        }
        q_scale = rotated
            .iter()
            .map(|x| x.v.abs_max())
            .fold(q_scale, f64::max);
        g_scale = g.iter().map(|x| x.v.abs_max()).fold(g_scale, f64::max);
        let thickness = Interval::from(input.thickness_m);
        positive(thickness)?;
        h = h + thickness;
        plies.push(Ply {
            h: thickness,
            q: rotated,
            g,
        });
    }
    positive(h)?;
    positive(q_scale.into())?;
    positive(g_scale.into())?;
    let q_inverse = Interval::from(1.0) / q_scale.into();
    let g_inverse = Interval::from(1.0) / g_scale.into();
    for ply in &mut plies {
        ply.q = ply.q.map(|x| scale_jet(x, q_inverse));
        ply.g = ply.g.map(|x| scale_jet(x, g_inverse));
    }
    let zero = || Envelope::constant(0.0, n);
    let mut a: [Envelope; 9] = core::array::from_fn(|_| zero());
    let (mut b, mut d) = (a.clone(), a.clone());
    let mut bottom = Interval::from(-0.5);
    for (i, ply) in plies.iter().enumerate() {
        check(cancelled)?;
        let dz = ply.h / h;
        let top = bottom + dz;
        let moments = [
            dz,
            dz * (bottom + top) / 2.0.into(),
            dz * (square(bottom) + bottom * top + square(top)) / 3.0.into(),
        ];
        for k in 0..9 {
            let q = Envelope::source(ply.q[k], i, n);
            a[k] = a[k].add(&q.scale(moments[0]));
            b[k] = b[k].add(&q.scale(moments[1]));
            d[k] = d[k].add(&q.scale(moments[2]));
        }
        bottom = top;
    }
    let neutral = neutral_matrix(a, &b)?;
    for entry in &neutral {
        entry.admit()?;
    }
    let relaxed = core::array::from_fn(|ij| {
        d[ij].sub(&(0..3).fold(zero(), |sum, k| {
            sum.add(&b[3 * (ij / 3) + k].mul(&neutral[3 * k + ij % 3]))
        }))
    });
    let (sin, cos) = angle.sin_cos();
    let axes = [[cos, sin], [-sin, cos]];
    let inverse = probe_inverse(divergence(&relaxed, axes))?;
    for entry in &inverse {
        entry.admit()?;
    }

    let mut primitive: [Envelope; 9] = core::array::from_fn(|_| zero());
    let mut compliance: [Envelope; 4] = core::array::from_fn(|_| zero());
    let mut raw = [zero(), zero()];
    bottom = Interval::from(-0.5);
    let root = Interval::from(3.0) / 5.0.into();
    let root = Interval::range(root.lo.sqrt(), root.hi.sqrt());
    let points = [
        (Interval::from(1.0) - root) / 2.0.into(),
        Interval::from(0.5),
        (Interval::from(1.0) + root) / 2.0.into(),
    ];
    let weights = [
        Interval::from(5.0) / 18.0.into(),
        Interval::from(4.0) / 9.0.into(),
        Interval::from(5.0) / 18.0.into(),
    ];
    for (i, ply) in plies.iter().enumerate() {
        check(cancelled)?;
        let dz = ply.h / h;
        let q = ply.q.map(|x| Envelope::source(x, i, n));
        let g = ply.g.map(|x| Envelope::source(x, i, n));
        raw[0] = raw[0].add(&g[0].scale(dz));
        raw[1] = raw[1].add(&g[3].scale(dz));
        let l00 = g[0].sqrt()?;
        let l00_inverse = l00.reciprocal()?;
        let l10 = g[2].mul(&l00_inverse);
        let l11_inverse = g[3].sub(&l10.square()).sqrt()?.reciprocal()?;
        // Once the actual membrane/probe/shear inverses have been enclosed,
        // one homogeneous layer has exactly H=(5h/6)G, independent of Q.
        // In particular, affine source G has exactly zero thermal curvature.
        if n == 1 {
            let factor = Interval::from(5.0) / 6.0.into() * h * g_scale.into();
            let stiffness = [g[0].scale(factor), g[1].scale(factor), g[3].scale(factor)];
            for entry in &stiffness {
                entry.admit()?;
            }
            positive(stiffness[0].value)?;
            positive(stiffness[2].value)?;
            return Ok(vec![upper_add(
                upper_add(stiffness[0].rows[0], stiffness[2].rows[0]),
                upper_mul(2.0, stiffness[1].rows[0]),
            )]);
        }
        let qn: [Envelope; 9] = core::array::from_fn(|ij| {
            (0..3).fold(zero(), |sum, k| {
                sum.add(&q[3 * (ij / 3) + k].mul(&neutral[3 * k + ij % 3]))
            })
        });
        let integral = |delta: Interval| -> [Envelope; 9] {
            core::array::from_fn(|k| {
                primitive[k].add(
                    &q[k]
                        .scale(bottom + delta / 2.0.into())
                        .sub(&qn[k])
                        .scale(delta),
                )
            })
        };
        // Quadratic equilibrium primitives give quartic energy. These positive
        // weights integrate it exactly; interval nodes enclose algebraic roots.
        for (point, weight) in points.into_iter().zip(weights) {
            let traction = divergence(&integral(dz * point), axes).map(|x| x.scale((-1.0).into()));
            let mut whitened: [Envelope; 4] = core::array::from_fn(|_| zero());
            for j in 0..2 {
                let tx = traction[0]
                    .mul(&inverse[j])
                    .add(&traction[1].mul(&inverse[2 + j]));
                let ty = traction[2]
                    .mul(&inverse[j])
                    .add(&traction[3].mul(&inverse[2 + j]));
                whitened[j] = tx.mul(&l00_inverse);
                whitened[2 + j] = ty.sub(&l10.mul(&whitened[j])).mul(&l11_inverse);
            }
            for j in 0..2 {
                compliance[3 * j] = compliance[3 * j].add(
                    &whitened[j]
                        .square()
                        .add(&whitened[2 + j].square())
                        .scale(dz * weight),
                );
            }
            compliance[1] = compliance[1].add(
                &whitened[0]
                    .mul(&whitened[1])
                    .add(&whitened[2].mul(&whitened[3]))
                    .scale(dz * weight),
            );
        }
        primitive = integral(dz);
        bottom = bottom + dz;
    }
    // Normalize before the final determinant too; both scales are constant in
    // temperature and disappear from H=g_scale*h*S^-1.
    let scale = compliance
        .iter()
        .map(|x| x.value.abs_max())
        .fold(0.0, f64::max);
    positive(scale.into())?;
    let c = compliance.map(|x| x.scale(Interval::from(1.0) / scale.into()));
    positive(c[0].value)?;
    positive(c[3].value)?;
    let det_inverse = c[0].mul(&c[3]).sub(&c[1].square()).reciprocal()?;
    let factor = Interval::from(g_scale) * h / scale.into();
    let stiffness = [
        c[3].mul(&det_inverse).scale(factor),
        c[1].mul(&det_inverse).scale(-factor),
        c[0].mul(&det_inverse).scale(factor),
    ];
    for entry in &stiffness {
        entry.admit()?;
    }
    for (diagonal, raw) in [(&stiffness[0], &raw[0]), (&stiffness[2], &raw[1])] {
        positive(diagonal.value)?;
        positive(diagonal.value / (raw.value * g_scale.into() * h))?;
    }
    check(cancelled)?;
    Ok((0..n)
        .map(|i| {
            upper_add(
                upper_add(stiffness[0].rows[i], stiffness[2].rows[i]),
                upper_mul(2.0, stiffness[1].rows[i]),
            )
        })
        .collect())
}

fn divergence(matrix: &[Envelope; 9], axes: [[f64; 2]; 2]) -> [Envelope; 4] {
    let n = matrix[0].rows.len();
    let mut out = core::array::from_fn(|_| Envelope::constant(0.0, n));
    for (column, axis) in axes.iter().enumerate() {
        let (x, y) = (Interval::from(axis[0]), Interval::from(axis[1]));
        let curvature = [square(x), square(y), Interval::from(2.0) * x * y];
        let stress: [Envelope; 3] = core::array::from_fn(|i| {
            (0..3).fold(Envelope::constant(0.0, n), |sum, j| {
                sum.add(&matrix[3 * i + j].scale(curvature[j]))
            })
        });
        out[column] = stress[0].scale(x).add(&stress[2].scale(y));
        out[2 + column] = stress[2].scale(x).add(&stress[1].scale(y));
    }
    out
}

fn probe_inverse(v: [Envelope; 4]) -> Result<[Envelope; 4], LaminateShearThermoelasticError> {
    let scale = v.iter().map(|x| x.value.abs_max()).fold(0.0, f64::max);
    positive(scale.into())?;
    let v = v.map(|x| x.scale(Interval::from(1.0) / scale.into()));
    let det = v[0].mul(&v[3]).sub(&v[1].mul(&v[2]));
    // A physical probe basis may have either orientation. Prove one sign over
    // the box, with at least the native normalized rank margin at every state.
    let sign = if det.value.lo > 128.0 * f64::EPSILON {
        1.0
    } else if det.value.hi < -128.0 * f64::EPSILON {
        -1.0
    } else {
        return Err(invalid(
            "whole temperature box cannot separate coupled bending-probe rank",
        ));
    };
    let inverse = det
        .scale(sign.into())
        .reciprocal()?
        .scale(Interval::from(sign) / scale.into());
    Ok([
        v[3].mul(&inverse),
        v[1].mul(&inverse).scale((-1.0).into()),
        v[2].mul(&inverse).scale((-1.0).into()),
        v[0].mul(&inverse),
    ])
}

// Fixed 3x3 Cholesky plus three solves, with every pivot enclosed above zero.
// This propagates the complete A^-1 B dependence, including mixed temperatures.
fn neutral_matrix(
    mut lower: [Envelope; 9],
    b: &[Envelope; 9],
) -> Result<[Envelope; 9], LaminateShearThermoelasticError> {
    for j in 0..3 {
        let mut diagonal = lower[3 * j + j].clone();
        for k in 0..j {
            diagonal = diagonal.sub(&lower[3 * j + k].square());
        }
        let pivot = diagonal.sqrt()?;
        let inverse = pivot.reciprocal()?;
        lower[3 * j + j] = pivot;
        for i in j + 1..3 {
            let mut value = lower[3 * i + j].clone();
            for k in 0..j {
                value = value.sub(&lower[3 * i + k].mul(&lower[3 * j + k]));
            }
            lower[3 * i + j] = value.mul(&inverse);
        }
    }
    let mut neutral = (*b).clone();
    for j in 0..3 {
        let mut rhs = [b[j].clone(), b[3 + j].clone(), b[6 + j].clone()];
        for i in 0..3 {
            let mut value = rhs[i].clone();
            for k in 0..i {
                value = value.sub(&lower[3 * i + k].mul(&rhs[k]));
            }
            rhs[i] = value.mul(&lower[3 * i + i].reciprocal()?);
        }
        for i in (0..3).rev() {
            let mut value = rhs[i].clone();
            for k in i + 1..3 {
                value = value.sub(&lower[3 * k + i].mul(&rhs[k]));
            }
            rhs[i] = value.mul(&lower[3 * i + i].reciprocal()?);
        }
        for i in 0..3 {
            neutral[3 * i + j] = rhs[i].clone();
        }
    }
    Ok(neutral)
}
