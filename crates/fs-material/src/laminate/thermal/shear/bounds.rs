//! Whole-temperature-box admission for the nonlocal shear entropy operator.
//!
//! Interval values and analytic absolute first/second derivative bounds follow
//! the same Whitney potential. A Hessian ROW SUM is propagated, so preparation
//! needs O(n²) scalar work and O(n) storage, without an n-by-n Hessian. Positive
//! square integrals keep the value enclosure useful without discarding any
//! derivative term. A failed enclosure refuses; it never becomes a pointwise
//! capacity test or an extrapolated material response.

use super::{LaminateShearThermoelasticError, ResolvedLaminateShearThermoelastic, check, invalid};
use crate::laminate::LaminateShearModel;
use crate::laminate::thermal::source::ResolvedLaminateThermal;
use crate::laminate::thermal::thermoelastic::{AffineLaminaThermoelastic, Interval, Jet};

/// Conservative fixed-total-shear entropy Hessian bounds over the Cartesian
/// product of all original temperature intervals and principal strain limits.
///
/// The i-th value bounds sum_j |dS_i/dT_j| PER reference area. Multiply by the
/// physical area, then add the original caloric and in-plane entropy derivatives
/// before admitting a positive thermal operator (for example, by Gershgorin).
/// This bounds the existing directional Whitney approximation, not a new 3-D
/// warping model, and does not itself certify a coupled thermal time step.
#[derive(Clone, Debug, PartialEq)]
pub struct LaminateShearCapacityBound {
    entropy_rows: Vec<f64>,
    section_shear_limit: f64,
}

impl LaminateShearCapacityBound {
    /// All laws must retain the exact source, reference geometry and order of
    /// `current`. The current state is checked against the native smooth branch;
    /// interval poles, nonfinite bounds and an unresolved positive shear-energy
    /// denominator anywhere in the box refuse before publishing a bound.
    /// Coupled equilibrated closures refuse: their point derivatives do not
    /// supply the Whitney interval enclosure used by this capacity proof.
    pub fn try_new(
        laws: &[AffineLaminaThermoelastic],
        current: &ResolvedLaminateThermal,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LaminateShearThermoelasticError> {
        check(&mut cancelled)?;
        // Coupled point derivatives do not certify this distinct Whitney
        // interval expression over a temperature box.
        if current.prepared().section().shear_model() != LaminateShearModel::DirectionalWhitney {
            return Err(invalid(
                "whole-temperature-box shear capacity bounds require the directional Whitney closure",
            ));
        }
        // Reuse source-identity, current-temperature and native-branch admission.
        ResolvedLaminateShearThermoelastic::try_new(laws, current, &mut cancelled)?;
        let n = laws.len();
        let geometry = current.prepared().section().layer_integrals();
        let mut layers = Vec::with_capacity(n);
        let mut limit = f64::INFINITY;
        for (law, geometry) in laws.iter().zip(geometry) {
            check(&mut cancelled)?;
            let (q, g) = law.stiffness_derivative_bounds();
            let (sin, cos) = law.source().input().orientation_rad.sin_cos();
            let (s, c) = (Interval::from(sin), Interval::from(cos));
            let (ss, cc, cs) = (square(s), square(c), c * s);
            let t = [
                [cc, ss],
                [ss, cc],
                [Interval::from(-2.0) * cs, Interval::from(2.0) * cs],
            ];
            let rotated_q = core::array::from_fn(|axis| {
                let mut sum = zero_jet();
                for a in 0..3 {
                    for b in 0..3 {
                        sum = sum + scale_jet(q[3 * a + b], t[a][axis] * t[b][axis]);
                    }
                }
                sum
            });
            let rotated_g = [
                scale_jet(g[0], cc) + scale_jet(g[1], ss),
                scale_jet(g[0] - g[1], cs),
                scale_jet(g[0], ss) + scale_jet(g[1], cc),
            ];
            for coefficient in rotated_q.iter().chain(rotated_g.iter()) {
                if [coefficient.v, coefficient.d, coefficient.dd]
                    .iter()
                    .any(|v| !v.finite())
                {
                    return Err(invalid(
                        "whole temperature box cannot bound the rotated source coefficients",
                    ));
                }
            }
            for value in [
                rotated_q[0].v,
                rotated_q[1].v,
                rotated_g[0].v,
                rotated_g[2].v,
            ] {
                positive(value)?;
            }
            // Invert the actual floating rotation coefficients: their squared
            // norm is not silently assumed to be exactly one.
            let rotation_norm = cc + ss;
            positive(rotation_norm)?;
            let envelope = Interval::from(law.maximum_absolute_engineering_strain())
                * (Interval::from(cos.abs()) + Interval::from(sin.abs()))
                / rotation_norm;
            if !envelope.finite() || envelope.hi <= 0.0 {
                return Err(invalid(
                    "whole-stack section shear envelope is unrepresentable",
                ));
            }
            limit = limit.min(envelope.hi);
            layers.push(BoundLayer {
                bottom: geometry.z_bottom_m,
                top: geometry.z_top_m,
                q: rotated_q,
                g: rotated_g,
            });
        }
        let hessian_rows = if n == 1 {
            // For ONE homogeneous layer the native correction is exactly 5/6,
            // including anisotropic off-axis G. Eliminate Q algebraically;
            // independently affine G then has exactly zero thermal curvature.
            let h = Interval::from(layers[0].top) - Interval::from(layers[0].bottom);
            positive(h)?;
            admit_homogeneous_native_scale(&layers[0], h)?;
            let factor = Interval::from(5.0) / Interval::from(6.0) * h;
            for g in &layers[0].g {
                if !(factor * g.v).finite() {
                    return Err(invalid(
                        "whole-box homogeneous shear stiffness is unrepresentable",
                    ));
                }
            }
            let cross = (layers[0].g[1].dd * factor).abs_max();
            vec![upper_add(
                upper_add(
                    (layers[0].g[0].dd * factor).abs_max(),
                    (layers[0].g[2].dd * factor).abs_max(),
                ),
                upper_mul(2.0, cross),
            )]
        } else {
            let x = direction(&layers, 0, 0, &mut cancelled)?;
            let y = direction(&layers, 1, 2, &mut cancelled)?;
            let mut raw_cross = Envelope::constant(0.0, n);
            for (i, layer) in layers.iter().enumerate() {
                check(&mut cancelled)?;
                raw_cross =
                    raw_cross.add(&Envelope::source(layer.g[1], i, n).scale(layer.thickness()?));
            }
            let cross = raw_cross.mul(&x.1.sqrt()?).mul(&y.1.sqrt()?);
            cross.admit()?;
            (0..n)
                .map(|i| {
                    upper_add(
                        upper_add(x.0.rows[i], y.0.rows[i]),
                        upper_mul(2.0, cross.rows[i]),
                    )
                })
                .collect()
        };
        let multiplier = upper_mul(0.5, upper_mul(limit, limit));
        let mut entropy_rows = Vec::with_capacity(n);
        for row in hessian_rows {
            check(&mut cancelled)?;
            // The original source affinity and equivalent native matrix
            // arithmetic have the same representation allowance as the plane
            // potential, in addition to all outward operations above.
            let value = upper_mul(upper_mul(multiplier, row), 1.0 + 4096.0 * f64::EPSILON);
            if !value.is_finite() || value < 0.0 {
                return Err(invalid(
                    "whole-stack entropy Hessian bound is unrepresentable",
                ));
            }
            entropy_rows.push(value);
        }
        check(&mut cancelled)?;
        Ok(Self {
            entropy_rows,
            section_shear_limit: limit,
        })
    }

    /// Bounds sum_j |dS_i/dT_j| / reference_area [J/(m² K²)].
    pub fn entropy_hessian_absolute_row_sums_j_m2_k2(&self) -> &[f64] {
        &self.entropy_rows
    }

    /// Conservative section-frame component envelope implied by ALL original
    /// principal shear envelopes. It is an enclosure, not permission to accept
    /// a state without checking those original rotated principal components.
    pub const fn section_shear_component_enclosure(&self) -> f64 {
        self.section_shear_limit
    }
}

struct BoundLayer {
    bottom: f64,
    top: f64,
    q: [Jet<Interval>; 2],
    g: [Jet<Interval>; 3],
}
impl BoundLayer {
    fn thickness(&self) -> Result<Interval, LaminateShearThermoelasticError> {
        let h = Interval::from(self.top) - Interval::from(self.bottom);
        positive(h)?;
        Ok(h)
    }
}
fn zero_jet() -> Jet<Interval> {
    Jet {
        v: Interval::from(0.0),
        d: Interval::from(0.0),
        dd: Interval::from(0.0),
    }
}
fn scale_jet(a: Jet<Interval>, scale: Interval) -> Jet<Interval> {
    Jet {
        v: a.v * scale,
        d: a.d * scale,
        dd: a.dd * scale,
    }
}

// Algebraic simplification must not hide a native overflow/underflow fallback
// elsewhere in the interval. For one layer, D=Qh³/12 and integral(r²/G)=Q²h⁵/(120G).
// Check the ORIGINAL squared numerator and product denominator before using
// the cancelled 5/6 result. Intermediate polynomial scales are checked too.
fn admit_homogeneous_native_scale(
    layer: &BoundLayer,
    h: Interval,
) -> Result<(), LaminateShearThermoelasticError> {
    let h2 = square(h);
    let h3 = h2 * h;
    let h4 = square(h2);
    let h5 = h4 * h;
    for (q, g) in [(layer.q[0].v, layer.g[0].v), (layer.q[1].v, layer.g[2].v)] {
        let q2 = square(q);
        for scale in [q * h, q * h2, q2 * h4, q2 * h5] {
            positive(scale)?;
        }
        let bending = q * h3 / Interval::from(12.0);
        let energy = q2 * h5 / (Interval::from(120.0) * g);
        let numerator = square(bending);
        let denominator = (g * h) * energy;
        positive(energy)?;
        positive(numerator)?;
        positive(denominator)?;
        positive(numerator / denominator)?;
    }
    Ok(())
}

fn positive(v: Interval) -> Result<(), LaminateShearThermoelasticError> {
    if v.finite() && v.lo > 0.0 {
        Ok(())
    } else {
        Err(invalid(
            "whole temperature box cannot enclose a positive smooth Whitney denominator",
        ))
    }
}
fn square(v: Interval) -> Interval {
    if !v.finite() {
        return Interval::range(f64::NEG_INFINITY, f64::INFINITY);
    }
    let upper = v.abs_max() * v.abs_max();
    let lower = if v.lo <= 0.0 && v.hi >= 0.0 {
        0.0
    } else {
        v.lo.abs().min(v.hi.abs()).powi(2)
    };
    let mut result = Interval::range(lower, upper);
    // Zero is an exact lower bound of a square, including underflow.
    if result.lo < 0.0 {
        result.lo = 0.0;
    }
    result
}
fn upper_add(a: f64, b: f64) -> f64 {
    if a == 0.0 {
        b
    } else if b == 0.0 {
        a
    } else {
        (a + b).next_up()
    }
}
fn upper_mul(a: f64, b: f64) -> f64 {
    if a == 0.0 || b == 0.0 {
        0.0
    } else {
        (a * b).next_up()
    }
}

/// Analytic norm enclosure, not a sampled Hessian or a derivative approximation.
/// `first[i] >= |f_i|`, `rows[i] >= sum_j |f_ij|` everywhere in the box.
#[derive(Clone)]
struct Envelope {
    value: Interval,
    first: Vec<f64>,
    rows: Vec<f64>,
}
impl Envelope {
    fn admit(&self) -> Result<(), LaminateShearThermoelasticError> {
        if self.value.finite()
            && self
                .first
                .iter()
                .chain(&self.rows)
                .all(|v| v.is_finite() && *v >= 0.0)
        {
            Ok(())
        } else {
            Err(invalid(
                "whole-stack derivative enclosure is unrepresentable",
            ))
        }
    }
    fn constant(v: f64, n: usize) -> Self {
        Self::interval(Interval::from(v), n)
    }
    fn interval(value: Interval, n: usize) -> Self {
        Self {
            value,
            first: vec![0.0; n],
            rows: vec![0.0; n],
        }
    }
    fn source(v: Jet<Interval>, coordinate: usize, n: usize) -> Self {
        let mut out = Self::interval(v.v, n);
        out.first[coordinate] = v.d.abs_max();
        out.rows[coordinate] = v.dd.abs_max();
        out
    }
    fn add(&self, b: &Self) -> Self {
        Self {
            value: self.value + b.value,
            first: self
                .first
                .iter()
                .zip(&b.first)
                .map(|(&a, &b)| upper_add(a, b))
                .collect(),
            rows: self
                .rows
                .iter()
                .zip(&b.rows)
                .map(|(&a, &b)| upper_add(a, b))
                .collect(),
        }
    }
    fn sub(&self, b: &Self) -> Self {
        let mut out = self.add(b);
        out.value = self.value - b.value;
        out
    }
    fn scale(&self, scale: Interval) -> Self {
        let bound = scale.abs_max();
        Self {
            value: self.value * scale,
            first: self.first.iter().map(|&a| upper_mul(a, bound)).collect(),
            rows: self.rows.iter().map(|&a| upper_mul(a, bound)).collect(),
        }
    }
    fn mul(&self, b: &Self) -> Self {
        let (a_abs, b_abs) = (self.value.abs_max(), b.value.abs_max());
        let a_sum = self.first.iter().copied().fold(0.0, upper_add);
        let b_sum = b.first.iter().copied().fold(0.0, upper_add);
        Self {
            value: self.value * b.value,
            first: self
                .first
                .iter()
                .zip(&b.first)
                .map(|(&a, &b)| upper_add(upper_mul(a, b_abs), upper_mul(b, a_abs)))
                .collect(),
            rows: (0..self.rows.len())
                .map(|i| {
                    upper_add(
                        upper_add(upper_mul(self.rows[i], b_abs), upper_mul(b.rows[i], a_abs)),
                        upper_add(
                            upper_mul(self.first[i], b_sum),
                            upper_mul(b.first[i], a_sum),
                        ),
                    )
                })
                .collect(),
        }
    }
    fn square(&self) -> Self {
        let mut out = self.mul(self);
        out.value = square(self.value);
        out
    }
    fn unary(&self, value: Interval, first: f64, second: f64) -> Self {
        let sum = self.first.iter().copied().fold(0.0, upper_add);
        Self {
            value,
            first: self.first.iter().map(|&d| upper_mul(first, d)).collect(),
            rows: self
                .rows
                .iter()
                .zip(&self.first)
                .map(|(&h, &d)| {
                    upper_add(upper_mul(first, h), upper_mul(second, upper_mul(d, sum)))
                })
                .collect(),
        }
    }
    fn reciprocal(&self) -> Result<Self, LaminateShearThermoelasticError> {
        positive(self.value)?;
        let inverse = Interval::from(1.0) / self.value;
        let inverse_square = square(inverse);
        let second = Interval::from(2.0) * inverse_square * inverse;
        if !inverse.finite() || !second.finite() {
            return Err(invalid(
                "whole-stack reciprocal derivative bound is unrepresentable",
            ));
        }
        Ok(self.unary(inverse, inverse_square.hi, second.hi))
    }
    fn sqrt(&self) -> Result<Self, LaminateShearThermoelasticError> {
        positive(self.value)?;
        let root = Interval::range(self.value.lo.sqrt(), self.value.hi.sqrt());
        let first = Interval::from(0.5) / root;
        let second = Interval::from(0.25) / (self.value * root);
        if !first.finite() || !second.finite() {
            return Err(invalid(
                "whole-stack square-root derivative bound is unrepresentable",
            ));
        }
        Ok(self.unary(root, first.hi, second.hi))
    }
}

/// H_axis = D_axis² / integral(r_axis² / G_axis), plus its correction H/S.
/// Centered second moments and positive three-point Gauss square integrals
/// are algebraically identical to the native polynomial integrals. The rule
/// integrates degree four exactly; interval nodes/weights enclose their exact
/// algebraic values. No quadrature-error or sampled-temperature claim is used.
fn direction(
    layers: &[BoundLayer],
    q_index: usize,
    g_index: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Envelope, Envelope), LaminateShearThermoelasticError> {
    let n = layers.len();
    let mut a = Envelope::constant(0.0, n);
    let mut first = a.clone();
    let mut shear = a.clone();
    // Fixed origin shift reduces interval dependency on an arbitrary reference
    // surface while retaining the same physical z integrals.
    let origin = 0.5 * layers[0].bottom + 0.5 * layers[n - 1].top;
    for (i, layer) in layers.iter().enumerate() {
        check(cancelled)?;
        let h = layer.thickness()?;
        let center = (Interval::from(layer.bottom) + Interval::from(layer.top))
            * Interval::from(0.5)
            - Interval::from(origin);
        let qh = Envelope::source(layer.q[q_index], i, n).scale(h);
        a = a.add(&qh);
        first = first.add(&qh.scale(center));
        shear = shear.add(&Envelope::source(layer.g[g_index], i, n).scale(h));
    }
    let neutral = first.mul(&a.reciprocal()?);
    let mut bending = Envelope::constant(0.0, n);
    let mut energy = bending.clone();
    let mut running = bending.clone();
    let root = Interval::from(3.0) / Interval::from(5.0);
    let root = Interval::range(root.lo.sqrt(), root.hi.sqrt());
    let points = [
        (Interval::from(1.0) - root) * Interval::from(0.5),
        Interval::from(0.5),
        (Interval::from(1.0) + root) * Interval::from(0.5),
    ];
    let weights = [
        Interval::from(5.0) / Interval::from(18.0),
        Interval::from(4.0) / Interval::from(9.0),
        Interval::from(5.0) / Interval::from(18.0),
    ];
    for (i, layer) in layers.iter().enumerate() {
        check(cancelled)?;
        let h = layer.thickness()?;
        let q = Envelope::source(layer.q[q_index], i, n);
        let g = Envelope::source(layer.g[g_index], i, n);
        let bottom = Envelope::interval(Interval::from(layer.bottom) - Interval::from(origin), n)
            .sub(&neutral);
        let center = bottom.add(&Envelope::interval(h * Interval::from(0.5), n));
        let second = center
            .square()
            .add(&Envelope::interval(square(h) / Interval::from(12.0), n));
        bending = bending.add(&q.mul(&second).scale(h));
        let inverse_g = g.reciprocal()?;
        for (point, weight) in points.into_iter().zip(weights) {
            let offset = h * point;
            let integral = bottom
                .scale(offset)
                .add(&Envelope::interval(square(offset) * Interval::from(0.5), n));
            let stress_integral = running.add(&q.mul(&integral));
            energy = energy.add(&stress_integral.square().mul(&inverse_g).scale(h * weight));
        }
        running = running.add(&q.mul(&center).scale(h));
    }
    positive(bending.value)?;
    positive(energy.value)?;
    positive(shear.value)?;
    let bending_squared = bending.square();
    // The native correction must retain a finite positive denominator too;
    // algebraic cancellation does not excuse an unrepresentable native branch.
    positive(shear.value * energy.value)?;
    let h = bending_squared.mul(&energy.reciprocal()?);
    positive(h.value)?;
    let correction = h.mul(&shear.reciprocal()?);
    positive(correction.value)?;
    h.admit()?;
    correction.admit()?;
    Ok((h, correction))
}
