//! Fixed-size compliance jets and a sufficient whole-interval SPD proof.
//! Engineering constants are original affine source paths, never fitted tensors.
use super::*;
use fs_math::det;

type Matrix3 = [[f64; 3]; 3];

/// (a+b*dt)/(1+d*dt). The constant cross product avoids subtracting
/// temperature-dependent products when differentiating a nearly constant ratio.
#[derive(Clone, Copy, Debug)]
struct AffineRatio {
    a: f64,
    b: f64,
    d: f64,
    cross: f64,
}
impl AffineRatio {
    fn new(a: f64, b: f64, d: f64) -> Self {
        Self {
            a,
            b,
            d,
            cross: (-a).mul_add(d, b),
        }
    }
    fn jet(self, dt: f64) -> Result<[f64; 3]> {
        let denominator = self.d.mul_add(dt, 1.0);
        let value = self.b.mul_add(dt, self.a) / denominator;
        let first = (self.cross / denominator) / denominator;
        let second = -2.0 * (self.d / denominator) * first;
        if !denominator.is_finite()
            || denominator <= 0.0
            || [value, first, second].iter().any(|v| !v.is_finite())
        {
            return Err(TensorSourceError::Invalid(
                "thermal compliance ratio is unrepresentable or leaves positive moduli".into(),
            ));
        }
        Ok([value, first, second])
    }
}

/// D=sqrt(E_reference), A=D S D. Normalized source ratios avoid pressure
/// products in inversion. Shear slots are the native Mandel values 2G.
#[derive(Clone, Copy, Debug)]
pub(super) struct ComplianceThermal {
    reference_temperature_k: f64,
    roots: [f64; 3],
    normal: [[AffineRatio; 3]; 3],
    shear: [AffineRatio; 3],
    isotropic: bool,
    rotation: Option<Tangent6>,
}
impl ComplianceThermal {
    pub(super) fn admit(
        source: &ResolvedTensorMaterial,
        slopes: &BTreeMap<String, f64>,
        reference_temperature_k: f64,
        interval: [f64; 2],
    ) -> Result<(Self, TensorTemperatureBounds)> {
        let isotropic = source.properties.elasticity == TensorElasticity::Isotropic;
        let value = |name: &str| {
            source
                .state
                .property(name)
                .expect("original elastic source")
                .value_si()
        };
        let young_names = if isotropic {
            ["young_modulus"; 3]
        } else {
            ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES
        };
        let poisson_names = if isotropic {
            ["poisson_ratio"; 3]
        } else {
            ORTHOTROPIC_POISSON_RATIO_PROPERTIES
        };
        let young = young_names.map(value);
        let roots = young.map(det::sqrt);
        let rates = young_names.map(|name| slopes[name] / value(name));
        if roots.iter().chain(rates.iter()).any(|v| !v.is_finite())
            || roots.iter().any(|v| *v <= 0.0)
        {
            return Err(TensorSourceError::Invalid(
                "thermal compliance scaling is unrepresentable".into(),
            ));
        }
        let zero = AffineRatio::new(0.0, 0.0, 0.0);
        let mut normal = [[zero; 3]; 3];
        for i in 0..3 {
            normal[i][i] = AffineRatio::new(1.0, 0.0, rates[i]);
        }
        // Native major ratios are (nu12,nu13,nu23), with denominator E1,E1,E2.
        for (k, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
            let name = poisson_names[k];
            let scale = roots[j] / roots[i];
            let ratio = AffineRatio::new(-value(name) * scale, -slopes[name] * scale, rates[i]);
            normal[i][j] = ratio;
            normal[j][i] = ratio;
        }
        let shear = if isotropic {
            let denominator = 1.0 + value("poisson_ratio");
            [AffineRatio::new(
                young[0] / denominator,
                slopes["young_modulus"] / denominator,
                slopes["poisson_ratio"] / denominator,
            ); 3]
        } else {
            ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES
                .map(|name| AffineRatio::new(2.0 * value(name), 2.0 * slopes[name], 0.0))
        };
        let model = Self {
            reference_temperature_k,
            roots,
            normal,
            shear,
            isotropic,
            rotation: None,
        };
        let bounds = model.interval_bounds(interval)?;
        Ok((model, bounds))
    }

    fn interval_bounds(&self, interval: [f64; 2]) -> Result<TensorTemperatureBounds> {
        let mut lower = f64::INFINITY;
        let mut first_norm = 0.0_f64;
        let mut second_norm = 0.0_f64;
        let guard = 8192.0 * f64::EPSILON;
        // Each ratio is monotone, and its derivative magnitudes are constant
        // numerators divided by positive affine denominators squared/cubed.
        // Endpoint absolute bounds therefore cover EVERY interior temperature.
        if self.isotropic {
            let diagonal = self.normal[0][0];
            let off = self.normal[0][1];
            // Fixed volumetric/deviatoric eigenspaces give exact scalar
            // bounds, including auxetic nu<-0.5 where Gershgorin is insufficient.
            for weight in [2.0_f64, -1.0] {
                let eigenvalue = AffineRatio::new(
                    weight.mul_add(off.a, diagonal.a),
                    weight.mul_add(off.b, diagonal.b),
                    diagonal.d,
                );
                for t in interval {
                    let jet = eigenvalue.jet(t - self.reference_temperature_k)?;
                    lower = lower.min((1.0 - guard) * jet[0]);
                    first_norm = first_norm.max((1.0 + guard) * jet[1].abs());
                    second_norm = second_norm.max((1.0 + guard) * jet[2].abs());
                }
            }
        } else {
            for i in 0..3 {
                let mut diagonal = 0.0;
                let mut radius = 0.0;
                let mut first_row = 0.0;
                let mut second_row = 0.0;
                for j in 0..3 {
                    let a = self.normal[i][j].jet(interval[0] - self.reference_temperature_k)?;
                    let b = self.normal[i][j].jet(interval[1] - self.reference_temperature_k)?;
                    if i == j {
                        diagonal = a[0].min(b[0]);
                    } else {
                        radius += a[0].abs().max(b[0].abs());
                    }
                    first_row += a[1].abs().max(b[1].abs());
                    second_row += a[2].abs().max(b[2].abs());
                }
                lower = lower.min((1.0 - guard) * diagonal - (1.0 + guard) * radius);
                first_norm = first_norm.max(first_row * (1.0 + guard));
                second_norm = second_norm.max(second_row * (1.0 + guard));
            }
        }
        // Either exact isotropic eigenvalues or symmetric diagonal dominance
        // proves SPD and ||A^-1||<=1/lower. The orthotropic sufficient proof
        // intentionally refuses some valid non-dominant paths.
        if !lower.is_finite() || lower <= 0.0 {
            return Err(TensorSourceError::Invalid(
                "complete thermal interval lacks a positive scaled-compliance bound".into(),
            ));
        }
        let inverse = 1.0 / lower;
        let scale = self.roots.iter().map(|v| v * v).fold(0.0_f64, f64::max);
        let rate = first_norm * inverse;
        let mut bounds = [
            scale * inverse,
            scale * inverse * rate,
            scale * inverse * (2.0 * rate * rate + second_norm * inverse),
        ];
        if bounds.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(TensorSourceError::Invalid(
                "thermal compliance interval bounds overflow".into(),
            ));
        }
        for shear in self.shear {
            for t in interval {
                let jet = shear.jet(t - self.reference_temperature_k)?;
                if jet[0] <= 0.0 {
                    return Err(TensorSourceError::Invalid(
                        "thermal shear modulus must stay positive over its original interval"
                            .into(),
                    ));
                }
                for k in 0..3 {
                    bounds[k] = bounds[k].max(jet[k].abs());
                }
            }
        }
        for bound in &mut bounds {
            *bound *= 1.0 + guard;
        }
        if bounds.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(TensorSourceError::Invalid(
                "thermal compliance interval bounds overflow".into(),
            ));
        }
        Ok(TensorTemperatureBounds {
            stiffness_norm_pa: bounds[0],
            slope_norm_pa_k: bounds[1],
            curvature_norm_pa_k2: bounds[2],
        })
    }

    pub(super) fn principal_jet(&self, temperature_k: f64) -> Result<[Tangent6; 3]> {
        let dt = temperature_k - self.reference_temperature_k;
        let mut a = [[[0.0; 3]; 3]; 3];
        for i in 0..3 {
            for j in i..3 {
                let jet = self.normal[i][j].jet(dt)?;
                for k in 0..3 {
                    a[k][i][j] = jet[k];
                    a[k][j][i] = jet[k];
                }
            }
        }
        let inverse = inverse_spd(a[0])?;
        let product = multiply(inverse, a[1]);
        let first = multiply(product, inverse);
        let repeated = multiply(product, first);
        let second = multiply(multiply(inverse, a[2]), inverse);
        let mut out = [[[0.0; 6]; 6]; 3];
        for i in 0..3 {
            for j in i..3 {
                let scale = self.roots[i] * self.roots[j];
                let values = [
                    inverse[i][j] * scale,
                    -first[i][j] * scale,
                    (2.0 * repeated[i][j] - second[i][j]) * scale,
                ];
                for k in 0..3 {
                    out[k][i][j] = values[k];
                    out[k][j][i] = values[k];
                }
            }
        }
        for i in 0..3 {
            let jet = self.shear[i].jet(dt)?;
            for k in 0..3 {
                out[k][i + 3][i + 3] = jet[k];
            }
        }
        if out.iter().flatten().flatten().any(|v| !v.is_finite()) {
            return Err(TensorSourceError::Invalid(
                "thermal inverse-compliance jet overflow".into(),
            ));
        }
        Ok(out)
    }

    pub(super) fn in_frame(mut self, rotation: [[f64; 3]; 3]) -> Self {
        if !self.isotropic {
            self.rotation = Some(crate::tensor::mandel_rotation(rotation));
        }
        self
    }

    /// The Bond matrix is cached once; each snapshot uses two fixed-size products.
    pub(super) fn rotate(&self, principal: Tangent6) -> Tangent6 {
        let Some(q) = self.rotation else {
            return principal;
        };
        let mut left = [[0.0; 6]; 6];
        let mut out = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                for k in 0..6 {
                    left[i][j] = q[i][k].mul_add(principal[k][j], left[i][j]);
                }
            }
        }
        for i in 0..6 {
            for j in i..6 {
                for k in 0..6 {
                    out[i][j] = left[i][k].mul_add(q[j][k], out[i][j]);
                }
                out[j][i] = out[i][j];
            }
        }
        out
    }
}

fn multiply(a: Matrix3, b: Matrix3) -> Matrix3 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| (0..3).fold(0.0, |value, k| a[i][k].mul_add(b[k][j], value)))
    })
}

/// Cholesky solves on the dimensionless SPD block, without pressure cofactors.
fn inverse_spd(a: Matrix3) -> Result<Matrix3> {
    let mut lower = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..=i {
            let mut value = a[i][j];
            for k in 0..j {
                value -= lower[i][k] * lower[j][k];
            }
            lower[i][j] = if i == j {
                if !value.is_finite() || value <= 0.0 {
                    return Err(TensorSourceError::Invalid(
                        "thermal normal compliance lost positive definiteness".into(),
                    ));
                }
                det::sqrt(value)
            } else {
                value / lower[j][j]
            };
        }
    }
    let mut inverse = [[0.0; 3]; 3];
    for column in 0..3 {
        let mut y = [0.0; 3];
        for i in 0..3 {
            let mut value = if i == column { 1.0 } else { 0.0 };
            for k in 0..i {
                value -= lower[i][k] * y[k];
            }
            y[i] = value / lower[i][i];
        }
        for i in (0..3).rev() {
            let mut value = y[i];
            for k in i + 1..3 {
                value -= lower[k][i] * inverse[k][column];
            }
            inverse[i][column] = value / lower[i][i];
        }
    }
    for i in 0..3 {
        for j in 0..i {
            inverse[i][j] = inverse[j][i];
        }
    }
    Ok(inverse)
}
