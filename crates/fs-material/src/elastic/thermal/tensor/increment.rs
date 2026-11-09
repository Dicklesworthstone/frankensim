//! Stable finite differences of the SAME admitted affine engineering sources.
//! These are coefficient increments, not fitted secants or a second material law.

use super::*;
use fs_math::dd::Dd;

type Matrix = [[Dd; 6]; 6];

pub(in crate::elastic::thermal) enum ThermalCoefficientIncrement {
    Isotropic { moduli: [Dd; 2], first: [Dd; 2] },
    Tensor { stiffness: Matrix, first: Matrix },
}

fn dd(x: f64) -> Dd {
    Dd::from_f64(x)
}

#[derive(Clone, Copy)]
struct Linear {
    values: [Dd; 2],
    slope: Dd,
}

impl Linear {
    fn at(reference: f64, slope: f64, reference_t: f64, t0: f64, dt: Dd) -> Self {
        let slope = dd(slope);
        let before = dd(reference) + slope * (dd(t0) - dd(reference_t));
        Self { values: [before, before + slope * dt], slope }
    }

    fn constant(value: f64) -> Self {
        Self { values: [dd(value); 2], slope: Dd::ZERO }
    }

    fn scaled(self, factor: f64) -> Self {
        Self {
            values: self.values.map(|x| x * dd(factor)),
            slope: self.slope * dd(factor),
        }
    }

    fn minus(self, other: Self) -> Self {
        Self {
            values: core::array::from_fn(|i| self.values[i] - other.values[i]),
            slope: self.slope - other.slope,
        }
    }
}

struct Ratio {
    values: [Dd; 2],
    first: [Dd; 2],
    delta: Dd,
    delta_first: Dd,
}

// For affine a,b, a'b-ab' is constant. Form increments before dividing,
// avoiding subtraction of nearly identical ratios or derivative endpoints.
fn ratio(a: Linear, b: Linear, dt: Dd) -> Ratio {
    let numerator = a.slope * b.values[0] - a.values[0] * b.slope;
    let delta = (numerator / b.values[0]) * (dt / b.values[1]);
    Ratio {
        values: core::array::from_fn(|i| a.values[i] / b.values[i]),
        first: core::array::from_fn(|i| (numerator / b.values[i]) / b.values[i]),
        delta,
        delta_first: -delta * (b.slope / b.values[0] + b.slope / b.values[1]),
    }
}

fn isotropic(e: Linear, nu: Linear, dt: Dd) -> ThermalCoefficientIncrement {
    let bulk = ratio(e, Linear::constant(3.0).minus(nu.scaled(6.0)), dt);
    let shear = ratio(e, Linear::constant(2.0).minus(nu.scaled(-2.0)), dt);
    ThermalCoefficientIncrement::Isotropic {
        moduli: [bulk.delta, shear.delta],
        first: [bulk.delta_first, shear.delta_first],
    }
}

fn product(a: &Matrix, b: &Matrix) -> Matrix {
    core::array::from_fn(|i| core::array::from_fn(|j| {
        (0..6).fold(Dd::ZERO, |sum, k| sum + a[i][k] * b[k][j])
    }))
}

fn rotate(a: &Matrix, axes: [[f64; 3]; 3]) -> Matrix {
    if axes == IDENTITY {
        return *a;
    }
    let q = mandel_rotation(axes).map(|row| row.map(dd));
    let qt = core::array::from_fn(|i| core::array::from_fn(|j| q[j][i]));
    product(&product(&q, a), &qt)
}

// Use the same native principal inverse as the original endpoint. Forming a
// global compliance inverse from an approximately orthogonal represented
// rotation would introduce a different tensor; rotate the increments last.
fn principal(
    law: &AffineTensorThermoelastic,
    source: &ResolvedThermoelastic,
) -> Result<Matrix, ThermoelasticError> {
    let mut values = [0.0; 9];
    for (i, coefficient) in law.elastic_coefficients.iter().enumerate() {
        values[i] = value(&source.elastic, coefficient.name, coefficient.dims)?;
    }
    let material = OrthotropicElastic::new(
        [values[0], values[1], values[2]],
        [values[3], values[4], values[5]],
        [values[6], values[7], values[8]],
        law.maximum_absolute_mandel_strain,
    )?;
    Ok(material.stiffness().map(|row| row.map(dd)))
}

pub(in crate::elastic::thermal) fn coefficient_increments(
    law: &AffineThermoelastic,
    before: &ResolvedThermoelastic,
    after: &ResolvedThermoelastic,
) -> Result<ThermalCoefficientIncrement, ThermoelasticError> {
    let t0 = before.temperature_k();
    let dt = dd(after.temperature_k()) - dd(t0);
    if dt == Dd::ZERO {
        let tensor = matches!(law, AffineThermoelastic::Tensor(source)
            if matches!(source.stiffness, ThermoelasticSymmetry::OrientedOrthotropic { .. }));
        return Ok(if tensor {
            ThermalCoefficientIncrement::Tensor {
                stiffness: [[Dd::ZERO; 6]; 6], first: [[Dd::ZERO; 6]; 6],
            }
        } else {
            ThermalCoefficientIncrement::Isotropic { moduli: [Dd::ZERO; 2], first: [Dd::ZERO; 2] }
        });
    }
    if let AffineThermoelastic::Isotropic(law) = law {
        return Ok(isotropic(
            Linear::at(
                value(&law.elastic_reference, YOUNG_MODULUS_PROPERTY, Pressure::DIMS)?,
                law.youngs_slope_pa_k, law.reference_temperature_k, t0, dt,
            ),
            Linear::at(
                law.poisson, law.poisson_slope_per_k, law.reference_temperature_k, t0, dt,
            ),
            dt,
        ));
    }
    let AffineThermoelastic::Tensor(law) = law else { unreachable!() };
    let coefficient = |i: usize| {
        let source = &law.elastic_coefficients[i];
        Linear::at(source.reference, source.slope, law.reference_temperature_k, t0, dt)
    };
    match law.stiffness {
        ThermoelasticSymmetry::Isotropic => {
            return Ok(isotropic(coefficient(0), coefficient(1), dt));
        }
        ThermoelasticSymmetry::IsotropicBulkShear => {
            return Ok(ThermalCoefficientIncrement::Isotropic {
                moduli: [coefficient(0).slope * dt, coefficient(1).slope * dt],
                first: [Dd::ZERO; 2],
            });
        }
        ThermoelasticSymmetry::IsotropicYoungShear => {
            let (e, g) = (coefficient(0), coefficient(1));
            // K = E/3 * [G/(3G-E)]: a product of one affine path and one
            // affine ratio. Polarization avoids products of squared moduli.
            let r = ratio(g, g.scaled(3.0).minus(e), dt);
            let mean_e = (e.values[0] + e.values[1]) * dd(0.5);
            let mean_r = (r.values[0] + r.values[1]) * dd(0.5);
            let mean_first = (r.first[0] + r.first[1]) * dd(0.5);
            let delta_e = e.slope * dt;
            return Ok(ThermalCoefficientIncrement::Isotropic {
                moduli: [(mean_e * r.delta + mean_r * delta_e) / dd(3.0), g.slope * dt],
                first: [(e.slope * r.delta + mean_e * r.delta_first
                    + mean_first * delta_e) / dd(3.0), Dd::ZERO],
            });
        }
        ThermoelasticSymmetry::OrientedOrthotropic { .. } => {}
    }

    let (mut ds, mut st0, mut st1, mut dst) = (
        [[Dd::ZERO; 6]; 6], [[Dd::ZERO; 6]; 6],
        [[Dd::ZERO; 6]; 6], [[Dd::ZERO; 6]; 6],
    );
    let mut insert = |i: usize, j: usize, r: Ratio| {
        for (row, col) in [(i, j), (j, i)] {
            ds[row][col] = r.delta;
            st0[row][col] = r.first[0];
            st1[row][col] = r.first[1];
            dst[row][col] = r.delta_first;
        }
    };
    for i in 0..3 {
        insert(i, i, ratio(Linear::constant(1.0), coefficient(i), dt));
        insert(i + 3, i + 3, ratio(Linear::constant(0.5), coefficient(i + 6), dt));
    }
    for (k, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
        insert(i, j, ratio(coefficient(k + 3).scaled(-1.0), coefficient(i), dt));
    }
    let c0 = principal(law, before)?;
    let c1 = principal(law, after)?;
    // Delta C = -C1 Delta S C0; every term retains Delta T before
    // multiplication. Delta C_T differentiates the same inverse identity.
    let dc = product(&product(&c1, &ds), &c0).map(|row| row.map(|x| -x));
    let a = product(&product(&dc, &st1), &c1);
    let b = product(&product(&c0, &dst), &c1);
    let c = product(&product(&c0, &st0), &dc);
    let dct = core::array::from_fn(|i| core::array::from_fn(|j| {
        -(a[i][j] + b[i][j] + c[i][j])
    }));
    Ok(ThermalCoefficientIncrement::Tensor {
        stiffness: rotate(&dc, law.stiffness.rotation()),
        first: rotate(&dct, law.stiffness.rotation()),
    })
}
