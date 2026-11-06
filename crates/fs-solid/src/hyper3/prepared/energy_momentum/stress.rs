//! Objective strain-space discrete gradient using the original material law.
//!
//! Cbar=(C0+C1)/2 is evaluated through F*=chol(Cbar)^T. It is NOT C(Fbar).
//! Symmetric algorithmic second Piola stress supplies potential work within
//! arithmetic/continuous-limit tolerance and zero midpoint internal torque. The analytic
//! tangent differentiates this same map with respect to endpoint F1.

use fs_la::factor::{Cholesky, Lu, cholesky, lu};
use fs_material::MaterialError;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_math::dd::Dd;

use crate::hyper3::pressure::determinant_increment;

pub(super) struct AlgorithmicStress {
    pub(super) piola: [f64; 9],
    pub(super) tangent_after: Option<[[f64; 9]; 9]>,
}

fn invalid(what: &str) -> MaterialError {
    MaterialError::State { what: what.into() }
}
fn dot(a: &[f64; 9], b: &[f64; 9]) -> f64 {
    a.iter().zip(b).fold(0.0, |s, (x, y)| x.mul_add(*y, s))
}
fn dot_dd(a: &[Dd; 9], b: &[Dd; 9]) -> Dd {
    a.iter().zip(b).fold(Dd::ZERO, |sum, (x, y)| sum + *x * *y)
}
fn product_dd(a: &[Dd; 9], b: &[Dd; 9]) -> [Dd; 9] {
    core::array::from_fn(|ij| {
        (0..3).fold(Dd::ZERO, |sum, k| {
            sum + a[3 * (ij / 3) + k] * b[3 * k + ij % 3]
        })
    })
}
fn transpose(a: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| a[3 * (ij % 3) + ij / 3])
}
fn product(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3).fold(0.0, |s, k| {
            a[3 * (ij / 3) + k].mul_add(b[3 * k + ij % 3], s)
        })
    })
}
fn symmetric(a: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| 0.5 * a[ij] + 0.5 * a[3 * (ij % 3) + ij / 3])
}
fn left_solve(factor: &Lu, rhs: &[f64; 9]) -> [f64; 9] {
    let mut result = [0.0; 9];
    for column in 0..3 {
        let mut x = core::array::from_fn::<_, 3, _>(|row| rhs[3 * row + column]);
        factor.solve(&mut x);
        for row in 0..3 {
            result[3 * row + column] = x[row];
        }
    }
    result
}

// Differentiate the native factor C=L L^T; return d(F*)=d(L^T).
fn cholesky_direction(factor: &Cholesky, dc: &[f64; 9]) -> [f64; 9] {
    let mut dl = [0.0; 9];
    for column in 0..3 {
        let mut diagonal = dc[3 * column + column];
        for k in 0..column {
            diagonal -= 2.0 * factor.l(column, k) * dl[3 * column + k];
        }
        dl[3 * column + column] = diagonal / (2.0 * factor.l(column, column));
        for row in column + 1..3 {
            let mut value = dc[3 * row + column];
            for k in 0..column {
                value -=
                    dl[3 * row + k] * factor.l(column, k) + factor.l(row, k) * dl[3 * column + k];
            }
            value -= factor.l(row, column) * dl[3 * column + column];
            dl[3 * row + column] = value / factor.l(column, column);
        }
    }
    transpose(&dl)
}

pub(super) fn evaluate<L: HyperelasticMaterial>(
    material: &L,
    before: &[f64; 9],
    after: &[f64; 9],
    tangent: bool,
) -> Result<AlgorithmicStress, MaterialError> {
    evaluate_with_increment(material, before, after, tangent, || {
        material.energy_density_increment(before, after)
    })
}

// Keep the endpoint work target separate from absolute storage. The material
// owner still supplies its own energy difference; pressure adds p*Delta J
// before either endpoint's large p*J value rounds away that material work.
// Evaluate the increment only after the original endpoint domain admissions.
fn evaluate_with_increment<L: HyperelasticMaterial>(
    material: &L,
    before: &[f64; 9],
    after: &[f64; 9],
    tangent: bool,
    increment: impl FnOnce() -> f64,
) -> Result<AlgorithmicStress, MaterialError> {
    // Both endpoints retain the original source and principal-stretch gates.
    material.first_piola(before)?;
    let p_after = material.first_piola(after)?;
    let energy_change = increment();
    let midpoint = core::array::from_fn(|i| {
        Dd::from_f64(0.5) * (Dd::from_f64(before[i]) + Dd::from_f64(after[i]))
    });
    let delta_f: [Dd; 9] =
        core::array::from_fn(|i| Dd::from_f64(after[i]) - Dd::from_f64(before[i]));
    let c0 = product(&transpose(before), before);
    let c1 = product(&transpose(after), after);
    let cbar = core::array::from_fn(|i| 0.5 * c0[i] + 0.5 * c1[i]);
    // Keep the bilinear metric increment and midpoint tails through the work
    // correction AND pullback. A large hydrostatic stress otherwise cancels
    // weak isochoric work even when the material's energy increment is exact.
    let delta_c = core::array::from_fn(|ij| {
        let (row, column) = (ij / 3, ij % 3);
        let (i, j) = (row.min(column), row.max(column));
        (0..3).fold(Dd::ZERO, |sum, k| {
            sum + midpoint[3 * k + i] * delta_f[3 * k + j]
                + delta_f[3 * k + i] * midpoint[3 * k + j]
        })
    });
    if !energy_change.is_finite()
        || cbar.iter().any(|v| !v.is_finite())
        || delta_c.iter().any(|v| !v.to_f64().is_finite())
    {
        return Err(invalid(
            "unrepresentable discrete-gradient metric or energy",
        ));
    }
    let factor =
        cholesky(&cbar, 3).map_err(|_| invalid("positive mean material metric required"))?;
    let representative = core::array::from_fn(|ij| factor.l(ij % 3, ij / 3));
    let inverse =
        lu(&representative, 3).map_err(|_| invalid("singular material metric representative"))?;
    let p = material.first_piola(&representative)?;
    let raw_s = left_solve(&inverse, &p);
    let s = symmetric(&raw_s);
    let denominator = dot_dd(&delta_c, &delta_c);
    let metric_norm = dot(&cbar, &cbar);
    if !denominator.to_f64().is_finite()
        || !metric_norm.is_finite()
        || s.iter().any(|v| !v.is_finite())
    {
        return Err(invalid("unrepresentable discrete-gradient stress"));
    }
    // Below this relative strain threshold the energy quotient loses digits.
    // Use its continuous limit; the caller still gates measured energy error.
    let resolved = denominator.to_f64() > 64.0 * f64::EPSILON * metric_norm;
    let s_dd = s.map(Dd::from_f64);
    let correction = if resolved {
        (Dd::from_f64(2.0) * Dd::from_f64(energy_change) - dot_dd(&s_dd, &delta_c)) / denominator
    } else {
        Dd::ZERO
    };
    let algorithmic_s = core::array::from_fn(|i| correction * delta_c[i] + s_dd[i]);
    let piola = product_dd(&midpoint, &algorithmic_s).map(Dd::to_f64);
    if piola.iter().any(|v| !v.is_finite()) {
        return Err(invalid("unrepresentable algorithmic Piola stress"));
    }
    let tangent_after = if tangent {
        let material_tangent = material.consistent_tangent(&representative)?;
        let mut result = [[0.0; 9]; 9];
        for column in 0..9 {
            let mut direction = [0.0; 9];
            direction[column] = 1.0;
            let term = product(&transpose(after), &direction);
            let d_delta_c = core::array::from_fn(|ij| term[ij] + term[3 * (ij % 3) + ij / 3]);
            let d_cbar = d_delta_c.map(|x| 0.5 * x);
            let d_representative = cholesky_direction(&factor, &d_cbar);
            let d_p: [f64; 9] =
                core::array::from_fn(|row| dot(&material_tangent[row], &d_representative));
            let frame_term = product(&d_representative, &raw_s);
            let d_raw_s = left_solve(&inverse, &core::array::from_fn(|i| d_p[i] - frame_term[i]));
            let d_s = symmetric(&d_raw_s);
            let d_s_dd = d_s.map(Dd::from_f64);
            let d_delta_c_dd = d_delta_c.map(Dd::from_f64);
            let d_correction = if resolved {
                let numerator = Dd::from_f64(2.0) * Dd::from_f64(p_after[column])
                    - dot_dd(&d_s_dd, &delta_c)
                    - dot_dd(&s_dd, &d_delta_c_dd);
                (numerator - correction * Dd::from_f64(2.0) * dot_dd(&delta_c, &d_delta_c_dd))
                    / denominator
            } else {
                Dd::ZERO
            };
            let d_algorithmic_s = core::array::from_fn(|i| {
                d_correction * delta_c[i] + correction * d_delta_c_dd[i] + d_s_dd[i]
            });
            let first = product_dd(&direction.map(Dd::from_f64), &algorithmic_s);
            let second = product_dd(&midpoint, &d_algorithmic_s);
            for row in 0..9 {
                result[row][column] = (Dd::from_f64(0.5) * first[row] + second[row]).to_f64();
            }
        }
        if result.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable algorithmic material tangent"));
        }
        Some(result)
    } else {
        None
    };
    Ok(AlgorithmicStress {
        piola,
        tangent_after,
    })
}

/// Add the conservative closed-boundary pressure potential to the original
/// frozen material. Positive pressure is inward, so INTERNAL storage is +p J;
/// the caller reports external pressure work separately as -p Delta V.
/// The same symmetric strain-space discrete gradient supplies its force and
/// endpoint tangent. Integrating face forces along the displacement segment
/// would close work alone, without ensuring zero torque about midpoint nodes.
pub(super) fn evaluate_with_uniform_pressure<L: HyperelasticMaterial>(
    material: &L,
    before: &[f64; 9],
    after: &[f64; 9],
    pressure_pa: f64,
    tangent: bool,
) -> Result<AlgorithmicStress, MaterialError> {
    if !pressure_pa.is_finite() {
        return Err(invalid("finite uniform pressure required"));
    }
    if pressure_pa == 0.0 {
        return evaluate(material, before, after, tangent);
    }
    evaluate_with_increment(
        &Pressurized {
            material,
            pressure_pa,
        },
        before,
        after,
        tangent,
        || {
            let material_change = material.energy_density_increment(before, after);
            let pressure_change = Dd::from_f64(pressure_pa) * determinant_increment(before, after);
            (Dd::from_f64(material_change) + pressure_change).to_f64()
        },
    )
}

/// Coefficient derivative of the SAME objective strain-space pressure term.
/// The caller has already admitted both endpoints and the representative with
/// the original material response. This is the geometric potential det(F),
/// not another constitutive law or a replacement material-domain admission.
pub(super) fn volume_stress<L: HyperelasticMaterial>(
    material: &L,
    before: &[f64; 9],
    after: &[f64; 9],
) -> Result<[f64; 9], MaterialError> {
    struct Volume<'a>(&'a fs_material::hyper::Hyperelastic);
    impl HyperelasticMaterial for Volume<'_> {
        fn hyperelastic(&self) -> &fs_material::hyper::Hyperelastic {
            self.0
        }
        fn energy_density(&self, f: &[f64; 9]) -> f64 {
            let cof = cofactor(f);
            f[0] * cof[0] + f[1] * cof[1] + f[2] * cof[2]
        }
        fn first_piola(&self, f: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
            Ok(cofactor(f))
        }
        fn consistent_tangent(&self, f: &[f64; 9]) -> Result<[[f64; 9]; 9], MaterialError> {
            let mut result = [[0.0; 9]; 9];
            for column in 0..9 {
                let mut direction = [0.0; 9];
                direction[column] = 1.0;
                let derivative = cofactor_direction(f, &direction);
                for row in 0..9 {
                    result[row][column] = derivative[row];
                }
            }
            Ok(result)
        }
    }
    Ok(evaluate_with_increment(
        &Volume(material.hyperelastic()),
        before,
        after,
        false,
        || determinant_increment(before, after).to_f64(),
    )?
    .piola)
}

struct Pressurized<'a, L> {
    material: &'a L,
    pressure_pa: f64,
}

// Analytic first/second derivatives of det(F), the material-volume density.
// Row-major cofactors avoid inverse-F cancellation in the pressure Hessian.
pub(super) fn cofactor(f: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        let (a, b) = ((ij / 3 + 1) % 3, (ij / 3 + 2) % 3);
        let (c, d) = ((ij % 3 + 1) % 3, (ij % 3 + 2) % 3);
        f[3 * a + c] * f[3 * b + d] - f[3 * a + d] * f[3 * b + c]
    })
}
fn cofactor_direction(f: &[f64; 9], direction: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        let (a, b) = ((ij / 3 + 1) % 3, (ij / 3 + 2) % 3);
        let (c, d) = ((ij % 3 + 1) % 3, (ij % 3 + 2) % 3);
        direction[3 * a + c] * f[3 * b + d] + f[3 * a + c] * direction[3 * b + d]
            - direction[3 * a + d] * f[3 * b + c]
            - f[3 * a + d] * direction[3 * b + c]
    })
}

impl<L: HyperelasticMaterial> HyperelasticMaterial for Pressurized<'_, L> {
    fn hyperelastic(&self) -> &fs_material::hyper::Hyperelastic {
        self.material.hyperelastic()
    }

    fn energy_density(&self, f: &[f64; 9]) -> f64 {
        let cof = cofactor(f);
        let jacobian = f[0] * cof[0] + f[1] * cof[1] + f[2] * cof[2];
        self.pressure_pa
            .mul_add(jacobian, self.material.energy_density(f))
    }

    fn first_piola(&self, f: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
        let original = self.material.first_piola(f)?;
        let cof = cofactor(f);
        Ok(core::array::from_fn(|i| {
            self.pressure_pa.mul_add(cof[i], original[i])
        }))
    }

    fn consistent_tangent(&self, f: &[f64; 9]) -> Result<[[f64; 9]; 9], MaterialError> {
        let mut tangent = self.material.consistent_tangent(f)?;
        for column in 0..9 {
            let mut direction = [0.0; 9];
            direction[column] = 1.0;
            let derivative = cofactor_direction(f, &direction);
            for row in 0..9 {
                tangent[row][column] = self
                    .pressure_pa
                    .mul_add(derivative[row], tangent[row][column]);
            }
        }
        Ok(tangent)
    }
}

#[cfg(test)]
mod tests;
