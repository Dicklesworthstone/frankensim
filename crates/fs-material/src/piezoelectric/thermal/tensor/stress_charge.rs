//! Exact thermal reduction of affine stress-charge tensors.
//!
//! In the fixed frame, write r=mandel(n n^T), h=e^T n, beta=n.epsilonS.n
//! and U=[r,h]. The reduced compliance is A=U^T C^-1 U+diag(0,beta).
//! For v=C^-1 r, w=C^-1 h and J=[-C'v,h'-C'w]:
//! A'=[[-v.C'v, h'.v-v.C'w],
//!     [h'.v-v.C'w, beta'+2h'.w-w.C'w]], A''=2 J^T C^-1 J.
//! These are derivatives of the complete tensor inverse; no eigenbasis,
//! diagonal truncation, finite differences or affine inverse assumption.
//!
//! The native scaled Cholesky factor supplies B with A''=2 B B^T. Its
//! nonnegative capacity contribution may be omitted in a lower capacity bound.
//! The bound below controls C^-1, h and beta throughout the affine SOURCE
//! interval; endpoint positivity of the REDUCED non-affine A alone is not used.
//!
//! Source convention and tensor conversion:
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.sme/sme_ug_theory.06.043.html>.
//! Helmholtz thermal derivatives:
//! <https://web.mit.edu/beretta/www/AdvancedThermodynamics_2_43_Lecture6.htm>.

use super::*;
use crate::tensor::{mandel_rotation, rotate_mandel};

fn dot<const N: usize>(a: &[f64; N], b: &[f64; N]) -> f64 {
    a.iter().zip(b).fold(0.0, |sum, (a, b)| a.mul_add(*b, sum))
}
fn apply<const R: usize, const C: usize>(matrix: &[[f64; C]; R], x: &[f64; C]) -> [f64; R] {
    core::array::from_fn(|i| dot(&matrix[i], x))
}
fn norm<const N: usize>(x: &[f64; N]) -> Result<f64, PiezoelectricError> {
    finite(x.iter().fold(0.0_f64, |value, x| value.hypot(*x)))
}
fn infinity_norm<const N: usize>(matrix: &[[f64; N]; N]) -> Result<f64, PiezoelectricError> {
    let mut result = 0.0_f64;
    for row in matrix {
        let sum = finite(row.iter().map(|value| value.abs()).sum())?;
        result = result.max(sum);
    }
    Ok(result)
}
fn admit_same_source_coefficients<const R: usize, const C: usize>(
    original: &[[f64; C]; R],
    refreshed: &[[f64; C]; R],
) -> Result<(), PiezoelectricError> {
    let scale = original
        .iter()
        .flatten()
        .chain(refreshed.iter().flatten())
        .fold(0.0_f64, |scale, value| scale.max(value.abs()));
    if original
        .iter()
        .flatten()
        .zip(refreshed.iter().flatten())
        .any(|(a, b)| (a - b).abs() > 1024.0 * f64::EPSILON * scale)
    {
        return Err(PiezoelectricError::Invalid(
            "thermal stress-charge source does not reproduce the original complete tensor",
        ));
    }
    Ok(())
}
fn mechanical_direction(n: [f64; 3]) -> [f64; 6] {
    let basis = mandel_basis();
    core::array::from_fn(|a| {
        (0..3).fold(0.0, |sum, i| {
            (0..3).fold(sum, |sum, j| (basis[a][i][j] * n[i]).mul_add(n[j], sum))
        })
    })
}
fn electric_stress(e: &[[f64; 6]; 3], n: &[f64; 3]) -> [f64; 6] {
    core::array::from_fn(|a| (0..3).fold(0.0, |sum, i| e[i][a].mul_add(n[i], sum)))
}

impl StressChargeThermalSource {
    fn rotate(&mut self, q: [[f64; 3]; 3]) -> Result<(), PiezoelectricError> {
        // Match the native tensor owner's frame transformation, including exact
        // symmetry of the rotated diagonal blocks. The source roles stay pinned.
        self.stiffness_slope = rotate_mandel(self.stiffness_slope, q);
        let m = mandel_rotation(q);
        let mut coupling = [[0.0; 6]; 3];
        for (i, row) in coupling.iter_mut().enumerate() {
            for (a, value) in row.iter_mut().enumerate() {
                *value = (0..3).fold(0.0, |sum, j| {
                    (0..6).fold(sum, |sum, b| {
                        (q[i][j] * self.coupling_slope[j][b]).mul_add(m[a][b], sum)
                    })
                });
            }
        }
        let mut permittivity = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in i..3 {
                let value = (0..3).fold(0.0, |sum, a| {
                    (0..3).fold(sum, |sum, b| {
                        (q[i][a] * self.permittivity_slope[a][b]).mul_add(q[j][b], sum)
                    })
                });
                permittivity[i][j] = value;
                permittivity[j][i] = value;
            }
        }
        if self
            .stiffness_slope
            .iter()
            .flatten()
            .chain(coupling.iter().flatten())
            .chain(permittivity.iter().flatten())
            .any(|x| !x.is_finite())
        {
            return Err(PiezoelectricError::Invalid(
                "rotated piezoelectric source derivative is unrepresentable",
            ));
        }
        self.coupling_slope = coupling;
        self.permittivity_slope = permittivity;
        Ok(())
    }

    pub(super) fn derivatives(
        &self,
        tensor: &TensorPiezoelectric,
        n: [f64; 3],
    ) -> Result<TensorThermalDerivatives, PiezoelectricError> {
        let r = mechanical_direction(n);
        let h = electric_stress(tensor.coupling_c_m2(), &n);
        let h_prime = electric_stress(&self.coupling_slope, &n);
        let v = tensor.solve_elastic(r)?;
        let w = tensor.solve_elastic(h)?;
        let cv = apply(&self.stiffness_slope, &v);
        let cw = apply(&self.stiffness_slope, &w);
        let beta_prime = dot(&n, &apply(&self.permittivity_slope, &n));
        let slopes = [
            finite(-dot(&v, &cv))?,
            finite(dot(&h_prime, &v) - dot(&v, &cw))?,
            finite(2.0_f64.mul_add(dot(&h_prime, &w), beta_prime - dot(&w, &cw)))?,
        ];
        let j = [
            cv.map(|value| -value),
            core::array::from_fn(|i| h_prime[i] - cw[i]),
        ];
        let curvature_factor = [
            tensor.elastic_inverse_factor_action(j[0])?,
            tensor.elastic_inverse_factor_action(j[1])?,
        ];
        Ok(TensorThermalDerivatives {
            slopes,
            curvature_factor,
        })
    }
}

impl ThermalLongitudinalPiezoelectric {
    /// Bind complete affine stress-charge sources C^E(T), e(T), epsilon^S(T)
    /// to the existing longitudinal rod/circuit thermal law.
    ///
    /// Every original tensor role is admitted over the entire interval, and
    /// both full C^E and epsilon^S endpoint tensors must be positive definite.
    /// Their affine interpolation then stays stable for every temperature.
    /// The resulting s^E, d and epsilon^T are generally NON-affine: native
    /// factor solves produce their actual first and second derivatives at each
    /// trial temperature. No coefficient or slope is frozen at the reference.
    ///
    /// The fixed homogeneous constraint is T=t n n^T and E=E_n n, including
    /// transverse relaxation. Original sources, frame, electrode direction and
    /// other query coordinates survive refresh. Repeated caloric responses
    /// use a scalar inverse and six stack-sized curvature actions.
    ///
    /// This reversible small-signal law has no pyroelectric offset, thermal
    /// expansion, phase/depoling model or empirical accuracy claim. Its interval
    /// capacity bound is conservative model admission, not a numerical certificate.
    pub fn from_stress_charge_tensor(
        card: &MaterialCard,
        original: &TensorPiezoelectric,
        direction: [f64; 3],
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, PiezoelectricError> {
        let longitudinal = original.longitudinal(direction)?;
        let direction = longitudinal
            .tensor_reduction()
            .expect("native reduction receipt")
            .direction;
        let mut source = StressChargeThermalSource {
            stiffness_slope: [[0.0; 6]; 6],
            coupling_slope: [[0.0; 6]; 3],
            permittivity_slope: [[0.0; 3]; 3],
            reference_curvature_factor: [[0.0; 6]; 2],
        };
        let mut reference_temperature_k = 0.0;
        for (prefix, rows, columns, block) in [
            ("piezoelectric_stiffness_e_mandel", 6, 6, 0),
            ("piezoelectric_stress_charge_mandel", 3, 6, 1),
            ("piezoelectric_permittivity_s", 3, 3, 2),
        ] {
            for i in 0..rows {
                for j in 0..columns {
                    let name = format!("{prefix}_{}{}", i + 1, j + 1);
                    let (reference, slope, _) = affine_source_path(
                        card,
                        original.material(),
                        &name,
                        temperature_axis,
                        temperature_range_k,
                        block != 1 && i == j,
                    )
                    .map_err(field_error)?;
                    reference_temperature_k = reference;
                    match block {
                        0 => source.stiffness_slope[i][j] = slope,
                        1 => source.coupling_slope[i][j] = slope,
                        _ => source.permittivity_slope[i][j] = slope,
                    }
                }
            }
        }
        let rotation = *original.source_to_frame_rotation();
        if rotation != [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
            source.rotate(rotation)?;
        }
        let reference = source.derivatives(original, direction)?;
        source.reference_curvature_factor = reference.curvature_factor;
        let reduction = TensorThermalReduction {
            source_frame: original.source_frame_identity(),
            frame: original.frame_identity(),
            source_to_frame: rotation,
            direction,
            coefficient_roundoff_scales: [0.0; 3],
            stress_charge: Some(source),
        };
        // A retained bundle can contain both tensor source schemas. Verify the
        // requested stress-charge roles reproduce the ORIGINAL full tensor,
        // including unexcited directions, before binding their derivatives.
        let refreshed = reduction.tensor(original.material().clone())?;
        admit_same_source_coefficients(
            original.stiffness_mandel_pa(),
            refreshed.stiffness_mandel_pa(),
        )?;
        admit_same_source_coefficients(original.coupling_c_m2(), refreshed.coupling_c_m2())?;
        admit_same_source_coefficients(original.permittivity_f_m(), refreshed.permittivity_f_m())?;
        Self::prepare(
            card,
            &longitudinal,
            temperature_axis,
            temperature_range_k,
            reference_temperature_k,
            reference.slopes,
            Some(reduction),
        )
    }
}

impl TensorThermalReduction {
    pub(in crate::piezoelectric::thermal) fn stress_charge_capacity_bounds(
        &self,
        law: &ThermalLongitudinalPiezoelectric,
    ) -> Result<(f64, f64), PiezoelectricError> {
        let source = self
            .stress_charge
            .as_ref()
            .expect("stress-charge admission");
        let n = self.direction;
        let r = mechanical_direction(n);
        let mut c_min = f64::INFINITY;
        let mut c_max = 0.0_f64;
        let mut beta_min = f64::INFINITY;
        let mut beta_max = 0.0_f64;
        let mut h_max = 0.0_f64;
        for temperature in law.temperature_range_k {
            let point = law
                .original
                .material()
                .query()
                .with(&law.temperature_axis, temperature)
                .map_err(|_| {
                    PiezoelectricError::Invalid("piezoelectric thermal query is unrepresentable")
                })?;
            let tensor = self.tensor(law.original.material().at_query_point(&law.card, &point)?)?;
            // Full native SPD admission precedes any longitudinal projection.
            tensor.longitudinal(n)?;
            let mut inverse = [[0.0; 6]; 6];
            for j in 0..6 {
                let column = tensor
                    .solve_elastic(core::array::from_fn(|i| if i == j { 1.0 } else { 0.0 }))?;
                for i in 0..6 {
                    inverse[i][j] = column[i];
                }
            }
            let lower = finite(1.0 / infinity_norm(&inverse)?)?;
            if lower <= 0.0 {
                return Err(PiezoelectricError::Invalid(
                    "piezoelectric stiffness lower bound underflows",
                ));
            }
            c_min = c_min.min(lower);
            c_max = c_max.max(infinity_norm(tensor.stiffness_mandel_pa())?);
            // The native positive storage evaluates beta robustly even when
            // signed off-diagonal dielectric terms nearly cancel.
            let beta = finite(2.0 * tensor.field_response([0.0; 6], n)?.stored_energy_j_m3)?;
            beta_min = beta_min.min(beta);
            beta_max = beta_max.max(beta);
            h_max = h_max.max(norm(&electric_stress(tensor.coupling_c_m2(), &n))?);
        }
        let down = 1.0 - 4096.0 * f64::EPSILON;
        let up = 1.0 + 4096.0 * f64::EPSILON;
        c_min = finite(c_min * down)?;
        c_max = finite(c_max * up)?;
        beta_min = finite(beta_min * down)?;
        beta_max = finite(beta_max * up)?;
        h_max = finite(h_max * up)?;
        // Affine SPD C lies above c_min I and below c_max I by endpoint
        // convexity. beta is affine and ||h|| is convex. Hence throughout:
        // s >= ||r||^2/c_max, s <= ||r||^2/c_min,
        // p <= beta_max+||h||_max^2/c_min, det(A) >= s_min*beta_min.
        // For G=S^-1 A S^-1, lambda_min(G) >= det(G)/trace(G).
        let [q0, q1] = law.scales;
        let r_norm = norm(&r)?;
        let r_squared = finite(r_norm * r_norm)?;
        let a_min = finite(r_squared / c_max / q0 / q0)?;
        let a_max = finite(r_squared / c_min / q0 / q0)?;
        let b_min = finite(beta_min / q1 / q1)?;
        let p_max = finite(h_max.mul_add(h_max / c_min, beta_max))?;
        let trace_max = finite(a_max + p_max / q1 / q1)?;
        let minimum = finite((a_min / trace_max) * b_min * down)?;
        if minimum <= 0.0 || c_min <= 0.0 || beta_min <= 0.0 {
            return Err(PiezoelectricError::Invalid(
                "complete piezoelectric thermal interval bound is unrepresentable",
            ));
        }
        // For Uhat=[r/q0,h/q1], differentiate Uhat^T C^-1 Uhat:
        // ||G'|| <= 2||Uhat'|| ||Uhat||/c_min
        //          + ||Uhat||^2 ||C'||/c_min^2 + |beta'|/q1^2.
        // Frobenius norms bound Uhat; infinity norm bounds symmetric C'.
        let h_prime_norm = norm(&electric_stress(&source.coupling_slope, &n))?;
        let c_prime_norm = infinity_norm(&source.stiffness_slope)?;
        let beta_prime = finite(dot(&n, &apply(&source.permittivity_slope, &n)))?;
        let u = finite((r_norm / q0).hypot(h_max / q1))?;
        let u_prime = finite(h_prime_norm / q1)?;
        let derivative = finite(
            (2.0 * u_prime) * (u / c_min)
                + (u / c_min) * u * (c_prime_norm / c_min)
                + beta_prime.abs() / q1 / q1,
        )?;
        let derivative = finite(derivative * up)?;
        if derivative == 0.0 && (h_prime_norm != 0.0 || c_prime_norm != 0.0 || beta_prime != 0.0) {
            return Err(PiezoelectricError::Invalid(
                "piezoelectric interval derivative bound underflows",
            ));
        }
        Ok((minimum, derivative))
    }
}
