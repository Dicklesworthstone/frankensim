//! Local physical parameter tangents of a simple mechanical mode.
//!
//! The binary coordinate/rate maps stay fixed while differentiating the
//! physical polynomial. Shared bordered solves provide the full left/right
//! coupling. Physical recovery differentiates phi^H M phi=1, including dM,
//! and keeps one ORIGINAL displacement component real. No moving maximum,
//! crossing tracker, finite-difference step or material law is inferred.

use super::{
    C64, PreparedQuadraticModes, PreparedQuadraticResponse, QuadraticMode, QuadraticResponseError,
    QuadraticResponseOptions, component_exponent, divide, finite, invalid, poll, product,
    quadratic_fraction, scale, scale_complex,
};
use fs_math::det;

/// One signed direction of ORIGINAL row-major M/C/K per caller-declared unit
/// of theta. All three arrays must be complete and finite. dM is symmetric;
/// it need not be positive definite. dC/dK may be nonsymmetric or indefinite,
/// retaining gyroscopic, circulatory, active and nonproportional directions.
#[derive(Clone, Copy, Debug)]
pub struct QuadraticParameterDirection<'a> {
    /// dM/dtheta, in the original generalized-coordinate mass units.
    pub mass: &'a [f64],
    /// dC/dtheta, including the signed derivative of gyroscopic coefficients.
    pub damping: &'a [f64],
    /// dK/dtheta, including full mutual and nonsymmetric stiffness directions.
    pub stiffness: &'a [f64],
}
impl QuadraticParameterDirection<'_> {
    fn arrays(&self) -> [&[f64]; 3] {
        [self.mass, self.damping, self.stiffness]
    }
}

/// Directional derivative of the published physical mode per unit theta.
/// This differentiates a locally isolated simple pole, not its ordering in a
/// spectrum. The fixed physical phase coordinate belongs to the prepared owner.
#[derive(Clone, Debug, PartialEq)]
pub struct QuadraticModeDerivative {
    /// ds/dtheta in original physical Laplace units; active signs are retained.
    pub laplace_frequency_per_s: C64,
    /// dphi/dtheta with 2 Re(phi^H M dphi)+phi^H dM phi=0 and a fixed phase.
    pub displacement: Vec<C64>,
    /// d(s*phi)/dtheta = ds*phi+s*dphi, in original physical velocity units.
    pub velocity: Vec<C64>,
    /// Shared solver's backward residual of the original bordered rows.
    pub linear_residual: f64,
    /// Independent differentiated ORIGINAL-polynomial row residual, recomputed
    /// from source M/C/K, the supplied direction and published physical outputs.
    pub polynomial_relative_residual: f64,
    /// Relative residual of the differentiated physical unit-mass constraint.
    pub normalization_relative_residual: f64,
}

/// One immutable response factor borrowed alongside its original source pencil.
/// Force residues and multiple sensitivity batches reuse the SAME factor via
/// response(). No source matrices or mass factor are copied. No derivative of
/// the preparation's numerical coordinate/rate scaling is taken.
#[derive(Debug)]
pub struct PreparedQuadraticSensitivities<'model, 'source> {
    source: &'model PreparedQuadraticModes<'source>,
    response: PreparedQuadraticResponse,
    phase_coordinate: usize,
    phase: C64,
    norm: f64,
    right_per_norm: Vec<C64>,
    mass_image_per_norm: Vec<C64>,
}

fn admit(n: usize, options: &QuadraticResponseOptions) -> Result<(), QuadraticResponseError> {
    super::admit(n, options)?;
    let (n, k) = (n as u128, options.numerical.max_vectors as u128);
    super::super::budget(
        "physical sensitivity bytes",
        16 * (2 * n * n + 8 * k * n + 20 * n) + 512 * k,
        options.max_physical_bytes as u128,
    )?;
    super::super::budget(
        "physical sensitivity work",
        256 * (k + 1) * n * n + 512 * k * n,
        options.max_physical_work as u128,
    )?;
    Ok(())
}

impl<'source> PreparedQuadraticModes<'source> {
    /// Prepare one rechecked simple mode for full physical M/C/K sensitivities.
    /// The existing response factor supplies both parameter tangents and force
    /// residues. Original source coefficients and the mass Cholesky are borrowed.
    ///
    /// In addition to the response admission, the same physical ceilings admit
    /// 16*(2*n^2+8*k*n+20*n)+512*k bytes and 256*(k+1)*n^2+512*k*n work units,
    /// with k=max_vectors. These cover preparation and each complete batch;
    /// caller direction arrays and the separately admitted source/factor are
    /// excluded. Cancellation publishes no partial batch or changed source.
    pub fn prepare_sensitivities<'model>(
        &'model self,
        mode: &QuadraticMode,
        options: QuadraticResponseOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedQuadraticSensitivities<'model, 'source>, QuadraticResponseError> {
        poll(&mut cancelled)?;
        admit(self.n, &options)?;
        let response = self.prepare_response(mode, options, &mut cancelled)?;
        let phase_coordinate = (1..self.n).fold(0, |best, i| {
            if response.displacement[i].abs() > response.displacement[best].abs() {
                i
            } else {
                best
            }
        });
        let vector = response.numerical.right_vector();
        let phase = divide(
            vector[phase_coordinate].conj(),
            vector[phase_coordinate].abs(),
        );
        let norm = self.mass_norm(vector, &mut cancelled)?;
        let (mantissa, exponent) = quadratic_fraction(norm);
        let right_per_norm = vector
            .iter()
            .map(|&v| product(v, C64::from_re(mantissa.recip()), -exponent))
            .collect::<Result<Vec<_>, _>>()?;
        let mass_image_per_norm = mass_image(self, &right_per_norm, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(PreparedQuadraticSensitivities {
            source: self,
            response,
            phase_coordinate,
            phase,
            norm,
            right_per_norm,
            mass_image_per_norm,
        })
    }
}

impl PreparedQuadraticSensitivities<'_, '_> {
    /// Original canonical physical mode and shared force-residue factor.
    pub fn response(&self) -> &PreparedQuadraticResponse {
        &self.response
    }

    /// This original displacement component is held real while differentiating.
    /// It is the largest baseline component (lowest-index tie), held FIXED for
    /// all directions. No differentiability is claimed across pivot switches,
    /// repeated poles, crossings or a finite parameter step.
    pub fn phase_reference_coordinate(&self) -> usize {
        self.phase_coordinate
    }

    /// Analytic derivatives for a complete nonempty batch of signed dM/dC/dK.
    /// Each direction uses the existing factor and O(n^2) source actions/audits.
    /// All outputs use the caller's parameter units, including changes of modal
    /// mass. A zero direction gives zero; finite material/geometry changes still
    /// require a new physical model and a new simple-mode preparation.
    pub fn parameter_derivatives(
        &self,
        directions: &[QuadraticParameterDirection<'_>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<QuadraticModeDerivative>, QuadraticResponseError> {
        poll(&mut cancelled)?;
        if directions.is_empty() || directions.len() > self.response.options.numerical.max_vectors {
            return Err(invalid(
                "complete nonempty parameter batch within max_vectors required",
            ));
        }
        let mut actions = Vec::with_capacity(directions.len());
        let mut exponents = Vec::with_capacity(directions.len());
        for direction in directions {
            let exponent = self.admit_direction(direction, &mut cancelled)?;
            actions.push(self.direction_action(direction, exponent, &mut cancelled)?);
            exponents.push(exponent);
        }
        let references: Vec<_> = actions.iter().map(Vec::as_slice).collect();
        let numerical = self
            .response
            .numerical
            .parameter_derivatives(&references, &mut cancelled)?;
        let (norm_fraction, norm_exponent) = quadratic_fraction(self.norm);
        let map = divide(self.phase, norm_fraction);
        let vector = self.response.numerical.right_vector();
        let mut results = Vec::with_capacity(directions.len());
        for ((direction, exponent), tangent) in directions.iter().zip(exponents).zip(numerical) {
            poll(&mut cancelled)?;
            let normalized_slope = tangent
                .vector
                .iter()
                .map(|&v| product(v, C64::from_re(norm_fraction.recip()), -norm_exponent))
                .collect::<Result<Vec<_>, _>>()?;
            let image = mass_image(self.source, &normalized_slope, &mut cancelled)?;
            let mut motion = RealSum::default();
            for (&a, &b) in self.mass_image_per_norm.iter().zip(&image) {
                motion.add(product(a.conj(), b, 0)?.re);
            }
            let (moving_mass, _) = self.mass_direction_form(
                direction,
                exponent,
                &self.right_per_norm,
                &mut cancelled,
            )?;
            let phase_slope =
                (tangent.vector[self.phase_coordinate] / vector[self.phase_coordinate]).im;
            let correction = C64::new(motion.value() + 0.5 * moving_mass, phase_slope);
            if !finite(correction) {
                return Err(invalid(
                    "differentiated mass or fixed phase gauge is not representable",
                ));
            }
            let ds = scale_complex(tangent.value, exponent + self.source.time_exponent)?;
            let mut displacement = Vec::with_capacity(self.source.n);
            for (i, (&v, &dv)) in vector.iter().zip(&tangent.vector).enumerate() {
                poll(&mut cancelled)?;
                let fixed = dv - product(correction, v, 0)?;
                displacement.push(product(
                    fixed,
                    map,
                    exponent - self.source.coordinate_exponents[i] - norm_exponent,
                )?);
            }
            // The declared phase constraint is exact in the published coordinate.
            // The independent original-polynomial audit below checks this recovery.
            displacement[self.phase_coordinate].im = 0.0;
            let velocity = displacement
                .iter()
                .zip(&self.response.displacement)
                .map(|(&dq, &q)| sum_products(self.response.pole, dq, ds, q))
                .collect::<Result<Vec<_>, _>>()?;
            let (polynomial_relative_residual, normalization_relative_residual) =
                self.original_residual(direction, exponent, ds, &displacement, &mut cancelled)?;
            let limit = self.response.options.numerical.linear_residual_tolerance;
            if polynomial_relative_residual > limit || normalization_relative_residual > limit {
                return Err(QuadraticResponseError::SensitivityResidual {
                    polynomial: polynomial_relative_residual,
                    normalization: normalization_relative_residual,
                    limit,
                });
            }
            results.push(QuadraticModeDerivative {
                laplace_frequency_per_s: ds,
                displacement,
                velocity,
                linear_residual: tangent.linear_residual,
                polynomial_relative_residual,
                normalization_relative_residual,
            });
        }
        poll(&mut cancelled)?;
        Ok(results)
    }

    fn coefficient(
        &self,
        value: f64,
        i: usize,
        j: usize,
        order: usize,
        exponent: i32,
    ) -> Result<f64, QuadraticResponseError> {
        Ok(scale(
            value,
            -self.source.coordinate_exponents[i]
                - self.source.coordinate_exponents[j]
                - order as i32 * self.source.time_exponent
                - exponent,
        )?)
    }

    fn admit_direction(
        &self,
        direction: &QuadraticParameterDirection<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<i32, QuadraticResponseError> {
        let n = self.source.n;
        if direction.arrays().iter().any(|a| a.len() != n * n) {
            return Err(invalid(
                "every dM/dC/dK direction must retain the full original matrix",
            ));
        }
        let mut largest = i32::MIN;
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..n {
                if direction.mass[i * n + j] != direction.mass[j * n + i] {
                    return Err(invalid("physical mass direction must be symmetric"));
                }
                for (order, array) in direction.arrays().into_iter().enumerate() {
                    let value = array[i * n + j];
                    if !value.is_finite() {
                        return Err(invalid("finite physical coefficient directions required"));
                    }
                    if value != 0.0 {
                        largest = largest.max(
                            quadratic_fraction(value).1
                                - self.source.coordinate_exponents[i]
                                - self.source.coordinate_exponents[j]
                                - order as i32 * self.source.time_exponent,
                        );
                    }
                }
            }
        }
        Ok(if largest == i32::MIN { 0 } else { largest })
    }

    fn direction_action(
        &self,
        direction: &QuadraticParameterDirection<'_>,
        exponent: i32,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, QuadraticResponseError> {
        let n = self.source.n;
        let s = self.response.numerical.value();
        let mut action = Vec::with_capacity(n);
        for i in 0..n {
            poll(cancelled)?;
            let mut sums = [ComplexSum::default(); 3];
            for (j, &v) in self.response.numerical.right_vector().iter().enumerate() {
                for (order, array) in direction.arrays().into_iter().enumerate() {
                    let a = self.coefficient(array[i * n + j], i, j, order, exponent)?;
                    sums[order].add(v.scale(a));
                }
            }
            let value = s * (s * sums[0].value() + sums[1].value()) + sums[2].value();
            if !finite(value) {
                return Err(invalid("normalized physical parameter action overflows"));
            }
            action.push(value);
        }
        Ok(action)
    }

    fn mass_direction_form(
        &self,
        direction: &QuadraticParameterDirection<'_>,
        exponent: i32,
        vector: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), QuadraticResponseError> {
        let mut sum = RealSum::default();
        let mut absolute = 0.0;
        for (i, &v) in vector.iter().enumerate() {
            poll(cancelled)?;
            for (j, &w) in vector[..=i].iter().enumerate() {
                let a =
                    self.coefficient(direction.mass[i * self.source.n + j], i, j, 0, exponent)?;
                let doubled = i32::from(i != j);
                sum.add(weighted_product(v.conj(), w, a, doubled)?.re);
                absolute += weighted_product(
                    C64::from_re(v.abs()),
                    C64::from_re(w.abs()),
                    a.abs(),
                    doubled,
                )?
                .re;
            }
        }
        if !sum.value().is_finite() || !absolute.is_finite() {
            return Err(invalid("physical mass-direction quadratic form overflows"));
        }
        Ok((sum.value(), absolute))
    }

    fn original_residual(
        &self,
        direction: &QuadraticParameterDirection<'_>,
        exponent: i32,
        physical_ds: C64,
        physical_dq: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), QuadraticResponseError> {
        let n = self.source.n;
        // Decode the PUBLISHED physical values with one common amplitude scale
        // before multiplying coordinate powers. Divide the tangent/direction by
        // its original exact binary normalization, leaving a homogeneous audit.
        let common = self
            .response
            .displacement
            .iter()
            .enumerate()
            .filter_map(|(i, &q)| {
                component_exponent(q).map(|e| e + self.source.coordinate_exponents[i])
            })
            .chain(physical_dq.iter().enumerate().filter_map(|(i, &q)| {
                component_exponent(q).map(|e| e + self.source.coordinate_exponents[i] - exponent)
            }))
            .max()
            .ok_or_else(|| invalid("nonzero prepared physical mode required"))?;
        let q = self
            .response
            .displacement
            .iter()
            .enumerate()
            .map(|(i, &v)| scale_complex(v, self.source.coordinate_exponents[i] - common))
            .collect::<Result<Vec<_>, _>>()?;
        let dq = physical_dq
            .iter()
            .enumerate()
            .map(|(i, &v)| {
                scale_complex(v, self.source.coordinate_exponents[i] - common - exponent)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let ds = scale_complex(physical_ds, -self.source.time_exponent - exponent)?;
        let s = scale_complex(self.response.pole, -self.source.time_exponent)?;
        let radius = s.abs().max(self.source.options.spectral_scale);
        let q_norm = q.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        let dq_norm = dq.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        let mut worst = 0.0_f64;
        for i in 0..n {
            poll(cancelled)?;
            let mut source_q = [ComplexSum::default(); 3];
            let mut source_dq = [ComplexSum::default(); 3];
            let mut parameter_q = [ComplexSum::default(); 3];
            let mut envelopes = [0.0; 3];
            let mut parameter_envelopes = [0.0; 3];
            for j in 0..n {
                for (order, array) in direction.arrays().into_iter().enumerate() {
                    let a =
                        self.coefficient(self.source.source[order][i * n + j], i, j, order, 0)?;
                    let da = self.coefficient(array[i * n + j], i, j, order, exponent)?;
                    source_q[order].add(q[j].scale(a));
                    source_dq[order].add(dq[j].scale(a));
                    parameter_q[order].add(q[j].scale(da));
                    envelopes[order] += a.abs();
                    parameter_envelopes[order] += da.abs();
                }
            }
            // Independent split ORIGINAL source contractions, not the cached
            // numerical T or its bordered rows. This includes ds*T_s*phi.
            let residual = s * (s * source_dq[0].value() + source_dq[1].value())
                + source_dq[2].value()
                + (s * source_q[0].value().scale(2.0) + source_q[1].value()) * ds
                + s * (s * parameter_q[0].value() + parameter_q[1].value())
                + parameter_q[2].value();
            let polynomial = (envelopes[0] * radius + envelopes[1]) * radius + envelopes[2];
            let spectral = 2.0 * radius * envelopes[0] + envelopes[1];
            let parameter = (parameter_envelopes[0] * radius + parameter_envelopes[1]) * radius
                + parameter_envelopes[2];
            let bound = polynomial * dq_norm + (spectral * ds.abs() + parameter) * q_norm;
            worst = worst.max(relative(residual.abs(), bound)?);
        }
        let image = mass_image(self.source, &q, cancelled)?;
        let d_image = mass_image(self.source, &dq, cancelled)?;
        let mut motion = RealSum::default();
        for (&a, &b) in image.iter().zip(&d_image) {
            motion.add(product(a.conj(), b, 1)?.re);
        }
        let (mass, mass_bound) = self.mass_direction_form(direction, exponent, &q, cancelled)?;
        let bound = 2.0 * norm(&image) * norm(&d_image) + mass_bound;
        let normalization = relative((motion.value() + mass).abs(), bound)?;
        Ok((worst, normalization))
    }
}

#[derive(Clone, Copy, Default)]
struct RealSum {
    sum: f64,
    correction: f64,
}
impl RealSum {
    fn add(&mut self, term: f64) {
        let adjusted = term - self.correction;
        let next = self.sum + adjusted;
        self.correction = (next - self.sum) - adjusted;
        self.sum = next;
    }
    fn value(self) -> f64 {
        self.sum
    }
}
#[derive(Clone, Copy, Default)]
struct ComplexSum {
    re: RealSum,
    im: RealSum,
}
impl ComplexSum {
    fn add(&mut self, z: C64) {
        self.re.add(z.re);
        self.im.add(z.im);
    }
    fn value(self) -> C64 {
        C64::new(self.re.value(), self.im.value())
    }
}
fn norm(values: &[C64]) -> f64 {
    values.iter().fold(0.0, |n, v| det::hypot(n, v.abs()))
}
fn mass_image(
    source: &PreparedQuadraticModes<'_>,
    vector: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<C64>, QuadraticResponseError> {
    let mut image = Vec::with_capacity(source.n);
    for i in 0..source.n {
        poll(cancelled)?;
        let mut sum = ComplexSum::default();
        for (j, &v) in vector.iter().enumerate().skip(i) {
            sum.add(v.scale(source.mass_factor.l(j, i)));
        }
        let value = sum.value();
        if !finite(value) {
            return Err(invalid("physical mass image overflows"));
        }
        image.push(value);
    }
    Ok(image)
}
fn relative(residual: f64, bound: f64) -> Result<f64, QuadraticResponseError> {
    if !residual.is_finite() || !bound.is_finite() || bound < 0.0 {
        return Err(invalid(
            "original sensitivity residual is not representable",
        ));
    }
    if bound == 0.0 {
        if residual == 0.0 {
            Ok(0.0)
        } else {
            Err(invalid(
                "nonzero sensitivity residual without an equation scale",
            ))
        }
    } else {
        Ok(residual / bound)
    }
}
fn weighted_product(
    a: C64,
    b: C64,
    weight: f64,
    exponent: i32,
) -> Result<C64, QuadraticResponseError> {
    let (Some(ea), Some(eb)) = (component_exponent(a), component_exponent(b)) else {
        return Ok(C64::ZERO);
    };
    if weight == 0.0 {
        return Ok(C64::ZERO);
    }
    let (fraction, ew) = quadratic_fraction(weight);
    let a = scale_complex(a, -ea)?;
    let b = scale_complex(b, -eb)?;
    Ok(scale_complex(
        (a * b).scale(fraction.copysign(weight)),
        exponent + ea + eb + ew,
    )?)
}
fn sum_products(a: C64, b: C64, c: C64, d: C64) -> Result<C64, QuadraticResponseError> {
    let first = component_exponent(a)
        .zip(component_exponent(b))
        .map(|(x, y)| x + y);
    let second = component_exponent(c)
        .zip(component_exponent(d))
        .map(|(x, y)| x + y);
    let Some(common) = first.into_iter().chain(second).max() else {
        return Ok(C64::ZERO);
    };
    Ok(scale_complex(
        product(a, b, -common)? + product(c, d, -common)?,
        common,
    )?)
}
