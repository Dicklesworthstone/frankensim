//! Rotor resonance with original generalized Maxwell bearing memory.
//!
//! M q'' + (C + Omega G) q' + K_relaxed q + sum(b_j sigma_j) = f,
//! sigma_j' + sigma_j/tau_j = k_j b_j^T q'.
//! Positive branch stiffness/time describe passive hereditary strain channels;
//! the original bearing C/K may remain nonsymmetric or active. The unrelaxed
//! branch stiffness MUST NOT also be included in the supplied K.
//!
//! Harmonics use exp(+i omega t). Exact branch elimination keeps the solve at
//! the mechanical dimension, while every stress, heat rate and constitutive
//! residual is recovered. The full-state spectrum owner retains all 2n+b
//! poles, including pure material modes. No loss-factor fit or sampled damping
//! ratio replaces the original memory law.

use super::resonance::{RotorHarmonicDirection, RotorObservation};
use super::{
    MAX_ROTOR_DOFS, ROTOR_QR_SWEEPS_PER_WINDOW, RotorBearingModel, RotorError,
    RotorHarmonicResponse, RotorOperatingPoint, action, finite, finite_complex, norm, poll, power,
};
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel};
use fs_math::{c64::C64, dd::Dd, det};
use fs_modal::relaxation_modes::maxwell_eigenvalues_with_cancel;
pub use fs_modal::relaxation_modes::{
    MaxwellDampingBranch, MaxwellModeError, PoleOptions, PreparedMaxwellModes,
};
use fs_phs::harmonic::resonance::search_resonance;
pub use fs_phs::harmonic::resonance::{DrivenResonance, ResonanceControls, ResonanceSearchError};

/// Native dense spectrum has at most 2*MAX_ROTOR_DOFS + this many states.
pub const MAX_ROTOR_MAXWELL_BRANCHES: usize = 64;

/// Original rotor operating point and borrowed, individually retained memories.
#[derive(Debug)]
pub struct RotorMaxwellOperatingPoint<'a> {
    point: RotorOperatingPoint<'a>,
    branches: &'a [MaxwellDampingBranch],
}

/// Complete physical Laplace spectrum, with every original material state.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorMaxwellPoleReport {
    /// Signed shaft speed used in C + Omega G [rad/s].
    pub spin_rad_s: f64,
    /// All 2n+b native poles [1/s], sorted by real then imaginary part.
    /// Positive growth is retained. Dense QR does not certify tiny decay rates.
    pub poles_per_s: Vec<C64>,
    /// Number of independently retained original branch states.
    pub memory_states: usize,
}

/// Original displacement/velocity, branch stress and separately measured work.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorMaxwellHarmonicResponse {
    /// Mechanical fields and ORIGINAL M/C/G/K work. Its power remainder is
    /// input minus those four terms, before subtracting material branch heat.
    /// Its equilibrium residual already includes the full branch forces.
    pub motion: RotorHarmonicResponse,
    /// Stress conjugate to each original scalar strain b_j^T q.
    pub branch_stress: Vec<C64>,
    /// Independent cycle-averaged stress work Re(sigma_j^H b_j^T v)/2 [W].
    pub branch_mean_work_w: Vec<f64>,
    /// Physical passive heat |sigma_j|^2/(2 k_j tau_j) [W], in source order.
    pub branch_mean_heat_w: Vec<f64>,
    /// Largest normalized residual of the original branch constitutive rows.
    pub relaxation_relative_residual: f64,
    /// Supplied work minus original mechanical work and all branch heat [W].
    pub mean_power_residual_w: f64,
}

/// Tangent through frequency, spin and source with the original memories fixed.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorMaxwellHarmonicDerivative {
    /// Original physical carrier and source-equation checks.
    pub response: RotorMaxwellHarmonicResponse,
    /// Original displacement derivative per caller parameter.
    pub displacement_derivative: Vec<C64>,
    /// Complete derivative of i*omega*q, including the moving frequency.
    pub velocity_derivative: Vec<C64>,
    /// Each original stress derivative, including its own relaxation dispersion.
    pub branch_stress_derivative: Vec<C64>,
    /// Differentiated original mechanical equilibrium residual.
    pub derivative_relative_residual: f64,
    /// Differentiated original constitutive residual.
    pub relaxation_derivative_relative_residual: f64,
}

impl RotorBearingModel {
    /// Attach caller-supplied physical strain channels at a signed shaft speed.
    /// The existing stiffness is RELAXED; no branch is merged or fitted.
    /// A zero strain channel is retained because its free material pole exists.
    /// Branch coefficients are fixed at this supplied operating state; no
    /// preload, temperature or speed-dependent material law is inferred.
    pub fn at_speed_with_maxwell_memory<'a>(
        &'a self,
        spin_rad_s: f64,
        branches: &'a [MaxwellDampingBranch],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorMaxwellOperatingPoint<'a>, RotorError> {
        poll(&mut cancelled)?;
        if branches.len() > MAX_ROTOR_MAXWELL_BRANCHES {
            return Err(invalid(
                "Maxwell branch count exceeds the native rotor bound",
            ));
        }
        for branch in branches {
            poll(&mut cancelled)?;
            if branch.strain.len() != self.n
                || branch.strain.iter().any(|value| !value.is_finite())
                || !branch.stiffness.is_finite()
                || branch.stiffness <= 0.0
                || !branch.relaxation_time_s.is_finite()
                || branch.relaxation_time_s <= 0.0
            {
                return Err(invalid(
                    "complete finite strain channels and positive stiffness/time required",
                ));
            }
        }
        let point = self.at_speed(spin_rad_s, &mut cancelled)?;
        Ok(RotorMaxwellOperatingPoint { point, branches })
    }
}

impl RotorMaxwellOperatingPoint<'_> {
    /// Signed frozen shaft speed [rad/s].
    pub fn spin_rad_s(&self) -> f64 {
        self.point.spin_rad_s
    }

    /// All mechanical AND memory poles from the existing scaled full-state owner.
    /// Repeated equal-rate branches and mechanically dark states are retained.
    /// Native QR has bounded work and no relative tiny-decay-rate certificate.
    pub fn poles(
        &self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorMaxwellPoleReport, RotorError> {
        let model = self.point.model;
        let poles_per_s = maxwell_eigenvalues_with_cancel(
            &model.mass,
            &self.point.velocity_operator,
            &model.stiffness,
            model.n,
            self.branches,
            ROTOR_QR_SWEEPS_PER_WINDOW,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(RotorMaxwellPoleReport {
            spin_rad_s: self.point.spin_rad_s,
            poles_per_s,
            memory_states: self.branches.len(),
        })
    }

    /// Prepare physical local mechanical mode shapes, stresses and residues.
    /// The native owner preserves full C/K and all original branch laws.
    /// Its displacement-based rational domain excludes s=-1/tau and pure
    /// material modes; poles() above still retains their full-state values.
    pub fn prepare_modes(
        &self,
        options: PoleOptions,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMaxwellModes<'_>, MaxwellModeError> {
        let model = self.point.model;
        PreparedMaxwellModes::new(
            &model.mass,
            &self.point.velocity_operator,
            &model.stiffness,
            model.n,
            self.branches,
            options,
            cancelled,
        )
    }

    /// Particular +i harmonic response of the complete original memory law.
    /// Costs O(n^3+b*n^2), factors only n mechanical unknowns and recovers every
    /// branch. Active baseline bearings do not imply decaying transients.
    /// Singular, nonrepresentable and cancelled queries publish no partial state.
    pub fn harmonic_response(
        &self,
        frequency_rad_s: f64,
        force: &[C64],
        relative_tolerance: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<RotorMaxwellHarmonicResponse, RotorError> {
        self.factored_response(frequency_rad_s, force, 0.0, relative_tolerance, cancelled)
            .map(|(response, _, _, _)| response)
    }

    /// Analytic full source/frequency/spin tangent using the SAME native LU.
    /// Every branch contributes dH/domega = i*k*tau/(1+i*omega*tau)^2.
    /// Original branch constitutive rows and mechanical equations are audited
    /// after recovery. No finite difference or frozen frequency-sampled loss is used.
    pub fn harmonic_response_with_derivative(
        &self,
        frequency_rad_s: f64,
        force: &[C64],
        direction: RotorHarmonicDirection<'_>,
        relative_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorMaxwellHarmonicDerivative, RotorError> {
        poll(&mut cancelled)?;
        let model = self.point.model;
        if direction.force_derivative.len() != model.n
            || !finite_complex(direction.force_derivative)
            || !direction.frequency_derivative.is_finite()
            || !direction.spin_derivative.is_finite()
        {
            return Err(invalid(
                "finite matching source and frequency/spin derivatives required",
            ));
        }
        let w = frequency_rad_s;
        let dw = direction.frequency_derivative;
        let (response, factor, scales, weights) =
            self.factored_response(w, force, dw, relative_tolerance, &mut cancelled)?;
        let q = &response.motion.displacement;
        let mut matrix_slope = action(&model.mass, q, &mut cancelled)?;
        let velocity_slope = action(&self.point.velocity_operator, q, &mut cancelled)?;
        let gyro_slope = action(&model.gyro, q, &mut cancelled)?;
        for i in 0..model.n {
            matrix_slope[i] = sum([
                matrix_slope[i].scale(-w).scale(dw).scale(2.0),
                times_i(velocity_slope[i].scale(dw)),
                times_i(gyro_slope[i].scale(w).scale(direction.spin_derivative)),
            ])?;
        }
        for (branch, weight) in self.branches.iter().zip(&weights) {
            poll(&mut cancelled)?;
            let stress_slope = weight.slope * strain(&branch.strain, q)?;
            for (value, &b) in matrix_slope.iter_mut().zip(&branch.strain) {
                *value = sum([*value, stress_slope.scale(b)])?;
            }
        }
        let mut dq: Vec<_> = direction
            .force_derivative
            .iter()
            .zip(&matrix_slope)
            .zip(&scales[..model.n])
            .map(|((&force, &slope), &scale)| divided(force - slope, scale))
            .collect();
        if !finite_complex(&dq) {
            return Err(invalid("scaled Maxwell tangent source is unrepresentable"));
        }
        factor.solve_with_cancel(&mut dq, &mut cancelled)?;
        let dv: Vec<_> = dq
            .iter()
            .zip(q)
            .map(|(&dq, &q)| sum([dq.scale(w), q.scale(dw)]).map(times_i))
            .collect::<Result<_, _>>()?;
        let mut stresses = Vec::with_capacity(self.branches.len());
        let mut material = vec![C64::ZERO; model.n];
        let mut relaxation_residual: f64 = 0.0;
        for ((branch, weight), &stress) in self
            .branches
            .iter()
            .zip(&weights)
            .zip(&response.branch_stress)
        {
            poll(&mut cancelled)?;
            let e = strain(&branch.strain, q)?;
            let de = strain(&branch.strain, &dq)?;
            let ds = sum([weight.transfer * de, weight.slope * e])?;
            let terms = [
                ds.scale(weight.row_real),
                times_i(ds.scale(weight.row_imaginary)),
                -times_i(de.scale(branch.stiffness).scale(weight.row_imaginary)),
                times_i((stress - e.scale(branch.stiffness)).scale(weight.row_slope)),
            ];
            relaxation_residual = relaxation_residual.max(relative_row(&terms)?);
            for (value, &b) in material.iter_mut().zip(&branch.strain) {
                *value = sum([*value, ds.scale(b)])?;
            }
            stresses.push(ds);
        }
        let mut inertia = action(&model.mass, &dq, &mut cancelled)?;
        let moving_inertia = action(&model.mass, q, &mut cancelled)?;
        for i in 0..model.n {
            inertia[i] = sum([
                inertia[i].scale(-w).scale(w),
                moving_inertia[i].scale(-w).scale(dw).scale(2.0),
            ])?;
        }
        let damping = action(&model.damping, &dv, &mut cancelled)?;
        let mut gyro = action(&model.gyro, &dv, &mut cancelled)?;
        let moving_gyro = action(&model.gyro, &response.motion.velocity, &mut cancelled)?;
        for i in 0..model.n {
            gyro[i] = sum([
                gyro[i].scale(self.point.spin_rad_s),
                moving_gyro[i].scale(direction.spin_derivative),
            ])?;
        }
        let stiffness = action(&model.stiffness, &dq, &mut cancelled)?;
        let residual = equilibrium_residual(
            [&inertia, &damping, &gyro, &stiffness, &material],
            direction.force_derivative,
        )?;
        admit_residual(residual.max(relaxation_residual), relative_tolerance)?;
        poll(&mut cancelled)?;
        Ok(RotorMaxwellHarmonicDerivative {
            response,
            displacement_derivative: dq,
            velocity_derivative: dv,
            branch_stress_derivative: stresses,
            derivative_relative_residual: residual,
            relaxation_derivative_relative_residual: relaxation_residual,
        })
    }

    /// Actual fixed-spin, fixed-source driven peak with full material dispersion.
    /// The selected-amplitude half-height width is not inferred from a free pole.
    /// The caller brackets one isolated maximum; shared query budgets apply.
    pub fn driven_resonance(
        &self,
        interval_rad_s: (f64, f64),
        force: &[C64],
        observation: RotorObservation<'_>,
        relative_tolerance: f64,
        controls: ResonanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenResonance, ResonanceSearchError<RotorError>> {
        poll(&mut cancelled).map_err(ResonanceSearchError::Source)?;
        observation
            .validate(self.point.model.n)
            .map_err(ResonanceSearchError::Source)?;
        let zero = [C64::ZERO; MAX_ROTOR_DOFS];
        search_resonance(interval_rad_s, controls, cancelled, |w, _, cancel| {
            let result = self.harmonic_response_with_derivative(
                w,
                force,
                RotorHarmonicDirection {
                    frequency_derivative: 1.0,
                    spin_derivative: 0.0,
                    force_derivative: &zero[..self.point.model.n],
                },
                relative_tolerance,
                cancel,
            )?;
            let (weights, values, slopes) = match observation {
                RotorObservation::Displacement(weights) => (
                    weights,
                    &result.response.motion.displacement,
                    &result.displacement_derivative,
                ),
                RotorObservation::Velocity(weights) => (
                    weights,
                    &result.response.motion.velocity,
                    &result.velocity_derivative,
                ),
            };
            Ok((
                strain(weights, values)?,
                strain(weights, slopes)?,
                result
                    .response
                    .motion
                    .relative_residual
                    .max(result.response.relaxation_relative_residual)
                    .max(result.derivative_relative_residual)
                    .max(result.relaxation_derivative_relative_residual),
            ))
        })
    }

    fn factored_response(
        &self,
        w: f64,
        force: &[C64],
        dw: f64,
        tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        (
            RotorMaxwellHarmonicResponse,
            LuComplex,
            [f64; MAX_ROTOR_DOFS],
            Vec<BranchWeights>,
        ),
        RotorError,
    > {
        poll(&mut cancelled)?;
        let model = self.point.model;
        let n = model.n;
        if !w.is_finite()
            || force.len() != n
            || !finite_complex(force)
            || !tolerance.is_finite()
            || tolerance <= 0.0
            || tolerance >= 1.0
        {
            return Err(invalid(
                "finite frequency/force and relative tolerance in (0,1) required",
            ));
        }
        let mut weights = Vec::with_capacity(self.branches.len());
        for branch in self.branches {
            poll(&mut cancelled)?;
            weights.push(branch_weights(branch, w, dw)?);
        }
        let mut matrix = vec![C64::ZERO; n * n];
        let mut scales = [0.0; MAX_ROTOR_DOFS];
        let mut q = Vec::with_capacity(n);
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut row_scale: f64 = 0.0;
            for j in 0..n {
                let ij = i * n + j;
                // Match the native rotor's compensated inertial detuning.
                let product = finite(w * model.mass[ij])?;
                let error = w.mul_add(model.mass[ij], -product);
                let mut value = ComplexSum::default();
                value.add(C64::new(
                    finite((-w).mul_add(product, model.stiffness[ij]) - w * error)?,
                    finite(w * self.point.velocity_operator[ij])?,
                ));
                for (branch, weight) in self.branches.iter().zip(&weights) {
                    value.add(
                        weight
                            .transfer
                            .scale(branch.strain[i])
                            .scale(branch.strain[j]),
                    );
                }
                matrix[ij] = value.value()?;
                row_scale = row_scale.max(finite(matrix[ij].abs())?);
            }
            if row_scale == 0.0 {
                return Err(LuComplexError::Singular { index: i }.into());
            }
            scales[i] = row_scale;
            for value in &mut matrix[i * n..(i + 1) * n] {
                *value = divided(*value, row_scale);
            }
            q.push(divided(force[i], row_scale));
        }
        if !finite_complex(&q) {
            return Err(invalid("scaled Maxwell forcing is unrepresentable"));
        }
        let factor = lu_complex_with_cancel(&matrix, n, &mut cancelled)?;
        factor.solve_with_cancel(&mut q, &mut cancelled)?;
        let v: Vec<_> = q.iter().map(|&q| times_i(q.scale(w))).collect();
        if !finite_complex(&q) || !finite_complex(&v) {
            return Err(invalid("Maxwell harmonic motion is unrepresentable"));
        }
        let mut stresses = Vec::with_capacity(self.branches.len());
        let mut work = Vec::with_capacity(self.branches.len());
        let mut heat = Vec::with_capacity(self.branches.len());
        let mut material = vec![C64::ZERO; n];
        let mut relaxation_residual: f64 = 0.0;
        for (branch, weight) in self.branches.iter().zip(&weights) {
            poll(&mut cancelled)?;
            let e = strain(&branch.strain, &q)?;
            let stress = weight.transfer * e;
            let terms = [
                stress.scale(weight.row_real),
                times_i(stress.scale(weight.row_imaginary)),
                -times_i(e.scale(branch.stiffness).scale(weight.row_imaginary)),
            ];
            relaxation_residual = relaxation_residual.max(relative_row(&terms)?);
            let strain_velocity = strain(&branch.strain, &v)?;
            work.push(finite(0.5 * (stress.conj() * strain_velocity).re)?);
            let amplitude = finite(
                stress.abs() / (det::sqrt(branch.stiffness) * det::sqrt(branch.relaxation_time_s)),
            )?;
            heat.push(finite((0.5 * amplitude) * amplitude)?);
            for (value, &b) in material.iter_mut().zip(&branch.strain) {
                *value = sum([*value, stress.scale(b)])?;
            }
            stresses.push(stress);
        }
        let mut inertia = action(&model.mass, &q, &mut cancelled)?;
        for value in &mut inertia {
            *value = value.scale(-w).scale(w);
        }
        let damping = action(&model.damping, &v, &mut cancelled)?;
        let mut gyro = action(&model.gyro, &v, &mut cancelled)?;
        for value in &mut gyro {
            *value = value.scale(self.point.spin_rad_s);
        }
        let stiffness = action(&model.stiffness, &q, &mut cancelled)?;
        let residual =
            equilibrium_residual([&inertia, &damping, &gyro, &stiffness, &material], force)?;
        admit_residual(residual.max(relaxation_residual), tolerance)?;
        let input_power = power(force, &v)?;
        let damping_power = power(&damping, &v)?;
        let stiffness_power = power(&stiffness, &v)?;
        let gyro_power = power(&gyro, &v)?;
        let inertia_power = power(&inertia, &v)?;
        let mechanical_remainder = sum([
            C64::from_re(input_power),
            C64::from_re(-damping_power),
            C64::from_re(-stiffness_power),
            C64::from_re(-gyro_power),
            C64::from_re(-inertia_power),
        ])?
        .re;
        let mut balance = ComplexSum::default();
        balance.add(C64::from_re(mechanical_remainder));
        for &value in &heat {
            balance.add(C64::from_re(-value));
        }
        let response = RotorMaxwellHarmonicResponse {
            motion: RotorHarmonicResponse {
                spin_rad_s: self.point.spin_rad_s,
                frequency_rad_s: w,
                displacement: q,
                velocity: v,
                relative_residual: residual,
                mean_input_power_w: input_power,
                mean_velocity_force_power_w: damping_power,
                mean_stiffness_force_power_w: stiffness_power,
                mean_gyroscopic_power_w: gyro_power,
                mean_inertial_power_w: inertia_power,
                mean_power_residual_w: mechanical_remainder,
            },
            branch_stress: stresses,
            branch_mean_work_w: work,
            branch_mean_heat_w: heat,
            relaxation_relative_residual: relaxation_residual,
            mean_power_residual_w: balance.value()?.re,
        };
        poll(&mut cancelled)?;
        Ok((response, factor, scales, weights))
    }
}

#[derive(Clone, Copy, Debug)]
struct BranchWeights {
    transfer: C64,
    slope: C64,
    row_real: f64,
    row_imaginary: f64,
    row_slope: f64,
}

fn branch_weights(
    branch: &MaxwellDampingBranch,
    w: f64,
    dw: f64,
) -> Result<BranchWeights, RotorError> {
    let k = branch.stiffness;
    let tau = branch.relaxation_time_s;
    let x = w * tau;
    let (transfer, slope, row_real, row_imaginary, row_slope) = if x.abs() <= 1.0 {
        let denominator = 1.0 + x * x;
        let imaginary = (k * x) / denominator;
        let transfer = C64::new(imaginary * x, imaginary);
        let slope = if dw == 0.0 {
            C64::ZERO
        } else {
            C64::new(
                (2.0 * x / denominator) / denominator,
                ((1.0 - x * x) / denominator) / denominator,
            )
            .scale(k)
            .scale(tau)
            .scale(dw)
        };
        (
            transfer,
            slope,
            1.0,
            x,
            if dw == 0.0 { 0.0 } else { tau * dw },
        )
    } else {
        // Reciprocal form avoids squaring or even forming a large omega*tau.
        let reciprocal = (1.0 / w) / tau;
        let denominator = 1.0 + reciprocal * reciprocal;
        let real = k / denominator;
        let transfer = C64::new(real, real * reciprocal);
        let slope = if dw == 0.0 {
            C64::ZERO
        } else {
            divided(
                C64::new(
                    (2.0 * reciprocal * reciprocal / denominator) / denominator,
                    (-reciprocal * (1.0 - reciprocal * reciprocal) / denominator) / denominator,
                )
                .scale(k),
                w,
            )
            .scale(dw)
        };
        (
            transfer,
            slope,
            reciprocal,
            1.0,
            if dw == 0.0 { 0.0 } else { dw / w },
        )
    };
    if !finite_complex(&[transfer, slope]) || !row_slope.is_finite() {
        return Err(invalid(
            "Maxwell transfer or directional derivative is unrepresentable",
        ));
    }
    if w != 0.0 && transfer.im == 0.0 && branch.strain.iter().any(|&value| value != 0.0) {
        // A positive finite branch always has nonzero signed loss away from
        // DC. Dropping it can erase representable work after multiplication
        // by frequency, even while the separately recovered heat remains finite.
        return Err(invalid("Maxwell loss quadrature underflowed"));
    }
    Ok(BranchWeights {
        transfer,
        slope,
        row_real,
        row_imaginary,
        row_slope,
    })
}

fn equilibrium_residual(terms: [&[C64]; 5], force: &[C64]) -> Result<f64, RotorError> {
    let mut residual = Vec::with_capacity(force.len());
    for (i, &f) in force.iter().enumerate() {
        residual.push(sum([
            terms[0][i],
            terms[1][i],
            terms[2][i],
            terms[3][i],
            terms[4][i],
            -f,
        ])?);
    }
    let mut norms = [0.0; 6];
    for (target, source) in norms.iter_mut().zip(terms.into_iter().chain([force])) {
        *target = norm(source)?;
    }
    let scale = norms.into_iter().fold(0.0_f64, f64::max);
    if scale == 0.0 {
        norm(&residual)
    } else {
        finite((norm(&residual)? / scale) / norms.into_iter().map(|v| v / scale).sum::<f64>())
    }
}
fn relative_row(terms: &[C64]) -> Result<f64, RotorError> {
    let scale = terms.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
    finite(scale)?;
    if scale == 0.0 {
        return Ok(0.0);
    }
    let mut residual = ComplexSum::default();
    let mut denominator = 0.0;
    for &term in terms {
        let value = divided(term, scale);
        residual.add(value);
        denominator += value.abs();
    }
    finite(residual.value()?.abs() / denominator)
}
fn admit_residual(measured: f64, limit: f64) -> Result<(), RotorError> {
    if measured > limit {
        Err(RotorError::Residual {
            relative: measured,
            limit,
        })
    } else {
        Ok(())
    }
}
fn invalid(what: &'static str) -> RotorError {
    RotorError::Invalid(what)
}
fn divided(value: C64, scale: f64) -> C64 {
    C64::new(value.re / scale, value.im / scale)
}
fn times_i(value: C64) -> C64 {
    C64::new(-value.im, value.re)
}
#[derive(Default)]
struct ComplexSum {
    real: Dd,
    imaginary: Dd,
}
impl ComplexSum {
    fn add(&mut self, value: C64) {
        self.real = self.real + Dd::from_f64(value.re);
        self.imaginary = self.imaginary + Dd::from_f64(value.im);
    }
    fn value(self) -> Result<C64, RotorError> {
        Ok(C64::new(
            finite(self.real.to_f64())?,
            finite(self.imaginary.to_f64())?,
        ))
    }
}
fn sum<const N: usize>(values: [C64; N]) -> Result<C64, RotorError> {
    let mut result = ComplexSum::default();
    for value in values {
        result.add(value);
    }
    result.value()
}
fn strain(direction: &[f64], q: &[C64]) -> Result<C64, RotorError> {
    let mut real = Dd::ZERO;
    let mut imaginary = Dd::ZERO;
    for (&b, &q) in direction.iter().zip(q) {
        real = real + Dd::from_f64(b) * Dd::from_f64(q.re);
        imaginary = imaginary + Dd::from_f64(b) * Dd::from_f64(q.im);
    }
    Ok(C64::new(
        finite(real.to_f64())?,
        finite(imaginary.to_f64())?,
    ))
}
