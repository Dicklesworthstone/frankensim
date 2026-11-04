//! Exact conservative motion in a complete, admitted generalized eigenbasis.
//!
//! All n modes of the original n-coordinate pencil are mandatory. This is a
//! change of coordinates, not modal reduction. The caller supplies the basis;
//! no unpolled eigensolve is hidden in preparation. Fixed-clock oscillator
//! coefficients replace the large augmented exponential. Preparation costs
//! O(n^3), retained storage and warm intervals O(n^2), with explicit capacities.
//! Force work is integrated in original coordinates independently of energy.

use super::{
    ExactLinearSolidEnergy, ExactLinearSolidStep, LinearSolidImpulseError,
    LinearSolidImpulseReport, impulse,
};
use crate::linear3::TetElasticAssembly;
use fs_exec::Cx;
use fs_math::det;
use fs_sparse::Csr;

/// Exact cubic forces and continuous displacement moments for moving supports.
pub mod polynomial;
mod sustained;

/// Required capacities and numerical allowances; no implicit defaults.
#[derive(Clone, Copy, Debug)]
pub struct CompleteSpectralSolidControls {
    pub maximum_free_dofs: usize,
    /// Conservative simultaneous retained/scratch byte envelope.
    pub maximum_preparation_bytes: usize,
    /// Admission uses 4*n^3 + 32*n^2 + 512*n scalar product terms.
    pub maximum_preparation_multiply_terms: u64,
    /// Each interval admits 16*n^2 + 512*n scalar product terms.
    pub maximum_step_multiply_terms: u64,
    /// Relative matrix symmetry, eigen residual, and infinity norm of
    /// Phi^T M Phi - I. Must be strictly less than one.
    pub maximum_relative_operator_error: f64,
    /// Absolute allowance for K*phi - lambda*M*phi [sqrt(kg)/s^2].
    /// Explicitly controls roundoff in rigid eigenvectors; zero is permitted.
    /// No eigenvalue is rounded, clipped, or replaced by this allowance.
    pub maximum_absolute_eigen_residual: f64,
    /// Initial u/v reconstruction error relative to each field's own norm.
    /// The zero field must reconstruct exactly; no unit-valued floor is used.
    pub maximum_relative_state_error: f64,
    pub maximum_work_error_j: f64,
    pub maximum_steps: u64,
}

#[derive(Debug)]
pub enum CompleteSpectralSolidError {
    Invalid(&'static str),
    Budget { required: u128, maximum: u128 },
    Operator { defect: f64, allowance: f64 },
    Energy { defect_j: f64, allowance_j: f64 },
    Cancelled,
}
impl core::fmt::Display for CompleteSpectralSolidError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "complete spectral solid refused: {self:?}")
    }
}
impl core::error::Error for CompleteSpectralSolidError {}
type Result<T> = core::result::Result<T, CompleteSpectralSolidError>;
fn invalid(what: &'static str) -> CompleteSpectralSolidError {
    CompleteSpectralSolidError::Invalid(what)
}
fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint()
        .map_err(|_| CompleteSpectralSolidError::Cancelled)
}
fn bound(required: u128, maximum: u128) -> Result<()> {
    if required > maximum {
        Err(CompleteSpectralSolidError::Budget { required, maximum })
    } else {
        Ok(())
    }
}
fn zeros(n: usize) -> Result<Vec<f64>> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| invalid("allocation refused"))?;
    values.resize(n, 0.0);
    Ok(values)
}
fn copy(values: &[f64]) -> Result<Vec<f64>> {
    let mut result = zeros(values.len())?;
    result.copy_from_slice(values);
    Ok(result)
}
fn add(sum: &mut f64, correction: &mut f64, value: f64) {
    let increment = value - *correction;
    let next = *sum + increment;
    *correction = (next - *sum) - increment;
    *sum = next;
}

#[derive(Clone, Copy)]
struct Oscillator {
    cosine: f64,
    h_sinc: f64,
    h2_cosc: f64,
    h2_sinc3: f64,
    minus_h_lambda_sinc: f64,
    h_cosc: f64,
    sinc: f64,
    h2_cosc4: f64,
}
impl Oscillator {
    fn prepare(lambda: f64, h: f64) -> Result<Self> {
        let z = lambda.sqrt() * h;
        let z2 = z * z;
        if !z2.is_finite() {
            return Err(invalid("unrepresentable oscillator phase"));
        }
        let (sinc, cosc, sinc3, cosc4, cosine) = if z.abs() <= 0.5 {
            // Separate entire-function series avoid subtraction of nearly
            // equal values. At z=0 these are exact rigid-mode coefficients.
            let series = |offset: usize| {
                let mut factorial = 1.0;
                for k in 2..=offset {
                    factorial *= k as f64;
                }
                let mut term = 1.0 / factorial;
                let mut sum = term;
                for k in 1..=9 {
                    let a = (2 * k + offset - 1) as f64;
                    let b = (2 * k + offset) as f64;
                    term *= -z2 / (a * b);
                    sum += term;
                }
                sum
            };
            let cosc = series(2);
            (
                series(1),
                cosc,
                series(3),
                series(4),
                (-z2).mul_add(cosc, 1.0),
            )
        } else {
            let (sine, cosine) = det::sin_cos(z);
            let sinc = sine / z;
            let cosc = (1.0 - cosine) / z2;
            (sinc, cosc, (1.0 - sinc) / z2, (0.5 - cosc) / z2, cosine)
        };
        let result = Self {
            cosine,
            h_sinc: h * sinc,
            h2_cosc: h * (h * cosc),
            h2_sinc3: h * (h * sinc3),
            minus_h_lambda_sinc: -(h * lambda) * sinc,
            h_cosc: h * cosc,
            sinc,
            h2_cosc4: h * (h * cosc4),
        };
        if [
            result.cosine,
            result.h_sinc,
            result.h2_cosc,
            result.h2_sinc3,
            result.minus_h_lambda_sinc,
            result.h_cosc,
            result.sinc,
            result.h2_cosc4,
        ]
        .iter()
        .any(|x| !x.is_finite())
        {
            return Err(invalid("unrepresentable oscillator coefficients"));
        }
        Ok(result)
    }
}

/// Full-coordinate conservative solid with every supplied basis column kept.
/// The basis is approximate only within its explicitly admitted residuals.
/// Original K/M energy checks remain active on every accepted interval.
/// Supports are homogeneous/stationary, as in ExactLinearSolid. There is no
/// constitutive damping, nonlinear finite-rotation, or sparse-scale claim.
pub struct CompleteSpectralSolid {
    stiffness: Vec<f64>,
    mass: Vec<f64>,
    basis: Vec<f64>,
    free_dofs: Vec<usize>,
    oscillators: Vec<Oscillator>,
    polynomial_oscillators: Option<Vec<polynomial::PolynomialOscillator>>,
    displacement: Vec<f64>,
    velocity: Vec<f64>,
    modal_displacement: Vec<f64>,
    modal_velocity: Vec<f64>,
    dt_s: f64,
    controls: CompleteSpectralSolidControls,
    steps: u64,
    stored_j: f64,
    kinetic_j: f64,
    initial_energy_j: f64,
    cumulative_work_j: f64,
    work_correction_j: f64,
}

/// Unpublished native impulse response for a surrounding full-coordinate
/// support or coupling audit. Prepared matrices and oscillators are borrowed.
pub(crate) struct SpectralImpulseCandidate {
    pub(crate) velocity_mps: Vec<f64>,
    pub(crate) report: LinearSolidImpulseReport,
    modal_velocity: Vec<f64>,
    cumulative_work_j: f64,
    work_correction_j: f64,
}

impl CompleteSpectralSolid {
    /// Admit every generalized eigenpair against the original pencil.
    /// Basis storage is row-major, with eigenvector j in column j:
    /// basis_columns[i*n+j] = phi_j[i]. Ordering need not be spectral.
    /// Exactly n finite nonnegative eigenvalues and n*n entries are required.
    /// The row-sum norm of Phi^T M Phi - I must be below one. Consequently the
    /// admitted square basis is nonsingular and the symmetric mass is positive
    /// definite to the declared numerical allowance; a duplicate/missing mode
    /// cannot pass as a complete basis.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_with_basis(
        assembly: &TetElasticAssembly,
        eigenvalues_rad2_s2: &[f64],
        basis_columns: &[f64],
        dt_s: f64,
        initial_displacement_m: &[f64],
        initial_velocity_mps: &[f64],
        controls: CompleteSpectralSolidControls,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        let n = assembly.free_dofs.len();
        let square = n.checked_mul(n).ok_or_else(|| invalid("shape overflow"))?;
        let good_relative = |x: f64| x.is_finite() && (0.0..1.0).contains(&x);
        if n == 0
            || !dt_s.is_finite()
            || dt_s <= 0.0
            || controls.maximum_steps == 0
            || controls.maximum_steps > (1_u64 << 53)
            || !good_relative(controls.maximum_relative_operator_error)
            || !good_relative(controls.maximum_relative_state_error)
            || !controls.maximum_absolute_eigen_residual.is_finite()
            || controls.maximum_absolute_eigen_residual < 0.0
            || !controls.maximum_work_error_j.is_finite()
            || controls.maximum_work_error_j < 0.0
            || eigenvalues_rad2_s2.len() != n
            || basis_columns.len() != square
            || initial_displacement_m.len() != n
            || initial_velocity_mps.len() != n
        {
            return Err(invalid(
                "complete finite basis, state, clock and controls required",
            ));
        }
        bound(n as u128, controls.maximum_free_dofs as u128)?;
        let nn = n as u128;
        // This bounds copied matrices, simultaneous dual/physical/modal
        // scratch, coefficients, free map, and accepted/output vectors.
        bound(
            64 * nn * nn + 1024 * nn,
            controls.maximum_preparation_bytes as u128,
        )?;
        bound(
            4 * nn * nn * nn + 32 * nn * nn + 512 * nn,
            controls.maximum_preparation_multiply_terms as u128,
        )?;
        bound(
            16 * nn * nn + 512 * nn,
            controls.maximum_step_multiply_terms as u128,
        )?;
        for i in 0..n {
            poll(cx)?;
            if (i > 0 && assembly.free_dofs[i - 1] >= assembly.free_dofs[i])
                || !eigenvalues_rad2_s2[i].is_finite()
                || eigenvalues_rad2_s2[i] < 0.0
                || !initial_displacement_m[i].is_finite()
                || !initial_velocity_mps[i].is_finite()
                || basis_columns[i * n..(i + 1) * n]
                    .iter()
                    .any(|x| !x.is_finite())
            {
                return Err(invalid(
                    "finite complete basis, state and sorted free map required",
                ));
            }
        }
        let tolerance = controls.maximum_relative_operator_error;
        let stiffness = symmetric(&assembly.stiffness, n, tolerance, cx)?;
        let mass = symmetric(&assembly.mass, n, tolerance, cx)?;
        let mut basis = zeros(square)?;
        for i in 0..n {
            poll(cx)?;
            basis[i * n..(i + 1) * n].copy_from_slice(&basis_columns[i * n..(i + 1) * n]);
        }
        let mut dual = zeros(square)?;
        for i in 0..n {
            for j in 0..n {
                poll(cx)?;
                dual[i * n + j] =
                    (0..n).fold(0.0, |s, k| mass[i * n + k].mul_add(basis[k * n + j], s));
            }
        }
        for j in 0..n {
            let (mut residual, mut scale) = (0.0_f64, 0.0_f64);
            for i in 0..n {
                poll(cx)?;
                let kphi = (0..n).fold(0.0, |s, k| {
                    stiffness[i * n + k].mul_add(basis[k * n + j], s)
                });
                let lambda_mphi = eigenvalues_rad2_s2[j] * dual[i * n + j];
                residual = residual.max((kphi - lambda_mphi).abs());
                scale = scale.max(kphi.abs()).max(lambda_mphi.abs());
                if !kphi.is_finite() || !lambda_mphi.is_finite() {
                    return Err(invalid("unrepresentable eigen residual"));
                }
            }
            let allowance = controls.maximum_absolute_eigen_residual + tolerance * scale;
            check_operator(residual, allowance)?;
        }
        for i in 0..n {
            let (mut row_error, mut correction) = (0.0, 0.0);
            for j in 0..n {
                poll(cx)?;
                let gram = (0..n).fold(0.0, |s, k| basis[k * n + i].mul_add(dual[k * n + j], s));
                add(
                    &mut row_error,
                    &mut correction,
                    (gram - if i == j { 1.0 } else { 0.0 }).abs(),
                );
            }
            check_operator(row_error, tolerance)?;
        }
        let modal_displacement = project(&dual, initial_displacement_m, cx)?;
        let modal_velocity = project(&dual, initial_velocity_mps, cx)?;
        for (field, modal) in [
            (initial_displacement_m, &modal_displacement),
            (initial_velocity_mps, &modal_velocity),
        ] {
            let reconstructed = reconstruct(&basis, modal, cx)?;
            let scale = field.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
            let error = field
                .iter()
                .zip(reconstructed)
                .fold(0.0_f64, |s, (a, b)| s.max((a - b).abs()));
            check_operator(error, controls.maximum_relative_state_error * scale)?;
        }
        let stored_j = quadratic(&stiffness, initial_displacement_m, cx)?;
        let kinetic_j = quadratic(&mass, initial_velocity_mps, cx)?;
        if !(stored_j + kinetic_j).is_finite() {
            return Err(invalid("unrepresentable initial energy"));
        }
        let mut oscillators = Vec::new();
        oscillators
            .try_reserve_exact(n)
            .map_err(|_| invalid("coefficient allocation refused"))?;
        for &lambda in eigenvalues_rad2_s2 {
            poll(cx)?;
            oscillators.push(Oscillator::prepare(lambda, dt_s)?);
        }
        let mut free_dofs = Vec::new();
        free_dofs
            .try_reserve_exact(n)
            .map_err(|_| invalid("free map allocation refused"))?;
        free_dofs.extend_from_slice(&assembly.free_dofs);
        poll(cx)?;
        Ok(Self {
            stiffness,
            mass,
            basis,
            free_dofs,
            oscillators,
            polynomial_oscillators: None,
            displacement: copy(initial_displacement_m)?,
            velocity: copy(initial_velocity_mps)?,
            modal_displacement,
            modal_velocity,
            dt_s,
            controls,
            steps: 0,
            stored_j,
            kinetic_j,
            initial_energy_j: stored_j + kinetic_j,
            cumulative_work_j: 0.0,
            work_correction_j: 0.0,
        })
    }

    pub fn displacement_m(&self) -> &[f64] {
        &self.displacement
    }
    pub fn velocity_mps(&self) -> &[f64] {
        &self.velocity
    }
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }
    pub fn steps(&self) -> u64 {
        self.steps
    }
    pub fn time_s(&self) -> f64 {
        self.steps as f64 * self.dt_s
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.stored_j
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        self.kinetic_j
    }
    pub fn cumulative_work_j(&self) -> f64 {
        self.cumulative_work_j
    }

    /// Apply an instantaneous impulse in the original `free_dofs()` order.
    /// All complete-basis modes receive Phi^T J, including admitted rigid
    /// modes. No stiffness inverse, mode cutoff or new factorization is used.
    /// The original M audits the actual represented physical velocity change;
    /// an insufficiently accurate basis or an unrepresentable kick refuses.
    ///
    /// Displacement, clock, step count and prepared oscillator/source maps stay
    /// fixed. J dot average(v) enters the existing cumulative work exactly once.
    /// Zero impulses preserve both physical and modal state bits. This reduced
    /// pencil has no discarded support reactions or independent contact law.
    pub fn apply_momentum_impulse(
        &mut self,
        impulse_kg_m_s: &[f64],
        maximum_impulse_error_kg_m_s: f64,
        cx: &Cx<'_>,
    ) -> core::result::Result<LinearSolidImpulseReport, LinearSolidImpulseError> {
        let candidate =
            self.trial_momentum_impulse(impulse_kg_m_s, maximum_impulse_error_kg_m_s, cx)?;
        Ok(self.accept_momentum_impulse(candidate))
    }

    /// The same complete basis used by native impulses, without publishing a
    /// state or changing modal coordinates. Contact-space mobility uses every
    /// column, including rigid modes; the final impulse still receives the
    /// original-M and actual-velocity audit before acceptance.
    pub(crate) fn momentum_response(
        &self,
        impulse_kg_m_s: &[f64],
        cx: &Cx<'_>,
    ) -> core::result::Result<Vec<f64>, LinearSolidImpulseError> {
        impulse::admit(impulse_kg_m_s, self.free_dofs.len(), 0.0, cx)?;
        let map_error = |error| match error {
            CompleteSpectralSolidError::Cancelled => LinearSolidImpulseError::Cancelled,
            _ => {
                LinearSolidImpulseError::Invalid("unrepresentable complete-basis impulse response")
            }
        };
        let modal = project(&self.basis, impulse_kg_m_s, cx).map_err(map_error)?;
        reconstruct(&self.basis, &modal, cx).map_err(map_error)
    }

    pub(crate) fn trial_momentum_impulse(
        &self,
        impulse_kg_m_s: &[f64],
        maximum_impulse_error_kg_m_s: f64,
        cx: &Cx<'_>,
    ) -> core::result::Result<SpectralImpulseCandidate, LinearSolidImpulseError> {
        let n = self.free_dofs.len();
        impulse::admit(impulse_kg_m_s, n, maximum_impulse_error_kg_m_s, cx)?;
        let map_error = |error| match error {
            CompleteSpectralSolidError::Cancelled => LinearSolidImpulseError::Cancelled,
            _ => {
                LinearSolidImpulseError::Invalid("unrepresentable complete-basis impulse response")
            }
        };
        let mut modal_velocity = self.modal_velocity.clone();
        let zero = impulse_kg_m_s.iter().all(|j| *j == 0.0);
        let velocity = if zero {
            self.velocity.clone()
        } else {
            let response = project(&self.basis, impulse_kg_m_s, cx).map_err(map_error)?;
            for (value, change) in modal_velocity.iter_mut().zip(response) {
                if change != 0.0 {
                    *value += change;
                }
                if !value.is_finite() {
                    return Err(LinearSolidImpulseError::Invalid(
                        "unrepresentable complete-basis velocity",
                    ));
                }
            }
            reconstruct(&self.basis, &modal_velocity, cx).map_err(map_error)?
        };
        let kinetic_after_j = if zero {
            self.kinetic_j
        } else {
            quadratic(&self.mass, &velocity, cx).map_err(map_error)?
        };
        let candidate = impulse::audit(
            &self.mass,
            &self.velocity,
            &velocity,
            impulse_kg_m_s,
            maximum_impulse_error_kg_m_s,
            kinetic_after_j,
            impulse::History {
                time_s: self.time_s(),
                stored_j: self.stored_j,
                kinetic_before_j: self.kinetic_j,
                initial_energy_j: self.initial_energy_j,
                cumulative_work_j: self.cumulative_work_j,
                work_correction_j: self.work_correction_j,
                cumulative_dissipated_j: 0.0,
                maximum_work_error_j: self.controls.maximum_work_error_j,
            },
            cx,
        )?;
        Ok(SpectralImpulseCandidate {
            velocity_mps: velocity,
            report: candidate.report,
            modal_velocity,
            cumulative_work_j: candidate.cumulative_work_j,
            work_correction_j: candidate.work_correction_j,
        })
    }

    pub(crate) fn accept_momentum_impulse(
        &mut self,
        candidate: SpectralImpulseCandidate,
    ) -> LinearSolidImpulseReport {
        self.velocity = candidate.velocity_mps;
        self.modal_velocity = candidate.modal_velocity;
        self.kinetic_j = candidate.report.energy.kinetic_after_j;
        self.cumulative_work_j = candidate.cumulative_work_j;
        self.work_correction_j = candidate.work_correction_j;
        candidate.report
    }

    /// Exact held-force motion, in complete free-DOF order.
    pub fn advance(&mut self, force_n: &[f64], cx: &Cx<'_>) -> Result<ExactLinearSolidStep> {
        self.advance_affine(force_n, force_n, cx)
    }

    /// Exact linear interpolation of the two force endpoints. A jump between
    /// intervals changes only the source and creates no velocity impulse.
    /// Mean kinematic velocity is Delta u / h; the separately integrated work
    /// below is needed for a varying force. Accepted state changes only after
    /// original-coordinate interval and cumulative energy closure both pass.
    pub fn advance_affine(
        &mut self,
        force_start_n: &[f64],
        force_end_n: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ExactLinearSolidStep> {
        poll(cx)?;
        let n = self.free_dofs.len();
        if force_start_n.len() != n
            || force_end_n.len() != n
            || force_start_n
                .iter()
                .chain(force_end_n)
                .any(|x| !x.is_finite())
        {
            return Err(invalid("complete finite force endpoints required"));
        }
        let (step, time_s) = self.next_interval()?;
        let mut increment = zeros(n)?;
        for i in 0..n {
            increment[i] = force_end_n[i] - force_start_n[i];
            if !increment[i].is_finite() {
                return Err(invalid("force increment overflow"));
            }
        }
        let force = project(&self.basis, force_start_n, cx)?;
        let delta = project(&self.basis, &increment, cx)?;
        let (mut q, mut v, mut mean_q) = (zeros(n)?, zeros(n)?, zeros(n)?);
        for i in 0..n {
            poll(cx)?;
            let [position, velocity, mean] = self.affine_mode(i, force[i], delta[i]);
            q[i] = position;
            v[i] = velocity;
            mean_q[i] = mean;
        }
        let displacement = reconstruct(&self.basis, &q, cx)?;
        let velocity = reconstruct(&self.basis, &v, cx)?;
        let mean_displacement = reconstruct(&self.basis, &mean_q, cx)?;
        let (mut supplied_work_j, mut work_roundoff) = (0.0, 0.0);
        let mut mean_velocity_mps = zeros(n)?;
        for i in 0..n {
            poll(cx)?;
            let du = displacement[i] - self.displacement[i];
            // Integration by parts: integral f du = f0*Delta u +
            // Delta f*(u1 - mean(u)). No energy difference enters this work.
            add(
                &mut supplied_work_j,
                &mut work_roundoff,
                force_start_n[i]
                    .mul_add(du, increment[i] * (displacement[i] - mean_displacement[i])),
            );
            mean_velocity_mps[i] = du / self.dt_s;
        }
        if q.iter()
            .chain(&v)
            .chain(&mean_velocity_mps)
            .any(|x| !x.is_finite())
        {
            return Err(invalid("unrepresentable propagated state"));
        }
        let stored_after_j = quadratic(&self.stiffness, &displacement, cx)?;
        let kinetic_after_j = quadratic(&self.mass, &velocity, cx)?;
        let change = quadratic_change(&self.stiffness, &self.displacement, &displacement, cx)?
            + quadratic_change(&self.mass, &self.velocity, &velocity, cx)?;
        let (mut cumulative_work_j, mut correction) =
            (self.cumulative_work_j, self.work_correction_j);
        add(&mut cumulative_work_j, &mut correction, supplied_work_j);
        let energy = ExactLinearSolidEnergy {
            stored_before_j: self.stored_j,
            stored_after_j,
            kinetic_before_j: self.kinetic_j,
            kinetic_after_j,
            supplied_work_j,
            dissipated_j: 0.0,
            cumulative_dissipated_j: 0.0,
            balance_error_j: supplied_work_j - change,
            cumulative_balance_error_j: cumulative_work_j
                - (stored_after_j + kinetic_after_j - self.initial_energy_j),
        };
        for defect_j in [energy.balance_error_j, energy.cumulative_balance_error_j] {
            if !defect_j.is_finite() || defect_j.abs() > self.controls.maximum_work_error_j {
                return Err(CompleteSpectralSolidError::Energy {
                    defect_j,
                    allowance_j: self.controls.maximum_work_error_j,
                });
            }
        }
        if !cumulative_work_j.is_finite() || !correction.is_finite() {
            return Err(invalid("unrepresentable accumulated work"));
        }
        poll(cx)?;
        self.displacement = displacement;
        self.velocity = velocity;
        self.modal_displacement = q;
        self.modal_velocity = v;
        self.steps = step;
        self.stored_j = stored_after_j;
        self.kinetic_j = kinetic_after_j;
        self.cumulative_work_j = cumulative_work_j;
        self.work_correction_j = correction;
        Ok(ExactLinearSolidStep {
            step,
            time_s,
            mean_velocity_mps,
            energy,
        })
    }

    // Same cached complete-mode propagation for ordinary affine intervals and
    // read-only held-force contact responses. This is the only oscillator law.
    fn affine_mode(&self, i: usize, g: f64, dg: f64) -> [f64; 3] {
        let a = self.oscillators[i];
        let (q0, v0) = (self.modal_displacement[i], self.modal_velocity[i]);
        let q = a.cosine.mul_add(
            q0,
            a.h_sinc.mul_add(v0, a.h2_cosc.mul_add(g, a.h2_sinc3 * dg)),
        );
        let v = a
            .minus_h_lambda_sinc
            .mul_add(q0, a.cosine.mul_add(v0, a.h_sinc.mul_add(g, a.h_cosc * dg)));
        let mean = a.sinc.mul_add(
            q0,
            a.h_cosc.mul_add(v0, a.h2_sinc3.mul_add(g, a.h2_cosc4 * dg)),
        );
        [q, v, mean]
    }

    fn next_interval(&self) -> Result<(u64, f64)> {
        let step = self
            .steps
            .checked_add(1)
            .filter(|s| *s <= self.controls.maximum_steps)
            .ok_or_else(|| invalid("accepted step budget exhausted"))?;
        let time_s = step as f64 * self.dt_s;
        if !time_s.is_finite() || time_s <= self.time_s() {
            return Err(invalid("representable later clock required"));
        }
        Ok((step, time_s))
    }
}

fn check_operator(defect: f64, allowance: f64) -> Result<()> {
    if !defect.is_finite() || !allowance.is_finite() || defect > allowance {
        Err(CompleteSpectralSolidError::Operator { defect, allowance })
    } else {
        Ok(())
    }
}
fn symmetric(matrix: &Csr, n: usize, tolerance: f64, cx: &Cx<'_>) -> Result<Vec<f64>> {
    if matrix.nrows() != n || matrix.ncols() != n {
        return Err(invalid("source matrix shape differs from complete basis"));
    }
    let mut dense = zeros(n * n)?;
    for i in 0..n {
        poll(cx)?;
        let (columns, values) = matrix.row(i);
        for (&j, &value) in columns.iter().zip(values) {
            if !value.is_finite() {
                return Err(invalid("source matrix must be finite"));
            }
            dense[i * n + j] = value;
        }
    }
    let scale = dense.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
    for i in 0..n {
        poll(cx)?;
        for j in 0..i {
            let defect = if scale == 0.0 {
                0.0
            } else {
                (dense[i * n + j] / scale - dense[j * n + i] / scale).abs()
            };
            check_operator(defect, tolerance)?;
            let mean = 0.5 * dense[i * n + j] + 0.5 * dense[j * n + i];
            dense[i * n + j] = mean;
            dense[j * n + i] = mean;
        }
    }
    Ok(dense)
}
fn project(matrix: &[f64], x: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>> {
    let mut result = zeros(x.len())?;
    project_into(matrix, x, &mut result, cx)?;
    Ok(result)
}
fn project_into(matrix: &[f64], x: &[f64], result: &mut [f64], cx: &Cx<'_>) -> Result<()> {
    let n = x.len();
    for j in 0..n {
        poll(cx)?;
        result[j] = (0..n).fold(0.0, |s, i| matrix[i * n + j].mul_add(x[i], s));
        if !result[j].is_finite() {
            return Err(invalid("projection overflow"));
        }
    }
    Ok(())
}
fn reconstruct(matrix: &[f64], x: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>> {
    let mut result = zeros(x.len())?;
    reconstruct_into(matrix, x, &mut result, cx)?;
    Ok(result)
}
fn reconstruct_into(matrix: &[f64], x: &[f64], result: &mut [f64], cx: &Cx<'_>) -> Result<()> {
    let n = x.len();
    for i in 0..n {
        poll(cx)?;
        result[i] = (0..n).fold(0.0, |s, j| matrix[i * n + j].mul_add(x[j], s));
        if !result[i].is_finite() {
            return Err(invalid("reconstruction overflow"));
        }
    }
    Ok(())
}
fn quadratic(matrix: &[f64], x: &[f64], cx: &Cx<'_>) -> Result<f64> {
    let n = x.len();
    let (mut result, mut correction) = (0.0, 0.0);
    for i in 0..n {
        poll(cx)?;
        let row = (0..n).fold(0.0, |s, j| matrix[i * n + j].mul_add(x[j], s));
        add(&mut result, &mut correction, 0.5 * x[i] * row);
    }
    if !result.is_finite() {
        return Err(invalid("physical energy overflow"));
    }
    Ok(result)
}
fn quadratic_change(matrix: &[f64], before: &[f64], after: &[f64], cx: &Cx<'_>) -> Result<f64> {
    let n = before.len();
    let (mut result, mut correction) = (0.0, 0.0);
    for i in 0..n {
        poll(cx)?;
        let row = (0..n).fold(0.0, |s, j| {
            matrix[i * n + j].mul_add(after[j] - before[j], s)
        });
        add(
            &mut result,
            &mut correction,
            (0.5 * before[i] + 0.5 * after[i]) * row,
        );
    }
    if !result.is_finite() {
        return Err(invalid("physical energy change overflow"));
    }
    Ok(result)
}
