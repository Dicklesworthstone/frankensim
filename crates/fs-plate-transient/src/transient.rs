//! Spatial DKT plate dynamics through the shared generalized-alpha integrator.
//!
//! For one fixed assembled mesh/support state, `M = s_m M0`, `K = s_k K0`
//! and `C = alpha M + beta K`. All actions use the existing sparse matrices;
//! no dense operator, modal truncation or second time integrator is introduced.
//! The four model parameters and any applied-load parameters share one discrete
//! adjoint. Scaling the whole K also scales prestress and stiffeners, if present;
//! it is a Young-modulus derivative only when those contributions scale with it.
//!
//! Geometry, support elimination and reference matrices remain fixed. This
//! inherits the supplied spatial DKT/inertia approximation and small-deflection plate
//! assumptions. Nonnegative Rayleigh coefficients dissipate energy when the
//! supplied stiffness is positive semidefinite; this adapter does not establish
//! that spectral property or identify physical damping from material data.

use crate::PlateModel;
use fs_solver::FlexiblePreconditioner;
use fs_sparse::{Csr, DirectOrdering, LdltError, LdltOptions, SymbolicLdlt};
use fs_time::galpha::second_order_adjoint::SecondOrderVjp;
use fs_time::galpha::second_order_adjoint::trajectory::StructuralTrajectoryModel;
use fs_time::galpha::{SecondOrderOperatorWeights, SecondOrderProblem};

/// Applied generalized forces in the reduced `(w, wx, wy)` coordinates.
/// Callbacks are pure at a fixed parameter point, overwrite all outputs, and
/// bound their own work/storage. Time may change the load location and value.
/// Model-dependent locations/loads must include their complete derivative here.
pub trait PlateLoad {
    /// Parameters appended after the four dynamics parameters.
    fn parameter_count(&self) -> usize;
    /// Applied load at the requested physical time, in reduced DOF order.
    fn forcing(&self, time: f64, output: &mut [f64]) -> Result<(), String>;
    /// Load-parameter VJP at fixed time, with the supplied reduced force seed.
    /// Outputs cover only this load's parameters, not the four dynamics entries.
    fn forcing_vjp(&self, time: f64, seed: &[f64], parameter_bar: &mut [f64])
    -> Result<(), String>;
}

/// Parameter-independent zero load for free vibration.
#[derive(Debug, Clone, Copy)]
pub struct NoPlateLoad;
impl PlateLoad for NoPlateLoad {
    fn parameter_count(&self) -> usize {
        0
    }
    fn forcing(&self, _: f64, output: &mut [f64]) -> Result<(), String> {
        output.fill(0.0);
        Ok(())
    }
    fn forcing_vjp(&self, _: f64, _: &[f64], output: &mut [f64]) -> Result<(), String> {
        output.fill(0.0);
        Ok(())
    }
}

/// Fixed model point; adjoint entries follow this field order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateDynamicsParameters {
    /// Positive multiplier of the complete reference stiffness.
    pub stiffness_scale: f64,
    /// Positive multiplier of the complete reference mass, including cross terms.
    pub mass_scale: f64,
    /// Nonnegative mass-proportional damping coefficient, 1/s.
    pub mass_damping_per_s: f64,
    /// Nonnegative stiffness-proportional damping coefficient, s.
    pub stiffness_damping_s: f64,
}

/// Admission ceilings checked before scanning operators or allocating scratch.
/// The model borrows existing matrices. During each action it needs no scratch;
/// initialization and integration retain their own explicit workspace limits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PlateDynamicsBudget {
    /// Maximum number of reduced displacement/slope degrees of freedom.
    pub max_dofs: usize,
    /// Sum of stored entries in the reference K and M matrices.
    pub max_nonzeros: usize,
    /// Maximum number of entries in the full-to-reduced support map.
    pub max_full_dofs: usize,
    /// Maximum number of parameters owned by the applied load.
    pub max_load_parameters: usize,
}

/// Admission, preconditioning or cancellation refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum PlateDynamicsError {
    /// Incompatible model, nonfinite inputs or inadmissible physical controls.
    InvalidInput(&'static str),
    /// The model exceeded an explicit pre-scan admission ceiling.
    Budget(&'static str),
    /// The original sparse mass factorization refused, including its byte cap.
    MassFactor(LdltError),
    /// Cancellation observed before publishing a model or preconditioner.
    Cancelled,
}
impl std::fmt::Display for PlateDynamicsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "plate dynamics failed: {self:?}")
    }
}
impl std::error::Error for PlateDynamicsError {}

fn poll<C: FnMut() -> bool>(cancelled: &mut C) -> Result<(), PlateDynamicsError> {
    if cancelled() {
        Err(PlateDynamicsError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(values: &[f64]) -> bool {
    values.iter().all(|value| value.is_finite())
}
fn dot_row(matrix: &Csr, row: usize, input: &[f64]) -> f64 {
    let (columns, values) = matrix.row(row);
    columns
        .iter()
        .zip(values)
        .fold(0.0, |sum, (&column, &value)| {
            value.mul_add(input[column], sum)
        })
}

/// One immutable spatial plate and load parameter point.
///
/// Implements the existing structural residual, explicit transposed actions,
/// parameter VJP and time-dependent trajectory loading. The supplied mass may
/// be lumped or consistent; positive definiteness is checked once at admission.
/// K may include prestress and stiffeners. Transposed actions use the original
/// stored entries, including all mass cross terms, without assuming bitwise
/// symmetry of floating-point assembly.
pub struct PlateDynamics<'a, L: ?Sized> {
    model: &'a PlateModel,
    parameters: PlateDynamicsParameters,
    load: &'a L,
    load_parameters: usize,
}

impl<L: ?Sized> Copy for PlateDynamics<'_, L> {}
impl<L: ?Sized> Clone for PlateDynamics<'_, L> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, L: PlateLoad + ?Sized> PlateDynamics<'a, L> {
    pub(crate) fn same_model_and_load(&self, other: &Self) -> bool {
        core::ptr::eq(self.model, other.model)
            && core::ptr::eq(self.load, other.load)
            && self.load_parameters == other.load_parameters
            && self.parameters == other.parameters
    }

    /// Bind admitted sparse operators and pure load callbacks.
    pub fn new<C: FnMut() -> bool>(
        model: &'a PlateModel,
        parameters: PlateDynamicsParameters,
        load: &'a L,
        budget: PlateDynamicsBudget,
        cancelled: &mut C,
    ) -> Result<Self, PlateDynamicsError> {
        Self::new_with_mass_options(
            model,
            parameters,
            load,
            budget,
            &LdltOptions::default(),
            cancelled,
        )
    }

    /// Bind the same physical model with explicit consistent-mass factor limits.
    /// Diagonal mass uses the original linear-time admission and needs no factor.
    /// For consistent mass, the existing sparse LDLT service checks positive
    /// definiteness once; its factors are then discarded. No mass lumping or
    /// modified physical operator is introduced. Set `max_factor_bytes` to bound
    /// the numeric factor and scratch separately from the source operator budget.
    /// Cancellation is checked around the bounded synchronous factorization.
    #[allow(clippy::too_many_lines)]
    pub fn new_with_mass_options<C: FnMut() -> bool>(
        model: &'a PlateModel,
        parameters: PlateDynamicsParameters,
        load: &'a L,
        budget: PlateDynamicsBudget,
        mass_options: &LdltOptions,
        cancelled: &mut C,
    ) -> Result<Self, PlateDynamicsError> {
        poll(cancelled)?;
        let n = model.free;
        let load_parameters = load.parameter_count();
        if n > budget.max_dofs
            || model.dof_map.len() > budget.max_full_dofs
            || load_parameters > budget.max_load_parameters
            || load_parameters.checked_add(4).is_none()
            || model
                .k
                .nnz()
                .checked_add(model.m.nnz())
                .is_none_or(|nnz| nnz > budget.max_nonzeros)
        {
            return Err(PlateDynamicsError::Budget(
                "model exceeds operator, DOF or parameter cap",
            ));
        }
        let p = parameters;
        if !finite(&[
            p.stiffness_scale,
            p.mass_scale,
            p.mass_damping_per_s,
            p.stiffness_damping_s,
            p.mass_scale * p.mass_damping_per_s,
            p.stiffness_scale * p.stiffness_damping_s,
        ]) || p.stiffness_scale <= 0.0
            || p.mass_scale <= 0.0
            || p.mass_damping_per_s < 0.0
            || p.stiffness_damping_s < 0.0
        {
            return Err(PlateDynamicsError::InvalidInput(
                "positive finite scales and nonnegative finite Rayleigh coefficients required",
            ));
        }
        if n == 0
            || !model.dof_map.len().is_multiple_of(3)
            || model.k.nrows() != n
            || model.k.ncols() != n
            || model.m.nrows() != n
            || model.m.ncols() != n
        {
            return Err(PlateDynamicsError::InvalidInput(
                "compatible nonempty reduced plate matrices required",
            ));
        }
        let mut seen = vec![false; n];
        for chunk in model.dof_map.chunks(256) {
            poll(cancelled)?;
            for &entry in chunk {
                if let Some(index) = entry {
                    if index >= n || seen[index] {
                        return Err(PlateDynamicsError::InvalidInput(
                            "invalid or duplicate reduced DOF mapping",
                        ));
                    }
                    seen[index] = true;
                }
            }
        }
        if seen.iter().any(|present| !present) {
            return Err(PlateDynamicsError::InvalidInput(
                "incomplete reduced DOF mapping",
            ));
        }
        let mut diagonal_mass = true;
        for row in 0..n {
            poll(cancelled)?;
            let (columns, values) = model.m.row(row);
            let diagonal = model.m.get(row, row);
            if !diagonal.is_finite()
                || diagonal <= 0.0
                || p.mass_scale * diagonal <= 0.0
                || !(p.mass_scale * diagonal).is_finite()
            {
                return Err(PlateDynamicsError::InvalidInput(
                    "finite positive scaled plate mass diagonal required",
                ));
            }
            for (column_chunk, value_chunk) in columns.chunks(256).zip(values.chunks(256)) {
                poll(cancelled)?;
                for (&column, &value) in column_chunk.iter().zip(value_chunk) {
                    if !value.is_finite()
                        || !(p.mass_scale * value).is_finite()
                        || !(p.mass_scale * p.mass_damping_per_s * value).is_finite()
                    {
                        return Err(PlateDynamicsError::InvalidInput(
                            "finite representable mass and damping entries required",
                        ));
                    }
                    if column != row && value != 0.0 {
                        diagonal_mass = false;
                        let mirror = model.m.get(column, row);
                        // Admit assembly roundoff only; do not silently use one
                        // triangle of a genuinely nonsymmetric inertia operator.
                        let scale = diagonal.sqrt() * model.m.get(column, column).sqrt();
                        if !scale.is_finite()
                            || !mirror.is_finite()
                            || (value - mirror).abs() > 64.0 * f64::EPSILON * scale
                        {
                            return Err(PlateDynamicsError::InvalidInput(
                                "consistent plate mass must be symmetric to assembly roundoff",
                            ));
                        }
                    }
                }
            }
            for chunk in model.k.row(row).1.chunks(256) {
                poll(cancelled)?;
                if chunk.iter().any(|&value| {
                    !value.is_finite()
                        || !(p.stiffness_scale * value).is_finite()
                        || !(p.stiffness_scale * p.stiffness_damping_s * value).is_finite()
                }) {
                    return Err(PlateDynamicsError::InvalidInput(
                        "finite representable stiffness and damping entries required",
                    ));
                }
            }
        }
        if !diagonal_mass {
            let symbolic = SymbolicLdlt::analyze(&model.m, DirectOrdering::Amd)
                .map_err(PlateDynamicsError::MassFactor)?;
            poll(cancelled)?;
            let factor = symbolic
                .factor(&model.m, mass_options)
                .map_err(PlateDynamicsError::MassFactor)?;
            if factor.inertia().negative != 0 {
                return Err(PlateDynamicsError::InvalidInput(
                    "consistent plate mass must be positive definite",
                ));
            }
        }
        poll(cancelled)?;
        Ok(Self {
            model,
            parameters,
            load,
            load_parameters,
        })
    }

    /// The immutable reference spatial pencil and support map.
    #[must_use]
    pub fn model(&self) -> &PlateModel {
        self.model
    }

    /// Jacobi inverse of the scaled mass diagonal, for both initialization
    /// primal and transposed mass solves. It is exact for lumped mass only;
    /// consistent mass retains its full sparse action in the iterative solve.
    /// Supply sufficient Krylov restart/cycle limits for that physical operator.
    #[must_use]
    pub fn mass_preconditioner(&self) -> PlateMassPreconditioner<'_> {
        PlateMassPreconditioner {
            mass: &self.model.m,
            scale: self.parameters.mass_scale,
        }
    }

    /// Jacobi preconditioner for the displacement-based generalized-alpha
    /// effective operator (and its transpose) at fixed step and spectral radius.
    /// Pass this to the trajectory's explicit adjoint argument. Forward Newton
    /// solves automatically use the same diagonal through the model hook, with
    /// the weights supplied by the integrator at its current step size.
    pub fn effective_preconditioner<C: FnMut() -> bool>(
        &self,
        step: f64,
        rho_inf: f64,
        cancelled: &mut C,
    ) -> Result<PlateJacobi, PlateDynamicsError> {
        poll(cancelled)?;
        if !step.is_finite() || step <= 0.0 || !(0.0..=1.0).contains(&rho_inf) {
            return Err(PlateDynamicsError::InvalidInput(
                "positive finite step and spectral radius in [0,1] required",
            ));
        }
        let alpha_m = (2.0 * rho_inf - 1.0) / (rho_inf + 1.0);
        let alpha_f = rho_inf / (rho_inf + 1.0);
        let gamma = 0.5 - alpha_m + alpha_f;
        let beta = 0.25 * (1.0 - alpha_m + alpha_f) * (1.0 - alpha_m + alpha_f);
        let (m_coefficient, k_coefficient) =
            self.effective_coefficients(SecondOrderOperatorWeights {
                mass: (1.0 - alpha_m) / (beta * step * step),
                damping: (1.0 - alpha_f) * gamma / (beta * step),
                tangent: 1.0 - alpha_f,
            });
        if !finite(&[m_coefficient, k_coefficient]) {
            return Err(PlateDynamicsError::InvalidInput(
                "unrepresentable effective operator coefficients",
            ));
        }
        let mut inverse = Vec::with_capacity(self.model.free);
        for row in 0..self.model.free {
            poll(cancelled)?;
            let diagonal = m_coefficient * self.model.m.get(row, row)
                + k_coefficient * self.model.k.get(row, row);
            if !diagonal.is_finite() || diagonal <= 0.0 || !(1.0 / diagonal).is_finite() {
                return Err(PlateDynamicsError::InvalidInput(
                    "positive finite invertible effective diagonal required",
                ));
            }
            inverse.push(1.0 / diagonal);
        }
        poll(cancelled)?;
        Ok(PlateJacobi { inverse })
    }

    fn effective_coefficients(&self, weights: SecondOrderOperatorWeights) -> (f64, f64) {
        let p = self.parameters;
        (
            p.mass_scale * (weights.mass + weights.damping * p.mass_damping_per_s),
            p.stiffness_scale * (weights.tangent + weights.damping * p.stiffness_damping_s),
        )
    }

    fn action(&self, input: &[f64], output: &mut [f64], mass: f64, stiffness: f64) {
        if input.len() != self.model.free || output.len() != self.model.free {
            output.fill(f64::NAN);
            return;
        }
        for (row, out) in output.iter_mut().enumerate() {
            let inertia = if mass == 0.0 {
                0.0
            } else {
                mass * dot_row(&self.model.m, row, input)
            };
            let resistance = if stiffness == 0.0 {
                0.0
            } else {
                stiffness * dot_row(&self.model.k, row, input)
            };
            *out = inertia + resistance;
        }
    }

    fn transpose(
        &self,
        input: &[f64],
        output: &mut [f64],
        mass: f64,
        stiffness: f64,
    ) -> Result<(), String> {
        if input.len() != self.model.free || output.len() != self.model.free || !finite(input) {
            return Err("incompatible plate transpose vectors".into());
        }
        output.fill(0.0);
        for (matrix, scale) in [(&self.model.m, mass), (&self.model.k, stiffness)] {
            if scale == 0.0 {
                continue;
            }
            for (row, &seed) in input.iter().enumerate() {
                let (columns, values) = matrix.row(row);
                for (&column, &value) in columns.iter().zip(values) {
                    output[column] += scale * value * seed;
                }
            }
        }
        if finite(output) {
            Ok(())
        } else {
            Err("nonfinite plate transposed action".into())
        }
    }
}

impl<L: PlateLoad + ?Sized> SecondOrderProblem for PlateDynamics<'_, L> {
    fn dimension(&self) -> usize {
        self.model.free
    }
    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        self.action(input, output, self.parameters.mass_scale, 0.0);
    }
    fn damping_apply(&self, input: &[f64], output: &mut [f64]) {
        let p = self.parameters;
        self.action(
            input,
            output,
            p.mass_scale * p.mass_damping_per_s,
            p.stiffness_scale * p.stiffness_damping_s,
        );
    }
    fn internal_force(&self, q: &[f64], output: &mut [f64]) {
        self.action(q, output, 0.0, self.parameters.stiffness_scale);
    }
    fn tangent_apply(&self, _: &[f64], input: &[f64], output: &mut [f64]) {
        self.internal_force(input, output);
    }
    fn preconditioner_apply(
        &self,
        _: &[f64],
        weights: SecondOrderOperatorWeights,
        _: usize,
        _: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        if residual.len() != self.model.free || output.len() != self.model.free || !finite(residual)
        {
            output.fill(f64::NAN);
            return;
        }
        let (mass, stiffness) = self.effective_coefficients(weights);
        if !finite(&[mass, stiffness]) {
            output.fill(f64::NAN);
            return;
        }
        for (row, (out, &value)) in output.iter_mut().zip(residual).enumerate() {
            let diagonal =
                mass * self.model.m.get(row, row) + stiffness * self.model.k.get(row, row);
            let inverse = 1.0 / diagonal;
            if !diagonal.is_finite() || diagonal <= 0.0 || !inverse.is_finite() {
                output.fill(f64::NAN);
                return;
            }
            *out = value * inverse;
        }
        if !finite(output) {
            output.fill(f64::NAN);
        }
    }
}

impl<L: PlateLoad + ?Sized> SecondOrderVjp for PlateDynamics<'_, L> {
    fn parameter_count(&self) -> usize {
        4 + self.load_parameters
    }
    fn mass_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        self.transpose(seed, output, self.parameters.mass_scale, 0.0)
    }
    fn damping_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        let p = self.parameters;
        self.transpose(
            seed,
            output,
            p.mass_scale * p.mass_damping_per_s,
            p.stiffness_scale * p.stiffness_damping_s,
        )
    }
    fn tangent_transpose_apply(
        &self,
        _: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        self.transpose(seed, output, 0.0, self.parameters.stiffness_scale)
    }
    fn residual_parameter_vjp(
        &self,
        q: &[f64],
        v: &[f64],
        a: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        if [q, v, a, seed]
            .iter()
            .any(|x| x.len() != self.model.free || !finite(x))
            || output.len() != self.parameter_count()
        {
            return Err("incompatible plate parameter derivative vectors".into());
        }
        output.fill(0.0);
        let p = self.parameters;
        for (row, &seed) in seed.iter().enumerate() {
            let kq = dot_row(&self.model.k, row, q);
            let kv = dot_row(&self.model.k, row, v);
            let ma = dot_row(&self.model.m, row, a);
            let mv = dot_row(&self.model.m, row, v);
            output[0] += seed * (kq + p.stiffness_damping_s * kv);
            output[1] += seed * (ma + p.mass_damping_per_s * mv);
            output[2] += seed * p.mass_scale * mv;
            output[3] += seed * p.stiffness_scale * kv;
        }
        if finite(output) {
            Ok(())
        } else {
            Err("nonfinite plate parameter derivative".into())
        }
    }
}

impl<L: PlateLoad + ?Sized> StructuralTrajectoryModel for PlateDynamics<'_, L> {
    fn forcing(&self, time: f64, output: &mut [f64]) -> Result<(), String> {
        if !time.is_finite()
            || output.len() != self.model.free
            || self.load.parameter_count() != self.load_parameters
        {
            return Err("incompatible plate load time, dimension or parameter count".into());
        }
        output.fill(f64::NAN);
        self.load.forcing(time, output)?;
        if finite(output) {
            Ok(())
        } else {
            Err("nonfinite or unwritten plate load".into())
        }
    }
    fn forcing_vjp(&self, time: f64, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        if !time.is_finite()
            || seed.len() != self.model.free
            || !finite(seed)
            || output.len() != self.parameter_count()
            || self.load.parameter_count() != self.load_parameters
        {
            return Err("incompatible plate load derivative vectors".into());
        }
        output[..4].fill(0.0);
        output[4..].fill(f64::NAN);
        self.load.forcing_vjp(time, seed, &mut output[4..])?;
        if finite(output) {
            Ok(())
        } else {
            Err("nonfinite or unwritten plate load derivative".into())
        }
    }
}

/// Allocation-free Jacobi mass inverse at one fixed parameter point.
pub struct PlateMassPreconditioner<'a> {
    mass: &'a Csr,
    scale: f64,
}
impl FlexiblePreconditioner for PlateMassPreconditioner<'_> {
    fn apply(&self, _: usize, residual: &[f64], output: &mut [f64]) {
        if residual.len() != self.mass.nrows() || output.len() != self.mass.nrows() {
            output.fill(f64::NAN);
            return;
        }
        for (row, (out, &value)) in output.iter_mut().zip(residual).enumerate() {
            *out = value / (self.scale * self.mass.get(row, row));
        }
    }
}

/// Fixed Jacobi inverse for the admitted effective plate operator.
#[derive(Debug, Clone)]
pub struct PlateJacobi {
    inverse: Vec<f64>,
}
impl FlexiblePreconditioner for PlateJacobi {
    fn apply(&self, _: usize, residual: &[f64], output: &mut [f64]) {
        if residual.len() != self.inverse.len() || output.len() != self.inverse.len() {
            output.fill(f64::NAN);
            return;
        }
        for ((out, &value), &inverse) in output.iter_mut().zip(residual).zip(&self.inverse) {
            *out = value * inverse;
        }
    }
}
