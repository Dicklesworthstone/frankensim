//! Exact stress-ramp creep of the original anisotropic Maxwell history.
//!
//! For original Gram rows B, branch coordinates q and rates R, let
//! x = sqrt(R) q and H = sqrt(R) B. Eliminating strain rate gives
//! x' = -(R - H C0^-1 H^T) x + H C0^-1 stress', with
//! C0 = C_eq + B^T B. Preparation uses a positive compliance factor and the
//! native QR/SVD owners, avoiding subtraction of nearly equal stiffnesses.
//! Each mode then uses the existing exact Prony ramp. Heat integrates |x|²
//! independently; work integrates the prescribed stress against strain rate.
//!
//! [`PreparedTensorMixedControl`] uses the same full-memory evolution with any
//! combination of prescribed stress and strain components. In particular,
//! plane-stress histories retain all original three-dimensional branch memory.
//! Coefficients are fixed. Nonzero Kelvin viscosity is a different differential
//! problem and refuses here. No thermal evolution or certified forward-error
//! bound is inferred.

use super::{
    TensorProny, TensorPronyError, TensorPronyHistory, Voigt, apply, dot, finite, invalid, mandel,
    poll, tensor, transpose_add,
};
use crate::Tangent6;
use crate::visco::{PronyRamp, stress_control::modal_integrals};
use fs_la::factor::Cholesky;
use fs_math::dd::Dd;
use std::sync::Arc;

mod mixed;
mod prepare;
pub use mixed::PreparedTensorMixedControl;

type Result<T> = core::result::Result<T, TensorPronyError>;

/// Explicit ceilings for the dense, reusable constitutive preparation.
#[derive(Clone, Copy, Debug)]
pub struct TensorStressControlBudget {
    /// Six memory coordinates per original branch, including null directions.
    pub max_memory_coordinates: usize,
    /// Conservative peak scalar storage, including native factor workspaces.
    pub max_dense_scalars: usize,
    /// Conservative scalar-work ceiling for QR/SVD preparation and its checks.
    /// Native noninterruptible factors also have a fixed bounded-work ceiling.
    pub max_preparation_work: usize,
}

impl Default for TensorStressControlBudget {
    fn default() -> Self {
        Self {
            max_memory_coordinates: 48,
            max_dense_scalars: 100_000,
            max_preparation_work: fs_la::eigen::MAX_EIGEN_UNPOLLED_SCALAR_WORK,
        }
    }
}

/// One accepted stress ramp; all energies are densities in J/m³.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorStressRampStep {
    /// Endpoint total stress in tensor-shear Voigt order [Pa]. Under mixed
    /// control, this includes the reactions in strain-controlled components.
    pub stress: Voigt,
    /// Resulting endpoint strain in tensor-shear Voigt order.
    pub strain: Voigt,
    /// Integrated strain increment before adding a possibly large old strain.
    pub strain_increment: Voigt,
    /// Original equilibrium and branch recoverable energy at the endpoint.
    pub stored_energy: f64,
    /// Independently evaluated local quadratic storage change.
    pub stored_energy_change: f64,
    /// Signed integral of stress : strain_rate over this ramp.
    pub work: f64,
    /// Nonnegative independent integral of the original dashpot heat rate.
    pub dissipated: f64,
    /// Work minus local storage change minus heat; never charged as heat.
    pub balance_residual: f64,
    /// Conservative componentwise whole-ramp strain bounds, tensor-shear order.
    /// These numerical bounds are not outward-rounded interval certificates.
    pub strain_bounds: [Voigt; 2],
}

/// Prepared full-tensor creep, bound to an immutable original material model.
///
/// Prepare from any existing history, then alternate stress and strain control
/// without translating or resetting its branch memory. Clones of that history
/// retain the binding. Independently created histories refuse, even if their
/// public material coefficients happen to compare equal.
///
/// Preparation is O((6*branches)^3); each ramp is O((6*branches)^2).
/// Original zero/rank-deficient branch factors are retained. Cancellation polls
/// surround bounded native factors and occur within assembly and every mode;
/// a refused ramp leaves every history field unchanged.
#[derive(Clone, Debug)]
pub struct PreparedTensorStressControl {
    model: Arc<TensorProny>,
    instantaneous: Option<Cholesky>,
    instantaneous_stiffness: Tangent6,
    stress_components: [bool; 6],
    stiffness_scale: f64,
    branch_rate_roots: Vec<f64>,
    rates: Vec<f64>,
    modes: Vec<f64>,
    forcing: Vec<Voigt>,
}

fn checked(value: f64, what: &'static str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(what))
    }
}

fn zeros(n: usize) -> Result<Vec<f64>> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| invalid("stress-control allocation refused"))?;
    values.resize(n, 0.0);
    Ok(values)
}

fn dd_dot(a: &Voigt, b: &Voigt) -> Dd {
    a.iter().zip(b).fold(Dd::ZERO, |s, (&a, &b)| {
        s + Dd::from_f64(a) * Dd::from_f64(b)
    })
}

impl PreparedTensorStressControl {
    fn compliance(&self, value: Voigt) -> Result<Voigt> {
        let mut compressed = [0.0; 6];
        let mut count = 0;
        for (axis, &controlled) in self.stress_components.iter().enumerate() {
            if controlled {
                compressed[count] = value[axis] / self.stiffness_scale;
                count += 1;
            }
        }
        if let Some(factor) = &self.instantaneous {
            factor.solve(&mut compressed[..count]);
        }
        let mut answer = [0.0; 6];
        count = 0;
        for (axis, &controlled) in self.stress_components.iter().enumerate() {
            if controlled {
                answer[axis] = compressed[count];
                count += 1;
            }
        }
        if !finite(&answer) {
            return Err(invalid("unrepresentable instantaneous compliance response"));
        }
        Ok(answer)
    }

    /// Number of original retained six-coordinate memory channels.
    #[must_use]
    pub fn memory_coordinates(&self) -> usize {
        self.rates.len()
    }

    /// Advance a piecewise-linear prescribed stress, including a hold, reversal
    /// or instantaneous elastic jump. Stress uses tensor-shear Voigt order.
    ///
    /// The complete continuous strain path must fit the original small-strain
    /// limit, using conservative sums of monotone modal contributions. A bound
    /// can refuse a path whose exact maximum would fit. Finite rounded source
    /// closure and a local work/storage/heat balance are required before commit.
    ///
    /// # Errors
    /// Model mismatch, invalid/unrepresentable data, an unresolved physical
    /// balance, an exceeded whole-path strain bound, or cancellation.
    #[allow(clippy::too_many_lines)] // One atomic source-bound constitutive transaction.
    pub fn advance(
        &self,
        history: &mut TensorPronyHistory,
        stress: Voigt,
        duration_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<TensorStressRampStep> {
        self.advance_controlled(history, stress, duration_s, cancelled)
    }

    #[allow(clippy::too_many_lines)] // Shared full-memory transaction for every control partition.
    fn advance_controlled(
        &self,
        history: &mut TensorPronyHistory,
        prescribed: Voigt,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TensorStressRampStep> {
        poll(&mut cancelled)?;
        if !Arc::ptr_eq(&self.model, &history.model) {
            return Err(invalid(
                "stress controller and history have different material bindings",
            ));
        }
        if !duration_s.is_finite() || duration_s < 0.0 || !finite(&prescribed) {
            return Err(invalid("finite controls and nonnegative duration required"));
        }
        let target = mandel(&prescribed);
        let delta = core::array::from_fn(|i| {
            target[i]
                - if self.stress_components[i] {
                    history.stress[i]
                } else {
                    history.strain[i]
                }
        });
        if !finite(&target) || !finite(&delta) {
            return Err(invalid("unrepresentable prescribed stress increment"));
        }
        if duration_s == 0.0 && delta.iter().all(|&v| v == 0.0) {
            return Ok(TensorStressRampStep {
                stress: history.stress(),
                strain: history.strain(),
                strain_increment: [0.0; 6],
                stored_energy: history.stored,
                stored_energy_change: 0.0,
                work: 0.0,
                dissipated: 0.0,
                balance_residual: 0.0,
                strain_bounds: [history.strain(); 2],
            });
        }
        let strain_delta = core::array::from_fn(|i| {
            if self.stress_components[i] {
                0.0
            } else {
                delta[i]
            }
        });
        let strain_loading = apply(&self.instantaneous_stiffness, &strain_delta);
        let rhs = core::array::from_fn(|i| {
            if self.stress_components[i] {
                delta[i] - strain_loading[i]
            } else {
                0.0
            }
        });
        let response = self.compliance(rhs)?;
        let elastic = core::array::from_fn(|i| response[i] + strain_delta[i]);
        let mut increment = elastic.map(Dd::from_f64);
        let mut lower = elastic.map(|v| v.min(0.0));
        let mut upper = elastic.map(|v| v.max(0.0));
        let base_stress = apply(&self.instantaneous_stiffness, &elastic);
        let mean_stress = core::array::from_fn(|i| {
            if self.stress_components[i] {
                0.5 * history.stress[i] + 0.5 * target[i]
            } else {
                history.stress[i] + 0.5 * base_stress[i]
            }
        });
        let mut work = dd_dot(&mean_stress, &elastic);
        // Stress ports integrate sigma_F * epsilon_F'; strain ports integrate
        // their evolving reaction sigma_K against the imposed epsilon_K'.
        let work_initial = core::array::from_fn(|i| {
            if self.stress_components[i] {
                history.stress[i]
            } else {
                -delta[i]
            }
        });
        let mut heat = Dd::ZERO;
        let n = self.rates.len();
        let mut responses = Vec::new();
        responses
            .try_reserve_exact(n)
            .map_err(|_| invalid("stress-control workspace allocation refused"))?;
        if duration_s > 0.0 {
            for mode in 0..n {
                poll(&mut cancelled)?;
                let mut y0 = Dd::ZERO;
                for row in 0..n {
                    let value = self.branch_rate_roots[row / 6] * history.q[row / 6][row % 6];
                    y0 = y0 + Dd::from_f64(self.modes[row * n + mode]) * Dd::from_f64(value);
                }
                let initial = checked(y0.to_f64(), "unrepresentable creep modal memory")?;
                let forcing = checked(
                    dd_dot(&self.forcing[mode], &delta).to_f64(),
                    "unrepresentable creep modal forcing",
                )?;
                let rate = self.rates[mode];
                let r = duration_s * rate;
                if (r == 0.0 && duration_s > 0.0) || (r.is_infinite() && forcing != 0.0) {
                    return Err(invalid("unrepresentable stress-ramp relaxation number"));
                }
                let ramp = PronyRamp::new(r);
                let (_, mean, _) = ramp.apply(initial, forcing);
                let (zeroth, first, parts) =
                    modal_integrals(initial, forcing, rate, duration_s, mean, &ramp);
                let mode_heat = ramp.squared_integral(initial, forcing, duration_s, rate);
                if ![zeroth, first, mode_heat].iter().all(|v| v.is_finite()) || mode_heat < 0.0 {
                    return Err(invalid("unrepresentable integrated creep response"));
                }
                heat = heat + Dd::from_f64(mode_heat);
                for axis in 0..6 {
                    if !self.stress_components[axis] {
                        continue;
                    }
                    let direction = self.forcing[mode][axis];
                    increment[axis] =
                        increment[axis] + Dd::from_f64(direction) * Dd::from_f64(zeroth);
                    for part in parts {
                        let contribution = direction * part;
                        lower[axis] += contribution.min(0.0);
                        upper[axis] += contribution.max(0.0);
                    }
                }
                work = work
                    + dd_dot(&work_initial, &self.forcing[mode]) * Dd::from_f64(zeroth)
                    + Dd::from_f64(forcing) * Dd::from_f64(first);
                responses.push((ramp, initial, forcing));
            }
        }
        let mut de = increment.map(Dd::to_f64);
        let mut dq = vec![[0.0; 6]; history.q.len()];
        for row in 0..n {
            if duration_s == 0.0 {
                let branch = &self.model.branches()[row / 6];
                let b = branch.resolved_gram_rows(&self.model.factors[row / 6]);
                dq[row / 6][row % 6] = dd_dot(&b[row % 6], &de).to_f64();
            } else {
                let mut change = Dd::ZERO;
                for (mode, &(ramp, initial, forcing)) in responses.iter().enumerate() {
                    // Scale amplitudes BEFORE the tiny decay decrement, so a
                    // weak physical memory change is not lost in heat units.
                    let coefficient = self.modes[row * n + mode] / self.branch_rate_roots[row / 6];
                    change = change
                        + Dd::from_f64(
                            ramp.spring_increment(coefficient * initial, coefficient * forcing),
                        );
                }
                dq[row / 6][row % 6] = change.to_f64();
            }
        }
        let mut strain = core::array::from_fn(|i| {
            if self.stress_components[i] {
                history.strain[i] + de[i]
            } else {
                target[i]
            }
        });
        let mut q = history.q.clone();
        for (next, change) in q.iter_mut().zip(&dq) {
            for axis in 0..6 {
                next[axis] += change[axis];
            }
        }
        let (calculated, scale) = self.source_stress(&strain, &q)?;
        let residual = core::array::from_fn(|i| {
            if self.stress_components[i] {
                target[i] - calculated[i]
            } else {
                0.0
            }
        });
        let tolerance = 8192.0 * f64::EPSILON * (n + 6) as f64;
        for axis in 0..6 {
            if residual[axis].abs() > tolerance * (scale[axis] + target[axis].abs()) {
                return Err(invalid(
                    "stress-control source reconstruction is unresolved",
                ));
            }
        }
        // Close only admitted arithmetic roundoff through C0. This is not a
        // physical extra stress jump and its defect is visible in the balance.
        let correction = self.compliance(residual)?;
        for axis in 0..6 {
            de[axis] += correction[axis];
            if self.stress_components[axis] {
                strain[axis] = history.strain[axis] + de[axis];
            }
        }
        for (j, next) in q.iter_mut().enumerate() {
            let b = self.model.branches()[j].resolved_gram_rows(&self.model.factors[j]);
            let dc = apply(&b, &correction);
            for axis in 0..6 {
                dq[j][axis] += dc[axis];
                next[axis] = history.q[j][axis] + dq[j][axis];
            }
        }
        let mut bounds = [[0.0; 6]; 2];
        let mut maximum = [0.0; 6];
        for axis in 0..6 {
            bounds[0][axis] = (history.strain[axis] + lower[axis]).min(strain[axis]);
            bounds[1][axis] = (history.strain[axis] + upper[axis]).max(strain[axis]);
            maximum[axis] = bounds[0][axis].abs().max(bounds[1][axis].abs());
        }
        self.model
            .validate_mandel_strain(&maximum)
            .map_err(|_| invalid("whole stress-ramp path exceeds the material strain bound"))?;
        let eq0 = apply(&self.model.equilibrium_gram, &history.strain);
        let eq_delta = apply(&self.model.equilibrium_gram, &de);
        let eq1 = apply(&self.model.equilibrium_gram, &strain);
        let mut storage_change = Dd::ZERO;
        let mut exchange_scale = 0.0;
        let mut stored = 0.5 * dot(&eq1, &eq1);
        for (old, change) in core::iter::once((&eq0, &eq_delta)).chain(history.q.iter().zip(&dq)) {
            for axis in 0..6 {
                let term = (Dd::from_f64(old[axis]) + Dd::from_f64(0.5 * change[axis]))
                    * Dd::from_f64(change[axis]);
                storage_change = storage_change + term;
                exchange_scale += term.abs().to_f64();
            }
        }
        for memory in &q {
            stored += 0.5 * dot(memory, memory);
        }
        let balance = (work - storage_change - heat).to_f64();
        let work = work.to_f64();
        let heat = heat.to_f64();
        let scale = work.abs() + heat + exchange_scale;
        let cumulative_work = history.work + work;
        let cumulative_heat = history.heat + heat;
        if !finite(&strain)
            || q.iter().any(|q| !finite(q))
            || ![
                stored,
                work,
                heat,
                balance,
                scale,
                cumulative_work,
                cumulative_heat,
            ]
            .iter()
            .all(|v| v.is_finite())
            || heat < 0.0
            || balance.abs() > tolerance * scale
        {
            return Err(invalid(
                "stress-control local work/storage/heat balance is unresolved",
            ));
        }
        let (closed, scale) = self.source_stress(&strain, &q)?;
        for axis in 0..6 {
            if self.stress_components[axis]
                && (closed[axis] - target[axis]).abs()
                    > tolerance * (scale[axis] + target[axis].abs())
            {
                return Err(invalid(
                    "stress-control rounded source closure is unresolved",
                ));
            }
        }
        poll(&mut cancelled)?;
        let next_stress = core::array::from_fn(|i| {
            if self.stress_components[i] {
                target[i]
            } else {
                closed[i]
            }
        });
        let physical_stress = tensor(&next_stress);
        let report = TensorStressRampStep {
            stress: core::array::from_fn(|i| {
                if self.stress_components[i] {
                    prescribed[i]
                } else {
                    physical_stress[i]
                }
            }),
            strain: tensor(&strain),
            strain_increment: tensor(&de),
            stored_energy: stored,
            stored_energy_change: storage_change.to_f64(),
            work,
            dissipated: heat,
            balance_residual: balance,
            strain_bounds: [tensor(&bounds[0]), tensor(&bounds[1])],
        };
        history.strain = strain;
        history.stress = next_stress;
        history.q = q;
        history.stored = stored;
        history.work = cumulative_work;
        history.heat = cumulative_heat;
        Ok(report)
    }

    fn source_stress(&self, strain: &Voigt, q: &[Voigt]) -> Result<(Voigt, Voigt)> {
        let eq = &self.model.equilibrium_gram;
        let mut stress = [0.0; 6];
        let mut magnitude = [0.0; 6];
        let eq_strain = apply(eq, strain);
        transpose_add(eq, &eq_strain, &mut stress);
        for axis in 0..6 {
            for row in 0..6 {
                magnitude[axis] += (eq[row][axis] * eq_strain[row]).abs();
            }
        }
        for (j, memory) in q.iter().enumerate() {
            let b = self.model.branches()[j].resolved_gram_rows(&self.model.factors[j]);
            transpose_add(&b, memory, &mut stress);
            for axis in 0..6 {
                for row in 0..6 {
                    magnitude[axis] += (b[row][axis] * memory[row]).abs();
                }
            }
        }
        if !finite(&stress) || !finite(&magnitude) {
            return Err(invalid("unrepresentable original constitutive stress"));
        }
        Ok((stress, magnitude))
    }
}
