//! Stress ramps with the original full Maxwell memory and SPD Kelvin viscosity.
//!
//! In energy coordinates z = [G_eq strain, q], let C = [G_eq; B] and
//! D = R_D^T R_D be the Kelvin tensor. The positive factor
//! F = [R_D^-T C^T; 0 sqrt(R)] satisfies z' = -F^T F z + C D^-1 stress.
//! Native QR/SVD prepares F = U sqrt(Lambda) V^T. Its loss coordinates
//! g = U^T [R_D^-T(C^T z - stress); sqrt(R) q] obey independent Prony
//! ramps, with strain' = -sum k_i g_i and k_i = R_D^-1 U_i[0..6].
//! Heat is the independent positive integral sum integral(g_i^2).
//!
//! This retains the original history and full noncommuting tensors. Singular
//! Kelvin tensors require differential/algebraic constraints and refuse here.
//! Fixed-temperature ramps are exact. Prescribed temperature ramps reuse this
//! controller in reduced time with bounded midpoint kinetics; elastic moduli
//! remain fixed. No thermal-state evolution, mixed loading, derivative or
//! certified forward-error bound is inferred.

use super::{
    TensorProny, TensorPronyError, TensorPronyHistory, Voigt, apply, dot, finite, invalid, mandel,
    poll, tensor,
};
use crate::visco::{PronyRamp, stress_control::modal_integrals};
use fs_la::factor::Qr;
use fs_math::dd::Dd;
use std::sync::Arc;

mod prepare;
mod temperature;
pub use temperature::TensorKelvinStressTemperatureStep;

type Result<T> = core::result::Result<T, TensorPronyError>;

/// Explicit ceilings on reusable dense constitutive preparation.
#[derive(Clone, Copy, Debug)]
pub struct TensorKelvinStressControlBudget {
    /// Six equilibrium coordinates plus six coordinates per original branch.
    pub max_state_coordinates: usize,
    /// Conservative peak scalar storage, including native factor workspaces.
    pub max_dense_scalars: usize,
    /// Scalar work ceiling, also bounded by the native unpolled-factor limit.
    pub max_preparation_work: usize,
}

impl Default for TensorKelvinStressControlBudget {
    fn default() -> Self {
        Self {
            max_state_coordinates: 48,
            max_dense_scalars: 100_000,
            max_preparation_work: fs_la::eigen::MAX_EIGEN_UNPOLLED_SCALAR_WORK,
        }
    }
}

/// Accepted ramp response; energies are densities in J/m³.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorKelvinStressRampStep {
    /// Prescribed endpoint stress, tensor-shear Voigt order [Pa].
    pub stress: Voigt,
    /// Resulting endpoint strain, tensor-shear Voigt order.
    pub strain: Voigt,
    /// Strain increment before addition to a possibly large previous strain.
    pub strain_increment: Voigt,
    /// Endpoint strain rate from the original constitutive stress equation [1/s].
    pub strain_rate: Voigt,
    /// Original equilibrium and Maxwell spring storage at the endpoint.
    pub stored_energy: f64,
    /// Local quadratic storage change, using unrounded state increments.
    pub stored_energy_change: f64,
    /// Signed integral of prescribed stress against the actual strain rate.
    pub work: f64,
    /// Independent nonnegative total Kelvin and Maxwell heat integral.
    pub dissipated: f64,
    /// Work minus local storage change minus physical heat.
    pub balance_residual: f64,
    /// Conservative componentwise whole-ramp strain bounds, tensor-shear order.
    /// These numerical bounds are not outward-rounded interval certificates.
    pub strain_bounds: [Voigt; 2],
}

/// Exact piecewise-linear stress loading of a history with SPD Kelvin viscosity.
///
/// Preparation retains all original equilibrium and branch coordinates, including
/// null branch directions, and uses O(n³) work and O(n²) storage, n=6+6*branches.
/// Each ramp uses O(n²) arithmetic and n scalar Prony preparations; it performs
/// no dense factorization. Histories can alternate between this controller and
/// original strain ramps without resetting their memory. Independently created
/// histories refuse even if their public coefficients compare equal.
#[derive(Clone, Debug)]
pub struct PreparedTensorKelvinStressControl {
    model: Arc<TensorProny>,
    viscosity_qr: Qr,
    viscosity_scale: f64,
    branch_rate_roots: Vec<f64>,
    rates: Vec<f64>,
    left_modes: Vec<f64>,
    strain_modes: Vec<Voigt>,
    memory_modes: Vec<f64>,
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
        .map_err(|_| invalid("Kelvin controller allocation refused"))?;
    values.resize(n, 0.0);
    Ok(values)
}

fn dd_dot(a: &Voigt, b: &Voigt) -> Dd {
    a.iter().zip(b).fold(Dd::ZERO, |sum, (&a, &b)| {
        sum + Dd::from_f64(a) * Dd::from_f64(b)
    })
}

struct Response {
    ramp: PronyRamp,
    initial: f64,
    forcing: f64,
    rate: f64,
    duration: f64,
    end: f64,
    integrals: (f64, f64, [f64; 2]),
}

impl Response {
    // Weight amplitudes before integrating, preserving physical increments when
    // the unweighted loss-coordinate integral is too small to represent.
    fn scaled_integrals(&self, coefficient: f64) -> Result<(f64, f64, [f64; 2])> {
        if coefficient == 0.0 {
            return Ok((0.0, 0.0, [0.0; 2]));
        }
        let (zeroth, first, parts) = self.integrals;
        let cached = (
            coefficient * zeroth,
            coefficient * first,
            parts.map(|v| coefficient * v),
        );
        // Ordinary coordinates reuse the two already prepared moments. Only
        // an underflowed/subnormal primitive needs the amplitude-first path.
        let resolved = (self.initial == 0.0 || parts[0].is_normal())
            && (self.forcing == 0.0 || parts[1].is_normal())
            && [zeroth, first].iter().all(|v| *v == 0.0 || v.is_normal());
        if resolved
            && [cached.0, cached.1, cached.2[0], cached.2[1]]
                .iter()
                .all(|v| v.is_finite())
        {
            return Ok(cached);
        }
        let initial = coefficient * self.initial;
        let forcing = coefficient * self.forcing;
        let result = if initial.is_finite() && forcing.is_finite() {
            let (_, mean, _) = self.ramp.apply(initial, forcing);
            modal_integrals(initial, forcing, self.rate, self.duration, mean, &self.ramp)
        } else {
            cached
        };
        if ![result.0, result.1, result.2[0], result.2[1]]
            .iter()
            .all(|v| v.is_finite())
        {
            return Err(invalid("unrepresentable Kelvin stress-ramp integral"));
        }
        Ok(result)
    }
}

impl PreparedTensorKelvinStressControl {
    /// Number of retained equilibrium and original branch energy coordinates.
    #[must_use]
    pub fn state_coordinates(&self) -> usize {
        self.rates.len()
    }

    fn solve_viscosity_root(&self, rhs: Voigt, transpose: bool) -> Result<Voigt> {
        let mut answer = [0.0; 6];
        for index in 0..6 {
            let i = if transpose { index } else { 5 - index };
            let mut value = Dd::from_f64(rhs[i] / self.viscosity_scale);
            if transpose {
                for (j, &v) in answer.iter().take(i).enumerate() {
                    value = value - Dd::from_f64(self.viscosity_qr.r(j, i)) * Dd::from_f64(v);
                }
            } else {
                for (j, &v) in answer.iter().enumerate().skip(i + 1) {
                    value = value - Dd::from_f64(self.viscosity_qr.r(i, j)) * Dd::from_f64(v);
                }
            }
            answer[i] = checked(
                value.to_f64() / self.viscosity_qr.r(i, i),
                "unrepresentable Kelvin compliance response",
            )?;
        }
        Ok(answer)
    }

    fn original_loss(
        &self,
        strain: &Voigt,
        q: &[Voigt],
        stress: &Voigt,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        let eq = apply(&self.model.equilibrium_gram, strain);
        let mut residual = stress.map(|v| -Dd::from_f64(v));
        let mut magnitude = stress.map(f64::abs);
        for axis in 0..6 {
            for (row, &memory) in eq.iter().enumerate() {
                let value =
                    Dd::from_f64(self.model.equilibrium_gram[row][axis]) * Dd::from_f64(memory);
                residual[axis] = residual[axis] + value;
                magnitude[axis] += value.abs().to_f64();
            }
        }
        for (j, branch) in self.model.branches().iter().enumerate() {
            let factor = branch.resolved_gram_rows(&self.model.factors[j]);
            for axis in 0..6 {
                for row in 0..6 {
                    let value = Dd::from_f64(factor[row][axis]) * Dd::from_f64(q[j][row]);
                    residual[axis] = residual[axis] + value;
                    magnitude[axis] += value.abs().to_f64();
                }
            }
        }
        let physical = self.solve_viscosity_root(residual.map(Dd::to_f64), true)?;
        let mut values = zeros(self.rates.len())?;
        let mut scale = zeros(self.rates.len())?;
        values[..6].copy_from_slice(&physical);
        for i in 0..6 {
            let mut bound = magnitude[i] / self.viscosity_scale;
            for (j, &previous) in scale.iter().take(i).enumerate() {
                bound += self.viscosity_qr.r(j, i).abs() * previous;
            }
            scale[i] = checked(
                bound / self.viscosity_qr.r(i, i).abs(),
                "unrepresentable Kelvin source scale",
            )?;
        }
        for (j, memory) in q.iter().enumerate() {
            for axis in 0..6 {
                let value = checked(
                    self.branch_rate_roots[j] * memory[axis],
                    "unrepresentable Maxwell heat coordinate",
                )?;
                values[6 + 6 * j + axis] = value;
                scale[6 + 6 * j + axis] = value.abs();
            }
        }
        Ok((values, scale))
    }

    fn endpoint_rate(
        &self,
        loss: &[f64],
        strain: &Voigt,
        q: &[Voigt],
        stress: &Voigt,
        tolerance: f64,
    ) -> Result<Voigt> {
        let rhs = core::array::from_fn(|i| -loss[i]);
        let rate = self.solve_viscosity_root(rhs, false)?;
        let kelvin = self
            .model
            .kv_pa_s()
            .ok_or_else(|| invalid("missing Kelvin tensor"))?;
        let eq = apply(&self.model.equilibrium_gram, strain);
        for axis in 0..6 {
            let mut actual = dd_dot(&kelvin[axis], &rate);
            let mut scale = stress[axis].abs();
            for (j, &velocity) in rate.iter().enumerate() {
                scale += (kelvin[axis][j] * velocity).abs();
            }
            for (row, &value) in eq.iter().enumerate() {
                let contribution =
                    Dd::from_f64(self.model.equilibrium_gram[row][axis]) * Dd::from_f64(value);
                actual = actual + contribution;
                scale += contribution.abs().to_f64();
            }
            for (j, branch) in self.model.branches().iter().enumerate() {
                let factor = branch.resolved_gram_rows(&self.model.factors[j]);
                for (row, &value) in q[j].iter().enumerate() {
                    let contribution = Dd::from_f64(factor[row][axis]) * Dd::from_f64(value);
                    actual = actual + contribution;
                    scale += contribution.abs().to_f64();
                }
            }
            let residual = (actual - Dd::from_f64(stress[axis])).to_f64();
            if !scale.is_finite() || !residual.is_finite() || residual.abs() > tolerance * scale {
                return Err(invalid(
                    "original Kelvin stress reconstruction is unresolved",
                ));
            }
        }
        Ok(rate)
    }

    /// Advance an exact stress ramp, hold, unloading or recovery interval.
    ///
    /// SPD Kelvin viscosity makes strain continuous: a zero-duration stress
    /// jump changes only stress and its compatible strain rate, with zero work
    /// or heat. Positive-duration ramps integrate continuous work and total
    /// Kelvin/Maxwell heat independently. Conservative whole-path strain bounds
    /// can refuse a path whose exact maximum would fit the material limit.
    ///
    /// # Errors
    /// Model mismatch, nonfinite or unrepresentable data, unresolved original
    /// source/energy balances, an exceeded strain bound, or cancellation. Every
    /// history field is unchanged on error.
    #[allow(clippy::too_many_lines)] // One source-bound atomic constitutive transaction.
    pub fn advance(
        &self,
        history: &mut TensorPronyHistory,
        stress: Voigt,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TensorKelvinStressRampStep> {
        poll(&mut cancelled)?;
        if !Arc::ptr_eq(&self.model, &history.model) {
            return Err(invalid(
                "Kelvin controller and history have different material bindings",
            ));
        }
        if !duration_s.is_finite() || duration_s < 0.0 || !finite(&stress) {
            return Err(invalid("finite stress and nonnegative duration required"));
        }
        let target = mandel(&stress);
        let delta = core::array::from_fn(|i| target[i] - history.stress[i]);
        if !finite(&target) || !finite(&delta) {
            return Err(invalid("unrepresentable prescribed stress increment"));
        }
        let n = self.rates.len();
        let tolerance = 8192.0 * f64::EPSILON * n as f64;
        if duration_s == 0.0 {
            let (loss, _) = self.original_loss(&history.strain, &history.q, &target)?;
            let rate =
                self.endpoint_rate(&loss, &history.strain, &history.q, &target, tolerance)?;
            poll(&mut cancelled)?;
            history.stress = target;
            return Ok(TensorKelvinStressRampStep {
                stress,
                strain: history.strain(),
                strain_increment: [0.0; 6],
                strain_rate: tensor(&rate),
                stored_energy: history.stored,
                stored_energy_change: 0.0,
                work: 0.0,
                dissipated: 0.0,
                balance_residual: 0.0,
                strain_bounds: [history.strain(); 2],
            });
        }
        let (initial_loss, initial_scale) =
            self.original_loss(&history.strain, &history.q, &history.stress)?;
        let mut responses = Vec::new();
        responses
            .try_reserve_exact(n)
            .map_err(|_| invalid("Kelvin ramp allocation refused"))?;
        let mut increment = [Dd::ZERO; 6];
        let mut first_increment = [Dd::ZERO; 6];
        let mut lower = [0.0; 6];
        let mut upper = [0.0; 6];
        let mut heat = Dd::ZERO;
        for mode in 0..n {
            poll(&mut cancelled)?;
            let initial = checked(
                initial_loss
                    .iter()
                    .enumerate()
                    .fold(Dd::ZERO, |sum, (row, &v)| {
                        sum + Dd::from_f64(self.left_modes[row * n + mode]) * Dd::from_f64(v)
                    })
                    .to_f64(),
                "unrepresentable Kelvin modal loss coordinate",
            )?;
            let forcing = checked(
                -dd_dot(&self.strain_modes[mode], &delta).to_f64(),
                "unrepresentable Kelvin modal stress forcing",
            )?;
            let rate = self.rates[mode];
            let r = duration_s * rate;
            if r.is_infinite() && forcing != 0.0 {
                return Err(invalid(
                    "unrepresentable Kelvin stress-ramp relaxation number",
                ));
            }
            let ramp = PronyRamp::new(r);
            let (end, mean, _) = ramp.apply(initial, forcing);
            let integrals = modal_integrals(initial, forcing, rate, duration_s, mean, &ramp);
            let mode_heat = ramp.squared_integral(initial, forcing, duration_s, rate);
            if !end.is_finite() || !mode_heat.is_finite() || mode_heat < 0.0 {
                return Err(invalid("unrepresentable Kelvin modal response or heat"));
            }
            heat = heat + Dd::from_f64(mode_heat);
            let response = Response {
                ramp,
                initial,
                forcing,
                rate,
                duration: duration_s,
                end,
                integrals,
            };
            for axis in 0..6 {
                let (zeroth, first, parts) =
                    response.scaled_integrals(-self.strain_modes[mode][axis])?;
                increment[axis] = increment[axis] + Dd::from_f64(zeroth);
                first_increment[axis] = first_increment[axis] + Dd::from_f64(first);
                for part in parts {
                    lower[axis] += part.min(0.0);
                    upper[axis] += part.max(0.0);
                }
            }
            responses.push(response);
        }
        let de = increment.map(Dd::to_f64);
        let strain = core::array::from_fn(|i| history.strain[i] + de[i]);
        let mut dq = vec![[0.0; 6]; history.q.len()];
        let mut q = history.q.clone();
        for row in 0..n - 6 {
            if row % 6 == 0 {
                poll(&mut cancelled)?;
            }
            let mut change = Dd::ZERO;
            for (mode, response) in responses.iter().enumerate() {
                change = change
                    + Dd::from_f64(
                        response
                            .scaled_integrals(-self.memory_modes[row * n + mode])?
                            .0,
                    );
            }
            dq[row / 6][row % 6] = change.to_f64();
            q[row / 6][row % 6] += change.to_f64();
        }
        let mut bounds = [[0.0; 6]; 2];
        let mut maximum = [0.0; 6];
        for axis in 0..6 {
            bounds[0][axis] = (history.strain[axis] + lower[axis]).min(strain[axis]);
            bounds[1][axis] = (history.strain[axis] + upper[axis]).max(strain[axis]);
            maximum[axis] = bounds[0][axis].abs().max(bounds[1][axis].abs());
        }
        self.model.validate_mandel_strain(&maximum).map_err(|_| {
            invalid("whole Kelvin stress-ramp path exceeds the material strain bound")
        })?;
        let mut work = Dd::ZERO;
        let mut exchange_scale = 0.0;
        for axis in 0..6 {
            let initial = Dd::from_f64(history.stress[axis]) * increment[axis];
            let change = Dd::from_f64(delta[axis]) * first_increment[axis];
            work = work + initial + change;
            exchange_scale += initial.abs().to_f64() + change.abs().to_f64();
        }
        let eq0 = apply(&self.model.equilibrium_gram, &history.strain);
        let eq_delta = apply(&self.model.equilibrium_gram, &de);
        let eq1 = apply(&self.model.equilibrium_gram, &strain);
        let mut storage_change = Dd::ZERO;
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
        let cumulative_work = history.work + work;
        let cumulative_heat = history.heat + heat;
        let scale = exchange_scale + heat;
        if !finite(&strain)
            || q.iter().any(|v| !finite(v))
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
                "Kelvin stress-control work/storage/heat balance is unresolved",
            ));
        }
        let (actual_loss, actual_scale) = self.original_loss(&strain, &q, &target)?;
        for row in 0..n {
            let mut predicted = Dd::ZERO;
            let mut scale = actual_scale[row] + initial_scale[row];
            for (mode, response) in responses.iter().enumerate() {
                let coefficient = self.left_modes[row * n + mode];
                predicted = predicted + Dd::from_f64(coefficient) * Dd::from_f64(response.end);
                scale += coefficient.abs()
                    * (response.end.abs() + response.initial.abs() + response.forcing.abs());
            }
            if !scale.is_finite()
                || (predicted - Dd::from_f64(actual_loss[row])).abs().to_f64() > tolerance * scale
            {
                return Err(invalid(
                    "original Kelvin/Maxwell loss reconstruction is unresolved",
                ));
            }
        }
        let rate = self.endpoint_rate(&actual_loss, &strain, &q, &target, tolerance)?;
        poll(&mut cancelled)?;
        let report = TensorKelvinStressRampStep {
            stress,
            strain: tensor(&strain),
            strain_increment: tensor(&de),
            strain_rate: tensor(&rate),
            stored_energy: stored,
            stored_energy_change: storage_change.to_f64(),
            work,
            dissipated: heat,
            balance_residual: balance,
            strain_bounds: [tensor(&bounds[0]), tensor(&bounds[1])],
        };
        history.strain = strain;
        history.stress = target;
        history.q = q;
        history.stored = stored;
        history.work = cumulative_work;
        history.heat = cumulative_heat;
        Ok(report)
    }
}
