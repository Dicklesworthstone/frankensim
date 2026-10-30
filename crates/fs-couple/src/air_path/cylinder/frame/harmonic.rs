//! Forced response of the SAME frozen frame and causal exterior-fluid states.
//!
//! Peak phasors use exp(+i omega t). Eliminating a fluid momentum gives
//! v_j=lambda_j/(lambda_j+i*omega) b^T*v. Only the structural dynamic stiffness
//! is factored; no branch, cross-mode term, material loss or added mass is fitted
//! or dropped. Every fluid momentum and its slip loss is reconstructed afterward.
//! This is a steady particular solution, NOT a warm-state time advance, acoustic
//! radiation model, continuum convergence certificate or live thermal update.

use super::MaterialFluidFrame;
use fs_exec::CancelGate;
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel};
use fs_math::c64::C64;

/// Separate limits for the reduced solve and its COMPLETE reconstructed result.
#[derive(Clone, Copy, Debug)]
pub struct FrameHarmonicBudget {
    /// Entries in each structural m-by-m matrix, not the original full state.
    pub max_modal_entries: usize,
    /// m^3 + m^2*(1 + material_ports + fluid_stations) + fluid_memories.
    /// Arithmetic work admission, not elapsed-time or total byte certification.
    pub max_setup_terms: usize,
    /// Full state, modal displacement/velocity and all physical loss channels.
    pub max_response_values: usize,
    /// Scaled infinity-norm dynamic-equation residual, strictly between 0 and 1.
    pub equation_relative_tolerance: f64,
    pub power_absolute_tolerance_w: f64,
    pub power_relative_tolerance: f64,
}

#[derive(Debug)]
pub enum FrameHarmonicError {
    Invalid(&'static str),
    Linear(LuComplexError),
    Equation { relative_residual: f64, tolerance: f64 },
    Power { residual_w: f64, tolerance_w: f64 },
    Cancelled,
}
impl core::fmt::Display for FrameHarmonicError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "frame/fluid harmonic response: {self:?}")
    }
}
impl std::error::Error for FrameHarmonicError {}
impl From<LuComplexError> for FrameHarmonicError {
    fn from(error: LuComplexError) -> Self {
        match error {
            LuComplexError::Cancelled => Self::Cancelled,
            other => Self::Linear(other),
        }
    }
}
type Result<T> = core::result::Result<T, FrameHarmonicError>;
fn invalid(message: &'static str) -> FrameHarmonicError { FrameHarmonicError::Invalid(message) }
fn poll(gate: &CancelGate) -> Result<()> {
    if gate.is_requested() { Err(FrameHarmonicError::Cancelled) } else { Ok(()) }
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() { Ok(value) } else { Err(invalid("unrepresentable harmonic arithmetic")) }
}
fn complex(value: C64) -> Result<C64> {
    finite(value.re)?;
    finite(value.im)?;
    Ok(value)
}
fn magnitude(value: C64) -> Result<f64> { finite(value.re.hypot(value.im)) }
fn squared(value: C64) -> Result<f64> {
    let norm = magnitude(value)?;
    let result = finite(norm * norm)?;
    if norm != 0.0 && result == 0.0 { return Err(invalid("harmonic quadratic account underflow")); }
    Ok(result)
}
fn projected(column: &[f64], values: &[C64]) -> Result<C64> {
    let mut value = C64::ZERO;
    for (&b, &x) in column.iter().zip(values) { value = complex(value + x.scale(b))?; }
    Ok(value)
}

/// All amplitudes are peak (not RMS). State order matches the original model:
/// [dry modal q, loaded modal momentum M*v, every original fluid momentum].
/// Losses are cycle-mean watts and stored energy is cycle-mean joules.
#[derive(Clone, Debug)]
pub struct FrameHarmonicResponse {
    pub omega_rad_s: f64,
    pub modal_displacement: Vec<C64>,
    pub modal_velocity: Vec<C64>,
    pub state: Vec<C64>,
    /// Original material-damping port order; never silently diagonalized.
    pub material_port_loss_w: Vec<f64>,
    /// Original loaded-region order, even when element IDs repeat.
    pub fluid_region_loss_w: Vec<f64>,
    pub supplied_power_w: f64,
    pub mean_stored_energy_j: f64,
    pub equation_relative_residual: f64,
    pub power_residual_w: f64,
}

/// One immutable source model and one reusable factorization at one frequency.
/// Preparing/solving never evaluates a replacement material or changes history.
/// Frozen source/geometry, modal truncation and spectral-cutoff limits are those
/// of the original model. A small residual is NOT a condition-number estimate.
pub struct PreparedFrameHarmonic<'m, 'f, F> {
    model: &'m MaterialFluidFrame<'f, F>,
    omega: f64,
    dynamic: Vec<C64>,
    factor: LuComplex,
    follow: Vec<C64>,
    slip: Vec<C64>,
    budget: FrameHarmonicBudget,
}

impl<'f, F> MaterialFluidFrame<'f, F> {
    /// Prepare forced vibration without factoring the complete fluid-memory
    /// state. The full loaded mass is included ONCE; the remaining impedance
    /// contains only the original low tail and causal memory resistances.
    ///
    /// DC and unsupported solid-viscosity bands refuse. This does not claim
    /// that every structural mode or every continuum frequency is resolved.
    pub fn prepare_harmonic(&self, omega: f64, budget: FrameHarmonicBudget,
        gate: &CancelGate) -> Result<PreparedFrameHarmonic<'_, 'f, F>> {
        poll(gate)?;
        if !omega.is_finite() || omega <= 0.0
            || !budget.equation_relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&budget.equation_relative_tolerance)
            || budget.equation_relative_tolerance == 0.0
            || !budget.power_absolute_tolerance_w.is_finite() || budget.power_absolute_tolerance_w <= 0.0
            || !budget.power_relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&budget.power_relative_tolerance) {
            return Err(invalid("positive frequency and explicit finite equation/power tolerances required"));
        }
        let omega2 = finite(omega * omega)?;
        if omega2 == 0.0 { return Err(invalid("squared harmonic frequency underflow")); }
        let m = self.basis.modes().len();
        let k = self.memory_masses.len();
        let square = m.checked_mul(m).ok_or_else(|| invalid("harmonic matrix size overflow"))?;
        let terms = self.material_ports.len().checked_add(self.stations.len())
            .and_then(|n| n.checked_add(1)).and_then(|n| n.checked_add(m))
            .and_then(|n| n.checked_mul(square)).and_then(|n| n.checked_add(k));
        let values = m.checked_mul(2).and_then(|n| n.checked_add(self.system.state_dim()))
            .and_then(|n| n.checked_add(self.material_ports.len()))
            .and_then(|n| n.checked_add(self.members.len()));
        if square > budget.max_modal_entries
            || terms.is_none_or(|n| n > budget.max_setup_terms)
            || values.is_none_or(|n| n > budget.max_response_values) {
            return Err(invalid("complete harmonic setup or response exceeds its budget"));
        }
        if let Some(damping) = self.basis.material_damping() {
            for (_, section) in damping.sections() {
                poll(gate)?;
                let (low, high) = section.omega_band_rad_s();
                if omega < low || omega > high {
                    return Err(invalid("harmonic frequency is outside original material-viscosity support"));
                }
            }
        }
        let mut dynamic = vec![C64::ZERO; square];
        for i in 0..m {
            poll(gate)?;
            for j in 0..m {
                dynamic[i*m+j].re = finite(-omega2 * self.modal_mass[i*m+j])?;
            }
            dynamic[i*m+i].re = finite(dynamic[i*m+i].re + self.basis.modes()[i].eigenvalue_per_s2)?;
        }
        for port in &self.material_ports {
            poll(gate)?;
            outer(&mut dynamic, m, &port.projection, C64::new(0.0, finite(omega*port.resistance)?))?;
        }
        let mut follow = vec![C64::ZERO; k];
        let mut slip = vec![C64::ZERO; k];
        for station in &self.stations {
            poll(gate)?;
            let flow = &self.members[station.member].flow;
            let mut impedance = C64::new(finite(station.weight_m*flow.low_tail_resistance_per_length_n_s_m2())?, 0.0);
            for (j, branch) in flow.branches().iter().enumerate() {
                poll(gate)?;
                let index = station.first_memory+j;
                let c = finite(station.weight_m*branch.resistance_per_length_n_s_m2)?;
                let rate = finite(c/self.memory_masses[index])?;
                if c <= 0.0 || rate <= 0.0 { return Err(invalid("unrepresentable physical fluid resistance/rate")); }
                let (a, b) = memory_transfer(rate, omega)?;
                follow[index] = a;
                slip[index] = b;
                impedance = complex(impedance + b.scale(c))?;
            }
            outer(&mut dynamic, m, &station.projection,
                complex(C64::new(-omega*impedance.im, omega*impedance.re))?)?;
        }
        let factor = lu_complex_with_cancel(&dynamic, m, || gate.is_requested())?;
        poll(gate)?;
        Ok(PreparedFrameHarmonic { model: self, omega, dynamic, factor, follow, slip, budget })
    }
}

impl<'m, 'f, F> PreparedFrameHarmonic<'m, 'f, F> {
    pub const fn model(&self) -> &'m MaterialFluidFrame<'f, F> { self.model }
    pub const fn omega_rad_s(&self) -> f64 { self.omega }
    /// Original-basis K - omega^2*M + i*omega*C with EXACT branch elimination.
    pub fn dynamic_stiffness(&self) -> &[C64] { &self.dynamic }
    pub fn factor_dimension(&self) -> usize { self.model.basis.modes().len() }

    /// Solve a complete original-basis generalized-force phasor. Reuses this
    /// frequency's LU for arbitrary phases/loads, then independently accounts
    /// each physical slip loss. No response escapes cancellation or refusal.
    pub fn solve_modal(&self, force: &[C64], gate: &CancelGate) -> Result<FrameHarmonicResponse> {
        poll(gate)?;
        let m = self.factor_dimension();
        if force.len() != m { return Err(invalid("one force phasor per original structural mode required")); }
        for &f in force { complex(f)?; }
        let mut q = force.to_vec();
        self.factor.solve_with_cancel(&mut q, || gate.is_requested())?;
        for &x in &q { complex(x)?; }
        let residual = equation_residual(&self.dynamic, &q, force, gate)?;
        if residual > self.budget.equation_relative_tolerance {
            return Err(FrameHarmonicError::Equation { relative_residual: residual,
                tolerance: self.budget.equation_relative_tolerance });
        }
        let v = q.iter().map(|q| complex(C64::new(-self.omega*q.im, self.omega*q.re)))
            .collect::<Result<Vec<_>>>()?;
        let mut state = vec![C64::ZERO; self.model.system.state_dim()];
        state[..m].copy_from_slice(&q);
        let mut stored = 0.0;
        for i in 0..m {
            poll(gate)?;
            state[m+i] = projected(&self.model.modal_mass[i*m..(i+1)*m], &v)?;
            stored = finite(stored + 0.25*self.model.basis.modes()[i].eigenvalue_per_s2*squared(q[i])?
                + 0.25*complex(v[i].conj()*state[m+i])?.re)?;
        }
        let mut material_loss = Vec::with_capacity(self.model.material_ports.len());
        for port in &self.model.material_ports {
            poll(gate)?;
            material_loss.push(finite(0.5*port.resistance*squared(projected(&port.projection, &v)?)?)?);
        }
        let mut fluid_loss = vec![0.0; self.model.members.len()];
        for station in &self.model.stations {
            poll(gate)?;
            let flow = &self.model.members[station.member].flow;
            let speed = projected(&station.projection, &v)?;
            let mut loss = finite(0.5*station.weight_m*flow.low_tail_resistance_per_length_n_s_m2()*squared(speed)?)?;
            for (j, branch) in flow.branches().iter().enumerate() {
                poll(gate)?;
                let index = station.first_memory+j;
                let fluid_v = complex(self.follow[index]*speed)?;
                // Direct i*omega/(lambda+i*omega) avoids subtracting nearly
                // equal velocities when a fast memory follows the body.
                let slip_v = complex(self.slip[index]*speed)?;
                let mass = self.model.memory_masses[index];
                state[2*m+index] = complex(fluid_v.scale(mass))?;
                stored = finite(stored + 0.25*mass*squared(fluid_v)?)?;
                loss = finite(loss + 0.5*station.weight_m*branch.resistance_per_length_n_s_m2*squared(slip_v)?)?;
            }
            fluid_loss[station.member] = finite(fluid_loss[station.member]+loss)?;
        }
        let mut supplied = 0.0;
        for (&f, &velocity) in force.iter().zip(&v) {
            supplied = finite(supplied + 0.5*complex(f*velocity.conj())?.re)?;
        }
        let dissipated = material_loss.iter().chain(&fluid_loss)
            .try_fold(0.0, |sum, value| finite(sum+value))?;
        let power_residual = finite(supplied-dissipated)?;
        let tolerance = finite(self.budget.power_absolute_tolerance_w
            + self.budget.power_relative_tolerance*supplied.abs().max(dissipated))?;
        if power_residual.abs() > tolerance {
            return Err(FrameHarmonicError::Power { residual_w: power_residual, tolerance_w: tolerance });
        }
        if stored < 0.0 { return Err(invalid("negative harmonic stored-energy account")); }
        poll(gate)?;
        Ok(FrameHarmonicResponse { omega_rad_s: self.omega, modal_displacement: q,
            modal_velocity: v, state, material_port_loss_w: material_loss,
            fluid_region_loss_w: fluid_loss, supplied_power_w: supplied,
            mean_stored_energy_j: stored, equation_relative_residual: residual,
            power_residual_w: power_residual })
    }
}

fn memory_transfer(rate: f64, omega: f64) -> Result<(C64, C64)> {
    let scale = rate.max(omega);
    let a = rate/scale;
    let b = omega/scale;
    if a == 0.0 || b == 0.0 { return Err(invalid("unrepresentable fluid-memory frequency ratio")); }
    let denominator = a*a+b*b;
    Ok((complex(C64::new(a*a/denominator, -a*b/denominator))?,
        complex(C64::new(b*b/denominator, a*b/denominator))?))
}
fn outer(matrix: &mut [C64], m: usize, column: &[f64], coefficient: C64) -> Result<()> {
    for i in 0..m {
        for j in 0..=i {
            let value = complex(coefficient.scale(finite(column[i]*column[j])?))?;
            let next = complex(matrix[i*m+j]+value)?;
            matrix[i*m+j] = next;
            matrix[j*m+i] = next;
        }
    }
    Ok(())
}
fn equation_residual(matrix: &[C64], q: &[C64], force: &[C64], gate: &CancelGate) -> Result<f64> {
    let m = q.len();
    let mut matrix_norm = 0.0_f64;
    let mut residual = 0.0_f64;
    let mut q_norm = 0.0_f64;
    let mut f_norm = 0.0_f64;
    for i in 0..m {
        poll(gate)?;
        let mut action = C64::ZERO;
        let mut row_norm = 0.0;
        for j in 0..m {
            action = complex(action+matrix[i*m+j]*q[j])?;
            row_norm = finite(row_norm+magnitude(matrix[i*m+j])?)?;
        }
        matrix_norm = matrix_norm.max(row_norm);
        residual = residual.max(magnitude(complex(action-force[i])?)?);
        q_norm = q_norm.max(magnitude(q[i])?);
        f_norm = f_norm.max(magnitude(force[i])?);
    }
    let scale = finite(matrix_norm*q_norm+f_norm)?;
    if scale == 0.0 {
        if residual != 0.0 { return Err(invalid("nonzero equation residual without a scale")); }
        Ok(0.0)
    } else { finite(residual/scale) }
}

#[cfg(test)]
mod tests;
