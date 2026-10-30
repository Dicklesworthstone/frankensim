//! Continuous forced response of the original coupled linear system.
//! Peak phasors use exp(-i omega t); averages therefore carry a factor 1/2.
use super::*;
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel};
use fs_math::{c64::C64, det};

pub mod dg;

/// Dense response capacity; no retained mode or acoustic pole is truncated.
pub const MAX_HARMONIC_UNKNOWNS: usize = 1024;
/// A batch retains complete responses, including every material/air history.
pub const MAX_HARMONIC_RIGHT_HAND_SIDES: usize = 64;
/// Retained response vectors and their struct headers, checked before assembly.
/// Operator/factor workspace remains bounded by MAX_HARMONIC_UNKNOWNS and the
/// original construction budget; caller-owned input vectors are not included.
pub const MAX_HARMONIC_BATCH_RESPONSE_BYTES: usize = 8 * 1_048_576;

/// AC particular solution, independent of accepted histories and DC preload.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalHarmonicResponse {
    pub omega_rad_s: f64,
    /// Original flattened modal coordinates, including any internal components.
    pub displacement_m_sqrt_kg: Vec<C64>,
    pub velocity_m_sqrt_kg_per_s: Vec<C64>,
    /// Energy-normalized Maxwell internal strains, in original branch order.
    pub material_memory: Vec<C64>,
    /// Acoustic displacement/momentum pairs, in original radiation-pole order.
    pub radiation_state: Vec<C64>,
    /// Sum of each connection's signed spring and damping efforts [N].
    /// For an ordinary link both act along its original column. A prepared
    /// constant-loss link has DISTINCT spring and damping columns; this sum
    /// is then only a channel diagnostic, not a force along either column.
    pub connection_forces_n: Vec<C64>,
    /// Spring effort `-k b^T q`, in original connection order [N].
    pub connection_spring_forces_n: Vec<C64>,
    /// Dissipative effort `-c d^T v`, in original connection order [N].
    /// `d = b` for Kelvin-Voigt; constant-loss links use their prepared `g(K)b`.
    pub connection_damping_forces_n: Vec<C64>,
    /// Sum of `b * spring_effort + d * damping_effort` on each ORIGINAL
    /// mass-normalized mode [N/sqrt(kg)], retaining both reaction directions.
    pub generalized_connection_forces_n_per_sqrt_kg: Vec<C64>,
    /// Mean quadratic storage of the AC solution only [J].
    pub mean_stored_energy_j: f64,
    /// Mean work rate of the supplied peak-force phasors [W].
    pub external_power_w: f64,
    pub component_dissipation_w: f64,
    pub connection_dissipation_w: f64,
    pub material_dissipation_w: f64,
    /// Mean heat rate of each original Maxwell branch [W], in the same order
    /// as `material_memory`. Source adapters may sum their retained ranges;
    /// a factor that combines several source films is not an individual film.
    pub material_branch_dissipation_w: Vec<f64>,
    /// Mean stored strain energy of each original Maxwell branch [J].
    /// For peak mismatch r these are |r|²/4 and heat is |r|²/(2 tau).
    pub material_branch_mean_storage_j: Vec<f64>,
    pub radiation_dissipation_w: f64,
    /// Total physical mean loss minus external mean power [W].
    pub power_residual_w: f64,
    /// Infinity-norm backward residual of the ORIGINAL equations.
    pub solve_relative_residual: f64,
}
impl ModalHarmonicResponse {
    pub fn total_dissipation_w(&self) -> f64 {
        self.component_dissipation_w
            + self.connection_dissipation_w
            + self.material_dissipation_w
            + self.radiation_dissipation_w
    }
}

fn complex_finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn lu_error(error: LuComplexError) -> ModalCouplingError {
    match error {
        LuComplexError::Cancelled => ModalCouplingError::Cancelled,
        LuComplexError::Singular { .. } => invalid("harmonic dynamic stiffness is singular"),
    }
}
fn projection(column: &[f64], q: &[C64]) -> C64 {
    column
        .iter()
        .zip(q)
        .fold(C64::ZERO, |sum, (&b, &x)| sum + x.scale(b))
}
fn add_rank(
    matrix: &mut [C64],
    dimension: usize,
    column: &[f64],
    gain: C64,
    support: &mut Vec<(usize, f64)>,
    gate: &CancelGate,
) -> Result<(), ModalCouplingError> {
    poll(Some(gate))?;
    if !complex_finite(gain) {
        return Err(invalid("harmonic operator overflow"));
    }
    // A material or connection may touch only one body inside a much larger
    // coupled instrument. Keep every nonzero coefficient in its original order
    // and reuse this support buffer across ranks; exact zeros need no matrix
    // update. No threshold, retained coordinate or history is removed.
    support.clear();
    for (index, &value) in column.iter().enumerate() {
        if !value.is_finite() {
            return Err(invalid("harmonic rank column is not finite"));
        }
        if value != 0.0 {
            support.push((index, value));
        }
    }
    for &(i, a) in support.iter() {
        poll(Some(gate))?;
        for &(j, b) in support.iter() {
            matrix[i * dimension + j] = matrix[i * dimension + j] + gain.scale(a * b);
        }
    }
    Ok(())
}

// The original mechanical/air operator may reserve additional trailing states
// for a source-owned continuous load. Every original row keeps its same basis.
struct HarmonicOperator {
    values: Vec<C64>,
    dimension: usize,
    original_dimension: usize,
}

/// Borrow the original equations while sharing only their scaled LU. The
/// borrowing lifetime prevents a factor from outliving its source operator.
struct HarmonicFactor<'a> {
    matrix: &'a [C64],
    dimension: usize,
    scale: f64,
    factor: LuComplex,
}
impl<'a> HarmonicFactor<'a> {
    fn new(matrix: &'a [C64], dimension: usize, gate: &CancelGate)
        -> Result<Self, ModalCouplingError>
    {
        poll(Some(gate))?;
        if matrix.iter().any(|&a| !complex_finite(a)) {
            return Err(invalid("harmonic operator overflow"));
        }
        let scale = matrix.iter().map(|a| a.abs()).fold(0.0_f64, f64::max);
        if !scale.is_finite() || scale == 0.0 {
            return Err(invalid("harmonic operator has no finite scale"));
        }
        let scaled: Vec<_> = matrix.iter()
            .map(|a| C64::new(a.re / scale, a.im / scale)).collect();
        let factor = lu_complex_with_cancel(&scaled, dimension, || gate.is_requested())
            .map_err(lu_error)?;
        Ok(Self { matrix, dimension, scale, factor })
    }

    fn solve(&self, force: &[C64], tolerance: f64, gate: &CancelGate)
        -> Result<(Vec<C64>, f64), ModalCouplingError>
    {
        poll(Some(gate))?;
        let Self { matrix, dimension, scale, factor } = self;
        let dimension = *dimension;
        let mut state = vec![C64::ZERO; dimension];
        for (x, f) in state.iter_mut().zip(force) {
            *x = C64::new(f.re / scale, f.im / scale);
        }
        if state.iter().any(|&x| !complex_finite(x)) {
            return Err(invalid("harmonic force scaling overflow"));
        }
        factor.solve_with_cancel(&mut state, || gate.is_requested()).map_err(lu_error)?;
        if state.iter().any(|&x| !complex_finite(x)) {
            return Err(invalid("harmonic response is not representable"));
        }
        let mut residual_norm = 0.0_f64;
        let mut operator_norm = 0.0_f64;
        let state_norm = state.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        let force_norm = force.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        for i in 0..dimension {
            poll(Some(gate))?;
            let rhs = force.get(i).copied().unwrap_or(C64::ZERO);
            let mut residual = -rhs;
            let mut row_norm = 0.0;
            for (j, &x) in state.iter().enumerate() {
                residual = residual + matrix[i * dimension + j] * x;
                row_norm += matrix[i * dimension + j].abs();
            }
            if !complex_finite(residual) || !row_norm.is_finite() {
                return Err(invalid("harmonic residual overflow"));
            }
            residual_norm = residual_norm.max(residual.abs());
            operator_norm = operator_norm.max(row_norm);
        }
        let denominator = operator_norm * state_norm + force_norm;
        if !denominator.is_finite() {
            return Err(invalid("harmonic residual scale overflow"));
        }
        let relative = if denominator == 0.0 { 0.0 } else { residual_norm / denominator };
        if relative > tolerance {
            return Err(ModalCouplingError::SolveResidual { relative, tolerance });
        }
        Ok((state, relative))
    }
}

fn solve_harmonic(matrix: &[C64], dimension: usize, force: &[C64], tolerance: f64,
    gate: &CancelGate) -> Result<(Vec<C64>, f64), ModalCouplingError>
{
    HarmonicFactor::new(matrix, dimension, gate)?.solve(force, tolerance, gate)
}

fn harmonic_work(dimension: usize, connections: usize, constant_loss: usize,
    memories: usize, allowance: usize, responses: usize) -> Result<usize, ModalCouplingError>
{
    // The first response keeps the original complete-operator allowance. Each
    // later response pays for triangular solves, residuals and all physical
    // observations; only factorization/assembly are shared. Full-support costs
    // remain conservative even when rank assembly touches a small local body.
    let per_response = connections.checked_add(memories)
        .and_then(|v| v.checked_add(allowance));
    dimension.checked_mul(dimension)
        .and_then(|square| {
            per_response.and_then(|per| per.checked_mul(responses))
                .and_then(|v| v.checked_add(dimension))
                .and_then(|v| v.checked_add(constant_loss))
                .and_then(|v| square.checked_mul(v))
        })
        .ok_or_else(|| invalid("harmonic work overflow"))
}

fn harmonic_response_bytes(modes: usize, connections: usize, memories: usize,
    poles: usize, responses: usize) -> Result<usize, ModalCouplingError>
{
    // Count every retained output array: q, v and generalized reaction; three
    // connection efforts; Maxwell histories and heat/storage; paired air states.
    // Struct size includes Vec headers, scalar diagnostics and alignment.
    modes.checked_mul(3)
        .and_then(|v| v.checked_add(connections.checked_mul(3)?))
        .and_then(|v| v.checked_add(memories))
        .and_then(|v| v.checked_add(poles.checked_mul(2)?))
        .and_then(|v| v.checked_mul(core::mem::size_of::<C64>()))
        .and_then(|v| v.checked_add(memories.checked_mul(2 * core::mem::size_of::<f64>())?))
        .and_then(|v| v.checked_add(core::mem::size_of::<ModalHarmonicResponse>()))
        .and_then(|v| v.checked_mul(responses))
        .ok_or_else(|| invalid("harmonic response storage overflow"))
}

impl CoupledModalSystem {
    /// Continuous linear response to peak generalized-force phasors in
    /// exp(-i omega t), with ALL original connections, Maxwell arms and air.
    /// This is the AC particular solution: initial transients, DC rest lengths,
    /// nonlinear contacts and the discrete held-force transfer are not included.
    /// Material/fit frequency validity remains with the coefficient producer.
    ///
    /// Maxwell states are eliminated algebraically, without approximating their
    /// frequency dependence. Acoustic displacement unknowns remain explicit, so
    /// a bare lossless air pole can produce its physical antiresonance rather
    /// than a spurious refusal from dividing by an isolated-pole impedance.
    ///
    /// `maximum_work` bounds d^3+d^2*(connections+constant_loss_links+material_arms+3), with
    /// d=mechanical_modes+acoustic_poles, before allocation; d is also capped at
    /// MAX_HARMONIC_UNKNOWNS. The existing solve-relative tolerance gates the
    /// original equations. Constant-loss connections reuse the prepared owner's
    /// SAME coupled-basis damping columns; the spring retains its original column.
    /// Rank assembly visits only exact nonzero column support, while the stated
    /// work allowance remains a conservative full-support bound.
    /// This is the continuous response of that fixed finite-dimensional image,
    /// not a new frequency-dependent loss law evaluated at the drive frequency.
    /// No accepted state, material history, sample clock or ledger is changed.
    pub fn harmonic_response(
        &self,
        omega_rad_s: f64,
        force: &[C64],
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<ModalHarmonicResponse, ModalCouplingError> {
        let operator = self.harmonic_operator(omega_rad_s, force, 0, maximum_work, gate)?;
        let (state, relative) = solve_harmonic(
            &operator.values,
            operator.dimension,
            force,
            self.config.solve_relative_tolerance,
            gate,
        )?;
        self.observe_harmonic_state(omega_rad_s, force, &state, relative, gate)
    }

    /// Complete responses to several independent peak-force vectors at ONE
    /// frequency, in the caller's order. The original operator is assembled and
    /// factored once; each RHS retains its own original-equation residual,
    /// Maxwell/air histories and full physical power/storage accounting.
    ///
    /// For d unknowns, C connections, G anchored-loss links, M Maxwell branches
    /// and r force vectors, `maximum_work` bounds d²*(d+G+r*(C+M+3)). This is the
    /// existing single-RHS allowance when r=1. Retained result storage is checked
    /// before assembly against MAX_HARMONIC_BATCH_RESPONSE_BYTES; 1..=64 RHS are
    /// admitted. Input vectors contain exactly one finite force per original
    /// mechanical mode. No partial result or accepted-state change is published.
    ///
    /// Responses are independent load cases, NOT incoherent power components.
    /// For coherent excitation, sum force phasors before solving: powers and
    /// stored energies are quadratic and must not be linearly superposed.
    pub fn harmonic_responses(
        &self, omega_rad_s: f64, forces: &[&[C64]], maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<Vec<ModalHarmonicResponse>, ModalCouplingError> {
        poll(Some(gate))?;
        if forces.is_empty() || forces.len() > MAX_HARMONIC_RIGHT_HAND_SIDES {
            return Err(invalid("harmonic batch requires 1..=64 complete force vectors"));
        }
        let n = self.mode_count();
        let poles = self.radiation.as_ref().map_or(0, |r| r.model.poles.len());
        let dimension = n.checked_add(poles).ok_or_else(|| invalid("harmonic dimension overflow"))?;
        let work = harmonic_work(dimension, self.connections.len(),
            self.connections.iter().filter(|c| c.constant_loss_reference_omega_rad_s.is_some()).count(),
            self.relaxation_branches.len(), 3, forces.len())?;
        let bytes = harmonic_response_bytes(n, self.connections.len(),
            self.relaxation_branches.len(), poles, forces.len())?;
        if dimension > MAX_HARMONIC_UNKNOWNS || work > maximum_work
            || bytes > MAX_HARMONIC_BATCH_RESPONSE_BYTES
        {
            return Err(invalid("complete harmonic batch exceeds its work or retained-response storage bound"));
        }
        for force in forces {
            poll(Some(gate))?;
            if force.len() != n || force.iter().any(|&f| !complex_finite(f)) {
                return Err(invalid("harmonic batch requires one finite force per original mode in every RHS"));
            }
        }
        let operator = self.harmonic_operator(omega_rad_s, forces[0], 0, maximum_work, gate)?;
        let factor = HarmonicFactor::new(&operator.values, operator.dimension, gate)?;
        let mut responses = Vec::with_capacity(forces.len());
        for &force in forces {
            let (state, relative) = factor.solve(force, self.config.solve_relative_tolerance, gate)?;
            responses.push(self.observe_harmonic_state(omega_rad_s, force, &state, relative, gate)?);
        }
        poll(Some(gate))?;
        Ok(responses)
    }

    fn harmonic_operator(
        &self,
        omega_rad_s: f64,
        force: &[C64],
        additional: usize,
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<HarmonicOperator, ModalCouplingError> {
        poll(Some(gate))?;
        if !self.plane_stress_material_points.is_empty() {
            return Err(invalid(
                "harmonic response does not support bound local plane-stress strain checks",
            ));
        }
        let n = self.mode_count();
        let poles = self.radiation.as_ref().map_or(0, |r| r.model.poles.len());
        let original_dimension = n
            .checked_add(poles)
            .ok_or_else(|| invalid("harmonic dimension overflow"))?;
        let dimension = original_dimension
            .checked_add(additional)
            .ok_or_else(|| invalid("harmonic dimension overflow"))?;
        let square = dimension
            .checked_mul(dimension)
            .ok_or_else(|| invalid("harmonic matrix overflow"))?;
        let work = harmonic_work(dimension, self.connections.len(),
            self.connections.iter().filter(|c| c.constant_loss_reference_omega_rad_s.is_some()).count(),
            self.relaxation_branches.len(), if additional == 0 { 3 } else { 5 }, 1)?;
        if dimension > MAX_HARMONIC_UNKNOWNS || work > maximum_work {
            return Err(invalid(
                "complete harmonic solve exceeds its explicit dense work budget",
            ));
        }
        let omega_sq = omega_rad_s * omega_rad_s;
        if !omega_rad_s.is_finite()
            || omega_rad_s <= 0.0
            || !omega_sq.is_finite()
            || omega_sq == 0.0
            || force.len() != n
            || force.iter().any(|&f| !complex_finite(f))
        {
            return Err(invalid(
                "harmonic response requires positive representable frequency and one finite force per mode",
            ));
        }
        let s = C64::new(0.0, -omega_rad_s);
        let modes: Vec<_> = self.models.iter().flat_map(|m| m.modes()).collect();
        let mut matrix = vec![C64::ZERO; square];
        let mut support = Vec::with_capacity(n);
        for (i, mode) in modes.iter().enumerate() {
            poll(Some(gate))?;
            let w = mode.angular_frequency_rad_s;
            matrix[i * dimension + i] = C64::new(
                w * w - omega_sq,
                -omega_rad_s * (2.0 * mode.damping_ratio * w),
            );
        }
        for (index, (column, law)) in self.columns.iter().zip(&self.connections).enumerate() {
            if law.constant_loss_reference_omega_rad_s.is_some() {
                let damping = self.harmonic_damping_column(index)?;
                add_rank(
                    &mut matrix,
                    dimension,
                    column,
                    C64::from_re(law.stiffness_n_m),
                    &mut support,
                    gate,
                )?;
                add_rank(
                    &mut matrix,
                    dimension,
                    damping,
                    C64::new(0.0, -omega_rad_s * law.damping_n_s_m),
                    &mut support,
                    gate,
                )?;
            } else {
                // Keep the original single-rank arithmetic for ordinary links.
                add_rank(
                    &mut matrix,
                    dimension,
                    column,
                    C64::new(law.stiffness_n_m, -omega_rad_s * law.damping_n_s_m),
                    &mut support,
                    gate,
                )?;
            }
        }
        let mut column = vec![0.0; n];
        for arm in &self.relaxation_branches {
            let root = det::sqrt(arm.stiffness);
            for (b, &p) in column.iter_mut().zip(arm.projection.iter().step_by(2)) {
                *b = root * p;
            }
            let gain = s / (s + C64::from_re(arm.relaxation_time_s.recip()));
            add_rank(&mut matrix, dimension, &column, gain, &mut support, gate)?;
        }
        if let Some(air) = &self.radiation {
            for (j, pole) in air.model.poles.iter().enumerate() {
                poll(Some(gate))?;
                let index = n + j;
                matrix[index * dimension + index] = C64::new(
                    pole.omega * pole.omega - omega_sq,
                    -omega_rad_s * (2.0 * pole.zeta * pole.omega),
                );
                for (&i, &g) in air.source_modes.iter().zip(&pole.coupling) {
                    matrix[i * dimension + index] = s.scale(g);
                    matrix[index * dimension + i] = -s.scale(g);
                }
            }
        }
        Ok(HarmonicOperator {
            values: matrix,
            dimension,
            original_dimension,
        })
    }

    fn harmonic_damping_column(&self, index: usize) -> Result<&[f64], ModalCouplingError> {
        if self.connections[index]
            .constant_loss_reference_omega_rad_s
            .is_some()
        {
            self.connection_damping_column(index)
                .ok_or_else(|| {
                    invalid("constant-loss harmonic response requires its prepared damping column")
                })
        } else {
            Ok(&self.columns[index])
        }
    }

    fn observe_harmonic_state(
        &self,
        omega_rad_s: f64,
        force: &[C64],
        state: &[C64],
        relative: f64,
        gate: &CancelGate,
    ) -> Result<ModalHarmonicResponse, ModalCouplingError> {
        let n = self.mode_count();
        let poles = self.radiation.as_ref().map_or(0, |r| r.model.poles.len());
        let s = C64::new(0.0, -omega_rad_s);
        let modes: Vec<_> = self.models.iter().flat_map(|m| m.modes()).collect();
        let mut column = vec![0.0; n];
        let mut result = ModalHarmonicResponse {
            omega_rad_s,
            displacement_m_sqrt_kg: state[..n].to_vec(),
            velocity_m_sqrt_kg_per_s: state[..n].iter().map(|&q| s * q).collect(),
            material_memory: Vec::with_capacity(self.relaxation_branches.len()),
            radiation_state: Vec::with_capacity(2 * poles),
            connection_forces_n: Vec::with_capacity(self.connections.len()),
            connection_spring_forces_n: Vec::with_capacity(self.connections.len()),
            connection_damping_forces_n: Vec::with_capacity(self.connections.len()),
            generalized_connection_forces_n_per_sqrt_kg: vec![C64::ZERO; n],
            mean_stored_energy_j: 0.0,
            external_power_w: 0.0,
            component_dissipation_w: 0.0,
            connection_dissipation_w: 0.0,
            material_dissipation_w: 0.0,
            material_branch_dissipation_w: Vec::with_capacity(self.relaxation_branches.len()),
            material_branch_mean_storage_j: Vec::with_capacity(self.relaxation_branches.len()),
            radiation_dissipation_w: 0.0,
            power_residual_w: 0.0,
            solve_relative_residual: relative,
        };
        for (i, mode) in modes.iter().enumerate() {
            let q = state[i];
            let v = result.velocity_m_sqrt_kg_per_s[i];
            let w = mode.angular_frequency_rad_s;
            result.mean_stored_energy_j += 0.25 * (q.scale(w).norm_sq() + v.norm_sq());
            result.component_dissipation_w += mode.damping_ratio * w * v.norm_sq();
            result.external_power_w += 0.5 * (force[i] * v.conj()).re;
        }
        for (index, (column, law)) in self.columns.iter().zip(&self.connections).enumerate() {
            poll(Some(gate))?;
            let damping_column = self.harmonic_damping_column(index)?;
            let extension = projection(column, &state[..n]);
            let velocity = if law.constant_loss_reference_omega_rad_s.is_some() {
                s * projection(damping_column, &state[..n])
            } else {
                s * extension
            };
            let spring = -extension.scale(law.stiffness_n_m);
            let damping = -velocity.scale(law.damping_n_s_m);
            result
                .connection_forces_n
                .push(-(extension.scale(law.stiffness_n_m) + velocity.scale(law.damping_n_s_m)));
            result.connection_spring_forces_n.push(spring);
            result.connection_damping_forces_n.push(damping);
            for ((force, &b), &d) in result
                .generalized_connection_forces_n_per_sqrt_kg
                .iter_mut()
                .zip(column)
                .zip(damping_column)
            {
                *force = *force + spring.scale(b) + damping.scale(d);
            }
            result.mean_stored_energy_j += 0.25 * law.stiffness_n_m * extension.norm_sq();
            result.connection_dissipation_w += 0.5 * law.damping_n_s_m * velocity.norm_sq();
        }
        for arm in &self.relaxation_branches {
            poll(Some(gate))?;
            let root = det::sqrt(arm.stiffness);
            let rate = arm.relaxation_time_s.recip();
            for (b, &p) in column.iter_mut().zip(arm.projection.iter().step_by(2)) {
                *b = root * p;
            }
            let strain = projection(&column, &state[..n]);
            let denominator = s + C64::from_re(rate);
            let memory = (C64::from_re(rate) / denominator) * strain;
            let mismatch = (s / denominator) * strain;
            result.material_memory.push(memory);
            let stored = 0.25 * mismatch.norm_sq();
            let heat = 0.5 * rate * mismatch.norm_sq();
            result.mean_stored_energy_j += stored;
            result.material_dissipation_w += heat;
            result.material_branch_mean_storage_j.push(stored);
            result.material_branch_dissipation_w.push(heat);
        }
        if let Some(air) = &self.radiation {
            for (j, pole) in air.model.poles.iter().enumerate() {
                poll(Some(gate))?;
                let q = state[n + j];
                let p = s * q;
                result.radiation_state.extend([q, p]);
                result.mean_stored_energy_j += 0.25 * (q.scale(pole.omega).norm_sq() + p.norm_sq());
                result.radiation_dissipation_w += pole.zeta * pole.omega * p.norm_sq();
            }
        }
        result.power_residual_w = result.total_dissipation_w() - result.external_power_w;
        if result
            .velocity_m_sqrt_kg_per_s
            .iter()
            .chain(&result.material_memory)
            .chain(&result.radiation_state)
            .chain(&result.connection_forces_n)
            .chain(&result.connection_spring_forces_n)
            .chain(&result.connection_damping_forces_n)
            .chain(&result.generalized_connection_forces_n_per_sqrt_kg)
            .any(|&z| !complex_finite(z))
            || [
                result.mean_stored_energy_j,
                result.external_power_w,
                result.total_dissipation_w(),
                result.power_residual_w,
            ]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(invalid("harmonic physical power or storage overflow"));
        }
        poll(Some(gate))?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
