//! Material changes at one accepted instant without resetting strain history.

use super::*;

/// Explicit allowance for externally prescribed coefficient work at fixed
/// displacement, momentum and viscous strain. This is not a heat allowance.
#[derive(Clone, Copy, Debug)]
pub struct PreparedTensorMaterialChangeBudget {
    pub maximum_absolute_parameter_work_j: f64,
}

/// Complete proposed continuation. The original history remains untouched on
/// success and refusal; publishing the returned history is the caller's act.
#[derive(Clone)]
pub struct PreparedTensorMaterialChange<'a> {
    pub history: PreparedTensorHistory<'a>,
    /// Signed equilibrium plus Maxwell storage change [J].
    pub parameter_work_j: f64,
    pub equilibrium_change_j: f64,
    pub memory_change_j: f64,
}

impl PreparedTensorHistory<'_> {
    /// Continue with changed elastic tensors, Kelvin viscosity and relaxation
    /// times while retaining every original branch's complete spring strain.
    ///
    /// Points and branches must correspond in their ORIGINAL physical order.
    /// This numerical owner checks reference volume, strain columns, mass,
    /// point heat ownership and branch shape. The source-binding caller must
    /// additionally establish the identity of those internal variables and the
    /// applicability of the prescribed material change. Relabeling a foreign
    /// branch as an existing one is not a valid continuation.
    ///
    /// An old time group may split by copying its exact displacement history.
    /// Several old groups may acquire the same new time only when their entire
    /// histories are exactly equal. Otherwise the new body's shared realization
    /// cannot represent the state and refuses; no strain averaging or numerical
    /// tolerance erases memory. Complete retained strain also permits changed
    /// factor rank without reconstruction, inversion or spectral truncation.
    ///
    /// Geometry, basis, reference mass and positive step duration stay fixed.
    /// Both bodies must have zero eigenstrain. Held/affine preparation and its
    /// conditioning scale may change because they carry no physical history.
    /// No time, relaxation or heat occurs here. The signed storage change has
    /// a distinct compensated parameter-work account; force work, heat, initial
    /// energy, accepted steps and momentum all survive unchanged. This is an
    /// external parameter operation, not a solved thermal/entropy feedback.
    pub fn rebind_materials<'b>(
        &self,
        next: &'b PreparedTensorBody,
        budget: PreparedTensorMaterialChangeBudget,
        gate: &CancelGate,
    ) -> Result<PreparedTensorMaterialChange<'b>> {
        poll(gate)?;
        let old = self.body;
        let n = old.coordinates();
        if !budget.maximum_absolute_parameter_work_j.is_finite()
            || budget.maximum_absolute_parameter_work_j < 0.0
            || next.coordinates() != n
            || next.mass_matrix_kg() != old.mass_matrix_kg()
            || next.dt_s.to_bits() != old.dt_s.to_bits()
            || next.points.len() != old.points.len()
            || next.heat_groups != old.heat_groups
            || self.steps > next.controls.maximum_steps
        {
            return Err(invalid(
                "material continuation requires fixed mass, point ownership, clock and bounded work/steps",
            ));
        }

        let mut parents: Vec<Option<usize>> = Vec::new();
        parents
            .try_reserve_exact(next.times.len())
            .map_err(|_| invalid("material history mapping allocation refused"))?;
        parents.resize(next.times.len(), None);
        let (mut memory_change_j, mut memory_roundoff) = (0.0, 0.0);
        for (index, (before, after)) in old.points.iter().zip(&next.points).enumerate() {
            poll(gate)?;
            if before.volume_m3().to_bits() != after.volume_m3().to_bits()
                || before.strain_columns() != after.strain_columns()
                || before.compatible_strain_columns() != after.compatible_strain_columns()
                || before.branches().len() != after.branches().len()
                || before.free_strain_mandel().iter().any(|v| *v != 0.0)
                || after.free_strain_mandel().iter().any(|v| *v != 0.0)
            {
                return Err(invalid(
                    "material continuation requires original strain maps, volumes, ordered branches and zero eigenstrain",
                ));
            }
            for (branch, (previous, replacement)) in
                before.branches().iter().zip(after.branches()).enumerate()
            {
                poll(gate)?;
                let source = old.point_times[index][branch];
                let target = next.point_times[index][branch];
                let history = &self.state[(2 + source) * n..(3 + source) * n];
                if let Some(first) = parents[target] {
                    if self.state[(2 + first) * n..(3 + first) * n] != *history {
                        return Err(PreparedTensorError::MemoryCollision {
                            relaxation_time_s: next.times[target],
                        });
                    }
                } else {
                    parents[target] = Some(source);
                }
                // .5*(z_new^2-z_old^2) = mean(z)*delta(z). Form each
                // difference from original factor rows before contraction;
                // unchanged factors give exact zero despite large storage.
                for (old_row, new_row) in previous.rows().iter().zip(replacement.rows()) {
                    let (mut mean, mut increment) = (0.0, 0.0);
                    for j in 0..n {
                        mean = (0.5 * old_row[j] + 0.5 * new_row[j]).mul_add(history[j], mean);
                        increment = (new_row[j] - old_row[j]).mul_add(history[j], increment);
                    }
                    add(&mut memory_change_j, &mut memory_roundoff, mean * increment);
                }
            }
        }

        let mut state = zeros(next.states)?;
        state[..2 * n].copy_from_slice(&self.state[..2 * n]);
        for (target, source) in parents.into_iter().enumerate() {
            poll(gate)?;
            let source = source.ok_or_else(|| invalid("unmapped replacement material history"))?;
            state[(2 + target) * n..(3 + target) * n]
                .copy_from_slice(&self.state[(2 + source) * n..(3 + source) * n]);
        }
        next.admit_displacement(&state[..n], gate)?;
        let energy = next.energy(&state, &self.velocity, gate)?;
        let equilibrium_change_j =
            coefficient_change(&old.equilibrium, &next.equilibrium, &self.state[..n], gate)?;
        // Kinetic storage is unchanged. Never subtract two total energies:
        // large kinetic energy must not erase the actual coefficient work.
        let parameter_work_j = equilibrium_change_j + memory_change_j;
        if !finite(&[equilibrium_change_j, memory_change_j, parameter_work_j]) {
            return Err(invalid("material coefficient work overflow"));
        }
        if parameter_work_j.abs() > budget.maximum_absolute_parameter_work_j {
            return Err(PreparedTensorError::ParameterWork {
                work_j: parameter_work_j,
                allowance_j: budget.maximum_absolute_parameter_work_j,
            });
        }
        let (mut accumulated_work, mut work_roundoff) =
            (self.parameter_work_j, self.parameter_work_correction_j);
        add(&mut accumulated_work, &mut work_roundoff, parameter_work_j);
        if !finite(&[accumulated_work, work_roundoff]) {
            return Err(invalid("material coefficient work ledger overflow"));
        }
        let defect_j =
            self.work_j + accumulated_work - self.heat_j - (energy.total() - self.initial_energy_j);
        if !defect_j.is_finite() || defect_j.abs() > next.controls.maximum_work_error_j {
            return Err(PreparedTensorError::Energy {
                defect_j,
                allowance_j: next.controls.maximum_work_error_j,
            });
        }
        let scratch = Scratch::new(next.states, n, next.affine.is_some())?;
        let mut velocity = zeros(n)?;
        velocity.copy_from_slice(&self.velocity);
        poll(gate)?;
        Ok(PreparedTensorMaterialChange {
            history: PreparedTensorHistory {
                body: next,
                state,
                velocity,
                scratch,
                steps: self.steps,
                energy,
                initial_energy_j: self.initial_energy_j,
                work_j: self.work_j,
                work_correction_j: self.work_correction_j,
                parameter_work_j: accumulated_work,
                parameter_work_correction_j: work_roundoff,
                heat_j: self.heat_j,
                heat_correction_j: self.heat_correction_j,
            },
            parameter_work_j,
            equilibrium_change_j,
            memory_change_j,
        })
    }
}

/// Fixed-strain .5*u^T*(K_new-K_old)*u, excluding unchanged kinetic energy.
fn coefficient_change(
    before: &[f64],
    after: &[f64],
    displacement: &[f64],
    gate: &CancelGate,
) -> Result<f64> {
    let n = displacement.len();
    let (mut change, mut roundoff) = (0.0, 0.0);
    for i in 0..n {
        poll(gate)?;
        let row = (0..n).fold(0.0, |sum, j| {
            (after[i * n + j] - before[i * n + j]).mul_add(displacement[j], sum)
        });
        add(&mut change, &mut roundoff, 0.5 * displacement[i] * row);
    }
    if !change.is_finite() {
        return Err(invalid("equilibrium coefficient work overflow"));
    }
    Ok(change)
}
