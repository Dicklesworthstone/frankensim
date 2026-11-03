//! Zero-duration momentum transfer on the original sparse transient operators.
//!
//! The event changes velocity through M_ff delta-v = J_f. Prescribed support
//! velocities remain fixed; their full consistent-mass impulse and work are
//! retained. Finite stiffness, damping and thermal forces carry no impulse,
//! but the new velocity changes the subsequent acceleration and reactions.

use super::super::*;

/// Physical residual tolerances for an instantaneous impulse [G0/G1].
/// Iteration and storage budgets remain those of the admitted transient run.
#[derive(Debug, Clone, Copy)]
pub struct NodalImpulseControls3 {
    pub absolute_momentum_tolerance_kg_m_s: f64,
    pub absolute_work_tolerance_j: f64,
    /// Applies to the impulse norm and to the independent impulse-work check.
    pub relative_tolerance: f64,
}

impl Default for NodalImpulseControls3 {
    fn default() -> Self {
        Self {
            absolute_momentum_tolerance_kg_m_s: 1e-12,
            absolute_work_tolerance_j: 1e-12,
            relative_tolerance: 1e-9,
        }
    }
}

/// Physical impulse and work at unchanged accepted time and displacement.
#[derive(Debug, Clone, PartialEq)]
pub struct NodalImpulse3 {
    pub time_s: f64,
    pub step: usize,
    pub impulse_solve_iterations: usize,
    pub momentum_residual_kg_m_s: f64,
    /// Full original-node momentum increment M*(v_after-v_before).
    /// Includes restrained rows and their off-diagonal inertia coupling.
    pub nodal_momentum_change_kg_m_s: Vec<[f64; 3]>,
    /// Impulse exerted BY the supports ON the solid, zero on free components.
    /// On restrained rows this is M*delta-v - J_applied.
    pub support_impulse_kg_m_s: Vec<[f64; 3]>,
    /// Physical full-field kinetic energies, ordered [before, after].
    pub kinetic_energy_j: [f64; 2],
    /// Stable average(v) dot M*delta-v, avoiding subtraction of large energies.
    pub kinetic_change_j: f64,
    pub applied_impulse_work_j: f64,
    /// average(v_support) dot J_support; moving supports can exchange energy.
    pub support_impulse_work_j: f64,
    pub work_residual_j: f64,
    /// Subsequent finite-force equilibrium, including changed damping forces.
    pub acceleration_solve_iterations: usize,
    pub force_residual_n: f64,
}

impl TetTransientRun3<'_, '_, '_, '_> {
    /// Apply explicit sorted, unique full-DOF impulses [N s = kg m/s].
    ///
    /// The same sparse consistent M used by Newmark supplies the velocity jump;
    /// no mass lumping or finite-duration force approximation enters. Applied
    /// impulses on fixed DOFs are retained and balanced by support impulses.
    /// Prescribed supports keep their displacement, velocity and acceleration.
    /// Their impulse work is included when their held velocity is nonzero.
    ///
    /// The event retains time, step count, displacement, elastic/thermal state,
    /// source materials and finite external loads. It then re-equilibrates
    /// acceleration and support reactions under those same finite loads and the
    /// new damping force. No strain, damping heat or thermal energy is generated
    /// during the zero-duration event. Impulses are external inputs; this method
    /// does not select a collision law, contact location or event time.
    ///
    /// Publication is transactional. Momentum is checked using the ACTUAL
    /// representable velocity change, and impulse work includes full consistent
    /// support coupling. Caller-owned contact partners receive the opposite
    /// applied impulse; constraints receive the opposite reported support one.
    ///
    /// # Errors
    /// Invalid/duplicate/out-of-range impulse, invalid controls, unrepresentable
    /// motion/work, native solve failure, momentum/work refusal or cancellation.
    pub fn apply_nodal_impulse(
        &mut self,
        impulse_kg_m_s: &[(usize, f64)],
        controls: NodalImpulseControls3,
        cx: &Cx<'_>,
    ) -> Result<NodalImpulse3, TransientError3> {
        check(cx)?;
        if !controls.absolute_momentum_tolerance_kg_m_s.is_finite()
            || controls.absolute_momentum_tolerance_kg_m_s <= 0.0
            || !controls.absolute_work_tolerance_j.is_finite()
            || controls.absolute_work_tolerance_j <= 0.0
            || !controls.relative_tolerance.is_finite()
            || controls.relative_tolerance <= 0.0
            || controls.relative_tolerance >= 1.0
        {
            return Err(TransientError3::Invalid(
                "positive finite impulse/work tolerances and relative tolerance below one required",
            ));
        }
        let full = self.geometry.nodes.len() * 3;
        if impulse_kg_m_s.len() > full {
            return Err(TransientError3::Budget);
        }
        let mut applied = vec![0.0; full];
        let mut previous = None;
        for &(dof, impulse) in impulse_kg_m_s {
            check(cx)?;
            if dof >= full || previous.is_some_and(|old| dof <= old) || !impulse.is_finite() {
                return Err(TransientError3::Invalid(
                    "impulse DOFs must be sorted, unique, in range and finite",
                ));
            }
            applied[dof] = impulse;
            previous = Some(dof);
        }
        let free_impulse: Vec<_> = self.free_dofs().iter().map(|&d| applied[d]).collect();
        let impulse_scale = number(length(&free_impulse))?;
        let (velocity_increment, impulse_solve_iterations) = if impulse_scale == 0.0 {
            (vec![0.0; self.free_dofs().len()], 0)
        } else {
            // A force-unit stopping floor must not erase a small physical
            // impulse. Normalize the equation, then check its SI residual below.
            let normalized: Vec<_> = free_impulse
                .iter()
                .map(|value| value / impulse_scale)
                .collect();
            let mut solve_controls = self.config;
            solve_controls.absolute_force_tolerance_n = f64::MIN_POSITIVE;
            solve_controls.solve.relative_residual_tolerance = solve_controls
                .solve
                .relative_residual_tolerance
                .min(controls.relative_tolerance);
            let (mut increment, iterations, _) =
                match solve(&self.operator(1.0, 0.0, cx), &normalized, solve_controls) {
                    Ok(result) => result,
                    Err(TransientError3::Solve {
                        residual_n,
                        allowance_n,
                        ..
                    }) => {
                        return Err(TransientError3::Momentum {
                            residual_kg_m_s: number(impulse_scale * residual_n)?,
                            allowance_kg_m_s: number(impulse_scale * allowance_n)?,
                        });
                    }
                    Err(error) => return Err(error),
                };
            for value in &mut increment {
                check(cx)?;
                *value = number(*value * impulse_scale)?;
            }
            (increment, iterations)
        };
        let before = self.nodal_motion(cx)?;
        let mut candidate = self.clone();
        let mut full_increment = std::array::from_fn(|_| vec![0.0; full]);
        let mut average_velocity: Vec<_> = before[1].iter().flatten().copied().collect();
        let mut changed = false;
        for (i, (&dof, &increment)) in self.free_dofs().iter().zip(&velocity_increment).enumerate()
        {
            if i % 256 == 0 {
                check(cx)?;
            }
            // Preserve exact no-op bits, including signed zero.
            let next = if increment == 0.0 {
                self.state.velocity[i]
            } else {
                number(self.state.velocity[i] + increment)?
            };
            candidate.state.velocity[i] = next;
            let actual = number(next - self.state.velocity[i])?;
            full_increment[1][dof] = actual;
            average_velocity[dof] = number(0.5 * next + 0.5 * self.state.velocity[i])?;
            changed |= actual != 0.0;
        }
        // The original full element action retains M_cf, which a reduced
        // solve alone cannot supply for physical support impulse accounting.
        let momentum = self.problem.apply_full_kinematics(&full_increment, cx)?.mv;
        let mut residual = 0.0_f64;
        for &dof in self.free_dofs() {
            check(cx)?;
            residual = number(residual.hypot(number(momentum[dof] - applied[dof])?))?;
        }
        let allowance = number(
            controls.absolute_momentum_tolerance_kg_m_s
                + controls.relative_tolerance * impulse_scale,
        )?;
        if residual > allowance {
            return Err(TransientError3::Momentum {
                residual_kg_m_s: residual,
                allowance_kg_m_s: allowance,
            });
        }
        let mut support = vec![0.0; full];
        for &dof in self.problem.fixed_dofs {
            check(cx)?;
            support[dof] = number(momentum[dof] - applied[dof])?;
        }
        let kinetic_change = number(inner(&average_velocity, &momentum))?;
        let applied_work = number(inner(&average_velocity, &applied))?;
        let support_work = number(inner(&average_velocity, &support))?;
        let work_residual = number(number(kinetic_change - applied_work)? - support_work)?;
        let work_allowance = number(
            controls.absolute_work_tolerance_j
                + controls.relative_tolerance
                    * kinetic_change
                        .abs()
                        .max(applied_work.abs())
                        .max(support_work.abs()),
        )?;
        if work_residual.abs() > work_allowance {
            return Err(TransientError3::Energy {
                defect_j: work_residual,
                allowance_j: work_allowance,
            });
        }
        let (acceleration_solve_iterations, force_residual_n) = if changed {
            let ku = candidate.apply(0.0, 1.0, &candidate.state.displacement, cx)?;
            let cv = candidate.apply(
                candidate.config.mass_damping_per_s,
                candidate.config.stiffness_damping_s,
                &candidate.state.velocity,
                cx,
            )?;
            let mut rhs = Vec::with_capacity(candidate.free_dofs().len());
            for i in 0..candidate.free_dofs().len() {
                if i % 256 == 0 {
                    check(cx)?;
                }
                rhs.push(number(
                    candidate.force[i] + candidate.thermal[i] - ku[i] - cv[i],
                )?);
            }
            if let Some(supports) = &candidate.prescribed {
                candidate.subtract_support_force(supports, &mut rhs, cx)?;
            }
            let (acceleration, iterations, force_residual) =
                solve(&candidate.operator(1.0, 0.0, cx), &rhs, candidate.config)?;
            candidate.state.acceleration = acceleration;
            candidate.state.kinetic_energy_j = candidate
                .energy(&candidate.state.displacement, &candidate.state.velocity, cx)?
                .0;
            candidate.state.total_energy_j =
                number(candidate.state.kinetic_energy_j + candidate.state.strain_energy_j)?;
            // The existing owner restores full supported kinetic energy and
            // finite-force reactions without changing any support kinematics.
            candidate.initialize_support_energy(cx)?;
            (iterations, force_residual)
        } else {
            (0, 0.0)
        };
        let report = NodalImpulse3 {
            time_s: self.state.time_s,
            step: self.state.step,
            impulse_solve_iterations,
            momentum_residual_kg_m_s: residual,
            nodal_momentum_change_kg_m_s: momentum.as_chunks::<3>().0.to_vec(),
            support_impulse_kg_m_s: support.as_chunks::<3>().0.to_vec(),
            kinetic_energy_j: [
                self.state.kinetic_energy_j,
                candidate.state.kinetic_energy_j,
            ],
            kinetic_change_j: kinetic_change,
            applied_impulse_work_j: applied_work,
            support_impulse_work_j: support_work,
            work_residual_j: work_residual,
            acceleration_solve_iterations,
            force_residual_n,
        };
        check(cx)?;
        *self = candidate;
        Ok(report)
    }
}
