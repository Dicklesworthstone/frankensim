//! Prescribed pin motion on the same condensed shell trajectory.
//!
//! Free rows receive f_f - M_fc a_c - C_fc v_c - K_fc q_c. Prescribed rows
//! retain their actual inertia, reactions, work and storage in the full model.
//! Endpoint acceleration drives the same average-acceleration kinematics as
//! the native free solve, avoiding inconsistent independent q/v/a samples.

use super::{
    DynamicState, MovingLoadError, MovingShellModel, MovingShellMotion, MovingShellStep,
    ShellPointWrench, ShellTimeConfig, finite, invalid,
};

/// Prescribed Cartesian translations of the original three pinned nodes.
///
/// Rows follow `ShellPlate::support.node_indices` in the supplied order; each
/// row has Cartesian x/y/z components in the model's fixed SI frame. Rotations
/// remain native free DOFs. The initial q/v/a may be nonzero; later calls supply
/// only endpoint acceleration, and the trajectory owns compatible q and v.
/// Support membership and the linear reference geometry stay fixed.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ShellSupportMotion {
    /// Prescribed displacement from each original node, m.
    pub displacement_m: [[f64; 3]; 3],
    /// Prescribed velocity in the model frame, m/s.
    pub velocity_m_s: [[f64; 3]; 3],
    /// Prescribed acceleration in the model frame, m/s².
    pub acceleration_m_s2: [[f64; 3]; 3],
}

impl ShellSupportMotion {
    pub(super) fn advance(
        &self,
        dt: f64,
        acceleration_m_s2: [[f64; 3]; 3],
    ) -> Result<Self, MovingLoadError> {
        for row in &acceleration_m_s2 {
            finite(row)?;
        }
        let mut next = *self;
        next.acceleration_m_s2 = acceleration_m_s2;
        for i in 0..3 {
            for j in 0..3 {
                let increment =
                    dt * f64::midpoint(self.acceleration_m_s2[i][j], acceleration_m_s2[i][j]);
                next.velocity_m_s[i][j] = self.velocity_m_s[i][j] + increment;
                next.displacement_m[i][j] =
                    self.displacement_m[i][j] + dt * (self.velocity_m_s[i][j] + 0.5 * increment);
            }
        }
        next.check()?;
        Ok(next)
    }

    fn check(&self) -> Result<(), MovingLoadError> {
        for field in [
            &self.displacement_m,
            &self.velocity_m_s,
            &self.acceleration_m_s2,
        ] {
            for row in field {
                finite(row)?;
            }
        }
        Ok(())
    }
}

impl MovingShellModel {
    /// Start native shell motion with prescribed motion at its three pins.
    ///
    /// Initial free acceleration is solved from the ORIGINAL full mechanical
    /// equation, including prescribed inertial, elastic and viscous terms.
    /// Full mechanical energy includes the moving pin masses and cross terms.
    /// Free displacement/velocity use the existing assembly's reduced order;
    /// `supports` uses the original support-node order. No operator is rebuilt.
    ///
    /// This is prescribed excitation by an external actuator, not a finite
    /// support-body model. Reaction and signed work expose that actuator's
    /// exchange. Call [`MovingShellMotion::advance_driven`] for every next step.
    pub fn start_driven<'a>(
        &'a self,
        displacement: &[f64],
        velocity: &[f64],
        supports: ShellSupportMotion,
        initial_loads: &[ShellPointWrench],
        config: ShellTimeConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MovingShellMotion<'a>, MovingLoadError> {
        self.start_with_supports(
            displacement,
            velocity,
            initial_loads,
            config,
            Some(supports),
            cancelled,
        )
    }

    pub(super) fn check_supports(
        &self,
        motion: &ShellSupportMotion,
    ) -> Result<(), MovingLoadError> {
        if self.plate.support.is_none() {
            return Err(invalid(
                "driven motion requires the original three pinned supports",
            ));
        }
        motion.check()
    }

    pub(super) fn support_value(&self, dof: usize, field: &[[f64; 3]; 3]) -> f64 {
        let supports = self.plate.support.expect("validated support membership");
        let point = supports
            .node_indices
            .iter()
            .position(|&node| 6 * node <= dof && dof < 6 * node + 3)
            .expect("original restrained translation");
        field[point][dof % 6]
    }

    pub(super) fn support_row(&self, matrix: &[f64], row: usize, field: &[[f64; 3]; 3]) -> f64 {
        let supports = self.plate.support.expect("validated support membership");
        let width = self.free_index.len();
        let mut sum = 0.0;
        for (&node, value) in supports.node_indices.iter().zip(field) {
            for axis in 0..3 {
                sum += matrix[row * width + 6 * node + axis] * value[axis];
            }
        }
        sum
    }

    pub(super) fn effective_load(
        &self,
        mut load: Vec<f64>,
        supports: Option<&ShellSupportMotion>,
    ) -> Result<Vec<f64>, MovingLoadError> {
        if let Some(motion) = supports {
            for (value, &row) in load.iter_mut().zip(&self.assembly.free_dofs) {
                *value -= self.support_row(
                    self.assembly.full_mass.values(),
                    row,
                    &motion.acceleration_m_s2,
                ) + self.support_row(
                    self.assembly.full_stiffness.values(),
                    row,
                    &motion.displacement_m,
                ) + self.assembly.full_damping.as_ref().map_or(0.0, |c| {
                    self.support_row(c.values(), row, &motion.velocity_m_s)
                });
            }
            finite(&load)?;
        }
        Ok(load)
    }

    pub(super) fn expand_with_supports(
        &self,
        free: &[f64],
        field: Option<&[[f64; 3]; 3]>,
    ) -> Vec<f64> {
        let mut full = self.expand(free);
        if let Some(field) = field {
            let supports = self.plate.support.expect("validated support membership");
            for (&node, value) in supports.node_indices.iter().zip(field) {
                full[6 * node..6 * node + 3].copy_from_slice(value);
            }
        }
        full
    }

    pub(super) fn refresh_full_energy(
        &self,
        state: &mut DynamicState,
        motion: &ShellSupportMotion,
    ) -> Result<(), MovingLoadError> {
        let velocity = self.expand_with_supports(&state.velocity, Some(&motion.velocity_m_s));
        let displacement =
            self.expand_with_supports(&state.displacement, Some(&motion.displacement_m));
        state.kinetic_energy_j = self.assembly.full_mass.quadratic_energy(&velocity);
        state.strain_energy_j = self.assembly.full_stiffness.quadratic_energy(&displacement);
        state.total_energy_j = state.kinetic_energy_j + state.strain_energy_j;
        finite(&[
            state.kinetic_energy_j,
            state.strain_energy_j,
            state.total_energy_j,
        ])
    }
}

impl MovingShellMotion<'_> {
    /// Accepted prescribed pin motion, or `None` on the stationary-support API.
    pub fn support_motion(&self) -> Option<&ShellSupportMotion> {
        self.supports.as_ref()
    }

    /// Advance the original free dynamics with explicit next pin acceleration.
    ///
    /// The pinned translations use v_next=v+dt*(a+a_next)/2 and
    /// q_next=q+dt*(v+v_next)/2. The complete-increment form retains small
    /// displacements beside the reference datum. The original free solver
    /// receives the equivalent endpoint load including full M/C/K coupling;
    /// its prepared factorization is reused unchanged. Reaction work and
    /// prescribed momentum/storage are included in the same returned account.
    /// Resolve time-varying acceleration through timestep refinement.
    ///
    /// The original timestamp, horizon and cancellation rules apply. A call
    /// without initialized driven supports or with any nonfinite acceleration
    /// refuses without changing the accepted free or prescribed motion.
    pub fn advance_driven(
        &mut self,
        time_s: f64,
        loads: &[ShellPointWrench],
        acceleration_m_s2: [[f64; 3]; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MovingShellStep, MovingLoadError> {
        self.advance_with_supports(time_s, loads, Some(acceleration_m_s2), cancelled)
    }
}
