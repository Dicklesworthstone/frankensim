//! Native backward-Euler inertia and the shared physical stress quadratures.
//!
//! Dotting M(u1-u0)/dt + A(u1) + grad(p1) = f by u1 gives
//! W_wall + W_body + W_pressure = Q_physical + (K1-K0) + N_BE,
//! with K = sum(M*u^2)/2 and N_BE = sum(M*(u1-u0)^2)/2.
//! M uses the SAME staggered half-cell volumes as momentum(), including graded
//! grids. N_BE is numerical loss; it is reported and NEVER deposited into h.
//! W_body includes material buoyancy and external forces. Buoyancy's signed
//! subset is attributed separately, without a gravitational-potential store.

use super::*;

/// One accepted-duration native unsteady-Stokes mechanical account.
///
/// All powers use the projected endpoint faces and the source-start material
/// coefficients. The original thermal owner deposits only power.cell_heat_w.
/// K is the native face mass norm, distinct from FlowStepRecord's cell-centred
/// diagnostic. This is a first-order discrete balance, not an error certificate.
#[derive(Debug, Clone, PartialEq)]
pub struct UnsteadyMechanicalStep {
    pub start_s: f64,
    pub duration_s: f64,
    pub end_s: f64,
    /// Body power in this receipt is the TOTAL of buoyancy and external force.
    pub power: ViscousPower,
    /// Signed buoyancy work rate per original voxel [W], zero on fixed solids.
    /// This is a subset of power.cell_body_power_w, NEVER additional heat.
    /// None when no material buoyancy law is attached.
    pub cell_buoyancy_power_w: Option<Vec<f64>>,
    /// Signed subset of power.supplied_body_power_w on the same native faces.
    /// External actuator power is total body power minus this value. There is
    /// no resolved gravitational-potential reservoir debited by this account.
    pub supplied_buoyancy_power_w: f64,
    pub initial_kinetic_energy_j: f64,
    pub final_kinetic_energy_j: f64,
    /// Nonnegative backward-Euler loss, outside physical heat.
    pub backward_euler_loss_j: f64,
    /// Supply + pressure - heat - kinetic change - numerical loss [J].
    pub residual_j: f64,
}

/// One accepted thermal advance's mechanical storage, supply and heat account.
/// Physical heat is already in the thermal EnergyAudit.external_heat_j.
/// Neither actuator work nor numerical loss is an additional thermal source.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct UnsteadyHeatingAudit {
    pub duration_s: f64,
    pub substeps: usize,
    pub initial_kinetic_energy_j: f64,
    pub final_kinetic_energy_j: f64,
    pub wall_work_j: [f64; 6],
    pub supplied_wall_work_j: f64,
    /// TOTAL signed work of buoyancy and external body forces.
    pub supplied_body_work_j: f64,
    /// Signed buoyancy subset of supplied_body_work_j, not additional supply
    /// or heat. Subtract it from total body work to recover actuator work.
    pub supplied_buoyancy_work_j: f64,
    pub generated_heat_j: f64,
    pub molecular_heat_j: f64,
    pub darcy_heat_j: f64,
    pub forchheimer_heat_j: f64,
    /// Inelastic kinetic capture by an explicitly stationary phase support.
    /// Separate from molecular/Darcy viscosity and numerical time loss.
    pub support_capture_heat_j: f64,
    pub backward_euler_loss_j: f64,
    /// Bounded pressure-constraint defect; never heat.
    pub pressure_constraint_work_j: f64,
    pub residual_j: f64,
}

impl UnsteadyHeatingAudit {
    pub(in crate::conjugate) fn from_kinetic_energy(kinetic_energy_j: f64) -> Self {
        Self {
            initial_kinetic_energy_j: kinetic_energy_j,
            final_kinetic_energy_j: kinetic_energy_j,
            ..Self::default()
        }
    }

    pub(in crate::conjugate) fn accept(
        &mut self,
        step: &UnsteadyMechanicalStep,
        controls: ViscousHeatControls,
    ) -> Result<(), ChtError> {
        if self.final_kinetic_energy_j.to_bits() != step.initial_kinetic_energy_j.to_bits() {
            return Err(invalid_power(
                "transient heat must retain the same native kinetic history",
            ));
        }
        let mut next = *self;
        let add = |sum: &mut f64, value: f64| -> Result<(), ChtError> {
            *sum = finite_value("unsteady_heat.integral", *sum + value)?;
            Ok(())
        };
        add(&mut next.duration_s, step.duration_s)?;
        next.substeps = next
            .substeps
            .checked_add(1)
            .ok_or_else(|| invalid_power("transient heat substep count overflow"))?;
        next.final_kinetic_energy_j = step.final_kinetic_energy_j;
        for (sum, power) in next.wall_work_j.iter_mut().zip(step.power.wall_power_w) {
            add(sum, step.duration_s * power)?;
        }
        for (sum, power) in [
            (
                &mut next.supplied_wall_work_j,
                step.power.supplied_wall_power_w,
            ),
            (
                &mut next.supplied_body_work_j,
                step.power.supplied_body_power_w,
            ),
            (
                &mut next.supplied_buoyancy_work_j,
                step.supplied_buoyancy_power_w,
            ),
            (&mut next.generated_heat_j, step.power.heat_w),
            (&mut next.molecular_heat_j, step.power.molecular_heat_w),
            (&mut next.darcy_heat_j, step.power.darcy_heat_w),
            (&mut next.forchheimer_heat_j, step.power.forchheimer_heat_w),
            (
                &mut next.pressure_constraint_work_j,
                step.power.pressure_constraint_power_w,
            ),
        ] {
            add(sum, step.duration_s * power)?;
        }
        add(&mut next.backward_euler_loss_j, step.backward_euler_loss_j)?;
        // Recompute the total from independent extensive accounts. Never
        // adjust a heat entry to make this residual vanish.
        next.residual_j = finite_value(
            "unsteady_heat.residual",
            next.supplied_wall_work_j + next.supplied_body_work_j + next.pressure_constraint_work_j
                - next.generated_heat_j
                - (next.final_kinetic_energy_j - next.initial_kinetic_energy_j)
                - next.backward_euler_loss_j,
        )?;
        let allowance = finite_value(
            "unsteady_heat.allowance",
            controls.power_tolerance_w * next.duration_s,
        )?;
        if next.residual_j.abs() > allowance
            || next.pressure_constraint_work_j.abs() > allowance
            || (next.generated_heat_j - next.molecular_heat_j - next.darcy_heat_j
                - next.forchheimer_heat_j
                - next.support_capture_heat_j).abs() > allowance
        {
            return Err(invalid_power(
                "integrated transient mechanical and physical heat accounts do not close",
            ));
        }
        *self = next;
        Ok(())
    }

    #[cfg(feature = "enthalpy-coupling")]
    pub(in crate::conjugate) fn accept_support_capture(
        &mut self,
        before_j: f64,
        after_j: f64,
        heat_j: f64,
        controls: ViscousHeatControls,
    ) -> Result<(), ChtError> {
        if before_j.to_bits() != self.final_kinetic_energy_j.to_bits()
            || !after_j.is_finite() || after_j < 0.0 || !heat_j.is_finite() || heat_j < 0.0
        {
            return Err(invalid_power("phase capture must retain the actual native kinetic history"));
        }
        let mut next = *self;
        next.final_kinetic_energy_j = after_j;
        next.generated_heat_j = finite_value("capture.heat", next.generated_heat_j + heat_j)?;
        next.support_capture_heat_j = finite_value("capture.heat", next.support_capture_heat_j + heat_j)?;
        next.residual_j = finite_value("capture.residual",
            next.supplied_wall_work_j + next.supplied_body_work_j + next.pressure_constraint_work_j
                - next.generated_heat_j
                - (next.final_kinetic_energy_j - next.initial_kinetic_energy_j)
                - next.backward_euler_loss_j)?;
        let allowance = finite_value("capture.allowance", controls.power_tolerance_w * next.duration_s)?;
        if (before_j - after_j - heat_j).abs() > allowance
            || next.residual_j.abs() > allowance
        {
            return Err(invalid_power("supported capture does not close native kinetic and heat storage"));
        }
        *self = next;
        Ok(())
    }
}

impl Solver<'_> {
    pub(in crate::conjugate) fn enable_unsteady_stokes_heat(
        &mut self,
        controls: ViscousHeatControls,
    ) -> Result<(), ChtError> {
        admit_heating(self.domain, self.config, controls)?;
        // Only momentum convection is omitted. Native rho*V/dt, pressure,
        // full symmetric stress and stationary Darcy coefficients remain.
        self.creeping = true;
        Ok(())
    }

    pub(in crate::conjugate) fn unsteady_mechanical_step(
        &self,
        start_s: f64,
        end_s: f64,
        buoyancy_force_n_m3: Option<&[[f64; 3]]>,
        controls: ViscousHeatControls,
        gate: &CancelGate,
    ) -> Result<UnsteadyMechanicalStep, ChtError> {
        let time = self
            .time
            .as_ref()
            .ok_or_else(|| invalid_power("transient heat requires native time levels"))?;
        if !self.creeping || time.scheme != TimeScheme::BackwardEuler || time.older.is_some() {
            return Err(invalid_power(
                "transient heat requires one native backward-Euler Stokes step",
            ));
        }
        let mut power = self.sample_viscous_power(&self.field(), &self.pressure, controls, gate)?;
        let (cell_buoyancy_power_w, supplied_buoyancy_power_w) = match buoyancy_force_n_m3 {
            Some(force) => {
                let (cells, power) = self.body_force_power(force, controls, gate)?;
                (Some(cells), power)
            }
            None => (None, 0.0),
        };
        power.reynolds = power.reynolds.max(self.source_reynolds(&time.previous)?);
        if power.reynolds > controls.maximum_reynolds {
            return Err(invalid_power(
                "source viscosity and retained starting velocity exceed the unsteady Stokes bound",
            ));
        }
        let mut initial = 0.0;
        let mut final_energy = 0.0;
        let mut numerical_loss = 0.0;
        for a in 0..3 {
            for &index in &self.comps[a].unknowns {
                poll(gate)?;
                // Exactly rho*face_area*length from the native inertia row.
                let mass = self.interior_face_mass_kg(a, index);
                let old = time.previous[a][index];
                let new = self.vel[a][index];
                initial = finite_value(
                    "unsteady_heat.initial_kinetic",
                    initial + 0.5 * mass * old * old,
                )?;
                final_energy = finite_value(
                    "unsteady_heat.final_kinetic",
                    final_energy + 0.5 * mass * new * new,
                )?;
                let difference = new - old;
                numerical_loss = finite_value(
                    "unsteady_heat.numerical_loss",
                    numerical_loss + 0.5 * mass * difference * difference,
                )?;
            }
        }
        let residual_j = finite_value(
            "unsteady_heat.step_residual",
            time.dt * power.residual_w - (final_energy - initial) - numerical_loss,
        )?;
        let allowance = finite_value(
            "unsteady_heat.step_allowance",
            time.dt * controls.power_tolerance_w,
        )?;
        if residual_j.abs() > allowance {
            return Err(invalid_power(
                "transient wall/body work, native kinetic storage and physical/numerical losses do not close",
            ));
        }
        poll(gate)?;
        Ok(UnsteadyMechanicalStep {
            start_s,
            duration_s: time.dt,
            end_s,
            power,
            cell_buoyancy_power_w,
            supplied_buoyancy_power_w,
            initial_kinetic_energy_j: initial,
            final_kinetic_energy_j: final_energy,
            backward_euler_loss_j: numerical_loss,
            residual_j,
        })
    }
}
