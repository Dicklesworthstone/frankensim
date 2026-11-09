//! Native backward-Euler inertia and the shared physical stress quadratures.
//!
//! Dotting M(u1-u0)/dt + A(u1) + grad(p1) = f by u1 gives
//! W_wall + W_body + W_pressure = Q_physical + (K1-K0) + N_BE,
//! with K = sum(M*u^2)/2 and N_BE = sum(M*(u1-u0)^2)/2.
//! M uses the SAME physical staggered volumes as momentum(), including graded
//! grids and the HALF cells at normal-traction reservoir mouths. W_pressure
//! includes their signed mechanical supply and the bounded continuity defect.
//! N_BE is numerical loss; it is reported and NEVER deposited into h.
//! This is unsteady Stokes: momentum convection and its kinetic-energy flux
//! are omitted together. The existing low-Reynolds applicability gate remains.
//! W_body includes material buoyancy and external forces. Buoyancy's signed
//! subset is attributed separately, without a gravitational-potential store.

use super::*;
use fs_math::dd::Dd;

/// One accepted-duration native unsteady-Stokes mechanical account.
///
/// All powers use the projected endpoint faces and the selected original
/// material epoch. Explicit coupling selects the source-start epoch; implicit
/// material coupling independently closes the actual source endpoint epoch.
/// The original thermal owner deposits only power.cell_heat_w.
/// K uses the native face mass norm. Its increment is retained directly from
/// the accepted face speeds, independently of rounded absolute diagnostics.
/// This is a first-order discrete balance, not an error certificate.
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
    /// Actual accepted native kinetic increment, without total subtraction [J].
    pub kinetic_energy_change_j: f64,
    // Carry face-reduction tails into the accepted multi-interval account.
    pub(in crate::conjugate) kinetic_energy_change: Dd,
    /// Nonnegative backward-Euler loss, outside physical heat.
    pub backward_euler_loss_j: f64,
    pub(in crate::conjugate) backward_euler_loss: Dd,
    /// Supply + pressure - heat - kinetic change - numerical loss [J].
    pub residual_j: f64,
}

impl UnsteadyMechanicalStep {
    pub(in crate::conjugate) fn admit_balance(
        &self,
        controls: ViscousHeatControls,
    ) -> Result<(), ChtError> {
        let allowance = finite_value(
            "unsteady_heat.step_allowance",
            self.duration_s * controls.power_tolerance_w,
        )?;
        if self.residual_j.abs() > allowance {
            return Err(invalid_power(
                "transient wall/body work, native kinetic storage and physical/numerical losses do not close",
            ));
        }
        Ok(())
    }
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
    /// Sum of actual native increments, including supported phase capture [J].
    /// Rounded absolute endpoint energies need not resolve this value.
    pub kinetic_energy_change_j: f64,
    kinetic_energy_change: Dd,
    pub wall_work_j: [f64; 6],
    pub supplied_wall_work_j: f64,
    /// Signed reservoir work by original plane, in Face3::ALL order [J].
    /// Positive supplies mechanical energy; negative returns it to a reservoir.
    /// Each component uses the same current reference as native momentum.
    pub pressure_work_j: [f64; 6],
    /// Sum of pressure_work_j, separate from the continuity-constraint defect.
    /// This is mechanical supply, never an additional thermal source.
    pub supplied_pressure_work_j: f64,
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
    backward_euler_loss: Dd,
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
        next.kinetic_energy_change = next.kinetic_energy_change + step.kinetic_energy_change;
        next.kinetic_energy_change_j = finite_value(
            "unsteady_heat.kinetic_change",
            next.kinetic_energy_change.to_f64(),
        )?;
        for (sum, power) in next.wall_work_j.iter_mut().zip(step.power.wall_power_w) {
            add(sum, step.duration_s * power)?;
        }
        for (sum, power) in next
            .pressure_work_j
            .iter_mut()
            .zip(step.power.pressure_power_w)
        {
            add(sum, step.duration_s * power)?;
        }
        for (sum, power) in [
            (
                &mut next.supplied_wall_work_j,
                step.power.supplied_wall_power_w,
            ),
            (
                &mut next.supplied_pressure_work_j,
                step.power.supplied_pressure_power_w,
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
        next.backward_euler_loss = next.backward_euler_loss + step.backward_euler_loss;
        next.backward_euler_loss_j = finite_value(
            "unsteady_heat.numerical_loss", next.backward_euler_loss.to_f64(),
        )?;
        // Recompute the total from independent extensive accounts. Never
        // adjust a heat entry to make this residual vanish.
        next.residual_j = finite_value(
            "unsteady_heat.residual",
            (Dd::from_f64(next.supplied_wall_work_j)
                + Dd::from_f64(next.supplied_body_work_j)
                + Dd::from_f64(next.supplied_pressure_work_j)
                + Dd::from_f64(next.pressure_constraint_work_j)
                - Dd::from_f64(next.generated_heat_j)
                - next.kinetic_energy_change
                - next.backward_euler_loss)
            .to_f64(),
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
        kinetic_energy_change: Dd,
        heat_j: f64,
        controls: ViscousHeatControls,
    ) -> Result<(), ChtError> {
        if before_j.to_bits() != self.final_kinetic_energy_j.to_bits()
            || !after_j.is_finite() || after_j < 0.0 || !heat_j.is_finite() || heat_j < 0.0
            || !kinetic_energy_change.to_f64().is_finite()
            || kinetic_energy_change.to_f64() > 0.0
        {
            return Err(invalid_power("phase capture must retain the actual native kinetic history"));
        }
        let mut next = *self;
        next.final_kinetic_energy_j = after_j;
        next.kinetic_energy_change = next.kinetic_energy_change + kinetic_energy_change;
        next.kinetic_energy_change_j = finite_value(
            "capture.kinetic_change", next.kinetic_energy_change.to_f64(),
        )?;
        next.generated_heat_j = finite_value("capture.heat", next.generated_heat_j + heat_j)?;
        next.support_capture_heat_j = finite_value("capture.heat", next.support_capture_heat_j + heat_j)?;
        next.residual_j = finite_value("capture.residual",
            (Dd::from_f64(next.supplied_wall_work_j)
                + Dd::from_f64(next.supplied_body_work_j)
                + Dd::from_f64(next.supplied_pressure_work_j)
                + Dd::from_f64(next.pressure_constraint_work_j)
                - Dd::from_f64(next.generated_heat_j)
                - next.kinetic_energy_change
                - next.backward_euler_loss).to_f64())?;
        let allowance = finite_value("capture.allowance", controls.power_tolerance_w * next.duration_s)?;
        if (kinetic_energy_change + Dd::from_f64(heat_j)).to_f64().abs() > allowance
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
        let pressure_work_axis = if self
            .config
            .faces
            .iter()
            .any(|face| face.prescribed_pressure().is_some())
        {
            Some(pressure::admit_pressure_work_heating(
                self.domain,
                self.config,
                controls,
            )?)
        } else {
            admit_heating(self.domain, self.config, controls)?;
            None
        };
        // Only momentum convection is omitted. Native rho*V/dt, pressure,
        // full symmetric stress and stationary Darcy coefficients remain.
        // The existing heat opt-in selects normal traction and zero tangential
        // trace at an opposed reservoir pair, exactly as in steady Stokes.
        self.creeping = true;
        self.pressure_work_axis = pressure_work_axis;
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
        let step = self.sample_unsteady_mechanical_step(
            start_s, end_s, buoyancy_force_n_m3, controls, gate,
        )?;
        step.admit_balance(controls)?;
        Ok(step)
    }

    /// Reassemble the original frozen native equations without solving,
    /// projecting or changing accepted time levels. The native shared scale
    /// prevents component cancellation from hiding an endpoint force defect.
    pub(in crate::conjugate) fn unsteady_material_momentum_residual(
        &mut self,
        gate: &CancelGate,
    ) -> Result<f64, ChtError> {
        self.coupled_darcy_residual(gate)
    }

    /// Sample original source coefficients at an unaccepted nonlinear
    /// endpoint. Physical Reynolds/partition gates still apply immediately;
    /// its kinetic balance is admitted before a receipt can be published.
    pub(in crate::conjugate) fn sample_unsteady_mechanical_step(
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
        if !self.creeping
            || time.scheme != TimeScheme::BackwardEuler
            || time.older.is_some()
            || time.face_mass_kg.is_some()
        {
            return Err(invalid_power(
                "transient heat requires one native constant-density backward-Euler Stokes step",
            ));
        }
        // Native pressure is already centered; never subtract the physical
        // component datum a second time or reconstruct a large absolute field.
        let mut power =
            self.sample_viscous_power::<false>(&self.field(), &self.pressure, controls, gate)?;
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
        let mut numerical_loss = Dd::ZERO;
        let mut kinetic_energy_change = Dd::ZERO;
        for a in 0..3 {
            for &index in &self.comps[a].unknowns {
                poll(gate)?;
                // Interior rows retain their exact original mass arithmetic.
                // A traction-mouth row mirrors one cell and is weighted by
                // one half in physical work, so its kinetic mass must likewise
                // use half of that row's rho*area*width inertia.
                let mass = if matches!(self.comps[a].kind[index], Kind::Outlet(_)) {
                    let (cell, _, weight) = self.power_face_cells(a, index);
                    weight
                        * (self.rho * self.comp_face_area(a, index) * self.domain.width(a, cell[a]))
                } else {
                    self.interior_face_mass_kg(a, index)
                };
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
                let increment = FaceKineticIncrement::new(mass, old, new);
                kinetic_energy_change = kinetic_energy_change + increment.stored_j();
                numerical_loss = numerical_loss + increment.backward_euler_loss_j();
            }
        }
        let kinetic_energy_change_j = finite_value(
            "unsteady_heat.kinetic_change", kinetic_energy_change.to_f64(),
        )?;
        let backward_euler_loss_j = finite_value(
            "unsteady_heat.numerical_loss", numerical_loss.to_f64(),
        )?;
        let residual_j = finite_value(
            "unsteady_heat.step_residual",
            (Dd::from_f64(time.dt) * Dd::from_f64(power.residual_w)
                - kinetic_energy_change - numerical_loss).to_f64(),
        )?;
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
            kinetic_energy_change_j,
            kinetic_energy_change,
            backward_euler_loss_j,
            backward_euler_loss: numerical_loss,
            residual_j,
        })
    }
}

#[cfg(test)]
mod audit_tests {
    use super::*;

    // Exact one-degree-of-freedom BE balances exercise the production
    // extensive account. Native mouth geometry is tested separately below.
    fn pressure_step(before: f64, force: f64, start: f64) -> UnsteadyMechanicalStep {
        let mass = 2.0;
        let dt = 0.25;
        let resistance = 0.5;
        let after = (mass * before + dt * force) / (mass + dt * resistance);
        let increment = FaceKineticIncrement::new(mass, before, after);
        let heat = resistance * after * after;
        let supplied = force * after;
        let kinetic = increment.stored_j();
        let numerical = increment.backward_euler_loss_j();
        let power = ViscousPower {
            cell_heat_w: vec![heat],
            heat_w: heat,
            molecular_heat_w: heat,
            darcy_heat_w: 0.0,
            forchheimer_heat_w: 0.0,
            cell_darcy_heat_w: None,
            cell_forchheimer_heat_w: None,
            wall_power_w: [0.0; 6],
            wall_tangential_force_n: [[0.0; 3]; 6],
            supplied_wall_power_w: 0.0,
            pressure_power_w: [0.0, supplied, 0.0, 0.0, 0.0, 0.0],
            supplied_pressure_power_w: supplied,
            cell_body_power_w: None,
            supplied_body_power_w: 0.0,
            pressure_constraint_power_w: 0.0,
            residual_w: supplied - heat,
            reynolds: 0.01,
        };
        UnsteadyMechanicalStep {
            start_s: start,
            duration_s: dt,
            end_s: start + dt,
            residual_j: (Dd::from_f64(dt) * Dd::from_f64(power.residual_w) - kinetic - numerical)
                .to_f64(),
            power,
            cell_buoyancy_power_w: None,
            supplied_buoyancy_power_w: 0.0,
            initial_kinetic_energy_j: 0.5 * mass * before * before,
            final_kinetic_energy_j: 0.5 * mass * after * after,
            kinetic_energy_change_j: kinetic.to_f64(),
            kinetic_energy_change: kinetic,
            backward_euler_loss_j: numerical.to_f64(),
            backward_euler_loss: numerical,
        }
    }

    #[test]
    fn signed_reservoir_work_closes_acceleration_braking_and_physical_heat() {
        let controls = ViscousHeatControls {
            maximum_reynolds: 0.1,
            power_tolerance_w: 1e-12,
        };
        let mut audit = UnsteadyHeatingAudit::default();
        let mut before = 0.0;
        let mut returned_work = false;
        let mut supplied = 0.0;
        for (step, force) in [2.0, -1.0, -2.0, 0.0].into_iter().enumerate() {
            let receipt = pressure_step(before, force, step as f64 * 0.25);
            receipt.admit_balance(controls).unwrap();
            returned_work |= receipt.power.supplied_pressure_power_w < 0.0;
            supplied += receipt.duration_s * receipt.power.supplied_pressure_power_w;
            audit.accept(&receipt, controls).unwrap();
            // Recover the signed accepted endpoint from its independent BE row.
            before = (2.0 * before + 0.25 * force) / 2.125;
        }
        assert!(returned_work);
        assert!(audit.generated_heat_j > 0.0);
        assert!(audit.backward_euler_loss_j > 0.0);
        assert_eq!(audit.supplied_pressure_work_j, supplied);
        assert_eq!(audit.pressure_work_j[1], supplied);
        assert_eq!(audit.supplied_wall_work_j, 0.0);
        assert_eq!(audit.supplied_body_work_j, 0.0);
        assert_eq!(audit.pressure_constraint_work_j, 0.0);
        assert!(
            (supplied
                - audit.generated_heat_j
                - audit.kinetic_energy_change_j
                - audit.backward_euler_loss_j)
                .abs()
                < 1e-14
        );
        assert_eq!(audit.generated_heat_j, audit.molecular_heat_j);
    }

    #[test]
    fn rejected_pressure_work_preserves_the_accepted_extensive_account() {
        let controls = ViscousHeatControls {
            maximum_reynolds: 0.1,
            power_tolerance_w: 1e-12,
        };
        let mut audit = UnsteadyHeatingAudit::default();
        let first = pressure_step(0.0, 2.0, 0.0);
        audit.accept(&first, controls).unwrap();
        let saved = audit;
        let mut wrong = pressure_step(0.5 / 2.125, -1.0, 0.25);
        wrong.power.supplied_pressure_power_w += 1.0;
        assert!(audit.accept(&wrong, controls).is_err());
        assert_eq!(audit, saved);

        let mut overflow = pressure_step(0.5 / 2.125, -1.0, 0.25);
        overflow.duration_s = 2.0;
        overflow.power.pressure_power_w[0] = f64::MAX;
        assert!(audit.accept(&overflow, controls).is_err());
        assert_eq!(audit, saved);
    }
}

#[cfg(test)]
mod pressure_mouth_tests {
    use super::*;

    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 2e-10,
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn native_pressure_mouth_inertia_and_regenerative_work_use_physical_half_cells() {
        let gate = CancelGate::new_clock_free();
        let fluid = FluidProperties {
            density_kg_m3: 2.0,
            specific_heat_j_kg_k: 2.0,
            conductivity_w_m_k: 0.01,
            kinematic_viscosity_m2_s: 1.0,
        };
        let controls = ViscousHeatControls {
            maximum_reynolds: 0.1,
            power_tolerance_w: 1e-9,
        };
        for axis in 0..3 {
            let mut widths = [vec![0.5], vec![0.75], vec![0.5]];
            widths[axis] = vec![0.5, 1.0, 1.5];
            let domain = VoxelDomain::graded(widths).unwrap();
            let area = domain.volume(0) / domain.width(axis, 0);
            let mass = fluid.density_kg_m3 * 3.0 * area;
            let mut faces = [FvBoundary::Symmetry; 6];
            faces[2 * axis] = FvBoundary::Outlet;
            faces[2 * axis + 1] = FvBoundary::Outlet;
            let mut config = SimpleConfig::new(faces);
            config.tolerance = 1e-11;
            config.pressure_tolerance = 1e-13;
            config.momentum_tolerance = 1e-12;
            let viscosity = vec![
                fluid.density_kg_m3 * fluid.kinematic_viscosity_m2_s;
                domain.cell_count()
            ];
            let run = |common: f64| {
                let mut solver =
                    material_solver(&domain, &fluid, &viscosity, None, None, &config, &gate).unwrap();
                solver.enable_unsteady_stokes_heat(controls).unwrap();
                let mut audit = UnsteadyHeatingAudit::default();
                let mut records = Vec::new();
                let mut expected = 0.0;
                for (step, head) in [0.0625, 0.0625, -0.0625, -0.125, -0.125]
                    .into_iter()
                    .enumerate()
                {
                    let previous = expected;
                    expected += 0.125 * head / (fluid.density_kg_m3 * 3.0);
                    let datum = common + step as f64 * (common / 8.0);
                    let mut offsets = [0.0; 6];
                    offsets[2 * axis] = datum + head;
                    offsets[2 * axis + 1] = datum;
                    solver.set_pressure_loads(1.0, offsets, &gate).unwrap();
                    solver.start_time(0.125, TimeScheme::BackwardEuler);
                    let mut settled = false;
                    for _ in 0..1000 {
                        let (mass, momentum) = solver.sweep(&gate).unwrap();
                        if mass <= 1e-11 && momentum <= 1e-11 {
                            settled = true;
                            break;
                        }
                    }
                    assert!(settled);
                    solver.project(&gate).unwrap();
                    let receipt = solver
                        .unsteady_mechanical_step(
                            step as f64 * 0.125,
                            (step + 1) as f64 * 0.125,
                            None,
                            controls,
                            &gate,
                        )
                        .unwrap();
                    for &face in &solver.comps[axis].unknowns {
                        near(solver.vel[axis][face], expected);
                    }
                    near(
                        receipt.initial_kinetic_energy_j,
                        0.5 * mass * previous * previous,
                    );
                    near(
                        receipt.final_kinetic_energy_j,
                        0.5 * mass * expected * expected,
                    );
                    near(
                        receipt.kinetic_energy_change_j,
                        0.5 * mass * (expected - previous) * (expected + previous),
                    );
                    near(
                        receipt.backward_euler_loss_j,
                        0.5 * mass * (expected - previous).powi(2),
                    );
                    near(
                        receipt.power.supplied_pressure_power_w,
                        head * area * expected,
                    );
                    near(receipt.power.heat_w, 0.0);
                    audit.accept(&receipt, controls).unwrap();
                    records.push(receipt);
                }
                assert!(records[2].power.supplied_pressure_power_w < 0.0);
                assert!(solver.vel[axis][solver.comps[axis].unknowns[0]] < 0.0);
                assert!(audit.backward_euler_loss_j > 0.0);
                near(audit.generated_heat_j, 0.0);
                (solver.vel, audit, records)
            };
            let baseline = run(0.0);
            let shifted = run(1_099_511_627_776.0);
            assert_eq!(baseline, shifted);
        }
    }
}
