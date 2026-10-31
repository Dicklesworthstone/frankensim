//! Simultaneous original mean-port closure for several finite DG ducts.
//!
//! The full FE cross-port matrix is retained. Only the small constant interface
//! matrix is factored cold; each interval uses the original FE trajectory and
//! the original DG private trials, including all finite return histories.

use super::*;
use fs_la::factor::{Lu, lu};

pub(super) struct PreparedMultipleDg {
    mobility: Vec<f64>,
    impedance: Vec<f64>,
    matrix: Vec<f64>,
    factor: Lu,
}

impl PreparedMultipleDg {
    fn new(mobility: Vec<f64>, impedance: Vec<f64>, gate: &CancelGate)
        -> Result<Self, ThermoviscousDgError>
    {
        poll(gate)?;
        let n = impedance.len();
        if n < 2 || n.checked_mul(n) != Some(mobility.len())
            || mobility.iter().any(|value| !value.is_finite())
            || impedance.iter().any(|value| !value.is_finite() || *value <= 0.0)
            || (0..n).any(|i| mobility[i * n + i] <= 0.0)
        {
            return Err(invalid("finite original multiport mobility and positive DG impedances required"));
        }
        // Q_FE = f + B p, Q_DG = -Q_FE, p = h + Z Q_DG.
        // Thus (I + Z B) p = h - Z f. B need not be symmetric.
        let mut matrix = mobility.clone();
        for (i, row) in matrix.chunks_exact_mut(n).enumerate() {
            poll(gate)?;
            for value in row.iter_mut() { *value *= impedance[i]; }
            row[i] += 1.0;
        }
        if matrix.iter().any(|value| !value.is_finite()) {
            return Err(invalid("original field/DG interface matrix overflowed"));
        }
        let factor = lu(&matrix, n).map_err(ThermoviscousTimeError::from)?;
        if !factor.growth.is_finite() {
            return Err(invalid("original field/DG interface factor is nonfinite"));
        }
        poll(gate)?;
        Ok(Self { mobility, impedance, matrix, factor })
    }

    fn pressures(&self, free: &[f64], mobility: &[f64], relations: &[(f64, f64)],
        tolerance: f64, gate: &CancelGate) -> Result<Vec<f64>, ThermoviscousDgError>
    {
        let n = self.impedance.len();
        if free.len() != n || relations.len() != n || mobility.len() != self.mobility.len()
            || mobility.iter().zip(&self.mobility).any(|(a, b)| a.to_bits() != b.to_bits())
            || relations.iter().zip(&self.impedance).any(|((h, z), cached)|
                !h.is_finite() || z.to_bits() != cached.to_bits())
        {
            return Err(invalid("original field/DG multiport response changed"));
        }
        let rhs: Vec<_> = relations.iter().zip(free).map(|((h, z), f)| h - z * f).collect();
        if rhs.iter().any(|value| !value.is_finite()) {
            return Err(invalid("original field/DG multiport drive is nonfinite"));
        }
        let mut pressure = rhs.clone();
        poll(gate)?;
        self.factor.solve(&mut pressure);
        if pressure.iter().any(|value| !value.is_finite()) {
            return Err(invalid("original field/DG multiport pressure is nonfinite"));
        }
        for (row, rhs) in self.matrix.chunks_exact(n).zip(rhs) {
            poll(gate)?;
            let value: f64 = row.iter().zip(&pressure).map(|(a, p)| a * p).sum();
            let scale = rhs.abs() + row.iter().zip(&pressure).map(|(a, p)| (a * p).abs()).sum::<f64>();
            if !value.is_finite() || !scale.is_finite()
                || (value - rhs).abs() > tolerance * scale.max(1e-300)
            {
                return Err(invalid("original field/DG multiport solve failed its residual admission"));
            }
        }
        Ok(pressure)
    }
}

impl<'a, 'c> ThermoviscousDg<'a, 'c> {
    /// Attach another original finite duct to an exclusive planar FE port.
    /// Attachment is cold: all ducts keep their authored initial states and
    /// no accepted trajectory may be reconfigured. The complete off-diagonal
    /// FE response couples every inlet, so a cork or airway return affects
    /// the bore through the same resolved field and its retained memories.
    pub fn with_additional_duct(mut self, port: usize, load: DgPortLoad, gate: &CancelGate)
        -> Result<Self, ThermoviscousDgError>
    {
        poll(gate)?;
        let area = self.source.inner().plane_wave_port_area_m2(port)
            .ok_or_else(|| invalid("additional DG connection requires an exclusive planar unit-pressure FE port"))?;
        if self.accepted_steps != 0 || self.ports.contains(&port)
            || load.state().step_count() != 0 || load.has_driven_mouth()
            || self.source.gas() != load.duct().inlet_gas()
            || self.prepared.duration_s().to_bits() != load.propagator().dt().to_bits()
            || (area - load.duct().inlet_area_m2()).abs() > 64.0 * f64::EPSILON * area
        {
            return Err(invalid("additional field/DG connection requires a distinct cold ordinary inlet and equal gas/area/clock"));
        }
        self.ports.push(port);
        self.additional.push(load);
        let zero = self.prepared.zero_state();
        let traction = vec![0.0; self.source.inner().port_mode_count()];
        let forces = vec![0.0; zero.wall_velocity_m_sqrt_kg_per_s.len()];
        let (free, mobility) = self.prepared.average_relations(
            &zero, &traction, &forces, &self.ports, &mut || gate.is_requested(),
        )?;
        if free.iter().any(|value| *value != 0.0) {
            return Err(invalid("unforced original field has nonzero mean port flow"));
        }
        let impedance = std::iter::once(&self.load).chain(&self.additional)
            .map(|load| load.port_response().map(|relation| relation.impedance_pa_s_m3))
            .collect::<Result<Vec<_>, _>>()?;
        self.multiple = Some(PreparedMultipleDg::new(mobility, impedance, gate)?);
        self.admit_state(self.state(), self.duct_storage_j())?;
        poll(gate)?;
        Ok(self)
    }

    pub(super) fn step_multiple(&mut self, traction_pa: &[f64], forces: &[f64], gate: &CancelGate)
        -> Result<&ThermoviscousDgFrame, ThermoviscousDgError>
    {
        poll(gate)?;
        if self.accepted_steps >= self.limits.maximum_steps
            || traction_pa.len() != self.source.inner().port_mode_count()
            || self.ports.iter().any(|&port| traction_pa.get(port).copied() != Some(0.0))
            || traction_pa.iter().any(|p| !p.is_finite() || p.abs() > self.limits.maximum_abs_pressure_pa)
        {
            return Err(invalid("bounded coupled clock and complete external tractions with zero internal pressure slots required"));
        }
        let (free, mobility) = self.prepared.average_relations(
            self.state(), traction_pa, forces, &self.ports, &mut || gate.is_requested(),
        )?;
        let relations = std::iter::once(&self.load).chain(&self.additional)
            .map(|load| load.port_response().map(|r| (r.history_pressure_pa, r.impedance_pa_s_m3)))
            .collect::<Result<Vec<_>, _>>()?;
        let pressure = self.multiple.as_ref()
            .ok_or_else(|| invalid("original field/DG multiport preparation is absent"))?
            .pressures(&free, &mobility, &relations, self.limits.interface_relative_tolerance, gate)?;
        if pressure.iter().any(|p| p.abs() > self.limits.maximum_abs_pressure_pa) {
            return Err(invalid("common resolved field/DG pressures exceed the physical budget"));
        }
        let mut applied = traction_pa.to_vec();
        for (&port, &pressure) in self.ports.iter().zip(&pressure) { applied[port] = pressure; }
        let field = self.prepared.step(self.state(), &applied, forces, &mut || gate.is_requested())?;
        let n = self.ports.len();
        let actual_flow: Vec<_> = self.ports.iter().map(|&port| -field.mean_inward_flow_m3_s[port]).collect();
        for (i, row) in mobility.chunks_exact(n).enumerate() {
            let flow = -free[i] - row.iter().zip(&pressure).map(|(b, p)| b * p).sum::<f64>();
            let scale = free[i].abs() + actual_flow[i].abs()
                + row.iter().zip(&pressure).map(|(b, p)| (b * p).abs()).sum::<f64>();
            if !flow.is_finite() || !actual_flow[i].is_finite() || !scale.is_finite()
                || (actual_flow[i] - flow).abs() > self.limits.interface_relative_tolerance * scale.max(1e-300)
            {
                return Err(invalid("resolved field and original multiport map disagree on swept volume"));
            }
        }
        let dt = self.prepared.duration_s();
        let authored_external_work_j = traction_pa.iter().zip(&field.mean_inward_flow_m3_s)
            .map(|(p, q)| p * (q * dt)).sum::<f64>() + field.work.structural_work_j;
        let inlets: Vec<_> = std::iter::once(&self.load).chain(&self.additional).map(DgPortLoad::inlet).collect();
        // Every trial borrows a distinct original owner. An error while making
        // or checking any trial drops all of them without publishing a state.
        let trials = std::iter::once(&mut self.load).chain(self.additional.iter_mut())
            .zip(&actual_flow).map(|(load, &flow)| {
                poll(gate)?;
                Ok(load.preview_step(flow)?)
            }).collect::<Result<Vec<_>, ThermoviscousDgError>>()?;
        let mut stored_energy_j = field.work.final_energy_j;
        let mut storage_change_j = field.work.delta_energy_j;
        let mut material_dissipation_j = field.work.dissipated_j();
        let mut terminal_work_j = 0.0;
        let mut interfaces = Vec::with_capacity(n);
        let mut scale = field.work.initial_energy_j.abs() + field.work.final_energy_j.abs()
            + authored_external_work_j.abs() + material_dissipation_j.abs();
        for (i, trial) in trials.iter().enumerate() {
            poll(gate)?;
            let duct = trial.output();
            let port = duct.port_outputs.iter().find(|port| port.id == inlets[i])
                .ok_or_else(|| invalid("original DG candidate omitted its connected inlet"))?;
            let (history, impedance) = relations[i];
            let pressure_scale = pressure[i].abs() + history.abs() + (impedance * actual_flow[i]).abs();
            if !port.conjugate_rate.is_finite() || !pressure_scale.is_finite()
                || (port.conjugate_rate - pressure[i]).abs()
                    > self.limits.interface_relative_tolerance * pressure_scale.max(1e-300)
            {
                return Err(invalid("resolved field and DG candidates disagree on common multiport pressure"));
            }
            let field_work = pressure[i] * (field.mean_inward_flow_m3_s[self.ports[i]] * dt);
            let interface = trial.inlet_work_j() + field_work;
            let interface_tolerance = self.limits.energy_absolute_tolerance_j
                + self.limits.energy_relative_tolerance * (trial.inlet_work_j().abs() + field_work.abs());
            if !interface.is_finite() || !interface_tolerance.is_finite() || interface.abs() > interface_tolerance {
                return Err(invalid("one original field/DG interface failed its independent work cancellation"));
            }
            interfaces.push(interface);
            let previous = trial.previous_stored_energy_j();
            let stored = duct.ledger.stored_energy_j;
            let loss = duct.ledger.cumulative_dissipated_j - trial.previous_cumulative_dissipated_j();
            let terminal = sum_terminal_work(trial.terminal_frames());
            stored_energy_j += stored;
            storage_change_j += stored - previous;
            material_dissipation_j += loss;
            terminal_work_j += terminal;
            scale += previous.abs() + stored.abs() + loss.abs() + terminal.abs();
        }
        let interface_work_residual_j = interfaces.iter().sum::<f64>();
        let balance_residual_j = storage_change_j + material_dissipation_j
            + terminal_work_j - authored_external_work_j;
        let balance_tolerance_j = self.limits.energy_absolute_tolerance_j
            + self.limits.energy_relative_tolerance * scale;
        if [stored_energy_j, storage_change_j, material_dissipation_j, terminal_work_j,
            authored_external_work_j, interface_work_residual_j, balance_residual_j, scale, balance_tolerance_j]
            .iter().any(|value| !value.is_finite())
            || stored_energy_j > self.limits.maximum_total_energy_j
            || field.state.pressure_pa.iter().any(|p| p.abs() > self.limits.maximum_abs_pressure_pa)
            || material_dissipation_j < -balance_tolerance_j
            || interface_work_residual_j.abs() > balance_tolerance_j
            || balance_residual_j.abs() > balance_tolerance_j
        {
            return Err(invalid("original resolved field/multiple DG pressure, storage or same-window work admission refused"));
        }
        let mut additional_ducts = Vec::with_capacity(n - 1);
        poll(gate)?;
        // All remaining publications are infallible. There is no cancellation
        // point or fallible check after the first original owner is accepted.
        let mut trials = trials.into_iter();
        let duct = trials.next().expect("primary original DG trial").accept();
        for (i, trial) in trials.enumerate() {
            additional_ducts.push(ThermoviscousDgPortFrame {
                port: self.ports[i + 1], duct: trial.accept(),
                mean_pressure_pa: pressure[i + 1], mean_flow_m3_s: actual_flow[i + 1],
                interface_work_residual_j: interfaces[i + 1],
            });
        }
        self.accepted_steps += 1;
        self.frame = Some(ThermoviscousDgFrame { step: self.accepted_steps, field, duct,
            additional_ducts, mean_pressure_pa: pressure[0], mean_flow_m3_s: actual_flow[0],
            stored_energy_j, storage_change_j, authored_external_work_j, material_dissipation_j,
            terminal_work_j, interface_work_residual_j, balance_residual_j, balance_tolerance_j });
        Ok(self.frame.as_ref().expect("joint field/multiple DG step accepted"))
    }
}
