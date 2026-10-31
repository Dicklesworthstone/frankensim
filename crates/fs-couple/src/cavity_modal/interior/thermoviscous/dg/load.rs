//! Original complete field and finite ducts behind one external acoustic port,
//! optionally sharing the primary duct's independently driven exterior mouths.
//!
//! Internal duct pressures are eliminated through the full FE mean-flow map.
//! The resulting `p = H + Z Q` is a same-window port relation, not a fitted
//! resonator: every original field and duct coordinate remains in its owner.
//! Positive Q enters the external FE port. All other external tractions and
//! structural forces are zero. The original mean-port partition still needs
//! independent interface, spatial and coupling-clock refinement.

use super::*;
use crate::bernoulli_aperture::radiation::environment::multiport::bake::mixed::DrivenMouthSurface;
use crate::dg_port::{DgIndependentMouthsFrame, DgIndependentMouthsResponse, DgPortTrial};
use fs_duct::dg::DgPortResponse;
use fs_la::factor::{Lu, lu};

/// A spatial thermoviscous region and all of its original attached ducts,
/// exposed at a distinct planar unit-pressure FE opening. This is an acoustic
/// load; it does not accept mean material mass or enthalpy transport.
pub struct ThermoviscousDgLoad<'a, 'c> {
    system: ThermoviscousDg<'a, 'c>,
    port: usize,
    area_m2: f64,
    selected: Vec<usize>,
    mobility: Vec<f64>,
    impedance: Vec<f64>,
    matrix: Vec<f64>,
    factor: Lu,
    /// Internal pressures at zero mouth flow are
    /// `free_duct_pressure - pressure_slope * p_ext`.
    pressure_slope: Vec<f64>,
    effective_mobility: f64,
    zero_traction: Vec<f64>,
    zero_forces: Vec<f64>,
    mouths: Option<PreparedMouths>,
}

/// Original primary-DG mouth coordinates after eliminating every FE/DG
/// interface. Only fixed mean maps are cached; all histories stay original.
struct PreparedMouths {
    original_impedance: Vec<f64>,
    /// Row-major internal held pressures per original outward mouth flow.
    internal_pressure_per_mouth: Vec<f64>,
    /// Complete [external FE flow, mouth flows] mean effort response.
    impedance: Vec<f64>,
    count: usize,
}

struct PortHistory {
    free_flow: Vec<f64>,
    ducts: Vec<DgPortResponse>,
    free_duct_pressure: Vec<f64>,
    response: DgPortResponse,
    mouths: Option<DgIndependentMouthsResponse>,
}

fn duct_responses(system: &ThermoviscousDg<'_, '_>)
    -> Result<(Vec<DgPortResponse>, Option<DgIndependentMouthsResponse>), ThermoviscousDgError>
{
    let mouths = system.load.has_driven_mouth()
        .then(|| system.load.independent_mouths_response()).transpose()?;
    let primary = if let Some(response) = &mouths {
        let width = response.history_pressure_pa.len();
        if width < 2 || width.checked_mul(width) != Some(response.impedance_pa_s_m3.len())
            || response.history_pressure_pa.iter().chain(&response.impedance_pa_s_m3)
                .any(|value| !value.is_finite())
        {
            return Err(invalid("resolved acoustic mouths lost their original complete pressure basis"));
        }
        DgPortResponse { history_pressure_pa: response.history_pressure_pa[0],
            impedance_pa_s_m3: response.impedance_pa_s_m3[0] }
    } else { system.load.port_response()? };
    let mut ducts = vec![primary];
    for load in &system.additional {
        if load.has_driven_mouth() {
            return Err(invalid("resolved mouths must share the primary duct's one exterior owner"));
        }
        ducts.push(load.port_response()?);
    }
    Ok((ducts, mouths))
}

fn solve_checked(matrix: &[f64], factor: &Lu, rhs: &[f64], tolerance: f64,
    gate: &CancelGate) -> Result<Vec<f64>, ThermoviscousDgError>
{
    poll(gate)?;
    if rhs.iter().any(|value| !value.is_finite()) {
        return Err(invalid("resolved acoustic port has a nonfinite original drive"));
    }
    let mut solution = rhs.to_vec();
    factor.solve(&mut solution);
    if solution.iter().any(|value| !value.is_finite()) {
        return Err(invalid("resolved acoustic port has a nonfinite internal pressure"));
    }
    for (row, &rhs) in matrix.chunks_exact(rhs.len()).zip(rhs) {
        poll(gate)?;
        let value: f64 = row.iter().zip(&solution).map(|(a, p)| a * p).sum();
        let scale = rhs.abs() + row.iter().zip(&solution).map(|(a, p)| (a * p).abs()).sum::<f64>();
        if !value.is_finite() || !scale.is_finite()
            || (value - rhs).abs() > tolerance * scale.max(1e-300)
        {
            return Err(invalid("resolved acoustic port failed its original internal solve residual"));
        }
    }
    Ok(solution)
}

impl<'a, 'c> ThermoviscousDgLoad<'a, 'c> {
    /// Select one original external opening before the first accepted interval.
    /// All attached ducts, authored initial fields and original material/thermal
    /// histories remain owned by `system`. Other unconnected pressures and all
    /// structural forces are held at zero by this one-port load.
    pub fn new(system: ThermoviscousDg<'a, 'c>, port: usize, gate: &CancelGate)
        -> Result<Self, ThermoviscousDgError>
    {
        poll(gate)?;
        if system.load.has_driven_mouth() {
            return Err(invalid("attach resolved exterior mouths explicitly after the complete field load"));
        }
        Self::prepare(system, port, gate)
    }

    fn prepare(system: ThermoviscousDg<'a, 'c>, port: usize, gate: &CancelGate)
        -> Result<Self, ThermoviscousDgError>
    {
        poll(gate)?;
        let area_m2 = system.source.inner().plane_wave_port_area_m2(port)
            .ok_or_else(|| invalid("resolved acoustic load requires an exclusive planar unit-pressure FE port"))?;
        if system.accepted_steps != 0 || system.ports.contains(&port) {
            return Err(invalid("resolved acoustic load requires a cold, unconnected external FE port"));
        }
        let mut selected = vec![port];
        selected.extend_from_slice(&system.ports);
        let zero = system.prepared.zero_state();
        let zero_traction = vec![0.0; system.source.inner().port_mode_count()];
        let zero_forces = vec![0.0; zero.wall_velocity_m_sqrt_kg_per_s.len()];
        let (free, mobility) = system.prepared.average_relations(
            &zero, &zero_traction, &zero_forces, &selected, &mut || gate.is_requested(),
        )?;
        let n = system.ports.len();
        let width = selected.len();
        let (ducts, mouth_response) = duct_responses(&system)?;
        let impedance: Vec<_> = ducts.iter().map(|response| response.impedance_pa_s_m3).collect();
        if free.iter().any(|value| *value != 0.0)
            || mobility.iter().any(|value| !value.is_finite())
            || impedance.iter().any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(invalid("resolved acoustic load requires finite original mobility and positive duct impedances"));
        }
        // Q_K = f_K + B_KK p_K + B_Ks p_s and p_K = h_K - Z_K Q_K.
        // Hence M p_K = h_K - Z_K f_K - Z_K B_Ks p_s, M = I + Z_K B_KK.
        let mut matrix = vec![0.0; n * n];
        let mut rhs = vec![0.0; n];
        for i in 0..n {
            poll(gate)?;
            for j in 0..n {
                matrix[i * n + j] = impedance[i] * mobility[(i + 1) * width + j + 1];
            }
            matrix[i * n + i] += 1.0;
            rhs[i] = impedance[i] * mobility[(i + 1) * width];
        }
        if matrix.iter().any(|value| !value.is_finite()) {
            return Err(invalid("resolved acoustic load interface matrix overflowed"));
        }
        let factor = lu(&matrix, n).map_err(ThermoviscousTimeError::from)?;
        if !factor.growth.is_finite() {
            return Err(invalid("resolved acoustic load has a nonfinite interface factor"));
        }
        let pressure_slope = solve_checked(&matrix, &factor, &rhs,
            system.limits.interface_relative_tolerance, gate)?;
        // The Schur complement retains every cross-port reaction. It is the
        // actual mean external flow per held external pressure, not B_ss alone.
        let effective_mobility = mobility[0] - mobility[1..width].iter()
            .zip(&pressure_slope).map(|(b, p)| b * p).sum::<f64>();
        if !effective_mobility.is_finite() || effective_mobility <= 0.0
            || !(1.0 / effective_mobility).is_finite()
        {
            return Err(invalid("resolved acoustic load requires strictly positive finite effective mobility"));
        }
        let mouths = if let Some(original) = mouth_response {
            let mouth_count = original.history_pressure_pa.len() - 1;
            let response_width = mouth_count + 1;
            let mut complete = vec![0.0; response_width * response_width];
            let mut internal_pressure_per_mouth = vec![0.0; n * mouth_count];
            complete[0] = 1.0 / effective_mobility;
            // The primary DG's actual inward flow is -Q_FE,K0. Derive its
            // coefficients from the complete FE row, preserving cancellation
            // and every finite passage's reaction rather than fitting a load.
            let primary_row = &mobility[width..2 * width];
            let primary_per_external = -(primary_row[0] - primary_row[1..].iter()
                .zip(&pressure_slope).map(|(b, slope)| b * slope).sum::<f64>())
                / effective_mobility;
            for row in 1..response_width {
                complete[row * response_width] = original.impedance_pa_s_m3[row * response_width]
                    * primary_per_external;
            }
            for mouth in 0..mouth_count {
                let mut rhs = vec![0.0; n];
                rhs[0] = original.impedance_pa_s_m3[mouth + 1];
                let column = solve_checked(&matrix, &factor, &rhs,
                    system.limits.interface_relative_tolerance, gate)?;
                let external = -mobility[1..width].iter().zip(&column)
                    .map(|(b, value)| b * value).sum::<f64>() / effective_mobility;
                complete[mouth + 1] = external;
                let mut primary_per_mouth = -primary_row[0] * external;
                for i in 0..n {
                    internal_pressure_per_mouth[i * mouth_count + mouth] = column[i];
                    primary_per_mouth -= primary_row[i + 1]
                        * (column[i] - pressure_slope[i] * external);
                }
                for row in 1..response_width {
                    complete[row * response_width + mouth + 1] =
                        original.impedance_pa_s_m3[row * response_width + mouth + 1]
                        + original.impedance_pa_s_m3[row * response_width] * primary_per_mouth;
                }
            }
            if complete.iter().chain(&internal_pressure_per_mouth).any(|value| !value.is_finite()) {
                return Err(invalid("resolved acoustic mouths overflowed the complete original response"));
            }
            Some(PreparedMouths { original_impedance: original.impedance_pa_s_m3,
                internal_pressure_per_mouth, impedance: complete, count: mouth_count })
        } else { None };
        poll(gate)?;
        Ok(Self { system, port, area_m2, selected, mobility, impedance, matrix,
            factor, pressure_slope, effective_mobility, zero_traction, zero_forces, mouths })
    }

    /// Attach independently driven rectangular mouths to the primary duct's
    /// original computed exterior, retaining its circular outlets and every
    /// self/mutual term. All FE/DG states survive this cold preparation. The
    /// source adapter must admit aperture geometry, gas and work projection.
    pub fn with_independent_mouths(self, surfaces: &[DrivenMouthSurface],
        mechanics_rate_hz: u32, gate: &CancelGate) -> Result<Self, ThermoviscousDgError>
    {
        poll(gate)?;
        if self.system.accepted_steps != 0 || self.mouths.is_some() || surfaces.is_empty() {
            return Err(invalid("resolved mouths require one nonempty cold exterior attachment"));
        }
        let mut system = self.system;
        system.load = system.load.with_independent_mouths(surfaces, mechanics_rate_hz, gate)
            .map_err(|error| if gate.is_requested() { ThermoviscousDgError::Cancelled }
                else { ThermoviscousDgError::Duct(AcousticRealizeError::Nonlinear(error.to_string())) })?;
        Self::prepare(system, self.port, gate)
    }

    pub const fn system(&self) -> &ThermoviscousDg<'a, 'c> { &self.system }
    pub const fn port(&self) -> usize { self.port }
    pub const fn area_m2(&self) -> f64 { self.area_m2 }
    pub fn duration_s(&self) -> f64 { self.system.prepared.duration_s() }
    /// Number of independently driven mouths sharing the original exterior.
    pub const fn mouth_count(&self) -> usize {
        match &self.mouths { Some(mouths) => mouths.count, None => 0 }
    }

    /// Exact affine relation of the original mean-port partition at this state.
    /// Positive flow enters the external FE opening. No owner is advanced.
    pub fn port_response(&self, gate: &CancelGate) -> Result<DgPortResponse, ThermoviscousDgError> {
        poll(gate)?;
        if self.mouths.is_some() {
            return Err(invalid("resolved exterior mouths require every independent flow coordinate"));
        }
        Ok(self.history(gate)?.response)
    }

    /// Mean efforts for [flow into the external FE port, outward mouth flows].
    /// Mouth order is the original shared radiation owner's attachment order.
    /// This preserves the full field/duct/exterior history and all cross terms.
    pub fn independent_mouths_response(&self, gate: &CancelGate)
        -> Result<DgIndependentMouthsResponse, ThermoviscousDgError>
    {
        poll(gate)?;
        let prepared = self.mouths.as_ref()
            .ok_or_else(|| invalid("resolved exterior mouth coordinates are not attached"))?;
        let history = self.history(gate)?;
        let original = history.mouths.as_ref().expect("admitted mouth history");
        let pressure = history.response.history_pressure_pa;
        let internal: Vec<_> = history.free_duct_pressure.iter().zip(&self.pressure_slope)
            .map(|(free, slope)| free - slope * pressure).collect();
        let width = self.selected.len();
        let row = &self.mobility[width..2 * width];
        let primary_flow = -history.free_flow[1] - row[0] * pressure
            - row[1..].iter().zip(&internal).map(|(b, p)| b * p).sum::<f64>();
        let mut history_pressure_pa = vec![pressure];
        let response_width = prepared.count + 1;
        for mouth in 0..prepared.count {
            history_pressure_pa.push(original.history_pressure_pa[mouth + 1]
                + original.impedance_pa_s_m3[(mouth + 1) * response_width] * primary_flow);
        }
        if history_pressure_pa.iter().any(|value| !value.is_finite()) {
            return Err(invalid("resolved acoustic mouth history pressure overflowed"));
        }
        Ok(DgIndependentMouthsResponse { history_pressure_pa,
            impedance_pa_s_m3: prepared.impedance.clone() })
    }

    fn history(&self, gate: &CancelGate) -> Result<PortHistory, ThermoviscousDgError> {
        let (free_flow, mobility) = self.system.prepared.average_relations(
            self.system.state(), &self.zero_traction, &self.zero_forces, &self.selected,
            &mut || gate.is_requested(),
        )?;
        let (ducts, mouths) = duct_responses(&self.system)?;
        let mouths_match = match (&self.mouths, &mouths) {
            (None, None) => true,
            (Some(prepared), Some(original)) => original.history_pressure_pa.len() == prepared.count + 1
                && original.impedance_pa_s_m3.len() == prepared.original_impedance.len()
                && original.impedance_pa_s_m3.iter().zip(&prepared.original_impedance)
                    .all(|(a, b)| a.to_bits() == b.to_bits()),
            _ => false,
        };
        if mobility.iter().zip(&self.mobility).any(|(a, b)| a.to_bits() != b.to_bits())
            || !mouths_match
            || ducts.iter().zip(&self.impedance).any(|(response, z)|
                !response.history_pressure_pa.is_finite() || response.impedance_pa_s_m3.to_bits() != z.to_bits())
        {
            return Err(invalid("resolved acoustic load original response changed after cold preparation"));
        }
        let rhs: Vec<_> = ducts.iter().zip(&free_flow[1..]).map(|(response, f)|
            response.history_pressure_pa - response.impedance_pa_s_m3 * f).collect();
        let free_duct_pressure = solve_checked(&self.matrix, &self.factor, &rhs,
            self.system.limits.interface_relative_tolerance, gate)?;
        let free = free_flow[0] + self.mobility[1..self.selected.len()].iter()
            .zip(&free_duct_pressure).map(|(b, p)| b * p).sum::<f64>();
        let response = DgPortResponse {
            history_pressure_pa: -free / self.effective_mobility,
            impedance_pa_s_m3: 1.0 / self.effective_mobility,
        };
        if !response.history_pressure_pa.is_finite() {
            return Err(invalid("resolved acoustic load history pressure overflowed"));
        }
        poll(gate)?;
        Ok(PortHistory { free_flow, ducts, free_duct_pressure, response, mouths })
    }

    /// Propagate private candidates at the requested external mean inward flow.
    /// The surrounding source owner must check its own interface projection and
    /// cancellation before calling the trial's infallible `accept`.
    pub(crate) fn preview_step(&mut self, flow_m3_s: f64, gate: &CancelGate)
        -> Result<ThermoviscousDgLoadTrial<'_>, ThermoviscousDgError>
    {
        self.preview_step_inner(flow_m3_s, None, gate)
    }

    /// Stage the complete field, every finite duct and the original shared
    /// exterior at independent external-FE and outward-mouth mean flows.
    pub(crate) fn preview_step_with_mouth_flows(&mut self, flow_m3_s: f64,
        mouth_flows_m3_s: &[f64], gate: &CancelGate)
        -> Result<ThermoviscousDgLoadTrial<'_>, ThermoviscousDgError>
    {
        self.preview_step_inner(flow_m3_s, Some(mouth_flows_m3_s), gate)
    }

    fn preview_step_inner(&mut self, flow_m3_s: f64, mouths: Option<&[f64]>, gate: &CancelGate)
        -> Result<ThermoviscousDgLoadTrial<'_>, ThermoviscousDgError>
    {
        poll(gate)?;
        let mouth_flows = match (&self.mouths, mouths) {
            (None, None) => &[][..],
            (Some(prepared), Some(flows)) if flows.len() == prepared.count
                && flows.iter().all(|value| value.is_finite()) => flows,
            _ => return Err(invalid("resolved exterior mouths require every finite independent flow")),
        };
        let limits = self.system.limits;
        if !flow_m3_s.is_finite() || self.system.accepted_steps >= limits.maximum_steps {
            return Err(invalid("resolved acoustic load requires finite flow and an unexhausted clock"));
        }
        let history = self.history(gate)?;
        let mut pressure = history.response.history_pressure_pa
            + history.response.impedance_pa_s_m3 * flow_m3_s;
        if let Some(prepared) = &self.mouths {
            for (mouth, flow) in mouth_flows.iter().enumerate() {
                pressure += prepared.impedance[mouth + 1] * flow;
            }
        }
        let internal: Vec<_> = history.free_duct_pressure.iter().zip(&self.pressure_slope)
            .enumerate().map(|(i, (p, slope))| {
                let mut value = p - slope * pressure;
                if let Some(prepared) = &self.mouths {
                    for (mouth, flow) in mouth_flows.iter().enumerate() {
                        value += prepared.internal_pressure_per_mouth[i * prepared.count + mouth] * flow;
                    }
                }
                value
            }).collect();
        if !pressure.is_finite() || pressure.abs() > limits.maximum_abs_pressure_pa
            || internal.iter().any(|p| !p.is_finite() || p.abs() > limits.maximum_abs_pressure_pa)
        {
            return Err(invalid("resolved acoustic load pressures exceed the original physical budget"));
        }
        let mut applied = self.zero_traction.clone();
        applied[self.port] = pressure;
        for (&port, &p) in self.system.ports.iter().zip(&internal) { applied[port] = p; }
        let field = self.system.prepared.step(self.system.state(), &applied, &self.zero_forces,
            &mut || gate.is_requested())?;
        let pressures: Vec<_> = std::iter::once(pressure).chain(internal.iter().copied()).collect();
        for (i, &port) in self.selected.iter().enumerate() {
            let row = &self.mobility[i * pressures.len()..(i + 1) * pressures.len()];
            let expected = history.free_flow[i] + row.iter().zip(&pressures).map(|(b, p)| b * p).sum::<f64>();
            let actual = field.mean_inward_flow_m3_s[port];
            let scale = history.free_flow[i].abs() + actual.abs() + flow_m3_s.abs()
                + row.iter().zip(&pressures).map(|(b, p)| (b * p).abs()).sum::<f64>();
            if !expected.is_finite() || !actual.is_finite() || !scale.is_finite()
                || (actual - expected).abs() > limits.interface_relative_tolerance * scale.max(1e-300)
                || (i == 0 && (actual - flow_m3_s).abs() > limits.interface_relative_tolerance * scale.max(1e-300))
            {
                return Err(invalid("resolved acoustic load and original field disagree on mean swept volume"));
            }
        }
        let dt = self.duration_s();
        let inlet_work_j = pressure * (field.mean_inward_flow_m3_s[self.port] * dt);
        let actual_flow: Vec<_> = self.system.ports.iter()
            .map(|&port| -field.mean_inward_flow_m3_s[port]).collect();
        let inlets: Vec<_> = std::iter::once(&self.system.load).chain(&self.system.additional)
            .map(DgPortLoad::inlet).collect();
        let trials = std::iter::once(&mut self.system.load).chain(self.system.additional.iter_mut())
            .zip(&actual_flow).enumerate().map(|(i, (load, &flow))| {
                poll(gate)?;
                Ok(if i == 0 && !mouth_flows.is_empty() {
                    load.preview_step_with_mouth_flows(flow, mouth_flows)?
                } else { load.preview_step(flow)? })
            }).collect::<Result<Vec<_>, ThermoviscousDgError>>()?;
        let mut stored_energy_j = field.work.final_energy_j;
        let mut storage_change_j = field.work.delta_energy_j;
        let mut material_dissipation_j = field.work.dissipated_j();
        let mut terminal_work_j = 0.0;
        let mut interfaces = Vec::with_capacity(trials.len());
        let mut scale = field.work.initial_energy_j.abs() + field.work.final_energy_j.abs()
            + inlet_work_j.abs() + material_dissipation_j.abs();
        for (i, trial) in trials.iter().enumerate() {
            poll(gate)?;
            let output = trial.output();
            let actual = output.port_outputs.iter().find(|port| port.id == inlets[i])
                .ok_or_else(|| invalid("original DG candidate lost its resolved-field inlet"))?.conjugate_rate;
            let relation = history.ducts[i];
            let mut pressure_scale = internal[i].abs() + relation.history_pressure_pa.abs()
                + (relation.impedance_pa_s_m3 * actual_flow[i]).abs();
            if i == 0 {
                if let Some(original) = &history.mouths {
                    pressure_scale += original.impedance_pa_s_m3[1..mouth_flows.len() + 1]
                        .iter().zip(mouth_flows).map(|(z, q)| (z * q).abs()).sum::<f64>();
                    let receipt = trial.independent_mouths_frame().ok_or_else(||
                        invalid("resolved acoustic mouth candidate lost its original exterior receipt"))?;
                    if receipt.mouths.len() != mouth_flows.len()
                        || receipt.mouths.iter().any(|mouth| !mouth.mean_pressure_pa.is_finite()
                            || mouth.mean_pressure_pa.abs() > limits.maximum_abs_pressure_pa)
                    {
                        return Err(invalid("resolved acoustic mouth pressure exceeds the original physical budget"));
                    }
                }
            }
            if !actual.is_finite() || !pressure_scale.is_finite()
                || (actual - internal[i]).abs() > limits.interface_relative_tolerance * pressure_scale.max(1e-300)
            {
                return Err(invalid("resolved acoustic load and original DG disagree on internal pressure"));
            }
            let field_work = internal[i] * (field.mean_inward_flow_m3_s[self.system.ports[i]] * dt);
            let residual = trial.inlet_work_j() + field_work;
            let allowance = limits.energy_absolute_tolerance_j
                + limits.energy_relative_tolerance * (trial.inlet_work_j().abs() + field_work.abs());
            if !residual.is_finite() || !allowance.is_finite() || residual.abs() > allowance {
                return Err(invalid("resolved acoustic load failed an internal field/duct work exchange"));
            }
            interfaces.push(residual);
            let previous = trial.previous_stored_energy_j();
            let stored = output.ledger.stored_energy_j;
            let loss = output.ledger.cumulative_dissipated_j - trial.previous_cumulative_dissipated_j();
            let terminal = sum_terminal_work(trial.terminal_frames());
            stored_energy_j += stored;
            storage_change_j += stored - previous;
            material_dissipation_j += loss;
            terminal_work_j += terminal;
            scale += previous.abs() + stored.abs() + loss.abs() + terminal.abs();
        }
        let interface_work_residual_j = interfaces.iter().sum::<f64>();
        let balance_residual_j = storage_change_j + material_dissipation_j + terminal_work_j - inlet_work_j;
        let balance_tolerance_j = limits.energy_absolute_tolerance_j + limits.energy_relative_tolerance * scale;
        if [stored_energy_j, storage_change_j, material_dissipation_j, terminal_work_j,
            inlet_work_j, interface_work_residual_j, balance_residual_j, scale, balance_tolerance_j]
            .iter().any(|value| !value.is_finite())
            || stored_energy_j > limits.maximum_total_energy_j
            || field.state.pressure_pa.iter().any(|p| p.abs() > limits.maximum_abs_pressure_pa)
            || material_dissipation_j < -balance_tolerance_j
            || interface_work_residual_j.abs() > balance_tolerance_j
            || balance_residual_j.abs() > balance_tolerance_j
        {
            return Err(invalid("resolved acoustic load pressure, storage or complete work admission refused"));
        }
        let additional_ducts = trials.iter().enumerate().skip(1).map(|(i, trial)|
            ThermoviscousDgPortFrame { port: self.system.ports[i], duct: trial.output().clone(),
                mean_pressure_pa: internal[i], mean_flow_m3_s: actual_flow[i],
                interface_work_residual_j: interfaces[i] }).collect();
        let frame = ThermoviscousDgFrame {
            step: self.system.accepted_steps + 1, field, duct: trials[0].output().clone(),
            additional_ducts, mean_pressure_pa: internal[0], mean_flow_m3_s: actual_flow[0],
            stored_energy_j, storage_change_j, authored_external_work_j: inlet_work_j,
            material_dissipation_j, terminal_work_j, interface_work_residual_j,
            balance_residual_j, balance_tolerance_j,
        };
        poll(gate)?;
        Ok(ThermoviscousDgLoadTrial { accepted_steps: &mut self.system.accepted_steps,
            accepted_frame: &mut self.system.frame, trials, frame, pressure, inlet_work_j })
    }
}

/// Private complete-field and duct candidates. Drop leaves all owners untouched.
pub(crate) struct ThermoviscousDgLoadTrial<'a> {
    accepted_steps: &'a mut u64,
    accepted_frame: &'a mut Option<ThermoviscousDgFrame>,
    trials: Vec<DgPortTrial<'a>>,
    frame: ThermoviscousDgFrame,
    pressure: f64,
    inlet_work_j: f64,
}
impl ThermoviscousDgLoadTrial<'_> {
    pub(crate) const fn mean_pressure_pa(&self) -> f64 { self.pressure }
    pub(crate) const fn inlet_work_j(&self) -> f64 { self.inlet_work_j }
    pub(crate) const fn frame(&self) -> &ThermoviscousDgFrame { &self.frame }
    /// Original primary-DG/exterior receipt. Its inlet pressure/work belong to
    /// that internal DG inlet, while `mean_pressure_pa`/`inlet_work_j` above
    /// belong to the external FE opening. Mouth work enters the exterior
    /// directly; it is not FE/DG storage or an extra debit of FE inlet work.
    pub(crate) fn independent_mouths_frame(&self) -> Option<&DgIndependentMouthsFrame> {
        self.trials[0].independent_mouths_frame()
    }

    /// All checks, allocation and cancellation belong before this publication.
    pub(crate) fn accept(self) {
        for trial in self.trials { trial.accept(); }
        *self.accepted_steps = self.frame.step;
        *self.accepted_frame = Some(self.frame);
    }
}
