use super::super::evolve::{forward, operator, totals};
use super::super::{
    Boundaries, Conserved, Error, Reconstruction, StepControl, StepReport, Totals, poll, positive,
};
use super::{Accumulator, BranchEnd, GasNetwork1d, NetworkStepReport, NetworkTotals, PortExchange};
use fs_qty::Time;

impl GasNetwork1d {
    /// Combined mass and energy in deterministic branch-then-plenum order.
    pub fn totals(&self) -> NetworkTotals {
        self.current_totals(&mut || false)
            .expect("admitted network storage and an infallible callback cannot refuse")
    }

    fn current_totals(&self, cancelled: &mut impl FnMut() -> bool) -> Result<NetworkTotals, Error> {
        let mut sum = Accumulator::default();
        for branch in &self.branches {
            let value = totals(&branch.states, branch.mesh.volumes(), cancelled)?;
            sum.add(value.mass, value.energy)?;
        }
        for plenum in &self.plenums {
            poll(cancelled)?;
            sum.add(
                plenum.volume_m3 * plenum.state.density,
                plenum.volume_m3 * plenum.state.energy,
            )?;
        }
        Ok(sum.totals())
    }

    /// Advance all owners by one accepted common SSP RK2 step.
    ///
    /// Actual port wave speeds limit both branch and compact-volume CFL. A trial
    /// failure halves the whole network step; no branch may publish independently.
    /// Numerical/source failures use the same first-order retry policy as Euler1d.
    /// State, clock and accepted report slices remain unchanged on refusal or
    /// cancellation. Successful stepping allocates no numerical workspace.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, Error> {
        poll(&mut cancelled)?;
        if !positive(maximum_duration.value())
            || !positive(controls.cfl)
            || controls.cfl > 0.45
            || controls.maximum_halvings > 32
        {
            return Err(Error::Invalid("network time step, CFL or retry controls"));
        }
        let mut rate = 0.0_f64;
        let mut before = Accumulator::default();
        for (index, branch) in self.branches.iter_mut().enumerate() {
            for (i, &state) in branch.states.iter().enumerate() {
                let p = branch.gas.decode(
                    state,
                    Some(branch.work.primitive[i].temperature),
                    &mut cancelled,
                )?;
                branch.work.primitive[i] = p;
                let area =
                    0.5 * branch.mesh.face_areas()[i] + 0.5 * branch.mesh.face_areas()[i + 1];
                let local = (p.velocity.abs() + p.sound_speed) * (area / branch.mesh.volumes()[i]);
                if !positive(local) {
                    return Err(Error::Invalid("unrepresentable initial network CFL"));
                }
                rate = rate.max(local);
            }
            let value = totals(&branch.states, branch.mesh.volumes(), &mut cancelled)?;
            self.work.branch_before[index] = value;
            before.add(value.mass, value.energy)?;
        }
        self.work.plenum_signal_rates.fill(0.0);
        for (i, plenum) in self.plenums.iter().enumerate() {
            self.work.plenum_primitive[i] = self.branches[0].gas.decode(
                plenum.state,
                Some(self.work.plenum_primitive[i].temperature),
                &mut cancelled,
            )?;
            before.add(
                plenum.volume_m3 * plenum.state.density,
                plenum.volume_m3 * plenum.state.energy,
            )?;
        }
        for connection in &self.connections {
            poll(&mut cancelled)?;
            let branch = &self.branches[connection.branch];
            let face = match connection.end {
                BranchEnd::Left => 0,
                BranchEnd::Right => branch.states.len(),
            };
            let volume = self.plenums[connection.plenum].volume_m3;
            let speed = self.work.plenum_primitive[connection.plenum].sound_speed;
            let local = 0.5 * branch.mesh.face_areas()[face] * (speed / volume);
            if !positive(local) {
                return Err(Error::Invalid("unrepresentable initial plenum CFL"));
            }
            self.work.plenum_signal_rates[connection.plenum] += local;
            let sum = self.work.plenum_signal_rates[connection.plenum];
            if !sum.is_finite() {
                return Err(Error::Invalid("unrepresentable aggregate plenum CFL"));
            }
            rate = rate.max(sum);
        }
        let before = before.totals();
        let mut dt = maximum_duration.value().min(0.9 * controls.cfl / rate);
        let mut last = Error::Invalid("no admissible network time increment");
        let mut first_order_retry = false;
        for halvings in 0..=controls.maximum_halvings {
            poll(&mut cancelled)?;
            let new_time = self.time().value() + dt;
            if !positive(dt) || !positive(new_time - self.time().value()) || !new_time.is_finite() {
                return Err(Error::Invalid(
                    "network time step cannot advance the represented clock",
                ));
            }
            let reconstruction = if first_order_retry {
                Reconstruction::FirstOrder
            } else {
                controls.reconstruction
            };
            match self.attempt(dt, controls.cfl, reconstruction, &mut cancelled) {
                Ok(maximum_rate) => {
                    let report = self.prepare_report(
                        before,
                        dt,
                        maximum_rate,
                        halvings,
                        first_order_retry && controls.reconstruction == Reconstruction::Muscl,
                        &mut cancelled,
                    )?;
                    // Everything after this final cancellation point is infallible.
                    // No accepted branch, compact state, clock or report is changed
                    // until the complete network candidate has been admitted.
                    poll(&mut cancelled)?;
                    for branch in &mut self.branches {
                        core::mem::swap(&mut branch.states, &mut branch.work.candidate);
                        branch.time = new_time;
                    }
                    for (plenum, &candidate) in
                        self.plenums.iter_mut().zip(&self.work.plenum_candidate)
                    {
                        plenum.state = candidate;
                    }
                    core::mem::swap(
                        &mut self.port_exchanges,
                        &mut self.work.port_exchange_candidate,
                    );
                    core::mem::swap(
                        &mut self.branch_reports,
                        &mut self.work.branch_report_candidate,
                    );
                    self.last_report = Some(report);
                    return Ok(report);
                }
                Err(Error::Cancelled) => return Err(Error::Cancelled),
                Err(error) => {
                    if !matches!(error, Error::CflExceeded { .. }) {
                        first_order_retry = true;
                    }
                    last = error;
                    dt *= 0.5;
                }
            }
        }
        Err(Error::StepRetries {
            attempts: controls.maximum_halvings + 1,
            last: Box::new(last),
        })
    }

    fn attempt(
        &mut self,
        dt: f64,
        cfl: f64,
        reconstruction: Reconstruction,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let a = self.evaluate(false, reconstruction, cancelled)?;
        admit_cfl(a, dt, cfl)?;
        self.forward_stage(false, dt, cancelled)?;
        let b = self.evaluate(true, reconstruction, cancelled)?;
        admit_cfl(b, dt, cfl)?;
        self.forward_stage(true, dt, cancelled)?;
        Ok(a.max(b))
    }

    fn evaluate(
        &mut self,
        second: bool,
        reconstruction: Reconstruction,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        self.work.endpoint_fluxes.fill([None; 2]);
        self.work.plenum_rhs.fill([0.0; 2]);
        self.work.plenum_compensation.fill([0.0; 2]);
        self.work.plenum_signal_rates.fill(0.0);
        for i in 0..self.plenums.len() {
            let q = if second {
                self.work.plenum_stage[i]
            } else {
                self.plenums[i].state
            };
            self.work.plenum_primitive[i] = self.branches[0].gas.decode(
                q,
                Some(self.work.plenum_primitive[i].temperature),
                cancelled,
            )?;
        }
        for connection in &self.connections {
            poll(cancelled)?;
            let branch = &mut self.branches[connection.branch];
            let (cell, face) = match connection.end {
                BranchEnd::Left => (0, 0),
                BranchEnd::Right => (branch.states.len() - 1, branch.states.len()),
            };
            let q = if second {
                branch.work.stage[cell]
            } else {
                branch.states[cell]
            };
            let primitive =
                branch
                    .gas
                    .decode(q, Some(branch.work.primitive[cell].temperature), cancelled)?;
            branch.work.primitive[cell] = primitive;
            let plenum_state = if second {
                self.work.plenum_stage[connection.plenum]
            } else {
                self.plenums[connection.plenum].state
            };
            let flux = branch.gas.plenum_flux(
                q,
                primitive,
                plenum_state,
                self.work.plenum_primitive[connection.plenum],
                connection.end.normal(),
                cancelled,
            )?;
            self.work.endpoint_fluxes[connection.branch][connection.end.index()] = Some(flux);
            let area = branch.mesh.face_areas()[face];
            let normal_area = connection.end.normal() * area;
            let value = Totals {
                mass: normal_area * flux.mass,
                momentum: normal_area * flux.momentum,
                energy: normal_area * flux.energy,
            };
            if !value.array().iter().all(|x| x.is_finite()) {
                return Err(Error::Invalid("unrepresentable network port exchange"));
            }
            for (k, entry) in [value.mass, value.energy].into_iter().enumerate() {
                let y = entry - self.work.plenum_compensation[connection.plenum][k];
                let next = self.work.plenum_rhs[connection.plenum][k] + y;
                if !next.is_finite() {
                    return Err(Error::Invalid("unrepresentable plenum port sum"));
                }
                self.work.plenum_compensation[connection.plenum][k] =
                    (next - self.work.plenum_rhs[connection.plenum][k]) - y;
                self.work.plenum_rhs[connection.plenum][k] = next;
            }
            let rate = 0.5 * area * (flux.signal_speed / self.plenums[connection.plenum].volume_m3);
            if !positive(rate) {
                return Err(Error::Invalid("unrepresentable plenum face CFL"));
            }
            self.work.plenum_signal_rates[connection.plenum] += rate;
            if !self.work.plenum_signal_rates[connection.plenum].is_finite() {
                return Err(Error::Invalid("unrepresentable plenum face-speed sum"));
            }
        }
        let mut maximum_rate = self
            .work
            .plenum_signal_rates
            .iter()
            .copied()
            .fold(0.0, f64::max);
        for (index, branch) in self.branches.iter_mut().enumerate() {
            let current = if second {
                branch.work.stage.as_slice()
            } else {
                branch.states.as_slice()
            };
            let stage = operator(
                &branch.gas,
                current,
                &branch.mesh,
                branch.boundaries,
                reconstruction,
                None,
                &mut branch.work.primitive,
                &mut branch.work.slopes,
                &mut branch.work.faces,
                self.work.endpoint_fluxes[index],
                cancelled,
            )?;
            maximum_rate = maximum_rate.max(stage.signal_rate);
            if second {
                self.work.stages_b[index] = stage;
            } else {
                self.work.stages_a[index] = stage;
            }
        }
        if !positive(maximum_rate) {
            return Err(Error::Invalid("unrepresentable network stage CFL"));
        }
        Ok(maximum_rate)
    }

    fn forward_stage(
        &mut self,
        second: bool,
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        for branch in &mut self.branches {
            for i in 0..branch.states.len() {
                poll(cancelled)?;
                let q = if second {
                    branch.work.stage[i]
                } else {
                    branch.states[i]
                };
                let advanced = forward(
                    q,
                    i,
                    dt,
                    &branch.mesh,
                    &branch.work.primitive,
                    &branch.work.faces,
                    None,
                );
                if !branch.gas.admissible_state(advanced) {
                    return Err(Error::Invalid(
                        "network branch RK stage outside source domain",
                    ));
                }
                if second {
                    let candidate = Conserved::from_array(core::array::from_fn(|k| {
                        0.5 * branch.states[i].array()[k] + 0.5 * advanced.array()[k]
                    }));
                    branch.gas.decode(
                        candidate,
                        Some(branch.work.primitive[i].temperature),
                        cancelled,
                    )?;
                    branch.work.candidate[i] = candidate;
                } else {
                    branch.work.stage[i] = advanced;
                }
            }
        }
        for (i, plenum) in self.plenums.iter().enumerate() {
            poll(cancelled)?;
            let q = if second {
                self.work.plenum_stage[i]
            } else {
                plenum.state
            };
            let ratio = dt / plenum.volume_m3;
            let advanced = Conserved {
                density: q.density + ratio * self.work.plenum_rhs[i][0],
                momentum: 0.0,
                energy: q.energy + ratio * self.work.plenum_rhs[i][1],
            };
            if !self.branches[0].gas.admissible_state(advanced) {
                return Err(Error::Invalid(
                    "network plenum RK stage outside source domain",
                ));
            }
            if second {
                let candidate = Conserved {
                    density: 0.5 * plenum.state.density + 0.5 * advanced.density,
                    momentum: 0.0,
                    energy: 0.5 * plenum.state.energy + 0.5 * advanced.energy,
                };
                self.branches[0].gas.decode(
                    candidate,
                    Some(self.work.plenum_primitive[i].temperature),
                    cancelled,
                )?;
                self.work.plenum_candidate[i] = candidate;
            } else {
                self.work.plenum_stage[i] = advanced;
            }
        }
        Ok(())
    }

    fn prepare_report(
        &mut self,
        before: NetworkTotals,
        dt: f64,
        maximum_rate: f64,
        halvings: usize,
        first_order_retry: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, Error> {
        let factor = 0.5 * dt;
        let mut after = Accumulator::default();
        let mut boundary = Accumulator::default();
        let (mut fallbacks, mut limited) = (0, 0);
        for (index, branch) in self.branches.iter().enumerate() {
            poll(cancelled)?;
            let a = self.work.stages_a[index];
            let b = self.work.stages_b[index];
            let left_area = branch.mesh.face_areas()[0];
            let right_area = branch.mesh.face_areas()[branch.states.len()];
            let end_change = Totals::from_array(core::array::from_fn(|k| {
                factor
                    * (left_area * (a.left.array()[k] + b.left.array()[k])
                        - right_area * (a.right.array()[k] + b.right.array()[k]))
            }));
            let left_boundary_change = Totals::from_array(core::array::from_fn(|k| {
                factor * left_area * (a.left.array()[k] + b.left.array()[k])
            }));
            let right_boundary_change = Totals::from_array(core::array::from_fn(|k| {
                -factor * right_area * (a.right.array()[k] + b.right.array()[k])
            }));
            let branch_after = totals(&branch.work.candidate, branch.mesh.volumes(), cancelled)?;
            let branch_before = self.work.branch_before[index];
            let impulse = factor * (a.area_force + b.area_force);
            let residual = Totals::from_array(core::array::from_fn(|k| {
                branch_after.array()[k]
                    - branch_before.array()[k]
                    - end_change.array()[k]
                    - if k == 1 { impulse } else { 0.0 }
            }));
            if !branch_after
                .array()
                .iter()
                .chain(&end_change.array())
                .chain(&left_boundary_change.array())
                .chain(&right_boundary_change.array())
                .chain(&residual.array())
                .all(|x| x.is_finite())
            {
                return Err(Error::Invalid("non-finite branch balance in network"));
            }
            self.work.branch_report_candidate[index] = StepReport {
                elapsed: Time::new(dt),
                maximum_cfl: dt * a.signal_rate.max(b.signal_rate),
                halvings,
                hll_fallback_faces: a.fallbacks + b.fallbacks,
                limited_cells: a.limited + b.limited,
                first_order_retry,
                before: branch_before,
                after: branch_after,
                boundary_change: end_change,
                left_boundary_change,
                right_boundary_change,
                area_momentum_impulse: impulse,
                gravity_momentum_impulse: 0.0,
                gravity_work: 0.0,
                conservation_residual: residual,
            };
            after.add(branch_after.mass, branch_after.energy)?;
            // Periodic ends are one internal face, including in a network with
            // otherwise disconnected branches. Do not accumulate their roundoff
            // as two independent external boundary exchanges.
            if branch.boundaries != Boundaries::Periodic {
                if self.endpoint_connections[index][0].is_none() {
                    boundary.add(
                        factor * left_area * (a.left.mass + b.left.mass),
                        factor * left_area * (a.left.energy + b.left.energy),
                    )?;
                }
                if self.endpoint_connections[index][1].is_none() {
                    boundary.add(
                        -factor * right_area * (a.right.mass + b.right.mass),
                        -factor * right_area * (a.right.energy + b.right.energy),
                    )?;
                }
            }
            fallbacks += a.fallbacks + b.fallbacks;
            limited += a.limited + b.limited;
        }
        for (plenum, q) in self.plenums.iter().zip(&self.work.plenum_candidate) {
            poll(cancelled)?;
            after.add(plenum.volume_m3 * q.density, plenum.volume_m3 * q.energy)?;
        }
        for (port, connection) in self.connections.iter().enumerate() {
            poll(cancelled)?;
            // Retain the exact opposite of the branch's accepted end exchange.
            // This also avoids storing another copy of both stages' port fluxes.
            let branch = self.work.branch_report_candidate[connection.branch];
            let into_branch = match connection.end {
                BranchEnd::Left => branch.left_boundary_change,
                BranchEnd::Right => branch.right_boundary_change,
            };
            let exchange = PortExchange {
                mass_to_plenum: -into_branch.mass,
                momentum_to_plenum: -into_branch.momentum,
                energy_to_plenum: -into_branch.energy,
            };
            if ![
                exchange.mass_to_plenum,
                exchange.momentum_to_plenum,
                exchange.energy_to_plenum,
            ]
            .iter()
            .all(|x| x.is_finite())
            {
                return Err(Error::Invalid("non-finite accepted port balance"));
            }
            self.work.port_exchange_candidate[port] = exchange;
        }
        let (after, boundary_change) = (after.totals(), boundary.totals());
        let conservation_residual = NetworkTotals {
            mass: after.mass - before.mass - boundary_change.mass,
            energy: after.energy - before.energy - boundary_change.energy,
        };
        if !conservation_residual.mass.is_finite() || !conservation_residual.energy.is_finite() {
            return Err(Error::Invalid("non-finite complete network balance"));
        }
        Ok(NetworkStepReport {
            elapsed: Time::new(dt),
            maximum_cfl: maximum_rate * dt,
            halvings,
            hll_fallback_faces: fallbacks,
            limited_cells: limited,
            first_order_retry,
            before,
            after,
            boundary_change,
            conservation_residual,
        })
    }
}

fn admit_cfl(rate: f64, dt: f64, maximum: f64) -> Result<(), Error> {
    let observed = rate * dt;
    if !observed.is_finite() || observed > maximum * (1.0 + 64.0 * f64::EPSILON) {
        Err(Error::CflExceeded { observed, maximum })
    } else {
        Ok(())
    }
}
