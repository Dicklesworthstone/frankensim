use super::super::evolve::{forward, operator, totals};
use super::super::{
    Boundaries, Conserved, Error, Reconstruction, StepControl, StepReport, Totals, poll, positive,
};
use super::{
    Accumulator, BranchEnd, GasNetwork1d, NetworkGravityReport, NetworkStepReport,
    NetworkThermalReport, NetworkTotals, NetworkWallOperator, PlenumBoundaryWork,
    PlenumThermalOperator, PortExchange,
};
use super::thermal::BranchSource;
use crate::multicomponent::network::wall_exchange::NetworkWallReport;
use crate::multicomponent::wall_exchange::WallControl;
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, Error> {
        self.step_inner(maximum_duration, controls, None, None, cancelled)
    }

    /// Advance prescribed chamber motion and port transport on the same SSP RK2
    /// stages. Rates are signed dV/dt [m³/s], one per retained plenum, and are
    /// constant over the accepted interval. Negative rates compress the gas.
    ///
    /// Mass and internal energy advance extensively on the represented volumes.
    /// The original variable-Cp source supplies pressure and compression heating;
    /// `-p dV/dt` is mechanical work, not an added heat source. Accepted volume
    /// endpoints and work are returned by [`Self::accepted_plenum_boundary_work`].
    /// Acoustic and volume-change CFL limits may shorten the requested duration.
    /// Failure publishes no state, geometry, clock or receipt. All-zero rates
    /// use the original rigid arithmetic. Compact chambers remain perfectly
    /// mixed: this does not resolve their internal piston shocks or inertia.
    pub fn step_with_volume_rates(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        volume_rates_m3_s: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, Error> {
        let moving = self.admit_volume_rates(volume_rates_m3_s, &mut cancelled)?;
        self.step_inner(
            maximum_duration,
            controls,
            moving.then_some(volume_rates_m3_s),
            None,
            cancelled,
        )
    }

    /// Advance transport, then original wall sources on plenums and any
    /// explicitly bound branch cells. The original plenum-only call is unchanged.
    ///
    /// Heat uses the actual accepted wave duration after all CFL retries. The
    /// complete gas and wall candidates publish together after the final poll.
    /// Splitting is first order even though the inner wave step is SSP RK2;
    /// original backward-Euler thermal work and source-domain caps still apply.
    pub fn step_with_thermal(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        thermal: &mut PlenumThermalOperator<'_>,
        thermal_controls: WallControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, Error> {
        self.step_inner(
            maximum_duration,
            controls,
            None,
            Some((thermal, thermal_controls)),
            cancelled,
        )
    }

    /// Advance the same transactional wave/wall step with distributed branch
    /// heat and momentum loss. Stationary-wall drag retains total gas energy;
    /// its lost kinetic energy is already available to the thermal solve.
    /// Different branch impulses remain in their original axial coordinates.
    pub fn step_with_walls(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        walls: &mut NetworkWallOperator<'_>,
        wall_controls: WallControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, Error> {
        self.step_inner(
            maximum_duration,
            controls,
            None,
            Some((walls, wall_controls)),
            cancelled,
        )
    }

    /// Compose chamber compression/expansion with the original finite wall
    /// inventories, shared walls and branch losses in one atomic interval.
    ///
    /// Transport and pressure work use SSP RK2; the existing implicit wall
    /// source follows at the actual accepted endpoint volume and duration.
    /// This composition remains first order in thermal/transport splitting.
    /// Accepted wall geometry is rebound only after the whole step succeeds.
    /// Authored conductances, wall masses and port areas remain explicit inputs:
    /// a volume rate alone does not determine their geometric evolution.
    pub fn step_with_walls_and_volume_rates(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        volume_rates_m3_s: &[f64],
        walls: &mut NetworkWallOperator<'_>,
        wall_controls: WallControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, Error> {
        let moving = self.admit_volume_rates(volume_rates_m3_s, &mut cancelled)?;
        self.step_inner(
            maximum_duration,
            controls,
            moving.then_some(volume_rates_m3_s),
            Some((walls, wall_controls)),
            cancelled,
        )
    }

    fn admit_volume_rates(
        &self,
        volume_rates_m3_s: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, Error> {
        poll(cancelled)?;
        if volume_rates_m3_s.len() != self.plenums.len() {
            return Err(Error::Invalid("one volume rate is required for every plenum"));
        }
        let mut moving = false;
        for &rate in volume_rates_m3_s {
            poll(cancelled)?;
            if !rate.is_finite() {
                return Err(Error::Invalid("plenum volume rates must be finite"));
            }
            moving |= rate != 0.0;
        }
        Ok(moving)
    }

    fn step_inner(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        volume_rates: Option<&[f64]>,
        mut thermal: Option<(&mut PlenumThermalOperator<'_>, WallControl)>,
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
        if let Some((operator, controls)) = &thermal {
            controls.validate()?;
            operator.validate_binding(self)?;
        }
        let mut rate = 0.0_f64;
        let mut before = Accumulator::default();
        for (index, branch) in self.branches.iter_mut().enumerate() {
            for (i, &state) in branch.states.iter().enumerate() {
                let p = if let Some(balance) = branch
                    .hydrostatic
                    .as_ref()
                    .filter(|b| b.cells[i].0 == state)
                {
                    poll(&mut cancelled)?;
                    balance.cells[i].1
                } else {
                    branch.gas.decode(
                        state,
                        Some(branch.work.primitive[i].temperature),
                        &mut cancelled,
                    )?
                };
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
            if let Some(rates) = volume_rates {
                let local = rates[i].abs() / plenum.volume_m3;
                if !local.is_finite() {
                    return Err(Error::Invalid("unrepresentable initial plenum volume rate"));
                }
                rate = rate.max(local);
            }
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
        let potential_before = self
            .centered_potential_storage(false, &mut cancelled)?
            .energy;
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
            let trial = self
                .attempt(dt, controls.cfl, reconstruction, volume_rates, &mut cancelled)
                .and_then(|maximum_rate| {
                    let heat = match &mut thermal {
                        Some((operator, controls)) => Some(operator.prepare(
                            &mut self.branches,
                            &mut self.work.plenum_candidate,
                            &self.work.plenum_work_candidate,
                            Time::new(dt),
                            *controls,
                            &mut cancelled,
                        )?),
                        None => None,
                    };
                    Ok((maximum_rate, heat))
                });
            match trial {
                Ok((maximum_rate, heat)) => {
                    let branch_sources = thermal
                        .as_ref()
                        .map(|(operator, _)| operator.branch_source_candidates());
                    let report = self.prepare_report(
                        (before, potential_before),
                        dt,
                        maximum_rate,
                        halvings,
                        first_order_retry && controls.reconstruction == Reconstruction::Muscl,
                        heat,
                        branch_sources,
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
                    for (i, plenum) in self.plenums.iter_mut().enumerate() {
                        plenum.state = self.work.plenum_candidate[i];
                        plenum.volume_m3 = self.work.plenum_work_candidate[i].final_volume_m3;
                    }
                    core::mem::swap(
                        &mut self.port_exchanges,
                        &mut self.work.port_exchange_candidate,
                    );
                    core::mem::swap(
                        &mut self.branch_reports,
                        &mut self.work.branch_report_candidate,
                    );
                    if let (Some((operator, _)), Some(heat)) = (&mut thermal, heat) {
                        operator.commit(heat, &self.work.plenum_work_candidate);
                    }
                    core::mem::swap(&mut self.plenum_work, &mut self.work.plenum_work_candidate);
                    self.last_report = Some(report);
                    return Ok(report);
                }
                Err(Error::Cancelled) => return Err(Error::Cancelled),
                Err(error) => {
                    if !matches!(error, Error::CflExceeded { .. } | Error::Thermal(_)) {
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
        volume_rates: Option<&[f64]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        // Every retry starts from accepted geometry. A nonzero rate must move
        // the represented volume rather than spend work on a rounded-away step.
        for (i, plenum) in self.plenums.iter().enumerate() {
            poll(cancelled)?;
            let rate = volume_rates.map_or(0.0, |rates| rates[i]);
            let volume = plenum.volume_m3 + dt * rate;
            if !positive(volume)
                || (rate != 0.0 && volume == plenum.volume_m3)
                || (rate != 0.0
                    && (volume - plenum.volume_m3).is_sign_positive() != rate.is_sign_positive())
            {
                return Err(Error::Invalid(
                    "plenum motion has no positive representable endpoint volume",
                ));
            }
            self.work.plenum_work_candidate[i] = PlenumBoundaryWork {
                initial_volume_m3: plenum.volume_m3,
                final_volume_m3: volume,
                work_on_gas_j: 0.0,
            };
        }
        let a = self.evaluate(false, reconstruction, volume_rates, cancelled)?;
        admit_cfl(a, dt, cfl)?;
        self.forward_stage(false, dt, volume_rates, cancelled)?;
        let b = self.evaluate(true, reconstruction, volume_rates, cancelled)?;
        admit_cfl(b, dt, cfl)?;
        self.forward_stage(true, dt, volume_rates, cancelled)?;
        Ok(a.max(b))
    }

    fn evaluate(
        &mut self,
        second: bool,
        reconstruction: Reconstruction,
        volume_rates: Option<&[f64]>,
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
            if let Some(rates) = volume_rates {
                let volume = if second {
                    self.work.plenum_work_candidate[i].final_volume_m3
                } else {
                    self.plenums[i].volume_m3
                };
                let local = rates[i].abs() / volume;
                if !local.is_finite() {
                    return Err(Error::Invalid("unrepresentable plenum stage volume rate"));
                }
                self.work.plenum_signal_rates[i] = local;
            }
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
            let primitive = if let Some(balance) =
                branch.hydrostatic.as_ref().filter(|b| b.cells[cell].0 == q)
            {
                balance.cells[cell].1
            } else {
                branch
                    .gas
                    .decode(q, Some(branch.work.primitive[cell].temperature), cancelled)?
            };
            branch.work.primitive[cell] = primitive;
            let (q, primitive) = if let Some(balance) = &branch.hydrostatic {
                balance.end_trace(&branch.gas, q, primitive, cell, face)?
            } else {
                (q, primitive)
            };
            let plenum_state = if second {
                self.work.plenum_stage[connection.plenum]
            } else {
                self.plenums[connection.plenum].state
            };
            // Identical admitted conservative states have the same original
            // thermodynamics. Reuse the face's immutable evaluation instead of
            // introducing an inversion-roundoff pressure difference at rest.
            let plenum_primitive = if q == plenum_state {
                primitive
            } else {
                self.work.plenum_primitive[connection.plenum]
            };
            let flux = branch.gas.plenum_flux(
                q,
                primitive,
                plenum_state,
                plenum_primitive,
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
            let volume = if second {
                self.work.plenum_work_candidate[connection.plenum].final_volume_m3
            } else {
                self.plenums[connection.plenum].volume_m3
            };
            let rate = 0.5 * area * (flux.signal_speed / volume);
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
                branch.hydrostatic.as_ref(),
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
        volume_rates: Option<&[f64]>,
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
                    branch.hydrostatic.as_ref(),
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
                    if !branch
                        .hydrostatic
                        .as_ref()
                        .is_some_and(|b| b.cells[i].0 == candidate)
                    {
                        branch.gas.decode(
                            candidate,
                            Some(branch.work.primitive[i].temperature),
                            cancelled,
                        )?;
                    }
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
            let volume_rate = volume_rates.map_or(0.0, |rates| rates[i]);
            if volume_rate != 0.0 {
                // RK averages extensive storage, never densities on different
                // volumes. Both stages use the same source-owned pressure work
                // recorded below; formation references cancel at fixed mass.
                let initial_volume = plenum.volume_m3;
                let final_volume = self.work.plenum_work_candidate[i].final_volume_m3;
                let stage_volume = if second { final_volume } else { initial_volume };
                let pressure_power = -self.work.plenum_primitive[i].pressure * volume_rate;
                if !pressure_power.is_finite() {
                    return Err(Error::Invalid("unrepresentable mechanical plenum power"));
                }
                let mass = stage_volume * q.density + dt * self.work.plenum_rhs[i][0];
                let energy = stage_volume * q.energy
                    + dt * (self.work.plenum_rhs[i][1] + pressure_power);
                let next = Conserved {
                    density: if second {
                        (0.5 * initial_volume * plenum.state.density + 0.5 * mass) / final_volume
                    } else {
                        mass / final_volume
                    },
                    momentum: 0.0,
                    energy: if second {
                        (0.5 * initial_volume * plenum.state.energy + 0.5 * energy) / final_volume
                    } else {
                        energy / final_volume
                    },
                };
                self.branches[0].gas.decode(
                    next,
                    Some(self.work.plenum_primitive[i].temperature),
                    cancelled,
                )?;
                let work = self.work.plenum_work_candidate[i].work_on_gas_j
                    + (0.5 * dt) * pressure_power;
                if !work.is_finite() {
                    return Err(Error::Invalid("unrepresentable mechanical plenum work"));
                }
                self.work.plenum_work_candidate[i].work_on_gas_j = work;
                if second {
                    self.work.plenum_candidate[i] = next;
                } else {
                    self.work.plenum_stage[i] = next;
                }
                continue;
            }
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
        before: (NetworkTotals, f64),
        dt: f64,
        maximum_rate: f64,
        halvings: usize,
        first_order_retry: bool,
        heat: Option<NetworkWallReport>,
        branch_sources: Option<&[BranchSource]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, Error> {
        let factor = 0.5 * dt;
        let (before, potential_before) = before;
        let mut after = Accumulator::default();
        let mut boundary = Accumulator::default();
        let mut gravity_work_sum = Accumulator::default();
        let mut potential_boundary = Accumulator::default();
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
            let gravity_impulse = factor * (a.gravity_force + b.gravity_force);
            let gravity_work = factor * (a.gravity_power + b.gravity_power);
            gravity_work_sum.add(0.0, gravity_work)?;
            let wall = branch_sources.map_or(BranchSource::default(), |sources| sources[index]);
            let residual = Totals::from_array(core::array::from_fn(|k| {
                branch_after.array()[k]
                    - branch_before.array()[k]
                    - end_change.array()[k]
                    - match k {
                        1 => impulse + gravity_impulse + wall.momentum,
                        2 => gravity_work + wall.heat + wall.work,
                        _ => 0.0,
                    }
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
                gravity_momentum_impulse: gravity_impulse,
                gravity_work,
                wall_momentum_impulse: wall.momentum,
                wall_heat: wall.heat,
                wall_work: wall.work,
                source_change: crate::body_force::ForceChange::default(),
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
                    if let Some(gravity) = &self.gravity {
                        potential_boundary.add(
                            0.0,
                            left_boundary_change.mass
                                * gravity.centered_branch_potential(
                                    index,
                                    branch,
                                    branch.mesh.positions()[0],
                                ),
                        )?;
                    }
                }
                if self.endpoint_connections[index][1].is_none() {
                    boundary.add(
                        -factor * right_area * (a.right.mass + b.right.mass),
                        -factor * right_area * (a.right.energy + b.right.energy),
                    )?;
                    if let Some(gravity) = &self.gravity {
                        potential_boundary.add(
                            0.0,
                            right_boundary_change.mass
                                * gravity.centered_branch_potential(
                                    index,
                                    branch,
                                    branch.mesh.positions()[branch.states.len()],
                                ),
                        )?;
                    }
                }
            }
            fallbacks += a.fallbacks + b.fallbacks;
            limited += a.limited + b.limited;
        }
        for (i, q) in self.work.plenum_candidate.iter().enumerate() {
            poll(cancelled)?;
            let volume = self.work.plenum_work_candidate[i].final_volume_m3;
            after.add(volume * q.density, volume * q.energy)?;
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
        let gravity_work = gravity_work_sum.totals().energy;
        let heat_to_gas = heat.map_or(0.0, |report| report.heat_to_gas_j);
        let wall_work = heat.map_or(0.0, |report| report.mechanical_work_to_gas_j);
        let mut boundary_work = Accumulator::default();
        for motion in &self.work.plenum_work_candidate {
            poll(cancelled)?;
            boundary_work.add(0.0, motion.work_on_gas_j)?;
        }
        let boundary_work_on_gas_j = boundary_work.totals().energy;
        let conservation_residual = NetworkTotals {
            mass: after.mass - before.mass - boundary_change.mass,
            energy: after.energy
                - before.energy
                - boundary_change.energy
                - gravity_work
                - heat_to_gas
                - wall_work
                - boundary_work_on_gas_j,
        };
        if !conservation_residual.mass.is_finite() || !conservation_residual.energy.is_finite() {
            return Err(Error::Invalid("non-finite complete network balance"));
        }
        let gravity = if let Some(frame) = &self.gravity {
            let after_centered_potential_energy =
                self.centered_potential_storage(true, cancelled)?.energy;
            let boundary_centered_potential_energy_change = potential_boundary.totals().energy;
            // Every term was accumulated AFTER subtracting the same potential
            // reference. Merely differencing absolute potential storage would
            // discard small resolved gravity work under a large common datum.
            let total_residual =
                (after.energy - before.energy - boundary_change.energy - heat_to_gas - wall_work
                    - boundary_work_on_gas_j)
                    + (after_centered_potential_energy
                        - potential_before
                        - boundary_centered_potential_energy_change);
            let absolute_conservation_residual =
                total_residual + frame.reference() * conservation_residual.mass;
            if !total_residual.is_finite()
                || !absolute_conservation_residual.is_finite()
                || !(after_centered_potential_energy + frame.reference() * after.mass).is_finite()
            {
                return Err(Error::Invalid("non-finite gas-plus-potential balance"));
            }
            Some(NetworkGravityReport {
                reference_potential_j_kg: frame.reference(),
                before_centered_potential_energy: potential_before,
                after_centered_potential_energy,
                boundary_centered_potential_energy_change,
                gravity_work,
                conservation_residual: total_residual,
                absolute_conservation_residual,
            })
        } else {
            None
        };
        let thermal = heat.map(|source| NetworkThermalReport {
            source,
            gas_wall_conservation_residual_j: conservation_residual.energy
                + source.energy_residual_j,
            gas_wall_potential_conservation_residual_j: gravity
                .map(|g| g.conservation_residual + source.energy_residual_j),
            absolute_gas_wall_potential_conservation_residual_j: gravity.map(|g| {
                (g.conservation_residual + source.energy_residual_j)
                    + g.reference_potential_j_kg * conservation_residual.mass
            }),
        });
        if thermal.is_some_and(|report| {
            !report.gas_wall_conservation_residual_j.is_finite()
                || report
                    .gas_wall_potential_conservation_residual_j
                    .is_some_and(|r| !r.is_finite())
                || report
                    .absolute_gas_wall_potential_conservation_residual_j
                    .is_some_and(|r| !r.is_finite())
        }) {
            return Err(Error::Invalid(
                "non-finite complete gas/wall/potential balance",
            ));
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
            boundary_work_on_gas_j,
            conservation_residual,
            gravity,
            thermal,
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
