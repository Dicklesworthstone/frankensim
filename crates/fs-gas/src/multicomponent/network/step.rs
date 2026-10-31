use super::super::evolve::{EndpointFlux, forward, operator, totals};
use super::super::{
    Boundaries, Cell, Error, MAX_COMPONENTS, MAX_SPECIES, Reconstruction, StepControl, StepReport,
    Totals, cell, poll, positive, sum,
};
use super::{
    Accumulator, BranchEnd, MixtureNetwork1d, NetworkStepReport, NetworkTotals, PlenumBoundaryWork,
    PortExchange,
};
use fs_qty::Time;

impl MixtureNetwork1d {
    /// Combined species masses and energy in deterministic branch-then-plenum order.
    pub fn totals(&self) -> NetworkTotals {
        self.current_totals(&mut || false)
            .expect("admitted mixture network storage and an infallible callback cannot refuse")
    }
    fn current_totals(&self, cancelled: &mut impl FnMut() -> bool) -> Result<NetworkTotals, Error> {
        let n = self.gas().species_count();
        let mut sum = Accumulator::new(n);
        for branch in &self.branches {
            let value = totals(&branch.states, branch.mesh.volumes(), n, cancelled)?;
            sum.add(&value.species_masses[..n], value.energy)?;
        }
        for plenum in &self.plenums {
            poll(cancelled)?;
            sum.add_cell(plenum.state(), plenum.volume_m3)?;
        }
        Ok(sum.totals())
    }

    /// Advance every branch and plenum by one globally accepted SSP RK2 step.
    ///
    /// Branch and finite-plenum characteristic speeds limit both stages. All
    /// retries restart from original accepted fields. Cancellation or refusal
    /// publishes no state, clock, or accepted report. Numerical buffers are reused.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, Error> {
        self.step_impl(maximum_duration, controls, None, None, cancelled)
    }

    /// Advance prescribed compact-volume motion and branch flow simultaneously.
    ///
    /// `volume_rates_m3_s` contains one finite signed dV/dt in retained plenum
    /// order: negative compresses, positive expands. Rates stay constant over
    /// the actual accepted duration, which may be shortened by acoustic or
    /// volume-rate limits. Each SSP RK2 stage advances extensive species masses
    /// and referenced internal energy, with the mechanical term `-p*dV/dt`
    /// evaluated through the original local composition and NASA-9 calorics.
    /// Shared port enthalpy fluxes remain equal and opposite at both stages.
    ///
    /// All-zero rates use the original rigid-volume arithmetic. Successful
    /// steps allocate no new workspace. Refusal/cancellation retains every
    /// state, volume, clock and accepted report. The per-plenum accepted volumes
    /// and pressure work are available through `accepted_plenum_boundary_work`.
    ///
    /// This is a finite-time perfectly mixed chamber model coupled to resolved
    /// branch waves. Chamber inertia, spatial piston shocks and wall motion of
    /// the branches are not resolved; the caller supplies the mechanical motion.
    pub fn step_with_volume_rates(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        volume_rates_m3_s: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, Error> {
        let moving = self.admit_volume_rates(volume_rates_m3_s, &mut cancelled)?;
        self.step_impl(
            maximum_duration,
            controls,
            moving.then_some(volume_rates_m3_s),
            None,
            cancelled,
        )
    }

    /// Advance the complete network against prescribed species-bearing reservoirs.
    ///
    /// Supply exactly one `[left, right]` pair per branch in retained order.
    /// `Some` overrides an unconnected physical end for this call; `None` keeps
    /// its original boundary or finite-plenum connection. An all-`None` slice
    /// executes the original step exactly. Periodic ends and ends already owned
    /// by a finite plenum cannot receive an external reservoir.
    ///
    /// Each stationary borrowed cell must use this network's original species
    /// bank, canonical order and per-species energy references, as required by
    /// [`super::super::MixtureGas::reservoir_flux`]. Raw cells carry no source
    /// identity tag; that compatibility is the caller's responsibility. Species
    /// dimensions, stationarity and the complete source state are admitted here.
    ///
    /// Pressure and temperature are stagnation conditions for inflow; outflow
    /// sees the reservoir backpressure. The original characteristic closure
    /// transports reservoir species and stagnation enthalpy into a branch and
    /// interior species and total enthalpy out. Every RK stage and halved retry
    /// recomputes each flux from that stage's actual branch end cell. External
    /// exchange enters the network boundary report; finite-plenum exchanges
    /// remain internal, with exact opposite accepted port reports.
    ///
    /// Reservoirs are prescribed throughout the call and are not depleted.
    /// Their states may change between accepted calls. No borrowed state is
    /// retained and no numerical heap workspace is allocated per step. All
    /// branches, finite plenums, clocks and accepted reports publish together;
    /// refusal or cancellation publishes none of them.
    pub fn step_with_reservoirs(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        reservoirs: &[[Option<Cell<'_>>; 2]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, Error> {
        let active = self.admit_reservoirs(reservoirs, &mut cancelled)?;
        self.step_impl(
            maximum_duration,
            controls,
            None,
            active.then_some(reservoirs),
            cancelled,
        )
    }

    /// Advance prescribed end reservoirs and moving finite plenums in one step.
    ///
    /// Source pairs and signed plenum rates have the same admission and retained
    /// order as [`Self::step_with_reservoirs`] and [`Self::step_with_volume_rates`].
    /// Both drives act in each shared SSP RK2 stage: finite-port states use the
    /// evolving plenum volume, species and energy, exterior fluxes use that
    /// stage's branch cells, and `-p*dV/dt` is included in extensive storage.
    /// Reports keep transported external energy separate from mechanical work.
    /// No sequential volume/transport split or intermediate publication occurs.
    /// All-None sources retain the volume-rate path; all-zero rates retain the
    /// rigid reservoir arithmetic. Every retry replays both drives together.
    pub fn step_with_reservoirs_and_volume_rates(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        reservoirs: &[[Option<Cell<'_>>; 2]],
        volume_rates_m3_s: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, Error> {
        let active = self.admit_reservoirs(reservoirs, &mut cancelled)?;
        let moving = self.admit_volume_rates(volume_rates_m3_s, &mut cancelled)?;
        self.step_impl(
            maximum_duration,
            controls,
            moving.then_some(volume_rates_m3_s),
            active.then_some(reservoirs),
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
            return Err(Error::Invalid(
                "one volume rate is required for every plenum",
            ));
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

    fn admit_reservoirs(
        &mut self,
        reservoirs: &[[Option<Cell<'_>>; 2]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, Error> {
        if reservoirs.len() != self.branches.len() {
            poll(cancelled)?;
            return Err(Error::Invalid(
                "one reservoir pair per mixture network branch required",
            ));
        }
        if reservoirs.iter().flatten().all(Option::is_none) {
            return Ok(false);
        }
        poll(cancelled)?;
        for (index, (branch, ends)) in self.branches.iter().zip(reservoirs).enumerate() {
            for (end, reservoir) in ends.iter().copied().enumerate() {
                if let Some(state) = reservoir {
                    poll(cancelled)?;
                    if branch.boundaries == Boundaries::Periodic
                        || self.endpoint_connections[index][end].is_some()
                    {
                        return Err(Error::Invalid(
                            "external mixture reservoir requires an unconnected physical end",
                        ));
                    }
                    if state.momentum != 0.0 {
                        return Err(Error::Invalid(
                            "stationary mixture reservoir has nonzero momentum",
                        ));
                    }
                    self.work.reservoir_primitive[index][end] =
                        branch.gas.primitive(state, &mut *cancelled)?;
                }
            }
        }
        Ok(true)
    }

    fn step_impl(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        volume_rates: Option<&[f64]>,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, Error> {
        poll(&mut cancelled)?;
        if !positive(maximum_duration.value())
            || !positive(controls.cfl)
            || controls.cfl > 0.45
            || controls.maximum_halvings > 32
        {
            return Err(Error::Invalid(
                "mixture network time step, CFL or retry controls",
            ));
        }
        let n = self.gas().species_count();
        let m = n + 2;
        let mut rate = 0.0_f64;
        let mut before = Accumulator::new(n);
        for (index, branch) in self.branches.iter_mut().enumerate() {
            for (i, q) in branch.states.chunks_exact(m).enumerate() {
                let p = branch.gas.decode(
                    q,
                    Some(branch.work.primitive[i].temperature),
                    &mut cancelled,
                )?;
                branch.work.primitive[i] = p;
                let area =
                    0.5 * branch.mesh.face_areas()[i] + 0.5 * branch.mesh.face_areas()[i + 1];
                let local = (p.velocity.abs() + p.sound_speed) * (area / branch.mesh.volumes()[i]);
                if !positive(local) {
                    return Err(Error::Invalid(
                        "unrepresentable initial mixture network CFL",
                    ));
                }
                rate = rate.max(local);
            }
            let value = totals(&branch.states, branch.mesh.volumes(), n, &mut cancelled)?;
            self.work.branch_before[index] = value;
            before.add(&value.species_masses[..n], value.energy)?;
        }
        self.work.plenum_signal_rates.fill(0.0);
        for (i, plenum) in self.plenums.iter().enumerate() {
            self.work.plenum_primitive[i] = self.branches[0].gas.decode_cell(
                plenum.state(),
                Some(self.work.plenum_primitive[i].temperature),
                &mut cancelled,
            )?;
            before.add_cell(plenum.state(), plenum.volume_m3)?;
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
                BranchEnd::Right => branch.cell_count(),
            };
            let local = 0.5
                * branch.mesh.face_areas()[face]
                * (self.work.plenum_primitive[connection.plenum].sound_speed
                    / self.plenums[connection.plenum].volume_m3);
            if !positive(local) {
                return Err(Error::Invalid("unrepresentable initial mixture plenum CFL"));
            }
            self.work.plenum_signal_rates[connection.plenum] += local;
            let total = self.work.plenum_signal_rates[connection.plenum];
            if !total.is_finite() {
                return Err(Error::Invalid(
                    "unrepresentable aggregate mixture plenum CFL",
                ));
            }
            rate = rate.max(total);
        }
        let before = before.totals();
        let mut dt = maximum_duration.value().min(0.9 * controls.cfl / rate);
        let mut last = Error::Invalid("no admissible mixture network time increment");
        let mut first_order_retry = false;
        for halvings in 0..=controls.maximum_halvings {
            poll(&mut cancelled)?;
            let new_time = self.time().value() + dt;
            if !positive(dt) || !positive(new_time - self.time().value()) || !new_time.is_finite() {
                return Err(Error::Invalid(
                    "mixture network time step cannot advance the represented clock",
                ));
            }
            let reconstruction = if first_order_retry {
                Reconstruction::FirstOrder
            } else {
                controls.reconstruction
            };
            match self.attempt(
                dt,
                controls.cfl,
                reconstruction,
                volume_rates,
                reservoirs,
                &mut cancelled,
            ) {
                Ok(maximum_rate) => {
                    let report = self.prepare_report(
                        before,
                        dt,
                        maximum_rate,
                        halvings,
                        first_order_retry && controls.reconstruction == Reconstruction::Muscl,
                        &mut cancelled,
                    )?;
                    // The complete candidate, balances and reports are admitted before
                    // the only publication point. Everything after this poll is infallible.
                    poll(&mut cancelled)?;
                    for branch in &mut self.branches {
                        core::mem::swap(&mut branch.states, &mut branch.work.candidate);
                        branch.time = new_time;
                    }
                    for (i, plenum) in self.plenums.iter_mut().enumerate() {
                        let q = &self.work.plenum_candidate[i * m..(i + 1) * m];
                        plenum.state.partial_densities.copy_from_slice(&q[..n]);
                        plenum.state.momentum = 0.0;
                        plenum.state.energy = q[n + 1];
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
                    core::mem::swap(&mut self.plenum_work, &mut self.work.plenum_work_candidate);
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
        volume_rates: Option<&[f64]>,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        // Admit actual, represented geometry before evaluating either stage.
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
        let a = self.evaluate(false, reconstruction, volume_rates, reservoirs, cancelled)?;
        admit_cfl(a, dt, cfl)?;
        self.forward_stage(false, dt, volume_rates, cancelled)?;
        let b = self.evaluate(true, reconstruction, volume_rates, reservoirs, cancelled)?;
        admit_cfl(b, dt, cfl)?;
        self.forward_stage(true, dt, volume_rates, cancelled)?;
        Ok(a.max(b))
    }

    fn evaluate(
        &mut self,
        second: bool,
        reconstruction: Reconstruction,
        volume_rates: Option<&[f64]>,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let n = self.gas().species_count();
        let m = n + 2;
        self.work.endpoint_data.fill([None; 2]);
        self.work.plenum_rhs.fill(0.0);
        self.work.plenum_compensation.fill(0.0);
        self.work.plenum_signal_rates.fill(0.0);
        for (i, plenum) in self.plenums.iter().enumerate() {
            let q = if second {
                cell(&self.work.plenum_stage[i * m..(i + 1) * m], n)
            } else {
                plenum.state()
            };
            self.work.plenum_primitive[i] = self.branches[0].gas.decode_cell(
                q,
                Some(self.work.plenum_primitive[i].temperature),
                cancelled,
            )?;
            if let Some(rates) = volume_rates {
                let volume = if second {
                    self.work.plenum_work_candidate[i].final_volume_m3
                } else {
                    plenum.volume_m3
                };
                let rate = rates[i].abs() / volume;
                if !rate.is_finite() {
                    return Err(Error::Invalid("unrepresentable plenum stage volume rate"));
                }
                self.work.plenum_signal_rates[i] = rate;
            }
        }
        for connection in &self.connections {
            poll(cancelled)?;
            let branch = &mut self.branches[connection.branch];
            let count = branch.cell_count();
            let (index, face) = match connection.end {
                BranchEnd::Left => (0, 0),
                BranchEnd::Right => (count - 1, count),
            };
            let q = if second {
                &branch.work.stage[index * m..(index + 1) * m]
            } else {
                &branch.states[index * m..(index + 1) * m]
            };
            let primitive =
                branch
                    .gas
                    .decode(q, Some(branch.work.primitive[index].temperature), cancelled)?;
            branch.work.primitive[index] = primitive;
            let p = connection.plenum;
            let plenum_state = if second {
                cell(&self.work.plenum_stage[p * m..(p + 1) * m], n)
            } else {
                self.plenums[p].state()
            };
            let offset = (2 * connection.branch + connection.end.index()) * m;
            let values = &mut self.work.endpoint_values[offset..offset + m];
            let data = branch.gas.plenum_flux(
                cell(q, n),
                primitive,
                plenum_state,
                self.work.plenum_primitive[p],
                connection.end.normal(),
                values,
                cancelled,
            )?;
            self.work.endpoint_data[connection.branch][connection.end.index()] = Some(data);
            let area = branch.mesh.face_areas()[face];
            let normal_area = connection.end.normal() * area;
            // Only species and energy enter compact storage. Momentum is a fixed
            // support impulse, retained later from the branch's accepted end report.
            for k in (0..n).chain(core::iter::once(n + 1)) {
                let value = normal_area * values[k];
                let slot = p * m + k;
                let y = value - self.work.plenum_compensation[slot];
                let next = self.work.plenum_rhs[slot] + y;
                if !value.is_finite() || !next.is_finite() {
                    return Err(Error::Invalid("unrepresentable mixture plenum port sum"));
                }
                self.work.plenum_compensation[slot] = (next - self.work.plenum_rhs[slot]) - y;
                self.work.plenum_rhs[slot] = next;
            }
            let volume = if second {
                self.work.plenum_work_candidate[p].final_volume_m3
            } else {
                self.plenums[p].volume_m3
            };
            let rate = 0.5 * area * (data.signal_speed / volume);
            if !positive(rate) {
                return Err(Error::Invalid("unrepresentable mixture plenum face CFL"));
            }
            self.work.plenum_signal_rates[p] += rate;
            if !self.work.plenum_signal_rates[p].is_finite() {
                return Err(Error::Invalid(
                    "unrepresentable mixture plenum face-speed sum",
                ));
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
            let offset = 2 * index * m;
            if let Some(reservoirs) = reservoirs {
                for (end, reservoir) in reservoirs[index].iter().copied().enumerate() {
                    if let Some(outside) = reservoir {
                        let i = if end == 0 { 0 } else { branch.cell_count() - 1 };
                        let inside = cell(&current[i * m..(i + 1) * m], n);
                        let primitive = branch.gas.decode_cell(
                            inside,
                            Some(branch.work.primitive[i].temperature),
                            cancelled,
                        )?;
                        branch.work.primitive[i] = primitive;
                        let start = offset + end * m;
                        let data = branch.gas.plenum_flux(
                            inside,
                            primitive,
                            outside,
                            self.work.reservoir_primitive[index][end],
                            if end == 0 { -1.0 } else { 1.0 },
                            &mut self.work.endpoint_values[start..start + m],
                            cancelled,
                        )?;
                        self.work.endpoint_data[index][end] = Some(data);
                    }
                }
            }
            let left = self.work.endpoint_data[index][0].map(|data| EndpointFlux {
                values: &self.work.endpoint_values[offset..offset + m],
                data,
            });
            let right = self.work.endpoint_data[index][1].map(|data| EndpointFlux {
                values: &self.work.endpoint_values[offset + m..offset + 2 * m],
                data,
            });
            let stage = operator(
                &branch.gas,
                current,
                &branch.mesh,
                branch.boundaries,
                reconstruction,
                &mut branch.work.primitive,
                &mut branch.work.slopes,
                &mut branch.work.faces,
                &mut branch.work.face_data,
                [left, right],
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
            return Err(Error::Invalid("unrepresentable mixture network stage CFL"));
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
        let n = self.gas().species_count();
        let m = n + 2;
        let mut advanced = [0.0; MAX_COMPONENTS];
        for branch in &mut self.branches {
            for i in 0..branch.cell_count() {
                poll(cancelled)?;
                let q = if second {
                    &branch.work.stage[i * m..(i + 1) * m]
                } else {
                    &branch.states[i * m..(i + 1) * m]
                };
                forward(
                    q,
                    i,
                    dt,
                    &branch.mesh,
                    &branch.work.primitive,
                    &branch.work.faces,
                    &mut advanced[..m],
                    n,
                );
                if !branch.gas.admissible(&advanced[..m]) {
                    return Err(Error::Invalid(
                        "mixture network branch RK stage outside source domain",
                    ));
                }
                if second {
                    let candidate = &mut branch.work.candidate[i * m..(i + 1) * m];
                    for k in 0..m {
                        candidate[k] = 0.5 * branch.states[i * m + k] + 0.5 * advanced[k];
                    }
                    branch.gas.decode(
                        candidate,
                        Some(branch.work.primitive[i].temperature),
                        cancelled,
                    )?;
                } else {
                    branch.work.stage[i * m..(i + 1) * m].copy_from_slice(&advanced[..m]);
                }
            }
        }
        for (i, plenum) in self.plenums.iter().enumerate() {
            poll(cancelled)?;
            let q = if second {
                cell(&self.work.plenum_stage[i * m..(i + 1) * m], n)
            } else {
                plenum.state()
            };
            let volume_rate = volume_rates.map_or(0.0, |rates| rates[i]);
            if volume_rate != 0.0 {
                // Advance the moving control volume's extensive state. A
                // density-only RK average would violate the geometric balance
                // and create/destroy species during volume motion.
                let initial_volume = plenum.volume_m3;
                let final_volume = self.work.plenum_work_candidate[i].final_volume_m3;
                let stage_volume = if second { final_volume } else { initial_volume };
                let pressure_power = -self.work.plenum_primitive[i].pressure * volume_rate;
                if !pressure_power.is_finite() {
                    return Err(Error::Invalid("unrepresentable mechanical plenum power"));
                }
                for k in 0..n {
                    let extensive = stage_volume * q.partial_densities[k]
                        + dt * self.work.plenum_rhs[i * m + k];
                    advanced[k] = if second {
                        (0.5 * initial_volume * plenum.state.partial_densities[k] + 0.5 * extensive)
                            / final_volume
                    } else {
                        extensive / final_volume
                    };
                }
                advanced[n] = 0.0;
                let extensive_energy = stage_volume * q.energy
                    + dt * (self.work.plenum_rhs[i * m + n + 1] + pressure_power);
                advanced[n + 1] = if second {
                    (0.5 * initial_volume * plenum.state.energy + 0.5 * extensive_energy)
                        / final_volume
                } else {
                    extensive_energy / final_volume
                };
                self.branches[0].gas.decode(
                    &advanced[..m],
                    Some(self.work.plenum_primitive[i].temperature),
                    cancelled,
                )?;
                let work =
                    self.work.plenum_work_candidate[i].work_on_gas_j + (0.5 * dt) * pressure_power;
                if !work.is_finite() {
                    return Err(Error::Invalid("unrepresentable mechanical plenum work"));
                }
                self.work.plenum_work_candidate[i].work_on_gas_j = work;
                let destination = if second {
                    &mut self.work.plenum_candidate[i * m..(i + 1) * m]
                } else {
                    &mut self.work.plenum_stage[i * m..(i + 1) * m]
                };
                destination.copy_from_slice(&advanced[..m]);
                continue;
            }
            let ratio = dt / plenum.volume_m3;
            for k in 0..n {
                advanced[k] = q.partial_densities[k] + ratio * self.work.plenum_rhs[i * m + k];
            }
            advanced[n] = 0.0;
            advanced[n + 1] = q.energy + ratio * self.work.plenum_rhs[i * m + n + 1];
            if !self.branches[0].gas.admissible(&advanced[..m]) {
                return Err(Error::Invalid(
                    "mixture network plenum RK stage outside source domain",
                ));
            }
            if second {
                let candidate = &mut self.work.plenum_candidate[i * m..(i + 1) * m];
                for k in 0..n {
                    candidate[k] = 0.5 * plenum.state.partial_densities[k] + 0.5 * advanced[k];
                }
                candidate[n] = 0.0;
                candidate[n + 1] = 0.5 * plenum.state.energy + 0.5 * advanced[n + 1];
                self.branches[0].gas.decode(
                    candidate,
                    Some(self.work.plenum_primitive[i].temperature),
                    cancelled,
                )?;
            } else {
                self.work.plenum_stage[i * m..(i + 1) * m].copy_from_slice(&advanced[..m]);
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
        let n = self.gas().species_count();
        let m = n + 2;
        let factor = 0.5 * dt;
        let mut after = Accumulator::new(n);
        let mut boundary = Accumulator::new(n);
        let (mut fallbacks, mut limited) = (0, 0);
        for (index, branch) in self.branches.iter().enumerate() {
            poll(cancelled)?;
            let a = &self.work.stages_a[index];
            let b = &self.work.stages_b[index];
            let al = branch.mesh.face_areas()[0];
            let ar = branch.mesh.face_areas()[branch.cell_count()];
            let mut end = [0.0; MAX_COMPONENTS];
            let mut left = [0.0; MAX_COMPONENTS];
            let mut right = [0.0; MAX_COMPONENTS];
            let mut residual = [0.0; MAX_COMPONENTS];
            let branch_after = totals(&branch.work.candidate, branch.mesh.volumes(), n, cancelled)?;
            let branch_before = self.work.branch_before[index];
            let impulse = factor * (a.area_force + b.area_force);
            for k in 0..m {
                left[k] = factor * al * (a.left[k] + b.left[k]);
                right[k] = -factor * ar * (a.right[k] + b.right[k]);
                end[k] = factor * (al * (a.left[k] + b.left[k]) - ar * (a.right[k] + b.right[k]));
                residual[k] = branch_after.get(k)
                    - branch_before.get(k)
                    - end[k]
                    - if k == n { impulse } else { 0.0 };
            }
            if !end[..m]
                .iter()
                .chain(&left[..m])
                .chain(&right[..m])
                .chain(&residual[..m])
                .all(|x| x.is_finite())
                || !impulse.is_finite()
            {
                return Err(Error::Invalid(
                    "non-finite mixture branch balance in network",
                ));
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
                boundary_change: Totals::from_values(&end[..m], n),
                left_boundary_change: Totals::from_values(&left[..m], n),
                right_boundary_change: Totals::from_values(&right[..m], n),
                area_momentum_impulse: impulse,
                conservation_residual: Totals::from_values(&residual[..m], n),
            };
            after.add(&branch_after.species_masses[..n], branch_after.energy)?;
            // A periodic end pair represents one internal face, never two
            // separately rounded contributions to the external network balance.
            if branch.boundaries != Boundaries::Periodic {
                if self.endpoint_connections[index][0].is_none() {
                    boundary.add(&left[..n], left[n + 1])?;
                }
                if self.endpoint_connections[index][1].is_none() {
                    boundary.add(&right[..n], right[n + 1])?;
                }
            }
            fallbacks += a.fallbacks + b.fallbacks;
            limited += a.limited + b.limited;
        }
        for i in 0..self.plenums.len() {
            poll(cancelled)?;
            after.add_cell(
                cell(&self.work.plenum_candidate[i * m..(i + 1) * m], n),
                self.work.plenum_work_candidate[i].final_volume_m3,
            )?;
        }
        for (port, connection) in self.connections.iter().enumerate() {
            poll(cancelled)?;
            let branch = &self.work.branch_report_candidate[connection.branch];
            let into_branch = match connection.end {
                BranchEnd::Left => &branch.left_boundary_change,
                BranchEnd::Right => &branch.right_boundary_change,
            };
            let mut species_masses_to_plenum = [0.0; MAX_SPECIES];
            for (value, &mass) in species_masses_to_plenum[..n]
                .iter_mut()
                .zip(&into_branch.species_masses[..n])
            {
                *value = -mass;
            }
            // Exact opposite of the accepted branch-end exchange, including
            // every species formation gauge and the local support impulse.
            self.work.port_exchange_candidate[port] = PortExchange {
                species_masses_to_plenum,
                species_count: n,
                momentum_to_plenum: -into_branch.momentum,
                energy_to_plenum: -into_branch.energy,
            };
        }
        let after = after.totals();
        let boundary_change = boundary.totals();
        let boundary_work_on_gas_j = sum(self
            .work
            .plenum_work_candidate
            .iter()
            .map(|x| x.work_on_gas_j));
        let mut species_masses = [0.0; MAX_SPECIES];
        for (k, value) in species_masses[..n].iter_mut().enumerate() {
            *value = after.species_masses[k]
                - before.species_masses[k]
                - boundary_change.species_masses[k];
        }
        let conservation_residual = NetworkTotals {
            species_masses,
            species_count: n,
            energy: after.energy - before.energy - boundary_change.energy - boundary_work_on_gas_j,
        };
        if !conservation_residual.energy.is_finite()
            || !conservation_residual.species_masses[..n]
                .iter()
                .all(|x| x.is_finite())
        {
            return Err(Error::Invalid(
                "non-finite complete mixture network balance",
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
