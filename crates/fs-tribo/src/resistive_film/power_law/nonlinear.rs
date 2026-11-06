//! Shared bounded Reynolds pressure solve and work-conjugate tangent.
//!
//! Constitutive leaves are statically dispatched. Geometry, pressure scaling,
//! Newton/backtracking, factorization, force projection and atomic publication
//! retain the original power-law implementation and numerical work bounds.

use super::{
    Channel, FilmError, FlowBalance, Geometry, MAX_CELLS, MAX_CHANNELS, MAX_PORTS,
    PowerLawControls, ReservoirPressures, ResistiveFilm, State, dot, fail, shear, solve_factored,
};

pub(super) trait ReynoldsLaw {
    fn pressure_operator(&self) -> &ResistiveFilm;
    fn pressure_controls(&self) -> PowerLawControls;
    fn native_newtonian(&self) -> bool;
    fn core_requires_positive_conductance(&self, index: usize) -> bool;
    #[allow(clippy::too_many_arguments)]
    fn channel_response(
        &self,
        index: usize,
        gap: f64,
        core_conductance: f64,
        pressure_drop: f64,
        motion: shear::Motion,
        moving: bool,
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<Channel, FilmError>;

    fn geometry(
        &self,
        q: &[f64],
        velocity: &[f64],
        walls: Option<&[shear::Coefficients]>,
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<Geometry, FilmError> {
        check()?;
        let minimum_gap = self.pressure_operator().validate_configuration(q)?;
        if velocity.len() != self.pressure_operator().ports
            || velocity.iter().any(|x| !x.is_finite())
        {
            return fail("invalid power-law film velocity/effort");
        }
        let mut geometry = Geometry {
            gaps: [0.0; MAX_CHANNELS],
            core_conductance: [0.0; MAX_CHANNELS],
            source: [0.0; MAX_CELLS],
            pressure_source: [0.0; MAX_CELLS],
            source_scale: 0.0,
            minimum_gap,
            motion: [shear::Motion::default(); MAX_CHANNELS],
        };
        let mut entrainment = [FlowBalance::default(); MAX_CELLS];
        for (k, edge) in self.pressure_operator().channels.iter().enumerate() {
            if k.is_multiple_of(32) {
                check()?;
            }
            let gap = edge.gap.reference_m - dot(&edge.gap.closure, q);
            geometry.gaps[k] = gap;
            let g = if gap > 0.0 {
                self.pressure_operator().factors[k] * gap * gap * gap
            } else {
                0.0
            };
            if !g.is_finite()
                || (gap > 0.0 && g == 0.0 && self.core_requires_positive_conductance(k))
            {
                return fail("power-law core differential conductance is unrepresentable");
            }
            geometry.core_conductance[k] = g;
            if let Some(wall) = walls.map(|walls| &walls[k]) {
                geometry.motion[k] = wall.motion(gap, edge.width_m, velocity)?;
                let flow = geometry.motion[k].entrainment;
                entrainment[edge.from].add(flow);
                if let Some(j) = edge.to {
                    entrainment[j].add(-flow);
                }
            }
        }
        for (i, cell) in self.pressure_operator().cells.iter().enumerate() {
            let source = cell.area_m2 * dot(&cell.gap.closure, velocity);
            if !source.is_finite() {
                return fail("power-law swept-volume rate is unrepresentable");
            }
            geometry.source[i] = source;
            let pressure_source = if walls.is_some() {
                source - entrainment[i].value()
            } else {
                source
            };
            if !pressure_source.is_finite() {
                return fail("power-law reduced pressure forcing is unrepresentable");
            }
            geometry.pressure_source[i] = pressure_source;
            geometry.source_scale = geometry.source_scale.max(pressure_source.abs());
        }
        Ok(geometry)
    }

    #[allow(clippy::too_many_arguments)]
    fn at_pressure(
        &self,
        geometry: &Geometry,
        pressure: &[f64; MAX_CELLS],
        reservoirs: Option<&ReservoirPressures>,
        walls: Option<&[shear::Coefficients]>,
        state: &mut State,
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<(), FilmError> {
        let n = self.pressure_operator().cell_count();
        state.system.lower[..n * n].fill(0.0);
        state.residual = geometry.pressure_source.map(|q| -q);
        state.maximum_residual = 0.0;
        state.power = 0.0;
        state.source_work = 0.0;
        state.wall_work = 0.0;
        state.reservoir_power = 0.0;
        state.supply_magnitude = 0.0;
        state.flow_scale = geometry.source_scale;
        state.system.pressure = *pressure;
        state.system.potential = *pressure;
        if let Some(boundary) = reservoirs {
            for p in &mut state.system.pressure[..n] {
                *p += boundary.datum.pressure_pa;
            }
        }
        if state.system.pressure[..n].iter().any(|p| {
            !p.is_finite() || p.abs() > self.pressure_operator().limits.maximum_pressure_pa
        }) {
            return fail("power-law pressure exceeds the declared incompressible domain");
        }
        for (k, edge) in self.pressure_operator().channels.iter().enumerate() {
            if k.is_multiple_of(32) {
                check()?;
            }
            let dp = pressure[edge.from]
                - edge.to.map_or_else(
                    || reservoirs.map_or(0.0, |r| r.potential[k]),
                    |j| pressure[j],
                );
            if !dp.is_finite() {
                return fail("power-law pressure drop is unrepresentable");
            }
            let moving = walls.is_some_and(|walls| walls[k].moving);
            let channel = self.channel_response(
                k,
                geometry.gaps[k],
                geometry.core_conductance[k],
                dp,
                geometry.motion[k],
                moving,
                check,
            )?;
            let pressure_flow = if moving {
                channel.wall.pressure_flow
            } else {
                channel.flow
            };
            let loss = if moving {
                let motion = geometry.motion[k];
                state.wall_work += channel.wall.mean_effort * motion.mean
                    + channel.wall.relative_effort * motion.relative;
                channel.wall.pressure_flow * dp + channel.wall.relative_effort * motion.relative
            } else {
                channel.flow * dp
            };
            // A source-owned profile may retain a finite local wall-speed
            // residual. Its accepted heat is physical; requested port work
            // remains separate and must satisfy the original work budget.
            let loss = channel.work.map_or(loss, |work| work.dissipated_power_w);
            if !loss.is_finite() || loss < 0.0 {
                return fail("nonfinite or negative power-law channel work");
            }
            state.power += loss;
            state.residual[edge.from] += pressure_flow;
            let g = channel.differential;
            let i = edge.from;
            state.system.lower[i * n + i] += g;
            if let Some(j) = edge.to {
                state.residual[j] -= pressure_flow;
                state.system.lower[j * n + j] += g;
                state.system.lower[i * n + j] -= g;
                state.system.lower[j * n + i] -= g;
            } else if let Some(boundary) = reservoirs {
                let supplied =
                    -boundary.reservoirs_pa[k].expect("admitted boundary") * channel.flow;
                state.reservoir_power += supplied;
                state.supply_magnitude += supplied.abs();
            }
            if reservoirs.is_some() || walls.is_some() {
                state.flow_scale = state.flow_scale.max(pressure_flow.abs());
            }
            state.channels[k] = channel;
        }
        for (i, &p) in state.system.pressure[..n].iter().enumerate() {
            state.source_work += p * geometry.source[i];
            state.maximum_residual = state.maximum_residual.max(state.residual[i].abs());
        }
        if !state.power.is_finite()
            || !state.source_work.is_finite()
            || !state.wall_work.is_finite()
            || !state.reservoir_power.is_finite()
            || !state.supply_magnitude.is_finite()
            || state.residual[..n].iter().any(|r| !r.is_finite())
        {
            return fail("power-law film volume or work sum is unrepresentable");
        }
        Ok(())
    }

    // Preserve the actual work-conjugate port projection for the force and
    // tangent gates, including its original floating-point accumulation.
    fn admit_port_work(
        &self,
        state: &mut State,
        velocity: &[f64],
        walls: Option<&[shear::Coefficients]>,
        supplied: bool,
    ) -> Result<(), FilmError> {
        self.pressure_operator()
            .project(&state.system.pressure, &mut state.force);
        if let Some(walls) = walls {
            for (k, coefficients) in walls.iter().enumerate() {
                if !coefficients.moving {
                    continue;
                }
                let wall = state.channels[k].wall;
                for (j, force) in state.force[..self.pressure_operator().ports]
                    .iter_mut()
                    .enumerate()
                {
                    *force += coefficients.mean[j] * wall.mean_effort
                        + coefficients.relative[j] * wall.relative_effort;
                }
            }
        }
        state.port_work = dot(velocity, &state.force[..self.pressure_operator().ports]);
        if state.force[..self.pressure_operator().ports]
            .iter()
            .any(|x| !x.is_finite())
            || !self
                .pressure_controls()
                .admits_supplied_work(state.port_work, state, supplied)
        {
            return fail("power-law film failed force and flow work balance");
        }
        Ok(())
    }

    fn flow_tolerance(&self, state: &State) -> Result<f64, FilmError> {
        let tolerance = self.pressure_controls().absolute_flow_tolerance_m3_s
            + self.pressure_controls().relative_flow_tolerance * state.flow_scale;
        if !tolerance.is_finite() {
            return fail("power-law flow tolerance is unrepresentable");
        }
        Ok(tolerance)
    }

    fn solve_with_reservoirs(
        &self,
        q: &[f64],
        velocity: &[f64],
        reservoirs: Option<&ReservoirPressures>,
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<State, FilmError> {
        self.solve_with_motion(q, velocity, reservoirs, None, check)
    }

    fn solve_with_motion(
        &self,
        q: &[f64],
        velocity: &[f64],
        reservoirs: Option<&ReservoirPressures>,
        walls: Option<&[shear::Coefficients]>,
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<State, FilmError> {
        let geometry = self.geometry(q, velocity, walls, check)?;
        let mut state = State::new(self.pressure_operator().cell_count(), geometry.minimum_gap);
        if self.native_newtonian() {
            let native = if walls.is_some() {
                let entrainment = geometry.motion.map(|motion| motion.entrainment);
                if let Some(boundary) = reservoirs {
                    self.pressure_operator()
                        .solve_with_boundaries::<true, true>(
                            q,
                            velocity,
                            &entrainment[..self.pressure_operator().channel_count()],
                            &boundary.potential,
                            boundary.datum,
                        )?
                } else {
                    self.pressure_operator().solve_with_entrainment::<true>(
                        q,
                        velocity,
                        &entrainment[..self.pressure_operator().channel_count()],
                    )?
                }
            } else if let Some(boundary) = reservoirs {
                self.pressure_operator()
                    .solve_with_boundaries::<false, true>(
                        q,
                        velocity,
                        &[],
                        &boundary.potential,
                        boundary.datum,
                    )?
            } else {
                self.pressure_operator().solve(q, velocity)?
            };
            check()?;
            self.at_pressure(
                &geometry,
                &native.potential,
                reservoirs,
                walls,
                &mut state,
                check,
            )?;
            state.system = native;
            if state.maximum_residual > self.flow_tolerance(&state)?
                || !self.pressure_controls().admits_supplied_work(
                    state.source_work + state.wall_work,
                    &state,
                    reservoirs.is_some(),
                )
            {
                return fail("native Newtonian limit exceeds power-law residual budgets");
            }
            self.admit_port_work(&mut state, velocity, walls, reservoirs.is_some())?;
            return Ok(state);
        }
        self.at_pressure(
            &geometry,
            &[0.0; MAX_CELLS],
            reservoirs,
            walls,
            &mut state,
            check,
        )?;
        // Two bounded workspaces avoid returning/copying the whole pressure
        // factorization and channel state for every line-search candidate.
        let mut trial = State::new(self.pressure_operator().cell_count(), geometry.minimum_gap);
        let mut backtracks = 0;
        for iteration in 0..=self.pressure_controls().maximum_iterations {
            check()?;
            // Connectivity is checked even for zero load: a sealed pocket does
            // not acquire an invented pressure datum just because its RHS is zero.
            state.system.factor()?;
            // A nonzero load always gets a pressure correction; an absolute
            // tolerance is not used as an artificial zero-force dead band.
            if (iteration > 0 || state.maximum_residual == 0.0)
                && state.maximum_residual <= self.flow_tolerance(&state)?
                && self.pressure_controls().admits_supplied_work(
                    state.source_work + state.wall_work,
                    &state,
                    reservoirs.is_some(),
                )
            {
                state.iterations = iteration;
                state.backtracks = backtracks;
                self.admit_port_work(&mut state, velocity, walls, reservoirs.is_some())?;
                return Ok(state);
            }
            if iteration == self.pressure_controls().maximum_iterations {
                return fail("power-law film pressure iteration budget exhausted");
            }
            let mut direction = state.residual.map(|r| -r);
            solve_factored(
                &state.system.lower,
                &state.system.scale,
                self.pressure_operator().cell_count(),
                &mut direction,
            )?;
            let mut fraction = 1.0_f64;
            let limit = self.pressure_operator().limits.maximum_pressure_pa;
            for (&p, &delta) in state.system.pressure[..self.pressure_operator().cell_count()]
                .iter()
                .zip(&direction)
            {
                let bound = if delta > 0.0 {
                    (limit - p) / delta
                } else if delta < 0.0 {
                    (-limit - p) / delta
                } else {
                    1.0
                };
                fraction = fraction.min(bound);
            }
            if fraction < 1.0 {
                fraction *= 0.99;
            }
            if !fraction.is_finite() || fraction <= 0.0 {
                return fail("power-law Newton step cannot remain inside the pressure domain");
            }
            let mut accepted = false;
            for _ in 0..self.pressure_controls().maximum_line_search_steps {
                check()?;
                let mut pressure = state.system.potential;
                for i in 0..self.pressure_operator().cell_count() {
                    pressure[i] += fraction * direction[i];
                }
                // A nonrepresentable physical trial can shorten the numerical
                // Newton step. Callback refusals must never be swallowed here.
                let mut callback_failure = None;
                let trial_result = self.at_pressure(
                    &geometry,
                    &pressure,
                    reservoirs,
                    walls,
                    &mut trial,
                    &mut || {
                        let result = check();
                        if let Err(error) = result {
                            callback_failure = Some(error);
                        }
                        result
                    },
                );
                if let Some(error) = callback_failure {
                    return Err(error);
                }
                check()?;
                if trial_result.is_ok() {
                    let decreased =
                        trial.maximum_residual <= state.maximum_residual * (1.0 - 1e-4 * fraction);
                    let converged = trial.maximum_residual <= self.flow_tolerance(&trial)?
                        && self.pressure_controls().admits_supplied_work(
                            trial.source_work + trial.wall_work,
                            &trial,
                            reservoirs.is_some(),
                        );
                    if decreased || converged {
                        core::mem::swap(&mut state, &mut trial);
                        accepted = true;
                        break;
                    }
                }
                backtracks += 1;
                fraction *= 0.5;
            }
            if !accepted {
                return fail("power-law film line-search budget exhausted");
            }
        }
        unreachable!("bounded pressure loop returns an accepted state or refusal")
    }

    #[allow(clippy::too_many_arguments)]
    fn tangent_with_walls(
        &self,
        q: &[f64],
        velocity: &[f64],
        dq: &[f64],
        dv: &[f64],
        reservoir_directions_pa: &[f64],
        out: &mut [f64],
        reservoirs: Option<&ReservoirPressures>,
        walls: &[shear::Coefficients],
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<(), FilmError> {
        let geometry = self.pressure_operator();
        if dq.len() != geometry.port_count()
            || dv.len() != geometry.port_count()
            || out.len() != geometry.port_count()
            || dq.iter().chain(dv).any(|x| !x.is_finite())
            || reservoir_directions_pa.len() != geometry.channel_count()
            || reservoir_directions_pa.iter().any(|x| !x.is_finite())
            || geometry
                .channels
                .iter()
                .zip(reservoir_directions_pa)
                .any(|(edge, &dp)| edge.to.is_some() && dp != 0.0)
        {
            return fail("invalid moving power-law film tangent direction");
        }
        let state = self.solve_with_motion(q, velocity, reservoirs, Some(walls), check)?;
        if state.channels[..geometry.channel_count()]
            .iter()
            .any(|c| c.wall.tangent_ambiguous)
        {
            return fail("power-law wall tangent is undefined at the pure-Couette core corner");
        }
        let mut rhs = [0.0; MAX_CELLS];
        let mut entrainment_change = [FlowBalance::default(); MAX_CELLS];
        for (i, cell) in geometry.cells.iter().enumerate() {
            rhs[i] = cell.area_m2 * dot(&cell.gap.closure, dv);
        }
        for (k, edge) in geometry.channels.iter().enumerate() {
            if k.is_multiple_of(32) {
                check()?;
            }
            let coefficients = &walls[k];
            let h = edge.gap.reference_m - dot(&edge.gap.closure, q);
            let dh = -dot(&edge.gap.closure, dq);
            let local = state.channels[k];
            let change = local.aperture_derivative * dh
                + local.wall.relative_flow_derivative * dot(&coefficients.relative, dv);
            rhs[edge.from] -= change;
            if let Some(j) = edge.to {
                rhs[j] += change;
            } else {
                rhs[edge.from] += local.differential * reservoir_directions_pa[k];
            }
            if coefficients.moving {
                let transported = edge.width_m
                    * (dh * dot(&coefficients.mean, velocity) + h * dot(&coefficients.mean, dv));
                entrainment_change[edge.from].add(transported);
                if let Some(j) = edge.to {
                    entrainment_change[j].add(-transported);
                }
            }
        }
        for (i, value) in rhs[..geometry.cell_count()].iter_mut().enumerate() {
            *value -= entrainment_change[i].value();
        }
        solve_factored(
            &state.system.lower,
            &state.system.scale,
            geometry.cell_count(),
            &mut rhs,
        )?;
        let mut candidate = [0.0; MAX_PORTS];
        geometry.project(&rhs, &mut candidate);
        for (k, edge) in geometry.channels.iter().enumerate() {
            let coefficients = &walls[k];
            if !coefficients.moving {
                continue;
            }
            let h = edge.gap.reference_m - dot(&edge.gap.closure, q);
            let dh = -dot(&edge.gap.closure, dq);
            let dp = state.system.potential[edge.from]
                - edge.to.map_or_else(
                    || reservoirs.map_or(0.0, |r| r.potential[k]),
                    |j| state.system.potential[j],
                );
            let ddp = rhs[edge.from] - edge.to.map_or(reservoir_directions_pa[k], |j| rhs[j]);
            let wall = state.channels[k].wall;
            let mean_change = -edge.width_m * (h * ddp + dp * dh);
            // Implicit stress derivatives share the cross coefficient Q_r:
            // (w L c)_p = -Q_r and (w L c)_h = -(dp Q_r + D r)/h.
            let relative_change = -wall.relative_flow_derivative * ddp
                + wall.relative_effort_derivative * dot(&coefficients.relative, dv)
                - (dp * wall.relative_flow_derivative
                    + wall.relative_effort_derivative * dot(&coefficients.relative, velocity))
                    * (dh / h);
            for (j, output) in candidate[..geometry.port_count()].iter_mut().enumerate() {
                *output +=
                    coefficients.mean[j] * mean_change + coefficients.relative[j] * relative_change;
            }
        }
        if candidate[..geometry.port_count()]
            .iter()
            .any(|x| !x.is_finite())
        {
            return fail("moving power-law film tangent is unrepresentable");
        }
        check()?;
        out.copy_from_slice(&candidate[..geometry.port_count()]);
        Ok(())
    }
}
