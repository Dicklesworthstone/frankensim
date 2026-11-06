//! Navier-slip walls on the original incompressible Reynolds pressure graph.
//!
//! The boundary law is `u(0)-Ulower=b_lower*u'(0)` and
//! `u(h)-Uupper=-b_upper*u'(h)`, with prescribed nonnegative slip lengths.
//! For `D=h+b_lower+b_upper`, the pressure-shear neutral fraction is
//! `a=(h/2+b_upper)/D`, entrainment is `w*h*(a*Ulower+(1-a)*Uupper)`,
//! and relative-wall resistance is `w*L*mu/D`. Pressure mobility, both wall
//! tractions and their geometry derivatives follow from this SAME profile.
//! The scalar parallel-plate law is described by Bazant and Vinogradova,
//! J. Fluid Mech. 613 (2008), <https://doi.org/10.1017/S002211200800356X>.
//!
//! Irreversible work includes bulk viscous heat and two distinct interfacial
//! losses `area*b*tau^2/mu`. It is not all fluid heat. A thermal owner chooses
//! where each interfacial loss goes; this operator owns no thermal stores.
//! The common Barus viscosity law is retained with FIXED slip lengths, so its
//! exact length-average multiplier applies to stress and each heat integral.
//! This is a quasistatic scalar slip model, not a rarefied-gas, textured-wall
//! tensor, non-Newtonian, stratified-viscosity or cavitating film model.

use super::{
    FilmError, FilmReport, MAX_CELLS, MAX_CHANNELS, MAX_PORTS, PressureDatum, PressureDrivenReport,
    ReservoirPressures, ShearChannelReport, ShearedFilm, Solution, Sum, barus_mean_at_datum, dot,
    fail, pressure_drop_at_datum, solve_factored,
};

/// Prescribed scalar Navier slip lengths for one original channel [m].
///
/// Zero imposes no slip at that wall. Values must come from the admitted
/// interface model; material identity or roughness does not supply a default.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct WallSlipLengths {
    /// Slip length at the lower wall [m].
    pub lower_m: f64,
    /// Slip length at the upper wall [m].
    pub upper_m: f64,
}

/// Flow, traction and physically distinct heat destinations in one passage.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SlipChannelReport {
    /// Complete pressure/Couette flow, resisting wall traction and total work.
    /// Its irreversible powers include bulk AND interfacial dissipation.
    pub film: ShearChannelReport,
    /// Viscous dissipation inside the fluid volume [W].
    pub bulk_dissipated_power_w: f64,
    /// Nonnegative dissipation at the lower sliding interface [W].
    pub lower_slip_dissipated_power_w: f64,
    /// Nonnegative dissipation at the upper sliding interface [W].
    pub upper_slip_dissipated_power_w: f64,
    /// Fluid velocity minus lower-wall velocity, along the channel [m/s].
    pub lower_slip_velocity_m_s: f64,
    /// Fluid velocity minus upper-wall velocity, along the channel [m/s].
    pub upper_slip_velocity_m_s: f64,
}

/// Pressure elimination with unequal Navier-slip walls and a complete tangent.
///
/// Each channel has a uniform cross-gap Newtonian viscosity. Geometry,
/// original mechanical ports and pressure limits remain owned by the input
/// film. Evaluation and tangents use bounded stack scratch and one native
/// pressure factorization. Slip lengths stay fixed during each evaluation.
#[derive(Clone, Debug)]
pub struct NavierSlipFilm {
    film: ShearedFilm,
    slips: Vec<WallSlipLengths>,
    reservoirs: Option<ReservoirPressures>,
}

impl NavierSlipFilm {
    /// Bind one explicit pair of slip lengths to every original channel.
    ///
    /// A hydraulic-equivalent viscosity cannot specify wall slip in an
    /// unresolved layered film; retained stratified inputs therefore refuse.
    pub fn new(film: ShearedFilm, slips: Vec<WallSlipLengths>) -> Result<Self, FilmError> {
        if film.stratified_profiles().is_some() {
            return fail("Navier-slip film requires uniform cross-gap channel viscosity");
        }
        if slips.len() != film.film.channel_count()
            || slips.iter().any(|slip| {
                !slip.lower_m.is_finite()
                    || !slip.upper_m.is_finite()
                    || slip.lower_m < 0.0
                    || slip.upper_m < 0.0
                    || !(slip.lower_m + slip.upper_m).is_finite()
            })
        {
            return fail(
                "Navier-slip film needs finite nonnegative slip lengths for every channel",
            );
        }
        Ok(Self {
            film,
            slips,
            reservoirs: None,
        })
    }

    /// Original geometry and wall-port binding, before the slip modification.
    #[must_use]
    pub const fn sheared_film(&self) -> &ShearedFilm {
        &self.film
    }

    /// Slip lengths in original channel order.
    #[must_use]
    pub fn wall_slip_lengths(&self) -> &[WallSlipLengths] {
        &self.slips
    }

    /// Bind prescribed gauge pressures to every original reservoir channel.
    /// Boundary channels need `Some(p)` and internal channels need `None`.
    pub fn with_reservoir_pressures(
        &self,
        pressures_pa: &[Option<f64>],
    ) -> Result<Self, FilmError> {
        let reservoirs = ReservoirPressures::new(&self.film.film, pressures_pa)?;
        Ok(Self {
            reservoirs: Some(reservoirs),
            ..self.clone()
        })
    }

    /// Rebind original channel viscosities while retaining geometry and slip.
    /// Slip LENGTHS are held fixed; interfacial mobility `b/mu` then changes.
    pub fn with_channel_viscosities(&self, viscosities_pa_s: &[f64]) -> Result<Self, FilmError> {
        let rebound = Self::new(
            self.film.with_layer_viscosities(viscosities_pa_s)?,
            self.slips.clone(),
        )?;
        match &self.reservoirs {
            Some(boundary) => rebound.with_reservoir_pressures(&boundary.reservoirs_pa),
            None => Ok(rebound),
        }
    }

    /// Evaluate force, pressure, wall slip and the three physical heat shares.
    /// Every output remains untouched on refusal.
    pub fn evaluate_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [SlipChannelReport],
    ) -> Result<PressureDrivenReport, FilmError> {
        self.evaluate::<true>(q, velocity, force, pressure, channels)
    }

    /// The same work-conjugate force without channel-report computation.
    pub fn evaluate_force_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
    ) -> Result<PressureDrivenReport, FilmError> {
        self.evaluate::<false>(q, velocity, force, pressure, &mut [])
    }

    fn evaluate<const REPORT: bool>(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [SlipChannelReport],
    ) -> Result<PressureDrivenReport, FilmError> {
        let film = &self.film.film;
        if force.len() != film.ports
            || pressure.len() != film.cell_count()
            || (REPORT && channels.len() != film.channel_count())
        {
            return fail("invalid Navier-slip film output shape");
        }
        let state = self.state(q, velocity)?;
        let mut candidate = [0.0; MAX_PORTS];
        let mut reports = [SlipChannelReport::default(); MAX_CHANNELS];
        film.project(&state.pressure.pressure, &mut candidate);
        let mut total = 0.0;
        let mut supply = Sum::default();
        let mut supply_magnitude = 0.0;
        for (k, edge) in film.channels.iter().enumerate() {
            let local = state.channels[k];
            let potential_drop = self.potential_drop(&state.pressure, k);
            let pressure_flow = state.pressure.conductance[k] * potential_drop;
            let pressure_heat = pressure_flow * local.pressure_drop;
            let shear_force = local.drag * local.relative_speed;
            let shear_heat = shear_force * local.relative_speed;
            let flow = pressure_flow + local.entrainment;
            if !flow.is_finite()
                || !pressure_heat.is_finite()
                || pressure_heat < 0.0
                || !shear_heat.is_finite()
                || shear_heat < 0.0
            {
                return fail("Navier-slip film has nonfinite or negative channel work");
            }
            total += pressure_heat + shear_heat;
            if self.film.coefficients[k].moving {
                let pressure_force = edge.width_m * local.gap * local.pressure_drop;
                for (j, output) in candidate[..film.ports].iter_mut().enumerate() {
                    *output += -self.mean_row(k, j, local) * pressure_force
                        + self.film.coefficients[k].relative_velocity[j] * shear_force;
                }
            }
            if let Some(p) = self.reservoirs.as_ref().and_then(|r| r.reservoirs_pa[k]) {
                let incoming = -p * flow;
                supply.add(incoming);
                supply_magnitude += incoming.abs();
            }
            if REPORT && local.gap > 0.0 {
                let area = edge.width_m * edge.length_m;
                let mu = film.channel_viscosities_pa_s[k] * local.multiplier;
                let pressure_stress = local.gap / edge.length_m * local.pressure_drop;
                let couette_stress = (shear_force / edge.width_m) / edge.length_m;
                let lower = -local.lower_fraction * pressure_stress - couette_stress;
                let upper = -local.upper_fraction * pressure_stress + couette_stress;
                let lower_slip = -(self.slips[k].lower_m / mu) * lower;
                let upper_slip = -(self.slips[k].upper_m / mu) * upper;
                let mid_stress = couette_stress + (local.lower_fraction - 0.5) * pressure_stress;
                let bulk = (area * local.gap / mu)
                    * (mid_stress * mid_stress + pressure_stress * pressure_stress / 12.0);
                let lower_heat = -area * lower * lower_slip;
                let upper_heat = -area * upper * upper_slip;
                let heat = pressure_heat + shear_heat;
                let split = bulk + lower_heat + upper_heat;
                let wall_power = shear_heat - local.entrainment * local.pressure_drop;
                if [
                    lower, upper, lower_slip, upper_slip, bulk, lower_heat, upper_heat, wall_power,
                ]
                .iter()
                .any(|x| !x.is_finite())
                    || bulk < 0.0
                    || lower_heat < 0.0
                    || upper_heat < 0.0
                    || !split.is_finite()
                    || (split - heat).abs() > 1e-9 * (split.abs() + heat.abs()).max(1e-30)
                {
                    return fail("Navier-slip film bulk and interface work do not close");
                }
                reports[k] = SlipChannelReport {
                    film: ShearChannelReport {
                        volume_flow_m3_s: flow,
                        entrained_flow_m3_s: local.entrainment,
                        lower_wall_traction_pa: lower,
                        upper_wall_traction_pa: upper,
                        wall_power_w: wall_power,
                        pressure_dissipated_power_w: pressure_heat,
                        shear_dissipated_power_w: shear_heat,
                    },
                    bulk_dissipated_power_w: bulk,
                    lower_slip_dissipated_power_w: lower_heat,
                    upper_slip_dissipated_power_w: upper_heat,
                    lower_slip_velocity_m_s: lower_slip,
                    upper_slip_velocity_m_s: upper_slip,
                };
            }
        }
        let work = dot(velocity, &candidate[..film.ports]);
        let supply = supply.value();
        let residual = work + supply - total;
        let magnitude = work.abs() + supply_magnitude + total.abs();
        if candidate[..film.ports].iter().any(|x| !x.is_finite())
            || !total.is_finite()
            || !work.is_finite()
            || !magnitude.is_finite()
            || !residual.is_finite()
            || residual.abs() > 1e-9 * magnitude.max(1e-30)
        {
            return fail("Navier-slip film failed pressure and wall work balance");
        }
        force.copy_from_slice(&candidate[..film.ports]);
        pressure.copy_from_slice(&state.pressure.pressure[..film.cell_count()]);
        if REPORT {
            channels.copy_from_slice(&reports[..film.channel_count()]);
        }
        Ok(PressureDrivenReport {
            film: FilmReport {
                minimum_gap_m: state.pressure.minimum_gap,
                peak_pressure_pa: pressure.iter().fold(0.0_f64, |m, p| m.max(p.abs())),
                dissipated_power_w: total,
            },
            mechanical_power_w: work,
            reservoir_power_w: supply,
            power_residual_w: residual,
        })
    }

    /// Full directional force derivative at fixed slip lengths and supplies.
    ///
    /// Includes the aperture dependence of hydraulic mobility, asymmetric
    /// entrainment, wall pressure traction, shear resistance and Barus viscosity.
    /// The accepted pressure factorization is reused for the tangent solve.
    pub fn tangent_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        dq: &[f64],
        dv: &[f64],
        output: &mut [f64],
    ) -> Result<(), FilmError> {
        let film = &self.film.film;
        if dq.len() != film.ports
            || dv.len() != film.ports
            || output.len() != film.ports
            || dq.iter().chain(dv).any(|x| !x.is_finite())
        {
            return fail("invalid Navier-slip film tangent direction");
        }
        let state = self.state(q, velocity)?;
        let mut rhs = [0.0; MAX_CELLS];
        for (i, cell) in film.cells.iter().enumerate() {
            rhs[i] = cell.area_m2 * dot(&cell.gap.closure, dv);
        }
        for (k, edge) in film.channels.iter().enumerate() {
            let local = state.channels[k];
            if local.gap <= 0.0 {
                continue;
            }
            let dh = -dot(&edge.gap.closure, dq);
            let dmean = (0..film.ports)
                .map(|j| self.mean_row(k, j, local) * dv[j])
                .sum::<f64>()
                - local.lower_fraction_dh * dh * local.relative_speed;
            let dq_couette = edge.width_m * (dh * local.mean_speed + local.gap * dmean);
            let change =
                local.conductance_dh * dh * self.potential_drop(&state.pressure, k) + dq_couette;
            rhs[edge.from] -= change;
            if let Some(j) = edge.to {
                rhs[j] += change;
            }
        }
        solve_factored(
            &state.pressure.lower,
            &state.pressure.scale,
            film.cell_count(),
            &mut rhs,
        )?;
        let alpha = film.pressure_viscosity_coefficient_per_pa;
        let datum = self.datum();
        let mut pressure_direction = [0.0; MAX_CELLS];
        for i in 0..film.cell_count() {
            pressure_direction[i] = rhs[i] / (datum.fluidity - alpha * state.pressure.potential[i]);
        }
        let mut candidate = [0.0; MAX_PORTS];
        film.project(&pressure_direction, &mut candidate);
        for (k, edge) in film.channels.iter().enumerate() {
            let local = state.channels[k];
            if !self.film.coefficients[k].moving {
                continue;
            }
            let dh = -dot(&edge.gap.closure, dq);
            let dp = pressure_direction[edge.from] - edge.to.map_or(0.0, |j| pressure_direction[j]);
            let drelative = dot(&self.film.coefficients[k].relative_velocity, dv);
            let (_, dmultiplier) = barus_mean_at_datum(
                alpha,
                datum.fluidity,
                state.pressure.potential[edge.from],
                self.to_potential(&state.pressure, k),
                rhs[edge.from],
                edge.to.map_or(0.0, |j| rhs[j]),
            )?;
            let dshear = local.drag * (drelative - local.relative_speed * dh / local.separation)
                + local.base_drag * local.relative_speed * dmultiplier;
            let pressure_force = edge.width_m * local.gap * local.pressure_drop;
            let dpressure_force = edge.width_m * (dh * local.pressure_drop + local.gap * dp);
            for (j, value) in candidate[..film.ports].iter_mut().enumerate() {
                let relative_row = self.film.coefficients[k].relative_velocity[j];
                *value += -self.mean_row(k, j, local) * dpressure_force
                    + relative_row * (dshear + local.lower_fraction_dh * dh * pressure_force);
            }
        }
        if candidate[..film.ports].iter().any(|x| !x.is_finite()) {
            return fail("Navier-slip film tangent is not representable");
        }
        output.copy_from_slice(&candidate[..film.ports]);
        Ok(())
    }

    fn state(&self, q: &[f64], velocity: &[f64]) -> Result<State, FilmError> {
        let film = &self.film.film;
        if q.len() != film.ports
            || velocity.len() != film.ports
            || q.iter().chain(velocity).any(|x| !x.is_finite())
        {
            return fail("invalid Navier-slip film configuration or effort");
        }
        let mut channels = [Passage::default(); MAX_CHANNELS];
        let mut conductance = [0.0; MAX_CHANNELS];
        let mut entrainment = [0.0; MAX_CHANNELS];
        for (k, edge) in film.channels.iter().enumerate() {
            let h = edge.gap.reference_m - dot(&edge.gap.closure, q);
            if !h.is_finite() || (self.film.coefficients[k].moving && h <= 0.0) {
                return fail("moving Navier-slip walls require a positive film aperture");
            }
            if h <= 0.0 {
                continue;
            }
            let slip = self.slips[k];
            let separation = h + (slip.lower_m + slip.upper_m);
            let lower_fraction = (0.5 * h + slip.upper_m) / separation;
            let upper_fraction = (0.5 * h + slip.lower_m) / separation;
            let offset = 0.5 * (slip.lower_m - slip.upper_m) / separation;
            // Positive second moments of the fluid and two wall mobilities.
            // Their centroid derivative cancels, so dK/dh is the fluid moment.
            let fluid_moment = offset * offset + 1.0 / 12.0;
            let moment = h * fluid_moment
                + slip.lower_m * lower_fraction * lower_fraction
                + slip.upper_m * upper_fraction * upper_fraction;
            let g = if slip.lower_m == 0.0 && slip.upper_m == 0.0 {
                film.factors[k] * h * h * h
            } else {
                (12.0 * film.factors[k]) * h * h * moment
            };
            let dg = g * (2.0 / h + fluid_moment / moment);
            let relative = dot(&self.film.coefficients[k].relative_velocity, velocity);
            let base_drag = self.film.coefficients[k].shear_scale / separation;
            let mut local = Passage {
                gap: h,
                separation,
                lower_fraction,
                upper_fraction,
                lower_fraction_dh: offset / separation,
                relative_speed: relative,
                base_drag,
                conductance_dh: dg,
                ..Passage::default()
            };
            let mean = (0..film.ports)
                .map(|j| self.mean_row(k, j, local) * velocity[j])
                .sum::<f64>();
            let flow = edge.width_m * h * mean;
            if [
                separation,
                lower_fraction,
                upper_fraction,
                moment,
                g,
                dg,
                relative,
                base_drag,
                mean,
                flow,
            ]
            .iter()
            .any(|x| !x.is_finite())
                || g <= 0.0
                || moment <= 0.0
            {
                return fail("Navier-slip film mobility or wall motion is not representable");
            }
            local.mean_speed = mean;
            local.entrainment = flow;
            channels[k] = local;
            conductance[k] = g;
            entrainment[k] = flow;
        }
        let pressure = match &self.reservoirs {
            Some(boundary) => film.solve_with_mobility::<true, true, true>(
                q,
                velocity,
                &entrainment[..film.channel_count()],
                &boundary.potential,
                boundary.datum,
                &conductance[..film.channel_count()],
            )?,
            None => film.solve_with_mobility::<true, false, true>(
                q,
                velocity,
                &entrainment[..film.channel_count()],
                &[],
                PressureDatum::ZERO,
                &conductance[..film.channel_count()],
            )?,
        };
        let alpha = film.pressure_viscosity_coefficient_per_pa;
        for (k, edge) in film.channels.iter().enumerate() {
            if channels[k].gap <= 0.0 {
                continue;
            }
            let from = pressure.potential[edge.from];
            let to = self.to_potential(&pressure, k);
            let multiplier =
                barus_mean_at_datum(alpha, self.datum().fluidity, from, to, 0.0, 0.0)?.0;
            let dp = pressure_drop_at_datum(alpha, self.datum().fluidity, from, to)?;
            let drag = channels[k].base_drag * multiplier;
            if !drag.is_finite() || drag < 0.0 {
                return fail("Navier-slip film wall resistance is not representable");
            }
            channels[k].pressure_drop = dp;
            channels[k].multiplier = multiplier;
            channels[k].drag = drag;
        }
        Ok(State { pressure, channels })
    }

    fn mean_row(&self, channel: usize, port: usize, local: Passage) -> f64 {
        let wall = &self.film.walls[channel];
        let lower = wall.lower[port];
        let upper = wall.upper[port];
        let relative = self.film.coefficients[channel].relative_velocity[port];
        let combined = local.lower_fraction + local.upper_fraction;
        if lower.abs() <= upper.abs() {
            local.upper_fraction.mul_add(relative, combined * lower)
        } else {
            (-local.lower_fraction).mul_add(relative, combined * upper)
        }
    }

    fn datum(&self) -> PressureDatum {
        self.reservoirs
            .as_ref()
            .map_or(PressureDatum::ZERO, |r| r.datum)
    }

    fn to_potential(&self, pressure: &Solution, channel: usize) -> f64 {
        self.film.film.channels[channel].to.map_or_else(
            || {
                self.reservoirs
                    .as_ref()
                    .map_or(0.0, |r| r.potential[channel])
            },
            |j| pressure.potential[j],
        )
    }

    fn potential_drop(&self, pressure: &Solution, channel: usize) -> f64 {
        pressure.potential[self.film.film.channels[channel].from]
            - self.to_potential(pressure, channel)
    }
}

#[derive(Clone, Copy, Default)]
struct Passage {
    gap: f64,
    separation: f64,
    lower_fraction: f64,
    upper_fraction: f64,
    lower_fraction_dh: f64,
    relative_speed: f64,
    mean_speed: f64,
    entrainment: f64,
    base_drag: f64,
    drag: f64,
    conductance_dh: f64,
    pressure_drop: f64,
    multiplier: f64,
}

struct State {
    pressure: Solution,
    channels: [Passage; MAX_CHANNELS],
}

#[cfg(test)]
mod tests;
