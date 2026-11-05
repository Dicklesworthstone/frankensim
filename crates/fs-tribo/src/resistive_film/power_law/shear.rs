//! Simultaneous Couette--Poiseuille flow of the original power-law liquid.
//!
//! With x in [-1,1], the shear stress is c - delta_p*h*x/(2*length).
//! Its offset c is solved so the integrated inverse stress law gives the
//! prescribed upper-minus-lower wall speed. Pressure flow and sliding drag
//! therefore sample ONE nonlinear velocity profile. Adding a pressure-only
//! power-law flux to a separately evaluated Couette drag would be incorrect.
//!
//! Mean translation contributes width*h*Umean to flow and the reciprocal
//! pressure traction to wall force. It contributes no heat. The positive
//! weighted variance of the inverse law's differential gives dQ/d(delta_p);
//! the existing pressure graph, bounded Newton solve and factorization remain
//! the numerical owners. All velocity rows use the native mechanical ports.

pub use super::super::shear::WallVelocityPorts;
use super::{
    FilmError, FilmReport, FlowBalance, MAX_CELLS, MAX_CHANNELS, MAX_PORTS, PowerLawFilmReport,
    PressureDrivenPowerLawFilm, PressureDrivenPowerLawReport, RegularizedPowerLaw, dot, fail,
    solve_factored,
};

pub(super) mod channel;

/// Flow, wall traction and irreversible power of the same nonlinear passage.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ShearedPowerLawChannelReport {
    /// Signed total flow, positive from the channel's original `from` cell.
    pub volume_flow_m3_s: f64,
    /// Signed excess over mean translation from the same nonlinear profile.
    /// Retained separately so large entrainment cannot erase this small flow.
    pub pressure_flow_m3_s: f64,
    /// Mean-wall translation alone, width*h*(Ulower+Uupper)/2 [m^3/s].
    pub entrained_flow_m3_s: f64,
    /// Lower-wall resisting traction along the channel orientation [Pa].
    /// Multiplying by width*length gives its work-conjugate force.
    pub lower_wall_traction_pa: f64,
    /// Upper-wall resisting traction along the channel orientation [Pa].
    pub upper_wall_traction_pa: f64,
    /// Signed input from both moving walls, including pressure traction [W].
    pub wall_power_w: f64,
    /// Nonnegative (Q-Qentrained)*delta_p of the COMBINED profile [W].
    pub pressure_dissipated_power_w: f64,
    /// Nonnegative width*length*c*(Uupper-Ulower) of that same profile [W].
    /// This is not a separately evaluated pure-Couette constitutive law.
    pub shear_dissipated_power_w: f64,
}

impl ShearedPowerLawChannelReport {
    /// Physical viscous heat, counted once; signed supply/wall work is separate.
    #[must_use]
    pub fn dissipated_power_w(self) -> f64 {
        self.pressure_dissipated_power_w + self.shear_dissipated_power_w
    }
}

/// A pressure-fed power-law film with explicitly authored lower/upper walls.
///
/// Geometry, supplies, laws and wall rows are immutable during evaluation.
/// Both pressure and local stress corrections use the original declared
/// iteration ceiling. Local stress roots resolve wall velocity to roundoff;
/// original cell-flow and complete port-work budgets still gate publication.
/// Calls allocate no heap scratch and poll the supplied check inside both
/// nonlinear solves. A failed solve never publishes partial caller outputs.
///
/// This remains a homogeneous, incompressible, full-film Reynolds reduction.
/// It adds no slip, layer evolution, thermal clock, cavitation or fluid inertia.
#[derive(Debug, Clone)]
pub struct ShearedPowerLawFilm {
    film: PressureDrivenPowerLawFilm,
    walls: Vec<WallVelocityPorts>,
    coefficients: Vec<Coefficients>,
    moving: bool,
}

#[derive(Debug, Clone)]
pub(super) struct Coefficients {
    pub(super) mean: Vec<f64>,
    pub(super) relative: Vec<f64>,
    pub(super) moving: bool,
}

#[derive(Clone, Copy, Default)]
pub(super) struct Motion {
    pub(super) mean: f64,
    pub(super) relative: f64,
    pub(super) entrainment: f64,
}

#[derive(Clone, Copy, Default)]
pub(super) struct WallResponse {
    pub(super) mean_effort: f64,
    pub(super) relative_effort: f64,
    pub(super) relative_flow_derivative: f64,
    pub(super) relative_effort_derivative: f64,
    pub(super) pressure_flow: f64,
    pub(super) mid_stress: f64,
    pub(super) tangent_ambiguous: bool,
}

impl Coefficients {
    pub(super) fn motion(
        &self,
        gap: f64,
        width: f64,
        velocity: &[f64],
    ) -> Result<Motion, FilmError> {
        if !self.moving {
            return Ok(Motion::default());
        }
        if gap <= 0.0 {
            return fail("moving power-law passage needs a positive gap");
        }
        let mean = dot(&self.mean, velocity);
        let relative = dot(&self.relative, velocity);
        let entrainment = (width * gap) * mean;
        if [mean, relative, entrainment].iter().any(|x| !x.is_finite()) {
            return fail("power-law wall motion is unrepresentable");
        }
        Ok(Motion {
            mean,
            relative,
            entrainment,
        })
    }
}

impl ShearedPowerLawFilm {
    /// Bind native wall-velocity rows in original hydraulic-channel order.
    /// Zero rows retain fixed walls; a nonzero row requires an open passage
    /// even when the particular evaluated wall velocity happens to be zero.
    pub fn new(
        film: PressureDrivenPowerLawFilm,
        walls: Vec<WallVelocityPorts>,
    ) -> Result<Self, FilmError> {
        let geometry = film.power_law_film().pressure_operator();
        if walls.len() != geometry.channel_count()
            || walls.iter().any(|wall| {
                wall.lower.len() != geometry.port_count()
                    || wall.upper.len() != geometry.port_count()
                    || wall.lower.iter().chain(&wall.upper).any(|x| !x.is_finite())
            })
        {
            return fail("power-law sliding film needs finite original-order wall rows");
        }
        let mut coefficients = Vec::with_capacity(walls.len());
        let mut moving = false;
        for wall in &walls {
            let mut mean = Vec::with_capacity(geometry.port_count());
            let mut relative = Vec::with_capacity(geometry.port_count());
            for (&lower, &upper) in wall.lower.iter().zip(&wall.upper) {
                let difference = upper - lower;
                let average = if lower.abs() <= upper.abs() {
                    0.5_f64.mul_add(difference, lower)
                } else {
                    (-0.5_f64).mul_add(difference, upper)
                };
                if !difference.is_finite() || !average.is_finite() {
                    return fail("power-law wall mean or relative row is unrepresentable");
                }
                mean.push(average);
                relative.push(difference);
            }
            let active = wall.lower.iter().chain(&wall.upper).any(|&x| x != 0.0);
            moving |= active;
            coefficients.push(Coefficients {
                mean,
                relative,
                moving: active,
            });
        }
        Ok(Self {
            film,
            walls,
            coefficients,
            moving,
        })
    }

    /// Constitutive owner with its authored pressure supplies; its evaluations
    /// omit moving walls, so they are not substitutes for this binding.
    #[must_use]
    pub const fn pressure_film(&self) -> &PressureDrivenPowerLawFilm {
        &self.film
    }

    /// Exact lower/upper velocity rows, with no velocity inferred from a gap.
    #[must_use]
    pub fn wall_ports(&self) -> &[WallVelocityPorts] {
        &self.walls
    }

    /// Refresh supplies while retaining all wall rows, laws and controls.
    pub fn with_reservoir_pressures(
        &self,
        pressures_pa: &[Option<f64>],
    ) -> Result<Self, FilmError> {
        Self::new(
            self.film.with_reservoir_pressures(pressures_pa)?,
            self.walls.clone(),
        )
    }

    /// Refresh explicit channel laws without advancing a thermal state.
    pub fn with_laws(&self, laws: Vec<RegularizedPowerLaw>) -> Result<Self, FilmError> {
        Self::new(self.film.with_laws(laws)?, self.walls.clone())
    }

    /// Solve the complete normal and tangential mechanical resistance.
    /// The report's mechanical power includes BOTH normal and wall ports.
    pub fn evaluate_force_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
    ) -> Result<PressureDrivenPowerLawReport, FilmError> {
        self.evaluate::<false>(q, velocity, force, pressure, &mut [], &mut || Ok(()))
    }

    /// Publish force, physical pressure, passage flow, traction and heat together.
    pub fn evaluate_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [ShearedPowerLawChannelReport],
    ) -> Result<PressureDrivenPowerLawReport, FilmError> {
        self.evaluate_with_check_into(q, velocity, force, pressure, channels, || Ok(()))
    }

    /// The same complete solve with cancellation checks and atomic publication.
    pub fn evaluate_with_check_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [ShearedPowerLawChannelReport],
        mut check: impl FnMut() -> Result<(), FilmError>,
    ) -> Result<PressureDrivenPowerLawReport, FilmError> {
        self.evaluate::<true>(q, velocity, force, pressure, channels, &mut check)
    }

    fn evaluate<const REPORT: bool>(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [ShearedPowerLawChannelReport],
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<PressureDrivenPowerLawReport, FilmError> {
        let original = self.film.power_law_film();
        let geometry = original.pressure_operator();
        if force.len() != geometry.port_count()
            || pressure.len() != geometry.cell_count()
            || (REPORT && channels.len() != geometry.channel_count())
        {
            return fail("invalid moving power-law film output shape");
        }
        let reservoirs = self.film.evaluation_reservoirs();
        let state = original.solve_with_motion(
            q,
            velocity,
            reservoirs,
            self.moving.then_some(&self.coefficients),
            check,
        )?;
        let mut candidate = [ShearedPowerLawChannelReport::default(); MAX_CHANNELS];
        if REPORT {
            for (k, edge) in geometry.channels.iter().enumerate() {
                let local = state.channels[k];
                let wall = local.wall;
                let dp = state.system.potential[edge.from]
                    - edge.to.map_or_else(
                        || reservoirs.map_or(0.0, |r| r.potential[k]),
                        |j| state.system.potential[j],
                    );
                let h = (edge.gap.reference_m - dot(&edge.gap.closure, q)).max(0.0);
                let motion = self.coefficients[k].motion(h, edge.width_m, velocity)?;
                let s = dp * (0.5 * h / edge.length_m);
                candidate[k] = ShearedPowerLawChannelReport {
                    volume_flow_m3_s: local.flow,
                    pressure_flow_m3_s: if self.coefficients[k].moving {
                        wall.pressure_flow
                    } else {
                        local.flow
                    },
                    entrained_flow_m3_s: motion.entrainment,
                    lower_wall_traction_pa: -(wall.mid_stress + s),
                    upper_wall_traction_pa: wall.mid_stress - s,
                    wall_power_w: wall.mean_effort * motion.mean
                        + wall.relative_effort * motion.relative,
                    pressure_dissipated_power_w: if self.coefficients[k].moving {
                        wall.pressure_flow * dp
                    } else {
                        local.flow * dp
                    },
                    shear_dissipated_power_w: wall.relative_effort * motion.relative,
                };
                if [
                    candidate[k].lower_wall_traction_pa,
                    candidate[k].upper_wall_traction_pa,
                    candidate[k].wall_power_w,
                ]
                .iter()
                .any(|x| !x.is_finite())
                {
                    return fail("power-law wall traction or diagnostic work is unrepresentable");
                }
            }
        }
        let report = PressureDrivenPowerLawReport {
            power_law: PowerLawFilmReport {
                film: FilmReport {
                    minimum_gap_m: state.system.minimum_gap,
                    peak_pressure_pa: state.system.pressure[..geometry.cell_count()]
                        .iter()
                        .fold(0.0_f64, |m, p| m.max(p.abs())),
                    dissipated_power_w: state.power,
                },
                iterations: state.iterations,
                backtracks: state.backtracks,
                maximum_flow_defect_m3_s: state.maximum_residual,
                work_defect_w: state.port_work + state.reservoir_power - state.power,
            },
            mechanical_power_w: state.port_work,
            reservoir_power_w: state.reservoir_power,
        };
        check()?;
        force.copy_from_slice(&state.force[..geometry.port_count()]);
        pressure.copy_from_slice(&state.system.pressure[..geometry.cell_count()]);
        if REPORT {
            channels.copy_from_slice(&candidate[..geometry.channel_count()]);
        }
        Ok(report)
    }

    /// Differentiate geometry, all wall/normal velocities and pressure supplies.
    /// Implicit stress-offset derivatives and the accepted pressure Jacobian are
    /// reused; no velocity-profile or pressure finite differences are taken.
    /// At the exact uniform pure-Couette core transition the authored stress
    /// law has no unique full tangent; this method refuses without publication.
    #[allow(clippy::too_many_arguments)]
    pub fn tangent_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        dq: &[f64],
        dv: &[f64],
        reservoir_directions_pa: &[f64],
        out: &mut [f64],
    ) -> Result<(), FilmError> {
        self.tangent_with_check_into(q, velocity, dq, dv, reservoir_directions_pa, out, || Ok(()))
    }

    /// Same analytic tangent, including checks before publishing caller output.
    #[allow(clippy::too_many_arguments)]
    pub fn tangent_with_check_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        dq: &[f64],
        dv: &[f64],
        reservoir_directions_pa: &[f64],
        out: &mut [f64],
        mut check: impl FnMut() -> Result<(), FilmError>,
    ) -> Result<(), FilmError> {
        if !self.moving {
            return self.film.tangent_with_check_into(
                q,
                velocity,
                dq,
                dv,
                reservoir_directions_pa,
                out,
                check,
            );
        }
        let original = self.film.power_law_film();
        let geometry = original.pressure_operator();
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
        let reservoirs = self.film.evaluation_reservoirs();
        let state = original.solve_with_motion(
            q,
            velocity,
            reservoirs,
            Some(&self.coefficients),
            &mut check,
        )?;
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
            let coefficients = &self.coefficients[k];
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
            let coefficients = &self.coefficients[k];
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

#[cfg(test)]
mod tests;
