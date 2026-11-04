//! Prescribed reservoir pressures on the original Reynolds-film graph.
//!
//! Externally pressurized lubrication supports load at zero wall speed. The
//! original elimination now solves `L psi = A B v + sum_boundary G psi_b`.
//! Its force is still `B^T A p`, but it is no longer a passive resistance:
//! `F dot v - sum_boundary p_b Q_out = sum_channels Q delta_p >= 0`.
//! Reservoir supply must therefore enter the surrounding mechanical/thermal
//! account; neither `F dot v` nor reservoir power alone is the viscous heat.
//!
//! Pressures share the original gauge datum, also the viscosity reference and
//! opposing mechanical pressure. This is the classical incompressible Reynolds
//! reduction (NASA, *Fundamentals of Fluid Film Lubrication*, chapter 14,
//! <https://ntrs.nasa.gov/citations/19910021217>). It supplies no reservoir
//! inventory, supply-line dynamics, cavitation, fluid inertia or thermal owner.

use fs_math::det;

use super::{
    FilmChannelReport, FilmError, FilmReport, MAX_CELLS, MAX_CHANNELS, MAX_PORTS, ResistiveFilm,
    dot, fail, solve_factored,
};

#[derive(Clone, Copy, Debug)]
pub(super) struct PressureDatum {
    pub(super) pressure_pa: f64,
    /// Original zero-gauge fluidity multiplier exp(-alpha*p_ref).
    pub(super) fluidity: f64,
}

impl PressureDatum {
    pub(super) const ZERO: Self = Self {
        pressure_pa: 0.0,
        fluidity: 1.0,
    };
}

/// Native film with prescribed pressure at EVERY original reservoir channel.
///
/// `Some(p)` is required on a boundary (`to == None`); internal channels require
/// `None`. The mapping retains original channel order, including shut passages.
/// All pressures are gauge Pa relative to the SAME external mechanical datum.
/// A uniform pressure offset changes the physical normal load. With a Barus
/// law it also changes viscosity, which remains referenced to zero gauge.
///
/// Construction compiles pressure potentials once. Evaluation and the full
/// force tangent use the native factorization and bounded stack storage.
#[derive(Clone, Debug)]
pub struct PressureDrivenFilm {
    film: ResistiveFilm,
    reservoirs: ReservoirPressures,
    ambient: bool,
}

/// One compiled pressure boundary shared by normal and sliding-wall films.
#[derive(Clone, Debug)]
pub(super) struct ReservoirPressures {
    pub(super) reservoirs_pa: Vec<Option<f64>>,
    pub(super) potential: Vec<f64>,
    pub(super) datum: PressureDatum,
}

/// Work from one simultaneous pressure/flow/force evaluation, in watts.
/// These are instantaneous powers, not a time-integrated passivity receipt.
#[derive(Clone, Copy, Debug)]
pub struct PressureDrivenReport {
    /// Native geometry, actual pressure bound and nonnegative viscous heat.
    pub film: FilmReport,
    /// `F dot v`, positive INTO the film from mechanical motion. May be negative
    /// when reservoir pressure drives the bodies apart.
    pub mechanical_power_w: f64,
    /// `-sum p_b Q_out`, positive INTO the film from pressure supplies.
    pub reservoir_power_w: f64,
    /// Mechanical plus reservoir input minus viscous heat.
    pub power_residual_w: f64,
}

impl ReservoirPressures {
    pub(super) fn new(
        film: &ResistiveFilm,
        reservoirs_pa: &[Option<f64>],
    ) -> Result<Self, FilmError> {
        if reservoirs_pa.len() != film.channel_count() {
            return fail("film reservoir pressure count must match original channels");
        }
        for (edge, pressure) in film.channels.iter().zip(reservoirs_pa) {
            if edge.to.is_none() != pressure.is_some()
                || pressure
                    .is_some_and(|p| !p.is_finite() || p.abs() > film.limits.maximum_pressure_pa)
            {
                return fail("film reservoir pressure topology or gauge domain");
            }
        }
        let pressure_pa = reservoirs_pa
            .iter()
            .flatten()
            .copied()
            .next()
            .unwrap_or(0.0);
        let alpha = film.pressure_viscosity_coefficient_per_pa;
        let fluidity = if alpha == 0.0 {
            1.0
        } else {
            (-alpha * pressure_pa).exp()
        };
        if !fluidity.is_finite() || fluidity <= 0.0 {
            return fail("unrepresentable film reservoir fluidity");
        }
        let datum = PressureDatum {
            pressure_pa,
            fluidity,
        };
        let mut potential = Vec::with_capacity(reservoirs_pa.len());
        for (k, pressure) in reservoirs_pa.iter().enumerate() {
            let value = if let Some(p) = pressure {
                let difference = p - pressure_pa;
                let scaled = alpha * difference;
                let value = if scaled == 0.0 {
                    fluidity * difference
                } else {
                    fluidity * difference * (-(-scaled).exp_m1() / scaled)
                };
                let endpoint_fluidity = fluidity - alpha * value;
                let viscosity = film.channel_viscosities_pa_s[k] / endpoint_fluidity;
                if !value.is_finite() || !viscosity.is_finite() || viscosity <= 0.0 {
                    return fail("unrepresentable film reservoir pressure or viscosity");
                }
                value
            } else {
                0.0
            };
            potential.push(value);
        }
        Ok(Self {
            reservoirs_pa: reservoirs_pa.to_vec(),
            potential,
            datum,
        })
    }
}

impl PressureDrivenFilm {
    /// Bind fixed reservoir pressures without changing geometry or transport.
    /// The same pressure limit and Barus source law apply to reservoir and cell
    /// endpoints. No omitted or duplicate reservoir can silently become ambient.
    pub fn new(film: ResistiveFilm, reservoirs_pa: &[Option<f64>]) -> Result<Self, FilmError> {
        let reservoirs = ReservoirPressures::new(&film, reservoirs_pa)?;
        let ambient = reservoirs_pa.iter().flatten().all(|p| *p == 0.0);
        Ok(Self { film, reservoirs, ambient })
    }

    /// Original channel-indexed boundary mapping, in gauge Pa.
    #[must_use]
    pub fn reservoir_pressures_pa(&self) -> &[Option<f64>] {
        &self.reservoirs.reservoirs_pa
    }

    /// Rebind the prescribed supplies while retaining the original graph and
    /// viscosity law. This stateless owner has no pressure history to reset.
    pub fn with_reservoir_pressures(
        &self,
        pressures_pa: &[Option<f64>],
    ) -> Result<Self, FilmError> {
        Self::new(self.film.clone(), pressures_pa)
    }

    #[must_use]
    pub fn port_count(&self) -> usize {
        self.film.port_count()
    }
    #[must_use]
    pub fn cell_count(&self) -> usize {
        self.film.cell_count()
    }
    #[must_use]
    pub fn channel_count(&self) -> usize {
        self.film.channel_count()
    }

    /// Solve the original film and publish force, physical cell gauge pressures,
    /// passage flows/heat, and complete mechanical/reservoir power together.
    /// The force has the original RESISTING sign; the body's force is its negative.
    /// A failed solve or work check leaves all output slices untouched.
    pub fn evaluate_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [FilmChannelReport],
    ) -> Result<PressureDrivenReport, FilmError> {
        let film = &self.film;
        if force.len() != film.ports
            || pressure.len() != film.cell_count()
            || channels.len() != film.channel_count()
        {
            return fail("invalid pressure-driven film output shape");
        }
        if self.ambient {
            // Preserve the original ambient pressure, channel work and force
            // arithmetic exactly, including internal-channel cancellation.
            let report = film.evaluate_with_channels_into(q, velocity, force, pressure, channels)?;
            let work = dot(velocity, force);
            return Ok(PressureDrivenReport {
                film: report,
                mechanical_power_w: work,
                reservoir_power_w: 0.0,
                power_residual_w: work - report.dissipated_power_w,
            });
        }
        let s = film.solve_with_boundaries::<false, true>(
            q,
            velocity,
            &[],
            &self.reservoirs.potential,
            self.reservoirs.datum,
        )?;
        let mut candidate = [0.0; MAX_PORTS];
        let mut reports = [FilmChannelReport::default(); MAX_CHANNELS];
        film.project(&s.pressure, &mut candidate);
        let alpha = film.pressure_viscosity_coefficient_per_pa;
        let mut loss = 0.0;
        let mut supply = 0.0;
        let mut supply_magnitude = 0.0;
        for (k, edge) in film.channels.iter().enumerate() {
            let to = edge.to.map_or(self.reservoirs.potential[k], |j| s.potential[j]);
            let delta = s.potential[edge.from] - to;
            let flow = s.conductance[k] * delta;
            // Recover a pressure DIFFERENCE without subtracting two biased
            // physical pressures. This is also exact for the constant-mu law.
            let dp = pressure_drop_at_datum(
                alpha,
                self.reservoirs.datum.fluidity,
                s.potential[edge.from],
                to,
            )?;
            let power = flow * dp;
            if !flow.is_finite() || !power.is_finite() || power < 0.0 {
                return fail("nonfinite or negative pressure-driven film channel work");
            }
            reports[k] = FilmChannelReport {
                volume_flow_m3_s: flow,
                dissipated_power_w: power,
            };
            loss += power;
            if let Some(p) = self.reservoirs.reservoirs_pa[k] {
                let incoming = -p * flow;
                supply += incoming;
                supply_magnitude += incoming.abs();
            }
        }
        let work = dot(velocity, &candidate[..film.ports]);
        let residual = work + supply - loss;
        let magnitude = work.abs() + supply_magnitude + loss;
        if candidate[..film.ports].iter().any(|f| !f.is_finite())
            || !magnitude.is_finite()
            || !residual.is_finite()
            || residual.abs() > 1e-9 * magnitude.max(1e-30)
        {
            return fail("pressure-driven film failed mechanical/reservoir work balance");
        }
        force.copy_from_slice(&candidate[..film.ports]);
        pressure.copy_from_slice(&s.pressure[..film.cell_count()]);
        channels.copy_from_slice(&reports[..film.channel_count()]);
        Ok(PressureDrivenReport {
            film: FilmReport {
                minimum_gap_m: s.minimum_gap,
                peak_pressure_pa: pressure.iter().fold(0.0_f64, |m, p| m.max(p.abs())),
                dissipated_power_w: loss,
            },
            mechanical_power_w: work,
            reservoir_power_w: supply,
            power_residual_w: residual,
        })
    }

    /// Exact directional force derivative including geometry, normal velocity
    /// AND supply pressures. Pressure directions follow original channel order;
    /// internal entries must be zero. Supply directions do not move the common
    /// mechanical/viscosity datum. Uses one native solve and its same factorization.
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
        let film = &self.film;
        if dq.len() != film.ports
            || dv.len() != film.ports
            || out.len() != film.ports
            || reservoir_directions_pa.len() != film.channel_count()
            || dq
                .iter()
                .chain(dv)
                .chain(reservoir_directions_pa)
                .any(|v| !v.is_finite())
            || film
                .channels
                .iter()
                .zip(reservoir_directions_pa)
                .any(|(e, d)| e.to.is_some() && *d != 0.0)
        {
            return fail("invalid pressure-driven film tangent direction or shape");
        }
        if self.ambient && reservoir_directions_pa.iter().all(|d| *d == 0.0) {
            return film.tangent_into(q, velocity, dq, dv, out);
        }
        let s = film.solve_with_boundaries::<false, true>(
            q,
            velocity,
            &[],
            &self.reservoirs.potential,
            self.reservoirs.datum,
        )?;
        let alpha = film.pressure_viscosity_coefficient_per_pa;
        let mut rhs = [0.0; MAX_CELLS];
        for (i, cell) in film.cells.iter().enumerate() {
            rhs[i] = cell.area_m2 * dot(&cell.gap.closure, dv);
        }
        for (k, (edge, factor)) in film.channels.iter().zip(&film.factors).enumerate() {
            let h = edge.gap.reference_m - dot(&edge.gap.closure, q);
            if h <= 0.0 {
                continue;
            }
            let dh = -dot(&edge.gap.closure, dq);
            let dg = 3.0 * factor * h * h * dh;
            let to = edge.to.map_or(self.reservoirs.potential[k], |j| s.potential[j]);
            let change = dg * (s.potential[edge.from] - to);
            rhs[edge.from] -= change;
            if let Some(j) = edge.to {
                rhs[j] += change;
            } else {
                let dpsi = (self.reservoirs.datum.fluidity - alpha * to) * reservoir_directions_pa[k];
                rhs[edge.from] += s.conductance[k] * dpsi;
            }
        }
        solve_factored(&s.lower, &s.scale, film.cell_count(), &mut rhs)?;
        if alpha != 0.0 {
            for (dp, psi) in rhs[..film.cell_count()].iter_mut().zip(&s.potential) {
                *dp /= self.reservoirs.datum.fluidity - alpha * psi;
            }
        }
        let mut candidate = [0.0; MAX_PORTS];
        film.project(&rhs, &mut candidate);
        if candidate[..film.ports].iter().any(|f| !f.is_finite()) {
            return fail("nonfinite pressure-driven film force tangent");
        }
        out.copy_from_slice(&candidate[..film.ports]);
        Ok(())
    }
}

/// Physical pressure difference without requiring its viscosity multiplier
/// to be representable. Sorted fluidities avoid the false negative-log endpoint
/// at large contrast; zero drop remains exact even at subnormal fluidity.
pub(super) fn pressure_drop_at_datum(
    alpha: f64,
    fluidity: f64,
    from: f64,
    to: f64,
) -> Result<f64, FilmError> {
    let signed_delta = from - to;
    if alpha == 0.0 || signed_delta == 0.0 {
        return if signed_delta.is_finite() {
            Ok(signed_delta)
        } else {
            fail("film pressure difference is not representable")
        };
    }
    let (high, low, sign) = if from > to { (from, to, 1.0) } else { (to, from, -1.0) };
    let r_high = fluidity - alpha * high;
    let r_low = fluidity - alpha * low;
    if !r_high.is_finite() || !r_low.is_finite() || r_high <= 0.0 || r_low <= 0.0 {
        return fail("film pressure difference needs positive finite fluidities");
    }
    let delta = high - low;
    let scaled = alpha * delta / r_high;
    let magnitude = if scaled == 0.0 {
        delta / r_high
    } else if scaled.is_finite() {
        scaled.ln_1p() / alpha
    } else {
        (det::ln(r_low) - det::ln(r_high)) / alpha
    };
    if !magnitude.is_finite() || magnitude < 0.0 {
        return fail("film pressure difference is not representable");
    }
    Ok(sign * magnitude)
}

// Integrate exp(alpha*p(x)) along the passage, where psi(x) is affine.
// The coincident-pressure divided difference and its derivative are evaluated
// by a fixed convergent series; no pressure subtraction or quadrature is used.
pub(super) fn barus_mean_at_datum(
    alpha: f64,
    fluidity: f64,
    from: f64,
    to: f64,
    from_change: f64,
    to_change: f64,
) -> Result<(f64, f64), FilmError> {
    if alpha == 0.0 {
        return Ok((1.0, 0.0));
    }
    let (high, low, high_change, low_change) = if from >= to {
        (from, to, from_change, to_change)
    } else {
        (to, from, to_change, from_change)
    };
    let r = fluidity - alpha * high;
    let delta = high - low;
    let x = alpha * delta / r;
    let (mean, change) = if x <= 0.125 {
        let mut logarithmic = 1.0;
        let mut high_weight = 0.5;
        let mut low_weight = 0.5;
        let mut power = 1.0;
        for n in 1..=20 {
            power *= -x;
            logarithmic += power / f64::from(n + 1);
            high_weight += power / f64::from(n + 2);
            low_weight += power * f64::from(n + 1) / f64::from(n + 2);
        }
        let mean = logarithmic / r;
        let change = if high_change == 0.0 && low_change == 0.0 {
            0.0
        } else {
            (alpha / r) * ((high_weight * high_change + low_weight * low_change) / r)
        };
        (mean, change)
    } else {
        let r_low = fluidity - alpha * low;
        let log_ratio = if x.is_finite() {
            det::ln(1.0 + x)
        } else {
            det::ln(r_low) - det::ln(r)
        };
        let mean = (log_ratio / alpha) / delta;
        let change = if high_change == 0.0 && low_change == 0.0 {
            0.0
        } else {
            ((high_change / r - low_change / r_low) - mean * (high_change - low_change)) / delta
        };
        (mean, change)
    };
    if !mean.is_finite() || mean <= 0.0 || !change.is_finite() {
        return fail("Barus sliding viscosity average or tangent is not representable");
    }
    Ok((mean, change))
}

#[cfg(test)]
mod tests;
