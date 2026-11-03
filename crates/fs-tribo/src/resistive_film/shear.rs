//! Sliding no-slip walls on the original incompressible Reynolds pressure graph.
//!
//! For normalized gap fluidity A = integral(1/mu dz), z0 = integral(z/mu dz)/A,
//! entrainment is Qc = width*h*(z0*Ulower + (1-z0)*Uupper), and relative-wall
//! resistance is S = width*length/(h*A). The native graph solves
//! L*psi = swept_volume_rate - div(Qc). A common Barus multiplier
//! mu(z,p)=mu_0(z)*exp(alpha*p) leaves Qc unchanged and gives
//! Qp=G*delta_psi. Its exact length-averaged shear multiplier is
//! M=delta_p/delta_psi, with M=exp(alpha*p) at equal endpoint pressure.
//! Thus S=width*length*M/(h*A); at alpha=0, psi=p and M=1 exactly.
//!
//! Normal pressure force alone does not account for a moving wall's work.
//! Positive RESISTING wall forces along the oriented channel are
//! F_lower = -width*h*z0*delta_p - S*(Uupper-Ulower),
//! F_upper = -width*h*(1-z0)*delta_p + S*(Uupper-Ulower).
//! Projecting these through the SAME velocity rows gives
//! F dot v = sum(Qp*delta_p + S*(Uupper-Ulower)^2) >= 0.
//! The pressure-dependent pump contribution is not counted as extra heat.
//! This is the stationary-profile Couette/Poiseuille specialization of the
//! generalized Reynolds reduction (Dowson 1962,
//! <https://doi.org/10.1016/S0020-7403(62)80038-1>).
//!
//! Profiles and linear wall rows are prescribed and held fixed within a solve.
//! The model supplies no fluid inertia, thermal/advection time owner, layer
//! inventory, interface motion, slip or dry-contact friction. Barus support
//! remains the native classical Reynolds reduction: a prescribed common alpha
//! and admitted pressure interval, not a general high-pressure or EHL claim.
//! Cross-gap pressure corrections are outside this model; see Gustafsson et al.
//! (2015), <https://doi.org/10.1016/j.apm.2015.03.028>.

use fs_math::det;

mod layer_work;
use layer_work::LayerWork;
pub use layer_work::MAX_LAYER_REPORTS;

use super::stratified::{GapViscosityProfile, StratifiedFilm};
use super::{
    FilmError, FilmReport, MAX_CELLS, MAX_CHANNELS, MAX_PORTS, ResistiveFilm, Solution, dot, fail,
    solve_factored,
};

/// Tangential wall velocities in one original channel's from-to orientation.
///
/// Each row maps the SAME work-conjugate mechanical rate passed to the film
/// into a physical wall speed [m/s]. Its transpose maps resisting wall force
/// [N] back to that mechanical port. Supply zero rows for fixed walls. Reverse
/// row signs when the same physical wall runs against a channel's orientation.
#[derive(Debug, Clone, PartialEq)]
pub struct WallVelocityPorts {
    pub lower: Vec<f64>,
    pub upper: Vec<f64>,
}

/// Flow and irreversible work from the simultaneous pressure/wall solve.
///
/// Total flow times pressure drop is not the irreversible power of a sliding
/// passage: a moving wall can pump fluid or oppose pressure-driven flow.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ShearChannelReport {
    /// Signed TOTAL volume flow, pressure-driven plus entrained [m^3/s].
    pub volume_flow_m3_s: f64,
    /// Signed wall-driven entrainment in the same orientation [m^3/s].
    pub entrained_flow_m3_s: f64,
    /// Nonnegative pressure-driven Qp*delta_p [W].
    pub pressure_dissipated_power_w: f64,
    /// Nonnegative S*(Uupper-Ulower)^2 [W].
    pub shear_dissipated_power_w: f64,
}

impl ShearChannelReport {
    /// Total physical viscous power in this passage [W], counted once.
    ///
    /// In a stratified passage the pressure/Couette shear-stress cross terms
    /// cancel over the complete gap, but generally not within each layer.
    /// Pressure-only layer fractions must not be applied to this combined heat;
    /// use [ShearedFilm::evaluate_with_layers_into] for its local distribution.
    #[must_use]
    pub fn dissipated_power_w(&self) -> f64 {
        self.pressure_dissipated_power_w + self.shear_dissipated_power_w
    }
}

#[derive(Debug, Clone)]
struct Coefficients {
    mean_velocity: Vec<f64>,
    relative_velocity: Vec<f64>,
    shear_scale: f64,
    moving: bool,
}

/// Native pressure elimination with complete work-conjugate sliding-wall force.
///
/// The original graph and profiles are retained. Force and tangent evaluations
/// allocate no heap storage and use the existing bounded pressure factorization.
/// At fixed geometry the constant-viscosity resistance is reciprocal. Barus
/// pressure dependence keeps nonnegative work but generally has a nonsymmetric
/// tangent. A normal-pressure force without the wall tractions is not the
/// operator represented here.
///
/// Mean and relative wall rows are compiled before evaluating velocity. Small
/// sliding coordinates therefore survive a large shared translation, and
/// pressure/shear tractions use those same rows for reciprocal projection.
/// Unrepresentable derived rows refuse at construction.
#[derive(Debug, Clone)]
pub struct ShearedFilm {
    film: ResistiveFilm,
    profiles: Option<Vec<GapViscosityProfile>>,
    walls: Vec<WallVelocityPorts>,
    coefficients: Vec<Coefficients>,
    layer_work: LayerWork,
}

impl ShearedFilm {
    /// Bind walls when each native channel viscosity is UNIFORM across its gap.
    ///
    /// Channel viscosities may differ. A hydraulic equivalent of an unresolved
    /// profile is insufficient for shear drag; use from_stratified to retain
    /// that profile and its distinct harmonic shear viscosity.
    pub fn new(film: ResistiveFilm, walls: Vec<WallVelocityPorts>) -> Result<Self, FilmError> {
        Self::compile(film, None, walls)
    }

    /// Bind walls to the exact retained layer profiles and their pressure owner.
    ///
    /// Both entrainment weights are integrated from positive distances to
    /// their respective walls. The upper weight is not obtained by subtracting
    /// a rounded near-wall neutral plane from one. Shear viscosity is 1/A,
    /// independently of the hydraulic-equivalent pressure coefficient.
    pub fn from_stratified(
        stratified: StratifiedFilm,
        walls: Vec<WallVelocityPorts>,
    ) -> Result<Self, FilmError> {
        let (film, profiles) = stratified.into_parts();
        Self::compile(film, Some(profiles), walls)
    }

    fn compile(
        film: ResistiveFilm,
        profiles: Option<Vec<GapViscosityProfile>>,
        walls: Vec<WallVelocityPorts>,
    ) -> Result<Self, FilmError> {
        if walls.len() != film.channel_count()
            || walls.iter().any(|wall| {
                wall.lower.len() != film.port_count()
                    || wall.upper.len() != film.port_count()
                    || wall.lower.iter().chain(&wall.upper).any(|x| !x.is_finite())
            })
        {
            return fail("sliding film needs finite wall rows for every original channel and port");
        }
        let layer_work = LayerWork::compile(film.channel_count(), profiles.as_deref())?;
        let mut coefficients = Vec::with_capacity(walls.len());
        for (i, (edge, wall)) in film.channels.iter().zip(&walls).enumerate() {
            let (lower_fraction, upper_fraction, shear_viscosity) = match &profiles {
                Some(profiles) => profile_coefficients(&profiles[i])?,
                None => (0.5, 0.5, film.channel_viscosities_pa_s()[i]),
            };
            let mut mean_velocity = Vec::with_capacity(film.ports);
            let mut relative_velocity = Vec::with_capacity(film.ports);
            let combined_weight = lower_fraction + upper_fraction;
            for (&lower, &upper) in wall.lower.iter().zip(&wall.upper) {
                let relative = upper - lower;
                if !relative.is_finite() {
                    return fail("relative wall velocity row is not representable");
                }
                // A fused interpolation from the smaller-magnitude endpoint
                // avoids losing equal subnormal rows or overflowing a finite
                // weighted mean through a larger intermediate endpoint.
                let mean = if lower.abs() <= upper.abs() {
                    upper_fraction.mul_add(relative, combined_weight * lower)
                } else {
                    (-lower_fraction).mul_add(relative, combined_weight * upper)
                };
                if !mean.is_finite() {
                    return fail("mean wall velocity row is not representable");
                }
                mean_velocity.push(mean);
                relative_velocity.push(relative);
            }
            let moving = wall.lower.iter().chain(&wall.upper).any(|&x| x != 0.0);
            let shear_scale = if moving {
                // Pair the smallest/largest positive factors first so a finite
                // product is not lost to an avoidable intermediate overflow.
                let mut terms = [edge.width_m, edge.length_m, shear_viscosity];
                terms.sort_by(f64::total_cmp);
                let scale = (terms[0] * terms[2]) * terms[1];
                if !scale.is_finite() || scale <= 0.0 {
                    return fail("sliding film shear coefficient is not representable");
                }
                scale
            } else {
                0.0
            };
            coefficients.push(Coefficients {
                mean_velocity,
                relative_velocity,
                shear_scale,
                moving,
            });
        }
        Ok(Self {
            film,
            profiles,
            walls,
            coefficients,
            layer_work,
        })
    }

    /// Original channel wall rows; no velocity is inferred from a gap closure.
    #[must_use]
    pub fn wall_ports(&self) -> &[WallVelocityPorts] {
        &self.walls
    }

    /// Exact original layers when constructed through from_stratified.
    #[must_use]
    pub fn stratified_profiles(&self) -> Option<&[GapViscosityProfile]> {
        self.profiles.as_deref()
    }

    /// Number of layer powers returned by [Self::evaluate_with_layers_into].
    ///
    /// A uniform channel has one implicit layer. Explicit profiles retain every
    /// authored layer in lower-to-upper order; shut channels retain their slots.
    #[must_use]
    pub fn layer_count(&self) -> usize {
        self.layer_work.len()
    }

    /// Flat output range for one original channel, or None outside the graph.
    #[must_use]
    pub fn channel_layer_range(&self, channel: usize) -> Option<core::ops::Range<usize>> {
        self.layer_work.range(channel)
    }

    /// Retained native geometry and pressure coefficients.
    ///
    /// This pressure-only operator omits moving-wall entrainment and tractions.
    /// Use this binding's evaluate/tangent methods for the complete mechanics.
    #[must_use]
    pub const fn pressure_operator(&self) -> &ResistiveFilm {
        &self.film
    }

    /// Combined resisting force and gauge pressure, evaluated simultaneously.
    ///
    /// A channel carrying any nonzero wall row requires a strictly positive
    /// hydraulic gap. A collapsed sliding passage refuses; no gap floor, dry
    /// contact law or clipped wall drag is substituted. Entirely fixed-wall
    /// passages retain native exact closure. Every output is atomic on refusal.
    pub fn evaluate_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
    ) -> Result<FilmReport, FilmError> {
        self.evaluate::<false, false>(q, velocity, force, pressure, &mut [], &mut [])
    }

    /// Return total oriented flow and physical channel work in the SAME solve.
    /// Reports stay in original channel order; no separate heat solve is run.
    pub fn evaluate_with_channels_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [ShearChannelReport],
    ) -> Result<FilmReport, FilmError> {
        self.evaluate::<true, false>(q, velocity, force, pressure, channels, &mut [])
    }

    /// Return physical work in every retained layer from the SAME native solve.
    ///
    /// Layer power [W] integrates the complete pressure/Couette shear-stress
    /// square, including layerwise cross terms. It is nonnegative and closes
    /// against the original per-channel dissipation to floating-point roundoff.
    /// No normalization, inferred wall heat split, mass or thermal clock is
    /// introduced. Use [Self::channel_layer_range] for channel-local destinations.
    ///
    /// All force, pressure, channel and layer outputs remain unchanged on any
    /// refusal, including an unrepresentably small positive layer contribution.
    pub fn evaluate_with_layers_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [ShearChannelReport],
        layer_power_w: &mut [f64],
    ) -> Result<FilmReport, FilmError> {
        self.evaluate::<true, true>(q, velocity, force, pressure, channels, layer_power_w)
    }

    fn evaluate<const REPORT_CHANNELS: bool, const REPORT_LAYERS: bool>(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [ShearChannelReport],
        layer_power_w: &mut [f64],
    ) -> Result<FilmReport, FilmError> {
        if force.len() != self.film.ports
            || pressure.len() != self.film.cell_count()
            || (REPORT_CHANNELS && channels.len() != self.film.channel_count())
            || (REPORT_LAYERS && layer_power_w.len() != self.layer_count())
        {
            return fail("invalid sliding film output shape");
        }
        let state = self.state(q, velocity)?;
        let mut candidate = [0.0; MAX_PORTS];
        let mut channel_candidate = [ShearChannelReport::default(); MAX_CHANNELS];
        self.film.project(&state.pressure.pressure, &mut candidate);
        let mut pressure_drops = [0.0; MAX_CHANNELS];
        let mut power = 0.0;
        for (k, edge) in self.film.channels.iter().enumerate() {
            let local = state.channels[k];
            let coefficient = &self.coefficients[k];
            let dp = state.pressure.pressure[edge.from]
                - edge.to.map_or(0.0, |j| state.pressure.pressure[j]);
            let potential_drop = state.pressure.potential[edge.from]
                - edge.to.map_or(0.0, |j| state.pressure.potential[j]);
            let pressure_flow = state.pressure.conductance[k] * potential_drop;
            let pressure_loss = pressure_flow * dp;
            let shear_force = local.drag * local.relative_speed;
            let shear_loss = shear_force * local.relative_speed;
            let flow = pressure_flow + local.entrainment;
            if coefficient.moving {
                let pressure_force = (edge.width_m * local.gap) * dp;
                for (j, output) in candidate[..self.film.ports].iter_mut().enumerate() {
                    *output += -coefficient.mean_velocity[j] * pressure_force
                        + coefficient.relative_velocity[j] * shear_force;
                }
            }
            if !flow.is_finite()
                || !pressure_loss.is_finite()
                || pressure_loss < 0.0
                || !shear_loss.is_finite()
                || shear_loss < 0.0
            {
                return fail("nonfinite or negative sliding film channel work");
            }
            power += pressure_loss + shear_loss;
            if REPORT_LAYERS {
                pressure_drops[k] = dp;
            }
            if REPORT_CHANNELS || REPORT_LAYERS {
                channel_candidate[k] = ShearChannelReport {
                    volume_flow_m3_s: flow,
                    entrained_flow_m3_s: local.entrainment,
                    pressure_dissipated_power_w: pressure_loss,
                    shear_dissipated_power_w: shear_loss,
                };
            }
        }
        let work = dot(velocity, &candidate[..self.film.ports]);
        if candidate[..self.film.ports].iter().any(|x| !x.is_finite())
            || !power.is_finite()
            || power < 0.0
            || !work.is_finite()
            || (work - power).abs() > 1e-9 * (work.abs() + power.abs()).max(1e-30)
        {
            return fail("sliding film failed reciprocal pressure and wall work balance");
        }
        if REPORT_LAYERS {
            self.layer_work.evaluate_into(
                &state,
                &channel_candidate[..self.film.channel_count()],
                &pressure_drops[..self.film.channel_count()],
                layer_power_w,
            )?;
        }
        force.copy_from_slice(&candidate[..self.film.ports]);
        pressure.copy_from_slice(&state.pressure.pressure[..self.film.cell_count()]);
        if REPORT_CHANNELS {
            channels.copy_from_slice(&channel_candidate[..self.film.channel_count()]);
        }
        Ok(FilmReport {
            minimum_gap_m: state.pressure.minimum_gap,
            peak_pressure_pa: pressure.iter().fold(0.0_f64, |m, p| m.max(p.abs())),
            dissipated_power_w: power,
        })
    }

    /// Exact directional derivative with independent configuration/rate directions.
    ///
    /// The native factorization is reused for
    /// L*dpsi = A*B*dv - (dL)*psi - div(dQc). Both wall pressure tractions and
    /// the M/h shear resistance are differentiated before their transpose
    /// projection, including the Barus inverse and exact viscosity average.
    /// Profiles and wall rows are fixed data during differentiation.
    /// Admission and output publication have the same atomic semantics as force.
    pub fn tangent_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        dq: &[f64],
        dv: &[f64],
        out: &mut [f64],
    ) -> Result<(), FilmError> {
        if dq.len() != self.film.ports
            || dv.len() != self.film.ports
            || out.len() != self.film.ports
            || dq.iter().chain(dv).any(|x| !x.is_finite())
        {
            return fail("invalid sliding film tangent shape or direction");
        }
        let state = self.state(q, velocity)?;
        let mut rhs = [0.0; MAX_CELLS];
        for (i, cell) in self.film.cells.iter().enumerate() {
            rhs[i] = cell.area_m2 * dot(&cell.gap.closure, dv);
        }
        for (k, (edge, factor)) in self
            .film
            .channels
            .iter()
            .zip(&self.film.factors)
            .enumerate()
        {
            let local = state.channels[k];
            let coefficient = &self.coefficients[k];
            let dh = -dot(&edge.gap.closure, dq);
            let dg = if local.gap > 0.0 {
                3.0 * factor * local.gap * local.gap * dh
            } else {
                0.0
            };
            let dpsi = state.pressure.potential[edge.from]
                - edge.to.map_or(0.0, |j| state.pressure.potential[j]);
            let entrainment_change = if coefficient.moving {
                let mean_change = dot(&coefficient.mean_velocity, dv);
                edge.width_m * (dh * local.mean_speed + local.gap * mean_change)
            } else {
                0.0
            };
            let flux_change = dg * dpsi + entrainment_change;
            rhs[edge.from] -= flux_change;
            if let Some(j) = edge.to {
                rhs[j] += flux_change;
            }
        }
        solve_factored(
            &state.pressure.lower,
            &state.pressure.scale,
            self.film.cell_count(),
            &mut rhs,
        )?;
        let potential_change = rhs;
        let alpha = self.film.pressure_viscosity_coefficient_per_pa();
        if alpha != 0.0 {
            for (dp, psi) in rhs[..self.film.cell_count()]
                .iter_mut()
                .zip(&state.pressure.potential)
            {
                *dp /= 1.0 - alpha * psi;
            }
        }
        let mut candidate = [0.0; MAX_PORTS];
        self.film.project(&rhs, &mut candidate);
        for (k, edge) in self.film.channels.iter().enumerate() {
            let coefficient = &self.coefficients[k];
            if !coefficient.moving {
                continue;
            }
            let local = state.channels[k];
            let dh = -dot(&edge.gap.closure, dq);
            let dp = state.pressure.pressure[edge.from]
                - edge.to.map_or(0.0, |j| state.pressure.pressure[j]);
            let pressure_change = rhs[edge.from] - edge.to.map_or(0.0, |j| rhs[j]);
            let relative_change = dot(&coefficient.relative_velocity, dv);
            let mut shear_change = (-local.drag * (dh / local.gap)) * local.relative_speed
                + local.drag * relative_change;
            if alpha != 0.0 {
                let (_, multiplier_change) = barus_mean(
                    alpha,
                    state.pressure.potential[edge.from],
                    edge.to.map_or(0.0, |j| state.pressure.potential[j]),
                    potential_change[edge.from],
                    edge.to.map_or(0.0, |j| potential_change[j]),
                )?;
                shear_change += (coefficient.shear_scale / local.gap)
                    * multiplier_change
                    * local.relative_speed;
            }
            let pressure_force_change = edge.width_m * (dh * dp + local.gap * pressure_change);
            for (j, output) in candidate[..self.film.ports].iter_mut().enumerate() {
                *output += -coefficient.mean_velocity[j] * pressure_force_change
                    + coefficient.relative_velocity[j] * shear_change;
            }
        }
        if candidate[..self.film.ports].iter().any(|x| !x.is_finite()) {
            return fail("nonfinite sliding film force tangent");
        }
        out.copy_from_slice(&candidate[..self.film.ports]);
        Ok(())
    }

    fn state(&self, q: &[f64], velocity: &[f64]) -> Result<State, FilmError> {
        if q.len() != self.film.ports
            || velocity.len() != self.film.ports
            || q.iter().chain(velocity).any(|x| !x.is_finite())
        {
            return fail("invalid sliding film configuration or effort");
        }
        let mut channels = [ChannelState::default(); MAX_CHANNELS];
        let mut entrainment = [0.0; MAX_CHANNELS];
        for (k, edge) in self.film.channels.iter().enumerate() {
            let coefficient = &self.coefficients[k];
            let h = edge.gap.reference_m - dot(&edge.gap.closure, q);
            if !h.is_finite() || (coefficient.moving && h <= 0.0) {
                return fail("sliding passage needs a finite positive gap");
            }
            if h <= 0.0 {
                continue;
            }
            let (mean_speed, relative_speed, drag) = if coefficient.moving {
                let mean = dot(&coefficient.mean_velocity, velocity);
                let relative = dot(&coefficient.relative_velocity, velocity);
                let drag = coefficient.shear_scale / h;
                if !mean.is_finite() || !relative.is_finite() || !drag.is_finite() || drag <= 0.0 {
                    return fail("sliding velocity or shear resistance is not representable");
                }
                (mean, relative, drag)
            } else {
                (0.0, 0.0, 0.0)
            };
            let flow = if coefficient.moving {
                (edge.width_m * h) * mean_speed
            } else {
                0.0
            };
            if !flow.is_finite() {
                return fail("sliding entrainment is not representable");
            }
            entrainment[k] = flow;
            channels[k] = ChannelState {
                gap: h,
                mean_speed,
                relative_speed,
                drag,
                entrainment: flow,
            };
        }
        let pressure = self.film.solve_with_entrainment::<true>(
            q,
            velocity,
            &entrainment[..self.film.channel_count()],
        )?;
        let alpha = self.film.pressure_viscosity_coefficient_per_pa();
        if alpha != 0.0 {
            for (k, edge) in self.film.channels.iter().enumerate() {
                if self.coefficients[k].moving {
                    let (multiplier, _) = barus_mean(
                        alpha,
                        pressure.potential[edge.from],
                        edge.to.map_or(0.0, |j| pressure.potential[j]),
                        0.0,
                        0.0,
                    )?;
                    channels[k].drag *= multiplier;
                    if !channels[k].drag.is_finite() || channels[k].drag <= 0.0 {
                        return fail("pressure-dependent sliding resistance is not representable");
                    }
                }
            }
        }
        Ok(State { pressure, channels })
    }
}

// Integrate exp(alpha*p(x)) along the passage, where psi(x) is affine.
// The coincident-pressure divided difference and its derivative are evaluated
// by a fixed convergent series; no pressure subtraction or quadrature is used.
fn barus_mean(
    alpha: f64,
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
    let r = 1.0 - alpha * high;
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
        let r_low = 1.0 - alpha * low;
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

#[derive(Clone, Copy, Default)]
struct ChannelState {
    gap: f64,
    mean_speed: f64,
    relative_speed: f64,
    drag: f64,
    entrainment: f64,
}

struct State {
    pressure: Solution,
    channels: [ChannelState; MAX_CHANNELS],
}

fn profile_coefficients(profile: &GapViscosityProfile) -> Result<(f64, f64, f64), FilmError> {
    let layers = profile.layers();
    if layers
        .iter()
        .all(|layer| layer.dynamic_viscosity_pa_s == layers[0].dynamic_viscosity_pa_s)
    {
        return Ok((0.5, 0.5, layers[0].dynamic_viscosity_pa_s));
    }
    let minimum = layers
        .iter()
        .map(|layer| layer.dynamic_viscosity_pa_s)
        .fold(f64::INFINITY, f64::min);
    let mut total = Sum::default();
    let mut lower = Sum::default();
    let mut upper = Sum::default();
    let mut previous = 0.0;
    for layer in layers {
        let width = layer.upper_fraction - previous;
        let weight = width * (minimum / layer.dynamic_viscosity_pa_s);
        total.add(weight);
        lower.add(weight * (previous + 0.5 * width));
        upper.add(weight * ((1.0 - layer.upper_fraction) + 0.5 * width));
        previous = layer.upper_fraction;
    }
    let sum = total.value();
    let fractions = [lower.value() / sum, upper.value() / sum];
    let shear_viscosity = minimum / sum;
    if fractions
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0 || *x > 1.0)
        || !shear_viscosity.is_finite()
        || shear_viscosity <= 0.0
    {
        return fail("stratified sliding coefficients are not representable");
    }
    Ok((fractions[0], fractions[1], shear_viscosity))
}

#[derive(Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }
    fn value(self) -> f64 {
        self.sum + self.correction
    }
}

#[cfg(test)]
mod tests;
