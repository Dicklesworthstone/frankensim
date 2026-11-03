//! Simultaneous incoming and outgoing particle streams on the native body.
//!
//! With held body-relative moments, `M'=qi-qo` and
//! `M*v'=F+Pi-Po+M*g`. Gross streams remain separate even at zero net mass.
//! Gain uses the time reverse of the original loss integrals: its mean log
//! coefficient is `A-C`, and its velocity variance is the SAME loss variance.
//! No incoming particle is assigned an inertial velocity independent of the
//! body: such a source has changing relative moments and needs time refinement.

use super::ExhaustEnvironment;
use super::forced::velocity_moments;
use crate::Vec3;
use crate::rocket::{RocketError, RocketForcing, RocketState, mass_integrals, poll};

/// One gross physical particle stream, held relative to the translating body.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParticleStreamMoments {
    /// Nonnegative incoming OR outgoing mass rate [kg/s].
    pub mass_flow_kg_s: f64,
    /// Particle momentum rate relative to the body [N], in its inertial axes.
    /// This is not a continuum pressure force or effective nozzle thrust.
    pub relative_momentum_rate_n: Vec3,
    /// Complete relative particle energy [W], including relative kinetic and
    /// internal/random energy in one declared datum. A formation datum may
    /// make this signed; mechanics cannot infer thermal admissibility from it.
    pub relative_energy_rate_w: f64,
}

impl ParticleStreamMoments {
    /// A physically closed stream, including all its moments.
    pub const ZERO: Self = Self {
        mass_flow_kg_s: 0.0,
        relative_momentum_rate_n: Vec3::ZERO,
        relative_energy_rate_w: 0.0,
    };

    fn validate(self) -> Result<(), RocketError> {
        if !self.mass_flow_kg_s.is_finite()
            || !self.relative_momentum_rate_n.is_finite()
            || !self.relative_energy_rate_w.is_finite()
        {
            return Err(RocketError::NonFiniteInput("particle stream moments"));
        }
        if self.mass_flow_kg_s < 0.0
            || (self.mass_flow_kg_s == 0.0
                && (self.relative_momentum_rate_n != Vec3::ZERO
                    || self.relative_energy_rate_w != 0.0))
        {
            return Err(RocketError::InvalidInput(
                "particle stream mass and moments",
            ));
        }
        Ok(())
    }
}

/// Separate gross streams; neither is erased by a zero signed mass rate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParticleExchangeForcing {
    /// Gross particles entering the retained body.
    pub incoming: ParticleStreamMoments,
    /// Gross particles leaving the retained body.
    pub outgoing: ParticleStreamMoments,
}

/// One gross stream integrated along the same changing body velocity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParticleStreamIntegral {
    /// Nonnegative gross transferred mass [kg].
    pub mass_kg: f64,
    /// Integrated momentum relative to the translating body [N s].
    pub relative_momentum_n_s: Vec3,
    /// Complete original relative energy, including its source datum [J].
    pub relative_energy_j: f64,
    /// Physical momentum in the body's declared inertial frame [kg m/s].
    pub momentum_world_kg_m_s: Vec3,
    /// Complete energy in that same inertial frame [J].
    pub energy_world_j: f64,
    /// Inertial correction retained apart from the original caloric datum [J].
    pub kinetic_transport_energy_j: f64,
}

impl ParticleStreamIntegral {
    const ZERO: Self = Self {
        mass_kg: 0.0,
        relative_momentum_n_s: Vec3::ZERO,
        relative_energy_j: 0.0,
        momentum_world_kg_m_s: Vec3::ZERO,
        energy_world_j: 0.0,
        kinetic_transport_energy_j: 0.0,
    };
}

/// Atomic native endpoint and complete gross stream and load accounting.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParticleExchangeStep {
    /// Original native state after the complete admitted interval.
    pub state_after: RocketState,
    /// Entire held interval [s], without an inferred burnout/coast split.
    pub duration_s: f64,
    /// Original gross incoming stream, transformed along the native trajectory.
    pub incoming: ParticleStreamIntegral,
    /// Original gross outgoing stream, transformed along the same trajectory.
    pub outgoing: ParticleStreamIntegral,
    /// Zero if neither stream carries mass; otherwise the common time mean.
    pub mean_velocity_world_m_s: Vec3,
    /// Time mean of squared speed, retaining its nonnegative variance [m2/s2].
    /// Zero if neither stream carries mass.
    pub mean_speed_squared_m2_s2: f64,
    /// Actual retained kinetic change, evaluated in the initial-velocity chart [J].
    pub kinetic_energy_change_j: f64,
    /// Required caloric credit: incoming minus outgoing complete relative
    /// energy. Incoming relative kinetic energy thermalizes here ONCE. The
    /// composing material owner must apply this; no second caloric model lives
    /// in this mechanics owner.
    pub retained_relative_energy_j: f64,
    /// Applied-force and changing retained-mass gravity impulse [N s].
    pub external_impulse_n_s: Vec3,
    /// Their complete work along the same body trajectory [J].
    pub external_work_j: f64,
    /// Actual retained inventory change + outgoing mass - incoming mass [kg].
    pub mass_residual_kg: f64,
    /// Actual retained momentum change + outgoing - incoming - external [N s].
    pub momentum_residual_kg_m_s: Vec3,
    /// Actual kinetic change + outgoing/incoming inertial corrections - work.
    /// Caloric reference energies are excluded from this mechanical check.
    pub kinetic_energy_residual_j: f64,
}

fn stream_integral(
    stream: ParticleStreamMoments,
    duration: f64,
    velocity: Vec3,
    mean_increment: Vec3,
    mean_squared_increment: f64,
) -> Result<ParticleStreamIntegral, RocketError> {
    let mass = stream.mass_flow_kg_s * duration;
    if !mass.is_finite() || (stream.mass_flow_kg_s > 0.0 && duration > 0.0 && mass == 0.0) {
        return Err(RocketError::Unrepresentable("gross particle mass"));
    }
    let relative_momentum = stream.relative_momentum_rate_n.scale(duration);
    let relative_energy = stream.relative_energy_rate_w * duration;
    if duration > 0.0
        && [
            (stream.relative_momentum_rate_n.x, relative_momentum.x),
            (stream.relative_momentum_rate_n.y, relative_momentum.y),
            (stream.relative_momentum_rate_n.z, relative_momentum.z),
            (stream.relative_energy_rate_w, relative_energy),
        ]
        .into_iter()
        .any(|(rate, integral)| rate != 0.0 && integral == 0.0)
    {
        return Err(RocketError::Unrepresentable(
            "gross particle moment underflow",
        ));
    }
    let momentum_relative = mean_increment.scale(mass).add(relative_momentum);
    let kinetic_relative =
        relative_momentum.dot(mean_increment) + (0.5 * mass) * mean_squared_increment;
    let kinetic = if mass == 0.0 {
        0.0
    } else {
        (0.5 * mass) * velocity.norm_squared() + velocity.dot(momentum_relative) + kinetic_relative
    };
    let result = ParticleStreamIntegral {
        mass_kg: mass,
        relative_momentum_n_s: relative_momentum,
        relative_energy_j: relative_energy,
        momentum_world_kg_m_s: velocity.scale(mass).add(momentum_relative),
        energy_world_j: relative_energy + kinetic,
        kinetic_transport_energy_j: kinetic,
    };
    if !result.relative_momentum_n_s.is_finite()
        || !result.momentum_world_kg_m_s.is_finite()
        || ![
            result.relative_energy_j,
            result.energy_world_j,
            result.kinetic_transport_energy_j,
        ]
        .into_iter()
        .all(f64::is_finite)
    {
        return Err(RocketError::Unrepresentable("gross particle integral"));
    }
    Ok(result)
}

// Form force*time/mass without first underflowing time/mass. Dividing the
// larger numerator first preserves a normal intermediate whenever possible;
// the other order covers very small masses. Refuse an unresolved increment.
fn normalized_impulse(force: Vec3, duration: f64, mass: f64) -> Result<Vec3, RocketError> {
    let component = |value: f64| {
        if value == 0.0 {
            return Ok(0.0);
        }
        let large = value.abs().max(duration);
        let small = value.abs().min(duration);
        for candidate in [
            (large / mass) * small,
            (small / mass) * large,
            (large * small) / mass,
        ] {
            if candidate.is_finite() && candidate > 0.0 {
                return Ok(candidate.copysign(value));
            }
        }
        Err(RocketError::Unrepresentable("particle force increment"))
    };
    if !force.is_finite() {
        return Err(RocketError::Unrepresentable("combined particle force"));
    }
    Ok(Vec3::new(
        component(force.x)?,
        component(force.y)?,
        component(force.z)?,
    ))
}

impl RocketState {
    /// Advance a complete held interval with simultaneous physical intake and
    /// outflow. The original transferable inventory may grow, shrink or stay
    /// fixed. Insufficient inventory refuses; no intake or burnout switch is
    /// invented. Species withdrawal and source qualification remain with the
    /// composing material owner.
    ///
    /// Exact constant relative fluxes and loads use the original native mass
    /// integrals, with gain obtained by time reversal. Both gross streams are
    /// transformed along the same trajectory. Particle kinetic thermalization
    /// is included in `retained_relative_energy_j`, never as an extra heat term.
    ///
    /// MIT 16.07 L14, eq. (11), covers capture and ejection momentum:
    /// <https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/98d19cb2b9983c14c6907052499e7a98_MIT16_07F09_Lec14.pdf>.
    /// Rotation, COM shifts, continuum pressure work and a stationary ambient
    /// particle distribution are outside this held body-relative model.
    pub fn advance_with_particle_exchange(
        &self,
        forcing: ParticleExchangeForcing,
        environment: ExhaustEnvironment,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ParticleExchangeStep, RocketError> {
        poll(&mut cancelled)?;
        forcing.incoming.validate()?;
        forcing.outgoing.validate()?;
        if !duration_s.is_finite() {
            return Err(RocketError::NonFiniteInput("particle exchange duration"));
        }
        if duration_s < 0.0 {
            return Err(RocketError::InvalidInput("particle exchange duration"));
        }
        // Original native admission owns the state's mass, coordinates and
        // velocity. A zero-duration identity also validates the declared loads.
        self.advance(
            RocketForcing {
                mass_flow_kg_s: 0.0,
                effective_exhaust_velocity_m_s: Vec3::ZERO,
                gravity_m_s2: environment.gravity_m_s2,
                external_force_n: environment.external_force_n,
            },
            0.0,
            &mut cancelled,
        )?;
        if duration_s == 0.0 {
            poll(&mut cancelled)?;
            return Ok(ParticleExchangeStep {
                state_after: *self,
                duration_s,
                incoming: ParticleStreamIntegral::ZERO,
                outgoing: ParticleStreamIntegral::ZERO,
                mean_velocity_world_m_s: Vec3::ZERO,
                mean_speed_squared_m2_s2: 0.0,
                kinetic_energy_change_j: 0.0,
                retained_relative_energy_j: 0.0,
                external_impulse_n_s: Vec3::ZERO,
                external_work_j: 0.0,
                mass_residual_kg: 0.0,
                momentum_residual_kg_m_s: Vec3::ZERO,
                kinetic_energy_residual_j: 0.0,
            });
        }
        let rate = forcing.incoming.mass_flow_kg_s - forcing.outgoing.mass_flow_kg_s;
        let signed_mass = rate * duration_s;
        let mut next = *self;
        next.propellant_mass_kg = rate.mul_add(duration_s, self.propellant_mass_kg);
        if !signed_mass.is_finite() || !next.propellant_mass_kg.is_finite() {
            return Err(RocketError::Unrepresentable("particle exchange inventory"));
        }
        if next.propellant_mass_kg < 0.0 {
            return Err(RocketError::InvalidInput("insufficient particle inventory"));
        }
        if duration_s > 0.0
            && rate != 0.0
            && (signed_mass == 0.0 || next.propellant_mass_kg == self.propellant_mass_kg)
        {
            return Err(RocketError::Unrepresentable("particle inventory change"));
        }
        let m0 = self.dry_mass_kg + self.propellant_mass_kg;
        let m1 = self.dry_mass_kg + next.propellant_mass_kg;
        let gain = rate > 0.0;
        let (high, low) = if gain { (m1, m0) } else { (m0, m1) };
        let x = signed_mass.abs() / high;
        if !high.is_finite() || !x.is_finite() || x > 1.0 || (signed_mass != 0.0 && x == 0.0) {
            return Err(RocketError::Unrepresentable("particle mass integrals"));
        }
        let (_, a, c) = if x > 0.0 {
            mass_integrals(x, high, low)
        } else {
            (0.0, 1.0, 0.5)
        };
        if !a.is_finite() || a <= 0.0 || !c.is_finite() || c <= 0.0 {
            return Err(RocketError::Unrepresentable("particle mass integrals"));
        }
        let force = environment.external_force_n;
        let gravity = environment.gravity_m_s2;
        let particle_force = forcing
            .incoming
            .relative_momentum_rate_n
            .sub(forcing.outgoing.relative_momentum_rate_n);
        let d = normalized_impulse(force.add(particle_force), duration_s, high)?;
        let g = gravity.scale(duration_s);
        let (loss_mean, variance, time_log) = velocity_moments(x, a, c, low / high, d, g)?;
        let mean = if gain && x < 0.125 {
            let mut tail = 0.0_f64;
            for n in (0..=16).rev() {
                tail = tail.mul_add(x, 1.0 / (f64::from(n) + 3.0));
            }
            d.add(g).scale(0.5).add(d.scale(x * tail))
        } else if gain {
            d.scale(a - c).add(g.scale(0.5))
        } else {
            loss_mean
        };
        let mean_squared = mean.norm_squared() + variance;
        let ideal_increment = if x < 0.125 {
            let mut tail = 0.0_f64;
            for n in (0..=16).rev() {
                tail = tail.mul_add(x, 1.0 / (f64::from(n) + 2.0));
            }
            d.add(g).add(d.scale(x * tail))
        } else {
            d.scale(a).add(g)
        };
        next.velocity_m_s = self.velocity_m_s.add(ideal_increment);
        next.position_m = self
            .position_m
            .add(self.velocity_m_s.scale(duration_s))
            .add(mean.scale(duration_s));
        poll(&mut cancelled)?;

        // <(M/high)*dv>-<dv>. Subtract the common constant-mass term
        // analytically; opposing force/gravity can leave only O(x^2) work.
        let weighted_difference = if x < 0.125 {
            let mut tail = 0.0_f64;
            for n in (0..=16).rev() {
                let k = f64::from(n);
                let coefficient = if gain {
                    0.5 / (k + 4.0)
                } else {
                    1.0 / ((k + 2.0) * (k + 4.0))
                };
                tail = tail.mul_add(x, coefficient);
            }
            d.add(g)
                .scale(if gain { 1.0 / 6.0 } else { 1.0 / 3.0 })
                .add(d.scale(x * tail))
                .scale(-x)
        } else if gain {
            d.scale(0.25 - 0.5 * (a - c)).add(g.scale(-x / 6.0))
        } else {
            d.scale(time_log).add(g.scale(1.0 / 3.0)).scale(-x)
        };
        let external_impulse = force.scale(duration_s).add(g.scale(high * (1.0 - 0.5 * x)));
        let work_relative = force.scale(duration_s).add(g.scale(high)).dot(mean)
            + g.dot(weighted_difference.scale(high));
        let v0 = self.velocity_m_s;
        let incoming = stream_integral(forcing.incoming, duration_s, v0, mean, mean_squared)?;
        let outgoing = stream_integral(forcing.outgoing, duration_s, v0, mean, mean_squared)?;
        let delta = next.velocity_m_s.sub(v0);
        let actual_mass = next.propellant_mass_kg - self.propellant_mass_kg;
        let gross_loss = outgoing.mass_kg - incoming.mass_kg;
        let mass_residual = actual_mass + gross_loss;
        // Retain the separately stored dry and transferable masses here: their
        // sum used in the analytic trajectory can round at either endpoint.
        let retained_increment = delta
            .scale(self.dry_mass_kg)
            .add(delta.scale(next.propellant_mass_kg));
        let relative_outflow = outgoing
            .relative_momentum_n_s
            .sub(incoming.relative_momentum_n_s);
        let momentum_relative = retained_increment
            .add(mean.scale(gross_loss))
            .add(relative_outflow)
            .sub(external_impulse);
        let kinetic_relative = (0.5 * self.dry_mass_kg) * delta.norm_squared()
            + (0.5 * next.propellant_mass_kg) * delta.norm_squared();
        let transport_relative = relative_outflow.dot(mean) + (0.5 * gross_loss) * mean_squared;
        let initial_speed_squared = if actual_mass == 0.0 && mass_residual == 0.0 {
            0.0
        } else {
            v0.norm_squared()
        };
        let kinetic_change = v0.dot(retained_increment)
            + kinetic_relative
            + (0.5 * actual_mass) * initial_speed_squared;
        let active = incoming.mass_kg > 0.0 || outgoing.mass_kg > 0.0;
        let mean_velocity = if active { v0.add(mean) } else { Vec3::ZERO };
        let result = ParticleExchangeStep {
            state_after: next,
            duration_s,
            incoming,
            outgoing,
            mean_velocity_world_m_s: mean_velocity,
            mean_speed_squared_m2_s2: if active {
                mean_velocity.norm_squared() + variance
            } else {
                0.0
            },
            kinetic_energy_change_j: kinetic_change,
            retained_relative_energy_j: incoming.relative_energy_j - outgoing.relative_energy_j,
            external_impulse_n_s: external_impulse,
            external_work_j: v0.dot(external_impulse) + work_relative,
            mass_residual_kg: mass_residual,
            momentum_residual_kg_m_s: momentum_relative.add(v0.scale(mass_residual)),
            kinetic_energy_residual_j: v0.dot(momentum_relative)
                + kinetic_relative
                + transport_relative
                - work_relative
                + (0.5 * mass_residual) * initial_speed_squared,
        };
        if !result.state_after.position_m.is_finite()
            || !result.state_after.velocity_m_s.is_finite()
            || !result.mean_velocity_world_m_s.is_finite()
            || !result.external_impulse_n_s.is_finite()
            || !result.momentum_residual_kg_m_s.is_finite()
            || ![
                result.mean_speed_squared_m2_s2,
                result.kinetic_energy_change_j,
                result.retained_relative_energy_j,
                result.external_work_j,
                result.mass_residual_kg,
                result.kinetic_energy_residual_j,
            ]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(RocketError::Unrepresentable("particle exchange balance"));
        }
        poll(&mut cancelled)?;
        Ok(result)
    }
}

#[cfg(test)]
#[path = "exchange/tests.rs"]
mod tests;
