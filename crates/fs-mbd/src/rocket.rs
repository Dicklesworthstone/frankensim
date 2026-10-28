//! Variable-mass point-rocket translation with an exact burnout event.
//!
//! During a powered interval, `m(t)=m0-mdot*t` and
//! `dv/dt = mdot*ve/m(t) + F/m(t) + g`. The effective-exhaust vector `ve`
//! points along thrust and already includes any nozzle pressure-thrust term.
//! Its integral is the Tsiolkovsky logarithm; the position integral is evaluated
//! separately, not by an average-mass or constant-acceleration approximation.
//! After propellant exhaustion, thrust vanishes and the vehicle coasts under
//! the same explicit gravity and external force at its dry mass.
//!
//! These inputs are piecewise constant over the requested step. A caller must
//! refine changing nozzle states, gravity, drag, guidance or environment. This
//! module supplies no atmosphere, propellant chemistry, attitude, staging,
//! rotating-frame force or empirical vehicle parameters.

use crate::Vec3;
use fs_math::det;

/// Inertial-frame point-vehicle state in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RocketState {
    /// Position in the declared inertial frame [m].
    pub position_m: Vec3,
    /// Velocity in that same frame [m/s].
    pub velocity_m_s: Vec3,
    /// Strictly positive vehicle mass that this engine may not consume [kg].
    pub dry_mass_kg: f64,
    /// Nonnegative remaining propellant inventory [kg].
    pub propellant_mass_kg: f64,
}

/// Explicit forcing, held fixed over one step including a possible coast.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RocketForcing {
    /// Nonnegative propellant consumption rate while fuel remains [kg/s].
    pub mass_flow_kg_s: f64,
    /// Thrust-direction vector with magnitude `total_nozzle_thrust/mdot` [m/s].
    /// It is the negative of the rearward exhaust-relative-velocity direction.
    /// Includes pressure thrust exactly once; ignored when the engine is off.
    pub effective_exhaust_velocity_m_s: Vec3,
    /// Explicit inertial-frame gravitational acceleration [m/s²].
    pub gravity_m_s2: Vec3,
    /// Other inertial-frame force [N], excluding propulsion and gravity.
    /// It remains active after burnout; its acceleration changes with mass.
    pub external_force_n: Vec3,
}

/// One complete accepted step; a refusal never returns a partially moved body.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RocketStep {
    /// State after both the powered interval and any remaining coast.
    pub state_after: RocketState,
    /// Propellant consumed during the powered interval [kg].
    pub burned_propellant_mass_kg: f64,
    /// Duration for which the specified mass flow operated [s].
    pub powered_duration_s: f64,
    /// Rest of the requested duration, with engine thrust exactly zero [s].
    pub coast_duration_s: f64,
    /// Engine thrust integrated over the powered duration [N s].
    /// This is not the vehicle momentum change: exhaust also carries momentum.
    pub thrust_impulse_n_s: Vec3,
}

/// Input, arithmetic, or cancellation refusal with no state mutation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RocketError {
    /// A named scalar or vector contains NaN or infinity.
    NonFiniteInput(&'static str),
    /// A named mass, rate or duration violates its physical domain.
    InvalidInput(&'static str),
    /// A required positive quantity or updated state is not representable.
    Unrepresentable(&'static str),
    /// Caller cancelled before the complete result was published.
    Cancelled,
}

impl core::fmt::Display for RocketError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "variable-mass rocket translation: {self:?}")
    }
}
impl std::error::Error for RocketError {}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RocketError> {
    if cancelled() {
        Err(RocketError::Cancelled)
    } else {
        Ok(())
    }
}

// Return L=ln(m0/m1), A=L/x, C=[x-(1-x)L]/x², x=burned/m0.
// Analytically cancelling the leading terms retains position and velocity
// even when m0 and m1 round to the same f64. C=1/2+x/6+x²/12+... .
fn mass_integrals(x: f64, initial_mass: f64, final_mass: f64) -> (f64, f64, f64) {
    if x < 0.125 {
        let (mut a, mut c) = (0.0_f64, 0.0_f64);
        for n in (0..=16).rev() {
            let k = f64::from(n + 1);
            a = a.mul_add(x, 1.0 / k);
            c = c.mul_add(x, 1.0 / (k * (k + 1.0)));
        }
        (x * a, a, c)
    } else {
        let l = if x < 0.5 {
            let remaining = 1.0 - x;
            // Stable -log1p(-x), including rounding in the subtraction.
            -det::ln(remaining) + (x + (remaining - 1.0)) / remaining
        } else {
            // The ratio can overflow or x can round to one while the positive
            // dry mass remains meaningful. Logarithms of each mass do not.
            det::ln(initial_mass) - det::ln(final_mass)
        };
        let a = l / x;
        let c = (-(final_mass / initial_mass)).mul_add(l, x) / (x * x);
        (l, a, c)
    }
}

impl RocketState {
    /// Advance exactly for piecewise-constant rate, effective exhaust vector,
    /// gravity and external force; split burnout inside the step and coast.
    ///
    /// Gravity is an acceleration, whereas the external input is a force and
    /// is integrated against the changing mass. At zero mass flow or empty
    /// propellant, the engine contributes no force regardless of its declared
    /// exhaust vector. A zero duration returns the unchanged admitted state.
    ///
    /// With burned fraction `x`, powered time `b`, and coast time `c`, the
    /// engine's displacement multiplier is `b*x*C + c*L` and the external
    /// force's is `b²*C/m0 + b*c*A/m0 + c²/(2*m1)`, where `A=L/x` and
    /// `C=[x-(1-x)*L]/x²`. These are exact integrals within floating roundoff.
    /// No mass is clipped to hide an overlong step: exhaustion is a physical
    /// event that sets remaining propellant to zero and ends propulsion.
    pub fn advance(
        &self,
        forcing: RocketForcing,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RocketStep, RocketError> {
        poll(&mut cancelled)?;
        for (name, value) in [
            ("position_m", self.position_m),
            ("velocity_m_s", self.velocity_m_s),
            (
                "effective_exhaust_velocity_m_s",
                forcing.effective_exhaust_velocity_m_s,
            ),
            ("gravity_m_s2", forcing.gravity_m_s2),
            ("external_force_n", forcing.external_force_n),
        ] {
            if !value.is_finite() {
                return Err(RocketError::NonFiniteInput(name));
            }
        }
        for (name, value) in [
            ("dry_mass_kg", self.dry_mass_kg),
            ("propellant_mass_kg", self.propellant_mass_kg),
            ("mass_flow_kg_s", forcing.mass_flow_kg_s),
            ("duration_s", duration_s),
        ] {
            if !value.is_finite() {
                return Err(RocketError::NonFiniteInput(name));
            }
            if value < 0.0 || (name == "dry_mass_kg" && value == 0.0) {
                return Err(RocketError::InvalidInput(name));
            }
        }
        let initial_mass = self.dry_mass_kg + self.propellant_mass_kg;
        if !initial_mass.is_finite() {
            return Err(RocketError::Unrepresentable("wet mass"));
        }
        let mut next = *self;
        let (mut burned, mut powered, mut coast) = (0.0, 0.0, duration_s);
        let (mut log_ratio, mut a, mut c) = (0.0, 1.0, 0.5);
        let mut fraction = 0.0;
        if duration_s > 0.0 && forcing.mass_flow_kg_s > 0.0 && self.propellant_mass_kg > 0.0 {
            let burnout = self.propellant_mass_kg / forcing.mass_flow_kg_s;
            if burnout == 0.0 {
                return Err(RocketError::Unrepresentable("burnout duration"));
            }
            if burnout <= duration_s {
                powered = burnout;
                coast = duration_s - powered;
                burned = self.propellant_mass_kg;
                next.propellant_mass_kg = 0.0;
            } else {
                powered = duration_s;
                coast = 0.0;
                burned = forcing.mass_flow_kg_s * powered;
                next.propellant_mass_kg -= burned;
                if next.propellant_mass_kg <= 0.0
                    || next.propellant_mass_kg == self.propellant_mass_kg
                {
                    return Err(RocketError::Unrepresentable("remaining propellant"));
                }
            }
            fraction = burned / initial_mass;
            if !fraction.is_finite() || fraction <= 0.0 || fraction > 1.0 {
                return Err(RocketError::Unrepresentable("burned mass fraction"));
            }
            (log_ratio, a, c) = mass_integrals(
                fraction,
                initial_mass,
                self.dry_mass_kg + next.propellant_mass_kg,
            );
            if !log_ratio.is_finite()
                || log_ratio <= 0.0
                || !a.is_finite()
                || a <= 0.0
                || !c.is_finite()
                || c <= 0.0
            {
                return Err(RocketError::Unrepresentable("variable-mass integrals"));
            }
        }
        poll(&mut cancelled)?;
        let final_mass = self.dry_mass_kg + next.propellant_mass_kg;
        let engine_position = powered * fraction * c + coast * log_ratio;
        let force_velocity = (powered / initial_mass) * a + coast / final_mass;
        let force_position = (powered / initial_mass) * (powered * c + coast * a)
            + (0.5 * coast / final_mass) * coast;
        next.velocity_m_s = self
            .velocity_m_s
            .add(forcing.effective_exhaust_velocity_m_s.scale(log_ratio))
            .add(forcing.external_force_n.scale(force_velocity))
            .add(forcing.gravity_m_s2.scale(duration_s));
        next.position_m = self
            .position_m
            .add(self.velocity_m_s.scale(duration_s))
            .add(
                forcing
                    .effective_exhaust_velocity_m_s
                    .scale(engine_position),
            )
            .add(forcing.external_force_n.scale(force_position))
            .add(
                forcing
                    .gravity_m_s2
                    .scale(0.5 * duration_s)
                    .scale(duration_s),
            );
        let impulse = forcing.effective_exhaust_velocity_m_s.scale(burned);
        if !next.position_m.is_finite() || !next.velocity_m_s.is_finite() || !impulse.is_finite() {
            return Err(RocketError::Unrepresentable(
                "rocket state or thrust impulse",
            ));
        }
        poll(&mut cancelled)?;
        Ok(RocketStep {
            state_after: next,
            burned_propellant_mass_kg: burned,
            powered_duration_s: powered,
            coast_duration_s: coast,
            thrust_impulse_n_s: impulse,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn state() -> RocketState {
        RocketState {
            position_m: Vec3::ZERO,
            velocity_m_s: Vec3::ZERO,
            dry_mass_kg: 100.0,
            propellant_mass_kg: 100.0,
        }
    }
    fn engine() -> RocketForcing {
        RocketForcing {
            mass_flow_kg_s: 10.0,
            effective_exhaust_velocity_m_s: Vec3::new(1000.0, 0.0, 0.0),
            gravity_m_s2: Vec3::ZERO,
            external_force_n: Vec3::ZERO,
        }
    }
    fn near(a: f64, b: f64) {
        assert!((a - b).abs() <= 2e-12 * b.abs().max(1e-300), "{a} != {b}");
    }
    fn vector_near(a: Vec3, b: Vec3) {
        for (u, v) in [(a.x, b.x), (a.y, b.y), (a.z, b.z)] {
            assert!((u - v).abs() <= 2e-11 * v.abs().max(1.0), "{a:?} != {b:?}");
        }
    }

    #[test]
    fn exact_vacuum_tsiolkovsky_velocity_and_position() {
        let step = state().advance(engine(), 5.0, || false).unwrap();
        let log = (4.0_f64 / 3.0).ln();
        near(step.state_after.velocity_m_s.x, 1000.0 * log);
        near(step.state_after.position_m.x, 1000.0 * (5.0 - 15.0 * log));
        assert_eq!(step.state_after.propellant_mass_kg, 50.0);
        assert_eq!(step.state_after.dry_mass_kg, 100.0);
        assert_eq!(step.thrust_impulse_n_s.x, 50000.0);
        assert_eq!(step.powered_duration_s, 5.0);
        assert_eq!(step.coast_duration_s, 0.0);
        assert_eq!(state().advance(engine(), 5.0, || false).unwrap(), step);
    }

    #[test]
    fn burnout_inside_step_coasts_without_consuming_dry_mass_or_retaining_thrust() {
        let s = RocketState {
            propellant_mass_kg: 20.0,
            ..state()
        };
        let step = s.advance(engine(), 4.0, || false).unwrap();
        let log = 1.2_f64.ln();
        assert_eq!(step.powered_duration_s, 2.0);
        assert_eq!(step.coast_duration_s, 2.0);
        assert_eq!(step.burned_propellant_mass_kg, 20.0);
        assert_eq!(step.state_after.propellant_mass_kg, 0.0);
        near(step.state_after.velocity_m_s.x, 1000.0 * log);
        near(step.state_after.position_m.x, 1000.0 * (2.0 - 8.0 * log));
        let empty = step.state_after.advance(engine(), 3.0, || false).unwrap();
        assert_eq!(empty.thrust_impulse_n_s, Vec3::ZERO);
        assert_eq!(
            empty.state_after.velocity_m_s,
            step.state_after.velocity_m_s
        );
        near(
            empty.state_after.position_m.x,
            step.state_after.position_m.x + 3.0 * step.state_after.velocity_m_s.x,
        );
    }

    #[test]
    fn zero_flow_coasts_under_explicit_gravity_and_nonpropulsive_force() {
        let f = RocketForcing {
            mass_flow_kg_s: 0.0,
            gravity_m_s2: Vec3::new(0.0, -9.81, 0.0),
            external_force_n: Vec3::new(200.0, 0.0, 0.0),
            ..engine()
        };
        let step = state().advance(f, 3.0, || false).unwrap();
        vector_near(step.state_after.velocity_m_s, Vec3::new(3.0, -29.43, 0.0));
        vector_near(step.state_after.position_m, Vec3::new(4.5, -44.145, 0.0));
        assert_eq!(
            step.state_after.propellant_mass_kg,
            state().propellant_mass_kg
        );
        assert_eq!(step.thrust_impulse_n_s, Vec3::ZERO);
        assert_eq!(
            state().advance(f, 0.0, || false).unwrap().state_after,
            state()
        );
    }

    #[test]
    fn external_force_uses_the_changing_mass_in_both_integrals() {
        let f = RocketForcing {
            effective_exhaust_velocity_m_s: Vec3::ZERO,
            external_force_n: Vec3::new(500.0, 250.0, -100.0),
            ..engine()
        };
        let result = state().advance(f, 5.0, || false).unwrap();
        let log = (4.0_f64 / 3.0).ln();
        vector_near(
            result.state_after.velocity_m_s,
            f.external_force_n.scale(log / 10.0),
        );
        vector_near(
            result.state_after.position_m,
            f.external_force_n.scale((5.0 - 15.0 * log) / 10.0),
        );
        assert_eq!(result.thrust_impulse_n_s, Vec3::ZERO);
    }

    #[test]
    fn split_steps_match_one_step_even_across_burnout_with_external_force() {
        let s = RocketState {
            propellant_mass_kg: 25.0,
            position_m: Vec3::new(3.0, -4.0, 5.0),
            velocity_m_s: Vec3::new(10.0, 20.0, -30.0),
            ..state()
        };
        let f = RocketForcing {
            effective_exhaust_velocity_m_s: Vec3::new(800.0, 200.0, -100.0),
            external_force_n: Vec3::new(30.0, -40.0, 50.0),
            gravity_m_s2: Vec3::new(0.0, -9.81, 0.0),
            ..engine()
        };
        let whole = s.advance(f, 6.0, || false).unwrap();
        for first_time in [0.01, 1.0, 2.5, 4.0] {
            let first = s.advance(f, first_time, || false).unwrap();
            let second = first
                .state_after
                .advance(f, 6.0 - first_time, || false)
                .unwrap();
            vector_near(whole.state_after.position_m, second.state_after.position_m);
            vector_near(
                whole.state_after.velocity_m_s,
                second.state_after.velocity_m_s,
            );
            assert_eq!(second.state_after.propellant_mass_kg, 0.0);
            assert_eq!(second.state_after.dry_mass_kg, s.dry_mass_kg);
        }
    }

    #[test]
    fn tiny_burn_keeps_velocity_and_position_when_total_mass_rounds_unchanged() {
        let s = RocketState {
            dry_mass_kg: 1e16,
            propellant_mass_kg: 1.0,
            ..state()
        };
        let f = RocketForcing {
            mass_flow_kg_s: 1e-4,
            ..engine()
        };
        let step = s.advance(f, 0.01, || false).unwrap();
        assert_eq!(
            s.dry_mass_kg + s.propellant_mass_kg,
            step.state_after.dry_mass_kg + step.state_after.propellant_mass_kg
        );
        near(step.state_after.velocity_m_s.x, 1e-19);
        near(step.state_after.position_m.x, 5e-22);
        assert!(step.state_after.propellant_mass_kg < s.propellant_mass_kg);
    }

    #[test]
    fn translation_rotation_and_galilean_changes_preserve_the_solution() {
        fn rotate(v: Vec3) -> Vec3 {
            Vec3::new(v.z, v.x, v.y)
        }
        let s = state();
        let f = RocketForcing {
            effective_exhaust_velocity_m_s: Vec3::new(800.0, 200.0, -100.0),
            external_force_n: Vec3::new(30.0, -40.0, 50.0),
            gravity_m_s2: Vec3::new(0.0, -9.81, 0.0),
            ..engine()
        };
        let result = s.advance(f, 4.0, || false).unwrap().state_after;
        let rotated = RocketState {
            position_m: rotate(s.position_m),
            velocity_m_s: rotate(s.velocity_m_s),
            ..s
        }
        .advance(
            RocketForcing {
                effective_exhaust_velocity_m_s: rotate(f.effective_exhaust_velocity_m_s),
                external_force_n: rotate(f.external_force_n),
                gravity_m_s2: rotate(f.gravity_m_s2),
                ..f
            },
            4.0,
            || false,
        )
        .unwrap()
        .state_after;
        assert_eq!(rotated.position_m, rotate(result.position_m));
        assert_eq!(rotated.velocity_m_s, rotate(result.velocity_m_s));
        let translation = Vec3::new(500.0, -20.0, 30.0);
        let boost = Vec3::new(7.0, 8.0, 9.0);
        let moved = RocketState {
            position_m: translation,
            velocity_m_s: boost,
            ..s
        }
        .advance(f, 4.0, || false)
        .unwrap()
        .state_after;
        vector_near(
            moved.position_m,
            result.position_m.add(translation).add(boost.scale(4.0)),
        );
        vector_near(moved.velocity_m_s, result.velocity_m_s.add(boost));
    }

    #[test]
    fn invalid_inputs_and_every_cancellation_boundary_leave_the_source_untouched() {
        let s = state();
        let f = engine();
        for bad in [
            RocketState {
                dry_mass_kg: 0.0,
                ..s
            },
            RocketState {
                propellant_mass_kg: -1.0,
                ..s
            },
            RocketState {
                position_m: Vec3::new(f64::NAN, 0.0, 0.0),
                ..s
            },
        ] {
            assert!(bad.advance(f, 1.0, || false).is_err());
        }
        for bad in [
            RocketForcing {
                mass_flow_kg_s: -1.0,
                ..f
            },
            RocketForcing {
                effective_exhaust_velocity_m_s: Vec3::new(f64::INFINITY, 0.0, 0.0),
                ..f
            },
        ] {
            assert!(s.advance(bad, 1.0, || false).is_err());
        }
        for dt in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(s.advance(f, dt, || false).is_err());
        }
        for stop in 1..=3 {
            let mut calls = 0;
            assert_eq!(
                s.advance(f, 20.0, || {
                    calls += 1;
                    calls == stop
                }),
                Err(RocketError::Cancelled)
            );
            assert_eq!(s, state());
        }
    }
}
