//! Passive recoil under the original local-speed-dependent Hertz traction.
//!
//! The original pressure integral is evaluated at every positive RK4 stage;
//! neither a centre-speed coefficient nor a homogeneous Coulomb ray replaces
//! it. With `q=[u_x,u_y,a*omega]`, `r=|W^-1/2 q|` and `k=3*pi/16`, a finite
//! local-speed interval with `0<mu_min<=mu(v)<=mu_max` gives
//! `N*mu_min*k/sqrt(1+k*k)*|q| <= H(q)` and
//! `|[F_x,F_y,M/a]| <= N*mu_max*sqrt(1+k*k)`.
//!
//! The initial energy, admitted numerical work defect, and original SPD
//! mobility bound the endpoint energy ball. Every endpoint/tail is checked
//! against that ball. Positive RK stages lie within 1.25 times its mass radius;
//! their local speeds and computed wrench norms are checked too. These bounds
//! support the shared finite-extinction tail without assuming degree-zero
//! traction. A positive finite-domain minimum suffices even when a weakening
//! law tends to zero at infinite speed. No coefficient floor is introduced.
//!
//! Without a global positive lower bound, the same original-law integrator
//! advances the finite horizon without assuming a stopping event. A conditional
//! extinction bound activates only inside a proved all-local-speed domain.
//! For clipped affine friction, speeds below half the cutoff give at least
//! half the original zero-speed coefficient. Pure viscosity instead uses its
//! exact Hertz second moment to bound the local linear decay rate: its safe
//! timestep does not collapse with speed, and no finite arrest is reported.
//! Proven zero traction is exact coasting, never arrest.
//! Spatial quadrature remains caller-owned: temporal error estimates do not
//! include its error. Refine pressure panels separately for nonlinear laws.
//! Load, radius, mobility, source curve and environment are held fixed. There
//! is no external drive, normal evolution, static branch selection or heating
//! feedback in this owner.

use core::f64::consts::PI;

use super::{HertzFrictionPatch, HertzFrictionQuadrature};
use crate::adhesive::spinning::transient::{
    CircularTransientExtinction, CircularTransientLaw, CircularTransientQuadrature,
    CircularTransientTraction, TransientError, advance_circular_transient,
    circular_transient_rate_domain,
};
use crate::sliding_patch::contact::relaxation::{CircularRelaxationInput, HertzRelaxationInput};
use crate::sliding_patch::contact::transient::{
    HertzTransientControls, HertzTransientError, HertzTransientReport,
};
use crate::sliding_patch::contact::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchError, HertzWrenchParameters,
};
use crate::{FrictionLaw, TriboError};
use fs_math::det;

const SPIN_RADIUS: f64 = 3.0 * PI / 16.0;
const BOUND_ROUNDOFF: f64 = 128.0 * f64::EPSILON;

/// Integrated original local-speed traction and its admitted whole-path bounds.
#[derive(Clone, Debug, PartialEq)]
pub struct HertzRateTransientReport {
    /// Original force/moment impulse, physical heat, travel and finite tail.
    pub trajectory: HertzTransientReport,
    /// Finite local-speed domain used for the nonlinear curve [m/s]. Genuine
    /// Coulomb uses infinity because its coefficient bounds hold globally.
    pub local_speed_bound_m_s: f64,
    /// Conservative original coefficient bounds on `[0,local_speed_bound]`.
    /// These are bounds, not an effective coefficient used by the trajectory.
    pub kinetic_coefficient_bounds: [f64; 2],
    /// Whole-path bound on `|[F_x,F_y,M/a]|` [N], including any stopping tail.
    /// Exact zero-traction coasting/rest returns zero.
    pub maximum_scaled_wrench_norm_n: f64,
}

/// Refusal leaves the original curve, rates, geometry and mechanics unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum HertzRateTransientError {
    /// Original shared friction curve or spatial pressure quadrature failed.
    Law(TriboError),
    /// Original circular mechanics, controls, work, effort or cancellation gate.
    Transient(HertzTransientError),
    /// The admitted speed domain has no representable positive traction bound
    /// and is not a proved exact zero-traction branch.
    NonCoercive {
        local_speed_bound_m_s: f64,
        kinetic_coefficient_bounds: [f64; 2],
    },
}
impl From<TriboError> for HertzRateTransientError {
    fn from(error: TriboError) -> Self {
        Self::Law(error)
    }
}
impl From<TransientError> for HertzRateTransientError {
    fn from(error: TransientError) -> Self {
        Self::Transient(error.with_exact())
    }
}
impl From<HertzWrenchError> for HertzRateTransientError {
    fn from(error: HertzWrenchError) -> Self {
        Self::Transient(HertzTransientError::Contact(error))
    }
}
impl From<HertzTransientError> for HertzRateTransientError {
    fn from(error: HertzTransientError) -> Self {
        Self::Transient(error)
    }
}
impl core::fmt::Display for HertzRateTransientError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rate-dependent Hertz trajectory: {self:?}")
    }
}
impl std::error::Error for HertzRateTransientError {}

impl HertzFrictionPatch {
    /// Integrate the original local-speed law under the complete body mobility.
    /// Positive pressure quadrature and actual physical power remain separate
    /// from temporal error and numerical work residuals. Each evaluation is a
    /// bounded original patch evaluation at the supplied spatial resolution.
    #[allow(clippy::too_many_lines)]
    pub fn advance_transient(
        &self,
        input: &HertzRelaxationInput,
        config: HertzWrenchConfig,
        quadrature: HertzFrictionQuadrature,
        controls: HertzTransientControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HertzRateTransientReport, HertzRateTransientError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        quadrature.validate()?;
        config.validate()?;
        let radius = self.contact_radius_m();
        let normal = self.normal_force_n();
        let upper_norm = det::sqrt(1.0 + SPIN_RADIUS * SPIN_RADIUS);
        if let FrictionLaw::Coulomb {
            static_mu,
            kinetic_mu,
        } = self.law
        {
            // Only this genuinely homogeneous original law may use the
            // existing invariant-ray admission and exact finite-stop shortcut.
            let law = HertzWrenchContact::new(HertzWrenchParameters {
                normal_force_n: normal,
                contact_radius_m: radius,
                static_mu,
                kinetic_mu,
            })?;
            let trajectory = law.advance_transient(input, config, controls, &mut cancelled)?;
            return Ok(HertzRateTransientReport {
                trajectory,
                local_speed_bound_m_s: f64::INFINITY,
                kinetic_coefficient_bounds: [kinetic_mu; 2],
                maximum_scaled_wrench_norm_n: checked(normal * kinetic_mu * upper_norm)?,
            });
        }
        if !input.dt_s.is_finite()
            || input.dt_s <= 0.0
            || input.velocity_m_s.into_iter().any(|v| !v.is_finite())
            || !input.spin_rate_rad_s.is_finite()
        {
            return Err(invalid("finite rates and positive duration required"));
        }
        let initial = [
            input.velocity_m_s[0],
            input.velocity_m_s[1],
            checked(radius * input.spin_rate_rad_s)?,
        ];
        if input.spin_rate_rad_s != 0.0 && initial[2] == 0.0 {
            return Err(invalid("unrepresentable radius-scaled spin rate"));
        }
        let (maximum_mass_radius, local_speed_bound_m_s) =
            circular_transient_rate_domain(input.scaled_mobility_per_kg, initial, config)?;
        // One actual admission evaluation, counted by the common integrator.
        let native = self.evaluate(input.velocity_m_s, input.spin_rate_rad_s, quadrature)?;
        poll(&mut cancelled)?;
        if det::hypot(native.force_n[0], native.force_n[1]) > config.maximum_force_n
            || native.spin_moment_nm.abs() > config.maximum_spin_moment_nm
        {
            return Err(HertzTransientError::EffortBudget.into());
        }
        let bounds = coefficient_bounds(&self.law, local_speed_bound_m_s)?;
        if self.proven_zero_traction(initial)? {
            if native.force_n != [0.0; 2]
                || native.spin_moment_nm != 0.0
                || native.dissipated_power_w != 0.0
            {
                return Err(invalid("unresolved exact zero-traction response"));
            }
            let trajectory = coast(input, initial == [0.0; 3])?;
            poll(&mut cancelled)?;
            return Ok(HertzRateTransientReport {
                trajectory,
                local_speed_bound_m_s,
                kinetic_coefficient_bounds: bounds,
                maximum_scaled_wrench_norm_n: 0.0,
            });
        }
        if bounds[1] <= 0.0 {
            return Err(HertzRateTransientError::NonCoercive {
                local_speed_bound_m_s,
                kinetic_coefficient_bounds: bounds,
            });
        }
        if native.force_n == [0.0; 2]
            && native.spin_moment_nm == 0.0
            && native.dissipated_power_w == 0.0
        {
            return Err(invalid("unresolved positive rate-dependent traction"));
        }
        let force_bound = checked(normal * bounds[1])?;
        let moment_bound = checked(force_bound * (SPIN_RADIUS * radius))?;
        let maximum_scaled_wrench_norm_n = checked(force_bound * upper_norm)?;
        let conditional_extinction = self.conditional_extinction(
            bounds,
            maximum_mass_radius,
            local_speed_bound_m_s,
            upper_norm,
        )?;
        let linear_wrench_norm_bound = if self.law.zero_speed_kinetic_coefficient()? == 0.0 {
            self.affine_slope()
                .filter(|slope| *slope > 0.0)
                .map(|slope| checked(normal * slope * (1.0 + BOUND_ROUNDOFF)))
                .transpose()?
        } else {
            None
        };
        let trajectory = advance_circular_transient(
            &CircularRelaxationInput {
                velocity_m_s: input.velocity_m_s,
                spin_rate_rad_s: input.spin_rate_rad_s,
                scaled_mobility_per_kg: input.scaled_mobility_per_kg,
                dt_s: input.dt_s,
            },
            config,
            controls,
            CircularTransientLaw {
                radius,
                sliding_capacity_n: force_bound,
                spinning_capacity_nm: moment_bound,
                power_norm_bounds: [
                    (bounds[0] / bounds[1]) * (SPIN_RADIUS / upper_norm) * (1.0 - BOUND_ROUNDOFF),
                    upper_norm,
                ],
                maximum_mass_radius,
                conditional_extinction,
                linear_wrench_norm_bound,
                quadrature: CircularTransientQuadrature::Rk4Doubling,
                evaluate: |velocity: [f64; 2], spin: f64| {
                    let speed = det::hypot(velocity[0], velocity[1]) + radius * spin.abs();
                    if !speed.is_finite() || speed > local_speed_bound_m_s {
                        return Err(HertzWrenchError::Invalid(
                            "rate-dependent stage left its admitted local-speed domain",
                        ));
                    }
                    let response = self
                        .evaluate(velocity, spin, quadrature)
                        .map_err(HertzWrenchError::Law)?;
                    let wrench_norm = det::hypot(
                        det::hypot(response.force_n[0], response.force_n[1]),
                        response.spin_moment_nm / radius,
                    );
                    if !wrench_norm.is_finite() || wrench_norm > maximum_scaled_wrench_norm_n {
                        return Err(HertzWrenchError::Invalid(
                            "spatial Hertz response exceeded its admitted wrench bound",
                        ));
                    }
                    if let Some(bound) = linear_wrench_norm_bound {
                        let rate = det::hypot(det::hypot(velocity[0], velocity[1]), radius * spin);
                        if wrench_norm > bound * rate {
                            return Err(HertzWrenchError::Invalid(
                                "viscous Hertz response exceeded its original linear bound",
                            ));
                        }
                    }
                    Ok(CircularTransientTraction {
                        force_n: response.force_n,
                        spin_moment_nm: response.spin_moment_nm,
                        dissipated_power_w: response.dissipated_power_w,
                        at_rest: response.at_rest,
                    })
                },
            },
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(HertzRateTransientReport {
            trajectory,
            local_speed_bound_m_s,
            kinetic_coefficient_bounds: bounds,
            maximum_scaled_wrench_norm_n,
        })
    }

    fn conditional_extinction(
        &self,
        bounds: [f64; 2],
        mass_radius: f64,
        speed_bound: f64,
        upper_norm: f64,
    ) -> Result<Option<CircularTransientExtinction>, HertzRateTransientError> {
        if bounds[0] > 0.0 {
            return Ok(None);
        }
        // These are domains for evaluating the ORIGINAL lower bound, not
        // replacement coefficients. Energy keeps the subsequent physical
        // tail inside the admitted ball; numerical moving stages still use
        // the original instantaneous clipped/nonlinear pressure integral.
        let speed = match self.law {
            FrictionLaw::VelocityDependent {
                mu_zero,
                slope_per_speed,
                ..
            } if mu_zero > 0.0 && slope_per_speed < 0.0 => 0.5 * (mu_zero / -slope_per_speed),
            FrictionLaw::Stribeck {
                characteristic_speed,
                ..
            }
            | FrictionLaw::Hyperbolic {
                characteristic_speed,
                ..
            } => characteristic_speed,
            FrictionLaw::ShiftedPower {
                characteristic_speed,
                exponent,
                ..
            } => characteristic_speed / exponent.max(1.0),
            _ => return Ok(None),
        };
        if !speed.is_finite() || speed <= 0.0 {
            return Ok(None);
        }
        let lower = coefficient_bounds(&self.law, speed)?[0];
        if lower <= 0.0 {
            return Ok(None);
        }
        let maximum_mass_radius =
            checked((speed / speed_bound) * mass_radius * (1.0 - 4.0 * BOUND_ROUNDOFF))?;
        let minimum_power_norm =
            checked((lower / bounds[1]) * (SPIN_RADIUS / upper_norm) * (1.0 - BOUND_ROUNDOFF))?;
        if maximum_mass_radius <= 0.0 || minimum_power_norm <= 0.0 {
            return Ok(None);
        }
        Ok(Some(CircularTransientExtinction {
            maximum_mass_radius,
            minimum_power_norm,
        }))
    }

    // Only exact autonomous zero-traction branches coast. A small quadrature
    // result or an underflowed positive response is not a proof of zero drag.
    fn proven_zero_traction(&self, q: [f64; 3]) -> Result<bool, TriboError> {
        if self.normal_force_n() == 0.0 || q == [0.0; 3] {
            return Ok(true);
        }
        if self.affine_slope() == Some(0.0) && self.law.zero_speed_kinetic_coefficient()? == 0.0 {
            return Ok(true);
        }
        if let FrictionLaw::VelocityDependent {
            slope_per_speed, ..
        } = self.law
        {
            if slope_per_speed < 0.0 {
                let minimum_speed = (det::hypot(q[0], q[1]) - q[2].abs()).max(0.0);
                return Ok(self.law.kinetic_mu(minimum_speed)? == 0.0);
            }
        }
        Ok(false)
    }
}

// Monotonic pieces use the shared original coefficient evaluation. For
// Stribeck, remove only the nonnegative viscous term to obtain a lower bound;
// its decreasing dry part is evaluated by the same original curve owner.
fn coefficient_bounds(law: &FrictionLaw, speed: f64) -> Result<[f64; 2], HertzRateTransientError> {
    let zero = law.zero_speed_kinetic_coefficient()?;
    let (lower, upper) = match *law {
        FrictionLaw::Stribeck {
            static_mu,
            kinetic_mu,
            characteristic_speed,
            viscous_per_speed,
        } => (
            FrictionLaw::Stribeck {
                static_mu,
                kinetic_mu,
                characteristic_speed,
                viscous_per_speed: 0.0,
            }
            .kinetic_mu(speed)?,
            checked(static_mu + checked(viscous_per_speed * speed)?)?,
        ),
        FrictionLaw::VelocityDependent {
            slope_per_speed, ..
        } if slope_per_speed > 0.0 => (zero, law.kinetic_mu(speed)?),
        _ => (law.kinetic_mu(speed)?, zero),
    };
    Ok([
        checked(lower * (1.0 - BOUND_ROUNDOFF))?,
        checked(upper * (1.0 + BOUND_ROUNDOFF))?,
    ])
}

fn coast(
    input: &HertzRelaxationInput,
    at_rest: bool,
) -> Result<HertzTransientReport, HertzRateTransientError> {
    Ok(HertzTransientReport {
        velocity_after_m_s: input.velocity_m_s,
        spin_rate_after_rad_s: input.spin_rate_rad_s,
        force_impulse_ns: [0.0; 2],
        spin_angular_impulse_nms: 0.0,
        slip_increment_m: [
            checked(input.dt_s * input.velocity_m_s[0])?,
            checked(input.dt_s * input.velocity_m_s[1])?,
        ],
        spin_increment_rad: checked(input.dt_s * input.spin_rate_rad_s)?,
        dissipated_energy_j: 0.0,
        prefix_dissipated_energy_j: 0.0,
        stopping_tail_heat_j: 0.0,
        contact_work_j: 0.0,
        work_residual_j: 0.0,
        prefix_work_residual_j: 0.0,
        estimated_scaled_rate_error_m_s: 0.0,
        estimated_scaled_impulse_error_ns: 0.0,
        estimated_scaled_travel_error_m: 0.0,
        estimated_heat_error_j: 0.0,
        tail_slip_bound_m: [0.0; 2],
        tail_spin_bound_rad: 0.0,
        stopping_time_bracket_s: at_rest.then_some([0.0; 2]),
        integrated_until_s: if at_rest { 0.0 } else { input.dt_s },
        scaled_velocity_residual_m_s: [0.0; 3],
        exact_ray: true,
        accepted_steps: 0,
        rejected_steps: 0,
        evaluations: 1,
    })
}

fn invalid(reason: &'static str) -> HertzRateTransientError {
    HertzTransientError::Invalid(reason).into()
}
fn checked(value: f64) -> Result<f64, HertzRateTransientError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable rate-dependent trajectory"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), HertzRateTransientError> {
    if cancelled() {
        Err(HertzTransientError::Cancelled.into())
    } else {
        Ok(())
    }
}
