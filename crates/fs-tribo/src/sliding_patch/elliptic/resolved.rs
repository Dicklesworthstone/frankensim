//! Local kinetic ellipses on the original elliptical Hertz pressure measure.
//!
//! The resolver supplies the original law at each actual rigid slip vector.
//! Geometry, pressure weights, associated traction, torque and heat remain here.
//! Signed rays through zero slip remove the traction-direction singularity.
//! A sine substitution on each chord removes the Hertz rim square root. The
//! exterior formula uses inverse distance, retaining finite coordinates when
//! zero slip is far from the patch. This is bounded quadrature, not an error
//! certificate; refine both resolutions for sharp source variation or nearly
//! rank-one friction.

use core::f64::consts::PI;

use super::*;

/// Explicit spatial work budget. The native quadrature allocates no memory.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ResolvedHertzQuadrature {
    /// Four nodes per panel in each of at most two radial intervals.
    pub radial_panels: usize,
    /// Four nodes per panel in the angular ray coordinate.
    pub angular_panels: usize,
    /// Actual resolver calls, admitted before the first call; at most 65,536.
    pub maximum_evaluations: usize,
}

impl ResolvedHertzQuadrature {
    pub fn validate(self) -> Result<(), TriboError> {
        if !(1..=256).contains(&self.radial_panels)
            || !(1..=256).contains(&self.angular_panels)
            || !(1..=65_536).contains(&self.maximum_evaluations)
        {
            return Err(invalid("resolved_hertz_quadrature"));
        }
        Ok(())
    }
}

/// One pressure-weighted node; velocity and position share the supplied frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzLocalPoint {
    pub position_m: [f64; 2],
    pub slip_velocity_m_s: [f64; 2],
    pub normal_force_weight_n: f64,
}

/// Original force/torque/heat, plus the independently integrated pressure load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResolvedHertzWrench {
    pub wrench: HertzSlidingWrench,
    pub evaluations: usize,
    pub integrated_normal_force_n: f64,
    /// No renormalization conceals pressure quadrature error.
    pub normal_force_residual_n: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ResolvedHertzError<E> {
    Native(TriboError),
    Source(E),
    EvaluationBudget { required: usize, maximum: usize },
    Cancelled,
}

impl<E> From<TriboError> for ResolvedHertzError<E> {
    fn from(error: TriboError) -> Self {
        Self::Native(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for ResolvedHertzError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "resolved elliptical Hertz friction: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for ResolvedHertzError<E> {}

/// Immutable pressure geometry, with no spatially uniform coefficient imposed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResolvedEllipticHertzPatch {
    parameters: EllipticHertzPatchParameters,
}

impl ResolvedEllipticHertzPatch {
    pub fn new(parameters: EllipticHertzPatchParameters) -> Result<Self, TriboError> {
        Ok(Self {
            parameters: normalize_parameters(parameters)?,
        })
    }

    #[must_use]
    pub const fn parameters(&self) -> EllipticHertzPatchParameters {
        self.parameters
    }

    /// Resolve and integrate one original kinetic ellipse per pressure node.
    ///
    /// Returned laws already carry their local material direction in the current
    /// tangent frame. No shear history or physical time advances. Resolver side
    /// effects must remain private candidates until this whole call succeeds.
    /// Cancellation is checked before and after every resolver invocation.
    pub fn evaluate<E, C: FnMut() -> bool>(
        &self,
        slip: [f64; 2],
        spin: f64,
        controls: ResolvedHertzQuadrature,
        mut resolve: impl FnMut(HertzLocalPoint, &mut C) -> Result<EllipticCoulomb, E>,
        cancelled: &mut C,
    ) -> Result<ResolvedHertzWrench, ResolvedHertzError<E>> {
        controls.validate()?;
        for x in slip {
            finite(x, "resolved_hertz_slip")?;
        }
        finite(spin, "resolved_hertz_spin")?;
        poll_local(cancelled)?;
        let p = self.parameters;
        let mut q = [0.0; 2];
        if spin != 0.0 {
            q = to_principal(p.footprint_axis, turn(slip))?;
            for (i, x) in q.iter_mut().enumerate() {
                *x = finite(*x / p.semi_axes_m[i], "resolved_hertz_inverse_slip")?;
            }
        }
        let qnorm = finite(det::hypot(q[0], q[1]), "resolved_hertz_inverse_speed")?;
        let direction = if qnorm == 0.0 {
            [1.0, 0.0]
        } else {
            q.map(|x| spin.signum() * (x / qnorm))
        };
        let cross_direction = turn(direction);
        let inside = spin != 0.0 && qnorm < spin.abs();
        let ratio = if spin == 0.0 {
            0.0
        } else if inside {
            qnorm / spin.abs()
        } else {
            spin.abs() / qnorm
        };
        if spin != 0.0 && qnorm > 0.0 && ratio == 0.0 {
            return Err(invalid("resolved_hertz_unrepresentable_motion_ratio").into());
        }
        let intervals = if inside { 2 } else { 1 };
        let required = 16 * controls.radial_panels * controls.angular_panels * intervals;
        if required > controls.maximum_evaluations {
            return Err(ResolvedHertzError::EvaluationBudget {
                required,
                maximum: controls.maximum_evaluations,
            });
        }
        // Fixed normalized radial nodes are shared by every shifted chord.
        let mut radial = [(0.0, 0.0); 4 * 256];
        let half = 0.5 / controls.radial_panels as f64;
        for panel in 0..controls.radial_panels {
            poll_local(cancelled)?;
            for (i, (node, weight)) in GAUSS.into_iter().enumerate() {
                radial[4 * panel + i] = ((2.0 * panel as f64 + 1.0 + node) * half, half * weight);
            }
        }
        let mut sum = [Dd::ZERO; 4];
        let mut load = Dd::ZERO;
        let angular_range = if spin == 0.0 { 2.0 * PI } else { PI };
        let angular_start = if spin != 0.0 && !inside {
            -FRAC_PI_2
        } else {
            0.0
        };
        let angular_half = angular_range / (2.0 * controls.angular_panels as f64);
        for panel in 0..controls.angular_panels {
            for (node, angular_weight) in GAUSS {
                let theta = angular_start + (2.0 * panel as f64 + 1.0 + node) * angular_half;
                let (s, c) = det::sin_cos(theta);
                // Inside: z = d e1 + t n, t = R sin(psi) - d cos(theta).
                // The pressure measure is R^2 cos(psi)^2 |t| dpsi dtheta.
                // Split at t=0 so traction direction is constant on each ray.
                let radius = if inside {
                    det::sqrt((1.0 - ratio) * (1.0 + ratio) + (ratio * c) * (ratio * c))
                } else {
                    1.0
                };
                let centre = ratio * c;
                let split = if inside {
                    det::asin(centre / radius)
                } else {
                    FRAC_PI_2
                };
                let boundaries = if spin == 0.0 {
                    [0.0, FRAC_PI_2, FRAC_PI_2]
                } else {
                    [-FRAC_PI_2, split, FRAC_PI_2]
                };
                // Outside: sin(actual ray angle)=k sin(theta), k=1/d.
                // Never form d or subtract remote O(d) coordinates.
                let ray_cos = if !inside && spin != 0.0 {
                    det::sqrt((1.0 - ratio) * (1.0 + ratio) + (ratio * c) * (ratio * c))
                } else {
                    c
                };
                let n = [0, 1].map(|i| {
                    if inside {
                        c * direction[i] + s * cross_direction[i]
                    } else {
                        -ray_cos * direction[i] + ratio * s * cross_direction[i]
                    }
                });
                let jn = turn(from_principal(
                    p.footprint_axis,
                    [p.semi_axes_m[0] * n[0], p.semi_axes_m[1] * n[1]],
                )?);
                for interval in boundaries.windows(2) {
                    let width = interval[1] - interval[0];
                    if width == 0.0 {
                        continue;
                    }
                    for &(fraction, radial_weight) in &radial[..4 * controls.radial_panels] {
                        poll_local(cancelled)?;
                        let psi = interval[0] + width * fraction;
                        let (rs, rc) = det::sin_cos(psi);
                        let (z, velocity, measure) = if spin == 0.0 {
                            ([rs * c, rs * s], slip, rs * rc * rc)
                        } else if inside {
                            let t = radius * rs - centre;
                            (
                                [0, 1].map(|i| ratio * direction[i] + t * n[i]),
                                jn.map(|x| (spin * t) * x),
                                radius * radius * rc * rc * t.abs(),
                            )
                        } else {
                            let x = ratio * s * s - c * ray_cos * rs;
                            let y = s * (ray_cos + ratio * c * rs);
                            let rate = spin.signum() * (qnorm * ray_cos + spin.abs() * c * rs);
                            (
                                [0, 1].map(|i| x * direction[i] + y * cross_direction[i]),
                                jn.map(|x| rate * x),
                                c * c * c * rc * rc * (1.0 + ratio * c * rs / ray_cos),
                            )
                        };
                        let position = from_principal(
                            p.footprint_axis,
                            [p.semi_axes_m[0] * z[0], p.semi_axes_m[1] * z[1]],
                        )?;
                        for x in velocity {
                            finite(x, "resolved_hertz_local_slip")?;
                        }
                        let weight = (3.0 / (2.0 * PI))
                            * measure
                            * width
                            * radial_weight
                            * angular_half
                            * angular_weight;
                        nonnegative_finite(weight, "resolved_hertz_pressure_weight")?;
                        let point = HertzLocalPoint {
                            position_m: position,
                            slip_velocity_m_s: velocity,
                            normal_force_weight_n: finite(
                                p.normal_force_n * weight,
                                "resolved_hertz_pressure_weight",
                            )?,
                        };
                        let law = resolve(point, cancelled).map_err(ResolvedHertzError::Source)?;
                        poll_local(cancelled)?;
                        let texture = law.initial_state().principal_axis();
                        let local = to_principal(texture, velocity)?;
                        let capacities = law
                            .kinetic_coefficients()
                            .map(|mu| mu * point.normal_force_weight_n);
                        for cap in capacities {
                            finite(cap, "resolved_hertz_force_capacity")?;
                        }
                        let resistance = support(capacities, local)?;
                        let force = from_principal(texture, resistance.map(|x| -x))?;
                        let values = [
                            Dd::from_f64(force[0]),
                            Dd::from_f64(force[1]),
                            dot(turn(position), force),
                            dot(resistance, local),
                        ];
                        for (total, value) in sum.iter_mut().zip(values) {
                            *total = *total + value;
                        }
                        load = load + Dd::from_f64(point.normal_force_weight_n);
                    }
                }
            }
        }
        let integrated_normal_force_n = finite(load.to_f64(), "resolved_hertz_integrated_load")?;
        let wrench = finish(sum, slip, spin, slip == [0.0; 2] && spin == 0.0)?;
        poll_local(cancelled)?;
        Ok(ResolvedHertzWrench {
            wrench,
            evaluations: required,
            integrated_normal_force_n,
            normal_force_residual_n: finite(
                (load - Dd::from_f64(p.normal_force_n)).to_f64(),
                "resolved_hertz_load_residual",
            )?,
        })
    }
}

fn poll_local<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), ResolvedHertzError<E>> {
    if cancelled() {
        Err(ResolvedHertzError::Cancelled)
    } else {
        Ok(())
    }
}
