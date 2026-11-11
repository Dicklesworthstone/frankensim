//! Fully developed axial Couette--Poiseuille flow in a concentric annulus.
//!
//! Cylindrical momentum fixes `tau(r)=C/r-delta_p*r/(2L)`. ONE bounded
//! constitutive solve enforces both axial wall velocities using the original
//! thermal Carreau law. Logarithmic radius resolves the inner cylinder; gap
//! and normalized annular area avoid cancellation between squared radii.
//!
//! This is a fixed-geometry, no-slip, incompressible laminar reduction at one
//! uniform frozen temperature. The caller owns entrance/transition admission,
//! pressure limits and thermal evolution. Pure cylinder rotation is provided
//! separately by [`rotation`], and [`helical`] resolves their coupled invariant
//! shear. Nonlinear profiles cannot be superposed. Eccentricity, cavitation, slip and density transport are not
//! supplied. A zero inner radius is a pipe with an axis, not a material wall.

use super::pipe::{CarreauPipeControls, PreparedCarreau, RateState, integrate, product_ratio};
use super::slit::{Moments, integrate_moments};
use super::{MaterialError, ThermalCarreauYasuda, parameters, positive, state_error};
use fs_math::det;

/// Pure circular Couette flow on the same source and concentric geometry.
pub mod rotation;

/// Coupled pressure, axial wall motion and rotation using one local viscosity.
pub mod helical;

#[cfg(test)]
mod tests;

/// Explicit local work bounds and no-slip tolerance, not material uncertainty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauAnnularControls {
    pub local: CarreauPipeControls,
    /// Relative wall-speed and constitutive correction tolerance. Axial flow
    /// uses max(|Uouter-Uinner|, gap*max wall rate); pure rotation uses the
    /// nonzero angular-speed difference and the inner-wall shear rate.
    pub relative_wall_velocity_tolerance: f64,
    /// From one through 128 bounded stress corrections.
    pub maximum_wall_stress_iterations: u32,
}

impl Default for CarreauAnnularControls {
    fn default() -> Self {
        Self {
            local: CarreauPipeControls::default(),
            relative_wall_velocity_tolerance: 1e-9,
            maximum_wall_stress_iterations: 96,
        }
    }
}

impl CarreauAnnularControls {
    fn admit(self) -> Result<(), MaterialError> {
        self.local.admit()?;
        if !positive(self.relative_wall_velocity_tolerance)
            || self.relative_wall_velocity_tolerance >= 1.0
            || !(1..=128).contains(&self.maximum_wall_stress_iterations)
        {
            return Err(parameters("invalid Carreau annular wall-speed budgets"));
        }
        Ok(())
    }
}

/// Original source, fixed temperature, curved geometry and supported rate span.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminarCarreauAnnulus {
    prepared: PreparedCarreau,
    inner_radius_m: f64,
    outer_radius_m: f64,
    length_m: f64,
    gap_m: f64,
    mean_radius_m: f64,
    log_ratio: f64,
    log_inner: f64,
    langevin: f64,
    area_m2: f64,
}

/// Accepted pressure, wall motion and one actual curved velocity profile.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauAnnularFlow {
    pub pressure_drop_pa: f64,
    /// Inner/outer AXIAL wall speeds [m/s].
    pub requested_wall_velocities_m_s: [f64; 2],
    pub accepted_wall_velocities_m_s: [f64; 2],
    /// Accepted minus requested outer-minus-inner speed [m/s].
    pub wall_velocity_residual_m_s: f64,
    pub volumetric_flow_m3_s: f64,
    pub mean_wall_entrainment_m3_s: f64,
    /// Flow relative to the arithmetic mean wall speed. Curved Couette flow
    /// contributes here even at zero pressure; this is not just pressure flow.
    pub relative_flow_m3_s: f64,
    /// Fluid-domain outward-normal stress tractions, in inner/outer order.
    pub wall_tractions_pa: [f64; 2],
    /// Axial forces on the fluid [N], including each actual cylindrical area.
    pub wall_forces_n: [f64; 2],
    pub pressure_power_w: f64,
    /// Work at the REQUESTED wall speeds; individual input ports may absorb work.
    pub wall_power_w: f64,
    /// Positive volume integral of tau*du/dr, never a fitted loss coefficient.
    pub dissipated_power_w: f64,
    /// Requested pressure/wall work minus actual-profile heat, without the
    /// cancellation of equal opposite common-translation work.
    pub work_residual_w: f64,
    /// dQ/d(delta_p) at fixed wall speeds [m^3/(s Pa)], strictly positive.
    pub pressure_conductance_m3_s_pa: f64,
    /// dQ/d(Uouter-Uinner) at fixed mean speed and pressure [m²].
    /// Its negative is d(relative resisting force)/d(delta_p).
    pub flow_relative_wall_velocity_derivative_m2: f64,
    pub relative_wall_resistance_n_s_m: f64,
    pub maximum_shear_rate_per_s: f64,
    /// Largest observed positive-integral refinement change; not an error bound.
    pub quadrature_relative_change: f64,
    pub wall_stress_iterations: u32,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauAnnularProfile {
    annulus: LaminarCarreauAnnulus,
    controls: CarreauAnnularControls,
    flow: CarreauAnnularFlow,
    offset_pa: f64,
    pressure_stress_pa: f64,
    cuts: Cuts,
    rate_scale: f64,
    stress_scale: f64,
}

/// Prescribed signed net flow and axial wall motion resolved by one profile.
/// Requested flow is retained separately; physical work uses accepted flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauAnnularFlowSolution {
    profile: CarreauAnnularProfile,
    requested_volumetric_flow_m3_s: f64,
    flow_residual_m3_s: f64,
    flow_iterations: u32,
}

impl CarreauAnnularFlowSolution {
    pub const fn profile(&self) -> &CarreauAnnularProfile {
        &self.profile
    }
    pub const fn requested_volumetric_flow_m3_s(&self) -> f64 {
        self.requested_volumetric_flow_m3_s
    }
    /// Accepted minus requested flow after removing common-wall entrainment.
    /// Large common translation therefore does not mask a hydraulic residual.
    pub const fn flow_residual_m3_s(&self) -> f64 {
        self.flow_residual_m3_s
    }
    pub const fn flow_iterations(&self) -> u32 {
        self.flow_iterations
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauAnnularSample {
    pub axial_velocity_m_s: f64,
    /// tau_zr, positive for velocity increasing with radius.
    pub shear_stress_pa: f64,
    pub radial_velocity_gradient_per_s: f64,
    pub dynamic_viscosity_pa_s: f64,
    pub dissipation_w_m3: f64,
    pub quadrature_relative_change: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Cuts {
    points: [f64; 12],
    count: usize,
}

impl Cuts {
    fn push(&mut self, point: f64) -> Result<(), MaterialError> {
        if point > 0.0 && point < 1.0 {
            // One stress zero OR extremum, one area midpoint, and at most four
            // source-transition crossings fit in this fixed radial scratch.
            if self.count == self.points.len() {
                return Err(state_error("Carreau annular radial split budget exhausted"));
            }
            self.points[self.count] = point;
            self.count += 1;
        }
        Ok(())
    }
    fn sort(&mut self) {
        self.points[..self.count].sort_unstable_by(f64::total_cmp);
    }
}

#[derive(Clone, Copy)]
enum Integral {
    Velocity,
    Flow,
    Heat,
}

struct Evaluation {
    relative_speed: f64,
    derivative: f64,
    moments: Moments,
    tangent_scale: f64,
    rate_scale: f64,
    stress_scale: f64,
    cuts: Cuts,
    change: f64,
}

impl LaminarCarreauAnnulus {
    pub fn new(
        law: ThermalCarreauYasuda,
        temperature_k: f64,
        inner_radius_m: f64,
        outer_radius_m: f64,
        length_m: f64,
        maximum_laminar_shear_rate_per_s: f64,
    ) -> Result<Self, MaterialError> {
        if ![inner_radius_m, outer_radius_m, length_m]
            .into_iter()
            .all(positive)
            || outer_radius_m <= inner_radius_m
        {
            return Err(parameters(
                "Carreau annulus needs positive ordered radii and length",
            ));
        }
        let gap_m = outer_radius_m - inner_radius_m;
        let mean_radius_m = inner_radius_m + 0.5 * gap_m;
        let log_ratio = log_ratio_from_gap(inner_radius_m, gap_m);
        let area_m2 = product_ratio(&[core::f64::consts::TAU, mean_radius_m, gap_m], &[])?;
        let prepared = PreparedCarreau::new(law, temperature_k, maximum_laminar_shear_rate_per_s)?;
        if ![
            gap_m,
            mean_radius_m,
            log_ratio,
            prepared.ceiling.stress(),
            mean_radius_m / inner_radius_m,
            inner_radius_m / outer_radius_m,
        ]
        .into_iter()
        .all(positive)
        {
            return Err(parameters(
                "Carreau annular geometry or stress envelope is unrepresentable",
            ));
        }
        Ok(Self {
            prepared,
            inner_radius_m,
            outer_radius_m,
            length_m,
            gap_m,
            mean_radius_m,
            log_ratio,
            log_inner: det::ln(inner_radius_m),
            langevin: langevin(log_ratio),
            area_m2,
        })
    }

    pub const fn law(&self) -> &ThermalCarreauYasuda {
        &self.prepared.law
    }
    pub const fn temperature_k(&self) -> f64 {
        self.prepared.temperature.temperature_k()
    }
    pub const fn inner_radius_m(&self) -> f64 {
        self.inner_radius_m
    }
    pub const fn outer_radius_m(&self) -> f64 {
        self.outer_radius_m
    }
    pub const fn gap_m(&self) -> f64 {
        self.gap_m
    }
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub const fn maximum_laminar_shear_rate_per_s(&self) -> f64 {
        self.prepared.ceiling.rate
    }

    pub fn profile_at_pressure_and_walls(
        &self,
        pressure_drop_pa: f64,
        wall_velocities_m_s: [f64; 2],
        controls: CarreauAnnularControls,
    ) -> Result<CarreauAnnularProfile, MaterialError> {
        self.profile_with_check(pressure_drop_pa, wall_velocities_m_s, controls, || Ok(()))
    }

    /// Resolve the pressure required for signed net flow and both axial wall
    /// speeds, including blockage and reverse pumping. Every accepted force,
    /// velocity and heat value belongs to the same original nonlinear profile.
    ///
    /// The source wall-stress envelope bounds pressure corrections. A supported
    /// equal-wall-stress anchor also admits annular cases whose wall motion
    /// would exceed the inner-wall rate ceiling at zero pressure. Unsupported
    /// pressure trials only tighten the bracket; material support is never
    /// extrapolated. Newton corrections use the original positive conductance.
    /// Constant-viscosity laws use the direct analytical inverse.
    ///
    /// The flow residual is relative to max(|Q - area*mean_wall_speed|,
    /// area*|outer_wall_speed-inner_wall_speed|/2); common translation is removed.
    /// The pressure-stress correction is relative to max(|pressure_stress|,
    /// |stress_offset|), retaining a meaningful budget at zero net flow or zero
    /// pressure. Wall residual and quadrature budgets remain independent.
    ///
    /// `maximum_flow_iterations` is in 1..=128 and bounds pressure trials,
    /// in addition to one source-admitted anchor solve bounded by `controls`.
    pub fn profile_at_flow_and_walls(
        &self,
        volumetric_flow_m3_s: f64,
        wall_velocities_m_s: [f64; 2],
        controls: CarreauAnnularControls,
        relative_flow_tolerance: f64,
        maximum_flow_iterations: u32,
    ) -> Result<CarreauAnnularFlowSolution, MaterialError> {
        self.profile_at_flow_and_walls_with_check(
            volumetric_flow_m3_s,
            wall_velocities_m_s,
            controls,
            relative_flow_tolerance,
            maximum_flow_iterations,
            || Ok(()),
        )
    }

    /// Same inverse with caller-owned checks in the source admission, bounded
    /// wall/pressure corrections, quadrature and before publication. Exhaustion
    /// or cancellation returns no partial result and leaves the owner unchanged.
    pub fn profile_at_flow_and_walls_with_check(
        &self,
        volumetric_flow_m3_s: f64,
        wall_velocities_m_s: [f64; 2],
        controls: CarreauAnnularControls,
        relative_flow_tolerance: f64,
        maximum_flow_iterations: u32,
        mut check: impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<CarreauAnnularFlowSolution, MaterialError> {
        controls.admit()?;
        check()?;
        let [inner, outer] = wall_velocities_m_s;
        let relative = outer - inner;
        let mean = if inner.abs() <= outer.abs() {
            0.5_f64.mul_add(relative, inner)
        } else {
            (-0.5_f64).mul_add(relative, outer)
        };
        if [volumetric_flow_m3_s, inner, outer, relative, mean]
            .into_iter()
            .any(|x| !x.is_finite())
            || !positive(relative_flow_tolerance)
            || relative_flow_tolerance >= 1.0
            || !(1..=128).contains(&maximum_flow_iterations)
        {
            return Err(parameters(
                "invalid prescribed Carreau annular flow or inverse controls",
            ));
        }
        let entrained = self.area_m2 * mean;
        let target = volumetric_flow_m3_s - entrained;
        let sliding_scale = if relative == 0.0 {
            0.0
        } else {
            product_ratio(&[self.area_m2, relative.abs()], &[2.0])?
        };
        let flow_scale = target.abs().max(sliding_scale);
        if !entrained.is_finite() || !target.is_finite() {
            return Err(state_error(
                "prescribed Carreau annular relative flow is unrepresentable",
            ));
        }
        let accepted =
            |profile: CarreauAnnularProfile, flow_iterations| CarreauAnnularFlowSolution {
                flow_residual_m3_s: profile.flow.relative_flow_m3_s - target,
                profile,
                requested_volumetric_flow_m3_s: volumetric_flow_m3_s,
                flow_iterations,
            };
        let assess = |profile: &CarreauAnnularProfile| -> Result<(f64, f64, bool), MaterialError> {
            let difference = profile.flow.relative_flow_m3_s - target;
            if !difference.is_finite() {
                return Err(state_error(
                    "Carreau annular inverse residual is unrepresentable",
                ));
            }
            let correction = if difference == 0.0 {
                0.0
            } else {
                // An unrepresentable Newton predictor is not an extrapolated
                // material state. It falls back to the finite source bracket.
                product_ratio(
                    &[difference.abs(), self.gap_m],
                    &[
                        2.0,
                        self.length_m,
                        profile.flow.pressure_conductance_m3_s_pa,
                    ],
                )
                .unwrap_or(f64::INFINITY)
                .copysign(difference)
            };
            let stress_scale = profile
                .offset_pa
                .abs()
                .max(profile.pressure_stress_pa.abs());
            let converged = (difference == 0.0
                || (difference / flow_scale).abs() <= relative_flow_tolerance)
                && (correction == 0.0
                    || (correction / stress_scale).abs() <= relative_flow_tolerance);
            Ok((difference, correction, converged))
        };
        if self.prepared.constant() {
            let conductance = product_ratio(
                &[
                    core::f64::consts::PI,
                    self.mean_radius_m,
                    self.mean_radius_m,
                    self.gap_m,
                    self.gap_m,
                    self.langevin,
                ],
                &[2.0, self.length_m, self.prepared.rest.viscosity],
            )?;
            let couette = if relative == 0.0 {
                0.0
            } else {
                product_ratio(&[self.area_m2, self.langevin, relative.abs()], &[2.0])?
                    .copysign(relative)
            };
            let pressure_flow = target - couette;
            let pressure = if pressure_flow == 0.0 {
                0.0
            } else {
                product_ratio(&[pressure_flow.abs()], &[conductance])?.copysign(pressure_flow)
            };
            let profile =
                self.profile_with_check(pressure, wall_velocities_m_s, controls, &mut check)?;
            if !assess(&profile)?.2 {
                return Err(state_error(
                    "Carreau annular analytical flow inverse cannot resolve its tolerance",
                ));
            }
            check()?;
            return Ok(accepted(profile, u32::from(pressure != 0.0)));
        }

        let mut anchor =
            self.flow_inverse_anchor(wall_velocities_m_s, mean, relative, controls, &mut check)?;
        let (difference, correction, converged) = assess(&anchor)?;
        if converged {
            check()?;
            return Ok(accepted(anchor, 0));
        }
        let ceiling = self.prepared.ceiling.stress();
        let (mut low, mut high) = (-ceiling, ceiling);
        if difference < 0.0 {
            low = anchor.pressure_stress_pa;
        } else {
            high = anchor.pressure_stress_pa;
        }
        let bounded = |candidate: f64, low: f64, high: f64| {
            if candidate.is_finite() && candidate > low && candidate < high {
                candidate
            } else {
                0.5 * low + 0.5 * high
            }
        };
        let mut stress = bounded(anchor.pressure_stress_pa - correction, low, high);
        for iteration in 1..=maximum_flow_iterations {
            check()?;
            if stress <= low || stress >= high {
                return Err(state_error(
                    "Carreau annular inverse cannot resolve a flow inside its source envelope",
                ));
            }
            let candidate = if self
                .inverse_pressure_supported(stress, relative, &anchor, controls, &mut check)?
            {
                let pressure = if stress == 0.0 {
                    0.0
                } else {
                    product_ratio(&[2.0, self.length_m, stress.abs()], &[self.gap_m])?
                        .copysign(stress)
                };
                let profile =
                    self.profile_with_check(pressure, wall_velocities_m_s, controls, &mut check)?;
                let (difference, correction, converged) = assess(&profile)?;
                if converged {
                    check()?;
                    return Ok(accepted(profile, iteration));
                }
                if difference < 0.0 {
                    low = stress;
                } else {
                    high = stress;
                }
                anchor = profile;
                stress - correction
            } else {
                // At fixed wall speeds the admitted pressure set is an
                // interval: the two wall stresses vary monotonically in
                // opposite directions as pressure increases. Keep its known
                // supported anchor while contracting the unsupported side.
                if stress > anchor.pressure_stress_pa {
                    high = stress;
                } else {
                    low = stress;
                }
                f64::NAN
            };
            let next = bounded(candidate, low, high);
            if next == stress {
                return Err(state_error(
                    "Carreau annular inverse lost pressure-stress resolution",
                ));
            }
            stress = next;
        }
        Err(state_error(
            "Carreau annular flow inverse exhausted its iteration budget",
        ))
    }

    /// Bounded stress solve and positive quadrature with checks before every
    /// capped constitutive inverse and before publication. No allocation or
    /// mutable material history occurs; refusal leaves the owner retryable.
    pub fn profile_with_check(
        &self,
        pressure_drop_pa: f64,
        wall_velocities_m_s: [f64; 2],
        controls: CarreauAnnularControls,
        mut check: impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<CarreauAnnularProfile, MaterialError> {
        controls.admit()?;
        check()?;
        let [inner, outer] = wall_velocities_m_s;
        let relative = outer - inner;
        let mean = if inner.abs() <= outer.abs() {
            0.5_f64.mul_add(relative, inner)
        } else {
            (-0.5_f64).mul_add(relative, outer)
        };
        if [pressure_drop_pa, inner, outer, relative, mean]
            .into_iter()
            .any(|x| !x.is_finite())
        {
            return Err(parameters(
                "Carreau annular pressure or wall speeds are nonfinite",
            ));
        }
        let s = if pressure_drop_pa == 0.0 {
            0.0
        } else {
            product_ratio(&[pressure_drop_pa.abs(), self.gap_m], &[2.0, self.length_m])?
                .copysign(pressure_drop_pa)
        };
        let (minimum, maximum) = self.stress_span(s)?;
        let predicted = if relative == 0.0 {
            Ok(0.0)
        } else {
            product_ratio(
                &[self.prepared.rest.viscosity, relative.abs()],
                &[self.mean_radius_m, self.log_ratio],
            )
            .map(|x| x.copysign(relative))
        };
        let mut c = predicted
            .map_or(0.5 * minimum + 0.5 * maximum, |uniform| {
                uniform - s * self.langevin
            })
            .clamp(minimum, maximum);
        let (mut low, mut high) = (minimum, maximum);
        for iteration in 1..=controls.maximum_wall_stress_iterations {
            check()?;
            let evaluation = self.evaluate(c, s, controls, &mut check)?;
            let residual = evaluation.relative_speed - relative;
            let velocity_scale = (self.gap_m * evaluation.rate_scale).max(relative.abs());
            let correction = residual / evaluation.derivative;
            let stress_scale = c.abs().max(s.abs());
            if velocity_scale.is_finite()
                && (residual == 0.0
                    || (residual / velocity_scale).abs()
                        <= controls.relative_wall_velocity_tolerance)
                && (correction == 0.0
                    || (correction / stress_scale).abs()
                        <= controls.relative_wall_velocity_tolerance)
            {
                let iterations = if s == 0.0 && relative == 0.0 {
                    0
                } else {
                    iteration
                };
                return self.finish(
                    pressure_drop_pa,
                    wall_velocities_m_s,
                    mean,
                    relative,
                    c,
                    s,
                    evaluation,
                    iterations,
                    controls,
                    &mut check,
                );
            }
            if !residual.is_finite()
                || !velocity_scale.is_finite()
                || residual > 0.0 && c == minimum
                || residual < 0.0 && c == maximum
            {
                return Err(state_error(
                    "Carreau annular wall motion exceeds the supported stress span",
                ));
            }
            if residual < 0.0 {
                low = c;
            } else {
                high = c;
            }
            let candidate = c - correction;
            let next = if candidate.is_finite() && candidate > low && candidate < high {
                candidate
            } else {
                0.5 * low + 0.5 * high
            };
            if next == c {
                return Err(state_error(
                    "Carreau annular no-slip solve lost floating-point resolution",
                ));
            }
            c = next;
        }
        Err(state_error(
            "Carreau annular no-slip solve exhausted its iteration budget",
        ))
    }

    // The equal-wall stress amplitude c at which BOTH reconstructed wall
    // stresses stay inside the ceiling. With s=-z*c and z=h/(2*Rm) the walls
    // are (Rm/a)*(c-z*c) and (Rm/b)*(c+z*c), analytically both exactly c, so
    // the ceiling itself is the supported amplitude; the reconstruction can
    // still round one wall above it. The ceiling is taken when its walls
    // reconstruct inside it. Otherwise c is lowered ONCE by kappa*eps, where
    // kappa bounds the first-order relative error of the worse (inner) wall:
    // z carries 3 roundings (gap=b-a, mean=a+gap/2, their quotient) and z*c one,
    // and 1-z is a cancellation that turns each of those 4 into z/(1-z) of
    // them in the wall; the prefactor Rm/a, the subtraction and the final
    // product add 3 more, and the lowering itself one: kappa = 4*z/(1-z) + 4.
    // The outer wall has no cancellation (4*z/(1+z) + 3 < kappa). The check on
    // the lowered amplitude stays as verification and refuses what the bound
    // does not cover. No profile is clipped: the lowering is about 4e-14 of c
    // for the annulus a=0.001, b=0.1 (z=0.98), far below the wall-velocity
    // tolerance that closes the solve below.
    fn equal_wall_limit(&self, ceiling: f64, z: f64) -> Result<f64, MaterialError> {
        let inside = |limit: f64| {
            self.wall_stresses(limit, -z * limit)
                .into_iter()
                .all(|tau| tau.is_finite() && tau.abs() <= ceiling)
        };
        let mut limit = ceiling;
        if !inside(limit) {
            limit = ceiling * (1.0 - (4.0 * z / (1.0 - z) + 4.0) * f64::EPSILON);
        }
        if !limit.is_finite() || limit <= 0.0 || !inside(limit) {
            return Err(state_error(
                "Carreau annular equal-wall stress span is unrepresentable",
            ));
        }
        Ok(limit)
    }

    // Equal actual wall stresses t give c=t and s=-h*t/(2*Rm). The
    // cylindrical stress is positive in BOTH boundary stresses. Consequently
    // this path connects zero relative motion to the largest admissible
    // positive relative speed; sign reversal gives the negative endpoint.
    // Unlike a zero-pressure start, it cannot discard wall motion that pressure
    // can admit by redistributing the inner-cylinder shear.
    fn flow_inverse_anchor(
        &self,
        walls: [f64; 2],
        mean: f64,
        relative: f64,
        controls: CarreauAnnularControls,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<CarreauAnnularProfile, MaterialError> {
        if relative == 0.0 {
            return self.profile_with_check(0.0, walls, controls, check);
        }
        let target = relative.abs();
        let mean_rate = target / self.gap_m;
        if !positive(mean_rate) || mean_rate > self.prepared.ceiling.rate {
            return Err(state_error(
                "Carreau annular wall motion exceeds every supported profile",
            ));
        }
        let z = 0.5 * (self.gap_m / self.mean_radius_m);
        let limit = self.equal_wall_limit(self.prepared.ceiling.stress(), z)?;
        check()?;
        let mut stress = self.prepared.sample_rate(mean_rate)?.stress().min(limit);
        let (mut low, mut high) = (0.0, limit);
        for iteration in 1..=controls.maximum_wall_stress_iterations {
            check()?;
            let mut evaluation = self.evaluate(stress, -z * stress, controls, check)?;
            let difference = evaluation.relative_speed - target;
            let derivative = product_ratio(
                &[
                    evaluation.derivative,
                    1.0 - z * (1.0 - 2.0 * evaluation.moments.mean),
                ],
                &[],
            )?;
            let correction = difference / derivative;
            let scale = target.max(self.gap_m * evaluation.rate_scale);
            if scale.is_finite()
                && (difference == 0.0
                    || (difference / scale).abs() <= controls.relative_wall_velocity_tolerance)
                && (correction == 0.0
                    || (correction / stress).abs() <= controls.relative_wall_velocity_tolerance)
            {
                let c = stress.copysign(relative);
                let s = -z * c;
                let pressure =
                    product_ratio(&[2.0, self.length_m, s.abs()], &[self.gap_m])?.copysign(s);
                evaluation.relative_speed = evaluation.relative_speed.copysign(relative);
                return self.finish(
                    pressure, walls, mean, relative, c, s, evaluation, iteration, controls, check,
                );
            }
            if !difference.is_finite() || !scale.is_finite() || difference < 0.0 && stress == limit
            {
                return Err(state_error(
                    "Carreau annular wall motion exceeds every supported profile",
                ));
            }
            if difference < 0.0 {
                low = stress;
            } else {
                high = stress;
            }
            let candidate = stress - correction;
            let next = if candidate.is_finite() && candidate > low && candidate < high {
                candidate
            } else {
                0.5 * low + 0.5 * high
            };
            if next == stress || next <= low || next >= high {
                return Err(state_error(
                    "Carreau annular source anchor lost stress resolution",
                ));
            }
            stress = next;
        }
        Err(state_error(
            "Carreau annular source anchor exhausted its iteration budget",
        ))
    }

    fn inverse_pressure_supported(
        &self,
        s: f64,
        relative: f64,
        anchor: &CarreauAnnularProfile,
        controls: CarreauAnnularControls,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<bool, MaterialError> {
        let (minimum, maximum) = self.stress_span(s)?;
        let distance = (s - anchor.pressure_stress_pa).abs();
        // dc/ds=-(1-2*weighted_mean_area) lies in [-1,1]. These bounds
        // avoid endpoint quadrature when the previous accepted stress and its
        // signed wall residual already prove admission. Otherwise evaluate the
        // original source boundary; no failed forward solve is swallowed.
        let lower_proven =
            minimum <= anchor.offset_pa - distance && anchor.flow.wall_velocity_residual_m_s <= 0.0;
        if !lower_proven && self.evaluate(minimum, s, controls, check)?.relative_speed > relative {
            return Ok(false);
        }
        let upper_proven =
            maximum >= anchor.offset_pa + distance && anchor.flow.wall_velocity_residual_m_s >= 0.0;
        if !upper_proven && self.evaluate(maximum, s, controls, check)?.relative_speed < relative {
            return Ok(false);
        }
        Ok(true)
    }

    // c=(a*tau(a)+b*tau(b))/(2 Rm), s=delta_p*h/(2L). Maximum |tau|
    // occurs at a wall: an interior extremum of C/r-G*r is a minimum of |tau|.
    fn stress_span(&self, s: f64) -> Result<(f64, f64), MaterialError> {
        let limit = self.prepared.ceiling.stress();
        if !s.is_finite() || s.abs() > limit {
            return Err(state_error(
                "Carreau annular pressure exceeds the supported stress span",
            ));
        }
        let inner = limit * (self.inner_radius_m / self.mean_radius_m);
        let outer = limit * (self.outer_radius_m / self.mean_radius_m);
        let mut low = (-inner - s).max(-outer + s);
        let mut high = (inner - s).min(outer + s);
        // Round only inward if reconstructing an endpoint stress would cross
        // the original envelope. No material sample is clipped or extrapolated.
        for _ in 0..8 {
            if self
                .wall_stresses(low, s)
                .into_iter()
                .all(|x| x.is_finite() && x.abs() <= limit)
            {
                break;
            }
            low = low.next_up();
        }
        for _ in 0..8 {
            if self
                .wall_stresses(high, s)
                .into_iter()
                .all(|x| x.is_finite() && x.abs() <= limit)
            {
                break;
            }
            high = high.next_down();
        }
        if !low.is_finite()
            || !high.is_finite()
            || low > high
            || [low, high].into_iter().any(|c| {
                self.wall_stresses(c, s)
                    .into_iter()
                    .any(|x| !x.is_finite() || x.abs() > limit)
            })
        {
            return Err(state_error(
                "Carreau annulus has no representable supported wall-stress span",
            ));
        }
        Ok((low, high))
    }

    fn wall_stresses(&self, c: f64, s: f64) -> [f64; 2] {
        [
            (self.mean_radius_m / self.inner_radius_m) * (c + s),
            (self.mean_radius_m / self.outer_radius_m) * (c - s),
        ]
    }

    // x=ln(r/a)/ln(b/a); y=(r²-a²)/(b²-a²). The latter is evaluated
    // from the physical gap, retaining the slit limit without squared radii.
    fn point(&self, x: f64) -> (f64, f64) {
        if x == 0.0 {
            return (self.inner_radius_m, 0.0);
        }
        if x == 1.0 {
            return (self.outer_radius_m, 1.0);
        }
        let t = if self.log_ratio < 0.5 {
            det::expm1(self.log_ratio * x) / (self.gap_m / self.inner_radius_m)
        } else {
            let r = det::exp(self.log_inner + self.log_ratio * x)
                .clamp(self.inner_radius_m, self.outer_radius_m);
            (r - self.inner_radius_m) / self.gap_m
        };
        let radius = self.gap_m.mul_add(t, self.inner_radius_m);
        let y = t
            * (self.inner_radius_m / self.mean_radius_m
                + 0.5 * (self.gap_m / self.mean_radius_m) * t);
        (radius, y.clamp(0.0, 1.0))
    }

    fn stress(&self, x: f64, c: f64, s: f64) -> f64 {
        let (radius, y) = self.point(x);
        (self.mean_radius_m / radius) * s.mul_add(1.0 - 2.0 * y, c)
    }

    fn coordinate_at_area(&self, y: f64) -> f64 {
        let a = self.inner_radius_m / self.mean_radius_m;
        let e = self.gap_m / self.mean_radius_m;
        let t = 2.0 * y / (a + det::sqrt(a * a + 2.0 * e * y));
        log_ratio_from_gap(self.inner_radius_m, self.gap_m * t) / self.log_ratio
    }

    fn cuts(
        &self,
        c: f64,
        s: f64,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<Cuts, MaterialError> {
        let mut cuts = Cuts {
            points: [0.0; 12],
            count: 2,
        };
        cuts.points[1] = 1.0;
        cuts.push(self.coordinate_at_area(0.5))?;
        if s != 0.0 {
            let zero = 0.5 * ((c + s) / s);
            let a = self.inner_radius_m / self.mean_radius_m;
            let extremum = -zero - a * (self.inner_radius_m / self.gap_m);
            for y in [zero, extremum] {
                if y > 0.0 && y < 1.0 {
                    cuts.push(self.coordinate_at_area(y))?;
                }
            }
        }
        cuts.sort();
        if let Some(rate) = self.prepared.transition_rate(self.prepared.ceiling.rate) {
            let transition = self.prepared.sample_rate(rate)?.stress();
            let original = cuts;
            for pair in original.points[..original.count].windows(2) {
                let left_value = self.stress(pair[0], c, s);
                let right_value = self.stress(pair[1], c, s);
                for target in [-transition, transition] {
                    if target <= left_value.min(right_value)
                        || target >= left_value.max(right_value)
                    {
                        continue;
                    }
                    let (mut low, mut high) = (pair[0], pair[1]);
                    for _ in 0..64 {
                        check()?;
                        let mid = 0.5 * low + 0.5 * high;
                        if mid == low || mid == high {
                            break;
                        }
                        if (self.stress(mid, c, s) < target) == (left_value < right_value) {
                            low = mid;
                        } else {
                            high = mid;
                        }
                    }
                    cuts.push(0.5 * low + 0.5 * high)?;
                }
            }
            cuts.sort();
        }
        Ok(cuts)
    }

    fn local(
        &self,
        stress: f64,
        controls: CarreauPipeControls,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<RateState, MaterialError> {
        check()?;
        Ok(self
            .prepared
            .inverse(stress.abs(), self.prepared.ceiling, controls)?
            .0)
    }

    // Every part has constant sign. Positive canonical quadrature therefore
    // controls the underlying magnitudes even when the total signed load is zero.
    #[allow(clippy::too_many_arguments)]
    fn integral(
        &self,
        c: f64,
        s: f64,
        cuts: Cuts,
        rate_scale: f64,
        stress_scale: f64,
        kind: Integral,
        lower: f64,
        upper: f64,
        controls: CarreauPipeControls,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<(f64, f64), MaterialError> {
        if rate_scale == 0.0 || lower == upper {
            return Ok((0.0, 0.0));
        }
        let mut total = 0.0;
        let mut compensation = 0.0;
        let mut change = 0.0_f64;
        for pair in cuts.points[..cuts.count].windows(2) {
            let (lo, hi) = (pair[0].max(lower), pair[1].min(upper));
            if lo >= hi {
                continue;
            }
            let mid = 0.5 * lo + 0.5 * hi;
            let sign = match kind {
                Integral::Heat => 1.0,
                Integral::Velocity => self.stress(mid, c, s).signum(),
                Integral::Flow => {
                    self.stress(mid, c, s).signum() * (1.0 - 2.0 * self.point(mid).1).signum()
                }
            };
            let (value, estimate) = integrate(
                &mut |t| {
                    let x = (hi - lo).mul_add(t, lo);
                    let (r, y) = self.point(x);
                    let stress = self.stress(x, c, s);
                    let rate = self.local(stress, controls, check)?.rate / rate_scale;
                    let radius = r / self.outer_radius_m;
                    let factor = match kind {
                        Integral::Velocity => radius,
                        Integral::Flow => radius * (1.0 - 2.0 * y).abs(),
                        Integral::Heat => radius * radius * (stress.abs() / stress_scale),
                    };
                    Ok(factor * rate)
                },
                None,
                controls,
            )?;
            let corrected = sign * (hi - lo) * value - compensation;
            let next = total + corrected;
            compensation = (next - total) - corrected;
            total = next;
            change = change.max(estimate);
        }
        Ok((total, change))
    }

    fn evaluate(
        &self,
        c: f64,
        s: f64,
        controls: CarreauAnnularControls,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<Evaluation, MaterialError> {
        let stresses = self.wall_stresses(c, s);
        let inner = self.local(stresses[0], controls.local, check)?;
        let outer = self.local(stresses[1], controls.local, check)?;
        let rate_scale = inner.rate.max(outer.rate);
        let stress_scale = stresses[0].abs().max(stresses[1].abs());
        let tangent_scale = inner.tangent.min(outer.tangent);
        let cuts = self.cuts(c, s, check)?;
        let (relative_speed, moments, change) = if self.prepared.constant() || rate_scale == 0.0 {
            let numerator = s.mul_add(self.langevin, c);
            let relative = if numerator == 0.0 {
                0.0
            } else {
                product_ratio(
                    &[self.mean_radius_m, self.log_ratio, numerator.abs()],
                    &[self.prepared.rest.viscosity],
                )?
                .copysign(numerator)
            };
            let mean = 0.5 * (1.0 - self.langevin);
            let variance = self.langevin / (4.0 * self.log_ratio);
            (
                relative,
                Moments {
                    weight: 1.0,
                    mean,
                    variance,
                    pressure: mean * (1.0 - mean) - variance,
                },
                0.0,
            )
        } else {
            let (value, mut change) = self.integral(
                c,
                s,
                cuts,
                rate_scale,
                stress_scale,
                Integral::Velocity,
                0.0,
                1.0,
                controls.local,
                check,
            )?;
            let relative = if value == 0.0 {
                0.0
            } else {
                product_ratio(
                    &[self.log_ratio, self.outer_radius_m, rate_scale, value.abs()],
                    &[],
                )?
                .copysign(value)
            };
            let mut moments = Moments::default();
            for pair in cuts.points[..cuts.count].windows(2) {
                if pair[0] == pair[1] {
                    continue;
                }
                let (part, estimate) = integrate_moments(
                    &mut |x| {
                        let y = self.point(x).1;
                        let local = self.local(self.stress(x, c, s), controls.local, check)?;
                        let weight = tangent_scale / local.tangent;
                        Ok(Moments {
                            weight,
                            mean: y,
                            variance: 0.0,
                            pressure: y * (1.0 - y) * weight,
                        })
                    },
                    pair[0],
                    pair[1],
                    controls.local,
                )?;
                moments = moments.combine(part);
                change = change.max(estimate);
            }
            (relative, moments, change)
        };
        let derivative = product_ratio(
            &[self.mean_radius_m, self.log_ratio, moments.weight],
            &[tangent_scale],
        )?;
        Ok(Evaluation {
            relative_speed,
            derivative,
            moments,
            tangent_scale,
            rate_scale,
            stress_scale,
            cuts,
            change,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        pressure: f64,
        requested_walls: [f64; 2],
        mean: f64,
        requested_relative: f64,
        c: f64,
        s: f64,
        evaluation: Evaluation,
        iterations: u32,
        controls: CarreauAnnularControls,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<CarreauAnnularProfile, MaterialError> {
        let Evaluation {
            relative_speed,
            derivative,
            moments,
            tangent_scale,
            rate_scale,
            stress_scale,
            cuts,
            mut change,
        } = evaluation;
        let conductance = product_ratio(
            &[
                core::f64::consts::TAU,
                self.mean_radius_m,
                self.mean_radius_m,
                self.gap_m,
                self.gap_m,
                self.log_ratio,
                moments.variance,
            ],
            &[self.length_m, tangent_scale],
        )?;
        let (relative_flow, heat) = if self.prepared.constant() || rate_scale == 0.0 {
            let pressure_flow = if pressure == 0.0 {
                0.0
            } else {
                product_ratio(&[conductance, pressure.abs()], &[])?.copysign(pressure)
            };
            let sliding_flow = if relative_speed == 0.0 {
                0.0
            } else {
                product_ratio(&[self.area_m2, self.langevin, relative_speed.abs()], &[2.0])?
                    .copysign(relative_speed)
            };
            let sliding_heat = if relative_speed == 0.0 {
                0.0
            } else {
                product_ratio(
                    &[
                        core::f64::consts::TAU,
                        self.length_m,
                        self.prepared.rest.viscosity,
                        relative_speed.abs(),
                        relative_speed.abs(),
                    ],
                    &[self.log_ratio],
                )?
            };
            (
                pressure_flow + sliding_flow,
                pressure * pressure_flow + sliding_heat,
            )
        } else {
            let (value, flow_change) = self.integral(
                c,
                s,
                cuts,
                rate_scale,
                stress_scale,
                Integral::Flow,
                0.0,
                1.0,
                controls.local,
                check,
            )?;
            let flow = if value == 0.0 {
                0.0
            } else {
                product_ratio(
                    &[
                        core::f64::consts::PI,
                        self.mean_radius_m,
                        self.gap_m,
                        self.log_ratio,
                        self.outer_radius_m,
                        rate_scale,
                        value.abs(),
                    ],
                    &[],
                )?
                .copysign(value)
            };
            let (value, heat_change) = self.integral(
                c,
                s,
                cuts,
                rate_scale,
                stress_scale,
                Integral::Heat,
                0.0,
                1.0,
                controls.local,
                check,
            )?;
            let heat = product_ratio(
                &[
                    core::f64::consts::TAU,
                    self.length_m,
                    self.log_ratio,
                    self.outer_radius_m,
                    self.outer_radius_m,
                    stress_scale,
                    rate_scale,
                    value,
                ],
                &[],
            )?;
            change = change.max(flow_change).max(heat_change);
            (flow, heat)
        };
        let force = |stress: f64| -> Result<f64, MaterialError> {
            if stress == 0.0 {
                Ok(0.0)
            } else {
                Ok(product_ratio(
                    &[
                        core::f64::consts::TAU,
                        self.length_m,
                        self.mean_radius_m,
                        stress.abs(),
                    ],
                    &[],
                )?
                .copysign(stress))
            }
        };
        let relative_force = force(c)?;
        let entrained = self.area_m2 * mean;
        let flow = entrained + relative_flow;
        let wall_power = relative_force * requested_relative - pressure * entrained;
        let work_residual =
            pressure.mul_add(relative_flow, relative_force * requested_relative) - heat;
        let stresses = self.wall_stresses(c, s);
        let wall_forces = [-force(c + s)?, force(c - s)?];
        let accepted_walls = [mean - 0.5 * relative_speed, mean + 0.5 * relative_speed];
        let residual = relative_speed - requested_relative;
        let pressure_power = pressure * flow;
        let cross = 0.5
            * self.area_m2
            * if self.prepared.constant() || rate_scale == 0.0 {
                self.langevin
            } else {
                1.0 - 2.0 * moments.mean
            };
        let resistance = product_ratio(
            &[core::f64::consts::TAU, self.length_m, self.mean_radius_m],
            &[derivative],
        )?;
        if [
            entrained,
            flow,
            relative_flow,
            heat,
            wall_power,
            work_residual,
            residual,
            pressure_power,
            cross,
            accepted_walls[0],
            accepted_walls[1],
            wall_forces[0],
            wall_forces[1],
        ]
        .into_iter()
        .any(|x| !x.is_finite())
            || heat < 0.0
        {
            return Err(state_error(
                "Carreau annular flux, force or physical work is unrepresentable",
            ));
        }
        check()?;
        Ok(CarreauAnnularProfile {
            annulus: *self,
            controls,
            offset_pa: c,
            pressure_stress_pa: s,
            cuts,
            rate_scale,
            stress_scale,
            flow: CarreauAnnularFlow {
                pressure_drop_pa: pressure,
                requested_wall_velocities_m_s: requested_walls,
                accepted_wall_velocities_m_s: accepted_walls,
                wall_velocity_residual_m_s: residual,
                volumetric_flow_m3_s: flow,
                mean_wall_entrainment_m3_s: entrained,
                relative_flow_m3_s: relative_flow,
                wall_tractions_pa: [-stresses[0], stresses[1]],
                wall_forces_n: wall_forces,
                pressure_power_w: pressure_power,
                wall_power_w: wall_power,
                dissipated_power_w: heat,
                work_residual_w: work_residual,
                pressure_conductance_m3_s_pa: conductance,
                flow_relative_wall_velocity_derivative_m2: cross,
                relative_wall_resistance_n_s_m: resistance,
                maximum_shear_rate_per_s: rate_scale,
                quadrature_relative_change: change,
                wall_stress_iterations: iterations,
            },
        })
    }
}

impl CarreauAnnularProfile {
    pub const fn annulus(&self) -> &LaminarCarreauAnnulus {
        &self.annulus
    }
    pub const fn flow(&self) -> &CarreauAnnularFlow {
        &self.flow
    }
    /// Mean radius-normalized weighted wall stress, not a planar midpoint stress.
    pub const fn stress_offset_pa(&self) -> f64 {
        self.offset_pa
    }

    pub fn sample(&self, radius_m: f64) -> Result<CarreauAnnularSample, MaterialError> {
        if !radius_m.is_finite()
            || radius_m < self.annulus.inner_radius_m
            || radius_m > self.annulus.outer_radius_m
        {
            return Err(parameters(
                "Carreau annular sample must lie between the physical cylinders",
            ));
        }
        self.sample_gap_fraction((radius_m - self.annulus.inner_radius_m) / self.annulus.gap_m)
    }

    /// Fractional gap position in [0,1], retaining thin-gap locations even when
    /// adding their distance to a large absolute radius would round it away.
    pub fn sample_gap_fraction(
        &self,
        fraction: f64,
    ) -> Result<CarreauAnnularSample, MaterialError> {
        self.sample_gap_fraction_with_check(fraction, || Ok(()))
    }

    pub fn sample_gap_fraction_with_check(
        &self,
        fraction: f64,
        mut check: impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<CarreauAnnularSample, MaterialError> {
        check()?;
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err(parameters("Carreau annular gap fraction must lie in [0,1]"));
        }
        let owner = &self.annulus;
        let x = if fraction == 0.0 {
            0.0
        } else if fraction == 1.0 {
            1.0
        } else {
            log_ratio_from_gap(owner.inner_radius_m, owner.gap_m * fraction) / owner.log_ratio
        };
        let (radius, y) = if fraction == 0.0 {
            (owner.inner_radius_m, 0.0)
        } else if fraction == 1.0 {
            (owner.outer_radius_m, 1.0)
        } else {
            (
                owner.gap_m.mul_add(fraction, owner.inner_radius_m),
                fraction
                    * (owner.inner_radius_m / owner.mean_radius_m
                        + 0.5 * (owner.gap_m / owner.mean_radius_m) * fraction),
            )
        };
        let stress = (owner.mean_radius_m / radius)
            * self
                .pressure_stress_pa
                .mul_add(1.0 - 2.0 * y, self.offset_pa);
        let local = owner.local(stress, self.controls.local, &mut check)?;
        let (velocity, change) = if fraction == 0.0 {
            (self.flow.accepted_wall_velocities_m_s[0], 0.0)
        } else if fraction == 1.0 {
            (self.flow.accepted_wall_velocities_m_s[1], 0.0)
        } else {
            let (lo, hi, wall, direction) = if x <= 0.5 {
                (0.0, x, self.flow.accepted_wall_velocities_m_s[0], 1.0)
            } else {
                (x, 1.0, self.flow.accepted_wall_velocities_m_s[1], -1.0)
            };
            let (value, change) = owner.integral(
                self.offset_pa,
                self.pressure_stress_pa,
                self.cuts,
                self.rate_scale,
                self.stress_scale,
                Integral::Velocity,
                lo,
                hi,
                self.controls.local,
                &mut check,
            )?;
            let relative = if value == 0.0 {
                0.0
            } else {
                product_ratio(
                    &[
                        owner.log_ratio,
                        owner.outer_radius_m,
                        self.rate_scale,
                        value.abs(),
                    ],
                    &[],
                )?
                .copysign(value)
            };
            (wall + direction * relative, change)
        };
        let heat = if local.rate == 0.0 {
            0.0
        } else {
            product_ratio(&[stress.abs(), local.rate], &[])?
        };
        if !velocity.is_finite() {
            return Err(state_error(
                "Carreau annular sampled velocity is unrepresentable",
            ));
        }
        check()?;
        Ok(CarreauAnnularSample {
            axial_velocity_m_s: velocity,
            shear_stress_pa: stress,
            radial_velocity_gradient_per_s: local.rate.copysign(stress),
            dynamic_viscosity_pa_s: local.viscosity,
            dissipation_w_m3: heat,
            quadrature_relative_change: change,
        })
    }
}

fn log_ratio_from_gap(radius: f64, gap: f64) -> f64 {
    let relative = gap / radius;
    if relative.is_finite() {
        let sum = 1.0 + relative;
        if sum == 1.0 {
            relative
        } else if relative < 0.5 {
            det::ln(sum) * (relative / (sum - 1.0))
        } else {
            det::ln(sum)
        }
    } else {
        det::ln(radius + gap) - det::ln(radius)
    }
}

// coth(x)-1/x, the annular Couette volume asymmetry and Poiseuille factor.
// Its small-gap series avoids subtracting two terms of order 1/x.
fn langevin(x: f64) -> f64 {
    if x < 0.125 {
        let z = x * x;
        x * (1.0 / 3.0
            + z * (-1.0 / 45.0
                + z * (2.0 / 945.0
                    + z * (-1.0 / 4725.0 + z * (2.0 / 93555.0 - z * (1382.0 / 638512875.0))))))
    } else {
        1.0 - x.recip() + 2.0 / det::expm1(2.0 * x)
    }
}
