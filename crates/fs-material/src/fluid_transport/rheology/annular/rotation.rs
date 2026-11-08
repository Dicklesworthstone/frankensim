//! Steady circular Couette flow of the original thermal Carreau-Yasuda law.
//!
//! For azimuthal velocity v=r*omega(r), the physical shear is
//! gamma=r*domega/dr=dv/dr-v/r, not dv/dr. Angular momentum fixes
//! r^2*tau=C, with tau=mu(|gamma|)*gamma. Thus rigid rotation has no
//! viscous loss and both cylindrical walls exchange the same torque.
//! The Newtonian solution and torque follow Fitzpatrick, Fluid Mechanics,
//! section 10.4: <https://farside.ph.utexas.edu/teaching/336L/Fluidhtml/node137.html>.
//!
//! The original calibrated source is inverted at the actual curved shear
//! stress. Positive quadrature in logarithmic radius resolves the speed
//! difference; the local differential stress supplies the torque derivative.
//! Preparation and bounded queries allocate no scratch. Constant viscosity
//! takes the analytic path, including the narrow-gap limit.
//!
//! This is the fully developed laminar base state at fixed temperature, not
//! Taylor-vortex stability, turbulent torque, axial/rotational superposition,
//! end effects, slip, or a thermal evolution. Density and the radial pressure
//! balance are not supplied. A caller must qualify the laminar regime and
//! deposit the reported heat through its thermal owner.

use super::super::pipe::{RateState, integrate, product_ratio};
use super::{CarreauAnnularControls, LaminarCarreauAnnulus, MaterialError};
use super::{log_ratio_from_gap, parameters, state_error};
use fs_math::det;

/// Wall ports and irreversible power for one accepted circular Couette state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauRotation {
    /// Inner/outer prescribed angular speeds about the same axis [rad/s].
    pub requested_wall_angular_velocities_rad_s: [f64; 2],
    /// Inner speed and the outer speed reached by the accepted profile.
    pub accepted_wall_angular_velocities_rad_s: [f64; 2],
    /// Actual outer-minus-inner speed, retained before adding common rotation.
    pub relative_angular_velocity_rad_s: f64,
    /// Accepted minus requested outer-minus-inner speed [rad/s].
    pub angular_velocity_residual_rad_s: f64,
    /// Azimuthal tractions on the fluid, using its inner/outer outward normals.
    pub wall_tractions_pa: [f64; 2],
    /// Torques exerted by the inner/outer walls ON THE FLUID [N m].
    pub wall_torques_n_m: [f64; 2],
    /// Work rate at the requested wall speeds, formed without common-spin cancellation.
    pub wall_power_w: f64,
    /// Integral of the original positive tau*gamma over the annular volume.
    pub dissipated_power_w: f64,
    /// Requested wall work minus accepted-profile heat [W].
    pub work_residual_w: f64,
    /// d(outer wall torque)/d(outer-minus-inner angular speed) [N m s].
    pub differential_resistance_n_m_s: f64,
    /// Largest local shear magnitude, attained at the inner cylinder [1/s].
    pub maximum_shear_rate_per_s: f64,
    /// Observed positive-integral refinement change; not a certified error bound.
    pub quadrature_relative_change: f64,
    /// Bounded outer solve iterations; zero for exact constant/rest limits.
    pub wall_stress_iterations: u32,
}

/// Accepted source, geometry and rotation. Sampling retains the same work bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauRotationProfile {
    annulus: LaminarCarreauAnnulus,
    controls: CarreauAnnularControls,
    rotation: CarreauRotation,
    inner: RateState,
    sign: f64,
    transition: Option<f64>,
}

/// Actual cylindrical shear and velocity at a physical radial position.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauRotationSample {
    /// Angular velocity [rad/s].
    pub angular_velocity_rad_s: f64,
    /// Azimuthal velocity r*omega [m/s].
    pub azimuthal_velocity_m_s: f64,
    /// Signed r*domega/dr [1/s], excluding rigid spin.
    pub shear_rate_per_s: f64,
    /// tau_rtheta on the increasing-radius face [Pa].
    pub radial_shear_stress_pa: f64,
    /// Original constitutive viscosity at the actual local rate [Pa s].
    pub dynamic_viscosity_pa_s: f64,
    /// Positive local irreversible power per volume [W/m^3].
    pub dissipation_w_m3: f64,
    /// Observed angular-profile quadrature change.
    pub quadrature_relative_change: f64,
}

struct Evaluation {
    speed: f64,
    // d(outer-minus-inner angular speed)/d(inner-wall stress).
    mobility: f64,
    change: f64,
    transition: Option<f64>,
}

impl LaminarCarreauAnnulus {
    /// Resolve differential rotation without substituting planar wall shear.
    /// Axial flow and rotation are separate reductions; nonlinear profiles
    /// cannot be superposed. Work and heat use the same accepted local law.
    pub fn rotation_at_wall_speeds(
        &self,
        wall_angular_velocities_rad_s: [f64; 2],
        controls: CarreauAnnularControls,
    ) -> Result<CarreauRotationProfile, MaterialError> {
        self.rotation_with_check(wall_angular_velocities_rad_s, controls, || Ok(()))
    }

    /// Checked counterpart, polling before local inverses and quadrature nodes
    /// and before publishing a result. Cancellation never returns a partial profile.
    pub fn rotation_with_check(
        &self,
        walls: [f64; 2],
        controls: CarreauAnnularControls,
        mut check: impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<CarreauRotationProfile, MaterialError> {
        controls.admit()?;
        check()?;
        let relative = walls[1] - walls[0];
        if walls.iter().any(|v| !v.is_finite()) || !relative.is_finite() {
            return Err(parameters(
                "Carreau rotation requires finite wall angular speeds",
            ));
        }
        let target = relative.abs();
        if target == 0.0 {
            let inner = self.prepared.rest;
            let evaluation = self.evaluate_rotation(inner, controls, &mut check)?;
            return self.finish_rotation(walls, inner, evaluation, controls, 0, &mut check);
        }
        let maximum = self.evaluate_rotation(self.prepared.ceiling, controls, &mut check)?;
        if target > maximum.speed {
            return Err(state_error(
                "Carreau rotation exceeds the supported inner-wall shear rate",
            ));
        }
        if target == maximum.speed {
            return self.finish_rotation(
                walls,
                self.prepared.ceiling,
                maximum,
                controls,
                0,
                &mut check,
            );
        }
        // Admit the original ceiling before inverting even a constant law:
        // a supported endpoint can otherwise round-trip one ULP above its rate.
        // 1-(a/b)^2 is formed from the gap, never by subtracting squared radii.
        let span = self.rotation_stress_span();
        let seed = product_ratio(&[2.0, target], &[span])?;
        if self.prepared.constant() {
            let inner = self
                .prepared
                .sample_rate(seed.min(self.prepared.ceiling.rate))?;
            let evaluation = self.evaluate_rotation(inner, controls, &mut check)?;
            return self.finish_rotation(walls, inner, evaluation, controls, 0, &mut check);
        }
        let (mut low, mut high) = (0.0, self.prepared.ceiling.rate);
        let mut rate = if seed < high { seed } else { 0.5 * high };
        for iteration in 1..=controls.maximum_wall_stress_iterations {
            check()?;
            let inner = self.prepared.sample_rate(rate)?;
            let evaluation = self.evaluate_rotation(inner, controls, &mut check)?;
            let residual = evaluation.speed - target;
            let correction = (residual / evaluation.mobility) / inner.tangent;
            let tolerance = controls.relative_wall_velocity_tolerance;
            if (residual / target).abs() <= tolerance && (correction / rate).abs() <= tolerance {
                return self
                    .finish_rotation(walls, inner, evaluation, controls, iteration, &mut check);
            }
            if residual < 0.0 {
                low = rate;
            } else {
                high = rate;
            }
            let newton = rate - correction;
            rate = if newton > low && newton < high {
                newton
            } else {
                0.5 * low + 0.5 * high
            };
            if rate == low || rate == high {
                return Err(state_error(
                    "Carreau rotation lost inner-wall rate resolution",
                ));
            }
        }
        Err(state_error(
            "Carreau rotation exhausted its wall-stress iteration budget",
        ))
    }

    fn rotation_stress_span(&self) -> f64 {
        (self.gap_m / self.outer_radius_m) * (1.0 + self.inner_radius_m / self.outer_radius_m)
    }

    fn rotation_local(
        &self,
        x: f64,
        inner: RateState,
        controls: CarreauAnnularControls,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<RateState, MaterialError> {
        check()?;
        if x == 0.0 || inner.rate == 0.0 {
            return Ok(inner);
        }
        let stress = self.rotation_stress(x, inner)?;
        Ok(self.prepared.inverse(stress, inner, controls.local)?.0)
    }

    fn rotation_stress(&self, x: f64, inner: RateState) -> Result<f64, MaterialError> {
        if x == 0.0 || inner.rate == 0.0 {
            return Ok(inner.stress());
        }
        let r = self.point(x).0;
        let stress = product_ratio(
            &[inner.stress(), self.inner_radius_m, self.inner_radius_m],
            &[r, r],
        )?;
        // The exact cylindrical stress decreases monotonically from the inner
        // wall. Bound only a last-bit product roundoff at that endpoint.
        Ok(stress.min(inner.stress()))
    }

    fn rotation_integral(
        &self,
        inner: RateState,
        bounds: [f64; 2],
        transition: Option<f64>,
        controls: CarreauAnnularControls,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<(f64, f64), MaterialError> {
        let [lower, upper] = bounds;
        if inner.rate == 0.0 || lower == upper {
            return Ok((0.0, 0.0));
        }
        let width = upper - lower;
        if self.prepared.constant() {
            // Integral gamma_i*exp(-2*l*x)*l dx; expm1 retains tiny gaps.
            let fraction = -det::expm1(-2.0 * self.log_ratio * width);
            let value = product_ratio(
                &[
                    0.5,
                    inner.rate,
                    det::exp(-2.0 * self.log_ratio * lower),
                    fraction,
                ],
                &[],
            )?;
            return Ok((value, 0.0));
        }
        let split = transition.map(|x| (x - lower) / width);
        let (value, change) = integrate(
            &mut |t| {
                let x = width.mul_add(t, lower);
                Ok(self.rotation_local(x, inner, controls, check)?.rate / inner.rate)
            },
            split,
            controls.local,
        )?;
        Ok((
            product_ratio(&[self.log_ratio, width, inner.rate, value], &[])?,
            change,
        ))
    }

    fn evaluate_rotation(
        &self,
        inner: RateState,
        controls: CarreauAnnularControls,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<Evaluation, MaterialError> {
        if inner.rate == 0.0 || self.prepared.constant() {
            let span = self.rotation_stress_span();
            let speed = if inner.rate == 0.0 {
                0.0
            } else {
                product_ratio(&[0.5, span, inner.rate], &[])?
            };
            return Ok(Evaluation {
                speed,
                mobility: product_ratio(&[0.5, span], &[inner.viscosity])?,
                change: 0.0,
                transition: None,
            });
        }
        let transition = if let Some(rate) = self.prepared.transition_rate(inner.rate) {
            let stress = self.prepared.sample_rate(rate)?.stress();
            let log = det::ln(inner.stress()) - det::ln(stress);
            Some(0.5 * log / self.log_ratio).filter(|x| *x > 0.0 && *x < 1.0)
        } else {
            None
        };
        let (speed, mut change) =
            self.rotation_integral(inner, [0.0, 1.0], transition, controls, check)?;
        let outer = self.rotation_local(1.0, inner, controls, check)?;
        let difference = inner.rate - outer.rate;
        // Leibniz on 1/2 int gamma(tau)/tau d(tau) avoids a second integral
        // except in narrow gaps, where subtracting nearly equal rates is ill-conditioned.
        let mobility = if difference > 0.01 * inner.rate {
            product_ratio(&[0.5, difference], &[inner.stress()])?
        } else {
            let (value, estimate) = integrate(
                &mut |x| {
                    let local = self.rotation_local(x, inner, controls, check)?;
                    let r = self.point(x).0;
                    product_ratio(
                        &[self.inner_radius_m, self.inner_radius_m, inner.tangent],
                        &[r, r, local.tangent],
                    )
                },
                transition,
                controls.local,
            )?;
            change = change.max(estimate);
            product_ratio(&[self.log_ratio, value], &[inner.tangent])?
        };
        Ok(Evaluation {
            speed,
            mobility,
            change,
            transition,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn finish_rotation(
        &self,
        walls: [f64; 2],
        inner: RateState,
        evaluation: Evaluation,
        controls: CarreauAnnularControls,
        iterations: u32,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<CarreauRotationProfile, MaterialError> {
        let requested = walls[1] - walls[0];
        let sign = requested.signum();
        let relative = evaluation.speed.copysign(requested);
        let residual = relative - requested;
        let outer_stress = self.rotation_stress(1.0, inner)?;
        let torque_scale = [
            core::f64::consts::TAU,
            self.length_m,
            self.inner_radius_m,
            self.inner_radius_m,
        ];
        let resistance = product_ratio(&torque_scale, &[evaluation.mobility])?;
        let torque = if inner.rate == 0.0 {
            0.0
        } else {
            product_ratio(
                &[
                    torque_scale[0],
                    torque_scale[1],
                    torque_scale[2],
                    torque_scale[3],
                    inner.stress(),
                ],
                &[],
            )?
            .copysign(requested)
        };
        let heat = if torque == 0.0 {
            0.0
        } else {
            product_ratio(&[torque.abs(), evaluation.speed], &[])?
        };
        let work = if torque == 0.0 {
            0.0
        } else {
            product_ratio(&[torque.abs(), requested.abs()], &[])?
        };
        let work_residual = -torque * residual;
        let outer_speed = walls[0] + relative;
        if !outer_speed.is_finite() || !work_residual.is_finite() {
            return Err(state_error(
                "Carreau rotation wall speed or work is unrepresentable",
            ));
        }
        check()?;
        Ok(CarreauRotationProfile {
            annulus: *self,
            controls,
            inner,
            sign,
            transition: evaluation.transition,
            rotation: CarreauRotation {
                requested_wall_angular_velocities_rad_s: walls,
                accepted_wall_angular_velocities_rad_s: [walls[0], outer_speed],
                relative_angular_velocity_rad_s: relative,
                angular_velocity_residual_rad_s: residual,
                wall_tractions_pa: [-sign * inner.stress(), sign * outer_stress],
                wall_torques_n_m: [-torque, torque],
                wall_power_w: work,
                dissipated_power_w: heat,
                work_residual_w: work_residual,
                differential_resistance_n_m_s: resistance,
                maximum_shear_rate_per_s: inner.rate,
                quadrature_relative_change: evaluation.change,
                wall_stress_iterations: iterations,
            },
        })
    }
}

impl CarreauRotationProfile {
    /// Original source and geometry, shared with axial annular flow.
    pub const fn annulus(&self) -> &LaminarCarreauAnnulus {
        &self.annulus
    }
    /// Accepted physical torque, heating and boundary residuals.
    pub const fn rotation(&self) -> &CarreauRotation {
        &self.rotation
    }

    /// Sample a radius between the original material walls.
    pub fn sample(&self, radius_m: f64) -> Result<CarreauRotationSample, MaterialError> {
        if !radius_m.is_finite()
            || radius_m < self.annulus.inner_radius_m
            || radius_m > self.annulus.outer_radius_m
        {
            return Err(parameters(
                "Carreau rotation sample must lie between its cylinders",
            ));
        }
        self.sample_gap_fraction((radius_m - self.annulus.inner_radius_m) / self.annulus.gap_m)
    }

    /// Sample a physical gap fraction, preserving resolution in narrow gaps.
    pub fn sample_gap_fraction(
        &self,
        fraction: f64,
    ) -> Result<CarreauRotationSample, MaterialError> {
        self.sample_gap_fraction_with_check(fraction, || Ok(()))
    }

    /// Cancellable sampling; the accepted profile is immutable.
    pub fn sample_gap_fraction_with_check(
        &self,
        fraction: f64,
        mut check: impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<CarreauRotationSample, MaterialError> {
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err(parameters(
                "Carreau rotation gap fraction must lie in [0,1]",
            ));
        }
        check()?;
        let a = &self.annulus;
        let x = if fraction == 1.0 {
            1.0
        } else {
            (log_ratio_from_gap(a.inner_radius_m, fraction * a.gap_m) / a.log_ratio).clamp(0.0, 1.0)
        };
        let local = a.rotation_local(x, self.inner, self.controls, &mut check)?;
        let stress = a.rotation_stress(x, self.inner)?;
        let (angular, change) = if fraction == 0.0 {
            (self.rotation.accepted_wall_angular_velocities_rad_s[0], 0.0)
        } else if fraction == 1.0 {
            (self.rotation.accepted_wall_angular_velocities_rad_s[1], 0.0)
        } else if x <= 0.5 {
            let (relative, change) = a.rotation_integral(
                self.inner,
                [0.0, x],
                self.transition,
                self.controls,
                &mut check,
            )?;
            (
                self.rotation.accepted_wall_angular_velocities_rad_s[0] + self.sign * relative,
                change,
            )
        } else {
            let (relative, change) = a.rotation_integral(
                self.inner,
                [x, 1.0],
                self.transition,
                self.controls,
                &mut check,
            )?;
            (
                self.rotation.accepted_wall_angular_velocities_rad_s[1] - self.sign * relative,
                change,
            )
        };
        let radius = a.gap_m.mul_add(fraction, a.inner_radius_m);
        let velocity = radius * angular;
        let heat = if local.rate == 0.0 {
            0.0
        } else {
            product_ratio(&[stress, local.rate], &[])?
        };
        if !angular.is_finite() || !velocity.is_finite() {
            return Err(state_error(
                "Carreau rotational velocity is unrepresentable",
            ));
        }
        check()?;
        Ok(CarreauRotationSample {
            angular_velocity_rad_s: angular,
            azimuthal_velocity_m_s: velocity,
            shear_rate_per_s: self.sign * local.rate,
            radial_shear_stress_pa: self.sign * stress,
            dynamic_viscosity_pa_s: local.viscosity,
            dissipation_w_m3: heat,
            quadrature_relative_change: change,
        })
    }
}
