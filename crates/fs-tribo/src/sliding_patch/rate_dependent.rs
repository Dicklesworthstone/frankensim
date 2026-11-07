//! Original dry-friction curves resolved over a sliding and spinning Hertz disk.
//!
//! A spinning patch has a continuum of local slip speeds. Evaluating the law
//! only at its centre (or at an RMS speed) cannot recover its force, moment, or
//! heat. This owner integrates the shared [`FrictionLaw`] at each local speed,
//! retaining the original Hertz pressure and rigid patch velocity field.
//!
//! Coulomb and unclipped affine speed laws use exact pressure moments. Other
//! curves integrate their difference from the exact zero-speed Coulomb
//! limit. That difference removes the direction singularity at a stationary
//! point, without introducing a speed regularizer. Deterministic composite
//! four-point Gauss quadrature has caller-owned radial and angular resolution;
//! it is a convergent spatial discretization, not an error certificate.
//! Clipped affine motion integrates its actual positive-friction region,
//! including an active disk crossing the Hertz boundary or lying outside it.
//! A shrinking disk or crescent therefore cannot vanish between full-patch
//! nodes or through subtraction of Coulomb responses.
//!
//! This is instantaneous prescribed gross sliding. It does not choose a static
//! reaction, solve velocity-weakening branch selection, evolve the contact,
//! partition heat, or infer an interface curve from a material name.

use core::f64::consts::PI;

use fs_math::det;

use super::{HertzSlidingPatch, HertzSlidingPatchParameters, HertzSlidingWrench, finite};
use crate::{FrictionLaw, TriboError};

mod clipped;

// Four-point Gauss--Legendre rule on [-1, 1].
const GAUSS: [(f64, f64); 4] = [
    (-0.861_136_311_594_052_6, 0.347_854_845_137_453_85),
    (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.861_136_311_594_052_6, 0.347_854_845_137_453_85),
];
const MAXIMUM_PANELS: usize = 256;

/// Explicit spatial resolution for the local, non-affine friction response.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HertzFrictionQuadrature {
    /// Panels in the radial transform; four nodes per panel. Full-patch rules
    /// use `r = a sin(theta)`; clipped rules follow the actual active span.
    pub radial_panels: usize,
    /// Panels per reflection-symmetric angular interval; four nodes per panel.
    /// An intersecting active disk splits at its circle intersection. Pure spin
    /// uses radial quadrature only.
    pub angular_panels: usize,
}

impl HertzFrictionQuadrature {
    fn validate(self) -> Result<(), TriboError> {
        // At most 1,048,576 original-law evaluations; no data-dependent growth.
        if self.radial_panels == 0
            || self.angular_panels == 0
            || self.radial_panels > MAXIMUM_PANELS
            || self.angular_panels > MAXIMUM_PANELS
        {
            return Err(TriboError::InvalidInput {
                field: "hertz_friction_quadrature_panels",
            });
        }
        Ok(())
    }
}

/// An explicit shared friction law and prescribed circular Hertz geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct HertzFrictionPatch {
    law: FrictionLaw,
    zero_speed_patch: HertzSlidingPatch,
}

impl HertzFrictionPatch {
    /// Validate the shared law and original normal-contact inputs once.
    ///
    /// These explicit coefficients retain the caller's source authority. A
    /// material adapter remains responsible for its source-card applicability.
    pub fn new(
        normal_force_n: f64,
        contact_radius_m: f64,
        law: FrictionLaw,
    ) -> Result<Self, TriboError> {
        let zero_speed_patch = HertzSlidingPatch::new(HertzSlidingPatchParameters {
            normal_force_n,
            contact_radius_m,
            friction_coefficient: law.zero_speed_kinetic_coefficient()?,
        })?;
        Ok(Self {
            law,
            zero_speed_patch,
        })
    }

    /// Original shared curve; no patch-effective coefficient replaces it.
    #[must_use]
    pub const fn friction_law(&self) -> &FrictionLaw {
        &self.law
    }

    /// Force, spin moment, and local-slip heat from one pressure distribution.
    ///
    /// Positive local velocity gives opposing local traction. With a strongly
    /// weakening curve the resultant components need not each oppose their
    /// corresponding generalized velocity; their *total* work is dissipative.
    /// Complete rest selects zero, without claiming a sticking solution.
    ///
    /// Refine both quadrature panel counts for mixed slip/spin, or the radial
    /// count for pure spin, until the requested physical outputs converge. The
    /// exact pure-translation, Coulomb, and unclipped affine paths do not depend
    /// on these counts. Invalid controls are refused even on an exact path.
    pub fn evaluate(
        &self,
        slip_velocity_m_s: [f64; 2],
        spin_rad_s: f64,
        quadrature: HertzFrictionQuadrature,
    ) -> Result<HertzSlidingWrench, TriboError> {
        quadrature.validate()?;
        // Also validates the velocities and preserves the original exact
        // Coulomb degeneracies, including an unloaded patch.
        let baseline = self
            .zero_speed_patch
            .evaluate(slip_velocity_m_s, spin_rad_s)?;
        let parameters = self.zero_speed_patch.parameters();
        if baseline.at_rest || parameters.normal_force_n == 0.0 {
            return Ok(baseline);
        }
        let affine_slope = self.affine_slope();
        if affine_slope == Some(0.0) {
            return Ok(baseline);
        }
        let speed = finite(
            det::hypot(slip_velocity_m_s[0], slip_velocity_m_s[1]),
            "hertz_friction_slip_speed",
        )?;
        if spin_rad_s == 0.0 {
            return HertzSlidingPatch::new(HertzSlidingPatchParameters {
                friction_coefficient: self.law.kinetic_mu(speed)?,
                ..parameters
            })?
            .evaluate(slip_velocity_m_s, spin_rad_s);
        }
        let edge_speed = finite(
            parameters.contact_radius_m * spin_rad_s.abs(),
            "hertz_friction_edge_speed",
        )?;
        if let Some(slope) = affine_slope {
            let unclipped = if slope >= 0.0 {
                true
            } else {
                let maximum_speed = finite(speed + edge_speed, "hertz_friction_maximum_speed")?;
                finite(
                    parameters.friction_coefficient + slope * maximum_speed,
                    "hertz_affine_coefficient",
                )? >= 0.0
            };
            if unclipped {
                return self.affine(
                    baseline,
                    slip_velocity_m_s,
                    spin_rad_s,
                    speed,
                    edge_speed,
                    slope,
                );
            }
        }
        if let FrictionLaw::VelocityDependent {
            mu_zero,
            slope_per_speed,
            ..
        } = self.law
        {
            // The lower clamp can make the entire moving patch frictionless.
            // Retain that exact limit rather than subtracting two quadratures.
            if slope_per_speed < 0.0 && self.law.kinetic_mu((speed - edge_speed).max(0.0))? == 0.0 {
                return finish([0.0; 2], 0.0, 0.0, slip_velocity_m_s, spin_rad_s);
            }
            if speed == 0.0 {
                let core_radius = mu_zero / (-slope_per_speed * edge_speed);
                return clipped_affine_spin(
                    baseline,
                    core_radius,
                    spin_rad_s,
                    quadrature.radial_panels,
                );
            }
            let cutoff_speed = mu_zero / -slope_per_speed;
            if speed < edge_speed && cutoff_speed <= edge_speed - speed {
                return clipped_affine_contained(
                    parameters,
                    slip_velocity_m_s,
                    spin_rad_s,
                    speed,
                    edge_speed,
                    cutoff_speed,
                    quadrature,
                );
            }
            return clipped::intersecting_disk(
                parameters,
                slip_velocity_m_s,
                spin_rad_s,
                speed,
                edge_speed,
                cutoff_speed,
                quadrature,
            );
        }
        finite(speed + edge_speed, "hertz_friction_maximum_speed")?;
        let mut correction = Compensated::default();
        let mut sampled_traction = false;
        let radial_half = (PI / 4.0) / quadrature.radial_panels as f64;
        let angular_half = (PI / 2.0) / quadrature.angular_panels as f64;
        // Cache the radial geometry once, not once per angular node. Fixed
        // stack storage keeps evaluation allocation-free and bounded.
        let mut radial_nodes = [(0.0, 0.0); 4 * MAXIMUM_PANELS];
        for radial in 0..quadrature.radial_panels {
            for (index, (node, weight)) in GAUSS.into_iter().enumerate() {
                let theta = (2.0 * radial as f64 + 1.0 + node) * radial_half;
                let (r, c) = det::sin_cos(theta);
                // dN/N = 3 sin(theta) cos(theta)^2 dtheta dphi/(2 pi).
                radial_nodes[4 * radial + index] = (r, 3.0 * r * c * c * radial_half * weight);
            }
        }
        let radial_nodes = &radial_nodes[..4 * quadrature.radial_panels];
        if speed == 0.0 {
            for &(r, weight) in radial_nodes {
                let local_speed = edge_speed * r;
                let coefficient = self.law.kinetic_mu(local_speed)?;
                sampled_traction |= coefficient > 0.0;
                let delta = coefficient - parameters.friction_coefficient;
                correction.add([0.0, weight * delta * r, weight * delta * local_speed]);
            }
        } else {
            for angular in 0..quadrature.angular_panels {
                for (node, weight) in GAUSS {
                    let phi = (2.0 * angular as f64 + 1.0 + node) * angular_half;
                    let (sine, cosine) = det::sin_cos(phi);
                    let angular_weight = (angular_half / PI) * weight;
                    for &(r, pressure_weight) in radial_nodes {
                        let local = [speed - edge_speed * r * cosine, edge_speed * r * sine];
                        let local_speed =
                            finite(det::hypot(local[0], local[1]), "hertz_friction_local_speed")?;
                        // A zero-speed point carries no area and zero
                        // correction: mu(v)-mu(0) tends to zero there.
                        if local_speed > 0.0 {
                            let coefficient = self.law.kinetic_mu(local_speed)?;
                            sampled_traction |= coefficient > 0.0;
                            let delta = coefficient - parameters.friction_coefficient;
                            let load = pressure_weight * angular_weight * delta;
                            let direction = local.map(|v| v / local_speed);
                            correction.add([
                                load * direction[0],
                                load * r * (sine * direction[1] - cosine * direction[0]),
                                load * local_speed,
                            ]);
                        }
                    }
                }
            }
        }
        if !sampled_traction {
            // Exact all-zero laws and the fully clipped region exited above.
            // A remaining active region was missed at this spatial resolution.
            return Err(TriboError::InvalidInput {
                field: "hertz_friction_unresolved_active_region",
            });
        }
        let correction = correction.total;
        let normal = parameters.normal_force_n;
        let force_correction = normal * correction[0];
        let force_n = if speed == 0.0 {
            [0.0; 2]
        } else {
            [0, 1].map(|i| baseline.force_n[i] - force_correction * (slip_velocity_m_s[i] / speed))
        };
        finish(
            force_n,
            baseline.spin_moment_nm
                - (normal * correction[1] * parameters.contact_radius_m) * spin_rad_s.signum(),
            baseline.dissipated_power_w + normal * correction[2],
            slip_velocity_m_s,
            spin_rad_s,
        )
    }

    fn affine_slope(&self) -> Option<f64> {
        match self.law {
            FrictionLaw::Coulomb { .. } => Some(0.0),
            FrictionLaw::VelocityDependent {
                slope_per_speed, ..
            } => Some(slope_per_speed),
            FrictionLaw::Stribeck {
                static_mu,
                kinetic_mu,
                viscous_per_speed,
                ..
            } if static_mu == kinetic_mu => Some(viscous_per_speed),
            FrictionLaw::Hyperbolic {
                static_mu,
                kinetic_mu,
                ..
            }
            | FrictionLaw::ShiftedPower {
                static_mu,
                kinetic_mu,
                ..
            } if static_mu == kinetic_mu => Some(0.0),
            _ => None,
        }
    }

    fn affine(
        &self,
        baseline: HertzSlidingWrench,
        slip: [f64; 2],
        spin: f64,
        speed: f64,
        edge_speed: f64,
        slope: f64,
    ) -> Result<HertzSlidingWrench, TriboError> {
        if slope == 0.0 {
            return Ok(baseline);
        }
        let parameters = self.zero_speed_patch.parameters();
        let translation_drag = (slope * speed) * parameters.normal_force_n;
        // Integral of r^2 p dA = (2/5) N a^2 for Hertz pressure.
        let scaled_spin_drag = ((slope * edge_speed) * parameters.normal_force_n) * (2.0 / 5.0);
        let force_n = if speed == 0.0 {
            [0.0; 2]
        } else {
            [0, 1].map(|i| baseline.force_n[i] - translation_drag * (slip[i] / speed))
        };
        finish(
            force_n,
            baseline.spin_moment_nm
                - scaled_spin_drag * parameters.contact_radius_m * spin.signum(),
            baseline.dissipated_power_w + translation_drag * speed + scaled_spin_drag * edge_speed,
            slip,
            spin,
        )
    }
}

fn clipped_affine_contained(
    parameters: HertzSlidingPatchParameters,
    slip: [f64; 2],
    spin: f64,
    speed: f64,
    edge_speed: f64,
    cutoff_speed: f64,
    quadrature: HertzFrictionQuadrature,
) -> Result<HertzSlidingWrench, TriboError> {
    // In unit-patch coordinates the zero-slip point is distance d from the
    // centre. The positive-friction disk has radius k = v_cut/(a*|omega|).
    // Use local-slip polar coordinates q=k*t, mu=mu_zero*(1-t), rather than
    // a grid on the full Hertz disk. The caller proved d+k <= 1.
    let d = speed / edge_speed;
    let k = cutoff_speed / edge_speed;
    let radial_half = 0.5 / quadrature.radial_panels as f64;
    let angular_half = (PI / 4.0) / quadrature.angular_panels as f64;
    let mut radial_nodes = [(0.0, 0.0); 4 * MAXIMUM_PANELS];
    for panel in 0..quadrature.radial_panels {
        for (i, (node, weight)) in GAUSS.into_iter().enumerate() {
            let t = (2.0 * panel as f64 + 1.0 + node) * radial_half;
            radial_nodes[4 * panel + i] = (t, radial_half * weight * t * (1.0 - t));
        }
    }
    let mut integral = Compensated::default();
    for panel in 0..quadrature.angular_panels {
        for (node, weight) in GAUSS {
            let phi = (2.0 * panel as f64 + 1.0 + node) * angular_half;
            let (sine, cosine) = det::sin_cos(phi);
            for &(t, radial_weight) in &radial_nodes[..4 * quadrature.radial_panels] {
                let q = k * t;
                let y = q * sine;
                let near_x = d - q * cosine;
                let far_x = d + q * cosine;
                // Both points are inside the admitted disk. Only endpoint
                // roundoff can make these pressure-square arguments negative.
                let near = det::sqrt(((1.0 - near_x) * (1.0 + near_x) - y * y).max(0.0));
                let far = det::sqrt(((1.0 - far_x) * (1.0 + far_x) - y * y).max(0.0));
                let sum = near + far;
                // Pair opposite local-slip directions analytically. Subtracting
                // nearly equal pressures loses the small physical net force as
                // the active disk shrinks; the rationalized difference does not.
                let difference = (4.0 * d * q * cosine) / sum;
                let weight = radial_weight * angular_half * weight;
                integral.add([
                    weight * cosine * difference,
                    weight * (q * sum - d * cosine * difference),
                    weight * t * sum,
                ]);
            }
        }
    }
    // The reflected half-disk supplies 3/pi times the Hertz load. Multiply the
    // dimensional scale before the small radius factors to retain tiny cores.
    let scale =
        (((3.0 / PI) * parameters.normal_force_n * parameters.friction_coefficient) * k) * k;
    let force = -scale * integral.total[0];
    finish(
        slip.map(|v| force * (v / speed)),
        -(scale * integral.total[1] * parameters.contact_radius_m) * spin.signum(),
        scale * cutoff_speed * integral.total[2],
        slip,
        spin,
    )
}

fn clipped_affine_spin(
    baseline: HertzSlidingWrench,
    core_radius: f64,
    spin: f64,
    radial_panels: usize,
) -> Result<HertzSlidingWrench, TriboError> {
    // r/a = core_radius*t and mu = mu_zero*(1-t), 0 <= t <= 1.
    // M = -3 N a mu_zero core_radius^3 integral t^2(1-t)
    //          sqrt(1-core_radius^2*t^2) dt.
    // Integrate only positive local friction, without subtracting a large
    // Coulomb moment. Scaling the known baseline retains its original units.
    let half = 0.5 / radial_panels as f64;
    let mut integral = Compensated::default();
    for panel in 0..radial_panels {
        for (node, weight) in GAUSS {
            let t = (2.0 * panel as f64 + 1.0 + node) * half;
            let r = core_radius * t;
            let pressure = det::sqrt((1.0 - r) * (1.0 + r));
            integral.add([half * weight * t * t * (1.0 - t) * pressure, 0.0, 0.0]);
        }
    }
    let factor = (16.0 / PI) * integral.total[0];
    // Multiply the dimensional result before each small radius factor, so a
    // representable result does not vanish when core_radius^3 underflows.
    let scale = |value| ((value * core_radius) * core_radius) * core_radius * factor;
    finish(
        [0.0; 2],
        scale(baseline.spin_moment_nm),
        scale(baseline.dissipated_power_w),
        [0.0; 2],
        spin,
    )
}

fn finish(
    force_n: [f64; 2],
    spin_moment_nm: f64,
    dissipated_power_w: f64,
    slip: [f64; 2],
    spin: f64,
) -> Result<HertzSlidingWrench, TriboError> {
    for force in force_n {
        finite(force, "hertz_rate_dependent_force")?;
    }
    finite(spin_moment_nm, "hertz_rate_dependent_moment")?;
    finite(dissipated_power_w, "hertz_rate_dependent_power")?;
    if dissipated_power_w < 0.0 {
        return Err(TriboError::InvalidInput {
            field: "hertz_rate_dependent_passivity",
        });
    }
    let power_residual_w = finite(
        force_n[0] * slip[0] + force_n[1] * slip[1] + spin_moment_nm * spin + dissipated_power_w,
        "hertz_rate_dependent_power_residual",
    )?;
    Ok(HertzSlidingWrench {
        force_n,
        spin_moment_nm,
        dissipated_power_w,
        power_residual_w,
        at_rest: false,
    })
}

#[derive(Default)]
struct Compensated {
    total: [f64; 3],
    error: [f64; 3],
}

impl Compensated {
    fn add(&mut self, value: [f64; 3]) {
        for (i, value) in value.into_iter().enumerate() {
            let next = value - self.error[i];
            let sum = self.total[i] + next;
            self.error[i] = (sum - self.total[i]) - next;
            self.total[i] = sum;
        }
    }
}
