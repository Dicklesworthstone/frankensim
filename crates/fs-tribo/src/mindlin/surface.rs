//! Radial fields of the original decoupled, proportional Hertz--Mindlin model.
//!
//! For the retained spring strain `z(x)`, the line load is `q(x)=k0*z(x)/(2a)`.
//! Its inverse Abel transform gives `tau(r)=-integral_r^a q'(x)/sqrt(x*x-r*r)/pi`
//! and `u(r)=2*integral_0^min(r,a) z(x)/sqrt(r*r-x*x)/pi`. Each stored strain
//! is affine in `1-x*x/a*a`, so both integrals are analytic. Pressure and shear
//! vanish outside the patch; elastic displacement there generally does not.
//!
//! These are the scalar radial fields of the parent's decoupled approximation,
//! not a general vector elastic half-space solution. `u` is the relative elastic
//! surface displacement of both bodies, not plastic slip, accumulated sliding
//! distance, indentation or a displacement of either body alone. A local slip
//! history additionally needs the macroscopic motion and contact entry/exit.
//!
//! Reference: Popov, Hess and Willert (2019), section 4.4, Eq. 4.25 and the
//! Cattaneo superposition in Eqs. 4.34 and 4.58:
//! <https://doi.org/10.1007/978-3-662-58709-6_4>.

use super::{MemoryInterval, MindlinDeresiewiczPoint, MindlinState};
use crate::{TriboError, finite, nonnegative_finite};
use core::f64::consts::PI;
use fs_math::det;

/// Irreversible local sliding and friction work during one fixed-load segment.
pub mod sliding;

/// Local accepted-state fields in the original fixed first tangent direction.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MindlinSurfaceResponse {
    /// Compressive Hertz pressure [Pa]; zero at and outside the contact edge.
    pub normal_pressure_pa: f64,
    /// Signed tangential traction entering the interface [Pa]. Its area integral
    /// is the retained endpoint force, not the last step's time-average force.
    pub tangential_traction_pa: f64,
    /// Signed relative elastic surface displacement along the first tangent [m].
    /// The exterior field decays to zero; its reference is the remote surfaces.
    pub elastic_displacement_m: f64,
}

impl MindlinDeresiewiczPoint {
    /// Reconstruct a radial point from the original accepted history.
    ///
    /// Any finite nonnegative `radius_m` is admitted, including the center,
    /// contact edge and exterior. The history must belong to this exact law.
    /// Evaluation takes O(retained intervals), allocates nothing and never
    /// changes the checkpoint or introduces a radial discretization.
    pub fn surface_response(
        &self,
        state: &MindlinState,
        radius_m: f64,
    ) -> Result<MindlinSurfaceResponse, TriboError> {
        if state.owner != *self {
            return Err(TriboError::InvalidInput {
                field: "mindlin_checkpoint_owner",
            });
        }
        nonnegative_finite(radius_m, "mindlin_surface_radius_m")?;
        let a = self.parameters.contact_radius_m;
        let mut result = MindlinSurfaceResponse::default();
        if radius_m < a {
            let radius = radius_m / a;
            let cutoff = (1.0 - radius) * (1.0 + radius);
            let pressure_scale = (self.parameters.normal_force_n / a) * (1.5 / PI) / a;
            result.normal_pressure_pa = pressure_scale * det::sqrt(cutoff);
            let mut lower = 0.0;
            let mut traction = 0.0;
            for interval in &state.memory {
                if lower >= cutoff {
                    break;
                }
                let upper = interval.upper.min(cutoff);
                let outer = det::sqrt(cutoff - lower);
                let inner = det::sqrt(cutoff - upper);
                // Rationalize the square-root difference for a thin slip ring.
                traction += interval.slope * ((upper - lower) / (outer + inner));
                lower = interval.upper;
            }
            result.tangential_traction_pa = ((self.capacity_n / a) * (1.5 / PI) / a) * traction;
        }
        result.elastic_displacement_m =
            self.critical_slip_m * displacement(&state.memory, radius_m, a);
        finite(result.normal_pressure_pa, "mindlin_surface_pressure_pa")?;
        finite(result.tangential_traction_pa, "mindlin_surface_traction_pa")?;
        finite(
            result.elastic_displacement_m,
            "mindlin_surface_displacement_m",
        )?;
        Ok(result)
    }
}

// Integrate in theta=asin(x/radius). All weights are positive, so a small
// displacement is not obtained by subtracting large exterior-field primitives.
fn displacement(memory: &[MemoryInterval], radius_m: f64, a: f64) -> f64 {
    let Some(last) = memory.last() else {
        return 0.0;
    };
    let inside = radius_m < a;
    let normalized = radius_m / a;
    let radius_squared = normalized * normalized;
    let last_lower = if memory.len() > 1 {
        memory[memory.len() - 2].upper
    } else {
        0.0
    };
    // The entire integration disk lies in one quadratic. This also evaluates
    // the exact center limit without dividing by its radius or losing 1-r^2.
    if inside && radius_squared <= 1.0 - last_lower {
        return last.slope + last.offset - 0.5 * last.slope * radius_squared;
    }
    let scale = a / radius_m;
    if scale == 0.0 {
        return 0.0;
    }
    let cutoff = if inside {
        (1.0 - normalized) * (1.0 + normalized)
    } else {
        0.0
    };
    let mut lower: f64 = 0.0;
    let mut integral = 0.0;
    for interval in memory {
        let low_r = lower.max(cutoff);
        let high_r = interval.upper;
        lower = high_r;
        if high_r <= low_r {
            continue;
        }
        let low_x = det::sqrt(1.0 - high_r);
        let high_x = det::sqrt(1.0 - low_r);
        let low_sine = (scale * low_x).min(1.0);
        let high_sine = (scale * high_x).min(1.0);
        let (low_cosine, high_cosine) = if inside {
            (
                det::sqrt((high_r - cutoff) / radius_squared),
                det::sqrt((low_r - cutoff) / radius_squared),
            )
        } else {
            // Preserve a thin edge interval's original radial coordinate;
            // squaring its rounded sine would lose tiny microslip fields.
            let gap = (1.0 - scale) * (1.0 + scale);
            (
                det::sqrt(gap + scale * scale * high_r),
                det::sqrt(gap + scale * scale * low_r),
            )
        };
        let width = high_r - low_r;
        let denominator = high_x * low_cosine + low_x * high_cosine;
        let radial_sine = width / denominator;
        let sine = (scale * radial_sine).min(1.0);
        let cosine = low_cosine * high_cosine + low_sine * high_sine;
        let angle = det::atan2(sine, cosine);
        let cosine_sum = low_cosine * high_cosine - low_sine * high_sine;
        // Exact mean x^2 in this angular interval, as a fraction of its range.
        // g(theta)=(theta*cos(theta)-sin(theta))/(2*theta*sin(theta)^2).
        // Its removable zero limit prevents far-field and thin-band cancellation.
        let g = if angle < 0.01 {
            let t = angle * angle;
            -1.0 / 6.0 - t * (7.0 / 180.0 + t * (31.0 / 5040.0 + t * 127.0 / 151200.0))
        } else {
            (angle * cosine - sine) / (2.0 * angle * sine * sine)
        };
        let fraction = 0.5 + cosine_sum * (radial_sine / denominator) * g;
        let strain_low = interval.slope * high_r + interval.offset;
        let strain_high = interval.slope * low_r + interval.offset;
        integral += angle * ((1.0 - fraction) * strain_low + fraction * strain_high);
    }
    (2.0 / PI) * integral
}
