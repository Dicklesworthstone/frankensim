//! Pressure-resolved Coulomb sliding and spin on a circular Hertz patch.
//!
//! Local relative velocity is `u + omega * [-y, x]` and local traction is
//! `-mu * p(r) * velocity / |velocity|`, with
//! `p(r) = 3 P / (2 pi a^2) * sqrt(1-r^2/a^2)`. Its exact surface integrals
//! give a coupled two-axis force and normal-axis moment. They are not separate
//! saturated force/torque limits and do not use an RMS patch lever radius.
//!
//! The closed forms are Antali (2024), *Nonlinear Dynamics* 112, 6785--6811,
//! DOI: 10.1007/s11071-024-09462-6, equations (62)--(63). The paper normalizes
//! moment by the pressure first-moment radius, which here is `3 pi a / 16`.
//! The small-spin branch evaluates the same analytic function with a convergent
//! series to avoid cancellation. Work is bounded and allocation-free.
//!
//! This is an instantaneous gross-sliding constitutive response under supplied
//! rigid patch kinematics. An isolated zero-slip point has zero area. At complete
//! rest the zero wrench is returned as one admissible selection; this function
//! does not solve set-valued static friction or decide slip onset. There is no
//! retained elastic shear, rolling-deformation loss, normal-contact evolution,
//! adhesion, material admission, finite-step motion or thermal partition. In
//! particular this is not a nonproportional Mindlin partial-slip approximation.

use core::f64::consts::PI;

use fs_math::det;

use crate::{TriboError, nonnegative_finite, positive_finite};

/// Coupled static/sliding contact with full translational and torsional recoil.
pub mod contact;

fn finite(value: f64, field: &'static str) -> Result<f64, TriboError> {
    crate::finite(value, field)?;
    Ok(value)
}

// The two series are paired by the dissipation identity, not independently fit.
// Compile-time coefficients keep repeated constitutive evaluations to Horner
// multiplication/addition without per-term division or any stored grid.
const SMALL_SPIN_COEFFICIENTS: [[f64; 2]; 32] = small_spin_coefficients();

const fn small_spin_coefficients() -> [[f64; 2]; 32] {
    let mut coefficients = [[0.0; 2]; 32];
    let mut asin_coefficient = 1.0 / 6.0;
    let mut sqrt_coefficient = -0.5;
    let mut index = 0;
    while index < 32 {
        let k = (index + 1) as f64;
        let next_asin = asin_coefficient * ((2.0 * k + 1.0) * (2.0 * k + 1.0))
            / ((2.0 * k + 2.0) * (2.0 * k + 3.0));
        let next_sqrt = sqrt_coefficient * ((k - 0.5) / (k + 1.0));
        let coefficient = (3.0 / 16.0)
            * (4.0 * asin_coefficient - next_asin + next_sqrt + 2.0 * sqrt_coefficient);
        coefficients[index] = [coefficient, -(2.0 * k) / (2.0 * k - 1.0) * coefficient];
        asin_coefficient = next_asin;
        sqrt_coefficient = next_sqrt;
        index += 1;
    }
    coefficients
}

/// Physical data supplied by the original normal and interface owners.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzSlidingPatchParameters {
    /// Prescribed compressive force [N]; zero gives zero traction.
    pub normal_force_n: f64,
    /// Circular contact radius [m]. No radius is inferred from a bulk name.
    pub contact_radius_m: f64,
    /// Explicit local kinetic Coulomb coefficient, uniform over the patch.
    pub friction_coefficient: f64,
}

/// Immutable exact Hertz-pressure sliding/spinning response.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzSlidingPatch {
    parameters: HertzSlidingPatchParameters,
    capacity_n: f64,
    spin_capacity_nm: f64,
}

/// Resultant acting on the side whose relative velocity was supplied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzSlidingWrench {
    /// Tangential force opposing the centre slip velocity [N].
    pub force_n: [f64; 2],
    /// Moment about the patch centre, opposing relative normal-axis spin [N m].
    pub spin_moment_nm: f64,
    /// Integral of `mu p |u + omega cross r|` [W]. No partition is assigned.
    pub dissipated_power_w: f64,
    /// `force dot u + moment * omega + dissipated_power` [W].
    pub power_residual_w: f64,
    /// Both supplied velocities are exactly zero. This is not a stiction solve.
    pub at_rest: bool,
}

impl HertzSlidingPatch {
    /// Admit explicit finite patch inputs without changing their source authority.
    pub fn new(parameters: HertzSlidingPatchParameters) -> Result<Self, TriboError> {
        nonnegative_finite(parameters.normal_force_n, "normal_force_n")?;
        positive_finite(parameters.contact_radius_m, "contact_radius_m")?;
        nonnegative_finite(parameters.friction_coefficient, "friction_coefficient")?;
        let capacity_n = finite(
            parameters.normal_force_n * parameters.friction_coefficient,
            "hertz_sliding_capacity_n",
        )?;
        let spin_capacity_nm = finite(
            capacity_n * ((3.0 * PI / 16.0) * parameters.contact_radius_m),
            "hertz_spinning_capacity_nm",
        )?;
        Ok(Self {
            parameters,
            capacity_n,
            spin_capacity_nm,
        })
    }

    /// Original admitted load, radius and local coefficient, unchanged.
    #[must_use]
    pub const fn parameters(&self) -> HertzSlidingPatchParameters {
        self.parameters
    }

    /// Pure translational sliding force `mu P` [N].
    #[must_use]
    pub const fn sliding_capacity_n(&self) -> f64 {
        self.capacity_n
    }

    /// Pure spinning moment `3 pi mu P a / 16` [N m].
    #[must_use]
    pub const fn spinning_capacity_nm(&self) -> f64 {
        self.spin_capacity_nm
    }

    /// Exact coupled wrench for arbitrary two-axis slip [m/s] and spin [rad/s].
    ///
    /// Both components use one pressure distribution and the same local Coulomb
    /// law. Superposing pure sliding and pure spinning capacities is incorrect.
    /// No speed regularizer is used; positive rescaling of all velocities leaves
    /// the wrench unchanged and rescales physical power by the same factor.
    pub fn evaluate(
        &self,
        slip_velocity_m_s: [f64; 2],
        spin_rad_s: f64,
    ) -> Result<HertzSlidingWrench, TriboError> {
        for value in slip_velocity_m_s {
            finite(value, "slip_velocity_m_s")?;
        }
        finite(spin_rad_s, "spin_rad_s")?;
        let speed = finite(
            det::hypot(slip_velocity_m_s[0], slip_velocity_m_s[1]),
            "slip_speed_m_s",
        )?;
        let at_rest = speed == 0.0 && spin_rad_s == 0.0;
        if at_rest || self.capacity_n == 0.0 {
            return Ok(HertzSlidingWrench {
                force_n: [0.0; 2],
                spin_moment_nm: 0.0,
                dissipated_power_w: 0.0,
                power_residual_w: 0.0,
                at_rest,
            });
        }
        let (force_ratio, spin_ratio) = if spin_rad_s == 0.0 {
            (1.0, 0.0)
        } else if speed == 0.0 {
            (0.0, 1.0)
        } else {
            let edge_speed = finite(
                self.parameters.contact_radius_m * spin_rad_s.abs(),
                "patch_edge_spin_speed_m_s",
            )?;
            shape_factors(speed, edge_speed)
        };
        let force_magnitude = finite(self.capacity_n * force_ratio, "sliding_force_n")?;
        let moment_magnitude = finite(self.spin_capacity_nm * spin_ratio, "spinning_moment_nm")?;
        let force_n = if speed == 0.0 {
            [0.0; 2]
        } else {
            slip_velocity_m_s.map(|v| -force_magnitude * (v / speed))
        };
        let spin_moment_nm = if spin_rad_s == 0.0 {
            0.0
        } else {
            -moment_magnitude.copysign(spin_rad_s)
        };
        // Euler's identity for the homogeneous convex dissipation integral:
        // D = |F| |u| + |M| |omega|. This is independently tested by integrating
        // local pressure times local slip speed over the full Hertz disk.
        let dissipated_power_w = finite(
            force_magnitude * speed + moment_magnitude * spin_rad_s.abs(),
            "sliding_patch_power_w",
        )?;
        let mechanical_power = finite(
            force_n[0] * slip_velocity_m_s[0]
                + force_n[1] * slip_velocity_m_s[1]
                + spin_moment_nm * spin_rad_s,
            "sliding_patch_mechanical_power_w",
        )?;
        let power_residual_w = finite(
            mechanical_power + dissipated_power_w,
            "sliding_patch_power_residual_w",
        )?;
        if force_ratio < 0.0 || spin_ratio < 0.0 || dissipated_power_w < 0.0 {
            return Err(TriboError::InvalidInput {
                field: "hertz_sliding_passivity",
            });
        }
        Ok(HertzSlidingWrench {
            force_n,
            spin_moment_nm,
            dissipated_power_w,
            power_residual_w,
            at_rest,
        })
    }
}

// Force normalized by mu P, moment normalized by 3 pi mu P a / 16.
fn shape_factors(speed: f64, edge_speed: f64) -> (f64, f64) {
    if speed <= edge_speed {
        let epsilon = speed / edge_speed;
        let square = epsilon * epsilon;
        return (
            (3.0 * PI / 32.0) * epsilon * (4.0 - square),
            (3.0 * square * square - 8.0 * square + 8.0) / 8.0,
        );
    }
    let z = edge_speed / speed;
    if z <= 0.5 {
        // If Q(z)=1+sum q_k z^(2k), homogeneity of the dissipation integral
        // gives M/(mu P a)=-sum 2k*q_k*z^(2k-1)/(2k-1). Derive q_k from
        // asin(z) and sqrt(1-z^2), so this is the SAME closed form, not a fit.
        // Thirty-two terms at z<=1/2 make the omitted tail < 2^-64 before
        // the decreasing coefficient factor; work is fixed and allocation-free.
        let square = z * z;
        let mut force = 0.0;
        let mut moment = 0.0;
        for coefficient in SMALL_SPIN_COEFFICIENTS.into_iter().rev() {
            force = force * square + coefficient[0];
            moment = moment * square + coefficient[1];
        }
        (1.0 + force * square, moment * z / (3.0 * PI / 16.0))
    } else {
        // Here 1 < epsilon < 2, so cancellation in the direct expression is
        // modest. The far-slip region uses the series above.
        let epsilon = 1.0 / z;
        let square = epsilon * epsilon;
        let angle = det::asin(z);
        let root = det::sqrt(square - 1.0);
        (
            (3.0 / (16.0 * epsilon)) * (square * (4.0 - square) * angle + (square + 2.0) * root),
            ((3.0 * square * square - 8.0 * square + 8.0) * angle + (6.0 - 3.0 * square) * root)
                / (4.0 * PI),
        )
    }
}
