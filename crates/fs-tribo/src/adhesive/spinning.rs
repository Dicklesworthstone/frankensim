//! Uniform junction shear over the actual circular solid-contact patch.
//!
//! At each point, `v = u + omega * [-y, x]` and the traction is
//! `-tau * v / |v|`. The resultant force and spin moment are integrals of the
//! SAME slip field; separately saturating force and torque is not this law.
//! A finite attached patch therefore resists spin even at tensile net load.
//!
//! The uniform-disk integrals are Antali (2024), Nonlinear Dynamics 112,
//! 6785--6811, DOI 10.1007/s11071-024-09462-6, equations (59)--(60).
//! His torque normalization is `2 tau A a / 3`. Complete elliptic integrals
//! use the parameter m, not its square root. Convergent paired series avoid
//! cancellation at small slip/spin ratios; AGM work is bounded near the rim.
//!
//! This is gross sliding on a supplied, frozen real area. It introduces no
//! pressure contribution, empirical friction coefficient, elastic shear,
//! peeling, normal evolution or heat allocation. At complete rest zero wrench
//! is one admissible selection, not a solution of the static reaction problem.
//! Evaluation is allocation-free and bounded by 32 series terms or 16 AGM steps.

use core::f64::consts::PI;

use fs_math::det;

use crate::{TriboError, finite, nonnegative_finite};

/// Static/sliding response under the actual three-coordinate body compliance.
pub mod contact;

/// Independently supplied solid-contact geometry and uniform junction strength.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveSlidingPatchParameters {
    /// True solid-contact radius, m. Zero describes a detached junction.
    pub contact_radius_m: f64,
    /// Uniform local shear capacity, Pa; independent of work of adhesion.
    pub shear_strength_pa: f64,
}

/// One circular uniform-shear patch. No net normal load is required.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveSlidingPatch {
    parameters: AdhesiveSlidingPatchParameters,
    capacity_n: f64,
    spin_capacity_nm: f64,
}

/// Actual integrated traction and spin moment on the supplied moving side.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveSlidingWrench {
    /// Tangential force on the side whose relative velocity was supplied, N.
    pub force_n: [f64; 2],
    /// Moment about the contact normal opposing relative spin, N m.
    pub spin_moment_nm: f64,
    /// Integral of `tau * |u + omega cross r|`, W; allocate at most once.
    pub dissipated_power_w: f64,
    /// `force dot u + moment * omega + dissipated_power`, W.
    pub power_residual_w: f64,
    /// Both physical rates are exactly zero; no static reaction is inferred.
    pub at_rest: bool,
}

impl AdhesiveSlidingPatch {
    /// Admit explicit real geometry and strength, including detached zero area.
    pub fn new(parameters: AdhesiveSlidingPatchParameters) -> Result<Self, TriboError> {
        nonnegative_finite(parameters.contact_radius_m, "adhesive patch radius")?;
        nonnegative_finite(
            parameters.shear_strength_pa,
            "adhesive patch shear strength",
        )?;
        let a = parameters.contact_radius_m;
        let capacity_n = parameters.shear_strength_pa * a * (PI * a);
        let spin_capacity_nm = capacity_n * ((2.0 / 3.0) * a);
        finite(capacity_n, "adhesive sliding capacity")?;
        finite(spin_capacity_nm, "adhesive spinning capacity")?;
        if a > 0.0
            && parameters.shear_strength_pa > 0.0
            && (capacity_n <= 0.0 || spin_capacity_nm <= 0.0)
        {
            return Err(invalid("unrepresentable adhesive patch capacity"));
        }
        Ok(Self {
            parameters,
            capacity_n,
            spin_capacity_nm,
        })
    }

    /// Original immutable circular geometry and local junction strength.
    #[must_use]
    pub const fn parameters(&self) -> AdhesiveSlidingPatchParameters {
        self.parameters
    }

    /// Pure translational capacity, `tau * pi * a^2`, N.
    #[must_use]
    pub const fn sliding_capacity_n(&self) -> f64 {
        self.capacity_n
    }

    /// Pure spinning capacity, `2 * tau * pi * a^3 / 3`, N m.
    #[must_use]
    pub const fn spinning_capacity_nm(&self) -> f64 {
        self.spin_capacity_nm
    }

    /// Coupled uniform-shear force and moment, with no velocity regularization.
    /// A positive rescaling of both rates leaves the wrench unchanged.
    pub fn evaluate(
        &self,
        slip_velocity_m_s: [f64; 2],
        spin_rad_s: f64,
    ) -> Result<AdhesiveSlidingWrench, TriboError> {
        for x in slip_velocity_m_s {
            finite(x, "adhesive slip velocity")?;
        }
        finite(spin_rad_s, "adhesive spin velocity")?;
        let speed = det::hypot(slip_velocity_m_s[0], slip_velocity_m_s[1]);
        finite(speed, "adhesive slip speed")?;
        let at_rest = speed == 0.0 && spin_rad_s == 0.0;
        if at_rest || self.capacity_n == 0.0 {
            return Ok(AdhesiveSlidingWrench {
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
            let edge_speed = self.parameters.contact_radius_m * spin_rad_s.abs();
            finite(edge_speed, "adhesive rim speed")?;
            shape_factors(speed, edge_speed)?
        };
        let force = self.capacity_n * force_ratio;
        let moment = self.spin_capacity_nm * spin_ratio;
        let force_n = if speed == 0.0 {
            [0.0; 2]
        } else {
            slip_velocity_m_s.map(|v| -force * (v / speed))
        };
        let spin_moment_nm = if spin_rad_s == 0.0 {
            0.0
        } else {
            -moment.copysign(spin_rad_s)
        };
        // Euler's identity for the one homogeneous dissipation integral.
        // Tests independently integrate the pointwise positive slip work.
        let dissipated_power_w = force * speed + moment * spin_rad_s.abs();
        let power_residual_w = force_n[0] * slip_velocity_m_s[0]
            + force_n[1] * slip_velocity_m_s[1]
            + spin_moment_nm * spin_rad_s
            + dissipated_power_w;
        for value in
            force_n
                .into_iter()
                .chain([spin_moment_nm, dissipated_power_w, power_residual_w])
        {
            finite(value, "adhesive sliding wrench")?;
        }
        if dissipated_power_w < 0.0 {
            return Err(invalid("adhesive sliding passivity"));
        }
        Ok(AdhesiveSlidingWrench {
            force_n,
            spin_moment_nm,
            dissipated_power_w,
            power_residual_w,
            at_rest,
        })
    }
}

// Q(z) = 1 + sum q_n z^(2n), q_n = -[(2n choose n)/4^n]^2 / ((2n-1)(n+1)).
// The same coefficients determine force AND torque by homogeneous dissipation.
const SERIES: [[f64; 3]; 32] = series_coefficients();

const fn series_coefficients() -> [[f64; 3]; 32] {
    let mut out = [[0.0; 3]; 32];
    let mut central = 1.0;
    let mut index = 0;
    while index < 32 {
        let n = (index + 1) as f64;
        let ratio = (2.0 * n - 1.0) / (2.0 * n);
        central *= ratio * ratio;
        let q = -central / ((2.0 * n - 1.0) * (n + 1.0));
        out[index] = [
            q,
            -3.0 * n * q / (2.0 * n - 1.0),
            -1.5 * (2.0 * n + 1.0) * q / (2.0 * n + 2.0),
        ];
        index += 1;
    }
    out
}

/// Force normalized by tau*A, moment by 2*tau*A*a/3.
pub(crate) fn shape_factors(speed: f64, edge_speed: f64) -> Result<(f64, f64), TriboError> {
    nonnegative_finite(speed, "adhesive shape slip speed")?;
    nonnegative_finite(edge_speed, "adhesive shape rim speed")?;
    if speed == 0.0 {
        return Ok((0.0, if edge_speed == 0.0 { 0.0 } else { 1.0 }));
    }
    if edge_speed == 0.0 {
        return Ok((1.0, 0.0));
    }
    let spin_dominant = speed <= edge_speed;
    let z = if spin_dominant {
        speed / edge_speed
    } else {
        edge_speed / speed
    };
    let m = z * z;
    let result = if z <= 0.5 {
        let mut q = 0.0;
        let mut small_spin = 0.0;
        let mut small_slip = 0.0;
        for c in SERIES.into_iter().rev() {
            q = q * m + c[0];
            small_spin = small_spin * m + c[1];
            small_slip = small_slip * m + c[2];
        }
        if spin_dominant {
            (z * (1.0 + m * q), 1.0 + m * (-0.75 + m * small_slip))
        } else {
            (1.0 + m * q, z * small_spin)
        }
    } else if z == 1.0 {
        // K(1) diverges, but (1-m)K(m) vanishes. Never form 0 * infinity.
        (8.0 / (3.0 * PI), 4.0 / (3.0 * PI))
    } else {
        let (k, e) = elliptic_ke(m)?;
        let one_minus_m = (1.0 - z) * (1.0 + z);
        let q_numerator = (1.0 + m) * e - one_minus_m * k;
        if spin_dominant {
            (
                (4.0 / (3.0 * PI)) * q_numerator / z,
                (2.0 / (3.0 * PI)) * ((4.0 - 2.0 * m) * e - one_minus_m * k),
            )
        } else {
            // Factor the vanishing K coefficient for accuracy near z=1.
            (
                (4.0 / (3.0 * PI)) * q_numerator / m,
                (2.0 / (3.0 * PI)) * ((4.0 * m - 2.0) * e + one_minus_m * (2.0 - 3.0 * m) * k)
                    / (z * m),
            )
        }
    };
    if !result.0.is_finite()
        || !result.1.is_finite()
        || result.0 < 0.0
        || result.1 < 0.0
        || result.0 > 1.0 + 32.0 * f64::EPSILON
        || result.1 > 1.0 + 32.0 * f64::EPSILON
    {
        return Err(invalid("unresolved adhesive slip-spin integral"));
    }
    Ok(result)
}

// Complete integrals K(m), E(m), with m in (1/4,1). The classical AGM
// converges quadratically; its complementary parameter remains representable.
fn elliptic_ke(m: f64) -> Result<(f64, f64), TriboError> {
    let mut a = 1.0;
    let mut b = det::sqrt(1.0 - m);
    let mut sum = 0.5 * m;
    let mut weight = 1.0;
    for _ in 0..16 {
        let c = 0.5 * (a - b);
        sum += weight * c * c;
        let next_a = 0.5 * (a + b);
        b = det::sqrt(a * b);
        a = next_a;
        if (a - b).abs() <= 2.0 * f64::EPSILON * a {
            let k = (0.5 * PI) / a;
            return Ok((k, k * (1.0 - sum)));
        }
        weight *= 2.0;
    }
    Err(invalid("adhesive elliptic iteration budget"))
}

fn invalid(field: &'static str) -> TriboError {
    TriboError::InvalidInput { field }
}
