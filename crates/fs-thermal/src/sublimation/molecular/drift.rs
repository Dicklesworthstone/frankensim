//! Incoming moments of a drifting Maxwellian at one local planar surface.
//!
//! The normal points FROM the solid INTO the gas. Velocities are relative to
//! the local surface, so positive normal drift moves away from it. Integrating
//! the incoming half-space changes mass, momentum AND crossing energy; adding
//! bulk kinetic energy to the stationary effusive energy is not sufficient.
//! This module neither advances a gas/body nor chooses a caloric reference.
//!
//! For the incoming mass integral see Ek et al. (2026), equation (1),
//! <https://doi.org/10.1111/jace.70896>. The higher moments below follow by
//! integrating the SAME Maxwellian with velocity and kinetic-energy weights.
//! Far-tail positive integrals obey the parabolic-cylinder recurrences in
//! <https://dlmf.nist.gov/12.5.E1> and <https://dlmf.nist.gov/12.8.E1>.
//! Local Maxwellian translational and equilibrated internal modes, independent
//! of velocity, are required. No Knudsen layer, multiple reflection, continuum
//! stress or velocity-dependent sticking law is inferred.

use crate::sublimation::{SublimationError, nonnegative, positive};
use fs_math::{
    dd::Dd,
    det::{erfc, exp, expm1, ln, sqrt},
};

const INV_SQRT_TAU: f64 = 0.398_942_280_401_432_7;
const LOG_SQRT_TAU: f64 = 0.918_938_533_204_672_7;
const NORMAL_TOLERANCE: f64 = 64.0 * f64::EPSILON;

/// Complete incoming species flux, BEFORE sticking or reflection selection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaxwellianIncidence {
    pub incoming_mass_flux_kg_m2_s: f64,
    /// Signed change from the same species' zero-normal-drift impingement.
    /// Retained independently even when the two gross rates round identically.
    pub incoming_mass_flux_change_from_zero_drift_kg_m2_s: f64,
    /// A genuinely empty species has no conditional molecular moments.
    pub moments: Option<IncidentMolecularMoments>,
}

/// Conditional moments of molecules actually crossing into the surface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IncidentMolecularMoments {
    /// Multiply by the ORIGINAL accepted incident mass to obtain its relative
    /// momentum; its component along the outward normal is negative.
    pub mean_relative_velocity_m_s: [f64; 3],
    /// Add to the ORIGINAL species enthalpy h_i(T_g), including its original
    /// formation/reference datum, to obtain relative-frame crossing energy.
    /// At zero normal drift this is |v_t|²/2 - R_i*T_g/2. It is not extra heat
    /// to add after transferring that complete molecular crossing energy.
    pub crossing_energy_minus_enthalpy_j_kg: f64,
}

fn invalid(message: &'static str) -> SublimationError {
    SublimationError::Invalid(message)
}

fn represented(value: f64) -> Result<f64, SublimationError> {
    if positive(value) {
        Ok(value)
    } else {
        Err(invalid(
            "nonzero drifting molecular moment is unrepresentable",
        ))
    }
}

// A power-of-two scale makes divisions exact except unavoidable subnormal
// rounding. Scaling BEFORE dot products avoids overflow followed by cancellation.
fn binary_scale(value: f64) -> f64 {
    f64::from_bits(value.to_bits() & 0x7ff0_0000_0000_0000).max(f64::MIN_POSITIVE)
}

fn dot(a: [f64; 3], b: [f64; 3]) -> Result<f64, SublimationError> {
    let largest = a.into_iter().map(f64::abs).fold(0.0, f64::max);
    if largest == 0.0 {
        return Ok(0.0);
    }
    let scale = binary_scale(largest);
    let sum = a.into_iter().zip(b).fold(Dd::ZERO, |sum, (x, y)| {
        sum + Dd::from_f64(x / scale) * Dd::from_f64(y)
    });
    let value = (sum * Dd::from_f64(scale)).to_f64();
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "relative molecular velocity projection is unrepresentable",
        ))
    }
}

fn from_log(value: Dd) -> Result<f64, SublimationError> {
    // Splitting the exponent also handles representable subnormals without an
    // underflowing exponential that a later physical scale could have rescued.
    let half = exp(0.5 * value.hi);
    represented((half * exp(value.lo)) * half)
}

// A normal a*b with a normal positive result keeps (a*b)*c: two roundings,
// relative error at most 2u + u^2 (u = 2^-53), bit-equal to the earlier code.
// Every other case is correctly rounded by exact_product, so a nonzero exact
// product that rounds to 0 refuses through represented(). A factor that is not
// finite and positive has no such product and refuses as well.
fn product(a: f64, b: f64, c: f64) -> Result<f64, SublimationError> {
    let ab = a * b;
    let direct = ab * c;
    if ab.is_normal() && positive(direct) && direct.is_normal() {
        return Ok(direct);
    }
    if !(positive(a) && positive(b) && positive(c)) {
        return represented(f64::NAN);
    }
    represented(exact_product(a, b, c))
}

// A finite positive value as m * 2^e exactly, with an integer m < 2^53.
fn binary_parts(value: f64) -> (u64, i32) {
    let bits = value.to_bits();
    let field = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1 << 52) - 1);
    if field == 0 {
        (fraction, -1074)
    } else {
        (fraction | (1 << 52), field - 1075)
    }
}

// 2^k for -1074 <= k <= 1023.
fn power_of_two(k: i32) -> f64 {
    if k >= -1022 {
        f64::from_bits(((k + 1023) as u64) << 52)
    } else {
        f64::from_bits(1 << (k + 1074))
    }
}

// The exact product of three finite positive values, rounded once to nearest
// with ties to even, on the f64 grid down to the subnormal quantum 2^-1074:
// the result is correctly rounded. Each factor is m * 2^e with an integer
// m < 2^53, so the product is the integer ma*mb*mc < 2^159 times
// 2^(ea+eb+ec). An fma two-product holds only the 106 bits of a*b, so this
// integer is formed exactly in 192 bits. Bits past the top 128 fold into one
// sticky bit at least 74 places below the half-way bit, so the comparison
// with half the last place stays exact. The last place follows from the
// exact exponent: 53 significant bits for a normal result, 2^-1074 below
// 2^-1022. No intermediate is rounded, so a subnormal result is not rounded
// twice, and the final scaling by a power of two is exact.
fn exact_product(a: f64, b: f64, c: f64) -> f64 {
    let (ma, ea) = binary_parts(a);
    let (mb, eb) = binary_parts(b);
    let (mc, ec) = binary_parts(c);
    let mut exponent = ea + eb + ec;
    let ab = u128::from(ma) * u128::from(mb);
    let low = (ab & u128::from(u64::MAX)) * u128::from(mc);
    // The product is high * 2^64 + (low mod 2^64), with high < 2^96.
    let high = (ab >> 64) * u128::from(mc) + (low >> 64);
    let low = low as u64;
    let high_bits = 128 - high.leading_zeros() as i32;
    let m = if high_bits <= 64 {
        (high << 64) | u128::from(low)
    } else {
        let k = high_bits - 64;
        exponent += k;
        let sticky = (low & ((1 << k) - 1)) != 0;
        (high << (64 - k)) | u128::from(low >> k) | u128::from(sticky)
    };
    let width = 128 - m.leading_zeros() as i32;
    // The exact product lies in [2^top, 2^(top + 1)).
    let top = exponent + width - 1;
    if top > 1023 {
        return f64::INFINITY;
    }
    let quantum = (top - 52).max(-1074);
    let shift = quantum - exponent;
    let q = if shift <= 0 {
        m << -shift
    } else if shift > width {
        0
    } else {
        let kept = m.checked_shr(shift as u32).unwrap_or(0);
        let rest = m & (u128::MAX >> (128 - shift));
        let half = 1 << (shift - 1);
        kept + u128::from(rest > half || (rest == half && (kept & 1) == 1))
    };
    // q <= 2^53 converts exactly.
    q as f64 * power_of_two(quantum)
}

// For z >= 2, I_k(z) = integral_0^inf y^k exp(-z*y-y²/2)dy and
// r_k = I_k/I_(k-1) = k/(z+r_(k+1)). The positive terminal bounds
// 0 <= r_(N+1) <= (N+1)/z enclose each exact backward convergent. Both bounds
// must agree for ALL three ratios. At most 496 backward iterations are used;
// nonconvergence refuses instead of silently supplying an arbitrary truncation.
fn tail_ratios(z: f64) -> Result<[f64; 3], SublimationError> {
    for depth in [16, 32, 64, 128, 256] {
        let (mut lo, mut hi) = (0.0, (depth + 1) as f64 / z);
        let mut lower = [0.0; 3];
        let mut upper = [0.0; 3];
        for k in (1..=depth).rev() {
            let next_lo = k as f64 / (z + hi);
            hi = k as f64 / (z + lo);
            lo = next_lo;
            if k <= 3 {
                lower[k - 1] = lo;
                upper[k - 1] = hi;
            }
        }
        if lower
            .iter()
            .zip(upper)
            .all(|(lo, hi)| positive(*lo) && hi - lo <= 8.0 * f64::EPSILON * lo)
        {
            return Ok(core::array::from_fn(|i| 0.5 * (lower[i] + upper[i])));
        }
    }
    Err(invalid("drifting Maxwellian tail moments did not converge"))
}

struct NormalIncidence {
    mass: f64,
    mass_change: f64,
    mean_inward: f64,
    // Normal contribution to the enthalpy correction, factored to avoid
    // squaring a tiny/huge thermal speed before combining physical terms.
    energy_scale: f64,
    energy_factor: Dd,
}

fn normal_incidence(
    rho: f64,
    sigma: f64,
    inward: f64,
) -> Result<NormalIncidence, SublimationError> {
    let a = inward / sigma;
    if a >= 16.0 {
        // The missing Gaussian tail is < 1e-56 relative here, far below one
        // f64 ulp in these positive moments. Use dimensional full moments so
        // a itself may overflow for a perfectly representable cold beam.
        let ratio = sigma / inward;
        return Ok(NormalIncidence {
            mass: product(rho, inward, 1.0)?,
            mass_change: product(rho, inward, 1.0 - INV_SQRT_TAU * ratio)?,
            mean_inward: represented(inward + sigma * ratio)?,
            energy_scale: inward,
            energy_factor: Dd::from_f64(0.5),
        });
    }
    if a <= -2.0 {
        let z = -a;
        if !(0.5 * z * z).is_finite() {
            return Err(invalid(
                "nonzero incoming Maxwellian tail is unrepresentable",
            ));
        }
        let [r1, r2, r3] = tail_ratios(z)?;
        let i1 = represented(r1 / (z + r1))?;
        let log_mass = Dd::from_f64(ln(rho)) + Dd::from_f64(ln(sigma))
            - Dd::from_f64(LOG_SQRT_TAU)
            - Dd::from_f64(0.5 * z) * Dd::from_f64(z)
            + Dd::from_f64(ln(i1));
        let mass = from_log(log_mass)?;
        let baseline = product(rho, sigma, INV_SQRT_TAU)?;
        return Ok(NormalIncidence {
            mass,
            mass_change: mass - baseline,
            mean_inward: represented(sigma * r2)?,
            energy_scale: sigma,
            energy_factor: (Dd::from_f64(r2) * Dd::from_f64(r3) - Dd::from_f64(3.0))
                * Dd::from_f64(0.5),
        });
    }
    let cdf = 0.5 * erfc(-a * core::f64::consts::FRAC_1_SQRT_2);
    let phi = INV_SQRT_TAU * exp(-0.5 * a * a);
    let ad = Dd::from_f64(a);
    let a1 = Dd::from_f64(phi) + ad * Dd::from_f64(cdf);
    let a2 = ad * a1 + Dd::from_f64(cdf);
    let a3 = ad * a2 + Dd::from_f64(2.0) * a1;
    let mass_change = if inward == 0.0 {
        0.0
    } else {
        // A1(a)-A1(0) = a*[Phi(a)+phi(0)*expm1(-a²/2)/a].
        // Multiplying rho*inward retains the small linear change even when a
        // underflows or adding it to the gross rate would round it away.
        let factor = if a == 0.0 {
            0.5
        } else {
            (Dd::from_f64(cdf) + Dd::from_f64(INV_SQRT_TAU) * Dd::from_f64(expm1(-0.5 * a * a) / a))
                .to_f64()
        };
        inward.signum() * product(rho, inward.abs(), represented(factor)?)?
    };
    Ok(NormalIncidence {
        mass: product(rho, sigma, represented(a1.to_f64())?)?,
        mass_change,
        mean_inward: represented(sigma * (a2 / a1).to_f64())?,
        energy_scale: sigma,
        energy_factor: (a3 / a1 - Dd::from_f64(3.0)) * Dd::from_f64(0.5),
    })
}

/// Exact local incoming mass and conditional molecular moments for one species.
///
/// Supply the original species' partial density, specific gas constant and
/// temperature. `relative_bulk_velocity_m_s` is gas velocity MINUS local
/// surface velocity; all three translational degrees of freedom are retained.
/// The normal must be unit length within 64 machine epsilons. Inputs are fully
/// validated even in vacuum; only zero partial density yields `moments=None`.
/// Nonzero unrepresentable public moments (including the separately retained
/// drift mass increment) refuse. No arbitrary tail cutoff or Mach limit is used.
///
/// If the resulting material owner reports relative gas-side mass j, momentum
/// P and complete energy q, a surface velocity w gives laboratory-frame
/// momentum P+w*j and energy q+w·P+|w|²*j/2. These already include wall work;
/// a separate pressure-work term for the same molecules would double count it.
pub fn incident_maxwellian(
    partial_density_kg_m3: f64,
    specific_gas_constant_j_kg_k: f64,
    gas_temperature_k: f64,
    relative_bulk_velocity_m_s: [f64; 3],
    outward_unit_normal: [f64; 3],
) -> Result<MaxwellianIncidence, SublimationError> {
    if !nonnegative(partial_density_kg_m3)
        || !positive(specific_gas_constant_j_kg_k)
        || !positive(gas_temperature_k)
        || relative_bulk_velocity_m_s.iter().any(|x| !x.is_finite())
        || outward_unit_normal
            .iter()
            .any(|x| !x.is_finite() || x.abs() > 1.0 + NORMAL_TOLERANCE)
        || (dot(outward_unit_normal, outward_unit_normal)? - 1.0).abs() > NORMAL_TOLERANCE
    {
        return Err(invalid(
            "finite species state, relative velocity and unit outward normal required",
        ));
    }
    if partial_density_kg_m3 == 0.0 {
        return Ok(MaxwellianIncidence {
            incoming_mass_flux_kg_m2_s: 0.0,
            incoming_mass_flux_change_from_zero_drift_kg_m2_s: 0.0,
            moments: None,
        });
    }
    let sigma = represented(sqrt(specific_gas_constant_j_kg_k) * sqrt(gas_temperature_k))?;
    let normal_velocity = dot(relative_bulk_velocity_m_s, outward_unit_normal)?;
    let mut tangent = [0.0; 3];
    for i in 0..3 {
        tangent[i] = (Dd::from_f64(relative_bulk_velocity_m_s[i])
            - Dd::from_f64(normal_velocity) * Dd::from_f64(outward_unit_normal[i]))
        .to_f64();
    }
    let normal = normal_incidence(partial_density_kg_m3, sigma, -normal_velocity)?;
    let mut velocity = [0.0; 3];
    for i in 0..3 {
        velocity[i] = (Dd::from_f64(tangent[i])
            - Dd::from_f64(normal.mean_inward) * Dd::from_f64(outward_unit_normal[i]))
        .to_f64();
        if !tangent[i].is_finite() || !velocity[i].is_finite() {
            return Err(invalid(
                "conditional incoming molecular velocity is unrepresentable",
            ));
        }
    }
    let largest = tangent
        .into_iter()
        .map(f64::abs)
        .fold(normal.energy_scale, f64::max);
    let scale = binary_scale(largest);
    let ratio = Dd::from_f64(normal.energy_scale / scale);
    let mut energy = normal.energy_factor * ratio * ratio;
    for value in tangent {
        let value = Dd::from_f64(value / scale);
        energy = energy + Dd::from_f64(0.5) * value * value;
    }
    let correction = (energy * Dd::from_f64(scale) * Dd::from_f64(scale)).to_f64();
    if !correction.is_finite() || (energy != Dd::ZERO && correction == 0.0) {
        return Err(invalid(
            "incoming molecular crossing energy correction is unrepresentable",
        ));
    }
    Ok(MaxwellianIncidence {
        incoming_mass_flux_kg_m2_s: normal.mass,
        incoming_mass_flux_change_from_zero_drift_kg_m2_s: normal.mass_change,
        moments: Some(IncidentMolecularMoments {
            mean_relative_velocity_m_s: velocity,
            crossing_energy_minus_enthalpy_j_kg: correction,
        }),
    })
}

#[cfg(test)]
mod tests;
