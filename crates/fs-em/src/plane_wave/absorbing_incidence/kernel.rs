//! Fixed-size Maxwell boundary solve with explicit incident interference flux.
use super::super::{
    ReflectionError, WaveMedium, checked_outgoing_normal, finite, normalized_interface_pair,
    positive,
};
use fs_math::{c64::C64, dd::Dd};

type Result<T> = core::result::Result<T, ReflectionError>;

/// Full Cartesian complex peak field at the interface, with ky=0.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InterfaceBoundaryField {
    pub electric_v_m: [C64; 3],
    pub magnetic_a_m: [C64; 3],
}

/// Actual interface fields and unnormalized cycle-mean normal power [W/m²].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbsorbingInterfaceReflection {
    pub transverse_wavenumber_rad_m: f64,
    pub incident_normal_wavenumber_rad_m: C64,
    pub transmitted_normal_wavenumber_rad_m: C64,
    pub incident_permittivity_f_m: C64,
    pub transmitted_permittivity_f_m: C64,
    pub incident_permeability_h_m: f64,
    pub transmitted_permeability_h_m: f64,
    /// Fixed tangential Ex (TM), then Ey (TE), on both sides.
    pub reflected_tangential_e_gain: [C64; 2],
    pub transmitted_tangential_e_gain: [C64; 2],
    pub incident: InterfaceBoundaryField,
    pub reflected: InterfaceBoundaryField,
    pub transmitted: InterfaceBoundaryField,
    /// Positive into +z, or zero for zero illumination.
    pub incident_normal_power_w_m2: f64,
    /// Signed along +z: outgoing reflected power is nonpositive.
    pub reflected_normal_power_w_m2: f64,
    /// Signed cross term between the incident and reflected fields. It may
    /// have either sign; it is neither interface absorption nor an extra ray.
    pub interference_normal_power_w_m2: f64,
    /// Flux entering z>0. Subsequent bulk absorption is separate.
    pub transmitted_normal_power_w_m2: f64,
}

#[derive(Clone, Copy)]
pub(super) struct Input {
    pub(super) omega: f64,
    pub(super) kx: f64,
    pub(super) epsilon_i: C64,
    pub(super) epsilon_t: C64,
    pub(super) mu_i: f64,
    pub(super) mu_t: f64,
    pub(super) electric: [C64; 2],
}

fn invalid(message: &'static str) -> ReflectionError {
    ReflectionError::Invalid(message)
}
fn magnitude(z: C64) -> f64 {
    z.re.abs().max(z.im.abs())
}
fn divide_real(z: C64, divisor: f64) -> C64 {
    C64::new(z.re / divisor, z.im / divisor)
}
fn contrast(z: C64, divisor: f64) -> Result<C64> {
    let value = divide_real(z, divisor);
    if !finite(value) || (z.re != 0.0 && value.re == 0.0) || (z.im != 0.0 && value.im == 0.0) {
        return Err(invalid("interface scattering contrast is unrepresentable"));
    }
    Ok(value)
}
fn power_of_two(exponent: i32) -> f64 {
    f64::from_bits(((exponent + 1023) as u64) << 52)
}
fn exponent(value: f64) -> i32 {
    (((value.to_bits() >> 52) & 0x7ff) as i32 - 1023).max(-1022)
}

// Real, exact binary reference coordinates, not a substitute material law.
// The product and ratio of these constitutive scales have exact square roots.
pub(in crate::plane_wave) fn reference_scales(
    omega: f64,
    epsilon: C64,
    mu: f64,
) -> Result<[f64; 4]> {
    let mut e = exponent(magnitude(epsilon));
    let mut m = exponent(mu);
    if (e + m) & 1 != 0 {
        if e > -1022 {
            e -= 1;
        } else {
            m -= 1;
        }
    }
    Ok([
        positive(power_of_two(e))?,
        positive(power_of_two(m))?,
        positive(omega * power_of_two((e + m) / 2))?,
        positive(power_of_two((e - m) / 2))?,
    ])
}

#[derive(Clone, Copy)]
pub(in crate::plane_wave) struct Medium {
    pub(in crate::plane_wave) normal: C64,
    pub(in crate::plane_wave) omega_epsilon: C64,
    pub(in crate::plane_wave) omega_mu: f64,
    pub(in crate::plane_wave) te_admittance: C64,
    // Unlike TM admittance, this remains finite at a critical exit.
    pub(in crate::plane_wave) tm_impedance: C64,
}
impl Medium {
    pub(in crate::plane_wave) fn new(omega: f64, kx: f64, epsilon: C64, mu: f64) -> Result<Self> {
        positive(mu)?;
        if !finite(epsilon) || epsilon == C64::ZERO || epsilon.im > 0.0 {
            return Err(invalid(
                "finite passive nonzero interface permittivity required",
            ));
        }
        // Binary constitutive references make their geometric mean exact.
        // No rounded refractive index/angle precedes the k_z discriminant;
        // omega² and the unscaled epsilon*mu product are never formed.
        let [epsilon_ref, mu_ref, wave_scale, _] = reference_scales(omega, epsilon, mu)?;
        let medium = WaveMedium::new(epsilon, mu, epsilon_ref, mu_ref)?;
        if (epsilon.re != 0.0 && medium.index_squared_ratio.re == 0.0)
            || (epsilon.im != 0.0 && medium.index_squared_ratio.im == 0.0)
        {
            return Err(invalid("interface constitutive contrast underflows"));
        }
        let transverse = Dd::from_f64(kx) / Dd::from_f64(wave_scale);
        let square = Dd::from_pair(
            medium.index_squared_ratio.re,
            medium.index_squared_real_residual,
        ) - transverse * transverse;
        if !square.hi.is_finite() || !square.lo.is_finite() {
            return Err(invalid("interface normal discriminant is unrepresentable"));
        }
        let normalized = checked_outgoing_normal(
            C64::new(square.to_f64(), medium.index_squared_ratio.im).sqrt(),
        )?;
        let normal = normalized.scale(wave_scale);
        let omega_mu = positive(omega * mu)?;
        let omega_epsilon = epsilon.scale(omega);
        if !finite(normal)
            || !finite(omega_epsilon)
            || omega_epsilon == C64::ZERO
            || (normalized.re != 0.0 && normal.re == 0.0)
            || (normalized.im != 0.0 && normal.im == 0.0)
            || (epsilon.re != 0.0 && omega_epsilon.re == 0.0)
            || (epsilon.im != 0.0 && omega_epsilon.im == 0.0)
        {
            return Err(invalid(
                "dimensional interface wave data are unrepresentable",
            ));
        }
        let te_admittance = contrast(normal, omega_mu)?;
        let tm_impedance = normal / omega_epsilon;
        if !finite(tm_impedance) || (normal.re > 0.0 && tm_impedance.re == 0.0) {
            return Err(invalid("interface wave impedance is unrepresentable"));
        }
        Ok(Self {
            normal,
            omega_epsilon,
            omega_mu,
            te_admittance,
            tm_impedance,
        })
    }

    pub(in crate::plane_wave) fn field(
        self,
        kx: f64,
        e: [C64; 2],
        h: [C64; 2],
    ) -> Result<InterfaceBoundaryField> {
        // h=[Hy,-Hx] is power-conjugate to e=[Ex,Ey]. Maxwell-Ampere
        // determines Ez without dividing by a vanishing transmitted k_z.
        let field = InterfaceBoundaryField {
            electric_v_m: [e[0], e[1], -(h[0] / self.omega_epsilon).scale(kx)],
            magnetic_a_m: [-h[1], h[0], e[1].scale(kx / self.omega_mu)],
        };
        if field
            .electric_v_m
            .iter()
            .chain(&field.magnetic_a_m)
            .any(|&z| !finite(z))
        {
            return Err(invalid("interface field is unrepresentable"));
        }
        Ok(field)
    }
}

fn gains(a: C64, b: C64) -> Result<[C64; 3]> {
    let (a, b, denominator, _) = normalized_interface_pair(a, b)?;
    let result = [
        (a - b) / denominator,
        a.scale(2.0) / denominator,
        b.scale(2.0) / denominator,
    ];
    if result.iter().any(|&z| !finite(z)) {
        return Err(invalid("interface scattering amplitude is unrepresentable"));
    }
    Ok(result)
}

// Apply the weight before squaring an amplitude. This avoids needlessly
// overflowing |E|² when its physical power is representable, and keeps a
// weak interference term independent of subtraction of the larger powers.
pub(in crate::plane_wave) fn quadratic(weight: Dd, amplitude: C64) -> Dd {
    let re = Dd::from_f64(amplitude.re);
    let im = Dd::from_f64(amplitude.im);
    (weight * re) * re + (weight * im) * im
}

pub(in crate::plane_wave) fn entrance_power(
    yi: [C64; 2],
    reflected_gain: [C64; 2],
    incident: [C64; 2],
    reflected: [C64; 2],
) -> [Dd; 3] {
    let mut pi = Dd::ZERO;
    let mut pr = Dd::ZERO;
    let mut interference = Dd::ZERO;
    for c in 0..2 {
        let weight = Dd::from_f64(yi[c].re) * Dd::from_f64(0.5);
        pi = pi + quadratic(weight, incident[c]);
        pr = pr - quadratic(weight, reflected[c]);
        interference = interference
            + quadratic(
                Dd::from_f64(yi[c].im) * Dd::from_f64(reflected_gain[c].im),
                incident[c],
            );
    }
    [pi, pr, interference]
}

pub(in crate::plane_wave) fn exit_power(medium: Medium, e: [C64; 2], h: [C64; 2]) -> Dd {
    quadratic(
        Dd::from_f64(medium.te_admittance.re) * Dd::from_f64(0.5),
        e[1],
    ) + quadratic(
        Dd::from_f64(medium.tm_impedance.re) * Dd::from_f64(0.5),
        h[0],
    )
}

pub(super) fn solve(data: Input) -> Result<AbsorbingInterfaceReflection> {
    let omega = positive(data.omega)?;
    if !data.kx.is_finite() || data.electric.iter().any(|&z| !finite(z)) {
        return Err(invalid(
            "finite transverse wave number and incident electric field required",
        ));
    }
    let incident = Medium::new(omega, data.kx, data.epsilon_i, data.mu_i)?;
    let transmitted = Medium::new(omega, data.kx, data.epsilon_t, data.mu_t)?;
    if incident.normal == C64::ZERO || incident.tm_impedance == C64::ZERO {
        return Err(invalid(
            "incident interface wave must carry positive normal power",
        ));
    }
    let yi = [C64::ONE / incident.tm_impedance, incident.te_admittance];
    if yi.iter().any(|&y| !finite(y) || y.re <= 0.0) {
        return Err(invalid(
            "incident interface wave must carry positive normal power",
        ));
    }
    if transmitted.te_admittance.re < 0.0 || transmitted.tm_impedance.re < 0.0 {
        return Err(invalid("transmitted interface wave is not passive"));
    }
    let q_scale = positive(magnitude(incident.normal).max(magnitude(transmitted.normal)))?;
    let qi = contrast(incident.normal, q_scale)?;
    let qt = contrast(transmitted.normal, q_scale)?;
    let mu_scale = positive(data.mu_i.max(data.mu_t))?;
    let mu_i = positive(data.mu_i / mu_scale)?;
    let mu_t = positive(data.mu_t / mu_scale)?;
    let epsilon_scale = positive(magnitude(data.epsilon_i).max(magnitude(data.epsilon_t)))?;
    let eps_i = contrast(data.epsilon_i, epsilon_scale)?;
    let eps_t = contrast(data.epsilon_t, epsilon_scale)?;
    // a/b=Y_i/Y_t. Multiplication through by q_t retains exact critical TM
    // fields; scaling both wave numbers and constitutive values bounds a,b.
    let te = gains(qi.scale(mu_t), qt.scale(mu_i))?;
    let tm = gains(eps_i * qt, eps_t * qi)?;
    let r = [tm[0], te[0]];
    let t = [tm[1], te[1]];
    let h_gain = [tm[2], te[2]];
    let ei = data.electric;
    let hi = [yi[0] * ei[0], yi[1] * ei[1]];
    let er = [r[0] * ei[0], r[1] * ei[1]];
    let hr = [-yi[0] * er[0], -yi[1] * er[1]];
    let et = [t[0] * ei[0], t[1] * ei[1]];
    // The independent magnetic numerator survives even when 1-r would round
    // to zero; no transmitted critical TM admittance is evaluated.
    let ht = [h_gain[0] * hi[0], h_gain[1] * hi[1]];
    let incident_field = incident.field(data.kx, ei, hi)?;
    let reflected_field = incident.field(data.kx, er, hr)?;
    let transmitted_field = transmitted.field(data.kx, et, ht)?;
    let [pi, pr, interference] = entrance_power(yi, r, ei, er);
    let pt = exit_power(transmitted, et, ht);
    let powers = [pi, pr, interference, pt].map(Dd::to_f64);
    if powers.iter().any(|p| !p.is_finite())
        || powers[0] < 0.0
        || powers[1] > 0.0
        || powers[3] < 0.0
    {
        return Err(invalid("interface power is unrepresentable or nonpassive"));
    }
    let illuminated = ei.iter().any(|&e| e != C64::ZERO);
    if illuminated && powers[0] == 0.0 {
        return Err(invalid("incident interface power underflows"));
    }
    let scale = powers.iter().fold(0.0_f64, |s, p| s.max(p.abs()));
    if scale > 0.0 {
        let scale = Dd::from_f64(scale);
        let defect = (pi / scale + pr / scale + interference / scale - pt / scale).to_f64();
        if !defect.is_finite() || defect.abs() > 256.0 * f64::EPSILON {
            return Err(invalid(
                "interface Poynting balance with interference failed",
            ));
        }
    }
    Ok(AbsorbingInterfaceReflection {
        transverse_wavenumber_rad_m: data.kx,
        incident_normal_wavenumber_rad_m: incident.normal,
        transmitted_normal_wavenumber_rad_m: transmitted.normal,
        incident_permittivity_f_m: data.epsilon_i,
        transmitted_permittivity_f_m: data.epsilon_t,
        incident_permeability_h_m: data.mu_i,
        transmitted_permeability_h_m: data.mu_t,
        reflected_tangential_e_gain: r,
        transmitted_tangential_e_gain: t,
        incident: incident_field,
        reflected: reflected_field,
        transmitted: transmitted_field,
        incident_normal_power_w_m2: powers[0],
        reflected_normal_power_w_m2: powers[1],
        interference_normal_power_w_m2: powers[2],
        transmitted_normal_power_w_m2: powers[3],
    })
}

#[cfg(test)]
mod tests;
