//! Principal-plane uniaxial Maxwell interface, without material-source policy.
//!
//! Time convention: exp(+i omega t), forward exp(-i k.z). The optic axis
//! lies in the x-z incidence plane; z points into the crystal. TE excites
//! the ordinary mode, TM the extraordinary mode. Tangential E uses the same
//! fixed x/y axes on both sides, including reflection.
//!
//! For a unit optic axis a and q=k/k0, the extraordinary dispersion is
//! q^T M q=1, M=(I-aa^T)/ne^2+aa^T/no^2. Its energy direction is Mq,
//! not q. The outgoing root is selected by positive normal energy flux.
//! The evanescent root has Im(kz)<0 and therefore decays into z>0.
//! Reference: F. Kaertner, MIT 6.974, section 2.8, equations 2.410--2.419:
//! https://ocw.mit.edu/courses/6-974-fundamentals-of-photonics-quantum-electronics-spring-2006/8e54d1625e9e71eb5b4e4998e6967e4e_chapter2.pdf

use fs_math::{c64::C64, dd::Dd, det};

pub(crate) const LIGHT_SPEED_M_S: f64 = 299_792_458.0;

/// Decoupled transmitted eigenmode for an optic axis in the incidence plane.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UniaxialMode {
    /// Electric field perpendicular to the incidence plane (TE).
    Ordinary,
    /// Magnetic field perpendicular to the incidence plane (TM).
    Extraordinary,
}

/// A bulk ray, a grazing critical field, or a decaying near field.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UniaxialWaveKind {
    Propagating,
    /// Exactly zero normal transmitted power. No ray enters the half-space.
    Critical,
    /// Normal decay with zero transmitted normal power; not a bulk ray.
    Evanescent,
}

/// Signed fields, normal Poynting fractions and separate phase/energy directions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniaxialInterfaceResponse {
    pub mode: UniaxialMode,
    pub wave_kind: UniaxialWaveKind,
    /// Reflected tangential electric field / incident tangential electric field.
    pub reflected_tangential_e: C64,
    /// Transmitted tangential E / incident tangential E, formed independently.
    pub transmitted_tangential_e: C64,
    /// Transmitted power-conjugate H / incident H, formed independently.
    /// Critical TM keeps finite H even when its tangential E is exactly zero.
    pub transmitted_conjugate_h: C64,
    pub reflected_power_fraction: f64,
    pub transmitted_power_fraction: f64,
    /// Conserved signed tangential wavenumber [rad/m].
    pub tangential_wavenumber_rad_m: f64,
    /// Positive real incident normal wavenumber [rad/m].
    pub incident_normal_wavenumber_rad_m: f64,
    /// Outgoing/decaying crystal normal wavenumber [rad/m].
    /// Its real part can be negative: energy flux selects the outgoing branch.
    pub transmitted_normal_wavenumber_rad_m: C64,
    /// Unit [x,z] wave-normal direction for a propagating or critical mode.
    /// Evanescent fields are inhomogeneous waves and return None.
    pub phase_direction: Option<[f64; 2]>,
    /// Unit [x,z] time-averaged Poynting direction, including crystal walk-off.
    /// Critical modes are tangent; evanescent fields do not supply a bulk ray.
    pub energy_direction: Option<[f64; 2]>,
}

/// Preserve the exact supplied coordinate: grazing cosines must not be lost
/// by rounding sin(theta) to one, and exact tangential critical data must not
/// be perturbed by conversion to an angle.
#[derive(Clone, Copy)]
pub(crate) enum Incidence {
    Cosine(f64),
    /// Signed kx/k0, with absolute value strictly below the incident phase index.
    TangentialIndex(f64),
}

#[derive(Clone, Copy)]
pub(crate) struct Input {
    pub omega_rad_s: f64,
    pub incident_index: f64,
    pub incident_mu_h_m: f64,
    pub ordinary_index: f64,
    pub extraordinary_index: f64,
    pub crystal_mu_h_m: f64,
    pub optic_axis: [f64; 2],
    pub incidence: Incidence,
    pub mode: UniaxialMode,
}

pub(crate) fn unit_axis(axis: [f64; 2]) -> Result<[f64; 2], &'static str> {
    if axis.iter().any(|v| !v.is_finite()) {
        return Err("principal-plane optic axis must be finite and nonzero");
    }
    let scale = axis[0].abs().max(axis[1].abs());
    if scale == 0.0 {
        return Err("principal-plane optic axis must be finite and nonzero");
    }
    let x = axis[0] / scale;
    let z = axis[1] / scale;
    let length = det::hypot(x, z);
    Ok([canonical_zero(x / length), canonical_zero(z / length)])
}

fn canonical_zero(x: f64) -> f64 {
    if x == 0.0 { 0.0 } else { x }
}

fn positive(x: f64) -> Result<f64, &'static str> {
    if x.is_finite() && x > 0.0 {
        Ok(x)
    } else {
        Err("uniaxial wave data must be positive, finite and representable")
    }
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

// A power-of-two scale does not round normal index values. Its square is
// never formed: all squared-index work happens at order-one magnitude.
pub(super) fn binary_scale(x: f64) -> f64 {
    let exponent = x.to_bits() & 0x7ff0_0000_0000_0000;
    if exponent == 0 {
        f64::MIN_POSITIVE
    } else {
        f64::from_bits(exponent)
    }
}

pub(super) fn outgoing_root(value: Dd) -> Result<(C64, UniaxialWaveKind), &'static str> {
    if !value.hi.is_finite() || !value.lo.is_finite() {
        return Err("uniaxial normal discriminant is unrepresentable");
    }
    if value == Dd::ZERO {
        return Ok((C64::ZERO, UniaxialWaveKind::Critical));
    }
    let magnitude = value.abs().sqrt().to_f64();
    positive(magnitude)?;
    if value.lt(Dd::ZERO) {
        Ok((C64::new(0.0, -magnitude), UniaxialWaveKind::Evanescent))
    } else {
        Ok((C64::from_re(magnitude), UniaxialWaveKind::Propagating))
    }
}

fn direction([x, z]: [f64; 2]) -> Result<[f64; 2], &'static str> {
    unit_axis([x, z]).map_err(|_| "uniaxial phase or energy direction is unrepresentable")
}

pub(crate) fn solve(input: Input) -> Result<UniaxialInterfaceResponse, &'static str> {
    let omega = positive(input.omega_rad_s)?;
    let ni = positive(input.incident_index)?;
    let no = positive(input.ordinary_index)?;
    let ne = positive(input.extraordinary_index)?;
    let mui = positive(input.incident_mu_h_m)?;
    let mut_ = positive(input.crystal_mu_h_m)?;
    let axis = unit_axis(input.optic_axis)?;

    let scale = binary_scale(ni.max(no).max(ne));
    let i = positive(ni / scale)?;
    let o = positive(no / scale)?;
    let e = positive(ne / scale)?;
    let id = Dd::from_f64(i);
    let od = Dd::from_f64(o);
    let ed = Dd::from_f64(e);
    let i2 = id * id;
    let o2 = od * od;
    let e2 = ed * ed;
    positive(i2.hi)?;
    positive(o2.hi)?;
    positive(e2.hi)?;
    let (p, p2, incident_normal, cosine) = match input.incidence {
        Incidence::Cosine(c) => {
            positive(c)?;
            if c > 1.0 {
                return Err("uniaxial incidence cosine must lie in (0,1]");
            }
            let cd = Dd::from_f64(c);
            let s2 = (Dd::ONE - cd) * (Dd::ONE + cd);
            let p2 = i2 * s2;
            let p = i * det::sqrt(s2.to_f64());
            (p, p2, positive(i * c)?, c)
        }
        Incidence::TangentialIndex(value) => {
            if !value.is_finite() || value.abs() >= ni {
                return Err("tangential index needs a propagating incident port");
            }
            let p = value / scale;
            if value != 0.0 && p == 0.0 {
                return Err("uniaxial tangential index is unrepresentable");
            }
            let pd = Dd::from_f64(p);
            let p2 = pd * pd;
            let q = positive((i2 - p2).sqrt().to_f64())?;
            (p, p2, q, positive(q / i)?)
        }
    };

    // Re-normalize squared direction weights in double-double arithmetic.
    // This keeps the isotropic limit exact even if the unit vector's rounded
    // squared length differs from one by an ulp.
    let ax = Dd::from_f64(axis[0]);
    let az = Dd::from_f64(axis[1]);
    let axis_norm2 = ax * ax + az * az;
    let xx = ax * ax / axis_norm2;
    let zz = az * az / axis_norm2;
    let xz = ax * az / axis_norm2;
    let f = if no == ne { o2 } else { o2 * xx + e2 * zz };
    let h = if no == ne { o2 } else { o2 * zz + e2 * xx };
    let g = (e2 - o2) * xz;
    let oe_dd = od * ed;
    let oe = positive(oe_dd.to_f64())?;

    let discriminant = match input.mode {
        UniaxialMode::Ordinary => o2 - p2,
        UniaxialMode::Extraordinary => f - p2,
    };
    // A matched incident/crystal dispersion has qz=ni*cos(theta) exactly.
    // Retain it even when cos² underflows or 1-cos² rounds to one.
    let exactly_matched = match input.mode {
        UniaxialMode::Ordinary => ni == no,
        UniaxialMode::Extraordinary => ni == no && no == ne,
    };
    let matched_cutoff = match input.mode {
        UniaxialMode::Ordinary => ni == no,
        UniaxialMode::Extraordinary => f == i2,
    };
    let (root, kind) = if matched_cutoff {
        (C64::from_re(incident_normal), UniaxialWaveKind::Propagating)
    } else {
        outgoing_root(discriminant)?
    };

    let (normal, energy) = match input.mode {
        UniaxialMode::Ordinary => (root, [p, root.re]),
        UniaxialMode::Extraordinary if no == ne => (root, [p, root.re]),
        UniaxialMode::Extraordinary => {
            let fd = positive(f.to_f64())?;
            let pd = Dd::from_f64(p);
            let shift_numerator = g * pd;
            let term = oe_dd * Dd::from_f64(root.re);
            let real = if kind == UniaxialWaveKind::Propagating && Dd::ZERO.lt(shift_numerator) {
                // Stable outgoing quadratic root when -Gp+oe*sqrt(D)
                // would subtract nearly equal numbers. The other root is
                // never used as an outgoing-wave criterion.
                let constant = h * p2 - o2 * e2;
                (constant / (-shift_numerator - term)).to_f64()
            } else {
                ((-shift_numerator + term) / f).to_f64()
            };
            let normal = C64::new(real, (oe / fd) * root.im);
            // Equivalent to Mq up to a positive common scale. The normal
            // component is independently sqrt(D); no cancellation can flip
            // its sign at grazing or for a negative normal phase velocity.
            let tangent_energy = ((oe_dd * pd + g * Dd::from_f64(root.re)) / f).to_f64();
            (normal, [tangent_energy, root.re])
        }
    };
    if !finite(normal) || normal.im > 0.0 {
        return Err("outgoing uniaxial normal wavenumber is unrepresentable");
    }

    let mu_scale = binary_scale(mui.max(mut_));
    let mi = positive(mui / mu_scale)?;
    let mt = positive(mut_ / mu_scale)?;
    // a/b = Yi/Yt, with singular admittances multiplied out. Ordinary is TE;
    // extraordinary is TM. Both field numerators survive critical incidence.
    let (a, b) = if exactly_matched {
        (C64::from_re(mt), C64::from_re(mi))
    } else {
        match input.mode {
            UniaxialMode::Ordinary => (C64::from_re(mt * incident_normal), root.scale(mi)),
            UniaxialMode::Extraordinary => (root.scale(mt * i), C64::from_re(mi * oe * cosine)),
        }
    };
    let field_scale = positive(a.re.abs().max(a.im.abs()).max(b.re.abs()).max(b.im.abs()))?;
    let a = C64::new(a.re / field_scale, a.im / field_scale);
    let b = C64::new(b.re / field_scale, b.im / field_scale);
    if kind == UniaxialWaveKind::Propagating && (a.re <= 0.0 || b.re <= 0.0) {
        return Err("positive uniaxial transmitted flux is unrepresentable");
    }
    let denominator = a + b;
    let denominator_power = positive(denominator.norm_sq())?;
    let reflected = (a - b) / denominator;
    let transmitted = a.scale(2.0) / denominator;
    let transmitted_h = b.scale(2.0) / denominator;
    let reflected_power = reflected.norm_sq();
    let transmitted_power = if kind == UniaxialWaveKind::Propagating {
        // Independent positive normal Poynting flux, not 1-|r|². After the
        // scale one real coefficient is one, preserving weak transmission.
        (4.0 * a.re * b.re) / denominator_power
    } else {
        0.0
    };
    if !finite(reflected)
        || !finite(transmitted)
        || !finite(transmitted_h)
        || !reflected_power.is_finite()
        || !transmitted_power.is_finite()
        || transmitted_power < 0.0
        || (reflected_power + transmitted_power - 1.0).abs() > 256.0 * f64::EPSILON
    {
        return Err("uniaxial Maxwell normal-power balance failed");
    }

    let wave_scale = positive((omega / LIGHT_SPEED_M_S) * scale)?;
    let tangential = wave_scale * p;
    let incident_normal = positive(wave_scale * incident_normal)?;
    let transmitted_normal = normal.scale(wave_scale);
    if !tangential.is_finite()
        || !finite(transmitted_normal)
        || p != 0.0 && tangential == 0.0
        || normal.re != 0.0 && transmitted_normal.re == 0.0
        || normal.im != 0.0 && transmitted_normal.im == 0.0
    {
        return Err("uniaxial SI wavenumber is unrepresentable");
    }
    let (phase_direction, energy_direction) = if kind == UniaxialWaveKind::Evanescent {
        (None, None)
    } else {
        (Some(direction([p, normal.re])?), Some(direction(energy)?))
    };
    Ok(UniaxialInterfaceResponse {
        mode: input.mode,
        wave_kind: kind,
        reflected_tangential_e: reflected,
        transmitted_tangential_e: transmitted,
        transmitted_conjugate_h: transmitted_h,
        reflected_power_fraction: reflected_power,
        transmitted_power_fraction: transmitted_power,
        tangential_wavenumber_rad_m: canonical_zero(tangential),
        incident_normal_wavenumber_rad_m: incident_normal,
        transmitted_normal_wavenumber_rad_m: transmitted_normal,
        phase_direction,
        energy_direction,
    })
}


/// Shared principal/arbitrary-axis contraction of the original spectral slopes.
/// Weights are |q_parallel/no|² and |q_perpendicular/ne|². The ng/n ratios
/// must be divided in double-double arithmetic before the signed sum: f64
/// division can erase a finite difference between nearly opposite slopes.
pub(super) fn dispersive_group_denominator(
    weights: [Dd; 2],
    indices: [f64; 2],
    group_indices: [f64; 2],
) -> Result<Option<f64>, &'static str> {
    let ordinary = Dd::from_f64(group_indices[0]) / Dd::from_f64(indices[0]);
    let extraordinary = Dd::from_f64(group_indices[1]) / Dd::from_f64(indices[1]);
    let denominator = weights[0] * ordinary + weights[1] * extraordinary;
    if denominator == Dd::ZERO {
        return Ok(None);
    }
    let denominator = denominator.to_f64();
    if !denominator.is_finite() || denominator == 0.0 {
        return Err("uniaxial dispersive denominator is unrepresentable");
    }
    Ok(Some(denominator))
}

/// Differentiate the ORIGINAL ordinary/extraordinary bulk dispersion surface.
/// The carrier remains the native outgoing interface mode; there is no second
/// root solve or angular finite difference here.
///
/// For the extraordinary wave, F=k^T M(omega) k-omega²/c0²=0, with
/// M=(I-aa^T)/ne²+aa^T/no². Thus v_group=-F_k/F_omega. Write
/// u=(q.a)/no, v=(q.b)/ne, q=c0*k/omega and b perpendicular to a.
/// Then v_group=c0*M*q/[u²*ng_o/no+v²*ng_e/ne].
/// Native Poynting direction supplies the robust orientation at critical
/// incidence and negative normal phase velocity. Positive dispersive denominator
/// aligns group and Poynting directions; an anomalous signed slope is not
/// promoted to a positive energy or signal velocity.
/// MIT 6.974 section 2.8.3, equations 2.412--2.414, plus differentiation of
/// the original frequency-dependent principal indices.
pub(crate) fn group_velocity(
    carrier: UniaxialInterfaceResponse,
    omega_rad_s: f64,
    indices: [f64; 2],
    group_indices: [f64; 2],
    optic_axis: [f64; 2],
) -> Result<Option<[f64; 2]>, &'static str> {
    let omega = positive(omega_rad_s)?;
    let [no, ne] = indices.map(positive);
    let (no, ne) = (no?, ne?);
    if group_indices.iter().any(|value| !value.is_finite()) {
        return Err("uniaxial group indices must be finite");
    }
    let Some(direction) = carrier.energy_direction else {
        // An evanescent field is not a bulk wave packet travelling into z>0.
        return Ok(None);
    };
    let speed = if carrier.mode == UniaxialMode::Ordinary {
        if group_indices[0] == 0.0 {
            return Ok(None);
        }
        LIGHT_SPEED_M_S / group_indices[0]
    } else {
        let axis = unit_axis(optic_axis)?;
        let k0 = positive(omega / LIGHT_SPEED_M_S)?;
        let q = [
            carrier.tangential_wavenumber_rad_m / k0,
            carrier.transmitted_normal_wavenumber_rad_m.re / k0,
        ];
        let along = (axis[0] * q[0] + axis[1] * q[1]) / no;
        let across = (-axis[1] * q[0] + axis[0] * q[1]) / ne;
        if !q.into_iter().chain([along, across]).all(f64::is_finite) {
            return Err("uniaxial group-wave coordinates are unrepresentable");
        }
        let weights = [
            Dd::from_f64(along) * Dd::from_f64(along),
            Dd::from_f64(across) * Dd::from_f64(across),
        ];
        let Some(denominator) =
            dispersive_group_denominator(weights, [no, ne], group_indices)?
        else {
            return Ok(None);
        };
        let metric_q = [
            axis[0] * (along / no) - axis[1] * (across / ne),
            axis[1] * (along / no) + axis[0] * (across / ne),
        ];
        let magnitude = positive(det::hypot(metric_q[0], metric_q[1]))?;
        (magnitude / denominator) * LIGHT_SPEED_M_S
    };
    let velocity = direction.map(|component| component * speed);
    if !speed.is_finite()
        || speed == 0.0
        || velocity.iter().any(|value| !value.is_finite())
        || velocity == [0.0; 2]
    {
        return Err("uniaxial group velocity is unrepresentable");
    }
    Ok(Some(velocity))
}

#[cfg(test)]
#[path = "kernel_tests.rs"]
mod tests;
