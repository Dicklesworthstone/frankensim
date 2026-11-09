//! Analytic normal-wavevector jets of the original uniaxial quadratic.
//!
//! All index work uses one fixed binary scale. Coefficients and contractions
//! retain double-double residuals; frequency derivatives are scaled by omega
//! until the final dimensional divisions. No neighboring carrier is solved.

use super::super::kernel::{Incidence, UniaxialMode, UniaxialWaveKind, binary_scale, unit_axis};
use crate::plane_wave::dispersion::IncidenceFrequencyConstraint;
use fs_math::{c64::C64, dd::Dd, det};

const C0: f64 = 299_792_458.0;
type Result<T> = core::result::Result<T, &'static str>;

/// Already source-admitted frequency data. The first spectral array is
/// [incident, ordinary, extraordinary] [omega*epsilon'/epsilon,
///  omega^2*epsilon''/epsilon]. The second keeps the original principal
/// [K', omega*K''] combinations, so zero bulk GVD is not reconstructed by
/// subtracting large index-derivative terms.
pub(in crate::plane_wave::uniaxial) struct Input {
    pub omega: f64,
    pub indices: [f64; 3],
    pub relative_permittivity_derivatives: [[f64; 2]; 3],
    pub principal_wavenumber_derivatives: [[f64; 2]; 2],
    pub incidence: Incidence,
    pub axis: [f64; 3],
    pub mode: UniaxialMode,
    pub kind: UniaxialWaveKind,
    pub normal_wavenumber: C64,
    pub constraint: IncidenceFrequencyConstraint,
}

pub(in crate::plane_wave::uniaxial) fn normal_derivatives(data: Input) -> Result<[C64; 2]> {
    let omega = positive(data.omega)?;
    let [ni, no, ne] = data.indices;
    for value in data.indices {
        positive(value)?;
    }
    if data
        .relative_permittivity_derivatives
        .iter()
        .flatten()
        .chain(data.principal_wavenumber_derivatives.iter().flatten())
        .any(|value| !value.is_finite())
        || data.axis.iter().any(|value| !value.is_finite())
        || !finite(data.normal_wavenumber)
    {
        return Err("finite original crystal frequency jets required");
    }
    // A tilted critical mode may have nonzero normal PHASE. Conversely a
    // propagating mode may have kz=0 while its normal energy remains positive.
    if data.kind == UniaxialWaveKind::Critical {
        return Err("critical crystal modes have no admitted normal frequency jet");
    }
    let scale = binary_scale(ni.max(no).max(ne));
    let kscale = positive((omega / C0) * scale)?;
    let q = Z::from_c64(data.normal_wavenumber).scale(Dd::ONE / Dd::from_f64(kscale));
    let i = Dd::from_f64(positive(ni / scale)?);
    let o = Dd::from_f64(positive(no / scale)?);
    let e = Dd::from_f64(positive(ne / scale)?);
    let (p, p2, cosine, incident_normal) = match data.incidence {
        Incidence::Cosine(c) => {
            positive(c)?;
            if c > 1.0 {
                return Err("crystal incidence cosine must lie in (0,1]");
            }
            let cdd = Dd::from_f64(c);
            let sine2 = (Dd::ONE - cdd) * (Dd::ONE + cdd);
            let p2 = i * i * sine2;
            (i * Dd::from_f64(det::sqrt(sine2.to_f64())), p2, c, i * cdd)
        }
        Incidence::TangentialIndex(value) => {
            if !value.is_finite() || value.abs() >= ni {
                return Err("crystal frequency jet requires a propagating incident port");
            }
            let p = Dd::from_f64(value / scale);
            if value != 0.0 && p == Dd::ZERO {
                return Err("crystal tangential coordinate is unrepresentable");
            }
            let p2 = p * p;
            let incident_normal = (i * i - p2).sqrt();
            let cosine = positive((incident_normal / i).to_f64())?;
            (p, p2, cosine, incident_normal)
        }
    };
    let [ri, ro, re] = data.relative_permittivity_derivatives;
    // The exactly isotropic carrier uses the original scalar root. When both
    // principal local jets also coincide, its two spectral branches agree
    // through this order. Reuse that scalar identity instead of inventing a
    // grazing discriminant from normalized anisotropic coefficient residuals.
    // Equal phase indices with distinct local jets remain separate branches.
    let mode = if no == ne
        && ro == re
        && data.principal_wavenumber_derivatives[0] == data.principal_wavenumber_derivatives[1]
    {
        UniaxialMode::Ordinary
    } else {
        data.mode
    };
    let axis = if data.axis[1] == 0.0 {
        let [x, z] = unit_axis([data.axis[0], data.axis[2]])?;
        [x, 0.0, z]
    } else {
        // Match the vector carrier's original robust normalization.
        let axis_scale = positive(data.axis.iter().fold(0.0_f64, |m, x| m.max(x.abs())))?;
        let v = data.axis.map(|x| x / axis_scale);
        let length = positive(det::hypot(det::hypot(v[0], v[1]), v[2]))?;
        v.map(|x| x / length)
    };
    let [ax, ay, az] = axis.map(Dd::from_f64);
    let axis2 = ax * ax + ay * ay + az * az;
    positive(axis2.to_f64())?;
    let [xx, yy, zz, xz] = [ax * ax, ay * ay, az * az, ax * az].map(|v| v / axis2);
    let fixed_angle = data.constraint == IncidenceFrequencyConstraint::FixedAngle;

    // Preserve the original bulk jets in exactly axial principal limits.
    // These identities also avoid rebuilding the zero GVD of linear n(lambda).
    if p == Dd::ZERO {
        let principal = match mode {
            UniaxialMode::Ordinary => Some(0),
            UniaxialMode::Extraordinary if xx == Dd::ZERO && yy == Dd::ZERO => Some(0),
            UniaxialMode::Extraordinary if zz == Dd::ZERO => Some(1),
            _ => None,
        };
        if let Some(principal) = principal {
            return bulk_multiple(
                Z::real(Dd::ONE),
                data.principal_wavenumber_derivatives[principal],
                omega,
            );
        }
    }

    if fixed_angle {
        let matched = match mode {
            UniaxialMode::Ordinary if ni == no && ri == ro => Some(0),
            UniaxialMode::Extraordinary
                if xx == Dd::ZERO && yy == Dd::ZERO && ni == ne && ri == re =>
            {
                Some(0)
            }
            UniaxialMode::Extraordinary
                if yy == Dd::ZERO && zz == Dd::ZERO && ni == no && ri == ro =>
            {
                Some(1)
            }
            UniaxialMode::Extraordinary
                if xx == Dd::ZERO && zz == Dd::ZERO && ni == ne && ri == re =>
            {
                Some(1)
            }
            _ => None,
        };
        if let Some(principal) = matched {
            // The matched cutoff is exact in the original local jets. Keep
            // the authored cosine even when its square rounds out of 1-c^2.
            // Equal phase indices with DISTINCT slopes never enter here.
            return bulk_multiple(
                Z::real(Dd::from_f64(cosine)),
                data.principal_wavenumber_derivatives[principal],
                omega,
            );
        }
        let common = match mode {
            UniaxialMode::Ordinary => ri == ro,
            UniaxialMode::Extraordinary => ri == ro && ro == re,
        };
        if common {
            // The entire local index ellipsoid scales with the incident
            // index. Its original normal carrier is a constant multiple of Ko.
            return bulk_multiple(
                q.scale(Dd::ONE / o),
                data.principal_wavenumber_derivatives[0],
                omega,
            );
        }
    }

    if !fixed_angle {
        if mode == UniaxialMode::Ordinary {
            // At fixed kx, use the ORIGINAL bulk K jets directly:
            // kz'=K*K'/kz and kz''=K*K''/kz-K'^2*kx^2/kz^3.
            // This preserves exactly stationary/anomalous bulk slopes and
            // avoids rebuilding them by subtracting index derivative terms.
            let conversion = Dd::from_f64(scale) / Dd::from_f64(C0);
            let [first, second_scaled] =
                data.principal_wavenumber_derivatives[0].map(|v| Dd::from_f64(v) / conversion);
            let first_over_q = Z::real(first).divide(q)?;
            let first_normal = first_over_q.scale(o * conversion).finish()?;
            let second_normal = Z::real(o * second_scaled)
                .sub(first_over_q.mul(first_over_q).scale(p2))
                .divide(q)?
                .scale(conversion / Dd::from_f64(omega))
                .finish()?;
            return Ok([first_normal, second_normal]);
        }
        if data.principal_wavenumber_derivatives == [[0.0; 2]; 2] {
            // Both original principal K are locally stationary through second
            // order, so the whole fixed-kx crystal carrier is stationary too.
            return Ok([C64::ZERO; 2]);
        }
    }

    let o2 = index_square(o, ro);
    let e2 = index_square(e, re);
    let (p1, p_second, p2_derivatives) = if fixed_angle {
        let r1 = Dd::from_f64(ri[0]);
        let r2 = Dd::from_f64(ri[1]);
        (
            p * r1 * Dd::from_f64(0.5),
            p * (r2 * Dd::from_f64(0.5) - r1 * r1 * Dd::from_f64(0.25)),
            [p2 * r1, p2 * r2],
        )
    } else {
        // kx is fixed, while p=c0*kx/(omega*scale) varies as 1/omega.
        (
            -p,
            p * Dd::from_f64(2.0),
            [-p2 * Dd::from_f64(2.0), p2 * Dd::from_f64(6.0)],
        )
    };
    let two = Dd::from_f64(2.0);
    let (q1, q2) = match mode {
        UniaxialMode::Ordinary => {
            let derivatives = root_jet(q, o2, p2_derivatives, ni == no, i * i, ri, fixed_angle)?;
            (derivatives[0], derivatives[1])
        }
        UniaxialMode::Extraordinary => {
            // Retain the factored outgoing root instead of reconstructing
            // F_q=2(f*q+g*p) from the rounded normal phase q. At tilted
            // grazing incidence that subtraction can erase normal energy.
            let f = combine(o2, xx + yy, e2, zz);
            let g = combine(e2, xz, o2, -xz);
            let (cutoff, factor, root, matched) = if axis[1] == 0.0 {
                let matched = f[0] == i * i;
                let root = if matched {
                    Z::real(incident_normal)
                } else {
                    outgoing(f[0] - p2, data.kind)?
                };
                let r_o = ro.map(Dd::from_f64);
                let r_e = re.map(Dd::from_f64);
                let first = (r_o[0] + r_e[0]) * Dd::from_f64(0.5);
                let difference = r_o[0] - r_e[0];
                let second = (r_o[1] + r_e[1]) * Dd::from_f64(0.5)
                    - difference * difference * Dd::from_f64(0.25);
                (f, [o * e, o * e * first, o * e * second], root, matched)
            } else {
                let determinant = product(o2, combine(o2, yy, e2, xx + zz));
                let numerator = product(f, product(o2, e2));
                let matched = numerator[0] == determinant[0] * i * i;
                let cutoff = quotient(numerator, determinant, matched.then_some(i * i));
                let factor = real_root_jet(determinant);
                let root = if matched {
                    Z::real(incident_normal)
                } else {
                    outgoing(numerator[0] - determinant[0] * p2, data.kind)?
                        .scale(Dd::ONE / factor[0])
                };
                (cutoff, factor, root, matched)
            };
            let r = root_jet(
                root,
                cutoff,
                p2_derivatives,
                matched,
                i * i,
                ri,
                fixed_angle,
            )?;
            // q=(-g*p+factor*root)/f. The unchanged native q enters only
            // quotient terms divided by f, never the vanishing root gradient.
            let numerator_first = Z::real(-(g[1] * p + g[0] * p1))
                .add(root.scale(factor[1]))
                .add(r[0].scale(factor[0]));
            let first = numerator_first.sub(q.scale(f[1])).scale(Dd::ONE / f[0]);
            let numerator_second = Z::real(-(g[2] * p + two * g[1] * p1 + g[0] * p_second))
                .add(root.scale(factor[2]))
                .add(r[0].scale(two * factor[1]))
                .add(r[1].scale(factor[0]));
            let second = numerator_second
                .sub(q.scale(f[2]))
                .sub(first.scale(two * f[1]))
                .scale(Dd::ONE / f[0]);
            (first, second)
        }
    };
    let conversion = Dd::from_f64(scale) / Dd::from_f64(C0);
    let first = q.add(q1).scale(conversion).finish()?;
    let second = q1
        .scale(two)
        .add(q2)
        .scale(conversion / Dd::from_f64(omega))
        .finish()?;
    Ok([first, second])
}

// Jets of r^2=cutoff-p^2. At an exactly matched cutoff, keep r=i*cos
// explicitly, and separate the incident-root terms before dividing by r.
// This retains both derivatives when cos^2 is below subtraction precision.
fn root_jet(
    root: Z,
    cutoff: [Dd; 3],
    p2_derivatives: [Dd; 2],
    matched: bool,
    i2: Dd,
    incident_relative: [f64; 2],
    fixed_angle: bool,
) -> Result<[Z; 2]> {
    let half = Dd::from_f64(0.5);
    let (first, leading_second) = if matched {
        let relative = if fixed_angle {
            incident_relative
        } else {
            [-2.0, 6.0]
        };
        let [a, b] = relative.map(Dd::from_f64);
        (
            Z::real((cutoff[1] - i2 * a) * half)
                .divide(root)?
                .add(root.scale(a * half)),
            Z::real((cutoff[2] - i2 * b) * half)
                .divide(root)?
                .add(root.scale(b * half)),
        )
    } else {
        (
            Z::real((cutoff[1] - p2_derivatives[0]) * half).divide(root)?,
            Z::real((cutoff[2] - p2_derivatives[1]) * half).divide(root)?,
        )
    };
    let second = leading_second.sub(first.mul(first.divide(root)?));
    Ok([first, second])
}

fn outgoing(value: Dd, kind: UniaxialWaveKind) -> Result<Z> {
    if ![value.hi, value.lo].into_iter().all(f64::is_finite) || value == Dd::ZERO {
        return Err("crystal normal-frequency discriminant is unresolved");
    }
    let root = value.abs().sqrt();
    match kind {
        UniaxialWaveKind::Propagating if Dd::ZERO.lt(value) => Ok(Z::real(root)),
        UniaxialWaveKind::Evanescent if value.lt(Dd::ZERO) => Ok(Z {
            re: Dd::ZERO,
            im: -root,
        }),
        _ => Err("crystal normal-frequency branch disagrees with the native carrier"),
    }
}

fn quotient(numerator: [Dd; 3], denominator: [Dd; 3], exact_value: Option<Dd>) -> [Dd; 3] {
    let value = exact_value.unwrap_or_else(|| numerator[0] / denominator[0]);
    let first = (numerator[1] - value * denominator[1]) / denominator[0];
    let second =
        (numerator[2] - value * denominator[2] - first * denominator[1] * Dd::from_f64(2.0))
            / denominator[0];
    [value, first, second]
}

fn real_root_jet(value: [Dd; 3]) -> [Dd; 3] {
    let root = value[0].sqrt();
    let first = (value[1] / root) * Dd::from_f64(0.5);
    let second = (value[2] / root) * Dd::from_f64(0.5) - first * (first / root);
    [root, first, second]
}

fn index_square(index: Dd, relative: [f64; 2]) -> [Dd; 3] {
    let square = index * index;
    [
        square,
        square * Dd::from_f64(relative[0]),
        square * Dd::from_f64(relative[1]),
    ]
}

fn combine(a: [Dd; 3], wa: Dd, b: [Dd; 3], wb: Dd) -> [Dd; 3] {
    core::array::from_fn(|j| a[j] * wa + b[j] * wb)
}

fn product(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    [
        a[0] * b[0],
        a[1] * b[0] + a[0] * b[1],
        a[2] * b[0] + a[1] * b[1] * Dd::from_f64(2.0) + a[0] * b[2],
    ]
}

fn bulk_multiple(factor: Z, [first, second_scaled]: [f64; 2], omega: f64) -> Result<[C64; 2]> {
    Ok([
        factor.scale(Dd::from_f64(first)).finish()?,
        factor
            .scale(Dd::from_f64(second_scaled) / Dd::from_f64(omega))
            .finish()?,
    ])
}

fn positive(value: f64) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err("positive representable crystal frequency data required")
    }
}

fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

// Only fixed-size complex DD contractions are needed. Smith's ratio division
// avoids squaring a weak normal gradient or overflowing a large one.
#[derive(Clone, Copy)]
struct Z {
    re: Dd,
    im: Dd,
}
impl Z {
    fn real(re: Dd) -> Self {
        Self { re, im: Dd::ZERO }
    }
    fn from_c64(value: C64) -> Self {
        Self {
            re: Dd::from_f64(value.re),
            im: Dd::from_f64(value.im),
        }
    }
    fn add(self, other: Self) -> Self {
        Self {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }
    fn sub(self, other: Self) -> Self {
        Self {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }
    fn mul(self, other: Self) -> Self {
        Self {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }
    fn scale(self, value: Dd) -> Self {
        Self {
            re: self.re * value,
            im: self.im * value,
        }
    }
    fn divide(self, denominator: Self) -> Result<Self> {
        if denominator.re == Dd::ZERO && denominator.im == Dd::ZERO {
            return Err("crystal normal-frequency gradient is unresolved");
        }
        let (re, im) = if denominator.im.abs().lt(denominator.re.abs()) {
            let ratio = denominator.im / denominator.re;
            let divisor = denominator.re + denominator.im * ratio;
            (
                (self.re + self.im * ratio) / divisor,
                (self.im - self.re * ratio) / divisor,
            )
        } else {
            let ratio = denominator.re / denominator.im;
            let divisor = denominator.im + denominator.re * ratio;
            (
                (self.re * ratio + self.im) / divisor,
                (self.im * ratio - self.re) / divisor,
            )
        };
        if ![re.hi, re.lo, im.hi, im.lo].into_iter().all(f64::is_finite) {
            return Err("crystal normal-frequency derivatives are unrepresentable");
        }
        Ok(Self { re, im })
    }
    fn finish(self) -> Result<C64> {
        let value = C64::new(self.re.to_f64(), self.im.to_f64());
        if !finite(value)
            || (value.re == 0.0 && self.re != Dd::ZERO)
            || (value.im == 0.0 && self.im != Dd::ZERO)
        {
            return Err("crystal normal-frequency derivatives are unrepresentable");
        }
        Ok(value)
    }
}

#[cfg(test)]
mod tests;
