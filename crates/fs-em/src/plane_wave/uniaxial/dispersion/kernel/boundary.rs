//! Fresnel frequency jets from the original principal-plane power ports.
//!
//! A noncritical transparent port has a real logarithmic admittance slope,
//! including TIR: the outgoing gradient root keeps a fixed real/imaginary
//! phase on its admitted local branch. The native complex fields retain that
//! phase and every permeability/normalization factor. No new root is solved.

use super::{Z, positive};
use crate::plane_wave::dispersion::IncidenceFrequencyConstraint;
use crate::plane_wave::uniaxial::kernel::{
    Incidence, UniaxialInterfaceResponse, UniaxialMode, UniaxialWaveKind, binary_scale, unit_axis,
};
use fs_math::{c64::C64, dd::Dd};

type Result<T> = core::result::Result<T, &'static str>;

pub(in crate::plane_wave::uniaxial) struct Input {
    pub indices: [f64; 3],
    /// Original [incident, ordinary, extraordinary] scaled epsilon derivatives
    /// divided by each original epsilon, through ordinary second order.
    pub relative_permittivity_derivatives: [[f64; 2]; 3],
    pub incidence: Incidence,
    pub axis: [f64; 2],
    pub constraint: IncidenceFrequencyConstraint,
    pub carrier: UniaxialInterfaceResponse,
}

/// [omega*r', omega²*r'']; tE has the same jets and normalized tH the opposite.
pub(in crate::plane_wave::uniaxial) fn amplitude_derivatives(data: Input) -> Result<[C64; 2]> {
    if data.carrier.wave_kind == UniaxialWaveKind::Critical {
        return Err("critical crystal boundaries have no admitted frequency jet");
    }
    for value in data.indices {
        positive(value)?;
    }
    if data
        .relative_permittivity_derivatives
        .iter()
        .flatten()
        .any(|value| !value.is_finite())
    {
        return Err("finite original crystal boundary frequency jets required");
    }
    let [ni, no, ne] = data.indices;
    let scale = binary_scale(ni.max(no).max(ne));
    let [i, o, e] = data.indices.map(|value| Dd::from_f64(value / scale));
    let [i2, o2, e2] = [i * i, o * o, e * e];
    positive(i2.to_f64())?;
    positive(o2.to_f64())?;
    positive(e2.to_f64())?;
    let (p2, incident_square) = match data.incidence {
        Incidence::Cosine(cosine) => {
            positive(cosine)?;
            if cosine > 1.0 {
                return Err("crystal boundary incidence cosine must lie in (0,1]");
            }
            let c = Dd::from_f64(cosine);
            let normal = i * c;
            let sine_squared = (Dd::ONE - c) * (Dd::ONE + c);
            (i2 * sine_squared, normal * normal)
        }
        Incidence::TangentialIndex(value) => {
            if !value.is_finite() || value.abs() >= ni {
                return Err("crystal boundary jet requires a propagating incident port");
            }
            let p = Dd::from_f64(value / scale);
            if value != 0.0 && p == Dd::ZERO {
                return Err("crystal boundary tangential coordinate is unrepresentable");
            }
            let p2 = p * p;
            (p2, i2 - p2)
        }
    };
    let [ax, az] = unit_axis(data.axis)?.map(Dd::from_f64);
    let norm = ax * ax + az * az;
    let xx = ax * ax / norm;
    let zz = az * az / norm;
    let [ri, ro, re] = data
        .relative_permittivity_derivatives
        .map(|pair| pair.map(Dd::from_f64));
    // Match the original native exactly-isotropic coefficient identity.
    let cutoff = match data.carrier.mode {
        UniaxialMode::Ordinary => o2,
        UniaxialMode::Extraordinary if no == ne => o2,
        UniaxialMode::Extraordinary => o2 * xx + e2 * zz,
    };
    let difference = cutoff - i2;
    // Factor each principal derivative against the incident derivative.
    // Matched local cutoffs retain exact zero even at authored tiny cosines.
    let difference_derivatives = core::array::from_fn(|j| {
        let contrast = match data.carrier.mode {
            UniaxialMode::Ordinary => o2 * (ro[j] - ri[j]),
            UniaxialMode::Extraordinary => o2 * xx * (ro[j] - ri[j]) + e2 * zz * (re[j] - ri[j]),
        };
        difference * ri[j] + contrast
    });
    let [root_first, root_linear_second, incident_first] =
        if difference == Dd::ZERO && difference_derivatives == [Dd::ZERO; 2] {
            [Dd::ZERO; 3]
        } else {
            let incident_relative = match data.constraint {
                IncidenceFrequencyConstraint::FixedAngle => ri,
                IncidenceFrequencyConstraint::FixedTransverseWavenumber => {
                    positive(incident_square.to_f64())?;
                    [
                        (i2 * ri[0] + p2 * Dd::from_f64(2.0)) / incident_square,
                        (i2 * ri[1] - p2 * Dd::from_f64(6.0)) / incident_square,
                    ]
                }
            };
            let numerators = core::array::from_fn::<_, 2, _>(|j| {
                difference_derivatives[j] - difference * incident_relative[j]
            });
            if numerators == [Dd::ZERO; 2] {
                [Dd::ZERO; 3]
            } else {
                // Preserve the same authored cutoff identity as the native carrier.
                let transmitted_square = if difference == Dd::ZERO {
                    incident_square
                } else {
                    cutoff - p2
                };
                if transmitted_square == Dd::ZERO {
                    return Err("crystal boundary gradient root is unrepresentable");
                }
                [
                    numerators[0] / transmitted_square,
                    numerators[1] / transmitted_square,
                    incident_relative[0],
                ]
            }
        };
    // L=log(Yi/Yt). TM uses the outgoing gradient root, never signed kz.
    // Permeability is the original fixed scalar snapshot, so its log jet is zero.
    let half = Dd::from_f64(0.5);
    let (sign, material_log) = match data.carrier.mode {
        UniaxialMode::Ordinary => (-half, [Dd::ZERO; 2]),
        UniaxialMode::Extraordinary if ri == ro && ro == re => (half, [Dd::ZERO; 2]),
        UniaxialMode::Extraordinary => (
            half,
            [
                ri[0] - (ro[0] + re[0]) * half,
                ri[1] - ri[0] * ri[0] - (ro[1] - ro[0] * ro[0] + re[1] - re[0] * re[0]) * half,
            ],
        ),
    };
    if material_log == [Dd::ZERO; 2] && root_first == Dd::ZERO && root_linear_second == Dd::ZERO {
        return Ok([C64::ZERO; 2]);
    }
    // 2ab/(a+b)^2 = tE*tH/2 uses the two independently formed native field
    // numerators. It never divides by r, including an exact Brewster zero.
    let reflected = Z::from_c64(data.carrier.reflected_tangential_e);
    let factor = Z::from_c64(data.carrier.transmitted_tangential_e)
        .mul(Z::from_c64(data.carrier.transmitted_conjugate_h))
        .scale(half);
    let logarithmic_first = material_log[0] + root_first * sign;
    let first = factor.scale(logarithmic_first);
    // Multiply the physical prefactor BEFORE the large grazing products.
    // Unmatched fixed-kx grazing has L''~cos^-4 but r''~cos^-3: forming the
    // naked log curvature or L'^2 would overflow while r'' is still finite.
    let second = factor
        .scale(material_log[1])
        .add(factor.scale(sign * root_linear_second))
        .sub(
            factor
                .scale(sign * root_first)
                .scale(root_first + incident_first * Dd::from_f64(2.0)),
        )
        .sub(reflected.mul(first).scale(logarithmic_first));
    Ok([first.finish()?, second.finish()?])
}

#[cfg(test)]
mod tests;
