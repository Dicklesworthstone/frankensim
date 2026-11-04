//! Exact uniform-medium ellipsoid resistance, prepared once from positive
//! geometric integrals. Anand & Narsimhan, JFM 983 (2024), Appendix B,
//! equations B1-B7, DOI 10.1017/jfm.2024.136. Their integrals include a*b*c;
//! the common factor cancels in the unscaled integrals used here.
//!
//! With Delta(t)=sqrt(product(a_i^2+t)), chi=int(dt/Delta) and
//! alpha_i=int(dt/((a_i^2+t)*Delta)), translation is
//! 16*pi*mu/(chi+a_i^2*alpha_i), and rotation about i is
//! (16*pi*mu/3)*(a_j^2+a_k^2)/(a_j^2*alpha_j+a_k^2*alpha_k).
use super::{Error, finite, poll, positive};

/// Cold-path quadrature work and refinement controls. The observed change
/// between successive meshes is an error estimate, not a certified enclosure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipsoidQuadrature {
    /// Initial number of equal panels in the compact integration coordinate.
    pub initial_panels: usize,
    /// Maximum number of panels, at most 4096 and at least twice the initial.
    pub maximum_panels: usize,
    /// Relative change allowed in each of the six physical coefficients,
    /// within [64*epsilon,1e-3]. No coefficient is accepted without refinement.
    pub relative_tolerance: f64,
}

/// Resistance in the boundary's authored ellipsoid axes, never sorted into
/// inertia order. All six eigenvalues are finite and strictly positive.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipsoidResistance {
    /// Translational resistance along each semiaxis [kg/s].
    pub translation_kg_s: [f64; 3],
    /// Rotational resistance about each semiaxis [kg m2/s].
    pub rotation_kg_m2_s: [f64; 3],
    /// Accepted panel count; zero denotes the exact sphere limit.
    pub quadrature_panels: usize,
    /// Largest relative change in a coefficient on the final refinement.
    pub estimated_relative_change: f64,
}

pub(super) fn prepare(
    axes: [f64; 3],
    viscosity: f64,
    controls: EllipsoidQuadrature,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<EllipsoidResistance, Error> {
    poll(cancelled)?;
    for a in axes {
        positive(a)?;
    }
    positive(viscosity)?;
    if controls.initial_panels == 0
        || controls.maximum_panels > 4096
        || controls.initial_panels > controls.maximum_panels / 2
        || !controls.relative_tolerance.is_finite()
        || !(64.0 * f64::EPSILON..=1e-3).contains(&controls.relative_tolerance)
    {
        return Err(Error::Invalid("invalid ellipsoid quadrature controls"));
    }
    let length = axes[0].min(axes[1]).min(axes[2]);
    if axes[0] == axes[1] && axes[1] == axes[2] {
        let translation = positive(6.0 * core::f64::consts::PI * viscosity * length)?;
        let rotation = positive((4.0 / 3.0) * translation * length * length)?;
        poll(cancelled)?;
        return Ok(EllipsoidResistance {
            translation_kg_s: [translation; 3],
            rotation_kg_m2_s: [rotation; 3],
            quadrature_panels: 0,
            estimated_relative_change: 0.0,
        });
    }
    let mut ratios_squared = [0.0; 3];
    for i in 0..3 {
        let ratio = positive(axes[i] / length)?;
        ratios_squared[i] = positive(ratio * ratio)?;
    }
    let translation_scale = positive(16.0 * core::f64::consts::PI * viscosity * length)?;
    let rotation_scale = positive((translation_scale / 3.0) * length * length)?;
    let coefficients = |integrals: [f64; 4]| -> Result<EllipsoidResistance, Error> {
        let mut translation = [0.0; 3];
        let mut rotation = [0.0; 3];
        for i in 0..3 {
            translation[i] = positive(
                translation_scale / positive(integrals[0] + ratios_squared[i] * integrals[i + 1])?,
            )?;
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            rotation[i] = positive(
                rotation_scale * positive(ratios_squared[j] + ratios_squared[k])?
                    / positive(
                        ratios_squared[j] * integrals[j + 1] + ratios_squared[k] * integrals[k + 1],
                    )?,
            )?;
        }
        Ok(EllipsoidResistance {
            translation_kg_s: translation,
            rotation_kg_m2_s: rotation,
            quadrature_panels: 0,
            estimated_relative_change: 0.0,
        })
    };
    let mut panels = controls.initial_panels;
    let mut previous = coefficients(integrate(ratios_squared, panels, cancelled)?)?;
    while panels <= controls.maximum_panels / 2 {
        panels *= 2;
        let mut current = coefficients(integrate(ratios_squared, panels, cancelled)?)?;
        let mut change = 0.0_f64;
        for i in 0..3 {
            change = change
                .max(finite(
                    (current.translation_kg_s[i] / previous.translation_kg_s[i] - 1.0).abs(),
                )?)
                .max(finite(
                    (current.rotation_kg_m2_s[i] / previous.rotation_kg_m2_s[i] - 1.0).abs(),
                )?);
        }
        current.quadrature_panels = panels;
        current.estimated_relative_change = change;
        if change <= controls.relative_tolerance {
            poll(cancelled)?;
            return Ok(current);
        }
        previous = current;
    }
    Err(Error::Invalid(
        "ellipsoid quadrature did not converge within the declared panel budget",
    ))
}

// t=length^2*(u^-2-1) maps the complete semi-infinite integral onto [0,1].
// Every quadrature contribution is positive. No tail cutoff or aspect-ratio
// fit is introduced. Compensated accumulation avoids summation-order loss.
fn integrate(
    ratios_squared: [f64; 3],
    panels: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 4], Error> {
    const NODES: [f64; 4] = [
        0.1834346424956498,
        0.525532409916329,
        0.7966664774136267,
        0.9602898564975363,
    ];
    const WEIGHTS: [f64; 4] = [
        0.362683783378362,
        0.3137066458778873,
        0.2223810344533745,
        0.1012285362903763,
    ];
    let width = 1.0 / panels as f64;
    let mut sums = [0.0; 4];
    let mut corrections = [0.0; 4];
    for panel in 0..panels {
        poll(cancelled)?;
        let center = (panel as f64 + 0.5) * width;
        for n in 0..4 {
            for sign in [-1.0, 1.0] {
                let u = center + sign * 0.5 * width * NODES[n];
                let u2 = u * u;
                let factors = ratios_squared.map(|r| 1.0 + (r - 1.0) * u2);
                let denominator = positive(
                    fs_math::det::sqrt(factors[0])
                        * fs_math::det::sqrt(factors[1])
                        * fs_math::det::sqrt(factors[2]),
                )?;
                let base = positive(width * WEIGHTS[n] / denominator)?;
                let terms = [
                    base,
                    positive(base * u2 / factors[0])?,
                    positive(base * u2 / factors[1])?,
                    positive(base * u2 / factors[2])?,
                ];
                for i in 0..4 {
                    let increment = terms[i] - corrections[i];
                    let next = finite(sums[i] + increment)?;
                    corrections[i] = (next - sums[i]) - increment;
                    sums[i] = next;
                }
            }
        }
    }
    for value in sums {
        positive(value)?;
    }
    Ok(sums)
}
