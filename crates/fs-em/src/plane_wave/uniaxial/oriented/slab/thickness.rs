//! Differentiate the original two-face Maxwell system at fixed frequency.
//!
//! If A(d)c(d)=b, then A c'=-A'c and A c''=-A''c-2A'c'. Only the
//! entrance-to-exit propagation factors depend on d. All three right-hand
//! sides use one native pivoted factorization; no neighboring slab or material
//! response is evaluated. Derivatives first use a dimensionless thickness
//! coordinate, keeping wave-number scales out of the boundary solve.

use super::kernel::{BoundaryFactorization, SlabWave, UniaxialSlabResponse, tangential};
use crate::plane_wave::uniaxial::dispersion::kernel::Z;
use fs_math::{c64::C64, dd::Dd};

type Vector = [C64; 3];
type Matrix = [[C64; 8]; 8];
type Result<T> = core::result::Result<T, &'static str>;

/// Ordinary (not Taylor-coefficient) thickness derivatives of a Cartesian
/// complex peak field. E has units V/m and H has units A/m before
/// differentiation. The field's native entrance/exit anchor is retained.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabFieldThicknessJet {
    /// dE/dd [V/m²].
    pub electric_first_v_m2: Vector,
    /// d²E/dd² [V/m³].
    pub electric_second_v_m3: Vector,
    /// dH/dd [A/m²].
    pub magnetic_first_a_m2: Vector,
    /// d²H/dd² [A/m³].
    pub magnetic_second_a_m3: Vector,
}

/// Native coherent slab response and its analytic thickness sensitivity.
/// Frequency, source fields, materials, optic axis and incidence are fixed.
/// The entrance remains z=0 and the exit moves with d. Consequently the
/// transmitted and backward fields are derivatives at that moving face;
/// an observation at fixed global z needs its separate propagation derivative.
/// No scalar pulse delay is assigned to a coherent vector field.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniaxialSlabThicknessDerivatives {
    /// Unchanged original carrier, including the original power/field gates.
    pub response: UniaxialSlabResponse,
    /// Reflected field anchored at the fixed entrance.
    pub reflected: SlabFieldThicknessJet,
    /// Transmitted field anchored at the moving exit.
    pub transmitted: SlabFieldThicknessJet,
    /// Ordinary then extraordinary, anchored at the fixed entrance.
    pub forward: [SlabFieldThicknessJet; 2],
    /// Ordinary then extraordinary, anchored at the moving exit.
    pub backward: [SlabFieldThicknessJet; 2],
    /// dR/dd [1/m], from the signed normal E/H flux.
    pub reflected_power_first_per_m: f64,
    /// d²R/dd² [1/m²], including the first-field cross product.
    pub reflected_power_second_per_m2: f64,
    /// dT/dd [1/m], from the original isotropic exit E/H flux.
    pub transmitted_power_first_per_m: f64,
    /// d²T/dd² [1/m²]. Critical/evanescent exit power is identically zero.
    pub transmitted_power_second_per_m2: f64,
}

pub(super) struct BoundaryWitness {
    pub factorization: BoundaryFactorization<8>,
    pub coefficients: [C64; 8],
    pub reflected: [SlabWave; 2],
    pub forward: [SlabWave; 2],
    pub backward: [SlabWave; 2],
    pub exit: [SlabWave; 2],
    pub phase_forward: [C64; 2],
    pub phase_backward: [C64; 2],
    pub amplitude: f64,
    pub incoming: f64,
}

pub(super) fn differentiate(
    response: UniaxialSlabResponse,
    w: BoundaryWitness,
) -> Result<UniaxialSlabThicknessDerivatives> {
    let wave_scale = w
        .forward
        .iter()
        .chain(&w.backward)
        .fold(0.0_f64, |scale, wave| {
            scale.max(wave.normal_wavenumber_rad_m.abs())
        });
    if !wave_scale.is_finite() || wave_scale <= 0.0 {
        return Err("finite positive slab thickness wave-number scale required");
    }
    let mut matrix_first = [[C64::ZERO; 8]; 8];
    let mut matrix_second = [[C64::ZERO; 8]; 8];
    for mode in 0..2 {
        let forward = phase_jet(
            w.forward[mode].normal_wavenumber_rad_m,
            w.phase_forward[mode],
            wave_scale,
        )?;
        let backward = phase_jet(
            -w.backward[mode].normal_wavenumber_rad_m,
            w.phase_backward[mode],
            wave_scale,
        )?;
        let a = tangential(w.forward[mode].electric_v_m, w.forward[mode].magnetic_a_m);
        let b = tangential(w.backward[mode].electric_v_m, w.backward[mode].magnetic_a_m);
        for row in 0..4 {
            matrix_first[row][mode + 4] = -b[row] * backward[0];
            matrix_second[row][mode + 4] = -b[row] * backward[1];
            matrix_first[row + 4][mode + 2] = a[row] * forward[0];
            matrix_second[row + 4][mode + 2] = a[row] * forward[1];
        }
    }
    let first = w
        .factorization
        .solve(rhs(&matrix_first, &w.coefficients, None)?)?;
    let second = w.factorization.solve(rhs(
        &matrix_second,
        &w.coefficients,
        Some((&matrix_first, &first)),
    )?)?;
    let coefficients = [w.coefficients, first, second];
    let reflected = fields(w.reflected, coefficients.map(|c| [c[0], c[1]]));
    let transmitted = fields(w.exit, coefficients.map(|c| [c[6], c[7]]));
    let rp = power_derivatives(&reflected, -w.incoming, wave_scale)?;
    let tp = if response.exit_wave_kind == super::super::super::UniaxialWaveKind::Propagating {
        power_derivatives(&transmitted, w.incoming, wave_scale)?
    } else {
        // One original isotropic critical/evanescent exit carries no normal
        // flux for ANY thickness. Avoid differentiating cancellation noise.
        [0.0; 2]
    };
    let internal = |waves: [SlabWave; 2], start: usize| -> Result<[SlabFieldThicknessJet; 2]> {
        let one = |mode| {
            let zero = SlabWave {
                electric_v_m: [C64::ZERO; 3],
                magnetic_a_m: [C64::ZERO; 3],
                ..waves[mode]
            };
            let fields = fields(
                [waves[mode], zero],
                coefficients.map(|c| [c[start + mode], C64::ZERO]),
            );
            dimensional_jet(fields, w.amplitude, wave_scale)
        };
        Ok([one(0)?, one(1)?])
    };
    Ok(UniaxialSlabThicknessDerivatives {
        response,
        reflected: dimensional_jet(reflected, w.amplitude, wave_scale)?,
        transmitted: dimensional_jet(transmitted, w.amplitude, wave_scale)?,
        forward: internal(w.forward, 2)?,
        backward: internal(w.backward, 4)?,
        reflected_power_first_per_m: rp[0],
        reflected_power_second_per_m2: rp[1],
        transmitted_power_first_per_m: tp[0],
        transmitted_power_second_per_m2: tp[1],
    })
}

fn phase_jet(k: C64, phase: C64, scale: f64) -> Result<[C64; 2]> {
    if phase == C64::ZERO {
        // A rounded zero transfer can still have a representable dimensional
        // derivative after multiplication by a large wave number. The carrier
        // remains valid, but do not publish fabricated zero sensitivity.
        return Err("underflowed crystal propagation has no admitted thickness jet");
    }
    let generator = C64::new(k.im / scale, -k.re / scale);
    if k.re != 0.0 && generator.im == 0.0 || k.im != 0.0 && generator.re == 0.0 {
        return Err("crystal thickness phase scale is unrepresentable");
    }
    let first = Z::from_c64(phase).mul(Z::from_c64(generator));
    let second = first.mul(Z::from_c64(generator));
    Ok([first.finish()?, second.finish()?])
}

fn rhs(
    matrix: &Matrix,
    carrier: &[C64; 8],
    cross: Option<(&Matrix, &[C64; 8])>,
) -> Result<[C64; 8]> {
    let mut result = [C64::ZERO; 8];
    for row in 0..8 {
        let mut value = Z::real(Dd::ZERO);
        for column in 0..8 {
            value = value.sub(Z::from_c64(matrix[row][column]).mul(Z::from_c64(carrier[column])));
            if let Some((matrix, first)) = cross {
                value = value.sub(
                    Z::from_c64(matrix[row][column])
                        .mul(Z::from_c64(first[column]))
                        .scale(Dd::from_f64(2.0)),
                );
            }
        }
        result[row] = value.finish()?;
    }
    Ok(result)
}

/// Each order stores [Ex,Ey,Ez,Hx,Hy,Hz] before incident amplitude restoration.
fn fields(waves: [SlabWave; 2], coefficients: [[C64; 2]; 3]) -> [[Z; 6]; 3] {
    coefficients.map(|c| {
        core::array::from_fn(|j| {
            let value = |m: usize| {
                if j < 3 {
                    waves[m].electric_v_m[j]
                } else {
                    waves[m].magnetic_a_m[j - 3]
                }
            };
            Z::from_c64(value(0))
                .mul(Z::from_c64(c[0]))
                .add(Z::from_c64(value(1)).mul(Z::from_c64(c[1])))
        })
    })
}

fn dimensional_jet(
    fields: [[Z; 6]; 3],
    amplitude: f64,
    wave_scale: f64,
) -> Result<SlabFieldThicknessJet> {
    let mut result = [[C64::ZERO; 6]; 2];
    for order in 0..2 {
        for component in 0..6 {
            let mut value = fields[order + 1][component].scale(Dd::from_f64(amplitude));
            for _ in 0..=order {
                value = value.scale(Dd::from_f64(wave_scale));
            }
            result[order][component] = value.finish()?;
        }
    }
    Ok(SlabFieldThicknessJet {
        electric_first_v_m2: [result[0][0], result[0][1], result[0][2]],
        electric_second_v_m3: [result[1][0], result[1][1], result[1][2]],
        magnetic_first_a_m2: [result[0][3], result[0][4], result[0][5]],
        magnetic_second_a_m3: [result[1][3], result[1][4], result[1][5]],
    })
}

fn power_derivatives(fields: &[[Z; 6]; 3], incident: f64, wave_scale: f64) -> Result<[f64; 2]> {
    let dot = |a: Z, b: Z| a.re * b.re + a.im * b.im;
    let flux = |e: usize, h: usize| {
        (dot(fields[e][0], fields[h][4]) - dot(fields[e][1], fields[h][3])) * Dd::from_f64(0.5)
    };
    let raw = [
        flux(1, 0) + flux(0, 1),
        flux(2, 0) + flux(1, 1) * Dd::from_f64(2.0) + flux(0, 2),
    ];
    let mut result = [0.0; 2];
    for order in 0..2 {
        let mut value = raw[order] / Dd::from_f64(incident);
        for _ in 0..=order {
            value = value * Dd::from_f64(wave_scale);
        }
        result[order] = value.to_f64();
        if !result[order].is_finite() || result[order] == 0.0 && value != Dd::ZERO {
            return Err("slab thickness power derivative is unrepresentable");
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests;
