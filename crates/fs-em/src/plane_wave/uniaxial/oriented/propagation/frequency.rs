//! Coherent spatial frequency jets in the original phasor/decay accumulator.
//!
//! Common phase AND its two derivatives are factored before modal interference.
//! A common evanescent exposure is applied only after the weighted derivative
//! polynomials, retaining derivatives whose carrier has already underflowed.

use super::super::boundary::FieldDerivatives;
use super::{
    Carrier, Complex, FieldSum, Fields, OrientedCrystalFieldSide, Result, Vector, ZERO,
    attenuated_fields, attenuation_exposure, multiply, relative_rotation,
};
use crate::plane_wave::uniaxial::dispersion::kernel::Z;
use fs_math::{c64::C64, dd::Dd};

#[derive(Clone, Copy)]
pub(in crate::plane_wave::uniaxial::oriented) struct FrequencyWave<'a> {
    pub normal: C64,
    /// Original dimensional [kz', kz''] in [s/m, s²/m].
    pub normal_derivatives: [Z; 2],
    pub electric: &'a Vector,
    pub magnetic: &'a Vector,
    /// Original [omega F', omega² F''] at the boundary.
    pub derivatives: &'a FieldDerivatives,
}

pub(in crate::plane_wave::uniaxial::oriented) fn field_frequency_derivatives(
    omega: f64,
    tangential: f64,
    tangential_derivatives: [Dd; 2],
    side: OrientedCrystalFieldSide,
    point: [f64; 3],
    waves: [Option<FrequencyWave<'_>>; 2],
) -> Result<FieldDerivatives> {
    // The original accumulator owns coordinate and half-space admission.
    let mut sums = [
        FieldSum::new(side, point)?,
        FieldSum::new(side, point)?,
        FieldSum::new(side, point)?,
    ];
    let depth = sums[0].depth;
    let mut reference: Option<(FrequencyWave<'_>, f64)> = None;
    for wave in waves.into_iter().flatten().filter(active) {
        let exposure = attenuation_exposure(wave.normal.im, depth)?;
        // Compare the original roots, not rounded large exposures: two
        // exposures may round equal while their differential decay is finite.
        let less_decay = reference.is_none_or(|(previous, _)| {
            if point[2] >= 0.0 {
                wave.normal.im > previous.normal.im
            } else {
                wave.normal.im < previous.normal.im
            }
        });
        if less_decay {
            reference = Some((wave, exposure));
        }
    }
    let Some((reference_wave, common_exposure)) = reference else {
        return Ok(zero_derivatives());
    };
    if common_exposure.is_infinite() {
        // Every selected wave is completely extinguished. Finite-degree
        // frequency polynomials cannot overcome infinite decay exposure.
        return Ok(zero_derivatives());
    }
    let carrier = Carrier {
        tangential,
        axis: [1.0, 0.0],
        normal: reference_wave.normal.re,
        anchor_z: 0.0,
    };
    for sum in &mut sums {
        sum.reference = Some(carrier);
    }
    let common_factors = core::array::from_fn(|order| {
        phase_factor(
            [
                reference_wave.normal_derivatives[order].re,
                reference_wave.normal_derivatives[order].im,
            ],
            tangential_derivatives[order],
            depth,
            point[0],
            omega,
            order,
        )
    });
    for wave in waves.into_iter().flatten().filter(active) {
        let values = boundary_fields(wave);
        let factors = core::array::from_fn(|order| {
            let normal = difference(
                wave.normal_derivatives[order],
                reference_wave.normal_derivatives[order],
            );
            phase_factor(normal, Dd::ZERO, depth, point[0], omega, order)
        });
        let values = product_jet(values, factors)?;
        // Subtract attenuation roots BEFORE multiplying by distance, retaining
        // weak differential decay under a large common evanescent exposure.
        let relative_exposure = attenuation_exposure_dd(
            Dd::from_f64(wave.normal.im) - Dd::from_f64(reference_wave.normal.im),
            depth,
        )?;
        let geometry = Carrier {
            normal: wave.normal.re,
            ..carrier
        };
        for order in 0..3 {
            let attenuated = attenuated_fields(values[order], relative_exposure);
            let needed_by_later_product = match order {
                0 => common_factors != [ZERO; 2],
                1 => common_factors[0] != ZERO,
                _ => false,
            };
            if needed_by_later_product
                && values[order]
                    .iter()
                    .flatten()
                    .zip(attenuated.iter().flatten())
                    .any(|(&before, &after)| before != Dd::ZERO && after == Dd::ZERO)
            {
                return Err(
                    "relative decay underflowed a field needed by a common frequency product",
                );
            }
            sums[order].add_fields(attenuated, geometry)?;
        }
    }
    let values = product_jet(sums.map(|sum| sum.fields), common_factors)?;
    let derivatives = [
        attenuated_fields(values[1], common_exposure),
        attenuated_fields(values[2], common_exposure),
    ];
    if derivatives.iter().flatten().all(|&value| value == ZERO) {
        return Ok(zero_derivatives());
    }
    let rotation = relative_rotation(
        carrier,
        Carrier {
            tangential: 0.0,
            axis: [0.0; 2],
            normal: 0.0,
            anchor_z: 0.0,
        },
        [point[0], point[1]],
        depth,
    )?;
    let rounded = derivatives.map(|values| {
        values.map(|value| {
            let value = multiply(value, rotation);
            C64::new(value[0].to_f64(), value[1].to_f64())
        })
    });
    if rounded
        .iter()
        .flatten()
        .any(|&value| !super::super::finite(value))
    {
        return Err("crystal spatial frequency field derivatives are unrepresentable");
    }
    Ok(FieldDerivatives {
        electric: rounded.map(|values| [values[0], values[1], values[2]]),
        magnetic: rounded.map(|values| [values[3], values[4], values[5]]),
    })
}

fn active(wave: &FrequencyWave<'_>) -> bool {
    wave.electric
        .iter()
        .chain(wave.magnetic)
        .any(|&value| value != C64::ZERO)
        || wave
            .derivatives
            .electric
            .iter()
            .chain(&wave.derivatives.magnetic)
            .flatten()
            .any(|&value| value != C64::ZERO)
}

fn boundary_fields(wave: FrequencyWave<'_>) -> [Fields; 3] {
    core::array::from_fn(|order| {
        core::array::from_fn(|component| {
            let value = if order == 0 {
                if component < 3 {
                    wave.electric[component]
                } else {
                    wave.magnetic[component - 3]
                }
            } else if component < 3 {
                wave.derivatives.electric[order - 1][component]
            } else {
                wave.derivatives.magnetic[order - 1][component - 3]
            };
            [value.re, value.im].map(Dd::from_f64)
        })
    })
}

fn difference(a: Z, b: Z) -> Complex {
    [a.re - b.re, a.im - b.im]
}

// Minus i times the ordinary phase derivative, scaled by the fixed omega.
// Neither omega² nor a squared phase slope is formed as an f64 intermediate.
fn phase_factor(
    normal: Complex,
    tangential: Dd,
    depth: Dd,
    x: f64,
    omega: f64,
    order: usize,
) -> Complex {
    let mut phase = [
        tangential * Dd::from_f64(x) + normal[0] * depth,
        normal[1] * depth,
    ];
    for _ in 0..=order {
        phase = phase.map(|value| value * Dd::from_f64(omega));
    }
    [phase[1], -phase[0]]
}

fn product_jet(values: [Fields; 3], factors: [Complex; 2]) -> Result<[Fields; 3]> {
    let product = |value: Complex, factor: Complex| {
        if value == ZERO {
            ZERO
        } else {
            multiply(value, factor)
        }
    };
    let add = |a: Complex, b: Complex| [a[0] + b[0], a[1] + b[1]];
    let mut result = values;
    for component in 0..6 {
        let first_product = product(values[0][component], factors[0]);
        result[1][component] = add(values[1][component], first_product);
        let cross = product(values[1][component], factors[0]).map(|v| v * Dd::from_f64(2.0));
        result[2][component] = add(
            add(values[2][component], cross),
            add(
                product(values[0][component], factors[1]),
                product(first_product, factors[0]),
            ),
        );
    }
    if result
        .iter()
        .flatten()
        .flatten()
        .any(|value| !value.hi.is_finite() || !value.lo.is_finite())
    {
        return Err("crystal spatial frequency polynomial is unrepresentable before attenuation");
    }
    Ok(result)
}

fn attenuation_exposure_dd(normal_imaginary: Dd, distance: Dd) -> Result<f64> {
    let leading = -normal_imaginary.hi * distance.hi;
    let exposure = if leading.is_infinite() {
        leading
    } else {
        (-normal_imaginary * distance).to_f64()
    };
    if !(exposure >= 0.0) {
        return Err("crystal relative evanescent frequency field must decay");
    }
    Ok(exposure)
}

fn zero_derivatives() -> FieldDerivatives {
    FieldDerivatives {
        electric: [[C64::ZERO; 3]; 2],
        magnetic: [[C64::ZERO; 3]; 2],
    }
}
