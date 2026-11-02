//! Positive per-film work on the original coherent field solution.
//!
//! A symmetric film between real reference-admittance ports has independent
//! even/odd dissipative modes. Their losses multiply |a+b|²/2 and |a-b|²/2,
//! where a,b are incident port waves reconstructed from the original F/G
//! boundary fields. This is the same positive slab-work decomposition as
//! fs-em::plane_wave::layers, in this owner's exp(-i omega t) convention.
//! Poynting/field conventions: Byrnes, sections 4.3--4.5,
//! <https://arxiv.org/abs/1603.02720>. No heat comes from a difference of nearly
//! equal interface fluxes, and no growing propagation exponential is formed.

use super::{
    AbsorbedPower, Complex, ComplexIor, ConductorError, ConductorPowerBalance, Geometry,
    LayerTransfer, MAX_COHERENT_STACK_LAYERS, Reflection, ResponseKind, det, response_as,
};
use crate::conductor::ConductorLayerPowerBalance;

#[derive(Clone, Copy)]
pub(super) struct Layer {
    transfer: LayerTransfer,
    // Polarization -> symmetric/antisymmetric incident-wave loss.
    loss: [[f64; 2]; 2],
}

impl Layer {
    pub(super) const IDENTITY: Self = Self {
        transfer: LayerTransfer::IDENTITY,
        loss: [[0.0; 2]; 2],
    };

    pub(super) fn new(
        ior: ComplexIor,
        thickness: f64,
        geometry: Geometry,
    ) -> Result<Self, ConductorError> {
        let transfer = LayerTransfer::new::<false>(ior, thickness, geometry);
        if ior.k == 0.0 || thickness == 0.0 || geometry.cosine == 0.0 {
            return Ok(Self {
                transfer,
                loss: [[0.0; 2]; 2],
            });
        }
        let n_squared = Complex::new(ior.eta, ior.k).square();
        let q = crate::conductor::thin_film::forward_sqrt(
            n_squared - Complex::new(geometry.transverse_squared, 0.0),
        );
        let depth = core::f64::consts::TAU * (thickness * 1e9) / geometry.wavelength_nm;
        let phase = q.re * depth;
        let decay = q.im * depth;
        let attenuation = det::exp(-decay);
        let (sine, cosine) = det::sin_cos(phase);
        let propagation = Complex::new(attenuation * cosine, attenuation * sine);
        let removed = -det::expm1(-2.0 * decay);
        let sum = Complex::new(1.0, 0.0) + propagation;
        let half_sine = det::sin(0.5 * phase);
        let difference = Complex::new(
            -det::expm1(-decay) + 2.0 * attenuation * half_sine * half_sine,
            -attenuation * sine,
        );
        let incident = [
            geometry.incident_eta * geometry.cosine,
            geometry.cosine / geometry.incident_eta,
        ];
        let mut loss = [[0.0; 2]; 2];
        for pol in 0..2 {
            // Y_layer=q (s), q/n² (p magnetic-field basis). Multiplying
            // numerator/denominator by n² avoids an extra complex division.
            let a = if pol == 0 {
                Complex::new(incident[pol], 0.0)
            } else {
                n_squared * incident[pol]
            };
            let scale = a.re.abs().max(a.im.abs()).max(q.re.abs()).max(q.im.abs());
            if !scale.is_finite() || scale == 0.0 {
                return Err(ConductorError::NumericalFailure);
            }
            let a = Complex::new(a.re / scale, a.im / scale);
            let b = Complex::new(q.re / scale, q.im / scale);
            let (even_work, odd_work) = if phase.abs().max(decay.abs()) < 0.1 {
                // sinhc(y) +/- sinc(x), with the O(d³) odd loss retained.
                // The omitted series term is below binary64 roundoff over
                // this interval. Im(n²)=2 eta k comes from the original law.
                let x = phase * phase;
                let y = decay * decay;
                let sinc_defect =
                    x * (1.0 / 6.0 - x * (1.0 / 120.0 - x * (1.0 / 5040.0 - x / 362880.0)));
                let sinhc_excess =
                    y * (1.0 / 6.0 + y * (1.0 / 120.0 + y * (1.0 / 5040.0 + y / 362880.0)));
                let plus = 2.0 + sinhc_excess - sinc_defect;
                let minus = sinhc_excess + sinc_defect;
                let factor =
                    attenuation * (incident[pol] / scale) * (2.0 * ior.eta * ior.k / scale) * depth;
                if pol == 0 {
                    (factor * plus, factor * minus)
                } else {
                    let normal = q.norm_squared();
                    let lateral = geometry.transverse_squared;
                    (
                        factor * (normal * minus + lateral * plus),
                        factor * (normal * plus + lateral * minus),
                    )
                }
            } else {
                let cross = a * Complex::new(b.re, -b.im);
                (
                    cross.re * removed - 2.0 * cross.im * propagation.im,
                    cross.re * removed + 2.0 * cross.im * propagation.im,
                )
            };
            let denominators = [
                (a * sum + b * difference).norm_squared(),
                (a * difference + b * sum).norm_squared(),
            ];
            for (parity, work) in [even_work, odd_work].into_iter().enumerate() {
                let value = 4.0 * work / denominators[parity];
                if !value.is_finite() || !(0.0..=1.0 + 256.0 * f64::EPSILON).contains(&value) {
                    return Err(ConductorError::NumericalFailure);
                }
                loss[pol][parity] = value.min(1.0);
            }
        }
        Ok(Self { transfer, loss })
    }
}

#[derive(Clone, Copy)]
struct Boundary {
    field: [Complex; 2],
    flux: [Complex; 2],
    log_scale: [f64; 2],
}

struct BoundaryResponse;
impl ResponseKind for BoundaryResponse {
    type Output = (ConductorPowerBalance, Self::History);
    type History = [Boundary; MAX_COHERENT_STACK_LAYERS + 1];
    const BOUNDARY_SCALES: bool = true;

    fn history(_: usize) -> Self::History {
        [Boundary {
            field: [Complex::new(0.0, 0.0); 2],
            flux: [Complex::new(0.0, 0.0); 2],
            log_scale: [0.0; 2],
        }; MAX_COHERENT_STACK_LAYERS + 1]
    }

    fn record(
        history: &mut Self::History,
        index: usize,
        field: [Complex; 2],
        flux: [Complex; 2],
        log_scale: [f64; 2],
    ) {
        history[index] = Boundary {
            field,
            flux,
            log_scale,
        };
    }

    fn finish(
        geometry: Geometry,
        reflected: Reflection,
        field: [Complex; 2],
        flux: [Complex; 2],
        history: Self::History,
    ) -> Result<Self::Output, ConductorError> {
        Ok((
            AbsorbedPower::finish(geometry, reflected, field, flux, ())?,
            history,
        ))
    }
}

pub(super) fn partition(
    geometry: Geometry,
    count: usize,
    mut layer: impl FnMut(usize) -> Result<Layer, ConductorError>,
) -> Result<ConductorLayerPowerBalance, ConductorError> {
    let mut losses = [[[0.0; 2]; 2]; MAX_COHERENT_STACK_LAYERS];
    let (total, boundaries) = response_as::<false, BoundaryResponse>(geometry, count, |index| {
        let prepared = layer(index)?;
        losses[index] = prepared.loss;
        Ok(prepared.transfer)
    })?;
    let mut result = ConductorLayerPowerBalance {
        total,
        substrate_absorbed: [0.0; 2],
        layers: [[0.0; 2]; MAX_COHERENT_STACK_LAYERS],
        layer_count: count,
    };
    if count == 0 || geometry.cosine == 0.0 {
        result.substrate_absorbed = total.absorbed;
        return Ok(result);
    }
    let incident = [
        geometry.incident_eta * geometry.cosine,
        geometry.cosine / geometry.incident_eta,
    ];
    let polarization_count = if geometry.cosine == 1.0 { 1 } else { 2 };
    for pol in 0..polarization_count {
        let y = incident[pol];
        let entrance = boundaries[0].field[pol] * y + boundaries[0].flux[pol];
        let mut log_gain = -det::ln(det::hypot(entrance.re, entrance.im));
        let mut absorbed = 0.0;
        for index in 0..count {
            let left = boundaries[index];
            let right = boundaries[index + 1];
            let next_gain = log_gain - left.log_scale[pol];
            let scale = log_gain.max(next_gain);
            // One common complex phase cancels from the quadratic work.
            // Recover only relative real scales, keeping opaque stacks bounded.
            let a = (left.field[pol] * y + left.flux[pol]) * det::exp(log_gain - scale);
            let b = (right.field[pol] * y - right.flux[pol]) * det::exp(next_gain - scale);
            let even = a + b;
            let odd = a - b;
            let even_norm = det::hypot(even.re, even.im) * det::sqrt(0.5 * losses[index][pol][0]);
            let odd_norm = det::hypot(odd.re, odd.im) * det::sqrt(0.5 * losses[index][pol][1]);
            let norm = det::hypot(even_norm, odd_norm);
            let value = if norm == 0.0 {
                0.0
            } else {
                det::exp(2.0 * (scale + det::ln(norm)))
            };
            if !value.is_finite() || value < 0.0 {
                return Err(ConductorError::NumericalFailure);
            }
            result.layers[index][pol] = value;
            absorbed += value;
            log_gain = next_gain;
        }
        let exit_flux = boundaries[count].flux[pol].re;
        let substrate = det::exp(det::ln(4.0 * y * exit_flux) + 2.0 * log_gain);
        absorbed += substrate;
        let allowance = 512.0 * f64::EPSILON * (count + 1) as f64;
        if !substrate.is_finite()
            || substrate < 0.0
            || !absorbed.is_finite()
            || (absorbed - total.absorbed[pol]).abs()
                > allowance * absorbed.max(total.absorbed[pol])
        {
            return Err(ConductorError::NumericalFailure);
        }
        result.substrate_absorbed[pol] = substrate;
    }
    if polarization_count == 1 {
        result.substrate_absorbed[1] = result.substrate_absorbed[0];
        for layer in &mut result.layers[..count] {
            layer[1] = layer[0];
        }
    }
    Ok(result)
}
