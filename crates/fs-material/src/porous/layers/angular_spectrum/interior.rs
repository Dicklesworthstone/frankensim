//! Interior fields of the same coherent refracting layer stack.
//!
//! A point uses the original backward layer matrix and its retained complex
//! normalization. Taking the ratio of the full and remaining layer scalings
//! leaves exp(-i k_n z), never a growing exponential or an underflowed exit
//! amplitude multiplied by one. The original critical-angle matrix limit
//! remains valid without dividing pressure by the normal wavenumber.

use super::{
    AcousticSpectrumError, AcousticSpectrumSample, FieldSum, PlaneField, PreparedAcousticSpectrum,
    poll, propagation,
};
use crate::porous::layers::transmission::{LogAmplitude, log_scaled};
use crate::porous::layers::{Layer, LayerTransfer, MAX_LAYERS, PreparedLayer, finite};
use fs_math::{c64::C64, dd::Dd, det};

struct LayerWave {
    transfer: LayerTransfer,
    front_amplitude: LogAmplitude,
    tangential_velocity_per_pa: [C64; 2],
}

/// Coherent pressure, velocity and intensity inside one original finite layer.
///
/// Optional preparation retains O(components) storage for this selected layer.
/// Exterior-only spectra retain their original storage and work. Preparation
/// repeats the original bounded layer propagation once per component, with no
/// material queries. Each observation is allocation-free O(components) and
/// uses the original local layer matrix; no new boundary solve is performed.
/// All parent source, planar, frozen-state and phasor assumptions remain.
pub struct PreparedLayerAcousticSpectrum<'s, 'p, 'a> {
    spectrum: &'s PreparedAcousticSpectrum<'p, 'a>,
    layer_index: usize,
    waves: Vec<LayerWave>,
}

impl<'p, 'a> PreparedAcousticSpectrum<'p, 'a> {
    /// Prepare fields in one declared layer, counting zero-thickness layers.
    /// Each side of a material interface is observed through its own layer;
    /// pressure and normal velocity are continuous, tangential velocity need
    /// not be. In porous media velocity is the original equivalent-fluid
    /// volume-flow velocity, not a reconstructed microscopic pore velocity.
    pub fn prepare_layer_field(
        &self,
        layer_index: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedLayerAcousticSpectrum<'_, 'p, 'a>, AcousticSpectrumError> {
        poll(cancelled)?;
        let prepared = self.stack.prepared_layers();
        let layer = prepared
            .layers
            .get(layer_index)
            .ok_or(AcousticSpectrumError::Invalid(
                "interior observation needs an original layer index",
            ))?;
        let density = if let Some(response) = layer.porous {
            response.dynamic_density_kg_m3()
        } else {
            C64::from_re(
                layer
                    .fluid
                    .unwrap_or(prepared.fluid.acoustic())
                    .density_kg_m3(),
            )
        };
        let mut waves = Vec::with_capacity(self.components.len());
        for wave in &self.components {
            poll(cancelled)?;
            let mut transfers = [LayerTransfer::default(); MAX_LAYERS];
            self.stack.scattering_with_layers::<AcousticSpectrumError>(
                wave.cosine,
                |index, step| {
                    poll(cancelled)?;
                    transfers[index] = step;
                    Ok(())
                },
            )?;
            let first = transfers[0];
            let mut amplitude = LogAmplitude::new(0.0);
            amplitude
                .multiply(C64::from_re(2.0) / (first.front_pressure + first.front_velocity))?;
            amplitude.multiply(wave.incident.pressure)?;
            for (source, step) in prepared.layers[..layer_index]
                .iter()
                .zip(&transfers[..layer_index])
            {
                poll(cancelled)?;
                let thickness = source.source.thickness_m();
                if thickness == 0.0 {
                    continue;
                }
                if step.normal_wavenumber.im == 0.0 {
                    amplitude.add(-log_scaled(16.0, step.normalization));
                } else {
                    amplitude.add(
                        step.normal_wavenumber.im * thickness - log_scaled(8.0, step.normalization),
                    );
                    amplitude.rotate(-step.normal_wavenumber.re * thickness);
                }
            }
            let tangential_velocity_per_pa = wave
                .tangential_wavenumber
                .map(|k| C64::from_re(k / prepared.omega) / density);
            if tangential_velocity_per_pa.iter().any(|&v| !finite(v)) {
                return Err(AcousticSpectrumError::Invalid(
                    "interior tangential acoustic momentum is unrepresentable",
                ));
            }
            waves.push(LayerWave {
                transfer: transfers[layer_index],
                front_amplitude: amplitude,
                tangential_velocity_per_pa,
            });
        }
        poll(cancelled)?;
        Ok(PreparedLayerAcousticSpectrum {
            spectrum: self,
            layer_index,
            waves,
        })
    }
}

impl<'s, 'p, 'a> PreparedLayerAcousticSpectrum<'s, 'p, 'a> {
    /// The original prepared coherent incident/reflected/transmitted spectrum.
    pub const fn spectrum(&self) -> &'s PreparedAcousticSpectrum<'p, 'a> {
        self.spectrum
    }

    /// Selected layer in the original exposed-to-rear ordering.
    pub const fn layer_index(&self) -> usize {
        self.layer_index
    }

    /// Pressure and all velocity components at `[x,y,z]` in metres, where
    /// `z` is depth from THIS layer's front face and `0 <= z <= thickness`.
    /// The transverse origin is unchanged from the original entrance plane.
    /// Both faces are included; a zero-thickness identity admits only z=0.
    /// The complete complex fields are summed before forming intensity,
    /// retaining internal standing waves and evanescent tunnelling.
    pub fn field(
        &self,
        position_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcousticSpectrumSample, AcousticSpectrumError> {
        poll(cancelled)?;
        let prepared = self.spectrum.stack.prepared_layers();
        let layer = &prepared.layers[self.layer_index];
        let thickness = layer.source.thickness_m();
        let depth = position_m[2];
        if position_m.iter().any(|x| !x.is_finite()) || depth < 0.0 || depth > thickness {
            return Err(AcousticSpectrumError::Invalid(
                "interior observation needs finite coordinates within its selected layer",
            ));
        }
        let mut sum = FieldSum([Dd::ZERO; 8]);
        for (component, local) in self.spectrum.components.iter().zip(&self.waves) {
            poll(cancelled)?;
            if local.front_amplitude.log == f64::NEG_INFINITY {
                continue;
            }
            let step = local.transfer;
            let mut amplitude = local.front_amplitude.clone();
            let (pressure, normal_velocity) = if depth == 0.0 {
                (step.front_pressure, step.front_velocity)
            } else {
                // An inverse local matrix is well conditioned over an
                // optically short distance. It also preserves a small depth
                // that would round away in thickness-depth for a long layer.
                let forward = -step.normal_wavenumber.im * depth <= 0.5;
                let distance = if forward { depth } else { thickness - depth };
                let partial = PreparedLayer {
                    source: match layer.source {
                        Layer::Porous { law, .. } => Layer::Porous {
                            law,
                            thickness_m: distance,
                        },
                        Layer::FluidGap { .. } => Layer::FluidGap {
                            thickness_m: distance,
                        },
                        Layer::FluidLayer { fluid, .. } => Layer::FluidLayer {
                            fluid,
                            thickness_m: distance,
                        },
                    },
                    ..*layer
                };
                // J T(d) J = T(-d), J=diag(1,-1), so the existing bounded
                // backward matrix also supplies the forward solution.
                let (p0, v0) = if forward {
                    (step.front_pressure, -step.front_velocity)
                } else {
                    (step.back_pressure, step.back_velocity)
                };
                let (p, v, scale, kn) =
                    prepared.propagate_layer(&partial, component.cosine, p0, v0)?;
                if forward {
                    amplitude.add(log_scaled(if kn.im == 0.0 { 16.0 } else { 8.0 }, scale));
                    if kn.im != 0.0 {
                        amplitude.add(-kn.im * depth);
                        amplitude.rotate(kn.re * depth);
                    }
                } else {
                    // Cancel full-minus-remaining distances symbolically.
                    // Even an opaque layer can have a finite entrance field.
                    amplitude.add(det::ln(scale) - det::ln(step.normalization));
                    amplitude.add(kn.im * depth);
                    amplitude.rotate(-kn.re * depth);
                }
                (p, if forward { -v } else { v })
            };
            let field = PlaneField {
                pressure: amplitude.coefficient(pressure, 0.0)?,
                velocity: [
                    amplitude.coefficient(pressure * local.tangential_velocity_per_pa[0], 0.0)?,
                    amplitude.coefficient(pressure * local.tangential_velocity_per_pa[1], 0.0)?,
                    amplitude.coefficient(
                        normal_velocity,
                        det::ln(component.cosine) - det::ln(prepared.z0),
                    )?,
                ],
            };
            sum.add(
                field,
                propagation(component, [position_m[0], position_m[1], 0.0], C64::ZERO)?,
            )?;
        }
        let sample = sum.finish()?;
        poll(cancelled)?;
        Ok(sample)
    }
}
