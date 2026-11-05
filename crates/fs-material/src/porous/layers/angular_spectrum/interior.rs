//! Interior fields of the same coherent refracting layer stack.
//!
//! A point uses the original backward layer matrix and its retained complex
//! normalization. Taking the ratio of the full and remaining layer scalings
//! leaves exp(-i k_n z), never a growing exponential or an underflowed exit
//! amplitude multiplied by one. The original critical-angle matrix limit
//! remains valid without dividing pressure by the normal wavenumber.

use super::{
    AcousticSpectrumError, AcousticSpectrumSample, FieldSum, PlaneField, PreparedAcousticSpectrum,
    PreparedPlaneWave, poll, propagation,
};
use crate::porous::layers::transmission::{LogAmplitude, log_scaled};
use crate::porous::layers::{
    Layer, LayerTransfer, MAX_LAYERS, PreparedLayer, RigidBackedLayers, finite,
};
use fs_math::{c64::C64, dd::Dd, det};

struct LayerWave {
    transfer: LayerTransfer,
    front_amplitude: LogAmplitude,
    tangential_velocity_per_pa: [C64; 2],
}

/// Coherent fields and separate local heat mechanisms in an original layer.
/// Heat is evaluated after summing the fields, not by adding component powers.
/// These equivalent-fluid mechanisms are not a fluid/frame temperature split.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticLayerSample {
    pub field: AcousticSpectrumSample,
    /// `-omega Im(rho) |v|^2 / 2` [W/m³], including tangential velocity.
    pub viscous_heating_w_m3: f64,
    /// `-omega Im(1/K) |p|^2 / 2` [W/m³].
    pub thermal_heating_w_m3: f64,
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
    fields: PreparedLayerFields,
}

// Shared local field/heat owner. Stack-specific boundary solves supply the
// original two face states and complex entrance amplitude, never a fitted loss.
pub(in crate::porous::layers) struct PreparedLayerFields {
    waves: Vec<LayerWave>,
    heat_scales: [f64; 2],
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
        let fields = PreparedLayerFields::prepare(
            prepared,
            layer,
            &self.components,
            cancelled,
            |wave, cancelled| {
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
                            step.normal_wavenumber.im * thickness
                                - log_scaled(8.0, step.normalization),
                        );
                        amplitude.rotate(-step.normal_wavenumber.re * thickness);
                    }
                }
                Ok((transfers[layer_index], amplitude))
            },
        )?;
        Ok(PreparedLayerAcousticSpectrum {
            spectrum: self,
            layer_index,
            fields,
        })
    }
}

impl<'s, 'p, 'a> PreparedLayerAcousticSpectrum<'s, 'p, 'a> {
    pub(super) fn heat_scales(&self) -> [f64; 2] {
        self.fields.heat_scales
    }

    pub(super) fn component_field(
        &self,
        index: usize,
        depth: f64,
    ) -> Result<PlaneField, AcousticSpectrumError> {
        let prepared = self.spectrum.stack.prepared_layers();
        PreparedLayerFields::component_field(
            prepared,
            &prepared.layers[self.layer_index],
            &self.spectrum.components[index],
            &self.fields.waves[index],
            depth,
        )
    }

    /// The original prepared coherent incident/reflected/transmitted spectrum.
    pub const fn spectrum(&self) -> &'s PreparedAcousticSpectrum<'p, 'a> {
        self.spectrum
    }

    /// Selected layer in the original exposed-to-rear ordering.
    pub const fn layer_index(&self) -> usize {
        self.layer_index
    }

    /// Observe the same original coherent field and its local physical heat.
    /// Pressure nodes may still have viscous heat through transverse velocity.
    /// Lossless layers have exactly zero heat; evanescent decay and reactive
    /// intensity are never assigned heat. This does not evolve temperature.
    /// A zero-thickness identity has no material volume and refuses heat, while
    /// [`Self::field`] still admits its original interface observation.
    /// Additional heat work is constant and allocation-free after the field sum.
    pub fn field_with_heat(
        &self,
        position_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcousticLayerSample, AcousticSpectrumError> {
        let prepared = self.spectrum.stack.prepared_layers();
        self.fields.field_with_heat(
            prepared,
            &prepared.layers[self.layer_index],
            &self.spectrum.components,
            position_m,
            cancelled,
        )
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
        let prepared = self.spectrum.stack.prepared_layers();
        self.fields.field(
            prepared,
            &prepared.layers[self.layer_index],
            &self.spectrum.components,
            position_m,
            cancelled,
        )
    }
}

impl PreparedLayerFields {
    pub(in crate::porous::layers) fn prepare<C: FnMut() -> bool>(
        prepared: &RigidBackedLayers<'_>,
        layer: &PreparedLayer<'_>,
        components: &[PreparedPlaneWave],
        cancelled: &mut C,
        mut resolve: impl FnMut(
            &PreparedPlaneWave,
            &mut C,
        ) -> Result<(LayerTransfer, LogAmplitude), AcousticSpectrumError>,
    ) -> Result<Self, AcousticSpectrumError> {
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
        let compliance = layer.porous.map_or(C64::ZERO, |response| {
            response.dynamic_bulk_modulus_pa().recip()
        });
        // Weight before squaring; a weak field can have representable heat.
        let heat_scales = [-density.im, -compliance.im].map(|loss| {
            if loss == 0.0 {
                0.0
            } else {
                det::exp(0.5 * (det::ln(prepared.omega) + det::ln(loss) - core::f64::consts::LN_2))
            }
        });
        let mut waves = Vec::with_capacity(components.len());
        for wave in components {
            poll(cancelled)?;
            let (transfer, front_amplitude) = resolve(wave, cancelled)?;
            let tangential_velocity_per_pa = wave
                .tangential_wavenumber
                .map(|k| C64::from_re(k / prepared.omega) / density);
            if tangential_velocity_per_pa.iter().any(|&v| !finite(v)) {
                return Err(AcousticSpectrumError::Invalid(
                    "interior tangential acoustic momentum is unrepresentable",
                ));
            }
            waves.push(LayerWave {
                transfer,
                front_amplitude,
                tangential_velocity_per_pa,
            });
        }
        poll(cancelled)?;
        Ok(Self { waves, heat_scales })
    }

    pub(in crate::porous::layers) fn field_with_heat(
        &self,
        prepared: &RigidBackedLayers<'_>,
        layer: &PreparedLayer<'_>,
        components: &[PreparedPlaneWave],
        position_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcousticLayerSample, AcousticSpectrumError> {
        poll(cancelled)?;
        if layer.source.thickness_m() == 0.0 {
            return Err(AcousticSpectrumError::Invalid(
                "a zero-thickness layer has no interior heat density",
            ));
        }
        if self.heat_scales.iter().any(|s| !s.is_finite() || *s < 0.0) {
            return Err(AcousticSpectrumError::Invalid(
                "interior passive constitutive heat is unrepresentable",
            ));
        }
        let field = self.field(prepared, layer, components, position_m, cancelled)?;
        let viscous = field
            .particle_velocity_m_s
            .into_iter()
            .fold(Dd::ZERO, |sum, v| {
                sum + Dd::from_f64(v.scale(self.heat_scales[0]).norm_sq())
            })
            .to_f64();
        let thermal = field.pressure_pa.scale(self.heat_scales[1]).norm_sq();
        if !viscous.is_finite() || !thermal.is_finite() {
            return Err(AcousticSpectrumError::Invalid(
                "interior acoustic heat density is unrepresentable",
            ));
        }
        poll(cancelled)?;
        Ok(AcousticLayerSample {
            field,
            viscous_heating_w_m3: viscous,
            thermal_heating_w_m3: thermal,
        })
    }

    pub(in crate::porous::layers) fn field(
        &self,
        prepared: &RigidBackedLayers<'_>,
        layer: &PreparedLayer<'_>,
        components: &[PreparedPlaneWave],
        position_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcousticSpectrumSample, AcousticSpectrumError> {
        poll(cancelled)?;
        let thickness = layer.source.thickness_m();
        let depth = position_m[2];
        if position_m.iter().any(|x| !x.is_finite()) || depth < 0.0 || depth > thickness {
            return Err(AcousticSpectrumError::Invalid(
                "interior observation needs finite coordinates within its selected layer",
            ));
        }
        let mut sum = FieldSum([Dd::ZERO; 8]);
        for (component, local) in components.iter().zip(&self.waves) {
            poll(cancelled)?;
            if local.front_amplitude.log == f64::NEG_INFINITY {
                continue;
            }
            sum.add(
                Self::component_field(prepared, layer, component, local, depth)?,
                propagation(component, [position_m[0], position_m[1], 0.0], C64::ZERO)?,
            )?;
        }
        let sample = sum.finish()?;
        poll(cancelled)?;
        Ok(sample)
    }

    // The same original local propagation serves point and finite-window
    // observations. This is at x=y=0; callers apply the transverse phase.
    fn component_field(
        prepared: &RigidBackedLayers<'_>,
        layer: &PreparedLayer<'_>,
        component: &PreparedPlaneWave,
        local: &LayerWave,
        depth: f64,
    ) -> Result<PlaneField, AcousticSpectrumError> {
        let thickness = layer.source.thickness_m();
        if local.front_amplitude.log == f64::NEG_INFINITY {
            return Ok(PlaneField {
                pressure: C64::ZERO,
                velocity: [C64::ZERO; 3],
            });
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
            let (p, v, scale, kn) = prepared.propagate_layer(&partial, component.cosine, p0, v0)?;
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
        Ok(PlaneField {
            pressure: amplitude.coefficient(pressure, 0.0)?,
            velocity: [
                amplitude.coefficient(pressure * local.tangential_velocity_per_pa[0], 0.0)?,
                amplitude.coefficient(pressure * local.tangential_velocity_per_pa[1], 0.0)?,
                amplitude.coefficient(
                    normal_velocity,
                    det::ln(component.cosine) - det::ln(prepared.z0),
                )?,
            ],
        })
    }
}
