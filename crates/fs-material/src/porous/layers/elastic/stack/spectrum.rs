//! Coherent fluid fields and porous heat in the original elastic plate stack.
//!
//! The original P/S boundary solve supplies section face states. The shared
//! angular-spectrum owner supplies exterior momentum, local fluid propagation,
//! coherent summation and irreversible heat. No solid displacement field or
//! constitutive damping is inferred from the lossless plate sources.

use super::{ElasticPlateStack, PreparedSection, SectionTransfer};
use crate::porous::layers::MAX_LAYERS;
use crate::porous::layers::angular_spectrum::{
    AcousticLayerSample, AcousticSpectrumError, AcousticSpectrumSample, AngularSpectrumBudget,
    IncidentPlaneWave, PreparedLayerFields, PreparedPlaneWave, poll, prepare_components,
    sample_components,
};
use crate::porous::layers::transmission::LogAmplitude;
use fs_math::c64::C64;

/// Coherent exterior waves through one original fluid-coupled elastic stack.
/// Preparation retains O(components) storage and performs one original section
/// solve per component. Observations allocate nothing and re-query no source.
/// The same frequency, planar geometry and lossless solid assumptions remain.
/// Pressure weights scale the infinitesimal response; fluid observations do
/// not establish admission of a finite strain in the original solid sources.
pub struct PreparedElasticAcousticSpectrum<'p, 'a> {
    stack: &'p ElasticPlateStack<'a>,
    components: Vec<PreparedPlaneWave>,
}

/// Coherent pressure, three-component fluid velocity and local porous heat in
/// one selected original fluid section, before, between or after elastic plates.
/// Optional preparation costs O(components × sections) and stores O(components).
/// Original section indexing includes plates and zero-thickness identities.
pub struct PreparedElasticLayerAcousticSpectrum<'s, 'p, 'a> {
    spectrum: &'s PreparedElasticAcousticSpectrum<'p, 'a>,
    layer_index: usize,
    fluid_index: usize,
    fields: PreparedLayerFields,
}

impl<'a> ElasticPlateStack<'a> {
    /// Retain coherent incident/reflected/transmitted pressure and velocity from
    /// this stack's original oblique P/S solution. Incident weights are complex
    /// peak pressures, never incoherent powers. Every section transfer polls
    /// cancellation, including plates and zero-thickness sections.
    pub fn prepare_angular_spectrum<'p>(
        &'p self,
        incident: &[IncidentPlaneWave],
        budget: AngularSpectrumBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedElasticAcousticSpectrum<'p, 'a>, AcousticSpectrumError> {
        let components = prepare_components(
            &self.prepared,
            self.receiver.density_kg_m3(),
            incident,
            budget,
            cancelled,
            |cosine, cancelled| {
                self.scattering_with_sections::<AcousticSpectrumError>(cosine, |_, _| {
                    poll(cancelled)
                })
            },
        )?;
        Ok(PreparedElasticAcousticSpectrum {
            stack: self,
            components,
        })
    }
}

impl<'p, 'a> PreparedElasticAcousticSpectrum<'p, 'a> {
    pub const fn stack(&self) -> &'p ElasticPlateStack<'a> {
        self.stack
    }

    pub fn component_count(&self) -> usize {
        self.components.len()
    }

    /// Total incident plus reflected field; z <= 0 is relative to the entrance.
    pub fn incident_field(
        &self,
        position_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcousticSpectrumSample, AcousticSpectrumError> {
        sample_components(&self.components, position_m, false, cancelled)
    }

    /// Receiving-fluid field; z >= 0 is relative to the final exit face.
    /// Evanescent transmission remains a reactive field, not absorbed heat.
    pub fn transmitted_field(
        &self,
        position_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcousticSpectrumSample, AcousticSpectrumError> {
        sample_components(&self.components, position_m, true, cancelled)
    }

    /// Select an original fluid/porous section. Plate indices refuse because
    /// fluid pressure/velocity and pore heat do not describe a solid interior.
    /// Each preparation reuses the original bounded section solver and cached
    /// material snapshots. Both passes poll cancellation per section.
    pub fn prepare_layer_field(
        &self,
        layer_index: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedElasticLayerAcousticSpectrum<'_, 'p, 'a>, AcousticSpectrumError> {
        poll(cancelled)?;
        let Some(PreparedSection::Fluid {
            index: fluid_index, ..
        }) = self.stack.sections.get(layer_index)
        else {
            return Err(AcousticSpectrumError::Invalid(
                "elastic-stack interior observation needs an original fluid section",
            ));
        };
        let prepared = &self.stack.prepared;
        let fields = PreparedLayerFields::prepare(
            prepared,
            &prepared.layers[*fluid_index],
            &self.components,
            cancelled,
            |wave, cancelled| {
                let mut transfers = [SectionTransfer::default(); MAX_LAYERS];
                self.stack
                    .scattering_with_sections::<AcousticSpectrumError>(
                        wave.cosine,
                        |index, step| {
                            poll(cancelled)?;
                            transfers[index] = step;
                            Ok(())
                        },
                    )?;
                let first = transfers[0].fields;
                let mut amplitude = LogAmplitude::new(0.0);
                amplitude
                    .multiply(C64::from_re(2.0) / (first.front_pressure + first.front_velocity))?;
                amplitude.multiply(wave.incident_pressure())?;
                for (section, step) in self.stack.sections[..layer_index]
                    .iter()
                    .zip(&transfers[..layer_index])
                {
                    poll(cancelled)?;
                    for factor in step.log_factors {
                        amplitude.add(factor);
                    }
                    if step.negative_transmission {
                        amplitude.multiply(C64::from_re(-1.0))?;
                    }
                    if let PreparedSection::Fluid { source, .. } = section
                        && step.fields.normal_wavenumber.im != 0.0
                    {
                        amplitude.rotate(-step.fields.normal_wavenumber.re * source.thickness_m());
                    }
                }
                Ok((transfers[layer_index].fields, amplitude))
            },
        )?;
        Ok(PreparedElasticLayerAcousticSpectrum {
            spectrum: self,
            layer_index,
            fluid_index: *fluid_index,
            fields,
        })
    }
}

impl<'s, 'p, 'a> PreparedElasticLayerAcousticSpectrum<'s, 'p, 'a> {
    pub const fn spectrum(&self) -> &'s PreparedElasticAcousticSpectrum<'p, 'a> {
        self.spectrum
    }

    pub const fn layer_index(&self) -> usize {
        self.layer_index
    }

    /// Depth z is measured from this section's front face. Both faces are
    /// included; x/y keep the original entrance-plane transverse origin.
    /// Zero-thickness fluid identities admit only z=0. Queries are O(components)
    /// and allocation-free, retaining standing waves and evanescent tunneling.
    pub fn field(
        &self,
        position_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcousticSpectrumSample, AcousticSpectrumError> {
        let prepared = &self.spectrum.stack.prepared;
        self.fields.field(
            prepared,
            &prepared.layers[self.fluid_index],
            &self.spectrum.components,
            position_m,
            cancelled,
        )
    }

    /// Sum complete complex fields before forming local viscous/thermal heat
    /// [W/m³]. Interference can suppress pressure while transverse pore velocity
    /// still heats. Lossless fluid sections have exactly zero heat; identities
    /// have no volume and refuse. This does not evolve temperature or assign a
    /// fluid/frame thermal partition. Transmitted and reactive power are not heat.
    pub fn field_with_heat(
        &self,
        position_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcousticLayerSample, AcousticSpectrumError> {
        let prepared = &self.spectrum.stack.prepared;
        self.fields.field_with_heat(
            prepared,
            &prepared.layers[self.fluid_index],
            &self.spectrum.components,
            position_m,
            cancelled,
        )
    }
}
