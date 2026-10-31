//! Material dispersion and particle populations in native spectral transport.
//!
//! The accepted material owns bulk density and optical data. The caller owns
//! current particle radii and dispersed mass density. Preparation solves Mie
//! at every wavelength, including independent spectral-midpoint probes; it
//! never freezes an index at a representative colour. Particles are passive,
//! homogeneous spheres with independent, unpolarized scattering.

mod indexed;
pub use indexed::MaterialParticleCellGrid;

use super::bulk_optics::BulkOpticalLaw;
use super::{MaterialOpticalBinding, MaterialOpticalError};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_matdb::PropertyKey;
use fs_material::state_point::{
    ResolvedMaterialStatePoint, VisibleDielectricAbsorption, VisibleOpticalStatePoint,
};
use fs_math::c64::C64;
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{Density, QuantitySpec, Temperature};

use crate::volumes::VolumeGrid;
use crate::volumes::transport::scattering::ScatteringGrid;
use crate::volumes::transport::scattering::mie::{
    MieError, MieOpticalSpectrum, MiePopulation, MiePopulationBudget, MiePopulationSolution,
    MieSphere, ParticleSizeDistribution, prepare_visible_spectrum,
};
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, VolumeTransportError};

/// A frozen, source-bound material and its current particle size distribution.
///
/// Rebind after a change in material state or particle geometry. Radius means
/// the physical current radius, not a reference radius before thermal expansion.
/// No surface roughness, oxide coating, material-name preset or elastic modulus
/// is invented to supply bulk scattering data.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialParticleScattering {
    accepted_state_identity: ContentHash,
    optical_state_identity: ContentHash,
    density_kg_m3: f64,
    temperature_k: Option<f64>,
    optics: BulkOpticalLaw,
    sizes: ParticleSizeDistribution,
    host_index_n: f64,
    identity: ContentHash,
}

impl MaterialParticleScattering {
    /// Select the complete density property key from the ORIGINAL accepted
    /// material bundle. Density and optics must share the card, physical
    /// condition and quantity conventions. The density requirement is a plain
    /// scalar in kg/m³; no mechanics-specific property set is required.
    ///
    /// Dielectric absorption must be a spectrum of bulk INTENSITY absorption,
    /// excluding interface and scattering losses. RGB-derived absorption is
    /// refused. The conversion is `k = alpha * lambda_vacuum / (4*pi)`.
    pub fn try_new(
        accepted: &ResolvedMaterialStatePoint,
        optical: &VisibleOpticalStatePoint,
        density_property: &PropertyKey,
        sizes: ParticleSizeDistribution,
        host_index_n: f64,
    ) -> Result<Self, ParticleScatteringError> {
        MaterialOpticalBinding::validate_bulk_state(accepted, optical.resolved())?;
        let property = accepted
            .property_by_key(density_property)
            .ok_or(ParticleScatteringError::Density)?;
        let requirement = property.requirement();
        let density_kg_m3 = property.value_si();
        if requirement.quantity() != QuantitySpec::dimensional(Density::DIMS)
            || requirement.elastic_component().is_some()
            || requirement.strain_component().is_some()
            || requirement.stress_component().is_some()
            || requirement.hardness_test().is_some()
            || !density_kg_m3.is_finite()
            || density_kg_m3 <= 0.0
        {
            return Err(ParticleScatteringError::Density);
        }
        let absorption = match optical {
            VisibleOpticalStatePoint::Dielectric(value) => Some(value.absorption()),
            VisibleOpticalStatePoint::DielectricSellmeier(value) => Some(value.absorption()),
            _ => None,
        };
        if matches!(
            absorption,
            Some(VisibleDielectricAbsorption::ReferenceRgb { .. })
        ) {
            return Err(ParticleScatteringError::SpectralAbsorptionRequired);
        }
        if !host_index_n.is_finite() || host_index_n <= 0.0 {
            return Err(ParticleScatteringError::HostIndex);
        }
        let optics = BulkOpticalLaw::from_state(optical)?;
        let temperature_k = unique_temperature(accepted);
        let mut hash = DomainHasher::new("org.frankensim.render.material-particles.v1");
        hash.update(accepted.identity().as_bytes());
        hash.update(optical.resolved().identity().as_bytes());
        hash.update(property.answer().receipt.content_hash().as_bytes());
        hash.update(&density_kg_m3.to_bits().to_le_bytes());
        hash.update(&host_index_n.to_bits().to_le_bytes());
        hash.update(&(sizes.bins().len() as u64).to_le_bytes());
        for &(radius, weight) in sizes.bins() {
            hash.update(&radius.to_bits().to_le_bytes());
            hash.update(&weight.to_bits().to_le_bytes());
        }
        Ok(Self {
            accepted_state_identity: accepted.identity(),
            optical_state_identity: optical.resolved().identity(),
            density_kg_m3,
            temperature_k,
            optics,
            sizes,
            host_index_n,
            identity: hash.finalize(),
        })
    }

    /// Wavelength-resolved absolute `n + i*k`, under `exp(-i*omega*t)`.
    pub fn complex_index(&self, vacuum_nm: f64) -> Result<C64, ParticleScatteringError> {
        Ok(self.optics.complex_index(vacuum_nm)?)
    }

    /// Solve the declared population with the ORIGINAL material evaluated at
    /// this vacuum wavelength. Preparation allocation is bounded by the number
    /// of admitted bins; no population loop is introduced at render vertices.
    pub fn at_wavelength(
        &self,
        cx: &Cx<'_>,
        vacuum_nm: f64,
        budget: MiePopulationBudget,
    ) -> Result<MiePopulationSolution, ParticleScatteringError> {
        cx.checkpoint().map_err(|_| MieError::Execution)?;
        if self.sizes.bins().len() > budget.max_bins {
            return Err(MieError::Budget.into());
        }
        let index = self.complex_index(vacuum_nm)?;
        let mut particles = Vec::with_capacity(self.sizes.bins().len());
        for (ordinal, &(radius, number_weight)) in self.sizes.bins().iter().enumerate() {
            if ordinal % 32 == 0 {
                cx.checkpoint().map_err(|_| MieError::Execution)?;
            }
            particles.push((
                MieSphere::try_new(
                    radius,
                    index.re,
                    index.im,
                    self.host_index_n,
                    self.density_kg_m3,
                )?,
                number_weight,
            ));
        }
        Ok(MiePopulation::try_new(&particles)?.at_wavelength(cx, vacuum_nm, budget)?)
    }

    /// Prepare native transport and measured spectral/angular discrepancies.
    /// The shared Mie spectral owner performs every declared knot and fresh
    /// midpoint solve using this material's wavelength-dependent index.
    pub fn visible_spectrum(
        &self,
        cx: &Cx<'_>,
        wavelengths_nm: &[f64],
        budget: MiePopulationBudget,
        angular_tolerance: f64,
        max_phase_knots: usize,
    ) -> Result<MaterialParticleSpectrum, ParticleScatteringError> {
        let spectrum = prepare_visible_spectrum(
            cx,
            self.host_index_n,
            wavelengths_nm,
            angular_tolerance,
            max_phase_knots,
            |wavelength| self.at_wavelength(cx, wavelength, budget),
        )?;
        Ok(MaterialParticleSpectrum {
            binding_identity: self.identity,
            density_kg_m3: self.density_kg_m3,
            temperature_k: self.temperature_k,
            spectrum,
        })
    }

    /// Current particle material density; distinct from dispersed mass per m³.
    #[must_use]
    pub const fn particle_density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }
    /// Original accepted material endpoint.
    #[must_use]
    pub const fn accepted_state_identity(&self) -> ContentHash {
        self.accepted_state_identity
    }
    /// Original optical source bundle at that same endpoint.
    #[must_use]
    pub const fn optical_state_identity(&self) -> ContentHash {
        self.optical_state_identity
    }
    /// Source, state, density, host and complete radius-number distribution.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }
    /// Admitted current radius/number distribution.
    #[must_use]
    pub const fn sizes(&self) -> &ParticleSizeDistribution {
        &self.sizes
    }
}

/// Material-bound native spectra, including the original Mie resolution reports.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialParticleSpectrum {
    binding_identity: ContentHash,
    density_kg_m3: f64,
    temperature_k: Option<f64>,
    spectrum: MieOpticalSpectrum,
}

impl MaterialParticleSpectrum {
    /// Coefficients and angular law for the existing native transport owners.
    /// Nonunit hosts require a matching dielectric-interior binding; this
    /// spectrum alone cannot create an interface or bend rays.
    #[must_use]
    pub const fn spectrum(&self) -> &MieOpticalSpectrum {
        &self.spectrum
    }

    /// Bind dispersed particle mass in kg/m³ to index-one ambient transport.
    ///
    /// `maximum_volume_fraction` is an explicit caller-selected dilute-model
    /// admission limit in (0,1]. It is checked in every cell, but does not prove
    /// independent scattering or exclude correlations between particles.
    /// Optional LTE temperatures must equal the one unambiguous temperature
    /// of the accepted material. Omitting them declares a nonemitting model.
    /// A new composition, material state or radius distribution needs rebinding.
    pub fn bind_ambient_grid<'a>(
        &self,
        cx: &Cx<'_>,
        grid: VolumeGrid<'a>,
        temperature_k: Option<&'a [f64]>,
        maximum_volume_fraction: f64,
        budget: AbsorbingGridBudget,
    ) -> Result<MaterialParticleGrid<'a>, ParticleScatteringError> {
        cx.checkpoint()
            .map_err(|_| VolumeTransportError::Cancelled)?;
        if self.spectrum.host_index_n != 1.0 {
            return Err(ParticleScatteringError::HostIndex);
        }
        if !maximum_volume_fraction.is_finite()
            || maximum_volume_fraction <= 0.0
            || maximum_volume_fraction > 1.0
        {
            return Err(ParticleScatteringError::VolumeFractionLimit);
        }
        if grid.data.len() > budget.max_grid_cells {
            return Err(VolumeTransportError::Budget.into());
        }
        if let Some(temperatures) = temperature_k {
            if temperatures.len() != grid.data.len() {
                return Err(VolumeTransportError::Temperature.into());
            }
            let accepted = self
                .temperature_k
                .ok_or(ParticleScatteringError::ThermalState)?;
            for (chunk, values) in temperatures.chunks(1024).enumerate() {
                cx.checkpoint()
                    .map_err(|_| VolumeTransportError::Cancelled)?;
                if let Some(offset) = values.iter().position(|&value| value != accepted) {
                    return Err(ParticleScatteringError::TemperatureMismatch {
                        index: chunk * 1024 + offset,
                    });
                }
            }
        }
        let mut maximum_observed_volume_fraction = 0.0_f64;
        for (chunk, values) in grid.data.chunks(1024).enumerate() {
            cx.checkpoint()
                .map_err(|_| VolumeTransportError::Cancelled)?;
            for (offset, &density) in values.iter().enumerate() {
                let index = chunk * 1024 + offset;
                if !density.is_finite() || density < 0.0 {
                    return Err(VolumeTransportError::Density { index }.into());
                }
                let fraction = density / self.density_kg_m3;
                if fraction > maximum_volume_fraction {
                    return Err(ParticleScatteringError::VolumeFraction { index });
                }
                maximum_observed_volume_fraction = maximum_observed_volume_fraction.max(fraction);
            }
        }
        let absorption =
            AbsorbingGrid::try_new(cx, grid, self.spectrum.absorption, temperature_k, budget)?;
        let medium = ScatteringGrid::try_new(
            cx,
            absorption,
            self.spectrum.scattering,
            self.spectrum.phase.clone(),
        )?;
        let mut hash = DomainHasher::new("org.frankensim.render.material-particle-grid.v1");
        hash.update(self.binding_identity.as_bytes());
        hash.update(medium.content_identity().as_bytes());
        hash.update(&maximum_volume_fraction.to_bits().to_le_bytes());
        Ok(MaterialParticleGrid {
            medium,
            identity: hash.finalize(),
            maximum_observed_volume_fraction,
        })
    }
}

/// Borrowed native scattering/absorption and optional LTE emission.
pub struct MaterialParticleGrid<'a> {
    medium: ScatteringGrid<'a>,
    identity: ContentHash,
    maximum_observed_volume_fraction: f64,
}

impl<'a> MaterialParticleGrid<'a> {
    /// Direct input to the existing renderer, photometry and beam deposition.
    #[must_use]
    pub const fn medium(&self) -> &ScatteringGrid<'a> {
        &self.medium
    }
    /// Material source and complete numerical field identity.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }
    /// Largest admitted dispersed solid/liquid fraction across original cells.
    #[must_use]
    pub const fn maximum_volume_fraction(&self) -> f64 {
        self.maximum_observed_volume_fraction
    }
}

fn unique_temperature(state: &ResolvedMaterialStatePoint) -> Option<f64> {
    let mut found = None;
    for (name, value) in state.query_point() {
        let value = *value;
        let quantity = state.axis_quantities().get(name);
        let absolute = quantity.and_then(|q| q.semantic_type()).is_some_and(|q| {
            q.kind() == QuantityKind::AbsoluteTemperature && q.form() == ValueForm::Static
        });
        let legacy = name == "T"
            && quantity.is_none_or(|q| *q == QuantitySpec::dimensional(Temperature::DIMS));
        if absolute || legacy {
            if found.is_some() || !value.is_finite() || value <= 0.0 {
                return None;
            }
            found = Some(value);
        }
    }
    found
}

/// Atomic refusal before an inconsistent optical material reaches transport.
#[derive(Clone, Debug, PartialEq)]
pub enum ParticleScatteringError {
    /// Original source material, state, quantities or optical law do not agree.
    Material(MaterialOpticalError),
    /// Missing, nonpositive or incorrectly typed bulk density property.
    Density,
    /// A complex particle index requires spectral bulk absorption, not RGB.
    SpectralAbsorptionRequired,
    /// Invalid host index, or a nonunit host passed to the ambient grid binder.
    HostIndex,
    /// No single positive absolute temperature is defined by the material state.
    ThermalState,
    /// LTE cell temperature differs from the frozen material endpoint.
    TemperatureMismatch { index: usize },
    /// Caller-selected volume fraction ceiling is not in (0,1].
    VolumeFractionLimit,
    /// Dispersed material exceeds the declared volume fraction ceiling.
    VolumeFraction { index: usize },
    /// Original Mie domain, arithmetic, resolution or work-budget refusal.
    Mie(MieError),
    /// Original native transport geometry, coefficient or cancellation refusal.
    Transport(VolumeTransportError),
}

impl core::fmt::Display for ParticleScatteringError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material particle scattering refused: {self:?}")
    }
}
impl core::error::Error for ParticleScatteringError {}
impl From<MaterialOpticalError> for ParticleScatteringError {
    fn from(value: MaterialOpticalError) -> Self {
        Self::Material(value)
    }
}
impl From<MieError> for ParticleScatteringError {
    fn from(value: MieError) -> Self {
        Self::Mie(value)
    }
}
impl From<VolumeTransportError> for ParticleScatteringError {
    fn from(value: VolumeTransportError) -> Self {
        Self::Transport(value)
    }
}

#[cfg(test)]
#[path = "particle_scattering/tests.rs"]
mod tests;
