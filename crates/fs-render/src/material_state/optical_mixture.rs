//! Coexisting accepted gas and particle populations on one physical grid.
//!
//! Fractions describe component MASS, while the borrowed grid contains TOTAL
//! mass per bulk volume in kg/m³. The exact native sums give
//! `sigma_a/s(lambda) = rho_total * sum(w_i * kappa_a/s_i(lambda))`.
//! Angular response is weighted by each component's scattering at the queried
//! wavelength; particle indices, asymmetry parameters and endpoint mixtures
//! are never averaged to invent a replacement optical material.

/// True incident-beam absorption assigned to the original constituents.
pub mod deposition;

use super::gas_scattering::MaterialGasScatteringBinding;
use super::particle_scattering::MaterialParticleSpectrum;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::scattering::{
    ScatteringGrid, ScatteringMixtureComponent, ScatteringOptics,
};
use crate::volumes::transport::{
    AbsorbingGrid, AbsorbingGridBudget, CellPalette, CellPaletteBudget,
    MAX_OPTICAL_MIXTURE_COMPONENTS, MassAbsorptionSpectrum, OpticalMixtureBudget,
    ThermalEmissionComponent, ThermalEmissionModel, VolumeTransportError,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_material::gas_scattering::RayleighGasState;
use fs_qty::semantic::{Composition, CompositionBasis};

/// An original accepted optical owner; no source data or preparation reports
/// are copied, replaced, or inferred from a material name.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MaterialOpticalComponent<'a> {
    /// One complete molecular composition. At most one gas owner per recipe;
    /// molecular species mixing remains with the original gas material owner.
    Gas(&'a MaterialGasScatteringBinding),
    /// One prepared material and particle size distribution.
    Particles(&'a MaterialParticleSpectrum),
}

/// Explicit thermal-source choice for the original accepted mixture.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MaterialMixtureEmission<'a> {
    /// Absorption and scattering only; no thermal radiation is added.
    None,
    /// One borrowed positive temperature per native cell. Every original gas
    /// and particle endpoint must agree with that common LTE temperature.
    CommonLte(&'a [f64]),
    /// Each internally equilibrated particle population emits at its OWN
    /// accepted material temperature. Populations need not be in mutual thermal
    /// equilibrium. The original gas temperature remains independent; its
    /// current nonabsorbing molecular law contributes no thermal source.
    /// An absorbing particle with no accepted temperature refuses admission.
    ConstituentTemperatures,
}

/// Fixed mass fractions and original component references, prepared once.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialOpticalMixture<'a> {
    components: &'a [MaterialOpticalComponent<'a>],
    composition: &'a Composition,
    absorption: MassAbsorptionSpectrum,
    optics: ScatteringOptics,
    source_identity: ContentHash,
}

impl<'a> MaterialOpticalMixture<'a> {
    /// Prepare exact weighted sums of the admitted component laws.
    ///
    /// `composition` must use `MassFraction`, have one positive entry per
    /// component, and close to unit mass within floating-point roundoff.
    /// Values are retained exactly, never silently renormalized. Different
    /// wavelength grids are merged without resampling; the native limit of
    /// 16 distinct spectral knots is an explicit refusal, not a lossy fit.
    /// Molecular inverse-fourth scattering remains analytic.
    ///
    /// This is the existing index-one, independent-scattering model. Every
    /// particle solve must have an index-one host. The mixture neither changes
    /// that host nor supplies coherent refraction, particle correlations, or
    /// dependent-scattering corrections. Grid admission accounts for the
    /// carrier volume displaced by particles using additive phase volumes.
    pub fn try_new(
        cx: &Cx<'_>,
        components: &'a [MaterialOpticalComponent<'a>],
        composition: &'a Composition,
        budget: OpticalMixtureBudget,
    ) -> Result<Self, MaterialOpticalMixtureError> {
        checkpoint(cx)?;
        if components.len() > budget.max_components
            || components.len() > MAX_OPTICAL_MIXTURE_COMPONENTS
        {
            return Err(VolumeTransportError::Budget.into());
        }
        let fractions = composition.fractions();
        if components.is_empty()
            || composition.basis() != CompositionBasis::MassFraction
            || fractions.len() != components.len()
        {
            return Err(MaterialOpticalMixtureError::Composition);
        }
        let sum: f64 = fractions.iter().sum();
        if fractions.iter().any(|&value| value <= 0.0)
            || (sum - 1.0).abs() > 32.0 * f64::EPSILON * components.len() as f64
        {
            return Err(MaterialOpticalMixtureError::Composition);
        }
        let mut absorption = Vec::with_capacity(components.len());
        let mut scattering = Vec::with_capacity(components.len());
        let mut gas_component = None;
        let mut hash = DomainHasher::new("org.frankensim.render.optical-mixture-sources.v1");
        hash.update(&(components.len() as u64).to_le_bytes());
        for (index, (component, &mass_weight)) in components.iter().zip(fractions).enumerate() {
            checkpoint(cx)?;
            hash.update(&mass_weight.to_bits().to_le_bytes());
            let (a, s, phase) = match component {
                MaterialOpticalComponent::Gas(gas) => {
                    if gas_component.replace(index).is_some() {
                        return Err(MaterialOpticalMixtureError::MultipleGases);
                    }
                    hash.update(&[0]);
                    hash.update(gas.state().content_identity().as_bytes());
                    (
                        MassAbsorptionSpectrum::constant(0.0)?,
                        gas.scattering_spectrum(),
                        gas.phase().clone(),
                    )
                }
                MaterialOpticalComponent::Particles(particles) => {
                    if particles.spectrum().host_index_n != 1.0 {
                        return Err(MaterialOpticalMixtureError::HostIndex { component: index });
                    }
                    hash.update(&[1]);
                    hash.update(particles.binding_identity().as_bytes());
                    (
                        particles.spectrum().absorption,
                        particles.spectrum().scattering,
                        particles.spectrum().phase.clone(),
                    )
                }
            };
            absorption.push((mass_weight, a));
            scattering.push(ScatteringMixtureComponent {
                mass_weight,
                scattering: s,
                phase,
            });
        }
        let absorption = MassAbsorptionSpectrum::weighted_sum(cx, &absorption, budget)?;
        let optics = ScatteringOptics::weighted_sum(cx, &scattering, budget)?;
        Ok(Self {
            components,
            composition,
            absorption,
            optics,
            source_identity: hash.finalize(),
        })
    }

    /// Original component slice, including each particle's resolution reports.
    #[must_use]
    pub const fn components(&self) -> &'a [MaterialOpticalComponent<'a>] {
        self.components
    }

    /// Original typed mass composition, with unchanged fraction bits.
    #[must_use]
    pub const fn composition(&self) -> &'a Composition {
        self.composition
    }

    /// Original component bindings and declared mass fractions. Spectral
    /// preparation resolution is additionally bound by the native grid identity.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    /// Exact sum of the original admitted mass absorption laws, in m²/kg total.
    #[must_use]
    pub const fn absorption_spectrum(&self) -> MassAbsorptionSpectrum {
        self.absorption
    }

    /// Exact summed scattering and its wavelength-weighted normalized phase.
    #[must_use]
    pub const fn scattering_optics(&self) -> &ScatteringOptics {
        &self.optics
    }

    fn thermal_emission_model(
        &self,
        cx: &Cx<'_>,
    ) -> Result<ThermalEmissionModel, MaterialOpticalMixtureError> {
        let mut components = Vec::with_capacity(self.components.len());
        for (component, &mass_weight) in self.components.iter().zip(self.composition.fractions()) {
            checkpoint(cx)?;
            let (absorption, temperature_k) = match component {
                MaterialOpticalComponent::Gas(_) => (MassAbsorptionSpectrum::constant(0.0)?, None),
                MaterialOpticalComponent::Particles(particles) => (
                    particles.spectrum().absorption,
                    particles.accepted_temperature_k(),
                ),
            };
            components.push(ThermalEmissionComponent {
                mass_weight,
                absorption,
                temperature_k,
            });
        }
        Ok(ThermalEmissionModel::try_new(
            cx,
            &components,
            OpticalMixtureBudget {
                max_components: self.components.len(),
            },
        )?)
    }
}

/// Coexisting materials in native cells, retaining the original source palette
/// and borrowing the simulation's total density and accepted gas states.
pub struct MaterialOpticalMixtureGrid<'a> {
    mixtures: CellPalette<MaterialOpticalMixture<'a>>,
    density_kg_m3: &'a [f64],
    gas_states: Option<&'a [RayleighGasState]>,
    medium: ScatteringGrid<'a>,
    identity: ContentHash,
    maximum_volume_fraction: f64,
}

impl<'a> MaterialOpticalMixtureGrid<'a> {
    /// Bind recipes to TOTAL mass density in the grid's original x-major order.
    /// Absorption and scattering share the same immutable cell-index map.
    ///
    /// Each mapped gas recipe requires the corresponding original accepted gas
    /// state with the same fixed model. Its EOS density is intrinsic gas mass
    /// per GAS volume. With aggregate particle volume fraction `phi`, admission
    /// requires positive gas volume and
    /// `rho_total * w_gas = rho_gas_intrinsic * (1 - phi)` within eight relative
    /// machine epsilons. Neither original state nor fraction is modified.
    /// Consistent additive volumes give
    /// `rho_total = 1 / (w_gas/rho_gas + sum(w_i/rho_bulk_i))`.
    /// Gas-free recipes remain dispersed particles in index-one ambient/vacuum;
    /// they do not consume an optional gas state. A supplied state buffer must
    /// match the whole grid.
    ///
    /// Particle volume fractions ADD: `sum(rho_total * w_i / rho_bulk_i)` must
    /// stay below the declared limit in (0,1]. This admission does not prove
    /// independent-scattering accuracy. `temperature_k = None` explicitly
    /// declares a nonemitting medium. Otherwise every present particle endpoint
    /// and original gas state must have that same positive LTE temperature;
    /// different-temperature sources are not collapsed into an invented mean.
    pub fn try_new(
        cx: &Cx<'_>,
        grid: VolumeGrid<'a>,
        mixtures: CellPalette<MaterialOpticalMixture<'a>>,
        gas_states: Option<&'a [RayleighGasState]>,
        temperature_k: Option<&'a [f64]>,
        maximum_volume_fraction: f64,
        budget: AbsorbingGridBudget,
    ) -> Result<Self, MaterialOpticalMixtureError> {
        Self::try_new_with_emission(
            cx,
            grid,
            mixtures,
            gas_states,
            temperature_k.map_or(
                MaterialMixtureEmission::None,
                MaterialMixtureEmission::CommonLte,
            ),
            maximum_volume_fraction,
            budget,
        )
    }

    /// Bind the original mixture with an explicit thermal-source model.
    ///
    /// Density, phase-volume closure, material endpoints and dilution obey the
    /// same admission as [`Self::try_new`]. `None` and `CommonLte` retain that
    /// constructor's native paths and temperature checks.
    ///
    /// `ConstituentTemperatures` prepares each recipe's native source directly
    /// from its original particle temperatures and mass-weighted absorption:
    /// `j_lambda = rho_total * sum(w_i * kappa_a_i(lambda) * B_lambda(T_i))`.
    /// No averaged temperature, additional emissivity, gas heating, or spectral
    /// continuation is inferred. Sources share the original cell map; no
    /// per-cell temperature field is created. A changed particle temperature
    /// requires a freshly accepted material spectrum, as for its optical law.
    pub fn try_new_with_emission(
        cx: &Cx<'_>,
        grid: VolumeGrid<'a>,
        mixtures: CellPalette<MaterialOpticalMixture<'a>>,
        gas_states: Option<&'a [RayleighGasState]>,
        emission: MaterialMixtureEmission<'a>,
        maximum_volume_fraction: f64,
        budget: AbsorbingGridBudget,
    ) -> Result<Self, MaterialOpticalMixtureError> {
        checkpoint(cx)?;
        let temperature_k = match emission {
            MaterialMixtureEmission::CommonLte(temperatures) => Some(temperatures),
            _ => None,
        };
        let constituent_temperatures =
            matches!(emission, MaterialMixtureEmission::ConstituentTemperatures);
        if grid.data.len() > budget.max_grid_cells || budget.max_steps_per_segment == 0 {
            return Err(VolumeTransportError::Budget.into());
        }
        if mixtures.cell_indices().len() != grid.data.len() {
            return Err(VolumeTransportError::Palette.into());
        }
        if gas_states.is_some_and(|states| states.len() != grid.data.len()) {
            return Err(MaterialOpticalMixtureError::StateCount);
        }
        if temperature_k.is_some_and(|values| values.len() != grid.data.len()) {
            return Err(VolumeTransportError::Temperature.into());
        }
        if !maximum_volume_fraction.is_finite()
            || maximum_volume_fraction <= 0.0
            || maximum_volume_fraction > 1.0
        {
            return Err(MaterialOpticalMixtureError::VolumeFractionLimit);
        }
        let density_kg_m3 = grid.data;
        let mut maximum_observed_fraction = 0.0_f64;
        for (cell, &density) in density_kg_m3.iter().enumerate() {
            if cell % 1024 == 0 {
                checkpoint(cx)?;
            }
            if !density.is_finite() || density < 0.0 {
                return Err(VolumeTransportError::Density { index: cell }.into());
            }
            let mixture = mixtures
                .at_cell(cell)
                .ok_or(VolumeTransportError::Palette)?;
            let temperature = temperature_k.map(|values| values[cell]);
            if temperature.is_some_and(|value| !value.is_finite() || value <= 0.0) {
                return Err(VolumeTransportError::Temperature.into());
            }
            let mut volume_fraction = 0.0;
            let mut gas_density = None;
            for (component, (material, &fraction)) in mixture
                .components
                .iter()
                .zip(mixture.composition.fractions())
                .enumerate()
            {
                let partial_density = checked_product(density, fraction)?;
                match material {
                    MaterialOpticalComponent::Gas(binding) => {
                        let state = gas_states
                            .ok_or(MaterialOpticalMixtureError::MissingGasStates)?
                            .get(cell)
                            .ok_or(MaterialOpticalMixtureError::StateCount)?;
                        if state.model_identity() != binding.state().model_identity() {
                            return Err(MaterialOpticalMixtureError::GasModelMismatch { cell });
                        }
                        // The EOS density is intrinsic. Compare its BULK-volume
                        // contribution only after every particle volume is known.
                        gas_density = Some((partial_density, state.gas().density));
                        if temperature.is_some_and(|value| value != state.gas().temperature) {
                            return Err(MaterialOpticalMixtureError::TemperatureMismatch {
                                cell,
                                component,
                            });
                        }
                    }
                    MaterialOpticalComponent::Particles(particles) => {
                        let fraction = partial_density / particles.particle_density_kg_m3();
                        if !fraction.is_finite() || (partial_density > 0.0 && fraction == 0.0) {
                            return Err(VolumeTransportError::Unrepresentable.into());
                        }
                        volume_fraction += fraction;
                        if let Some(temperature) = temperature {
                            let accepted = particles.accepted_temperature_k().ok_or(
                                MaterialOpticalMixtureError::ThermalState { cell, component },
                            )?;
                            if temperature != accepted {
                                return Err(MaterialOpticalMixtureError::TemperatureMismatch {
                                    cell,
                                    component,
                                });
                            }
                        }
                    }
                }
            }
            if volume_fraction > maximum_volume_fraction {
                return Err(MaterialOpticalMixtureError::VolumeFraction { cell });
            }
            if let Some((partial_density, intrinsic_density)) = gas_density {
                let gas_volume = 1.0 - volume_fraction;
                if gas_volume <= 0.0 {
                    return Err(MaterialOpticalMixtureError::GasFreeVolume { cell });
                }
                let accepted = checked_product(intrinsic_density, gas_volume)?;
                if (partial_density - accepted).abs()
                    > 8.0 * f64::EPSILON * partial_density.max(accepted)
                {
                    return Err(MaterialOpticalMixtureError::GasDensityMismatch { cell });
                }
            }
            maximum_observed_fraction = maximum_observed_fraction.max(volume_fraction);
        }
        let palette_budget = CellPaletteBudget {
            max_models: mixtures.values().len(),
            max_cells: mixtures.cell_indices().len(),
        };
        let mut absorption = Vec::new();
        let mut thermal_sources = Vec::new();
        if constituent_temperatures {
            thermal_sources
                .try_reserve_exact(palette_budget.max_models)
                .map_err(|_| VolumeTransportError::Budget)?;
        } else {
            absorption
                .try_reserve_exact(palette_budget.max_models)
                .map_err(|_| VolumeTransportError::Budget)?;
        }
        let mut scattering = Vec::with_capacity(palette_budget.max_models);
        let mut hash = DomainHasher::new("org.frankensim.render.optical-mixture-grid.v1");
        hash.update(&(mixtures.values().len() as u64).to_le_bytes());
        for mixture in mixtures.values() {
            checkpoint(cx)?;
            if constituent_temperatures {
                thermal_sources.push(mixture.thermal_emission_model(cx)?);
            } else {
                absorption.push(mixture.absorption);
            }
            scattering.push(mixture.optics.clone());
            hash.update(mixture.source_identity.as_bytes());
        }
        hash.update(&[u8::from(gas_states.is_some())]);
        if let Some(states) = gas_states {
            for chunk in states.chunks(1024) {
                checkpoint(cx)?;
                for state in chunk {
                    hash.update(state.content_identity().as_bytes());
                }
            }
        }
        let scattering = mixtures.with_values(cx, scattering, palette_budget)?;
        let absorption = if constituent_temperatures {
            hash.update(b"constituent-temperatures-v1");
            let sources = mixtures.with_values(cx, thermal_sources, palette_budget)?;
            AbsorbingGrid::try_new_thermal_sources(cx, grid, sources, budget)?
        } else {
            let absorption = mixtures.with_values(cx, absorption, palette_budget)?;
            AbsorbingGrid::try_new_indexed(cx, grid, absorption, temperature_k, budget)?
        };
        let medium = ScatteringGrid::try_new_indexed(cx, absorption, scattering)?;
        hash.update(medium.content_identity().as_bytes());
        hash.update(&maximum_volume_fraction.to_bits().to_le_bytes());
        Ok(Self {
            mixtures,
            density_kg_m3,
            gas_states,
            medium,
            identity: hash.finalize(),
            maximum_volume_fraction: maximum_observed_fraction,
        })
    }

    /// Original source recipes and shared physical-cell map.
    #[must_use]
    pub const fn mixtures(&self) -> &CellPalette<MaterialOpticalMixture<'a>> {
        &self.mixtures
    }

    /// Original simulation total-density buffer, kg/m³.
    #[must_use]
    pub const fn density_kg_m3(&self) -> &'a [f64] {
        self.density_kg_m3
    }

    /// Original accepted gas snapshot, when supplied.
    #[must_use]
    pub const fn gas_states(&self) -> Option<&'a [RayleighGasState]> {
        self.gas_states
    }

    /// Direct input to native rendering, photometry and power deposition.
    #[must_use]
    pub const fn medium(&self) -> &ScatteringGrid<'a> {
        &self.medium
    }

    /// Sources, exact native spectra/phases, physical fields and dilution limit.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Largest admitted SUM of particle volume fractions over original cells.
    #[must_use]
    pub const fn maximum_volume_fraction(&self) -> f64 {
        self.maximum_volume_fraction
    }
}

/// Refusal before inconsistent coexisting material data reaches transport.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialOpticalMixtureError {
    /// Wrong basis/count, zero active fraction, or mass closure beyond roundoff.
    Composition,
    /// Multiple gas owners must first be combined by the gas material owner.
    MultipleGases,
    /// A particle population was prepared for a nonunit optical host.
    HostIndex { component: usize },
    /// A mapped gas component has no original accepted cell-state field.
    MissingGasStates,
    /// Gas states do not use the same complete cell layout as density.
    StateCount,
    /// Fixed molecular composition or optical sources changed in this cell.
    GasModelMismatch { cell: usize },
    /// Declared gas partial density disagrees with intrinsic density times its
    /// unoccupied gas-volume fraction.
    GasDensityMismatch { cell: usize },
    /// A positive gas mass fraction has no remaining physical gas volume.
    GasFreeVolume { cell: usize },
    /// A particle endpoint has no unambiguous accepted absolute temperature.
    ThermalState { cell: usize, component: usize },
    /// A constituent endpoint disagrees with the declared common LTE temperature.
    TemperatureMismatch { cell: usize, component: usize },
    /// Declared aggregate particle-volume ceiling is not in (0,1].
    VolumeFractionLimit,
    /// The sum of all particle volume fractions exceeds the admission ceiling.
    VolumeFraction { cell: usize },
    /// Original native spectral, phase, geometry, budget or cancellation refusal.
    Transport(VolumeTransportError),
}

impl core::fmt::Display for MaterialOpticalMixtureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material optical mixture refused: {self:?}")
    }
}
impl core::error::Error for MaterialOpticalMixtureError {}
impl From<VolumeTransportError> for MaterialOpticalMixtureError {
    fn from(value: VolumeTransportError) -> Self {
        Self::Transport(value)
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), MaterialOpticalMixtureError> {
    cx.checkpoint()
        .map_err(|_| VolumeTransportError::Cancelled.into())
}

fn checked_product(a: f64, b: f64) -> Result<f64, MaterialOpticalMixtureError> {
    let result = a * b;
    if !result.is_finite() || (a > 0.0 && b > 0.0 && result == 0.0) {
        return Err(VolumeTransportError::Unrepresentable.into());
    }
    Ok(result)
}

#[cfg(test)]
mod tests;
