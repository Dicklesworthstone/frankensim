//! Original fluid and porous material laws in the same resolved 3-D pressure field.
//!
//! Each original tetrahedron names a source region. Free-fluid regions consume
//! FluidAcoustics; rigid porous regions consume the ORIGINAL JCA/JCAL owner,
//! including its fluid substitution, frequency support and optional molecular
//! calorics. Density and bulk modulus already include porosity. No conjugation,
//! additional porosity factor or independently added air attenuation is applied:
//! both original owners use exp(+i omega t).
//!
//! This is the equivalent-fluid rigid-frame approximation, not elastic Biot
//! waves, anisotropic pores, evolving saturation or a new material qualification.
//! Pressure and superficial normal velocity cross conforming material interfaces
//! through the same full FE weak form. No homogenized boundary impedance replaces
//! a resolved porous volume, and heat remains separate from exterior wall losses.

use fs_feec::acoustics::helmholtz::HelmholtzFe;
use fs_feec::acoustics::helmholtz::heterogeneous::{
    BoundaryAdmittance, EquivalentFluidCoefficients, EquivalentFluidControls, HeterogeneousError,
    PreparedEquivalentFluid,
};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_material::porous::{PorousFluidLaw, PorousFluidResponse};
use fs_material::relaxation::RelaxationState;
use fs_math::c64::C64;

/// Borrowed original sources; every material domain and phase policy stays there.
#[derive(Clone, Copy)]
pub enum AcousticVolumeRegion<'s> {
    /// A lossless adiabatic gas/liquid region with its original EOS closure.
    Fluid(&'s dyn FluidAcoustics),
    /// One rigid porous skeleton filled by one coherent original fluid.
    Porous {
        law: &'s PorousFluidLaw,
        fluid: &'s dyn ThermoviscousFluid,
    },
    /// The same rigid skeleton with joint molecular and pore-wall relaxation.
    /// No free-gas molecular attenuation is added to this effective response.
    RelaxingPorous {
        law: &'s PorousFluidLaw,
        relaxation: &'s RelaxationState,
    },
}

/// Exact frozen source response used throughout a conforming region.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AcousticVolumeRegionState {
    Fluid(FluidAcousticState),
    Porous(PorousFluidResponse),
}

impl AcousticVolumeRegionState {
    /// Original numerical coefficients, in the FE owner's same convention.
    #[must_use]
    pub fn coefficients(self) -> EquivalentFluidCoefficients {
        match self {
            Self::Fluid(state) => EquivalentFluidCoefficients {
                dynamic_density_kg_m3: C64::from_re(state.density_kg_m3()),
                dynamic_bulk_modulus_pa: C64::from_re(state.adiabatic_bulk_modulus_pa()),
            },
            Self::Porous(response) => EquivalentFluidCoefficients {
                dynamic_density_kg_m3: response.dynamic_density_kg_m3(),
                dynamic_bulk_modulus_pa: response.dynamic_bulk_modulus_pa(),
            },
        }
    }
}

/// Original material errors remain attributable to the source region.
#[derive(Debug)]
pub enum PorousRegionsInteriorError {
    ElementCount {
        expected: usize,
        supplied: usize,
    },
    UnknownRegion {
        tetrahedron: usize,
        region: usize,
        source_count: usize,
    },
    Material {
        region: usize,
        source: MaterialError,
    },
    Interior(HeterogeneousError),
}

impl core::fmt::Display for PorousRegionsInteriorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ElementCount { expected, supplied } => write!(f,
                "acoustic material regions need {expected} tetrahedron assignments, received {supplied}"),
            Self::UnknownRegion { tetrahedron, region, source_count } => write!(f,
                "acoustic tetrahedron {tetrahedron} names region {region}, with {source_count} sources supplied"),
            Self::Material { region, source } => write!(f, "acoustic volume region {region}: {source}"),
            Self::Interior(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for PorousRegionsInteriorError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Material { source, .. } => Some(source),
            Self::Interior(error) => Some(error),
            _ => None,
        }
    }
}

impl From<HeterogeneousError> for PorousRegionsInteriorError {
    fn from(error: HeterogeneousError) -> Self {
        Self::Interior(error)
    }
}

/// Source-retaining complete pressure solve for free-fluid/porous region maps.
///
/// Original sources are sampled once per USED region, after complete mapping
/// and FE work/storage admission. Every later drive reuses that exact response
/// and factorization. Reprepare for another frequency or equilibrium state.
pub struct PorousRegionsInterior<'f, 'c, 's> {
    sources: Vec<AcousticVolumeRegion<'s>>,
    element_regions: Vec<usize>,
    states: Vec<Option<AcousticVolumeRegionState>>,
    operator: PreparedEquivalentFluid<'f, 'c>,
}

impl<'f, 'c, 's> PorousRegionsInterior<'f, 'c, 's> {
    pub fn new(
        geometry: &'f HelmholtzFe<'c>,
        frequency_hz: f64,
        sources: &[AcousticVolumeRegion<'s>],
        element_regions: &[usize],
        admittances: &[BoundaryAdmittance],
        pressure_release_faces: &[[u32; 3]],
        controls: EquivalentFluidControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, PorousRegionsInteriorError> {
        if cancelled() {
            return Err(HeterogeneousError::Cancelled.into());
        }
        if element_regions.len() != geometry.tetrahedra().len() {
            return Err(PorousRegionsInteriorError::ElementCount {
                expected: geometry.tetrahedra().len(),
                supplied: element_regions.len(),
            });
        }
        for (tetrahedron, &region) in element_regions.iter().enumerate() {
            if cancelled() {
                return Err(HeterogeneousError::Cancelled.into());
            }
            if region >= sources.len() {
                return Err(PorousRegionsInteriorError::UnknownRegion {
                    tetrahedron,
                    region,
                    source_count: sources.len(),
                });
            }
        }
        let omega = core::f64::consts::TAU * frequency_hz;
        let mut states = vec![None; sources.len()];
        let operator = PreparedEquivalentFluid::new(
            geometry,
            frequency_hz,
            admittances,
            pressure_release_faces,
            controls,
            &mut |element, _| {
                let region = element_regions[element];
                let state = match states[region] {
                    Some(state) => state,
                    None => {
                        let evaluated = match sources[region] {
                            AcousticVolumeRegion::Fluid(fluid) => fluid.acoustic_state()
                                .map(AcousticVolumeRegionState::Fluid),
                            AcousticVolumeRegion::Porous { law, fluid } => law.effective_fluid(fluid, omega)
                                .map(AcousticVolumeRegionState::Porous),
                            AcousticVolumeRegion::RelaxingPorous { law, relaxation } => law.effective_relaxing(relaxation, omega)
                                .map(AcousticVolumeRegionState::Porous),
                        };
                        let state = evaluated.map_err(|source| PorousRegionsInteriorError::Material {
                            region,
                            source,
                        })?;
                        states[region] = Some(state);
                        state
                    }
                };
                Ok::<_, PorousRegionsInteriorError>(state.coefficients())
            },
            cancelled,
        )?;
        Ok(Self {
            sources: sources.to_vec(),
            element_regions: element_regions.to_vec(),
            states,
            operator,
        })
    }

    #[must_use]
    pub fn sources(&self) -> &[AcousticVolumeRegion<'s>] {
        &self.sources
    }

    #[must_use]
    pub fn element_regions(&self) -> &[usize] {
        &self.element_regions
    }

    /// Frozen original state; unknown and unused regions return None.
    #[must_use]
    pub fn region_state(&self, region: usize) -> Option<AcousticVolumeRegionState> {
        self.states.get(region).copied().flatten()
    }

    /// Native complete-field solves and original per-element heat forms.
    #[must_use]
    pub const fn operator(&self) -> &PreparedEquivalentFluid<'f, 'c> {
        &self.operator
    }
}
