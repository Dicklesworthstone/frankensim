//! Material-state ingress shared by all specimen and scene owners.
//!
//! Bulk optics and mechanics must come from one immutable card at one query
//! point. Surface finish is an independent input, not inferred from chemistry.
//! Conversion retains the existing visible-band conductor/Cauchy/Beer-Lambert
//! approximations. The separately declared gray/diffuse thermal binding reuses
//! a conduction emissivity sample at the same complete accepted material state.
//! Neither binding authenticates measurements or supplies experimental validation.

use fs_blake3::{ContentHash, DomainHasher};
use fs_conduction::radiation::SurfaceEmissivity;
use fs_material::state_point::{
    ResolvedMaterialStatePoint, VisibleConductorStatePoint, VisibleDielectricStatePoint,
    VisibleOpticalStatePoint,
};

use crate::conductor::{
    ConductorDataStatus, ConductorError, ConductorIorSample, ConductorOptics, ConductorSource,
    ConductorSurface,
};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricError, DielectricGlass, DielectricSurface,
    GlassProvenance,
};
use crate::thermal_surface::{GrayDiffuseThermalSurface, ThermalSurfaceError};
use crate::tracer::Material;

/// Explicit opaque gray/diffuse observation of a material state. The original
/// emissivity sample is retained, so cooling and rendering can share the exact
/// source and temperature. Gray/diffuse behavior is an additional declared
/// approximation; it is not inferred from a total emissivity measurement.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialThermalBinding {
    accepted_state_identity: ContentHash,
    emissivity: SurfaceEmissivity,
    surface: GrayDiffuseThermalSurface,
}

impl MaterialThermalBinding {
    /// Resolve a rendered thermal material from the same card, complete query
    /// point and quantity conventions as the accepted mechanical/thermal state.
    /// The radiation owner supplies the physical temperature; its coordinate
    /// need not be named `T`. Pressure, phase and other source conditions cannot
    /// be omitted or changed between cooling and observation.
    ///
    /// The selected model emits Planck radiance and reflects `1-epsilon`.
    /// It replaces visible-band conductor/dielectric scattering; combining
    /// those independent models would require a matched spectral/directional
    /// absorptivity law. A later thermal endpoint requires a fresh binding.
    pub fn try_gray_diffuse(
        accepted: &ResolvedMaterialStatePoint,
        emissivity: &SurfaceEmissivity,
        surface_state_identity: ContentHash,
    ) -> Result<Self, MaterialThermalError> {
        let receipt = emissivity.receipt();
        if accepted.card_identity() != emissivity.card_identity()
            || accepted.query_point() != receipt.query_point.as_slice()
            || accepted.axis_quantities() != &receipt.axis_quantities
        {
            return Err(MaterialThermalError::StateMismatch);
        }
        if surface_state_identity == ContentHash([0; 32]) {
            return Err(MaterialThermalError::MissingSurfaceIdentity);
        }
        let mut h = DomainHasher::new("org.frankensim.render.material-thermal-binding.v1");
        h.update(accepted.identity().as_bytes());
        h.update(emissivity.receipt().content_hash().as_bytes());
        h.update(surface_state_identity.as_bytes());
        Ok(Self {
            accepted_state_identity: accepted.identity(),
            emissivity: emissivity.clone(),
            surface: GrayDiffuseThermalSurface::try_new(
                emissivity.temperature_k(),
                emissivity.value(),
                h.finalize(),
            )
            .map_err(MaterialThermalError::Thermal)?,
        })
    }

    /// Direct input to the original tracer, for the actual accepted geometry.
    #[must_use]
    pub const fn material(&self) -> Material {
        Material::GrayDiffuseThermal {
            surface: self.surface,
        }
    }

    /// Exact accepted constitutive state observed by this binding.
    /// A later state requires rebinding before attaching it to physical geometry.
    #[must_use]
    pub const fn accepted_state_identity(&self) -> ContentHash {
        self.accepted_state_identity
    }

    /// Original emissivity source, uncertainty and selected query receipt.
    #[must_use]
    pub const fn emissivity(&self) -> &SurfaceEmissivity {
        &self.emissivity
    }

    /// Physical spectral radiance and the declared gray/diffuse assumption.
    #[must_use]
    pub const fn surface(&self) -> GrayDiffuseThermalSurface {
        self.surface
    }
}

/// Refusal before thermal appearance can alter a scene.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaterialThermalError {
    /// Accepted state and emissivity differ in card, coordinates or quantities.
    StateMismatch,
    /// The caller did not identify the independently declared surface state.
    MissingSurfaceIdentity,
    /// The gray/diffuse constitutive model refused its physical inputs.
    Thermal(ThermalSurfaceError),
}

impl core::fmt::Display for MaterialThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::StateMismatch => {
                f.write_str("thermal appearance requires the accepted card, complete query point and quantities")
            }
            Self::MissingSurfaceIdentity => {
                f.write_str("surface-state identity must not be zero")
            }
            Self::Thermal(error) => error.fmt(f),
        }
    }
}

impl core::error::Error for MaterialThermalError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Thermal(error) => Some(error),
            _ => None,
        }
    }
}

/// Supported bulk optical models, without a fallback to an artistic preset.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum OpticalAppearance {
    /// Opaque complex-index conductor with isotropic GGX surface scattering.
    Conductor {
        /// Admitted wavelength table and source assertion.
        optics: ConductorOptics,
        /// Explicit surface finish.
        surface: ConductorSurface,
    },
    /// Homogeneous Cauchy dispersion and Beer-Lambert absorption.
    Dielectric {
        /// Admitted bulk optical law.
        glass: DielectricGlass,
        /// Explicit smooth or rough surface finish.
        surface: DielectricSurface,
    },
}

impl OpticalAppearance {
    /// Direct input to the existing spectral path tracer.
    #[must_use]
    pub const fn material(self) -> Material {
        match self {
            Self::Conductor { optics, surface } => Material::Conductor { optics, surface },
            Self::Dielectric { glass, surface } => Material::Dielectric { glass, surface },
        }
    }
}

/// Optical observation of the same card/state consumed by a mechanical model.
///
/// The mechanical model owns its requirements: a string need not acquire a
/// Poisson ratio or yield stress merely to be rendered. Identities bind inputs
/// for replay and cannot confer measurement or physical-validation authority.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialOpticalBinding {
    mechanical_state_identity: ContentHash,
    optical_state_identity: ContentHash,
    material_card_identity: ContentHash,
    surface_state_identity: ContentHash,
    appearance: OpticalAppearance,
}

impl MaterialOpticalBinding {
    /// Bind the family already selected by the material resolver, with an
    /// explicit positive isotropic GGX roughness for either family.
    pub fn try_visible(
        mechanical: &ResolvedMaterialStatePoint,
        optical: &VisibleOpticalStatePoint,
        roughness_alpha: f64,
        surface_state_identity: ContentHash,
    ) -> Result<Self, MaterialOpticalError> {
        match optical {
            VisibleOpticalStatePoint::Conductor(state) => {
                Self::try_conductor(mechanical, state, roughness_alpha, surface_state_identity)
            }
            VisibleOpticalStatePoint::Dielectric(state) => Self::try_dielectric(
                mechanical,
                state,
                Some(roughness_alpha),
                surface_state_identity,
            ),
        }
    }

    /// Bind the card's nine complex-index samples to the conductor model.
    pub fn try_conductor(
        mechanical: &ResolvedMaterialStatePoint,
        optical: &VisibleConductorStatePoint,
        roughness_alpha: f64,
        surface_state_identity: ContentHash,
    ) -> Result<Self, MaterialOpticalError> {
        Self::validate_state(mechanical, optical.resolved(), surface_state_identity)?;
        let mut samples = [ConductorIorSample::try_new(380.0, 1.0, 1.0)?; 9];
        for (slot, sample) in samples.iter_mut().zip(optical.samples()) {
            *slot = ConductorIorSample::try_new(sample.wavelength_nm, sample.eta, sample.k)?;
        }
        let source = ConductorSource::try_new(
            optical.resolved().identity(),
            ConductorDataStatus::MaterialStateClaims,
        )?;
        Ok(Self::new(
            mechanical,
            optical.resolved(),
            surface_state_identity,
            OpticalAppearance::Conductor {
                optics: ConductorOptics::try_new(samples, source)?,
                surface: ConductorSurface::try_rough(roughness_alpha)?,
            },
        ))
    }

    /// Bind dielectric optics, converting SI Cauchy coefficients to the
    /// renderer's micrometre coefficient convention. Evaluation wavelengths
    /// remain nanometres. `None` selects a smooth boundary.
    pub fn try_dielectric(
        mechanical: &ResolvedMaterialStatePoint,
        optical: &VisibleDielectricStatePoint,
        roughness_alpha: Option<f64>,
        surface_state_identity: ContentHash,
    ) -> Result<Self, MaterialOpticalError> {
        Self::validate_state(mechanical, optical.resolved(), surface_state_identity)?;
        let [a, b_m2, c_m4] = optical.cauchy_coefficients_si();
        let glass = DielectricGlass::new(
            CauchyIor::try_new(a, b_m2 * 1.0e12, c_m4 * 1.0e24)?,
            BeerLambertAbsorption::try_from_rgb_transmittance(
                optical.reference_transmittance_linear_rgb(),
                optical.reference_distance_m(),
            )?,
            GlassProvenance::Custom,
        );
        let surface = match roughness_alpha {
            Some(alpha) => DielectricSurface::try_rough(alpha)?,
            None => DielectricSurface::SMOOTH,
        };
        Ok(Self::new(
            mechanical,
            optical.resolved(),
            surface_state_identity,
            OpticalAppearance::Dielectric { glass, surface },
        ))
    }

    fn validate_state(
        mechanical: &ResolvedMaterialStatePoint,
        optical: &ResolvedMaterialStatePoint,
        surface: ContentHash,
    ) -> Result<(), MaterialOpticalError> {
        if mechanical.material() != optical.material()
            || mechanical.card_identity() != optical.card_identity()
            || mechanical.query_point() != optical.query_point()
            || mechanical.axis_quantities() != optical.axis_quantities()
        {
            return Err(MaterialOpticalError::StateMismatch);
        }
        if surface == ContentHash([0; 32]) {
            return Err(MaterialOpticalError::MissingSurfaceIdentity);
        }
        Ok(())
    }

    fn new(
        mechanical: &ResolvedMaterialStatePoint,
        optical: &ResolvedMaterialStatePoint,
        surface_state_identity: ContentHash,
        appearance: OpticalAppearance,
    ) -> Self {
        Self {
            mechanical_state_identity: mechanical.identity(),
            optical_state_identity: optical.identity(),
            material_card_identity: mechanical.card_identity(),
            surface_state_identity,
            appearance,
        }
    }

    /// Validated bulk optical model and independent surface finish.
    #[must_use]
    pub const fn appearance(self) -> OpticalAppearance {
        self.appearance
    }

    /// Exact accepted mechanical state used to admit the optical observation.
    #[must_use]
    pub const fn accepted_state_identity(self) -> ContentHash {
        self.mechanical_state_identity
    }

    /// Direct input to any existing tracer primitive, independent of specimen.
    #[must_use]
    pub const fn material(self) -> Material {
        self.appearance.material()
    }

    /// Identity of both state bundles, the surface input and rendered material.
    #[must_use]
    pub fn identity(self) -> ContentHash {
        let mut hasher = DomainHasher::new("org.frankensim.render.material-optical-binding.v1");
        for identity in [
            self.mechanical_state_identity,
            self.optical_state_identity,
            self.material_card_identity,
            self.surface_state_identity,
            self.material().content_identity(),
        ] {
            hasher.update(identity.as_bytes());
        }
        hasher.finalize()
    }
}

/// Named refusal before any scene or accepted physical state is changed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MaterialOpticalError {
    /// Bulk mechanics and optics refer to different material cards or states.
    StateMismatch,
    /// Surface finish has no nonzero caller identity.
    MissingSurfaceIdentity,
    /// The existing conductor model refused its parameters.
    Conductor(ConductorError),
    /// The existing dielectric model refused its parameters.
    Dielectric(DielectricError),
}

impl core::fmt::Display for MaterialOpticalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::StateMismatch => write!(
                f,
                "mechanics and optics require the same card, query point and quantities"
            ),
            Self::MissingSurfaceIdentity => write!(f, "surface-state identity must not be zero"),
            Self::Conductor(error) => error.fmt(f),
            Self::Dielectric(error) => error.fmt(f),
        }
    }
}

impl core::error::Error for MaterialOpticalError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Conductor(error) => Some(error),
            Self::Dielectric(error) => Some(error),
            _ => None,
        }
    }
}

impl From<ConductorError> for MaterialOpticalError {
    fn from(error: ConductorError) -> Self {
        Self::Conductor(error)
    }
}

impl From<DielectricError> for MaterialOpticalError {
    fn from(error: DielectricError) -> Self {
        Self::Dielectric(error)
    }
}
