//! Coherent finite films illuminated from an original absorbing medium.
//!
//! Specify signed real kx and the actual complex peak [Ex,Ey] at the first
//! face. The original source laws and explicit permeabilities own every q,
//! field and dissipation term. A real numerical E/H reference makes bounded
//! slab composition possible without replacing the incident material law.
//! Positive film loss is independent of the entrance interference flux.
//!
//! Preparation is O(films), bounded to 64. Illumination and depth queries
//! allocate nothing and never reevaluate a source law or rebuild the stack.
//! All faces are stationary, planar, isotropic and coherent at one frequency.
//! Source amplitude is local to the entrance, not an inferred distant launch.
//! The angular-spectrum and lossless-port owners remain separate.
//!
//! Poynting and film conventions: S. J. Byrnes, "Multilayer optical
//! calculations", §§4.1, 4.3–4.5 and Appendix B, arXiv:1603.02720v5,
//! <https://arxiv.org/pdf/1603.02720>, with the conjugate time convention.
use super::{DielectricLayer, MAX_LAYERS, bad};
use crate::plane_wave::{DielectricResponse, ReflectionError, positive};
use fs_material::field_storage::linear_magnetic_permeability;
use fs_math::c64::C64;

mod kernel;
pub use kernel::{AbsorbingLayerPoint, AbsorbingLayerReflection};

mod two_sided;
pub use kernel::two_sided::{AbsorbingPortPower, AbsorbingTwoSidedReflection};
pub use two_sided::{AbsorbingTwoSidedIllumination, PreparedAbsorbingTwoSidedLayers};

/// Original sources, frequency, real transverse wave number and prepared
/// coherent fields. Lossy, evanescent and exact-critical films/exits are allowed;
/// the incident wave must carry positive normal power at unit amplitude.
pub struct PlanarAbsorbingDielectricLayers<'a> {
    incident: DielectricResponse<'a>,
    exit: DielectricResponse<'a>,
    layers: Vec<DielectricLayer<'a>>,
    prepared: kernel::Prepared,
}

/// One actual coherent illumination of the retained source-bound stack.
pub struct AbsorbingLayerIllumination<'s, 'a> {
    source: &'s PlanarAbsorbingDielectricLayers<'a>,
    electric: [C64; 2],
    solved: kernel::Illumination,
}

impl<'a> PlanarAbsorbingDielectricLayers<'a> {
    pub fn try_new(
        incident: DielectricResponse<'a>,
        exit: DielectricResponse<'a>,
        omega_rad_s: f64,
        transverse_wavenumber_rad_m: f64,
        layers: &[DielectricLayer<'a>],
    ) -> Result<Self, ReflectionError> {
        if layers.len() > MAX_LAYERS
            || layers
                .iter()
                .any(|l| !l.thickness_m.is_finite() || l.thickness_m < 0.0)
        {
            return Err(bad(
                "finite layers require <=64 finite nonnegative thicknesses",
            ));
        }
        if matches!(incident, DielectricResponse::IonizedGas(_))
            || matches!(exit, DielectricResponse::IonizedGas(_))
            || layers
                .iter()
                .any(|l| matches!(l.response, DielectricResponse::IonizedGas(_)))
        {
            return Err(bad(
                "absorbing-incidence films require original whole-field plasma drift admission",
            ));
        }
        let omega = positive(omega_rad_s)?;
        let input = kernel::Input {
            omega,
            kx: transverse_wavenumber_rad_m,
            epsilon_i: incident.permittivity(omega)?,
            epsilon_t: exit.permittivity(omega)?,
            mu_i: linear_magnetic_permeability(incident.material())?,
            mu_t: linear_magnetic_permeability(exit.material())?,
        };
        let mut resolved = Vec::with_capacity(layers.len());
        for layer in layers {
            resolved.push(kernel::FilmInput {
                epsilon: layer.response.permittivity(omega)?,
                mu: linear_magnetic_permeability(layer.response.material())?,
                thickness: layer.thickness_m,
            });
        }
        let prepared = kernel::Prepared::new(input, &resolved)?;
        Ok(Self {
            incident,
            exit,
            layers: layers.to_vec(),
            prepared,
        })
    }

    pub const fn incident_response(&self) -> DielectricResponse<'a> {
        self.incident
    }
    pub const fn exit_response(&self) -> DielectricResponse<'a> {
        self.exit
    }
    pub fn layer_count(&self) -> usize {
        self.layers.len()
    }
    pub fn layer(&self, index: usize) -> Option<DielectricLayer<'a>> {
        self.layers.get(index).copied()
    }
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.prepared.input.omega
    }
    pub fn transverse_wavenumber_rad_m(&self) -> f64 {
        self.prepared.input.kx
    }
    pub fn incident_permittivity_f_m(&self) -> C64 {
        self.prepared.input.epsilon_i
    }
    pub fn exit_permittivity_f_m(&self) -> C64 {
        self.prepared.input.epsilon_t
    }

    /// Actual peak Ex (TM), Ey (TE) [V/m] at z=0. Relative phase is preserved
    /// in all six field components. The powers include incident/reflected
    /// interference; normalizing them into additive R/T fractions is invalid
    /// in general. Zero illumination gives zero fields and heat.
    pub fn illuminate(
        &self,
        incident_tangential_e_v_m: [C64; 2],
    ) -> Result<AbsorbingLayerIllumination<'_, 'a>, ReflectionError> {
        Ok(AbsorbingLayerIllumination {
            source: self,
            electric: incident_tangential_e_v_m,
            solved: self.prepared.illuminate(incident_tangential_e_v_m)?,
        })
    }
}

impl<'s, 'a> AbsorbingLayerIllumination<'s, 'a> {
    pub fn source(&self) -> &'s PlanarAbsorbingDielectricLayers<'a> {
        self.source
    }
    pub fn response(&self) -> &AbsorbingLayerReflection {
        &self.solved.response
    }
    /// Integrated heat in each original finite film [W/m²], in source order.
    /// A thermal consumer supplies area and duration and deposits each term once.
    pub fn absorbed_power_per_layer_w_m2(&self) -> &[f64] {
        &self.solved.absorbed[..self.source.layer_count()]
    }
    /// Depth [m] from the selected film's first face, inclusive of both faces.
    /// Full E/H and the local dielectric heating include longitudinal TM E.
    pub fn layer_field(
        &self,
        layer_index: usize,
        depth_m: f64,
    ) -> Result<AbsorbingLayerPoint, ReflectionError> {
        self.source
            .prepared
            .field(self.electric, layer_index, depth_m)
    }
}
