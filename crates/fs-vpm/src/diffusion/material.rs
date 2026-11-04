//! Material viscosity for the existing particle-strength-exchange operator.
//!
//! One original [`FluidTransport`] evaluation supplies the checked kinematic
//! viscosity `nu = mu / rho` in SI. The binding borrows the source and retains
//! its complete transport snapshot; it neither copies source authority into
//! a raw coefficient nor invents a new material model.

use super::{PseDiffusion, PseKernel};
use crate::VpmError;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

/// An original material refusal or an invalid native diffusion length scale.
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialPseError {
    /// The original source could not supply a coherent transport state.
    Material(MaterialError),
    /// The original PSE owner rejected the resulting numerical inputs.
    Diffusion(VpmError),
}

impl core::fmt::Display for MaterialPseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "PSE material transport: {error}"),
            Self::Diffusion(error) => write!(f, "material PSE diffusion: {error}"),
        }
    }
}

impl std::error::Error for MaterialPseError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Material(error) => Some(error),
            Self::Diffusion(error) => Some(error),
        }
    }
}

/// Source-bound homogeneous Newtonian diffusion at one evaluated state point.
///
/// Use [`Self::law`] directly with the existing checked diffusion, viscous
/// advection, remeshing or prepared-graph entry points. Their kernels, budgets,
/// conservation and cancellation behavior remain the original owners'. No
/// source callback or extra allocation is added to any pair sweep.
///
/// Density and dynamic viscosity are frozen in this binding. Prepare another
/// binding from the next evaluated source state when temperature, pressure or
/// composition changes. A snapshot alone does not establish that one constant
/// viscosity is valid throughout an evolution. Spatially varying viscosity,
/// non-Newtonian vorticity transport, thermal feedback and phase changes are
/// not supplied. The borrowed source retains its original receipts and domain;
/// a raw [`FluidTransportState`] remains raw, without material-card authority.
#[derive(Debug)]
pub struct MaterialPseDiffusion<'m, F: FluidTransport + ?Sized> {
    source: &'m F,
    transport: FluidTransportState,
    law: PseDiffusion,
}

impl<'m, F: FluidTransport + ?Sized> MaterialPseDiffusion<'m, F> {
    /// Bind the original infinite-support Gaussian kernel at width [m].
    /// The source is evaluated exactly once; its checked `mu/rho` is used as-is.
    ///
    /// # Errors
    /// Retains the original material refusal or native PSE input refusal.
    pub fn new(source: &'m F, kernel_width_m: f64) -> Result<Self, MaterialPseError> {
        Self::prepare(source, kernel_width_m, PseKernel::Gaussian)
    }

    /// Bind the original compact Wendland kernel at support radius [m].
    /// This radius is independent of the Biot–Savart regularization core.
    ///
    /// # Errors
    /// Retains the original material refusal or native PSE input refusal.
    pub fn compact_wendland(
        source: &'m F,
        support_radius_m: f64,
    ) -> Result<Self, MaterialPseError> {
        Self::prepare(source, support_radius_m, PseKernel::CompactWendland)
    }

    fn prepare(source: &'m F, width: f64, kernel: PseKernel) -> Result<Self, MaterialPseError> {
        let transport = source
            .transport_state()
            .map_err(MaterialPseError::Material)?;
        let law = PseDiffusion::with_kernel(transport.kinematic_viscosity_m2_s(), width, kernel)
            .map_err(MaterialPseError::Diffusion)?;
        Ok(Self {
            source,
            transport,
            law,
        })
    }

    /// Original source owner, including any receipts and applicability limits.
    #[must_use]
    pub const fn source(&self) -> &'m F {
        self.source
    }

    /// Immutable coherent transport snapshot actually used for this law.
    #[must_use]
    pub const fn transport(&self) -> &FluidTransportState {
        &self.transport
    }

    /// Native constant-viscosity law for every existing PSE consumer.
    /// Keep this binding alongside the run when retaining source provenance.
    #[must_use]
    pub const fn law(&self) -> PseDiffusion {
        self.law
    }
}
