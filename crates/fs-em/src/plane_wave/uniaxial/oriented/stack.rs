//! Coherent transparent crystal stacks with independent fixed global axes.
//!
//! Each finite layer retains its original principal-index laws and scalar
//! permeability. Native Maxwell modes are matched at every parallel face.
//! A right-to-left reflection/Schur sweep and forward field recovery use
//! only decaying-side propagation, with O(layers) preparation and storage.
//! The absorbing extension retains these same original sources and placements;
//! see [`super::absorbing::stack`]. No effective medium or causal reservoir is inferred.

use super::{DielectricResponse, Incidence, ReflectionError, UniaxialOpticalResponse};
use crate::plane_wave::angular_spectrum::AngularSpectrumError;
use fs_material::field_storage::linear_magnetic_permeability;
use fs_math::{c64::C64, det};

mod kernel;
pub use kernel::{StackRegion, UniaxialLayerField, UniaxialStackResponse};
mod spectrum;
pub use spectrum::PreparedUniaxialStackSpectrum;

const C0: f64 = 299_792_458.0;
pub const MAX_UNIAXIAL_STACK_LAYERS: usize = 256;
type Result<T> = core::result::Result<T, AngularSpectrumError>;

/// One homogeneous original crystal and its authored placement in the stack.
#[derive(Clone, Copy)]
pub struct UniaxialStackLayer<'m> {
    crystal: UniaxialOpticalResponse<'m>,
    optic_axis: [f64; 3],
    thickness_m: f64,
}

impl<'m> UniaxialStackLayer<'m> {
    pub fn try_new(
        crystal: UniaxialOpticalResponse<'m>,
        optic_axis: [f64; 3],
        thickness_m: f64,
    ) -> core::result::Result<Self, ReflectionError> {
        Ok(Self {
            crystal,
            optic_axis: super::kernel::unit(optic_axis).map_err(ReflectionError::Invalid)?,
            thickness_m: super::positive(thickness_m)?,
        })
    }
    pub const fn crystal_response(&self) -> UniaxialOpticalResponse<'m> { self.crystal }
    pub const fn optic_axis(&self) -> [f64; 3] { self.optic_axis }
    pub const fn thickness_m(&self) -> f64 { self.thickness_m }
}

/// Checked before allocating frequency data or solving any Maxwell modes.
#[derive(Clone, Copy, Debug)]
pub struct UniaxialStackLimits {
    pub maximum_layers: usize,
    /// Peak logical buffers for ONE plane-wave response, including resolved
    /// source data, elimination workspace, retained fields and accumulator.
    /// Caller-owned source data, fixed stack temporaries and allocator overhead
    /// are excluded. Prepared spectra additionally bound their aggregate storage.
    pub maximum_workspace_bytes: usize,
}

/// Original entrance/exit half-spaces and an ordered borrowed crystal slice.
/// Layer zero starts at z=0; all layer normals point along global +z.
pub struct OrientedUniaxialStack<'s, 'm> {
    incident: DielectricResponse<'m>,
    exit: DielectricResponse<'m>,
    layers: &'s [UniaxialStackLayer<'m>],
    incident_mu: f64,
    exit_mu: f64,
    thickness_m: f64,
    limits: UniaxialStackLimits,
}

impl<'s, 'm> OrientedUniaxialStack<'s, 'm> {
    pub fn try_new(
        incident: DielectricResponse<'m>,
        exit: DielectricResponse<'m>,
        layers: &'s [UniaxialStackLayer<'m>],
        limits: UniaxialStackLimits,
    ) -> Result<Self> {
        if layers.is_empty() || layers.len() > MAX_UNIAXIAL_STACK_LAYERS
            || layers.len() > limits.maximum_layers
            || kernel::workspace_bytes(layers.len(), 1)? > limits.maximum_workspace_bytes {
            return Err(bad("uniaxial stack exceeds its layer/workspace budget"));
        }
        let mut depth = 0.0;
        for layer in layers {
            let next = depth + layer.thickness_m;
            if !next.is_finite() || next <= depth {
                return Err(bad("uniaxial stack face positions are unrepresentable"));
            }
            depth = next;
        }
        Ok(Self {
            incident, exit, layers,
            incident_mu: linear_magnetic_permeability(incident.material()).map_err(ReflectionError::from)?,
            exit_mu: linear_magnetic_permeability(exit.material()).map_err(ReflectionError::from)?,
            thickness_m: depth, limits,
        })
    }
    pub const fn incident_response(&self) -> DielectricResponse<'m> { self.incident }
    pub const fn exit_response(&self) -> DielectricResponse<'m> { self.exit }
    pub const fn layers(&self) -> &'s [UniaxialStackLayer<'m>] { self.layers }
    pub const fn thickness_m(&self) -> f64 { self.thickness_m }

    pub(super) const fn preparation_limits(&self) -> UniaxialStackLimits { self.limits }

    pub(super) const fn exterior_permeabilities_h_m(&self) -> [f64; 2] {
        [self.incident_mu, self.exit_mu]
    }

    /// Source-supported positive frequency and original complex tangential
    /// Ex/Ey [V/m]. All internal reflections and polarization conversion remain
    /// coherent. A cancellation or failed gate returns no partial snapshot.
    pub fn response(
        &self, omega: f64, cosine: f64, electric: [C64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<UniaxialStackResponse> {
        poll(&mut cancelled)?;
        let data = self.resolve_frequency(omega, &mut cancelled)?;
        kernel::solve(&data, Incidence::Cosine(cosine), electric, [1.0,0.0], &mut cancelled)
    }

    /// Signed conserved kx/k0; the original entrance must propagate.
    pub fn response_at_tangential_index(
        &self, omega: f64, index: f64, electric: [C64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<UniaxialStackResponse> {
        poll(&mut cancelled)?;
        let data = self.resolve_frequency(omega, &mut cancelled)?;
        kernel::solve(&data, Incidence::TangentialIndex(index), electric, [1.0,0.0], &mut cancelled)
    }

    fn resolve_frequency(&self, omega: f64, cancelled: &mut impl FnMut() -> bool) -> Result<Frequency> {
        let omega = super::positive(omega)?;
        let index = |law: DielectricResponse<'m>, mu| -> Result<f64> {
            Ok(super::positive((C0 * det::sqrt(mu)) * det::sqrt(law.lossless_permittivity(omega)?))?)
        };
        let ni = index(self.incident, self.incident_mu)?;
        let nt = index(self.exit, self.exit_mu)?;
        let mut layers = Vec::new();
        layers.try_reserve_exact(self.layers.len()).map_err(|_| bad("uniaxial stack allocation failed"))?;
        let mut entrance = 0.0;
        for source in self.layers {
            poll(cancelled)?;
            let [no, ne] = source.crystal.phase_indices(omega)?;
            let exit = entrance + source.thickness_m;
            layers.push(ResolvedLayer {
                no, ne, mu: source.crystal.permeability_h_m(),
                axis: source.optic_axis, thickness: source.thickness_m, entrance, exit,
            });
            entrance = exit;
        }
        Ok(Frequency { omega, ni, nt, mui: self.incident_mu, mut_: self.exit_mu, layers })
    }
}

#[derive(Clone, Copy)]
struct ResolvedLayer {
    no: f64, ne: f64, mu: f64, axis: [f64; 3],
    thickness: f64, entrance: f64, exit: f64,
}
struct Frequency {
    omega: f64, ni: f64, nt: f64, mui: f64, mut_: f64,
    layers: Vec<ResolvedLayer>,
}
fn bad(message: &'static str) -> AngularSpectrumError { AngularSpectrumError::Invalid(message) }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() { Err(AngularSpectrumError::Cancelled) } else { Ok(()) }
}
