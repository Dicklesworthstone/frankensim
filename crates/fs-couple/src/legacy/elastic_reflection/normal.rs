//! Original oriented pin mechanics on source-bound anisotropic wave modes.

mod coating;
pub use coating::{NormalOscillatorBackedElasticCoating, NormalSupportedElasticCoatingResponse};

use super::oriented::support_motion;
use super::{ElasticSupportError, OrientedOscillatorSupport};
use fs_material::elastic::acoustic::interface::normal::{
    NormalElasticCharacteristic, NormalElasticSurfaceResponse, PreparedNormalElasticSurface,
};
use fs_material::elastic::acoustic::interface::surface::ElasticSurfaceCondition;
use fs_math::c64::C64;

/// Coherent normal-incidence reflection from the original anisotropic solid
/// and the original principal support masses, springs and resistances.
///
/// Support directions use the characteristic's physical surface frame. Native
/// channels keep ascending Christoffel speed and the same polarization in both
/// directions. A quasi-mode is not relabeled as an isotropic P or shear wave.
#[derive(Debug)]
pub struct NormalOscillatorSupportedSurface<'a> {
    support: OrientedOscillatorSupport,
    omega_rad_s: f64,
    rotation: [f64; 9],
    reflection: PreparedNormalElasticSurface<'a>,
}

/// Native quasi-mode reflection and storage in the actual support directions.
#[derive(Debug)]
pub struct NormalSupportedElasticResponse {
    /// Original material fields, strain admission and sole positive support work.
    pub reflection: NormalElasticSurfaceResponse,
    pub principal_velocity_m_s: [C64; 3],
    pub principal_displacement_m: [C64; 3],
    /// Peak support displacement in the characteristic's surface frame [m].
    pub support_displacement_m: [C64; 3],
    /// Original cycle-mean pin kinetic and spring energies [J/m²].
    pub mean_kinetic_energy_j_m2: [f64; 3],
    pub mean_elastic_energy_j_m2: [f64; 3],
    pub mean_stored_energy_j_m2: f64,
}

impl<'a> NormalOscillatorSupportedSurface<'a> {
    pub fn prepare(
        characteristic: &NormalElasticCharacteristic<'a>,
        support: OrientedOscillatorSupport,
        omega_rad_s: f64,
    ) -> Result<Self, ElasticSupportError> {
        Self::prepare_controlled(characteristic, support, omega_rad_s, &mut || false)
    }

    pub fn prepare_controlled(
        characteristic: &NormalElasticCharacteristic<'a>,
        support: OrientedOscillatorSupport,
        omega_rad_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ElasticSupportError> {
        if cancelled() {
            return Err(ElasticSupportError::Cancelled);
        }
        let impedance = support.impedance_at_frequency(omega_rad_s)?;
        let reflection = characteristic.prepare_surface_cancellable(
            ElasticSurfaceCondition::TensorImpedance(impedance),
            cancelled,
        )?;
        Ok(Self {
            support,
            omega_rad_s,
            rotation: support.principal_to_surface().matrix().m,
            reflection,
        })
    }

    pub const fn support(&self) -> &OrientedOscillatorSupport {
        &self.support
    }

    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    pub const fn reflection(&self) -> &PreparedNormalElasticSurface<'a> {
        &self.reflection
    }

    /// Drive original ascending-speed quasi-mode peak velocities [m/s]. The
    /// prepared native characteristic/support solve and shared pin potentials
    /// are reused without allocation, source re-evaluation or a new solve.
    pub fn drive(
        &self,
        incoming_velocity_m_s: [C64; 3],
    ) -> Result<NormalSupportedElasticResponse, ElasticSupportError> {
        let reflection = self.reflection.drive(incoming_velocity_m_s)?;
        let motion = support_motion(
            *self.support.principal_supports(),
            self.rotation,
            reflection.boundary_velocity_m_s,
            self.omega_rad_s,
        )?;
        Ok(NormalSupportedElasticResponse {
            reflection,
            principal_velocity_m_s: motion.principal_velocity,
            principal_displacement_m: motion.storage.displacement,
            support_displacement_m: motion.surface_displacement,
            mean_kinetic_energy_j_m2: motion.storage.kinetic,
            mean_elastic_energy_j_m2: motion.storage.elastic,
            mean_stored_energy_j_m2: motion.storage.total,
        })
    }
}
