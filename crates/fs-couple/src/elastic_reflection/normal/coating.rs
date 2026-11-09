//! Resolved anisotropic coating waves and original rear support mechanics.

use super::super::oriented::support_motion;
use super::super::{ElasticSupportError, OrientedOscillatorSupport};
use fs_material::elastic::acoustic::interface::normal::NormalElasticCharacteristic;
use fs_material::elastic::acoustic::interface::normal::coating::{
    NormalElasticCoating, NormalElasticCoatingResponse, PreparedNormalElasticCoating,
};
use fs_material::elastic::acoustic::interface::surface::ElasticSurfaceCondition;
use fs_math::c64::C64;

/// A finite original anisotropic solid coating on three original principal
/// mass/spring/damper supports, illuminated from an original elastic solid.
///
/// Both characteristics describe the same physical surface axes and positive
/// normal, each expressed in its own original material frame. Native modes
/// retain ascending speed and the same polarization in both directions.
/// The material owner resolves normal-wave propagation, interference, vector
/// conversion and both original strain domains; the support remains local.
///
/// Both solids are lossless. Original support resistance is the sole heat
/// account. This harmonic owner does not advance causal support memory or
/// infer bulk coating loss, oblique branches or finite-aperture transport.
#[derive(Debug)]
pub struct NormalOscillatorBackedElasticCoating<'a> {
    support: OrientedOscillatorSupport,
    omega_rad_s: f64,
    rotation: [f64; 9],
    coating: PreparedNormalElasticCoating<'a>,
}

/// Original coating fields plus storage at the actual rear mechanical support.
#[derive(Debug)]
pub struct NormalSupportedElasticCoatingResponse<'p, 'a> {
    /// Reflection and original interior fields, including material-frame strain.
    pub coating: NormalElasticCoatingResponse<'p, 'a>,
    pub principal_velocity_m_s: [C64; 3],
    pub principal_displacement_m: [C64; 3],
    /// Peak rear-support displacement in the common surface axes [m].
    pub support_displacement_m: [C64; 3],
    /// Original cycle-mean pin kinetic and spring storage [J/m²].
    pub mean_kinetic_energy_j_m2: [f64; 3],
    pub mean_elastic_energy_j_m2: [f64; 3],
    /// Total pin storage, excluding the coating's bulk elastic energy [J/m²].
    pub mean_support_stored_energy_j_m2: f64,
}

impl<'a> NormalOscillatorBackedElasticCoating<'a> {
    /// Bind one carrier frequency in Hz and one finite original coating.
    /// Wall mechanics and native propagation use the same TAU*frequency.
    /// Exactly zero thickness recovers the original direct supported surface.
    pub fn prepare(
        incident: &NormalElasticCharacteristic<'a>,
        coating: &NormalElasticCharacteristic<'a>,
        thickness_m: f64,
        support: OrientedOscillatorSupport,
        frequency_hz: f64,
    ) -> Result<Self, ElasticSupportError> {
        Self::prepare_controlled(
            incident,
            coating,
            thickness_m,
            support,
            frequency_hz,
            &mut || false,
        )
    }

    pub fn prepare_controlled(
        incident: &NormalElasticCharacteristic<'a>,
        coating: &NormalElasticCharacteristic<'a>,
        thickness_m: f64,
        support: OrientedOscillatorSupport,
        frequency_hz: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ElasticSupportError> {
        if cancelled() {
            return Err(ElasticSupportError::Cancelled);
        }
        let omega_rad_s = core::f64::consts::TAU * frequency_hz;
        let impedance = support.impedance_at_frequency(omega_rad_s)?;
        let coating = NormalElasticCoating::new(
            incident,
            coating,
            thickness_m,
            ElasticSurfaceCondition::TensorImpedance(impedance),
        )?
        .prepare_cancellable(frequency_hz, cancelled)?;
        Ok(Self {
            support,
            omega_rad_s,
            rotation: support.principal_to_surface().matrix().m,
            coating,
        })
    }

    pub const fn support(&self) -> &OrientedOscillatorSupport {
        &self.support
    }

    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    pub const fn coating(&self) -> &PreparedNormalElasticCoating<'a> {
        &self.coating
    }

    /// Apply incoming native quasi-mode peak velocities [m/s]. The original
    /// prepared wave maps and shared pin-storage evaluator are reused, with no
    /// allocation, factorization, material query or transient clock advance.
    pub fn drive(
        &self,
        incoming_velocity_m_s: [C64; 3],
    ) -> Result<NormalSupportedElasticCoatingResponse<'_, 'a>, ElasticSupportError> {
        self.with_storage(self.coating.drive(incoming_velocity_m_s)?)
    }

    /// Apply physical incoming velocity in the common surface Cartesian axes.
    pub fn drive_cartesian(
        &self,
        incoming_velocity_m_s: [C64; 3],
    ) -> Result<NormalSupportedElasticCoatingResponse<'_, 'a>, ElasticSupportError> {
        self.with_storage(self.coating.drive_cartesian(incoming_velocity_m_s)?)
    }

    fn with_storage<'p>(
        &'p self,
        coating: NormalElasticCoatingResponse<'p, 'a>,
    ) -> Result<NormalSupportedElasticCoatingResponse<'p, 'a>, ElasticSupportError> {
        let motion = support_motion(
            *self.support.principal_supports(),
            self.rotation,
            coating.support_velocity_m_s,
            self.omega_rad_s,
        )?;
        Ok(NormalSupportedElasticCoatingResponse {
            coating,
            principal_velocity_m_s: motion.principal_velocity,
            principal_displacement_m: motion.storage.displacement,
            support_displacement_m: motion.surface_displacement,
            mean_kinetic_energy_j_m2: motion.storage.kinetic,
            mean_elastic_energy_j_m2: motion.storage.elastic,
            mean_support_stored_energy_j_m2: motion.storage.total,
        })
    }
}
