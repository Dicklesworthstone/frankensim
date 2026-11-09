//! A resolved elastic coating backed by original oriented support mechanics.

use super::oriented::support_motion;
use super::{ElasticSupportError, OrientedOscillatorSupport};
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::elastic::acoustic::interface::ElasticMode;
use fs_material::elastic::acoustic::interface::layer::supported::{
    ElasticCoatingResponse, IsotropicElasticCoating, PreparedElasticCoating,
};
use fs_material::elastic::acoustic::interface::surface::ElasticSurfaceCondition;
use fs_math::c64::C64;

/// Coherent reflection through an original finite solid coating into three
/// original principal mass/spring/damper supports.
///
/// The incident solid occupies z < 0, the coating 0 <= z <= h, and the support
/// acts on the rear +z face. The original material tensors govern all P/SV/SH
/// propagation, conversion, interference and interior strain. Only the backing
/// is locally reacting; the coating is a resolved slab, including evanescence.
///
/// Both solids are lossless and explicitly isotropic. Support resistance is the
/// sole dissipative account; no coating loss is inferred from standing waves.
/// The original stationary, linear pin model and common principal frame apply.
#[derive(Debug)]
pub struct OscillatorBackedElasticCoating<'a> {
    support: OrientedOscillatorSupport,
    omega_rad_s: f64,
    rotation: [f64; 9],
    coating: PreparedElasticCoating<'a>,
}

/// Resolved coating fields and mechanical storage at the actual rear support.
#[derive(Debug)]
pub struct SupportedElasticCoatingResponse<'p, 'a> {
    /// Native reflection, interior fields, material strain domains and positive
    /// support work. Interior fields remain queryable through `field_at`.
    pub coating: ElasticCoatingResponse<'p, 'a>,
    /// Peak velocity in each original principal support direction [m/s].
    pub principal_velocity_m_s: [C64; 3],
    /// Peak displacement in each original principal support direction [m].
    pub principal_displacement_m: [C64; 3],
    /// Peak displacement of the REAR support in the surface frame [m].
    pub support_displacement_m: [C64; 3],
    /// Cycle-mean kinetic storage of each original pin [J/m²].
    pub mean_kinetic_energy_j_m2: [f64; 3],
    /// Cycle-mean spring storage of each original pin [J/m²].
    pub mean_elastic_energy_j_m2: [f64; 3],
    /// Total pin storage [J/m²], excluding the resolved coating's bulk energy.
    pub mean_support_stored_energy_j_m2: f64,
}

impl<'a> OscillatorBackedElasticCoating<'a> {
    /// Prepare one frequency in Hz. The support angular frequency is derived
    /// once using the SAME `TAU * frequency_hz` as the native layer propagator.
    /// At zero thickness the original direct supported surface is recovered.
    pub fn prepare(
        incident: &'a ElasticAcousticState,
        coating: &'a ElasticAcousticState,
        thickness_m: f64,
        support: OrientedOscillatorSupport,
        frequency_hz: f64,
        mode: ElasticMode,
        incidence_cosine: f64,
    ) -> Result<Self, ElasticSupportError> {
        Self::prepare_controlled(
            incident,
            coating,
            thickness_m,
            support,
            frequency_hz,
            mode,
            incidence_cosine,
            &mut || false,
        )
    }

    /// Bounded, transactional native preparation. It retains the full coating
    /// state rather than replacing the slab by a scalar effective impedance.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_controlled(
        incident: &'a ElasticAcousticState,
        coating: &'a ElasticAcousticState,
        thickness_m: f64,
        support: OrientedOscillatorSupport,
        frequency_hz: f64,
        mode: ElasticMode,
        incidence_cosine: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ElasticSupportError> {
        if cancelled() {
            return Err(ElasticSupportError::Cancelled);
        }
        let omega_rad_s = core::f64::consts::TAU * frequency_hz;
        let impedance = support.impedance_at_frequency(omega_rad_s)?;
        let coating = IsotropicElasticCoating::new(
            incident,
            coating,
            thickness_m,
            ElasticSurfaceCondition::TensorImpedance(impedance),
        )?
        .prepare_incidence_cancellable(mode, incidence_cosine, frequency_hz, cancelled)?;
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

    pub const fn coating(&self) -> &PreparedElasticCoating<'a> {
        &self.coating
    }

    /// Apply peak incident P/SV/SH velocities [m/s]. Prepared field maps and
    /// the shared original pin-storage evaluator are reused without allocation
    /// or factorization. Front-face motion is never substituted for rear motion.
    pub fn drive(
        &self,
        incoming_velocity_m_s: [C64; 3],
    ) -> Result<SupportedElasticCoatingResponse<'_, 'a>, ElasticSupportError> {
        let coating = self.coating.drive(incoming_velocity_m_s)?;
        let motion = support_motion(
            *self.support.principal_supports(),
            self.rotation,
            coating.support_velocity_m_s,
            self.omega_rad_s,
        )?;
        Ok(SupportedElasticCoatingResponse {
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
