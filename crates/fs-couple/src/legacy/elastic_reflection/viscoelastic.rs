//! Original dissipative layers backed by the original normal wall oscillator.
//!
//! The native viscoelastic stack uses exp(+i omega t), pressure=-sigma_zz and
//! positive velocity into the backing. The original wall impedance uses the
//! opposite phasor sign, so its conjugate is supplied explicitly. The same
//! phase conversion gives x=v/(i omega) for the returned backing displacement.

mod heating;
pub use heating::OscillatorBackedViscoelasticHeating;

use super::{ElasticSupportError, displacement, mean_energy};
use fs_material::elastic::acoustic::viscoelastic::layers::{
    DrivenSupportedViscoelasticLayers, SupportedViscoelasticLayers, ViscoelasticBacking,
    ViscoelasticLayer, ViscoelasticLayerBudget,
};
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_math::{c64::C64, dd::Dd};
use fs_phs::{WallPin, wall_specific_impedance};

/// Normal fluid-pressure reflection by original viscoelastic solids and a
/// physical mass/spring/damper backing. Original bulk/shear heat and backing
/// resistance remain separate; no effective material loss replaces either.
pub struct OscillatorBackedViscoelasticLayers<'a> {
    backing: WallPin,
    omega_rad_s: f64,
    layers: SupportedViscoelasticLayers<'a>,
}

/// Native resolved pressure/strain/heat fields and actual rear pin storage.
pub struct SupportedViscoelasticResponse<'s, 'a> {
    /// Original material-layer heat and backing heat are separate TOTAL watts.
    pub layers: DrivenSupportedViscoelasticLayers<'s, 'a>,
    /// Peak displacement under the native exp(+i omega t) convention [m].
    pub support_displacement_m: C64,
    pub mean_kinetic_energy_j_m2: f64,
    pub mean_elastic_energy_j_m2: f64,
    /// Cycle-mean pin storage per area; no coating bulk storage is included.
    pub mean_support_stored_energy_j_m2: f64,
    /// Same original pin storage over the drive's actual positive area [J].
    pub mean_support_stored_energy_j: f64,
}

impl<'a> OscillatorBackedViscoelasticLayers<'a> {
    pub fn prepare(
        layers: &[ViscoelasticLayer<'a>],
        incident_source: &'a dyn FluidAcoustics,
        backing: WallPin,
        omega_rad_s: f64,
        budget: ViscoelasticLayerBudget,
    ) -> Result<Self, ElasticSupportError> {
        Self::prepare_controlled(
            layers,
            incident_source,
            backing,
            omega_rad_s,
            budget,
            &mut || false,
        )
    }

    pub fn prepare_controlled(
        layers: &[ViscoelasticLayer<'a>],
        incident_source: &'a dyn FluidAcoustics,
        backing: WallPin,
        omega_rad_s: f64,
        budget: ViscoelasticLayerBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ElasticSupportError> {
        if cancelled() {
            return Err(ElasticSupportError::Cancelled);
        }
        let impedance = wall_specific_impedance(&backing, omega_rad_s)?.conj();
        let layers = SupportedViscoelasticLayers::try_new(
            layers,
            incident_source,
            ViscoelasticBacking::Impedance(impedance),
            omega_rad_s,
            budget,
            cancelled,
        )?;
        Ok(Self {
            backing,
            omega_rad_s,
            layers,
        })
    }

    pub const fn backing(&self) -> &WallPin {
        &self.backing
    }

    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    pub const fn layers(&self) -> &SupportedViscoelasticLayers<'a> {
        &self.layers
    }

    /// Reuse the native pressure solve and heat integrals, then evaluate the
    /// original pin potentials from actual rear motion. This advances no clock.
    pub fn drive(
        &self,
        incident_peak_pressure_pa: C64,
        area_m2: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SupportedViscoelasticResponse<'_, 'a>, ElasticSupportError> {
        let layers = self
            .layers
            .drive(incident_peak_pressure_pa, area_m2, cancelled)?;
        let velocity = layers.support_normal_velocity_m_s();
        // Reuse the original displacement evaluator with an explicit phasor
        // conjugation; silently copying its -i*omega convention would flip x.
        let support_displacement_m = displacement(velocity.conj(), self.omega_rad_s)?.conj();
        let mean_kinetic_energy_j_m2 = mean_energy(self.backing.surface_density, velocity.abs())?;
        let mean_elastic_energy_j_m2 = mean_energy(
            self.backing.stiffness_per_area,
            support_displacement_m.abs(),
        )?;
        let mean_support_stored_energy_j_m2 = (Dd::from_f64(mean_kinetic_energy_j_m2)
            + Dd::from_f64(mean_elastic_energy_j_m2))
        .to_f64();
        let mean_support_stored_energy_j = mean_support_stored_energy_j_m2 * area_m2;
        if !mean_support_stored_energy_j_m2.is_finite()
            || !mean_support_stored_energy_j.is_finite()
            || (mean_support_stored_energy_j_m2 > 0.0 && mean_support_stored_energy_j == 0.0)
        {
            return Err(ElasticSupportError::NumericalRange);
        }
        if cancelled() {
            return Err(ElasticSupportError::Cancelled);
        }
        Ok(SupportedViscoelasticResponse {
            layers,
            support_displacement_m,
            mean_kinetic_energy_j_m2,
            mean_elastic_energy_j_m2,
            mean_support_stored_energy_j_m2,
            mean_support_stored_energy_j,
        })
    }
}
