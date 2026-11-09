//! Original dissipative solid coatings ending on a passive mechanical backing.
//!
//! A lossless incident fluid occupies z<0. Perfectly bonded original isotropic
//! solids follow in their authored order, and the backing acts at the last
//! face. The native exp(+i omega t) convention retains positive loss moduli,
//! inward-decaying waves and longitudinal strain. Compression p=-sigma_zz and
//! velocity toward +z obey p=Z*v at the backing. An impedance authored for
//! exp(-i omega t), such as a mechanical pin snapshot, must be CONJUGATED by its
//! source adapter. A single complex impedance does not supply a causal law.
//!
//! The original scalar bulk/shear Maxwell sources own every film's propagation,
//! strain ceiling and positive local heat. Backing resistance owns separate
//! positive work at its actual rear velocity. No physical rear fluid, averaged
//! attenuation law, flux-deficit heat, or fictitious interface loss is added.
//! This rung is normal-incidence longitudinal motion; a fluid cannot supply a
//! shear incident wave. It introduces no oblique conversion, slip, anisotropy,
//! finite aperture, evolving temperature or experimental source qualification.

use super::super::{ViscoelasticSlabError, ViscoelasticSlabField, kernel, poll};
use super::{
    ViscoelasticLayer, ViscoelasticLayerAbsorption, ViscoelasticLayerBudget, error,
    prepare_materials, prepare_stack,
};
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::{c64::C64, dd::Dd};
use kernel::layers::supported::{Termination, traveling_power};

type Result<T> = core::result::Result<T, ViscoelasticSlabError>;

/// Rear normal boundary in the original exp(+i omega t) convention.
/// Impedances are specific pressure/velocity ratios [Pa s/m], evaluated by the
/// caller's original support law at the stack's angular frequency.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ViscoelasticBacking {
    /// Exactly zero compressive traction, with unrestricted normal motion.
    TractionFree,
    /// Exactly zero normal velocity, with unrestricted reaction traction.
    Clamped,
    /// Finite passive p=Z*v; Re(Z)>=0. Reactance of either sign is admitted.
    Impedance(C64),
}

impl ViscoelasticBacking {
    fn impedance(self) -> Option<C64> {
        match self {
            Self::TractionFree => Some(C64::ZERO),
            Self::Clamped => None,
            Self::Impedance(value) => Some(value),
        }
    }
}

/// Independent physical powers over the common positive illuminated area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SupportedViscoelasticPower {
    pub incident_power_w: f64,
    pub reflected_power_w: f64,
    /// Positive original bulk constitutive work in the physical layers [W].
    pub bulk_absorbed_power_w: f64,
    /// Positive original shear constitutive work in the physical layers [W].
    pub shear_absorbed_power_w: f64,
    /// Positive Re(Z)*|v_rear|²*area/2 [W], separate from all coating heat.
    pub backing_absorbed_power_w: f64,
    /// Incident minus reflection, intrinsic heat and backing work [W].
    pub balance_residual_w: f64,
}

/// Prepared original stack with one physical incident fluid and one backing.
/// An empty stack is the direct fluid/backing boundary; zero-thickness layers
/// still validate their original material support but supply no field or heat.
pub struct SupportedViscoelasticLayers<'a> {
    layers: Vec<ViscoelasticLayer<'a>>,
    incident_source: &'a dyn FluidAcoustics,
    incident_state: FluidAcousticState,
    backing: ViscoelasticBacking,
    omega: f64,
    wave: Option<kernel::layers::Stack>,
    termination: Termination,
    last_interior: Option<usize>,
}

/// Admitted reflection and full original interior response to one incident wave.
/// Repeated observations allocate nothing and never re-query any source law.
pub struct DrivenSupportedViscoelasticLayers<'s, 'a> {
    prepared: &'s SupportedViscoelasticLayers<'a>,
    input: [C64; 2],
    area: f64,
    reflected: C64,
    support_pressure: C64,
    support_velocity: C64,
    power: SupportedViscoelasticPower,
}

impl<'a> SupportedViscoelasticLayers<'a> {
    /// Prepare original sources once, with the same layer/branch caps and
    /// cancellation as the native two-fluid stack. Internal and rear reference
    /// impedances are mathematical wave coordinates at the same physical face.
    /// No reference-port fluid source, material or absorbed power is exposed.
    pub fn try_new(
        layers: &[ViscoelasticLayer<'a>],
        incident_source: &'a dyn FluidAcoustics,
        backing: ViscoelasticBacking,
        angular_frequency_rad_s: f64,
        budget: ViscoelasticLayerBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        if !angular_frequency_rad_s.is_finite() || angular_frequency_rad_s <= 0.0 {
            return Err(ViscoelasticSlabError::OutsideSourceSupport);
        }
        if backing
            .impedance()
            .is_some_and(|z| !kernel::finite(z) || z.re < 0.0)
        {
            return Err(ViscoelasticSlabError::Invalid(
                "backing impedance must be finite and passive",
            ));
        }
        let material = prepare_materials(layers, angular_frequency_rad_s, budget, cancelled)?;
        let incident_state = incident_source.acoustic_state()?;
        let left = incident_state.characteristic_impedance_pa_s_m();
        let right = material.reference.last().copied().unwrap_or(left);
        let wave = if layers.is_empty() {
            None
        } else {
            Some(prepare_stack(
                layers,
                &material,
                angular_frequency_rad_s,
                [left, right],
                cancelled,
            )?)
        };
        let scattering = wave
            .as_ref()
            .map_or([C64::ZERO, C64::ONE, C64::ONE, C64::ZERO], |value| {
                value.scattering
            });
        let termination = Termination::new(scattering, right, backing.impedance())
            .map_err(ViscoelasticSlabError::Invalid)?;
        let last_interior = layers.iter().rposition(|layer| layer.thickness_m > 0.0);
        poll(cancelled)?;
        Ok(Self {
            layers: layers.to_vec(),
            incident_source,
            incident_state,
            backing,
            omega: angular_frequency_rad_s,
            wave,
            termination,
            last_interior,
        })
    }

    pub fn layer_count(&self) -> usize {
        self.layers.len()
    }
    pub fn layer(&self, index: usize) -> Option<ViscoelasticLayer<'a>> {
        self.layers.get(index).copied()
    }
    pub fn layer_moduli_pa(&self, index: usize) -> Option<[C64; 2]> {
        self.wave.as_ref()?.layer(index).map(|value| value.moduli)
    }
    pub const fn incident_source(&self) -> &'a dyn FluidAcoustics {
        self.incident_source
    }
    pub const fn incident_state(&self) -> FluidAcousticState {
        self.incident_state
    }
    pub const fn backing(&self) -> ViscoelasticBacking {
        self.backing
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    /// Original reflected pressure divided by incident pressure, at z=0.
    pub const fn pressure_reflection_coefficient(&self) -> C64 {
        self.termination.reflected
    }

    /// Drive the same positive area with one peak incident pressure phasor.
    /// Every original layer's whole-depth strain ceiling remains binding. Work
    /// is O(layers), with no drive allocation or constitutive/source query.
    pub fn drive<'s>(
        &'s self,
        incident_peak_pressure_pa: C64,
        area_m2: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DrivenSupportedViscoelasticLayers<'s, 'a>> {
        poll(cancelled)?;
        if !kernel::finite(incident_peak_pressure_pa)
            || !incident_peak_pressure_pa.abs().is_finite()
            || !area_m2.is_finite()
            || area_m2 <= 0.0
        {
            return Err(ViscoelasticSlabError::Invalid(
                "supported drive needs finite pressure and positive area",
            ));
        }
        let input = [
            incident_peak_pressure_pa,
            self.termination.incoming_right * incident_peak_pressure_pa,
        ];
        let support_pressure = self.termination.pressure * incident_peak_pressure_pa;
        let support_velocity = self.termination.velocity * incident_peak_pressure_pa;
        if input
            .into_iter()
            .chain([support_pressure, support_velocity])
            .any(|v| !kernel::finite(v) || !v.abs().is_finite())
        {
            return Err(ViscoelasticSlabError::Invalid(
                "supported boundary field is unrepresentable",
            ));
        }
        let (incident, reflected, reflected_power, heat) = if let Some(wave) = &self.wave {
            let power = wave.drive(input, area_m2, cancelled).map_err(error)?;
            (
                power.incident[0],
                power.outgoing_pressure[0],
                power.outgoing[0],
                power.heat,
            )
        } else {
            let reflected = self.termination.reflected * incident_peak_pressure_pa;
            let impedance = self.incident_state.characteristic_impedance_pa_s_m();
            (
                traveling_power(incident_peak_pressure_pa, impedance, area_m2)
                    .map_err(ViscoelasticSlabError::Invalid)?,
                reflected,
                traveling_power(reflected, impedance, area_m2)
                    .map_err(ViscoelasticSlabError::Invalid)?,
                [0.0; 2],
            )
        };
        let backing_heat = self
            .termination
            .absorbed_power(support_velocity, area_m2)
            .map_err(ViscoelasticSlabError::Invalid)?;
        let sum = [reflected_power, heat[0], heat[1], backing_heat]
            .into_iter()
            .fold(Dd::ZERO, |sum, value| sum + Dd::from_f64(value));
        let residual = (Dd::from_f64(incident) - sum).to_f64();
        let scale = (Dd::from_f64(incident) + sum).to_f64();
        if !kernel::finite(reflected)
            || !scale.is_finite()
            || !residual.is_finite()
            || residual.abs()
                > 8192.0
                    * (self.layers.len() + 1) as f64
                    * f64::EPSILON
                    * scale.max(f64::MIN_POSITIVE)
        {
            return Err(ViscoelasticSlabError::Invalid(
                "reflection, intrinsic work and backing work disagree",
            ));
        }
        poll(cancelled)?;
        Ok(DrivenSupportedViscoelasticLayers {
            prepared: self,
            input,
            area: area_m2,
            reflected,
            support_pressure,
            support_velocity,
            power: SupportedViscoelasticPower {
                incident_power_w: incident,
                reflected_power_w: reflected_power,
                bulk_absorbed_power_w: heat[0],
                shear_absorbed_power_w: heat[1],
                backing_absorbed_power_w: backing_heat,
                balance_residual_w: residual,
            },
        })
    }
}

impl<'s, 'a> DrivenSupportedViscoelasticLayers<'s, 'a> {
    pub const fn prepared(&self) -> &'s SupportedViscoelasticLayers<'a> {
        self.prepared
    }
    pub const fn area_m2(&self) -> f64 {
        self.area
    }
    pub const fn power(&self) -> SupportedViscoelasticPower {
        self.power
    }
    pub const fn reflected_pressure_pa(&self) -> C64 {
        self.reflected
    }
    /// Actual rear compression p=-sigma_zz, acting into the backing [Pa].
    pub const fn support_compressive_traction_pa(&self) -> C64 {
        self.support_pressure
    }
    /// Actual rear velocity toward +z [m/s], retained even for a very stiff
    /// finite support. It is not a subtraction of rounded traveling waves.
    pub const fn support_normal_velocity_m_s(&self) -> C64 {
        self.support_velocity
    }

    pub fn layer_absorption(
        &self,
        index: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ViscoelasticLayerAbsorption> {
        poll(cancelled)?;
        let wave = self
            .prepared
            .wave
            .as_ref()
            .ok_or(ViscoelasticSlabError::Invalid(
                "empty coating has no physical layer",
            ))?;
        let local = wave
            .layer_drive(index, self.input, self.area)
            .map_err(error)?;
        poll(cancelled)?;
        Ok(ViscoelasticLayerAbsorption {
            source_identity: self.prepared.layers[index].source.source_identity,
            bulk_absorbed_power_w: local.heat[0],
            shear_absorbed_power_w: local.heat[1],
            maximum_absolute_strain_bound: local.strain_bound,
        })
    }

    /// Direct positive original bulk/shear work on an actual layer interval.
    /// No subtraction of cumulative heats or interface fluxes supplies it.
    pub fn absorption_between(
        &self,
        index: usize,
        depths_m: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ViscoelasticLayerAbsorption> {
        poll(cancelled)?;
        let wave = self
            .prepared
            .wave
            .as_ref()
            .ok_or(ViscoelasticSlabError::Invalid(
                "empty coating has no physical layer",
            ))?;
        let local = wave
            .layer_drive(index, self.input, self.area)
            .map_err(error)?;
        let heat = wave
            .layer(index)
            .expect("validated coating layer")
            .interval_heat(local.inward, local.parity, depths_m, self.area)
            .map_err(ViscoelasticSlabError::Invalid)?;
        poll(cancelled)?;
        Ok(ViscoelasticLayerAbsorption {
            source_identity: self.prepared.layers[index].source.source_identity,
            bulk_absorbed_power_w: heat[0],
            shear_absorbed_power_w: heat[1],
            maximum_absolute_strain_bound: local.strain_bound,
        })
    }

    /// The original normal field at depth from this layer's first face.
    /// At the physical rear face the directly retained support velocity and
    /// pressure also feed the same original constitutive strain/heat owner.
    pub fn field(
        &self,
        index: usize,
        depth_m: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ViscoelasticSlabField> {
        poll(cancelled)?;
        let wave = self
            .prepared
            .wave
            .as_ref()
            .ok_or(ViscoelasticSlabError::Invalid(
                "empty coating has no physical layer",
            ))?;
        let layer = wave.layer(index).ok_or(ViscoelasticSlabError::Invalid(
            "physical layer index is outside the coating",
        ))?;
        let field = if self.prepared.last_interior == Some(index) && depth_m == layer.length {
            layer
                .field_from_trace(self.support_pressure, self.support_velocity)
                .map_err(ViscoelasticSlabError::Invalid)?
        } else {
            wave.field(index, self.input, depth_m).map_err(error)?
        };
        poll(cancelled)?;
        Ok(ViscoelasticSlabField {
            compressive_traction_pa: field.pressure,
            normal_velocity_m_s: field.velocity,
            longitudinal_strain: field.strain,
            bulk_heat_density_w_m3: field.heat_density[0],
            shear_heat_density_w_m3: field.heat_density[1],
        })
    }
}

#[cfg(test)]
mod tests;
