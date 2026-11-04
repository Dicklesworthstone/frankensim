//! Coherent acoustic beams through the original planar layer solver.
//!
//! The caller supplies a finite angular sum of complex PEAK pressures at the
//! entrance plane. Each component uses the original angle-resolved R/T solve,
//! conserving both tangential wavenumbers. With `exp(+i omega t)`, a forward
//! component is `p = A exp(-i k dot x)` and momentum gives `v = k p/(rho omega)`.
//! The receiving normal wavenumber may be imaginary: the original outgoing
//! branch then decays with depth while retaining reactive normal intensity.
//!
//! Fields are superposed before forming `p conj(v)/2`. The incoherent power
//! average in [`super::angular`] cannot replace this interference calculation.
//! All components have the SAME original frequency and frozen layer/material
//! state. The transverse planes are infinite; angular discretization and any
//! compact-aperture approximation belong to the caller. No lateral edge
//! scattering, mean flow, exterior bulk loss, or angular error bound is added.

use super::{LayeredTransmission, RigidBackedLayers, TransmittingLayers, finite};
use crate::MaterialError;
use fs_math::{c64::C64, dd::Dd, det};

mod interior;
pub use super::elastic::stack::{
    PreparedElasticAcousticSpectrum, PreparedElasticLayerAcousticSpectrum,
};
pub(in crate::porous::layers) use interior::PreparedLayerFields;
pub use interior::{AcousticLayerSample, PreparedLayerAcousticSpectrum};

/// Hard ceiling on cached coherent plane-wave components.
pub const MAX_ANGULAR_SPECTRUM_COMPONENTS: usize = 4096;

/// Caller-owned allowance for preparation storage and work per observation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AngularSpectrumBudget {
    /// Nonempty component count, bounded by [`MAX_ANGULAR_SPECTRUM_COMPONENTS`].
    /// Zero-amplitude components also count against this allowance.
    pub max_components: usize,
}

/// One coherent incident plane wave at the original stack's frequency.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IncidentPlaneWave {
    /// Unit propagation vector in the interface frame; positive z enters the
    /// layers. Its norm must equal one to floating-point rounding accuracy.
    pub direction: [f64; 3],
    /// Complex peak pressure [Pa] at the entrance origin. Include any angular
    /// quadrature weight here; no area, solid-angle or square-root weight is
    /// inferred. Components with the same direction interfere coherently.
    pub pressure_pa: C64,
}

/// Coherent exterior fields and their time-averaged energy flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticSpectrumSample {
    pub pressure_pa: C64,
    pub particle_velocity_m_s: [C64; 3],
    /// `Re(p conj(v))/2` [W/m²], formed after coherent superposition.
    pub active_intensity_w_m2: [f64; 3],
    /// `Im(p conj(v))/2` [W/m²]; reactive flow is not absorbed heat.
    pub reactive_intensity_w_m2: [f64; 3],
}

/// Original material failure, explicit input/arithmetic refusal, or cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum AcousticSpectrumError {
    Material(MaterialError),
    Invalid(&'static str),
    Cancelled,
}
impl From<MaterialError> for AcousticSpectrumError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl core::fmt::Display for AcousticSpectrumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coherent acoustic spectrum refused: {self:?}")
    }
}
impl std::error::Error for AcousticSpectrumError {}

#[derive(Clone, Copy)]
struct PlaneField {
    pressure: C64,
    velocity: [C64; 3],
}
impl PlaneField {
    fn new(pressure: C64, velocity_per_pa: [C64; 3]) -> Result<Self, AcousticSpectrumError> {
        let velocity = velocity_per_pa.map(|value| value * pressure);
        if !finite(pressure) || velocity.iter().any(|&v| !finite(v)) {
            return Err(AcousticSpectrumError::Invalid(
                "unrepresentable acoustic component pressure or velocity",
            ));
        }
        Ok(Self { pressure, velocity })
    }
}

pub(in crate::porous::layers) struct PreparedPlaneWave {
    pub(in crate::porous::layers) cosine: f64,
    tangential_wavenumber: [f64; 2],
    incident_normal_wavenumber: f64,
    receiving_normal_wavenumber: C64,
    incident: PlaneField,
    reflected: PlaneField,
    transmitted: PlaneField,
}
impl PreparedPlaneWave {
    pub(in crate::porous::layers) fn incident_pressure(&self) -> C64 {
        self.incident.pressure
    }
}

/// Reusable coherent field on both exterior sides of ONE original layer stack.
///
/// Preparation caches each component's original transmission, reflection and
/// pressure/velocity coefficients. Observations cost O(components), use fixed
/// local accumulation storage, allocate nothing and query no material source.
/// Cancellation is checked per component and before publication; preparation
/// additionally checks at each original layer transfer. The borrowed stack
/// retains the source identities and all constitutive/no-claim boundaries.
pub struct PreparedAcousticSpectrum<'p, 'a> {
    stack: &'p TransmittingLayers<'a>,
    components: Vec<PreparedPlaneWave>,
}

pub(in crate::porous::layers) fn poll(
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AcousticSpectrumError> {
    if cancelled() {
        Err(AcousticSpectrumError::Cancelled)
    } else {
        Ok(())
    }
}

impl<'a> TransmittingLayers<'a> {
    /// Prepare coherent exterior waves using this stack's original R/T solver.
    /// No original incident, layer or receiving material is sampled again.
    /// The supplied component order is retained for deterministic accumulation.
    pub fn prepare_angular_spectrum<'p>(
        &'p self,
        incident: &[IncidentPlaneWave],
        budget: AngularSpectrumBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedAcousticSpectrum<'p, 'a>, AcousticSpectrumError> {
        let components = prepare_components(
            self.prepared_layers(),
            self.receiving_state().density_kg_m3(),
            incident,
            budget,
            cancelled,
            |cosine, cancelled| {
                self.scattering_with_layers::<AcousticSpectrumError>(cosine, |_, _| poll(cancelled))
            },
        )?;
        Ok(PreparedAcousticSpectrum {
            stack: self,
            components,
        })
    }
}

// Fluid and elastic stacks retain their own boundary solves while sharing the
// exact exterior momentum, phase convention, admission and coherent summation.
pub(in crate::porous::layers) fn prepare_components<C: FnMut() -> bool>(
    prepared: &RigidBackedLayers<'_>,
    receiving_density: f64,
    incident: &[IncidentPlaneWave],
    budget: AngularSpectrumBudget,
    cancelled: &mut C,
    mut scattering: impl FnMut(f64, &mut C) -> Result<LayeredTransmission, AcousticSpectrumError>,
) -> Result<Vec<PreparedPlaneWave>, AcousticSpectrumError> {
    poll(cancelled)?;
    if budget.max_components == 0
        || budget.max_components > MAX_ANGULAR_SPECTRUM_COMPONENTS
        || incident.is_empty()
        || incident.len() > budget.max_components
    {
        return Err(AcousticSpectrumError::Invalid(
            "acoustic spectrum needs nonempty components within explicit and hard budgets",
        ));
    }
    let k0 = prepared.exterior_wavenumber_rad_m();
    let omega = prepared.angular_frequency_rad_s();
    let impedance = prepared
        .fluid()
        .acoustic()
        .characteristic_impedance_pa_s_m();
    let mut components = Vec::with_capacity(incident.len());
    for wave in incident {
        poll(cancelled)?;
        let direction = wave.direction;
        let norm_squared = direction.iter().map(|d| d * d).sum::<f64>();
        if direction.iter().any(|d| !d.is_finite())
            || direction[2] <= 0.0
            || direction[2] > 1.0
            || (norm_squared - 1.0).abs() > 64.0 * f64::EPSILON
            || !finite(wave.pressure_pa)
        {
            return Err(AcousticSpectrumError::Invalid(
                "incident wave needs a unit direction with positive normal component and finite pressure",
            ));
        }
        let response = scattering(direction[2], cancelled)?;
        let tangential_wavenumber = [k0 * direction[0], k0 * direction[1]];
        let incident_normal_wavenumber = k0 * direction[2];
        let receiving_normal_wavenumber = response.receiving_normal_wavenumber_rad_m;
        let incident_velocity = direction.map(|d| C64::from_re(d / impedance));
        let reflected_velocity = [
            incident_velocity[0],
            incident_velocity[1],
            -incident_velocity[2],
        ];
        let receiving_velocity = [
            C64::from_re((tangential_wavenumber[0] / omega) / receiving_density),
            C64::from_re((tangential_wavenumber[1] / omega) / receiving_density),
            C64::new(
                (receiving_normal_wavenumber.re / omega) / receiving_density,
                (receiving_normal_wavenumber.im / omega) / receiving_density,
            ),
        ];
        components.push(PreparedPlaneWave {
            cosine: direction[2],
            tangential_wavenumber,
            incident_normal_wavenumber,
            receiving_normal_wavenumber,
            incident: PlaneField::new(wave.pressure_pa, incident_velocity)?,
            reflected: PlaneField::new(
                wave.pressure_pa * response.pressure_reflection,
                reflected_velocity,
            )?,
            transmitted: PlaneField::new(
                wave.pressure_pa * response.pressure_transmission,
                receiving_velocity,
            )?,
        });
    }
    poll(cancelled)?;
    Ok(components)
}

impl<'p, 'a> PreparedAcousticSpectrum<'p, 'a> {
    /// Original prepared layer stack, including its original source owners.
    pub const fn stack(&self) -> &'p TransmittingLayers<'a> {
        self.stack
    }

    pub fn component_count(&self) -> usize {
        self.components.len()
    }

    /// TOTAL incident plus reflected field at `[x,y,z]` metres, with `z<=0`
    /// measured from the entrance face. The reflected normal velocity has
    /// the outgoing sign; interference is retained in the resulting intensity.
    pub fn incident_field(
        &self,
        position_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcousticSpectrumSample, AcousticSpectrumError> {
        sample_components(&self.components, position_m, false, cancelled)
    }

    /// Transmitted field at `[x,y,z]` metres, with `z>=0` measured from the
    /// EXIT face. Its transverse origin is aligned with the entrance origin.
    /// Critical waves have zero normal wavenumber; evanescent waves decay
    /// into this half-space without being relabeled as absorbed power.
    pub fn transmitted_field(
        &self,
        position_m: [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AcousticSpectrumSample, AcousticSpectrumError> {
        sample_components(&self.components, position_m, true, cancelled)
    }
}

pub(in crate::porous::layers) fn sample_components(
    components: &[PreparedPlaneWave],
    position: [f64; 3],
    transmitted: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AcousticSpectrumSample, AcousticSpectrumError> {
    poll(cancelled)?;
    if position.iter().any(|x| !x.is_finite())
        || (transmitted && position[2] < 0.0)
        || (!transmitted && position[2] > 0.0)
    {
        return Err(AcousticSpectrumError::Invalid(
            "observation needs finite coordinates on the selected exterior side",
        ));
    }
    let mut sum = FieldSum([Dd::ZERO; 8]);
    for wave in components {
        poll(cancelled)?;
        if transmitted {
            sum.add(
                wave.transmitted,
                propagation(wave, position, wave.receiving_normal_wavenumber)?,
            )?;
        } else {
            sum.add(
                wave.incident,
                propagation(
                    wave,
                    position,
                    C64::from_re(wave.incident_normal_wavenumber),
                )?,
            )?;
            sum.add(
                wave.reflected,
                propagation(
                    wave,
                    position,
                    C64::from_re(-wave.incident_normal_wavenumber),
                )?,
            )?;
        }
    }
    let sample = sum.finish()?;
    poll(cancelled)?;
    Ok(sample)
}

fn propagation(
    wave: &PreparedPlaneWave,
    position: [f64; 3],
    normal: C64,
) -> Result<C64, AcousticSpectrumError> {
    let decay = det::exp(normal.im * position[2]);
    if decay == 0.0 {
        return Ok(C64::ZERO);
    }
    let phase = -(wave.tangential_wavenumber[0] * position[0]
        + wave.tangential_wavenumber[1] * position[1]
        + normal.re * position[2]);
    if !phase.is_finite() || !decay.is_finite() {
        return Err(AcousticSpectrumError::Invalid(
            "unrepresentable acoustic propagation phase or decay",
        ));
    }
    let (sin, cos) = det::sin_cos(phase);
    Ok(C64::new(cos, sin).scale(decay))
}

// Fixed storage and compensated coherent summation retain interference nodes
// without accumulating powers separately or allocating a sample-sized buffer.
struct FieldSum([Dd; 8]);
impl FieldSum {
    fn add(&mut self, field: PlaneField, propagation: C64) -> Result<(), AcousticSpectrumError> {
        for (index, value) in core::iter::once(field.pressure)
            .chain(field.velocity)
            .enumerate()
        {
            let value = value * propagation;
            if !finite(value) {
                return Err(AcousticSpectrumError::Invalid(
                    "nonfinite acoustic component field",
                ));
            }
            self.0[2 * index] = self.0[2 * index] + Dd::from_f64(value.re);
            self.0[2 * index + 1] = self.0[2 * index + 1] + Dd::from_f64(value.im);
        }
        Ok(())
    }

    fn finish(self) -> Result<AcousticSpectrumSample, AcousticSpectrumError> {
        let fields: [C64; 4] =
            core::array::from_fn(|i| C64::new(self.0[2 * i].to_f64(), self.0[2 * i + 1].to_f64()));
        let pressure = fields[0];
        let velocity = [fields[1], fields[2], fields[3]];
        let intensity = velocity.map(|v| pressure.scale(0.5) * v.conj());
        if fields.iter().chain(&intensity).any(|&value| !finite(value)) {
            return Err(AcousticSpectrumError::Invalid(
                "nonfinite coherent acoustic field or intensity",
            ));
        }
        Ok(AcousticSpectrumSample {
            pressure_pa: pressure,
            particle_velocity_m_s: velocity,
            active_intensity_w_m2: intensity.map(|i| i.re),
            reactive_intensity_w_m2: intensity.map(|i| i.im),
        })
    }
}
