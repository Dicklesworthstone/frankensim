//! Coherent refraction through fluid-coupled elastic plates and original gaps.
//!
//! Every nonzero plate is homogeneous, explicitly isotropic and lossless. Its
//! two faces retain the full P/S potential solution of the single-plate owner:
//! normal traction/velocity are continuous and shear traction is zero. Positive
//! fluid or equivalent-fluid separation is required between nonzero plates.
//! Adjacent solids are NOT treated as a bonded laminate or sliding contact.
//!
//! The same tangential wavenumber crosses every section, including evanescent
//! gaps and plates. A homogeneous pressure/velocity pair and compensated log
//! amplitude propagate backwards through the original bounded kernels; no
//! admittance division, growing transfer matrix, fitted plate theory or
//! multiplication of incoherent power transmissions is introduced.
//!
//! The plate boundary equations are those of Almeida et al.,
//! <https://arxiv.org/abs/2302.08826>. Each separating fluid uses the original
//! pressure/normal-velocity transfer described in the parent layers module.

use super::super::absorption::{
    LayerAbsorbedPower, LayerAbsorption, LayerAbsorptionBudget, LayerAbsorptionError, admit,
    finish_partition, layer_loss, poll,
};
use super::super::transmission::{
    LayeredTransmission, LogAmplitude, TransmittedState, log_scaled, positive_product_ratio,
};
use super::super::{
    Layer, LayerTransfer, MAX_LAYERS, PreparedLayer, RigidBackedLayers, bad, divided, finite,
    magnitude_scale, normal_wavenumber_ratio, retained_fluid_snapshot,
};
use super::plate::{PlateBoundary, bounded_wave, parity_boundary_pair};
use crate::MaterialError;
use crate::elastic::acoustic::ElasticAcousticState;
use crate::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::{c64::C64, dd::Dd, det};

/// Ordered physical sections from the incident fluid to the receiving fluid.
/// Zero thickness is an identity but does not waive original source admission.
#[derive(Clone, Copy, Debug)]
pub enum ElasticPlateStackLayer<'a> {
    /// Original isotropic solid with independent shear-free faces.
    Plate {
        state: &'a ElasticAcousticState,
        thickness_m: f64,
    },
    /// Existing lossless fluid or rigid-frame porous equivalent-fluid layer.
    Fluid(Layer<'a>),
}

impl ElasticPlateStackLayer<'_> {
    fn thickness_m(self) -> f64 {
        match self {
            Self::Plate { thickness_m, .. } => thickness_m,
            Self::Fluid(layer) => layer.thickness_m(),
        }
    }
}

struct PreparedPlate<'a> {
    source: &'a ElasticAcousticState,
    thickness: f64,
    longitudinal: FluidAcousticState,
    shear_speed: f64,
    longitudinal_ratio: f64,
    shear_ratio: f64,
    impedance_ratio: f64,
}

enum PreparedSection<'a> {
    Fluid { source: Layer<'a>, index: usize },
    Plate(PreparedPlate<'a>),
}

// A section's original bounded pressure/velocity pair and the two separately
// accumulated log factors relating its faces. A plate's signed transmission
// only changes their common phase, which cancels in every local work integral.
#[derive(Clone, Copy, Default)]
struct SectionTransfer {
    fields: LayerTransfer,
    log_factors: [f64; 2],
}

/// At most [`MAX_LAYERS`] fluid-coupled sections at one original frequency.
/// Preparation allocates bounded storage and resolves each distinct fluid once.
/// Incidence queries allocate nothing and re-query no material source.
///
/// The original static elastic response regime and strain ceiling are retained;
/// these infinitesimal coefficients do not admit a finite driving amplitude.
/// Absorption is solely in the original porous layers. Solid damping, bonded
/// laminates, lateral boundaries and reconstructed interior fields are absent.
pub struct ElasticPlateStack<'a> {
    prepared: RigidBackedLayers<'a>,
    sections: Vec<PreparedSection<'a>>,
    receiver_source: &'a dyn FluidAcoustics,
    receiver: FluidAcousticState,
    speed_ratio: f64,
    density_ratio: f64,
    lossless: bool,
}

impl<'a> ElasticPlateStack<'a> {
    pub fn new(
        incident: &'a dyn ThermoviscousFluid,
        receiver_source: &'a dyn FluidAcoustics,
        omega: f64,
        layers: &[ElasticPlateStackLayer<'a>],
    ) -> Result<Self, MaterialError> {
        if layers.len() > MAX_LAYERS
            || layers
                .iter()
                .any(|layer| !layer.thickness_m().is_finite() || layer.thickness_m() < 0.0)
        {
            return Err(bad("elastic stack needs <=64 finite nonnegative sections"));
        }
        let mut fluids = Vec::with_capacity(layers.len());
        let mut separated = true;
        for &layer in layers {
            match layer {
                ElasticPlateStackLayer::Plate { state, thickness_m } => {
                    if !state.is_explicitly_isotropic() {
                        return Err(bad(
                            "elastic stack requires explicitly isotropic plate sources",
                        ));
                    }
                    if thickness_m > 0.0 {
                        if !separated {
                            return Err(bad(
                                "nonzero plates require positive fluid separation; bonded solids are unsupported",
                            ));
                        }
                        separated = false;
                    }
                }
                ElasticPlateStackLayer::Fluid(fluid) => {
                    separated |= fluid.thickness_m() > 0.0;
                    fluids.push(fluid);
                }
            }
        }
        // Compress only the preparation list. The complete physical order is
        // retained below, and all pores share this ONE incident fluid snapshot.
        let prepared = RigidBackedLayers::new(incident, omega, &fluids)?;
        let exterior = prepared.fluid().acoustic();
        let mut sections = Vec::with_capacity(layers.len());
        let mut fluid_index = 0;
        for &layer in layers {
            match layer {
                ElasticPlateStackLayer::Fluid(source) => {
                    sections.push(PreparedSection::Fluid {
                        source,
                        index: fluid_index,
                    });
                    fluid_index += 1;
                }
                ElasticPlateStackLayer::Plate { state, thickness_m } => {
                    let modes = state.directional_waves([0.0, 0.0, 1.0])?.modes;
                    let longitudinal = FluidAcousticState::try_new(
                        state.density_kg_m3(),
                        modes[2].phase_speed_m_s,
                    )?;
                    let shear_speed = modes[0].phase_speed_m_s;
                    let longitudinal_ratio =
                        exterior.sound_speed_m_s() / longitudinal.sound_speed_m_s();
                    let shear_ratio = exterior.sound_speed_m_s() / shear_speed;
                    let impedance_ratio =
                        exterior.characteristic_impedance_pa_s_m() / modes[2].impedance_pa_s_m;
                    if [longitudinal_ratio, shear_ratio, impedance_ratio]
                        .iter()
                        .any(|value| !value.is_finite() || *value <= 0.0)
                    {
                        return Err(bad("elastic stack contrast is unrepresentable"));
                    }
                    sections.push(PreparedSection::Plate(PreparedPlate {
                        source: state,
                        thickness: thickness_m,
                        longitudinal,
                        shear_speed,
                        longitudinal_ratio,
                        shear_ratio,
                        impedance_ratio,
                    }));
                }
            }
        }
        let original_incident: &dyn FluidAcoustics = incident;
        let receiver = match retained_fluid_snapshot(
            receiver_source,
            original_incident,
            exterior,
            &prepared.layers,
        ) {
            Some(state) => state,
            None => receiver_source.acoustic_state()?,
        };
        let speed_ratio = exterior.sound_speed_m_s() / receiver.sound_speed_m_s();
        let density_ratio = exterior.density_kg_m3() / receiver.density_kg_m3();
        if [speed_ratio, density_ratio]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(bad("elastic stack receiving contrast is unrepresentable"));
        }
        let lossless = prepared
            .layers
            .iter()
            .all(|layer| layer.porous.is_none() || layer.source.thickness_m() == 0.0);
        Ok(Self {
            prepared,
            sections,
            receiver_source,
            receiver,
            speed_ratio,
            density_ratio,
            lossless,
        })
    }

    /// Original ordered section sources, including zero-thickness entries.
    pub fn layers(&self) -> impl ExactSizeIterator<Item = ElasticPlateStackLayer<'a>> + '_ {
        self.sections.iter().map(|section| match section {
            PreparedSection::Fluid { source, .. } => ElasticPlateStackLayer::Fluid(*source),
            PreparedSection::Plate(plate) => ElasticPlateStackLayer::Plate {
                state: plate.source,
                thickness_m: plate.thickness,
            },
        })
    }

    /// Original physical sections, including lossless plates and identities.
    pub fn layer_count(&self) -> usize {
        self.sections.len()
    }

    pub const fn incident_source(&self) -> &'a dyn ThermoviscousFluid {
        self.prepared.source
    }

    pub const fn incident_state(&self) -> &ThermoviscousFluidState {
        &self.prepared.fluid
    }

    pub const fn receiving_source(&self) -> &'a dyn FluidAcoustics {
        self.receiver_source
    }

    pub const fn receiving_state(&self) -> &FluidAcousticState {
        &self.receiver
    }

    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.prepared.omega
    }

    /// Coherent pressure R at the entrance and T at the final receiving face.
    /// Original fluid/solid critical angles and evanescent gaps are retained.
    pub fn scattering(&self, cosine: f64) -> Result<LayeredTransmission, MaterialError> {
        self.scattering_with_sections(cosine, |_, _| Ok(()))
    }

    /// Resolve acoustic heat in each original porous section, including layers
    /// before, between and after elastic plates. Rows and the work budget count
    /// ALL declared sections in their physical order. Lossless plates, fluid
    /// gaps and zero-thickness identities retain exact-zero heat rows.
    ///
    /// The existing coherent P/S plate and equivalent-fluid propagation supplies
    /// both faces of every layer. Positive analytic viscous/thermal work retains
    /// losses below the resolution of `1-R-T`, including oblique tangential
    /// pore velocity. Transmitted sound and reactive evanescence are not heat.
    /// No material is re-queried, no additional field solve or spatial quadrature
    /// occurs, and no temperature or fluid/frame heat partition is inferred.
    /// Fixed storage and two linear passes are bounded by `max_layers`, with
    /// cancellation per section in both passes and before publication.
    pub fn layer_absorption(
        &self,
        cosine: f64,
        incident_power_w: f64,
        budget: LayerAbsorptionBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LayerAbsorption, LayerAbsorptionError> {
        admit(
            self.layer_count(),
            cosine,
            incident_power_w,
            budget,
            &mut cancelled,
        )?;
        let mut transfers = [SectionTransfer::default(); MAX_LAYERS];
        let response =
            self.scattering_with_sections::<LayerAbsorptionError>(cosine, |i, step| {
                poll(&mut cancelled)?;
                transfers[i] = step;
                Ok(())
            })?;
        // Unit incident pressure means p_front + Z_ref*v_front = 2. Absolute
        // phase is unnecessary: the existing integral needs only the common
        // face magnitudes and each layer's retained internal phase relation.
        let front = transfers[0].fields;
        let mut amplitude = LogAmplitude::new(if self.sections.is_empty() {
            0.0
        } else {
            core::f64::consts::LN_2 - det::ln((front.front_pressure + front.front_velocity).abs())
        });
        let mut layers = [LayerAbsorbedPower::default(); MAX_LAYERS];
        let mut total = Dd::ZERO;
        for (i, section) in self.sections.iter().enumerate() {
            poll(&mut cancelled)?;
            let step = transfers[i];
            let front_log = amplitude.log;
            for factor in step.log_factors {
                amplitude.add(factor);
            }
            if let PreparedSection::Fluid { source, index } = section
                && source.thickness_m() > 0.0
                && let Some(porous) = self.prepared.layers[*index].porous
            {
                let [viscous, thermal] = layer_loss(
                    &self.prepared,
                    cosine,
                    source.thickness_m(),
                    &porous,
                    step.fields,
                    [front_log, amplitude.log],
                )?;
                let power = incident_power_w * (viscous + thermal);
                if !power.is_finite() {
                    return Err(LayerAbsorptionError::Invalid(
                        "elastic stack layer absorbed power is unrepresentable",
                    ));
                }
                layers[i] = LayerAbsorbedPower {
                    viscous_absorbed_fraction: viscous,
                    thermal_absorbed_fraction: thermal,
                    absorbed_power_w: power,
                };
                total = total + Dd::from_f64(viscous) + Dd::from_f64(thermal);
            }
        }
        finish_partition(
            incident_power_w,
            [response.reflected_fraction, response.transmitted_fraction],
            layers,
            self.layer_count(),
            total,
            &mut cancelled,
        )
    }

    fn scattering_with_sections<E: From<MaterialError>>(
        &self,
        cosine: f64,
        mut observe: impl FnMut(usize, SectionTransfer) -> Result<(), E>,
    ) -> Result<LayeredTransmission, E> {
        if !cosine.is_finite() || cosine <= 0.0 || cosine > 1.0 {
            return Err(bad("elastic stack needs incidence cosine in (0,1]").into());
        }
        let q = normal_wavenumber_ratio(
            self.speed_ratio,
            cosine,
            self.prepared.fluid().acoustic().sound_speed_m_s(),
            self.receiver.sound_speed_m_s(),
        );
        let receiving_kn = q.scale(self.prepared.k0);
        if !finite(receiving_kn) || (q != C64::ZERO && receiving_kn == C64::ZERO) {
            return Err(bad("elastic stack receiving wavenumber is unrepresentable").into());
        }
        let load = if q == C64::ZERO {
            C64::ZERO
        } else {
            let value = positive_product_ratio(self.density_ratio, magnitude_scale(q), cosine)?;
            if q.re > 0.0 {
                C64::from_re(value)
            } else {
                C64::new(0.0, -value)
            }
        };
        let initial_scale = magnitude_scale(load).max(1.0);
        let mut pressure = C64::from_re(1.0 / initial_scale);
        let mut velocity = divided(load, initial_scale);
        let mut amplitude = LogAmplitude::new(-det::ln(initial_scale));
        let mut sign = 1.0;
        for (section_index, section) in self.sections.iter().enumerate().rev() {
            let mut step = SectionTransfer {
                fields: LayerTransfer::identity(pressure, velocity),
                log_factors: [0.0; 2],
            };
            match section {
                PreparedSection::Fluid { source, index } => {
                    let thickness = source.thickness_m();
                    if thickness == 0.0 {
                        observe(section_index, step)?;
                        continue;
                    }
                    let normalization;
                    let kn;
                    (pressure, velocity, normalization, kn) = self.prepared.propagate_layer(
                        &self.prepared.layers[*index],
                        cosine,
                        pressure,
                        velocity,
                    )?;
                    step.fields.front_pressure = pressure;
                    step.fields.front_velocity = velocity;
                    step.fields.normalization = normalization;
                    step.fields.normal_wavenumber = kn;
                    if kn.im == 0.0 {
                        step.log_factors[0] = -log_scaled(16.0, normalization);
                        amplitude.add(step.log_factors[0]);
                    } else {
                        step.log_factors[0] = kn.im * thickness - log_scaled(8.0, normalization);
                        amplitude.add(step.log_factors[0]);
                        amplitude.rotate(-kn.re * thickness);
                    }
                }
                PreparedSection::Plate(plate) => {
                    if plate.thickness == 0.0 {
                        observe(section_index, step)?;
                        continue;
                    }
                    let next = plate.propagate(&self.prepared, cosine, pressure, velocity)?;
                    let scale = magnitude_scale(next.pressure).max(magnitude_scale(next.velocity));
                    if !finite(next.pressure) || !finite(next.velocity) || scale == 0.0 {
                        return Err(bad("elastic stack boundary pair is unrepresentable").into());
                    }
                    pressure = divided(next.pressure, scale);
                    velocity = divided(next.velocity, scale);
                    step.fields.front_pressure = pressure;
                    step.fields.front_velocity = velocity;
                    // Keep every signed transmission factor in log magnitude;
                    // an exact zero remains zero without an admittance pole.
                    if next.transmission < 0.0 {
                        sign = -sign;
                    }
                    step.log_factors = [
                        det::ln(next.transmission.abs()),
                        next.transmission_log_scale - det::ln(scale),
                    ];
                    amplitude.add(step.log_factors[0]);
                    amplitude.add(step.log_factors[1]);
                }
            }
            observe(section_index, step)?;
        }
        let state = TransmittedState::from_pair(pressure, velocity, amplitude)?;
        let transmission = state.coefficient(C64::from_re(sign))?;
        let transmitted = if load.re == 0.0 {
            0.0
        } else {
            state.power(C64::ONE, det::ln(load.re))
        };
        let reflected = state.reflection.norm_sq();
        let absorbed = if self.lossless {
            0.0
        } else {
            state.entering - transmitted
        };
        let tolerance = 4096.0 * f64::EPSILON;
        if !transmitted.is_finite()
            || transmitted < 0.0
            || !absorbed.is_finite()
            || absorbed < -tolerance
            || (reflected + transmitted + absorbed - 1.0).abs() > tolerance
        {
            return Err(bad("elastic stack failed passive normal-power balance").into());
        }
        Ok(LayeredTransmission {
            pressure_reflection: state.reflection,
            pressure_transmission: transmission,
            receiving_normal_wavenumber_rad_m: receiving_kn,
            reflected_fraction: reflected,
            transmitted_fraction: transmitted,
            absorbed_fraction: absorbed,
        })
    }
}

impl PreparedPlate<'_> {
    fn propagate(
        &self,
        prepared: &RigidBackedLayers<'_>,
        cosine: f64,
        pressure: C64,
        velocity: C64,
    ) -> Result<PlateBoundary, MaterialError> {
        let cp = self.longitudinal.sound_speed_m_s();
        if cosine == 1.0 {
            // The exact longitudinal limit uses the original fluid kernel;
            // it does not approximate an oblique elastic plate by a fluid.
            let initial = magnitude_scale(pressure).max(magnitude_scale(velocity));
            let layer = PreparedLayer {
                source: Layer::FluidLayer {
                    fluid: &self.longitudinal,
                    thickness_m: self.thickness,
                },
                porous: None,
                fluid: Some(self.longitudinal),
            };
            let (pressure, velocity, scale, _) = prepared.propagate_layer(
                &layer,
                1.0,
                divided(pressure, initial),
                divided(velocity, initial),
            )?;
            return Ok(PlateBoundary {
                pressure,
                velocity,
                transmission: 1.0,
                transmission_log_scale: -det::ln(initial) - log_scaled(16.0, scale),
            });
        }
        let ci = prepared.fluid().acoustic().sound_speed_m_s();
        let qp = normal_wavenumber_ratio(self.longitudinal_ratio, cosine, ci, cp);
        let qs = normal_wavenumber_ratio(self.shear_ratio, cosine, ci, self.shear_speed);
        let p_kn = qp.scale(prepared.k0);
        let s_kn = qs.scale(prepared.k0);
        for (q, kn) in [(qp, p_kn), (qs, s_kn)] {
            if !finite(kn) || (q != C64::ZERO && kn == C64::ZERO) {
                return Err(bad("elastic stack plate wavenumber is unrepresentable"));
            }
        }
        let p_cosine = divided(qp, self.longitudinal_ratio);
        let s_cosine = divided(qs, self.shear_ratio);
        let sine = det::sqrt((1.0 - cosine) * (1.0 + cosine)) / self.shear_ratio;
        let zeta = self.impedance_ratio / cosine;
        if !finite(p_cosine)
            || !finite(s_cosine)
            || !sine.is_finite()
            || !zeta.is_finite()
            || zeta <= 0.0
        {
            return Err(bad("elastic stack plate refraction is unrepresentable"));
        }
        parity_boundary_pair(
            bounded_wave(p_cosine, p_kn, prepared.omega / cp, self.thickness)?,
            bounded_wave(
                s_cosine,
                s_kn,
                prepared.omega / self.shear_speed,
                self.thickness,
            )?,
            sine,
            self.shear_speed / cp,
            zeta,
            pressure,
            velocity,
        )
    }
}
