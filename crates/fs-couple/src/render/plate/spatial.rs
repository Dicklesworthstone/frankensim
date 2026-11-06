//! Stationary finite-wavelength pressure from a plate's original signed modes.
//!
//! This adapter retains the existing exact-ZOH body mechanics and the shared
//! differentiating Rayleigh receiver. It supplies neither a second radiation
//! kernel nor a reaction load. The model is a flat panel in an infinite rigid
//! baffle, observed in its +z half-space. A zero net monopole area does not erase
//! its moving surface or its off-axis radiation.
use super::{CompactPlateVoice, PlateVoiceConfig, RenderError, admit_radiators, invalid};
use crate::acoustic_realize::AcousticRealizeError;
use crate::pcm_wav::baffled::{BaffledPressure, RayleighMedium};
use crate::pcm_wav::baffled::gas::{GasPropagation, GasPropagationReport};
use crate::thin_plate::PlateSurfaceRadiators;

/// Actual path geometry, admitted band and numerical work of a surface observer.
/// Spatial refinement controls phase variation and near-field quadrature; these
/// controls are not a rigorous integration-error or modal-completeness bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlateSpatialObserverReceipt {
    /// Receiver position in the original chart's flat baffle frame [m].
    pub position_m: [f64; 3],
    /// Number of signed surface quadrature points retained from the reduction.
    pub surface_samples: usize,
    /// Diameter of the largest retained quadrature cell [m].
    pub maximum_cell_diameter_m: f64,
    /// Declared spatially resolved upper frequency [Hz], including every mode.
    /// This is an admission band, not a filter on a subsequently applied force.
    pub maximum_frequency_hz: f64,
    /// Minimum and maximum physical travel times from the actual samples [s].
    /// Reconstruction latency is separate and is never counted as gas travel.
    pub propagation_delay_s: [f64; 2],
    /// Additional bandlimited derivative/interpolation delay, in output samples.
    pub filter_delay_samples: f64,
    /// Receiver multiply-add estimate per sample, excluding body mechanics.
    pub multiply_adds_per_sample: usize,
}

pub(super) struct SpatialPlateObserver {
    receiver: BaffledPressure,
    velocities: Vec<f64>,
    receipt: PlateSpatialObserverReceipt,
    gas_propagation: Option<GasPropagationReport>,
}

impl CompactPlateVoice {
    /// Host one geometry-derived modal reduction with its signed physical
    /// surface and observe it through the original retarded Rayleigh owner.
    /// All modal coordinates, masses, forces and optional reaction loads remain
    /// unchanged. `config.listener_m` is only a reference range in this image;
    /// each actual surface-to-receiver path supplies its own attenuation and
    /// time of flight through the original gas.
    ///
    /// Before constructing the receiver, the largest cell must satisfy
    /// `k_max * diameter <= pi/4` and `diameter <= receiver_z/4`. The declared
    /// band must contain every retained natural frequency and be at most 0.45
    /// of the sample rate. These are explicit resolution controls; convergence
    /// still requires spatial/temporal refinement and enough structural modes.
    /// The force is not silently bandlimited. Discontinuous controls can excite
    /// content outside the declared band, where no accuracy claim is made.
    ///
    /// The differentiating/interpolating receiver adds its reported numerical
    /// latency to the physical path delays. No compact delay or postmixed
    /// single-range atmosphere filter is applied a second time.
    ///
    /// # Errors
    /// Invalid body/configuration, underresolved geometry, insufficient band,
    /// invalid half-space receiver, or the shared receiver's bounded admission.
    pub fn from_surface_radiators(
        reduction: PlateSurfaceRadiators,
        initial_force_n: f64,
        config: PlateVoiceConfig,
        receiver_position_m: [f64; 3],
        maximum_frequency_hz: f64,
    ) -> Result<Self, RenderError> {
        Self::from_surface_radiators_with_propagation(reduction, initial_force_n, config,
            receiver_position_m, maximum_frequency_hz, None, &mut || false)
    }

    /// Observe the same signed plate surface through the shared material-gas
    /// Rayleigh owner. Each quadrature path retains its own absorption,
    /// dispersion and travel time; the declared gas is also the radiation
    /// medium. No mean-range filter or additional compact delay is applied.
    ///
    /// Original modal motion, force participation and reaction loads remain
    /// unchanged. Preparation polls cancellation and honors the supplied work
    /// limits. The returned report separates common numerical padding from
    /// physical propagation. This observation adds no mechanical loss or heat.
    #[allow(clippy::too_many_arguments)]
    pub fn from_surface_radiators_with_gas_propagation(
        reduction: PlateSurfaceRadiators,
        initial_force_n: f64,
        config: PlateVoiceConfig,
        receiver_position_m: [f64; 3],
        maximum_frequency_hz: f64,
        propagation: &GasPropagation,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RenderError> {
        Self::from_surface_radiators_with_propagation(reduction, initial_force_n, config,
            receiver_position_m, maximum_frequency_hz, Some(propagation), cancelled)
    }

    #[allow(clippy::too_many_arguments)]
    fn from_surface_radiators_with_propagation(
        reduction: PlateSurfaceRadiators,
        initial_force_n: f64,
        config: PlateVoiceConfig,
        receiver_position_m: [f64; 3],
        maximum_frequency_hz: f64,
        propagation: Option<&GasPropagation>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RenderError> {
        if cancelled() { return Err(RenderError::Voice(AcousticRealizeError::Cancelled)); }
        admit_radiators(reduction.bodies(), config)?;
        if !initial_force_n.is_finite() || initial_force_n.abs() > config.maximum_abs_force_n
            || reduction.bodies().iter()
                .any(|b| !(initial_force_n * b.drive_participation).is_finite())
        {
            return Err(invalid("plate footprint force exceeds its finite newton budget"));
        }
        let diameter = reduction.maximum_cell_diameter_m();
        if !maximum_frequency_hz.is_finite() || maximum_frequency_hz <= 0.0
            || maximum_frequency_hz > 0.45 * f64::from(config.sample_rate_hz)
            || reduction.bodies().iter()
                .any(|b| b.omega / core::f64::consts::TAU > maximum_frequency_hz)
        {
            return Err(invalid("spatial plate band must contain every retained mode and remain at or below 0.45 of the sample rate"));
        }
        if receiver_position_m.iter().any(|v| !v.is_finite()) || receiver_position_m[2] < 0.05
            || !diameter.is_finite() || diameter <= 0.0
        {
            return Err(invalid("spatial plate needs a finite positive cell diameter and a receiver at least 5 cm above the baffle"));
        }
        let phase = core::f64::consts::TAU * maximum_frequency_hz * diameter
            / config.gas.sound_speed;
        if !phase.is_finite() || phase > core::f64::consts::FRAC_PI_4
            || diameter > 0.25 * receiver_position_m[2]
        {
            return Err(invalid("plate surface quadrature underresolves the declared wavelength or receiver clearance; refine its subdivisions"));
        }
        let mut propagation_delay_s = [f64::INFINITY, 0.0_f64];
        for sample in reduction.surface() {
            if cancelled() { return Err(RenderError::Voice(AcousticRealizeError::Cancelled)); }
            // Use the original receiver's metric and gas, with no reference
            // range substitution or modal-area averaging of path lengths.
            let radius = fs_math::det::sqrt((0..3)
                .map(|i| (receiver_position_m[i] - sample.position_m[i]).powi(2)).sum());
            let delay = radius / config.gas.sound_speed;
            if !delay.is_finite() || delay <= 0.0 {
                return Err(invalid("plate surface path must have a finite positive travel time"));
            }
            propagation_delay_s[0] = propagation_delay_s[0].min(delay);
            propagation_delay_s[1] = propagation_delay_s[1].max(delay);
        }
        let (mut receiver, gas_propagation) = match propagation {
            Some(propagation) => {
                let mut stopped = false;
                let result = BaffledPressure::new_multirate_gas_differentiating(
                    reduction.surface(), reduction.bodies().len(), config.sample_rate_hz,
                    config.sample_rate_hz, receiver_position_m, &config.gas, propagation,
                    &mut || { stopped = cancelled(); stopped },
                );
                let (receiver, report) = result.map_err(|error| if stopped {
                    RenderError::Voice(AcousticRealizeError::Cancelled)
                } else { receiver_error(error) })?;
                (receiver, Some(report))
            }
            None => (BaffledPressure::new_differentiating(
                reduction.surface(), reduction.bodies().len(), config.sample_rate_hz,
                receiver_position_m, RayleighMedium::from(&config.gas),
            ).map_err(receiver_error)?, None),
        };
        // At this one-to-one clock, differentiation commutes with the fixed
        // Rayleigh sum. Its shared cold optimization runs the same derivative
        // once per output rather than once per retained structural mode.
        receiver.differentiate_pressure().map_err(receiver_error)?;
        let receipt = PlateSpatialObserverReceipt {
            position_m: receiver_position_m,
            surface_samples: reduction.surface().len(),
            maximum_cell_diameter_m: diameter,
            maximum_frequency_hz,
            propagation_delay_s,
            filter_delay_samples: receiver.filter_delay_samples(),
            multiply_adds_per_sample: receiver.multiply_adds_per_sample(),
        };
        let velocities = vec![0.0; reduction.bodies().len()];
        let (bodies, _) = reduction.into_parts();
        if cancelled() { return Err(RenderError::Voice(AcousticRealizeError::Cancelled)); }
        Ok(Self {
            bodies, config, force_n: initial_force_n, sample_index: 0, poisoned: false,
            delay_whole: 0, retarded_elapsed_s: 0.0, history: Vec::new(), observer: None,
            spatial: Some(SpatialPlateObserver { receiver, velocities, receipt, gas_propagation }),
        })
    }

    /// Geometry, band, physical delay interval and numerical latency of the
    /// explicit surface observer. Compact observation returns `None`.
    #[must_use]
    pub fn spatial_observer_receipt(&self) -> Option<&PlateSpatialObserverReceipt> {
        self.spatial.as_ref().map(|spatial| &spatial.receipt)
    }

    /// Original per-path material propagation receipts, when explicitly selected.
    /// Absent for the original lossless spatial image and for compact observers.
    #[must_use]
    pub fn spatial_gas_propagation_report(&self) -> Option<&GasPropagationReport> {
        self.spatial.as_ref().and_then(|spatial| spatial.gas_propagation.as_ref())
    }

    // Public step_block has already admitted count, clock, poison and nonempty
    // output. Every buffer below belongs to the cold, fixed receiver image.
    pub(super) fn step_spatial_block(&mut self, out: &mut [f64]) -> Result<(), RenderError> {
        let dt = self.sample_period_s();
        let spatial = self.spatial.as_mut().expect("surface image selected by step_block");
        for slot in out {
            for (body, velocity) in self.bodies.iter_mut().zip(&mut spatial.velocities) {
                let force = self.force_n * body.drive_participation;
                if let Err(error) = body.drive(force, dt) {
                    self.poisoned = true;
                    return Err(RenderError::Voice(error));
                }
                *velocity = body.modal_velocity_m_s();
            }
            let pressure = match spatial.receiver.step_trace(&spatial.velocities) {
                Ok(pressure) => pressure,
                Err(error) => {
                    self.poisoned = true;
                    return Err(receiver_error(error));
                }
            };
            if !pressure.is_finite() || pressure.abs() > self.config.maximum_abs_pressure_pa {
                self.poisoned = true;
                return Err(invalid("spatial plate pressure exceeds its declared finite output budget"));
            }
            *slot = pressure;
            self.sample_index += 1;
        }
        Ok(())
    }
}

fn receiver_error(error: String) -> RenderError {
    RenderError::Voice(AcousticRealizeError::Nonlinear(error))
}

#[cfg(test)]
mod tests;
