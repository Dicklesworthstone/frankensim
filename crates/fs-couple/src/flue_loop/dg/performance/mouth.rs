//! Original rectangular Rayleigh observation and compact self reaction at the
//! actual DG flue mouth. The accepted inlet flow supplies its outward sign.

use super::{BaffledPressure, FlueDgPerformance, FlueError, FlueMouthRadiation};
use super::{GestureCompileError, RenderError, invalid, source_error};
use crate::acoustic_realize::AcousticRealizeError;
use crate::flue_loop::mouth_radiation::RectangularMouthLoad;
use crate::pcm_wav::baffled::openings::RectangularOutletReceiver;
use fs_exec::CancelGate;
use super::{GasPropagation, GasPropagationReport};

pub(super) struct MouthPressure {
    input: FlueMouthRadiation,
    receiver: BaffledPressure,
    observer_gas: Option<GasPropagationReport>,
}

impl MouthPressure {
    pub(super) fn push(&mut self, inward_flow_m3_s: f64) -> Result<f64, FlueError> {
        // DG inlet flow is positive into the bore. Both observed openings face
        // the same exterior, so mouth radiation takes its negative once.
        self.receiver.step_trace(&[-inward_flow_m3_s])
            .map_err(|error| FlueError::Physics(AcousticRealizeError::Nonlinear(error)))
    }
}

impl FlueDgPerformance {
    /// Add the solved mouth flow to this performance's exterior microphone.
    ///
    /// The acoustic rectangle is explicit, or defaults to the source width and
    /// cut-up. The complete exterior gas comes from the already computed terminal.
    /// The jet geometry remains independent. The original rectangular
    /// Rayleigh owner retains distance, phase, surface quadrature and numerical
    /// filter delay; no gain or duplicate pressure source is introduced.
    ///
    /// This explicitly describes disjoint openings in one infinite rigid baffle.
    /// Selecting observation alone retains the existing mouth reaction. The
    /// compact self load and complete mutual load are explicit separate choices
    /// through with_mouth_reaction and with_mutual_mouth_reaction. None resolves
    /// a jet or free-field scattering by an instrument body.
    pub fn with_mouth_radiation(mut self, input: FlueMouthRadiation)
        -> Result<Self, GestureCompileError>
    {
        if self.poisoned || self.voice.samples_rendered() != 0 || self.mouth.is_some()
            || self.observer_gas.is_some() {
            return Err(invalid("DG mouth observation requires one cold attachment"));
        }
        let outlet = self.outlet.as_ref()
            .ok_or_else(|| invalid("DG mouth observation requires the bound baffled outlet microphone"))?;
        let source = self.voice.dg().expect("admitted original DG voice").load()
            .radiation_record().expect("admitted computed DG microphone source");
        let size = input.size_m(&self.voice.geometry);
        if source.driven_mouth().is_some_and(|mouth|
            mouth.center_m != input.center_m || mouth.cells != input.cells || mouth.size_m != size)
        {
            return Err(invalid("joint DG mouth observation must retain its bound rectangular source"));
        }
        let exterior = source.gas_record().state();
        if self.voice.gas().pressure.to_bits() != exterior.pressure.to_bits() {
            return Err(invalid("DG mouth and exterior require equal static pressure"));
        }
        for opening in source.apertures() {
            let dx = ((opening.center_m[0] - input.center_m[0]).abs() - 0.5 * size[0]).max(0.0);
            let dy = ((opening.center_m[1] - input.center_m[1]).abs() - 0.5 * size[1]).max(0.0);
            if !dx.is_finite() || !dy.is_finite() || fs_math::det::hypot(dx, dy) < opening.radius_m {
                return Err(invalid("DG mouth rectangle and radiating terminal must not overlap"));
            }
        }
        let receiver = BaffledPressure::rectangular_outlet_multirate_differentiating(
            size, self.config.sample_rate_hz, self.config.sample_rate_hz,
            RectangularOutletReceiver {
                center_m: input.center_m,
                position_m: outlet.config.position_m,
                cells: input.cells,
                maximum_frequency_hz: outlet.config.maximum_frequency_hz,
            }, exterior.into(),
        ).map_err(|error| GestureCompileError::Render(RenderError::Voice(
            AcousticRealizeError::Nonlinear(error),
        )))?;
        self.mouth = Some(MouthPressure { input, receiver, observer_gas: None });
        Ok(self)
    }

    /// Apply the same physical exterior gas law to ALL observed openings.
    /// Attach the optional mouth first. Terminal and mouth preparations retain
    /// their actual paths and receive the same reported numerical padding before
    /// publication, so summing them preserves physical phase and flight time.
    /// This changes only the observer, not the jet, bore or radiation reaction.
    pub fn with_observer_gas(
        mut self, propagation: &GasPropagation, gate: &CancelGate,
    ) -> Result<Self, GestureCompileError> {
        if self.poisoned || self.voice.samples_rendered() != 0 || self.observer_gas.is_some() {
            return Err(invalid("DG gas observation requires one cold attachment"));
        }
        let outlet = self.outlet.as_ref()
            .ok_or_else(|| invalid("gas propagation requires the bound exterior microphone"))?;
        let load = self.voice.dg().expect("admitted DG source").load();
        let source = load.radiation_record().expect("admitted computed exterior");
        let (mut receiver, mut report) = load.baffled_receiver_with_gas(
            self.config.sample_rate_hz, outlet.config, propagation, gate,
        ).map_err(|error| GestureCompileError::Render(RenderError::Voice(error)))?;
        let mut mouth_receiver = if let Some(mouth) = &self.mouth {
            Some(BaffledPressure::rectangular_outlet_multirate_gas_differentiating(
                mouth.input.size_m(&self.voice.geometry),
                self.config.sample_rate_hz, self.config.sample_rate_hz,
                RectangularOutletReceiver {
                    center_m: mouth.input.center_m, position_m: outlet.config.position_m,
                    cells: mouth.input.cells,
                    maximum_frequency_hz: outlet.config.maximum_frequency_hz,
                }, source.gas_record().state(), propagation, &mut || gate.is_requested(),
            ).map_err(|error| GestureCompileError::Render(RenderError::Voice(
                AcousticRealizeError::Nonlinear(error),
            )))?)
        } else { None };
        if let Some((mouth, mouth_report)) = &mut mouth_receiver {
            let latency = report.additional_latency_samples.max(mouth_report.additional_latency_samples);
            let align_error = |error| GestureCompileError::Render(RenderError::Voice(
                AcousticRealizeError::Nonlinear(error)));
            receiver.align_gas_latency(&mut report, latency).map_err(align_error)?;
            mouth.align_gas_latency(mouth_report, latency).map_err(align_error)?;
        }
        if gate.is_requested() { return Err(invalid("DG material observation cancelled")); }
        self.outlet.as_mut().expect("admitted exterior observer").receiver = receiver;
        self.observer_gas = Some(report);
        if let Some((receiver, report)) = mouth_receiver {
            let mouth = self.mouth.as_mut().expect("original physical mouth");
            mouth.receiver = receiver;
            mouth.observer_gas = Some(report);
        }
        Ok(self)
    }

    /// Source-derived gas propagation on every circular terminal path.
    pub fn observer_gas(&self) -> Option<&GasPropagationReport> { self.observer_gas.as_ref() }

    /// Same exterior law and common numerical padding on the original mouth.
    pub fn mouth_observer_gas(&self) -> Option<&GasPropagationReport> {
        self.mouth.as_ref().and_then(|mouth| mouth.observer_gas.as_ref())
    }

    /// Replace the authored mouth loss and end correction with the existing
    /// compact self reaction of this exact rectangle and exterior material.
    ///
    /// The retained passive load participates in the current jet/DG solve and
    /// owns its near-field storage and escaped energy. It is independent of the
    /// computed outlet load; mutual reaction between openings is not included.
    /// Original compact-band, quadrature, cancellation and cold-source checks
    /// apply. The flue's legacy `mouth_loss_rel` must be zero.
    pub fn with_mouth_reaction(mut self, gate: &CancelGate)
        -> Result<Self, GestureCompileError>
    {
        if self.poisoned || self.voice.samples_rendered() != 0 {
            return Err(invalid("DG mouth reaction requires a cold performance"));
        }
        let mouth = self.mouth.as_ref()
            .ok_or_else(|| invalid("DG mouth reaction requires its actual rectangular receiver"))?;
        let outlet = self.outlet.as_ref().expect("mouth requires the original outlet observer");
        let source = self.voice.dg().expect("admitted original DG voice").load()
            .radiation_record().expect("admitted computed DG microphone source");
        let size = mouth.input.size_m(&self.voice.geometry);
        let model = RectangularMouthLoad::new(
            size,
            mouth.input.cells, *source.gas_record().state(), self.voice.sample_period_s(),
            outlet.config.maximum_frequency_hz, gate,
        ).map_err(|error| GestureCompileError::Render(source_error(error)))?;
        self.voice.attach_mouth_radiation(model, size)
            .map_err(|error| GestureCompileError::Render(source_error(error)))?;
        Ok(self)
    }

    /// Bind one reciprocal mouth/terminal exterior to the current source solve.
    /// Every source pressure and flow uses the original mixed Rayleigh matrix
    /// and one causal rational memory. Both directions of the mutual coupling
    /// act in this sample; the original signed mouth flow is -Qin.
    ///
    /// The shared geometry resolution must cover the full mechanical Nyquist
    /// band even when the microphone observes a narrower band. The original
    /// jet remains an exploration model, and the baffle is not a free-field
    /// instrument body. Independent compact mouth reaction cannot also attach.
    pub fn with_mutual_mouth_reaction(mut self, gate: &CancelGate)
        -> Result<Self, GestureCompileError>
    {
        if self.poisoned || self.voice.samples_rendered() != 0 {
            return Err(invalid("joint DG mouth reaction requires a cold performance"));
        }
        let mouth = self.mouth.as_ref()
            .ok_or_else(|| invalid("joint DG mouth reaction requires its actual rectangular receiver"))?;
        let outlet = self.outlet.as_ref().expect("mouth requires the original outlet observer");
        let surface = crate::bernoulli_aperture::radiation::environment::multiport::bake::mixed::DrivenMouthSurface {
            center_m: mouth.input.center_m,
            size_m: mouth.input.size_m(&self.voice.geometry),
            cells: mouth.input.cells,
            terminal_radial_rings: outlet.config.radial_rings,
            terminal_angular_points: outlet.config.angular_points,
        };
        self.voice.attach_joint_mouth_radiation(surface, self.config.sample_rate_hz, gate)
            .map_err(|error| GestureCompileError::Render(source_error(error)))?;
        Ok(self)
    }

    /// Actual rectangular mouth receiver, including its own propagation support.
    pub fn mouth_receiver(&self) -> Option<&BaffledPressure> {
        self.mouth.as_ref().map(|mouth| &mouth.receiver)
    }

    /// Original placement and surface resolution of the optional mouth.
    pub fn mouth_radiation(&self) -> Option<FlueMouthRadiation> {
        self.mouth.as_ref().map(|mouth| mouth.input)
    }

    /// Complete propagation support of every selected exterior contribution.
    /// Numerical filter latency remains separately exposed by each receiver.
    /// Internal mean-inlet observation has no exterior propagation kernel.
    pub fn propagation_delay_samples(&self) -> Option<(usize, usize)> {
        let terminal = self.baffled_receiver()?.delay_samples;
        Some(self.mouth_receiver().map_or(terminal, |mouth| (
            terminal.0.min(mouth.delay_samples.0),
            terminal.1.max(mouth.delay_samples.1),
        )))
    }
}
