//! Actual power resonances of original and wave-loaded port-Hamiltonian models.
//!
//! This composition binds complete native physics to the shared signed-power
//! query owner. It does not add a root algorithm, fit a mode, or infer material
//! heat from a wave deficit. Every trial closes physical feedback before one
//! factorization; only final extrema/edge calls request higher response jets.

use std::cell::RefCell;

use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{self as query, ResonanceError};
use fs_phs::harmonic::scattering::ScatteringModel;
use fs_phs::harmonic::{HarmonicError, HarmonicModel, HarmonicParameterDirection};

pub use fs_modal::harmonic::resonance::{
    HalfPowerBandwidth, HalfPowerBandwidthParameterResponse, HarmonicPowerCurvatureSample,
    HarmonicPowerParameterJet, HarmonicPowerSample, PowerExtremum, PowerExtremumParameterResponse,
    PowerParameterSurvey, PowerSurvey, ResonanceControls,
};

mod observe;

/// Original physical preparation. Scattering has already closed all selected
/// references while retaining every original state, constraint and memory arm.
#[derive(Clone, Copy)]
pub enum PhsPowerModel<'m> {
    Original(&'m HarmonicModel),
    Scattering(&'m ScatteringModel),
}

/// Observe the original signed quantity itself, in watts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhsPowerObservation {
    /// Independent original R/T work. For a linear passive image this is heat;
    /// for an equilibrium linearization it retains signed incremental semantics.
    /// Selected reference loading is never counted as material dissipation.
    Dissipation,
    /// Aggregate original conjugate input/output work, including actual loaded
    /// inputs rather than effective incident-wave forcing.
    Supplied,
    /// Signed physical work at one ORIGINAL input/output slot.
    Port(usize),
    /// Incident flux at one selected port in ScatteringModel::wave_ports order.
    IncidentWave(usize),
    /// Outgoing flux at one selected port in ScatteringModel::wave_ports order.
    OutgoingWave(usize),
}

/// Derivative orders requested from the same coherent original source law.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhsPowerSourceOrder {
    Frequency,
    Curvature,
    Parameter,
}

/// Original input jets. Selected scattering slots contain incident output waves;
/// unselected slots contain ordinary conjugate inputs. None authors zero.
#[derive(Clone, Copy, Debug)]
pub struct PhsPowerSourceJet<'s> {
    pub input: &'s [C64],
    pub angular_derivative: Option<&'s [C64]>,
    pub angular_second_derivative: Option<&'s [C64]>,
    pub parameter_derivative: Option<&'s [C64]>,
    pub angular_parameter_derivative: Option<&'s [C64]>,
}

impl<'s> PhsPowerSourceJet<'s> {
    /// Frequency- and parameter-independent original source phasors.
    pub fn fixed(input: &'s [C64]) -> Self {
        Self {
            input,
            angular_derivative: None,
            angular_second_derivative: None,
            parameter_derivative: None,
            angular_parameter_derivative: None,
        }
    }
}

/// A lending source law may reuse its own buffers at every frequency.
///
/// All orders must describe the SAME continuous source and parameter. A higher
/// order may calculate additional derivatives, but must reproduce its input
/// and first tangent. Angular derivatives are per rad/s; parameter derivatives
/// use the same physical parameter as PhsPowerParameterDirection. Source-owned
/// computation/storage limits remain with the law. The query caps call count.
pub trait PhsPowerSource {
    fn jet(
        &mut self,
        angular_frequency_rad_s: f64,
        order: PhsPowerSourceOrder,
    ) -> Result<PhsPowerSourceJet<'_>, HarmonicError>;
}

impl PhsPowerSource for &[C64] {
    fn jet(
        &mut self,
        _: f64,
        _: PhsPowerSourceOrder,
    ) -> Result<PhsPowerSourceJet<'_>, HarmonicError> {
        Ok(PhsPowerSourceJet::fixed(self))
    }
}

/// One local original physical direction, reused throughout a query.
#[derive(Clone, Copy, Debug, Default)]
pub struct PhsPowerParameterDirection<'d> {
    pub physical: HarmonicParameterDirection<'d>,
    /// dZ/dtheta in selected wave-port order. References remain independent of
    /// frequency; dispersive physical loads belong in the native memory states.
    pub references: Option<&'d [f64]>,
}

/// Source-independent binding of one physical model and observed power channel.
///
/// Native preparations/factors retain their own memory, work and residual caps.
/// The shared query controls bound the complete refinement or survey call count.
/// Returned bandwidth Q describes this observed positive power maximum, not a
/// free-mode decay, material loss factor or assumed Lorentzian line shape.
pub struct PhsPowerQuery<'m> {
    model: PhsPowerModel<'m>,
    observation: PhsPowerObservation,
}

impl<'m> PhsPowerQuery<'m> {
    pub fn new(
        model: PhsPowerModel<'m>,
        observation: PhsPowerObservation,
    ) -> Result<Self, HarmonicError> {
        if let PhsPowerObservation::IncidentWave(p) | PhsPowerObservation::OutgoingWave(p) =
            observation
        {
            match model {
                PhsPowerModel::Scattering(m) if p < m.wave_ports().len() => {}
                _ => {
                    return Err(invalid(
                        "wave power requires an admitted selected wave port",
                    ));
                }
            }
        }
        Ok(Self { model, observation })
    }

    /// Sample signed P and its complete frequency tangent with two native solves.
    pub fn sample(
        &self,
        omega: f64,
        source: &mut (impl PhsPowerSource + ?Sized),
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicPowerSample, HarmonicError> {
        poll(&mut cancelled)?;
        let input = source.jet(omega, PhsPowerSourceOrder::Frequency)?;
        poll(&mut cancelled)?;
        match self.model {
            PhsPowerModel::Original(model) => observe::original_sample(
                self.observation,
                &model
                    .factor(omega, &mut cancelled)?
                    .response_with_frequency_derivative(
                        input.input,
                        input.angular_derivative,
                        &mut cancelled,
                    )?,
            ),
            PhsPowerModel::Scattering(model) => observe::loaded_sample(
                self.observation,
                &model
                    .factor(omega, &mut cancelled)?
                    .response_with_frequency_derivative(
                        input.input,
                        input.angular_derivative,
                        &mut cancelled,
                    )?,
            ),
        }
    }

    /// Complete frequency curvature with three native solves and one factor.
    pub fn curvature(
        &self,
        omega: f64,
        source: &mut (impl PhsPowerSource + ?Sized),
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicPowerCurvatureSample, HarmonicError> {
        poll(&mut cancelled)?;
        let input = source.jet(omega, PhsPowerSourceOrder::Curvature)?;
        poll(&mut cancelled)?;
        match self.model {
            PhsPowerModel::Original(model) => observe::original_curvature(
                self.observation,
                &model
                    .factor(omega, &mut cancelled)?
                    .response_with_frequency_curvature(
                        input.input,
                        input.angular_derivative,
                        input.angular_second_derivative,
                        &mut cancelled,
                    )?,
            ),
            PhsPowerModel::Scattering(model) => observe::loaded_curvature(
                self.observation,
                &model
                    .factor(omega, &mut cancelled)?
                    .response_with_frequency_curvature(
                        input.input,
                        input.angular_derivative,
                        input.angular_second_derivative,
                        &mut cancelled,
                    )?,
            ),
        }
    }

    /// Complete original five-order power jet with five solves on one factor.
    /// Original J/R/Q/G/T and changing-reference terms retain native semantics.
    pub fn parameter_jet(
        &self,
        omega: f64,
        source: &mut (impl PhsPowerSource + ?Sized),
        direction: PhsPowerParameterDirection<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicPowerParameterJet, HarmonicError> {
        poll(&mut cancelled)?;
        self.check_direction(direction)?;
        let input = source.jet(omega, PhsPowerSourceOrder::Parameter)?;
        poll(&mut cancelled)?;
        match self.model {
            PhsPowerModel::Original(model) => observe::original_parameter(
                self.observation,
                &model
                    .factor(omega, &mut cancelled)?
                    .response_with_frequency_parameter_derivative(
                        input.input,
                        direction.physical,
                        input.angular_derivative,
                        input.angular_second_derivative,
                        input.parameter_derivative,
                        input.angular_parameter_derivative,
                        &mut cancelled,
                    )?,
            ),
            PhsPowerModel::Scattering(model) => observe::loaded_parameter(
                self.observation,
                &model
                    .factor(omega, &mut cancelled)?
                    .response_with_frequency_parameter_derivative(
                        input.input,
                        direction.physical,
                        input.angular_derivative,
                        input.angular_second_derivative,
                        input.parameter_derivative,
                        input.angular_parameter_derivative,
                        direction.references,
                        &mut cancelled,
                    )?,
            ),
        }
    }

    /// Resolve one isolated signed-power extremum. The higher source/physics
    /// order is requested only once, at the returned final frequency.
    pub fn extremum(
        &self,
        bracket_rad_s: (f64, f64),
        source: &mut (impl PhsPowerSource + ?Sized),
        controls: ResonanceControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PowerExtremum, ResonanceError<HarmonicError>> {
        let source = RefCell::new(source);
        let cancel = RefCell::new(cancelled);
        query::refine_power_extremum_with_cancel(
            bracket_rad_s,
            controls,
            || (*cancel.borrow_mut())(),
            |w| self.sample(w, &mut **source.borrow_mut(), &mut *cancel.borrow_mut()),
            |w| self.curvature(w, &mut **source.borrow_mut(), &mut *cancel.borrow_mut()),
        )
    }

    /// Local continuation of this same physical extremum, including source-law,
    /// material and reference directions. No perturbed model is factored.
    pub fn extremum_with_parameter_derivative(
        &self,
        bracket_rad_s: (f64, f64),
        source: &mut (impl PhsPowerSource + ?Sized),
        direction: PhsPowerParameterDirection<'_>,
        controls: ResonanceControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PowerExtremumParameterResponse, ResonanceError<HarmonicError>> {
        self.check_direction(direction)
            .map_err(ResonanceError::Evaluation)?;
        let source = RefCell::new(source);
        let cancel = RefCell::new(cancelled);
        query::refine_power_extremum_with_parameter_derivative_with_cancel(
            bracket_rad_s,
            controls,
            || (*cancel.borrow_mut())(),
            |w| self.sample(w, &mut **source.borrow_mut(), &mut *cancel.borrow_mut()),
            |w| {
                self.parameter_jet(
                    w,
                    &mut **source.borrow_mut(),
                    direction,
                    &mut *cancel.borrow_mut(),
                )
            },
        )
    }

    /// Actual half-power crossings around the same positive maximum. Brackets
    /// isolate its two flanks; the supplied peak must come from this source,
    /// model and channel. Signed or absent crossings refuse, without abs(P).
    pub fn bandwidth(
        &self,
        peak: &PowerExtremum,
        lower: (f64, f64),
        upper: (f64, f64),
        source: &mut (impl PhsPowerSource + ?Sized),
        controls: ResonanceControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HalfPowerBandwidth, ResonanceError<HarmonicError>> {
        let cancel = RefCell::new(cancelled);
        query::refine_power_bandwidth_with_cancel(
            peak,
            lower,
            upper,
            controls,
            || (*cancel.borrow_mut())(),
            |w| self.sample(w, source, &mut *cancel.borrow_mut()),
        )
    }

    /// Both half-power edge motions, measured width and Q for the same direction.
    /// Higher native jets are evaluated only at the two final crossing samples.
    #[allow(clippy::too_many_arguments)]
    pub fn bandwidth_with_parameter_derivative(
        &self,
        peak: &PowerExtremumParameterResponse,
        lower: (f64, f64),
        upper: (f64, f64),
        source: &mut (impl PhsPowerSource + ?Sized),
        direction: PhsPowerParameterDirection<'_>,
        controls: ResonanceControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HalfPowerBandwidthParameterResponse, ResonanceError<HarmonicError>> {
        self.check_direction(direction)
            .map_err(ResonanceError::Evaluation)?;
        let source = RefCell::new(source);
        let cancel = RefCell::new(cancelled);
        query::refine_power_bandwidth_with_parameter_derivative_with_cancel(
            peak,
            lower,
            upper,
            controls,
            || (*cancel.borrow_mut())(),
            |w| self.sample(w, &mut **source.borrow_mut(), &mut *cancel.borrow_mut()),
            |w| {
                self.parameter_jet(
                    w,
                    &mut **source.borrow_mut(),
                    direction,
                    &mut *cancel.borrow_mut(),
                )
            },
        )
    }

    /// Resolve both maxima and minima supported by the authored frequency grid.
    /// Every native trial retains all states; grid resolution is not a claim of
    /// spectral completeness. Global call and extremum caps refuse truncation.
    pub fn survey(
        &self,
        frequencies_rad_s: &[f64],
        source: &mut (impl PhsPowerSource + ?Sized),
        controls: ResonanceControls,
        maximum_extrema: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PowerSurvey, ResonanceError<HarmonicError>> {
        let source = RefCell::new(source);
        let cancel = RefCell::new(cancelled);
        query::survey_power_extrema_with_cancel(
            frequencies_rad_s,
            controls,
            maximum_extrema,
            || (*cancel.borrow_mut())(),
            |w| self.sample(w, &mut **source.borrow_mut(), &mut *cancel.borrow_mut()),
            |w| self.curvature(w, &mut **source.borrow_mut(), &mut *cancel.borrow_mut()),
        )
    }

    /// Local physical continuation of every grid-resolved signed-power extremum.
    /// One final full jet per feature; first-order solves serve all grid/refinement trials.
    #[allow(clippy::too_many_arguments)]
    pub fn survey_with_parameter_derivative(
        &self,
        frequencies_rad_s: &[f64],
        source: &mut (impl PhsPowerSource + ?Sized),
        direction: PhsPowerParameterDirection<'_>,
        controls: ResonanceControls,
        maximum_extrema: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PowerParameterSurvey, ResonanceError<HarmonicError>> {
        self.check_direction(direction)
            .map_err(ResonanceError::Evaluation)?;
        let source = RefCell::new(source);
        let cancel = RefCell::new(cancelled);
        query::survey_power_extrema_with_parameter_derivative_with_cancel(
            frequencies_rad_s,
            controls,
            maximum_extrema,
            || (*cancel.borrow_mut())(),
            |w| self.sample(w, &mut **source.borrow_mut(), &mut *cancel.borrow_mut()),
            |w| {
                self.parameter_jet(
                    w,
                    &mut **source.borrow_mut(),
                    direction,
                    &mut *cancel.borrow_mut(),
                )
            },
        )
    }

    fn check_direction(
        &self,
        direction: PhsPowerParameterDirection<'_>,
    ) -> Result<(), HarmonicError> {
        if let Some(reference) = direction.references {
            match self.model {
                PhsPowerModel::Scattering(m)
                    if reference.len() == m.wave_ports().len()
                        && reference.iter().all(|v| v.is_finite()) => {}
                _ => {
                    return Err(invalid(
                        "reference direction must follow the selected loaded ports",
                    ));
                }
            }
        }
        Ok(())
    }
}

fn invalid(what: &'static str) -> HarmonicError {
    HarmonicError::Invalid { what }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), HarmonicError> {
    if cancelled() {
        Err(HarmonicError::Cancelled)
    } else {
        Ok(())
    }
}
