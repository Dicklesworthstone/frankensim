//! Complete physical receiver peaks, cancellation notches and design linewidths.

use super::fields::transfer::{StructuralInteriorRelaxationReceiverSpectrum, evaluate};
use super::{
    Error, PreparedStructuralInteriorRelaxation, PreparedStructuralInteriorRelaxationDirection,
    StructuralInteriorRelaxationSourceOrder as Order, StructuralInteriorRelaxationSourceSpectrum,
};
use core::cell::RefCell;
use fs_exec::Cx;
use fs_modal::harmonic::resonance::{
    self as query, HalfPowerBandwidth, HalfPowerBandwidthParameterResponse, ResonanceControls,
    ResonanceError, TransferExtremum, TransferExtremumParameterResponse, TransferParameterSurvey,
    TransferSurvey,
};

impl PreparedStructuralInteriorRelaxation<'_, '_, '_, '_> {
    /// Resolve an isolated peak or cancellation notch of one original field.
    /// Coherent source and receiver spectra retain their complete frequency
    /// slopes. Two state solves share one factor at each sampled frequency;
    /// all material and constant-potential coordinates remain present.
    pub fn transfer_extremum<
        S: StructuralInteriorRelaxationSourceSpectrum + ?Sized,
        R: StructuralInteriorRelaxationReceiverSpectrum + ?Sized,
    >(
        &self,
        sources: &mut S,
        receiver: &mut R,
        bracket_rad_s: (f64, f64),
        controls: ResonanceControls,
        cx: &Cx<'_>,
    ) -> Result<TransferExtremum, ResonanceError<Error>> {
        let (sources, receiver) = (RefCell::new(sources), RefCell::new(receiver));
        query::refine_transfer_extremum(bracket_rad_s, controls, cx, |w| {
            evaluate(
                &self.inner,
                &sources,
                &receiver,
                w,
                Order::Frequency,
                None,
                &mut || cx.checkpoint().is_err(),
            )
            .map(|jet| jet.sample)
        })
    }

    /// Differentiate the same simple peak or notch under complete material,
    /// geometry, source and receiver changes. The original direction is mapped
    /// once; only the final sample needs all five same-factor response orders.
    /// An exact simple transfer zero retains its frequency and complex-amplitude
    /// derivatives; phase and logarithmic amplitude are absent at that zero.
    pub fn transfer_extremum_with_parameter_derivative<
        S: StructuralInteriorRelaxationSourceSpectrum + ?Sized,
        R: StructuralInteriorRelaxationReceiverSpectrum + ?Sized,
    >(
        &self,
        sources: &mut S,
        receiver: &mut R,
        direction: &PreparedStructuralInteriorRelaxationDirection<'_, '_>,
        bracket_rad_s: (f64, f64),
        controls: ResonanceControls,
        cx: &Cx<'_>,
    ) -> Result<TransferExtremumParameterResponse, ResonanceError<Error>> {
        let (sources, receiver) = (RefCell::new(sources), RefCell::new(receiver));
        query::refine_transfer_extremum_with_parameter_derivative(
            bracket_rad_s,
            controls,
            cx,
            |w| {
                evaluate(
                    &self.inner,
                    &sources,
                    &receiver,
                    w,
                    Order::Frequency,
                    None,
                    &mut || cx.checkpoint().is_err(),
                )
                .map(|jet| jet.sample)
            },
            |w| {
                evaluate(
                    &self.inner,
                    &sources,
                    &receiver,
                    w,
                    Order::Parameter,
                    Some(&direction.inner),
                    &mut || cx.checkpoint().is_err(),
                )
            },
        )
    }

    /// Resolve the two half-|H_peak|^2 crossings of the same positive response
    /// maximum. Q is the measured receiver peak frequency divided by its actual
    /// width, not a free-pole decay or material absorption Q. The source, receiver
    /// and model must agree with the supplied peak; each flank must be isolated.
    #[allow(clippy::too_many_arguments)]
    pub fn transfer_bandwidth<
        S: StructuralInteriorRelaxationSourceSpectrum + ?Sized,
        R: StructuralInteriorRelaxationReceiverSpectrum + ?Sized,
    >(
        &self,
        sources: &mut S,
        receiver: &mut R,
        peak: &TransferExtremum,
        lower_bracket_rad_s: (f64, f64),
        upper_bracket_rad_s: (f64, f64),
        controls: ResonanceControls,
        cx: &Cx<'_>,
    ) -> Result<HalfPowerBandwidth, ResonanceError<Error>> {
        let (sources, receiver) = (RefCell::new(sources), RefCell::new(receiver));
        query::refine_half_power_bandwidth(
            peak,
            lower_bracket_rad_s,
            upper_bracket_rad_s,
            controls,
            cx,
            |w| {
                evaluate(
                    &self.inner,
                    &sources,
                    &receiver,
                    w,
                    Order::Frequency,
                    None,
                    &mut || cx.checkpoint().is_err(),
                )
                .map(|jet| jet.sample)
            },
        )
    }

    /// Continue both actual half-power crossings, linewidth and observed Q.
    /// The peak-height change, coherent source/receiver changes and every
    /// physical direction enter the crossing derivative. No Lorentzian shape
    /// or finite parameter step is inferred; flat crossings and switches refuse.
    #[allow(clippy::too_many_arguments)]
    pub fn transfer_bandwidth_with_parameter_derivative<
        S: StructuralInteriorRelaxationSourceSpectrum + ?Sized,
        R: StructuralInteriorRelaxationReceiverSpectrum + ?Sized,
    >(
        &self,
        sources: &mut S,
        receiver: &mut R,
        direction: &PreparedStructuralInteriorRelaxationDirection<'_, '_>,
        peak: &TransferExtremumParameterResponse,
        lower_bracket_rad_s: (f64, f64),
        upper_bracket_rad_s: (f64, f64),
        controls: ResonanceControls,
        cx: &Cx<'_>,
    ) -> Result<HalfPowerBandwidthParameterResponse, ResonanceError<Error>> {
        let (sources, receiver) = (RefCell::new(sources), RefCell::new(receiver));
        query::refine_half_power_bandwidth_with_parameter_derivative(
            peak,
            lower_bracket_rad_s,
            upper_bracket_rad_s,
            controls,
            cx,
            |w| {
                evaluate(
                    &self.inner,
                    &sources,
                    &receiver,
                    w,
                    Order::Frequency,
                    None,
                    &mut || cx.checkpoint().is_err(),
                )
                .map(|jet| jet.sample)
            },
            |w| {
                evaluate(
                    &self.inner,
                    &sources,
                    &receiver,
                    w,
                    Order::Parameter,
                    Some(&direction.inner),
                    &mut || cx.checkpoint().is_err(),
                )
            },
        )
    }

    /// Survey the peaks and notches resolved by the authored positive rad/s
    /// grid, with one shared evaluation and extremum ceiling. This is a driven
    /// receiver survey; unresolved pairs and a complete free-mode census are
    /// outside its claim. A failure publishes no partial survey.
    pub fn survey_transfer_extrema<
        S: StructuralInteriorRelaxationSourceSpectrum + ?Sized,
        R: StructuralInteriorRelaxationReceiverSpectrum + ?Sized,
    >(
        &self,
        sources: &mut S,
        receiver: &mut R,
        frequencies_rad_s: &[f64],
        controls: ResonanceControls,
        maximum_extrema: usize,
        cx: &Cx<'_>,
    ) -> Result<TransferSurvey, ResonanceError<Error>> {
        let (sources, receiver) = (RefCell::new(sources), RefCell::new(receiver));
        query::survey_transfer_extrema(frequencies_rad_s, controls, maximum_extrema, cx, |w| {
            evaluate(
                &self.inner,
                &sources,
                &receiver,
                w,
                Order::Frequency,
                None,
                &mut || cx.checkpoint().is_err(),
            )
            .map(|jet| jet.sample)
        })
    }

    /// Resolve and locally differentiate each grid-separated receiver feature.
    /// Cached first-order samples serve discovery/refinement; one complete
    /// original jet per returned feature is reserved within the same budget.
    #[allow(clippy::too_many_arguments)]
    pub fn survey_transfer_extrema_with_parameter_derivative<
        S: StructuralInteriorRelaxationSourceSpectrum + ?Sized,
        R: StructuralInteriorRelaxationReceiverSpectrum + ?Sized,
    >(
        &self,
        sources: &mut S,
        receiver: &mut R,
        direction: &PreparedStructuralInteriorRelaxationDirection<'_, '_>,
        frequencies_rad_s: &[f64],
        controls: ResonanceControls,
        maximum_extrema: usize,
        cx: &Cx<'_>,
    ) -> Result<TransferParameterSurvey, ResonanceError<Error>> {
        let (sources, receiver) = (RefCell::new(sources), RefCell::new(receiver));
        query::survey_transfer_extrema_with_parameter_derivative(
            frequencies_rad_s,
            controls,
            maximum_extrema,
            cx,
            |w| {
                evaluate(
                    &self.inner,
                    &sources,
                    &receiver,
                    w,
                    Order::Frequency,
                    None,
                    &mut || cx.checkpoint().is_err(),
                )
                .map(|jet| jet.sample)
            },
            |w| {
                evaluate(
                    &self.inner,
                    &sources,
                    &receiver,
                    w,
                    Order::Parameter,
                    Some(&direction.inner),
                    &mut || cx.checkpoint().is_err(),
                )
            },
        )
    }
}
