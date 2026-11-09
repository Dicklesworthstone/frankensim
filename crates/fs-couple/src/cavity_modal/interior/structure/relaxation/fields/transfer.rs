//! Original coherent field observations for complete material-memory queries.

use super::super::sources::{
    StructuralInteriorRelaxationSourceOrder as Order, StructuralInteriorRelaxationSourceSpectrum,
};
use super::{
    Error, Prepared, PreparedDirection, StructuralInteriorRelaxationFrequencyResponse, poll,
};
use crate::cavity_modal::interior::structure::resonance::StructuralInteriorResponseField;
use core::cell::RefCell;
use fs_math::c64::C64;
use fs_modal::harmonic::HarmonicResponseError;
use fs_modal::harmonic::resonance::{
    HarmonicReceiverParameterJet, HarmonicTransferParameterJet, HarmonicTransferSample,
};

/// One physical field and its coherent receiver law at the current frequency.
/// Weights form a linear transfer row, without complex conjugation. Original
/// field units and coordinate order survive; pressure potential is not pressure.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorRelaxationReceiverJet<'r> {
    pub field: StructuralInteriorResponseField,
    pub receiver: HarmonicReceiverParameterJet<'r>,
}

/// A receiver may vary with frequency and the same physical design parameter.
///
/// Supply every source-independent receiver term, including the mixed one for
/// a moving sensor or changing propagation law. A provider can reuse its own
/// buffers. None declares a zero derivative, never an unknown derivative.
/// All requested orders must describe one smooth observation on the interval;
/// a final complete jet must reproduce the original first-order sample.
/// Domain, work and storage bounds of the authored receiver belong to its
/// provider. The physical query polls cancellation around every provider call.
pub trait StructuralInteriorRelaxationReceiverSpectrum {
    fn jet(
        &mut self,
        angular_frequency_rad_s: f64,
        order: Order,
    ) -> Result<StructuralInteriorRelaxationReceiverJet<'_>, Error>;
}

/// Fixed complex observation in the original field basis, with no allocation.
/// For a pressure point, weights are the original pressure-basis evaluation
/// row. A geometry-dependent row instead supplies its full receiver jet.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorRelaxationFixedReceiver<'r> {
    pub field: StructuralInteriorResponseField,
    pub weights: &'r [C64],
}

impl StructuralInteriorRelaxationReceiverSpectrum
    for StructuralInteriorRelaxationFixedReceiver<'_>
{
    fn jet(
        &mut self,
        _angular_frequency_rad_s: f64,
        _order: Order,
    ) -> Result<StructuralInteriorRelaxationReceiverJet<'_>, Error> {
        Ok(StructuralInteriorRelaxationReceiverJet {
            field: self.field,
            receiver: HarmonicReceiverParameterJet {
                weights: self.weights,
                angular_derivative: None,
                angular_second_derivative: None,
                parameter_derivative: None,
                angular_parameter_derivative: None,
            },
        })
    }
}

/// The shared numerical query owns its callback ceiling. This physical
/// evaluator owns the complete same-factor response and authored observations.
/// Only Frequency and Parameter are needed by transfer searches.
#[allow(clippy::too_many_arguments)]
pub(in super::super) fn evaluate<
    S: StructuralInteriorRelaxationSourceSpectrum + ?Sized,
    R: StructuralInteriorRelaxationReceiverSpectrum + ?Sized,
>(
    prepared: &Prepared<'_>,
    sources: &RefCell<&mut S>,
    receivers: &RefCell<&mut R>,
    omega: f64,
    order: Order,
    direction: Option<&PreparedDirection<'_, '_>>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HarmonicTransferParameterJet, Error> {
    poll(cancelled)?;
    let mut source_owner = sources.borrow_mut();
    let source = source_owner.jet(omega, order)?;
    poll(cancelled)?;
    let mut receiver_owner = receivers.borrow_mut();
    let receiver = receiver_owner.jet(omega, order)?;
    poll(cancelled)?;
    let result = match order {
        Order::Frequency => {
            let response = prepared
                .factor(omega, cancelled)?
                .solve_with_frequency_derivative(
                    source.values.structural_force,
                    source.values.inward_flow_load,
                    source
                        .angular_derivative
                        .map(|d| (d.structural_force, d.inward_flow_load)),
                    cancelled,
                )?;
            let (values, first) = frequency_field(&response, receiver.field);
            let sample = match receiver.receiver.angular_derivative {
                Some(slope) => HarmonicTransferSample::project_complex_with_derivative(
                    receiver.receiver.weights,
                    slope,
                    values,
                    first,
                ),
                None => HarmonicTransferSample::project_complex(
                    receiver.receiver.weights,
                    values,
                    first,
                ),
            }
            .map_err(observation)?;
            HarmonicTransferParameterJet {
                sample,
                angular_second_derivative: C64::ZERO,
                parameter_derivative: C64::ZERO,
                angular_parameter_derivative: C64::ZERO,
            }
        }
        Order::Parameter => {
            let direction = direction.ok_or(Error::Invalid(
                "a prepared original direction is required for a complete receiver query",
            ))?;
            let response = prepared
                .factor(omega, cancelled)?
                .solve_with_frequency_parameter_derivative(source, direction, cancelled)?;
            HarmonicTransferParameterJet::project_complex_with_derivatives(
                receiver.receiver,
                response.field_parameter_jet(receiver.field),
            )
            .map_err(observation)?
        }
        Order::Curvature => {
            return Err(Error::Invalid(
                "receiver queries request first-order samples or complete parameter jets",
            ));
        }
    };
    poll(cancelled)?;
    Ok(result)
}

fn frequency_field(
    result: &StructuralInteriorRelaxationFrequencyResponse,
    field: StructuralInteriorResponseField,
) -> (&[C64], &[C64]) {
    let original = &result.response;
    match field {
        StructuralInteriorResponseField::StructuralDisplacement => (
            &original.structural_displacement,
            &result.structural_displacement_derivative_per_rad_s,
        ),
        StructuralInteriorResponseField::StructuralVelocity => (
            &original.structural_velocity,
            &result.structural_velocity_derivative_per_rad_s,
        ),
        StructuralInteriorResponseField::PressurePotential => (
            &original.pressure_potential_pa_s,
            &result.pressure_potential_derivative_pa_s_per_rad_s,
        ),
        StructuralInteriorResponseField::Pressure => (
            &original.pressure_pa,
            &result.pressure_derivative_pa_per_rad_s,
        ),
        StructuralInteriorResponseField::FluidForce => (
            &original.fluid_force,
            &result.fluid_force_derivative_per_rad_s,
        ),
        StructuralInteriorResponseField::DisplacedFlow => (
            &original.displaced_flow_m3_s,
            &result.displaced_flow_derivative_m3_s_per_rad_s,
        ),
    }
}

fn observation(error: HarmonicResponseError) -> Error {
    match error {
        HarmonicResponseError::Cancelled => Error::Cancelled,
        _ => Error::Invalid("invalid or unrepresentable original coherent receiver jet"),
    }
}
