//! Uniform, source-bound Navier slip on the original spherical boundary.
//!
//! Slip length is an interface property, independent of bulk viscosity and
//! body density. The boundary is impermeable; tangential velocity jump equals
//! b/mu times tangential viscous traction, including spherical curvature.
//! This is continuum interfacial slip, not a gas rarefaction correction.

mod resistance;
pub(super) use resistance::ResistanceFactors;

use super::CreepingFlowError as Error;
use fs_matdb::EvaluationDecision;
use fs_material::state_point::{ResolvedInterfaceStatePoint, ResolvedMaterialStatePoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Length, QuantitySpec, Temperature};

/// Physical Navier slip length [m] on the ordered solid/fluid interface.
pub const NAVIER_SLIP_LENGTH_PROPERTY: &str = "navier_slip_length";

/// Borrowed original interface source for a homogeneous spherical boundary.
///
/// Surface A is the sphere's solid; surface B is the surrounding fluid. The
/// original resolver owns source selection, units, validity and uncertainty.
/// The sphere checks the solid assignment; the caller must assign surface B
/// to its actual `FluidTransport` source and hold other external conditions.
/// No-slip remains an explicit boundary choice when this view is absent;
/// missing slip evidence is never resolved as a zero-length measurement.
#[derive(Clone, Copy, Debug)]
pub struct NavierSlipBoundary<'a> {
    source: &'a ResolvedInterfaceStatePoint,
    slip_length_m: f64,
}

impl<'a> NavierSlipBoundary<'a> {
    /// Consume an already-resolved interface bundle without re-querying it.
    /// Resolve `navier_slip_length` with `Length::DIMS` and nonnegative scalar
    /// admissibility using the existing interface-state resolver. The selected
    /// source must be a constant within its validity, with only fixed external
    /// conditions. Rate-, stress-, history- or temperature-dependent slip needs
    /// an evolving interface owner and is refused here.
    ///
    /// The canonical surface-temperature coordinate is `T` [absolute K]. A
    /// bounded constant remains usable during thermal expansion only while
    /// the actual solid temperature stays inside that original source bound.
    pub fn new(source: &'a ResolvedInterfaceStatePoint) -> Result<Self, Error> {
        let property = source
            .property(NAVIER_SLIP_LENGTH_PROPERTY)
            .ok_or(Error::Invalid(
                "Navier sphere requires an unambiguous source-resolved interface slip length",
            ))?;
        let slip_length_m = property.value_si();
        if property.requirement().quantity() != QuantitySpec::dimensional(Length::DIMS)
            || !slip_length_m.is_finite()
            || slip_length_m < 0.0
            || !property.source_axes().is_empty()
            || property.answer().receipt.decision != EvaluationDecision::ConstantWithinValidity
            || !property.has_only_fixed_external_conditions(Some("T"), false)
        {
            return Err(Error::Invalid(
                "Navier sphere needs finite nonnegative metres and constant static interface support",
            ));
        }
        Ok(Self {
            source,
            slip_length_m,
        })
    }

    /// Complete original interface identity, ordered materials and receipts.
    pub const fn source(&self) -> &'a ResolvedInterfaceStatePoint {
        self.source
    }

    /// Fixed physical slip length [m], never a dimensionless drag multiplier.
    pub const fn slip_length_m(&self) -> f64 {
        self.slip_length_m
    }

    pub(super) fn require_solid(&self, solid: &ResolvedMaterialStatePoint) -> Result<(), Error> {
        if &self.source.surface_materials()[0] != solid.material() {
            return Err(Error::Invalid(
                "Navier interface surface A must be the sphere's original solid",
            ));
        }
        let property = self
            .source
            .property(NAVIER_SLIP_LENGTH_PROPERTY)
            .expect("immutable admitted slip source");
        if let Some((lower, upper)) = property.answer().evidence.model.validity.bound("T") {
            let absolute = QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            ));
            let wrong_quantity = solid.axis_quantities().get("T").is_some_and(|quantity| {
                *quantity != absolute && *quantity != QuantitySpec::dimensional(Temperature::DIMS)
            });
            if wrong_quantity
                || !solid.state_coordinate("T").is_some_and(|temperature| {
                    temperature.is_finite()
                        && temperature > 0.0
                        && temperature >= lower
                        && temperature <= upper
                })
            {
                return Err(Error::Invalid(
                    "sphere surface temperature needs absolute kelvin within its Navier slip source support",
                ));
            }
        }
        Ok(())
    }
}
