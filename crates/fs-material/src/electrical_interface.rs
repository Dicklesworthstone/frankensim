//! Source-resolved, zero-thickness Ohmic interface response.
//!
//! Specific electrical contact resistance belongs to an ordered surface system,
//! not to either bulk material. This module consumes the existing interface-card
//! resolver; it neither fits a contact law nor chooses a default metal, coating,
//! pressure, temperature, surface finish, or heat-partition model.

use fs_matdb::{InterfaceSystemCard, MaterialStateId, QueryPoint, SurfaceSpec};
use fs_qty::{Dims, QuantitySpec};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedInterfaceStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_interface_state_point,
};

/// Canonical specific electrical contact resistance R'' = voltage jump / J_n.
pub const ELECTRICAL_CONTACT_RESISTANCE_PROPERTY: &str = "electrical_contact_resistance";
/// Dimensions of ohm m2, NOT bulk resistivity (ohm m) or thermal resistance.
pub const ELECTRICAL_CONTACT_RESISTANCE_DIMS: Dims = Dims([4, 1, -3, 0, -2, 0]);

/// One source-supported Ohmic contact evaluated at an explicit interface state.
/// The complete original receipt and ordered surface descriptors are retained.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedElectricalInterface {
    resolved: ResolvedInterfaceStatePoint,
    surfaces: [SurfaceSpec; 2],
    resistance_ohm_m2: f64,
}

impl ResolvedElectricalInterface {
    /// Adopt an already resolved multiphysics interface bundle without querying
    /// or choosing a source again. Other interface properties and their receipts
    /// remain in the same bundle; the Ohmic law consumes only its own coefficient.
    pub fn try_from_resolved(
        card: &InterfaceSystemCard,
        resolved: ResolvedInterfaceStatePoint,
    ) -> Result<Self, MaterialStatePointError> {
        if resolved.card_identity() != card.content_hash() {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "electrical interface requires its original ordered surface card",
            });
        }
        let resistance = resolved.property(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY)
            .ok_or_else(|| MaterialStatePointError::InvalidRequirement {
                property: ELECTRICAL_CONTACT_RESISTANCE_PROPERTY.into(),
                reason: "missing or ambiguous electrical interface resistance",
            })?;
        let resistance_ohm_m2 = resistance.value_si();
        if resistance.requirement().quantity()
            != QuantitySpec::dimensional(ELECTRICAL_CONTACT_RESISTANCE_DIMS)
            || !resistance_ohm_m2.is_finite() || resistance_ohm_m2 < 0.0
        {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "electrical interface requires finite nonnegative specific resistance in ohm m2",
            });
        }
        Ok(Self {
            resolved,
            surfaces: [card.surface_a().clone(), card.surface_b().clone()],
            resistance_ohm_m2,
        })
    }

    /// Refresh the complete original interface bundle, then rebuild the Ohmic
    /// response from its pinned source. No earlier resistance survives a state
    /// change, and failure leaves this interface unchanged.
    pub fn at_query_point(
        &self,
        card: &InterfaceSystemCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_resolved(card, self.resolved.at_query_point(card, point)?)
    }

    /// Complete source resolution, including pressure/history/temperature axes.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedInterfaceStatePoint { &self.resolved }
    /// Specific contact resistance in ohm m2. Divide by the actual contact area
    /// to obtain resistance; never add it directly to bulk resistivity.
    #[must_use]
    pub const fn resistance_ohm_m2(&self) -> f64 { self.resistance_ohm_m2 }
    /// Ordered source surface states, including their exact texture frames.
    #[must_use]
    pub const fn surfaces(&self) -> &[SurfaceSpec; 2] { &self.surfaces }

    /// Check the material assignments at the two receiving surfaces in order.
    /// Equal chemistry with a different process, phase or revision is not equal.
    /// Geometry and texture assignment remain separate caller responsibilities.
    pub fn require_surface_materials(
        &self, first: &MaterialStateId, second: &MaterialStateId,
    ) -> Result<(), MaterialStatePointError> {
        if &self.surfaces[0].material != first || &self.surfaces[1].material != second {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "electrical interface surface material assignment mismatch",
            });
        }
        Ok(())
    }
}

/// Resolve one passive Ohmic interface using the canonical source/query owner.
///
/// All source conditions are explicit in `point`. An interface-temperature model
/// (first surface, second surface, film temperature, etc.) must be declared by
/// the caller; this function never infers it from bulk temperatures. Finite
/// nonnegative R'' is admitted, including an explicitly sourced zero. Absence
/// of a claim is not perfect contact. Re-query after any state change.
///
/// Heat partition and actual contact area are independent physical inputs to
/// the consumer. This is not a breakdown, tunnelling, thermoelectric, thermal
/// contact, wear, or contact-mechanics law, and conveys no empirical validation.
///
/// # Errors
/// Missing/ambiguous claims, unsupported coordinates, wrong dimensions or
/// quantity semantics, and negative/nonfinite resistance refuse atomically.
pub fn resolve_electrical_interface_state(
    card: &InterfaceSystemCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<ResolvedElectricalInterface, MaterialStatePointError> {
    let requirement = ScalarPropertyRequirement::try_new(
        ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
        ELECTRICAL_CONTACT_RESISTANCE_DIMS,
        ScalarAdmissibility::NonNegative,
    )?;
    let resolved = resolve_interface_state_point(card, point, &[requirement], selection)?;
    ResolvedElectricalInterface::try_from_resolved(card, resolved)
}
