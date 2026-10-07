//! Manufactured joint coordinates for the native thermal material query.
//!
//! These exact legacy-axis spellings are coherent SI: `normal_pressure` (Pa),
//! `thickness` and `gap` (m), `torque` (N m), and `bolt_count` (dimensionless).
//! No axis aliases, pressure-to-torque relation, contact law or probability
//! distribution is inferred. The selected interface card still owns R''.
//! Categorical finish, fluid and pattern declarations remain card/state
//! provenance, not numeric surrogates for an unprovided constitutive law.

use fs_matdb::{ClaimSelection, ClaimSet, EnvelopeAnswer, PropertySupportError, QueryPoint};
use fs_qty::{Dims, QtyAny};
use crate::{InterfaceState, ProjectError};

fn invalid(detail: impl Into<String>) -> ProjectError {
    ProjectError { code: "project-interface-query", detail: detail.into(),
        hint: "provide a positive finite joint-state band in coherent SI and a card qualified on the exact declared coordinate names".into() }
}

/// Which location of the explicitly declared manufacturing band to evaluate.
#[derive(Debug, Clone, Copy)]
pub enum StatePoint { Nominal, Lower, Upper }

fn coordinate(value: QtyAny, half_width: QtyAny, dims: Dims, at: StatePoint) -> Result<f64, ProjectError> {
    let low = value.value - half_width.value;
    let high = value.value + half_width.value;
    if value.dims != dims || half_width.dims != dims
        || !value.value.is_finite() || !half_width.value.is_finite()
        || half_width.value < 0.0 || !low.is_finite() || low <= 0.0 || !high.is_finite()
    { return Err(invalid("joint-state value and half-width must form a positive finite band with the declared dimensions")); }
    Ok(match at { StatePoint::Nominal => value.value, StatePoint::Lower => low, StatePoint::Upper => high })
}

/// Form the exact physical query used by binding and native contact lowering.
///
/// The half-widths are support requirements, never silently converted to
/// random inputs. Nominal lowering uses `StatePoint::Nominal`.
pub fn query_point(state: &InterfaceState, temperature_axis: &str, temperature: f64, at: StatePoint)
    -> Result<QueryPoint, ProjectError>
{
    if temperature_axis.is_empty() || !temperature.is_finite() || temperature <= 0.0 {
        return Err(invalid("interface temperature requires a named axis and positive finite kelvin"));
    }
    let (axis, value, count) = match state {
        InterfaceState::DryContact { pressure, pressure_half_width, .. } =>
            ("normal_pressure", coordinate(*pressure, *pressure_half_width, crate::spec::dims::PRESSURE, at)?, None),
        InterfaceState::Adhesive { thickness, thickness_half_width }
        | InterfaceState::Tim { thickness, thickness_half_width } =>
            ("thickness", coordinate(*thickness, *thickness_half_width, crate::spec::dims::LENGTH, at)?, None),
        InterfaceState::GapWithFluid { gap, gap_half_width, .. } =>
            ("gap", coordinate(*gap, *gap_half_width, crate::spec::dims::LENGTH, at)?, None),
        InterfaceState::BoltedWithPattern { bolt_count, torque, torque_half_width, .. } => {
            if *bolt_count == 0 { return Err(invalid("a bolted interface must declare at least one fastener")); }
            ("torque", coordinate(*torque, *torque_half_width, Dims([2, 1, -2, 0, 0, 0]), at)?, Some(*bolt_count))
        }
    };
    if temperature_axis == axis || (count.is_some() && temperature_axis == "bolt_count") {
        return Err(invalid("temperature and manufactured-state axes must be distinct"));
    }
    let mut point = QueryPoint::new().with(temperature_axis, temperature)
        .and_then(|p| p.with(axis, value)).map_err(|e| invalid(e.to_string()))?;
    if let Some(count) = count {
        point = point.with("bolt_count", f64::from(count)).map_err(|e| invalid(e.to_string()))?;
    }
    Ok(point)
}

/// Admit the WHOLE temperature x manufacturing-band box through matdb's
/// continuous-support owner, then return the NOMINAL-state temperature ends.
/// Checking only two diagonal samples would miss an interior source conflict.
/// Returning the manufacturing corners as nominal values would incorrectly
/// turn pressure/thickness dependence into temperature dependence.
pub fn query_envelope(
    claims: &ClaimSet, property: &str, state: &InterfaceState,
    temperature_axis: &str, low_k: f64, high_k: f64, selection: ClaimSelection,
) -> Result<EnvelopeAnswer, PropertySupportError> {
    let point = |t, at| query_point(state, temperature_axis, t, at)
        .map_err(|e| PropertySupportError::InvalidEnvelope { reason: e.detail });
    let support = claims.query_envelope(property,
        &point(low_k, StatePoint::Lower)?, &point(high_k, StatePoint::Upper)?, selection)?;
    let nominal = claims.query_envelope(property,
        &point(low_k, StatePoint::Nominal)?, &point(high_k, StatePoint::Nominal)?, selection)?;
    if support.lower.receipt.selected != nominal.lower.receipt.selected {
        return Err(PropertySupportError::InvalidEnvelope {
            reason: "nominal joint state differs from the claim admitted over the complete state band".into(),
        });
    }
    Ok(nominal)
}

#[cfg(test)]
mod tests;
