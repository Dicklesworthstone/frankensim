//! Passive live constitutive response for the shared dissipative element.
//!
//! The physics owner resolves its material law. This boundary supplies no
//! material fit, interpolation, storage, or independent time integration.

use core::fmt;

use super::{ElementError, temperature};

/// How a damper's reference resistance follows its live temperature and flow.
///
/// The signed relative flow is `f_a - f_b` in the element's power-domain SI
/// units: m/s, rad/s, A, or the corresponding domain flow. A nonnegative
/// factor keeps `R f^2` passive, but does not require positive differential
/// resistance. The resulting effort and heat maps must be differentiable
/// where a smooth coupled solve uses them; uniqueness is not implied.
///
/// Implementations must be pure in these inputs. Material-card validity and
/// source lookup belong to the owner; an unsupported state returns its
/// reason and is never clamped or extrapolated by the element.
pub trait DissipationLaw: fmt::Debug + Send + Sync {
    /// The dimensionless factor and its partial derivatives, in the order
    /// `[d_factor/d_kelvin, d_factor/d_relative_flow]`.
    ///
    /// The factor must be finite and nonnegative; zero is an explicitly
    /// lossless state. Partials must be finite. The element multiplies the
    /// factor by its original coefficient in effort per flow. At zero flow,
    /// a bounded one-sided factor derivative is sufficient when effort and
    /// heat have unique derivatives there: `1 + abs(flow)/scale`, for
    /// example, may return a zero flow partial at the origin.
    ///
    /// # Errors
    /// The source law's reason for refusing this temperature or flow.
    fn factor(&self, kelvin: f64, relative_flow: f64) -> Result<(f64, [f64; 2]), String>;

    /// Complete stable identity of the source law, its parameters and support.
    /// These words join the element's durable problem identity.
    fn identity_words(&self, words: &mut Vec<u64>);
}

/// One evaluated law and the exact derivatives used by both tangent directions.
pub(super) struct Response {
    pub effort: f64,
    pub effort_by_temperature: f64,
    pub effort_by_flow: f64,
    pub power_by_temperature: f64,
    pub power_by_flow: f64,
}

/// Resolve the live constitutive response before any caller output is written.
pub(super) fn response(
    reference: f64,
    law: &dyn DissipationLaw,
    kelvin: f64,
    relative_flow: f64,
) -> Result<Response, ElementError> {
    let kelvin = temperature(kelvin)?;
    if !relative_flow.is_finite() {
        return Err(ElementError::Input {
            value: relative_flow,
        });
    }
    let (factor, [by_temperature, by_flow]) = law
        .factor(kelvin, relative_flow)
        .map_err(|reason| ElementError::Law { reason })?;
    if !(factor.is_finite() && factor >= 0.0 && by_temperature.is_finite() && by_flow.is_finite()) {
        return Err(ElementError::Law {
            reason: "dissipation factor must be nonnegative and finite, with finite partials"
                .into(),
        });
    }
    let resistance = reference * factor;
    if !resistance.is_finite() || (factor > 0.0 && resistance == 0.0) {
        return Err(ElementError::Law {
            reason: "scaled resistance is not representable".into(),
        });
    }
    let effort = resistance * relative_flow;
    let power = effort * relative_flow;
    // F = R0 factor(T, f) f; P = F f. The heat output is -P, so
    // its tangent carries the negative of the two power partials below.
    let effort_by_temperature = reference * (relative_flow * by_temperature);
    let effort_by_flow = reference * (factor + relative_flow * by_flow);
    let power_by_temperature = effort_by_temperature * relative_flow;
    let power_by_flow = effort + effort_by_flow * relative_flow;
    if ![
        effort,
        power,
        effort_by_temperature,
        effort_by_flow,
        power_by_temperature,
        power_by_flow,
    ]
    .iter()
    .all(|value| value.is_finite())
    {
        return Err(ElementError::Law {
            reason: "dissipative effort, power or tangent is not finite".into(),
        });
    }
    Ok(Response {
        effort,
        effort_by_temperature,
        effort_by_flow,
        power_by_temperature,
        power_by_flow,
    })
}
