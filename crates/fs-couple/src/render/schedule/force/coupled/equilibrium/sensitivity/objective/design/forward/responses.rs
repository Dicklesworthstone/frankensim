//! Values of declared physical response constraints, not their derivatives.
use super::*;
use super::super::constraints::{ConstraintSense, ResponseConstraint, ResponseQuantity};
use fs_dcontact::ContactStorage;
use fs_phs::Storage;

/// One original, unpenalized response and its signed normalized residual.
#[derive(Clone, Debug, PartialEq)]
pub struct ForwardConstraintResult {
    /// Metres or newtons, according to the original ResponseQuantity.
    pub value: f64,
    /// (value-bound)/scale, negated for AtLeast. Inequalities require <= 0;
    /// equality acceptance needs an explicitly supplied physical tolerance.
    pub residual: f64,
}

/// Complete primal observations plus all response rows in declaration order.
/// No derivative, multiplier, confidence interval or physical validation claim.
#[derive(Clone, Debug, PartialEq)]
pub struct ForwardDesignResponses {
    /// Same observations as evaluate_forward, with zero unassessed constraints.
    pub observations: ForwardDesignEvaluation,
    /// Every authored constraint, including violated ones. Empty only when the
    /// problem declares no constraints; never a truncated or failed-case family.
    pub constraints: Vec<ForwardConstraintResult>,
}

// The same zero base Hamiltonian used by contact preload: ONLY the original
// contact potential contributes, with its original weights and exponent.
struct NoOtherEnergy;
impl Storage for NoOtherEnergy {
    fn hamiltonian(&self, _: &[f64]) -> f64 { 0.0 }
    fn gradient(&self, _: &[f64], out: &mut [f64]) { out.fill(0.0); }
}

pub(super) fn observe(constraint: &ResponseConstraint, network: &CoupledModalSystem,
    q: &[f64], contacts: &[(ModalContact, ModalContactConfig)], gate: &CancelGate)
    -> Result<ForwardConstraintResult, ModalCouplingError>
{
    poll(Some(gate))?;
    let value = match &constraint.quantity {
        ResponseQuantity::Displacement(map) => {
            let start = network.offsets[map.component];
            dot(&map.shapes, &q[start..start + map.shapes.len()])?
        }
        ResponseQuantity::SpringForce(i) => {
            let spring = &network.connections[*i];
            let extension = finite(dot(&network.columns[*i], q)? - spring.rest_extension_m)?;
            finite(-spring.stiffness_n_m * extension)?
        }
        ResponseQuantity::ContactForce(i) | ResponseQuantity::ContactPenetration(i) => {
            let (contact, config) = &contacts[*i];
            let column = contact_column(network, contact, *config)?;
            let closure = dot(&column, q)?;
            if matches!(&constraint.quantity, ResponseQuantity::ContactPenetration(_)) {
                finite(closure - contact.law.gaps()[0])?.max(0.0)
            } else {
                // Exactly the preload's potential gradient, including at a kink.
                // static_differential would incorrectly impose a derivative gate.
                let storage = ContactStorage::new(Box::new(NoOtherEnergy), 1, vec![contact.law.clone()])
                    .map_err(ModalCouplingError::ContactLaw)?;
                let mut gradient = [0.0; 2];
                storage.gradient(&[-closure, 0.0], &mut gradient);
                let force = finite(-gradient[0])?;
                if force < 0.0 { return Err(invalid("forward normal reaction cannot be attractive")); }
                force
            }
        }
    };
    let sign = if constraint.sense == ConstraintSense::AtLeast { -1.0 } else { 1.0 };
    let residual = finite(sign * finite(value - constraint.bound)? / constraint.scale)?;
    poll(Some(gate))?;
    Ok(ForwardConstraintResult { value, residual })
}
