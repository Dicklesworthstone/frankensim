//! Physical response constraints over the SAME settled state and tangent factor.
//! Violated design constraints are finite optimizer residuals, not simulation
//! failures. Original physical ceilings and contact switching exclusions remain
//! hard refusals. Every extra row uses an analytic adjoint, never perturbations.
use super::*;

/// A scalar response of one solved load case, with its declared physical units.
#[derive(Clone, Debug)]
pub enum ResponseQuantity {
    /// Signed displacement of a fixed mass-normalized attachment [m].
    Displacement(ModalAttachment),
    /// Signed spring reaction on its LEFT side, -k*(left-right-rest) [N].
    /// Static dashpots carry no force.
    SpringForce(usize),
    /// Compressive normal reaction, nonnegative [N].
    ContactForce(usize),
    /// Positive contact penetration, zero when separated [m].
    ContactPenetration(usize),
}

/// The residual supplied to the optimizer, without clipping or squaring.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConstraintSense {
    /// (response - bound)/scale <= 0.
    AtMost,
    /// (bound - response)/scale <= 0.
    AtLeast,
    /// (response - bound)/scale == 0.
    Equal,
}

/// A physical requirement on one independent experiment, not a penalty term.
#[derive(Clone, Debug)]
pub struct ResponseConstraint {
    /// Unique nonempty name, at most 128 bytes.
    pub name: String,
    /// Index of the load case in construction order.
    pub case: usize,
    /// Observation or reaction in that case's actual equilibrium.
    pub quantity: ResponseQuantity,
    /// Upper bound, lower bound or equality.
    pub sense: ConstraintSense,
    /// Bound in metres or newtons, according to quantity.
    pub bound: f64,
    /// Strictly positive normalization in the SAME units as the bound.
    /// This is not a physical tolerance or an estimated noise level.
    pub scale: f64,
}

/// Complete value/Jacobian row, in the original constraint declaration order.
#[derive(Clone, Debug, PartialEq)]
pub struct ResponseConstraintResult {
    /// Unnormalized observation [m] or reaction [N].
    pub value: f64,
    /// Signed dimensionless residual; inequalities are feasible at <= 0.
    pub residual: f64,
    /// d(residual)/d(decision), including explicit constitutive parameter terms.
    pub gradient: Vec<f64>,
    /// Recomputed full-coordinate residual of this row's adjoint solve.
    pub adjoint_relative_residual: f64,
}

impl EquilibriumDesign {
    /// Attach a complete immutable family of physical constraints. Each may
    /// require one extra adjoint per evaluation, reusing its case's primal and
    /// prepared tangent. Count is explicitly bounded and hard-capped at 64;
    /// every solve/pullback retains the original sensitivity query work cap.
    /// No objective penalty is added and no violation is converted to an error.
    /// Consuming self prevents changing requirements on an existing study.
    pub fn with_constraints(mut self, constraints: Vec<ResponseConstraint>, maximum_constraints: usize,
        gate: &CancelGate) -> Result<Self, DesignError>
    {
        checkpoint(gate)?;
        if maximum_constraints > 64 || constraints.len() > maximum_constraints || !self.constraints.is_empty() {
            return Err(bad("physical constraints need one bounded complete family, at most 64 rows"));
        }
        for (i, constraint) in constraints.iter().enumerate() {
            checkpoint(gate)?;
            if !name_ok(&constraint.name) || constraints[..i].iter().any(|c| c.name == constraint.name)
                || constraint.case >= self.cases.len() || !constraint.bound.is_finite()
                || !constraint.scale.is_finite() || constraint.scale <= 0.0 {
                return Err(bad("physical constraint needs a unique name, existing case, finite bound and positive scale"));
            }
            match &constraint.quantity {
                ResponseQuantity::Displacement(map) => {
                    let model = self.models.get(map.component).ok_or_else(|| bad("constraint names an unknown component"))?;
                    if map.shapes.len() != model.modes().len() || map.shapes.iter().any(|x| !x.is_finite()) {
                        return Err(bad("constraint displacement map must match its complete finite modal basis"));
                    }
                }
                ResponseQuantity::SpringForce(j) if *j < self.springs.len() => {}
                ResponseQuantity::ContactForce(j) | ResponseQuantity::ContactPenetration(j)
                    if *j < self.contacts.len() => {}
                _ => return Err(bad("constraint names an unknown spring or normal contact")),
            }
        }
        self.constraints = constraints;
        checkpoint(gate)?;
        Ok(self)
    }

    /// Frozen physical requirements. Parameter box bounds remain separate.
    #[must_use]
    pub fn constraints(&self) -> &[ResponseConstraint] { &self.constraints }
}

pub(super) fn evaluate(constraint: &ResponseConstraint, linear: &EquilibriumLinearization<'_>,
    contacts: &[(ModalContact, ModalContactConfig)], case: &DesignLoadCase, variables: &[DesignVariable],
    gate: &CancelGate) -> Result<ResponseConstraintResult, ModalCouplingError>
{
    poll(Some(gate))?;
    let mut state_gradient = vec![0.0; linear.q.len()];
    // Direct d(response)/d(parameter) terms MUST accompany the indirect state
    // derivative. A reaction constraint cannot use only -R_p^T lambda.
    let mut explicit = Vec::new();
    let value = match &constraint.quantity {
        ResponseQuantity::Displacement(map) => {
            let start = linear.network.offsets[map.component];
            state_gradient[start..start + map.shapes.len()].copy_from_slice(&map.shapes);
            dot(&map.shapes, &linear.q[start..start + map.shapes.len()])?
        }
        ResponseQuantity::SpringForce(i) => {
            let spring = &linear.network.connections[*i];
            let column = &linear.network.columns[*i];
            let x = finite(dot(column, &linear.q)? - spring.rest_extension_m)?;
            for (g, b) in state_gradient.iter_mut().zip(column) { *g = finite(-spring.stiffness_n_m*b)?; }
            explicit.push((DesignField::SpringStiffness(*i), -x));
            explicit.push((DesignField::SpringRest(*i), spring.stiffness_n_m));
            finite(-spring.stiffness_n_m*x)?
        }
        ResponseQuantity::ContactForce(i) => {
            let point = &linear.points[*i];
            for (g, b) in state_gradient.iter_mut().zip(&point.column) { *g = finite(point.tangent*b)?; }
            explicit.push((DesignField::ContactStiffness(*i), point.force_per_stiffness));
            explicit.push((DesignField::ContactGap(*i), -point.tangent));
            explicit.push((DesignField::ContactWeight(*i), point.force_per_weight));
            point.force
        }
        ResponseQuantity::ContactPenetration(i) => {
            let point = &linear.points[*i];
            let penetration = finite(dot(&point.column, &linear.q)? - contacts[*i].0.law.gaps()[0])?;
            // The linearization already refuses at or near the activity switch.
            if penetration > 0.0 {
                state_gradient.copy_from_slice(&point.column);
                explicit.push((DesignField::ContactGap(*i), -1.0));
                penetration
            } else { 0.0 }
        }
    };
    let adjoint = linear.solve(&state_gradient, gate)?;
    let pullback = linear.parameter_pullback(&adjoint.values, gate)?;
    let mut physical_forces = Vec::with_capacity(case.loads.len());
    for load in &case.loads {
        poll(Some(gate))?;
        let start = linear.network.offsets[load.attachment.component];
        physical_forces.push(dot(&load.attachment.shapes, &adjoint.values[start..start + load.attachment.shapes.len()])?);
    }
    let sign = if constraint.sense == ConstraintSense::AtLeast { -1.0 } else { 1.0 };
    let residual = finite(sign*finite(value - constraint.bound)?/constraint.scale)?;
    let mut gradient = Vec::with_capacity(variables.len());
    for variable in variables {
        poll(Some(gate))?;
        let mut physical = 0.0;
        for &field in &variable.fields {
            let mut derivative = field_derivative(field, constraint.case, &pullback, &physical_forces);
            for &(parameter, direct) in &explicit {
                if parameter == field { derivative = finite(derivative + direct)?; }
            }
            physical = finite(physical + derivative)?;
        }
        gradient.push(finite(sign*finite(physical*variable.scale)?/constraint.scale)?);
    }
    poll(Some(gate))?;
    Ok(ResponseConstraintResult { value, residual, gradient, adjoint_relative_residual: adjoint.relative_residual })
}
