//! Physical displacement targets and actuator gradients over one frozen solve.
//! This is an objective adapter, not another optimizer or a fit-validity claim.
use super::*;

/// A physical displacement observation at an explicitly declared attachment.
#[derive(Clone, Debug)]
pub struct DisplacementTarget {
    /// Existing component-major mass-normalized observation map [1/sqrt(kg)].
    pub attachment: ModalAttachment,
    /// Target physical displacement [m].
    pub target_m: f64,
    /// Positive residual normalization [m], supplied rather than inferred.
    pub scale_m: f64,
    /// Nonnegative dimensionless objective weight, not a confidence claim.
    pub weight: f64,
}

/// J = sum 0.5*weight*((observation-target)/scale)^2 at the frozen equilibrium.
/// All results are complete or the entire evaluation refuses. There is ONE
/// adjoint solve regardless of the number of contacts, targets or actuators.
#[derive(Clone, Debug, PartialEq)]
pub struct DisplacementObjective {
    /// Dimensionless weighted least-squares value.
    pub value: f64,
    /// Actual predictions [m], in target input order.
    pub observations_m: Vec<f64>,
    /// Total dJ/dF [1/N] for additive physical forces on the supplied actuators.
    /// The objective has no explicit dependence on those force amplitudes.
    pub physical_force_gradient: Vec<f64>,
    /// Explicit dJ/d(target_m), in target order.
    pub target_gradient: Vec<f64>,
    /// Explicit dJ/d(scale_m), in target order.
    pub scale_gradient: Vec<f64>,
    /// Explicit dJ/d(weight), in target order.
    pub weight_gradient: Vec<f64>,
    /// Explicit derivatives of each observation's own shape entries.
    /// When an observation shape is also a mechanical parameter, ADD these to
    /// the negative mechanical residual pullback; neither contribution replaces
    /// the other. Coordinates/bases are fixed, not re-solved geometric modes.
    pub observation_shape_gradient: Vec<Vec<f64>>,
    /// Full-coordinate adjoint solution and recomputed residual.
    pub adjoint: EquilibriumTangentSolution,
    /// (dR/dp)^T lambda. Total mechanical-parameter gradients are its NEGATIVE
    /// here, plus any shared observation-shape partial explicitly given above.
    pub residual_pullback: EquilibriumParameterPullback,
}

impl EquilibriumLinearization<'_> {
    /// Evaluate physical displacement targets and total physical-actuator
    /// gradients. All maps, units, weights and scales must be explicit.
    /// `max_ports` bounds targets+actuators and is hard-capped at 1024; their
    /// visits are additionally admitted against the original query work budget.
    /// No source state or parameter is changed, including on cancellation/error.
    pub fn displacement_objective(&self, targets:&[DisplacementTarget], actuators:&[ModalAttachment],
        max_ports:usize, gate:&CancelGate)->Result<DisplacementObjective,ModalCouplingError>
    {
        poll(Some(gate))?;
        let ports=targets.len().checked_add(actuators.len()).ok_or_else(||invalid("displacement objective port count overflow"))?;
        if targets.is_empty() || max_ports>1024 || ports>max_ports {
            return Err(invalid("displacement objective requires targets and a bounded complete port family"));
        }
        let width=self.network.columns.len()+self.points.len()+1;
        let terms=self.q.len().checked_mul(width+ports).and_then(|v|v.checked_add(width*width))
            .ok_or_else(||invalid("displacement objective work overflow"))?;
        if terms>self.budget.max_query_terms {return Err(invalid("displacement objective exceeds max_query_terms"));}
        for target in targets {
            self.check_attachment(&target.attachment)?;
            if !target.target_m.is_finite() || !target.scale_m.is_finite() || target.scale_m<=0.0
                || !target.weight.is_finite() || target.weight<0.0 {
                return Err(invalid("displacement targets require finite units, positive scales and nonnegative weights"));
            }
        }
        for actuator in actuators {self.check_attachment(actuator)?;}
        let mut value=0.0;let mut state_gradient=vec![0.0;self.q.len()];
        let mut observations_m=Vec::with_capacity(targets.len());
        let mut target_gradient=Vec::with_capacity(targets.len());
        let mut scale_gradient=Vec::with_capacity(targets.len());
        let mut weight_gradient=Vec::with_capacity(targets.len());
        let mut observation_shape_gradient=Vec::with_capacity(targets.len());
        for target in targets {
            poll(Some(gate))?;
            let start=self.network.offsets[target.attachment.component];
            let q=&self.q[start..start+target.attachment.shapes.len()];
            let observed=dot(&target.attachment.shapes,q)?;
            let normalized=finite(finite(observed-target.target_m)?/target.scale_m)?;
            let squared=finite(normalized*normalized)?;
            value=finite(value+0.5*target.weight*squared)?;
            let derivative=finite(target.weight*normalized/target.scale_m)?;
            for (j,shape) in target.attachment.shapes.iter().enumerate() {
                state_gradient[start+j]=finite(state_gradient[start+j]+derivative*shape)?;
            }
            observations_m.push(observed);target_gradient.push(-derivative);
            scale_gradient.push(finite(-target.weight*squared/target.scale_m)?);
            weight_gradient.push(0.5*squared);
            observation_shape_gradient.push(q.iter().map(|q|finite(derivative*q)).collect::<Result<_,_>>()?);
        }
        let adjoint=self.solve(&state_gradient,gate)?;
        let residual_pullback=self.parameter_pullback(&adjoint.values,gate)?;
        let mut physical_force_gradient=Vec::with_capacity(actuators.len());
        for actuator in actuators {
            poll(Some(gate))?;
            let start=self.network.offsets[actuator.component];
            physical_force_gradient.push(dot(&actuator.shapes,&adjoint.values[start..start+actuator.shapes.len()])?);
        }
        poll(Some(gate))?;
        Ok(DisplacementObjective {value,observations_m,physical_force_gradient,target_gradient,scale_gradient,
            weight_gradient,observation_shape_gradient,adjoint,residual_pullback})
    }
    fn check_attachment(&self,attachment:&ModalAttachment)->Result<(),ModalCouplingError> {
        let model=self.network.models.get(attachment.component).ok_or_else(||invalid("objective attachment names an unknown component"))?;
        if attachment.shapes.len()!=model.modes().len() || attachment.shapes.iter().any(|v|!v.is_finite()) {
            return Err(invalid("objective attachment must match its complete finite modal basis"));
        }
        Ok(())
    }
}
