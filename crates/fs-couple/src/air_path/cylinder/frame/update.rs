//! Prescribed exterior-material changes on the original frame/fluid history.
//!
//! The dimensionless Stokes spectral points are spatial coordinates, fixed by
//! the admitted radius and quadrature. Their masses scale with density, while
//! their resistances scale with viscosity. Changing the latter changes diffusion
//! rates, not the meaning or ordering of a fluid momentum coordinate.

use super::*;

/// A material intervention at a fixed physical time, not a mechanical step.
/// The owner prescribing the intervention owes `parameter_work_j` exactly once.
#[derive(Clone, Copy, Debug)]
pub struct FrameFluidUpdate {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Signed work at unchanged structural displacement and ALL canonical
    /// momenta: one half p^T (M_new^-1-M_old^-1) p, including fluid memory.
    pub parameter_work_j: f64,
    /// Difference of the separately evaluated Hamiltonians minus parameter
    /// work; only arithmetic roundoff remains, not unaccounted material heat.
    pub energy_residual_j: f64,
    /// Source metadata is replaced even if every physical coefficient agrees.
    pub operator_changed: bool,
}

impl<F: FluidTransport> MaterialFluidFrame<'_, F> {
    /// Bind complete replacement exterior sources without erasing motion.
    ///
    /// One source is required for EACH original region, in declaration order.
    /// Original geometry, dry modes, supports, solid damping, Gauss traces,
    /// spectral points and coordinate ordering remain fixed. The full coupled
    /// added-mass matrix and every fluid memory mass/resistance are assembled
    /// through the same original constructor, with the original inverse-mass
    /// tolerance. No eigenmode retuning or diagonal damping approximation occurs.
    ///
    /// `state` stays bit-for-bit unchanged: structural q, loaded structural
    /// momentum M*qdot, and each fluid momentum m_j*v_j. If density changes,
    /// velocities change through the new inverse mass; the signed change in
    /// kinetic storage is explicit external parameter work. A viscosity-only
    /// change preserves all inertia and instantaneous stored energy exactly.
    ///
    /// This is an externally prescribed material intervention, not a finite
    /// gas exchange, thermal, pressure-load, entrainment or phase-change solve.
    /// No background enthalpy, buoyancy or mean-flow work is inferred. Native
    /// prepared steppers borrow the old model and must be dropped before this
    /// mutable call, then prepared again for the continuing unchanged state.
    ///
    /// All source evaluation, state/operator admission, work-budget checks and
    /// cancellation happen before publication. Refusal changes neither model
    /// nor state. The original preparation budgets bound the unchanged topology.
    pub fn rebind_fluids(&mut self, state: &[f64], fluids: Vec<F>,
        maximum_absolute_parameter_work_j: f64, gate: &CancelGate) -> Result<FrameFluidUpdate> {
        poll(gate)?;
        self.admit_state(state)?;
        if fluids.len() != self.members.len() {
            return Err(invalid("one replacement fluid source per original region required"));
        }
        if !maximum_absolute_parameter_work_j.is_finite() || maximum_absolute_parameter_work_j < 0.0 {
            return Err(invalid("finite nonnegative material intervention work budget required"));
        }

        let mut members = Vec::with_capacity(self.members.len());
        for (original, fluid) in self.members.iter().zip(fluids) {
            poll(gate)?;
            let flow = CylinderFlow::from_fluid(original.flow.radius_m(), fluid,
                original.flow.resolution(), gate).map_err(|error| match error {
                    AcousticRealizeError::Cancelled => FrameFluidError::Cancelled,
                    other => FrameFluidError::Fluid(other),
                })?;
            members.push(LoadedFrameMember { element: original.element, span: original.span,
                length_m: original.length_m, flow });
        }
        let operators = regions::assemble(&self.basis, &members, &self.stations,
            &self.material_ports, self.mass_inverse_tolerance, gate)?;
        if operators.system.state_dim() != self.system.state_dim() {
            return Err(invalid("material intervention changed original coordinate topology"));
        }
        let before = finite(self.system.hamiltonian(state))?;
        let after = finite(operators.system.hamiltonian(state))?;
        if operators.system.effort(state).iter().any(|v| !v.is_finite()) {
            return Err(invalid("replacement fluid makes original momentum effort unrepresentable"));
        }
        let old_q = self.system.storage().exact_quadratic_matrix()
            .ok_or_else(|| invalid("original frame/fluid quadratic storage is absent"))?;
        let next_q = operators.system.storage().exact_quadratic_matrix()
            .ok_or_else(|| invalid("replacement frame/fluid quadratic storage is absent"))?;
        let parameter_work_j = momentum_work(old_q, next_q, self.basis.modes().len(), state, gate)?;
        if parameter_work_j.abs() > maximum_absolute_parameter_work_j {
            return Err(FrameFluidError::ParameterWork { parameter_work_j,
                maximum_absolute_parameter_work_j });
        }
        let update = FrameFluidUpdate {
            stored_before_j: before,
            stored_after_j: after,
            parameter_work_j,
            energy_residual_j: finite(after-before-parameter_work_j)?,
            operator_changed: old_q != next_q
                || self.system.structure().1 != operators.system.structure().1,
        };
        poll(gate)?;
        // No fallible operation follows publication. State is borrowed
        // immutably, and all basis/trace/structural owners remain untouched.
        self.members = members;
        self.modal_mass = operators.modal_mass;
        self.memory_masses = operators.memory_masses;
        self.inverse_residual = operators.inverse_residual;
        self.system = operators.system;
        Ok(update)
    }
}

/// Compare only changed kinetic metrics: subtracting two complete Hamiltonians
/// can hide small intervention work underneath a large unchanged elastic energy.
fn momentum_work(old_q: &[f64], next_q: &[f64], modes: usize, state: &[f64],
    gate: &CancelGate) -> Result<f64> {
    let n = state.len();
    let mut work = 0.0;
    for i in 0..modes {
        poll(gate)?;
        let row = modes+i;
        let mut velocity_change = 0.0;
        for k in 0..modes {
            let index = row*n+modes+k;
            velocity_change = finite(velocity_change
                + finite(next_q[index]-old_q[index])?*state[modes+k])?;
        }
        work = finite(work+0.5*state[row]*velocity_change)?;
    }
    for i in 2*modes..n {
        if i%128 == 0 { poll(gate)?; }
        work = finite(work+0.5*state[i]*finite(next_q[i*n+i]-old_q[i*n+i])?*state[i])?;
    }
    Ok(work)
}

