//! Physical characteristic-boundary work and numerical upwind loss.
//!
//! At an outward-oriented boundary, let v be the interior volume-flow trace,
//! z its local characteristic impedance, a=R/(R+z), and b=z/(R+z).
//! The original flux is p_hat=a*(p+z*v), v_hat=b*(p/z+v). Its actual load
//! power is a*b*(p/sqrt(z)+sqrt(z)*v)^2. The remaining DG loss is the
//! upwind mismatch (b*p/sqrt(z)-a*sqrt(z)*v)^2. Their sum is the original
//! boundary operator's b*p^2/z+a*z*v^2, but only the first is load work.
//! These two source squares are integrated on the same complete held-input
//! trajectory; no stored-energy deficit is relabeled as physical absorption.

use super::*;
use fs_phs::PreparedStepError;
use fs_phs::linear_exact::{DenseExponentialBudget, PreparedSourceLoss};

/// Independent boundary energies over one original held-input interval.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DgBoundaryEnergy {
    /// Work delivered to physical matched/resistive loads [J]. A matched
    /// exterior carries this energy away; a resistive termination absorbs it.
    pub load_absorbed_j: f64,
    /// Nonnegative loss from the DG characteristic boundary mismatch [J].
    /// This is a discretization error, never a physical heat source.
    pub numerical_loss_j: f64,
}

impl DgBoundaryEnergy {
    /// Total loss represented by the original boundary operator [J].
    pub fn total_j(self) -> f64 {
        self.load_absorbed_j + self.numerical_loss_j
    }
}

/// Original step with independently integrated physical and numerical loss.
#[derive(Clone, Debug, PartialEq)]
pub struct DgBoundaryStepOutput {
    pub step: DgStepOutput,
    pub boundary: DgBoundaryEnergy,
}

/// Source-bound boundary-energy observations of the complete native trajectory.
/// Warm evaluation performs O(n^2+n*ports) work with no exponential or solve.
#[derive(Clone, Debug)]
pub struct PreparedDgBoundaryEnergy {
    operator_identity: ContentHash,
    ports: Vec<PortId>,
    dimension: usize,
    physical: PreparedSourceLoss,
    numerical: PreparedSourceLoss,
}

fn invalid(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled { completed_steps: 0 })
    } else {
        Ok(())
    }
}

fn integral_error(error: PreparedStepError) -> DuctError {
    match error {
        PreparedStepError::Cancelled => DuctError::Cancelled { completed_steps: 0 },
        PreparedStepError::Solver(_) => invalid("DG boundary source-energy integration unresolved"),
    }
}

fn source_coefficients(
    boundary: DgBoundary,
    zc: f64,
    normal: f64,
) -> Result<Option<([f64; 2], [f64; 2])>, DuctError> {
    let Some((alpha, beta)) = boundary.characteristic_weights(zc) else {
        return Ok(None);
    };
    let root_zc = det::sqrt(zc);
    let root_ab = det::sqrt(alpha) * det::sqrt(beta);
    let load = [root_ab / root_zc, normal * root_ab * root_zc];
    let mismatch = [beta / root_zc, -normal * alpha * root_zc];
    if !zc.is_finite() || zc <= 0.0 || load.iter().chain(&mismatch).any(|v| !v.is_finite()) {
        return Err(invalid(
            "DG characteristic boundary source is unrepresentable",
        ));
    }
    Ok(Some((load, mismatch)))
}

impl PreparedDgBoundaryEnergy {
    pub const fn operator_identity(&self) -> ContentHash {
        self.operator_identity
    }

    /// Physical held-input interval [s], identical to the original propagator.
    pub fn duration(&self) -> f64 {
        self.physical.duration()
    }

    /// Exact caller-owned scratch length for allocation-free warm observation.
    pub fn scratch_len(&self) -> usize {
        self.ports.len() + 2 * self.dimension
    }

    // Preserve the native ordered sum of duplicate inputs, polling before
    // each entry. The resulting vector is bounded by the declared port count.
    fn resolve_inputs_into(
        &self,
        inputs: &[PortInput],
        efforts: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DuctError> {
        efforts.fill(0.0);
        for input in inputs {
            poll(cancelled)?;
            if !input.effort.is_finite() {
                return Err(invalid("DG boundary energy port effort must be finite"));
            }
            let slot = self
                .ports
                .iter()
                .position(|id| *id == input.id)
                .ok_or_else(|| invalid("DG boundary energy references an undeclared port"))?;
            efforts[slot] += input.effort;
            if !efforts[slot].is_finite() {
                return Err(invalid("DG boundary energy summed port effort overflows"));
            }
        }
        Ok(())
    }

    /// Integrate boundary source squares from the INITIAL complete state and
    /// the same held inputs as the native propagator. No state advances here.
    pub fn integrated_energy_j(
        &self,
        state: &DgState,
        inputs: &[PortInput],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgBoundaryEnergy, DuctError> {
        let mut scratch = vec![0.0; self.scratch_len()];
        self.integrated_energy_into_j(state, inputs, &mut scratch, cancelled)
    }

    /// Same positive integrals without warm allocation. Scratch may change on
    /// refusal; the original physical state and its clock never do.
    pub fn integrated_energy_into_j(
        &self,
        state: &DgState,
        inputs: &[PortInput],
        scratch: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgBoundaryEnergy, DuctError> {
        poll(cancelled)?;
        if state.operator_identity != self.operator_identity
            || state.tilde.len() != self.dimension
            || scratch.len() != self.scratch_len()
        {
            return Err(invalid(
                "DG boundary energy requires its original state and exact scratch",
            ));
        }
        let (efforts, coordinates) = scratch.split_at_mut(self.ports.len());
        self.resolve_inputs_into(inputs, efforts, cancelled)?;
        let load_absorbed_j = self
            .physical
            .evaluate_into_controlled(&state.tilde, efforts, coordinates, &mut *cancelled)
            .map_err(integral_error)?;
        let numerical_loss_j = self
            .numerical
            .evaluate_into_controlled(&state.tilde, efforts, coordinates, &mut *cancelled)
            .map_err(integral_error)?;
        let energy = DgBoundaryEnergy {
            load_absorbed_j,
            numerical_loss_j,
        };
        if !energy.total_j().is_finite() {
            return Err(invalid("DG total boundary energy overflows"));
        }
        poll(cancelled)?;
        Ok(energy)
    }
}

impl<M: DgMedium> DgDuct<M> {
    /// Prepare actual load-work and numerical-loss integrals on this duct's
    /// unchanged complete operator and the original propagator's interval.
    /// Local material and area set each boundary's impedance; the inlet's
    /// outward flow is negative axial flow. Both load and mismatch remain
    /// positive squares, including every held-input/history cross term.
    ///
    /// All attached memories participate in the trajectory, but their own
    /// dissipation belongs to their original heat owners. This cache reports
    /// boundary loss only. The native source-integral owner bounds the complete
    /// state to 128 coordinates; no spatial coordinate or memory is discarded.
    /// Its two caches share the supplied byte/work allowance. Closed, open and
    /// prescribed-flow boundaries contribute no characteristic loss rows.
    pub fn prepare_boundary_energy(
        &self,
        propagator: &DgPropagator,
        budget: DenseExponentialBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedDgBoundaryEnergy, DuctError> {
        poll(cancelled)?;
        let n = self.total_dofs();
        let ports = self.ports.len();
        if propagator.operator_identity != self.operator_identity
            || n > 128
            || 2 * n > budget.max_dimension
        {
            return Err(invalid(
                "DG boundary energy needs its native propagator and complete bounded state",
            ));
        }
        let bytes = ports
            .checked_add(4)
            .and_then(|rows| rows.checked_mul(n))
            .and_then(|entries| entries.checked_mul(8))
            .ok_or_else(|| invalid("DG boundary energy source extent overflows"))?;
        if bytes > budget.max_dense_bytes {
            return Err(invalid("DG boundary sources exceed the byte allowance"));
        }
        let channel_budget = DenseExponentialBudget {
            max_dimension: budget.max_dimension,
            max_dense_bytes: (budget.max_dense_bytes - bytes) / 2,
            max_multiply_terms: budget.max_multiply_terms / 2,
        };
        let local = self.mesh.order + 1;
        let last = self.elements.len() - 1;
        let mut physical = Vec::with_capacity(2 * n);
        let mut numerical = Vec::with_capacity(2 * n);
        for (boundary, element, node, medium, normal) in [
            (self.inlet, 0, 0, &self.gas, -1.0),
            (self.outlet, last, local - 1, &self.outlet_gas, 1.0),
        ] {
            poll(cancelled)?;
            let elem = &self.elements[element];
            let zc = medium.density_kg_m3() * medium.sound_speed_m_s() / elem.areas[node];
            let Some((load, mismatch)) = source_coefficients(boundary, zc, normal)? else {
                continue;
            };
            let start = physical.len();
            physical.resize(start + n, 0.0);
            numerical.resize(start + n, 0.0);
            let base = element * local;
            for j in 0..local {
                let p = elem.mp_inv_sqrt[node * local + j];
                let q = elem.mq_inv_sqrt[node * local + j];
                physical[start + base + j] = load[0] * p;
                physical[start + self.field_dofs + base + j] = load[1] * q;
                numerical[start + base + j] = mismatch[0] * p;
                numerical[start + self.field_dofs + base + j] = mismatch[1] * q;
            }
        }
        let mut input_map = vec![0.0; n * ports];
        for (port, (_, column)) in self.ports.iter().enumerate() {
            poll(cancelled)?;
            for i in 0..n {
                input_map[i * ports + port] = column[i];
            }
        }
        let physical = PreparedSourceLoss::new_controlled(
            n,
            ports,
            &self.l_tilde,
            &input_map,
            &physical,
            propagator.dt,
            channel_budget,
            &mut *cancelled,
        )
        .map_err(integral_error)?;
        let numerical = PreparedSourceLoss::new_controlled(
            n,
            ports,
            &self.l_tilde,
            &input_map,
            &numerical,
            propagator.dt,
            channel_budget,
            &mut *cancelled,
        )
        .map_err(integral_error)?;
        Ok(PreparedDgBoundaryEnergy {
            operator_identity: self.operator_identity,
            ports: self.ports.iter().map(|(id, _)| *id).collect(),
            dimension: n,
            physical,
            numerical,
        })
    }

    /// Advance the ORIGINAL propagator transactionally and accept the step
    /// against independently integrated load work plus numerical boundary loss.
    /// This closes the current interval without using its energy deficit as a
    /// loss estimate. Only load_absorbed_j may feed a receiving physical owner.
    ///
    /// This complete balance is restricted to boundary-only dissipation: ducts
    /// carrying porous, viscothermal, molecular or compliant-wall attachments
    /// refuse until their independent losses join the same accounting window.
    /// Their boundary-energy queries remain supported by the preparation above.
    /// Original state, accepted clock and cumulative work survive any refusal
    /// or cancellation. Previous history retains its original accounting scope.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_boundary_energy(
        &self,
        propagator: &DgPropagator,
        boundary: &PreparedDgBoundaryEnergy,
        state: &mut DgState,
        inputs: &[PortInput],
        absolute_tolerance_j: f64,
        relative_tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgBoundaryStepOutput, DuctError> {
        poll(cancelled)?;
        if boundary.operator_identity != self.operator_identity
            || propagator.operator_identity != self.operator_identity
            || boundary.duration().to_bits() != propagator.dt.to_bits()
            || self.viscothermal.is_some()
            || self.molecular.is_some()
            || self.compliant_wall.is_some()
            || self.gas.porous_diffusive_model().is_some()
            || !absolute_tolerance_j.is_finite()
            || absolute_tolerance_j < 0.0
            || !relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&relative_tolerance)
            || state.step_count.checked_add(1).is_none()
        {
            return Err(invalid(
                "checked DG boundary step needs matching owners, interval, boundary-only loss and finite tolerances",
            ));
        }
        let mut efforts = vec![0.0; boundary.ports.len()];
        boundary.resolve_inputs_into(inputs, &mut efforts, cancelled)?;
        let canonical_inputs: Vec<_> = boundary
            .ports
            .iter()
            .zip(&efforts)
            .map(|(&id, &effort)| PortInput { id, effort })
            .collect();
        let energy = boundary.integrated_energy_j(state, &canonical_inputs, cancelled)?;
        let mut candidate = state.clone();
        let original = propagator.step(&mut candidate, &canonical_inputs)?;
        poll(cancelled)?;
        let exchanges: Vec<_> = original
            .port_outputs
            .iter()
            .zip(&efforts)
            .map(|(port, effort)| continuous::state::DgTrajectoryPort {
                id: port.id,
                pressure_impulse_pa_s: port.delta_x,
                work_into_duct_j: effort * port.delta_x,
            })
            .collect();
        let (mut accepted, step) = self.continuous_operator().candidate_from_trajectory(
            state,
            candidate.tilde,
            propagator.dt,
            &exchanges,
            energy.total_j(),
            absolute_tolerance_j,
            relative_tolerance,
        )?;
        // Adoption of an externally composed trajectory clears these caches.
        // Here the ORIGINAL time owner supplied the endpoint and its anchors:
        // keep them so checking lossless steps does not restart modal phase or
        // introduce extra projections/roundoff at every accounting interval.
        accepted.modal_anchor_p = candidate.modal_anchor_p;
        accepted.modal_anchor_q = candidate.modal_anchor_q;
        accepted.modal_anchor_dt_bits = candidate.modal_anchor_dt_bits;
        accepted.unforced_steps = candidate.unforced_steps;
        poll(cancelled)?;
        *state = accepted;
        Ok(DgBoundaryStepOutput {
            step,
            boundary: energy,
        })
    }
}

#[cfg(test)]
mod tests;
