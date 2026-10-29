//! Reciprocal acoustic storage in the original complete modal propagator.
//! The shared passive model is also used by piano and impact consumers.
use super::*;
use crate::render::plate::impact::radiation::Model;
use fs_math::det;

#[derive(Clone, Debug)]
pub(super) struct Binding {
    pub(super) model: Model,
    pub(super) source_modes: Vec<usize>,
    pub(super) maximum_poles: usize,
    pub(super) mechanical_bound: f64,
    coupling_bound: f64,
    pole_bound: f64,
}
impl Binding {
    pub(super) fn frequency_bound(&self, mechanical: f64) -> Result<f64, ModalCouplingError> {
        finite(mechanical.max(self.pole_bound) + self.coupling_bound)
    }
    pub(super) fn energy(&self, state: &[f64]) -> Result<f64, ModalCouplingError> {
        if state.len() != 2*self.model.poles.len() || state.iter().any(|v| !v.is_finite()) {
            return Err(invalid("complete finite acoustic state required"));
        }
        self.model.poles.iter().zip(state.chunks_exact(2)).try_fold(0.0, |sum, (p, x)|
            finite(sum + 0.5*((p.omega*x[0]).powi(2) + x[1]*x[1])))
    }
}

impl CoupledModalSystem {
    /// Cold-attach the shared positive-real exterior impedance without changing
    /// mechanical force addresses, material histories or the integration owner.
    /// Every signed model column maps one distinct original velocity; all mutual
    /// reactions survive. Acoustic coordinates start at rest and are never
    /// supplied to a mechanical contact or pressure observer as extra modes.
    ///
    /// The producer owns geometry, gas, basis and fit validity. This method
    /// proves composition, not acoustic fidelity of an authored pole model.
    /// Complete acoustic/material states share the ordinary exact capacity and
    /// setup limits, up to 256 poles. A prepared owner keeps the poles as its
    /// air bank beside the transition (`exact::prepared::air`), up to the shared
    /// model's own pole limit, within the prepared byte and work budget.
    /// A completely zero coupling preserves the original system exactly.
    pub fn with_radiation_load(mut self, model: &Model, source_modes: &[usize],
        maximum_poles: usize, gate: &CancelGate) -> Result<Self, ModalCouplingError>
    {
        poll(Some(gate))?;
        model.validate().map_err(|_| invalid("invalid complete passive acoustic model"))?;
        let ceiling = if self.prepared_budget.is_some() {
            crate::render::plate::impact::radiation::MAX_POLES
        } else { 256 };
        if self.samples_rendered() != 0 || self.radiation.is_some()
            || maximum_poles > ceiling || model.poles.len() > maximum_poles
            || source_modes.len() != model.ports
            || source_modes.iter().enumerate().any(|(i, &k)|
                k >= self.mode_count() || source_modes[..i].contains(&k)) {
            return Err(invalid("radiation requires one cold complete mapping within its pole allowance"));
        }
        let mut column_sums = vec![0.0; model.ports];
        let mut row_max = 0.0_f64;
        let mut pole_bound = 0.0_f64;
        for p in &model.poles {
            poll(Some(gate))?;
            finite(p.omega*p.omega)?;
            finite(2.0*p.zeta*p.omega)?;
            pole_bound = pole_bound.max(p.omega);
            let mut sum = 0.0;
            for (column, &value) in column_sums.iter_mut().zip(&p.coupling) {
                *column = finite(*column + value.abs())?;
                sum = finite(sum + value.abs())?;
            }
            row_max = row_max.max(sum);
        }
        let column_max = column_sums.into_iter().fold(0.0_f64, f64::max);
        let coupling_bound = finite(det::sqrt(row_max)*det::sqrt(column_max))?;
        if coupling_bound == 0.0 { return Ok(self); }
        if row_max.max(column_max)*self.dt > 0.25 {
            return Err(invalid("acoustic coupling exceeds the declared mechanical clock resolution"));
        }
        let binding = Binding { model: model.clone(), source_modes: source_modes.to_vec(),
            maximum_poles, mechanical_bound: self.linear_frequency_upper_bound_rad_s,
            coupling_bound, pole_bound };
        let bound = binding.frequency_bound(binding.mechanical_bound)?;
        if bound > core::f64::consts::PI/self.dt*self.config.nyquist_guard_fraction {
            return Err(invalid("complete acoustic/mechanical frequency bound exceeds Nyquist admission"));
        }
        let exact = if let Some(budget) = self.prepared_budget {
            exact::ExactNetwork::new_prepared_with_relaxation(&self.models, &self.columns, &self.connections,
                self.mode_count(), self.dt, budget, &self.relaxation_branches, Some(&binding), gate)?
        } else {
            self.admit_augmented_setup(self.relaxation_branches.len(), 2*model.poles.len())?;
            exact::ExactNetwork::new_with_extensions(&self.models, &self.columns,
                &self.connections, self.mode_count(), &self.relaxation_branches, Some(&binding), gate)?
        };
        poll(Some(gate))?;
        self.radiation_state = vec![0.0; 2*model.poles.len()];
        self.candidate_radiation_state = self.radiation_state.clone();
        self.radiation = Some(binding);
        self.exact = Some(exact);
        self.linear_frequency_upper_bound_rad_s = bound;
        Ok(self)
    }

    pub(super) fn admit_augmented_setup(&self, material: usize, acoustic: usize)
        -> Result<(), ModalCouplingError>
    {
        let base = self.mode_count().checked_mul(2).ok_or_else(|| invalid("augmented base extent"))?;
        let state = base.checked_add(material).and_then(|n| n.checked_add(acoustic))
            .ok_or_else(|| invalid("complete augmented state extent"))?;
        let links = self.connections.len()+1;
        let terms = self.mode_count().checked_mul(links).and_then(|n| n.checked_mul(links))
            .and_then(|n| n.checked_add(state.checked_mul(state)?))
            .and_then(|n| n.checked_add(base.checked_mul(material.checked_add(acoustic)?)?))
            .ok_or_else(|| invalid("complete acoustic/material setup overflow"))?;
        if terms > self.config.max_setup_terms {
            return Err(invalid("complete acoustic/material setup exceeds max_setup_terms"));
        }
        Ok(())
    }

    /// Accepted acoustic q/p pairs, in original model pole order.
    #[must_use]
    pub fn radiation_state(&self) -> &[f64] { &self.radiation_state }
    /// The original complete passive load, when one is attached.
    #[must_use]
    pub fn radiation_model(&self) -> Option<&Model> { self.radiation.as_ref().map(|r| &r.model) }
    /// Acoustic storage, already included in the complete network energy.
    pub fn radiation_energy_j(&self) -> Result<f64, ModalCouplingError> {
        self.radiation.as_ref().map_or(Ok(0.0), |r| r.energy(&self.radiation_state))
    }

    pub(super) fn stage_radiation_state(&mut self, endpoint: &[f64]) -> Result<f64, ModalCouplingError> {
        let start = 2*self.mode_count() + self.relaxation_branches.len();
        if endpoint.len() != start + self.radiation_state.len() {
            return Err(invalid("complete acoustic endpoint dimension"));
        }
        let Some(binding) = &self.radiation else { return Ok(0.0); };
        self.candidate_radiation_state.copy_from_slice(&endpoint[start..]);
        let energy = binding.energy(&self.candidate_radiation_state)?;
        for (pole, state) in binding.model.poles.iter().zip(self.candidate_radiation_state.chunks_exact(2)) {
            for (&mode, &coupling) in binding.source_modes.iter().zip(&pole.coupling) {
                // Original parallel interconnection: F_host = -L^T p_air.
                self.forces[mode] = finite(self.forces[mode] - coupling*state[1])?;
            }
        }
        Ok(energy)
    }
}

#[cfg(test)]
mod tests;
