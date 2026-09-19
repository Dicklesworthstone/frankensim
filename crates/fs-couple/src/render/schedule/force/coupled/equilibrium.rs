//! Full-network static preload; never independent equilibria after coupling.
use super::*;

impl CoupledModalSystem {
    /// Establish a settled equilibrium of the COMPLETE spring-connected network.
    /// All components must still be unexcited at sample zero. Dashpots carry no
    /// static force. The supplied load is assumed to have settled before this
    /// window; returned storage is not a measured work history. This offline
    /// initialization allocates, uses the same setup-size cap as construction,
    /// and publishes no states on a cancellation, residual or budget refusal.
    pub fn initialize_static_equilibrium(&mut self, external: &[f64], gate: &CancelGate)
        -> Result<f64, ModalCouplingError>
    {
        poll(Some(gate))?;
        if self.samples_rendered() != 0 || external.len() != self.mode_count()
            || external.iter().any(|x| !x.is_finite())
            || self.models.iter().flat_map(|m| m.states()).any(|s|
                s.displacement_m_sqrt_kg != 0.0 || s.velocity_m_sqrt_kg_per_s != 0.0) {
            return Err(invalid("coupled preload requires zero initial vibration, sample zero and finite complete forces"));
        }
        let mut compliance = Vec::with_capacity(self.mode_count());
        for model in &self.models {
            for mode in model.modes() {
                let d = finite((mode.angular_frequency_rad_s * mode.angular_frequency_rad_s).recip())?;
                if d <= 0.0 { return Err(invalid("positive static modal compliance is not representable")); }
                compliance.push(d);
            }
        }
        let roots: Vec<f64> = self.connections.iter().map(|c| c.stiffness_n_m.sqrt()).collect();
        let matrix = connection_matrix(&self.columns, &compliance, &roots, Some(gate))?;
        let factor = cholesky(&matrix, roots.len()).map_err(ModalCouplingError::Factor)?;
        let free: Vec<f64> = compliance.iter().zip(external).map(|(d,g)| finite(d*g)).collect::<Result<_,_>>()?;
        let mut rhs = Vec::with_capacity(roots.len());
        for ((column, link), root) in self.columns.iter().zip(&self.connections).zip(&roots) {
            rhs.push(finite(root * (dot(column, &free)? - link.rest_extension_m))?);
        }
        let mut solution = rhs.clone();
        factor.solve(&mut solution);
        check_solve(&matrix, &solution, &rhs, self.config.solve_relative_tolerance)?;
        let mut total = external.to_vec();
        for (j, column) in self.columns.iter().enumerate() {
            poll(Some(gate))?;
            let reaction = finite(-roots[j] * solution[j])?;
            limit("static connection force", reaction.abs(), self.config.maximum_abs_connection_force_n)?;
            for (g,b) in total.iter_mut().zip(column) { *g = finite(*g + b*reaction)?; }
        }
        for i in 0..self.models.len() {
            poll(Some(gate))?;
            let states: Vec<_> = (self.offsets[i]..self.offsets[i+1]).map(|j|
                Ok(ModalAcousticState { displacement_m_sqrt_kg: finite(compliance[j]*total[j])?,
                    velocity_m_sqrt_kg_per_s: 0.0 })).collect::<Result<_,ModalCouplingError>>()?;
            self.candidates[i].restore_states(&states)
                .map_err(|source| ModalCouplingError::Component { component: i, source })?;
        }
        // Recompute the stationary force balance from the ACTUAL rounded states
        // and spring extensions, not just the solved reaction variables.
        let mut implied = external.to_vec();
        let mut scale: Vec<f64> = external.iter().map(|g| g.abs()).collect();
        let mut energy = component_energy(&self.candidates)?;
        for (column, link) in self.columns.iter().zip(&self.connections) {
            poll(Some(gate))?;
            let x = extension(&self.candidates, column, link.rest_extension_m)?;
            let reaction = finite(-link.stiffness_n_m*x)?;
            limit("static connection force", reaction.abs(), self.config.maximum_abs_connection_force_n)?;
            energy = finite(energy + 0.5*link.stiffness_n_m*x*x)?;
            for j in 0..implied.len() {
                implied[j] = finite(implied[j]+column[j]*reaction)?;
                scale[j] = finite(scale[j]+(column[j]*reaction).abs())?;
            }
        }
        let mut relative = 0.0_f64;
        let mut j = 0;
        for model in &self.candidates {
            for (mode,state) in model.modes().iter().zip(model.states()) {
                let actual = finite(mode.angular_frequency_rad_s*mode.angular_frequency_rad_s*state.displacement_m_sqrt_kg)?;
                let denominator = finite(scale[j]+actual.abs())?;
                let residual = finite(actual-implied[j])?.abs();
                relative = relative.max(if denominator == 0.0 { 0.0 } else { residual/denominator });
                j += 1;
            }
        }
        if relative > self.config.solve_relative_tolerance {
            return Err(ModalCouplingError::SolveResidual { relative, tolerance: self.config.solve_relative_tolerance });
        }
        limit("static total energy", energy, self.config.maximum_total_energy_j)?;
        poll(Some(gate))?;
        std::mem::swap(&mut self.models, &mut self.candidates);
        Ok(energy)
    }
}
