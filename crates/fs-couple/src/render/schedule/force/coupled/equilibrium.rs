//! Full-network static preload; never independent equilibria after coupling.
use super::*;

mod preload;
mod normal;
mod free;
use free::FreeResponse;

impl CoupledModalSystem {
    /// Establish a settled equilibrium of the COMPLETE spring-connected network.
    /// All components must still be unexcited at sample zero. Dashpots carry no
    /// static force. The supplied load is assumed to have settled before this
    /// window; returned storage is not a measured work history. This offline
    /// initialization allocates, uses the same setup-size cap as construction,
    /// and publishes no states on a cancellation, residual or budget refusal.
    /// Explicit free coordinates are legal only when the declared bilateral
    /// springs constrain all of them. No artificial tether or pose is supplied.
    /// The additional bounded support solve uses at most max_connections rows.
    pub fn initialize_static_equilibrium(&mut self, external: &[f64], gate: &CancelGate)
        -> Result<f64, ModalCouplingError>
    {
        poll(Some(gate))?;
        self.require_static_initialization(external)?;
        let response = StaticResponse::new(self, gate)?;
        let q = response.solve(self, external, true, true, gate)?;
        let energy = self.stage_static_equilibrium(&q, external, gate)?;
        poll(Some(gate))?;
        std::mem::swap(&mut self.models, &mut self.candidates);
        Ok(energy)
    }

    /// Settle held forces against bilateral springs AND all declared contacts.
    /// Uses the same mass-normalized coordinates and existing fs-dcontact law.
    /// Every component must be at zero Q/V before sample zero. Contact activity
    /// is solved, not supplied; separated contacts have zero reaction. Contact
    /// and bilateral dashpots exert no static force. Returns total stored energy
    /// including contacts, not a reconstructed pre-window actuator-work history.
    ///
    /// The nonlinear sweep/root/setup caps are explicit. Final contact forces,
    /// penetrations, actual modal force balance and total energy must all pass
    /// before any component changes. A failure/cancellation leaves the complete
    /// accepted network unchanged. This method does not attach the contacts:
    /// pass these same descriptions to ContactModalSystem/MultiContactModalSystem
    /// to continue their dynamics, and retain the held forces until release.
    /// Free coordinates require an independently supported bilateral network;
    /// a mass supported ONLY by unilateral contacts remains outside this solve.
    pub fn initialize_contact_equilibrium(
        &mut self,
        external: &[f64],
        contacts: &[(contact::ModalContact, contact::ModalContactConfig)],
        config: contact::multiple::MultiContactConfig,
        gate: &CancelGate,
    ) -> Result<f64, ModalCouplingError> {
        preload::initialize(self, external, contacts, config, gate)
    }

    fn require_static_initialization(&self, external: &[f64]) -> Result<(), ModalCouplingError> {
        if self.samples_rendered() != 0 || external.len() != self.mode_count()
            || external.iter().any(|x| !x.is_finite())
            || self.models.iter().flat_map(|m| m.states()).any(|s|
                s.displacement_m_sqrt_kg != 0.0 || s.velocity_m_sqrt_kg_per_s != 0.0) {
            return Err(invalid("coupled preload requires zero initial vibration, sample zero and finite complete forces"));
        }
        Ok(())
    }

    // Stage only final displacements, never the unrestrained free prediction.
    // The same component and bilateral-force gates serve both preload paths.
    fn stage_static_equilibrium(&mut self, q: &[f64], external: &[f64], gate: &CancelGate)
        -> Result<f64, ModalCouplingError>
    {
        for i in 0..self.models.len() {
            poll(Some(gate))?;
            let states: Vec<_> = (self.offsets[i]..self.offsets[i+1]).map(|j|
                ModalAcousticState { displacement_m_sqrt_kg: q[j], velocity_m_sqrt_kg_per_s: 0.0 }).collect();
            self.candidates[i].restore_states(&states)
                .map_err(|source| ModalCouplingError::Component { component: i, source })?;
        }
        // Actual rounded states and spring extensions, not solved variables.
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
        Ok(energy)
    }
}

// One static inverse-stiffness owner, also used for contact cross-compliance.
// The dynamic exact-ZOH response is NOT the static inverse stiffness.
struct StaticResponse {
    compliance: Vec<f64>,
    roots: Vec<f64>,
    matrix: Vec<f64>,
    factor: Cholesky,
    free: Option<FreeResponse>,
}
impl StaticResponse {
    fn new(network: &CoupledModalSystem, gate: &CancelGate) -> Result<Self, ModalCouplingError> {
        let free_count = network.models.iter().flat_map(|m| m.modes())
            .filter(|m| m.angular_frequency_rad_s == 0.0).count();
        FreeResponse::admit(network, free_count)?;
        let mut indices = Vec::with_capacity(free_count);
        let mut compliance = Vec::with_capacity(network.mode_count());
        for model in &network.models {
            poll(Some(gate))?;
            for mode in model.modes() {
                if mode.angular_frequency_rad_s == 0.0 {
                    indices.push(compliance.len());
                    compliance.push(0.0); // placeholder; free positions come from the support solve
                } else {
                    let d = finite((mode.angular_frequency_rad_s * mode.angular_frequency_rad_s).recip())?;
                    if d <= 0.0 { return Err(invalid("positive static modal compliance is not representable")); }
                    compliance.push(d);
                }
            }
        }
        let roots: Vec<f64> = network.connections.iter().map(|c| c.stiffness_n_m.sqrt()).collect();
        let matrix = connection_matrix(&network.columns, &compliance, &roots, Some(gate))?;
        let factor = cholesky(&matrix, roots.len()).map_err(ModalCouplingError::Factor)?;
        let free = if indices.is_empty() { None } else {
            Some(FreeResponse::new(network, indices, &roots, &matrix, &factor, gate)?)
        };
        Ok(Self { compliance, roots, matrix, factor, free })
    }

    fn solve(&self, network: &CoupledModalSystem, external: &[f64], include_rest: bool,
        enforce_limits: bool, gate: &CancelGate) -> Result<Vec<f64>, ModalCouplingError>
    {
        let free: Vec<f64> = self.compliance.iter().zip(external).map(|(d,g)| finite(d*g)).collect::<Result<_,_>>()?;
        let mut rhs = Vec::with_capacity(self.roots.len());
        for ((column, link), root) in network.columns.iter().zip(&network.connections).zip(&self.roots) {
            let rest = if include_rest { link.rest_extension_m } else { 0.0 };
            rhs.push(finite(root * (dot(column, &free)? - rest))?);
        }
        let mut solution = rhs.clone();
        self.factor.solve(&mut solution);
        check_solve(&self.matrix, &solution, &rhs, network.config.solve_relative_tolerance)?;
        let free_q = if let Some(free) = &self.free {
            let q = free.complete(external, &mut rhs, &solution, network.config.solve_relative_tolerance, gate)?;
            solution.copy_from_slice(&rhs);
            self.factor.solve(&mut solution);
            check_solve(&self.matrix, &solution, &rhs, network.config.solve_relative_tolerance)?;
            Some(q)
        } else { None };
        let mut total = external.to_vec();
        for (j, column) in network.columns.iter().enumerate() {
            poll(Some(gate))?;
            let reaction = finite(-self.roots[j] * solution[j])?;
            if enforce_limits {
                limit("static connection force", reaction.abs(), network.config.maximum_abs_connection_force_n)?;
            }
            for (g,b) in total.iter_mut().zip(column) { *g = finite(*g + b*reaction)?; }
        }
        let mut q: Vec<f64> = self.compliance.iter().zip(total).map(|(d,g)| finite(d*g)).collect::<Result<_,_>>()?;
        if let (Some(free), Some(free_q)) = (&self.free, &free_q) {
            free.check_and_insert(external, &solution, free_q, &mut q, network.config.solve_relative_tolerance, gate)?;
        }
        Ok(q)
    }
}
