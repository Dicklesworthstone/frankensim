//! Original physical-state conversion and native frame/fluid step admission.
use super::*;

impl<F> MaterialFluidFrame<'_, F> {
    /// Initialize original modal displacements/velocities plus ALL physical
    /// quadrature fluid velocities. No fluid-rest assumption or history reset.
    pub fn state_from_motion(&self, modes: &[ModalAcousticState], fluid_velocities: &[f64]) -> Result<Vec<f64>> {
        let m = self.basis.modes().len();
        if modes.len() != m || fluid_velocities.len() != self.memory_masses.len()
            || modes.iter().any(|s| !s.displacement_m_sqrt_kg.is_finite() || !s.velocity_m_sqrt_kg_per_s.is_finite())
            || fluid_velocities.iter().any(|v| !v.is_finite()) {
            return Err(invalid("complete finite modal motion and all fluid velocities required"));
        }
        let velocity: Vec<_> = modes.iter().map(|s| s.velocity_m_sqrt_kg_per_s).collect();
        let mut state = vec![0.0; self.system.state_dim()];
        for i in 0..m {
            state[i] = modes[i].displacement_m_sqrt_kg;
            state[m+i] = dot(&self.modal_mass[i*m..(i+1)*m], &velocity)?;
        }
        for (k, (&mass, &v)) in self.memory_masses.iter().zip(fluid_velocities).enumerate() {
            state[2*m+k] = finite(mass*v)?;
            if v != 0.0 && state[2*m+k] == 0.0 { return Err(invalid("fluid momentum underflow")); }
        }
        self.admit_state(&state)?;
        Ok(state)
    }

    /// Recover original-basis q and qdot=M_inverse*p, not momentum-as-velocity.
    /// Feed these to basis().reconstruct_into for physical nodal observations.
    pub fn modal_states(&self, state: &[f64]) -> Result<Vec<ModalAcousticState>> {
        self.admit_state(state)?;
        let m = self.basis.modes().len();
        let effort = self.system.effort(state);
        Ok((0..m).map(|i| ModalAcousticState { displacement_m_sqrt_kg: state[i],
            velocity_m_sqrt_kg_per_s: effort[m+i] }).collect())
    }

    /// Advance held incremental nodal forces/moments using the ORIGINAL pHS
    /// stepper. Fixed supports and stationary far-field fluid are explicit.
    /// Do not apply the frozen equilibrium's preload twice. Cut force events at
    /// interval boundaries. Refine time, dry modes and fluid spectrum separately.
    ///
    /// The caller's state/source frame is never mutated. Publication follows
    /// original material clock admission, native solve, finite-state admission,
    /// and independent nodal-work / separate solid-and-fluid loss checks.
    pub fn step_candidate(&self, state: &[f64], nodal_force: &[f64], dt: f64,
        energy_tolerance_j: f64, gate: &CancelGate) -> Result<FrameFluidStep> {
        poll(gate)?;
        positive(dt)?;
        positive(energy_tolerance_j)?;
        self.admit_state(state)?;
        if let Some(damping) = self.basis.material_damping() { damping.admit_timestep(dt)?; }
        let m = self.basis.modes().len();
        let mut input = vec![0.0; m];
        self.basis.project_force_into(nodal_force, [0.0; 2], &mut input, || gate.is_requested())?;
        poll(gate)?;
        let native = fs_phs::step(&self.system, state, &input, dt)?;
        poll(gate)?;
        self.admit_state(&native.x)?;
        let midpoint: Vec<_> = state.iter().zip(&native.x).map(|(a,b)| 0.5*a+0.5*b).collect();
        let effort = self.system.effort(&midpoint);
        let velocity = &effort[m..2*m];
        let mut material_port_dissipation_j = Vec::with_capacity(self.material_ports.len());
        for port in &self.material_ports {
            poll(gate)?;
            let speed = dot(&port.projection, velocity)?;
            material_port_dissipation_j.push(finite(dt*port.resistance*speed*speed)?);
        }
        let mut fluid_member_dissipation_j = vec![0.0; self.members.len()];
        for station in &self.stations {
            poll(gate)?;
            let speed = dot(&station.projection, velocity)?;
            let flow = &self.members[station.member].flow;
            let mut power = station.weight_m*flow.low_tail_resistance_per_length_n_s_m2()*speed*speed;
            for (k, branch) in flow.branches().iter().enumerate() {
                let slip = speed-effort[2*m+station.first_memory+k];
                power += station.weight_m*branch.resistance_per_length_n_s_m2*slip*slip;
            }
            fluid_member_dissipation_j[station.member] = finite(
                fluid_member_dissipation_j[station.member]+dt*power)?;
        }
        let mut nodal_work_j = 0.0;
        for (dof, force) in nodal_force.iter().enumerate() {
            if dof % 128 == 0 { poll(gate)?; }
            let increment: f64 = self.basis.modes().iter().enumerate().map(|(i,mode)|
                mode.displacement[dof]*(native.x[i]-state[i])).sum();
            nodal_work_j = finite(nodal_work_j+force*increment)?;
        }
        let stored_before_j = self.system.hamiltonian(state);
        let stored_after_j = self.system.hamiltonian(&native.x);
        let loss: f64 = material_port_dissipation_j.iter().chain(&fluid_member_dissipation_j).sum();
        let energy_residual_j = finite(stored_after_j-stored_before_j+loss-nodal_work_j)?;
        let port_residual = finite(native.supplied-nodal_work_j)?;
        if energy_residual_j.abs() > energy_tolerance_j || port_residual.abs() > energy_tolerance_j {
            return Err(FrameFluidError::Balance {
                residual_j: if energy_residual_j.abs() >= port_residual.abs() { energy_residual_j } else { port_residual },
                tolerance_j: energy_tolerance_j,
            });
        }
        poll(gate)?;
        Ok(FrameFluidStep { native, stored_before_j, stored_after_j, material_port_dissipation_j,
            fluid_member_dissipation_j, nodal_work_j, energy_residual_j })
    }

    pub(super) fn admit_state(&self, state: &[f64]) -> Result<()> {
        if state.len() != self.system.state_dim() || state.iter().any(|v| !v.is_finite()) {
            return Err(invalid("complete finite frame/fluid state required"));
        }
        finite(self.system.hamiltonian(state))?;
        if self.system.effort(state).iter().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable frame/fluid effort"));
        }
        Ok(())
    }
}
