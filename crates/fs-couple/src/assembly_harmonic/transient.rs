//! Exact continuous loading of the accepted assembly, including load memory.
//!
//! The held quantity is the applied incremental force. Attachment flow evolves
//! continuously with the structure and the original rational states in ONE
//! generator. No held-flow split, pole dropping, retuning or direct-term inverse
//! is used. The existing fs-phs trajectory/Gramian owner performs propagation.

use fs_exec::Cx;
use fs_phs::linear_exact::DenseExponentialBudget;
use fs_plate::shell::preload::assembly::vibration::harmonic::ports::AssemblyHarmonicPorts;
use fs_plate::shell::preload::assembly::vibration::response::{
    AssemblyResponseBudget, AssemblyResponseSource,
};
use fs_vfit::mimo::MimoPassivityReport;
use fs_vfit::mimo::load::{MimoLoadError, PreparedMimoLoad};

mod joint;
use joint::{PreparedJointFlow, check, dot, invalid, multiply};

/// Independent source, connection and original exact-owner resource ceilings.
#[derive(Clone, Copy, Debug)]
pub struct AssemblyMimoTransientBudget {
    /// Original complete assembly projection and canonical-source budget.
    pub structure: AssemblyResponseBudget,
    /// Original exponential/Gramian cache budget, including all observations.
    pub exact: DenseExponentialBudget,
    /// Every physical attachment is retained, including dark/fixed ports.
    pub maximum_ports: usize,
    /// Complete [q,p,original MIMO states] dimension, without truncation.
    pub maximum_real_states: usize,
    /// Connection matrices, candidate/source clone envelope and physical output
    /// payload. The borrowed body/ports and separately bounded source/cache are
    /// excluded. Includes two original source memory allowances conservatively.
    pub maximum_connection_bytes: usize,
    /// Conservative dense connection construction and direct physical-port
    /// contraction cap. Assembly transforms retain their original source owners.
    pub maximum_connection_terms: u64,
    /// Relative gate for independent work, flow and energy checks. No absolute
    /// power/energy floor is inserted. Must be finite and in (0, 1e-8].
    pub maximum_relative_balance_error: f64,
}

/// Source, physical-map, numerical or cooperative cancellation refusal.
#[derive(Debug)]
pub enum AssemblyMimoTransientError {
    Cancelled,
    Source(MimoLoadError),
    Structure(String),
    Flow(String),
}
impl core::fmt::Display for AssemblyMimoTransientError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "assembled MIMO transient refused: {self:?}")
    }
}
impl std::error::Error for AssemblyMimoTransientError {}
type Result<T> = core::result::Result<T, AssemblyMimoTransientError>;

/// One accepted interval in the original full component and port ordering.
#[derive(Debug)]
pub struct AssemblyMimoTransientStep {
    /// Endpoint displacement/rotation, including the already loaded mean shape.
    pub displacements: Vec<Vec<f64>>,
    pub velocities: Vec<Vec<f64>>,
    /// Physical displacement increment divided by the unchanged source clock.
    pub mean_velocities: Vec<Vec<f64>>,
    /// Resisting load force/moment ON the structure, -W*effort, at the endpoint.
    pub reaction_forces: Vec<Vec<f64>>,
    /// Exact interval load impulse, -W*integral(effort dt).
    pub reaction_impulses: Vec<Vec<f64>>,
    /// Physical W^T*v at the endpoint, independently checked against B^T*v_r.
    pub port_flows: Vec<f64>,
    pub port_efforts: Vec<f64>,
    /// Exact trajectory integrals; their product is NOT exchanged work.
    pub integrated_port_flows: Vec<f64>,
    pub integrated_port_efforts: Vec<f64>,
    pub perturbation_energy_j: f64,
    pub energy_change_j: f64,
    /// Original full component force/moment dotted with displacement increment.
    pub supplied_work_j: f64,
    /// Independent exact integral of the original v^T C_glue v, counted once.
    pub glue_dissipation_j: f64,
    /// Signed integral of effort^T flow. This is external boundary work, NOT a
    /// fitted load's invented storage or heat. Existing history may release it.
    pub external_load_work_j: f64,
    /// Delta H_structure + glue heat + external load work - supplied work.
    pub energy_balance_residual_j: f64,
}

/// Fixed-clock exact motion of one immutable accepted body and original load.
///
/// The original source's declared band and continuous/sampled report classes
/// remain unchanged. Exact evolution of its supplied rational equations does
/// not certify that physical source outside its band, a broadband excitation,
/// omitted structural modes, nonlinear motion or a changed contact branch.
/// Source memory has no inferred positive storage/heat certificate here.
///
/// Preparation and warm contractions have the original bounded atomic units.
/// Structural state AND original MIMO history publish only after all physical
/// maps, independent accounts and the final cancellation check succeed.
pub struct PreparedAssemblyMimoTransient<'p, 'm, 'v, 'body> {
    ports: &'p AssemblyHarmonicPorts<'m, 'v, 'body>,
    structure: AssemblyResponseSource<'v, 'body>,
    source: PreparedMimoLoad,
    joint: PreparedJointFlow,
    state: Vec<f64>,
    tolerance: f64,
}

impl<'p, 'm, 'v, 'body> PreparedAssemblyMimoTransient<'p, 'm, 'v, 'body> {
    /// Consume the already admitted proper MIMO source at its original clock.
    /// Its original history is retained; structural perturbation starts at rest.
    /// Port units/order must match the source. The same accepted assembly that
    /// owns the cached physical ports supplies M, C, K and the glue witness.
    pub fn new(
        ports: &'p AssemblyHarmonicPorts<'m, 'v, 'body>,
        source: PreparedMimoLoad,
        budget: AssemblyMimoTransientBudget,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        let mut cancelled = || cx.checkpoint().is_err();
        check(&mut cancelled)?;
        let shapes = &ports.physical_ports()[0].component_force_shapes;
        let physical = shapes.iter().map(|v| v.len() as u128).sum();
        joint::preflight(
            ports.vibration().mode_count(),
            ports.port_count(),
            &source,
            physical,
            shapes.len(),
            budget,
        )?;
        let structure = AssemblyResponseSource::prepare(ports.vibration(), budget.structure, cx)
            .map_err(AssemblyMimoTransientError::Structure)?;
        let joint = PreparedJointFlow::new(
            structure.system(),
            structure.stiffness_matrix(),
            structure.inverse_mass_matrix(),
            structure.damping_matrix(),
            ports.force_map(),
            &source,
            budget,
            &mut cancelled,
        )?;
        let state = vec![0.; structure.system().state_dim()];
        check(&mut cancelled)?;
        Ok(Self {
            ports,
            structure,
            source,
            joint,
            state,
            tolerance: budget.maximum_relative_balance_error,
        })
    }

    pub fn duration_s(&self) -> f64 {
        self.source.time_step_s()
    }

    /// Original immutable source, including accepted pole history and band.
    pub fn source(&self) -> &PreparedMimoLoad {
        &self.source
    }

    /// Continuous and sampled reports; GridOnly remains GridOnly.
    pub fn source_reports(&self) -> (&MimoPassivityReport, &MimoPassivityReport) {
        self.source.reports()
    }

    pub fn structure(&self) -> &AssemblyResponseSource<'v, 'body> {
        &self.structure
    }

    /// Structural [q,p] with p=M*v. Copy together with source real coordinates
    /// for an exact checkpoint on THIS body, mode order, load and clock.
    pub fn canonical_state(&self) -> &[f64] {
        &self.state
    }

    /// Restore the complete joint state atomically in the unchanged source chart.
    /// Pair coordinates follow MimoContinuous::real_coordinates (Re, -Im).
    pub fn restore_state(&mut self, structural: &[f64], load: &[f64], cx: &Cx<'_>) -> Result<()> {
        let mut cancelled = || cx.checkpoint().is_err();
        check(&mut cancelled)?;
        if structural.len() != self.state.len()
            || load.len() != self.source.continuous_operator().dimension()
            || structural.iter().any(|v| !v.is_finite())
            || load.iter().any(|v| !v.is_finite())
            || !self.structure.system().hamiltonian(structural).is_finite()
        {
            return Err(invalid("structural checkpoint must be complete and finite"));
        }
        let mut state = structural.to_vec();
        state.extend_from_slice(load);
        let flow = self.joint.endpoint_flows(&state)?;
        let candidate = self
            .source
            .continuous_operator()
            .candidate_from_coordinates(load, &flow)
            .map_err(AssemblyMimoTransientError::Source)?;
        state.truncate(structural.len());
        check(&mut cancelled)?;
        self.state = state;
        self.source = candidate;
        Ok(())
    }

    /// Initialize structural perturbation while retaining every original load
    /// memory coordinate. No existing radiation/material history is zeroed.
    pub fn initialize(&mut self, q: &[f64], v: &[f64], cx: &Cx<'_>) -> Result<()> {
        check(&mut || cx.checkpoint().is_err())?;
        let n = self.ports.vibration().mode_count();
        if q.len() != n || v.len() != n || q.iter().chain(v).any(|v| !v.is_finite()) {
            return Err(invalid(
                "initial perturbation must match the finite structural basis",
            ));
        }
        let mut state = q.to_vec();
        state.extend(multiply(self.structure.mass_matrix(), n, v)?);
        let load = self.source.continuous_operator().real_coordinates();
        self.restore_state(&state, &load, cx)
    }

    /// Advance one held physical incremental force/moment vector. Static loads
    /// are already balanced; the continuous load sees the evolving actual flow.
    pub fn step(&mut self, loads: &[&[f64]], cx: &Cx<'_>) -> Result<AssemblyMimoTransientStep> {
        self.step_with_checkpoint(loads, cx, &mut || cx.checkpoint().is_err())
    }

    fn step_with_checkpoint(
        &mut self,
        loads: &[&[f64]],
        cx: &Cx<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AssemblyMimoTransientStep> {
        check(cancelled)?;
        let vibration = self.ports.vibration();
        let n = vibration.mode_count();
        let force = vibration
            .project_loads(loads, cx)
            .map_err(AssemblyMimoTransientError::Structure)?;
        let mut initial = self.state.clone();
        initial.extend(self.source.continuous_operator().real_coordinates());
        let candidate = self
            .joint
            .flow
            .evaluate(&initial, &force)
            .map_err(|e| AssemblyMimoTransientError::Flow(e.to_string()))?;
        check(cancelled)?;
        let next = &candidate.state[..2 * n];
        let velocity = multiply(self.structure.inverse_mass_matrix(), n, &next[n..])?;
        let delta: Vec<_> = next[..n]
            .iter()
            .zip(&self.state[..n])
            .map(|(a, b)| a - b)
            .collect();
        let mean: Vec<_> = delta.iter().map(|v| v / self.duration_s()).collect();
        let port_efforts = self.joint.endpoint_efforts(&candidate.state)?;
        let retained_flows = self.joint.endpoint_flows(&candidate.state)?;
        let staged_source = self
            .source
            .continuous_operator()
            .candidate_from_coordinates(&candidate.state[2 * n..], &retained_flows)
            .map_err(AssemblyMimoTransientError::Source)?;
        let integral_velocity = &candidate.integrals[2..2 + n];
        let integrated_port_efforts = candidate.integrals[2 + n..].to_vec();
        // The exact integral of structural velocity is its displacement
        // increment. This also avoids subtractive Gramian/endpoint comparison
        // at a complete oscillation, where both integrals vanish physically.
        let integrated_port_flows = self.joint.project_velocities(&delta)?;
        let glue = candidate.integrals[0];
        let load_work = candidate.integrals[1];
        let modal_work = dot(&force, &delta)?;
        let integrated_work = dot(&force, integral_velocity)?;
        let energy = self.structure.system().hamiltonian(next);
        let delta_energy = joint::energy_change(
            self.structure.stiffness_matrix(),
            self.structure.inverse_mass_matrix(),
            &self.state,
            next,
        )?;
        let displacements = vibration
            .displacements(&next[..n], cx)
            .map_err(AssemblyMimoTransientError::Structure)?;
        let velocities = vibration
            .perturbation_displacements(&velocity, cx)
            .map_err(AssemblyMimoTransientError::Structure)?;
        let mean_velocities = vibration
            .perturbation_displacements(&mean, cx)
            .map_err(AssemblyMimoTransientError::Structure)?;
        let physical_delta = vibration
            .perturbation_displacements(&delta, cx)
            .map_err(AssemblyMimoTransientError::Structure)?;
        let mut supplied = 0.;
        let mut work_scale = 0.;
        for (load, movement) in loads.iter().zip(&physical_delta) {
            for (index, (&f, &dq)) in load.iter().zip(movement).enumerate() {
                if index % 1024 == 0 {
                    check(cancelled)?;
                }
                let term = f * dq;
                supplied += term;
                work_scale += term.abs();
            }
        }
        let mut reactions: Vec<_> = velocities.iter().map(|v| vec![0.; v.len()]).collect();
        let mut impulses = reactions.clone();
        let mut port_flows = vec![0.; self.ports.port_count()];
        for (a, port) in self.ports.physical_ports().iter().enumerate() {
            check(cancelled)?;
            let (mut flow_integral, mut flow_scale, mut integral_scale) = (0., 0., 0.);
            for (component, shape) in port.component_force_shapes.iter().enumerate() {
                for (j, &w) in shape.iter().enumerate() {
                    if j % 1024 == 0 {
                        check(cancelled)?;
                    }
                    if w == 0. {
                        continue;
                    }
                    let flow = w * velocities[component][j];
                    let integral = w * physical_delta[component][j];
                    port_flows[a] += flow;
                    flow_scale += flow.abs();
                    flow_integral += integral;
                    integral_scale += integral.abs();
                    reactions[component][j] -= w * port_efforts[a];
                    impulses[component][j] -= w * integrated_port_efforts[a];
                }
            }
            if !flow_scale.is_finite()
                || !integral_scale.is_finite()
                || !port_flows[a].is_finite()
                || !flow_integral.is_finite()
                || (port_flows[a] - retained_flows[a]).abs()
                    > self.tolerance * flow_scale.max(retained_flows[a].abs())
                || (flow_integral - integrated_port_flows[a]).abs()
                    > self.tolerance * integral_scale.max(integrated_port_flows[a].abs())
            {
                return Err(invalid(
                    "original physical attachment flow/work map refused",
                ));
            }
        }
        let scale = energy
            .abs()
            .max(self.structure.system().hamiltonian(&self.state).abs())
            .max(work_scale)
            .max(glue.abs())
            .max(load_work.abs())
            .max(modal_work.abs())
            .max(integrated_work.abs());
        let residual = delta_energy + glue + load_work - supplied;
        if !scale.is_finite()
            || !energy.is_finite()
            || !residual.is_finite()
            || !supplied.is_finite()
            || glue < 0.
            || residual.abs() > self.tolerance * scale
            || (supplied - modal_work).abs() > self.tolerance * scale
            || (supplied - integrated_work).abs() > self.tolerance * scale
            || displacements
                .iter()
                .chain(&velocities)
                .chain(&mean_velocities)
                .chain(&reactions)
                .chain(&impulses)
                .flatten()
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "original physical force, glue and load work balance refused",
            ));
        }
        let state = next.to_vec();
        check(cancelled)?;
        self.state = state;
        self.source = staged_source;
        Ok(AssemblyMimoTransientStep {
            displacements,
            velocities,
            mean_velocities,
            reaction_forces: reactions,
            reaction_impulses: impulses,
            port_flows,
            port_efforts,
            integrated_port_flows,
            integrated_port_efforts,
            perturbation_energy_j: energy,
            energy_change_j: delta_energy,
            supplied_work_j: supplied,
            glue_dissipation_j: glue,
            external_load_work_j: load_work,
            energy_balance_residual_j: residual,
        })
    }
}

#[cfg(test)]
mod tests;
