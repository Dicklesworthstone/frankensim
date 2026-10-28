//! Source viscosity on the original simultaneous modal/contact/radiation owner.
//!
//! Each physical member port contributes `c b b^T` to damping. Project its SAME
//! work-conjugate displacement map through the retained frame basis. The original
//! coupled modal propagator then carries all cross terms and independently
//! integrates `c (b^T v)^2`; no diagonal loss fit or second integrator is used.

use super::{FrameVibrationBasis, FrameVibrationError, checked_sum, poll};
use crate::modal_acoustic_time::{
    ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use crate::render::schedule::force::coupled::{
    CoupledModalSystem, MAX_MODAL_CONNECTIONS, ModalAttachment, ModalConnection,
    ModalCouplingConfig, PreparedModalBudget,
};
use fs_exec::CancelGate;

/// One complete frame component: its unchanged modal basis, original sources,
/// oscillator state, and EVERY projected material loss port travel together.
/// Additional authored per-mode damping remains a distinct loss mechanism.
pub struct FrameVibrationComponent {
    basis: FrameVibrationBasis,
    model: ModalAcousticTimeModel,
    connections: Vec<ModalConnection>,
}

impl FrameVibrationBasis {
    /// Bind the same physical frame to the shared simultaneous modal owner.
    ///
    /// Every member's original axial and bending viscosity enters as a positive
    /// physical resistance port, retaining directional and inter-mode terms.
    /// The source states/geometry stay frozen; material loss is exported by the
    /// original connection dissipation account, not deposited into a heat cell.
    /// `additional_damping_ratios` describes a separate authored mechanism; pass
    /// zeros when all loss is supplied by the retained material sections.
    ///
    /// The declared `max_projection_terms` bounds six scalar products per
    /// retained mode per material port BEFORE allocating projected columns.
    /// Native modal and coupled setup budgets still apply. Explicit truncation
    /// produces the Galerkin operator `Phi^T C Phi`, retaining its cross terms;
    /// it does not recover the dynamics or dissipation of omitted coordinates.
    ///
    /// No pressure transfer is invented. This component can join the existing
    /// contact, friction, radiation, and other modal components through its
    /// ordinary reciprocal attachments. Use `into_system` for one frame, or
    /// include ALL `into_parts` connections in a larger original coupled system.
    pub fn coupled_component(
        &self,
        sample_rate_hz: u32,
        additional_damping_ratios: &[f64],
        modal_budget: ModalAcousticTimeBudget,
        max_projection_terms: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FrameVibrationComponent, FrameVibrationError> {
        poll(&mut cancelled)?;
        let model =
            self.oscillator_model(sample_rate_hz, additional_damping_ratios, modal_budget)?;
        let ports = self
            .material_damping
            .as_ref()
            .map_or(0, |d| d.ports().len());
        if ports > MAX_MODAL_CONNECTIONS
            || ports
                .checked_mul(self.modes.len())
                .and_then(|n| n.checked_mul(6))
                .is_none_or(|n| n > max_projection_terms)
        {
            return Err(FrameVibrationError::Invalid(
                "frame material-port projection budget",
            ));
        }
        let mut connections = Vec::with_capacity(ports);
        if let Some(operator) = &self.material_damping {
            operator.admit_timestep(1.0 / f64::from(sample_rate_hz))?;
            for port in operator.ports() {
                poll(&mut cancelled)?;
                let mut shapes = Vec::with_capacity(self.modes.len());
                for mode in &self.modes {
                    poll(&mut cancelled)?;
                    shapes.push(checked_sum(
                        port.dofs
                            .iter()
                            .zip(&port.displacement_weights)
                            .map(|(&dof, &weight)| weight * mode.displacement[dof]),
                    )?);
                }
                connections.push(ModalConnection {
                    left: ModalAttachment {
                        component: 0,
                        shapes,
                    },
                    right: ModalAttachment {
                        component: 0,
                        shapes: vec![0.0; self.modes.len()],
                    },
                    stiffness_n_m: 0.0,
                    damping_n_s_m: port.damping_n_s_m,
                    rest_extension_m: 0.0,
                });
            }
        }
        let basis = self.clone();
        poll(&mut cancelled)?;
        Ok(FrameVibrationComponent {
            basis,
            model,
            connections,
        })
    }
}

impl FrameVibrationComponent {
    /// Frozen reference, reciprocal physical maps and original material sources.
    pub fn basis(&self) -> &FrameVibrationBasis {
        &self.basis
    }

    /// Original accepted coordinates, in this retained mass-normalized basis.
    pub fn states(&self) -> &[ModalAcousticState] {
        self.model.states()
    }

    /// Restore an explicit modal state before network preparation. The original
    /// modal owner validates it transactionally; sources and ports are unchanged.
    pub fn restore_states(
        &mut self,
        states: &[ModalAcousticState],
    ) -> Result<(), FrameVibrationError> {
        self.model.restore_states(states)?;
        Ok(())
    }

    /// Original physical-port order, three channels per source-bound member.
    /// Indices correspond to `basis().material_damping().ports()`. All shapes
    /// currently refer to component zero; `into_parts` remaps them together.
    pub fn material_connections(&self) -> &[ModalConnection] {
        &self.connections
    }

    /// Bind into a larger shared modal/contact/radiation assembly. Append the
    /// model at `component`, include EVERY returned connection, and retain the
    /// basis for physical force/observation maps and source provenance. Discarding
    /// a returned material connection changes the physical model.
    pub fn into_parts(
        mut self,
        component: usize,
    ) -> (
        FrameVibrationBasis,
        ModalAcousticTimeModel,
        Vec<ModalConnection>,
    ) {
        for connection in &mut self.connections {
            connection.left.component = component;
            connection.right.component = component;
        }
        (self.basis, self.model, self.connections)
    }

    /// Prepare the ORIGINAL simultaneous coupled owner with every material port.
    /// The returned basis retains source records and maps physical loads and
    /// observations of `system.components()[0]`. Nonzero initial states survive.
    ///
    /// With `prepared`, use its existing cached exact propagator; otherwise use
    /// its ordinary bounded complete propagator. Material loss appears in
    /// `connection_dissipation_j`, additional per-mode loss in
    /// `component_dissipation_j`. These are parts of the existing total account.
    pub fn into_system(
        self,
        config: ModalCouplingConfig,
        prepared: Option<PreparedModalBudget>,
        gate: &CancelGate,
    ) -> Result<(FrameVibrationBasis, CoupledModalSystem), FrameVibrationError> {
        let (basis, model, connections) = self.into_parts(0);
        let system = match prepared {
            Some(budget) => {
                CoupledModalSystem::new_prepared(vec![model], connections, config, budget, gate)?
            }
            None => CoupledModalSystem::new(vec![model], connections, config, gate)?,
        };
        Ok((basis, system))
    }
}
